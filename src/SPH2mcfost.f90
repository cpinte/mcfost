module SPH2mcfost

  use parametres
  use constantes
  use utils
  use density, only : normalize_dust_density, reduce_density
  use read_phantom, only : read_phantom_bin_files, read_phantom_hdf_files

  implicit none

  procedure(read_phantom_bin_files), pointer :: read_phantom_files => null()

contains

  subroutine setup_SPH2mcfost(SPH_file,SPH_limits_file, n_SPH, extra_heating)

    use read_gadget2, only : read_gadget2_file
    use dump_utils, only : get_error_text
    use utils, only : read_comments

    character(len=512), intent(in) :: SPH_file, SPH_limits_file
    integer, intent(out) :: n_SPH

    integer, parameter :: iunit = 1

    real(dp), allocatable, dimension(:) :: x,y,z,h,vx,vy,vz,rho,massgas,SPH_grainsizes,T_gas
    integer,  allocatable, dimension(:) :: particle_id
    real(dp), allocatable, dimension(:,:) :: rhodust, massdust
    real, allocatable, dimension(:) :: extra_heating
    logical, allocatable, dimension(:) :: mask ! size == np, not n_SPH, index is original SPH id

    real(dp), dimension(6) :: SPH_limits
    real :: factor
    integer :: ndusttypes, ierr, i, ilen
    character(len=100) :: line_buffer
    logical :: check_previous_tesselation

    if (lphantom_file) then
       write(*,*) "Performing phantom2mcfost setup"
       if (n_phantom_files==1) then
          write(*,*) "Reading phantom density file: "//trim(density_files(1))  ! todo : update: we do not use SPH_file anymore
       else
          write(*,*) "Reading phantom density files: "
          do i=1,n_phantom_files
             write(*,*) " - "//trim(density_files(i))  ! todo : update: we do not use SPH_file anymore
          enddo
       endif

       ! Are we reading phantom binary or HDF5 files ?
       ilen = index(density_files(1),'.',back=.true.) ! last position of the '.' character
       if (ilen > 0) then
          if (density_files(1)(ilen:ilen+3) == ".h5") then
             read_phantom_files => read_phantom_hdf_files
          else
             read_phantom_files => read_phantom_bin_files
          endif
       else
          read_phantom_files => read_phantom_bin_files
       endif

       call read_phantom_files(iunit,n_phantom_files,density_files, x,y,z,h,vx,vy,vz,T_gas, &
            particle_id, massgas,massdust,rho,rhodust,extra_heating,ndusttypes, &
            SPH_grainsizes,mask,n_SPH,ierr)

       if (lphantom_avg) then ! We are averaging the dump
          factor = 1.0/n_phantom_files
          massgas = massgas * factor
          massdust = massdust * factor
          rho = rho * factor
          rhodust = rhodust * factor
          h = h/(n_phantom_files**(1./3.))
       endif

       ! Todo : extra heating must be passed to mcfost
       if (ierr /=0) then
          write(*,*) "Error code =", ierr,  get_error_text(ierr)
          call exit(1)
       endif
    else if (lgadget2_file) then
       write(*,*) "Performing gadget2mcfost setup"
       write(*,*) "Reading Gadget-2 density file: "//trim(SPH_file)
       call read_gadget2_file(iunit,SPH_file, x,y,z,h,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
       T_gas = x       !to allocate memory
       T_gas = 2.74
    else if (lascii_SPH_file) then
       write(*,*) "Performing SPH2mcfost setup"
       write(*,*) "Reading SPH density file: "//trim(SPH_file)
       call read_ascii_SPH_file(iunit,SPH_file, x,y,z,h,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
       T_gas = x       !to allocate memory
       T_gas = 2.74
    endif
    write(*,*) "Done"

    if (lignore_dust) then
       ndusttypes = 0
       if (allocated(rhodust)) deallocate(rhodust,massdust)
    endif

    if ((.not.lfix_star).and.(lphantom_file .or. lgadget2_file)) call compute_stellar_parameters()

    ! Model limits
    call read_SPH_limits_file(SPH_limits_file, SPH_limits)

    ! Voronoi tesselation
    check_previous_tesselation = (.not. lrandomize_Voronoi)
    call SPH_to_Voronoi(n_SPH, ndusttypes, particle_id, x,y,z,h, vx,vy,vz, T_gas, massgas,massdust,rho,rhodust,SPH_grainsizes, SPH_limits, check_previous_tesselation, mask=mask)

    deallocate(x,y,z,h)
    if (allocated(vx)) deallocate(vx,vy,vz)
    deallocate(massgas,rho)
    if (allocated(rhodust)) deallocate(rhodust,massdust)

    ! Deleting partcles/cells in masked arreas (Hill sphere, etc)
    if (allocated(mask)) call delete_masked_particles()

    return

  end subroutine setup_SPH2mcfost

  !*********************************************************

  subroutine read_SPH_limits_file(SPH_limits_file, SPH_limits)

    character(len=512), intent(in) :: SPH_limits_file
    real(dp), dimension(6), intent(out) :: SPH_limits

    integer :: ios

    write(*,*) " "
    if (llimits_file) then
       write(*,*) "Reading limits file: "//trim(SPH_limits_file)
       open(unit=1, file=SPH_limits_file, status='old', iostat=ios)
       if (ios/=0) call error("cannot open "//trim(SPH_limits_file))
       call read_comments(1)
       read(1,*) SPH_limits(1), SPH_limits(3), SPH_limits(5)
       read(1,*) SPH_limits(2), SPH_limits(4), SPH_limits(6)
       close(unit=1)
    else
       SPH_limits(:) = 0
    endif

    return

  end subroutine read_SPH_limits_file

  !*********************************************************

  subroutine SPH_to_Voronoi(n_SPH, ndusttypes, particle_id, x,y,z,h, vx,vy,vz, T_gas, massgas,massdust,rho,rhodust,&
       SPH_grainsizes, SPH_limits, check_previous_tesselation, mask)

    use Voronoi_grid
    use density, only : densite_gaz, masse_gaz, densite_pouss, masse
    use grains, only : n_grains_tot, M_grain
    use disk_physics, only : compute_othin_sublimation_radius
    use mem

    integer, intent(in) :: n_SPH, ndusttypes
    real(dp), dimension(n_SPH), intent(inout) :: x,y,z,h,rho,massgas
    real(dp), dimension(:), allocatable, intent(inout) :: vx,vy,vz ! dimension n_SPH or 0
    real(dp), dimension(:), allocatable, intent(in) :: T_gas
    integer, dimension(n_SPH), intent(in) :: particle_id
    real(dp), dimension(ndusttypes,n_SPH), intent(in) :: rhodust, massdust
    real(dp), dimension(ndusttypes), intent(in) :: SPH_grainsizes
    real(dp), dimension(6), intent(in) :: SPH_limits
    logical, intent(in) :: check_previous_tesselation
    logical, dimension(:), allocatable, intent(in), optional :: mask

    logical :: lwrite_ASCII = .false. ! produce an ASCII file for yorick

    real, allocatable, dimension(:) :: a_SPH, log_a_SPH, rho_dust
    real(dp) :: mass, somme, Mtot, Mtot_dust
    real :: f, limit_threshold, density_factor, destruction_factor
    integer :: icell, l, k, iSPH, n_force_empty, i, id_n

    real(dp), dimension(6) :: limits


    if (lcorrect_density_elongated_cells) then
       density_factor = correct_density_factor_elongated_cells
    else
       density_factor = 1
    endif

    limit_threshold = (1.0 - SPH_keep_particles) * 0.5 ;

    icell_ref = 1

    write(*,*) "# Farthest particules :"
    write(*,*) "x =", minval(x), maxval(x)
    write(*,*) "y =", minval(y), maxval(y)
    write(*,*) "z =", minval(z), maxval(z)

    write(*,*) "Found", n_SPH, " SPH particles with ", ndusttypes, "dust grains"

    if (lwrite_ASCII) then
       !  Write the file for the grid version of mcfost
       !- N_part: total number of particles
       !  - r_in: disk inner edge in AU
       !  - r_out: disk outer edge in AU
       !  - p: surface density exponent, Sigma=Sigma_0*(r/r_0)^(-p), p>0
       !  - q: temperature exponent, T=T_0*(r/r_0)^(-q), q>0
       !  - m_star: star mass in solar masses
       !  - m_disk: disk mass in solar masses (99% gas + 1% dust)
       !  - H_0: disk scale height at 100 AU, in AU
       !  - rho_d: dust density in g.cm^-3
       !  - flag_ggrowth: T with grain growth, F without
       !
       !
       !    N_part lines containing:
       !  - x,y,z: coordinates of each particle in AU
       !  - h: smoothing length of each particle in AU
       !  - s: grain size of each particle in µm
       !
       !  Without grain growth: 2 lines containing:
       !  - n_sizes: number of grain sizes
       !  - (s(i),i=1,n_sizes): grain sizes in µm
       !  OR
       !  With grain growth: 1 line containing:
       !  - s_min,s_max: smallest and largest grain size in µm

       open(unit=1,file="SPH_phantom.txt",status="replace")
       write(1,*) size(x)
       write(1,*) minval(sqrt(x**2 + y**2))
       write(1,*) maxval(sqrt(x**2 + y**2))
       write(1,*) 1 ! p
       write(1,*) 0.5 ! q
       write(1,*) 1.0 ! mstar
       write(1,*) 1.e-3 !mdisk
       write(1,*) 10 ! h0
       write(1,*) 3.5 ! rhod
       write(1,*) .false.
       !rhoi = massoftype(itypei)*(hfact/hi)**3  * udens ! g/cm**3

       do icell=1,size(x)
          write(1,*) x(icell), y(icell), z(icell), 1.0, 1.0
       enddo

       write(1,*) 1
       write(1,*) 1.0
       close(unit=1)
    endif

    if (abs(maxval(SPH_limits)) < tiny_real) then
       write(*,*) "Selecting spatial range which contains"
       write(*,*) SPH_keep_particles*100, "% of particles in each dimension"

       k = max(int(limit_threshold * n_SPH),1)
       limits(1) = select_inplace(k,real(x))
       limits(3) = select_inplace(k,real(y))
       limits(5) = select_inplace(k,real(z))

       k = int((1.0-limit_threshold) * n_SPH)
       limits(2) = select_inplace(k,real(x))
       limits(4) = select_inplace(k,real(y))
       limits(6) = select_inplace(k,real(z))
    else
       limits(:) = SPH_limits(:)
    endif

    write(*,*) "# Model limits :"
    write(*,*) "x =", limits(1), limits(2)
    write(*,*) "y =", limits(3), limits(4)
    write(*,*) "z =", limits(5), limits(6)

    !*******************************
    ! Voronoi tesselation
    !*******************************
    ! Make the Voronoi tesselation on the SPH particles ---> define_Voronoi_grid : volume
    !call Voronoi_tesselation_cmd_line(n_SPH, x,y,z, limits)

    call Voronoi_tesselation(n_SPH, particle_id, x,y,z,h,vx,vy,vz, limits, check_previous_tesselation)
    !deallocate(x,y,z)
    write(*,*) "Using n_cells =", n_cells


    !*************************
    ! Densities
    !*************************
    call allocate_densities(n_cells_max = n_SPH + n_etoiles) ! we allocate all the SPH particule for libmcfost
    ! Tableau de densite et masse de gaz
    !do icell=1,n_cells
    !   densite_gaz(icell) = rho(icell) / masse_mol_gaz * m3_to_cm3 ! rho is in g/cm^3 --> part.m^3
    !   masse_gaz(icell) =  densite_gaz(icell) * masse_mol_gaz * volume(icell)
    !enddo
    !masse_gaz(:) = masse_gaz(:) * AU3_to_cm3

    ! I need to work with masses, as Voronoi and SPH volume could be different
    do icell=1,n_cells
       iSPH = Voronoi(icell)%id
       if (iSPH > 0) then
          masse_gaz(icell)    = massgas(iSPH) * Msun_to_g ! en g
          densite_gaz(icell)  = masse_gaz(icell) /  (masse_mol_gaz * volume(icell) * AU3_to_m3)
       else ! star
          masse_gaz(icell)    = 0.
          densite_gaz(icell)  = 0.
       endif
    enddo

    ! Tableau de densite et masse de poussiere
    ! interpolation en taille
    if (ndusttypes >= 1) then
       lvariable_dust = .true.

       ! mcfost adds an extra grain size follwing the gas
       ! this is required if ndusttypes == 1, but I do it in any case
       allocate(a_SPH(ndusttypes+1),log_a_SPH(ndusttypes+1),rho_dust(ndusttypes+1))

       write(*,*) "Found the following grain sizes in SPH calculation:"
       do l=1, ndusttypes
          if (dust_pop(1)%porosity > tiny_real) then
             if (l==1) call warning("Grain sizes are adjusted for porosity")
             ! Stokes number is reduced by porosity
             ! We shift the mcfost grain sizes relative to phantom grain sizes
             a_SPH(l+1) = SPH_grainsizes(l) / (1.-dust_pop(1)%porosity)
          else
             a_SPH(l+1) = SPH_grainsizes(l)
          endif
          if (lfluffy) then
             if (l==1) call warning("Grain sizes are adjusted for fluffyness")
             a_SPH(l+1) = a_SPH(l+1) / fluffyness
          endif
          write(*,*) real(SPH_grainsizes(l)), "microns   --->  ",  real(a_SPH(l+1)), "microns"
       enddo

       ! Adding smallest grain size following the gas
       if (a_SPH(2) <= 1.0+1e-6) then
          write(*,*) "WARNING: assuming dust grains smaller than 0.1mum are following the gas"
          a_SPH(1) = 0.1 ;
       else
          write(*,*) "WARNING: assuming dust grains smaller than 1mum are following the gas"
          a_SPH(1) = 1. ;
       endif

       ! temporary for old phantom dumps were grain sizes were not defined
       if (SPH_grainsizes(1) < tiny_real) then
          a_SPH(1) = 1. ;
          a_SPH(2) = 1000. ;
          write(*,*) "WARNING: SPH dust grain size  not found"
          write(*,*) "Forcing small & big grains to be ", a_SPH(1), "and", a_SPH(2), "mum"
       endif

       log_a_SPH(:) = 0.
       rho_dust(:) = 0.

       ! We are making sure that the grain sizes are sorted
       do l=1, ndusttypes+1
          if (l<=ndusttypes) then
             if (a_sph(l) >= a_sph(l+1)) then
                write(*,*) "ERROR : grains must be ordered from small to large"
                do k=1, ndusttypes
                   write(*,*) k, a_sph(k)
                enddo
                write(*,*) "Exiting"
                call exit(1)
             endif
          endif

          if (a_SPH(l) > 0.) then
             log_a_SPH(l) = log(a_SPH(l))
          else
             call error("grains size must be > 0 in SPH file")
          endif
       enddo

       Mtot = 0. ; Mtot_dust = 0.
       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) then
             Mtot = Mtot +  massgas(icell)
             !Mtot_dust = Mtot_dust + masse_gaz(icell)* (sum(dustfrac(:,iSPH)))
             Mtot_dust = Mtot_dust + sum(massdust(:,icell))
          endif
       enddo

       !Mtot_dust = Mtot_dust * Msun_to_g

       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) then
             ! mass & density indices are shifted by 1
             do l=1, ndusttypes+1
                if (l==1) then
                   ! small grains follow the gas, we do not care about normalization here
                   rho_dust(l) = massgas(iSPH) / volume(icell)
                else
                   rho_dust(l) = massdust(l-1,iSPH) / volume(icell)
                endif
             enddo ! l

             l=1
             do k=1,n_grains_tot
                if (r_grain(k) < a_SPH(1)) then ! small grains
                   densite_pouss(k,icell) = rho_dust(1)
                else if (r_grain(k) > a_SPH(ndusttypes+1)) then ! large grains
                   densite_pouss(k,icell) = rho_dust(ndusttypes+1)
                else ! interpolation
                   if (r_grain(k) > a_sph(l+1)) l = l+1
                   f = (log(r_grain(k))-log_a_sph(l))/(log_a_sph(l+1)-log_a_sph(l))
                   densite_pouss(k,icell) = rho_dust(l) + f * (rho_dust(l+1)  - rho_dust(l))
                endif
             enddo !k
          else ! iSPH == 0, star
             densite_pouss(:,icell) = 0.
          endif
       enddo ! icell

       ! Using the parameter file gas-to-dust ratio for now
       ! until phantom provides a proper grain size distribution
       call normalize_dust_density( sum(masse_gaz) * g_to_Msun / disk_zone(1)%gas_to_dust)
    else ! ndusttypes = 0 : using the gas density

       lvariable_dust = .false.
       write(*,*) "Using gas-to-dust ratio in mcfost parameter file"

       destruction_factor = 1.
       do icell=1,n_cells
          if (T_gas(icell) >= 1500. .and. .not. lturn_off_dust_subl) destruction_factor = 0.0001
          masse(icell) = 0.
          do k=1,n_grains_tot
             densite_pouss(k,icell) = densite_gaz(icell) * nbre_grains(k) * destruction_factor
             masse(icell) = masse(icell) + densite_pouss(k,icell) * M_grain(k) * volume(icell)
          enddo
       enddo
       masse(:) = masse(:) * AU3_to_cm3
       f = 1./disk_zone(1)%gas_to_dust * sum(masse_gaz)/sum(masse)
       densite_pouss(:,:) = densite_pouss(:,:) * f
       masse(:) = masse(:) * f
    endif ! ndusttypes == 0

    !*************************
    ! Mask
    !*************************
    if (present(mask)) then
       if (allocated(mask)) then
          do icell=1,n_cells
             iSPH = Voronoi(icell)%original_id
             if (iSPH > 0) then
                Voronoi(icell)%masked = mask(iSPH)
             else
                Voronoi(icell)%masked = .false.
             endif
          enddo
       else
          do icell=1,n_cells
             Voronoi(icell)%masked = .false.
          enddo
       endif
    else
       do icell=1,n_cells
          Voronoi(icell)%masked = .false.
       enddo
    endif

    ! We eventually reduce density to avoid artefacts: superseeded by cell cutting
    if (density_factor < 1.-1e-6) then
       ! Removing cells at the "surface" of the SPH model:
       ! density is reduced so that they do not appear in images or cast artificial shadows,
       ! but we can still compute a temperature (forcing them to be optically thin)
       n_force_empty = 0
       cell_loop : do icell=1,n_cells
          ! We reduce the density on cells that are very elongated
          if (Voronoi(icell)%was_cut) then
             n_force_empty = n_force_empty + 1
             call reduce_density(icell, density_factor)
             cycle cell_loop
          endif

          ! We reduce the density on cells that are touching a wall
          do i=Voronoi(icell)%first_neighbour, Voronoi(icell)%last_neighbour
             id_n = neighbours_list(i) ! id du voisin
             if (id_n < 0) then
                n_force_empty = n_force_empty + 1
                call reduce_density(icell, density_factor)
                cycle cell_loop
             endif
          enddo
       enddo cell_loop
       write(*,*) "Density was reduced by", density_factor, "in", n_force_empty, "cells surrounding the model, ie", (1.0*n_force_empty)/n_cells * 100, "% of cells"
    endif

    write(*,*) 'Total  gas mass in model :',  real(sum(masse_gaz) * g_to_Msun),' Msun'
    write(*,*) 'Total dust mass in model :', real(sum(masse) * g_to_Msun),' Msun'

    search_not_empty : do k=1,n_grains_tot
       do icell=1, n_cells
          if (densite_pouss(k,icell) > 0.0_sp) then
             icell_not_empty = icell
             exit search_not_empty
          endif
       enddo !icell
    enddo search_not_empty

    if (ndusttypes >= 1) deallocate(a_SPH,log_a_SPH,rho_dust)

    return

  end subroutine SPH_to_Voronoi

  !*********************************************************

  subroutine compute_stellar_parameters()

    integer :: i

    character(len=512) :: isochrone_file, filename
    character(len=100) :: line_buffer
    character(len=1)   :: s_age

    character(len=2) :: SpT
    real :: L, R, T, M, minM, maxM, logg, minM_Allard, maxM_Allard, minM_Siess, maxM_Siess
    real(kind=dp) :: Gcm_to_Rsun
    integer :: age, k

    logical :: lread_Siess, lread_Allard
    integer, parameter :: nSpT_Siess = 29
    real, dimension(nSpT_Siess) :: logR_Siess, logTeff_Siess, logM_Siess
    integer, parameter :: nSpT_Allard = 50
    real, dimension(nSpT_Allard) :: logR_Allard, logTeff_Allard, logM_Allard

    if (n_etoiles < 1) return

    minM_Siess = 0.1
    maxM_siess = 7.0

    minM_Allard = 0.0005
    maxM_Allard = minM_Siess

    ! Which models do we need to read ?
    lread_Siess = .False. ; lread_Allard = .False.
    do i=1, n_etoiles
       M = etoile(i)%M
       if (M > minM_Siess)  then
          lread_Siess = .True.
       else
          lread_Allard = .True.
       endif
    enddo

    if (lread_Siess) then
       ! Siess models
       isochrone_file = "Siess/isochrone_"//trim(system_age)//".txt"
       write(*,*) "Reading isochrone file: "//trim(isochrone_file)
       filename = trim(mcfost_utils)//"/Isochrones/"//trim(isochrone_file)

       open(unit=1,file=filename,status="old")
       do i=1,3
          read(1,*) line_buffer
       enddo
       minM_Siess = 1.e30 ; maxM = 0 ;
       do i=1, nSpT_Siess
          read(1,*) SpT, L, r, T, M
          logR_Siess(i) = log(r) ; logTeff_Siess(i) = log(T) ; logM_Siess(i) = log(M)
       enddo
       close(unit=1)
    endif

    if (lread_Allard) then
       ! Allard models if mass < 0.1Msun
       isochrone_file = "Allard/model.AMES-Cond-2000.M-0.0.2MASS.AB"
       write(*,*) "Reading isochrone file: "//trim(isochrone_file)
       filename = trim(mcfost_utils)//"/Isochrones/"//trim(isochrone_file)

       s_age = system_age(1:1)
       read(s_age,*) age ! age is an int with age in Myr

       open(unit=1,file=filename,status="old")
       Gcm_to_Rsun = 1e9 * cm_to_m/Rsun
       ! Skipping age block
       do k=1,age-1
          ! header
          do i=1,4
             read(1,*) line_buffer
          enddo
          ! data
          line_loop : do i=1,nSpT_Allard
             read(1,'(A)') line_buffer
             if (line_buffer(1:1) == "-") exit line_loop
          enddo line_loop
       enddo

       ! header
       do i=1,4
          read(1,*) line_buffer
       enddo
       ! data
       k=0
       line_loop2 : do i=1,nSpT_Allard
          read(1,'(A)') line_buffer
          if (line_buffer(1:1) == "-") exit line_loop2
          k = k+1
          read(line_buffer,*) M, T, L, logg, R
          logR_Allard(i) = log(r * Gcm_to_Rsun) ; logTeff_Allard(i) = log(T) ; logM_Allard(i) = log(M)
       enddo line_loop2
       close(unit=1)
    endif

    ! interpoler L et T, les fonctions sont plus smooth
    write(*,*) ""
    write(*,*) "New stellar parameters:"
    do i=1, n_etoiles
       if (lturn_off_planets .and. i>1) then
          write(*,*) " "
          write(*,*) "*** WARNING : turning off emission fron sink particle"
          write(*,*) "*** object #", i, "M=", etoile(i)%M, "Msun"
          write(*,*) "*** The object will not radiate"
          etoile(i)%T = 3.
          etoile(i)%r = 1e-4
       else if (etoile(i)%M < minM_Allard) then
          write(*,*) " "
          write(*,*) "*** WARNING : stellar object mass is below isochrone range"
          write(*,*) "*** object #", i, "M=", etoile(i)%M, "Msun"
          write(*,*) "*** The object will not radiate"
          etoile(i)%T = 3.
          etoile(i)%r = 0.01
       else if (etoile(i)%M < maxM_Allard) then
          etoile(i)%T = exp(interp(logTeff_Allard(1:k), logM_Allard(1:k), log(etoile(i)%M)))
          etoile(i)%r = exp(interp(logR_Allard(1:k), logM_Allard(1:k), log(etoile(i)%M)))
       else ! using Siess' models
          if (etoile(i)%M > maxM_Siess) then
             write(*,*) " "
             write(*,*) "*** WARNING : stellar object mass is above in isochrone range"
             write(*,*) "*** object #", i, "M=", etoile(i)%M, "Msun"
             write(*,*) "*** Stellar properties are extrapolated"
          endif
          etoile(i)%T = exp(interp(logTeff_Siess, logM_Siess, log(etoile(i)%M)))
          etoile(i)%r = exp(interp(logR_Siess, logM_Siess, log(etoile(i)%M)))
       endif

       ! Pas de fUV et pas de spectre stellaire pour le moment
       etoile(i)%fUV = 0.0 ; etoile(i)%slope_UV = 0.0 ;
       etoile(i)%lb_body = .false.
       etoile(i)%spectre = "None"

       write(*,*) "Star #",i,"  Teff=", etoile(i)%T, "K, r=", etoile(i)%r, "Rsun"
    enddo

    ! Passage rayon en AU
    etoile(:)%r = etoile(:)%r * Rsun_to_AU

    return

  end subroutine compute_stellar_parameters

  !*********************************************************

  subroutine delete_masked_particles()

    use Voronoi_grid
    use density, only : densite_gaz, masse_gaz, densite_pouss, masse

    integer :: icell, k

    k=0
    do icell=1, n_cells
       if (Voronoi(icell)%masked) then
          k=k+1
          masse_gaz(icell)       = 0.
          densite_gaz(icell)     = 0.
          masse(icell)           = 0.
          densite_pouss(:,icell) = 0.
       endif
    enddo

    write(*,*) k, "masked cells have been made transparent"

    return

  end subroutine delete_masked_particles

  !*********************************************************

  subroutine delete_Hill_sphere()

    use Voronoi_grid
    use density, only : densite_gaz, masse_gaz, densite_pouss, masse

    integer :: istar, icell, n_delete
    real(kind=dp) :: d2, r_Hill2, r_hill, dx, dy, dz

    ! We assume that the central star is the actual star
    ! and following sink particles are planets
    do istar=2, n_etoiles
       n_delete = 0

       d2 = (etoile(istar)%x - etoile(1)%x)**2 + &
            (etoile(istar)%y - etoile(1)%y)**2 + &
            (etoile(istar)%z - etoile(1)%z)**2
       r_Hill2 = d2 * (etoile(istar)%m / (3*etoile(1)%m))**(2./3)
       r_Hill = sqrt(r_Hill2)

       write(*,*) "Sink particle #", istar, "Hill radius =", r_Hill, "au"

       cell_loop : do icell=1, n_cells
          if (icell == istar) cycle cell_loop
          dx = Voronoi(icell)%xyz(1) - etoile(istar)%x
          if (dx > r_Hill) cycle cell_loop
          dy = Voronoi(icell)%xyz(2) - etoile(istar)%y
          if (dy > r_Hill) cycle cell_loop
          dz = Voronoi(icell)%xyz(3) - etoile(istar)%z
          if (dz > r_Hill) cycle cell_loop

          d2 = dx**2 + dy**2 + dz**2
          if (d2 < r_Hill2) then ! particle is in Hill sphere
             masse_gaz(icell)    = 0.
             densite_gaz(icell) = 0.
             masse(icell) = 0.
             densite_pouss(:,icell) = 0.
             n_delete = n_delete + 1
          endif
       enddo cell_loop

       write(*,*) "Deleting", n_delete, "particles in Hill sphere of sink particle #", istar
    enddo

    return

  end subroutine delete_Hill_sphere

  !*********************************************************

  subroutine randomize_azimuth(n_points, x,y, vx,vy)

    use naleat, only : seed, stream, gtype
#include "sprng_f.h"


    integer, intent(in) :: n_points
    real(kind=dp), dimension(n_points), intent(inout) :: x, y, vx,vy

    integer, parameter :: nb_proc = 1
    integer :: i, id, istar

    real(kind=dp) :: cos_phi, sin_phi, phi, x_tmp, y_tmp

    particle_loop : do i=1, n_points
       ! We do not touch the sink particles
       do istar=1, n_etoiles
          if (i == etoile(istar)%icell) cycle particle_loop
       enddo

       call random_number(phi)
       phi = 2*pi*phi
       cos_phi = cos(phi) ; sin_phi = sin(phi)

       !-- position
       x_tmp = x(i) * cos_phi + y(i) * sin_phi
       y_tmp = -x(i) * sin_phi + y(i) * cos_phi
       x(i) = x_tmp ; y(i) = y_tmp


       !-- velocities
       x_tmp = vx(i) * cos_phi + vy(i) * sin_phi
       y_tmp = -vx(i) * sin_phi + vy(i) * cos_phi
       vx(i) = x_tmp ; vy(i) = y_tmp
    enddo particle_loop

    return

  end subroutine randomize_azimuth

  !*********************************************************

  subroutine read_ascii_SPH_file(iunit,filename,x,y,z,h,massgas,rhogas,rhodust,ndusttypes,n_SPH,ierr)

    integer,               intent(in) :: iunit
    character(len=*),      intent(in) :: filename
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,rhogas,massgas
    real(dp), intent(out), dimension(:,:), allocatable :: rhodust
    integer, intent(out) :: ndusttypes, n_SPH,ierr

    integer :: syst_status, alloc_status, ios, i
    character(len=512) :: cmd

    ierr = 0

    cmd = "wc -l "//trim(filename)//" > ntest.txt"
    call appel_syst(cmd,syst_status)
    open(unit=1,file="ntest.txt",status="old")
    read(1,*) n_SPH
    close(unit=1)
    ndusttypes =1

    write(*,*) "n_SPH read_test_ascii_file = ", n_SPH

    alloc_status = 0
    allocate(x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),massgas(n_SPH),rhogas(n_SPH),rhodust(ndusttypes,n_SPH), stat=alloc_status)
    if (alloc_status /=0) call error("Allocation error in phanton_2_mcfost")

    open(unit=1, file=filename, status='old', iostat=ios)
    do i=1, n_SPH
       read(1,*) x(i), y(i), z(i), h(i), massgas(i)
       rhogas(i) = massgas(i)
    enddo

    write(*,*) "MinMax=", minval(massgas), maxval(massgas)

    write(*,*) "Using stars from mcfost parameter file"

    return

  end subroutine read_ascii_SPH_file

  !*********************************************************

!  subroutine read_Judith_file(iunit,filename,x,y,z,h,massgas,rhogas,rhodust,ndusttypes,n_SPH,ierr)
!
!    !The 3D cube is 20 cells in  colatitude direction, 215 radially, and 680 azimuthally.
!    ! So when you read in the 1D  columns of arrays, you should maybe reform the arrays as [20,215,680] to get the 3D cube again
!
!    integer,               intent(in) :: iunit
!    character(len=*),      intent(in) :: filename
!    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,rhogas,massgas
!    real(dp), intent(out), dimension(:,:), allocatable :: rhodust
!    integer, intent(out) :: ndusttypes, n_SPH,ierr
!
!    integer :: syst_status, alloc_status, ios, i
!    character(len=512) :: cmd
!
!    ierr = 0
!
!    cmd = "wc -l "//trim(filename)//" > ntest.txt"
!    call appel_syst(cmd,syst_status)
!    open(unit=1,file="ntest.txt",status="old")
!    read(1,*) n_SPH
!    n_SPH = n_SPH - 1 ! removing 1 line of comments
!    close(unit=1)
!    ndusttypes =1
!
!    write(*,*) "n_SPH = ", n_SPH
!
!    alloc_status = 0
!    allocate(x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),massgas(n_SPH),rhogas(n_SPH),rhodust(ndusttypes,n_SPH), stat=alloc_status)
!    if (alloc_status /=0) call error("Allocation error in phanton_2_mcfost")
!
!    open(unit=1, file=filename, status='old', iostat=ios)
!    do i=1, n_SPH
!       read(1,*) x(i), y(i), z(i), rhogas(i), T(i), vx(i), vy(i), vz(i)
!    enddo
!
!    ! Correcting units : positions in au and velocities in m/s
!    x(:) = x(:) * cm_to_au ; y(:) = y(:) * cm_to_au ; z(:) = z(:) * cm_to_au
!    vx(:) = vx(:) * cm_to_m ; vy(:) = vy(:) * cm_to_m ; vz(:) = vz(:) * cm_to_m
!
!    write(*,*) "Using stars from mcfost parameter file"
!
!    return
!
!  end subroutine read_Judith_file

end module SPH2mcfost
