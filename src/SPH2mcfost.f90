module SPH2mcfost

  use parametres
  use constantes
  use utils
  use density, only : normalize_dust_density, reduce_density

  implicit none

contains

  subroutine setup_SPH2mcfost(SPH_file,SPH_limits_file, n_SPH, extra_heating)

    use read_phantom, only : read_phantom_file
    use read_gadget2, only : read_gadget2_file
    use dump_utils, only : get_error_text
    use utils, only : read_comments
    use prop_star, only : etoile
    use disk, only : lscale_SPH, scale_SPH, lfix_star


    character(len=512), intent(in) :: SPH_file, SPH_limits_file
    integer, intent(out) :: n_SPH

    integer, parameter :: iunit = 1

    real(dp), allocatable, dimension(:) :: x,y,z,h,vx,vy,vz,rho,massgas,SPH_grainsizes
    integer,  allocatable, dimension(:) :: particle_id
    real(dp), allocatable, dimension(:,:) :: rhodust, massdust
    real, allocatable, dimension(:) :: extra_heating

    real(dp), dimension(6) :: SPH_limits
    integer :: ndusttypes, ierr
    character(len=100) :: line_buffer

    if (lphantom_file) then
       write(*,*) "Performing phantom2mcfost setup"
       write(*,*) "Reading phantom density file: "//trim(SPH_file)
       call read_phantom_file(iunit,SPH_file,x,y,z,h,vx,vy,vz, &
            particle_id,massgas,massdust,rho,rhodust,extra_heating,ndusttypes,SPH_grainsizes,n_SPH,ierr)
       ! Todo : extra heating must be passed to mcfost
       if (ierr /=0) then
          write(*,*) "Error code =", ierr,  get_error_text(ierr)
          call exit(1)
       endif
    else if (lgadget2_file) then
       write(*,*) "Performing gadget2mcfost setup"
       write(*,*) "Reading Gadget-2 density file: "//trim(SPH_file)
       call read_gadget2_file(iunit,SPH_file, x,y,z,h,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    else if (lascii_SPH_file) then
       write(*,*) "Performing SPH2mcfost setup"
       write(*,*) "Reading SPH density file: "//trim(SPH_file)
       call read_ascii_SPH_file(iunit,SPH_file, x,y,z,h,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    endif
    write(*,*) "Done"

    if ((.not.lfix_star).and.(lphantom_file .or. lgadget2_file)) call compute_stellar_parameters()

    if (lscale_SPH) then
       write(*,*) "**************************************************"
       write(*,*) "WARNING : rescaling SPH simulation by:", scale_SPH
       write(*,*) "**************************************************"
       x = x * scale_SPH ; y = y * scale_SPH ; z = z * scale_SPH
       etoile(:)%x = etoile(:)%x * scale_SPH ; etoile(:)%y = etoile(:)%y * scale_SPH ; etoile(:)%z = etoile(:)%z * scale_SPH
    endif

    ! Model limits
    call read_SPH_limits_file(SPH_limits_file, SPH_limits)

    ! Voronoi tesselation
    call SPH_to_Voronoi(n_SPH, ndusttypes, x,y,z,h, vx,vy,vz, massgas,massdust,rho,rhodust,SPH_grainsizes, SPH_limits, .true.)

    deallocate(massgas,rho)
    if (ndusttypes > 0) then
       deallocate(rhodust,massdust)
    endif

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

  subroutine SPH_to_Voronoi(n_SPH, ndusttypes, x,y,z,h, vx,vy,vz, massgas,massdust,rho,rhodust,SPH_grainsizes, &
       SPH_limits, check_previous_tesselation)

    use Voronoi_grid
    use opacity, only : densite_pouss, masse
    use molecular_emission, only : densite_gaz, masse_gaz
    use grains, only : n_grains_tot, M_grain
    use disk_physics, only : compute_othin_sublimation_radius
    use mem
    use disk, only : SPH_keep_particles, disk_zone

    integer, intent(in) :: n_SPH, ndusttypes
    real(dp), dimension(n_SPH), intent(in) :: x,y,z,h,rho,massgas
    real(dp), dimension(:), allocatable, intent(in) :: vx,vy,vz ! dimension n_SPH or 0
    real(dp), dimension(ndusttypes,n_SPH), intent(in) :: rhodust, massdust
    real(dp), dimension(ndusttypes), intent(in) :: SPH_grainsizes
    real(dp), dimension(6), intent(in) :: SPH_limits
    logical, intent(in) :: check_previous_tesselation

    real, parameter :: density_factor = 1 !e-6
    logical :: lwrite_ASCII = .true. ! produce an ASCII file for yorick

    real, allocatable, dimension(:) :: a_SPH, log_a_SPH, rho_dust
    real(dp) :: mass, somme, Mtot, Mtot_dust
    real :: f, limit_threshold
    integer :: icell, l, k, iSPH, n_force_empty, i, id_n

    real(dp), dimension(6) :: limits

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
    call Voronoi_tesselation(n_SPH, x,y,z,h, limits, check_previous_tesselation)
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
          a_SPH(l+1) = SPH_grainsizes(l)
          write(*,*) real(SPH_grainsizes(l)), "microns"
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
          a_SPH(2) = 1. ;
          a_SPH(2) = 1000. ;
          write(*,*) "WARNING: SPH dust grain size  not found"
          write(*,*) "forcing small & big grains to be 1mum & 1mm"
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

       do icell=1,n_cells
          masse(icell) = 0.
          do k=1,n_grains_tot
             densite_pouss(k,icell) = densite_gaz(icell) * nbre_grains(k)
             masse(icell) = masse(icell) + densite_pouss(k,icell) * M_grain(k) * volume(icell)
          enddo
       enddo
       masse(:) = masse(:) * AU3_to_cm3
       f = 1./disk_zone(1)%gas_to_dust * sum(masse_gaz)/sum(masse)
       densite_pouss(:,:) = densite_pouss(:,:) * f
       masse(:) = masse(:) * f
    endif ! ndusttypes == 0


    !*************************
    ! Velocities
    !*************************
    if (lemission_mol) then
       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) then
             Voronoi(icell)%vxyz(1) = vx(iSPH)
             Voronoi(icell)%vxyz(2) = vy(iSPH)
             Voronoi(icell)%vxyz(3) = vz(iSPH)
          endif
       enddo
    endif

    ! We eventually reduce density to avoid artefacts: superseeded by cell cutting
    if (density_factor < 1.-1e-6) then
       ! Removing cells at the "surface" of the SPH model:
       ! density is reduced so that they do not appear in images or cast artificial shadows,
       ! but we can still compute a temperature (forcing them to be optically thin)
       n_force_empty = 0.0
       cell_loop : do icell=1,n_cells
          ! We reduce the density on cells that are very elongated
          if (Voronoi(icell)%delta_edge > 3 * Voronoi(icell)%h) then
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

    use prop_star

    integer :: i

    character(len=512) :: isochrone_file, filename
    character(len=100) :: line_buffer

    integer, parameter :: nSpT = 29
    character(len=2), dimension(nSpT) :: SpT
    real :: L, R, T, M, minM, maxM
    real, dimension(nSpT) :: logL, logR, logTeff, logM

    if (n_etoiles < 1) return

    isochrone_file = "Siess/isochrone_"//trim(system_age)//".txt"

    write(*,*) ""
    write(*,*) "Reading isochrone file: "//trim(isochrone_file)
    filename = trim(mcfost_utils)//"/Isochrones/"//trim(isochrone_file)
    open(unit=1,file=filename,status="old")
    do i=1,3
       read(1,*) line_buffer
    enddo
    minM = 1.e30 ; maxM = 0 ;
    do i=1, nSpT
       read(1,*) SpT(i), L, r, T, M
       logL(i) = log(L) ; logR(i) = log(r) ; logTeff(i) = log(T) ; logM(i) = log(M)
       if (M < minM) minM = M
       if (M > maxM) maxM = M
    enddo
    close(unit=1)

    ! interpoler L et T, les fonctions sont plus smooth
    write(*,*) "New stellar parameters:"
    do i=1, n_etoiles
       if ((etoile(i)%M < minM) .or. (etoile(i)%M > maxM))  then
          write(*,*) "   *** WARNING : stellar object mass not in isochrone range"
          write(*,*) "   *** object #", i, "M=", etoile(i)%M, "Msun"
       endif
       etoile(i)%T = exp(interp(logTeff, logM, log(etoile(i)%M)))
       etoile(i)%r = exp(interp(logR, logM, log(etoile(i)%M)))

       ! Pas de fUV et pas de spectre stellaire pour le moment
       etoile(i)%fUV = 0.0 ; etoile(i)%slope_UV = 0.0 ;
       etoile(i)%lb_body = .true. ; etoile(i)%spectre = "None"

       write(*,*) "Star #",i,"  Teff=", etoile(i)%T, "K, r=", etoile(i)%r, "Rsun"
    enddo
    write(*,*) ""

    ! Passage rayon en AU
    etoile(:)%r = etoile(:)%r * Rsun_to_AU

    return

  end subroutine compute_stellar_parameters

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

end module SPH2mcfost
