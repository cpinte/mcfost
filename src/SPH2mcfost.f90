module SPH2mcfost

  use parametres
  use constantes
  use utils
  use sort, only : find_kth_smallest_inplace
  use density, only : normalize_dust_density, reduce_density, read_Voronoi_fits_file
  use read_phantom, only : read_phantom_bin_files, read_phantom_hdf_files
  use sort, only : index_quicksort
  use stars, only : compute_stellar_parameters

  implicit none

  procedure(read_phantom_bin_files), pointer :: read_phantom_files => null()

contains

  subroutine setup_SPH2mcfost(extra_heating)

    use read_gadget2, only : read_gadget2_file
    use io_phantom_utils, only : get_error_text
    use utils, only : read_comments


    integer, parameter :: iunit = 1

    real(dp), allocatable, dimension(:) :: x,y,z,h,vx,vy,vz,rho,massgas,SPH_grainsizes,T_gas
    real(dp), allocatable, dimension(:) :: vturb,mass_ne_on_massgas,atomic_mask
    integer,  allocatable, dimension(:) :: particle_id
    real(dp), allocatable, dimension(:,:) :: rhodust, massdust, dust_moments
    real, allocatable, dimension(:) :: extra_heating
    integer, allocatable, dimension(:) :: mask ! size == np, not n_SPH, index is original SPH id (update, I think it is id)

    integer, dimension(:), allocatable :: is_ghost
    real(dp), dimension(6) :: SPH_limits
    real :: factor
    integer :: ndusttypes, ierr, i, ilen, n_SPH
    logical :: check_previous_tesselation, ldust_moments
    real(dp) :: mass_per_H

    ldust_moments = .false.
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
            SPH_grainsizes,mask,n_SPH,ldust_moments,dust_moments,mass_per_H,ierr)

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
       write(*,*) "Reading Gadget-2 density file: "//trim(density_files(1))
       call read_gadget2_file(iunit,density_files(1), x,y,z,h,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    else if (ldensity_file) then
       call read_Voronoi_fits_file(density_files(1), x,y,z,h,vx,vy,vz,particle_id,massgas,n_SPH)
       ldust_moments = .false.
       ndusttypes=0
    else
       call error("Unknown SPH structure.")
    endif
    write(*,*) "Done"

    if (lignore_dust) then
       ndusttypes = 0
       if (allocated(rhodust)) deallocate(rhodust)
       if (allocated(massdust)) deallocate(massdust)
    endif

    if (.not.lfix_star) call compute_stellar_parameters()

    ! Model limits
    call read_SPH_limits_file(limits_file, SPH_limits)

    ! Voronoi tesselation
    check_previous_tesselation = (.not. lrandomize_Voronoi)
    call SPH_to_Voronoi(n_SPH, ndusttypes, particle_id, x,y,z,h, vx,vy,vz, &
         T_gas, massgas,massdust,rho,rhodust,SPH_grainsizes, SPH_limits, check_previous_tesselation, &
         is_ghost, ldust_moments, dust_moments, mass_per_H, mask=mask)

    deallocate(x,y,z,h)
    if (allocated(vx)) deallocate(vx,vy,vz)
    if (allocated(rhodust)) deallocate(rhodust,massdust)

    ! setup needed for Atomic line transfer
    if (lemission_atom) then
       call hydro_to_Voronoi_atomic(n_SPH,T_gas,vturb,massgas,mass_ne_on_massgas,atomic_mask)
    endif
    deallocate(massgas)
    if (allocated(rho)) deallocate(rho)
    if (allocated(vturb)) deallocate(vturb)
    if (allocated(mass_ne_on_massgas)) deallocate (mass_ne_on_massgas)
    if (allocated(atomic_mask)) deallocate(atomic_mask)

    ! Deleting particles/cells in masked areas (Hill sphere, etc)
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

  subroutine print_moments(ki,a0,rhoi)

    real(dp), intent(in) :: ki(0:3),a0,rhoi

    print "(2x,a,1pg0.4,a)",         'a) number density of dust = ',ki(0),' cm^-3'
    print "(2x,a,1pg0.4,a,1pg0.2,a)",'b) average grain radius   = ',ki(1)/ki(0),' times a0, or ',ki(1)/ki(0)*a0, ' '
    print "(2x,a,1pg0.4,a,1pg0.2,a)",'c) average grain area     = ',4.*pi*a0**2*ki(2)/ki(0),' '
    print "(2x,a,1pg0.4,a,1pg0.2,a)",'d) average particle size  = ',ki(3)/ki(0),' monomers'
    print "(2x,a,1pg0.4,a,1pg0.2,a)",'e) opacity at Td=100      = ',ki(3)*pi*a0**3/rhoi*6.7*100. / mum_to_cm**3  ! Eq (42-43) of Siess et al. 2022

  end subroutine print_moments

  !*********************************************************

  subroutine SPH_to_Voronoi(n_SPH, ndusttypes, particle_id, x,y,z,h, vx,vy,vz, T_gas, massgas,massdust,rho,rhodust,&
       SPH_grainsizes, SPH_limits, check_previous_tesselation, is_ghost, ldust_moments, dust_moments, mass_per_H, mask)

    ! ************************************************************************************ !
    ! n_sph : number of points in the input model
    ! ndusttypes : number of dust species in the input model
    ! particle_id : index of the particle / cell centre in the input model
    ! x,y,z : coordinates of the particle / cell
    ! h : smoothing length to cut elongated cells
    ! vx,vy,vz : vector components of the velocity fields (cartesian) in the input model
    ! massgas : total mass of the gas
    ! massdust : total mass of the dust
    ! rho : gas density
    ! rhodust : dust density
    ! sph_grainsizes : size of grains for sph models
    ! sph_limits : limit of the input model box
    ! check_previous_tesselation :
    ! mask : integer array, 1 if a SPH particle will be made transparent, 2 if deleted before tesselation
    ! ************************************************************************************ !
    use Voronoi_grid
    use density, only : densite_gaz, masse_gaz, densite_pouss, masse
    use grains, only : n_grains_tot, M_grain
    use disk_physics, only : compute_othin_sublimation_radius
    use mem
    use reconstruct_from_moments
    use, intrinsic :: ieee_arithmetic

    integer, intent(in) :: n_SPH, ndusttypes
    real(dp), dimension(:), allocatable, intent(inout) :: x,y,z,h,massgas!,rho, !move rho to allocatable, assuming not always allocated
    real(dp), dimension(:), allocatable, intent(inout) :: rho
    real(dp), dimension(:), allocatable, intent(inout) :: vx,vy,vz ! dimension n_SPH or 0
    real(dp), dimension(:), allocatable, intent(in) :: T_gas
    integer, dimension(n_SPH), intent(in) :: particle_id
    real(dp), dimension(:,:), allocatable, intent(inout) :: rhodust, massdust, dust_moments ! ndusttypes,n_SPH
    real(dp), dimension(:), allocatable, intent(in) :: SPH_grainsizes ! ndusttypes
    real(dp), dimension(6), intent(in) :: SPH_limits
    logical, intent(in) :: check_previous_tesselation
    integer, dimension(:), allocatable, intent(in), optional :: mask
    logical, intent(in) :: ldust_moments
    real(dp), intent(in) :: mass_per_H

    integer, dimension(:), allocatable, intent(out) :: is_ghost

    logical :: use_single_grain

    real, allocatable, dimension(:) :: a_SPH, log_a_SPH, rho_dust
    real(dp), allocatable, dimension(:) :: gsize, grainsize_f, dN_ds, N_monomers, rho_monomers
    real(dp), dimension(4) :: lambsol, lambguess

    real(dp) :: mass, somme, Mtot, Mtot_dust, facteur, a, mass_factor
    real :: f, limit_threshold, density_factor
    integer :: icell, l, k, iSPH, n_force_empty, i, id_n, ierr, N_pb

    real(dp), dimension(6) :: limits

    real(dp), parameter :: a0 = 1.28e-4 ! microns. Siess et al 2022

    real(dp), dimension(4) :: ki, err
    real(dp) ::rhoi, norm, mdust, mdust_tot, factor


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

    if (ndusttypes >= 1) then
       write(*,*) "Found", n_SPH, " hydro sites with ", ndusttypes, "dust grains."
    else
       write(*,*) "Found", n_SPH, " hydro sites."
    endif

    if (abs(maxval(SPH_limits)) < tiny_real) then
       write(*,*) "Selecting spatial range which contains"
       write(*,*) SPH_keep_particles*100, "% of particles in each dimension"

       k = max(int(limit_threshold * n_SPH),1)
       limits(1) = find_kth_smallest_inplace(k,real(x))
       limits(3) = find_kth_smallest_inplace(k,real(y))
       limits(5) = find_kth_smallest_inplace(k,real(z))

       k = int((1.0-limit_threshold) * n_SPH)
       limits(2) = find_kth_smallest_inplace(k,real(x))
       limits(4) = find_kth_smallest_inplace(k,real(y))
       limits(6) = find_kth_smallest_inplace(k,real(z))
    else
       limits(:) = SPH_limits(:)
    endif

    if (ldelete_outside_rSPH) then
       limits(1) = max(-rsph_max,limits(1))
       limits(3) = max(-rsph_max,limits(3))
       limits(5) = max(-rsph_max,limits(5))

       limits(2) = min(rsph_max,limits(2))
       limits(4) = min(rsph_max,limits(4))
       limits(6) = min(rsph_max,limits(6))
    endif

    write(*,*) "# Model limits :"
    write(*,*) "x =", limits(1), limits(2)
    write(*,*) "y =", limits(3), limits(4)
    write(*,*) "z =", limits(5), limits(6)

    call test_duplicate_particles(n_SPH, particle_id, x,y,z, massgas,massdust,rho,rhodust, is_ghost)

    !*******************************
    ! Voronoi tesselation
    !*******************************
    call Voronoi_tesselation(n_SPH, particle_id, x,y,z,h,vx,vy,vz, is_ghost, mask, limits, check_previous_tesselation)
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

    mass = 0.
    do icell=1,n_cells
       mass = mass + masse_gaz(icell)
    enddo !icell
    mass =  mass * g_to_Msun

    if (lforce_Mgas) then ! Todo : testing for now, use a routine to avoid repeting code
       write(*,*) "Forcing gas mass to value in parameter file rather"
       ! Normalisation
       if (mass > 0.0) then ! pour le cas ou gas_to_dust = 0.
          facteur = disk_zone(1)%diskmass * disk_zone(1)%gas_to_dust / mass

          ! Somme sur les zones pour densite finale
          do icell=1,n_cells
             densite_gaz(icell) = densite_gaz(icell) * facteur
             masse_gaz(icell) = masse_gaz(icell) * facteur
          enddo ! icell
       else
          call error('Gas mass is 0 in hydero file')
       endif
    endif

    ! Tableau de densite et masse de poussiere
    ! interpolation en taille
    if (ldust_moments) then
       write(*,*) "Reconstructing grain size distribution from moments ..."

       mass = 0.0
       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) mass = mass +  massgas(iSPH) * dust_moments(3,iSPH) * 12.*amu/mass_per_H
       enddo
       write(*,*) "Dust mass in hydro model is ", real(mass), "Msun"

       lvariable_dust = .true.
       allocate(grainsize_f(n_grains_tot),gsize(n_grains_tot),dN_ds(n_grains_tot),&
            N_monomers(n_grains_tot),rho_monomers(n_grains_tot))

       N_monomers(:) = (r_grain(:)/a0)**3 ! number of monomers
       dN_ds(:) = 3*r_grain(:)**2/a0**3

       mass_factor = 12.*amu/mass_per_H

       N_pb = 0
       !$omp parallel &
       !$omp default(none) &
       !$omp private(icell,iSPH,use_single_grain,rhoi,ki,err,ierr,lambsol,rho_monomers,a,i,norm,mdust) &
       !$omp shared(n_cells,Voronoi,densite_gaz,n_grains_tot,r_grain,mass_factor,dN_ds,N_monomers) &
       !$omp shared(densite_pouss,dust_moments,masse_gaz,volume) &
       !$omp reduction(+:N_pb)
       !$omp do schedule(dynamic,1)
       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) then
             use_single_grain = .false.

             !rhoi =  densite_gaz(icell) * cm_to_m**3 ! check that it is in cgs
             ki = dust_moments(:,iSPH) !/ mass_per_H

             call reconstruct_gamma_dist(ki,lambsol,err,ierr)
             if (ierr/=1) N_pb = N_pb+1 ! 1 means no error!!!

             do i=1,n_grains_tot
                rho_monomers(i) = gamma_func_from_moments(N_monomers(i),ki(1:2),lambsol(1:2))
                if (.not. ieee_is_finite(rho_monomers(i))) then
                   use_single_grain = .true.
                   exit
                endif
             enddo

             densite_pouss(:,icell) = rho_monomers(:) * dN_ds(:)

             ! Simple approximation : we assume 1 single grain size
             if (use_single_grain) then
                densite_pouss(:,icell) = 0._dp
                if (dust_moments(1,iSPH) > tiny_dp) then
                   a = a0 * dust_moments(2,iSPH)/dust_moments(1,iSPH)
                   i = locate(1.0_dp*r_grain(:),a)
                   densite_pouss(i,icell) = 1.0
                endif
             endif
          else ! iSPH == 0, star
             densite_pouss(:,icell) = 0._dp
          endif
       enddo ! icell
       !$omp end do
       !$omp end parallel

       write(*,*) "Done"
       if (N_pb > 0) write(*,*) "N cells with pb", N_pb, "/", n_cells

       ! We normalise to the total dust mass in the dump
       mdust_tot = 0.0_dp
       do icell=1,n_cells
          iSPH = Voronoi(icell)%id

          mass = 0.0_dp
          do l=1,n_grains_tot
             mass=mass + (densite_pouss(l,icell) *1.0_dp) * M_grain(l)
          enddo !l
          mass = mass * volume(icell)

          if (iSPH > 0) mdust =  massgas(iSPH) * dust_moments(3,iSPH) * 12.*amu/mass_per_H
          mdust_tot = mdust_tot + mdust

          if (mass > tiny_dp) then
             factor = mdust/ mass
             densite_pouss(:,icell) = densite_pouss(:,icell) * factor
          endif
       enddo !icell
       call normalize_dust_density(mdust_tot) ! this should only calculates the array masse

    elseif (ndusttypes >= 1) then
       lvariable_dust = .true.

       ! mcfost adds an extra grain size follwing the gas
       ! this is required if ndusttypes == 1, but I do it in any case
       allocate(a_SPH(ndusttypes+1),log_a_SPH(ndusttypes+1),rho_dust(ndusttypes+1))

       write(*,*) "Found the following grain sizes in hydro calculations:"
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

       if (lforce_SPH_amin) a_SPH(1) = SPH_amin
       if (lforce_SPH_amax) then
          if (ndusttypes == 1) then
             a_SPH(2) = SPH_amax
          else
             call error("Can only froce SPH grain size is ndusttypes == 1")
          endif
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
    ! Mask
    !*************************
    if (present(mask)) then
       if (allocated(mask)) then
          do icell=1,n_cells
             iSPH = Voronoi(icell)%id
             if (iSPH > 0) then
                Voronoi(icell)%masked = mask(iSPH)
             else
                Voronoi(icell)%masked = 0
             endif
          enddo
       else
          do icell=1,n_cells
             Voronoi(icell)%masked = 0
          enddo
       endif
    else
       do icell=1,n_cells
          Voronoi(icell)%masked = 0
       enddo
    endif

    ! We eventually reduce density to avoid artefacts: superseeded by cell cutting
    !massgas is also modified with the density factor for compatiblity with atomictransfer
    !since massgas is inout.
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
             massgas(Voronoi(icell)%id) = massgas(Voronoi(icell)%id) * density_factor
             cycle cell_loop
          endif

          ! We reduce the density on cells that are touching a wall
          do i=Voronoi(icell)%first_neighbour, Voronoi(icell)%last_neighbour
             id_n = neighbours_list(i) ! id du voisin
             if (id_n < 0) then
                n_force_empty = n_force_empty + 1
                call reduce_density(icell, density_factor)
                massgas(Voronoi(icell)%id) = massgas(Voronoi(icell)%id) * density_factor
                cycle cell_loop
             endif
          enddo
       enddo cell_loop
       write(*,*) "Density was reduced by", density_factor, "in", n_force_empty, "cells surrounding the model, ie",&
            (1.0*n_force_empty)/n_cells * 100, "% of cells"
    endif

    search_not_empty : do k=1,n_grains_tot
       do icell=1, n_cells
          if (densite_pouss(k,icell) > 0.0_sp) then
             icell_not_empty = icell
             exit search_not_empty
          endif
       enddo !icell
    enddo search_not_empty

    if (ndusttypes >= 1) deallocate(a_SPH,log_a_SPH,rho_dust)

    write(*,*) 'Total  gas mass in model :',  real(sum(masse_gaz) * g_to_Msun),' Msun'
    write(*,*) 'Total dust mass in model :', real(sum(masse) * g_to_Msun),' Msun'

    return

  end subroutine SPH_to_Voronoi

  !****************************************************************************

  subroutine Hydro_to_Voronoi_atomic(n_SPH,T_tmp,vt_tmp,mass_gas,mass_ne_on_massgas,mask)
    ! copy additional parameters from hydro code needed for atomic line transfer
    ! ************************************************************************************ !
    ! n_sph : number of points in the input model
    ! particle_id : index of the particle / cell centre in the input model
    ! T_tmp : temperature of the particle/ cell
    ! vt_tmp : turbulent velocity in the particle / cell
    ! massgas : total mass of the gas
    ! rho : gas density
    ! mask : -1 means skip, 0 means transparent, 1 means compute atomic transfer
    ! ************************************************************************************ !
    use parametres
    use constantes,   only : masseH
    use Voronoi_grid, only : Voronoi, volume
    use disk_physics, only : compute_othin_sublimation_radius
    use mem
    use elements_type, only : wght_per_H, read_abundance
    use grid, only : alloc_atomrt_grid, nHtot, ne, v_char, lmagnetized, vturb, T, icompute_atomRT, lcalc_ne, &
         check_for_zero_electronic_density

    integer, intent(in) :: n_SPH
    real(dp), dimension(n_SPH), intent(in) :: T_tmp,mass_gas
    real(dp), dimension(:), allocatable, intent(in) :: mass_ne_on_massgas ! not always allocated
    real(dp), dimension(:), allocatable, intent(in) :: vt_tmp,mask

    integer :: icell,voroindex
    real(kind=dp) :: rho_to_nH, vxmax, vxmin, vymax, vymin, vzmax, vzmin, vmax
    real, parameter :: Lextent = 1.01

    !*******************************
    ! Accomodating model
    !*******************************
    !alloc space for all physical quantities.
    !Velocity field arrays not allocated because lvoronoi is true
    !-> fills element abundances structures for elements
    call alloc_atomrt_grid
    call read_abundance
    rho_to_nH = 1d3 / masseH / wght_per_H !convert from density to number density nHtot

    Vxmax = 0
    Vymax = 0
    Vzmax = 0
    Vxmin = 1d100
    Vymin = 1d100
    Vzmin = 1d00
    Vmax = 0
    icompute_atomRT(:) = 0
    do icell=1,n_cells
       voroindex = Voronoi(icell)%id
       if (voroindex > 0) then
          nHtot(icell)  = Msun_to_kg * rho_to_nH * mass_gas(voroindex) /  (volume(icell) * AU3_to_m3)

          !store the ratio of mass electron on mass hydrogen
          !to do the tesselation only on gas particle
          !-> if ne is not present in the model, it is simply 0.
          if (allocated(mass_ne_on_massgas)) then
             ne(icell) = Msun_to_kg * mass_ne_on_massgas(voroindex) * mass_gas(voroindex) / (volume(icell) * AU3_to_m3)
          endif

          T(icell) = T_tmp(voroindex)

          if (allocated(vt_tmp)) vturb(icell) = vt_tmp(voroindex)

          if (allocated(mask)) then
             icompute_atomRT(icell) = mask(voroindex)
          else
             icompute_atomRT(icell) = 1
          endif
          vxmax = max(vxmax, abs(Voronoi(icell)%vxyz(1)))
          vxmin = min(vxmin, max(abs(Voronoi(icell)%vxyz(1)),0.0))
          vymax = max(vymax, abs(Voronoi(icell)%vxyz(2)))
          vymin = min(vymin,  max(abs(Voronoi(icell)%vxyz(2)),0.0))
          vzmax = max(vzmax, abs(Voronoi(icell)%vxyz(3)))
          vzmin = min(vzmin,  max(abs(Voronoi(icell)%vxyz(3)),0.0))

          Vmax = max( vmax, sqrt(Voronoi(icell)%vxyz(1)**2 + Voronoi(icell)%vxyz(2)**2 + Voronoi(icell)%vxyz(3)**2) )
       endif
    end do

    v_char = Vmax * Lextent

    write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT>0)), " density zones"
    write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT==0)), " transparent zones"
    write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT<0)), " dark zones"

    call check_for_zero_electronic_density()

    write(*,*) "Maximum/minimum velocities in the model (km/s):"
    write(*,*) "|Vx|", vxmax*1d-3, vxmin*1d-3
    write(*,*) "|Vy|", vymax*1d-3, vymin*1d-3
    write(*,*) "|Vz|", vzmax*1d-3, vzmin*1d-3

    write(*,*) "Typical line extent due to V fields (km/s):"
    write(*,*) v_char/1d3

    write(*,*) "Maximum/minimum turbulent velocity (km/s):"
    write(*,*) maxval(vturb)/1d3, minval(vturb, mask=icompute_atomRT>0)/1d3

    write(*,*) "Maximum/minimum Temperature in the model (K):"
    write(*,*) maxval(T), minval(T,mask=icompute_atomRT>0)
    write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
    write(*,*) maxval(nHtot), minval(nHtot,mask=icompute_atomRT>0)
    if (.not.lcalc_ne) then
       write(*,*) "Maximum/minimum ne density in the model (m^-3):"
       write(*,*) maxval(ne), minval(ne,mask=icompute_atomRT>0)
    endif

    return

  end subroutine hydro_to_Voronoi_atomic

  !*********************************************************

  subroutine test_duplicate_particles(n_SPH, particle_id, x,y,z, massgas,massdust,rho,rhodust, is_ghost)
    ! Filtering particles at the same position
    ! They are defined as ghost of the main one. The main gets the total mass and denity.
    ! We keep the main particle_id to be able to give the ghost particle the same temperature as the main particle

    integer, intent(in) :: n_SPH
    integer, dimension(n_SPH), intent(in) :: particle_id
    real(kind=dp), dimension(n_SPH), intent(inout) :: x, y, z, massgas
    real(kind=dp), dimension(:), allocatable, intent(inout) :: rho
    real(dp), dimension(:,:), allocatable, intent(inout) :: rhodust, massdust

    integer, dimension(:), allocatable, intent(out) :: is_ghost

    real, parameter :: prec = 1e-14

    real(kind=dp), dimension(:), allocatable :: x2
    integer, dimension(:), allocatable :: order
    real :: dr2
    integer :: i, j, ii, jj, nkill, alloc_status

    alloc_status=0
    allocate(x2(n_SPH), order(n_SPH), is_ghost(n_SPH), stat=alloc_status)
    if (alloc_status /=0) call error("Allocation error test_duplicate_particles")

    x2(:) = x(:)**2 + y(:)**2 + z(:)**2
    order = index_quicksort(x2)

    is_ghost(:) = 0
    nkill = 0
    do i=1, n_SPH
       ii = order(i)
       loop2 : do j=i+1, n_SPH
          jj = order(j)
          dr2 = (x(ii)-x(jj))**2 + (y(ii)-y(jj))**2 + (z(ii)-z(jj))**2
          if (dr2 < prec * x2(ii)) then
             is_ghost(jj) = particle_id(ii)
             nkill = nkill+1

             ! Adding ghost particle to main one
             massgas(ii) = massgas(ii) + massgas(jj)
             if (.not.lignore_dust) then
               rho(ii) = rho(ii) + rho(jj)
               rhodust(:,ii) = rhodust(:,ii) + rhodust(:,jj)
               massdust(:,ii) = massdust(:,ii) + massdust(:,jj)
             endif
          else
             exit loop2
          endif
       enddo loop2 ! j
    enddo ! i
    deallocate(order,x2)

    if (nkill > 0) then
       write(*,*)
       write(*,*) nkill, "hydro sites were flagged as ghosts and merged"
       write(*,*)
    endif

    return

  end subroutine test_duplicate_particles

  !*********************************************************

  subroutine delete_masked_particles()

    use Voronoi_grid
    use density, only : densite_gaz, masse_gaz, densite_pouss, masse

    integer :: icell, k

    k=0
    do icell=1, n_cells
       if (Voronoi(icell)%masked == 1) then
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
