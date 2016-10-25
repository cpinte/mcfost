module SPH2mcfost

  use parametres
  use constantes
  use utils

  implicit none

contains

  subroutine setup_SPH2mcfost(SPH_file,SPH_limits_file, n_SPH)

    use read_phantom, only : read_phantom_file, read_phantom_input_file
    use read_gadget2, only : read_gadget2_file
    use dump_utils, only : get_error_text

    character(len=512), intent(in) :: SPH_file, SPH_limits_file
    integer, intent(out) :: n_SPH

    integer, parameter :: iunit = 1

    real(dp), allocatable, dimension(:) :: x,y,z,rho,massgas,grainsize
    real(dp), allocatable, dimension(:,:) :: rhodust, massdust

    real(dp), dimension(6) :: SPH_limits
    integer :: ndusttypes, ierr, ios
    character(len=100) :: line_buffer

    if (lphantom_file) then
       write(*,*) "Performing phantom2mcfost setup"
       write(*,*) "Reading phantom density file: "//trim(SPH_file)
       call read_phantom_file(iunit,SPH_file, x,y,z,massgas,massdust,rho,rhodust,ndusttypes,grainsize,n_SPH,ierr)
       if (ierr /=0) then
          write(*,*) "Error code =", ierr,  get_error_text(ierr)
          stop
       endif
    else if (lgadget2_file) then
       write(*,*) "Performing gadget2mcfost setup"
       write(*,*) "Reading Gadget-2 density file: "//trim(SPH_file)
       call read_gadget2_file(iunit,SPH_file, x,y,z,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    else if (lascii_SPH_file) then
       write(*,*) "Performing SPH2mcfost setup"
       write(*,*) "Reading SPH density file: "//trim(SPH_file)
       call read_ascii_SPH_file(iunit,SPH_file, x,y,z,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    endif
    write(*,*) "Done"

    if (lphantom_file .or. lgadget2_file) call compute_stellar_parameters()

    !*******************************
    ! Model limits
    !*******************************
    write(*,*) " "
    if (llimits_file) then
       write(*,*) "Reading limits file: "//trim(SPH_limits_file)
       open(unit=1, file=SPH_limits_file, status='old', iostat=ios)
       if (ios/=0) then
          write(*,*) "ERROR : cannot open "//trim(SPH_limits_file)
          write(*,*) "Exiting"
          stop
       endif
       read(1,*) line_buffer
       read(1,*) SPH_limits(1), SPH_limits(3), SPH_limits(5)
       read(1,*) SPH_limits(2), SPH_limits(4), SPH_limits(6)
       close(unit=1)
    else
       SPH_limits(:) = 0
    endif

    call SPH_to_Voronoi(n_SPH, ndusttypes, x,y,z,massgas,massdust,rho,rhodust,grainsize, SPH_limits, .true.)

    deallocate(massgas,rho)
    if (ndusttypes > 0) then
       deallocate(rhodust,massdust)
    endif

    return

  end subroutine setup_SPH2mcfost

  !*********************************************************

  subroutine SPH_to_Voronoi(n_SPH, ndusttypes, x,y,z,massgas,massdust,rho,rhodust,grainsize, &
       SPH_limits, check_previous_tesselation)

    use Voronoi_grid
    use opacity, only : densite_pouss, masse
    use molecular_emission, only : densite_gaz, masse_gaz
    use grains, only : n_grains_tot, M_grain
    use disk_physics, only : compute_othin_sublimation_radius
    use mem

    integer, intent(in) :: n_SPH, ndusttypes
    real(dp), dimension(n_SPH), intent(in) :: x,y,z,rho,massgas
    real(dp), dimension(ndusttypes,n_SPH), intent(in) :: rhodust, massdust
    real(dp), dimension(ndusttypes), intent(in) :: grainsize
    real(dp), dimension(6), intent(in) :: SPH_limits
    logical, intent(in) :: check_previous_tesselation

    real, parameter :: limit_threshold = 0.01

    logical :: lwrite_ASCII = .false. ! produce an ASCII file for yorick

    real, allocatable, dimension(:) :: a_SPH, log_a_SPH, rho_dust
    real(dp) :: mass, somme, Mtot, Mtot_dust, dust_to_gas
    real :: graindens, f
    integer :: ierr, icell, l, k, iSPH

    real(dp), dimension(6) :: limits


    icell_ref = 1

    write(*,*) "# Farthest particules :"
    write(*,*) "x =", minval(x), maxval(x)
    write(*,*) "y =", minval(y), maxval(y)
    write(*,*) "z =", minval(z), maxval(z)

    !if (ndusttypes==1) then
    !   call read_phantom_input_file("hltau.in",iunit,grainsize,graindens,ierr)
    !   write(*,*) grainsize,graindens
    !endif
    write(*,*) "Found", n_SPH, " SPH particles with ", ndusttypes, "dust grains"

    if (lwrite_ASCII) then
       ! Write the file for the grid version of mcfost
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
       write(*,*) 1.0-2*limit_threshold, "particles in each dimension"

       k = int(limit_threshold * n_SPH)
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
    call Voronoi_tesselation(n_SPH, x,y,z, limits, check_previous_tesselation)
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

    do icell=1,n_cells
       iSPH = Voronoi(icell)%id
       if (iSPH > 0) then
          masse_gaz(icell)    = massgas(iSPH) /  g_to_Msun
          densite_gaz(icell)  = masse_gaz(icell) /  (AU3_to_cm3 * masse_mol_gaz * volume(icell))
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
       do l=1, ndusttypes
          a_SPH(l+1) = grainsize(l)
       enddo
       ! temorary for old phantom dumps were grain sizes were not defined
       if (grainsize(1) < 1) then
          a_SPH(2) = 1000. ;
          write(*,*) "WARNING: forcing big grains to be 1mm"
       endif
       write(*,*) "WARNING: HARD-CODED: assuming dust grains smaller than 1mum are following the gas"
       a_SPH(1) = 1. ;

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
                stop
             endif
          endif

          if (a_SPH(l) > 0.) then
             log_a_SPH(l) = log(a_SPH(l))
          else
             write(*,*) "ERROR : grains size must be > 0 in SPH file"
             write(*,*) "Exiting"
             stop
          endif
       enddo

       ! I need to work with masses
       Mtot = 0. ;
       Mtot_dust = 0. ;
       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) then
             masse_gaz(icell)    = massgas(iSPH) /  g_to_Msun
             densite_gaz(icell)  = masse_gaz(icell) /  (AU3_to_cm3 * masse_mol_gaz * volume(icell))

             Mtot = Mtot +  masse_gaz(icell)
             !Mtot_dust = Mtot_dust + masse_gaz(icell)* (sum(dustfrac(:,iSPH)))
             Mtot_dust = Mtot_dust + sum(massdust(:,icell))
          else ! star
             masse_gaz(icell)    = 0.
             densite_gaz(icell)  = 0.
          endif
       enddo
       Mtot_dust = Mtot_dust / g_to_Msun

       ! Remark : you may want to had a factor here if there are only a few grainsize
       ! in the SPH dump
       dust_to_gas = Mtot_dust/Mtot

       do icell=1,n_cells
          iSPH = Voronoi(icell)%id
          if (iSPH > 0) then
             ! mass & density indices are shifted by 1
             do l=1, ndusttypes+1
                if (l==1) then
                   ! small grains follow the gas, we do not care about normalization here
                   rho_dust(l) = massgas(iSPH) / (100. * volume(icell))
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
             !stop
          else ! iSPH == 0, star
             densite_pouss(:,icell) = 0.
          endif
       enddo ! icell

       ! Normalisation en taille de grain: on a 1 grain de chaque taille dans le disque,
       ! puis 1 grain en tout dans le disque
       do l=1,n_grains_tot
          somme=0.0
          do icell=1,n_cells
             if (densite_pouss(l,icell) <= 0.0) densite_pouss(l,icell) = tiny_dp
             somme=somme+densite_pouss(l,icell)*volume(icell)
          enddo !icell
          densite_pouss(l,:) = densite_pouss(l,:) * (nbre_grains(l)/somme)
       enddo !l

       ! Normalisation : Calcul masse totale
       mass = 0.0
       do icell=1,n_cells
          do l=1,n_grains_tot
             mass=mass + densite_pouss(l,icell) * M_grain(l) * (volume(icell) * AU3_to_cm3)
          enddo !l
       enddo !icell
       mass =  mass !/Msun_to_g
       densite_pouss(:,:) = densite_pouss(:,:) * sum(masse_gaz)/mass * dust_to_gas

       do icell=1,n_cells
          masse(icell) = 0.
          do l=1,n_grains_tot
             masse(icell) = masse(icell) + densite_pouss(l,icell) * M_grain(l) * volume(icell)
          enddo !l
       enddo ! icell
       masse(:) = masse(:) * AU3_to_cm3

    else ! ndusttypes = 0 : using the gas density

       lvariable_dust = .false.
       write(*,*) "Using gas to dust ratio in mcfost parameter file"

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
    endif

    write(*,*) 'Total  gas mass in model:',  real(sum(masse_gaz) * g_to_Msun),' Msun'
    write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'

    search_not_empty : do k=1,n_grains_tot
       do icell=1, n_cells
          if (densite_pouss(k,icell) > 0.0_dp) then
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

    isochrone_file = "Siess/isochrone_3Myr.txt"

    write(*,*) ""
    write(*,*) "Reading isochrone file: "//trim(isochrone_file)
    filename = trim(mcfost_utils)//"Isochrones/"//trim(isochrone_file)
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
       etoile(i)%lb_body = .true.
       write(*,*) "Star #",i,"  Teff=", etoile(i)%T, "K, r=", etoile(i)%r, "Rsun"
    enddo
    write(*,*) ""

    ! Passage rayon en AU
    etoile(:)%r = etoile(:)%r * Rsun_to_AU

    return

  end subroutine compute_stellar_parameters

  !*********************************************************


  subroutine read_ascii_SPH_file(iunit,filename,x,y,z,massgas,rhogas,rhodust,ndusttypes,n_SPH,ierr)

    integer,               intent(in) :: iunit
    character(len=*),      intent(in) :: filename
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,rhogas,massgas
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
    allocate(x(n_SPH),y(n_SPH),z(n_SPH),massgas(n_SPH),rhogas(n_SPH),rhodust(ndusttypes,n_SPH), stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error in phanton_2_mcfost"
       write(*,*) "Exiting"
       stop
    endif

    open(unit=1, file=filename, status='old', iostat=ios)
    do i=1, n_SPH
       read(1,*) x(i), y(i), z(i), massgas(i)
       rhogas(i) = massgas(i)
    enddo

    write(*,*) "MinMax=", minval(massgas), maxval(massgas)

    write(*,*) "Using stars from mcfost parameter file"

    return

  end subroutine read_ascii_SPH_file

end module SPH2mcfost
