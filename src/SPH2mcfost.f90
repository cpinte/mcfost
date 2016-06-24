module SPH2mcfost

  use parametres
  use constantes
  use utils

  implicit none

contains


  subroutine setup_SPH2mcfost(phantom_file, phantom_limits_file)

    use io_phantom, only : read_phantom_file, read_phantom_input_file
    use dump_utils, only : get_error_text
    use Voronoi_grid
    use opacity, only : densite_pouss, masse
    use molecular_emission, only : densite_gaz, masse_gaz
    use grains, only : n_grains_tot, M_grain
    use mem

    character(len=512), intent(in) :: phantom_file, phantom_limits_file

    integer, parameter :: iunit = 1
    real(db), allocatable, dimension(:) :: x,y,z,rho,massgas
    real(db), allocatable, dimension(:,:) :: rhodust
    real, allocatable, dimension(:) :: a_SPH
    real :: grainsize,graindens, f
    integer :: ierr, n_SPH, n_Voronoi, ndusttypes, alloc_status, icell, l, k, ios

    logical :: lwrite_ASCII = .true. ! produce an ASCII file for yorick

    character(len=100) :: line_buffer
    real(db), dimension(6) :: limits

    icell_ref = 1

    write(*,*) "Performing phantom2mcfost setup"
    write(*,*) "Reading phantom density file: "//trim(phantom_file)

    if (lphantom_test) then
       write(*,*) "USING PHANTOM TEST MODE"
       call read_test_ascii_file(iunit,phantom_file, x,y,z,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    else
       call read_phantom_file(iunit,phantom_file, x,y,z,massgas,rho,rhodust,ndusttypes,n_SPH,ierr)
    endif
    if (ierr /=0) then
       write(*,*) "Error code =", ierr,  get_error_text(ierr)
       stop
    endif
    write(*,*) "# Farthest particules :"
    write(*,*) "x =", minval(x), maxval(x)
    write(*,*) "y =", minval(y), maxval(y)
    write(*,*) "z =", minval(z), maxval(z)

    !if (ndusttypes==1) then
    !   call read_phantom_input_file("hltau.in",iunit,grainsize,graindens,ierr)
    !   write(*,*) grainsize,graindens
    !endif
    write(*,*) "Found", n_SPH, " SPH particles with ", ndusttypes, "dust grains"
    allocate(a_SPH(ndusttypes))

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


    !*******************************
    ! Model limits
    !*******************************
    write(*,*) " "
    write(*,*) "Reading phantom limits file: "//trim(phantom_limits_file)
    open(unit=1, file=phantom_limits_file, status='old', iostat=ios)
    if (ios/=0) then
       write(*,*) "ERROR : cannot open "//trim(phantom_limits_file)
       write(*,*) "Exiting"
       stop
    endif
    read(1,*) line_buffer
    read(1,*) limits(1), limits(3), limits(5)
    read(1,*) limits(2), limits(4), limits(6)
    write(*,*) "# Model limits :"
    write(*,*) "x =", limits(1), limits(2)
    write(*,*) "y =", limits(3), limits(4)
    write(*,*) "z =", limits(5), limits(6)
    close(unit=1)

    !*******************************
    ! Voronoi tesselation
    !*******************************
    ! Make the Voronoi tesselation on the SPH particles ---> define_Voronoi_grid : volume
    !call Voronoi_tesselation_cmd_line(n_SPH, x,y,z, limits, n_Voronoi)
    call Voronoi_tesselation(n_SPH, x,y,z, limits, n_Voronoi)
    deallocate(x,y,z)
    write(*,*) "Using n_cells =", n_cells
    !call test_walls()

    !*************************
    ! Densities
    !*************************
    call allocate_densities()
    ! Tableau de densite et masse de gaz
    !do icell=1,n_cells
    !   densite_gaz(icell) = rho(icell) / masse_mol_gaz * m3_to_cm3 ! rho is in g/cm^3 --> part.m^3
    !   masse_gaz(icell) =  densite_gaz(icell) * masse_mol_gaz * volume(icell)
    !enddo
    !masse_gaz(:) = masse_gaz(:) * AU3_to_cm3

    do icell=1,n_cells
       masse_gaz(icell) = massgas(icell) /  g_to_Msun
       densite_gaz(icell)  = masse_gaz(icell) /  (AU3_to_cm3 * masse_mol_gaz * volume(icell))
    enddo

    ! Tableau de densite et masse de poussiere
    ! interpolation en taille
    if (ndusttypes > 1) then
       lvariable_dust = .true.
       write(*,*) "*********************************************"
       write(*,*) "This part has not been tested"
       write(*,*) "Dust mass is going to incorrect !!!"
       write(*,*) "rhodust is not calibrated for mcfost yet"
       write(*,*) "*********************************************"
       l=1
       do icell=1,n_cells
          do k=1,n_grains_tot
             if (r_grain(l) < a_SPH(1)) then ! small grains
                densite_pouss(k,icell) = rhodust(1,icell)
             else if (r_grain(k) < a_SPH(ndusttypes)) then ! large grains
                densite_pouss(k,icell) = rhodust(ndusttypes,icell)
             else ! interpolation
                if (r_grain(k) > a_sph(l+1)) l = l+1
                f = (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))

                densite_pouss(k,icell) = rhodust(l,icell) + f * (rhodust(l+1,icell)  - rhodust(l,icell))
             endif
             !write(*,*) "Todo : densite_pouss : missing factor"
             !          stop
             masse(icell) = masse(icell) + densite_pouss(k,icell) * M_grain(k) * volume(icell)
          enddo !l
       enddo ! icell
       masse(:) = masse(:) * AU3_to_cm3
    else ! using the gas density
       lvariable_dust = .false.
       write(*,*) "Forcing gas/dust == 100"
       do icell=1,n_cells
          do k=1,n_grains_tot
             densite_pouss(k,icell) = densite_gaz(icell) * nbre_grains(k)
             masse(icell) = masse(icell) + densite_pouss(k,icell) * M_grain(k) * volume(icell)
          enddo
       enddo
       masse(:) = masse(:) * AU3_to_cm3
       f = 0.01 * sum(masse_gaz)/sum(masse)
       densite_pouss(:,:) = densite_pouss(:,:) * f
       masse(:) = masse(:) * f
    endif

    write(*,*) 'Total  gas mass in model:', real(sum(masse_gaz) * g_to_Msun),' Msun'
    write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'
    deallocate(massgas,rho,rhodust,a_SPH)

    search_not_empty : do k=1,n_grains_tot
       do icell=1, n_cells
          if (densite_pouss(k,icell) > 0.0_db) then
             icell_not_empty = icell
             exit search_not_empty
          endif
       enddo !icell
    enddo search_not_empty

    write(*,*) "End setup_phantom2mcfost"

    return

  end subroutine setup_SPH2mcfost

  !************

  subroutine read_test_ascii_file(iunit,filename,x,y,z,massgas,rhogas,rhodust,ndusttypes,n_SPH,ierr)

    integer,               intent(in) :: iunit
    character(len=*),      intent(in) :: filename
    real(db), intent(out), dimension(:),   allocatable :: x,y,z,rhogas,massgas
    real(db), intent(out), dimension(:,:), allocatable :: rhodust
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
    endif

    write(*,*) shape(x)

    open(unit=1, file=filename, status='old', iostat=ios)
    do i=1, n_SPH
       read(1,*) x(i), y(i), z(i), massgas(i)
       rhogas(i) = massgas(i)
    enddo

    write(*,*) "MinMax=", minval(massgas), maxval(massgas)

    return

  end subroutine read_test_ascii_file

end module SPH2mcfost
