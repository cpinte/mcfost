module pluto_mod
  use parametres
  use messages
  use getline
  use atmos_type, only : lmagnetized, laccretion_shock, Taccretion, empty_cells, thetao, thetai

  use utils
  use sph2mcfost, only : SPH_to_Voronoi, Hydro_to_Voronoi_atomic

  implicit none

contains


  subroutine setup_mhd_to_mcfost()
    !here lmhd_voronoi ist true even with model ascii !
    integer                                  :: n_points,k,icell
    integer                                  :: Nread
    real(kind=dp), allocatable, dimension(:) :: x,y,z,h!allocated in reading
    real(kind=dp), allocatable, dimension(:) :: vx,vy,vz,mass_gas, mass_ne_on_massgas, T_tmp, vt_tmp, dz
    real(kind=dp), allocatable, dimension(:) :: rho, hydro_grainsizes
    real(kind=dp), allocatable, dimension(:,:) :: rhodust, massdust

    real(kind=dp)                            :: tilt
    integer, parameter                       :: Nheader = 3 !Add more, for ascii file
    character(len=MAX_LENGTH)                :: rotation_law
    integer                                  :: syst_status, acspot, alloc_status
    character(len=MAX_LENGTH)                :: inputline, FormatLine, cmd
    logical                                  :: check_previous_tesselation
    integer,  allocatable, dimension(:)      :: particle_id
    real(kind=dp), dimension(6)              :: hydro_limits
    integer                                  :: ndusttypes !, voroindex, N_fixed_ne = 0
    real, parameter                          :: limit_factor = 1.005!, Lextent = 1.01

    write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH

    !There will be an error if lphantom_file is true. Because density_files and density_file
    !stores pluto's model name. But also the filename from phantom.. So to date, the two
    !codes cannot be merged.
    !This is to be able to use read previous tesselation


    lfix_star = .true.
    lphantom_file = .false.
    lignore_dust = .true.
    lrandomize_voronoi = .false.
    check_previous_tesselation = (.not.lrandomize_voronoi)
    n_points = 0 ! to avoid compiler warning

    if (lignore_dust) then
       ndusttypes = 0
       if (allocated(rhodust)) deallocate(rhodust,massdust)

    else
       call error("Dust not handled yet for pluto models!")
    endif

    !lvoronoi = .true. -> set to .true. if lmhd_voronoi
    lmagnetized = lzeeman_polarisation !if Zeeman polarisation read magnetic field!
    if (lmagnetized) then
       call warning("Magnetic field not read yet with Voronoi grid. lzeeman_polarisation set to .false.!")
       lzeeman_polarisation = .false.
       lmagnetized = .false.
       !!call alloc_magnetic_field()
    endif

    if (lpluto_file) then
       write(*,*) "Voronoi tesselation on Pluto model..."
       !read and convert to mcfost units
       call read_pluto() ! to date empty
    elseif (lmodel_ascii) then
       !read density_file
       !Reading ascii file in the same format as atmos_type.f90 / readAtmos_ascii()
       !but I intend to deprecate readAtmos_ascii.

       cmd = "wc -l "//trim(density_file)//" > ntest.txt"
       call appel_syst(cmd,syst_status)
       open(unit=1,file="ntest.txt",status="old")
       read(1,*) N_points
       close(unit=1)
       N_points = N_points - Nheader
       write(*,*) " Input model has ", N_points," grid points !"

       N_points = N_points + n_etoiles

       open(unit=1,file=density_file, status="old")
       call getnextline(1, "#", FormatLine, inputline, Nread)
       read(inputline(1:Nread),*) rotation_law

       lmagnetoaccr = .false.
       lspherical_velocity = .false.
       !-> .false. with Voronoi
       !read T shock and if accretion spots
       call getnextline(1, "#", FormatLine, inputline, Nread)
       read(inputline(1:Nread),*) Taccretion, acspot
       laccretion_shock =  (acspot == 1)
       if (Taccretion==0.0_dp) Taccretion = -1.0_dp

       !unused her, but still present
       call getnextline(1, "#", FormatLine, inputline, Nread)
       read(inputline(1:Nread),*) thetai, thetao, tilt
       !for compatibility
       tilt = tilt * pi / 180.0
       thetai = thetai * pi / 180.0
       thetao = thetao * pi / 180.

       allocate(h(n_points), stat=alloc_status)
       if (alloc_status > 0) then
          call error("Allocation error smoothing length h")
       endif
       !cut cells larger than 3*h

       allocate(particle_id(n_points), stat=alloc_status)
       if (alloc_status > 0) then
          call error("Allocation error particle_id (mhd_to_mcfost)")
       endif
       particle_id(:) = 0


       allocate(x(n_points), y(n_points), z(n_points), stat=alloc_status)
       if (alloc_status > 0) then
          call error("Allocation error x, y, z")
       endif

       allocate(vx(n_points), vy(n_points), vz(n_points), stat=alloc_status)
       if (alloc_status > 0) then
          call error("Allocation error vx, vy, vz")
       endif

       allocate(T_tmp(n_points), vt_tmp(n_points), dz(n_points), mass_gas(n_points), mass_ne_on_massgas(n_points), &
            stat=alloc_status)
       if (alloc_status > 0) then
          call error("Allocation error T_tmp")
       endif

       icell = 0
       Nread = 0
       do icell=1, n_points
          particle_id(icell) = icell
          if (lmagnetized) then
             call error("Magnetic field not available with Voronoi!")
          else
             call getnextLine(1, "#", FormatLine, inputline, Nread)
             read(inputline(1:Nread),*) x(icell), y(icell), z(icell), T_tmp(icell), mass_gas(icell),&
                  mass_ne_on_massgas(icell), vx(icell), vy(icell), vz(icell), vt_tmp(icell), dz(icell), h(icell)
          endif
       enddo
       !density_file

		!Included in the input model now.
!        h(:) = maxval(sqrt(x**2+y**2+z**2))

    else
       call error("lpluto_file or lmodel_ascii required for lmhd_voronoi!")
    end if


    !???
    ! 		if ((.not.lfix_star).and.(lpluto_file.or.lmodel_ascii)) then
    ! ! 			write(*,*) "Error, lfix_star has to be .true. using pluto models !"
    ! ! 			stop
    ! 			!call compute_stellar_parameters()
    ! 		end if

    ! Setting limits explicitely here !
    !call read_SPH_limits_file(" ", .false.)
    hydro_limits(:) = 0

    k = 1
    hydro_limits(1) = find_kth_smallest_inplace(k,real(x))*limit_factor
    hydro_limits(3) = find_kth_smallest_inplace(k,real(y))*limit_factor
    hydro_limits(5) = find_kth_smallest_inplace(k,real(z))*limit_factor

    k = n_points
    hydro_limits(2) = find_kth_smallest_inplace(k,real(x))*limit_factor
    hydro_limits(4) = find_kth_smallest_inplace(k,real(y))*limit_factor
    hydro_limits(6) = find_kth_smallest_inplace(k,real(z))*limit_factor
    if (n_etoiles > 0) then
    	write(*,*) "# Model limits (Rstar) :"
    	write(*,*) "x =", hydro_limits(1)/etoile(1)%r, hydro_limits(2)/etoile(1)%r
    	write(*,*) "y =", hydro_limits(3)/etoile(1)%r, hydro_limits(4)/etoile(1)%r
    	write(*,*) "z =", hydro_limits(5)/etoile(1)%r, hydro_limits(6)/etoile(1)%r
    endif

    !also work with grid-based code
    !massdust, rhodust, hydro_grainsizes not allocated if ndusttypes = 0 !
    call sph_to_voronoi(n_points-n_etoiles, ndusttypes, particle_id, x, y, z, h, vx, vy, vz, &
         T_tmp, mass_gas, massdust, rho, rhodust, hydro_grainsizes, hydro_limits, check_previous_tesselation)
    ! -> correction for small density applied on mass_gas directly inside

    call hydro_to_Voronoi_atomic(n_points,T_tmp,vt_tmp,mass_gas,mass_ne_on_massgas,dz)
! 	call empty_cells

    !deallocating temporary variables from input file.
    deallocate(h,vx,vy,vz,mass_gas, mass_ne_on_massgas, x,y,z,T_tmp, vt_tmp, dz)

    return
  end subroutine setup_mhd_to_mcfost


  subroutine read_pluto()
    !read pluto format in hdf5

    call error("Pluto interface not available yet!")

    return
  end subroutine read_pluto

end module pluto_mod
