module pluto_mod
  !use density, only : densite_gaz, masse_gaz !used to build nHtot
  use parametres
  use messages
  use constantes, only : masseH
  use getline
  use atmos_type, only : alloc_atomic_atmos, alloc_magnetic_field, wght_per_H, nHtot, ne, &
       v_char, lmagnetized, vturb, T, icompute_atomRT, calc_ne, laccretion_shock, &
       Taccretion, write_atmos_domain

  use utils
  use Voronoi_grid, only : Voronoi, volume, neighbours_list
  use sph2mcfost, only : SPH_to_Voronoi
  ! 	use grid,

  implicit none

contains


  subroutine setup_mhd_to_mcfost()
    !here lmhd_voronoi ist true even with model ascii !
    integer										:: n_points, i, j, k, icell
    integer										:: Nread
    real(kind=dp), allocatable, dimension(:)	:: x,y,z,h!allocated in reading
    real(kind=dp), allocatable, dimension(:)	:: vx,vy,vz,mass_gas, mass_ne_on_massgas, T_tmp, vt_tmp, dz
    real(kind=dp), allocatable, dimension(:)	:: rho, hydro_grainsizes
    real(kind=dp), allocatable, dimension(:,:)	:: rhodust, massdust
    logical, allocatable, dimension(:)			:: mask

    real(kind=dp)								:: thetai, thetao, tilt
    integer, parameter							:: Nheader = 3 !Add more, for ascii file
    character(len=MAX_LENGTH)					:: rotation_law
    integer										:: syst_status, acspot, alloc_status
    character(len=MAX_LENGTH)					:: inputline, FormatLine, cmd
    real(kind=dp)								:: rho_to_nH, vxmax, vxmin, vymax, vymin, vzmax, vzmin, vmax
    logical										:: check_previous_tesselation
    integer,  allocatable, dimension(:)			:: particle_id
    real(kind=dp), dimension(6)					:: hydro_limits
    integer										:: ndusttypes, voroindex, N_fixed_ne = 0
    real, parameter								:: limit_factor = 1.005, Lextent = 1.01

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

       !unused her, but still present
       call getnextline(1, "#", FormatLine, inputline, Nread)
       read(inputline(1:Nread),*) thetai, thetao, tilt

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

       allocate(T_tmp(n_points), vt_tmp(n_points), dz(n_points), mass_gas(n_points), mass_ne_on_massgas(n_points), stat=alloc_status)
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
             read(inputline(1:Nread),*) x(icell), y(icell), z(icell), T_tmp(icell), mass_gas(icell), mass_ne_on_massgas(icell), vx(icell), vy(icell), vz(icell), vt_tmp(icell), dz(icell)
          endif
       enddo
       !density_file

       h(:) = maxval(sqrt(x**2+y**2+z**2))

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
    hydro_limits(1) = select_inplace(k,real(x))*limit_factor
    hydro_limits(3) = select_inplace(k,real(y))*limit_factor
    hydro_limits(5) = select_inplace(k,real(z))*limit_factor

    k = n_points
    hydro_limits(2) = select_inplace(k,real(x))*limit_factor
    hydro_limits(4) = select_inplace(k,real(y))*limit_factor
    hydro_limits(6) = select_inplace(k,real(z))*limit_factor
    write(*,*) "# Model limits (Rstar) :"
    write(*,*) "x =", hydro_limits(1)/etoile(1)%r, hydro_limits(2)/etoile(1)%r
    write(*,*) "y =", hydro_limits(3)/etoile(1)%r, hydro_limits(4)/etoile(1)%r
    write(*,*) "z =", hydro_limits(5)/etoile(1)%r, hydro_limits(6)/etoile(1)%r

    !also work with grid-based code
    !massdust, rhodust, hydro_grainsizes not allocated if ndusttypes = 0 !
    !** bug handling
    !     	ndusttypes = 1 !oterwise bug
    !     	allocate(massdust(1, n_points),rhodust(1, n_points),hydro_grainsizes(1))
    !** bug handling
    call sph_to_voronoi(n_points, ndusttypes, particle_id, x, y, z, h, vx, vy, vz, mass_gas, massdust, rho, rhodust, hydro_grainsizes, hydro_limits, check_previous_tesselation)
    ! -> correction for small density applied on mass_gas directly inside

    !** bug handling
    !     	deallocate(massdust,rhodust,hydro_grainsizes)
    !** bug handling

    !*******************************
    ! Accomodating model
    !*******************************
    !alloc space for all physical quantities.
    !Velocity field arrays not allocated because lvoronoi is true
    !-> fills element abundances structures for elements
    call alloc_atomic_atmos
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
          ne(icell) = Msun_to_kg * mass_ne_on_massgas(voroindex) * mass_gas(voroindex) / (volume(icell) * AU3_to_m3)

          T(icell) = T_tmp(voroindex)

          vturb(icell) = vt_tmp(voroindex)

          icompute_atomRT(icell) = dz(voroindex)

          vxmax = max(vxmax, abs(Voronoi(icell)%vxyz(1)))
          vxmin = min(vxmin, max(abs(Voronoi(icell)%vxyz(1)),0.0))
          vymax = max(vymax, abs(Voronoi(icell)%vxyz(2)))
          vymin = min(vymin,  max(abs(Voronoi(icell)%vxyz(2)),0.0))
          vzmax = max(vzmax, abs(Voronoi(icell)%vxyz(3)))
          vzmin = min(vzmin,  max(abs(Voronoi(icell)%vxyz(3)),0.0))

          Vmax = max( vmax, sqrt(Voronoi(icell)%vxyz(1)**2 + Voronoi(icell)%vxyz(2)**2 + Voronoi(icell)%vxyz(3)**2) )
       endif
    end do

    !deallocating temporary variables from input file.
    deallocate(h,vx,vy,vz,mass_gas, mass_ne_on_massgas, x,y,z,T_tmp, vt_tmp, dz)

    v_char = Vmax * Lextent


    write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT>0)), " density zones"
    write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT==0)), " transparent zones"
    write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT<0)), " dark zones"



    ! 		if (lmagnetized) then
    !
    ! 		endif


    calc_ne = .false.
    icell_loop : do icell=1,n_cells
       !check that in filled cells there is electron density otherwise we need to compute it
       !from scratch.
       if (icompute_atomRT(icell) > 0) then

          if (ne(icell) <= 0.0_dp) then
             write(*,*) "  ** No electron density found in the model! ** "
             calc_ne = .true.
             exit icell_loop
          endif

       endif
    enddo icell_loop
    N_fixed_ne = size(pack(icompute_atomRT,mask=(icompute_atomRT==2)))
    if (N_fixed_ne > 0) then
       write(*,'("Found "(1I5)" cells with fixed electron density values! ("(1I3)" %)")') N_fixed_ne, nint(real(N_fixed_ne) / real(n_cells) * 100)
    endif

    call write_atmos_domain() !but for consistency with the plot functions in python


    write(*,*) "Maximum/minimum velocities in the model (km/s):"
    write(*,*) "|Vx|", vxmax*1d-3, vxmin*1d-3
    write(*,*) "|Vy|", vymax*1d-3, vymin*1d-3
    write(*,*) "|Vz|", vzmax*1d-3, vzmin*1d-3



    write(*,*) "Typical line extent due to V fields (km/s):"
    write(*,*) v_char/1d3

    write(*,*) "Maximum/minimum turbulent velocity (km/s):"
    write(*,*) maxval(vturb)/1d3, minval(vturb, mask=icompute_atomRT>0)/1d3


    write(*,*) "Maximum/minimum Temperature in the model (K):"
    write(*,*) MAXVAL(T), MINVAL(T,mask=icompute_atomRT>0)
    write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
    write(*,*) MAXVAL(nHtot), MINVAL(nHtot,mask=icompute_atomRT>0)
    if (.not.calc_ne) then
       write(*,*) "Maximum/minimum ne density in the model (m^-3):"
       write(*,*) MAXVAL(ne), MINVAL(ne,mask=icompute_atomRT>0)
    endif


    return
  end subroutine setup_mhd_to_mcfost


  subroutine read_pluto()
    !read pluto format in hdf5

    call error("Pluto interface not available yet!")

    return
  end subroutine read_pluto

end module pluto_mod
