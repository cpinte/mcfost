module read_idefix

  use parametres
  use messages
  use mcfost_env
  use VTKtools
  use grid
  use cylindrical_grid
  use density
  use stars, only : compute_stellar_parameters
  use read_fargo3d, only : convert_planets, n_planets_max

  implicit none

contains

  subroutine read_idefix_parameters(filename)

    character(len=*), intent(in) :: filename

    integer :: geometry, time
    integer, dimension(3) :: dimensions
    character(:), allocatable :: origin
    real, dimension(:), allocatable :: x1, x2, x3
    real(dp) :: dx

    idefix%filename = filename
    call readVTK_header(filename, idefix%iunit, idefix%position, idefix%geometry, dimensions, idefix%time, origin, x1, x2, x3)

    idefix%origin = origin
    idefix%dimensions = dimensions
    idefix%nx1 = dimensions(1) ! r
    idefix%nx2 = dimensions(2) ! theta
    idefix%nx3 = dimensions(3) ! phi

    ! Checking is the grid is linear
    idefix%log_spacing = (x1(idefix%nx1) - x1(idefix%nx1-1)) > 1.01 * (x1(2) - x1(1))

    ! Model boundaries
    idefix%x1_min = x1(1)           ! r
    idefix%x1_max = x1(idefix%nx1)
    idefix%x2_min = x2(1)           ! phi in cylindrical, theta in spherical
    idefix%x2_max = x2(idefix%nx2)
    idefix%x3_min = x3(1)           ! z in cylindrical, phi in spherical
    idefix%x3_max = x3(idefix%nx3)
    ! phi is 0 to 2*pi, while fargo3d is -pi to pi

    allocate(idefix%x1(idefix%nx1), idefix%x2(idefix%nx2), idefix%x3(idefix%nx3))
    idefix%x1 = x1 ; idefix%x2 = x2 ; idefix%x3 = x3

    ! Updating mcfost parameters
    grid_type = idefix%geometry
    n_rad = idefix%nx1-1
    n_rad_in = 1

    if (idefix%geometry == 1) then
       call warning("idefix : using cylindrical grid")
       nz = (idefix%nx3-1)/2
       n_az = idefix%nx2-1
    else if (idefix%geometry ==2) then
       call warning("idefix : using spherical grid")
       nz = (idefix%nx2-1)/2+1
       n_az = idefix%nx3-1
       lregular_theta = .true.
       theta_max = 0.5 * pi - idefix%x2_min
    else
       call error("idefix: unknow geometry")
    endif

    if (lscale_length_units) then
       write(*,*) 'Lengths are rescaled by ', real(scale_length_units_factor)
    else
       scale_length_units_factor = 1.0
    endif

    disk_zone(1)%rin  = idefix%x1_min * scale_length_units_factor
    disk_zone(1)%edge=0.0
    disk_zone(1)%rmin = disk_zone(1)%rin
    disk_zone(1)%rout = idefix%x1_max * scale_length_units_factor
    disk_zone(1)%rmax = disk_zone(1)%rout

    write(*,*) "n_rad=", n_rad, "nz=", nz, "n_az=", n_az
    write(*,*) "rin=", real(disk_zone(1)%rin), "rout=", real(disk_zone(1)%rout)

    return

  end subroutine read_idefix_parameters

  !---------------------------------------------

  subroutine read_idefix_model()

    ! idefix data is ordered in x1 = r, x2 = theta, x3 = phi

    real, dimension(:,:,:), allocatable  :: rho, vx1, vx2, vx3
    integer :: ios, iunit, alloc_status, l, recl, i,j, jj, phik, icell, id, n_planets, n_etoiles_old, i2, i3

    character(len=128) :: filename
    character(len=16), dimension(4) :: file_types

    real(dp) :: Ggrav_idefix, umass, usolarmass, ulength, utime, udens, uvelocity, ulength_au, mass, facteur, omega
    type(star_type), dimension(:), allocatable :: etoile_old

    real(dp), dimension(n_planets_max) :: x, y, z, vx, vy, vz, Mp, Omega_p, time

    ! Planet properties hard coded for now
    !real, parameter :: Mp = 1e-3
    !real, parameter :: Omega_p = 1.0
    !real, parameter :: x = 1.0, y=0.0, z=0.0
    !real, parameter :: vx=0.0, vy=1.0, vz=1.0

    usolarmass = 1.0_dp
    ulength_au = 1.0_dp
    Ggrav_idefix = 1.0_dp

    if (lscale_length_units) then
       ulength_au = ulength_au * scale_length_units_factor
    else
       scale_length_units_factor = 1.0
    endif

    if (lscale_mass_units) then
       write(*,*) 'Mass are rescaled by ', real(scale_mass_units_factor)
       usolarmass = usolarmass * scale_mass_units_factor
    else
       scale_mass_units_factor = 1.0
    endif

    umass = usolarmass *  Msun_to_kg
    ulength = ulength_au * AU_to_m
    utime = sqrt(ulength**3/((Ggrav/Ggrav_idefix)*umass))

    udens = umass / ulength**3
    uvelocity = ulength / utime

    call read_idefix_planets("./", n_planets,x,y,z,vx,vy,vz,Mp,time)
    ! Omega from .ini
    call convert_planets(n_planets, x,y,z,vx,vy,vz,Mp,time,Omega_p,ulength_au,uvelocity,usolarmass,utime)
    idefix%corrotating_frame = (n_planets > 0) ! todo : from ini

    if (idefix%corrotating_frame) then
       Omega = Omega_p(which_planet)
    else
       Omega = 0.0_dp
    endif

    ! reading data
    write(*,*) "Reading Idefix data ..."
    call readVTK_data(idefix%iunit, idefix%position, idefix%dimensions, rho, vx1, vx2, vx3)
    write(*,*) "Done"

    !-----------------------------------
    ! Passing data to mcfost
    !-----------------------------------

    ! copied from read_athena++

    write(*,*) "Converting idefix files to mcfost ..."
    lvelocity_file = .true.
    vfield_coord = 3 ! spherical

    allocate(vfield3d(n_cells,3), stat=alloc_status)
    if (alloc_status /= 0) call error("memory allocation error idefix vfield3d")

    write(*,*) "Constant spatial distribution"

    do i=1, n_rad
       jj= 0
       bz : do j=j_start+1,nz-1 ! 1 extra empty cell in theta on each side
          if (j==0) cycle bz
          jj = jj + 1
          do phik=1, n_az
             ! Order of dimensions is not the same as mcfost depending on geometry
             if (lcylindrical) then
                i2 = phik
                i3 = jj
             else
                i2 = jj
                i3 = phik
             endif


             icell = cell_map(i,j,phik)

             densite_gaz(icell) = rho(i,i2,i3) * udens
             densite_pouss(:,icell) = rho(i,i2,i3) * udens

             ! todo : check in cyl
             vfield3d(icell,1)  = vx1(i,i2,i3) * uvelocity ! vr
             vfield3d(icell,2)  = (vx3(i,i2,i3) + r_grid(icell)/ulength_au * Omega) * uvelocity ! vphi : planet at r=1
             vfield3d(icell,3)  = vx2(i,i2,i3) * uvelocity ! vtheta
          enddo ! k
       enddo bz
    enddo ! i
    deallocate(rho,vx1,vx2,vx3)

    ! -- another copy and paste from read_fargo3d
    ! Normalisation density : copy and paste from read_density_file for now : needs to go in subroutine

    ! Calcul de la masse de gaz de la zone
    mass = 0.
    do icell=1,n_cells
       mass = mass + densite_gaz(icell) *  masse_mol_gaz * volume(icell)
    enddo !icell
    mass =  mass * AU3_to_m3 * g_to_Msun

    ! Normalisation
    if (mass > 0.0) then ! pour le cas ou gas_to_dust = 0.
       facteur = disk_zone(1)%diskmass * disk_zone(1)%gas_to_dust / mass

       ! Somme sur les zones pour densite finale
       do icell=1,n_cells
          densite_gaz(icell) = densite_gaz(icell) * facteur
          masse_gaz(icell) = densite_gaz(icell) * masse_mol_gaz * volume(icell) * AU3_to_m3
       enddo ! icell
    else
       call error('Gas mass is 0')
    endif

    write(*,*) 'Total  gas mass in model:', real(sum(masse_gaz) * g_to_Msun),' Msun'
    call normalize_dust_density()
    ! -- end copy and paste from read_fargo3d

    write(*,*) "Done"
    return

  end subroutine read_idefix_model

  !---------------------------------------------

  subroutine read_idefix_planets(dir, n_planets,x,y,z,vx,vy,vz,Mp,time)

    character(len=*), intent(in) :: dir
    integer, intent(out) :: n_planets
    real(dp), dimension(n_planets_max), intent(out) :: x, y, z, vx, vy, vz, Mp, time

    integer :: n_etoiles_old, iunit, ios, n_etoile_old, i, i_planet, id

    character(len=1) :: s
    character(len=128) :: filename

    ios = 0
    iunit = 1

    n_planets = 0
    planet_loop : do i_planet=1, n_planets_max
       write(s,"(I1)") i_planet-1
       filename=dir//"/planet"//s//".dat"
       open(unit=iunit, file=filename, status="old", form="formatted", iostat=ios)
       if (ios /= 0) exit planet_loop
       n_planets = n_planets+1
       write(*,*) "Reading "//trim(filename)
       read(fargo3d%id,*) id
       do while(ios==0)
          read(iunit,*) time(i_planet), x(i_planet), y(i_planet), z(i_planet), vx(i_planet), vy(i_planet), vz(i_planet), &
               Mp(i_planet)
          if (i==id) exit
       enddo
       close(iunit)
    enddo planet_loop

    write(s,"(I1)") n_planets
    write(*,*) "Found "//s// " planets"

    return

  end subroutine read_idefix_planets

  !---------------------------------------------

end module read_idefix
