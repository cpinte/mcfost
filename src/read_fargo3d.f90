module read_fargo3d

  use parametres
  use messages
  use mcfost_env
  use grid
  use cylindrical_grid
  use density
  use stars, only : compute_stellar_parameters

  implicit none

contains

  subroutine read_fargo3d_parameters(dir,id)

    character(len=*), intent(in) :: dir, id

    character(len=128) :: key, val, filename
    integer :: ios, n_lines, iunit

    fargo3d%dir = dir
    fargo3d%id = id

    iunit = 1
    ios = 0

    ! Reading parameter files and extracting key parameters
    filename = trim(dir)//"/variables.par"
    open(unit=iunit,file=trim(filename),status='old',form='formatted',iostat=ios)
    if (ios /= 0) call error("opening "//trim(filename))

    ! On compte les lignes avec des donnees
    n_lines = 0
    do while(ios==0)
       n_lines=n_lines+1
       read(1,*,iostat=ios) key, val

       select case(trim(key))
          ! reading grid parameters
       case("NX")
          read(val,*,iostat=ios) fargo3d%nx
       case("NY")
          read(val,*,iostat=ios) fargo3d%ny
       case("NZ")
          read(val,*,iostat=ios) fargo3d%nz
       case("XMIN")
          read(val,*,iostat=ios) fargo3d%xmin
       case("XMAX")
          read(val,*,iostat=ios) fargo3d%xmax
       case("YMIN")
          read(val,*,iostat=ios) fargo3d%ymin
       case("YMAX")
          read(val,*,iostat=ios) fargo3d%ymax
       case("ZMIN")
          read(val,*,iostat=ios) fargo3d%zmin
       case("ZMAX")
          read(val,*,iostat=ios) fargo3d%zmax
       case("SPACING")
          if ((val == "LOG").or.(val == "log").or.(val == "Log")) then
             fargo3d%log_spacing = .true.
          else
             fargo3d%log_spacing = .false.
             llinear_rgrid = .true.
          endif
       case("FRAME")
          if ((val == "C").or.(val == "c")) then
             fargo3d%corrotating_frame = .true.
          else
             fargo3d%corrotating_frame = .false.
          endif
       case("REALTYPE")
          if ((val == "float64").or.(val == "FLOAT64")) then
             fargo3d%realtype = dp
          else
             fargo3d%realtype = sp
          endif

          ! Additional parameters
       case("DT")
          read(val,*,iostat=ios) fargo3d%dt
       case("ASPECTRATIO")
          read(val,*,iostat=ios) fargo3d%aspect_ratio
       case("NU")
          read(val,*,iostat=ios) fargo3d%nu
       case("PLANETCONFIG")
          fargo3d%planetconfig = val
       case("CS")
          read(val,*,iostat=ios) fargo3d%cs
       case("GAMMA")
          read(val,*,iostat=ios) fargo3d%gamma
       end select
    end do

    ! Updating mcfost parameters
    grid_type = 2
    n_rad = fargo3d%ny
    n_rad_in = 1
    n_az = fargo3d%nx
    nz = fargo3d%nz/2 + 1
    lregular_theta = .true.
    theta_max = 0.5 * pi - fargo3d%zmin

    if (lscale_length_units) then
       write(*,*) 'Lengths are rescaled by ', real(scale_length_units_factor)
    else
       scale_length_units_factor = 1.0
    endif

    disk_zone(1)%rin  = fargo3d%ymin * scale_length_units_factor
    disk_zone(1)%edge=0.0
    disk_zone(1)%rmin = disk_zone(1)%rin
    disk_zone(1)%rout = fargo3d%ymax * scale_length_units_factor
    disk_zone(1)%rmax = disk_zone(1)%rout

    write(*,*) "n_rad=", n_rad, "nz=", nz, "n_az=", n_az
    write(*,*) "rin=", real(disk_zone(1)%rin), "rout=", real(disk_zone(1)%rout)

    return

  end subroutine read_fargo3d_parameters

  !---------------------------------------------

  subroutine read_fargo3d_files()

    ! fargo3d data is ordered in x = phi, y = r, z = theta

    real(dp), dimension(:,:,:), allocatable  :: fargo3d_density, fargo3d_vx, fargo3d_vy, fargo3d_vz
    integer :: ios, iunit, alloc_status, l, recl, i,j, jj, phik, icell, id, n_etoiles_old
    real(dp) :: x, y, z, vx, vy, vz, Mp, Omega_p, time

    character(len=128) :: filename
    character(len=16), dimension(4) :: file_types

    real(dp) :: Ggrav_fargo3d, umass, usolarmass, ulength, utime, udens, uvelocity, ulength_au, mass, facteur
    type(star_type), dimension(:), allocatable :: etoile_old

    ! Todo : add option to skip velocity files if only continuum is needed

    ! Todo : offset in azimuth between fargo3d and mcfost --> just rotate planet coordinates ?

    ! Todo : correct : Omega_p and vphi (r_grid) for unit scaling


    usolarmass = 1.0_dp
    ulength_au = 1.0_dp
    Ggrav_fargo3d = 1.0_dp

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
    utime = sqrt(ulength**3/((Ggrav/Ggrav_fargo3d)*umass))

    udens = umass / ulength**3
    uvelocity = ulength / utime

    ! dimensions are az, r, theta
    allocate(fargo3d_density(fargo3d%nx,fargo3d%ny,fargo3d%nz),fargo3d_vx(fargo3d%nx,fargo3d%ny,fargo3d%nz), &
         fargo3d_vy(fargo3d%nx,fargo3d%ny,fargo3d%nz),fargo3d_vz(fargo3d%nx,fargo3d%ny,fargo3d%nz),stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error when reading fargo3d files')
    fargo3d_density = 0.0_dp ; fargo3d_vx  = 0.0_dp ; fargo3d_vy = 0.0_dp ; fargo3d_vz = 0.0_dp

    ios = 0
    iunit = 1

    ! Reading planet properties
    filename = trim(fargo3d%dir)//"/planet0.dat"
    write(*,*) "Reading "//trim(filename)
    open(unit=iunit, file=filename, status="old", form="formatted", iostat=ios)
    if (ios /= 0) call error("opening fargo3d file:"//trim(filename))

    read(fargo3d%id,*) id
    do while(ios==0)
       read(iunit,*) i, x, y, z, vx, vy, vz, Mp, time, Omega_p
       if (i==id) exit
    enddo
    simu_time = time * utime

    ! Checking units :
    ! Omega_p * uvelocity == 29.78 km/s : OK

    n_etoiles_old = n_etoiles
    n_etoiles = 2 ! Hard coded for new

    if (lfix_star) then
       write(*,*) ""
       write(*,*) "Stellar parameters will not be updated, only the star positions, velocities and masses"
       if (n_etoiles /= n_etoiles_old) call error("Wrong number of stars in mcfost parameter files")
    else
       write(*,*) ""
       write(*,*) "Updating the stellar properties:"
       write(*,*) "There are now", n_etoiles, "stars in the model"

       ! Saving if the accretion rate was forced
       allocate(etoile_old(n_etoiles_old))
       if (allocated(etoile)) then
          etoile_old(:) = etoile(:)
          deallocate(etoile)
       endif
       allocate(etoile(n_etoiles))
       do i=1, min(n_etoiles, n_etoiles_old)
          etoile(i)%force_Mdot = etoile_old(i)%force_Mdot
          etoile(i)%Mdot = etoile_old(i)%Mdot
       enddo
       ! If we have new stars
       do i=n_etoiles_old,n_etoiles
          etoile(i)%force_Mdot = .false.
          etoile(i)%Mdot = 0.
       enddo
       deallocate(etoile_old)
       etoile(:)%find_spectrum = .true.
    endif

    etoile(1)%x = 0_dp ; etoile(1)%y = 0_dp ; etoile(1)%z = 0_dp
    etoile(1)%vx = 0_dp ; etoile(1)%vy = 0_dp ; etoile(1)%vz = 0_dp
    etoile(1)%M = 1_dp * usolarmass

    ! -x and -y as phi is defined differently in fargo3d
    etoile(2)%x = -x * ulength_au
    etoile(2)%y = -y * ulength_au
    etoile(2)%z =  z * ulength_au

    ! -vx and -y as phi is defined differently in fargo3d
    etoile(2)%vx = -vx * uvelocity
    etoile(2)%vy = -vy * uvelocity
    etoile(2)%vz =  vz * uvelocity

    etoile(2)%M = Mp * usolarmass

    do i=1,n_etoiles
       if (etoile(i)%M > 0.013) then
          write(*,*) "Star   #", i, "xyz=", real(etoile(i)%x), real(etoile(i)%y), real(etoile(i)%z), "au, M=", &
               real(etoile(i)%M), "Msun, Mdot=", real(etoile(i)%Mdot), "Msun/yr"
       else
          write(*,*) "Planet #", i, "xyz=", real(etoile(i)%x), real(etoile(i)%y), real(etoile(i)%z), "au, M=", &
               real(etoile(i)%M * GxMsun/GxMjup), "MJup, Mdot=", real(etoile(i)%Mdot), "Msun/yr"
       endif
       if (i>1) write(*,*)  "       distance=", real(sqrt((etoile(i)%x - etoile(1)%x)**2 + &
            (etoile(i)%y - etoile(1)%y)**2 + (etoile(i)%z - etoile(1)%z)**2)), "au"
    enddo
    if (.not.lfix_star) call compute_stellar_parameters()

    if (lplanet_az) then
       which_planet = 1 ! 1 planet for now
       RT_n_az = 1
       RT_az_min = planet_az + atan2(-y, -x) / deg_to_rad
       RT_az_max = RT_az_min
       write(*,*) "Moving planet #", which_planet, "to azimuth =", planet_az
       write(*,*) "WARNING: updating the azimuth to:", RT_az_min
    endif


    !-----------------------------------
    ! Passing data to mcfost
    !-----------------------------------
    write(*,*) "Converting fargo3d files to mcfost ..."
    lvelocity_file = .true.
    vfield_coord = 3 ! spherical

    recl = dp*fargo3d%nx*fargo3d%ny*fargo3d%nz

    file_types(1) = "gasdens"
    file_types(2) = "gasvx"
    file_types(3) = "gasvy"
    file_types(4) = "gasvz"
    do l=1, 4
       filename = trim(fargo3d%dir)//"/"//trim(file_types(l))//trim(trim(fargo3d%id))//".dat"
       write(*,*) "Reading "//trim(filename)
       open(unit=iunit, file=filename, status="old", form="unformatted", iostat=ios, access="direct" , recl=recl)
       if (ios /= 0) call error("opening fargo3d file:"//trim(filename))
       select case(l)
       case(1)
          read(iunit, rec=1, iostat=ios) fargo3d_density
       case(2)
          read(iunit, rec=1, iostat=ios) fargo3d_vx ! vphi
       case(3)
          read(iunit, rec=1, iostat=ios) fargo3d_vy ! vr
       case(4)
          read(iunit, rec=1, iostat=ios) fargo3d_vz ! vtheta
       end select
       if (ios /= 0) call error("reading fargo3d file:"//trim(filename))
       close(iunit)
    enddo

    allocate(vfield3d(n_cells,3), stat=alloc_status)
    if (alloc_status /= 0) call error("memory allocation error fargo3d vfield3d")

    write(*,*) "Constant spatial distribution"

    do i=1, n_rad
       jj= 0
       bz : do j=j_start+1,nz-1 ! 1 extra empty cell in theta on each side
          if (j==0) cycle bz
          jj = jj + 1
          do phik=1, n_az
             icell = cell_map(i,j,phik)

             densite_gaz(icell) = fargo3d_density(phik,i,jj) * udens
             densite_pouss(:,icell) = fargo3d_density(phik,i,jj) * udens

             vfield3d(icell,1)  = fargo3d_vy(phik,i,jj) * uvelocity! vr
             vfield3d(icell,2)  = (fargo3d_vx(phik,i,jj) + r_grid(icell)/ulength_au * Omega_p) * uvelocity ! vphi : planet at r=1
             vfield3d(icell,3)  = fargo3d_vz(phik,i,jj) * uvelocity! vtheta
          enddo ! k
       enddo bz
    enddo ! i
    deallocate(fargo3d_density,fargo3d_vx,fargo3d_vy,fargo3d_vz)

    !write(*,*) maxval(vfield3d(:,1))/1000., maxval(vfield3d(:,2))/1000., maxval(vfield3d(:,3))/1000.
    !write(*,*) etoile(2)%vx/1000., etoile(2)%vy/1000., etoile(2)%vz/1000.
    !stop

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

    write(*,*) "Done"

  end subroutine read_fargo3d_files

end module read_fargo3d
