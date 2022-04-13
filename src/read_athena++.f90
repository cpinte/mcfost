module read_athena

  use parametres
  use messages
  use mcfost_env
  use grid
  use cylindrical_grid
  use density
  use stars, only : compute_stellar_parameters
  use hdf5
  use utils_hdf5 ! mentiplay
  use hdf5_utils ! needed to read attributes, need to merge everything
  !use h5lt ! we can use this instead, official hdf5 high level library

  implicit none

contains

  subroutine read_athena_parameters(filename)

    character(len=*), intent(in) :: filename

    integer :: ierr, alloc_status
    integer(HID_T) :: file_id, group_id

    character(len=32) :: coord, DatasetNames
    integer :: MaxLevel
    real(sp) :: time
    integer, dimension(3) :: RootGridSize
    real, dimension(3) :: RootGridX1, RootGridX2, RootGridX3

    ierr = 0
    alloc_status=0

    athena%filename = filename

    call open_hdf5file(filename,file_id,ierr)
    if (ierr /= 0) call error("cannot open athena HDF file "//trim(filename))

    ! Open main root group "/"
    call open_hdf5group(file_id,'/',group_id,ierr)
    if (ierr /= 0) call error("cannot open athena '/' group")

    ! Reading attributes :
    call hdf_read_attribute(group_id,"", "Coordinates",coord)
    if (trim(coord) /= "spherical_polar") call error("mcfost can only read spherical_polar athena++ grids for now")

    call hdf_read_attribute(group_id,"", "DatasetNames",DatasetNames)
    if (trim(DatasetNames) /= "prim") call error("mcfost can only read prim athena++ grids for now")

    call hdf_read_attribute(group_id,"", "MaxLevel",Maxlevel)
    if (MaxLevel > 0) call error("mcfost can only read athena++ grids with MaxLevel=0 for now")

    call hdf_read_attribute(group_id,"", "RootGridSize",RootGridSize) ! nr, ntheta, naz
    athena%nx1=RootGridSize(1)
    athena%nx2=RootGridSize(2)
    athena%nx3=RootGridSize(3)


    call hdf_read_attribute(group_id,"", "RootGridX1",RootGridX1) ! r
    call hdf_read_attribute(group_id,"", "RootGridX2",RootGridX2) ! theta
    call hdf_read_attribute(group_id,"", "RootGridX3",RootGridX3) ! phi

    athena%x1_min=RootGridX1(1)
    athena%x2_min=RootGridX2(1)
    athena%x3_min=RootGridX3(1)

    athena%x1_max=RootGridX1(2)
    athena%x2_max=RootGridX2(2)
    athena%x3_max=RootGridX3(2)

    athena%log_spacing = (RootGridX1(3) > 1.0)
    athena%corrotating_frame = .true.

    call hdf_read_attribute(group_id,"", "Time",time)
    athena%time = time

!    call hdf_read_attribute(group_id,"", "VariableNames",VariableNames) ! No interface yet for array of strings
!    if (trim(VariableNames) /= "rho") call error("mcfost cannot read VariableNames in athena++ dump")

    ! Updating mcfost parameters
    grid_type = 2
    n_rad = athena%nx1
    n_rad_in = 1
    n_az = athena%nx3
    nz = athena%nx2/2 + 1
    lregular_theta = .true.
    theta_max = 0.5 * pi - athena%x2_min

    if (lscale_length_units) then
       write(*,*) 'Lengths are rescaled by ', real(scale_length_units_factor)
    else
       scale_length_units_factor = 1.0
    endif

    disk_zone(1)%rin  = athena%x1_min * scale_length_units_factor
    disk_zone(1)%edge=0.0
    disk_zone(1)%rmin = disk_zone(1)%rin
    disk_zone(1)%rout = athena%x1_max * scale_length_units_factor
    disk_zone(1)%rmax = disk_zone(1)%rout

    write(*,*) "n_rad=", n_rad, "nz=", nz, "n_az=", n_az
    write(*,*) "rin=", real(disk_zone(1)%rin), "rout=", real(disk_zone(1)%rout)

    return

  end subroutine read_athena_parameters

  !---------------------------------------------

  subroutine read_athena_model()

    ! athena data is ordered in x1 = r, x2 = theta, x3 = phi

    integer :: n_blocks

    integer :: ierr, alloc_status, n_variables, n_etoiles_old
    integer(HID_T) :: file_id, group_id

    integer, dimension(3) :: block_size

    integer, dimension(:,:), allocatable :: logical_locations
    real, dimension(:,:,:,:,:), allocatable :: data

    real, dimension(:,:,:), allocatable :: rho, vx1, vx2, vx3 ! todo : we can save memory and only data to directly pass it to mcfost

    integer :: nx1, nx2, nx3, bs1, bs2, bs3
    integer :: i, iblock, il, jl, kl, iu, ju, ku, j, jj, phik, icell

    real(dp) :: Ggrav_athena, umass, usolarmass, ulength, utime, udens, uvelocity, ulength_au, mass, facteur
    type(star_type), dimension(:), allocatable :: etoile_old

    ! Planet properties hard coded for now
    real, parameter :: Mp = 1e-3
    real, parameter :: Omega_p = 1.0
    real, parameter :: x = 1.0, y=0.0, z=0.0
    real, parameter :: vx=0.0, vy=1.0, vz=1.0

    usolarmass = 1.0_dp
    ulength_au = 1.0_dp
    Ggrav_athena = 1.0_dp

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
    utime = sqrt(ulength**3/((Ggrav/Ggrav_athena)*umass))

    udens = umass / ulength**3
    uvelocity = ulength / utime

    simu_time = athena%time * utime

    ! --- copy/paste from read_fargo3d
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
    ! ---- end copy/paste from read_fargo3d

    write(*,*) "Reading "//trim(athena%filename)
    call open_hdf5file(athena%filename,file_id,ierr)
    if (ierr /= 0) call error("cannot open athena HDF file "//trim(athena%filename))

    ! Open main root group "/"
    call open_hdf5group(file_id,'/',group_id,ierr)
    if (ierr /= 0) call error("cannot open athena '/' group")

    call hdf_read_attribute(group_id,"", "NumVariables",n_variables) ! rho, press, vel1, vel2, vel3
    if (n_variables /= 5) call error("mcfost can only read athena++ grids with NumVariables=5 for now")

    call hdf_read_attribute(group_id,"", "MeshBlockSize",block_size)
    bs1 = block_size(1) ; bs2 = block_size(2) ; bs3 = block_size(3)

    call hdf_read_attribute(group_id,"", "NumMeshBlocks",n_blocks)

    ! Allocating Memory
    allocate(logical_locations(3,n_blocks), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ logical_locations')
    allocate(data(bs1,bs2,bs3,n_blocks,n_variables), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ data')

    nx1 = athena%nx1 ; nx2 = athena%nx2 ; nx3 = athena%nx3
    allocate(rho(nx1,nx2,nx3), vx1(nx1,nx2,nx3), vx2(nx1,nx2,nx3), vx3(nx1,nx2,nx3), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ rho')

    ! Levels is always if MaxLevel=0, we do not read it
    call hdf_read_dataset(group_id,"LogicalLocations",logical_locations)

    ! compare with python --> ok
    !write(*,*) logical_locations(:,1)
    !write(*,*) logical_locations(:,2)
    !write(*,*) logical_locations(:,3)
    !write(*,*) "..."
    !write(*,*) logical_locations(:,n_blocks-2)
    !write(*,*) logical_locations(:,n_blocks-1)
    !write(*,*) logical_locations(:,n_blocks)

    call hdf_read_dataset(group_id,"prim",data)
    ! compare with python raw data --> ok
    !write(*,*) data(1,1,1,1,:)
    !write(*,*) data(2,1,1,1,:)
    !write(*,*) data(1,1,1,2,:)

    call close_hdf5file(file_id,ierr)
    write(*,*) "Athena++ file read sucessfully"

    ! Convert to non-raw data, ie merge all blocks
    do iblock=1, n_blocks
       ! Calculate destination indices
       il = logical_locations(1,iblock) * bs1
       jl = logical_locations(2,iblock) * bs2
       kl = logical_locations(3,iblock) * bs3

       iu = il + bs1
       ju = jl + bs2
       ku = kl + bs3

       rho(il+1:iu,jl+1:ju,kl+1:ku) = data(:,:,:,iblock,1)
       vx1(il+1:iu,jl+1:ju,kl+1:ku) = data(:,:,:,iblock,3) ! vr
       vx2(il+1:iu,jl+1:ju,kl+1:ku) = data(:,:,:,iblock,4) ! vtheta
       vx3(il+1:iu,jl+1:ju,kl+1:ku) = data(:,:,:,iblock,5) ! vphi
    enddo
    deallocate(data)

    !write(*,*) rho(1,1,1), rho(1,1,2), rho(1,2,1), rho(2,1,1) ! test ok with python
    !write(*,*) rho(101,101,101) ! test ok with python

    !-----------------------------------
    ! Passing data to mcfost
    !-----------------------------------
    write(*,*) "Converting athena++ model to mcfost ..."
    lvelocity_file = .true.
    vfield_coord = 3 ! spherical

    allocate(vfield3d(n_cells,3), stat=alloc_status)
    if (alloc_status /= 0) call error("memory allocation error athena++ vfield3d")

    write(*,*) "Constant spatial distribution"

    do i=1, n_rad
       jj= 0
       bz : do j=j_start+1,nz-1 ! 1 extra empty cell in theta on each side
          if (j==0) cycle bz
          jj = jj + 1
          do phik=1, n_az
             icell = cell_map(i,j,phik)

             densite_gaz(icell) =  rho(i,jj,phik) * udens
             densite_pouss(:,icell) = rho(i,jj,phik) * udens

             vfield3d(icell,1)  = vx1(i,jj,phik) * uvelocity! vr
             vfield3d(icell,2)  = (vx3(i,jj,phik) + r_grid(icell)/ulength_au * Omega_p) * uvelocity ! vphi : planet at r=1
             vfield3d(icell,3)  = vx2(i,jj,phik) * uvelocity! vtheta
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

  end subroutine read_athena_model


end module read_athena
