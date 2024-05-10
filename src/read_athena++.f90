module read_athena

  use mem, only : allocate_densities
  use parametres
  use messages
  use mcfost_env
  use grid
  ! use cylindrical_grid
  ! use spherical_grid
  use density
  use stars, only : compute_stellar_parameters, stars_cell_indices
  use hdf5
  use utils, only : meshgrid_3d, volumegrid_3d, to_cartesian
  use utils_hdf5 ! mentiplay
  use hdf5_utils ! needed to read attributes, need to merge everything
  use arb2mcfost
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
    athena%arb_grid = .true.

    call open_hdf5file(filename,file_id,ierr)
    if (ierr /= 0) call error("cannot open athena HDF file "//trim(filename))

    ! Open main root group "/"
    call open_hdf5group(file_id,'/',group_id,ierr)
    if (ierr /= 0) call error("cannot open athena '/' group")

    ! Reading attributes :
    call hdf_read_attribute(group_id,"", "Coordinates",coord)
    write(*,*) "Coordinates is ", coord

    if (trim(coord) == "cartesian") then
      athena%coord = 0
    else if (trim(coord) == "cylindrical") then
      athena%coord = 1
    else if (trim(coord) == "spherical_polar") then
      athena%coord = 2
    else
      call error("Unknown athena++ grid. Exiting")
    endif


    ! call hdf_read_attribute(group_id,"", "DatasetNames",DatasetNames)
    ! if (trim(DatasetNames) /= "prim") call error("mcfost can only read prim athena++ grids for now")

    athena%maxlevel = 0
    call hdf_read_attribute(group_id,"", "MaxLevel",Maxlevel)
    if (MaxLevel > 0)  then
      athena%maxlevel = MaxLevel
      write(*,*) "Mesh refinment detected in athena++ model. "
      ! call error("mcfost can only read athena++ grids with MaxLevel=0 for now")
    else
      if (athena%coord == 1 .or. athena%coord == 2) then
        athena%arb_grid = .false.
      endif
    endif

    ! nr, ntheta, naz
    ! x, y, z
    ! nr, ntheta, nphi
    call hdf_read_attribute(group_id,"", "RootGridSize",RootGridSize)
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

    if (athena%arb_grid) then
      return

    else
      ! Updating mcfost parameters
      ! call warning("athena : forcing spherical grid") ! only spherical grid is implemented for now
      grid_type = athena%coord
      n_rad = athena%nx1
      n_rad_in = 1
      if (grid_type == 3) then
        n_az = athena%nx3
        nz = athena%nx2/2 + 1
      else if (grid_type == 2) then
        n_az = athena%nx2
        nz = athena%nx3/2 + 1
      endif

      lregular_theta = .true.
      theta_max = 0.5 * pi - athena%x2_min

      if (lscale_length_units) then
         write(*,*) 'Lengths are rescaled by ', real(scale_length_units_factor)
      else
         scale_length_units_factor = 1.0
      endif

      disk_zone(1)%geometry = grid_type
      disk_zone(1)%rin  = athena%x1_min * scale_length_units_factor
      disk_zone(1)%edge = 0.0
      disk_zone(1)%rmin = disk_zone(1)%rin
      disk_zone(1)%rout = athena%x1_max * scale_length_units_factor
      disk_zone(1)%rmax = disk_zone(1)%rout

      write(*,*) "n_rad=", n_rad, "nz=", nz, "n_az=", n_az
      write(*,*) "rin=", real(disk_zone(1)%rin), "rout=", real(disk_zone(1)%rout)

      return
    end if

  end subroutine read_athena_parameters

  !---------------------------------------------

  subroutine read_athena_model()

    ! athena data is ordered in x1 = r, x2 = theta, x3 = phi

    integer :: n_blocks

    integer :: ierr, alloc_status, n_variables, n_etoiles_old
    integer(HID_T) :: file_id, group_id

    integer, dimension(3) :: block_size

    integer, dimension(:), allocatable :: n_variables_1
    integer, dimension(:,:), allocatable :: logical_locations
    real, dimension(:,:,:,:,:), allocatable :: data

    real(kind=dp), dimension(:,:,:), allocatable :: rho, vx1, vx2, vx3 ! todo : we can save memory and only data to directly pass it to mcfost
    real(kind=dp), dimension(:,:,:), allocatable :: rho_tmp, vx1_tmp, vx2_tmp, vx3_tmp, x1_tmp, x2_tmp, x3_tmp, v_tmp
    real(kind=dp), dimension(:), allocatable :: rho_a, vx1_a, vx2_a, vx3_a, x1_a, x2_a, x3_a, v_a, x1f_tmp, x2f_tmp, x3f_tmp ! For the arbitrary grids where we need position to be passed to voronoi
    real(kind=dp), dimension(:), allocatable :: xx, yy, zz, vxx, vyy, vzz, mass_gas, h
    integer, dimension(:), allocatable :: particle_id
    real, dimension(:,:), allocatable :: x1f, x2f, x3f, x1v, x2v, x3v

    integer :: nx1, nx2, nx3, bs1, bs2, bs3
    integer :: i, iblock, il, jl, kl, iu, ju, ku, j, jj, phik, icell, it, jt, kt

    real(dp) :: Ggrav_athena, umass, usolarmass, ulength, utime, udens, uvelocity, ulength_au, mass, facteur
    type(star_type), dimension(:), allocatable :: etoile_old

    ! Planet properties hard coded for now
    real, parameter :: Mp = 1e-3
    real, parameter :: Omega_p = 1.0
    real, parameter :: x = 6.0, y=0.0, z=0.0
    real, parameter :: vx=0.0, vy=1.0, vz=1.0
    logical :: print_messages

    ! print_messages = .true.
    ! call hdf_set_print_messages(print_messages)

    usolarmass = 1.0_dp
    ulength_au = 1.0_dp
    Ggrav_athena = 1.0_dp

    if (lscale_length_units) then
       write(*,*) 'Lengths are rescaled by ', real(scale_length_units_factor)
       ulength_au = ulength_au * scale_length_units_factor
    endif

    if (lscale_mass_units) then
       write(*,*) 'Masses are rescaled by ', real(scale_mass_units_factor)
       usolarmass = usolarmass * scale_mass_units_factor
    endif

    umass = usolarmass *  Msun_to_kg
    ulength = ulength_au * AU_to_m
    utime = sqrt(ulength**3/((Ggrav/Ggrav_athena)*umass))

    udens = umass / ulength**3
    uvelocity = ulength / utime

    simu_time = athena%time * utime

    ! --- copy/paste from read_fargo3d
    n_etoiles_old = n_etoiles
    n_etoiles = 1 ! Hard coded for now

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

    ! ! -x and -y as phi is defined differently in fargo3d
    ! etoile(2)%x = -x * ulength_au
    ! etoile(2)%y = -y * ulength_au
    ! etoile(2)%z =  z * ulength_au
    !
    ! ! -vx and -y as phi is defined differently in fargo3d
    ! etoile(2)%vx = -vx * uvelocity
    ! etoile(2)%vy = -vy * uvelocity
    ! etoile(2)%vz =  vz * uvelocity
    !
    ! etoile(2)%M = Mp * usolarmass

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

    if (athena%arb_grid) then
      lVoronoi = .true.
      nx1 = n_blocks * bs1 * bs2 * bs3
      allocate(rho_a(nx1), vx1_a(nx1), vx2_a(nx1), vx3_a(nx1), stat=alloc_status)
      allocate(x1_a(nx1), x2_a(nx1), x3_a(nx1), v_a(nx1), stat=alloc_status)
      ! allocate(x1f(n_blocks, bs1+1), x2f(n_blocks, bs2+1), x3f(n_blocks, bs3+1), \
      !          x1v(n_blocks, bs1), x2v(n_blocks, bs2), x3v(n_blocks, bs3), stat=alloc_status)
      allocate(x1f(bs1+1, n_blocks), x2f(bs2+1, n_blocks), x3f(bs3+1, n_blocks), stat=alloc_status)
      allocate(x1v(bs1, n_blocks), x2v(bs2, n_blocks), x3v(bs3, n_blocks), stat=alloc_status)
      allocate(x1_tmp(bs1, bs2, bs3), x2_tmp(bs1, bs2, bs3), x3_tmp(bs1, bs2, bs3), v_tmp(bs1, bs2, bs3), stat=alloc_status)
      allocate(x1f_tmp(bs1), x2f_tmp(bs2), x3f_tmp(bs3), stat=alloc_status)
      allocate(particle_id(nx1), stat=alloc_status)
      particle_id(:) = 0

      call hdf_read_dataset(group_id,"x1f",x1f)
      call hdf_read_dataset(group_id,"x2f",x2f)
      call hdf_read_dataset(group_id,"x3f",x3f)
      call hdf_read_dataset(group_id,"x1v",x1v)
      call hdf_read_dataset(group_id,"x2v",x2v)
      call hdf_read_dataset(group_id,"x3v",x3v)
    else
      nx1 = athena%nx1 ; nx2 = athena%nx2 ; nx3 = athena%nx3
      allocate(rho(nx1,nx2,nx3), vx1(nx1,nx2,nx3), vx2(nx1,nx2,nx3), vx3(nx1,nx2,nx3), stat=alloc_status)
    endif
    if (alloc_status > 0) call error('Allocation error athena++ rho')


    ! Levels is always if MaxLevel=0, we do not read it
    call hdf_read_dataset(group_id,"LogicalLocations",logical_locations)

    call hdf_read_dataset(group_id,"prim",data)

    call close_hdf5file(file_id,ierr)
    write(*,*) "Athena++ file read sucessfully"

    write(*,*) "logical_locations", shape(logical_locations)

    ! allocate(tmp_flatten(bs1*bs2*bs3), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ rho')

    if (athena%arb_grid) then
      ! Now... What do we do here...
      i = 1
      do iblock=1, n_blocks
         iu = iblock*bs1*bs2*bs3
         il = iu - bs1*bs2*bs3 + 1

         rho_a(il:iu) = reshape(data(:,:,:,iblock,1), (/size(data(:,:,:,iblock,1))/) )
         vx1_a(il:iu) = reshape(data(:,:,:,iblock,2), (/size(data(:,:,:,iblock,2))/) )
         vx2_a(il:iu) = reshape(data(:,:,:,iblock,3), (/size(data(:,:,:,iblock,3))/) )
         vx3_a(il:iu) = reshape(data(:,:,:,iblock,4), (/size(data(:,:,:,iblock,4))/) )

         ! call meshgrid_3d(x1v(iblock, :), x2v(iblock, :), x3v(iblock, :), x1_tmp, x2_tmp, x3_tmp)  ! (x, y, z, xx, yy, zz)
         call meshgrid_3d(x1v(:, iblock), x2v(:, iblock), x3v(:, iblock), x1_tmp, x2_tmp, x3_tmp)  ! (x, y, z, xx, yy, zz)

         ! write(*,*) "outside of meshgrid_3d"

         x1_a(il:iu) = reshape(x1_tmp, (/size(x1_tmp)/) )
         x2_a(il:iu) = reshape(x2_tmp, (/size(x2_tmp)/) )
         x3_a(il:iu) = reshape(x3_tmp, (/size(x3_tmp)/) )

         call volumegrid_3d(x1f(:, iblock), x2f(:, iblock), x3f(:, iblock), v_tmp, coord=athena%coord)

         v_a(il:iu) = reshape(v_tmp, (/size(v_tmp)/) )

      enddo
      write(*,*) "Data successfully reshaped and read "
      deallocate(data, x1v, x2v, x3v, x1_tmp, x2_tmp, x3_tmp, v_tmp, x1f, x2f, x3f)

      ! Need to convert from density to mass
      mass_gas = rho_a*udens*v_a !* AU3_to_m3  * g_to_Msun
      ! write(*,*) "AU3_to_m3 * g_to_Msun", AU3_to_m3 * g_to_Msun
      ! write(*,*) "masse_mol_gaz", masse_mol_gaz
      write(*,*) "udens", udens
      write(*,*) "uvelocity", uvelocity
      write(*,*) "ulength", ulength


      ! now that v_a is no longer needed, it will become h
      h = v_a**1/3

      ! Convert coordinates and velocities to Cartesian if necessary
      if (.not. athena%coord==0) then
        ! coordinates
        allocate(xx(size(x1_a)), yy(size(x1_a)), zz(size(x1_a)))
        call to_cartesian(x1_a, x2_a, x3_a, xx, yy, zz, athena%coord)

        ! velocites
        allocate(vxx(size(vx1_a)), vyy(size(vx1_a)), vzz(size(vx1_a)))
        call to_cartesian_velocities(vx1_a, vx2_a, vx3_a, vxx, vyy, vzz, x1_a, x2_a, x3_a, athena%coord)
        vfield_coord = 1
        ! call to_cartesian(vx1_a, vx2_a, vx3_a, vx, vy, vz, athena%coord)
        write(*,*) "Data successfully converted to cartesian "
      else
        ! Already in cartesian coordinates
        xx = x1_a
        yy = x2_a
        zz = x3_a

        vxx = vx1_a!*uvelocity
        vyy = vx2_a!*uvelocity
        vzz = vx3_a!*uvelocity
      endif

      deallocate(x1_a, x2_a, x3_a, vx1_a, vx2_a, vx3_a, rho_a, v_a)

      ! do i=1, 10
      !   write(*,*) "h", h(i), "x", xx(i), "y", yy(i), "z", zz(i)
      ! enddo
      !
      ! write(*,*) "second lot"

      ! do i=115480, 115490
      !   write(*,*) "h", h(i), "x", xx(i), "y", yy(i), "z", zz(i), "vxx", vxx(i), "vyy", vyy(i),  "vzz", vzz(i)
      ! enddo

      ! write(*,*) "third lot"
      !
      ! do i=3900000, 3900010
      !   write(*,*) "h", h(i), "x", xx(i), "y", yy(i), "z", zz(i)
      ! enddo

      vxx = vxx*uvelocity
      vyy = vyy*uvelocity
      vzz = vzz*uvelocity

      call setup_arb_to_mcfost(xx, yy, zz, h, vxx, vyy, vzz, mass_gas, particle_id)

      return
    else
      call setup_grid()
      call define_grid()
      call stars_cell_indices()

      call allocate_densities()

      ! Convert to non-raw data, ie merge all blocks
      do iblock=1, n_blocks
         ! Calculate destination indices
         il = logical_locations(1,iblock) * bs1
         jl = logical_locations(2,iblock) * bs2
         kl = logical_locations(3,iblock) * bs3
         write(*,*) il, jl, kl

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

      allocate(vfield3d(n_cells,3), stat=alloc_status)
      if (alloc_status /= 0) call error("memory allocation error athena++ vfield3d")

      write(*,*) "Constant spatial distribution"

      ! Athena cartesian: x, y, z,
      ! Athena cylindrical: r, phi, z
      ! Athena spherical: r, theta, phi
      ! MCFOST: 1 = cartesian, 2 = cylindrical, 3 = spherical
      ! TO DO: Add parameter for co-rotating frame
      vfield_coord = disk_zone(1)%geometry + 1
      if (vfield_coord == 3) then
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
                 ! I guess the below line is only true if the simulation is done in a co-rotating frame
                 vfield3d(icell,2)  = (vx3(i,jj,phik) + r_grid(icell)/ulength_au * Omega_p) * uvelocity ! vphi : planet at r=1
                 ! vfield3d(icell,2)  = vx3(i,jj,phik) * uvelocity ! vphi : planet at r=1
                 vfield3d(icell,3)  = vx2(i,jj,phik) * uvelocity! vtheta
              enddo ! k
           enddo bz
        enddo ! i
      else if (vfield_coord == 2) then
        do i=1, n_rad
           ! jj= 0
           do j=j_start+1,nz-1 ! 1 extra empty cell in theta on each side
              jj = jj + 1
              do phik=1, n_az
                 icell = cell_map(i,j,phik)

                 densite_gaz(icell) =  rho(i,jj,phik) * udens
                 densite_pouss(:,icell) = rho(i,jj,phik) * udens

                 vfield3d(icell,1)  = vx1(i,jj,phik) * uvelocity ! vr
                 vfield3d(icell,2)  = vx2(i,jj,phik) * uvelocity ! vphi
                 vfield3d(icell,3)  = vx3(i,jj,phik) * uvelocity ! vz
              enddo ! k
           enddo
        enddo ! i
      endif

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
    endif

  end subroutine read_athena_model


end module read_athena
