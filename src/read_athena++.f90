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
    integer(HID_T) :: file_id, group_id, a_id, type_id

    character(len=32) :: coord, DatasetNames, VariableNames
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

    call hdf_read_attribute(group_id,"", "VariableNames",VariableNames) ! No interface yet for array of strings
    if (trim(VariableNames) /= "rho") call error("mcfost cannot read VariableNames in athena++ dump")

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

    integer :: ierr, alloc_status, n_variables
    integer(HID_T) :: file_id, group_id

    integer, dimension(3) :: block_size

    integer, dimension(:,:), allocatable :: logical_locations
    real, dimension(:,:,:,:,:), allocatable :: data

    real, dimension(:,:,:), allocatable :: rho ! tmp

    integer :: nx1, nx2, nx3, bs1, bs2, bs3
    integer :: iblock, il, jl, kl, iu, ju, ku

    call open_hdf5file(athena%filename,file_id,ierr)
    if (ierr /= 0) call error("cannot open athena HDF file "//trim(athena%filename))

    ! Open main root group "/"
    call open_hdf5group(file_id,'/',group_id,ierr)
    if (ierr /= 0) call error("cannot open athena '/' group")

    call hdf_read_attribute(group_id,"", "NumVariables",n_variables) ! rho, press, vel1, vel2, vel3
    if (n_variables /= 5) call error("mcfost can only read athena++ grids with NumVariables=5 for now")

    call hdf_read_attribute(group_id,"", "MeshBlockSize",block_size)
    bs1 = block_size(1) ; bs2 = block_size(2) ; bs3 = block_size(3)

    write(*,*) bs1, bs2, bs3

    call hdf_read_attribute(group_id,"", "NumMeshBlocks",n_blocks)

    !write(*,*) nx1, nx2, nx3, lx1, lx2, lx3

    ! Allocating Memory
    allocate(logical_locations(3,n_blocks), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ logical_locations')
    allocate(data(bs1,bs2,bs3,n_blocks,n_variables), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ data')

    ! tmp : rho
    allocate(rho(nx1,nx2,nx3), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error athena++ rho')

    !call hdf_read_attribute(group_id,"", "NumCycles",NumCycles) ! I am not sure what that one is



    !write(*,*) lx1,lx2,lx3 ! ok with python


    !write(*,*) "test"

    ! Levels is always if MaxLevel=0, we do not read it
    call hdf_read_dataset(group_id,"LogicalLocations",logical_locations)

    ! compare with python --> ok
    write(*,*) logical_locations(:,1)
    write(*,*) logical_locations(:,2)
    write(*,*) logical_locations(:,3)
    write(*,*) "..."
    write(*,*) logical_locations(:,n_blocks-2)
    write(*,*) logical_locations(:,n_blocks-1)
    write(*,*) logical_locations(:,n_blocks)

    call hdf_read_dataset(group_id,"prim",data)
    ! compare with python raw data --> ok
    !write(*,*) data(1,1,1,1,:)
    !write(*,*) data(2,1,1,1,:)
    !write(*,*) data(1,1,1,2,:)

    call close_hdf5file(file_id,ierr)
    write(*,*) "Athena++ file read sucessfully", ierr

    read(*,*)

    ! Convert to non-raw data
    do iblock=1, n_blocks
       ! Calculate destination indices
       il = logical_locations(1,iblock) * bs1
       jl = logical_locations(2,iblock) * bs2
       kl = logical_locations(3,iblock) * bs3

       iu = il + bs1
       ju = jl + bs2
       ku = kl + bs3

       !write(*,*) iblock, il, jl, kl


       ! Calculate source indices
       ! We do not filter, so it is always 1 to blocksize
       rho(il+1:iu,jl+1:ju,kl+1:ku) = data(:,:,:,iblock,1)


       ! convert python index to fortran
       !rho = data(:,:,:,iblock,1)
    enddo

    write(*,*) rho(1,1,1), rho(1,1,2), rho(1,2,1), rho(2,1,1)


    return

  end subroutine read_athena_model


end module read_athena
