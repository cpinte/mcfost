module read_fargo3d

  use parametres
  use messages
  use mcfost_env

  implicit none

  type fargo3d_model
     integer :: nx, ny, nz, realtype
     real(kind=dp) :: xmin,xmax, ymin,ymax, zmin,zmax
     logical :: log_spacing, corrotating_frame

     real(kind=dp) :: dt, aspect_ratio, nu, gamma, cs
     character(len=128) :: dir, id, planetconfig

  end type fargo3d_model

  type(fargo3d_model) :: fargo3d

  contains

subroutine read_fargo3d_parameters(dir,id)

  character(len=*), intent(in) :: dir, id

  character(len=128) :: s, key, val, filename
  integer :: ios, n_lines, iunit

  fargo3d%dir = dir
  fargo3d%id = id

  write(*,*) fargo3d%dir

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

     write(*,*) trim(key), " --> ", trim(val)

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
           call error("MCFOST linear grid for fargo3d not implemented yet")
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
  write(*,*) "Reading FARGO3D model"
  write(*,*) "Forcing spherical grid and updating dimension"

  grid_type = 2
  n_rad = fargo3d%ny
  n_rad_in = 1
  n_az = fargo3d%nx
  nz = fargo3d%nz/2 + 1

  disk_zone(1)%rin  = fargo3d%ymin
  disk_zone(1)%edge=0.0
  disk_zone(1)%rmin = disk_zone(1)%rin
  disk_zone(1)%rout = fargo3d%ymax
  disk_zone(1)%rmax = fargo3d%ymax

  write(*,*) "n_rad=", n_rad, "nz=", nz, "n_az=", n_az
  write(*,*) "rin=", real(disk_zone(1)%rin), "rout=", real(disk_zone(1)%rout)

  return

end subroutine read_fargo3d_parameters


!---------------------------------------------

subroutine check_fargo3d_grid(r,theta,phi)

  real(dp), dimension(*) :: r, theta, phi

  return

end subroutine check_fargo3d_grid

!---------------------------------------------


subroutine read_fargo3d_files()

  ! Reading dump file
  !real*8  :: data(nx*ny*nz)
  !open(unit=100, status="old", file=filename, form="unformatted", access="direct", recl = NX*NY*NZ*8)
  !read(100,rec=1) data

  real(dp), dimension(:,:,:), allocatable  :: data
  integer :: ios, iunit, alloc_status, i

  character(len=128) :: filename
  character(len=16), dimension(4) :: file_types

  allocate(data(fargo3d%nx,fargo3d%ny,fargo3d%nz), stat=alloc_status) ! ie --> az, r, theta
  if (alloc_status > 0) call error('Allocation error when reading fargo3d files')
  data = 0.0_dp

  ios = 0
  iunit = 1

  file_types(1) = "gasdens"
  file_types(2) = "gasvx"
  file_types(3) = "gasvy"
  file_types(4) = "gasvz"
  do i=1, 4
     filename = trim(fargo3d%dir)//"/"//trim(file_types(i))//trim(trim(fargo3d%id))//".dat"
     write(*,*) "Reading "//trim(filename)
     open(unit=iunit, status="old", file=filename, form="unformatted", iostat=ios, access="direct", recl=8*fargo3d%nx*fargo3d%ny*fargo3d%nz)
     if (ios /= 0) call error("opening fargo3d file:"//trim(filename))
     read(iunit, rec=1, iostat=ios) data
     if (ios /= 0) call error("reading fargo3d file:"//trim(filename))
     close(iunit)


  enddo

  deallocate(data)
  stop

end subroutine read_fargo3d_files



end module read_fargo3d
