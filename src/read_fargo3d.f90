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

  stop

  return

end subroutine read_fargo3d_parameters


subroutine read_fargo3_files()

  ! Reading grid

  ! Reading dump file

end subroutine read_fargo3_files

end module read_fargo3d
