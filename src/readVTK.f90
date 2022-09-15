module VTKtools

  use messages, only : error

  implicit none

contains

  function getLineSize(file, position)

    integer,intent(in) :: file, position
    integer :: getLineSize

    character :: c
    integer :: loc, s

    loc = 0
    do
       read(file,pos=position+loc,iostat=s) c
       if(s<0) then
          !EOF (presumably)
          loc = -1
          exit
       endif
       ! end of line!
       if(iachar(c) == 10) then
          exit
       endif
       loc = loc+1
    enddo
    getLineSize = loc

    return

  end function getLineSize

  subroutine readLine(file, position, line, newposition)

    integer, intent(in) :: file, position
    character (:), intent(out), allocatable :: line
    integer, intent(out) :: newposition

    integer size
    size = getLineSize(file, position)
    if(size<0) then
       ! eof
       newposition = -1
       return
    endif
    allocate(character(size) :: line)
    read(file, pos=position) line
    newposition = position+size+1

    return

  end subroutine readLine

  subroutine readField(file, oldposition, geometry, periodicity, time, newposition)

    integer, intent(in) :: file, oldposition
    integer, intent(out) :: geometry
    integer, dimension(3), intent(out) :: periodicity
    integer geo1, geo2
    real, intent(out) :: time
    integer, intent(out) :: newposition

    integer position, n
    character (:), allocatable :: line
    integer nfield, lineSize

    position = oldposition
    ! Read the number of FIELD in the header
    call readLine(file, position, line, newposition)
    lineSize = newposition-position-1
    read(line(lineSize-1:lineSize),*) nfield

    position = newposition
    ! read each entry
    do n = 1, nfield
       call readLine(file, position, line, newposition)
       if(line(1:8) == "GEOMETRY") then
          read(file,pos=newposition) geometry
          ! update file pointer, and add carriage return
          inquire(file, pos=position)
          position=position+1
       else if(line(1:11) == "PERIODICITY") then
          read(file,pos=newposition) periodicity
          ! update file pointer, and add carriage return
          inquire(file, pos=position)
          position=position+1
       else if(line(1:4) == "TIME") then
          read(file,pos=newposition) time
          ! update file pointer, and add carriage return
          inquire(file, pos=position)
          position=position+1
       else
          print *, "Unknown entry: ", line
       endif
    enddo
    newposition = position

    return

  end subroutine readField

  subroutine readPoints(file, oldposition, geometry, dimensions, x1, x2, x3, newposition)

    integer, intent(in) :: file, oldposition, geometry
    integer, dimension(3), intent(in) :: dimensions
    real, intent(out), allocatable :: x1(:), x2(:), x3(:)
    integer, intent(out) :: newposition
    character (:), allocatable :: line

    real, allocatable :: pts(:)
    real, dimension(:,:,:), allocatable :: x, y, z
    real q1, q2;

    integer position, lineSize, nPoints
    integer i,j,k,n, nx1mid, nx2mid


    position = oldposition

    call readLine(file, position, line, newposition)

    lineSize = newposition-position-1
    read(line(7:lineSize),*) nPoints

    position=newposition
    allocate(pts(3*npoints))
    allocate(x(dimensions(1),dimensions(2),dimensions(3)))
    allocate(y(dimensions(1),dimensions(2),dimensions(3)))
    allocate(z(dimensions(1),dimensions(2),dimensions(3)))

    read(file,pos=newposition) pts
    ! update file pointer, and add carriage return
    inquire(file, pos=position)
    position=position+1
    newposition=position

    ! load cartesian points
    do k = 1, dimensions(3)
       do j = 1, dimensions(2)
          do i = 1, dimensions(1)
             n= i+(j-1)*dimensions(1)+(k-1)*dimensions(1)*dimensions(2)
             x(i,j,k) = pts(3*(n-1)+1)
             y(i,j,k) = pts(3*(n-1)+2)
             z(i,j,k) = pts(3*(n-1)+3)

             !print *,"n=",n,"i=",i,"j=",j,"k=",k,"x=",x(i,j,k),"y=",y(i,j,k), "z=",z(i,j,k)
          enddo
       enddo
    enddo

    ! convert cartesian to dedicated geometry
    if(geometry==2) then
       ! spherical geometry
       allocate(x1(dimensions(1)))
       do i = 1, dimensions(1)
          x1(i) = sqrt(x(i,1,1)**2 + y(i,1,1)**2 + z(i,1,1)**2)
       enddo

       nx1mid = (dimensions(1)-1)/2
       nx2mid = (dimensions(2)-1)/2
       allocate(x2(dimensions(2)))
       do i = 1, dimensions(2)
          x2(i) = acos( z(1,i,1) / sqrt(x(1,i,1)**2+y(1,i,1)**2+z(1,i,1)**2))
       enddo

       allocate(x3(dimensions(3)))
       do i = 1, dimensions(3)
          x3(i) = atan2(y(nx1mid,nx2mid,i), x(nx1mid,nx2mid,i))
          if(x3(i)< 0.0) then
             x3(i) = x3(i) + 8.0*atan(1.0)
          endif
       enddo
    else
       print *, "Unable to transform this geometry"
    endif

    return

  end subroutine readPoints

  subroutine readScalars(file, oldposition, dimensions, array, name, newposition)

    integer, intent(in) :: file, oldposition
    integer, dimension(3), intent(in) :: dimensions
    real, dimension(:,:,:), intent(out), allocatable :: array
    character (3), intent(out) :: name
    integer, intent(out) :: newposition


    character (:), allocatable :: line
    integer position, lineSize

    position = oldposition

    ! read the name of the scalar
    call readLine(file, position, line, newposition)
    lineSize = newposition-position-1
    read(line(9:11),*) name
    position=newposition
    ! read the lookup table
    call readLine(file, position, line, position)
    ! read the meat
    allocate(array(dimensions(1)-1,dimensions(2)-1,dimensions(3)-1))
    read(file,pos=position) array
    ! update file pointer, and add carriage return
    inquire(file, pos=position)
    position=position+1
    newposition=position
  end subroutine readScalars

  subroutine readVTK_header(filename, unit, position, dimensions, time, origin, x1, x2, x3)
    ! read only the header of a VTK file and keepthe unit open

    character(len=*), intent(in) :: filename
    integer, intent(out) :: unit
    integer, dimension(3), intent(out) :: dimensions
    real, intent(out) :: time
    character(:), allocatable, intent(out) :: origin
    real, dimension(:), allocatable, intent(out) :: x1, x2, x3

    character (:), allocatable :: line
    character (3) :: varName
    integer :: position, newposition
    integer :: lineSize, nPoints
    integer :: geometry
    integer, dimension(3) :: periodicity

    open(newunit=unit, file=filename,form="unformatted",access="stream",CONVERT='BIG_ENDIAN',status="old")
    position = 1

    ! 1st line: 1 vtk Data......
    call readLine(unit, position, line, position)

    ! 2nd Line : Idefix xxxxx
    call readLine(unit, position, line, position)
    if (line(1:6) /= "Idefix") call error(trim(filename)//" does not seem to be an Idefix file: "//line(1:6))
    origin = line

    ! 3rd Line : BINARY
    call readLine(unit, position, line, position)
    if (line(1:6) /= "BINARY") call error(trim(filename)//" does not seem to be an Idefix fileis not a binary Idefix file")

    ! 4th Line : DATASET STRUCTURED (or other)
    call readLine(unit, position, line, position)

    ! Field (+ geometry, periodicity and time))
    call readLine(unit, position, line, newposition)
    if (line(1:5) /= "FIELD") call error("FIELD not found in vtk file: "//trim(filename))
    call readField(unit, position, geometry, periodicity, time, position)

    ! Dimensions
    call readLine(unit, position, line, newposition)
    if (line(1:10) /= "DIMENSIONS") call error("reading DIMENSIONS in vtk file: "//trim(filename))
    lineSize = newposition-position-1
    read(line(11:lineSize),*) dimensions
    position=newposition

    ! Points
    call readLine(unit, position, line, newposition)
    if (line(1:6) /= "POINTS") call error("reading POINTS in vtk file: "//trim(filename))
    call readPoints(unit, position, geometry, dimensions, x1, x2, x3, position)

    return

  end subroutine readVTK_header


  subroutine readVTK_data(unit, old_position, dimensions, rho, vx1, vx2, vx3)
    ! read only the data of a VTK file (after readVTK_header)

    integer, intent(in) :: unit, old_position
    integer, dimension(3), intent(in) :: dimensions
    real, allocatable, dimension(:,:,:), intent(out) :: rho, vx1, vx2, vx3


    integer :: position, newposition
    character (:), allocatable :: line
    character (3) :: varName
    real, allocatable, dimension(:,:,:) :: array

    position = old_position
    do
       call readLine(unit, position, line, newposition)
       if(newposition<0) exit
       if(line(1:9) == "CELL_DATA") then
          ! we're entering the meat
          !skip the extra line feed
          position=newposition+1
       else if(line(1:7) == "SCALARS") then
          call readScalars(unit, position, dimensions, array, varName, position)
          if(varName=="RHO") then
             rho = array
          else if(varName=="VX1") then
             vx1 = array
          else if(varName=="VX2") then
             vx2 = array
          else if(varName=="VX3") then
             vx3 = array
          endif
       else
          write(*,*) "Unknown line in file:", line
          exit
       endif
    enddo

    close(unit=unit)

  end subroutine readVTK_data


!  subroutine readVTK(filename, dimensions, time, origin, x1, x2, x3, rho, vx1, vx2, vx3)
!    implicit none
!    character(len=*), intent(in) :: filename
!    integer, dimension(3), intent(out) :: dimensions
!    real, intent(out) :: time
!    character (:), allocatable, intent(out) :: origin
!    real, allocatable, intent(out) :: x1(:), x2(:), x3(:)
!    real, allocatable, dimension(:,:,:), intent(out) :: rho, vx1, vx2, vx3
!
!    character (:), allocatable :: line
!    character (3) :: varName
!    integer :: position, newposition
!    integer :: lineSize, nPoints
!    integer :: geometry
!    integer, dimension(3) :: periodicity
!
!    real, allocatable, dimension(:,:,:) :: array
!
!
!    ! This is a comment line; it is ignored by the compiler
!    open(unit = 1, file= filename,form="unformatted",access="stream",CONVERT='BIG_ENDIAN')
!    position = 1
!
!    ! 1st line: 1 vtk Data......
!    call readLine(1, position, line, position)
!
!    ! 2nd Line : Idefix xxxxx
!    call readLine(1, position, line, position)
!    origin = line
!
!    ! 3rd Line : BINARY (should be checked!)
!    call readLine(1, position, line, position)
!
!    ! 4th Line : DATASET STRUCTURED (or other)
!    call readLine(1, position, line, position)
!
!    !Field data and others
!    do
!       call readLine(1, position, line, newposition)
!       if(newposition<0) exit
!       if(line(1:5) == "FIELD") then
!          call readField(1, position, geometry, periodicity, time, position)
!       else if(line(1:10) == "DIMENSIONS") then
!          lineSize = newposition-position-1
!          read(line(11:lineSize),*) dimensions
!          position=newposition
!       else if(line(1:6) == "POINTS") then
!          call readPoints(1, position, geometry, dimensions, x1, x2, x3, position)
!       else if(line(1:9) == "CELL_DATA") then
!          ! we're entering the meat
!          !skip the extra line feed
!          position=newposition+1
!       else if(line(1:7) == "SCALARS") then
!          call readScalars(1, position, dimensions, array, varName, position)
!          if(varName=="RHO") then
!             rho = array
!          else if(varName=="VX1") then
!             vx1 = array
!          else if(varName=="VX2") then
!             vx2 = array
!          else if(varName=="VX3") then
!             vx3 = array
!          endif
!       else
!          print *, "Unknown line in file:", line
!          exit
!       endif
!    enddo
!    close(1)
!
!    return
!
!  end subroutine readVTK

end module VTKtools


!program loadVTK
!  use VTKtools
!  implicit none
!  ! NB: this is the dimension of the grid, which define the center of cell faces
!  integer, dimension(3) :: dimensions
!  real :: time
!  ! who created the file
!  character (:), allocatable :: origin
!  ! in spherical coordinates x1=r, x2=theta, x3=phi (coordinates of cell faces)
!  real, allocatable :: x1(:), x2(:), x3(:)
!  ! in spherical coordinates: vx1= vr, vx2=vtheta, vx3=vphi
!  ! note that the field are at cell centers (hence the size is reduced by 1 in each direction)
!  real, allocatable, dimension(:,:,:) :: rho, vx1, vx2, vx3
!
!  !! some variables to play around
!  integer nx2mid, i,j,k
!
!  character(len=13) :: filename = "data.0000.vtk"
!
!  call readVTK(filename,dimensions, time, origin, x1,x2,x3, rho, vx1,vx2,vx3)
!
!  print *,"File created by ",origin
!
!  print *,"Midplane density:"
!  print *,"-------------------------"
!
!  nx2mid = (dimensions(2)-1)/2
!  do i = 1, dimensions(1)-1
!    print *, "r=",x1(i), "rho=",rho(i,nx2mid,1)
!  enddo
!  print *,"-------------------------"
!
!end
