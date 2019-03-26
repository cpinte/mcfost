!
! First attempt to read pluto hdf5 files into MCFOST
!
!
!
!
!
!
module READ_PLUTO

  use parametres
  use messages
  use constantes
  use getline
  
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE read_hdF5(filename, N_points, x, y, z, h, rho,vx, vy, vz)
  ! ---------------------------------------------------------- !
   ! Read HDF5 files, formatted by C. Zanni and G. Pantolmos
   ! Presenlty.
  ! ---------------------------------------------------------- !
   character(len=*), intent(in)	:: filename
   integer :: icell, Nread, syst_status
   character(len=MAX_LENGTH) :: inputline, FormatLine, cmd
   integer, intent(inout)										:: N_points
   double precision, dimension(:), allocatable, intent(out) :: x, y, z, h
   double precision, dimension(:), allocatable, intent(out) :: rho, vx, vy, vz
  
   !density_files(1) assumes only one model for pluto
   !read here ....
   ! h = 3d0 * Rsun !assuming a large h to now
  
   !convert from pluto to mcfost units
   CALL pluto_to_mcfost(N_points,x,y,z,h)!,vx,vy,vz)

  RETURN
  END SUBROUTINE read_hdf5


  SUBROUTINE pluto_to_mcfost(N, x,y,z,h)
   integer, intent(in) 										  :: N
   double precision, intent(inout), allocatable, dimension(:) :: x,y,z,h

  ! Convert PLUTO quantities & units to mcfost quantities & units
  ! x,y,z are in au
  ! densities assumed in kg/m3 here and velocities in m/s
  ! mcfost densities are in g/cm3 for dust and mol RT. But atom RT
  ! uses SI.
  ! Defines rings here ? and take into account them later?
  ! before tesselation?

  
  !Force using MCFOST stars
  lfix_star = .true.
  if (lfix_star) then
   write(*,*) "Using ", n_etoiles, "MCFOST stars"
   !adding some rings/spots at the surface
  else
   write(*,*) "Stars other than MCFOST defined are not allowed, exiting"
   stop
  end if

  !density in PLUTO units to kg/m3
  !velocites --> m/s
  !positions read in m converted in AU
  
  if (lascii_pluto_file) then
  ! grid written to ascii is in meters
   x = x * m_to_AU
   y = y * m_to_AU
   z = z * m_to_AU
   h = h * m_to_AU
  else !HDF5 true pluto's model
   !...
  end if

 RETURN
 END SUBROUTINE pluto_to_mcfost

  SUBROUTINE readAtmos_ascii(filename, N_points, x, y, z, h, rho,vx, vy, vz)
  ! ------------------------------------------ !
   ! read atmos ascii file in GridType atmos.
   ! ultimately this model will be mapped
   ! onto a Voronoi mesh.
   ! Note: The temperature is not written !
  ! ------------------------------------------ !
   use getline
   use utils, only : appel_syst
   character(len=*), intent(in)	:: filename
   integer :: icell, Nread, syst_status
   character(len=MAX_LENGTH) :: inputline, FormatLine, cmd
   integer, intent(out)										:: N_points
   double precision, dimension(:), allocatable, intent(out) :: x, y, z, h
   double precision, dimension(:), allocatable, intent(out) :: rho, vx, vy, vz
   
   write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
   
   
   cmd = "wc -l "//trim(filename)//" > ntest.txt"
   call appel_syst(cmd,syst_status)
   open(unit=1,file="ntest.txt",status="old")
   read(1,*) N_points
   close(unit=1)
   write(*,*) "Found ", N_points, "voronoi points"
   
   allocate(x(N_points), y(N_points), z(N_points), h(N_points), rho(N_points))
   allocate(vx(N_points), vy(N_points), vz(N_points))

   
   open(unit=1,file=trim(filename))
   do icell=1, N_points
    CALL getnextline(1, "#", FormatLine, inputline, Nread)
    READ(inputLine(1:Nread), *) x(icell), y(icell), z(icell), rho(icell), &
    	vx(icell), vy(icell), vz(icell), h(icell)

   end do
   
   close(unit=1)
   
   !also define stellar rings in this ?
   CALL pluto_to_mcfost(N_points, x,y,z,h)

  RETURN
  END SUBROUTINE readAtmos_ascii

END MODULE READ_PLUTO
