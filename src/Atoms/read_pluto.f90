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
  stop
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
  
   x = x * m_to_AU
   y = y * m_to_AU
   z = z * m_to_AU
   h = h * m_to_AU


 RETURN
 END SUBROUTINE pluto_to_mcfost


END MODULE READ_PLUTO
