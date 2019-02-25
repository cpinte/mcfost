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
  use input
  use atmos_type, only : atmos, init_atomic_atmos
  use getline
  
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE read_hdF5(filename, x, y, z)
  ! ---------------------------------------------------------- !
   ! Read HDF5 files, formatted by C. Zanni and G. Pantolmos
   ! Presenlty.
  ! ---------------------------------------------------------- !
   character(len=*), intent(in)	:: filename
   integer :: icell, Nread, syst_status
   character(len=MAX_LENGTH) :: inputline, FormatLine
   double precision, dimension(:), allocatable, intent(out) :: x, y, z
  
  !CALL init_atomic_atmos()
  !CALL pluto_to_mcfost(Nmodel,x,y,z,h,vx,vy,vz)

  RETURN
  END SUBROUTINE read_hdf5


  SUBROUTINE pluto_to_mcfost(x,y,z)
   double precision, intent(inout), allocatable, dimension(:) :: x,y,z
  ! Convert phantom quantities & units to mcfost quantities & units
  ! x,y,z are in au
  ! densities assumed in kg/m3 here and velocities in m/s

  
  !Force using MCFOST stars
  lfix_star = .true.
  if (lfix_star) then
   write(*,*) "Using MCFOST stars"
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


 RETURN
 END SUBROUTINE pluto_to_mcfost

  SUBROUTINE readAtmos_ascii(filename, x, y, z)
  ! ------------------------------------------ !
   ! read atmos ascii file in GridType atmos.
   ! ultimately this model will be mapped
   ! onto a Voronoi mesh.
   ! Note: The temperature is not written !
  ! ------------------------------------------ !
   use getline
   character(len=*), intent(in)	:: filename
   integer :: icell, Nread, syst_status
   character(len=MAX_LENGTH) :: inputline, FormatLine, cmd
   double precision, dimension(:), allocatable, intent(out) :: x, y, z
   
   write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
   
   
   cmd = "wc -l "//trim(filename)//" > ntest.txt"
   call appel_syst(cmd,syst_status)
   open(unit=1,file="ntest.txt",status="old")
   read(1,*) n_cells
   close(unit=1)
   write(*,*) "Found ", n_cells, "voronoi points"
   
   CALL init_atomic_atmos() !init atmos structure for atomline RT
   allocate(x(n_cells), y(n_Cells), z(n_cells))
   
   open(unit=1,file=trim(filename))
   do icell=1, atmos%Nspace
    CALL getnextline(1, "#", FormatLine, inputline, Nread)
    READ(inputLine(1:Nread), *) x(icell), y(icell), z(icell), atmos%nHtot(icell), &
    	atmos%Vxyz(icell,1),atmos%Vxyz(icell,2), atmos%Vxyz(icell,3)

   end do
   
   close(unit=1)
   
   !also define stellar rings in this ?
   CALL pluto_to_mcfost(x,y,z)

  RETURN
  END SUBROUTINE readAtmos_ascii

END MODULE READ_PLUTO
