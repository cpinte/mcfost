! This module contains functions to project the velocity
! and the magnetic field according the parameters of the cell
! in which a (local) line profile is computed.
! Use MCFOST routines.


MODULE project

 use grid_type, only : atmos, cell_point

 IMPLICIT NONE

 CONTAINS

! Taken from MCFOST
 FUNCTION vproject(icell) result(y)
  ! velocity is taken from the model
  integer :: icell
  double precision :: y

  y = atmos%ux(icell) * cell_point%u + atmos%uy(icell) * cell_point%v + &
      atmos%uz(icell) * cell_point%w

 RETURN
 END FUNCTION

 FUNCTION Bproject(icell) result(y)
  integer :: icell
  double precision :: y
  
  y = 0.

 END FUNCTION

END MODULE project
