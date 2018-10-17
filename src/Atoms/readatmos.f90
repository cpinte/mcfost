!
!
! This module reads external model "atmosphere" file from MHD simulations.
! Standard fits format is used.
!
! atmos structure is allocated here for the RT Calculations
!
! Note: The word "atmosphere" is a general designation for accretion disk models,
! or whatever. It only requires the total Hydrogen densities (nHtot), eventually 
! electron densities , Temperature, velocity fields, magnetic fields etc. Further
! the conditions of atomic line transfer should apply.
!
MODULE READATMOS

 use grid_type, only : atmos
 use project
 use solvene, only : SolveElectronDensity
 use constant

 IMPLICIT NONE 
 
 CONTAINS
 
 SUBROUTINE fillAtmos()
 
 
 END SUBROUTINE fillAtmos

END MODULE READATMOS