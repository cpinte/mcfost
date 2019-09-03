MODULE thomson_scattering

 use atmos_type, only : atmos
 use constant
 use mcfost_env, only : dp

 IMPLICIT NONE

 CONTAINS

 FUNCTION Thomson(icell)
 ! ------------------------------------------------------------- !
  ! Thomson scattering cross-section in non relativistic limit.
  ! (i.e., wavelength independent)
  ! unit m^2
 ! ------------------------------------------------------------- !
  integer 			:: icell
  real(kind=dp)     :: Thomson


  Thomson = atmos%ne(icell) * sigma_e

 RETURN
 END FUNCTION Thomson


END MODULE thomson_scattering
