MODULE thomson_scattering

 use grid_type, only: atmos
 use constant
 use math, only : dpow

 IMPLICIT NONE

 CONTAINS

 FUNCTION Thomson(icell) result (scatt)
 ! Thomson scattering cross-section in non relativistic limit.
 ! (i.e., wavelength independent)
 ! unit m^2
  integer :: icell
  double precision :: sigma, scatt

  sigma = 8.*PI/3. * dpow(Q_ELECTRON/(sqrt(4.0*PI*EPSILON_0)*&
          (sqrt(M_ELECTRON)*CLIGHT)), 4d0)

  scatt = atmos%ne(icell) * sigma

 RETURN
 END FUNCTION Thomson


END MODULE thomson_scattering
