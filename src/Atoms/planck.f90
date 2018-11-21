MODULE Planck

 use constant

 IMPLICIT NONE

 real(8), parameter :: MAX_EXPONENT = 150.0

 CONTAINS

 SUBROUTINE Bplanck (T, Nlam, lambda, Bnu)
 ! -----------------------------------------------------
 ! Return an array of planck functions at all wavelengths
 ! for the cell temperature.
 ! Bnu in J/s/m2/Hz/sr
 ! lambda in nm
 ! -----------------------------------------------------
  integer, intent(in) :: Nlam
  double precision, intent(in) :: T
  double precision, dimension(Nlam), intent(out) :: Bnu
  double precision, intent(in), dimension(Nlam) :: lambda
  double precision, dimension(Nlam) :: hnu_kT, twohnu3_c2

   twohnu3_c2 = 2.*HPLANCK*CLIGHT / (NM_TO_M * lambda)**3

   hnu_kT = (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*lambda*T)
   where(hnu_kT.lt.MAX_EXPONENT)
     Bnu = twohnu3_c2 / (dexp(hnu_kT)-1.)
   else where
     Bnu = 0. ! exponential is infinite, Bnu goes to zero
   end where

 RETURN
 END SUBROUTINE Bplanck

END MODULE Planck
