MODULE Planck

 use constant
 use spectrum_type, only : NLTEspec

 IMPLICIT NONE

 real(8), parameter :: MAX_EXPONENT = 150.0

 CONTAINS

 SUBROUTINE Bplanck (T, Bnu)
 ! -----------------------------------------------------
 ! Return an array of planck functions at all wavelengths
 ! for the cell temperature.
 ! Bnu in J/s/m2/Hz/sr
 ! lambda in nm
 ! -----------------------------------------------------
  double precision, intent(in) :: T
  double precision, dimension(NLTEspec%Nwaves), intent(out) :: Bnu
  double precision, dimension(NLTEspec%Nwaves) :: hnu_kT, twohnu3_c2

   twohnu3_c2 = 2.*HPLANCK*CLIGHT / (NM_TO_M * NLTEspec%lambda)**3

   hnu_kT = (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*NLTEspec%lambda*T)
   where(hnu_kT.lt.MAX_EXPONENT)
     Bnu = twohnu3_c2 / (dexp(hnu_kT)-1.)
   else where
     Bnu = 0. ! exponential is infinite, Bnu goes to zero
   end where

 RETURN
 END SUBROUTINE Bplanck

END MODULE Planck
