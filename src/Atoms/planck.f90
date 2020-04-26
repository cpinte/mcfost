module Planck

 use constant
 !!use spectrum_type, only : NLTEspec
 use mcfost_env, only : dp

 implicit none

 real(kind=dp), parameter :: MAX_EXPONENT = 600.0

 contains

!  subroutine Bplanck (T, Bnu)
!  ! -----------------------------------------------------
!  ! Return an array of planck functions at all wavelengths
!  ! for the cell temperature.
!  ! Bnu in J/s/m2/Hz/sr
!  ! lambda in nm
!  ! -----------------------------------------------------
!   real(kind=dp), intent(in) :: T
!   real(kind=dp), dimension(NLTEspec%Nwaves), intent(out) :: Bnu
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: hnu_kT, twohnu3_c2
! 
!    twohnu3_c2 = 2.*HPLANCK*CLIGHT / (NM_TO_M * NLTEspec%lambda)**3
! 
!    hnu_kT = (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*NLTEspec%lambda*T)
!    where(hnu_kT < MAX_EXPONENT)
!      Bnu = twohnu3_c2 / (exp(hnu_kT)-1.)
!    else where
!      Bnu = 0. ! exponential is infinite, Bnu goes to zero
!    end where
! 
!  return
!  end subroutine Bplanck
 
 elemental Function Bpnu (T, lambda)
 ! -----------------------------------------------------
 ! Return an array of planck functions at all wavelengths
 ! for the cell temperature.
 ! Bnu in J/s/m2/Hz/sr
 ! lambda in nm
 ! -----------------------------------------------------
  real(kind=dp), intent(in) :: T, lambda
  !real(kind=dp), intent(out) :: Bpnu
  real(kind=dp) :: hnu_kT, twohnu3_c2, Bpnu

   twohnu3_c2 = 2.*HPLANCK*CLIGHT / (NM_TO_M * lambda)**3
   hnu_kT = (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*lambda*T)
   
   if (hnu_kT < MAX_EXPONENT) then 
     Bpnu = twohnu3_c2 / (exp(hnu_kT)-1.)
   else
     Bpnu = 0. ! exponential is infinite, Bnu goes to zero
   end if

 return
 end function bpnu
 
!  subroutine dBnu_dT (T, dBnu)
!  ! -----------------------------------------------------
!  ! Return an array of derivatives of planck functions at all wavelengths
!  ! for temperature T
!  ! -----------------------------------------------------
!   real(kind=dp), intent(in) :: T
!   real(kind=dp), dimension(NLTEspec%Nwaves), intent(out) :: dBnu
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: hnu_kT, Bnu
! 
!    hnu_kT= (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*NLTEspec%lambda*T)
!    CALL Bplanck(T, Bnu)
!    
!    dBnu(:) = Bnu(:)*hnu_kT / T * exp(hnu_kT) / (exp(hnu_kT)-1)
! 
!  return
!  end subroutine dBnu_dT
 
!  subroutine TdBnu_dT_Bnu (T, dBnu)
!  ! -----------------------------------------------------
!  ! Return an array of derivatives of planck functions at all wavelengths
!  ! for temperature T divided by Bnu times T
!  ! -----------------------------------------------------
!   real(kind=dp), intent(in) :: T
!   real(kind=dp), dimension(NLTEspec%Nwaves), intent(out) :: dBnu
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: hnu_kT
! 
!    hnu_kT= (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*NLTEspec%lambda*T)
!    
!    dBnu(:) = hnu_kT * exp(hnu_kT) / (exp(hnu_kT)-1)
! 
!  return
!  end subroutine TdBnu_dT_Bnu
 
!  function uLD (T) result(ulimb)
!  ! -----------------------------------------------------
!  ! -----------------------------------------------------
!   real(kind=dp), intent(in) :: T
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: ulimb
! 
!    CALL TdBnu_dT_Bnu(T, ulimb)
! 
!  return
!  end function uLD

end module Planck
