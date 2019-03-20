! --------------------------------------------------------------------------- !
 ! Implements Voigt function generators
! --------------------------------------------------------------------------- !
MODULE voigtfunctions

 use constant

 IMPLICIT NONE

 CONTAINS

 SUBROUTINE Humlicek(N, a, v, W)
  ! W = L + i*F
  integer, intent(in)              :: N
  double precision, intent(in)     :: a
  double precision, intent(in)     :: v(N)
  complex(kind=16), intent(out)    :: W(N)
  double precision                 :: s(N)
  complex(kind=16)                 :: z(N)
  complex(kind=16)                 :: u(N)

  z = cmplx(a,-v)
  s = abs(v) + a

  ! Approximation in region IV
  u = z * z
  W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
   u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
(32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
   u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))

  where(s >= 15)
   ! Approximation in region I
   W = (z * 0.5641896) / (0.5 + (z * z))
  else where (s >= 5.5)
   ! Approximation in region II
   u = z * z
   W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
  else where (a >= 0.195*abs(v) - 0.176)
  ! Approximation in region III
   W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
   0.5642236*z)))) / &
   (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
   z*(6.699398 + z)))))
  end where

 RETURN
 END SUBROUTINE Humlicek

 FUNCTION VoigtHumlicek(N, a,v, F) result(L)
 ! --------------------------------------------------------------------------- !
  ! Humlicek 1982, JQSRT 27, p. 437
  ! Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
  ! function needed in Stokes radiative transfer.
  ! W = Voigt + i*Faraday-Voigt (dispersion profile)
  !
  ! Adapted from rh code (Uitenbroek 2001, ApJ 389-398)
 ! --------------------------------------------------------------------------- !
  integer, intent(in)           :: N
  double precision, intent(in)  :: a
  double precision, intent(in)  :: v(N)
  double precision              :: L(N)
  double precision, intent(out) :: F(N)
  complex(kind=16)              :: W(N)

  ! real of W is L, and imag is F
  CALL Humlicek(N, a, v, W)

  F = imag(W)
  L = real(W) !factor 2 here or not ?
 RETURN
 END FUNCTION VoigtHumlicek
 
 FUNCTION VoigtAbrarov(N, a, v, F) result(L)
 ! --------------------------------------------------------------------------- !
  ! S.M. Abrarov et al. JQSRT (2010) 372–375
  ! Rapidly convergent series for high-accuracy calculation of the
  ! Voigt function
  !
  ! Accuracy up to 1d-13 ?
  ! TBD: Faraday-Voigt profile
 ! --------------------------------------------------------------------------- !
  double precision              :: a
  double precision				:: v(N)
  integer                       :: N
  double precision              :: L(N)
  double precision              :: x(N)
  double precision              :: y
  double precision              :: K
  double precision, intent(out) :: F(N)
  integer                       :: j
  integer                       :: tM
  double precision              :: cosx(N)
  double precision              :: sinx(N)
  double precision              :: xpy(N)
  double precision              :: ymx(N)
  double precision              :: exptMy
  double precision              :: exptMyp
  double precision              :: tMx(N)
  double precision              :: tMy
  
  x = dsqrt(log(2d0)) * v
  y = dsqrt(log(2d0)) * a
  tM = 12
  
  xpy = (x**2 + y**2)
  ymx = (y**2 - x**2)
  cosx = cos(tM*x)
  sinx = sin(tM*x)
  tMy = tM * y
  exptMy = dexp(-tMy)
  exptMyp = dexp(tMy)
  tMx = tM * x
  
  L = - exptMy * (exptMyp*y - y*cosx+x*sinx) / (tM * xpy)
  do j=0,2*tM - 1
   K = 2d0 / tM * dexp(-(j*PI/tM)**2)
   L = L + K * tM**2 * y * (j**2 * PI**2 + tM**2 * xpy) & 
      / (j**4 * PI**4 + 2*j**2 * PI**2 * tM**2 * ymx + tM**4 *xpy**2) + &
      0.5*tM * exptMy * (&
       ((j*PI-tMx)*sin(j*PI-tMx) - tMy*cos(j*PI-tMx))/&
        (j**2 * PI**2 - 2*j*PI*tMx + tM**2*xpy) + &
       ((j*PI + tMx)*sin(j*PI+tMx)-tMy*cos(j*PI+tMx)) /& 
        (j**2 * PI**2 + 2*j*PI*tMx + tM**2*xpy) &
      )*K
  end do
  
  F = 0d0

 RETURN
 END FUNCTION VoigtAbrarov

 FUNCTION Voigt(N, a, v, F, VoigtAlgorithm) result (L)
 ! --------------------------------------------------------------------------- !
  ! RETURN THE VOIGT FUNCTION FOR A LINE PROFILE
  ! ACCORDING TO A DETERMINED ALGORITHM
  ! L-> Voigt function, F->Dispersion profile
  ! if desired (when polarisation is desired in fact)
  ! VoigtAlgorithm-> aglorithm for the voigt function
 ! --------------------------------------------------------------------------- !
  integer, intent(in)                       :: N !size of the final voigt
  double precision, intent(in)              :: a
  double precision, intent(in)				:: v(N)
  double precision, intent(out)				:: F(N)
  double precision			                :: L(N)
  character(len=*), optional    :: VoigtAlgorithm

  if (.not. present(VoigtAlgorithm)) VoigtAlgorithm="HUMLICEK"

  SELECT CASE (VoigtAlgorithm)
   CASE ("HUMLICEK")
    L = VoigtHumlicek(N, a, v, F)
   CASE ("ABRAROV2010")
    L = VoigtAbrarov(N, a, v, F)
   CASE DEFAULT
    L = VoigtHumlicek(N, a, v, F) !using Humlicek by default
 END SELECT

 RETURN
 END FUNCTION Voigt


END MODULE voigtfunctions
