! --------------------------------------------------------------------------- !
 ! Implements Voigt function procedures
! --------------------------------------------------------------------------- !
MODULE voigtfunctions

 use constant
 use mcfost_env, only : dp
 use math, only : locate

 IMPLICIT NONE
 
 PROCEDURE(VoigtHumlicek), pointer :: Voigt => null()

 CONTAINS
 
  FUNCTION VoigtHumlicek(N, a, v, F) result(L)
  !--------------------------------------------------------------
  ! Generates Voigt and anomalous dispersion profiles
  ! See Humlicek (1982) JQSRT 27, 437
  ! W4 = Voigt + i*Faraday-Voigt (dispersion profile)
  !--------------------------------------------------------------
	real(kind=dp) :: a
	real(kind=dp) :: v(N), L(N)
	real(kind=dp), intent(out), optional :: F(N)
	complex(kind=dp) :: w4, z, u, v4!,t
	real(kind=dp) :: s
	integer :: i, N
		
    do i = 1, n
		!!!z = cmplx(v(i), a)
		!t = cmplx(a, -v(i))
		z = cmplx(a,-v(i))
		s = dabs(v(i)) + a
		u = z*z

		if (s >= 15.d0) then
		! Region I
			w4 = z * 0.5641896d0 / (0.5d0+u)
		elseif (s >= 5.5) then
		! Region II
			w4 = z*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
		elseif (a >= 0.195d0*dabs(v(i))-0.176d0) then
		! Region III
			w4 = (16.4955d0+z*(20.20933d0+z*(11.96482d0+z*(3.778987d0+z*0.5642236d0)))) / &
					(16.4955d0+z*(38.82363d0+z*(39.27121d0+z*(21.69274d0+z*(6.699398d0+z)))))
		else 
		! Region IV
			w4 = z*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
					u*(1.320522d0-u*0.56419d0))))))
			v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
					u*(61.57037d0-u*(1.841439d0-u)))))))
			w4 = cdexp(u) - w4/v4
		endif
			L(i) = dble(w4)
			if (present(F)) F(i) = aimag(w4)
	enddo  

  RETURN
  END FUNCTION VoigtHumlicek
  
  FUNCTION dirac_line(N, a, v, F) result(L)
	real(kind=dp) :: a
	real(kind=dp) :: v(N), L(N)
	real(kind=dp), intent(out), optional :: F(N)
	integer :: i, N
	
	L(:) = 0d0
	i = locate(v, 0d0)
	L(i) = 1d0

  RETURN
  END FUNCTION dirac_line
  
  FUNCTION VoigtXXXX(N, a, v, F) result(L) !for unpolarised profile
	real(kind=dp) :: a
	real(kind=dp) :: v(N), L(N)
	real(kind=dp), intent(out), optional :: F(N)
	integer :: i, N
  !TBD
  RETURN
  END FUNCTION VoigtXXXX

END MODULE voigtfunctions


!  SUBROUTINE Humlicek(N, a, v, W)
!   ! W = L + i*F
!   integer, intent(in)              :: N
!   double precision, intent(in)     :: a
!   double precision, intent(in)     :: v(N)
!   complex(kind=16), intent(out)    :: W(N)
!   double precision                 :: s(N)
!   complex(kind=16)                 :: z(N)
!   complex(kind=16)                 :: u(N)
! 
!   z = cmplx(a,-v)
!   s = abs(v) + a
! 
!   ! Approximation in region IV
!   u = z * z
!   W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
!    u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
! (32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
!    u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))
! 
!   where(s >= 15)
!    ! Approximation in region I
!    W = (z * 0.5641896) / (0.5 + (z * z))
!   else where (s >= 5.5)
!    ! Approximation in region II
!    u = z * z
!    W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
!   else where (a >= 0.195*abs(v) - 0.176)
!   ! Approximation in region III
!    W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
!    0.5642236*z)))) / &
!    (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
!    z*(6.699398 + z)))))
!   end where
! 
!  RETURN
!  END SUBROUTINE Humlicek

!  FUNCTION VoigtHumlicek(N, a,v, F) result(L)
!  ! --------------------------------------------------------------------------- !
!   ! Humlicek 1982, JQSRT 27, p. 437
!   ! Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
!   ! function needed in Stokes radiative transfer.
!   ! W = Voigt + i*Faraday-Voigt (dispersion profile)
!   !
!   ! Adapted from rh code (Uitenbroek 2001, ApJ 389-398)
!  ! --------------------------------------------------------------------------- !
!   integer, intent(in)           :: N
!   double precision, intent(in)  :: a
!   double precision, intent(in)  :: v(N)
!   double precision              :: L(N)
!   double precision, intent(out) :: F(N)
!   complex(kind=16)              :: W(N)
! 
!   ! real of W is L, and imag is F
!   CALL Humlicek(N, a, v, W)
! 
!   F = imag(W)
!   L = real(W) !factor 2 here or not ?
!  RETURN
!  END FUNCTION VoigtHumlicek

  
!  double precision, dimension(0:23) :: an
!  data an /0.2954089751509193d0, 0.27584023329217705d0, 0.22457395522461585d0,0.1594149382739117d0,&
!     0.09866576641545417d0, 0.05324414078763942d0, 0.025052150005393646d0, 0.010277465670539542d0, &
!     0.003676164332844841d0, 0.0011464936412422387d0, 0.0003117570150461923d0, 0.00007391433429603544d0, &
!     0.00001527949342800332d0, 2.7539566082259447d-6, 4.327858781854298d-7, 5.9300304091976984d-8 , &
!     7.084490303405148d-9 , 7.379520677523531d-10 , 6.702171205812766d-11 , 5.30726900090507d-12 , &
!     3.6642873348622373d-13 , 2.2062473036786997d-14 , 1.1544519576969374d-15 ,5.62191094758381d-17 /
!  
!  Function an_abra(j, tM) result(y)
!   integer, intent(in) :: j, tM
!   double precision :: y
!   
!    y = 2d0 * dsqrt(PI) / real(tM) * dexp(-j*j * PI*PI / real(tM) / real(tM)) 
!  return
!  end function an_abra
!  
!  FUNCTION VoigtAbrarov(N, a, v, F) result(L)
!  ! --------------------------------------------------------------------------- !
!   ! S.M. Abrarov et al. JCPC (2010) 876–882
!   ! Rapidly convergent series for high-accuracy calculation of the
!   ! Voigt function
!   !
!   ! Accuracy up to 1d-13 ?
!  ! --------------------------------------------------------------------------- !
!   double precision              :: a
!   double precision				:: v(N)
!   integer						:: N
!   double precision              :: L(N)
!   !double precision              :: x(N)
!   !double precision              :: y
!   double precision, intent(out) :: F(N)
!   integer                       :: j
!   integer, parameter 			:: M = 24
!   integer, parameter            :: tM = 12
!   !integer						:: tM2
!   !double precision              :: xpy(N)
!   complex(kind=16)				:: jcmplx, z(N), w(N), A1(N), B(N), C(N) !,K(2, N)
!   
!   jcmplx = cmplx(0d0, 1d0)
! 
!   !Slower then humliceck why ?
! ! --> Fatser, less accurate
!   z = cmplx(v, a)    																		  
!   A1 = tM*z; B = exp(jcmplx*A1); C = A1*A1
!   w = jcmplx * (1.-B) / A1
!   do j=1,M-1!M with an_abra(j, tM)
!    w = w + jcmplx * A1 * tM / dsqrt(PI) * an(j) * ((-1)**j * B - 1) / (j*j * PI*PI - C)
!   end do
!   L = real(w); F = imag(w)
!   
! ! --> to debug, it is more accurate
! !   K(:,:) = 0d0  
! !   x = dsqrt(log(2d0)) * v
! !   y = dsqrt(log(2d0)) * a
! !   xpy = (x**2 + y**2)
! !   an(:) = an(:)/(2d0*dsqrt(PI))
! ! 
! !   tM2 = tM*tM
! !   do j=0, M-1
! !   
! !    K(1,:) = K(1,:) + an(j) * (((jcmplx*j*PI*tM + tM2*Y)*(1-exp(-jcmplx*j*PI-tM*Y)*cos(tM*x))+exp(-jcmplx*j*PI-tM*y)*tM2*x*sin(tM*x))/&
! !    (tM2*x**2 - (j*PI-jcmplx*tM*y)**2) - &
! !    ((jcmplx*j*PI*tM-tM2*y)*(1-exp(jcmplx*j*PI-tM*y)*cos(tM*x))-exp(jcmplx*j*PI-tM*y)*tM2*x*sin(tM*x)) /&
! !    (tM2*x**2 - (j*PI+jcmplx*tM*y)**2))
! !   
! !    K(2,:) = K(2,:) + an(j) * ((tM2*x - exp(jcmplx * j * PI + tM*y)*(tM2*x*cos(tM*x)+(jcmplx*j*PI*tM+tM2*y)*sin(tM*x)))/&
! !    (tM2*x**2 - (j*PI - jcmplx*tM*y)**2)+&
! !        (tM2*x - exp(jcmplx*j*PI -tM*y)*(tM2*x*cos(tM*x)-(jcmplx*j*PI*tM - tM2*y)*sin(tM*x)))/&
! !    (tM2*x**2 - (j*PI + jcmplx*tM*y)**2))
! ! 
! !   end do
! !   L = real(K(1,:)) - an(0) * (y-exp(-tM*y)*(y*cos(tM*x)-x*sin(tM*x))) / xpy
! !   F = imag(K(2,:)) - an(0)*(x-exp(-tM*y)*(x*cos(tM*x)+y*sin(tM*x))) / xpy
!   
! 
!  RETURN
!  END FUNCTION VoigtAbrarov

!! Deprecated I use PROCEDURE NOW
!  FUNCTION Voigt(N, a, v, F, VoigtAlgorithm) result (L)
!  ! --------------------------------------------------------------------------- !
!   ! RETURN THE VOIGT FUNCTION FOR A LINE PROFILE
!   ! ACCORDING TO A DETERMINED ALGORITHM
!   ! L-> Voigt function, F->Dispersion profile
!   ! if desired (when polarisation is desired in fact)
!   ! VoigtAlgorithm-> aglorithm for the voigt function
!  ! --------------------------------------------------------------------------- !
!   integer, intent(in)                       :: N !size of the final voigt
!   double precision, intent(in)              :: a
!   double precision, intent(in)				:: v(N)
!   double precision, intent(out)				:: F(N)
!   double precision			                :: L(N)
!   character(len=*), optional    :: VoigtAlgorithm
! 
!   if (.not. present(VoigtAlgorithm)) VoigtAlgorithm="HUMLICEK"
! 
!   SELECT CASE (VoigtAlgorithm)
!    CASE ("HUMLICEK")
!     L = VoigtHumlicek(N, a, v, F)
!    CASE ("ABRAROV")
!     L = VoigtAbrarov(N, a, v, F)
!    CASE DEFAULT
!     L = VoigtHumlicek(N, a, v, F) !using Humlicek by default
!  END SELECT
! 
!  RETURN
!  END FUNCTION Voigt