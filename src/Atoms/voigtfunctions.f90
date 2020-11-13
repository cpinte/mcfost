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
	complex(kind=dp) :: w4, z, u, v4
	real(kind=dp) :: s
	integer :: i, N
		
    do i = 1, n

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
			w4 = exp(u) - w4/v4 !cdexp, can be optimized by splitting the exp on real and imag parts
		endif
		L(i) = real(w4,kind=dp)!w4%re!!dble(w4)
		if (present(F)) F(i) = aimag(w4)!w4%im!
		
	enddo  

  RETURN
  END FUNCTION VoigtHumlicek
  
  FUNCTION dirac_line(N, v) result(L)
	real(kind=dp) :: v(N), L(N), minv = 1d20
	integer :: i, N, i1

! 	L(:) = 0.0_dp
! 	i = locate(v, 0.0_dp)
! 	L(i) = 1.0_dp
	
	i1 = N
	do i=1, N
		minv = min(minv, abs(v(i)))
		if (abs(v(i)) == minv) then
			v(i) = 1.0
			v(i1) = 0.0
			i1 = i
		else
			v(i) = 0.0
		endif
	
	enddo


  RETURN
  END FUNCTION dirac_line
  
  
   FUNCTION gate_line(N, v, vlim) result(L)
	real(kind=dp) :: vlim
	real(kind=dp) :: v(N), L(N)
	integer :: i, N
	
	do i=1, N
	
		if (abs(v(i)) > vlim ) then
			L(i) = 0.0_dp
		else
			L(i) = 1.0_dp
		endif
	
	enddo

  RETURN
  END FUNCTION gate_line 
  
  FUNCTION Voigt_thomson(N, a, vbroad, v) result(L) !for unpolarised profile
	real(kind=dp) :: a, vbroad
	real(kind=dp) :: v(N), L(N)
	integer :: i, N
	real(kind=dp) :: aeff, ratio, eta, al

	al = a * vbroad
		
	aeff = (vbroad**5. + 2.69269*vbroad**4. * aL + 2.42843*vbroad**3. * aL**2. + &
					4.47163*vbroad**2. *aL**3. + 0.07842*vbroad*aL**4. + aL**5.)**(0.2)
          		

          	
	ratio = aL/aeff
	eta = 1.36603*ratio - 0.47719*(ratio*ratio) + 0.11116*(ratio*ratio*ratio)
			
         
	L(:) = eta * ( aeff/pi * (v**2 + aeff**2)**(-1.0) ) + &
			(1.0_dp - eta) * exp(-(v/aeff)**2) / aeff / sqrtpi
					

  RETURN
  END FUNCTION Voigt_thomson

END MODULE voigtfunctions


  
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
