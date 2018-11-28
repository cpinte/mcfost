! --------------------------------------------------------------------------- !
 ! Implements Voigt functions generator for atomic line radiative transfer
! --------------------------------------------------------------------------- !
MODULE voigtfunctions

 use constant

 IMPLICIT NONE

!  integer, parameter :: NC=34, NT=10, NGR=31
!  double precision, parameter :: EXPMAX=70.,TWOSQRTPI=1.12837917,&
!      C0=0.0897935610625833,C1=29.608813203268, THREEPI=9.42477796076938,&
!      INVSQRTPI=0.564189583547756
!  double precision, parameter :: TINY=1d-8
!  double precision, dimension(NT) :: t, w
! 
!  data t / 0.2453407083, 0.7374737285, 1.2340762153, &
!          1.7385377121, 2.2549740020, 2.7888060584, &
!          3.3478545673, 3.9447640401, 4.6036824495, &
!          5.3874808900 /
! 
!  data w / 4.6224366960d-01, 2.8667550536d-01, &
!          1.0901720602d-01, 2.4810520887d-02, &
!          3.2437733422d-03, 2.2833863601d-04, &
!          7.8025564785d-06, 1.0860693707d-07, &
!          4.3993409922d-10, 2.2293936455d-13 /


 CONTAINS


!  FUNCTION VoigtK1(N, a,v) result(L)
!   double precision :: a, v(:), c(NC)
!   integer :: N
!   double precision :: L(N)
!   integer :: nn
!   double precision :: a2, v2(N), u1(N), dn01(N), dn02(N), dn(N)
!   double precision :: v2i(N), funct(N), an, q, g(N), coef(N), bn01(N), bn02(N), bn(N), v1(N)
! 
! 
! 
!   DATA c / 0.1999999999972224, -0.1840000000029998,            &
!            0.1558399999965025, -0.1216640000043988,            &
!            0.0877081599940391, -0.0585141248086907,            &
!            0.0362157301623914, -0.0208497654398036,            &
!            0.0111960116346270, -0.56231896167109d-02,          &
!            0.26487634172265d-02, -0.11732670757704d-02,        &
!            0.4899519978088d-03, -0.1933630801528d-03,          &
!            0.722877446788d-04, -0.256555124979d-04,            &
!            0.86620736841d-05, -0.27876379719d-05,              &
!            0.8566873627d-06, -0.2518433784d-06,                &
!            0.709360221d-07, -0.191732257d-07, 0.49801256d-08,  &
!            -0.12447734d-08, 0.2997777d-09, -0.696450d-10,      &
!            0.156262d-10, -0.33897d-11, 0.7116d-12, -0.1447d-12,&
!            0.285d-13, -0.55d-14, 0.10d-14, -0.2d-15 /
! 
!   a2 = a*a
!   v2 = v*v
! 
!   u1 = dexp(a2-v2) * cos(2.*v*a)
!   where ((v2-a2).gt.EXPMAX)
!    u1 = 0.
!   end where
! 
!    bn01 = 0.0
!    bn02 = 0.0
!    v1   = v / 5.0
!    coef = 4.0 * v1*v1 - 2.0
!    do nn = NC, 1, -1 ! for(n=NC-1,n>=0,n--)
!     bn = coef*bn01 - bn02 + c(nn)
!     bn02 = bn01
!     bn01 = bn
!    end do
!    dn02 = v1*(bn - bn02)
!    dn01 = 1.0 - 2.0*v*dn02
! 
!   where (v.gt.5.)
!    v2i = 1./v2
!    dn01 = -v2i * (0.5 + v2i*(0.75 + v2i*(1.875 + &
!               v2i*(6.5625 + v2i*(29.53125 + &
!               v2i*(1162.4218 + v2i*1055.7421))))))
!    dn02 = (1.0 - dn01) / (2.0 * v)
!   end where
! 
! 
!   funct = a*dn01
!   if (a.gt.TINY) then
!    q = 1.
!    an = a
!    do nn=3,50
!     dn = (v*dn01 + dn02) * (-2./nn)
!     dn02 = dn01
!     dn01 = dn
!     if (mod(nn, 1).gt.0) then ! if (n%2)
!      q = -q
!      an = an*a2
!      g = dn*an
!      funct = funct + q*g
! !      if (dabs(g/funct).le.TINY) then
! !       L = (u1 - TWOSQRTPI*funct)
! !       RETURN
! !      end if
!      if (MAXVAL(dabs(g/funct)) <= TINY) then
!       L = (u1 - TWOSQRTPI*funct)
!       RETURN
!      end if
!     end if
!    end do
!   end if
!   L = (u1 - TWOSQRTPI*funct)
! 
!  RETURN
!  END FUNCTION VoigtK1
! 
!  FUNCTION VoigtK2(N, a,v) result(L)
!   double precision :: a, v(:)
!   integer :: N
!   double precision :: L(N)
!   double precision :: g(N), r(N), s(N), a2
!   integer :: nn
!   a2 = a*a
!   g = 0.
!   do nn=1,NT
!    r = t(nn) - v
!    s = t(nn) + v
!    ! beware here, by atan I mean atan(x)
!    ! in the interval [-pi/2,+pi/2] radians.
!    ! and log is natural (base e)
!    g = g+(4.0*t(nn)*t(nn) - 2.0) * (r*atan(r/a) + s*atan(s/a) -&
!           0.5*a*(log(a2 + r*r) + log(a2 + s*s))) * w(nn)
!   end do
!   L = g/PI
!  RETURN
!  END FUNCTION VoigtK2
! 
!  FUNCTION VoigtK3(N,a,v) result(L)
!   double precision :: a, v(:)
!   integer N
!   double precision :: L(N)
!   integer :: nn
!   double precision :: g(N), a2
!   g = 0.
!   a2 = a*a
! 
!   do nn=1,NT
!    g = g+(1.0/((v - t(nn))**2 + a2) + 1.0/((v + t(nn))**2 + a2)) * w(nn)
!   end do
!   L = (a*g) / PI
!  RETURN
!  END FUNCTION VoigtK3
! 
!  FUNCTION VoigtArmstrong(N, a, v) result(L)
!  ! --------------------------------------------------------------------------- !
!   ! Armstrong 1967, JQSRT 7, pp. 61-88
!   ! (slow for damping parameters larger than 1.5, accurate to
!   ! 6 significant figures).
!   !
!   ! Written by Han. Utienbroek from rh code (Uitenbroek 2001, ApJ 389-398)
!  ! --------------------------------------------------------------------------- !
!   integer :: N
!   double precision :: L(N)
!   double precision :: a, v(:)
!   logical :: times_minus1(N)
! 
!   times_minus1 = .false.
!   !!!!!!!!!!!!!!!!!!!
!   where (v < 0.0)
!      v=-v
!      times_minus1=.true. !trick here
!   end where
!   ! A CHANGER CAR ICI LA VALEUR DE V est transformé si v < 0
!   ! peut casser la grille d'entrée si v est une variable d'un objet
!   ! J'utilise un petit trick avec times_minus1, mais je comprends
!   ! pas trop ce qu'il se passe. Surement une connerie de Fortran...
!   !!!!!!!!!!!!!!!!!!!
! 
!   L = VoigtK3(N, a,v)
! 
!   where (((a.lt.1.).and.(v.lt.4.0)).or.((a.lt.1.8/(v+1.))))
!    L = VoigtK1(N, a,v)
!   else where ((a.lt.2.5).and.(v.lt.4.))
!    L = VoigtK2(N, a,v)
!   end where
! 
!   where (times_minus1)
!    v = -v
!    times_minus1 = .false. !trick here
!   end where
! 
!  RETURN
!  END FUNCTION VoigtArmstrong


 SUBROUTINE Humlicek(N, a, v, W)
  ! W = L + i*F
  double precision, intent(in) :: a, v(:)
  integer, intent(in) :: N
  complex(kind=8), intent(out) :: W(N)
  double precision :: s(N)
  complex(kind=8) :: z(N), u(N)

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
  ! Written by Han. Utienbroek from rh code (Uitenbroek 2001, ApJ 389-398)
 ! --------------------------------------------------------------------------- !
  double precision :: a, v(:)
  integer :: N
  double precision :: L(N)
  double precision, intent(out) :: F(N)
  complex(kind=8) :: W(N)

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
  ! TBD: Faraday-Voigt profile
 ! --------------------------------------------------------------------------- !
  double precision              :: a
  double precision              :: v(:)
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
 
 FUNCTION VoigtLookup (N, a, v, F) result(L)
  !beta test
  use spectrum_type, only : NLTEspec
  use math, only : interp2Darr
 
  integer :: N
  double precision :: a, v(:)
  double precision :: L(N)
  double precision, intent(out) :: F(N)
  double precision :: y(N,1), xa(1)

  
  if (.not.NLTEspec%precomputed_voigt) then
   write(*,*) "Error: Trying to use precomputed voigt table that does not exist"
   stop
  end if
  
  xa(1) = a
  y = interp2Darr(size(NLTEspec%v_values), NLTEspec%v_values, &
          size(NLTEspec%a_values), NLTEspec%a_values, NLTEspec%VoigtTable, &
          N, v, 1, xa)
  L = y(:,1)
  F = 0
  
 RETURN
 END FUNCTION VoigtLookup

 FUNCTION Voigt(N, a, v, F, VoigtAlgorithm) result (L)
 ! --------------------------------------------------------------------------- !
  ! RETURN THE VOIGT FUNCTION FOR A LINE PROFILE
  ! ACCORDING TO A DETERMINED ALGORITHM
  ! L-> Voigt function, F->Dispersion profile
  ! if desired (when polarisation is desired in fact)
  ! VoigtAlgorithm-> aglorithm for the voigt function
 ! --------------------------------------------------------------------------- !
  double precision :: a, v(:)
  integer :: N !size of the final voigt
  double precision, intent(out) :: F(:)
  double precision :: L(N)
  character(len=*), optional :: VoigtAlgorithm

  if (.not. present(VoigtAlgorithm)) VoigtAlgorithm="ARMSTRONG"

  SELECT CASE (VoigtAlgorithm)
!   CASE ("ARMSTRONG")
!    L = VoigtArmstrong(N, a, v)
  CASE ("HUMLICEK")
   L = VoigtHumlicek(N, a, v, F)
  CASE ("ABRAROV2010")
   L = VoigtAbrarov(N, a, v, F)
  CASE ("LookupTable")
   L = VoigtLookup(N, a, v, F)
  CASE DEFAULT
   L = VoigtHumlicek(N, a, v, F) !using Humlicek by default
  END SELECT

 RETURN
 END FUNCTION Voigt


END MODULE voigtfunctions
