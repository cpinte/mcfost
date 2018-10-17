MODULE voigtfunctions

 use constant
 use math, only : SQ

 IMPLICIT NONE

 integer, parameter :: NC=34, NT=10, NGR=31
 double precision, parameter :: EXPMAX=70.,TWOSQRTPI=1.12837917,&
     C0=0.0897935610625833,C1=29.608813203268, THREEPI=9.42477796076938,&
     INVSQRTPI=0.564189583547756
 double precision, parameter :: TINY=1d-8
 double precision, dimension(NT) :: t, w

 data t / 0.2453407083, 0.7374737285, 1.2340762153, &
         1.7385377121, 2.2549740020, 2.7888060584, &
         3.3478545673, 3.9447640401, 4.6036824495, &
         5.3874808900 /

 data w / 4.6224366960d-01, 2.8667550536d-01, &
         1.0901720602d-01, 2.4810520887d-02, &
         3.2437733422d-03, 2.2833863601d-04, &
         7.8025564785d-06, 1.0860693707d-07, &
         4.3993409922d-10, 2.2293936455d-13 /


 CONTAINS

 ! VoigtAtmostrong function generators
 FUNCTION VoigtK1(a,v) result(L)
  double precision :: a, v, L, c(NC)
  integer :: n
  double precision :: a2, v2, u1, dn01, dn02, dn
  double precision :: v2i, funct, an, q, g, coef, bn01, bn02, bn, v1



  DATA c / 0.1999999999972224, -0.1840000000029998,            &
           0.1558399999965025, -0.1216640000043988,            &
           0.0877081599940391, -0.0585141248086907,            &
           0.0362157301623914, -0.0208497654398036,            &
           0.0111960116346270, -0.56231896167109d-02,          &
           0.26487634172265d-02, -0.11732670757704d-02,        &
           0.4899519978088d-03, -0.1933630801528d-03,          &
           0.722877446788d-04, -0.256555124979d-04,            &
           0.86620736841d-05, -0.27876379719d-05,              &
           0.8566873627d-06, -0.2518433784d-06,                &
           0.709360221d-07, -0.191732257d-07, 0.49801256d-08,  &
           -0.12447734d-08, 0.2997777d-09, -0.696450d-10,      &
           0.156262d-10, -0.33897d-11, 0.7116d-12, -0.1447d-12,&
           0.285d-13, -0.55d-14, 0.10d-14, -0.2d-15 /

  a2 = a*a
  v2 = v*v

  if ((v2-a2).gt.EXPMAX) then
   u1 = 0.
  else
   u1 = dexp(a2-v2) * cos(2.*v*a)
  end if

  if (v.gt.5.) then
   v2i = 1./v2
   dn01 = -v2i * (0.5 + v2i*(0.75 + v2i*(1.875 + &
              v2i*(6.5625 + v2i*(29.53125 + &
              v2i*(1162.4218 + v2i*1055.7421))))))
   dn02 = (1.0 - dn01) / (2.0 * v)
  else
   bn01 = 0.0
   bn02 = 0.0
   v1   = v / 5.0
   coef = 4.0 * v1*v1 - 2.0
   do n = NC, 1, -1 ! for(n=NC-1,n>=0,n--)
    bn = coef*bn01 - bn02 + c(n)
    bn02 = bn01
    bn01 = bn
   end do
   dn02 = v1*(bn - bn02)
   dn01 = 1.0 - 2.0*v*dn02;
  end if
  funct = a*dn01
  if (a.gt.TINY) then
   q = 1.
   an = a
   do n=3,50
    dn = (v*dn01 + dn02) * (-2./n)
    dn02 = dn01
    dn01 = dn
    if (mod(n, 1).gt.0) then ! if (n%2)
     q = -q
     an = an*a2
     g = dn*an
     funct = funct + q*g
     if (dabs(g/funct).le.TINY) then
      L = (u1 - TWOSQRTPI*funct)
      RETURN
     end if
    end if
   end do
  end if

  L = (u1 - TWOSQRTPI*funct)
 RETURN
 END FUNCTION VoigtK1

 FUNCTION VoigtK2(a,v) result(L)
  double precision :: a, v, L
  double precision :: g, r, s, a2
  integer :: n
  a2 = a*a
  g = 0.
  do n=1,NT
   r = t(n) - v
   s = t(n) + v
   ! beware here, by atan I mean atan(x)
   ! in the interval [-pi/2,+pi/2] radians.
   ! and log is natural (base e)
   g = g+(4.0*t(n)*t(n) - 2.0) * (r*atan(r/a) + s*atan(s/a) -&
          0.5*a*(log(a2 + r*r) + log(a2 + s*s))) * w(n)
  end do
  L = g/PI
 RETURN
 END FUNCTION VoigtK2

 FUNCTION VoigtK3(a,v) result(L)
  double precision :: a, v, L
  integer :: n
  double precision :: g, a2
  g = 0.
  a2 = a*a

  do n=1,NT
   g = g+(1.0/(SQ(v - t(n)) + a2) + 1.0/(SQ(v + t(n)) + a2)) * w(n)
  end do
  L = (a*g) / PI
 RETURN
 END FUNCTION VoigtK3

 FUNCTION VoigtArmstrong(a, v) result(L)
 ! Armstrong 1967, JQSRT 7, pp. 61-88
 ! (slow for damping parameters larger than 1.5, accurate to
 ! 6 significant figures).
  double precision :: L
  double precision :: a, v
  logical :: times_minus1=.false.

  !!!!!!!!!!!!!!!!!!!
  if (v.lt.0.0) then
     v=-v
     times_minus1=.true. !trick here
  end if
  ! A CHANGER CAR ICI LA VALEUR DE V est transformé si v < 0
  ! peut casser la grille d'entrée si v est une variable d'un objet
  ! J'utilise un petit trick avec times_minus1, mais je comprends
  ! pas trop ce qu'il se passe. Surement une connerie de Fortran...
  !!!!!!!!!!!!!!!!!!!

  if (((a.lt.1.).and.(v.lt.4.0)).or.((a.lt.1.8/(v+1.)))) then
  !write(*,*) "VoigtK1 ","a=",a,"v=",v
   L = VoigtK1(a,v)
  else if ((a.lt.2.5).and.(v.lt.4.)) then
  !write(*,*) "VoigtK2 ","a=",a,"v=",v
   L = VoigtK2(a,v)
  else
  !write(*,*) "VoigtK3 ","a=",a,"v=",v
   L = VoigtK3(a,v)
  end if

  if (times_minus1) then
   v = -v
   times_minus1 = .false. !trick here
  end if

 RETURN
 END FUNCTION VoigtArmstrong

 FUNCTION VoigtHui(a, v, F) result(L)
 ! Hui, Armstrong & Wray 1978, JQSRT 19, pp. 509-516
 ! (same speed in whole parameter space, faster than Armstrong for
 ! large a, otherwise comparable, only 1% accurate for larger v).
  double precision :: a, v, L
  double precision, intent(out) :: F
  
  L = 0d0
  F = 0d0
  write(*,*) "Not implemented yet, exiting..."
  stop
 RETURN
 END FUNCTION VoigtHui

 FUNCTION VoigtRybicki(a, v) result(L)
 ! George Rybicki's accurate (Ref)
 ! (accurate to at least 8 significant figures, but a factor 2
 ! slower than Armstrong and Hui et al.)
  double precision :: a, v, L, c(NGR)
  integer :: m, n
  double precision :: a1, a2, b1, b2, e, s, t, zi, zr 
  
  m = -15
  do n=1,NGR
    c(n) = C0 * dexp(-(real(m*m,kind=8))/9.)
    m = m+1  
  end do 
  if (a.eq.0.0) then 
   L = dexp(-v*v) !Doppler profile
   RETURN
  end if 
  
  a1 = 3. * a
  a2 = a*a
  e = dexp(-THREEPI*a)
  if (a.lt.0.) then
   zr = 0.5*(e+1./e)*cos(THREEPI*v)
   zi = 0.5*(e-1./e)*sin(THREEPI*v)
   L = INVSQRTPI*dexp(a2-v*v)*cos(2*a*v)
  else
   zr = e*cos(THREEPI*v)
   zi = e*sin(THREEPI*v)
   L = 0.  
  end if
  b1 = (1-zr)*a*1.5
  b2 = -zi
  s = -8. - 1.5*v 
  t = s*s + 2.25 * a2
  do n=1,NGR
   t = t+s+0.25
   s = s+0.5
   b1 = a1-b1 
   b2 = -b2
   if (t.gt.2.5d-12) then
     L = L + c(n)*(b1+b2*s)/t
   else 
     L = L - c(n)*s*C1
   end if 
  end do
  L = L * SQRTPI
 RETURN
 END FUNCTION VoigtRybicki

 SUBROUTINE Humlicek(a, v, W)
  ! W = L + i*F
  double precision, intent(in) :: a, v
  complex(kind=8), intent(out) :: W
  double precision :: s
  complex(kind=8) :: z, u

  z = cmplx(a,-v)
  s = abs(v) + a
  if (s >= 15.0) then
   ! Approximation in region I
   W = (z * 0.5641896) / (0.5 + (z * z))
  else if (s >= 5.5) then
   ! Approximation in region II
   u = z * z
   W = (z * (1.410474 + u*0.5641896)) / (0.75 + (u*(3.0 + u)))
 else if (a >= 0.195*abs(v) - 0.176) then
  ! Approximation in region III
  W = (16.4955 + z*(20.20933 + z*(11.96482 + z*(3.778987 + &
  0.5642236*z)))) / &
  (16.4955 + z*(38.82363 + z*(39.27121 + z*(21.69274 + &
  z*(6.699398 + z)))))
 else
  ! Approximation in region IV
  u = z * z
  W = exp(u) - (z*(36183.31 - u*(3321.99 - u*(1540.787 - &
   u*(219.031 - u*(35.7668 - u*(1.320522 - u*0.56419)))))) / &
(32066.6 - u*(24322.84 - u*(9022.228 - u*(2186.181 - &
   u*(364.2191 - u*(61.57037 - u*(1.841439 - u))))))))
  endif
 RETURN
 END SUBROUTINE Humlicek

 FUNCTION VoigtHumlicek(a,v, F) result(L)
 ! Humlicek 1982, JQSRT 27, p. 437
 ! Relative accuracy 1.0E-04. Also calculates Faraday-Voigt
 ! function needed in Stokes radiative transfer.
  double precision :: a, v, L
  double precision, intent(out) :: F
  complex(kind=8) :: W

  ! real of W is L, and imag is F
  CALL Humlicek(a, v, W)

  F = imag(W)
  L = real(W) !factor 2 here or not ?
 RETURN
 END FUNCTION VoigtHumlicek

 FUNCTION VoigtLookup(a, v) result(L)
 ! Voigt function from file or lookup table
  double precision :: a, v, L
  
  L = 0d0
  write(*,*) "Not implemented yet, exiting..."
  stop
 RETURN
 END FUNCTION VoigtLookup

 FUNCTION VoigtLookup2(a, v, F) result(L)
 ! Voigt function from file or lookup table
 ! including dispersion profile
  double precision :: a, v, L
  double precision, intent(out) :: F
  
  L = 0d0
  F = 0d0
  write(*,*) "Not implemented yet, exiting..."
  stop
 RETURN
 END FUNCTION VoigtLookup2

 ! FUNCTION Voigt[...] (a, v, F) result(L)
 !  double precision :: a, v, F, L
 !  ! ...
 ! RETURN
 ! END FUNCTION Voigt[...]

 FUNCTION Voigt(a, v, F, VoigtAlgorithm) result (L)
  ! RETURN THE VOIGT FUNCTION FOR A LINE PROFILE
  ! ACCORDING TO A DETERMINED ALGORITHM
  ! L-> Voigt function, F->Dispersion profile
  ! if desired (when polarisation is desired in fact)
  ! VoigtAlgorithm-> aglorithm for the voigt function
  double precision :: a, v, L
  double precision, intent(out) :: F
  character(len=10), optional :: VoigtAlgorithm
  character(len=10) :: Algorithm

  if (.not. present(VoigtAlgorithm)) then
   Algorithm="ARMSTRONG "
  else
   Algorithm = VoigtAlgorithm
  end if

  SELECT CASE (Algorithm)
  CASE ("ARMSTRONG ")
   L = VoigtArmstrong(a, v)
  CASE ("HUMLICEK  ")
   !Zeeman Case
   L = VoigtHumlicek(a,v, F)
  CASE ("RYBICKI   ")
   L = VoigtRybicki(a, v)
   ! -- Not implemented yet --!
  CASE ("HUI_ETAL  ")
   L = VoigtHui(a, v, F)
  CASE ("LOOKUP    ")
   L = VoigtLookup(a, v)
  CASE ("LOOKUP_2  ")
   L = VoigtLookup2(a, v, F)
  ! -- Not implemented yet --!
  CASE DEFAULT
   write(*,*) "Algorithm ",Algorithm, &
       " for voigt function unkown!"
   write(*,*) "len=",len(Algorithm)
   write(*,*) "Exciting..."
   stop
  END SELECT

 RETURN
 END FUNCTION Voigt


END MODULE voigtfunctions
