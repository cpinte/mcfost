MODULE math

  IMPLICIT NONE


  CONTAINS
  
   FUNCTION Integrate_x(N, x, y) result(integ)
    integer :: N
    double precision, dimension(N) :: x, y
    double precision :: integ
    integer  :: k
    
    integ = 0.5 * dabs(x(2)-x(1)) * y(1)
    do k=2,N
     integ = integ + dabs(x(k)-x(k-1)) * 0.5 * (y(k)+y(k-1))
    end do
   
   RETURN
   END FUNCTION Integrate_x
   
   FUNCTION Integrate_dx(N, dx, y) result(integ)
    integer :: N
    double precision, dimension(N) :: dx, y
    double precision :: integ
    integer  :: k
    
    integ = 0.5 * dx(1) * y(1)
    do k=2,N
     integ = integ + dx(k) * 0.5 * (y(k)+y(k-1))
    end do
   
   RETURN
   END FUNCTION Integrate_dx

   FUNCTION SQ(x) result(y)
    double precision :: x, y
    y = x*x
   RETURN
   END FUNCTION SQ


   FUNCTION CUBE(x) result(y)
    double precision :: x, y
    y = x*x*x
   RETURN
   END FUNCTION CUBE


   FUNCTION dPOW(x,a) result(y)
    ! ------------------------------
    ! double precision pow function
    ! ------------------------------
    double precision :: x, a, y
    integer :: i, one
    i = 1
    y = 1.0

    if (a.eq.0.) RETURN ! x**0 = 1.
    if (a.eq.1.) then
     y = x
     RETURN
    end if

    if (abs(mod(a,1.)).gt.0.) then ! a is not an integer
     y = x**a ! a can be a negative real
     RETURN
    end if
    one = a/abs(a)
    do while (i.le.abs(int(a)))
     if (one.gt.0) y = y*x
     if (one.lt.0) y = y/(x)
     i = i+1
    end do

   RETURN
   END FUNCTION dPOW

  FUNCTION E1 (x) result (y)
   !First exponential integral
   real, dimension(6) :: a6
   real, dimension(4) :: a4, b4
   real(8) :: y, x

   data a6  /-0.57721566,  0.99999193, -0.24991055,&
             0.05519968, -0.00976004,  0.00107857 /
   data a4  / 8.5733287401, 18.0590169730, &
               8.6347608925,  0.2677737343 /
   data b4  /9.5733223454, 25.6329561486, &
              21.0996530827,  3.9584969228/

   if (x<=0.0) then
    write(*,*) "Arg of Exponential integral has to be > 0"
    write(*,*) "exiting..."
    stop
   else if (x.gt.0 .and. x.le.1.) then
    y = -log(x) + a6(1) + x*(a6(2) + x*(a6(3) + &
       x*(a6(4) + x*(a6(5) + x*a6(6)))))
   else if (x.gt.1 .and. x.le.80.) then
    y  = a4(4)/x +  a4(3) + x*(a4(2) + x*(a4(1) + x))
    y = y / (&
       b4(4) + x*(b4(3) + x*(b4(2) + x*(b4(1) + x))) &
        )
    y = y * exp(-x);
   else
    y = 0.0
   end if

  RETURN
  END FUNCTION E1

  FUNCTION E2(x) result (y)
  ! second exponential integral, using recurrence relation
  ! En+1 = 1/n * ( exp(-x) -xEn(x))
   real(8) :: x, y

   if (x.le.0.) then
    write(*,*) "Arg of Exponential integral has to be > 0"
    write(*,*) "exiting..."
    stop
   else if (x.gt.0 .and. x.le.80.) then
    y = 1./1. * (exp(-x) - x* E1(x))
   else
    y = 0.
   end if

  RETURN
  END FUNCTION
  
 FUNCTION Gaunt_bf(Nl,u, n_eff) result(GII)
  ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
  ! See also Menzel & Pekeris 1935
  integer :: Nl
  double precision :: n_eff, GII(Nl)
  double precision :: u(Nl), n23

  n23 = n_eff**(-6.6666666666667d-1) !n^-2/3

  GII = 1.+0.1728 * n23 * (u+1)**(-6.6666666666667d-1) * &
      (u-1) - 0.0496*n23*n23 * (u+1)**(-4d0/3d0) * (u*u + 4./3. * u + 1) !+...

  if (MINVAL(GII) < 0d0) then
   !write(*,*) "Warning, Gaunt factor is less than 0"
   !write(*,*) Nl, minval(u), maxval(u), n_eff, minval(GII)
   GII = dabs(GII)
  end if

  if (MAXVAL(GII) > 2d0) then
    write(*,*) "Bound-free Gaunt factor gII = ", MAXVAL(GII)
  end if

 RETURN
 END FUNCTION Gaunt_bf

 FUNCTION Gaunt_ff(N, lambda, charge, T) result(GIII)
  use constant, only : HPLANCK, CLIGHT, NM_TO_M, KBOLTZMANN, E_RYDBERG
 ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
 !
 ! Note: There is a problem with this expansion at higher temperatures
 ! (T > 3.0E4 and longer wavelengths (lambda > 2000 nm). Set to
 ! 1.0 when the value goes below 1.0
  integer, intent(in) :: N
  double precision, intent(in) :: lambda(N)
  double precision, dimension(N) :: x, x3, y, GIII
  double precision :: charge, T

  x = ((HPLANCK * CLIGHT)/(lambda * NM_TO_M)) / &
       (E_RYDBERG * (charge)**(2d0))
  x3 = (x**(3.3333333d-1))
  y  = (2.0 * lambda * NM_TO_M * KBOLTZMANN*T) / &
       (HPLANCK*CLIGHT)

  gIII = 1.0 + 0.1728*x3 * (1.0 + y) - &
        0.0496*(x3*x3) * (1.0 + (1.0 + y)*0.33333333*y)

  where (GIII <= 1d0) gIII = 1.

  if (MAXVAL(GIII) > 2d0) write(*,*) "free-free Gaunt factor gIII = ", &
      MAXVAL(gIII)

 RETURN
 END FUNCTION Gaunt_ff

  SUBROUTINE cent_deriv(n,x,y,yp)
  integer :: n, k
  Real(8) :: x(n), y(n)
  real(8), dimension(n), intent(out) ::  yp
  real(8) :: der, der1, lambda, dx , dx1

  do k = 2, n - 1
  dx = x(k) - x(k-1)
  der = (y(k) - y(k-1)) / dx

  dx1 = x(k+1) - x(k)
  der1 = (y(k+1) - y(k)) / dx1

  if(der*der1 .gt. 0.0d0) then
  lambda = (1.d0 + dx1 / (dx1 + dx)) / 3.d0
   yp(k) = (der / (lambda * der1 + (1.d0 - lambda) * der)) *der1
  Else
  yp(k) = 0.0d0
  end if
 enddo

 yp(1) =  (y(1) - y(2)) / (x(1) - x(2))
 yp(n) = (y(n-1) - y(n)) / (x(n-1) - x(n))

 RETURN
 END SUBROUTINE cent_deriv


 SUBROUTINE bezier2_interp(n, x, y, np, xp, yp)
 !
 !
 ! POINTS OUTSIDE THE DOMAIN ARE ATTRIBUTED TO THE OUTER INTIAL VALUES
 Implicit None
 Integer :: n, np, k
 Real(8), dimension(n), intent(in) :: x, y
 Real(8), dimension(np), intent(in) :: xp
 Real(8), dimension(np), intent(out) :: yp
 Real(8) :: cntrl, dx, yprime(n), lambda, u(np)


 !
 ! Compute derivatives
 !
 cntrl = 0
 yprime = 0
 call cent_deriv(n, x, y, yprime)

 do k = 2, n
  dx =  x(k) - x(k-1)


 cntrl = 0.5d0 * (y(k) - 0.5d0*dx * yprime(k) + y(k-1) + 0.5d0*dx * yprime(k-1))
 if(cntrl .GT. max(y(k),y(k-1)) .OR. cntrl .LT. min(y(k),y(k-1))) cntrl = y(k-1)

 where(xp .LT. x(k) .AND. xp .GE. x(k-1))
 u = (xp - x(k-1)) / dx
 yp = y(k-1) * (1.d0 - u)**2 + y(k) * u**2 + 2.d0 * cntrl * u * (1.d0 - u)
 End where
 end do

 !
 ! Points outside the domain
 !
 where(xp .GE. x(n))
 yp = y(n)
 end where
 where(xp .LE. x(1))
 yp = y(1)
 end where

 RETURN
 END SUBROUTINE bezier2_interp

 SUBROUTINE bezier3_interp(n, x, y, np, xp, yp)
  Integer :: n, np, k
  Real(8), dimension(n), intent(in) :: x, y
  Real(8), dimension(np), intent(in) :: xp
  Real(8), dimension(np), intent(out) :: yp
  Real(8) :: c1, c2, yprime(n), dx, u(np), mmi, mma


  c1 = 0
  c2 = 0
  yprime = 0
  !
  ! Compute derivatives
  !
  call cent_deriv(n, x, y, yprime)
  !
  do k = 2, n
   dx =  x(k) - x(k-1)

   c1 = y(k-1) + dx * yprime(k-1) / 3.0d0
   c2 = y(k) - dx * yprime(k) / 3.0d0

   mmi = min(y(k),y(k-1))
   mma = max(y(k),y(k-1))
   if(c1 .LT. mmi .OR. c1 .GT. mma) c1 = y(k-1)
   if(c2 .LT. mmi .OR. c2 .GT. mma) c2 = y(k)

   where(xp .LT. x(k) .AND. xp .GE. x(k-1))
    u = (xp - x(k-1)) / dx
    yp = y(k-1) * (1.d0 - u)**3 + y(k) * u**3 + &
    3.d0 * c1 * u * (1.d0 - u)**2 + 3.0d0 * c2 * u**2 * (1.0d0 - u)
   End where

  end do

  !
  ! Points outside the domain
  !
  where(xp .GE. x(n))
  yp = y(n)
  end where
  where(xp .LE. x(1))
  yp = y(1)
  end where
  RETURN
  END SUBROUTINE bezier3_interp

  !One point Bezier interpolation, wrapper around BezierN_interp
  ! return scallar
  FUNCTION interp1D(x1a,ya,x1) result(y)
   REAL(8), dimension(:), intent(in) :: x1a, ya
   REAL(8), intent(in) :: x1
   REAL(8) :: y
   real(8), dimension(1) :: tmp1, tmp3

   tmp1(1) = x1
   call bezier3_interp(size(x1a),x1a,ya,1, tmp1,tmp3)
   !call bezier2_interp(size(x1a),x1a,ya,1, tmp1,tmp3)

   y = tmp3(1)

  RETURN
  END FUNCTION

  FUNCTION interp1Darr(x1a,ya,x1) result(y)
   REAL(8), dimension(:), intent(in) :: x1a, ya, x1
   REAL(8), dimension(size(x1)) :: y

   call bezier3_interp(size(x1a),x1a,ya,size(x1), x1,y)
   !call bezier2_interp(size(x1a),x1a,ya,1, tmp1,tmp3)

  RETURN
  END FUNCTION

  FUNCTION interp2D(x1a,x2a,ya,x1,x2) result(y)
   ! interpolate at points x1 and x2
   ! inside the grid vectors x1a and x2a.
   ! special case for Barklem data.
   ! Using 2 x 1D Bezier cubic splines (in the x and y direction)
   ! return scalar
   REAL(8), dimension(:), intent(in) :: x1a, x2a
   REAL(8), DIMENSION(:,:), intent(in) :: ya
   REAL(8), intent(in) :: x1,x2
   REAL(8) :: y
   real(8), dimension(1) :: tmp1, tmp2, tmp3

   INTEGER :: j
   REAL(8), DIMENSION(size(x1a)) :: ymtmp
   real(8), dimension(size(x2a)) :: yb

   tmp1(1) = x1
   tmp2(1) = x2

   do j=1,size(x1a) !y-axis interpolation
    yb = ya(j,:)
    call bezier3_interp(size(x2a),x2a,yb,1, tmp2, tmp3)
    ymtmp(j) = tmp3(1) !tmp3 needed to extract the scalar value
   end do
   call bezier3_interp(size(x1a), x1a, ymtmp, 1, tmp1, tmp3)

   y = tmp3(1)

  RETURN
  END FUNCTION interp2D

  FUNCTION interp2Darr(N1a, x1a, N2a, x2a,ya,N1, x1, N2, x2) result(y)
   integer, intent(in) :: N1a, N2a, N1, N2
   double precision, intent(in) :: x1a(N1a), x2a(N2a)
   double precision, intent(in) :: x1(N1), x2(N2)
   double precision, DIMENSION(N1a,N2a), intent(in) :: ya
   double precision  :: y(N1,N2)
   integer :: i, j
   
   do i=1,N1
    do j=1,N2
     y(i,j) = Interp2D(x1a, x2a, ya,x1(i),x2(j))
    end do
   end do
   
   RETURN
  END FUNCTION interp2Darr

  FUNCTION gammln(xx)
   IMPLICIT NONE
   INTEGER :: i
   REAl(8), INTENT(IN) :: xx
   REAL(8) :: gammln
   REAL(8) :: ser,tmp,x,y
   REAL(8) :: stp = 2.5066282746310005
   REAL(8), DIMENSION(6) :: coef = &
      (/76.18009172947146,&
       -86.50532032941677,24.01409824083091,&
       -1.231739572450155,0.1208650973866179e-2,&
      -0.5395239384953e-5/)
   x=xx
   tmp=x+5.5
   tmp=(x+0.5)*log(tmp)-tmp
   ser=1.000000000190015
   y=x
   do i=1,size(coef)
    y=y+1.0
    ser=ser+coef(i)/y
   end do
   gammln=tmp+log(stp*ser/x)
  RETURN
  END FUNCTION gammln
  
  FUNCTION locate(xx,x) result(y)

  double precision, dimension(:), intent(in) :: xx
  double precision, intent(in) :: x
  integer :: y

  
  ! 1D array
  y = minloc((xx-x)**2,1) !(xx(:)-x)*(xx(:)-x)

  RETURN
  END FUNCTION locate


  FUNCTION fact(N) result (f)
   integer :: i, N
   double precision :: f
   f = 1d0
   if (N.eq.0) RETURN
   do i=1,N
    f = f * real(i,kind=8)
    if (N.gt.301) exit
   end do
  RETURN
  END FUNCTION


  ! 3j, 6j and 9j symbols
  ! From A. A Ramos, Hazel
!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
  FUNCTION w3js(j1,j2,j3,m1,m2,m3)
   integer :: m1, m2, m3, j1, j2, j3
   integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
   double precision :: w3js, cc, denom, cc1, cc2


   w3js = 0.d0
   if (m1+m2+m3 /= 0) RETURN
   ia = j1 + j2
   if (j3 > ia) RETURN
   ib = j1 - j2
   if (j3 < abs(ib)) RETURN
   if (abs(m1) > j1) RETURN
   if (abs(m2) > j2) RETURN
   if (abs(m3) > j3) RETURN

   jsum = j3 + ia
   ic = j1 - m1
   id = j2 - m2

   ie = j3 - j2 + m1
   im = j3 - j1 - m2
   zmin = max0(0,-ie,-im)
   ig = ia - j3
   ih = j2 + m2
   zmax = min0(ig,ih,ic)
   cc = 0.d0

   do z = zmin, zmax, 2
    denom = fact(z/2)*fact((ig-z)/2)*&
          fact((ic-z)/2)*fact((ih-z)/2)*&
          fact((ie+z)/2)*fact((im+z)/2)
    if (mod(z,4) /= 0) denom = -denom
    cc = cc + 1.d0/denom
   end do

   cc1 = fact(ig/2)*fact((j3+ib)/2)*&
        fact((j3-ib)/2)/fact((jsum+2)/2)
   cc2 = fact((j1+m1)/2)*fact(ic/2)*&
        fact(ih/2)*fact(id/2)*fact((j3-m3)/2)*&
        fact((j3+m3)/2)
   cc = cc * sqrt(1.d0*cc1*cc2)
   if (mod(ib-m3,4) /= 0) cc = -cc
   w3js = cc
   if (abs(w3js) < 1.d-8) w3js = 0.d0
1000        RETURN
   END FUNCTION w3js


!----------------------------------------------------------------
! This function calculates the 6-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
   FUNCTION w6js(j1,j2,j3,l1,l2,l3)
    integer :: j1,j2,j3,l1,l2,l3
    integer :: ia, ib, ic, id, ie, iif, ig, ih, sum1, sum2, sum3, sum4
    integer :: w, wmin, wmax, ii, ij, ik
    real(kind=8) :: w6js, omega, theta1, theta2, theta3, theta4, theta, denom

    w6js = 0.d0
    ia = j1 + j2
    if (ia < j3) RETURN
    ib = j1 - j2
    if (abs(ib) > j3) RETURN
    ic = j1 + l2
    if (ic < l3) RETURN
    id = j1 - l2
    if (abs(id) > l3) RETURN
    ie = l1 + j2
    if (ie < l3) RETURN
    iif = l1 - j2
    if (abs(iif) > l3) RETURN
    ig = l1 + l2
    if (ig < j3) RETURN
    ih = l1 - l2
    if (abs(ih) > j3) RETURN
    sum1=ia + j3
    sum2=ic + l3
    sum3=ie + l3
    sum4=ig + j3
    wmin = max0(sum1, sum2, sum3, sum4)
    ii = ia + ig
    ij = j2 + j3 + l2 + l3
    ik = j3 + j1 + l3 + l1
    wmax = min0(ii,ij,ik)
    omega = 0.d0
    do w = wmin, wmax, 2
     denom = fact((w-sum1)/2)*fact((w-sum2)/2)*fact((w-sum3)/2)&
       *fact((w-sum4)/2)*fact((ii-w)/2)*fact((ij-w)/2)&
       *fact((ik-w)/2)
     if (mod(w,4) /= 0) denom = -denom
     omega = omega + fact(w/2+1) / denom
    end do
    theta1 = fact((ia-j3)/2)*fact((j3+ib)/2)*&
             fact((j3-ib)/2)/fact(sum1/2+1)
    theta2 = fact((ic-l3)/2)*fact((l3+id)/2)*&
             fact((l3-id)/2)/fact(sum2/2+1)
    theta3 = fact((ie-l3)/2)*fact((l3+iif)/2)*&
             fact((l3-iiF)/2)/fact(sum3/2+1)
    theta4 = fact((ig-j3)/2)*fact((j3+ih)/2)*&
             fact((j3-ih)/2)/fact(sum4/2+1)
    theta = theta1 * theta2 * theta3 * theta4
    w6js = omega * sqrt(theta)
    write(*,*) "theta", theta
    if (abs(w6js) < 1.d-8) w6js = 0.d0

1000        RETURN
  END FUNCTION w6js

!----------------------------------------------------------------
! This function calculates the 9-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
  FUNCTION w9js(j1,j2,j3,j4,j5,j6,j7,j8,j9)
  integer :: j1,j2,j3,j4,j5,j6,j7,j8,j9
  integer :: i, kmin, kmax, k
  real(kind=8) :: x, s, x1, x2, x3, w9js

  kmin = abs(j1-j9)
  kmax = j1 + j9
  i = abs(j4-j8)
  if (i > kmin) kmin = i
  i = j4 + j8
  if (i < kmax) kmax = i
  i = abs(j2-j6)
  if (i > kmin) kmin = i
  i = j2 + j6
  if (i < kmax) kmax = i
  x = 0.d0
  do k = kmin, kmax, 2
   s = 1.d0
   if (mod(k,2) /= 0) s = -1.d0
   x1 = w6js(j1,j9,k,j8,j4,j7)
   x2 = w6js(j2,j6,k,j4,j8,j5)
   x3 = w6js(j1,j9,k,j6,j2,j3)
   x = x + s*x1*x2*x3*dfloat(k+1)
  end do
  w9js = x
 RETURN
 END FUNCTION w9js



END MODULE math
