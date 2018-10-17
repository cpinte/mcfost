MODULE math

  IMPLICIT NONE

  CONTAINS

   FUNCTION SQ(x) result(y)
    real(8) :: x, y
    y = x*x
   RETURN
   END FUNCTION SQ

   FUNCTION SQarr(N, x) result(y)
    integer :: N
    double precision, dimension(N) :: x
    double precision :: y(N)
    integer :: la

    do la=1,N
     y(la) = SQ(x(la))
    end do
    return
   END FUNCTION SQarr


   FUNCTION CUBE(x) result(y)
    real(8) :: x, y
    y = x*x*x
   RETURN
   END FUNCTION CUBE

   FUNCTION CUBEarr(N, x) result(y)
    integer :: N
    real(8), dimension(N) :: x
    integer :: la
    real(8), dimension(N) :: y

    do la=1, N
     y(la) = CUBE(x(la))
    end do
   RETURN
   END FUNCTION CUBEarr

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

   FUNCTION dPOWarr(N, x,a) result(y)
    ! ------------------------------
    ! double precision pow function
    ! ------------------------------
    integer :: N
    double precision :: x(N), a
    double precision, dimension(N) :: y
    integer :: la

    do la=1,N
     y(la) = dPOW(x(la),a)
    end do

   RETURN
   END FUNCTION dPOWarr

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
 Real(8), dimension(n) :: x, y
 Real(8), dimension(np) :: xp, yp
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
  Real(8), dimension(n) :: x, y
  Real(8), dimension(np) :: xp, yp
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
   REAL(8), dimension(:) :: x1a
   REAL(8), DIMENSION(:) :: ya
   REAL(8) :: x1
   REAL(8) :: y
   real(8), dimension(1) :: tmp1, tmp2, tmp3

   tmp1(1) = x1
   call bezier3_interp(size(x1a),x1a,ya,1, tmp1,tmp3)
   !call bezier2_interp(size(x1a),x1a,ya,1, tmp1,tmp3)

   y = tmp3(1)

  RETURN
  END FUNCTION

  FUNCTION interp2D(x1a,x2a,ya,x1,x2) result(y)
   ! interpolate at points x1 and x2
   ! inside the grid vectors x1a and x2a.
   ! special case for Barklem data.
   ! Using 2 x 1D Bezier cubic splines (in the x and y direction)
   ! return scalar
   REAL(8), dimension(:) :: x1a, x2a
   REAL(8), DIMENSION(:,:) :: ya
   REAL(8) :: x1,x2
   REAL(8) :: y
   real(8), dimension(1) :: tmp1, tmp2, tmp3

   INTEGER :: j
   REAL(8), DIMENSION(size(x1a)) :: ymtmp

   tmp1(1) = x1
   tmp2(1) = x2

   do j=1,size(x1a) !y-axis interpolation
    call bezier3_interp(size(x2a),&
           x2a,ya(j,:),1, tmp2, tmp3)
    ymtmp(j) = tmp3(1) !tmp3 needed to extract the scalar value
   end do
  call bezier3_interp(size(x1a), &
           x1a, ymtmp, 1, tmp1, tmp3)

   y = tmp3(1)

  RETURN
  END FUNCTION interp2D

  FUNCTION interp2Darr(N1a, x1a, N2a, x2a,ya,N1, x1, N2, x2) result(y)
   ! interpolate at vectorised points x1 and x2
   ! inside the grid vectors x1a and x2a.
   ! return vector
   integer :: N1a, N2a, N1, N2
   REAL(8) :: x1a(N1a), x2a(N2a)
   REAL(8), DIMENSION(N1a,N2a) :: ya
   REAL(8) :: x1(N1),x2(N2)
   REAL(8), dimension(N1,N2) :: y

   INTEGER :: m, n
   REAL(8), DIMENSION(N1a, N2) :: ymtmp


   do n=1,N1a ! y axis interpolation
    call bezier3_interp(N2a,&
           x2a,ya(n,:),N2, x2, ymtmp(n,:))
   end do
   do m =1,N2 ! interpolate on the x direction
   call bezier3_interp(N1a, &
           x1a, ymtmp(:,m), N1,x1, y(:,m))
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


  SUBROUTINE SORT(X,N)
    ! Modified NetlibLapack sort function
    integer, intent(in) :: N
    double precision, intent(inout) :: X(N)
    double precision :: S,T, Y(N)
    INTEGER I,J,K,L,M

    Y(:) = X(:)
    I = 1
10  K = I
20  J = I
    I = I + 1
    IF ( J .EQ. N ) GOTO 30
    IF ( X(I) .GE. X(J) ) GOTO 20
    Y(K) = I
    GOTO 10
30  IF ( K .EQ. 1 ) RETURN
    Y(K) = N + 1
40  M = 1
    L = 1
50  I = L
    IF ( I .GT. N ) GOTO 120
    S = X(I)
    J = Y(I)
    K = J
    IF ( J .GT. N ) GOTO 100
    T = X(J)
    L = Y(J)
    X(I) = L
60  IF ( S .GT. T ) GOTO 70
    Y(M) = S
    M = M + 1
    I = I + 1
    IF ( I .EQ. K ) GOTO 80
    S = X(I)
    GOTO 60
70  Y(M)= T
    M = M + 1
    J = J + 1
    IF ( J .EQ. L ) GOTO 110
    T = X(J)
    GOTO 60
80  Y(M) = T
    K = M + L - J
    I = J - M
90  M = M + 1
    IF ( M .EQ. K ) GOTO 50
    Y(M) = X(M+I)
    GOTO 90
100 X(I) = J
    L = J
110 Y(M) = S
    K = M + K - I
    I = I - M
    GOTO 90
120 I = 1
130 K = I
    J = X(I)
140 X(I) = Y(I)
    I = I + 1
    IF ( I .LT. J ) GOTO 140
    Y(K) = I
    IF ( I .LE. N ) GOTO 130
    IF ( K .EQ. 1 ) RETURN
    GOTO 40
  END SUBROUTINE SORT


  FUNCTION fact(N) result (f)
   ! Factorial function up to N = 101
   integer :: i, N
   integer(8) :: f
   f = 1
   if (N.eq.0) RETURN
   do i=1,N
    f = f * dble(i)
    if (N.gt.101) exit
   end do
  RETURN
  END FUNCTION


 SUBROUTINE SolveLinearEq(N, A, b, improve_sol)
 ! -----------------------------------------------------------
 ! Solves for a linear system equation from a squared matrix !
 ! the system to solve is Ax=b for instance A=A(Nlevel, Nlevel)
 ! and b(Nlevel) the population of each level
 ! A the matrix to decompose, note that ludcmp modify A
 ! b the right hand side vector. Note that b is modified
 ! in lubksb and contains the solution
 ! -----------------------------------------------------------
  integer, intent(in) :: N
  integer :: index(N), error, i, j
  double precision, intent(inout) :: A(N,N), b(N)
  logical, intent(in) :: improve_sol
  double precision :: A_copy(N,N), b_copy(N), residual(N), d

  if (improve_sol) then
   do i=1,N
    b_copy(i) = b(i)
    do j=1,N
     A_copy(i,j) = A(i,j)
    end do
   end do
  end if

  ! - initial solution
  CALL ludcmp(A,index,d,error)
  CALL lubksb(A,index,B)
  ! leave here if not improve_sol
  if (improve_sol) then
   do i=1,N
    residual(i) = b_copy(i)
    do j=1,N
     residual(i) = residual(i) - A_copy(i,j)*b(j)
    end do
   end do
   CALL lubksb(A,index,B)
   ! - correct the initial solution
   do i=1,N
    b(i) = b(i) + residual(i)
   end do
  end if


 RETURN
 END SUBROUTINE

  ! 3j, 6j and 9j symbols
  ! From A. A Ramos, Hazel
!----------------------------------------------------------------
! This function calculates the 3-j symbol
! J_i and M_i have to be twice the actual value of J and M
!----------------------------------------------------------------
  FUNCTION w3js(j1,j2,j3,m1,m2,m3)
   integer :: m1, m2, m3, j1, j2, j3
   integer :: ia, ib, ic, id, ie, im, ig, ih, z, zmin, zmax, jsum
   real(kind=8) :: w3js, cc, denom, cc1, cc2


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

    RETURN
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
