MODULE math

  use mcfost_env, only : dp
  use constantes, only : tiny_dp, huge_dp

  IMPLICIT NONE


  CONTAINS

!   building
!   FUNCTION overlapping_transitions(lambda, Nl, Nblue, Nl2, Nblue2) result(overlap)
!   Find where lambda(Nblue:Nblue+Nl-1)==lambda(Nblue2:Nblue2+Nl2-1)
!   meaning where two transitions overlap.
!   All transitions share the same wavelength grid, so they have the same value of lambda
!   for the same index. But they have differents Nblue, Nl
!    integer :: Nblue, Nl, Nl2, Nblue2, overlap(Nl)
!    real(kind=dp) :: lambda(:)
!    integer :: i, j
!
!    overlap(:) = 0
!
!    i_loop : do i=1,Nl
!     do j=1,Nl2
!      if (lambda(Nblue+i-1) == lambda(Nblue2+j-1)) then
!        overlap(i) = Nblue+i-1 !index on lambda grid
!        cycle i_loop
!      end if
!
!     enddo
!    enddo i_loop
!
!   RETURN
!   END FUNCTION overlapping_transitions


  FUNCTION flatten(n1, n2, M)
  !for each i in n1 write the n2 values of n1(i,:)
     		     !In the flattened array, there are:
     		     ! ilvl=1
     		     !    icell=1->ncells
     		     !                   ilvl=2
     		     !                   icell=1->Ncells
   integer :: n1, n2, i, j
   real(kind=dp) :: flatten(n1*n2), M(n1,n2)

   do i=1, n1
    !write(*,*) i
    do j=1, n2
    !write(*,*) j, n2*(i-1)+j, n2
     flatten(n2*(i-1)+j) = M(i,j)
     !n2 values per i
    enddo
   enddo

  RETURN
  END FUNCTION flatten

  FUNCTION reform(n1, n2, F)
  !reform the (n1, n2) matrix from flatten(n1,n2,F)
   integer :: n1, n2, i, j
   real(kind=dp) :: F(n1*n2), reform(n1,n2)

   do i=1, n1
    do j=1, n2
     reform(i,j) = F(n2*(i-1)+j)
    enddo
   enddo

  RETURN
  END FUNCTION reform

  FUNCTION flatten2(n1, n2, M)
  !for each j in n2 write the n1 values of n1(:,j)
     		     !In the flattened array, there are:
     		     ! icell=1
     		     !    ilvl=1->Nl
     		     !               icell=2
     		     !                   ilvl=1->Nl
   integer :: n1, n2, i, j
   real(kind=dp) :: flatten2(n1*n2), M(n1,n2)

   do i=1, n1 !Nlevel
    !write(*,*) i
    do j=1, n2 !Ncells
     !write(*,*) j, n1*(j-1)+i, n2
     flatten2(n1*(j-1)+i) = M(i,j)
     !n1 values per j
    enddo
   enddo

  RETURN
  END FUNCTION flatten2

  FUNCTION reform2(n1, n2, F)
  !reform the (n1, n2) matrix from flatten2(n1,n2,F)
   integer :: n1, n2, i, j
   real(kind=dp) :: F(n1*n2), reform2(n1,n2)

   do i=1, n1
    do j=1, n2
     reform2(i,j) = F(n1*(j-1)+i)
    enddo
   enddo

  RETURN
  END FUNCTION reform2

   FUNCTION is_nan_infinity(y) result(val)
    real(kind=dp) :: y, val

     val = 0
    if (y /= y) then
     write(*,*) "(Nan):", y
     val = 1
     return
    else if (y > 0 .and. (y==y*10)) then
     write(*,*) "(infinity):", y
     val = 2
     return
    end if

   RETURN
   END FUNCTION is_nan_infinity

   FUNCTION any_nan_infinity_matrix(y) result(val)
    real(kind=dp) :: y(:,:)
    integer :: val, i, j

     val = 0
     do i=1,size(y(:,1))
      do j=1, size(y(1,:))
       if (y(i,j) /= y(i,j)) then
        write(*,*) "(Nan):", y(i,j)
        val = 1
        return
       else if (y(i,j) > 0 .and. (y(i,j)==y(i,j)*10)) then
        write(*,*) "(infinity):", y(i,j), y(i,j)*(1+0.1)
        val = 2
        return
       end if
      end do
     end do
   RETURN
   END FUNCTION any_nan_infinity_matrix

   FUNCTION any_nan_infinity_vector(y) result(val)
    real(kind=dp) :: y(:)
    integer :: val, i

     val = 0
     do i=1,size(y)
       if (y(i) /= y(i)) then
        write(*,*) "(Nan):", y(i)
        val = 1
        return
       else if (y(i)>0 .and. (y(i)==y(i)*10)) then
        write(*,*) "(infinity):", y(i), y(i)*(1+0.1)
        val = 2
        return
       end if
     end do
   RETURN
   END FUNCTION any_nan_infinity_vector


   function cmf_to_of (Nx, y, dk)
   !Shift a function y, centered at 0 in the velocity space
   !by a delta in index of dk.
   !dk is an integer positive or negative.
   !The function y has to be linearly spaced such that y'(1) = y(1 + dk)
   !where y' is the function projected onto the observer's frame.
   ! dk is int(dv * di)
   ! with di = (index(y[1]) - index(y[0]))/(y[1]-y[0]) and dv a velocity shift.
    integer, intent(in) :: Nx, dk
    real(kind=dp), intent(in), dimension(Nx) :: y
    real(kind=dp), dimension(Nx) :: cmf_to_of
    integer :: j

    do j=1, Nx

     if ( (j + dk < 1) .or. (j + dk > Nx) ) then
     	cmf_to_of(j) = 0.0
     else
     	cmf_to_of(j) = y(j+dk)
     endif

    enddo


   return
   end function cmf_to_of


   !c'est bourrin, il n y a pas de test
   function linear_1D(N,x,y,Np,xp)
     real(kind=dp) :: x(N),y(N),linear_1D(Np), xp(Np), t
     integer :: N, Np, i, j

     do i=1,N-1
        do j=1,Np
           if (xp(j)>=x(i) .and. xp(j)<=x(i+1)) then
              t = (xp(j) - x(i)) / (x(i+1)-x(i))
              linear_1D(j) = (1.0_dp - t) * y(i)  + t * y(i+1)
           endif
        enddo
     enddo

     return

   end function linear_1D

   function linear_1D_sorted(n,x,y, np,xp)
     ! assumes that both x and xp are in increasing order
     ! We only loop once over the initial array, and we only perform 1 test per element

     integer, intent(in)                      :: n, np
     real(kind=dp), dimension(n),  intent(in) :: x,y
     real(kind=dp), dimension(np), intent(in) :: xp
     real(kind=dp), dimension(np)             :: linear_1D_sorted

     real(kind=dp) :: t
     integer :: i, j, i0, j0

     linear_1D_sorted(:) = 0._dp

     ! We do a first pass, to find the 1st index to interpolate
     ! Below x(1), we keep the values to 0 (ie no extrapolation)
     j0=np+1
     do j=1, np
        if (xp(j) > x(1)) then
           j0 = j
           exit
        endif
     enddo

     ! We perform the 2nd pass where we do the actual interpolation
     ! For points larger than x(n), value will stay at 0
     i0 = 2
     do j=j0, np
        loop_i : do i=i0, n
           if (x(i) > xp(j)) then
              t = (xp(j) - x(i-1)) / (x(i) - x(i-1))
              linear_1D_sorted(j) = (1.0_dp - t) * y(i-1)  + t * y(i)
              i0 = i
              exit loop_i
           endif
        enddo loop_i
     enddo

     return

   end function linear_1D_sorted

   function convolve(x, y, K)
   !x, y, and K have the same dimension
   !using trapezoidal rule
    real(kind=dp), intent(in) :: x(:), y(:), K(:)
    real(kind=dp) :: convolve(size(x)), dx
    integer :: i, j, Nx, shift

    Nx = size(x)
    convolve(:) = 0.0_dp


    do j=1,Nx
    	if (j==1) then
    		dx = dabs(x(j+1)-x(j))
    	else
    		dx = dabs(x(j)-x(j-1))
    	endif

    	do i=1,Nx-1 !Should not happen but here. Work in python because negative index exists
    		if (j-i < 1 .or. j-(i+1) < 1) then
    			convolve(i) = 0.0
    		else
    			convolve(i) = convolve(i) + 0.5 * dx * (y(i)*K(j-i) + y(i+1)*K(j-(i+1)))
			endif
    	enddo


    enddo


  return
  end function convolve

   FUNCTION Integrate_x(N, x, y) result(integ)
    integer :: N
    real(kind=dp), dimension(N) :: x, y
    real(kind=dp) :: integ
    integer  :: k

    integ = y(1) !0.5 * ( x(2)-x(1) ) * y(1)
    do k=2,N
     integ = integ + (x(k)-x(k-1)) * 0.5 * (y(k)+y(k-1))
    end do

   RETURN
   END FUNCTION Integrate_x

   FUNCTION Integrate_nu(N, x, y)
    integer :: N
    real(kind=dp), dimension(N) :: x, y
    real(kind=dp) :: Integrate_nu, sum1, sum2
    integer  :: k

    Integrate_nu = y(1) + y(N)
    sum1 = 0d0; sum2 = 0d0
    do k=2,N-1
     if (mod(k,2) /= 0) then
      sum1 = sum1 + y(k)
     else
      sum2 = sum2 + y(k)
     endif
    end do

    Integrate_nu = (Integrate_nu + 4*sum1 + 2*sum2) * (x(N)-x(1))/(N*3.)

   RETURN
   END FUNCTION Integrate_nu


   FUNCTION SQ(x) result(y)
    real(kind=dp) :: x, y
    y = x*x
   RETURN
   END FUNCTION SQ


   FUNCTION CUBE(x) result(y)
    real(kind=dp) :: x, y
    y = x*x*x
   RETURN
   END FUNCTION CUBE


   FUNCTION dPOW(x,a) result(y)
    ! ------------------------------
    ! real(kind=dp) pow function
    ! ------------------------------
    real(kind=dp) :: x, a, y
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
   real(kind=dp) :: y, x

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
   real(kind=dp) :: x, y

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
  real(kind=dp) :: x(n), y(n)
  real(kind=dp), dimension(n), intent(out) ::  yp
  real(kind=dp) :: der, der1, lambda, dx , dx1

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
 real(kind=dp), dimension(n), intent(in) :: x, y
 real(kind=dp), dimension(np), intent(in) :: xp
 real(kind=dp), dimension(np), intent(out) :: yp
 real(kind=dp) :: cntrl, dx, yprime(n), lambda, u(np)


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
  real(kind=dp), dimension(n), intent(in) :: x, y
  real(kind=dp), dimension(np), intent(in) :: xp
  real(kind=dp), dimension(np), intent(out) :: yp
  real(kind=dp) :: c1, c2, yprime(n), dx, u(np), mmi, mma


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
   real(kind=dp), dimension(:), intent(in) :: x1a, ya
   real(kind=dp), intent(in) :: x1
   real(kind=dp) :: y
   real(kind=dp), dimension(1) :: tmp1, tmp3

   tmp1(1) = x1
   call bezier3_interp(size(x1a),x1a,ya,1, tmp1,tmp3)
   !call bezier2_interp(size(x1a),x1a,ya,1, tmp1,tmp3)

   y = tmp3(1)

  RETURN
  END FUNCTION

  FUNCTION interp1Darr(x1a,ya,x1) result(y)
   real(kind=dp), dimension(:), intent(in) :: x1a, ya, x1
   real(kind=dp), dimension(size(x1)) :: y

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
   real(kind=dp), dimension(:), intent(in) :: x1a, x2a
   real(kind=dp), DIMENSION(:,:), intent(in) :: ya
   real(kind=dp), intent(in) :: x1,x2
   real(kind=dp) :: y
   real(kind=dp), dimension(1) :: tmp1, tmp2, tmp3

   INTEGER :: j
   real(kind=dp), DIMENSION(size(x1a)) :: ymtmp
   real(kind=dp), dimension(size(x2a)) :: yb

   tmp1(1) = x1
   tmp2(1) = x2

   do j=1,size(x1a) !y-axis interpolation
    yb = ya(j,:)
    call bezier2_interp(size(x2a),x2a,yb,1, tmp2, tmp3)
    ymtmp(j) = tmp3(1) !tmp3 needed to extract the scalar value
   end do
   call bezier2_interp(size(x1a), x1a, ymtmp, 1, tmp1, tmp3)

   y = tmp3(1)

  RETURN
  END FUNCTION interp2D

  FUNCTION interp2Darr(N1a, x1a, N2a, x2a,ya,N1, x1, N2, x2) result(y)
   integer, intent(in) :: N1a, N2a, N1, N2
   real(kind=dp), intent(in) :: x1a(N1a), x2a(N2a)
   real(kind=dp), intent(in) :: x1(N1), x2(N2)
   real(kind=dp), DIMENSION(N1a,N2a), intent(in) :: ya
   real(kind=dp)  :: y(N1,N2)
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
   real(kind=dp), INTENT(IN) :: xx
   real(kind=dp) :: gammln
   real(kind=dp) :: ser,tmp,x,y
   real(kind=dp) :: stp = 2.5066282746310005
   real(kind=dp), DIMENSION(6) :: coef = &
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

  FUNCTION locate(xx,x,mask) result(y)

  real(kind=dp), dimension(:), intent(in) :: xx
  real(kind=dp), intent(in) :: x
  logical, intent(in), dimension(:), optional :: mask
  integer :: y

  if (present(mask)) then
   y = minloc((xx-x)**2,1,mask=mask)
  else
  ! 1D array
   y = minloc((xx-x)**2,1) !(xx(:)-x)*(xx(:)-x)
  end if

  RETURN
  END FUNCTION locate


  FUNCTION fact(N) result (f)
   integer :: i, N
   real(kind=dp) :: f
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
   real(kind=dp) :: w3js, cc, denom, cc1, cc2


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
