module integrate

  use constantes, only: dp

  implicit none

 public :: integrate_trap,integrate_trap_log,integrate_trap_log_uniform
 public :: check_log_grid_is_uniform
 public :: gauss_legendre_nodes_weights,gauss_legendre_nodes_weights_log
 public :: integrate_gauss_legendre

 private

 abstract interface
    pure real function func(x)
      use constantes, only: dp
      real(kind=dp), intent(in) :: x
    end function func
 end interface

 interface integrate_trap
  module procedure integrate_trap_func,integrate_trap_array
 end interface

 interface integrate_trap_log
  module procedure integrate_trap_log,integrate_trap_log_func
 end interface

contains
!-------------------------------------------------------------------
! helper routine to integrate a function using the trapezoidal rule
!-------------------------------------------------------------------
pure real function integrate_trap_func(f,xmin,xmax) result(g)
 real(dp), intent(in) :: xmin,xmax
 procedure(func) :: f
 real(dp) :: fx,fprev,dx,x
 integer, parameter :: npts = 4096
 integer :: i

 g = 0.
 dx = (xmax-xmin)/(npts)
 fprev = f(xmin)
 do i=2,npts
    x = xmin + i*dx
    fx = f(x)
    g = g + 0.5*dx*(fx + fprev)
    fprev = fx
 enddo

end function integrate_trap_func

!-------------------------------------------------------------------
! helper routine to integrate a function using the trapezoidal rule
!-------------------------------------------------------------------
pure real function integrate_trap_array(n,x,f) result(g)
 integer, intent(in) :: n
 real(dp), intent(in) :: x(n),f(n)
 integer :: i

 g = 0.
 do i=2,n
    g = g + 0.5*(x(i)-x(i-1))*(f(i) + f(i-1))
 enddo

end function integrate_trap_array

!--------------------------------------------------------
!+
!  Integrate function on logarithmic grid
!  i.e. \int f(x) dx = \int x f(x) d(ln x)
!  using trapezoidal rule
!+
!--------------------------------------------------------
real function integrate_trap_log(n,x,f) result(fint)
 integer, intent(in) :: n
 real(dp), intent(in) :: x(:),f(:)
 real(dp) :: dlogx
 integer :: i

 fint = 0.
 if (n < 2) return
 do i=2,n
    dlogx = log(x(i)/x(i-1))
    fint = fint + 0.5*(f(i)*x(i) + f(i-1)*x(i-1))*dlogx
 enddo

end function integrate_trap_log

!--------------------------------------------------------
!+
!  Integrate function on evenly spaced logarithmic grid
!  i.e. \int f(x) dx = \int x f(x) d(ln x)
!  using trapezoidal rule
!+
!--------------------------------------------------------
real function integrate_trap_log_uniform(n,x,f) result(fint)
 integer, intent(in) :: n
 real(dp), intent(in) :: x(:),f(:)
 real(dp) :: dlogx
 integer :: i

 fint = 0.
 if (n < 2) return
 dlogx = log(x(2)/x(1))
 do i=2,n
    fint = fint + 0.5*(f(i)*x(i) + f(i-1)*x(i-1))*dlogx
 enddo

end function integrate_trap_log_uniform

!--------------------------------------------------------
!+
!  Integrate function on logarithmic grid
!  with a function call
!+
!--------------------------------------------------------
real function integrate_trap_log_func(n,x,f) result(fint)
 integer, intent(in) :: n
 real(dp), intent(in) :: x(:)
 procedure(func) :: f
 real :: dlogx
 integer :: i

 fint = 0.
 if (n < 2) return
 do i=2,n
    dlogx = log(x(i)/x(i-1))
    fint = fint + 0.5*(f(x(i))*x(i) + f(x(i-1))*x(i-1))*dlogx
 enddo

end function integrate_trap_log_func

!--------------------------------------------------------
!+
!  check logarithmic grid is uniform
!+
!--------------------------------------------------------
logical function check_log_grid_is_uniform(x) result(uniform)
 real(dp), intent(in) :: x(:)
 real :: dx,dxprev
 integer :: n,i

 uniform = .true.
 n = size(x)
 if (n < 2) then
    uniform = .false.
    return
 endif

 dxprev = x(2)/x(1)
 do i=3,n
    dx = x(i)/x(i-1)
    !print*,dx,dxprev,abs(dx-dxprev),epsilon(0.)
    if (abs(dx - dxprev) > 1.e-13) then
       uniform = .false.
       return
    endif
    dxprev = dx
 enddo

end function check_log_grid_is_uniform

!--------------------------------------------------------
!+
!  integrate function with Gauss-Legendre quadrature
!+
!--------------------------------------------------------
real function integrate_gauss_legendre(n,w,f) result(fint)
 integer, intent(in) :: n
 real(dp), dimension(n), intent(in) :: w, f
 integer :: i

 fint = 0.
 do i=1,n
    fint = fint + w(i)*f(i)
 enddo

end function integrate_gauss_legendre

!--------------------------------------------------------
!+
!  compute nodes and weights for Gauss-Legendre quadrature
!+
!--------------------------------------------------------
subroutine gauss_legendre_nodes_weights(x, xmin, xmax, w, use_log)
 real, parameter :: pi = 4.*atan(1.)
 real, intent(in) :: xmin, xmax
 real, dimension(:), intent(out) :: x, w
 logical, intent(in), optional :: use_log
 integer :: i,j,m,n
 real :: z, z1, pp, p1, p2, p3
 logical :: uselog

 uselog = .false.
 if (present(use_log)) uselog = use_log

 n = size(x)
 m = (n + 1) / 2
 do i = 1, m
    z = cos(pi * (i - 0.25) / (n + 0.5))
    do
       p1 = 1.0
       p2 = 0.0
       do j=1,n
          p3 = p2
          p2 = p1
          p1 = ((2.0 * j - 1.0) * z * p2 - (j - 1.0) * p3) / j
       enddo
       pp = n * (z * p1 - p2) / (z * z - 1.0)
       z1 = z
       z = z1 - p1 / pp
       if (abs(z - z1) <= 1.0e-14) exit
    enddo
    if (uselog) then
       x(i) = exp(0.5 * ((log(xmax) - log(xmin)) * (-z) + log(xmin) + log(xmax)))
       x(n + 1 - i) = exp(0.5 * ((log(xmax) - log(xmin)) * z + log(xmin) + log(xmax)))
       w(i) = 2.0 / ((1.0 - z * z) * pp * pp) * 0.5 * (log(xmax) - log(xmin)) * x(i)
    else
       x(i) = 0.5 * ((xmax - xmin) * (-z) + xmin + xmax)
       x(n + 1 - i) = 0.5 * ((xmax - xmin) * z + xmin + xmax)
       w(i) = 2.0 / ((1.0 - z * z) * pp * pp) * 0.5 * (xmax - xmin)
    endif
    w(n + 1 - i) = w(i)
 enddo

end subroutine gauss_legendre_nodes_weights

!--------------------------------------------------------
!+
!  interface to above enforcing logarithmic grid
!+
!--------------------------------------------------------
subroutine gauss_legendre_nodes_weights_log(x, xmin, xmax, w)
 real, intent(in) :: xmin, xmax
 real, dimension(:), intent(out) :: x, w

 call gauss_legendre_nodes_weights(x, xmin, xmax, w, use_log=.true.)

end subroutine gauss_legendre_nodes_weights_log

end module integrate
