module reconstruct_from_moments

  use mcfost_env, only : dp

  implicit none
  real, parameter :: pi = 4.*atan(1.)

  public :: reconstruct,reconstruct_maxent,fsolve_error,integrand

  private

contains
!
! reconstruct a function from a limited number
! of specified moments mu_n, corresponding to
!
!  mu_n  = \int f(x) x^n dx
!
! Arguments:
!  IN:
!    mu(:) : array of moments of arbitrary length
!    x(:)  : grid of positions on which to evaluate function
!  OUT:
!    f(size(x)) : reconstructed function evaluated
!                 on grid positions
!    lambsol(size(mu)) : parameters used in max entropy reconstruction
!    ierr : error code (>0 if something went wrong)
!
subroutine reconstruct(mu,x,f,lambsol,ierr,lambguess)
 real(dp), intent(in) :: mu(:)
 real(dp), intent(in) :: x(:)
 real(dp), intent(out) :: f(size(x))
 real(dp), intent(out) :: lambsol(size(mu))
 integer, intent(out) :: ierr
 real(dp), intent(in), optional :: lambguess(size(mu))

 ! in principle can choose method here
 if (present(lambguess)) then
    call reconstruct_maxent(mu,x,f,lambsol,ierr,lambguess)
 else
    call reconstruct_maxent(mu,x,f,lambsol,ierr,lambguess)
 endif

end subroutine reconstruct

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_maxent(mu,x,f,lambsol,ierr,lambguess)
 real(dp), intent(in) :: mu(:)
 real(dp), intent(in) :: x(:)
 real(dp), intent(out) :: f(size(x))
 real(dp), intent(out) :: lambsol(size(mu))
 integer, intent(out) :: ierr
 real(dp), intent(in), optional :: lambguess(size(mu))
 integer :: n_moments
 real(dp), parameter :: tol = 1.e-6
 real(dp) :: lsum(size(mu))

 lambsol = 0.
 ierr = 0
 f = 0.
 n_moments = size(mu)

 ! initial guesses for Lagrange multipliers
 if (present(lambguess)) then
    ! use guesses passed as arguments (e.g. from previous attempt)
    lambsol = lambguess
 else
    ! use predefined guess
    lambsol(1) = -log(sqrt(2.*pi))
 endif

 call fsolve(residual,n_moments,lambsol,lsum,tol,ierr)
 f = integrand(x,lambsol,n_moments)

contains
!
!  residual of the moment approximation function
!
!    Calculates the residual of the moment approximation function.
!
!    Parameters:
!        lamb (array): an array of Lagrange constants used to approximate the distribution
!        x (array):
!        k (integer): order of the moment
!        mu (array): an array of the known moments needed to approximate the distribution function
!
!    Returns:
!        rhs: the integrated right hand side of the moment approximation function
!
subroutine residual(n,lamb,l_sum)
 use utils, only:integrate_trap
 integer, intent(in) :: n
 real(dp), intent(in)  :: lamb(n) ! guess for  solution
 real(dp), intent(out) :: l_sum(n) ! function evaluated for given lamb

 integer :: k

 do k=1,n
    !l_sum(k) = sum(integrand(x,lamb,k-1)) - mu(k)
    l_sum(k) = integrate_trap(size(x),x,integrand(x,lamb,k-1)) - mu(k)
 enddo
 !print*,' residuals are: ',l_sum(:)

end subroutine residual

end subroutine reconstruct_maxent

!
! compute the integrand (bit inside the integral) of the kth moment
!
!  IN:
!    x   : grid on which to evaluate the integrand
!    lamb: array of Lagrange multipliers
!    k   : order of moment to calculate (e.g. 0, 1, 2...)
!
function integrand(x,lamb,k) result(y)
 real(dp), intent(in) :: x(:)
 real(dp), intent(in) :: lamb(:)
 real(dp) :: y(size(x))  ! result
 integer, intent(in) :: k
 real(dp) :: xi(size(lamb))
 integer :: i,j

 do i=1,size(x)
    do j=1,size(lamb)
       xi(j) = x(i)**(j-1)
    enddo
    y(i) = x(i)**k * exp(dot_product(lamb,xi))
 enddo

end function integrand

!
! print the error message corresponding to the error code
!
function fsolve_error(ierr) result(string)
 integer, intent(in) :: ierr
 character(len=62) :: string

 select case(ierr)
 case(0)
     string = 'improper input parameters'
 case(1)
     string = 'relative error between X and solution is at most tol'
 case(2)
     string = 'number of calls to residual function exceeded 200*(N+1)'
 case(3)
     string = 'TOL is too small. No further improvement in solution possible'
 case(4)
     string = 'the iteration is not making good progress'
 case default
    string = ''
 end select

end function fsolve_error

end module reconstruct_from_moments
