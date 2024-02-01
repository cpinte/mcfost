module reconstruct_from_moments

  use constantes

 implicit none

 public :: reconstruct,reconstruct_maxent,fsolve_error,integrand
 public :: reconstruct_gamma_dist,gamma_func,gamma_func_from_moments,gamma_func_moment

 logical, private :: use_log_grid = .false.
 logical, private :: log_grid_is_uniform = .false.
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
 real, intent(in) :: mu(:)
 real(dp), intent(in) :: x(:)
 real, intent(out) :: f(size(x))
 real(dp), intent(out) :: lambsol(size(mu))
 integer, intent(out) :: ierr
 real(dp), intent(in), optional :: lambguess(size(mu))
 real(dp) :: lsum(size(mu))

 ! in principle can choose method here
 if (present(lambguess)) then
    call reconstruct_maxent(mu,x,f,lambsol,lsum,ierr,lambguess)
 else
    call reconstruct_maxent(mu,x,f,lambsol,lsum,ierr)
 endif

end subroutine reconstruct

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_maxent(mu,x,f,lambsol,lsum,ierr,lambguess,use_log,weights)
 use integrate, only:check_log_grid_is_uniform
 real, intent(in) :: mu(:)
 real(dp), intent(in) :: x(:)
 real, intent(out) :: f(size(x))
 real(dp), intent(out) :: lambsol(size(mu))
 real(dp), intent(out) :: lsum(size(mu))
 integer, intent(out) :: ierr
 real(dp), intent(in), optional :: lambguess(size(mu))
 real(dp), intent(in), optional :: weights(size(x))
 logical, intent(in), optional :: use_log
 integer :: n_moments
 real, parameter :: tol = 1.e-10
 real(dp), allocatable :: logx(:)
 logical :: use_gauss_legendre = .true.


 lambsol = 0.
 ierr = 0
 f = 0.
 n_moments = size(mu)

 ! initial guesses for Lagrange multipliers
 if (present(lambguess)) then
    ! use guesses passed as arguments (e.g. from previous attempt)
    lambsol = lambguess
    !print*,' using lambguess = ',lambguess
 else
    ! use predefined guess
    lambsol(1) = -log(sqrt(2.*pi))
 endif

 !
 ! if we are using a log grid, precalculate the grid
 ! of logx (this is an optimisation)
 !
 use_log_grid = .false.
 if (present(use_log)) use_log_grid = use_log
 if (use_log_grid) then
    log_grid_is_uniform = check_log_grid_is_uniform(x)
    allocate(logx(size(x)))
    logx = log(x)
 endif
 !
 ! if pre-calculated weights are sent, use Gauss-Legendre quadrature
 !
 use_gauss_legendre = present(weights)

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
 use integrate, only:integrate_trap,integrate_trap_log,integrate_trap_log_uniform,integrate_gauss_legendre
 integer, intent(in) :: n
 real (dp), intent(in)  :: lamb(n) ! guess for  solution
 real (dp), intent(out) :: l_sum(n) ! function evaluated for given lamb
 integer :: k

 if (use_gauss_legendre) then
    do k=1,n
       l_sum(k) = integrate_gauss_legendre(size(x),weights,integrand(x,lamb,k-1,logx=logx)) - mu(k)
    enddo
 else
    if (use_log_grid) then
       if (log_grid_is_uniform) then
          do k=1,n
             l_sum(k) = integrate_trap_log_uniform(size(x),x,integrand(x,lamb,k-1,logx=logx)) - mu(k)
          enddo
       else
          do k=1,n
             l_sum(k) = integrate_trap_log(size(x),x,integrand(x,lamb,k-1,logx=logx)) - mu(k)
          enddo
       endif
    else
       do k=1,n
          l_sum(k) = integrate_trap(size(x),x,integrand(x,lamb,k-1)) - mu(k)
       enddo
    endif
 endif
 !print*,' moments are: ',l_sum(:) + mu(:)
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
function integrand(x,lamb,k,logx) result(y)
 real(dp), intent(in) :: x(:)
 real(dp), intent(in) :: lamb(:)
 real(dp) :: y(size(x))  ! result
 real(dp), intent(in), optional :: logx(:)
 integer, intent(in) :: k
 real(dp) :: xi(size(lamb))
 integer :: i,j
 !real, parameter :: e100 = exp(-100.)

 do i=1,size(x)
    if (use_log_grid) then
       if (present(logx)) then
          do j=1,size(lamb)
             xi(j) = logx(i)**(j-1)
          enddo
       else
          do j=1,size(lamb)
             xi(j) = log(x(i))**(j-1)
          enddo
       endif
    else
       do j=1,size(lamb)
          xi(j) = x(i)**(j-1)
       enddo
    endif
    ! function is exp( lam_0 + lam_1 x + lam_2 x^2 + ...)
    ! or if log_grid is set exp( lam_0 + lam_1 log(x) + lam_2 log(x)^2 + ...)
    !term = dot_product(lamb,xi)
    !if (term < -100.) then ! avoid wasted effort on exp(-100) etc
    !   y(i) = x(i)**(k/3.)*e100
    !else
    !   y(i) = x(i)**(k/3.)*exp(term)
    !endif

    ! GENERALIZED GAMMA DISTRIBUTION (standard Gamma distribution if p=1)
    y(i) = x(i)**(k/3.)*gamma_func(x(i),lamb)

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
 case(5)
     string = 'gave up on k_3'
 case default
    string = ''
 end select

end function fsolve_error

!
! maximum entropy reconstruction of function given moments
!
subroutine reconstruct_gamma_dist(mu,lambsol,lsum,ierr,lambguess,verbose)
 real(dp), intent(in) :: mu(:)
 real(dp), intent(out) :: lambsol(size(mu))
 real(dp), intent(out) :: lsum(size(mu))
 integer, intent(out) :: ierr
 real, intent(in), optional :: lambguess(size(mu))
 logical, intent(in), optional :: verbose
 integer :: n_moments,k
 real, parameter :: tol = 1.e-2
 logical :: debug

 lambsol = 0.
 ierr = 0
 n_moments = 2 !size(mu)
 if (mu(1) < tiny(0.)) then
    lambsol = 0.
    lsum = 0.
    ierr = 1
    return
 endif
 debug = .false.
 if (present(verbose)) debug = verbose

 ! initial guesses for Lagrange multipliers
 if (present(lambguess)) then
    lambsol = lambguess
 else
    lambsol(1) = 2.
    if (n_moments > 1) lambsol(2) = 0.5
 endif

 if (debug) print*,'INPUT  moments = ',mu,'guess = ',lambsol(1:2)
 call fsolve(residual_fit_gamma,n_moments,lambsol,lsum,tol,ierr)
 lambsol = abs(lambsol)

 if ((ierr /= 1 .or. any(abs(lsum(1:n_moments)) > 0.1))) then
    if (debug) then
       print*,'err=',ierr,'2 parameter moments = ',(gamma_func_moment(n_moments,lambsol,mu,k),k=0,3),&
              'd_on_p,p=',lambsol(1:2),'err=',lsum(1:2)
    endif

    ! try two parameter solve with a different initial guess
    lambsol(1) = 1.1
    if (n_moments > 1) lambsol(2) = 2.
    call fsolve(residual_fit_gamma,n_moments,lambsol,lsum,tol,ierr)
    lambsol = abs(lambsol)

    if (debug) then
       print*,'err=',ierr,'2 parameter moments = ',(gamma_func_moment(n_moments,lambsol,mu,k),k=0,3),&
           'd_on_p,p=',lambsol(1:2),'err=',lsum(1:2)
    endif

    !
    ! if the error is > 15% or the parameters are negative, or d/p > 30, then
    ! revert to a one parameter solve with p fixed to 1
    !
    if (any(abs(lsum(1:n_moments)) > 0.15)) then! .or. any(lambsol(1:n_moments) < 0) .or. lambsol(1)*lambsol(2) > 30.) then
       lambsol(1) = 1.5
       if (n_moments > 1) lambsol(2) = 1.0
       call fsolve(residual_fit_gamma,1,lambsol,lsum,tol,ierr)

       ! since we take the abs, must be able to get the same solution with d_on_p > 0
       if (lambsol(1) < 0.) then
          lambsol(1) = abs(lambsol(1))
          call fsolve(residual_fit_gamma,1,lambsol,lsum,tol,ierr)
       endif

       ierr = 5  ! report that we gave up on k_3
       lsum(2) = gamma_func_moment(n_moments,lambsol,mu,3)/mu(4) - 1.

       if (debug) then
          print*,'err=',ierr,'1 parameter moments = ',(gamma_func_moment(n_moments,lambsol,mu,k),k=0,3),&
              'd_on_p,p=',lambsol(1:2),'err=',(gamma_func_moment(n_moments,lambsol,mu,k+1)/mu(k+2) - 1.,k=1,2)
       endif
    endif
 endif

contains

subroutine residual_fit_gamma(n,lamb,l_sum)
 integer, intent(in) :: n
 real(dp), intent(in)  :: lamb(n) ! guess for  solution
 real(dp), intent(out) :: l_sum(n) ! function evaluated for given lamb
 integer :: k

 ! solve for moments 3 and 4 (k=2,3)
 do k=1,n
    ! lsum is the error between the desired moments and the moments of the distribution
    l_sum(k) = gamma_func_moment(n,lamb,mu,k+1)/mu(k+2) - 1.
 enddo

end subroutine residual_fit_gamma

end subroutine reconstruct_gamma_dist

!--------------------------------------------------------------------
! GENERALIZED GAMMA DISTRIBUTION
!
!  f(x) = beta * p / theta * (x/theta)^(d-1) * exp(-(x/theta)^p) / Gamma(d/p)
!
! input parameters are beta, theta, d_on_p and p
! for p=1 this is just the Gamma distribution
!
! see: https://en.wikipedia.org/wiki/Gamma_distribution
!--------------------------------------------------------------------
real function gamma_func(x,params)
 real(dp), intent(in) :: x
 real(dp), intent(in) :: params(4)
 real(dp) :: beta,theta,d_on_p,p,d,expterm

 beta = params(1)   ! overall normalisation factor
 if (beta < tiny(0.)) then
    gamma_func = 0.
    return
 endif
 theta = params(2)  ! scale parameter theta on wikipedia
 d_on_p = abs(params(3))     ! shape parameter k on wikipedia
 p = abs(params(4))     ! shape parameter
 d = d_on_p * p

 ! use expression below to avoid NaNs with large numbers
 expterm = exp(-(x/theta)**p)
 if (expterm < tiny(0.)) then
    gamma_func = 0.
 else
    !expterm = (theta)**(-d) * (expterm)
    !if (expterm < tiny(0.)) then
    !   gamma_func = 0.
    !else
    !   gamma_func = beta * p / Gamma(d_on_p) * x**(d-1) * (expterm)
    !endif
    !gamma_func = beta * p / theta / Gamma(d_on_p) * (x/theta)**(d-1) * (expterm)
    gamma_func = beta * p / Gamma(d_on_p) * x**(d-1) * ((theta)**(-d) * (expterm))
 endif
 !if (isnan(gamma_func)) print*,gamma_func,beta,p, Gamma(d_on_p), x**(d-1), ((theta)**(-d) * (expterm))

end function gamma_func

!----------------------------------------------------------------------
! Evaluate the generalised gamma function but with the parameters
! given in terms of the first two moments k_0 and k_1 (mu(1) and mu(2))
! and the fitted parameters d_on_p and p (params(1) and params(2))
!----------------------------------------------------------------------
real function gamma_func_from_moments(x,mu,params)
 real(dp), intent(in) :: x
 real(dp), intent(in) :: mu(2),params(2)
 real(dp) :: beta,theta,d_on_p,p,all_params(4)

 d_on_p = abs(params(1))
 p = abs(params(2))
 theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./(3.*p)))**3
 beta = mu(1)   ! overall normalisation factor

 all_params = [beta,theta,d_on_p,p]
 gamma_func_from_moments = gamma_func(x,all_params)

end function gamma_func_from_moments

!----------------------------------------------------
! analytic moments of generalised gamma distribution
! given desired moments mu(1) and mu(2)
! and parameters d_on_p and p (lambsol)
!
! can either give d_on_p and p (n=2),
! or d alone (n=1) in which case p=1
!----------------------------------------------------
real function gamma_func_moment(n,lambsol,mu,k)
 integer, intent(in) :: n,k
 real(dp), intent(in) :: lambsol(n)
 real(dp), intent(in) :: mu(2)
 real(dp) :: d_on_p,theta,p

 d_on_p = abs(lambsol(1))
 if (n > 1) then
    p = abs(lambsol(2))
 else
    p = 1.
 endif
 theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./(3.*p)))**3

 gamma_func_moment = mu(1)*theta**(k/3.)*Gamma(d_on_p + k/(3.*p))/Gamma(d_on_p)

end function gamma_func_moment

end module reconstruct_from_moments
