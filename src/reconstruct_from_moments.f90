module reconstruct_from_moments

  use constantes

  implicit none

  public :: reconstruct_gamma_dist,gamma_func,gamma_func_from_moments,gamma_func_moment,fsolve_error

  private

contains

  function fsolve_error(ierr) result(string)
    ! print the error message corresponding to the error code
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

  !*********************************************************

  subroutine reconstruct_gamma_dist(mu,lambsol,lsum,ierr,lambguess,verbose)
    ! reconstruction of function given moments
    real(dp), intent(in) :: mu(:)
    real(dp), intent(out) :: lambsol(size(mu))
    real(dp), intent(out) :: lsum(size(mu))
    integer, intent(out) :: ierr
    real(dp), intent(in), optional :: lambguess(size(mu))
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

       ! if the error is > 15% or the parameters are negative, or d/p > 30, then
       ! revert to a one parameter solve with p fixed to 1
       if (any(abs(lsum(1:n_moments)) > 0.15) .or. abs(lambsol(1)*lambsol(2)) > 30.) then
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

  !*********************************************************

  real(dp) function gamma_func(x,params)
    ! GENERALIZED GAMMA DISTRIBUTION
    !
    !  f(x) = beta * p / theta * (x/theta)^(d-1) * exp(-(x/theta)^p) / Gamma(d/p)
    !
    ! input parameters are beta, theta, d_on_p and p
    ! for p=1 this is just the Gamma distribution
    !
    ! see: https://en.wikipedia.org/wiki/Gamma_distribution

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
    if (expterm < 1e6*tiny_real) then
       gamma_func = 0.
    else
       !expterm = (theta)**(-d) * (expterm)
       !if (expterm < tiny(0.)) then
       !   gamma_func = 0.
       !else
       !   gamma_func = beta * p / Gamma(d_on_p) * x**(d-1) * (expterm)
       !endif
       !gamma_func = beta * p / theta / Gamma(d_on_p) * (x/theta)**(d-1) * (expterm)
       !write(*,*) "test", beta, p, d_on_p, d, expterm
       !write(*,*) beta * p
       !write(*,*)  beta * p / Gamma(d_on_p)
       !write(*,*)  x**(d-1)
       !write(*,*)  ((theta)**(-d) * (expterm))
       !write(*,*) "TEST", beta * p / Gamma(d_on_p) * x**(d-1) * ((theta)**(-d) * (expterm))

       gamma_func = beta * p / Gamma(d_on_p) * x**(d-1) * ((theta)**(-d) * (expterm))
    endif
    !if (isnan(gamma_func)) write(*,*) gamma_func,beta,p, Gamma(d_on_p), x**(d-1), ((theta)**(-d) * (expterm))

  end function gamma_func

  !*********************************************************

  real(dp) function gamma_func_from_moments(x,mu,params)
    ! Evaluate the generalised gamma function but with the parameters
    ! given in terms of the first two moments k_0 and k_1 (mu(1) and mu(2))
    ! and the fitted parameters d_on_p and p (params(1) and params(2))

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

  !*********************************************************

  real(dp) function gamma_func_moment(n,lambsol,mu,k)
    ! analytic moments of generalised gamma distribution
    ! given desired moments mu(1) and mu(2)
    ! and parameters d_on_p and p (lambsol)
    !
    ! can either give d_on_p and p (n=2),
    ! or d alone (n=1) in which case p=1

    integer, intent(in) :: n,k
    real(dp), intent(in) :: lambsol(n)
    real(dp), intent(in) :: mu(2)
    real(dp) :: d_on_p,theta,p, ratio

    d_on_p = abs(lambsol(1))
    if (n > 1) then
       p = abs(lambsol(2))
    else
       p = 1.
    endif
    theta = (mu(2)/mu(1) * Gamma(d_on_p) / Gamma(d_on_p + 1./(3.*p)))**3

    !write(*,*) "test", d_on_p + k/(3.*p), d_on_p
    !write(*,*) Gamma(d_on_p)
    !write(*,*) "BBB", mu(1)*theta**(k/3.), Gamma(d_on_p + k/(3.*p)), Gamma(d_on_p)

    ratio  = Gamma(d_on_p + k/(3.*p))/Gamma(d_on_p)

    gamma_func_moment = mu(1)*theta**(k/3.)*ratio

    !if (isnan(gamma_func_moment)) write(*,*)  "test", mu(1)*theta**(k/3.), Gamma(d_on_p + k/(3.*p)), Gamma(d_on_p)

  end function gamma_func_moment

end module reconstruct_from_moments
