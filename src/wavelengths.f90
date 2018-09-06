module wavelengths

  use mcfost_env, only : dp, band
  use parametres, only : n_pop
  use messages

  implicit none

  logical :: lmono0, lmono

  ! Gamme de longueurs d'onde utilisees
  real :: lambda_min, lambda_max
  integer :: n_lambda, n_lambda2
  real(kind=dp) :: delta_lambda
  real(kind=dp), dimension(:), allocatable :: tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda !n_lambda

  real, dimension (:,:), allocatable :: tab_amu1, tab_amu2 !n_lambda,n_pop
  real, dimension (:,:), allocatable :: tab_amu1_coating, tab_amu2_coating !n_lambda,n_pop
  real, dimension(:), allocatable :: tab_lambda2, tab_lambda2_inf, tab_lambda2_sup, tab_delta_lambda2

contains


subroutine init_lambda()
  ! Initialisation table de longueurs d'onde

  integer :: i, alloc_status

  allocate(tab_lambda(n_lambda), tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_lambda')
  tab_lambda=0.0 ;  tab_lambda_inf=0.0 ; tab_lambda_sup=0.0 ; tab_delta_lambda=0.0

  allocate(tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_amu')
  tab_amu1=0.0; tab_amu2=0.0;

  allocate(tab_amu1_coating(n_lambda, n_pop), tab_amu2_coating(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_amu_coating')
  tab_amu1_coating=0.0; tab_amu2_coating=0.0;

  if (lmono0) then
     ! Lecture longueur d'onde
     read(band,*) tab_lambda(1)
     tab_delta_lambda(1) = 1.0
     tab_lambda_inf(1) = tab_lambda(1)
     tab_lambda_sup(1) = tab_lambda(1)
  else
     ! Initialisation longueurs d'onde
     !delta_lambda = (lambda_max/lambda_min)**(1.0/real(n_lambda))
     delta_lambda =  exp( (1.0_dp/real(n_lambda,kind=dp)) * log(lambda_max/lambda_min) )

     tab_lambda_inf(1) = lambda_min
     tab_lambda(1)=lambda_min*sqrt(delta_lambda)
     tab_lambda_sup(1) = lambda_min*delta_lambda
     do i=2, n_lambda
        tab_lambda(i)= tab_lambda(i-1)*delta_lambda
        tab_lambda_sup(i)= tab_lambda_sup(i-1)*delta_lambda
        tab_lambda_inf(i)= tab_lambda_sup(i-1)
     enddo

     do i=1, n_lambda
        tab_delta_lambda(i) = tab_lambda_sup(i) - tab_lambda_inf(i)
     enddo

  endif

end subroutine init_lambda

!**********************************************************************

subroutine init_lambda2()
  ! Initialisation table en lambda sed

  implicit none

  integer :: i

  n_lambda=n_lambda2
  do i=1, n_lambda2
     tab_lambda(i)= tab_lambda2(i)
     tab_lambda_inf(i)= tab_lambda2_inf(i)
     tab_lambda_sup(i)= tab_lambda2_sup(i)
     tab_delta_lambda(i)= tab_delta_lambda2(i)
  enddo

  return

end subroutine init_lambda2

!**********************************************************************


end module wavelengths
