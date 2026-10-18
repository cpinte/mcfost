program test_nan_checks
  ! Unit test: NaN and infinity detection in utils.f90 must work with the
  ! release compiler flags (-ffinite-math-only, -fp-model fast=2), which
  ! remove comparisons such as y /= y. Built by "make unit_tests" in src/
  ! and run by test_suite/test_unit.py. Exits with status 1 on failure.

  use mcfost_env, only : dp
  use utils, only : is_nan_infinity, is_nan_infinity_vector, is_nan_infinity_matrix, &
       is_nan_dp, is_inf_dp, is_finite_dp
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf

  implicit none

  real(dp) :: nan, pinf, minf, x(3), m(2,2)
  integer :: n_fail, n_test

  nan  = ieee_value(1.0_dp, ieee_quiet_nan)
  pinf = ieee_value(1.0_dp, ieee_positive_inf)
  minf = ieee_value(1.0_dp, ieee_negative_inf)
  n_fail = 0 ; n_test = 0

  ! is_nan_infinity : 0 = finite, 1 = NaN, 2 = +/- infinity
  call check("is_nan_infinity(1.0)",    is_nan_infinity(1.0_dp), 0)
  call check("is_nan_infinity(-0.0)",   is_nan_infinity(-0.0_dp), 0)
  call check("is_nan_infinity(huge)",   is_nan_infinity(huge(1.0_dp)), 0)
  call check("is_nan_infinity(-huge)",  is_nan_infinity(-huge(1.0_dp)), 0)
  call check("is_nan_infinity(subnormal)", is_nan_infinity(tiny(1.0_dp)/4), 0)
  call check("is_nan_infinity(NaN)",    is_nan_infinity(nan), 1)
  call check("is_nan_infinity(+Inf)",   is_nan_infinity(pinf), 2)
  call check("is_nan_infinity(-Inf)",   is_nan_infinity(minf), 2)

  x = [1.0_dp, 2.0_dp, 3.0_dp] ; call check("vector finite", is_nan_infinity_vector(x), 0)
  x = [1.0_dp, nan, 2.0_dp]    ; call check("vector NaN",    is_nan_infinity_vector(x), 1)
  x = [1.0_dp, 2.0_dp, pinf]   ; call check("vector +Inf",   is_nan_infinity_vector(x), 2)

  m = 1.0_dp                   ; call check("matrix finite", is_nan_infinity_matrix(m), 0)
  m = 1.0_dp ; m(2,1) = nan    ; call check("matrix NaN",    is_nan_infinity_matrix(m), 1)
  m = 1.0_dp ; m(1,2) = minf   ; call check("matrix -Inf",   is_nan_infinity_matrix(m), 2)

  ! elemental helpers
  call check_l("is_nan_dp(NaN)",       is_nan_dp(nan), .true.)
  call check_l("is_nan_dp(+Inf)",      is_nan_dp(pinf), .false.)
  call check_l("is_inf_dp(-Inf)",      is_inf_dp(minf), .true.)
  call check_l("is_inf_dp(NaN)",       is_inf_dp(nan), .false.)
  call check_l("is_finite_dp(1e300)",  is_finite_dp(1.0e300_dp), .true.)
  call check_l("is_finite_dp(NaN)",    is_finite_dp(nan), .false.)
  call check_l("is_finite_dp(+Inf)",   is_finite_dp(pinf), .false.)

  write(*,'(a,i0,a,i0)') "test_nan_checks: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  subroutine check(name, got, expected)
    character(len=*), intent(in) :: name
    integer, intent(in) :: got, expected

    n_test = n_test + 1
    if (got /= expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,i0,a,i0)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check

  subroutine check_l(name, got, expected)
    character(len=*), intent(in) :: name
    logical, intent(in) :: got, expected

    n_test = n_test + 1
    if (got .neqv. expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,l1,a,l1)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check_l

end program test_nan_checks
