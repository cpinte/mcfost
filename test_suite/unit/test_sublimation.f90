program test_sublimation
  ! Unit test: sublimation decision used in the physics loop (disk_physics)
  !  - sublimation_band : hysteresis half-width from the Monte Carlo noise
  !  - sublimation_change : sublimation, re-condensation, frozen grains
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use disk_physics, only : sublimation_band, sublimation_change, sublimation_precision, &
       sublimation_band_max, sublimation_n_sigma, sublimation_max_flips, &
       sublimation_no_change, sublimation_sublimates, sublimation_recondenses, sublimation_frozen

  implicit none

  integer :: n_fail, n_test
  real, parameter :: T_sub = 1500.0, delta = 0.02

  n_fail = 0 ; n_test = 0

  ! hysteresis band
  call check_r("band, no packet",           sublimation_band(0.0), sublimation_band_max)
  call check_r("band, below min packets",   sublimation_band(5.0), sublimation_band_max)
  call check_r("band, 100 packets",         sublimation_band(100.0), sublimation_n_sigma * 0.25 / 10.0)
  call check_r("band, many packets (floor)", sublimation_band(1.0e8), sublimation_precision)
  call check_l("band monotonic", sublimation_band(1.0e3) <= sublimation_band(1.0e2), .true.)

  ! dust present
  call check_i("present, T well above T_sub", sublimation_change(.false., 1600.0, T_sub, delta, 0), &
       sublimation_sublimates)
  call check_i("present, T just above T_sub (in band)", sublimation_change(.false., 1520.0, T_sub, delta, 0), &
       sublimation_no_change)
  call check_i("present, T below T_sub", sublimation_change(.false., 1000.0, T_sub, delta, 0), &
       sublimation_no_change)

  ! dust sublimated (test-particle temperature)
  call check_i("sublimated, T well below T_sub", sublimation_change(.true., 1400.0, T_sub, delta, 0), &
       sublimation_recondenses)
  call check_i("sublimated, T just below T_sub (in band)", sublimation_change(.true., 1480.0, T_sub, delta, 0), &
       sublimation_no_change)
  call check_i("sublimated, T above T_sub", sublimation_change(.true., 2000.0, T_sub, delta, 0), &
       sublimation_no_change)

  ! frozen after too many changes, only when a change would happen
  call check_i("frozen, would sublimate", sublimation_change(.false., 1600.0, T_sub, delta, sublimation_max_flips), &
       sublimation_frozen)
  call check_i("frozen, would re-condense", sublimation_change(.true., 1400.0, T_sub, delta, sublimation_max_flips), &
       sublimation_frozen)
  call check_i("not frozen below the limit", sublimation_change(.true., 1400.0, T_sub, delta, sublimation_max_flips-1), &
       sublimation_recondenses)
  call check_i("frozen but no change", sublimation_change(.false., 1000.0, T_sub, delta, sublimation_max_flips), &
       sublimation_no_change)

  ! snow line : low T_sub
  call check_i("ice, T = 150 K > T_sub = 100 K", sublimation_change(.false., 150.0, 100.0, delta, 0), &
       sublimation_sublimates)

  write(*,'(a,i0,a,i0)') "test_sublimation: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  subroutine check_i(name, got, expected)
    character(len=*), intent(in) :: name
    integer, intent(in) :: got, expected

    n_test = n_test + 1
    if (got /= expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,i0,a,i0)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check_i

  subroutine check_r(name, got, expected)
    character(len=*), intent(in) :: name
    real, intent(in) :: got, expected

    n_test = n_test + 1
    if (abs(got - expected) > 1.0e-6 * abs(expected)) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,es12.5,a,es12.5)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check_r

  subroutine check_l(name, got, expected)
    character(len=*), intent(in) :: name
    logical, intent(in) :: got, expected

    n_test = n_test + 1
    if (got .neqv. expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,l1,a,l1)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check_l

end program test_sublimation
