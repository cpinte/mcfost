program test_sublimation_law
  ! Unit test: density-dependent sublimation laws
  !  - sublimation_temperature (disk_physics) : constant and power law
  !  - read_sublimation_law_file (dust_prop) : parsing of the .sub files
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use grains, only : sublimation_law_constant, sublimation_law_power_law
  use disk_physics, only : sublimation_temperature, sublimation_rho_floor
  use dust_prop, only : read_sublimation_law_file

  implicit none

  integer :: n_fail, n_test, law
  real :: T0, rho0, gamma
  character(len=*), parameter :: power_law_file = "test_sublimation_law_power_law.sub", &
       constant_file = "test_sublimation_law_constant.sub"

  n_fail = 0 ; n_test = 0

  ! constant law : independent of the density, power-law parameters ignored
  call check_r("constant law", sublimation_temperature(sublimation_law_constant, 1500.0, 2000.0, 1.0, 0.0195, 1e-10), &
       1500.0)

  ! power law, Pollack et al. 1994 coefficients
  call check_r("power law at rho0", sublimation_temperature(sublimation_law_power_law, 1500.0, 2000.0, 1.0, 0.0195, 1.0), &
       2000.0)
  call check_r("power law at 1e-10 g/cm3", &
       sublimation_temperature(sublimation_law_power_law, 1500.0, 2000.0, 1.0, 0.0195, 1e-10), &
       2000.0 * 1e-10**0.0195)
  call check_r("power law with rho0 /= 1", &
       sublimation_temperature(sublimation_law_power_law, 1500.0, 2000.0, 1e-8, 0.0195, 1e-10), &
       2000.0 * 1e-2**0.0195)
  call check_l("power law increases with density", &
       sublimation_temperature(sublimation_law_power_law, 0.0, 2000.0, 1.0, 0.0195, 1e-12) < &
       sublimation_temperature(sublimation_law_power_law, 0.0, 2000.0, 1.0, 0.0195, 1e-8), .true.)
  ! empty cells : floored density, finite temperature
  call check_r("power law, zero density (floor)", &
       sublimation_temperature(sublimation_law_power_law, 1500.0, 2000.0, 1.0, 0.0195, 0.0), &
       2000.0 * sublimation_rho_floor**0.0195)

  ! .sub file : power law, with comments, blank lines, mixed case keys and values
  open(unit=10, file=power_law_file, status="replace")
  write(10,'(a)') "# Sublimation properties for a test material"
  write(10,'(a)') "format_version = 1"
  write(10,'(a)') ""
  write(10,'(a)') "Law = Power_Law      # constant | power_law"
  write(10,'(a)') "T0 = 2000.           # K"
  write(10,'(a)') "  rho0=1.0e-2"
  write(10,'(a)') "gamma = 0.0195"
  write(10,'(a)') "reference = Pollack et al. 1994; Isella & Natta 2005"
  close(10)
  call read_sublimation_law_file(power_law_file, law, T0, rho0, gamma)
  call check_i("file, power law", law, sublimation_law_power_law)
  call check_r("file, T0", T0, 2000.0)
  call check_r("file, rho0", rho0, 1.0e-2)
  call check_r("file, gamma", gamma, 0.0195)

  ! .sub file : constant law, no parameter needed
  open(unit=10, file=constant_file, status="replace")
  write(10,'(a)') "format_version = 1"
  write(10,'(a)') "law = constant"
  close(10)
  call read_sublimation_law_file(constant_file, law, T0, rho0, gamma)
  call check_i("file, constant law", law, sublimation_law_constant)

  call delete_file(power_law_file)
  call delete_file(constant_file)

  write(*,'(a,i0,a,i0)') "test_sublimation_law: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  subroutine delete_file(filename)
    character(len=*), intent(in) :: filename

    open(unit=10, file=filename, status="old")
    close(10, status="delete")
  end subroutine delete_file

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
    if (abs(got - expected) > 1.0e-5 * abs(expected)) then
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

end program test_sublimation_law
