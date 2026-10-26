program test_version
  ! Unit test: version comparison used by the automatic update (utils)
  !  - is_newer_version : a binary is only updated to a newer release, never
  !    downgraded (e.g. a development version not released yet)
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use utils, only : is_newer_version

  implicit none

  integer :: n_fail, n_test

  n_fail = 0 ; n_test = 0

  call check("newer patch",               is_newer_version("4.1.16", "4.1.15"), .true.)
  call check("older patch (no downgrade)", is_newer_version("4.1.14", "4.1.15"), .false.)
  call check("same version",              is_newer_version("4.1.15", "4.1.15"), .false.)
  call check("numeric, not lexical",      is_newer_version("4.1.10", "4.1.9"),  .true.)
  call check("newer minor",               is_newer_version("4.2.0", "4.1.15"),  .true.)
  call check("newer major",               is_newer_version("5.0.0", "4.9.99"),  .true.)
  call check("missing patch number",      is_newer_version("4.2", "4.1.15"),    .true.)
  call check("missing patch, same",       is_newer_version("4.1", "4.1.0"),     .false.)

  write(*,'(a,i0,a,i0)') "test_version: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  subroutine check(name, got, expected)
    character(len=*), intent(in) :: name
    logical, intent(in) :: got, expected

    n_test = n_test + 1
    if (got .neqv. expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,l1,a,l1)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check

end program test_version
