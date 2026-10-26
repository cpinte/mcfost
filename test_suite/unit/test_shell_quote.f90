program test_shell_quote
  ! Unit test: shell commands built from user strings (utils, H5)
  !  - shell_quote : a path with spaces, $(...), backticks, double and single
  !    quotes or ; is passed to the shell as one word and never executed
  !  - system_call : exit status of the command
  ! Works in a scratch directory created in the current directory and removed at the end.
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use utils, only : shell_quote, system_call, is_dir, is_file

  implicit none

  character(len=*), parameter :: base = "test_shell_quote.tmp"
  character(len=*), parameter :: weird = "My Models $(touch PWNED1) `touch PWNED2` ""q"" it's; touch PWNED3"

  character(len=:), allocatable :: dir
  character(len=256) :: line
  integer :: n_fail, n_test, status, unit

  n_fail = 0 ; n_test = 0

  call system_call("rm -rf "//shell_quote(base), status)
  call system_call("mkdir -p "//shell_quote(base), status)
  dir = base//"/"//weird

  call check_int("exit status of true", status_of("true"), 0)
  call check_int("exit status of false", merge(1, 0, status_of("false") /= 0), 1)

  call check_str("no special character", shell_quote("data_th"), "'data_th'")
  call check_str("single quote", shell_quote("it's"), "'it'\''s'")

  ! mkdir creates exactly the requested directory, without running anything
  call system_call("mkdir -p "//shell_quote(dir), status)
  call check_int("mkdir status", status, 0)
  call check_log("directory created", is_dir(dir))
  call check_log("$(...) not executed", .not.is_file("PWNED1") .and. .not.is_file(base//"/PWNED1"))
  call check_log("backticks not executed", .not.is_file("PWNED2") .and. .not.is_file(base//"/PWNED2"))
  call check_log("; not executed", .not.is_file("PWNED3") .and. .not.is_file(base//"/PWNED3"))

  ! the string reaches the command unchanged
  call system_call("printf '%s\n' "//shell_quote(weird)//" > "//shell_quote(dir//"/out.txt"), status)
  call check_int("printf status", status, 0)
  line = ""
  open(newunit=unit, file=dir//"/out.txt", status="old", action="read", iostat=status)
  if (status == 0) then
     read(unit,'(a)', iostat=status) line
     close(unit)
  endif
  call check_str("string passed unchanged", trim(line), weird)

  call system_call("rm -rf "//shell_quote(base), status)
  call check_log("scratch directory removed", .not.is_dir(base))

  write(*,'(a,i0,a,i0)') "test_shell_quote: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  integer function status_of(cmd)
    character(len=*), intent(in) :: cmd
    call system_call(cmd, status_of)
  end function status_of

  subroutine check_log(name, ok)
    character(len=*), intent(in) :: name
    logical, intent(in) :: ok

    n_test = n_test + 1
    if (.not.ok) then
       n_fail = n_fail + 1
       write(*,'(a,a)') "FAIL ", name
    endif
  end subroutine check_log

  subroutine check_int(name, got, expected)
    character(len=*), intent(in) :: name
    integer, intent(in) :: got, expected

    n_test = n_test + 1
    if (got /= expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,i0,a,i0)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check_int

  subroutine check_str(name, got, expected)
    character(len=*), intent(in) :: name, got, expected

    n_test = n_test + 1
    if (got /= expected) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,a,a,a)') "FAIL ", name, ": got ", got, ", expected ", expected
    endif
  end subroutine check_str

end program test_shell_quote
