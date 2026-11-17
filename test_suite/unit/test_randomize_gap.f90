program test_randomize_gap
  ! Unit test: options -random_gap and -random_outside_gap (mess_up_SPH): the particles
  ! inside (resp. outside) the gap of a planet are rotated around the z axis, the others
  ! do not move. Star at the origin, planet of 1e-3 Msun at 50 au (Hill radius 3.5 au),
  ! gap of +/- 2 Hill radii; one particle in the gap, one outside.
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use mcfost_env, only : dp
  use mess_up_SPH, only : randomize_gap

  implicit none

  integer :: n_fail
  real(dp), parameter :: udist = 1.495978707e13_dp, factor = 2.0_dp

  n_fail = 0
  call run(.true.,  moved_in=.true.,  moved_out=.false., label="-random_gap")
  call run(.false., moved_in=.false., moved_out=.true.,  label="-random_outside_gap")

  if (n_fail > 0) then
     write(*,*) n_fail, " failure(s)"
     call exit(1)
  endif
  write(*,*) "test_randomize_gap: all checks passed"

contains

  subroutine run(inside, moved_in, moved_out, label)
    logical, intent(in) :: inside, moved_in, moved_out
    character(len=*), intent(in) :: label

    real(dp) :: xyzh(4,2), vxyzu(4,2), xyzmh_ptmass(18,2), xyzh0(4,2)
    logical :: moved(2)

    xyzmh_ptmass = 0.0_dp
    xyzmh_ptmass(4,1) = 1.0_dp                                ! star
    xyzmh_ptmass(1,2) = 50.0_dp ; xyzmh_ptmass(4,2) = 1.0e-3_dp ! planet
    xyzh(:,1) = [0.0_dp, 51.0_dp, 0.0_dp, 1.0_dp] ! in the gap (|r - 50| < 7 au)
    xyzh(:,2) = [0.0_dp, 20.0_dp, 0.0_dp, 1.0_dp] ! outside the gap
    vxyzu = 0.0_dp ; vxyzu(1,:) = 1.0_dp
    xyzh0 = xyzh

    call randomize_gap(2, 2, xyzh, vxyzu, xyzmh_ptmass, udist, factor, inside)

    moved(:) = any(abs(xyzh(1:2,:) - xyzh0(1:2,:)) > 1e-10_dp, dim=1)
    call check(moved(1) .eqv. moved_in, label//": particle in the gap")
    call check(moved(2) .eqv. moved_out, label//": particle outside the gap")
    ! a rotation around z keeps the cylindrical radius and the height
    call check(all(abs(hypot(xyzh(1,:),xyzh(2,:)) - hypot(xyzh0(1,:),xyzh0(2,:))) < 1e-10_dp) .and. &
         all(xyzh(3:4,:) == xyzh0(3:4,:)), label//": rotation around the z axis")

  end subroutine run

  subroutine check(ok, label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if (.not.ok) then
       write(*,*) "FAIL: ", label
       n_fail = n_fail + 1
    endif
  end subroutine check

end program test_randomize_gap
