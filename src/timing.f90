module timing
  ! Wall-clock time per phase of a run, accumulated over all the calls of the
  ! phase (e.g. every hydrostatic iteration or wavelength). Timers are started
  ! and stopped outside parallel regions. Printed at the end of the run with
  ! -timing, one "TIMING <phase> <seconds> <calls>" line per phase, read by the
  ! benchmark scripts (performance_suite). Phases can be nested: the times of
  ! the inner phases are included in the outer ones.

  use, intrinsic :: iso_fortran_env, only : int64

  implicit none

  integer, parameter, public :: t_total = 1, t_setup = 2, t_init = 3, t_dust_prop = 4, t_opacity = 5, &
       t_mc_temperature = 6, t_temperature = 7, t_temperature_nRE = 8, t_step2_setup = 9, &
       t_mc_step2 = 10, t_source_fct = 11, t_ray_tracing = 12, t_output = 13, t_mol = 14, t_atom = 15
  integer, parameter :: n_timers = 15

  character(len=16), dimension(n_timers), parameter :: timer_name = [character(len=16) :: &
       "total", "setup", "init", "dust_prop", "opacity", "mc_temperature", "temperature", &
       "temperature_nRE", "step2_setup", "mc_step2", "source_fct", "ray_tracing", "output", "mol", "atom"]

  logical, public :: ltiming = .false.

  integer(int64), dimension(n_timers) :: count_start = 0, count_sum = 0
  integer, dimension(n_timers) :: n_calls = 0

  public :: timer_start, timer_stop, write_timing

  private

contains

  subroutine timer_start(i)

    integer, intent(in) :: i

    call system_clock(count_start(i))

  end subroutine timer_start

  !---------------------------------

  subroutine timer_stop(i)

    integer, intent(in) :: i
    integer(int64) :: count

    call system_clock(count)
    count_sum(i) = count_sum(i) + (count - count_start(i))
    n_calls(i) = n_calls(i) + 1

  end subroutine timer_stop

  !---------------------------------

  subroutine write_timing()

    integer(int64) :: count_rate
    integer :: i

    if (.not.ltiming) return

    call system_clock(count_rate=count_rate)
    do i=1, n_timers
       if (n_calls(i) > 0) write(*,'(" TIMING ", A16, F12.3, I8)') timer_name(i), &
            real(count_sum(i),kind=8)/real(count_rate,kind=8), n_calls(i)
    enddo

  end subroutine write_timing

end module timing
