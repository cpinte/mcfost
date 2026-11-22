module numa
  ! Memory placement on NUMA nodes (multi-socket nodes, or sockets split in several domains).
  ! A memory page is placed on the NUMA domain of the thread that first writes it ("first
  ! touch"). The routines below fill an array in parallel, in equal contiguous blocks per
  ! thread (static schedule), so that:
  !  - arrays indexed by cell and read by every thread (opacities, emission CDFs, Voronoi data)
  !    have their pages spread over the domains, as with numactl --interleave;
  !  - per-thread arrays with the thread index last (xKJ_abs(:,id), xT_ech(:,id)) have the block
  !    of thread id-1 equal to its own slice, which stays local to it.
  ! Same values as a serial fill (bit-for-bit). Called with the array and its size (sequence
  ! association: any rank), e.g. call first_touch_dp(kappa, size(kappa,kind=int64), 0.0_dp).
  ! Measured on AMD EPYC 7543 (8 NUMA domains, 64 threads): -32 % on the Monte Carlo of a
  ! Phantom/Voronoi temperature run, -33 % on the Monte Carlo of its images (review tracker R4).

  use mcfost_env, only : dp
  use iso_fortran_env, only : int64

  implicit none

  public :: first_touch_dp, first_touch_sp, first_touch_int, int64
  private

contains

  subroutine first_touch_dp(a, n, val)
    integer(kind=int64), intent(in) :: n
    real(kind=dp), dimension(n), intent(out) :: a
    real(kind=dp), intent(in) :: val
    integer(kind=int64) :: i
    !$omp parallel do schedule(static) default(none) shared(a,n,val) private(i)
    do i=1, n
       a(i) = val
    enddo
    !$omp end parallel do
  end subroutine first_touch_dp

  subroutine first_touch_sp(a, n, val)
    integer(kind=int64), intent(in) :: n
    real, dimension(n), intent(out) :: a
    real, intent(in) :: val
    integer(kind=int64) :: i
    !$omp parallel do schedule(static) default(none) shared(a,n,val) private(i)
    do i=1, n
       a(i) = val
    enddo
    !$omp end parallel do
  end subroutine first_touch_sp

  subroutine first_touch_int(a, n, val)
    integer(kind=int64), intent(in) :: n
    integer, dimension(n), intent(out) :: a
    integer, intent(in) :: val
    integer(kind=int64) :: i
    !$omp parallel do schedule(static) default(none) shared(a,n,val) private(i)
    do i=1, n
       a(i) = val
    enddo
    !$omp end parallel do
  end subroutine first_touch_int

end module numa
