module MRW
  use mcfost_env, only : dp
  use constants, only : pi

  implicit none

  integer, parameter :: n = 10000
  real(dp), dimension(n), save :: zeta, y_MRW

  real(dp) :: gamma_MRW = 2.0_dp
  real(dp), parameter :: cst_ct = 3.0_dp / (pi * pi)
  real(dp), parameter :: mrw_wall_safety = 1.0e-5_dp
  integer, parameter :: n_min_interactions_MRW = 3
  integer, parameter :: n_mrw_step_max = 100000

contains

  pure real(dp) function cumulative_zeta_from_y(y_val) result(z_val)
    ! Evaluate Min et al. (2009) Eq. (7):
    !   P(t) = 2 * Sum_{j=1}^{infinity} (-1)^(j+1) * y^(j^2)
    ! For y <= 0.5, the direct series converges in <= 8 terms without cancellation.
    ! For y > 0.5, use the Jacobi theta_4 Poisson summation identity:
    !   P(t) = 1 - 2 * sqrt(pi / (-ln y)) * Sum_{k=0}^{infinity} exp(-pi^2 * (k + 1/2)^2 / (-ln y))
    ! which converges in <= 3 positive terms with zero catastrophic cancellation as y -> 1.
    real(dp), intent(in) :: y_val

    integer :: j
    real(dp) :: term, tau, s, expo

    if (y_val <= 0.0_dp) then
       z_val = 0.0_dp
    else if (y_val >= 1.0_dp) then
       z_val = 1.0_dp
    else if (y_val <= 0.5_dp) then
       s = 0.0_dp
       j = 0
       do
          j = j + 1
          term = y_val**(j * j)
          if (term < 1.0e-18_dp) exit
          if (mod(j, 2) == 0) then
             s = s - term
          else
             s = s + term
          endif
       enddo
       z_val = 2.0_dp * s
    else
       tau = -log(y_val)
       s = 0.0_dp
       j = 0
       do
          expo = (pi * pi) * ((real(j, dp) + 0.5_dp)**2) / tau
          if (expo > 100.0_dp) exit
          term = exp(-expo)
          s = s + term
          if (term < 1.0e-18_dp * s) exit
          j = j + 1
       enddo
       z_val = max(0.0_dp, min(1.0_dp, 1.0_dp - 2.0_dp * sqrt(pi / tau) * s))
    endif

  end function cumulative_zeta_from_y

  !----------------------------------------------

  subroutine initialize_cumulative_zeta()
    ! Pre-compute the inverse table y_MRW(zeta) on a uniform zeta grid in [0, 1]
    ! by bisection on Min et al. (2009) Eq. (7), enabling O(1) runtime sampling.

    implicit none

    logical, save :: lzeta_initialized = .false.
    integer :: i, iter
    real(dp) :: y_lo, y_hi, y_mid

    if (lzeta_initialized) return

    zeta(1) = 0.0_dp
    y_MRW(1) = 0.0_dp

    do i = 2, n - 1
       zeta(i) = real(i - 1, dp) / real(n - 1, dp)
       y_lo = y_MRW(i - 1)
       y_hi = 1.0_dp
       do iter = 1, 60
          y_mid = 0.5_dp * (y_lo + y_hi)
          if (cumulative_zeta_from_y(y_mid) < zeta(i)) then
             y_lo = y_mid
          else
             y_hi = y_mid
          endif
       enddo
       y_MRW(i) = 0.5_dp * (y_lo + y_hi)
    enddo

    zeta(n) = 1.0_dp
    y_MRW(n) = 1.0_dp
    lzeta_initialized = .true.

    return

  end subroutine initialize_cumulative_zeta

  !----------------------------------------------

  pure real(dp) function interp_y_from_zeta(zeta_val) result(y_val)
    ! O(1) lookup of y in (0, 1] from uniform random variate zeta_val in [0, 1).
    ! Clamps zeta_val away from 0 so -log(y_val) is always finite.
    real(dp), intent(in) :: zeta_val

    real(dp) :: xi, pos, frac
    integer :: idx

    xi = min(max(zeta_val, 1.0e-30_dp), 1.0_dp - epsilon(1.0_dp))
    pos = xi * real(n - 1, dp)
    idx = min(int(pos) + 1, n - 1)
    frac = min(max(pos - real(idx - 1, dp), 0.0_dp), 1.0_dp)
    y_val = y_MRW(idx) * (1.0_dp - frac) + y_MRW(idx + 1) * frac

  end function interp_y_from_zeta

  !----------------------------------------------

  real(dp) function sample_zeta(id)

    use naleat, only : stream
#include "sprng_f.h"

    integer, intent(in) :: id

    real(dp) :: zeta_random

    call initialize_cumulative_zeta()
    zeta_random = real(sprng(stream(id)), kind=dp)
    sample_zeta = interp_y_from_zeta(zeta_random)

  end function sample_zeta

  !----------------------------------------------

  subroutine sample_lambertian_direction(id, un, vn, wn, u, v, w)
    ! Sample an outward Lambertian direction (cos(psi) = sqrt(xi_1), phi = 2*pi*xi_2)
    ! around the unit outward sphere normal (un, vn, wn) (Min et al. 2009, Section 2.2).
    use naleat, only : stream
#include "sprng_f.h"

    integer, intent(in) :: id
    real(dp), intent(in) :: un, vn, wn
    real(dp), intent(out) :: u, v, w

    real(dp) :: rand1, rand2, cpsi, spsi, phi, a, b, c, inv_c

    rand1 = real(sprng(stream(id)), kind=dp)
    rand2 = real(sprng(stream(id)), kind=dp)

    cpsi = sqrt(rand1)
    spsi = sqrt(max(0.0_dp, 1.0_dp - rand1))
    phi = 2.0_dp * pi * rand2
    a = spsi * cos(phi)
    b = spsi * sin(phi)

    if (abs(wn) <= 0.9_dp) then
       c = sqrt(1.0_dp - wn * wn)
       inv_c = 1.0_dp / c
       u = (a * wn * un - b * vn) * inv_c + cpsi * un
       v = (a * wn * vn + b * un) * inv_c + cpsi * vn
       w = cpsi * wn - a * c
    else
       c = sqrt(1.0_dp - un * un)
       inv_c = 1.0_dp / c
       u = cpsi * un - a * c
       v = (a * un * vn - b * wn) * inv_c + cpsi * vn
       w = (a * un * wn + b * vn) * inv_c + cpsi * wn
    endif

  end subroutine sample_lambertian_direction

  !----------------------------------------------

  subroutine make_MRW_step(id, x,y,z, u,v,w, d, chi_trans, ct)

    use naleat, only : random_isotropic_direction

    integer, intent(in) :: id
    real(kind=dp), intent(inout) :: x, y, z
    real(kind=dp), intent(out) :: u, v, w, ct
    real(kind=dp), intent(in) :: d, chi_trans

    real(dp) :: un, vn, wn, y_sample

    ! Place photon randomly on sphere of radius d (R0) around current position
    call random_isotropic_direction(id, un, vn, wn)
    x = x + un * d
    y = y + vn * d
    z = z + wn * d

    ! Select outgoing Lambertian direction on the emitting MRW sphere surface
    call sample_lambertian_direction(id, un, vn, wn, u, v, w)

    ! Sample y0 from cumulative exit distribution P(t) and solve Eq. (8) of Min et al. (2009):
    ! D = 1 / (3 * chi_trans), cst_ct = 3 / pi^2, where chi_trans [AU^-1] is the local
    ! Rosseland transport mean extinction (including kappa_factor) in B&W mode.
    y_sample = sample_zeta(id)
    ct = -log(y_sample) * cst_ct * chi_trans * (d * d)

    return

  end subroutine make_MRW_step

end module MRW
