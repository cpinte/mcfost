module MRW
  use mcfost_env, only : dp
  use constants, only : pi

  implicit none

  integer, parameter :: n = 10000
  real(dp), dimension(n), save :: zeta, y_MRW

  ! gamma: minimum sphere radius in transport mean free paths. In the thick midplane cells just
  ! inside the dark zone (Pinte 2009 tau = 1e5-1e6), T is low by 7 % (gamma = 2) and 2-3 % (5)
  ! without the extrapolation length and exit law below, 3 % and <= 1 % with them
  real(dp) :: gamma_MRW = 5.0_dp
  real(dp), parameter :: cst_ct = 3.0_dp / (pi * pi)
  ! Extrapolation length of the diffusion solution at a free boundary, in mean free paths:
  ! Hopf constant q(inf) = 0.7104 of the Milne problem (Chandrasekhar 1950; 2/3 in the
  ! Eddington approximation). The exact transport walk leaves a sphere of radius R0 as the
  ! diffusion solution with an absorbing boundary at R0 + z_extrap / chi: mean path
  ! (gamma + 0.71)^2 / 2 instead of gamma^2 / 2 mean free paths. Using R0 underestimates the
  ! time spent in the sphere, hence the deposited energy (cause given by Min et al. 2009, Sect. 3.2)
  real(dp), parameter :: z_extrap_MRW = 0.7104_dp
  ! Chaining of the MRW steps while d * chi_trans > gamma. A chained step restarts from the exit
  ! point as if the packet had just been scattered isotropically and loses its outward
  ! direction: with z_extrap on every step the chain diffuses too slowly (D = gamma^2 / (3
  ! (gamma + z_extrap)^2) instead of 1/3) and over-deposits, +16 % (gamma 2) to +6 % (gamma 5)
  ! in a grey slab (mrw_study/test_slab_hybrid.py). Without chaining, every walk starts after
  ! n_min_interactions_MRW real interactions and the MRW reproduces exact transport to 0.5 %.
  !   0: no chaining (default), 1: chained, z_extrap on every step,
  !   2: chained, z_extrap on the first step only (the chained steps are then consistent, D = 1/3)
  integer :: mrw_chain_mode = 0
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
    use messages, only : error
#include "sprng_f.h"

    integer, intent(in) :: id

    real(dp) :: zeta_random

    ! the table is initialised in init_MRW_cells, before the transport (y_MRW(n) = 1 once set)
    if (y_MRW(n) < 1.0_dp) call error("MRW: zeta table not initialised")
    zeta_random = real(sprng(stream(id)), kind=dp)
    sample_zeta = interp_y_from_zeta(zeta_random)

  end function sample_zeta

  !----------------------------------------------

  pure real(dp) function exit_peak_fraction(gam) result(f)
    ! Fraction of 9 mu^8 in the exit law of the walks that interact in a sphere of gam transport
    ! mean free paths, fitted to <mu> of the exact walk for gam = 2-20 (mrw_study: <mu> within
    ! 0.0015, <mu^2> within 0.005); tends to 1.06 / gam
    real(dp), intent(in) :: gam

    f = min(1.063_dp * gam / (gam * gam - 1.284_dp * gam + 2.702_dp), 1.0_dp)

  end function exit_peak_fraction

  !----------------------------------------------

  subroutine sample_exit_direction(id, un, vn, wn, u, v, w, gam)
    ! Sample the direction of a packet leaving the MRW sphere, around the unit outward normal
    ! (un, vn, wn), phi = 2*pi*xi_2. Large spheres: mu = cos(psi) from p(mu) = mu (1 + 1.5 mu)
    ! on [0,1], the emergent flux of an optically thick medium in the Eddington approximation
    ! (limb darkening, <mu> = 0.71). A Lambertian exit (p = 2 mu, <mu> = 2/3, Min et al. 2009) is
    ! too oblique. Sphere of gam transport mean free paths (optional): the exact walk also
    ! escapes without interaction (probability exp(-gam), mu = 1), and its diffuse part is more
    ! forward-peaked: mu (1 + 1.5 mu) mixed with 9 mu^8 in proportion exit_peak_fraction(gam)
    ! (<mu> = 0.83 at gam 2, 0.76 at 5, 0.73 at 10)
    use naleat, only : stream
#include "sprng_f.h"

    integer, intent(in) :: id
    real(dp), intent(in) :: un, vn, wn
    real(dp), intent(out) :: u, v, w
    real(dp), intent(in), optional :: gam

    real(dp) :: rand1, rand2, rand3, p_ball, p_peak, cpsi, spsi, phi, a, b, c, inv_c, f
    integer :: k

    rand1 = real(sprng(stream(id)), kind=dp)
    rand2 = real(sprng(stream(id)), kind=dp)

    rand3 = 1.0_dp
    p_ball = 0.0_dp
    p_peak = 0.0_dp
    if (present(gam)) then
       rand3 = real(sprng(stream(id)), kind=dp)
       p_ball = exp(-gam)
       p_peak = p_ball + (1.0_dp - p_ball) * exit_peak_fraction(gam)
    endif

    if (rand3 < p_ball) then
       cpsi = 1.0_dp
    else if (rand3 < p_peak) then
       cpsi = rand1**(1.0_dp / 9.0_dp)
    else if (rand1 <= 0.0_dp) then
       cpsi = 0.0_dp
    else
       ! Solve the CDF (mu^2 + mu^3) / 2 = xi_1 by Newton: the function is convex and the start
       ! point is above the root, so the iterations decrease monotonically
       cpsi = min(sqrt(2.0_dp * rand1), 1.0_dp)
       do k = 1, 50
          f = cpsi * cpsi * (1.0_dp + cpsi) - 2.0_dp * rand1
          if (f <= 1.0e-15_dp) exit
          cpsi = cpsi - f / (cpsi * (2.0_dp + 3.0_dp * cpsi))
       enddo
       cpsi = min(max(cpsi, 0.0_dp), 1.0_dp)
    endif
    spsi = sqrt(max(0.0_dp, 1.0_dp - cpsi * cpsi))
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

  end subroutine sample_exit_direction

  !----------------------------------------------

  subroutine make_MRW_step(id, x,y,z, u,v,w, d, chi_trans, ct, z_extrap)
    ! z_extrap: extrapolation length [transport mean free paths], default z_extrap_MRW

    use naleat, only : random_isotropic_direction

    integer, intent(in) :: id
    real(kind=dp), intent(inout) :: x, y, z
    real(kind=dp), intent(out) :: u, v, w, ct
    real(kind=dp), intent(in) :: d, chi_trans
    real(kind=dp), intent(in), optional :: z_extrap

    real(dp) :: un, vn, wn, y_sample, z_ex

    ! Place photon randomly on sphere of radius d (R0) around current position
    call random_isotropic_direction(id, un, vn, wn)
    x = x + un * d
    y = y + vn * d
    z = z + wn * d

    ! Outgoing direction on the sphere surface
    call sample_exit_direction(id, un, vn, wn, u, v, w, d * chi_trans)

    ! Sample y0 from cumulative exit distribution P(t) and solve Eq. (8) of Min et al. (2009):
    ! D = 1 / (3 * chi_trans), cst_ct = 3 / pi^2, where chi_trans [AU^-1] is the local
    ! Rosseland transport mean extinction (including kappa_factor) in B&W mode. The absorbing
    ! boundary of the diffusion solution is at the extrapolated radius d + z_extrap / chi_trans
    ! (the packet is still moved by d): reproduces the path distribution of the exact walk
    z_ex = z_extrap_MRW
    if (present(z_extrap)) z_ex = z_extrap
    y_sample = sample_zeta(id)
    ct = -log(y_sample) * cst_ct * chi_trans * (d + z_ex / chi_trans)**2

    return

  end subroutine make_MRW_step

end module MRW
