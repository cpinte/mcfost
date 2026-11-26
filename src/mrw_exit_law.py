#!/usr/bin/env python3
"""Exit direction of the MRW walk: exact law, and the table of MRW_exit_table.f90.

A packet starting at the centre of a sphere of radius gamma (transport mean free paths) with an
isotropic direction, conservative isotropic scattering, leaves it with an exit cosine mu (to the
outward normal) distributed as
    p(mu) = exp(-gamma) delta(mu - 1) + (1 - exp(-gamma)) p_s(mu, gamma),
where the first term is the packets that leave without interaction, and p_s the others. p_s is the
emergent intensity of a point source at the centre of the sphere, times mu:
    p_s(mu) ~ mu I_s(mu),  I_s(mu) = int_0^{2 gamma mu} S(r(s)) exp(-s) ds,  r^2 = gamma^2 - 2 gamma mu s + s^2,
with S = J the mean intensity (albedo 1), J = J_unc + J_c, J_unc = exp(-r) / (4 pi r^2), and u = r J_c
solution of the radial transport equation (integral equation with an E1 kernel):
    u(r) = 1/2 int_0^gamma [exp(-r') / (4 pi r') + u(r')] [E1(|r - r'|) - E1(r + r')] dr'.
Large gamma: p_s tends to the Milne law mu H(mu) / alpha_1, with a 1/gamma curvature correction.
Small gamma: a peak at mu -> 1, ~ (pi/2) gamma exp(-gamma) / ((1 - exp(-gamma)) sqrt(1 - mu^2)), from the
light scattered once close to the centre. The law agrees with exact walks within their noise
(mrw_study, 1e6 walks per gamma, gamma = 2-30).

Table: mu as a function of (u, log gamma), with xi = sin^2(pi u / 2) the cumulative probability of
p_s (uniform u resolves both ends: 1 - mu ~ (1 - xi)^2 at 1; at 0, grazing chords are short, I_s -> 0
and mu ~ xi^(1/3)). Bilinear interpolation error at gamma between the grid points: mean 6e-5 in mu,
max 1e-3 for xi > 1e-3, up to 6e-3 for xi < 1e-3 (mu < 0.08, 0.1 % of the packets).
Usage: python3 mrw_exit_law.py            writes MRW_exit_table.f90 (next to this script)
       python3 mrw_exit_law.py --check    interpolation error at gamma between the grid points
"""
import os
import subprocess
import sys
import numpy as np
import scipy
from scipy.special import expn
from scipy.integrate import quad

GAM_MIN, GAM_MAX = 2.0, 100.0
N_GAM, N_U = 33, 129


def _kernel_cell(r, a, b):
    """1/2 int_a^b [E1(|r - t|) - E1(r + t)] dt (antiderivative of E1: -E2)."""
    if b <= r:
        i1 = expn(2, r - b) - expn(2, r - a)
    elif a >= r:
        i1 = expn(2, a - r) - expn(2, b - r)
    else:
        i1 = (1.0 - expn(2, r - a)) + (1.0 - expn(2, b - r))
    return 0.5 * (i1 - (expn(2, r + a) - expn(2, r + b)))


def source_function(R, n=600):
    """Collided mean intensity J_c at the cell centres (cells refined at the centre and surface)."""
    e = R * (0.5 - 0.5 * np.cos(np.pi * np.linspace(0, 1, n + 1)))
    c = 0.5 * (e[1:] + e[:-1])
    K = np.array([[_kernel_cell(ri, a, b) for a, b in zip(e[:-1], e[1:])] for ri in c])

    def src(r):
        f = lambda t: np.exp(-t) / (4 * np.pi * t) * 0.5 * (expn(1, abs(r - t)) - expn(1, r + t)) if t > 0 else 0.0
        return quad(f, 0, R, points=[r], limit=400)[0]

    u = np.linalg.solve(np.eye(n) - K, np.array([src(ri) for ri in c]))
    return c, u / c


def emergent(R, mu, c, Jc):
    """mu I_s(mu): emergent intensity of the scattered light times mu (unnormalised p_s)."""
    def f(s, m):
        r = np.sqrt(max(R * R - 2 * R * m * s + s * s, 1e-300))
        return (np.interp(r, c, Jc) + np.exp(-r) / (4 * np.pi * r * r)) * np.exp(-s)
    out = []
    for m in np.atleast_1d(mu):
        smax = 2 * R * m
        pts = [R * m] if 0 < R * m < smax else None
        out.append(m * quad(f, 0, smax, args=(m,), points=pts, limit=400, epsrel=1e-9)[0])
    return np.array(out)


def inverse_cdf(R):
    """mu(xi) for p_s at gamma = R, as an interpolating function."""
    c, Jc = source_function(R, 600 if R <= 30 else 900)
    mu = np.concatenate([np.linspace(1e-4, 0.9, 400), 1 - np.geomspace(0.1, 1e-8, 160)])
    p = emergent(R, mu, c, Jc)
    # p ~ mu at 0 (linear), integrable (1 - mu)^(-1/2) at 1: trapezoid, plus the analytic end interval
    x = np.concatenate([[0.0], mu])
    pp = np.concatenate([[0.0], p])
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (pp[1:] + pp[:-1]) * np.diff(x))])
    tail = 2.0 * p[-1] * (1.0 - mu[-1])          # int of a (1 - mu)^(-1/2) over the last interval
    cdf = np.concatenate([cdf, [cdf[-1] + tail]]) / (cdf[-1] + tail)
    x = np.concatenate([x, [1.0]])
    return lambda xi: np.interp(xi, cdf, x)


def gamma_grid():
    return np.geomspace(GAM_MIN, GAM_MAX, N_GAM)


def u_grid():
    return np.linspace(0.0, 1.0, N_U)


def table():
    xi = np.sin(0.5 * np.pi * u_grid()) ** 2
    T = np.array([inverse_cdf(g)(xi) for g in gamma_grid()])
    T[:, 0], T[:, -1] = 0.0, 1.0
    return T


def script_hash():
    """git blob hash of this script (git hash-object), identifies the version that wrote the table."""
    try:
        return subprocess.run(["git", "hash-object", os.path.abspath(__file__)], capture_output=True,
                              text=True, check=True).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"


def write_fortran(T, path=os.path.join(os.path.dirname(os.path.abspath(__file__)), "MRW_exit_table.f90")):
    lines = ["! Generated by mrw_exit_law.py: do not edit. Regenerate with: python3 mrw_exit_law.py",
             f"! mrw_exit_law.py blob {script_hash()} (git hash-object), numpy {np.__version__}, scipy {scipy.__version__}",
             f"! Grid: {N_U} values of u in [0, 1], {N_GAM} values of gamma in [{GAM_MIN:g}, {GAM_MAX:g}] (log spaced);",
             "! bilinear interpolation error in mu: mean 6e-5, max 1e-3 for xi > 1e-3 (mrw_exit_law.py --check)",
             "! Exit cosine of the packets that interact in the MRW sphere (exact law, see mrw_exit_law.py):",
             "! mu_exit_table(i_u, i_gam) at u = (i_u - 1) / (n_u_exit - 1), xi = sin^2(pi u / 2) the cumulative",
             "! probability, and gamma = gam_min_exit (gam_max_exit / gam_min_exit)^((i_gam - 1) / (n_gam_exit - 1))",
             "module MRW_exit_table",
             "  use mcfost_env, only : dp",
             "  implicit none",
             f"  integer, parameter :: n_u_exit = {N_U}, n_gam_exit = {N_GAM}",
             f"  real(dp), parameter :: gam_min_exit = {GAM_MIN:.1f}_dp, gam_max_exit = {GAM_MAX:.1f}_dp",
             "  real(dp), dimension(n_u_exit, n_gam_exit), save, protected :: mu_exit_table",
             "  ! one DATA statement per gamma (continuation lines below the limit of the standard)"]
    for j in range(T.shape[0]):
        vals = [f"{v:.10f}_dp" for v in T[j]]
        lines.append(f"  data mu_exit_table(:, {j + 1}) / &")
        for i in range(0, len(vals), 5):
            end = ", &" if i + 5 < len(vals) else " /"
            lines.append("       " + ", ".join(vals[i:i + 5]) + end)
    lines += ["end module MRW_exit_table", ""]
    open(path, "w").write("\n".join(lines))


def interp(T, gam, xi):
    """Bilinear interpolation as in MRW.f90 (u from xi, log gamma)."""
    gg = gamma_grid()
    g = np.clip(gam, GAM_MIN, GAM_MAX)
    t = np.log(g / GAM_MIN) / np.log(GAM_MAX / GAM_MIN) * (N_GAM - 1)
    j = min(int(t), N_GAM - 2); w = t - j
    u = 2.0 / np.pi * np.arcsin(np.sqrt(xi)) * (N_U - 1)
    i = np.minimum(u.astype(int), N_U - 2); a = u - i
    row = lambda k: (1 - a) * T[k, i] + a * T[k, i + 1]
    return (1 - w) * row(j) + w * row(j + 1)


if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    T = table()
    if "--check" in sys.argv:
        gg = gamma_grid()
        xs = np.linspace(0, 1, 20001)[1:-1]
        for g in np.sqrt(gg[1:] * gg[:-1])[[0, 1, 3, 8, 16, 31]]:
            err = np.abs(interp(T, g, xs) - inverse_cdf(g)(xs))
            print(f"gamma {g:7.3f}: max |mu error| {err.max():.1e}, mean {err.mean():.1e}")
    else:
        write_fortran(T)
        print("wrote MRW_exit_table.f90")
