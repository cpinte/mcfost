# Pinte et al. 2009 benchmark

Parameter files for the continuum radiative transfer benchmark of
Pinte et al. 2009, A&A 498, 967 (arXiv:0903.1231): a very optically thick
disc with anisotropic scattering and polarisation, compared between MCFOST,
TORUS, MCMax and Pinball (and ProDiMo, RADMC and RADICAL with isotropic or
no scattering).

| File | Dust mass (Msun) | Midplane optical depth at 0.81 um |
|---|---|---|
| `Pinte2009_tau=1e3.para` | 3e-8 | 1.22e3 |
| `Pinte2009_tau=1e4.para` | 3e-7 | 1.22e4 |
| `Pinte2009_tau=1e5.para` | 3e-6 | 1.22e5 |
| `Pinte2009_tau=1e6.para` | 3e-5 | 1.22e6 |

Model (section 3 of the paper), common to the four files:

- star: uniformly radiating blackbody, 4000 K, 2 Rsun, at 140 pc;
- disc: Rin = 0.1 au, Rout = 400 au with sharp edges, surface density
  r^-1.5, Gaussian vertical profile with h = 10 au (r/100 au)^1.125;
- dust: 1 um astronomical silicate spheres (Weingartner & Draine 2001,
  `Draine_Si_sUV.dat`), 3.5 g/cm3, Mie theory;
- images: 1 um, 251 x 251 pixels over 900 au (25.61 mas/pixel), Stokes
  parameters, 10 inclinations equally spaced in cos(i) = 0.05, ..., 0.95
  (i = 18.2 to 87.1 deg).

The star must be a blackbody: the parameter file can only select a stellar
spectrum, so use the `-star_bb` option for every run:

    mcfost Pinte2009_tau=1e3.para -star_bb                # temperature and SED
    mcfost Pinte2009_tau=1e3.para -star_bb -img 1.0       # images and polarisation maps (ray-tracing)
    mcfost Pinte2009_tau=1e3.para -star_bb -img 1.0 -mc   # also Monte Carlo images

`mcfost Pinte2009_tau=1e3.para -star_bb -op 0.81` prints the midplane
optical depth ("Integ tau in midplane"), which must match the table above.

The benchmark compared temperatures, SEDs, scattered light images and linear
polarisation maps (circular polarisation was not compared). The description of
the benchmark and the results of all codes (MCFOST, TORUS, MCMax, Pinball,
RADMC, ProDiMo, RADICAL; anisotropic, isotropic and no scattering cases) were
published on the benchmark webpage, first at Exeter
(http://www.astro.ex.ac.uk/people/cpinte/benchmark/, paper section 3), now at
https://ipag.osug.fr/~pintec/benchmark/index.shtml. A copy is distributed with
the MCFOST utils, in `$MCFOST_UTILS/benchmarks/Pinte2009/` (run
`mcfost -update_utils` if the directory is missing).

The photon numbers are those used for MCFOST tests, not the values of the
paper (not specified there): increase `nbr_photons_image` for Monte Carlo
images and polarisation at high inclinations.

## Tests

`test_mcfost.py::test_pinte2009_benchmark` (long test suite, `pytest --run-slow`)
runs the tau = 1e3 and 1e6 models and compares the SEDs at the 10 inclinations
with MCMax and MCFOST 2009, the midplane temperatures with MCMax, TORUS and
MCFOST 2009, and the vertical temperatures at 0.2 and 200 au with MCMax and
MCFOST 2009. The tau = 1e6 model is run with `-diff_approx`: the temperature of
the optically thick midplane, not reached by the packets, is computed with the
diffusion approximation (it is set to 1 K otherwise). Largest differences
(MCFOST 4.1.16, October 2026):

| Optical depth | SED per wavelength (i <= 70-76 deg) | SED median, all inclinations | Temperatures |
|---|---|---|---|
| 1e3 | 16 % | 8 % (81.4 deg) | 7.6 % |
| 1e6 | 32 % (short wavelengths, 18 deg) | 3 % | 19 % (inner au, MCMax and TORUS) |

For tau = 1e6, the test also computes the 1 um image and compares it with
MCMax and MCFOST 2009 at the 10 inclinations (and TORUS at 70 and 87 deg),
excluding the stellar pixels: scattered flux (largest difference 7.5 %, MCMax),
cuts along the major and minor axes (15 % up to 64 deg, medians 14 %), and the
flux-weighted linear polarisation degree (24 %: MCMax is 10-20 % lower, the
MCFOST 2009 maps, computed before the 2026 fix of multiple scattering, about
8 % higher).

`test_mcfost.py::test_pinte2009_temperature` (long test suite) runs the
tau = 1e3 model with isotropic scattering (`-isotropic`) and without scattering
(`-no_scattering`), and compares the midplane temperatures and the vertical cuts
at 0.2 and 200 au with RADMC, ProDiMo, TORUS and MCFOST 2009 (isotropic) and
RADICAL, ProDiMo and MCFOST 2009 (no scattering): medians within 3 %, largest
differences 14 % (isotropic) and 18 % (no scattering, near the edges of the disc
where the RADICAL and ProDiMo grids are coarser).
