# Pascucci et al. 2004 benchmark

Parameter files for the 2D continuum radiative transfer
benchmark of Pascucci et al. 2004, A&A 417, 793 (arXiv:astro-ph/0402357):
temperature structure and SED of a passive disc with isotropic scattering, for
visual optical depths from 0.1 to 100. The parameter files are the version 3.0
files used for MCFOST previously (`src/Pascucci_3.0*.para`), converted to the
4.1 format.

| File | Dust mass (Msun) | Packets for the temperature |
|---|---|---|
| `Pascucci_tau=0.1.para` | 1.1e-7 | 1.28e6 |
| `Pascucci_tau=1.para`   | 1.1e-6 | 1.28e6 |
| `Pascucci_tau=10.para`  | 1.1e-5 | 1.28e7 |
| `Pascucci_tau=100.para` | 1.1e-4 | 1.28e8 |

Model, common to the four files:

- star: blackbody, 5800 K, 1 Rsun;
- disc: rho(r,z) = rho0 (r/500 au)^-1 exp(-pi/4 (z/h(r))^2), h(r) = 125 au
  (r/500 au)^1.125, from 1 to 1000 au. In MCFOST terms: Gaussian scale height
  125 sqrt(2/pi) = 99.736 au at 500 au, flaring exponent 1.125, surface density
  exponent +0.125;
- dust: 0.12 um astronomical silicate grains, 3.6 g/cm3, isotropic scattering;
- SEDs at the 61 wavelengths of `Pascucci.lambda` (MCFOST utils), at
  inclinations 12.5, 42.5 and 77.5 deg.

Running:

    mcfost Pascucci_tau=1.para -star_bb

- `-star_bb`: the star of the benchmark is a blackbody, the parameter file can
  only select a stellar spectrum.
- MCFOST recognises the benchmark from the name of the parameter file (it must
  start with `Pascucci`) and then forces isotropic scattering and a grain
  density of 3.6 g/cm3.
- The parameter files compute the SEDs at 12.5 and 77.5 deg. For 42.5 deg, set
  the inclination line to `42.5  42.5  1  F` (inclinations of a multi-inclination
  run are equally spaced in cos(i)).
- A spherical grid can be tested by changing the grid lines to `2` and
  `100 50 1 20`, as in `src/Pascucci_spherical_3.0.para`.

Dust properties: in benchmark mode, MCFOST replaces the Mie cross-sections by
those tabulated for the benchmark (`Pascucci_optSi.dat` in the MCFOST utils,
interpolated in log-log in wavelength), with isotropic scattering. The optical
indices file of the parameter files (`Draine_Si.dat`) is still read, but its
cross-sections are not used.

Optical depth: `mcfost Pascucci_tau=1.para -star_bb -op 0.55` prints the
midplane optical depth ("Integ tau in midplane"): 0.1013, 1.013, 10.13 and
101.3 for the four files (analytic value with the tabulated cross-sections:
0.987 for 1.1e-6 Msun).

## Reference results

The results of the five codes of the paper are distributed with the MCFOST
utils, in `$MCFOST_UTILS/benchmarks/Pascucci2004/` (added in October 2026:
run `mcfost -update_utils` if the directory is missing)
(MC3D,
MCTRANSF, RADICAL, RADMC, STEINRAY), for the four optical depths (`01tau`,
`1tau`, `10tau`, `100tau`), in the native formats in which they were
distributed for the benchmark: temperature structures (e.g. `dtemp.dat`,
`tm.dat`, `bench_tm.*`, `bencT*`) and SEDs at the three inclinations
(e.g. `allseds.dat`, `sp`, `bench_sp.*`, `sp_ray_new`). The formats differ
between codes; see the headers of the files and the paper for units.

The SEDs in `allseds.dat` (RADMC, RADICAL) are lambda L_lambda in erg/s:
MCFOST `sed_rt.fits.gz` (lambda F_lambda in W/m2 at 140 pc) multiplied by
1e3 x 4 pi d^2. Check with MCFOST 4.1.16 and the tabulated cross-sections (`Pascucci_tau=1.para`,
27/10/2026),
median ratio MCFOST / reference:

| Inclination | < 1 um | 1-10 um | 10-100 um | > 100 um |
|---|---|---|---|---|
| 12.5 deg (RADMC, RADICAL) | 0.90 | 0.98 | 1.02 | 1.03 |
| 77.5 deg (RADMC, RADICAL) | 1.06 | 0.98 | 1.02 | 1.05 |

The difference below 1 um is not explained yet: it is the same with the Mie
cross-sections, and the direct stellar flux of MCFOST is within 2 % of the
blackbody at these wavelengths.

## Tests

`test_mcfost.py::test_pascucci_benchmark` runs these models and compares the
SEDs above 1 um and the midplane temperatures with RADMC and RADICAL:
tau = 0.1 and 1 in the default test suite (pull requests), tau = 10 and 100
with `pytest --run-slow`. The tests are skipped if the results are not in the
utils.
