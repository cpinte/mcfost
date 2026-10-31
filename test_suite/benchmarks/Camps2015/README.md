# Camps et al. 2015 benchmark: stochastically heated grains

Driver for the benchmark of Camps et al. 2015, A&A 580, A87
(arXiv:1506.05304), https://shg.ugent.be: emissivity of a mixture of
silicate, graphite and PAH grains (BARE-GR-S model of Zubko et al. 2004)
exposed to imposed radiation fields, with stochastic heating. The first time
the driver runs, it extracts the dust model and the results of the paper into
`data/` (not versioned, 51 MB) from the copy of the benchmark data in the
MCFOST utils (`$MCFOST_UTILS/benchmarks/Camps2015`, archives of the site
unchanged, with their SHA-256), or downloads them from the site if the utils
do not have them.

Input fields (the mean intensity of the dust emission itself is ignored):

- the Mathis et al. 1983 interstellar radiation field, Eq. 11 of the paper,
  scaled by U = 1e-4 ... 1e6;
- diluted blackbodies, 3000 to 18000 K, dilution factors 8.28e-12 to 9.42e-16.

Running, from this directory or any other:

    shg_benchmark.py run -m <mcfost binary> [-w workdir] [-f fields] [-g grains] [-t n_T] [-j jobs]
    shg_benchmark.py compare [-w workdir]

`run` writes one mcfost model per grain type and field (51 runs, about 20 s
on 8 cores), `compare` prints, for each field and grain type, the ratio of
the emitted to the absorbed power per grain size, and the ratio of the mcfost
emissivity to DustEM, SKIRT and CRT from 3 to 1000 um (where the reference is
above 1e-3 of its peak): median [min, max].

## Model

- one thin spherical shell (1 cell) at distance R of a 1 Rsun blackbody star,
  1e-10 Msun of dust, a single grain type with heating method 3 (nRE),
  64, 62 and 24 log-spaced sizes for silicates, graphite and PAHs (as the
  MCFOST results of the paper), 300 temperature bins (`-nT`);
- optical properties: the tabulated files of the benchmark in the Misselt
  format, read by mcfost when the optical indices file name starts with
  `Misselt` (first line added by the driver: name of the calorimetry file);
- `-benchmark_SHG` replaces the Monte Carlo radiation field by the field of
  the benchmark in `Temp_nRE`: the Mathis field scaled by `-Mathis_field U`,
  or the star diluted at the distance of the cell, W = (R_star / 2R)^2 (the
  driver places the shell at R = R_star / (2 sqrt(D)) for a dilution factor D).
  The nRE grains are not re-emitted (no iteration);
- wavelength grid of the temperature calculation: the 1201 log-spaced
  wavelengths of the benchmark, from 0.02 um (the stars have no flux below,
  the fields are negligible there);
- emissivity: computed by the driver from the equilibrium temperatures and
  temperature probabilities of `data_th/Temperature_nRE.fits.gz`, with the
  cross-sections interpolated as mcfost does, and the size distributions of
  the benchmark (`SizeInputs`) integrated over the mcfost size bins.

## Results

October 2026 (branch `pah_nre`), 300 temperature bins, ratio of the mcfost
emissivity to DustEM from 3 to 1000 um, over the 17 fields:

| Grain | mcfost / DustEM | Median over the fields | Emitted / absorbed power per size | Other codes / DustEM |
|---|---|---|---|---|
| Silicate | 0.93-1.04 | 1.005-1.011 | 0.97-1.00 | 0.48-1.32 |
| Graphite | 0.91-1.03 | 0.996-1.001 | 0.99-1.00 | 0.89-1.15 |
| PAH      | 0.88-1.05 | 0.984-1.020 | 0.99-1.00 | 0.62-1.42 |

With 1000 temperature bins, the PAH results change by less than 2 %.
Individual large PAHs (25-50 A, single-grain solution of the paper) need
1000 bins: the 3.3 um emission of a 50 A PAH in the Mathis field is 19
times DustEM with 100 bins (the default), 2.8 with 300 bins, 1.4 with 1000
bins (SKIRT: 1.8, CRT: 1.2), as the fixed 1-3000 K grid is coarse compared
with the energy of one photon.

The MCFOST results distributed with the paper (2015) predate these
corrections and are not used as a reference.

## Test

`test_mcfost.py::test_camps2015_SHG_benchmark` (long test, local runs only)
runs the 51 models and checks the energy balance of each grain size (within
6 %), and the emissivity against DustEM (median within 3 %, every wavelength
within 15 %).
