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
polarisation maps (circular polarisation was not compared). The results of
all codes were published at http://www.astro.ex.ac.uk/people/cpinte/benchmark/
(paper, section 3); that page may no longer be online.

The photon numbers are those used for MCFOST tests, not the values of the
paper (not specified there): increase `nbr_photons_image` for Monte Carlo
images and polarisation at high inclinations.
