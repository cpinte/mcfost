


The MCFOST radiative transfer code
==================================

- Code: <https://github.com/cpinte/mcfost>
- Docs: <https://mcfost.readthedocs.io/>

[![build-gfortran-macos](https://github.com/cpinte/mcfost/actions/workflows/build-gfortran-macos.yml/badge.svg)](https://github.com/cpinte/mcfost/actions/workflows/build-gfortran-macos.yml) 
[![Documentation Status](https://readthedocs.org/projects/mcfost/badge/?version=stable)](https://mcfost.readthedocs.io/en/stable/?badge=stable)

About
-----

MCFOST is a 3D continuum and line radiative transfer code based on an hybrid Monte Carlo and ray-tracing method. It is mainly designed to study the circumstellar environment of young stellar objects, but has been used for a wide range of astrophysical problems. The calculations are done exactly within the limitations of the Monte Carlo noise and machine precision, i.e. no approximation are used in the calculations. The code has been strongly optimized for speed.


Code Papers
-----------
Core papers :
- Pinte et al. 2006:  https://ui.adsabs.harvard.edu/abs/2006A%26A...459..797P/abstract
- Pinte et al. 2009 : https://ui.adsabs.harvard.edu/abs/2009A%26A...498..967P/abstract

Radiative transfer in atomic lines:
- Tessore et al. 2021 : https://ui.adsabs.harvard.edu/abs/2021A%26A...647A..27T/abstract


Licence
-------

See LICENCE file for usage and distribution conditions.

The code is open source under GPLv3.
We also kindly ask you cite the code papers in scientific publications if you are using the code in your research.

If the code is useful for your research, please get in touch with us.
We welcome scientific collaborations, and hopefully we will be able to help.


Contributing
------------
We welcome contributions, including (but not limited to):

1. Code, via [pull request](https://github.com/cpinte/mcfost/pulls). Please read developer section of user guide for guidelines.
2. Documentation, also by [pull request](https://github.com/cpinte/mcfost/pulls). Docs can be edited in the docs/ directory of the main code.
3. Suggestions for features or bug reports, via the [issue tracker](https://github.com/cpinte/mcfost/issues/new). Please file bugs via github rather than by email.

Visualising your results
------------------------
We suggest to use [pymcfost](https://github.com/cpinte/pymcfost). Alternative options are discussed in the [documentation](https://mcfost.readthedocs.io/en/latest/tools.html).
