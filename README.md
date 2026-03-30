


The MCFOST radiative transfer code
==================================

<p align='center'>
  <br/>
  <img src="https://github.com/cpinte/mcfost/blob/main/logo/mcfost_logo.png" width="300" height="300">
  <br/>
</p>

- Code: <https://github.com/cpinte/mcfost>
- Docs: <https://mcfost.readthedocs.io>

[![test-suite](https://github.com/cpinte/mcfost/actions/workflows/test-suite.yml/badge.svg)](https://github.com/cpinte/mcfost/actions/workflows/test-suite.yml?query=branch%3Amain++)
[![Documentation Status](https://readthedocs.org/projects/mcfost/badge/?version=latest)](https://mcfost.readthedocs.io/en/latest/)

About
-----

MCFOST is a 3D continuum and line radiative transfer code based on an hybrid Monte Carlo and ray-tracing method. It is mainly designed to study the circumstellar environment of young stellar objects, but has been used for a wide range of astrophysical problems. The calculations are done exactly within the limitations of the Monte Carlo noise and machine precision, i.e. no approximation are used in the calculations. The code has been strongly optimized for speed.

Code of conduct
---------------
If you wish to use the code, please make sure you agree to adhere to the [code of conduct](https://github.com/cpinte/mcfost?tab=coc-ov-file).

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

1. Code, via [pull request](https://github.com/cpinte/mcfost/pulls). Please read the developer section of the user guide for guidelines. 
We use the [pre-commit](https://pre-commit.com) framework to automatically fix
some coding bad practices.
It is recommended to install pre-commit by running the following commands from the top level of the repo
```shell
python3 -m pip install pre-commit
pre-commit install
```
2. Documentation, also by [pull request](https://github.com/cpinte/mcfost/pulls). Docs can be edited in the docs/ directory of the main code. 
3. Suggestions for features or bug reports, via the [issue tracker](https://github.com/cpinte/mcfost/issues/new). Please file bugs via github rather than by email.

Questions?
----------

Discussions about the code and its use have moved [here](https://github.com/cpinte/mcfost/discussions).

Visualising your results
------------------------
We suggest to use [pymcfost](https://github.com/cpinte/pymcfost). Alternative options are discussed in the [documentation](https://mcfost.readthedocs.io/en/latest/tools.html).
