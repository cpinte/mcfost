MCFOST test suite
=================


A test suite is available to ensure that mcfost behaves as expected. The released binary passes the test suite. This test suite is only useful if you compile mcfost yourself.

So far, functional tests only have been implemented, no unit tests, ie
the goal is to ensure that mcfost calculations remain consistent with results from previous versions (which were benchmarked against other codes).

The test suite is performed in 2 steps.

1. Downloading the the reference solution::

     $ cd ~/mcfost/test_suite
     $ ./get_test_data.sh

2. Calculations of the same models with the current version and validation of results::

     $ cd ../src
     $ git checkout master
     $ make
     $ cd ../test_suite
     $ pytest -v

Running the reference and test solution takes ~ 14 minutes each on a macbook pro (Intel i7 4 cores at 2.8Ghz)

 .. note:: Comparisons with version earlier than v3.0.34 will fail as some definitions regarging the orientation of the images, and the Stokes vectors were changed. Those changes were validated manually.

Benchmarks
----------

The parameter files of the benchmark of `Pinte et al. 2009 <https://ui.adsabs.harvard.edu/abs/2009A&A...498..967P>`_
(very optically thick disc, anisotropic scattering and polarisation, midplane optical depths from 10\ :sup:`3`
to 10\ :sup:`6`) are in ``test_suite/benchmarks/Pinte2009``. The star of the benchmark is a blackbody: run them
with the ``-star_bb`` option. See the README in that directory for the model and the comparisons.

The parameter files of the benchmark of `Pascucci et al. 2004 <https://ui.adsabs.harvard.edu/abs/2004A&A...417..793P>`_
(2D disc with isotropic scattering, visual optical depths from 0.1 to 100), with the results of the five codes
of the paper, are in ``test_suite/benchmarks/Pascucci2004``. They also need the ``-star_bb`` option.
