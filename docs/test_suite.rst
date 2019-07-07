MCFOST test suite
=================


A test suite is available to ensure that mcfost behaves as expected. The released binary passes the test suite. This test suite is only useful if you compile mcfost yourself.

So far, functional tests only have been implemented, no unit tests, ie
the goal is to ensure that mcfost calculations remain consistent with results from previous versions (which were benchmarked against other codes).

The test suite is performed in 2 steps.

1. Calculations of the reference solution::

     $ cd ~/mcfost/src
     $ git checkout ef6c1060d779
     $ make
     $ cd ../test_suite
     $ ./compute_test_suite.sh

2. Calculations of the same models with the current version and validation of results::

     $ cd ../src
     $ git checkout master
     $ make
     $ cd ../test_suite
     $ pytest -v

Running the reference and test solution takes 10 to 12 minutes each on a macbook pro (Intel i7 4 cores at 2.8Ghz)

 .. note:: Comparisons with version earlier than v3.0.44 will fail as some definitions regarging the orientation of the images, and the Stokes vectors were changed. Those changes were validated manually.
           Some combinations of parameters were not possible before commit ef6c1060d779, which is now the reference for the calculations. Comparison with earlier versions were performed manually.
