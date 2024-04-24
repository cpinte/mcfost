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
