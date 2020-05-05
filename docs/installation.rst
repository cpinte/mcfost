Requirements and Installation
=============================

MCFOST is written in Fortran 90 (with a few functions in C++) and parallelized with the Open MP
protocol, i.e. it can use several processors and/or cores on a single
machine. For simplicity, we provide pre-compiled binaries which are statically linked and do not
have any dependencies. Source code is available on github and acces can also be provided if needed.

The following environment is required for the binaries:

-  a 64-bits Unix system, including Linux or MacOS X,
-  any data analysis software capable of handling FITS files.

MCFOST is also available for Xeon Phi but performance is not optimal
yet. Stay tuned.



Binary Installation Procedure
-----------------------------

Unless you need to dig into the source, this is the recommended installation
procedure. The binaries have been compiled with the Intel compiler and should
be optmized for most architectures. The binary also updates itself regularly (unless this options is truened off) making sure you have an up-to-date version.

1. Download the tar ball with the pre-compiled binary:

http://ipag.osug.fr/public/pintec/mcfost/linux/mcfost.tgz
http://ipag.osug.fr/public/pintec/mcfost/mac/mcfost.tgz

2. Extract the file and make it executable::

     $ tar -xvzf mcfost.tgz
     $ chmod +x mcfost

3. Make a directory where you move mcfost (or move mcfost in any directory defined in your shell path)::

     $ mkdir -p ~/mcfost/bin
     $ mv mcfost ~/mcfost/bin

   You can then add this directory to your path with::

   $ setenv PATH ${HOME}/mcfost/bin:${PATH}

   for tcsh or csh shell, or::

   $ export PATH=${HOME}/mcfost/bin:${PATH}

   for bash or bash-like shell.

4. Set the environment variable ``MCFOST_UTILS`` to point to a directory
   where you want mcfost to store its data files.
   E.g. edit your shell startup files to include either::

   $ setenv MCFOST_UTILS ~/mcfost/utils

   for tcsh or csh shell, or::

   $ export MCFOST_UTILS=~/mcfost/utils

   for bash or bash-like shell.

5. You should now be able to run::

     $ mcfost --help


 to get a short list of the available options.

6. Download mcfost data files, current reference parameter files with::

      $ mcfost -setup



Installation from source
------------------------

mcfost source code is hosted on github (as a private directory, public release will be made soon).


1. Clone the repository (ask Christophe for access).

2. Set the following environment variables::

     $ export MCFOST_INSTALL=/my/install/dir
     $ export MCFOST_GIT=1
     $ export MCFOST_AUTO_UPDATE=0

(The last line will prevent the mcfost binary to try to update itself, which is probably not desired if you compile it yourself).

3. Select your compiler::

     $ export SYSTEM=ifort

   or::

     $ export SYSTEM=gfortran

.. important:: If you use gfortran on a mac, please make sure you actually use the GNU gcc/g++/gfortran compiler suite. mcfost has been tested with the homebrew version of the gcc suite. In particular, the default mac gcc is aliased to clang, which has issues compiling some of the libraries.

   To test if you use the system compiler, simply run::

     $ which gcc

   The homebrew version of gcc should be ``/usr/local/bin/gcc``.
   If the output is ``/usr/bin/gcc``, it is likely that you will have issues at the linking stage.

   The solution is to run::

     $ brew link gcc

   which should result in::

     $ which gcc
     /usr/local/bin/gcc


4. Change directory to ``mcfost/lib`` and run the installation script::

   $ ./install.sh

   This should install the required files to ``/my/install/dir/lib`` and
   ``/my/install/dir/include``.
5. Enter the src directory and compile with::

     $ make

6. If you plan to use ``mcfost+phantom`` to perform live rdaiation hydrodynamics calculations, you can compile the mcfost library with::

     $ make libmcfost.a

7. Set the environment variable ``MCFOST_UTILS`` to point to a directory
   where you want mcfost to store its data files.
   E.g. edit your shell startup files to include either::

   $ setenv MCFOST_UTILS ~/mcfost/utils

   for tcsh or csh shell, or::

   $ export MCFOST_UTILS=~/mcfost/utils

   for bash or bash-like shell.

8. You should now be able to run::

     $ mcfost --help


 to get a short list of the available options.

9. Download mcfost data files with::

      $ mcfost -setup


.. note:: mcfost uses the xgboost machine learning library to predict chemical abundances. This features is experimental and xgboost is sometimes tricky to compile with the intel compiler. You can turn the feature off by seting the environement variable `MCFOST_NO_XGBOOST` to yes.


MCFOST_UTILS Environment variable
----------------------------------

MCFOST uses a database of stellar spectra, optical properties and atomic and
molecular data. These files are generally put in a directory named
mcfost/utils, although any name can be used. The environment variable
``MCFOST_UTILS`` must be set to the path name of this directory.

An additional (optional) environment variable ``MY_MCFOST_UTILS`` can be
defined by the user to add his own data files. This has an advantage to
ensure that no personal data files will be overwritten during an update of
the utils directory.


Parallelization
---------------

By default, MCFOST will parallelize itself across all available cpu/cores.
If you want to restrict it to a subset, you can specify the
number of cores to use with the environment variable ``OMP_NUM_THREADS``::

$ setenv OMP_NUM_THREADS <n_cores>

If you wish to disable parallelization entirely, you can use ::

$ setenv OMP_NUM_THREADS 1

Here are scaling results from testing on a 2014 Mac Pro (3 GHz, 8 core Intel
Xeon E5 with 32 GB DDR3 RAM) by Marshall Perrin. This is for calculating
the SED for one particular model file (chosen arbitrarily) from an MCMC
chain prepared by Schuyler Wolff. The scaling is not quite 1/N, but it's
pretty good up to 8 threads, which is the # of true CPU cores this
computer has.
`Hyperthreading* <http://en.wikipedia.org/wiki/Hyper-threading>`__
results in the computer appearing to have 16 virtual cores, but the
performance gain from trying to use these all is marginal.


+----------------+---------------------+--------------------------+
| # of threads   | CPU time used [s]   | Total elapsed time [s]   |
+================+=====================+==========================+
| 1              | 141                 | 141                      |
+----------------+---------------------+--------------------------+
| 2              | 159                 | 79                       |
+----------------+---------------------+--------------------------+
| 4              | 160                 | 40                       |
+----------------+---------------------+--------------------------+
| 8 = n_cores    | 186                 | 23                       |
+----------------+---------------------+--------------------------+
| 16             | 276                 | 18                       |
+----------------+---------------------+--------------------------+

Here are similar results for the ref2.19.para reference parameter file:

+----------------+---------------------+--------------------------+
| # of threads   | CPU time used [s]   | Total elapsed time [s]   |
+================+=====================+==========================+
| 1              | 22                  | 22                       |
+----------------+---------------------+--------------------------+
| 2              | 21                  | 10                       |
+----------------+---------------------+--------------------------+
| 3              | 22                  | 7                        |
+----------------+---------------------+--------------------------+
| 4              | 24                  | 6                        |
+----------------+---------------------+--------------------------+
| 8 = n_cores    | 30                  | 3                        |
+----------------+---------------------+--------------------------+
| 12             | 40                  | 3                        |
+----------------+---------------------+--------------------------+
| 16             | 47                  | 3                        |
+----------------+---------------------+--------------------------+


Setting the stacksize
---------------------

To speed up the calculations, MCFOST stores some arrays privately for each
thread. This means that storage can exceed the default OpenMP stacksize. To
avoid this, include those commeand in your ``.bashrc`` or equivalent::

$ export OMP_STACKSIZE=512M
$ ulimit -s unlimited


Upgrading to New Versions
-------------------------

The mcfost binary will try to update itself every week. An update can be manual
performed via the command ``mcfost -u``. If you wish to update to new binary
version between releases, you can do so by forcing the update via ``mcfost -fu``.

The ``MCFOST_UTILS`` data can updated via ``mcfost -update-utils``.

MCFOST will check for updates automatically at
start-up if the last update is older than 7 days (this should take less
than 1 second). This behaviour can be changed by setting the environment
variable ``MCFOST_AUTO_UPDATE`` to an integer defining the number of days
between which mcfost will check for updates. If ``MCFOST_AUTO_UPDATE`` is
set to 0, mcfost will not check for updates automatically (this is the
recommended behaviour is you are using the source code).


If you are using the source code, MCFOST can be updated via::

    $ git pull
    $ make
