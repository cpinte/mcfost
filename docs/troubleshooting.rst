Troubleshooting
===================


Most common error messages
--------------------------

Error: " WARNING : first cell is in diffusion approximation zone. Increase spatial grid resolution"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

It is an optical depth problem, which require a proper sampling of the grid.
The solution is
to increase n_rad_in which is the number of subdivision of the first
cell
The solution is to increase ``n_rad_in``. You may need to try a few values
empirically to find one that works. Around 30 or 40 is recommended, ie
something like that:
``150 60 1 30 n_rad, nz (or n_theta), n_az, n_rad_in``
to define the grid geometry.


When I try to run mcfost, I get the error message "Exec format error. Binary file not executable". What's going on?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You're trying to run a Linux version of MCFOST on a Mac, or vice versa,
or some other incorrect operating system combination. You probably
downloaded or copied the wrong file by accident, and should just obtain
a new copy from the download site.


I get random segmentation faults
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you have increased the stacksize::

  $ export OMP_STACKSIZE=512M
  $ ulimit -s unlimited

.. note::

   There is a hard limit at 64kB for the stacksize on MacOS. This will prevent post-processing of large phantom dumps (more than 2 millions particles). The solution is either to move to linux, or to re-compile the code with gfortran (which relies less on the stack than ifort), but this will result in some speed penalties.
