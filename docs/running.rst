Running MCFOST
=================

The operation of MCFOST is controlled in two ways: a parameter file that
defines most aspects of the physical system to be simulated and the
simulation options, and a series of command line parameters that can
override the defaults or adjust various optional elements. As this is
research-grade code rather than a public project, the interface has
grown over time as new features get added on, somewhat organically and
haphazardly.

.. note:: you can type ``mcfost --help`` to get a basic help message listing
          the various command line arguments and options available.


Typical organisation of a calculation
-----------------------------------------

In case you are interested in what MCFOST is doing and not just the
results, here are the main steps in the continuum and line calculations.

Continuum calculations generally follow this scheme:

1. Setup the density structure and optical properties.
2. Perform temperature and scattering source function calculations (Monte Carlo).
3. SED and/or images calculations (ray-tracing).

These 3 steps are transparent for the user and done in a single run of
MCFOST.

NLTE line calculations generally follow this scheme:

1. Setup density structure and continuum optical properties
2. Perform temperature and radiation field calculations (Monte Carlo)
3. Setup velocity field, external radiation field and line optical properties
4. (optional) chemistry calculations via external code to calculate the abundances
5. 1+1D NLTE level population calculations (grid based, short characteristics, iterative scheme)
6. 2D or 3D NLTE level population calculations (Monte Carlo, long characteristics, iterative scheme)
7. line emission calculations (ray-tracing)

If chemical abundances are assumed to be constant, step 4 is not
required. In this case, all the calculations are performed in a single
run of MCFOST. Steps 5 and 6 are skipped if LTE is assumed (valid for
low level molecular rotational transitions).



SED Calculation
---------------

The basic usage is just::

$ mcfost <parameter_file>

That will compute the thermal structure of the disk, and its SED for the
disk structure and viewing geometry described in the named parameter
file. These files will be saved in a new subdirectory ``data_th`` under
the current directory.

If you run the code, and the desired output directory exists, it will be
moved to ``data_th_old``, and the new directory created again for the
new output. If you try this again, the command will fail. It will never
overwrite any existing files; you have to explicitly delete your old
outputs.

Image and Polarization Map Calculations
---------------------------------------

To compute an image at a certain wavelength, run::

$ mcfost <parameter_file> -img 2.0

where ``2.0`` is the desired wavelength in microns. This will produce
output in a new subdirectory ``data_2.0`` under the current directory.

Molecular Line Calculations
---------------------------

A basic line calculation is performed by using::

  $ mcfost <parameter_file> -mol

This will compute the temperature structure, then the population
level and finally the line emission map using ray-tracing.
