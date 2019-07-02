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



Temperature and SED Calculation
-------------------------------

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

mcfost will search for a previously calculated temperature structure in ``data_th``. Note that the temperature must have been calculated with the physical parameters (and command line option) for the final image to make any sense. Only parameters relative to images (number of pixels, map and pixel size, orientation) can be modified between the calculation of the temperature and of an image.

By default mcfost needs a previously computed temperature structure to make an image. At short wavelength, you can decide that thermal emission is negligible. In that case, you can run mcfost with the `-only_scatt` option, and mcfost will not try to read a temperature structure.

Molecular Line Calculations
---------------------------

A basic line calculation is performed by using::

  $ mcfost <parameter_file> -mol

This will compute the temperature structure (ie you do not need to run `mcfost <parameter_file>` in advance to compute the temperature as in the case of images), then the population level and finally the line emission map using ray-tracing.


To calculate, for example, CO J=3-2 (870 micron), for radial velocities from -10
km/s to +10 km/s in delta v of 0.05 km/s (200 channels):

::

 #Molecular RT settings
   T T T 15.                    # lpop, laccurate_pop, LTE, profile width (km.s^-1)
   0.05                         # v_turb (km.s^-1)
   1                            # nmol
   co@xpol.dat 6                # molecular data filename, level_max
   10 100                       # vmax (km.s^-1), n_speed
   T 1.e-4 abundance.fits.gz    # cst molecule abundance ?, abundance, abundance file
   T  1                         # ray tracing ?,  number of lines in ray-tracing
   2                            # transition numbers

If you need a molecular map + corresponding continuum, the recommended workflow is then::

    $ mcfost <parameter_file> -mol

which will compute the temperature, SED and molecular maps, followed by::

  $ mcfost <parameter_file> -img 1300

which will compute a continuum image at 1.3mm.
