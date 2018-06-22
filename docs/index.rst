Welcome to MCFOST's documentation
=================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. role:: math(raw)
   :format: html latex
..

About this documentation
========================

This manual for MCFOST is collectively by Marshall Perrin, Gaspard
Duchene, and Christophe Pinte. It is based in part upon an earlier
incomplete manual draft started by Christophe, and additional
independent documentation by Sebastien Perez. No guarantees implied;
this is likely to be outdated and certainly is incomplete!

Contributions of improvements are welcomed.

1. MCFOST Overview
==================

MCFOST is a Monte Carlo Radiative Transfer (plus ray-tracing) code for
simulating circumstellar disks. subsequently updated and described by
Pinte et al. (2006). However the code has been substantially further
enhanced since then, with features that are not well documented in any
one place (see Pinte et al. 2009 for a more recent paper describing the
code, though again, that is far from up-to-date). See Pinte et al.
(2008), Duchene et al. (2010) for recent applications.

MCFOST is a 3D continuum and line radiative transfer code based on the
Monte Carlo method. It is mainly designed to study the circumstellar
environment of young stellar objects. The calculations are done exactly
within the limitations of the Monte Carlo noise and machine precision,
i.e. no approximation are used in the calculations. The code has been
strongly optimized for speed. *The code is not public domain but
available on request on a collaborative basis, with explicit permissions
of the authors.* MCFOST is in constant development and this
documentation is very likely to be lagging. Please contact Christophe
Pinte for the latest updates.

The algorithms used for the continuum transfer are described in detail
in Pinte et al. (2006, 2009), following an earlier description by Menard
(1989). In short, the code computes the temperature and scattering
source function everywhere in the disk via a Monte Carlo method: photon
packets are propagated stochastically through the model volume following
the equations of radiative transfer, and information on their properties
is retained along their path. The radiation field, and quantities
derived from it (for instance temperature, radiation pressure, etc) are
obtained by averaging this “Monte Carlo” information. Observables
quantities (SEDs and images) are then obtained via a ray-tracing method,
which calculates the output intensities by integrating formally the
source function estimated by the Monte Carlo calculations. Full
calculations of the polarization are included using the Stokes
formalism.

MCFOST also includes a non-LTE line transfer module. The adopted schemes
are not yet described in the literature, but the code uses an improved
version of the algorithm presented in Hogerheijde & van der Tak (2000).
NLTE level population are obtained via iterations between Monte Carlo
radiative transfer calculations and statistical equilibrium. To speed up
calculations, an initial guess is computed using a 1+1D scheme with
short characteristic before the full 2D or 3D calculation. Output
spectra and channel maps are calculated via a ray-tracing procedure.

Associated with MCFOST, but separate software, are tools for
visualization, grid calculations and model fitting tools for SEDs,
images and line emission. These tools allows one to draw robust
constraints on the derived parameters.

1.2 A Caution
-------------

*Warning: Use at your own risk!!! The author does not take any
responsibility for the use (or misuse of the code). There might be
bugs.* *This documentation is not complete, and may in places be out of
date with respect to the latest changes in the code. *

1.3 Disclaimers
---------------

MCFOST is available on a **collaborative basis**.

Using MCFOST implies that you agree to :

-  offer us (C. Pinte, F. Ménard, G. Duchêne) co-author right on any
       resulting publication.

-  NOT distribute MCFOST without our explicit agreement.

-  contact us if you initiate a new scientific project with MCFOST.

Please also acknowledge funding from the European Commission's seventh
Framework Program

(contract PERG06-GA-2009-256513) and from Agence Nationale pour la
Recherche (ANR) of France under contract ANR-2010-JCJC-0504-01.

The IDL code MCRE (MCFOST Results Explorer) is also available on a
collaborative basis. Using MCRE implies to offer M. Perrin co-author
right. However its use is deprecated. Instead users are encouraged to
use the Python mcfost package, available from
`*https://github.com/cpinte/mcfost-python* <https://github.com/cpinte/mcfost-python>`__,
by M. Perrin, C. Pinte, and S. Wolff.

1.4 What can be modelled using MCFOST?
--------------------------------------

MCFOST is primarily designed to study protoplanetary disks. The code can
reproduce most of the observations of disks:

-  SEDs

-  scattered light images

-  IR and mm visibilities

-  atomic and molecular line maps

The Monte Carlo method being generic, any complex structure can be
handled by MCFOST and its use can be extended to other astrophysical
objects. For instance, tests have been performed on infalling envelopes
and AGB stars.

1.5 Typical organisation of a calculation
-----------------------------------------

In case you are interested in what MCFOST is doing and not just the
results, here are the main steps in the continuum and line calculations.

Continuum calculations generally follow this scheme:

1. Setup the density structure and optical properties.

2. Perform temperature and scattering source function calculations
       (Monte Carlo)

3. SED and/or images calculations (ray-tracing)

These 3 steps are transparent for the user and done in a single run of
MCFOST.

NLTE line calculations generally follow this scheme:

1. Setup density structure and optical properties

2. Perform temperature and radiation field calculations (Monte Carlo)

3. Setup velocity field, external radiation field and line optical
       properties

4. (optional) chemistry calculations via external code to calculate the
       chemical abundances

5. 1+1D NLTE level population calculations (grid based, short
       characteristics, iterative scheme)

6. 2D or 3D NLTE level population calculations (Monte Carlo, long
       characteristics, iterative scheme)

7. line emission calculations (ray-tracing)

If chemical abundances are assumed to be constant, step 4 is not
required. In this case, all the calculations are performed in a single
run of MCFOST. Steps 5 and 6 are skipped if LTE is assumed (valid for
low level molecular rotational transitions).

1.6 Mailing List
----------------

| If you decide to use MCFOST, you should sign up to the mailing list in
  order to be informed when new versions are available. Please send an
  email to `*sympa@ujf-grenoble.fr* <mailto:sympa@ujf-grenoble.fr>`__
  with header :
| subscribe mcfost <First Name> <Last Name>

First name and Last Name are optional.

2. Requirements and Installation
================================

MCFOST is written in Fortran 90 and parallelized with the Open MP
protocol, i.e. it can use several processors and/or cores on a single
machine. For simplicity, we provide pre-compiled binaries which do not
have any dependencies. Source code can also be provided if needed.

The following environment is required for the binaries:

-  a 64-bits Unix system, including Linux or MacOS X,

-  any data analysis software capable of handling FITS files.

MCFOST is also available for Xeon Phi but performance is not optimal
yet. Stay tuned.

For visualization and fitting tools, 3 packages are available so far for
MCFOST:

-  an IDL GUI interface (lead author: Perrin)

-  a more simple, command-line based IDL interface (lead author:
       Duchêne)

-  a command-line based, Yorick interface (lead author: Pinte).

-  [new] a Python package for display and manipulation of models. (Work
       in progress; Perrin & Pinte)

Only the analysis of the results are done with IDL or Yorick or Python.
They are not required for the use of MCFOST, but make its use easier,
especially when dealing with large numbers of models.

MCFOST uses a database of stellar spectra, optical properties and atomic
and molecular data. These files are generally put in a directory named
mcfost/utils, although any name can be used. The environment variable
MCFOST\_UTILS must be set to the path name of this directory.

An additional environment variable MY\_MCFOST\_UTILS can be defined by
the user to add his own data files. This has an advantage to ensure that
no personal data files will be overwritten during an update of the utils
directory. See Section 6 below for more details on how this is used.

**Installation Procedure:**

1. Copy the file mcfost to somewhere in your shell $PATH.

2. chmod 755 mcfost; rehash

3. You should now be able to run mcfost --help

4. | Set the environment variable MCFOST\_UTILS to point to a directory
         where you want mcfost to store its data files. E.g. edit your
         shell startup files to include either
       | setenv MCFOST\_UTILS /path/you/put/the/files [tcsh or csh
         shell], or
       | export MCFOST\_UTILS=/path/you/put/the/files [bash or bash-like
         shell]

5. run mcfost -setup to download mcfost data files, current reference
       parameter files and manual.

3. Running MCFOST
=================

The operation of MCFOST is controlled in two ways: a parameter file that
defines most aspects of the physical system to be simulated and the
simulation options, and a series of command line parameters that can
override the defaults or adjust various optional elements. As this is
research-grade code rather than a public project, the interface has
grown over time as new features get added on, somewhat organically and
haphazardly.

Note that you can type

mcfost --help

to get a basic help message listing the various command line arguments
and options available.

SED Calculation
---------------

The basic usage is just

mcfost “parameter\_file”

That will compute the thermal structure of the disk, and its SED for the
disk structure and viewing geometry described in the named parameter
file. These files will be saved in a new subdirectory “data\_th” under
the current directory. For historical reasons, this command will produce
a purely Monte Carlo-based SED, which are generally noisy in the near-
and far-IR unless a large number of photons are used. To compute a ray
tracing SED instead, type the following command:

mcfost parameter\_file -rt

This is particularly useful for edge-on optically thick disks. For close
to pole-on disks, the overhead may be too high and the Monte Carlo alone
approach is probably preferable unless very high signal to noise is
required.

If you run the code, and the desired output directory exists, it will be
moved to e.g. data\_th\_old, and the new directory created again for the
new output. If you try this again, the command will fail. It will never
overwrite any existing files; you have to explicitly delete your old
outputs.

Note: When using the -rt option, polarisation will be turned off
automatically in SED mode. If you really need SED+RT+pola you can have
it using the image ray-tracing mode, ie using -rt2 instead of -rt but
that may take a while. SEDs can also be computed with polarization in
Monte Carlo.

Image and Polarization Map Calculations
---------------------------------------

To compute an image at a certain wavelength, run

mcfost parameter\_file -img 2.0 -rt

where ‘2.0’ is the desired wavelength in microns. This will produce
output in a new subdirectory “data\_2.0” under the current directory.
The ray tracing mode should always be used for image calculations, and
will likely become the default soon.

Molecular Line Calculations
---------------------------

A basic line calculation is performed by using: mcfost <parameter\_file>
-mol. This will compute the temperature structure, then the population
level and finally the line emission map using ray-tracing.

Parallelization Options
-----------------------

By default, MCFOST will parallelize itself across all available CPU
cores. If you want to restrict it to a subset, you can specify the
number of cores to use with the environment variable OMP\_NUM\_THREADS.
In some tests on a circa 2010-era Mac Pro, performance degrades above 4
cores in ray-tracing mode, so we recommend

setenv OMP\_NUM\_THREADS 4

If you wish to disable parallelization entirely, you can use

setenv OMP\_NUM\_THREADS 1

We recommend carrying out your own speed tests to see if there is an
optimal setting on your specific computers.

Here are results from testing on a 2014 Mac Pro (3 GHz, 8 core Intel
Xeon E5 with 32 GB DDR3 RAM) by Marshall Perrin. This is for calculating
the SED for one particular model file (chosen arbitrarily) from an MCMC
chain prepared by Schuyler Wolff. The scaling is not quite 1/N, but it’s
pretty good up to 8 threads, which is the # of true CPU cores this
computer has.
`*Hyperthreading* <http://en.wikipedia.org/wiki/Hyper-threading>`__
results in the computer appearing to have 16 virtual cores, but the
performance gain from trying to use these all is marginal.

|image0|

+----------------+---------------------+--------------------------+
| # of threads   | CPU time used [s]   | Total elapsed time [s]   |
+================+=====================+==========================+
| 1              | 141                 | 141                      |
+----------------+---------------------+--------------------------+
| 2              | 159                 | 79                       |
+----------------+---------------------+--------------------------+
| 4              | 160                 | 40                       |
+----------------+---------------------+--------------------------+
| 8              | 186                 | 23                       |
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
| 8              | 30                  | 3                        |
+----------------+---------------------+--------------------------+
| 12             | 40                  | 3                        |
+----------------+---------------------+--------------------------+
| 16             | 47                  | 3                        |
+----------------+---------------------+--------------------------+

4. Outputs from MCFOST
======================

Outputs from SED mode calculations
----------------------------------

A calculation in SED mode will output the following files:

-  A copy of the input parameter file. A new line will be added at the
       end of this file, recording any command line options given in the
       call to MCFOST.

-  **sed\_mc.fits.gz**

-  | **Temperature.fits.gz** (and Temperature\_DA.fits.gz if diffusion
         approximation is used)
       | + additional temperature files of out-of-equilibrium grains are
         present

-  **sed\_rt.fits.gz** [only if the -rt command line option has been
       used]

Units for all SED files are W.m\ :sup:`-2`. The temperature map is in K,
with values of 1K for regions where the temperature was not computed and
is later estimated via the diffusion approximation (this occurs in
extremely deep regions of disks).

Output SED files typically are 4 dimensional. The dimensions are
[wavelength, inclination, azimuth, stokes]. Dimensions may be of length
1; for instance typically the azimuth axis has only 1 value for the case
of axisymmetric models.

An additional .sed\_th.fits.gz is also created in the data\_th
directory. It records the SED produced by the Bjorkman & Wood algorithm
during the calculation of the temperature. This file is generally noisy
at wavelengths where the flux is low. It is mainly used for testing
purposes. In particular, each re-emission takes place at the exact same
location of the absorption, while in sed\_mc.fits.gz and
sed\_rt.fits.gz, the emissivity is assumed to be constant within a cell.
The comparison of .sed\_th.fits.gz and sed\_mc.fits.gz can be useful to
detect inadequate spatial resolution (ie cells that are too big to
accurately sample the local emissivity)

Should I use Monte Carlo (MC) or ray-traced (RT) SED calculations?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sed\_rt offers several advantages : it correspond to exact inclinations
(instead of bins for sed\_mc) and is converges almost equally fast at
all inclinations. The ray-tracing as the default mode for SEDs as there
is an extra overhead time (and it uses extra memory) to compute the
scattered specific intensities and it is usually faster to compute the
SED in MC mode for non edge-on disks. Parallelization is also limited in
ray-tracing mode, due to high memory sharing between threads.

In short, if the disk inclination is low, MC SED is generally faster
unless you need very high signal-to-noise. For high inclinations, RT SED
is always preferable (in particular because the SED is very sensitive to
small changes of inclinations). It is also preferable to use the RT SED
if you wish to separate the various contributions.

Outputs from Image calculations
-------------------------------

-  A copy of the input parameter file. A new line will be added at the
       end of this file, recording any command line options given in the
       call to MCFOST.

-  **RT.fits.gz** Raytraced output image (if -rt option is used)

-  **MC.fits.gz** Monte-Carlo output image (otherwise)

Units for the images are W.m\ :sup:`-2`.pixel\ :sup:`-1`. These images
have at least 3 dimensions, two spatial dimensions and one for
inclination. Up to 5 dimensions may be present, in which case the 4th
dimension is for the azimuthal viewing angle in case of a
non-axisymmetric disk, and 5th dimension is the Stokes parameters for
polarized imaging (I,Q,U and V) and different contributions if the
contributions are separated (star alone, scattered light from the star,
thermal emission, scattered light from the dust thermal emission).

If the Stokes parameter option is turned on, the output is Stokes
polarized intensities [I, Q, U, V] and *not* fractional polarization [I,
Q/I, U/I, V/I]. That is to say, all the polarizations are written with
the same dimensionality as intensities and units of
W.m\ :sup:`-2`.pixel\ :sup:`-1`.

Other kinds of output [optional]
--------------------------------

**Disk Density Files**

You can request an output of the computed disk density. This is done
with the -disk\_struct (or -output\_density\_grid) command line option.
This will output files in the *current* directory

-  gas\_density.fits.gz [map of the gas density - in g/cm\ :sup:`3`]

-  dust\_mass\_density.fits.gz [map of the dust density - in
       g/cm\ :sup:`3`]

-  dust\_particule\_density.fits.gz [map of the dust density - in
       part/cm\ :sup:`3`]

-  grain\_sizes.fits.gz [grain size bins - in μm]

-  grid.fits.gz [map of the radii and height of the individual grid
       cells - in AU]

-  volume.fits.gz [map of the volume of the individual grid cells - in
       AU\ :sup:`3`]

The various density files will have 3d arrays of dimension [n\_rad, n\_z
or n\_theta, n\_az], where those parameters are as defined below in the
parameter file geometry. The grid.fits file will have those same
dimensions for its first 3 axes, then the 4th axis will have dimension
2. The first plane in that axis will give the radius for each grid cell,
the second will give the height.

**Dust Property Files**

You can request an output of the computed dust propertis. This is done
with the -dust\_prop command line option. This will generate a
data\_dust directory that will contain the following output files (along
with a copy of the parameter file):

-  albedo.fits.gz [table of the dust albedo]

-  g.fits.gz [table of the scattering phase function asymmetry
       parameter]

-  kappa.fits.gz [table of dust opacities - in cm\ :sup:`2`/g]

-  kappa\_grain.fits.gz [table of dust opacities as a function of grain
       size - in cm\ :sup:`2`/g]

-  lambda.fits.gz [table of wavelength at which the dust properties are
       computed - in micron]

-  phase\_function.fits.gz [table of scattering phase function from 0 to
       181 degrees, versus wavelength]

-  polarizability.fits.gz [table of scattering polarizability from 0 to
       181 degrees, versus wavelength]

The lambda.fits.gz file also provides the wavelength array for an SED
computation when lsed\_complete is set to T (see below).

If the dust properties are only needed at a specific wavelength (say 2.2
micron), the user can instead use the “-op 2.2” command line option. In
this case, only the initialization sequence is performed, and the albedo
and asymmetry parameter are given on screen at the end (as well as the
total optical depth along the equatorial plane and along the viewing
“inclination of interest” indicated in the parameter file). Note that
these quantities are always provided in the normal run of MCFOST (for
both SED and image calculations); this option merely stops the
calculation after this stage. Note that a data\_2.2 directory is created
containing only a copy of the parameter file.

**other...**

5. Parameter File Contents
==========================

The MCFOST input parameter file is a fixed format text file, that has
gained more features over time as the code has developed. Each line
consists of one or more values separated by whitespace. Anything after
the last value is interpreted as a comment. The number of lines in the
file and the number of items on each line must be exactly right or
MCFOST will not run. The parsing code is flexible enough to accept
different amounts of white space, though. Extra empty lines can also be
added.

The file format changes slightly with different versions of MCFOST. The
first line of the parameter file gives a version number that specifies
how that file should be interpreted. MCFOST tries to remain compatible
with old versions of the parameter files.

The following description of the many parameters is very very
incomplete... See Pinte et al. 2006 for most of the relevant equations.

Photon count setup
------------------

**nbr\_photons\_eq\_th**: Number of “photons” launched when computing
the temperature structure of the disk (typically 100 000, unless large
noise is seen in temperature map)

**nbr\_photons\_lambda**: Number of “photons” *per wavelength* that must
reach the observer when computing the output SED, for the inclination
bin (“captor”) with highest inclination. Typically 10 000, unless
anomalously large noise is seen in the SED. With 10 000 packets
received, the MC noise will be ≈ 1 %.

**nbr\_photons\_image**: Number of “photons” launched when computing a
disk image. This number must be adjusted depending on the spatial and
inclination resolution. Low values, typically 100 000, are enough when
using the ray-tracing mode.

NB: if the ray tracing option is used (-rt) then the Monte Carlo process
is followed by a ray tracing calculation so the number above are
generally adequate. For a pure Monte Carlo calculations, much larger
values are generally needed.

Wavelength setup
----------------

**nlambda, lambda\_min, lambda\_max:** Number of wavelength bins, and
minimum and maximum wavelengths in units of microns. Note that
lambda\_min and lambad\_max are the extreme wavelengths at the edges of
the bins; the central wavelengths of the bins will be a little bit
inside those extreme values. If you need to force specific wavelengths,
use the lsed\_complete flag.

**ltemp:** do you want to compute the temperature map? Most of the time
it's T (unless you have already previously computed the temperature but
not the SED for some reason)

**lsed:** do you want to compute an SED? Most of the time it's T (unless
you care only about the temperature structure of the disk)

**lsed\_complete**: are you happy with the default evenly log-spaced
sampling in wavelength (T) or do you want to use the wavelength file
indicated in the following line (F)

**wavelength file**: a text file containing the wavelengths at which the
SED is to be sampled. One entry per line, all entries in microns. The
file should be located in $MCFOST\_UTILS/Lambda/ or in the local
directory where MCFOST is executed, and a copy will be placed in the
data\_th directory

**separate contributions:** Save output file that tracks the various
different contributions (star vs disk, scattered vs direct). Applicable
to both SED and image calculations, works in RT mode. In MC SED mode,
this setting is *always* on, as it is almost free to track this.

**output stokes parameters**: Keep track of the Stokes parameters for
each photon. Output images will include polarization axis. Only
applicable to image calculations, works in RT mode.

Grid geometry
-------------

| **grid\_type**: selecting a cylindrical, spherical grid geometry or a
  Voronoi tesselation of a set of pre-defined points. The spherical grid
  is mandatory if an envelope is present; for disks, the cylindrical
  grid is more appropriate.
| The Voronoi tesselation is particularly suited to use the results of
  hydrodynamical simulations, in particular SPH, where 1 SPH particle
  can be directly matched to 1 cell in MCFOST.

**n\_rad**: Number of grid cells along the radial direction

**nz (n\_theta)**: Number of grid cell along the vertical direction
(from the midplane up, for a cylindrical grid), or the “latitude”
direction (for a spherical grid)

**n\_az**: Number of grid cells along the azimuthal direction (set to 1
for an axisymmetric structure)

**n\_rad\_in**: Number of sub-cells in which to split the first cell
along the radial axis. This should be in the 10-30 range for an
optically thick disk so as to deal accurately with the highest density
regions of the disk. If the temperature map shows T=1K in the second
pixel from the disk inner radius, you need to increase n\_rad\_in

Images
------

**grid (nx, ny)**: Number of pixels for the output images (applies to
both Monte Carlo and RayTracing images). It must be odd so that the star
is exactly in the central pixel.

**map\_size**: physical size of the produced maps (in AU). Applied to
the largest axis if nx ≠ ny

**zoom factor**: can only be accessed as an option from version 2.18.
Multiplicative factor to obtain a zoomed in image centered on the star
instead of an image of the whole extent of the computing grid. Can
larger (zoom in) or smaller (zoom out) than 1. The effective map size is
map\_size / zoom

***MC\_n\_incl**: Number of inclination “bins” for the Monte Carlo
images and SEDs; the whole [0,90] range is sampled and divided evenly in
cos(i)*

***MC\_n\_phi**: Number of azimuthal viewing angles (set to 1 for an
axisymmetric disk)*

*Note: these MC parameters are no longer included in version 3.0+
parameter files.*

**RT\_imin**: Minimum inclination to compute the Ray Tracing images
(does not need to be 0\ :sup:`o`). 0\ :sup:`o` is pole-on and
90\ :sup:`o` is edge-on.

**RT\_imax**: Maximum inclination to compute the Ray Tracing images
(does not need to be 90\ :sup:`o`)

**RT\_n\_incl**: Number of inclinations for which Ray Tracing images
should be computed in the [imin,imax] range, evenly spaced in cos(i)

**RT\_centered**: Set to F to force images to be computed exactly at
imin and imax; if set to T, the sampled inclination values will be
computed at the midway point of the cos(i) bins

**RT\_az\_min**: Minimum azimuthal angle to compute the Ray Tracing
images (does not need to be 0\ :sup:`o`).

For a pole-on model, 0\ :sup:`o` means that the x-axis of the model is
towards the observer and 90\ :sup:`o` means the y-axis is towards the
observer.

**RT\_az\_max**: Maximum azimuthal angle to compute the Ray Tracing
images.

**RT\_n\_az**: Number of azimuthal angles for which Ray Tracing images
should be computed in the [RT\_az\_min,RT\_az\_max] range, linearly
spaced.

**distance**: Distance to the object in pc

**disk PA**: Position angle of the semi-minor axis of the disk, measured
counter-clockwise. If disk PA is not set to 0\ :sup:`o`, some of the
image “symmetries” (see below) will automatically be set to F

Scattering Method
-----------------

**scattering method**: compute the average dust properties per grain
size (value: 1), per grid cell (value: 2), or the most appropriate for
the ongoing calculation (value: 0, preferred in virtually all cases).
Value 2 can not be used if the Stokes parameters are required.

**Mie/hg**: choice of scattering phase function; Mie theory (value: 1)
is strongly preferred over the Henyey-Greenstein parametric description
(value: 2) since it is physically grounded and allows calculations in
full-Stokes mode. If the HG phase function is selected, Mie theory will
first be used to compute the effective *g* value, which will then be
used to randomly select scattering angles.

    *Optical properties calculations are much faster when the option is
    set to 2, reducing significantly MCFOST’s initialization time. It
    may be useful to compute quickly temperature structures and SEDs
    when very large grains are present and/or a large number of
    wavelengths are used (the overhead due to Mie theory is almost
    always negligeable for monochromatic images).*

Symmetries
----------

**image symmetry**: Is the image left/right symmetric?

**central symmetry**: Is the image point symmetric?

**axial symmetry**: Is the disk structure axisymmetric?

In most cases, all three symmetries should be set to T. If you use a
non-axisymmetric disk structure, or if there is more than one star
illuminating the disk, then they should be set to F. If a disk PA
different than 0\ :sup:`o` has been set, then first two symmetries
should be set to F and the last one to T.

Disk physics
------------

**dust\_settling:** no settling (0), parametric (1), following
Dubrulle’s (2) or Fromang’s prescription (3).

    *Parametric settling is independent of radius, it is simply a
    scaling of the scale height as a function of the grain size.*

    *Dubrulle’s and Fromang’s settling assume a diffusion equation
    depending on the α viscosity. Dubrulle’s vertical profile remains
    Gaussian, while Fromang’s presciption is more realistic and departs
    from the Gaussian profile at high altitude. (see Dubrulle et al 1995
    and Fromang et al 2009).*

**exp\_strat**: power law describing the parametric settling, where H(a)
decreases as a\ :sup:`-exp\_strat` for a > a\_strat

**a\_strat**: minimum size [microns] for grains affected by vertical
settling (for parametric settling)

**dust radial\_migration**: simple prescription for radial migration of
dust grains (TBW)

**hydrostatic equilibrium**: (*work in progress*) compute the
hydrostatic equilibrium assuming Tgas = Tdust

**sublimate\_dust**: this option will iteratively remove the dust if the
temperature reaches the dust sublimation temperature

**viscous heating**: (*work in progress*) includes additional heating
source due to viscous accretion (note that it does not account for the
accretion shock on the star)

**viscosity**: alpha parameter describing the strength of the viscosity
(used for settling (mode 2 and 3), and viscous heating).

Number of Zones
---------------

The following sections density structure--grain properties will be
repeated n times, depending on the number of zones set here. This lets
you describe a complex multi-component system.

MCFOST can use as many zones as required but the memory usage and cpu
time for the initialization will increase with the number of zones. If
you use a large number of zones you might also need to ensure that the
resolution of the spatial grid is high enough.

Density structure
-----------------

| **zone type**: disk (value: 1), disk with outer tapered-edge (value:
  2), spherically symmetric envelope (value: 3), debris disk (value : 4)
  or an azimuthally asymetric wall (value : 5).
| If at least one of the zones is described as an envelope, the
  computing grid must be spherical

**disk dust mass**: in units of M\ :sub:`sun`

**gas-to-dust ratio**: quantity only used in molecular emission
calculations

**scale height H0:** value of the disk scale height (technically, sigma
of the Gaussian vertical density profile) at the reference radius, in AU

**reference radius R0**: in AU

**γ\ :sub:`vert` vertical profile exponent :** exponent of vertical
density profile (only relevant for type 4, debris disk)

**Rin**: inner radius of disk or envelope in AU (assuming sublimation
calculation is not enabled)

**Rout or Rc for tapered-edge**: outer radius of disk or envelope in AU

**edge**: make nonzero for gradual rather than abrupt falloff inside
rin/outside rout.

    *This is to set a smooth decline in surface density at the inner and
    outer edges of the disk. Inside of r\_in (and outside of r\_out),
    the density drops following a Gaussian whose sigma is the "edge"
    parameter.*

**β :** flaring exponent = power law index of the H(r) scale height
function, typically in the [1.0-1.25] range

**p1 and p2 : surface density exponent, or -gamma for tapered-edge:**
power law index of the surface density profile (generally <0).
**-gamma\_exp** for tapered-edge. Defines p\_in and p\_out in case of a
debris disk.

The disk density structures are defined as :

1. :math:`\Sigma(r)\ \alpha{{\text{\ r}{}^{p}}^{1}}^{}`

2. :math:`\ \Sigma(r)\ \alpha{{\text{\ r}{}^{p1\ }}^{}}^{}exp( - (\frac{\text{\ r}}{\text{Rc}}){}^{2 + p2}\ )\ \alpha{{\text{\ r}{}^{- \gamma\ \ }}^{}}^{}exp( - (\frac{\text{\ r}}{\text{Rc}}){}^{2 - \gamma_{\exp}}\ )`

3. :math:`\rho(r)\ \alpha\ \text{r\ }^{p1}`

4. :math:`\varrho(r,\ z)\ \alpha{{\ \ (\ (\frac{\text{\ r}}{\text{Rc}}){}^{- 2p_{\text{in}}}\  + \ (\frac{\text{\ r}}{\text{Rc}}){}^{- 2p_{\text{out}}})^{- 1/2}\ }^{}}^{}\  \times \ exp( - (\frac{\ \left| z \right|}{h(r)}\ )^{\gamma_{\text{vert}}})`
       see Augereau et al, 1999, A&A, 348, 557. p\ :sub:`in` = p1 > 0
       and p\ :sub:`out` = p2 < 0

If Rout is set to 0, it is automatically set Rout to 8 Rc in the case of
a disk with tapered-edge

For types 1,2, 4 and 5, the local scale height is defined as h(r) =
H\ :sub:`0` x (r/R:sub:`0`)^β

*Cavity*
--------

***cavity**: Is there a cavity?*

***height**: height at the reference radius above which density is set
to zero, in AU*

***reference radius**: in AU*

***flaring exponent**: power law index of the z\ :sub:`max`\ (r)
function describing the cavity wall*

*Note: all this is only relevant if density structure is an envelope. *

*Note 2: This section removed in version 3.0+ parameter files. *

Grain properties
----------------

NB: there should be as many blocks containing the following parameters
as there are zones in the disk. The code will crash otherwise.

**N\_species**: Number of dust populations present in the disk zone; if
N\_species > 1, the dust grains of different species are assumed to be
physically disjoint, but distributed in the same manner through the
disk. *All the following lines in the block must be duplicated
N\_species times if N\_species > 1.*

**Grain\_type:** spherical grains (Mie) or distribution of hollow
spheres (DHS)

**N\_components**: Number of materials that make up a given “specie”;
these materials are assumed to be physically joint within each dust
grain. *The line with the optical indices and volume fraction must be
duplicated N\_components times if N\_components > 1.*

**mixing rule:** are the components randomly mixed within the volume of
a grain (value: 1, effective mixing theory following Bruggeman rule) or
is the second component forming a coating on top of the first one
(value: 2). The effective optical index of the new "mixed" grain is
computed before any Mie theory computation. Does not apply is
N\_components = 1.

*Coating can only be used with 2 components (the 1st one is the core,
the 2nd one the shell).*

**porosity**: porosity of the dust grains (in the [0,1] range, 0 for
compact grains, near 1 for porous ones)

**mass fraction**: fraction of the mass contained in this specie (the
sum of the N\_species “mass fractions” should be equal to 1, MCFOST will
renormalize the values so that the sum is 1)

**DHS\_V\ :sub:`max`** maximum void fraction for DHS calculation

**optical indices**: file containing the optical index of the material
as a function of wavelength (files must be located in
$MCFOST\_UTILS/Dust/)

**volume\_fraction**: fraction of the volume of a grain contained in
this component (the sum of the N\_components “volume fractions” should
be equal to 1, MCFOST will renormalize the values so that the sum is 1)

**heating method**: indicated whether radiative equilibrium and local
thermal equilibrium are assumed: for an optically thick disk, both
should be true (value: 1); for an optically thin disk, only the RE is
assumed (value: 2, will yield a temperature map that has a third
dimension spanning the grain size distribution); for out-of-equilibrium
grains, none of them is true (value: 3, typical for PAH grains)

**amin**: minimum grain size, in microns

**amax**: maximum grain size, in microns

**aexp**: power law index of the grain size distribution dN(a)/da

**n\_grains**: number of grain size bins to the sample [evenly in
log(a)] the grain size distribution; dust properties are only computed
for these grain sizes and subsequent interpolations are used whenever
necessary. Typical value is in the 50-100 range, less in case of
multiple dust specie/disk zones to limit RAM requirement (and
computation time of the Mie theory).

Molecular RT settings
---------------------

**lpop:** do you wish to compute the level populations (this might not
be the case if you use populations from an external code, ProDiMo for
instance)

**lpop\_accurate:** if the variable is set to false, mcfost will just
perform a 1+1D (or 1+1+1D) line transfer. If set to true, the result of
the 1+1D line transfer will be used as a starting point for the full 2D
or 3D line transfer.

**LTE:** assume LTE level populations

**profile width:** internal line width used for the line transfer
calculation. Bascically, it means that cells with relative projected
velocities that exceed this value will not see each other during the
transfer. The value needs to be larger than the local line width.

**molecular data file:** LAMBDA data file used for the line transfer

**level max:** maximum level up to which the line transfer will be
performed. Level above “level max” will not be populated

**vmax [km/s]:** maximum velocity is the produced channel maps

**n\_speed:** number of velocity points between 0 and vmax

**cst\_abundance:** do you wish to use a constant molecular abundance
over the disk. If true, use the provided abundance, if false read the
abundance from the following fits file. The resolution of the fits file
must be the same as the grid used by mcfost.

**ray-tracing:** produce or not a ray-traced data cube of the molecule

**number of lines:** number of transition ray-traced, the indices given
in the next lines correspond to the transition indices in the LAMBDA
files. For instance the J=1-0 transition is usually #0

Star properties
---------------

**num\_stars:** number of stars illuminating the disk. *If num\_stars >
1, the following group of two lines must be duplicated num\_stars times*

**Temp**: effective temperature of the star (only used to compute the
stellar luminosity), in K

**radius**: stellar radius (used to compute the stellar luminosity and
set the spatial origin of the photon packets: the star is assumed to be
a uniformly radiating sphere), in R\ :sub:`sun`

**M**: stellar mass (only used for molecular line calculations,
hudrostatic equilibrium and viscous heating via accretion, in
M\ :sub:`sun`

**x,y,z,**: position of the star, in AU. If this is not (0,0,0), the
image symmetries must be set to F

**is a blackbody**: is the stellar emission approximated by a blackbody
of temperature Temp? If not, the stellar spectrum indicated in the
following line is used instead

**fUV, slope\_fUV:** photospheric UV excess and its slope in
*F\ :sub:`ν`*

    *The basic underlying assumption is that the total emission from the
    star can be approximated by whatever stellar model you use plus a
    single power law UV excess. That UV excess is characterized by two
    parameters: fUV is the scaling factor while slope\_fUV is the power
    law index. The definition of fUV can be found in Woitke et al.
    `2010 <http://adsabs.harvard.edu/abs/2010MNRAS.405L..26W>`__, and is
    the ratio of the UV flux between 92 and 250nm to the total
    photospheric luminosity. If fUV=0, then there's no excess on top of
    the phtosphere. The UV is only important for chemistry purposes in
    practice, unless you use some crazy high value of fUV.*

**T\ :sub:`eff` vs stellar atmosphere model:**

    *The effective temperature is only used to compute the total stellar
    luminosity (via 4π R\ :sub:`star`\ :sup:`2`* x
    *σT\ :sub:`eff`\ :sup:`4`) and it doesn't need to match the
    atmosphere model. In practice, the atmosphere model shoud be chosen
    with the effective temperature closest to the target's, use the same
    Teff and adjust Rstar to get the right Lstar. However, one could use
    more precise numbers for both Teff and Rstar (when available in the
    literature, for instance). Technically, so long as you get the right
    luminosity and right atmospheric model, it doesn't matter if you use
    a "wrong" Teff.*

Stellar spot
------------

A spot can be added to photosphere using the option:

mcfost <para> -spot <T\_spot> <spot\_surface\_fraction> <i\_spot>
<phi\_spot>

The spot surface fraction is defined between 0 and 1. The spot
inclination is defined between 0 and 180 degrees: 0 degree for a spot on
the pole and 90 degrees for a spot on the equator.

The option currently only works if the photosphere has no extra UV and
is only available in MC mode so far.

If i is defined between 0 and 90 degrees, the spot is facing the
observer when using the central azimuthal bin.

6. Input Data Files
===================

All the data used by mcfost is generally read from the directory
$MCFOST\_UTILS. An additional environment variable MY\_MCFOST\_UTILS can
be defined (giving the path to another directory) by the user to add
her/his own data file. This has to advantage to ensure that no personal
data file will be overwritten during an update.

When looking for data files, mcfost will first look in the local working
directory, then in MY\_MCFOST\_UTILS if it is defined and finally in
MY\_MCFOST\_UTILS. If the data file used is in the local directory, it
will be copied into the output directory as well.

The data directory also has a version number which correponds to the
minimum MCFOST version which can use this data set (ie MCFOST version ≥
MCFOST\_UTILS version). MCFOST will ask the user to update the data
directory if needed (note this also means that you won’t be able to use
an old MCFOST version if the data directory is updated).

Dust Optical Constants
----------------------

Filenames that are underlined in the following list indicate
well-identified dust species (i.e., the origin of the file can be
tracked in the literature with reasonable certainty). Comments at the
top the files must start with a #. The first line contains the material
density in g/cm\ :sup:`3` and the sublimation temperature. It is
followed by an empty line and then 3 columns containing the wavelengths
(in microns) and the real and imaginary parts of the optical indices.
Wavelengths can be in ascending or descending order but must be ordered.

-  **CO2\_50K-reduit.dat:** CO2 ice at 50K (origin?)

-  **Draine\_Si.dat: “**\ Astronomical silicates” from Bruce Draine
       *(very similar though not exactly the same as those shown in Fig
       3 of Draine & Lee 1984, ApJ, 285, 89 - use dlsi\_opct.dat
       instead);* file **Draine\_Si\_sUV.dat** is the version which has
       been smoothed in the UV.

-  **Graphite\_para.dat:** graphite with a “parallel” orientation (from
       Draine & Lee 1984, ApJ, 285, 89?)

-  **Graphite\_perp.dat:** graphite with a “perpendicular” orientation
       (from Draine & Lee 1984, ApJ, 285, 89?)

-  ***MW89.dat*:** Highly porous composite interstellar grains. From
       `*Mathis and Whiffen
       1989* <http://adsabs.harvard.edu/cgi-bin/bib_query?1989ApJ...341..808M>`__.

-  ***Olivine\_modif.dat*:** Amorphous MgFeSiO4 olivine. From Dorschner,
       Begemann, et al. 1995. A&A 300 503, see their Fig 5). File
       **Olivine.dat** contains the same optical indices but has no
       information for lambda > 500 micron, forcing the code to
       extrapolate. File **Olivine\_no\_SI.dat** contains the same
       optical indices as **Olivine\_modif.dat**, except for an
       artificial linear interpolation across the 10 and 20 micron
       features to suppress them (no obvious reason why that would be
       useful...).

-  ***ac\_opct.dat*:** (amorphous) carbonaceous dust, Rouleau & Martin
       1991, ApJ, 377, 526

-  ***crstsi\_opct.dat*:** crystalline silicate dust (olivine), optical
       constants derived by Li, A., & Draine, B. T. 2001, ApJ, 550,
       L213, no corresponding optical constants figure in the literature

-  ***dlsi\_opct.dat*:** amorphous silicate dust ("astronomical
       silicates") from Draine, B. T., & Lee, H. M. 1984, ApJ, 285, 89

-  ***ice\_opct.dat*:** H2O-dominated "dirty" ice (mostly crystalline)
       compiled from several sources described in Li, A., & Greenberg,
       J. M. 1998, A&A, 331, 291. no corresponding optical constants
       figure in the literature

-  **indices\_acar.dat:** “amorphous carbon” optical indices file used
       by Sebastian Wolf with MC3D (origin?). *Judging from the shape of
       the plots, this material also contain some silicates, though.*

-  **indices\_acar2.dat:** “amorphous carbon” with high porosity (~85%)
       optical indices file used by Sebastian Wolf with MC3D (origin?).
       *Judging from the shape of the plots, this material also contain
       some silicates, though.*

-  **indices\_mc3d.dat:** “silicates” optical indices file used by
       Sebastian Wolf with MC3D; identical to the indices in file
       **dlsi\_opct.dat**, except that this file has no coverage for
       lambda < 0.2 micron. File **indices\_mc3d\_sil.dat** is
       identical.

-  ***lgsi\_opct.dat*:** amorphous olivine (MgFeSiO\_4), based on lab
       data and synthesized by Li, A., & Greenberg, J. M. 1997, A&A,
       323, 566, optical constants displayed as solid line in Figure 1
       of Li, A., & Greenberg, J. M. 1997, A&A, 323, 566

-  **porousice\_james.dat:** Water ice with optical indices only
       provided at two wavelengths (very limited usefulness)

-  ***waterice210K.dat*:** Water ice (at 210K?), as compiled in Pollack
       et al. 1994, ApJ, 421, 615

-  **PAHion.dat**: ionized PAHs computed by B. T. Draine (compiled from
       Laor, A., & Draine, B.T. 1993, ApJ 402:441, Draine, B.T., & Lee,
       H.M. 1984, ApJ 285:89, and Li, A., & Draine, B.T. 2001, ApJ
       554:778)

-  **PAHneu.dat**: neutral PAHs computed by B. T. Draine (compiled from
       Laor, A., & Draine, B.T. 1993, ApJ 402:441, Draine, B.T., & Lee,
       H.M. 1984, ApJ 285:89, and Li, A., & Draine, B.T. 2001, ApJ
       554:778)

Stellar Models
--------------

We should also add references for these.

7. Running a 3D model
=====================

MCFOST is completely 3D and can be used with any density structure.
Interfaces for hydrodynamics codes can be built on demand if there are
specific needs. There is also a default interface using FITS file to
input MCFOST with any arbitrary density structure. This interface is
likely to fulfill most of the needs.

It can be used with the command:

mcfost <parameter\_file> -density\_file <your\_density\_file.fits.gz>
-3D (+ any other option)

The density array in the FITS file must have 4 dimensions :

density(1:n\_rad, 1:nz, 1:n\_az, 1:n\_grains)

(it then has n\_rad x nz x n\_az x n\_grains elements).

The option -density\_file also works in the 2D version of the code, ie
without -3D. In that case the density file still requires 4 dimensions,
but with n\_az = 1 :

density(1:n\_rad, 1:nz, 1, 1:n\_grains)

| Note that the spatial grid is not currently passed to MCFOST in the
  FITS file. Instead, n\_rad, nz and n\_az must match the values
  provided in the parameter file. MCFOST will build the spatial grid
  from the parameter file (the grid can be cylindrical or spherical) and
  will assume the density is provided on this grid. The safest way to
  create the density file is then to first run mcfost <parameter\_file>
  -density\_struct and read the file data\_disk/grid.fits.gz to compute
  the density at the exact same points. As far as possible, you should
  try to match the MCFOST grid (ie, Rin, Rout, H\ :sub:`0`, p) to the
  hydrodynamics simulations, this is will ensure a best use of the
  MCFOST cells.
| The value n\_grains provided is the FITS file is independent of the
  one given in the parameter file on the other hand. MCFOST will assume
  that the n\_grains 3D density structures provided in the FITS file are
  independent and will interpolate (in log of the grain size) between
  these densities for all the grain sizes requested in the parameter
  file. MCFOST does not extrapolate, ie all the grains smaller (larger)
  than the minimum (maximum) grain size provided in the FITS file will
  have the same spatial distribution as the minimum (maximum) grain
  size.

The gas is assumed to follow the density of the smallest grain size
provided (you can however create an artificially small grain size for
the gas if you wish).

There are 2 options to normalize the different spatial distribution of
the grains relative to each other. These 2 options can be selected using
the keyword read\_n\_a in the header of the 1st HDU of the FITS file.

If read\_n\_a is set to 0, the density structure does not need to be
normalized. MCFOST will use the dust mass, grain size power law
distribution and gas-to-dust mass ratio provided in the parameter file.
In other words, the fits file is only used to pass the spatial
distribution of each grain sizes and MCFOST will take care of all the
required normalization.

If there are more than 1 grain size, the values of the n\_grains grain
sizes used in the FITS file in a second HDU with 1 dimension (NAXIS=1) :
NAXIS1 = n\_grains (unit : micron) :

grain\_sizes(1:n\_grains)

These grain sizes must be given in increasing order (and of course match
the order of the data cube).

If read\_n\_a is set to 1, the relative density in each grain size bin
needs passed to MCFOST and you must provide the grain size distribution
integrated over the whole disk volume. This is done by using a 3st HDU
where you pass an array of dimension(1:n\_grains)providing the number of
grains in each grain size bin : dn(a)/da

Note that the total normalization of the distribution is not required as
MCFOST will use the total dust mass as set in the parameter file. The
array n\_a is only used to define the relative number of the grains of
each size.

The gas density can be passed by setting the keyword read\_gas\_density
to 1. In that case, an extra keyword gas\_to\_dust must also be passed
to give the integrated gas to dust mass ratio.

| An extra HDU must also be passed, giving the gas density with 3
  dimensions:
| gas\_density(1:n\_rad, 1:nz, 1:n\_az).
| As for the dust density array, the array has 3 dimensions also in 2D,
  but with n\_az=1. The total gas mass is set by the mcfost parameter
  file in any case. If the keyword read\_gas\_density is set to 0 or
  absent, mcfost will assume that the gas has the same spatial
  distribution as the smallest grains in the FITS file, and will use the
  gas-to-dust mass ratio provided in the parameter file.

| The gas velocity field can also be passed by setting
  read\_gas\_velocity to 1. The fits file must have an extra HDU with 4
  dimensions:
| gas\_velocity(1:n\_rad, 1:nz, 1:n\_az, 3).
| The last index correspond to v\_x, v\_y, v\_z and the velocity is
  given in km/s.

If the velocity is not passed via the fits interface, mcfost will assume
that the velocity field is Keplerian (unless modified by command line
options).

8. Tools for Creating Parameter Files
=====================================

**From IDL:** There is a script ‘crea\_grid\_mcfost’ originally by M.
Perrin / G. Duchêne and recently updated by Elodie Choquet.

1. Edit a parameter file by hand, to set all the parameters you want
       constant.

2. | Run the IDL command crea\_grid\_mcfost, with keywords setting the
         parameters you want to vary. e.g.
       | crea\_grid\_mcfost, mdust=[1e-9, 1e-8, 1e-7], amax=[1,10,100]
       | would create a set of 9 models with varying dust mass and max
         grain size.

3. If you want to add additional elements onto an existing grid, call
       the same function again with the new desired parameters, and add
       the init\_counter= parameter set to the starting number to use
       (one greater than the number of models output in the prior run)

**From Yorick:** This tool is maintained by C. Pinte. Usage is similar
to the IDL interface. It also includes grid fitting and Bayesian
inference tools, a genetic algorithm, an ant colony algorithm and
parallel MCMC sampler. Default tool for Grenoble people ;-)

The IDL and yorick interfaces will be progressively deprecated and be
replaced by the python interface.

**From Python:** There are Python tools under development by Marshall
Perrin, Christophe Pinte, Schuyler Wolff et al., including both grid
creation, model display, and model fitting.

The code is available here :

`*https://github.com/cpinte/mcfost-python* <https://github.com/cpinte/mcfost-python>`__

Note that every time the parameter file format changes in a new version
of MCFOST, all of the above scripts must be updated to be compatible
with the changes..

9. Tools for Computing a Model Grid
===================================

The right way to do this depends on what machine you’re running on.

-  To run on a single computer, you can just use a for loop in a shell
       script, or the equivalent

-  On the FOSTINO cluster, use distrib\_mcfost.py or
       distrib\_addimages.py (ask Marshall) or the new yorick interface
       (ask Christophe)

-  On a Mac with xgrid set up, use mcfost\_runxgrid.py (Marshall)

10. IDL Tools for Visualizing MCFOST Output
===========================================

*Note - this tool is now deprecated. Use the Python tools instead. *

MCRE (MCFOST Results Explorer, author : Marshall) is a graphical
application for browsing and analyzing large (10k-100k+) model grids
computed using MCFOST. It has its own manual; please see that document.

There are a number of stand-alone simple tools that are bundled with
MCRE, and provide a simple IDL command line interface to plotting
various outputs.

-  **mcseds** Plot SED for all inclinations provided by a model

-  **mcsed** Plot SED for online one inclination

-  **mcimgs** Display images at all inclinations (for one wavelength)

-  **mcimg** Display image at one inclination and one wavelength

-  **mc\_rt** Replaced mcimg for raytraced image - should merge code.

-  **mcwaveimgs** Display images at all wavelengths (broken?)

-  **mcdiff** Overplot 2 SEDs for comparison

-  **mcrim** Plot puffed-up inner rim graphically

-  **mcgrid** Compute coordinates for each grid cell in a computation
       (broken?)

-  **mctemp** Display temperature plot

11. Python Tools for Visualizing MCFOST Output
==============================================

Marshall and Christophe have started developing a python module for
manipulating MCFOST output. Thus far it is very incomplete compared to
the IDL functionality (but does have nicer looking plots in general,
since it’s Python). This is available from github, including its own
documentation. Very much a work in progress, with contributions
welcomed.

The python-code is available via github :

https://github.com/cpinte/mcfost-python

12. Yorick tools for Vizualizing MCFOST Output
==============================================

The Yorick interface (author : Christophe) includes a number of command
line functions to display and manipulate various outputs. Detailed
documentation of the function is included in the yorick source code.
Here is a list of the most useful commands :

-  **open\_mcfost** read a mcfost model and put the parameter file and
       resulating calculations in a yorick structure

-  **read\_params** and **write\_params** read a parameter file into a
       yorick structure (and write a yorick structure to parameter file)

-  **plot\_sed** plot SED with extinction, can separate the various flux
       contribution

-  **interfero (and calc\_vis)** plot (or output) visibility curves

-  **plot\_temp** plot temperature structure

-  **plot\_dens** plot density structure

-  **surface\_density** compute and plot the surface density

-  **plot\_contrib** plot the various flux contribution in an image

-  **plot\_pola** compute and plot polarization maps from the Stokes
       maps

-  **sum\_mcfost** merge several Monte-Carlo models

-  **plot\_line** plot line data (integrated map and intergrated line
       profile)

-  **spectral\_index** compute the spectral index of a model between 2
       wavelengths

-  **partial\_mass** compute the model mass between 2 radii

-  **ProDiMo\_star\_2\_mcfost\_star** convert a ProDiMo ASCII input
       stellar spectrum to a MCFOST fits input stellar spectrum

-  **create\_grid** create a uniform grid of models

-  **sedfit** fit and vis a SED or visibility data set with a grid of
       models and compute the Bayesian inference

-  **start\_grid** set up and launch a grid on FOSTINO

-  **stat\_grid** check the status of a running grid on FOSTINO

-  **mcfost\_grid2ascii** output MCFOST grid fitting results to an ASCII
       file

-  **mcfost\_genetic** run a genetic model fitting (to be executed on
       FOSTINO)

-  **mcfost\_eMCMC** run a parallel MCMC (to be executed on FOSTINO)

13. Upgrading to New Versions
=============================

You can run mcfost -v to see the current version, and if a new version
is available

mcfost -u to automatically download & install the new version.

mcfost -update\_utils will similarly update the utilities directory.

mcfost -fu to force the update to whatever newest executable is
available, independent of version number. In other words, this overrides
the version number checking, precisely for instances like a small change
to fix a bug shortly after a new version is out.

From version 2.18.5, mcfost will check for updates automatically at
start-up if the last update is older than 7 days (this should take less
than 1 second). This behaviour can be changed by setting the environment
variable MCFOST\_AUTO\_UPDATE to an integer defining the number of days
between which mcfost will check for updates. If MCFOST\_AUTO\_UPDATE is
set to 0, mcfost will not check for updates automatically.

14. Troubleshooting
===================

Error: “ WARNING : first cell is in diffusion approximation zone. Increase spatial grid resolution”
---------------------------------------------------------------------------------------------------

The solution is to increase n\_rad\_in. You may need to try a few values
empirically to find one that works. Around 30 or 40 is recommended.

Christophe says:

*“It is an optical depth problem, you have more than 10^9 in the
midplane, which require a proper sampling of the grid. The solution is
to increase n\_rad\_in which is the number of subdivision of the first
cell.I would suggest to use :*

*150 60 1 30 n\_rad, nz (or n\_theta), n\_az, n\_rad\_in*

*to define the grid geometry.*

*I use (n\_rad - n\_rad\_in + 1) to build a log spaced radial grid but
by doing so, the optical depth of the first cell is often way larger
than 1 in the optical which prevent proper calculations.*

*So, I use n\_rad\_in to subdivide the first cell in sub-cells with a
sampling much concentrating towards the inner edge.*

*In your case, MCFOST is complaining because the optical depth of the
first cell is larger than the threshold chosen to define the diffusion
approximation zone, ie 1000, which means I cannot run the diffusion
approximation calculation (the direct flux from the star will make the
calculation to crash + most important, with such an opacity in the
first, the results will be completely wrong).*

*I generally use something like 30 or 40 for n\_rad\_in and between 100
and 200 for n\_rad and it seems to work fine with the models I have done
so far.”*

When I try to run mcfost, I get the error message “Exec format error. Binary file not executable”. What’s going on?
--------------------------------------------------------------------------------------------------------------------

You’re trying to run a Linux version of MCFOST on a Mac, or vice versa,
or some other incorrect operating system combination. You probably
downloaded or copied the wrong file by accident, and should just obtain
a new copy from the download site.

During SED calculations, I sometimes see warning printouts “PB r 0 1”. What does this mean?
--------------------------------------------------------------------------------------------

For instance, here is a portion of the MCFOST output:

127 / 128 2

128 / 128 4

WARNING : temperature > sublimation temperature

WARNING : temperature = 1680.505

Temperature calculation complete in 0h 0m 34s

0.1190514 1.000000

0.1687343 1.000000

0.2391511 1.000000

0.3389544 0.9999878

0.4804080 0.9984285

0.6808934 0.9520234

PB r 0 1

0.9650461 0.5454955

PB r 0 1

1.367782 0.1429528

1.938590 3.5336383E-02

PB r 0 1

2.747609 2.5369067E-02

3.894250 1.7464584E-02

PB r 0 1

5.519411 6.8406179E-03

7.822790 2.5935131E-03

etc.

*This indicates a resolution problem in the inner grid cells, i.e., the
gradient of temperature and/or density is too steep within a given cell,
leading to suspicious (at least) results. You should increase the total
number of pixels along the radial direction and/or increase n\_rad\_in
and that should solve it.*

*The same kind of messages with* PB z *indicate a resolution issue at
the disk surface. It probably means that the scale height is changing
too rapidly between 2 zones. Is that really what you wanted to do ?*

References
==========

Hogerheijde, M. R. & van der Tak, F. F. S. 2000,
`A&A <http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?2000A%26A...362..697H&amp;data_type=PDF_HIGH&amp;whole_paper=YES&amp;type=PRINTER&amp;filetype=.pdf>`__,
`362, 697 <http://adsabs.harvard.edu/abs/2000A%26A...362..697H>`__

Pinte, C., Harries, T. J., Min, M., et al. 2009,
`A&A <http://www.aanda.org/articles/aa/pdf/2009/18/aa11555-08.pdf>`__,
`498, 967 <http://adsabs.harvard.edu/abs/2009A%26A...498..967P>`__

Pinte, C., Ménard, F., Duchêne, G., & Bastien, P. 2006,
`A&A <http://www.aanda.org/10.1051/0004-6361:20053275/pdf>`__, `459,
797 <http://adsabs.harvard.edu/abs/2006A%26A...459..797P>`__

.. |image0| image:: media/image2.png
   :width: 2.68229in
   :height: 2.17076in
