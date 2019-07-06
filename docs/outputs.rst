Outputs from MCFOST
===================

Outputs from SED mode calculations
----------------------------------

A calculation in SED mode will create a direcatory named ``data_th`` with the following files:

-  A copy of the input parameter file. A new line will be added at the end of
   this file, recording any command line options given in the call to MCFOST.

-  ``Temperature.fits.gz`` (and ``Temperature_DA.fits.gz`` if diffusion
   approximation is used) + additional temperature files of out-of-equilibrium grains are
   present

-  ``sed_rt.fits.gz``. This is the main output with the SED computed with a
   ray-tracing method. Ray-ytraced SEDs offer several advantages : they correspond to exact inclinations
   (instead of bins for ``sed_mc.fits.gz``) and converge almost equally fast at
   all inclinations.

-  ``sed_mc.fits.gz``. This file contains the SED computed with the MC
   method. It is mainly kept because it is small and can be used to
   double-check that the calculation went allright.

Units for all SED files are W.m\ :sup:`-2`. The temperature map is in K,
with values of 1K for regions where the temperature was not computed and
is later estimated via the diffusion approximation (this occurs in
extremely deep regions of disks).

Output SED files typically are 4 dimensional. The dimensions are
[wavelength, inclination, azimuth, contributions]. Dimensions may be of length
1; for instance typically the azimuth axis has only 1 value for the case
of axisymmetric models.

An additional ``.sed_th.fits.gz`` is also created in the ``data_th``
directory. It records the SED produced by the Bjorkman & Wood algorithm
during the calculation of the temperature. This file is generally noisy
at wavelengths where the flux is low. It is mainly used for testing
purposes. In particular, each re-emission takes place at the exact same
location of the absorption, while in ``sed_mc.fits.gz`` and
``sed_rt.fits.gz``, the emissivity is assumed to be constant within a cell.
The comparison of ``.sed_th.fits.gz`` and ``sed_mc.fits.gz`` can be useful to
detect inadequate spatial resolution (ie cells that are too big to
accurately sample the local emissivity).

Outputs from Image calculations
-------------------------------

-  A copy of the input parameter file. A new line will be added at the
   end of this file, recording any command line options given in the
   call to MCFOST.

-  ``RT.fits.gz``: ray-traced output image

-  ``MC.fits.gz``: Monte-Carlo output image (if -mc option is used)

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

.. note:: the -casa option makes the fits file compliant with casa and reduces the number of dimensions. Only the first inclination and azimuth are considered. Additional keywords are added for casa, units are changed to Jy/pixel.


Outputs from Molecular Line calculations
----------------------------------------

-  A copy of the input parameter file. A new line will be added at the
   end of this file, recording any command line options given in the
   call to MCFOST.

-  ``lines.fits.gz``: ray-traced data cube

Units are in  W.m\ :sup:`-2`.pixel\ :sup:`-1`. These images have 6 dimensions, 2 for spatial dimensions,
1 for velocity, 1 for the transition number, and 2 for the system inclination and azimuth.

The files contains multiple HDUs.
The 2nd HDU contains the continuum images, it has the same dimension as the main cube without the velocity.
The 3rd and 4th HDUs contain the transition numbers and frequencies respectively.


.. note:: the -casa option makes the fits file compliant with casa and reduces the number of dimensions. Only the first inclination and azimuth are considered.
          The fits file also only contains the main HDU in that case.


Additional optional outputs
---------------------------

Disk Density Files
^^^^^^^^^^^^^^^^^^

You can request an output of the computed disk density. This is done
with the ``-disk_struct`` (or ``-output_density_grid``) command line option.
This will output files in a ``data_disk`` directory:

-  ``gas_density.fits.gz`` [map of the gas density - in g/cm\ :sup:`3`]

-  ``dust_mass_density.fits.gz`` [map of the integrated dust density - in g/cm\ :sup:`3`]

-  ``dust_particule_density.fits.gz`` [map of the dust density - in part/cm\ :sup:`3`]

-  ``grain_sizes.fits.gz`` [grain size bins - in micron]

-  ``grid.fits.gz`` [map of the radii and height of the individual grid cells - in au]

-  ``volume.fits.gz`` [map of the volume of the individual grid cells - in au\ :sup:`3`]

The various density files will have 3d arrays of dimension [``n_rad``, ``n_z``
or ``n_theta``, ``n_az``], where those parameters are as defined below in the
parameter file geometry. The ``grid.fits.gz`` file will have those same
dimensions for its first 3 axes, then the 4th axis will have dimension
2. The first plane in that axis will give the radius for each grid cell,
the second will give the height.

Dust Property Files
^^^^^^^^^^^^^^^^^^^

You can request an output of the computed dust propertis. Thies is done
with the ``-dust_prop`` command line option. This will generate a
``data_dust`` directory that will contain the following output files (along
with a copy of the parameter file):

-  ``albedo.fits.gz`` [table of the dust albedo]

-  ``g.fits.gz`` [table of the scattering phase function asymmetry parameter]

-  ``kappa.fits.gz`` [table of dust opacities - in cm\ :sup:`2`/g]

-  ``kappa\_grain.fits.gz`` [table of dust opacities as a function of grain size - in cm\ :sup:`2`/g]

-  ``lambda.fits.gz`` [table of wavelength at which the dust properties are computed - in micron]

-  ``phase\_function.fits.gz`` [table of scattering phase function from 0 to 180 degrees, versus wavelength]

-  ``polarizability.fits.gz`` [table of scattering polarizability from 0 to 180 degrees, versus wavelength]

The ``lambda.fits.gz`` file also provides the wavelength array for an SED
computation when ``lsed_complete`` is set to True (see below).

If the dust properties are only needed at a specific wavelength (say 2.2
microns), the user can instead use the ``-op 2.2`` command line option. In
this case, only the initialization sequence is performed, and the albedo
and asymmetry parameter are given on screen at the end (as well as the
total optical depth along the equatorial plane and along the viewing
"inclination of interest" indicated in the parameter file). Note that
these quantities are always provided in the normal run of MCFOST (for
both SED and image calculations); this option merely stops the
calculation after this stage. Note that a ``data_2.2`` directory is created
containing only a copy of the parameter file.

Optical depth maps
^^^^^^^^^^^^^^^^^^

tbw

tau=1 surface
^^^^^^^^^^^^^

tbw

Tracking origin of packets
^^^^^^^^^^^^^^^^^^^^^^^^^^

tbw

Radiation field
^^^^^^^^^^^^^^^
tbw

Files for ProDiMo
^^^^^^^^^^^^^^^^^

tbw

Files for Astochem
^^^^^^^^^^^^^^^^^^

tbw
