MCFOST + Phantom
================


MCFOST can read both the standard and hdf5 phantom dump formats.
Typical syntax for reading a phantom file is:``mcfost <para_file> -phantom <phantom_dump>``
Any additional option can be used in conjonction with the `-phantom` option.

When reading a phantom dump, sections of the parameter file will be ignored, in particular any section describing the disc density, model grid and stellar properties (unless the stellar properties are forced).


Here are some important sections to consider for post-processing Phantom dump files.

Number of photons for temperature and image calculation
#######################################################

You may need to use more photons to compute the temperature. As a rule of thumb, you should use at least 100 more packets than phantom SPH particles, ie the section is the parameter file should be something similar to

::

 #Number of photon packages
   1.28e8                  nbr_photons_eq_th  : T computation
   1.28e3                  nbr_photons_lambda : SED computation
   1.28e6                  nbr_photons_image  : images computation



Turn off SED computation
########################

You probably don't want to compute an SED if you have a Phantom dump file.

::

   T F T                   compute temperature?, compute sed?, use default wavelength grid for ouput ?


Turn off image symmetries
#########################

You don't want image symmetries if you have a Phantom dump file.

::

 #Symetries
   F                       image symmetry
   F                       central symmetry
   F                       axial symmetry (important only if N_phi > 1)

Can ignore the following, they come from the phantom dump
#########################################################

::

 #Grid geometry and size
 #Disk physics
 #Number of zones : 1 zone = 1 density structure + corresponding grain properties

Density structure
#################

Density comes from phantom, so can ignore "#Density structure" however you need
to set "gas-to-dust mass ratio":

::

   1.e-3    100.           dust mass,  gas-to-dust mass ratio


Star properties
###############

By default, mcfost will use the properties (ie mass) of the sink particles to estimate the brightness temperature and luminosity.
mcfost will do so by interpolating through isochrones, assuming an age of 3Myr by default. The age can be change with the `-age` option.
The stellar properties can be fixed to the values in the parameter files with `-fix_stars`. Note that in that case, the number of stars in the mcfost parameter file must matches the number of sink particles in the phantom dump.
