Command line options
====================

Most of this section remains to be written.

Basic options
-------------

``-h``, ``-help``: displays the help message, with a summary of the main
options. Running `mcfost` withour any argument displays the same message.

``-v`` : displays version number, and available updates

``-setup`` : performs mcfost initial setup (downloads the data and parameter files)

``-get_para``: downloads the current version of the parameter file

``-u``: updates MCFOST to most recent version

``-fu``: forces MCFOST to update to the most recent version, even if there is
no new release available

``-update_utils``: updates MCFOST_UTILS to most recent version

``-history``: displays full MCFOST history since v2.12.9

``-max_mem <value>`` [GB]: maximum memory that MCFOST can use (approx), default 8

Main options
------------

``-img <wavelength>`` [microns]: computes image at specified wavelength

``-mol``: calculates molecular maps




Stellar Properties
-------------------

``-spot <T_spot> <spot_surface_fraction> <i_spot> <phi_spot>``: adds a cold or
hot spot on the photosphare. The spot surface fraction is defined between 0 and 1. The spot
inclination is defined between 0 and 180 degrees: 0 degree for a spot on
the pole and 90 degrees for a spot on the equator. The option currently only
works if the photosphere has no extra UV and is only available in MC mode so far.
If i is defined between 0 and 90 degrees, the spot is facing the
observer when using the central azimuthal bin.

``-limb_darkening <ld_filename>`` adds limb darkening (including polarized limb
darkening) on the stellar photosphere. Exemples of limb darkening files can be found in
``$MCFOST_UTILS/Stellar_Polarization``.

``-age <age>`` (default: 3) [Ma]. When using results from hydrodynamics
simulations (e.g., SPH), mcfost  will assume an age to determine the stellar
luminosity and temperature from the mass. The age can be selected using::
Isochrones are found in ``$MCFOST_UTILS/Stellar_Polarization/Siess1``.


Coupling with other codes
-------------------------

``-prodimo``: creates required files for ProDiMo.

``-p2m``: reads the results from ProDiMo.

``-astrochem``: creates required files for Astrochem.

``-phamtom`` : reads a phantom dump file.

``-gadget`` : reads a gadget-2 dump file.

``-limits <limit-file>`` : x,y,z values used for the Voronoi tesselation.

``-keep_particles <fraction>`` (default: 0.99): fraction of SPH particles to
keep for the Voronoi tesselation.
