Command line options
====================

Most of this section remains to be written.

Stellar spot
------------

A spot can be added to photosphere using the option::

$ mcfost <para> -spot <T_spot> <spot_surface_fraction> <i_spot> <phi_spot>

The spot surface fraction is defined between 0 and 1. The spot
inclination is defined between 0 and 180 degrees: 0 degree for a spot on
the pole and 90 degrees for a spot on the equator.

The option currently only works if the photosphere has no extra UV and
is only available in MC mode so far.

If i is defined between 0 and 90 degrees, the spot is facing the
observer when using the central azimuthal bin.
