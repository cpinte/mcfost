MCFOST Overview
==================

MCFOST is a 3D continuum and line radiative transfer code based on the
Monte Carlo method (plus ray-tracing). It is mainly designed to study the circumstellar
environment of young stellar objects. The calculations are done exactly
within the limitations of the Monte Carlo noise and machine precision,
i.e. no approximation are used in the calculations. The code has been
strongly optimized for speed.

MCFOST is primarily designed to study protoplanetary disks. The code can
reproduce most of the observations of disks:

* SEDs
* scattered light images
* IR and mm visibilities
* atomic and molecular line maps

The Monte Carlo method being generic, any complex structure can be
handled by MCFOST and its use can be extended to other astrophysical
objects. For instance, tests have been performed on infalling envelopes
and AGB stars.

.. important:: The code is not public domain yet but available on request on a
               collaborative basis, with explicit permissions of the authors
               (see disclaimers below).

The core of the algorithms are described in
`Pinte et al. (2006) <http://adsabs.harvard.edu/abs/2006A%26A...459..797P>`__
and `Pinte et al. 2009
<http://adsabs.harvard.edu/abs/2009A%26A...498..967P>`__.
However the code has been substantially further
enhanced, with features that are mainly documented in this documentation.
In short, the code computes the temperature and scattering
source function everywhere in the disk via a Monte Carlo method: photon
packets are propagated stochastically through the model volume following
the equations of radiative transfer, and information on their properties
is retained along their path. The radiation field, and quantities
derived from it (for instance temperature, radiation pressure, etc) are
obtained by averaging this "Monte Carlo" information. Observables
quantities (SEDs and images) are then obtained via a ray-tracing method,
which calculates the output intensities by integrating formally the
source function estimated by the Monte Carlo calculations. Full
calculations of the polarization are included using the Stokes
formalism.


MCFOST also includes a non-LTE line transfer module. The adopted schemes
are not yet described in the literature, but the code uses an improved
version of the algorithm presented in `Hogerheijde & van der Tak (2000)
<http://adsabs.harvard.edu/abs/2000A%26A...362..697H>`__.
NLTE level population are obtained via iterations between Monte Carlo
radiative transfer calculations and statistical equilibrium. To speed up
calculations, an initial guess is computed using a 1+1D scheme with
short characteristic before the full 2D or 3D calculation. Output
spectra and channel maps are calculated via a ray-tracing procedure.


.. note:: MCFOST is in constant development and this documentation is very
          likely to be lagging. Please contact Christophe Pinte for the latest updates.


Associated with MCFOST, but separate software, are tools for
visualization, grid calculations and model fitting tools for SEDs,
images and line emission. These tools allows one to draw robust
constraints on the derived parameters.


.. warning:: Use at your own risk!!! The author does not take any
             responsibility for the use (or misuse of the code). There might be
             bugs.


Disclaimers
---------------

MCFOST is available on a **collaborative basis**. Using MCFOST implies that you agree to :

*  offer us (C. Pinte, F. Menard, G. Duchene) co-author right on any resulting publication.
*  NOT distribute MCFOST without our explicit agreement.
*  contact us if you initiate a new scientific project with MCFOST.

Please also acknowledge funding from the Australian Research Council
under contracts FT170100040 and DP180104235,
from Agence Nationale pour la
Recherche (ANR) of France under contract ANR-16-CE31-0013.

The IDL code MCRE (MCFOST Results Explorer) is also available on a
collaborative basis. Using MCRE implies to offer M. Perrin co-author
right. However its use is deprecated. Instead users are encouraged to
use the Python mcfost package, available from
`github <https://github.com/cpinte/mcfost-python>`__,
by M. Perrin, C. Pinte, and S. Wolff.


Slack channel
-------------

If you decide to use MCFOST, you should sign up to the `slack channel <https://mcfost.slack.com/>`__ in
order to be informed of new releases, bugs. This is also a collaborative space
to discuss issues and future developments of MCFOST.


..  If you decide to use MCFOST, you should sign up to the mailing list in
    order to be informed when new versions are available. Please send an
    email to `*sympa@ujf-grenoble.fr* <mailto:sympa@ujf-grenoble.fr>`__
    with header :
    subscribe mcfost <First Name> <Last Name>
    First name and Last Name are optional.
