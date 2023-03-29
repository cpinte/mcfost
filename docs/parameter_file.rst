Parameter File Contents
=======================

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

* ``nbr_photons_eq_th``: Number of of photon packets launched when computing
  the temperature structure of the disk (typically 100 000, unless large
  noise is seen in temperature map),

* ``nbr_photons_lambda``: Number of photons packets **per wavelength** that must
  reach the observer when computing the output SED, for the inclination
  bin with highest inclination. Typically 10 000, unless
  anomalously large noise is seen in the SED. With 10 000 packets
  received, the MC noise will be of the order of 1%.

* ``nbr_photons_image``: Number of of photon packets launched when computing a
  disk image. This number must be adjusted depending on the spatial and
  inclination resolution. Low values, typically 100 000, are enough when
  using the ray-tracing mode.


Wavelength setup
----------------

* ``nlambda``, ``lambda_min``, ``lambda_max``: Number of wavelength bins, and
  minimum and maximum wavelengths in units of microns.

.. important::  ``lambda_min`` and  ``lambda_max`` are defining the
                range of wavelengths used to compute the temperature structure. They must
                include all wavelengths were there is significant emission from the star and
                the disk.

.. note:: ``lambda_min`` and ``lambad_max`` are the extreme wavelengths at the edges of
          the bins; the central wavelengths of the bins will be a little bit
          inside those extreme values. If you need to force specific wavelengths,
          use the ``use_default_wavelength_grid`` flag.

* ``ltemp:`` do you want to compute the temperature map? Most of the time
  it's T (unless you have already previously computed the temperature but
  not the SED for some reason).

* ``lsed:`` do you want to compute an SED? Most of the time it's T (unless
  you care only about the temperature structure of the disk)

* ``use_default_wavelength_grid``: are you happy with the default evenly log-spaced
  sampling in wavelength (T) or do you want to use the wavelength file
  indicated in the following line (F)

* ``wavelength_file``: a text file containing the wavelengths at which the
  SED is to be sampled. One entry per line, all entries in microns. The
  file should be located in ``$MCFOST_UTILS/Lambda/`` or in the local
  directory where MCFOST is executed, and a copy will be placed in the
  ``data_th``  directory

* ``separate_contributions:`` Save output file that tracks the various
  different contributions (star vs disk, scattered vs direct). Applicable
  to both SED and image calculations, works in RT mode. In MC SED mode,
  this setting is *always* on, as it is almost free to track this.

* ``output_stokes_parameters``: Keep track of the Stokes parameters for
  each photon. Output images will include polarization axis. Only
  applicable to image calculations, works in RT mode.

Grid geometry
-------------

* ``grid_type``: selecting a cylindrical, spherical grid geometry.
  The spherical grid
  is mandatory if an envelope is present; for disks, the cylindrical
  grid is more appropriate, and benefits from more optimizations.

.. note:: A Voronoi tesselation is also available for a set of pre-defined points.
           The Voronoi tesselation is  automatically used when using some density input
           files. It is particularly suited to use the results of
           hydrodynamical simulations, in particular SPH, where 1 SPH particle
           can be directly matched to 1 cell in MCFOST.

* ``n_rad``: Number of grid cells along the radial direction

* ``nz (n_theta)``: Number of grid cell along the vertical direction
  (from the midplane up, for a cylindrical grid), or the latitude
  direction (for a spherical grid)

* ``n_az``: Number of grid cells along the azimuthal direction (set to 1
  for an axisymmetric structure)

* ``n_rad_in``: Number of sub-cells in which to split the first cell
  along the radial axis. This should be in the 10-30 range for an
  optically thick disk so as to deal accurately with the highest density
  regions of the disk. If the temperature map shows T=1K in the second
  pixel from the disk inner radius, you need to increase ``n_rad_in``.


Images
------

* ``grid (nx, ny)``: Number of pixels for the output images (applies to
  both Monte Carlo and RayTracing images).

* ``map_size``: physical size of the produced maps (in au). Applied to
  the largest axis if nx is different from ny

* ``RT_imin``: Minimum inclination to compute the Ray Tracing images
  (does not need to be 0\ :sup:`o`). 0\ :sup:`o` is pole-on and
  90\ :sup:`o` is edge-on.

* ``RT_imax``: Maximum inclination to compute the Ray Tracing images
  (does not need to be 90\ :sup:`o`)

* ``RT_n_incl``: Number of inclinations for which Ray Tracing images
  should be computed in the [imin,imax] range, evenly spaced in cos(i)

* ``RT_centered``: Set to F to force images to be computed exactly at
  imin and imax; if set to T, the sampled inclination values will be
  computed at the midway point of the cos(i) bins

* ``RT_az_min``: Minimum azimuthal angle to compute the Ray Tracing
  images (does not need to be 0\ :sup:`o`).
  For a pole-on model, 0\ :sup:`o` means that the x and y axes of the model
  correspond to the x and y axes if the map.

* ``RT_az_max``: Maximum azimuthal angle to compute the Ray Tracing
  images.

* ``RT_n_az``: Number of azimuthal angles for which Ray Tracing images
  should be computed in the [``RT_az_min``, ``RT_az_max``] range, linearly
  spaced.

* ``distance``: Distance to the object in pc

* ``disk_PA``: position angle of the disk major axis from North to East.
  PA = means that the red-shifted side of the disk is towards North.

.. important:: For version below 3.0.43, PA was defined as the
  position angle of the semi-minor axis of the disk, measured
  counter-clockwise, with red-shifted side towards West. This
  convention can be recovered with comman line option -old_PA
  (this removes 90 degrees from the PA set in the parameter file).


  If disk PA is not set to +/-90\ :sup:`o`, some of the
  image symmetries (see below) will automatically be set to F

Scattering Method
-----------------

``scattering_method``: compute the average dust properties per grain
size (value: 1), per grid cell (value: 2), or the most appropriate for
the ongoing calculation (value: 0, preferred in virtually all cases).
Value 2 can not be used if the Stokes parameters are required.

``Mie/hg``: choice of scattering phase function; Mie theory (value: 1)
is strongly preferred over the Henyey-Greenstein parametric description
(value: 2) since it is physically grounded and allows calculations in
full-Stokes mode. If the HG phase function is selected, Mie theory will
first be used to compute the effective *g* value, which will then be
used to randomly select scattering angles.

.. note:: Optical properties calculations are much faster when the option is
          set to 2, reducing significantly MCFOST's initialization time. It
          may be useful to compute quickly temperature structures and SEDs
          when very large grains are present and/or a large number of
          wavelengths are used (the overhead due to Mie theory is almost
          always negligeable for monochromatic images).

Symmetries
----------

* ``image_symmetry``: Is the image left/right symmetric?

* ``central_symmetry``: Is the model structure symmetric relative to the origin?

* ``plane_symmetry``: Is the model structure symmetric relative a vertical plane?

In most cases for a 2D model, all three symmetries should be set to T. If you
use an asymmetric disk structure, or if there is more than one star
illuminating the disk, then they should be set to F. If a disk PA
different than 0\ :sup:`o` has been set, then first two symmetries
should be set to F and the last one to T.

Disk physics
------------

* ``dust_settling:`` no settling (0), parametric (1), following
  Dubrulle's (2) or Fromang's prescription (3).

.. note::
    * Parametric settling is independent of radius, it is simply a
      scaling of the scale height as a function of the grain size.

    * Dubrulle's and Fromang's settling assume a diffusion equation
      depending on the viscosity. Dubrulle's vertical profile remains
      Gaussian, while Fromang's presciption is more realistic and departs
      from the Gaussian profile at high altitude. (see Dubrulle et al 1995
      and Fromang et al 2009).

* ``exp_strat``: power law describing the parametric settling, where H(a)
  decreases as a\ :sup:`-exp_strat` for a > a_strat

* ``a_strat``: minimum size [microns] for grains affected by vertical
  settling (for parametric settling)

* ``dust_radial_migration``: simple prescription for radial migration of
  dust grains (TBW)

* ``hydrostatic_equilibrium``: (*work in progress*) compute the
  hydrostatic equilibrium assuming Tgas = Tdust

* ``sublimate_dust``: this option will iteratively remove the dust if the
  temperature reaches the dust sublimation temperature

* ``viscous_heating``: (*work in progress*) includes additional heating
  source due to viscous accretion (note that it does not account for the
  accretion shock on the star)

* ``viscosity``: alpha parameter describing the strength of the viscosity
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

* ``zone_type``: disk (value: 1), disk with outer tapered-edge (value:
  2), spherically symmetric envelope (value: 3), debris disk (value : 4)
  or an azimuthally asymetric wall (value : 5).
  If at least one of the zones is described as an envelope, the
  computing grid must be spherical.

* ``disk_dust_mass``: in units of M\ :sub:`sun`

* ``gas_to_dust_ratio``: quantity only used in molecular emission
  calculations

* ``scale_height_H0:`` value of the disk scale height (technically, sigma
  of the Gaussian vertical density profile) at the reference radius, in au

* ``reference_radius_R0``: in au

* ``vertical_profile_exponent``: exponent of vertical
  density profile (only relevant for type 4, debris disk)

* ``Rin``: inner radius of disk or envelope in AU (assuming sublimation
  calculation is not enabled)

* ``Rout``: outer radius of disk or envelope [au]

* ``Rc``: critical radius for tapered-edge disk model [au]

* ``edge``: make nonzero for gradual rather than abrupt falloff inside
  rin/outside rout. This is to set a smooth decline in surface density at the inner and
  outer edges of the disk. Inside of r_in (and outside of r_out),
  the density drops following a Gaussian whose sigma is the "edge"
  parameter.

* ``\beta``: flaring exponent = power law index of the H(r) scale height
  function, typically in the [1.0-1.25] range

* ``p1`` and ``p2`` : surface density exponent, or ``-gamma`` for tapered-edge.
  Power law index of the surface density profile (generally <0).
  ``-gamma_exp`` for tapered-edge. Defines ``p_in`` and ``p_out`` in case of a
  debris disk.

The disk density structures are defined as :

1. :math:`\Sigma(r)\ \alpha \ r^{p1}`

2. :math:`\ \Sigma(r)\ \alpha \ r^{p1}\ \exp( - (\frac{r}{Rc})^{2 + p2})\
   \alpha \  r^{-\gamma}\ \exp( - (\frac{r}{Rc})^{2 - \gamma_{\exp}})`

3. :math:`\rho(r)\ \alpha\ r^{p1}`

4. :math:`\rho(r,z)\ \alpha \ ((\frac{r}{Rc})^{- 2p_{\text{in}}} +
   (\frac{r}{Rc})^{- 2p_{\text{out}}})^{- 1/2} \times \exp( -(\frac{\left| z \right|}{h(r)})^{\gamma_{\text{vert}}})`
   see Augereau et al, 1999, A&A, 348, 557. p\ :sub:`in` = p1 > 0
   and p\ :sub:`out` = p2 < 0

If ``Rout`` is set to 0, it is automatically set ``Rout`` to ``8 Rc`` in the case of
a disk with tapered-edge.

For types 1,2, 4 and 5, the local scale height is defined as :math:`h(r) =
h_0  (\frac{r}{r_0})^{\beta}`


Grain properties
----------------

.. note:: there should be as many blocks containing the following parameters
          as there are zones in the disk. The code will crash otherwise.

* ``n_species``: Number of dust populations present in the disk zone; if
  N_species > 1, the dust grains of different species are assumed to be
  physically disjoint, but distributed in the same manner through the
  disk. *All the following lines in the block must be duplicated
  N_species times if N_species > 1.*

* ``Grain_type:`` spherical grains (Mie) or distribution of hollow
  spheres (DHS)

* ``n_components``: Number of materials that make up a given specie;
  these materials are assumed to be physically joint within each dust
  grain. *The line with the optical indices and volume fraction must be
  duplicated N_components times if N_components > 1.*

* ``mixing_rule:`` are the components randomly mixed within the volume of
  a grain (value: 1, effective mixing theory following Bruggeman rule) or
  is the second component forming a coating on top of the first one
  (value: 2). The effective optical index of the new "mixed" grain is
  computed before any Mie theory computation. Does not apply is
  N_components = 1.
  *Coating can only be used with 2 components (the 1st one is the core,
  the 2nd one the shell).*

* ``porosity``: porosity of the dust grains (in the [0,1] range, 0 for
  compact grains, near 1 for porous ones)

* ``mass_fraction``: fraction of the mass contained in this specie (the
  sum of the N_species mass fractions should be equal to 1, MCFOST will
  renormalize the values so that the sum is 1)

* ``DHS_Vmax``: maximum void fraction for DHS calculation

* ``optical_indices_file``: file containing the optical index of the material
  as a function of wavelength (files must be located in
  ``$MCFOST_UTILS/Dust/``)

* ``volume_fraction``: fraction of the volume of a grain contained in
  this component (the sum of the N_components volume fractions should
  be equal to 1, MCFOST will renormalize the values so that the sum is 1)

* ``heating_method``: indicated whether radiative equilibrium and local
  thermal equilibrium are assumed: for an optically thick disk, both
  should be true (value: 1); for an optically thin disk, only the RE is
  assumed (value: 2, will yield a temperature map that has a third
  dimension spanning the grain size distribution); for out-of-equilibrium
  grains, none of them is true (value: 3, typical for PAH grains)

* ``amin``: minimum grain radius, in microns

* ``amax``: maximum grain radius, in microns

* ``aexp``: power law index of the grain size distribution dN(a)/da

* ``n_grains``: number of grain size bins to the sample [evenly in
  log(a)] the grain size distribution; dust properties are only computed
  for these grain sizes and subsequent interpolations are used whenever
  necessary. Typical value is in the 50-100 range, less in case of
  multiple dust specie/disk zones to limit RAM requirement (and
  computation time of the Mie theory).

Molecular RT settings
---------------------

* ``lpop:`` do you wish to compute the level populations (this might not
  be the case if you use populations from an external code, ProDiMo for
  instance)

* ``lpop_accurate:`` if the variable is set to false, mcfost will just
  perform a 1+1D (or 1+1+1D) line transfer. If set to true, the result of
  the 1+1D line transfer will be used as a starting point for the full 2D
  or 3D line transfer.

* ``LTE:`` assume LTE level populations

* ``profile_width:`` internal line width used for the line transfer
  calculation. Bascically, it means that cells with relative projected
  velocities that exceed this value will not see each other during the
  transfer. The value needs to be larger than the local line width.

* ``molecular_data_file:`` LAMBDA data file used for the line transfer

* ``level_max:`` maximum level up to which the line transfer will be
  performed. Level above ``level_max`` will not be populated

* ``vmax [km/s]:`` maximum velocity is the produced channel maps

* ``n_speed:`` number of velocity points between 0 and vmax

* ``cst_abundance:`` do you wish to use a constant molecular abundance
  over the disk. If true, use the provided abundance, if false read the
  abundance from the following fits file. The resolution of the fits file
  must be the same as the grid used by mcfost.

* ``ray_tracing:`` produce or not a ray-traced data cube of the molecule

* ``n_lines:`` number of transition ray-traced, the indices given
  in the next lines correspond to the transition indices in the LAMBDA
  files. For instance the J=1-0 transition is usually #0

Star properties
---------------

* ``n_stars:`` number of stars illuminating the disk. *If num_stars >
  1, the following group of two lines must be duplicated num_stars times*

* ``Teff``: effective temperature of the star (only used to compute the
  stellar luminosity), in K

* ``Rstar``: stellar radius (used to compute the stellar luminosity and
  set the spatial origin of the photon packets: the star is assumed to be
  a uniformly radiating sphere), in R\ :sub:`sun`

* ``Mstar``: stellar mass (only used for molecular line calculations,
  hudrostatic equilibrium and viscous heating via accretion, in
  M\ :sub:`sun`

* ``x``, ``y``, ``z``: position of the star, in au. If this is not (0,0,0), the
  image symmetries must be set to F

* ``automatic_spectrum?``: should the stellar spectrum estimated automatically from the effective temperature and radius ?
  If not, the stellar spectrum indicated in the following line is used instead

  .. note:: ``Teff`` vs stellar atmosphere model. The effective
          temperature is only used to compute the total stellar
          luminosity and it doesn't need to match the atmosphere
          model. In practice, the atmosphere model will be chosen
          with the effective temperature closest to the target's. In any case, the luminosity will be set by
          ``Teff`` and ``Rstar`` (assuming there is no UV excess).


* ``fUV, slope_fUV:`` photospheric UV excess and its slope in Fnu.

  .. note::  The basic underlying assumption is that the total emission from the
             star can be approximated by whatever stellar model you use plus a
             single power law UV excess. That UV excess is characterized by two
             parameters: fUV is the scaling factor while slope_fUV is the power
             law index. The definition of fUV can be found in `Woitke et al.
             2010 <http://adsabs.harvard.edu/abs/2010MNRAS.405L..26W>`__, and is
             the ratio of the UV flux between 92 and 250nm to the total
             photospheric luminosity. If fUV=0, then there's no excess on top of
             the phtosphere. The UV is only important for chemistry purposes in
             practice, unless you use some crazy high value of fUV.
