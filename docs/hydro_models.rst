
Running mcfost on a hydro model
===============================

MCFOST is completely 3D and can be used with any density structure.
Interfaces for hydrodynamics codes can be built on demand if there are
specific needs. Current interfaces where mcfost can natively read dumps exist for

* phantom
* fargo-3d
* athena++
* idefix
* pluto


Generic interface (via fits files)
##################################

There is also a default interface using FITS file to
input MCFOST with any arbitrary density structure. This interface is
likely to fulfill most of the needs if you are not using one the hydro codes above.

It can be used with either a structured or unstructed grid.

Structured grid
---------------

For a structured grid, the fits interface can be used as::

    $ mcfost <parameter_file> -density_file <your_density_file.fits.gz> -3D (+ any other option)

or::

    $ mcfost <parameter_file> -density_file <your_surface_density_file.fits.gz> (+ any other option)

The `-sigma_file` option works in a similar way as the `-density_file` option but without the vertical dimension. mcfost use the scale height and flaring index provided in the parameter file to reconstruct the 3D density structure.


.. important:: When using the `-density_file` or `-sigma_file` options, the number of zones must be set to 1 in mcfost


The density array in the FITS file must have 4 dimensions : ``density(1:n_rad, 1:nz, 1:n_az, 1:n_grains)``
(it then has n_rad x nz x n_az x n_grains elements).

The option ``-density_file`` also works in the 2D version of the code, ie
without ``-3D``. In that case the density file still requires 4 dimensions,
but with n_az = 1 : ``density(1:n_rad, 1:nz, 1, 1:n_grains)``

.. important:: the spatial grid is not currently passed to MCFOST in the
          FITS file. Instead, ``n_rad``, ``nz`` and ``n_az`` must match the values
          provided in the parameter file. MCFOST will build the spatial grid
          from the parameter file (the grid can be cylindrical or spherical) and
          will assume the density is provided on this grid. The safest way to
          create the density file is then to first run::

            $ mcfost <parameter_file> -density_struct

          and read the file ``data_disk/grid.fits.gz`` to compute
          the density at the exact same points. As far as possible, you should
          try to match the MCFOST grid (ie, Rin, Rout, H\ :sub:`0`, p) to the
          hydrodynamics simulations, this is will ensure a best use of the
          MCFOST cells.

.. note:: The value ``n_grains`` provided is the FITS file is independent of the
          one given in the parameter file on the other hand. MCFOST will assume
          that the ``n_grains`` 3D density structures provided in the FITS file are
          independent and will interpolate (in log of the grain size) between
          these densities for all the grain sizes requested in the parameter
          file. MCFOST does not extrapolate, ie all the grains smaller (larger)
          than the minimum (maximum) grain size provided in the FITS file will
          have the same spatial distribution as the minimum (maximum) grain
          size.

.. important:: If not provided, the gas is assumed to follow the density of the smallest grain size
               provided (you can however create an artificially small grain size for
               the gas if you wish).

There are 2 options to normalize the different spatial distribution of
the grains relative to each other. These 2 options can be selected using
the keyword ``read_n_a`` in the header of the 1st HDU of the FITS file.

If ``read_n_a`` is set to 0, the density structure does not need to be
normalized. MCFOST will use the dust mass, grain size power law
distribution and gas-to-dust mass ratio provided in the parameter file.
In other words, the fits file is only used to pass the spatial
distribution of each grain sizes and MCFOST will take care of all the
required normalization.

If there are more than 1 grain size, the values of the ``n_grains`` grain
sizes used in the FITS file in a second HDU with 1 dimension (NAXIS=1) :
NAXIS1 = n_grains (unit : micron) : ``grain_sizes(1:n_grains)``.
These grain sizes must be given in increasing order (and of course match
the order of the data cube).

If ``read_n_a`` is set to 1, the relative density in each grain size bin
needs passed to MCFOST and you must provide the grain size distribution
integrated over the whole disk volume. This is done by using a 3st HDU
where you pass an array of dimension (1:n_grains) providing the number of
grains in each grain size bin : dn(a)/da

.. note:: the total normalization of the distribution is not required as
          MCFOST will use the total dust mass as set in the parameter file. The
          array n_a is only used to define the relative number of the grains of
          each size.

The gas density can be passed by setting the keyword ``read_gas_density``
to 1. In that case, an extra keyword ``gas_to_dust`` must also be passed
to give the integrated gas to dust mass ratio.
An extra HDU must also be passed, giving the gas density with 3
dimensions: ``gas_density(1:n_rad, 1:nz, 1:n_az)``.
As for the dust density array, the array has 3 dimensions also in 2D,
but with n_az=1. The total gas mass is set by the mcfost parameter
file in any case. If the keyword ``read_gas_density`` is set to 0 or
absent, mcfost will assume that the gas has the same spatial
distribution as the smallest grains in the FITS file, and will use the
gas-to-dust mass ratio provided in the parameter file.

.. note:: some parameters, such as ``read_gas_density`` are longer than 8 characters, which is the fits standard. They are however handled correctly by mcfost. In general, you can put a `HIERARCH` card in front of this parameter so that it respect fits standard. See for example `Astopy's documentation <http://docs.astropy.org/en/stable/io/fits/usage/headers.html#hierarch-cards>`_.


The gas velocity field, in Cartesian or cylindrical coordinates, can
also be passed by setting ``read_gas_velocity`` to 1 (Cartesian
coordinates) or 2 (cylindrical coordinates). The fits file must have
an extra HDU with 4 dimensions: ``gas_velocity(1:n_rad, 1:nz, 1:n_az, 3)``.
The last index correspond to ``v_x``, ``v_y``, ``v_z`` in Cartesian
coordinates, or ``v_r``, ``v_phi``, ``v_z`` in cylindrical
coordinates. The velocity is given in m/s.

.. note:: If the velocity is not passed via the fits interface, mcfost will assume
          that the velocity field is Keplerian (unless modified by command line
          options).


Unstructured set of points
--------------------------

mcfost can also read an arbitrary set of points, in which case mcfost will perform a Voronoi tesselation on the provided points (as for SPH calculations, for instance with the `-phantom` option)
The syntax is the same as above::

  $ mcfost <parameter_file> -density_file <your_file.fits.gz>

mcfost will automatically detect if the data is structured or unstructured.

The fits file needs to have at least 2 hdus:

* the first will contain the *mass* of each "particle" as a 1D array of size ``n_points``.

.. note:: because the volume of a Voronoi cell can be different from the initial volume in the hydro model, the mass and not the density of each particle is required.

* the 2nd hdu must be 2D with dimensions of ``3 x n_points``, and provide the x,y,z coordinates for all the data points



Phantom interface
#################

MCFOST can read both the standard and hdf5 phantom dump formats.
Typical syntax for reading a phantom file is:``mcfost <para_file> -phantom <phantom_dump>``
Any additional option can be used in conjonction with the `-phantom` option.

When reading a phantom dump, sections of the parameter file will be ignored, in particular any section describing the disc density, model grid and stellar properties (unless the stellar properties are forced).

Idefix, pluto, fargo3d and athena++ interfaces
##############################################

MCFOST can also read native dumps from idefix, pluto, fargo3d, athena++ if they are on a regular grid. The mcfost grid will be matched eactly to the hydro grid, which means mcfost can used the density in each cell, rather than the mass (as for a set of points or a SPH model). Syntax is the similar to the phantom case.

For idefix::

  $ mcfost <para_file> -idefix <idefix_vyk_file>

For pluto::

   $ mcfost <para_file> -pluto <pluto_dump>

For athena++::

   $ mcfost <para_file> -athena++ <athena_adhf_file>

For fargo3d::

   $ mcfost <para_file> -fargo3d <fargo3d_directory> <dump_number>




Scaling units
#############

Values in hydro simulations might be provided in code units rather than physical units. They can be rescaled with the options ``-scale_length_units <factor>`` and  ``-scale_mass_units <factor>``.




Important sections in parameter file
####################################

Here are some important sections to consider for post-processing hydro dump files.

Number of photons for temperature and image calculation
-------------------------------------------------------

You may need to use more photons than for a 2D model to compute the temperature. As a rule of thumb, you should use at least 100 more packets than SPH particles or grid cells, ie the section is the parameter file should be something similar to

::

 #Number of photon packages
   1.28e8                  nbr_photons_eq_th  : T computation
   1.28e3                  nbr_photons_lambda : SED computation
   1.28e6                  nbr_photons_image  : images computation



Turn off SED computation
------------------------

You probably don't want to compute an SED if you have a hydro file. Turning it off will save computtation time.

::

   T F T                   compute temperature?, compute sed?, use default wavelength grid for ouput ?


Turn off image symmetries
-------------------------

You don't want image symmetries if you are processing an hydrodynamical dump.

::

 #Symetries
   F                       image symmetry
   F                       central symmetry
   F                       axial symmetry (important only if N_phi > 1)

Can ignore the following, they come from the hydro dump
-------------------------------------------------------

::

 #Grid geometry and size
 #Disk physics
 #Number of zones : 1 zone = 1 density structure + corresponding grain properties

Density structure
-----------------

The density structure comes from the hydro dump, so can ignore "#Density structure" however you need
to set "gas-to-dust mass ratio" in your model does not have dust:

::

   1.e-3    100.           dust mass,  gas-to-dust mass ratio


.. note:: You can force the gas mass to the value in the parameter file using the option ``-force_Mgas``.


Star properties
---------------

By default, mcfost will use the properties (ie mass) of the sink particles to estimate the brightness temperature and luminosity.
mcfost will do so by interpolating through isochrones, assuming an age of 3Myr by default. The age can be change with the ``-age`` option.
The stellar properties can be fixed to the values in the parameter files with ``-fix_stars``. Note that in that case, the number of stars in the mcfost parameter file must matches the number of sink particles in the hydro dump.
