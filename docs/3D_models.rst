
Running a 3D model
==================

MCFOST is completely 3D and can be used with any density structure.
Interfaces for hydrodynamics codes can be built on demand if there are
specific needs. There is also a default interface using FITS file to
input MCFOST with any arbitrary density structure. This interface is
likely to fulfill most of the needs.

It can be used with the command::

$ mcfost <parameter_file> -density_file <your_density_file.fits.gz> -3D (+ any other option)


.. important:: When using the `-density_file` or `-sigma_file` options, the number of zones must be set to 1 in mcfost


The density array in the FITS file must have 4 dimensions : ``density(1:n_rad, 1:nz, 1:n_az, 1:n_grains)``
(it then has n_rad x nz x n_az x n_grains elements).

The option ``-density_file`` also works in the 2D version of the code, ie
without ``-3D``. In that case the density file still requires 4 dimensions,
but with n_az = 1 : ``density(1:n_rad, 1:nz, 1, 1:n_grains)``

.. note:: the spatial grid is not currently passed to MCFOST in the
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
coordinates. The velocity is given in km/s.

.. note:: If the velocity is not passed via the fits interface, mcfost will assume
          that the velocity field is Keplerian (unless modified by command line
          options).


An alternative option `-sigma_file` works in a similar way with the vertical dimension and passes only the surface density. mcfost use the scale height and flaring index provided in the parameter file to reconstruct the 3D density structure.
