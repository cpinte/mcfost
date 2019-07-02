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

File organisation
-----------------
``-seed <seed>``: modifies seed for random number generator ; results stored in ``seed=XXX`` directory

``-root_dir <root_dir>``: results stored in ``root_dir`` directory

``-no_backup`` : stops if directory ``data_XX`` already exists without attempting to backup existing directory

``-prodimo_input_dir <dir>``: input directory for ProDiMo

Options related to images
-------------------------

``-zoom <zoom>``: overrides value in parameter file

``-resol <nx> <ny>``: overrides value in parameter file

``-PA``: override value in parameter file

``-only_scatt`` : ignore dust thermal emission

``-casa`` : write an image ready for CASA

``-nphot_img`` : overwrite the value in the parameter file

``-rt``: use ray-tracing method to compute images or SEDs (on by default)

``-rt1`` or ``-rt2``: use ray-tracing method and force ray-tracing method

``-no-rt``: do not output the ray-tracing results

``-mc``:  keep Monte-Carlo output in ray-tracing mode

``-n_MC_bins <n_inclinations> <n_azimuth>`` (default : 10 1)

``-planet_az <angle>`` [deg]: adjust the model azimuth so that the planet is at
desired azimuth in the map

Options related to temperature equilibrium
``-no_T``: skip temperature calculations, force ltemp to F

``-diff_approx``: enforce computation of T structure with diff approx.

``-no_diff_approx``: compute T structure with only MC method

``-only_diff_approx``: only compute the diffusion approx

``-tau_dark_zone_obs <tau_dark_zone>`` (default : 100)

``-tau_dark_zone_eq_th <tau_dark_zone>`` (default : 1500)

``-origin``: save origin of packets received the interest bin

``-rs (remove specie) <specie_number> <Temperature>``

``-reemission_stats``

``-weight_emission``: weight emission towards disk surface

``-force_PAH_equilibrium``: mainly for testing purposes

``-force_PAH_out_equilibrium``: mainly for testing purposes

``-Tmax_PAH <T>``: changes the maximum temperature allowed for PAH (default: 2500)

``-ISM_heating``: includes heating by ISM radiation

``-chi_ISM <chi>``: changes the chi of ISM radiation (default: 1)

``-no_internal_energy``: ignoring internal energy in Phantom dumps

Disk structure
--------------

``-disk_struct``: computes the density structure and stops:
gas_density.fits.gz and dust_density.fits.gz -> density map
grid.fits.gz -> radii and height in the grid
volume.fits.gz -> volume per cell at each radius

``-r_subdivide <radius>``: forces cell subdivision elsewhere
than at inner radius

``-3D``: 3D geometrical grid

``-warp`` <h_warp>: defined at reference radius

``-tilt``: <angle> [degrees]

``-cavity <h0> <r0> <flaring exponent>``: carves a cavity in the the disc. Defines a surface, and empties the density above this surface.

``-output_J``: outputs the radiation field in each cell

``-output_UV_field``: output the UV radation field between 912 and 2000 Angströms in each cell. units: Habings

``-puffed_up_rim  <h_rim / h0> <r_rim> <delta_r>``: add a puffed up inner rim to the disc, with an increased scale height by a factor  <h rim / h0> and up to a radius r. The width over the scale goes back to normal is delta_r

The updated scale height is

.. math::
   h(r) = h_0(r) \times \left(1.0 + \frac{\frac{h_{rim}}{h_0} - 1.0}{\exp(\frac{r- r_{rim}}{\Delta r}) + 1.0} \right)


``-density_file or -df <density_file>``: reads a fits file with the density (gas + dust + velocity, see section "Running a 3D model")

``-sigma_file or -sigma <surface_density_file>``: reads a fits file with the surface density

``-correct_density <factor> <Rmin> <Rmax>``: applies a correcting factir to the density between Rmin and Rmax.  Can be used to generate rings or gaps.

``-gap <depth> <R> <sigma>``: creates a Gaussian gap in the disc [depth is between 0 and 1, R and Sigma in au]

``-Seb_F <number>``: select the dust diffusion method frollowing Sebastian Fromang's prescriptions.  1 = gaussian, 2 = cst diffusion coeff

``-cutoff <number>``: upper limit of the grid [scale height] default = 7

``-n_rad``: overwrite value in parameter file

``-nz``: overwrite value in parameter file

``-z_scaling_env <scaling_factor>``: scale a spherical envelope along the z-axis


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
Isochrones are found in ``$MCFOST_UTILS/Stellar_Polarization/Siess``.

Dust properties
---------------

``-dust_prop``: computes opacity, albedo, asymmetry parameter,
polarizability and saves results in ``data_dust``

``-op <wavelength>`` [microns]: computes dust properties at
specified wavelength and stops

``-aggregate <GMM_input_file> <GMM_output_file>``

``-optical_depth_map``, ``-od``: generates a map of integrated optical depth
along radial and vertical directions and stops;
results stored in ``optical_depth_map.fits.gz``

``-average_grain_size``: computes average grain size in each cell,
weighted by their geometrical cross-section;
results stored in ``average_grain_size.fits.gz``

``-HG``: uses an Heynyey-Greenstein function

``-force_HG <g>``: uses an Heynyey-Greenstein function and forces the g value

``-isotropic``: forces isotropic scattering

``-no_scattering``: forces albedo = 0

``-qsca=qabs``: forces albedo = 0.5

``-phase-function <s11.fits>``: uses a tabulated phase function (rt2 only)

``-tau=1_surface``: when computing an image, generates a fits file with the coordinates x,y,z of the tau=1 surface for each pixel. The dimension of the fits file is the same as the image fits file + an extra dimension with the 3 values x,y and z.

Options related to molecular emission
-------------------------------------

``-freeze-out <T>``: freeze-out molecules (abundance = 0) in regions below this temperature

``-freeze-out_depletion <relative depletion>``: between 0 and 1, depletion factor in the freeze-out region (default: 0)

``-photo-dissociation``: photo-dissociate molecules at high-UV (see Pinte et al 2018)

``-photo-desorption``: photo-desorb molecules at high-UV (see Pinte et al 2018)

``-prodimo``: generates a a fits file forProDiMo.fits.gz with the required information to run a thermo-chemical model with ProDiMo

``-prodimo_fPAH``: force a fPAH value for ProDiMo

``-only_top``: molecular emssion from the top half of the disk

``-only_bottom``: molecular emssion from the bottom half of the disk

``-correct_Tgas <factor>``: applies a factor to the gas temperature

``-chi_infall <value>``: v_infall/v_kepler

``-cylindrical_rotation``: forces Keplerian velocity of independent of z
