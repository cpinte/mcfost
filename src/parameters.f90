module parameters
!****************************
! general parameters
!****************************

  use mcfost_env, only : dp

  implicit none
  save

  real :: para_version
  real(kind=dp) :: simu_time = 0.0

  logical :: lpara, lstop_after_init
  integer :: step_index, etape_i, etape_f

  ! Number of photons launched
  logical :: ldust_transfer
  integer :: n_photons_loop, n_photons_eq_th, n_photons_lambda, n_photons_image, n_photons_spectrum
  real :: n_photons_lim = 1.e4 ! how many times more we would have received without the Disk
  integer :: nnfot1
  real :: n_photons_total
  real(kind=dp) :: E_paquet
  integer :: n_dif_max_eq_th = 100000 ! Max number of scattering allowed in thermal eq. calculation OUTDATED
  real :: tau_dark_zone_eq_th = 1500 !1500.   15000 pour benchmark at tau=1e6
  real :: tau_dark_zone_obs = 100 ! same as 1000 if 1500. above
  integer :: n_Stokes

  ! Number of angles to sample the phase function
  integer, parameter :: nang_scatt = 180  ! TODO : it bugs if it is not 180

  ! lvariable_dust = true if the dust properties are variable in each cell
  logical :: lvariable_dust, lmigration, lhydrostatic, ldust_sublimation
  integer :: settling_type ! 1 = Parametric, 2 = Dubrulle or 3 = Fromang

  logical :: lRE_LTE, lRE_nLTE, lnRE, lonly_LTE, lonly_nLTE
  logical :: loutput_J, loutput_J_step1, loutput_UV_field, lxJ_abs, lxJ_abs_step1

  ! Scattering calculation method: to choose to optimize memory size and CPU time
  ! 0 -> automatique
  ! 1 -> draw scattering grain size + Mueller matrix per grain
  ! 2 -> average Mueller matrix per cell (benchmark)
  integer :: scattering_method, scattering_method0
  logical :: lscattering_method1

  ! Mie or HG theory
  integer :: aniso_method ! 1 = full phase function, 2 = HG
  logical :: lmethod_aniso1

  integer :: RT_sed_method ! cf routine dust_map for definition

  ! Thermal emission steps
  logical :: ltemp, lsed, lsed_complete, l_em_disk_image, lchauff_int, lextra_heating, lno_internal_energy
  character(len=512), dimension(:), allocatable :: indices
  character(len=512) :: tab_wavelength

  !gas transfer here
  logical :: ldust_gas

  ! Emission moleculaire
  logical :: lemission_mol,  lpop, lprecise_pop, lmol_LTE, ldust_mol, lonly_top, lonly_bottom

  ! Atomic line radiative transfer
  logical :: lexit_after_nonlte_loop, lstop_after_jnu
  logical :: lemission_atom, lelectron_scattering, lforce_lte
  logical :: ldissolve, loutput_rates, lzeeman_polarisation, ldust_atom
  integer :: N_rayons_mc, istep_start, istep_end

  !HEALpix
  integer :: healpix_lorder, healpix_lmin, healpix_lmax !lmin and lmax not yet (for local evaluation)

  logical :: lcheckpoint, lsafe_stop
  !Convergence relative errors
  real :: dpops_max_error, art_hv, safe_stop_time
  integer :: checkpoint_period

  !Ng's acceleration
  logical :: lng_acceleration
  integer :: Ng_Norder, Ng_Nperiod

  !electron density
  logical :: lsolve_for_ne, lescape_prob
  integer :: ndelay_iterate_ne, n_iterate_ne !0 means once SEE is solved. Otherwise, > 1, iterated every n_iterate_ne during the nlte_loop

  logical :: lmhd_voronoi
  integer :: limit_mem !if limit_mem == 0: try to keep a maximum of quantities on ram (faster but very ram-consuming)
  					   !                    currently it means that background and continua are stored on the full wavelength grid for non-LTE transfer.
  					   !if limit_mem == 1: the continua are stored on a small frequency grid and interpolated locally on the full grid.
  					   !					This approach is faster than computing the continua for each wavelength point, and is relatively cheap in ram.
  					   !					-> good trade-off between 0 and 2.
  					   !if limit_mem == 2: everything is computed locally on the full grid. Slow but cheap in memory.

  ! Image decomposition
  logical :: lsepar_contrib, lsepar_pola, lonly_capt_interet, lsepar_ori
  integer :: N_type_flux
  ! the fluxes are I, (Q,U,V), (star, scatt, disk th, disk th scatt.)

  ! rotation of the Disk plane in deg., counter-clockwise.
  real(kind=dp) :: ang_disque

  ! Production of symmetric images
  ! The symmetry is performed before choosing the pixels
  ! Is the system centrosymmetric?
  ! does the system have axial symmetry (only counts if N_phi > 1)
  logical :: l_sym_ima, l_sym_centrale, l_sym_axiale

  ! maps parameters
  integer :: N_thet, N_incl, N_phi, capt_interet, delta_capt, capt_inf, capt_sup, capt_debut, capt_fin
  integer ::  npix_x, npix_y, npix_x_save, npix_y_save
  real :: angle_interet, zoom, tau_seuil, wl_seuil, image_offset_centre(3)

  real  :: cutoff = 7.0

  !must be initialized to 0
  integer(kind=8) :: mem_alloc_tot !total memory allocated dynamically or not in bytes

  ! Density grid resolution
  ! Number of cells in the r direction (log sampling)
  integer :: grid_type ! 1 = cylindrical, 2 = spherical
  integer :: n_rad, n_rad_in  ! subdivision of the first cell
  ! Number of vertical layers (= stratifications)
  integer :: nz, p_n_rad, p_nz, p_n_az, p_n_lambda_pos, p_n_lambda_grain
  ! Number of azimuthal cells
  integer :: n_az, j_start, pj_start
  ! Total number of cells
  integer :: n_cells, nrz, p_n_cells
  integer, target :: icell1 = 1
  logical :: lregular_theta
  real :: theta_max, theta_mask_max

  logical :: letape_th, limg, lorigine, laggregate, lFresnel, lFresnel_per_size, l3D, lremove, lwarp, lcavity, ltilt, lwall
  logical :: lopacite_only, lseed, ldust_prop, ldisk_struct, lwrite_velocity, loptical_depth_to_cell, ltau_map, lreemission_stats
  logical :: lapprox_diffusion, lcylindrical, lspherical, llinear_rgrid, lVoronoi, is_there_disk, lno_backup
  logical :: laverage_grain_size, lisotropic, lno_scattering, lqsca_equal_qabs, lonly_diff_approx, lforce_diff_approx
  logical :: ldensity_file, lsigma_file, lvelocity_file, lphantom_file, lphantom_multi, lphantom_avg
  logical :: lgadget2_file, llimits_file, lforce_SPH_amin, lforce_SPH_amax, lmcfost_lib
  logical :: lweight_emission, lcorrect_density, lProDiMo2mcfost, lProDiMo2mcfost_test, lastrochem, lML
  logical :: lspot, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium, lchange_Tmax_PAH, lISM_heating, lcasa, lJy, lforce_Mgas
  integer :: ISR_model ! 0 : no ISM radiation field, 1 : ProDiMo, 2 : Bate & Keto
  integer :: vfield_coord ! 1 : Cartesian, 2 : cylindrical, 3 : spherical

  logical :: lfargo3d, lathena, lidefix, lpluto, lsphere_model, lmodel_1d, lheader_only !future lsymspheric

  ! benchmarks
  logical :: lbenchmark_Pascucci, lbenchmark_vanZadelhoff1, lbenchmark_vanZadelhoff2, lDutrey94, lHH30mol
  logical :: lbenchmark_water1, lbenchmark_water2, lbenchmark_water3, lbenchmark_SHG, lMathis_field
  real :: Mathis_field

  ! Prodimo
  logical :: lprodimo, lprodimo_input_dir, lforce_ProDiMo_PAH

  logical :: lSeb_Charnoz, lread_Seb_Charnoz, lread_Seb_Charnoz2, lread_Misselt, lread_DustEM
  logical :: lread_grain_size_distrib, lphase_function_file,ltau_surface, lflux_fraction_surface
  logical :: lwrite_column_density, lwrite_mol_column_density, lwrite_abundance
  character(len=8) :: sflux_fraction, stau_surface
  real(kind=dp) :: flux_fraction
  real :: tau_surface

  ! Phantom
  logical :: ldudt_implicit, lscale_length_units, lscale_mass_units, lignore_dust
  logical :: ldelete_Hill_sphere, lrandomize_Voronoi, lrandomize_azimuth, lrandomize_gap, lrandomize_outside_gap, lcentre_on_sink
  logical :: lmask_inside_rsph, ldelete_outside_rsph, ldelete_above_theta, lmask_outside_rsph, lmask_above_theta, lexpand_z
  real(kind=dp) :: ufac_implicit,scale_length_units_factor,scale_mass_units_factor,correct_density_factor_elongated_cells
  real(kind=dp) :: SPH_amin, SPH_amax, fluffyness, gap_factor, rsph_min, rsph_max, rsph_mask_max, expand_z_factor
  logical :: lupdate_velocities, lno_vr, lno_vz, lvphi_Kep, lfluffy, lnot_random_Voronoi, lignore_sink
  integer :: isink_centre

  ! Disk parameters
  real :: distance ! Distance of the Disk in pc
  real(kind=dp) :: map_size

  ! Polarisation
  logical :: loverwrite_s12
  real :: Pmax

  integer :: n_zones, n_regions

  type disk_zone_type
     real(kind=dp) :: Rin, Rmin, Rc, Rout, Rmax, Rref, edge, exp_beta, surf
     real(kind=dp) :: moins_gamma_exp, sclht, diskmass, gas_to_dust, vert_exponent
     integer :: geometry ! 1=disk, 2=tappered-disk, 3=envelope
     integer :: region
  end type disk_zone_type

  type disk_region_type
     integer :: n_zones
     real(kind=dp) :: Rmin, Rmax
     integer :: iRmin, iRmax
     integer, dimension(:), allocatable :: zones
  end type disk_region_type

  type cavity_type
     real(kind=dp) ::  exp_beta, sclht, Rref
  end type cavity_type

  type(disk_zone_type), dimension(:), allocatable, target :: disk_zone
  type(disk_region_type), dimension(:), allocatable, target :: regions
  type(cavity_type) :: cavity

  real(kind=dp) :: Rmin, Rmax, Rmax2, diskmass, correct_Rsub

  real :: exp_strat, a_strat
  real :: alpha

  ! Analytical description of the puffed-up inner rim
  logical :: lpuffed_rim
  real :: puffed_rim_h, puffed_rim_r, puffed_rim_delta_r

  real :: z_warp, tilt_angle

  ! SPH
  real :: SPH_keep_particles, planet_az, delta_planet_az
  logical :: lplanet_az, lfix_star, lcorrect_density_elongated_cells, lturn_off_planets, lturn_off_Lacc, lforce_Mdot
  integer :: which_planet, idelta_planet_az

  logical :: lgap_Gaussian
  real :: f_gap_Gaussian, r_gap_Gaussian, sigma_gap_Gaussian

  ! Local density correction (in a ring)
  real :: correct_density_factor, correct_density_Rin, correct_density_Rout

  ! Vertical scaling of the envelope
  real :: z_scaling_env

  character(len=512) :: sigma_file, grain_size_file, limits_file
  character(len=512), dimension(:), allocatable :: density_files
  integer :: n_phantom_files

  ! Stars
  type star_type
     ! todo : indicate all units
     real :: T, M, fUV, slope_UV, othin_sublimation_radius
     real :: Mdot ! Msun
     real(kind=dp) :: r, x,y,z, vx,vy,vz
     logical :: lb_body, out_model, find_spectrum, force_Mdot
     character(len=512) :: spectrum
     integer :: icell
  end type star_type
  logical :: laccretion_shock

  integer :: n_stars
  type(star_type), dimension(:), allocatable :: star
  logical :: lstar_bb


  ! Spot
  real :: T_spot, surf_fraction_spot, theta_spot, phi_spot

  ! Limb darkening
  logical :: llimb_darkening
  character(len=512) :: limb_darkening_file
  real, dimension(:), allocatable :: mu_limb_darkening, limb_darkening, pola_limb_darkening

  ! ISM radiation following ProDiMo's definitions
  real, parameter :: Wdil =  9.85357e-17
  real, parameter :: TCmb = 2.73
  real, parameter :: T_ISM_stars = 20000.
  real :: chi_ISM = 1.0

  character(len=8) :: system_age

  logical :: lscatt_ray_tracing, lscatt_ray_tracing1, lscatt_ray_tracing2, loutput_mc
  ! ray-tracing 1: saves the diffuse radiation field
  ! ray-tracing 2: saves the specific intensity

  integer :: n_pop


  real :: T_max, T_min, Tmax_PAH ! Sublimation Temp and cloud Temp
  integer  :: n_T

  real, dimension(:), allocatable :: s11_file

  ! inclinaisons
  real :: RT_imin, RT_imax, RT_az_min, RT_az_max
  integer ::  RT_n_incl, RT_n_az
  logical :: lRT_i_centered

  integer :: nLevels
  real(kind=dp) :: profile_width

  type fargo3d_model
     integer :: nx, ny, nz, realtype
     real(kind=dp) :: xmin,xmax, ymin,ymax, zmin,zmax
     logical :: log_spacing, corrotating_frame

     real(kind=dp) :: dt, aspect_ratio, nu, gamma, cs
     character(len=128) :: dir, id, planetconfig
  end type fargo3d_model

  type(fargo3d_model) :: fargo3d

  type athena_model
     integer :: nx1, nx2, nx3
     real(kind=dp) :: x1_min,x1_max, x2_min,x2_max, x3_min,x3_max
     logical :: log_spacing, corrotating_frame

     real(kind=dp) :: time
     character(len=128) :: filename
  end type athena_model

  type(athena_model) :: athena

  type idefix_model
     integer :: nx1, nx2, nx3, iunit, position, geometry, id
     integer, dimension(3) :: dimensions
     real(kind=dp) :: x1_min,x1_max, x2_min,x2_max, x3_min,x3_max
     logical :: log_spacing, corrotating_frame

     real :: time
     character(len=128) :: filename
     character(len=32) :: origin

     real, dimension(:), allocatable :: x1, x2, x3
  end type idefix_model

  type(idefix_model) :: idefix

  type pluto_model
     integer :: nx1, nx2, nx3, iunit, position, geometry
     integer, dimension(3) :: dimensions
     real(kind=dp) :: x1_min,x1_max, x2_min,x2_max, x3_min,x3_max
     logical :: log_spacing, corrotating_frame

     real :: time
     character(len=128) :: dir, id
     character(len=32) :: origin

     real, dimension(:), allocatable :: x1, x2, x3
  end type pluto_model

  type(pluto_model) :: pluto

  type symspheric_model
   integer :: nr
   real(kind=dp), dimension(:), allocatable :: r, T, rho, vt, ne
   real(kind=dp), dimension(:,:), allocatable :: v
   integer, dimension(:), allocatable :: iz

   integer :: Ncorona
   logical :: lcoronal_illumination
   real(kind=dp), allocatable :: I_coro(:,:), x_coro(:), m(:)

   real(kind=dp) :: rbot, rtop, rstar, s, E_corona

  end type symspheric_model

  type(symspheric_model) :: atmos_1d

end module parameters
