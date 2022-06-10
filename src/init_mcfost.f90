module init_mcfost

  use parametres
  use naleat
  use grains, only : aggregate_file, mueller_aggregate_file
  use density, only : species_removed, T_rm
  use molecular_emission
  !$ use omp_lib
  use benchmarks
  use read_params
  use input, only : Tfile, lect_lambda, read_phase_function, read_molecules_names
  use ProdiMo
  use utils
  use read_fargo3d, only : read_fargo3d_parameters
  use read_athena, only : read_athena_parameters
  use read1d_models, only : read_grid_1d

  implicit none

contains

subroutine set_default_variables()

  ! Pour code sequentiel
  nb_proc=1 ; lpara=.false.

  n_zones=1
  lmcfost_lib = .false.

  ! Pour code parallel
  !$omp parallel default(none) &
  !$omp shared(nb_proc, lpara)
  !$ nb_proc=omp_get_num_threads() ; lpara=.true.
  !$omp end parallel

  ldisk_struct = .false.
  ldust_prop = .false.
  lstop_after_init= .false.
  lwall=.false.
  lpah=.false.
  limg=.false.
  lorigine=.false.
  laggregate=.false.
  l3D=.false.
  lopacite_only=.false.
  lseed=.false.
  loptical_depth_map=.false.
  lreemission_stats=.false.
  n_az = 1
  root_dir = "."
  seed_dir = "."
  lapprox_diffusion = .true.
  lonly_diff_approx = .false.
  lbenchmark_Pascucci = .false.
  lbenchmark_vanZadelhoff1 = .false.
  lbenchmark_vanZadelhoff2 = .false.
  lbenchmark_water1 = .false.
  lbenchmark_water2 = .false.
  lbenchmark_water3 = .false.
  lbenchmark_SHG = .false.
  lDutrey94 = .false.
  lHH30mol = .false.
  lemission_mol=.false.
  lcheckpoint = .false.
  checkpoint_period = 15
  !HEALpix
  healpix_lorder = 1
  healpix_lmin = 0
  healpix_lmax = 7 !6 !5
  ! Atomic lines Radiative Transfer (AL-RT)
  lsafe_stop = .false.
  safe_stop_time = 155520.0!1.8days in seconds, default
  lemission_atom = .false.
  lelectron_scattering = .false.
  lstop_after_jnu = .false.
  lsolve_for_ne = .false.
  n_iterate_ne = -1 !negative means never updated after/during non-LTE loop.
  ndelay_iterate_ne = 0
  lmodel_1d = .false.
  lmodel_ascii = .false.
  lmhd_voronoi = .false.
  lzeeman_polarisation = .false.
  lforce_lte = .false.
  ldissolve = .false.
  lng_acceleration = .false.
  Ng_Norder = 0
  Ng_Nperiod = 5
  istep_start = 1
  ! AL-RT
  laccurate_integ = .false.
  N_rayons_mc = 100
  loutput_rates = .false.
  !
  ! Max relative error in transfer. ATM only atomic line transfer
  dpops_max_error = 1e-1
  dpops_sub_max_error = 1e-2
  art_hv = 0.0 !default is frac * min(vD)
  art_hv_nlte = 0.0 !default is frac * min(vD)
  !
  lpuffed_rim = .false.
  lno_backup = .false.
  loutput_UV_field = .false.
  loutput_J = .false.
  loutput_J_step1 = .false.
  laverage_grain_size = .false.
  lfreeze_out = .false.
  lphoto_dissociation = .false.
  lphoto_desorption = .false.
  l_em_disk_image = .true.
  lisotropic = .false.
  lno_scattering = .false.
  lqsca_equal_qabs = .false.
  lscatt_ray_tracing=.true.
  lscatt_ray_tracing1=.false.
  lscatt_ray_tracing2=.false.
  loutput_mc=.false.
  ldensity_file=.false.
  lvelocity_file=.false.
  lphantom_file=.false.
  lphantom_multi = .false.
  lphantom_avg = .false.
  lforce_Mgas = .false.
  lforce_SPH_amin = .false.
  lforce_SPH_amax = .false.
  lascii_SPH_file = .false.
  lgadget2_file=.false.
  llimits_file = .false.
  lsigma_file = .false.
  lweight_emission=.false.
  lprodimo=.false.
  lProDiMo_input_dir=.false.
  lProDiMo2mcfost=.false.
  lProDiMo2mcfost_test = .false.
  lforce_ProDiMo_PAH = .false.
  lforce_diff_approx = .false.
  lgap_Gaussian=.false.
  lspot = .false.
  lSeb_Charnoz = .false.
  lread_Seb_Charnoz = .false.
  lread_Seb_Charnoz2 = .false.
  lonly_bottom = .false.
  lonly_top = .false.
  lcorrect_Tgas = .false.
  lcorrect_density=.false.
  lremove = .false.
  lforce_PAH_equilibrium=.false.
  lforce_PAH_out_equilibrium=.false.
  lread_grain_size_distrib=.false.
  lMathis_field = .false.
  lchange_Tmax_PAH=.false.
  lISM_heating = .false. ; ISR_model = 0
  llimb_darkening = .false.
  lVoronoi = .false.
  lcavity = .false.
  freeze_out_depletion = 0.
  lmono0 = .false.
  lmono = .false.
  lextra_heating = .false.
  lno_internal_energy = .false.
  ldudt_implicit = .false.
  linfall = .false.
  lcylindrical_rotation = .false.
  lforce_HG = .false.
  lphase_function_file = .false.
  ltau1_surface=.false.
  lcasa=.false.
  lplanet_az = .false.
  which_planet = 0
  lML = .false.
  lcorrect_density_elongated_cells=.false.
  lfix_star = .false.
  lscale_length_units = .false.
  lscale_mass_units = .false.
  lignore_dust = .false.
  lupdate_velocities = .false.
  lno_vr = .false.
  lno_vz = .false.
  lvphi_Kep = .false.
  lfluffy = .false.
  ldelete_hill_sphere = .false.
  ldelete_inside_rsph = .false.
  ldelete_outside_rsph = .false.
  lrandomize_Voronoi = .false.
  lrandomize_azimuth = .false.
  lrandomize_gap = .false.
  lrandomize_outside_gap = .false.
  lcentre_on_sink = .false.
  lwrite_column_density = .false.
  lwrite_mol_column_density = .false.
  lturn_off_planets = .false.
  lturn_off_Lacc = .false.
  lforce_Mdot = .false.
  lregular_theta = .false.
  theta_max = 0.5*pi

  tmp_dir = "./"

  ! Geometrie Grille
  z_scaling_env = 1.0

  ! Methodes par defaut
  RT_sed_method = 1

  system_age = "3Myr"

  SPH_keep_particles = 0.999

  vfield_coord = 0

  return

end subroutine set_default_variables

!**********************************************

subroutine get_mcfost_utils_dir()

  integer :: i, n_dir

 ! Test if MCFOST_UTILS is defined
  call get_environment_variable('MCFOST_UTILS',mcfost_utils)
  if (mcfost_utils == "") call error("environnement variable MCFOST_UTILS is not defined.")
  call get_environment_variable('MY_MCFOST_UTILS',my_mcfost_utils)

  ! Directories to search (ordered)
  if (my_mcfost_utils == "") then
     n_dir = 2
     allocate(search_dir(n_dir))
     search_dir(1) = "." ; search_dir(2) = mcfost_utils ;
  else
     n_dir = 3
     allocate(search_dir(n_dir))
     search_dir(1) = "." ; search_dir(2) = my_mcfost_utils ; search_dir(3) = mcfost_utils ;
  endif

  allocate(dust_dir(n_dir),mol_dir(n_dir),star_dir(n_dir),lambda_dir(n_dir),ML_dir(n_dir))
  dust_dir(1) = "./" ; mol_dir(1) = "./" ; star_dir(1) = "./" ; lambda_dir(1) = "./" ; ML_dir(1) = "./"
  do i=2,n_dir
     dust_dir(i)   = trim(search_dir(i))//"/Dust/"
     mol_dir(i)    = trim(search_dir(i))//"/Molecules/"
     star_dir(i)   = trim(search_dir(i))//"/Stellar_Spectra/"
     lambda_dir(i) = trim(search_dir(i))//"/Lambda/"
     ML_dir(i) = trim(search_dir(i))//"/ML/"
  enddo

  return

end subroutine get_mcfost_utils_dir

!**********************************************

subroutine initialisation_mcfost()

  implicit none

  integer :: ios, nbr_arg, i_arg, nx, ny, syst_status, imol, i
  integer :: current_date, update_date, mcfost_auto_update, ntheta, nazimuth, ilen
  real(kind=dp) :: wvl
  real :: opt_zoom, utils_version, PA

  character(len=512) :: cmd, s, str_seed, para, base_para
  character(len=4) :: n_chiffres
  character(len=128)  :: fmt1, fargo3d_dir, fargo3d_id, athena_file

  logical :: lresol, lMC_bins, lPA, lzoom, lmc, lHG, lonly_scatt, lupdate, lno_T, lno_SED, lpola, lstar_bb

  real :: nphot_img = 0.0, n_rad_opt = 0, nz_opt = 0, n_T_opt = 0

  integer, parameter :: nmax_stars = 100
  real, dimension(nmax_stars) :: star_Mdot
  logical, dimension(nmax_stars) :: star_force_Mdot
  real :: Mdot
  integer :: istar_Mdot


  write(*,*) "You are running MCFOST "//trim(mcfost_release)
  write(*,*) "Git SHA = ", sha_id

  ! Local logical variables
  lzoom=.false.
  lresol=.false.
  lPA = .false.
  lmc = .false.
  lMC_bins = .false.
  lHG = .false.
  lonly_scatt = .false.
  lno_T = .false.
  lno_SED = .false.
  lpola = .false.
  lstar_bb = .false.
  star_force_Mdot(:) = .false.
  star_Mdot(:) = 0.0

  ! Global logical variables
  call set_default_variables()

  ! Setting up web-server
  call get_environment_variable('MCFOST_WEB_SERVER',s) ; if  (s/="") web_server = s
  webpage = trim(web_server)//trim(webpage)
  utils_webpage = trim(web_server)//trim(utils_webpage)

  ! Looking for the mcfost utils directory
  call get_mcfost_utils_dir()

  ! Do we need to search for an update ?
  mcfost_auto_update = 7
  call get_environment_variable('MCFOST_AUTO_UPDATE',s)
  if (s/="") read(s,*) mcfost_auto_update

  if (mcfost_auto_update > 0) then
     cmd = 'date +%s > date.tmp' ;   call appel_syst(cmd, syst_status)
     open(unit=1,file="date.tmp",status="old")
     read(1,*) current_date
     close(unit=1,status="delete")

     if (is_file(trim(mcfost_utils)//"/.last_update")) then
        open(unit=1,file=trim(mcfost_utils)//"/.last_update",status="old")
        read(1,*) update_date
        close(unit=1)
     else
        update_date = 0 ! if the file .last_update does not exist, I assume mcfost is old
     endif

     if ( (current_date-update_date) > mcfost_auto_update*86400) then
        write(n_chiffres,fmt="(i4)") ceiling(log10(1.0*mcfost_auto_update))
        fmt1 = '(" Your version of mcfost is more than ", i'//adjustl(trim(n_chiffres))//'," days old")'
        write(*,fmt=fmt1) mcfost_auto_update
        write(*,*) "Checking for update ..."
        lupdate = mcfost_update(.false.,.false., mcfost_auto_update)
        if (lupdate) then ! On redemarre mcfost avec la meme ligne de commande
           write(*,*) "Restarting MCFOST ..."
           write(*,*) ""
           call appel_syst(cmd_opt,syst_status)
           if (syst_status /=0) call error("MCFOST did not manage to restart")
        endif ! lupdate
     endif
  endif

  ! Ligne de commande
  call get_command(cmd_opt)

  ! Nbre d'arguments
  nbr_arg = command_argument_count()
  if (nbr_arg < 1) call display_help()

  call get_command_argument(1,para)

  ! Basic options
  if (para(1:1)=="-") then
     if (para(2:2)=="v") then ! mcfost version
        call mcfost_v()
     else if (para(2:8)=="history") then ! mcfost history
        call mcfost_history()
     else if (para(2:2)=="h") then
        call display_help()
     else if (para(2:6)=="setup") then ! download the utils and para file the 1st time the code is used
        call mcfost_setup()
     else if (para(2:9)=="get_para") then ! download current reference file
        call mcfost_get_ref_para()
     else if (para(2:14)=="get_yorick") then ! force update utils
        call mcfost_get_yorick()
     else if (para(2:4)=="url") then
        write(*,*) trim(webpage)//"linux/mcfost.tgz"
        write(*,*) trim(webpage)//"mac/mcfost.tgz"
        call exit(0)
     else  if (para(2:13)=="update_utils") then ! update utils
        call update_utils(.false.)
     else if (para(2:14)=="fupdate_utils") then ! force update utils
        call update_utils(.true.)
     else if (para(2:2)=="u") then ! update binary
        lupdate = mcfost_update(.false.,.true.)
     else if (para(2:3)=="fu") then ! force update binary
        lupdate =  mcfost_update(.true.,.true.)
     else
        call display_help()
     endif
     call exit(0)
  endif

  utils_version =  get_mcfost_utils_version()
  if (utils_version /= required_utils_version) then
     write(*,*) "ERROR: wrong version of the MCFOST_UTILS database"
     write(*,*) "Utils:", utils_version, "required:",required_utils_version
     write(*,*) "Please update with mcfost -update_utils"
     write(*,*) "Exiting."
     call exit(1)
  endif

  ilen = index(para,'/',back=.true.) ! last position of the '/' character
  base_para = para(ilen+1:)

  ! Les benchmarks
  if (para(1:8)=="Pascucci") then
     lbenchmark_Pascucci = .true.
     write(*,*) "Running Pascucci et al 2004 benchmark"
  else if (para(1:13)=="vanZadelhoff1") then
     lbenchmark_vanZadelhoff1 = .true.
     write(*,*) "Running van Zadelhoff et al 2002 benchmark 1"
  else if (para(1:13)=="vanZadelhoff2") then
     lbenchmark_vanZadelhoff2 = .true.
     write(*,*) "Running van Zadelhoff et al 2002 benchmark 2"
  else if (para(1:6)=="water1") then
     lbenchmark_water1 = .true.
     write(*,*) "Running water benchmark 1"
  else if (para(1:6)=="water2") then
     lbenchmark_water2 = .true.
     write(*,*) "Running water benchmark 2"
  else if (para(1:6)=="water3") then
     lbenchmark_water3 = .true.
     write(*,*) "Running water benchmark 3"
  else if (para(1:8)=="Dutrey94") then
     lDutrey94 = .true.
     write(*,*) "Running comparison with Dutrey et al 1994"
  else if (para(1:7)=="HH30mol") then
     lHH30mol = .true.
     write(*,*) "Running comparison with Pety et al 2006"
  endif

  i_arg=2

  ! Options ligne de commande
  do while (i_arg <= nbr_arg)
     call get_command_argument(i_arg,s)
     select case(trim(s))
     case("-v")
        call mcfost_v()
        call exit(0)
     case("-seed")
        lseed=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("seed needed")
        call get_command_argument(i_arg,str_seed)
        read(str_seed,*,iostat=ios) seed
        if (ios/=0) call error("seed needed")
        write(*,*) "Updating seed =", seed
        seed_dir = "seed="//trim(str_seed)
        i_arg = i_arg+1
     case("-img")
        limg=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("wavelength needed for -img. Error #1")
        call get_command_argument(i_arg,band)
        read(band,*,iostat=ios) wvl
        if (ios/=0) call error("wavelength needed for -img. Error #2")
        i_arg = i_arg+1
     case("-op")
        limg=.true.
        lopacite_only=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("wavelength needed for -op. Error #1")
        call get_command_argument(i_arg,band)
        read(band,*,iostat=ios) wvl
        if (ios/=0) call error("wavelength needed for -op. Error #2")
        i_arg = i_arg+1
     case("-origin")
        lorigine=.true.
        i_arg = i_arg+1
     case("-zoom")
        i_arg = i_arg+1
        lzoom = .true.
        if (i_arg > nbr_arg) call error("zoom needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) opt_zoom
        if (ios/=0) call error("zoom needed")
        i_arg = i_arg+1
     case("-pah")
        lpah=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("emmisivity needed")
        call get_command_argument(i_arg,model_pah)
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("grain type needed")
        call get_command_argument(i_arg,pah_grain)
        i_arg = i_arg+1
     case("-aggregate")
        laggregate=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("GMM input file needed")
        call get_command_argument(i_arg,aggregate_file)
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("GMM input file needed")
        call get_command_argument(i_arg,mueller_aggregate_file)
     case("-3D")
        l3D=.true.
        i_arg = i_arg+1
     case("-warp")
        if (.not.l3D) then
           write(*,*) "WARNING : forcing 3D mode"
           l3D=.true.
        endif
        lwarp=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("h_warp needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) z_warp
        i_arg= i_arg+1
     case("-tilt")
        if (.not.l3D) then
           call warning("forcing 3D mode")
           l3D=.true.
        endif
        call warning("tilt will only be applied to 1st zone")
        ltilt=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("tilt angle needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) tilt_angle
        i_arg= i_arg+1
     case("-rs")
        lremove=.true.
        i_arg= i_arg+1
        if (i_arg > nbr_arg) call error("species number needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) species_removed
        i_arg= i_arg+1
        if (i_arg > nbr_arg) call error("Temperature needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) T_rm
        i_arg= i_arg+1
     case("-resol")
        lresol=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("resolution needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) nx
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("resolution needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) ny
        i_arg= i_arg+1
     case("-n_MC_bins")
        lMC_bins=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("MC bins needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) ntheta
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("MC bins needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) nazimuth
        i_arg= i_arg+1
     case("-PA")
        lPA = .true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("PA needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) PA
        i_arg= i_arg+1
     case("-disk_struct","-output_density_grid","-ds")
        ldisk_struct=.true.
        i_arg = i_arg+1
        lstop_after_init= .true.
     case("+disk_struct","+ds")
        ldisk_struct=.true.
        i_arg = i_arg+1
        lstop_after_init= .false.
     case("-output_UV","-output_UV_field")
        loutput_UV_field=.true.
        i_arg = i_arg+1
     case("-output_J")
        loutput_J=.true.
        i_arg = i_arg+1
     case("-output_J1","-output_J_step1","-output_J_step_th")
        loutput_J_step1=.true.
        i_arg = i_arg+1
     case("-average_grain_size")
        laverage_grain_size=.true.
        i_arg = i_arg+1
     case("-killing_level")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("killing_level needed")
        call get_command_argument(i_arg,s)
        read(s,*) n_dif_max_eq_th
        i_arg = i_arg+1
     case("-tau_dark_zone_eq_th")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("tau_dark_zone needed")
        call get_command_argument(i_arg,s)
        read(s,*) tau_dark_zone_eq_th
        i_arg = i_arg+1
     case("-tau_dark_zone_obs")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("tau_dark_zone needed")
        call get_command_argument(i_arg,s)
        read(s,*) tau_dark_zone_obs
        i_arg = i_arg+1
     case("-root_dir")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("root_dir needed")
        call get_command_argument(i_arg,s)
        root_dir=trim(root_dir)//"/"//s
        i_arg = i_arg+1
     case("-tmp_dir")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) call error("root_dir needed")
        call get_command_argument(i_arg,tmp_dir)
        i_arg = i_arg+1
     case("-dust_prop")
        i_arg = i_arg+1
        ldust_prop=.true.
        lstop_after_init= .true.
     case("+dust_prop")
        i_arg = i_arg+1
        ldust_prop=.true.
        lstop_after_init= .false.
     case("-optical_depth_map","-od")
        i_arg = i_arg+1
        loptical_depth_map=.true.
     case("-reemission_stats")
        i_arg = i_arg+1
        lreemission_stats=.true.
     case("-no_diff_approx")
        i_arg = i_arg+1
        lapprox_diffusion = .false.
     case("-only_diff_approx")
        i_arg = i_arg+1
        lonly_diff_approx=.true.
     case("-diff_approx")
        i_arg = i_arg+1
        lforce_diff_approx=.true.
     case("-mol")
        i_arg = i_arg+1
        lemission_mol=.true.
     case("-safe_stop")
        i_arg = i_arg + 1
        lsafe_stop = .true.
     case("-safe_stop_time")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("time needed (safe_stop)")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) safe_stop_time
        safe_stop_time = safe_stop_time * 86400.!convert in sec
        i_arg= i_arg+1
     case("-checkpoint")
     	call error("option -checkpoint not yet")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("Period needed for checkpoint!")
        lcheckpoint = .true.
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) checkpoint_period !in iterations
        i_arg= i_arg+1
     case("-atom")
        ! Option to solve for the RTE for atoms
        i_arg = i_arg+1
        lemission_atom=.true.
     case("-output_rates")
        i_arg = i_arg + 1
        loutput_rates = .true.
     case("-electron_scatt") !force solving ne density even if provided in the model
     	call error("option -electron_scatt not yet")
        i_arg = i_arg + 1
        lelectron_scattering = .true.
     case("-solve_ne") !force solving ne density even if provided in the model
        i_arg = i_arg + 1
        lsolve_for_ne = .true.
     case("-iterate_ne") !Iterate electron density in the NLTE loop
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("Ne period needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) n_iterate_ne
        i_arg= i_arg+1
     case("-Ndelay_iterate_ne")!number of iteration before solving for electron density
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("Ne delay needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) ndelay_iterate_ne
        i_arg= i_arg+1
     case("-calc_jnu_atom")
     	call error("option -calc_jnu_atom not yet")
     	i_arg = i_arg + 1
     	lstop_after_jnu = .true.
     case("-puffed_up_rim")
        lpuffed_rim = .true.
        if (i_arg + 3 > nbr_arg) call error("rim parameters needed")
        i_arg = i_arg+1
        call get_command_argument(i_arg,s)
        read(s,*) puffed_rim_h
        i_arg = i_arg+1
        call get_command_argument(i_arg,s)
        read(s,*) puffed_rim_r
        i_arg = i_arg+1
        call get_command_argument(i_arg,s)
        read(s,*) puffed_rim_delta_r
        i_arg = i_arg+1
     case("-no_backup")
        lno_backup=.true.
        i_arg = i_arg + 1
     case("-Tfile")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,Tfile)
        i_arg = i_arg + 1
     case("-freeze_out","-freeze-out")
        i_arg = i_arg + 1
        lfreeze_out=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) T_freeze_out
     case("-freeze_out_depletion","-freeze-out_depletion")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) freeze_out_depletion
     case("-photodissociation","-photo_dissociation","-photo-dissociation")
        i_arg = i_arg + 1
        lphoto_dissociation=.true.
     case("-photodesorption","-photo_desorption","-photo-desorption")
        i_arg = i_arg + 1
        lphoto_desorption=.true.
     case("-isotropic","-iso")
        i_arg = i_arg + 1
        lisotropic=.true.
     case("-no_scattering","-no_scatt")
        i_arg = i_arg + 1
        lno_scattering=.true.
     case("-qsca=qabs")
        i_arg = i_arg + 1
        lqsca_equal_qabs=.true.
     case("-rt")
        i_arg = i_arg + 1
        lscatt_ray_tracing=.true.
        lscatt_ray_tracing1=.false.
        lscatt_ray_tracing2=.false.
        if (.not.lmc) loutput_mc=.false.
     case("-rt1")
        i_arg = i_arg + 1
        lscatt_ray_tracing=.true.
        lscatt_ray_tracing1=.true.
        lscatt_ray_tracing2=.false.
        if (.not.lmc) loutput_mc=.false.
     case("-rt2")
        i_arg = i_arg + 1
        lscatt_ray_tracing=.true.
        lscatt_ray_tracing1=.false.
        lscatt_ray_tracing2=.true.
        if (.not.lmc) loutput_mc=.false.
     case("-no-rt")
        i_arg = i_arg + 1
        lscatt_ray_tracing=.false.
        lscatt_ray_tracing1=.false.
        lscatt_ray_tracing2=.false.
        loutput_mc=.true.
     case("-mc")
        i_arg = i_arg + 1
        lmc=.true.
        loutput_mc=.true.
     case("-density_file","-df")
        i_arg = i_arg + 1
        ldensity_file=.true.
        call get_command_argument(i_arg,s)
        density_file = s
        i_arg = i_arg + 1
     case("-start_step")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("1 or 2 needed for -start_step")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) istep_start
        i_arg= i_arg+1
     case("-mhd_voronoi")
        i_arg = i_arg + 1
        lmhd_voronoi = .true.
        lVoronoi = .true.
        l3D = .true.
        ldensity_file = .false.
        !later lpluto_file = .true.
     case ("-model_1d")
     	i_arg = i_arg + 1
     	lmodel_1d = .true.
     	lmodel_ascii = .false.
      ldensity_file = .false.
     case("-model_ascii")
        i_arg = i_arg + 1
        lmodel_ascii = .true.
        ldensity_file = .false.
     case("-zeeman_polarisation")
     	call error("Zeeman polarisation not yet!")
        i_arg = i_arg + 1
        lzeeman_polarisation=.true.
     case("-accurate_integ")
        i_arg = i_arg + 1
        laccurate_integ = .true.
     case("-art_line_resol")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("resolution (km/s) needed with -art_line_resol !")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) art_hv
        i_arg= i_arg+1
     case("-healpix_lorder")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("l value needed for healpix !")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) healpix_lorder
        i_arg= i_arg+1
        if ( (healpix_lorder < 0).or.(healpix_lorder > 28) ) then
         call error ("healpix l must be positive, below 28!")
        endif
        if ( (healpix_lorder > 7) ) then
        	call warning("healpix l order > 7 resulting in a high number of rays!")
        endif
     case("-Ng_Norder")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("Ng'acc order needed with -Ng_Norder !")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) Ng_Norder
        i_arg= i_arg+1
        if (Ng_Norder <= 1) then
         call warning ("Ng Norder <= 1, turning off Ng'acceleration!")
         lng_acceleration = .false.
        else
         lng_acceleration = .true.
         write(*,*) " Ng's acceleration activated"
        endif
     case("-Ng_Nperiod")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("Ng'acc period needed with -Ng_Nperiod !")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) Ng_Nperiod
        i_arg= i_arg+1
        if (Ng_Nperiod < 0) then
         call warning ("Ng Nperiod < 0 setting to abs!")
         Ng_Nperiod = abs(Ng_Nperiod)
        endif
     case("-Nrays_mc_step")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("Number of rays needed with -Nrays_mc_step !")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) N_rayons_mc
        i_arg= i_arg+1
        if (N_rayons_mc <= 0) then
         call error ("N_rayons_mc must be > 0")
        endif
     case("-max_err")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("relative error needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) dpops_max_error
        i_arg= i_arg+1
        if (dpops_max_error <= 0) then
         call error ("MAx relative error has to be > 0")
        endif
     case("-max_err_sub")
        i_arg = i_arg + 1
        if (i_arg > nbr_arg) call error("relative error sub needed")
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) dpops_sub_max_error
        i_arg= i_arg+1
        if (dpops_sub_max_error <= 0) then
         call error ("Max sub relative error has to be > 0")
        endif
     case("-see_lte")
        i_arg = i_arg + 1
        lforce_lte=.true.
     case("-level_dissolution")
     	call error("Level's dissolution not yet!")
        i_arg = i_arg + 1
        ldissolve =.true.
     case("-phantom")
        i_arg = i_arg + 1
        lphantom_file=.true.
        lVoronoi = .true.
        l3D = .true.
        call get_command_argument(i_arg,s)
        n_phantom_files = 1
        allocate(density_files(n_phantom_files))
        density_files(1) = s
        i_arg = i_arg + 1
        if (.not.llimits_file) limits_file = "phantom.limits"
     case("-phantom-multi","phantom_multi","-phantom-add","-phantom-avg")
        if (s == "-phantom-avg") lphantom_avg = .true.
        i_arg = i_arg + 1
        lphantom_file=.true.
        lphantom_multi = .true. ! not used in practise
        lVoronoi = .true.
        l3D = .true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) n_phantom_files
        allocate(density_files(n_phantom_files))
        do i=1, n_phantom_files
           call get_command_argument(i_arg,s)
           i_arg = i_arg + 1
           density_files(i) = s
        enddo
        if (.not.llimits_file) limits_file = "phantom.limits"
     case("-ascii_SPH")
        i_arg = i_arg + 1
        lascii_SPH_file = .true.
        lVoronoi = .true.
        l3D = .true.
        call get_command_argument(i_arg,s)
        density_file = s
        i_arg = i_arg + 1
        if (.not.llimits_file) limits_file = "phantom.limits"
     case("-SPH_amin")
        lforce_SPH_amin = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) SPH_amin
        i_arg = i_arg + 1
     case("-SPH_amax")
        lforce_SPH_amax = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) SPH_amax
        i_arg = i_arg + 1
     case("-force_Mgas")
        lforce_Mgas = .true.
        i_arg = i_arg + 1
     case("-gadget","-gadget2")
        i_arg = i_arg + 1
        lgadget2_file=.true.
        lVoronoi = .true.
        l3D = .true.
        call get_command_argument(i_arg,s)
        density_file = s
        i_arg = i_arg + 1
        if (.not.llimits_file) limits_file = "gadget2.limits"
     case("-limits_file","-limits")
        i_arg = i_arg + 1
        llimits_file = .true.
        call get_command_argument(i_arg,s)
        limits_file = s
        i_arg = i_arg + 1
     case("-keep_particles")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) SPH_keep_particles
        if ((SPH_keep_particles < 0) .or. (SPH_keep_particles > 1)) then
           call error("keep_particles value must between 0 and 1")
        endif
     case("-sigma_file","-sigma")
        i_arg = i_arg + 1
        lsigma_file=.true.
        call get_command_argument(i_arg,s)
        sigma_file = s
        i_arg = i_arg + 1
     case("-weight_emission")
        i_arg = i_arg+1
        lweight_emission=.true.
     case("-correct_density")
        i_arg = i_arg+1
        lcorrect_density=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) correct_density_factor
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) correct_density_Rin
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) correct_density_Rout
     case("-prodimo")
        i_arg = i_arg + 1
        lprodimo = .true.
        mcfost2ProDiMo_version = 5
        lISM_heating=.true.
        ISR_model = 1
     case("-astrochem","-Astrochem","-AstroChem")
        i_arg = i_arg + 1
        lastrochem = .true.
        ldisk_struct=.true.
        lstop_after_init=.false.
        loutput_UV_field=.true.
     case("-prodimo4")
        i_arg = i_arg + 1
        lprodimo = .true.
        mcfost2ProDiMo_version = 4
     case("-prodimo3")
        i_arg = i_arg + 1
        lprodimo = .true.
        mcfost2ProDiMo_version = 3
     case("-prodimo2")
        i_arg = i_arg + 1
        lprodimo = .true.
        mcfost2ProDiMo_version = 2
     case("-prodimo1")
        i_arg = i_arg + 1
        lprodimo = .true.
        mcfost2ProDiMo_version = 1
     case("-prodimo_input_dir")
        i_arg = i_arg + 1
        lprodimo_input_dir=.true.
        call get_command_argument(i_arg,ProDiMo_input_dir)
        i_arg = i_arg + 1
     case("-prodimo_fPAH")
        i_arg = i_arg + 1
        lforce_ProDiMo_PAH = .true.
        call get_command_argument(i_arg,sProDiMo_fPAH)
        i_arg = i_arg + 1
        read(sProDiMo_fPAH,*) ProDiMo_fPAH
     case("-gap")
        i_arg = i_arg+1
        lgap_Gaussian=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) f_gap_Gaussian
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) r_gap_Gaussian
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) sigma_gap_Gaussian
     case("-only_scatt")
        i_arg = i_arg+1
        lonly_scatt = .true.
     case("-HG","-hg")
        i_arg = i_arg+1
        lHG=.true.
     case("-force_HG","-force_hg")
        i_arg = i_arg+1
        lHG=.true.
        lforce_HG=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) forced_g
     case("-p2m")
        i_arg = i_arg+1
        lProDiMo2mcfost=.true.
     case("-prodimo2mcfost")
        i_arg = i_arg+1
        lProDiMo2mcfost_test=.true.
     case("-spot")
        i_arg = i_arg+1
        lspot=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) T_spot
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) surf_fraction_spot
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) theta_spot
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) phi_spot
     case("-Seb_C")
        i_arg = i_arg+1
        lSeb_Charnoz=.true.
     case("-read_Seb_C")
        i_arg = i_arg+1
        lread_Seb_Charnoz=.true.
     case("-read_Seb_C2")
        i_arg = i_arg + 1
        lread_Seb_Charnoz2=.true.
        call get_command_argument(i_arg,s)
        density_file = s
        i_arg = i_arg + 1
     case("-only_top")
        i_arg = i_arg+1
        lonly_top=.true.
     case("-only_bottom")
        i_arg = i_arg+1
        lonly_bottom=.true.
     case("-correct_Tgas")
        i_arg = i_arg+1
        lcorrect_Tgas=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) correct_Tgas
     case("-force_PAH_equilibrium")
        lforce_PAH_equilibrium=.true.
        i_arg = i_arg+1
     case("-force_PAH_out_equilibrium")
        lforce_PAH_out_equilibrium=.true.
        i_arg = i_arg+1
        if (lforce_PAH_equilibrium) then
           write(*,*) "ERROR: cannot force eq. and out eq."
           write(*,*) "Exiting"
        endif
     case("-grain_size_distrib_file")
        i_arg = i_arg + 1
        lread_grain_size_distrib=.true.
        call get_command_argument(i_arg,s)
        grain_size_file = s
        i_arg = i_arg+1
     case("-Tmax_PAH")
        i_arg = i_arg + 1
        lchange_Tmax_PAH=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) Tmax_PAH
     case("-benchmark_SHG")
        i_arg = i_arg + 1
        lbenchmark_SHG=.true.
     case("-Mathis_field")
        if (.not.lbenchmark_SHG) then
           call error("Mathis field can only be used with the SHG benchmark")
        endif
        i_arg = i_arg + 1
        lMathis_field=.true.
        call get_command_argument(i_arg,s)
        read(s,*) Mathis_field
        i_arg = i_arg+1
     case("-no_T")
        i_arg = i_arg + 1
        lno_T=.true.
     case("-no_SED")
        i_arg = i_arg + 1
        lno_SED=.true.
     case("-ISM_heating")
        i_arg = i_arg + 1
        lISM_heating=.true.
        ISR_model = 1
     case("-chi_ISM")
        i_arg = i_arg + 1
        lISM_heating=.true.
        call get_command_argument(i_arg,s)
        read(s,*) chi_ISM
        i_arg = i_arg + 1
        ISR_model = 1
     case("-ISM_heating_Bate")
        i_arg = i_arg + 1
        lISM_heating=.true.
        ISR_model = 2
     case("-casa")
        i_arg = i_arg + 1
        lcasa=.true.
     case("-cutoff")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) cutoff
        i_arg = i_arg + 1
     case("-nphot_img")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) nphot_img
        i_arg = i_arg + 1
     case("-n_rad")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) n_rad_opt
        i_arg = i_arg + 1
     case("-nz")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) nz_opt
        i_arg = i_arg + 1
     case("-nT")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) n_T_opt
        i_arg = i_arg + 1
     case("-limb_darkening")
        i_arg = i_arg + 1
        llimb_darkening = .true.
        call get_command_argument(i_arg,limb_darkening_file)
        i_arg = i_arg + 1
     case("-max_mem")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) max_mem
        max_mem = max_mem/2. ! facteur a la louche
        i_arg = i_arg + 1
     case("-cavity")
        lcavity = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) cavity%sclht
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) cavity%rref
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) cavity%exp_beta
        i_arg = i_arg + 1
     case("-age")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,system_age)
        i_arg = i_arg + 1
     case("-no_internal_energy")
        i_arg = i_arg + 1
        lno_internal_energy = .true.
     case("-chi_infall")
        i_arg = i_arg + 1
        linfall = .true.
        call get_command_argument(i_arg,s)
        read(s,*) chi_infall
        i_arg = i_arg + 1
     case("-cylindrical_rotation","-cyl_rotation","-cyl_rot")
        i_arg = i_arg + 1
        lcylindrical_rotation = .true.
     case("-phase_function","-phase-function","-phase_function_file","-phase-function-file")
        i_arg = i_arg + 1
        lphase_function_file = .true.
        call get_command_argument(i_arg,s) ; call read_phase_function(s)
        i_arg = i_arg + 1
     case("-tau=1_surface")
        i_arg = i_arg + 1
        ltau1_surface=.true.
     case("-z_scaling_env")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) z_scaling_env
        i_arg = i_arg + 1
     case("-planet_az")
        lplanet_az = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) planet_az
        i_arg = i_arg + 1
     case("-planet")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) which_planet
        write(*,*) "PLANET", which_planet, "AZ=", planet_az
        i_arg = i_arg + 1
     case("-turn-off_planets")
        i_arg = i_arg + 1
        lturn_off_planets = .true.
     case("-turn-off_Lacc")
        i_arg = i_arg + 1
        lturn_off_Lacc = .true.
     case("-correct_density_elongated_cells")
        i_arg = i_arg+1
        lcorrect_density_elongated_cells=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,* ) correct_density_factor_elongated_cells
     case("-pola")
        i_arg = i_arg + 1
        lpola=.true.
     case("-ML","-ml")
        i_arg = i_arg + 1
        lML=.true.
     case("-fix_star","-fix_stars")
        i_arg = i_arg + 1
        lfix_star=.true.
     case("-scale_length_units")
        lscale_length_units = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) scale_length_units_factor
        i_arg = i_arg + 1
     case("-scale_mass_units")
        lscale_mass_units = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) scale_mass_units_factor
        i_arg = i_arg + 1
     case("-ignore_dust")
        i_arg = i_arg + 1
        lignore_dust=.true.
     case("-no_vr")
        i_arg = i_arg + 1
        lupdate_velocities = .true.
        lno_vr = .true.
     case("-no_vz")
        i_arg = i_arg + 1
        lupdate_velocities = .true.
        lno_vz = .true.
     case("-vphi_Kep","-vphi_kep")
        i_arg = i_arg + 1
        lupdate_velocities = .true.
        lvphi_Kep = .true.
     case("-fluffy","-fluffyness")
        lfluffy = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) fluffyness
        i_arg = i_arg + 1
     case("-delete_Hill_sphere")
        i_arg = i_arg + 1
        ldelete_Hill_sphere = .true.
     case("-delete_inside_rsph")
        i_arg = i_arg + 1
        ldelete_inside_rsph = .true.
        call get_command_argument(i_arg,s)
        read(s,*) rsph_min
        i_arg = i_arg + 1
     case("-delete_outside_rsph")
        i_arg = i_arg + 1
        ldelete_outside_rsph = .true.
        call get_command_argument(i_arg,s)
        read(s,*) rsph_max
        i_arg = i_arg + 1
     case("-random_az")
        i_arg = i_arg + 1
        lrandomize_azimuth = .true.
        lrandomize_Voronoi = .true.
     case("-random_gap")
        i_arg = i_arg + 1
        lrandomize_gap = .true.
        lrandomize_Voronoi = .true.
        call get_command_argument(i_arg,s)
        read(s,*) gap_factor
        i_arg = i_arg + 1
     case("-random_outside_gap")
        i_arg = i_arg + 1
        lrandomize_outside_gap = .true.
        lrandomize_Voronoi = .true.
        call get_command_argument(i_arg,s)
        read(s,*) gap_factor
        i_arg = i_arg + 1
     case("-cd","-column_density")
        i_arg = i_arg + 1
        lwrite_column_density = .true.
     case("-mol_cd","-mol_column_density")
        i_arg = i_arg + 1
        lwrite_mol_column_density = .true.
     case("-centre_on_sink")
        i_arg = i_arg + 1
        lcentre_on_sink = .true.
        call get_command_argument(i_arg,s)
        read(s,*) isink_centre
        i_arg = i_arg + 1
     case("-star_bb")
        i_arg = i_arg + 1
        lstar_bb = .true.
     case("-Mdot")
        i_arg = i_arg + 1
        lforce_Mdot = .true.
        call get_command_argument(i_arg,s)
        read(s,*) istar_Mdot
        star_force_Mdot(istar_Mdot) = .true.
        i_arg = i_arg + 1
        call get_command_argument(i_arg,s)
        read(s,*) Mdot
        star_Mdot(istar_Mdot) = Mdot
        i_arg = i_arg + 1
     case("-fargo3d","-fargo")
        i_arg = i_arg + 1
        lfargo3d = .true.
        call get_command_argument(i_arg,fargo3d_dir)
        i_arg = i_arg + 1
        call get_command_argument(i_arg,fargo3d_id)
        i_arg = i_arg + 1
        read(fargo3d_id,*,iostat=ios) i
        if (ios/=0) call error("fargo3d dump number needed")
     case("-athena++","-athena")
        i_arg = i_arg + 1
        lathena = .true.
        call get_command_argument(i_arg,athena_file)
         i_arg = i_arg + 1
      case default
        write(*,*) "Error: unknown option: "//trim(s)
        write(*,*) "Use 'mcfost -h' to get list of available options"
        call exit(0)
     end select
  enddo ! while

  ! Lecture du fichier de parametres
  if (lProDiMo2mcfost) then
     call read_mcfost2ProDiMo(para)
  else
     call read_para(para)
  endif

  if (lfargo3d) then
     l3D = .true.
     if (n_zones > 1) call error("fargo3d mode only work with 1 zone")
     call warning("fargo3d : forcing spherical grid") ! only spherical grid is implemented for now
     disk_zone(1)%geometry = 2
     call read_fargo3d_parameters(fargo3d_dir, fargo3d_id)
  endif
  if (lathena) then
     l3D = .true.
     if (n_zones > 1) call error("athena mode only work with 1 zone")
     call warning("athena : forcing spherical grid") ! only spherical grid is implemented for now
     disk_zone(1)%geometry = 2
     call read_athena_parameters(athena_file)
  endif
  if (lmodel_1d) then
   l3d = .false.
   n_zones = 1
   disk_zone(1)%geometry = 2
   call warning("model_1d : reading 1d  stellar atmosphere model")
   call read_grid_1d(density_file)
  endif
  if (lmodel_ascii .or. lmodel_1d) ldensity_file = .false. !because -df is used to 
                                                           !the pass the model to mcfost!

  if (n_zones > 1) lvariable_dust=.true.

  if (lemission_mol.and.para_version < 2.11) call error("parameter version must be larger than 2.10")
  if (lemission_atom.and.para_version < 2.11) call error("Atomic line RT only available for latest versions")

  if (lno_T) ltemp = .false.
  if (lno_SED) then
     lsed = .false.
     lsed_complete = .false.
  endif

  if (map_size < tiny_real) call error("map size is set to 0")

  if (((.not.limg).and.(.not.ldust_prop)).and.lsepar_pola.and.lscatt_ray_tracing.and.(.not.lscatt_ray_tracing2)) then
     call warning("polarization is turned off in ray-traced SEDs", &
          msg2="it can be turned back on with -rt2")
     lsepar_pola = .false.
  endif

  if (lpola) lsepar_pola = .true.

  if (lsepar_pola) then
     n_Stokes = 4
     if (lsepar_contrib) then
        N_type_flux = 8
     else
        N_type_flux = 4
     endif
  else
     n_Stokes = 1
     if (lsepar_contrib) then
        N_type_flux = 5
     else
        N_type_flux = 1
     endif
  endif

  if (lwrite_column_density .and. .not. ldisk_struct) call error("-cd option requires - or +disk_struct option")

  if (lstar_bb) etoile(:)%lb_body = .true.

  if (lforce_Mdot) then
     do i=1, n_etoiles
        if (i > nmax_stars) call error("Maximal number of stars is hard-coded to 100")
        if (star_force_Mdot(i)) then
           etoile(i)%force_Mdot = .true.
           etoile(i)%Mdot = star_Mdot(i)
        endif
     enddo
  endif

  write(*,*) 'Input file read successfully'

  if ((lmodel_ascii.or.lmodel_1d).and.(lascii_sph_file.or.lphantom_file)) then
   call error("Cannot use Phantom and MHD files at the same time presently.")
  end if

  ! Correction sur les valeurs du .para
  if (lProDiMo) then
     if ((lsed_complete).or.(tab_wavelength(1:7) /= "ProDiMo")) then
        write(*,*) "WARNING: ProDiMo mode, forcing the wavelength grid using ProDiMo_UV3_9.lambda"
        lsed_complete = .false.
        tab_wavelength = prodimo_tab_wavelength
     endif
  endif

  if (lProDiMo2mcfost_test) lProDiMo2mcfost = .true.

  if (ldust_prop) then
     write(*,*) "Computation of dust properties as a function of wavelength"
     ! Make sure the full wavelength range is used to compute
     ! the dust properties

     if (lstop_after_init) then
        ! on change les parametres par default pour gagner du temps
        ! et pour avoir des quantites integrees !!!
        ! BUG ici : +dust_prop renvoie les prop de la 1ere cellule
        !limg=.false.
        !lmono=.false.
        !lmono0=.false.
        lvariable_dust=.false.
        scattering_method=2 ; scattering_method0=2
     else
        if (lvariable_dust) call error("Cannot use +dust_prop when settling is on, use -dust_prop instead")
     endif

     basename_data_dir = "data_dust"
     data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
     call save_data_prop(para,base_para)
  endif

  if (ldisk_struct) then
     write(*,*) "Computation of disk structure"
     basename_data_dir = "data_disk"
     data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
     call save_data_prop(para,base_para)
  endif

  if (lread_Seb_Charnoz) lvariable_dust = .true.

  if (lonly_scatt) l_em_disk_image=.false.
  if (lHG.or.lisotropic) aniso_method=2

  if (nphot_img > tiny_real) nbre_photons_image = max(nphot_img / nbre_photons_loop,1.)
  if (n_rad_opt > 0) n_rad = n_rad_opt
  if (nz_opt > 0) nz = nz_opt
  if (n_T_opt > 0) n_T = n_T_opt

  ! Defining pixel values
  if (lresol) then
     npix_x = nx
     npix_y = ny
  endif

  if (limg) then
     if (l_em_disk_image) then
        write(*,*) "Scattered light + thermal emission map calculation"
     else
        write(*,*) "Scattered light map calculation"
     endif
     ltemp=.false.
     lsed=.false.
     lsed_complete=.false.
     lmono0=.true.
     n_lambda=1
     if (wvl <= 1.0) then ! mcfost fait meme les accords grammaticaux 8-)
        write(*,*) "Calculating image at wavelength =", real(wvl)," micron"
     else
        write(*,*) "Calculating image at wavelength =", real(wvl)," microns"
     endif
     basename_data_dir = "data_"//trim(band)
     if (lreemission_stats) then
        write(*,*) "The [-reemission_stats] option is not relevant for image calculation"
        write(*,*) "It is therefore discarded here"
     endif
  else ! SED : we do not need pixels
     basename_data_dir = "data_th"
     npix_x_save = npix_x ; npix_y_save = npix_y
     npix_x = 1 ; npix_y = 1
     lmono0 = .false.
  endif
  data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
  lmono = lmono0

  if (lProDiMo .and. (mcfost2ProDiMo_version == 1)) then
     ! Version 1: lambda : 13 bins entre 0.0912 et 3410.85 microns donne les 11 bins de ProDiMo
     write(*,*) "***************************"
     write(*,*) "* Modelling for ProDiMo   *"
     write(*,*) "* Forcing wavelength grid *"
     write(*,*) "***************************"
     n_lambda = 39 ; lambda_min = 0.0912 ; lambda_max = 3410.85
     lsed=.true. ; lsed_complete = .true.
  endif

  if (lML) then ! We use the same wavelength in step2 as in the DENT grid
     write(*,*) "***************************************************"
     write(*,*) "* Forcing wavelength grid for xgboost predictions *"
     write(*,*) "***************************************************"
     lsed=.true. ; lsed_complete = .false.
     tab_wavelength = DENT_tab_wavelength
     lscatt_ray_tracing = .false. ; call warning("turning ray-tracing off")
  endif

  ! Discrimination type de run (image vs SED/Temp)
  !                       et
  ! verification coherence du fichier de parametres
  n_lambda2 = n_lambda
  if ((.not.lmono0).and.(lsed).and.(.not.lsed_complete)) call lect_lambda()

  if ((.not.lmono0).and.llimb_darkening) call error("Limb darkening only implementing in imaging mode")

  if ((ltemp.or.lsed.or.lsed_complete).and.(.not.lstop_after_init)) then
     write(*,*) "Thermal equilibrium calculation"
     if (lmono) then
        call error("thermal equilibrium cannot be calculated with only 1 wavelength!", &
             msg2="Set n_lambda to a value higher than 1 in parameter file")
     endif
  endif

  if ((lTemp).and.(.not.(lstop_after_init))) then
     if (lRE_LTE) then
        write(*,*) "Temperature calculation under LTE approximation"
     else
        write(*,*) "Temperature calculation without LTE approximation"
     endif
  endif

  if ((.not.limg).and.(.not.lsed).and.(.not.ltemp).and.(.not.lemission_mol).and.&
       (.not.lemission_atom).and.(.not.lstop_after_init)) then
     write(*,*) "ERROR: Nothing to calculate!"
     write(*,*) "You can: "
     write(*,*) "- use the [-img wavelength] option to calculate an image"
     write(*,*) "- set ltemp, lsed and/or lsed_complete to T"
     write(*,*) "  for thermal equilibrium calculations"
     write(*,*) "- use the [-dust_prop] option to compute global dust properties"
     write(*,*) 'Exiting'
     call exit(1)
  endif

  if (.not.limg) loutput_mc = .true.

  if (lPA) ang_disque = PA
  ang_disque = - ang_disque ! rotation North vers East

  if (lzoom) then
     zoom = opt_zoom
     write(*,*) "Updating zoom =", zoom
  endif

  if (lMC_bins) then
     N_thet = ntheta
     N_phi = nazimuth
  endif

  if (lemission_mol) then
     do imol=1,n_molecules
        call read_molecules_names(imol)
        basename_data_dir2(imol) = "data_"//trim(mol(imol)%name)
     enddo

     if (lmol_LTE) then
        write(*,*) "Molecular line transfer under LTE approximation"
     else
        write(*,*) "NLTE molecular line transfer"
     endif

  endif

  if (lspot) then
     if (lscatt_ray_tracing) call error("stellar spots are not implemented in ray-tracing mode yet")
     if (etoile(1)%fUV > 0.) call error("stellar spots are not implemented if the star has a fUV")
     write(*,*) "Spot will be applied on star #1"
  endif


  if (lpara) then
     if (nb_proc==1) then
        write (*,'(" Parallelized code on ", I3, " processor")')  nb_proc
     else
        write (*,'(" Parallelized code on ", I3, " processors")')  nb_proc
     endif
     write(*,*) " "
  else
     write (*,'(" Sequential code")')
  endif

  if ((l_sym_ima).and.(abs(ang_disque) > 1e-6)) then
     call warning("PA different from zero: removing image symetry")
     l_sym_ima=.false.
  endif

  cmd = 'date'
  call appel_syst(cmd, syst_status)

  data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
  do imol=1, n_molecules
     data_dir2(imol) = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir2(imol))
  enddo

  if (.not.lstop_after_init) call save_data(para,base_para)

  if ((l3D).and.(n_az==1).and.(.not.lVoronoi)) then
     write(*,*) "WARNING: using 3D version of MCFOST with a 2D grid"
  endif
  if (n_az > 1) then
     write(*,*) "Forcing 3D mode"
     l3D = .true.
  endif

  if (linfall) l_sym_ima = .false.

  if (lscatt_ray_tracing) then
     if ((.not. lscatt_ray_tracing1) .and. (.not. lscatt_ray_tracing2)) then
        if (lmono0.and.(.not.l3D)) then
           lscatt_ray_tracing2 = .true.
           write(*,*) "Using ray-tracing method 2"
        else
           lscatt_ray_tracing1 = .true.
           write(*,*) "Using ray-tracing method 1"
        endif
     endif
  else
     lscatt_ray_tracing1 = .false.
     lscatt_ray_tracing2 = .false.
  endif

  ! Forcing azimuthal angle to 0 in rt2
  ! as angles are defined with phi=0 in mind
  if (lscatt_ray_tracing2) then
     if (RT_n_az > 1) then
        write(*,*) "There can be only 1 azimuthal angle in RT method 2"
        RT_n_az = 1
     endif
     RT_az_min = 0
     RT_az_max = 0
  endif

  lonly_LTE = .false.
  lonly_nLTE = .false.
  if (lRE_LTE .and. .not.lRE_nLTE .and. .not. lnRE) lonly_LTE = .true.
  if (lRE_nLTE .and. .not.lRE_LTE .and. .not. lnRE) lonly_nLTE = .true.

  ! Signal handler
  ! do i=1,17
  !    sig=signal(i,sig_handler,-1)
  ! enddo

  return

end subroutine initialisation_mcfost

!********************************************************************

subroutine display_help()

  implicit none

  write(*,*) "usage : mcfost parameter_file [options]"
  write(*,*)
  write(*,*) "mcfost -h, -help : displays this help message"
  write(*,*) "       -v : displays version number, and available updates"
  write(*,*) "       -get_para : downloads the current version of the parameter file"
  write(*,*) "       -get_yorick : downloads the current version of yorick scripts"
  write(*,*) "       -u : updates MCFOST to most recent version"
  write(*,*) "       -update_utils : updates MCFOST_UTILS to most recent version"
  write(*,*) "       -history : displays full MCFOST history since v2.12.9"
  write(*,*) "       -max_mem [GB] : maximum memory that MCFOST can use (approx), default 8"
  write(*,*) " "
  write(*,*) " Main mcfost options"
  write(*,*) "        : -img <wavelength> (microns) : computes image at specified wavelength"
  write(*,*) "        : -mol : calculates molecular emission"
  write(*,*) "        : -atom : calculates atomic lines emission"
  write(*,*) " "
  write(*,*) " Coupling with other codes"
  write(*,*) "        : -prodimo : creates required files for ProDiMo"
  write(*,*) "        : -p2m : reads the results from ProDiMo"
  write(*,*) "        : -astrochem : creates the files for astrochem"
  write(*,*) "        : -phantom <dump> : reads a phantom dump file"
  write(*,*) "        : -gadget : reads a gadget-2 dump file"
  write(*,*) "        : -fargo3d <dir> <id> : reads a fargo3d model"
  write(*,*) "        : -model_ascii_atom <file> : read the <file> from ascii file"
  write(*,*) "        : -mhd_voronoi : interface between grid-based code and Voronoi mesh."
  write(*,*) "        : -model_1d : interface with stellar atmosphere models."
  write(*,*) "        : -athena++ <dump> : reads an athena++ athdf file"
  write(*,*) " "
  write(*,*) " Options related to data file organisation"
  write(*,*) "        : -seed <seed> : modifies seed for random number generator;"
  write(*,*) "                         results stored in 'seed=XXX' directory"
  write(*,*) "        : -root_dir <dir> : results stored in 'root_dir' directory"
  write(*,*) '        : -tmp_dir <dir>  : redults where previous cartial calculations are store (default: ".")'
  write(*,*) "        : -no_backup  : stops if directory data_XX already exists"
  write(*,*) "                        without attempting to backup existing directory"
  write(*,*) "        : -prodimo_input_dir <dir> : input files for ProDiMo"
  write(*,*) " "
  write(*,*) " Options related to images"
  write(*,*) "        : -zoom <zoom> (overrides value in parameter file)"
  write(*,*) "        : -resol <nx> <ny> (overrides value in parameter file)"
  write(*,*) "        : -PA (override value in parameter file)"
  write(*,*) "        : -only_scatt : ignore dust thermal emission"
  write(*,*) "        : -casa : write an image ready for CASA"
  write(*,*) "        : -nphot_img : overwrite the value in the parameter file"
  write(*,*) "        : -rt : use ray-tracing method to compute images or SEDs (on by default)"
  write(*,*) "        : -rt1 or -rt2 : use ray-tracing method and force ray-tracing method"
  write(*,*) "        : -no-rt : do not output the ray-tracing results"
  write(*,*) "        : -mc :  keep Monte-Carlo output in ray-tracing mode"
  write(*,*) "        : -n_MC_bins <n_inclinations> <n_azimuth> (default : 10 1)"
  write(*,*) "        : -planet_az <angle> [deg] : adjust the model azimuth so that the planet is at"
  write(*,*) "                                     desired azimuth in the map"
  write(*,*) "        : -planet <sink_particle_number> : select the sink particle used to"
  write(*,*) "                                           perform the dump rotation"
  write(*,*) "        : -turn-off_planets : sink particles with id > 1 will not emit"
  write(*,*) "        : -turn-off_Lacc : ignore accretion on sink particles"
  write(*,*) "        : -turn-off_dust_subl : ignore dust sublimation"
  write(*,*) " "
  write(*,*) " Options related to temperature equilibrium"
  write(*,*) "        : -no_T : skip temperature calculations, force ltemp to F"
  write(*,*) "        : -Tfile <file> : read a given temperature file"
  write(*,*) "        : -diff_approx : enforce computation of T structure with diff approx."
  write(*,*) "        : -no_diff_approx : compute T structure with only MC method"
  write(*,*) "        : -only_diff_approx : only compute the diffusion approx"
  write(*,*) "        : -tau_dark_zone_obs <tau_dark_zone> (default : 100)"
  write(*,*) "        : -tau_dark_zone_eq_th <tau_dark_zone> (default : 1500)"
  write(*,*) "        : -origin : save origin of packets received the interest bin"
  write(*,*) "        : -rs (remove species) <species_number> <Temperature>"
  write(*,*) "        : -reemission_stats"
  write(*,*) "        : -weight_emission  : weight emission towards disk surface"
  write(*,*) "        : -force_PAH_equilibrium : mainly for testing purposes"
  write(*,*) "        : -force_PAH_out_equilibrium : mainly for testing purposes"
  write(*,*) "        : -Tmax_PAH <T> : changes the maximum temperature allowed for PAH (default: 2500)"
  write(*,*) "        : -ISM_heating : includes heating by ISM radiation"
  write(*,*) "        : -chi_ISM <chi> : changes the chi of ISM radiation (default: 1)"
  write(*,*) "        : -no_internal_energy : ignoring internal energy in Phantom dumps"
  write(*,*) " "
  write(*,*) " Options related to disk structure"
  write(*,*) "        : -disk_struct : computes the density structure and stops:"
  write(*,*) "                         gas_density.fits.gz and dust_density.fits.gz -> density map"
  write(*,*) "                         grid.fits.gz -> radii and height in the grid"
  write(*,*) "                         volume.fits.gz -> volume per cell at each radius"
  write(*,*) "        : -3D : 3D geometrical grid"
  write(*,*) "        : -warp : <h_warp> @ reference radius"
  write(*,*) "        : -tilt : <angle> [degrees]"
  write(*,*) "        : -cavity <h0> <r0> <flaring exponent>"
  write(*,*) "        : -output_J"
  write(*,*) "        : -output_UV_field"
  write(*,*) "        : -puffed_up_rim  <h rim / h0> <r> <delta_r>"
  write(*,*) "        : -density_file or -df <density_file>"
  write(*,*) "        : -sigma_file or -sigma <surface_density_file>"
  write(*,*) "        : -correct_density <factor> <Rmin> <Rmax>"
  write(*,*) "        : -gap <depth> <R> <sigma> [depth is between 0 and 1, R and Sigma in au]"
  write(*,*) "        : -Seb_F <number>  1 = gaussian, 2 = cst diffusion coeff"
  write(*,*) "        : -cutoff <number>, upper limit of the grid [scale height] default = 7"
  write(*,*) "        : -n_rad : overwrite value in parameter file"
  write(*,*) "        : -nz : overwrite value in parameter file"
  write(*,*) "        : -z_scaling_env <scaling_factor> : scale a spherical envelope along the z-axis"
  write(*,*) "        : -correct_density_elongated_cells <factor> : apply a density correction to elongated Voronoi cells"
  write(*,*) "        : -column_density or -cd : generates a fits file with the column densities from each cell"
  write(*,*) " "
  write(*,*) " Options related to star properties"
  write(*,*) "        : -star_bb : forces a black-body for the stellar spectrum"
  write(*,*) "        : -spot <T_spot> <surface_fraction> <theta> <phi>, T_spot in K, theta & phi in degrees"
  write(*,*) "        : -limb_darkening <filename>"
  write(*,*) "        : -Mdot <star_id> <Mdot> : forces accretion rate on star (Msun/yr)"
  write(*,*) " "
  write(*,*) " Options related to dust properties"
  write(*,*) "        : -dust_prop : computes opacity, albedo, asymmetry parameter,"
  write(*,*) "                       polarizability and saves results in data_dust"
  write(*,*) "        : -op <wavelength> (microns) : computes dust properties at"
  write(*,*) "                                    specified wavelength and stops"
  write(*,*) "        : -aggregate <GMM_input_file> <GMM_output_file>"
  write(*,*) "        : -optical_depth_map : generates a map of integrated optical depth"
  write(*,*) "          -od                along radial and vertical directions and stops;"
  write(*,*) "                             results stored in optical_depth_map.fits.gz"
  write(*,*) "        : -average_grain_size : computes average grain size in each cell,"
  write(*,*) "                             weighted by their geometrical cross-section;"
  write(*,*) "                             results stored in average_grain_size.fits.gz"
  write(*,*) "        : -HG : uses an Heynyey-Greenstein function"
  write(*,*) "        : -force_HG <g> : uses an Heynyey-Greenstein function and forces the g value"
  write(*,*) "        : -isotropic : forces isotropic scattering"
  write(*,*) "        : -no_scattering : forces albedo = 0"
  write(*,*) "        : -qsca=qabs : forces albedo = 0.5"
  write(*,*) "        : -phase-function <s11.fits> : uses a tabulated phase function (rt2 only)"
  write(*,*) "        : -tau=1_surface"
  write(*,*) " "
  write(*,*) "        : -Nrays_mc_step <Nray> : Number of rays for angular quadrature in Monte Carlo step"
  write(*,*) " Options related to molecular emission"
  write(*,*) "        : -freeze-out <T>"
  write(*,*) "        : -freeze-out_depletion <relative depletion> between 0 and 1"
  write(*,*) "        : -photo-dissociation"
  write(*,*) "        : -photo-desorption"
  write(*,*) "        : -prodimo"
  write(*,*) "        : -prodimo_fPAH : force a fPAH value for ProDiMo"
  write(*,*) "        : -only_top : molecular emssion from the top half of the disk"
  write(*,*) "        : -only_bottom : molecular emssion from the bottom half of the disk"
  write(*,*) "        : -correct_Tgas <factor> : applies a factor to the gas temperature"
  write(*,*) "        : -chi_infall <value> : v_infall/v_kepler"
  write(*,*) "        : -cylindrical_rotation : forces Keplerian velocity independent of z"
  write(*,*) " "
  write(*,*) " Options related to atomic lines emission"
  !healpix options missing. Waiting finle adaptive scheme.
  write(*,*) "        : -start_step <int> : Select the first step for non-LTE loop (default 1)"
!   write(*,*) "        : -checkpoint <int> : activate checkpointing of non-LTE populations every <int> iterations"
  write(*,*) "        : -solve_ne : force the calculation of electron density"
  write(*,*) "        : -iterate_ne <Nperiod> : Iterate ne with populations every Nperiod"
  write(*,*) "        : -Ndelay_iterate_ne <Ndelay> : Iterate ne with populations after Ndelay"
  write(*,*) "        : -see_lte : Force rate matrix to be at LTE"
!   write(*,*) "        : -level_dissolution : Level's dissolution of hydrogenic ions"
  write(*,*) "        : -accurate_integ : increase the accuracy of the monte carlo angular integration"
  write(*,*) "        : -art_line_resol <v> : resolution of the non-LTE grid of art in km/s"
  write(*,*) "        : -output_rates : write radiative rates, rate matrix and full opacities"
!   write(*,*) "        : -electron_scatt : Lambda-iterate the mean intensity with SEE"
!   write(*,*) "        : -calc_jnu_atom : Stop the code after Jnu_scattering has been computed and written. "
!-> this one could also be use for mol tranfer, like Ng or tab_wavelength
  write(*,*) "        : -max_err <max_err> : max relative error"
  write(*,*) "        : -max_err_sub <max_err_sub> : max relative error for sub-iterations"
  write(*,*) "        : -Ng_Norder <Norder> : Order of Ng's acceleration"
  write(*,*) "        : -Ng_Nperiod <Nperiod> : Cycle of Ng's iteration"
!   write(*,*) "        : -zeeman_polarisation : Stokes profiles Zeeman."
  write(*,*) "        : -safe_stop : stop calculation if time > calc_time_limit"
  write(*,*) "        : -safe_stop_time <real> : calc_time_limit in days "

  write(*,*) " "
  write(*,*) " Options related to phantom"
  write(*,*) "        : -limits_file or limits <limit-file> : x,y,z values used for the Voronoi tesselation"
  write(*,*) "        : -keep_particles <fraction> : fraction of SPH particles to keep for"
  write(*,*) "                                       the Voronoi tesselation (default : 0.99)"
  write(*,*) "        : -ignore_dust : ignore the dust fraction in a phantom dump"
  write(*,*) "        : -age <age> : age used to compute stellar parameters from mass of sink particles"
  write(*,*) "        : -fix_star : do not compute stellar parameters from sink particle, use values in para file"
  write(*,*) "        : -scale_length_units <scaling_factor> : over-ride the length units read in by this factor"
  write(*,*) "        : -scale_mass_units <scaling_factor> : over-ride the mass units read in by this factor"
  write(*,*) "        : -fluffyness <factor> : shift grain sizes between phantom and mcfost"
  write(*,*) "        : -delete_Hill_sphere : delete SPH particles inside Hill spheres of planets"
  write(*,*) "        : -delete_inside_rsph <r> : delete SPH particles inside spherical radius"
  write(*,*) "        : -delete_outside_rsph <r> : delete SPH particles outside spherical radius"
  write(*,*) "        : -no_vr : force the radial velocities to be 0"
  write(*,*) "        : -no_vz : force the vertical velocities to be 0"
  write(*,*) "        : -vphi_Kep : force the azimuthal velocities to be Keplerian"
  write(*,*) "        : -centre_on_sink <number> : centre the model on the sink particle"
  write(*,*) "        : -SPH_amin <size> [mum] : force the grain size that follow the gas"
  write(*,*) "        : -SPH_amax <size> [mum] : force the grain size that follow the dust"
  write(*,*) "                                   (only works with 1 grain size dump)"
  write(*,*) "        : -force_Mgas : force the gas mass to be the value given the mcfost parameter file"
  write(*,*) ""
  write(*,*) "You can find the full documentation at:"
  write(*,*) trim(doc_webpage)

  call exit(0)

end subroutine display_help

!********************************************************************

subroutine save_data_prop(para, base_para)

  character(len=*), intent(in) :: para, base_para

  character(len=1024) :: cmd
  integer ::  syst_status

  if (is_dir(trim(data_dir))) then
     write(*,*) "Directory "//trim(data_dir)//" already exists! Erasing it..."
     cmd = 'rm -Rf '//trim(data_dir)
     call appel_syst(cmd,syst_status)
  endif

  ! Cree le dossier data
  write (*,*) 'Creating directory '//trim(data_dir)
  cmd = 'mkdir -p '//trim(data_dir)//" ; "// &
       ! copie le fichier de parametres
       'cp '//trim(para)//' '//trim(data_dir)//" ; "// &
       ! options de la ligne de commande
       'echo " " >>  '//trim(data_dir)//'/'//trim(base_para)//" ; "// &
       'echo "Executed command line : '//trim(cmd_opt)//'" >> '//trim(data_dir)//'/'//trim(base_para)//" ; "// &
       ! date du calcul
       'date >> '//trim(data_dir)//'/'//trim(base_para)//" ; "// &
       ! machine de calcul
       'uname -a >> '//trim(data_dir)//'/'//trim(base_para)//" ; "// &
       ! id SHA
       'echo sha = '//sha_id//' >> '//trim(data_dir)//'/'//trim(base_para)
  call appel_syst(cmd,syst_status)

  return

end subroutine save_data_prop

!********************************************************************

subroutine save_data(para,base_para)
  !*************************************************
  ! Si le dossier data existe on le sauve
  !*************************************************

  implicit none

  character(len=*), intent(in) :: para, base_para

  integer :: syst_status
  character(len=1024) :: cmd

  logical :: lnew_run, lmove_data
  integer :: etape, etape_start, etape_end
  character(len=512) :: local_data_dir, local_basename_data_dir

  lmove_data = .false.

  if (lonly_diff_approx) return

  if (lmono0) then
     etape_start=1
     etape_end=1
  else
     if (ltemp) then
        etape_start=1
     else
        etape_start=2
     endif
     if (lemission_mol) then
        etape_end=1+n_molecules
     else
        etape_end=1
     endif
  endif

  do etape=etape_start,etape_end

     if (etape == 1) then
        local_data_dir = data_dir
        local_basename_data_dir = basename_data_dir
     else
        local_data_dir = data_dir2(etape-1)
        local_basename_data_dir = basename_data_dir2(etape-1)
     endif


     if (is_dir(trim(root_dir)//"/"//trim(seed_dir))) then
        if (is_dir(trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(local_basename_data_dir))) then
           ! le dossier data existe
           lnew_run = .true.
           lmove_data=.true.
        else ! le dossier data n'existe pas
           lnew_run=.true.
           lmove_data=.false.
        endif ! if y a un dossier data
     else
        lnew_run = .true. ! le dossier data n'existe pas
     endif

     if (lmove_data) then
        if (lno_backup) then
           call error('Directory '//trim(local_data_dir)//' already exists : exiting!')
        else
           write (*,*) 'Directory '//trim(local_data_dir)//' already exists : backing it up in '//trim(local_data_dir)//'_old'
           if (is_dir(trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(local_basename_data_dir)//'_old')) then
              call error('Directory '//trim(local_data_dir)//'_old already exists : exiting!')
           endif
           cmd = 'mv '//trim(local_data_dir)//' '//trim(local_data_dir)//'_old'
           call appel_syst(cmd,syst_status)
        endif
     endif

     if (lnew_run) then
        ! Cree le dossier data
        write (*,*) 'Creating directory '//trim(local_data_dir)
        cmd = 'mkdir -p '//trim(local_data_dir)//" ; "// &
             ! copie le fichier de parametres
             'cp '//trim(para)//' '//trim(local_data_dir)//" ; "// &
             ! options de la ligne de commande
             'echo " " >>  '//trim(local_data_dir)//'/'//trim(base_para)//" ; "// &
             'echo "Executed command line : '//trim(cmd_opt)//'" >> '//trim(local_data_dir)//'/'//trim(base_para)//" ; "// &
             ! date du calcul
             'date >> '//trim(local_data_dir)//'/'//trim(base_para)//" ; "// &
             ! machine de calcul
             'uname -a >> '//trim(local_data_dir)//'/'//trim(base_para)//" ; "// &
             ! id SHA
             'echo sha = '//sha_id//' >> '//trim(local_data_dir)//'/'//trim(base_para)
        ! Copie du fichier lambda si besoin
        if (lsed.and.(.not.lsed_complete).and.(.not.lmono0).and.(etape==1)) then
           cmd = trim(cmd)//' ; cp '//trim(lambda_filename)//' '//trim(local_data_dir)
        endif
        call appel_syst(cmd,syst_status)
     endif
  enddo ! etape

  return

end subroutine save_data

!********************************************************************

end module init_mcfost
