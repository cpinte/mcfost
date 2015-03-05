module init_mcfost

  use parametres
  use disk
  use naleat
  use grains, only : aggregate_file, mueller_aggregate_file
  use em_th, only : specie_removed, T_rm, Tfile
  use wall
  use molecular_emission
  use ray_tracing
  !$ use omp_lib
  use benchmarks
  use read_params
  use input
  use mem
  use ProdiMo
  use utils
  use filter

  implicit none

  contains

subroutine initialisation_mcfost()

  implicit none

  integer :: ios, nbr_arg, i_arg, nx, ny, syst_status, imol, mcfost_no_disclaimer, n_dir, i
  integer :: current_date, update_date, mcfost_auto_update
  real(kind=db) :: wvl
  real :: opt_zoom, utils_version, PA

  character(len=512) :: cmd, s, str_seed
  character(len=4) :: n_chiffres
  character(len=128)  :: fmt1

  logical :: lresol, lPA, lzoom, lmc, ln_zone, lHG, lonly_scatt, lupdate, lno_T

  write(*,*) "You are running MCFOST "//trim(mcfost_release)
  write(*,*) "Git SHA = ", sha_id

  lmc = .false.

  call get_environment_variable('HOME',home)
  if (home == "") then
     home="./"
  else
     home=trim(home)//"/"
  endif

  ! Pour code sequentiel
  nb_proc=1 ; lpara=.false.

  ! Pour code parallel
  !$omp parallel default(none) &
  !$omp shared(nb_proc, lpara)
  !$ nb_proc=omp_get_num_threads() ; lpara=.true.
  !$omp end parallel

  ldisk_struct = .false.
  ldust_prop = .false.
  ldust_gas_ratio = .false.
  lstop_after_init= .false.
  lwall=.false.
  lgap=.false.
  lpah=.false.
  ln_zone=.false. ; n_zones=1
  limg=.false.
  lorigine=.false.
  laggregate=.false.
  l3D=.false.
  lresol=.false.
  lPA = .false.
  lzoom=.false.
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
  lpuffed_rim = .false.
  lopacity_wall = .false.
  lno_backup = .false.
  loutput_J = .false.
  loutput_UV_field = .false.
  laverage_grain_size = .false.
  llinear_grid=.false.
  lr_subdivide=.false.
  lfreeze_out = .false.
  l_em_disk_image = .true.
  lisotropic = .false.
  lno_scattering = .false.
  lqsca_equal_qabs = .false.
  lscatt_ray_tracing=.false.
  lscatt_ray_tracing1=.false.
  lscatt_ray_tracing2=.false.
  loutput_mc=.true.
  ldensity_file=.false.
  lsigma_file = .false.
  ldebris=.false.
  lkappa_abs_grain=.false.
  lweight_emission=.false.
  lprodimo=.false.
  lProDiMo_input_dir=.false.
  lProDiMo2mcfost=.false.
  lProDiMo2mcfost_test = .false.
  lforce_ProDiMo_PAH = .false.
  lforce_diff_approx = .false.
  lstrat_SPH = .false. ; lno_strat_SPH=.true.
  lstrat_SPH_bin = .false. ; lno_strat_SPH_bin=.true.
  lgap_ELT=.false.
  lLaure_SED=.false.
  lforce_T_Laure_SED = .false.
  lspot = .false.
  lSeb_Charnoz = .false.
  lread_Seb_Charnoz = .false.
  lread_Seb_Charnoz2 = .false.
  lforce_1st_scatt = .false.
  lold_grid = .false.
  lonly_bottom = .false.
  lonly_top = .false.
  lcorrect_Tgas = .false.
  lcorrect_density=.false.
  lremove = .false.
  lonly_scatt = .false.
  lHG = .false.
  lforce_PAH_equilibrium=.false.
  lforce_PAH_out_equilibrium=.false.
  lread_grain_size_distrib=.false.
  lMathis_field = .false.
  lchange_Tmax_PAH=.false.
  lno_T = .false.
  lISM_heating = .false.

  ! Geometrie Grille
  lcylindrical=.true.
  lspherical=.not.lcylindrical

  ! Methodes par defaut
  RT_sed_method = 1


  ! Test if MCFOST_UTILS is defined
  call get_environment_variable('MCFOST_UTILS',mcfost_utils)
  if (mcfost_utils == "") then
     write(*,*) "ERROR: environnement variable MCFOST_UTILS is not defined."
     write(*,*) "Exiting."
     stop
  endif
  call get_environment_variable('MY_MCFOST_UTILS',my_mcfost_utils)

  ! Ligne de commande
  call get_command(cmd_opt)

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
           if (syst_status /=0) then
              write(*,*) "ERROR : MCFOST did not manage to restart itself"
              write(*,*) "Exiting"
           endif
           stop
        endif ! lupdate
     endif
  endif

  ! Directories to search (ordered)
  if (my_mcfost_utils == "") then
     !write(*,*) "WARNING: environnement variable MY_MCFOST_UTILS is not defined."
     allocate(search_dir(2)) ; n_dir = 2
     search_dir(1) = "." ; search_dir(2) = mcfost_utils ;
  else
     allocate(search_dir(3)) ; n_dir = 3
     search_dir(1) = "." ; search_dir(2) = my_mcfost_utils ; search_dir(3) = mcfost_utils ;
  endif

  allocate(dust_dir(n_dir),mol_dir(n_dir),star_dir(n_dir),lambda_dir(n_dir))
  dust_dir(1) = "./" ; mol_dir(1) = "./" ; star_dir(1) = "./" ; lambda_dir(1) = "./" ;
  do i=2,n_dir
     dust_dir(i)   = trim(search_dir(i))//"/Dust/"
     mol_dir(i)    = trim(search_dir(i))//"/Molecules/"
     star_dir(i)   = trim(search_dir(i))//"/Stellar_Spectra/"
     lambda_dir(i) = trim(search_dir(i))//"/Lambda/"
  enddo

  ! Nbre d'arguments
  nbr_arg = command_argument_count()
  if (nbr_arg < 1) call display_help()

  call get_command_argument(1,para)

  ! Basic options
  if (para(1:1)=="-") then
     if (para(2:2)=="v") then ! mcfost version
        call mcfost_v()
     else if (para(2:2)=="h") then ! mcfost history
        call mcfost_history()
     else if (para(2:6)=="setup") then ! download the utils and para file the 1st time the code is used
        call mcfost_setup()
     else if (para(2:9)=="get_para") then ! download current reference file
        call mcfost_get_ref_para()
     else if ((para(2:11)=="get_manual") .or. (para(2:8)=="get_doc")) then ! download current manual
        call mcfost_get_manual()
     else if (para(2:14)=="get_yorick") then ! force update utils
        call get_yorick()
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
     stop
  endif

  utils_version =  get_mcfost_utils_version()
  if (utils_version /= required_utils_version) then
     write(*,*) "ERROR: wrong version of the MCFOST_UTILS database"
     write(*,*) "Utils:", utils_version, "required:",required_utils_version
     write(*,*) "Please update with mcfost -update_utils"
     write(*,*) "Exiting."
     stop
  endif


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
        stop
     case("-seed")
        lseed=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : seed needed"
           stop
        endif
        call get_command_argument(i_arg,str_seed)
        read(str_seed,*,iostat=ios) seed
        if (ios/=0) then
           write(*,*) "Error : seed needed"
           stop
        endif
        write(*,*) "Updating seed =", seed
        seed_dir = "seed="//trim(str_seed)
        i_arg = i_arg+1
     case("-img")
        limg=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : wavelength needed"
           stop
        endif
        call get_command_argument(i_arg,band)
        read(band,*,iostat=ios) wvl
        if (ios/=0) then
           write(*,*) "Error : wavelength needed"
           stop
        endif
        i_arg = i_arg+1
     case("-op")
        limg=.true.
        lopacite_only=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : wavelength needed"
           stop
        endif
        call get_command_argument(i_arg,band)
        read(band,*,iostat=ios) wvl
        if (ios/=0) then
           write(*,*) "Error : wavelength needed"
           stop
        endif
        i_arg = i_arg+1
     case("-gap")
        lgap=.true.
        write(*,*) "Gap"
        i_arg = i_arg+1
     case("-n-zone")
        ln_zone=.true.
        write(*,*) "n-zone disk calculation"
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : Number of zones needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) n_zones
        if (ios/=0) then
           write(*,*) "Error : Number of zones needed"
           stop
        endif
        i_arg = i_arg+1
     case("-origin")
        lorigine=.true.
        i_arg = i_arg+1
     case("-zoom")
        i_arg = i_arg+1
        lzoom = .true.
        if (i_arg > nbr_arg) then
           write(*,*) "Error : zoom needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) opt_zoom
        if (ios/=0) then
           write(*,*) "Error : zoom needed"
           stop
        endif
        i_arg = i_arg+1
     case("-pah")
        lpah=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : emmisivity needed"
           stop
        endif
        call get_command_argument(i_arg,model_pah)
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : grain type needed"
           stop
        endif
        call get_command_argument(i_arg,pah_grain)
        i_arg = i_arg+1
     case("-aggregate")
        laggregate=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : GMM input file needed"
           stop
        endif
        call get_command_argument(i_arg,aggregate_file)
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : GMM output file needed"
           stop
        endif
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
        if (i_arg > nbr_arg) then
           write(*,*) "Error : h_warp needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) z_warp
        i_arg= i_arg+1
     case("-rs")
        lremove=.true.
        i_arg= i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : specie number needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) specie_removed
        i_arg= i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : Temperature needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) T_rm
        i_arg= i_arg+1
     case("-strat_SPH")
        lstrat_SPH=.true.
        lno_strat_SPH=.false.
        i_arg = i_arg+1
     case("-no_strat_SPH")
        lstrat_SPH=.true.
        lno_strat_SPH=.true.
        i_arg = i_arg+1
     case("-strat_SPH_bin")
        lstrat_SPH_bin=.true.
        lno_strat_SPH_bin=.false.
        i_arg = i_arg+1
     case("-no_strat_SPH_bin")
        lstrat_SPH_bin=.true.
        lno_strat_SPH_bin=.true.
        i_arg = i_arg+1
     case("-resol")
        lresol=.true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : resolution needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) nx
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : resolution needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*,iostat=ios) ny
        i_arg= i_arg+1
     case("-PA")
        lPA = .true.
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : PA needed"
           stop
        endif
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
     case("-output_J")
        loutput_J=.true.
        i_arg = i_arg+1
     case("-output_UV_field")
        loutput_UV_field=.true.
        i_arg = i_arg+1
     case("-average_grain_size")
        laverage_grain_size=.true.
        i_arg = i_arg+1
     case("-killing_level")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : killing_level needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*) n_dif_max_eq_th
        i_arg = i_arg+1
     case("-tau_dark_zone_eq_th")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : tau_dark_zone needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*) tau_dark_zone_eq_th
        i_arg = i_arg+1
     case("-tau_dark_zone_obs")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : tau_dark_zone needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        read(s,*) tau_dark_zone_obs
        i_arg = i_arg+1
     case("-root_dir")
        i_arg = i_arg+1
        if (i_arg > nbr_arg) then
           write(*,*) "Error : root_dir needed"
           stop
        endif
        call get_command_argument(i_arg,s)
        root_dir=trim(root_dir)//"/"//s
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
     case("-puffed_up_rim")
        lpuffed_rim = .true.
        if (i_arg + 3 > nbr_arg) then
           write(*,*) "Error : rim parameters needed"
           stop
        endif
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
     case("-opacity_wall")
        lopacity_wall = .true.
        if (i_arg + 2 > nbr_arg) then
           write(*,*) "Error : wall parameters needed"
           stop
        endif
        i_arg = i_arg+1
        call get_command_argument(i_arg,s)
        read(s,*) h_wall
        i_arg = i_arg+1
        call get_command_argument(i_arg,s)
        read(s,*) tau_wall
        i_arg = i_arg+1
     case("-spherical")
        lcylindrical=.false.
        lspherical=.true.
        i_arg = i_arg + 1
     case("-no_backup")
        lno_backup=.true.
        i_arg = i_arg + 1
     case("-Tfile")
        i_arg = i_arg + 1
        call get_command_argument(i_arg,Tfile)
        i_arg = i_arg + 1
     case("-linear_grid")
        i_arg = i_arg + 1
        llinear_grid=.true.
     case("-r_subdivide")
        i_arg = i_arg + 1
        lr_subdivide=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) r_subdivide
     case("-freeze_out","-freeze-out")
        i_arg = i_arg + 1
        lfreeze_out=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) T_freeze_out
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
     case("-sigma_file","-sigma")
        i_arg = i_arg + 1
        lsigma_file=.true.
        call get_command_argument(i_arg,s)
        sigma_file = s
        i_arg = i_arg + 1
     case("-debris")
        i_arg = i_arg+1
        ldebris=.true.
        if (i_arg > nbr_arg) then
           write(*,*) "Error : debris disk structure file needed"
           stop
        endif
        call get_command_argument(i_arg,debris_file)
        i_arg = i_arg+1
     case("-kappa_abs_grain")
        i_arg = i_arg+1
        lkappa_abs_grain = .true.
     case("-dg_ratio")
        i_arg = i_arg+1
        ldust_gas_ratio = .true.
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
     case("-gap_ELT")
        i_arg = i_arg+1
        lgap_ELT=.true.
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) r_gap_ELT
        call get_command_argument(i_arg,s)
        i_arg = i_arg + 1
        read(s,*) sigma_gap_ELT
     case("-only_scatt")
        i_arg = i_arg+1
        lonly_scatt = .true.
     case("-HG","-hg")
        i_arg = i_arg+1
        lHG=.true.
     case("-p2m")
        i_arg = i_arg+1
        lProDiMo2mcfost=.true.
     case("-prodimo2mcfost")
        i_arg = i_arg+1
        lProDiMo2mcfost_test=.true.
     case("-Laure_SED")
        i_arg = i_arg+1
        lLaure_SED = .true.
         call get_command_argument(i_arg,Laure_SED_filename)
        i_arg = i_arg + 1
     case("-Laure_SED_force_T")
        i_arg = i_arg+1
        lLaure_SED = .true.
         call get_command_argument(i_arg,Laure_SED_filename)
        i_arg = i_arg + 1
        lforce_T_Laure_SED = .true.
     case("-spot")
        i_arg = i_arg+1
        lspot=.true.
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
     case("-force_1st_scatt")
        i_arg = i_arg+1
        lforce_1st_scatt=.true.
        write(*,*) "WARNING: forcing 1st scattering event when tau < 10"
     case("-old_grid")
        i_arg = i_arg+1
        lold_grid=.true.
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
           write(*,*) "ERROR: Mathis field can only be used with the SHG benchmark"
           write(*,*) "Exiting"
           stop
        endif
        i_arg = i_arg + 1
        lMathis_field=.true.
        call get_command_argument(i_arg,s)
        read(s,*) Mathis_field
        i_arg = i_arg+1
     case("-no_T")
        i_arg = i_arg + 1
        lno_T=.true.
     case("-ISM_heating")
        i_arg = i_arg + 1
        lISM_heating=.true.
     case("-casa")
        i_arg = i_arg + 1
        lcasa=.true.
     case default
        call display_help()
     end select
  enddo ! while

  ! Display the disclaimer if needed
  mcfost_no_disclaimer = 0
  call get_environment_variable('MCFOST_NO_DISCLAIMER',s)
  if (s/="") read(s,*) mcfost_no_disclaimer
  if (mcfost_no_disclaimer == 0) call display_disclaimer()

  ! Lecture du fichier de parametres
  if (lProDiMo2mcfost) then
     call read_mcfost2ProDiMo()
  else
     call read_para()
  endif
  call check_init()

  if (lemission_mol.and.para_version < 2.11) then
     write(*,*) "ERROR: parameter version must be larger than 2.10"
     write(*,*) "line transfer"
     write(*,*) "Exiting."
     stop
  endif

  if (map_size < tiny_real) then
     write(*,*) "ERROR: map size is set to 0"
     write(*,*) "Exiting"
     stop
  endif

  if ((.not.limg).and.lsepar_pola.and.lscatt_ray_tracing.and.(.not.lscatt_ray_tracing2)) then
     write(*,*) "WARNING: polarization is turned off in ray-traced SEDs"
     write(*,*) "         it can be turned back on with -rt2"
     lsepar_pola = .false.
  endif

  if (lno_T) ltemp = .false.

  write(*,*) 'Input file read successfully'


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
        lstrat=.false.
        scattering_method=2
     endif
     basename_data_dir = "data_dust"
     data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
     call save_data()
  endif

  if (ldisk_struct) then
     write(*,*) "Computation of disk structure"
     basename_data_dir = "data_disk"
     data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
     call save_data()
  endif

  if (ln_zone) then
     !lstrat=.true.
     write(*,*) "WARNING: lstrat is not set automatically to TRUE !!!!!!!!!"

     if (exp_strat > 1.0e-6) then
        write(*,*) "!!!  WARNING  !!!"
        write(*,*) "Using the [-2zone] option automatically turns on the disk settling "
        write(*,*) "To avoid this, you should set [exp_strat] to 0.0 in the parameter file"
        write(*,*) "Continuing anyway"
     endif
  endif

  if (lread_Seb_Charnoz) lstrat = .true.


  if (lonly_scatt) l_em_disk_image=.false.
  if (lHG.or.lisotropic) aniso_method=2


  ! Discrimination type de run (image vs SED/Temp)
  !                       et
  ! verification cohérence du fichier de paramètres
  n_lambda2 = n_lambda
  if ((.not.lmono0).and.(lsed).and.(.not.lsed_complete)) call lect_lambda()

  if (limg) then
     if (l_em_disk_image) then
        write(*,*) "Scattered light + thermal emission map calculation"
     else
        write(*,*) "Scattered light map calculation"
     endif
     ltemp=.false.
     lsed=.false.
     lsed_complete=.false.
     lmono=.true.
     lmono0=.true.
     n_lambda=1
     if (wvl <= 1.0) then ! mcfost fait meme les accords grammaticaux 8-)
        write(*,*) "Calculating image at wavelength =", wvl," micron"
     else
        write(*,*) "Calculating image at wavelength =", wvl," microns"
     endif
     basename_data_dir = "data_"//trim(band)
     if (lreemission_stats) then
        write(*,*) "The [-reemission_stats] option is not relevant for image calculation"
        write(*,*) "It is therefore discarded here"
     endif
  endif

  if ((ltemp.or.lsed.or.lsed_complete).and.(.not.lstop_after_init)) then
     write(*,*) "Thermal equilibrium calculation"
     if (lforce_1st_scatt) then
        write(*,*) "The [-force_1st_scatt] option is not relevant for SED calculation"
        write(*,*) "It is therefore discarded here"
        lforce_1st_scatt = .false.
     endif

     if (lmono) then
        write(*,*) "Error : thermal equilibrium cannot be calculated with only 1 wavelength!"
        write(*,*) "      Set n_lambda to a value higher than 1 in parameter file"
        write(*,*) 'Exiting'
        stop
     endif
     if ((.not.lstop_after_init)) basename_data_dir = "data_th"
  endif

  if ((lTemp).and.(.not.(lstop_after_init))) then
     if (lRE_LTE) then
        write(*,*) "Temperature calculation under LTE approximation"
     else
        write(*,*) "Temperature calculation without LTE approximation"
     endif
  endif

  if ((.not.limg).and.(.not.lsed).and.(.not.ltemp).and.(.not.lemission_mol).and.(.not.lstop_after_init)) then
     write(*,*) "Error : Nothing to calculate!"
     write(*,*) "You can: "
     write(*,*) "- use the [-img wavelength] option to calculate an image"
     write(*,*) "- set ltemp, lsed and/or lsed_complete to T"
     write(*,*) "  for thermal equilibrium calculations"
     write(*,*) "- use the [-dust_prop] option to compute global dust properties"
     write(*,*) 'Exiting'
     stop
  endif

  if ( (lforce_1st_scatt).and.(lscatt_ray_tracing) ) then
     write(*,*) "ERROR: force_1st_scatt is not compatible with rt"
     write(*,*) "Exiting"
     stop
  endif

  if (lresol) then
     igridx = nx
     igridy = ny

      maxigrid = max(igridx, igridy)

      if (igridx == igridy) then
         deltapix_x = 1
         deltapix_y = 1
      else if (igridx > igridy) then
         deltapix_x = 1
         deltapix_y = 1 - (igridx/2) + (igridy/2)
      else
         deltapix_x = 1 - (igridy/2) + (igridx/2)
         deltapix_y = 1
      endif
      size_pix=maxigrid/(map_size)
  endif

  if (lPA) ang_disque = PA

  if (lzoom) then
     zoom = opt_zoom
     write(*,*) "Updating zoom =", zoom
  endif

! BUG dans version methode 2 de dust_map (ray-tracing)
!  if (lsed) then
!     igridx = 1
!     igridy = 1
!  endif

  if (lstrat_SPH) lstrat=.true.

  if (lemission_mol) then
     lscatt_ray_tracing = .false. ! tmp : scatt ray-tracing has no sense yet for mol emssion
     do imol=1,n_molecules
        call read_molecules_names(imol)
        basename_data_dir2(imol) = "data_"//trim(mol(imol)%name)
     enddo

     if (para_version < 2.07) then
        nTrans = nTrans_tot
        mol(1)%nTrans_raytracing = nTrans_tot
     endif

     if (lmol_LTE) then
        write(*,*) "Molecular line transfer under LTE approximation"
     else
        write(*,*) "NLTE molecular line transfer"
     endif

  endif

  if (lProDiMo .and. (mcfost2ProDiMo_version == 1) ) then
     ! Version 1: lambda : 13 bins entre 0.0912 et 3410.85 microns donne les 11 bins de ProDiMo
     write(*,*) "***************************"
     write(*,*) "* Modelling for ProDiMo   *"
     write(*,*) "* Forcing wavelength grid *"
     write(*,*) "***************************"
     n_lambda = 39
     lambda_min = 0.0912
     lambda_max = 3410.85
     lsed_complete = .true.
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
     write(*,*) "WARNING: PA different from zero: removing image symetry"
     l_sym_ima=.false.
  endif

  cmd = 'date'
  call appel_syst(cmd, syst_status)


  data_dir = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir)
  do imol=1, n_molecules
     data_dir2(imol) = trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(basename_data_dir2(imol))
  enddo

  call save_data()

  if ((l3D).and.(n_az==1)) then
     write(*,*) "WARNING: using 3D version of MCFOST with a 2D grid"
  endif


  if (lscatt_ray_tracing .and. (.not. lscatt_ray_tracing1) .and. (.not. lscatt_ray_tracing2)) then
     if (lmono0.and.(.not.l3D)) then
        lscatt_ray_tracing2 = .true.
        write(*,*) "Using ray-tracing method 2"
     else
        lscatt_ray_tracing1 = .true.
        write(*,*) "Using ray-tracing method 1"
     endif
  endif

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
  write(*,*) "mcfost -help : displays this help message"
  write(*,*) "       -v : displays version number, and available updates"
  write(*,*) "       -get_para or get_doc : downloads the current version of the parameter file"
  write(*,*) "       -get_manual : downloads the current version of manual"
  write(*,*) "       -get_yorick : downloads the current version of yorick scripts"
  write(*,*) "       -u : updates MCFOST to most recent version"
  write(*,*) "       -update_utils : updates MCFOST_UTILS to most recent version"
  write(*,*) "       -h : displays full MCFOST history since v2.12.9"
  write(*,*) " "
  write(*,*) " Main mcfost options"
  write(*,*) "        : -img <wavelength> (microns) : computes image at specified wavelength"
  write(*,*) "        : -rt : use ray-tracing method to compute images or SEDs"
  write(*,*) "        : -mol : calculates molecular emission"
  write(*,*) "        : -prodimo : creates the files for ProDiMo"
  write(*,*) " "
  write(*,*) " Options related to data file organisation"
  write(*,*) "        : -seed <seed> : modifies seed for random number generator;"
  write(*,*) "                         results stored in 'seed=XXX' directory"
  write(*,*) "        : -root_dir <root_dir> : results stored in 'root_dir' directory"
  write(*,*) "        : -no_backup  : stops if directory data_XX already exists"
  write(*,*) "                        without attempting to backup existing directory"
  write(*,*) "        : -prodimo_input_dir <dir> : input files for ProDiMo"
  write(*,*) " "
  write(*,*) " Options related to images"
  write(*,*) "        : -zoom <zoom> (overrides value in parameter file)"
  write(*,*) "        : -resol <nx> <ny> (overrides value in parameter file)"
  write(*,*) "        : -PA (override value in parameter file)"
  write(*,*) "        : -only_scatt : ignore dust thermal emission"
  write(*,*) "        : -force_1st_scatt : uses forced scattering in image calculation;"
  write(*,*) "                             useful for optically thin disk in MC mode"
!  write(*,*) "        : -rt1 : use ray-tracing method 1 (SED calculation)"
!  write(*,*) "        : -rt2 : use ray-tracing method 2 (image calculation)"
  write(*,*) "        : -mc  : keep Monte-Carlo output in ray-tracing mode"
  write(*,*) "        : -casa : write an image ready for CASA"
  write(*,*) " "
  write(*,*) " Options related to temperature equilibrium"
  write(*,*) "        : -no_T : skip temperature calculations, force ltemp to F"
  write(*,*) "        : -diff_approx : enforce computation of T structure with diff approx."
  write(*,*) "        : -no_diff_approx : compute T structure with only MC method"
  write(*,*) "        : -only_diff_approx : only compute the diffusion approx"
  write(*,*) "        : -tau_dark_zone_obs <tau_dark_zone> (default : 100)"
  write(*,*) "        : -tau_dark_zone_eq_th <tau_dark_zone> (default : 1500)"
  write(*,*) "        : -origin : save origin of packets received the interest bin"
  write(*,*) "        : -rs (remove specie) <specie_number> <Temperature>"
  write(*,*) "        : -reemission_stats"
  write(*,*) "        : -weight_emission  : weight emission towards disk surface"
  write(*,*) "        : -force_PAH_equilibrium : mainly for testing purposes"
  write(*,*) "        : -force_PAH_out_equilibrium : mainly for testing purposes"
  write(*,*) "        : -Tmax_PAH <T> : changes the maximum temperature allowed for PAH (default: 2500)"
  write(*,*) "        : -ISM_heating : includes heating by ISM radiation"
  write(*,*) " "
  write(*,*) " Options related to disk structure"
  write(*,*) "        : -disk_struct : computes the density structure and stops:"
  write(*,*) "                         gas_density.fits.gz and dust_density.fits.gz -> density map"
  write(*,*) "                         grid.fits.gz -> radii and height in the grid"
  write(*,*) "                         volume.fits.gz -> volume per cell at each radius"
  write(*,*) "        : -r_subdivide <radius>  : forces cell subdivision elsewhere"
  write(*,*) "                                   than at inner radius"
  write(*,*) "        : -3D : 3D geometrical grid"
  write(*,*) "        : -warp : <h_warp> @ reference radius"
  write(*,*) "        : -strat_SPH"
  write(*,*) "        : -no_strat_SPH"
  write(*,*) "        : -output_J"
  write(*,*) "        : -output_UV_field"
  write(*,*) "        : -puffed_up_rim  <h rim / h0> <r> <delta_r>"
!  write(*,*) "        : -wall <h_wall> <tau_wall>, implies 3D, density wall"
  write(*,*) "        : -opacity_wall <h_wall> <tau_wall>, ONLY an opacity wall in MC,"
  write(*,*) "                            NOT a density wall"
  write(*,*) "        : -linear_grid : linearly spaced grid"
  write(*,*) "        : -density_file or -df <density_file>"
  write(*,*) "        : -sigma_file or -sigma <surface_density_file>"
  write(*,*) "        : -debris <debris_disk_structure_file>"
  write(*,*) "        : -correct_density <factor> <Rmin> <Rmax>"
  write(*,*) "        : -gap_ELT <R> <sigma>"
  write(*,*) "        : -Laure_SED <file>"
  write(*,*) "        : -Laure_SED_force_T <file>"
  write(*,*) "        : -Seb_F <number>  1 = gaussian, 2 = cst diffusion coeff"
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
  write(*,*) "        : -isotropic : forces isotropic scattering"
  write(*,*) "        : -no_scattering : forces albedo = 0"
  write(*,*) "        : -qsca=qabs : forces albedo = 0.5"
  write(*,*) " "
  write(*,*) " Options related to molecular emission"
  write(*,*) "        : -freeze_out <T>"
  write(*,*) "        : -prodimo"
  write(*,*) "        : -prodimo_fPAH : force a fPAH value for ProDiMo"
  write(*,*) "        : -only_top : molecular emssion from the top half of the disk"
  write(*,*) "        : -only_bottom : molecular emssion from the bottom half of the disk"
  write(*,*) "        : -correct_Tgas <factor> : applies a factor to the gas temperature"

  stop

end subroutine display_help

!********************************************************************

subroutine display_disclaimer()

  character(len=10) :: accept
  character(len=512) :: cmd
  integer :: syst_status

  write(*,*) "*******************************************"
  write(*,*) "*          MCFOST DISCLAIMER              *"
  !write(*,*) "*     You are running MCFOST "//trim(mcfost_release)//"       *"
  write(*,*) "*    @ C. Pinte, F. Menard, G. Duchene    *"
  write(*,*) "*                                         *"
  write(*,*) "* MCFOST is available on a collaborative  *"
  write(*,*) "* basis. Using MCFOST implies that you    *"
  write(*,*) "* agree to :                              *"
  write(*,*) "*  - offer us co-author right on any      *"
  write(*,*) "* resulting publication.                  *"
  write(*,*) "*  - NOT distribute MCFOST without our    *"
  write(*,*) "* explicit agreement.                     *"
  write(*,*) "*  - contact us if you initiate a new     *"
  write(*,*) "* scientific project with MCFOST.         *"

  if (.not.is_file(trim(home)//"/.mcfost/accept_disclaimer_"//mcfost_release)) then
     write(*,*) "*                                         *"
     write(*,*) "* Do you accept ? (yes/no)"
     read(*,*) accept

     if ( (accept(1:3) == "yes").or.(accept(1:3) == "Yes").or.(accept(1:3) == "YES") &
          .or.(accept(1:1) == "Y").or.(accept(1:1) == "y") ) then
        cmd = 'mkdir -p ~/.mcfost'
        call appel_syst(cmd,syst_status)
        open(unit=1,file=trim(home)//"/.mcfost/accept_disclaimer_"//mcfost_release,status="new")
        close(unit=1)
        write(*,*) "* Thank you !                             *"
        !write(*,*) "* This screen will not appear again       *"
     else
        write(*,*) "* Exiting MCFOST                          *"
        write(*,*) "*******************************************"
        stop
     endif
  endif ! accept disclaimer

  write(*,*) "*******************************************"

  return

end subroutine display_disclaimer


!********************************************************************

subroutine save_data
  !*************************************************
  ! Si le dossier data existe on le sauve
  !*************************************************

  implicit none
  integer :: syst_status
  character(len=1024) :: cmd

  logical :: lnew_run, lmove_data
  integer :: etape, etape_start, etape_end
  character(len=512) :: local_data_dir, local_basename_data_dir

  lmove_data = .false.

  if (lonly_diff_approx) return
  if (lsed.and.(.not.ltemp)) return

  if (ldust_prop .or. ldisk_struct) then
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
          'echo " " >>  '//trim(data_dir)//'/'//trim(para)//" ; "// &
          'echo "Executed command line : '//trim(cmd_opt)//'" >> '//trim(data_dir)//'/'//trim(para)//" ; "// &
          ! date du calcul
          'date >> '//trim(data_dir)//'/'//trim(para)//" ; "// &
          ! machine de calcul
          'uname -a >> '//trim(data_dir)//'/'//trim(para)//" ; "// &
          ! id SHA
          'echo sha = '//sha_id//' >> '//trim(data_dir)//'/'//trim(para)
     call appel_syst(cmd,syst_status)

  else

     if (lmono0) then
        etape_start=1
        etape_end=1
     else
        if (ltemp.or.lsed) then
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
              write (*,*) 'Directory '//trim(local_data_dir)//' already exists : exiting!'
              stop
           else
              write (*,*) 'Directory '//trim(local_data_dir)//' already exists : backing it up in '//trim(local_data_dir)//'_old'
              if (is_dir(trim(root_dir)//"/"//trim(seed_dir)//"/"//trim(local_basename_data_dir)//'_old')) then
                 write (*,*) 'Directory '//trim(local_data_dir)//'_old already exists : exiting!'
                 stop
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
                'echo " " >>  '//trim(local_data_dir)//'/'//trim(para)//" ; "// &
                'echo "Executed command line : '//trim(cmd_opt)//'" >> '//trim(local_data_dir)//'/'//trim(para)//" ; "// &
                ! date du calcul
                'date >> '//trim(local_data_dir)//'/'//trim(para)//" ; "// &
                ! machine de calcul
                'uname -a >> '//trim(local_data_dir)//'/'//trim(para)//" ; "// &
                ! id SHA
                'echo sha = '//sha_id//' >> '//trim(local_data_dir)//'/'//trim(para)
           ! Copie du fichier lambda si besoin
           if (lsed.and.(.not.lsed_complete).and.(.not.lmono0).and.(etape==1)) then
              cmd = trim(cmd)//' ; cp '//trim(lambda_filename)//' '//trim(local_data_dir)
           endif
           call appel_syst(cmd,syst_status)
        endif
      enddo ! etape

  endif ! ldust_prop

  return

end subroutine save_data

!********************************************************************

end module init_mcfost
