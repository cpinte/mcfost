module dust_transfer

  use parameters
  use grains
  use naleat, only : seed, stream, gtype
  use dust_prop
  use temperature
  use thermal_emission
  use constants
  use scattering
  use grid
  use optical_depth
  use density
  use PAH
  use thermal_emission
  use disk_physics
  use output
  use input
  use benchmarks
  use diffusion
  use stars
  use mem
  use utils
  use ProDiMo
  use init_mcfost
  use SPH2mcfost
  use ML_ProDiMo
  use read_1d_models, only : setup_model1d_to_mcfost !to check
  use mhd2mcfost, only : setup_mhd_to_mcfost !to check
  use read_fargo3d, only : read_fargo3d_files
  use read_athena, only : read_athena_model
  use read_idefix, only : read_idefix_model
  use read_pluto, only : read_pluto_files
  use read_spherical_grid, only : read_spherical_model
  !$ use omp_lib

  implicit none

  contains

subroutine init_dust_transfer()
! Added the case where Mueller matrices are given as inputs
! 20/04/2023

  use MRW, only : initialize_cumulative_zeta
  use utils
  implicit none

#include "sprng_f.h"

  integer :: lambda_threshold
  integer, target :: lambda, lambda0
  integer, pointer :: p_lambda
  integer :: i, p_icell
  logical :: laffichage, lcompute_dust_prop
  real, allocatable, dimension(:) :: extra_heating

  ! Packet energy mise a 1
  E_paquet = 1.0_dp

  ! Compute equation (7) of Min et al (2009) for Modified Random Walk
  call initialize_cumulative_zeta()

  ! Building the wavelength & basic dust properties grid
  call init_lambda()

  if (lbenchmark_Pascucci) call init_Pascucci_benchmark()
  call init_optical_indices()

  ! Building the model volume and corresponding grid
  call order_zones()
  call define_physical_zones()

  ! Building the dust grain population
  call build_grain_size_distribution()

  laffichage=.true.

  if (lVoronoi) then
     if (lmhd_voronoi) then
        call setup_mhd_to_mcfost() !uses sph_to_voronoi
     else
        call setup_SPH2mcfost(extra_heating)
     endif
     call setup_grid()
  else ! Setting up a regular grid
     call setup_grid()
     call define_grid() ! included in setup_phantom2mcfost
     call stars_cell_indices()

     call allocate_densities()
     if (ldensity_file) then
        call read_density_file()
     else if (lread_Seb_Charnoz) then
        call densite_Seb_Charnoz()
     else if (lread_Seb_Charnoz2) then
        call densite_Seb_Charnoz2()
     else if (lfargo3d) then
        call read_fargo3d_files()
     else if (lathena) then
        call read_athena_model()
     else if (lsphere_model) then
        !on a structured spherical grid
        call read_spherical_model(density_files(1))
     else if (lmodel_1d) then !1d spherically symmetric "stellar atmosphere" models
        call setup_model1d_to_mcfost()
     else if (lidefix) then
        call read_idefix_model()
     else if (lpluto) then
        call read_pluto_files()
     else
        if (lsigma_file) call read_sigma_file()
        call define_density()
     endif

     if (lwall) call define_density_wall3D()
  endif

  call setup_scattering()
  ! Dynamic allocation of all other tables
  call alloc_dust_prop()

  call alloc_dynamique()

  stream = 0.0
  do i=1, nb_proc
     stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
  enddo

  if (lProDiMo) call setup_ProDiMo()
  if (lML) call init_ML()

  if ((ldisk_struct).and.(.not. ldust_sublimation)) then
     ! We write it later if there is sublimation
     if (lastrochem) then
        call write_disk_struct(.true.,.true.,.false.)
     else
        if (n_cells <= 1000000) then
           call write_disk_struct(.true.,lwrite_column_density,lwrite_velocity)
        else ! We do not write the density as the file is big
           call write_disk_struct(.false.,lwrite_column_density,lwrite_velocity)
        endif
     endif
  endif

  if (lmono) then ! monochromatic code: image or SED
     lambda=1
     letape_th = .false.

     if (aniso_method==1) then
        lmethod_aniso1=.true.
     else
        lmethod_aniso1=.false.
        if (laggregate) call error("you must use scattering method 1 when grains are aggregates")
     endif
     call star_energy_distribution()

     if (llimb_darkening) call read_limb_darkening_file(1)

     if (ldust_sublimation) then
        call read_sublimation_radius()
        if (.not.lVoronoi) then
           call define_grid()
           call define_dust_density()
        endif
     endif

     call prop_grains(1)
     if (lscatt_ray_tracing) then
        call alloc_ray_tracing()
        call init_directions_ray_tracing()
     endif
     call opacity(1)
     call integ_tau(1)

     if (loptical_depth_to_cell) call write_optical_depth_to_cell(1)

     write(*,*) ""
     write(*,*) "Dust properties in cell #", icell_not_empty
     p_icell = icell_not_empty
     if (aniso_method==2) write(*,*) "g             ", tab_g_pos(p_icell,1)
     write(*,*) "albedo        ", tab_albedo_pos(p_icell,1)
     if (lsepar_pola.and.(scattering_method == 2)) write(*,*) "polarisability", maxval(-tab_s12_o_s11_pos(:,p_icell,1))
     if (lopacite_only) call exit(0)

     if (l_em_disk_image) then ! the disk is emitting
        if (.not.(ldust_prop.and.lstop_after_init)) then ! we do not need the temperature if we only compute the dust prop
           call lect_Temperature()
        endif
     else ! Only the star is emitting
        Tdust=0.0
     endif !l_em_disk_image

     ! Re-apply T-based sublimation so -img matches the thermal-run dust
     if (ldust_sublimation) then
        call sublimate_dust()
        call opacity(1)
     endif

  else ! not lmono

     if (aniso_method==1) then
        lmethod_aniso1 = .true.
     else
        lmethod_aniso1 = .false.
     endif

     ! Number of steps to determine for the thermal calculation
     if (lTemp) then
        letape_th=.true.
     else
        letape_th=.false.
        if (.not.(ldust_prop.and.lstop_after_init)) then ! we do not need the temperature if we only compute the dust prop
           call lect_Temperature()
        endif
     endif
     if (lsed) then
        if (lsed_complete) then
           n_lambda2 = n_lambda
        endif
     endif

     if (lTemp.or.lsed_complete) then
        call star_energy_distribution()
        if (lISM_heating) then
           call ism_energy_distribution(ISR_model)
        else
           E_ISM = 0.0 ;
        endif

        if (lscatt_ray_tracing.and.lsed_complete) then
           call alloc_ray_tracing()
           call init_directions_ray_tracing()
        endif

        ! Try to restore dust calculation from previous run
        call read_saved_dust_prop(letape_th, lcompute_dust_prop)
        if (lcompute_dust_prop) then
           write(*,'(a30, $)') "Computing dust properties ..."
        else
           write(*,'(a46, $)') "Reading dust properties from previous run ..."
        endif

        do lambda=1,n_lambda
           if (lcompute_dust_prop) call prop_grains(lambda)
           call opacity(lambda)!_eqdiff!_data  ! ~ takes 2 seconds  PB : takes a long time in RT as using method 2 for scattering
        enddo !n
        if (lcompute_dust_prop) call save_dust_prop(letape_th)
        write(*,*) "Done"

        if (ldust_sublimation)  then
           call compute_othin_sublimation_radius()
           if (.not.lVoronoi) then
              call define_grid()
              call define_dust_density()
           endif

           if (ldisk_struct) call write_disk_struct(.false.,lwrite_column_density,lwrite_velocity)

           do lambda=1,n_lambda
              ! recalculation for opacity 2 :peut etre eviter mais implique + meme : garder tab_s11 en mem
              call prop_grains(lambda)
              call opacity(lambda)
           enddo
        endif ! ldust_sublimation

        test_tau : do lambda=1,n_lambda
           if (tab_lambda(lambda) > wl_seuil) then
              lambda_threshold=lambda
              exit test_tau
           endif
        enddo test_tau
        write(*,*) "lambda =", tab_lambda(lambda_threshold)
        call integ_tau(lambda_threshold)
        if (loptical_depth_to_cell) call write_optical_depth_to_cell(lambda_threshold)

        if (lspherical.or.l3D) then
           write(*,*) "No dark zone"
           call no_dark_zone()
           lapprox_diffusion=.false.
        else
           if (lapprox_diffusion) then
              if (lcylindrical) then
                 call define_dark_zone(lambda_threshold,p_lambda,tau_dark_zone_eq_th,.true.) ! BUG avec 1 cell
              else
                 write(*,*) "No dark zone"
                 call no_dark_zone()
              endif
           else
              write(*,*) "No dark zone"
              call no_dark_zone()
           endif
        endif

        if (lscatt_ray_tracing) then
           if (lspherical.or.l3D) then
              call no_dark_zone()
           endif
        endif

        if (lonly_diff_approx) then
           call lect_temperature()
           call Temp_approx_diffusion_vertical()
           call diffusion_approx_nLTE_nRE()
           call ecriture_temperature(2)
           return
        endif

        if (lTemp) call init_reemission(lextra_heating,extra_heating)

        !$omp parallel default(none) private(lambda) shared(n_lambda)
        !$omp do schedule(static,1)
        do lambda=1, n_lambda
           call repartition_energie(lambda)
        enddo
        !$omp end do
        !$omp end parallel

        call repartition_wl_em()
     endif ! lTemp.or.lsed_complete

     if (lTemp.and.lnRE) call init_emissivite_nRE()
  endif ! lmono

  if (laverage_grain_size) call taille_moyenne_grains()

end subroutine init_dust_transfer


!***********************************************************

subroutine recompute_opacities()
  ! Recompute grain properties and opacity tables across all wavelengths.
  ! Called when density structure or grain distribution changes during physics iterations.
  ! C. Pinte
  ! 07/2026

  implicit none

  integer :: lambda

  write(*,'(a30, $)') "Computing dust properties ..."
  do lambda = 1, n_lambda
     call prop_grains(lambda)
     call opacity(lambda)
  enddo
  write(*,*) "Done"

  return

end subroutine recompute_opacities

!***********************************************************

subroutine dust_transfer_sub()
  ! Main entry point for radiative transfer calculation
  ! C. Pinte

  implicit none

  integer :: iter, n_iter_physics
  logical :: lconverged

  call init_dust_transfer()
  if (lonly_diff_approx) return

  ! ── Thermal structure & Physics Iteration Loop ──
  if (letape_th) then
     n_iter_physics = 1
     if (lhydrostatic) n_iter_physics = 10 ! Max iterations for physics convergence

     physics_loop : do iter = 1, n_iter_physics
        call run_thermal_mc()

        lconverged = .true.

        ! Physics iteration hook: Hydrostatic equilibrium
        if (lhydrostatic) then
           write(*,*) "Updating hydrostatic equilibrium (iteration", iter, ")"
           call equilibre_hydrostatique()
           call recompute_opacities()
           lconverged = .false. ! Architectural hook for convergence criteria
        endif

        ! Additional physics iteration hooks can be added here
        ! e.g., dust sublimation iteration, grain growth, etc.

        if (lconverged) exit physics_loop
     enddo physics_loop
  endif

  ! ── Dust Species Removal ──
  if (lremove) then
     call remove_species()
     if (lTemp .and. lsed_complete) call recompute_opacities()
  endif

  ! ── Observables (Image / SED) ──
  if (lmono0) then
     call run_image_mc()
  else if (lsed) then
     call run_sed_mc()
  endif

  if (lscatt_ray_tracing) call dealloc_ray_tracing()

  return

end subroutine dust_transfer_sub

!***********************************************************

subroutine mc_photon_loop(lambda_in, p_lambda_in, n_photons2, n_phot_lim, nnfot1_start, laffichage)
  ! Monte Carlo photon loop — shared kernel used by thermal, image, and SED modes.
  ! The mode-dependent behaviour is controlled by module-level flags:
  !   letape_th  -> thermal structure (multi-wavelength, re-emission)
  !   lmono      -> monochromatic (forced scattering, image/SED)
  !   lmono0     -> image mode (fixed photon count)
  ! C. Pinte
  ! 07/2026 — Extracted from dust_transfer_sub

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: lambda_in, p_lambda_in, n_photons2, nnfot1_start
  real, intent(in) :: n_phot_lim
  logical, intent(in) :: laffichage

  ! Local variables
  real(kind=dp), dimension(4) :: Stokes
  integer :: id, icell, capt, ibar, n_photons1_cumul, n_photons2_local, nnfot1
  real :: rand
  logical :: lpacket_alive, lintersect, flag_star, flag_scatt, flag_ISM
  real(kind=dp) :: x, y, z, u, v, w
  real(kind=dp), target :: nnfot2, n_phot_sed2
  real(kind=dp), pointer :: p_nnfot2
  real(kind=dp) :: n_phot_envoyes_in_loop
  integer, target :: lambda_local
  integer, target :: p_lambda_local
  integer, pointer :: p_lambda_ptr

  ! Copy intent(in) arguments to local targets for OMP firstprivate
  lambda_local = lambda_in
  n_photons2_local = n_photons2
  ibar = 1
  n_photons1_cumul = 0
  ! Note: p_lambda_local and p_lambda_ptr are private (not firstprivate) because
  ! pointer association must be established inside the parallel region on each
  ! thread's own private targets. A firstprivate pointer would still point to the
  ! master thread's stack frame, causing cross-thread data races.

  ! OMP parallel photon transport loop
  !$omp parallel &
  !$omp default(none) &
  !$omp firstprivate(lambda_local, n_photons2_local) &
  !$omp private(p_lambda_local, p_lambda_ptr) &
  !$omp private(id,icell,lpacket_alive,lintersect,p_nnfot2,nnfot2,n_phot_sed2,n_phot_envoyes_in_loop,rand) &
  !$omp private(x,y,z,u,v,w,Stokes,flag_star,flag_ISM,flag_scatt,capt,nnfot1) &
  !$omp shared(lambda_in, p_lambda_in, nnfot1_start,n_photons_loop,capt_sup,n_phot_lim,lscatt_ray_tracing1) &
  !$omp shared(n_phot_envoyes,nb_proc,p_n_lambda_pos,n_lambda) &
  !$omp shared(stream,laffichage,lmono,lmono0,lProDiMo,lML,letape_th,tab_lambda,n_photons_lambda, n_photons1_cumul,ibar) &
  !$omp reduction(+:E_abs_nRE)
  ! Establish thread-private pointer association inside the parallel region.
  ! p_lambda_ptr is always fixed at p_lambda_in (the scattering grid index set by the
  ! caller). This matches main's effective behavior: in main, firstprivate(p_lambda)
  ! copied the pointer association to each thread, but the target (master's lambda)
  ! was never modified inside the parallel region, so p_lambda was effectively
  ! constant at its pre-parallel value in all threads.
  ! - Thermal MC: p_lambda_in = 1 (from lambda0), stays fixed even as lambda_local
  !   changes per photon via select_wl_em. Correct for scattering method 2 (p_lambda
  !   indexes the wavelength-averaged scattering matrix).
  ! - SED/image (lmono=.true.): lambda_local never changes, so p_lambda_ptr stays
  !   correctly at p_lambda_in (either current lambda or 1 for sub-sampled grid).
  p_lambda_local = p_lambda_in
  p_lambda_ptr => p_lambda_local
  if (letape_th) then
     p_nnfot2 => nnfot2
     E_abs_nRE = 0.0
  else
     if (lmono0) then
        p_nnfot2 => nnfot2
     else
        p_nnfot2 => n_phot_sed2

        if (lProDiMo.or.lML)  then
           p_nnfot2 => nnfot2  ! Constant number of packets per wavelength
           ! Increase number of packets in the UV
           if (tab_lambda(lambda_local) < 0.5) n_photons2_local = n_photons_lambda * 10
        endif
     endif
  endif

  id = 1 ! For sequential code
  !$ id = omp_get_thread_num() + 1

  !$omp do schedule(dynamic,1)
  do nnfot1=nnfot1_start,n_photons_loop
     nnfot2 = 0.0_dp
     p_nnfot2 = 0.0_dp
     n_phot_envoyes_in_loop = 0.0_dp
     n_phot_sed2 = 0.0_dp
     photon : do while ((p_nnfot2 < n_photons2_local).and.(n_phot_envoyes_in_loop < n_phot_lim))
        nnfot2=nnfot2+1.0_dp
        n_phot_envoyes(lambda_local,id) = n_phot_envoyes(lambda_local,id) + 1.0_dp
        n_phot_envoyes_in_loop = n_phot_envoyes_in_loop + 1.0_dp

        ! Wavelength choice
        if (.not.lmono) then
           rand = sprng(stream(id))
           call select_wl_em(rand,lambda_local)
        endif

        ! Packet emission
        call emit_packet(id,lambda_local, icell,x,y,z,u,v,w,stokes,flag_star,flag_ISM,lintersect)
        lpacket_alive = .true.

        ! Packet propagation
        if (lintersect) call propagate_packet(id,lambda_local,p_lambda_ptr,icell,x,y,z,u,v,w,stokes, &
             flag_star,flag_ISM,flag_scatt,lpacket_alive)

        ! The packet has now exited : on le met dans le bon capteur
        if (lpacket_alive.and.(.not.flag_ISM)) then
           call capteur(id,lambda_local,icell,x,y,z,u,v,w,Stokes,flag_star,flag_scatt,capt)
           if (capt == capt_sup) n_phot_sed2 = n_phot_sed2 + 1.0_dp ! number of photons received for step 2
        endif
     enddo photon !nnfot2

     ! Progress bar
     !$omp atomic
     n_photons1_cumul = n_photons1_cumul+1
     if (laffichage) then
        if (real(n_photons1_cumul) > 0.02*ibar * real(n_photons_loop)) then
           call progress_bar(ibar)
           !$omp atomic
           ibar = ibar+1
        endif
     endif
  enddo !nnfot1
  !$omp end do
  !$omp end parallel
  if (laffichage) call progress_bar(50)

  return

end subroutine mc_photon_loop

!***********************************************************

subroutine run_thermal_mc()
  ! Runs the thermal structure Monte Carlo calculation
  ! Iterates over non-equilibrium grains (nRE) if needed
  ! Calculates temperature structure and writes output
  ! C. Pinte
  ! 07/2026 — Extracted from dust_transfer_sub

  use radiation_field, only : reset_radiation_field
  implicit none

  integer :: n_photons2, nnfot1_start, n_iter, itime
  real :: n_phot_lim, time
  logical :: flag_em_nRE, laffichage
  integer, target :: lambda, lambda0
  integer, pointer :: p_lambda
  real(kind=dp) :: n_phot_sed2

  ! Reset energy and temperature arrays for clean iteration state
  call reset_radiation_field()
  call reset_temperature()

  laffichage = .true.
  n_photons2 = n_photons_eq_th
  n_phot_lim = 1.0e30 ! packets are not killed
  nnfot1_start = 1
  n_phot_sed2 = 0.0_dp

  lambda = 1
  lambda0 = 1
  p_lambda => lambda

  letape_th = .true.

  write(*,*) "Computing temperature structure ..."
  if (laffichage) call progress_bar(0)

  n_iter = 0
  do
     n_photons2 = n_photons_eq_th

     ! Todo : we do no need to pass lambda here
     call mc_photon_loop(lambda, p_lambda, n_photons2, n_phot_lim, nnfot1_start, laffichage)

     letape_th = .false. ! A priori, on a calcule la temperature
     if (lRE_LTE) then
        call Temp_finale()
        if (lreemission_stats) call reemission_stats()
     endif
     if (lRE_nLTE) then
        call Temp_finale_nLTE()
     endif
     if (lnRE) then
        call Temp_nRE(flag_em_nRE)
        call update_proba_abs_nRE()
        if (n_iter > 10) then
           flag_em_nRE = .true.
           write(*,*) "WARNING: Reaching the maximum number of iterations"
           write(*,*) "radiation field may not be converged"
        endif

        if (.not.flag_em_nRE) then ! need to iterate
           call emission_nRE()
           letape_th = .true.
           n_iter = n_iter + 1
           write(*,*) "Starting iteration", n_iter
        endif
     endif

     ! todo_physics : this should be moved into the dust_transfer_sub and compared with remove_species
     if (ldust_sublimation) then
        call sublimate_dust()
     endif

     if (.not.letape_th) exit
  enddo

  if (.not. lmcfost_lib) then
     call ecriture_temperature(1)
     call ecriture_sed(1)
  endif

  if (lapprox_diffusion.and.l_is_dark_zone.and.&
   (lemission_mol.or.lprodimo.or.lML.or.lforce_diff_approx.or.lemission_atom)) then
     call Temp_approx_diffusion_vertical()
     ! call Temp_approx_diffusion()
     call diffusion_approx_nLTE_nRE()
     if (.not. lmcfost_lib) call ecriture_temperature(2)
  endif

  ! Reset for the next step
  sed=0.0; sed_q=0.0 ; sed_u=0.0 ; sed_v=0.0
  n_phot_sed=0.0;  n_phot_sed2=0.0; n_phot_envoyes=0.0
  sed_star=0.0 ; sed_star_scat=0.0 ; sed_disk=0.0 ; sed_disk_scat=0.0

  call system_clock(time_end)
  if (time_end < time_begin) then
     time=(time_end + (1.0 * time_max)- time_begin)/real(time_tick)
  else
     time=(time_end - time_begin)/real(time_tick)
  endif
  if (time > 60) then
     itime = int(time)
     write (*,'(" Temperature calculation complete in ", I3, "h", I3, "m", I3, "s")')  &
          itime/3600, mod(itime/60,60), mod(itime,60)
  else
     itime = int(time)
     write (*,'(" Temperature calculation complete in ", F5.2, "s")')  time
  endif
  if (loutput_J_step1) call ecriture_J(1)

  return

end subroutine run_thermal_mc

!***********************************************************

subroutine run_image_mc()
  ! Runs monochromatic image calculation (MC photon loop + ray-tracing maps)
  ! C. Pinte
  ! 07/2026 — Extracted from dust_transfer_sub

  implicit none

  integer, target :: lambda, lambda0
  integer :: ibin, iaz, itime, i
  integer, pointer :: p_lambda
  integer :: nnfot1_start, n_photons2
  real :: n_phot_lim, time
  integer :: time_1, time_2, time_RT, time_source_fct
  logical :: laffichage

  time_source_fct = 0 ; time_RT = 0
  lmono = .true.
  E_paquet = 1.0_dp
  laffichage = .true.
  nnfot1_start = 1
  n_photons2 = n_photons_image
  n_phot_lim = 1.0e30

  if (lscattering_method1) then
     lambda = 1
     p_lambda => lambda
  else
     if (p_n_lambda_pos == n_lambda) then
        lambda = 1
        p_lambda => lambda
     else
        lambda0 = 1
        p_lambda => lambda0
     endif
  endif

  ! lambda and p_lambda are always 1 here

  if (.not.lMueller_pos_multi .and. lscatt_ray_tracing) then
     call calc_local_scattering_matrices(lambda)
  endif

  if (lspherical.or.l3D) then
     call no_dark_zone()
  else
     if (lcylindrical) call define_dark_zone(lambda,p_lambda,tau_dark_zone_obs,.false.)
  endif

  if (lweight_emission) call define_proba_weight_emission(lambda)

  call repartition_energie(lambda)

  write(*,*) "frac. energy emitted by star(s) : ", real(frac_E_stars(1))
  if (n_stars > 1) then
     write(*,*) "Relative fraction of energy emitted by each star:"
     do i=1, n_stars
        write(*,*) "Star #", i, "-->", real(prob_E_star(1,i))
     enddo
  endif

  write(*,*) "Computing MC radiation field ..."
  if (laffichage) call progress_bar(0)

  ! MC photon loop

  ! todo : lambda and p_lambda are always 1 here : it does not look like we need to pass lambda as input here
  call mc_photon_loop(lambda, p_lambda, n_photons2, n_phot_lim, nnfot1_start, laffichage)

  ! Stokes FITS output from MC packets
  if (loutput_mc) call write_stokes_fits()

  ! Ray-tracing map
  if (lscatt_ray_tracing) then

     call system_clock(time_end)
     time=(time_end - time_begin)/real(time_tick)
     if (time > 60) then
        itime = int(time)
        write (*,'(" Time = ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
     else
        write (*,'(" Time = ", F5.2, "s")')  time
     endif

     do ibin=1,RT_n_incl
        if (lscatt_ray_tracing1) then
           do iaz=1, RT_n_az
              call system_clock(time_1)
              call init_dust_source_fct1(lambda,ibin,iaz)
              call system_clock(time_2)
              time_source_fct = time_source_fct + (time_2 - time_1)

              time_1 = time_2
              call dust_map(lambda,ibin,iaz)
              if (ltau_surface) call compute_tau_surface_map(lambda,tau_surface, ibin,iaz)
              if (ltau_map) call compute_tau_map(lambda, ibin, iaz)
              call system_clock(time_2)
              time_RT = time_RT + (time_2 - time_1)
           enddo
        else
           iaz=1
           call system_clock(time_1)
           call init_dust_source_fct2(lambda,p_lambda,ibin)
           call system_clock(time_2)
           time_source_fct = time_source_fct + (time_2 - time_1)

           time_1 = time_2
           call dust_map(lambda,ibin,iaz)
           if (ltau_surface) call compute_tau_surface_map(lambda,tau_surface, ibin,iaz)
           if (ltau_map) call compute_tau_map(lambda, ibin, iaz)
           call system_clock(time_2)
           time_RT = time_RT + (time_2 - time_1)
        endif

        call system_clock(time_end)
        time=(time_end - time_begin)/real(time_tick)
        if (time > 60) then
           itime = int(time)
           write (*,'(" Time = ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
        else
           write (*,'(" Time = ", F5.2, "s")')  time
        endif
     enddo

     call ecriture_map_ray_tracing()
     if (ltau_surface) call write_tau_surface(0) ! 0 for continuum
     if (ltau_map) call write_tau_map(0)
     write(*,*) "Source fct time", time_source_fct/real(time_tick), "s"
     write(*,*) "RT time        ", time_RT/real(time_tick), "s"
  endif

  return

end subroutine run_image_mc

!***********************************************************

subroutine run_sed_mc()
  ! Runs multi-wavelength SED calculation (MC photon loop + ray-tracing + output)
  ! Handles both lsed_complete (.true. and .false.) paths internally.
  ! C. Pinte
  ! 07/2026 — Extracted from dust_transfer_sub

  implicit none

  integer, target :: lambda_target, lambda0_target, p_lambda_local
  integer, pointer :: p_lambda, p_lambda_ptr
  integer :: lambda, lambda_local, ibin, iaz, itime, nnfot1_start, n_photons2, nnfot1, n_lambda_sed, id, icell
  real :: n_phot_lim, time, tau
  real(kind=dp) :: x, y, z, u, v, w, nnfot2, lmin, lmax
  real(kind=dp), dimension(4) :: Stokes
  integer :: time_1, time_2, time_RT, time_source_fct
  logical :: laffichage, lscatt_ray_tracing1_save, lscatt_ray_tracing2_save
  logical :: flag_star, flag_ISM, flag_scatt, lpacket_alive, lintersect, lcompute_dust_prop

  time_source_fct = 0 ; time_RT = 0
  lmono = .true.
  E_paquet = 1.0_dp
  laffichage = .false.
  n_photons2 = n_photons_lambda
  n_phot_lim = n_photons_lim
  nnfot1_start = 1

  ! Setup step 2 for non-complete SED (separate wavelength grid)
  if (.not.lsed_complete) then
     lambda_target = 1
     if (ldust_prop) then
        p_lambda => lambda_target
     else
        lambda0_target = 1 ; p_lambda => lambda0_target
     endif

     call setup_scattering()
     call realloc_step2()
     call init_lambda2()
     call init_optical_indices()
     call star_energy_distribution()
     E_ISM = 0.0

     if (lscatt_ray_tracing) then
        call alloc_ray_tracing()
        call init_directions_ray_tracing()
     endif

     call read_saved_dust_prop(letape_th, lcompute_dust_prop)
     if (lcompute_dust_prop) then
        write(*,'(a30, $)') "Computing dust properties ..."
     else
        write(*,'(a46, $)') "Reading dust properties from previous run ..."
     endif
     do lambda=1,n_lambda2
        if (lcompute_dust_prop) call prop_grains(lambda)
        lambda_target = lambda
        call opacity(lambda)
     enddo
     if (lcompute_dust_prop) call save_dust_prop(letape_th)
     write(*,*) "Done"

     n_lambda_sed = n_lambda2
  else
     if (.not.lMueller_pos_multi .and. lscatt_ray_tracing) then
        call realloc_ray_tracing_scattering_matrix()
     endif

     n_lambda_sed = n_lambda
  endif

  ! Loop over SED wavelengths
  do lambda = 1, n_lambda_sed
     lambda_target = lambda
     if (lscattering_method1) then
        p_lambda => lambda_target
     else
        if (p_n_lambda_pos == n_lambda) then
           p_lambda => lambda_target
        else
           lambda0_target = 1
           p_lambda => lambda0_target
        endif
     endif

     if (.not.lMueller_pos_multi .and. lscatt_ray_tracing) then
        call calc_local_scattering_matrices(lambda)
     endif

     if (lspherical.or.l3D) then
        call no_dark_zone()
     else
        if (lcylindrical) call define_dark_zone(lambda,p_lambda,tau_dark_zone_obs,.false.)
     endif

     if (lweight_emission) call define_proba_weight_emission(lambda)

     call repartition_energie(lambda)

     if (lambda == 1) then
        write(*,*) "# Wavelength [mum]  frac. E star     tau midplane"
     endif

     ! Optical depth along midplane
     x=0.0 ; y=0.0 ; z=0.0
     Stokes = 0.0_dp ; Stokes(1) = 1.0_dp
     w = 0.0 ; u = 1.0 ; v = 0.0
     call index_cell(x,y,z, icell)
     call optical_length_tot(1,lambda,Stokes,icell,x,y,z,u,v,w,tau,lmin,lmax)
     write(*,*) "", real(tab_lambda(lambda)) ,"  ", real(frac_E_stars(lambda)), "  ", tau

     ! MC photon loop
     call mc_photon_loop(lambda, p_lambda, n_photons2, n_phot_lim, nnfot1_start, laffichage)

     ! Interstellar radiation field
     if ((.not.letape_th).and.(lProDiMo.or.lML)) then
        lscatt_ray_tracing1_save = lscatt_ray_tracing1
        lscatt_ray_tracing2_save = lscatt_ray_tracing2
        lscatt_ray_tracing1 = .false.
        lscatt_ray_tracing2 = .false.

        if (lProDiMo) call save_J_ProDiMo(lambda)
        if (lML)      call save_J_ML(lambda,.false.)

        !$omp parallel &
        !$omp default(none) &
        !$omp shared(lambda,p_lambda,n_photons_lambda,n_photons_loop,n_phot_envoyes_ISM) &
        !$omp private(id, flag_star,flag_ISM,flag_scatt,nnfot1,x,y,z,u,v,w,stokes) &
        !$omp private(lintersect,icell,lpacket_alive,nnfot2,lambda_local,p_lambda_local,p_lambda_ptr)

        p_lambda_local = p_lambda
        p_lambda_ptr => p_lambda_local

        !$omp do schedule(dynamic,1)
        do nnfot1=1,n_photons_loop
           id = 1
           !$ id = omp_get_thread_num() + 1
           nnfot2 = 0.0_dp
           photon_ISM : do while (nnfot2 < n_photons_lambda)
              n_phot_envoyes_ISM(lambda,id) = n_phot_envoyes_ISM(lambda,id) + 1.0_dp

              call emit_packet_ISM(id, icell,x,y,z,u,v,w,stokes,lintersect)
              flag_ISM = .true.
              flag_star = .false.

              if (.not.lintersect) then
                 cycle photon_ISM
              else
                 nnfot2 = nnfot2 + 1.0_dp
                 lambda_local = lambda
                 call propagate_packet(id,lambda_local,p_lambda_ptr,icell,x,y,z,u,v,w,stokes, &
                      flag_star,flag_ISM,flag_scatt,lpacket_alive)
              endif
           enddo photon_ISM
        enddo
        !$omp end do
        !$omp end parallel

        lscatt_ray_tracing1 = lscatt_ray_tracing1_save
        lscatt_ray_tracing2 = lscatt_ray_tracing2_save

        if (lML) call save_J_ML(lambda,.true.)
     endif

     ! SED ray-tracing
     if (lscatt_ray_tracing) then
        do ibin=1,RT_n_incl
           if (lscatt_ray_tracing1) then
              do iaz=1, RT_n_az
                 call system_clock(time_1)
                 call init_dust_source_fct1(lambda,ibin,iaz)
                 call system_clock(time_2)
                 time_source_fct = time_source_fct + (time_2 - time_1)

                 time_1 = time_2
                 call dust_map(lambda,ibin,iaz)
                 call system_clock(time_2)
                 time_RT = time_RT + (time_2 - time_1)
              enddo
           else
              iaz=1
              call system_clock(time_1)
              call init_dust_source_fct2(lambda,p_lambda,ibin)
              call system_clock(time_2)
              time_source_fct = time_source_fct + (time_2 - time_1)

              time_1 = time_2
              call dust_map(lambda,ibin,iaz)
              call system_clock(time_2)
              time_RT = time_RT + (time_2 - time_1)
           endif
        enddo

        if (lscatt_ray_tracing1) then
           xI_scatt = 0.0_dp
        else
           I_spec = 0.0_dp ; I_spec_star = 0.0_dp
        endif
     endif

  enddo ! lambda

  ! Write outputs after completing all wavelengths
  call ecriture_sed(2)
  if (lscatt_ray_tracing) then
     call ecriture_sed_ray_tracing()
     write(*,*) "Source fct time", time_source_fct/real(time_tick), "s"
     write(*,*) "RT time        ", time_RT/real(time_tick), "s"
  endif
  if (lProDiMo) call mcfost2ProDiMo()
  if (loutput_UV_field) call ecriture_UV_field()
  if (loutput_J) call ecriture_J(2)

  return

end subroutine run_sed_mc


!***********************************************************

subroutine emit_packet(id,lambda, icell,x0,y0,z0,u0,v0,w0,stokes,flag_star,flag_ISM,lintersect)
  ! C. Pinte
  ! 27/05/09

  integer, intent(in) :: id, lambda

  ! Packet position and direction
  integer, intent(out) :: icell
  real(kind=dp), intent(out) :: x0,y0,z0,u0,v0,w0
  real(kind=dp), dimension(4), intent(out) :: Stokes
  logical, intent(out) :: lintersect

  ! Packet properties
  logical, intent(out) :: flag_star, flag_ISM
  real :: rand, rand2, rand3, rand4
  integer :: i_star

  real(kind=dp) :: w02

  real :: hc_lk, correct_spot, cos_thet_spot, x_spot, y_spot, z_spot


  ! TODO : flag_scat et flag_direct_star, id en argument ??

  lintersect = .true.

  rand = sprng(stream(id))
  if (rand <= frac_E_stars(lambda)) then ! Emission depuis �toile
     flag_star=.true.
     flag_ISM=.false.

     rand = sprng(stream(id))
     ! Choix de l'�toile
     call select_star(lambda,rand,i_star)
     ! Emission depuis l'�toile
     rand  = sprng(stream(id))
     rand2 = sprng(stream(id))
     rand3 = sprng(stream(id))
     rand4 = sprng(stream(id))
     call emit_packet_uniform_sphere(id, i_star,rand,rand2,rand3,rand4, icell,x0,y0,z0,u0,v0,w0,w02,lintersect)
     ! Unpolarized light emanant de l'star
     Stokes(1) = E_paquet ; Stokes(2) = 0.0 ; Stokes(3) = 0.0 ; Stokes(4) = 0.0

     !********************************************************
     ! Hot spot parameters
     !********************************************************

     if (lspot) then
        !write(*,*) "*******************"
        !write(*,*) "*  Adding a spot  *"
        !write(*,*) "*******************"
        ! Not very smart, it calculates for each packet

        ! Position
        z_spot = cos(theta_spot/180.*pi)
        x_spot = sin(theta_spot/180.*pi) * cos(phi_spot/180.*pi)
        y_spot = sin(theta_spot/180.*pi) * sin(phi_spot/180.*pi)

        ! Angle subtended by the spot
        cos_thet_spot = sqrt(1.0 - surf_fraction_spot)

        ! If the photon is in the spot, the intensity is corrected
        ! Multiply by r_star car x0, y0, et z0 ont ete multiplies par r_star
        !write(*,*) "test"
        if (x_spot*x0+y_spot*y0+z_spot*z0  > cos_thet_spot * star(1)%r) then
           !  Rapport des intensites point chaud / star
           hc_lk = hp * c_light / (tab_lambda(lambda)*1e-6 * kb)
           correct_spot = (exp(hc_lk/star(1)%T) - 1)/(exp(hc_lk/T_spot) - 1)

           ! Packet energy correction
           Stokes(:) = Stokes(:) * correct_spot
        endif
     endif ! lspot

  else  if (rand <= frac_E_disk(lambda)) then! Emission from the Disk
     flag_star=.false.
     flag_ISM=.false.

     ! Initial position
     rand = sprng(stream(id))
     call select_cellule(lambda,rand, icell)

     rand  = sprng(stream(id))
     rand2 = sprng(stream(id))
     rand3 = sprng(stream(id))
     call  pos_em_cell(icell, rand,rand2,rand3,x0,y0,z0)

     ! Flight direction (uniform)
     call random_isotropic_direction(id, u0,v0,w0)

     ! Stokes parameters : lumi�re non polaris�e
     Stokes(1) = E_paquet ; Stokes(2) = 0.0 ; Stokes(3) = 0.0 ; Stokes(4) = 0.0

     if (lweight_emission) then
        Stokes(1) = Stokes(1) * correct_E_emission(icell)
     endif
  else ! Emission ISM
     flag_star=.false.
     flag_ISM=.true.
     call emit_packet_ISM(id,icell,x0,y0,z0,u0,v0,w0,stokes,lintersect)
  endif !(rand < prob_E_star)

  return

end subroutine emit_packet

!***********************************************************

subroutine propagate_packet(id,lambda,p_lambda,icell,x,y,z,u,v,w,stokes,flag_star,flag_ISM,flag_scatt,lpacket_alive)
  ! C. Pinte
  ! 27/05/09
  ! Added the case where Mueller matrices are given as inputs
  ! 20/04/2023

  ! - flag_star_direct and flag_scatt to initialize : on a besoin des 2
  ! - separate 1st scattering and the rest
  ! - lom deleted !

  use MRW, only : make_MRW_step, gamma_MRW

  integer, intent(in) :: id
  integer, intent(inout) :: lambda, p_lambda
  integer, target, intent(inout) :: icell
  real(kind=dp), intent(inout) :: x,y,z,u,v,w
  real(kind=dp), dimension(4), intent(inout) :: stokes

  logical, intent(inout) :: flag_star, flag_ISM, lpacket_alive
  logical, intent(out) :: flag_scatt

  real(kind=dp), dimension(4,4) :: M
  real(kind=dp) :: u1,v1,w1, phi, cospsi
  !real(kind=dp) :: Planck_opacity, rec_Planck_opacity, d, diff_coeff
  integer :: igrain, itheta
  integer :: n_iteractions_in_cell, icell_old
  integer, pointer :: p_icell
  real :: rand, rand2, tau, dvol

  logical :: flag_direct_star, flag_sortie

  flag_scatt = .false.
  flag_sortie = .false.

  if (flag_star) then
     flag_direct_star = .true.
  else
     flag_direct_star = .false.
  endif

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell1
  endif

  ! Loop over packet interactions:
  ! - advance the packet
  ! - make it interact with dust if needed
  n_iteractions_in_cell = 0
  infinie : do

     ! Flight length
     rand = sprng(stream(id))
     if (rand == 1.0) then
        tau=1.0e30
     else if (rand > 1.0e-6) then
        tau = -log(1.0-rand)
     else
        tau = rand
     endif

     ! Packet propagation until next interaction
     !if (.not.letape_th) then
     !   if (.not.flag_star) Stokes=0.
     !endif

     ! Do we need to perform a MRW ?
     !if ((n_iteractions_in_cell > 5)) then
     !   d = distance_to_closest_wall(icell,x,y,z)
     !   call compute_Planck_opacities(icell, Planck_opacity,rec_Planck_opacity)
     !
     !   !write(*,*) icell, d, rec_Planck_opacity,  d * rec_Planck_opacity, gamma_MRW
     !   !read(*,*)
     !
     !   do while  (d * rec_Planck_opacity > gamma_MRW )
     !      write(*,*) "MRW"
     !      call make_MRW_step(id,icell, x,y,z,Stokes(1), d, rec_Planck_opacity)
     !      d = distance_to_closest_wall(icell,x,y,z)
     !      ! todo : do we update the rec_plack_opacity ???
     !      !call compute_Planck_opacities(icell, Planck_opacity,rec_Planck_opacity)
     !   enddo
     !
     !   ! MRW is finished, we choose wl and direction
     !endif

     ! now that MRW is done, we perform the normal propagation
     icell_old = icell
     call physical_length(id,lambda,p_lambda,Stokes,icell,x,y,z,u,v,w,flag_star,flag_direct_star,tau,dvol,flag_sortie,lpacket_alive)
     if (icell == icell_old) then
        n_iteractions_in_cell = n_iteractions_in_cell + 1
        !write(*,*) "icell=", icell, n_iteractions_in_cell
     else
        n_iteractions_in_cell = 0
     endif

     if (flag_sortie) return ! Photon life finished

!     if ((icell>n_cells).and.(.not.flag_sortie)) then
!        write(*,*) "*********************"
!        write(*,*) "PB cell", icell, id, x,y,z,u,v,w
!        write(*,*) flag_star,flag_direct_star,tau,dvol,flag_sortie
!        write(*,*) "*********************"
!     endif

     ! Otherwise the photon life continues : il y a interaction
     ! Scattering or absorption
     flag_direct_star = .false.
     if (lmono) then   ! Forced scattering: multiply packet energy by albedo
        ! dark zone test
        if (l_dark_zone(icell)) then ! skip the photon
           lpacket_alive = .false.
           return
        endif

        ! Multiply by albedo
        Stokes(:)=Stokes(:)*tab_albedo_pos(p_icell,lambda)
        if (Stokes(1) < tiny_real_x1e6)then ! skip the photon
           lpacket_alive = .false.
           return
        endif

        ! Diffusion forcee: rand < albedo
        rand = -1.0
     else ! Choose absorption or scattering
        rand = sprng(stream(id))
     endif ! lmono


     if (rand < tab_albedo_pos(p_icell,lambda)) then ! Scatttering event
        flag_scatt=.true.
        flag_direct_star = .false.

        if (lscattering_method1) then ! method 1: choice of scattering grain
           rand = sprng(stream(id))
           igrain = select_scattering_grain(lambda,p_icell, rand) ! ok, not too bad, not much smaller

           rand = sprng(stream(id))
           rand2 = sprng(stream(id))
           if (lmethod_aniso1) then ! Mie phase function
              call angle_diff_theta(lambda,igrain,rand,rand2,itheta,cospsi)
              rand = sprng(stream(id))
              !  call angle_diff_phi(l,Stokes(1),Stokes(2),Stokes(3),itheta,rand,phi)
              PHI = PI * ( 2.0 * rand - 1.0 )
              ! propagation direction after scattering
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
              if (lsepar_pola) then
                 call get_Mueller_matrix_per_grain(lambda,itheta,rand2,igrain, M)
                 call update_Stokes(Stokes,u,v,w,u1,v1,w1,M)
              endif
           else ! HG phase function
              call hg(tab_g(igrain,lambda),rand, itheta, COSPSI) !HG
              if (lisotropic) then ! Isotropic scattering
                 itheta=1
                 cospsi=2.0*rand-1.0
              endif
              rand = sprng(stream(id))
              !  call angle_diff_phi(l,Stokes(1),Stokes(2),Stokes(3),itheta,rand,phi)
              PHI = PI * ( 2.0 * rand - 1.0 )
              ! propagation direction after scattering
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
           endif

        else ! method 2: scattering on the grain population
           rand = sprng(stream(id))
           rand2= sprng(stream(id))
           if (lmethod_aniso1) then ! Mie phase function
              call angle_diff_theta_pos(p_lambda,p_icell, rand, rand2, itheta, cospsi)
              if (lisotropic) then ! Isotropic scattering
                 itheta=1
                 cospsi=2.0*rand-1.0
              endif
              rand = sprng(stream(id))
              ! call angle_diff_phi(l,Stokes(1),Stokes(2),Stokes(3),itheta,rand,phi)
              PHI = PI * ( 2.0 * rand - 1.0 )
              ! propagation direction after scattering
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
              if (lsepar_pola) then
                 call get_Mueller_matrix_per_cell(lambda,itheta,rand2,p_icell, M)
                 call update_Stokes(Stokes,u,v,w,u1,v1,w1,M)
	      endif
           else ! HG phase function
              call hg(tab_g_pos(p_icell,lambda),rand, itheta, cospsi) !HG
              if (lisotropic)  then ! Isotropic scattering
                 itheta=1
                 cospsi=2.0*rand-1.0
              endif
              rand = sprng(stream(id))
              ! call angle_diff_phi(l,STOKES(1),STOKES(2),STOKES(3),itheta,rand,phi)
              phi = pi * ( 2.0 * rand - 1.0 )
              ! propagation direction after scattering
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
           endif
        endif

        ! Flight direction update
        u = u1 ; v = v1 ; w = w1

     else ! Absorption + eventual re-emission

        if ((.not.lmono).and.lnRE) then
           ! energy fraction absorbed by non-equilibrium grains
           E_abs_nRE = E_abs_nRE + Stokes(1) * (1.0_dp - proba_abs_RE(icell,lambda))
           ! Multiply by probability of absorption by grain in radiative equilibrium
           Stokes = Stokes * proba_abs_RE(icell,lambda)

           if (Stokes(1) < tiny_real)  then ! skip the photon
              lpacket_alive = .false.
              return
           endif
        endif ! lnRE

        flag_star=.false.
        flag_scatt=.false.
        flag_direct_star = .false.
        flag_ISM=.false.

        ! Wavelength choice
        if (lonly_LTE) then
           rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
           call im_reemission_LTE(id,icell,p_icell,rand,rand2,lambda)
        else if (lonly_NLTE) then
           rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
           call im_reemission_NLTE(id,icell,p_icell,rand,rand2,lambda)
        else
           ! We need to select which type of dust grain will re-emit
           rand = sprng(stream(id))
           if (rand <= Proba_abs_RE_LTE(icell,lambda)) then
              ! Case RE - LTE
              rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
              call im_reemission_LTE(id,icell,p_icell,rand,rand2,lambda)
           else if (rand <= Proba_abs_RE_LTE_p_nLTE(icell,lambda)) then
              ! Case RE - nLTE
              rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
              call im_reemission_NLTE(id,icell,p_icell,rand,rand2,lambda)
           else
              ! Case nRE - qRE
              rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
              call im_reemission_qRE(id,icell,p_icell,rand,rand2,lambda)
           endif
        endif ! only_LTE

        ! New packet direction: isotropic emission
        call random_isotropic_direction(id, u, v, w)

        ! Dust emission is unpolarised, resetimg Stokes vector
        Stokes(2)=0.0 ; Stokes(3)=0.0 ; Stokes(4)=0.0
     endif ! tab_albedo_pos

  enddo infinie

  write(*,*) "BUG propagate_packet"
  return

end subroutine propagate_packet

!***********************************************************

subroutine dust_map(lambda,ibin,iaz)
  ! Creation of the dust emission map
  ! by ray-tracing in a given direction
  ! C. Pinte
  ! 24/01/08

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: lambda, ibin, iaz
  real(kind=dp) :: u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcorner

  real(kind=dp) :: taille_pix, l, x0, y0, z0
  integer :: i,j, id, npix_x_max, n_iter_max, n_iter_min, ri_RT, phi_RT, ech_method


  integer, parameter :: n_rad_RT = 128, n_phi_RT = 30  ! OK, it works with n_rad_RT = 1000
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp) :: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  logical :: lresolved

  ! Viewing direction for ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition of image plane basis vectors in the universal frame

  ! Image x-vector without PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Image x-vector with PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo: can be simpler because rotation axis is perpendicular to x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Image y-vector with PA: orthogonal to x_plan_image and uvw
  y_plan_image = -cross_product(x_plan_image, uvw)


  if (lmono0) then
     write(*,*) "Ray-tracing ..."
     write(*,*) "i=", tab_RT_incl(ibin), "az=", tab_RT_az(iaz)
     write(*,*) "Vector to observer =", real(u),real(v),real(w)
     write(*,*) "x-image =           ", real(x_plan_image(:))
     write(*,*) "y-image =           ", real(y_plan_image(:))
  endif

  ! Initial position outside model (observer side)
  ! = image center
  l = 10.*Rmax  ! stay far away

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Method 1 = log sampling in r and uniform in phi
  ! Method 2 = linear sampling of pixels (square) with iteration on sub-pixels
  if (lsed) then
     ech_method = RT_sed_method
  else ! image
     ech_method = 2
  endif

  if (ech_method==1) then
     ! No sub-pixels because pixels are not square
     n_iter_min = 1
     n_iter_max = 1

     ! dx and dy are only required for stellar map here
     taille_pix = (map_size/zoom)  ! en AU
     dx(:) = 0.
     dy(:) = 0.

     i = 1
     j = 1
     lresolved = .false.

     rmin_RT = 0.01_dp * Rmin
     rmax_RT = 2.0_dp * Rmax

     tab_r(1) = rmin_RT
     fact_r = exp( (1.0_dp/(real(n_rad_RT,kind=dp) -1))*log(rmax_RT/rmin_RT) )

     do ri_RT = 2, n_rad_RT
        tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
     enddo

     fact_A = sqrt(pi * (fact_r - 1.0_dp/fact_r)  / n_phi_RT )

     if (l_sym_ima) then
        cst_phi = pi  / real(n_phi_RT,kind=dp)
     else
        cst_phi = two_pi  / real(n_phi_RT,kind=dp)
     endif

     ! Loop over sampling radii
     !$omp parallel &
     !$omp default(none) &
     !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
     !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j,ibin,iaz) &
     !$omp shared(n_iter_min,n_iter_max,lambda,l_sym_ima,cst_phi)
     id = 1 ! For sequential code

     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1

        r = tab_r(ri_RT)
        taille_pix =  fact_A * r ! square root of pixel area

        do phi_RT=1,n_phi_RT ! from 0 to pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)

           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image ! This is actually the center because dx = dy = 0.
           ! this is of course the expensive line:
           call intensite_pixel_dust(id,ibin,iaz,n_iter_min,n_iter_max,lambda,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

     ! We need dx and dy /= 0 for star_map now
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix
  else ! method 2: linear sampling with sub-pixels
     lresolved = .true.

     ! Vectors defining pixels (dx,dy) in the universal frame
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     ! Bottom-left corner of the image
     Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

     if (l_sym_ima) then
        npix_x_max = npix_x/2 + modulo(npix_x,2)
     else
        npix_x_max = npix_x
     endif

     ! Loop over image pixels
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,lambda,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y,n_iter_min,n_iter_max,ibin,iaz)
     id =1 ! For sequential code
     n_iter_min = 2
     n_iter_max = 6

     !$omp do schedule(dynamic,1)
     do i = 1, npix_x_max ! We only compute half the map if it is symmetric
        !$ id = omp_get_thread_num() + 1
        do j = 1,npix_y
           ! Bottom-left corner of the pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           call intensite_pixel_dust(id,ibin,iaz,n_iter_min,n_iter_max,lambda,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel
  endif ! method

  ! Adding stellar contribution
  call compute_stars_map(lambda, ibin, iaz, u,v,w, taille_pix,dx,dy, lresolved)

  id = 1 ! We add the map on the first cpu id
  Stokes_ray_tracing(lambda,:,:,ibin,iaz,1,id) = Stokes_ray_tracing(lambda,:,:,ibin,iaz,1,id) + stars_map(:,:,1)
  if (lsepar_contrib) then
     Stokes_ray_tracing(lambda,:,:,ibin,iaz,n_Stokes+1,id) = Stokes_ray_tracing(lambda,:,:,ibin,iaz,n_Stokes+1,id) &
          + stars_map(:,:,1)
  endif
  if (lsepar_pola.and.llimb_darkening) then
     Stokes_ray_tracing(lambda,:,:,ibin,iaz,2:3,id) = Stokes_ray_tracing(lambda,:,:,ibin,iaz,2:3,id) &
          + stars_map(:,:,2:3)
  endif

  if (lmono0) write(*,*) "Done"

  return

end subroutine dust_map

!***********************************************************

subroutine compute_stars_map(lambda, ibin, iaz, u,v,w, taille_pix, dx_map, dy_map, lresolved)
  ! Make a ray-traced map of the stars
  ! Also save the projected location of the stars in the map (in arcsec)

  use utils, only : interp

  integer, intent(in) :: lambda, ibin, iaz
  real(kind=dp), intent(in) :: u,v,w, taille_pix
  real(kind=dp), dimension(3), intent(in) :: dx_map, dy_map ! normalized to taille_pix
  logical, intent(in) :: lresolved

  integer, parameter :: n_ray_star_SED = 1024

  real(kind=dp), dimension(4) :: Stokes
  real(kind=dp), dimension(3) :: dx_screen, dy_screen, vec, xyz
  real(kind=dp) :: factor, factor2, lmin, lmax, norm, x, y, z, argmt, srw02, tau_avg
  real(kind=dp) :: delta, norm_screen2, offset_x, offset_y, fx, fy
  real :: cos_thet, rand, rand2, tau, pix_size, LimbDarkening, Pola_LimbDarkening, P, phi, factor_pix
  integer, dimension(n_stars) :: n_ray_star
  integer :: id, icell, iray, istar, i,j, x_center, y_center, alloc_status
  logical :: in_map, lpola, is_in_image

  integer, parameter :: nx_screen = 10
  real(kind=dp), dimension(-nx_screen:nx_screen,-nx_screen:nx_screen) :: tau_screen

  ! ToDo : this is not optimum as there can be many pixels & most of them do not contain a star
  ! allacatable array as it can be big and not fit in stack memory
  real, dimension(:,:,:), allocatable :: map_1star, Q_1star, U_1star

  stars_map(:,:,:) = 0.0
  if (n_stars < 1) return

  factor_pix = 1.0 / (taille_pix*distance)

  alloc_status = 0
  allocate(map_1star(npix_x,npix_y,nb_proc),stat=alloc_status)
  map_1star = 0.0

  alloc_status = 0

  lpola = .false.
  if (lsepar_pola.and.llimb_darkening) then
     lpola = .true.
     allocate(Q_1star(npix_x,npix_y,nb_proc),U_1star(npix_x,npix_y,nb_proc),stat=alloc_status)
     Q_1star = 0.0 ; U_1star = 0.0
  endif

  if (alloc_status/=0) call error("Allocation error in compute_stars_map")

  x_center = npix_x/2 + 1
  y_center = npix_y/2 + 1

  ! Energy
  factor = E_stars(lambda) * tab_lambda(lambda) * 1.0e-6 &
       / (distance*pc_to_AU*AU_to_Rsun)**2 * 1.35e-12

  ! Test if star is resolved
  n_ray_star(:) = max(n_ray_star_SED / n_stars,1)

  if (lresolved) then
     pix_size = map_size/zoom / max(npix_x,npix_y)
     do istar=1, n_stars
        if (2*star(istar)%r > pix_size) then
           ! on average 100 rays per pixels
           n_ray_star(istar) = max(100 * int(4*pi*(star(istar)%r/pix_size)**2), n_ray_star_SED)
           if (istar==1) write(*,*) ""
           write(*,*) "Star #",istar,"is resolved, using",n_ray_star(istar),"rays for the stellar disk"
        endif
     enddo
  endif

  do istar=1, n_stars
     ! if (star(istar)%icell == 0) cycle ! star is not in the grid ! We don't need to skip those stars anymore

     ! Compute optical depth screen in front of the star at limited resolution, e.g. 10x10
     delta = star(istar)%r / nx_screen
     norm_screen2 = 1./delta**2

     dx_screen(:) = delta * dx_map(:)/norm2(dx_map(:))
     dy_screen(:) = delta * dy_map(:)/norm2(dy_map(:))

     ! Todo : this could be made parallel, but it is very fast so far
     id = 1
     do j=-nx_screen, nx_screen
        do i=-nx_screen, nx_screen
           x = star(istar)%x + dx_screen(1) * i +  dy_screen(1) * j
           y = star(istar)%y + dx_screen(2) * i +  dy_screen(2) * j
           z = star(istar)%z + dx_screen(3) * i +  dy_screen(3) * j

           icell = star(istar)%icell
           call index_cell(x,y,z, icell)

           Stokes = 0.0_dp
           call optical_length_tot(id,lambda,Stokes,icell,x,y,z,u,v,w,tau,lmin,lmax)

           tau_screen(i,j) = tau
        enddo ! j
     enddo ! i

     map_1star(:,:,:) = 0.0
     if (lpola) then
        Q_1star(:,:,:) = 0.0
        U_1star(:,:,:) = 0.0
     endif

     norm = 0.0_dp
     tau_avg = 0.0_dp

     ! star ponctuelle
     !  x0=0.0_dp ;  y0= 0.0_dp ; z0= 0.0_dp
     !  Stokes = 0.0_dp
     !  call optical_length_tot(1,lambda,Stokes,i,j,x0,y0,z0,u,v,w,tau,lmin,lmax)
     !  Flux_etoile =  exp(-tau)
     !  write(*,*)  "F0", Flux_etoile

     ! star non ponctuelle
     !$omp parallel &
     !$omp default(none) &
     !$omp shared(stream,istar,n_ray_star,llimb_darkening,limb_darkening,mu_limb_darkening,lsepar_pola) &
     !$omp shared(pola_limb_darkening,lambda,u,v,w,tab_RT_az,lsed,star,l3D,RT_sed_method,lpola,lmono0) &
     !$omp shared(x_center,y_center,taille_pix,dx_map,dy_map,nb_proc,map_1star,Q_1star,U_1star,lresolved) &
     !$omp shared(tau_screen,dx_screen,dy_screen,norm_screen2) &
     !$omp private(id,i,j,iray,rand,rand2,x,y,z,srw02,argmt,cos_thet,LimbDarkening,Stokes,fx,fy,offset_x,offset_y,vec) &
     !$omp private(Pola_LimbDarkening,icell,tau,lmin,lmax,in_map,P,phi,is_in_image) &
     !$omp reduction(+:norm,tau_avg)
     in_map = .true. ! for SED
     LimbDarkening = 1.0

     id = 1 ! For sequential code
     !$ id = omp_get_thread_num() + 1

     !$omp do schedule(static,n_ray_star(istar)/nb_proc)
     do iray=1,n_ray_star(istar)
        ! Random position on the stellar Disk
        rand  = sprng(stream(id))
        rand2 = sprng(stream(id))

        ! Random starting position on a sphere of radius 1
        z = 2.0_dp * rand - 1.0_dp
        srw02 = sqrt(1.0-z*z)
        argmt = pi*(2.0_dp*rand2-1.0_dp)
        x = srw02 * cos(argmt)
        y = srw02 * sin(argmt)

        cos_thet = abs(x*u + y*v + z*w) ;

        if (llimb_darkening) then
           LimbDarkening = interp(limb_darkening, mu_limb_darkening, cos_thet)
           if (lsepar_pola) then
              Pola_LimbDarkening = interp(pola_limb_darkening, mu_limb_darkening, cos_thet)
           endif
        endif

        ! Random starting position on a sphere of radius r_star
        vec = (/x,y,z/) * star(istar)%r ! offset vector from center of star
        x = star(istar)%x + vec(1)
        y = star(istar)%y + vec(2)
        z = star(istar)%z + vec(3)

        ! Compute exact optical depth for each point
        !icell = star(istar)%icell
        !Stokes = 0.0_dp
        !call optical_length_tot(id,lambda,Stokes,icell,x,y,z,u,v,w,tau,lmin,lmax)

        ! Todo : we need to add a test to check if that ray will intersect another star


        ! Interpolation of optical depth : bilinear interpolation on the precomputed screen
        ! offset in in # of screen pixels (dx_screen is propto delta, so there is a delta**2 normlization)
        is_in_image = .true.
        offset_x = dot_product(vec,dx_screen) * norm_screen2 ; i = floor(offset_x) ; fx = offset_x - i
        if ((i < -nx_screen).or.(i >= nx_screen)) is_in_image = .false.

        offset_y = dot_product(vec,dy_screen) * norm_screen2 ; j = floor(offset_y) ; fy = offset_y - j
        if ((j < -nx_screen).or.(j >= nx_screen)) is_in_image = .false.

        if (is_in_image) then
           tau =  tau_screen(i,j)     * (1-fx) * (1-fy) &
                + tau_screen(i+1,j)   * fx * (1-fy) &
                + tau_screen(i,j+1)   * (1-fx) * fy &
                + tau_screen(i+1,j+1) * fx * fy
        else
           tau = 0.0_dp
        endif


        ! Average optical depth to the star
        if (lmono0) tau_avg = tau_avg + tau

        ! Coordonnees pixel
        if (lresolved) then
           call find_pixel(x,y,z, taille_pix, dx_map, dy_map, i,j,in_map)
        else
           i=1 ; j=1
        endif

        if (in_map) then
           map_1star(i,j,id) = map_1star(i,j,id) + exp(-tau) * cos_thet * LimbDarkening
           if (lpola) then ! Average polarisation in the pixel
              P = exp(-tau) * cos_thet * LimbDarkening * Pola_LimbDarkening

              ! Todo : this is temporary, only works for a star centered
              ! to be fixed : phi needs to be calculated properly
              phi = atan2((j-y_center)*1.0,(i-x_center)*1.0)

              Q_1star(i,j,id) = Q_1star(i,j,id) + P * cos(2*phi)
              U_1star(i,j,id) = U_1star(i,j,id) + P * sin(2*phi)
           endif
        endif
        norm = norm + cos_thet * LimbDarkening
     enddo ! iray
     !$omp end do
     !$omp end parallel

     ! Normalizing map and adding all the stars
     factor2 =  (factor * prob_E_star(lambda,istar)) / norm

     do id=1, nb_proc
        stars_map(:,:,1) = stars_map(:,:,1) + map_1star(:,:,id) * factor2
     enddo

     if (lpola) then
        ! Normalizing maps and adding all the stars
        stars_map(:,:,2) = stars_map(:,:,2) + sum(Q_1star(:,:,:),dim=3) * factor2
        stars_map(:,:,3) = stars_map(:,:,3) + sum(U_1star(:,:,:),dim=3) * factor2
     endif

     if (lmono0) then
        tau_avg = tau_avg/n_ray_star(istar)
        write(*,fmt='(" Optical depth from star #", i2, " is ", E12.5)') istar, tau_avg
     endif


     !---  Projected position of centres of each star
     xyz(1) = star(istar)%x ; xyz(2) = star(istar)%y ; xyz(3) = star(istar)%z

     ! Offset from map center in arcsec
     star_position(istar,ibin,iaz,1) = - dot_product(xyz, dx_map) * factor_pix ! RA negative axis
     star_position(istar,ibin,iaz,2) = dot_product(xyz, dy_map) * factor_pix

     ! Radial velocities
     star_vr(istar,ibin,iaz) = star(istar)%vx * u + star(istar)%vy * v + star(istar)%vz * w

  enddo ! n_stars

  deallocate(map_1star)
  if (lpola) deallocate(Q_1star, U_1star)

  return

end subroutine compute_stars_map

!***********************************************************

subroutine find_pixel(x,y,z,taille_pix, dx_map,dy_map, i, j, in_map)

  real(kind=dp), intent(in) :: x,y,z, taille_pix
  real(kind=dp), dimension(3), intent(in) :: dx_map, dy_map ! normalized to taille_pix
  integer, intent(out) :: i,j
  logical, intent(out) :: in_map

  real(kind=dp), dimension(3) :: xyz
  real(kind=dp) :: x_map, y_map, factor

  xyz(1) = x ; xyz(2) = y ; xyz(3) = z

  factor = 1.0 / taille_pix**2 ! dx_map length is taille_pix and we want the p[rojection in unit of taille_pix

  ! Offset from map center in units of pixel size
  x_map = dot_product(xyz, dx_map) * factor !/ taille_pix**2
  y_map = dot_product(xyz, dy_map) * factor !/ taille_pix**2

  if (modulo(npix_x,2) == 1) then
     i = nint(x_map) + npix_x/2 + 1
  else
     i = nint(x_map + 0.5) + npix_x/2
  endif
  if (modulo(npix_y,2) == 1) then
     j = nint(y_map) + npix_y/2 + 1
  else
     j = nint(y_map + 0.5) + npix_y/2
  endif

  if ((i<1).or.(i>npix_x).or.(j<1).or.(j>npix_y)) then
     in_map = .false.
  else
     in_map = .true.
  endif

  return

end subroutine find_pixel

!***********************************************************

subroutine intensite_pixel_dust(id,ibin,iaz,n_iter_min,n_iter_max,lambda,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! Calculates the intensity of a square pixel of arbitrary size, position, and orientation
  ! par une methode de Ray-tracing
  ! (u,v,w) points towards the observer
  ! TODO: Integration by Romberg method to determine the number of sub-pixels
  ! necessary
  ! Unit: W.m-2: nu.F_nu
  ! C. Pinte
  ! 12/04/07

  implicit none

  integer, intent(in) :: lambda, ibin, iaz,ipix,jpix,id, n_iter_min, n_iter_max
  real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
  real(kind=dp), intent(in) :: pixelsize,u,v,w

  real(kind=dp), dimension(N_type_flux) :: Stokes, Stokes_old

  real(kind=dp) :: x0,y0,z0,u0,v0,w0, npix2
  real(kind=dp), dimension(3) :: sdx, sdy

  real(kind=dp), parameter :: precision = 1.e-2_dp
  integer :: i, j, subpixels, ri, zj, phik, iter, icell

  logical :: lintersect

  ! TODO: there is something weird in this routine !!!

  ! Ray tracing: propagating in the other direction
  u0 = -u ; v0 = -v ; w0 = -w

  ! the number of sub-pixels in x is 2^(iter-1)
  subpixels = 1
  iter = 1

  infinie : do ! loop until convergence

     npix2 =  real(subpixels)**2
     Stokes_old(:) = Stokes(:)
     Stokes(:) = 0.0_dp

     ! Vectors defining sub-pixels
     sdx(:) = dx(:) / real(subpixels,kind=dp)
     sdy(:) = dy(:) / real(subpixels,kind=dp)

     ! The observer is outside the grid
     ri = 2*n_rad ; zj=1 ; phik=1

     ! Loop over sub-pixels that calculates the intensity at the center
     ! of each sub-pixel
     do i = 1,subpixels
        do j = 1,subpixels
           ! Center of sub-pixel
           x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
           y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
           z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)

           ! Go to the grid edge: reverse propagation
           call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)  !BUG

           if (lintersect) then ! If hitting the grid, potentially have flux
              ! Flux received in the pixel
              Stokes(:) = Stokes(:) + integ_ray_dust(lambda,icell,x0,y0,z0,u0,v0,w0)
           endif
        enddo !j
     enddo !i
     Stokes(:) = Stokes(:) / npix2

     if (iter < n_iter_min) then
        ! Iterate by default
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! Stop to avoid infinite loop
        ! write(*,*) "Warning : converging pb in ray-tracing"
        ! write(*,*) " Pixel", ipix, jpix
        exit infinie
     else
        ! Test on the difference
        if (abs(Stokes(1) - Stokes_old(1)) > precision * Stokes_old(1)) then
           ! Not converged
           subpixels = subpixels * 2
        else
           ! Converged
           exit infinie
        endif
     endif ! iter

     iter = iter + 1

  enddo infinie

  ! Take pixel surface into account (in sr)
  Stokes = Stokes * (pixelsize / (distance*pc_to_AU) )**2

  if (lsed) then
     ! Implicit summation over pixels
     Stokes_ray_tracing(lambda,ipix,jpix,ibin,iaz,:,id) = Stokes_ray_tracing(lambda,ipix,jpix,ibin,iaz,:,id) + Stokes(:)
  else
     Stokes_ray_tracing(lambda,ipix,jpix,ibin,iaz,:,id) = Stokes(:)
  endif

  return

end subroutine intensite_pixel_dust

!***********************************************************

subroutine compute_tau_surface_map(lambda,tau,ibin,iaz)

  real, intent(in) :: tau
  integer, intent(in) :: lambda, ibin, iaz
  real(kind=dp) :: u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcenter

  integer :: i,j, id, p_lambda, icell

  real :: ltot
  real(kind=dp) :: l, taille_pix, x0, y0, z0, u0, v0, w0
  logical :: lintersect, flag_star, flag_direct_star, flag_sortie, lpacket_alive

  real(kind=dp), dimension(4) :: Stokes

  p_lambda=lambda
  Stokes(1) = 1 ; Stokes(2:4) = 0.

  ! Viewing direction for ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition of image plane basis vectors in the universal frame

  ! Image x-vector without PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Image x-vector with PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo: can be simpler because rotation axis is perpendicular to x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Image y-vector with PA: orthogonal to x_plan_image and uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! Initial position outside model (observer side)
  ! = image center
  l = 10.*Rmax  ! stay far away

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0


  ! Vectors defining pixels (dx,dy) in the universal frame
  taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
  dx(:) = x_plan_image * taille_pix
  dy(:) = y_plan_image * taille_pix

  ! Bottom-left corner of the image
  Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

  ! Loop over image pixels
  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,id,Stokes,icell,lintersect,x0,y0,z0,u0,v0,w0) &
  !$omp private(flag_star,flag_direct_star,ltot,flag_sortie,lpacket_alive) &
  !$omp shared(tau,Icorner,lambda,P_lambda,pixelcenter,dx,dy,u,v,w) &
  !$omp shared(taille_pix,npix_x,npix_y,ibin,iaz,tau_surface_map,move_to_grid)
  id = 1 ! For sequential code

  !$omp do schedule(dynamic,1)
  do i = 1, npix_x
     !$ id = omp_get_thread_num() + 1
     do j = 1,npix_y
        ! Bottom-left corner of the pixel
        pixelcenter(:,id) = Icorner(:) + (i-0.5_dp) * dx(:) + (j-0.5_dp) * dy(:)

        x0 = pixelcenter(1,id)
        y0 = pixelcenter(2,id)
        z0 = pixelcenter(3,id)

        ! Ray tracing: propagating in the other direction
        u0 = -u ; v0 = -v ; w0 = -w

        ! Go to the grid edge: reverse propagation
        call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

        if (lintersect) then ! If hitting the grid, potentially have flux
           lpacket_alive = .true.
           call physical_length(id,lambda,p_lambda,Stokes,icell,x0,y0,z0,u0,v0,w0, &
                flag_star,flag_direct_star,tau,ltot,flag_sortie,lpacket_alive)
           if (flag_sortie) then ! We do not reach the surface tau=1
              tau_surface_map(i,j,ibin,iaz,:,id) = 0.0
           else
              tau_surface_map(i,j,ibin,iaz,1,id) = x0
              tau_surface_map(i,j,ibin,iaz,2,id) = y0
              tau_surface_map(i,j,ibin,iaz,3,id) = z0
           endif
        else ! We do not reach the disk
           tau_surface_map(i,j,ibin,iaz,:,id) = 0.0
        endif

     enddo !j
  enddo !i
  !$omp end do
  !$omp end parallel

  return

end subroutine compute_tau_surface_map

!***********************************************************

subroutine compute_tau_map(lambda,ibin,iaz)

  integer, intent(in) :: lambda, ibin, iaz
  real(kind=dp) :: u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcenter

  integer :: i,j, id, p_lambda, icell

  real(kind=dp) :: lmin, lmax
  real :: tau
  real(kind=dp) :: l, taille_pix, x0, y0, z0, u0, v0, w0
  logical :: lintersect, flag_star, flag_direct_star, flag_sortie, lpacket_alive

  real(kind=dp), dimension(4) :: Stokes

  p_lambda=lambda
  Stokes(1) = 1 ; Stokes(2:4) = 0.

  ! Viewing direction for ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition of image plane basis vectors in the universal frame

  ! Image x-vector without PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Image x-vector with PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo: can be simpler because rotation axis is perpendicular to x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Image y-vector with PA: orthogonal to x_plan_image and uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! Initial position outside model (observer side)
  ! = image center
  l = 10.*Rmax  ! stay far away

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0


  ! Vectors defining pixels (dx,dy) in the universal frame
  taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
  dx(:) = x_plan_image * taille_pix
  dy(:) = y_plan_image * taille_pix

  ! Bottom-left corner of the image
  Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

  ! Loop over image pixels
  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,id,Stokes,icell,lintersect,x0,y0,z0,u0,v0,w0) &
  !$omp private(flag_star,flag_direct_star,flag_sortie,lpacket_alive,tau,lmin,lmax) &
  !$omp shared(Icorner,lambda,P_lambda,pixelcenter,dx,dy,u,v,w) &
  !$omp shared(taille_pix,npix_x,npix_y,ibin,iaz,tau_map,move_to_grid)
  id = 1 ! For sequential code

  !$omp do schedule(dynamic,1)
  do i = 1, npix_x
     !$ id = omp_get_thread_num() + 1
     do j = 1,npix_y
        ! Bottom-left corner of the pixel
        pixelcenter(:,id) = Icorner(:) + (i-0.5_dp) * dx(:) + (j-0.5_dp) * dy(:)

        x0 = pixelcenter(1,id)
        y0 = pixelcenter(2,id)
        z0 = pixelcenter(3,id)

        ! Ray tracing: propagating in the other direction
        u0 = -u ; v0 = -v ; w0 = -w

        ! Go to the grid edge: reverse propagation
        call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

        if (lintersect) then ! If hitting the grid, potentially have flux
           call optical_length_tot(id,lambda,Stokes,icell,x0,y0,z0,u0,v0,w0,tau,lmin,lmax)
           tau_map(i,j,ibin,iaz,id) = tau
        else ! We do not reach the disk
           tau_map(i,j,ibin,iaz,id) = 0.0
        endif

     enddo !j
  enddo !i
  !$omp end do
  !$omp end parallel

  return

end subroutine compute_tau_map

end module dust_transfer
