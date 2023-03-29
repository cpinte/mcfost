module dust_transfer

  use parametres
  use grains
  use naleat, only : seed, stream, gtype
  use dust_prop
  use temperature
  use thermal_emission
  use constantes
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
  use read1d_models, only : setup_model1d_to_mcfost !to check
  use mhd2mcfost, only : setup_mhd_to_mcfost !to check
  use read_fargo3d, only : read_fargo3d_files
  use read_athena, only : read_athena_model
  use read_idefix, only : read_idefix_model
  use read_pluto, only : read_pluto_files
  use read_spherical_grid, only : read_spherical_model
  !$ use omp_lib

  implicit none

  contains

subroutine transfert_poussiere()

  use thermal_emission, only : frac_E_stars, frac_E_disk

  implicit none

#include "sprng_f.h"

  ! Energie des paquets
  real(kind=dp), dimension(4) :: Stokes

  ! Parametres simu
  integer :: itime, lambda_seuil, nbre_phot2
  integer :: ind_etape, first_etape_obs
  integer :: etape_start, nnfot1_start, n_iter, ibin, iaz, ibar, nnfot1_cumul

  real :: time, n_phot_lim
  logical :: lpacket_alive, lintersect

  logical :: lscatt_ray_tracing1_save, lscatt_ray_tracing2_save

  integer, target :: lambda, lambda0
  integer, pointer, save :: p_lambda
  integer :: capt

  real(kind=dp) :: x,y,z, u,v,w, lmin, lmax
  real :: rand, tau
  integer :: i, icell, p_icell, n_SPH
  logical :: flag_star, flag_scatt, flag_ISM

  logical :: laffichage, flag_em_nRE, lcompute_dust_prop

  ! Param�tres parallelisation
  integer :: id=1

  real(kind=dp), target :: nnfot2, n_phot_sed2
  real(kind=dp), pointer :: p_nnfot2
  real(kind=dp) :: n_phot_envoyes_in_loop

  integer :: time_1, time_2, time_RT, time_source_fct

  real, allocatable, dimension(:) :: extra_heating

  time_source_fct = 0 ; time_RT = 0

  lambda0 = -99 ; nnfot2=0.0_dp ; n_phot_sed2 = 0.0_dp

  ! Energie des paquets mise a 1
  E_paquet = 1.0_dp

  ! Building the wavelength & basic dust properties grid
  call init_lambda()

  if (lbenchmark_Pascucci) call init_Pascucci_benchmark()
  call init_indices_optiques()

  ! Building the model volume and corresponding grid
  call order_zones()
  call define_physical_zones()

  ! Building the dust grain population
  call build_grain_size_distribution()

  if (lphantom_file .or. lgadget2_file .or. lascii_SPH_file) then
     call setup_SPH2mcfost(density_file, limits_file, n_SPH, extra_heating)
     call setup_grid()
  else if (lmhd_voronoi) then
     call setup_mhd_to_mcfost() !uses sph_to_voronoi
     call setup_grid()
  else
     call setup_grid()
     call define_grid() ! included in setup_phantom2mcfost
     call stars_cell_indices()
  endif

  laffichage=.true.

  if (.not.(lphantom_file .or. lgadget2_file .or. lascii_SPH_file .or. lmhd_voronoi)) then ! already done by setup_SPH2mcfost
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
        call read_spherical_model(density_file)
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
  ! Allocation dynamique de tous les autres tableaux
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

  if (lmono) then ! code monochromatique
     lambda=1
     etape_i=1
     etape_f=1
     letape_th = .false.
     first_etape_obs=1

     n_phot_lim=1.0e30

     if (aniso_method==1) then
        lmethod_aniso1=.true.
     else
        lmethod_aniso1=.false.
        if (laggregate) call error("you must use scattering method 1 when grains are aggregates")
     endif
     call repartition_energie_etoiles()

     if (llimb_darkening) call read_limb_darkening_file(1)

     if (ldust_sublimation) then
        call read_sublimation_radius()
        call define_grid()
        call define_dust_density()
     endif

     call prop_grains(1)
     if (lscatt_ray_tracing) then
        call alloc_ray_tracing()
        call init_directions_ray_tracing()
     endif
     call opacite(1,1)
     call integ_tau(1)

     if (loptical_depth_to_cell) call write_optical_depth_to_cell(1)

     write(*,*) ""
     write(*,*) "Dust properties in cell #", icell_ref
     p_icell = icell_ref
     if (aniso_method==2) write(*,*) "g             ", tab_g_pos(p_icell,1)
     write(*,*) "albedo        ", tab_albedo_pos(p_icell,1)
     if (lsepar_pola.and.(scattering_method == 2)) write(*,*) "polarisability", maxval(-tab_s12_o_s11_pos(:,p_icell,1))

     if (lopacite_only) call exit(0)

     if (l_em_disk_image) then ! le disque �met
        if (.not.(ldust_prop.and.lstop_after_init)) then ! we do not need the temperature if we only compute the dust prop
           call lect_Temperature()
        endif
     else ! Seule l'�toile �met
        Tdust=0.0
     endif !l_em_disk_image

  else ! not lmono

     if (aniso_method==1) then
        lmethod_aniso1 = .true.
     else
        lmethod_aniso1 = .false.
     endif

     first_etape_obs=2
     ! Nbre d'�tapes � d�terminer pour code thermique
     if (lTemp) then
        etape_i=1
        letape_th=.true.
     else
        etape_i=2
        letape_th=.false.
        if (.not.(ldust_prop.and.lstop_after_init)) then ! we do not need the temperature if we only compute the dust prop
           call lect_Temperature()
        endif
     endif
     if (lsed) then
        if (lsed_complete) then
           etape_f=1+n_lambda
           n_lambda2 = n_lambda
        else
           etape_f=1+n_lambda2 ! modif nombre �tape
        endif
     else
        etape_f=1
     endif


     if (lTemp.or.lsed_complete) then
        call repartition_energie_etoiles()
        if (lISM_heating) then
           call repartition_energie_ISM(ISR_model)
        else
           E_ISM = 0.0 ;
        endif

        if (lscatt_ray_tracing.and.lsed_complete) then
           call alloc_ray_tracing()
           call init_directions_ray_tracing()
        endif

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

        ! Try to restore dust calculation from previous run
        call read_saved_dust_prop(letape_th, lcompute_dust_prop)
        if (lcompute_dust_prop) then
           write(*,'(a30, $)') "Computing dust properties ..."
        else
           write(*,'(a46, $)') "Reading dust properties from previous run ..."
        endif

        do lambda=1,n_lambda
           if (lcompute_dust_prop) call prop_grains(lambda)
           call opacite(lambda, p_lambda)!_eqdiff!_data  ! ~ takes 2 seconds  PB : takes a long time in RT as using method 2 for scattering
        enddo !n
        if (lcompute_dust_prop) call save_dust_prop(letape_th)
        write(*,*) "Done"

        if (ldust_sublimation)  then
           call compute_othin_sublimation_radius()
           call define_grid()
           call define_dust_density()

           if (ldisk_struct) call write_disk_struct(.false.,lwrite_column_density,lwrite_velocity)

           do lambda=1,n_lambda
              ! recalcul pour opacite 2 :peut etre eviter mais implique + meme : garder tab_s11 en mem
              call prop_grains(lambda)
              call opacite(lambda, p_lambda)
           enddo
        endif ! ldust_sublimation

        test_tau : do lambda=1,n_lambda
           if (tab_lambda(lambda) > wl_seuil) then
              lambda_seuil=lambda
              exit test_tau
           endif
        enddo test_tau
        write(*,*) "lambda =", tab_lambda(lambda_seuil)
        call integ_tau(lambda_seuil)
        if (loptical_depth_to_cell) call write_optical_depth_to_cell(lambda_seuil)

        if (lspherical.or.l3D) then
           write(*,*) "No dark zone"
           call no_dark_zone()
           lapprox_diffusion=.false.
        else
           if (lapprox_diffusion) then
              if (lcylindrical) then
                 call define_dark_zone(lambda_seuil,p_lambda,tau_dark_zone_eq_th,.true.) ! BUG avec 1 cellule
              else
                 write(*,*) "No dark zone"
                 call no_dark_zone()
              endif
           else
              write(*,*) "No dark zone"
              call no_dark_zone()
           endif
        endif

        if (lonly_diff_approx) then
           call lect_temperature()
           call Temp_approx_diffusion_vertical()
           ! call Temp_approx_diffusion()
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

  etape_start=etape_i
  nnfot1_start=1
  lambda=1 ! pour eviter depassement tab a l'initialisation
  ind_etape = etape_start

  !************************************************************
  !  Boucle principale sur les �tapes du calcul
  !************************************************************
  n_iter = 0 ! Nbre iteration grains hors equilibre
  do while (ind_etape <= etape_f)
     indice_etape=ind_etape

     if (letape_th) then ! Calcul des temperatures
        nbre_phot2 = nbre_photons_eq_th
        n_phot_lim = 1.0e30 ! on ne tue pas les paquets
     else ! calcul des observables
        ! on devient monochromatique
        lmono=.true.
        E_paquet = 1.0_dp

        !Todo: maybe we can use a variable lscatt_method2_mono
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

        if (lmono0) then ! image
           laffichage=.true.
           nbre_phot2 = nbre_photons_image
           n_phot_lim = 1.0e30 ! On ne limite pas le nbre de photons
        else ! SED
           lambda=1
           laffichage=.false.
           nbre_phot2 = nbre_photons_lambda
           n_phot_lim = nbre_photons_lim
        endif

        if ((ind_etape==first_etape_obs).and.lremove) then
           call remove_species()
           if (lTemp.and.lsed_complete) then
              write(*,'(a30, $)') "Computing dust properties ..."
              do lambda=1, n_lambda
                 call prop_grains(lambda) ! recalcul pour opacite
                 call opacite(lambda, p_lambda)
              enddo
              write(*,*) "Done"
           endif
        endif

        if ((ind_etape==first_etape_obs).and.(lsed_complete).and.(.not.lmono0)) then
           if (.not.lMueller_pos_multi .and. lscatt_ray_tracing) call realloc_ray_tracing_scattering_matrix()
        endif

        if ((ind_etape==first_etape_obs).and.(.not.lsed_complete).and.(.not.lmono0)) then ! Changement des lambda
           ! if we reallocate, we are now in monochromatic
           ! except if we want to save the dust properties
           if (ldust_prop) then
              p_lambda => lambda
           else
              lambda0 = 1 ; p_lambda => lambda0
           endif

           ! reorganisation memoire
           call realloc_step2()

           call init_lambda2()
           call init_indices_optiques()

           call repartition_energie_etoiles()
           E_ISM = 0.0 ! ISM done a second step in SED step2 calculation

           if (lscatt_ray_tracing) then
              call alloc_ray_tracing()
              call init_directions_ray_tracing()
           endif

           ! Recalcul des propri�t�s optiques
           ! Try to restore dust calculation from previous run
           call read_saved_dust_prop(letape_th, lcompute_dust_prop)
           if (lcompute_dust_prop) then
              write(*,'(a30, $)') "Computing dust properties ..."
           else
              write(*,'(a46, $)') "Reading dust properties from previous run ..."
           endif
           do lambda=1,n_lambda2
              if (lcompute_dust_prop) call prop_grains(lambda)
              call opacite(lambda, p_lambda)!_eqdiff!_data  ! ~ takes 2 seconds
           enddo !n
           if (lcompute_dust_prop) call save_dust_prop(letape_th)
           write(*,*) "Done"
        endif ! ind_etape==first_etape_obs

        lambda = ind_etape - first_etape_obs + 1

        if (.not.lMueller_pos_multi .and. lscatt_ray_tracing) call calc_local_scattering_matrices(lambda, p_lambda)

        if (lspherical.or.l3D) then
           call no_dark_zone()
        else
           if (lcylindrical) call define_dark_zone(lambda,p_lambda,tau_dark_zone_obs,.false.)
        endif
        !call no_dark_zone()
        ! n_dif_max = seuil_n_dif(lambda)

        if (lweight_emission) call define_proba_weight_emission(lambda)

        call repartition_energie(lambda)
        if (lmono0) then
           write(*,*) "frac. energy emitted by star(s) : ", real(frac_E_stars(1))
           if (n_etoiles > 1) then
              write(*,*) "Relative fraction of energy emitted by each star:"
              do i=1, n_etoiles
                 write(*,*) "Star #", i, "-->", real(prob_E_star(1,i))
              enddo
           endif
        endif

     endif !letape_th

     if (ind_etape==etape_start) then
        call system_clock(time_end)

        time=(time_end - time_begin)/real(time_tick)
         if (time > 60) then
            itime = int(time)
            write (*,'(" Initialization complete in ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
         else
            write (*,'(" Initialization complete in ", F5.2, "s")')  time
         endif
      endif

     if (letape_th) write(*,*) "Computing temperature structure ..."
     if (lmono0)    write(*,*) "Computing MC radiation field ..."

     if (laffichage) call progress_bar(0)

     if ((ind_etape >= first_etape_obs).and.(.not.lmono0)) then
        if (ind_etape == first_etape_obs) write(*,*) "# Wavelength [mum]  frac. E star     tau midplane"
        ! Optical depth along midplane
        x=0.0 ; y=0.0 ; z=0.0
        Stokes = 0.0_dp ; Stokes(1) = 1.0_dp
        w = 0.0 ; u = 1.0 ; v = 0.0
        call indice_cellule(x,y,z, icell)
        call optical_length_tot(1,lambda,Stokes,icell,x,y,y,u,v,w,tau,lmin,lmax)
        write(*,*) "", real(tab_lambda(lambda)) ,"  ", real(frac_E_stars(lambda)), "  ", tau
     endif


     ! Les pointeurs (meme les tab) doivent �tre priv�s !!! COPYIN
     !$omp parallel &
     !$omp default(none) &
     !$omp firstprivate(lambda,p_lambda) &
     !$omp private(id,icell,lpacket_alive,lintersect,p_nnfot2,nnfot2,n_phot_envoyes_in_loop,rand) &
     !$omp private(x,y,z,u,v,w,Stokes,flag_star,flag_ISM,flag_scatt,n_phot_sed2,capt) &
     !$omp shared(nnfot1_start,nbre_photons_loop,capt_sup,n_phot_lim,lscatt_ray_tracing1) &
     !$omp shared(nbre_phot2,n_phot_envoyes,nb_proc) &
     !$omp shared(stream,laffichage,lmono,lmono0,lProDiMo,lML,letape_th,tab_lambda,nbre_photons_lambda, nnfot1_cumul,ibar) &
     !$omp reduction(+:E_abs_nRE)
     if (letape_th) then
        p_nnfot2 => nnfot2
        E_abs_nRE = 0.0
     else
        if (lmono0) then
           p_nnfot2 => nnfot2
        else
           p_nnfot2 => n_phot_sed2

           if (lProDiMo.or.lML)  then
              p_nnfot2 => nnfot2  ! Nbre de paquet cst par lambda
              ! Augmentation du nbre de paquets dans UV
              if (tab_lambda(lambda) < 0.5) nbre_phot2 = nbre_photons_lambda * 10
           endif
        endif
     endif

     id = 1 ! Pour code sequentiel
     !$ id = omp_get_thread_num() + 1
     ibar=1 ;  nnfot1_cumul = 0

     !$omp do schedule(dynamic,1)
     do nnfot1=nnfot1_start,nbre_photons_loop
        p_nnfot2 = 0.0_dp
        n_phot_envoyes_in_loop = 0.0_dp
        photon : do while ((p_nnfot2 < nbre_phot2).and.(n_phot_envoyes_in_loop < n_phot_lim))
           nnfot2=nnfot2+1.0_dp
           n_phot_envoyes(lambda,id) = n_phot_envoyes(lambda,id) + 1.0_dp
           n_phot_envoyes_in_loop = n_phot_envoyes_in_loop + 1.0_dp

           ! Choix longueur d'onde
           if (.not.lmono) then
              rand = sprng(stream(id))
              call select_wl_em(rand,lambda)
           endif

           ! Emission du paquet
           call emit_packet(id,lambda, icell,x,y,z,u,v,w,stokes,flag_star,flag_ISM,lintersect)
           lpacket_alive = .true.

           ! Propagation du packet
           if (lintersect) call propagate_packet(id,lambda,p_lambda,icell,x,y,z,u,v,w,stokes, &
                flag_star,flag_ISM,flag_scatt,lpacket_alive)

           ! La paquet est maintenant sorti : on le met dans le bon capteur
           if (lpacket_alive.and.(.not.flag_ISM)) then
              call capteur(id,lambda,icell,x,y,z,u,v,w,Stokes,flag_star,flag_scatt,capt)
              if (capt == capt_sup) n_phot_sed2 = n_phot_sed2 + 1.0_dp ! nbre de photons recus pour etape 2
           endif
        enddo photon !nnfot2

        ! Progress bar
        !$omp atomic
        nnfot1_cumul = nnfot1_cumul+1
        if (laffichage) then
           if (real(nnfot1_cumul) > 0.02*ibar * real(nbre_photons_loop)) then
              call progress_bar(ibar)
              !$omp atomic
              ibar = ibar+1
           endif
        endif
     enddo !nnfot1
     !$omp end do
     !$omp end parallel
     if (laffichage) call progress_bar(50)

     ! Champ de radiation interstellaire
     if ((.not.letape_th).and.(lProDiMo.or.lML)) then
        ! Pas de ray-tracing avec les packets ISM
        lscatt_ray_tracing1_save = lscatt_ray_tracing1
        lscatt_ray_tracing2_save = lscatt_ray_tracing2
        lscatt_ray_tracing1 = .false.
        lscatt_ray_tracing2 = .false.

        ! Sauvegarde champ stellaire et th separement
        if (lProDiMo) call save_J_ProDiMo(lambda)
        if (lML)      call save_J_ML(lambda,.false.)

        !$omp parallel &
        !$omp default(none) &
        !$omp shared(lambda,p_lambda,nbre_photons_lambda,nbre_photons_loop,n_phot_envoyes_ISM) &
        !$omp private(id, flag_star,flag_ISM,flag_scatt,nnfot1,x,y,z,u,v,w,stokes,lintersect,icell,lpacket_alive,nnfot2)

        flag_star = .false.

        !$omp do schedule(dynamic,1)
        do nnfot1=1,nbre_photons_loop
           !$ id = omp_get_thread_num() + 1
           nnfot2 = 0.0_dp
           photon_ISM : do while (nnfot2 < nbre_photons_lambda)
              n_phot_envoyes_ISM(lambda,id) = n_phot_envoyes_ISM(lambda,id) + 1.0_dp

              ! Emission du paquet
              call emit_packet_ISM(id, icell,x,y,z,u,v,w,stokes,lintersect)
              flag_ISM = .true.

              ! Le photon sert a quelquechose ou pas ??
              if (.not.lintersect) then
                 cycle photon_ISM
              else
                 nnfot2 = nnfot2 + 1.0_dp
                 ! Propagation du packet
                 call propagate_packet(id,lambda,p_lambda,icell,x,y,z,u,v,w,stokes,flag_star,flag_ISM,flag_scatt,lpacket_alive)
              endif
           enddo photon_ISM ! nnfot2
        enddo ! nnfot1
        !$omp end do
        !$omp end parallel

        lscatt_ray_tracing1 = lscatt_ray_tracing1_save
        lscatt_ray_tracing2 = lscatt_ray_tracing2_save

        if (lML) call save_J_ML(lambda,.true.)
     endif ! champ ISM

     !----------------------------------------------------
     if (lmono0) then ! Creation image
        if (loutput_mc) call write_stokes_fits()

        ! Carte ray-tracing
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
                    call dust_map(lambda,ibin,iaz) ! Ne prend pas de temps en SED
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
                 call dust_map(lambda,ibin,iaz) ! Ne prend pas de temps en SED
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
        endif

     elseif (letape_th) then ! Calcul de la structure en temperature

        letape_th=.false. ! A priori, on a calcule la temperature
        if (lRE_LTE) then
           call Temp_finale()
           if (lreemission_stats) call reemission_stats()
        end if
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

           if (.not.flag_em_nRE) then ! il faut iterer
              call emission_nRE()
              letape_th=.true.
              first_etape_obs = first_etape_obs + 1
              etape_f = etape_f + 1
              n_iter = n_iter + 1
              write(*,*) "Starting iteration", n_iter
           endif
        endif

        if (ldust_sublimation) then
           call sublimate_dust()
        endif

        ! A-t-on fini le calcul des grains hors eq ?
        if (.not.letape_th) then ! oui, on passe a la suite
           call ecriture_temperature(1)
           call ecriture_sed(1)

           if (lapprox_diffusion.and.l_is_dark_zone.and.(lemission_mol.or.lprodimo.or.lML.or.lforce_diff_approx)) then
              call Temp_approx_diffusion_vertical()
              ! call Temp_approx_diffusion()
              call diffusion_approx_nLTE_nRE()
              call ecriture_temperature(2)
           endif

           ! Remise a zero pour etape suivante
           sed=0.0; sed_q=0.0 ; sed_u=0.0 ; sed_v=0.0
           n_phot_sed=0.0;  n_phot_sed2=0.0; n_phot_envoyes=0.0
           sed_star=0.0 ; sed_star_scat=0.0 ; sed_disk=0.0 ; sed_disk_scat=0.0
        endif ! .not.letape_th

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
           write (*,'(" Temperature calculation complete in ", F5.2, "s")')  time
        endif
        if (loutput_J_step1) call ecriture_J(1)

     else ! Etape 2 SED

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
                    call dust_map(lambda,ibin,iaz) ! Ne prend pas de temps en SED
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
                 call dust_map(lambda,ibin,iaz) ! Ne prend pas de temps en SED
                 call system_clock(time_2)
                 time_RT = time_RT + (time_2 - time_1)
              endif
           enddo

           ! Pour longeur d'onde suivante
           if (lscatt_ray_tracing1) then
              xI_scatt = 0.0_dp
           else
              I_spec = 0.0_dp ; I_spec_star = 0.0_dp
           endif
        endif

        if (ind_etape==etape_f) then ! Ecriture SED ou spectre
           call ecriture_sed(2)
           if (lscatt_ray_tracing) call ecriture_sed_ray_tracing()
           if (lProDiMo) call mcfost2ProDiMo()
           if (loutput_UV_field) call ecriture_UV_field()
           if (loutput_J) call ecriture_J(2)
        endif

     endif

     ind_etape = ind_etape + 1
  enddo ! nbre_etapes

  if (lscatt_ray_tracing) then
     call dealloc_ray_tracing()
     write(*,*) "Source fct time", time_source_fct/real(time_tick), "s"
     write(*,*) "RT time        ", time_RT/real(time_tick), "s"
  endif

  return

end subroutine transfert_poussiere

!***********************************************************

subroutine emit_packet(id,lambda, icell,x0,y0,z0,u0,v0,w0,stokes,flag_star,flag_ISM,lintersect)
  ! C. Pinte
  ! 27/05/09

  integer, intent(in) :: id, lambda

  ! Position et direction du packet
  integer, intent(out) :: icell
  real(kind=dp), intent(out) :: x0,y0,z0,u0,v0,w0
  real(kind=dp), dimension(4), intent(out) :: Stokes
  logical, intent(out) :: lintersect

  ! Proprietes du packet
  logical, intent(out) :: flag_star, flag_ISM
  real :: rand, rand2, rand3, rand4
  integer :: i_star

  real(kind=dp) :: w02, srw02
  real :: argmt

  real :: hc_lk, correct_spot, cos_thet_spot, x_spot, y_spot, z_spot


  ! TODO : flag_scat et flag_direct_star, id en argument ??

  lintersect = .true.

  rand = sprng(stream(id))
  if (rand <= frac_E_stars(lambda)) then ! Emission depuis �toile
     flag_star=.true.
     flag_ISM=.false.

     rand = sprng(stream(id))
     ! Choix de l'�toile
     call select_etoile(lambda,rand,i_star)
     ! Emission depuis l'�toile
     rand  = sprng(stream(id))
     rand2 = sprng(stream(id))
     rand3 = sprng(stream(id))
     rand4 = sprng(stream(id))
     call em_sphere_uniforme(id, i_star,rand,rand2,rand3,rand4, icell,x0,y0,z0,u0,v0,w0,w02,lintersect)
     ! Lumiere non polarisee emanant de l'etoile
     Stokes(1) = E_paquet ; Stokes(2) = 0.0 ; Stokes(3) = 0.0 ; Stokes(4) = 0.0

     !********************************************************
     ! Parametres du point chaud
     !********************************************************

     if (lspot) then
        !write(*,*) "*******************"
        !write(*,*) "*  Adding a spot  *"
        !write(*,*) "*******************"
        ! Pas tres malin ca, ca fait les calculs a chaque paquet

        ! Position
        z_spot = cos(theta_spot/180.*pi)
        x_spot = sin(theta_spot/180.*pi) * cos(phi_spot/180.*pi)
        y_spot = sin(theta_spot/180.*pi) * sin(phi_spot/180.*pi)

        ! Angle sous-tendu par le spot
        cos_thet_spot = sqrt(1.0 - surf_fraction_spot)

        ! Si le photon est dans le spot, on corrige l'intensite
        ! On multiplis par r_star car x0, y0, et z0 ont ete multiplies par r_star
        !write(*,*) "test"
        if (x_spot*x0+y_spot*y0+z_spot*z0  > cos_thet_spot * etoile(1)%r) then
           !  Rapport des intensites point chaud / etoile
           hc_lk = hp * c_light / (tab_lambda(lambda)*1e-6 * kb)
           correct_spot = (exp(hc_lk/etoile(1)%T) - 1)/(exp(hc_lk/T_spot) - 1)

           ! Correction energy packet
           Stokes(:) = Stokes(:) * correct_spot
        endif
     endif ! lspot

  else  if (rand <= frac_E_disk(lambda)) then! Emission depuis le disque
     flag_star=.false.
     flag_ISM=.false.

     ! Position initiale
     rand = sprng(stream(id))
     call select_cellule(lambda,rand, icell)

     rand  = sprng(stream(id))
     rand2 = sprng(stream(id))
     rand3 = sprng(stream(id))
     call  pos_em_cellule(icell, rand,rand2,rand3,x0,y0,z0)

     ! Direction de vol (uniforme)
     rand = sprng(stream(id))
     W0 = 2.0 * rand - 1.0
     W02 =  1.0 - W0*W0
     SRW02 = sqrt (  W02 )
     rand = sprng(stream(id))
     ARGMT = PI * ( 2.0 * rand - 1.0 )
     U0 = SRW02 * cos(ARGMT)
     V0 = SRW02 * sin(ARGMT)

     ! Parametres de stokes : lumi�re non polaris�e
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

  ! - flag_star_direct et flag_scatt a initialiser : on a besoin des 2
  ! - separer 1ere diffusion et reste
  ! - lom supprime !

  integer, intent(in) :: id
  integer, intent(inout) :: lambda, p_lambda
  integer, target, intent(inout) :: icell
  real(kind=dp), intent(inout) :: x,y,z,u,v,w
  real(kind=dp), dimension(4), intent(inout) :: stokes

  logical, intent(inout) :: flag_star, flag_ISM, lpacket_alive
  logical, intent(out) :: flag_scatt

  real(kind=dp) :: u1,v1,w1, phi, cospsi, w02, srw02, argmt
  integer :: taille_grain, itheta
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
     p_icell => icell_ref
  endif

  ! Boucle sur les interactions du paquets:
  ! - on avance le paquet
  ! - on le fait interagir avec la poussiere si besoin
  infinie : do

     ! Longueur de vol
     rand = sprng(stream(id))
     if (rand == 1.0) then
        tau=1.0e30
     else if (rand > 1.0e-6) then
        tau = -log(1.0-rand)
     else
        tau = rand
     endif

     ! Propagation du packet jusqu'a la prochaine interaction
     !if (.not.letape_th) then
     !   if (.not.flag_star) Stokes=0.
     !endif
     call physical_length(id,lambda,p_lambda,Stokes,icell,x,y,z,u,v,w,flag_star,flag_direct_star,tau,dvol,flag_sortie,lpacket_alive)
     if (flag_sortie) return ! Vie du photon terminee

!     if ((icell>n_cells).and.(.not.flag_sortie)) then
!        write(*,*) "*********************"
!        write(*,*) "PB cell", icell, id, x,y,z,u,v,w
!        write(*,*) flag_star,flag_direct_star,tau,dvol,flag_sortie
!        write(*,*) "*********************"
!     endif

     ! Sinon la vie du photon continue : il y a interaction
     ! Diffusion ou absorption
     flag_direct_star = .false.
     if (lmono) then   ! Diffusion forcee : on multiplie l'energie du packet par l'albedo
        ! test zone noire
        if (l_dark_zone(icell)) then ! on saute le photon
           lpacket_alive = .false.
           return
        endif

        ! Multiplication par albedo
        Stokes(:)=Stokes(:)*tab_albedo_pos(p_icell,lambda)
        if (Stokes(1) < tiny_real_x1e6)then ! on saute le photon
           lpacket_alive = .false.
           return
        endif

        ! Diffusion forcee: rand < albedo
        rand = -1.0
     else ! Choix absorption ou diffusion
        rand = sprng(stream(id))
     endif ! lmono


     if (rand < tab_albedo_pos(p_icell,lambda)) then ! Diffusion
        flag_scatt=.true.
        flag_direct_star = .false.

        if (lscattering_method1) then ! methode 1 : choix du grain diffuseur
           rand = sprng(stream(id))
           taille_grain = select_scattering_grain(lambda,p_icell, rand) ! ok, not too bad, not much smaller

           rand = sprng(stream(id))
           rand2 = sprng(stream(id))
           if (lmethod_aniso1) then ! fonction de phase de Mie
              call angle_diff_theta(lambda,taille_grain,rand,rand2,itheta,cospsi)
              rand = sprng(stream(id))
              !  call angle_diff_phi(l,Stokes(1),Stokes(2),Stokes(3),itheta,rand,phi)
              PHI = PI * ( 2.0 * rand - 1.0 )
              ! direction de propagation apres diffusion
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
              if (lsepar_pola) then
                 ! Nouveaux param�tres de Stokes
                 if (laggregate) then
                    call new_stokes_gmm(lambda,itheta,rand2,taille_grain,u,v,w,u1,v1,w1,stokes)
                 else
                    call new_stokes(lambda,itheta,rand2,taille_grain,u,v,w,u1,v1,w1,stokes)
                 endif
              endif
           else ! fonction de phase HG
              call hg(tab_g(taille_grain,lambda),rand, itheta, COSPSI) !HG
              if (lisotropic) then ! Diffusion isotrope
                 itheta=1
                 cospsi=2.0*rand-1.0
              endif
              rand = sprng(stream(id))
              !  call angle_diff_phi(l,Stokes(1),Stokes(2),Stokes(3),itheta,rand,phi)
              PHI = PI * ( 2.0 * rand - 1.0 )
              ! direction de propagation apres diffusion
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
              ! Param�tres de Stokes non modifi�s
           endif

        else ! methode 2 : diffusion sur la population de grains
           rand = sprng(stream(id))
           rand2= sprng(stream(id))
           if (lmethod_aniso1) then ! fonction de phase de Mie
              call angle_diff_theta_pos(p_lambda,p_icell, rand, rand2, itheta, cospsi)
              if (lisotropic) then ! Diffusion isotrope
                 itheta=1
                 cospsi=2.0*rand-1.0
              endif
              rand = sprng(stream(id))
              ! call angle_diff_phi(l,Stokes(1),Stokes(2),Stokes(3),itheta,rand,phi)
              PHI = PI * ( 2.0 * rand - 1.0 )
              ! direction de propagation apres diffusion
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
              ! Nouveaux param�tres de Stokes
              if (lsepar_pola) call new_stokes_pos(p_lambda,itheta,rand2,p_icell,u,v,w,u1,v1,w1,Stokes)
           else ! fonction de phase HG
              call hg(tab_g_pos(p_icell,lambda),rand, itheta, cospsi) !HG
              if (lisotropic)  then ! Diffusion isotrope
                 itheta=1
                 cospsi=2.0*rand-1.0
              endif
              rand = sprng(stream(id))
              ! call angle_diff_phi(l,STOKES(1),STOKES(2),STOKES(3),itheta,rand,phi)
              phi = pi * ( 2.0 * rand - 1.0 )
              ! direction de propagation apres diffusion
              call cdapres(cospsi, phi, u, v, w, u1, v1, w1)
              ! Param�tres de Stokes non modifi�s
           endif
        endif

        ! Mise a jour direction de vol
        u = u1 ; v = v1 ; w = w1

     else ! Absorption

        if ((.not.lmono).and.lnRE) then
           ! fraction d'energie absorbee par les grains hors equilibre
           E_abs_nRE = E_abs_nRE + Stokes(1) * (1.0_dp - proba_abs_RE(icell,lambda))
           ! Multiplication par proba abs sur grain en eq. radiatif
           Stokes = Stokes * proba_abs_RE(icell,lambda)

           if (Stokes(1) < tiny_real)  then ! on saute le photon
              lpacket_alive = .false.
              return
           endif
        endif ! lnRE

        flag_star=.false.
        flag_scatt=.false.
        flag_direct_star = .false.
        flag_ISM=.false.

        ! Choix longueur d'onde
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
              ! Cas RE - LTE
              rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
              call im_reemission_LTE(id,icell,p_icell,rand,rand2,lambda)
           else if (rand <= Proba_abs_RE_LTE_p_nLTE(icell,lambda)) then
              ! Cas RE - nLTE
              rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
              call im_reemission_NLTE(id,icell,p_icell,rand,rand2,lambda)
           else
              ! Cas nRE - qRE
              rand = sprng(stream(id)) ; rand2 = sprng(stream(id))
              call im_reemission_qRE(id,icell,p_icell,rand,rand2,lambda)
           endif
        endif ! only_LTE

        ! Nouvelle direction de vol : emission isotrope
        rand = sprng(stream(id))
        w = 2.0 * rand - 1.0
        w02 =  1.0 - w*w
        srw02 = sqrt (w02)
        rand = sprng(stream(id))
        argmt = pi * ( 2.0 * rand - 1.0 )
        u = srw02 * cos(argmt)
        v = srw02 * sin(argmt)

        ! Emission non polaris�e : remise � 0 des parametres de Stokes
        Stokes(2)=0.0 ; Stokes(3)=0.0 ; Stokes(4)=0.0
     endif ! tab_albedo_pos

  enddo infinie

  write(*,*) "BUG propagate_packet"
  return

end subroutine propagate_packet

!***********************************************************

subroutine dust_map(lambda,ibin,iaz)
  ! Creation de la carte d'emission de la poussiere
  ! par ray-tracing dans une direction donnee
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


  integer, parameter :: n_rad_RT = 128, n_phi_RT = 30  ! OK, ca marche avec n_rad_RT = 1000
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp) :: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  logical :: lresolved

  ! Direction de visee pour le ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition des vecteurs de base du plan image dans le repere universel

  ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Vecteur x image avec PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
  y_plan_image = -cross_product(x_plan_image, uvw)


  if (lmono0) then
     write(*,*) "Ray-tracing ..."
     write(*,*) "i=", tab_RT_incl(ibin), "az=", tab_RT_az(iaz)
     write(*,*) "Vector to observer =", real(u),real(v),real(w)
     write(*,*) "x-image =           ", real(x_plan_image(:))
     write(*,*) "y-image =           ", real(y_plan_image(:))
  endif

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Methode 1 = echantillonage log en r et uniforme en phi
  ! Methode 2 = echantillonage lineaire des pixels (carres donc) avec iteration sur les sous-pixels
  if (lsed) then
     ech_method = RT_sed_method
  else ! image
     ech_method = 2
  endif

  if (ech_method==1) then
     ! Pas de sous-pixel car les pixels ne sont pas carres
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
        cst_phi = deux_pi  / real(n_phi_RT,kind=dp)
     endif

     ! Boucle sur les rayons d'echantillonnage
     !$omp parallel &
     !$omp default(none) &
     !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
     !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j,ibin,iaz) &
     !$omp shared(n_iter_min,n_iter_max,lambda,l_sym_ima,cst_phi)
     id = 1 ! pour code sequentiel

     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1

        r = tab_r(ri_RT)
        taille_pix =  fact_A * r ! racine carree de l'aire du pixel

        do phi_RT=1,n_phi_RT ! de 0 a pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)

           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image ! C'est le centre en fait car dx = dy = 0.
           ! this is of course the expensive line:
           call intensite_pixel_dust(id,ibin,iaz,n_iter_min,n_iter_max,lambda,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

     ! We need dx and dy /= 0 for star_map now
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix
  else ! method 2 : echantillonnage lineaire avec sous-pixels
     lresolved = .true.

     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     ! Coin en bas gauche de l'image
     Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

     if (l_sym_ima) then
        npix_x_max = npix_x/2 + modulo(npix_x,2)
     else
        npix_x_max = npix_x
     endif

     ! Boucle sur les pixels de l'image
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,lambda,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y,n_iter_min,n_iter_max,ibin,iaz)
     id =1 ! pour code sequentiel
     n_iter_min = 2
     n_iter_max = 6

     !$omp do schedule(dynamic,1)
     do i = 1, npix_x_max ! We only compute half the map if it is symmetric
        !$ id = omp_get_thread_num() + 1
        do j = 1,npix_y
           ! Coin en bas gauche du pixel
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
  real(kind=dp) :: facteur, facteur2, lmin, lmax, norme, x, y, z, argmt, srw02, tau_avg
  real(kind=dp) :: delta, norm_screen2, offset_x, offset_y, fx, fy
  real :: cos_thet, rand, rand2, tau, pix_size, LimbDarkening, Pola_LimbDarkening, P, phi, factor_pix
  integer, dimension(n_etoiles) :: n_ray_star
  integer :: id, icell, iray, istar, i,j, x_center, y_center, alloc_status
  logical :: in_map, lpola, is_in_image

  integer, parameter :: nx_screen = 10
  real(kind=dp), dimension(-nx_screen:nx_screen,-nx_screen:nx_screen) :: tau_screen

  ! ToDo : this is not optimum as there can be many pixels & most of them do not contain a star
  ! allacatable array as it can be big and not fit in stack memory
  real, dimension(:,:,:), allocatable :: map_1star, Q_1star, U_1star

  stars_map(:,:,:) = 0.0
  if (n_etoiles < 1) return

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

  ! Energie
  facteur = E_stars(lambda) * tab_lambda(lambda) * 1.0e-6 &
       / (distance*pc_to_AU*AU_to_Rsun)**2 * 1.35e-12

  ! Test si etoile est resolue
  n_ray_star(:) = max(n_ray_star_SED / n_etoiles,1)

  if (lresolved) then
     pix_size = map_size/zoom / max(npix_x,npix_y)
     do istar=1, n_etoiles
        if (2*etoile(istar)%r > pix_size) then
           ! on average 100 rays per pixels
           n_ray_star(istar) = max(100 * int(4*pi*(etoile(istar)%r/pix_size)**2), n_ray_star_SED)
           if (istar==1) write(*,*) ""
           write(*,*) "Star is resolved, using",n_ray_star,"rays for the stellar disk"
        endif
     enddo
  endif

  do istar=1, n_etoiles
     ! if (etoile(istar)%icell == 0) cycle ! star is not in the grid ! We don't need to skip those stars anymore

     ! Compute optical depth screen in front of the star at limited resolution, e.g. 10x10
     delta = etoile(istar)%r / nx_screen
     norm_screen2 = 1./delta**2

     dx_screen(:) = delta * dx_map(:)/norm2(dx_map(:))
     dy_screen(:) = delta * dy_map(:)/norm2(dy_map(:))

     ! Todo : this could be made parallel, but it is very fast so far
     id = 1
     do j=-nx_screen, nx_screen
        do i=-nx_screen, nx_screen
           x = etoile(istar)%x + dx_screen(1) * i +  dy_screen(1) * j
           y = etoile(istar)%y + dx_screen(2) * i +  dy_screen(2) * j
           z = etoile(istar)%z + dx_screen(3) * i +  dy_screen(3) * j

           icell = etoile(istar)%icell
           call indice_cellule(x,y,z, icell)

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

     norme = 0.0_dp
     tau_avg = 0.0_dp

     ! Etoile ponctuelle
     !  x0=0.0_dp ;  y0= 0.0_dp ; z0= 0.0_dp
     !  Stokes = 0.0_dp
     !  call optical_length_tot(1,lambda,Stokes,i,j,x0,y0,z0,u,v,w,tau,lmin,lmax)
     !  Flux_etoile =  exp(-tau)
     !  write(*,*)  "F0", Flux_etoile

     ! Etoile non ponctuelle
     !$omp parallel &
     !$omp default(none) &
     !$omp shared(stream,istar,n_ray_star,llimb_darkening,limb_darkening,mu_limb_darkening,lsepar_pola) &
     !$omp shared(pola_limb_darkening,lambda,u,v,w,tab_RT_az,lsed,etoile,l3D,RT_sed_method,lpola,lmono0) &
     !$omp shared(x_center,y_center,taille_pix,dx_map,dy_map,nb_proc,map_1star,Q_1star,U_1star,lresolved) &
     !$omp shared(tau_screen,dx_screen,dy_screen,norm_screen2) &
     !$omp private(id,i,j,iray,rand,rand2,x,y,z,srw02,argmt,cos_thet,LimbDarkening,Stokes,fx,fy,offset_x,offset_y,vec) &
     !$omp private(Pola_LimbDarkening,icell,tau,lmin,lmax,in_map,P,phi,is_in_image) &
     !$omp reduction(+:norme,tau_avg)
     in_map = .true. ! for SED
     LimbDarkening = 1.0

     id = 1 ! Pour code sequentiel
     !$ id = omp_get_thread_num() + 1

     !$omp do schedule(static,n_ray_star(istar)/nb_proc)
     do iray=1,n_ray_star(istar)
        ! Position aleatoire sur la disque stellaire
        rand  = sprng(stream(id))
        rand2 = sprng(stream(id))

        ! Position de depart aleatoire sur une sphere de rayon 1
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

        ! Position de depart aleatoire sur une sphere de rayon r_etoile
        vec = (/x,y,z/) * etoile(istar)%r ! offset vector from center of star
        x = etoile(istar)%x + vec(1)
        y = etoile(istar)%y + vec(2)
        z = etoile(istar)%z + vec(3)

        ! Compute exact optical depth for each point
        !icell = etoile(istar)%icell
        !Stokes = 0.0_dp
        !call optical_length_tot(id,lambda,Stokes,icell,x,y,z,u,v,w,tau,lmin,lmax)

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
        norme = norme + cos_thet * LimbDarkening
     enddo ! iray
     !$omp end do
     !$omp end parallel

     ! Normalizing map and Adding all the stars
     facteur2 =  (facteur * prob_E_star(lambda,istar)) / norme

     do id=1, nb_proc
        stars_map(:,:,1) = stars_map(:,:,1) + map_1star(:,:,id) * facteur2
     enddo

     if (lpola) then
        ! Normalizing maps and adding all the stars
        stars_map(:,:,2) = stars_map(:,:,2) + sum(Q_1star(:,:,:),dim=3) * facteur2
        stars_map(:,:,3) = stars_map(:,:,3) + sum(U_1star(:,:,:),dim=3) * facteur2
     endif

     if (lmono0) then
        tau_avg = tau_avg/n_ray_star(istar)
        write(*,fmt='(" Optical depth from star #", i2, " is ", E12.5)') istar, tau_avg
     endif


     !---  Projected position of centres of each star
     xyz(1) = etoile(istar)%x ; xyz(2) = etoile(istar)%y ; xyz(3) = etoile(istar)%z

     ! Offset from map center in arcsec
     star_position(istar,ibin,iaz,1) = - dot_product(xyz, dx_map) * factor_pix ! RA negative axis
     star_position(istar,ibin,iaz,2) = dot_product(xyz, dy_map) * factor_pix

     ! Radial velocities
     star_vr(istar,ibin,iaz) = etoile(istar)%vx * u + etoile(istar)%vy * v + etoile(istar)%vz * w

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
  ! Calcule l'intensite d'un pixel carre de taille, position et orientation arbitaires
  ! par une methode de Ray-tracing
  ! (u,v,w) pointe vers l'observateur
  ! TODO : Integration par methode de Romberg pour determiner le nbre de sous-pixel
  ! necessaire
  ! Unite : W.m-2 : nu.F_nu
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

  ! TODO : il y a un truc bizarre dans cette routine !!!

  ! Ray tracing : on se propage dans l'autre sens
  u0 = -u ; v0 = -v ; w0 = -w

  ! le nbre de subpixel en x est 2^(iter-1)
  subpixels = 1
  iter = 1

  infinie : do ! boucle jusqu'a convergence

     npix2 =  real(subpixels)**2
     Stokes_old(:) = Stokes(:)
     Stokes(:) = 0.0_dp

     ! Vecteurs definissant les sous-pixels
     sdx(:) = dx(:) / real(subpixels,kind=dp)
     sdy(:) = dy(:) / real(subpixels,kind=dp)

     ! L'obs est en dehors de la grille
     ri = 2*n_rad ; zj=1 ; phik=1

     ! Boucle sur les sous-pixels qui calcule l'intensite au centre
     ! de chaque sous pixel
     do i = 1,subpixels
        do j = 1,subpixels
           ! Centre du sous-pixel
           x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
           y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
           z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)

           ! On se met au bord de la grille : propagation a l'envers
           call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)  !BUG

           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
              ! Flux recu dans le pixel
              Stokes(:) = Stokes(:) + integ_ray_dust(lambda,icell,x0,y0,z0,u0,v0,w0)
           endif
        enddo !j
     enddo !i
     Stokes(:) = Stokes(:) / npix2

     if (iter < n_iter_min) then
        ! On itere par defaut
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! On arrete pour pas tourner dans le vide
        ! write(*,*) "Warning : converging pb in ray-tracing"
        ! write(*,*) " Pixel", ipix, jpix
        exit infinie
     else
        ! On fait le test sur a difference
        if (abs(Stokes(1) - Stokes_old(1)) > precision * Stokes_old(1)) then
           ! On est pas converge
           subpixels = subpixels * 2
        else
           ! On est converge
           exit infinie
        endif
     endif ! iter

     iter = iter + 1

  enddo infinie

  ! Prise en compte de la surface du pixel (en sr)
  Stokes = Stokes * (pixelsize / (distance*pc_to_AU) )**2

  if (lsed) then
     ! Sommation sur les pixels implicite
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

  ! Direction de visee pour le ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition des vecteurs de base du plan image dans le repere universel

  ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Vecteur x image avec PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0


  ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
  taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
  dx(:) = x_plan_image * taille_pix
  dy(:) = y_plan_image * taille_pix

  ! Coin en bas gauche de l'image
  Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

  ! Boucle sur les pixels de l'image
  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,id,Stokes,icell,lintersect,x0,y0,z0,u0,v0,w0) &
  !$omp private(flag_star,flag_direct_star,ltot,flag_sortie,lpacket_alive) &
  !$omp shared(tau,Icorner,lambda,P_lambda,pixelcenter,dx,dy,u,v,w) &
  !$omp shared(taille_pix,npix_x,npix_y,ibin,iaz,tau_surface_map,move_to_grid)
  id = 1 ! pour code sequentiel

  !$omp do schedule(dynamic,1)
  do i = 1, npix_x
     !$ id = omp_get_thread_num() + 1
     do j = 1,npix_y
        ! Coin en bas gauche du pixel
        pixelcenter(:,id) = Icorner(:) + (i-0.5_dp) * dx(:) + (j-0.5_dp) * dy(:)

        x0 = pixelcenter(1,id)
        y0 = pixelcenter(2,id)
        z0 = pixelcenter(3,id)

        ! Ray tracing : on se propage dans l'autre sens
        u0 = -u ; v0 = -v ; w0 = -w

        ! On se met au bord de la grille : propagation a l'envers
        call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

        if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
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

  ! Direction de visee pour le ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition des vecteurs de base du plan image dans le repere universel

  ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Vecteur x image avec PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0


  ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
  taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
  dx(:) = x_plan_image * taille_pix
  dy(:) = y_plan_image * taille_pix

  ! Coin en bas gauche de l'image
  Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

  ! Boucle sur les pixels de l'image
  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,id,Stokes,icell,lintersect,x0,y0,z0,u0,v0,w0) &
  !$omp private(flag_star,flag_direct_star,flag_sortie,lpacket_alive,tau,lmin,lmax) &
  !$omp shared(Icorner,lambda,P_lambda,pixelcenter,dx,dy,u,v,w) &
  !$omp shared(taille_pix,npix_x,npix_y,ibin,iaz,tau_map,move_to_grid)
  id = 1 ! pour code sequentiel

  !$omp do schedule(dynamic,1)
  do i = 1, npix_x
     !$ id = omp_get_thread_num() + 1
     do j = 1,npix_y
        ! Coin en bas gauche du pixel
        pixelcenter(:,id) = Icorner(:) + (i-0.5_dp) * dx(:) + (j-0.5_dp) * dy(:)

        x0 = pixelcenter(1,id)
        y0 = pixelcenter(2,id)
        z0 = pixelcenter(3,id)

        ! Ray tracing : on se propage dans l'autre sens
        u0 = -u ; v0 = -v ; w0 = -w

        ! On se met au bord de la grille : propagation a l'envers
        call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

        if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
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
