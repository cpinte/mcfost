module mem

  use parametres
  use grains
  use em_th
  use opacity
  use resultats
  use prop_star
  use naleat
  use molecular_emission
  use ray_tracing
  use utils

  implicit none

  contains

subroutine alloc_dynamique()
  ! Alloue les tableaux dynamiquement en fonction de l'organisation
  ! des calculs, de la presence de stratification ...
  ! Permet d'optimiser la taille mémoire
  ! C. Pinte
  ! 12/05/05

  integer ::  alloc_status
  real :: mem_size

  nrz = n_rad * nz
  if (l3D) then
     n_cells = 2*nrz*n_az
  else
     n_cells = nrz
  endif

  if (n_cells < 1e6) then
     write(*,*) "Using", n_cells, "cells"
  else
     write(*,*) "Using", real(n_cells)/1e6, "million cells"
  endif

  if (lvariable_dust) then
     p_n_rad=n_rad ; p_nz = nz ; p_n_cells = n_cells
  else
     p_n_rad=1 ;  p_nz=1 ; p_n_cells = 1
  endif

  lonly_LTE = .false.
  lonly_nLTE = .false.
  if (lRE_LTE .and. .not.lRE_nLTE .and. .not. lnRE) lonly_LTE = .true.
  if (lRE_nLTE .and. .not.lRE_LTE .and. .not. lnRE) lonly_nLTE = .true.

  if (l3D) then
     j_start = -nz
     if (lvariable_dust) then
        p_n_az = n_az
     else
        p_n_az = 1
     endif
  else
     j_start = 1
     p_n_az = 1
  endif

  if ((p_nz /= 1).and.l3D) then
     pj_start = -nz
  else
     pj_start = 1
  endif

  ! parametrage methode de diffusion
  ! 1 : per dust grain
  ! 2 : per cell
  if (scattering_method == 0) then
     if (.not.lmono) then
        mem_size = (1.0*p_n_cells) * (nang_scatt+1) * n_lambda * 4 / 1024**3
        if (mem_size > max_mem) then
           scattering_method = 1
        else
           scattering_method = 2
        endif
     else
        if (lscatt_ray_tracing) then
           scattering_method = 2
        else
           scattering_method = 2
        endif
     endif
  endif

  lMueller_pos_multi = .false.
  if (lmono) then
     p_n_lambda_pos = 1
  else
     if (scattering_method==1) then
        p_n_lambda_pos = 1
     else
        p_n_lambda_pos = n_lambda
        lMueller_pos_multi = .true.
     endif
  endif

  write(*,fmt='(" Using scattering method ",i1)') scattering_method
  lscattering_method1 = (scattering_method==1)

  allocate(stream(nb_proc), gauss_random_saved(nb_proc), lgauss_random_saved(nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error random number stream'
     stop
  endif
  gauss_random_saved = 0.0_db
  lgauss_random_saved = .false.

  allocate(n_phot_envoyes(n_lambda,nb_proc),  stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_phot_envoyes'
     stop
  endif
  n_phot_envoyes = 0.0

  ! **************************************************
  ! Tableaux relatifs a la grille
  ! **************************************************
  allocate(r_lim(0:n_rad), r_lim_2(0:n_rad), r_lim_3(0:n_rad), &
  delta_z(n_rad), dr2_grid(n_rad), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_lim'
     stop
  endif
  r_lim = 0.0 ; r_lim_2=0.0; r_lim_3=0.0 ; delta_z=0.0 ; dr2_grid=0.0

  allocate(r_grid(n_cells), z_grid(n_cells), phi_grid(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_lim'
     stop
  endif
  r_grid=0.0; z_grid=0.0 ; phi_grid = 0.0

  allocate(tab_region(n_rad), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_region'
     stop
  endif
  tab_region = 0

  allocate(z_lim(n_rad,nz+2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error z_lim'
     stop
  endif
  z_lim = 0.0

  allocate(w_lim(0:nz),  theta_lim(0:nz),tan_theta_lim(0:nz),tan_phi_lim(n_az), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tan_phi_lim'
     stop
  endif
  w_lim = 0.0
  theta_lim=0.0
  tan_theta_lim = 0.0
  tan_phi_lim = 0.0

  allocate(zmax(n_rad),volume(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error zmax, volume'
     stop
  endif
  zmax = 0.0 ; volume=0.0

  allocate(masse(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error masse'
     stop
  endif
  masse = 0.0

  allocate(densite_pouss(n_grains_tot,n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error densite_pouss'
     stop
  endif
  densite_pouss = 0.0

  if (lSigma_file) then
     allocate(Surface_density(n_rad), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Sigma'
        stop
     endif
     Surface_density = 0.0_db
  endif

  allocate(l_dark_zone(n_cells), ri_in_dark_zone(n_az), ri_out_dark_zone(n_az),&
          zj_sup_dark_zone(n_rad,n_az), zj_inf_dark_zone(n_rad,n_az), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error l_dark_zone'
     stop
  endif
  l_is_dark_zone = .false.
  l_dark_zone = .false.
  ri_in_dark_zone=0
  ri_out_dark_zone=0
  zj_sup_dark_zone=0
  if (l3D) zj_inf_dark_zone=0

  if (l3D) then
     allocate(r_in_opacite(-nz-1:nz+1,n_az), r_in_opacite2(-nz-1:nz+1,n_az), stat=alloc_status)
  else
     allocate(r_in_opacite(nz+1,1), r_in_opacite2(nz+1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_in_opacite'
     stop
  endif
  r_in_opacite=0.0 ; r_in_opacite2=0.0

  ! **************************************************
  ! Tableaux relatifs aux grains
  ! **************************************************
  allocate(nbre_grains(n_grains_tot), r_grain(n_grains_tot),  r_grain_min(n_grains_tot), r_grain_max(n_grains_tot), &
       S_grain(n_grains_tot), M_grain(n_grains_tot), r_core(n_grains_tot), grain(n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_grain'
     stop
  endif
  nbre_grains = 0.0   ; r_core=0.0
  r_grain=0.0 ; r_grain_min=0.0 ; r_grain_max=0.0 ; S_grain=0.0 ; M_grain=0.0

  allocate(tab_albedo(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo'
     stop
  endif
  tab_albedo = 0

  allocate(C_ext(n_grains_tot,n_lambda), C_sca(n_grains_tot,n_lambda), &
       C_abs(n_grains_tot,n_lambda), C_abs_norm(n_grains_tot,n_lambda), stat=alloc_status) !  C_geo(n_grains_tot)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error q_ext'
     stop
  endif
  C_ext = 0 ; C_sca = 0 ; C_abs = 0 ; C_abs_norm =0


  allocate(tab_g(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_g'
     stop
  endif
  tab_g = 0

  ! **************************************************
  ! Tableaux relatifs aux prop en fct de lambda
  ! **************************************************
  allocate(E_stars(n_lambda), E_disk(n_lambda), E_ISM(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error E_stars'
     stop
  endif
  E_stars = 0.0
  E_disk = 0.0
  E_ISM = 0.0

  allocate(frac_E_stars(n_lambda), frac_E_disk(n_lambda), E_totale(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error frac_E_stars'
     stop
  endif
  frac_E_stars = 0.0 ; frac_E_disk = 0.0 ; E_totale = 0.0

  allocate(spectre_etoiles_cumul(0:n_lambda), spectre_etoiles(n_lambda), spectre_emission_cumul(0:n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre_etoile'
     stop
  endif
  spectre_etoiles_cumul = 0.0
  spectre_etoiles = 0.0
  spectre_emission_cumul = 0.0

  allocate(tab_lambda(n_lambda), tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_lambda'
     stop
  endif
  tab_lambda=0.0 ;  tab_lambda_inf=0.0 ; tab_lambda_sup=0.0 ; tab_delta_lambda=0.0

  allocate(tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_amu'
     stop
  endif
  tab_amu1=0.0; tab_amu2=0.0;

  allocate(tab_amu1_coating(n_lambda, n_pop), tab_amu2_coating(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_amu_coating'
     stop
  endif
  tab_amu1_coating=0.0; tab_amu2_coating=0.0;


  ! Tableaux relatifs aux prop optiques des cellules
  allocate(kappa(n_cells,n_lambda), kappa_sca(n_cells,n_lambda), kappa_abs_LTE(n_cells,n_lambda), stat=alloc_status)

  if (.not.lonly_LTE .or. .not.lonly_nLTE) then
     allocate(proba_abs_RE_LTE(n_cells,n_lambda),  stat=alloc_status)
     if (lRE_nLTE)  allocate(kappa_abs_nLTE(n_cells,n_lambda), stat=alloc_status)
     if (lRE_nLTE.or.lnRE) allocate(proba_abs_RE_LTE_p_nLTE(n_cells,n_lambda), stat=alloc_status)
     if (lnRE) allocate(proba_abs_RE(n_cells,n_lambda), kappa_abs_RE(n_cells,n_lambda), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error kappa'
        stop
     endif
     kappa=0.0 ; kappa_abs_LTE=0.0
     proba_abs_RE_LTE=0.0
     if (lRE_nLTE) kappa_abs_nLTE=0.0
     if (lRE_nLTE.or.lnRE) proba_abs_RE_LTE_p_nLTE=0.0
     if (lnRE)  then
        proba_abs_RE=0.0
        kappa_abs_RE=0.0
     endif
  endif

  ! todo : could be p_n_cells ...
  allocate(tab_albedo_pos(n_cells,n_lambda), tab_g_pos(n_cells,n_lambda),stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos, tab_g_pos'
     stop
  endif
  tab_albedo_pos = 0 ; tab_g_pos = 0.0

  ! **************************************************
  ! Tableaux relatifs aux prop optiques des cellules ou des grains
  ! **************************************************
  if (scattering_method == 2) then ! prop par cellule
     if (lsepar_pola) then
        mem_size = (5. * nang_scatt) * p_n_cells * p_n_lambda_pos * 4. / 1024.**3
     else
        mem_size = (2. * nang_scatt) * p_n_cells * p_n_lambda_pos * 4. / 1024.**3
     endif
     if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering matrices"

     allocate(tab_s11_pos(0:nang_scatt, p_n_cells, p_n_lambda_pos), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_s11_pos'
        stop
     endif
     tab_s11_pos = 0

     allocate(prob_s11_pos(0:nang_scatt, p_n_cells, p_n_lambda_pos), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_s11_pos'
        stop
     endif
     prob_s11_pos = 0

     if (lsepar_pola) then
        allocate(tab_s12_o_s11_pos(0:nang_scatt, p_n_cells, p_n_lambda_pos), &
             tab_s33_o_s11_pos(0:nang_scatt, p_n_cells, p_n_lambda_pos), &
             tab_s34_o_s11_pos(0:nang_scatt, p_n_cells, p_n_lambda_pos), &
             stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s12_o_s11_pos'
           stop
        endif
        tab_s12_o_s11_pos = 0
        tab_s33_o_s11_pos = 0
        tab_s34_o_s11_pos = 0
     endif
  else ! scattering method==1 --> prop par grains
     mem_size = (1.0 * n_grains_tot) * p_n_cells * n_lambda * 4. / 1024.**3
     if (mem_size < max_mem) then
        low_mem_scattering = .false.
        if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"

        allocate(ksca_CDF(0:n_grains_tot,p_n_cells,n_lambda), stat=alloc_status)

        if (alloc_status > 0) then
           write(*,*) 'Allocation error ksca_CDF'
           stop
        endif
        ksca_CDF = 0
     else ! Array is to big, we will recompute ksca_CDF on the fly
        low_mem_scattering = .true.
        write(*,*) "Using low memory mode for scattering properties"
     endif
  endif ! method

  ! **************************************************
  ! tableaux relatifs aux prop optiques des grains
  ! **************************************************
  allocate(tab_s11(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11'
     stop
  endif
  tab_s11 = 0

  allocate(tab_s12(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s12'
     stop
  endif
  tab_s12 = 0

  allocate(tab_s33(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s33'
     stop
  endif
  tab_s33 = 0

  allocate(tab_s34(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s34'
     stop
  endif
  tab_s34 = 0

  if (laggregate) then
     allocate(tab_mueller(4,4,0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_mueller'
        stop
     endif
     tab_mueller = 0
  endif

  allocate(prob_s11(n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_s11'
     stop
  endif
  prob_s11 = 0

  ! **************************************************
  ! Tableaux relatifs aux prop d'emission des cellules
  ! **************************************************
  allocate(l_emission_pah(0:n_rad+1,0:nz+1), stat=alloc_status) ! OUTDATED
  if (alloc_status > 0) then
     write(*,*) 'Allocation error l_emission_pah'
     stop
  endif
  l_emission_pah = .false.


  allocate(prob_E_cell(0:n_cells,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_cell'
     stop
  endif
  prob_E_cell = 0.0

  if (lweight_emission) then
     allocate(weight_proba_emission(n_cells), correct_E_emission(n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_E_cell'
        stop
     endif
     weight_proba_emission = 1.0 ;
     correct_E_emission = 1.0 ;
  endif


  if (lorigine) then
     allocate(disk_origin(n_lambda, n_rad, nz, nb_proc), star_origin(n_lambda, nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error disk_origin'
        stop
     endif
     disk_origin = 0.0
     star_origin = 0.0
  endif

  ! **************************************************
  ! Tableaux de temperature
  ! **************************************************
  allocate(Temperature(n_cells), Temperature_old(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error Temperature'
     stop
  endif
  Temperature = 0.0 ; Temperature_old=0.0


  if (lRE_nLTE) then
     allocate(Temperature_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_cells),stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Temperature_1grain'
        stop
     endif
     Temperature_1grain = 0.0
  endif

  if (lnRE) then
     if ( (.not.ltemp).and.(lsed.or.lmono0.or.lProDiMo.or.lProDiMo2mcfost) ) then ! si ltemp --> tableau alloue ci-dessous
        allocate(tab_Temp(n_T), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_Temp'
           stop
        endif
        tab_Temp = 0.0
     endif

     allocate(Proba_Temperature(n_T,grain_nRE_start:grain_nRE_end,n_cells), &
          Temperature_1grain_nRE(grain_nRE_start:grain_nRE_end,n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Proba_Temperature'
        stop
     endif
     Proba_Temperature=0.0
     Temperature_1grain_nRE=0.0

     allocate(l_RE(grain_nRE_start:grain_nRE_end,n_cells), lchange_nRE(grain_nRE_start:grain_nRE_end,n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error l_RE'
        stop
     endif
     l_RE=.false. ; lchange_nRE = .false.
  endif

  ! **************************************************
  ! Tableaux relatifs au *calcul* de la temperature
  ! **************************************************
  if (lTemp) then
     allocate(tab_Temp(n_T), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_Temp'
        stop
     endif
     tab_Temp = 0.0

     allocate(log_frac_E_em(n_T,n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error log_frac_E_em'
        stop
     endif
     log_frac_E_em = 0.0

     allocate(DensE(n_rad,0:nz,n_az), DensE_m1(n_rad,0:nz,n_az), Dcoeff(n_rad,0:nz,n_az), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error kappa_abs_1grain'
        stop
     endif
     DensE = 0.0 ; DensE_m1 = 0.0

     if (lRE_LTE) then
        allocate(xKJ_abs(n_cells,nb_proc), E0(n_cells), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error xKJ_abs'
           stop
        endif
        xKJ_abs = 0.0 ; E0 = 0.0
     endif

     lxJ_abs = loutput_J .or. lRE_nLTE .or. lnRE .or. loutput_UV_field
     if (lxJ_abs) then
        allocate(xJ_abs(n_cells,n_lambda,nb_proc), J0(n_cells,n_lambda), stat=alloc_status) ! BIG array
        if (alloc_status > 0) then
           write(*,*) 'Allocation error xJ_abs'
           stop
        endif
        xJ_abs=0.0 ; J0 = 0.0
     endif

     allocate(xT_ech(n_cells,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error xT_ech'
        stop
     endif
     xT_ech = 2

     if (lreemission_stats) then
        allocate(nbre_reemission(n_cells,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error nbre_reemission'
           stop
        endif
        nbre_reemission = 0.0
     endif

     mem_size = (1.0 * p_n_cells) * n_T * n_lambda * 4. / 1024.**3
     if (mem_size < max_mem) then
        low_mem_th_emission = .false.
        if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for temperature calculation"
        allocate(kdB_dT_CDF(n_lambda,n_T,p_n_cells), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error kdB_dT_CDF'
           stop
        endif
        kdB_dT_CDF = 0
     else
        low_mem_th_emission = .true.
        write(*,*) "Using low memory mode for thermal emission"
        allocate(kdB_dT_1grain_LTE_CDF(n_lambda,grain_RE_LTE_start:grain_RE_LTE_end,n_T), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error kdB_dT_1grain_LTE_CDF'
           stop
        endif
        kdB_dT_1grain_LTE_CDF = 0
     endif

     if (lRE_nLTE) then

        mem_size = (1.0 * (grain_RE_nLTE_end-grain_RE_nLTE_start+2)) * p_n_cells * n_lambda * 4. / 1024.**3
        if (mem_size < max_mem) then
           low_mem_th_emission_nLTE = .false.
           if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"
           allocate(kabs_nLTE_CDF(grain_RE_nLTE_start-1:grain_RE_nLTE_end,n_cells,n_lambda),stat=alloc_status)
           if (alloc_status > 0) then
              write(*,*) 'Allocation error kabs_nLTE_CDF'
              stop
           endif
           kabs_nLTE_CDF = 0.0
        else
           low_mem_th_emission_nLTE = .true.
           write(*,*) "Using low memory mode for nLTE thermal emission"
        endif

        allocate(kdB_dT_1grain_nLTE_CDF(n_lambda,grain_RE_nLTE_start:grain_RE_nLTE_end,n_T),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error kdB_dT_1grain_nLTE_CDF'
           stop
        endif
        kdB_dT_1grain_nLTE_CDF=0.0

        allocate(log_frac_E_em_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_T),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error log_frac_E_em_1grain'
           stop
        endif
        log_frac_E_em_1grain=0.0

        allocate(xT_ech_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_cells,nb_proc),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error xT_ech_1grain'
           stop
        endif
        xT_ech_1grain = 2
     endif


     if (lnRE) then
        allocate(frac_E_em_1grain_nRE(grain_nRE_start:grain_nRE_end,n_T),&
             log_frac_E_em_1grain_nRE(grain_nRE_start:grain_nRE_end,n_T), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error frac_E_em_1grain'
           stop
        endif
        frac_E_em_1grain_nRE=0.0
        log_frac_E_em_1grain_nRE=0.0

        allocate(Temperature_1grain_nRE_old(grain_nRE_start:grain_nRE_end,n_cells), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error Temperature_1grain_nRE_old'
           stop
        endif
        Temperature_1grain_nRE_old =0.0

        allocate(Tpeak_old(grain_nRE_start:grain_nRE_end,n_cells), &
             maxP_old(grain_nRE_start:grain_nRE_end,n_cells), &
             stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error Tpeak'
           stop
        endif
        Tpeak_old=0
        maxP_old=0.

        allocate(xT_ech_1grain_nRE(grain_nRE_start:grain_nRE_end,n_cells,nb_proc),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error xT_ech_1grain_nRE'
           stop
        endif
        xT_ech_1grain_nRE = 2

        allocate(kdB_dT_1grain_nRE_CDF(n_lambda,grain_nRE_start:grain_nRE_end,n_T),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error kdB_dT_1grain_nRE_CDF'
           stop
        endif
        kdB_dT_1grain_nRE_CDF=0.0

        if (lRE_nlTE) then
           allocate(Temperature_1grain_old(grain_RE_nLTE_start:grain_RE_nLTE_end,n_cells),stat=alloc_status)
           if (alloc_status > 0) then
              write(*,*) 'Allocation error Temperature_1grain_old'
              stop
           endif
           Temperature_old=0.
        endif
     endif

  endif ! lTemp

  if (lnRE) then
     allocate(Emissivite_nRE_old(n_cells,n_lambda), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Emissivite_nRE_old'
        stop
     endif
     Emissivite_nRE_old = 0.0
  endif

  ! **************************************************
  ! Tableaux relatifs aux SEDs
  ! **************************************************
  if (lTemp.or.lsed) then
     allocate(sed(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed'
        stop
     endif
     sed = 0.0

     allocate(sed_q(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_q'
        stop
     endif
     sed_q = 0.0

     allocate(sed_u(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_u'
        stop
     endif
     sed_u = 0.0

     allocate(sed_v(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_v'
        stop
     endif
     sed_v = 0.0

     allocate(sed_star(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_star'
        stop
     endif
     sed_star = 0.0

     allocate(sed_star_scat(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_star_scat'
        stop
     endif
     sed_star_scat = 0.0

     allocate(sed_disk(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_disk'
        stop
     endif
     sed_disk = 0.0

     allocate(sed_disk_scat(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_disk_scat'
        stop
     endif
     sed_disk_scat = 0.0

     allocate(n_phot_sed(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error n_phot_sed'
        stop
     endif
     n_phot_sed = 0.0

     allocate(sed1_io(n_lambda,N_thet,N_phi),sed2_io(n_lambda,N_thet,N_phi,9),wave2_io(n_lambda,2), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_io'
        stop
     endif
     sed1_io=0.0
     sed2_io=0.0
     wave2_io=0.0
  endif ! ltemp.or.lSED


  ! **************************************************
  ! Tableaux relatifs aux images
  ! **************************************************
  if (lmono0.and.loutput_mc) then
     allocate(STOKE_io(igridx,igridy,capt_debut:capt_fin,N_phi,N_type_flux), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error STOKE_io'
        stop
     endif
     STOKE_io = 0.0

     allocate(STOKEI(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error STOKEI'
        stop
     endif
     STOKEI = 0.0

     if (lsepar_pola) then
        allocate(STOKEQ(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEQ'
           stop
        endif
        STOKEQ = 0.0

        allocate(STOKEU(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEU'
           stop
        endif
        STOKEU=0.0

        allocate(STOKEV(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEV'
           stop
        endif
        STOKEV = 0.0
     endif

     if (lsepar_contrib) then
        allocate(STOKEI_star(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_star'
           stop
        endif
        STOKEI_star = 0.0

        allocate(STOKEI_star_scat(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_star_scat'
           stop
        endif
        STOKEI_star_scat = 0.0

        allocate(STOKEI_disk(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_disk'
           stop
        endif
        STOKEI_disk = 0.0

        allocate(STOKEI_disk_scat(n_lambda,IGRIDX,IGRIDY,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_disk_scat'
           stop
        endif
        STOKEI_disk_scat = 0.0
     endif !lsepar

  endif ! lmono0

  ! **************************************************
  ! Tableaux relatifs a l'emission moleculaire
  ! **************************************************
  if (lemission_mol) then
     allocate(tab_abundance(n_cells), Tcin(n_cells), lcompute_molRT(n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Tcin & tab_abundance'
        stop
     endif
     tab_abundance = 0.0
     lcompute_molRT = .true.
     Tcin=0.0

     allocate(vfield(n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error vfield'
        stop
     endif
     vfield=0.0 !; vx=0.0 ; vy=0.0

     allocate(v_turb(n_cells), v_line(n_cells), deltaVmax(n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sigma2'
        stop
     endif
     v_turb = 0.0 ; v_line = 0.0 ;   deltaVmax = 0.0

     allocate(tab_dnu_o_freq(n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_dnu'
        stop
     endif
     tab_dnu_o_freq=0.0

     allocate(norme_phiProf_m1(n_cells), sigma2_phiProf_m1(n_cells), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error norme_phiProf_m1'
        stop
     endif
     norme_phiProf_m1 = 0.0 ; sigma2_phiProf_m1 = 0.0

  endif ! lemission_mol

  allocate(densite_gaz(n_cells), densite_gaz_midplane(n_rad), masse_gaz(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error densite_gaz'
     stop
  endif
  densite_gaz = 0.0 ; densite_gaz_midplane = 0.0 ; masse_gaz = 0.0

  return

end subroutine alloc_dynamique

!**********************************************************************

subroutine dealloc_em_th()

  deallocate(n_phot_envoyes)


  deallocate(tab_albedo,C_ext,C_sca,C_abs,C_abs_norm,tab_g) ! q_geo

  !deallocate(E_stars,E_disk,frac_E_stars,E_totale)

  deallocate(spectre_etoiles_cumul,spectre_etoiles, spectre_emission_cumul, CDF_E_star, prob_E_star, E_stars)

  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2,tab_amu1_coating,tab_amu2_coating)

  deallocate(kappa,kappa_abs_LTE,kappa_sca)
  if (allocated(proba_abs_RE_LTE)) then
     deallocate(proba_abs_RE_LTE)
     if (lRE_nLTE) deallocate(kappa_abs_nLTE)
     if (lRE_nLTE.or.lnRE) deallocate(proba_abs_RE_LTE_p_nLTE)
     if (lnRE) deallocate(proba_abs_RE,kappa_abs_RE)
  endif

  deallocate(tab_albedo_pos,tab_g_pos)

  if (scattering_method == 2) then ! prop par cellule
     deallocate(tab_s11_pos,prob_s11_pos)
     if (lsepar_pola) deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
  else ! prop par grains
     deallocate(ksca_CDF)
  endif ! method

  deallocate(tab_s11,tab_s12,tab_s33,tab_s34,prob_s11)

  deallocate(l_emission_pah) ! OUTDATED
  deallocate(prob_E_cell)
  if (lorigine) deallocate(disk_origin,star_origin)

  if (lTemp) then
     deallocate(tab_Temp)
     if (lsed_complete) then
        deallocate(log_frac_E_em)
        deallocate(DensE, DensE_m1, Dcoeff)
        if (allocated(xKJ_abs)) deallocate(xKJ_abs,E0)
        if (allocated(xJ_abs)) deallocate(xJ_abs,J0)
        if (allocated(nbre_reemission)) deallocate(nbre_reemission)
        if (allocated(kdB_dT_CDF)) deallocate(xT_ech,kdB_dT_CDF)
     endif

     if (lRE_nLTE) then
        deallocate(kabs_nLTE_CDF,kdB_dT_1grain_nLTE_CDF,log_frac_E_em_1grain)
        deallocate(xT_ech_1grain)
     endif

     if (lnRE) then
        deallocate(frac_E_em_1grain_nRE,log_frac_E_em_1grain_nRE)
        deallocate(Temperature_1grain_nRE_old,kdB_dT_1grain_nRE_CDF,xT_ech_1grain_nRE)
        deallocate(Emissivite_nRE_old)
        deallocate(Tpeak_old)
        if (lRE_nlTE) deallocate(Temperature_1grain_old)
     endif
  endif ! lTemp

  if (lTemp.or.lsed) then
     deallocate(sed,sed_q,sed_u,sed_v)
     deallocate(sed_star,sed_star_scat,sed_disk,sed_disk_scat,n_phot_sed)
     deallocate(sed1_io,sed2_io,wave2_io)
  endif ! ltemp.or.lSED


  return

end subroutine dealloc_em_th

!******************************************************************************

subroutine realloc_dust_mol()

  integer :: alloc_status, mem_size

  allocate(tab_lambda(n_lambda), tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), &
       tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), &
       tab_amu1_coating(n_lambda, n_pop), tab_amu2_coating(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_lambda (realloc)'
     stop
  endif
  tab_lambda=0.0 ; tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda= 0.0 ;
  tab_amu1=0.0 ; tab_amu2=0.0 ; tab_amu1_coating=0.0 ; tab_amu2_coating=0.0

  allocate(tab_albedo(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo (realloc)'
     stop
  endif
  tab_albedo = 0

  allocate(C_ext(n_grains_tot,n_lambda), C_sca(n_grains_tot,n_lambda), &
       C_abs(n_grains_tot,n_lambda),  C_abs_norm(n_grains_tot,n_lambda), tab_g(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error C_ext (realloc)'
     stop
  endif
  C_ext = 0 ; C_sca = 0 ; C_abs = 0 ; C_abs_norm = 0 ; tab_g = 0


  allocate(prob_s11(n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_s11 (realloc)'
     stop
  endif
  prob_s11 = 0

  allocate(tab_s11(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11 (realloc)'
     stop
  endif
  tab_s11 = 0

  allocate(tab_s12(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s12 (realloc)'
     stop
  endif
  tab_s12 = 0

  allocate(tab_s33(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s33 (realloc)'
     stop
  endif
  tab_s33 = 0

  allocate(tab_s34(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s34 (realloc)'
     stop
  endif
  tab_s34 = 0

  ! TODO : cette partie prend bcp de memoire
  low_mem_scattering = .false.
  allocate(ksca_CDF(0:n_grains_tot,p_n_cells,n_lambda), stat=alloc_status)
  mem_size = n_grains_tot * p_n_cells * n_lambda * 4. / 1024.**3
  if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"
  if (alloc_status > 0) then
     write(*,*) 'Allocation error ksca_CDF (realloc)'
     stop
  endif
  ksca_CDF = 0

  ! Tableaux relatifs aux prop optiques des cellules
  allocate(kappa(n_cells,n_lambda),kappa_abs_LTE(n_cells,n_lambda), kappa_sca(n_cells,n_lambda), &
       emissivite_dust(n_cells,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error emissivite_dust (realloc)'
     stop
  endif
  kappa = 0.0 ; kappa_abs_LTE = 0.0 ; kappa_sca = 0.0 ; emissivite_dust = 0.0

  if (.not.lonly_LTE .or. .not.lonly_nLTE) then
     allocate(proba_abs_RE_LTE(n_cells,n_lambda), stat=alloc_status)
     if (lRE_nLTE) allocate(kappa_abs_nLTE(n_cells,n_lambda), stat=alloc_status)
     if (lRE_nLTE.or.lnRE) allocate(proba_abs_RE_LTE_p_nLTE(n_cells,n_lambda), stat=alloc_status)
     if (lnRE) allocate(proba_abs_RE(n_cells,n_lambda), kappa_abs_RE(n_cells,n_lambda), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error kappa (realloc)'
        stop
     endif
     proba_abs_RE_LTE=0.0
     if (lRE_nLTE) kappa_abs_nLTE = 0.0
     if (lRE_nLTE.or.lnRE) proba_abs_RE_LTE_p_nLTE=0.0
     if (lnRE)  then
        proba_abs_RE=0.0
        kappa_abs_RE=0.0
     endif
  endif

  ! todo : could be p_n_cells
  allocate(tab_albedo_pos(n_cells,n_lambda), tab_g_pos(n_cells,n_lambda),stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos, tab_g_pos (realloc)'
     stop
  endif
  tab_albedo_pos = 0
  tab_g_pos = 0.0

  allocate(spectre_etoiles_cumul(0:n_lambda),spectre_etoiles(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre_etoile'
     stop
  endif
  spectre_etoiles_cumul = 0.0
  spectre_etoiles = 0.0

  allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), E_stars(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_star'
     stop
  endif
  CDF_E_star = 0.0
  prob_E_star = 0.0
  E_stars = 0.0


  return

end subroutine realloc_dust_mol

!******************************************************************************

subroutine clean_mem_dust_mol()

  ! Ne reste que tab_lambda, tab_delta_lambda, tab_lambda_inf, tab_lambda_sup, kappa, kappa_sca, emissivite_dust
  ! et spectre_etoiles, spectre_etoiles_cumul
  deallocate(tab_amu1, tab_amu2,tab_amu1_coating, tab_amu2_coating)
  deallocate(tab_albedo)
  deallocate(C_ext, C_sca, C_abs, C_abs_norm, tab_g)
  deallocate(prob_s11,tab_s11,tab_s12,tab_s33,tab_s34,ksca_CDF)
  deallocate(kappa_abs_LTE)
  if (allocated(proba_abs_RE_LTE)) then
     deallocate(proba_abs_RE_LTE)
     if (lRE_nLTE) deallocate(kappa_abs_nLTE)
     if (lRE_nLTE.or.lnRE) deallocate(proba_abs_RE_LTE_p_nLTE)
     if (lnRE) deallocate(proba_abs_RE,kappa_abs_RE)
  endif
  deallocate(tab_albedo_pos, tab_g_pos)

  return

end subroutine clean_mem_dust_mol

!******************************************************************************

subroutine realloc_step2()

  integer :: alloc_status, mem_size
  integer :: p_n_lambda2_pos = 1

  p_n_lambda_pos = p_n_lambda2_pos ! just in case

  ! parametrage methode de diffusion
  if (scattering_method == 0) then
     if ((lvariable_dust).and.(.not.lmono).and.(.not.lscatt_ray_tracing)) then
        scattering_method = 1
     else
        scattering_method = 2
     endif
  endif
  write(*,fmt='(" Using scattering method ",i1)') scattering_method
  lscattering_method1 = (scattering_method==1)

  ! Liberation memoire
  if (ltemp) then
     if (lRE_LTE)  deallocate(kdB_dT_CDF, log_frac_E_em, xT_ech,xKJ_abs)
     if (lRE_nLTE) deallocate(kabs_nLTE_CDF, kdB_dT_1grain_nLTE_CDF, log_frac_E_em_1grain,xT_ech_1grain)
     if (lreemission_stats) deallocate(nbre_reemission)
     if (lnRE) deallocate(kdB_dT_1grain_nRE_CDF,frac_E_em_1grain_nRE,log_frac_E_em_1grain_nRE,xT_ech_1grain_nRE)
     if (lxJ_abs) deallocate(xJ_abs,J0)
  endif

  if (lProDiMo.or.loutput_UV_field) then
     allocate(xJ_abs(n_cells,n_lambda2,nb_proc), J0(n_cells,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error xJ_abs in realloc_step2'
        stop
     endif
     xJ_abs = 0.0
  endif

  ! Liberation memoire step1 et reallocation step 2
  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2,tab_amu1_coating,tab_amu2_coating)
  allocate(tab_lambda(n_lambda2),tab_lambda_inf(n_lambda2),tab_lambda_sup(n_lambda2),tab_delta_lambda(n_lambda2),&
       tab_amu1(n_lambda2,n_pop),tab_amu2(n_lambda2,n_pop), &
       tab_amu1_coating(n_lambda2,n_pop),tab_amu2_coating(n_lambda2,n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_lambda in init_lambda2'
     stop
  endif
  tab_lambda=0.0 ; tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda=0.0
  tab_amu1=0.0 ; tab_amu2=0.0 ;  tab_amu1_coating=0.0 ; tab_amu2_coating=0.0

  deallocate(sed)
  allocate(sed(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed'
     stop
  endif
  sed = 0.0


  deallocate(sed_q)
  allocate(sed_q(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_q'
     stop
  endif
  sed_q = 0.0

  deallocate(sed_u)
  allocate(sed_u(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_u'
     stop
  endif
  sed_u = 0.0

  deallocate(sed_v)
  allocate(sed_v(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_v'
     stop
  endif
  sed_v = 0.0

  deallocate(sed_star)
  allocate(sed_star(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_star'
     stop
  endif
  sed_star = 0.0

  deallocate(sed_star_scat)
  allocate(sed_star_scat(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_star_scat'
     stop
  endif
  sed_star_scat = 0.0

  deallocate(sed_disk)
  allocate(sed_disk(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_disk'
     stop
  endif
  sed_disk = 0.0

  deallocate(sed_disk_scat)
  allocate(sed_disk_scat(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_disk_scat'
     stop
  endif
  sed_disk_scat = 0.0

  deallocate(n_phot_sed)
  allocate(n_phot_sed(n_lambda2,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_phot_sed'
     stop
  endif
  n_phot_sed = 0.0

  deallocate(sed2_io, wave2_io)
  allocate(sed2_io(n_lambda2,N_thet,N_phi,9),wave2_io(n_lambda2,2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed2_io'
     stop
  endif
  sed2_io=0.0
  wave2_io=0.0

  deallocate(n_phot_envoyes)
  allocate(n_phot_envoyes(n_lambda2,nb_proc),  stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_phot_envoyes'
     stop
  endif
  n_phot_envoyes = 0.0

  deallocate(prob_E_cell)
  allocate(prob_E_cell(0:n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_cell'
     stop
  endif
  prob_E_cell = 0.0

  deallocate(CDF_E_star,prob_E_star)
  allocate(CDF_E_star(n_lambda2,0:n_etoiles), prob_E_star(n_lambda2,n_etoiles), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_star'
     stop
  endif
  CDF_E_star = 0.0 ; prob_E_star = 0.0

  deallocate(spectre_etoiles_cumul, spectre_etoiles, spectre_emission_cumul)
  allocate(spectre_etoiles_cumul(0:n_lambda2),spectre_etoiles(n_lambda2),spectre_emission_cumul(0:n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre_etoile'
     stop
  endif
  spectre_etoiles_cumul = 0.0
  spectre_etoiles = 0.0
  spectre_emission_cumul = 0.0

  deallocate(E_stars, E_disk,E_ISM)
  allocate(E_stars(n_lambda2), E_disk(n_lambda2), E_ISM(n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error E_stars'
     stop
  endif
  E_stars = 0.0
  E_disk = 0.0
  E_ISM = 0.0

  deallocate(frac_E_stars, frac_E_disk, E_totale)
  allocate(frac_E_stars(n_lambda2), frac_E_disk(n_lambda2), E_totale(n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error frac_E_stars'
     stop
  endif
  frac_E_stars = 0.0 ; frac_E_disk = 0.0 ; E_totale = 0.0

  deallocate(tab_albedo)
  allocate(tab_albedo(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo'
     stop
  endif
  tab_albedo = 0

  deallocate(C_ext)
  allocate(C_ext(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error C_ext'
     stop
  endif
  C_ext = 0

  deallocate(C_sca)
  allocate(C_sca(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error C_sca'
     stop
  endif
  C_sca = 0

  deallocate(C_abs,C_abs_norm)
  allocate(C_abs(n_grains_tot,n_lambda2), C_abs_norm(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error C_abs'
     stop
  endif
  C_abs = 0 ; C_abs_norm = 0

  deallocate(tab_g)
  allocate(tab_g(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_g'
     stop
  endif
  tab_g = 0

  deallocate(tab_albedo_pos,tab_g_pos)
  ! todo : could be p_n_cells
  allocate(tab_albedo_pos(n_cells,n_lambda2), tab_g_pos(n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos, tab_g_pos'
     stop
  endif
  tab_albedo_pos = 0
  tab_g_pos = 0

  deallocate(tab_s11)
  allocate(tab_s11(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11'
     stop
  endif
  tab_s11 = 0

  deallocate(tab_s12)
  allocate(tab_s12(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s12'
     stop
  endif
  tab_s12 = 0

  deallocate(tab_s33)
  allocate(tab_s33(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s33'
     stop
  endif
  tab_s33 = 0

  deallocate(tab_s34)
  allocate(tab_s34(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s34'
     stop
  endif
  tab_s34 = 0

  if (laggregate) then
     deallocate(tab_mueller)
     allocate(tab_mueller(4,4,0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_mueller'
        stop
     endif
     tab_mueller = 0
  endif

  deallocate(prob_s11)
  allocate(prob_s11(n_lambda2,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_s11'
     stop
  endif
  prob_s11 = 0

  if (scattering_method == 2) then
     deallocate(tab_s11_pos)
     allocate(tab_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_s11_pos'
        stop
     endif
     tab_s11_pos = 0

     if (lsepar_pola) then
        deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
        allocate(tab_s12_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
             tab_s33_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
             tab_s34_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
             stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s12_o_s11_pos'
           stop
        endif
        tab_s12_o_s11_pos = 0
        tab_s33_o_s11_pos = 0
        tab_s34_o_s11_pos = 0
     endif

     deallocate(prob_s11_pos)
     allocate(prob_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_s11_pos'
        stop
     endif
     prob_s11_pos = 0
  else
     deallocate(ksca_CDF)
     low_mem_scattering = .false.
     allocate(ksca_CDF(0:n_grains_tot,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     mem_size = n_grains_tot * p_n_cells * p_n_lambda2_pos * 4. / 1024.**3
     if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"
     if (alloc_status > 0) then
        write(*,*) 'Allocation error ksca_CDF'
        stop
     endif
     ksca_CDF = 0
  endif ! method

  deallocate(kappa, kappa_sca, kappa_abs_LTE)
  deallocate(proba_abs_RE_LTE)
  if (lRE_nLTE) deallocate(kappa_abs_nLTE)
  if (lRE_nLTE.or.lnRE) deallocate(proba_abs_RE_LTE_p_nLTE)
  if (lnRE) deallocate(proba_abs_RE,kappa_abs_RE)
  allocate(kappa(n_cells,n_lambda2), kappa_sca(n_cells,n_lambda2), kappa_abs_LTE(n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa'
     stop
  endif
  kappa=0.0 ; kappa_sca = 0.0 ; kappa_abs_LTE=0.0
  if (lRE_nLTE) then
     allocate(kappa_abs_nLTE(n_cells,n_lambda2))
     if (alloc_status > 0) then
        write(*,*) 'Allocation error kappa_abs_nLTE'
        stop
     endif
     kappa_abs_nLTE=0.0
  endif

  if (lorigine) then
     deallocate(disk_origin, star_origin)
     allocate(disk_origin(n_lambda2, n_rad, nz, nb_proc), star_origin(n_lambda2, nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error disk_origin'
        stop
     endif
     disk_origin = 0.0
     star_origin = 0.0
  endif

  return

end subroutine realloc_step2

!**********************************************************************

subroutine realloc_ray_tracing_scattering_matrix()

  integer, parameter :: p_n_lambda2_pos = 1
  integer :: alloc_status

  ! parametrage methode de diffusion
  scattering_method = 2
  write(*,fmt='(" Using scattering method ",i1)') scattering_method
  lscattering_method1 = (scattering_method==1)

  if (allocated(tab_s11_pos)) deallocate(tab_s11_pos)
  allocate(tab_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11_pos'
     stop
  endif
  tab_s11_pos = 0

  if (lsepar_pola) then
     if (allocated(tab_s12_o_s11_pos)) deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
     allocate(tab_s12_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
          tab_s33_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
          tab_s34_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
          stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_s12_o_s11_pos'
        stop
     endif
     tab_s12_o_s11_pos = 0
     tab_s33_o_s11_pos = 0
     tab_s34_o_s11_pos = 0
  endif

  if (allocated(prob_s11_pos)) deallocate(prob_s11_pos)
  allocate(prob_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_s11_pos'
     stop
  endif
  prob_s11_pos = 0

  return

end subroutine realloc_ray_tracing_scattering_matrix

!**********************************************************************


subroutine alloc_emission_mol(imol)

  integer, intent(in) :: imol
  integer :: alloc_status, n_speed, n_speed_rt, nTrans_raytracing


  n_speed = mol(imol)%n_speed_rt ! I use the same now
  n_speed_rt = mol(imol)%n_speed_rt
  nTrans_raytracing = mol(imol)%nTrans_raytracing

  allocate(kappa_mol_o_freq(n_cells,nTrans_tot), emissivite_mol_o_freq(n_cells,nTrans_tot), &
       stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa_mol_o_freq'
     stop
  endif
  kappa_mol_o_freq=0.0
  emissivite_mol_o_freq = 0.0

  allocate(tab_nLevel(n_cells,nLevels), tab_nLevel_old(n_cells,nLevels), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_nLevel'
     stop
  endif
  tab_nLevel = 0.0
  tab_nLevel_old = 0.0

  allocate(maser_map(n_cells,nTrans_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error maser_map'
     stop
  endif
  maser_map = 0.0

  allocate(tab_v(-n_largeur_Doppler*n_speed:n_largeur_Doppler*n_speed), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_v'
     stop
  endif
  tab_v=0.0

  allocate(tab_deltaV(-n_speed:n_speed,n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_deltaV'
     stop
  endif
  tab_deltaV = 0.0

   allocate(tab_Cmb_mol(nTrans_tot), stat=alloc_status)
   if (alloc_status > 0) then
      write(*,*) 'Allocation error tab_Cmb_mol'
      stop
   endif
   tab_Cmb_mol = 0.0

   allocate(Jmol(nTrans_tot,nb_proc), stat=alloc_status)
   if (alloc_status > 0) then
      write(*,*) 'Allocation error Jmol'
      stop
   endif
   Jmol = 0.0_db

  if (ldouble_RT) then
     allocate(kappa_mol_o_freq2(n_cells,nTrans_tot), emissivite_mol_o_freq2(n_cells,nTrans_tot),&
          stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error kappa_mol2'
        stop
     endif
     kappa_mol_o_freq2=0.0
     emissivite_mol_o_freq2=0.0

     allocate(tab_nLevel2(n_cells,nLevels), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_nLevel2'
        stop
     endif
     tab_nLevel2 = 0.0

     allocate(Jmol2(nTrans_tot,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Jmol2'
        stop
     endif
     Jmol2 = 0.0_db
  endif

  ! Methode d'echantillonnage
  if (igridx > 1) then
     RT_line_method = 2 ! creation d'une carte avec pixels carres

     write(*,*) "WARNING : memory size if lots of pixels"
     allocate(spectre(igridx,igridy,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_incl,RT_n_az), &
          continu(igridx,igridy,nTrans_raytracing,RT_n_incl,RT_n_az), stars_map(igridx,igridy,1), stat=alloc_status)
  else
     RT_line_method = 1 ! utilisation de pixels circulaires
     allocate(spectre(1,1,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_incl,RT_n_az), &
          continu(1,1,nTrans_raytracing,RT_n_incl,RT_n_az), stars_map(1,1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre'
     stop
  endif
  spectre=0.0
  continu=0.0
  stars_map=0.0

  return

end subroutine alloc_emission_mol

!**********************************************************************

subroutine dealloc_emission_mol()

  ! Dealloue ce qui n'a pas ete libere par  clean_mem_dust_mol
  deallocate(tab_lambda, tab_delta_lambda, tab_lambda_inf, tab_lambda_sup)
  deallocate(kappa, kappa_sca, emissivite_dust)
  if (allocated(spectre_etoiles)) deallocate(spectre_etoiles)
  if (allocated(spectre_etoiles_cumul)) deallocate(spectre_etoiles_cumul)
  if (allocated(CDF_E_star)) deallocate(CDF_E_star)
  if (allocated(prob_E_star)) deallocate(prob_E_star)
  if (allocated(E_stars)) deallocate(E_stars)
  if (allocated(E_ISM)) deallocate(E_ISM)

  deallocate(Level_energy,poids_stat_g,j_qnb,Aul,fAul,Bul,fBul,Blu,fBlu,transfreq, &
       itransUpper,itransLower,nCollTrans,nCollTemps,collTemps,collBetween, &
       iCollUpper,iCollLower,indice_Trans)

  deallocate(kappa_mol_o_freq, emissivite_mol_o_freq, tab_nLevel, tab_nLevel_old, &
       tab_v, tab_deltaV, spectre,continu, stars_map, tab_Cmb_mol, Jmol, maser_map)

  if (ldouble_RT) deallocate(kappa_mol_o_freq2, emissivite_mol_o_freq2, tab_nLevel2, Jmol2)

  deallocate(I0, I0c, tab_speed_rt) ! besoin de dealoue tab_speed_rt pour plusieyrs ray
  if (lorigine) deallocate(origine_mol)

  return

end subroutine dealloc_emission_mol

!**********************************************************************

subroutine sig_handler(sig)

  integer, intent(in) ::  sig

  select case(sig)
  case(2)
     write (*,*) 'mcfost : SIGINT Caught'
     stop
  case(15)
     write (*,*) 'mcfost : SIGTERM Caught'
     stop
  case default
     write (*,*) 'mcfost : signal caught :', sig
     read(*,*)
     return
  end select

  !call exit(1)
  stop

end subroutine sig_handler

!*************************************************

end module mem
