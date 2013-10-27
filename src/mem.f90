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

  allocate(stream(nb_proc), gauss_random_saved(nb_proc), lgauss_random_saved(nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  gauss_random_saved = 0.0_db
  lgauss_random_saved = .false.

  allocate(n_phot_envoyes(n_lambda), n_phot_envoyes_loc(n_lambda),  stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_phot_envoyes'
     stop
  endif
  n_phot_envoyes = 0.0 ; n_phot_envoyes_loc = 0.0

  !***************************************************
  ! Tableaux relatifs a la propagation des paquets : othin ...
  !***************************************************
  allocate(n_cell_traversees(nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_cell_traversees'
     stop
  endif
  n_cell_traversees = 0

  allocate(tab_cell_r(nb_proc,n_cell_max), tab_cell_z(nb_proc,n_cell_max),  stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_cell_r'
     stop
  endif
  tab_cell_r = 0 ;  tab_cell_z = 0


  allocate(tab_length(nb_proc,n_cell_max), tab_tau(nb_proc,0:n_cell_max),&
       tab_length_tot(nb_proc,0:n_cell_max), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_length'
     stop
  endif
  tab_length = 0.0 ; tab_tau = 0.0 ; tab_length_tot=0.0!tab_tau(id:0)=0.0


  allocate(tab_x0(nb_proc,n_cell_max), tab_y0(nb_proc,n_cell_max), tab_z0(nb_proc,n_cell_max),stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_x0'
     stop
  endif
  tab_x0 = 0.0 ; tab_y0 = 0.0 ; tab_z0 = 0.0

  ! ... ou integration dichotomique
  linteg_dic = .false.
  if (linteg_dic) then
     allocate(delta0(0:n_rad), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error delta0'
        stop
     endif
     delta0 = 0.0
  endif


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

  if (l3D) then
     allocate(r_grid(n_rad,-nz:nz), z_grid(n_rad,-nz:nz), phi_grid(n_az), stat=alloc_status)
  else
     allocate(r_grid(n_rad,nz), z_grid(n_rad,nz), phi_grid(n_az), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_lim'
     stop
  endif
  r_grid=0.0; z_grid=0.0 ; phi_grid = 0.0

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

  allocate(zmax(n_rad),volume(n_rad), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error zmax, volume'
     stop
  endif
  zmax = 0.0 ; volume=0.0

  if (l3D) then
     allocate(masse(n_rad,-nz:nz,n_az), stat=alloc_status)
  else
     allocate(masse(n_rad,nz,1), stat=alloc_status)
  endif

  if (alloc_status > 0) then
     write(*,*) 'Allocation error masse'
     stop
  endif
  masse = 0.0

  if (l3D) then
     allocate(densite_pouss(n_rad,-nz-1:nz+1,n_az,n_grains_tot), stat=alloc_status)
  else
     allocate(densite_pouss(n_rad,nz+1,1,n_grains_tot), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error densite_pouss'
     stop
  endif
  densite_pouss = 0.0

  if (l3D) then
     allocate(l_dark_zone(0:n_rad+1,-nz-1:nz+1,n_az), ri_in_dark_zone(n_az), ri_out_dark_zone(n_az),&
          zj_sup_dark_zone(n_rad,n_az), zj_inf_dark_zone(n_rad,n_az), stat=alloc_status)
  else
     allocate(l_dark_zone(0:n_rad+1,0:nz+1,1), ri_in_dark_zone(n_az), ri_out_dark_zone(n_az),&
          zj_sup_dark_zone(n_rad,n_az), stat=alloc_status)
  endif
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
  allocate(nbre_grains(n_grains_tot), r_grain(n_grains_tot),  S_grain(n_grains_tot), M_grain(n_grains_tot), &
       r_core(n_grains_tot), is_grain_PAH(n_grains_tot), grain(n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_grain'
     stop
  endif
  nbre_grains = 0.0   ; r_core=0.0
  r_grain=0.0 ; S_grain=0.0 ; M_grain=0.0
  is_grain_PAH=.false.

  allocate(tab_albedo(n_lambda,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo'
     stop
  endif
  tab_albedo = 0

  allocate(q_ext(n_lambda,n_grains_tot), q_sca(n_lambda,n_grains_tot), &
       q_abs(n_lambda,n_grains_tot), q_geo(n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error q_ext'
     stop
  endif
  q_ext = 0 ; q_sca = 0 ; q_abs = 0 ; q_geo =0


  allocate(tab_g(n_lambda,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_g'
     stop
  endif
  tab_g = 0


  ! **************************************************
  ! Tableaux relatifs aux prop en fct de lambda
  ! **************************************************
  allocate(E_stars(n_lambda), E_disk(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error E_stars'
     stop
  endif
  E_stars = 0.0
  E_disk = 0.0

  allocate(frac_E_stars(n_lambda), E_totale(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error frac_E_stars'
     stop
  endif
  frac_E_stars = 0.0 ; E_totale = 0.0

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
  if (l3D) then
     allocate(amax_reel(n_lambda,n_rad,-nz-1:nz+1,n_az), kappa(n_lambda,n_rad,-nz-1:nz+1,n_az), &
          kappa_abs_eg(n_lambda,n_rad,-nz-1:nz+1,n_az), proba_abs_RE(n_lambda,n_rad,-nz-1:nz+1,n_az), &
          stat=alloc_status)
  else
     allocate(amax_reel(n_lambda,n_rad,nz+1,1), kappa(n_lambda,n_rad,nz+1,1), &
          kappa_abs_eg(n_lambda,n_rad,nz+1,1), proba_abs_RE(n_lambda,n_rad,nz+1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa'
     stop
  endif
  amax_reel = 0.0 ; kappa=0.0 ; kappa_abs_eg=0.0 ; proba_abs_RE=0.0


  if (l3D) then
     allocate(tab_albedo_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az), tab_g_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az),&
          stat=alloc_status)
  else
     allocate(tab_albedo_pos(n_lambda,p_n_rad,p_nz,1), tab_g_pos(n_lambda,p_n_rad,p_nz,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos, tab_g_pos'
     stop
  endif
  tab_albedo_pos = 0
  tab_g_pos = 0.0


  ! **************************************************
  ! Tableaux relatifs aux prop optiques des cellules ou des grains
  ! **************************************************
  if (scattering_method == 2) then ! prop par cellule
     if (l3D) then
        allocate(tab_s11_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
     else
        allocate(tab_s11_pos(n_lambda,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_s11_pos'
        stop
     endif
     tab_s11_pos = 0

     if (lsepar_pola) then
        if (l3D) then
           allocate(tab_s12_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
        else
           allocate(tab_s12_pos(n_lambda,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
        endif
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s12_pos'
           stop
        endif
        tab_s12_pos = 0

        if (l3D) then
           allocate(tab_s33_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
        else
           allocate(tab_s33_pos(n_lambda,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
        endif
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s33_pos'
           stop
        endif
        tab_s33_pos = 0

        if (l3D) then
           allocate(tab_s34_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
        else
           allocate(tab_s34_pos(n_lambda,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
        endif
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s34_pos'
           stop
        endif
        tab_s34_pos = 0
     endif

     if (l3D) then
        allocate(prob_s11_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
     else
        allocate(prob_s11_pos(n_lambda,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_s11_pos'
        stop
     endif
     prob_s11_pos = 0
     p_n_lambda = 1
  else ! prop par grains
     p_n_lambda = n_lambda

     if (l3D) then
        allocate(probsizecumul(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:n_grains_tot), stat=alloc_status)
     else
        allocate(probsizecumul(n_lambda,p_n_rad,p_nz,1,0:n_grains_tot), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error probsizecumul'
        stop
     endif
     probsizecumul = 0

     if (l3D) then
        allocate(ech_prob(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:n_prob+1), stat=alloc_status)
     else
        allocate(ech_prob(n_lambda,p_n_rad,p_nz,1,0:n_prob+1), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error ech_prob'
        stop
     endif
     ech_prob = 0

     if (l3D) then
        allocate(valeur_prob(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:n_prob+1), stat=alloc_status)
     else
        allocate(valeur_prob(n_lambda,p_n_rad,p_nz,1,0:n_prob+1), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error valeur_prob'
        stop
     endif
     valeur_prob = 0

     allocate(xspline(n_lambda,p_n_rad,p_nz,0:n_prob+1), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error xspline'
        stop
     endif
     xspline = 0
  endif ! method


  ! **************************************************
  ! tableaux relatifs aux prop optiques des grains
  ! **************************************************
  allocate(tab_s11(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11'
     stop
  endif
  tab_s11 = 0

  allocate(tab_s12(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s12'
     stop
  endif
  tab_s12 = 0

  allocate(tab_s33(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s33'
     stop
  endif
  tab_s33 = 0

  allocate(tab_s34(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s34'
     stop
  endif
  tab_s34 = 0

  if (laggregate) then
     allocate(tab_mueller(p_n_lambda,n_grains_tot,4,4,0:nang_scatt), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_mueller'
        stop
     endif
     tab_mueller = 0
  endif

  allocate(prob_s11(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
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


  if (l3D) then
     allocate(prob_E_cell(n_lambda,0:n_rad*2*nz*n_az), stat=alloc_status)
  else
     allocate(prob_E_cell(n_lambda,0:n_rad*nz), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_cell'
     stop
  endif
  prob_E_cell = 0.0

  if (lweight_emission) then
     allocate(weight_proba_emission(n_rad,nz), correct_E_emission(n_rad,nz), stat=alloc_status)
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
  if (l3D) then
     allocate(Temperature(n_rad,-nz:nz,n_az), Temperature_old(n_rad,-nz:nz,n_az), stat=alloc_status)
  else
     allocate(Temperature(n_rad,nz,1), Temperature_old(n_rad,nz,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error Temperature'
     stop
  endif
  Temperature = 0.0 ; Temperature_old=0.0


  if (lRE_nLTE) then
     allocate(Temperature_1grain(n_rad,nz,grain_RE_nLTE_start:grain_RE_nLTE_end),stat=alloc_status)
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

     allocate(Proba_Temperature(n_T,n_rad,nz,grain_nRE_start:grain_nRE_end), &
          Temperature_1grain_nRE(n_rad,nz,grain_nRE_start:grain_nRE_end), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Proba_Temperature'
        stop
     endif
     Proba_Temperature=0.0
     Temperature_1grain_nRE=0.0

     allocate(l_RE(n_rad,nz,grain_nRE_start:grain_nRE_end),stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error l_RE'
        stop
     endif
     l_RE=.false.
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

     if (l3D) then
        allocate(log_frac_E_em(n_rad,-nz:nz,n_az,n_T), stat=alloc_status)
     else
        allocate(log_frac_E_em(n_rad,nz,1,n_T), stat=alloc_status)
     endif
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


     if (l3D) then
        allocate(xKJ_abs(n_rad,-nz:nz,n_az,nb_proc), xE_abs(n_rad,-nz:nz,n_az,nb_proc), &
             xJ_abs(n_lambda,n_rad,nz,nb_proc), nbre_reemission(n_rad,-nz:nz,n_az,nb_proc),&
             stat=alloc_status)
     else
        allocate(xKJ_abs(n_rad,nz,1,nb_proc), xE_abs(n_rad,nz,1,nb_proc), &
             xJ_abs(n_lambda,n_rad,nz,nb_proc), nbre_reemission(n_rad,nz,1,nb_proc), &
             stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error xKJ_abs'
        stop
     endif
     xKJ_abs = 0.0 ; xE_abs = 0.0 ; xJ_abs=0.0 ; nbre_reemission = 0.0

     if (l3D) then
        allocate(E0(n_rad,-nz:nz,n_az), J0(n_lambda,n_rad,-nz:nz,n_az),  stat=alloc_status)
     else
        allocate(E0(n_rad,nz,1),J0(n_lambda,n_rad,nz,1),  stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error E0'
        stop
     endif
     E0 = 0.0
     J0 = 0.0

     if (l3D) then
        allocate(xT_ech(nb_proc,n_rad,-nz:nz,n_az), stat=alloc_status)
     else
        allocate(xT_ech(nb_proc,n_rad,nz,1), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error xT_ech'
        stop
     endif
     xT_ech = 2

     if (l3D) then
        allocate(prob_delta_T(p_n_rad,-p_nz:p_nz,p_n_az,n_T,n_lambda), stat=alloc_status)
     else
        allocate(prob_delta_T(p_n_rad,p_nz,1,n_T,n_lambda), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_delta_T'
        stop
     endif
     prob_delta_T = 0


     if (lRE_nLTE) then
        allocate(prob_kappa_abs_1grain(n_lambda,n_rad,nz+1,0:n_grains_RE_nLTE),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error kappa_abs_1grain'
           stop
        endif
        prob_kappa_abs_1grain=0.0

        allocate(prob_delta_T_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_T,n_lambda),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error prob_delta_T_1grain'
           stop
        endif
        prob_delta_T_1grain=0.0

        allocate(log_frac_E_em_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_T),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error log_frac_E_em_1grain'
           stop
        endif
        log_frac_E_em_1grain=0.0

        allocate(xE_abs_1grain(n_rad,nz,grain_RE_nLTE_start:grain_RE_nLTE_end,nb_proc),stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error xE_abs_1grain'
           stop
        endif
        xE_abs_1grain = 0.0

        allocate(xT_ech_1grain(nb_proc,n_rad,nz,grain_RE_nLTE_start:grain_RE_nLTE_end),stat=alloc_status)
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

        !allocate(log_frac_E_em_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_T),stat=alloc_status)
        !if (alloc_status > 0) then
        !   write(*,*) 'Allocation error log_frac_E_em_1grain'
        !   stop
        !endif
        !log_frac_E_em_1grain=0.0

        allocate(Temperature_1grain_nRE_old(n_rad,nz,grain_nRE_start:grain_nRE_end), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error Proba_Temperature'
           stop
        endif
        Temperature_1grain_nRE_old =0.0

        allocate(Tpeak_old(n_rad,nz,grain_nRE_start:grain_nRE_end), &
             maxP_old(n_rad,nz,grain_nRE_start:grain_nRE_end), &
             stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error Tpeak'
           stop
        endif
        Tpeak_old=0
        maxP_old=0.

        if (lRE_nlTE) then
           allocate(Temperature_1grain_old(n_rad,nz,grain_RE_nLTE_start:grain_RE_nLTE_end),stat=alloc_status)
           if (alloc_status > 0) then
              write(*,*) 'Allocation error Temperature_old'
              stop
           endif
           Temperature_old=0.
        endif
     endif

  endif ! lTemp

  if (lnRE) then
     allocate(Emissivite_nRE_old(n_lambda,n_rad,nz,1), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Proba_Temperature'
        stop
     endif
     Emissivite_nRE_old = 0.0
  endif



 ! if (lProDiMo) then
 !    allocate(J_prodimo(n_lambda,n_rad,nz),stat=alloc_status)
 !    if (alloc_status > 0) then
 !       write(*,*) 'Allocation error J_prodimo'
 !       stop
 !    endif
 !    J_prodimo = 0.
 ! endif

  ! **************************************************
  ! Tableaux relatifs aux SEDs
  ! **************************************************
  if (lTemp.or.lsed) then
     allocate(sed(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed'
        stop
     endif
     sed = 0.0

     allocate(sed_q(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_q'
        stop
     endif
     sed_q = 0.0

     allocate(sed_u(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_u'
        stop
     endif
     sed_u = 0.0

     allocate(sed_v(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_v'
        stop
     endif
     sed_v = 0.0

     allocate(sed_star(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_star'
        stop
     endif
     sed_star = 0.0

     allocate(sed_star_scat(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_star_scat'
        stop
     endif
     sed_star_scat = 0.0

     allocate(sed_disk(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_disk'
        stop
     endif
     sed_disk = 0.0

     allocate(sed_disk_scat(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sed_disk_scat'
        stop
     endif
     sed_disk_scat = 0.0

     allocate(n_phot_sed(nb_proc,n_lambda,N_thet,N_phi), stat=alloc_status)
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
  if (lmono0) then
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

     if (n_cartes > 1) then
        allocate(STOKEI1(n_lambda,IGRIDX2(1),IGRIDY2(1),capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI1'
           stop
        endif
        STOKEI1 = 0.0
     endif
     if (n_cartes > 2) then
        allocate(STOKEI2(n_lambda,IGRIDX2(2),IGRIDY2(2),capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI1'
           stop
        endif
        STOKEI2 = 0.0
     endif
     if (n_cartes > 3) then
        allocate(STOKEI3(n_lambda,IGRIDX2(3),IGRIDY2(3),capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI1'
           stop
        endif
        STOKEI3 = 0.0
     endif

  endif ! lmono0

  ! **************************************************
  ! Tableaux relatifs a l'emission moleculaire
  ! **************************************************
  if (lemission_mol) then
     allocate(Tcin(n_rad,nz,1), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Tcin'
        stop
     endif
     Tcin=0.0

     allocate(tab_abundance(n_rad,nz), lcompute_molRT(n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_abundance'
        stop
     endif
     tab_abundance = 0.0
     lcompute_molRT = .true.

     allocate(vfield(n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error vfield'
        stop
     endif
     vfield=0.0 !; vx=0.0 ; vy=0.0

     allocate(v_turb(n_rad, nz), v_line(n_rad, nz), deltaVmax(n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sigma2'
        stop
     endif
     v_turb = 0.0 ; v_line = 0.0 ;   deltaVmax = 0.0

     allocate(tab_dnu_o_freq(n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_dnu'
        stop
     endif
     tab_dnu_o_freq=0.0

     allocate(norme_phiProf_m1(n_rad,nz), sigma2_phiProf_m1(n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error norme_phiProf_m1'
        stop
     endif
     norme_phiProf_m1 = 0.0 ; sigma2_phiProf_m1 = 0.0
  endif ! lemission_mol

  if (l3D) then
     allocate(densite_gaz(n_rad,-nz:nz,n_az), masse_gaz(n_rad,-nz:nz,n_az), stat=alloc_status)
  else
     allocate(densite_gaz(n_rad,-nz:nz,n_az), masse_gaz(n_rad,-nz:nz,n_az), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error densite_gaz'
     stop
  endif
  densite_gaz = 0.0 ; masse_gaz = 0.0

  return

end subroutine alloc_dynamique

!**********************************************************************

subroutine dealloc_em_th()

  deallocate(n_phot_envoyes,n_phot_envoyes_loc)

  deallocate(n_cell_traversees,tab_cell_r,tab_cell_z,tab_length,tab_tau,tab_length_tot)

  deallocate(tab_x0,tab_y0,tab_z0)

  deallocate(tab_albedo,q_ext,q_sca,q_abs,q_geo,tab_g)

  deallocate(E_stars,E_disk,frac_E_stars,E_totale)

  deallocate(spectre_etoiles_cumul,spectre_etoiles, spectre_emission_cumul)

  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2)

  deallocate(amax_reel,kappa,kappa_abs_eg,proba_abs_RE)

  deallocate(tab_albedo_pos,tab_g_pos)

  if (scattering_method == 2) then ! prop par cellule
     deallocate(tab_s11_pos,prob_s11_pos)
     if (lsepar_pola) deallocate(tab_s12_pos,tab_s33_pos,tab_s34_pos)
  else ! prop par grains
     deallocate(probsizecumul,ech_prob,valeur_prob,xspline)
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
        deallocate(xKJ_abs,xE_abs,xJ_abs,nbre_reemission)
        deallocate(E0,J0,xT_ech,prob_delta_T)
     endif

     if (lRE_nLTE) then
        deallocate(prob_kappa_abs_1grain,prob_delta_T_1grain,log_frac_E_em_1grain)
        deallocate(xE_abs_1grain,xT_ech_1grain)
     endif

     if (lnRE) then
        deallocate(frac_E_em_1grain_nRE,log_frac_E_em_1grain_nRE)
        deallocate(Temperature_1grain_nRE_old)
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

  integer :: alloc_status

  allocate(tab_lambda(n_lambda), tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), &
       tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_lambda (realloc)'
     stop
  endif
  tab_lambda=0.0 ; tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda= 0.0 ; tab_amu1=0.0 ; tab_amu2=0.0

  allocate(tab_albedo(n_lambda,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo (realloc)'
     stop
  endif
  tab_albedo = 0

  allocate(q_ext(n_lambda,n_grains_tot), q_sca(n_lambda,n_grains_tot), &
       q_abs(n_lambda,n_grains_tot), tab_g(n_lambda,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error q_ext (realloc)'
     stop
  endif
  q_ext = 0 ; q_sca = 0 ; q_abs = 0 ; tab_g = 0


  allocate(prob_s11(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_s11 (realloc)'
     stop
  endif
  prob_s11 = 0

  allocate(tab_s11(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11 (realloc)'
     stop
  endif
  tab_s11 = 0

  allocate(tab_s12(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s12 (realloc)'
     stop
  endif
  tab_s12 = 0

  allocate(tab_s33(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s33 (realloc)'
     stop
  endif
  tab_s33 = 0

  allocate(tab_s34(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s34 (realloc)'
     stop
  endif
  tab_s34 = 0

  ! TODO : cette partie prend bcp de memoire
  if (l3D) then
     allocate(probsizecumul(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:n_grains_tot), stat=alloc_status)
  else
     allocate(probsizecumul(n_lambda,p_n_rad,p_nz,1,0:n_grains_tot), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error probsizecumul (realloc)'
     stop
  endif
  probsizecumul = 0

  ! Tableaux relatifs aux prop optiques des cellules
  if (l3D) then
     allocate(kappa(n_lambda,n_rad,-nz-1:nz+1,n_az),kappa_abs_eg(n_lambda,n_rad,-nz-1:nz+1,n_az), &
          kappa_sca(n_lambda,n_rad,-nz-1:nz+1,n_az), &
          emissivite_dust(n_lambda,n_rad,-nz-1:nz+1,n_az),proba_abs_RE(n_lambda,n_rad,-nz-1:nz+1,n_az), &
          amax_reel(n_lambda,n_rad,-nz-1:nz+1,n_az), stat=alloc_status)
  else
     allocate(kappa(n_lambda,n_rad,nz+1,1),kappa_abs_eg(n_lambda,n_rad,nz+1,1), &
          kappa_sca(n_lambda,n_rad,nz+1,1), &
          emissivite_dust(n_lambda,n_rad,nz+1,1),proba_abs_RE(n_lambda,n_rad,nz+1,1),&
          amax_reel(n_lambda,n_rad,nz+1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa (realloc)'
     stop
  endif
  kappa = 0.0 ; kappa_abs_eg = 0.0 ; kappa_sca = 0.0 ; emissivite_dust = 0.0

  if (l3D) then
     allocate(tab_albedo_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az), tab_g_pos(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az),&
          stat=alloc_status)
  else
     allocate(tab_albedo_pos(n_lambda,p_n_rad,p_nz,1), tab_g_pos(n_lambda,p_n_rad,p_nz,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos, tab_g_pos (realloc)'
     stop
  endif
  tab_albedo_pos = 0
  tab_g_pos = 0.0

  if (l3D) then
     allocate(ech_prob(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:n_prob+1), stat=alloc_status)
  else
     allocate(ech_prob(n_lambda,p_n_rad,p_nz,1,0:n_prob+1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error ech_prob (realloc)'
     stop
  endif
  ech_prob = 0

  if (l3D) then
     allocate(valeur_prob(n_lambda,p_n_rad,-p_nz:p_nz,p_n_az,0:n_prob+1), stat=alloc_status)
  else
     allocate(valeur_prob(n_lambda,p_n_rad,p_nz,1,0:n_prob+1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error valeur_prob (realloc)'
     stop
  endif
  valeur_prob = 0

  return

end subroutine realloc_dust_mol

!******************************************************************************

subroutine clean_mem_dust_mol()

  integer :: alloc_status

  ! Ne reste que kappa et emissivite_dust
  deallocate(tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda, tab_amu1, tab_amu2)
  deallocate(tab_albedo)
  deallocate(q_ext, q_sca, q_abs, tab_g)
  deallocate(prob_s11,tab_s11,tab_s12,tab_s33,tab_s34,probsizecumul)
  deallocate(kappa_abs_eg,proba_abs_RE,amax_reel)
  deallocate(tab_albedo_pos, tab_g_pos)
  deallocate(ech_prob,valeur_prob)

  return

end subroutine clean_mem_dust_mol

!******************************************************************************

subroutine realloc_step2()

  integer :: alloc_status


  if (scattering_method == 2) then ! prop par cellule
     p_n_lambda = 1
  else ! prop par grains
     p_n_lambda = n_lambda2
  endif

  ! Liberation memoire
  if (ltemp) then
     if (lRE_LTE)  deallocate(prob_delta_T, log_frac_E_em, xT_ech)
     if (lRE_nLTE) deallocate(prob_kappa_abs_1grain, prob_delta_T_1grain, log_frac_E_em_1grain,xT_ech_1grain)
     deallocate(xJ_abs, xKJ_abs, nbre_reemission)
  endif

  if (lProDiMo) then
     allocate(xJ_abs(n_lambda2,n_rad,nz,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error xJ_abs in realloc_step2'
        stop
     endif
     xJ_abs = 0.0
  endif

  ! Liberation memoire step1 et reallocation step 2
  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2)
  allocate(tab_lambda(n_lambda2),tab_lambda_inf(n_lambda2),tab_lambda_sup(n_lambda2),tab_delta_lambda(n_lambda2),&
       tab_amu1(n_lambda2,n_pop),tab_amu2(n_lambda2,n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_lambda in init_lambda2'
     stop
  endif
  tab_lambda=0.0 ; tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda=0.0
  tab_amu1=0.0 ; tab_amu2=0.0

  deallocate(sed)
  allocate(sed(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed'
     stop
  endif
  sed = 0.0



  deallocate(sed_q)
  allocate(sed_q(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_q'
     stop
  endif
  sed_q = 0.0

  deallocate(sed_u)
  allocate(sed_u(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_u'
     stop
  endif
  sed_u = 0.0

  deallocate(sed_v)
  allocate(sed_v(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_v'
     stop
  endif
  sed_v = 0.0

  deallocate(sed_star)
  allocate(sed_star(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_star'
     stop
  endif
  sed_star = 0.0

  deallocate(sed_star_scat)
  allocate(sed_star_scat(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_star_scat'
     stop
  endif
  sed_star_scat = 0.0

  deallocate(sed_disk)
  allocate(sed_disk(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_disk'
     stop
  endif
  sed_disk = 0.0

  deallocate(sed_disk_scat)
  allocate(sed_disk_scat(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sed_disk_scat'
     stop
  endif
  sed_disk_scat = 0.0

  deallocate(n_phot_sed)
  allocate(n_phot_sed(nb_proc,n_lambda2,N_thet,N_phi), stat=alloc_status)
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

  deallocate(n_phot_envoyes, n_phot_envoyes_loc)
  allocate(n_phot_envoyes(n_lambda2), n_phot_envoyes_loc(n_lambda2),  stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_phot_envoyes'
     stop
  endif
  n_phot_envoyes = 0.0 ; n_phot_envoyes_loc = 0.0

  deallocate(prob_E_cell)
  if (l3D) then
     allocate(prob_E_cell(n_lambda2,0:n_rad*2*nz*n_az), stat=alloc_status)
  else
     allocate(prob_E_cell(n_lambda2,0:n_rad*nz), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_cell'
     stop
  endif
  prob_E_cell = 0.0

  deallocate(prob_E_star)
  allocate(prob_E_star(n_lambda2,0:n_etoiles), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error prob_E_star'
     stop
  endif
  prob_E_star = 0.0

  deallocate(spectre_etoiles_cumul, spectre_etoiles, spectre_emission_cumul)
  allocate(spectre_etoiles_cumul(0:n_lambda2),spectre_etoiles(n_lambda2),spectre_emission_cumul(0:n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre_etoile'
     stop
  endif
  spectre_etoiles_cumul = 0.0
  spectre_etoiles = 0.0
  spectre_emission_cumul = 0.0

  deallocate(E_stars, E_disk)
  allocate(E_stars(n_lambda2), E_disk(n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error E_stars'
     stop
  endif
  E_stars = 0.0
  E_disk = 0.0

  deallocate(frac_E_stars, E_totale)
  allocate(frac_E_stars(n_lambda2), E_totale(n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error frac_E_stars'
     stop
  endif
  frac_E_stars = 0.0 ; E_totale = 0.0

  deallocate(tab_albedo)
  allocate(tab_albedo(n_lambda2,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo'
     stop
  endif
  tab_albedo = 0

  deallocate(q_ext)
  allocate(q_ext(n_lambda2,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error q_ext'
     stop
  endif
  q_ext = 0

  deallocate(q_sca)
  allocate(q_sca(n_lambda2,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error q_sca'
     stop
  endif
  q_sca = 0

   deallocate(q_abs)
  allocate(q_abs(n_lambda2,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error q_abs'
     stop
  endif
  q_abs = 0

  deallocate(tab_g)
  allocate(tab_g(n_lambda2,n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_g'
     stop
  endif
  tab_g = 0


  deallocate(tab_albedo_pos)
  if (l3D) then
     allocate(tab_albedo_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az), stat=alloc_status)
  else
     allocate(tab_albedo_pos(n_lambda2,p_n_rad,p_nz,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos'
     stop
  endif
  tab_albedo_pos = 0

  deallocate(tab_g_pos)
  if (l3D) then
     allocate(tab_g_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az), stat=alloc_status)
  else
     allocate(tab_g_pos(n_lambda2,p_n_rad,p_nz,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_g_pos'
     stop
  endif
  tab_g_pos = 0

  deallocate(tab_s11)
  allocate(tab_s11(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s11'
     stop
  endif
  tab_s11 = 0

  deallocate(tab_s12)
  allocate(tab_s12(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s12'
     stop
  endif
  tab_s12 = 0

  deallocate(tab_s33)
  allocate(tab_s33(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s33'
     stop
  endif
  tab_s33 = 0

  deallocate(tab_s34)
  allocate(tab_s34(p_n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_s34'
     stop
  endif
  tab_s34 = 0

  if (laggregate) then
     deallocate(tab_mueller)
     allocate(tab_mueller(n_lambda2,n_grains_tot,4,4,0:nang_scatt), stat=alloc_status)
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
     if (l3D) then
        allocate(tab_s11_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
     else
        allocate(tab_s11_pos(n_lambda2,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_s11_pos'
        stop
     endif
     tab_s11_pos = 0

     if (lsepar_pola) then
        deallocate(tab_s12_pos)
        if (l3D) then
           allocate(tab_s12_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
        else
           allocate(tab_s12_pos(n_lambda2,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
        endif
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s12_pos'
           stop
        endif
        tab_s12_pos = 0

        deallocate(tab_s33_pos)
        if (l3D) then
           allocate(tab_s33_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
        else
           allocate(tab_s33_pos(n_lambda2,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
        endif
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s33_pos'
           stop
        endif
        tab_s33_pos = 0

        deallocate(tab_s34_pos)
        if (l3D) then
           allocate(tab_s34_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
        else
           allocate(tab_s34_pos(n_lambda2,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
        endif
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_s34_pos'
           stop
        endif
        tab_s34_pos = 0
     endif

     deallocate(prob_s11_pos)
     if (l3D) then
        allocate(prob_s11_pos(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:nang_scatt), stat=alloc_status)
     else
        allocate(prob_s11_pos(n_lambda2,p_n_rad,p_nz,1,0:nang_scatt), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_s11_pos'
        stop
     endif
     prob_s11_pos = 0
  else
     deallocate(probsizecumul)
     if (l3D) then
        allocate(probsizecumul(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:n_grains_tot), stat=alloc_status)
     else
        allocate(probsizecumul(n_lambda2,p_n_rad,p_nz,1,0:n_grains_tot), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error probsizecumul'
        stop
     endif
     probsizecumul = 0

     deallocate(ech_prob)
     if (l3D) then
        allocate(ech_prob(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:n_prob+1), stat=alloc_status)
     else
        allocate(ech_prob(n_lambda2,p_n_rad,p_nz,1,0:n_prob+1), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error ech_prob'
        stop
     endif
     ech_prob = 0

     deallocate(valeur_prob)
     if (l3D) then
        allocate(valeur_prob(n_lambda2,p_n_rad,-p_nz:p_nz,p_n_az,0:n_prob+1), stat=alloc_status)
     else
        allocate(valeur_prob(n_lambda2,p_n_rad,p_nz,1,0:n_prob+1), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error valeur_prob'
        stop
     endif
     valeur_prob = 0
  endif ! method

  deallocate(amax_reel, kappa, kappa_abs_eg, proba_abs_RE)
  if (l3D) then
     allocate(amax_reel(n_lambda2,n_rad,-nz-1:nz+1,n_az), kappa(n_lambda2,n_rad,-nz-1:nz+1,n_az), &
          kappa_abs_eg(n_lambda2,n_rad,-nz-1:nz+1,n_az), proba_abs_RE(n_lambda2,n_rad,-nz-1:nz+1,n_az), &
          stat=alloc_status)
 else
     allocate(amax_reel(n_lambda2,n_rad,nz+1,1), kappa(n_lambda2,n_rad,nz+1,1), &
          kappa_abs_eg(n_lambda2,n_rad,nz+1,1), proba_abs_RE(n_lambda2,n_rad,nz+1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa'
     stop
  endif
  amax_reel = 0.0 ; kappa=0.0 ; kappa_abs_eg=0.0

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

subroutine alloc_emission_mol(imol)

  integer, intent(in) :: imol
  integer :: alloc_status, n_speed, n_speed_rt, nTrans_raytracing


  n_speed = mol(imol)%n_speed
  n_speed_rt = mol(imol)%n_speed_rt
  nTrans_raytracing = mol(imol)%nTrans_raytracing

  allocate(kappa_mol_o_freq(n_rad,nz+1,nTrans_tot), emissivite_mol_o_freq(n_rad,nz+1,nTrans_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa_mol_o_freq'
     stop
  endif
  kappa_mol_o_freq=0.0
  emissivite_mol_o_freq = 0.0

  allocate(tab_nLevel(n_rad,nz,nLevels), tab_nLevel_old(n_rad,nz,nLevels), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_nLevel'
     stop
  endif
  tab_nLevel = 0.0
  tab_nLevel_old = 0.0


  allocate(tab_v(-n_largeur_Doppler*n_speed:n_largeur_Doppler*n_speed), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_v'
     stop
  endif
  tab_v=0.0

  allocate(tab_deltaV(-n_speed:n_speed,n_rad,nz), stat=alloc_status)
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
     allocate(kappa_mol_o_freq2(n_rad,nz+1,nTrans_tot), emissivite_mol_o_freq2(n_rad,nz+1,nTrans_tot), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error kappa_mol2'
        stop
     endif
     kappa_mol_o_freq2=0.0
     emissivite_mol_o_freq2=0.0

     allocate(tab_nLevel2(n_rad,nz,nLevels), stat=alloc_status)
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
     allocate(spectre(igridx,igridy,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_ibin), &
          continu(igridx,igridy,nTrans_raytracing,RT_n_ibin), stat=alloc_status)
  else
     RT_line_method = 1 ! utilisation de pixels circulaires
     allocate(spectre(1,1,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_ibin), &
          continu(1,1,nTrans_raytracing,RT_n_ibin), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre'
     stop
  endif
  spectre=0.0
  continu=0.0

  allocate(maser_map(n_rad,nz,nTrans_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error maser_map'
     stop
  endif
  maser_map = 0.0

  return

end subroutine alloc_emission_mol

!**********************************************************************

subroutine dealloc_emission_mol()

  deallocate(kappa,emissivite_dust)

  deallocate(Level_energy,poids_stat_g,j_qnb,Aul,fAul,Bul,fBul,Blu,fBlu,transfreq, &
       itransUpper,itransLower,nCollTrans,nCollTemps,collTemps,collBetween, &
       iCollUpper,iCollLower,indice_Trans)

  deallocate(kappa_mol_o_freq, emissivite_mol_o_freq, tab_nLevel, tab_nLevel_old, &
       tab_v, tab_deltaV, spectre,continu, tab_Cmb_mol, Jmol, maser_map)

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
