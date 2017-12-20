module mem

  use parametres
  use grains
  use opacity
  use resultats
  use prop_star
  use naleat
  use molecular_emission
  use ray_tracing
  use utils

  implicit none

  contains

subroutine allocate_densities(n_cells_max)

  integer, intent(in), optional :: n_cells_max

  integer ::  alloc_status, Nc

  if (present(n_cells_max)) then
     Nc = n_cells_max
  else
     Nc = n_cells
  endif

  allocate(masse(Nc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error mass'
     stop
  endif
  masse = 0.0

  allocate(densite_pouss(n_grains_tot,Nc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error densite_pouss'
     stop
  endif
  densite_pouss = 0.0

  allocate(densite_gaz(Nc), densite_gaz_midplane(n_rad), masse_gaz(Nc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error densite_gaz'
     stop
  endif
  densite_gaz = 0.0 ; densite_gaz_midplane = 0.0 ; masse_gaz = 0.0

end subroutine allocate_densities

!*****************************************************************

subroutine deallocate_densities

  deallocate(masse, densite_pouss, densite_gaz, densite_gaz_midplane, masse_gaz)

  return

end subroutine deallocate_densities

!*****************************************************************

subroutine alloc_dust_prop()

  integer ::  alloc_status

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

  return

end subroutine alloc_dust_prop

!*****************************************************************

subroutine alloc_dynamique(n_cells_max)
  ! Alloue les tableaux dynamiquement en fonction de l'organisation
  ! des calculs, de la presence de stratification ...
  ! Permet d'optimiser la taille mémoire
  ! C. Pinte
  ! 12/05/05

  use stars, only : allocate_stellar_spectra
  use thermal_emission, only : allocate_temperature, allocate_thermal_emission, &
       allocate_weight_proba_emission, allocate_thermal_energy

  integer, intent(in), optional :: n_cells_max

  integer ::  alloc_status, Nc, p_Nc
  real :: mem_size

  if (present(n_cells_max)) then
     if (n_cells_max < n_cells) then
        write(*,*) "ERROR in alloc_dynamique : n_cells_max must be larger than n_cells"
        write(*,*) n_cells_max, n_cells
        write(*,*) "Exiting"
        stop
     endif

     Nc = n_cells_max
     if (p_n_cells == n_cells) then
        p_Nc = Nc
     else
        p_Nc = 1
     endif
  else
     Nc = n_cells
     p_Nc = p_n_cells
  endif

  allocate(stream(nb_proc), gauss_random_saved(nb_proc), lgauss_random_saved(nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error random number stream'
     stop
  endif
  gauss_random_saved = 0.0_dp
  lgauss_random_saved = .false.

  allocate(n_phot_envoyes(n_lambda,nb_proc),  stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error n_phot_envoyes'
     stop
  endif
  n_phot_envoyes = 0.0

  if (lSigma_file) then
     allocate(Surface_density(n_rad), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Sigma'
        stop
     endif
     Surface_density = 0.0_dp
  endif

  allocate(l_dark_zone(Nc), ri_in_dark_zone(n_az), ri_out_dark_zone(n_az),&
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

  ! **************************************************
  ! Tableaux relatifs aux prop en fct de lambda
  ! **************************************************
  call allocate_stellar_spectra(n_lambda)

  call allocate_thermal_energy(Nc)

  ! Tableaux relatifs aux prop optiques des cellules
  allocate(kappa(Nc,n_lambda), kappa_abs_LTE(Nc,n_lambda), stat=alloc_status)

  if (.not.lonly_LTE .or. .not.lonly_nLTE) then
     allocate(proba_abs_RE_LTE(Nc,n_lambda),  stat=alloc_status)
     if (lRE_nLTE)  allocate(kappa_abs_nLTE(Nc,n_lambda), stat=alloc_status)
     if (lRE_nLTE.or.lnRE) allocate(proba_abs_RE_LTE_p_nLTE(Nc,n_lambda), stat=alloc_status)
     if (lnRE) allocate(proba_abs_RE(Nc,n_lambda), kappa_abs_RE(Nc,n_lambda), stat=alloc_status)
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

  ! todo : could be p_Nc ...
  allocate(tab_albedo_pos(Nc,n_lambda),stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos, tab_albedo_pos'
     stop
  endif
  tab_albedo_pos = 0

  if (aniso_method==2) then
     allocate(tab_g_pos(Nc,n_lambda),stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_albedo_pos, tab_g_pos'
        stop
     endif
     tab_g_pos = 0.0
  endif

  ! **************************************************
  ! Tableaux relatifs aux prop optiques des cellules ou des grains
  ! **************************************************
  if (scattering_method == 2) then ! prop par cellule
     if (lsepar_pola) then
        mem_size = (5. * nang_scatt) * p_Nc * p_n_lambda_pos * 4. / 1024.**3
     else
        mem_size = (2. * nang_scatt) * p_Nc * p_n_lambda_pos * 4. / 1024.**3
     endif
     if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering matrices"

     allocate(tab_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_s11_pos'
        stop
     endif
     tab_s11_pos = 0

     allocate(prob_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error prob_s11_pos'
        stop
     endif
     prob_s11_pos = 0

     if (lsepar_pola) then
        allocate(tab_s12_o_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), &
             tab_s33_o_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), &
             tab_s34_o_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), &
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
     mem_size = (1.0 * n_grains_tot) * p_Nc * n_lambda * 4. / 1024.**3
     if (mem_size < max_mem) then
        low_mem_scattering = .false.
        if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"

        allocate(ksca_CDF(0:n_grains_tot,p_Nc,n_lambda), stat=alloc_status)

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
  ! Tableaux relatifs aux prop d'emission des cellules
  ! **************************************************
  allocate(l_emission_pah(0:n_rad+1,0:nz+1), stat=alloc_status) ! OUTDATED
  if (alloc_status > 0) then
     write(*,*) 'Allocation error l_emission_pah'
     stop
  endif
  l_emission_pah = .false.

  if (lweight_emission) then
     call allocate_weight_proba_emission(Nc)
  endif


  if (lorigine) then
     allocate(disk_origin(n_lambda, Nc, nb_proc), star_origin(n_lambda, nb_proc), stat=alloc_status)
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
  call allocate_temperature(Nc)

  ! **************************************************
  ! Tableaux relatifs au *calcul* de la temperature
  ! **************************************************
  if (lTemp) call allocate_thermal_emission(Nc, p_Nc)

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
     allocate(STOKE_io(npix_x,npix_y,capt_debut:capt_fin,N_phi,N_type_flux), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error STOKE_io'
        stop
     endif
     STOKE_io = 0.0

     allocate(STOKEI(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error STOKEI'
        stop
     endif
     STOKEI = 0.0

     if (lsepar_pola) then
        allocate(STOKEQ(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEQ'
           stop
        endif
        STOKEQ = 0.0

        allocate(STOKEU(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEU'
           stop
        endif
        STOKEU=0.0

        allocate(STOKEV(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEV'
           stop
        endif
        STOKEV = 0.0
     endif

     if (lsepar_contrib) then
        allocate(STOKEI_star(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_star'
           stop
        endif
        STOKEI_star = 0.0

        allocate(STOKEI_star_scat(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_star_scat'
           stop
        endif
        STOKEI_star_scat = 0.0

        allocate(STOKEI_disk(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error STOKEI_disk'
           stop
        endif
        STOKEI_disk = 0.0

        allocate(STOKEI_disk_scat(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
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
     allocate(tab_abundance(Nc), Tcin(Nc), lcompute_molRT(Nc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error Tcin & tab_abundance'
        stop
     endif
     tab_abundance = 0.0
     lcompute_molRT = .true.
     Tcin=0.0

     allocate(vfield(Nc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error vfield'
        stop
     endif
     vfield=0.0 !; vx=0.0 ; vy=0.0

     allocate(v_turb(Nc), v_line(Nc), deltaVmax(Nc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sigma2'
        stop
     endif
     v_turb = 0.0 ; v_line = 0.0 ;   deltaVmax = 0.0

     allocate(tab_dnu_o_freq(Nc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_dnu'
        stop
     endif
     tab_dnu_o_freq=0.0

     allocate(norme_phiProf_m1(Nc), sigma2_phiProf_m1(Nc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error norme_phiProf_m1'
        stop
     endif
     norme_phiProf_m1 = 0.0 ; sigma2_phiProf_m1 = 0.0

  endif ! lemission_mol

  return

end subroutine alloc_dynamique

!**********************************************************************

subroutine dealloc_em_th()

  use thermal_emission, only : deallocate_thermal_emission
  use stars, only : deallocate_stellar_spectra

  deallocate(n_phot_envoyes)

  deallocate(tab_albedo,C_ext,C_sca,C_abs,C_abs_norm,tab_g) ! q_geo

  !deallocate(E_stars,E_disk,frac_E_stars,E_totale)
  call deallocate_stellar_spectra()

  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2,tab_amu1_coating,tab_amu2_coating)

  deallocate(kappa,kappa_abs_LTE)
  if (allocated(proba_abs_RE_LTE)) then
     deallocate(proba_abs_RE_LTE)
     if (lRE_nLTE) deallocate(kappa_abs_nLTE)
     if (lRE_nLTE.or.lnRE) deallocate(proba_abs_RE_LTE_p_nLTE)
     if (lnRE) deallocate(proba_abs_RE,kappa_abs_RE)
  endif

  deallocate(tab_albedo_pos)
  if (allocated(tab_g_pos)) deallocate(tab_g_pos)

  if (scattering_method == 2) then ! prop par cellule
     deallocate(tab_s11_pos,prob_s11_pos)
     if (lsepar_pola) deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
  else ! prop par grains
     deallocate(ksca_CDF)
  endif ! method

  deallocate(tab_s11,tab_s12,tab_s33,tab_s34,prob_s11)

  deallocate(l_emission_pah) ! OUTDATED
  if (lorigine) deallocate(disk_origin,star_origin)

  if (lTemp) call deallocate_thermal_emission()

  if (lTemp.or.lsed) then
     deallocate(sed,sed_q,sed_u,sed_v)
     deallocate(sed_star,sed_star_scat,sed_disk,sed_disk_scat,n_phot_sed)
     deallocate(sed1_io,sed2_io,wave2_io)
  endif ! ltemp.or.lSED


  return

end subroutine dealloc_em_th

!******************************************************************************

subroutine realloc_dust_mol()

  use stars, only : allocate_stellar_spectra

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
  allocate(kappa(n_cells,n_lambda),kappa_abs_LTE(n_cells,n_lambda), &
       emissivite_dust(n_cells,n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error emissivite_dust (realloc)'
     stop
  endif
  kappa = 0.0 ; kappa_abs_LTE = 0.0 ; emissivite_dust = 0.0

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
  allocate(tab_albedo_pos(n_cells,n_lambda),stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos (realloc)'
     stop
  endif
  tab_albedo_pos = 0

  if (aniso_method==2) then
     allocate(tab_g_pos(n_cells,n_lambda),stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_g_pos (realloc)'
        stop
     endif
     tab_g_pos = 0.0
  endif

  call allocate_stellar_spectra(n_lambda)

  return

end subroutine realloc_dust_mol

!******************************************************************************

subroutine clean_mem_dust_mol()

  ! Ne reste que tab_lambda, tab_delta_lambda, tab_lambda_inf, tab_lambda_sup, kappa, emissivite_dust
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
  deallocate(tab_albedo_pos)
  if (allocated(tab_g_pos)) deallocate(tab_g_pos)

  return

end subroutine clean_mem_dust_mol

!******************************************************************************

subroutine realloc_step2()

  use radiation_field, only : allocate_radiation_field_step2
  use stars, only : allocate_stellar_spectra, deallocate_stellar_spectra
  use thermal_emission, only : deallocate_temperature_calculation, realloc_emitting_fractions

  integer :: alloc_status, mem_size, p_n_lambda2_pos

  n_lambda = n_lambda2
  if (ldust_prop) then
     p_n_lambda2_pos = n_lambda2
  else
     p_n_lambda2_pos = 1
  endif
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
  if (lTemp) then

     call deallocate_temperature_calculation()
  endif

  call allocate_radiation_field_step2()

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

  call deallocate_stellar_spectra()
  call allocate_stellar_spectra(n_lambda2)

  call realloc_emitting_fractions()

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

  deallocate(tab_albedo_pos)
  ! todo : could be p_n_cells
  allocate(tab_albedo_pos(n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_albedo_pos'
     stop
  endif
  tab_albedo_pos = 0

  if (allocated(tab_g_pos)) then
     deallocate(tab_g_pos)
     ! todo : could be p_n_cells
     allocate(tab_g_pos(n_cells,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tab_g_pos'
        stop
     endif
     tab_g_pos = 0
  endif

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

  deallocate(kappa,kappa_abs_LTE)
  deallocate(proba_abs_RE_LTE)
  if (lRE_nLTE) deallocate(kappa_abs_nLTE)
  if (lRE_nLTE.or.lnRE) deallocate(proba_abs_RE_LTE_p_nLTE)
  if (lnRE) deallocate(proba_abs_RE,kappa_abs_RE)
  allocate(kappa(n_cells,n_lambda2), kappa_abs_LTE(n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa'
     stop
  endif
  kappa=0.0 ; kappa_abs_LTE=0.0
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
     allocate(disk_origin(n_lambda2, n_cells, nb_proc), star_origin(n_lambda2, nb_proc), stat=alloc_status)
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
   Jmol = 0.0_dp

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
     Jmol2 = 0.0_dp
  endif

  ! Methode d'echantillonnage
  if (npix_x_save > 1) then
     RT_line_method = 2 ! creation d'une carte avec pixels carres
     npix_x = npix_x_save ; npix_y = npix_y_save ! we update the value after the SED calculation

     write(*,*) "WARNING : memory size if lots of pixels"
     allocate(spectre(npix_x,npix_y,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_incl,RT_n_az), &
          continu(npix_x,npix_y,nTrans_raytracing,RT_n_incl,RT_n_az), stars_map(npix_x,npix_y,1), stat=alloc_status)
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

  use stars, only : deallocate_stellar_spectra

  ! Dealloue ce qui n'a pas ete libere par  clean_mem_dust_mol
  deallocate(tab_lambda, tab_delta_lambda, tab_lambda_inf, tab_lambda_sup)
  deallocate(kappa,emissivite_dust)

  call deallocate_stellar_spectra()

  deallocate(Level_energy,poids_stat_g,j_qnb,Aul,fAul,Bul,fBul,Blu,fBlu,transfreq, &
       itransUpper,itransLower,nCollTrans,nCollTemps,collTemps,collBetween, &
       iCollUpper,iCollLower)

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
