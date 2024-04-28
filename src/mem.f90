module mem

  use parametres
  use grains
  use dust_prop
  use naleat
  use molecular_emission
  use utils
  use messages
  use wavelengths
  use output
  use density
  use cylindrical_grid

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
  if (alloc_status > 0) call error('Allocation error mass')
  masse = 0.0

  allocate(densite_pouss(n_grains_tot,Nc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error densite_pouss')
  densite_pouss = 0.0

  allocate(densite_gaz(Nc), densite_gaz_midplane(n_rad), masse_gaz(Nc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error densite_gaz')
  densite_gaz = 0.0 ; densite_gaz_midplane = 0.0 ; masse_gaz = 0.0

end subroutine allocate_densities

!*****************************************************************

subroutine deallocate_densities

  if (allocated(masse)) deallocate(masse,densite_pouss,densite_gaz,densite_gaz_midplane,masse_gaz)

  return

end subroutine deallocate_densities

!*****************************************************************

subroutine alloc_dust_prop()
 ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
 ! 20/04/2023

  integer ::  alloc_status

  allocate(tab_albedo(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo')
  tab_albedo = 0

  allocate(C_ext(n_grains_tot,n_lambda), C_sca(n_grains_tot,n_lambda), &
       C_abs(n_grains_tot,n_lambda), C_abs_norm(n_grains_tot,n_lambda), stat=alloc_status) !  C_geo(n_grains_tot)
  if (alloc_status > 0) call error('Allocation error q_ext')
  C_ext = 0 ; C_sca = 0 ; C_abs = 0 ; C_abs_norm =0


  allocate(tab_g(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_g')
  tab_g = 0

  ! **************************************************
  ! tableaux relatifs aux prop optiques des grains
  ! **************************************************
  if (laggregate.or.lmueller) then
     allocate(tab_mueller(4,4,0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_mueller')
     tab_mueller = 0
  else
     allocate(tab_s11(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s11')
     tab_s11 = 0

     allocate(tab_s12(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s12')
     tab_s12 = 0

     allocate(tab_s33(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s33')
     tab_s33 = 0

     allocate(tab_s34(0:nang_scatt,n_grains_tot,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s34')
     tab_s34 = 0
  endif

  allocate(prob_s11(n_lambda,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error prob_s11')
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
  ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
  ! 20/04/2023

  use stars, only : allocate_stellar_spectra
  use thermal_emission, only : allocate_temperature, allocate_thermal_emission, &
       allocate_weight_proba_emission, allocate_thermal_energy
  use radiation_field, only :  allocate_radiation_field_step1
  use output, only : allocate_origin

  integer, intent(in), optional :: n_cells_max

  integer ::  alloc_status, Nc, p_Nc
  real :: mem_size

  if (present(n_cells_max)) then
     if (n_cells_max < n_cells) then
        write(*,*) "n_cells_max=", n_cells_max, "<", n_cells
        call error("n_cells_max must be larger than n_cells")
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

  allocate(stream(nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error random number stream')

  allocate(n_phot_envoyes(n_lambda,nb_proc),  stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error n_phot_envoyes')
  n_phot_envoyes = 0.0

  allocate(l_dark_zone(Nc), ri_in_dark_zone(n_az), ri_out_dark_zone(n_az),&
          zj_sup_dark_zone(n_rad,n_az), zj_inf_dark_zone(n_rad,n_az), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error l_dark_zone')
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
  allocate(kappa(p_Nc,n_lambda), kappa_abs_LTE(p_Nc,n_lambda), tab_albedo_pos(p_Nc,n_lambda), kappa_factor(Nc), stat=alloc_status)
  kappa=0.0 ; kappa_abs_LTE=0.0 ; tab_albedo_pos = 0
  if (alloc_status > 0) then
     write(*,*) 'Allocation error kappa and albedo'
     stop
  endif
  if (lRE_nLTE)  then
     allocate(kappa_abs_nLTE(p_Nc,n_lambda), stat=alloc_status)
     kappa_abs_nLTE=0.0
  endif

  ! the following arrays are updated live per cell, so we cannot use a pointer
  if (.not.(lonly_LTE.or.lonly_nLTE)) then
     allocate(proba_abs_RE_LTE(Nc,n_lambda),  stat=alloc_status)
     proba_abs_RE_LTE=0.0
  endif
  if (lRE_nLTE.or.lnRE) then
     allocate(proba_abs_RE_LTE_p_nLTE(Nc,n_lambda), stat=alloc_status)
     proba_abs_RE_LTE_p_nLTE=0.0
  endif
  if (lnRE) then
     allocate(kappa_abs_RE(Nc,n_lambda), proba_abs_RE(Nc,n_lambda), stat=alloc_status)
     kappa_abs_RE=0.0 ; proba_abs_RE=0.0
  endif
  if (alloc_status > 0) call error('Allocation error kappa_abs')

  if (aniso_method==2) then
     allocate(tab_g_pos(p_Nc,n_lambda),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_albedo_pos, tab_g_pos')
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

     allocate(prob_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error prob_s11_pos')
     prob_s11_pos = 0

     if (lmueller) then
     	allocate(tab_mueller_pos(4,4,0:nang_scatt, p_Nc, p_n_lambda_pos), stat=alloc_status)
     	if (alloc_status > 0) call error('Allocation error tab_mueller_pos')
     	tab_mueller_pos = 0
     else
	allocate(tab_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), stat=alloc_status)
	if (alloc_status > 0) call error('Allocation error tab_s11_pos')
	tab_s11_pos = 0

	if (lsepar_pola) then
	   allocate(tab_s12_o_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), &
		    tab_s33_o_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), &
		    tab_s34_o_s11_pos(0:nang_scatt, p_Nc, p_n_lambda_pos), &
		    stat=alloc_status)
	   if (alloc_status > 0) call error('Allocation error tab_s12_o_s11_pos')
	   tab_s12_o_s11_pos = 0
	   tab_s33_o_s11_pos = 0
	   tab_s34_o_s11_pos = 0
	endif
     endif
  else ! scattering method==1 --> prop par grains
     mem_size = (1.0 * n_grains_tot) * p_Nc * n_lambda * 4. / 1024.**3
     if (mem_size < max_mem) then
        low_mem_scattering = .false.
        if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"

        allocate(ksca_CDF(0:n_grains_tot,p_Nc,n_lambda), stat=alloc_status)

        if (alloc_status > 0) call error('Allocation error ksca_CDF')
        ksca_CDF = 0
     else ! Array is to big, we will recompute ksca_CDF on the fly
        low_mem_scattering = .true.
        write(*,*) "Using low memory mode for scattering properties"
     endif
  endif ! method

  ! **************************************************
  ! Tableaux relatifs aux prop d'emission des cellules
  ! **************************************************
  if (lweight_emission) call allocate_weight_proba_emission(Nc)
  if (lorigine) call allocate_origin()

  ! **************************************************
  ! Tableaux de temperature
  ! **************************************************
  call allocate_temperature(Nc)

  ! **************************************************
  ! Tableaux relatifs au *calcul* de la temperature
  ! **************************************************
  if (lTemp) call allocate_thermal_emission(Nc, p_Nc)

  call allocate_radiation_field_step1(Nc)

  ! **************************************************
  ! Tableaux relatifs aux SEDs
  ! **************************************************
  if (lTemp.or.lsed) call allocate_sed()

  ! **************************************************
  ! Tableaux relatifs aux images
  ! **************************************************
  if (lmono0.and.loutput_mc) call allocate_mc_images()

  ! **************************************************
  ! Tableaux relatifs a l'emission moleculaire
  ! **************************************************
  if (lemission_mol) then
     allocate(tab_abundance(Nc), Tcin(Nc), lcompute_molRT(Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Tcin & tab_abundance')
     tab_abundance = 0.0
     lcompute_molRT = .true.
     Tcin=0.0

     allocate(vfield(Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error vfield')
     vfield=0.0 !; vx=0.0 ; vy=0.0

     allocate(v_turb(Nc), v_line(Nc), deltaVmax(Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error sigma2')
     v_turb = 0.0 ; v_line = 0.0 ;   deltaVmax = 0.0

     allocate(tab_dnu_o_freq(Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_dnu')
     tab_dnu_o_freq=0.0

     allocate(norme_phiProf_m1(Nc), sigma2_phiProf_m1(Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error norme_phiProf_m1')
     norme_phiProf_m1 = 0.0 ; sigma2_phiProf_m1 = 0.0

  endif ! lemission_mol

  return

end subroutine alloc_dynamique

!**********************************************************************

subroutine deallocate_em_th_mol()
 ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
 ! 20/04/2023

  use thermal_emission, only : deallocate_thermal_emission
  use stars, only : deallocate_stellar_spectra
   use output, only : deallocate_origin

  deallocate(n_phot_envoyes)

  deallocate(tab_albedo,C_ext,C_sca,C_abs,C_abs_norm,tab_g) ! q_geo

  !deallocate(E_stars,E_disk,frac_E_stars,E_totale)
  call deallocate_stellar_spectra()

  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2,tab_amu1_coating,tab_amu2_coating)

  deallocate(kappa,kappa_abs_LTE,kappa_factor)
  if (allocated(proba_abs_RE_LTE)) then
     deallocate(proba_abs_RE_LTE)
  endif
  if (lRE_nLTE) deallocate(kappa_abs_nLTE)
  if (allocated(proba_abs_RE_LTE_p_nLTE)) deallocate(proba_abs_RE_LTE_p_nLTE)
  if (allocated(proba_abs_RE)) deallocate(proba_abs_RE)
  if (allocated(kappa_abs_RE)) deallocate(kappa_abs_RE)

  deallocate(tab_albedo_pos)
  if (allocated(tab_g_pos)) deallocate(tab_g_pos)

  if (scattering_method == 2) then ! prop par cellule
     if (lmueller) then
        deallocate(tab_mueller_pos,prob_s11_pos)
     else
        deallocate(tab_s11_pos,prob_s11_pos)
        if (lsepar_pola) deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
     endif
  else ! prop par grains
     if (allocated(ksca_CDF)) deallocate(ksca_CDF)
  endif ! method

  if (lmueller.or.laggregate) then
     deallocate (tab_mueller, prob_s11)
  else
     deallocate(tab_s11,tab_s12,tab_s33,tab_s34,prob_s11)
  endif

  if (lorigine) call deallocate_origin()

  if (lTemp) call deallocate_thermal_emission()

  if (lTemp.or.lsed) call deallocate_sed()

  return

end subroutine deallocate_em_th_mol

!******************************************************************************

subroutine realloc_dust_mol(imol)

  use stars, only : allocate_stellar_spectra

  integer, intent(in) :: imol

  integer :: alloc_status, iTrans_min, iTrans_max

  iTrans_min = mol(imol)%iTrans_min
  iTrans_max = mol(imol)%iTrans_max

  allocate(tab_lambda(n_lambda), tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), &
       tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), &
       tab_amu1_coating(n_lambda, n_pop), tab_amu2_coating(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_lambda (realloc)')
  tab_lambda=0.0 ; tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda= 0.0 ;
  tab_amu1=0.0 ; tab_amu2=0.0 ; tab_amu1_coating=0.0 ; tab_amu2_coating=0.0

  allocate(tab_albedo(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo (realloc)')
  tab_albedo = 0

  allocate(C_ext(n_grains_tot,n_lambda), C_sca(n_grains_tot,n_lambda), &
       C_abs(n_grains_tot,n_lambda),  C_abs_norm(n_grains_tot,n_lambda), tab_g(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error C_ext (realloc)')
  C_ext = 0 ; C_sca = 0 ; C_abs = 0 ; C_abs_norm = 0 ; tab_g = 0

  ! Tableaux relatifs aux prop optiques des cellules
  allocate(kappa(p_n_cells,iTrans_min:iTrans_max),kappa_abs_LTE(p_n_cells,iTrans_min:iTrans_max), kappa_factor(n_cells), &
       emissivite_dust(n_cells,iTrans_min:iTrans_max), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error emissivite_dust (realloc)')
  kappa = 0.0 ; kappa_abs_LTE = 0.0 ; emissivite_dust = 0.0

  if (lRE_nLTE) then
     allocate(kappa_abs_nLTE(p_n_cells,iTrans_min:iTrans_max), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kappa_abs_nLTE (realloc)')
     kappa_abs_nLTE = 0.0
  endif

  allocate(tab_albedo_pos(p_n_cells,iTrans_min:iTrans_max),stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo_pos (realloc)')
  tab_albedo_pos = 0

  if (aniso_method==2) then
     allocate(tab_g_pos(p_n_cells,iTrans_min:iTrans_max),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_g_pos (realloc)')
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
  !deallocate(kappa_abs_LTE)
  deallocate(tab_albedo_pos)
  if (allocated(tab_g_pos)) deallocate(tab_g_pos)

  return

end subroutine clean_mem_dust_mol

!******************************************************************************

subroutine realloc_dust_atom()
  ! This routine should be the mirror of the mol one, except for the tab_lambda which is allocated
  ! when the gas atom RT grid is defined (from reading the atomic species).
  ! Still, the test on the allocation and the call of clean_mem_dust_mol in init_dust_atom means
  ! that they are not deallocated after temperature calculation. What am I missing ?

  !   use stars, only : allocate_stellar_spectra !-> not use yet in atom transfer.

  integer :: alloc_status

  if (lvariable_dust) then
     write(*,*) " WARNING: sizes of dust transfer could be very big !"
     !TO DO: better storing of quantities / recuction of n_lambda
  endif

  !Note: tab_lambda(n_lambda) is already allocated in atomic_transfer
  !	the tab_lambda_* or tab_delta_lambda should be de-allocated when tab_lambda is allocated in atom_rt.
  allocate(tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), &
       tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), &
       tab_amu1_coating(n_lambda, n_pop), tab_amu2_coating(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_lambda (realloc)')
!   tab_lambda=0.0
  tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda= 0.0 ;
  tab_amu1=0.0 ; tab_amu2=0.0 ; tab_amu1_coating=0.0 ; tab_amu2_coating=0.0

  allocate(tab_albedo(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo (realloc)')
  tab_albedo = 0

  allocate(C_ext(n_grains_tot,n_lambda), C_sca(n_grains_tot,n_lambda), &
       C_abs(n_grains_tot,n_lambda),  C_abs_norm(n_grains_tot,n_lambda), tab_g(n_grains_tot,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error C_ext (realloc)')
  C_ext = 0 ; C_sca = 0 ; C_abs = 0 ; C_abs_norm = 0 ; tab_g = 0

  ! Tableaux relatifs aux prop optiques des cellules
  if (allocated(kappa)) deallocate(kappa)
  if (allocated(kappa_abs_LTE)) deallocate(kappa_abs_LTE)
  if (allocated(kappa_factor)) deallocate(kappa_factor)
  allocate(kappa(p_n_cells,n_lambda),kappa_abs_LTE(p_n_cells,n_lambda), kappa_factor(n_cells), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error emissivite_dust (realloc atom)')
  kappa = 0.0 ; kappa_abs_LTE = 0.0
  !mind the shape of the array compared to the others.
  if (allocated(emissivite_dust)) deallocate(emissivite_dust)
  allocate(emissivite_dust(n_lambda,n_cells),stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error emissivite_dust (realloc atom)')
  if (lvariable_dust.or.(sizeof(emissivite_dust)/1024.**3 > 5)) then
     write(*,*) "  *** WARNING: using ", sizeof(emissivite_dust)/1024.**3, " GB for emissivite_dust"
  endif
  emissivite_dust = 0.0

  if (lRE_nLTE) then
     if (allocated(kappa_abs_nlte)) deallocate(kappa_abs_nlte)
     allocate(kappa_abs_nLTE(p_n_cells,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kappa_abs_nLTE (realloc atom)')
     kappa_abs_nLTE = 0.0
  endif

  allocate(tab_albedo_pos(p_n_cells,n_lambda),stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo_pos (realloc atom)')
  tab_albedo_pos = 0

  if (aniso_method==2) then
     allocate(tab_g_pos(p_n_cells,n_lambda),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_g_pos (realloc atom)')
     tab_g_pos = 0.0
  endif

!   call allocate_stellar_spectra(n_lambda)

  return

end subroutine realloc_dust_atom

!******************************************************************************

subroutine realloc_step2()
! Ajout du cas ou les matrices de Mueller sont donnees en entrees
! 20/04/2023

  use dust_ray_tracing, only : select_scattering_method
  use radiation_field, only : allocate_radiation_field_step2
  use stars, only : allocate_stellar_spectra, deallocate_stellar_spectra
  use thermal_emission, only : deallocate_temperature_calculation, realloc_emitting_fractions
  use output, only : allocate_origin

  integer :: alloc_status, mem_size, p_n_lambda2_pos

  n_lambda = n_lambda2
  if (ldust_prop) then
     p_n_lambda2_pos = n_lambda2
  else
     p_n_lambda2_pos = 1
  endif
  p_n_lambda_pos = p_n_lambda2_pos ! just in case

  ! parametrage methode de diffusion
  call select_scattering_method(p_n_cells)

  ! Liberation memoire
  if (lTemp) call deallocate_temperature_calculation()

  call allocate_radiation_field_step2()

  ! Liberation memoire step1 et reallocation step 2
  deallocate(tab_lambda,tab_lambda_inf,tab_lambda_sup,tab_delta_lambda,tab_amu1,tab_amu2,tab_amu1_coating,tab_amu2_coating)
  allocate(tab_lambda(n_lambda2),tab_lambda_inf(n_lambda2),tab_lambda_sup(n_lambda2),tab_delta_lambda(n_lambda2),&
       tab_amu1(n_lambda2,n_pop),tab_amu2(n_lambda2,n_pop), &
       tab_amu1_coating(n_lambda2,n_pop),tab_amu2_coating(n_lambda2,n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_lambda in init_lambda2')
  tab_lambda=0.0 ; tab_lambda_inf = 0.0 ; tab_lambda_sup = 0.0 ; tab_delta_lambda=0.0
  tab_amu1=0.0 ; tab_amu2=0.0 ;  tab_amu1_coating=0.0 ; tab_amu2_coating=0.0

  deallocate(n_phot_envoyes)
  allocate(n_phot_envoyes(n_lambda2,nb_proc),  stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error n_phot_envoyes')
  n_phot_envoyes = 0.0

  call deallocate_sed()
  call allocate_sed()

  call deallocate_stellar_spectra()
  call allocate_stellar_spectra(n_lambda2)

  call realloc_emitting_fractions()

  deallocate(tab_albedo)
  allocate(tab_albedo(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo')
  tab_albedo = 0

  deallocate(C_ext)
  allocate(C_ext(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error C_ext')
  C_ext = 0

  deallocate(C_sca)
  allocate(C_sca(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error C_sca')
  C_sca = 0

  deallocate(C_abs,C_abs_norm)
  allocate(C_abs(n_grains_tot,n_lambda2), C_abs_norm(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error C_abs')
  C_abs = 0 ; C_abs_norm = 0

  deallocate(tab_g)
  allocate(tab_g(n_grains_tot,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_g')
  tab_g = 0

  deallocate(tab_albedo_pos)
  ! todo : could be p_n_cells
  allocate(tab_albedo_pos(n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_albedo_pos')
  tab_albedo_pos = 0

  if (allocated(tab_g_pos)) then
     deallocate(tab_g_pos)
     ! todo : could be p_n_cells
     allocate(tab_g_pos(n_cells,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_g_pos')
     tab_g_pos = 0
  endif

  if (laggregate.or.lmueller) then
     deallocate(tab_mueller)
     allocate(tab_mueller(4,4,0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_mueller')
     tab_mueller = 0
  else
     deallocate(tab_s11)
     allocate(tab_s11(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s11')
     tab_s11 = 0

     deallocate(tab_s12)
     allocate(tab_s12(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s12')
     tab_s12 = 0

     deallocate(tab_s33)
     allocate(tab_s33(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s33')
     tab_s33 = 0

     deallocate(tab_s34)
     allocate(tab_s34(0:nang_scatt,n_grains_tot,n_lambda2), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s34')
     tab_s34 = 0
  endif

  deallocate(prob_s11)
  allocate(prob_s11(n_lambda2,n_grains_tot,0:nang_scatt), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error prob_s11')
  prob_s11 = 0

  if (scattering_method == 2) then
     if (lmueller) then
        if (allocated(tab_mueller_pos)) deallocate(tab_mueller_pos)
        allocate(tab_mueller_pos(4,4,0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_mueller_pos')
        tab_mueller_pos = 0
     else
        if (allocated(tab_s11_pos)) deallocate(tab_s11_pos)
        allocate(tab_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_s11_pos')

        if (lsepar_pola) then
           if (allocated(tab_s12_o_s11_pos)) deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
           allocate(tab_s12_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
                tab_s33_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
                tab_s34_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
                stat=alloc_status)
           if (alloc_status > 0) call error('Allocation error tab_s12_o_s11_pos')
           tab_s12_o_s11_pos = 0
           tab_s33_o_s11_pos = 0
           tab_s34_o_s11_pos = 0
        endif
     endif

     if (allocated(prob_s11_pos)) deallocate(prob_s11_pos)
     allocate(prob_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error prob_s11_pos')
     prob_s11_pos = 0
  else
     if (allocated(ksca_CDF)) deallocate(ksca_CDF)
     low_mem_scattering = .false.
     mem_size = n_grains_tot * p_n_cells * p_n_lambda2_pos * 4. / 1024.**3
     if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"
     allocate(ksca_CDF(0:n_grains_tot,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error ksca_CDF')
     ksca_CDF = 0
  endif ! method

  deallocate(kappa,kappa_abs_LTE,kappa_factor)
  if (allocated(proba_abs_RE_LTE)) deallocate(proba_abs_RE_LTE)
  if (lRE_nLTE) deallocate(kappa_abs_nLTE)
  if (lRE_nLTE.or.lnRE) deallocate(proba_abs_RE_LTE_p_nLTE)
  if (lnRE) deallocate(proba_abs_RE,kappa_abs_RE)
  allocate(kappa(p_n_cells,n_lambda2), kappa_abs_LTE(p_n_cells,n_lambda2), kappa_factor(n_cells), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error kappa')
  kappa=0.0 ; kappa_abs_LTE=0.0
  if (lRE_nLTE) then
     allocate(kappa_abs_nLTE(p_n_cells,n_lambda2))
     if (alloc_status > 0) call error('Allocation error kappa_abs_nLTE')
     kappa_abs_nLTE=0.0
  endif

  if (lorigine) call allocate_origin()

  return

end subroutine realloc_step2

!**********************************************************************

subroutine realloc_ray_tracing_scattering_matrix()
! Ajout du cas ou les matrices de Mueller sont donnees en entrees
! 20/04/2023

  integer, parameter :: p_n_lambda2_pos = 1
  integer :: alloc_status

  ! parametrage methode de diffusion
  scattering_method = 2
  write(*,fmt='(" Using scattering method ",i1)') scattering_method
  lscattering_method1 = (scattering_method==1)

  if (lmueller) then
     if (allocated(tab_mueller_pos)) deallocate(tab_mueller_pos)
     allocate(tab_mueller_pos(4,4,0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_mueller_pos')
     tab_mueller_pos = 0
  else
     if (allocated(tab_s11_pos)) deallocate(tab_s11_pos)
     allocate(tab_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tab_s11_pos')
     tab_s11_pos = 0

     if (lsepar_pola) then
        if (allocated(tab_s12_o_s11_pos)) deallocate(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos)
        allocate(tab_s12_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
             tab_s33_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
             tab_s34_o_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), &
             stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_s12_o_s11_pos')
        tab_s12_o_s11_pos = 0
        tab_s33_o_s11_pos = 0
        tab_s34_o_s11_pos = 0
     endif
  endif

  if (allocated(prob_s11_pos)) deallocate(prob_s11_pos)
  allocate(prob_s11_pos(0:nang_scatt,p_n_cells,p_n_lambda2_pos), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error prob_s11_pos')
  prob_s11_pos = 0

  return

end subroutine realloc_ray_tracing_scattering_matrix

!**********************************************************************


subroutine alloc_emission_mol(imol)

  integer, intent(in) :: imol
  integer :: alloc_status, n_speed
  integer :: iTrans_min, iTrans_max, level_min, level_max

  alloc_status = 0

  n_speed = mol(imol)%n_speed_rt ! I use the same now

  if (lmol_LTE) then ! Reducing memory in LTE mode
     iTrans_min = minval(mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing))
     iTrans_max = maxval(mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing))

     level_min = iTransLower(iTrans_min)
     level_max = iTransUpper(iTrans_max)
  else ! Todo : we can reduce memory here too, up to level_max
     iTrans_min = 1
     iTrans_max = nTrans_tot

     level_min=1
     level_max=nLevels
  endif

  mol(imol)%iTrans_min = iTrans_min
  mol(imol)%iTrans_max = iTrans_max
  mol(imol)%level_min = level_min
  mol(imol)%level_max = level_max

  allocate(kappa_mol_o_freq(n_cells,iTrans_min:iTrans_max), emissivite_mol_o_freq(n_cells,iTrans_min:iTrans_max), &
       stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error kappa_mol_o_freq')
  kappa_mol_o_freq=0.0
  emissivite_mol_o_freq = 0.0

  allocate(tab_nLevel(n_cells,level_min:level_max), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_nLevel')
  tab_nLevel = 0.0

  ! Todo : we don't need this most of the time
  !allocate(maser_map(n_cells,nTrans_tot), stat=alloc_status)
  !if (alloc_status > 0) call error('Allocation error maser_map')
  !maser_map = 0.0

  allocate(tab_v(-n_largeur_Doppler*n_speed:n_largeur_Doppler*n_speed), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_v')
  tab_v=0.0

  ! We can skip that in LTE too
  allocate(tab_Cmb_mol(nTrans_tot), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_Cmb_mol')
  tab_Cmb_mol = 0.0

  if (.not.lmol_LTE) then
     allocate(Jmol(iTrans_min:iTrans_max,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Jmol')
     Jmol = 0.0_dp

     if (ldouble_RT) then ! Todo : this seems to be turned off at the moment
        allocate(kappa_mol_o_freq2(n_cells,nTrans_tot), emissivite_mol_o_freq2(n_cells,nTrans_tot),&
             stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error kappa_mol2')
        kappa_mol_o_freq2=0.0
        emissivite_mol_o_freq2=0.0

        allocate(tab_nLevel2(n_cells,nLevels), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_nLevel2')
        tab_nLevel2 = 0.0

        allocate(Jmol2(nTrans_tot,nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error Jmol2')
        Jmol2 = 0.0_dp
     endif
  endif

  call allocate_mol_maps(imol)

  if (ltau_surface.or.lflux_fraction_surface) then
     if (.not.allocated(tau_surface_map)) then
        allocate(tau_surface_map(npix_x,npix_y,RT_n_incl,RT_n_az,3,nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tau_surface')
        tau_surface_map = 0.0
     endif
  endif

  return

end subroutine alloc_emission_mol

!**********************************************************************

subroutine dealloc_emission_mol()

  use stars, only : deallocate_stellar_spectra

  ! Dealloue ce qui n'a pas ete libere par  clean_mem_dust_mol
  deallocate(tab_lambda, tab_delta_lambda, tab_lambda_inf, tab_lambda_sup)
  deallocate(kappa,kappa_abs_LTE,kappa_factor,emissivite_dust)

  call deallocate_stellar_spectra()

  deallocate(Level_energy,poids_stat_g,j_qnb,Aul,fAul,Bul,fBul,Blu,fBlu,transfreq, &
       itransUpper,itransLower,nCollTrans,nCollTemps,collTemps,collBetween, &
       iCollUpper,iCollLower)

  deallocate(kappa_mol_o_freq, emissivite_mol_o_freq, tab_nLevel, &
       tab_v, stars_map, tab_Cmb_mol)
  if (allocated(Jmol)) deallocate(Jmol)

  call deallocate_mol_maps()

  if (ldouble_RT) deallocate(kappa_mol_o_freq2, emissivite_mol_o_freq2, tab_nLevel2, Jmol2)

  deallocate(I0, I0c, tab_speed_rt) ! besoin de dealouer tab_speed_rt pour plusieurs rayons
  if (lorigine) deallocate(origine_mol)

  return

end subroutine dealloc_emission_mol

!**********************************************************************

subroutine sig_handler(sig)

  integer, intent(in) ::  sig

  select case(sig)
  case(2)
     write (*,*) 'mcfost : SIGINT Caught'
     call exit(1)
  case(15)
     write (*,*) 'mcfost : SIGTERM Caught'
     call exit(1)
  case default
     write (*,*) 'mcfost : signal caught :', sig
     read(*,*)
     return
  end select

end subroutine sig_handler

!*************************************************

end module mem
