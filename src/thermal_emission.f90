module thermal_emission

  use parametres
  use constantes
  use grains
  use dust_prop
  !$ use omp_lib
  use temperature
  use utils
  use PAH
  use grid
  use stars, only : spectre_etoiles, E_ISM, E_stars
  use input
  use density
  use cylindrical_grid

  implicit none
  save

  public :: prob_E_cell, E_totale, nbre_reemission, frac_E_stars, frac_E_disk, L_packet_th, E_abs_nRE, correct_E_emission, &
       DensE, DensE_m1, Dcoeff

  public :: allocate_temperature, deallocate_temperature_calculation, realloc_emitting_fractions, &
       allocate_thermal_emission, deallocate_thermal_emission, allocate_weight_proba_emission, &
       define_proba_weight_emission, emission_nre, im_reemission_lte, im_reemission_nlte, &
       im_reemission_qre, init_emissivite_nre, init_reemission, select_wl_em, &
       select_cellule, temp_finale, temp_finale_nlte, temp_nre, update_proba_abs_nre, &
       repartition_wl_em, repartition_energie, allocate_thermal_energy, set_min_temperature, &
       reset_temperature

  private

  ! Choix cellule d'emission pour cas monochromatique
  real(kind=dp), dimension(:,:), allocatable :: prob_E_cell !n_lambda,0:n_rad*nz
  real(kind=dp), dimension(:), allocatable :: frac_E_stars, frac_E_disk, E_totale !n_lambda

  ! fraction d'energie reemise sur energie etoile
  ! (Opacite moyenne de Planck * coeff)
  real(kind=dp), dimension(:,:), allocatable :: log_Qcool_minus_extra_heating ! 0:n_T, n_cells
  real(kind=dp), dimension(:,:), allocatable :: log_E_em_1grain  !n_grains,0:n_T
  real(kind=dp), dimension(:,:), allocatable :: E_em_1grain_nRE, log_E_em_1grain_nRE !n_grains,0:n_T

  ! Probabilite cumulee en lambda d'emissivite de la poussiere
  ! avec correction de temperature (dp/dT)
  ! (Bjorkman & Wood 2001, A&A 554-615 -- eq 9)
  real(kind=dp), dimension(:,:,:), allocatable :: kdB_dT_CDF ! 0:n_T,n_cells,n_lambda
  real(kind=dp), dimension(:,:,:), allocatable :: kdB_dT_1grain_LTE_CDF, kdB_dT_1grain_nLTE_CDF, kdB_dT_1grain_nRE_CDF ! n_grains,0:n_T,n_lambda

  integer, dimension(:,:), allocatable :: xT_ech ! n_cells, id
  integer, dimension(:,:,:), allocatable :: xT_ech_1grain, xT_ech_1grain_nRE ! n_grains, n_cells, id

  ! Proba cumulee en lambda d'emettre selon un corps noir
  real(kind=dp), dimension(:), allocatable :: spectre_emission_cumul !(0:n_lambda)

  ! emissivite en unite qq (manque une cst mais travail en relatif)
  real(kind=dp), dimension(:,:), allocatable :: Emissivite_nRE_old ! n_lambda, n_rad, nz, n_az

  real(kind=dp), dimension(:,:), allocatable :: nbre_reemission ! n_cells, id

  real(kind=dp) :: L_packet_th, E_abs_nRE

  ! Biais de l'emission vers la surface du disque
  real(kind=dp), dimension(:), allocatable :: weight_proba_emission, correct_E_emission

  real(kind=dp), dimension(:,:,:), allocatable :: DensE, DensE_m1 !n_rad, 0:n_z, n_az
  real(kind=dp), dimension(:,:,:), allocatable :: Dcoeff !n_rad, n_z, n_az



  contains

subroutine allocate_thermal_energy(Nc)

  integer, intent(in) :: Nc

  integer :: alloc_status

  allocate(E_disk(n_lambda), frac_E_stars(n_lambda), frac_E_disk(n_lambda), E_totale(n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error frac_E_stars')
  E_disk = 0.0 ; frac_E_stars = 0.0 ; frac_E_disk = 0.0 ; E_totale = 0.0

  allocate(prob_E_cell(0:Nc,n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error( 'Allocation error prob_E_cell 1')
  prob_E_cell = 0.0

  return

end subroutine allocate_thermal_energy

!***************************************************

subroutine allocate_thermal_emission(Nc,p_Nc)

  use radiation_field, only : allocate_radiation_field_step1

  integer, intent(in) :: Nc, p_Nc

  integer :: alloc_status
  real :: mem_size

  allocate(spectre_emission_cumul(0:n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error spectre_etoile')
  spectre_emission_cumul = 0.0

  allocate(tab_Temp(n_T), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_Temp')
  tab_Temp = 0.0

  allocate(log_Qcool_minus_extra_heating(n_T,p_Nc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error log_E_em')
  log_Qcool_minus_extra_heating = 0.0

  allocate(DensE(n_rad,0:nz,n_az), DensE_m1(n_rad,0:nz,n_az), Dcoeff(n_rad,0:nz,n_az), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error kappa_abs_1grain')
  DensE = 0.0 ; DensE_m1 = 0.0

  call allocate_radiation_field_step1(Nc)

  allocate(xT_ech(Nc,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error xT_ech')
  xT_ech = 2

  if (lreemission_stats) then
     allocate(nbre_reemission(Nc,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error nbre_reemission')
     nbre_reemission = 0.0
  endif

  mem_size = (1.0 * p_Nc) * n_T * n_lambda * 4. / 1024.**3
  if (mem_size < max_mem) then
     low_mem_th_emission = .false.
     if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for temperature calculation"
     allocate(kdB_dT_CDF(n_lambda,n_T,p_Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kdB_dT_CDF')
     kdB_dT_CDF = 0
  else
     low_mem_th_emission = .true.
     write(*,*) "Using low memory mode for thermal emission"
     allocate(kdB_dT_1grain_LTE_CDF(n_lambda,grain_RE_LTE_start:grain_RE_LTE_end,n_T), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kdB_dT_1grain_LTE_CDF')
     kdB_dT_1grain_LTE_CDF = 0
  endif

  if (lRE_nLTE) then
     mem_size = (1.0 * (grain_RE_nLTE_end-grain_RE_nLTE_start+2)) * p_Nc * n_lambda * 4. / 1024.**3
     if (mem_size < max_mem) then
        low_mem_th_emission_nLTE = .false.
        if (mem_size > 1) write(*,*) "Trying to allocate", mem_size, "GB for scattering probability"
        allocate(kabs_nLTE_CDF(grain_RE_nLTE_start-1:grain_RE_nLTE_end,Nc,n_lambda),stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error kabs_nLTE_CDF')
        kabs_nLTE_CDF = 0.0
     else
        low_mem_th_emission_nLTE = .true.
        write(*,*) "Using low memory mode for nLTE thermal emission"
     endif

     allocate(kdB_dT_1grain_nLTE_CDF(n_lambda,grain_RE_nLTE_start:grain_RE_nLTE_end,n_T),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kdB_dT_1grain_nLTE_CDF')
     kdB_dT_1grain_nLTE_CDF=0.0

     allocate(log_E_em_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,n_T),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error log_E_em_1grain')
     log_E_em_1grain=0.0

     allocate(xT_ech_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,Nc,nb_proc),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error xT_ech_1grain')
     xT_ech_1grain = 2
  endif

  if (lnRE) then
     allocate(E_em_1grain_nRE(grain_nRE_start:grain_nRE_end,n_T),&
          log_E_em_1grain_nRE(grain_nRE_start:grain_nRE_end,n_T), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error E_em_1grain')
     E_em_1grain_nRE=0.0
     log_E_em_1grain_nRE=0.0

     allocate(Tdust_1grain_nRE_old(grain_nRE_start:grain_nRE_end,Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Tdust_1grain_nRE_old')
     Tdust_1grain_nRE_old =0.0

     allocate(Tpeak_old(grain_nRE_start:grain_nRE_end,Nc), &
          maxP_old(grain_nRE_start:grain_nRE_end,Nc), &
          stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Tpeak')
     Tpeak_old=0
     maxP_old=0.

     allocate(xT_ech_1grain_nRE(grain_nRE_start:grain_nRE_end,Nc,nb_proc),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error xT_ech_1grain_nRE')
     xT_ech_1grain_nRE = 2

     allocate(kdB_dT_1grain_nRE_CDF(n_lambda,grain_nRE_start:grain_nRE_end,n_T),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kdB_dT_1grain_nRE_CDF')
     kdB_dT_1grain_nRE_CDF=0.0

     if (lRE_nlTE) then
        allocate(Tdust_1grain_old(grain_RE_nLTE_start:grain_RE_nLTE_end,Nc),stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error Tdust_1grain_old')
        Tdust_old=0.
     endif

     allocate(Emissivite_nRE_old(Nc,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Emissivite_nRE_old')
     Emissivite_nRE_old = 0.0
  endif

  return

end subroutine allocate_thermal_emission

!***************************************************

subroutine deallocate_thermal_emission()

  use radiation_field, only : deallocate_radiation_field

  deallocate(tab_Temp)
  if (lsed_complete) then
     deallocate(log_Qcool_minus_extra_heating)
     deallocate(DensE, DensE_m1, Dcoeff)
     call deallocate_radiation_field()
     if (allocated(nbre_reemission)) deallocate(nbre_reemission)
     if (allocated(kdB_dT_CDF)) deallocate(xT_ech,kdB_dT_CDF)
  endif

  return

end subroutine deallocate_thermal_emission

!***************************************************

subroutine deallocate_temperature_calculation()

  use radiation_field, only : deallocate_radiation_field

  if (lRE_LTE) then
     if (allocated(kdB_dT_CDF)) deallocate(xT_ech,kdB_dT_CDF)
     deallocate(log_Qcool_minus_extra_heating)
  endif
  if (lRE_nLTE) deallocate(kabs_nLTE_CDF, kdB_dT_1grain_nLTE_CDF, log_E_em_1grain,xT_ech_1grain)
  if (lreemission_stats) deallocate(nbre_reemission)
  if (lnRE) then
     deallocate(kdB_dT_1grain_nRE_CDF,E_em_1grain_nRE,log_E_em_1grain_nRE,xT_ech_1grain_nRE)
     deallocate(Tdust_1grain_nRE_old,Emissivite_nRE_old,Tpeak_old)
     if (lRE_nlTE) deallocate(Tdust_1grain_old)
  endif
  call deallocate_radiation_field()

  return

end subroutine deallocate_temperature_calculation

!***************************************************

subroutine realloc_emitting_fractions()

  integer :: alloc_status

  alloc_status = 0
  deallocate(E_disk, frac_E_stars, frac_E_disk, E_totale)
  allocate(E_disk(n_lambda2), frac_E_stars(n_lambda2), frac_E_disk(n_lambda2), E_totale(n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error frac_E_stars')
  E_disk = 0.0 ; frac_E_stars = 0.0 ; frac_E_disk = 0.0 ; E_totale = 0.0

  deallocate(prob_E_cell)
  allocate(prob_E_cell(0:n_cells,n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error prob_E_cell 2')
  prob_E_cell = 0.0

  return

end subroutine realloc_emitting_fractions

!***************************************************


subroutine allocate_temperature(Nc)

  integer, intent(in) :: Nc

  integer :: alloc_status

  allocate(Tdust(Nc), Tdust_old(Nc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error Tdust')
  Tdust = 0.0 ; Tdust_old=0.0

  if (lRE_nLTE) then
     allocate(Tdust_1grain(grain_RE_nLTE_start:grain_RE_nLTE_end,Nc),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Tdust_1grain')
     Tdust_1grain = 0.0
  endif

  if (lnRE) then
     if ( (.not.ltemp).and.(lsed.or.lmono0.or.lProDiMo.or.lProDiMo2mcfost) ) then ! si ltemp --> tableau alloue ci-dessous
        allocate(tab_Temp(n_T), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_Temp')
        tab_Temp = 0.0
     endif

     allocate(Proba_Tdust(n_T,grain_nRE_start:grain_nRE_end,Nc), &
          Tdust_1grain_nRE(grain_nRE_start:grain_nRE_end,Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error Proba_Tdust')
     Proba_Tdust=0.0
     Tdust_1grain_nRE=0.0

     allocate(l_RE(grain_nRE_start:grain_nRE_end,Nc), lchange_nRE(grain_nRE_start:grain_nRE_end,Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error l_RE')
     l_RE=.false. ; lchange_nRE = .false.
  endif

  return

end subroutine allocate_temperature

!***************************************************

subroutine repartition_wl_em()
! Pour step 1
! Corrige fonction de repartition de l'emission en prenant en compte
! le chauffage initial du disque (pour select_wl_em)
! Corrige energie des paquets du step 1
! Doit venir apres repartition energie_etoiles et repartition_energie
! C. Pinte
! 18/05/06

  implicit none

  integer :: lambda
  real(kind=dp) :: E_star_tot, E_disk_tot, E_ISM_tot, delta_wl, L_tot

  if (lTemp) then
     spectre_emission_cumul(0) = 0.0
     ! Fonction de répartition émssion
     do lambda=1,n_lambda
        delta_wl=tab_delta_lambda(lambda)*1.e-6
        spectre_emission_cumul(lambda)=spectre_emission_cumul(lambda-1) + &
             (E_stars(lambda) + E_disk(lambda) + E_ISM(lambda)) * delta_wl
     enddo

     ! Normalisation
     do lambda=1,n_lambda
        spectre_emission_cumul(lambda)=spectre_emission_cumul(lambda)/spectre_emission_cumul(n_lambda)
     enddo
  endif

  ! Energie des paquets pour step 1
  E_star_tot = 0.0
  E_disk_tot = 0.0
  E_ISM_tot = 0.0
  do lambda=1, n_lambda
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     E_star_tot = E_star_tot + E_stars(lambda) * delta_wl
     E_disk_tot = E_disk_tot + E_disk(lambda)  * delta_wl
     E_ISM_tot  = E_ISM_tot  + E_ISM(lambda)   * delta_wl
  enddo

  L_tot = 2.0*pi*hp*c_light**2 * (E_star_tot + E_disk_tot + E_ISM_tot)
  L_packet_th = L_tot/nbre_photons_tot

  return

end subroutine repartition_wl_em

!***************************************************

subroutine select_wl_em(aleat,lambda)
! Choix de la longueur d'onde dans le corps noir precedemment cree
! Dichotomie
! lambda est l'indice de la longueur d'onde
! Utilise les résultats de  repartition_wl_em
! C. Pinte

  implicit none

  real, intent(in) :: aleat
  integer, intent(out) :: lambda ! indice de longueur d'onde

  integer :: k, kmin, kmax

  kmin=0
  kmax=n_lambda

  k=(kmin + kmax)/2

  ! Dichotomie
  do while (spectre_emission_cumul(k) /= aleat)
     if (spectre_emission_cumul(k) < aleat) then
        kmin = k
     else
        kmax = k
     endif

     k = (kmin + kmax)/2

     if ((kmax-kmin) <= 1) then
        exit ! Sortie du while
     endif
  enddo   ! while

  lambda=kmax

end subroutine select_wl_em

!***************************************************

subroutine init_reemission(lheating,dudt)
! Pretabule les opacites de Planck en fonction de T et de la position
! + les termes k dB/dT pour tirage de la longueur d'onde du photon emis
! avec correction de temperature
! C. Pinte 17/01/05

  use radiation_field, only : J0

  implicit none

  logical, intent(in) :: lheating
  real, dimension(:), allocatable, intent(in), optional :: dudt

  integer :: k,lambda,t, pop, id, icell
  real(kind=dp) :: integ, coeff_exp, cst_wl, cst, wl
  real(kind=dp) ::  temp, cst_E, delta_wl
  real(kind=dp), dimension(0:n_lambda) :: integ3
  real(kind=dp), dimension(n_lambda,n_T) :: B, dB_dT

  real(kind=dp) :: Qcool, extra_heating, u_o_dt, Qcool_minus_extra_heating, Qcool0

  write(*,'(a36, $)') " Initializing thermal properties ..."

  cst_E = 2.0*hp*c_light**2 * quatre_pi

  call init_tab_temp()

  do t=1,n_T
     ! Calcul du corps noir et de sa derivee
     ! A une cst pres (pas la meme ds les 2 cas!!)

     Temp=tab_Temp(t)
     cst=cst_th/Temp
     do lambda=1, n_lambda
        ! longueur d'onde en metre
        wl = tab_lambda(lambda)*1.e-6
        delta_wl=tab_delta_lambda(lambda)*1.e-6
        cst_wl=cst/wl
        if (cst_wl < 500.0) then
           coeff_exp=exp(cst_wl)
           ! Les formules prennent en compte le delta_wl de l'integration
           B(lambda,T) = 1.0/((wl**5)*(coeff_exp-1.0))*delta_wl
           dB_dT(lambda,T) = B(lambda,T) * cst_wl*coeff_exp/(coeff_exp-1.0) !/Temp * temp a cause de dT mais ca change rien en pratique
        else
           B(lambda,T)=0.0
           dB_dT(lambda,T)=0.0
        endif
     enddo !lambda
  enddo ! T

  ! produit par opacite (abs seule) massique
  ! Todo : this loop is OK in 2D but takes ~ 5sec for 0.5 million cells in 3D
  !$omp parallel default(none) &
  !$omp private(id,icell,T,lambda,integ, Qcool,Qcool0,extra_heating,Qcool_minus_extra_heating,Temp,u_o_dt) &
  !$omp shared(cst_E,kappa_abs_LTE,kappa_factor,volume,B,lextra_heating,xT_ech,log_Qcool_minus_extra_heating,J0) &
  !$omp shared(n_T,n_cells,p_n_cells,n_lambda,tab_Temp,ldudt_implicit,ufac_implicit,dudt,lRE_nLTE,lvariable_dust,icell_ref)
  id = 1 ! Pour code sequentiel
  !$ id = omp_get_thread_num() + 1

  !$omp do
  do icell=1,p_n_cells
     do T=1, n_T
        Temp = tab_Temp(T)
        integ=0.0
        do lambda=1, n_lambda
           ! kappa en Au-1    \
           ! volume en AU3     >  pas de cst pour avoir E_em en SI
           ! B * cst_E en SI = W.m-2.sr-1 (.m-1 * m) cat delta_wl inclus
           integ = integ + kappa_abs_LTE(icell,lambda) * B(lambda,T)  ! kappa_factor, and volume are not included here
        enddo !lambda

        ! Proper normalization
        ! At this stage, this is Q- in W / (AU/m)^2
        Qcool = integ * cst_E
        if (T==1) Qcool0 = Qcool

        ! Non radiative heating :
        ! We solve Q+ = int kappa.Jnu.dnu = Q- - extra_heating = int kappa.Bnu.dnu - extra_heating
        ! Here we conpute the extra_heating term
        if (.not.lextra_heating) then
           ! Energie venant de l'equilibre avec nuage à T_min
           extra_heating = Qcool0
        else
           if (ldudt_implicit) then
              u_o_dt = ufac_implicit * Temp ! u(T)/dt
              ! dudt is meant to be u^n/dt here
              extra_heating = max(Qcool0, (u_o_dt - dudt(icell)) / (AU_to_m**2 * volume(icell) * kappa_factor(icell)))
           else
              extra_heating = max(Qcool0, dudt(icell) / (AU_to_m**2 * volume(icell) * kappa_factor(icell)) )
           endif
        endif

        Qcool_minus_extra_heating = Qcool - extra_heating
        if (Qcool_minus_extra_heating > tiny_dp) then
           ! ToDO : I will need to correct the log interpolation as it should now broken to to the shift in energy
           log_Qcool_minus_extra_heating(T,icell)=log(Qcool_minus_extra_heating)
        else
           ! This sets the initial disk temperature
           xT_ech(icell,id)=T+1
           log_Qcool_minus_extra_heating(T,icell)=-1000.
        endif

        if (lRE_nLTE.and.(T==1)) then
           do lambda=1, n_lambda
              J0(icell,lambda) =  volume(icell) * B(lambda,T) * cst_E
           enddo
        endif
     enddo ! T

  enddo !icell
  !$omp enddo
  !$omp end parallel

  if (low_mem_th_emission) then
     do T=1, n_T
        do k=grain_RE_LTE_start,grain_RE_LTE_end
           integ3(1) = 0.0
           do lambda=2, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1), ni densite
              integ3(lambda) = integ3(lambda-1) + C_abs_norm(k,lambda) * dB_dT(lambda,T)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_1grain_LTE_CDF(lambda,k,T) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     enddo ! T
  else  ! .not. low_mem_th_emission
     do icell=1,p_n_cells
        do T=1, n_T
           integ3(0) = 0.0
           do lambda=1, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1)
              integ3(lambda) = integ3(lambda-1) + kappa_abs_LTE(icell,lambda) * kappa_factor(icell) * dB_dT(lambda,T)
           enddo !l

           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_CDF(lambda,T,icell) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo ! T
     enddo !icell
  endif ! low_mem_th_emission

  if (lRE_nLTE) then
     ! produit par opacite (abs seule) massique d'une seule taille de grain
     do T=1, n_T
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ=0.0
           do lambda=1, n_lambda
              ! WARNING : il manque la densite et le volume par rapport au cas LTE !!!
              integ = integ + C_abs_norm(k,lambda) * B(lambda,T)
           enddo !lambda
           ! Le coeff qui va bien
           if (integ > tiny_dp) then
              log_E_em_1grain(k,T)=log(integ*cst_E)
           else
              log_E_em_1grain(k,T)=-1000.
           endif
        enddo ! k

        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ3(1) = 0.0
           do lambda=2, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1), ni densite
              integ3(lambda) = integ3(lambda-1) + C_abs_norm(k,lambda) * dB_dT(lambda,T)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_1grain_nLTE_CDF(lambda,k,T) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     enddo ! T
  endif ! lnLTE

  if (lnRE) then
     do T=1, n_T
        do k=grain_nRE_start,grain_nRE_end
           integ=0.0
           do lambda=1, n_lambda
              integ = integ + C_abs_norm(k,lambda) * B(lambda,T)
           enddo !lambda
           integ = integ * cst_E
           E_em_1grain_nRE(k,T) = integ
           if (integ > 0.0_dp) then
              log_E_em_1grain_nRE(k,T) = log(integ)
           else
              write(*,*) "Error in init_reemission"
              write(*,*) "Pb in opacity of non-equilibrium grains"
              write(*,*) "grain", k, "C_abs_norm seems to be incorrect"
              write(*,*) C_abs_norm(k,:)
              write(*,*) "Exiting"
              call exit(1)
           endif
        enddo

        do k=grain_nRE_start,grain_nRE_end
           integ3(1) = 0.0
           do lambda=2, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1), ni densite
              integ3(lambda) = integ3(lambda-1) + C_abs_norm(k,lambda) * dB_dT(lambda,T)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_1grain_nRE_CDF(lambda,k,T) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     enddo ! T
  endif !lnRE

  if (lextra_heating .and. ldudt_implicit) then
     ! as u(T) depends on T, we need to check that int kappa_nu.Bnu(T).dnu - (u(T) - u_n)/dt
     ! is still an increasing function of T
     do icell=1, p_n_cells
        do T=1,n_T-1
           if (log_Qcool_minus_extra_heating(T+1,icell) < log_Qcool_minus_extra_heating(T,icell)) then
              call error("Qrad_minus_dudt is not an increasing function of T")
           endif
        enddo
     enddo
  endif

  do pop=1, n_pop
     if (dust_pop(pop)%methode_chauffage == 3) then
        if (dust_pop(pop)%is_Misselt_opacity_file) call read_Misselt_specific_heat(pop)
        if (dust_pop(pop)%is_DustEM_opacity_file)  call read_DustEM_specific_heat(pop)
     endif
  enddo

  write(*,*) "Done"
  return

end subroutine init_reemission


!********************************************************************

subroutine Temp_LTE(icell, Ti, Temp)

  use radiation_field, only : xKJ_abs

  integer, intent(in) :: icell
  integer, intent(out) :: Ti
  real, intent(out) :: Temp

  real(kind=dp) :: Qheat, log_Qheat, frac
  integer :: p_icell

  if (lvariable_dust) then
     p_icell = icell
  else
     p_icell = icell_ref
  endif

  Qheat=sum(xKJ_abs(icell,:)) * L_packet_th / volume(icell)  ! does not include kappa_factor, same as log_Qcool
  if (Qheat < tiny_dp) then
     Temp = T_min ; Ti = 2
  else
     log_Qheat = log(Qheat)

     if (log_Qheat <  log_Qcool_minus_extra_heating(1,p_icell)) then
        Temp = T_min ; Ti = 2
     else
        ! Temperature echantillonee juste sup. a la temperature de la cellule
        Ti = maxval(xT_ech(icell,:))

        ! On incremente eventuellement la zone de temperature
        do while((log_Qcool_minus_extra_heating(Ti,p_icell) < log_Qheat).and.(Ti < n_T))
           Ti=Ti+1
        enddo  ! LIMITE MAX

        ! Interpolation lineaire entre energies emises pour des
        ! temperatures echantillonees voisines
        frac=(log_Qheat-log_Qcool_minus_extra_heating(Ti-1,p_icell)) / &
             (log_Qcool_minus_extra_heating(Ti,p_icell)-log_Qcool_minus_extra_heating(Ti-1,p_icell))
        Temp=exp(log(tab_Temp(Ti))*frac+log(tab_Temp(Ti-1))*(1.0-frac))
     endif
  endif

  return

end subroutine Temp_LTE

!********************************************************************

subroutine im_reemission_LTE(id,icell,p_icell,aleat1,aleat2,lambda)
! Calcul de la temperature de la cellule et stokage energie recue + T
! Reemission d'un photon a la bonne longeur d'onde

  use radiation_field, only : xKJ_abs, E0

  implicit none

  integer, intent(in) :: id, icell, p_icell
  real, intent(in) ::  aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, Ti, T1, T2, k, heating_method
  real :: Temp
  real(kind=dp) :: frac_T1, frac_T2, proba

  if (lreemission_stats) nbre_reemission(icell,id) = nbre_reemission(icell,id) + 1.0_dp

  call Temp_LTE(icell, Ti, Temp)

  ! Save pour prochaine reemission et/ou T finale
  xT_ech(icell,id) = Ti

  !**********************************************************************
  ! Choix de la longeur d'onde de reemission
  ! Dichotomie, la loi de proba est obtenue par interpolation lineaire en T
  T2 = Ti ; T1 = Ti-1
  Frac_T2=(Temp-tab_Temp(T1))/(tab_temp(T2)-tab_Temp(T1))
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  if (low_mem_th_emission) then
     ! Select absorbing grain
     heating_method = 1 ! LTE
     k = select_absorbing_grain(lambda,icell, aleat1, heating_method)

     ! Select wavelength
     do while((l2-l1) > 1)
        proba=frac_T1*kdB_dT_1grain_LTE_CDF(l,k,T1)+frac_T2*kdB_dT_1grain_LTE_CDF(l,k,T2)
        if(aleat2 > proba) then
           l1=l
        else
           l2=l
        endif
        l=(l1+l2)/2
     enddo
  else
     ! Select wavelength
     do while((l2-l1) > 1)
        proba=frac_T1*kdB_dT_CDF(l,T1,p_icell)+frac_T2*kdB_dT_CDF(l,T2,p_icell)
        if(aleat2 > proba) then
           l1=l
        else
           l2=l
        endif
        l=(l1+l2)/2
     enddo
  endif

  lambda=l+1

  return

end subroutine im_reemission_LTE

!********************************************************************

subroutine im_reemission_NLTE(id,icell,p_icell,aleat1,aleat2,lambda)
! Calcul de la temperature de la cellule
! Reemission d'un photon a la bonne longeur d'onde

  use radiation_field, only : xJ_abs, J0

  implicit none

  integer, intent(in) :: id,icell, p_icell
  real, intent(in) :: aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, kmin, kmax, lambda0, ilambda, heating_method
  real(kind=dp) :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_E_abs, J_abs

  lambda0=lambda

  ! Selection du grain qui absorbe le photon
  if (low_mem_th_emission_nLTE) then
     ! Select absorbing grain
     heating_method = 2 ! nLTE
     k = select_absorbing_grain(lambda0,icell, aleat1, heating_method)
  else
     kmin=grain_RE_nLTE_start
     kmax=grain_RE_nLTE_end
     k=(kmin+kmax)/2

     do while((kmax-kmin) > 1)
        if (kabs_nLTE_CDF(k,icell,lambda0) < aleat1) then
           kmin = k
        else
           kmax = k
        endif
        k = (kmin + kmax)/2
     enddo   ! while
     k=kmax
  endif

  ! Mean intensity
  ! Somme sur differents processeurs
  J_abs=0.0
  do ilambda=1, n_lambda
     J_abs =  J_abs + C_abs_norm(k,ilambda)  * (sum(xJ_abs(icell,ilambda,:)) + J0(icell,ilambda))
  enddo ! lambda
 ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans E_em_1grain
  log_E_abs=log(J_abs*L_packet_th/volume(icell) )

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech_1grain(k,icell,:))

  ! On incremente eventuellement la zone de temperature
  do while((log_E_em_1grain(k,T_int) < log_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! LIMITE MAX

  ! Save pour prochaine reemission et/ou T finale
  xT_ech_1grain(k,icell,id)=T_int

  ! Interpolation lineaire entre energies emises pour des
  ! temperatures echantillonees voisines
  T2=T_int
  Temp2=tab_Temp(T2)
  T1=T_int-1
  Temp1=tab_Temp(T1)
  frac=(log_E_abs-log_E_em_1grain(k,T1))/(log_E_em_1grain(k,T2)-log_E_em_1grain(k,T1))
  Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))


  !**********************************************************************
  ! Choix de la longeur d'onde de reemission
  ! Dichotomie, la loi de proba est obtenue par interpolation lineaire en T
  frac_T2=(Temp-Temp1)/(Temp2-Temp1)
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  do while((l2-l1) > 1)
     proba=frac_T1*kdB_dT_1grain_nLTE_CDF(l,k,T1)+frac_T2*kdB_dT_1grain_nLTE_CDF(l,k,T2)
     if(aleat2 > proba) then
        l1=l
     else
        l2=l
     endif
     l=(l1+l2)/2
  enddo
  lambda=l+1

  return

end subroutine im_reemission_NLTE

!********************************************************************

subroutine Temp_finale()
! Calcule la temperature finale du disque dans chaque cellule
! Cas de l'"absoprtion continue"
! C. Pinte
! 24/01/05

  use radiation_field, only : xKJ_abs, E0

  implicit none

  real :: Temp
  integer :: Ti, icell

  !$omp parallel &
  !$omp default(none) &
  !$omp private(icell,Temp,Ti) &
  !$omp shared(Tdust,n_cells)
  !$omp do schedule(dynamic,10)
  do icell=1,n_cells
     call Temp_LTE(icell, Ti, Temp)
     Tdust(icell) = Temp
  enddo !icell
  !$omp enddo
  !$omp end parallel

  if (maxval(Tdust) > T_max) then
     write(*,*) "WARNING : temperature > sublimation temperature"
     write(*,*) "WARNING : temperature = ", maxval(Tdust)
  else
     write(*,*) "Max. temperature = ", maxval(Tdust)
  endif

  return

end subroutine Temp_finale

!**********************************************************************

subroutine set_min_temperature(Tmin)

  real, intent(in) :: Tmin

  integer :: icell, i

  i=0
  do icell=1, n_cells
     if (Tdust(icell) < Tmin) then
        Tdust(icell) = Tmin
        i = i+1
     endif
  end do

  write(*,*) "Tdust was set to ", Tmin, "in", i, "cells"

  return

end subroutine set_min_temperature

!**********************************************************************

subroutine Temp_finale_nLTE()
! Calcule la temperature finale du disque dans chaque cellule
! Cas de l'"absoprtion continue"
! C. Pinte
! 24/01/05

  use radiation_field, only : xJ_abs, J0

  implicit none

  integer :: T_int, T1, T2,icell
  real(kind=dp) :: Temp, Temp1, Temp2, frac, log_E_abs, J_absorbe

  integer :: k, lambda

  ! Calcul de la temperature de la cellule et stokage energie recue + T
  ! Utilisation temperature precedente

  !$omp parallel &
  !$omp default(none) &
  !$omp private(log_E_abs,T_int,T1,T2,Temp1,Temp2,Temp,frac,icell,J_absorbe) &
  !$omp shared(L_packet_th,Tdust,tab_Temp,n_cells,n_lambda) &
  !$omp shared(xJ_abs,densite_pouss,Tdust_1grain, xT_ech_1grain,log_E_em_1grain) &
  !$omp shared(C_abs_norm,volume, grain_RE_nLTE_start, grain_RE_nLTE_end, n_T, T_min, J0)
  !$omp do schedule(dynamic,10)
  do icell=1,n_cells
     do k=grain_RE_nLTE_start, grain_RE_nLTE_end
        if (densite_pouss(k,icell) > tiny_dp) then
           J_absorbe=0.0
           do lambda=1, n_lambda
              J_absorbe =  J_absorbe + C_abs_norm(k,lambda)  * (sum(xJ_abs(icell,lambda,:)) + J0(icell,lambda))
           enddo ! lambda

           ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans E_em_1grain
           J_absorbe = J_absorbe*L_packet_th/volume(icell)
           if (J_absorbe < tiny_dp) then
              Tdust_1grain(k,icell) = T_min
           else
              log_E_abs=log(J_absorbe)

              if (log_E_abs <  log_E_em_1grain(k,1)) then
                 Tdust_1grain(k,icell) = T_min
              else
                 ! Temperature echantillonee juste sup. a la temperature de la cellule
                 T_int=maxval(xT_ech_1grain(k,icell,:))

                 ! On incremente eventuellement la zone de temperature
                 do while((log_E_em_1grain(k,T_int) < log_E_abs).and.(T_int < n_T))
                    T_int=T_int+1
                 enddo  ! LIMITE MAX

                 ! Interpolation lineaire entre energies emises pour des
                 ! temperatures echantillonees voisines
                 T2=T_int
                 Temp2=tab_Temp(T2)
                 T1=T_int-1
                 Temp1=tab_Temp(T1)
                 frac=(log_E_abs-log_E_em_1grain(k,T1))/(log_E_em_1grain(k,T2)-log_E_em_1grain(k,T1))
                 Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

                 ! Save
                 Tdust_1grain(k,icell)=Temp
              endif
           endif
        else
           Tdust_1grain(k,icell)=0.0
        endif
     enddo !k
  enddo !icell
  !$omp enddo
  !$omp end parallel

  if (maxval(Tdust_1grain) > T_max) then
     write(*,*) "WARNING : temperature > sublimation temperature"
     write(*,*) "WARNING : temperature = ", maxval(Tdust_1grain)
  else
     write(*,*) "Max. temperature NLTE = ", maxval(Tdust_1grain)
  endif

  return

end subroutine Temp_finale_nLTE

!**********************************************************************

subroutine Temp_nRE(lconverged)
! C. Pinte
! 26/01/07

  use radiation_field, only : xJ_abs, J0

  implicit none

!  logical, intent(in) :: lfinal ! TODO
   logical, intent(out) :: lconverged

  real, parameter :: precision = 5.e-2

  real(kind=dp), dimension(:,:,:), allocatable :: A, B
  real(kind=dp), dimension(:,:), allocatable :: X
  real(kind=dp), dimension(n_T) :: nu_bin, delta_nu_bin
  real(kind=dp), dimension(0:n_T) :: T_lim, U_lim

  real(kind=dp), dimension(n_lambda) :: C_abs_norm_o_dnu, tab_nu

  real(kind=dp), dimension(:,:), allocatable :: kJnu, lambda_Jlambda

  real(kind=dp) :: delta_T, kJnu_interp, t_cool, t_abs, mean_abs_E, kTu, mean_abs_nu, cst_t_cool, Int_k_lambda_Jlambda

  integer :: l, T, T1, T2, lambda, alloc_status, id, T_int, icell

  integer, parameter :: n_cooling_time = 100

  real(kind=dp), dimension(n_cooling_time+1) :: en_lim
  real(kind=dp), dimension(n_cooling_time) :: en, delta_en, Cabs

  real :: E_max
  real(kind=dp) :: Temp1, Temp2, frac, log_E_abs, maxP
  real(kind=dp) :: somme1, somme2, wl, wl_mum, delta_wl, Temp, cst_wl ! pour test
  integer :: Tpeak

  write(*,*) "Calculating temperature probabilities of non equilibrium grains ..."
  if (lforce_PAH_equilibrium) write(*,*) "Forcing equilibrium " ;

  allocate(A(n_T,n_T,nb_proc), B(n_T,n_T,nb_proc), X(n_T,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error A')
  A=0.0 ; B=0.0 ;

  allocate(kJnu(n_lambda,nb_proc), lambda_Jlambda(n_lambda,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error kJnu')
  kJnu = 0.0 ;  lambda_Jlambda=0.0

  delta_T=exp((1.0_dp/(real(n_T,kind=dp)))*log(T_max/T_min))
  T_lim(0) = T_min
  T_lim(1:n_T) = tab_Temp(:)*sqrt(delta_T)

  tab_nu = c_light/(tab_lambda*1.0e-6)

  ! l est l'indice de taille de grains
  do l=grain_nRE_start, grain_nRE_end
     ! bin centers [Hz]
     nu_bin(:) = specific_heat(real(tab_Temp,kind=dp),l) * tab_Temp(:)/hp ! evaluate specific heat & enthalpy at bin centers

     ! bin widths [Hz]
     U_lim(:) = specific_heat(T_lim(:),l)*T_lim(:) ! enthalpy of bin boundaries [erg]
     delta_nu_bin(:) = (U_lim(1:n_T)-U_lim(0:n_T-1))/hp  ! bin width [Hz]

     ! compute transition matrix  A   [Hz]

     !! cooling terms :  select one-off upper diagonal, on peut les calculer a l'avance
     ! Afi with final < inital
     A(:,:,:)=0.
     do T=2,n_T
        A(T-1,T,:) = E_em_1grain_nRE(l,T) / (nu_bin(T)-nu_bin(T-1)) ! OK, il manque un 1/h par rapport a Krugel p265
        !if (tab_Temp(T) < 10.)  A(T-1,T,:) =  A(T-1,T,:) * 1e5
     enddo

     ! Pour kJnu
     do lambda=1, n_lambda
        ! C_abs_norm / dnu  : dnu = c dwl/wl^2    1/dnu = wl^2/(c.dwl)
        C_abs_norm_o_dnu(lambda) = C_abs_norm(l,lambda) *  tab_lambda(lambda)**2 / (c_light*tab_delta_lambda(lambda)) * 1.0e-6
        !C_abs_norm_o_dnu(lambda) = C_abs_norm(l,lambda) *  tab_lambda(lambda) / (c_light*tab_delta_lambda(lambda)) * 1.0e-6
        !C_abs_norm_o_dnu(lambda) = C_abs_norm(l,lambda)  * (tab_lambda(lambda))/c_light * 1.0e-6
     enddo

     ! Pour cooling time
     E_max = real(maxval(U_lim))
     en_lim = spanl(1.0e-10*E_max,E_max, n_cooling_time+1)
     do T=1,n_cooling_time
        en(T) = 0.5 * (en_lim(T+1) + en_lim(T))
        delta_en(T) = en_lim(T+1) - en_lim(T)
        Cabs(T) = interp(real(C_abs_norm(l,:),kind=dp),real(tab_lambda,kind=dp),en(T)/hp)
     enddo

     cst_t_cool =(real(hp,kind=dp)**3 * real(c_light,kind=dp)**2) / (8*pi)

     ! Boucle sur les cellules

     !$omp parallel &
     !$omp default(none) &
     !$omp private(Int_k_lambda_Jlambda, lambda, wl, wl_mum, T, T2) &
     !$omp private(kJnu_interp,id,t_cool,t_abs,mean_abs_E,mean_abs_nu,kTu) &
     !$omp private(frac,T1,Temp1,Temp2,T_int,log_E_abs,icell) &
     !$omp shared(l,kJnu, lambda_Jlambda, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium) &
     !$omp shared(n_cells, C_abs_norm_o_dnu, xJ_abs, J0, L_packet_th, volume, n_T, disk_zone,etoile) &
     !$omp shared(tab_nu, n_lambda, tab_delta_lambda, tab_lambda,en,delta_en,Cabs) &
     !$omp shared(delta_nu_bin,Proba_Tdust, A,B,X,nu_bin,tab_Temp,T_min,T_max,lbenchmark_SHG,lMathis_field,Mathis_field) &
     !$omp shared(Tdust_1grain_nRE,log_E_em_1grain_nRE,cst_t_cool,C_abs_norm,l_RE,r_grid) &
     !$omp shared(densite_pouss,l_dark_zone,Tdust,lchange_nRE)

     id = 1 ! pour code sequentiel
     ! ganulation faible car le temps calcul depend fortement des cellules
     !$omp do schedule(dynamic,1)
     do icell=n_cells, 1, -1
        !$ id = omp_get_thread_num() + 1
        if (l_dark_zone(icell)) then
           l_RE(:,icell) = .true.
        else
           if (densite_pouss(l,icell) > tiny_dp) then
              ! Champ de radiation
              Int_k_lambda_Jlambda=0.0
              do lambda=1, n_lambda
                 ! conversion lambda et delta_lambda de micron -> m
                 lambda_Jlambda(lambda,id) = (sum(xJ_abs(icell,lambda,:)) + J0(icell,lambda))
                 Int_k_lambda_Jlambda = Int_k_lambda_Jlambda + C_abs_norm(l,lambda) * lambda_Jlambda(lambda,id)
                 kJnu(lambda,id) =   C_abs_norm_o_dnu(lambda)  * lambda_Jlambda(lambda,id)
              enddo ! lambda

              lambda_Jlambda(:,id) =  lambda_Jlambda(:,id)*L_packet_th * (1.0/volume(icell))
              Int_k_lambda_Jlambda = Int_k_lambda_Jlambda*L_packet_th * (1.0/volume(icell))
              kJnu(:,id) = kJnu(:,id)*L_packet_th * (1.0/volume(icell))

              if (lbenchmark_SHG) then ! Adding TRUST radiation field
                 Int_k_lambda_Jlambda = 0.0
                 if (lMathis_field) then
                    id = 1
                    do lambda=1, n_lambda
                       wl_mum =  tab_lambda(lambda)
                       wl = tab_lambda(lambda) * 1e-6

                       if (wl_mum < 0.0912) then
                          lambda_Jlambda(lambda,id) = 0.0_dp
                       else if (wl_mum < 0.110) then
                          lambda_Jlambda(lambda,id) = 3069 * wl_mum**3.4172
                       else if (wl_mum <  0.134) then
                          lambda_Jlambda(lambda,id) = 1.627
                       else if (wl_mum <  0.250) then
                          lambda_Jlambda(lambda,id) = 0.0566 * wl_mum**(-1.6678)
                       else
                          lambda_Jlambda(lambda,id) = 1e-14 * Blambda_dp(wl,7500.) &
                               + 1e-13 * Blambda_dp(wl,4000.) &
                               + 4e-13 * Blambda_dp(wl,3000.)
                       endif
                       lambda_Jlambda(lambda,id) = lambda_Jlambda(lambda,id)* wl * Mathis_field * 1.3e-2

                       Int_k_lambda_Jlambda = Int_k_lambda_Jlambda + C_abs_norm(l,lambda) * lambda_Jlambda(lambda,id)
                       kJnu(lambda,id) =   C_abs_norm_o_dnu(lambda)  * lambda_Jlambda(lambda,id)
                    enddo ! lambda
                 else ! BB field
                    do lambda=1, n_lambda
                       wl = tab_lambda(lambda)*1e-6
                       lambda_Jlambda(lambda,id) = (Blambda_dp(wl,etoile(1)%T))*wl / disk_zone(1)%Rin**2 * 7.e-8 !* 7.25965e-08

                       Int_k_lambda_Jlambda = Int_k_lambda_Jlambda + C_abs_norm(l,lambda) * lambda_Jlambda(lambda,id)
                       kJnu(lambda,id) =   C_abs_norm_o_dnu(lambda)  * lambda_Jlambda(lambda,id)
                    enddo ! lambda
                 endif ! lMathis
              endif ! End adding TRUST radiation field

              ! decide whether we really need to use this model, instead of calculating the equilibrium temperature

              ! time interval for photon absorption
              t_abs = sum(C_abs_norm(l,:)*lambda_Jlambda(:,id)*tab_lambda(:)) * 1.0e-6/c_light
              if (t_abs > tiny_dp) then
                 t_abs = hp/(4.*pi*t_abs)
              else
                 t_abs = huge_real
              endif

              ! mean absorbed photon energy
              mean_abs_E = t_abs * 4*pi * Int_k_lambda_Jlambda ! D01 eq 46
              mean_abs_nu = mean_abs_E/hp

              ! Inversion rapide Energie en temperature
              ! On prend energie juste en dessus energie moyenne -> t_cool un peu trop court
              kTu = kb * tab_Temp(1)
              search : do T=2, n_T
                 if (nu_bin(T) > mean_abs_nu) then
                    kTu = kb * tab_Temp(T)
                    exit search
                 endif
              enddo search

              ! radiative cooling time at mean photon energy
              t_cool = 0.0
              if (T <= n_T) then
                 integ : do T=1, n_cooling_time
                    if (en(T) > mean_abs_E) exit integ
                    if (en(T)/kTu < 500.) then
                       t_cool = t_cool + en(T)**3 * Cabs(T)/(exp(en(T)/kTu)-1.0_dp) * delta_en(T)
                    endif
                 enddo integ
                 if (t_cool > tiny_dp) then
                    t_cool = cst_t_cool*mean_abs_E/t_cool * 1.0e-4
                 else
                    t_cool = huge_real
                 endif
              endif

              ! Calcul Temperature equilibre
              if (Int_k_lambda_Jlambda < tiny_dp) then
                 Tdust_1grain_nRE(l,icell) = T_min
              else
                 log_E_abs=log(Int_k_lambda_Jlambda)

                 if (log_E_abs <  log_E_em_1grain_nRE(l,1)) then
                    Tdust_1grain_nRE(l,icell) = T_min
                 else
                    T_int=2
                    do while((log_E_em_1grain_nRE(l,T_int) < log_E_abs).and.(T_int < n_T))
                       T_int=T_int+1
                    enddo  ! LIMITE MAX

                    ! Interpolation lineaire entre energies emises pour des
                    ! temperatures echantillonees voisines
                    T2=T_int
                    Temp2=tab_Temp(T2)
                    T1=T_int-1
                    Temp1=tab_Temp(T1)
                    frac=(log_E_abs-log_E_em_1grain_nRE(l,T1))/&
                         (log_E_em_1grain_nRE(l,T2)-log_E_em_1grain_nRE(l,T1))
                    Tdust_1grain_nRE(l,icell)=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))
                 endif
              endif

              if (Tdust_1grain_nRE(l,icell) > T_max) then
                 ! Impossible de definir proba de temperature
                 t_cool = 1.0 ; t_abs = 0.0
                 write(*,*) "ERROR : temperature of non equilibrium grains is larger than", T_max
                 write(*,*) "cell", icell, "R=", real(r_grid(icell)), real(densite_pouss(l,icell)), &
                      real(Tdust_1grain_nRE(l,icell))
                 write(*,*) "Exiting"
                 call exit(1)
              endif

              ! Forcing equilibrium or no equilibrium
              if (lforce_PAH_equilibrium) then
                 t_cool=1.0 ; t_abs = 0.0
              else if (lforce_PAH_out_equilibrium) then
                 t_cool = 0.0 ; t_abs = 1.0
              endif

              !t_cool = 0.0

              if (t_cool < t_abs) then !  calcul proba temperature
                 l_RE(l,icell) = .false.
                 ! heating terms : select lower triangle
                 ! Afi with with f > i
                 do T=1,n_T
                    do T2=1,T-1
                       ! Tres largement plus stable en interpolation lineaire que log !!!
                       ! Ca fait n'importe quoi en log en cas de flux faible
                       !KJ_abs_interp = exp(interp(log(kJnu(:,id)), log(tab_nu), log(nu_bin(T) - nu_bin(T2))))
                       !A(T,T2,id) = 0.5 * KJ_abs_interp * delta_nu_bin(T)/ (nu_bin(T) - nu_bin(T2))
                       kJnu_interp = interp(kJnu(:,id), tab_nu, nu_bin(T) - nu_bin(T2))
                       A(T,T2,id) = kJnu_interp * delta_nu_bin(T)/ (nu_bin(T) - nu_bin(T2))
                    end do
                 end do

                 !! diagonal terms
                 !do T=1,n_T
                 !   A(T,T,id) = -sum(A(T,:,id)) ! sum along k direction (diag. should be zero initially)
                 !enddo

                 !! compute probabilities P(T) using iterative formulation of GD89

                 !! compute B array
                 B(:,:,id) = 0.
                 B(n_T,1:n_T-1,id) = A(n_T,1:n_T-1,id)
                 do T=n_T-1,1,-1
                    do T2=1,T-1
                       B(T,T2,id) = A(T,T2,id) + B(T+1,T2,id)
                    enddo
                 enddo

                 !! compute P array
                 X(1,id) = 1.e-250_dp
                 do T=2, n_T
                    if (X(T-1,id) <  1.e-300_dp) X(T-1,id) = 1.e-300_dp
                    if (X(T-1,id) >  1.e250_dp) X(1:T-1,id) = X(1:T-1,id) * 1.0e-50_dp ! permet de stabiliser en cas d'erreur d'arrondis
                    X(T,id) =  sum(B(T,1:T-1,id)*X(1:T-1,id)) / max(A(T-1, T,id),tiny_dp)
                    if (A(T-1, T,id) < tiny_dp) then
                       write(*,*) l, T, id, "normalization error"
                       call exit(1)
                    endif

                 enddo


                 !X(1,id) = 1.e-250_dp
                 !do T=2, n_T
                 !   if (X(T-1,id) <  1.e-300_dp) X(T-1,id) = 1.e-300_dp
                 !   if (X(T-1,id) >  1.e250_dp) X(1:T-1,id) = X(1:T-1,id) * 1.0e-50_dp ! permet de stabiliser en cas d'erreur d'arrondis
                 !   X(T,id) =  abs(sum(A(1:T-1,T-1,id)*X(1:T-1,id)) / max(A(T-1, T,id),tiny_dp))
                 !enddo

                 !! Normalisation
                 X(1,id) = X(2,id)
                 X(:,id)=X(:,id)/sum(X(:,id))

                 do T=1, n_T ! boucle car passage tableau double -> simple bug avec ifort
                    Proba_Tdust(T,l,icell) = X(T,id) ! P(T).dT, probability of being in bin of temperature T
                 enddo

              else  ! test : t_cool > t_abs : on est a a l'equilibre

                 if (.not.l_RE(l,icell)) then ! on n'etait pas a l'equilibre a l'iteration precedente
                    l_RE(l,icell) = .true.
                    lchange_nRE(l,icell) = .true.
                 else
                    lchange_nRE(l,icell) = .false. ! on etait deja a l'equilibre
                 endif

              endif ! test : t_cool vs t_abs
           endif ! densite_pouss > 0.
        endif ! l_dark_zone
     enddo !icell
     !$omp end do
     !$omp end parallel
  enddo !l

  deallocate(A, B, X)

  write(*,*) "Done"

  lconverged = .true.
  ! Verification convergence sur la temperature
  if (lRE_LTE) then
     do icell=1,n_cells
        delta_T = Tdust(icell) - Tdust_old(icell)
        if (delta_T > precision * Tdust_old(icell)) lconverged = .false.
     enddo
     Tdust_old = Tdust
  endif

  if (lRE_nLTE) then
     do icell=1, n_cells
        do l=grain_RE_nLTE_start,grain_RE_nLTE_end
           delta_T = Tdust_1grain(l,icell) - Tdust_1grain_old(l,icell)
           if (delta_T > precision * Tdust_1grain_old(l,icell)) lconverged = .false.
        enddo
     enddo
     Tdust_1grain_old = Tdust_1grain
  endif

  do icell=1, n_cells
     do l=grain_nRE_start, grain_nRE_end
        if (l_RE(l,icell)) then
           delta_T = Tdust_1grain_nRE(l,icell) - Tdust_1grain_nRE_old(l,icell)
           if (delta_T > precision * Tdust_1grain_nRE_old(l,icell)) lconverged = .false.
           Tdust_1grain_nRE_old(l,icell) =  Tdust_1grain_nRE(l,icell)
        else
           Tpeak = 1
           maxP =  Proba_Tdust(1,l,icell)
           do T=2, n_T
              if (Proba_Tdust(T,l,icell) > maxP) then
                 maxP=Proba_Tdust(T,l,icell)
                 Tpeak=T
              endif
           enddo

           if (Tpeak /= Tpeak_old(l,icell)) lconverged=.false.

           if ( abs(maxP-maxP_old(l,icell)) > precision * maxP_old(l,icell) ) lconverged=.false.

           Tpeak_old(l,icell) = Tpeak
           maxP_old(l,icell) = maxP
        endif
     enddo !l
  enddo !icell


  if (lbenchmark_SHG) then ! renormalisation des proba, probleme en transfert complet ...
     do l=grain_nRE_start,grain_nRE_end ! TEST pour 1 taille de grains
        do icell=1, n_cells
           somme1=0.0
           somme2=0.0
           !temp=Tdust_1grain(i,j,l) ! Tab_Temp(100)
           do lambda=1, n_lambda
              somme1 = somme1 + C_abs_norm(l,lambda)  * lambda_Jlambda(lambda,1)
              wl = tab_lambda(lambda)*1.0e-6
              delta_wl = tab_delta_lambda(lambda)*1.0e-6
              if (l_RE(l,icell)) then
                 Temp = Tdust_1grain_nRE(l,icell)
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < 500.) then
                    somme2 = somme2 + C_abs_norm(l,lambda)* 1.0/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
                 endif
              else
                 do T=1, n_T
                    Temp = tab_Temp(T)
                    cst_wl=cst_th/(Temp*wl)
                    if (cst_wl < 500.) then
                       somme2 = somme2 + C_abs_norm(l,lambda)* 1.0/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl * &
                            Proba_Tdust(T,l,icell)
                    endif !cst_wl
                 enddo
              endif
           enddo ! lambda

           somme2=somme2*2.0*hp*c_light**2
           !write(*,*) i,j,l,l_RE(i,j,l), Tdust_1grain_nRE(i,j,l), real(somme1), real(somme2), real(somme1/somme2)
           if (somme2 > tiny_dp) then
              Proba_Tdust(:,l,icell)  = Proba_Tdust(:,l,icell) * real(somme1/somme2)
           endif
        enddo !icell
     enddo ! l
  endif ! lbenchmark_SHG

!  lconverged = .true. ! TMP

  return

end subroutine Temp_nRE

!********************************************************************

subroutine im_reemission_qRE(id,icell,p_icell,aleat1,aleat2,lambda)
! Calcul de la temperature de la cellule
! Reemission d'un photon a la bonne longeur d'onde

  use radiation_field, only : xJ_abs, J0

  implicit none

  integer, intent(in) :: id, icell, p_icell ! p_icell usused
  real, intent(in) :: aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, lambda0, ilambda
  real(kind=dp) :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_E_abs, J_abs

  integer, parameter :: heating_method = 3

  lambda0=lambda

  k = select_absorbing_grain(lambda0,icell, aleat1, heating_method)

  ! Mean intensity
  ! Somme sur differents processeurs
  J_abs=0.0
  do ilambda=1, n_lambda
     J_abs =  J_abs + C_abs_norm(k,ilambda)  * (sum(xJ_abs(icell,ilambda,:)) + J0(icell,lambda))
  enddo ! ilambda
  ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans E_em_1grain
  log_E_abs=log(J_abs*L_packet_th/volume(icell))

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech_1grain_nRE(k,icell,:))

  ! On incremente eventuellement la zone de temperature
  do while((log_E_em_1grain_nRE(k,T_int) < log_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! LIMITE MAX

  ! Save pour prochaine reemission et/ou T finale
  xT_ech_1grain_nRE(k,icell,id)=T_int

  ! Interpolation lineaire entre energies emises pour des
  ! temperatures echantillonees voisines
  T2=T_int
  Temp2=tab_Temp(T2)
  T1=T_int-1
  Temp1=tab_Temp(T1)
  frac=(log_E_abs-log_E_em_1grain_nRE(k,T1))/(log_E_em_1grain_nRE(k,T2)-log_E_em_1grain_nRE(k,T1))
  Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

  !**********************************************************************
  ! Choix de la longeur d'onde de reemission
  ! Dichotomie, la loi de proba est obtenue par interpolation lineaire en T
  frac_T2=(Temp-Temp1)/(Temp2-Temp1)
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  do while((l2-l1) > 1)
     proba=frac_T1*kdB_dT_1grain_nRE_CDF(l,k,T1)+frac_T2*kdB_dT_1grain_nRE_CDF(l,k,T2)
     if(aleat2 > proba) then
        l1=l
     else
        l2=l
     endif
     l=(l1+l2)/2
  enddo
  lambda=l+1

  return

end subroutine im_reemission_qRE

!**********************************************************************

subroutine update_proba_abs_nRE()
  ! Update the various probabilities of absorption for the immediate re-emission
  ! to take into account nRE grains which are actually at equilibrium
  ! C. Pinte
  ! 19/12/2014

  ! TODO :
  ! reemission immadiate : 3 types : RE_LTE, RE_nLTE, nRE_qRE
  ! reemission decalee : nRE

  ! besoin de kappa_abs_tot et de kappa_abs_RE_tot : pour mettre a jour les proba relatives

  ! P_LTE = k_LTE / (k_LTE + k_nLTE + k_qRE) = k_LTE / k_RE
  ! P_LTE_p_nLTE = (k_LTE + k_nLTE) / (k_LTE + k_nLTE + k_qRE) = k_LTE_p_nLTE / k_RE
  ! P_abs_RE =  (k_LTE + k_nLTE + k_qRE) /  (k_LTE + k_nLTE + k_qRE + k_nRE) = k_RE / k_tot
  ! kRE est modifie cat k_qRE l'est

  ! je peux m'en sortir en ne sauvegardant que 1) k_LTE + k_nLTE = k_LTE_p_nLTE et 2) kRE iteration precedente
  ! --> donne k_RE = k_LTE_p_nLTE + k_qRE
  ! puis P_LTE = P_LTE_old * k_RE_old / k_RE_new et P_abs_RE = k_RE_new / k_RE_old * P_abs_RE_old

  real(kind=dp) :: correct, correct_m1,  kappa_abs_RE_new,  kappa_abs_RE_old, delta_kappa_abs_qRE
  integer :: l, lambda, icell
  logical :: lall_grains_eq

  write(*,*) "Setting grains at qRE and updating nRE pobabilities"

  do lambda=1, n_lambda
     do icell=1,n_cells
        delta_kappa_abs_qRE = 0.0_dp
        lall_grains_eq = .true.
        do l=grain_nRE_start,grain_nRE_end
           if (lchange_nRE(l,icell)) then ! 1 grain a change de status a cette iteration
              delta_kappa_abs_qRE =  C_abs_norm(l,lambda) * densite_pouss(l,icell)
           else
              if (.not.l_RE(l,icell)) lall_grains_eq = .false. ! il reste des grains qui ne sont pas a l'equilibre
           endif
        enddo !l

        if (delta_kappa_abs_qRE > tiny_dp) then ! au moins 1 grain a change de status, on met a jour les differentes probabilites
           kappa_abs_RE_old = kappa_abs_RE(icell,lambda)
           kappa_abs_RE_new = kappa_abs_RE_old + delta_kappa_abs_qRE
           kappa_abs_RE(icell,lambda) = kappa_abs_RE_new

           if (kappa_abs_RE_old < tiny_dp) then
              write(*,*) "Oops, opacity of equilibrium grains is 0, cannot perform correction"
              write(*,*) "Something went wrong."
              write(*,*) "Having only grains out of equilibrium is not implemented yet."
              write(*,*) "Cell #", icell, " lambda #", lambda
              write(*,*) "Exiting"
              call exit(1)
           else
              correct = kappa_abs_RE_new / kappa_abs_RE_old ! > 1
              correct_m1 = 1.0_dp/correct ! < 1

              ! Proba d'abs sur un grain en equilibre (pour reemission immediate)
              if (lall_grains_eq) then
                 proba_abs_RE(icell,lambda) = 1.0_dp ! 1 ecrit en dur pour eviter erreur d'arrondis
              else
                 proba_abs_RE(icell,lambda) = proba_abs_RE(icell,lambda) * correct
              endif

              ! Parmis les grains a eq, proba d'absorbe sur un grain a LTE ou sur un grain a LTE ou nLTE
              Proba_abs_RE_LTE(icell,lambda) =  Proba_abs_RE_LTE(icell,lambda) * correct_m1
              Proba_abs_RE_LTE_p_nLTE(icell,lambda) =  Proba_abs_RE_LTE_p_nLTE(icell,lambda) * correct_m1
           endif

        endif
     enddo ! icell
  enddo ! lambda

  write(*,*) "Done"

  return

end subroutine update_proba_abs_nRE

!**********************************************************************

subroutine emission_nRE()
  ! calcule la reemission des grains hors equilibre :
  ! - prob_E_cell
  ! - frac_E_stars a 0
  ! C. Pinte
  ! 05/02/07

  implicit none

  integer :: k, T, lambda, icell
  real(kind=dp) :: Temp, cst_wl, cst_wl_max, wl, delta_wl, fn
  real(kind=dp) :: E_emise, frac_E_abs_nRE, Delta_E
  real(kind=dp), dimension(n_cells) :: E_cell, E_cell_old

  ! proba emission en fct lambda puis pour chaque lambda proba en fct position
  ! proba differentielle en fonction proba emission precedente

  !E_abs_nRE est en nbre de paquets !!!
  frac_E_abs_nRE = E_abs_nRE / nbre_photons_tot
  write(*,*) "Non equilibrium grains have absorbed", real(frac_E_abs_nRE), "of emitted energy"
  write(*,*) "Re-emitting this amount of energy"

  ! Quantite energie a reemettre

  ! energie des paquets en supposant que ce soit le meme nombre de paquets que pour la premiere etape
  !E_paquet = frac_E_abs_nRE

  ! Nouvelle methode : E_paquet reste constante et on modifie le nombre de paquets
  ! + Modification E_paquet pour compenser erreur arrondi du au passage a un entier
  fn = (nbre_photons_tot * frac_E_abs_nRE) / nbre_photons_loop ! nbre de photons "reel"
  nbre_photons_eq_th = int( max( fn, 1.0) ) ! nbre de phtons entier
  ! fn * 1.0 = nbre_photons_eq_th * E_paquet
  E_paquet = fn /  nbre_photons_eq_th

  ! Pour iteration suivante
  E_abs_nRE = 0.0

  ! Repartion spatiale energie
  cst_wl_max = log(huge_real)-1.0e-4

  spectre_emission_cumul(0) = 0.0
  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl = tab_delta_lambda(lambda)*1.e-6
     E_cell=0.0

     ! Emission par cellule des PAHs

     !$omp parallel default(none) &
     !$omp private(k,E_emise,Temp,cst_wl,T,icell) &
     !$omp shared(lambda,wl,delta_wl,E_cell,E_cell_old,tab_lambda,tab_delta_lambda,grain_nRE_start,grain_nRE_end) &
     !$omp shared(n_cells,l_RE, Tdust_1grain_nRE,n_T,C_abs_norm,densite_pouss,volume,tab_Temp,Proba_Tdust) &
     !$omp shared(Emissivite_nRE_old,cst_wl_max,lchange_nRE)
     !$omp do
     do icell=1,n_cells
        E_emise = 0.0
        do k=grain_nRE_start,grain_nRE_end
           if (l_RE(k,icell)) then ! le grain a une temperature
              if (lchange_nRE(k,icell)) then ! la grain passe en qRE a cette iteration : il faut le compter
                 Temp = Tdust_1grain_nRE(k,icell)
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* &
                         volume(icell)/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
                 endif
              endif ! le grain etait en qRE avant, il est traite en re-emission immediate
           else ! densite de proba de Temperature
              do T=1,n_T
                 temp=tab_Temp(T)
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* volume(icell)/ &
                         ((wl**5)*(exp(cst_wl)-1.0)) * Proba_Tdust(T,k,icell) * delta_wl
                 endif !cst_wl
              enddo !T
           endif
        enddo !k
        E_cell(icell) =   E_emise
        ! Recup et mise a jour de l'ancienne emmissivite
        E_cell_old(icell) = Emissivite_nRE_old(icell,lambda)
        if (E_cell(icell) >   E_cell_old(icell)) then
           Emissivite_nRE_old(icell,lambda) = E_emise
        else ! Test en cas de pb d'arrondis
           ! On ne met pas a jour
           E_cell(icell) = E_cell_old(icell)
        endif
     enddo !icell
     !$omp end do
     !$omp end parallel

     prob_E_cell(0,lambda)=0.0
     do icell=1, n_cells
        Delta_E = E_cell(icell) - E_cell_old(icell)
        prob_E_cell(icell,lambda) = prob_E_cell(icell-1,lambda) + Delta_E
     enddo

     spectre_emission_cumul(lambda) =  spectre_emission_cumul(lambda-1) + prob_E_cell(n_cells,lambda)

     ! Les etoiles ont deja emis
     frac_E_stars(lambda) = 0.0

     if (prob_E_cell(n_cells,lambda) > tiny_dp) then
        prob_E_cell(1:n_cells,lambda)=prob_E_cell(1:n_cells,lambda)/prob_E_cell(n_cells,lambda)
     else
        prob_E_cell(1:n_cells,lambda)=0.0
     endif
  enddo ! lambda


  ! Normalisation de la proba d'emission / fct lambda des grains hors equilibre
  do lambda=1,n_lambda
     spectre_emission_cumul(lambda)=spectre_emission_cumul(lambda)/spectre_emission_cumul(n_lambda)
  enddo

  return

end subroutine emission_nRE

!**********************************************************************

subroutine init_emissivite_nRE()
  ! initialise le tableau Emissivite_nRE_old avec l'emissivite a Tmin
  ! C. Pinte
  ! 05/02/07

  implicit none

  integer :: lambda, k, icell
  real(kind=dp) :: E_emise, facteur, cst_wl, wl
  real(kind=dp) :: Temp, cst_wl_max, delta_wl

  cst_wl_max = log(huge_real)-1.0e-4
  temp=tab_Temp(1)

  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl = tab_delta_lambda(lambda)*1.e-6
     cst_wl=cst_th/(Temp*wl)

     if (cst_wl < 100.0) then
        facteur = 1.0_dp/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
     else
        facteur = 0.0_dp
     endif

     ! Emission par cellule des PAHs
     E_emise = 0.0_dp
     do icell=1,n_cells
        do k=grain_nRE_start,grain_nRE_end
           E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* volume(icell) * facteur !* Proba_Tdust = 1 pour Tmin
        enddo !k
        Emissivite_nRE_old(icell,lambda) = E_emise
     enddo !icell

  enddo ! lambda

  return

end subroutine init_emissivite_nRE

!**********************************************************************

subroutine repartition_energie(lambda)
! Calcule la repartition de l'energie emise a la longuer d'onde consideree
! entre l'etoile et les differentes cellules du disque
!  - frac_E_star donne fraction emise par etoile
!  - prob_E_cell donne la proba d'emission cumulée des cellules
!  - E_totale donne energie totale emise (pour calibration des images)
! Utilise une table de temperature pretabulee
! Pour version du code monochromatique avec scattering + em th
! C. Pinte
! 04/02/05
! 01/02/07 : regroupement des cas LTE, nLTE et nRE pour faire une normalisation
! globale en cas de type de grains diffrents

  implicit none

  integer, intent(in) :: lambda

  integer :: k, T, alloc_status
  integer, target :: icell
  integer, pointer :: p_icell
  real(kind=dp) :: Temp, wl, cst_wl, E_star, surface, E_emise, cst_wl_max
  real(kind=dp) :: delta_T
  real(kind=dp), dimension(:), allocatable :: E_cell, E_cell_corrected

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif

  cst_wl_max = log(huge_real)-1.0e-4

  wl = tab_lambda(lambda)*1.e-6

  ! Emission totale des etoiles
  E_star = E_stars(lambda)

  allocate(E_cell(n_cells), E_cell_corrected(n_cells), stat = alloc_status)
  if (alloc_status /= 0) call error("Allocation error in repartition_energie")
  E_cell(:) = 0.0_dp
  E_cell_corrected(:) = 0.0_dp

  ! Cas LTE
  if (lRE_LTE) then
     do icell=1, n_cells
        if (.not.l_dark_zone(icell)) then
           Temp=Tdust(icell)
           if (Temp < tiny_real) then
              ! E_cell(icell) = 0.0_dp
           else
              cst_wl=cst_th/(Temp*wl)
              if (cst_wl < cst_wl_max) then
                 E_cell(icell) = 4.0*kappa_abs_LTE(p_icell,lambda)*kappa_factor(icell)*volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
              endif !cst_wl
           endif ! Temp==0.0
        else ! dark_zone
           E_cell(icell) = 0.0
        endif
     enddo !icell
  endif

  ! Cas nLTE
  if (lRE_nLTE) then
     do icell=1,n_cells
        E_emise = 0.0
        if (.not.l_dark_zone(icell)) then
           do k=grain_RE_nLTE_start,grain_RE_nLTE_end
              Temp=Tdust_1grain(k,icell)
              if (Temp < tiny_real) then
                 !E_emise = E_emise + 0.0
              else
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then

                    E_emise = E_emise +   4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* &
                         volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
                 endif !cst_wl
              endif ! Temp==0.0
           enddo !k
        endif
        E_cell(icell) = E_cell(icell) + E_emise
     enddo !icell
  endif

  ! Cas nRE
  if (lnRE) then

     ! Initialisation du tableau de temperature pour images avec grains
     ! hors equilibre
     if (lmono0) then
        delta_T=exp((1.0_dp/(real(n_T,kind=dp)))*log(T_max/T_min))
        tab_Temp(1)=T_min*sqrt(delta_T)
        do t=2,n_T
           tab_Temp(t)=delta_T*tab_Temp(t-1)
        enddo
     endif

     do icell=1,n_cells
        E_emise = 0.0
        if (.not.l_dark_zone(icell)) then
           do k=grain_nRE_start,grain_nRE_end
              if (l_RE(k,icell)) then ! le grain a une temperature
                 temp=Tdust_1grain_nRE(k,icell)
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* &
                         volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
                 endif !cst_wl
              else ! la grain a une proba de T
                 do T=1,n_T
                    temp=tab_Temp(T)
                    cst_wl=cst_th/(Temp*wl)
                    if (cst_wl < cst_wl_max) then
                       E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* &
                            volume(icell)/((wl**5)*(exp(cst_wl)-1.0)) * Proba_Tdust(T,k,icell)
                    endif !cst_wl
                 enddo !T
              endif ! l_RE
           enddo !k
        endif
        E_cell(icell) =  E_cell(icell) + E_emise
     enddo !icell

  endif

  if (lweight_emission) then
     do icell=1,n_cells
        E_cell_corrected(icell) = E_cell(icell) * weight_proba_emission(icell)
     enddo
  else ! no weight_emission
     E_cell_corrected(:) = E_cell(:)
  endif

  E_disk(lambda) = sum(E_cell)

  if (E_star+E_disk(lambda)+E_ISM(lambda) < tiny_dp) then
     write(*,*) "Error: wl #", lambda, " No energy"
     write(*,*) "Exiting"
     call exit(1)
  endif

  frac_E_stars(lambda)=E_star/(E_star+E_disk(lambda)+E_ISM(lambda))
  frac_E_disk(lambda)=(E_star+E_disk(lambda))/(E_star+E_disk(lambda)+E_ISM(lambda))

  ! Energie totale emise a une distance emise egale a la distance terre-etoile
  ! on chosit cette distance pour calibrer le flux / pi*B(lambda)

  ! Normalisation
  ! Pas de besoin de cst : r_etoile**2 en Au2,
  ! kappa_abs en AU-1, volume en Au**3
  surface=4*pi*(pc_to_AU*distance)**2
  if (l_sym_centrale) then
     E_totale(lambda) = 2.0*pi*hp*c_light**2/surface * (E_star+E_disk(lambda)+E_ISM(lambda)) * real(N_thet)*real(N_phi)
  else
     E_totale(lambda) = 2.0*pi*hp*c_light**2/surface * (E_star+E_disk(lambda)+E_ISM(lambda)) * real(2*N_thet)*real(N_phi)
  endif

  ! Distribution (spatiale) cumulee d'energie
  prob_E_cell(0,lambda)=0.0
  do icell=1, n_cells
     prob_E_cell(icell,lambda) = prob_E_cell(icell-1,lambda) + E_cell_corrected(icell)
  enddo

  ! Normalisation du facteur de correction
  if (lweight_emission) then
     correct_E_emission = correct_E_emission * prob_E_cell(n_cells,lambda) / E_disk(lambda)
  endif

  if (prob_E_cell(n_cells,lambda) > tiny_dp) then
     prob_E_cell(:,lambda)=prob_E_cell(:,lambda)/prob_E_cell(n_cells,lambda)
  else
     prob_E_cell(:,lambda)=0.0
  endif

  deallocate(E_cell,E_cell_corrected)

  return

end subroutine repartition_energie

!**********************************************************************

integer function select_absorbing_grain(lambda,icell, aleat, heating_method) result(k)
  ! This routine will select randomly the absorbing/emitting grain
  ! from the CDF of kabs
  ! Because we cannot store all the CDF for all cells
  ! (n_grains x ncells x n_lambda)
  ! The CDF is recomputed on the fly here
  ! The normalization is saved via kappa_abs, so we do not need to compute
  ! all the CDF

  implicit none

  integer, intent(in) :: lambda, heating_method, icell
  integer :: p_icell
  real, intent(in) :: aleat
  real(kind=dp) :: prob, CDF, norm
  integer :: kstart, kend

  if (lvariable_dust) then
     p_icell = icell
  else
     p_icell = icell_ref
  endif

  ! We scale the random number so that it is between 0 and kappa_sca (= last value of CDF)
  if (heating_method == 1) then
     norm =  kappa_abs_LTE(p_icell,lambda) * kappa_factor(icell) / ( AU_to_cm * mum_to_cm**2 )
     kstart = grain_RE_LTE_start ; kend = grain_RE_LTE_end
  else if (heating_method == 2) then
     norm =  kappa_abs_nLTE(p_icell,lambda) * kappa_factor(icell) / ( AU_to_cm * mum_to_cm**2 )
     kstart = grain_RE_nLTE_start ; kend = grain_RE_nLTE_end
  else
     ! todo : maybe update with an extra variable kappa_abs_qRE
     if (lRE_nLTE) then
        norm =  (kappa_abs_RE(icell, lambda) - &
        (kappa_abs_LTE(p_icell,lambda) + kappa_abs_nLTE(p_icell,lambda)) * kappa_factor(icell)) &
             / ( AU_to_cm * mum_to_cm**2 )
     else
        norm =  (kappa_abs_RE(icell, lambda) -  kappa_abs_LTE(p_icell,lambda) * kappa_factor(icell)) &
             / ( AU_to_cm * mum_to_cm**2 )
     endif
     kstart = grain_nRE_start ; kend = grain_nRE_end
  endif

  if (heating_method <= 2) then
     if (aleat < 0.5) then ! We start from first grain
        prob = aleat * norm
        CDF = 0.0
        do k=kstart, kend
           CDF = CDF + C_abs(k,lambda) * densite_pouss(k,icell)
           if (CDF > prob) exit
        enddo
     else ! We start from the end of the grain size distribution
        prob = (1.0-aleat) * norm
        CDF = 0.0
        do k=kend, kstart, -1
           CDF = CDF + C_abs(k,lambda) * densite_pouss(k,icell)
           if (CDF > prob) exit
        enddo
     endif
  else ! Same thing but with lRE to only include dust grains that are at qRE
     if (aleat < 0.5) then ! We start from first grain
        prob = aleat * norm
        CDF = 0.0
        do k=kstart, kend
           if (l_RE(k,icell)) CDF = CDF + C_abs(k,lambda) * densite_pouss(k,icell)
           if (CDF > prob) exit
        enddo
     else ! We start from the end of the grain size distribution
        prob = (1.0-aleat) * norm
        CDF = 0.0
        do k=kend, kstart, -1
           if (l_RE(k,icell)) CDF = CDF + C_abs(k,lambda) * densite_pouss(k,icell)
           if (CDF > prob) exit
        enddo
     endif
  endif
  !if (k < kstart) k=kstart
  !if (k > kend)   k=kend


  return

end function select_absorbing_grain

!**********************************************************************

subroutine select_cellule(lambda,aleat, icell)
  ! Sélection de la cellule qui va émettre le photon
  ! C. Pinte
  ! 04/02/05
  ! Modif 3D 10/06/05

  implicit none

  integer, intent(in) :: lambda
  real, intent(in) :: aleat
  integer, intent(out) :: icell
  integer :: k, kmin, kmax

  ! Dichotomie
  kmin=0
  kmax=n_cells
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_E_cell(k,lambda) < aleat) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   icell=kmax

   return

end subroutine select_cellule

!***********************************************************

subroutine define_proba_weight_emission(lambda)
  ! Augmente le poids des cellules pres de la surface
  ! par exp(-tau)
  ! Le poids est applique par weight_repartion_energie
  ! C. Pinte
  ! 19/11/08

  integer, intent(in) :: lambda

!--  use optical_depth, only : optical_length_tot
!--
!--  implicit none
!--
!--
!--
!--  real, dimension(n_cells) :: tau_min
!--  real(kind=dp), dimension(4) :: Stokes
!--  real(kind=dp) :: x0, y0, z0, u0, v0, w0, angle, lmin, lmax
!--  real :: tau
!--  integer :: i, j, n, id, icell
!--  integer, parameter :: nbre_angle = 101
!--
!--  tau_min(:) = 1.e30 ;
!--
!--  do icell=1,n_cells
!--     do n=1,nbre_angle
!--        id=1
!--        ! position et direction vol
!--        angle= pi * real(n)/real(nbre_angle+1)! entre 0 et pi
!--        i = cell_map_i(icell)
!--        j = cell_map_j(icell)
!--        x0=1.00001*r_lim(i-1) ! cellule 1 traitee a part
!--        y0=0.0
!--        z0=0.99999*z_lim(i,j+1)
!--        u0=cos(angle)
!--        v0=0.0
!--        w0=sin(angle)
!--
!--        Stokes(:) = 0.0_dp ;
!--        call optical_length_tot(id,lambda,Stokes,icell,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)
!--        if (tau < tau_min(icell)) tau_min(icell) = tau
!--
!--        x0 = 0.99999*r_lim(i)
!--        call optical_length_tot(id,lambda,Stokes,icell,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)
!--        if (tau < tau_min(icell)) tau_min(icell) = tau
!--
!--     enddo
!--  enddo ! icell
!--
!--
!--  weight_proba_emission(1:n_cells) =  exp(-tau_min(:))
!--
!--  ! correct_E_emission sera normalise dans repartition energie
!--  correct_E_emission(1:n_cells) = 1.0_dp / weight_proba_emission(1:n_cells)
!--
!--  return

end subroutine define_proba_weight_emission

!***********************************************************

subroutine allocate_weight_proba_emission(Nc)

  integer, intent(in) :: Nc

  integer ::  alloc_status

  allocate(weight_proba_emission(Nc), correct_E_emission(Nc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error prob_E_cell 3')
  weight_proba_emission = 1.0 ;
  correct_E_emission = 1.0 ;

  return

end subroutine allocate_weight_proba_emission

!***********************************************************

subroutine reset_temperature()

  use radiation_field, only : xKJ_abs, xJ_abs

  if (lRE_LTE) then
     xKJ_abs(:,:) = 0.0_dp
  endif
  if (lRE_nLTE .or. lnRE) xJ_abs(:,:,:) = 0.0_dp
  xT_ech = 2
  Tdust = 1.0

  return

end subroutine reset_temperature

end module thermal_emission
