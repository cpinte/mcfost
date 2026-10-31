module thermal_emission

  use parameters
  use constants
  use grains
  use dust_prop
  !$ use omp_lib
  use temperature
  use utils
  use PAH
  use grid
  use stars, only : E_ISM, E_stars
  use input
  use density
  use cylindrical_grid

  implicit none
  save

  public :: prob_E_cell, E_totale, nbre_reemission, frac_E_stars, frac_E_disk, L_packet_th, E_abs_nRE, E_abs_nRE_id, &
       correct_E_emission, &
       DensE, DensE_m1, Dcoeff

  public :: allocate_temperature, deallocate_temperature_calculation, realloc_emitting_fractions, &
       allocate_thermal_emission, deallocate_thermal_emission, allocate_weight_proba_emission, &
       define_proba_weight_emission, emission_nre, im_reemission_lte, im_reemission_nlte, &
       im_reemission_qre, init_emissivite_nre, init_reemission, select_wl_em, &
       select_cellule, temp_finale, temp_finale_nlte, temp_nre, update_proba_abs_nre, &
       repartition_wl_em, repartition_energie, allocate_thermal_energy, set_min_temperature, &
       reset_temperature, reset_thermal_state

  private

  ! Choice of emission cell for monochromatic case
  real(kind=dp), dimension(:,:), allocatable :: prob_E_cell !n_lambda,0:n_rad*nz
  real(kind=dp), dimension(:), allocatable :: frac_E_stars, frac_E_disk, E_totale !n_lambda

  ! fraction of energy re-emitted over star energy
  ! (average Planck opacity * coeff)
  real(kind=dp), dimension(:,:), allocatable :: log_Qcool_minus_extra_heating ! 0:n_T, n_cells
  real(kind=dp), dimension(:,:), allocatable :: log_E_em_1grain  !n_grains,0:n_T
  real(kind=dp), dimension(:,:), allocatable :: E_em_1grain_nRE, log_E_em_1grain_nRE !n_grains,0:n_T
  ! enthalpy of non equilibrium grains, at the temperature bin centers [Hz] and edges [J]
  real(kind=dp), dimension(:,:), allocatable :: nu_bin_nRE, U_lim_nRE ! n_T and 0:n_T, n_grains

  ! Cumulative probability in lambda of dust emissivity
  ! with temperature correction (dp/dT)
  ! (Bjorkman & Wood 2001, A&A 554-615 -- eq 9)
  real(kind=dp), dimension(:,:,:), allocatable :: kdB_dT_CDF ! 0:n_T,n_cells,n_lambda
  real(kind=dp), dimension(:,:,:), allocatable :: kdB_dT_1grain_LTE_CDF, kdB_dT_1grain_nLTE_CDF, kdB_dT_1grain_nRE_CDF ! n_grains,0:n_T,n_lambda

  integer, dimension(:,:), allocatable :: xT_ech ! n_cells, id
  integer, dimension(:,:,:), allocatable :: xT_ech_1grain, xT_ech_1grain_nRE ! n_grains, n_cells, id

  ! Cumulative probability in lambda to emit as a blackbody
  real(kind=dp), dimension(:), allocatable :: spectre_emission_cumul !(0:n_lambda)

  ! emissivity in some units (missing a constant but working in relative)
  real(kind=dp), dimension(:,:), allocatable :: Emissivite_nRE_old ! n_lambda, n_rad, nz, n_az

  real(kind=dp), dimension(:,:), allocatable :: nbre_reemission ! n_cells, id

  real(kind=dp) :: L_packet_th, E_abs_nRE
  real(kind=dp), dimension(:), allocatable :: E_abs_nRE_id ! energy absorbed by nRE grains, per thread

  ! Emission bias towards the disk surface
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

  allocate(xT_ech(Nc,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error xT_ech')
  xT_ech = 2

  if (lreemission_stats) then
     allocate(nbre_reemission(Nc,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error nbre_reemission')
     nbre_reemission = 0.0
  endif

  mem_size = (1.0 * p_Nc) * n_T * n_lambda * 8. / 1024.**3 ! double precision
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
     mem_size = (1.0 * (grain_RE_nLTE_end-grain_RE_nLTE_start+2)) * p_Nc * n_lambda * 8. / 1024.**3 ! double precision
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

     allocate(nu_bin_nRE(n_T,grain_nRE_start:grain_nRE_end), U_lim_nRE(0:n_T,grain_nRE_start:grain_nRE_end), &
          stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error nu_bin_nRE')
     nu_bin_nRE = 0.0_dp ; U_lim_nRE = 0.0_dp

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

     allocate(E_abs_nRE_id(nb_proc),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error E_abs_nRE_id')
     E_abs_nRE_id = 0.0_dp ; E_abs_nRE = 0.0_dp

     allocate(kdB_dT_1grain_nRE_CDF(n_lambda,grain_nRE_start:grain_nRE_end,n_T),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error kdB_dT_1grain_nRE_CDF')
     kdB_dT_1grain_nRE_CDF=0.0

     if (lRE_nlTE) then
        allocate(Tdust_1grain_old(grain_RE_nLTE_start:grain_RE_nLTE_end,Nc),stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error Tdust_1grain_old')
        Tdust_1grain_old=0.
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
     deallocate(kdB_dT_1grain_nRE_CDF,E_em_1grain_nRE,log_E_em_1grain_nRE,xT_ech_1grain_nRE,nu_bin_nRE,U_lim_nRE)
     deallocate(Tdust_1grain_nRE_old,Emissivite_nRE_old,Tpeak_old,E_abs_nRE_id)
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
! For step 1
! Corrects the emission distribution function by accounting for
! the initial disk heating (for select_wl_em)
! Corrects the packet energy for step 1
! Must come after repartition_energie_etoiles and repartition_energie
! C. Pinte
! 18/05/06

  implicit none

  integer :: lambda
  real(kind=dp) :: E_star_tot, E_disk_tot, E_ISM_tot, delta_wl, L_tot

  if (lTemp) then
     spectre_emission_cumul(0) = 0.0
     ! Emission distribution function
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

  ! Packet energy for step 1
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
  L_packet_th = L_tot/n_photons_total

  return

end subroutine repartition_wl_em

!***************************************************

subroutine select_wl_em(rand,lambda)
! Select wavelength from the previously computed blackbody
! Dichotomie
! lambda is the wavelength index
! Uses the results of repartition_wl_em
! C. Pinte

  implicit none

  real, intent(in) :: rand
  integer, intent(out) :: lambda ! wavelength index

  integer :: k, kmin, kmax

  kmin=0
  kmax=n_lambda

  k=(kmin + kmax)/2

  ! Dichotomie
  do while (spectre_emission_cumul(k) /= rand)
     if (spectre_emission_cumul(k) < rand) then
        kmin = k
     else
        kmax = k
     endif

     k = (kmin + kmax)/2

     if ((kmax-kmin) <= 1) then
        exit ! Exit the while loop
     endif
  enddo   ! while

  lambda=kmax

end subroutine select_wl_em

!***************************************************

subroutine init_reemission(lheating,dudt)
! Pre-tabulates Planck opacities as a function of T and position
! + the k dB/dT terms for wavelength sampling of the emitted photon
! with temperature correction
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
  real(kind=dp), dimension(0:n_T) :: T_lim

  real(kind=dp) :: Qcool, extra_heating, u_o_dt, Qcool_minus_extra_heating, Qcool0

  write(*,'(a36, $)') " Initializing thermal properties ..."

  cst_E = 2.0*hp*c_light**2 * four_pi

  call init_tab_temp()

  do t=1,n_T
     ! Calculation of the blackbody and its derivative
     ! Within a constant (not the same in the two cases!!)

     Temp=tab_Temp(t)
     cst=thermal_const/Temp
     do lambda=1, n_lambda
        ! wavelength in meters
        wl = tab_lambda(lambda)*1.e-6
        delta_wl=tab_delta_lambda(lambda)*1.e-6
        cst_wl=cst/wl
        if (cst_wl < 500.0) then
           coeff_exp=exp(cst_wl)
           ! The formulas take into account the delta_wl of the integration
           B(lambda,T) = 1.0/((wl**5)*(coeff_exp-1.0))*delta_wl
           dB_dT(lambda,T) = B(lambda,T) * cst_wl*coeff_exp/(coeff_exp-1.0) !/Temp * temp because of dT but it changes nothing in practice
        else
           B(lambda,T)=0.0
           dB_dT(lambda,T)=0.0
        endif
     enddo !lambda
  enddo ! T

  ! multiplied by opacity (abs seule) massique
  ! Todo : this loop is OK in 2D but takes ~ 5sec for 0.5 million cells in 3D
  !$omp parallel default(none) &
  !$omp private(id,icell,T,lambda,integ, Qcool,Qcool0,extra_heating,Qcool_minus_extra_heating,Temp,u_o_dt) &
  !$omp shared(cst_E,kappa_abs_LTE,kappa_factor,volume,B,lextra_heating,xT_ech,log_Qcool_minus_extra_heating,J0) &
  !$omp shared(n_T,n_cells,p_n_cells,n_lambda,tab_Temp,ldudt_implicit,ufac_implicit,dudt,lRE_nLTE,lvariable_dust,icell1)
  id = 1 ! for sequential code
  !$ id = omp_get_thread_num() + 1

  !$omp do
  do icell=1,p_n_cells
     do T=1, n_T
        Temp = tab_Temp(T)
        integ=0.0
        do lambda=1, n_lambda
           ! kappa en Au-1    \
           ! volume in AU3  --> no constant needed to have E_em in SI
           ! B * cst_E en SI = W.m-2.sr-1 (.m-1 * m) because delta_wl included
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
           ! Energy coming from equilibrium with cloud à T_min
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
           integ3(0) = 0.0
           do lambda=1, n_lambda ! from the first wavelength bin, as kdB_dT_CDF (B13)
              ! No need for constant - nor volume (normalization to 1), nor density
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
              ! No need for constant - ni du volume, nor kappa factor (normalisation a 1)
              integ3(lambda) = integ3(lambda-1) + kappa_abs_LTE(icell,lambda)  * dB_dT(lambda,T)
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
     ! multiplied by opacity (abs seule) mass opacity of a single grain size
     do T=1, n_T
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ=0.0
           do lambda=1, n_lambda
              ! WARNING:  density and volume are missing compared to the LTE case !!!
              integ = integ + C_abs_norm(k,lambda) * B(lambda,T)
           enddo !lambda
           ! The appropriate coefficient
           if (integ > tiny_dp) then
              log_E_em_1grain(k,T)=log(integ*cst_E)
           else
              log_E_em_1grain(k,T)=-1000.
           endif
        enddo ! k

        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ3(0) = 0.0
           do lambda=1, n_lambda ! from the first wavelength bin, as kdB_dT_CDF (B13)
              ! No need for constant - nor volume (normalization to 1), nor density
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
           integ3(0) = 0.0
           do lambda=1, n_lambda ! from the first wavelength bin, as kdB_dT_CDF (B13)
              ! No need for constant - nor volume (normalization to 1), nor density
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

  if (lnRE) then
     ! Enthalpy of the temperature bins (after reading the heat capacity files), used by
     ! Temp_nRE at each iteration
     T_lim(0) = T_min
     T_lim(1:n_T) = tab_Temp(:) * sqrt(exp((1.0_dp/(real(n_T,kind=dp)))*log(T_max/T_min)))
     do k=grain_nRE_start,grain_nRE_end
        nu_bin_nRE(:,k) = enthalpy(real(tab_Temp,kind=dp),k)/hp ! enthalpy at bin centers [Hz]
        U_lim_nRE(:,k) = enthalpy(T_lim(:),k) ! enthalpy of bin boundaries [J]
     enddo
  endif

  write(*,*) "Done"
  return

end subroutine init_reemission


!********************************************************************

subroutine Temp_LTE(id, icell, Ti, Temp, frac)

  use radiation_field, only : xKJ_abs

  integer, intent(in) :: id, icell
  integer, intent(out) :: Ti
  real, intent(out) :: Temp
  real(kind=dp), intent(out) :: frac

  real(kind=dp) :: Qheat, log_Qheat
  integer :: p_icell

  if (lvariable_dust) then
     p_icell = icell
  else
     p_icell = icell1
  endif

  if (id==0) then
     Qheat=sum(xKJ_abs(icell,:)) * L_packet_th / volume(icell)
  else
     Qheat=xKJ_abs(icell,id) * nb_proc * L_packet_th / volume(icell)  ! does not include kappa_factor, same as log_Qcool
  endif

  if (Qheat < tiny_dp) then
     Temp = T_min ; Ti = 2
  else
     log_Qheat = log(Qheat)

     if (log_Qheat <  log_Qcool_minus_extra_heating(1,p_icell)) then
        Temp = T_min ; Ti = 2
     else
        ! Sampled temperature just above the cell temperature
        if (id==0) then
           Ti = minval(xT_ech(icell,:))
        else
           Ti = xT_ech(icell,id)
        endif

        ! Eventually increment the temperature zone
        do while((log_Qcool_minus_extra_heating(Ti,p_icell) < log_Qheat).and.(Ti < n_T))
           Ti=Ti+1
        enddo  ! LIMITE MAX

        ! Linear interpolation between energies emitted for
        ! neighboring sampled temperatures
        frac=(log_Qheat-log_Qcool_minus_extra_heating(Ti-1,p_icell)) / &
             (log_Qcool_minus_extra_heating(Ti,p_icell)-log_Qcool_minus_extra_heating(Ti-1,p_icell))
        Temp=exp(log(tab_Temp(Ti))*frac+log(tab_Temp(Ti-1))*(1.0-frac))
     endif
  endif

  ! Save for next re-emission and/or final T
  if (id>0) xT_ech(icell,id) = Ti

  return

end subroutine Temp_LTE

!********************************************************************

subroutine im_reemission_LTE(id,icell,p_icell,rand1,rand2,lambda)
! Calculation of the cell temperature and storage of received energy + T
! Re-emission of a photon at the right wavelength

  implicit none

  integer, intent(in) :: id, icell, p_icell
  real, intent(in) ::  rand1, rand2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, Ti, T1, T2, k, heating_method, p_k
  real :: Temp
  real(kind=dp) :: frac_T1, frac_T2, proba

  if (lreemission_stats) nbre_reemission(icell,id) = nbre_reemission(icell,id) + 1.0_dp

  call Temp_LTE(id, icell, Ti, Temp, frac_T2)

  !**********************************************************************
  ! Wavelength choice for re-emission
  ! Bisection, the probability law is obtained by linear interpolation in T
  T2 = Ti ; T1 = Ti-1
  !Frac_T2=(Temp-tab_Temp(T1))/(tab_temp(T2)-tab_Temp(T1))
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  if (low_mem_th_emission) then
     ! Select absorbing grain
     heating_method = 1 ! LTE
     k = select_absorbing_grain(lambda,icell, rand1, heating_method)

     ! Select wavelength
     do while((l2-l1) > 1)
        proba=frac_T1*kdB_dT_1grain_LTE_CDF(l,k,T1)+frac_T2*kdB_dT_1grain_LTE_CDF(l,k,T2)
        if(rand2 > proba) then
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
        if(rand2 > proba) then
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

subroutine im_reemission_NLTE(id,icell,p_icell,rand1,rand2,lambda)
! Calculation of the cell temperature
! Re-emission of a photon at the right wavelength

  use radiation_field, only : xJ_abs, J0

  implicit none

  integer, intent(in) :: id,icell, p_icell
  real, intent(in) :: rand1, rand2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, kmin, kmax, lambda0, ilambda, heating_method, p_k
  real(kind=dp) :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_E_abs, J_abs

  lambda0=lambda

  ! Selection of the grain that absorbs the photon
  if (low_mem_th_emission_nLTE) then
     ! Select absorbing grain
     heating_method = 2 ! nLTE
     k = select_absorbing_grain(lambda0,icell, rand1, heating_method)
  else
     ! kabs_nLTE_CDF(grain_RE_nLTE_start-1) = 0: first grain k with CDF(k) >= rand1
     kmin=grain_RE_nLTE_start-1
     kmax=grain_RE_nLTE_end
     k=(kmin+kmax)/2

     do while((kmax-kmin) > 1)
        if (kabs_nLTE_CDF(k,icell,lambda0) < rand1) then
           kmin = k
        else
           kmax = k
        endif
        k = (kmin + kmax)/2
     enddo   ! while
     k=kmax
  endif

  ! Mean intensity: estimate of this thread x nb_proc, as for LTE grains (Temp_LTE), instead
  ! of summing the other threads at each event; the temperature also follows this thread
  J_abs=0.0
  do ilambda=1, n_lambda
     J_abs =  J_abs + C_abs_norm(k,ilambda)  * (xJ_abs(icell,ilambda,id) * nb_proc + J0(icell,ilambda))
  enddo ! lambda
  ! WARNING: must divide by dust_density_o_n_grains as it is not taken into account in E_em_1grain
  log_E_abs=log(J_abs*L_packet_th/volume(icell) )

  ! Sampled temperature just above the cell temperature
  T_int=xT_ech_1grain(k,icell,id)

  ! Eventually increment the temperature zone
  do while((log_E_em_1grain(k,T_int) < log_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! LIMITE MAX

  ! Save for next re-emission and/or final T
  xT_ech_1grain(k,icell,id)=T_int

  ! Linear interpolation between energies emitted for
  ! neighboring sampled temperatures
  T2=T_int
  Temp2=tab_Temp(T2)
  T1=T_int-1
  Temp1=tab_Temp(T1)
  frac=(log_E_abs-log_E_em_1grain(k,T1))/(log_E_em_1grain(k,T2)-log_E_em_1grain(k,T1))
  Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))


  !**********************************************************************
  ! Wavelength choice for re-emission
  ! Bisection, the probability law is obtained by linear interpolation in T
  frac_T2=(Temp-Temp1)/(Temp2-Temp1)
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  do while((l2-l1) > 1)
     proba=frac_T1*kdB_dT_1grain_nLTE_CDF(l,k,T1)+frac_T2*kdB_dT_1grain_nLTE_CDF(l,k,T2)
     if(rand2 > proba) then
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
! Calculates the final disk temperature dans chaque cell
! Case of "continuous absorption"
! C. Pinte
! 24/01/05

  implicit none

  real :: Temp
  integer :: Ti, icell, id

  real(dp) :: frac

  id=0 ! as we will sum the results from all threads

  !$omp parallel &
  !$omp default(none) &
  !$omp private(icell,Temp,Ti,frac) &
  !$omp shared(id,Tdust,n_cells)
  !$omp do schedule(dynamic,10)
  do icell=1,n_cells
     call Temp_LTE(id,icell, Ti, Temp, frac)
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
! Calculates the final disk temperature dans chaque cell
! Case of "continuous absorption"
! C. Pinte
! 24/01/05

  use radiation_field, only : xJ_abs, J0

  implicit none

  integer :: T_int, T1, T2,icell
  real(kind=dp) :: Temp, Temp1, Temp2, frac, log_E_abs, J_absorbe

  integer :: k, lambda, p_k

  ! Calculation of the cell temperature and storage of received energy + T
  ! Using previous temperature

  !$omp parallel &
  !$omp default(none) &
  !$omp private(J_absorbe,log_E_abs,T_int,T1,T2,Temp1,Temp2,Temp,frac,icell,p_k) &
  !$omp shared(L_packet_th,Tdust,tab_Temp,n_cells,n_lambda,kappa_abs_LTE,lvariable_dust) &
  !$omp shared(xJ_abs,dust_density_o_n_grains,n_grains,Tdust_1grain, xT_ech_1grain,log_E_em_1grain) &
  !$omp shared(C_abs_norm,volume, grain_RE_nLTE_start, grain_RE_nLTE_end, n_T, T_min, J0, grain)
  !$omp do schedule(dynamic,10)
  do icell=1,n_cells
     do k=grain_RE_nLTE_start, grain_RE_nLTE_end
        p_k = merge(k,grain(k)%zone,lvariable_dust)
        if (dust_density_o_n_grains(p_k,icell) > tiny_dp) then
           J_absorbe=0.0
           do lambda=1, n_lambda
              J_absorbe =  J_absorbe + C_abs_norm(k,lambda)  * (sum(xJ_abs(icell,lambda,:)) + J0(icell,lambda))
           enddo ! lambda

           ! WARNING: must divide by dust_density_o_n_grains as it is not taken into account in E_em_1grain
           J_absorbe = J_absorbe*L_packet_th/volume(icell)
           if (J_absorbe < tiny_dp) then
              Tdust_1grain(k,icell) = T_min
           else
              log_E_abs=log(J_absorbe)

              if (log_E_abs <  log_E_em_1grain(k,1)) then
                 Tdust_1grain(k,icell) = T_min
              else
                 ! Sampled temperature just above the cell temperature: lowest over the threads,
                 ! each thread follows its own estimate of the field (as Temp_LTE with id = 0)
                 T_int=max(minval(xT_ech_1grain(k,icell,:)),2)

                 ! Eventually increment the temperature zone
                 do while((log_E_em_1grain(k,T_int) < log_E_abs).and.(T_int < n_T))
                    T_int=T_int+1
                 enddo

                 ! Linear interpolation between energies emitted for
                 ! neighboring sampled temperatures
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
  real(kind=dp), dimension(0:n_T) :: U_lim

  real(kind=dp), dimension(n_lambda) :: C_abs_norm_o_dnu, tab_nu
  real(kind=dp), dimension(0:n_lambda) :: nu_edge

  real(kind=dp), dimension(:,:), allocatable :: kJnu, lambda_Jlambda, kJnu_o_nu, N_above
  integer, dimension(:,:), allocatable :: k_edge
  real(kind=dp), dimension(:,:), allocatable :: dnu_edge

  real(kind=dp) :: delta_T, N_lo, N_hi, t_cool, t_abs, mean_abs_E, kTu, mean_abs_nu, cst_t_cool, Int_k_lambda_Jlambda

  integer :: l, T, T1, T2, lambda, alloc_status, id, T_int, icell, p_l, k

  integer, parameter :: n_cooling_time = 100

  real(kind=dp), dimension(n_cooling_time+1) :: en_lim
  real(kind=dp), dimension(n_cooling_time) :: en, delta_en, Cabs

  real :: E_max
  real(kind=dp) :: Temp1, Temp2, frac, log_E_abs, maxP
  real(kind=dp) :: wl, wl_mum, J_lambda
  integer :: Tpeak

  write(*,*) "Calculating temperature probabilities of non equilibrium grains ..."
  if (lforce_PAH_equilibrium) write(*,*) "Forcing equilibrium " ;

  allocate(A(n_T,n_T,nb_proc), B(n_T,n_T,nb_proc), X(n_T,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error A')
  A=0.0 ; B=0.0 ;

  allocate(kJnu(n_lambda,nb_proc), lambda_Jlambda(n_lambda,nb_proc), kJnu_o_nu(n_lambda,nb_proc), &
       N_above(0:n_lambda,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error kJnu')
  kJnu = 0.0 ;  lambda_Jlambda=0.0 ; kJnu_o_nu = 0.0 ; N_above = 0.0

  allocate(k_edge(n_T,n_T), dnu_edge(n_T,n_T), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error k_edge')
  k_edge = 0 ; dnu_edge = 0.0_dp

  tab_nu = c_light/(tab_lambda*1.0e-6)
  ! frequency edges of the wavelength bins (decreasing)
  nu_edge(0) = c_light/(tab_lambda_inf(1)*1.0e-6)
  nu_edge(1:n_lambda) = c_light/(tab_lambda_sup(:)*1.0e-6)

  ! l est l'indice de taille de grains
  do l=grain_nRE_start, grain_nRE_end
     ! bin centers [Hz]
     nu_bin(:) = nu_bin_nRE(:,l) ! enthalpy at bin centers [Hz]

     ! bin widths [Hz]
     U_lim(:) = U_lim_nRE(:,l) ! enthalpy of bin boundaries [J]
     delta_nu_bin(:) = (U_lim(1:n_T)-U_lim(0:n_T-1))/hp  ! bin width [Hz]

     ! compute transition matrix  A   [Hz]

     !! cooling terms :  select one-off upper diagonal, they can be calculated in advance
     ! Afi with final < inital
     A(:,:,:)=0.
     do T=2,n_T
        A(T-1,T,:) = E_em_1grain_nRE(l,T) / (nu_bin(T)-nu_bin(T-1)) ! OK, missing a 1/h factor compared to Krugel p265
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
     ! t_abs and t_cool are both computed with C_abs_norm (C_abs x 1e4 AU_to_cm in m2):
     ! the common factor cancels in their comparison
     E_max = real(maxval(U_lim))
     en_lim = spanl(1.0e-10*E_max,E_max, n_cooling_time+1)
     do T=1,n_cooling_time
        en(T) = 0.5 * (en_lim(T+1) + en_lim(T))
        delta_en(T) = en_lim(T+1) - en_lim(T)
        Cabs(T) = interp(real(C_abs_norm(l,:),kind=dp),real(tab_lambda,kind=dp),hp*c_light/en(T) * 1.0e6_dp) ! lambda [mum]
     enddo

     cst_t_cool =(real(hp,kind=dp)**3 * real(c_light,kind=dp)**2) / (8*pi)

     ! Heating terms: the photon frequencies U_lim(T)/h - nu_bin(T2) at the edges of the
     ! transitions T2 -> T depend only on the grain: their wavelength bins are found once here
     do T2=1,n_T-1
        do T=T2,n_T-1
           call photon_edge(U_lim(T)/hp - nu_bin(T2), nu_edge, k_edge(T,T2), dnu_edge(T,T2))
        enddo
     enddo

     ! Loop over cells

     !$omp parallel &
     !$omp default(none) &
     !$omp private(Int_k_lambda_Jlambda, lambda, wl, wl_mum, T, T2) &
     !$omp private(N_lo,N_hi,id,t_cool,t_abs,mean_abs_E,mean_abs_nu,kTu,k) &
     !$omp private(frac,T1,Temp1,Temp2,T_int,log_E_abs,icell,p_l,J_lambda) &
     !$omp shared(l,kJnu, lambda_Jlambda, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium) &
     !$omp shared(n_cells, C_abs_norm_o_dnu, xJ_abs, J0, L_packet_th, volume, n_T, star) &
     !$omp shared(tab_nu, nu_edge, n_lambda, tab_delta_lambda, tab_lambda,en,delta_en,Cabs,kJnu_o_nu,N_above,k_edge,dnu_edge) &
     !$omp shared(delta_nu_bin,Proba_Tdust, A,B,X,nu_bin,tab_Temp,T_min,T_max,lbenchmark_SHG,lMathis_field,Mathis_field) &
     !$omp shared(Tdust_1grain_nRE,log_E_em_1grain_nRE,cst_t_cool,C_abs_norm,l_RE,r_grid,z_grid,n_grains) &
     !$omp shared(dust_density_o_n_grains,l_dark_zone,Tdust,lchange_nRE,lvariable_dust,grain)

     id = 1 ! pour code sequentiel
     ! low granularity because calculation time depends strongly on cells
     p_l = merge(l, grain(l)%zone, lvariable_dust)
     !$omp do schedule(dynamic,1)
     do icell=n_cells, 1, -1
        !$ id = omp_get_thread_num() + 1
        if (l_dark_zone(icell)) then
           lchange_nRE(l,icell) = .not.l_RE(l,icell)
           l_RE(l,icell) = .true.
        else
           if (dust_density_o_n_grains(p_l,icell) > tiny_dp) then
              ! Champ de radiation
              Int_k_lambda_Jlambda=0.0
              do lambda=1, n_lambda
                 ! conversion of lambda and delta_lambda from micron -> m
                 lambda_Jlambda(lambda,id) = (sum(xJ_abs(icell,lambda,:)) + J0(icell,lambda))
                 Int_k_lambda_Jlambda = Int_k_lambda_Jlambda + C_abs_norm(l,lambda) * lambda_Jlambda(lambda,id)
                 kJnu(lambda,id) =   C_abs_norm_o_dnu(lambda)  * lambda_Jlambda(lambda,id)
              enddo ! lambda

              lambda_Jlambda(:,id) =  lambda_Jlambda(:,id)*L_packet_th * (1.0/volume(icell))
              Int_k_lambda_Jlambda = Int_k_lambda_Jlambda*L_packet_th * (1.0/volume(icell))
              kJnu(:,id) = kJnu(:,id)*L_packet_th * (1.0/volume(icell))

              if (lbenchmark_SHG) then
                 ! Imposed radiation field of the SHG benchmark (Camps et al. 2015) instead of the MC
                 ! field: lambda_Jlambda is 4pi J_lambda Delta_lambda [W.m-2], as the MC estimate
                 ! above and the emission E_em_1grain_nRE (cst_E includes 4pi)
                 Int_k_lambda_Jlambda = 0.0
                 do lambda=1, n_lambda
                    wl_mum = tab_lambda(lambda)
                    wl = wl_mum * 1e-6
                    if (lMathis_field) then
                       ! Mathis et al. (1983) ISRF, Eq. 11 of Camps et al. 2015 [W.m-3.sr-1] x U
                       if (wl_mum < 0.0912) then
                          J_lambda = 0.0_dp
                       else if (wl_mum < 0.110) then
                          J_lambda = 3069 * wl_mum**3.4172
                       else if (wl_mum <  0.134) then
                          J_lambda = 1.627
                       else if (wl_mum <  0.250) then
                          J_lambda = 0.0566 * wl_mum**(-1.6678)
                       else
                          J_lambda = 1e-14 * Blambda_dp(wl,7500.) + 1e-13 * Blambda_dp(wl,4000.) &
                               + 4e-13 * Blambda_dp(wl,3000.)
                       endif
                       J_lambda = J_lambda * Mathis_field
                    else
                       ! star diluted at the distance of the cell: W = (R_star/2r)^2
                       J_lambda = star(1)%r**2 / (4.0_dp * (r_grid(icell)**2 + z_grid(icell)**2)) * Blambda_dp(wl,star(1)%T)
                    endif
                    lambda_Jlambda(lambda,id) = four_pi * J_lambda * tab_delta_lambda(lambda) * 1e-6

                    Int_k_lambda_Jlambda = Int_k_lambda_Jlambda + C_abs_norm(l,lambda) * lambda_Jlambda(lambda,id)
                    kJnu(lambda,id) =   C_abs_norm_o_dnu(lambda)  * lambda_Jlambda(lambda,id)
                 enddo ! lambda
              endif ! lbenchmark_SHG

              ! decide whether we really need to use this model, instead of calculating the equilibrium temperature

              ! time interval for photon absorption: lambda_Jlambda = 4pi J_lambda Delta_lambda,
              ! absorption rate = sum C_abs 4pi J_lambda Delta_lambda / (hc/lambda)
              t_abs = sum(C_abs_norm(l,:)*lambda_Jlambda(:,id)*tab_lambda(:)) * 1.0e-6/c_light
              if (t_abs > tiny_dp) then
                 t_abs = hp/t_abs
              else
                 t_abs = huge_real
              endif

              ! mean absorbed photon energy
              mean_abs_E = t_abs * Int_k_lambda_Jlambda ! D01 eq 46
              mean_abs_nu = mean_abs_E/hp

              ! Fast energy to temperature inversion
              ! Take energy just above average energy -> t_cool slightly too short
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
                    t_cool = cst_t_cool*mean_abs_E/t_cool
                 else
                    t_cool = huge_real
                 endif
              endif

              ! Equilibrium temperature calculation
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

                    ! Linear interpolation between energies emitted for
                    ! neighboring sampled temperatures
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
                 ! Impossible to define temperature probability
                 t_cool = 1.0 ; t_abs = 0.0
                 write(*,*) "ERROR : temperature of non equilibrium grains is larger than", T_max
                 write(*,*) "cell", icell, "R=", real(r_grid(icell)), real(dust_density_o_n_grains(p_l,icell) * n_grains(l)) , &
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
                 ! status change if the grain was at equilibrium at the previous iteration
                 lchange_nRE(l,icell) = l_RE(l,icell)
                 l_RE(l,icell) = .false.
                 ! heating terms : select lower triangle
                 ! Afi with with f > i: photons with an energy between the edges of bin f minus U_i,
                 ! integral of kJnu/nu dnu (number of photons), from the field averaged over each
                 ! wavelength bin, so that the heating rates add up to the absorbed power on any
                 ! wavelength grid. The top bin receives all the photons that would heat the grain
                 ! above it (Draine & Li 2001)
                 kJnu_o_nu(:,id) = kJnu(:,id) / tab_nu(:)
                 N_above(0,id) = 0.0_dp
                 do lambda=1, n_lambda
                    N_above(lambda,id) = N_above(lambda-1,id) + kJnu_o_nu(lambda,id) * (nu_edge(lambda-1) - nu_edge(lambda))
                 enddo
                 do T2=1,n_T-1
                    k = k_edge(T2,T2)
                    N_lo = N_above(k,id) + kJnu_o_nu(k+1,id) * dnu_edge(T2,T2)
                    do T=T2+1,n_T-1
                       k = k_edge(T,T2)
                       N_hi = N_above(k,id) + kJnu_o_nu(k+1,id) * dnu_edge(T,T2)
                       A(T,T2,id) = N_lo - N_hi
                       N_lo = N_hi
                    end do
                    A(n_T,T2,id) = N_lo
                 end do

                 !! diagonal terms
                 !do T=1,n_T
                 !   A(T,T,id) = -sum(A(T,:,id)) ! sum along k direction (diag. should be zero initially)
                 !enddo

                 !! compute probabilities P(T) using iterative formulation of GD89

                 !! compute B array, stored transposed (B(T2,T): heating rate from T2 to T or above,
                 !! contiguous in T2 for the P recursion). Every element used is set for each cell
                 B(1:n_T-1,n_T,id) = A(n_T,1:n_T-1,id)
                 do T=n_T-1,1,-1
                    do T2=1,T-1
                       B(T2,T,id) = A(T,T2,id) + B(T2,T+1,id)
                    enddo
                 enddo

                 !! compute P array
                 X(1,id) = 1.e-250_dp
                 do T=2, n_T
                    if (X(T-1,id) <  1.e-300_dp) X(T-1,id) = 1.e-300_dp
                    if (X(T-1,id) >  1.e250_dp) X(1:T-1,id) = X(1:T-1,id) * 1.0e-50_dp ! allows to stabilize in case of rounding errors
                    X(T,id) =  sum(B(1:T-1,T,id)*X(1:T-1,id)) / max(A(T-1, T,id),tiny_dp)
                    if (A(T-1, T,id) < tiny_dp) then
                       write(*,*) l, T, id, "normalization error"
                       call exit(1)
                    endif

                 enddo


                 !X(1,id) = 1.e-250_dp
                 !do T=2, n_T
                 !   if (X(T-1,id) <  1.e-300_dp) X(T-1,id) = 1.e-300_dp
                 !   if (X(T-1,id) >  1.e250_dp) X(1:T-1,id) = X(1:T-1,id) * 1.0e-50_dp ! allows to stabilize in case of rounding errors
                 !   X(T,id) =  abs(sum(A(1:T-1,T-1,id)*X(1:T-1,id)) / max(A(T-1, T,id),tiny_dp))
                 !enddo

                 !! Normalisation
                 X(:,id)=X(:,id)/sum(X(:,id))

                 do T=1, n_T ! boucle car double to single array assignment bug with ifort
                    Proba_Tdust(T,l,icell) = X(T,id) ! P(T).dT, probability of being in bin of temperature T
                 enddo

              else  ! test : t_cool > t_abs : on est a a l'equilibre

                 if (.not.l_RE(l,icell)) then ! not at equilibrium at previous iteration
                    l_RE(l,icell) = .true.
                    lchange_nRE(l,icell) = .true.
                 else
                    lchange_nRE(l,icell) = .false. ! already at equilibrium
                 endif

              endif ! test : t_cool vs t_abs
           endif ! dust_density_o_n_grains > 0.
        endif ! l_dark_zone
     enddo !icell
     !$omp end do
     !$omp end parallel
  enddo !l

  deallocate(A, B, X, kJnu, lambda_Jlambda, kJnu_o_nu, N_above, k_edge, dnu_edge)

  write(*,*) "Done"

  lconverged = .true.
  ! Verification of temperature convergence
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

!  lconverged = .true. ! TMP

  return

end subroutine Temp_nRE

!********************************************************************

pure subroutine photon_edge(nu, nu_edge, k, dnu)
  ! Position of frequency nu in bins of decreasing edges nu_edge(0:n), such that the
  ! integral above nu of a density constant in each bin (N_above(k) the integral above
  ! nu_edge(k)) is N_above(k) + density(k+1) * dnu, with 0 <= k < n and dnu >= 0:
  ! k = 0, dnu = 0 above nu_edge(0); k = n-1, dnu = full last bin below nu_edge(n)

  real(kind=dp), intent(in) :: nu
  real(kind=dp), dimension(0:), intent(in) :: nu_edge
  integer, intent(out) :: k
  real(kind=dp), intent(out) :: dnu

  integer :: n, kmid, kmax

  n = size(nu_edge) - 1
  if (nu >= nu_edge(0)) then
     k = 0 ; dnu = 0.0_dp
  else if (nu <= nu_edge(n)) then
     k = n-1 ; dnu = nu_edge(n-1) - nu_edge(n)
  else
     ! nu_edge(kmax) <= nu < nu_edge(k)
     k = 0 ; kmax = n
     do while (kmax - k > 1)
        kmid = (k + kmax) / 2
        if (nu_edge(kmid) > nu) then
           k = kmid
        else
           kmax = kmid
        endif
     enddo
     dnu = nu_edge(k) - nu
  endif

end subroutine photon_edge

!********************************************************************

subroutine im_reemission_qRE(id,icell,p_icell,rand1,rand2,lambda)
! Calculation of the cell temperature
! Re-emission of a photon at the right wavelength

  use radiation_field, only : xJ_abs, J0

  implicit none

  integer, intent(in) :: id, icell, p_icell ! p_icell usused
  real, intent(in) :: rand1, rand2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, lambda0, ilambda, p_k
  real(kind=dp) :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_E_abs, J_abs

  integer, parameter :: heating_method = 3

  lambda0=lambda

  k = select_absorbing_grain(lambda0,icell, rand1, heating_method)

  ! Mean intensity: estimate of this thread x nb_proc, as for LTE grains (Temp_LTE), instead
  ! of summing the other threads at each event; the temperature also follows this thread
  J_abs=0.0
  do ilambda=1, n_lambda
     J_abs =  J_abs + C_abs_norm(k,ilambda)  * (xJ_abs(icell,ilambda,id) * nb_proc + J0(icell,ilambda))
  enddo ! ilambda
  ! WARNING: must divide by dust_density_o_n_grains as it is not taken into account in E_em_1grain
  log_E_abs=log(J_abs*L_packet_th/volume(icell))

  ! Sampled temperature just above the cell temperature
  T_int=xT_ech_1grain_nRE(k,icell,id)

  ! Eventually increment the temperature zone
  do while((log_E_em_1grain_nRE(k,T_int) < log_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! LIMITE MAX

  ! Save for next re-emission and/or final T
  xT_ech_1grain_nRE(k,icell,id)=T_int

  ! Linear interpolation between energies emitted for
  ! neighboring sampled temperatures
  T2=T_int
  Temp2=tab_Temp(T2)
  T1=T_int-1
  Temp1=tab_Temp(T1)
  frac=(log_E_abs-log_E_em_1grain_nRE(k,T1))/(log_E_em_1grain_nRE(k,T2)-log_E_em_1grain_nRE(k,T1))
  Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

  !**********************************************************************
  ! Wavelength choice for re-emission
  ! Bisection, the probability law is obtained by linear interpolation in T
  frac_T2=(Temp-Temp1)/(Temp2-Temp1)
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  do while((l2-l1) > 1)
     proba=frac_T1*kdB_dT_1grain_nRE_CDF(l,k,T1)+frac_T2*kdB_dT_1grain_nRE_CDF(l,k,T2)
     if(rand2 > proba) then
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
  ! immediate re-emission : 3 types : RE_LTE, RE_nLTE, nRE_qRE
  ! delayed re-emission : nRE

  ! need kappa_abs_tot et de kappa_abs_RE_tot : to update relative probabilities

  ! P_LTE = k_LTE / (k_LTE + k_nLTE + k_qRE) = k_LTE / k_RE
  ! P_LTE_p_nLTE = (k_LTE + k_nLTE) / (k_LTE + k_nLTE + k_qRE) = k_LTE_p_nLTE / k_RE
  ! P_abs_RE =  (k_LTE + k_nLTE + k_qRE) /  (k_LTE + k_nLTE + k_qRE + k_nRE) = k_RE / k_tot
  ! kRE is modified because k_qRE is

  ! can manage by only saving 1) k_LTE + k_nLTE = k_LTE_p_nLTE et 2) kRE iteration precedente
  ! --> donne k_RE = k_LTE_p_nLTE + k_qRE
  ! puis P_LTE = P_LTE_old * k_RE_old / k_RE_new et P_abs_RE = k_RE_new / k_RE_old * P_abs_RE_old

  real(kind=dp) :: correct, correct_m1,  kappa_abs_RE_new,  kappa_abs_RE_old, delta_kappa_abs_qRE
  integer :: l, lambda, icell, p_l
  logical :: lall_grains_eq

  write(*,*) "Setting grains at qRE and updating nRE probabilities"

  do lambda=1, n_lambda
     do icell=1,n_cells
        delta_kappa_abs_qRE = 0.0_dp
        lall_grains_eq = .true.
        do l=grain_nRE_start,grain_nRE_end
           p_l = merge(l,grain(l)%zone,lvariable_dust)
           if (lchange_nRE(l,icell)) then ! grain changed status at this iteration: added to (qRE) or
              ! removed from (back to nRE) the grains with immediate re-emission
              if (l_RE(l,icell)) then
                 delta_kappa_abs_qRE = delta_kappa_abs_qRE + C_abs_norm(l,lambda) * dust_density_o_n_grains(p_l,icell) * n_grains(l)
              else
                 delta_kappa_abs_qRE = delta_kappa_abs_qRE - C_abs_norm(l,lambda) * dust_density_o_n_grains(p_l,icell) * n_grains(l)
              endif
           endif
           if (.not.l_RE(l,icell)) lall_grains_eq = .false. ! there are grains remaining not at equilibrium
        enddo !l

        if (abs(delta_kappa_abs_qRE) > tiny_dp) then ! at least 1 grain changed status, updating the various probabilites
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
              correct = kappa_abs_RE_new / kappa_abs_RE_old ! > 1 if grains reached qRE, < 1 if they left it
              correct_m1 = 1.0_dp/correct

              ! Proba d'abs sur un grain en equilibre (pour reemission immediate)
              if (lall_grains_eq) then
                 proba_abs_RE(icell,lambda) = 1.0_dp ! 1 hardcoded to avoid rounding errors
              else
                 proba_abs_RE(icell,lambda) = proba_abs_RE(icell,lambda) * correct
              endif

              ! Among atoms at eq. proba of absorbing on a grain at LTE ou sur un grain a LTE ou nLTE
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
  ! calculates re-emission of non-equilibrium grains :
  ! - prob_E_cell
  ! - frac_E_stars a 0
  ! C. Pinte
  ! 05/02/07

  implicit none

  integer :: k, T, lambda, icell, p_k
  real(kind=dp) :: Temp, cst_wl, cst_wl_max, wl, delta_wl, fn
  real(kind=dp) :: E_emise, frac_E_abs_nRE, Delta_E, d_k
  real(kind=dp), dimension(n_cells) :: E_cell, E_cell_old

  ! emission proba as a function of lambda then for each lambda proba as a function of position
  ! differential proba as a function of proba emission precedente

  !E_abs_nRE est en nbre de paquets !!!
  E_abs_nRE = sum(E_abs_nRE_id)
  frac_E_abs_nRE = E_abs_nRE / n_photons_total
  write(*,*) "Non equilibrium grains have absorbed", real(frac_E_abs_nRE), "of emitted energy"
  write(*,*) "Re-emitting this amount of energy"

  ! Quantite energie a reemettre

  ! energie des paquets en supposant que ce soit le meme nombre de paquets que pour la premiere etape
  !E_paquet = frac_E_abs_nRE

  ! Nouvelle methode : E_paquet reste constante et on modifie le nombre de paquets
  ! + Modification E_paquet pour compenser erreur arrondi du au passage a un entier
  fn = (n_photons_total * frac_E_abs_nRE) / n_photons_loop ! "real" photon number
  n_photons_eq_th = int( max( fn, 1.0) ) ! integer photon number
  ! fn * 1.0 = n_photons_eq_th * E_paquet
  E_paquet = fn /  n_photons_eq_th

  ! For next iteration
  E_abs_nRE = 0.0 ; E_abs_nRE_id = 0.0_dp

  ! Spatial energy distribution
  cst_wl_max = log(huge_real)-1.0e-4

  spectre_emission_cumul(0) = 0.0
  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl = tab_delta_lambda(lambda)*1.e-6
     E_cell=0.0

     ! PAH emission by cell

     !$omp parallel default(none) &
     !$omp private(k,E_emise,Temp,cst_wl,T,icell,p_k,d_k) &
     !$omp shared(lambda,wl,delta_wl,E_cell,E_cell_old,tab_lambda,tab_delta_lambda,grain_nRE_start,grain_nRE_end) &
     !$omp shared(n_cells,l_RE, Tdust_1grain_nRE,n_T,C_abs_norm,dust_density_o_n_grains,n_grains,volume,tab_Temp,Proba_Tdust) &
     !$omp shared(Emissivite_nRE_old,cst_wl_max,lchange_nRE,lvariable_dust,grain)
     !$omp do
     do icell=1,n_cells
        E_emise = 0.0

        do k=grain_nRE_start,grain_nRE_end
           p_k = merge(k,grain(k)%zone,lvariable_dust)
           if (l_RE(k,icell)) then ! the grain has a temperature
              if (lchange_nRE(k,icell)) then ! the grain switches to qRE at this iteration:  it must be counted
                 d_k = dust_density_o_n_grains(p_k,icell)
                 if (d_k > tiny_dp) then
                    Temp = Tdust_1grain_nRE(k,icell)
                    cst_wl=thermal_const/(Temp*wl)
                    if (cst_wl < cst_wl_max) then
                       E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*d_k*n_grains(k)* &
                            volume(icell)/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
                    endif
                 endif
              endif ! the grain was in qRE before, it is handled in immediate re-emission
           else ! Temperature probability density
              d_k = dust_density_o_n_grains(p_k,icell)
              if (d_k > tiny_dp) then
                 do T=1,n_T
                    temp=tab_Temp(T)
                    cst_wl=thermal_const/(Temp*wl)
                    if (cst_wl < cst_wl_max) then
                       E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*d_k*n_grains(k) &
                            * volume(icell)/ ((wl**5)*(exp(cst_wl)-1.0)) * Proba_Tdust(T,k,icell) * delta_wl
                    endif !cst_wl
                 enddo !T
              endif
           endif
        enddo !k

        E_cell(icell) =   E_emise
        ! Retrieve and update old emissivity
        E_cell_old(icell) = Emissivite_nRE_old(icell,lambda)
        if (E_cell(icell) >   E_cell_old(icell)) then
           Emissivite_nRE_old(icell,lambda) = E_emise
        else ! Test in case of rounding problems
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

     ! Stars have already emitted
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
  ! initializes table Emissivite_nRE_old with emissivity at Tmin
  ! C. Pinte
  ! 05/02/07

  implicit none

  integer :: lambda, k, icell, p_k
  real(kind=dp) :: E_emise, factor, cst_wl, wl
  real(kind=dp) :: Temp, cst_wl_max, delta_wl

  cst_wl_max = log(huge_real)-1.0e-4
  temp=tab_Temp(1)

  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl = tab_delta_lambda(lambda)*1.e-6
     cst_wl=thermal_const/(Temp*wl)

     if (cst_wl < 100.0) then
        factor = 1.0_dp/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
     else
        factor = 0.0_dp
     endif

     ! PAH emission by cell
     do icell=1,n_cells
        E_emise = 0.0_dp

        do k=grain_nRE_start,grain_nRE_end
           p_k = merge(k,grain(k)%zone,lvariable_dust)
           E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*dust_density_o_n_grains(p_k,icell)*n_grains(k)* volume(icell) * factor !* Proba_Tdust = 1 pour Tmin
        enddo !k
        Emissivite_nRE_old(icell,lambda) = E_emise
     enddo !icell

  enddo ! lambda

  return

end subroutine init_emissivite_nRE

!**********************************************************************

subroutine repartition_energie(lambda)
! Calculates the distribution of energy emitted at the considered wavelength
! between the star and the various disk cells
!  - frac_E_star donne fraction emise par star
!  - prob_E_cell donne la proba d'emission cumulée des cellules
!  - E_totale donne energie totale emise (for image calibration)
! Uses a pre-tabulated temperature table
! For monochromatic code version avec scattering + em th
! C. Pinte
! 04/02/05
! 01/02/07 : grouping LTE, nLTE, and nRE cases pour faire une normalisation
! in case of different grain types

  implicit none

  integer, intent(in) :: lambda

  integer :: k, T, alloc_status, p_k
  integer, target :: icell
  integer, pointer :: p_icell
  real(kind=dp) :: Temp, wl, cst_wl, E_star, surface, E_emise, cst_wl_max
  real(kind=dp) :: delta_T
  real(kind=dp), dimension(:), allocatable :: E_cell, E_cell_corrected
  real(kind=dp), dimension(n_T) :: den_T ! Planck denominator on the temperature grid (0: no emission)

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell1
  endif

  cst_wl_max = log(huge_real)-1.0e-4

  wl = tab_lambda(lambda)*1.e-6

  ! Total star emission
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
              cst_wl=thermal_const/(Temp*wl)
              if (cst_wl < cst_wl_max) then
                 E_cell(icell) = 4.0*kappa_abs_LTE(p_icell,lambda)*kappa_factor(icell)*volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
              endif !cst_wl
           endif ! Temp==0.0
        else ! dark_zone
           E_cell(icell) = 0.0
        endif
     enddo !icell
  endif

  ! nLTE case
  if (lRE_nLTE) then
     do icell=1,n_cells
        E_emise = 0.0
        if (.not.l_dark_zone(icell)) then
           do k=grain_RE_nLTE_start,grain_RE_nLTE_end
              p_k = merge(k,grain(k)%zone,lvariable_dust)
              Temp=Tdust_1grain(k,icell)
              if (Temp > tiny_real) then
                 cst_wl=thermal_const/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    E_emise = E_emise +   4.0*C_abs_norm(k,lambda) * dust_density_o_n_grains(p_k,icell)*n_grains(k)* &
                         volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
                 endif !cst_wl
              endif ! Temp==0.0
           enddo !k
        endif
        E_cell(icell) = E_cell(icell) + E_emise
     enddo !icell
  endif

  ! nRE case
  if (lnRE) then

     ! Initialization of temperature table for images avec grains
     ! non-equilibrium
     if (lmono0) then
        delta_T=exp((1.0_dp/(real(n_T,kind=dp)))*log(T_max/T_min))
        tab_Temp(1)=T_min*sqrt(delta_T)
        do t=2,n_T
           tab_Temp(t)=delta_T*tab_Temp(t-1)
        enddo
     endif

     do T=1,n_T
        temp=tab_Temp(T)
        cst_wl=thermal_const/(Temp*wl)
        if (cst_wl < cst_wl_max) then
           den_T(T) = (wl**5)*(exp(cst_wl)-1.0)
        else
           den_T(T) = 0.0_dp
        endif
     enddo

     !$omp parallel do &
     !$omp default(none) &
     !$omp shared(n_cells,l_dark_zone,grain_nRE_start,grain_nRE_end,grain,lvariable_dust,l_RE,Tdust_1grain_nRE) &
     !$omp shared(wl,cst_wl_max,C_abs_norm,lambda,dust_density_o_n_grains,n_grains,volume,n_T,den_T,Proba_Tdust,E_cell) &
     !$omp private(icell,E_emise,k,p_k,T,Temp,cst_wl) schedule(dynamic,100)
     do icell=1,n_cells
        E_emise = 0.0
        if (.not.l_dark_zone(icell)) then

           do k=grain_nRE_start,grain_nRE_end
              p_k = merge(k,grain(k)%zone,lvariable_dust)
              if (l_RE(k,icell)) then ! the grain has a temperature
                 temp=Tdust_1grain_nRE(k,icell)
                 cst_wl=thermal_const/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*dust_density_o_n_grains(p_k,icell)*n_grains(k)* &
                         volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
                 endif !cst_wl
              else ! the grain has a probability of T
                 do T=1,n_T
                    if (den_T(T) > 0.0_dp) then
                       E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*dust_density_o_n_grains(p_k,icell)*n_grains(k)* &
                            volume(icell)/den_T(T) * Proba_Tdust(T,k,icell)
                    endif !cst_wl
                 enddo !T
              endif ! l_RE
           enddo !k
        endif
        E_cell(icell) =  E_cell(icell) + E_emise
     enddo !icell
     !$omp end parallel do

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

  ! Total energy emitted at an emitted distance equal to the Earth-star distance
  ! choose this distance to calibrate flux! pi*B(lambda)

  ! Normalisation
  ! Pas de besoin de cst : r_etoile**2 en Au2,
  ! kappa_abs en AU-1, volume en Au**3
  surface=4*pi*(pc_to_AU*distance)**2
  if (l_sym_centrale) then
     E_totale(lambda) = 2.0*pi*hp*c_light**2/surface * (E_star+E_disk(lambda)+E_ISM(lambda)) * real(N_thet)*real(N_phi)
  else
     E_totale(lambda) = 2.0*pi*hp*c_light**2/surface * (E_star+E_disk(lambda)+E_ISM(lambda)) * real(2*N_thet)*real(N_phi)
  endif

  ! Cumulative (spatial) energy distribution
  prob_E_cell(0,lambda)=0.0
  do icell=1, n_cells
     prob_E_cell(icell,lambda) = prob_E_cell(icell-1,lambda) + E_cell_corrected(icell)
  enddo

  ! Correction factor normalization
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

integer function select_absorbing_grain(lambda,icell, rand, heating_method) result(l)
  ! This routine will select randomly the absorbing/emitting grain
  ! from the CDF of kabs
  ! Because we cannot store all the CDF for all cells
  ! (n_grains x ncells x n_lambda)
  ! The CDF is recomputed on the fly here
  ! The normalization is saved via kappa_abs, so we do not need to compute
  ! all the CDF

  implicit none

  integer, intent(in) :: lambda, heating_method, icell
  real, intent(in) :: rand
  real(kind=dp) :: prob, CDF, norm
  integer :: p_icell, kstart, kend
  integer, target :: k, k1
  integer, pointer :: p_k

  if (lvariable_dust) then
     p_icell = icell
     p_k => k
  else
     k1 = 1
     p_k => k1
     p_icell = icell1
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
     if (rand < 0.5) then ! We start from first grain
        prob = rand * norm
        CDF = 0.0
        do k=kstart, kend
           CDF = CDF + C_abs(k,lambda) * dust_density_o_n_grains(p_k,icell) * n_grains(k)
           if (CDF > prob) exit
        enddo
     else ! We start from the end of the grain size distribution
        prob = (1.0-rand) * norm
        CDF = 0.0
        do k=kend, kstart, -1
           CDF = CDF + C_abs(k,lambda) * dust_density_o_n_grains(p_k,icell) * n_grains(k)
           if (CDF > prob) exit
        enddo
     endif
  else ! Same thing but with lRE to only include dust grains that are at qRE
     if (rand < 0.5) then ! We start from first grain
        prob = rand * norm
        CDF = 0.0
        do k=kstart, kend
           if (l_RE(k,icell)) CDF = CDF + C_abs(k,lambda) * dust_density_o_n_grains(p_k,icell) * n_grains(k)
           if (CDF > prob) exit
        enddo
     else ! We start from the end of the grain size distribution
        prob = (1.0-rand) * norm
        CDF = 0.0
        do k=kend, kstart, -1
           if (l_RE(k,icell)) CDF = CDF + C_abs(k,lambda) * dust_density_o_n_grains(p_k,icell) * n_grains(k)
           if (CDF > prob) exit
        enddo
     endif
  endif

  !if (k < kstart) k=kstart
  !if (k > kend)   k=kend
  l = k

  return

end function select_absorbing_grain

!**********************************************************************

subroutine select_cellule(lambda,rand, icell)
  ! Sélection de la cell qui va émettre le photon
  ! C. Pinte
  ! 04/02/05
  ! Modif 3D 10/06/05

  implicit none

  integer, intent(in) :: lambda
  real, intent(in) :: rand
  integer, intent(out) :: icell
  integer :: k, kmin, kmax

  ! Dichotomie
  kmin=0
  kmax=n_cells
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_E_cell(k,lambda) < rand) then
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
  ! Increases weight of cells near the surface
  ! par exp(-tau)
  ! Weight is applied by weight_repartion_energie
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
!--        x0=1.00001*r_lim(i-1) ! cell 1 traitee a part
!--        y0=0.0
!--        z0=0.99999*z_lim(i,j+1)
!--        u0=cos(angle)
!--        v0=0.0
!--        w0=sin(angle)
!--
!--        Stokes(:) = 0.0_dp ;
!--        call optical_length_tot(id,lambda,Stokes,icell,x0,y0,z0,u0,v0,w0,tau,lmin,lmax)
!--        if (tau < tau_min(icell)) tau_min(icell) = tau
!--
!--        x0 = 0.99999*r_lim(i)
!--        call optical_length_tot(id,lambda,Stokes,icell,x0,y0,z0,u0,v0,w0,tau,lmin,lmax)
!--        if (tau < tau_min(icell)) tau_min(icell) = tau
!--
!--     enddo
!--  enddo ! icell
!--
!--
!--  weight_proba_emission(1:n_cells) =  exp(-tau_min(:))
!--
!--  ! correct_E_emission will be normalized in energy distribution
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
  ! energy absorbed by nRE grains: only reset by emission_nRE when re-emitted, the last
  ! (converged) pass and the SED step would otherwise carry over to the next calculation
  if (lnRE) then
     E_abs_nRE = 0.0_dp ; E_abs_nRE_id = 0.0_dp
  endif
  xT_ech = 2
  Tdust = 1.0

  return

end subroutine reset_temperature

!**********************************************************************

subroutine reset_thermal_state()
  ! Restores the temperatures and the equilibrium status of the non-LTE and
  ! non-equilibrium grains as at initialisation, before the emission tables are
  ! rebuilt for a new temperature calculation on a modified structure (physics loop).
  ! The packet budget and emissivities modified by emission_nRE are restored by the
  ! caller (update_thermal_emission_tables), and kappa_abs_RE / proba_abs_RE by
  ! recompute_opacities (opacity assumes that all nRE grains are out of equilibrium).

  if (lRE_nLTE) then
     Tdust_1grain = 0.0
     xT_ech_1grain = 2
     if (allocated(Tdust_1grain_old)) Tdust_1grain_old = 0.0
  endif

  if (lnRE) then
     Tdust_1grain_nRE = 0.0
     Proba_Tdust = 0.0
     Tdust_1grain_nRE_old = 0.0
     Tpeak_old = 0
     maxP_old = 0.0
     xT_ech_1grain_nRE = 2
     l_RE = .false.
     lchange_nRE = .false.
     E_abs_nRE = 0.0_dp ; E_abs_nRE_id = 0.0_dp
  endif

  return

end subroutine reset_thermal_state

end module thermal_emission
