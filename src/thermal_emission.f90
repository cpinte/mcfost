module thermal_emission

  use parametres
  use constantes
  use em_th
  use grains
  use opacity
  use prop_star
  use disk
  !$ use omp_lib

  use utils
  use PAH
  use grid

  implicit none

  contains

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
  real :: E_star_tot, E_disk_tot, delta_wl

  spectre_emission_cumul(0) = 0.0
  ! Fonction de répartition émssion
  do lambda=1,n_lambda
     spectre_emission_cumul(lambda)=spectre_emission_cumul(lambda-1) + spectre_etoiles(lambda)/frac_E_stars(lambda)
  enddo

  ! Normalisation
   do lambda=1,n_lambda
     spectre_emission_cumul(lambda)=spectre_emission_cumul(lambda)/spectre_emission_cumul(n_lambda)
  enddo


  ! Energie des paquets
  E_star_tot = 0.0
  E_disk_tot = 0.0
  do lambda=1, n_lambda
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     E_star_tot = E_star_tot + E_stars(lambda) * delta_wl
     E_disk_tot = E_disk_tot + E_disk(lambda) * delta_wl
  enddo

  L_tot_E = (E_disk_tot + E_star_tot) / E_star_tot
  L_tot = L_etoile *  L_tot_E

  if (l_sym_centrale) then
     E_photon = L_tot  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)*(distance*pc_to_AU)**2) * real(N_thet)*real(N_phi)
  else
     E_photon = L_tot  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)*(distance*pc_to_AU)**2) * real(2*N_thet)*real(N_phi)
  endif

  n_phot_L_tot = (1.0/nbre_photons_tot) * L_tot
  n_phot_L_tot0 = n_phot_L_tot

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

subroutine init_reemission()
! Pretabule les opacites de Planck en fonction de T et de la position
! + les termes k dB/dT pour tirage de la longueur d'onde du photon emis
! avec correction de temperature
! C. Pinte 17/01/05

  implicit none

  integer :: k,lambda,t, pop, icell, p_icell
  real(kind=dp) :: integ, coeff_exp, cst_wl, cst, wl
  real ::  temp, cst_E, delta_wl
  real(kind=dp), dimension(0:n_lambda) :: integ3
  real(kind=dp), dimension(n_lambda) :: B, dB_dT

  write(*,'(a36, $)') " Initializing thermal properties ..."

  !  cst_E=2.0*hp*c_light**2/L_etoile
  ! Depuis chauffage interne, L_etoile est sorti car il aurait du etre remplace par L_tot
  ! qui depend des temp et donc de log_frac_E_em
  ! La rourine reste indépendante des autres !
  cst_E=2.0*hp*c_light**2

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
           ! Les formules prennent en compte le * lambda de l'integ log
           ! edit : non, mais delta_wl
           B(lambda) = 1.0/((wl**5)*(coeff_exp-1.0))*delta_wl
           dB_dT(lambda) = B(lambda)*cst_wl*coeff_exp/(coeff_exp-1.0) !/Temp * temp a cause de dT mais ca change rien en pratique
        else
           B(lambda)=0.0
           dB_dT(lambda)=0.0
        endif
     enddo !lambda

     ! produit par opacite (abs seule) massique
     ! Todo : this loop is OK in 2D but takes ~ 5sec for 0.5 million cells in 3D
     do icell=1,n_cells
        integ=0.0

        ! this double loop is not needed when the dust grains are the same everywhere
        do lambda=1, n_lambda
           ! kappa en Au-1    \
           ! volume en AU3     >  pas de cst pour avoir frac_E_em en SI
           ! B en SI (cst_E)  /
           ! R*-2 en AU-2    /   --> dans cst_E
           integ = integ + kappa_abs_LTE(icell,lambda)* volume(icell) * B(lambda)
        enddo !lambda

        ! Le coeff qui va bien
        integ = integ*cst_E
        if (integ > tiny_dp) then
           log_frac_E_em(T,icell)=log(integ)
        else
           log_frac_E_em(T,icell)=-1000.
        endif

        if (lRE_nLTE.and.(T==1)) then
           do lambda=1, n_lambda
              J0(icell,lambda) =  volume(icell) * B(lambda) * cst_E
           enddo
        endif
     enddo !icell


     if (low_mem_th_emission) then
        do k=grain_RE_LTE_start,grain_RE_LTE_end
           integ3(1) = 0.0
           do lambda=2, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1), ni densite
              integ3(lambda) = integ3(lambda-1) + C_abs_norm(k,lambda) * dB_dT(lambda)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_1grain_LTE_CDF(lambda,k,T) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     else  ! .not. low_mem_th_emission
        if (lvariable_dust) then ! Calcul dans toutes les cellules
           do icell=1,p_n_cells
              integ3(0) = 0.0
              do lambda=1, n_lambda
                 ! Pas besoin de cst , ni du volume (normalisation a 1)
                 integ3(lambda) = integ3(lambda-1) + kappa_abs_LTE(icell,lambda) * dB_dT(lambda)
              enddo !l

              ! Normalisation a 1
              if (integ3(n_lambda) > tiny(0.0_dp)) then
                 do lambda=1, n_lambda
                    kdB_dT_CDF(lambda,T,icell) = integ3(lambda)/integ3(n_lambda)
                 enddo !l
              endif
           enddo !icell
        else ! Pas de strat : on calcule ds une cellule non vide et on dit que ca
           ! correspond a la cellule pour prob_delta_T (car idem pour toutes les cellules)
           icell = icell_not_empty
           integ3(0) = 0.0
           do lambda=1, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1)
              integ3(lambda) = integ3(lambda-1) + kappa_abs_LTE(icell,lambda) * dB_dT(lambda)
           enddo !l

           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              p_icell = icell_ref ! When lvariable is off --> values stored in 1st cell

              do lambda=1, n_lambda
                 kdB_dT_CDF(lambda,T,p_icell) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        endif !lvariable_dust
     endif ! low_mem_th_emission

     if (lRE_nLTE) then
        ! produit par opacite (abs seule) massique d'une seule taille de grain
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ=0.0
           do lambda=1, n_lambda
              ! WARNING : il manque la densite et le volume par rapport au cas LTE !!!
              integ = integ + C_abs_norm(k,lambda) * B(lambda)
           enddo !lambda
           ! Le coeff qui va bien
           if (integ > tiny_real) then
              log_frac_E_em_1grain(k,T)=log(integ*cst_E)
           else
              log_frac_E_em_1grain(k,T)=-1000.
           endif
        enddo


        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ3(1) = 0.0
           do lambda=2, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1), ni densite
              integ3(lambda) = integ3(lambda-1) + C_abs_norm(k,lambda) * dB_dT(lambda)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_1grain_nLTE_CDF(lambda,k,T) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     endif ! lnLTE

     if (lnRE) then
        do k=grain_nRE_start,grain_nRE_end
           integ=0.0
           do lambda=1, n_lambda
              integ = integ + C_abs_norm(k,lambda) * B(lambda)
           enddo !lambda
           integ = integ * cst_E
           frac_E_em_1grain_nRE(k,T) = integ
           if (integ > 0.0_dp) then
              log_frac_E_em_1grain_nRE(k,T) = log(integ)
           else
              write(*,*) "Error in init_reemission"
              write(*,*) "Pb in opacity of non-equilibrium grains"
              write(*,*) "grain", k, "C_abs_norm seems to be incorrect"
              write(*,*) C_abs_norm(k,:)
              write(*,*) "Exiting"
              stop
           endif
        enddo

        do k=grain_nRE_start,grain_nRE_end
           integ3(1) = 0.0
           do lambda=2, n_lambda
              ! Pas besoin de cst , ni du volume (normalisation a 1), ni densite
              integ3(lambda) = integ3(lambda-1) + C_abs_norm(k,lambda) * dB_dT(lambda)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_dp)) then
              do lambda=1, n_lambda
                 kdB_dT_1grain_nRE_CDF(lambda,k,T) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     endif !lnRE

  enddo !t

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

subroutine im_reemission_LTE(id,icell,p_icell,aleat1,aleat2,lambda)
! Calcul de la temperature de la cellule et stokage energie recue + T
! Reemission d'un photon a la bonne longeur d'onde

  implicit none

  integer, intent(in) :: id, icell, p_icell
  real, intent(in) ::  aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, heating_method
  real :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_frac_E_abs, J_abs

  ! Absorption d'un photon : on ajoute son energie dans la cellule
  !xE_abs(ri,zj,phik,id) = xE_abs(ri,zj,phik,id) + E
  !E_abs=sum(xE_abs(ri,zj,phik,:))
  !log_frac_E_abs=log(E_abs*n_phot_L_tot + E0(ri,zj,phik)) ! le E0 comprend le L_tot car il est calcule a partir de frac_E_em

  if (lreemission_stats) nbre_reemission(icell,id) = nbre_reemission(icell,id) + 1.0_dp

  J_abs=sum(xKJ_abs(icell,:)) ! plante avec sunf95 sur donald + ifort sur icluster2 car negatif (-> augmentation taille minimale des cellules dans define_grid3)
  if (J_abs > 0.) then
     log_frac_E_abs=log(J_abs*n_phot_L_tot + E0(icell)) ! le E0 comprend le L_tot car il est calcule a partir de frac_E_em
  else
     log_frac_E_abs = -300
  endif

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech(icell,:))

  ! On incremente eventuellement la zone de temperature
  do while((log_frac_E_em(T_int,icell) < log_frac_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! limite max

  ! Save pour prochaine reemission et/ou T finale
  xT_ech(icell,id)=T_int

  ! Interpolation lineaire entre energies emises pour des
  ! temperatures echantillonees voisines
  T2=T_int
  Temp2=tab_Temp(T2)
  T1=T_int-1
  Temp1=tab_Temp(T1)

  frac=(log_frac_E_abs-log_frac_E_em(T1,icell))/(log_frac_E_em(T2,icell)-log_frac_E_em(T1,icell))
  Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

  !**********************************************************************
  ! Choix de la longeur d'onde de reemission
  ! Dichotomie, la loi de proba est obtenue par interpolation lineaire en T
  frac_T2=(Temp-Temp1)/(Temp2-Temp1)
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

  implicit none

  integer, intent(in) :: id,icell, p_icell
  real, intent(in) :: aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, kmin, kmax, lambda0, ilambda, heating_method
  real :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_frac_E_abs, J_abs

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
 ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans frac_E_em_1grain
  log_frac_E_abs=log(J_abs*n_phot_L_tot/volume(icell) )

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech_1grain(k,icell,:))

  ! On incremente eventuellement la zone de temperature
  do while((log_frac_E_em_1grain(k,T_int) < log_frac_E_abs).and.(T_int < n_T))
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
  frac=(log_frac_E_abs-log_frac_E_em_1grain(k,T1))/(log_frac_E_em_1grain(k,T2)-log_frac_E_em_1grain(k,T1))
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

  implicit none

  real, dimension(n_cells) :: KJ_abs
  real :: Temp, Temp1, Temp2, frac, log_frac_E_abs
  integer :: T_int, T1, T2, icell, i

  ! Calcul de la temperature de la cellule et stokage energie recue + T
  ! Utilisation temperature precedente

  ! Somme sur differents processeurs
  ! boucle pour eviter des problemes d'allocation memoire en 3D
  !KJ_abs(:) = sum(xKJ_abs(:,:),dim=2)
  KJ_abs(:) = 0.0
  do i=1,nb_proc
     KJ_abs(:) = KJ_abs(:) + xKJ_abs(1:n_cells,i)
  enddo

  KJ_abs(:)= KJ_abs(:)*n_phot_L_tot + E0(1:n_cells) ! le E0 comprend le L_tot car il est calcule a partir de frac_E_em

  !$omp parallel &
  !$omp default(none) &
  !$omp private(log_frac_E_abs,T_int,T1,T2,Temp1,Temp2,Temp,frac,icell) &
  !$omp shared(KJ_abs,xT_ech,log_frac_E_em,Temperature,tab_Temp,n_cells,T_min,n_T)
  !$omp do schedule(dynamic,10)
  do icell=1,n_cells
     if (KJ_abs(icell) < tiny_real) then
        Temperature(icell) = T_min
     else
        log_frac_E_abs=log(KJ_abs(icell))
        if (log_frac_E_abs <  log_frac_E_em(1,icell)) then
           Temperature(icell) = T_min
        else
           ! Temperature echantillonee juste sup. a la temperature de la cellule
           !          xT_int(:)=xT_ech(:,i,j)
           !          T_int=maxval(xT_int)
           T_int=maxval(xT_ech(icell,:))

           ! On incremente eventuellement la zone de temperature
           do while((log_frac_E_em(T_int,icell) < log_frac_E_abs).and.(T_int < n_T))
              T_int=T_int+1
           enddo  ! LIMITE MAX

           ! Interpolation lineaire entre energies emises pour des
           ! temperatures echantillonees voisines
           T2=T_int
           Temp2=tab_Temp(T2)
           T1=T_int-1
           Temp1=tab_Temp(T1)
           frac=(log_frac_E_abs-log_frac_E_em(T1,icell))/(log_frac_E_em(T2,icell)-log_frac_E_em(T1,icell))
           Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

           ! Save
           Temperature(icell)=Temp
        endif
     endif
  enddo !icell
  !$omp enddo
  !$omp end parallel


! if (l3D)  Temperature(:,0,1) = T_min

  if (maxval(Temperature) > T_max) then
     write(*,*) "WARNING : temperature > sublimation temperature"
     write(*,*) "WARNING : temperature = ", maxval(Temperature)
  else
     write(*,*) "Max. temperature = ", maxval(Temperature)
  endif

  return

end subroutine Temp_finale

!**********************************************************************

subroutine set_min_temperature(Tmin)

  real, intent(in) :: Tmin

  integer :: icell, i

  i=0
  do icell=1, n_cells
     if (Temperature(icell) < Tmin) then
        Temperature(icell) = Tmin
        i = i+1
     endif
  end do

  write(*,*) "Temperature was set to ", Tmin, "in", i, "cells"

  return

end subroutine set_min_temperature

!**********************************************************************

subroutine Temp_finale_nLTE()
! Calcule la temperature finale du disque dans chaque cellule
! Cas de l'"absoprtion continue"
! C. Pinte
! 24/01/05

  implicit none

  integer :: T_int, T1, T2, icell
  real(kind=dp) :: Temp, Temp1, Temp2, frac, log_frac_E_abs, J_absorbe

  integer :: k, lambda

  ! Calcul de la temperature de la cellule et stokage energie recue + T
  ! Utilisation temperature precedente

  !$omp parallel &
  !$omp default(none) &
  !$omp private(log_frac_E_abs,T_int,T1,T2,Temp1,Temp2,Temp,frac,icell) &
  !$omp shared(J_absorbe,n_phot_L_tot,xT_ech,log_frac_E_em,Temperature,tab_Temp,n_cells,n_lambda,kappa_abs_LTE) &
  !$omp shared(xJ_abs,densite_pouss,Temperature_1grain, xT_ech_1grain,log_frac_E_em_1grain) &
  !$omp shared(C_abs_norm,volume, grain_RE_nLTE_start, grain_RE_nLTE_end, n_T, T_min, J0)
  !$omp do schedule(dynamic,10)
  do icell=1,n_cells
     do k=grain_RE_nLTE_start, grain_RE_nLTE_end
        if (densite_pouss(k,icell) > tiny_dp) then
           J_absorbe=0.0
           do lambda=1, n_lambda
              J_absorbe =  J_absorbe + C_abs_norm(k,lambda)  * (sum(xJ_abs(icell,lambda,:)) + J0(icell,lambda))
           enddo ! lambda

           ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans frac_E_em_1grain
           J_absorbe = J_absorbe*n_phot_L_tot/volume(icell)
           if (J_absorbe < tiny_dp) then
              Temperature_1grain(k,icell) = T_min
           else
              log_frac_E_abs=log(J_absorbe)

              if (log_frac_E_abs <  log_frac_E_em_1grain(k,1)) then
                 Temperature_1grain(k,icell) = T_min
              else
                 ! Temperature echantillonee juste sup. a la temperature de la cellule
                 T_int=maxval(xT_ech_1grain(k,icell,:))

                 ! On incremente eventuellement la zone de temperature
                 do while((log_frac_E_em_1grain(k,T_int) < log_frac_E_abs).and.(T_int < n_T))
                    T_int=T_int+1
                 enddo  ! LIMITE MAX

                 ! Interpolation lineaire entre energies emises pour des
                 ! temperatures echantillonees voisines
                 T2=T_int
                 Temp2=tab_Temp(T2)
                 T1=T_int-1
                 Temp1=tab_Temp(T1)
                 frac=(log_frac_E_abs-log_frac_E_em_1grain(k,T1))/(log_frac_E_em_1grain(k,T2)-log_frac_E_em_1grain(k,T1))
                 Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

                 ! Save
                 Temperature_1grain(k,icell)=Temp
              endif
           endif
        else
           Temperature_1grain(k,icell)=0.0
        endif
     enddo !k
  enddo !icell
  !$omp enddo
  !$omp end parallel

  if (maxval(Temperature_1grain) > T_max) then
     write(*,*) "WARNING : temperature > sublimation temperature"
     write(*,*) "WARNING : temperature = ", maxval(Temperature_1grain)
  else
     write(*,*) "Max. temperature NLTE = ", maxval(Temperature_1grain)
  endif

  return

end subroutine Temp_finale_nLTE

!**********************************************************************

subroutine Temp_nRE(lconverged)
! C. Pinte
! 26/01/07

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
  real :: E_max
  real(kind=dp), dimension(n_cooling_time+1) :: en_lim
  real(kind=dp), dimension(n_cooling_time) :: en, delta_en, Cabs

  real :: Temp1, Temp2, frac, log_frac_E_abs

  integer :: Tpeak
  real :: maxP

  real(kind=dp) :: somme1, somme2, wl, wl_mum, delta_wl, Temp, cst_wl ! pour test

  write(*,*) "Calculating temperature probabilities of non equilibrium grains ..."
  if (lforce_PAH_equilibrium) write(*,*) "Forcing equilibrium " ;

  allocate(A(n_T,n_T,nb_proc), B(n_T,n_T,nb_proc), X(n_T,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error A'
     stop
  endif
  A=0.0 ; B=0.0 ;

  allocate(kJnu(n_lambda,nb_proc), lambda_Jlambda(n_lambda,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error A'
     stop
  endif
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
        A(T-1,T,:) = frac_E_em_1grain_nRE(l,T) / (nu_bin(T)-nu_bin(T-1)) ! OK, il manque un 1/h par rapport a Krugel p265
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
     !$omp private(frac,T1,Temp1,Temp2,T_int,log_frac_E_abs,icell) &
     !$omp shared(l,kJnu, lambda_Jlambda, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium) &
     !$omp shared(n_cells, C_abs_norm_o_dnu, xJ_abs, J0, n_phot_L_tot, volume, n_T, disk_zone,etoile) &
     !$omp shared(tab_nu, n_lambda, tab_delta_lambda, tab_lambda,en,delta_en,Cabs) &
     !$omp shared(delta_nu_bin,Proba_temperature, A,B,X,nu_bin,tab_Temp,T_min,T_max,lbenchmark_SHG,lMathis_field,Mathis_field) &
     !$omp shared(Temperature_1grain_nRE,log_frac_E_em_1grain_nRE,cst_t_cool,C_abs_norm,l_RE,r_grid) &
     !$omp shared(densite_pouss,l_dark_zone,Temperature,lchange_nRE)

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

              lambda_Jlambda(:,id) =  lambda_Jlambda(:,id)*n_phot_L_tot * (1.0/volume(icell))
              Int_k_lambda_Jlambda = Int_k_lambda_Jlambda*n_phot_L_tot * (1.0/volume(icell))
              kJnu(:,id) = kJnu(:,id)*n_phot_L_tot * (1.0/volume(icell))

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
              if (t_abs > tiny_real) then
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
              if (Int_k_lambda_Jlambda < tiny_real) then
                 Temperature_1grain_nRE(l,icell) = T_min
              else
                 log_frac_E_abs=log(Int_k_lambda_Jlambda)

                 if (log_frac_E_abs <  log_frac_E_em_1grain_nRE(l,1)) then
                    Temperature_1grain_nRE(l,icell) = T_min
                 else
                    T_int=2
                    do while((log_frac_E_em_1grain_nRE(l,T_int) < log_frac_E_abs).and.(T_int < n_T))
                       T_int=T_int+1
                    enddo  ! LIMITE MAX

                    ! Interpolation lineaire entre energies emises pour des
                    ! temperatures echantillonees voisines
                    T2=T_int
                    Temp2=tab_Temp(T2)
                    T1=T_int-1
                    Temp1=tab_Temp(T1)
                    frac=(log_frac_E_abs-log_frac_E_em_1grain_nRE(l,T1))/&
                         (log_frac_E_em_1grain_nRE(l,T2)-log_frac_E_em_1grain_nRE(l,T1))
                    Temperature_1grain_nRE(l,icell)=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))
                 endif
              endif

              if (Temperature_1grain_nRE(l,icell) > T_max) then
                 ! Impossible de definir proba de temperature
                 t_cool = 1.0 ; t_abs = 0.0
                 write(*,*) "ERROR : temperature of non equilibrium grains is larger than", T_max
                 write(*,*) "cell", icell, "R=", real(r_grid(icell)), real(densite_pouss(l,icell)), &
                      real(Temperature_1grain_nRE(l,icell))
                 write(*,*) "Exiting"
                 stop
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
                       stop
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
                    Proba_Temperature(T,l,icell) = X(T,id) ! P(T).dT, probability of being in bin of temperature T
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
        delta_T = Temperature(icell) - Temperature_old(icell)
        if (delta_T > precision * Temperature_old(icell)) lconverged = .false.
     enddo
     Temperature_old = Temperature
  endif

  if (lRE_nLTE) then
     do icell=1, n_cells
        do l=grain_RE_nLTE_start,grain_RE_nLTE_end
           delta_T = Temperature_1grain(l,icell) - Temperature_1grain_old(l,icell)
           if (delta_T > precision * Temperature_1grain_old(l,icell)) lconverged = .false.
        enddo
     enddo
     Temperature_1grain_old = Temperature_1grain
  endif

  do icell=1, n_cells
     do l=grain_nRE_start, grain_nRE_end
        if (l_RE(l,icell)) then
           delta_T = Temperature_1grain_nRE(l,icell) - Temperature_1grain_nRE_old(l,icell)
           if (delta_T > precision * Temperature_1grain_nRE_old(l,icell)) lconverged = .false.
           Temperature_1grain_nRE_old(l,icell) =  Temperature_1grain_nRE(l,icell)
        else
           Tpeak = 1
           maxP =  Proba_Temperature(1,l,icell)
           do T=2, n_T
              if (Proba_Temperature(T,l,icell) > maxP) then
                 maxP=Proba_Temperature(T,l,icell)
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
           !temp=Temperature_1grain(i,j,l) ! Tab_Temp(100)
           do lambda=1, n_lambda
              somme1 = somme1 + C_abs_norm(l,lambda)  * lambda_Jlambda(lambda,1)
              wl = tab_lambda(lambda)*1.0e-6
              delta_wl = tab_delta_lambda(lambda)*1.0e-6
              if (l_RE(l,icell)) then
                 Temp = Temperature_1grain_nRE(l,icell)
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
                            Proba_Temperature(T,l,icell)
                    endif !cst_wl
                 enddo
              endif
           enddo ! lambda

           somme2=somme2*2.0*hp*c_light**2
           !write(*,*) i,j,l,l_RE(i,j,l), Temperature_1grain_nRE(i,j,l), real(somme1), real(somme2), real(somme1/somme2)
           if (somme2 > tiny_dp) then
              Proba_Temperature(:,l,icell)  = Proba_Temperature(:,l,icell) * real(somme1/somme2)
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

  implicit none

  integer, intent(in) :: id, icell, p_icell ! p_icell usused
  real, intent(in) :: aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, lambda0, ilambda
  real :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_frac_E_abs, J_abs

  integer, parameter :: heating_method = 3

  lambda0=lambda

  k = select_absorbing_grain(lambda0,icell, aleat1, heating_method)

  ! Mean intensity
  ! Somme sur differents processeurs
  J_abs=0.0
  do ilambda=1, n_lambda
     J_abs =  J_abs + C_abs_norm(k,ilambda)  * (sum(xJ_abs(icell,ilambda,:)) + J0(icell,lambda))
  enddo ! ilambda
  ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans frac_E_em_1grain
  log_frac_E_abs=log(J_abs*n_phot_L_tot/volume(icell))

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech_1grain_nRE(k,icell,:))

  ! On incremente eventuellement la zone de temperature
  do while((log_frac_E_em_1grain_nRE(k,T_int) < log_frac_E_abs).and.(T_int < n_T))
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
  frac=(log_frac_E_abs-log_frac_E_em_1grain_nRE(k,T1))/(log_frac_E_em_1grain_nRE(k,T2)-log_frac_E_em_1grain_nRE(k,T1))
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
              write(*,*) "Oups, opacity of equilibrium grains is 0, cannot perform correction"
              write(*,*) "Something went wrong."
              write(*,*) "Having only grains out of equilibrium is not implemented yet."
              write(*,*) "Cell #", icell, " lambda #", lambda
              write(*,*) "Exiting"
              stop
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
  real :: Temp, cst_wl, cst_wl_max, wl, delta_wl, fn
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
     !$omp shared(n_cells,l_RE, Temperature_1grain_nRE,n_T,C_abs_norm,densite_pouss,volume,tab_Temp,Proba_Temperature) &
     !$omp shared(Emissivite_nRE_old,cst_wl_max,lchange_nRE)
     !$omp do
     do icell=1,n_cells
        E_emise = 0.0
        do k=grain_nRE_start,grain_nRE_end
           if (l_RE(k,icell)) then ! le grain a une temperature
              if (lchange_nRE(k,icell)) then ! la grain passe en qRE a cette iteration : il faut le compter
                 Temp = Temperature_1grain_nRE(k,icell)
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
                         ((wl**5)*(exp(cst_wl)-1.0)) * Proba_Temperature(T,k,icell) * delta_wl
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

     if (prob_E_cell(n_cells,lambda) > tiny_real) then
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
  real :: Temp, cst_wl_max, delta_wl

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
           E_emise = E_emise + 4.0*C_abs_norm(k,lambda)*densite_pouss(k,icell)* volume(icell) * facteur !* Proba_Temperature = 1 pour Tmin
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

  integer :: k, T, icell, alloc_status
  real(kind=dp) :: Temp, wl, cst_wl, E_star, surface, E_emise, cst_wl_max
  real(kind=dp) :: delta_T
  real(kind=dp), dimension(:), allocatable :: E_cell, E_cell_corrected

  cst_wl_max = log(huge_real)-1.0e-4

  wl = tab_lambda(lambda)*1.e-6

  ! Emission totale des etoiles
  E_star = E_stars(lambda)

  allocate(E_cell(n_cells), E_cell_corrected(n_cells), stat = alloc_status)
  if (alloc_status /= 0) then
     write(*,*) "Allocation error in repartition_energie"
     write(*,*) "Exiting"
     stop
  endif
  E_cell(:) = 0.0_dp
  E_cell_corrected(:) = 0.0_dp

  ! Cas LTE
  if (lRE_LTE) then
     do icell=1, n_cells
        if (.not.l_dark_zone(icell)) then
           Temp=Temperature(icell)
           if (Temp < tiny_real) then
              ! E_cell(icell) = 0.0_dp
           else
              cst_wl=cst_th/(Temp*wl)
              if (cst_wl < cst_wl_max) then
                 E_cell(icell) = 4.0*kappa_abs_LTE(icell,lambda)*volume(icell)/((wl**5)*(exp(cst_wl)-1.0))
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
              Temp=Temperature_1grain(k,icell)
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
                 temp=Temperature_1grain_nRE(k,icell)
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
                            volume(icell)/((wl**5)*(exp(cst_wl)-1.0)) * Proba_Temperature(T,k,icell)
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

  if (prob_E_cell(n_cells,lambda) > tiny_real) then
     prob_E_cell(:,lambda)=prob_E_cell(:,lambda)/prob_E_cell(n_cells,lambda)
  else
     prob_E_cell(:,lambda)=0.0
  endif

  deallocate(E_cell,E_cell_corrected)

  return

end subroutine repartition_energie

!**********************************************************************

subroutine chauffage_interne()
  ! Calcule le temperature initiale du disque avant le step 1
  ! C. Pinte
  ! 27/02/06

  implicit none

  if (lRE_LTE) then
     ! Energie venant de l'equilibre avec nuage à T_min
     E0(1:n_cells) = exp(log_frac_E_em(1,1:n_cells))

!!$  ! Energie venant du chauffage visqueux
!!$  do i=1,n_rad
!!$     do j=1,nz
!!$        E0(i,j) = E0(i,j) !+ masse(i) * alpha *
!!$     enddo
!!$  enddo
!!$
!!$  ! Calcul temperature
!!$  if (lLTE) then
!!$     do i=1, nb_proc
!!$        xKJ_abs(:,:,i) = E0(:,:)/nb_proc
!!$     enddo
!!$     call Temp_finale
!!$  else
!!$     write(*,*) "Coding in progress ..."
!!$     stop
!!$  endif

     ! Energie emise aux differente longueurs d'onde : repartition_energie
     call Temp_finale()
  endif

  return

end subroutine chauffage_interne

!**********************************************************************

integer function select_absorbing_grain(lambda,icell, aleat, heating_method) result(k)
  ! This routine will select randomly the scattering grain
  ! from the CDF of kabs
  ! Because we cannot store all the CDF for all cells
  ! (n_grains x ncells x n_lambda)
  ! The CDF is recomputed on the fly here
  ! The normalization is saved via kappa_abs, so we do not need to compute
  ! all the CDF

  implicit none

  integer, intent(in) :: lambda, icell, heating_method
  real, intent(in) :: aleat
  real :: prob, CDF, norm
  integer :: kstart, kend

  ! We scale the random number so that it is between 0 and kappa_sca (= last value of CDF)
  if (heating_method == 1) then
     norm =  kappa_abs_LTE(icell,lambda) / ( AU_to_cm * mum_to_cm**2 )
     kstart = grain_RE_LTE_start ; kend = grain_RE_LTE_end
  else if (heating_method == 2) then
     norm =  kappa_abs_nLTE(icell,lambda) / ( AU_to_cm * mum_to_cm**2 )
     kstart = grain_RE_nLTE_start ; kend = grain_RE_nLTE_end
  else
     ! todo : maybe update with an extra variable kappa_abs_qRE
     if (lRE_nLTE) then
        norm =  (kappa_abs_RE(icell, lambda) -  kappa_abs_LTE(icell,lambda) - kappa_abs_nLTE(icell,lambda)) &
             / ( AU_to_cm * mum_to_cm**2 )
     else
        norm =  (kappa_abs_RE(icell, lambda) -  kappa_abs_LTE(icell,lambda)) &
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


  return

end function select_absorbing_grain

!**********************************************************************

subroutine reset_radiation_field()

  if (lRE_LTE) then
     xKJ_abs(:,:) = 0.0_dp
     E0 = 0.0_dp
  endif
  if (lRE_nLTE .or. lnRE) xJ_abs(:,:,:) = 0.0_dp
  xT_ech = 2
  Temperature = 1.0


!  E_stars = 0.0 ; E_disk = 0.0 ; E_ISM = 0.0
!
!  frac_E_stars = 0.0 ; frac_E_disk = 0.0 ; E_totale = 0.0
!
!  spectre_etoiles_cumul = 0.0
!  spectre_etoiles = 0.0
!  spectre_emission_cumul = 0.0
!
!  prob_E_cell = 0.0
!
!  log_frac_E_em = 0.0
!
!
!  kdB_dT_CDF = 0

  return

end subroutine reset_radiation_field

!**********************************************************************

end module thermal_emission
