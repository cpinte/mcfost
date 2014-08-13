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

  integer :: i,j,pk,k,lambda,t, pop
  real(kind=db) :: integ, integ2, coeff_exp, cst_wl, cst, wl
  real ::  temp, cst_E, delta_wl
  real(kind=db), dimension(0:n_lambda) :: integ3
  real(kind=db), dimension(n_lambda) :: B, dB_dT

  ! Pour spline
  real, dimension(1:n_T) :: xa,ya,y2
  real ::  yp1, ypn
  real, dimension(n_rad,nz) :: dfrac_E_em_1, dfrac_E_em_n

  lxJ_abs = loutput_J .or. loutput_UV_field .or. lRE_nLTE .or. lnRE !.or. lProDiMo

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
           !        write(*,*) lambda, wl, B(lambda)
           dB_dT(lambda) = B(lambda)*cst_wl*coeff_exp/(coeff_exp-1.0) !/Temp * temp a cause de dT mais ca change rien en pratique
        else
           B(lambda)=0.0
           dB_dT(lambda)=0.0
        endif
     enddo !lambda

     ! produit par opacite (abs seule) massique
     do i=1,n_rad
        bz : do j=j_start, nz
           if (j==0) cycle bz
           do pk=1, n_az
              integ=0.0

              do lambda=1, n_lambda
                 ! kappa en Au-1    \
                 ! volume en AU3     >  pas de cst pour avoir frac_E_em en SI
                 ! B en SI (cst_E)  /
                 ! R*-2 en AU-2    /   --> dans cst_E
                 integ = integ + kappa_abs_eg(lambda,i,j,pk)* volume(i) * B(lambda)
                 !              write(*,*) i,j,lambda,integ, kappa_abs(lambda,i,j),volume(i), B(lambda)
                 !              read(*,*)

              enddo !lambda


              ! Le coeff qui va bien
              integ = integ*cst_E
              if (integ > tiny_db) then
                 log_frac_E_em(i,j,pk,T)=log(integ)
              else
                 log_frac_E_em(i,j,pk,T)=-1000.
              endif

              if (T==1) then
                 do lambda=1, n_lambda
                    J0(lambda,i,j,pk) =  volume(i) * B(lambda) * cst_E
                 enddo
              endif

              !              write(*,*) log_frac_E_em(1,1,1,1), log_frac_E_em(1,-1,1,1)

           !!! Pour spline, calcul de la derivee / T de frac_E_em
!!$           if ((t==1).and.(j<=nz)) then
!!$              integ=0.0
!!$              do lambda=1, n_lambda
!!$                 integ = integ + kappa_abs(lambda,i,j)*volume(i)*dB_dT(lambda)
!!$              enddo !lambda
!!$              dfrac_E_em_1(i,j)=integ*cst_E
!!$           endif
!!$
!!$           if ((t==n_T).and.(j<=nz)) then
!!$              integ=0.0
!!$              do lambda=1, n_lambda
!!$                 integ = integ + kappa_abs(lambda,i,j)*volume(i)*dB_dT(lambda)
!!$              enddo !lambda
!!$              dfrac_E_em_n(i,j)=integ*cst_E
!!$           endif
           !!! Fin spline
           enddo !pk
        enddo bz !j
     enddo !i

     if (lstrat) then ! Calcul dans toutes les cellules
        do i=1,p_n_rad
           bz2 : do j=pj_start, p_nz
              if (j==0) cycle bz2
              do pk=1, p_n_az
                 integ3(0) = 0.0
                 do lambda=1, n_lambda
                    ! Pas besoin de cst , ni du volume (normalisation a 1)
                    integ3(lambda) = integ3(lambda-1) + kappa_abs_eg(lambda,i,j,pk) * dB_dT(lambda)
                 enddo !l

                 ! Normalisation a 1
                 if (integ3(n_lambda) > tiny(0.0_db)) then
                    do lambda=1, n_lambda
                       prob_delta_T(i,j,pk,T,lambda) = integ3(lambda)/integ3(n_lambda)
                    enddo !l
                 endif
              enddo !pk
           enddo bz2 !j
        enddo !i
     else ! Pas de strat : on calcule ds une cellule non vide et on dit que ca
        ! correspond a la cellule pour prob_delta_T (car idem pour toutes les cellules)
        i=ri_not_empty
        j=zj_not_empty
        pk=phik_not_empty
        integ3(0) = 0.0
        do lambda=1, n_lambda
           ! Pas besoin de cst , ni du volume (normalisation a 1)
           integ3(lambda) = integ3(lambda-1) + kappa_abs_eg(lambda,i,j,pk) * dB_dT(lambda)
        enddo !l

        ! Normalisation a 1
        if (integ3(n_lambda) > tiny(0.0_db)) then
           do lambda=1, n_lambda
              !      write(*,*) i,j,T,lambda, integ3(lambda), integ3(n_lambda)
              prob_delta_T(1,1,1,T,lambda) = integ3(lambda)/integ3(n_lambda)
           enddo !l
        endif

     endif !lstrat

     if (lRE_nLTE) then
        ! produit par opacite (abs seule) massique d'une seule taille de grain
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           integ=0.0
           do lambda=1, n_lambda
              ! WARNING : il manque la densite et le volume par rapport au cas LTE !!!
              integ = integ + q_abs(lambda,k) * B(lambda)
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
              integ3(lambda) = integ3(lambda-1) + q_abs(lambda,k) * dB_dT(lambda)
           enddo !l
           ! Normalisation a 1
           if (integ3(n_lambda) > tiny(0.0_db)) then
              do lambda=1, n_lambda
                 prob_delta_T_1grain(k,T,lambda) = integ3(lambda)/integ3(n_lambda)
              enddo !l
           endif
        enddo !k
     endif ! lnLTE

     if (lnRE) then
        do k=grain_nRE_start,grain_nRE_end
           integ=0.0
           do lambda=1, n_lambda
              integ = integ + q_abs(lambda,k) * B(lambda)
           enddo !lambda
           integ = integ * cst_E
           frac_E_em_1grain_nRE(k,T) = integ
           if (integ > 0.0_db) then
              log_frac_E_em_1grain_nRE(k,T) = log(integ)
           else
              write(*,*) "Error in init_reemission"
              write(*,*) "Pb in opacity of non-equilibrium grains"
              write(*,*) "Exiting"
              stop
           endif
           enddo
     endif !lnRE

  enddo !t

  do pop=1, n_pop
     if (lread_Misselt.and.dust_pop(pop)%is_opacity_file) call read_file_specific_heat(pop)
  enddo

  return

end subroutine init_reemission

!********************************************************************

subroutine reemission(id,ri,zj,phik,pri,pzj,pphik,E,aleat,lambda)
! Calcul de la temperature de la cellule et stokage energie recue + T
! Reemission d'un photon a la bonne longeur d'onde

  implicit none

  integer, intent(in) :: id,ri,zj,phik,pri,pzj, pphik
  real(kind=db), intent(in) :: E
  real, intent(in) ::  aleat
  integer, intent(out) :: lambda

  integer :: l, l1, l2, T_int, T1, T2
  real :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_frac_E_abs, J_abs, E_abs

  ! Absorption d'un photon : on ajoute son energie dans la cellule
  !xE_abs(ri,zj,phik,id) = xE_abs(ri,zj,phik,id) + E
  !E_abs=sum(xE_abs(ri,zj,phik,:))
  !log_frac_E_abs=log(E_abs*n_phot_L_tot + E0(ri,zj,phik)) ! le E0 comprend le L_tot car il est calcule a partir de frac_E_em

  nbre_reemission(ri,zj,phik,id) = nbre_reemission(ri,zj,phik,id) + 1.0_db

  J_abs=sum(xKJ_abs(ri,zj,phik,:)) ! plante avec sunf95 sur donald + ifort sur icluster2 car negatif (-> augmentation taille minimale des cellules dans define_grid3)
  if (J_abs > 0.) then
     log_frac_E_abs=log(J_abs*n_phot_L_tot + E0(ri,zj,phik)) ! le E0 comprend le L_tot car il est calcule a partir de frac_E_em
  else
     log_frac_E_abs = -300
  endif

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech(:,ri,zj,phik))

  ! On incremente eventuellement la zone de temperature
  do while((log_frac_E_em(ri,zj,phik,T_int) < log_frac_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! limite max

  ! Save pour prochaine reemission et/ou T finale
  xT_ech(id,ri,zj,phik)=T_int

  ! Interpolation lineaire entre energies emises pour des
  ! temperatures echantillonees voisines
  T2=T_int
  Temp2=tab_Temp(T2)
  T1=T_int-1
  Temp1=tab_Temp(T1)

  frac=(log_frac_E_abs-log_frac_E_em(ri,zj,phik,T1))/(log_frac_E_em(ri,zj,phik,T2)-log_frac_E_em(ri,zj,phik,T1))
  Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

  !**********************************************************************
  ! Choix de la longeur d'onde de reemission
  ! Dichotomie, la loi de proba est obtenue par interpolation lineaire en T
!  write(*,*) Temp2, Temp1
!   write(*,*) "t1", Temp2, Temp1
  frac_T2=(Temp-Temp1)/(Temp2-Temp1)
  frac_T1=1.0-frac_T2

  l1=0
  l2=n_lambda
  l=(l1+l2)/2

  do while((l2-l1) > 1)
     proba=frac_T1*prob_delta_T(pri,pzj,pphik,T1,l)+frac_T2*prob_delta_T(pri,pzj,pphik,T2,l)
     if(aleat > proba) then
        l1=l
     else
        l2=l
     endif
     l=(l1+l2)/2
  enddo
  lambda=l+1

  return

end subroutine reemission

!********************************************************************

subroutine reemission_NLTE(id,ri,zj,pri,pzj,E,aleat1,aleat2,lambda)
! Calcul de la temperature de la cellule
! Reemission d'un photon a la bonne longeur d'onde

  implicit none

  integer, intent(in) :: id,ri,zj,pri,pzj
  real(kind=db), intent(in) :: E
  real, intent(in) :: aleat1, aleat2
  integer, intent(inout) :: lambda

  integer :: l, l1, l2, T_int, T1, T2, k, kmin, kmax, lambda0, ilambda
  real :: Temp, Temp1, Temp2, frac_T1, frac_T2, proba, frac, log_frac_E_abs, J_abs, E_abs

  lambda0=lambda

  ! Selection du grain qui absorbe le photon
  kmin=0
  kmax=n_grains_RE_nLTE
  k=(kmin+kmax)/2

  do while((kmax-kmin) > 1)
     if (prob_kappa_abs_1grain(lambda0,pri,pzj,k) < aleat1) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
  enddo   ! while
  k=kmax


  ! Absorption d'un photon : on ajoute son energie dans la cellule
!  xE_abs(ri,zj,id) = xE_abs(ri,zj,id) + E
!  xE_abs_1grain(ri,zj,k,id) = xE_abs_1grain(ri,zj,k,id) + E
!  E_abs=sum(xE_abs_1grain(ri,zj,k,:))
!  frac_E_abs=E_abs/(nbre_photons_tot*densite_pouss(ri,zj,k)*volume(ri))

  ! Mean intensity
  ! Somme sur differents processeurs
  J_abs=0.0
  do ilambda=1, n_lambda
     J_abs =  J_abs + q_abs(ilambda,k)  * (sum(xJ_abs(ilambda,ri,zj,:)) + J0(ilambda,ri,zj,1))
  enddo ! ilambda
 ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans frac_E_em_1grain
  log_frac_E_abs=log(J_abs*n_phot_L_tot/volume(ri) )

  ! Temperature echantillonee juste sup. a la temperature de la cellule
  T_int=maxval(xT_ech_1grain(:,ri,zj,k))

  ! On incremente eventuellement la zone de temperature
  do while((log_frac_E_em_1grain(k,T_int) < log_frac_E_abs).and.(T_int < n_T))
     T_int=T_int+1
  enddo  ! LIMITE MAX

  ! Save pour prochaine reemission et/ou T finale
  xT_ech_1grain(id,ri,zj,k)=T_int


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
     proba=frac_T1*prob_delta_T_1grain(k,T1,l)+frac_T2*prob_delta_T_1grain(k,T2,l)
     if(aleat2 > proba) then
        l1=l
     else
        l2=l
     endif
     l=(l1+l2)/2
  enddo
  lambda=l+1

  return

end subroutine reemission_NLTE

!********************************************************************

subroutine Temp_finale()
! Calcule la temperature finale du disque dans chaque cellule
! Cas de l'"absoprtion continue"
! C. Pinte
! 24/01/05

  implicit none

  integer :: l, l1, l2, T_int, T1, T2
  real :: Temp, Temp1, Temp2, proba, frac, log_frac_E_abs

  real(kind=db), dimension(:,:,:), allocatable :: J_abs

  integer :: i,j, pk

  if (l3D) then
     allocate(J_abs(n_rad,-nz:nz,n_az))
  else
     allocate(J_abs(n_rad,nz,1))
  endif

  ! Calcul de la temperature de la cellule et stokage energie recue + T
  ! Utilisation temperature precedente

  ! Somme sur differents processeurs
  J_abs=sum(xKJ_abs,dim=4)
  J_abs(:,:,:)= J_abs(:,:,:)*n_phot_L_tot + E0(:,:,:) ! le E0 comprend le L_tot car il est calcule a partir de frac_E_em

  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,pk,log_frac_E_abs,T_int,T1,T2,Temp1,Temp2,Temp,frac) &
  !$omp shared(J_abs,xT_ech,log_frac_E_em,Temperature,tab_Temp,n_rad,nz,lstrat,T_min,n_T,j_start,n_az)
  !$omp do schedule(dynamic,10)
  do i=1,n_rad
     bz : do j=j_start,nz
        if (j==0) cycle bz
        do pk=1, n_az
           if (J_abs(i,j,pk) < tiny_db) then
              Temperature(i,j,pk) = T_min
           else
              log_frac_E_abs=log(J_abs(i,j,pk))
              if (log_frac_E_abs <  log_frac_E_em(i,j,pk,1)) then
                 Temperature(i,j,pk) = T_min
              else
                 ! Temperature echantillonee juste sup. a la temperature de la cellule
                 !          xT_int(:)=xT_ech(:,i,j)
                 !          T_int=maxval(xT_int)
                 T_int=maxval(xT_ech(:,i,j,pk))

                 ! On incremente eventuellement la zone de temperature
                 do while((log_frac_E_em(i,j,pk,T_int) < log_frac_E_abs).and.(T_int < n_T))
                    T_int=T_int+1
                 enddo  ! LIMITE MAX

                 ! Interpolation lineaire entre energies emises pour des
                 ! temperatures echantillonees voisines
                 T2=T_int
                 Temp2=tab_Temp(T2)
                 T1=T_int-1
                 Temp1=tab_Temp(T1)
                 frac=(log_frac_E_abs-log_frac_E_em(i,j,pk,T1))/(log_frac_E_em(i,j,pk,T2)-log_frac_E_em(i,j,pk,T1))
                 Temp=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))

                 ! Save
                 Temperature(i,j,pk)=Temp
              endif
           endif
        enddo !pk
     enddo bz !j
  enddo !i
  !$omp enddo
  !$omp end parallel


! if (l3D)  Temperature(:,0,1) = T_min

  if (maxval(Temperature) > T_max) then
     write(*,*) "WARNING : temperature > sublimation temperature"
     write(*,*) "WARNING : temperature = ", maxval(Temperature)
  else
     write(*,*) "Max. temperature = ", maxval(Temperature)
  endif

  deallocate(J_abs)

  return

end subroutine Temp_finale

!**********************************************************************

subroutine Temp_finale_nLTE()
! Calcule la temperature finale du disque dans chaque cellule
! Cas de l'"absoprtion continue"
! C. Pinte
! 24/01/05

  implicit none

  integer :: l, l1, l2, T_int, T1, T2
  real(kind=db) :: Temp, Temp1, Temp2, proba, frac, log_frac_E_abs, J_absorbe

  real :: rcyl, z

  integer :: i,j, k, lambda

  ! Calcul de la temperature de la cellule et stokage energie recue + T
  ! Utilisation temperature precedente

!!$ !Verifications
!!$  do i=1,n_rad
!!$     do j=1,nz
!!$        J_absorbe=0.0
!!$        do k=1, n_grains_tot
!!$           do lambda=1, n_lambda
!!$              J_absorbe =  J_absorbe + kappa_abs_1grain(lambda,i,j,k)  * sum(xJ_abs(lambda,i,j,:))
!!$           enddo ! lambda
!!$        enddo
!!$        write(*,*) i,j, J_absorbe, sum(xKJ_abs(i,j,:))
!!$     enddo
!!$  enddo
!!$  do i=1,n_rad
!!$     do j=1,nz
!!$        J_absorbe=0.0
!!$        do k=1, n_grains_tot
!!$           J_absorbe= J_absorbe + frac_E_em_1grain(k,5)*densite_pouss(i,j,k)*volume(i)
!!$        enddo
!!$        write(*,*) i,j, J_absorbe, frac_E_em(i,j,5)
!!$
!!$     enddo
!!$  enddo
!!$  stop

  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,log_frac_E_abs,T_int,T1,T2,Temp1,Temp2,Temp,frac) &
  !$omp shared(J_absorbe,n_phot_L_tot,xT_ech,log_frac_E_em,Temperature,tab_Temp,n_rad,nz,lstrat,n_lambda,kappa_abs_eg) &
  !$omp shared(prob_kappa_abs_1grain,xJ_abs,densite_pouss,Temperature_1grain, xT_ech_1grain,log_frac_E_em_1grain) &
  !$omp shared(q_abs,volume, grain_RE_nLTE_start, grain_RE_nLTE_end, n_T, T_min, J0)
  !$omp do schedule(dynamic,10)
  do i=1,n_rad
     do j=1,nz
        do k=grain_RE_nLTE_start, grain_RE_nLTE_end
           if (densite_pouss(i,j,1,k) > tiny_db) then
              J_absorbe=0.0
              do lambda=1, n_lambda
                 !J_absorbe =  J_absorbe + kappa_abs_1grain(lambda,i,j,k)  * sum(xJ_abs(lambda,i,j,:))
                 J_absorbe =  J_absorbe + q_abs(lambda,k)  * (sum(xJ_abs(lambda,i,j,:)) + J0(lambda,i,j,1))
              enddo ! lambda

              ! WARNING : il faut diviser par densite_pouss car il n'est pas pris en compte dans frac_E_em_1grain
              !J_absorbe = J_absorbe*n_phot_L_tot/(densite_pouss(i,j,1,k)*volume(i))
              J_absorbe = J_absorbe*n_phot_L_tot/volume(i)
              if (J_absorbe < tiny_db) then
                 Temperature_1grain(i,j,k) = T_min
              else
                 log_frac_E_abs=log(J_absorbe)

                 if (log_frac_E_abs <  log_frac_E_em_1grain(k,1)) then
                    Temperature_1grain(i,j,k) = T_min
                 else

                    ! Temperature echantillonee juste sup. a la temperature de la cellule
                    !          xT_int(:)=xT_ech(:,i,j)
                    !          T_int=maxval(xT_int)
                    T_int=maxval(xT_ech_1grain(:,i,j,k))

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
                    Temperature_1grain(i,j,k)=Temp!J_abs(i,j)
                 endif
              endif
           else
              Temperature_1grain(i,j,k)=0.0
           endif
        enddo!k
     enddo !j
  enddo !i
  !$omp enddo
  !$omp end parallel

  if (maxval(Temperature_1grain) > T_max) then
     write(*,*) "WARNING : temperature > sublimation temperature"
     write(*,*) "WARNING : temperature = ", maxval(Temperature_1grain)
  else
     write(*,*) "Max. temperature = ", maxval(Temperature_1grain)
  endif

  return

end subroutine Temp_finale_nLTE

!**********************************************************************

subroutine Temp_nRE(lconverged)
! C. Pinte
! 26/01/07

  implicit none

  logical, intent(out) :: lconverged

  real, parameter :: precision = 1.e-2

  real(kind=db), dimension(:,:,:), allocatable :: A, B
  real(kind=db), dimension(:,:), allocatable :: X
  real(kind=db), dimension(n_T) :: nu_bin, delta_nu_bin
  real(kind=db), dimension(0:n_T) :: T_lim, U_lim

  real(kind=db), dimension(n_lambda) :: log_tab_nu, q_abs_o_dnu, tab_nu

  real(kind=db), dimension(:,:), allocatable :: KJ_absorbe_nRE, log_KJ_absorbe, Jabs

  real(kind=db) :: delta_T, KJ_abs_interp, t_cool, t_abs, mean_abs_E, kTu, mean_abs_nu, cst_t_cool, KJ_absorbe

  integer :: l, i, j, k, T, T1, T2, lambda, alloc_status, id, T_int

  integer, parameter :: n_cooling_time = 100
  real :: E_max
  real(kind=db), dimension(n_cooling_time+1) :: en_lim
  real(kind=db), dimension(n_cooling_time) :: en, delta_en, Cabs

  real :: Temp1, Temp2, frac, log_frac_E_abs

  integer :: Tpeak
  real :: maxP

  real(kind=db) :: somme1, somme2, wl, delta_wl, Temp, cst_wl ! pour test

  write(*,*) "Calculating temperature probabilities of non equilibrium grains ..."
  if (lforce_PAH_equilibrium) write(*,*) "Forcing equilibrium " ;

  allocate(A(n_T,n_T,nb_proc), B(n_T,n_T,nb_proc), X(n_T,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error A'
     stop
  endif
  A=0.0 ; B=0.0 ;

  allocate(KJ_absorbe_nRE(n_lambda,nb_proc), log_KJ_absorbe(n_lambda,nb_proc), Jabs(n_lambda,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error A'
     stop
  endif
  KJ_absorbe_nRE = 0.0 ; log_KJ_absorbe = 0.0 ; Jabs=0.0

  delta_T=exp((1.0_db/(real(n_T,kind=db)))*log(T_max/T_min))
  T_lim(0) = T_min
  T_lim(1:n_T) = tab_Temp(:)*sqrt(delta_T)

  tab_nu = c_light/(tab_lambda*1.0e-6)

  ! l est l'indice de taille de grains
  do l=grain_nRE_start, grain_nRE_end
     ! bin centers [Hz]
     nu_bin(:) = specific_heat(real(tab_Temp,kind=db),l) * tab_Temp(:)/hp ! evaluate specific heat & enthalpy at bin centers

     ! bin widths [Hz]
     U_lim(:) = specific_heat(T_lim(:),l)*T_lim(:) ! enthalpy of bin boundaries [erg]
     delta_nu_bin(:) = (U_lim(1:n_T)-U_lim(0:n_T-1))/hp  ! bin width [Hz]

     ! compute transition matrix  A   [Hz]

     !! cooling terms :  select one-off upper diagonal, on peut les calculer a l'avance
     ! Afi with final < inital
     A(:,:,:)=0.
     do T=2,n_T
        A(T-1,T,:) = frac_E_em_1grain_nRE(l,T) / (nu_bin(T)-nu_bin(T-1)) ! OK, il manque un 1/h par rapport a Krugel p265
     enddo

     ! Pour KJ_absorbe
     do lambda=1, n_lambda
        ! q_abs / dnu  : dnu = c wl^2 / dwl
        q_abs_o_dnu(lambda) = q_abs(lambda,l)  / tab_delta_lambda(lambda) * (tab_lambda(lambda))**2/c_light * 1.0e-6
     enddo

     ! Pour cooling time
     E_max = real(maxval(U_lim))
     en_lim = spanl(1.0e-6*E_max,E_max, n_cooling_time+1)
     do T=1,n_cooling_time
        en(T) = 0.5 * (en_lim(T+1) + en_lim(T))
        delta_en(T) = en_lim(T+1) - en_lim(T)
        Cabs(T) = interp(real(q_abs(:,l),kind=db),real(tab_lambda,kind=db),en(T)/hp)
     enddo

     cst_t_cool =(real(hp,kind=db)**3 * real(c_light,kind=db)**2) / (8*pi)


     ! Boucle sur les cellules

     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,KJ_absorbe, lambda, T, T2) &
     !$omp private(KJ_abs_interp,id,t_cool,t_abs,mean_abs_E,mean_abs_nu,kTu) &
     !$omp private(frac,T1,Temp1,Temp2,T_int,k,log_frac_E_abs) &
     !$omp shared(l,KJ_absorbe_nRE, Jabs, log_KJ_absorbe, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium) &
     !$omp shared(n_rad, nz, q_abs_o_dnu, xJ_abs, J0, n_phot_L_tot, volume, n_T) &
     !$omp shared(log_tab_nu, tab_nu, n_lambda, tab_delta_lambda, tab_lambda,en,delta_en,Cabs) &
     !$omp shared(delta_nu_bin,Proba_temperature, A,B,X,nu_bin,tab_Temp,T_min,T_max) &
     !$omp shared(Temperature_1grain_nRE,log_frac_E_em_1grain_nRE,cst_t_cool,q_abs,l_RE,r_grid,densite_pouss)

     id =1 ! pour code sequentiel
     ! ganulation faible car le temps calcul depend fortement des cellules
     !$omp do schedule(dynamic,1)
     do i=n_rad, 1, -1
        !$ id = omp_get_thread_num() + 1
        do j=1,nz

           if (densite_pouss(i,j,1,l) > tiny_db) then
              ! Champ de radiation
              KJ_absorbe=0.0
              do lambda=1, n_lambda
                 ! conversion lambda et delta_lambda de micron -> m
                 Jabs(lambda,id) = (sum(xJ_abs(lambda,i,j,:)) + J0(lambda,i,j,1))
                 KJ_absorbe = KJ_absorbe + q_abs(lambda,l) * Jabs(lambda,id)
                 KJ_absorbe_nRE(lambda,id) =   q_abs_o_dnu(lambda)  * Jabs(lambda,id)
              enddo ! lambda

              Jabs(:,id) =  Jabs(:,id)*n_phot_L_tot * (1.0/volume(i))
              KJ_absorbe = KJ_absorbe*n_phot_L_tot * (1.0/volume(i))
              KJ_absorbe_nRE(:,id) = KJ_absorbe_nRE(:,id)*n_phot_L_tot * (1.0/volume(i))

              ! decide whether we really need to use this model, instead of calculating the equilibrium temperature

              ! time interval for photon absorption
              t_abs = sum(q_abs(:,l)*Jabs(:,id)*tab_lambda(:)) * 1.0e-6/c_light
              t_abs = hp/(4.*pi*t_abs)

              ! mean absorbed photon energy
              mean_abs_E = t_abs * 4*pi * KJ_absorbe ! D01 eq 46
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
                       t_cool = t_cool + en(T)**3 * Cabs(T)/(exp(en(T)/kTu)-1.0_db) * delta_en(T)
                    endif
                 enddo integ
                 if (t_cool > tiny_db) then
                    t_cool = cst_t_cool*mean_abs_E/t_cool * 1.0e-4
                 else
                    t_cool = huge_real
                 endif
              endif

              ! Calcul Temperature equilibre
              l_RE(i,j,l) = .true.
              if (KJ_absorbe < tiny_real) then
                 Temperature_1grain_nRE(i,j,l) = T_min
              else
                 log_frac_E_abs=log(KJ_absorbe)

                 if (log_frac_E_abs <  log_frac_E_em_1grain_nRE(l,1)) then
                    Temperature_1grain_nRE(i,j,l) = T_min
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
                    Temperature_1grain_nRE(i,j,l)=exp(log(Temp2)*frac+log(Temp1)*(1.0-frac))
                 endif
              endif


              if (Temperature_1grain_nRE(i,j,l) > T_max) then
                 ! Impossible de definir proba de temperature
                 t_cool = 1.0 ; t_abs = 0.0
                 write(*,*) "ERROR : temperature of non equilibrium grains is larger than", T_max
                 write(*,*) "cell", i, "R=", real(r_grid(i,1)), real(densite_pouss(i,j,1,l)), real(Temperature_1grain_nRE(i,j,l))
                 write(*,*) "Exiting"
                 stop
              endif

              ! Forcing equilibrium or no equilibrium
              if (lforce_PAH_equilibrium) then
                 t_cool=1.0 ; t_abs = 0.0
              else if (lforce_PAH_out_equilibrium) then
                 t_cool = 0.0 ; t_abs = 1.0
              endif

              if (t_cool < t_abs) then !  calcul proba temperature
                 l_RE(i,j,l) = .false.
                 log_KJ_absorbe(:,id) = log(KJ_absorbe_nRE(:,id)+tiny_db)

                 ! heating terms : select lower triangle
                 ! Afi with with f > i
                 do T=1,n_T
                    do T2=1,T-1
                       ! Tres largement plus stable en interpolation lineaire que log !!!
                       ! Ca fait n'importe quoi en log en cas de flux faible
                       KJ_abs_interp = interp(KJ_absorbe_nRE(:,id), tab_nu, nu_bin(T) - nu_bin(T2))
                       A(T,T2,id) = KJ_abs_interp * delta_nu_bin(T)/ (nu_bin(T) - nu_bin(T2))
                    end do
                 end do

                 !! diagonal terms
                 do T=1,n_T
                    A(T,T,id) = -sum(A(T,:,id)) ! sum along k direction (diag. should be zero initially)
                 enddo

                 !! compute probabilities P(T) using iterative formulation of GD89

                 !! compute B array
                 B(:,:,id) = 0.
                 B(n_T,:,id) = A(n_T,:,id)
                 do T=n_T-1,1,-1
                    B(T,:,id) = A(T,:,id) + B(T+1,:,id)
                 enddo

                 !! compute P array
                 X(1,id) = 1.e-250_db
                 do T=2, n_T
                    if (X(T-1,id) <  1.e-300_db) X(T-1,id) = 1.e-300_db
                    if (X(T-1,id) >  1.e250_db) X(1:T-1,id) = X(1:T-1,id) * 1.0e-50_db ! permet de stabiliser en cas d'erreur d'arrondis
                    X(T,id) =  sum(B(T,1:T-1,id)*X(1:T-1,id)) / max(A(T-1, T,id),tiny_db)
                 enddo

                 !! Normalisation
                 X(:,id)=X(:,id)/sum(X(:,id))

                 do T=1, n_T ! boucle car passage tableau double -> simple bug avec ifort
                    Proba_Temperature(T,i,j,l) = X(T,id) ! P(T).dT, probability of being in bin of temperature T
                 enddo

              endif ! test : t_cool > t_abs
           endif ! densite_pouss > 0.
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel
  enddo !l

  deallocate(A, B, X)

  write(*,*) "Done"

  lconverged = .true.
  ! Verification convergence sur la temperature
  if (lRE_LTE) then
     do i=1, n_rad
        do j=1, nz
           delta_T = Temperature(i,j,1) - Temperature_old(i,j,1)
           if (delta_T > precision * Temperature_old(i,j,1)) then
              lconverged = .false.
           endif
        enddo
     enddo
     Temperature_old = Temperature
  endif

  if (lRE_nLTE) then
     do i=1, n_rad
        do j=1, nz
           do l=grain_RE_nLTE_start,grain_RE_nLTE_end
              delta_T = Temperature_1grain(i,j,l) - Temperature_1grain_old(i,j,l)
              if (delta_T > precision * Temperature_1grain_old(i,j,l)) then
                 lconverged = .false.
              endif
           enddo
        enddo
     enddo
     Temperature_1grain_old = Temperature_1grain
  endif

  do i=1, n_rad
     do j=1, nz
        do l=grain_nRE_start, grain_nRE_end
           if (l_RE(i,j,l)) then
              delta_T = Temperature_1grain_nRE(i,j,l) - Temperature_1grain_nRE_old(i,j,l)
              if (delta_T > precision * Temperature_1grain_nRE_old(i,j,l)) then
                 lconverged = .false.
              endif
              Temperature_1grain_nRE_old(i,j,l) =  Temperature_1grain_nRE(i,j,l)
           else
              Tpeak = 1
              maxP =  Proba_Temperature(1,i,j,l)
              do T=2, n_T
                 if (Proba_Temperature(T,i,j,l) > maxP) then
                    maxP=Proba_Temperature(T,i,j,l)
                    Tpeak=T
                 endif
              enddo

              if (Tpeak /= Tpeak_old(i,j,l)) then
                 lconverged=.false.
              endif

              if ( abs(maxP-maxP_old(i,j,l)) > precision * maxP_old(i,j,l) ) then
                 lconverged=.false.
              endif

              Tpeak_old(i,j,l) = Tpeak
              maxP_old(i,j,l) = maxP
           endif
        enddo !l
     enddo !j
  enddo !i


  !! TEST
!!$  l=grain_nRE_start ! TEST pour 1 taille de grains
!!$  do i=1, n_rad
!!$     do j=1, nz
!!$        somme1=0.0
!!$        somme2=0.0
!!$        !temp=Temperature_1grain(i,j,l) ! Tab_Temp(100)
!!$        do lambda=1, n_lambda
!!$           somme1 =   somme1 + q_abs(lambda,l)  * (sum(xJ_abs(lambda,i,j,:)) + J0(lambda,i,j,1) )
!!$           wl = tab_lambda(lambda)*1.0e-6
!!$           delta_wl = tab_delta_lambda(lambda)*1.0e-6
!!$           if (l_RE(i,j,l)) then
!!$              Temp = Temperature_1grain_nRE(i,j,l)
!!$              cst_wl=cst_th/(Temp*wl)
!!$              if (cst_wl < 500.) then
!!$                 somme2 = somme2 + q_abs(lambda,l)* 1.0/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
!!$              endif
!!$           else
!!$              do T=1, n_T
!!$                 Temp = tab_Temp(T)
!!$                 cst_wl=cst_th/(Temp*wl)
!!$                 if (cst_wl < 500.) then
!!$                    somme2 = somme2 + q_abs(lambda,l)* 1.0/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl * Proba_Temperature(T,i,j,l)
!!$                 endif !cst_wl
!!$              enddo
!!$           endif
!!$        enddo ! lambda
!!$        somme1=somme1 *n_phot_L_tot/volume(i)
!!$        somme2=somme2*2.0*hp*c_light**2
!!$        write(*,*) i,j,l_RE(i,j,l),  real(somme1/somme2)
!!$       ! read(*,*)
!!$     enddo !j
!!$  enddo !i
  !! FIN TEST

  return

end subroutine Temp_nRE

!**********************************************************************

subroutine emission_nRE()
  ! calcule la reemission des grains hors equilibre :
  ! - prob_E_cell
  ! - frac_E_stars a 0
  ! C. Pinte
  ! 05/02/07

  implicit none

  integer :: i, j, k, T, l, lambda, n_max, alloc_status
  real :: Temp, cst_wl, cst_wl_max, wl, delta_wl, fn
  real(kind=db) :: E_emise, frac_E_abs_nRE, Delta_E
  real(kind=db), dimension(:), allocatable :: E_cell, E_cell_old

  ! proba emission en fct lambda puis pour chaque lambda proba en fct position
  ! proba differentielle en fonction proba emission precedente

  !E_abs_nRE est en nbre de paquets !!!
  frac_E_abs_nRE =  E_abs_nRE / nbre_photons_tot
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

  if (l3D) then
     n_max = n_rad*2*nz*n_az
  else
     n_max = n_rad*nz
  endif

  allocate(E_cell(n_max), E_cell_old(n_max), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif

  spectre_emission_cumul(0) = 0.0
  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl = tab_delta_lambda(lambda)*1.e-6
     E_cell=0.0

     ! Emission par cellule des PAHs

     !$omp parallel default(none) &
     !$omp private(i,j,l,k,E_emise,Temp,cst_wl,T) &
     !$omp shared(lambda,wl,delta_wl,E_cell,E_cell_old,tab_lambda,tab_delta_lambda,grain_nRE_start,grain_nRE_end,n_max) &
     !$omp shared(n_rad,nz,l_RE, Temperature_1grain_nRE,n_T,q_abs,densite_pouss,volume,tab_Temp,Proba_Temperature) &
     !$omp shared(Emissivite_nRE_old,cst_wl_max)
     !$omp do
     do i=1,n_rad
        do j=1,nz
           ! Combinaison des 2 indices pour dichotomie
           l=j+nz*(i-1)
           E_emise = 0.0
           do k=grain_nRE_start,grain_nRE_end
              if (l_RE(i,j,k)) then ! le grain a une temperature
                 Temp = Temperature_1grain_nRE(i,j,k)
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    E_emise = E_emise + 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)* &
                         volume(i)/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
                 endif
              else ! densite de proba de Temperature
                 do T=1,n_T
                    temp=tab_Temp(T)
                    cst_wl=cst_th/(Temp*wl)
                    if (cst_wl < cst_wl_max) then
                       E_emise = E_emise + 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)* &
                            volume(i)/((wl**5)*(exp(cst_wl)-1.0)) * Proba_Temperature(T,i,j,k) &
                            * delta_wl
                    endif !cst_wl
                 enddo !T
              endif
           enddo !k
           E_cell(l) =   E_emise
           ! Recup et mise a jour de l'ancienne emmissivite
           E_cell_old(l) = Emissivite_nRE_old(lambda,i,j,1)
           if (E_cell(l) >   E_cell_old(l)) then
              Emissivite_nRE_old(lambda,i,j,1) = E_emise
           else ! Test en cas de pb d'arrondis
              ! On ne met pas a jour
              E_cell(l) = E_cell_old(l)
           endif
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel


     prob_E_cell(lambda,0)=0.0
     do l=1, n_max
        Delta_E = E_cell(l) - E_cell_old(l)
        prob_E_cell(lambda,l) = prob_E_cell(lambda,l-1) + Delta_E
     enddo

     spectre_emission_cumul(lambda) =  spectre_emission_cumul(lambda-1) + prob_E_cell(lambda,n_max)

     ! Les etoiles ont deja emis
     frac_E_stars(lambda) = 0.0

     if (prob_E_cell(lambda,n_max) > tiny_real) then
        prob_E_cell(lambda,:)=prob_E_cell(lambda,:)/prob_E_cell(lambda,n_max)
     else
        prob_E_cell(lambda,:)=0.0
     endif
  enddo ! lambda


  ! Normalisation de la proba d'emission / fct lambda des grains hors equilibre
  do lambda=1,n_lambda
     spectre_emission_cumul(lambda)=spectre_emission_cumul(lambda)/spectre_emission_cumul(n_lambda)
  enddo

  deallocate(E_cell, E_cell_old)

  return

end subroutine emission_nRE

!**********************************************************************

subroutine init_emissivite_nRE()
  ! initialise le tableau Emissivite_nRE_old avec l'emissivite a Tmin
  ! C. Pinte
  ! 05/02/07

  implicit none


  integer :: lambda, i, j, k
  real(kind=db) :: E_emise, facteur, cst_wl, wl
  real :: Temp, cst_wl_max, delta_wl

  cst_wl_max = log(huge_real)-1.0e-4
  temp=tab_Temp(1)

  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl = tab_delta_lambda(lambda)*1.e-6
     cst_wl=cst_th/(Temp*wl)

     ! Emission par cellule des PAHs
     do i=1,n_rad
        if (cst_wl < 100.0) then
           facteur = volume(i)/((wl**5)*(exp(cst_wl)-1.0)) * delta_wl
        else
           facteur = 0.0_db
        endif

        do j=1,nz
           E_emise = 0.0_db
           do k=grain_nRE_start,grain_nRE_end
              E_emise = E_emise + 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)* facteur !* Proba_Temperature = 1 pour Tmin
           enddo !k
           Emissivite_nRE_old(lambda,i,j,1) = E_emise
        enddo !j
     enddo !i

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

  integer :: i, j, pk, jj, k, n_max, l, T, alloc_status
  real(kind=db) :: Temp, wl, cst_wl, E_star, surface, Ener, frac, E_emise, cst_wl_max
  real(kind=db) :: delta_T
  real(kind=db), dimension(:), allocatable :: E_cell, E_cell_corrected

  cst_wl_max = log(huge_real)-1.0e-4

  wl = tab_lambda(lambda)*1.e-6

  ! Emission totale des etoiles
  E_star = E_stars(lambda)

  if (l3D) then
     n_max = n_rad*2*nz*n_az
  else
     n_max = n_rad*nz
  endif

  allocate(E_cell(n_max), E_cell_corrected(n_max), stat=alloc_status)
   if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  E_cell=0.0
  E_cell_corrected=0.0

  ! Cas LTE
  if (lRE_LTE) then
     do i=1,n_rad
        bz : do j=j_start,nz
           if (j==0) cycle bz
           if (l3D) then
              if (j < 0) then
                 jj = j + nz + 1
              else
                 jj = j + nz
              endif
           endif
           do pk=1, n_az
              if (l3D) then
                 ! Combinaison des 3 indices pour dichotomie
                 l= pk+n_az*(jj+2*nz*(i-1)-1)

              else
                 ! Combinaison des 2 indices pour dichotomie
                 l= j+nz*(i-1)
              endif

              Temp=Temperature(i,j,pk)
              if (Temp < tiny_real) then
                 ! prob_E_cell(lambda,k) =  prob_E_cell(lambda,k-1)
              else
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    if (.not.test_dark_zone(i,j,pk,0.0_db,0.0_db)) then
                       E_cell(l) =  4.0*kappa_abs_eg(lambda,i,j,pk)*volume(i)/((wl**5)*(exp(cst_wl)-1.0))
                    endif
                 endif !cst_wl
              endif ! Temp==0.0

              if (lweight_emission) then
                 E_cell_corrected(l) = E_cell(l) * weight_proba_emission(i,j)
              else
                 E_cell_corrected(l) = E_cell(l)
              endif

           enddo !pk
        enddo bz !j
     enddo !i
  endif



  ! Cas nLTE
  if (lRE_nLTE) then
     do i=1,n_rad
        do j=1,nz
           ! Combinaison des 2 indices pour dichotomie
           l=j+nz*(i-1)
           E_emise = 0.0
           do k=grain_RE_nLTE_start,grain_RE_nLTE_end
              Temp=Temperature_1grain(i,j,k)
              if (Temp < tiny_real) then
                 !E_emise = E_emise + 0.0
              else
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    if (i==1) then
                       Ener = 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)*volume(i)/((wl**5)*(exp(cst_wl)-1.0))
                       frac = (r_in_opacite(j,1)-rmin)/(r_lim(1)-rmin)
                       E_emise = E_emise + Ener * frac
                    else if (.not.test_dark_zone(i,j,1,0.0_db,0.0_db)) then
                       E_emise = E_emise +   4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)*volume(i)/((wl**5)*(exp(cst_wl)-1.0))
                    endif
                 endif !cst_wl
              endif ! Temp==0.0
           enddo !k
           E_cell(l) = E_emise
           if (lweight_emission) then
              E_cell_corrected(l) = E_cell(l) * weight_proba_emission(i,j)
           else
              E_cell_corrected(l) = E_cell(l)
           endif

        enddo !j
     enddo !i
  endif

  ! Cas nRE
  if (lnRE) then

     ! Initialisation du tableau de temperature pour images avec grains
     ! hors equilibre
     if (lmono0) then
        delta_T=exp((1.0_db/(real(n_T,kind=db)))*log(T_max/T_min))
        tab_Temp(1)=T_min*sqrt(delta_T)
        do t=2,n_T
           tab_Temp(t)=delta_T*tab_Temp(t-1)
        enddo
     endif


     do i=1,n_rad
        do j=1,nz
           ! Combinaison des 2 indices pour dichotomie
           l=j+nz*(i-1)
           E_emise = 0.0
           do k=grain_nRE_start,grain_nRE_end
              if (l_RE(i,j,k)) then ! le grain a une temperature
                 temp=Temperature_1grain_nRE(i,j,k)
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < cst_wl_max) then
                    if (i==1) then
                       Ener = 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)*volume(i)/((wl**5)* &
                            (exp(cst_wl)-1.0))
                       frac = (r_in_opacite(j,1)-rmin)/(r_lim(1)-rmin)
                       E_emise = E_emise + Ener * frac
                    else if (.not.test_dark_zone(i,j,1,0.0_db,0.0_db)) then
                       E_emise = E_emise + 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)* &
                            volume(i)/((wl**5)*(exp(cst_wl)-1.0))
                    endif
                 endif !cst_wl
              else ! la grain a une proba de T
                 do T=1,n_T
                    temp=tab_Temp(T)
                    cst_wl=cst_th/(Temp*wl)
                    if (cst_wl < cst_wl_max) then
                       if (i==1) then
                          Ener = 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)*volume(i)/((wl**5)* &
                               (exp(cst_wl)-1.0)) * Proba_Temperature(T,i,j,k)
                          frac = (r_in_opacite(j,1)-rmin)/(r_lim(1)-rmin)
                          E_emise = E_emise + Ener * frac
                       else if (.not.test_dark_zone(i,j,1,0.0_db,0.0_db)) then
                          E_emise = E_emise + 4.0*q_abs(lambda,k)*densite_pouss(i,j,1,k)* &
                               volume(i)/((wl**5)*(exp(cst_wl)-1.0)) * Proba_Temperature(T,i,j,k)
                       endif
                    endif !cst_wl
                 enddo !T
              endif ! l_RE
           enddo !k
           E_cell(l) =  E_cell(l) + E_emise
           if (lweight_emission) then
              E_cell_corrected(l) = E_cell(l) * weight_proba_emission(i,j)
           else
              E_cell_corrected(l) = E_cell(l)
           endif
        enddo !j
     enddo !i

  endif


  E_disk(lambda) = sum(E_cell)
  frac_E_stars(lambda)=E_star/(E_star+E_disk(lambda))

  ! Energie totale emise a une distance emise egale a la distance terre-etoile
  ! on chosit cette distance pour calibrer le flux / pi*B(lambda)

  ! Normalisation
  ! Pas de besoin de cst : r_etoile**2 en Au2,
  ! kappa_abs en AU-1, volume en Au**3
  surface=4*pi*(pc_to_AU*distance)**2
  if (l_sym_centrale) then
     E_totale(lambda) = 2.0*pi*hp*c_light**2/surface * (E_star+E_disk(lambda)) * real(N_thet)*real(N_phi)
  else
     E_totale(lambda) = 2.0*pi*hp*c_light**2/surface * (E_star+E_disk(lambda)) * real(2*N_thet)*real(N_phi)
  endif


  ! Distribution (spatiale) cumulee d'energie
  prob_E_cell(lambda,0)=0.0
  do l=1, n_max
     prob_E_cell(lambda,l) = prob_E_cell(lambda,l-1) + E_cell_corrected(l)
  enddo


  ! Normalisation du facteur de correction
  if (lweight_emission) then
     correct_E_emission = correct_E_emission * prob_E_cell(lambda,n_max) / E_disk(lambda)
  endif

  if (prob_E_cell(lambda,n_max) > tiny_real) then
     prob_E_cell(lambda,:)=prob_E_cell(lambda,:)/prob_E_cell(lambda,n_max)
  else
     prob_E_cell(lambda,:)=0.0
  endif

  return

end subroutine repartition_energie

!**********************************************************************

subroutine chauffage_interne()
  ! Calcule le temperature initiale du disque avant le step 1
  ! C. Pinte
  ! 27/02/06

  implicit none

  integer :: lambda, i, j
  real :: wl, delta_wl, cst_wl, cst_wl_max


  ! Energie venant de l'equilibre avec nuage à T_min
  E0(:,:,:) = exp(log_frac_E_em(:,:,:,1))

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

  return

end subroutine chauffage_interne

!**********************************************************************

end module thermal_emission
