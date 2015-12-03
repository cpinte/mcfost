module diffusion

  use parametres
  use constantes
  use opacity
  use em_th

  implicit none

  contains

subroutine setDiffusion_coeff(i)
  ! calcule les coefficients de diffusion pour
  ! les cellules dans la zone ou est utilisee
  ! l'approximation de diffusion
  ! C. Pinte
  ! 15/02/07

  implicit none

  integer, intent(in) :: i

  real(kind=db) :: cst_Dcoeff, wl, delta_wl, cst, cst_wl, coeff_exp, dB_dT, Temp, somme
  integer :: j, k, lambda, icell

  real(kind=db), parameter :: precision = 1.0e-1_db ! Variation de temperature au dela de laquelle le coeff de diff et mis a jour
  ! le mettre a 0, evite le drole de BUG

  cst_Dcoeff = c_light*pi/(12.*sigma)

  do k=1,n_az
     do j=1, nz
        icell = cell_map(i,j,k)
        !Temp=Temperature(i,j,k)
        if (abs(DensE(i,j,k) - DensE_m1(i,j,k)) > precision * DensE_m1(i,j,k)) then
           ! On met a jour le coeff
       !    write(*,*) "update", i,j, DensE(i,j,k), DensE_m1(i,j,k)
           Temp = DensE(i,j,k)**0.25

           cst=cst_th/Temp
           somme=0.0_db
           do lambda=1, n_lambda
              ! longueur d'onde en metre
              wl = tab_lambda(lambda)*1.e-6
              delta_wl=tab_delta_lambda(lambda)*1.e-6
              cst_wl=cst/wl
              if (cst_wl < 200.0) then
                 coeff_exp=exp(cst_wl)
                 dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
              else
                 dB_dT = 0.0_db
              endif
              somme = somme + dB_dT/kappa(icell,lambda) * delta_wl
           enddo
           !kappa_R = 4.*sigma * Temp**3 / (pi * somme)
           ! Dcoeff = c_light/(3kappa_R) car kappa volumique
           Dcoeff(i,j,k) =  cst_Dcoeff * somme/Temp**3
        endif
     enddo
  enddo

  ! Condition limite
  Dcoeff(:,0,:) = Dcoeff(:,1,:)

  return

end subroutine setDiffusion_coeff

!************************************************************

subroutine setDiffusion_coeff0(i)
  ! calcule les coefficients de diffusion pour
  ! les cellules dans la zone ou est utilisee
  ! l'approximation de diffusion
  ! C. Pinte
  ! 15/02/07

  implicit none

  integer, intent(in) :: i

  real(kind=db) :: cst_Dcoeff, wl, delta_wl, cst, cst_wl, coeff_exp, dB_dT, Temp, somme
  integer :: j, k, lambda, icell

  cst_Dcoeff = c_light*pi/(12.*sigma)

  do k=1,n_az
     do j=1,nz
        icell = cell_map(i,j,k)
        Temp=Temperature(i,j,k)
        cst=cst_th/Temp
        somme=0.0_db
        do lambda=1, n_lambda
           ! longueur d'onde en metre
           wl = tab_lambda(lambda)*1.e-6
           delta_wl=tab_delta_lambda(lambda)*1.e-6
           cst_wl=cst/wl
           if (cst_wl < 200.0) then
              coeff_exp=exp(cst_wl)
              dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
           else
              dB_dT = 0.0_db
           endif
           somme = somme + dB_dT/kappa(icell,lambda) * delta_wl
        enddo
        ! kappa_R = 4.*sigma * Temp**3 / (pi * somme)
        ! Dcoeff = c_light/(3kappa_R) car kappa volumique
        Dcoeff(i,j,k) =  cst_Dcoeff * somme/Temp**3
     enddo
  enddo

  ! Condition limite
  ! Dcoeff(:,0,:) = Dcoeff(:,1,:)

  return

end subroutine setDiffusion_coeff0

!************************************************************

subroutine Temperature_to_DensE(ri)
  ! Calcule la densite d'energie des cellules a partir de leur
  ! temperature
  ! Converti toute la grille pour avoir les cellules sur le bord
  ! Execute une seule fois donc pas de soucis
  ! C. Pinte
  ! 15/02/07

  implicit none

  integer, intent(in) :: ri

  DensE(ri,1:nz,:) =  Temperature(ri,:,:)**4 ! On se tappe de la constante non ??? 4.*sigma/c

  ! Condition limite : pas de flux au niveau du plan median
  DensE(ri,0,:) = DensE(ri,1,:)

  DensE_m1(ri,:,:) = DensE(ri,:,:)

  return

end subroutine Temperature_to_DensE

!************************************************************

subroutine DensE_to_temperature()
  ! Calcule la temperature des cellules dans le zone de diffusion
  ! a partir de leur densite d'energie
  ! C. Pinte
  ! 15/02/07

  implicit none


  integer :: i, j, k

  do k=1,n_az
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1,zj_sup_dark_zone(i,1) + delta_cell_dark_zone
           Temperature(i,j,k) = DensE(i,j,k)**0.25
        enddo !j
     enddo !j
  enddo !k
  return

end subroutine DensE_to_temperature

!***********************************************************

subroutine clean_temperature()
  ! Met la temperature a Tmin dans la zone sombre
  ! C. Pinte
  ! 15/02/07

  implicit none

  integer :: i, j, k

  do k=1,n_az
     do i=ri_in_dark_zone(k), ri_out_dark_zone(k)
        do j=1,zj_sup_dark_zone(i,1)
           Temperature(i,j,k) = T_min
        enddo !j
     enddo !j
  enddo !k
  return

end subroutine clean_temperature

!***********************************************************

subroutine Temp_approx_diffusion()
  ! Calcul de la temperature dans les zones profondes du disque a l'aide
  ! d'une equation de diffusion (equation parabolique).
  ! Le coeff de diffusion depend de T donc pb non linaire
  ! => resolution par le schema explicite d'integration : la solution stationnaire
  ! est la limite a t --> inf de la solution non stationnaire
  ! C. Pinte
  ! 15/02/07

  implicit none

  real, dimension(n_rad,nz,1) :: Temp0
  real :: max_delta_E_r, stabilite, precision
  integer :: n_iter, i
  logical :: lconverged


  write(*,*) "Computing 2D diffusion approx. in central parts of the disk"

  ! Pour initier la premiere boucle
  lconverged=.false.

  precision = 1.0e-6
  stabilite=10.


  Temperature_old = T_min
  Temp0 = Temperature

  do while (.not.lconverged)

     ! Passage temperature -> densite d'energie
     do i=1, n_rad
        call temperature_to_DensE(i)
     enddo

     ! Calcul coeff de diffusion : condition limite
     do i=1,n_rad
        call setDiffusion_coeff0(i)
     enddo

     if (stabilite < 0.01) then
        write(*,*) "Error : diffusion approximation does not seem to converge"
        write(*,*) "Exiting"
        return
     endif

     ! Iterations
     n_iter = 0
     infinie : do
        n_iter = n_iter + 1

        ! Un pas de temps de l'equation de diffusion
        call iter_Temp_approx_diffusion(stabilite,max_delta_E_r,lconverged)

        !test divergence
        if (.not.lconverged) then
           stabilite = 0.5 * stabilite
           precision = 0.1 * precision
           Temperature = Temp0
           Temperature_old = T_min
           exit infinie
        endif

        ! C'est debile : temp pas definie
!        Temperature_old = Temperature

!        write(*,*) n_iter, max_delta_E_r, precision  !, maxval(DensE)

        ! Test_convergence
        if (max_delta_E_r < precision) exit infinie

        ! Calcul des nouveaux coeff de diffusion pour iteration suivante
        do i=1, n_rad
           call setDiffusion_coeff(i)
        enddo

     enddo infinie
  enddo

  ! Reppassage en temperature
  call DensE_to_temperature()

  write(*,*) "Temperature computed with diffusion approximation using", n_iter," iterations"

  return

end subroutine Temp_approx_diffusion

!************************************************************

subroutine Temp_approx_diffusion_vertical()
  ! Calcul de la temperature dans les zones profondes du disque a l'aide
  ! d'une equation de diffusion (equation parabolique).
  ! Le coeff de diffusion depend de T donc pb non linaire
  ! => resolution par le schema explicite d'integration : la solution stationnaire
  ! est la limite a t --> inf de la solution non stationnaire
  ! C. Pinte
  ! 15/02/07

  implicit none

  real :: max_delta_E_r, stabilite, precision
  integer :: n_iter, i, k
  logical :: lconverged


  write(*,*) "Computing 1+1D diffusion approx. in central parts of the disk"

  call clean_temperature()

  ! Pour initier la premiere boucle
  !-- Temperature_old = T_min
  !-- Temp0 = Temperature

  k=1 ! variable azimuth
  !$omp parallel &
  !$omp default(none) &
  !$omp shared(k,ri_in_dark_zone,ri_out_dark_zone,n_rad) &
  !$omp private(i,lconverged,n_iter,precision,stabilite,max_delta_E_r)
  !$omp do schedule(dynamic,1)
  do i=max(ri_in_dark_zone(k) - delta_cell_dark_zone,3) , min(ri_out_dark_zone(k) + delta_cell_dark_zone,n_rad-2)

     lconverged=.false.
     precision = 1.0e-6
     stabilite=2.

     !-- do while (.not.lconverged)
        ! Passage temperature -> densite d'energie
        call temperature_to_DensE(i)

        ! Calcul coeff de diffusion : condition limite
        call setDiffusion_coeff0(i)

        if (stabilite < 0.01) then
           write(*,*) "Error : diffusion approximation does not seem to converge"
           write(*,*) "Exiting"
           stop
        endif

        ! Iterations
        n_iter = 0
        infinie : do
           n_iter = n_iter + 1

           ! Un pas de temps de l'equation de diffusion
           call iter_Temp_approx_diffusion_vertical(i,stabilite,max_delta_E_r,lconverged)

           !test divergence
          !-- if (.not.lconverged) then
          !--    stabilite = 0.5 * stabilite
          !--    precision = 0.1 * precision
          !--    Temperature(i,:,:) = Temp0(i,:,:)
          !--    Temperature_old = T_min
          !--    exit infinie
          !-- endif

!           write(*,*) n_iter, max_delta_E_r, precision  !, maxval(DensE)

           ! Test_convergence
           if (max_delta_E_r < precision) exit infinie

           ! Calcul des nouveaux coeff de diffusion pour iteration suivante
           call setDiffusion_coeff(i)

        enddo infinie
     !-- enddo
     write(*,*) "Radius", i  ,"/", ri_out_dark_zone(k) + delta_cell_dark_zone, "   T computed using", n_iter," iterations"
  enddo !i
  !$omp end do
  !$omp end parallel

  write(*,*) "Done"

  ! Reppassage en temperature
  call DensE_to_temperature()

  return

end subroutine Temp_approx_diffusion_vertical

!***********************************************************

subroutine iter_Temp_approx_diffusion(stabilite,max_delta_E_r,lconverge)
  ! Realise une iteration du schema explicite d'integration
  ! de l'equation de diffusion sur la densite d'energie des cellules
  ! ie : un pas de temps de l'equation non stationnaire
  ! C. Pinte
  ! 15/02/07

  implicit none

  real, intent(in) :: stabilite
  real, intent(out) :: max_delta_E_r
  logical, intent(out) :: lconverge

  real(kind=db) :: dt, dE_dr_m1, dE_dr_p1, d2E_dr2, d2E_dz2
  real(kind=db) :: D_Laplacien_E, dr, dz, delta_E, delta_E_r
  real(kind=db), dimension(n_rad,nz,n_az) :: tab_dt
  integer :: i,j,k

  lconverge = .true.

  ! Determination du pas de temps pour respecter le critere de stabilite

  ! pas de temps pour chacune des cellules
  tab_dt = huge_db ! pour ne pas selectionner les cellules hors zone de diff
  do k=1, n_az
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1, zj_sup_dark_zone(i,k) + delta_cell_dark_zone
           dr = r_grid(i,j)-r_grid(i-1,j)
           dz = delta_z(i)
           ! tab_dt(i,j,k) = min(dr,dz)**2/Dcoeff(i,j,k)
           tab_dt(i,j,k) = 1.0_db/(Dcoeff(i,j,k)*(1.0_db/dr**2 + 1.0_db/dz**2))
        enddo !j
     enddo !i
  enddo !k

  ! On prend le mini + un facteur de securite
  dt = stabilite * 0.5 * minval(tab_dt)

 ! write(*,*) "dt", dt

  ! Sauvegarde densite d'energie
  DensE_m1 = DensE

  ! Boucle sur les celules de la zone de diffusion
  max_delta_E_r = 0.0

  ! TODO : Boucle croisee a chaque iteration pour mieux lisser le bruit ??


  do k=1, n_az
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,dE_dr_m1,dE_dr_p1,d2E_dr2,delta_E_r,D_Laplacien_E,delta_E,d2E_dz2) &
     !$omp shared(DensE_m1,r_grid,z_grid,Dcoeff,ri_in_dark_zone,ri_out_dark_zone,zj_sup_dark_zone,max_delta_E_r) &
     !$omp shared(DensE,dt,delta_z,k,n_rad)
     !$omp do schedule(dynamic,10)
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1, zj_sup_dark_zone(i,k) + delta_cell_dark_zone

           ! Calcul du Laplacien en cylindrique
           ! Attention D rentre dans le laplacien car il depend de laposition
           ! Ne marche qu'en 2D pour le moment
           dE_dr_m1 = (DensE_m1(i,j,k) - DensE_m1(i-1,j,k))/(r_grid(i,j)-r_grid(i-1,j))
           dE_dr_p1 = (DensE_m1(i+1,j,k) - DensE_m1(i,j,k))/(r_grid(i+1,j)-r_grid(i,j))

           !    frac=(log(r_lim(i))-log(r_grid(i)))/(log(r_grid(i+1))-log(r_grid(i)))
           !    Dcoeff_p=exp(log(Dcoeff(i,j,k))*frac+log(Dcoeff(i+1,j,k))*(1.0_db-frac))

           !   frac=(log(r_lim(i-1))-log(r_grid(i-1)))/(log(r_grid(i))-log(r_grid(i-1)))
           !   Dcoeff_m=exp(log(Dcoeff(i-1,j,k))*frac+log(Dcoeff(i,j,k))*(1.0_db-frac))


           !Dcoeff_p = 0.5_db * (Dcoeff(i,j,k) + Dcoeff(i+1,j,k))
           !Dcoeff_m = 0.5_db * (Dcoeff(i,j,k) + Dcoeff(i-1,j,k))

           !Dcoeff_p = 0.5_db * (Dcoeff(i,j,k) + interp(Dcoeff(i+1,:,k),z_grid(i+1,:),z_grid(i,j)))
           !Dcoeff_m = 0.5_db * (Dcoeff(i,j,k) + interp(Dcoeff(i-1,:,k),z_grid(i-1,:),z_grid(i,j)))

           !d2E_dr2  = (dE_dr_p1*Dcoeff_p  - dE_dr_m1*Dcoeff_m) / (2._db*(r_grid(i+1,j)-r_grid(i-1,j)))

           d2E_dr2  =  Dcoeff(i,j,k) * (dE_dr_p1 - dE_dr_m1) /(2._db*(r_grid(i+1,j)-r_grid(i-1,j)))

           !Dcoeff_p = 0.5_db * (Dcoeff(i,j,k) + Dcoeff(i,j+1,k))
           !Dcoeff_m = 0.5_db * (Dcoeff(i,j,k) + Dcoeff(i-1,j-1,k))

           !dE_dz_p1 = (DensE_m1(i,j+1,k) - DensE_m1(i,j,k))
           !dE_dz_m1 = (DensE_m1(i,j,k) - DensE_m1(i,j-1,k))

           !d2E_dz2  = (dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m) / (2.0_db * delta_z(i)**2)

           d2E_dz2  =   Dcoeff(i,j,k) * (DensE_m1(i,j+1,k) + DensE_m1(i,j-1,k) - 2.0 * DensE(i,j,k)) / (2.0 * delta_z(i)**2)


           ! Laplacien
           D_Laplacien_E = d2E_dr2 + d2E_dz2

           ! On avance d'un pas de temps
           delta_E =  D_Laplacien_E * dt

           DensE(i,j,k) = DensE_m1(i,j,k) + delta_E

           if (DensE(i,j,k) < 0.) write(*,*) DensE(i,j,k), DensE_m1(i,j,k), delta_E

           ! Augmentation relative de la densite d'energie
           delta_E_r = delta_E/DensE(i,j,k)
           if (delta_E_r > max_delta_E_r) max_delta_E_r = delta_E_r
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel
  enddo !k

  ! Condition limite : pas de flux au niveau du plan median
  DensE(:,0,:) = DensE(:,1,:)


  if (maxval(DensE) > 1.0e30) then
     write(*,*) "Diffusion approximation is diverging"
     write(*,*) "Reducing time step"
     lconverge = .false.
  endif

  return

end subroutine iter_Temp_approx_diffusion

!************************************************************

subroutine iter_Temp_approx_diffusion_vertical(ri,stabilite,max_delta_E_r,lconverge)
  ! Realise une iteration du schema explicite d'integration
  ! de l'equation de diffusion sur la densite d'energie des cellules
  ! ie : un pas de temps de l'equation non stationnaire
  ! C. Pinte
  ! 15/02/07

  implicit none

  integer, intent(in) :: ri
  real, intent(in) :: stabilite
  real, intent(out) :: max_delta_E_r
  logical, intent(out) :: lconverge

  real(kind=db) :: dt, dE_dz_m1, dE_dz_p1, d2E_dz2
  real(kind=db) :: D_Laplacien_E, dz, delta_E, delta_E_r
  real(kind=db), dimension(nz) :: tab_dt
  integer :: j,k

  lconverge = .true.
  k=1

  ! Determination du pas de temps pour respecter le critere de stabilite

  ! pas de temps pour chacune des cellules
  tab_dt = huge_db ! pour ne pas selectionner les cellules hors zone de diff
  do j=1, zj_sup_dark_zone(ri,k) + delta_cell_dark_zone
     dz = delta_z(ri)
     tab_dt(j) = dz**2/Dcoeff(ri,j,k)
  enddo !j

  ! On prend le mini + un facteur de securite
  dt = stabilite * 0.5 * minval(tab_dt)


!  write(*,*) "dt", dt


  ! Sauvegarde densite d'energie
  DensE_m1(ri,:,:) = DensE(ri,:,:)

  ! Boucle sur les celules de la zone de diffusion
  max_delta_E_r = 0.0

  do j=1, zj_sup_dark_zone(ri,k) + delta_cell_dark_zone

     dE_dz_p1 = (DensE_m1(ri,j+1,k) - DensE_m1(ri,j,k))
     dE_dz_m1 = (DensE_m1(ri,j,k) - DensE_m1(ri,j-1,k))


     ! Ca rend le truc instable
    ! Dcoeff_p = 0.5_db * (Dcoeff(ri,j,k) + Dcoeff(ri,j+1,k))
    ! Dcoeff_m = 0.5_db * (Dcoeff(ri,j,k) + Dcoeff(ri-1,j-1,k))
     ! Plus stable
  !   Dcoeff_p = Dcoeff(ri,j,k)
  !   Dcoeff_m =  Dcoeff(ri,j,k)
  !   d2E_dz2  = (dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m) / (2.0_db * delta_z(ri)**2)

     d2E_dz2  =  Dcoeff(ri,j,k) * (dE_dz_p1 - dE_dz_m1) / (2.0_db * delta_z(ri)**2)


  !   write(*,*) "****"
  !   write(*,*) dE_dz_p1, dE_dz_m1
  !   write(*,*) Dcoeff_p, Dcoeff_m
  !   write(*,*) j, dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m, dE_dz_p1*Dcoeff_p

     ! Laplacien
     D_Laplacien_E =  d2E_dz2

     ! On avance d'un pas de temps
     delta_E =  D_Laplacien_E * dt

   !  write(*,*) "DeltaE", j, delta_E

     DensE(ri,j,k) = DensE_m1(ri,j,k) + delta_E

     ! Augmentation relative de la densite d'energie
     delta_E_r = delta_E/DensE(ri,j,k)
     if (delta_E_r > max_delta_E_r) max_delta_E_r = delta_E_r
  enddo !j

  ! Condition limite : pas de flux au niveau du plan median
  DensE(ri,0,:) = DensE(ri,1,:)

  if (maxval(DensE) > 1.0e30) then
     write(*,*) "Diffusion approximation is diverging"
     write(*,*) "Reducing time step"
     lconverge = .false.
  endif

  return

end subroutine iter_Temp_approx_diffusion_vertical

!************************************************************

subroutine diffusion_approx_nLTE_nRE()

  integer :: i,j,k

  if (lRE_nLTE) then
     do k=1,n_az
        do i=ri_in_dark_zone(k), ri_out_dark_zone(k)
           do j=1,zj_sup_dark_zone(i,1)
              Temperature_1grain(i,j,:) = Temperature(i,j,k)
           enddo !j
        enddo !j
     enddo !k
  endif

  if (lnRE) then
      do k=1,n_az
        do i=ri_in_dark_zone(k), ri_out_dark_zone(k)
           do j=1,zj_sup_dark_zone(i,1)
              Temperature_1grain_nRE(i,j,:) = Temperature(i,j,k)
              l_RE(i,j,:) = .true.
           enddo !j
        enddo !j
     enddo !k
  endif

  return

end subroutine diffusion_approx_nLTE_nRE

end module diffusion
