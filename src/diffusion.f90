module diffusion

  use parametres
  use constantes
  use dust_prop
  use messages
  use wavelengths
  use thermal_emission
  use temperature
  use cylindrical_grid

  implicit none
  save

  contains

subroutine setDiffusion_coeff(i)
  ! calcule les coefficients de diffusion pour
  ! les cellules dans la zone ou est utilisee
  ! l'approximation de diffusion
  ! C. Pinte
  ! 15/02/07

  integer, intent(in) :: i

  real(kind=dp) :: cst_Dcoeff, wl, delta_wl, cst, cst_wl, coeff_exp, dB_dT, Temp, somme
  integer :: j, k, lambda, icell, p_icell

  real(kind=dp), parameter :: precision = 1.0e-1_dp ! Variation de temperature au dela de laquelle le coeff de diff et mis a jour
  ! le mettre a 0, evite le drole de BUG

  cst_Dcoeff = pi/(12.*sigma)

  p_icell = icell_ref

  do k=1,n_az
     do j=j_start, nz
        if (j==0) cycle
        icell = cell_map(i,j,k)
        if (lvariable_dust) p_icell = icell
        !Temp=Tdust(i,j,k)
        if (abs(DensE(i,j,k) - DensE_m1(i,j,k)) > precision * DensE_m1(i,j,k)) then
           ! On met a jour le coeff
       !    write(*,*) "update", i,j, DensE(i,j,k), DensE_m1(i,j,k)
           Temp = DensE(i,j,k)**0.25

           cst=cst_th/Temp
           somme=0.0_dp
           do lambda=1, n_lambda
              ! longueur d'onde en metre
              wl = tab_lambda(lambda)*1.e-6
              delta_wl=tab_delta_lambda(lambda)*1.e-6
              cst_wl=cst/wl
              if (cst_wl < 200.0) then
                 coeff_exp=exp(cst_wl)
                 dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
              else
                 dB_dT = 0.0_dp
              endif
              somme = somme + dB_dT/(kappa(p_icell,lambda) * kappa_factor(icell))  * delta_wl
           enddo
           ! kappa_R = 4.*sigma * Temp**3 / (pi * somme)
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


  integer, intent(in) :: i

  real(kind=dp) :: cst_Dcoeff, wl, delta_wl, cst, cst_wl, coeff_exp, dB_dT, Temp, somme
  integer :: j, k, lambda, icell, p_icell

  cst_Dcoeff = pi/(12.*sigma)

  p_icell = icell_ref

  do k=1,n_az
     do j=j_start,nz
        if (j==0) cycle
        icell = cell_map(i,j,k)
        if (lvariable_dust) p_icell = icell
        Temp=Tdust(icell)
        cst=cst_th/Temp
        somme=0.0_dp
        do lambda=1, n_lambda
           ! longueur d'onde en metre
           wl = tab_lambda(lambda)*1.e-6
           delta_wl=tab_delta_lambda(lambda)*1.e-6
           cst_wl=cst/wl
           if (cst_wl < 200.0) then
              coeff_exp=exp(cst_wl)
              dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
           else
              dB_dT = 0.0_dp
           endif
           somme = somme + dB_dT/(kappa(p_icell,lambda)*kappa_factor(icell)) * delta_wl
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

  integer, intent(in) :: ri

  integer :: j, k, icell

  do j=j_start,nz
     if (j==0) cycle
     do k=1,n_az
        icell = cell_map(ri,j,k)
        DensE(ri,j,k) =  Tdust(icell)**4 ! On se tappe de la constante non ??? 4.*sigma/c
     enddo
  enddo

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

  integer :: i, j, k

  do k=1,n_az
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1,zj_sup_dark_zone(i,1) + delta_cell_dark_zone
           Tdust(cell_map(i,j,k)) = DensE(i,j,k)**0.25
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

  integer :: i, j, k

  do k=1,n_az
     do i=ri_in_dark_zone(k), ri_out_dark_zone(k)
        do j=1,zj_sup_dark_zone(i,1)
           Tdust(cell_map(i,j,k)) = T_min
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

  real, dimension(n_cells) :: Temp0
  real :: max_delta_E_r, stabilite, precision
  integer :: n_iter, i
  logical :: lconverged


  write(*,*) "Computing 2D diffusion approx. in central parts of the disk"

  ! Pour initier la premiere boucle
  lconverged=.false.

  precision = 1.0e-6
  stabilite=10.


  Tdust_old = T_min
  Temp0 = Tdust

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
           Tdust = Temp0
           Tdust_old = T_min
           exit infinie
        endif

        ! C'est debile : temp pas definie
!        Tdust_old = Tdust

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

  real :: max_delta_E_r, stabilite, precision
  integer :: n_iter, i, k
  logical :: lconverged


  write(*,*) "Computing 1+1D diffusion approx. in central parts of the disk"

  call clean_temperature()

  ! Pour initier la premiere boucle
  !-- Tdust_old = T_min
  !-- Temp0 = Tdust

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

        if (stabilite < 0.01) call error("diffusion approximation does not seem to converge")

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
          !--    Tdust(i,:,:) = Temp0(i,:,:)
          !--    Tdust_old = T_min
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

  real, intent(in) :: stabilite
  real, intent(out) :: max_delta_E_r
  logical, intent(out) :: lconverge

  real(kind=dp) :: dt, dE_dr_m1, dE_dr_p1, d2E_dr2, d2E_dz2
  real(kind=dp) :: D_Laplacien_E, dr, dz, delta_E, delta_E_r
  real(kind=dp), dimension(n_rad,nz,n_az) :: tab_dt
  integer :: i,j,k

  lconverge = .true.

  ! Determination du pas de temps pour respecter le critere de stabilite

  ! pas de temps pour chacune des cellules
  tab_dt = huge_dp ! pour ne pas selectionner les cellules hors zone de diff
  do k=1, n_az
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1, zj_sup_dark_zone(i,k) + delta_cell_dark_zone
           dr = r_grid(cell_map(i,j,1))-r_grid(cell_map(i-1,j,1))
           dz = delta_z(i,j)
           ! tab_dt(i,j,k) = min(dr,dz)**2/Dcoeff(i,j,k)
           tab_dt(i,j,k) = 1.0_dp/(Dcoeff(i,j,k)*(1.0_dp/dr**2 + 1.0_dp/dz**2))
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
     !$omp shared(DensE,dt,delta_z,k,n_rad,cell_map)
     !$omp do schedule(dynamic,10)
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1, zj_sup_dark_zone(i,k) + delta_cell_dark_zone

           ! Calcul du Laplacien en cylindrique
           ! Attention D rentre dans le laplacien car il depend de laposition
           ! Ne marche qu'en 2D pour le moment
           dE_dr_m1 = (DensE_m1(i,j,k) - DensE_m1(i-1,j,k))/(r_grid(cell_map(i,j,k))-r_grid(cell_map(i-1,j,k)))
           dE_dr_p1 = (DensE_m1(i+1,j,k) - DensE_m1(i,j,k))/(r_grid(cell_map(i+1,j,k))-r_grid(cell_map(i,j,k)))

           !    frac=(log(r_lim(i))-log(r_grid(i)))/(log(r_grid(i+1))-log(r_grid(i)))
           !    Dcoeff_p=exp(log(Dcoeff(i,j,k))*frac+log(Dcoeff(i+1,j,k))*(1.0_dp-frac))

           !   frac=(log(r_lim(i-1))-log(r_grid(i-1)))/(log(r_grid(i))-log(r_grid(i-1)))
           !   Dcoeff_m=exp(log(Dcoeff(i-1,j,k))*frac+log(Dcoeff(i,j,k))*(1.0_dp-frac))


           !Dcoeff_p = 0.5_dp * (Dcoeff(i,j,k) + Dcoeff(i+1,j,k))
           !Dcoeff_m = 0.5_dp * (Dcoeff(i,j,k) + Dcoeff(i-1,j,k))

           !Dcoeff_p = 0.5_dp * (Dcoeff(i,j,k) + interp(Dcoeff(i+1,:,k),z_grid(i+1,:),z_grid(i,j)))
           !Dcoeff_m = 0.5_dp * (Dcoeff(i,j,k) + interp(Dcoeff(i-1,:,k),z_grid(i-1,:),z_grid(i,j)))

           !d2E_dr2  = (dE_dr_p1*Dcoeff_p  - dE_dr_m1*Dcoeff_m) / (2._dp*(r_grid(i+1,j)-r_grid(i-1,j)))

           d2E_dr2  =  Dcoeff(i,j,k) * (dE_dr_p1 - dE_dr_m1) /(2._dp*(r_grid(cell_map(i+1,j,k))-r_grid(cell_map(i-1,j,k))))

           !Dcoeff_p = 0.5_dp * (Dcoeff(i,j,k) + Dcoeff(i,j+1,k))
           !Dcoeff_m = 0.5_dp * (Dcoeff(i,j,k) + Dcoeff(i-1,j-1,k))

           !dE_dz_p1 = (DensE_m1(i,j+1,k) - DensE_m1(i,j,k))
           !dE_dz_m1 = (DensE_m1(i,j,k) - DensE_m1(i,j-1,k))

           !d2E_dz2  = (dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m) / (2.0_dp * delta_z(i)**2)

           d2E_dz2  =   Dcoeff(i,j,k) * (DensE_m1(i,j+1,k) + DensE_m1(i,j-1,k) - 2.0 * DensE(i,j,k)) / (2.0 * delta_z(i,j)**2)


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

  integer, intent(in) :: ri
  real, intent(in) :: stabilite
  real, intent(out) :: max_delta_E_r
  logical, intent(out) :: lconverge

  real(kind=dp) :: dt, dE_dz_m1, dE_dz_p1, d2E_dz2
  real(kind=dp) :: D_Laplacien_E, dz, delta_E, delta_E_r
  real(kind=dp), dimension(nz) :: tab_dt
  integer :: j,k

  lconverge = .true.
  k=1

  ! Determination du pas de temps pour respecter le critere de stabilite

  ! pas de temps pour chacune des cellules
  tab_dt = huge_dp ! pour ne pas selectionner les cellules hors zone de diff
  do j=1, zj_sup_dark_zone(ri,k) + delta_cell_dark_zone
     dz = delta_z(ri,j)
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
    ! Dcoeff_p = 0.5_dp * (Dcoeff(ri,j,k) + Dcoeff(ri,j+1,k))
    ! Dcoeff_m = 0.5_dp * (Dcoeff(ri,j,k) + Dcoeff(ri-1,j-1,k))
     ! Plus stable
  !   Dcoeff_p = Dcoeff(ri,j,k)
  !   Dcoeff_m =  Dcoeff(ri,j,k)
  !   d2E_dz2  = (dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m) / (2.0_dp * delta_z(ri)**2)

     d2E_dz2  =  Dcoeff(ri,j,k) * (dE_dz_p1 - dE_dz_m1) / (2.0_dp * delta_z(ri,j)**2)


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

  integer :: i,j,k, icell

  if (lRE_nLTE) then
     do k=1,n_az
        do i=ri_in_dark_zone(k), ri_out_dark_zone(k)
           do j=1,zj_sup_dark_zone(i,1)
              icell = cell_map(i,j,k)
              Tdust_1grain(:,icell) = Tdust(icell)
           enddo !j
        enddo !j
     enddo !k
  endif

  if (lnRE) then
      do k=1,n_az
        do i=ri_in_dark_zone(k), ri_out_dark_zone(k)
           do j=1,zj_sup_dark_zone(i,1)
              icell = cell_map(i,j,k)
              Tdust_1grain_nRE(:,icell) = Tdust(icell)
              l_RE(:,icell) = .true.
           enddo !j
        enddo !j
     enddo !k
  endif

  return

end subroutine diffusion_approx_nLTE_nRE

!************************************************************

subroutine diffusion_opacity(icell, Planck_opacity,rec_Planck_opacity)
  ! Compute current Planck reciprocal mean opacity for all cells
  ! (note : diffusion coefficient needs to be defined with Rosseland opacity
  ! in B&W mode)
  ! Diffusion coefficient is D = 1/(rho * opacity)
  ! This opacity/diffusion coefficient includes scattering
  ! See Min et al 2009 and Robitaille et al 2010
  use parametres
  use constantes
  use wavelengths, only : n_lambda, tab_lambda, tab_delta_lambda
  use Temperature, only : Tdust
  use dust_prop, only : kappa
  use Voronoi_grid, only : Voronoi
  use cylindrical_grid, only : volume
  use density, only : masse_gaz, densite_gaz

  integer,  intent(in)  :: icell
  real(dp), intent(out) :: Planck_opacity,rec_Planck_opacity ! cm2/g (ie per gram of gas)

  integer :: lambda
  real(dp) :: somme, somme2, cst, cst_wl, B, dB_dT, coeff_exp, wl, delta_wl, norm, Temp

  integer, pointer :: p_icell
  integer, target :: icell0

  icell0 = icell

  temp = Tdust(icell)
  if ((temp > 1) .and. (Voronoi(icell)%original_id > 0)) then

     if (lvariable_dust) then
        p_icell => icell0
     else
        p_icell => icell_ref
     endif

     somme  = 0.0_dp
     somme2 = 0.0_dp
     norm = 0.0_dp
     cst    = cst_th/temp
     do lambda = 1,n_lambda
        ! longueur d'onde en metre
        wl       = tab_lambda(lambda)*1.e-6
        delta_wl = tab_delta_lambda(lambda)*1.e-6
        cst_wl   = cst/wl
        if (cst_wl < 200.0) then
           coeff_exp = exp(cst_wl)
           B = 1.0_dp/((wl**5)*(coeff_exp-1.0))*delta_wl
           !dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
        else
           B = 0.0_dp
           !dB_dT = 0.0_dp
        endif
        somme  = somme  + B/(kappa(p_icell,lambda) * kappa_factor(icell0))*delta_wl
        somme2  = somme2  + B * (kappa(p_icell,lambda) * kappa_factor(icell0))*delta_wl
        norm = norm + B*delta_wl
     enddo
     rec_Planck_opacity = norm/somme
     Planck_opacity = somme2/norm
  else
     rec_Planck_opacity = 0.
     Planck_opacity = 0.
  endif

end subroutine diffusion_opacity

!************************************************************

end module diffusion
