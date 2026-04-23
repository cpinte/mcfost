module diffusion

  use parameters
  use constants
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
  ! Computes the diffusion coefficients for
  ! the cells in the zone where the
  ! diffusion approximation is used
  ! C. Pinte
  ! 15/02/07

  integer, intent(in) :: i

  real(kind=dp) :: cst_Dcoeff, wl, delta_wl, cst, cst_wl, coeff_exp, dB_dT, Temp, total_sum
  integer :: j, k, lambda, icell, p_icell

  real(kind=dp), parameter :: precision = 1.0e-1_dp ! Temperature variation beyond which the diffusion coefficient is updated
  ! Setting it to 0 avoids a peculiar BUG

  cst_Dcoeff = pi/(12.*sigma)

  p_icell = icell1

  do k=1,n_az
     do j=j_start, nz
        if (j==0) cycle
        icell = cell_map(i,j,k)
        if (lvariable_dust) p_icell = icell
        !Temp=Tdust(i,j,k)
        if (abs(DensE(i,j,k) - DensE_m1(i,j,k)) > precision * DensE_m1(i,j,k)) then
           ! Update the coefficient
       !    write(*,*) "update", i,j, DensE(i,j,k), DensE_m1(i,j,k)
           Temp = DensE(i,j,k)**0.25

           cst=thermal_const/Temp
           total_sum=0.0_dp
           do lambda=1, n_lambda
              ! wavelength in metres
              wl = tab_lambda(lambda)*1.e-6
              delta_wl=tab_delta_lambda(lambda)*1.e-6
              cst_wl=cst/wl
              if (cst_wl < 200.0) then
                 coeff_exp=exp(cst_wl)
                 dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
              else
                 dB_dT = 0.0_dp
              endif
              total_sum = total_sum + dB_dT/(kappa(p_icell,lambda) * kappa_factor(icell))  * delta_wl
           enddo
           ! kappa_R = 4.*sigma * Temp**3 / (pi * total_sum)
           ! Dcoeff = c_light/(3kappa_R) car kappa volumique
           Dcoeff(i,j,k) =  cst_Dcoeff * total_sum/Temp**3
        endif
     enddo
  enddo

  ! Boundary condition
  Dcoeff(:,0,:) = Dcoeff(:,1,:)

  return

end subroutine setDiffusion_coeff

!************************************************************

subroutine setDiffusion_coeff0(i)
  ! Computes the diffusion coefficients for
  ! the cells in the zone where the
  ! diffusion approximation is used
  ! C. Pinte
  ! 15/02/07


  integer, intent(in) :: i

  real(kind=dp) :: cst_Dcoeff, wl, delta_wl, cst, cst_wl, coeff_exp, dB_dT, Temp, total_sum
  integer :: j, k, lambda, icell, p_icell

  cst_Dcoeff = pi/(12.*sigma)

  p_icell = icell1

  do k=1,n_az
     do j=j_start,nz
        if (j==0) cycle
        icell = cell_map(i,j,k)
        if (lvariable_dust) p_icell = icell
        Temp=Tdust(icell)
        cst=thermal_const/Temp
        total_sum=0.0_dp
        do lambda=1, n_lambda
           ! wavelength in metres
           wl = tab_lambda(lambda)*1.e-6
           delta_wl=tab_delta_lambda(lambda)*1.e-6
           cst_wl=cst/wl
           if (cst_wl < 200.0) then
              coeff_exp=exp(cst_wl)
              dB_dT = cst_wl*coeff_exp/((wl**5)*(coeff_exp-1.0)**2)
           else
              dB_dT = 0.0_dp
           endif
           total_sum = total_sum + dB_dT/(kappa(p_icell,lambda)*kappa_factor(icell)) * delta_wl
        enddo
        ! kappa_R = 4.*sigma * Temp**3 / (pi * total_sum)
        ! Dcoeff = c_light/(3kappa_R) car kappa volumique
        Dcoeff(i,j,k) =  cst_Dcoeff * total_sum/Temp**3
     enddo
  enddo

  ! Boundary condition
  ! Dcoeff(:,0,:) = Dcoeff(:,1,:)

  return

end subroutine setDiffusion_coeff0

!************************************************************

subroutine Temperature_to_DensE(ri)
  ! Computes the energy density of cells from their
  ! temperature
  ! Converts the whole grid to have cells on the boundary
  ! Executed only once, so no performance concern
  ! C. Pinte
  ! 15/02/07

  integer, intent(in) :: ri

  integer :: j, k, icell

  do j=j_start,nz
     if (j==0) cycle
     do k=1,n_az
        icell = cell_map(ri,j,k)
        DensE(ri,j,k) =  Tdust(icell)**4 ! The constant does not matter here (4.*sigma/c)
     enddo
  enddo

  ! Boundary condition: no flux at the mid-plane
  DensE(ri,0,:) = DensE(ri,1,:)

  DensE_m1(ri,:,:) = DensE(ri,:,:)

  return

end subroutine Temperature_to_DensE

!************************************************************

subroutine DensE_to_temperature()
  ! Computes the temperature of cells in the diffusion zone
  ! from their energy density
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
  ! Set the temperature to Tmin in the dark zone
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
  ! Compute the temperature in the deep disk zones using
  ! a diffusion equation (parabolic equation).
  ! The diffusion coefficient depends on T, making this a nonlinear problem.
  ! => solved with an explicit integration scheme: the steady-state solution
  ! is the limit as t --> inf of the transient solution
  ! C. Pinte
  ! 15/02/07

  real, dimension(n_cells) :: Temp0
  real :: max_delta_E_r, stabilite, precision
  integer :: n_iter, i
  logical :: lconverged


  write(*,*) "Computing 2D diffusion approx. in central parts of the disk"

  ! To initialise the first loop iteration
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

     ! Compute diffusion coefficients: boundary condition
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

        ! One time step of the diffusion equation
        call iter_Temp_approx_diffusion(stabilite,max_delta_E_r,lconverged)

        !test divergence
        if (.not.lconverged) then
           stabilite = 0.5 * stabilite
           precision = 0.1 * precision
           Tdust = Temp0
           Tdust_old = T_min
           exit infinie
        endif

        ! This is wrong: temperature is undefined here
!        Tdust_old = Tdust

!        write(*,*) n_iter, max_delta_E_r, precision  !, maxval(DensE)

        ! Test_convergence
        if (max_delta_E_r < precision) exit infinie

        ! Computes new diffusion coefficients for the next iteration
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
  ! Compute the temperature in the deep disk zones using
  ! a diffusion equation (parabolic equation).
  ! The diffusion coefficient depends on T, making this a nonlinear problem.
  ! => solved with an explicit integration scheme: the steady-state solution
  ! is the limit as t --> inf of the transient solution
  ! C. Pinte
  ! 15/02/07

  real :: max_delta_E_r, stabilite, precision
  integer :: n_iter, i, k
  logical :: lconverged


  write(*,*) "Computing 1+1D diffusion approx. in central parts of the disk"

  call clean_temperature()

  ! To initialise the first loop iteration
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

        ! Compute diffusion coefficients: boundary condition
        call setDiffusion_coeff0(i)

        if (stabilite < 0.01) call error("diffusion approximation does not seem to converge")

        ! Iterations
        n_iter = 0
        infinie : do
           n_iter = n_iter + 1

           ! One time step of the diffusion equation
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

           ! Computes new diffusion coefficients for the next iteration
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
  ! Performs one iteration of the explicit integration scheme
  ! of the diffusion equation on the cell energy density
  ! ie: one time step of the non-stationary equation
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
           dz = cell_height(i,j)
           ! tab_dt(i,j,k) = min(dr,dz)**2/Dcoeff(i,j,k)
           tab_dt(i,j,k) = 1.0_dp/(Dcoeff(i,j,k)*(1.0_dp/dr**2 + 1.0_dp/dz**2))
        enddo !j
     enddo !i
  enddo !k

  ! We take the minimum + a safety factor
  dt = stabilite * 0.5 * minval(tab_dt)

 ! write(*,*) "dt", dt

  ! Saves the energy density
  DensE_m1 = DensE

  ! Loop over the cells in the diffusion zone
  max_delta_E_r = 0.0

  ! TODO : Cross loop at each iteration to better smooth the noise ??


  do k=1, n_az
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,dE_dr_m1,dE_dr_p1,d2E_dr2,delta_E_r,D_Laplacien_E,delta_E,d2E_dz2) &
     !$omp shared(DensE_m1,r_grid,z_grid,Dcoeff,ri_in_dark_zone,ri_out_dark_zone,zj_sup_dark_zone,max_delta_E_r) &
     !$omp shared(DensE,dt,cell_height,k,n_rad,cell_map)
     !$omp do schedule(dynamic,10)
     do i=max(ri_in_dark_zone(k) -delta_cell_dark_zone,3), min(ri_out_dark_zone(k)+ delta_cell_dark_zone,n_rad-2)
        do j=1, zj_sup_dark_zone(i,k) + delta_cell_dark_zone

           ! Computes the Laplacian in cylindrical coordinates
           ! Warning: D enters the Laplacian since it depends on the position
           ! Only works in 2D for now
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

           !d2E_dz2  = (dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m) / (2.0_dp * cell_height(i)**2)

           d2E_dz2  =   Dcoeff(i,j,k) * (DensE_m1(i,j+1,k) + DensE_m1(i,j-1,k) - 2.0 * DensE(i,j,k)) / (2.0 * cell_height(i,j)**2)


           ! Laplacien
           D_Laplacien_E = d2E_dr2 + d2E_dz2

           ! We advance by one time step
           delta_E =  D_Laplacien_E * dt

           DensE(i,j,k) = DensE_m1(i,j,k) + delta_E

           if (DensE(i,j,k) < 0.) write(*,*) DensE(i,j,k), DensE_m1(i,j,k), delta_E

           ! Relative energy density increase
           delta_E_r = delta_E/DensE(i,j,k)
           if (delta_E_r > max_delta_E_r) max_delta_E_r = delta_E_r
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel
  enddo !k

  ! Boundary condition: no flux at the mid-plane
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
  ! Performs one iteration of the explicit integration scheme
  ! of the diffusion equation on the cell energy density
  ! ie: one time step of the non-stationary equation
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
     dz = cell_height(ri,j)
     tab_dt(j) = dz**2/Dcoeff(ri,j,k)
  enddo !j

  ! We take the minimum + a safety factor
  dt = stabilite * 0.5 * minval(tab_dt)


!  write(*,*) "dt", dt


  ! Saves the energy density
  DensE_m1(ri,:,:) = DensE(ri,:,:)

  ! Loop over the cells in the diffusion zone
  max_delta_E_r = 0.0

  do j=1, zj_sup_dark_zone(ri,k) + delta_cell_dark_zone

     dE_dz_p1 = (DensE_m1(ri,j+1,k) - DensE_m1(ri,j,k))
     dE_dz_m1 = (DensE_m1(ri,j,k) - DensE_m1(ri,j-1,k))


     ! It makes things unstable
    ! Dcoeff_p = 0.5_dp * (Dcoeff(ri,j,k) + Dcoeff(ri,j+1,k))
    ! Dcoeff_m = 0.5_dp * (Dcoeff(ri,j,k) + Dcoeff(ri-1,j-1,k))
     ! Plus stable
  !   Dcoeff_p = Dcoeff(ri,j,k)
  !   Dcoeff_m =  Dcoeff(ri,j,k)
  !   d2E_dz2  = (dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m) / (2.0_dp * cell_height(ri)**2)

     d2E_dz2  =  Dcoeff(ri,j,k) * (dE_dz_p1 - dE_dz_m1) / (2.0_dp * cell_height(ri,j)**2)


  !   write(*,*) "****"
  !   write(*,*) dE_dz_p1, dE_dz_m1
  !   write(*,*) Dcoeff_p, Dcoeff_m
  !   write(*,*) j, dE_dz_p1*Dcoeff_p - dE_dz_m1*Dcoeff_m, dE_dz_p1*Dcoeff_p

     ! Laplacien
     D_Laplacien_E =  d2E_dz2

     ! We advance by one time step
     delta_E =  D_Laplacien_E * dt

   !  write(*,*) "DeltaE", j, delta_E

     DensE(ri,j,k) = DensE_m1(ri,j,k) + delta_E

     ! Relative energy density increase
     delta_E_r = delta_E/DensE(ri,j,k)
     if (delta_E_r > max_delta_E_r) max_delta_E_r = delta_E_r
  enddo !j

  ! Boundary condition: no flux at the mid-plane
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

subroutine compute_Planck_opacities(icell, Planck_opacity,rec_Planck_opacity)
  ! Compute current Planck reciprocal mean opacity for all cells
  ! (note : diffusion coefficient needs to be defined with Rosseland opacity
  ! in B&W mode)
  ! Diffusion coefficient is D = 1/(rho * opacity)
  ! This opacity/diffusion coefficient includes scattering
  ! See Min et al 2009 and Robitaille et al 2010
  use parameters
  use constants
  use wavelengths, only : n_lambda, tab_lambda, tab_delta_lambda
  use Temperature, only : Tdust
  use dust_prop, only : kappa

  integer,  intent(in)  :: icell
  real(dp), intent(out) :: Planck_opacity,rec_Planck_opacity ! cm2/g (ie per gram of gas)

  integer :: lambda
  real(dp) :: total_sum, somme2, cst, cst_wl, B, coeff_exp, wl, delta_wl, norm, T !dB_dT

  integer, pointer :: p_icell
  integer, target :: icell0

  icell0 = icell

  T = Tdust(icell) ! WARNING : this is 0 until Temp_finale
  T = 20
!  if ((T > 1) .and. (Voronoi(icell)%original_id > 0)) then
  if (T > 1) then
     if (lvariable_dust) then
        p_icell => icell0
     else
        p_icell => icell1
     endif

     total_sum  = 0.0_dp
     somme2 = 0.0_dp
     norm = 0.0_dp
     cst    = thermal_const/T
     do lambda = 1,n_lambda
        ! wavelength in metres
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
        total_sum  = total_sum  + B/kappa(p_icell,lambda)
        somme2  = somme2  + B * kappa(p_icell,lambda)
        norm = norm + B*delta_wl
     enddo
     rec_Planck_opacity = norm/total_sum * kappa_factor(icell)
     Planck_opacity = somme2/norm * kappa_factor(icell)
  else
     rec_Planck_opacity = 1e-30
     Planck_opacity = 1e-30
  endif

end subroutine compute_Planck_opacities

!************************************************************

end module diffusion
