module radiation_field

  use parametres
  use constantes
  use ray_tracing, only : lscatt_ray_tracing1, lscatt_ray_tracing2, n_theta_I, n_phi_I, n_az_rt, I_spec, I_spec_star
  use opacity, only : kappa_abs_LTE
  use messages

  implicit none

  public :: J0, xJ_abs, xN_abs, xKJ_abs, E0
  public :: save_radiation_field, allocate_radiation_field_step1, allocate_radiation_field_step2, &
       deallocate_radiation_field, reset_radiation_field

  private

  ! pour stockage des cellules par lequelles on passe
  ! longueur de vol cumulee dans la cellule
  real(kind=dp), dimension(:,:), allocatable :: xKJ_abs ! n_cells, id
  ! xJabs represente J_lambda dans le bin lambda -> bin en log : xJabs varie comme lambda.F_lambda
  real(kind=dp), dimension(:,:,:), allocatable :: xJ_abs ! n_cells, n_lambda, nb_proc
  real, dimension(:,:,:), allocatable :: xN_abs ! n_cells, n_lambda, nb_proc, for ProDiMo

  ! Initial energy and radiation field
  real(kind=dp), dimension(:), allocatable :: E0 ! n_cells
  real(kind=dp), dimension(:,:), allocatable :: J0 !n_cells, n_lambda, n_rad, nz, n_az


contains

subroutine save_radiation_field(id,lambda,p_lambda,icell0, Stokes, l,  x0,y0,z0, x1,y1,z1, u,v, w, flag_star, flag_direct_star)

  use dust_ray_tracing, only : calc_xI_scatt, calc_xi_scatt_pola

  integer, intent(in) :: id,lambda,p_lambda,icell0
  real(kind=dp), dimension(4), intent(in) :: Stokes
  real(kind=dp) :: l, x0,y0,z0, x1,y1,z1, u,v,w
  logical, intent(in) :: flag_star, flag_direct_star


  real(kind=dp) :: xm,ym,zm, phi_pos, phi_vol
  integer :: psup, phi_I, theta_I, phi_k


  if (letape_th) then
     if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + kappa_abs_LTE(icell0,lambda) * l * Stokes(1)
     if (lxJ_abs_step1) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
  else
     if (lxJ_abs) then ! loutput_UV_field .or. loutput_J .or. lprodimo
        xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
        ! Pour statistique: nbre de paquet contribuant a intensite specifique
        if (lProDiMo) xN_abs(icell0,lambda,id) = xN_abs(icell0,lambda,id) + 1.0
     endif ! lProDiMo

     if (lscatt_ray_tracing1) then
        xm = 0.5_dp * (x0 + x1)
        ym = 0.5_dp * (y0 + y1)
        zm = 0.5_dp * (z0 + z1)

        if (l3D) then ! phik & psup=1 in 3D
           phi_k = 1
           psup = 1
        else
           phi_pos = atan2(xm,ym)
           phi_k = floor(  modulo(phi_pos, deux_pi) / deux_pi * n_az_rt ) + 1
           if (phi_k > n_az_rt) phi_k=n_az_rt

           if (zm > 0.0_dp) then
              psup = 1
           else
              psup = 2
           endif
        endif

        if (lsepar_pola) then
           call calc_xI_scatt_pola(id,lambda,p_lambda,icell0,phi_k,psup,l,Stokes(:),flag_star)
        else
           ! ralentit d'un facteur 5 le calcul de SED
           ! facteur limitant
           call calc_xI_scatt(id,lambda,p_lambda,icell0,phi_k,psup,l,Stokes(1),flag_star)
        endif

     else if (lscatt_ray_tracing2) then ! only 2D
        if (flag_direct_star) then
           I_spec_star(icell0,id) = I_spec_star(icell0,id) + l * Stokes(1)
        else
           xm = 0.5_dp * (x0 + x1)
           ym = 0.5_dp * (y0 + y1)
           zm = 0.5_dp * (z0 + z1)
           phi_pos = atan2(xm,ym)

           phi_vol = atan2(-u,-v) + deux_pi ! deux_pi pour assurer diff avec phi_pos > 0

           !  if (l_sym_ima) then
           !     delta_phi = modulo(phi_vol - phi_pos, deux_pi)
           !     if (delta_phi > pi) delta_phi = deux_pi - delta_phi
           !     phi_I =  nint( delta_phi  / pi * (n_phi_I -1) ) + 1
           !     if (phi_I > n_phi_I) phi_I = n_phi_I
           !  else
           phi_I =  floor( modulo(phi_vol - phi_pos, deux_pi) / deux_pi * n_phi_I ) + 1
           if (phi_I > n_phi_I) phi_I = 1
           !  endif

           if (zm > 0.0_dp) then
              theta_I = floor(0.5_dp*( w + 1.0_dp) * n_theta_I) + 1
           else
              theta_I = floor(0.5_dp*(-w + 1.0_dp) * n_theta_I) + 1
           endif
           if (theta_I > n_theta_I) theta_I = n_theta_I

           I_spec(1:n_Stokes,theta_I,phi_I,icell0,id) = I_spec(1:n_Stokes,theta_I,phi_I,icell0,id) + l * Stokes(1:n_Stokes)

           if (lsepar_contrib) then
              if (flag_star) then
                 I_spec(n_Stokes+2,theta_I,phi_I,icell0,id) = I_spec(n_Stokes+2,theta_I,phi_I,icell0,id) + l * Stokes(1)
              else
                 I_spec(n_Stokes+4,theta_I,phi_I,icell0,id) = I_spec(n_Stokes+4,theta_I,phi_I,icell0,id) + l * Stokes(1)
              endif
           endif ! lsepar_contrib

        endif ! flag_direct_star
     endif !lscatt_ray_tracing
  endif !letape_th

  return

end subroutine save_radiation_field

!*************************************************************************************

subroutine allocate_radiation_field_step1(Nc)

  integer, intent(in) :: Nc
  integer :: alloc_status

  if (lRE_LTE) then
     allocate(xKJ_abs(Nc,nb_proc), E0(Nc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error xKJ_abs')
     xKJ_abs = 0.0 ; E0 = 0.0
  endif

  lxJ_abs = lProDiMo.or.lML.or.loutput_UV_field.or.loutput_J
  lxJ_abs_step1 = lRE_nLTE .or. lnRE .or. loutput_J_step1
  if (lxJ_abs_step1 .or. (lxJ_abs.and.lsed.and.lsed_complete)) then
     allocate(xJ_abs(Nc,n_lambda,nb_proc), J0(Nc,n_lambda), stat=alloc_status) ! BIG array
     if (alloc_status > 0) call error('Allocation error xJ_abs')
     xJ_abs=0.0 ; J0 = 0.0
  endif

  return

end subroutine allocate_radiation_field_step1

!*************************************************************************************

subroutine allocate_radiation_field_step2()

  integer :: alloc_status

  lxJ_abs = lProDiMo.or.loutput_UV_field.or.loutput_J
  if (lxJ_abs) then
     allocate(xJ_abs(n_cells,n_lambda,nb_proc), J0(n_cells,n_lambda), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error xJ_abs in realloc_step2')
     xJ_abs = 0.0 ; J0 = 0.0
  endif


  return

end subroutine allocate_radiation_field_step2

!*************************************************************************************

subroutine deallocate_radiation_field()

  if (allocated(xKJ_abs)) deallocate(xKJ_abs,E0)
  if (allocated(xJ_abs)) deallocate(xJ_abs,J0)

  return

end subroutine deallocate_radiation_field

!*************************************************************************************

subroutine reset_radiation_field()

  if (lRE_LTE) then
     xKJ_abs(:,:) = 0.0_dp
     E0 = 0.0_dp
  endif
  if (lRE_nLTE .or. lnRE) then
     xJ_abs(:,:,:) = 0.0_dp
     J0 = 0.0_dp
  endif

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
!  log_E_em = 0.0
!
!
!  kdB_dT_CDF = 0

  return

end subroutine reset_radiation_field

!**********************************************************************

end module radiation_field
