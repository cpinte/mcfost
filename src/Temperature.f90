module Temperature

  use mcfost_env
  use parametres

  implicit none
  save

  real, dimension(:), allocatable :: tab_Temp

  real, dimension(:), allocatable :: Tdust, Tdust_old !n_rad,nz,n_az
  real, dimension(:,:), allocatable :: Tdust_1grain, Tdust_1grain_nRE !n_rad,nz, n_grains
  real, dimension(:,:), allocatable :: Tdust_1grain_old, Tdust_1grain_nRE_old, maxP_old !n_rad,nz, n_grains
  integer, dimension(:,:), allocatable :: Tpeak_old
  real, dimension(:,:,:), allocatable :: Proba_Tdust !n_T, n_cells, n_grains

  real, dimension(:), allocatable :: E_disk !n_lambda

  logical, dimension(:,:), allocatable :: l_RE, lchange_nRE ! n_grains, n_cells

contains

  subroutine init_tab_Temp()

    real(kind=dp) :: delta_T
    integer :: t

    tab_Temp=0.0
    ! Echantillonage temperature
    !delta_T=(T_max)**(1.0/(n_T-1))
    delta_T=exp((1.0_dp/(real(n_T,kind=dp)))*log(T_max/T_min))
    tab_Temp(1)=T_min*sqrt(delta_T)
    do t=2,n_T
       tab_Temp(t)=delta_T*tab_Temp(t-1)
    enddo

    return

  end subroutine init_tab_Temp

end module Temperature
