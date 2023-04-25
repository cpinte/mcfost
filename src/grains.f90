module grains

  use mcfost_env, only : dp

  implicit none
  save

  type dust_grain_type
     integer :: zone, methode_chauffage, pop
     logical :: is_PAH
  end type dust_grain_type

  type(dust_grain_type), dimension(:), allocatable  :: grain

  type dust_pop_type
     integer :: n_grains, methode_chauffage, zone
     real :: amin, amax, aexp, frac_mass, porosity ! aexp is > 0

     real :: avg_grain_mass, masse, sblow, rho1g_avg, T_sub, dhs_maxf

     character(len=3) :: type ! Mie or DHS
     integer ::  mixing_rule, n_components ! en cas de coating, composant 1 est le coeur 2 est le manteau
     ! mixing rule : 1 = EMT, 2 = coating
     character(len=512), dimension(10) :: indices
     real, dimension(10) :: component_rho1g, component_volume_fraction, component_T_sub

     logical :: is_opacity_file, is_PAH, is_Misselt_opacity_file, is_DustEM_opacity_file, lcoating
     integer :: ind_debut, ind_fin
  end type dust_pop_type

  type(dust_pop_type), dimension(:), allocatable, target :: dust_pop

  ! Nombre de taille de grains
  ! le pas et (amax - amin)/n_grains
  integer :: n_grains_tot , n_grains_RE_LTE, n_grains_RE_nLTE, n_grains_nRE
  integer :: grain_RE_LTE_start, grain_RE_LTE_end, grain_RE_nLTE_start, grain_RE_nLTE_end, grain_nRE_start, grain_nRE_end

  real(kind=dp),dimension(:), allocatable :: nbre_grains !n_grains_tot
  real, dimension(:), allocatable :: r_grain, r_grain_min, r_grain_max, r_core, S_grain, M_grain !n_grains_tot
  real, dimension(:,:), allocatable :: frac_mass_pop !n_zones, n_pop

  ! Pour lecture des fichiers d'opacite, par exemple PAH de B.Draine
  integer, dimension(:), allocatable :: op_file_na,  op_file_n_lambda ! n_pop
  real, dimension(:,:,:), allocatable :: op_file_Qext, op_file_Qsca, op_file_g ! op_file_n_lambda,op_file_na, n_pop
  real, dimension(:,:), allocatable :: op_file_log_r_grain ! op_file_na, n_pop
  real, dimension(:,:), allocatable :: op_file_lambda, op_file_delta_lambda ! op_file_n_lambda, n_pop
  integer, dimension(:), allocatable :: file_sh_nT
  real(kind=dp), dimension(:,:), allocatable :: file_sh_T, file_sh ! nT, n_pop

  real, dimension(:,:), allocatable :: amin, amax, aexp ! n_zones, n_especes

  ! Parametres de diffusion des grains
  real, dimension(:,:,:), allocatable :: tab_s11, tab_s12, tab_s33, tab_s34, prob_s11 !n_lambda,n_grains,180
  real, dimension(:,:), allocatable :: tab_g, tab_albedo, C_ext, C_sca, C_abs, C_abs_norm !n_grains, n_lambda
  !real, dimension(:), allocatable :: q_geo ! n_grains section geometrique en m^2

  logical :: lforce_HG
  real :: forced_g

  ! aggregats
  real, dimension(:,:,:,:,:), allocatable :: tab_mueller !4,4, 180, n_grains,n_lambda

  ! Parametres de diffusion des cellules
  real, dimension(:,:), allocatable :: tab_albedo_pos, tab_g_pos ! n_cells,n_lambda
  real, dimension(:,:,:), allocatable :: tab_s11_pos, tab_s12_o_s11_pos, tab_s33_o_s11_pos, tab_s34_o_s11_pos, prob_s11_pos ! 0:180, n_cells,n_lambda
  real, dimension(:,:,:,:,:), allocatable :: tab_mueller_pos ! 4,4,0:180, n_cells,n_lambda

  character(len=512) :: aggregate_file, mueller_aggregate_file, mueller_file
  real :: R_sph_same_M

end module grains
