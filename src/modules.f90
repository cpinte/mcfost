module grains
!****************************
! Parametres des grains
!****************************
  use mcfost_env, only : dp

  implicit none
  save

  type dust_grain_type
     real :: nbre_grains !, r
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

  real,dimension(:), allocatable :: nbre_grains !n_grains_tot
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

  ! Tab de lambda

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

  character(len=512) :: aggregate_file, mueller_aggregate_file
  real :: R_sph_same_M

end module grains

!********************************************************************

module opacity

  use mcfost_env, only : dp

  implicit none
  save

  integer, dimension(:), allocatable :: lexit_cell
  real(kind=dp), dimension(:), allocatable :: zmax !n_rad
  real(kind=dp), dimension(:), allocatable :: volume !n_rad en AU^3
  real(kind=dp), dimension(:), allocatable :: masse  !en g ! n_cells
  real(kind=dp), dimension(:,:), allocatable :: masse_rayon ! en g!!!!  n_rad, n_az
  real(kind=dp), dimension(:), allocatable :: delta_z ! taille verticale des cellules cylindriques
  real(kind=dp), dimension(:), allocatable :: dr2_grid ! differentiel en r^2 des cellules
  real(kind=dp), dimension(:), allocatable :: r_grid, z_grid ! Position en cylindrique !!! des cellules
  real(kind=dp), dimension(:), allocatable :: phi_grid
  real(kind=dp), dimension(:), allocatable :: r_lim, r_lim_2, r_lim_3 ! lim rad sup de la cellule (**2) !0:n_rad
  real(kind=dp), dimension(:,:), allocatable :: z_lim ! lim vert inf de la cellule !n_rad,nz+1
  real(kind=dp), dimension(:), allocatable :: tan_phi_lim ! lim azimuthale de la cellule ! n_az
  real(kind=dp), dimension(:), allocatable :: w_lim, theta_lim, tan_theta_lim ! lim theta sup de la cellule ! 0:nz
  integer, dimension(:), allocatable :: tab_region ! n_rad : indice de region pour chaque cellule

  integer, dimension(:,:,:), allocatable :: cell_map
  integer, dimension(:), allocatable :: cell_map_i, cell_map_j, cell_map_k

  real, dimension(:,:), allocatable :: kappa !n_cells, n_lambda
  real, dimension(:,:), allocatable :: kappa_abs_LTE ! n_cells, n_lambda
  real, dimension(:,:), allocatable :: kappa_abs_nLTE, kappa_abs_RE ! n_cells, n_lambda
  real, dimension(:,:), allocatable :: proba_abs_RE, proba_abs_RE_LTE, proba_abs_RE_LTE_p_nLTE
  real, dimension(:,:,:), allocatable :: kabs_nLTE_CDF, kabs_nRE_CDF ! 0:n_grains, n_cells, n_lambda
  real(kind=dp), dimension(:,:), allocatable :: emissivite_dust ! emissivite en SI (pour mol)

  real, dimension(:,:), allocatable :: densite_pouss ! n_grains, n_cells en part.cm-3
  integer :: icell_not_empty

  real, dimension(:,:,:), allocatable :: ksca_CDF ! 0:n_grains, n_cells, n_lambda
  !* ksca_CDF(i) represente la probabilite cumulee en-dessous d'une
  !* certaine taille de grain. Ce tableau est utilise pour le tirage
  !* aleatoire de la taille du grain diffuseur, puisqu'elle doit prendre
  !* en compte le nombre de grains en meme temps que leur probabilite
  !* individuelle de diffuser (donnee par qsca*pi*a**2).

  logical :: l_is_dark_zone
  logical, dimension(:), allocatable :: l_dark_zone !n_cells
  integer, parameter :: delta_cell_dark_zone=3

  integer, dimension(:), allocatable :: ri_in_dark_zone, ri_out_dark_zone !n_az
  integer, dimension(:,:), allocatable :: zj_sup_dark_zone, zj_inf_dark_zone !n_rad, n_az

  logical, dimension(:,:), allocatable :: l_emission_pah !0:n_rad+1,0:nz+1

  real(kind=dp), dimension(:,:,:), allocatable :: Dcoeff !n_rad, n_z, n_az
  real, dimension(:,:,:), allocatable :: DensE, DensE_m1 !n_rad, 0:n_z, n_az

end module opacity


!********************************************************************

module molecular_emission

  use mcfost_env, only : dp

  implicit none
  save

  logical :: ldouble_RT

  ! Lecture du fichier de data mol
  real :: molecularweight
  real :: abundance
  integer :: nLevels
  real(kind=dp) :: largeur_profile
  real, dimension(:), allocatable :: Level_energy
  integer, dimension(:), allocatable ::  j_qnb
  ! g est dp car les calculs utilisant g sont en dp
  real(kind=dp), dimension(:), allocatable :: poids_stat_g
  integer :: nTrans_tot

  real(kind=dp), dimension(:), allocatable :: Aul, Blu, Bul, fAul, fBlu, fBul, transfreq
  integer, dimension(:), allocatable :: itransUpper, itransLower
  integer :: nCollPart
  character(len=512), dimension(:), allocatable :: collBetween
  integer, dimension(:), allocatable :: nCollTrans, nCollTemps
  real, dimension(:,:), allocatable :: collTemps

  integer, dimension(:,:), pointer :: iCollUpper, iCollLower
  real, dimension(:,:,:), pointer :: collRates

  real, dimension(:), allocatable :: Tcin ! Temperature cinetique
  real :: correct_Tgas
  logical :: lcorrect_Tgas

  real :: nH2, masse_mol
  ! masse_mol_gaz sert uniquement pour convertir masse disque en desnite de particule
  real(kind=dp), dimension(:,:), allocatable :: kappa_mol_o_freq, kappa_mol_o_freq2 ! n_cells, nTrans
  real(kind=dp), dimension(:,:), allocatable :: emissivite_mol_o_freq,  emissivite_mol_o_freq2 ! n_cells, nTrans
  real, dimension(:), allocatable :: vfield, vfield_x, vfield_y, vfield_z ! n_cells
  real, dimension(:,:), allocatable :: tab_nLevel, tab_nLevel2, tab_nLevel_old ! n_cells, nLevels

  real, dimension(:), allocatable :: v_turb, v_line ! n_cells

  real ::  vitesse_turb, dv, dnu
  integer, parameter :: n_largeur_Doppler = 15
  real(kind=dp), dimension(:), allocatable :: tab_v ! n_speed

  ! densite_gaz gives the midplane density for j=0
  real(kind=dp), dimension(:), allocatable :: densite_gaz, masse_gaz ! n_rad, nz, n_az, Unites: part.m-3 et g : H2
  real(kind=dp), dimension(:), allocatable :: densite_gaz_midplane
  real(kind=dp), dimension(:), allocatable :: Surface_density

  real(kind=dp), dimension(:,:), allocatable :: ds
  real(kind=dp), dimension(:,:,:,:), allocatable :: I0, I02 ! nSpeed,nTrans,iray,ncpus
  real(kind=dp), dimension(:,:,:), allocatable :: I0c ! Intensite dans le continu: nTrans,iray,ncpus
  real(kind=dp), dimension(:,:,:), allocatable :: Doppler_P_x_freq

  real(kind=dp), dimension(:,:), allocatable :: Jmol, Jmol2 ! nTrans, n_cpu
  real(kind=dp), dimension(:), allocatable :: tab_Cmb_mol ! nTrans

  logical :: linfall, lkeplerian, lcylindrical_rotation
  real :: chi_infall

  real(kind=dp), dimension(:), allocatable :: deltaVmax ! n_cells
  real(kind=dp), dimension(:), allocatable :: tab_dnu_o_freq ! n_cells
  real(kind=dp), dimension(:), allocatable :: norme_phiProf_m1, sigma2_phiProf_m1 ! n_cells

  real, dimension(:), allocatable :: tab_abundance ! n_cells
  logical, dimension(:), allocatable :: lcompute_molRT ! n_cells

  logical ::  lfreeze_out, lphoto_dissociation, lphoto_desorption
  real :: T_freeze_out, freeze_out_depletion

  real(kind=dp), dimension(:,:,:,:), allocatable ::  origine_mol ! nv, nTrans, n_cells, nb_proc

  integer :: RT_line_method, n_molecules

  type molecule
     integer :: n_speed_rt, n_speed_center_rt, n_extraV_rt, nTrans_raytracing, iLevel_max
     real :: vmax_center_rt, extra_deltaV_rt, abundance
     logical :: lcst_abundance, lline
     character(len=512) :: abundance_file, filename
     character(len=32) :: name
     integer, dimension(100) :: indice_Trans_rayTracing
  end type molecule

  type(molecule), dimension(:), allocatable :: mol

  real(kind=dp), dimension(:), allocatable :: tab_speed_rt

  real, dimension(:,:), allocatable :: maser_map ! n_cells, n_trans


end module molecular_emission

!***********************************************************

module ray_tracing

  use mcfost_env, only : dp
  use parametres, only : nang_scatt

  implicit none
  save


  real(kind=dp), dimension(:,:), allocatable :: n_phot_envoyes


  ! inclinaisons
  real :: RT_imin, RT_imax, RT_az_min, RT_az_max
  integer ::  RT_n_incl, RT_n_az
  logical :: lRT_i_centered

  real, dimension(:), allocatable :: tab_RT_incl, tab_RT_az
  real(kind=dp), dimension(:), allocatable :: tab_uv_rt, tab_w_rt
  real(kind=dp), dimension(:,:), allocatable :: tab_u_rt, tab_v_rt

  ! Sauvegarde champ de radiation pour rt2
  integer ::  n_phi_I,  n_theta_I ! 15 et 9 ok avec 30 et 30 en mode SED

  ! Pour rt 2 : nbre d'angle de visee en azimuth
  ! TODO : calculer automatiquement en fct de la fct de phase + interpolation
  integer :: nang_ray_tracing, nang_ray_tracing_star

  real, dimension(:,:,:), allocatable :: tab_s11_ray_tracing, tab_s12_ray_tracing, tab_s33_ray_tracing, tab_s34_ray_tracing ! 0:nang_scatt, n_cells, n_lambda
  real, dimension(:,:,:), allocatable ::  cos_thet_ray_tracing, omega_ray_tracing ! nang_ray_tracing, 2 (+z et -z), nb_proc
  real, dimension(:,:,:), allocatable ::  cos_thet_ray_tracing_star, omega_ray_tracing_star ! nang_ray_tracing, 2 (+z et -z), nb_proc

  real, dimension(0:nang_scatt) :: tab_cos_scatt

  ! intensite specifique
  real, dimension(:), allocatable :: J_th ! n_cells

  ! methode RT 1 : saving scattered specific intensity (SED + image 3D)
  ! todo faire sauter le 2 pour gagner une dimension et rester sous la limite de 7
  integer :: n_az_rt, n_theta_rt
  real, dimension(:,:,:,:,:,:), allocatable ::  xI_scatt ! 4, RT_n_incl * RT_n_az, n_cells, n_az_rt, n_theta_rt, ncpus
  real(kind=dp), dimension(:,:,:), allocatable ::  I_scatt ! 4, n_az_rt, 2
  integer, dimension(:,:,:), allocatable :: itheta_rt1 ! RT_n_incl,RT_n_az,nb_proc
  real(kind=dp), dimension(:,:,:), allocatable ::  sin_omega_rt1, cos_omega_rt1, sin_scatt_rt1 ! RT_n_incl,RT_n_az,nb_proc
  real(kind=dp), dimension(:,:,:,:), allocatable ::  eps_dust1 !N_type_flux, n_cells, n_az_rt,n_theta_rt

  ! methode RT 2 : saving specific intensity (image 2D)
  real, dimension(:,:,:,:,:), allocatable :: I_spec ! 4, n_theta_I, n_phi_I, n_cells, ncpus
  real, dimension(:,:), allocatable :: I_spec_star ! n_cells, ncpus

  ! Fonction source: Ok en simple
  real, dimension(:,:,:,:), allocatable ::  I_sca2 ! n_type_flux, nang_ray_tracing, 2, n_cells
  real, dimension(:,:,:,:), allocatable ::  eps_dust2 ! n_type_flux, nang_ray_tracing, 2, n_rad, nz
  real, dimension(:,:,:,:), allocatable ::  eps_dust2_star ! n_type_flux, nang_ray_tracing, 2, n_rad, nz

  real, dimension(:,:,:,:,:,:,:), allocatable :: Stokes_ray_tracing ! n_lambda, nx, ny, RT_n_incl, RT_n_az, n_type_flux, ncpus
  real, dimension(:,:,:,:,:,:), allocatable :: tau_surface ! nx, ny, RT_n_incl, RT_n_az, 3, ncpus
  real, dimension(:,:,:), allocatable :: stars_map ! nx, ny, 4

  real, dimension(:,:,:,:,:), allocatable :: weight_Inu_fct_phase ! n_rayon_rt, dir, n_theta_I, n_phi_I, nang_scatt

end module ray_tracing

!********************************************************************
