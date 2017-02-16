module parametres
!****************************
! Parametres generaux
!****************************

  implicit none
  save

  real, parameter :: mcfost_version = 3.0
  character(8), parameter :: mcfost_release = "3.0.10"
  real, parameter :: required_utils_version = 3.0

  character(len=128), parameter :: webpage=      "http://ipag.osug.fr/public/pintec/mcfost/"
  character(len=128), parameter :: utils_webpage="http://ipag.osug.fr/public/pintec/mcfost_utils/"

  real :: para_version

  ! Système
  integer :: nb_proc
  integer, parameter :: cache_line_size = 64 ! 64 bytes = 16 floats = 8 double, from Core 2 Duo to i7 + Xeon Phi
  logical :: lpara, lstop_after_init
  integer, parameter :: sp = selected_real_kind(p=6,r=37)
  integer, parameter :: dp = selected_real_kind(p=13,r=200)
  integer, parameter :: limite_stack = 5000000
  integer :: indice_etape, etape_i, etape_f
  integer :: time_begin, time_end, time_tick, time_max

  real :: max_mem = 4. ! GBytes maximum size for 1 array (mcfost can have 2 arrays of this size)
  logical :: low_mem_scattering, low_mem_th_emission, low_mem_th_emission_nLTE, lMueller_pos_multi

  ! Nombre de photons lances
  logical :: ldust_transfer
  integer :: nbre_photons_loop, nbre_photons_eq_th, nbre_photons_lambda, nbre_photons_image, nbre_photons_spectre
  real :: nbre_photons_lim = 1.e4 ! combien de fois plus on aurait recu sans disque
  integer :: nnfot1
  real(kind=dp) :: E_paquet
  integer :: n_dif_max_eq_th = 100000 ! Nbre max de dif autorises dans calcul eq. th OUTDATED
  real :: tau_dark_zone_eq_th = 1500 !1500.   15000 pour benchmark a tau=1e6
  real :: tau_dark_zone_obs = 100 ! idem que 1000 si 1500. ci-dessus
  integer :: n_Stokes

  ! Nbre d'angles pour echantillonner la fct de phase
  integer, parameter :: nang_scatt = 180  ! TODO : ca bug si c'est pas 180

  ! Nbre de longueurs d'onde utilisees
  integer :: n_lambda
  logical :: lmono0, lmono

  ! lvariable_dust = true si les proprites de la poussiere sont variables dans chaque cellule
  logical :: lvariable_dust, lmigration, lhydrostatic, ldust_sublimation
  integer :: settling_type ! 1 = Parametric, 2 = Dubrulle or 3 = Fromang

  logical :: lRE_LTE, lRE_nLTE, lnRE, lonly_LTE, lonly_nLTE, loutput_J, loutput_UV_field, lxJ_abs

  ! Methode de calcul de la diffusion : a choisir pour optimiser taille memoire et temps cpu
  ! 0 -> automatique
  ! 1 -> choix taille du grain diffuseur + matrice Mueller par grain
  ! 2 -> matrice de Mueller moyenne par cellule (benchmark)
  integer :: scattering_method
  logical :: lscattering_method1

  ! Theorie de Mie ou HG
  integer :: aniso_method ! 1 = full phase function, 2 = HG
  logical :: lmethod_aniso1

  integer :: RT_sed_method ! cf routine dust_map pour def

  ! Etapes de l'émission thermique
  logical :: ltemp, lsed, lsed_complete, l_em_disk_image, lchauff_int
  character(len=512), dimension(:), allocatable :: indices
  character(len=512) :: tab_wavelength

  ! Emission moleculaire
  logical :: lemission_mol,  lpop, lprecise_pop, lmol_LTE, ldust_mol, lonly_top, lonly_bottom

  ! Decomposition image
  logical :: lsepar_contrib, lsepar_pola, lonly_capt_interet
  integer :: N_type_flux
  ! les flux sont I, (Q,U,V), (star, scatt, disk th, disk th scatt.)

  ! rotation du plan du disque en deg., sens trigo.
  real(kind=dp) :: ang_disque
  real(kind=dp) :: sin_disk, cos_disk, cos_disk_x2, sin_disk_x2

  ! Production d'images symetriques
  ! La symetrie est effectuee avant de choisir les pixels
  ! le système est-il centrosymetrique
  ! le systeme a-t-il une symetrie axiale (ne compte que si N_phi > 1)
  logical :: l_sym_ima, l_sym_centrale, l_sym_axiale


  ! Parametres des cartes
  integer :: N_thet, N_incl, N_phi, capt_interet, delta_capt, capt_inf, capt_sup, capt_debut, capt_fin
  integer ::  igridx, igridy, deltapix_x, deltapix_y, maxigrid, npix_x, npix_y
  real :: angle_interet, zoom, size_pix, tau_seuil, wl_seuil

  real  :: cutoff = 7.0

  ! Résolution de la grille de densité
  ! Nombre de cellules dans la direction r (echantillonage log)
  integer :: grid_type ! 1 = cylindrical, 2 = spherical
  integer :: n_rad, n_rad_in  ! subdivision de la premiere cellule
  ! Nombre de couches verticales ( = de stratifications)
  integer :: nz, p_n_rad, p_nz, p_n_az, p_n_lambda_pos, p_n_lambda_grain
  ! Nombre de cellules azimuthales
  integer :: n_az, j_start, pj_start
  ! Nombre de cellules totale
  integer :: n_cells, nrz, p_n_cells, icell_ref

  integer :: n_lambda2

  logical :: letape_th, limg, lorigine, laggregate, l3D, lremove, lwarp, lcavity, ltilt, lwall
  logical :: lopacite_only, lseed, ldust_prop, ldisk_struct, loptical_depth_map, lreemission_stats
  logical :: lapprox_diffusion, lcylindrical, lspherical, lVoronoi, is_there_disk, lno_backup, lonly_diff_approx, lforce_diff_approx
  logical :: laverage_grain_size, lisotropic, lno_scattering, lqsca_equal_qabs
  logical :: ldensity_file, lsigma_file, lphantom_file, lgadget2_file, lascii_SPH_file, llimits_file
  logical :: lweight_emission, lcorrect_density, lProDiMo2mcfost, lProDiMo2mcfost_test, lastrochem
  logical :: lspot, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium, lchange_Tmax_PAH, lISM_heating, lcasa

  character(len=512) :: mcfost_utils, my_mcfost_utils, data_dir, root_dir, basename_data_dir, seed_dir
  character(len=512) :: lambda_filename, band, model_pah, pah_grain, cmd_opt
  character(len=512), dimension(100) :: data_dir2, basename_data_dir2
  character(len=512), dimension(:), allocatable :: search_dir, dust_dir, mol_dir, star_dir, lambda_dir

  ! benchmarks
  logical :: lbenchmark_Pascucci, lbenchmark_vanZadelhoff1, lbenchmark_vanZadelhoff2, lDutrey94, lHH30mol
  logical :: lbenchmark_water1, lbenchmark_water2, lbenchmark_water3, lbenchmark_SHG, lMathis_field
  real :: Mathis_field

  ! Prodimo
  logical :: lprodimo, lprodimo_input_dir, lforce_ProDiMo_PAH

  logical, parameter :: ltest_rt3 = .false. ! marche pas
  logical, parameter :: ltest_rt4 = .false.  ! marche pas non plus

  logical :: lSeb_Charnoz, lread_Seb_Charnoz, lread_Seb_Charnoz2, lread_Misselt, lread_DustEM, lread_grain_size_distrib

end module parametres

!********************************************************************

module disk
!****************************
! Parametres du disque
!****************************

  use parametres

  implicit none
  save

  real :: distance ! Distance du disque en pc
  real(kind=dp) :: map_size

  ! Disque decrit en une ou n zones
  integer :: n_zones, n_regions

  type disk_zone_type
     real(kind=dp) :: Rin, Rmin, Rc, Rout, Rmax, Rref, edge, exp_beta, surf
     real(kind=dp) :: moins_gamma_exp, sclht, diskmass, gas_to_dust, vert_exponent
     integer :: geometry ! 1=disk, 2=tappered-disk, 3=envelope
     integer :: region
  end type disk_zone_type

  type disk_region_type
     integer :: n_zones
     real(kind=dp) :: Rmin, Rmax
     integer :: iRmin, iRmax
     integer, dimension(:), allocatable :: zones
  end type disk_region_type

  type cavity_type
     real(kind=dp) ::  exp_beta, sclht, Rref
  end type cavity_type

  type(disk_zone_type), dimension(:), allocatable, target :: disk_zone
  type(disk_region_type), dimension(:), allocatable, target :: regions
  type(cavity_type) :: cavity

  real(kind=dp) :: Rmin, Rmax, Rmax2, diskmass, correct_Rsub

  real :: r_subdivide
  logical :: lr_subdivide

  ! La densite(z=0) evolue comme r**(surf-exp_beta)

  ! Dullemond et Dubrulle
  real :: alpha

  ! Vitesses (en m/s) (Dullemond et Dubrulle)
  real, parameter :: v_sound = 380.0
!  real, parameter :: dv_gaz = 0.02
  ! Rayon de definition des vitesses
  real, parameter :: rref_v = 50.

  !! Loi de puissance
  real :: exp_strat, a_strat! = 0.7
  !! Loi exp (Garaud , Barriere 2004)
!!  real, parameter :: fact_strat = 0.3

  ! Grille
  real(kind=dp), parameter :: prec_grille=1.0e-14_dp
  real(kind=dp), parameter :: prec_grille_sph=1.0e-10_dp

  real(kind=dp), dimension(:,:,:), allocatable :: disk_origin
  real(kind=dp), dimension(:,:), allocatable :: star_origin
  real(kind=dp) :: frac_star_origin

  real, dimension(:), allocatable :: E_disk !n_lambda

  real :: w0_sup, w0_inf

  real :: z_warp, tilt_angle

  ! Description analytique du puffed-up inner rim
  real :: puffed_rim_h, puffed_rim_r, puffed_rim_delta_r
  logical :: lpuffed_rim

  character(len=512) :: density_file, sigma_file, grain_size_file, limits_file
  character(len=512), dimension(:), allocatable :: sh_file

  ! Correction locale de la desnite (dans un anneau)
  real :: correct_density_factor, correct_density_Rin, correct_density_Rout

  logical :: lgap_Gaussian
  real :: f_gap_Gaussian, r_gap_Gaussian, sigma_gap_Gaussian

  ! Fichier de Gasp pour la structure du disque
  real :: struct_file_rin, struct_file_rout, struct_file_zmax, struct_file_beta
  real :: struct_file_amin, struct_file_amax,  struct_file_rref
  integer :: struct_file_n_grains, struct_file_nspecies

end module disk

!********************************************************************

module prop_star

  use parametres

  implicit none
  save

  integer :: n_etoiles

  type star_type
     real :: r, T, M, fUV, slope_UV, othin_sublimation_radius
     real(kind=dp) :: x,y,z
     logical :: lb_body, out_model
     character(len=512) :: spectre
     integer :: icell
  end type star_type

  type(star_type), dimension(:), allocatable :: etoile
  real, dimension(:,:), allocatable :: CDF_E_star, prob_E_star

  real, dimension(:), allocatable :: E_stars !n_lambda

  real :: L_etoile

  real, dimension(:), allocatable :: E_ISM
  real, parameter :: R_ISM = 1.5 ! rayon de la sphere d'ou est emis le champ ISM

  ! Spot
  real :: T_spot, surf_fraction_spot, theta_spot, phi_spot

  ! Limb darkening
  logical :: llimb_darkening
  character(len=512) :: limb_darkening_file
  real, dimension(:), allocatable :: mu_limb_darkening, limb_darkening, pola_limb_darkening

  character(len=8) :: system_age

end module prop_star

!********************************************************************

module grains
!****************************
! Parametres des grains
!****************************
  use parametres

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
  integer :: n_pop
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
  real(kind=dp), dimension(:), allocatable :: tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda !n_lambda
  real, dimension (:,:), allocatable :: tab_amu1, tab_amu2 !n_lambda,n_pop
  real, dimension (:,:), allocatable :: tab_amu1_coating, tab_amu2_coating !n_lambda,n_pop
  real, dimension(:), allocatable :: tab_lambda2, tab_lambda2_inf, tab_lambda2_sup, tab_delta_lambda2

  ! Parametres de diffusion des grains
  real, dimension(:,:,:), allocatable :: tab_s11, tab_s12, tab_s33, tab_s34, prob_s11 !n_lambda,n_grains,180
  real, dimension(:,:), allocatable :: tab_g, tab_albedo, C_ext, C_sca, C_abs, C_abs_norm !n_grains, n_lambda
  !real, dimension(:), allocatable :: q_geo ! n_grains section geometrique en m^2

  ! aggregats
  real, dimension(:,:,:,:,:), allocatable :: tab_mueller !4,4, 180, n_grains,n_lambda

  ! Parametres de diffusion des cellules
  real, dimension(:,:), allocatable :: tab_albedo_pos, tab_g_pos ! n_cells,n_lambda
  real, dimension(:,:,:), allocatable :: tab_s11_pos, tab_s12_o_s11_pos, tab_s33_o_s11_pos, tab_s34_o_s11_pos, prob_s11_pos ! 0:180, n_cells,n_lambda

  character(len=512) :: aggregate_file, mueller_aggregate_file
  real :: R_sph_same_M

end module grains

!********************************************************************

module naleat
!**********************************************************
! Parametres des nombres aleatoires
!**********************************************************
! SPRNG renvoie un nbre aleatoire dans [0,1[
! faire attention a conversion dp -> sp sur machine 64 bits
! qui peut donner 1.0 -> BUG
!**********************************************************
  use parametres, only : dp

  implicit none
  save

#include "sprng_f.h"

  integer :: seed = 269753

  integer :: gtype = 2
  ! 0 -> Additive Lagged Fibonacci Generator
  ! 1 -> 48 bit Linear Congruential Generator with Prime Addend
  ! 2 -> 64 bit Linear Congruential Generator with Prime Addend
  ! 3 -> Combined multiple recursive generator
  ! 4 -> Multiplicative Lagged Fibonacci Generator

  SPRNG_POINTER, dimension(:), allocatable :: stream

  ! Generateur gaussien
  real(kind=dp), dimension(:), allocatable:: gauss_random_saved
  logical, dimension(:), allocatable :: lgauss_random_saved

  interface
     function gauss_random(id)
       use parametres
       integer, intent(in) :: id
       real(kind=dp) :: gauss_random
     end function gauss_random
  end interface

end module naleat

!********************************************************************

module opacity

  use disk
  use parametres
  use grains

  implicit none
  save

  real(kind=dp) :: zmaxmax
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

  real(kind=dp), dimension(:,:), allocatable :: kappa !n_cells, n_lambda
  real, dimension(:,:), allocatable :: kappa_abs_LTE, kappa_abs_nLTE, kappa_sca, kappa_abs_RE ! n_cells, n_lambda
  real, dimension(:,:), allocatable :: proba_abs_RE, proba_abs_RE_LTE, Proba_abs_RE_LTE_p_nLTE
  real, dimension(:,:,:), allocatable :: kabs_nLTE_CDF, kabs_nRE_CDF ! 0:n_grains, n_cells, n_lambda
  real(kind=dp), dimension(:,:), allocatable :: emissivite_dust ! emissivite en SI (pour mol)

  real(kind=dp), dimension(:,:), allocatable :: densite_pouss ! n_grains, n_cells en part.cm-3
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

module resultats
  use parametres
  use parametres
  use opacity, only : zmax

  implicit none
  save

  integer :: ntot = 0, nkill = 0
  real :: flux_etoile, flux_direct

  real(kind=dp), dimension(:,:,:,:,:,:), allocatable :: STOKEI, STOKEQ, STOKEU, STOKEV !id,lambda,x,y,n_thet,n_phi
  real(kind=dp), dimension(:,:,:,:,:,:), allocatable :: STOKEI_star, STOKEI_star_scat, STOKEI_disk, STOKEI_disk_scat
  real, dimension(:,:,:,:,:), allocatable :: stoke_io ! x,y, theta, phi, type

  real(kind=dp), dimension(:,:,:,:), allocatable :: sed, n_phot_sed, sed_u, sed_q, sed_v
  real(kind=dp), dimension(:,:,:,:), allocatable, target :: sed_star, sed_star_scat, sed_disk, sed_disk_scat!id,lambda,n_thet,n_phi,x,y
  real, dimension(:,:,:), allocatable :: sed1_io
  real, dimension(:,:,:,:), allocatable :: sed2_io
  real, dimension(:,:), allocatable :: wave2_io
  real(kind=dp), dimension(:,:), allocatable :: n_phot_envoyes

  ! Line transfer
  real, dimension(:,:,:,:,:,:), allocatable :: spectre ! speed,trans,thetai,phi,x,y
  real, dimension(:,:,:,:,:), allocatable :: continu ! trans,thetai,phi,x,y

end module resultats

!********************************************************************

module em_th

  use parametres

  implicit none
  save

  real :: T_max, T_min, Tmax_PAH ! Temp_sublimation et Temp nuage

  ! Gamme de longueurs d'onde utilisees
  real :: lambda_min, lambda_max = 3000.0
  ! Echelle_log
  ! Nbre d'intervalle de longueur d'onde
  real(kind=dp) :: delta_lambda

  ! Proba cumulee en lambda d'emettre selon un coprs noir
  real, dimension(:), allocatable :: spectre_etoiles_cumul, spectre_etoiles !(0:n_lambda)
  real, dimension(:), allocatable :: spectre_emission_cumul !(0:n_lambda)

  ! Nbre de temperatures utilisee pour pretabule kplanck et prob_delta_T
  ! et temperature echantillonee
  integer  :: n_T
  real, dimension(:), allocatable :: tab_Temp
!  real, dimension(n_rad,nz,n_T) :: spline_Temp

  ! fraction d'energie reemise sur energie etoile
  ! (Opacite moyenne de Planck * coeff)
  real, dimension(:,:), allocatable :: log_frac_E_em ! 0:n_T, n_cells
  real, dimension(:,:), allocatable :: log_frac_E_em_1grain  !n_grains,0:n_T
  real, dimension(:,:), allocatable :: frac_E_em_1grain_nRE, log_frac_E_em_1grain_nRE !n_grains,0:n_T

  ! Probabilite cumulee en lambda d'emissivite de la poussiere
  ! avec correction de temperature (dp/dT)
  ! (Bjorkman & Wood 2001, A&A 554-615 -- eq 9)
  real, dimension(:,:,:), allocatable :: kdB_dT_CDF ! 0:n_T,n_cells,n_lambda
  real, dimension(:,:,:), allocatable :: kdB_dT_1grain_LTE_CDF, kdB_dT_1grain_nLTE_CDF, kdB_dT_1grain_nRE_CDF ! n_grains,0:n_T,n_lambda

  ! pour stockage des cellules par lequelles on passe
  ! longueur de vol cumulee dans la cellule
  real(kind=dp), dimension(:,:), allocatable :: xKJ_abs, nbre_reemission ! n_cells, id
  real(kind=dp), dimension(:), allocatable :: E0 ! n_cells
  real(kind=dp), dimension(:,:), allocatable :: J0 !n_cells, n_lambda, n_rad, nz, n_az
  ! xJabs represente J_lambda dans le bin lambda -> bin en log : xJabs varie comme lambda.F_lambda
  real(kind=dp), dimension(:,:,:), allocatable :: xJ_abs ! n_cells, n_lambda, nb_proc
  real, dimension(:,:,:), allocatable :: xN_abs ! n_cells, n_lambda, nb_proc

  integer, dimension(:,:), allocatable :: xT_ech ! n_cells, id
  integer, dimension(:,:,:), allocatable :: xT_ech_1grain, xT_ech_1grain_nRE ! n_grains, n_cells, id

  real(kind=dp) :: E_abs_nRE, E_abs_nREm1
  ! emissivite en unite qq (manque une cst mais travail en relatif)
  real(kind=dp), dimension(:,:), allocatable :: Emissivite_nRE_old ! n_lambda, n_rad, nz, n_az

  real, dimension(:), allocatable :: Temperature, Temperature_old !n_rad,nz,n_az
  real, dimension(:,:), allocatable :: Temperature_1grain, Temperature_1grain_nRE !n_rad,nz, n_grains
  real, dimension(:,:), allocatable :: Temperature_1grain_old, Temperature_1grain_nRE_old, maxP_old !n_rad,nz, n_grains
  integer, dimension(:,:), allocatable :: Tpeak_old
  real, dimension(:,:,:), allocatable :: Proba_Temperature !n_T, n_cells,, n_grains
  logical, dimension(:,:), allocatable :: l_RE, lchange_nRE ! n_grains, n_cells
  real :: nbre_photons_tot, L_packet_th


  ! Choix cellule d'emission pour cas monochromatique
  real(kind=dp), dimension(:,:), allocatable :: prob_E_cell !n_lambda,0:n_rad*nz
  real, dimension(:), allocatable :: frac_E_stars, frac_E_disk, E_totale !n_lambda

  ! Biais de l'emission vers la surface du disque
  real, dimension(:), allocatable :: weight_proba_emission, correct_E_emission

  ! Suppresion de grains
  integer :: specie_removed
  real :: T_rm

  character(len=512) :: Tfile = "./data_th/Temperature.fits.gz"
  character(len=512) :: Tfile_nLTE = "./data_th/Temperature_nLTE.fits.gz"
  character(len=512) :: Tfile_Diff_approx = "./data_th/Temperature_Diff_approx.fits.gz"
  character(len=512) :: Tfile_nRE = "./data_th/Temperature_nRE.fits.gz"

  real, dimension(:,:,:), allocatable :: kappa_abs_tmp

end module em_th

!********************************************************************

module constantes

  use parametres

  implicit none
  save

  ! Quelques reels utiles
  real(kind=dp), parameter :: pi =3.141592653589793238462643383279502884197_dp ! ca devrait etre bon la
  real(kind=dp), parameter :: deux_pi = 2.0_dp * pi
  real(kind=dp), parameter :: un_sur_deux_pi = 1.0_dp/deux_pi
  real(kind=dp), parameter :: quatre_pi = 4.0_dp * pi
  real(kind=dp), parameter :: quatre_tiers_pi = 4.0_dp/3.0_dp * pi
  real(kind=dp), parameter :: pi_sur_deux = 0.5_dp * pi
  real(kind=dp), parameter :: un_tiers = 1.0_dp / 3.0_dp

  ! Constantes en SI !!!!!!!!
  real, parameter :: hp = 6.6260693e-34  ! Planck (J.Hz-1)
  real, parameter :: kb = 1.3806505e-23  ! Boltzmann (J.K^-1)
  real, parameter :: c_light = 299792458. ! vitesse lumiere (m.s^-1)
  real, parameter :: cst_th=c_light*hp/kb   ! pour calcul de (h c)/(lambda k T)
  real, parameter :: sigma = 5.6697e-8 ! Stefan (en W/(m^2.K^4))
  real, parameter :: Ggrav = 6.672e-11 ! (m^3.s^-2.kg^-1)    e-8 en cgs

  real, parameter :: mole = 6.022e23   ! Nombre d'Avogadro
  real, parameter :: masseH = 1.0/mole ! masse d'un atome d'hydrogene en g
  real, parameter :: mu = 2.3 ! en g,  2.3 selon Walker 2004
  real, parameter :: masse_mol_gaz = mu * masseH
  real, parameter :: T_Cmb = 2.73

  ! Changements d'unites
  ! Angles
  real(kind=dp), parameter :: deg_to_rad = pi/180.0_dp
  real(kind=dp), parameter :: rad_to_deg = 1.0/deg_to_rad
  real(kind=dp), parameter :: arcsec_to_rad = deg_to_rad / 3600.
  real(kind=dp), parameter :: rad_to_arcsec = 1.0/arcsec_to_rad
  real(kind=dp), parameter :: arcsec_to_deg = 1. / 3600.
  real(kind=dp), parameter :: deg_to_arcsec = 3600.

  ! Longueurs
  real(kind=dp), parameter :: AU_to_m = 149597870700._dp ! IAU 2012 definition
  real(kind=dp), parameter :: m_to_AU = 1.0_dp/AU_to_m

  real(kind=dp), parameter :: AU_to_cm = AU_to_m * 100._dp
  real(kind=dp), parameter :: cm_to_AU = 1.0_dp/AU_to_cm

  real(kind=dp), parameter :: AU3_to_m3 = AU_to_m**3
  real(kind=dp), parameter :: AU3_to_cm3 = AU_to_cm**3

  real(kind=dp), parameter :: mum_to_m = 1.0e-6_dp
  real(kind=dp), parameter :: m_to_mum = 1.0e6_dp
  real(kind=dp), parameter :: mum_to_cm = 1.0e-4_dp
  real(kind=dp), parameter :: cm_to_mum = 1.0e4_dp

  real(kind=dp), parameter :: A_to_mum = 1.0e-4_dp

  real(kind=dp), parameter :: m_to_cm = 1.0e2_dp
  real(kind=dp), parameter :: m3_to_cm3 = m_to_cm**3
  real(kind=dp), parameter :: cm_to_m = 1.0e-2_dp
  real(kind=dp), parameter :: m_to_km = 1.0e-3_dp
  real(kind=dp), parameter :: km_to_m = 1.0e3_dp

  real(kind=dp), parameter :: Rsun_to_AU = 0.00466666666_dp
  real(kind=dp), parameter :: Au_to_Rsun = 1.0_dp/Rsun_to_AU

  real(kind=dp), parameter :: pc_to_AU = 206264.81_dp
  real(kind=dp), parameter :: rad_to_sec = pc_to_AU
  real(kind=dp), parameter :: AU_to_pc = 1.0/pc_to_AU
  real(kind=dp), parameter :: sec_to_rad = AU_to_pc

  ! Energies
  real(kind=dp), parameter :: eV_to_J = 1.60217653e-19_dp
  real(kind=dp), parameter :: erg_to_J = 1.0e-7_dp
  real, parameter :: jansky = 1.0e-26 ! W.m^-2.Hz-1 (F_nu en jansky)

  ! Masses
  real(kind=dp), parameter :: Msun_to_g = 1.9891e33_dp
  real(kind=dp), parameter :: g_to_Msun = 1.0_dp/Msun_to_g

  real(kind=dp), parameter :: Msun_to_kg = 1.9891e30_dp
  real(kind=dp), parameter :: kg_to_Msun = 1.0/Msun_to_kg

  real(kind=dp), parameter :: g_to_kg = 1.0e-3_dp
  real(kind=dp), parameter :: kg_to_g = 1.0e3_dp

  ! Limites de precision numerique
  real, parameter :: tiny_real = tiny(0.0)
  real, parameter :: huge_real = huge(1.0)
  real(kind=dp), parameter :: tiny_dp = tiny(0.0_dp)
  real(kind=dp), parameter :: huge_dp = huge(1.0_dp)

  real, parameter ::  tiny_real_x1e6 =  tiny_real * 1.0e6

  integer, parameter :: huge_integer = huge(1)
  real , parameter :: max_int = real(huge_integer) * (1.0-1.0e-5)

end module constantes

!***********************************************************

module molecular_emission

  use parametres
  use constantes

  implicit none
  save

  logical :: ldouble_RT

  ! Lecture du fichier de data mol
  real :: molecularweight
  real :: abundance
  integer :: nLevels
  real(kind=dp) :: largeur_profile
  real, dimension(:), allocatable :: Level_energy, j_qnb
  ! g est dp car les calculs utilisant g sont en dp
  real(kind=dp), dimension(:), allocatable :: poids_stat_g
  integer :: nTrans_tot, nTrans
  integer, dimension(:), allocatable :: indice_Trans

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
  real, dimension(:), allocatable :: vfield ! n_cells
!  real, dimension(:,:,:), allocatable :: vx, vy
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

  logical :: linfall, lkeplerian

  real(kind=dp), dimension(:), allocatable :: deltaVmax ! n_cells
  real(kind=dp), dimension(:,:), allocatable :: tab_deltaV ! n_speed, n_cells
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

  use parametres, only : dp
  use grains, only : nang_scatt

  implicit none
  save

  logical :: lscatt_ray_tracing, lscatt_ray_tracing1, lscatt_ray_tracing2, loutput_mc
  ! ray-tracing 1 : sauve le champ de radiation diffuse
  ! ray-tracing 2 : sauve l'intensite specifique


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
  real, dimension(:,:,:), allocatable :: stars_map ! nx, ny, 4

  real, dimension(:,:,:,:,:), allocatable :: weight_Inu_fct_phase ! n_rayon_rt, dir, n_theta_I, n_phi_I, nang_scatt

end module ray_tracing

!********************************************************************
