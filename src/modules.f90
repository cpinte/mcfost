module parametres
!****************************
! Parametres generaux
!****************************
  
  implicit none
  save

  real, parameter :: mcfost_version = 2.14
  character(8), parameter :: mcfost_release = "2.14.5"
  !character(len=128), parameter :: webpage="http://www-laog.obs.ujf-grenoble.fr/public/pintec/mcfost/"
  character(len=128), parameter :: webpage=      "http://ipag.osug.fr/public/pintec/mcfost/"
  character(len=128), parameter :: utils_webpage="http://ipag.osug.fr/public/pintec/mcfost_utils/"
  

  real :: para_version

  ! Système
  integer :: nb_proc, checkpoint, checkpoint_time, checkpoint_level
  logical :: lpara, lcheckpoint
  integer, parameter :: sl = selected_real_kind(p=6,r=37)
  integer, parameter :: db = selected_real_kind(p=13,r=200)
  integer, parameter :: limite_stack = 5000000
  integer :: indice_etape, etape_i, etape_f
  integer :: time_begin, time_end, time_tick, time_max
  
  ! Nombre de photons lances
  logical :: ldust_transfer
  integer :: nbre_photons_loop, nbre_photons_eq_th, nbre_photons_lambda, nbre_photons_image, nbre_photons_spectre
  real :: nbre_photons_lim = 1.e4 ! combien de fois plus on aurait recu sans disque
  integer :: nnfot1
  real(kind=db), dimension(:), allocatable, target :: nnfot2
  real(kind=db) :: E_paquet
  integer :: n_dif_max_eq_th = 100000 ! Nbre max de dif autorises dans calcul eq. th OUTDATED
  real :: tau_dark_zone_eq_th = 1500 !1500.   15000 pour benchmark a tau=1e6
  real :: tau_dark_zone_obs = 100 ! idem que 1000 si 1500. ci-dessus
  logical :: linteg_dic
  integer :: n_Stokes

  ! Nbre d'angles pour echantillonner la fct de phase 
  integer, parameter :: nang_scatt = 180  ! TODO : ca bug si c'est pas 180

  ! Nbre de longueurs d'onde utilisees
  integer :: n_lambda
  logical :: lmono0, lmono

  ! lstrat = true si on prend en compte la stratification
  logical :: lstrat, ldust_sublimation, lRE_LTE, lRE_nLTE, lnRE, loutput_J, loutput_UV_field, lxJ_abs, lxJ_abs2

  ! Methode de calcul de la diffusion : a choisir pour optimiser taille memoire et temps cpu
  ! 0 -> automatique
  ! 1 -> choix taille du grain diffuseur + matrice Mueller par grain
  ! 2 -> matrice de Mueller moyenne par cellule (benchmark)
  integer :: scattering_method
  logical :: lscattering_method1

  ! Theorie de Mie ou HG
  integer :: aniso_method 
  logical :: lmethod_aniso1 
  
  integer :: RT_sed_method ! cf routine dust_map pour def

  ! Etapes de l'émission thermique
  logical :: ltemp, lsed, lsed_complete, l_em_disk_image, lchauff_int
  character(len=512), dimension(:), allocatable :: indices
  character(len=512) :: tab_wavelength

  ! Emission moleculaire
  logical :: lemission_mol,  lpop, lprecise_pop, lmol_LTE, ldust_mol, lsetup_gas, lonly_top, lonly_bottom
  
  ! Decomposition image
  logical :: lsepar_contrib, lsepar_pola, lonly_capt_interet
  integer :: N_type_flux 
  ! les flux sont I, (Q,U,V), (star, scatt, disk th, disk th scatt.)

  ! rotation du plan du disque en deg., sens trigo.
  real(kind=db) :: ang_disque
  real(kind=db) :: sin_disk, cos_disk

  ! Production d'images symetriques
  ! La symetrie est effectuee avant de choisir les pixels
  ! le système est-il centrosymetrique
  ! le systeme a-t-il une symetrie axiale (ne compte que si N_phi > 1)
  logical :: l_sym_ima, l_sym_centrale, l_sym_axiale

  
  ! Parametres des cartes
  integer :: N_thet, N_incl, N_phi, capt_interet, delta_capt, capt_inf, capt_sup, capt_debut, capt_fin
  integer ::  igridx, igridy, deltapix_x, deltapix_y, maxigrid
  integer, dimension(:), allocatable ::  igridx2, igridy2, deltapix_x2, deltapix_y2, maxigrid2
  real :: angle_interet, zoom, size_pix, tau_seuil, wl_seuil

  real, dimension(:), allocatable :: size_pix2
  real, parameter :: cutoff = 7.0


  ! Résolution de la grille de densité
  ! Nombre de cellules dans la direction r (echantillonage log) 
  integer :: grid_type ! 1 = cylindrical, 2 = spherical
  integer :: n_rad, n_rad_in  ! subdivision de la premiere cellule
  ! Nombre de couches verticales ( = de stratifications)
  integer :: nz, p_n_rad, p_nz, p_n_az, p_n_lambda
  ! Nombre de cellules azimuthales
  integer :: n_az, j_start, pj_start
  logical :: llinear_grid

  integer :: n_lambda2

  integer :: n_cartes
  integer, parameter :: n_cartes_max = 4

  logical :: letape_th, limg, lgap, lpah, lorigine, laggregate, l3D, lremove, lwarp, lcavity
  logical :: lstrat_SPH, lno_strat_SPH, lstrat_SPH_bin, lno_strat_SPH_bin, loutput_density_grid
  logical :: lopacite_only, lseed, ldust_prop, lopacity_map, lreemission_stats
  logical :: lapprox_diffusion, lcylindrical, lspherical, is_there_disk, lno_backup, lonly_diff_approx, lforce_diff_approx
  logical :: laverage_grain_size, lisotropic, lno_scattering, lqsca_equal_qabs, lgap_laure
  logical :: lkappa_abs_grain, ldust_gas_ratio
  logical :: lweight_emission, lcorrect_density, lProDiMo2mcfost, lProDiMo2mcfost_test, lLaure_SED, lforce_T_Laure_SED
  logical :: lspot, lforce_1st_scatt

  character(len=512) :: mcfost_utils, home, data_dir, root_dir, basename_data_dir, seed_dir
  character(len=512) :: dust_dir, mol_dir, star_dir, lambda_dir, lambda_filename
  character(len=512) :: para, band, model_pah, pah_grain, cmd_opt
  character(len=512), dimension(100) :: data_dir2, basename_data_dir2

  ! benchmarks
  logical :: lbenchmark_Pascucci, lbenchmark_vanZadelhoff1, lbenchmark_vanZadelhoff2, lDutrey94, lHH30mol
  logical :: lbenchmark_water1, lbenchmark_water2, lbenchmark_water3

  ! Disque de debris calcule par code dynamique
  ! ldebris -> using a text file for input density table
  ! lfits -> using a FITS file for input density table
  ! Can't do both! This is checked within read_param
  logical :: ldebris, lfits
  character(len=512) :: debris_file, struct_fits_file

  character(len=512) :: Laure_SED_filename

  ! Prodimo
  logical :: lprodimo, lprodimo_input_dir, lforce_ProDiMo_PAH

  logical, parameter :: ltest_rt3 = .false. ! marche pas
  logical, parameter :: ltest_rt4 = .false.  ! marche pas non plus

  ! Simu Sebastien Fromang
  logical :: lSeb_Fromang
  integer :: Seb_Fromang_model
  logical :: lSeb_Charnoz, lread_Seb_Charnoz, lread_Seb_Charnoz2

end module parametres

!********************************************************************

module disk
!****************************
! Parametres du disque
!****************************

  use parametres

  implicit none
  save

!*     size_neb est la dimension maxi de la nebuleuse spherique
! Ne pas modifier par n_rad
  real :: distance ! Distance du disque en pc

!* Parametres du disque
!* --------------------
!*
!* toutes les distances sont donnees en AU et la masse du disque
!* en masse solaire
!*     exp_beta.......... exposant decrivant le "flaring"
!*     surf.......... exposant de la distribution SURFACIQUE de densite
!*                    surf .neq. -2.0 !!!!!
!*     rin, rout..... limites interne et externe du disque
!*     rref.......... rayon de reference ("r0")
!*     sclht......... echelle de hauteur a r0 ("h0")
!*
  real(kind=db) :: size_neb
  
  type disk_zone_type 
     real(kind=db) :: Rin, Rmin, Rc, Rout, Rref, edge, exp_beta, surf, sclht, diskmass
     integer :: geometry ! 1=disk, 2=tappered-disk, 3=envelope
  end type disk_zone_type

  type cavity_type
     real(kind=db) ::  exp_beta, sclht, rref
  end type cavity_type

  type(disk_zone_type), dimension(:), allocatable, target :: disk_zone
  type(cavity_type) :: cavity

  real(kind=db) :: rmin, rout, rout2, diskmass, correct_Rsub

  real :: r_subdivide
  logical :: lr_subdivide

  ! La densite(z=0) evolue comme r**(surf-exp_beta)

  ! Dullemond et Dubrulle
  real :: alpha

  ! Ratio gaz poussiere
  real :: gas_dust

  ! Vitesses (en m/s) (Dullemond et Dubrulle)
  real, parameter :: v_sound = 380.0
!  real, parameter :: dv_gaz = 0.02
  ! Rayon de definition des vitesses
  real, parameter :: rref_v = 50.

  !! Loi de puissance
  real :: exp_strat, a_strat! = 0.7
  !! Loi exp (Garaud , Barriere 2004)
!!  real, parameter :: fact_strat = 0.3
  
  ! angle sous-tendu par le disque
  real(kind=db) :: cos_max2, r_bord2
  
  ! Grille
  real(kind=db), parameter :: prec_grille=1.0e-14_db
  real(kind=db), parameter :: prec_grille_sph=1.0e-10_db
  
  ! Tableau des longueurs de vol
  integer :: n_cell_max
  integer, dimension(:), allocatable :: n_cell_traversees !n_cpu
  integer, dimension(:,:), allocatable :: tab_cell_r, tab_cell_z ! n_cpu, n_cellule 
  real(kind=db), dimension(:,:), allocatable :: tab_length, tab_tau, tab_length_tot ! n_cpu, n_cellule
  real, dimension(:,:), allocatable :: tab_x0, tab_y0, tab_z0 ! n_cpu, n_cellule

  ! Disque decrit en une ou n zones
    integer :: n_zones
  
  ! Definition des regions = zones deconnectees
  logical :: lold_grid
  integer :: n_regions
  integer, dimension(:), allocatable :: region ! n_zones
  real(kind=db), dimension(:), allocatable :: Rmin_region, Rmax_region ! n_regions
  

  real(kind=db), dimension(:,:,:,:), allocatable :: disk_origin
  real(kind=db), dimension(:,:), allocatable :: star_origin
  real(kind=db) :: frac_star_origin

  real, dimension(:), allocatable :: E_disk !n_lambda

  real :: w0_sup, w0_inf

  real :: z_warp

  ! Pas d'integration initial dans chaque cellule pour integration dichotomique
  real, dimension(:), allocatable :: delta0 !0:n_rad

  ! Description analytique du puffed-up inner rim
  real :: puffed_rim_h, puffed_rim_r, puffed_rim_delta_r
  logical :: lpuffed_rim

  character(len=512) :: density_file

  ! Correction locale de la desnite (dans un anneau)
  real :: correct_density_factor, correct_density_Rin, correct_density_Rout

  logical :: lgap_ELT
  real :: r_gap_ELT, sigma_gap_ELT

  ! Fichier de Gasp pour la structure du disque
  real :: struct_file_rin, struct_file_rout, struct_file_zmax, struct_file_beta
  real :: struct_file_amin, struct_file_amax,  struct_file_rref
  integer :: struct_file_n_grains, struct_file_nspecies
  
end module disk

!********************************************************************

module wall

  use parametres

  implicit none
  save

  logical :: lwall, lopacity_wall
  real(kind=db) :: h_wall, tau_wall, kappa_wall

end module wall

!********************************************************************

module prop_star
  
  use parametres

  implicit none 
  save

  integer :: n_etoiles

! Parametres de l'etoile
!  real, parameter :: r_etoile = 0.0086333333 ! 1.85*R_soleil
!  real, parameter :: r_etoile = 0.00466666666 ! R_soleil
!  real, parameter :: T_etoile = 4000.
!  real, parameter :: M_etoile = 0.5

! Parametres du point chaud
  ! rayon projete du spot
!  real, parameter :: r_spot = 0.00*r_etoile
!  real, parameter :: T_spot = 8000.  
!  real, parameter :: theta_spot = 20.*pi/180.
!  real, parameter :: phi_spot = 0.*pi/180.

  type star_type
     real :: r, T, M, x, y, z, fUV, slope_UV
     logical :: lb_body
     character(len=512) :: spectre
  end type star_type

  type(star_type), dimension(:), allocatable :: etoile
  real, dimension(:,:), allocatable :: prob_E_star

  real, dimension(:), allocatable :: E_stars !n_lambda

  real :: L_etoile

  real, dimension(:), allocatable :: I_ISM
  real, parameter :: R_ISM = 1.5 ! rayon de la sphere d'ou est emis le champ ISM

end module prop_star

!********************************************************************

module grains
!****************************
! Parametres des grains
!****************************
  use parametres

  implicit none
  save

!**********************************************************************
!*        INDICE REF.= N-IK (AMU1-IAMU2) ; A = RAYON GRAIN (MICRON)
!*        LA DISTRIBUTION EN TAILLE DES GRAINS EST PRISE SELON 
!*        MATHIS, RUMPL & NORDSIECK (1997):
!*                 N(A) = A^(-AEXP)    AVEC AEXP=3.5 A PRIORI
!*        LA DISTRIBUTION S'ETEND DE AMIN A AMAX (MICRONS).
!*
!*        RHO1G = densite grains  
!*        AMU1 ET AMU2 DIFFERENTS DE 0.0
!*
!* Parametres de la distribution de grains
!*
!* amin=0.005, amax=0.9, aexp=3.7 are ISM values from Mathis & Whiffen (1989)
!*

  type dust_grain_type
     real :: nbre_grains, a
     integer :: zone, methode_chauffage, pop
     logical :: is_PAH
  end type dust_grain_type

  type(dust_grain_type), dimension(:), allocatable  :: grain

  type dust_pop_type     
     integer :: n_grains, methode_chauffage, zone
     real :: amin, amax, aexp, frac_mass, porosity
     
     real :: avg_grain_mass, masse, sblow, rho1g_avg, T_sub
     
     integer ::  mixing_rule, n_components ! en cas de coating, composant 1 est le coeur 2 est le manteau
     ! mixing rule : 1 = EMT, 2 = coating
     character(len=512), dimension(10) :: indices
     real, dimension(10) :: component_rho1g, component_volume_fraction, component_T_sub
    
     logical :: is_PAH, lcoating
     integer :: ind_debut, ind_fin
     integer :: pop_geo
  end type dust_pop_type
  
  type(dust_pop_type), dimension(:), allocatable, target :: dust_pop

  ! Nombre de taille de grains 
  ! le pas et (amax - amin)/n_grains
  integer :: n_pop
  integer :: n_grains_tot , n_grains_RE_LTE, n_grains_RE_nLTE, n_grains_nRE
  integer :: grain_RE_LTE_start, grain_RE_LTE_end, grain_RE_nLTE_start, grain_RE_nLTE_end, grain_nRE_start, grain_nRE_end

  real,dimension(:), allocatable :: nbre_grains !n_grains_tot
  real, dimension(:), allocatable :: r_grain, r_core, S_grain, M_grain !n_grains_tot
  real, dimension(:,:), allocatable :: frac_mass_pop !n_zones, n_pop
  logical, dimension(:), allocatable :: is_pop_PAH, is_grain_PAH !n_pop et n_grains_tot
  
  ! Pour lecture des fichiers de PAH de B.Draine
  real, dimension(:,:,:), allocatable :: PAH_Q_ext, PAH_Q_abs, PAH_Q_sca, PAH_g
  integer, parameter :: PAH_n_rad = 30,  PAH_n_lambda = 1201
  real, dimension(PAH_n_rad) :: log_PAH_rad
  real, dimension(PAH_n_lambda) :: PAH_lambda, PAH_delta_lambda
  logical :: lread_PAH = .false.

  real, dimension(:,:), allocatable :: amin, amax, aexp ! n_zones, n_especes
  real :: wavel

  ! Tab de lambda
  real, dimension(:), allocatable :: tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda !n_lambda
  real, dimension (:,:), allocatable :: tab_amu1, tab_amu2 !n_lambda,n_pop
  real, dimension (:,:), allocatable :: tab_amu1_coating, tab_amu2_coating !n_lambda,n_pop
  real, dimension(:), allocatable :: tab_lambda2, tab_lambda2_inf, tab_lambda2_sup, tab_delta_lambda2

  ! Parametres de diffusion des grains
  real, dimension(:,:,:), allocatable :: tab_s11, tab_s12, tab_s33, tab_s34 !n_lambda,n_grains,180
  real, dimension(:,:), allocatable :: tab_g, tab_albedo, q_ext, q_sca, q_abs!n_lambda,n_grains
  real, dimension(:), allocatable :: q_geo ! n_grains section geometrique en m^2
  real, dimension(:,:,:), allocatable :: prob_s11 !n_lambda,n_grains,0:180
  
  ! aggregats
  real, dimension(:,:,:,:,:), allocatable :: tab_mueller !n_lambda,n_grains,4,4,180

  ! Parametres de diffusion des cellules
  real, dimension(:,:,:,:), allocatable :: tab_albedo_pos, tab_g_pos !n_lambda,n_rad,nz+1, (n_az)
  real, dimension(:,:,:,:,:), allocatable :: tab_s11_pos, tab_s12_pos!n_lambda,n_rad,nz+1,(n_az), 180
  real, dimension(:,:,:,:,:), allocatable ::  tab_s33_pos, tab_s34_pos !n_lambda,n_rad,nz+1,(n_az), 180
  real, dimension(:,:,:,:,:), allocatable :: prob_s11_pos !n_lambda,n_rad,nz+1,(n_az),0:180

  character(len=512) :: aggregate_file, mueller_aggregate_file
  real :: R_sph_same_M

end module grains

!********************************************************************

module naleat
!**********************************************************
! Parametres des nombres aleatoires
!**********************************************************
! SPRNG renvoie un nbre aleatoire dans [0,1[
! faire attention a conversion db -> sl sur machine 64 bits
! qui peut donner 1.0 -> BUG
!**********************************************************
  use parametres, only : db

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
  real(kind=db), dimension(:), allocatable:: gauss_random_saved
  logical, dimension(:), allocatable :: lgauss_random_saved 

  interface
     function gauss_random(id)
       use parametres
       integer, intent(in) :: id
       real(kind=db) :: gauss_random
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

  real(kind=db) :: zmaxmax
  real(kind=db), dimension(:), allocatable :: zmax !n_rad
  real(kind=db), dimension(:), allocatable :: volume !n_rad en AU^3
  real(kind=db), dimension(:,:,:), allocatable :: masse  !en g !!! n_rad, nz, n_az
  real, dimension(:,:), allocatable :: dust_gas_ratio
  real(kind=db), dimension(:,:), allocatable :: masse_rayon ! en g!!!!  n_rad, n_az
  real(kind=db), dimension(:), allocatable :: delta_z ! taille verticale des cellules cylindriques
  real(kind=db), dimension(:), allocatable :: dr2_grid ! differentiel en r^2 des cellules
  real(kind=db), dimension(:,:), allocatable :: r_grid, z_grid ! Position en cylindrique !!! des cellules
  real(kind=db), dimension(:), allocatable :: phi_grid
  real(kind=db), dimension(:), allocatable :: r_lim, r_lim_2, r_lim_3 ! lim rad sup de la cellule (**2) !0:n_rad
  real(kind=db), dimension(:,:), allocatable :: z_lim ! lim vert inf de la cellule !n_rad,nz+1
  real(kind=db), dimension(:), allocatable :: tan_phi_lim ! lim azimuthale de la cellule ! n_az
  real(kind=db), dimension(:), allocatable :: w_lim, theta_lim, tan_theta_lim ! lim theta sup de la cellule ! 0:nz
  

  real, dimension(:,:,:,:), allocatable :: amax_reel !n_lambda,n_rad,nz+1, (n_az)
  real(kind=db), dimension(:,:,:,:), allocatable :: kappa
  real(kind=db), dimension(:,:,:,:), allocatable :: kappa_abs_eg, kappa_sca !n_lambda,n_rad,nz+1, (n_az)
  real, dimension(:,:,:,:), allocatable :: proba_abs_RE !n_lambda,n_rad,nz+1, (n_az)
  real, dimension(:,:,:,:), allocatable :: prob_kappa_abs_1grain !n_lambda,n_rad,nz+1, 0:n_grains
  real(kind=db), dimension(:,:,:,:), allocatable :: emissivite_dust ! emissivite en SI (pour mol)

  real(kind=db), dimension(:,:,:,:), allocatable :: densite_pouss !n_rad,nz+1, (n_az), n_grains en part.cm-3
  integer :: ri_not_empty, zj_not_empty, phik_not_empty
!  real, dimension(n_lambda,n_rad,nz+1,0:n_grains) :: probsizecumul 
  real, dimension(:,:,:,:,:), allocatable :: probsizecumul !n_lambda,n_rad,nz+1,(n_az),)0:n_grains
  
  integer, parameter :: n_prob = 3
  ! proba_resol doit etre inferieur a 1/n_prob
  real, parameter :: proba_resol = 0.1/n_prob
  ! 0. et 1., proba_resol , 1.-proba_resol + n_prob-2 valeurs echantillonees regulierement
!  integer, dimension(n_lambda,n_rad,nz+1,0:n_prob+1) :: ech_prob 
!  real, dimension(n_lambda,n_rad,nz+1,0:n_prob+1) :: valeur_prob, xspline

  integer, dimension(:,:,:,:,:), allocatable :: ech_prob !n_lambda,n_rad,nz+1,0:n_prob+1
  real, dimension(:,:,:,:,:), allocatable :: valeur_prob !n_lambda,n_rad,nz+1,(n_az),0:n_prob+1
  real, dimension(:,:,:,:), allocatable :: xspline !n_lambda,n_rad,nz+1,0:n_prob+1

  logical :: l_is_dark_zone 
  logical, dimension(:,:,:), allocatable :: l_dark_zone !0:n_rad+1,0:nz+1, n_az
  real, dimension(:,:), allocatable :: r_in_opacite, r_in_opacite2 !nz+1, (n_az)
  integer, parameter :: delta_cell_dark_zone=3

  integer, dimension(:), allocatable :: ri_in_dark_zone, ri_out_dark_zone !n_az
  integer, dimension(:,:), allocatable :: zj_sup_dark_zone, zj_inf_dark_zone !n_rad, n_az

  logical, dimension(:,:), allocatable :: l_emission_pah !0:n_rad+1,0:nz+1

  real(kind=db), dimension(:,:,:), allocatable :: Dcoeff !n_rad, n_z, n_az
  real, dimension(:,:,:), allocatable :: DensE, DensE_m1 !n_rad, 0:n_z, n_az

  real, dimension(:), allocatable :: kappa_lambda,albedo_lambda,g_lambda
  real, dimension(:,:), allocatable :: pol_lambda_theta

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
  
  real(kind=db), dimension(:,:,:,:,:,:), allocatable :: STOKEI, STOKEQ, STOKEU, STOKEV !id,lambda,x,y,n_thet,n_phi
  real(kind=db), dimension(:,:,:,:,:,:), allocatable :: STOKEI_star, STOKEI_star_scat, STOKEI_disk, STOKEI_disk_scat
  real, dimension(:,:,:,:,:), allocatable :: stoke_io ! x,y, theta, phi, type
  real(kind=db), dimension(:,:,:,:,:,:), allocatable :: STOKEI1, STOKEI2, STOKEI3, STOKEI4 !id,lambda,x,y,n_thet,n_phi

  real(kind=db), dimension(:,:,:,:), allocatable :: sed, n_phot_sed, sed_u, sed_q, sed_v
  real(kind=db), dimension(:,:,:,:), allocatable, target :: sed_star, sed_star_scat, sed_disk, sed_disk_scat!id,lambda,n_thet,n_phi,x,y
  real, dimension(:,:,:), allocatable :: sed1_io
  real, dimension(:,:,:,:), allocatable :: sed2_io
  real, dimension(:,:), allocatable :: wave2_io
  real(kind=db), dimension(:,:,:,:), allocatable, target :: n_phot_sed2
  real(kind=db), dimension(:,:), allocatable :: n_phot_envoyes, n_phot_envoyes_loc

  ! Line transfer
  real, dimension(:,:,:,:,:), allocatable :: spectre ! speed,trans,thetai,x,y
  real, dimension(:,:,:,:), allocatable :: continu ! trans,thetai,x,y

end module resultats

!********************************************************************

module em_th
  use parametres

  implicit none 
  save
  
  real :: T_max, T_min ! Temp_sublimation et Temp nuage

  real :: L_bol0, L_bol1, L_bol2, L_tot, L_tot_E

  ! Energie d'un photon
  real :: E_photon

  ! Gamme de longueurs d'onde utilisees
  real :: lambda_min, lambda_max = 3000.0
  ! Echelle_log
  ! Nbre d'intervalle de longueur d'onde
  real(kind=db) :: delta_lambda

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
  real, dimension(:,:,:,:), allocatable :: log_frac_E_em !n_rad,nz+1,(n_az),0:n_T
  real, dimension(:,:), allocatable :: log_frac_E_em_1grain  !n_grains,0:n_T
  real, dimension(:,:), allocatable :: frac_E_em_1grain_nRE, log_frac_E_em_1grain_nRE !n_grains,0:n_T

  ! Probabilite cumulee en lambda d'emissivite de la poussiere
  ! avec correction de temperature (dB/dT)
  ! (Bjorkman & Wood 2001, A&A 554-615 -- eq 9)
  real, dimension(:,:,:,:,:), allocatable :: prob_delta_T !n_rad,nz+1,(n_az),0:n_T,n_lambda
  real, dimension(:,:,:), allocatable :: prob_delta_T_1grain !n_grains,0:n_T,n_lambda


  ! pour stockage des cellules par lequelles on passe
  ! longueur de vol cumulee dans la cellule
  real(kind=db), dimension(:,:,:,:), allocatable :: xKJ_abs, xE_abs, nbre_reemission !n_rad, nz, n_az, id
  real(kind=db), dimension(:,:,:), allocatable :: E0 !n_rad, nz, n_az
  real(kind=db), dimension(:,:,:,:), allocatable :: J0 !n_lambda, n_rad, nz, n_az
  real(kind=db), dimension(:,:,:,:), allocatable :: xE_abs_1grain !id, n_rad, nz, n_grains
  ! xJabs represente J_lambda dans le bin lambda -> bin en log : xJabs varie comme lambda.F_lambda 
  real(kind=db), dimension(:,:,:,:), allocatable :: xJ_abs !id, n_lambda, n_rad, nz 
  real, dimension(:,:,:,:), allocatable :: xN_abs !id, n_lambda, n_rad, nz 
  integer, dimension(:,:,:,:), allocatable :: xT_ech !id, n_rad, nz, n_az
  integer, dimension(:,:,:,:), allocatable :: xT_ech_1grain !id, n_rad, nz, n_grains

  real(kind=db) :: E_abs_nRE, E_abs_nREm1
  ! emissivite en unite qq (manque une cst mais travail en relatif)
  real(kind=db), dimension(:,:,:,:), allocatable :: Emissivite_nRE_old ! n_lambda, n_rad, nz, n_az

  real, dimension(:,:,:), allocatable :: Temperature, Temperature_old, Temperature_Laure_SED !n_rad,nz,n_az
  real, dimension(:,:,:), allocatable :: Temperature_1grain, Temperature_1grain_nRE !n_rad,nz, n_grains
  real, dimension(:,:,:), allocatable :: Temperature_1grain_old, Temperature_1grain_nRE_old, maxP_old !n_rad,nz, n_grains
  integer, dimension(:,:,:), allocatable :: Tpeak_old
  real, dimension(:,:,:,:), allocatable :: Proba_Temperature !n_T, n_rad,nz, n_grains
  logical, dimension(:,:,:), allocatable :: l_RE ! n_rad, nz, n_grains
  real :: nbre_photons_tot, n_phot_L_tot,  n_phot_L_tot0


  ! Choix cellule d'emission pour cas monochromatique
  real(kind=db), dimension(:,:), allocatable :: prob_E_cell !n_lambda,0:n_rad*nz
  real, dimension(:), allocatable :: frac_E_stars, E_totale !n_lambda
  
  ! Biais de l'emission vers la surface du disque
  real, dimension(:,:), allocatable :: weight_proba_emission, correct_E_emission

  ! Suppresion de grains
  integer :: specie_removed
  real :: T_rm

  character(len=512) :: Tfile = "./data_th/Temperature.fits.gz"
  character(len=512) :: Tfile_Diff_approx = "./data_th/Temperature_Diff_approx.fits.gz"
  character(len=512) :: Tfile_nRE = "./data_th/Temperature_nRE.fits.gz"

  real, dimension(:,:,:), allocatable :: kappa_abs_tmp 

end module em_th

!********************************************************************

module ki2

  use parametres
  implicit none
  save

  real :: sigma
  integer :: n_ker
  real, dimension(:,:), allocatable :: noyau, intensite
  
end module ki2

!***********************************************************

module ode_data

  implicit none
  save

  real :: coeff_exp, coeff1, coeff2, rcyl


end module ode_data

!***********************************************************

module constantes

  use parametres

  implicit none 
  save

  ! Quelques reels utiles
  real(kind=db), parameter :: pi =3.141592653589793238462643383279502884197_db ! ca devrait etre bon la
  real(kind=db), parameter :: deux_pi = 2.0_db * pi
  real(kind=db), parameter :: quatre_pi = 4.0_db * pi
  real(kind=db), parameter :: quatre_tiers_pi = 4.0_db/3.0_db * pi
  real(kind=db), parameter :: pi_sur_deux = 0.5_db * pi
  real(kind=db), parameter :: un_tiers = 1.0_db / 3.0_db

  ! Constantes en SI !!!!!!!!
  real, parameter :: hp = 6.6260693e-34  ! Planck (J.Hz-1)
  real, parameter :: kb = 1.3806505e-23  ! Boltzmann (J.K^-1)
  real, parameter :: c_light = 299792458. ! vitesse lumiere (m.s^-1)
  real, parameter :: cst_th=c_light*hp/kb   ! pour calcul de (h c)/(lambda k T)
  real, parameter :: sigma = 5.6697e-8 ! Stefan (en W/(m^2.K^4))
  real, parameter :: Ggrav = 6.672e-11 ! (m^3.s^-2.kg^-1)    e-8 en cgs

  real, parameter :: mole = 6.022e23   ! Nombre d'Avogadro
  real, parameter :: masseH = 1.0/mole ! masse d'un atome d'hydrogene en g
  real, parameter :: T_Cmb = 2.73

  ! Changements d'unites
  ! Angles
  real(kind=db), parameter :: deg_to_rad = pi/180.0_db
  real(kind=db), parameter :: rad_to_deg = 1.0/deg_to_rad
  real(kind=db), parameter :: arcsec_to_rad = deg_to_rad / 3600.
  real(kind=db), parameter :: rad_to_arcsec = 1.0/arcsec_to_rad
  real(kind=db), parameter :: arcsec_to_deg = 1. / 3600.
  real(kind=db), parameter :: deg_to_arcsec = 3600.

  ! Longueurs
  real(kind=db), parameter :: AU_to_m = 149597870700._db ! IAU 2012 definition
  real(kind=db), parameter :: m_to_AU = 1.0_db/AU_to_m

  real(kind=db), parameter :: AU_to_cm = AU_to_m * 100._db
  real(kind=db), parameter :: cm_to_AU = 1.0_db/AU_to_cm
  real(kind=db), parameter :: AU3_to_cm3 = AU_to_cm**3


  real(kind=db), parameter :: mum_to_m = 1.0e-6_db
  real(kind=db), parameter :: m_to_mum = 1.0e6_db
  real(kind=db), parameter :: mum_to_cm = 1.0e-4_db
  real(kind=db), parameter :: cm_to_mum = 1.0e4_db

  real(kind=db), parameter :: m_to_cm = 1.0e2_db
  real(kind=db), parameter :: m3_to_cm3 = m_to_cm**3
  real(kind=db), parameter :: cm_to_m = 1.0e-2_db
  real(kind=db), parameter :: m_to_km = 1.0e-3_db
  real(kind=db), parameter :: km_to_m = 1.0e3_db

  real(kind=db), parameter :: Rsun_to_AU = 0.00466666666_db
  real(kind=db), parameter :: Au_to_Rsun = 1.0_db/Rsun_to_AU

  real(kind=db), parameter :: pc_to_AU = 206264.81_db
  real(kind=db), parameter :: rad_to_sec = pc_to_AU
  real(kind=db), parameter :: AU_to_pc = 1.0/pc_to_AU
  real(kind=db), parameter :: sec_to_rad = AU_to_pc

  ! Energies
  real(kind=db), parameter :: eV_to_J = 1.60217653e-19_db 
  real(kind=db), parameter :: erg_to_J = 1.0e-7_db
  real, parameter :: jansky = 1.0e-26 ! W.m^-2.Hz-1 (F_nu en jansky)

  ! Masses
  real(kind=db), parameter :: Msun_to_g = 1.9891e33_db
  real(kind=db), parameter :: g_to_Msun = 1.0_db/Msun_to_g

  real(kind=db), parameter :: Msun_to_kg = 1.9891e30_db
  real(kind=db), parameter :: kg_to_Msun = 1.0/Msun_to_kg

  real(kind=db), parameter :: g_to_kg = 1.0e-3_db
  real(kind=db), parameter :: kg_to_g = 1.0e3_db

  ! Limites de precision numerique
  real, parameter :: tiny_real = tiny(0.0)
  real, parameter :: huge_real = huge(1.0)
  real(kind=db), parameter :: tiny_db = tiny(0.0_db)
  real(kind=db), parameter :: huge_db = huge(1.0_db)

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
  real(kind=db) :: largeur_profile
  real, dimension(:), allocatable :: Level_energy, j_qnb
  ! g est db car les calculs utilisant g sont en db
  real(kind=db), dimension(:), allocatable :: poids_stat_g 
  integer :: nTrans_tot, nTrans
  integer, dimension(:), allocatable :: indice_Trans
 
  real(kind=db), dimension(:), allocatable :: Aul, Blu, Bul, fAul, fBlu, fBul, transfreq
  integer, dimension(:), allocatable :: itransUpper, itransLower 
  integer :: nCollPart
  character(len=512), dimension(:), allocatable :: collBetween
  integer, dimension(:), allocatable :: nCollTrans, nCollTemps
  real, dimension(:,:), allocatable :: collTemps

  integer, dimension(:,:), pointer :: iCollUpper, iCollLower
  real, dimension(:,:,:), pointer :: collRates

  real, dimension(:,:,:), allocatable :: Tcin ! Temperature cinetique

  real :: nH2, masse_mol
  real, parameter :: masse_mol_gaz = 2. * masseH ! en g,  2.3 selon Walker 2004 
  ! masse_mol_gaz sert uniquement pour convertir masse disque en desnite de particule
  real(kind=db), dimension(:,:,:), allocatable :: frac_E_Trans ! n_rad, n_z, nTrans
  real(kind=db), dimension(:,:,:), allocatable :: kappa_mol_o_freq, kappa_mol_o_freq2 ! n_rad, nz, nTrans
  real(kind=db), dimension(:,:,:), allocatable :: emissivite_mol_o_freq,  emissivite_mol_o_freq2 ! n_rad, nz, nTrans
  real, dimension(:,:), allocatable :: vfield ! n_rad, nz
!  real, dimension(:,:,:), allocatable :: vx, vy
  real, dimension(:,:,:), allocatable :: tab_nLevel, tab_nLevel2, tab_nLevel_old ! n_rad, nz, nLevels 
 
  real, dimension(:,:), allocatable :: v_turb, v_line ! n_rad, nz

  real ::  vitesse_turb, dv, dnu
  integer, parameter :: n_largeur_Doppler = 15
  real(kind=db), dimension(:), allocatable :: tab_v ! n_speed

  real(kind=db), dimension(:,:,:), allocatable :: densite_gaz, masse_gaz ! n_rad, nz, n_az, Unites: part.m-3 et g : H2

  real(kind=db), dimension(:,:), allocatable :: ds 
  real(kind=db), dimension(:,:,:,:), allocatable :: I0, I02 ! nSpeed,nTrans,iray,ncpus
  real(kind=db), dimension(:,:,:), allocatable :: I0c ! Intensite dans le continu: nTrans,iray,ncpus
  real(kind=db), dimension(:,:,:), allocatable :: Doppler_P_x_freq

  real(kind=db), dimension(:,:), allocatable :: Jmol, Jmol2 ! nTrans, n_cpu
  real(kind=db), dimension(:), allocatable :: tab_Cmb_mol ! nTrans  

  logical :: linfall, lkeplerian

  real(kind=db), dimension(:,:), allocatable :: deltaVmax ! n_rad, nz
  real(kind=db), dimension(:,:,:), allocatable :: tab_deltaV ! n_speed, n_rad, nz
  real(kind=db), dimension(:,:), allocatable :: tab_dnu_o_freq ! n_rad, nz
  real(kind=db), dimension(:,:), allocatable :: norme_phiProf_m1, sigma2_phiProf_m1 ! n_rad, nz

  real, dimension(:,:), allocatable :: tab_abundance
  logical, dimension(:,:), allocatable :: lcompute_molRT

  logical ::  lfreeze_out
  real :: T_freeze_out

  real(kind=db), dimension(:,:,:,:,:), allocatable ::  origine_mol ! nv, nTrans, n_rad, nz, nb_proc

  integer :: RT_line_method, n_molecules

  type molecule
     integer :: n_speed, n_speed_rt, n_speed_center_rt, n_extraV_rt, nTrans_raytracing, iLevel_max
     real :: vmax_center_rt, extra_deltaV_rt, abundance
     logical :: lcst_abundance, lline
     character(len=512) :: abundance_file, filename
     character(len=32) :: name
     integer, dimension(100) :: indice_Trans_rayTracing
  end type molecule

  type(molecule), dimension(:), allocatable :: mol

  real(kind=db), dimension(:), allocatable :: tab_speed_rt 

  real, dimension(:,:,:), allocatable :: maser_map


end module molecular_emission

!***********************************************************

module ray_tracing

  use grains

  implicit none
  save

  logical :: lscatt_ray_tracing, lscatt_ray_tracing1, lscatt_ray_tracing2, loutput_mc

  ! inclinaisons
  real :: RT_imin, RT_imax
  integer ::  RT_n_ibin
  logical :: lRT_i_centered
  real, dimension(:), allocatable :: tab_RT_incl
  real(kind=db), dimension(:), allocatable :: tab_u_rt, tab_v_rt, tab_w_rt

  ! Sauvegarde champ de radiation pour rt2
  integer ::  n_phi_I,  n_theta_I ! 15 et 9 ok avec 30 et 30 en mode SED
 
  ! Pour rt 2 : nbre d'angle de visee en azimuth 
  ! TODO : calculer automatiquement en fct de la fct de phase + interpolation
  integer :: nang_ray_tracing, nang_ray_tracing_star    

  real, dimension(:,:,:,:), allocatable :: tab_s11_ray_tracing, tab_s12_ray_tracing, tab_s33_ray_tracing, tab_s34_ray_tracing ! n_lambda, n_rad, nz, nang_scatt
  real, dimension(:,:,:,:), allocatable :: tab_s12_o_s11_ray_tracing, tab_s33_o_s11_ray_tracing, tab_s34_o_s11_ray_tracing ! n_lambda, n_rad, nz, nang_scatt

  real, dimension(:,:,:), allocatable ::  cos_thet_ray_tracing, omega_ray_tracing ! nang_ray_tracing, 2 (+z et -z), nb_proc
  real, dimension(:,:,:), allocatable ::  cos_thet_ray_tracing_star, omega_ray_tracing_star ! nang_ray_tracing, 2 (+z et -z), nb_proc

  real, dimension(0:nang_scatt) :: tab_cos_scatt

  ! intensite specifique
  real, dimension(:,:), allocatable :: J_th ! n_rad, nz 

  ! methode RT 1
  integer, parameter :: n_az_rt = 45
  real, dimension(:,:,:,:,:,:,:), allocatable ::  xI_scatt ! 4, n_rad, nz, n_az_rt, 2, ncpus
  real, dimension(:,:,:,:,:,:), allocatable ::  xsin_scatt, xN_scatt ! n_rad, nz, n_az_rt, 2, ncpus
  real(kind=db), dimension(:,:,:), allocatable ::  I_scatt ! 4, n_az_rt, 2
  integer, dimension(:,:), allocatable :: itheta_rt1
  real(kind=db), dimension(:,:), allocatable ::  sin_omega_rt1, cos_omega_rt1, sin_scatt_rt1
  real(kind=db), dimension(:,:,:,:,:), allocatable ::  eps_dust1 !4,n_rad, nz, n_az_rt,0:1,4 
  
  ! methode RT 2
  real, dimension(:,:,:,:,:,:), allocatable :: xI ! 4, n_theta_I, n_phi_I, nrad, nz, ncpus
  real, dimension(:,:,:), allocatable :: xI_star, xw_star, xl_star ! nrad, nz, ncpus
  
  ! Fonction source: Ok en simple
  real, dimension(:,:,:,:,:,:), allocatable ::  I_sca2 ! 4, nang_ray_tracing, 2, n_rad, nz, ncpus
  real, dimension(:,:,:,:,:), allocatable ::  eps_dust2 ! 4, nang_ray_tracing, 2, n_rad, nz  
  real, dimension(:,:,:,:,:), allocatable ::  eps_dust2_star ! 4, nang_ray_tracing, 2, n_rad, nz  

  real, dimension(:,:,:,:,:,:), allocatable :: Stokes_ray_tracing ! n_lambda, nx, ny, RT_n_ibin, n_type_flux, ncpus

  real, dimension(:,:,:,:,:), allocatable :: weight_Inu_fct_phase ! n_rayon_rt, dir, n_theta_I, n_phi_I, nang_scatt

end module ray_tracing

!********************************************************************
