module parametres
!****************************
! Parametres generaux
!****************************

  use mcfost_env, only : dp

  implicit none
  save

  real :: para_version

  logical :: lpara, lstop_after_init
  integer :: indice_etape, etape_i, etape_f

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

  logical :: lRE_LTE, lRE_nLTE, lnRE, lonly_LTE, lonly_nLTE
  logical :: loutput_J, loutput_J_step1, loutput_UV_field, lxJ_abs, lxJ_abs_step1

  ! Methode de calcul de la diffusion : a choisir pour optimiser taille memoire et temps cpu
  ! 0 -> automatique
  ! 1 -> choix taille du grain diffuseur + matrice Mueller par grain
  ! 2 -> matrice de Mueller moyenne par cellule (benchmark)
  integer :: scattering_method, scattering_method0
  logical :: lscattering_method1

  ! Theorie de Mie ou HG
  integer :: aniso_method ! 1 = full phase function, 2 = HG
  logical :: lmethod_aniso1

  integer :: RT_sed_method ! cf routine dust_map pour def

  ! Etapes de l'émission thermique
  logical :: ltemp, lsed, lsed_complete, l_em_disk_image, lchauff_int, lextra_heating, lno_internal_energy
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
  integer ::  npix_x, npix_y, npix_x_save, npix_y_save
  real :: angle_interet, zoom, tau_seuil, wl_seuil

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
  logical :: ldensity_file, lsigma_file, lvelocity_file, lphantom_file, lgadget2_file, lascii_SPH_file, llimits_file
  logical :: lweight_emission, lcorrect_density, lProDiMo2mcfost, lProDiMo2mcfost_test, lastrochem
  logical :: lspot, lforce_PAH_equilibrium, lforce_PAH_out_equilibrium, lchange_Tmax_PAH, lISM_heating, lcasa
  integer :: ISR_model ! 0 : no ISM radiation field, 1 : ProDiMo, 2 : Bate & Keto

  ! benchmarks
  logical :: lbenchmark_Pascucci, lbenchmark_vanZadelhoff1, lbenchmark_vanZadelhoff2, lDutrey94, lHH30mol
  logical :: lbenchmark_water1, lbenchmark_water2, lbenchmark_water3, lbenchmark_SHG, lMathis_field
  real :: Mathis_field

  ! Prodimo
  logical :: lprodimo, lprodimo_input_dir, lforce_ProDiMo_PAH

  logical, parameter :: ltest_rt3 = .false. ! marche pas
  logical, parameter :: ltest_rt4 = .false.  ! marche pas non plus

  logical :: lSeb_Charnoz, lread_Seb_Charnoz, lread_Seb_Charnoz2, lread_Misselt, lread_DustEM
  logical :: lread_grain_size_distrib, lphase_function_file,ltau1_surface

  ! Phantom
  logical :: ldudt_implicit
  real(kind=dp) :: ufac_implicit

end module parametres
