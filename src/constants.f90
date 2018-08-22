module constantes

  use parametres

  implicit none
  save

  ! Quelques reels utiles
  real(kind=dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp ! ca devrait etre bon la
  real(kind=dp), parameter :: deux_pi = 2.0_dp * pi
  real(kind=dp), parameter :: un_sur_deux_pi = 1.0_dp/deux_pi
  real(kind=dp), parameter :: quatre_pi = 4.0_dp * pi
  real(kind=dp), parameter :: quatre_tiers_pi = 4.0_dp/3.0_dp * pi
  real(kind=dp), parameter :: pi_sur_deux = 0.5_dp * pi
  real(kind=dp), parameter :: un_tiers = 1.0_dp / 3.0_dp

  ! Constantes en SI !!!!!!!!
  real, parameter :: hp = 6.6260693e-34    ! Planck (J.Hz-1)
  real, parameter :: kb = 1.3806505e-23    ! Boltzmann (J.K^-1)
  real, parameter :: c_light = 299792458.  ! vitesse lumiere (m.s^-1)
  real, parameter :: cst_th=c_light*hp/kb  ! pour calcul de (h c)/(lambda k T)
  real, parameter :: sigma = 5.6697e-8     ! Stefan (en W/(m^2.K^4))
  real, parameter :: Ggrav = 6.672e-11     ! (m^3.s^-2.kg^-1)    e-8 en cgs
  real, parameter :: electron_charge = 1.6021766208e-19  ! Coulombs

  real, parameter :: Na     = 6.022140857e23    ! Nombre d'Avogadro
  real, parameter :: amu    = 1.660531000E-24  ! atomar mass unit
  real, parameter :: masseH = 1.00794 * amu    ! masse d'un atome d'hydrogene en g
  real, parameter :: mu = 2.3                  ! en g,  2.3 selon Walker 2004
  real, parameter :: masse_mol_gaz = mu * masseH
  real, parameter :: T_Cmb = 2.73              ! K

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
  real, parameter :: Lsun = 3.839e26 ! W

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
