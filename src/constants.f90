module constantes

  use parametres

  implicit none
  save

  ! Quelques reels utiles
  real(kind=dp), parameter :: pi = 3.141592653589793238462643383279502884197_dp
  real(kind=dp), parameter :: deux_pi = 2.0_dp * pi
  real(kind=dp), parameter :: un_sur_deux_pi = 1.0_dp/deux_pi
  real(kind=dp), parameter :: quatre_pi = 4.0_dp * pi
  real(kind=dp), parameter :: quatre_tiers_pi = 4.0_dp/3.0_dp * pi
  real(kind=dp), parameter :: pi_sur_deux = 0.5_dp * pi
  real(kind=dp), parameter :: un_tiers = 1.0_dp / 3.0_dp

  real(kind=dp), parameter ::  SQRTPI=sqrt(pi)!1.77245385090551
  ! --- 1/(2sqrt(2)), needed for anisotropy of radiation ---
  real(kind=dp), parameter ::  TWOSQRTTWO = 0.35355339059327

  ! Constantes en SI
  real(kind=dp), parameter :: hp = 6.626070040e-34_dp  ! Planck (J.s) CODATA 2014
  real(kind=dp), parameter :: kb = 1.38064852e-23_dp   ! Boltzmann (J.K^-1) CODATA 2014
  real(kind=dp), parameter :: c_light = 299792458._dp  ! vitesse lumiere (m.s^-1) CODATA 2014
  real, parameter :: cst_th=c_light*hp/kb  ! pour calcul de (h c)/(lambda k T)
  real, parameter :: sigma = 5.670367e-8     ! Stefan (en W/(m^2.K^4)) CODATA 2014
  real, parameter :: Ggrav = 6.67428e-11   ! (m^3.s^-2.kg^-1) e-8 en cgs, CODATA 2016, value recommended by IAU 2015 B3
  real(kind=dp), parameter :: electron_charge = 1.6021766208e-19_dp  ! Coulombs CODATA 2014
  real(kind=dp), parameter :: mel = 9.1093897d-31 !Electron mass [kg]
  real(kind=dp), parameter :: radconst = 4.*sigma/c_light       ! radiation constant

  real(kind=dp), parameter :: Na = 6.022140857e23_dp   ! Nombre d'Avogadro CODATA 2014
  real(kind=dp), parameter :: amu = 1.0_dp/Na          ! atomic mass unit [g]
  real(kind=dp), parameter :: masseH = 1.007825032231_dp * amu   ! H atom mass [g] CODATA 2014
  real, parameter :: mu = 2.3                          ! [g]  2.3 following Walker 2004
  real, parameter :: masse_mol_gaz = mu * masseH
  real, parameter :: T_Cmb = 2.7260                    ! K

  real(kind=dp), parameter ::  E_ION_HMIN = 0.754 * electron_charge !Ionization energy (affinity) Hmin in [J]
  real(kind=dp), parameter ::  epsilon_0=8.854187817d-12  !Vacuum permittivity [F/m]
  real(kind=dp), parameter ::  RBOHR=5.29177349d-11  !Bohr radius [m]
  real(kind=dp), parameter ::  E_RYDBERG=2.1798741d-18 !Ion. pot. Hydrogen [J]
  real(kind=dp), parameter ::  THETA0 = 5.03974756d+3!log10(e) * eV/k [K^-1]
  real(kind=dp), parameter ::  ABARH = 7.42d-41 !polarizabilty of Hydrogen in [Fm^2]
  real(kind=dp), parameter ::  pia0squarex2 = deux_pi * RBOHR**2 !constant for collision Cross-sections

  !there is an error in LL04
  real(kind=dp), parameter ::  LARMOR = (electron_charge / (quatre_pi*mel)) !s^-1 / T
  !LamB = nuL/nu0 * lambda0 = nuL/(c/lambda0) * lambda0 = nuL/c * lambda0**2
  !lamD = lambda0 * vbroad/c --> vB = lamB/lamD = nuL/c * lambda0**2 * c / lambda0 / vbroad
  !vB = nuL * lambda0 / vbroad in T^-1

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
  real(kind=dp), parameter ::  NM_TO_M =1.0d-9
  real(kind=dp), parameter ::  M_TO_NM =1.0d9
  !micron to nm is also km_to_m (and nm to micron m_to_km)

  real(kind=dp), parameter :: Rsun = 6.957e8_dp ! IAU 2015 B3 definition, https://arxiv.org/abs/1605.09788
  real, parameter :: Teff_Sun = 5772.0 ! IAU 2015
  real(kind=dp), parameter :: Rsun_to_AU = Rsun/AU_to_m
  real(kind=dp), parameter :: Au_to_Rsun = 1.0_dp/Rsun_to_AU

  real(kind=dp), parameter :: pc_to_AU = 648000_dp/pi ! IAU 2015 B2
  real(kind=dp), parameter :: rad_to_sec = pc_to_AU
  real(kind=dp), parameter :: AU_to_pc = 1.0/pc_to_AU
  real(kind=dp), parameter :: sec_to_rad = AU_to_pc

  ! Energies
  real(kind=dp), parameter :: eV_to_J = electron_charge
  real(kind=dp), parameter :: erg_to_J = 1.0e-7_dp
  real, parameter :: jansky = 1.0e-26 ! W.m^-2.Hz-1 (F_nu en jansky)
  real, parameter :: Lsun = 3.828e26 ! [W] IAU 2015
  real(kind=dp), parameter ::  j_to_ev=1/ev_to_j!6.241506363094028e+18
  real(kind=dp), parameter ::  j_to_cm1=5.040963080525957d+22
  real, parameter :: Tsun = 5777 ! [K] effective tenperature of the sun

  ! Masses
  real(kind=dp), parameter :: g_to_kg = 1.0e-3_dp
  real(kind=dp), parameter :: kg_to_g = 1.0e3_dp

  real(kind=dp), parameter :: GxMsun   = 1.3271244e20_dp ! IAU 2015
  real(kind=dp), parameter :: logg_Sun = log10(GxMsun/Rsun**2 * 100) ! Needs to be in cgs to match stellar atmosphere models
  real(kind=dp), parameter :: GxMearth = 3.986004e14_dp  ! IAU 2015
  real(kind=dp), parameter :: GxMJup   = 1.2668653e17_dp ! IAU 2015

  real(kind=dp), parameter :: Msun_to_kg = GxMsun/Ggrav
  real(kind=dp), parameter :: kg_to_Msun = 1.0/Msun_to_kg

  real(kind=dp), parameter :: Msun_to_g = Msun_to_kg * kg_to_g
  real(kind=dp), parameter :: g_to_Msun = 1.0_dp/Msun_to_g

  ! Time
  real(kind=dp), parameter :: year_to_s = 31557600 ! 365.25 days

  ! Mixed constants
  real(kind=dp), parameter :: kb_on_mH = kb/(amu*g_to_kg)  ! convert amu to kg
  real(kind=dp), parameter :: sigma_e = 8.0*pi/3.0 * (electron_charge/(sqrt(quatre_pi*EPSILON_0) *&
       (sqrt(mel)*c_light)))**4.d0 !Thomson cross-section
  !Photoionisation Xsection of Hydrogen, at nu0, alpha0 = sigma0_H * g_bg(0) * neff / Z/Z
  !Note an error in Hubeny Mihalas eq. 7.92. unit should be cm2 not cm^-2 !
  real(kind=dp), parameter    :: sigma0_H = (32d0)/(3.0*pi*sqrt(3d0)) * EPSILON_0 * &
       (HP**(3d0)) / (c_light * (Mel*electron_charge)**(2d0)) ! 7.904e-22 m^2
  !here I have a problem if I try to compute sigma0_H_ff using 7.100 of Hubeny Mihalas with SI units value
  !So I Take the cgs result and turn it to SI ...
  !we multiply sigma0_H_ff by nion (m^-3) * ne(m^-3) to have chi in 1d-10 m^5 * m-3 * m^-3 in m^-1
  real(kind=dp), parameter    :: sigma0_H_ff = 3.6923284d8 * 1d-10 ! cm^5 K^1/2 Hz^-3 -> m^5 K^1/2 Hz^-3
  !K0 = (Q_ELECTRON**2)/(4.0*PI*EPSILON_0) / sqrt(M_ELECTRON)
  !K0 = (K0**3) * 4./3. * sqrt(2*pi/3./KBOLTZMANN) / HPLANCK / CLIGHT
  !sigma0_H_ff = K0
  real(kind=dp), parameter :: amu_kg = amu * g_to_kg      ! atomic mass unit [kg]
  real(kind=dp), parameter    :: hc = HP * c_light
  real(kind=dp), parameter    :: hc_fourPI = hc/quatre_pi
  real(kind=dp), parameter    :: fourPI_hc = quatre_pi/hc
  real(kind=dp), parameter    :: twohc = (2.*HP * C_LIGHT) / (NM_TO_M)**(3d0)
  real(kind=dp), parameter    :: hc_k = cst_th / NM_TO_M
  real(kind=dp), parameter    :: fourPI_h = quatre_pi / HP
  real(kind=dp), parameter, private :: factor1D = 2.0, factor3D = 8.0/PI
  !previously: Vtherm = 2*KBOLTZMANN/AMU and v=sqrt(Vtherm * T / m + xit**2)
  real(kind=dp), parameter    :: Vtherm = kb/amu_kg * factor1D !m^2/s^2/K
  !1d3 because AMU is in g

  ! Limites de precision numerique
  real, parameter :: tiny_real = tiny(0.0)
  real, parameter :: huge_real = huge(1.0)
  real(kind=dp), parameter :: tiny_dp = tiny(0.0_dp)
  real(kind=dp), parameter :: huge_dp = huge(1.0_dp)

  real, parameter ::  tiny_real_x1e6 =  tiny_real * 1.0e6

  integer, parameter :: huge_integer = huge(1)
  real , parameter :: max_int = real(huge_integer) * (1.0-1.0e-5)

end module constantes
