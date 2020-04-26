MODULE constant

  ! ----------------------------------------------------------------- !
   ! To avoid duplicates constants, constant calls constantes from
   ! MCFOST and takes what it needs.
  ! ----------------------------------------------------------------- !

  use constantes, only : PI, KM_TO_M, masseH, c_light, hp, kb, cst_th, electron_charge
  use mcfost_env, only : dp

  integer, parameter :: NELEM_WEIGHTS = 99

  ! --- Physical constants ---

  real(kind=dp), parameter ::  CLIGHT=c_light!2.99792458d8 ! Speed of light [m/s]
  real(kind=dp), parameter ::  HPLANCK=hp!6.6260755d-34 !Planck's constant [Js]
  real(kind=dp), parameter ::  KBOLTZMANN=kb!1.380658d-23  !Boltzman's constant [J/K]
  real(kind=dp), parameter ::  AMU= masseH * 1d-3!1.6605402d-27 !Atomic mass unit [kg]
  real(kind=dp), parameter ::  M_ELECTRON=9.1093897d-31 !Electron mass [kg]
  real(kind=dp), parameter ::  Q_ELECTRON=electron_charge!1.60217733d-19  !Electron charge [C]
  real(kind=dp), parameter ::  EPSILON_0=8.854187817d-12  !Vacuum permittivity [F/m]
  !real(kind=dp), parameter ::  MU_0=1.2566370614d-6  !Magnetic induct. of vac.
  real(kind=dp), parameter ::  RBOHR=5.29177349d-11  !Bohr radius [m]
  real(kind=dp), parameter ::  E_RYDBERG=2.1798741d-18 !Ion. pot. Hydrogen [J]
  real(kind=dp), parameter ::  EV = electron_charge!1.60217733d-19 ! One electronVolt [J]
  real(kind=dp), parameter ::  THETA0 =5.03974756d+3!log10(e) * eV/k [K^-1]
  real(kind=dp), parameter ::  ABARH=7.42d-41 !polarizabilty of Hydrogen in [Fm^2]
  real(kind=dp), parameter    ::  pia0squarex2 = PI * 2d0 * RBOHR**2 !constant for collision Cross-sections

 ! --- Unit conversions ---

  real(kind=dp), parameter ::  NM_TO_M =1.0d-9
  real(kind=dp), parameter ::  M_TO_NM =1.0d9
  real(kind=dp), parameter ::  CM_TO_M=1.0d-02
  !real(kind=dp), parameter ::  KM_TO_M =1.0d+03
  real(kind=dp), parameter ::  ERG_TO_JOULE=1.0d-07
  real(kind=dp), parameter ::  JOULE_TO_EV=1/EV!6.241506363094028e+18
  real(kind=dp), parameter ::  JOULE_TO_CM1=5.040963080525957d+22
  real(kind=dp), parameter ::  G_TO_KG=1.0d-03
  real(kind=dp), parameter ::  MICRON_TO_NM=1.0d+03

  
  ! ------- Useful RT constants --------- !
  real(kind=dp), parameter :: sigma_e = 8.0*PI/3.0 * (Q_ELECTRON/(sqrt(4.0*PI*EPSILON_0) *&
                                       (sqrt(M_ELECTRON)*CLIGHT)))**4.d0 !Thomson cross-section

  real(kind=dp), parameter    :: hc = HPLANCK * CLIGHT
  real(kind=dp), parameter    :: fourPI = 4d0*PI
  real(kind=dp), parameter    :: hc_fourPI = hc/fourPI
  real(kind=dp), parameter    :: fourPI_hc = fourPI/hc
  real(kind=dp), parameter    :: twohc = (2.*HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
  real(kind=dp), parameter    :: hc_k = cst_th / NM_TO_M!(HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
  real(kind=dp), parameter    :: fourPI_h = fourPI / HPLANCK
  !Photoionisation Xsection of Hydrogen, at nu0, alpha0 = sigma0_H * g_bg(0) * neff / Z/Z
  !Note an error in Hubeny Mihalas eq. 7.92. unit should be cm2 not cm^-2 !
  real(kind=dp), parameter    :: sigma0_H = (32d0)/(PI*3.*sqrt(3d0)) * EPSILON_0 * &
          (HPLANCK**(3d0)) / (CLIGHT * (M_ELECTRON*Q_ELECTRON)**(2d0)) ! 7.904e-22 m^2

  !here I have a problem if I try to compute sigma0_H_ff using 7.100 of Hubeny Mihalas with SI units value
  !So I Take the cgs result and turn it to SI ...
  !we multiply sigma0_H_ff by nion (m^-3) * ne(m^-3) to have chi in 1d-10 m^5 * m-3 * m^-3 in m^-1
  real(kind=dp), parameter    :: sigma0_H_ff = 3.6923284d8 * 1d-10 ! cm^5 K^1/2 Hz^-3 -> m^5 K^1/2 Hz^3
   !K0 = (Q_ELECTRON**2)/(4.0*PI*EPSILON_0) / sqrt(M_ELECTRON)
   !K0 = (K0**3) * 4./3. * sqrt(2*pi/3./KBOLTZMANN) / HPLANCK / CLIGHT
   !sigma0_H_ff = K0
          
 ! --- Mathematical constants ---
  real(kind=dp), parameter ::  SQRTPI=sqrt(pi)!1.77245385090551


  ! --- 1/(2sqrt(2)), needed for anisotropy of radiation ---

  real(8), parameter ::  TWOSQRTTWO = 0.35355339059327

  !there is an error in LL04
  real(kind=dp), parameter ::  LARMOR = (Q_ELECTRON / (4.0*PI*M_ELECTRON)) !s^-1 / T
  !LamB = nuL/nu0 * lambda0 = nuL/(c/lambda0) * lambda0 = nuL/c * lambda0**2
  !lamD = lambda0 * vbroad/c --> vB = lamB/lamD = nuL/c * lambda0**2 * c / lambda0 / vbroad
  !vB = nuL * lambda0 / vbroad in T^-1


  ! --- Ionization energy Hmin in [J] ---
  !-> we call it affinity
  real(8), parameter ::  E_ION_HMIN = 0.754 * EV
  real(8), parameter, private :: factor1D = 2d0, factor3D = 8d0/PI
  !previously: Vtherm = 2*KBOLTZMANN/AMU and v=sqrt(Vtherm * T / m + xit**2)
  real(8), parameter :: Vtherm = KBOLTZMANN/AMU * factor1D !m^2/s^2/K
  
  real(8), dimension(NELEM_WEIGHTS) :: atomic_weights
  !starts at 1 for H, ends at NELEM_WEIGHTS
  DATA atomic_weights /1.008,4.003,6.939,9.013,10.81,12.01,  &
       14.01,16.0,19.0,20.18,22.990,24.310,26.98,28.09,30.98,  &
       32.07,35.45,39.95,39.1,40.08,44.96,47.9,50.94,52.0,   &
       54.94,55.85,58.94,58.71,63.55,65.37,69.72,72.6,74.92, &
       78.96,79.91,83.8,85.48,87.63,88.91,91.22,92.91,95.95, &
       99.0,101.1,102.9,106.4,107.9,112.4,114.8,118.7,121.8, &
       127.6,126.9,131.3,132.9,137.4,138.9,140.1,140.9,      &
       144.3,147.0,150.4,152.0,157.3,158.9,162.5,164.9,      &
       167.3,168.9,173.0,175.0,178.5,181.0,183.9,186.3,      &
       190.2,192.2,195.1,197.0,200.6,204.4,207.2,209.0,      &
       210.0,211.0,222.0,223.0,226.1,227.1,232.0,231.0,      &
       238.0,237.0,244.0,243.0,247.0,247.0,251.0, 254.0/

  character(len=2), dimension(NELEM_WEIGHTS) :: elemental_ID

  DATA elemental_ID /'H ','He','Li','Be','B ','C ','N ','O ',    &
                 'F ','Ne','Na','Mg','Al','Si','P ','S ','Cl', &
                 'Ar','K ','Ca','Sc','Ti','V','Cr','Mn','Fe',  &
                 'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br', &
                 'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru', &
                 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ', &
                 'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
                 'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
                 'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
                 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac', &
                 'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es'/

END MODULE constant
