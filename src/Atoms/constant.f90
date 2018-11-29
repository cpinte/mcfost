MODULE constant

  ! ----------------------------------------------------------------- !
   ! To avoid duplicates constants, constant calls constantes from
   ! MCFOST and takes what it needs.
  ! ----------------------------------------------------------------- !

  use constantes, only : PI, KM_TO_M! MCFOST

  integer, parameter :: NELEM_WEIGHTS = 99

  ! --- Physical constants ---

  double precision, parameter ::  CLIGHT=2.99792458d8 ! Speed of light [m/s]
  double precision, parameter ::  HPLANCK=6.6260755d-34 !Planck's constant [Js]
  double precision, parameter ::  KBOLTZMANN=1.380658d-23  !Boltzman's constant [J/K]
  double precision, parameter ::  AMU=1.6605402d-27 !Atomic mass unit [kg]
  double precision, parameter ::  M_ELECTRON=9.1093897d-31 !Electron mass [kg]
  double precision, parameter ::  Q_ELECTRON=1.60217733d-19  !Electron charge [C]
  double precision, parameter ::  EPSILON_0=8.854187817d-12  !Vacuum permittivity [F/m]
  !double precision, parameter ::  MU_0=1.2566370614d-6  !Magnetic induct. of vac.
  double precision, parameter ::  RBOHR=5.29177349d-11  !Bohr radius [m]
  double precision, parameter ::  E_RYDBERG=2.1798741d-18 !Ion. pot. Hydrogen [J]
  double precision, parameter ::  EV =1.60217733d-19 ! One electronVolt [J]
  double precision, parameter ::  THETA0 =5.03974756d+3!log10(e) * eV/k [K^-1]
  double precision, parameter ::  ABARH=7.42d-41 !polarizabilty of Hydrogen in [Fm^2]


 ! --- Unit conversions ---

  double precision, parameter ::  NM_TO_M =1.0d-9
  double precision, parameter ::  CM_TO_M=1.0d-02
  !double precision, parameter ::  KM_TO_M =1.0d+03
  double precision, parameter ::  ERG_TO_JOULE=1.0d-07
  double precision, parameter ::  JOULE_TO_EV=1/EV!6.241506363094028e+18
  double precision, parameter ::  JOULE_TO_CM1=5.040963080525957d+22
  double precision, parameter ::  G_TO_KG=1.0d-03
  double precision, parameter ::  MICRON_TO_NM=1.0d+03

 ! --- Mathematical constants ---


  !double precision, parameter ::  PI  =3.14159265358979
  double precision, parameter ::  SQRTPI=1.77245385090551


  ! --- 1/(2sqrt(2)), needed for anisotropy of radiation ---

  real(8), parameter ::  TWOSQRTTWO = 0.35355339059327

  real(8), parameter ::  LARMOR = (Q_ELECTRON / (4.0*PI*M_ELECTRON)) * NM_TO_M


  ! --- Ionization energy Hmin in [J] ---

  real(8), parameter ::  E_ION_HMIN = 0.754 * EV

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
