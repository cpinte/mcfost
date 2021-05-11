MODULE atom_type

  use math, only : w6js
  use mcfost_env, only : dp
  use constant, only : M_ELECTRON, AMU, E_RYDBERG

  IMPLICIT NONE


   integer, parameter :: ATOM_LABEL_WIDTH=20
   integer, parameter :: ATOM_ID_WIDTH=2, MAX_LENGTH=512
!    real(kind=dp), parameter :: BRK=4.0
!    integer, parameter :: MSHELL=5
!    character(len=20), dimension(16) :: col_recipes !16 recipes for collision rates
!    data col_recipes / "OMEGA", "CE", "CI", "CP", "CH0", "CH+", "CH", "CR",&
!                     "AR85-CHP", "AR85-CHH", "AR85-CEA", "BURGESS", "SHULL82",&
!                     "BADNELL", "SUMMERS", "AR85-CDI" /

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE ZeemanType
   integer :: Ncomponent
   integer, allocatable, dimension(:)  :: q
   real(kind=dp), allocatable, dimension(:)  :: shift, strength
  END TYPE ZeemanType

  !Sets at reading: First LINES THEN CONTINUA
  !line, 1->Nl, continuum, Nl+1->Nc+Nl
  TYPE AtomicTransition
   character(len=16) :: trtype
   integer :: ik
   logical :: lcontrib_to_opac =.true.
  END TYPE AtomicTransition

  TYPE AtomicLine
   logical           :: symmetric, polarizable
   !logical, dimension(:), allocatable :: negative_opacity
   logical           :: Voigt=.true., write_flux_map=.false., PFR=.false.,&
      damping_initialized=.false. !true if we store the damping on the whole grid for all lines.
   character(len=17) :: vdWaals
   character(len=20) :: trtype="ATOMIC_LINE", Coupling="LS"
   ! i, j start at 1 (not 0 like in C)
   integer :: i, j, Nxrd=0
   integer :: Nlambda ! Number of wavelength points between Nblue and Nred = Nred-Nblue + 1
   integer :: Nblue=0, Nred = 0, Nmid=0 !Index of the line on the whole grid in the absence of velocity shifts or magnetic fields
   real(kind=dp) :: lambda0, lambdamin, lambdamax ! boundary of the line in absence of velocity shifts or magnetic fields
   real(kind=dp) :: isotope_frac, g_Lande_eff, Aji, Bji, Bij, Grad, cStark, fosc
   real(kind=dp) :: twohnu3_c2, gij
   !qcore is not used yet: size of the line "core" if non uniform grid (for interp)
   !qwing is used to compute the maximum extension of the line: Vmax = qwing * max(vD)
   real :: qcore, qwing
   real(kind=dp) :: glande_i, glande_j
   real(kind=dp), dimension(4) :: cvdWaals
   !!integer, allocatable, dimension(:,:) :: dk!size (Nray, id) = index displacement of a line due to velocity
   !Nlambda,(Nproc or Nspace, depends on the choice For flux calculations. Always Nlambda, Nspace for NLTE
   real(kind=dp), allocatable, dimension(:,:)  :: phi
   !!real(kind=dp), allocatable, dimension(:,:,:) :: eta
   real(kind=dp), allocatable, dimension(:,:,:,:) :: mapc!lcont flux to be stored
   real(kind=dp), allocatable, dimension(:,:,:,:,:) :: map!line flux to be stored
   !used for wavelength integration
   real(kind=dp), allocatable, dimension(:,:,:) :: phi_loc!, phiZ, psi !3, Nlambda, Nray
   !wlam is the integration wavelenght weight = phi
   real(kind=dp), allocatable, dimension(:)  :: lambda, CoolRates_ij, w_lam, Rij, Rji, Jbar
   !real(kind=dp), allocatable, dimension(:) :: fomega !for Rayleigh scattering
   real(kind=dp) :: Qelast, adamp!, chi0, eta0 ! at a cell
   real(kind=dp), dimension(:), allocatable :: Tex, deta !differentiel of source function
   !keep CLIGHT * (nu0 - nu)/nu0 for lines												(method for Voigt)
   real(kind=dp), dimension(:), allocatable :: u, a, pvoigt_eta!aeff, r, r1
   !damping for all cells(a) and Thomson effective damping (eff)
   real(kind=dp), allocatable, dimension(:,:) :: rho_pfr
   character(len=ATOM_LABEL_WIDTH) :: name ! for instance Halpha, h, k, Hbeta, D1, D2 etc
   integer :: ZeemanPattern! 0 = WF, -1=EFFECTIVEE TRIPLET, 1=FULL
   type (AtomType), pointer :: atom => NULL()
   type (ZeemanType) :: zm
  END TYPE AtomicLine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE AtomicContinuum
   logical :: hydrogenic
   !logical, dimension(:), allocatable :: negative_opacity
   integer :: i, j, Nlambda, Nblue = 0, Nred = 0, Nmid = 0, N0 = 0, Nb, Nr
   real(kind=dp) :: lambda0, isotope_Frac, alpha0, lambdamin, lambdamax !continuum maximum frequency > frequency photoionisation
   real(kind=dp), allocatable, dimension(:)  :: lambda, alpha, twohnu3_c2, CoolRates_ij, w_lam, ehnukt
   real(kind=dp), allocatable, dimension(:)  :: lambda_file, alpha_file
   real(kind=dp), dimension(:), allocatable :: Tex
   character(len=ATOM_LABEL_WIDTH) :: name !read in the atomic file
   type (AtomType), pointer :: atom => NULL()
   real(kind=dp), allocatable, dimension(:) :: Rji, Rij
   real(kind=dp), allocatable, dimension(:,:)  :: gij, Jnu, alpha_nu
   character(len=20) :: trtype="ATOMIC_CONTINUUM"
  END TYPE AtomicContinuum

  TYPE AtomType
   character(len=ATOM_ID_WIDTH) :: ID
   character(len=28) :: inputFile
   character(len=ATOM_LABEL_WIDTH), allocatable, dimension(:)  :: label
   logical                :: NLTEpops, active, set_ltepops
   ! atom can be passive but NLTEpops true. This is the case of
   ! populations read from previous run
   character(len=15)      :: initial_solution
   integer                :: Nlevel, Nline, Ncont, Ntr, Npfr=0, Ntr_line
   integer                :: periodic_table, activeindex !order of the active atom in the
   														! active atoms array ActiveAtoms
   ! BY CONVENTION, stage=0 for neutrals, 1 for singly ionised
   ! ions etc ...
   integer, allocatable, dimension(:)  :: stage, Lorbit
   !integer(8)            :: offset_coll, colunit
   real(kind=dp) :: Rydberg, scatt_limit !minimum wavelength for Rayleigh scattering
   real(kind=dp)                :: cswitch = 1.0_dp, Abund, weight, massf !mass fraction
   real(kind=dp), allocatable, dimension(:) :: g, E, vbroad!, ntotal
   real(kind=dp), allocatable, dimension(:) :: qS, qJ
   ! allocated in readatom.f90, freed with freeAtoms()
   !futur deprecation, because I will stop using RH routine
   character(len=MAX_LENGTH), allocatable, dimension(:) :: collision_lines !to keep all remaning lines in atomic file
   !Nlevel * Nlevel * Nproc
   real(kind=dp), dimension(:,:,:), allocatable :: Gamma, C!, Rij
   real(kind=dp), dimension(:,:), pointer :: n, nstar
   real(kind=dp), dimension(:,:), allocatable :: b !not a pointer this one
   ! arrays of lines, continua containing different line, continuum each
   type (AtomicLine), allocatable, dimension(:)         :: lines
   type (AtomicContinuum) , allocatable, dimension(:)   :: continua
   type (AtomicTransition), allocatable, dimension(:)   :: at !Atomic transition, lines first in readatom
   !one emissivity per atom, used in the construction of the gamma matrix
   !where I have to distinguish between atom own opac and overlapping transitions
   real(kind=dp), allocatable, dimension(:,:) :: etac, chic
   real(kind=dp), allocatable, dimension(:,:,:) :: eta !Nwaves, Nrays, Nproc
  END TYPE AtomType
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE Element
   !Element is a structure similar to atom but for all the Nelem
   ! of the periodic table. Even if no atomic model (store in AtomType)
   ! exists for an Element, this Element is used for electron density
   ! calculations, overall abundances, eventually chemical equilibrium etc.
   ! Further, pure-LTE calculations can be performed with this Element
   ! if desired.
   character(len=ATOM_ID_WIDTH) :: ID
   logical :: abundance_set
   integer :: Nstage, Nmolecule ! umber of Molecules having an element
   integer, allocatable, dimension(:)  :: mol_index !track molecules in which Element is present
   real(kind=dp) :: weight, abund, massf
   real(kind=dp), allocatable, dimension(:)  :: ionpot
   real(kind=dp), allocatable, dimension(:,:)  :: pf!, n !LTE populations, not used nor allocated anymore
   !n is the population for each stage at a given grid point
   type (AtomType), pointer :: model => NULL()
 END TYPE Element

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 CONTAINS
 
 FUNCTION rydberg_atom(atom)!correct from electron mass
  type (AtomType), intent(in) :: atom
  real(kind=dp) :: rydberg_atom, deltam
  
   deltam = 1. + M_ELECTRON / (atom%weight * AMU)
 
   rydberg_atom = E_RYDBERG / deltam !in joules
 RETURN 
 END FUNCTION rydberg_atom
 
 !only for Hydrogen or approximate for n = n_eff
 function line_oscillator_strength(n, np)
 !From Johnson 1972
  real(kind=dp), intent(in) :: n, np !n_lower and n_upper
  real(kind=dp) :: line_oscillator_strength
  real(kind=dp) :: x, g
  
  x = 1.0 - (n/np)**2
  
  if (n < 2.0) then !1 if n integer
   g = 1.1330 -0.4059 / x + 0.07014 / x / x
  else if (n >= 2.0 .and. n < 3.0) then ! 2 if n integer
   g = 1.0785 - 0.2319 / x + 0.02947 / x / x
  else !n>=3
   g = 0.9935 + 0.2328 / n - 0.1296 / n / n - &
    1./n * (0.6282 - 0.5598 / n + 0.5299 / n / n) / x + &
    1./n/n * (0.3887 - 1.181 /n + 1.470 /n / n) / x / x
  endif
  
  !32/3/sqrt(3.)/pi
  line_oscillator_strength = 1.9603 * n * g / (np * x)**3.
 
 return
 end function line_oscillator_strength
 
 function n_eff(Erydb, Ej, Ei, Z)
  !effective (hydrogenic) quantum number n*
  !Ei is the energy of the level for which we compute n*
  !and Ej is the energy of the next continuum.
  !In case of a line beware, Ej is not E(line%j) but E(next_Contiuum(line%i))
  real(kind=dp) :: n_eff
  real(kind=dp), intent(in) :: Erydb, Ej, Ei
  integer, intent(in) :: Z
    
  n_eff = Z * sqrt(Erydb / (Ej-Ei) )
 
 return
 end function n_eff


 !for lines only, for continuum it is simply cont%j
 FUNCTION find_continuum(atom, l)
  integer :: l, find_continuum, l0
  type (AtomType) :: atom
  
  l0 = l

  find_continuum = l + 1
  do while ((atom%stage(find_continuum) < atom%stage(l)+1).and.(find_continuum <= atom%Nlevel))
   find_continuum = find_continuum + 1
  end do

  if (atom%stage(find_continuum) == atom%stage(l)+1) then
     return !continuum level reach
 else
     write(*,*) "Coundn't find continuum level for level ", l0, find_continuum, atom%Nlevel
     find_continuum = 0
     write(*,*) atom%ID, atom%stage(l), l0
     return
 end if

 RETURN
 END FUNCTION find_continuum
 
 FUNCTION find_ground_state_ion(Nlevel, stage, Nstage)
  !Search the index of each ground state of each ionic stage
  integer, intent(in) :: Nlevel, Nstage
  integer, intent(in), dimension(Nlevel) :: stage
  integer :: find_ground_state_ion(Nstage)
  integer :: k, j, m
  
  j = stage(1) !first stage to start, then the next index is necessarily at least j+1
  m = 1
  do k=1, Nlevel
  
   if (stage(k)==j) then
    !write(*,*) "k=", k, j, m
    find_ground_state_ion(m) = k
    j = j + 1
    m = m + 1 !index should start at 1, not at j if j > 1
   endif
  
  enddo
 
 
 RETURN
 END FUNCTION find_ground_state_ion

 FUNCTION atomic_orbital_radius(n, l, Z)
 !return atomic orbital radius wrt the Bohr radius
  real(kind=dp) :: atomic_orbital_radius
  real(kind=dp) :: n ! quantum principal number
  integer :: l ! orbital quantum number
  integer :: Z

  atomic_orbital_radius = n*n/real(Z) * (1d0 + 0.5*(1d0 - real(l)*(real(l)+1)/n/n))

  RETURN
 END FUNCTION atomic_orbital_radius

 FUNCTION atomic_orbital_sqradius(n, l, Z)
 !Bates-Damguard mean square radius
  real(kind=dp) :: atomic_orbital_sqradius
  !in a0**2 units
  real(kind=dp) :: n ! quantum principal number
  integer :: l ! orbital quantum number
  integer :: Z ! charge : Z = 1 for H I

  atomic_orbital_sqradius = 0.5*n*n/real(Z*Z) * (5.*n*n + 1 - 3.*real(l)*(real(l)+1))


  RETURN
 END FUNCTION atomic_orbital_sqradius


 FUNCTION getOrbital(orbit, unknown_orbital) result (L)
  integer :: L
  character(len=1), intent(in) :: orbit
  logical, intent(out) :: unknown_orbital
  
  unknown_orbital = .false. 
  
  L = -1
  
  SELECT CASE (orbit)
   CASE ('S')
    L = 0
   CASE ('P')
    L = 1
   CASE ('D')
    L = 2
   CASE ('F')
    L = 3
   CASE ('G')
    L = 4
   CASE ('H')
    L = 5
   CASE ('I')
    L = 6
   CASE ('J')
    L = 7
   CASE ('K')
    L = 8
   CASE ('L')
    L = 9
   CASE ('M')
    L = 10
   CASE ('N')
    L = 11
   CASE ('O')
    L = 12
   CASE ('Q')
    L = 13
   CASE ('R')
    L = 14
   CASE ('T')
    L = 15
   CASE ('U')
    L = 16
   CASE ('V')
    L = 17
   CASE ('W')
    L = 18
   CASE ('X')
    L = 19
   CASE DEFAULT
    write(*,'("Orbit " (1A1) "unknown")') orbit
    unknown_orbital = .true.
   END SELECT

  END FUNCTION getOrbital


  !Not working correctly, using effective quantum Number instead. n**2 = Z**2 * Rydberg/deltaE
  FUNCTION getPrincipal(label, n) result(determined)
  ! Read principal quantum number of H I
  ! (should work on other elements)
  ! max len of quantum number n is 2!
   logical :: determined
   character(len=ATOM_LABEL_WIDTH) :: label
   real(kind=dp), intent(out) :: n
   integer :: err, st
   determined=.false.
   !           n principal
   !       12345
   !label='H I 1S 2SE          '
   !       123456
   !label='H I 10S 2SE         '


   if (label(2:2)==" ") then !"H X ..." case
    if (label(4:4)=="I") then !case H II
     st = 6
    else
     st = 5
    end if !always at this position for HI because
          !label of HI is always 'H I NX.....'
          !or 'H I NNX....'!1234567     123456
   else !other than H e.g., CA II ...., CA I ...
    if (label(6:6)=="I") then
     write(*,*) "(1) Case not handled for parsing principal"
     determined = .false.
     RETURN
    else if (label(4:4)=="I") then
     if (label(5:5)=="I") then
      st = 7
     else if (label(5:5)==" ") then
       st = 7
     else
      write(*,*) "(2) Error parsing label for principal quantum number"
      determined = .false.
      RETURN
     end if
    else
      write(*,*) "(3) Error parsing label for principal quantum number"
      determined = .false.
      RETURN
    end if
   end if

   read(label(st:st+1),*,iostat=err) n
   ! if error, label(st+1) is a string,
   ! meaning that len(n) = 1 (one digit)
   if (err.gt.0) read( label(st:st),*) n
   ! but if err is 0, len(n)=2, label(st+1) is integer
   ! two digits n
   determined = .true.

  RETURN
  END FUNCTION getPrincipal

  !Should be more efficient ?
  ! easy to break
  SUBROUTINE determinate(label, g, S, L,J, determined)
  ! get principal quantum number from label
   logical, intent(out) :: determined
   integer, intent(out) :: L!, J
   logical :: unknown_orbital
   real(kind=dp) :: J
   real(kind=dp), intent(out) :: S
   real(kind=dp), intent(in) :: g
   character(len=ATOM_LABEL_WIDTH), intent(in) :: label
   character(len=ATOM_LABEL_WIDTH+1) :: multiplet
   character(len=1) :: orbit
   integer :: i, Nl, parity, dk
   multiplet = trim(label)
   Nl = len(multiplet)

   dk = -1
   i = Nl-1
   !Parse label to get Multiplicity = 2*S + 1 and Orbital
   !'NA I 2P6 3S 2SE     '
   !'NA I 2P6 3P 2PO 3   '
   !'H I 1S 2SE          '
   !2S+1L is after E or O from decreasing order. Parity is +2 after E,O
  do while (i.ne.1)
   if (multiplet(i:i).eq."O" .or. multiplet(i:i).eq."E") then
    if (multiplet(i+2:i+2).eq."1" .or. multiplet(i+2:i+2).eq."0" &
   .or. multiplet(i+2:i+2).eq."2" .or. multiplet(i+2:i+2).eq."3" .or. multiplet(i+2:i+2).eq."3") then
     read(multiplet(i+2:i+2),*) parity
     !write(*,'("Parity: " (1I1))') parity
    end if
    read(multiplet(i-2:i-2),*) S
    orbit = multiplet(i-1:i-1)
    i = 1
    determined = .true.
   else
    i = i + dk
   end if
  end do

  if (i.gt.1) then
   determined = .false.
   RETURN
  end if

  S = (real(S,kind=8) - 1.)/2.
  L = getOrbital(orbit, unknown_orbital)
  if (unknown_orbital) then
   determined = .false.
   return
  endif
   
  J = L + S!(g-1.)/2.
  if (J>(L+S)) then
   determined = .false.
  end if

  RETURN
  END SUBROUTINE determinate

  !LL04 eq. 10.12, p. 514
  FUNCTION wKul (atom, kr, K) result(wK)
   type (AtomType) :: atom
   integer :: kr, i, j
   integer :: K, twice
   real wK, Jj, Ji

   wK = 1.
   twice = 2 !
   if (K.eq.0) RETURN

   i = atom%lines(kr)%i
   j = atom%lines(kr)%j

   Jj = twice*atom%qJ(j)
   Ji = twice*atom%qJ(i)

   wK = w6js(twice, twice, twice*K, int(Jj), int(Jj), &
         int(Ji))/ w6js(twice, twice, 0, int(Jj), &
         int(Jj), int(Ji))

   !WK = wKul**2
  RETURN
  END FUNCTION wKul

   FUNCTION atomZnumber(atom)
   ! --------------------------------------
   ! return the atomic number of an atom
   ! with ID = atomID.
   ! Hydrogen is 1
   ! --------------------------------------
    type (AtomType), intent(in) :: atom
    integer :: atomZnumber

    atomZnumber = atom%periodic_table

   RETURN
   END FUNCTION atomZnumber

   SUBROUTINE PTrowcol(Z, row, col)
   ! ---------------------------------------------------
   ! Position in the periodic table of an atom according
   ! to this Z number.
   ! min row is 1, max row is 87 Francium
   ! No special ceses yet for Lanthanides and Actinides
   ! min col is 1, max col is 18
   ! beware He, is on the 18 column, because they are 16
   ! empty columns between H and He
   ! ---------------------------------------------------

   integer i
   integer, intent(in) :: Z
   integer, intent(out) :: row, col
   integer :: istart(7)

   row = 0
   col = 0

   istart = (/1, 3, 11, 19, 37, 55, 87/) !first Z of each rows

   ! -- find row
   do i=1,6
    if ((Z .ge. istart(i)) .and. (Z .lt. istart(i + 1))) then
     ! we are on the row istart(i)
     row=i
     exit
    end if
   end do
   ! find column on the row
   col = Z - istart(i) + 1

   ! for relative position just comment out the following lines
   ! or take col=col-10 (or-18 for He) for these cases.

   ! special case of Z=2 for Helium because there are 16
   ! empty columns between Hydrogen and Helium
   if (Z.eq.2) then
    col = col + 16
   ! ten empty lines between Be and B
   ! ten empty lines between Mg and Al
   else if (((istart(i).eq.3).or.(istart(i)).eq.11).and.&
      (Z.gt.istart(i)+1)) then
    !row = istart(i)
    col = col+10
   end if
  RETURN
  END SUBROUTINE PTrowcol

END MODULE atom_type
