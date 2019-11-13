MODULE atom_type

  use math, only : w6js
  use accelerate, only : Ng
  use mcfost_env, only : dp

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
   logical           :: symmetric, polarizable!!, lcontrib_to_opac !default is yes, set at reading
   logical           :: Voigt=.true., PFR=.false.,&
      damping_initialized=.false. !true if we store the damping on the whole grid for all lines.
   character(len=17) :: vdWaals
   character(len=20) :: trtype="ATOMIC_LINE", Coupling="LS"
   ! i, j start at 1 (not 0 like in C)
   integer :: i, j, Nlambda, Nblue=0, Nxrd=0, Nred = 0, Nmid=0
   real(kind=dp) :: lambda0, isotope_frac, g_Lande_eff, Aji, Bji, Bij, Grad, cStark, fosc
   real(kind=dp) :: twohnu3_c2, gij
   real(kind=dp) :: qcore, qwing, glande_i, glande_j
   real(kind=dp), dimension(4) :: cvdWaals
   !Nlambda,Nproc
   !for one cell
   real(kind=dp), allocatable, dimension(:,:)  :: phi, phi_loc !in a specific direction
   !used for wavelength integration
   real(kind=dp), allocatable, dimension(:,:,:) :: phi_ray !for one cell Nlambda, Nray, Nproc
   real(kind=dp), allocatable, dimension(:,:,:) :: phiZ, psi !3, Nlambda, Nray
   !wlam is the integration wavelenght weight = phi
   real(kind=dp), allocatable, dimension(:)  :: lambda, CoolRates_ij, w_lam
   real(kind=dp) :: wphi
   real(kind=dp) :: Qelast, Rij, Rji, adamp ! at a cell
   !keep CLIGHT * (nu0 - nu)/nu0 for lines
   real(kind=dp), dimension(:), allocatable :: u, a !damping for all cells, 
   real(kind=dp), allocatable, dimension(:,:) :: rho_pfr
   !!Nlevel, wavelength and proc
   !!Stores the information for that atom only, necessary to  construct the Gamma matrix
   !!and to compute the cross-coupling terms. We need to know for each wavelength and each
   !!proc what are the active transitions involved in the Gamma matrix.
   !!Nlambda, Ndep
   !!real(kind=dp), allocatable, dimension(:,:,:)  :: U, V
   character(len=ATOM_LABEL_WIDTH) :: name ! for instance Halpha, h, k, Hbeta, D1, D2 etc
   integer :: ZeemanPattern! 0 = WF, -1=EFFECTIVEE TRIPLET, 1=FULL
   type (AtomType), pointer :: atom => NULL()
   type (ZeemanType) :: zm
  END TYPE AtomicLine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE AtomicContinuum
   logical :: hydrogenic!!, lcontrib_to_opac
   integer :: i, j, Nlambda, Nblue = 0, Nred = 0, Nmid = 0
   real(kind=dp) :: lambda0, isotope_Frac, alpha0, lambdamin !continuum maximum frequency > frequency photoionisation
   real(kind=dp), allocatable, dimension(:)  :: lambda, alpha, twohnu3_c2, CoolRates_ij, w_lam
   real(kind=dp), allocatable, dimension(:)  :: lambda_file, alpha_file
   real(kind=dp) :: wmu
   real(kind=dp) :: Rji, Rij
   character(len=ATOM_LABEL_WIDTH) :: name !read in the atomic file
   type (AtomType), pointer :: atom => NULL()
   !!Nlambda, Ndep
   real(kind=dp), allocatable, dimension(:,:)  :: gij
   character(len=20) :: trtype="ATOMIC_CONTINUUM"
  END TYPE AtomicContinuum

  TYPE AtomType
   character(len=ATOM_ID_WIDTH) :: ID
   character (len=15)             :: dataFile!popsinFile, popsoutFile
   character(len=28) :: inputFile
   character(len=ATOM_LABEL_WIDTH), allocatable, dimension(:)  :: label
   logical                :: NLTEpops, active
   ! atom can be passive but NLTEpops true. This is the case of
   ! populations read from previous run
   character(len=15)      :: initial_solution
   integer                :: Nlevel, Nline, Ncont, Ntr, Npfr=0, Ntr_line
   integer                :: periodic_table, activeindex !order of the active atom in the
   														! active atoms array ActiveAtoms
   ! BY CONVENTION, stage=0 for neutrals, 1 for singly ionised
   ! ions etc ...
   integer, allocatable, dimension(:)  :: stage, Lorbit
   integer(8)            :: offset_coll, colunit
   real(kind=dp)                :: Abund, weight
   real(kind=dp), allocatable, dimension(:) :: g, E, vbroad!, ntotal
   real(kind=dp), allocatable, dimension(:) :: qS, qJ
   ! allocated in readatom.f90, freed with freeAtoms()
   
   !futur deprecation, because I will stop using RH routine
   character(len=MAX_LENGTH), allocatable, dimension(:) :: collision_lines !to keep all remaning lines in atomic file
   !Nlevel * Nlevel * Nproc
   real(kind=dp), dimension(:,:,:), allocatable :: Gamma, C
   real(kind=dp), dimension(:,:), pointer :: n, nstar
   ! arrays of lines, continua containing different line, continuum each
   type (AtomicLine), allocatable, dimension(:)         :: lines
   type (AtomicContinuum) , allocatable, dimension(:)   :: continua
   type (AtomicTransition), allocatable, dimension(:)   :: at !Atomic transition, lines first in readatom
   !one emissivity per atom, used in the construction of the gamma matrix
   !where I have to distinguish between atom own opac and overlapping transitions
   real(kind=dp), allocatable, dimension(:,:,:) :: eta !Nwaves, Nproc
!	_down = from j (upper) to l (lower); _up from i (lower) to lp (upper)
   real(kind=dp), allocatable, dimension(:,:,:,:) :: Uji_down, chi_up, chi_down
   type (Ng) :: Ngs
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
   real(kind=dp) :: weight, abund
   real(kind=dp), allocatable, dimension(:)  :: ionpot
   real(kind=dp), allocatable, dimension(:,:)  :: pf, n !LTE populations, not used nor allocated anymore
   !n is the population for each stage at a given grid point
   type (AtomType), pointer :: model => NULL()
 END TYPE Element

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 CONTAINS
 
 !for lines only, for continuum it is simply cont%j
 FUNCTION find_continuum(atom, l)
  integer :: l, find_continuum
  type (AtomType) :: atom
  
  find_continuum = l + 1
  do while ((atom%stage(find_continuum) < atom%stage(l)+1).and.(find_continuum <= atom%Nlevel))
   find_continuum = find_continuum + 1
  end do
  
  if (atom%stage(find_continuum) == atom%stage(l)+1) then
     return !continuum level reach
 else
     write(*,*) "Coundn't find continuum level for level ", l, find_continuum, atom%Nlevel
     find_continuum = 0
     return
 end if
  
 RETURN
 END FUNCTION find_continuum
 

 FUNCTION atomic_orbital_radius(n, l, Z)
 !return atomic orbital radius wrt the Bohr radius
  real(kind=dp) :: atomic_orbital_radius
  real(kind=dp) :: n ! quantum principal number
  integer :: l ! orbital quantum number
  integer :: Z ! charge ?
  
  atomic_orbital_radius = n*n/real(Z) * (1d0 + 0.5*(1d0 - real(l)*(real(l)+1)/n/n))
  
  RETURN
 END FUNCTION atomic_orbital_radius
 
 FUNCTION atomic_orbital_sqradius(n, l, Z)
 !Bates-Damguard mean suare radius 
  real(kind=dp) :: atomic_orbital_sqradius
  !in a0**2 units
  real(kind=dp) :: n ! quantum principal number
  integer :: l ! orbital quantum number
  integer :: Z ! charge ?
  
  atomic_orbital_sqradius = 0.5*n*n/real(Z*Z) * (5.*n*n + 1 - 3.*real(l)*(real(l)+1))
  
  
  RETURN
 END FUNCTION atomic_orbital_sqradius
 

 FUNCTION getOrbital(orbit) result (L)
  integer :: L
  character(len=1), intent(in) :: orbit
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
    write(*,*) "exiting..."
    stop
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
   integer, intent(out) :: L
   real(kind=dp), intent(out) :: S, J
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
   .or. multiplet(i+2:i+2).eq."2" .or. multiplet(i+2:i+2).eq."3") then
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
  L = getOrbital(orbit)
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
