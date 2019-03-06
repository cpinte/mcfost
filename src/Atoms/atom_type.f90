MODULE atom_type

  use math, only : w6js

  IMPLICIT NONE


   integer, parameter :: ATOM_LABEL_WIDTH=20
   integer, parameter :: ATOM_ID_WIDTH=2



  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE AtomicLine
   !type (AtomType), pointer   :: atom
   logical           :: symmetric, polarizable
   logical           :: Voigt=.true., PFR=.false.,&
      damping_initialized=.false. !true if we store the damping on the whole grid for all lines.
   character(len=17) :: vdWaals
   character(len=20) :: trtype="ATOMIC_LINE"
   ! i, j start at 1 (not 0 like in C)
   integer :: i, j, Nlambda, Nblue=0, Nxrd=0, Nred = 0, Nmid=0
   real(8) :: lambda0, isotope_frac, g_Lande_eff, Aji, Bji, Bij, Grad, cStark, fosc
   real(8) :: qcore, qwing
   real(8), dimension(4) :: cvdWaals
   real(8), allocatable, dimension(:,:)  :: phi, phi_Q, phi_U, phi_V, psi_Q, psi_U, psi_V
   double precision, allocatable, dimension(:)  :: lambda, CoolRates_ij, wphi!, c_shift, c_fraction
   double precision :: Qelast, adamp, Rij, Rji
   real(8), allocatable, dimension(:,:) :: rho_pfr
   !type (AtomType), pointer :: atom
  END TYPE AtomicLine
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE AtomicContinuum
   logical :: hydrogenic
   integer :: i, j, Nlambda, Nblue = 0, Nred = 0, Nmid = 0
   real(8) :: lambda0, isotope_Frac, alpha0
   real(8), allocatable, dimension(:)  :: lambda, alpha, CoolRates_ij
   double precision :: Rji, Rij
   !type (AtomType), pointer :: atom
   character(len=20) :: trtype="ATOMIC_CONTINUUM"
  END TYPE AtomicContinuum
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !TYPE FixedTransition
  ! fixed transitions are usefull to describe NLTE problem
  ! in the solar chromosphere in 2D,3D without not so much
  ! computational cost.
  ! They are not implemented because only relevent for solar
  ! application
  !END TYPE FixedTransition
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  TYPE AtomType
   character(len=ATOM_ID_WIDTH) :: ID
   character (len=15)             :: dataFile!popsinFile, popsoutFile
   character(len=28) :: inputFile
   character(len=ATOM_LABEL_WIDTH), allocatable, dimension(:)  :: label
   logical                :: NLTEpops, active
   ! atom can be passive but NLTEpops true. This is the case of
   ! populations read from previous run
   character(len=15)      :: initial_solution
   integer                :: Nlevel, Nline, Ncont, Nfixed, Npfr=0
   integer                :: periodic_table, activeindex
   ! BY CONVENTION, stage=0 for neutrals, 1 for singly ionised
   ! ions etc ...
   integer, allocatable, dimension(:)  :: stage, Lorbit
   integer(8)            :: offset_coll, colunit
   real(8)                :: Abund, weight
   real(8), allocatable, dimension(:) :: g, E, vbroad, ntotal
   real(8), allocatable, dimension(:) :: qS, qJ
   ! allocated in readatom.f90, freed with freeAtoms()
   real(8), dimension(:), allocatable :: C, Gamma !now depth dependence is dropped
   real(8), dimension(:,:), pointer :: n, nstar
   ! arrays of lines, continua containing different line, continuum each
   type (AtomicLine), allocatable, dimension(:)         :: lines
   type (AtomicContinuum) , allocatable, dimension(:)   :: continua
   !type (FixedTransition) :: fts
   ! normally in rhthreads
   real, allocatable, dimension(:)  :: gij, Vij, wla, chi_up, chi_down, Uji_down, eta
  END TYPE AtomType
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !building
  TYPE AtomicTransition
   ! a type storing all transitions from atom
   !independently of their type
   character(len=20) :: trtype !transition type
   !UNION transitions
    !type (AtomicContinuum) :: line
    !type (AtomicLine) :: continuum
   !END UNION transitions
  END TYPE AtomicTransition
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
   integer :: Nstage, Nmolecule
   integer, allocatable, dimension(:)  :: mol_index
   real(8) :: weight, abund
   real(8), allocatable, dimension(:)  :: ionpot
   real(8), allocatable, dimension(:,:)  :: pf, n
   !n is the population for each stage at a given grid point
   type (AtomType), pointer :: model => NULL()
 END TYPE Element

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 TYPE ZeemanType
  integer :: Ncomponent
  integer, allocatable, dimension(:)  :: q
  double precision, allocatable, dimension(:)  :: shift, strength
 END TYPE ZeemanType

 CONTAINS

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

  FUNCTION getPrincipal(label, n) result(determined)
  ! Read principal quantum number of H I
  ! (should work on other elements)
  ! max len of quantum number n is 2!
   logical :: determined
   character(len=ATOM_LABEL_WIDTH) :: label
   real(8), intent(out) :: n
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
   real(8), intent(out) :: S, J
   real(8), intent(in) :: g
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
  J = (g-1.)/2.
  if (J.gt.(L+S)) then
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

END MODULE atom_type
