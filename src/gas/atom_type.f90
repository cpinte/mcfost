module atom_type

   use mcfost_env, only : dp
   use constantes, only   : mel, AMU_KG, E_RYDBERG, vtherm
   use messages, only : error

   implicit none

   integer, parameter   :: ATOM_LABEL_WIDTH=20
   integer, parameter   :: ATOM_ID_WIDTH=2
   integer, parameter   :: Nmax_Trans_raytracing = 100

   type ZeemanType
      integer :: Ncomponent
      integer, allocatable, dimension(:)  :: q
      real(kind=dp), allocatable, dimension(:)  :: shift, strength
   end type ZeemanType


   type AtomicContinuum
      logical :: hydrogenic, lcontrib
      integer :: i, j
      integer :: Nb, Nr, Nlambda
      real(kind=dp) :: lambda0, isotope_Frac, alpha0, lambdamin, lambdamax !continuum maximum frequency > frequency photoionisation
      real(kind=dp), allocatable, dimension(:)  :: alpha, twohnu3_c2
      real(kind=dp), allocatable, dimension(:)  :: lambda_file, alpha_file
      type (AtomType), pointer :: atom => NULL()
   end type AtomicContinuum

   type AtomicLine
      logical  :: polarizable
      logical  :: Voigt, lcontrib
      character(len=17) :: vdWaals
      integer :: i, j
      integer :: Nb, Nr, Nlambda
      integer :: Nover_inf, Nover_sup !max index of overlap with another line
      real(kind=dp) :: lambda0, lambdamin, lambdamax
      real(kind=dp) :: Aji, Bji, Bij, Grad, cStark, fosc
      real(kind=dp) :: twohnu3_c2, gij, vmax !m/s
      real :: qwing
      real(kind=dp), dimension(4) :: cvdWaals
      real(kind=dp), dimension(:), allocatable :: a
      real(kind=dp), dimension(:,:,:,:,:), allocatable :: map !2d flux in the line
      ! real(kind=dp), dimension(:), allocatable :: u, pvoigt_eta, aeff, r, r1 !thomson approx. of Voigt.
      ! real(kind=dp), allocatable, dimension(:,:)  :: phi
      type (AtomType), pointer :: atom => NULL()
      integer :: ZeemanPattern
      real :: glande_i, glande_j, g_Lande_eff
      type (ZeemanType) :: zm
   end type AtomicLine

   type AtomType
      character(len=ATOM_ID_WIDTH) :: ID
      character(len=28) :: filename
      logical :: active, lline !non-lte ?, images ?
      integer :: initial, nTrans_raytracing
      !initial:
      !0->LTE; 1->OLD_POPULATIONS; 2->ZERO_RADIATION; 3->CSWITCH; 4->SOBOLEV/CEP 
      !only for lines !
      integer :: j_Trans_rayTracing(Nmax_Trans_raytracing), i_Trans_rayTracing(Nmax_Trans_raytracing)
      ! integer, allocatable, dimension(:) :: tab_trans -> not useful anymore (?)
      integer, allocatable, dimension(:) :: i_trans, j_trans !return the index of a transition from 1 to atom%Ntr
      integer, allocatable, dimension(:,:) :: ij_to_trans !from i and j return the index of a transiton
      !Compatibility with RH, stored the collision in character format!
      character(len=512), allocatable, dimension(:) :: collision_lines !to keep all remaning lines in atomic file
      real(kind=dp)                :: cswitch, Abund, weight, massf !mass fraction
      integer                      :: periodic_table, activeindex !order of the active atom in the
      character(len=ATOM_LABEL_WIDTH), allocatable, dimension(:)  :: label
      integer                :: Nlevel, Nline, Ncont, Ntr, Ntr_line
      ! BY CONVENTION, stage=0 for neutrals, 1 for singly ionised
      integer, allocatable, dimension(:)  :: stage, Lorbit
      real(kind=dp) :: Rydberg
      real(kind=dp), allocatable, dimension(:) :: g, E
      real, allocatable, dimension(:) :: qS, qJ
      ! atom can be passive but NLTEpops true. This is the case of
      ! populations read from previous run
      logical                :: NLTEpops, set_ltepops
      logical :: lgauss_prof = .false.
      !real(kind=dp), allocatable :: ug(:), phi_g(:,:)
      real(kind=dp), dimension(:,:), pointer :: n, nstar
      type (AtomicLine), allocatable, dimension(:)         :: lines
      type (AtomicContinuum) , allocatable, dimension(:)   :: continua
   end type AtomType


   type atomPointerArray
      type(AtomType), pointer :: p => NULL()
   end type atomPointerArray

   type (atomPointerArray), dimension(:), allocatable :: Atoms, ActiveAtoms, PassiveAtoms
   type (AtomType), pointer :: Hydrogen => NULL(), Helium => NULL()
   integer :: n_atoms, Nactiveatoms, NpassiveAtoms
   logical :: helium_is_active

   integer, parameter :: cwitch_niter = 12!nint(ceil(log(cswitch_val)/log(cswitch_down_scaling_factor)))
   real(kind=dp), parameter :: cswitch_val = 1d12
   real(kind=dp), parameter :: cswitch_down_scaling_factor = 10.0!ceil(exp(log(cswitch_val)/cswitch_niter))
   logical :: cswitch_enabled = .false.
   
   contains

   function trans_number(atom, ilevel,jlevel)
      !from lower level i and upper level j
      !return the transition number from the first line 1 
      !to the last continuum = atom%Ntrans
      !
      !
      ! TO DO: check that the returned number is consistent !
      !        define an independent function not relying on atom%ij_to_trans
      type (AtomType), intent(in) :: atom
      integer, intent(in) :: ilevel, jlevel
      integer :: i, j, k , trans_number

      trans_number = atom%ij_to_trans(ilevel,jlevel)

      return

      ! do i=1, atom%Ntr
      !    k = atom%tab_trans(i)
      !    if(k < atom%Ntr_line + 1) then !line
      !       do j=1, atom%Nline
      !          if ((atom%lines(j)%i==ilevel).and.(atom%lines(j)%j==jlevel)) then
      !             write(*,*) "found transition it is a line", k, atom%lines(j)%lambda0
      !             write(*,*) atom%lines(j)%i, atom%lines(j)%j
      !             kr = k
      !             return
      !          endif
      !       enddo
      !    else !continua
      !       do j=atom%Ntr_line+1, atom%Ntr
      !          if ((atom%continua(j)%i==ilevel).and.(atom%continua(j)%j==jlevel)) then
      !             write(*,*) "found transition it is a continuum", k, atom%continua(j)%lambda0
      !             write(*,*) atom%continua(j)%i, atom%continua(j)%j
      !             kr = k
      !             return
      !          endif
      !       enddo
      !    endif
      ! enddo

      ! call error("(get_trans_number) did not find the transition ! ")
      ! return

   end function trans_number

   function rydberg_atom(atom)!correct from electron mass
      type (AtomType), intent(in) :: atom
      real(kind=dp) :: rydberg_atom, deltam

      deltam = 1. + mel / (atom%weight * amu_kg)

      rydberg_atom = E_RYDBERG / deltam !in joules
      return
   end function rydberg_atom

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
   function find_continuum(atom, l)
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

      return
   end function find_continuum

   function find_ground_state_ion(Nlevel, stage, Nstage)
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


      return
   end function find_ground_state_ion

   function atomic_orbital_radius(n, l, Z)
      !return atomic orbital radius wrt the Bohr radius
      real(kind=dp) :: atomic_orbital_radius
      real(kind=dp) :: n ! quantum principal number
      integer :: l ! orbital quantum number
      integer :: Z

      atomic_orbital_radius = n*n/real(Z) * (1d0 + 0.5*(1d0 - real(l)*(real(l)+1)/n/n))

      return
   end function atomic_orbital_radius

   function atomic_orbital_sqradius(n, l, Z)
      !Bates-Damguard mean square radius
      real(kind=dp) :: atomic_orbital_sqradius
      !in a0**2 units
      real(kind=dp) :: n ! quantum principal number
      integer :: l ! orbital quantum number
      integer :: Z ! charge : Z = 1 for H I

      atomic_orbital_sqradius = 0.5*n*n/real(Z*Z) * (5.*n*n + 1 - 3.*real(l)*(real(l)+1))


      return
   end function atomic_orbital_sqradius


   function getOrbital(orbit, unknown_orbital) result (L)
      integer :: L
      character(len=1), intent(in) :: orbit
      logical, intent(out) :: unknown_orbital

      unknown_orbital = .false.

      L = -1

      select case (orbit)
         case ('S')
            L = 0
         case ('P')
            L = 1
         case ('D')
            L = 2
         case ('F')
            L = 3
         case ('G')
            L = 4
         case ('H')
            L = 5
         case ('I')
            L = 6
         case ('J')
            L = 7
         case ('K')
            L = 8
         case ('L')
            L = 9
         case ('M')
            L = 10
         case ('N')
            L = 11
         case ('O')
            L = 12
         case ('Q')
            L = 13
         case ('R')
            L = 14
         case ('T')
            L = 15
         case ('U')
            L = 16
         case ('V')
            L = 17
         case ('W')
            L = 18
         case ('X')
            L = 19
         case DEFAULT
            write(*,'("Orbit " (1A1) "unknown")') orbit
            unknown_orbital = .true.
      end select
      return
   end function getOrbital


   !Not working correctly, using effective quantum Number instead. n**2 = Z**2 * Rydberg/deltaE
   function getPrincipal(label, n) result(determined)
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
            return
         else if (label(4:4)=="I") then
            if (label(5:5)=="I") then
               st = 7
            else if (label(5:5)==" ") then
               st = 7
            else
               write(*,*) "(2) Error parsing label for principal quantum number"
               determined = .false.
               return
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

      return
   end function getPrincipal

   !Should be more explicit and easy to parse at reading of the model!
   ! easy to break
   subroutine parse_label(label, g, S, L,J, determined)
      ! get principal quantum number from label
      real(kind=dp), intent(in) :: g
      character(len=ATOM_LABEL_WIDTH), intent(in) :: label
      logical, intent(out) :: determined
      integer, intent(out) :: L
      real, intent(out) :: J, S
    
    
      logical :: unknown_orbital
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
         return
      end if

      S = (real(S) - 1.)/2.
      L = getOrbital(orbit, unknown_orbital)
      if (unknown_orbital) then
         determined = .false.
         return
      endif

      J = real(L) + S!(g-1.)/2.
      if (J>(L+S)) then
         determined = .false.
      end if

      return
   end subroutine parse_label


   function atomZnumber(atom)
   ! --------------------------------------
   ! return the atomic number of an atom
   ! with ID = atomID.
   ! Hydrogen is 1
   ! --------------------------------------
      type (AtomType), intent(in) :: atom
      integer :: atomZnumber

      atomZnumber = atom%periodic_table

      return
   end function atomZnumber

  
!But why, a cswitch per atom ? It is going at the same speed for all atoms right ?
   function maxval_cswitch_atoms ()
      !for all atoms, check the maximum value of the cswitch
      integer :: n
      real(kind=dp) ::  maxval_cswitch_atoms

      maxval_cswitch_atoms = 1.0_dp
      do n=1, N_atoms

         maxval_cswitch_atoms = max(maxval_cswitch_atoms, atoms(n)%p%cswitch)

      enddo

      return
   end function maxval_cswitch_atoms

   !could be done for only one common parameter by the way
   subroutine adjust_cswitch_atoms ()
      !for all active atoms, decreases the cswitch value from cswitch_down_scaling_factor
      integer :: n
      logical :: print_message = .false.
      real(kind=dp) :: new_cs

      print_message = .false.

      do n=1, NactiveAtoms

         if (activeatoms(n)%p%cswitch > 1.0) then
            activeatoms(n)%p%cswitch = max(1.0_dp, min(activeatoms(n)%p%cswitch, &
               activeatoms(n)%p%cswitch/cswitch_down_scaling_factor))
            if (.not. print_message) then
               print_message = .true.
               new_cs = activeatoms(n)%p%cswitch
            endif
         endif
      enddo

      if (print_message) write(*,'(" cswitch for next iteration: "(1ES17.8E3))') new_cs


      return
   end subroutine adjust_cswitch_atoms

   function vbroad(temp, w, xi)
      real(kind=dp), intent(in) :: temp, w, xi
      real(kind=dp) :: vbroad
  
      vbroad = sqrt( Vtherm*temp/w + xi*xi )
  
      return
    end function vbroad

end module atom_type