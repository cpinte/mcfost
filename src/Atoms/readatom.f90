MODULE readatom

  use atom_type, only : AtomicLine, AtomicContinuum, AtomType, Element, determinate, rydberg_atom, &
  						n_eff, find_continuum
  use atmos_type, only : Nelem, Hydrogen, Helium, Elements, T, ne, vturb, lmagnetized, icompute_atomRT, &
  							Natom, NpassiveAtoms, NactiveAtoms, Atoms, PassiveAtoms, ActiveAtoms
  use zeeman, only : Lande_eff, ZeemanMultiplet
  use getlambda
  use constant
  use uplow
  use getline
  use barklem, only : getBarklem
  use io_atomic_pops, only	: read_pops_atom
  use collision, only : read_collisions
  use broad, only : Damping
  use solvene, only : Max_ionisation_stage, get_max_nstage

  !$ use omp_lib

  !MCFOST's originals
  use messages
  use mcfost_env, only : mcfost_utils ! convert from the relative location of atomic data
                                      ! to mcfost's environnement folders.
  use parametres, only : art_hv

  IMPLICIT NONE

  !-> Futur parameter or define in the model atom for each line ?
  !real(kind=dp), parameter :: maxvel_atom_transfer = 50.0 !in km/s
  
  character, parameter :: COMMENT_CHAR="#"
  character(len=*), parameter :: ATOMS_INPUT = "./atoms.input"!"/Atoms/atoms.input"
  character(len=*), parameter :: path_to_atoms = "/Atoms/"
  !real, parameter :: MAX_ABUND_ERROR=0.001
  
  integer, parameter :: cwitch_niter = 12!nint(ceil(log(cswitch_val)/log(cswitch_down_scaling_factor)))
  real(kind=dp), parameter :: cswitch_val = 1d12
  real(kind=dp), parameter :: cswitch_down_scaling_factor = 10.0!ceil(exp(log(cswitch_val)/cswitch_niter))
  logical :: cswitch_enabled = .false.
  
  !Until an extrapolation routine is found, the lambdamax is found where Gaunt(lambda) < 0
  real, parameter :: fact_pseudo_cont = 100.0 ! the hydrogenic bound-free will be extrapolated
  						!beyond the edge limit (lambda0) up to lambda0 * fact_pseudo_cont
  						!if it is 1.0, no extrapolation i.e., lambdamax=lambda0

  CONTAINS

  SUBROUTINE readModelAtom(atomunit, atom, atom_file)
   !!!
   ! read independent atomic model
   ! Iinitialize the atom values
   ! Variables depending on the size of the grid
   ! for instance Rij/Rji at each grid points
   ! are allocated as 1 dim (flattened) arrays
   !!!
    integer, intent(in) :: atomunit

    integer :: kr, k, la, alloc_status
    type (AtomType), intent(inout), target :: atom
    character(len=*), intent(in) :: atom_file
    character(len=MAX_LENGTH) :: inputline, FormatLine
    integer :: Nread, i,j, EOF, nll, nc, Nfixed !deprecation future
    real, allocatable, dimension(:) :: levelNumber
    logical :: Debeye, match, res=.false.
    logical, dimension(:), allocatable :: determined, parse_label
    real(kind=dp), dimension(:), allocatable :: old_nHtot
    character(len=20) :: shapeChar, symmChar, optionChar, vdWChar, nuDepChar
    character(len=2) :: IDread
    real(kind=dp) :: C1, vDoppler, f, lambdaji
    real(kind=dp) :: lambdamin, geff, gamma_j, gamma_i
    EOF = 0
    !write(*,*) "Atom is part of the active atoms ?", atom%active

    !for Aji
    C1 = 2.*PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON / CLIGHT)    

    open(unit=atomunit,file=atom_file,status="old")
    !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
    write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH

    !Read ID and fill atom abundance and weight
    CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
    read(inputline,*) IDread

    IDread(2:2) = to_lower(IDread(2:2))
    atom%ID = IDread
    write(*,*) "Reading atomic model of atom ", atom%ID
    match = .false.
    do nll=1,Nelem
     if (Elements(nll)%ptr_elem%ID.eq.atom%ID) then
      atom%massf = Elements(nll)%ptr_elem%massf
      write(*,*) "Abundance of atom ",atom%ID,": A =",Elements(nll)%ptr_elem%Abund
      if (atom%ID=="H" .or. atom%ID=="He") then 
       write(*,*) " -> mass fraction (%) = ", 100.*real(atom%massf)
      else
       write(*,*) " -> mass fraction (m/m(Fe) %) = ", 100.*real(atom%massf/Elements(26)%ptr_elem%massf)
      endif
      if (Elements(nll)%ptr_elem%abundance_set) then
        atom%periodic_table=nll
        atom%Abund=Elements(nll)%ptr_elem%Abund
        atom%weight=Elements(nll)%ptr_elem%weight
        atom%Rydberg = rydberg_atom(atom)
        match=.true.
      end if
      exit
     end if
    end do

    if (match .neqv..true.) then
     write(*,*) "Error no abundance found for atom ", atom_file, " ", atom%ID
     stop
    end if

    !read Nlevel, Nline, Ncont and Nfixed transitions
    CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
    read(inputline,*) atom%Nlevel, atom%Nline, atom%Ncont, Nfixed
    write(*,*) "Nlevel=",atom%Nlevel," Nline=",atom%Nline,&
              " Ncont=", atom%Ncont!, " Nfixed=", atom%Nfixed
              !deprecated Nfixed will be removed

    !read for each level, Energie (%E), statistical weight (%g), label,
    ! stage and levelNo
    
    atom%cswitch = 1.0_dp

    allocate(atom%label(atom%Nlevel))
    allocate(atom%E(atom%Nlevel))
    allocate(atom%g(atom%Nlevel))
    allocate(atom%stage(atom%Nlevel)) ! note that stage is 0 for neutrals
    allocate(levelNumber(atom%Nlevel))

    allocate(atom%Lorbit(atom%Nlevel))
    allocate(atom%qS(atom%Nlevel))
    allocate(atom%qJ(atom%Nlevel))
    allocate(determined(atom%Nlevel))
    allocate(parse_label(atom%Nlevel))
    allocate(atom%nstar(atom%Nlevel,n_cells))

    atom%Ntr = atom%Nline + atom%Ncont
    allocate(atom%at(atom%Ntr))
	atom%Ntr_line = atom%Nline

    atom%nstar(:,:) = 0d0

    do i=1,atom%Nlevel
     CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
     read(inputline,*) atom%E(i), atom%g(i), atom%label(i), &
                       atom%stage(i), levelNumber(i)
     atom%E(i) = atom%E(i) * HPLANCK*CLIGHT / (CM_TO_M)
!      write(*,*) "E(eV) = ", atom%E(i) * JOULE_TO_EV, &
!                "g = ", atom%g(i), "label = ", atom%label(i), &
!                "stage = ", atom%stage(i)

      !Default value:  Note that -99 is used for
      !unknown values:
      !either Landé, or S, J, L etc.
      !If usefull, levels with -99 will be skiped
      !for instance for zeeman polar.
       atom%qJ(i) = -99d0
       atom%qS(i) = -99d0
       atom%Lorbit(i) = -99
       determined(i) = .false.
       parse_label(i) = .false.
    end do  ! Note that the levelNumber starts at 0, but indexing in
            ! fortran starts at 1 so that the good number to use
            ! here is +1 wrt the the number read in file.

    ! Check if there is at least one continuum transitions
    if (atom%stage(atom%Nlevel) /= atom%stage(atom%Nlevel-1)+1) then
     write(*,*) atom%stage
     write(*,*) atom%stage(atom%Nlevel), atom%stage(atom%Nlevel-1)+1
     write(*,*) "Atomic model does not have an overlying continuum"
     write(*,*) "exiting..."
     stop
    end if
    !! DO NOT USE i as a loop index here !!

    !deprecation, they are now computed cell by cell for saving memory
    !allocate(atom%ntotal(n_cells)) !-> but cheap in mem? so keep it ?
!-> futur deprecation compute localy ??? or not
    allocate(atom%vbroad(n_cells)) !-> if profile keep in memory needed ??
    !write(*,*)
    !VDoppler = KBOLTZMANN/(AMU * atom%weight) * 8d0/PI !* 2d0!m/s

!-> futur deprec, not need to store
    atom%vbroad = sqrt(vtherm/atom%weight * T + vturb**2) !vturb in m/s
    !atom%ntotal = atom%Abund * nHtot
    !note: for Hydrogen, ntotal is actually (nHtot - nHmin)

    VDoppler = sqrt(vtherm/atom%weight * maxval(T) + maxval(vturb)**2)
    

!     write(*,*) Vdoppler, sqrt(Vtherm*maxval(T)/atom%weight + maxval(vturb)**2)
! 	do k=n_cells, 1, -1
! 	 if (icompute_atomRT(k)>0) then
! 	  write(*,*) k, T(k), atom%ID, atom%vbroad(k)/1e3, vtherm/atom%weight / 1e3, vturb(k)/1e3
! 	 endif
! 	enddo
! stop
    !Now read all bound-bound transitions
    allocate(atom%lines(atom%Nline))

    do kr=1,atom%Nline
     atom%lines(kr)%atom => atom
     atom%at(kr)%trtype = atom%lines(kr)%trtype; atom%at(kr)%ik=kr
     atom%at(kr)%lcontrib_to_opac=.true.

     atom%lines(kr)%isotope_frac = 1.
     atom%lines(kr)%g_lande_eff = -99.0
     atom%lines(kr)%glande_i = -99; atom%lines(kr)%glande_j = -99
     !atom%lines(kr)%trtype="ATOMIC_LINE"
     CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
     Nread = len(trim(inputline)) ! because, if blanck
           ! beyond cStark it will be interpreted
           ! has "additional geff", but its not.
           !
     !futur implement: line%name
     atom%lines(kr)%ZeemanPattern = 1 !should be read in file
     if (.not.lmagnetized) atom%lines(kr)%ZeemanPattern = 0
     ! -1 = effective triplet, +1
!      write(*,*) "Reading line #", kr
     if (Nread.eq.112) then
         read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
         symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
         atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
         atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
         atom%lines(kr)%Grad, atom%lines(kr)%cStark
     else if (Nread.gt.112) then
       write(*,*) " ->Read aditional g_lande_eff for that line"
       read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
       symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
       atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
       atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
       atom%lines(kr)%Grad, atom%lines(kr)%cStark, &
       atom%lines(kr)%g_Lande_eff, atom%lines(kr)%glande_j, atom%lines(kr)%glande_i
       !if glande_eff given, we need to read a value for gi and gj even if it is 0.
       !if glane is <= -99, gi, and gj and geff are computed eventually.
       !landé upper / lower levels in case the coupling scheme is not accurate
       if (atom%lines(kr)%g_lande_eff > -99) atom%lines(kr)%ZeemanPattern = -1 !effective T assumed
       if (atom%lines(kr)%g_lande_eff <= -99) atom%lines(kr)%ZeemanPattern = 1 !Full using gi and gj
       if (atom%lines(kr)%glande_j <= -99 .or. atom%lines(kr)%glande_i<= -99) then
        CALL Warning("Unable to use read lande factors, try to compute them..")
        atom%lines(kr)%g_lande_eff = -99 !force calculation
       end if
     else
       write(*,*) inputline
       CALL error(" Unable to parse atomic file line")
     end if
      i = i + 1
      j = j + 1 !because in C, indexing starts at 0, but starts at 1 in fortran


!      write(*,*) "Reading line #", kr, 1d9 * (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i)), 'nm'

      !therefore, the first level is 1 (C=0), the second 2 (C=1) etc
      !Lymann series: 2->1, 3->1
      atom%lines(kr)%i = min(i,j)
      atom%lines(kr)%j = max(i,j)
         
      !tmp
      if (atom%ID=="H") then
      
      	!Ly alpha
      	if (atom%g(i)==2 .and. atom%g(j)==8) atom%lines(kr)%write_flux_map =.true. 
      	!H alpha
      	if (atom%g(i)==8 .and. atom%g(j)==18) atom%lines(kr)%write_flux_map=.true.
      	!H beta
      	if (atom%g(i)==8 .and. atom%g(j)==32) atom%lines(kr)%write_flux_map =.true. 
      	!H gamma
      	if (atom%g(i)==8 .and. atom%g(j)==50) atom%lines(kr)%write_flux_map =.true. 
      	!Pa beta
      	if (atom%g(i)==18 .and. atom%g(j)==50) atom%lines(kr)%write_flux_map =.true. 
      	!Br gamma
      	if (atom%g(i)==32 .and. atom%g(j)==98) atom%lines(kr)%write_flux_map =.true. 


	  elseif (atom%ID=="He") then
	  	if (kr >= 4 .and. kr <= 6) then
	  		atom%lines(kr)%write_flux_map =.true. 
	  	endif
      else
      
      	if (kr <= 4) then
      		atom%lines(kr)%write_flux_map =.true. 
      	endif

      endif
      
      if (atom%lines(kr)%qwing < 2.0) then
      
      	call Warning("qwing for line lower than 2! setting to 2")
      	atom%lines(kr)%qwing = 2.0
      	
      endif


!       write(*,*)  j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
!       symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
!       atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
!       atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
!       atom%lines(kr)%Grad, atom%lines(kr)%cStark, &
!       atom%lines(kr)%g_Lande_eff

      ! filling S (2S+1), J and Lorbit for each line level
      !determine only for lines
      ! because continuum levels are not used for polari-
      ! -sation yet.
      ! continuum levels parsing is usefull
      ! to obtain W2(Jj,Ji) appearing in
      ! resonant dichroism
      if (.not.parse_label(atom%lines(kr)%i)) then
      CALL determinate(atom%label(atom%lines(kr)%i),&
       atom%g(atom%lines(kr)%i),&
       atom%qS(atom%lines(kr)%i),&
       atom%Lorbit(atom%lines(kr)%i),&
       atom%qJ(atom%lines(kr)%i), &
       determined(atom%lines(kr)%i))
       parse_label(atom%lines(kr)%i) = .true.
      end if
      if (.not.parse_label(atom%lines(kr)%j)) then
      CALL determinate(atom%label(atom%lines(kr)%j),&
       atom%g(atom%lines(kr)%j),&
       atom%qS(atom%lines(kr)%j),&
       atom%Lorbit(atom%lines(kr)%j),&
       atom%qJ(atom%lines(kr)%j), &
       determined(atom%lines(kr)%j))
       parse_label(atom%lines(kr)%i) = .true. !even if determined is false.
      end if
      ! not that if J > L+S determined is FALSE
      ! just like if the term could not be parsed
      ! Because we parse only atomic lines
      !  continuum transitions are by default not determined
      ! even if they are parsable
      !if (.not. determined(i) .and. i.le.atom%Nline) then
      !write(*,*) "Could not parssed level ", i, atom%label(i)
      !if (atom%qS(i)+atom%Lorbit(i).lt.atom%qJ(i)) then
      !  write(*,*) "J > L+S: term not allowed!"
      !end if
      !end if
      !write(*,'("S="(1F2.2)", L="(1I2)", J="(1F2.2))') &
      !  atom%qS(i), atom%Lorbit(i), atom%qJ(i)
      !!
      !If Selection rule is OK and g_lande_eff not given from file and atomic label
      !correctly determined
      if ((abs(atom%qJ(atom%lines(kr)%i) - atom%qJ(atom%lines(kr)%j)) <= 1.) .and. &
           (atom%lines(kr)%g_Lande_eff <= -99) .and. &
           (determined(atom%lines(kr)%j)) .and. (determined(atom%lines(kr)%i))) then !
       ! do not compute geff if term is not
       ! determined or if the geff is read from file
       ! ie if g_lande_eff > -99
       ! fill lande_g_factor if determined and not given
       ! for b-b transitions
       CALL Lande_eff(atom, kr)
       !compute indiviual glande of levels also
       !!testing
       !!write(*,*) wKul(atom, kr, 0)**2
       !!write(*,*) wKul(atom, kr, 1)**2
       !!write(*,*) wKul(atom, kr, 2)**2
       !!stop
       !write(*,*) "geff = ", atom%lines(kr)%g_lande_eff
      end if
      !!atom%has_atomic_polarization = .false. !different criterion than polarizable
      !!atom%lines(kr)%has_alignement = .false. !depending on the polarisability factor
      											! if we neglect J-states coherences
      atom%lines(kr)%polarizable = .false. !check also deltaJ
      !write(*,*) "dJ=", abs(atom%qJ(atom%lines(kr)%i) - atom%qJ(atom%lines(kr)%j))
      atom%lines(kr)%polarizable = (lmagnetized) .and. &
      								(atom%lines(kr)%g_lande_eff > -99)! .and. (abs(atom%qJ(atom%lines(kr)%i) - atom%qJ(atom%lines(kr)%j)) <= 1.)
       !not need to be determined here. Because geff can be read from file and be > -99
       !even if the levels are not determined. In this case deltaJ = 0 (-99+99).
       !Exception if one of the level is determined but not the other, in this case
       !line is assumed to be not polarizable and you have to chandle that in the atomic file.
       ! Otherwise I assume you know what you do by providing a landé factor to a line.
       !In case of PRT but the line is not polarized, set a Zeeman structure with Nc=1,
       !S=0, q=0, shift=0
       ! line%ZP = 0 if WEAK_FIEKD solution; otherwise has to be -1, 1 or .not.polarizable
       if (atom%lines(kr)%ZeemanPattern /= 0) CALL ZeemanMultiplet(atom%lines(kr))

      ! oscillator strength saved
      atom%lines(kr)%fosc = f
      lambdaji = (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i))
      atom%lines(kr)%Aji = C1 / (lambdaji**2) * (atom%g(i) / atom%g(j)) * f
      atom%lines(kr)%Bji = (lambdaji**3) / (2.0 * HPLANCK * CLIGHT) &
                          *atom%lines(kr)%Aji
      atom%lines(kr)%Bij = (atom%g(j) / atom%g(i)) * atom%lines(kr)%Bji
      atom%lines(kr)%lambda0 = lambdaji / NM_TO_M

      !gi * Bij = gj * Bji
      !gi/gj = Bji/Bij so that niBij - njBji = Bij (ni - gij nj)
      atom%lines(kr)%gij = atom%lines(kr)%Bji / atom%lines(kr)%Bij !gi/gj
      atom%lines(kr)%twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji

!       write(*,*) " ->", " Aji (1e7 s^-1) = ", atom%lines(kr)%Aji/1d7,&
!         "Grad (1e7 s^-1) = ", atom%lines(kr)%Grad/1d7, &
!         "gj = ", atom%g(j)," gi = ",  atom%g(i)

      !write(*,*) "line ", atom%lines(kr)%j,'->',atom%lines(kr)%i, " @",&
      !           lambdaji/NM_TO_M," nm : Aji = ", &
      !           atom%lines(kr)%Aji, " Bji = ", atom%lines(kr)%Bji,&
      !           " Bij = ", atom%lines(kr)%Bij

     ! Now parse line string, used to construct the profile function
     if (trim(shapeChar).ne."PFR" .and. trim(shapeChar).ne."VOIGT" &
         .and. trim(shapeChar).ne. "GAUSS") then
         write(*,*) "Invalid value for line-shape string"
         write(*,*) "exiting..."
         stop
     else if (trim(shapeChar).eq."PFR") then
        write(*,*) "Presently PFR transitions not allowed"
        stop
        atom%Npfr = atom%Npfr  + 1
        atom%lines(kr)%PFR = .true.
     else if (trim(shapeChar).eq."GAUSS") then
        !write(*,*) "Using Gaussian profile for that line"
        atom%lines(kr)%Voigt = .false.
!      else if (trim(shapeChar).eq."COMPOSIT") then
!         !write(*,*) "Using a Multi-component profile for that line"
!         ! frist read the number of component
!         CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
!         read(inputline,*) atom%lines(kr)%Ncomponent
!         !write(*,*) "Line has ", atom%lines(kr)%Ncomponent, " components"
!         ! read shift and fraction
!         allocate(atom%lines(kr)%c_shift(atom%lines(kr)%Ncomponent))
!         allocate(atom%lines(kr)%c_fraction(atom%lines(kr)%Ncomponent))
!         c_sum = 0.
!         do nc=1, atom%lines(kr)%Ncomponent
!          CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
!          read(inputline, *) atom%lines(kr)%c_shift(nc), &
!                            atom%lines(kr)%c_fraction(nc)
!          c_sum = c_sum + atom%lines(kr)%c_fraction(nc);
!         end do
!         if (c_sum.gt.(1.+MAX_ABUND_ERROR) .or. &
!             c_sum.lt.(1.-MAX_ABUND_ERROR)) then
!             write(*,*) "Line component fractions do not add up to unity"
!             write(*,*) "exiting..."
!             stop
!         else
!             !write(*,*) "Line ",atom%lines(kr)%j, "->",atom%lines(kr)%i,&
!             !  " has ",  atom%lines(kr)%Ncomponent," components"
!         end if

     else !default is Voigt ! line%voigt is default .true.
      atom%lines(kr)%Voigt = .true.
        !write(*,*) "Using Voigt profile for that line"
        !allocate(atom%lines(kr)%c_shift(1))
        !allocate(atom%lines(kr)%c_fraction(1))
        !atom%lines(kr)%Ncomponent = 1
        !atom%lines(kr)%c_shift(1) = 0.0
        !atom%lines(kr)%c_fraction(1) = 1.
     end if
     
 !force Gaussian for test
! CALL Warning("USING GAUSSIAN LINE PROFILES")
!      atom%lines(kr)%Voigt = .false.

     !Now parse Broedening recipe
     if (trim(vdWChar).eq."PARAMTR") then
       atom%lines(kr)%vdWaals = "RIDDER_RENSBERGEN"
       CALL Error ("Not used anymore!")
     else if (trim(vdWChar).eq."UNSOLD") then
       !write(*,*) "Using UNSOLD recipe for broadening"
       atom%lines(kr)%vdWaals = "UNSOLD"
       atom%lines(kr)%cvdWaals(4) = 0.
       atom%lines(kr)%cvdWaals(2) = 0.
     else if (trim(vdWChar).eq."BARKLEM") then
      atom%lines(kr)%vdWaals = "BARKLEM"
      CALL getBarklem(atom, kr, res)
      !if (res) &
        write(*,*) "Using BARKLEM recipe for broadening"
      if (.not. res) then
       write(*,*) &
         "Line <atom%lines(kr)%j>->atom%lines(kr)%i>", &
         " cannot be treated with Barklem type", &
         " broadening."
       write(*,*) "using UNSOLD"
       atom%lines(kr)%vdWaals = "UNSOLD"
       atom%lines(kr)%cvdWaals(4) = 0.
       atom%lines(kr)%cvdWaals(2) = 0.
      end if
     else
       write(*,*) 'Invalid value for vdWaals broadening reicpe'
       write(*,*) "exiting..."
       stop
     end if

     !symmetric without magnetic field, means that we can compute locally a line profile
     !only for half of the profile. 
     !!Futur deprecation
     if (trim(symmChar).eq."ASYMM") then
      atom%lines(kr)%symmetric = .false.
     else
      atom%lines(kr)%symmetric = .true.
      !write(*,*) "Symmetric line profile"
     end if
     !Should be replaced by another flag: for instance with magnetic field we need to use
     !week field if we shift profiles otherwise the profile is not symm
     atom%lines(kr)%symmetric = .false.


!      if (atom%active) then !Should do it for passive atoms too
!      if (lmagnetized) then !.or.line%scattpol ...
!       if (atom%lines(kr)%g_Lande_eff.gt.-99 .or. &
!           determined(atom%lines(kr)%i) .and. &
!           determined(atom%lines(kr)%j).and. &
!           abs(atom%qJ(atom%lines(kr)%i) - &
!            atom%qJ(atom%lines(kr)%j)).le.1.) then
!
! !        if (atom%lines(kr)%Ncomponent.gt.1) then
! !            !write(*,*) &
! !            !"Cannot treat composite line with polar"
! !            atom%lines(kr)%polarizable=.false.
! !        else
!          atom%lines(kr)%polarizable=.true.
! !        end if
!       end if
!      else
!       !write(*,*) "Treating line ",atom%lines(kr)%j,&
!       !    "->",atom%lines(kr)%i,&
!       !    " without polarization"
!       atom%lines(kr)%polarizable=.false.
!      end if !not mag
!     end if ! end loop over active b-b transitions of atom
   end do !end loop over bound-bound transitions
   
    ! ----------------------------------------- !
    !starts reading bound-free transitions
    ! cross-sections allocated once the final
    ! wavelength grid is known.
    ! ----------------------------------------- !
    allocate(atom%continua(atom%Ncont))
    do kr=1,atom%Ncont
     atom%continua(kr)%isotope_frac=1.
     atom%continua(kr)%atom => atom
     atom%at(kr)%lcontrib_to_opac=.true.
     !write(*,*) atom%Ntr_line+kr, kr, atom%nline+kr
     atom%at(atom%Nline+kr)%trtype = atom%continua(kr)%trtype; atom%at(kr+atom%Nline)%ik=kr
     !write(*,*)  atom%at(kr+atom%Nline)%ik, kr

     CALL getnextline(atomunit, COMMENT_CHAR, &
          FormatLine, inputline, Nread)
     read(inputline, *) j, i, atom%continua(kr)%alpha0,&
      atom%continua(kr)%Nlambda, nuDepChar, atom%continua(kr)%lambdamin
     j = j + 1
     i = i + 1
!      write(*,*) "Reading continuum #", kr, atom%continua(kr)%lambdamin, "nm", &
!      					1d9 * (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i)), "nm"

     !because in C indexing starts at 0, but starts at 1 in fortran
     atom%continua(kr)%j = max(i,j)
     atom%continua(kr)%i = min(i,j)
     lambdaji = (HPLANCK*CLIGHT)/ (atom%E(j)-atom%E(i))
     atom%continua(kr)%lambda0 = lambdaji/NM_TO_M !nm
     atom%continua(kr)%lambdamax = atom%continua(kr)%lambda0
     !wavelength max for extrapolated bound-free cross-section
     !used only if continuum is hydrogenic

     !write(*,*) "continuum ", atom%continua(kr)%j,&
     !    '->',atom%continua(kr)%i, " @",&
     !    lambdaji/NM_TO_M," nm : alpha0 = ", &
     !    atom%continua(kr)%alpha0, &
     !    " Nlambda = ", atom%continua(kr)%Nlambda, &
     !    " type = ", nuDepChar," lambda min = ",&
     !     lambdamin," nm"

     if (trim(nuDepChar).eq."EXPLICIT") then
      ! Nlambda set in atomic file
      allocate(atom%continua(kr)%alpha_file(atom%continua(kr)%Nlambda))
      allocate(atom%continua(kr)%lambda_file(atom%continua(kr)%Nlambda))
      allocate(atom%continua(kr)%lambda(atom%continua(kr)%Nlambda))
      atom%continua(kr)%hydrogenic=.false.
      ! rearanging them in increasing order
      ! because in the atomic file they are
      ! given in decreasing order !
      do la=atom%continua(kr)%Nlambda,1,-1
       CALL getnextline(atomunit, COMMENT_CHAR, &
             FormatLine, inputline, Nread)
       read(inputline,*) atom%continua(kr)%lambda_file(la), &
          atom%continua(kr)%alpha_file(la)
       ! though they are printed in decreasing order
       !write(*,*) "l = ",atom%continua(kr)%lambda(la), &
       !    " a = ", atom%continua(kr)%alpha(la)
      end do
      atom%continua(kr)%lambda(:) = atom%continua(kr)%lambda_file(:)
      do la=2,atom%continua(kr)%Nlambda
        if (atom%continua(kr)%lambda_file(la).lt.&
           atom%continua(kr)%lambda_file(la-1)) then
          write(*,*) "continuum wavelength not monotonous"
          write(*,*) "exiting..."
          stop
        end if
      end do
      !not extrapolated if explicit ?
      !should consider neglecting occupation probability for that case
      atom%continua(kr)%lambdamax = maxval(atom%continua(kr)%lambda_file)
     else if (trim(nuDepChar).eq."HYDROGENIC") then
       atom%continua(kr)%hydrogenic=.true.
       !!tmp
       atom%continua(kr)%lambdamin = 5.0_dp
       !atom%continua(kr)%lambdamin = 0.05 * atom%continua(kr)%lambda0
       !atom%continua(kr)%lambdamin = max(5.0_dp, 1d-3 * atom%continua(kr)%lambda0)
       !!tmp
       write(*,'(" Continuum "(1I3)" -> "(1I3)" at "(1F12.5)" nm")') atom%continua(kr)%i, atom%continua(kr)%j, atom%continua(kr)%lambda0
       write(*,'(" -> lower edge cut at "(1F12.5)" nm !")'), atom%continua(kr)%lambdamin       
       
       if (atom%continua(kr)%lambdamin>=atom%continua(kr)%lambda0) then
        write(*,*) "Minimum wavelength for continuum is larger than continuum edge."
        write(*,*) kr, atom%continua(kr)%lambda0, atom%continua(kr)%lambdamin
        write(*,*) "exiting..."
        stop
       end if
       !used to extrapolate the b-f cross-section for level dissolution
       !if I go to far in lambda0 actually g_bf is < 0 and atm I set g_bf = 0 if g_bf < 0.
       !to go further in lambda (i.e. for lambda with g_bf(lambda)->0), I need to properly extrapolated g_bf.
       !!atom%continua(kr)%lambdamax = atom%continua(kr)%lambda0 * fact_pseudo_cont!7
       !Meanwhile, I search the lambda for which g_bf is < 0 (then 0)
       if (ldissolve) then !only if we actually want to extrapolate
       	if (atom%ID=="H") then ! .or. atom%ID=="He") then
        	CALL search_cont_lambdamax (atom%continua(kr), atom%Rydberg, atom%stage(i)+1,atom%E(j),atom%E(i)) 
            CALL make_sub_wavelength_grid_cont_linlog(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
			!!CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)   
		else
			!there is dissolve but not for this atom
            CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
! 			call make_sub_wavelength_grid_cont_log_nu(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
       	endif
       else !no dissolve
			CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)  
! 			call make_sub_wavelength_grid_cont_log_nu(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
       endif
       ! %lambda allocated inside the routines.
!        CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
       !make it log if not occupation probability ?
!        CALL make_sub_wavelength_grid_cont_log(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
	   !Can be done elsewhere, with atomic lines ? But unlike some lines grid, this does not depend on the local condition ATM
     else
      write(*,*) "Invalid continuum type : ", trim(nuDepChar)
      write(*,*) "exiting..."
      stop
     end if
    end do !end loop over bound-free transitions

    ! now fixed transitions
    ! fixed transitions are usefull to describe NLTE problem
    ! in the solar chromosphere in 2D,3D without not so much
    ! computational cost.
    ! They are not implemented because only relevent for solar
    ! application
    if (Nfixed.gt.0) then
     write(*,*) "Fixed transitions not implemented yet"
     write(*,*) "exiting..."
     stop
    end if !end reading fixed trans

    !now compute wavelengths grid for each line done elsewhere

   !Now even for passive atoms we write atomic data.
!     if (atom%ID(2:2) .eq." ") then
!       atom%dataFile = atom%ID(1:1)//".fits.gz" !.fits to be updated, .gz not
!     else
!       atom%dataFile = atom%ID(1:2)//".fits.gz"
!     end if
	

   atom%set_ltepops = .true. !by default compute lte populations
   atom%NLTEpops = .false.
   
    ! allocate some space
    if (atom%initial_solution.eq."ZERO_RADIATION") then
    	if (.not.atom%active) then
    		write(*,*) atom%ID, " is passive! cannot use ZERO_RADIATION solution, set to LTE."
    		atom%initial_solution="LTE_POPULATIONS"
    	endif
    endif

   if (atom%active) then

    !Not implemented, futur removal in atomic file
    if (atom%Npfr.gt.0) then
     write(*,*) "PFR not implemented yet, do not write file"
    end if

    ! reading collision rates of RH
    if (atom%ID /= "H") then
    	write(*,*) "  -> Reading collision data from RH for atom ", atom%ID
! 		atom%colunit = atom%periodic_table*2 + 1
     call read_collisions(atomunit, atom)
    endif

    !!allocate(atom%C(atom%Nlevel*atom%Nlevel,n_cells))
    !!now Collision matrix is constructed cell by cell, therefore allocated elsewhere
    allocate(atom%n(atom%Nlevel,n_cells))
    atom%n = 0d0
    if (atom%initial_solution .eq. "LTE_POPULATIONS") then
    	atom%n = atom%nstar !still need to be computed
    	atom%set_ltepops = .true.
    	write(*,*) " -> Setting initial solution to LTE "
    	
    else if (atom%initial_solution .eq. "CSWITCH") then
    	atom%n = atom%nstar !still need to be computed
    	atom%set_ltepops = .true.
    	if (.not. lforce_lte) then !otherwise cswitch is not LTE
    		write(*,*) " -> Setting initial solution to LTE with CSWITCH "
    		atom%cswitch = cswitch_val
    		if (cswitch_enabled == .false.) cswitch_enabled = .true.!we need at least one
    	endif
    else if (atom%initial_solution .eq. "ZERO_RADIATION") then
    
    	atom%n = atom%nstar
    	atom%set_ltepops = .true.
    	!nlte pops is false
    	    	
	else if (atom%initial_solution .eq. "OLD_POPULATIONS") then
	
	   write(*,*) " -> Reading (non-LTE AND LTE) populations from file..."
       CALL read_pops_atom(atom)
       atom%NLTEpops = .true.
       atom%set_ltepops = .false. !read and USE also LTE populations from file!!

    end if
   else !not active = PASSIVE
     if (atom%initial_solution .eq. "OLD_POPULATIONS") then
     
       allocate(atom%n(atom%Nlevel,n_cells)) !not allocated if passive, n->nstar
	   write(*,*) " -> Reading (non-LTE AND LTE) populations from file for passive atom..."
       CALL read_pops_atom(atom)
       atom%NLTEpops = .true.
       atom%set_ltepops = .false. !read and USE also LTE populations from file!!

       !atom%NLTEpops = .false. still false at this point as we need pops to do electron densities
     else !pure passive without nlte pops from previous run
!        atom%NLTEpops=.false.  !-> default values
!        atom%set_ltepops = .true.
       atom%n => atom%nstar !initialised to zero
     !atom%n is an alias for nstar in this case
     end if
   end if !end is active
   
   !check nHtot if we read NLTE populations
!    if (atom%NLTEpops .and. atom%ID=='H') then !NLTE pops is false if departure coefficients
!     write(*,*) "Using NLTE populations for total H density"
!     allocate(old_nHtot(n_cells)); old_nHtot = 0.0
!     old_nHtot = nHtot
!     write(*,*) " -> old max/min nHtot (m^-3)", maxval(nHtot), minval(nHtot)
!     !nHtot = 0.
!     !nHtot = sum(atom%n,dim=1)
!     
!     write(*,*) " -> new max/min nHtot (m^-3)", maxval(nHtot), minval(nHtot)
!     write(*,*) "    :: max/min ratio", maxval(nHtot)/maxval(old_nHtot,mask=old_nHtot>0), &
!       minval(nHtot,mask=nHtot>0)/minval(old_nHtot,mask=old_nHtot > 0)
!     deallocate(old_nHtot)
!    endif
   
!    atom%n = 0.
!    atom%n => atom%nstar
!    atom%NLTEpops = .false.
!    atom%set_ltepops = .false.
!    if (loc(atom%n) /= loc(atom%nstar)) then
!     call error ("pointers error")
!    endif

  deallocate(levelNumber)
  deallocate(determined)
  deallocate(parse_label)
    !close atomic file
  close(unit=atomunit) !later it will be open again for
        !reading collision
  RETURN
  END SUBROUTINE readModelAtom

  SUBROUTINE readAtomicModels(unit)
  !Read all atomic files present in atoms.input file
  !successive call of readModelAtom()
  integer :: EOF=0,Nread, nmet, mmet, nblancks, nact, npass
  integer :: kr, k, imax, ic
  real, parameter :: epsilon = 5e-3
  real :: eps
  real(kind=dp) :: epsilon_l_max !if epsilon > 1/pi/adamp, the value of xwing_lorentz is negative
  real(kind=dp) :: max_adamp, adamp, maxvel, vel, min_resol
  integer, intent(in) :: unit
  character(len=MAX_LENGTH) :: inputline
  character(len=15) :: FormatLine
  character(len=MAX_LENGTH) :: popsfile, filename
  character(len=MAX_KEYWORD_SIZE) :: actionKey, popsKey
  character(len=2) :: IDread
  type (AtomType), pointer :: atom

  !create formatline
  !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
  !write(FormatLine,'("("I3,A")")') MAX_LENGTH,"A"
  write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH

!   if (fact_pseudo_cont <= 0.0) then
!    write(*,*) "fact_pseudo_cont=", fact_pseudo_cont
!    call error("Wrong value for fact_pseudo_cont!")
!   endif
!   if (fact_pseudo_cont > 1.0) then
!    Write(*,*) " Hydrogenic continua extrapolated up to lambda0 x ", fact_pseudo_cont
!   endif

  Nactiveatoms = 0
  Npassiveatoms = 0
  open(unit=unit,file=TRIM(ATOMS_INPUT), status="old")!mcfost_utils)//TRIM(ATOMS_INPUT)

  !get number of atomic models to read
  CALL getnextline(unit, COMMENT_CHAR,FormatLine, &
       inputline, Nread)
  read(inputline,*) Natom
  write(*,*) "Reading ", Natom, " species"
  !Allocate sapace for Natom in Atoms
  allocate(Atoms(Natom))

  do nmet = 1, Natom
   CALL getnextline(unit, COMMENT_CHAR, FormatLine, inputline, Nread)

   allocate(Atoms(nmet)%ptr_atom)

   read(inputline,'(1A28, 1A7, 1A22, 1A20)') filename, actionKey, popsKey, popsFile
   !write(*,*) ".",trim(filename),"."
   !write(*,*) ".",adjustl(actionKey),"."
   !write(*,*) ".",adjustl(popsKey),"."
   !write(*,*) ".",trim(popsFile),"."

   Atoms(nmet)%ptr_atom%initial_solution=adjustl(popsKey)
   Atoms(nmet)%ptr_atom%inputFile=trim(filename)


   ! would be pssoible in the future to read all J value also
   if (Atoms(nmet)%ptr_atom%initial_solution.ne."OLD_POPULATIONS" &
      .and.Atoms(nmet)%ptr_atom%initial_solution.ne."LTE_POPULATIONS"&
      .and.Atoms(nmet)%ptr_atom%initial_solution.ne."ZERO_RADIATION"&
		.and.Atoms(nmet)%ptr_atom%initial_solution.ne."CSWITCH")then
	     write(*,*) "Initial solution ", Atoms(nmet)%ptr_atom%initial_solution,&
      " unkown!"
     write(*,*) "Exiting..."
     stop
   end if
   !!Atoms(nmet)%ptr_atom%dataFile = trim(popsFile) ! now the file atomID.fits.gz contains
                                               ! all informations including populations.

   !Active atoms are treated in NLTE
   if (adjustl(actionKey).eq."ACTIVE") then
     Atoms(nmet)%ptr_atom%active=.true.
     !write(*,*) "atom is active"
     Nactiveatoms = Nactiveatoms+1
   else
     Atoms(nmet)%ptr_atom%active=.false.
     NpassiveAtoms = NpassiveAtoms+1
     !write(*,*) "atom is passive"
   end if
!   ! just opoen the model to check that Hydrogen is the first model
!   !
   open(unit=unit+nmet,file=trim(mcfost_utils)//TRIM(path_to_atoms)//trim(filename),status="old")
    CALL getnextline(unit+nmet, COMMENT_CHAR, FormatLine, inputline, Nread)
    read(inputline,*) IDread
    if (nmet.eq.1 .and. IDread.ne."H ") then
      write(*,*) "First atomic model read has to be Hydrogen"
      write(*,*) "Exting..."
      stop
    end if
   close(unit+nmet)
   ! check duplicates
   IDread(2:2) = to_lower(IDread(2:2))
   do mmet = 1,nmet-1 !compare the actual (IDread) with previous
    !write(*,*) mmet, nmet
    if (Atoms(mmet)%ptr_atom%ID.eq.IDread) then
     write(*,*) "Already read a model for this atom ", IDread
     write(*,*) "exiting..."
     stop
    end if
   end do

   ! read and fill atomic structures saved in atoms
   ! create an alias for Hydrogen
   CALL readModelAtom(unit+nmet, Atoms(nmet)%ptr_atom, &
        trim(mcfost_utils)//TRIM(path_to_atoms)//trim(filename))
   !write(*,*) "IS ACTIVE = ", Atoms(nmet)%active
   !Atoms(nmet)%ptr_atom = atom
   !CALL freeAtom(atom)
   !write(*,*) nmet, Atoms(nmet)%ptr_atom%ID
   !if (nmet>1) write(*,*) nmet-1, Atoms(nmet-1)%ptr_atom%ID
  end do
  close(unit)
  
  if (lfix_backgrnd_opac) then
  	do nmet=1, Natom
  		if (atoms(nmet)%ptr_atom%initial_solution.eq."OLD_POPULATIONS") then
  			call warning(" Using previous NLTE populations with fixed bacgkrnd!!!!")
  			write(*,*) " ----->n check consistency of the result"
  			exit
  		endif
  	enddo
  endif

  ! Alias to the most importent one
  Hydrogen=>Atoms(1)%ptr_atom
  if (.not.associated(Hydrogen, Atoms(1)%ptr_atom)) CALL Error(" Hydrogen alias not associated to atomic model!")

  ! Aliases to active atoms
  nact = 0; npass = 0
  if (Nactiveatoms > 0) then
    allocate(ActiveAtoms(Nactiveatoms))
  end if
  if (Natom - Nactiveatoms == Npassiveatoms) then
   allocate(PassiveAtoms(Npassiveatoms))
  else
   write(*,*) "Error, number of passive atoms is not N Nactive"
   stop
  end if

  ! keep a duplicate in Elements
  write(*,*) "order#     ID   periodic-table#    ACTIVE    #lines   #continua"
  do nmet=1,Natom
   write(*,*) nmet, Atoms(nmet)%ptr_atom%ID, &
    Atoms(nmet)%ptr_atom%periodic_table, Atoms(nmet)%ptr_atom%active, &
    Atoms(nmet)%ptr_atom%Nline, Atoms(nmet)%ptr_atom%Ncont
   ! create alias in Elements for elements that have
   ! a model atom. It means all elements here.
   !atom => Atoms(nmet)
   Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%model => Atoms(nmet)%ptr_atom
   
   if (.not.associated(Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%model, &
    Atoms(nmet)%ptr_atom)) then
		write(*,*) Atoms(nmet)%ptr_atom%id 
    	CALL error(" Elemental model not associated to atomic model!")
    endif
    
	!Check Nstage
	if (Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%Nstage /= maxval(Atoms(nmet)%ptr_atom%stage) + 1) then
		write(*,*) Atoms(nmet)%ptr_atom%id, maxval(Atoms(nmet)%ptr_atom%stage) + 1
		write(*,*) "Ns pf = ", Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%Nstage
		call error("Model has more ionisation stages than the one in the partition function!")
	endif
    
   if (Atoms(nmet)%ptr_atom%periodic_table.eq.2)  then
           NULLIFY(Helium) !Because, it is associated to an Elem by default
           Helium => Atoms(nmet)%ptr_atom
     if (.not.associated(Helium,Atoms(nmet)%ptr_atom)) &
      CALL Warning(" Helium alias is not associated to an atomic model!")
   end if
   
   if (allocated(ActiveAtoms)) then
     if (Atoms(nmet)%ptr_atom%active) then
      nact = nact + 1 !got the next index of active atoms
      ActiveAtoms(nact)%ptr_atom => Atoms(nmet)%ptr_atom
      Atoms(nmet)%ptr_atom%activeindex = nact
!       atom2 => ActiveAtoms(nact)
!       atom2 = atom
     end if
   end if
   
   if (allocated(PassiveAtoms)) then
     if (.not.Atoms(nmet)%ptr_atom%active) then
      npass = npass+1
!       atom2 => PassiveAtoms(npass)
!       atom2 = atom
      PassiveAtoms(npass)%ptr_atom => Atoms(nmet)%ptr_atom
     end if
   end if
  end do
  
  min_Resol = 1d30
  do nmet=1,Natom
  	do k=1,n_cells
  		if (icompute_atomRT(k)>0) then
  			min_resol = min(Atoms(nmet)%ptr_atom%vbroad(k), min_resol)
  		endif
  	enddo
  enddo
  hv = 0.46 * real(min_resol) * 1e-3
  
  if (art_hv > 0.0) then
	hv = art_hv
  endif
  write(*,*) "resol(km/s)", hv, " min(vD) (km/s) = ", min_resol * 1d-3


  !Move after LTEpops for first estimates of damping
  !line wave grid define here to have the max damping
  write(*,*) " Generating sub wavelength grid and lines boundary for all atoms..."
  do nmet=1, Natom  
	atom => Atoms(nmet)%ptr_atom
   !!write(*,*) "ID:", Atoms(nmet)%ptr_atom%ID
   do kr=1, Atoms(nmet)%ptr_atom%Nline

		maxvel = Atoms(nmet)%ptr_atom%lines(kr)%qwing * maxval(atom%vbroad)

!      ic = find_continuum(atom,atom%lines(kr)%i)
! 
!      max_adamp = 1d100
!      maxvel = 0
!      adamp = 0.0
!      imax = 1
!      !!write(*,*) " line @ ", Atoms(nmet)%ptr_atom%lines(kr)%lambda0,'nm'
!      do k=1, n_cells
!        if (icompute_atomRT(k)>0) then
!        
!         vel = Atoms(nmet)%ptr_atom%vbroad(k)
! 
!         if (Atoms(nmet)%ptr_atom%lines(kr)%voigt) then 
!          !no damping if Gaussian
!          CALL Damping(k, Atoms(nmet)%ptr_atom, kr, adamp)
!          !-> LTE populations not known, we cannot use damping here
!          !Because some of the damping are zero if prop to populations
!          !Typically there are only Stark and Natural depending on electrons and lines.
!          !If only Adamp its to short
!          !adamp = atom%lines(kr)%Aji * (NM_TO_M*atom%lines(kr)%lambda0) / (4.*pi) / atom%vbroad(k)
!          !max_adamp = max(max_adamp, adamp)
!          max_adamp = min(max_adamp, adamp) ! min actually
!          
!          epsilon_l_max = 1.0/pi/adamp
! 
!          if (n_eff(atom%Rydberg, atom%E(ic), atom%E(atom%lines(kr)%i), atom%stage(ic)) <= 3.0) then 
!           eps = 1e-4
!          else
!           eps = epsilon
!          endif
! !-> with vD         
! !         maxvel = max(maxvel, &
! !          vel * sqrt(abs(-log(epsilon))), vel * sqrt(min(adamp, 0.9*1.0/pi/epsilon)/pi/min(0.9*epsilon_l_max, epsilon) - min(adamp,0.9*1.0/pi/epsilon)**2) )
! !-> only from damping
! ! 		 maxvel = max(maxvel, vel * sqrt(min(adamp, 0.9*1.0/pi/eps)/pi/min(0.9*epsilon_l_max, eps) - min(adamp,0.9*1.0/pi/eps)**2))	 
! 		else
! 		 !-> exact
! ! 		 maxvel = max(maxvel, vel * sqrt(abs(-log(eps))))
!         endif
!   		
!   		maxvel = max(maxvel, vel)      
! 
! 
!        endif
!        
!      enddo
    !! maxvel = maxvel_atom_transfer * 1e3 * ( 1e-3 * vel/1.0)

! 	 maxvel = maxval(atom%vbroad)
     !write(*,*) "maxvel for line ", kr, maxvel * 1e-3
     !!-> group of lines, linear in v
     call define_local_profile_grid (Atoms(nmet)%ptr_atom%lines(kr))

 
 !-> depends if we interpolate profile on finer grid !   
     !!-> linear
!      CALL make_sub_wavelength_grid_line_lin(Atoms(nmet)%ptr_atom%lines(kr),&
!                                         maxval(Atoms(nmet)%ptr_atom%vbroad), max_adamp)  
     !!-> logarithmic
!      CALL make_sub_wavelength_grid_line(Atoms(nmet)%ptr_atom%lines(kr),&
!                                         maxval(Atoms(nmet)%ptr_atom%vbroad), max_adamp)

   
   enddo !over lines
  
  enddo !over atoms
  write(*,*) "..done"
  Max_ionisation_stage = get_max_nstage()


!   write(*,*) Atoms(1)%ptr_Atom%ID,ActiveAtoms(1)%ptr_Atom%ID
!   write(*,*) Hydrogen%ID, elements(1)%ptr_elem%model%ID
!   stop

!  NULLIFY(atom, atom2)

!   do nmet=1, Npassiveatoms
!    write(*,*) nmet, Atoms(nmet)%ptr_atom%ID, PassiveAtoms(nmet)%ptr_atom%ID
!   end do

!   write(*,*) "  -> writing lines individual wavelength grid"
!   !write lines grid
!   CALL write_lines_grid()

  RETURN
  END SUBROUTINE readAtomicModels
  
 SUBROUTINE search_cont_lambdamax (cont, Rinf, Z, Ej, Ei)
 !Search the lambdamax = first lambda for which Gaunt < 0
 !Because I do not have an extrapolation routine, there is no need to extrapolate beyond
 !lambdamax if gaunt is negative.
 
 !See Gaunt_bf in Hydrogen.f90
  real(kind=dp), intent(in) :: Ej, Ei, Rinf
  integer, intent(in) :: Z
  type (AtomicContinuum), intent(inout) :: cont
  real, parameter :: l1 = 1e3, dlam = 1.0 !nm !0.01
  real(kind=dp) :: n_eff,  u ! = n_Eff**2 * eps = hnu/Z/Z/E_RYDBERG - 1
  real(kind=dp) :: Gaunt_bf, x
  
  x = l1 * cont%lambda0
  cont%lambdamax = cont%lambda0
  n_eff = Z * sqrt(Rinf / (Ej-Ei))
  Gaunt_bf = 1.0_dp
  do while(cont%lambdamax < x)
               
   u = n_eff**2 * HPLANCK*CLIGHT / (NM_TO_M * cont%lambdamax) / Z*Z / E_RYDBERG - 1
  
   Gaunt_bf = 1d0 + 0.1728 * (n_eff**(-2./3.)) * (u+1d0)**(-2./3.) * (u-1d0) &
             - 0.0496*(n_eff**(-4./3.)) * (u+1d0)**(-4./3.) * (u*u + 4./3. * u + 1d0)

   if (Gaunt_bf < 0) exit

   cont%lambdamax = cont%lambdamax + dlam

  enddo

  cont%lambdamax = min(cont%lambdamax,1d5)
  write(*,*) cont%atom%ID, " cont, n=",real(n_eff), cont%j, cont%i, real(cont%lambda0),"nm"
  write(*,*) " -> pseudo-continuum: ", real(cont%lambdamax), " nm: frac=", real(cont%lambdamax/cont%lambda0)
             

 RETURN
 END SUBROUTINE search_cont_lambdamax


 !But why, a cswitch per atom ? It is going at the same speed for all atoms right ?
 function maxval_cswitch_atoms ()
 !for all atoms, check the maximum value of the cswitch
 	integer :: n
 	real(kind=dp) ::  maxval_cswitch_atoms
 	
	maxval_cswitch_atoms = 1.0_dp
 	do n=1, Natom
 	
 		maxval_cswitch_atoms = max(maxval_cswitch_atoms, atoms(n)%ptr_atom%cswitch)
 	
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
 	
 		if (activeatoms(n)%ptr_atom%cswitch > 1.0) then
 			activeatoms(n)%ptr_atom%cswitch = max(1.0_dp, min(activeatoms(n)%ptr_atom%cswitch, activeatoms(n)%ptr_atom%cswitch/cswitch_down_scaling_factor))
 			if (print_message == .false.) then
 				print_message = .true.
 				new_cs = activeatoms(n)%ptr_atom%cswitch
 			endif
 		endif
 	enddo
 	
 	if (print_message) write(*,'(" cswitch for next iteration: "(1ES17.8E3))') new_cs

 
 return
 end subroutine adjust_cswitch_atoms

END MODULE readatom
