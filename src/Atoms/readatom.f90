MODULE readatom

  use atom_type, only : AtomicLine, AtomicContinuum, AtomType, Element, determinate
  use atmos_type, only : atmos, Nelem, Hydrogen, Helium
  use zeeman, only : Lande_eff, ZeemanMultiplet
  use getlambda
  use constant
  use uplow
  use getline
  use barklem, only : getBarklem
  use writeatom, only : readPops
  use collision, only : read_collisions

  !$ use omp_lib

  !MCFOST's originals
  use messages
  use mcfost_env, only : mcfost_utils ! convert from the relative location of atomic data
                                      ! to mcfost's environnement folders.

  IMPLICIT NONE

  character, parameter :: COMMENT_CHAR="#"
  character(len=*), parameter :: ATOMS_INPUT = "./atoms.input"!"/Atoms/atoms.input"
  character(len=*), parameter :: path_to_atoms = "/Atoms/"
  !real, parameter :: MAX_ABUND_ERROR=0.001


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
    character(len=20) :: shapeChar, symmChar, optionChar, vdWChar, nuDepChar
    character(len=2) :: IDread
    real(8) :: C1, vDoppler, f, lambdaji
    real(8) :: lambdamin, geff!, c_sum
    EOF = 0
    !write(*,*) "Atom is part of the active atoms ?", atom%active

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
     if (atmos%Elements(nll)%ptr_elem%ID.eq.atom%ID) then
      write(*,*) "Abundance of atom ",atom%ID,": A =",atmos%Elements(nll)%ptr_elem%Abund
      if (atmos%Elements(nll)%ptr_elem%abundance_set.eqv..true.) then
        atom%periodic_table=nll
        atom%Abund=atmos%Elements(nll)%ptr_elem%Abund
        atom%weight=atmos%Elements(nll)%ptr_elem%weight
        match=.true.
      end if
      exit
     end if
    end do

    if (match .neqv..true.) then
     write(*,*) "Error no abundance found for atom ", atom_file, atom%ID
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
    allocate(atom%nstar(atom%Nlevel,atmos%Nspace))

    atom%Ntr = atom%Nline + atom%Ncont
    allocate(atom%at(atom%Ntr))

    atom%nstar(:,:) = 0d0

    do i=1,atom%Nlevel
     CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
     read(inputline,*) atom%E(i), atom%g(i), atom%label(i), &
                       atom%stage(i), levelNumber(i)
     atom%E(i) = atom%E(i) * HPLANCK*CLIGHT / (CM_TO_M)
!      write(*,*) "E(cm^-1) = ", atom%E(i) * JOULE_TO_EV, &
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
    !allocate(atom%ntotal(atmos%Nspace))
    allocate(atom%vbroad(atmos%Nspace))
    !write(*,*)
    VDoppler = KBOLTZMANN/(AMU * atom%weight) * 8d0/PI !* 2d0!m/s
    !        = vtherm/atom%weight

    atom%vbroad = dsqrt(vtherm/atom%weight * atmos%T)! + atmos%vturb**2) !vturb in m/s
    !atom%ntotal = atom%Abund * atmos%nHtot

    VDoppler = dsqrt(Vdoppler*maxval(atmos%T))! + maxval(atmos%vturb)**2)
!     write(*,*) Vdoppler, dsqrt(Vtherm*maxval(atmos%T)/atom%weight + maxval(atmos%vturb)**2)
    !Now read all bound-bound transitions
    allocate(atom%lines(atom%Nline))

    do kr=1,atom%Nline
     atom%lines(kr)%atom => atom
     !!atom%lines(kr)%lcontrib_to_opac=.true. !init
     atom%at(kr)%trtype = atom%lines(kr)%trtype; atom%at(kr)%ik=kr
     atom%at(kr)%lcontrib_to_opac=.true.; atom%Ntr_line = atom%Nline

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
     if ((PRT_SOLUTION=="NO_STOKES").or.(.not.atmos%magnetized)) atom%lines(kr)%ZeemanPattern = 0
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

     write(*,*) "Reading line #", kr, 1d9 * (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i)), 'nm'

      !therefore, the first level is 1 (C=0), the second 2 (C=1) etc
      !Lymann series: 2->1, 3->1
      atom%lines(kr)%i = min(i,j)
      atom%lines(kr)%j = max(i,j)

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
      atom%lines(kr)%polarizable = (atmos%magnetized) .and. &
      								(atom%lines(kr)%g_lande_eff > -99) .and. &
      			(abs(atom%qJ(atom%lines(kr)%i) - atom%qJ(atom%lines(kr)%j)) <= 1.)
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
      
      
      atom%lines(kr)%gij = atom%lines(kr)%Bji / atom%lines(kr)%Bij !gi/gj
      atom%lines(kr)%twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji

      write(*,*) " ->", " Aji (1e7 s^-1) = ", atom%lines(kr)%Aji/1d7,&
        "Grad (1e7 s^-1) = ", atom%lines(kr)%Grad/1d7, &
        "gj = ", atom%g(j)," gi = ",  atom%g(i)

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
     
 !!force Gaussian for test
     !atom%lines(kr)%Voigt = .false.

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
     if (trim(symmChar).eq."ASYMM") then
      atom%lines(kr)%symmetric = .false.
     else
      atom%lines(kr)%symmetric = .true.
      !write(*,*) "Symmetric line profile"
     end if
     atom%lines(kr)%symmetric = .false.


!      if (atom%active) then !Should do it for passive atoms too
!      if (atmos%Magnetized) then !.or.line%scattpol ...
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
    ! ----------------------------------------- !
    allocate(atom%continua(atom%Ncont))
    do kr=1,atom%Ncont
     atom%continua(kr)%isotope_frac=1.
     atom%continua(kr)%atom => atom
     !!atom%continua(kr)%lcontrib_to_opac=.true.
     atom%at(kr)%lcontrib_to_opac=.true.
     atom%at(atom%Nline+kr)%trtype = atom%continua(kr)%trtype; atom%at(kr+atom%Nline)%ik=kr

     CALL getnextline(atomunit, COMMENT_CHAR, &
          FormatLine, inputline, Nread)
     read(inputline, *) j, i, atom%continua(kr)%alpha0,&
      atom%continua(kr)%Nlambda, nuDepChar, atom%continua(kr)%lambdamin
     j = j + 1
     i = i + 1
     write(*,*) "Reading continuum #", kr, atom%continua(kr)%lambdamin, "nm", &
     					1d9 * (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i)), "nm"

     !because in C indexing starts at 0, but starts at 1 in fortran
     atom%continua(kr)%j = max(i,j)
     atom%continua(kr)%i = min(i,j)
     lambdaji = (HPLANCK*CLIGHT)/ (atom%E(j)-atom%E(i))
     atom%continua(kr)%lambda0 = lambdaji/NM_TO_M !nm
     !lambdamin=atom%continua(kr)%lambdamin

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
     else if (trim(nuDepChar).eq."HYDROGENIC") then
       atom%continua(kr)%hydrogenic=.true.
       if (atom%continua(kr)%lambdamin>=atom%continua(kr)%lambda0) then
        write(*,*) "Minimum wavelength for continuum is "&
          "larger than continuum edge."
        write(*,*) "exiting..."
        stop
       end if
       !input and fill wavelength grid for HYRDROGENIC
       ! Matters only when the complete wavelength dependence
       ! is not read from file (HYDROGENIC).
       ! %lambda allocated inside the routines.
       CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin)
       !write(*,*) "lambdacontmin = ", atom%continua(kr)%lambda(1), &
       !" lambdacontmax = ", atom%continua(kr)%lambda(atom%continua(kr)%Nlambda)
     else
      write(*,*) "Invalid value for continuum chromatic ",&
         "dependence."
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

    !now compute wavelengths grid for
    ! each line
    !Unlike RH, I do it even for passive atoms
  do kr=1,atom%Nline !line%lambda allocated inside
     !-> Actually with this grid, all lines of an atom have the same grid
     !because it depends only on vD and v_char which depends on the atom and on the model.
     !This is because the Number of core/wing points are fixed.
     !CALL make_sub_wavelength_grid_line(atom%lines(kr),Vdoppler) !MAXVAL(atom%vbroad)
     CALL make_sub_wavelength_grid_line_lin(atom%lines(kr),Vdoppler)
     !CALL make_sub_wavelength_grid_asymm(atom%lines(kr),Vdoppler)
  end do
   !Now even for passive atoms we write atomic data.
   ! Unlike RH, all data are in the same fits file.
    if (atom%ID(2:2) .eq." ") then
      atom%dataFile = atom%ID(1:1)//".fits.gz" !.fits to be updated, .gz not
    else
      atom%dataFile = atom%ID(1:2)//".fits.gz"
    end if
    !! done at writing if we do not store them by iterations, but only at the end of the NLTE loop
    !!CALL create_pops_file(atom)
    !!write(*,*) "Populations file for writing: .",trim(atom%dataFile),"."

   if (atom%active) then

    !Not implemented, futur removal in atomic file
    if (atom%Npfr.gt.0) then
     write(*,*) "PFR not implemented yet, do not write file"
    end if
    if (atmos%XRD) then
     write(*,*) "Cross redistribution ", &
                "not implemented yet."
    end if

    ! reading collision rates
    call read_collisions(atomunit, atom)

    ! allocate some space
    atom%colunit = atom%periodic_table*2 + 1

    !!allocate(atom%C(atom%Nlevel*atom%Nlevel,atmos%Nspace))
    !!now Collision matrix is constructed cell by cell, therefore allocated elsewhere
    allocate(atom%n(atom%Nlevel,atmos%Nspace))
    atom%n = 0d0
    atom%NLTEpops = .false.
	if (atom%dataFile.ne."" .and. atom%initial_solution .eq. "OLD_POPULATIONS") then
	   write(*,*) " -> Reading populations from file..."
       CALL readPops(atom)
       atom%NLTEpops = .true.
       write(*,*) " min/max pops for each level:"
       do kr=1,atom%Nlevel
        write(*,*) "    ", kr, ">>", minval(atom%n(kr,:)), maxval(atom%n(kr,:))
       enddo
    end if
   else !not active
    if (atom%dataFile.ne.""  .and. atom%initial_solution .eq. "OLD_POPULATIONS") then
       allocate(atom%n(atom%Nlevel,atmos%Nspace)) !not allocated if passive, n->nstar
       CALL readPops(atom)
       atom%NLTEpops = .true.
    else !pure passive without nlte pops from previous run
     atom%NLTEpops=.false.
     atom%n => atom%nstar !initialised to zero
     !atom%n is an alias for nstar in this case
    end if
   end if !end is active

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
  integer, intent(in) :: unit
  character(len=MAX_LENGTH) :: inputline
  character(len=15) :: FormatLine
  character(len=MAX_LENGTH) :: popsfile, filename
  character(len=MAX_KEYWORD_SIZE) :: actionKey, popsKey
  character(len=2) :: IDread
!   type (AtomType), target :: atom

  !create formatline
  !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
  !write(FormatLine,'("("I3,A")")') MAX_LENGTH,"A"
  write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH


  atmos%Nactiveatoms = 0
  atmos%Npassiveatoms = 0
  open(unit=unit,file=TRIM(ATOMS_INPUT), status="old")!mcfost_utils)//TRIM(ATOMS_INPUT)

  !get number of atomic models to read
  CALL getnextline(unit, COMMENT_CHAR,FormatLine, &
       inputline, Nread)
  read(inputline,*) atmos%Natom
  write(*,*) "Reading ", atmos%Natom, " species"
  !Allocate sapace for Natom in atmos%Atoms
  allocate(atmos%Atoms(atmos%Natom))

  do nmet = 1, atmos%Natom
   CALL getnextline(unit, COMMENT_CHAR, FormatLine, inputline, Nread)

   allocate(atmos%Atoms(nmet)%ptr_atom)

   read(inputline,'(1A28, 1A7, 1A22, 1A20)') filename, actionKey, popsKey, popsFile
   !write(*,*) ".",trim(filename),"."
   !write(*,*) ".",adjustl(actionKey),"."
   !write(*,*) ".",adjustl(popsKey),"."
   !write(*,*) ".",trim(popsFile),"."

   atmos%Atoms(nmet)%ptr_atom%initial_solution=adjustl(popsKey)
   atmos%Atoms(nmet)%ptr_atom%inputFile=trim(filename)


   ! would be pssoible in the future to read all J value also
   if (atmos%Atoms(nmet)%ptr_atom%initial_solution.ne."OLD_POPULATIONS" &
      .and.atmos%Atoms(nmet)%ptr_atom%initial_solution.ne."LTE_POPULATIONS"&
      .and.atmos%Atoms(nmet)%ptr_atom%initial_solution.ne."ZERO_RADIATION"&
       .and.atmos%Atoms(nmet)%ptr_atom%initial_solution.ne."SOBOLEV") then
     write(*,*) "Initial solution ", atmos%Atoms(nmet)%ptr_atom%initial_solution,&
      " unkown!"
     write(*,*) "Exiting..."
     stop
   end if
   atmos%Atoms(nmet)%ptr_atom%dataFile = trim(popsFile) ! now the file atomID.fits.gz contains
                                               ! all informations including populations.

   !Active atoms are treated in NLTE
   if (adjustl(actionKey).eq."ACTIVE") then
     atmos%Atoms(nmet)%ptr_atom%active=.true.
     !write(*,*) "atom is active"
     atmos%Nactiveatoms = atmos%Nactiveatoms+1
   else
     atmos%Atoms(nmet)%ptr_atom%active=.false.
     atmos%NpassiveAtoms = atmos%NpassiveAtoms+1
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
    if (atmos%Atoms(nmet)%ptr_atom%ID.eq.IDread) then
     write(*,*) "Already read a model for this atom ", IDread
     write(*,*) "exiting..."
     stop
    end if
   end do

   ! read and fill atomic structures saved in atmos%atoms
   ! create an alias for Hydrogen
   CALL readModelAtom(unit+nmet, atmos%Atoms(nmet)%ptr_atom, &
        trim(mcfost_utils)//TRIM(path_to_atoms)//trim(filename))
   !write(*,*) "IS ACTIVE = ", atmos%Atoms(nmet)%active
   !atmos%Atoms(nmet)%ptr_atom = atom
   !CALL freeAtom(atom)
   !write(*,*) nmet, atmos%Atoms(nmet)%ptr_atom%ID
   !if (nmet>1) write(*,*) nmet-1, atmos%Atoms(nmet-1)%ptr_atom%ID
  end do
  close(unit)

  ! Alias to the most importent one
  Hydrogen=>atmos%Atoms(1)%ptr_atom
  if (.not.associated(Hydrogen, atmos%Atoms(1)%ptr_atom)) CALL Error(" Hydrogen alias not associated to atomic model!")

  ! Aliases to active atoms
  nact = 0; npass = 0
  if (atmos%Nactiveatoms > 0) then
    allocate(atmos%ActiveAtoms(atmos%Nactiveatoms))
  end if
  if (atmos%Natom - atmos%Nactiveatoms == atmos%Npassiveatoms) then
   allocate(atmos%PassiveAtoms(atmos%Npassiveatoms))
  else
   write(*,*) "Error, number of passive atoms is not N Nactive"
   stop
  end if

  ! keep a duplicate in Elements
  write(*,*) "order#     ID   periodic-table#    ACTIVE    #lines   #continua"
  do nmet=1,atmos%Natom
   write(*,*) nmet, atmos%Atoms(nmet)%ptr_atom%ID, &
    atmos%Atoms(nmet)%ptr_atom%periodic_table, atmos%Atoms(nmet)%ptr_atom%active, &
    atmos%Atoms(nmet)%ptr_atom%Nline, atmos%Atoms(nmet)%ptr_atom%Ncont
   ! create alias in atmos%Elements for elements that have
   ! a model atom. It means all elements here.
   !atom => atmos%Atoms(nmet)
   atmos%Elements(atmos%Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%model &
         => atmos%Atoms(nmet)%ptr_atom
   if (.not.associated(atmos%Elements(atmos%Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%model, &
    atmos%Atoms(nmet)%ptr_atom)) CALL Warning(" Elemental model not associated to atomic model!")
   if (atmos%Atoms(nmet)%ptr_atom%periodic_table.eq.2)  then
           NULLIFY(Helium) !Because, it is associated to an Elem by default
           Helium => atmos%Atoms(nmet)%ptr_atom
     if (.not.associated(Helium,atmos%Atoms(nmet)%ptr_atom)) &
      CALL Warning(" Helium alias not associated to atomic model!")
   end if
   if (allocated(atmos%ActiveAtoms)) then
     if (atmos%Atoms(nmet)%ptr_atom%active) then
      nact = nact + 1 !got the next index of active atoms
      atmos%ActiveAtoms(nact)%ptr_atom => atmos%Atoms(nmet)%ptr_atom
      atmos%Atoms(nmet)%ptr_atom%activeindex = nact
!       atom2 => atmos%ActiveAtoms(nact)
!       atom2 = atom
     end if
   end if
   if (allocated(atmos%PassiveAtoms)) then
     if (.not.atmos%Atoms(nmet)%ptr_atom%active) then
      npass = npass+1
!       atom2 => atmos%PassiveAtoms(npass)
!       atom2 = atom
      atmos%PassiveAtoms(npass)%ptr_atom => atmos%Atoms(nmet)%ptr_atom
     end if
   end if
  end do

!   write(*,*) atmos%Atoms(1)%ptr_Atom%ID,atmos%ActiveAtoms(1)%ptr_Atom%ID
!   write(*,*) Hydrogen%ID, atmos%elements(1)%ptr_elem%model%ID
!   stop

!  NULLIFY(atom, atom2)

!   do nmet=1, atmos%Npassiveatoms
!    write(*,*) nmet, atmos%Atoms(nmet)%ptr_atom%ID, atmos%PassiveAtoms(nmet)%ptr_atom%ID
!   end do

!   write(*,*) "  -> writing lines individual wavelength grid"
!   !write lines grid
!   CALL write_lines_grid()

  RETURN
  END SUBROUTINE readAtomicModels


END MODULE readatom
