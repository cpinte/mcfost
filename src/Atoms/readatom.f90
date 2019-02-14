MODULE readatom

  use atom_type, only : AtomicLine, AtomicContinuum, AtomType, Element, determinate
  use atmos_type, only : atmos, Nelem, Hydrogen, Helium
  use zeeman, only : Lande_eff
  use getlambda
  use constant
  use uplow
  use getline
  use barklem, only : getBarklem

  !$ use omp_lib

  !MCFOST's originals
  use mcfost_env, only : mcfost_utils ! convert from the relative location of atomic data
                                      ! to mcfost's environnement folders.

  IMPLICIT NONE

  character, parameter :: COMMENT_CHAR="#"
  character(len=*), parameter :: ATOMS_INPUT = "/Atoms/atoms.input"
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
    integer :: kr, k, la, ftell !ftell here for ifort ??
    type (AtomType), intent(inout) :: atom
    type (AtomicLine) :: line, line1
    type (AtomicContinuum) :: AtomicContinuum
    character(len=*), intent(in) :: atom_file
    character(len=MAX_LENGTH) :: inputline, FormatLine
    integer :: Nread, i,j, EOF, nll, nc
    real, allocatable, dimension(:) :: levelNumber
    logical :: Debeye, match, res=.false.
    logical, dimension(:), allocatable :: determined, parse_label
    character(len=20) :: shapeChar, symmChar, optionChar, vdWChar, nuDepChar
    character(len=2) :: IDread
    real(8) :: C1, vtherm, f, lambdaji
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
    !write(*,*) "AtomID = ", atom%ID

    match = .false.
    do nll=1,Nelem
     if (atmos%Elements(nll)%ID.eq.atom%ID) then
      !write(*,*) "Abundance found for atom ",atom%ID,atmos%Elements(nll)%Abund
      if (atmos%Elements(nll)%abundance_set.eqv..true.) then
        atom%periodic_table=nll
        atom%Abund=atmos%Elements(nll)%Abund
        atom%weight=atmos%Elements(nll)%weight
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
    read(inputline,*) atom%Nlevel, atom%Nline, atom%Ncont, atom%Nfixed
    write(*,*) "Nlevel=",atom%Nlevel," Nline=",atom%Nline,&
              " Ncont=", atom%Ncont, " Nfixed=", atom%Nfixed

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
    atom%nstar = 0d0

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
    if (atom%stage(atom%Nlevel).ne.&
           atom%stage(atom%Nlevel-1)+1) then
     write(*,*) "Atomic model does not have an overlying continuum"
     write(*,*) "exiting..."
     stop
    end if

    allocate(atom%ntotal(atmos%Nspace))
    allocate(atom%vbroad(atmos%Nspace))
    vtherm = 2.*KBOLTZMANN/(AMU * atom%weight) !m/s

    atom%vbroad = dsqrt(vtherm*atmos%T + atmos%vturb**2) !vturb in m/s
    atom%ntotal = atom%Abund * atmos%nHtot
!     write(*,*) atom%ID, " maxVD(km/s)=", maxval(atom%vbroad)/1d3,&
!                          " minVD(km/s)=", minval(atom%vbroad,mask=atom%vbroad>0)/1d3
    !! DO NOT USE i as a loop index here !!

    !Now read all bound-bound transitions
    allocate(atom%lines(atom%Nline))

    do kr=1,atom%Nline
     !atom%lines(kr)%atom => atom !hmm work this?
     atom%lines(kr)%isotope_frac = 1.
     atom%lines(kr)%g_lande_eff = -99.0
     !atom%lines(kr)%trtype="ATOMIC_LINE"
     CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
     Nread = len(trim(inputline)) ! because, if blanck
           ! beyond cStark it will be interpreted
           ! has "additional geff", but its not.
           !
     if (Nread.eq.112) then
         read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
         symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
         atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
         atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
         atom%lines(kr)%Grad, atom%lines(kr)%cStark
     else if (Nread.gt.112) then
       write(*,*) "Read aditional g_lande_eff for that line"
       read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
       symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
       atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
       atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
       atom%lines(kr)%Grad, atom%lines(kr)%cStark, &
       atom%lines(kr)%g_Lande_eff
     else
       write(*,*) "Error reading line"
       write(*,*) "exiting.."
       stop
     end if
      i = i + 1
      j = j + 1 !because in C, indexing starts at 0, but starts at 1 in fortran
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
      end if
      if (.not.parse_label(atom%lines(kr)%j)) then
      CALL determinate(atom%label(atom%lines(kr)%j),&
       atom%g(atom%lines(kr)%j),&
       atom%qS(atom%lines(kr)%j),&
       atom%Lorbit(atom%lines(kr)%j),&
       atom%qJ(atom%lines(kr)%j), &
       determined(atom%lines(kr)%j))
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

      if (((determined(atom%lines(kr)%j)).and.&
          (determined(atom%lines(kr)%i))).and. &
           atom%lines(kr)%g_Lande_eff.eq.-99) then
       ! do not compute geff if term is not
       ! determined or if the geff is read from file
       ! ie if g_lande_eff > -99
       ! fill lande_g_factor if determined and not given
       ! for b-b transitions
       CALL Lande_eff(atom, kr)
       !!testing
       !!write(*,*) wKul(atom, kr, 0)**2
       !!write(*,*) wKul(atom, kr, 1)**2
       !!write(*,*) wKul(atom, kr, 2)**2
       !!stop
       !write(*,*) "geff = ", atom%lines(kr)%g_lande_eff
      end if
      ! oscillator strength saved
      atom%lines(kr)%fosc = f
      lambdaji = (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i))
      atom%lines(kr)%Aji = C1 / (lambdaji**2) * (atom%g(i) / atom%g(j)) * f
      atom%lines(kr)%Bji = (lambdaji**3) / (2.0 * HPLANCK * CLIGHT) &
                          *atom%lines(kr)%Aji
      atom%lines(kr)%Bij = (atom%g(j) / atom%g(i)) * atom%lines(kr)%Bji
      atom%lines(kr)%lambda0 = lambdaji / NM_TO_M

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
      line%Voigt = .true.
        !write(*,*) "Using Voigt profile for that line"
        !allocate(atom%lines(kr)%c_shift(1))
        !allocate(atom%lines(kr)%c_fraction(1))
        !atom%lines(kr)%Ncomponent = 1
        !atom%lines(kr)%c_shift(1) = 0.0
        !atom%lines(kr)%c_fraction(1) = 1.
     end if

     !Now parse Broedening recipe
     if (trim(vdWChar).eq."PARAMTR") then
       atom%lines(kr)%vdWaals = "RIDDER_RENSBERGEN"
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

     if (trim(symmChar).eq."ASYMM") then !futur deprecation, but atmosphere will always move
      !write(*,*) "Line has an asymmetric profile."
      atom%lines(kr)%symmetric = .false.
     else
      atom%lines(kr)%symmetric = .true.
      !write(*,*) "Symmetric line profile"
     end if
     atom%lines(kr)%polarizable = .false.

     if (atom%active) then
      allocate(atom%lines(kr)%Rij(atmos%Nspace))
      allocate(atom%lines(kr)%Rji(atmos%Nspace))

     if (atmos%Magnetized) then !.or.line%scattpol ...
      if (atom%lines(kr)%g_Lande_eff.gt.-99 .or. &
          determined(atom%lines(kr)%i) .and. &
          determined(atom%lines(kr)%j).and. &
          abs(atom%qJ(atom%lines(kr)%i) - &
           atom%qJ(atom%lines(kr)%j)).le.1.) then

!        if (atom%lines(kr)%Ncomponent.gt.1) then
!            !write(*,*) &
!            !"Cannot treat composite line with polar"
!            atom%lines(kr)%polarizable=.false.
!        else
         atom%lines(kr)%polarizable=.true.
!        end if
      end if
     else
      !write(*,*) "Treating line ",atom%lines(kr)%j,&
      !    "->",atom%lines(kr)%i,&
      !    " without polarization"
      atom%lines(kr)%polarizable=.false.
     end if !not mag
    end if ! end loop over active b-b transitions of atom
   end do !end loop over bound-bound transitions

    ! ----------------------------------------- !
    !starts reading bound-free transitions
    ! ----------------------------------------- !
    allocate(atom%continua(atom%Ncont))
    do kr=1,atom%Ncont
     atom%continua(kr)%isotope_frac=1.
     !atom%continua(kr)%trtype="ATOMIC_CONTINUUM"

     CALL getnextline(atomunit, COMMENT_CHAR, &
          FormatLine, inputline, Nread)
     read(inputline, *) j, i, atom%continua(kr)%alpha0,&
      atom%continua(kr)%Nlambda, nuDepChar, lambdamin
     j = j + 1
     i = i + 1
     !because in C indexing starts at 0, but starts at 1 in fortran
     atom%continua(kr)%j = max(i,j)
     atom%continua(kr)%i = min(i,j)
     lambdaji = (HPLANCK*CLIGHT)/ (atom%E(j)-atom%E(i))
     atom%continua(kr)%lambda0 = lambdaji/NM_TO_M !nm

     !write(*,*) "continuum ", atom%continua(kr)%j,&
     !    '->',atom%continua(kr)%i, " @",&
     !    lambdaji/NM_TO_M," nm : alpha0 = ", &
     !    atom%continua(kr)%alpha0, &
     !    " Nlambda = ", atom%continua(kr)%Nlambda, &
     !    " type = ", nuDepChar," lambda min = ",&
     !     lambdamin," nm"

     if (trim(nuDepChar).eq."EXPLICIT") then
      ! Nlambda set in atomic file
      allocate(atom%continua(kr)%lambda(atom%continua(kr)%Nlambda))
      allocate(atom%continua(kr)%alpha(atom%continua(kr)%Nlambda))
      atom%continua(kr)%hydrogenic=.false.
      ! rearanging them in increasing order
      ! because in the atomic file they are
      ! given in decreasing order !
      do la=atom%continua(kr)%Nlambda,1,-1
       CALL getnextline(atomunit, COMMENT_CHAR, &
             FormatLine, inputline, Nread)
       read(inputline,*) atom%continua(kr)%lambda(la), &
          atom%continua(kr)%alpha(la)
       ! though they are printed in decreasing order
       !write(*,*) "l = ",atom%continua(kr)%lambda(la), &
       !    " a = ", atom%continua(kr)%alpha(la)
      end do
      do la=2,atom%continua(kr)%Nlambda
        if (atom%continua(kr)%lambda(la).lt.&
           atom%continua(kr)%lambda(la-1)) then
          write(*,*) "continuum wavelength not monotonous"
          write(*,*) "exiting..."
          stop
        end if
      end do
     else if (trim(nuDepChar).eq."HYDROGENIC") then
       atom%continua(kr)%hydrogenic=.true.
       if (lambdamin.ge.atom%continua(kr)%lambda0) then
        write(*,*) "Minimum wavelength for continuum is "&
          "larger than continuum edge."
        write(*,*) "exiting..."
        stop
       end if
       !input and fill wavelength grid for HYRDROGENIC
       ! Matters only when the complete wavelength dependence
       ! is not read from file (HYDROGENIC).
       ! %lambda allocated inside the routines.
       ! %alpha not needed in this case.
       CALL make_sub_wavelength_grid_cont(atom%continua(kr), lambdamin)
       !write(*,*) "lambdacontmin = ", atom%continua(kr)%lambda(1), &
       !" lambdacontmax = ", atom%continua(kr)%lambda(atom%continua(kr)%Nlambda)
     else
      write(*,*) "Invalid value for continuum chromatic ",&
         "dependence."
      write(*,*) "exiting..."
      stop
     end if

     if (atom%active) then
      allocate(atom%continua(kr)%Rij(atmos%Nspace))
      allocate(atom%continua(kr)%Rji(atmos%Nspace))
     end if
    end do !end loop over bound-free transitions

    ! now fixed transitions
    ! fixed transitions are usefull to describe NLTE problem
    ! in the solar chromosphere in 2D,3D without not so much
    ! computational cost.
    ! They are not implemented because only relevent for solar
    ! application
    if (atom%Nfixed.gt.0) then
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
     CALL make_sub_wavelength_grid_line(atom%lines(kr), &
                                        minval(atom%vbroad,mask=atom%vbroad>0)) !MAXVAL(atom%vbroad)
  end do
   !Now even for passive atoms we write atomic data.
   ! Unlike RH, all data are in the same fits file.
    if (atom%ID(2:2) .eq." ") then
      atom%dataFile = atom%ID(1:1)//".fits.gz"
    else
      atom%dataFile = atom%ID(1:2)//".fits.gz"
    end if
    !write(*,*) "Populations file for writing: .", &
    !          trim(atom%popsoutFile),"."

   if (atom%active) then
!    !create file name for writing populations
!    !for each atom. File will be created at the end
!    !of the NLTE calculations
!     if (atom%ID(2:2) .eq." ") then
!       atom%popsoutFile = "pops."//atom%ID(1:1)//".fits.gz"
!     else
!       atom%popsoutFile = "pops."//atom%ID(1:2)//".fits.gz"
!     end if
!     !write(*,*) "Populations file for writing: .", &
!     !         trim(atom%popsoutFile),"."

    if (atom%Npfr.gt.0) then
     write(*,*) "PFR not implemented yet, do not write file"
    end if
    if (atmos%XRD) then
     write(*,*) "Cross redistribution ", &
                "not implemented yet."
    end if

    ! allocate some space
    atom%offset_coll = ftell(atomunit)
    atom%colunit = atom%periodic_table*2 + 1
    write(*,*) "offset file for reading collision", atom%offset_coll, "unit = ", atom%colunit
    !!allocate(atom%C(atom%Nlevel*atom%Nlevel,atmos%Nspace))
    !!now Collision matrix is constructed cell by cell, therefore allocated elsewhere
    allocate(atom%n(atom%Nlevel,atmos%Nspace))
    atom%n = 0d0
    atom%NLTEpops = .true.
   else !not active
    if (atom%dataFile.ne."" & !atom%popsinFile.ne.""
       .and. atom%initial_solution &
        .eq. "OLD_POPULATIONS") then
       atom%NLTEpops = .true.
       write(*,*) "First define a function that writes"
       write(*,*) "populations as fits and read it"
       stop
       allocate(atom%n(atom%Nlevel,atmos%Nspace))
       !CALL readPopulations(atom)
    else
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
  type (AtomType) :: atom

  !create formatline
  !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
  !write(FormatLine,'("("I3,A")")') MAX_LENGTH,"A"
  write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH


  atmos%Nactiveatoms = 0
  atmos%Npassiveatoms = 0
  open(unit=unit,file=trim(mcfost_utils)//TRIM(ATOMS_INPUT), status="old")

  !get number of atomic models to read
  CALL getnextline(unit, COMMENT_CHAR,FormatLine, &
       inputline, Nread)
  read(inputline,*) atmos%Natom
  write(*,*) "Reading ", atmos%Natom, " species"
  !Allocate sapace for Natom in atmos%Atoms
  allocate(atmos%Atoms(atmos%Natom))

  do nmet = 1, atmos%Natom
   CALL getnextline(unit, COMMENT_CHAR, FormatLine, &
     inputline, Nread)
   !write(*,*) "nmet=",nmet, inputline
   read(inputline,'(1A28, 1A7, 1A22, 1A20)') &
        filename, actionKey, popsKey, popsFile
   !write(*,*) ".",trim(filename),"."
   !write(*,*) ".",adjustl(actionKey),"."
   !write(*,*) ".",adjustl(popsKey),"."
   !write(*,*) ".",trim(popsFile),"."
   atmos%Atoms(nmet)%initial_solution=adjustl(popsKey)
   atmos%Atoms(nmet)%inputFile=trim(filename)


   ! would be pssoible in the future to read all J value also
   if (atmos%Atoms(nmet)%initial_solution.ne."OLD_POPULATIONS" &
      .and.atmos%Atoms(nmet)%initial_solution.ne."LTE_POPULATIONS"&
      .and.atmos%Atoms(nmet)%initial_solution.ne."ZERO_RADIATION"&
       .and.atmos%Atoms(nmet)%initial_solution.ne."SOBOLEV") then
     write(*,*) "Initial solution ", atmos%Atoms(nmet)%initial_solution,&
      " unkown!"
     write(*,*) "Exiting..."
     stop
   end if
!   atmos%Atoms(nmet)%popsinFile = trim(popsFile)
   atmos%Atoms(nmet)%dataFile = trim(popsFile) ! now the file atomID.fits.gz contains
                                               ! all informations including populations.


   !Active atoms are treated in NLTE
   if (adjustl(actionKey).eq."ACTIVE") then
     atmos%atoms(nmet)%active=.true.
     !write(*,*) "atom is active"
     atmos%Nactiveatoms = atmos%Nactiveatoms+1
   else
     atmos%atoms(nmet)%active=.false.
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
    if (atmos%atoms(mmet-1)%ID.eq.IDread) then
     write(*,*) "Already read a model for this atom ", IDread
     write(*,*) "exiting..."
     stop
    end if
   end do

   ! read and fill atomic structures saved in atmos%atoms
   ! create an alias for Hydrogen
   CALL readModelAtom(unit+nmet, atmos%Atoms(nmet), &
        trim(mcfost_utils)//TRIM(path_to_atoms)//trim(filename))
   !write(*,*) "IS ACTIVE = ", atmos%Atoms(nmet)%active
  end do
  close(unit)
  ! Alias to the most importent one
  Hydrogen=>atmos%Atoms(1)

  ! Aliases to active atoms
  if (atmos%Nactiveatoms.gt.0) then
     nact = 1
    allocate(atmos%ActiveAtoms(atmos%Nactiveatoms))
  end if
  if (atmos%Natom - atmos%Nactiveatoms == atmos%Npassiveatoms) then
   allocate(atmos%PassiveAtoms(atmos%Npassiveatoms))
   npass = 1
  else
   write(*,*) "Error, number of passive atoms is not N-Nactive"
   stop
  end if

  ! keep a duplicate in Elements
  write(*,*) "order#     ID   periodic-table#    ACTIVE    #lines   #continua"
  do nmet=1,atmos%Natom
   write(*,*) nmet, atmos%Atoms(nmet)%ID, &
    atmos%Atoms(nmet)%periodic_table, atmos%Atoms(nmet)%active, &
    atmos%Atoms(nmet)%Nline, atmos%Atoms(nmet)%Ncont
   ! create alias in atmos%Elements for elements that have
   ! a model atom.
   atmos%Elements(atmos%Atoms(nmet)%periodic_table)%model &
         => atmos%Atoms(nmet)
   if (atmos%atoms(nmet)%periodic_table.eq.2)  &
           Helium => atmos%Atoms(nmet)
   if (associated(atmos%ActiveAtoms)) then
     if (atmos%Atoms(nmet)%active) then 
      atmos%ActiveAtoms(nact:nact) => atmos%Atoms(nmet:nmet)
      nact = nact + 1 !got the next index of active atoms
     end if
   end if
   if (associated(atmos%PassiveAtoms)) then
     if (.not.atmos%Atoms(nmet)%active) then
      atmos%PassiveAtoms(npass:npass) => atmos%Atoms(nmet:nmet)
      npass = npass + 1
     end if
   end if
  end do


  RETURN
  END SUBROUTINE readAtomicModels


END MODULE readatom
