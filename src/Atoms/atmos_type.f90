MODULE atmos_type

  use atom_type, only : AtomType, Element, ATOM_ID_WIDTH, &
    AtomicLine, AtomicContinuum
  use uplow
  use constant
  
  !$ use omp_lib
  
  !MCFOST's originals
  use mcfost_env, only : mcfost_utils ! convert from the relative location of atomic data
                                      ! to mcfost's environnement folders.
  use parametres
  
  IMPLICIT NONE


  integer, parameter :: Nelem=99

  ! relative location in mcfost/utils/ specified by the mcfost environnement
  ! variables.
  character(len=50), parameter :: ABUNDANCE_FILE="/Atoms/abundance.input"
  character(len=50), parameter :: KURUCZ_PF_FILE="/Atoms/pf_Kurucz.fits.gz"

  TYPE atomPointerArray
   type(AtomType), pointer :: ptr_atom => NULL()
  END TYPE atomPointerArray

  ! Initialized only once at the begening of the main iteration loop !
  TYPE GridType
   ! Nspace is the number of cells, n_cells meaning that all quantities are computed
   ! for each cell starting from 1 to n_cells
   integer :: Nspace, Npf, Nactiveatoms, Npassiveatoms, Natom, Nrays = 1
   double precision :: metallicity, totalAbund, avgWeight, wght_per_H
   double precision :: tau !a variable used to sum up contribution of opacities along a ray
   ! an arrray containing the project local velocity
   ! in the cell, to be used to compute the local profile
   ! T is a function of x,y,z or T is a function of the cell point (T(Nspace))
   ! Each component of the velocity is a function of x, y, z so
   ! ux(x,y,z) = ux(Nspace) and vel=(ux,uy,uz)
   ! Don't know yet if it is useful here
   double precision, allocatable, dimension(:,:) :: Vxyz, Bxyz
   type (Element), dimension(:), allocatable :: Elements
   !!type (AtomType), pointer, dimension(:) :: Atoms, ActiveAtoms, PassiveAtoms
   !type (AtomType), dimension(:), allocatable :: Atoms, ActiveAtoms, PassiveAtoms
   type (atomPointerArray), dimension(:), allocatable :: Atoms, ActiveAtoms, PassiveAtoms
   !type (AtomType), pointer :: Hydrogen => NULL(), Helium => NULL()
   double precision, dimension(:), allocatable :: nHtot, ne, Tpf, T, vturb
   double precision, dimension(:), allocatable :: nHmin !Hminus populations
   double precision :: B_char = 0d0, v_char=0d0
           !B_char in Tesla and v_char in m/s, default 0T and 1km/s
   logical :: Magnetized = .false., XRD=.false., calc_ne, &
     H_LTE=.false. ! force LTE populations of H for
                   ! background opacities
   character(len=20) :: B_SOLUTION
   logical, allocatable, dimension(:) :: lcompute_atomRT !where line RT is taken into account on the grid
  END TYPE GridType


  TYPE (GridType), target :: atmos
  TYPE (AtomType), pointer :: Hydrogen => NULL(), Helium => NULL()


  CONTAINS

  SUBROUTINE freeAtom(atom)
  ! This subroutine free atmos%atoms properly including lines
  ! and continua.

   integer :: n, l, c
   type (AtomType), intent(inout) :: atom
   type (AtomicLine) :: line
   type (AtomicContinuum) :: cont
   
    !why error with deallocate(lines)
    atom%ID="" !just for checking

    if (allocated(atom%label)) deallocate(atom%label)
    if (allocated(atom%stage)) deallocate(atom%stage)
    if (allocated(atom%Lorbit)) deallocate(atom%Lorbit)
    if (allocated(atom%g)) deallocate(atom%g)
    if (allocated(atom%E)) deallocate(atom%E)
    if (allocated(atom%vbroad)) deallocate(atom%vbroad)
    if (allocated(atom%ntotal)) deallocate(atom%ntotal)
    if (allocated(atom%qS)) deallocate(atom%qS)
    if (allocated(atom%qJ)) deallocate(atom%qJ)
    if (allocated(atom%C)) deallocate(atom%C)
    if (allocated(atom%Gamma)) deallocate(atom%Gamma)

    !! maintenant se sont des pointeurs n et nstar
    !! so we need ot check if %n is an alias or not.
    !! if it is an alias, when %nstar is deallocated,
    !! %n needs to be nullified, not deallocated as
    !! it is already delinked from %nstar
    if ((atom%active).or.(atom%NLTEpops)) then
     if (associated(atom%n)) deallocate(atom%n)
     if (associated(atom%nstar)) deallocate(atom%nstar)
    else !passive 
      NULLIFY(atom%n)
      if (associated(atom%nstar)) deallocate(atom%nstar)
    end if
!     if (allocated(atom%gij)) deallocate(atom%gij)
!     if (allocated(atom%Vij)) deallocate(atom%Vij)
    if (allocated(atom%chi_up)) deallocate(atom%chi_up)
    if (allocated(atom%chi_down)) deallocate(atom%chi_down)
    if (allocated(atom%Uji_down)) deallocate(atom%Uji_down)
    if (allocated(atom%eta)) deallocate(atom%eta)

!!!     for this atom, free lines if allocated
     if (allocated(atom%lines)) then
     do l=1,atom%Nline
      line = atom%lines(l)
      if (allocated(line%phi)) deallocate(line%phi)
      if (allocated(line%phi_Q)) deallocate(line%phi_Q)
      if (allocated(line%phi_U)) deallocate(line%phi_U)
      if (allocated(line%phi_V)) deallocate(line%phi_V)
      if (allocated(line%psi_Q)) deallocate(line%psi_Q)
      if (allocated(line%psi_U)) deallocate(line%psi_U)
      if (allocated(line%psi_V)) deallocate(line%psi_V)
      if (allocated(line%lambda)) deallocate(line%lambda)
      !The net cooling rates is stored instead of the radiative rates
      !which are now scalar.
      if (allocated(line%CoolRates_ij)) deallocate(line%CoolRates_ij)
!       if (allocated(line%Rij)) deallocate(line%Rij)
!       if (allocated(line%Rji)) deallocate(line%Rji)
      if (allocated(line%wlam)) deallocate(line%wlam)
      !if (allocated(line%Qelast)) deallocate(line%Qelast)
      !if (allocated(line%c_shift)) deallocate(line%c_shift)
      !if (allocated(line%c_fraction)) deallocate(line%c_fraction)
      !if (allocated(line%adamp)) deallocate(line%adamp)
      if (allocated(line%rho_pfr)) deallocate(line%rho_pfr)
      if (allocated(line%gij)) deallocate(line%gij)
      if (allocated(line%Vij)) deallocate(line%Vij)
!       if (allocated(line%chi_up)) deallocate(line%chi_up)
!       if (allocated(line%chi_down)) deallocate(line%chi_down)
!       if (allocated(line%Uji_down)) deallocate(line%Uji_down)
     end do
      deallocate(atom%lines)
     end if
!!!     now free continua if allocated
     if (allocated(atom%continua)) then
     do c=1, atom%Ncont
      cont = atom%continua(c)
      if (allocated(cont%lambda)) deallocate(cont%lambda)
      if (allocated(cont%alpha)) deallocate(cont%alpha)
      if (allocated(cont%CoolRates_ij)) deallocate(cont%CoolRates_ij)
      if (allocated(cont%wlam)) deallocate(cont%wlam)
      if (allocated(cont%gij)) deallocate(cont%gij)
      if (allocated(cont%Vij)) deallocate(cont%Vij)
!       if (allocated(cont%chi_up)) deallocate(cont%chi_up)
!       if (allocated(cont%chi_down)) deallocate(cont%chi_down)
!       if (allocated(cont%Uji_down)) deallocate(cont%Uji_down)
     end do
      deallocate(atom%continua)
     end if

  RETURN
  END SUBROUTINE freeAtom

  SUBROUTINE readKurucz_pf(code, Npf, Nstage, ionpot, pf)
   !return pf function for atom with Z = code
   ! there is an HDU for each element (+one for T grid already read)
   integer :: unit, EOF, blocksize, i, j
   integer, intent(in) :: code
   integer :: NAXISPF(2), naxis_found, hdutype, Z, shape_pf(2)
   character(len=256) :: some_comments
   integer,  intent(inout) :: Nstage
   integer,  intent(in) :: Npf
   real(8), allocatable, dimension(:,:) :: data_krz
   real(8), allocatable, intent(inout), dimension(:,:) :: pf
   real(8), allocatable, intent(inout), dimension(:) :: ionpot
   logical :: anynull

   EOF = 0
   !get unique unit number
   CALL ftgiou(unit,EOF)
   ! open fits file in readmode'
   CALL ftopen(unit, TRIM(mcfost_utils)//TRIM(KURUCZ_PF_FILE), 0, blocksize, EOF)

    !move to nexthdu, 1=Primary, already read
    !Hydrogen = 2 = code+1
    CALL FTMAHD(unit,code+1,hdutype,EOF)
    if (EOF.ne.0) then
     write(*,*) "error reading partition function"
     write(*,*) "error to change HDU at ", code+1
     write(*,*) "exiting..."
     stop
    end if
    CALL ftgkyj(unit, "NSTAGE", Nstage, some_comments, EOF) !replace j with d if read double instead of int
    CALL ftgkyj(unit, "Z", Z, some_comments, EOF)
    !j'arrive pas a lire ce string, but it is not necessary
    !but should asks christophe to know how to do
    !CALL ftgkyj(unit, "ELEM", AtomID, commentfits, EOF)

    !write(*,*) "Nstage", Nstage, "Z=", Z, "for elem:", elemental_ID(code)
    !write(*,*) "Npf points = ", Npf
    !now read partition function
    CALL ftgknj(unit, 'NAXIS', 1, 2, NAXISPF, naxis_found, EOF)
    if (NAXISPF(1).ne.Npf.and.NAXISPF(2).ne.Nstage) then
      write(*,*) "error reading partition function"
      write(*,*) "NAXISPF(len=2)=",NAXISPF," NPOINTS=", &
                 Npf, " NSTAGE=", Nstage
      write(*,*) "exiting"
      stop
    end if
    !data are actually transposed from my fits to fortran
    allocate(data_krz(Npf+1,Nstage))
    allocate(ionpot(Nstage))
    allocate(pf(Nstage, Npf))
    !now read the value of pf for that atom
    ! remember: pf(:,1) =  ion potentials
    !           pf(:,2:) = partition functions for all ion pots.
    !e = real*4
    !d = real*8            !shape(data_krz)
    shape_pf(1) = Npf+1; shape_pf(2) = Nstage
    CALL FTG2Dd(unit,1,-999,shape_pf,Npf+1,Nstage,data_krz,anynull,EOF)
    !do i=1,Nstage
    ! write(*,*) "potential in cm-1 for elem ", elemental_ID(code),":", &
    !            data_krz(1,i)
     !fill ionpot array and pf array for that elem
     ! store value in Joule
     ionpot(:) = data_krz(1,:) * HPLANCK*CLIGHT / (CM_TO_M)
     !write(*,*) "chi(i) = ", data_krz(1,i)
     !do not forget to write partition function has (Nstage,Npf)
     !instead of (Npf, Nstage) as read by FTG2D
     !store value in logarithm of partition function, for interpolation
     do i=1,Nstage
      do j=2,Npf+1
       !!if (code.eq.1) write(*,*) 'pf=', data_krz(j,i)
       pf(i,j-1) = LOG10(data_krz(j,i))
       !!if (code.eq.1) write(*,*) 'pf10powLog=', 10**(pf(i,j-1))
      end do
     end do
    !end do
    !check shapes
    !write(*,*) "Shapes ionpot and pf: ", shape(ionpot), shape(pf)


    !close file
    call ftclos(unit, EOF)
    ! free the unit number
    call ftfiou(unit, EOF)
    deallocate(data_krz)
  RETURN
  END SUBROUTINE readKurucz_pf

  SUBROUTINE free_atomic_atmos()
  integer :: n

  if (allocated(atmos%Elements)) then
    do n=1,Nelem
     if (allocated(atmos%Elements(n)%mol_index)) deallocate (atmos%Elements(n)%mol_index)
     if (allocated(atmos%Elements(n)%ionpot)) deallocate (atmos%Elements(n)%ionpot)
     if (allocated(atmos%Elements(n)%pf)) deallocate (atmos%Elements(n)%pf)
     if (allocated(atmos%Elements(n)%n))  deallocate (atmos%Elements(n)%n)
     !if alias exist for this element free it.
     if (associated(atmos%Elements(n)%model)) NULLIFY(atmos%Elements(n)%model)
!       ! de-alias %model to the proper atmos%Atoms(nmet), atmos%Atoms will be freed after.
    end do
   deallocate(atmos%Elements) !here it works
  end if
  
  if (allocated(atmos%nHtot)) deallocate(atmos%nHtot)
  if (allocated(atmos%ne)) deallocate(atmos%ne)
  if (allocated(atmos%Tpf)) deallocate(atmos%Tpf)
  if (allocated(atmos%T)) deallocate(atmos%T)
  if (allocated(atmos%vturb)) deallocate(atmos%vturb)
  if (allocated(atmos%nHmin)) deallocate(atmos%nHmin)
  if (allocated(atmos%lcompute_atomRT)) deallocate(atmos%lcompute_atomRT)
  if (allocated(atmos%Vxyz)) deallocate(atmos%Vxyz)
  if (allocated(atmos%Bxyz)) deallocate(atmos%Bxyz)

  !write(*,*) "Free Atoms"
  ! start freeing Atoms if previously allocated
  ! in readAtom()
   !if (allocated(atmos%Atoms)) then
  do n=1,atmos%Npassiveatoms
   if (associated(atmos%PassiveAtoms(n)%ptr_atom)) NULLIFY(atmos%PassiveAtoms(n)%ptr_atom)
  end do
  if (allocated(atmos%PassiveAtoms)) deallocate(atmos%PassiveAtoms)
  do n=1,atmos%Nactiveatoms
   if (associated(atmos%ActiveAtoms(n)%ptr_atom)) NULLIFY(atmos%ActiveAtoms(n)%ptr_atom)
  end do
  if (allocated(atmos%ActiveAtoms)) deallocate(atmos%ActiveAtoms)
  do n=1,atmos%Natom
   if (associated(atmos%Atoms(n)%ptr_atom)) then !now c'est un pointeur array
!    ! H is alway presents as it is required !
!    ! but He exists only if a model atom is given and so
!    ! Helium is allocated only if He exists i.e., if
!    ! atmos%Atoms(n)%periodic_tabel.eq.2 with n a number.
    if (associated(Hydrogen)) NULLIFY(Hydrogen) ! Hydrogen is always associated
    if (associated(Helium)) NULLIFY(Helium)
!    do n=1,atmos%Natom
     CALL freeAtom(atmos%Atoms(n)%ptr_atom)
!    end do
!    if (associated(atmos%ActiveAtoms)) deallocate(atmos%ActiveAtoms)
!   deallocate(atmos%Atoms)
    NULLIFY(atmos%Atoms(n)%ptr_atom)
   end if
  end do
  deallocate(atmos%Atoms)
  RETURN
  END SUBROUTINE free_atomic_atmos

  SUBROUTINE fillElements()
  !This routine read the abundance of elements listed in the
  ! abundance.input file, and keep the log partition function values
  ! for each elements in atmos%Elements
  !type (GridType) :: grid
  type (Element) :: elem
  integer :: EOF, n, blocksize, unit, i, j, k, ni
  integer :: NAXIST, naxis_found, hdutype, Z
  character(len=256) :: some_comments
  logical :: anynull
  character(len=4) :: charID
  real(8) :: A

  EOF = 0
  !Start reading partition function by first reading the grid
  !get unique unit number
  CALL ftgiou(unit,EOF)
  ! open fits file in readmode'
  CALL ftopen(unit, TRIM(mcfost_utils)//TRIM(KURUCZ_PF_FILE), 0, blocksize, EOF)
  !some check about the first axis!
  !get first axis
  CALL ftgknj(unit, 'NAXIS', 1, 1, NAXIST, naxis_found, EOF)
  !read dimension of partition function table
  CALL ftgkyj(unit, "NPOINTS", atmos%Npf, some_comments, EOF)
  if (NAXIST.ne.atmos%Npf) then
   write(*,*) "Error in reading pf data !"
   write(*,*) "NAXIST=",NAXIST," NPOINTS=", atmos%Npf
   write(*,*) "exiting... !"
   stop
  end if
  !write(*,*) NAXIST, atmos%Npf

  allocate(atmos%Tpf(atmos%Npf))
  !now read temperature grid
  CALL ftgpvd(unit,1,1,atmos%Npf,-999,atmos%Tpf,anynull,EOF) !d because double !!
!   do k=1,atmos%Npf
!    write(*,*) "T=", atmos%Tpf(k), "K"
!   end do
  !close file
  call ftclos(unit, EOF)
  ! free the unit number
  call ftfiou(unit, EOF)

  n = 1
  EOF = 0
  !write(*,*) "Metallicity = ", atmos%metallicity

  !read abundances
  open(unit=1, file=TRIM(mcfost_utils)//TRIM(ABUNDANCE_FILE),status="old")
  do while (EOF.eq.0)
   read(1, '(1A4, 1F5.3)', IOSTAT=EOF) charID, A
   if (charID(1:1) /= '#') then
    atmos%Elements(n)%weight=atomic_weights(n)
    atmos%Elements(n)%ID = charID(1:2) !supposed to be lowercases!!!!!
    atmos%Elements(n)%abund = 10.**(A-12.)
    if (A <= -99) atmos%Elements(n)%abund = 0d0
    !write(*,*) " abund=",atmos%Elements(n)%abund, A
    atmos%Elements(n)%abundance_set=.true.

    !write(*,*) "Element ", atmos%Elements(n)%ID," has an abundance of A=",&
    !           atmos%Elements(n)%abund," and a weight of w=",&
    !           atmos%Elements(n)%weight
    if (n.gt.1 .and. atmos%metallicity.gt.0.) then
     atmos%Elements(n)%abund = atmos%Elements(n)%abund * &
                             10.**(atmos%metallicity)
    end if
    atmos%totalAbund = atmos%totalAbund+atmos%Elements(n)%abund
    atmos%avgWeight = atmos%avgWeight+atmos%Elements(n)%abund* &
                     atmos%Elements(n)%weight
    !write(*,*) "updating total Abundance = ",atmos%totalAbund
    !write(*,*) "updating average weight = ", atmos%avgWeight

    !attribute partition function to each element
    !the temperature grid is shared and saved in atmos%Tpf
    !write(*,*) "Fill partition function for element ", atmos%Elements(n)%ID
    CALL readKurucz_pf(n, atmos%Npf, &
              atmos%Elements(n)%Nstage, atmos%Elements(n)%ionpot, &
              atmos%Elements(n)%pf)
    !do i=1,atmos%Elements(n)%Nstage !convert from J->eV
    ! write(*,*) "Inpot for stage ", i, &
    !            atmos%Elements(n)%ionpot(i)*JOULE_TO_EV, " eV"
    !end do
    !! allocate space for futur populations of each stage
    !!allocate(atmos%Elements(n)%n(atmos%Elements(n)%Nstage, atmos%Nspace))
    !write(*,*) "Elem pop" !BEWARE of n index and Elements(n)%n
    ! there is a conflict
!     do ni=1,Nelem
!      atmos%Elements(n)%n(:,:)=0d0
!     end do

    !write(*,*) "************************************"
    n = n + 1 !go to next element
   end if
   if (n.gt.Nelem) then
    exit
   end if
  end do
  atmos%wght_per_H = atmos%avgWeight
  atmos%avgWeight = atmos%avgWeight/atmos%totalAbund
  !write(*,*) "Final total A = ", atmos%totalAbund
  !write(*,*) "Final total average weight = ", atmos%avgWeight
  !write(*,*) "Weight per Hydrogen = ", atmos%wght_per_H
  close(unit=1)

  RETURN
  END SUBROUTINE fillElements

  SUBROUTINE init_atomic_atmos()!(Nspace)
   !integer, intent(in) :: Nspace
   atmos%Nspace = n_cells
   
   if (allocated(atmos%T).or.(allocated(atmos%nHtot)) & 
       .or.(allocated(atmos%ne))) then
    write(*,*) "A atomic atmosphere is already allocated, exiting..."
    stop
   end if
   
   if (.not.allocated(atmos%T)) allocate(atmos%T(atmos%Nspace))
   if (.not.allocated(atmos%Elements)) allocate(atmos%Elements(Nelem))
   
   if (.not.lstatic.and..not.lvoronoi) then 
    !Nrad, Nz, Naz-> Velocity vector cartesian components
    allocate(atmos%Vxyz(atmos%Nspace,3))
    atmos%Vxyz = 0d0
    !if lmagnetoaccr then atmos%Vxyz = (VRcyl, Vz, Vphi)
    !it is then projected like it is done in v_proj for keplerian rot
    !meaning that the velocity field has infinite resolution. 
    !Otherwise,
    !Vxyz=(Vx,Vy,Vz) and the velocity is constant in the cell just like in the
    !Vornoi case. Note that this second case of (Vx,Vy,Vz) is just a test case
    !For generating a model for Voronoi Tesselation.
    !The normal case is atmos%Vxyz beeing (VRcyl, Vz, Vphi) for
    !lmagnetoaccr case only, where both cylindrical and spherical velocities vectors
    !overlap (different than lkeplerian with only Vphi and linfall with only Vrsph.
    !If Lvoronoi, it is not allated since velocity is in Voronoi(:)%vxyz
   end if

   atmos%Natom = 0
   atmos%Nactiveatoms = 0
   atmos%totalAbund=0
   atmos%avgWeight=0
   atmos%wght_per_H=0
   atmos%metallicity = 0. !not used yet
   atmos%Npf = 0
   atmos%calc_ne = .false.
   atmos%T = 0d0


   if (.not.allocated(atmos%nHtot)) allocate(atmos%nHtot(atmos%Nspace))
   if (.not.allocated(atmos%vturb)) allocate(atmos%vturb(atmos%Nspace))
   if (.not.allocated(atmos%ne)) allocate(atmos%ne(atmos%Nspace))
   if (.not.allocated(atmos%nHmin)) allocate(atmos%nHmin(atmos%Nspace))

   
   atmos%ne = 0d0
   atmos%nHmin = 0d0
   atmos%nHtot = 0d0
   atmos%vturb = 0d0 !m/s

   CALL fillElements()
   
   allocate(atmos%lcompute_atomRT(atmos%Nspace))

   RETURN
   END SUBROUTINE init_atomic_atmos

   SUBROUTINE define_atomRT_domain(itiny_T, itiny_nH)
   ! Set where to solve for the RT equation: where a cell is not 
   ! transparent = where there is a significant density and temperature.
   ! Determines also if we have to force electron density calculation
    integer :: icell
    double precision, optional :: itiny_nH, itiny_T
    double precision :: Vchar=0d0, tiny_nH=1d0, tiny_T=5d2
    !with Very low value of densities or temperature it is possible to have
    !nan or infinity in the populations calculations because of numerical precisions
    
    if (present(itiny_T)) then
     write(*,*) "changing the value of tiny_T (K) = ", tiny_T," to", itiny_T
     tiny_T = itiny_T !K
    end if

    if (present(itiny_nH)) then
     write(*,*) "changing the value of tiny_nH (m^-3) = ", tiny_nH," to", itiny_nH
     tiny_nH = itiny_nH !m^-3
    end if
    
    if (maxval(atmos%ne) == 0d0) atmos%calc_ne = .true.

    if (tiny_T <= 0) then
     write(*,*) "changing the value of tiny_T = ", tiny_T," to", 10d0
     tiny_T = 10d0
    end if
    if (tiny_nH <= 0) then
     write(*,*) "changing the value of tiny_nH = ", tiny_nH," to", 1d0
     tiny_nH = 1d0
    end if
    
    !atmos%lcompute_atomRT = (atmos%nHtot > tiny_nH) .and. (atmos%T > tiny_T) 
    do icell=1,atmos%Nspace
     atmos%lcompute_atomRT(icell) = &
       (atmos%nHtot(icell) > tiny_nH) .and. (atmos%T(icell) > tiny_T)
    end do

   RETURN
   END SUBROUTINE define_atomRT_domain
   
   SUBROUTINE init_magnetic_field()
   ! ----------------------------------------- !
   ! Allocate space for magnetic field
   ! and solution for the Zeeman polarisation.
   ! ----------------------------------------- !
    if (.not.atmos%magnetized) RETURN
    allocate(atmos%Bxyz(atmos%Nspace, 3))
    atmos%Bxyz = 0.0
    atmos%B_SOLUTION = "WEAK_FIELD"
    !in the weak field regime
    ! V = -delaLambda_B * geff * cos(incl_mag, normal) * dI/dlambda
   RETURN
   END SUBROUTINE init_magnetic_field

   FUNCTION atomZnumber(atomID) result(Z)
   ! --------------------------------------
   ! return the atomic number of an atom
   ! with ID = atomID.
   ! Hydrogen is 1
   ! --------------------------------------
    character(len=ATOM_ID_WIDTH) :: atomID
    integer :: Z, i

    Z = 1
    !do i=1,Nelem
    !  if (atmos%Elements(i)%ID.eq.atomID) then
    !    Z = i
    !    exit
    !  end if
    do while (atmos%Elements(Z)%ID.ne.atomID)
     Z=Z+1
    end do

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

END MODULE atmos_type

