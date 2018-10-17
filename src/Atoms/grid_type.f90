MODULE grid_type

  use atom_type, only : AtomType, Element, ATOM_ID_WIDTH, &
    AtomicLine, AtomicContinuum
  use uplow
  use constant

  IMPLICIT NONE


  integer, parameter :: Nelem=99
  character(len=50), parameter :: ABUNDANCE_FILE="../../utils/Atoms/abundance.input"
  character(len=50), parameter :: KURUCZ_PF_FILE="../../utils/Atoms/pf_Kurucz.fits.gz"


  ! Initialized only once at the begening of the main iteration loop !
  TYPE GridType
   integer :: Nspace, Ncells, Npf, Nactiveatoms, Natom
   real(8) :: metallicity, totalAbund, avgWeight, wght_per_H, wavelength_ref
   ! an arrray containing the project local velocity
   ! in the cell, to be used to compute the local profile
   ! T is a function of x,y,z or Nspace=Nx*Ny*Nz
   ! Each component of the velocity is a function of x, y, z so
   ! ux(x,y,z) = ux(Nspace) and vel=(ux,uy,uz)
   double precision, allocatable, dimension(:) :: ux, uy, uz, Bx, By, Bz !from the model
   type (Element), dimension(:), allocatable :: Elements
   type (AtomType), pointer, dimension(:) :: Atoms, ActiveAtoms 
   type (AtomType), pointer :: Hydrogen => NULL(), Helium => NULL()
   real(8), dimension(:), allocatable :: nHtot, ne, Tpf, T, vturb
   real(8), dimension(:), allocatable :: nHmin !Hminus populations
   real :: B_char = 0., vmicro_char=1.
   logical :: Magnetized = .false., XRD=.false., moving=.true., calc_ne,&
     H_LTE=.false. ! force LTE populations of H for
                   ! background opacities
  END TYPE GridType
  
  TYPE CellPoint
   double precision :: u, v, w
  END TYPE CellPoint

  TYPE (GridType), target :: atmos
  !alias to atmos%Atoms(1)
  TYPE (AtomType), pointer :: Hydrogen, Helium
  
  TYPE (CellPoint), pointer :: cell_point 

  CONTAINS

  SUBROUTINE freeAtom(atom)
  ! This subroutine free atmos%atoms properly including lines
  ! and continua.

   integer :: n, l, c
   type (AtomType), intent(inout) :: atom
   type (AtomicLine) :: line
   type (AtomicContinuum) :: cont
   
    !why error with deallocate(lines)


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
     deallocate(atom%n)
     deallocate(atom%nstar)
    else !passive 
      NULLIFY(atom%n)
      deallocate(atom%nstar)
    end if
    if (allocated(atom%gij)) deallocate(atom%gij)
    if (allocated(atom%Vij)) deallocate(atom%Vij)
    if (allocated(atom%wla)) deallocate(atom%wla)
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
      if (allocated(line%Rij)) deallocate(line%Rij)
      if (allocated(line%Rji)) deallocate(line%Rji)
      if (allocated(line%wphi)) deallocate(line%wphi)
      !if (allocated(line%Qelast)) deallocate(line%Qelast)
      if (allocated(line%c_shift)) deallocate(line%c_shift)
      if (allocated(line%c_fraction)) deallocate(line%c_fraction)
      !if (allocated(line%adamp)) deallocate(line%adamp)
      if (allocated(line%rho_prd)) deallocate(line%rho_prd)
     end do
!      deallocate(atom%lines)
     end if
!!!     now free continua if allocated
     if (allocated(atom%continua)) then
     do c=1, atom%Ncont
      cont = atom%continua(c)
      if (allocated(cont%lambda)) deallocate(cont%lambda)
      if (allocated(cont%alpha)) deallocate(cont%alpha)
      if (allocated(cont%Rij)) deallocate(cont%Rij)
      if (allocated(cont%Rji)) deallocate(cont%Rji)
     end do
      deallocate(atom%continua)
     end if

  RETURN
  END SUBROUTINE freeAtom

  SUBROUTINE readKurucz_pf(code, Npf, Nstage, ionpot, pf)
   !return pf function for atom with Z = code
   integer :: unit, EOF, blocksize, i, j
   integer, intent(in) :: code
   integer :: NAXISPF(2), naxis_found, hdutype, Z
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
   CALL ftopen(unit, TRIM(KURUCZ_PF_FILE), 0, blocksize, EOF)

    !move to nexthdu, 1=Primary, already read
    !Hydrogen = 2 = code+1
    CALL FTMAHD(unit,code+1,hdutype,EOF)
    if (EOF.ne.0) then
     write(*,*) "error reading partition function"
     write(*,*) "error to change HDU at ", code+1
     write(*,*) "exiting..."
     stop
    end if
    CALL ftgkyj(unit, "NSTAGE", Nstage, some_comments, EOF)
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
    !d = real*8
    CALL FTG2Dd(unit,1,-999,shape(data_krz),Npf+1,Nstage,data_krz,anynull,EOF)
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

  SUBROUTINE freeGrid()
  integer :: n

  !Some error of deallocate, how it works with derived types with allocatable,
  ! derived types as attributes and pointers and pointers in derived types attributes...
  ! should not deallocate(atmos%Atoms) working if at least all pointers have been deallocated
  ! first ?
  write(*,*) "Warning: check deallocation of pointers and derived types in derived types.."

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
  if (allocated(atmos%ux)) deallocate(atmos%ux)
  if (allocated(atmos%uy)) deallocate(atmos%uy)
  if (allocated(atmos%uz)) deallocate(atmos%uz)
  if (allocated(atmos%Bx)) deallocate(atmos%Bx)
  if (allocated(atmos%By)) deallocate(atmos%By)
  if (allocated(atmos%Bz)) deallocate(atmos%Bz)

  !write(*,*) "Free Atoms"
  ! start freeing Atoms if previously allocated
  ! in readAtom()
   if (associated(atmos%Atoms)) then !now c'est un pointeur array
!    ! H is alway presents as it is required !
!    ! but He exists only if a model atom is given and so
!    ! Helium is allocated only if He exists i.e., if
!    ! atmos%Atoms(n)%periodic_tabel.eq.2 with n a number.
    NULLIFY(Hydrogen) ! Hydrogen is always associated
    if (associated(Helium)) NULLIFY(Helium)
    do n=1,atmos%Natom
     CALL freeAtom(atmos%Atoms(n))
    end do
!    if (associated(atmos%ActiveAtoms)) deallocate(atmos%ActiveAtoms)
!    deallocate(atmos%Atoms)
   end if
  !write(*,*) "After free atoms"
  RETURN
  END SUBROUTINE freeGrid

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
  CALL ftopen(unit, TRIM(KURUCZ_PF_FILE), 0, blocksize, EOF)
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
  CALL ftgpve(unit,1,1,atmos%Npf,-999,atmos%Tpf,anynull,EOF)
  !do k=1,atmos%Npf
  ! write(*,*) "T=", atmos%Tpf(k), "K"
  !end do
  !close file
  call ftclos(unit, EOF)
  ! free the unit number
  call ftfiou(unit, EOF)

  n = 1
  EOF = 0
  !write(*,*) "Metallicity = ", atmos%metallicity

  !read abundances
  open(unit=1, file=ABUNDANCE_FILE,status="old")
  do while (EOF.eq.0)
   read(1, '(1A4, 1F5.3)', IOSTAT=EOF) charID, A
   if (charID(1:1) /= '#') then
    atmos%Elements(n)%weight=atomic_weights(n)
    atmos%Elements(n)%ID = charID(1:2) !supposed to be lowercases!!!!!
    atmos%Elements(n)%abund = 10.**(A-12.)
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
    ! allocate space for futur populations of each stage
    allocate(atmos%Elements(n)%n(atmos%Elements(n)%Nstage, atmos%Nspace))
    !write(*,*) "Elem pop" !BEWARE of n index and Elements(n)%n
    ! there is a conflict
    do ni=1,Nelem
     atmos%Elements(n)%n(:,:)=0.
    end do

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


  SUBROUTINE initGrid(Nspace, T, ne, nHtot, vturb, vel_char)
   integer, intent(in) :: Nspace
   integer :: k, n
   real(8), intent(in), dimension(Nspace) :: T, ne, nHtot, vturb
   real(8), intent(in) :: vel_char !km/s
   real(8) :: maxNe = 0

   if (.not.allocated(atmos%T)) allocate(atmos%T(Nspace))
   if (.not.allocated(atmos%Elements)) allocate(atmos%Elements(Nelem))
   if (.not.allocated(atmos%ux)) allocate(atmos%ux(Nspace))
   if (.not.allocated(atmos%uy)) allocate(atmos%uy(Nspace))
   if (.not.allocated(atmos%uz)) allocate(atmos%uz(Nspace))


   atmos%Natom = 0
   atmos%Nactiveatoms = 0
   atmos%totalAbund=0
   atmos%avgWeight=0
   atmos%wght_per_H=0
   atmos%metallicity = 0. !not used yet
   atmos%Npf = 0
   atmos%Nspace = Nspace
   atmos%wavelength_ref=500. !nm optionnal
   atmos%calc_ne = .false.

   if (.not.allocated(atmos%nHtot)) allocate(atmos%nHtot(Nspace))
   if (.not.allocated(atmos%vturb)) allocate(atmos%vturb(Nspace))
   if (.not.allocated(atmos%ne)) allocate(atmos%ne(Nspace))
   if (.not.allocated(atmos%nHmin)) allocate(atmos%nHmin(Nspace))

   write(*,*) "Use the proper vmicro_char in case of large ", &
     "velocity gradient" ! it is use to sample the line domain
   atmos%vmicro_char = vel_char*1000. !km/s to m/s
   ! Note that if velocity field are important, this term
   ! should be use as a typical width


   atmos%T = T !T is a flattenned array in N Dimension
   atmos%nHtot = nHtot*1d6 !m-3
   atmos%ne = ne*1d6 !cm-3->m-3
   atmos%vturb = vturb * 1000.

   if (MAXVAL(atmos%ne).eq.0.0) then 
     write(*,*) "Solving for electron density"
     atmos%calc_ne=.true.
   end if 
    
   CALL fillElements()

   RETURN
   END SUBROUTINE initGrid

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

END MODULE grid_type

