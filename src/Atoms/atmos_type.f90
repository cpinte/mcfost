MODULE atmos_type

  use atom_type, only : AtomType, Element, ATOM_ID_WIDTH, &
    AtomicLine, AtomicContinuum, AtomicTransition
  use uplow
  use constant
  use messages
  use constantes, only : tiny_dp

  
  !$ use omp_lib
  
  !MCFOST's originals
  use mcfost_env, only : mcfost_utils ! convert from the relative location of atomic data
                                      ! to mcfost's environnement folders.
  use parametres
  use utils, only : appel_syst
  
  IMPLICIT NONE


  integer, parameter :: Nelem=99

  ! relative location in mcfost/utils/ specified by the mcfost environnement
  ! variables.
  character(len=50), parameter :: ABUNDANCE_FILE="./abundance.input"!"/Atoms/abundance.input"
  character(len=50), parameter :: KURUCZ_PF_FILE="/Atoms/pf_Kurucz.fits.gz"
  !path to your partition function data

  !Wrappers to create arrays of Pointers
!!-> move then to Atom_type.f90
  TYPE atomPointerArray
   type(AtomType), pointer :: ptr_atom => NULL()
  END TYPE atomPointerArray
  
  !Is this one really needed ? and why it doesn't work in solvene ?
  TYPE elemPointerArray
   type(Element), pointer :: ptr_elem => NULL()
  END TYPE elemPointerArray

  ! Initialized only once at the begening of the main iteration loop !
  TYPE GridType
   ! Nspace is the number of cells, n_cells meaning that all quantities are computed
   ! for each cell starting from 1 to n_cells
   integer :: Nspace, Npf, Nactiveatoms, Npassiveatoms, Natom, Nrays = 1
   double precision :: metallicity, totalAbund, avgWeight, wght_per_H
   
   !! will be deprecated as I do not need optical_length_tot for stellar map. Only for checking
   double precision :: tau !a variable used to sum up contribution of opacities along a ray
   !!
   ! an arrray containing the project local velocity
   ! in the cell, to be used to compute the local profile
   ! T is a function of x,y,z or T is a function of the cell point (T(Nspace))
   ! Each component of the velocity is a function of x, y, z so
   ! ux(x,y,z) = ux(Nspace) and vel=(ux,uy,uz)
   ! Don't know yet if it is useful here
   double precision, allocatable, dimension(:,:) :: Vxyz, Bxyz
   !type (Element), dimension(:), allocatable :: Elements
   type (elemPointerArray), dimension(:), allocatable :: Elements
   !!type (AtomType), pointer, dimension(:) :: Atoms, ActiveAtoms, PassiveAtoms
   !type (AtomType), dimension(:), allocatable :: Atoms, ActiveAtoms, PassiveAtoms
   type (atomPointerArray), dimension(:), allocatable :: Atoms, ActiveAtoms, PassiveAtoms
   !type (AtomType), pointer :: Hydrogen => NULL(), Helium => NULL()
   double precision, dimension(:), allocatable :: nHtot, ne, Tpf, T!, vturb
   double precision, dimension(:), allocatable :: nHmin !Hminus populations
   double precision :: B_char = 0d0, v_char=0d0
   logical :: include_xcoupling = .false., coherent_scattering =.false.
           !B_char in Tesla and v_char in m/s, default 0T and 1km/s
   logical :: Magnetized = .false., XRD=.false., calc_ne
   !dark zone are in lcompute_atomRT->icompute_atomRT == -1
   integer, allocatable, dimension(:) :: icompute_atomRT!, ldark_zone!where line RT is taken into account on the grid
   !logical, allocatable, dimension(:) :: lcompute_atomRT!, ldark_zone!where line RT is taken into account on the grid
  END TYPE GridType


  TYPE (GridType), target :: atmos
  TYPE (AtomType), pointer :: Hydrogen => NULL(), Helium => NULL()


  CONTAINS
 
 FUNCTION VBROAD_atom(icell, atom) result(vD)
  double precision :: vD
  integer :: icell
  type (AtomType) :: atom
  
   vD = dsqrt(Vtherm*atmos%T(icell)/atom%weight)! + atmos%vturb(icell)**2)
   
 RETURN
 END FUNCTION VBROAD_atom
 
 FUNCTION nTotal_atom(icell, atom) result(ntotal)
  double precision :: ntotal
  integer :: icell
  type (AtomType) :: atom
  
   ntotal =  atom%Abund * atmos%nHtot(icell)
  
 RETURN
 END FUNCTION nTotal_atom
  
  SUBROUTINE freeZeemanMultiplet(line)
   type(AtomicLine), intent(inout) :: line

   if (allocated(line%zm%strength)) deallocate(line%zm%strength)
   if (allocated(line%zm%q)) deallocate(line%zm%q)
   if (allocated(line%zm%shift)) deallocate(line%zm%shift)
   
  RETURN
  END SUBROUTINE freeZeemanMultiplet
 
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
    !if (allocated(atom%vbroad)) deallocate(atom%vbroad)
    !if (allocated(atom%ntotal)) deallocate(atom%ntotal)
    if (allocated(atom%qS)) deallocate(atom%qS)
    if (allocated(atom%qJ)) deallocate(atom%qJ)
    !if (allocated(atom%C)) deallocate(atom%C)
    !free collision here
    if (allocated(atom%Gamma)) deallocate(atom%Gamma)
    if (allocated(atom%collision_lines)) deallocate(atom%collision_lines)

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
    if (allocated(atom%Xc)) deallocate(atom%Xc)

    !if (allocated(atom%chi_up)) deallocate(atom%chi_up)
    !if (allocated(atom%chi_down)) deallocate(atom%chi_down)
    !if (allocated(atom%Uji_down)) deallocate(atom%Uji_down)
    if (allocated(atom%eta)) deallocate(atom%eta)

!!!     for this atom, free lines if allocated
     if (allocated(atom%lines)) then
     do l=1,atom%Nline
      line = atom%lines(l)
      if (allocated(line%phi)) deallocate(line%phi)
!       if (allocated(line%phi_Q)) deallocate(line%phi_Q)
!       if (allocated(line%phi_U)) deallocate(line%phi_U)
!       if (allocated(line%phi_V)) deallocate(line%phi_V)
!       if (allocated(line%psi_Q)) deallocate(line%psi_Q)
!       if (allocated(line%psi_U)) deallocate(line%psi_U)
!       if (allocated(line%psi_V)) deallocate(line%psi_V)
      if (allocated(line%lambda)) deallocate(line%lambda)
      !The net cooling rates is stored instead of the radiative rates
      !which are now scalar.
      if (allocated(line%CoolRates_ij)) deallocate(line%CoolRates_ij)
!       if (allocated(line%Rij)) deallocate(line%Rij)
!       if (allocated(line%Rji)) deallocate(line%Rji)
      !if (allocated(line%wlam)) deallocate(line%wlam) !over frequencies and angles and proc
      !if (allocated(line%Qelast)) deallocate(line%Qelast)
      !if (allocated(line%c_shift)) deallocate(line%c_shift)
      !if (allocated(line%c_fraction)) deallocate(line%c_fraction)
      !if (allocated(line%adamp)) deallocate(line%adamp)
      if (allocated(line%rho_pfr)) deallocate(line%rho_pfr)
      !if (allocated(line%gij)) deallocate(line%gij)
      !if (allocated(line%Vij)) deallocate(line%Vij)
!       if (allocated(line%chi_up)) deallocate(line%chi_up)
       !if (allocated(line%chi)) deallocate(line%chi)
       !if (allocated(line%U)) deallocate(line%U)
      !!if (allocated(line%Jbar)) deallocate(line%Jbar)
      if (associated(line%atom)) NULLIFY(line%atom)
      CALL freeZeemanMultiplet(line)
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
      !if (allocated(cont%wlam)) deallocate(cont%wlam)
      !if (allocated(cont%gij)) deallocate(cont%gij)
      !if (allocated(cont%Vij)) deallocate(cont%Vij)
!       if (allocated(cont%chi_up)) deallocate(cont%chi_up)
       !if (allocated(cont%chi)) deallocate(cont%chi)
       !if (allocated(cont%U)) deallocate(cont%U)
	  !!if (allocated(cont%Jbar)) deallocate(cont%Jbar)
      if (associated(cont%atom)) NULLIFY(cont%atom)
     end do
      deallocate(atom%continua)
     end if

  RETURN
  END SUBROUTINE freeAtom


  SUBROUTINE realloc_transitions(atom, Nt, mask)
  !Only kept indexes in atom%at if line contribute to opac
   integer, intent(in) :: Nt
   type (AtomType), intent(inout) :: atom
   logical, intent(in), dimension(:) :: mask
   integer :: k, k1
   type (AtomicTransition), dimension(atom%Ntr) :: trans
   
    !temporary storage
    trans = atom%at

    deallocate(atom%at)
    allocate(atom%at(Nt));atom%Ntr=Nt

    k1=1
    do k=1,atom%Nline
     !write(*,*) trans(k)%ik, trans(k)%trtype
     if (mask(k)) then 
      atom%at(k1) = trans(k)
      !write(*,*) atom%at(k1)%ik, atom%at(k1)%trtype
      k1 = k1 + 1
     end if
    end do
    
    atom%Ntr_line = k1-1
    
    k1 = atom%Ntr_line + 1
    do k=Atom%Nline+1,atom%Nline+atom%Ncont
     !write(*,*) trans(k)%ik, trans(k)%trtype
     if (mask(k)) then 
      atom%at(k1) = trans(k)
      !write(*,*) atom%at(k1)%ik, atom%at(k1)%trtype
      k1 = k1 + 1
     end if
    end do
    
    !write(*,*) atom%Ntr_line, atom%Ntr, size(trans), size(atom%lines), size(atom%continua)
   k1 = 0
  RETURN
  END SUBROUTINE realloc_transitions
  

  SUBROUTINE realloc_line_transitions_deprec(atom, Nl_new, mask)
  !basicazlly does what PAck does
   integer, intent(in) :: Nl_new
   type (AtomType), intent(inout) :: atom
   logical, intent(in), dimension(:) :: mask
   !integer, intent(in), dimension(:) :: indexes
   integer :: k, k1
   type (AtomicLine), dimension(:), allocatable :: lines
   
    !temporary storage
    allocate(lines(atom%Nline)); lines=atom%lines
     do k=1,atom%Nline 
      if (allocated(atom%lines(k)%phi)) &
      	deallocate(atom%lines(k)%phi)
      if (allocated(atom%lines(k)%lambda)) &
      	deallocate(atom%lines(k)%lambda)
      if (allocated(atom%lines(k)%CoolRates_ij)) &
      	deallocate(atom%lines(k)%CoolRates_ij)
      !if (allocated(atom%lines(k)%wlam)) &
      !	deallocate(atom%lines(k)%wlam)
      if (allocated(atom%lines(k)%phi)) &
      	deallocate(atom%lines(k)%phi)
      if (allocated(atom%lines(k)%rho_pfr)) &
      	deallocate(atom%lines(k)%rho_pfr)
      !if (allocated(atom%lines(k)%gij)) &
      !	deallocate(atom%lines(k)%gij)
      !if (allocated(atom%lines(k)%Vij)) &
      !	deallocate(atom%lines(k)%Vij)

      !!if (allocated(atom%lines(k)%Jbar)) &
      	!!deallocate(atom%lines(k)%Jbar)
      if (associated(atom%lines(k)%atom)) &
      	NULLIFY(atom%lines(k)%atom)
      CALL freeZeemanMultiplet(atom%lines(k))
     end do
    !Pack the transitions
    deallocate(atom%lines)
    allocate(atom%lines(Nl_new));atom%Nline=Nl_new
    k1=1
    do k=1,size(lines)
     if (mask(k)) then 
      atom%lines(k1) = lines(k)
      k1 = k1 + 1
     end if
    end do
    deallocate(lines)
  
  RETURN
  END SUBROUTINE realloc_line_transitions_deprec
  
  SUBROUTINE realloc_continuum_transitions_deprec(atom, Nc_new, mask)
   integer, intent(in) :: Nc_new
   logical, intent(in), dimension(:) :: mask
   type (AtomType), intent(inout) :: atom
   !integer, intent(in), dimension(:) :: indexes
   integer :: k, k1
   type (AtomicContinuum), dimension(:), allocatable :: conta


    allocate(conta(atom%Ncont)); conta=atom%continua
    do k=1,atom%Ncont
     if (allocated(atom%continua(k)%lambda)) &
     	deallocate(atom%continua(k)%lambda)
     if (allocated(atom%continua(k)%alpha)) &
     	deallocate(atom%continua(k)%alpha)
     if (allocated(atom%continua(k)%coolrates_ij)) &
     	deallocate(atom%continua(k)%coolrates_ij)
     !need to be reallocated in case of active atoms
     !if (allocated(atom%continua(k)%gij)) &
     !	deallocate(atom%continua(k)%gij)
     !if (allocated(atom%continua(k)%Vij)) &
     !	deallocate(atom%continua(k)%Vij)
	 !!if (allocated(atom%continua(k)%Jbar)) &
	    !!deallocate(atom%continua(k)%Jbar)
     if (associated(atom%continua(k)%atom)) &
     	NULLIFY(atom%continua(k)%atom)
    end do
    !Pack the transitions
    deallocate(atom%continua)
    allocate(atom%continua(Nc_new));atom%Ncont=Nc_new
    k1=1
    do k=1,size(conta)
     !write(*,*) size(conta), k, k1, mask(k)
     if (mask(k)) then 
      atom%continua(k1) = conta(k)
      !write(*,*) maxval(atmos%Atoms(n)%ptr_atom%continua(k1)%alpha), maxval(conta(k)%alpha)
      k1 = k1 + 1
     end if
    end do
    deallocate(conta)

  RETURN
  END SUBROUTINE realloc_continuum_transitions_deprec

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
     if (allocated(atmos%Elements(n)%ptr_elem%mol_index)) deallocate (atmos%Elements(n)%ptr_elem%mol_index)
     if (allocated(atmos%Elements(n)%ptr_elem%ionpot)) deallocate (atmos%Elements(n)%ptr_elem%ionpot)
     if (allocated(atmos%Elements(n)%ptr_elem%pf)) deallocate (atmos%Elements(n)%ptr_elem%pf)
     if (allocated(atmos%Elements(n)%ptr_elem%n))  deallocate (atmos%Elements(n)%ptr_elem%n)
     !if alias exist for this element free it.
     if (associated(atmos%Elements(n)%ptr_elem%model)) NULLIFY(atmos%Elements(n)%ptr_elem%model)
!       ! de-alias %model to the proper atmos%Atoms(nmet), atmos%Atoms will be freed after.
    end do
   deallocate(atmos%Elements) !here it works
  end if
    
  if (allocated(atmos%nHtot)) deallocate(atmos%nHtot)
  if (allocated(atmos%ne)) deallocate(atmos%ne)
  if (allocated(atmos%Tpf)) deallocate(atmos%Tpf)
  if (allocated(atmos%T)) deallocate(atmos%T)
  !if (allocated(atmos%vturb)) deallocate(atmos%vturb)
  if (allocated(atmos%nHmin)) deallocate(atmos%nHmin)
  if (allocated(atmos%icompute_atomRT)) deallocate(atmos%icompute_atomRT)
!   if (allocated(atmos%lcompute_atomRT)) deallocate(atmos%lcompute_atomRT)
!   if (allocated(atmos%ldark_zone)) deallocate(atmos%ldark_zone)
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
  open(unit=1, file=TRIM(ABUNDANCE_FILE),status="old")!trim(mcfost_utils)//TRIM(ABUNDANCE_FILE)
  do while (EOF.eq.0)
   read(1, '(1A4, 1F5.3)', IOSTAT=EOF) charID, A
   if (charID(1:1) /= '#') then
    allocate(atmos%Elements(n)%ptr_elem)
    atmos%Elements(n)%ptr_elem%weight=atomic_weights(n)
    atmos%Elements(n)%ptr_elem%ID = charID(1:2) !supposed to be lowercases!!!!!
    atmos%Elements(n)%ptr_elem%abund = 10.**(A-12.)
    if (A <= -99) atmos%Elements(n)%ptr_elem%abund = 0d0
    !write(*,*) " abund=",atmos%Elements(n)%abund, A
    atmos%Elements(n)%ptr_elem%abundance_set=.true.

!     write(*,*) "Element ", atmos%Elements(n)%ptr_elem%ID," has an abundance of A=",&
!               atmos%Elements(n)%ptr_elem%abund," and a weight of w=",&
!               atmos%Elements(n)%ptr_elem%weight
!     stop
    if (n.gt.1 .and. atmos%metallicity.gt.0.) then
     atmos%Elements(n)%ptr_elem%abund = atmos%Elements(n)%ptr_elem%abund * &
                             10.**(atmos%metallicity)
    end if
    atmos%totalAbund = atmos%totalAbund+atmos%Elements(n)%ptr_elem%abund
    atmos%avgWeight = atmos%avgWeight+atmos%Elements(n)%ptr_elem%abund* &
                     atmos%Elements(n)%ptr_elem%weight
    !write(*,*) "updating total Abundance = ",atmos%totalAbund
    !write(*,*) "updating average weight = ", atmos%avgWeight

    !attribute partition function to each element
    !the temperature grid is shared and saved in atmos%Tpf
    !write(*,*) "Fill partition function for element ", atmos%Elements(n)%ID
    CALL readKurucz_pf(n, atmos%Npf, &
              atmos%Elements(n)%ptr_elem%Nstage, atmos%Elements(n)%ptr_elem%ionpot, &
              atmos%Elements(n)%ptr_elem%pf)
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
  write(*,*) "Total Abundance in the atmosphere = ", atmos%totalAbund
  write(*,*) "Total average weight = ", atmos%avgWeight
  write(*,*) "Weight per Hydrogen = ", atmos%wght_per_H
  write(*,*) ""
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
    atmos%Vxyz(:,:) = 0d0
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
   atmos%T(:) = 0d0


   if (.not.allocated(atmos%nHtot)) allocate(atmos%nHtot(atmos%Nspace))
   !if (.not.allocated(atmos%vturb)) allocate(atmos%vturb(atmos%Nspace))
   if (.not.allocated(atmos%ne)) allocate(atmos%ne(atmos%Nspace))
   !if (.not.allocated(atmos%nHmin)) allocate(atmos%nHmin(atmos%Nspace)) temporary

   
   atmos%ne(:) = 0d0
   !atmos%nHmin(:) = 0d0 !temporary
   atmos%nHtot(:) = 0d0
   !atmos%vturb(:) = 0d0 !m/s

   CALL fillElements()
   
!    allocate(atmos%lcompute_atomRT(atmos%Nspace))
   allocate(atmos%icompute_atomRT(atmos%Nspace))
   atmos%icompute_atomRT(:) = 0 !everything transparent at init.

   RETURN
   END SUBROUTINE init_atomic_atmos

!to do, this routine should change according to mcfost denity array
   SUBROUTINE define_atomRT_domain(itiny_T, itiny_nH)
   ! Set where to solve for the RT equation: where a cell is not 
   ! transparent = where there is a significant density and temperature.
   ! Determines also if we have to force electron density calculation
   ! Might be used to tell where we solve the atomic line RT in case of a wind, magnetospheric
   !protoplanetary disk model.
   ! It is also set dark zones.
   ! icompute_atomRT = 0 -> transparent (rho, T <= tiny_H, tiny_T)
   ! icompute_atomRT = 1 -> filled (rho, T > tiny_H, tiny_T)
   ! icompute_atomRT = -1 -> dark (rho goes to infinity and T goes to 0)
    integer :: icell
    double precision, optional :: itiny_nH, itiny_T
    double precision :: Vchar=0d0, tiny_nH=1d0, tiny_T=1d2
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
     write(*,*) "changing the value of tiny_T = ", tiny_T," to", 2d2
     tiny_T = 2d0
    end if
    if (tiny_nH <= 0) then
     write(*,*) "changing the value of tiny_nH = ", tiny_nH," to", 1d0
     tiny_nH = 1d0
    end if
    
    !atmos%lcompute_atomRT = (atmos%nHtot > tiny_nH) .and. (atmos%T > tiny_T) 
    ! atmos%icompute_atomRT(:) = 0 !transparent !already init
    where((atmos%nHtot > tiny_nH) .and. (atmos%T > tiny_T))
    	atmos%icompute_atomRT = 1
    end where
!     do icell=1,atmos%Nspace
! !      atmos%lcompute_atomRT(icell) = &
! !        (atmos%nHtot(icell) > tiny_nH) .and. (atmos%T(icell) > tiny_T)
!     if ((atmos%nHtot(icell) > tiny_nH) .and. (atmos%T(icell) > tiny_T)) &
!      	atmos%icompute_atomRT(icell) = 1 !filled
!     end do

   RETURN
   END SUBROUTINE define_atomRT_domain
   
   SUBROUTINE init_magnetic_field()
   ! ----------------------------------------- !
   ! Allocate space for magnetic field
   ! ----------------------------------------- !
    !if (.not.atmos%magnetized) RETURN
    allocate(atmos%Bxyz(atmos%Nspace, 3))
    atmos%Bxyz = 0.0
   RETURN
   END SUBROUTINE init_magnetic_field
   

!building
   FUNCTION B_project(icell, x,y,z,u,v,w, gamma, chi) result(Bmodule)
   ! ------------------------------------------- !
   ! Returned the module of the magnetic field at
   ! one point of the model icell and the angles
   ! gamma and chi.
   ! Gamma is the angle between B and the los
   ! chi is the angle between Bperp and ex
   ! Bperp is the projection of B onto a plane
   ! perpendicular to the los
   ! ------------------------------------------- !
    integer :: icell
    double precision :: x, y, z, u, v, w, bproj, Bmodule
    double precision :: r, bx, by, bz
    double precision, intent(inout) :: gamma, chi
    
     bproj = 0d0 
     Bmodule = 0d0!Module of B.
     gamma = 0d0; chi = 0d0
     if (lVoronoi) then !Bxyz in cartesian
       bproj = atmos%Bxyz(icell,1)*u + atmos%Bxyz(icell,2)*v + &
       					   atmos%Bxyz(icell,3) * w
       Bmodule = dsqrt(sum(atmos%Bxyz(icell,:)**2))
       gamma = acos(bproj/Bmodule)
       chi = acos(u*sin(gamma))
     else
      if (lmagnetoaccr) then
       r = dsqrt(x**2 + y**2)
       if (r <= tiny_dp) RETURN
       bx = x/r * atmos%Bxyz(icell,1) !R
       by = y/r * atmos%Bxyz(icell,1) !theta
       bz = atmos%Bxyz(icell,2) !z
       if (z<0) bz = -bz
       !module
       Bmodule = dsqrt(bx**2 + by**2 +bz**2)
       bproj = bx * u + by * v * bz * w
       gamma = acos(bproj / Bmodule)
       chi = acos(u*sin(gamma))
      else
       !temporary, to work with uniform
        Bmodule = dsqrt(sum(atmos%Bxyz(icell,:)**2))
        r = dsqrt(x**2 + y**2 + z**2)
        bx = x/r * atmos%Bxyz(icell,1) !x
        by = y/r * atmos%Bxyz(icell,2) !y
        bz = z/r * atmos%Bxyz(icell,3) !z
        bproj = bx*u + by*v + bz*w
        gamma = acos(bproj / Bmodule)
        chi = acos(bx*sin(gamma))
       !CALL Error("Geometry for magnetic field projection unkown")  
      end if     
     end if
    
    RETURN
   END FUNCTION B_project
   
!    !!should be moved in Atom_type
!    FUNCTION atomZnumber_old(atomID) result(Z)
!    --------------------------------------
!    return the atomic number of an atom
!    with ID = atomID.
!    Hydrogen is 1
!    --------------------------------------
!     character(len=ATOM_ID_WIDTH) :: atomID
!     integer :: Z, i
! 
!     Z = 1
!     do i=1,Nelem
!      if (atmos%Elements(i)%ID.eq.atomID) then
!        Z = i
!        exit
!      end if
!     do while (atmos%Elements(Z)%ptr_elem%ID.ne.atomID)
!      Z=Z+1
!     end do
! 
!    RETURN
!    END FUNCTION atomZnumber_old
   
  SUBROUTINE write_atmos_domain()
 ! ------------------------------------ !
 ! ------------------------------------ !
  use fits_utils, only : print_error
  integer :: unit, EOF = 0, blocksize, naxes(3), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1

  CALL ftinit(unit,"atomRT_domain.fits.gz",blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .false.
  bitpix = 64

  if (lVoronoi) then
   naxis = 1
   naxes(1) = atmos%Nspace 
   nelements = naxes(1)
  else
   if (l3D) then
    naxis = 3
    naxes(1) = n_rad
    naxes(2) = 2*nz
    naxes(3) = n_az
    nelements = naxes(1) * naxes(2) * naxes(3)
   else
    naxis = 2
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2)
   end if
  end if

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "integer", "0=empty,1=filled;-1=dark", ' ', EOF)
  !write data
  CALL ftpprj(unit,group,fpixel,nelements,atmos%icompute_atomRT,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
   
  
  RETURN
  END SUBROUTINE write_atmos_domain
   
  SUBROUTINE writeVfield()
 ! ------------------------------------ !
 ! ------------------------------------ !
  use fits_utils, only : print_error
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,"Vfield.fits.gz",blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .false.
  bitpix = -64  

  if (lVoronoi) then
   naxis = 2
   naxes(1) = atmos%Nspace 
   naxes(2) = 3 
   nelements = naxes(1) * naxes(2)
  else
   if (l3D) then
    naxis = 4
    naxes(4) = 3
    naxes(1) = n_rad
    naxes(2) = 2*nz
    naxes(3) = n_az
    nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
   else
    naxis = 3 !2 + 1 for all compo
    naxes(3) = 3
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2) * naxes(3)
   end if
  end if

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "m.s^-1", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,atmos%Vxyz,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
   
  
  RETURN
  END SUBROUTINE writeVfield
  

  SUBROUTINE readAtmos_ascii(filename)
  ! ------------------------------------------- !
   ! Read from ascii file a model to be used
   ! with the spherical grid of mcfost.
   ! Suppose that the model read is computed on 
   ! a predefined grid by mcfost.
  ! ------------------------------------------- !
   use getline
   use constantes
   use grid, only : cell_map
   character(len=*), intent(in)	:: filename
   real(kind=dp) Tshk!, Oi, Oo
   real(kind=dp) :: Tring, thetai, thetao
   integer, parameter :: Nhead = 3 !Add more
   character(len=MAX_LENGTH) :: rotation_law
   integer :: icell, Nread, syst_status, N_points, k, i, j, acspot
   character(len=MAX_LENGTH) :: inputline, FormatLine, cmd
   logical :: accretion_spots
   real :: south
   real, parameter :: Lextent = 1.5 !vchar=Lextent * vchar
   real(kind=dp), dimension(n_cells) :: rr, zz, pp
   logical :: is_not_dark
   real(kind=dp) :: rho_to_nH, Lr, rmi, rmo, Mdot = 1d-7, tc, phic
   
   lmagnetoaccr = .false.
   lVoronoi = .false.
   !read from file the velocity law if any
   CALL init_atomic_atmos()
   !For now we do not necessarily need to compute B if no polarization at all
   atmos%magnetized = (lmagnetic_field) .and. (PRT_SOLUTION /= "NO_STOKES")
   if (atmos%magnetized) CALL init_magnetic_field()
   
   rho_to_nH = 1d3 /masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3

   
   write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
   
   !could add header with magnetic field and so on
   !location of spots + lmagnetoaccretion flags if other kind of models with the use of spots
   ! + Tschok
   
   cmd = "wc -l "//trim(filename)//" > ntest.txt"
   call appel_syst(cmd,syst_status)
   open(unit=1,file="ntest.txt",status="old")
   read(1,*) N_points
   close(unit=1)
   !-2 because 2 headers lines
   write(*,*) "Found ", N_points - Nhead, " points and grid has", n_cells, " points"
   if (N_points - Nhead/= n_cells) then
    write(*,*) "Should read a model for the exact same grid as mcfost !"
    stop
   end if

   open(unit=1,file=filename, status="old")
   CALL getnextline(1, "#", FormatLine, inputline, Nread)
   read(inputline(1:Nread),*) rotation_law
   
   if (.not.lstatic) then
    SELECT CASE (rotation_law)
     CASE ("magneto-accretion")
     	lmagnetoaccr = .true.
     	write(*,*) " Velocity law is ", trim(rotation_law), lmagnetoaccr
     CASE ("None")
     	write(*,*) " No interpolation of the velocity field", trim(rotation_law)
     	write(*,*) " Profile and v_proj not modified to use that with not Voronoi grid yet"
     	stop
     CASE DEFAULT
        write(*,*) " Velocity law ", rotation_law," not handled yet"
        stop
    END SELECT
   end if
   
   !read T shock and if accretion spots
   CALL getnextline(1, "#", FormatLine, inputline, Nread)
   read(inputline(1:Nread),*) Tshk, acspot
   accretion_spots = .false.
   if (acspot==1) accretion_spots = .true.
   !now read thetao and thetai
   CALL getnextline(1, "#", FormatLine, inputline, Nread)
   read(inputline(1:Nread),*) thetai, thetao
   thetai = thetai * pi/180.
   thetao = thetao * pi/180.
      
   do i=1, n_rad
     do j=j_start,nz !j_start = -nz in 3D
      do k=1, n_az
       if (j==0) then !midplane
        !icell = cell_map(i,1,k)
        cycle
       else
         icell = cell_map(i,j,k)
       end if
    Nread = 0
    if (atmos%magnetized) then
     write(*,*) "Magnetic atmos not ready yet"
     stop
     CALL getnextline(1, "#", FormatLine, inputline, Nread)
     read(inputline(1:Nread),*) rr(icell), zz(icell), pp(icell), atmos%T(icell), atmos%nHtot(icell), atmos%ne(icell), & 
         atmos%Vxyz(icell,1), atmos%Vxyz(icell,2), atmos%Vxyz(icell,3), &
         atmos%Bxyz(icell,1), atmos%Bxyz(icell,2), atmos%Bxyz(icell,3), atmos%icompute_atomRT(icell)
    else
     CALL getnextLine(1, "#", FormatLine, inputline, Nread)
     read(inputline(1:Nread),*) rr(icell), zz(icell), pp(icell), atmos%T(icell), atmos%nHtot(icell), atmos%ne(icell), & 
         atmos%Vxyz(icell,1), atmos%Vxyz(icell,2), atmos%Vxyz(icell,3), atmos%icompute_atomRT(icell)
    end if
     end do
    end do
   end do

   !rho -> nH
   atmos%nHtot = atmos%nHtot * rho_to_nH
   close(unit=1)
   write(*,*) "Read ", size(pack(atmos%icompute_atomRT,mask=atmos%icompute_atomRT>0)), " density zones"
   write(*,*) "Read ", size(pack(atmos%icompute_atomRT,mask=atmos%icompute_atomRT==0)), " transparent zones"
   write(*,*) "Read ", size(pack(atmos%icompute_atomRT,mask=atmos%icompute_atomRT<0)), " dark zones"

   if (thetai-thetao /= 0.) then
   rmi = 1d0 / sin(thetai)**2
   rmo = 1d0 / sin(thetao)**2
!    if ((present(Oi)) .and. (present(Oo))) then
!     rmi = 1d0 / sin(Oi * PI/180.)**2
!     rmo = 1d0 / sin(Oo * PI/180.)**2
!     thetai = Oi * PI/180.!asin(dsqrt(1d0/rmi))
!     thetao = Oo * PI/180 !asin(dsqrt(1d0/rmo))
!    else
!     rmi = 2d0
!     rmo = 3d0
!     thetai = asin(dsqrt(1d0/2d0))
!     thetao = asin(dsqrt(1d0/3d0))   
!    end if
    Lr = Ggrav * etoile(1)%M * Msun_to_kg * Mdot / (etoile(1)%r*au_to_m) * (1. - 2./(rmi + rmo))   
   
   
    write(*,*) "Angles at stellar surface (deg)", thetao*rad_to_deg, thetai*rad_to_deg
    Tring = Lr / (4d0 * PI * (etoile(1)%r*au_to_m)**2 * sigma * abs(cos(thetai)-cos(thetao)))
    Tring = Tring**0.25
   
   
    if (Tshk > 0) Tring = Tshk

   
   !2 is for two rings, 2Pi is dphi from 0 to 2pi / total sphere surface
    write(*,*) "Ring T ", Tring, "K"
    write(*,*) "Surface ", 100*(abs(cos(thetai)-cos(thetao))), "%" !couverte par les deux anneaux sur toute la sphere
    write(*,*) "Luminosity", Lr/Lsun, "Lsun"
    
    

    !Add ring
     if (accretion_spots) then
      etoile(1)%Nr = 2
   	  allocate(etoile(1)%SurfB(etoile(1)%Nr))
   	  do k=1,etoile(1)%Nr!should be 1 ring for testing
   	    south = 1.
    	tc = 0d0 * PI/180 !center of vector position
    	phic = 0d0 !vector position, pointing to the center of the spot
    	if (k==2) south = -1. !for the ring in the southern hemisphere
    	etoile(1)%SurfB(k)%T = Tring
        !center of the spot
   	    etoile(1)%SurfB(k)%r(1) = cos(phic)*sin(tc)
   	    etoile(1)%SurfB(k)%r(2) = sin(phic)*sin(tc)
   	    etoile(1)%SurfB(k)%r(3) = cos(tc)*south
    	etoile(1)%SurfB(k)%muo = cos(thetao); etoile(1)%SurfB(k)%mui = cos(thetai)
    	etoile(1)%SurfB(k)%phio = 2*PI; etoile(1)%SurfB(k)%phii = 0d0
  	  end do
  	 endif    
    
   end if !thetai-thetao /= 0
   write(*,*) "Mean molecular weight", atmos%avgWeight


   if (.not.lstatic) then
    atmos%v_char = atmos%v_char + maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)),&
    	   dim=1,mask=sum(atmos%Vxyz**2,dim=2)>0)
   end if
   
   if (atmos%magnetized) then
    atmos%B_char = maxval(sum(atmos%Bxyz(:,:)**2,dim=2)) !* Lextent
    write(*,*)  "Typical Magnetic field modulus (G)", atmos%B_char * 1d4
   end if

   !!here not used, done by model
  !!CALL define_atomRT_domain() !can be defined earlier so that we avoid create a dark_zone array?	
  							  !for instance, cell by cell init of icompute_atomRT, when density
  							  !is computed
  !but this one is needed
   if (maxval(atmos%ne) == 0d0) atmos%calc_ne = .true.

  !no need if we do not the dark_zones from input file.
   CALL write_atmos_domain() !but for consistency with the plot functions in python
   
   if (.not.lstatic) then
    write(*,*) "Maximum/minimum velocities in the model (km/s):"
    write(*,*) maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1)*1d-3,  &
     		  minval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1,&
     		  mask=sum(atmos%Vxyz**2,dim=2)>0)*1d-3
   end if
   write(*,*) "Typical velocity in the model (km/s):"
   atmos%v_char = Lextent * atmos%v_char
   write(*,*) atmos%v_char/1d3

   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%icompute_atomRT>0)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%icompute_atomRT>0)  

  RETURN
  END SUBROUTINE readAtmos_ascii
  
  SUBROUTINE writeAtmos_ascii()
  ! ------------------------------------- !
   ! Write GridType atmos to ascii file.
   ! ultimately this model will be mapped
   ! onto a Voronoi mesh.
  ! ------------------------------------- !
   use grid
   character(len=10)	:: filename="model.s"
   integer :: icell, i, j, k
   double precision, dimension(n_cells) :: XX, YY, ZZ
   double precision :: r, rcyl, phi, z, sinTheta, rho_to_nH, fact
   double precision :: dx, dy, dz, smooth_scale
   
   rho_to_nH = 1d0 / AMU /atmos%avgWeight
   
   if (n_az <= 1) then
    write(*,*) "Error, in this case, lmagnetoaccr is false, and", &
    	' n_az > 1'
    stop
   end if
   !X, Y, Z cartesian coordinates
   do i=1, n_rad
    do j=j_start,nz
     do k=1, n_az
      if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
      else
       icell = cell_map(i,j,k)
       rcyl = r_grid(icell)
       z = z_grid(icell)/z_scaling_env
      end if
      phi = phi_grid(icell)
      r = dsqrt(z**2 + rcyl**2)
      fact = 1d0
      if (z<0) fact = -1
      sinTheta = fact*dsqrt(1.-(z/r)**2)
      !XX(icell) = rcyl*cos(phi)*AU_to_m!r*cos(phi)*sinTheta * AU_to_m
      !YY(icell) = rcyl*sin(phi)*AU_to_m!r*sin(phi)*sinTheta * AU_to_m
       XX(icell) = r*cos(phi)*sinTheta * AU_to_m
       YY(icell) = r*sin(phi)*sinTheta * AU_to_m   
       ZZ(icell) = z*AU_to_m
     end do
    end do
   end do
   
   if (lstatic) then 
    allocate(atmos%Vxyz(n_cells,3))
    atmos%Vxyz = 0d0
   end if
   
   open(unit=1,file=trim(filename),status="replace")
   atmos%nHtot = atmos%nHtot / rho_to_nH !using total density instead of nHtot for writing
   !write Masses per cell instead of densities. Voronoi cells Volume is different than
   !grid cell right ? 
   write(*,*) "Maximum/minimum density in the model (kg/m^3):"
   write(*,*) maxval(atmos%nHtot), minval(atmos%nHtot,mask=atmos%icompute_atomRT>0)
   atmos%nHtot = atmos%nHtot * volume * AU3_to_m3  * kg_to_Msun!Msun
   do icell=1, atmos%Nspace
     dx = XX(icell) - etoile(1)%x*AU_to_m
     dy = YY(icell) - etoile(1)%y*AU_to_m
     dz = ZZ(icell) - etoile(1)%z*AU_to_m
     smooth_scale = 5. * Rmax * AU_to_m
     if (min(dx,dy,dz) < 2*etoile(1)%r*AU_to_m) then
      if ((dx**2+dy**2+dz**2) < (2*etoile(1)%r*AU_to_m)**2) then
         smooth_scale = dsqrt(dx**2+dy**2+dz**2)/3. * AU_to_m
         !=etoile(1)%r*AU_to_m!etoile(1)%r/3. * AU_to_m!m
      end if
     end if
        write(1,"(8E)") XX(icell), YY(icell), ZZ(icell), atmos%nHtot(icell), &
    	 atmos%Vxyz(icell,1), atmos%Vxyz(icell,2), atmos%Vxyz(icell,3), smooth_scale
   end do
   write(*,*) "Maximum/minimum mass in the model (Msun):"
   write(*,*) maxval(atmos%nHtot), minval(atmos%nHtot,mask=atmos%icompute_atomRT>0)
   write(*,*) "Maximum/minimum mass in the model (Kg):"
   write(*,*) maxval(atmos%nHtot)*Msun_to_kg, &
   				minval(atmos%nHtot,mask=atmos%icompute_atomRT>0)*Msun_to_kg
   write(*,*) "Maximum/minimum velocities in the model (km/s):"
   write(*,*) maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1)*1d-3,  &
     		  minval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1,&
     		  mask=sum(atmos%Vxyz**2,dim=2)>0)*1d-3
   close(unit=1)

  RETURN
  END SUBROUTINE writeAtmos_ascii

END MODULE atmos_type

