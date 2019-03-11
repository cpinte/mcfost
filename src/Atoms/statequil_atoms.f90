MODULE statequil_atoms

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type, only : NLTEspec

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE initGamma()
 ! ------------------------------------------------------ !
  !for each active atom, allocate and init Gamma matrix
  !deallocated in freeAtom()
  ! Collision matrix has to be allocated
  ! if already allocated, simply set Gamma to its new icell
  ! value.
 ! ------------------------------------------------------ !
  integer :: nact, Nlevel, i, j, ij
  type (AtomType), pointer :: atom
  
  do nact = 1, atmos%Nactiveatoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   Nlevel = atom%Nlevel
   ij = 0
   if (.not.allocated(atom%Gamma)) allocate(atom%Gamma(Nlevel, Nlevel))
   do i=1, Nlevel
    do j=1,Nlevel
      ij = ij + 1
      atom%Gamma(i,j) = atom%C(ij) 
    end do
   end do
   NULLIFY(atom)
  end do
 
 RETURN
 END SUBROUTINE initGamma
 

 SUBROUTINE allocNetCoolingRates(atom)
 ! ----------------------------------------------- !
  ! For each transition of AtomType atom
  ! allocate cooling rates.
  ! It avoids to store 2 x Ntransitons * Nspace
  ! arrays of radiative rates.
  ! Instead Ntranstions * Nspace arrays of cooling
  ! rates
 ! ----------------------------------------------- !
  integer :: kc, kr
  type (AtomType), intent(inout) :: atom
   write(*,*) "Allocating net cooling rates..."
   write(*,*) "  -> continuum cooling rates not implemented yet"
!    do kc=1,atom%Ncont
!     allocate(atom%continua(kc)%CoolRates_ij(atmos%Nspace))
!     atom%continua(kc)%CoolRates_ij(:) = 0d0
!    end do
   do kr=1,atom%Nline
    allocate(atom%lines(kr)%CoolRates_ij(atmos%Nspace))
    atom%lines(kr)%CoolRates_ij(:) = 0d0
   end do 

 RETURN
 END SUBROUTINE allocNetCoolingRates
 
 SUBROUTINE initRadiativeRates(atom)
 ! ---------------------------------------- !
  ! For each transition of AtomType atom
  ! set radiative rates to 0
 ! ---------------------------------------- !
  integer :: kc, kr
  type (AtomType), intent(inout) :: atom
  
   do kc=1,atom%Ncont
    atom%continua(kc)%Rij = 0d0
    atom%continua(kc)%Rji = 0d0
   end do
   do kr=1,atom%Nline
    atom%lines(kr)%Rij = 0d0
    atom%lines(kr)%Rji = 0d0
   end do 
   
 RETURN
 END SUBROUTINE initRadiativeRates
 
 SUBROUTINE fillCrossCoupling_terms()
 
 RETURN
 END SUBROUTINE fillCrossCoupling_terms
 
 SUBROUTINE addtoRadiativeRates()
 ! -------------------------------------------- !
  ! For each atoms and each transitions
  ! compute the radiative rates
  ! It is a integral over frequency and angles
  !
 ! -------------------------------------------- !
 
 
 RETURN
 END SUBROUTINE addtoRadiativeRates

 SUBROUTINE statequil()

 RETURN
 END SUBROUTINE statequil


END MODULE statequil_atoms
