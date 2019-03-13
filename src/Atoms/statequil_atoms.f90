MODULE statequil_atoms

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type, only : NLTEspec
 use constant

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


 
 SUBROUTINE fillCrossCoupling_terms(id, icell)
 
  !previously allocated for this thread
  integer, intent(in) :: id, icell
  integer nact, i, j, kl
  type (AtomType), pointer :: atom
  type (AtomicContinuum) :: cont
  type (AtomicLine) :: line
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2, chicc
  
  do nact=1,atmos%Nactiveatoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   do kl=1,atom%Ncont
    cont = atom%continua(kl)
    twohnu3_c2 = 2d0 * HPLANCK * CLIGHT / (NM_TO_M * NLTEspec%lambda)**(3d0)
    i = cont%i; j = cont%j
    chicc = cont%Vij(:,id) * (atom%n(i,icell) - cont%gij(:,id)*atom%n(j,icell))
    
    atom%Uji_down(j,:,id) = atom%Uji_down(j,:,id) + twohnu3_c2*cont%gij(:,id)*cont%Vij(:,id)
    atom%chi_up(i,:,id) = atom%chi_up(i,:,id) + chicc
    atom%chi_down(j,:,id) = atom%chi_down(j,:,id) + chicc
   end do

   do kl=1,atom%Nline
    line = atom%lines(kl)
    twohnu3_c2 = line%Aji/line%Bji
    i = line%i; j = line%j
    chicc = line%Vij(:,id) * (atom%n(i,icell) - line%gij(:,id)*atom%n(j,icell))

    atom%Uji_down(j,:,id) = atom%Uji_down(j,:,id) + twohnu3_c2*line%gij(:,id)*line%Vij(:,id)
    atom%chi_up(i,:,id) = atom%chi_up(i,:,id) + chicc
    atom%chi_down(j,:,id) = atom%chi_down(j,:,id) + chicc
   end do   
   NULLIFY(atom)
  end do
 RETURN
 END SUBROUTINE fillCrossCoupling_terms
 
 SUBROUTINE fillGamma(id,icell)
 !---------------------------------- !
  ! for all atoms
  !
 !---------------------------------- !
  integer, intent(in) :: id, icell
  integer :: nact, nc, kr
  type (AtomType), pointer :: atom

  

  NULLIFY(atom)
 RETURN
 END SUBROUTINE fillGamma
 
 SUBROUTINE statequil()

 RETURN
 END SUBROUTINE statequil
 
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
 
 SUBROUTINE addtoRadiativeRates()
 ! -------------------------------------------- !
  ! For each atoms and each transitions
  ! compute the radiative rates
  ! It is a integral over frequency and angles
  !
 ! -------------------------------------------- !
 
 
 RETURN
 END SUBROUTINE addtoRadiativeRates


END MODULE statequil_atoms
