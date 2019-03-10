MODULE Opacity

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type, only : NLTEspec

 IMPLICIT NONE

 !store the pure continuum NLTE opacities to be added to the total
 !continuum opacity after NLTEloop ends
 double precision, dimension(:,:), allocatable :: chic_nlte, etac_nlte
 !cross-coupling opacities
 double precision, dimension(:,:), allocatable :: chi_xc, eta_xc

 CONTAINS

 SUBROUTINE NLTEOpacity(id, icell, x, y, z, x1, y1, z1, u, v, w, l)
  !
  !
  ! chi = Vij * (ni - gij * nj)
  ! eta = twohnu3_c2 * gij * Vij * nj
  ! Continuum:
  ! Vij = alpha
  ! gij = nstari/nstarj * exp(-hc/kT/lamba)
  ! Lines:
  ! twoHnu3/c2 = Aji/Bji
  ! gij = Bji/Bij (*rho if exists)
  ! Vij = Bij * hc/4PI * phi
  !
  !
  integer, intent(in) :: id, icell
  double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
  integer :: nact, Nred, Nblue, kc, kr, i, j
  type(AtomicLine) :: line
  type(AtomicContinuum) :: cont
  type(AtomType) :: aatom
  double precision :: gij, Vij
  
  do nact = 1, atmos%Nactiveatoms
   aatom = atmos%ActiveAtoms(nact)%ptr_atom
   
   !Separated loop on b-f and b-b- transitions
   do kc = 1, aatom%Ncont
    cont = aatom%continua(kc)
    Nred = cont%Nred; Nblue = cont%Nblue
    i = cont%i; j=cont%j
    !gij for continuum is simply continuum%alpha
    !Vij =
    
    !NLTEspec%AtomOpac%chi(:,id) = 
    !Do not forget to add continuum opacities to the all continnum opacities
    !after all populations have been converged
    
   end do
   
   do kr = 1, aatom%Nline
    line = aatom%lines(kr)
    Nred = line%Nred; Nblue = line%Nblue
    i = line%i; j=line%j
   end do
  
  end do

 RETURN
 END SUBROUTINE NLTEOpacity


END MODULE Opacity
