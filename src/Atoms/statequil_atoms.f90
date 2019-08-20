MODULE statequil_atoms

 use atmos_type, only : atmos, nTotal_atom
 use atom_type
 use spectrum_type, only : NLTEspec
 use constant
 use opacity
 use math, only : locate, any_nan_infinity_matrix
 use utils, only : GaussSlv
 use parametres
 use accelerate
 use collision, only : CollisionRate, Collision_Hydrogen
 use metal, only : bound_free_Xsection

 IMPLICIT NONE

 PROCEDURE (FillGamma_bf_Hjde), pointer :: FillGamma_bf => NULL()
 PROCEDURE (FillGamma_bb_Hjde), pointer :: FillGamma_bb => NULL()

 CONTAINS
 
 SUBROUTINE initGamma_atom(id, icell, atom)
  integer, intent(in) :: icell, id
  type (AtomType), pointer, intent(inout) :: atom
  
	 if (atom%ID=="H") then
	  atom%Gamma(:,:,id) = Collision_Hydrogen(icell)
	 else
      atom%Gamma(:,:,id) = CollisionRate(icell, atom) 
     end if
  
  RETURN 
 END SUBROUTINE initGamma_atom
 
  
 SUBROUTINE initGamma(id, icell)
 ! ------------------------------------------------------ !
  !for each active atom, allocate and init Gamma matrix
  !deallocated in freeAtom()
  ! Collision matrix has to be allocated
  ! if already allocated, simply set Gamma to its new icell
  ! value.
  !
  ! n(l)*Cl->lp = n(lp)*Clp->l
  ! e.g.,
  ! Cij = Cji * (nj/ni)^star, with Cji = C(ij) = 
  ! colision rate from j to i.
  !
  ! (ij = j->i = (i-1)*N + j)
  ! (ji = i->j = (j-1)*N + i)
 ! ------------------------------------------------------ !
  integer :: nact, Nlevel, lp, l, ij, ji, nati, natf
  integer, intent(in) :: icell, id
  type (AtomType), pointer :: atom
  double precision :: c_lte = 1d0
  
!   nati = (1. * (id-1)) / NLTEspec%NPROC * atmos%Nactiveatoms + 1
!   natf = (1. * id) / NLTEspec%NPROC * atmos%Nactiveatoms
  nati = 1; natf = atmos%Nactiveatoms
  do nact = nati, natf
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   
!    write(*,*) "icell in initG", icell, "id=", id!, nact, atom%ID
   Nlevel = atom%Nlevel
!   open(unit=12, file="Cji_Cij_H4x4.dat",status="unknown")

!    do lp=1,Nlevel
!     do l=1,Nlevel
!       ij = (l-1)*nlevel + lp!lp->l
!       ji = (lp-1)*nlevel + l!l->lp
!       !          col/row
!       atom%Gamma(lp,l,id) =  c_lte*atom%Ckij(icell,ij)!Gamma_llp = Gamma(lp, l) rate from lp to l
!       !!                                                      Cul                  Clu=Cul*nu/nl with Cul=Gul here
! !      write(12,'(6E)') real(lp), real(l), real(ij), real(ji), atom%Ckij(icell,ij), atom%Gamma(lp,l) * atom%nstar(lp, icell)/atom%nstar(l,icell)
!     end do
!    end do
	 CALL initGamma_atom(id, icell, atom)
! 	 if (atom%ID=="H") then
! 	  atom%Gamma(:,:,id) = Collision_Hydrogen(icell)
! 	 else
!       atom%Gamma(:,:,id) = CollisionRate(icell, atom) 
!      end if
     !write(*,*) id, icell, 'max,min Cul for atom', atom%ID, maxval(atom%Gamma(:,:,id)), minval(atom%Gamma(:,:,id))
     
! close(12)
   NULLIFY(atom)
  end do
 RETURN
 END SUBROUTINE initGamma
 
 SUBROUTINE fillGamma(id, icell, n_rayons, switch_to_lte)
  ! ------------------------------------------------------------------------- !
   ! Fill the rate matrix Gamma, whose elements are Gamma(l',l) is the rate
   ! of transition from level l' to l.
   ! At initialisation, Gamma(l',l) = C(J,I), the collisional rates from upper
   ! level j to lower level i.
   !
   !
   ! The Sum_l" represents a sommation over each column for the lth row.
   !
  ! ------------------------------------------------------------------------- !

  integer, intent(in) :: id, n_rayons, icell !here for debug not used
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact
  type (AtomType), pointer :: atom
  
  !this pointers functions depends on the choice of the methode
  do nact=1, atmos%NactiveAtoms
   atom=>atmos%ActiveAtoms(nact)%ptr_atom
   CALL FillGamma_bf(id, icell, atom, n_rayons, switch_to_lte)
   CALL FillGamma_bb(id, icell, atom, n_rayons, switch_to_lte)
   atom=>NULL()
  enddo

 RETURN
 END SUBROUTINE fillGamma

 ! for that atom
 SUBROUTINE FillGamma_bf_hjde(id, icell, atom, n_rayons, switch_to_lte)
  integer, intent(in) :: id, n_rayons, icell
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact, kr, i, j, Nblue, Nred, l
  !type (AtomType), pointer :: atom
  type (AtomicContinuum) :: cont
  real(kind=dp) :: Jnu
  real(kind=dp), allocatable, dimension(:,:) :: Ieff
  real(kind=dp), dimension(:), allocatable :: twohnu3_c2k, weight, gijk

  if (present(switch_to_lte)) then
   if (switch_to_lte) RETURN !Gamma initialized to Cul
  end if
   
   do kr=1,atom%Ncont
    cont = atom%continua(kr)
    i = cont%i; j = cont%j
    Nblue = cont%Nblue; Nred = cont%Nred
    
    allocate(Ieff(cont%Nlambda,n_rayons), twohnu3_c2k(cont%Nlambda), weight(cont%Nlambda), gijk(cont%Nlambda))
    !Accelerated I
    Ieff(:,:) = NLTEspec%I(Nblue:Nred,1:n_rayons,id)*dexp(-NLTEspec%dtau(Nblue:Nred,1:n_rayons,id)) + \
             NLTEspec%Psi(Nblue:Nred,1:n_rayons,id) * atom%eta(Nblue:Nred,1:n_rayons,id)

             
    !nstar(i)/nstar(j)*exp(-hnu/kT)
    gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell) * &
    		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
    !2hnu3/c2 = 2hc/lambda3
    twohnu3_c2k(:) = twohc / NLTEspec%lambda(Nblue:Nred)**(3d0)
    
    weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
    ! alpha_nu * 4pi / h * dnu / nu = domega dnu/hnu alpha_nu

    !explicit do loop
    do l=1,cont%Nlambda
     Jnu = sum(Ieff(l,1:n_rayons))/n_rayons
     !write(*,*) j, i, NLTEspec%lambda(Nblue+l-1), Jnu(l)
     atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu*weight(l)
     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + (Jnu+twohnu3_c2k(l))*gijk(l)*weight(l)  
     
    enddo
    
    deallocate(Ieff,twohnu3_c2k, weight, gijk)!, Jnu)
   end do 
!    if (any_nan_infinity_matrix(atom%Gamma(:,:,id))>0) then
!     write(*,*) "Error in Gamma_bf"
!     write(*,*) atom%Gamma(:,:,id)
!    end if

  
 RETURN
 END SUBROUTINE FillGamma_bf_hjde
 
 SUBROUTINE FillGamma_bb_hjde(id, icell, atom, n_rayons, switch_to_lte)
 !here, I develop eq. 21 of Uitenbroek 2001, and I substitue I with I*dexp(-dtau) + Psi * eta
 ! "Hogereijde-like"
  integer, intent(in) :: icell, id, n_rayons !icell is not need, except if we introduce Xcoupling terms here
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact, kr, i, j, Nblue, Nred, l, iray
  !type (AtomType), pointer :: atom
  type (AtomicLine) :: line
  real(kind=dp), dimension(:), allocatable :: Ieff, weight
  real(kind=dp) :: factor, gij, twohnu3_c2, Vij, norm, norm_test, WTest(1)

  if (present(switch_to_lte)) then
   if (switch_to_lte) RETURN !Gamma initialized to Cul
  end if

   do kr=1,atom%Nline
    line = atom%lines(kr)
    i = line%i; j = line%j
    Nblue = line%Nblue; Nred = line%Nred
    
    allocate(Ieff(line%Nlambda), weight(line%Nlambda))

    twohnu3_c2 = line%Aji / line%Bji 

    norm = 0d0
    do iray=1, n_rayons
     norm = norm + sum(line%phi(:,iray,id)*line_wlam(line))!/n_rayons !angle and frequency integrated
    enddo
    !the n_rayons simplify
    weight(:) =  line_wlam(line) / norm! / factor
!     write(*,*) id, j, i, maxval(weight), norm
    
    gij = line%Bji/line%Bij

    do iray=1, n_rayons

      Ieff(:) = NLTEspec%I(Nblue:Nred,iray,id)*dexp(-NLTEspec%dtau(Nblue:Nred,iray,id)) + \
             NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id)


     !write(*,*) "id=",id, iray, maxval(atom%eta(Nblue:Nred,iray,id))
     do l=1,line%Nlambda
     !write(*,*) icell, id, iray, l, Ieff(l), line%Bij, line%phi(l, iray, id), weight(l)
      atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + line%Bij*line%phi(l, iray,id)*Ieff(l)*weight(l)
      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%phi(l, iray, id)*weight(l)*(line%Aji+line%Bji*Ieff(l))
 
      
     enddo
    enddo

    deallocate(Ieff, weight)
   enddo !over lines
  
!    if (any_nan_infinity_matrix(atom%Gamma(:,:,id))>0) then
!     write(*,*) "Error in Gamma_bb"
!     write(*,*) atom%Gamma(:,:,id)
!     stop
!    end if
 RETURN 
 END SUBROUTINE FillGamma_bb_hjde
 
 
 !!But there is a problem of construction here, because I have line1%phi and line2%phi, they may
 !overlap, but line1%Nlambda /= line2%Nlambda, therefore the product of the two is complicated.
 !I need either to give up about this for the moment, or to find the region where the two lines overlap
 !and make the product.
 !I cannot fix the number of points of line especially if it overlap, because of the addition of grid
 !points or I have to change that? The bad solution is to store Xcoupling on the whole grid, with zero
 !outside the transition, then the multiplication is easy.
 SUBROUTINE FillGamma_Xcoupling(id, j, i, kr, atom)
 !Attempt to add cross-coupling terms to Gamma.
 !Cross coupling are of the form:
 !Sum_all_active_transitions A*chi*Uji
 !with A = 1 for j->i transition and A = -1 for i->i" transitions
 !with i (lower level of j) beeing the upper level (if any) of the transition i->i"
 !Uji and chi are either continuum of line depending on if j->i is a cont or a line.
 
 !This routine handles the case of i is the upper level of another transition (i is the lower level
 !of either a line or cont, and the Xcoupling from j->i computed in FillGamma_bb/bf_MALI())
 !but we do not know what is the type of this transition.
  integer, intent(in) :: id, i, j, kr
  integer :: krr, Ntrans
  type(AtomType), intent(inout) :: atom
  type(AtomicContinuum) :: cont
  type(AtomicLine) :: line
  
  Ntrans = 	atom%Ncont + atom%Nline
  
  do krr=1, Ntrans
   if (krr<= atom%Ncont) then !continua
    cont = atom%continua(krr)
    if (cont%j==i) then
!      atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + sum(cont%chi(:,id)*fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)*&
!        sum(NLTEspec%Psi(cont%Nblue:cont%Nred,:,id))*cont%U(:,id)) 
     end if
   else !lines
    line = atom%lines(krr)
    if (line%j==i) then
!       atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + sum(cont%chi(:,id)*fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)*&
!        sum(NLTEspec%Psi(cont%Nblue:cont%Nred,:,id))*cont%U(:,id))    
    end if
   end if
  
  end do !over all trans
 
 RETURN
 END SUBROUTINE FillGamma_Xcoupling
 
!Presently, the cross-coupling is only for the line with itself
 SUBROUTINE FillGamma_bb_mali(id, icell, atom, n_rayons, switch_to_lte)
  integer, intent(in) :: icell, id, n_rayons
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact, kr, i, j, Nblue, Nred, l, iray, krr, jp
  !type (AtomType), pointer :: atom
  type (AtomicLine) :: line
  real(kind=dp), dimension(:), allocatable :: Ieff, weight, xcc
  real(kind=dp) :: factor, gij, twohnu3_c2, Vij, norm, norm_test, WTest(1)

  if (present(switch_to_lte)) then
   if (switch_to_lte) RETURN !Gamma initialized to Cul
  end if

   do kr=1,atom%Nline
    line = atom%lines(kr)
    i = line%i; j = line%j
    Nblue = line%Nblue; Nred = line%Nred
    
    allocate(Ieff(line%Nlambda), weight(line%Nlambda), xcc(line%Nlambda))

    twohnu3_c2 = line%Aji / line%Bji 

    norm = 0d0
    do iray=1, n_rayons
     norm = norm + sum(line%phi(:,iray,id)*line_wlam(line))!/n_rayons !angle and frequency integrated
    enddo
    !the n_rayons simplify
    weight(:) =  line_wlam(line) / norm! / factor
!     write(*,*) id, j, i, maxval(weight), norm
    
    gij = line%Bji/line%Bij


    do iray=1, n_rayons!* dexp(-NLTEspec%dtau(Nblue:Nred,iray,id))
      Ieff(:) = NLTEspec%I(Nblue:Nred,iray,id) - &
      NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id)
  
     do l=1,line%Nlambda
      atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + line%Bij*line%phi(l, iray,id)*Ieff(l)*weight(l) 
      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%phi(l, iray, id)*weight(l)*(line%Aji+line%Bji*Ieff(l))
     
      !X coupling j->i
      !Partial coupling within a line only
      atom%Gamma(j,i,id) = atom%Gamma(j, i, id) - line%chi(l,iray,id)*NLTEspec%Psi(Nblue+l-1,iray,id) * &
      		line%U(l, iray, id) * weight(l) / hc
     enddo
    enddo
    deallocate(Ieff, weight, xcc)
   enddo !over lines  
     
 RETURN
 END SUBROUTINE FillGamma_bb_mali
 
 SUBROUTINE FillGamma_bf_mali(id, icell, atom, n_rayons, switch_to_lte)
  integer, intent(in) :: id, n_rayons, icell
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact, kr, i, j, Nblue, Nred, l, krr, jp, iray
  !type (AtomType), pointer :: atom
  type (AtomicContinuum) :: cont
  real(kind=dp) :: Jnu
  real(kind=dp), allocatable, dimension(:,:) :: Ieff
  real(kind=dp), dimension(:), allocatable :: twohnu3_c2k, weight, gijk

  if (present(switch_to_lte)) then
   if (switch_to_lte) RETURN !Gamma initialized to Cul
  end if
   
   do kr=1,atom%Ncont
    cont = atom%continua(kr)
    i = cont%i; j = cont%j
    Nblue = cont%Nblue; Nred = cont%Nred
    
    allocate(Ieff(cont%Nlambda,n_rayons), twohnu3_c2k(cont%Nlambda), weight(cont%Nlambda), gijk(cont%Nlambda))
    !Accelerated I
    !* dexp(-NLTEspec%dtau(Nblue:Nred,:,id))
    Ieff(:,:) = NLTEspec%I(Nblue:Nred,1:n_rayons,id) - &
      	NLTEspec%Psi(Nblue:Nred,1:n_rayons,id) * atom%eta(Nblue:Nred,1:n_rayons,id)

             
    !nstar(i)/nstar(j)*exp(-hnu/kT)
    gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell) * &
    		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
    !2hnu3/c2 = 2hc/lambda3
    twohnu3_c2k(:) = twohc / NLTEspec%lambda(Nblue:Nred)**(3d0)
    
    weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
    ! alpha_nu * 4pi / h * dnu / nu = domega dnu/hnu alpha_nu

    !explicit do loop
    do l=1,cont%Nlambda
     Jnu = sum(Ieff(l,:))/n_rayons
     !write(*,*) j, i, NLTEspec%lambda(Nblue+l-1), Jnu(l)
     atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu*weight(l)
     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + (Jnu+twohnu3_c2k(l))*gijk(l)*weight(l)
     
           !here I can add the Xcoupling from j to i, easy its the same type of this routine
     
     !!X coupling j->i
      
      atom%Gamma(j,i,id) = atom%Gamma(j, i, id) - cont%chi(l,id)*sum(NLTEspec%Psi(Nblue+l-1,1:n_rayons,id))/n_rayons * &
      		cont%U(l,id) * weight(l)   
    enddo
    
    deallocate(Ieff,twohnu3_c2k, weight, gijk)!, Jnu)
   end do !over cont

   
 RETURN
 END SUBROUTINE FillGamma_bf_mali
 
 
 SUBROUTINE FillGamma_bf_zero_radiation(id, icell, atom, n_rayons, switch_to_lte)
  integer, intent(in) :: id, n_rayons, icell
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact, kr, i, j, Nblue, Nred, l
  !type (AtomType), pointer :: atom
  type (AtomicContinuum) :: cont
  real(kind=dp), dimension(:), allocatable :: twohnu3_c2k, weight, gijk

  if (present(switch_to_lte)) then
   if (switch_to_lte) RETURN !Gamma initialized to Cul
  end if
  
  !do nact=1,atmos%NactiveAtoms !loop over each active atoms
  
   !atom => atmos%ActiveAtoms(nact)%ptr_atom
   
   do kr=1,atom%Ncont
    cont = atom%continua(kr)
    i = cont%i; j = cont%j
    Nblue = cont%Nblue; Nred = cont%Nred
    
    allocate(twohnu3_c2k(cont%Nlambda), weight(cont%Nlambda), gijk(cont%Nlambda))
       
    !nstar(i)/nstar(j)*exp(-hnu/kT)
    gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell) * &
    		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
    !2hnu3/c2 = 2hc/lambda3
    twohnu3_c2k(:) = twohc / NLTEspec%lambda(Nblue:Nred)**(3d0)
    
    weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
    ! alpha_nu * 4pi / h / n_rayons * dnu / nu = domega dnu/hnu alpha_nu

    !explicit do loop
    do l=1,cont%Nlambda
     !!atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + 0d0
     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + twohnu3_c2k(l)*gijk(l)*weight(l)
    enddo
    
    deallocate(twohnu3_c2k, weight, gijk)
   end do !over cont
   
   !atom => NULL()
  !enddo !over atom
  
 RETURN
 END SUBROUTINE FillGamma_bf_zero_radiation
 
 SUBROUTINE FillGamma_bb_zero_radiation(id, icell, atom, n_rayons, switch_to_lte)
  integer, intent(in) :: icell, id, n_rayons
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact, kr, i, j
  !type (AtomType), pointer :: atom
  type (AtomicLine) :: line

  if (present(switch_to_lte)) then
   if (switch_to_lte) RETURN !Gamma initialized to Cul
  end if
  
  !do nact=1,atmos%NactiveAtoms !loop over each active atoms
  
   !atom => atmos%ActiveAtoms(nact)%ptr_atom

   do kr=1,atom%Nline
    line = atom%lines(kr)
    i = line%i; j = line%j

    !!atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + 0d0
    atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%Aji !integration over domega/4PI dv phi = 1

   enddo !over lines
   !NULLIFY(atom)
  !end do !loop over atoms    
  
 RETURN 
 END SUBROUTINE FillGamma_bb_zero_radiation


 SUBROUTINE SEE_atom(id, icell,atom)
 ! --------------------------------------------------------------------!
  ! For atom atom solves for the Statistical Equilibrium Equations (SEE) 
  !
  ! We solve for :
  !  Sum_l' Gamma_l'l n_l' = 0 (m^-3 s^-1)
  !
  ! which is equivalent in matrix notation to:
  !
  ! GG(l,l') dot n_l' = 0 with l' is on the columns and l on the rows
  !
  ! In particular for a 2-level atom the two equations are:
  !
  ! n1*G_11 + n2*G_21 = 0
  ! n1*G12 + n2*G22 = 0,
  ! with one of these equations has to be replaced by
  ! n1 + n2 = N
  !
  ! For numerical stability, the row with the largest populations is 
  ! removed (i.e., it is replaced by the mass conservation law).
  !
  ! This matrix Gamma for atoms (or later molecules) is compared to 
  ! molecules A matrix in mcfost:
  ! diag(Gamma) = -diag(A); offdiag(Gamma) = -offdiag(A)
  ! A = transpose(Gamma)
 ! --------------------------------------------------------------------!
  
  integer, intent(in) :: icell, id
  type(AtomType), intent(inout) :: atom
  integer :: lp, imaxpop, l
  double precision, dimension(atom%Nlevel, atom%Nlevel) :: Aij
  
  !we need to do that here as Gamma is updated for each rays.
!!atom%Gamma(:,:,id)=atom%Gamma(:,:,id)*-1d0
   do l = 1, atom%Nlevel
    atom%Gamma(l,l,id) = 0d0
    atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!!atom%Gamma(l,l,id) = sum(abs(atom%Gamma(l,:,id))) !sum over rows for this column
   end do
  !write(*,*) atom%ID, id, icell, maxval(atom%Gamma(:,:,id)), minval(atom%Gamma(:,:,id))
  !write(*,*) id, icell, atom%n(:,icell)

  imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
  !write(*,*) "imaxpop", imaxpop, atom%n(imaxpop,icell)
  atom%n(:,icell) = 0d0
  atom%n(imaxpop,icell) = nTotal_atom(icell, atom)
  
  !Sum_l'_imaxpop * n_l' = N
  atom%Gamma(:,imaxpop,id) = 1d0 !all columns of the last row for instance
  !(G11 G21)  (n1)  (0)
  !(       ) .(  ) =( )
  !(1    1 )  (n2)  (N)
  
  !Y a peut être un transpose ici  par rapport à MCFOST, atom%Gamma.T ?

  Aij = transpose(atom%Gamma(:,:,id))
  !!Aij = atom%Gamma(:,:,id)
  CALL GaussSlv(Aij, atom%n(:,icell),atom%Nlevel)
  !!atom%n(:,icell) = atom%n(:,icell) * nTotal_atom(icell, atom)
  if (any_nan_infinity_matrix(Aij)>0) then
    write(*,*) atom%Gamma(:,:,id)
    write(*,*) id, icell, atom%n(:,icell)
    stop
  end if

 RETURN
 END SUBROUTINE SEE_atom
 
 SUBROUTINE updatePopulations(id, icell)
 ! --------------------------------------------------------------------!
  ! Performs a solution of SEE for each atom.
  !
  ! to do: implements Ng's acceleration iteration. 
 ! --------------------------------------------------------------------!

  integer, intent(in) :: id, icell
  type(AtomType), pointer :: atom
  integer :: nat, nati, natf
  logical :: accelerate = .false.
  double precision :: dM
  
  !nati = (1. * (id-1)) / NLTEspec%NPROC * atmos%Nactiveatoms + 1
  !natf = (1. * id) / NLTEspec%NPROC * atmos%Nactiveatoms
  nati = 1; natf = atmos%Nactiveatoms
  do nat=nati,natf !loop over each active atoms
   atom => atmos%ActiveAtoms(nat)%ptr_atom
   CALL SEE_atom(id, icell, atom)
   atom => NULL()
  end do
 
 RETURN
 END SUBROUTINE updatePopulations
 
 SUBROUTINE Gamma_LTE(id,icell)
  ! ------------------------------------------------------------------------- !
   ! Fill the rate matrix Gamma, whose elements are Gamma(lp,l) is the rate
   ! of transition from level lp to l.
   ! At initialisation, Gamma(lp,l) = C(J,I), the collisional rates from upper
   ! level j to lower level i.
   !
   ! This is the LTE case where Gamma(l',l) = C(l'l)
   ! Gamma(l',l) = Cl'l - delta(l,l')Sum_l" (Cll").
   !
   ! This Gamma is frequency and angle independent. 
   ! FOR ALL ATOMS
  ! ------------------------------------------------------------------------- !
  integer, intent(in) :: id, icell
  integer :: nact, kr, l, lp, nati, natf
  type (AtomType), pointer :: atom

  !nati = (1. * (id-1)) / NLTEspec%NPROC * atmos%Nactiveatoms + 1
  !natf = (1. * id) / NLTEspec%NPROC * atmos%Nactiveatoms
  nati = 1; natf = atmos%Nactiveatoms
  do nact=nati,natf !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
!    do l=1,atom%Nlevel
!     do lp=1,atom%Nlevel   
!       atom%Gamma(lp, l) = atom%Ckij(icell, (l-1)*atom%Nlevel+lp) !lp->l; C_kij = C(j->i)
!     end do
!   end do

   !fill the diagonal here, delta(l',l)
   !because for each G(i,i); Cii, Rii is 0.
   !and Gii = -Sum_j Cij + Rij = -sum_j Gamma(i,j)
   !diagonal of Gamma, Gamma(col,col) is sum_row Gamma(col,row)
   !G(1,1) = - (G12 + G13 + G14 ...)
   !G(2,2) = - (G21 + G23 + G24 ..) first index is for column and second row
   do l = 1, atom%Nlevel
    atom%Gamma(l,l,id) = 0d0
    atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
   end do

   
   NULLIFY(atom)
  end do !loop over atoms

 RETURN
 END SUBROUTINE Gamma_LTE
 
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
