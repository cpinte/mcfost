MODULE statequil_atoms

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type, only : NLTEspec
 use constant
 use constantes, only : AU_to_m
 use opacity, only : cont_wlam, NLTEopacity, calc_psi_operator
 use math, only : locate
 use utils, only : GaussSlv
 use grid, only : cross_cell

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE init_local_field_atom(id, icell, iray, x0, y0, z0, u0, v0, w0)
  ! ------------------------------------------------------------------------- !
   ! Computes local radiation field, keeping External radiation constant
   ! in case of sub-iterations are turned on.
   ! The local radiation field is proportional to Snu = eta / chi, for each
   ! atom, for each cell and thread.
  ! ------------------------------------------------------------------------- !  
  
  double precision, intent(in) :: x0, y0, z0, u0, v0, w0
  integer, intent(in) :: id, icell, iray
  double precision :: l_dum, l_c_dum, l_v_dum, x1, x2, x3, ds
  integer 		   :: n_c_dum
  !recompute opacity of this cell., but I need angles and pos...
  !NLTEspec%I not changed
  ! move to the cell icell.
  CALL cross_cell(x0,y0,z0, u0,v0,w0, icell, &
       						n_c_dum, x1,x2,x3, n_c_dum, l_dum, l_c_dum, l_v_dum)
  NLTEspec%AtomOpac%chi(:,id) = 0d0
  NLTEspec%AtomOpac%eta(:,id) = 0d0
  !set atom%eta to zero also
  CALL initCrossCoupling(id)
  !NOTE Zeeman opacities are not re set to zero and are accumulated
  !change that or always use FIELD_FREE
  !is equivalent to P(icell,id, iray) ?
  !Compute opacity, eta and chi for this cell in the direction u0, v0, w0
  CALL NLTEOpacity(id, icell, x0, y0, z0, x1, x2, x3, u0, v0, w0, l_dum, .true.)
  !last .true. is to compute atom%gij, atom%Uji,atom%Vij 
  CALL fillCrossCoupling_terms(id, icell)
  !il faut recalculer Psi dans les sous-iterations et Ieff si Hogereijde.
  !Sinon, seulement Psi depend de chi a besoin d'être mis à jour.
  ds = l_dum * AU_to_m
  !recompute Psi and eventually Ieff.
  CALL calc_psi_operator(id, icell, iray, ds)
  
 RETURN
 END SUBROUTINE init_local_field_atom
 
 
 SUBROUTINE initGamma(icell) !icell temporary
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
  integer :: nact, Nlevel, lp, l, ij, ji
  integer, intent(in) :: icell
  type (AtomType), pointer :: atom
  double precision :: c_lte = 1d0
  
  do nact = 1, atmos%Nactiveatoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   Nlevel = atom%Nlevel
   if (.not.allocated(atom%Gamma)) allocate(atom%Gamma(Nlevel, Nlevel))
!   open(unit=12, file="Cji_Cij_H4x4.dat",status="unknown")

   do lp=1,Nlevel
    do l=1,Nlevel
      ij = (l-1)*nlevel + lp!lp->l
      ji = (lp-1)*nlevel + l!l->lp
      !          col/row
      atom%Gamma(lp,l) = c_lte*atom%Ckij(icell,ij) !Gamma_llp = Gamma(lp, l) rate from lp to l
      !!                                                      Cul                  Clu=Cul*nu/nl with Cul=Gul here
!      write(12,'(6E)') real(lp), real(l), real(ij), real(ji), atom%Ckij(icell,ij), atom%Gamma(lp,l) * atom%nstar(lp, icell)/atom%nstar(l,icell)
    end do
   end do
! close(12)
   NULLIFY(atom)
  end do
 
 RETURN
 END SUBROUTINE initGamma

 SUBROUTINE initCrossCoupling(id)
  ! --------------------------------------------------- !
   ! Init the X coupling terms for each active atoms.
   ! they are used in the Rybicki Hummer MALI scheme,
   ! with full preconditioning (overlapping transitions
   ! and background continua)
   !
   ! init a cell level (when eval_operator is true) for
   ! thread id.
   ! For all wavelength
  ! --------------------------------------------------- !

  integer, intent(in) :: id
  integer :: nact
  
  do nact=1,atmos%Nactiveatoms
   atmos%ActiveAtoms(nact)%ptr_atom%Uji_down(:,:,id) = 0d0!0 if j<i
   atmos%ActiveAtoms(nact)%ptr_atom%chi_up(:,:,id) = 0d0
   atmos%ActiveAtoms(nact)%ptr_atom%chi_down(:,:,id) = 0d0
   
   !init eta at the same time
   atmos%ActiveAtoms(nact)%ptr_atom%eta(:,id) = 0d0
  end do
 
 RETURN
 END SUBROUTINE initCrossCoupling
 
 SUBROUTINE fillCrossCoupling_terms(id, icell)
  ! -------------------------------------------- !
   ! fill X coupling terms for all transitions
   ! of all active atoms.
  ! -------------------------------------------- !
  !previously allocated for this thread
  integer, intent(in) :: id, icell
  integer nact, i, j, kl, Nblue, Nred
  type (AtomType), pointer :: atom
  type (AtomicContinuum) :: cont
  type (AtomicLine) :: line
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2, chicc
  
  do nact=1,atmos%Nactiveatoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom

   do kl=1,atom%Ncont
    cont = atom%continua(kl)
    Nblue = cont%Nblue; Nred = cont%Nred
    twohnu3_c2 = 2d0 * HPLANCK * CLIGHT / (NM_TO_M * NLTEspec%lambda)**(3d0)
    i = cont%i; j = cont%j
    chicc = 0d0
    chicc(Nblue:Nred) = cont%Vij(:,id) * (atom%n(i,icell) &
    					-cont%gij(:,id)*atom%n(j,icell))
    					
    atom%Uji_down(j,Nblue:Nred,id) = atom%Uji_down(j,Nblue:Nred,id) +&
    				 twohnu3_c2(Nblue:Nred)*cont%gij(:,id)*cont%Vij(:,id)
    				 
    atom%chi_up(i,Nblue:Nred,id) = atom%chi_up(i,Nblue:Nred,id) + chicc(Nblue:Nred)
    atom%chi_down(j,Nblue:Nred,id) = atom%chi_down(j,Nblue:Nred,id) + chicc(Nblue:Nred)
   end do

   do kl=1,atom%Nline
    line = atom%lines(kl)
    twohnu3_c2 = line%Aji/line%Bji
    i = line%i; j = line%j
    Nblue = line%Nblue; Nred=line%Nred
    chicc = 0d0
    chicc(Nblue:Nred) = line%Vij(:,id) * (atom%n(i,icell) &
    					- line%gij(:,id)*atom%n(j,icell))

    atom%Uji_down(j,Nblue:Nred,id) = atom%Uji_down(j,Nblue:Nred,id) + &
    				twohnu3_c2(Nblue:Nred)*line%gij(:,id)*line%Vij(:,id)
    				
    atom%chi_up(i,Nblue:Nred,id) = atom%chi_up(i,Nblue:Nred,id) + chicc(Nblue:Nred)
    atom%chi_down(j,Nblue:Nred,id) = atom%chi_down(j,Nblue:Nred,id) + chicc(Nblue:Nred)
    
   end do   
   NULLIFY(atom)
  end do
 RETURN
 END SUBROUTINE fillCrossCoupling_terms

 SUBROUTINE fillGamma_Hogereijde(id, icell, iray, n_rayons)
  ! ------------------------------------------------------------------------- !
   ! Fill the rate matrix Gamma, whose elements are Gamma(l',l) is the rate
   ! of transition from level l' to l.
   ! At initialisation, Gamma(l',l) = C(J,I), the collisional rates from upper
   ! level j to lower level i.
   !
   ! This is the case of Hogeriejde, similar to molecular lines in MCFOST.
   ! The intensity I in the rate matrix elements: 
   !  Cl'l + Rl'l - delta(l,l') Sum_l" (Cllp" + Rll"); appearing in Rll' and
   !  Rll" is directly replaced by Iexp(-dtau) + (1-exp(-dtau))*Snu; 
   ! with Iexp(-dtau) = Ieff and (1-exp(-dtau))/chi = Psi by definition.
   ! I = Ieff + Psi * eta
   !
   ! The Sum_l" represents a sommation over each column for the lth row.
   !
  ! ------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell, iray, n_rayons
  integer :: nact, nc, kr, switch,i, j, Nblue, Nred, l
  type (AtomType), pointer :: atom
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2, Ieff
  double precision :: norm = 0d0
  double precision :: c_nlte = 0d0
  
  Ieff(:) = 0d0
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   !loop over transitions, b-f and b-b for these atoms
   !To do; define a transition_type with either cont or line
   Ieff(:) = NLTEspec%Ieff(:,iray,id) + NLTEspec%Psi(:,iray,id)! * atom%eta(:,id)
   do kr=1,atom%Ncont
    norm = 4d0*PI / HPLANCK / n_rayons! ?? 1/4PI * 4PI or *4PI alone
    
    i = atom%continua(kr)%i; j = atom%continua(kr)%j
    Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred 
    twohnu3_c2 = 2d0 * HPLANCK * CLIGHT / (NLTEspec%lambda*NM_TO_M)**3.
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%continua(kr)%Vij(:,id)*Ieff(Nblue:Nred)*cont_wlam(atom%continua(kr))) * norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%continua(kr)%Vij(:,id)*atom%continua(kr)%gij(:,id)*cont_wlam(atom%continua(kr))) * norm

   end do

   do kr=1,atom%Nline
    norm =  4d0 * PI / (CLIGHT * HPLANCK) / n_rayons! ?? 1/4PI * 4PI or *4PI alone
    
    i = atom%lines(kr)%i; j = atom%lines(kr)%j
    Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred 
    twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji 
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%lines(kr)%Vij(:,id)*Ieff(Nblue:Nred)*atom%lines(kr)%wlam(:)) * norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:)) * norm
   end do
   
   if (iray==n_rayons) then !otherwise we remove several times GammaDiag
   	do l = 1, atom%Nlevel
    	atom%Gamma(l,l) = -sum(atom%Gamma(l,:)) !sum over rows for this column
   	end do
!    	write(*,*) icell, id, "Gamma (i, j, G(i,j), G(i,i), G(j,i), G(j,j)" !Remember that even if no continuum transitions, Gamma is init to Cji
!      do i=1,atom%Nlevel
!      do j=1,atom%Nlevel
!      write(*,*) i, j, atom%Gamma(i,j), atom%Gamma(i,i),  atom%Gamma(j,i), atom%Gamma(j,j)
!      end do
!     end do
! !     stop
   end if
   
   NULLIFY(atom)
  end do !loop over atoms  

 RETURN
 END SUBROUTINE fillGamma_Hogereijde
 
 SUBROUTINE FillGamma_ZeroRadiation(id, icell)
  ! ------------------------------------------------------------------------- !
   ! Fill the rate matrix Gamma, whose elements are Gamma(lp,l) is the rate
   ! of transition from level lp to l.
   ! At initialisation, Gamma(lp,l) = C(J,I), the collisional rates from upper
   ! level j to lower level i.
   !
   ! This is the case where the radiation field is O, hence
   ! Gamma(l',l) = Cl'l + Al'l - delta(l,l')Sum_l" (Cll" + All")
   !
   ! This Gamma is frequency and angle independent. 
   !
  ! ------------------------------------------------------------------------- !
  integer, intent(in) :: id, icell
  integer :: nact, kr, i, j, l
  type (AtomType), pointer :: atom
  integer :: c_nlte = 1
    
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom

   do kr=1,atom%Nline
    
    i = atom%lines(kr)%i; j = atom%lines(kr)%j

    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*atom%lines(kr)%Aji

   end do

   do l = 1, atom%Nlevel
    atom%Gamma(l,l) = -sum(atom%Gamma(l,:)) !sum over rows for this column
   end do

   NULLIFY(atom)
  end do !loop over atoms  

 RETURN
 END SUBROUTINE fillGamma_ZeroRadiation
 

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
   !
  ! ------------------------------------------------------------------------- !
  integer, intent(in) :: id, icell
  integer :: nact, kr, l, lp
  type (AtomType), pointer :: atom

  
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
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
    atom%Gamma(l,l) = 0d0
    atom%Gamma(l,l) = -sum(atom%Gamma(l,:)) !sum over rows for this column
   end do

!     do lp=1, atom%Nlevel
!     do l=1,atom%nlevel
!      write(*,*) lp, l, atom%Gamma(lp,l)
!     end do
!     end do
   
   NULLIFY(atom)
  end do !loop over atoms
 RETURN
 END SUBROUTINE Gamma_LTE
 
 SUBROUTINE fillGamma(id, icell, iray, n_rayons, Ieff)
  ! ------------------------------------------------------------------------- !
   ! Fill the rate matrix Gamma, whose elements are Gamma(lp,l) is the rate
   ! of transition from level lp to l.
   ! At initialisation, Gamma(lp,l) = C(J,I), the collisional rates from upper
   ! level j to lower level i.
   !
   ! It uses the fully preconditioned scheme of Rybciki and Hummer. 
   ! Switch between Hogereijde like radiation field and (M-)ALI like radiation
   ! field is done by changing Ieff.
   ! I = Ieff + Psi*eta
   ! Hogereijde: I = Idagexp(-dtau) + (1.-exp(-dtau))/chi * eta
   ! MALI      : I = Idag - Psi * etadag + eta * Psi
  ! ------------------------------------------------------------------------- !
  integer, intent(in) :: id, icell, iray, n_rayons
  integer :: nact, nc, kr, switch,i, j, Nblue, Nred, jp, krr
  type (AtomType), pointer :: atom
  double precision, dimension(NLTEspec%Nwaves), intent(inout) :: Ieff
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2
  double precision :: norm = 0d0, diag
  double precision :: c_nlte = 1d0, dummy = 0
  
!   if (iray==1) dummy = 0
  
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   if (atmos%nLTE_methode=="MALI") Ieff = &
   						NLTEspec%I(:,iray,id) - NLTEspec%Psi(:,iray,id) * atom%eta(:,id)
   !loop over transitions, b-f and b-b for these atoms
   !To do; define a transition_type with either cont or line
   do kr=1,atom%Ncont
    norm = 4d0*PI / HPLANCK / n_rayons
    
    i = atom%continua(kr)%i; j = atom%continua(kr)%j
    Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred 
    twohnu3_c2 = 2d0 * HPLANCK * CLIGHT / (NLTEspec%lambda*NM_TO_M)**3.
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%continua(kr)%Vij(:,id)*Ieff(Nblue:Nred)*cont_wlam(atom%continua(kr))) * norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%continua(kr)%Vij(:,id)*atom%continua(kr)%gij(:,id)*cont_wlam(atom%continua(kr))) * norm

    ! cross-coupling terms
     atom%Gamma(i,j) = atom%Gamma(i,j) -c_nlte*sum((atom%chi_up(i,:,id)*NLTEspec%Psi(:, iray, id)*&
     	atom%Uji_down(j,:,id))*cont_wlam(atom%continua(kr)))* norm

     ! check if i is an upper level of another transition
     do krr=1,atom%Ncont
      jp = atom%continua(krr)%j
      if (jp==i) then !i upper level of this transition
       atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((atom%chi_down(j,:,id)*nLTEspec%Psi(:, iray, id)*&
       	atom%Uji_down(i,:,id))*cont_wlam(atom%continua(kr))) * norm
      end if 
     end do
   NLTEspec%J(:,id) = NLTEspec%J(:,id) + 4d0 * PI / n_rayons * NLTEspec%I(:,iray,id)
   end do

   do kr=1,atom%Nline
    norm =  4d0 * PI / (CLIGHT * HPLANCK) / n_rayons !
    
    i = atom%lines(kr)%i; j = atom%lines(kr)%j
    Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred 
    twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji 
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%lines(kr)%Vij(:,id)*Ieff(Nblue:Nred)*atom%lines(kr)%wlam(:)) * norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:)) * norm
!     	
!     if (j==3 .and. i==2) then 
! dummy = dummy + sum(twohnu3_c2(Nblue:Nred)*atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:))*norm
!      write(*,*) atom%lines(kr)%Aji/1d7, dummy/1d7
!     end if
!      write(*,*) 'line', atom%ID, icell, id, i, j, c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
!     	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:)) * norm	
    	
    ! cross-coupling terms
     atom%Gamma(i,j) = atom%Gamma(i,j) -c_nlte*sum(atom%chi_up(i,Nblue:Nred,id)*NLTEspec%Psi(Nblue:Nred, iray, id)*&
     	atom%Uji_down(j,Nblue:Nred,id)*atom%lines(kr)%wlam(:)) * norm
     ! check if i is an upper level of another transition
     do krr=1,atom%Nline
      jp = atom%lines(krr)%j
      if (jp==i) then !i upper level of this transition
       atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum(atom%chi_down(j,Nblue:Nred,id)*nLTEspec%Psi(Nblue:Nred, iray, id)*&
       	atom%Uji_down(i,Nblue:Nred,id)*atom%lines(kr)%wlam(:)) * norm
      end if 
     end do
    !
    !NLTEspec%J(:,id) = NLTEspec%J(:,id) + 4d0 * PI / n_rayons * NLTEspec%I(:,iray,id)

   end do
   
   if (iray==n_rayons) then !otherwise we remove several times GammaDiag
   	do jp = 1, atom%Nlevel
    	atom%Gamma(jp,jp) = -sum(atom%Gamma(jp,:)) !sum over rows for this column
   	end do
   end if

   NULLIFY(atom)
  end do !loop over atoms  
 RETURN
 END SUBROUTINE fillGamma

 SUBROUTINE SEE_atom(id, icell, atom)
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
 ! --------------------------------------------------------------------!

  integer, intent(in) :: icell, id
  type(AtomType), intent(inout) :: atom
  integer :: lp, imaxpop
  double precision, dimension(atom%Nlevel, atom%Nlevel) :: Aij
  
  !we need to do that here as Gamma is updated for each rays.
!    do lp = 1, atom%Nlevel
!     write(*,*) atom%Gamma(lp,lp), sum(atom%Gamma(lp,:))
!     atom%Gamma(lp,lp) = -sum(atom%Gamma(lp,:)) !sum over rows for this column
!    end do
  

  imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
  atom%n(:,icell) = 0d0
  atom%n(imaxpop,icell) = atom%ntotal(icell)

  
  !Sum_l'_imaxpop * n_l' = N
  atom%Gamma(:,imaxpop) = 1d0 !all columns of the last row for instance
  !(G11 G21)  (n1)  (0)
  !(       ) .(  ) =( )
  !(1    1 )  (n2)  (N)
  
  !Y a peut être un transpose ici  par rapport à MCFOST, atom%Gamma.T ?

  Aij = transpose(atom%Gamma)
  CALL GaussSlv(Aij, atom%n(:,icell),atom%Nlevel)
  atom%Gamma(:,:) = -1d12

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
  integer :: nat
  
  do nat=1,atmos%NactiveAtoms
   atom => atmos%ActiveAtoms(nat)%ptr_atom
   CALL SEE_atom(id, icell, atom)
   !Ng's acceleration
   !print Ng acceleration pops.
   NULLIFY(atom)
  end do
 
 RETURN
 END SUBROUTINE updatePopulations
 
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
