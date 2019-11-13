MODULE statequil_atoms

 use atmos_type, only : atmos, nTotal_atom
 use atom_type
 use spectrum_type, only : NLTEspec
 use constant
 use opacity
 use math, only : locate, any_nan_infinity_matrix, any_nan_infinity_vector, is_nan_infinity, integrate_dx
 use utils, only : GaussSlv
 use parametres
 use accelerate
 use collision, only : CollisionRate !future deprecation
 use impact, only : Collision_Hydrogen
 use metal, only : bound_free_Xsection
 
 use mcfost_env, only : dp
 use constantes, only : tiny_dp

 IMPLICIT NONE
 
!  PROCEDURE(FillGamma_atom_hogerheijde), pointer :: FillGamma_atom => NULL()
 PROCEDURE(FillGamma_atom_mali_mu), pointer :: FillGamma_atom_mu => NULL()


 CONTAINS
 
 SUBROUTINE fill_collision_matrix(id, icell)
  integer, intent(in) :: icell, id
  type (AtomType), pointer :: atom
  integer :: nact
  
   do nact=1, atmos%NactiveAtoms
    atom => atmos%ActiveAtoms(nact)%ptr_atom
    
    CALL collision_matrix_atom(id, icell, atom)
  	 
    atom => NULL()
  enddo
  
  RETURN 
 END SUBROUTINE fill_collision_matrix
 
 SUBROUTINE collision_matrix_atom(id, icell, atom)
  integer, intent(in) :: icell, id
  type (AtomType), pointer, intent(inout) :: atom
  

    
    atom%C(:,:,id) = 0d0
	if (atom%ID=="H") then
	  atom%C(:,:,id) = Collision_Hydrogen(icell)
	else
      atom%C(:,:,id) = CollisionRate(icell, atom)
    end if
     
    if (any_nan_infinity_matrix(atom%C(:,:,id))>0) then
    	write(*,*) atom%C(:,:,id)
    	write(*,*) "Bug at compute_collision_matrix", id, icell, atom%ID, atom%n(:,icell)
    	stop
  	endif

  
  RETURN 
 END SUBROUTINE collision_matrix_atom
 
 SUBROUTINE initGamma_atom(id, atom)
  integer, intent(in) :: id
  type (AtomType), pointer, intent(inout) :: atom
  
   atom%Gamma(:,:,id) = atom%C(:,:,id)
  
  RETURN 
 END SUBROUTINE initGamma_atom
 
 SUBROUTINE initGamma(id)
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
  integer, intent(in) :: id
  type (AtomType), pointer :: atom
  
  nati = 1; natf = atmos%Nactiveatoms
  do nact = nati, natf
     atom => atmos%ActiveAtoms(nact)%ptr_atom

	 CALL initGamma_atom(id, atom)

   NULLIFY(atom)
  end do
 RETURN
 END SUBROUTINE initGamma
 
!  SUBROUTINE initGamma_atom(id, icell, atom)
!   integer, intent(in) :: icell, id
!   type (AtomType), pointer, intent(inout) :: atom
!   
!      atom%Gamma(:,:,id) = 0d0
! 	 if (atom%ID=="H") then
! 	  atom%Gamma(:,:,id) = Collision_Hydrogen(icell)
! 	  !atom%C(:,:,id) = Collision_Hydrogen(icell)
! 	 else
!       atom%Gamma(:,:,id) = CollisionRate(icell, atom)
!       !atom%C(:,:,id) = CollisionRate(icell,atom) 
!      end if
!      
!      !if (any_nan_infinity_matrix(atom%C(:,:,id))>0) then
!      if (any_nan_infinity_matrix(atom%Gamma(:,:,id))>0) then
!     	write(*,*) atom%Gamma(:,:,id)
!     	!write(*,*) atom%C(:,:,id)
!     	write(*,*) "Bug at init Gamma", id, icell, atom%ID, atom%n(:,icell)
!     	stop
!   	 endif
!   
!   RETURN 
!  END SUBROUTINE initGamma_atom
  
!  SUBROUTINE initGamma(id, icell)
!  ! ------------------------------------------------------ !
!   !for each active atom, allocate and init Gamma matrix
!   !deallocated in freeAtom()
!   ! Collision matrix has to be allocated
!   ! if already allocated, simply set Gamma to its new icell
!   ! value.
!   !
!   ! n(l)*Cl->lp = n(lp)*Clp->l
!   ! e.g.,
!   ! Cij = Cji * (nj/ni)^star, with Cji = C(ij) = 
!   ! colision rate from j to i.
!   !
!   ! (ij = j->i = (i-1)*N + j)
!   ! (ji = i->j = (j-1)*N + i)
!  ! ------------------------------------------------------ !
!   integer :: nact, Nlevel, lp, l, ij, ji, nati, natf
!   integer, intent(in) :: icell, id
!   type (AtomType), pointer :: atom
!   real(kind=dp) :: c_lte = 1d0
!   
! !   nati = (1. * (id-1)) / NLTEspec%NPROC * atmos%Nactiveatoms + 1
! !   natf = (1. * id) / NLTEspec%NPROC * atmos%Nactiveatoms
!   nati = 1; natf = atmos%Nactiveatoms
!   do nact = nati, natf
!    atom => atmos%ActiveAtoms(nact)%ptr_atom
!    
! !    write(*,*) "icell in initG", icell, "id=", id!, nact, atom%ID
!    Nlevel = atom%Nlevel
! !   open(unit=12, file="Cji_Cij_H4x4.dat",status="unknown")
! 
! !    do lp=1,Nlevel
! !     do l=1,Nlevel
! !       ij = (l-1)*nlevel + lp!lp->l
! !       ji = (lp-1)*nlevel + l!l->lp
! !       !          col/row
! !       atom%Gamma(lp,l,id) =  c_lte*atom%Ckij(icell,ij)!Gamma_llp = Gamma(lp, l) rate from lp to l
! !       !!                                                      Cul                  Clu=Cul*nu/nl with Cul=Gul here
! ! !      write(12,'(6E)') real(lp), real(l), real(ij), real(ji), atom%Ckij(icell,ij), atom%Gamma(lp,l) * atom%nstar(lp, icell)/atom%nstar(l,icell)
! !     end do
! !    end do
! 	 CALL initGamma_atom(id, icell, atom)
! ! 	 if (atom%ID=="H") then
! ! 	  atom%Gamma(:,:,id) = Collision_Hydrogen(icell)
! ! 	 else
! !       atom%Gamma(:,:,id) = CollisionRate(icell, atom) 
! !      end if
!      !write(*,*) id, icell, 'max,min Cul for atom', atom%ID, maxval(atom%Gamma(:,:,id)), minval(atom%Gamma(:,:,id))
! ! close(12)
!    NULLIFY(atom)
!   end do
!  RETURN
!  END SUBROUTINE initGamma

!  SUBROUTINE fillGamma(id, icell, n_rayons, switch_to_lte)
!   ! ------------------------------------------------------------------------- !
!    ! Fill the rate matrix Gamma, whose elements are Gamma(l',l) is the rate
!    ! of transition from level l' to l.
!    ! At initialisation, Gamma(l',l) = C(J,I), the collisional rates from upper
!    ! level j to lower level i.
!    !
!    !
!    ! The Sum_l" represents a sommation over each column for the lth row.
!    !
!   ! ------------------------------------------------------------------------- !
! 
!   integer, intent(in) :: id, n_rayons, icell !here for debug not used
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: nact
!   type (AtomType), pointer :: atom
!   
! 
!   do nact=1, atmos%NactiveAtoms
!    atom=>atmos%ActiveAtoms(nact)%ptr_atom
!    CALL FillGamma_atom(id, icell, atom, n_rayons, switch_to_lte)
!    atom=>NULL()
!   enddo
! 
!  RETURN
!  END SUBROUTINE fillGamma
!   
!  !check xcc
!  SUBROUTINE FillGamma_atom_mali(id, icell, atom, n_rayons, switch_to_lte)
! 
!   integer, intent(in) :: icell, id, n_rayons 
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: kr, kc, i, j, Nblue, Nred, l, la, lap, iray, lp, ip, jp, iray2
!   type (AtomicLine) :: line, other_line
!   type (AtomicContinuum) :: cont, other_cont
!   real(kind=dp), dimension(:), allocatable :: weight, gijk
!   real(kind=dp) :: factor, gij, twohnu3_c2, Vij, Jnu, Ieff
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) then
!     !Remove therm in delta(l,l')
!      do l = 1, atom%Nlevel
!       atom%Gamma(l,l,id) = 0d0
!       write(*,*) l, -sum(atom%Gamma(l,:,id))
!       atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!      enddo
!       write(*,*) atom%Gamma(:,:,id)
!      stop
!      RETURN !Gamma initialized to Cul
!    endif
!   end if 
!   
!   tr_loop : do kr=1, atom%Ntr
!    kc = atom%at(kr)%ik
!    
!    SELECT CASE (atom%at(kr)%trtype)
!    
!     CASE ("ATOMIC_LINE")
!      line = atom%lines(kc)
!      i = line%i; j = line%j
!      Nblue = line%Nblue; Nred = line%Nred
!     
!      allocate(weight(line%Nlambda))
! 
!      twohnu3_c2 = line%Aji / line%Bji 
! 
! 
!     !the n_rayons simplify
!      weight(:) =  line_wlam(line) * line%wphi(id)
!     
!      gij = line%Bji/line%Bij
! 
!      atom%Gamma(j,i,id) = line%Aji
!      do iray=1, n_rayons
!      
!       do l=1,line%Nlambda
!       	la = Nblue + l -1
! 
!           Ieff = NLTEspec%I(la,iray,id) - NLTEspec%Psi(la,iray,id) * atom%eta(la,iray,id)
! 
!           atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + line%Bij*line%phi(l, iray,id)*weight(l)*Ieff    
!           atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%phi(l, iray, id)*weight(l)*line%Bji*Ieff
!        
!          atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - NLTEspec%Psi(la,iray,id)*&
!          		atom%Uji_down(j,la,iray,id)*atom%chi_up(i,la,iray,id) !j>l
!          		
! !          write(*,*) "j->i", line%lambda0, j, i, NLTEspec%lambda(la)
! !          write(*,*) "Gji=", line%phi(l, iray, id)*weight(l)*line%Bji*Ieff, "Xc(j->i)=",NLTEspec%Psi(la,iray,id)*&
! !           	atom%Uji_down(j,la,iray,id)*atom%chi_up(i,la,iray,id) / norm
!          		
!          inner_tr_loop : do lp=1, atom%Ntr
!           SELECT CASE (atom%at(lp)%trtype)
!            CASE ("ATOMIC_LINE")
!             other_line = atom%lines(atom%at(lp)%ik)
!             ip = other_line%i
!             jp = other_line%j
!            CASE ("ATOMIC_CONTINUUM")
!             other_cont = atom%continua(atom%at(lp)%ik)
!             ip = other_cont%i
!             jp = other_cont%j
!           END SELECT
!           !Only if the lower level of j->i is an upper level of another transition
!           if (jp==i) then
! !             lap = overlap(l) !index on the global grid where they overlap
! !             if (lap <= 0) cycle inner_tr_loop
! 			lap = la
! ! 			write(*,*) "....>", atom%at(lp)%trtype,  "i->lower",i, ip
!              atom%Gamma(i,j,id) = atom%Gamma(i,j,id) +  NLTEspec%Psi(lap,iray,id)*&
!              	atom%Uji_down(i,lap,iray,id)*atom%chi_down(j, lap, iray, id) !j<l
! !             	
! !         write(*,*) "Xc(i->j)=", NLTEspec%Psi(lap,iray,id)*&
! !             	atom%Uji_down(i,lap,iray,id)*atom%chi_down(j, lap, iray, id) / norm2  
!           endif
!          enddo inner_tr_loop
! 
!       enddo
!      enddo
! 
! !      if (atmos%include_xcoupling) deallocate(overlap)
!      deallocate(weight)
!         
!     CASE ("ATOMIC_CONTINUUM")
!      cont = atom%continua(kc)
!      i = cont%i; j = cont%j
!      Nblue = cont%Nblue; Nred = cont%Nred
!     
!      allocate(weight(cont%Nlambda), gijk(cont%Nlambda))
! 
!      gijk(:) = 0d0 !avoids to check at each lambda the condition on n(j)
!      !nstar(i)/nstar(j)*exp(-hnu/kT)
!      if (atom%nstar(j, icell) >0) &
!      	gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell)* &
!     		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
!     
!      weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
!      
!      ! alpha_nu * 4pi / h * dnu / nu = domega dnu/hnu alpha_nu
! 
!      !explicit do loop
!      do l=1,cont%Nlambda
!       twohnu3_c2 = twohc / NLTEspec%lambda(Nblue-1+l)**(3d0)
!       Jnu = 0d0
!       la = Nblue + l -1
!       do iray=1, n_rayons
!       
!            Jnu = Jnu + NLTEspec%I(la,iray,id) - NLTEspec%Psi(la,iray,id) * atom%eta(la,iray,id)        
! 
!       enddo
!       Jnu = Jnu * cont%wmu(id)
! 
!       atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu*weight(l)
!       atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + (Jnu+twohnu3_c2)*gijk(l)*weight(l)
!       
! 
!        do iray=1,n_rayons
!         atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - atom%Uji_down(j,la,iray,id)*atom%chi_up(i,la,iray,id)*&
!         		NLTEspec%Psi(la,iray,id)
!        enddo
!        
!          inner_tr_loop2 : do lp=1, atom%Ntr
!           SELECT CASE (atom%at(lp)%trtype)
!            CASE ("ATOMIC_LINE")
!             other_line = atom%lines(atom%at(lp)%ik)
!             ip = other_line%i
!             jp = other_line%j
!            CASE ("ATOMIC_CONTINUUM")
!             other_cont = atom%continua(atom%at(lp)%ik)
!             ip = other_cont%i
!             jp = other_cont%j
!           END SELECT
!           if (jp==i) then
!           	do iray=1,n_rayons
! 			 lap = la
!              atom%Gamma(i,j,id) = atom%Gamma(i,j,id) +  NLTEspec%Psi(lap,iray,id)*&
!             	atom%Uji_down(i,lap,iray,id)*atom%chi_down(j, lap, iray, id)
!             enddo
!           endif
!           
!          enddo inner_tr_loop2
! 
!       
!      enddo
! 
! !      if (atmos%include_xcoupling) deallocate(overlap)
!      deallocate(weight, gijk)
!  
!     CASE DEFAULT
!     
!      CALL Error("Unkown transition type", atom%at(kr)%trtype)
!      
!    END SELECT
!   
!   end do tr_loop
!   
!   !Remove therm in delta(l,l')
!   !But cross-coupling are in ?? and should not
!   do l = 1, atom%Nlevel
!     atom%Gamma(l,l,id) = 0d0
!     atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!   end do
!  
!  RETURN
!  END SUBROUTINE FillGamma_atom_mali
!  
!  !--> some discrepancies with Hogereheijde not mu, should be resolved the problem was in I*exp(-tau) + Psi*eta, a minus was here before the +
!  SUBROUTINE FillGamma_atom_hogerheijde(id, icell, atom, n_rayons, switch_to_lte)
! 
!   integer, intent(in) :: icell, id, n_rayons 
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: kr, kc, i, j, Nblue, Nred, l, la, lap, iray, lp, ip, jp, iray2
!   type (AtomicLine) :: line, other_line
!   type (AtomicContinuum) :: cont, other_cont
!   real(kind=dp), dimension(:), allocatable :: weight, gijk
!   real(kind=dp) :: factor, gij, twohnu3_c2, Vij, Jnu, Ieff
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) then
!     !Remove therm in delta(l,l')
!      do l = 1, atom%Nlevel
!       atom%Gamma(l,l,id) = 0d0
!       write(*,*) l, -sum(atom%Gamma(l,:,id))
!       atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!      enddo
!       write(*,*) atom%Gamma(:,:,id)
!      stop
!      RETURN !Gamma initialized to Cul
!    endif
!   end if 
!   
!   tr_loop : do kr=1, atom%Ntr
!    kc = atom%at(kr)%ik
!    
!    SELECT CASE (atom%at(kr)%trtype)
!    
!     CASE ("ATOMIC_LINE")
!      line = atom%lines(kc)
!      i = line%i; j = line%j
!      Nblue = line%Nblue; Nred = line%Nred
!     
!      allocate(weight(line%Nlambda))
! 
!      twohnu3_c2 = line%Aji / line%Bji 
! 
! 
!     !the n_rayons simplify
!      weight(:) =  line_wlam(line) * line%wphi(id)
!     
!      gij = line%Bji/line%Bij
! 
!      atom%Gamma(j,i,id) = line%Aji !init
!      do iray=1, n_rayons
!      
!       do l=1,line%Nlambda
!       	la = Nblue + l -1
! 
!        Ieff = NLTEspec%I(la,iray,id)*dexp(-NLTEspec%tau(la,iray,id)) + &
!              NLTEspec%Psi(la,iray,id) * atom%eta(la,iray,id)
! 
! 
!        atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + line%Bij*line%phi(l, iray,id)*weight(l)*Ieff    
!        atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%phi(l, iray, id)*weight(l)*line%Bji*Ieff
! 
! 
!       enddo
!      enddo
! 
!      deallocate(weight)
!         
!     CASE ("ATOMIC_CONTINUUM")
!      cont = atom%continua(kc)
!      i = cont%i; j = cont%j
!      Nblue = cont%Nblue; Nred = cont%Nred
!     
!      allocate(weight(cont%Nlambda), gijk(cont%Nlambda))
! 
!      gijk(:) = 0d0
!      !nstar(i)/nstar(j)*exp(-hnu/kT)
!      if (atom%nstar(j, icell) >0) &
!      	gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell)* &
!     		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
!     
!      weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
!      
!      ! alpha_nu * 4pi / h * dnu / nu = domega dnu/hnu alpha_nu
! 
!      !explicit do loop
!      do l=1,cont%Nlambda
!       twohnu3_c2 = twohc / NLTEspec%lambda(Nblue-1+l)**(3d0)
!       Jnu = 0d0
!       la = Nblue + l -1
!       do iray=1, n_rayons
! 
!         Jnu = Jnu + NLTEspec%I(Nblue+l-1,iray,id)*dexp(-NLTEspec%tau(Nblue+l-1,iray,id)) + &
!              NLTEspec%Psi(Nblue+l-1,iray,id) * atom%eta(Nblue+l-1,iray,id)
! 
!       enddo
!       Jnu = Jnu * cont%wmu(id)
! 
!       atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu*weight(l)
!       atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + (Jnu+twohnu3_c2)*gijk(l)*weight(l)
!       
!      enddo
! 
! !      if (atmos%include_xcoupling) deallocate(overlap)
!      deallocate(weight, gijk)
!  
!     CASE DEFAULT
!     
!      CALL Error("Unkown transition type", atom%at(kr)%trtype)
!      
!    END SELECT
!   
!   end do tr_loop
! 
!   do l = 1, atom%Nlevel
!     atom%Gamma(l,l,id) = 0d0
!     atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!   end do
!  
!  RETURN
!  END SUBROUTINE FillGamma_atom_hogerheijde
 
 SUBROUTINE fillGamma_mu(id, icell, iray, n_rayons, switch_to_lte)
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

  integer, intent(in) :: id, n_rayons, icell, iray !here for debug not used
  logical, optional, intent(in) :: switch_to_lte
  integer :: nact
  type (AtomType), pointer :: atom
  

  do nact=1, atmos%NactiveAtoms
   atom=>atmos%ActiveAtoms(nact)%ptr_atom
   CALL FillGamma_atom_mu(id, icell, iray, atom, n_rayons, switch_to_lte)
   atom=>NULL()
  enddo

 RETURN
 END SUBROUTINE fillGamma_mu
  
 SUBROUTINE FillGamma_atom_mali_mu(id, icell, iray, atom, n_rayons, switch_to_lte)

  integer, intent(in) :: icell, id, n_rayons, iray
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: kr, kc, i, j, l, Nblue, Nred, ip, jp, Nl2, Nb2, Nr2, lp
  real(kind=dp), dimension(:), allocatable :: integ, lambda
  real(kind=dp) :: factor, Vij, Jnu, Ieff, Jeff
  !to test
  real(kind=dp) :: GammaRad(atom%Nlevel,atom%Nlevel)

  if (present(switch_to_lte)) then
   if (switch_to_lte) then
    !Remove therm in delta(l,l')
     do l = 1, atom%Nlevel
      atom%Gamma(l,l,id) = 0d0
      write(*,*) l, -sum(atom%Gamma(l,:,id))
      atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
     enddo
     RETURN !Gamma initialized to Cul
   endif
  end if 
  
  tr_loop : do kr=1, atom%Ntr
   kc = atom%at(kr)%ik
   
   SELECT CASE (atom%at(kr)%trtype)
   
    CASE ("ATOMIC_LINE")
     i = atom%lines(kc)%i; j = atom%lines(kc)%j
     Nblue = atom%lines(kc)%Nblue; Nred = atom%lines(kc)%Nred
    
     allocate(integ(atom%lines(kc)%Nlambda))!, lambda(line%Nlambda))
     !lambda = (NLTEspec%lambda(Nblue:Nred)-line%lambda0)*Clight/line%lambda0

    
     integ(:) = (NLTEspec%I(Nblue:Nred,iray,id) - NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id)) * atom%lines(kc)%phi_ray(:, iray,id)

     !Jeff = Integrate_x(line%Nlambda, lambda, integ) / n_rayons
     Jeff = sum(integ*atom%lines(kc)%w_lam(:))/n_rayons

     !if (iray == 1) atom%Gamma(j,i,id) = line%Aji

     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + Jeff*atom%lines(kc)%Bji + atom%lines(kc)%Aji/n_rayons
     atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jeff*atom%lines(kc)%Bij
     
     integ(:) = 0d0
     !include weights
     !integ = NLTEspec%Psi(Nblue:Nred,iray,id)*&
     !    		atom%Uji_down(j,Nblue:Nred,iray,id)*atom%chi_up(i,Nblue:Nred,iray,id)
     !atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - sum(integ)!Integrate_x(line%Nlambda, lambda, integ) / n_rayons !j>l
     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - &
     	sum(NLTEspec%Psi(:,iray,id)*atom%Uji_down(j,:,iray,id)*atom%chi_up(i,:,iray,id)) / n_rayons
         		
         inner_tr_loop : do lp=1, atom%Ntr
          SELECT CASE (atom%at(lp)%trtype)
           CASE ("ATOMIC_LINE")
            ip = atom%lines(atom%at(lp)%ik)%i
            jp = atom%lines(atom%at(lp)%ik)%j
            Nb2 = atom%lines(atom%at(lp)%ik)%Nblue
            Nr2 = atom%lines(atom%at(lp)%ik)%Nred
            Nl2 = atom%lines(atom%at(lp)%ik)%Nlambda
            !deallocate(lambda)
			!allocate(lambda(Nl2)); lambda = (NLTEspec%lambda(Nb2:Nr2)-other_line%lambda0)*CLIGHT/other_line%lambda0
           CASE ("ATOMIC_CONTINUUM")
            ip = atom%continua(atom%at(lp)%ik)%i
            jp = atom%continua(atom%at(lp)%ik)%j
            Nb2 = atom%continua(atom%at(lp)%ik)%Nblue
            Nr2 = atom%continua(atom%at(lp)%ik)%Nred
            Nl2 = atom%continua(atom%at(lp)%ik)%Nlambda
            !deallocate(lambda)
			!allocate(lambda(Nl2)); lambda = NLTEspec%lambda(Nb2:Nr2)

          END SELECT
          !Only if the lower level of j->i is an upper level of another transition
          if (jp==i) then              
          	 integ = 0d0
          	 !integ = NLTEspec%Psi(Nblue:Nred,iray,id)*&
             !	atom%Uji_down(i,Nblue:Nred,iray,id)*atom%chi_down(j,Nblue:Nred, iray, id)
             !atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + sum(integ)!Integrate_x(Nl2, lambda, integ) / n_rayons  !j<l
             atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + &
             	sum(NLTEspec%Psi(:,iray,id)*atom%Uji_down(i,:,iray,id)*atom%chi_down(j,:, iray, id))/n_rayons
          endif
         enddo inner_tr_loop

     deallocate(integ)!, lambda)
        
    CASE ("ATOMIC_CONTINUUM")
     i = atom%continua(kc)%i; j = atom%continua(kc)%j
     Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
    
     allocate(integ(atom%continua(kc)%Nlambda))!, lambda(cont%Nlambda))
     !lambda = cont_wlam(cont) !NLTEspec%lambda(Nblue:Nred)

    
     
      
      integ(:) = (NLTEspec%I(Nblue:Nred,iray,id)  - NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id))*atom%continua(kc)%alpha(:)
      
      Jnu = sum(atom%continua(kc)%w_lam(:)*integ)/n_rayons!Integrate_x(cont%Nlambda, lambda, integ) / n_rayons !Jbarcont
      
      integ(:) = ( (NLTEspec%I(Nblue:Nred,iray,id) - NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id)) + &
       atom%continua(kc)%twohnu3_c2(:))*atom%continua(kc)%gij(:,icell)*atom%continua(kc)%alpha(:)
       
      Jeff = sum(atom%continua(kc)%w_lam(:)*integ)/n_rayons!Integrate_x(cont%Nlambda, lambda, integ) / n_rayons

      atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu * fourPI_h
      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + Jeff * fourPI_h
      
      !integ(:) = atom%Uji_down(j,Nblue:Nred,iray,id)*atom%chi_up(i,Nblue:Nred,iray,id)*&
      !  		NLTEspec%Psi(Nblue:Nred,iray,id)
      !atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - sum(integ)!Integrate_x(cont%Nlambda, lambda, integ) / n_rayons
      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - &
      	sum(atom%Uji_down(j,:,iray,id)*atom%chi_up(i,:,iray,id)*NLTEspec%Psi(:,iray,id))/n_rayons
       
         inner_tr_loop2 : do lp=1, atom%Ntr
          SELECT CASE (atom%at(lp)%trtype)
           CASE ("ATOMIC_LINE")
            ip = atom%lines(atom%at(lp)%ik)%i
            jp = atom%lines(atom%at(lp)%ik)%j
            Nb2 = atom%lines(atom%at(lp)%ik)%Nblue
            Nr2 = atom%lines(atom%at(lp)%ik)%Nred
            Nl2 = atom%lines(atom%at(lp)%ik)%Nlambda
             !deallocate(lambda)
             !allocate(lambda(Nl2)); lambda = (NLTEspec%lambda(Nb2:Nr2)-other_line%lambda0)*CLIGHT/other_line%lambda0
           CASE ("ATOMIC_CONTINUUM")
            ip = atom%continua(atom%at(lp)%ik)%i
            jp = atom%continua(atom%at(lp)%ik)%j
            Nb2 = atom%continua(atom%at(lp)%ik)%Nblue
            Nr2 = atom%continua(atom%at(lp)%ik)%Nred
            Nl2 = atom%continua(atom%at(lp)%ik)%Nlambda
            ! deallocate(lambda)
            ! allocate(lambda(Nl2)); lambda = NLTEspec%lambda(Nb2:Nr2)
          END SELECT
          if (jp==i) then
          	 integ(:) = 0d0
			 !integ(:) = NLTEspec%Psi(Nblue:Nred,iray,id)*&
             !	atom%Uji_down(i,Nblue:Nred,iray,id)*atom%chi_down(j,Nblue:Nred, iray, id)
             !atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + sum(integ)!Integrate_x(Nl2, lambda, integ) / n_rayons
             atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + &
             	sum(NLTEspec%Psi(:,iray,id)*atom%Uji_down(i,:,iray,id)*atom%chi_down(j,:, iray, id))/n_rayons
          endif
          
         enddo inner_tr_loop2

     deallocate(integ)!, lambda)
 
    CASE DEFAULT
    
     CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
   END SELECT
  
  end do tr_loop
  
  !Remove therm in delta(l,l')
  !But cross-coupling are in ?? and should not
  do l = 1, atom%Nlevel
    atom%Gamma(l,l,id) = 0d0
    atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
  end do
  
  !invoked here to avoid the term delta
  !CALL xcc_mali_mu(id, iray, atom)
 
 RETURN
 END SUBROUTINE FillGamma_atom_mali_mu
 
 
!  SUBROUTINE xcc_mali_mu(id, iray, atom)
! 
!   integer, intent(in) :: id, iray
!   type (AtomType), intent(inout), pointer :: atom
!   integer :: kr, kc, i, j, l, Nblue, Nred, ip, jp, Nl2, Nb2, Nr2, lp
!   real(kind=dp), dimension(:), allocatable :: integ, lambda
! 
!   
!   tr_loop : do kr=1, atom%Ntr
!    kc = atom%at(kr)%ik
!    
!    SELECT CASE (atom%at(kr)%trtype)
!    
!     CASE ("ATOMIC_LINE")
!      line = atom%lines(kc)
!      i = line%i; j = line%j
!      Nblue = line%Nblue; Nred = line%Nred
!     
!      allocate(integ(line%Nlambda), lambda(line%Nlambda))
!      !lambda = (NLTEspec%lambda(Nblue:Nred)-line%lambda0)*Clight/line%lambda0
!      lambda = line_wlam(line)
! 
!  
!      integ(:) = 0d0
!      !include weights
!      !integ = NLTEspec%Psi(Nblue:Nred,iray,id)*&
!      !    		atom%Uji_down(j,Nblue:Nred,iray,id)*atom%chi_up(i,Nblue:Nred,iray,id)
!      !atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - sum(integ)!Integrate_x(line%Nlambda, lambda, integ) / n_rayons !j>l
!      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - &
!      	sum(NLTEspec%Psi(:,iray,id)*atom%Uji_down(j,:,iray,id)*atom%chi_up(i,:,iray,id))
!          		
!          inner_tr_loop : do lp=1, atom%Ntr
!           SELECT CASE (atom%at(lp)%trtype)
!            CASE ("ATOMIC_LINE")
!             other_line = atom%lines(atom%at(lp)%ik)
!             ip = other_line%i
!             jp = other_line%j
!             Nb2 = other_line%Nblue
!             Nr2 = other_line%Nred
!             Nl2 = other_line%Nlambda
!             !deallocate(lambda)
! 			!allocate(lambda(Nl2)); lambda = (NLTEspec%lambda(Nb2:Nr2)-other_line%lambda0)*CLIGHT/other_line%lambda0
!            CASE ("ATOMIC_CONTINUUM")
!             other_cont = atom%continua(atom%at(lp)%ik)
!             ip = other_cont%i
!             jp = other_cont%j
!             Nb2 = other_cont%Nblue
!             Nr2 = other_cont%Nred
!             Nl2 = other_cont%Nlambda
!             !deallocate(lambda)
! 			!allocate(lambda(Nl2)); lambda = NLTEspec%lambda(Nb2:Nr2)
! 
!           END SELECT
!           !Only if the lower level of j->i is an upper level of another transition
!           if (jp==i) then              
!           	 integ = 0d0
!           	 !integ = NLTEspec%Psi(Nblue:Nred,iray,id)*&
!              !	atom%Uji_down(i,Nblue:Nred,iray,id)*atom%chi_down(j,Nblue:Nred, iray, id)
!              !atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + sum(integ)!Integrate_x(Nl2, lambda, integ) / n_rayons  !j<l
!              atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + &
!              	sum(NLTEspec%Psi(:,iray,id)*atom%Uji_down(i,:,iray,id)*atom%chi_down(j,:, iray, id))
!           endif
!          enddo inner_tr_loop
! 
!      deallocate(integ, lambda)
!         
!     CASE ("ATOMIC_CONTINUUM")
!      cont = atom%continua(kc)
!      i = cont%i; j = cont%j
!      Nblue = cont%Nblue; Nred = cont%Nred
!     
!      allocate(integ(cont%Nlambda), lambda(cont%Nlambda))
!      lambda = cont_wlam(cont) !NLTEspec%lambda(Nblue:Nred)
!       
!       !integ(:) = atom%Uji_down(j,Nblue:Nred,iray,id)*atom%chi_up(i,Nblue:Nred,iray,id)*&
!       !  		NLTEspec%Psi(Nblue:Nred,iray,id)
!       !atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - sum(integ)!Integrate_x(cont%Nlambda, lambda, integ) / n_rayons
!       atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - &
!       	sum(atom%Uji_down(j,:,iray,id)*atom%chi_up(i,:,iray,id)*NLTEspec%Psi(:,iray,id))
!        
!          inner_tr_loop2 : do lp=1, atom%Ntr
!           SELECT CASE (atom%at(lp)%trtype)
!            CASE ("ATOMIC_LINE")
!             other_line = atom%lines(atom%at(lp)%ik)
!             ip = other_line%i
!             jp = other_line%j
!             Nb2 = other_line%Nblue
!             Nr2 = other_line%Nred
!             Nl2 = other_line%Nlambda
!              !deallocate(lambda)
!              !allocate(lambda(Nl2)); lambda = (NLTEspec%lambda(Nb2:Nr2)-other_line%lambda0)*CLIGHT/other_line%lambda0
!            CASE ("ATOMIC_CONTINUUM")
!             other_cont = atom%continua(atom%at(lp)%ik)
!             ip = other_cont%i
!             jp = other_cont%j
!             Nb2 = other_cont%Nblue
!             Nr2 = other_cont%Nred
!             Nl2 = other_cont%Nlambda
!             ! deallocate(lambda)
!             ! allocate(lambda(Nl2)); lambda = NLTEspec%lambda(Nb2:Nr2)
!           END SELECT
!           if (jp==i) then
!           	 integ(:) = 0d0
! 			 !integ(:) = NLTEspec%Psi(Nblue:Nred,iray,id)*&
!              !	atom%Uji_down(i,Nblue:Nred,iray,id)*atom%chi_down(j,Nblue:Nred, iray, id)
!              !atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + sum(integ)!Integrate_x(Nl2, lambda, integ) / n_rayons
!              atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + &
!              	sum(NLTEspec%Psi(:,iray,id)*atom%Uji_down(i,:,iray,id)*atom%chi_down(j,:, iray, id))
!           endif
!           
!          enddo inner_tr_loop2
! 
!      deallocate(integ, lambda)
!  
!     CASE DEFAULT
!     
!      CALL Error("Unkown transition type", atom%at(kr)%trtype)
!      
!    END SELECT
!   
!   end do tr_loop
! 
!  
!  RETURN
!  END SUBROUTINE xcc_mali_mu
 

!  SUBROUTINE FillGamma_atom_hogerheijde_mu(id, icell, iray, atom, n_rayons, switch_to_lte)
! 
!   integer, intent(in) :: icell, id, n_rayons, iray
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: kr, kc, i, j, l, Nblue, Nred, ip, jp, Nl2, Nb2, Nr2, lp
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: integ
!   real(kind=dp) :: factor, Vij, Jnu, Ieff, Jeff
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) then
!     !Remove therm in delta(l,l')
!      do l = 1, atom%Nlevel
!       atom%Gamma(l,l,id) = 0d0
!       write(*,*) l, -sum(atom%Gamma(l,:,id))
!       atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!      enddo
!      RETURN !Gamma initialized to Cul
!    endif
!   end if 
!   
!   tr_loop : do kr=1, atom%Ntr
!    kc = atom%at(kr)%ik
!    
!    SELECT CASE (atom%at(kr)%trtype)
!    
!     CASE ("ATOMIC_LINE")
!      i = atom%lines(kc)%i; j = atom%lines(kc)%j
!      Nblue = atom%lines(kc)%Nblue; Nred = atom%lines(kc)%Nred
! 
!      integ(Nblue:Nred) = (NLTEspec%I(Nblue:Nred,iray,id)*dexp(-NLTEspec%tau(Nblue:Nred,iray,id)) + NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id)) * &
!        atom%lines(kc)%phi_ray(:, iray,id)
!        
!      Jeff = sum(integ(Nblue:Nred)*atom%lines(kc)%w_lam(:))/n_rayons
! 
!      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + Jeff*atom%lines(kc)%Bji + atom%lines(kc)%Aji/n_rayons
!      atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jeff*atom%lines(kc)%Bij
!     
!         
!     CASE ("ATOMIC_CONTINUUM")
!      i = atom%continua(kc)%i; j = atom%continua(kc)%j
!      Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
! 
!      integ(Nblue:Nred) = (NLTEspec%I(Nblue:Nred,iray,id)*dexp(-NLTEspec%tau(Nblue:Nred,iray,id)) + NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id))*atom%continua(kc)%alpha(:)
!       
!       Jnu = sum(atom%continua(kc)%w_lam(:)*integ(Nblue:Nred))/n_rayons
!       
!       
!       integ(Nblue:Nred) = ( (NLTEspec%I(Nblue:Nred,iray,id)*dexp(-NLTEspec%tau(Nblue:Nred,iray,id)) + NLTEspec%Psi(Nblue:Nred,iray,id) * atom%eta(Nblue:Nred,iray,id)) + &
!        atom%continua(kc)%twohnu3_c2(:))*atom%continua(kc)%gij(:,icell)*atom%continua(kc)%alpha(:)
!       Jeff = sum(atom%continua(kc)%w_lam(:)*integ)/n_rayons
!       
!       atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu * fourPI_h
!       atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + Jeff * fourPI_h
!        
!     CASE DEFAULT
!     
!      CALL Error("Unkown transition type", atom%at(kr)%trtype)
!      
!    END SELECT
!   
!   end do tr_loop
!   
!   !Remove therm in delta(l,l')
!   !But cross-coupling are in ?? and should not
!   do l = 1, atom%Nlevel
!     atom%Gamma(l,l,id) = 0d0
!     atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
!   end do
!  
!  RETURN
!  END SUBROUTINE FillGamma_atom_hogerheijde_mu 

 SUBROUTINE FillGamma_atom_hogerheijde_mu(id, icell, iray, atom, n_rayons, switch_to_lte)

  integer, intent(in) :: icell, id, n_rayons, iray
  type (AtomType), intent(inout), pointer :: atom
  logical, optional, intent(in) :: switch_to_lte
  integer :: kr, kc, i, j, l, Nblue, Nred, Nl
  real(kind=dp) :: nu, Ieff, Jeff, Jeffp

  if (present(switch_to_lte)) then
   if (switch_to_lte) then
    !Remove therm in delta(l,l')
     do l = 1, atom%Nlevel
      atom%Gamma(l,l,id) = 0d0
      !write(*,*) l, -sum(atom%Gamma(l,:,id))
      atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
     enddo
     RETURN !Gamma initialized to Cul
   endif
  end if 
  
  tr_loop : do kr=1, atom%Ntr
   kc = atom%at(kr)%ik
   
   SELECT CASE (atom%at(kr)%trtype)
   
    CASE ("ATOMIC_LINE")
     i = atom%lines(kc)%i; j = atom%lines(kc)%j
     Nblue = atom%lines(kc)%Nblue; Nred = atom%lines(kc)%Nred
     Nl = atom%lines(kc)%Nlambda

       
     Jeff = 0.0_dp
     do l=1, Nl
      Ieff = NLTEspec%I(Nblue+l-1,iray,id)*dexp(-NLTEspec%tau(Nblue+l-1,iray,id)) + &
      		 NLTEspec%Psi(Nblue+l-1,iray,id) * atom%eta(Nblue+l-1,iray,id)
      		 
      Jeff = Jeff + Ieff * atom%lines(kc)%phi_ray(l,iray,id)*atom%lines(kc)%w_lam(l) / n_rayons      
     enddo

     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + Jeff*atom%lines(kc)%Bji + atom%lines(kc)%Aji/n_rayons
     atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jeff*atom%lines(kc)%Bij
    
        
    CASE ("ATOMIC_CONTINUUM")
     i = atom%continua(kc)%i; j = atom%continua(kc)%j
     Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
     Nl = atom%continua(kc)%Nlambda


     Jeff = 0.0_dp
     Jeffp = 0.0_dp
     do l=1, Nl
      Ieff = NLTEspec%I(Nblue+l-1,iray,id)*dexp(-NLTEspec%tau(Nblue+l-1,iray,id)) + &
      		 NLTEspec%Psi(Nblue+l-1,iray,id) * atom%eta(Nblue+l-1,iray,id)
      		 
      Jeff = Jeff + Ieff * atom%continua(kc)%alpha(l) * atom%continua(kc)%w_lam(l) / n_rayons
      Jeffp = Jeffp + (Ieff + atom%continua(kc)%twohnu3_c2(l)) * &
              atom%continua(kc)%gij(l,icell) * atom%continua(kc)%alpha(l) * atom%continua(kc)%w_lam(l) / n_rayons
     enddo
      
      atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jeff * fourPI_h
      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + Jeffp * fourPI_h
       
    CASE DEFAULT
    
     CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
   END SELECT
  
  end do tr_loop
  
  !Remove therm in delta(l,l')
  !But cross-coupling are in ?? and should not
  do l = 1, atom%Nlevel
    atom%Gamma(l,l,id) = 0.0_dp
    atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
  end do
 
 RETURN
 END SUBROUTINE FillGamma_atom_hogerheijde_mu 


 SUBROUTINE FillGamma_atom_zero_radiation(id, icell, atom)
 !here, I develop eq. 21 of Uitenbroek 2001, and I substitue I with I*dexp(-tau) + Psi * eta
 ! "Hogereijde-like". There is no expansion of eta in S_jS_i Uji * nj, hence no Xcoupling yet
 ! to do: Xcoupling
  integer, intent(in) :: icell, id
  type (AtomType), intent(inout), pointer :: atom
  integer :: kr, kc, i, j, Nblue, Nred, l, iray
  !no too much memory needed for that
  !type (AtomicLine) :: line
  !type (AtomicContinuum) :: cont

  tr_loop : do kr=1, atom%Ntr
   kc = atom%at(kr)%ik
   
   SELECT CASE (atom%at(kr)%trtype)
   
    CASE ("ATOMIC_LINE")
     !line = atom%lines(kc)
     i = atom%lines(kc)%i; j = atom%lines(kc)%j
     Nblue = atom%lines(kc)%Nblue; Nred = atom%lines(kc)%Nred


     atom%Gamma(j,i,id) = atom%lines(kc)%Aji !init

        
    CASE ("ATOMIC_CONTINUUM")
     !cont = atom%continua(kc)
     i = atom%continua(kc)%i; j = atom%continua(kc)%j
     Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
    
     !allocate(weight(cont%Nlambda), gijk(cont%Nlambda))
    
     !weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
     ! alpha_nu * 4pi / h * dnu / nu = domega dnu/hnu alpha_nu

     !explicit do loop
     do l=1,atom%continua(kc)%Nlambda

      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + &
       fourPI_h * atom%continua(kc)%alpha(l)*atom%continua(kc)%twohnu3_c2(l)*atom%continua(kc)%gij(l,icell)*atom%continua(kc)%w_lam(l)
         
     enddo
    
 
    CASE DEFAULT
    
     CALL Error("Unkown transition type", atom%at(kr)%trtype)
     
   END SELECT
  
  end do tr_loop
  
  do l = 1, atom%Nlevel
    atom%Gamma(l,l,id) = 0d0
    atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
  end do

  
 RETURN
 END SUBROUTINE FillGamma_atom_zero_radiation
 
 
 !This version assumes, that the term in delta(l,l') has already been removed
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
  real(kind=dp), dimension(atom%Nlevel, atom%Nlevel) :: Aij
  
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
  !Aij = atom%Gamma(:,:,id)
! write(*,*) id, icell
! write(*,*) "nafter=",atom%n(:,icell)
! write(*,*) "Aij-Gamma(id)", Aij-transpose(atom%Gamma(:,:,id))
  CALL GaussSlv(Aij, atom%n(:,icell),atom%Nlevel)
  !!atom%n(:,icell) = atom%n(:,icell) * nTotal_atom(icell, atom)
  if ((any_nan_infinity_matrix(atom%Gamma(:,:,id))>0).or.&
  				(any_nan_infinity_vector(atom%n(:,icell))>0)) then
    write(*,*) atom%Gamma(:,:,id)
    write(*,*) id, icell, atom%n(:,icell)
    write(*,*) Aij
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
  real(kind=dp) :: dM
  
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
 

!! --> OLD VERSIONS 
 !futur deprecation, FillGamma will direclty compute the full matrix for lines and continua
 !using the new atom%at informations and atom%Ntr
!  SUBROUTINE fillGamma_1(id, icell, n_rayons, switch_to_lte)
!   ! ------------------------------------------------------------------------- !
!    ! Fill the rate matrix Gamma, whose elements are Gamma(l',l) is the rate
!    ! of transition from level l' to l.
!    ! At initialisation, Gamma(l',l) = C(J,I), the collisional rates from upper
!    ! level j to lower level i.
!    !
!    !
!    ! The Sum_l" represents a sommation over each column for the lth row.
!    !
!   ! ------------------------------------------------------------------------- !
! 
!   integer, intent(in) :: id, n_rayons, icell !here for debug not used
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: nact
!   type (AtomType), pointer :: atom
!   
!   !Later I will replace by a unique function
!   do nact=1, atmos%NactiveAtoms
!    atom=>atmos%ActiveAtoms(nact)%ptr_atom
!    CALL FillGamma_bf_hjde(id, icell, atom, n_rayons, switch_to_lte)
!    CALL FillGamma_bb_hjde(id, icell, atom, n_rayons, switch_to_lte)
!    atom=>NULL()
!   enddo
! 
!  RETURN
!  END SUBROUTINE fillGamma_1
!  
!  SUBROUTINE FillGamma_bb_hjde(id, icell, atom, n_rayons, switch_to_lte)
!  !here, I develop eq. 21 of Uitenbroek 2001, and I substitue I with I*dexp(-tau) + Psi * eta
!  ! "Hogereijde-like"
!   integer, intent(in) :: icell, id, n_rayons !icell is not need, except if we introduce Xcoupling terms here
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: nact, kr, i, j, Nblue, Nred, l, iray
!   !type (AtomType), pointer :: atom
!   type (AtomicLine) :: line
!   real(kind=dp), dimension(:), allocatable :: weight
!   real(kind=dp) :: factor, gij, twohnu3_c2, Vij, norm, Ieff
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) RETURN !Gamma initialized to Cul
!   end if
! 
!    do kr=1,atom%Nline
!     line = atom%lines(kr)
!     i = line%i; j = line%j
!     Nblue = line%Nblue; Nred = line%Nred
!     
!     allocate(weight(line%Nlambda)); weight(:) = line_wlam(line)
! 
!     twohnu3_c2 = line%Aji / line%Bji 
! 
!     norm = 0d0
!     do iray=1, n_rayons
!      norm = norm + sum(line%phi(:,iray,id)*weight(:))!/n_rayons !angle and frequency integrated
!     enddo
!     !the n_rayons simplify
!     weight(:) =  weight(:) / norm! / factor
! !     write(*,*) id, j, i, maxval(weight), norm
!     
!     gij = line%Bji/line%Bij
! 
!     atom%Gamma(j,i,id) = line%Aji !init
!     do iray=1, n_rayons
!      !write(*,*) "id=",id, iray, maxval(atom%eta(Nblue:Nred,iray,id))
!      do l=1,line%Nlambda
! 
!       Ieff = NLTEspec%I(Nblue-1+l,iray,id)*dexp(-NLTEspec%tau(Nblue-1+l,iray,id)) + \
!              NLTEspec%Psi(Nblue-1+l,iray,id) * atom%eta(Nblue-1+l,iray,id)
!      !write(*,*) icell, id, iray, l, Ieff(l), line%Bij, line%phi(l, iray, id), weight(l)
!       atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + line%Bij*line%phi(l, iray,id)*Ieff*weight(l)
!       atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%phi(l, iray, id)*weight(l)*line%Bji*Ieff
!       !atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%phi(l, iray, id)*weight(l)*(line%Aji+line%Bji*Ieff(l))
!  
!       
!      enddo
!     enddo
! 
!     deallocate(weight)
!    enddo !over lines
!   
! !    if (any_nan_infinity_matrix(atom%Gamma(:,:,id))>0) then
! !     write(*,*) "Error in Gamma_bb"
! !     write(*,*) atom%Gamma(:,:,id)
! !     stop
! !    end if
!  RETURN 
!  END SUBROUTINE FillGamma_bb_hjde
!  
!  ! for that atom
!  SUBROUTINE FillGamma_bf_hjde(id, icell, atom, n_rayons, switch_to_lte)
!   integer, intent(in) :: id, n_rayons, icell
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: nact, kr, i, j, Nblue, Nred, l, iray
!   !type (AtomType), pointer :: atom
!   type (AtomicContinuum) :: cont
!   real(kind=dp) :: Jnu, twohnu3_c2
!   real(kind=dp), dimension(:), allocatable :: weight, gijk
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) RETURN !Gamma initialized to Cul
!   end if
!    
!    do kr=1,atom%Ncont
!     cont = atom%continua(kr)
!     i = cont%i; j = cont%j
!     Nblue = cont%Nblue; Nred = cont%Nred
!     
!     allocate(weight(cont%Nlambda), gijk(cont%Nlambda))
!              
!     !nstar(i)/nstar(j)*exp(-hnu/kT)
!     gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell) * &
!     		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
!     
!     weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
!     ! alpha_nu * 4pi / h * dnu / nu = domega dnu/hnu alpha_nu
! 
!     !explicit do loop
!     do l=1,cont%Nlambda
!       twohnu3_c2 = twohc / NLTEspec%lambda(Nblue-1+l)**(3d0)
!       Jnu = 0d0
!       do iray=1,n_rayons
!        Jnu = Jnu + NLTEspec%I(Nblue-1+l,iray,id)*dexp(-NLTEspec%tau(Nblue-1+l,iray,id)) + \
!              NLTEspec%Psi(Nblue-1+l,iray,id) * atom%eta(Nblue-1+l,iray,id)
!       enddo
!       Jnu = Jnu / n_rayons
!      !write(*,*) j, i, NLTEspec%lambda(Nblue+l-1), Jnu(l)
!       atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + Jnu*weight(l)
!       atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + (Jnu+twohnu3_c2)*gijk(l)*weight(l)  
!      
!     enddo
!     
!     deallocate(weight, gijk)
!    end do 
! !    if (any_nan_infinity_matrix(atom%Gamma(:,:,id))>0) then
! !     write(*,*) "Error in Gamma_bf"
! !     write(*,*) atom%Gamma(:,:,id)
! !    end if
! 
!   
!  RETURN
!  END SUBROUTINE FillGamma_bf_hjde
!  
!  SUBROUTINE FillGamma_bf_zero_radiation(id, icell, atom, n_rayons, switch_to_lte)
!   integer, intent(in) :: id, n_rayons, icell
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: nact, kr, i, j, Nblue, Nred, l
!   !type (AtomType), pointer :: atom
!   type (AtomicContinuum) :: cont
!   real(kind=dp), dimension(:), allocatable :: twohnu3_c2k, weight, gijk
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) RETURN !Gamma initialized to Cul
!   end if
!   
!   !do nact=1,atmos%NactiveAtoms !loop over each active atoms
!   
!    !atom => atmos%ActiveAtoms(nact)%ptr_atom
!    
!    do kr=1,atom%Ncont
!     cont = atom%continua(kr)
!     i = cont%i; j = cont%j
!     Nblue = cont%Nblue; Nred = cont%Nred
!     
!     allocate(twohnu3_c2k(cont%Nlambda), weight(cont%Nlambda), gijk(cont%Nlambda))
!        
!     !nstar(i)/nstar(j)*exp(-hnu/kT)
!     gijk(:) = atom%nstar(i, icell)/atom%nstar(j,icell) * &
!     		 dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
!     !2hnu3/c2 = 2hc/lambda3
!     twohnu3_c2k(:) = twohc / NLTEspec%lambda(Nblue:Nred)**(3d0)
!     
!     weight(:) = fourPI_h * cont_wlam(cont) * bound_free_Xsection(cont)
!     ! alpha_nu * 4pi / h / n_rayons * dnu / nu = domega dnu/hnu alpha_nu
! 
!     !explicit do loop
!     do l=1,cont%Nlambda
!      !!atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + 0d0
!      atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + twohnu3_c2k(l)*gijk(l)*weight(l)
!     enddo
!     
!     deallocate(twohnu3_c2k, weight, gijk)
!    end do !over cont
!    
!    !atom => NULL()
!   !enddo !over atom
!   
!  RETURN
!  END SUBROUTINE FillGamma_bf_zero_radiation
!  
!  SUBROUTINE FillGamma_bb_zero_radiation(id, icell, atom, n_rayons, switch_to_lte)
!   integer, intent(in) :: icell, id, n_rayons
!   type (AtomType), intent(inout), pointer :: atom
!   logical, optional, intent(in) :: switch_to_lte
!   integer :: nact, kr, i, j
!   !type (AtomType), pointer :: atom
!   type (AtomicLine) :: line
! 
!   if (present(switch_to_lte)) then
!    if (switch_to_lte) RETURN !Gamma initialized to Cul
!   end if
!   
!   !do nact=1,atmos%NactiveAtoms !loop over each active atoms
!   
!    !atom => atmos%ActiveAtoms(nact)%ptr_atom
! 
!    do kr=1,atom%Nline
!     line = atom%lines(kr)
!     i = line%i; j = line%j
! 
!     !!atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + 0d0
!     atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + line%Aji !integration over domega/4PI dv phi = 1
! 
!    enddo !over lines
!    !NULLIFY(atom)
!   !end do !loop over atoms   
!   
!  RETURN 
!  END SUBROUTINE FillGamma_bb_zero_radiation
! 
! 
!  SUBROUTINE SEE_atom_1(id, icell,atom)
!  ! --------------------------------------------------------------------!
!   ! For atom atom solves for the Statistical Equilibrium Equations (SEE) 
!   !
!   ! We solve for :
!   !  Sum_l' Gamma_l'l n_l' = 0 (m^-3 s^-1)
!   !
!   ! which is equivalent in matrix notation to:
!   !
!   ! GG(l,l') dot n_l' = 0 with l' is on the columns and l on the rows
!   !
!   ! In particular for a 2-level atom the two equations are:
!   !
!   ! n1*G_11 + n2*G_21 = 0
!   ! n1*G12 + n2*G22 = 0,
!   ! with one of these equations has to be replaced by
!   ! n1 + n2 = N
!   !
!   ! For numerical stability, the row with the largest populations is 
!   ! removed (i.e., it is replaced by the mass conservation law).
!   !
!   ! This matrix Gamma for atoms (or later molecules) is compared to 
!   ! molecules A matrix in mcfost:
!   ! diag(Gamma) = -diag(A); offdiag(Gamma) = -offdiag(A)
!   ! A = transpose(Gamma)
!  ! --------------------------------------------------------------------!
!   
!   integer, intent(in) :: icell, id
!   type(AtomType), intent(inout) :: atom
!   integer :: lp, imaxpop, l
!   real(kind=dp), dimension(atom%Nlevel, atom%Nlevel) :: Aij
! 
!   !-> futuru, will be done in FillGamma
!    do l = 1, atom%Nlevel
!     atom%Gamma(l,l,id) = 0d0
!     atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id)) !sum over rows for this column
! !!atom%Gamma(l,l,id) = sum(abs(atom%Gamma(l,:,id))) !sum over rows for this column
!    end do
!   !write(*,*) atom%ID, id, icell, maxval(atom%Gamma(:,:,id)), minval(atom%Gamma(:,:,id))
!   !write(*,*) id, icell, atom%n(:,icell)
! 
!   imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
!   !write(*,*) "imaxpop", imaxpop, atom%n(imaxpop,icell)
!   atom%n(:,icell) = 0d0
!   atom%n(imaxpop,icell) = nTotal_atom(icell, atom)
!   
!   !Sum_l'_imaxpop * n_l' = N
!   atom%Gamma(:,imaxpop,id) = 1d0 !all columns of the last row for instance
!   !(G11 G21)  (n1)  (0)
!   !(       ) .(  ) =( )
!   !(1    1 )  (n2)  (N)
!   
!   !Y a peut être un transpose ici  par rapport à MCFOST, atom%Gamma.T ?
! 
!   Aij = transpose(atom%Gamma(:,:,id))
!   !!Aij = atom%Gamma(:,:,id)
!   CALL GaussSlv(Aij, atom%n(:,icell),atom%Nlevel)
!   !!atom%n(:,icell) = atom%n(:,icell) * nTotal_atom(icell, atom)
!   if (any_nan_infinity_matrix(Aij)>0) then
!     write(*,*) atom%Gamma(:,:,id)
!     write(*,*) id, icell, atom%n(:,icell)
!     stop
!   end if
! 
!  RETURN
!  END SUBROUTINE SEE_atom_1
!  
!  SUBROUTINE updatePopulations_1(id, icell)
!  ! --------------------------------------------------------------------!
!   ! Performs a solution of SEE for each atom.
!   !
!   ! to do: implements Ng's acceleration iteration. 
!  ! --------------------------------------------------------------------!
! 
!   integer, intent(in) :: id, icell
!   type(AtomType), pointer :: atom
!   integer :: nat, nati, natf
!   logical :: accelerate = .false.
!   real(kind=dp) :: dM
!   
!   !nati = (1. * (id-1)) / NLTEspec%NPROC * atmos%Nactiveatoms + 1
!   !natf = (1. * id) / NLTEspec%NPROC * atmos%Nactiveatoms
!   nati = 1; natf = atmos%Nactiveatoms
!   do nat=nati,natf !loop over each active atoms
!    atom => atmos%ActiveAtoms(nat)%ptr_atom
!    CALL SEE_atom_1(id, icell, atom)
!    atom => NULL()
!   end do
!  
!  RETURN
!  END SUBROUTINE updatePopulations_1

END MODULE statequil_atoms
