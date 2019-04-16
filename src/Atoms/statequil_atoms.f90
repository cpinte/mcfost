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
  double precision, intent(in) :: x0, y0, z0, u0, v0, w0
  integer, intent(in) :: id, icell, iray
  double precision :: l_dum, l_c_dum, l_v_dum, x1, x2, x3, ds
  integer 		   :: n_c_dum
  !recompute opacity of this cell., but I need angles and pos...
  !NLTEspec%I not changed
  CALL cross_cell(x0,y0,z0, u0,v0,w0, icell, &
       						n_c_dum, x1,x2,x3, n_c_dum, l_dum, l_c_dum, l_v_dum)
  NLTEspec%AtomOpac%chi(:,id) = 0d0
  NLTEspec%AtomOpac%eta(:,id) = 0d0
  !set atom%eta to zero also
  CALL initCrossCoupling(id)
  !NOTE Zeeman opacities are not re set to zero and are accumulated
  !change that or always use FIELD_FREE
  !is equivalent to P(icell,id, iray) ?
  CALL NLTEOpacity(id, icell, x0, y0, z0, x1, x2, x3, u0, v0, w0, l_dum, .true.)
  CALL fillCrossCoupling_terms(id, icell)
  !il faut recalculer Psi dans les sous-iterations et Ieff if Hogeroijde
  ds = l_dum * AU_to_m
  CALL calc_psi_operator(id, icell, iray, ds)
  
 RETURN
 END SUBROUTINE init_local_field_atom
 
 SUBROUTINE fillGamma(id, icell, iray, n_rayons, Ieff)
 ! Should work with MALI or Hogereijde.
 ! Only Psi and Ieff change for these cases.
  integer, intent(in) :: id, icell, iray, n_rayons
  integer :: nact, nc, kr, switch,i, j, Nblue, Nred, jp, krr
  type (AtomType), pointer :: atom
  double precision, dimension(NLTEspec%Nwaves), intent(inout) :: Ieff
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2
  double precision :: norm = 0d0, hc_4PI, diag
  double precision :: c_nlte = 1d0
  
  hc_4PI = HPLANCK * CLIGHT / (4d0 * PI)

  
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
!     write(*,*) atom%ID, icell, id,atom%Gamma(i,j), atom%Gamma(j,i), c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
!     	atom%continua(kr)%Vij(:,id)*atom%continua(kr)%gij(:,id)*cont_wlam(atom%continua(kr))) * norm
    ! cross-coupling terms
     atom%Gamma(i,j) = atom%Gamma(i,j) -c_nlte*sum((atom%chi_up(i,:,id)*NLTEspec%Psi(:, iray, id)*&
     	atom%Uji_down(j,:,id))*cont_wlam(atom%continua(kr)))* norm
!      write(*,*) atom%ID, icell, id, i, j, sum((atom%chi_up(i,:,id)*NLTEspec%Psi(:, iray, id)*&
!      	atom%Uji_down(j,:,id))*cont_wlam(atom%continua(kr))) * norm

     ! check if i is an upper level of another transition
     do krr=1,atom%Ncont
      jp = atom%continua(krr)%j
      if (jp==i) then !i upper level of this transition
       atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((atom%chi_down(j,:,id)*nLTEspec%Psi(:, iray, id)*&
       	atom%Uji_down(i,:,id))*cont_wlam(atom%continua(kr))) * norm
      end if 
     end do
     write(*,*) 'cont', atom%ID, icell, id, i, j, c_nlte*sum((atom%chi_down(j,:,id)*nLTEspec%Psi(:, iray, id)*&
       	atom%Uji_down(i,:,id))*cont_wlam(atom%continua(kr))) * norm

   end do

   do kr=1,atom%Nline
    i = atom%lines(kr)%i; j = atom%lines(kr)%j
    norm =  4d0 * PI / CLIGHT / HPLANCK / n_rayons
    Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred 
    twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji 
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%lines(kr)%Vij(:,id)*Ieff(Nblue:Nred)*atom%lines(kr)%wlam(:)) * norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:)) * norm
     write(*,*) 'line', atom%ID, icell, id, i, j, c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:)) * norm	
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

   end do
   
   !remove diagonal here, delat(l',l)
   do i=1,atom%Nlevel
    atom%Gamma(i,i) = 0d0 !because i that case Gamma_l"l = Cii - Cii = 0d0
    diag = 0d0
    do j=1,atom%Nlevel
     diag = diag + atom%Gamma(j,i) 
    end do
    atom%Gamma(i,i) = -diag
   end do
stop
   NULLIFY(atom)
  end do !loop over atoms  
 RETURN
 END SUBROUTINE fillGamma
 
 SUBROUTINE initGamma(icell) !icell temporary
 ! ------------------------------------------------------ !
  !for each active atom, allocate and init Gamma matrix
  !deallocated in freeAtom()
  ! Collision matrix has to be allocated
  ! if already allocated, simply set Gamma to its new icell
  ! value.
 ! ------------------------------------------------------ !
  integer :: nact, Nlevel, i, j, ij
  integer, intent(in) :: icell
  type (AtomType), pointer :: atom
  
  do nact = 1, atmos%Nactiveatoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   Nlevel = atom%Nlevel
   ij = 0
   if (.not.allocated(atom%Gamma)) allocate(atom%Gamma(Nlevel, Nlevel))
   do i=1, Nlevel
    do j=1,Nlevel
      ij = ij + 1
      atom%Gamma(i,j) = atom%Ckij(icell,ij)!atom%C(ij) 
      !write(*,*) icell, i, j, ij, atom%Gamma(i,j), atom%Ckij(icell,ij)
    end do
   end do
   NULLIFY(atom)
  end do
 
 RETURN
 END SUBROUTINE initGamma

 SUBROUTINE initCrossCoupling(id)
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

 SUBROUTINE fillGamma_Hogereijde(id, icell, iray, n_rayons, Ieffi)
  integer, intent(in) :: id, icell, iray, n_rayons
  integer :: nact, nc, kr, switch,i, j, Nblue, Nred, jp, krr
  type (AtomType), pointer :: atom
  double precision, dimension(NLTEspec%Nwaves), intent(in) :: Ieffi
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2, Ieff
  double precision :: norm = 0d0, hc_4PI, diag
  integer :: c_nlte = 0
  
  hc_4PI = HPLANCK * CLIGHT / (4d0 * PI)
  
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   Ieff = Ieffi + NLTEspec%Psi(:,iray,id) * atom%eta(:,id)
  write(*,*) "Check Norm"
  stop
   !loop over transitions, b-f and b-b for these atoms
   !To do; define a transition_type with either cont or line
   do kr=1,atom%Ncont
    norm = sum(cont_wlam(atom%continua(kr))) * n_rayons
    i = atom%continua(kr)%i; j = atom%continua(kr)%j
    Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred 
    twohnu3_c2 = 2d0 * HPLANCK * CLIGHT / (NLTEspec%lambda*NM_TO_M)**3.
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%continua(kr)%Vij(:,id)*Ieff(Nblue:Nred)*cont_wlam(atom%continua(kr))) / norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%continua(kr)%Vij(:,id)*atom%continua(kr)%gij(:,id)*cont_wlam(atom%continua(kr))) / norm
     !write(*,*) atom%ID, icell, id, i, j, atom%Gamma(i, j), atom%Gamma(j,i)
   end do

   do kr=1,atom%Nline
    i = atom%lines(kr)%i; j = atom%lines(kr)%j
    norm = 1d0 * HC_4PI* n_rayons
    Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred 
    twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji 
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%lines(kr)%Vij(:,id)*Ieff(Nblue:Nred)*atom%lines(kr)%wlam(:)) / norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:))
   end do
   
   !remove diagonal here, delat(l',l)
   do i=1,atom%Nlevel
    atom%Gamma(i,i) = 0d0 !because i that case Gamma_l"l = Cii - Cii = 0d0
    diag = 0d0
    do j=1,atom%Nlevel
     diag = diag + atom%Gamma(j,i) 
    end do
    atom%Gamma(i,i) = -diag
   end do
   
   NULLIFY(atom)
  end do !loop over atoms  
 
 RETURN
 END SUBROUTINE fillGamma_Hogereijde
 
 SUBROUTINE fillGamma_MALI(id,icell,iray,n_rayons, switch_to_ALI)
 !------------------------------------------------- !
  ! for all atoms, performs wavelength
  ! and angle integration of rate matrix
  ! The angle integration is done by nray calls of 
  ! this routine
  ! What are the integration weights ??? 
  ! I cannot use monte-carlo for wavelength integration
  ! because lambda is fixed for each line.
  ! or generate a grid and interpolate the value on this
  ! random grid ?
 !------------------------------------------------- !
  integer, intent(in) :: id, icell, iray, n_rayons
  logical, intent(in), optional :: switch_to_ALI
  integer :: nact, nc, kr, switch,i, j, Nblue, Nred, jp, krr
  type (AtomType), pointer :: atom
  double precision, dimension(NLTEspec%Nwaves) :: Ieff
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2
  double precision :: norm = 0d0, hc_4PI, diag
  integer :: c_nlte = 0
    
  !FLAG to change from MALI to ALI
  if (present(switch_to_ALI)) then
   if (switch_to_ALI) switch = 0
  else
   switch = 1
  end if
  write(*,*) "Check Norm"
  stop
  Ieff = 0d0 !part of the wavelength-angle integration
  hc_4PI = HPLANCK * CLIGHT / 4d0 / PI
  
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   Ieff(:) = NLTEspec%I(:,iray,id) - NLTEspec%Psi(:,iray,id) * atom%eta(:,id)!at this iray, id
   !loop over transitions, b-f and b-b for these atoms
   !To do; define a transition_type with either cont or line
   do kr=1,atom%Ncont
    norm = sum(cont_wlam(atom%continua(kr))) * n_rayons
    i = atom%continua(kr)%i; j = atom%continua(kr)%j
    Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred 
    twohnu3_c2 = 2d0 * HPLANCK * CLIGHT / (NLTEspec%lambda*NM_TO_M)**3.
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%continua(kr)%Vij(:,id)*Ieff(Nblue:Nred)*cont_wlam(atom%continua(kr))) / norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%continua(kr)%Vij(:,id)*atom%continua(kr)%gij(:,id)*cont_wlam(atom%continua(kr))) / norm
     !write(*,*) atom%ID, icell, id, i, j, atom%Gamma(i, j), atom%Gamma(j,i)

    ! cross-coupling terms
     atom%Gamma(i,j) = atom%Gamma(i,j) -c_nlte*sum((atom%chi_up(i,:,id)*NLTEspec%Psi(:, iray, id)*&
     	atom%Uji_down(j,:,id))*cont_wlam(atom%continua(kr))) / norm
     !write(*,*) atom%ID, icell, id, i, j, atom%Gamma(i, j), atom%Gamma(j,i)

     ! check if i is an upper level of another transition
     do krr=1,atom%Ncont
      jp = atom%continua(krr)%j
      if (jp==i) then !i upper level of this transition
       atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum((atom%chi_down(j,:,id)*nLTEspec%Psi(:, iray, id)*&
       	atom%Uji_down(i,:,id))*cont_wlam(atom%continua(kr))) / norm
      end if 
     end do
     !write(*,*) atom%ID, icell, id, i, j, atom%Gamma(i, j), atom%Gamma(j,i)

   end do

   do kr=1,atom%Nline
    i = atom%lines(kr)%i; j = atom%lines(kr)%j
    norm = 1d0 * HC_4PI* n_rayons
    Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred 
    twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji 
    
    !Ieff (Uij = 0 'cause i<j)
    atom%Gamma(i,j) = atom%Gamma(i,j) + c_nlte*sum(atom%lines(kr)%Vij(:,id)*Ieff(Nblue:Nred)*atom%lines(kr)%wlam(:)) / norm
    !Uji + Vji*Ieff
    !Uji and Vji express with Vij
    atom%Gamma(j,i) = atom%Gamma(j,i) +  c_nlte*sum((Ieff(Nblue:Nred)+twohnu3_c2(Nblue:Nred))*&
    	atom%lines(kr)%Vij(:,id)*atom%lines(kr)%gij(:,id)*atom%lines(kr)%wlam(:))
    	
    ! cross-coupling terms
     atom%Gamma(i,j) = atom%Gamma(i,j) -c_nlte*sum(atom%chi_up(i,Nblue:Nred,id)*NLTEspec%Psi(Nblue:Nred, iray, id)*&
     	atom%Uji_down(j,Nblue:Nred,id)*atom%lines(kr)%wlam(:)) / norm
     ! check if i is an upper level of another transition
     do krr=1,atom%Nline
      jp = atom%lines(krr)%j
      if (jp==i) then !i upper level of this transition
       atom%Gamma(j,i) = atom%Gamma(j,i) + c_nlte*sum(atom%chi_down(j,Nblue:Nred,id)*nLTEspec%Psi(Nblue:Nred, iray, id)*&
       	atom%Uji_down(i,Nblue:Nred,id)*atom%lines(kr)%wlam(:)) / norm
      end if 
     end do
    ! 
   end do
   
   !remove diagonal here, delat(l',l)
   do i=1,atom%Nlevel
    atom%Gamma(i,i) = 0d0 !because i that case Gamma_l"l = Cii - Cii = 0d0
    diag = 0d0
    do j=1,atom%Nlevel
     diag = diag + atom%Gamma(j,i) 
    end do
    atom%Gamma(i,i) = -diag
   end do
   
   NULLIFY(atom)
  end do !loop over atoms
 RETURN
 END SUBROUTINE fillGamma_MALI

 SUBROUTINE Gamma_LTE(id,icell)
 !------------------------------------------------- !
 !------------------------------------------------- !
  integer, intent(in) :: id, icell
  integer :: nact, nc, kr, switch,i, j, Nblue, Nred, jp, krr
  type (AtomType), pointer :: atom
  double precision :: diag

  
  do nact=1,atmos%Nactiveatoms !loop over each active atoms
   atom => atmos%ActiveAtoms(nact)%ptr_atom
   !remove diagonal here, delat(l',l)
   do i=1,atom%Nlevel
    atom%Gamma(i,i) = 0d0 !because i that case Gamma_l"l = Cii - Cii = 0d0
    diag = 0d0
    do j=1,atom%Nlevel !Sum ober lprime prime for each l
     diag = diag + atom%Gamma(j,i) 
    end do
    atom%Gamma(i,i) = -diag !only if l=lprime
   end do
   
   NULLIFY(atom)
  end do !loop over atoms
 RETURN
 END SUBROUTINE Gamma_LTE

 SUBROUTINE SEE_atom(id, icell, atom)
  integer, intent(in) :: icell, id
  type(AtomType), intent(inout) :: atom
  integer :: j, i, imaxpop
  
  !search the row with maximum population and remove it
  imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
  atom%n(:,icell) = 0d0
  atom%n(imaxpop,icell) = atom%ntotal(icell)
  
  
  atom%Gamma(imaxpop,:) = 1d0
!   do i=1, atom%Nlevel
!    do j=1, atom%Nlevel
!    write(*,*) i, j, atom%gamma(i, j)
!    end do
!   end do
!   stop
  CALL GaussSlv(atom%Gamma, atom%n(:,icell),atom%Nlevel)
if (minval(atom%n(:,icell)) <= 0d0) then
 write(*,*) "toto"
 write(*,*) atom%n(:,icell)
 write(*,*) atom%Gamma
 stop
end if
 RETURN
 END SUBROUTINE SEE_atom
 
 SUBROUTINE updatePopulations(id, icell)
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
