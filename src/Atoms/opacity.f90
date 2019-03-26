MODULE Opacity

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type, only : NLTEspec
 use constant
 use constantes, only				 : tiny_dp, huge_dp
 use messages
 !!use voigtfunctions, only 			 : Voigt
 use broad, only 					 : Damping
 use parametres
 use profiles, only : Profile, Profile_lambda
 !!use molecular_emission, only : v_proj
 use math, only : locate


 IMPLICIT NONE

 !store the pure continuum NLTE opacities to be added to the total
 !continuum opacity after NLTEloop ends, if lstore_opac add them to Kc and jc
 !double precision, dimension(:,:), allocatable :: chic_nlte, etac_nlte


 CONTAINS
 
  SUBROUTINE alloc_wlambda()
  ! --------------------------------------------------- !
   ! Allocates wavelength integration
   ! weights.
   ! only from Nblue to Nred !
   ! weights for line are direction dependent because
   ! of the absorption profile phi
   ! For continua, the weights are computed on the fly
   ! reducing the amount of memory needed.
   ! For lines we compute them in NLTEopac() for icell
   ! and iray, because of phi.
  ! --------------------------------------------------- !
   use atmos_type, only : atmos
   integer :: kr, kc, nact, Nred, Nblue, Nlambda, la
   
   do nact=1,atmos%Nactiveatoms
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
      Nred = atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%Nred
      Nblue = atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%Nblue
      Nlambda = atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%Nlambda
      allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%wlam(Nlambda))
      atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%wlam(:) = 0d0
!       do la=1,Nlambda!la=Nblue, Nred 
!        atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%wlam(la) = & 
!         NLTEspec%lambda(la+Nblue-1)-NLTEspec%lambda(la+Nblue-1-1)
!       end do
    end do
!    do kc=1,atmos%ActiveAtoms(nact)%ptr_atom%Ncont
!      Nred = atmos%ActiveAtoms(nact)%ptr_atom%continua(kc)%Nred
!      Nblue = atmos%ActiveAtoms(nact)%ptr_atom%continua(kc)%Nblue
!      Nlambda = atmos%ActiveAtoms(nact)%ptr_atom%continua(kc)%Nlambda
!      allocate(atmos%ActiveAtoms(nact)%ptr_atom%continua(kc)%wlam(Nlambda))
!      atmos%ActiveAtoms(nact)%ptr_atom%continua(kc)%wlam(:) = 0d0
!!       do la=1,Nlambda!=Nblue, Nred
!!        atmos%ActiveAtoms(nact)%ptr_atom%continua(kc)%wlam(la) = & 
!!         NLTEspec%lambda(la+Nblue-1)-NLTEspec%lambda(la+Nblue-1-1)
!!       end do
!    end do  
   end do
 
  RETURN
  END SUBROUTINE alloc_wlambda  
  
 FUNCTION line_wlam(line) result(wlam)
 ! --------------------------------------------------------- !
  ! The integral of I(nu, omega) over frequencies and
  ! angles is:
  ! int_omega int_nu I * phi(nu, omega) * domega* dnu/nu
  ! or: (with integral boundaries reversed)
  ! int_omega domega int_lambda I*phi  * dlambda/lambda
  ! or:
  ! int_omega domega int_v I*phi * dv * lambda0/lambda /c
 ! --------------------------------------------------------- !
  type(AtomicLine), intent(in) :: line
  double precision, dimension(line%Nlambda) :: wlam
  integer :: la, Nblue, Nred, la_start, la_end
  double precision :: norm = CLIGHT/HPLANCK

  Nblue = line%Nblue; Nred = line%Nred
  la_start = 1; la_end = line%Nlambda
  !compute
  !check that we are not at the boundaries, otherwise we cannot go beyond
  if (Nblue==1) then 
   la_start = 2
   wlam(1) = line%lambda0*(NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) / NLTEspec%lambda(Nblue) * norm
  end if
  
  if (Nred==NLTEspec%Nwaves) then
   la_end = line%Nlambda-1
   wlam(line%Nlambda) = line%lambda0*(NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) / NLTEspec%lambda(Nred) * norm
  end if
  do la=1,line%Nlambda! !first is Nblue - (Nblue -1) which exists
  					!last is Nred - (Nred-1)
   wlam(la) = line%lambda0*(NLTEspec%lambda(la+Nblue-1)-NLTEspec%lambda(la+Nblue-1-1)) / NLTEspec%lambda(la+Nblue-1) * norm
  end do
 
 RETURN
 END FUNCTION line_wlam
 
 FUNCTION cont_wlam(cont) result(wlam)
 ! --------------------------------------------------------- !
 ! computes dlam/lam for a continnum 
 ! --------------------------------------------------------- !
  type(AtomicContinuum), intent(in) :: cont
  double precision, dimension(cont%Nlambda) :: wlam
  integer :: la, Nblue, Nred, la_start, la_end
  double precision :: norm = 1d0 / HPLANCK
  

  Nblue = cont%Nblue; Nred = cont%Nred
  la_start = 1; la_end = cont%Nlambda
  !compute
  !check that we are not at the boundaries, otherwise we cannot go beyond
  if (Nblue==1) then 
   la_start = 2
   wlam(1) = (NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) / NLTEspec%lambda(Nblue) * norm
  end if
  
  if (Nred==NLTEspec%Nwaves) then
   la_end = cont%Nlambda-1
   wlam(cont%Nlambda) = (NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) / NLTEspec%lambda(Nred) * norm
  end if
  do la=la_start,la_end !first is Nblue - (Nblue -1) which exists
  					!last is Nred - (Nred-1)
   
   wlam(la) = (NLTEspec%lambda(la+Nblue-1)-NLTEspec%lambda(la+Nblue-1-1)) / NLTEspec%lambda(la+Nblue-1) * norm
   
  end do
 
 RETURN
 END FUNCTION cont_wlam

 SUBROUTINE NLTEOpacity(id, icell, x, y, z, x1, y1, z1, u, v, w, l, get_weights)
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
  ! if get_weights, compute lines weight for this cell icell and rays
  !
  integer, intent(in) :: id, icell
  double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
  logical, intent(in) :: get_weights
  integer :: nact, Nred, Nblue, kc, kr, i, j, nk
  type(AtomicLine) :: line
  type(AtomicContinuum) :: cont
  type(AtomType), pointer :: aatom
  double precision, parameter :: twohc = (2. * HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
  double precision, parameter :: hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
  double precision, dimension(NLTEspec%Nwaves) :: Vij, exp_lambda, gij, twohnu3_c2
  double precision, allocatable :: phi(:), phiZ(:,:), psiZ(:,:)
  integer, parameter :: NvspaceMax = 100
  character(len=20) :: VoigtMethod="HUMLICEK"
  double precision :: hc_4PI
  
  exp_lambda = dexp(-hc_k / (NLTEspec%lambda * atmos%T(icell)))
  twohnu3_c2 = twohc / NLTEspec%lambda(:)**(3d0)
  hc_4PI = HPLANCK * CLIGHT / (4d0 * PI)
  
  do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom
   aatom%eta(:,id) = 0d0

   
   	do kc = 1, aatom%Ncont
    	cont = aatom%continua(kc)
    	Nred = cont%Nred; Nblue = cont%Nblue
    	i = cont%i; j=cont%j
        gij = 0d0
   	 	Vij(Nblue:Nred) = cont%alpha(Nblue:Nred)
    	if (aatom%n(j,icell) < tiny_dp) then
    	 write(*,*) aatom%n(j,icell)
     	 CALL error("too small populations")
     	 CYCLE
    	end if
    	gij(Nblue:Nred) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * exp_lambda(Nblue:Nred)
    	aatom%continua(kc)%gij(:,id) = gij(Nblue:Nred)
    	aatom%continua(kc)%Vij(:,id) = Vij(Nblue:Nred)

    
    !store total emissivities and opacities
       NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
       		
		NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
    	gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)
    	
    	aatom%eta(Nblue:Nred,id) = aatom%eta(Nblue:Nred,id) + &
    		gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)
        !! aatom%continua(kc)%wlam(Nblue:Nred) = hc_4PI
    !Do not forget to add continuum opacities to the all continnum opacities
    !after all populations have been converged    
   	end do

   do kr = 1, aatom%Nline
    line = aatom%lines(kr)
    Nred = line%Nred; Nblue = line%Nblue
    i = line%i; j=line%j
    
    if ((aatom%n(j,icell) <=tiny_dp).or.(aatom%n(i,icell) <=tiny_dp)) then !no transition
        write(*,*) aatom%n(:,icell)
     	CALL error("too small populations")
     	CYCLE
    end if 
    gij = 0d0
    Vij = 0d0

    gij(:) = line%Bji / line%Bij
    twohnu3_c2(Nblue:Nred) = line%Aji / line%Bji
    if (line%voigt)  CALL Damping(icell, aatom, kr, line%adamp)
    allocate(phi(line%Nlambda))
    if (PRT_SOLUTION=="FULL_STOKES") &
    	allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
    CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)


     Vij(Nblue:Nred) = hc_4PI * line%Bij * phi(:) !normalized in Profile()
                                                             ! / (SQRTPI * aatom%vbroad(icell))
    
     aatom%lines(kr)%gij(:,id) = gij(Nblue:Nred)
     aatom%lines(kr)%Vij(:,id) = Vij(Nblue:Nred)    
      
     NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
       		
     NLTEspec%AtomOpac%eta(Nblue:Nred,id)= &
     		NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
       		twohnu3_c2(Nblue:Nred) * gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell)
       		
    aatom%eta(Nblue:Nred,id) = aatom%eta(Nblue:Nred,id) + &
    	twohnu3_c2(Nblue:Nred) * gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell)
    	
     if (get_weights) then
      !normalized
      aatom%lines(kr)%wlam(:) = line_wlam(aatom%lines(kr)) / sum(phi*line_wlam(aatom%lines(kr)))
     end if
    
     if (line%polarizable .and. PRT_SOLUTION == "FULL_STOKES") then
       write(*,*) "Beware, NLTE part of Zeeman opac not set to 0 between iteration!"
       do nk = 1, 3
         !magneto-optical
         NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) + &
           hc_4PI * line%Bij * (aatom%n(i,icell)-gij*aatom%n(j,icell)) * psiZ(nk,:)
         !dichroism
         NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) + &
           hc_4PI * line%Bij * (aatom%n(i,icell)-gij*aatom%n(j,icell)) * psiZ(nk,:)
         !emissivity
         NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) + &
          twohnu3_c2 * gij * hc_4PI * line%Bij * aatom%n(j,icell) * phiZ(nk,:)
       end do 
     end if
    deallocate(phi)
    if (PRT_SOLUTION=="FULL_STOKES") deallocate(phiZ, psiZ)
   end do
  
  end do !over activeatoms

 RETURN
 END SUBROUTINE NLTEOpacity


END MODULE Opacity
