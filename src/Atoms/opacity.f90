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
 use profiles, only : Profile
 !!use molecular_emission, only : v_proj
 use math, only : locate, integrate_dx


 IMPLICIT NONE

 !store the pure continuum NLTE opacities to be added to the total
 !continuum opacity after NLTEloop ends, if lstore_opac add them to Kc and jc
 !double precision, dimension(:,:), allocatable :: chic_nlte, etac_nlte


 CONTAINS
 
  SUBROUTINE add_to_psi_operator(id, icell, iray, ds)
  ! ----------------------------------------------------------- !
   ! Computes Psi and Ieff at the cell icell, for the thread id
   ! in the direction iray, using ds path length of the ray.
  ! ----------------------------------------------------------- !
   integer, intent(in) :: iray, id, icell
   double precision, intent(in) :: ds
   double precision, dimension(NLTEspec%Nwaves) :: chi!, eta_loc
  
   if (lstore_opac) then
     chi(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1)
     !eta_loc(:) = NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id) + &
     !               NLTEspec%AtomOpac%jc(icell,:)
   else
     chi(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)
     !eta_loc(:) = NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)
   end if !store_opac
   !Choose eval operatore depending on the method

   SELECT CASE (atmos%nLTE_methode)
     CASE ("MALI")
       NLTEspec%dtau(:,iray,id) = 0d0 !not used in SEE
       NLTEspec%Psi(:,iray,id) = (1d0 - dexp(-chi(:)*ds)) / chi !in meter
     CASE ("HOGEREIJDE")
       NLTEspec%dtau(:,iray,id) = chi(:)*ds !J = Sum_ray I0*exp(-dtau) + Psi*S
       !NLTEspec%Psi(:,iray,id) = ((1d0 - dexp(-NLTEspec%dtau(:,iray,id))))* eta_loc/(chi+1d-300)
       NLTEspec%Psi(:,iray,id) = ((1d0 - dexp(-NLTEspec%dtau(:,iray,id))))/(chi+1d-300)
       !or do not use atom%eta, and Psi = Psi * eta_total
     CASE DEFAULT
       CALL error(" In evaluate Psi operators!", atmos%nLTE_methode)
   END SELECT

  RETURN
  END SUBROUTINE add_to_psi_operator

 
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
  ! gives dv/c = dlambda/lambda = dnu/nu
  ! times c: dv.
  ! the integral of the radiative rates is
  ! integ (domega) * integ(dv/ch)
  ! a 1/ch is missing. 
  !
  ! Vij = hnu/4pi * Bij * phi if integral is over dnu/hnu
  ! and Vij = hc/4pi * Bij * phi if integral is over (dv/hc).
  !
  ! phi in s in the former case, and phi in s/m in the later.
 ! --------------------------------------------------------- !
  type(AtomicLine), intent(in) :: line
  double precision, dimension(line%Nlambda) :: wlam
  integer :: la, Nblue, Nred, la_start, la_end, la0
  double precision :: norm !beware this is not the result of the integral 
  						   ! just a factor to convert dv/c from dlambda/lambda
   !la0: index of wavelengths on the frequency grid (size Nwaves). 
   !la:  index of wavelengths on the lambda grid of the line (size Nlambda).
   !la0 = Nblue - 1 + la; line expands from Nblue to Nred on the frequency grid.
   !la=1 <=> la0=Nblue; la=Nlambda <=> la0 = Nred = Nblue - 1 + Nlambda
   !dlambda = (lambda(la0 + 1) - lambda(la0 - 1)) * 0.5 <=> mean value.

  norm = 5d-1 * CLIGHT/line%lambda0
  Nblue = line%Nblue; Nred = line%Nred
  la_start = 1; la_end = line%Nlambda


  wlam(1) = (NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) * norm
  
  wlam(line%Nlambda) = (NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) * norm

  do la=2,line%Nlambda-1

   la0 = Nblue - 1 + la
   wlam(la) = (NLTEspec%lambda(la0 + 1)-NLTEspec%lambda(la0 - 1)) * norm

  end do

 RETURN
 END FUNCTION line_wlam
 
 FUNCTION cont_wlam(cont) result(wlam)
 ! --------------------------------------------------------- !
  ! computes dlam/lam for a continnum 
  ! dnu/nu = dlam/lam
  ! the integral of the radiative rates is
  ! integ (domega) * integ(dlam/hlam)
  ! a 1/h is missing
 ! --------------------------------------------------------- !
  type(AtomicContinuum), intent(in) :: cont
  double precision, dimension(cont%Nlambda) :: wlam
  integer :: la, Nblue, Nred, la_start, la_end , la0

  Nblue = cont%Nblue; Nred = cont%Nred
  la_start = 1; la_end = cont%Nlambda

  wlam(1) = 5d-1 * &
  	(NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) / NLTEspec%lambda(Nblue)
  	
  wlam(cont%Nlambda) = 5d-1 * & 
  	(NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) / NLTEspec%lambda(Nred)
  	
  do la=2,cont%Nlambda-1
  
   la0 = Nblue - 1 + la
   wlam(la) = 5d-1 * &
   	(NLTEspec%lambda(la0+1)-NLTEspec%lambda(la0-1)) / NLTEspec%lambda(la0)
   	
  end do

 RETURN
 END FUNCTION cont_wlam

 !add options to save %eta, %gij %Vij only during the NLTEloop. For the image we
 !do not need them as we do not update the pops
 !same condition as get_weights
 !Je n'ai besoin de gij, Vij qu'a l'indice de la cellule
 SUBROUTINE NLTEOpacity(id, icell, x, y, z, x1, y1, z1, u, v, w, l, iterate)
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
  ! if iterate, compute lines weight for this cell icell and rays and eta, Vij gij for that atom.
  ! if not iterate means that atom%gij atom%vij atom%eta are not allocated (after NLTE for image for instance)
  !
  integer, intent(in) :: id, icell
  double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
  logical, intent(in) :: iterate
  integer :: nact, Nred, Nblue, kc, kr, i, j, nk
  type(AtomicLine) :: line
  type(AtomicContinuum) :: cont
  type(AtomType), pointer :: aatom
  double precision, parameter :: twohc = (2. * HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
  double precision, parameter :: hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
  double precision, dimension(NLTEspec%Nwaves) :: Vij, exp_lambda, gij, twohnu3_c2
  double precision, allocatable :: phi(:), phiZ(:,:), psiZ(:,:)
  character(len=20) :: VoigtMethod="HUMLICEK"
  double precision :: hc_4PI
  
  exp_lambda = dexp(-hc_k / (NLTEspec%lambda * atmos%T(icell)))
  twohnu3_c2 = twohc / NLTEspec%lambda(:)**(3d0)
  hc_4PI = HPLANCK * CLIGHT / (4d0 * PI)
  
  do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom
   if (iterate) aatom%eta(:,id) = 0d0 !init Eta only if iterate, for this cell and thread

   
   	do kc = 1, aatom%Ncont
    	cont = aatom%continua(kc)
    	Nred = cont%Nred; Nblue = cont%Nblue
    	!Should never happen. Because for NLTEOpac all transitions are present.
    	!However, in the case of images, some transitions are removed. But anyway
    	!
        if (Nred == -99 .and. Nblue == -99) CALL ERROR("NLTEOPAC")

    	i = cont%i; j=cont%j
        gij = 0d0
   	 	Vij(Nblue:Nred) = cont%alpha(Nblue:Nred)
    	if (aatom%n(j,icell) < tiny_dp) then
    	 write(*,*) aatom%n(j,icell)
     	 CALL Warning("too small cont populations") !or Error()
     	 CYCLE
    	end if
    	
    	 gij(Nblue:Nred) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * exp_lambda(Nblue:Nred)

    
    !store total emissivities and opacities
       NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
       		
		NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
    	gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)
    	
    	if (iterate) then
    	 aatom%continua(kc)%gij(:,id) = gij(Nblue:Nred)
    	 aatom%continua(kc)%Vij(:,id) = Vij(Nblue:Nred)
    	 aatom%eta(Nblue:Nred,id) = aatom%eta(Nblue:Nred,id) + &
    		gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)
        !! aatom%continua(kc)%wlam(Nblue:Nred) = 
    	end if
    !Do not forget to add continuum opacities to the all continnum opacities
    !after all populations have been converged    
   	end do

   do kr = 1, aatom%Nline
    line = aatom%lines(kr)
    Nred = line%Nred; Nblue = line%Nblue
    !if (Nred == -99 .and. Nblue == -99) CYCLE

    i = line%i; j=line%j
    
    if ((aatom%n(j,icell) < tiny_dp).or.(aatom%n(i,icell) < tiny_dp)) then !no transition
    	write(*,*) tiny_dp
        write(*,*) aatom%n(:,icell)
     	CALL Warning("too small line populations") !or Error()
     	CYCLE
    end if 
    gij = 0d0
    Vij = 0d0

    gij(:) = line%Bji / line%Bij
    twohnu3_c2(Nblue:Nred) = line%Aji / line%Bji
    if (line%voigt)  CALL Damping(icell, aatom, kr, line%adamp)
    if (line%adamp>5.) write(*,*) " large damping for line", line%j, line%i, line%atom%ID, line%adamp
    allocate(phi(line%Nlambda))
    if (PRT_SOLUTION=="FULL_STOKES") &
    	allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
    !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
    !they are allocated in this case.
    CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)


     Vij(Nblue:Nred) = hc_4PI * line%Bij * phi(:) !normalized in Profile()
                                                             ! / (SQRTPI * aatom%vbroad(icell)) 
      
     NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
       		
     NLTEspec%AtomOpac%eta(Nblue:Nred,id)= &
     		NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
       		twohnu3_c2(Nblue:Nred) * gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell)
      
    if (iterate) then
     aatom%lines(kr)%gij(:,id) = gij(Nblue:Nred)
     aatom%lines(kr)%Vij(:,id) = Vij(Nblue:Nred)    		
     aatom%eta(Nblue:Nred,id) = aatom%eta(Nblue:Nred,id) + &
    	twohnu3_c2(Nblue:Nred) * gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell)

     ! unit is m/s / (m/s * s/m) = m/s !it is iray dependent as it is computed along each direction
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

 SUBROUTINE NLTEOpacity_lambda(la, id, icell, x, y, z, x1, y1, z1, u, v, w, l)
  !Even if it is NLTEOpac, it means Opacity for Active atoms.
  !It is called at the end of the NLTEloop
  ! no pol yet.
  integer, intent(in) :: id, icell, la
  double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
  integer :: nact, Nred, Nblue, kc, kr, i, j, nk, Nlambda_line
  type(AtomicLine) :: line
  type(AtomicContinuum) :: cont
  type(AtomType), pointer :: aatom
  double precision, parameter :: twohc = (2. * HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
  double precision, parameter :: hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
  double precision, dimension(1) :: Vij, exp_lambda, gij, twohnu3_c2
  double precision, allocatable :: phi(:), phiZ(:,:), psiZ(:,:)
  character(len=20) :: VoigtMethod="HUMLICEK"
  double precision :: hc_4PI
  
  exp_lambda(1) = dexp(-hc_k / (NLTEspec%lambda(la) * atmos%T(icell)))
  twohnu3_c2(1) = twohc / NLTEspec%lambda(1)**(3d0)
  hc_4PI = HPLANCK * CLIGHT / (4d0 * PI)
  
  do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom

   
   	do kc = 1, aatom%Ncont
    	cont = aatom%continua(kc)
    	Nred = cont%Nred; Nblue = cont%Nblue
     !wavelength does not fall inside line domaine? OK cylcle
     if ((NLTEspec%lambda(la) < NLTEspec%lambda(Nblue)).or.&
        (NLTEspec%lambda(la) > NLTEspec%lambda(Nred))) then
       CYCLE
     end if

    	i = cont%i; j=cont%j
        gij = 0d0
   	 	Vij(1) = cont%alpha(la)
    	if (aatom%n(j,icell) < tiny_dp) then
    	 write(*,*) aatom%n(j,icell)
     	 CALL Warning("too small cont populations") !or Error()
     	 CYCLE
    	end if
    	
    	 gij(1) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * exp_lambda(1)

    
    !store total emissivities and opacities
       NLTEspec%AtomOpac%chi(la,id) = &
     		NLTEspec%AtomOpac%chi(la,id) + &
       		Vij(1) * (aatom%n(i,icell)-gij(1)*aatom%n(j,icell))
       		
		NLTEspec%AtomOpac%eta(la,id) = NLTEspec%AtomOpac%eta(la,id) + &
    	gij(1) * Vij(1) * aatom%n(j,icell) * twohnu3_c2(1)
    	
    !Do not forget to add continuum opacities to the all continnum opacities
    !after all populations have been converged    
   	end do

   do kr = 1, aatom%Nline
    line = aatom%lines(kr)
    Nred = line%Nred; Nblue = line%Nblue

     if ((NLTEspec%lambda(la) < NLTEspec%lambda(Nblue)).or.&
        (NLTEspec%lambda(la) > NLTEspec%lambda(Nred))) then
       CYCLE
     end if

    i = line%i; j=line%j
    
    if ((aatom%n(j,icell) <=tiny_dp).or.(aatom%n(i,icell) <=tiny_dp)) then !no transition
        write(*,*) aatom%n(:,icell)
     	CALL Warning("too small line populations") !or Error()
     	CYCLE
    end if 
    gij = 0d0
    Vij = 0d0

    gij(1) = line%Bji / line%Bij
    twohnu3_c2(1) = line%Aji / line%Bji
    if (line%voigt)  CALL Damping(icell, aatom, kr, line%adamp)
    Nlambda_line = line%Nlambda
    line%Nlambda = 1
    allocate(phi(line%Nlambda))

    !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
    !they are allocated in this case.
    CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)
    line%Nlambda = Nlambda_line

     Vij(1) = hc_4PI * line%Bij * phi(1) !normalized in Profile()
                                                             ! / (SQRTPI * aatom%vbroad(icell)) 
      
     NLTEspec%AtomOpac%chi(la,id) = &
     		NLTEspec%AtomOpac%chi(la,id) + &
       		Vij(1) * (aatom%n(i,icell)-gij(1)*aatom%n(j,icell))
       		
     NLTEspec%AtomOpac%eta(la,id)= &
     		NLTEspec%AtomOpac%eta(la,id) + &
       		twohnu3_c2(1) * gij(1) * Vij(1) * aatom%n(j,icell)
      
    deallocate(phi)
   end do
  
  end do !over activeatoms

 RETURN
 END SUBROUTINE NLTEOpacity_lambda

END MODULE Opacity
