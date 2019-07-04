MODULE Opacity

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type, only : NLTEspec, initAtomOpac, init_Xcoupling 
 use constant
 use constantes, only				 : tiny_dp, huge_dp, AU_to_m
 use messages
 !!use voigtfunctions, only 			 : Voigt
 use broad, only 					 : Damping
 use parametres
 use profiles, only : Profile
 use metal, only : bound_free_Xsection, Background, BackgroundLines
 !!use molecular_emission, only : v_proj
 use math, only : locate, integrate_dx
 use grid, only : cross_cell



 IMPLICIT NONE

 !store the pure continuum NLTE opacities to be added to the total
 !continuum opacity after NLTEloop ends, if lstore_opac add them to Kc and jc
 !They are temporary stored in a large array (Nspace, Nlambda)

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
       CALL Error("Not implemented MALI") !should create an error before, in init_NLTE
       !!not allocated
       !!NLTEspec%dtau(:,iray,id) = 0d0 !not used in SEE
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
!   NLTEspec%AtomOpac%chi(:,id) = 0d0
!   NLTEspec%AtomOpac%eta(:,id) = 0d0
  !We need to recompute LTE opacities for this cell and ray
  CALL initAtomOpac(id, .true.)
!   NLTEspec%Psi(:,iray,id) = 0d0; NLTEspec%dtau(:,iray,id) = 0d0
  !set atom%eta to zero also
  !CALL initCrossCoupling(id)
  !NOTE Zeeman opacities are not re set to zero and are accumulated
  !change that or always use FIELD_FREE
  !is equivalent to P(icell,id, iray) ?
  !Compute opacity, eta and chi for this cell in the direction u0, v0, w0
  CALL NLTEOpacity(id, icell, iray, x0, y0, z0, x1, x2, x3, u0, v0, w0, l_dum, .true.)
  if (lstore_opac) then
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, x2, x3, u0, v0, w0, l_dum)
  else
      CALL Background(id, icell, x0, y0, z0, x1, x2, x3, u0, v0, w0, l_dum)
  end if
  !last .true. is to compute atom%gij, atom%Uji,atom%Vij 
  !CALL fillCrossCoupling_terms(id, icell)
  !il faut recalculer Psi dans les sous-iterations et Ieff si Hogereijde.
  !Sinon, seulement Psi depend de chi a besoin d'être mis à jour.
  ds = l_dum * AU_to_m
  !recompute Psi and eventually Ieff.
  CALL add_to_psi_operator(id, icell, iray, ds)
  
 RETURN
 END SUBROUTINE init_local_field_atom 
  
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

  norm = 5d-1 / line%lambda0 !* CLIGHT, removed because we want dv/c not dv
  Nblue = line%Nblue; Nred = line%Nred
  la_start = 1; la_end = line%Nlambda


  wlam(1) = (NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) * norm
  
  wlam(line%Nlambda) = (NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) * norm
  !write(*,*) 1, wlam(1)
  do la=2,line%Nlambda-1

   la0 = Nblue - 1 + la
   wlam(la) = (NLTEspec%lambda(la0 + 1)-NLTEspec%lambda(la0 - 1)) * norm
   !write(*,*) la, wlam(la)
  end do
  !write(*,*) line%Nlambda, wlam(line%Nlambda)

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

  !Nblue = cont%Nblue; Nred = cont%Nred
  Nblue = 1; Nred = NLTEspec%Nwaves !---> Because ATM cont are kept on the whole grid
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
 SUBROUTINE NLTEOpacity(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)
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
  integer, intent(in) :: id, icell, iray
  double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
  logical, intent(in) :: iterate
  integer :: nact, Nred, Nblue, kc, kr, i, j, nk
  type(AtomicLine) :: line
  type(AtomicContinuum) :: cont
  type(AtomType), pointer :: aatom
  real(kind=dp) :: gij, twohnu3_c2
  double precision, dimension(:), allocatable :: Vij, gijk, twohnu3_c2k
  double precision, allocatable :: phi(:), phiZ(:,:), psiZ(:,:)
  character(len=20) :: VoigtMethod="HUMLICEK"


  do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom
   if (iterate) aatom%eta(:,iray,id) = 0d0 !init Eta only if iterate, for this cell and thread

   
   	do kc = 1, aatom%Ncont
    	cont = aatom%continua(kc)
    	Nred = cont%Nred; Nblue = cont%Nblue    	
    	i = cont%i; j=cont%j

    	if (aatom%n(j,icell) <= tiny_dp .or. aatom%n(i,icell) <= tiny_dp) then
    	 write(*,*) aatom%n(j,icell), aatom%n(i,icell)
     	 CALL Warning("too small cont populations") !or Error()
     	 ! if (aatom%n(i, icell) <= tiny_dp) CYCLE !but do not skip if atom%j < tiny_dp
!      	 aatom%n(j,icell) = 0d0 !here because it is < tiny_dp, otherwise it is n(i) and we skip
    	end if
      	allocate(gijk(cont%Nlambda))
        gijk(:) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))


      !Cannot be negative because we alread tested if < tiny_dp
      if ((aatom%n(i,icell) <= minval(gijk)*aatom%n(j,icell)).or.&
        (aatom%n(i,icell) <= maxval(gijk)*aatom%n(j,icell))) then
          write(*,*) id, icell, aatom%ID, &
          	" ** Stimulated emission for continuum transition ",j,i," neglected"
!          write(*,*) cont%i, cont%j, aatom%ID
!          write(*,*) "at cell-1", atmos%nHtot(icell-1), atmos%nHtot(icell-1)*aatom%Abund, atmos%T(icell-1), atmos%ne(icell-1), atmos%nHmin(icell-1)
!          write(*,*) "at cell", atmos%nHtot(icell), atmos%nHtot(icell)*aatom%Abund, atmos%T(icell), atmos%ne(icell), atmos%nHmin(icell)
!          write(*,*) id, icell, i, j, aatom%n(i,icell), aatom%n(j,icell)
!          write(*,*) minval(gijk)*aatom%n(j,icell), maxval(gijk)*aatom%n(j,icell)
!          write(*,*) "nstar =", (aatom%nstar(nk,icell), nk=1,aatom%Nlevel)
!          write(*,*) "n     =", (aatom%n(nk,icell), nk=1,aatom%Nlevel)
!                   stop
!         deallocate(gijk)
!        CYCLE
      end if
	  allocate(Vij(cont%Nlambda), twohnu3_c2k(cont%Nlambda))
    	
      twohnu3_c2k(:) = twohc / NLTEspec%lambda(cont%Nblue:cont%Nred)**(3d0)
   	  Vij(:) = bound_free_Xsection(cont) 	

    
    !store total emissivities and opacities
        NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(:) * (aatom%n(i,icell)-gijk(:)*aatom%n(j,icell))
       		
		NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
    	gijk(:) * Vij(:) * aatom%n(j,icell) * twohnu3_c2k
    	
    	!They are allocated if nlte so if we enter here
!     	NLTEspec%AtomOpac%chic_nlte(icell, Nblue:Nred) = &
!     	 NLTEspec%AtomOpac%chic_nlte(icell, Nblue:Nred) + &
!     	 	Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
!     	NLTEspec%AtomOpac%etac_nlte(icell, Nblue:Nred) = &
!     	 NLTEspec%AtomOpac%etac_nlte(icell, Nblue:Nred) + &
!     	   gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)

    	
    	if (iterate) &
         aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
         							gijk(:) * Vij(:) * aatom%n(j,icell) * twohnu3_c2k

    !Do not forget to add continuum opacities to the all continnum opacities
    !after all populations have been converged    
     deallocate(Vij, gijk, twohnu3_c2k)
   	end do

   do kr = 1, aatom%Nline
    line = aatom%lines(kr)
    Nred = line%Nred; Nblue = line%Nblue
    !if (Nred == -99 .and. Nblue == -99) CYCLE
    i = line%i; j=line%j
    
    if ((aatom%n(j,icell) <= tiny_dp).or.(aatom%n(i,icell) <= tiny_dp)) then !no transition
    	write(*,*) tiny_dp, aatom%n(j, icell), aatom%n(i,icell)
        write(*,*) aatom%n(:,icell)
     	CALL Warning("too small line populations") !or Error()
    end if 

    gij = line%Bji / line%Bij !array of constant Bji/Bij
    !Cannot be negative because we alread tested if < tiny_dp
    if ((aatom%n(i,icell) <= gij*aatom%n(j,icell)).or.&
        (aatom%n(i,icell) <= gij*aatom%n(j,icell))) then
          write(*,*) id, icell, aatom%ID, &
          	" ** Stimulated emission for line transition ",j,i," neglected"
!          write(*,*) id, icell, i, j, aatom%n(i,icell), aatom%n(j,icell)
!          write(*,*) gij
!          write(*,*) "nstar =", (aatom%nstar(nk,icell), nk=1,aatom%Nlevel)
!          write(*,*) "n     =", (aatom%n(nk,icell), nk=1,aatom%Nlevel)
!         stop
!	      CYCLE
    end if
    
    twohnu3_c2 = line%Aji / line%Bji
    if (line%voigt)  CALL Damping(icell, aatom, kr, line%adamp)
    if (line%adamp>5.) write(*,*) " large damping for line", line%j, line%i, line%atom%ID, line%adamp
    allocate(phi(line%Nlambda),Vij(line%Nlambda))
    if (PRT_SOLUTION=="FULL_STOKES") allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
    !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
    !they are allocated in this case.
    CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)


     Vij(:) = hc_4PI * line%Bij * phi(:) !normalized in Profile()
                                                             ! / (SQRTPI * VBROAD_atom(icell,aatom)) 
      
     NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(:) * (aatom%n(i,icell)-gij*aatom%n(j,icell))
       		
     NLTEspec%AtomOpac%eta(Nblue:Nred,id)= &
     		NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
       		twohnu3_c2 * gij * Vij(:) * aatom%n(j,icell)
      
    !line and cont are not pointers. Modification of line does not affect atom%lines(kr)
    if (iterate) then
      aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
      								twohnu3_c2 * gij * Vij(:) * aatom%n(j,icell)
      aatom%lines(kr)%phi(:,iray,id) = phi(:)
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
     
    deallocate(phi, Vij)
    if (PRT_SOLUTION=="FULL_STOKES") deallocate(phiZ, psiZ)
   end do
  
  end do !over activeatoms

 RETURN
 END SUBROUTINE NLTEOpacity

!  SUBROUTINE NLTEOpacity_lambda(la, id, icell, x, y, z, x1, y1, z1, u, v, w, l)
!   !Even if it is NLTEOpac, it means Opacity for Active atoms.
!   !It is called at the end of the NLTEloop
!   ! no pol yet.
!   integer, intent(in) :: id, icell, la
!   double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
!   integer :: nact, Nred, Nblue, kc, kr, i, j, nk, Nlambda_line
!   type(AtomicLine) :: line
!   type(AtomicContinuum) :: cont
!   type(AtomType), pointer :: aatom
!   double precision, parameter :: twohc = (2. * HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
!   double precision, parameter :: hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
!   double precision, dimension(1) :: Vij, exp_lambda, gij, twohnu3_c2
!   double precision, allocatable :: phi(:), phiZ(:,:), psiZ(:,:)
!   character(len=20) :: VoigtMethod="HUMLICEK"
!   double precision :: hc_4PI
!   
!   exp_lambda(1) = dexp(-hc_k / (NLTEspec%lambda(la) * atmos%T(icell)))
!   twohnu3_c2(1) = twohc / NLTEspec%lambda(1)**(3d0)
!   hc_4PI = HPLANCK * CLIGHT / (4d0 * PI)
!   
!   do nact = 1, atmos%Nactiveatoms
!    aatom => atmos%ActiveAtoms(nact)%ptr_atom
! 
!    
!    	do kc = 1, aatom%Ncont
!     	cont = aatom%continua(kc)
!     	Nred = cont%Nred; Nblue = cont%Nblue
!      !wavelength does not fall inside line domaine? OK cylcle
!      if ((NLTEspec%lambda(la) < NLTEspec%lambda(Nblue)).or.&
!         (NLTEspec%lambda(la) > NLTEspec%lambda(Nred))) then
!        CYCLE
!      end if
! 
!     	i = cont%i; j=cont%j
!         gij = 0d0
!    	 	Vij(1) = cont%alpha(la)
!     	if (aatom%n(j,icell) < tiny_dp) then
!     	 write(*,*) aatom%n(j,icell)
!      	 CALL Warning("too small cont populations") !or Error()
!      	 CYCLE
!     	end if
!     	
!     	 gij(1) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * exp_lambda(1)
! 
!     
!     !store total emissivities and opacities
!        NLTEspec%AtomOpac%chi(la,id) = &
!      		NLTEspec%AtomOpac%chi(la,id) + &
!        		Vij(1) * (aatom%n(i,icell)-gij(1)*aatom%n(j,icell))
!        		
! 		NLTEspec%AtomOpac%eta(la,id) = NLTEspec%AtomOpac%eta(la,id) + &
!     	gij(1) * Vij(1) * aatom%n(j,icell) * twohnu3_c2(1)
!     	
!     !Do not forget to add continuum opacities to the all continnum opacities
!     !after all populations have been converged    
!    	end do
! 
!    do kr = 1, aatom%Nline
!     line = aatom%lines(kr)
!     Nred = line%Nred; Nblue = line%Nblue
! 
!      if ((NLTEspec%lambda(la) < NLTEspec%lambda(Nblue)).or.&
!         (NLTEspec%lambda(la) > NLTEspec%lambda(Nred))) then
!        CYCLE
!      end if
! 
!     i = line%i; j=line%j
!     
!     if ((aatom%n(j,icell) <=tiny_dp).or.(aatom%n(i,icell) <=tiny_dp)) then !no transition
!         write(*,*) aatom%n(:,icell)
!      	CALL Warning("too small line populations") !or Error()
!      	CYCLE
!     end if 
!     gij = 0d0
!     Vij = 0d0
! 
!     gij(1) = line%Bji / line%Bij
!     twohnu3_c2(1) = line%Aji / line%Bji
!     if (line%voigt)  CALL Damping(icell, aatom, kr, line%adamp)
!     Nlambda_line = line%Nlambda
!     line%Nlambda = 1
!     allocate(phi(line%Nlambda))
! 
!     !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
!     !they are allocated in this case.
!     CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)
!     line%Nlambda = Nlambda_line
! 
!      Vij(1) = hc_4PI * line%Bij * phi(1) !normalized in Profile()
!                                                              ! / (SQRTPI * VBROAD_atom(icell,aatom)) 
!       
!      NLTEspec%AtomOpac%chi(la,id) = &
!      		NLTEspec%AtomOpac%chi(la,id) + &
!        		Vij(1) * (aatom%n(i,icell)-gij(1)*aatom%n(j,icell))
!        		
!      NLTEspec%AtomOpac%eta(la,id)= &
!      		NLTEspec%AtomOpac%eta(la,id) + &
!        		twohnu3_c2(1) * gij(1) * Vij(1) * aatom%n(j,icell)
!       
!     deallocate(phi)
!    end do
!   
!   end do !over activeatoms
! 
!  RETURN
!  END SUBROUTINE NLTEOpacity_lambda
!  SUBROUTINE NLTEOpacity(id, icell, x, y, z, x1, y1, z1, u, v, w, l, iterate)
!   !
!   !
!   ! chi = Vij * (ni - gij * nj)
!   ! eta = twohnu3_c2 * gij * Vij * nj
!   ! Continuum:
!   ! Vij = alpha
!   ! gij = nstari/nstarj * exp(-hc/kT/lamba)
!   ! Lines:
!   ! twoHnu3/c2 = Aji/Bji
!   ! gij = Bji/Bij (*rho if exists)
!   ! Vij = Bij * hc/4PI * phi
!   !
!   ! if iterate, compute lines weight for this cell icell and rays and eta, Vij gij for that atom.
!   ! if not iterate means that atom%gij atom%vij atom%eta are not allocated (after NLTE for image for instance)
!   !
!   integer, intent(in) :: id, icell
!   double precision, intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
!   logical, intent(in) :: iterate
!   integer :: nact, Nred, Nblue, kc, kr, i, j, nk
!   type(AtomicLine) :: line
!   type(AtomicContinuum) :: cont
!   type(AtomType), pointer :: aatom
!   double precision, dimension(NLTEspec%Nwaves) :: Vij, exp_lambda, gij, twohnu3_c2
!   double precision, allocatable :: phi(:), phiZ(:,:), psiZ(:,:)
!   character(len=20) :: VoigtMethod="HUMLICEK"
!   
!   exp_lambda = dexp(-hc_k / (NLTEspec%lambda * atmos%T(icell)))
!   twohnu3_c2 = twohc / NLTEspec%lambda(:)**(3d0)
! 
!   
!   do nact = 1, atmos%Nactiveatoms
!    aatom => atmos%ActiveAtoms(nact)%ptr_atom
!    if (iterate) aatom%eta(:,id) = 0d0 !init Eta only if iterate, for this cell and thread
! 
!    
!    	do kc = 1, aatom%Ncont
!     	cont = aatom%continua(kc)
!     	Nred = cont%Nred; Nblue = cont%Nblue
!     	!Should never happen. Because for NLTEOpac all transitions are present.
!     	!However, in the case of images, some transitions are removed. But anyway
!     	!
!         !if (Nred == -99 .and. Nblue == -99) CALL ERROR("NLTEOPAC")
! 
!     	i = cont%i; j=cont%j
!         gij = 0d0
!         Vij = 0d0
!    	 	Vij(:) = cont%alpha(:) !Nblue:Nred
!     	if (aatom%n(j,icell) <= tiny_dp .or. aatom%n(i,icell) <= tiny_dp) then
!     	 write(*,*) aatom%n(j,icell), aatom%n(i,icell)
!      	 CALL Warning("too small cont populations") !or Error()
!      	 ! if (aatom%n(i, icell) <= tiny_dp) CYCLE !but do not skip if atom%j < tiny_dp
! !      	 aatom%n(j,icell) = 0d0 !here because it is < tiny_dp, otherwise it is n(i) and we skip
!     	end if
!     	
! !       if (aatom%nstar(j,icell) <= tiny_dp) then 
! !        gij(Nblue:Nred) = 0d0
! !       else	 
!        gij(Nblue:Nred) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * exp_lambda(Nblue:Nred)
! !       end if
! !       write(*,*) " --------------------------------------------------------- "
! !       write(*,*) "-->", id, icell,i, j, aatom%n(i,icell), aatom%n(j,icell)
! !       write(*,*) "nstar =", (aatom%nstar(nk,icell), nk=1,aatom%Nlevel)
! !       write(*,*) "n     =", (aatom%n(nk,icell), nk=1,aatom%Nlevel)
!       !write(*,*) minval(gij(Nblue:Nred)), maxval(gij(Nblue:Nred))
! 
!       !Cannot be negative because we alread tested if < tiny_dp
!       if ((aatom%n(i,icell) <= minval(gij(Nblue:Nred))*aatom%n(j,icell)).or.&
!         (aatom%n(i,icell) <= maxval(gij(Nblue:Nred))*aatom%n(j,icell))) then
!          write(*,*) " ** Stimulated emission for continuum transition, chi < 0", cont%i, cont%j, &
!           aatom%ID
!          write(*,*) id, icell, i, j, aatom%n(i,icell), aatom%n(j,icell)
!          write(*,*) minval(gij(Nblue:Nred)), maxval(gij(Nblue:Nred))
!          write(*,*) "nstar =", (aatom%nstar(nk,icell), nk=1,aatom%Nlevel)
!          write(*,*) "n     =", (aatom%n(nk,icell), nk=1,aatom%Nlevel)
!         CYCLE
!       end if
! !       write(*,*) " --------------------------------------------------------- "
! 
!     
!     !store total emissivities and opacities
!         NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
!      		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
!        		Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
!        		
! 		NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
!     	gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)
!     	
!     	!They are allocated if nlte so if we enter here
! !     	NLTEspec%AtomOpac%chic_nlte(icell, Nblue:Nred) = &
! !     	 NLTEspec%AtomOpac%chic_nlte(icell, Nblue:Nred) + &
! !     	 	Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
! !     	NLTEspec%AtomOpac%etac_nlte(icell, Nblue:Nred) = &
! !     	 NLTEspec%AtomOpac%etac_nlte(icell, Nblue:Nred) + &
! !     	   gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell) * twohnu3_c2(Nblue:Nred)
! 
!     	
!     	if (iterate) then
!     	 aatom%continua(kc)%gij(:,id) = gij!(Nblue:Nred)
!     	 aatom%continua(kc)%Vij(:,id) = Vij!(Nblue:Nred)
!     	 aatom%eta(:,id) = aatom%eta(:,id) + &!Nblue:Nred
!     		gij(:) * Vij(:) * aatom%n(j,icell) * twohnu3_c2(:)
!     	end if
!     !Do not forget to add continuum opacities to the all continnum opacities
!     !after all populations have been converged    
!    	end do
! 
!    do kr = 1, aatom%Nline
!     line = aatom%lines(kr)
!     Nred = line%Nred; Nblue = line%Nblue
!     !if (Nred == -99 .and. Nblue == -99) CYCLE
! 
!     i = line%i; j=line%j
!     
!     if ((aatom%n(j,icell) <= tiny_dp).or.(aatom%n(i,icell) <= tiny_dp)) then !no transition
!     	write(*,*) tiny_dp, aatom%n(j, icell), aatom%n(i,icell)
!         write(*,*) aatom%n(:,icell)
!      	CALL Warning("too small line populations") !or Error()
!     end if 
!     gij = 0d0
!     Vij = 0d0
! 
!     gij(:) = line%Bji / line%Bij !array of constant Bji/Bij
!     !Cannot be negative because we alread tested if < tiny_dp
!     if ((aatom%n(i,icell) <= minval(gij(Nblue:Nred))*aatom%n(j,icell)).or.&
!         (aatom%n(i,icell) <= maxval(gij(Nblue:Nred))*aatom%n(j,icell))) then
!          write(*,*) " ** Stimulated emission for line transition, chi < 0", line%i, line%j, &
!           aatom%ID
!          write(*,*) id, icell, i, j, aatom%n(i,icell), aatom%n(j,icell)
!          write(*,*) minval(gij(Nblue:Nred)), maxval(gij(Nblue:Nred))
!          write(*,*) "nstar =", (aatom%nstar(nk,icell), nk=1,aatom%Nlevel)
!          write(*,*) "n     =", (aatom%n(nk,icell), nk=1,aatom%Nlevel)
!         stop
!     end if
!     
!     twohnu3_c2(Nblue:Nred) = line%Aji / line%Bji
!     if (line%voigt)  CALL Damping(icell, aatom, kr, line%adamp)
!     if (line%adamp>5.) write(*,*) " large damping for line", line%j, line%i, line%atom%ID, line%adamp
!     allocate(phi(line%Nlambda))
!     if (PRT_SOLUTION=="FULL_STOKES") &
!     	allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
!     !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
!     !they are allocated in this case.
!     CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)
! 
! 
!      Vij(Nblue:Nred) = hc_4PI * line%Bij * phi(:) !normalized in Profile()
!                                                              ! / (SQRTPI * VBROAD_atom(icell,aatom)) 
!       
!      NLTEspec%AtomOpac%chi(Nblue:Nred,id) = &
!      		NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
!        		Vij(Nblue:Nred) * (aatom%n(i,icell)-gij(Nblue:Nred)*aatom%n(j,icell))
!        		
!      NLTEspec%AtomOpac%eta(Nblue:Nred,id)= &
!      		NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
!        		twohnu3_c2(Nblue:Nred) * gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell)
!       
!     if (iterate) then
!      aatom%lines(kr)%gij(:,id) = gij(Nblue:Nred)
!      aatom%lines(kr)%Vij(:,id) = Vij(Nblue:Nred)    		
!      aatom%eta(Nblue:Nred,id) = aatom%eta(Nblue:Nred,id) + &
!     	twohnu3_c2(Nblue:Nred) * gij(Nblue:Nred) * Vij(Nblue:Nred) * aatom%n(j,icell)
! 
!      ! unit is m/s / (m/s * s/m) = m/s !it is iray dependent as it is computed along each direction
!      !aatom%lines(kr)%wlam(:) = line_wlam(aatom%lines(kr)) / sum(phi*line_wlam(aatom%lines(kr)))
!      aatom%lines(kr)%wlam_norm(id) = 1d0 / sum(phi*line_wlam(aatom%lines(kr)))
!      end if
!     
!      if (line%polarizable .and. PRT_SOLUTION == "FULL_STOKES") then
!        write(*,*) "Beware, NLTE part of Zeeman opac not set to 0 between iteration!"
!        do nk = 1, 3
!          !magneto-optical
!          NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) + &
!            hc_4PI * line%Bij * (aatom%n(i,icell)-gij*aatom%n(j,icell)) * psiZ(nk,:)
!          !dichroism
!          NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) + &
!            hc_4PI * line%Bij * (aatom%n(i,icell)-gij*aatom%n(j,icell)) * psiZ(nk,:)
!          !emissivity
!          NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) + &
!           twohnu3_c2 * gij * hc_4PI * line%Bij * aatom%n(j,icell) * phiZ(nk,:)
!        end do 
!      end if
!     deallocate(phi)
!     if (PRT_SOLUTION=="FULL_STOKES") deallocate(phiZ, psiZ)
!    end do
!   
!   end do !over activeatoms
! 
!  RETURN
!  END SUBROUTINE NLTEOpacity
END MODULE Opacity
