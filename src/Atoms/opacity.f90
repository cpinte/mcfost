MODULE Opacity

 use atmos_type, only : atmos, hydrogen
 use atom_type
 use spectrum_type, only : NLTEspec, initAtomOpac, init_psi_operator!,init_Xcoupling
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

   NLTEspec%dtau(:,iray,id) = chi(:)*ds !J = Sum_ray I0*exp(-dtau) + Psi*S
   NLTEspec%Psi(:,iray,id) = ((1d0 - dexp(-NLTEspec%dtau(:,iray,id))))/(chi+1d-300)


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
  CALL initAtomOpac(id)
  CALL init_psi_operator(id, iray)
!   NLTEspec%Psi(:,iray,id) = 0d0; NLTEspec%dtau(:,iray,id) = 0d0
  !set atom%eta to zero also
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

  norm = 5d-1 / line%lambda0 * CLIGHT !because we want dv
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
 
 !building, not ready
 SUBROUTINE calc_J_coherent(id, icell, n_rayons)
  integer, intent(in) :: id, icell, n_rayons
  
  NLTEspec%J(:,icell) = sum(NLTEspec%I(:,1:n_rayons,id),dim=2)
  NLTEspec%Jc(:,icell) = sum(NLTEspec%Ic(:,1:n_rayons,id),dim=2)

 
 RETURN
 END SUBROUTINE calc_J_coherent
 
 SUBROUTINE NLTEOpacity(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)
  !
  !
  ! chi = Vij * (ni - gij * nj) = ni * Vij - nj * Vji
  ! eta = twohnu3_c2 * gij * Vij * nj = Uji * nj = 2hnu3/c2 * Vji * nj
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
  real(kind=dp) :: gij, twohnu3_c2, stm
  double precision, dimension(:), allocatable :: Vij, gijk, twohnu3_c2k
  double precision, allocatable :: phiZ(:,:), psiZ(:,:)
  
  
  atom_loop : do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom
   
   if (iterate) aatom%eta(:,iray,id) = 0d0
   
   	tr_loop : do kr = 1, aatom%Ntr
   	
   	    !!not useful now. Transition not appearing for the image are removed in aatom%at
   	    !!and for the NLTEloop all transitions are kept ATM.
   	    !!if (.not.aatom%at(kr)%lcontrib_to_opac) cycle
   	
        stm = 1d0 
        kc = aatom%at(kr)%ik !relative index of a transition among continua or lines
write(*,*) "kc=",kc, id, icell, iray, "stm=",stm,"trtype=",aatom%at(kr)%trtype
        !if (.not.aatom%at(kr)%contrib_to_opac) CYCLE
        
        SELECT CASE (aatom%at(kr)%trtype)
        
        CASE ('ATOMIC_CONTINUUM')
write(*,*) "cont loc", loc(cont)
    	  cont = aatom%continua(kc)
    	  !!if (.not.cont%lcontrib_to_opac) CYCLE tr_loop

    	  Nred = cont%Nred;Nblue = cont%Nblue    	
    	  i = cont%i;j = cont%j
if ((i>aatom%Nlevel .or. i < 1).or.(j>aatom%Nlevel.or.j<1).or.(id<1.or.id>NLTEspec%NPROC).or.&
	(icell<1.or.icell>atmos%Nspace)) then
 write(*,*) "bf", id, icell, i, j
 stop
end if
        
    	  if (aatom%n(j,icell) < tiny_dp .or. aatom%n(i,icell) < tiny_dp) then !test
    	   write(*,*) id, icell, loc(aatom), loc(cont)
    	   write(*,*) tiny_dp, aatom%n(j,icell), aatom%n(i,icell)
    	   write(*,*) aatom%n(:,icell)
    	 
          if (aatom%n(j,icell)==0d0 .or. aatom%n(i,icell)==0d0) then
           write(*,*) icell, iray, id, aatom%ID, aatom%Nlevel, kc, shape(aatom%n)
           write(*,*) i, cont%i, j, cont%j
           write(*,*) aatom%n(:,icell)
           write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%n(cont%i,icell), aatom%n(cont%j,icell)
           write(*,*) "1", aatom%n(1,icell), "2", aatom%n(2,icell), "3", aatom%n(3,icell),"4", aatom%n(4,icell)
           stop
          end if    	 
    	 
     	  CALL WARNING("too small cont populations")
     	  aatom%n(j,icell) = max(tiny_dp, aatom%n(j,icell))
     	  aatom%n(i,icell) = max(tiny_dp, aatom%n(i,icell))
    	  end if !end test
    	  
      	  allocate(gijk(cont%Nlambda))
          gijk(:) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
      
         !allocate Vij, to avoid computing bound_free_Xsection(cont) 3 times for a continuum
	     allocate(Vij(cont%Nlambda), twohnu3_c2k(cont%Nlambda))
	     
    	 Vij(:) = bound_free_Xsection(cont) 	
    	       !write(*,*) 'bf', icell, id, iray, maxval(Vij), minval(Vij)

         twohnu3_c2k(:) = twohc / NLTEspec%lambda(Nblue:Nred)**(3d0) * &
              Vij(:) * gijk(:) * aatom%n(j,icell) !eta
         
         Vij(:) = Vij(:) * (aatom%n(i,icell) - stm * gijk(:)*aatom%n(j,icell)) !chi


         !keep copy of continuum only
    	 NLTEspec%AtomOpac%chic_nlte(Nblue:Nred, id) = NLTEspec%AtomOpac%chic_nlte(Nblue:Nred, id) + Vij(:)
    	 NLTEspec%AtomOpac%etac_nlte(Nblue:Nred, id) = NLTEspec%AtomOpac%etac_nlte(Nblue:Nred, id) + twohnu3_c2k(:)

         !gather total NLTE
         NLTEspec%AtomOpac%chi(Nblue:Nred,id) = NLTEspec%AtomOpac%chi(Nblue:Nred,id) + Vij(:)
       		
         NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + twohnu3_c2k(:)
    	
    	 if (iterate) then
    	   aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + twohnu3_c2k(:)
         end if			
   
        deallocate(Vij, gijk, twohnu3_c2k)
   
       CASE ('ATOMIC_LINE')
write(*,*) 'line loc', loc(line)
        line = aatom%lines(kc)
        !!if (.not.line%lcontrib_to_opac) CYCLE tr_loop
        Nred = line%Nred;Nblue = line%Nblue
        i = line%i;j = line%j
if ((i>aatom%Nlevel .or. i < 1).or.(j>aatom%Nlevel.or.j<1).or.(id<1.or.id>NLTEspec%NPROC).or.&
	(icell<1.or.icell>atmos%Nspace)) then
 write(*,*) "bb", id, icell, iray,  id, icell, i, j
 stop
end if
    
    if ((aatom%n(j,icell) < tiny_dp).or.(aatom%n(i,icell) < tiny_dp)) then !no transition
        write(*,*) id, icell, loc(aatom), loc(line)
    	write(*,*) tiny_dp, aatom%n(j, icell), aatom%n(i,icell)
        write(*,*) aatom%n(:,icell)
        
        if (aatom%n(j,icell)==0d0 .or. aatom%n(i,icell)==0d0) then
         write(*,*) icell, iray, id, aatom%ID, aatom%Nlevel, kc, shape(aatom%n)
         write(*,*) i, line%i, j, line%j
         write(*,*) aatom%n(:,icell)
         write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%n(line%i,icell), aatom%n(line%j,icell)
         write(*,*) "1", aatom%n(1,icell), "2", aatom%n(2,icell), "3", aatom%n(3,icell),"4", aatom%n(4,icell)
         stop
        end if
     	CALL WARNING("too small line populations")
     	aatom%n(j,icell) = max(tiny_dp, aatom%n(j,icell))
     	aatom%n(i,icell) = max(tiny_dp, aatom%n(i,icell))
    end if 

        gij = line%Bji / line%Bij

        twohnu3_c2 = line%Aji / line%Bji
        if (line%voigt)  CALL Damping(icell, aatom, kc, line%adamp)
        if (line%adamp>5.) write(*,*) " large damping for line", line%j, line%i, line%atom%ID, line%adamp
    
        !allocate(phi(line%Nlambda),Vij(line%Nlambda))
        allocate(Vij(line%Nlambda))
    
        if (PRT_SOLUTION=="FULL_STOKES") allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
        !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
        !they are allocated in this case.
        CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l,Vij,phiZ,psiZ)!, phi, phiZ, psiZ)

        if (iterate) aatom%lines(kc)%phi(:,iray,id) = Vij(:)!phi(:)
         Vij(:) = hc_4PI * line%Bij * Vij(:)!phi(:) !normalized in Profile()
                                                             ! / (SQRTPI * VBROAD_atom(icell,aatom)) 
            !write(*,*) 'bb', icell, id, iray, maxval(Vij), minval(Vij)

         !opacity total
         NLTEspec%AtomOpac%chi(Nblue:Nred,id) = NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		Vij(:) * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell))
       		
         NLTEspec%AtomOpac%eta(Nblue:Nred,id)= NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
       		twohnu3_c2 * gij * Vij(:) * aatom%n(j,icell)
      
        if (iterate) then
           aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
              twohnu3_c2 * gij * Vij(:) * aatom%n(j,icell)
        end if	

    
        if (line%polarizable .and. PRT_SOLUTION == "FULL_STOKES") then
         write(*,*) "Beware, NLTE part of Zeeman opac not set to 0 between iteration!"
         do nk = 1, 3
          !magneto-optical
          NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) + &
           hc_4PI * line%Bij * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell)) * psiZ(nk,:)
          !dichroism
          NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) + &
           hc_4PI * line%Bij * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell)) * psiZ(nk,:)
          !emissivity
          NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) + &
          twohnu3_c2 * gij * hc_4PI * line%Bij * aatom%n(j,icell) * phiZ(nk,:)
         end do 
        end if
     
       deallocate(Vij)!, phi)
       if (PRT_SOLUTION=="FULL_STOKES") deallocate(phiZ, psiZ)
    
    CASE DEFAULT
     CALL Error("Transition type unknown", aatom%at(kr)%trtype)
    END SELECT
    
   end do tr_loop
     !write(*,*) icell, id, "eta=", maxval(aatom%eta(:,iray,id)), minval(aatom%eta(:,iray,id))

   aatom => NULL()
   !stop
  end do atom_loop

 RETURN
 END SUBROUTINE NLTEOpacity

END MODULE Opacity