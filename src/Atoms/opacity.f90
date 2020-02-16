MODULE Opacity

	use atmos_type, only								: atmos, hydrogen, ntotal_atom
	use atom_type
	use spectrum_type, only								: NLTEspec, initAtomOpac, init_psi_operator, initAtomOpac_nlte
	use constant
	use constantes, only								: tiny_dp, huge_dp, AU_to_m
	use messages
	use broad, only										: Damping
	use parametres
	use profiles, only									: Profile
	use getlambda, only									: hv

	use math, only										: locate, integrate_nu, integrate_x
	use mcfost_env, only								: dp
	use input, only										: ds
	use molecular_emission, only						: v_proj

 


	IMPLICIT NONE
 
	real(kind=dp), dimension(:,:,:), allocatable :: chi_loc, eta_loc
	!real(kind=dp), dimension(:,:), allocatable :: chic_nlte_loc, etac_nlte_loc
	real(kind=dp), private, parameter :: tiny_chi = 1d-50
 
	CONTAINS
 
	SUBROUTINE calc_psi_operator_m(id, icell, iray, chi, dtau)
	!Mali
	! I' = (I - eta) Psi + Psi* eta with constant background opacities during succesive iterations

		integer, intent(in) :: iray, id, icell
		real(kind=dp), dimension(NLTEspec%Nwaves), intent(in) :: chi, dtau
   
		NLTEspec%Psi(:,iray,id) = (1d0 - dexp(-dtau)) / chi
   
		NLTEspec%etau(:,iray,id) = dexp(-dtau)


	RETURN
	END SUBROUTINE calc_psi_operator_m
	
	
	SUBROUTINE calc_psi_operator(id, icell, iray, chi, dtau, S)
	
	! I' = I * etau + S * Psi

		integer, intent(in) :: iray, id, icell
		real(kind=dp), dimension(NLTEspec%Nwaves), intent(in) :: chi, dtau, S
   
		NLTEspec%Psi(:,iray,id) = (1d0 - dexp(-dtau)) 
   
		NLTEspec%etau(:,iray,id) = dexp(-dtau)
		
		NLTEspec%S(:,iray,id) = S(:)


	RETURN
	END SUBROUTINE calc_psi_operator


	SUBROUTINE calc_Jnu(id, icell, n_rayons)
		integer, intent(in) :: id, icell, n_rayons
		integer :: iray
		
		!TMP because %J is not used yet.
		NLTEspec%Jc(:,icell) = 0.0
		if (allocated(NLTEspec%J)) then
			NLTEspec%J(:,icell) = 0.0
			do iray=1,n_rayons	
				NLTEspec%J(:,icell) = NLTEspec%J(:,icell) + NLTEspec%I(:,iray,id) / n_rayons
				NLTEspec%Jc(:,icell) = NLTEspec%Jc(:,icell) + NLTEspec%Ic(:,iray,id) / n_rayons
			enddo
		else
			do iray=1,n_rayons	
				NLTEspec%Jc(:,icell) = NLTEspec%Jc(:,icell) + NLTEspec%Ic(:,iray,id) / n_rayons
			enddo
		endif
  

	RETURN
	END SUBROUTINE calc_Jnu  

 
   
	!Does not compute etac for id, because they are computed on the whole grid for each global 
	!iterations. They are updated locally only for sub iterations
	SUBROUTINE NLTE_bound_free(icell)
		integer, intent(in) :: icell
		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la
		type(AtomType), pointer :: aatom
  
		NLTEspec%AtomOpac%Kc_nlte(:,icell) = 0d0
		NLTEspec%AtomOpac%jc_nlte(:,icell) = 0d0

		atom_loop : do nact = 1, atmos%Nactiveatoms
			aatom => atmos%ActiveAtoms(nact)%ptr_atom
   
		!loop only on continua
			tr_loop : do kr = aatom%Ntr_line+1,aatom%Ntr


        		kc = aatom%at(kr)%ik !relative index of a transition among continua or lines


				Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
				i = aatom%continua(kc)%i; j = aatom%continua(kc)%j

          
         		do la=1,aatom%continua(kc)%Nlambda

					if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0) then
					
						NLTEspec%AtomOpac%Kc_nlte(Nblue+la-1,icell) = NLTEspec%AtomOpac%Kc_nlte(Nblue+la-1,icell) + &
            aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell))
            
						NLTEspec%AtomOpac%jc_nlte(Nblue+la-1,icell) = NLTEspec%AtomOpac%jc_nlte(Nblue+la-1,icell) + &
            aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
            
					endif
          
				enddo   


			end do tr_loop
   


			aatom => NULL()
		end do atom_loop
	
 
	RETURN
	END SUBROUTINE NLTE_bound_free
 
	SUBROUTINE compute_nlte_bound_free
	!Does not update eta_c
	!Special version of Flux/image
	!$ use omp_lib
		integer :: icell, id


		!$omp parallel &
		!$omp default(none) &
		!$omp private(icell,id) &
		!$omp shared(atmos)
		!$omp do schedule(dynamic,1)
		do icell=1,atmos%Nspace
			!$ id = omp_get_thread_num() + 1
			if (atmos%icompute_atomRT(icell) > 0) then
				CALL NLTE_bound_free(icell)
			endif
		end do
		!$omp end do
		!$omp end parallel
 
	RETURN
	END SUBROUTINE compute_nlte_bound_free
	
	!Done in the propagation, to have it at id (icell) for iray==1
	SUBROUTINE calc_etac_atom_loc(id, icell)
	!Store the continuous emissivity for the running cell, for all atoms
		integer, intent(in) :: icell, id
		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la
		type(AtomType), pointer :: aatom
  

		atom_loop : do nact = 1, atmos%Nactiveatoms
			aatom => atmos%ActiveAtoms(nact)%ptr_atom
   
		!loop only on continua
			tr_loop : do kr = aatom%Ntr_line+1,aatom%Ntr


        		kc = aatom%at(kr)%ik !relative index of a transition among continua or lines


				Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
				i = aatom%continua(kc)%i; j = aatom%continua(kc)%j

          
         		do la=1,aatom%continua(kc)%Nlambda

					if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0) then
            
						aatom%etac(Nblue+la-1,id) = aatom%etac(Nblue+la-1,id) + &
            aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
            
					endif
          
				enddo   


			end do tr_loop
   


			aatom => NULL()
		end do atom_loop

		
	RETURN
	END SUBROUTINE calc_etac_atom_loc

 
	SUBROUTINE calc_eta_atom_loc(id, icell, iray, initialize) !update eta and psi, iterate always true
	! ------------------------------------------------------------------------- !
	! ------------------------------------------------------------------------- !  
  
		integer, intent(in) :: id, icell, iray
		logical, intent(in) :: initialize
		integer :: kr, kc, i, j, Nred, Nblue, la, nact
		real(kind=dp), dimension(NLTEspec%Nwaves) :: dtau, chi
		type (AtomType), pointer :: aatom
		integer :: dk

		!LTE, unchanged opacities for local sub iterations
		chi(:) = chi_loc(:,iray,id)
  
		atom_loop : do nact = 1, atmos%Nactiveatoms
			aatom => atmos%ActiveAtoms(nact)%ptr_atom
			
			!before continua, because we need to add etac for all rays in eta.
			if (initialize) then
				aatom%chic(:,id) = 0.0 !local but we want to keep chic for all rays to add it to the total chi of each ray
			
				tr_loopc : do kr = aatom%Ntr_line+1,aatom%Ntr


        			kc = aatom%at(kr)%ik


					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j

          
         			do la=1,aatom%continua(kc)%Nlambda

						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0) then

							aatom%chic(Nblue+la-1,id) = aatom%chic(Nblue+la-1,id)+  aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell))
							aatom%etac(Nblue+la-1,id) = aatom%etac(Nblue+la-1, id) + aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
						
						endif
            
          
					enddo   


				end do tr_loopc
						
			endif
			
			!etac is zero at iray==1, computed and then added for each ray (at this point, non zero)
			aatom%eta(:,iray,id) = aatom%etac(:,id)
			chi(:) = chi(:) + aatom%chic(:,id) !constant over direction, added to chi_loc(:,iray,id)
   
			tr_loop : do kr = 1, aatom%Ntr_line

   	
				kc = aatom%at(kr)%ik

				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
				dk = aatom%lines(kr)%dk(iray,id)

				if ((aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0) cycle

         		!NLTE subroutine. uses dk instead of a profile
				do la=1, aatom%lines(kc)%Nlambda

					chi(Nblue+la-1-dk) = chi(Nblue+la-1-dk) + &
					hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell) - aatom%lines(kc)%gij*aatom%n(j,icell))
         
					aatom%eta(Nblue+la-1-dk,iray,id) = aatom%eta(Nblue+la-1-dk,iray,id) + &
					aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)

				enddo
    
			end do tr_loop

			aatom => NULL()
		end do atom_loop
		
		! Check that chi_loc is correct
        !opacity total = LTE (doesn't change during sub iterations) + NLTE
		dtau = chi(:) * ds(iray,id)

		CALL calc_psi_operator_m(id, icell, iray, chi, dtau)


	RETURN
	END SUBROUTINE calc_eta_atom_loc
 
	SUBROUTINE calc_total_source_loc(id, icell, iray, initialize) !update eta and psi, iterate always true
	! ------------------------------------------------------------------------- !
	! ------------------------------------------------------------------------- !  
  
		integer, intent(in) :: id, icell, iray
		logical, intent(in) :: initialize
		integer :: kr, kc, i, j, Nred, Nblue, la, nact
		!real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, dtau
		type (AtomType), pointer :: aatom
		integer :: dk

		!LTE, unchanged opacities for local sub iterations
		NLTEspec%chi(:,iray,id) = chi_loc(:,iray,id)
		NLTEspec%S(:,iray,id) = eta_loc(:,iray,id)!! + &
			! Is always new during the sub iterations, because computed after integrate_ray_line
			!!NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell)
			!-> But still, should use the old scattering opacity for local iterations ??
  
		atom_loop : do nact = 1, atmos%Nactiveatoms
			aatom => atmos%ActiveAtoms(nact)%ptr_atom
			
			!before continua, because we need to add etac for all rays in eta.
			!if it is too long I reintroduce etac and chic to compute only for one ray. (see calc_eta_atom_loc)
			!if (initialize) then
				!aatom%chic(:,id) = 0.0 !local but we want to keep chic for all rays to add it to the total chi of each ray
			
				tr_loopc : do kr = aatom%Ntr_line+1,aatom%Ntr


        			kc = aatom%at(kr)%ik


					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j

          
         			do la=1,aatom%continua(kc)%Nlambda

						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0) then

							!aatom%chic(Nblue+la-1,id) = aatom%chic(Nblue+la-1,id)+  aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell))
							!aatom%etac(Nblue+la-1,id) = aatom%etac(Nblue+la-1, id) + aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
							NLTEspec%chi(Nblue+la-1,iray,id) = NLTEspec%chi(Nblue+la-1,iray,id) + aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell))
							NLTEspec%S(Nblue+la-1,iray,id) = NLTEspec%S(Nblue+la-1,iray,id) + aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
						
						endif
            
          
					enddo   


				end do tr_loopc
						
			!endif
			
			!chi(:) = chi(:) + aatom%chic(:,id)
			!eta(:) = eta(:) + aatom%etac(:,id)
   
			tr_loop : do kr = 1, aatom%Ntr_line

   	
				kc = aatom%at(kr)%ik

				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
				dk = aatom%lines(kr)%dk(iray,id)

				if ((aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0) cycle

         		!NLTE subroutine. uses dk instead of a profile
				do la=1, aatom%lines(kc)%Nlambda

					NLTEspec%chi(Nblue+la-1-dk,iray,id) = NLTEspec%chi(Nblue+la-1-dk,iray,id) + &
					hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell) - aatom%lines(kc)%gij*aatom%n(j,icell))
         
					NLTEspec%S(Nblue+la-1-dk,iray,id) = NLTEspec%S(Nblue+la-1-dk,iray,id) + &
					aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)

				enddo
    
			end do tr_loop

			aatom => NULL()
		end do atom_loop
		
		! Check that chi_loc is correct
        !opacity total = LTE (doesn't change during sub iterations) + NLTE
		!dtau = chi(:) * ds(iray,id)
		NLTEspec%S(:,iray,id) = NLTEspec%S(:,iray,id)/NLTEspec%chi(:,iray,id)
		!!CALL calc_psi_operator(id, icell, iray, chi, dtau, eta)


	RETURN
	END SUBROUTINE calc_total_source_loc
 
	SUBROUTINE NLTE_bound_bound(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)
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
		real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
		logical, intent(in) :: iterate
		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la, Nvspace, dk
		type(AtomType), pointer :: aatom
		integer, parameter 										:: NvspaceMax = 1
		real(kind=dp), dimension(NvspaceMax)						:: Omegav
		real(kind=dp) :: v0
  
  
!-> using dk instead of a full profile ATM
  
!   In principle Nvspace should depend on atom because of atom%vbroad(icell)
!   v_proj in m/s at point icell
		Omegav = 0d0
		Nvspace = 1
! 
		v0 = v_proj(icell,x,y,z,u,v,w) !can be lVoronoi here; for projection
		omegav(1) = v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w)
		dk = nint(1e-3 * omegav(1) / hv)
  
! 		call compute_shift_index(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l)

  
		atom_loop : do nact = 1, atmos%Nactiveatoms
		aatom => atmos%ActiveAtoms(nact)%ptr_atom
		
!Mali
! 			if (iterate) then 
! 				aatom%eta(:,iray,id) = NLTEspec%AtomOpac%jc_nlte(:,icell)
! 			endif

			tr_loop : do kr = 1, aatom%Ntr_line

   	
				kc = aatom%at(kr)%ik !relative index of a transition among continua or lines

				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
				aatom%lines(kc)%dk(iray,id) = dk

				if ((aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0) cycle


				do la=1, aatom%lines(kc)%Nlambda
        
					NLTEspec%AtomOpac%chi(Nblue+la-1-dk,id) = NLTEspec%AtomOpac%chi(Nblue+la-1-dk,id) + &
				hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell) - aatom%lines(kc)%gij*aatom%n(j,icell))
       	 
					NLTEspec%AtomOpac%eta(Nblue+la-1-dk,id)= NLTEspec%AtomOpac%eta(Nblue+la-1-dk,id) + &
				aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)

		
! 					if (iterate) then
! 
! 					!Can do that here because b-b opacities are computed on the fly
!          
! !-> no, dk is the default mode for NLTE.
! !Although it is complicated to allow for variations of the velocity across a cell
!            !keep line profile over directions for this cell (threds) 
!            !for sub iterations and fillGamma()
! !            aatom%lines(kc)%phi_ray(:,iray,id) = aatom%lines(kc)%phi(la,icell)
! !-> MALI         
! 							aatom%eta(Nblue+la-1-dk,iray,id) = aatom%eta(Nblue+la-1-dk,iray,id) + &
!                aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 
!               
! 
! 					end if	
				
				enddo  
			    
    
			end do tr_loop

			aatom => NULL()

		end do atom_loop

	RETURN
	END SUBROUTINE NLTE_bound_bound

END MODULE Opacity

 
