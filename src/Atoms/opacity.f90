MODULE Opacity

 use atmos_type, only : atmos, hydrogen
 use atom_type
 use spectrum_type, only : NLTEspec, initAtomOpac, init_psi_operator, initAtomOpac_nlte
 use constant
 use constantes, only				 : tiny_dp, huge_dp, AU_to_m
 use messages
 !!use voigtfunctions, only 			 : Voigt
 use broad, only 					 : Damping
 use parametres
 use profiles, only : Profile
 !use metal!, only : bound_free_Xsection, Metal_bb_new
 !!use molecular_emission, only : v_proj
 use math, only : locate, integrate_dx, integrate_nu, integrate_x
 !use grid!, only : cross_cell, test_exit_grid
 use mcfost_env, only : dp
 use input, only : ds
 !use stars!, only : intersect_stars
 


 IMPLICIT NONE
 
 real(kind=dp), dimension(:,:,:), allocatable :: chi_loc
 
 CONTAINS
 
  SUBROUTINE calc_psi_operator(id, icell, iray, chi, dtau)
  ! ----------------------------------------------------------- !
   ! Computes Psi and Ieff at the cell icell, for the thread id
   ! in the direction iray, using ds path length of the ray.
  ! ----------------------------------------------------------- !
   integer, intent(in) :: iray, id, icell
   real(kind=dp), dimension(NLTEspec%Nwaves), intent(in) :: chi, dtau
   
   !build from all value
   NLTEspec%Psi(:,iray,id) = (1d0 - dexp(-dtau)) / (chi + tiny_dp)
   
   if (.not.atmos%include_xcoupling) & 
   									NLTEspec%tau(:,iray,id) = dtau


!    if (atmos%include_xcoupling) then
!      NLTEspec%Psi(:,iray,id) = (1d0 - dexp(-dtau)) / (chi + tiny_dp) * dexp(-tau+dtau)
!    else
!      NLTEspec%tau(:,iray,id) = dtau !only allocated for .not.lxcoupling
!      NLTEspec%Psi(:,iray,id) = (1d0 - dexp(-dtau))/(chi+tiny_dp)
!    end if


  RETURN
  END SUBROUTINE calc_psi_operator
  
!  FUNCTION opacity_atom_loc(id, icell, iray) !update eta and psi
!   ! ------------------------------------------------------------------------- !
!   ! ------------------------------------------------------------------------- !  
!   
!   integer, intent(in) :: id, icell, iray
!   real(kind=dp) :: stm
!   integer :: kr, kc, i, j, Nred, Nblue, la, nact
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: opacity_atom_loc
! 
!   
! 
!   opacity_atom_loc(:) = 0.0_dp
!   
!   atom_loop : do nact = 1, atmos%Nactiveatoms
!    aatom => atmos%ActiveAtoms(nact)%ptr_atom
!    
!    	tr_loop : do kr = 1, aatom%Ntr
! 
!    	
!         stm = 1d0 !reset for each transition.
!         kc = aatom%at(kr)%ik !relative index of a transition among continua or lines
! 
!         
!         SELECT CASE (aatom%at(kr)%trtype)
!         
!          CASE ('ATOMIC_LINE')
!          !line = aatom%lines(kc)
!           Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
!           i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
! 
!             
!           if (aatom%lines(kc)%gij*aatom%n(j,icell) > aatom%n(i,icell)) stm = 0d0
! 
!          !opacity total
!           opacity_atom_loc(Nblue:Nred) = opacity_atom_loc(Nblue:Nred) + &
!        		hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell))
!          
! 
!          
!            aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
!               aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * aatom%n(j,icell)
!               
! !            if (atmos%include_xcoupling) then
! !             aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id)
! !             aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + &
! !             	hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell)) * fourPI_hc * aatom%lines(kc)%w_lam(:)
! !             aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + &
! !                 hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell)) * fourPI_hc * aatom%lines(kc)%w_lam(:)
! !            end if
! 
!         
!         CASE ('ATOMIC_CONTINUUM')
! 
! 
!     	  Nred = aatom%continua(kc)%Nred;Nblue = aatom%continua(kc)%Nblue    	
!     	  i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
!     	
!          !I don't know how to handle negative chi due to stimulated emission...
!          !Problem of physics (model) or code ?
!          check_stm : do la=1,aatom%continua(kc)%Nlambda
!           if (aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell) > aatom%n(i,icell)) then 
!             stm = 0d0 
!             exit check_stm
!           endif
!          enddo check_stm   
! 
!          !gather total NLTE
!          opacity_atom_loc(Nblue:Nred) = opacity_atom_loc(Nblue:Nred) + &
!          aatom%continua(kc)%alpha(:) * (aatom%n(i,icell) - stm * aatom%continua(kc)%gij(:,icell)*aatom%n(j,icell))
! 
!     	   aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
!     	   	aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) * aatom%n(j,icell)
! !            if (atmos%include_xcoupling) then
! !             aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) !U
! !             aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + aatom%continua(kc)%alpha(:) * fourPI_h * aatom%continua(kc)%w_lam(:)
! !             aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + aatom%continua(kc)%alpha(:) *fourPI_h * aatom%continua(kc)%w_lam(:)
! !            end if
!    
!     
!     CASE DEFAULT
!      CALL Error("Transition type unknown", aatom%at(kr)%trtype)
!     END SELECT
!     
!    end do tr_loop
! 
!    aatom => NULL()
!   end do atom_loop
! 
! 
!  RETURN
!  END SUBROUTINE opacity_atom_loc
 
 SUBROUTINE init_eta_atom_loc(id, icell, iray) !update eta and psi
  ! ------------------------------------------------------------------------- !
  ! ------------------------------------------------------------------------- !  
  
  integer, intent(in) :: id, icell, iray
  real(kind=dp) :: stm
  integer :: kr, kc, i, j, Nred, Nblue, la, nact
  real(kind=dp), dimension(NLTEspec%Nwaves) :: dtau, chi
  type (AtomType), pointer :: aatom


  !CALL initAtomOpac(id)
  !CALL initAtomOpac_nlte(id) !not need
  CALL init_psi_operator(id,iray) !init also eta
  
  chi(:) = 0.0_dp
  
  atom_loop : do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom
   
   	tr_loop : do kr = 1, aatom%Ntr

   	
        stm = 1d0 !reset for each transition.
        kc = aatom%at(kr)%ik !relative index of a transition among continua or lines

        
        SELECT CASE (aatom%at(kr)%trtype)
        
         CASE ('ATOMIC_LINE')
         !line = aatom%lines(kc)
          Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
          i = aatom%lines(kc)%i;j = aatom%lines(kc)%j

            
          if (aatom%lines(kc)%gij*aatom%n(j,icell) > aatom%n(i,icell)) stm = 0d0

         !opacity total
          chi(Nblue:Nred) = chi(Nblue:Nred) + &
       		hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell))
         

         
           aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
              aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * aatom%n(j,icell)
              
!            if (atmos%include_xcoupling) then
!             aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id)
!             aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + &
!             	hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell)) * fourPI_hc * aatom%lines(kc)%w_lam(:)
!             aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + &
!                 hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_ray(:,iray,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell)) * fourPI_hc * aatom%lines(kc)%w_lam(:)
!            end if

        
        CASE ('ATOMIC_CONTINUUM')


    	  Nred = aatom%continua(kc)%Nred;Nblue = aatom%continua(kc)%Nblue    	
    	  i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
    	
         !I don't know how to handle negative chi due to stimulated emission...
         !Problem of physics (model) or code ?
         check_stm : do la=1,aatom%continua(kc)%Nlambda
          if (aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell) > aatom%n(i,icell)) then 
            stm = 0d0 
            exit check_stm
          endif
         enddo check_stm   

         !gather total NLTE
         chi(Nblue:Nred) = chi(Nblue:Nred) + &
         aatom%continua(kc)%alpha(:) * (aatom%n(i,icell) - stm * aatom%continua(kc)%gij(:,icell)*aatom%n(j,icell))

    	   aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
    	   	aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) * aatom%n(j,icell)
!            if (atmos%include_xcoupling) then
!             aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) !U
!             aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + aatom%continua(kc)%alpha(:) * fourPI_h * aatom%continua(kc)%w_lam(:)
!             aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + aatom%continua(kc)%alpha(:) *fourPI_h * aatom%continua(kc)%w_lam(:)
!            end if
   
    
    CASE DEFAULT
     CALL Error("Transition type unknown", aatom%at(kr)%trtype)
    END SELECT
    
   end do tr_loop

   aatom => NULL()
  end do atom_loop

  dtau = ( chi_loc(:,iray,id) + chi(:) ) * ds(iray,id)
!   if (iray==1) then 
!    write(*,*) ds
!    stop
!   endif
  !recompute Psi
  CALL calc_psi_operator(id, icell, iray, chi, dtau)


 RETURN
 END SUBROUTINE init_eta_atom_loc


 SUBROUTINE calc_J_coherent(id, icell, n_rayons)
  integer, intent(in) :: id, icell, n_rayons
  
  NLTEspec%J(:,icell) = sum(NLTEspec%I(:,1:n_rayons,id),dim=2) / n_rayons
  NLTEspec%Jc(:,icell) = sum(NLTEspec%Ic(:,1:n_rayons,id),dim=2) / n_rayons

 
 RETURN
 END SUBROUTINE calc_J_coherent
 
 
 !beware not test yet on populations
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
  real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
  logical, intent(in) :: iterate
  integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la
!too much time needed for that
!   type(AtomicLine) :: line
!   type(AtomicContinuum) :: cont
  type(AtomType), pointer :: aatom
  real(kind=dp) :: stm
  
  
  atom_loop : do nact = 1, atmos%Nactiveatoms
   aatom => atmos%ActiveAtoms(nact)%ptr_atom
   
   	tr_loop : do kr = 1, aatom%Ntr

   	
        stm = 1d0 !reset for each transition.
        kc = aatom%at(kr)%ik !relative index of a transition among continua or lines

        
        SELECT CASE (aatom%at(kr)%trtype)
        
         CASE ('ATOMIC_LINE')
         !line = aatom%lines(kc)
          Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
          i = aatom%lines(kc)%i;j = aatom%lines(kc)%j

          CALL Profile(aatom%lines(kc),icell,x,y,z,x1,y1,z1,u,v,w,l, id)
          !fill line%phi_loc(:,id)
            
          if (aatom%lines(kc)%gij*aatom%n(j,icell) > aatom%n(i,icell)) stm = 0d0

         !opacity total
         NLTEspec%AtomOpac%chi(Nblue:Nred,id) = NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
       		hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell))
       		
         NLTEspec%AtomOpac%eta(Nblue:Nred,id)= NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
       		aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,id) * aatom%n(j,icell)
      
         if (iterate) then
         
           !keep line profile over directions for this cell (threds) 
           !for sub iterations and fillGamma()
           aatom%lines(kc)%phi_ray(:,iray,id) = aatom%lines(kc)%phi_loc(:,id)
         
           aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
              aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,id) * aatom%n(j,icell)
              
           if (atmos%include_xcoupling) then
            aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,id)
            aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + &
            	hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell)) * fourPI_hc * aatom%lines(kc)%w_lam(:)
            aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + &
                hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,id) * (aatom%n(i,icell) - stm * aatom%lines(kc)%gij*aatom%n(j,icell)) * fourPI_hc * aatom%lines(kc)%w_lam(:)
           end if
         end if	

        
        CASE ('ATOMIC_CONTINUUM')
!write(*,*) "cont loc", loc(cont)
    	  !cont = aatom%continua(kc)
    	  !!if (.not.cont%lcontrib_to_opac) CYCLE tr_loop

    	  Nred = aatom%continua(kc)%Nred;Nblue = aatom%continua(kc)%Nblue    	
    	  i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
    	
         !I don't know how to handle negative chi due to stimulated emission...
         !Problem of physics (model) or code ?
         check_stm : do la=1,aatom%continua(kc)%Nlambda
          if (aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell) > aatom%n(i,icell)) then 
            stm = 0d0 
            exit check_stm
          endif
         enddo check_stm   

         !keep copy of continuum only
    	 NLTEspec%AtomOpac%chic_nlte(Nblue:Nred, id) = NLTEspec%AtomOpac%chic_nlte(Nblue:Nred, id) + &
    	          aatom%continua(kc)%alpha(:) * (aatom%n(i,icell) - stm * aatom%continua(kc)%gij(:,icell)*aatom%n(j,icell))

    	 NLTEspec%AtomOpac%etac_nlte(Nblue:Nred, id) = NLTEspec%AtomOpac%etac_nlte(Nblue:Nred, id) + &
    	 	      aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) * aatom%n(j,icell)

         !gather total NLTE
         NLTEspec%AtomOpac%chi(Nblue:Nred,id) = NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
         aatom%continua(kc)%alpha(:) * (aatom%n(i,icell) - stm * aatom%continua(kc)%gij(:,icell)*aatom%n(j,icell))
       		
         NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
         aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) * aatom%n(j,icell)
    	
    	 if (iterate) then
    	   aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
    	   	aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) * aatom%n(j,icell)
           if (atmos%include_xcoupling) then
            aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + aatom%continua(kc)%alpha(:) * aatom%continua(kc)%twohnu3_c2(:) * aatom%continua(kc)%gij(:,icell) !U
            aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + aatom%continua(kc)%alpha(:) * fourPI_h * aatom%continua(kc)%w_lam(:)
            aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + aatom%continua(kc)%alpha(:) *fourPI_h * aatom%continua(kc)%w_lam(:)
           end if
         end if			
   
    
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


 
  !renorm if not exactly 1 for phi, and integrate the 1/nrays factor.
 !The line and cont weights could eventually be stored on memory to avoid recompute them
 !at each step
!  SUBROUTINE compute_integral_weight(id, icell, iray, n_rayons, x0, y0, z0, u0, v0, w0)
!   integer, intent(in) :: id, n_rayons, icell, iray
!   real(kind=dp), intent(in) :: x0, y0, z0, u0, v0, w0
!   real(kind=dp) :: x1, y1, z1, l
!   real(kind=dp) :: l_c_dum, l_v_dum
!   integer 		   :: n_c_dum
!   type (AtomType) :: atom
!   type (AtomicLine) :: line !too much memory needed for that
!   type (AtomicContinuum) :: cont
!   type (AtomicTransition) :: at
!   real(kind=dp), allocatable, dimension(:) ::weight
!   integer :: nact, k, la
!  
!   CALL cross_cell(x0,y0,z0, u0,v0,w0, icell, &
!        						n_c_dum, x1,y1,z1, n_c_dum, l, l_c_dum, l_v_dum)
!   
!    do nact=1, atmos%Nactiveatoms
!    
!     do k=1, atmos%ActiveAtoms(nact)%ptr_atom%Ntr
!      
!      at = atmos%ActiveAtoms(nact)%ptr_atom%at(k)
!      SELECT CASE (at%trtype)
!       CASE ("ATOMIC_LINE")
!        line = atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)       
!       ! allocate(weight(line%Nlambda)); weight = (NLTEspec%lambda(line%Nblue:line%Nred)-line%lambda0)*CLIGHT/line%lambda0
!        if (iray==1) atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id) = 0d0
!        
!        if (line%voigt) CALL Damping(icell, atmos%ActiveAtoms(nact)%ptr_atom, at%ik, line%adamp)
!        CALL Profile(line,icell,x0,y0,z0,x1,y1,z1,u0,v0,w0,l,line%phi(:,iray,id))
!        !test
! !        if (line%j==3 .and. line%i==2) then
! !         write(*,*) id, iray, line%phi(:,iray,id)
! !        endif
! 
!        !line is not an alias
! !        atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id) = &
! !        		atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id) + &
! !        			Integrate_x(line%Nlambda, weight, line%phi(:,iray,id))
!        atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id) = &
!        		atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id) + &
!        			sum(line_wlam(line)*line%phi(:,iray,id))
!        
! !        write(*,*) iray, atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi
! !        open(1, file="toto", status="old")
! !        do la=1, line%Nlambda
! !        write(1, '(3F)') weight(la), NLTEspec%lambda(line%Nblue+la-1), line%phi(la,iray,id)
! !        enddo
! !        stop
!        
!        if (iray==n_rayons) then
!         !should be close to n_rayons, since the integral is done n_rayons times
!         !write(*,*) atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id)
!         atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id) = 1d0/atmos%ActiveAtoms(nact)%ptr_atom%lines(at%ik)%wphi(id)
!        endif
!        !deallocate(weight)
!        
!         
!       CASE ("ATOMIC_CONTINUUM") !ray integration weight is 1/nrayons
!        cont = atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)
!        
!        if (iray==1) atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)%wmu(id) = 0d0
!        atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)%wmu(id) = &
!         atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)%wmu(id) + 1d0 !should be n_rayons at the end
!        if (iray==n_rayons) then 
!         if (atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)%wmu(id) /= real(n_rayons,kind=dp)) CALL Error("error in cont weight")
!         atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)%wmu(id) = 1d0/atmos%ActiveAtoms(nact)%ptr_atom%continua(at%ik)%wmu(id)
!        endif
!       CASE DEFAULT
!        call error ("unknown, type", at%trtype)
!      END SELECT
!     enddo
!    enddo
!  
!  RETURN
!  END SUBROUTINE compute_integral_weight
 
!  SUBROUTINE NLTEOpacity(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)
!   !
!   !
!   ! chi = Vij * (ni - gij * nj) = ni * Vij - nj * Vji
!   ! eta = twohnu3_c2 * gij * Vij * nj = Uji * nj = 2hnu3/c2 * Vji * nj
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
!   integer, intent(in) :: id, icell, iray
!   real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
!   logical, intent(in) :: iterate
!   integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la
!   type(AtomicLine) :: line
!   type(AtomicContinuum) :: cont
!   type(AtomType), pointer :: aatom
!   real(kind=dp) :: gij, twohnu3_c2, stm, wphi
!   real(kind=dp), dimension(:), allocatable :: Vij, gijk, twohnu3_c2k
!   real(kind=dp), allocatable :: phiZ(:,:), psiZ(:,:)
!   
!   
!   atom_loop : do nact = 1, atmos%Nactiveatoms
!    aatom => atmos%ActiveAtoms(nact)%ptr_atom
!    
!    	tr_loop : do kr = 1, aatom%Ntr
!    	
!    	    !!not useful now. Transition not appearing for the image are removed in aatom%at
!    	    !!and for the NLTEloop all transitions are kept ATM.
!    	    !!if (.not.aatom%at(kr)%lcontrib_to_opac) cycle
!    	
!         stm = 1d0 !reset for each transition.
!         kc = aatom%at(kr)%ik !relative index of a transition among continua or lines
! !write(*,*) "kc=",kc, id, icell, iray, "stm=",stm,"trtype=",aatom%at(kr)%trtype
!         !if (.not.aatom%at(kr)%contrib_to_opac) CYCLE
!         
!         SELECT CASE (aatom%at(kr)%trtype)
!         
!          CASE ('ATOMIC_LINE')
! !write(*,*) 'line loc', loc(line)
!          line = aatom%lines(kc)
!         !!if (.not.line%lcontrib_to_opac) CYCLE tr_loop
!           Nred = line%Nred;Nblue = line%Nblue
!           i = line%i;j = line%j
!     
!     if ((aatom%n(j,icell) < tiny_dp).or.(aatom%n(i,icell) < tiny_dp)) then !no transition    
!         !probably False sharing
!         if (aatom%n(j,icell)==0d0 .or. aatom%n(i,icell)==0d0) then
!          write(*,*) " ***************************** "
!      	 CALL WARNING("False sharing")
! !          write(*,*) icell, iray, id, aatom%ID, aatom%Nlevel, kc, shape(aatom%n)
! !          write(*,*) i, line%i, j, line%j
! !          write(*,*) aatom%n(:,icell)
! !          write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%n(line%i,icell), aatom%n(line%j,icell)
! !          write(*,*) "1", aatom%n(1,icell), "2", aatom%n(2,icell), "3", aatom%n(3,icell),"4", aatom%n(4,icell)
!          !stop
!          write(*,*) " ***************************** "
!      	 cycle tr_loop !go to next transitions without computing opac
!      	else if (aatom%n(j,icell) < 0 .or. aatom%n(i,icell) < 0) then
!      	 CALL WARNING("False sharing ? negative pops")
!      	 cycle tr_loop
!         else
!          write(*,*) " ***************************** "
!      	 CALL WARNING("too small line populations")
!          write(*,*) icell, iray, id, aatom%ID, aatom%Nlevel, kc, shape(aatom%n)
!          write(*,*) i, line%i, j, line%j
!          write(*,*) aatom%n(:,icell)
!          write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%n(line%i,icell), aatom%n(line%j,icell)
!          write(*,*) "1", aatom%n(1,icell), "2", aatom%n(2,icell), "3", aatom%n(3,icell),"4", aatom%n(4,icell)
!          write(*,*) " ***************************** "
!      	endif
!      	!treating problem in LTE ?
! !          aatom%n(j,icell) = aatom%nstar(j,icell)
! !      	 aatom%n(i,icell) = aatom%nstar(i,icell)
! !          aatom%n(j,icell) = max(tiny_dp, aatom%n(j,icell))
! !      	 aatom%n(i,icell) = max(tiny_dp, aatom%n(i,icell))
!     end if 
! 
!         gij = line%Bji / line%Bij
! 
!         twohnu3_c2 = line%Aji / line%Bji
!         if (line%voigt)  CALL Damping(icell, aatom, kc, line%adamp)
!         !if (line%adamp>5.) write(*,*) " large damping for line", line%j, line%i, line%atom%ID, line%adamp
!     
!         !allocate(phi(line%Nlambda),Vij(line%Nlambda))
!         allocate(Vij(line%Nlambda))
!     
!         if (PRT_SOLUTION=="FULL_STOKES") allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
!         !phiZ and psiZ are used only if Zeeman polarisation, which means we care only if
!         !they are allocated in this case.
!         CALL Profile(line, icell,x,y,z,x1,y1,z1,u,v,w,l,Vij,phiZ,psiZ)!, phi, phiZ, psiZ)
! 
!         if (iterate) aatom%lines(kc)%phi(:,iray,id) = Vij(:)!phi(:)
! !attention, line ne pointe pas vers Vij, car line est une instance de atom%lines(kc) pas un pointeur
! !        if (line%j==3 .and. line%i==2) then
! !         write(*,*) Vij
! !        stop
! !        endif        
!         
!          Vij(:) = hc_fourPI * line%Bij * Vij(:)!phi(:) !normalized in Profile()
!                                                              ! / (SQRTPI * VBROAD_atom(icell,aatom)) 
!             !write(*,*) 'bb', icell, id, iray, maxval(Vij), minval(Vij)
!             
!          !I don't know how to handle negative chi due to stimulated emission...
!          !Problem of physics (model) or code ?
!         if (gij*aatom%n(j,icell) > aatom%n(i,icell)) stm = 0d0
! 
!          !opacity total
!          NLTEspec%AtomOpac%chi(Nblue:Nred,id) = NLTEspec%AtomOpac%chi(Nblue:Nred,id) + &
!        		Vij(:) * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell))
!        		
!          NLTEspec%AtomOpac%eta(Nblue:Nred,id)= NLTEspec%AtomOpac%eta(Nblue:Nred,id) + &
!        		twohnu3_c2 * gij * Vij(:) * aatom%n(j,icell)
!       
!         if (iterate) then
!            aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + &
!               twohnu3_c2 * gij * Vij(:) * aatom%n(j,icell)
!               
!            if (atmos%include_xcoupling) then
!             aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + twohnu3_c2 * gij * Vij(:)
!             aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + &
!             	Vij(:) * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell)) * fourPI_hc * line%wphi(id) * line_wlam(line) 
!             aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + &
!             	Vij(:) * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell)) * fourPI_hc * line%wphi(id) * line_wlam(line) 
!            end if
!         end if	
! 
!     
!         if (line%polarizable .and. PRT_SOLUTION == "FULL_STOKES") then
!          write(*,*) "Beware, NLTE part of Zeeman opac not set to 0 between iteration!"
!          do nk = 1, 3
!           !magneto-optical
!           NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) + &
!            hc_fourPI * line%Bij * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell)) * psiZ(nk,:)
!           !dichroism
!           NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) + &
!            hc_fourPI * line%Bij * (aatom%n(i,icell) - stm * gij*aatom%n(j,icell)) * psiZ(nk,:)
!           !emissivity
!           NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) + &
!           twohnu3_c2 * gij * hc_fourPI * line%Bij * aatom%n(j,icell) * phiZ(nk,:)
!          end do 
!         end if
!      
!        deallocate(Vij)!, phi)
!        if (PRT_SOLUTION=="FULL_STOKES") deallocate(phiZ, psiZ)
!         
!         CASE ('ATOMIC_CONTINUUM')
! !write(*,*) "cont loc", loc(cont)
!     	  cont = aatom%continua(kc)
!     	  !!if (.not.cont%lcontrib_to_opac) CYCLE tr_loop
! 
!     	  Nred = cont%Nred;Nblue = cont%Nblue    	
!     	  i = cont%i;j = cont%j
!         
!     	  if (aatom%n(j,icell) < tiny_dp .or. aatom%n(i,icell) < tiny_dp) then !test  
!            if (aatom%n(j,icell)==0d0 .or. aatom%n(i,icell)==0d0) then
!             write(*,*) " ***************************** "
!      	    CALL WARNING("False sharing")
! !             write(*,*) icell, iray, id, aatom%ID, aatom%Nlevel, kc, shape(aatom%n)
! !             write(*,*) i, cont%i, j, cont%j
! !             write(*,*) aatom%n(:,icell)
! !             write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%n(cont%i,icell), aatom%n(cont%j,icell)
! !             write(*,*) "1", aatom%n(1,icell), "2", aatom%n(2,icell), "3", aatom%n(3,icell),"4", aatom%n(4,icell)
!             write(*,*) " ***************************** "
!      	    cycle tr_loop !go to next transition
!      	   else if (aatom%n(j,icell) < 0 .or. aatom%n(i,icell) < 0) then
!      	    CALL WARNING("False sharing ? negative pops")
!      	    cycle tr_loop
!            else
!            write(*,*) " ***************************** "
!      	   CALL WARNING("too small cont populations")
!            write(*,*) icell, iray, id, aatom%ID, aatom%Nlevel, kc, shape(aatom%n)
!            write(*,*) i, cont%i, j, cont%j
!            write(*,*) aatom%n(:,icell)
!            write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%n(cont%i,icell), aatom%n(cont%j,icell)
!            write(*,*) "1", aatom%n(1,icell), "2", aatom%n(2,icell), "3", aatom%n(3,icell),"4", aatom%n(4,icell)
!            write(*,*) " ***************************** "
!      	   endif 	 
!      	  !treating pb in LTE ?
! !      	  aatom%n(j,icell) = aatom%nstar(j,icell)
! !      	  aatom%n(i,icell) = aatom%nstar(i,icell)
! !      	  aatom%n(j,icell) = max(tiny_dp, aatom%n(j,icell))
! !      	  aatom%n(i,icell) = max(tiny_dp, aatom%n(i,icell))
!     	  end if !end test
!     	  
!       	  allocate(gijk(cont%Nlambda)); gijk(:) = 0d0
!       	  if (aatom%nstar(j,icell)>0) &
!           	gijk(:) = aatom%nstar(i, icell)/aatom%nstar(j,icell) * dexp(-hc_k / (NLTEspec%lambda(Nblue:Nred) * atmos%T(icell)))
!       
!          !allocate Vij, to avoid computing bound_free_Xsection(cont) 3 times for a continuum
! 	     allocate(Vij(cont%Nlambda), twohnu3_c2k(cont%Nlambda))
! 	     
!     	 Vij(:) = bound_free_Xsection(cont) 	
!     	       !write(*,*) 'bf', icell, id, iray, maxval(Vij), minval(Vij)
! 
!          twohnu3_c2k(:) = twohc / NLTEspec%lambda(Nblue:Nred)**(3d0) * &
!               Vij(:) * gijk(:) * aatom%n(j,icell) !eta
!               
!          
!          !I don't know how to handle negative chi due to stimulated emission...
!          !Problem of physics (model) or code ?
!          check_stm : do la=1,cont%Nlambda
!           if (gijk(la)*aatom%n(j,icell) > aatom%n(i,icell)) then 
!             stm = 0d0 
!             exit check_stm
!           endif
!          enddo check_stm   
!          Vij(:) = Vij(:) * (aatom%n(i,icell) - stm * gijk(:)*aatom%n(j,icell)) !chi
! 
!          !keep copy of continuum only
!     	 NLTEspec%AtomOpac%chic_nlte(Nblue:Nred, id) = NLTEspec%AtomOpac%chic_nlte(Nblue:Nred, id) + Vij(:)
!     	 NLTEspec%AtomOpac%etac_nlte(Nblue:Nred, id) = NLTEspec%AtomOpac%etac_nlte(Nblue:Nred, id) + twohnu3_c2k(:)
! 
!          !gather total NLTE
!          NLTEspec%AtomOpac%chi(Nblue:Nred,id) = NLTEspec%AtomOpac%chi(Nblue:Nred,id) + Vij(:)
!        		
!          NLTEspec%AtomOpac%eta(Nblue:Nred,id) = NLTEspec%AtomOpac%eta(Nblue:Nred,id) + twohnu3_c2k(:)
!     	
!     	 if (iterate) then
!     	   aatom%eta(Nblue:Nred,iray,id) = aatom%eta(Nblue:Nred,iray,id) + twohnu3_c2k(:)
!            if (atmos%include_xcoupling) then
!             aatom%Uji_down(j,Nblue:Nred, iray, id) = aatom%Uji_down(j,Nblue:Nred, iray, id) + twohnu3_c2k(:)/aatom%n(j,icell) !U
!             aatom%chi_up(i,Nblue:Nred, iray, id) = aatom%chi_up(i,Nblue:Nred, iray, id) + Vij(:) * fourPI_h *cont_wlam(cont)*cont%wmu(id)
!             aatom%chi_down(j,Nblue:Nred,iray,id) = aatom%chi_down(j,Nblue:Nred,iray,id) + Vij(:) *fourPI_h *cont_wlam(cont)*cont%wmu(id)!chi
!            end if
!          end if			
!    
!         deallocate(Vij, gijk, twohnu3_c2k)
!     
!     CASE DEFAULT
!      CALL Error("Transition type unknown", aatom%at(kr)%trtype)
!     END SELECT
!     
!    end do tr_loop
!      !write(*,*) icell, id, "eta=", maxval(aatom%eta(:,iray,id)), minval(aatom%eta(:,iray,id))
! 
!    aatom => NULL()
!    !stop
!   end do atom_loop
! 
!  RETURN
!  END SUBROUTINE NLTEOpacity
 
