! This module computes the metals bound-free opacitiy according
! to Kramer's law if Hydrogenic or using detailed photoionisation
! cross-section read in the atomic model.
!
! Note that Hydrogen bound-free is computed in the hydrogen.f90 module
!
! The module also computes the bound-bound opacity of metals (including
! Hydrogen this time) treated as PASSIVE atoms (either LTE populations or
! NLTE pops for a previous run) and sums all PASSIVE opacities in NLTEspec%chi and %eta
!
! chi in m^-1, eta in J/m3/s/Hz/sr
MODULE metal

 use atmos_type, only                : atmos, Hydrogen, Helium, B_project, vbroad_atom, ntotal_atom
 use constant
 use math, only                      : bezier3_interp, cent_deriv, any_nan_infinity_vector, linear_1D
 use atom_type
 use spectrum_type, only			 : NLTEspec, initAtomOpac
 use hydrogen_opacities
 use voigtfunctions, only 			 : Voigt
 use Profiles, only					 : Profile !, line_profiles
 use broad, only 					 : Damping
 use thomson_scattering
 use Rayleigh_scattering
 use Planck

 ! MCFOST's original
 use mcfost_env, only : dp
 use molecular_emission, only		 : v_proj
 use parametres
 use input
 use constantes, only				 : tiny_dp, huge_dp

 IMPLICIT NONE

 CONTAINS
 
 !futur deprecation
 FUNCTION bound_free_Xsection(cont) result(alpha)
  type (AtomicContinuum) :: cont
  real(kind=dp) :: alpha(cont%Nlambda) 
  									 
   if (cont%Hydrogenic) then !Kramer's formula with quantum mechanical correction
     alpha = H_bf_Xsection(cont, NLTEspec%lambda(cont%Nblue:cont%Nred))
   else !interpolation of the read Cross-section
    !alpha = interp_dp(cont%alpha, cont%lambda, NLTEspec%lambda(cont%Nblue:cont%Nred))
    CALL bezier3_interp(size(cont%alpha), cont%lambda, cont%alpha, & !read values
     	cont%Nlambda, NLTEspec%lambda(cont%Nblue:cont%Nred), alpha) !interpolation grid

   endif
 
 RETURN
 END FUNCTION bound_free_Xsection
 
  
!   SUBROUTINE compute_photoionisation_xsections()
!    ! Computes and stores the photo-ionisation cross sections + cont%gij
!    ! for each contiuum of each atom
!     type(AtomType), pointer :: atom
!     integer :: n, kc, k, alloc_status
!     type(AtomicContinuum) :: cont
! !     real(kind=dp), allocatable :: dumm_alpha(:)
!     
!     do n=1,atmos%Natom
!      atom => atmos%Atoms(n)%ptr_atom
!      
!      !only over continua counted in opacities
!      do kr=atom%Ntr_line+1, atom%Ntr!1, atom%Ncont !loop over continua
!       kc = atom%at(kr)%ik 
!       cont = atom%continua(kc)
! 
!       if (cont%Hydrogenic) then !Kramer's formula with quantum mechanical correction
!         allocate(atom%continua(kc)%alpha(cont%Nlambda), stat=alloc_status)
!         if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")
!         
!         atom%continua(kc)%alpha(:) = 0d0
!         !computes only where it is not zero
!         atom%continua(kc)%alpha(cont%Nblue:cont%Nred) = & 
!         	H_bf_Xsection(cont, NLTEspec%lambda(cont%Nblue:cont%Nred))
!         	
!       else !interpolation of the read Cross-section
! !         allocate(dumm_alpha(size(cont%lambda)),stat=alloc_status)
! !         if (alloc_status > 0) CALL ERROR("Allocation error dumm_alpha")
! !         dumm_alpha(:) = cont%alpha(:)
! !       cont is not an alias but a copy of contiua(kc), so cont%alpha not deallocated
! 		deallocate(atom%continua(kc)%alpha)
!         allocate(atom%continua(kc)%alpha(cont%Nlambda), stat=alloc_status)
!         if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")
!         !alpha = interp_dp(cont%alpha, cont%lambda, NLTEspec%lambda(cont%Nblue:cont%Nred))
!         
!         CALL bezier3_interp(size(cont%alpha), cont%lambda, cont%alpha, & !read values
!      	cont%Nlambda, NLTEspec%lambda(cont%Nblue:cont%Nred), atom%continua(kc)%alpha) !interpolation grid
!      	!the real Nlambda from Nblue to Nred
!      	!deallocate(dumm_alpha)
!         deallocate(atom%continua(kc)%lambda) !not used anymore, we kept in memory the interpolated data
!       endif
!      enddo
!      atom => NULL()
!     enddo
!   
!   
!   RETURN
!   END SUBROUTINE compute_photoionisation_xsections

  SUBROUTINE dealloc_atom_quantities
  
    type(AtomType), pointer :: atom
    integer :: n, k, kc, j, alloc_status
    type(AtomicContinuum) :: cont
    type(AtomicLine) :: line

    
    do n=1,atmos%Natom
     atom => atmos%Atoms(n)%ptr_atom
        
     do k = 1, atom%Ntr   
     
       kc = atom%at(k)%ik 
        
       SELECT CASE (atom%at(k)%trtype)
        
        CASE ('ATOMIC_LINE')
          line = atom%lines(kc)
          
          !Damping value
          !no need for it if we do keep the profile on cell
          if (line%voigt) then
           deallocate(atom%lines(kc)%a) !futur deprec depending of the choice of Voigt
           deallocate(atom%lines(kc)%aeff,atom%lines(kc)%r,atom%lines(kc)%r1)
          endif
          !line profile for v = 0
          deallocate(atom%lines(kc)%phi)
          deallocate(atom%lines(kc)%phi_loc)
          
          deallocate(atom%lines(kc)%u)

        
        CASE ('ATOMIC_CONTINUUM')
     
         !do kc=1, atom%Ncont !loop over continua
          cont = atom%continua(kc)
      
          !gij for bound-free transitions
          deallocate(atom%continua(kc)%gij)
          deallocate(atom%continua(kc)%twohnu3_c2)
          deallocate(atom%continua(kc)%alpha)
		  deallocate(atom%continua(kc)%Vji, atom%continua(kc)%Jnu)

 
         !enddo
       CASE DEFAULT
        CALL Error("Transition type unknown", atom%at(k)%trtype)
       END SELECT
     enddo
     atom => NULL()
    enddo  
  
  RETURN
  END SUBROUTINE dealloc_atom_quantities
  
  SUBROUTINE alloc_atom_quantities
  
    type(AtomType), pointer :: atom
    integer :: n, k, kc, j, alloc_status, la
!     real(kind=dp), allocatable :: dumm_alpha(:)

    
    do n=1,atmos%Natom
     atom => atmos%Atoms(n)%ptr_atom
        
     do k = 1, atom%Ntr   
     
       kc = atom%at(k)%ik 
        
       SELECT CASE (atom%at(k)%trtype)
        
        CASE ('ATOMIC_LINE')          
          !Tex
          allocate(atom%lines(kc)%Tex(atmos%Nspace), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%Tex")
          !alloc at LTE
          atom%lines(kc)%Tex(:) = atmos%T(:) !0.0_dp          
          
          !some memory can be saved depending of the choice of the profile
          ! I have to do it for large simu
          
          !Damping value
          !no need for it if we do keep the profile on cell
          !except if updated, but if we chose thomson method it is need along to line%aeff
          if (atom%lines(kc)%voigt) then
           allocate(atom%lines(kc)%a(atmos%Nspace), stat=alloc_status)
           if (alloc_status > 0) CALL ERROR("Allocation error line%adamp")
           atom%lines(kc)%a(:) = 0d0
           allocate(atom%lines(kc)%aeff(atmos%Nspace), stat=alloc_status)
           if (alloc_status > 0) CALL ERROR("Allocation error line%eff")
           atom%lines(kc)%aeff(:) = 0d0
           allocate(atom%lines(kc)%r(atmos%Nspace), stat=alloc_status)
           if (alloc_status > 0) CALL ERROR("Allocation error line%r")
           atom%lines(kc)%r(:) = 0d0
           allocate(atom%lines(kc)%r1(atmos%Nspace), stat=alloc_status)
           if (alloc_status > 0) CALL ERROR("Allocation error line%r1")
           atom%lines(kc)%r1(:) = 0d0
          endif
          !line profile for v = 0
          !no need to keep it if not interp, no shift. No need also if Thomson
          allocate(atom%lines(kc)%phi(atom%lines(kc)%Nlambda, atmos%Nspace), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%phi")
          atom%lines(kc)%phi(:,:) = 0d0
          !-> use in NLTEloop only now
          !allocate(atom%lines(kc)%phi_ray(line%Nlambda, atmos%Nrays, NLTEspec%NPROC), stat=alloc_status)
          !if (alloc_status > 0) CALL ERROR("Allocation error line%phi_ray")
          allocate(atom%lines(kc)%phi_loc(atom%lines(kc)%Nlambda, NLTEspec%NPROC), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%phi_loc")
          
          
          allocate(atom%lines(kc)%u(atom%lines(kc)%Nlambda), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%u")
          !so that argument of profile is line%u / vbroad(icell)
          atom%lines(kc)%u(:) = (NLTEspec%lambda(atom%lines(kc)%Nblue:atom%lines(kc)%Nred)-atom%lines(kc)%lambda0) * &
           CLIGHT / atom%lines(kc)%lambda0 !m/s
           
          allocate(atom%lines(kc)%Jbar(NLTEspec%NPROC), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%Jbar")
          atom%lines(kc)%Jbar(:) = 0d0
        
        CASE ('ATOMIC_CONTINUUM')
     
         !do kc=1, atom%Ncont !loop over continua
          
          !Tex
          allocate(atom%continua(kc)%Tex(atmos%Nspace), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%Tex")
          !first LTE
          atom%continua(kc)%Tex(:) = atmos%T(:)
      
          !gij for bound-free transitions
          allocate(atom%continua(kc)%gij(atom%continua(kc)%Nlambda, atmos%Nspace), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%gij")
          atom%continua(kc)%gij(:,:) = 0d0
!either Vji or Jnu will be remove eventually 
          !gij * Vij(=alpha) = Vji
          allocate(atom%continua(kc)%Vji(atom%continua(kc)%Nlambda, NLTEspec%NPROC), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%Vji")
          atom%continua(kc)%Vji(:,:) = 0d0
          
          allocate(atom%continua(kc)%Jnu(atom%continua(kc)%Nlambda, NLTEspec%NPROC), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%Jnu")
          atom%continua(kc)%Jnu(:,:) = 0d0

          !twohnu3_c2 for continuum waves
          allocate(atom%continua(kc)%twohnu3_c2(atom%continua(kc)%Nlambda),stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%twohnu3_c2")
          atom%continua(kc)%twohnu3_c2(:) = twohc / NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred)**3.
          
          !special for bound-free cross-sections who do not depend on cell
    		if (atom%continua(kc)%Hydrogenic) then !Kramer's formula with quantum mechanical correction
        		allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
        		if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")

        		atom%continua(kc)%alpha(:) = 0d0
        		!computes only where it is not zero
        		atom%continua(kc)%alpha(:) = H_bf_Xsection(atom%continua(kc), NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))

      		else !interpolation of the read Cross-section
!       cont is not an alias but a copy of contiua(kc), so cont%alpha not deallocated
        		allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
        		if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")
         
         atom%continua(kc)%alpha(:) = linear_1D(size(atom%continua(kc)%lambda_file),atom%continua(kc)%lambda_file,&
            atom%continua(kc)%alpha_file, atom%continua(kc)%Nlambda,NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))
!         		CALL bezier3_interp(size(atom%continua(kc)%lambda_file), atom%continua(kc)%lambda_file, atom%continua(kc)%alpha_file, & !read values
!      	atom%continua(kc)%Nlambda, NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), atom%continua(kc)%alpha) !interpolation grid

      		endif
 
         !enddo
       CASE DEFAULT
        CALL Error("Transition type unknown", atom%at(k)%trtype)
       END SELECT
     enddo
     atom => NULL()
    enddo  

  RETURN
  END SUBROUTINE alloc_atom_quantities
  
  SUBROUTINE compute_atom_quantities(icell) !those quantities are to be updated, most likely during iterations
   ! Computes and stores the photo-ionisation cross sections + cont%gij
   ! for each contiuum of each atom
    integer, intent(in) :: icell
    type(AtomType), pointer :: atom
    integer :: n, k, kc, alloc_status, i, j, Nblue, Nred
    real(kind=dp) :: vbroad, aL, cte, cte2, adamp
    
    
    do n=1,atmos%Natom
     atom => atmos%Atoms(n)%ptr_atom
     do k = 1, atom%Ntr   
     
       kc = atom%at(k)%ik 
        
       SELECT CASE (atom%at(k)%trtype)
        
        CASE ('ATOMIC_LINE')
         vbroad = atom%vbroad(icell)
         
         i=atom%lines(kc)%i; j=atom%lines(kc)%j

         
         if (atom%lines(kc)%voigt) then
          CALL Damping(icell, atom, kc, adamp)
          !no need to keep it if we keep the profile for shift or interp
          atom%lines(kc)%a(icell) = adamp
          !line%u(:) = atom%lines(kc)%u / atom%vbroad(icell)
          !                             if full profile to be interpolated or shifter.
          !                             if we use another approx for Voigt, it is computed on the fly, not shifted nor interp
          atom%lines(kc)%phi(:,icell) = Voigt(atom%lines(kc)%Nlambda, adamp, atom%lines(kc)%u(:)/vbroad)
          !write(*,*) "d = ", line%adamp
          !write(*,*) maxval(atom%lines(kc)%phi(:,icell)), minval(atom%lines(kc)%phi(:,icell))
          !Thomson method effective damping for all cells
          aL = adamp * vbroad !(m/s), adamp in doppler units
          atom%lines(kc)%aeff(icell) = (vbroad**5. + 2.69269*vbroad**4. * aL + 2.42843*vbroad**3. * aL**2. + &
          		4.47163*vbroad**2.*aL**3. + 0.07842*vbroad*aL**4. + aL**5.)**(0.2) !1/5
          		
          !!there should have simplification here
          !!r = line%atom%vbroad(icell)*line%a(icell)/line%aeff(icell)
          !!eta = 1.36603*r - 0.47719*r*r + 0.11116*r*r*r
          cte = vbroad*adamp/atom%lines(kc)%aeff(icell)
          cte2 = 1.36603*cte - 0.47719*cte*cte + 0.11116*cte*cte*cte
          atom%lines(kc)%r(icell) = cte2/pi
          !for the lorentzian it is eta/pi and for the gaussian it is (1-eta)/sqrtpi/aeff
          atom%lines(kc)%r1(icell) = (1. - cte2)/sqrtpi/atom%lines(kc)%aeff(icell)
         else !Gaussian
          !no need if no interp nor shift or if we use Thomson (but keep adamp for thomson)
          atom%lines(kc)%phi(:,icell) = dexp(-(atom%lines(kc)%u(:)/atom%vbroad(icell))**2)
         endif
         !futur test to see if we keep it or not depending on the method of calculation of profiles
         atom%lines(kc)%phi(:,icell) = atom%lines(kc)%phi(:,icell) / (SQRTPI * atom%vbroad(icell))

         if ((atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij) <= 0 ) then
          !write(*,*) i, j, atom%n(i,icell), atom%n(j, icell), atom%lines(kc)%gij
          !write(*,*) atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij
          write(*,*) atom%ID, " nuij (10^15 Hz)", 1d-6 * CLIGHT / atom%lines(kc)%lambda0, " icell=", icell
          write(*,*) " lambdaij (nm)", atom%lines(kc)%lambda0 
          CALL warning ("background line: ni < njgij")
          atom%lines(kc)%neg_opac = .true. !should be an icell function ???
          									! to remove only on specific position
          									
         endif
    
        CASE ("ATOMIC_CONTINUUM")
        
         i=atom%continua(kc)%i; j=atom%continua(kc)%j
         Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
   
         if (atom%nstar(j,icell) < tiny_dp .or. atom%nstar(i,icell) < tiny_dp) then
          atom%continua(kc)%gij(:,icell) = 0d0
         else
          atom%continua(kc)%gij(:, icell) = atom%nstar(i,icell)/atom%nstar(j,icell) *  &
     			dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell))
     			
         endif
! write(*,*) "test"        
! write(*,*) icell, atmos%T(icell), atom%nstar(j, icell), atom%n(j, icell), atom%nstar(j, icell)/atom%n(j,icell) , minval(dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell))), maxval(dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell)))
! write(*,*) atom%n(i,icell) - atom%n(j, icell)*minval(atom%continua(kc)%gij(:, icell))
! write(*,*) atom%n(i,icell), atom%n(j, icell)*atom%nstar(i, icell)/atom%nstar(j,icell)

         if (atom%n(i,icell) - atom%n(j, icell)*minval(atom%continua(kc)%gij(:, icell)) <= 0 ) then
         !if (1d0 - minval(dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell))) <= 0 ) then
          !write(*,*) i, j, atom%n(i,icell), atom%n(j, icell), minval(atom%continua(kc)%gij(:, icell)), atom%nstar(i,icell), atom%nstar(j,icell)
          !write(*,*) atom%n(i,icell) - atom%n(j, icell)*minval(atom%continua(kc)%gij(:, icell))
          write(*,*) atom%ID, " nu+ (10^15 Hz)", 1d-6 * CLIGHT / atom%continua(kc)%lambda0, " icell=", icell
          write(*,*) " lambda+ (nm)", atom%continua(kc)%lambda0
          CALL warning ("background cont: ni < njgij")
          atom%continua(kc)%neg_opac = .true.
         endif
! write(*,*) "endtest"        
        
        CASE DEFAULT
         CALL Error("Transition type unknown", atom%at(k)%trtype)
       END SELECT
     enddo
     atom => NULL()
    enddo  
  
  
  RETURN
  END SUBROUTINE compute_atom_quantities
  
 !for each cell. Done and kept in memory, recomputed if LTE populations and electrons changed
 SUBROUTINE Metal_bf_new(icell)
 !cross-section in cm2 per particle is given by Kramers’ formula
  !with n the principal quantum number of the level i from
 !which the atom or ion is ionized, Z the ion charge, ν in Hz and gbf the dimensionless
 !Gaunt factor, a quantummechanical correction factor of order unity.
 !The Kramers cross-section decays ∼ ν−3 above the threshold (“edge”) frequency ν0,
 !being zero below it because the threshold energy is the required minimum. Think the
  ! inverse in terms of wavelengths
  integer, intent(in)							            :: icell
  integer                                                   :: m, kr, kc, i, j, Nblue, Nred, la
  type (AtomType), pointer                                  :: atom
  
  ! Go throught all bound-free transitions of each PASSIVE
  ! metal and add the opacity and emissivity if lambda
  ! is lower (greater) than the wavelength threshold lambdaEdge
  ! (the frequency threshold) and if greater (lower) than
  ! wavelength min (frequency max). See Hydrogen b-f for more
  ! informations, and Hubeny & Mihalas chap. 7

  do m=1,atmos%Npassiveatoms
  ! run over all passive atoms

   atom => atmos%PassiveAtoms(m)%ptr_atom!atmos%Atoms(m)
   !if (metal%ID == "H ") CYCLE !H cont is treated in Hydrogen_bf()
   
    do kc=atom%Ntr_line+1,atom%Ntr
     kr = atom%at(kc)%ik 
    !do kr=1,atomx%Ncont
    
     i = atom%continua(kr)%i
     j = atom%continua(kr)%j 
     Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred
     

     !inner wavelength loop to allow to handle negative continuum opacity
     do la=1, atom%continua(kr)%Nlambda
        if (atom%n(i,icell) - atom%n(j,icell)*atom%continua(kr)%gij(la,icell) < 0) cycle
        
		!if (atom%n(i,icell) - atom%n(j,icell)*atom%continua(kr)%gij(la,icell) >= 0.) then
    		 NLTEspec%AtomOpac%Kc(Nblue+la-1,icell,1) = NLTEspec%AtomOpac%Kc(Nblue+la-1,icell,1) + &
       	atom%continua(kr)%alpha(la) * (atom%n(i,icell)-atom%continua(kr)%gij(la,icell)*atom%n(j,icell))
        !endif

     	NLTEspec%AtomOpac%jc(Nblue+la-1,icell) = NLTEspec%AtomOpac%jc(Nblue+la-1,icell) + &
       	atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * atom%continua(kr)%gij(la,icell) * atom%n(j,icell)
     
     enddo

!      NLTEspec%AtomOpac%Kc(Nblue:Nred,icell,1) = NLTEspec%AtomOpac%Kc(Nblue:Nred,icell,1) + &
!        	kappa_neg * atom%continua(kr)%alpha(:) * (atom%n(i,icell)-atom%continua(kr)%gij(:,icell)*atom%n(j,icell))
! 
!      NLTEspec%AtomOpac%jc(Nblue:Nred,icell) = NLTEspec%AtomOpac%jc(Nblue:Nred,icell) + &
!        	atom%continua(kr)%alpha(:) * atom%continua(kr)%twohnu3_c2(:) * atom%continua(kr)%gij(:,icell) * atom%n(j,icell)

    end do ! loop over Ncont
    atom => NULL()
  end do !loop over metals
  
 RETURN
 END SUBROUTINE Metal_bf_new

!->futur deprec
 !I try to compute H_bf here also
 SUBROUTINE Metal_bf(id, icell)
 !cross-section in cm2 per particle is given by Kramers’ formula
  !with n the principal quantum number of the level i from
 !which the atom or ion is ionized, Z the ion charge, ν in Hz and gbf the dimensionless
 !Gaunt factor, a quantummechanical correction factor of order unity.
 !The Kramers cross-section decays ∼ ν−3 above the threshold (“edge”) frequency ν0,
 !being zero below it because the threshold energy is the required minimum. Think the
  ! inverse in terms of wavelengths
  integer, intent(in)							            :: icell, id
  logical 										            :: obtained_n
  integer                                                   :: m, kr, kc, i, j, Z, nc, Nblue, Nred
  type (AtomType)                                           :: metal
  type (AtomicContinuum)                                    :: continuum
  real(kind=dp)                                          :: lambdaEdge
  real(kind=dp), dimension(:), allocatable               :: twohnu3_c2, gijk, Vij
   !obtained_n = .false. !true if the routine to read principal quantum number is fine
   
   !if (atmos%Npassiveatoms==1 .and. atmos%PassiveAtoms(1)%ptr_atom%ID == "H ") RETURN


  ! Go throught all bound-free transitions of each PASSIVE
  ! metal and add the opacity and emissivity if lambda
  ! is lower (greater) than the wavelength threshold lambdaEdge
  ! (the frequency threshold) and if greater (lower) than
  ! wavelength min (frequency max). See Hydrogen b-f for more
  ! informations, and Hubeny & Mihalas chap. 7

  do m=1,atmos%Npassiveatoms
  ! run over all passive atoms
   metal = atmos%PassiveAtoms(m)%ptr_atom!atmos%Atoms(m)
   !if (metal%ID == "H ") CYCLE !H cont is treated in Hydrogen_bf()
    do kc=metal%Ntr_line+1,metal%Ntr
     kr = metal%at(kc)%ik 
    !do kr=1,metal%Ncont
     continuum = metal%continua(kr)
     i = continuum%i
     j = continuum%j !+1 wrt C indexing
     Nblue = continuum%Nblue; Nred = continuum%Nred
     !!if (.not.metal%at(metal%Nline+kr)%lcontrib_to_opac) CYCLE !avoid continua not defined on the grid
     
     allocate(twohnu3_c2(continuum%Nlambda), gijk(continuum%Nlambda), Vij(continuum%Nlambda))
     !hc_kla = hc_k/NLTEspec%lambda(Nblue:Nred) !factor 1/NM_TO_M in hc_k
     twohnu3_c2 = twohc / NLTEspec%lambda(Nblue:Nred)**3
     Vij(:) = bound_free_Xsection(continuum)

     lambdaEdge = continuum%lambda0! or ionisation wavelength or wavelength
               ! associated to the minimal frquency needed
               ! to unbound an electron
               
               
    !futur ndeprecation, if the photoionized level is too low, just neglect emission
   if (metal%n(j,icell) < tiny_dp) then !==0, but < tiny_dp is okay also
        if ((metal%n(j,icell) < 0) .or. (metal%n(i,icell) < 0)) then
        !-> Futur deprecation of this test, psoitivity is tested when populations are set
         write(*,*) "(Metal_bf) Warning at icell, negative pops", icell," T(K)=", atmos%T(icell), "ne=", atmos%ne(icell)
         write(*,*) metal%ID, j, metal%n(j,icell), i, metal%n(i,icell)
         write(*,*) "nstar=", metal%n(:,icell)
         write(*,*) " cycling"
         deallocate(gijk, twohnu3_c2, Vij)
         cycle
        end if
        gijk(:) = 0d0
        
    else
      if (metal%n(i,icell) < 0) then
        !-> Futur deprecation of this test, psoitivity is tested when populations are set
        write(*,*) metal%ID, icell, i, " negative pops =", metal%n(i,icell)
        deallocate(gijk, twohnu3_c2, Vij)
        write(*,*) " cycling"
        cycle        
      endif
        gijk(:) = metal%nstar(i,icell)/metal%nstar(j,icell) *  &
     			dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell))
    
    end if


!not storeopac here now
      NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
       				Vij(:) * (metal%n(i,icell)-gijk(:)*metal%n(j,icell))

      NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
      				twohnu3_c2(:) * gijk(:) * Vij(:)*metal%n(j,icell)
  

     deallocate(gijk, twohnu3_c2, Vij)
    end do ! loop over Ncont
  end do !loop over metals

 RETURN
 END SUBROUTINE Metal_bf

 !kept in memory, profile are computed at the centre of the cell for a 0 velocity
 !interpolation turns the CMF profile to the observer's frame profile
 !Zeeman polarisation not included yet
 SUBROUTINE Metal_bb_new (id, icell,x,y,z,x1,y1,z1,u,v,w,l)
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& 
                                				               x1,y1,z1,l
  integer, intent(in) 							            :: icell, id
  integer													:: Nred, Nblue, Nvspace, nv
  integer 													:: kr, m, i, j, nk, kc, la
  type (AtomType), pointer									:: atom
  integer, parameter 										:: NvspaceMax = 1
  real(kind=dp), dimension(NvspaceMax)						:: Omegav
!   real(kind=dp) 											:: delta_vol_phi, xphi, yphi, zphi,&
!   															   v0, v1, dv
  real														:: kappa_neg
  
  
  
  !In principle Nvspace should depend on atom because of atom%vbroad(icell)
  ! v_proj in m/s at point icell
  Omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   !v0 = v_proj(icell,x,y,z,u,v,w) !can be lVoronoi here; for projection
   omegav(1) = v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w)
  end if
  
  !project magnetic field if any
  

!   if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
!      v1 = v_proj(icell,x1,y1,z1,u,v,w)
!      dv = dabs(v1-v0)
!      Nvspace = max(2,nint(20*dv/line%atom%vbroad(icell)))
!      Nvspace = min(Nvspace,NvspaceMax)
!      omegav(Nvspace) = v1
!     do nv=2,Nvspace-1
!       delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
!       xphi=x+delta_vol_phi*u
!       yphi=y+delta_vol_phi*v
!       zphi=z+delta_vol_phi*w
!       omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
!     end do 
!   end if
! write(*,*) "Nvspace=", Nvspace

  do m=1,atmos%Npassiveatoms
  
     atom => atmos%PassiveAtoms(m)%ptr_atom
     
     do kc=1,atom%Ntr_line

     kr = atom%at(kc)%ik
     !line = atom%lines(kr) tà slow to copy an entire struct
     i = atom%lines(kr)%i; j = atom%lines(kr)%j
     Nred = atom%lines(kr)%Nred; Nblue = atom%lines(kr)%Nblue

     !correction factor for negative opacity
     kappa_neg = 1.0 !inner wavelength loop not needed
     if (atom%n(i,icell) - atom%lines(kr)%gij*atom%n(j,icell) < 0.) kappa_neg = 0.
     

     CALL Profile(atom%lines(kr),icell,x,y,z,x1,y1,z1,u,v,w,l,id, Nvspace, Omegav)

     
!         if (atom%lines(kr)%gij*atom%n(j,icell) > atom%n(i,icell)) then 
!            write(*,*) 'line', atom%n(i,icell)/ntotal_atom(icell,atom), &
!            atom%n(j,icell)*atom%lines(kc)%gij/ntotal_atom(icell,atom), atom%lines(kc)%gij
!         endif
     
!      write(*,*) maxval(atom%lines(kr)%phi_loc(:,id)), 1./(atom%vbroad(icell)*SQRTPI), atom%lines(kr)%twohnu3_c2, atom%lines(kr)%Bij
!      write(*,*) atom%lines(kr)%j, atom%lines(kr)%i, atom%n(j,icell), atom%n(i,icell), atom%vbroad(icell)

     NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
     	 kappa_neg * hc_fourPI * atom%lines(kr)%Bij * atom%lines(kr)%phi_loc(:,id) * (atom%n(i,icell)-atom%lines(kr)%gij*atom%n(j,icell))
     
     NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
     	kappa_neg * atom%lines(kr)%twohnu3_c2 * atom%lines(kr)%gij * hc_fourPI * atom%lines(kr)%Bij * atom%lines(kr)%phi_loc(:,id) * atom%n(j,icell)
!      write(*,*) "line ", j,i, " opac = ",hc_fourPI * atom%lines(kr)%Bij /(atom%vbroad(icell)*SQRTPI) * (atom%n(i,icell)-atom%lines(kr)%gij*atom%n(j,icell))
!      write(*,*) "line ", j,i, " eta = ",atom%lines(kr)%twohnu3_c2 * atom%lines(kr)%gij * hc_fourPI * atom%lines(kr)%Bij /(atom%vbroad(icell)*SQRTPI) * atom%n(j,icell)
!      write(*,*) "line ", j,i, " tauloc = ",hc_fourPI * atom%lines(kr)%Bij /(atom%vbroad(icell)*SQRTPI) * (atom%n(i,icell)-atom%lines(kr)%gij*atom%n(j,icell)) * l


!      do kc=1,atmos%PassiveAtoms(m)%ptr_atom%Ntr_line
! 
!      kr = atmos%PassiveAtoms(m)%ptr_atom%at(kc)%ik
!      !line = atom%lines(kr) tà slow to copy an entire struct
!      i = atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%i
!      j = atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%j
!      Nred = atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%Nred
!      Nblue = atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%Nblue
! 
! 
!      !CALL cmf_to_observer(line,icell,x,y,z,x1,y1,z1,u,v,w,l,id,1)
!      CALL Profile(atmos%PassiveAtoms(m)%ptr_atom%lines(kr),icell,x,y,z,x1,y1,z1,u,v,w,l,id,1)
! 
!      NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
!      	 hc_fourPI * atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%Bij * & 
!      	 atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%%phi_loc(:,id) * &
!      	 (atmos%PassiveAtoms(m)%ptr_atom%n(i,icell)-&
!      	 atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%gij*atmos%PassiveAtoms(m)%ptr_atom%n(j,icell))
! 
!      NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
!      	atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%twohnu3_c2 * &
!      	atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%gij * hc_fourPI * &
!      	atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%Bij * &
!      	atmos%PassiveAtoms(m)%ptr_atom%lines(kr)%%phi_loc(:,id) * atmos%PassiveAtoms(m)%ptr_atom%n(j,icell)
 
      
!      do la=1, atom%lines(kr)%Nlambda
!       NLTEspec%AtomOpac%chi_p(Nblue+la-1,id) = NLTEspec%AtomOpac%chi_p(Nblue+la-1,id) + &
!      	 hc_fourPI * atom%lines(kr)%Bij * atom%lines(kr)%%phi_loc(la,id) * (atom%n(i,icell)-atom%lines(kr)%gij*atom%n(j,icell))
! 
!       NLTEspec%AtomOpac%eta_p(Nblue+la-1,id) = NLTEspec%AtomOpac%eta_p(Nblue+la-1,id) + &
!      	atom%lines(kr)%twohnu3_c2 * atom%lines(kr)%gij * hc_fourPI * atom%lines(kr)%Bij * atom%lines(kr)%%phi_loc(la,id) * atom%n(j,icell)
!      enddo
     

    end do !end loop on lines for this atom
    atom => NULL()
  end do !end loop over Natom

 RETURN
 END SUBROUTINE Metal_bb_new

!->futur deprec
 SUBROUTINE Metal_bb (id, icell,x,y,z,x1,y1,z1,u,v,w,l)
  ! Computes the emissivity and extinction of passive lines.
  ! i.e., Atoms with detailed atomic structures read but
  ! not treated in NLTE.
  ! Because damping is wavelength independent and depend only on
  ! the grid (cell) points, here, if line%damping_initialized
  ! do not CALL Damping()
  ! the x,y,z and u,v,w quantities are used to compute the projected velocities at the
  ! cell point we are computing the opacities.
  ! Chip is only computed in Stokes transfer and contains the magneto-optical elements.
  integer 													:: kr, m, i, j, nk, kc
  integer, intent(in) 							            :: icell, id
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  real(kind=dp) 											:: twohnu3_c2, gij
  integer													:: Nred, Nblue, Nvspace, nv
  type (AtomicLine)										    :: line
  type (AtomType)											:: atom
  integer, parameter 										:: NvspaceMax = 1
  real(kind=dp), dimension(NvspaceMax)						:: Omegav

write(*,*) "Omegav not fill"
stop

  do m=1,atmos%Npassiveatoms
   atom = atmos%PassiveAtoms(m)%ptr_atom
     do kc=1,atom%Ntr_line
    !do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     kr = atom%at(kc)%ik
     line = atom%lines(kr)
     i = line%i; j = line%j
     Nred = line%Nred; Nblue = line%Nblue
     !!if (.not.atom%at(kr)%lcontrib_to_opac) CYCLE !avoid lines not defined on the grid


     !Still compute emission even if low, but neglect stm ???
     !if ((atom%n(j,icell) <tiny_dp).or.(atom%n(i,icell) <tiny_dp)) then !no transition
       !!but show the message only if pops is negative
!-> Futur deprecation of this test, psoitivity is tested when populations are set
      if ((atom%n(j,icell) < 0 ).or.(atom%n(i,icell) < 0)) then
        write(*,*) "(Metal_bb) Warning at icell=", icell," T(K)=", atmos%T(icell)
        write(*,*) atom%ID," density <= tiny dp ", i, j, line%lambda0, atom%n(i,icell), atom%n(j,icell)
        write(*,*) "skipping this level"
        write(*,*) "nstar=", atom%nstar(:,icell)
        write(*,*) "n = ", atom%n(:,icell)
        CYCLE
      endif
      !CYCLE
     !end if

     !allocate(Vij(line%Nlambda)); Vij = 0d0
     gij = line%Bji / line%Bij ! = gi/gj
     twohnu3_c2 = line%Aji / line%Bji
     
     !before it was afyer allocate phi
     if (line%voigt) CALL Damping(icell, atom, kr, line%adamp)

     CALL Profile (line, icell,x,y,z,x1,y1,z1,u,v,w,l,id, Nvspace, Omegav)


     !write(*,*) allocated(phiZ), allocated(psiZ), line%polarizable, PRT_SOLUTION

     !Sum up all contributions for this line with the other
     !Vij(:) = hc_fourPI * line%Bij * phi(:)!already normalized / (SQRTPI * VBROAD_atom(icell,atom))
      
     NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + hc_fourPI * line%Bij * line%phi_loc(:,id) * (atom%n(i,icell)-gij*atom%n(j,icell))

     NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + twohnu3_c2 * gij * hc_fourPI * line%Bij * line%phi_loc(:,id) * atom%n(j,icell)

!          write(*,*) "check"
!          write(*,*) gij,  atom%g(i)/atom%g(j)
!          write(*,*) atom%g(i)/atom%g(j)*atom%n(j,icell)/atom%n(i,icell), exp(-hc/KBOLTZMANN/atmos%T(icell)/NM_TO_M/line%lambda0)
!          stop

     if (line%polarizable .and. PRT_SOLUTION == "FULL_STOKES") then
       do nk = 1, 3
         !magneto-optical
         NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) + &
           hc_fourPI * line%Bij * (atom%n(i,icell)-gij*atom%n(j,icell)) * line%psi(nk,:,1)
         !dichroism
         NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) + &
           hc_fourPI * line%Bij * (atom%n(i,icell)-gij*atom%n(j,icell)) * line%psi(nk,:,1)
         !emissivity
         NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) + &
          twohnu3_c2 * gij * hc_fourPI * line%Bij * atom%n(j,icell) * line%phiZ(nk,:,1)
       end do 
     end if
     
    end do !end loop on lines for this atom
  end do !end loop over Natom

 RETURN
 END SUBROUTINE Metal_bb

 
 SUBROUTINE compute_opacities()
 !$ use omp_lib
  integer :: icell, id, icell0
   icell0 = 1

   !$omp parallel &
   !$omp default(none) &
   !$omp private(icell,id,icell0) &
   !$omp shared(atmos)
   !$omp do schedule(dynamic,1)
   do icell=1, atmos%Nspace
    !$ id = omp_get_thread_num() + 1
    if (atmos%icompute_atomRT(icell) > 0) then
     CALL compute_atom_quantities(icell) !need for BackgroundContinua
    									!and lines 
     CALL BackgroundContinua(icell)
          
    !lines profiles projected after
    endif
   end do
   !$omp end do
   !$omp end parallel

 RETURN
 END SUBROUTINE compute_opacities
 
 SUBROUTINE BackgroundContinua (icell)
  integer, intent(in) :: icell
  integer :: la
  real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, Bpnu!,sca

   Bpnu = 0d0
   CALL Bplanck(atmos%T(icell), Bpnu)

   chi = 0d0
   eta = 0d0
   NLTEspec%AtomOpac%Kc(:,icell,:) = 0.0_dp
   NLTEspec%AtomOpac%jc(:,icell) = 0.0_dp
   
   CALL Hminus_bf(icell, chi,eta)
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
   NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + eta

   CALL Hminus_ff(icell, chi)
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
   NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu
  
   CALL Hydrogen_ff(icell, chi)
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
   NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu
 

   if (associated(Helium)) then
    CALL atom_ff_transitions(Helium, icell, chi)
    NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
    NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu
   endif

   if (atmos%Npassiveatoms > 0) CALL Metal_bf_new(icell)
   

   NLTEspec%AtomOpac%Kc(:,icell,2) = Thomson(icell)
 
   CALL HI_Rayleigh(1, icell)
   !!CALL Rayleigh(1, icell, hydrogen)
   if (associated(Helium)) CALL HeI_Rayleigh(1, icell)
   
   !Total opac once source functions are known
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + &
                                     NLTEspec%AtomOpac%Kc(:,icell,2)


 RETURN
 END SUBROUTINE BackgroundContinua
 
!  SUBROUTINE BackgroundLines(id,icell, x, y, z, u, v, w, &
!                                   x1, y1, z1, l)
!   integer, intent(in) :: icell, id
!   real(kind=dp), intent(in) :: x, y, z, u, v, w, &
!                                   x1, y1, z1, l
! 
!    !no need to test, already done
!    !if (atmos%icompute_atomRT(icell)<1) RETURN 
! 
!    !if (atmos%Npassiveatoms == 0) RETURN !tested oustide now
!    CALL Metal_bb_new(id, icell, x, y, z, x1, y1, z1, u, v, w, l)
! 
! 
!  RETURN
!  END SUBROUTINE BackgroundLines
 

!futur deprecation
 SUBROUTINE Background(id,icell,x,y,z,x1,y1,z1,u,v,w,l)
  integer, intent(in) :: icell, id
  real(kind=dp), intent(in) :: x, y, z, u, v, w, &
                                  x1, y1, z1, l!only relevant for b-b when vector fields are present
  real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, Bpnu, chip!, sca

  !test already done
  !if (atmos%icompute_atomRT(icell)<1) RETURN !nH <= tiny_nH or T <= tiny_T == empty cell
  ! all opac are zero, return.

!   if ((atmos%nHtot(icell)==0d0).or.(atmos%T(icell)==0d0)) &
!     RETURN ! stoping for this cell,
            ! it is free of (significant) gas
            ! so no emission/absorption. Coefficients set to 0d0 for all wavelengths
  !Do not forget that it is still possible however, that lcompute_atomRT is TRUE,
  !but that the temperature is too low to have all the levels non zero.
  !The compute_atomRT ensures that at least one level has non-zero populations,
  !to avoid division by zero.

   CALL Bplanck(atmos%T(icell), Bpnu)

   chi = 0d0
   eta = 0d0
   chip = 0d0


   NLTEspec%AtomOpac%sca_c(:,id) = atmos%ne(icell) * sigma_e
   CALL HI_Rayleigh(id, icell)
   !CALL Rayleigh(id, icell, Hydrogen)
   !if (associated(Helium)) CALL Rayleigh(id, icell, Helium)


   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%sca_c(:,id)

   CALL Hydrogen_ff(icell, chi)
   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
   NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + chi * Bpnu

   CALL Hminus_bf(icell, chi, eta)
   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
   NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + eta

   CALL Hminus_ff(icell, chi)
   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
   NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + chi * Bpnu

    !I treat H bound-free in metal_bf. I'll change the name of the subroutine later
!    if (.not.Hydrogen%active) then !passive bound-free !do not enter if active !!!
!     CALL Hydrogen_bf(icell, chi, eta)
!     NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
!     NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + eta
!    end if
   
   if (atmos%Npassiveatoms == 0) RETURN !no passive bound-bound and bound-free

   						
   !--> at this point, eta_p and chi_p are 0 because of initAtomOpac(id), therefore
   !after metal_bf they only points to continuum bound-free.
    CALL Metal_bf(id, icell) !Return if Npassive=1 and PassiveAtoms(1)==" H"
!!     NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
!!     NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + eta

   !keep pure continuum opacities now
   if (MINVAL(NLTEspec%AtomOpac%eta_p(:,id))<0 .or. &
    MINVAL(NLTEspec%AtomOpac%chi_p(:,id)) < 0) then
    write(*,*) "Beware, negative contopac"
   end if
   NLTEspec%AtomOpac%eta_c(:,id) = NLTEspec%AtomOpac%eta_p(:,id)
   NLTEspec%AtomOpac%chi_c(:,id) = NLTEspec%AtomOpac%chi_p(:,id)

   ! we already RETURNs if no passive transitions (H included)
   CALL Metal_bb(id, icell, x, y, z, x1, y1, z1, u, v, w, l)


 RETURN
 END SUBROUTINE Background

END MODULE metal
