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

 use atmos_type, only                : atmos, Hydrogen, Helium, B_project, vbroad_atom
 use constant
 use math, only                      : bezier3_interp, cent_deriv, any_nan_infinity_vector
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
           deallocate(atom%lines(kc)%aeff)
          endif
          !line profile for v = 0
          deallocate(atom%lines(kc)%phi)
          !deallocate(atom%lines(kc)%phi_loc) -> used in NLTEloop only now
          
          deallocate(atom%lines(kc)%u)

        
        CASE ('ATOMIC_CONTINUUM')
     
         !do kc=1, atom%Ncont !loop over continua
          cont = atom%continua(kc)
      
          !gij for bound-free transitions
          deallocate(atom%continua(kc)%gij)
          deallocate(atom%continua(kc)%twohnu3_c2)
          deallocate(atom%continua(kc)%alpha)


 
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
    integer :: n, k, kc, j, alloc_status
    type(AtomicContinuum) :: cont
    type(AtomicLine) :: line
!     real(kind=dp), allocatable :: dumm_alpha(:)

    
    do n=1,atmos%Natom
     atom => atmos%Atoms(n)%ptr_atom
        
     do k = 1, atom%Ntr   
     
       kc = atom%at(k)%ik 
        
       SELECT CASE (atom%at(k)%trtype)
        
        CASE ('ATOMIC_LINE')
          line = atom%lines(kc)
          
          !some memory can be saved depending of the choice of the profile
          ! I have to do it for large simu
          
          !Damping value
          !no need for it if we do keep the profile on cell
          !except if updated, but if we chose thomson method it is need along to line%aeff
          if (line%voigt) then
           allocate(atom%lines(kc)%a(atmos%Nspace), stat=alloc_status)
           if (alloc_status > 0) CALL ERROR("Allocation error line%adamp")
           atom%lines(kc)%a(:) = 0d0
           allocate(atom%lines(kc)%aeff(atmos%Nspace), stat=alloc_status)
           if (alloc_status > 0) CALL ERROR("Allocation error line%eff")
           atom%lines(kc)%aeff(:) = 0d0
          endif
          !line profile for v = 0
          !no need to keep it if not interp, no shift. No need also if Thomson
          allocate(atom%lines(kc)%phi(line%Nlambda, atmos%Nspace), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%phi")
          atom%lines(kc)%phi(:,:) = 0d0
          !-> use in NLTEloop only now
          !allocate(atom%lines(kc)%phi_ray(line%Nlambda, atmos%Nrays, NLTEspec%NPROC), stat=alloc_status)
          !if (alloc_status > 0) CALL ERROR("Allocation error line%phi_ray")
          allocate(atom%lines(kc)%phi_loc(line%Nlambda, NLTEspec%NPROC), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%phi_loc")
          
          
          allocate(atom%lines(kc)%u(line%Nlambda), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error line%u")
          !so that argument of profile is line%u / vbroad(icell)
          atom%lines(kc)%u(:) = (NLTEspec%lambda(line%Nblue:line%Nred)-line%lambda0) * &
           CLIGHT / line%lambda0 !m/s
        
        CASE ('ATOMIC_CONTINUUM')
     
         !do kc=1, atom%Ncont !loop over continua
          cont = atom%continua(kc)
      
          !gij for bound-free transitions
          allocate(atom%continua(kc)%gij(cont%Nlambda, atmos%Nspace), stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%gij")
          atom%continua(kc)%gij(:,:) = 0d0
          !twohnu3_c2 for continuum waves
          allocate(atom%continua(kc)%twohnu3_c2(cont%Nlambda),stat=alloc_status)
          if (alloc_status > 0) CALL ERROR("Allocation error cont%twohnu3_c2")
          atom%continua(kc)%twohnu3_c2(:) = twohc / NLTEspec%lambda(cont%Nblue:cont%Nred)**3
          
          !special for bound-free cross-sections who do not depend on cell
    		if (cont%Hydrogenic) then !Kramer's formula with quantum mechanical correction
        		allocate(atom%continua(kc)%alpha(cont%Nlambda), stat=alloc_status)
        		if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")
        
        		atom%continua(kc)%alpha(:) = 0d0
        		!computes only where it is not zero
        		atom%continua(kc)%alpha(:) = H_bf_Xsection(cont, NLTEspec%lambda(cont%Nblue:cont%Nred))
        	
      		else !interpolation of the read Cross-section
!       cont is not an alias but a copy of contiua(kc), so cont%alpha not deallocated
        		allocate(atom%continua(kc)%alpha(cont%Nlambda), stat=alloc_status)
        		if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")
        !alpha = interp_dp(cont%alpha, cont%lambda, NLTEspec%lambda(cont%Nblue:cont%Nred))
        
        		CALL bezier3_interp(size(cont%lambda_file), cont%lambda_file, cont%alpha_file, & !read values
     	cont%Nlambda, NLTEspec%lambda(cont%Nblue:cont%Nred), atom%continua(kc)%alpha) !interpolation grid

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
    integer :: n, k, kc, alloc_status, i, j
    real(kind=dp) :: vbroad, aL
    type(AtomicContinuum) :: cont !takes time to allocate these cont and line
    type(AtomicLine) :: line      ! but we do it beforehand, before the actual transfer
    							  !so that's okay
    
    
    do n=1,atmos%Natom
     atom => atmos%Atoms(n)%ptr_atom
     do k = 1, atom%Ntr   
     
       kc = atom%at(k)%ik 
        
       SELECT CASE (atom%at(k)%trtype)
        
        CASE ('ATOMIC_LINE')
         line = atom%lines(kc)
         vbroad = atom%vbroad(icell)
         
         if (line%voigt) then
          CALL Damping(icell, atom, kc, line%adamp)
          !no need to keep it if we keep the profile for shift or interp
          atom%lines(kc)%a(icell) = line%adamp
          !line%u(:) = atom%lines(kc)%u / atom%vbroad(icell)
          !                             if full profile to be interpolated or shifter.
          !                             if we use another approx for Voigt, it is computed on the fly, not shifted nor interp
          atom%lines(kc)%phi(:,icell) = Voigt(line%Nlambda, line%adamp, line%u(:)/vbroad)!/atom%vbroad(icell))
          !write(*,*) "d = ", line%adamp
          !write(*,*) maxval(atom%lines(kc)%phi(:,icell)), minval(atom%lines(kc)%phi(:,icell))
          !Thomson method effective damping for all cells
          aL = line%adamp * vbroad !(m/s), adamp in doppler units
          atom%lines(kc)%aeff(icell) = (vbroad**5. + 2.69269*vbroad**4. * aL + 2.42843*vbroad**3. * aL**2. + &
          		4.47163*vbroad**2.*aL**3. + 0.07842*vbroad*aL**4. + aL**5.)**(0.2) !1/5
         else !Gaussian
          !no need if no interp nor shift or if we use Thomson (but keep adamp for thomson)
          atom%lines(kc)%phi(:,icell) = dexp(-(line%u(:)/atom%vbroad(icell))**2)
         endif
         !futur test to see if we keep it or not depending on the method of calculation of profiles
         atom%lines(kc)%phi(:,icell) = atom%lines(kc)%phi(:,icell) / (SQRTPI * atom%vbroad(icell))
        
        CASE ("ATOMIC_CONTINUUM")
        
         cont = atom%continua(kc)
         i=cont%i; j=cont%j
   
         if (atom%nstar(j,icell) < tiny_dp .or. atom%nstar(i,icell) < tiny_dp) then
          atom%continua(kc)%gij(:,icell) = 0d0
         else
          atom%continua(kc)%gij(:, icell) = atom%nstar(i,icell)/atom%nstar(j,icell) *  &
     			dexp(-hc_k/NLTEspec%lambda(cont%Nblue:cont%Nred)/atmos%T(icell))
     			
         endif
        
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
  integer                                                   :: m, kr, kc, i, j, Nblue, Nred
  type (AtomType), pointer                                  :: metal

  ! Go throught all bound-free transitions of each PASSIVE
  ! metal and add the opacity and emissivity if lambda
  ! is lower (greater) than the wavelength threshold lambdaEdge
  ! (the frequency threshold) and if greater (lower) than
  ! wavelength min (frequency max). See Hydrogen b-f for more
  ! informations, and Hubeny & Mihalas chap. 7

  do m=1,atmos%Npassiveatoms
  ! run over all passive atoms
   metal => atmos%PassiveAtoms(m)%ptr_atom!atmos%Atoms(m)
   !if (metal%ID == "H ") CYCLE !H cont is treated in Hydrogen_bf()
    do kc=metal%Ntr_line+1,metal%Ntr
     kr = metal%at(kc)%ik 
    !do kr=1,metal%Ncont
     i = metal%continua(kr)%i
     j = metal%continua(kr)%j 
     Nblue = metal%continua(kr)%Nblue; Nred = metal%continua(kr)%Nred



     NLTEspec%AtomOpac%Kc(Nblue:Nred,icell,1) = NLTEspec%AtomOpac%Kc(Nblue:Nred,icell,1) + &
       				metal%continua(kr)%alpha(:) * (metal%n(i,icell)-metal%continua(kr)%gij(:,icell)*metal%n(j,icell))
       				
     NLTEspec%AtomOpac%jc(Nblue:Nred,icell) = NLTEspec%AtomOpac%jc(Nblue:Nred,icell) + &
       				metal%continua(kr)%twohnu3_c2(:) * metal%continua(kr)%gij(:,icell) * metal%continua(kr)%alpha(:) * metal%n(j,icell)


    end do ! loop over Ncont
    metal => NULL()
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
  integer													:: Nred, Nblue
  integer 													:: kr, m, i, j, nk, kc, la
  type (AtomType), pointer									:: atom


  do m=1,atmos%Npassiveatoms
     atom => atmos%PassiveAtoms(m)%ptr_atom
     
     do kc=1,atom%Ntr_line

     kr = atom%at(kc)%ik
     !line = atom%lines(kr) tà slow to copy an entire struct
     i = atom%lines(kr)%i; j = atom%lines(kr)%j
     Nred = atom%lines(kr)%Nred; Nblue = atom%lines(kr)%Nblue


     CALL Profile(atom%lines(kr),icell,x,y,z,x1,y1,z1,u,v,w,l,id)
     
!      write(*,*) maxval(atom%lines(kr)%phi_loc(:,id)), 1./(atom%vbroad(icell)*SQRTPI), atom%lines(kr)%twohnu3_c2, atom%lines(kr)%Bij
!      write(*,*) atom%lines(kr)%j, atom%lines(kr)%i, atom%n(j,icell), atom%n(i,icell), atom%vbroad(icell)

     NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
     	 hc_fourPI * atom%lines(kr)%Bij * atom%lines(kr)%phi_loc(:,id) * (atom%n(i,icell)-atom%lines(kr)%gij*atom%n(j,icell))

     NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
     	atom%lines(kr)%twohnu3_c2 * atom%lines(kr)%gij * hc_fourPI * atom%lines(kr)%Bij * atom%lines(kr)%phi_loc(:,id) * atom%n(j,icell)
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
  integer													:: Nred, Nblue
  type (AtomicLine)										    :: line
  type (AtomType)											:: atom



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

     CALL Profile (line, icell,x,y,z,x1,y1,z1,u,v,w,l,id)


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
!    CALL Metal_bb_new(1, icell0, 0d0,0d0,0d0 ,0d0, 0d0, 0d0, 0d0, 0D0, 0d0, 566893968.924461_dp)
!    write(*,*) maxval(NLTEspec%AtomOpac%Kc(:,icell0,1)), maxval(NLTEspec%AtomOpac%jc(:,icell0)), maxval(NLTEspec%AtomOpac%Kc(:,icell0,2))
!    write(*,*) maxval(NLTEspec%AtomOpac%chi_p(:,1)), maxval(NLTEspec%AtomOpac%eta_p(:,1))
!   stop
  
 RETURN
 END SUBROUTINE compute_opacities
 
 SUBROUTINE BackgroundContinua (icell)
  integer, intent(in) :: icell
  real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, Bpnu!,sca

   Bpnu = 0d0
   CALL Bplanck(atmos%T(icell), Bpnu)

   chi = 0d0
   eta = 0d0

   NLTEspec%AtomOpac%Kc(:,icell,1) = atmos%ne(icell) * sigma_e !Thomson(icell)

!    !Not working properly if frequencies are removed
   !CALL Rayleigh(1, icell, Hydrogen)
   !if (associated(Helium)) CALL Rayleigh(1, icell, Helium)
   CALL HI_Rayleigh(1, icell)
   if (associated(Helium)) CALL HeI_Rayleigh(1, icell)

! 
   NLTEspec%AtomOpac%Kc(:,icell,2) = NLTEspec%AtomOpac%Kc(:,icell,1)

   CALL Hydrogen_ff(icell, chi)
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
   NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu
   
   CALL Hminus_bf(icell, chi,eta)
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
   NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + eta


   CALL Hminus_ff(icell, chi)
   NLTEspec%AtomOpac%Kc(:,icell,1) = NLTEspec%AtomOpac%Kc(:,icell,1) + chi
   NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu


   if (atmos%Npassiveatoms == 0) RETURN !no passive bound-bound and bound-free
   CALL Metal_bf_new(icell) !here we do not need id, because lstore_atom_opac=.true.


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
