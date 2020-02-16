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

	use atmos_type, only				: atmos, Hydrogen, Helium, B_project, vbroad_atom, ntotal_atom
	use constant
	use math, only 						: bezier3_interp, cent_deriv, any_nan_infinity_vector, linear_1D
	use atom_type
	use spectrum_type, only				: NLTEspec, initAtomOpac
	use hydrogen_opacities
	use voigtfunctions, only			: Voigt
	use Profiles!, only					 : Profile
	use broad, only 					: Damping
	use stark, only						: Stark_profile
	use thomson_scattering
	use Rayleigh_scattering
	use Planck
	use getlambda, only					: hv
	use occupation_probability, only	: wocc_n, D_i

	! MCFOST's original
	use mcfost_env, only				: dp
	use molecular_emission, only		: v_proj
	use parametres
	use input
	use constantes, only				: tiny_dp, huge_dp

	IMPLICIT NONE
 
	CONTAINS
 


	SUBROUTINE dealloc_atom_quantities
  
		type(AtomType), pointer :: atom
		integer :: n, k, kc, j, alloc_status

    
		do n=1,atmos%Natom
			atom => atmos%Atoms(n)%ptr_atom
        
			do k = 1, atom%Ntr   
     
				kc = atom%at(k)%ik 
        
				SELECT CASE (atom%at(k)%trtype)
        
					CASE ('ATOMIC_LINE')
						          
          !Damping value
          !no need for it if we do keep the profile on cell
						
						!tmp, need better version of Voigt, but the idea is here
						
						if (atom%lines(kc)%voigt) then !constant for thomson
							if (associated(Profile,Iprofile_thomson)) then
								deallocate(atom%lines(kc)%aeff,atom%lines(kc)%r,atom%lines(kc)%r1)
							endif
          !line profile for v = 0

!           deallocate(atom%lines(kc)%phi_loc)
          
							!Used by all if not kept in memory
							if (allocated(atom%lines(kc)%a)) deallocate(atom%lines(kc)%a) !futur deprec depending of the choice of Voigt

							if (allocated(atom%lines(kc)%u)) deallocate(atom%lines(kc)%u)
						endif
						deallocate(atom%lines(kc)%phi)
						
! 						if (allocated(atom%lines(kc)%Jbar)) deallocate(atom%lines(kc)%Jbar)
						if (allocated(atom%lines(kc)%Rij)) then
							deallocate(atom%lines(kc)%Rij,atom%lines(kc)%Rji)
						endif
						
						if (allocated(atom%lines(kc)%Tex)) deallocate(atom%lines(kc)%Tex)


 
					CASE ('ATOMIC_CONTINUUM')
     
         !do kc=1, atom%Ncont !loop over continua      
          !gij for bound-free transitions
						deallocate(atom%continua(kc)%gij)
						deallocate(atom%continua(kc)%twohnu3_c2)
						deallocate(atom%continua(kc)%alpha)
							
						if (allocated(atom%continua(kc)%Rij)) then
							deallocate(atom%continua(kc)%Rij)
							deallocate(atom%continua(kc)%Rji)
						endif
							
						if (allocated(atom%continua(kc)%Tex)) deallocate(atom%continua(kc)%Tex)


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
		real(kind=dp) :: n_eff, u

    
		do n=1,atmos%Natom
			atom => atmos%Atoms(n)%ptr_atom
        
			do k = 1, atom%Ntr   
     
			kc = atom%at(k)%ik 
        
				SELECT CASE (atom%at(k)%trtype)
        
				CASE ('ATOMIC_LINE')     
     
          !Tex
          !Can be NPROC instead of Nspace ? 
					allocate(atom%lines(kc)%Tex(atmos%Nspace), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error line%Tex")
					atom%lines(kc)%Tex(:) = atmos%T(:)

          
          !some memory can be saved depending of the choice of the profile
          ! I have to do it for large simu
          
          !Damping value
          !no need for it if we do keep the profile on cell
          !except if updated, but if we chose thomson method it is need along to line%aeff
					if (atom%lines(kc)%voigt) then
						if (associated(Profile, Iprofile_cmf_to_obs))then
							allocate(atom%lines(kc)%a(atmos%Nspace), stat=alloc_status)
							if (alloc_status > 0) CALL ERROR("Allocation error line%adamp")
							atom%lines(kc)%a(:) = 0d0
						else if (associated(Profile, Iprofile_thomson)) then
						
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
					endif
          !line profile for v = 0
          !no need to keep it if not interp, no shift. No need also if Thomson (approx), or Iprofile (exact)
!Keep it anyway, as ATM, dk is the default mode even at LTE
! 					if (associated(Voigt, Iprofile_cmf_to_obs).or.atmos%NactiveAtom>0) then
						allocate(atom%lines(kc)%phi(atom%lines(kc)%Nlambda, atmos%Nspace), stat=alloc_status)
						if (alloc_status > 0) CALL ERROR("Allocation error line%phi")
						atom%lines(kc)%phi(:,:) = 0d0
!same
						allocate(atom%lines(kc)%u(atom%lines(kc)%Nlambda), stat=alloc_status)
						if (alloc_status > 0) CALL ERROR("Allocation error line%u")
						!so that argument of profile is line%u / vbroad(icell)
						atom%lines(kc)%u(:) = (NLTEspec%lambda(atom%lines(kc)%Nblue:atom%lines(kc)%Nred)-atom%lines(kc)%lambda0) * &
						CLIGHT / atom%lines(kc)%lambda0 !m/s
! 					endif
          
          !line displacement on the global frequency arrays due to velocity fields.
          !Function of direction, and stored only once per proc if not sequencial 
          !(Other wise, integration of SEE is done cell by cell for all rays)
          !dk is actually compute for each cell during the propagation, but it needs
          !to be stored only for the cell we are solving for the SEE.
! 					if () then
						allocate(atom%lines(kc)%dk(atmos%Nrays, NLTEspec%Nproc), stat=alloc_status)
						if (alloc_status > 0) CALL error ("Allocation error line%dk")
						atom%lines(kc)%dk(:,:) = 0 !no shift
! 					endif
 
!->NO, dk is the default mode. Otherwise an interpolation should be used for accurate results
!But at this point I do not allow to compute exact or approximative (Thomson) Voigt functions
! for NLTE. Interp > Voigt for the same accuracy. Thomson is faster but less accurate. So dk is prefered.
!Can change relatively easily, especially if the profile varies across the cell due to velocity.         
!           allocate(atom%lines(kc)%phi_loc(atom%lines(kc)%Nlambda, NLTEspec%NPROC), stat=alloc_status)
!           if (alloc_status > 0) CALL ERROR("Allocation error line%phi_loc")
          
!-> Redundant with Rij: Rij = Jbar/Bij           
! 					allocate(atom%lines(kc)%Jbar(NLTEspec%NPROC), stat=alloc_status)
! 					if (alloc_status > 0) CALL ERROR("Allocation error line%Jbar")
! 					atom%lines(kc)%Jbar(:) = 0d0
					
					allocate(atom%lines(kc)%Rij(NLTEspec%NPROC), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error line%Rij")
					atom%lines(kc)%Rij(:) = 0d0
					allocate(atom%lines(kc)%Rji(NLTEspec%NPROC), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error line%Rji")
					atom%lines(kc)%Rji(:) = 0d0

                  
				CASE ('ATOMIC_CONTINUUM')

          !Tex
					allocate(atom%continua(kc)%Tex(atmos%Nspace), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error cont%Tex")
					atom%continua(kc)%Tex(:) = atmos%T(:)
      
          !gij for bound-free transitions
					allocate(atom%continua(kc)%gij(atom%continua(kc)%Nlambda, atmos%Nspace), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error cont%gij")
					atom%continua(kc)%gij(:,:) = 0d0
          
					allocate(atom%continua(kc)%Rij(NLTEspec%NPROC), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error cont%Rij")
					atom%continua(kc)%Rij(:) = 0d0
					allocate(atom%continua(kc)%Rji(NLTEspec%NPROC), stat=alloc_status)
					if (alloc_status > 0) CALL ERROR("Allocation error cont%Rji")
					atom%continua(kc)%Rji(:) = 0d0

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
						n_eff = atom%stage(atom%continua(kc)%j)*dsqrt(atom%Rydberg/(atom%E(atom%continua(kc)%j)-atom%E(atom%continua(kc)%i)))

					else !interpolation of the read Cross-section
!       cont is not an alias but a copy of contiua(kc), so cont%alpha not deallocated
        				allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
        				if (alloc_status > 0) CALL ERROR("Allocation error cont%alpha")
         
!          atom%continua(kc)%alpha(:) = linear_1D(size(atom%continua(kc)%lambda_file),atom%continua(kc)%lambda_file,&
!             atom%continua(kc)%alpha_file, atom%continua(kc)%Nlambda,NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))
						CALL bezier3_interp(size(atom%continua(kc)%lambda_file), atom%continua(kc)%lambda_file, atom%continua(kc)%alpha_file, & !read values
     	atom%continua(kc)%Nlambda, NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), atom%continua(kc)%alpha) !interpolation grid

					endif

!tests
!allocate(atom%continua(kc)%negative_opacity(n_cells))
!atom%continua(kc)%negative_opacity(:) = .false.
 
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
		integer :: n, k, kc, alloc_status, i, j, Nblue, Nred, la
		real(kind=dp) :: vbroad, aL, cte, cte2, adamp, wi, wj
    
    
		do n=1,atmos%Natom
			atom => atmos%Atoms(n)%ptr_atom
			do k = 1, atom%Ntr   
     
				kc = atom%at(k)%ik 
        
				SELECT CASE (atom%at(k)%trtype)
        
				CASE ('ATOMIC_LINE')
					vbroad = atom%vbroad(icell)
         
					i=atom%lines(kc)%i; j=atom%lines(kc)%j
					wj = 1.0; wi = 1.0
					if (ldissolve) then
						if (atom%ID=="H") then! .or. atom%ID=="He") then
							wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j))+1.0)
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i))+1.0)
						endif
					endif

 !There are duplicates space array to remove. We do not need aeff if not Thomson profile        
					if (atom%lines(kc)%voigt) then

						CALL Damping(icell, atom, kc, adamp)
          !no need to keep it if we keep the profile for shift or interp
          !But damping can depend on populations of atom and electrons so ...
						if (associated(Profile, Iprofile_cmf_to_obs)) atom%lines(kc)%a(icell) = adamp
          !line%u(:) = atom%lines(kc)%u / atom%vbroad(icell)
          !                             if full profile to be interpolated or shifter.
          !                             if we use another approx for Voigt, it is computed on the fly, not shifted nor interp
!keep it at the moment for all case
						atom%lines(kc)%phi(:,icell) = Voigt(atom%lines(kc)%Nlambda, adamp, atom%lines(kc)%u(:)/vbroad)
!-> do nothing at the moment, but write the futur profile
CALL Stark_profile(icell, atom%lines(kc))
          !write(*,*) "d = ", line%adamp
          !write(*,*) maxval(atom%lines(kc)%phi(:,icell)), minval(atom%lines(kc)%phi(:,icell))
          !Thomson method effective damping for all cells
						aL = adamp * vbroad !(m/s), adamp in doppler units
						if (associated(Profile, Iprofile_thomson)) then
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
						endif
					else !Gaussian
!no need if no interp nor shift or if we use Thomson (but keep adamp for thomson)
						atom%lines(kc)%phi(:,icell) = dexp(-(atom%lines(kc)%u(:)/atom%vbroad(icell))**2)
					endif
!futur test to see if we keep it or not depending on the method of calculation of profiles
					atom%lines(kc)%phi(:,icell) = atom%lines(kc)%phi(:,icell) / (SQRTPI * atom%vbroad(icell))

!Include occupa ?
					if ((atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij) <= 0 ) then
						write(*,*) atom%ID, " nuij (10^15 Hz)", 1d-6 * CLIGHT / atom%lines(kc)%lambda0, " icell=", icell
						write(*,*) " lambdaij (nm)", atom%lines(kc)%lambda0 
						CALL warning ("background line: ni < njgij")
          				write(*,*) atom%n(i,icell), atom%n(j,icell), atom%lines(kc)%gij
          	          				write(*,*) atom%nstar(i,icell), atom%nstar(j,icell), atom%lines(kc)%gij					
					
					endif
    
				CASE ("ATOMIC_CONTINUUM")
        
					i=atom%continua(kc)%i; j=atom%continua(kc)%j
					Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred

					wj = 1.0; wi = 1.0
					if (ldissolve) then
						if (atom%ID=="H") then! .or. atom%ID=="He") then
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(j)))
						endif
					endif
   
					if (atom%nstar(j,icell) < tiny_dp .or. atom%nstar(i,icell) < tiny_dp) then
						atom%continua(kc)%gij(:,icell) = 0d0
					else
						atom%continua(kc)%gij(:, icell) = atom%nstar(i,icell)/atom%nstar(j,icell) *  &
     			dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell))
     			
					endif

					do la=1, atom%continua(kc)%Nlambda

						if (atom%n(i,icell) - atom%n(j, icell)*atom%continua(kc)%gij(la, icell) <= 0 ) then

							write(*,*) atom%ID, " nu+ (10^15 Hz)", 1d-6 * CLIGHT / atom%continua(kc)%lambda0, " icell=", icell
							write(*,*) " lambda+ (nm)", atom%continua(kc)%lambda0
							CALL warning ("background cont: ni < njgij")

							exit
						endif
					enddo
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

	 !compute the dissolve fraction for pseudo -continua only in hydrogenic case.
	 !Actually, for explicit continua, because lambda <= lambda0, the dissolve fraction
	 !is  one. For hydrogenic with lambda > lambda0 it is < 1
  integer, intent(in)							            :: icell
  integer                                                   :: m, kr, kc, i, j, Nblue, Nred, la
  type (AtomType), pointer                                  :: atom
  real(kind=dp)                                             :: n_eff, wj, wi, chi
  real(kind=dp)                                             :: Diss, chi_ion
  
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
!     do kc=atom%Ntr_line+5,atom%Ntr_line+5
     kr = atom%at(kc)%ik 
     
    !do kr=1,atom%Ncont
     i = atom%continua(kr)%i
     j = atom%continua(kr)%j 
! write(*,*)"lam0=", atom%continua(kr)%lambda0, "chi0=",(atom%E(j) - atom%E(i))/EV
! write(*,*) atmos%T(icell), atom%n(i,icell), atom%n(j,icell)*atom%nstar(i,icell)/atom%nstar(j,icell)

     Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred
     n_eff = real(atom%stage(j)) * dsqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))

			wj = 1.0; wi = 1.0
			if (ldissolve) then
				if (atom%ID=="H") then! .or. atom%ID=="He") then
												!n_eff
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(j)))
				endif
			endif
     
     !get the ionisation potential for the ion ot be use in the dissolve fraction
     chi_ion = atmos%Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(i)+1)

    !explicit do loop to avoid a internal loop in D_i. Should be similar to the implicit do loop 
    do la=1, atom%continua(kr)%nlambda
! write(*,*) "lam=", NLTEspec%lambda(Nblue+la-1)
     if ( (atom%n(i,icell) - atom%n(j,icell) * wi * atom%continua(kr)%gij(la,icell)) <= 0.0) cycle
     
     Diss = D_i(icell, n_eff, real(atom%stage(i)), &
       		   1.0, NLTEspec%lambda(Nblue+la-1), atom%continua(kr)%lambda0, chi_ion)

! write(*,*) "np=", atom%n(j,icell), "sigma=", atom%continua(kr)%alpha(la), "gijk=", atom%continua(kr)%gij(la,icell), &
! "ni=", atom%n(i, icell), "2hnu3_c2=", atom%continua(kr)%twohnu3_c2(la)
! write(*,*) "icell=", icell, "chi=",      Diss * atom%continua(kr)%alpha(la) * (atom%n(i,icell) - wi * atom%continua(kr)%gij(la,icell)*atom%n(j,icell)) , &
!  "eta=", Diss * atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * atom%continua(kr)%gij(la,icell) * atom%n(j,icell) * wi				
!      
     NLTEspec%AtomOpac%Kc(Nblue+la-1,icell) = NLTEspec%AtomOpac%Kc(Nblue+la-1,icell) + &
        Diss * atom%continua(kr)%alpha(la) * (atom%n(i,icell) - wi * atom%continua(kr)%gij(la,icell)*atom%n(j,icell)) 
       				
     NLTEspec%AtomOpac%jc(Nblue+la-1,icell) = NLTEspec%AtomOpac%jc(Nblue+la-1,icell) + &
       Diss * atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * atom%continua(kr)%gij(la,icell) * atom%n(j,icell) * wi
!      NLTEspec%AtomOpac%jc(Nblue+la-1,icell) = NLTEspec%AtomOpac%jc(Nblue+la-1,icell) + &
!        atom%n(i,icell) * Diss * atom%continua(kr)%alpha(la) * bpnu(atmos%T(icell), NLTEspec%lambda(Nblue+la-1))
 
!     NLTEspec%AtomOpac%Kc(Nblue+la-1,icell) = NLTEspec%AtomOpac%Kc(Nblue+la-1,icell) + &
!        1d-22 * atom%n(i,icell)	
! 		
!      NLTEspec%AtomOpac%jc(Nblue+la-1,icell) = NLTEspec%AtomOpac%jc(Nblue+la-1,icell) + &
!        bpnu(atmos%T(icell), NLTEspec%lambda(Nblue+la-1))*1d-22 * atom%n(i,icell)
! if (icell==atmos%nspace-2) then
!     write(*,*) atmos%T(icell), atmos%ne(icell), atom%nstar(i,icell)
!     write(*,*) "n=", kr, "lambda=", NLTEspec%lambda(Nblue+la-1), "chi=", 1d-22 * atom%n(i,icell)
!     write(*,*) "chib=", atom%continua(kr)%alpha(la) * (atom%n(i,icell) - wi * atom%continua(kr)%gij(la,icell)*atom%n(j,icell))
!     write(*,*) "Xs=", atom%continua(kr)%alpha(la), "stm=",atom%nstar(j,icell)/atom%nstar(i,icell)*atom%continua(kr)%gij(la,icell), &
!     "stm2=", dexp(-hc_k/atmos%T(icell)/NLTEspec%lambda(Nblue+la-1))
! endif    				
    enddo

! if (icell==atmos%nspace-2) stop
    end do ! loop over Ncont
    atom => NULL()
  end do !loop over metals

 RETURN
 END SUBROUTINE Metal_bf_new

 SUBROUTINE Metal_bb_new (id, icell,x,y,z,x1,y1,z1,u,v,w,l)
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& 
                                				               x1,y1,z1,l
  integer, intent(in) 							            :: icell, id
  integer													:: Nred, Nblue, Nvspace, nv
  integer 													:: kr, m, i, j, nk, kc, la, dk
  type (AtomType), pointer									:: atom
  integer, parameter 										:: NvspaceMax = 1
  real(kind=dp), dimension(NvspaceMax)						:: Omegav
!   real(kind=dp) 											:: delta_vol_phi, xphi, yphi, zphi,&
!   															   v0, v1, dv
  real(kind=dp)												:: n_eff, wi, wj, phi

  !In principle Nvspace should depend on atom because of atom%vbroad(icell)
  ! v_proj in m/s at point icell
  Omegav = 0d0
  Nvspace = 1
  dk = 0

  !v0 = v_proj(icell,x,y,z,u,v,w) !can be lVoronoi here; for projection
  omegav(1) = v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w)

  
  dk = nint(1e-3 * omegav(1) / hv)


  do m=1,atmos%Npassiveatoms
  
     atom => atmos%PassiveAtoms(m)%ptr_atom
     
     do kc=1,atom%Ntr_line

     kr = atom%at(kc)%ik
     
     !check that dk is the good one for SEE calculation, needed to integrate
     !I(Nblue-dk:Nred+dk)*phi(:)
     !atom%lines(kr)%dk(1,id) = dk !there should be rays or not useful in dk ?? since integration is ray by ray
     
     !line = atom%lines(kr) tà slow to copy an entire struct
     i = atom%lines(kr)%i; j = atom%lines(kr)%j
     Nred = atom%lines(kr)%Nred; Nblue = atom%lines(kr)%Nblue

			wj = 1.0; wi = 1.0
			if (ldissolve) then
				if (atom%ID=="H") then! .or. atom%ID=="He") then
												!n_eff
					wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j))+1.0)
					wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i))+1.0)	
				endif
			endif

     if ( (wj*atom%n(i,icell) - wi*atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0 ) cycle


!      CALL Profile(atom%lines(kr),icell,x,y,z,x1,y1,z1,u,v,w,l,id, Nvspace, Omegav)

!!      NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
!!      	 hc_fourPI * atom%lines(kr)%phi_loc(:,id) * atom%lines(kr)%Bij * (atom%n(i,icell)-wi*atom%lines(kr)%gij*atom%n(j,icell) ) 
!!      
!!      NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
!!      	 hc_fourPI * atom%lines(kr)%Aji * atom%lines(kr)%phi_loc(:,id) * atom%n(j,icell) * wi

!      do la=1, atom%lines(kr)%Nlambda
!        NLTEspec%AtomOpac%chi_p(Nblue+la-1,id) = NLTEspec%AtomOpac%chi_p(Nblue+la-1,id) + &
!      	 hc_fourPI * atom%lines(kr)%phi_loc(la,id) * atom%lines(kr)%Bij * (wj * atom%n(i,icell)- wi * atom%lines(kr)%gij*atom%n(j,icell) ) 
!      
!        NLTEspec%AtomOpac%eta_p(Nblue+la-1,id) = NLTEspec%AtomOpac%eta_p(Nblue+la-1,id) + &
!      	hc_fourPI * atom%lines(kr)%Aji * atom%lines(kr)%phi_loc(la,id) * atom%n(j,icell) * wi
!      enddo
     
     do la=1, atom%lines(kr)%Nlambda
     
       NLTEspec%AtomOpac%chi_p(Nblue+la-1-dk,id) = NLTEspec%AtomOpac%chi_p(Nblue+la-1-dk,id) + &
     	hc_fourPI * atom%lines(kr)%phi(la,icell) * atom%lines(kr)%Bij * (wj * atom%n(i,icell)- wi * atom%lines(kr)%gij*atom%n(j,icell) ) 
     
       NLTEspec%AtomOpac%eta_p(Nblue+la-1-dk,id) = NLTEspec%AtomOpac%eta_p(Nblue+la-1-dk,id) + &
     	hc_fourPI * atom%lines(kr)%Aji * atom%lines(kr)%phi(la,icell) * atom%n(j,icell) * wi

     enddo


    end do !end loop on lines for this atom
    atom => NULL()
  end do !end loop over Natom

 RETURN
 END SUBROUTINE Metal_bb_new


 
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
		integer :: la, nat
		real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, Bpnu!,sca

		Bpnu = 0d0
		CALL Bplanck(atmos%T(icell), Bpnu)

		chi = 0d0
		eta = 0d0
		NLTEspec%AtomOpac%Kc(:,icell) = 0.0_dp
		NLTEspec%AtomOpac%jc(:,icell) = 0.0_dp
		NLTEspec%AtomOpac%sca_c(:,icell) = 0.0_dp
   
		CALL Hminus_bf_wishart(icell, chi,eta)
		NLTEspec%AtomOpac%Kc(:,icell) = NLTEspec%AtomOpac%Kc(:,icell) + chi
		NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + eta

		CALL Hminus_ff_bell_berr(icell, chi)
		NLTEspec%AtomOpac%Kc(:,icell) = NLTEspec%AtomOpac%Kc(:,icell) + chi
		NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu

		CALL Hydrogen_ff(icell, chi)
		NLTEspec%AtomOpac%Kc(:,icell) = NLTEspec%AtomOpac%Kc(:,icell) + chi
		NLTEspec%AtomOpac%jc(:,icell) = NLTEspec%AtomOpac%jc(:,icell) + chi * Bpnu

		if (atmos%Npassiveatoms > 0) CALL Metal_bf_new(icell)

		NLTEspec%AtomOpac%sca_c(:,icell) = Thomson(icell)

!    CALL HI_Rayleigh(1, icell)
!    !CALL Rayleigh(1, icell, hydrogen)
!    if (associated(Helium)) then
!     !write(*,*) " Rayleigh on neutral helium"
!     CALL HeI_Rayleigh(1, icell)
!    endif
   
		!Total opac once source functions are known
		NLTEspec%AtomOpac%Kc(:,icell) = NLTEspec%AtomOpac%Kc(:,icell) + &
                                     NLTEspec%AtomOpac%sca_c(:,icell)


	RETURN
	END SUBROUTINE BackgroundContinua
 

END MODULE metal
