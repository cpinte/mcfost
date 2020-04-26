module Opacity

	use atmos_type, only								: hydrogen, ntotal_atom, T, ne, NactiveAtoms, Natom, Npassiveatoms, ActiveAtoms, PassiveAtoms, Atoms, Elements, &
															icompute_atomRT
	use atom_type
	use spectrum_type, only								: Itot, Icont, lambda, nlambda, Nlambda_cont, lambda_cont, dk, dk_min, dk_max, &
															chi0_bb, eta0_bb, chi_c, sca_c, eta_c, chi, eta, chi_c_nlte, eta_c_nlte, Psi, Stot, chitot, Jnu, Jnu_cont
	use constant
	use constantes, only								: tiny_dp, huge_dp, AU_to_m
	use messages
	use broad, only										: Damping
	use parametres
	use voigtfunctions, only							: Voigt
	use stark, only										: Stark_profile
	use profiles!, only									: Profile
	use getlambda, only									: hv
	use planck, only 									: bpnu
	use occupation_probability, only					: wocc_n, D_i

	use background_opacity, only						: Thomson, Hydrogen_ff, Hminus_bf, Hminus_bf_geltman, &
															Hminus_bf_geltman, Hminus_bf_Wishart, Hminus_ff_john, Hminus_ff, Hminus_ff_bell_berr, lte_bound_free, H_bf_Xsection
	use mcfost_env, only								: dp
	use input, only										: ds
	use molecular_emission, only						: v_proj


 


	implicit none
 
	contains
	
	subroutine dealloc_atom_quantities
  
		type(AtomType), pointer :: atom
		integer :: n, k, kc, j, alloc_status

    
		do n=1,Natom
			atom => Atoms(n)%ptr_atom
        
			do k = 1, atom%Ntr   
     
				kc = atom%at(k)%ik 
        
				select case (atom%at(k)%trtype)
        
					case ('ATOMIC_LINE')
						
						if (atom%lines(kc)%voigt) then !constant for thomson
							if (associated(Profile,Iprofile_thomson)) then
								deallocate(atom%lines(kc)%aeff,atom%lines(kc)%r,atom%lines(kc)%r1)
							endif
          
							!Used by all if not kept in memory
							if (allocated(atom%lines(kc)%a)) deallocate(atom%lines(kc)%a) !futur deprec depending of the choice of Voigt

							if (allocated(atom%lines(kc)%u)) deallocate(atom%lines(kc)%u)
						endif
						if (allocated(atom%lines(kc)%phi)) deallocate(atom%lines(kc)%phi)
						if (allocated(atom%lines(kc)%phi_loc)) deallocate(atom%lines(kc)%phi_loc)
						
						if (allocated(atom%lines(kc)%Rij)) then
							deallocate(atom%lines(kc)%Rij,atom%lines(kc)%Rji)
						endif
						
						if (allocated(atom%lines(kc)%Tex)) deallocate(atom%lines(kc)%Tex)

 
					case ('ATOMIC_CONTINUUM')
     

						if (allocated(atom%continua(kc)%gij)) deallocate(atom%continua(kc)%gij)
						if (allocated(atom%continua(kc)%alpha_nu)) deallocate(atom%continua(kc)%alpha_nu)
						if (allocated(atom%continua(kc)%twohnu3_c2)) deallocate(atom%continua(kc)%twohnu3_c2)
						if (allocated(atom%continua(kc)%alpha)) deallocate(atom%continua(kc)%alpha)
							
						if (allocated(atom%continua(kc)%Rij)) then
							deallocate(atom%continua(kc)%Rij)
							deallocate(atom%continua(kc)%Rji)
						endif
							
						if (allocated(atom%continua(kc)%Tex)) deallocate(atom%continua(kc)%Tex)


					case default
						call Error("Transition type unknown", atom%at(k)%trtype)
					end select
				enddo
				atom => NULL()
			enddo  
  
	return
	end subroutine dealloc_atom_quantities
  
   !I will not allocate twohnu3_c2, gij for lte bound free but allocate alpha for all at the moment
	subroutine alloc_atom_quantities
  
		type(AtomType), pointer :: atom
		integer :: n, k, kc, j, alloc_status, la, n_cells_d, icell
		real(kind=dp) :: n_eff, u

    
		do n=1,Natom
			atom => Atoms(n)%ptr_atom
        
			do k = 1, atom%Ntr   
     
			kc = atom%at(k)%ik 
        
				select case (atom%at(k)%trtype)
        
				case ('ATOMIC_LINE')     

					if (atom%active) then
						allocate(atom%lines(kc)%Tex(n_cells), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%Tex")
						atom%lines(kc)%Tex(:) = T(:)
					endif
          
          !some memory can be saved depending of the choice of the profile
          ! I have to do it for large simu
          
          !Damping value
          !no need for it if we do keep the profile on cell
          !except if updated, but if we chose thomson method it is need along to line%aeff
					if (atom%lines(kc)%voigt) then
							allocate(atom%lines(kc)%a(n_cells), stat=alloc_status)
							if (alloc_status > 0) call ERROR("Allocation error line%adamp")
							atom%lines(kc)%a(:) = 0d0
						if (associated(Profile, Iprofile_thomson)) then
! 						else if (associated(Profile, Iprofile_thomson)) then
						
							allocate(atom%lines(kc)%aeff(n_cells), stat=alloc_status)
							if (alloc_status > 0) call ERROR("Allocation error line%eff")
							atom%lines(kc)%aeff(:) = 0d0
							allocate(atom%lines(kc)%r(n_cells), stat=alloc_status)
							if (alloc_status > 0) call ERROR("Allocation error line%r")
							atom%lines(kc)%r(:) = 0d0
							allocate(atom%lines(kc)%r1(n_cells), stat=alloc_status)
							if (alloc_status > 0) call ERROR("Allocation error line%r1")
							atom%lines(kc)%r1(:) = 0d0
						endif
					endif
					
          !line profile for v = 0
          !no need to keep it if not interp, no shift. No need also if Thomson (approx), or Iprofile (exact)
!Keep it anyway, as ATM, dk is the default mode even at LTE
! 					if (associated(Voigt, Iprofile_cmf_to_obs).or.atmos%NactiveAtom>0) then
						allocate(atom%lines(kc)%phi(atom%lines(kc)%Nlambda, n_cells), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%phi")
						atom%lines(kc)%phi(:,:) = 0d0
!same
						allocate(atom%lines(kc)%u(atom%lines(kc)%Nlambda), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%u")
						!so that argument of profile is line%u / vbroad(icell)
						atom%lines(kc)%u(:) = (lambda(atom%lines(kc)%Nblue:atom%lines(kc)%Nred)-atom%lines(kc)%lambda0) * &
						CLIGHT / atom%lines(kc)%lambda0 !m/s
! 					endif
          
          !line displacement on the global frequency arrays due to velocity fields.
          !Function of direction, and stored only once per proc if not sequencial 
          !(Other wise, integration of SEE is done cell by cell for all rays)
          !dk is actually compute for each cell during the propagation, but it needs
          !to be stored only for the cell we are solving for the SEE.
! 					if () then
! 						allocate(atom%lines(kc)%dk(atmos%Nrays, NLTEspec%Nproc), stat=alloc_status)
! 						if (alloc_status > 0) call error ("Allocation error line%dk")
! 						atom%lines(kc)%dk(:,:) = 0 !no shift
! 					endif
 
!->NO, dk is the default mode. Otherwise an interpolation should be used for accurate results
!But at this point I do not allow to compute exact or approximative (Thomson) Voigt functions
! for NLTE. Interp > Voigt for the same accuracy. Thomson is faster but less accurate. So dk is prefered.
!Can change relatively easily, especially if the profile varies across the cell due to velocity.         
!           allocate(atom%lines(kc)%phi_loc(atom%lines(kc)%Nlambda, NLTEspec%NPROC), stat=alloc_status)
!           if (alloc_status > 0) call ERROR("Allocation error line%phi_loc")
          
!-> Redundant with Rij: Rij = Jbar/Bij           
! 					allocate(atom%lines(kc)%Jbar(NLTEspec%NPROC), stat=alloc_status)
! 					if (alloc_status > 0) call ERROR("Allocation error line%Jbar")
! 					atom%lines(kc)%Jbar(:) = 0d0
					
					if (atom%active) then
						allocate(atom%lines(kc)%Rij(nb_proc), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%Rij")
						atom%lines(kc)%Rij(:) = 0d0
						allocate(atom%lines(kc)%Rji(nb_proc), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%Rji")
						atom%lines(kc)%Rji(:) = 0d0
					endif

                  
				case ('ATOMIC_CONTINUUM')

          !Tex
          			if (atom%active) then
						allocate(atom%continua(kc)%Tex(n_cells), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%Tex")
						atom%continua(kc)%Tex(:) = T(:)
					endif
      
          !gij for bound-free transitions
          			if (atom%active) then
! 						allocate(atom%continua(kc)%gij(atom%continua(kc)%Nlambda, n_cells), stat=alloc_status)
! 						if (alloc_status > 0) call ERROR("Allocation error cont%gij")
! 						atom%continua(kc)%gij(:,:) = 0d0

! 						allocate(atom%continua(kc)%sigma(atom%continua(kc)%Nlambda, n_cells), stat=alloc_status)
! 						if (alloc_status > 0) call ERROR("Allocation error cont%gij")
! 						atom%continua(kc)%sigma(:,:) = 0d0
          
						allocate(atom%continua(kc)%Rij(nb_proc), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%Rij")
						atom%continua(kc)%Rij(:) = 0d0
						allocate(atom%continua(kc)%Rji(nb_proc), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%Rji")
						atom%continua(kc)%Rji(:) = 0d0

          !twohnu3_c2 for continuum waves (Not needed ?)
						allocate(atom%continua(kc)%twohnu3_c2(atom%continua(kc)%Nlambda),stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%twohnu3_c2")
						atom%continua(kc)%twohnu3_c2(:) = twohc / lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred)**3.
          			endif
          			
          			
          !special for bound-free cross-sections who do not depend on cell
					if (atom%continua(kc)%Hydrogenic) then !Kramer's formula with quantum mechanical correction
						allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%alpha")

						atom%continua(kc)%alpha(:) = 0d0
        		!computes only where it is not zero
						atom%continua(kc)%alpha(:) = H_bf_Xsection(atom%continua(kc), lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))
						n_eff = atom%stage(atom%continua(kc)%j)*sqrt(atom%Rydberg/(atom%E(atom%continua(kc)%j)-atom%E(atom%continua(kc)%i)))

					else !interpolation of the read Cross-section
!       cont is not an alias but a copy of contiua(kc), so cont%alpha not deallocated
        				allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
        				if (alloc_status > 0) call ERROR("Allocation error cont%alpha")
         
!          atom%continua(kc)%alpha(:) = linear_1D(size(atom%continua(kc)%lambda_file),atom%continua(kc)%lambda_file,&
!             atom%continua(kc)%alpha_file, atom%continua(kc)%Nlambda,NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))
						call bezier3_interp(size(atom%continua(kc)%lambda_file), atom%continua(kc)%lambda_file, atom%continua(kc)%alpha_file, & !read values
     	atom%continua(kc)%Nlambda, lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), atom%continua(kc)%alpha) !interpolation grid

					endif
					
					if (atom%active) then
						n_cells_d = 1
						if (ldissolve) then
							if (atom%ID=="H") n_cells_d = n_cells
						endif
						allocate(atom%continua(kc)%alpha_nu(atom%continua(kc)%Nr-atom%continua(kc)%Nb+1,n_cells_d))
						atom%continua(kc)%alpha_nu(:,:) = 0.0_dp
						
						do icell=1, n_cells_d
							call bezier3_interp(size(atom%continua(kc)%alpha), lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), atom%continua(kc)%alpha, &
     						size(atom%continua(kc)%alpha_nu(:,icell)), lambda(atom%continua(kc)%Nb:atom%continua(kc)%Nr), atom%continua(kc)%alpha_nu(:,icell))						
						enddo
						
					endif
 
         !enddo
				case default
					call Error("Transition type unknown", atom%at(k)%trtype)
				end select
			enddo
			atom => NULL()
		enddo

	RETURN
	end SUBROUTINE alloc_atom_quantities
  
	SUBROUTINE compute_atom_quantities(icell)
		integer, intent(in) :: icell
		type(AtomType), pointer :: atom
		integer :: n, k, kc, alloc_status, i, j, Nblue, Nred, la, icell_d
		real(kind=dp) :: vbroad, aL, cte, cte2, adamp, wi, wj, gij, neff, chi_ion
    
    
		do n=1,Natom
			atom => Atoms(n)%ptr_atom
			do k = 1, atom%Ntr   
     
				kc = atom%at(k)%ik 
        
				select case (atom%at(k)%trtype)
        
				case ('ATOMIC_LINE')
					vbroad = atom%vbroad(icell)
         
					i=atom%lines(kc)%i; j=atom%lines(kc)%j
					wj = 1.0; wi = 1.0
					if (ldissolve) then
						if (atom%ID=="H") then! .or. atom%ID=="He") then
							wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1.0))
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1.0))
						endif
					endif

 !There are duplicates space array to remove. We do not need aeff if not Thomson profile        
					if (atom%lines(kc)%voigt) then

						call Damping(icell, atom, kc, adamp)
          !no need to keep it if we keep the profile for shift or interp
          !But damping can depend on populations of atom and electrons so ...
! 						if (associated(Profile, Iprofile_cmf_to_obs)) &
						if (allocated(atom%lines(kc)%a)) atom%lines(kc)%a(icell) = adamp
          !line%u(:) = atom%lines(kc)%u / atom%vbroad(icell)
          !                             if full profile to be interpolated or shifter.
          !                             if we use another approx for Voigt, it is computed on the fly, not shifted nor interp
!keep it at the moment for all case
						!-> replaced by call of profile()
						!!atom%lines(kc)%phi(:,icell) = Voigt(atom%lines(kc)%Nlambda, adamp, atom%lines(kc)%u(:)/vbroad)

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
					!!else !Gaussian
!no need if no interp nor shift or if we use Thomson (but keep adamp for thomson)
						!!atom%lines(kc)%phi(:,icell) = exp(-(atom%lines(kc)%u(:)/atom%vbroad(icell))**2)
					endif
!futur test to see if we keep it or not depending on the method of calculation of profiles
					!!atom%lines(kc)%phi(:,icell) = atom%lines(kc)%phi(:,icell) / (SQRTPI * atom%vbroad(icell))
					atom%lines(kc)%phi(:,icell) = local_profile_i(atom%lines(kc),icell,atom%lines(kc)%Nlambda,lambda(atom%lines(kc)%Nblue:atom%lines(kc)%Nred), &
																							0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp)

!Include occupa, cannot be lower, since ni =  niwi and nj = njwj
					if ((wj/wi * atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij) <= 0 ) then
						write(*,*) atom%ID, " nuij (10^15 Hz)", 1d-6 * CLIGHT / atom%lines(kc)%lambda0, " icell=", icell
						write(*,*) " lambdaij (nm)", atom%lines(kc)%lambda0 
						call warning ("background line: ni < njgij")
						write(*,*) "i = ", i, " j = ", j
						write(*,*) "w(i) = ", wi, " w(j) = ", wj
          				write(*,*) "ni=", atom%n(i,icell), " nj=", atom%n(j,icell), " gij=", atom%lines(kc)%gij
          				write(*,*) "nstari=", atom%n(i,icell), " njstar=", atom%n(j,icell)
          				stop				
					endif
    
				case ("ATOMIC_CONTINUUM")
        
					i=atom%continua(kc)%i; j=atom%continua(kc)%j
					Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
					chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
					neff = atom%stage(j) * sqrt(atom%Rydberg / (atom%E(j) - atom%E(i)) )

					icell_d = 1
					wj = 1.0; wi = 1.0
					if (ldissolve) then
						if (atom%ID=="H") then! .or. atom%ID=="He") then
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(j)))
							icell_d = icell
						endif
						if (atom%active) then
						
							do la=1, atom%continua(kc)%Nr-atom%continua(kc)%Nb+1
								atom%continua(kc)%alpha_nu(la,icell_d) = atom%continua(kc)%alpha_nu(la,icell_d) * &
								D_i(icell_d, neff, real(atom%stage(i)), 1.0, lambda(atom%continua(kc)%Nb+la-1), atom%continua(kc)%lambda0, chi_ion)
							enddo
						endif
					endif
   
					!if (atom%nstar(j,icell) < tiny_dp .or. atom%nstar(i,icell) < tiny_dp) then
					if (atom%active) then
! 						if (atom%nstar(i,icell) < tiny_dp) then
! 							atom%continua(kc)%gij(:,icell) = 0d0
! 						else
! 							do la=1, atom%continua(kc)%Nlambda
! 								if (exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)) * atom%nstar(i,icell) <= tiny_dp) then
! 									atom%continua(kc)%gij(la, icell) = 0.0_dp
! 								else
! 									atom%continua(kc)%gij(la, icell) = exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)) * atom%nstar(i,icell)/atom%nstar(j,icell)
! 								endif
! 								
!      						enddo
! 						endif
						
						gij = atom%nstar(i,icell)/atom%nstar(j,icell)
					
						do la=1, atom%continua(kc)%Nlambda

! 							if (atom%n(i,icell) - atom%n(j, icell)*atom%continua(kc)%gij(la, icell) <= 0 ) then
							if (atom%n(i,icell) - atom%n(j,icell) * gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)) <=0 ) then

								write(*,*) atom%ID, " nu+ (10^15 Hz)", 1d-6 * CLIGHT / atom%continua(kc)%lambda0, " icell=", icell
								write(*,*) " lambda+ (nm)", atom%continua(kc)%lambda0
								call warning ("background cont: ni < njgij")
								write(*,*) "i = ", i, " j = ", j
								write(*,*) "w(i) = ", wi, " w(j) = ", wj
          						write(*,*) atom%n(i,icell), wi*atom%n(j, icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))!atom%continua(kc)%gij(la, icell)
          						write(*,*) atom%nstar(i,icell), wi*atom%nstar(j, icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))!atom%continua(kc)%gij(la, icell)
								exit
							endif
						
						enddo
					endif !active
				case default
					call Error("Transition type unknown", atom%at(k)%trtype)
				end select
			enddo
			atom => NULL()
		enddo  
  
  
	RETURN
	end SUBROUTINE compute_atom_quantities
	
	!metal_bf typically
	subroutine opacity_atom_bf_lte(icell)

		integer, intent(in)											:: icell
		integer														:: m, kr, kc, i, j, Nblue, Nred, la
		type (AtomType), pointer									:: atom
		real(kind=dp)												:: n_eff, wj, wi
		real(kind=dp)												:: Diss, chi_ion, gij


		do m=1,Npassiveatoms

			atom => PassiveAtoms(m)%ptr_atom

			do kc=atom%Ntr_line+1,atom%Ntr
				kr = atom%at(kc)%ik 

				i = atom%continua(kr)%i
				j = atom%continua(kr)%j 

				Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred
				wj = 1.0; wi = 1.0

				if (ldissolve) then
					if (atom%ID=="H") then
						n_eff = real(i,kind=dp)
						wi = wocc_n(icell, n_eff, real(atom%stage(i)), real(atom%stage(j)))
					else
						n_eff = atom%stage(j)*sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))
					endif
				endif
				
				chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))

				do la=1, atom%continua(kr)%nlambda
				
					Diss = D_i(icell, n_eff, real(atom%stage(i)), 1.0, lambda_cont(Nblue+la-1), atom%continua(kc)%lambda0, chi_ion)
					gij = atom%nstar(i,icell)/atom%nstar(j,icell) * exp(-hc_k/T(icell)/lambda_cont(Nblue+la-1))

					
					if ( (atom%n(i,icell) - atom%n(j,icell) * gij) > 0.0) then
						
							chi_c(Nblue+la-1,icell) = chi_c(Nblue+la-1,icell) + &
							diss * atom%continua(kr)%alpha(la) * (atom%n(i,icell) - gij*atom%n(j,icell)) 
       				
							eta_c(Nblue+la-1,icell) = eta_c(Nblue+la-1,icell) + &
							diss * atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * gij * atom%n(j,icell)
					else
					
! 							eta_c(Nblue+la-1,icell) = eta_c(Nblue+la-1,icell) + &
! 							diss * atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * gij * atom%n(j,icell)   
							
							chi_c(Nblue+la-1,icell) = chi_c(Nblue+la-1,icell) - &
							diss * atom%continua(kr)%alpha(la) * (atom%n(i,icell) - gij*atom%n(j,icell)) 
							eta_c(Nblue+la-1,icell) = eta_c(Nblue+la-1,icell) + 0.0_dp

					endif
				enddo				

			end do ! loop over Ncont
			atom => NULL()
		end do !loop over metals

	return
	end subroutine opacity_atom_bf_lte
	
	subroutine background_continua (icell)
		integer, intent(in) :: icell
		integer :: la, nat
		real(kind=dp), dimension(Nlambda_cont) :: chi, eta, Bp

		Bp = Bpnu(T(icell), lambda_cont)

		chi_c(:,icell) = 0.0_dp
		eta_c(:,icell) = 0.0_dp
		sca_c(:,icell) = 0.0_dp

		call Hminus_bf_wishart(icell, Nlambda_cont, lambda_cont, chi, eta)
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + eta(:)

		call Hminus_ff_bell_berr(icell, Nlambda_cont, lambda_cont, chi)
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + chi(:) * Bp(:)
 
		call Hydrogen_ff(icell, Nlambda_cont, lambda_cont, chi)
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + chi(:) * Bp(:)

		!now atomic LTE bound-free
		call lte_bound_free(icell, Nlambda_cont, lambda_cont, chi, eta)
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + eta(:)

		sca_c(:,icell) = Thomson(icell)
		!Rayleigh scattering is included no more at the moment. Might be necessary for some temperature
		!below the Lyman limit

		!Total opac once source functions are known
		chi_c(:,icell) = chi_c(:,icell) +  sca_c(:,icell)

	return
	end subroutine background_continua
	
	subroutine compute_background_continua()
	!$ use omp_lib
		integer :: icell, id, icell0
		icell0 = 1
		!$omp parallel &
		!$omp default(none) &
		!$omp private(icell,id,icell0) &
		!$omp shared(n_cells, icompute_atomRT)
		!$omp do schedule(dynamic,1)
		do icell=1, n_cells
			!$ id = omp_get_thread_num() + 1
			if (icompute_atomRT(icell) > 0) then
				call compute_atom_quantities(icell) 
				!!need for BackgroundContinua
    			!!and lines 
				call background_continua(icell)	
			!else
				!Nothing to do here, opacity is zero		
			endif
		end do
		!$omp end do
		!$omp end parallel

	return
	end subroutine compute_background_continua
	
	subroutine interp_contopac()
	!$ use omp_lib
		integer :: icell, id, icell0
		icell0 = 1

		!$omp parallel &
		!$omp default(none) &
		!$omp private(icell,id,icell0) &
		!$omp shared(n_cells, icompute_atomRT,chi0_bb, eta0_bb)
		!$omp do schedule(dynamic,1)
		do icell=1, n_cells
			!$ id = omp_get_thread_num() + 1
			if (icompute_atomRT(icell) > 0) then			
				call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
				if (any_nan_infinity_vector(chi0_bb(:,icell)) /= 0) then
					write(*,*) chi0_bb(:,icell)
					call error("nan.infinity error in chi0_bb")
				endif
				if (any_nan_infinity_vector(eta0_bb(:,icell)) /= 0) then
					write(*,*) eta0_bb(:,icell)
					call error("nan.infinity error in eta0_bb")
				endif
			endif
		end do
		!$omp end do
		!$omp end parallel

	return
	end subroutine interp_contopac

		
	!scattering term is not taken as the one the continuum
	subroutine continuum_line (icell, lambda0, chi00, eta00)
		integer, intent(in) :: icell
		real(kind=dp), intent(in) :: lambda0
		real(kind=dp), intent(out) :: chi00, eta00
		integer :: la, nat
		real(kind=dp) :: chi(1), eta(1), l(1), Bp

		chi00 = 0.0
		eta00 = 0.0
		Bp = Bpnu(T(icell), lambda0)
		l(1) = lambda0

  
		call Hminus_bf_wishart(icell, 1, l, chi, eta)
		chi00 = chi00 + chi(1)
		eta00 = eta00 + eta(1)

		call Hminus_ff_bell_berr(icell, 1, l, chi)
		chi00 = chi00 + chi(1)
		eta00 = eta00 + chi(1) * Bp
 
		call Hydrogen_ff(icell, 1, l, chi)
		chi00 = chi00 + chi(1)
		eta00 = eta00 + chi(1) * Bp

		!now atomic LTE bound-free
		call lte_bound_free(icell, 1, l, chi, eta)
		chi00 = chi00 + chi(1)
		eta00 = eta00 + eta(1)

		chi00 = chi00 + Thomson(icell)
		!if using electron scattering from the continuum
		!if (lelectron_scattering) eta00 = eta00 + interp_dp(Jnu_cont(:,icell) * sca_c(:,icell), lambda_cont, lambda0) 
		
		chi00 = chi00 + interp_dp(chi_c_nlte(:,icell), lambda_cont, lambda0)
		eta00 = eta00 + interp_dp(eta_c_nlte(:,icell), lambda_cont, lambda0)

	return
	end subroutine continuum_line
   
	subroutine NLTE_bound_free(icell)
		integer, intent(in) :: icell
		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la
		type(AtomType), pointer :: aatom
		real(kind=dp) :: wi, wj, chi_ion, Diss, nn, gij
  
		chi_c_nlte(:,icell) = 0d0
		eta_c_nlte(:,icell) = 0d0

		atom_loop : do nact = 1, Nactiveatoms
			aatom => ActiveAtoms(nact)%ptr_atom
   
			tr_loop : do kr = aatom%Ntr_line+1,aatom%Ntr


        		kc = aatom%at(kr)%ik !relative index of a transition among continua or lines


				Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
				i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
				nn = n_eff(aatom%Rydberg, aatom%E(j), aatom%E(i), aatom%stage(j))!real(atom%stage(j)) * sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))

				
				wj = 1.0; wi = 1.0
				if (ldissolve) then
					if (aatom%ID=="H") then
												!nn
						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(j)))
					endif
				endif
     
				!get the ionisation potential for the ion ot be use in the dissolve fraction
				chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))

          		
         		do la=1,aatom%continua(kc)%Nlambda
				
					Diss = D_i(icell, nn, real(aatom%stage(i)), 1.0, lambda_cont(Nblue+la-1), aatom%continua(kc)%lambda0, chi_ion)
					gij = aatom%nstar(i,icell)/aatom%nstar(j,icell) * exp(-hc_k/T(icell)/lambda_cont(Nblue+la-1))

					if ((aatom%n(i,icell) - aatom%n(j,icell) * gij) > 0.0_dp) then
					
						chi_c_nlte(Nblue+la-1,icell) = chi_c_nlte(Nblue+la-1,icell) + &
            			Diss * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - gij * aatom%n(j,icell))
            
						eta_c_nlte(Nblue+la-1,icell) = eta_c_nlte(Nblue+la-1,icell) + &
            			Diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * gij * aatom%n(j,icell)

            
					else !neg or null
					
! 						eta_c_nlte(Nblue+la-1,icell) = eta_c_nlte(Nblue+la-1,icell) + &
!             			Diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * gij * aatom%n(j,icell)
 			
					!assuming small inversions
						eta_c_nlte(Nblue+la-1,icell) = eta_c_nlte(Nblue+la-1,icell) + 0.0_dp
						chi_c_nlte(Nblue+la-1,icell) = chi_c_nlte(Nblue+la-1,icell) - Diss * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - gij * aatom%n(j,icell))
 
					endif
          
				enddo   

			end do tr_loop
   


			aatom => NULL()
		end do atom_loop
	
 
	return
	end subroutine NLTE_bound_free
 
	subroutine compute_nlte_bound_free
	!$ use omp_lib
		integer :: icell, id


		!$omp parallel &
		!$omp default(none) &
		!$omp private(icell,id) &
		!$omp shared(n_cells,icompute_atomRT)
		!$omp do schedule(dynamic,1)
		do icell=1,n_cells
			!$ id = omp_get_thread_num() + 1
			if (icompute_atomRT(icell) > 0) then
				call NLTE_bound_free(icell)
			endif
		end do
		!$omp end do
		!$omp end parallel
 
	return
	end subroutine compute_nlte_bound_free
	
	subroutine interp_background_opacity(icell, chii, etai)
	!linear interpolation of continuum opacities from lambda_cont grid to lambda grid
	!lambda >> lambda_cont
		integer, intent(in) :: icell
		integer :: i, j, i0,j0, kc
		!it seems I'm having an interpolation error at the first point, so just avoinding it..
		real(kind=dp), dimension(Nlambda_cont) :: chi0, eta0
		real(kind=dp), dimension(Nlambda), intent(out) :: chii, etai
		real(kind=dp) :: u
		
		chii = 0.0_dp
		etai = 0.0_dp
     
     	chi0 = chi_c(:,icell)
     	eta0 = eta_c(:,icell)

!-> I use the table eta_es that is proportional to Jnu(Itot) instead
!      	if (lelectron_scattering) then
!      		eta0 = eta0 + sca_c(:,icell) * Jnu_cont(:,icell)
!      	endif
     	
     	if (NactiveAtoms>0) then
     		chi0 = chi0 + chi_c_nlte(:,icell)
     		eta0 = eta0 + eta_c_nlte(:,icell)
     	endif
     	
     	i0 = 1
     	j0 = 1

		!linear interpolation of both arrays
     	do i=i0,Nlambda_cont-1
        	do j=j0,Nlambda
           		if (lambda(j)>=lambda_cont(i) .and. lambda(j)<=lambda_cont(i+1)) then
            		u = (lambda(j) - lambda_cont(i)) / (lambda_cont(i+1) - lambda_cont(i))
!             		if (j<=3 .and. i<=3)write(*,*) "---"
!             		if (j<=3 .and. i<=3) write(*,*) icell, i,j, lambda_cont(i), lambda(j), " u=",u
              		chii(j) = (1.0_dp - u) * chi0(i) + u * chi0(i+1)
              		etai(j) = (1.0_dp - u) * eta0(i) + u * eta0(i+1)
!             		if (j<=3 .and. i<=3)write(*,*) "chii=",chii(j), etai(j)
!             		if (j<=3 .and. i<=3)write(*,*) "chi0i=",chi0(i), eta0(i)
!             		if (j<=3 .and. i<=3)write(*,*) "chi0i1=",chi0(i+1), eta0(i+1)
!             		if (j<=3 .and. i<=3)write(*,*) "---"
           		endif
        	enddo
     	enddo
     	
! 		j0=Nlambda+1
! 		do j=1, Nlambda
!         	if (lambda(j) > lambda_cont(1)) then
! 				j0 = j
!            		exit
!         	endif
!      	enddo
! 
! 		do j=j0, Nlambda
!         loop_i : do i=i0, Nlambda_cont
!         	if (lambda(i) > lambda_cont(j)) then
!             	u = (lambda(j) - lambda_cont(i-1)) / (lambda_cont(i) - lambda_cont(i-1))
!              	chii(j) = (1.0_dp - u) * chi0(i-1)  + u * eta0(i)
!              	etai(j) = (1.0_dp - u) * eta0(i-1)  + u * eta0(i)
!               	i0 = i
!               exit loop_i
!            endif
!         enddo loop_i
!      enddo

	return
	end subroutine interp_background_opacity

	subroutine opacity_atom_loc(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)

		integer, intent(in) :: id, icell, iray
		logical, intent(in) :: iterate
		real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la, Nvspace
		type(AtomType), pointer :: aatom
! 		integer :: nv
! 		integer, parameter 										:: NvspaceMax = 1
! 		real(kind=dp), dimension(NvspaceMax)					:: Omegav, dkv
! 		real(kind=dp) :: v0, v1, delta_vol_phi, xphi, yphi, zphi
		integer :: dk0
		real(kind=dp) :: wi, wj, nn, v0
		!tmp
		real(kind=dp), dimension(Nlambda) :: phi0

!-> chi and eta contains already chi_c, chi_c_nlte and eta_c, eta_c_nlte

!-> using dk instead of a full profile ATM
  
!   In principle Nvspace should depend on atom because of atom%vbroad(icell)
!   v_proj in m/s at point icell
! 		Omegav = 0d0
! 		Nvspace = 1
		v0 = v_proj(icell,x,y,z,u,v,w)
! 		dkv(1) = nint(1e-3 * v0 / hv)
! 		omegav(1) = v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w)
! 		v1 = v_proj(icell,x1,y1,z1,u,v,w)
! 		dkv(NvspaceMax) = nint(1e-3 * v1 / hv)
! 		dk0 = nint(1e-3 * v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w) / hv)
!		dk0 = dk(iray,id)
! 		write(*,*) 1, 1e-3 * v0, hv, nint(1e-3 * v0 / hv)
! 		do nv=2, NvspaceMax-1
!       		delta_vol_phi = (real(nv,kind=dp))/(real(NvspaceMax,kind=dp)) * l
! 			xphi=x+delta_vol_phi*u
! 			yphi=y+delta_vol_phi*v
! 			zphi=z+delta_vol_phi*w
! 			write(*,*) nv, 1e-3 * v_proj(icell,xphi,yphi,zphi,u,v,w), hv, nint(1e-3 * v_proj(icell,xphi,yphi,zphi,u,v,w) / hv)
! 			dkv(nv) = nint(1e-3 * v_proj(icell,xphi,yphi,zphi,u,v,w) / hv)
! 		enddo
! 		write(*,*) NvspaceMax, 1e-3 * v1, hv, nint(1e-3 * v1 / hv)
! 		stop
  
		atom_loop : do nact = 1, Natom!Nactiveatoms
			aatom => Atoms(nact)%ptr_atom!ActiveAtoms(nact)%ptr_atom

			tr_loop : do kr = 1,aatom%Ntr_line

				kc = aatom%at(kr)%ik !relative index of a transition among continua or lines

				Nred = aatom%lines(kc)%Nred; Nblue = aatom%lines(kc)%Nblue
				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
				
 				wj = 1.0; wi = 1.0
				if (ldissolve) then
					if (aatom%ID=="H") then
												!nn
						wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1)) !1 for H
						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1))
					endif
				endif 
				
! 				write(*,*) aatom%n(i,icell), aatom%n(j,icell), aatom%lines(kc)%gij
! 				write(*,*) aatom%lines(kc)%Aji, aatom%lines(kc)%Bij
! 				write(*,*) Nblue, Nred, Nblue+dk_min, Nred+dk_max
! 				write(*,*) lambda(Nblue), lambda(Nred), lambda(Nblue+dk_min), lambda(Nred+dk_max)
				
				!->not needed because phi is define only on its bound in local_profile_i
				!phi0(:) = 0.0_dp
							
! if (v0 <= hv) then
				phi0(Nblue+dk_min:Nred+dk_max) = local_profile_i(aatom%lines(kc),icell,Nred+dk_max-dk_min-Nblue+1,lambda(Nblue+dk_min:Nred+dk_max), x,y,z,x1,y1,z1,u,v,w,l)
! 				phi0(Nblue+dk_min:Nred+dk_max) = local_profile_interp(aatom%lines(kc),icell,Nred+dk_max-dk_min-Nblue+1,lambda(Nblue+dk_min:Nred+dk_max), x,y,z,x1,y1,z1,u,v,w,l)
! write(*,*) phi0(Nblue), phi0(Nred), phi0(Nblue+dk_min), phi0(Nred+dk_max)
! stop
! if (maxval(phi0) <= 0.0) then
! write(*,*) "line", kc, id, icell, Nred+dk_max-dk_min-Nblue+1, l
! write(*,*) "phi=", phi0(Nblue+dk_min:Nred+dk_max)
! write(*,*) "lam=", lambda(Nblue+dk_min:Nred+dk_max)
! stop
! endif

				if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then


					chi(Nblue+dk_min:Nred+dk_max,id) = chi(Nblue+dk_min:Nred+dk_max,id) + &
						hc_fourPI * aatom%lines(kc)%Bij * phi0(Nblue+dk_min:Nred+dk_max) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
! 					if (icell==66 .and. kc==2) then
! 						do v0=1, Nred+dk_max-dk_min-Nblue+1
! 							write(*,*) lambda(Nblue+dk_min+v0-1), hc_fourPI * aatom%lines(kc)%Bij * phi0(Nblue+dk_min+v0-1) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell)),hc_fourPI * aatom%lines(kc)%Aji * phi0(Nblue+dk_min+v0-1) * aatom%n(j,icell)  
! 						enddo 						
! 					stop
! 					endif
					eta(Nblue+dk_min:Nred+dk_max,id)= eta(Nblue+dk_min:Nred+dk_max,id) + &
						hc_fourPI * aatom%lines(kc)%Aji * phi0(Nblue+dk_min:Nred+dk_max) * aatom%n(j,icell)
! else
! 					chi(Nblue+dk0:Nred+dk0,id) = chi(Nblue+dk0:Nred+dk0,id) + &
! 					hc_fourPI * aatom%lines(kc)%Bij * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell)) * &
! 					aatom%lines(kc)%phi(:,icell)
! 					       	 
! 					eta(Nblue+dk0:Nred+dk0,id)= eta(Nblue+dk0:Nred+dk0,id) + &
! 					hc_fourPI * aatom%lines(kc)%Aji * aatom%n(j,icell) * &
! 					aatom%lines(kc)%phi(:,icell)
! endif
				else !neg or null
! 					eta(Nblue+dk_min:Nred+dk_max,id)= eta(Nblue+dk_min:Nred+dk_max,id) + &
! 						hc_fourPI * aatom%lines(kc)%Aji * aatom%n(j,icell) * phi0(Nblue+dk_min:Nred+dk_max)
				!small inversions
					eta(Nblue+dk_min:Nred+dk_max,id)= eta(Nblue+dk_min:Nred+dk_max,id) + 0.0_dp
					chi(Nblue+dk_min:Nred+dk_max,id) = chi(Nblue+dk_min:Nred+dk_max,id) - &
					hc_fourPI * aatom%lines(kc)%Bij * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell)) * phi0(Nblue+dk_min:Nred+dk_max)
				endif
				
				if (iterate) then
					aatom%lines(kc)%phi_loc(:,iray,id) = phi0(Nblue+dk_min:Nred+dk_max)
				endif

    
			end do tr_loop

			aatom => NULL()

		end do atom_loop
		
! 					do la=1, aatom%lines(kc)%Nlambda
! 						chi(Nblue+la-1+dk0,id) = chi(Nblue+la-1+dk0,id) + &
! 						hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
!        	 
! 						eta(Nblue+la-1+dk0,id)= eta(Nblue+la-1+dk0,id) + &
! 						hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 														
! 					enddo 
! 		if (any_nan_infinity_vector(eta(:,id)) /= 0) then
! 				write(*,*) eta(:,id)
! 				call error("nan.infinity error in opacity_atom_loc (eta)")
! 		endif
! 		if (any_nan_infinity_vector(chi(:,id)) /= 0) then
! 				write(*,*) chi(:,id)
! 				call error("nan.infinity error in opacity_atom_loc (chi)")
! 		endif
! 		write(*,*) "e=",eta(:,id)
! 		write(*,*) "c=",chi(:,id)
! 		stop

	return
	end subroutine opacity_atom_loc
	
	
end module Opacity

 
	
! 	subroutine opacity_atom_loc(id, icell, iray)
! 
! 		integer, intent(in) :: icell, iray, id
! 		!!real(kind=dp), intent(in) :: x0, y0, z0, x1, y1, z1, u, v, w, l
! 		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, Nvspace, dk0, lam
! 		type(AtomType), pointer :: aatom
! 		real(kind=dp) :: wi, wj, nn, chi_ion, Diss
! 		
! 		dk0 = dk(iray,id)
!   
! 		chi(:) = linear_1D(NLTEspec%Nwaves_cont, NLTEspec%lambda_cont, NLTEspec%Kc(:,icell) + NLTEspec%Kc_nlte(:,icell), NLTEspec%Nwaves, NLTEspec%lambda)
! 		!assumes sca_c * Jcont is the electron scattering emissivity
! 		eta(:) = linear_1D(NLTEspec%Nwaves_cont, NLTEspec%lambda_cont,  NLTEspec%Kc(:,icell) + NLTEspec%jc_nlte(:,icell), NLTEspec%Nwaves, NLTEspec%lambda)
! 	    
! 		atom_loop : do nact = 1, Natom
! 		aatom => Atoms(nact)%ptr_atom
! 
! 			tr_loop : do kr = 1, aatom%Ntr_line
! 			
! 				kc = aatom%at(kr)%ik 
! 
! 			
! ! 				select case (aatom%at(kr)%trtrype)
! ! 
! ! 				case ("ATOMIC LINE")
! 
! 					Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
! 					i = aatom%lines(kc)%i;j = aatom%lines(kc)%j				
! 				
!  					wj = 1.0; wi = 1.0
! 					if (ldissolve) then
! 						if (aatom%ID=="H") then
! 												!nn
! 							wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1)) !1 for H
! 							wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1))
! 						endif
! 					endif 
! 
! 
! 				
! 					if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then
! 				
! 
!         
! 						chi(Nblue+dk0:Nred+dk0) = chi(Nblue+dk0:Nred+dk0) + &
! 							hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(:,icell) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
!        	 
! 						eta(Nblue+dk0:Nred+dk0) = eta(Nblue+dk:Nred+dk) + &
! 							hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi(:,icell) * aatom%n(j,icell)
! 				
! 
! 				
! 					else !neg or null
! 				
! 
!        	 
! 						eta(Nblue+dk0:Nred+dk0) = eta(Nblue+dk0:Nred+dk0) + &
! 							hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi(:,icell) * aatom%n(j,icell)
! 				
! 				
! 
! 					endif
! ! 				case ("ATOMIC_CONTINUUM")
! ! 
! ! 					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
! ! 					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
! ! 					nn = n_eff(aatom%Rydberg, aatom%E(j), aatom%E(i), aatom%stage(j))!real(atom%stage(j)) * sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))
! ! 
! ! 				
! ! 					wj = 1.0; wi = 1.0
! ! 					if (ldissolve) then
! ! 						if (aatom%ID=="H") then
! ! 												!nn
! ! 							wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(j)))
! ! 						endif
! ! 					endif
! !      
! ! 				!get the ionisation potential for the ion ot be use in the dissolve fraction
! ! 					chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
! ! 
! !           		
! !           			do la=1, atom%continua(kc)%Nlambda
! ! 				
! ! 						Diss = D_i(icell, nn, real(aatom%stage(i)), 1.0, NLTEspec%lambda(la), aatom%continua(kc)%lambda0, chi_ion)
! ! 
! ! 						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0_dp) then
! ! 					
! ! 							chi(Nblue:Nred) = chi(Nblue:Nred) + &
! !             					Diss * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell))
! !             
! ! 							eta(Nblue:Nred) = eta(Nblue:Nred) + &
! !             					Diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
! ! 
! !             
! ! 						else !neg or null
! ! 					
! ! 							eta(Nblue:Nred) = eta(Nblue:Nred) + &
! !             					Diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
! ! 
! !  
! ! 						endif	
! ! 					enddo			
! 			    
! !     			case default
! !     				call error("transition unknown", atom%at(kr)%trtypr)
! !     			end select
!     			
! 			end do tr_loop
! 
! 			aatom => NULL()
! 
! 		end do atom_loop
! 
! 	return
! 	end subroutine opacity_atom_loc	
! 
! 
! 	!for all angle.
! 	!Total source fonction in the different radiative rates ?? For all atoms ??
! 	!Or the coupling with other atom should be only in chi and I ??
! 	subroutine calc_total_source_loc(id, icell, n_rayons)
! 	! ------------------------------------------------------------------------- !
! 	! ------------------------------------------------------------------------- !  
!   
! 		integer, intent(in) :: id, icell, n_rayons
! 		integer :: kr, kc, i, j, Nred, Nblue, la, nact
! 		!real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, dtau
! 		type (AtomType), pointer :: aatom
! 		integer :: dk, iray
! 		real(kind=dp) :: nn, wi, wj, chi_ion, dissolve
! 
! 		!LTE, unchanged opacities for local sub iterations
! 		do iray=1, n_rayons
! 			NLTEspec%chi(:,iray,id) = NLTEspec%AtomOpac%Kc(:,icell) 
! 			NLTEspec%S(:,iray,id) = NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%sca_c(:,icell)*NLTEspec%Jc(:,icell)
! 			if (allocated(chi_loc)) then
! 				NLTEspec%chi(:,iray,id) = NLTEspec%chi(:,iray,id) + chi_loc(:,iray,id)
! 				NLTEspec%S(:,iray,id) = NLTEspec%S(:,iray,id) + eta_loc(:,iray,id)
! 			endif
! 		enddo
! 
!   
! 		atom_loop : do nact = 1, Nactiveatoms
! 			aatom => ActiveAtoms(nact)%ptr_atom
! 			
! 
! 				tr_loopc : do kr = aatom%Ntr_line+1,aatom%Ntr
! 
! 
!         			kc = aatom%at(kr)%ik
! 
! 
! 					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
! 					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
! 					
! 					nn = n_eff(aatom%Rydberg, aatom%E(j), aatom%E(i), aatom%stage(j))!real(atom%stage(j)) * sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))
! 					chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
! 
! 				
! 					wj = 1.0; wi = 1.0
! 					if (ldissolve) then
! 						if (aatom%ID=="H") then
! 												!nn
! 							wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(j)))
! 						endif
! 					endif
! 
!           
!          			do la=1,aatom%continua(kc)%Nlambda
!          			
! 						dissolve = D_i(icell, nn, real(aatom%stage(i)), 1.0, NLTEspec%lambda(Nblue+la-1), aatom%continua(kc)%lambda0, chi_ion)
! 
! 						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0_dp) then
! 
! 
! 							NLTEspec%chi(Nblue+la-1,:,id) = NLTEspec%chi(Nblue+la-1,:,id) + dissolve * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell))
! 							NLTEspec%S(Nblue+la-1,:,id) = NLTEspec%S(Nblue+la-1,:,id) + dissolve * aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
! 							
! 						else !negative or null
! 
! 							NLTEspec%S(Nblue+la-1,:,id) = NLTEspec%S(Nblue+la-1,:,id) + dissolve * aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
! 
! 													
! 						endif
!             
!           
! 					enddo   
! 
! 
! 				end do tr_loopc
! 
!    
! 			tr_loop : do kr = 1, aatom%Ntr_line
! 
!    	
! 				kc = aatom%at(kr)%ik
! 
! 				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
! 				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
! 				
! 				wj = 1.0; wi = 1.0
! 				if (ldissolve) then
! 					if (aatom%ID=="H") then
! 												!nn
! 						wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1)) !1 for H
! 						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1))
! 					endif
! 				endif
! 
! 
!          		!NLTE subroutine. uses dk instead of a profile
! 				
! 				if ((wj/wi * aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then
! 
! 					do la=1, aatom%lines(kc)%Nlambda
! 
! 						do iray=1, n_rayons
! 							dk = aatom%lines(kr)%dk(iray,id)
! 
! 							NLTEspec%chi(Nblue+la-1+dk,iray,id) = NLTEspec%chi(Nblue+la-1+dk,iray,id) + &
! 							hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
!          
! 							NLTEspec%S(Nblue+la-1+dk,iray,id) = NLTEspec%S(Nblue+la-1+dk,iray,id) + &
! 							aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 						enddo
! 					enddo
! 					
! 				else !negative or null
! 				
! 					do la=1, aatom%lines(kc)%Nlambda
!          				do iray=1, n_rayons
!          					dk = aatom%lines(kr)%dk(iray,id)
! 
! 							NLTEspec%S(Nblue+la-1+dk,iray,id) = NLTEspec%S(Nblue+la-1+dk,iray,id) + &
! 							aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 						enddo	
! 					enddo
! 					
! 				endif
! 
!     
! 			end do tr_loop
! 
! 			aatom => NULL()
! 		end do atom_loop
! 		
! 		NLTEspec%S(:,:,id) = NLTEspec%S(:,:,id)/NLTEspec%chi(:,:,id)
! 
! 
! 	return
! 	end subroutine calc_total_source_loc
 	
 	!angle by angle
! 	subroutine calc_total_source_loc_o(id, icell, iray, initialize) !not used yet, avoid recomputing contopac for each dir
! 	! ------------------------------------------------------------------------- !
! 	! ------------------------------------------------------------------------- !  
!   
! 		integer, intent(in) :: id, icell, iray
! 		logical, intent(in) :: initialize
! 		integer :: kr, kc, i, j, Nred, Nblue, la, nact
! 		!real(kind=dp), dimension(NLTEspec%Nwaves) :: chi, eta, dtau
! 		type (AtomType), pointer :: aatom
! 		integer :: dk
! 		real(kind=dp) :: nn, wi, wj, chi_ion, dissolve
! 
! 		!LTE, unchanged opacities for local sub iterations
! 		NLTEspec%chi(:,iray,id) = NLTEspec%AtomOpac%Kc(:,icell) 
! 		!Jc should be updated or the old value should be used ? Jc is computed at the end of the cell loop
! 		!to avoid to update it
! 		NLTEspec%S(:,iray,id) = NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%sca_c(:,icell)*NLTEspec%Jc(:,icell)
! 		if (allocated(chi_loc)) then
! 			NLTEspec%chi(:,iray,id) = NLTEspec%chi(:,iray,id) + chi_loc(:,iray,id)
! 			NLTEspec%S(:,iray,id) = NLTEspec%S(:,iray,id) + eta_loc(:,iray,id)
! 		endif
! 
!   
! 		atom_loop : do nact = 1, Nactiveatoms
! 			aatom => ActiveAtoms(nact)%ptr_atom
! 			
! 			!before continua, because we need to add etac for all rays in eta.
! 			!if it is too long I reintroduce etac and chic to compute only for one ray. (see calc_eta_atom_loc)
! 			!if (initialize) then
! 				!aatom%chic(:,id) = 0.0 !local but we want to keep chic for all rays to add it to the total chi of each ray
! 			
! 				tr_loopc : do kr = aatom%Ntr_line+1,aatom%Ntr
! 
! 
!         			kc = aatom%at(kr)%ik
! 
! 
! 					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
! 					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
! 					
! 					nn = n_eff(aatom%Rydberg, aatom%E(j), aatom%E(i), aatom%stage(j))!real(atom%stage(j)) * sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))
! 					chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
! 
! 				
! 					wj = 1.0; wi = 1.0
! 					if (ldissolve) then
! 						if (aatom%ID=="H") then
! 												!nn
! 							wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(j)))
! 						endif
! 					endif
! 
!           
!          			do la=1,aatom%continua(kc)%Nlambda
!          			
! 						dissolve = D_i(icell, nn, real(aatom%stage(i)), 1.0, NLTEspec%lambda(Nblue+la-1), aatom%continua(kc)%lambda0, chi_ion)
! 
! 						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0_dp) then
! 
! 
! 							NLTEspec%chi(Nblue+la-1,iray,id) = NLTEspec%chi(Nblue+la-1,iray,id) + dissolve * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell))
! 							NLTEspec%S(Nblue+la-1,iray,id) = NLTEspec%S(Nblue+la-1,iray,id) + dissolve * aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
! 							
! 						else !negative or null
! 
! 							NLTEspec%S(Nblue+la-1,iray,id) = NLTEspec%S(Nblue+la-1,iray,id) + dissolve * aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
! 
! 													
! 						endif
!             
!           
! 					enddo   
! 
! 
! 				end do tr_loopc
! 
!    
! 			tr_loop : do kr = 1, aatom%Ntr_line
! 
!    	
! 				kc = aatom%at(kr)%ik
! 
! 				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
! 				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
! 				dk = aatom%lines(kr)%dk(iray,id)
! 				
! 				wj = 1.0; wi = 1.0
! 				if (ldissolve) then
! 					if (aatom%ID=="H") then
! 												!nn
! 						wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1)) !1 for H
! 						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1))
! 					endif
! 				endif
! 
! 
!          		!NLTE subroutine. uses dk instead of a profile
! 				
! 				if ((wj/wi * aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then
! 
! 					do la=1, aatom%lines(kc)%Nlambda
! 
! 						NLTEspec%chi(Nblue+la-1+dk,iray,id) = NLTEspec%chi(Nblue+la-1+dk,iray,id) + &
! 						hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
!          
! 						NLTEspec%S(Nblue+la-1+dk,iray,id) = NLTEspec%S(Nblue+la-1+dk,iray,id) + &
! 						aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 					enddo
! 					
! 				else !negative or null
! 				
! 					do la=1, aatom%lines(kc)%Nlambda
!          
! 						NLTEspec%S(Nblue+la-1+dk,iray,id) = NLTEspec%S(Nblue+la-1+dk,iray,id) + &
! 						aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 									
! 					enddo
! 					
! 				endif
! 
!     
! 			end do tr_loop
! 
! 			aatom => NULL()
! 		end do atom_loop
! 		
! 		NLTEspec%S(:,iray,id) = NLTEspec%S(:,iray,id)/NLTEspec%chi(:,iray,id)
! 
! 
! 	return
! 	end subroutine calc_total_source_loc_o
	
	!Not up to date, check with calc_total_source_loc
	!Done in the propagation, to have it at id (icell) for iray==1
! 	subroutine calc_etac_atom_loc(id, icell)
! 	!Store the continuous emissivity for the running cell, for all atoms
! 		integer, intent(in) :: icell, id
! 		integer :: nact, Nred, Nblue, kc, kr, i, j, nk, la
! 		type(AtomType), pointer :: aatom
!   
! 
! 		atom_loop : do nact = 1, Nactiveatoms
! 			aatom => ActiveAtoms(nact)%ptr_atom
!    
! 		!loop only on continua
! 			tr_loop : do kr = aatom%Ntr_line+1,aatom%Ntr
! 
! 
!         		kc = aatom%at(kr)%ik !relative index of a transition among continua or lines
! 
! 
! 				Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
! 				i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
! 
!           
!          		do la=1,aatom%continua(kc)%Nlambda
! 
! 					if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0) then
!             
! 						aatom%etac(Nblue+la-1,id) = aatom%etac(Nblue+la-1,id) + &
!             aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
!             
! 					endif
!           
! 				enddo   
! 
! 
! 			end do tr_loop
!    
! 
! 
! 			aatom => NULL()
! 		end do atom_loop
! 
! 		
! 	return
! 	end subroutine calc_etac_atom_loc

 
	!Not up to date, check with calc_total_source_loc
! 	subroutine calc_eta_atom_loc(id, icell, iray, initialize) !update eta and psi, iterate always true
! 	! ------------------------------------------------------------------------- !
! 	! ------------------------------------------------------------------------- !  
!   
! 		integer, intent(in) :: id, icell, iray
! 		logical, intent(in) :: initialize
! 		integer :: kr, kc, i, j, Nred, Nblue, la, nact
! 		real(kind=dp), dimension(NLTEspec%Nwaves) :: dtau, chi
! 		type (AtomType), pointer :: aatom
! 		integer :: dk
! 
! 		!LTE, unchanged opacities for local sub iterations
! 		chi(:) = chi_loc(:,iray,id)
!   
! 		atom_loop : do nact = 1, Nactiveatoms
! 			aatom => ActiveAtoms(nact)%ptr_atom
! 			
! 			!before continua, because we need to add etac for all rays in eta.
! 			if (initialize) then
! 				aatom%chic(:,id) = 0.0 !local but we want to keep chic for all rays to add it to the total chi of each ray
! 			
! 				tr_loopc : do kr = aatom%Ntr_line+1,aatom%Ntr
! 
! 
!         			kc = aatom%at(kr)%ik
! 
! 
! 					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
! 					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
! 
!           
!          			do la=1,aatom%continua(kc)%Nlambda
! 
! 						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0) then
! 
! 							aatom%chic(Nblue+la-1,id) = aatom%chic(Nblue+la-1,id)+  aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell)*aatom%n(j,icell))
! 							aatom%etac(Nblue+la-1,id) = aatom%etac(Nblue+la-1, id) + aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell)
! 						
! 						endif
!             
!           
! 					enddo   
! 
! 
! 				end do tr_loopc
! 						
! 			endif
! 			
! 			!etac is zero at iray==1, computed and then added for each ray (at this point, non zero)
! 			aatom%eta(:,iray,id) = aatom%etac(:,id)
! 			chi(:) = chi(:) + aatom%chic(:,id) !constant over direction, added to chi_loc(:,iray,id)
!    
! 			tr_loop : do kr = 1, aatom%Ntr_line
! 
!    	
! 				kc = aatom%at(kr)%ik
! 
! 				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
! 				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
! 				dk = aatom%lines(kr)%dk(iray,id)
! 
! 				if ((aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0) cycle
! 
!          		!NLTE subroutine. uses dk instead of a profile
! 				do la=1, aatom%lines(kc)%Nlambda
! 
! 					chi(Nblue+la-1-dk) = chi(Nblue+la-1-dk) + &
! 					hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell) - aatom%lines(kc)%gij*aatom%n(j,icell))
!          
! 					aatom%eta(Nblue+la-1-dk,iray,id) = aatom%eta(Nblue+la-1-dk,iray,id) + &
! 					aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 
! 				enddo
!     
! 			end do tr_loop
! 
! 			aatom => NULL()
! 		end do atom_loop
! 		
! 		! Check that chi_loc is correct
!         !opacity total = LTE (doesn't change during sub iterations) + NLTE
! 		dtau = chi(:) * ds(iray,id)
! 
! 		call calc_psi_operator_m(id, icell, iray, chi, dtau)
! 
! 
! 	return
! 	end subroutine calc_eta_atom_loc

	!? total absorption but only total emissivity per atom for SEE ? Or total Source function
	!at a given frequency ??
! 	subroutine opacity_atom_loc(id, icell, n_rayons)
! 	! ------------------------------------------------------------------------- !
! 	! ------------------------------------------------------------------------- !  
!   
! 		integer, intent(in) :: id, icell, n_rayons
! 		integer :: kr, kc, i, j, Nred, Nblue, la, nact
! 		type (AtomType), pointer :: aatom
! 		integer :: dk, iray
! 		real(kind=dp) :: nn, wi, wj, chi_ion, dissolve
! 
! 		NLTEspec%chi(:,:,id) = 0.0
! 		do iray=1, n_rayons
! 			NLTEspec%chi(:,iray,id) = NLTEspec%AtomOpac%Kc(:,icell) 
! 			if (allocated(chi_loc)) then
! 				NLTEspec%chi(:,iray,id) = NLTEspec%chi(:,iray,id) + chi_loc(:,iray,id)
! 			endif
! 		enddo
!   
! 		atom_loop : do nact = 1, Nactiveatoms
! 			aatom => ActiveAtoms(nact)%ptr_atom
! 			aatom%eta(:,:,id) = 0.0
! 
! 			do iray=1, n_rayons
! 				aatom%eta(:,iray,id) = aatom%eta(:,iray,id) +  NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%sca_c(:,icell)*NLTEspec%Jc(:,icell)
! 				if (allocated(eta_loc)) aatom%eta(:,iray,id) = aatom%eta(:,iray,id) + eta_loc(:,iray,id)
! 			enddo
! 
! 				tr_loopc : do kr = aatom%Ntr_line+1,aatom%Ntr
! 
! 
!         			kc = aatom%at(kr)%ik
! 
! 
! 					Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
! 					i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
! 					
! 					nn = n_eff(aatom%Rydberg, aatom%E(j), aatom%E(i), aatom%stage(j))!real(atom%stage(j)) * sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))
! 					chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
! 
! 				
! 					wj = 1.0; wi = 1.0
! 					if (ldissolve) then
! 						if (aatom%ID=="H") then
! 												!nn
! 							wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(j)))
! 						endif
! 					endif
! 
!           
!          			do la=1,aatom%continua(kc)%Nlambda
!          			
! 						dissolve = D_i(icell, nn, real(aatom%stage(i)), 1.0, NLTEspec%lambda(Nblue+la-1), aatom%continua(kc)%lambda0, chi_ion)
! 
! 						if ((aatom%n(i,icell) - aatom%n(j,icell) * aatom%continua(kc)%gij(la,icell)) > 0.0_dp) then
! 
! 							NLTEspec%chi(Nblue+la-1,:,id) = NLTEspec%chi(Nblue+la-1,:,id) + dissolve * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - aatom%continua(kc)%gij(la,icell) * aatom%n(j,icell))
! 							aatom%eta(Nblue+la-1,:,id) = aatom%eta(Nblue+la-1,:,id) + dissolve * aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
! 							
! 						else !negative or null
! 
! 							aatom%eta(Nblue+la-1,:,id) = aatom%eta(Nblue+la-1,:,id) + dissolve * aatom%continua(kc)%alpha(la) * aatom%n(j,icell) * aatom%continua(kc)%twohnu3_c2(la) * aatom%continua(kc)%gij(la,icell)
! 													
! 						endif
!             
!           
! 					enddo   
! 
! 
! 				end do tr_loopc
! 						
!    
! 			tr_loop : do kr = 1, aatom%Ntr_line
! 
!    	
! 				kc = aatom%at(kr)%ik
! 
! 				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
! 				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j
! 				
! 				wj = 1.0; wi = 1.0
! 				if (ldissolve) then
! 					if (aatom%ID=="H") then
! 												!nn
! 						wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1)) !1 for H
! 						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1))
! 					endif
! 				endif
! 				
! 
!          		!NLTE subroutine. uses dk instead of a profile
! 				
! 				if ((wj/wi * aatom%n(i,icell) - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then
! 
! 					do la=1, aatom%lines(kc)%Nlambda
! 					
! 						do iray=1, n_rayons
! 							dk = aatom%lines(kr)%dk(iray,id)
! 
! 
! 							NLTEspec%chi(Nblue+la-1-dk,iray,id) = NLTEspec%chi(Nblue+la-1-dk,iray,id) + &
! 							hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
!          
! 							aatom%eta(Nblue+la-1-dk,iray,id) = aatom%eta(Nblue+la-1-dk,iray,id) + &
! 							aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 						
! 						enddo
! 					enddo
! 					
! 				else !negative or null
! 				
! 					do la=1, aatom%lines(kc)%Nlambda
! 						do iray=1, n_rayons
!          
! 							aatom%eta(Nblue+la-1-dk,iray,id) = aatom%eta(Nblue+la-1-dk,iray,id) + &
! 							wi*aatom%lines(kc)%twohnu3_c2 * aatom%lines(kc)%gij * hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi(la,icell) * aatom%n(j,icell)
! 						enddo
! 					enddo
! 					
! 				endif
! 
!     
! 			end do tr_loop
! 
! 			aatom => NULL()
! 		end do atom_loop
! 
! 
! 	return
! 	end subroutine opacity_atom_loc

! 	subroutine calc_psi_operator_m(id, icell, iray, chi, dtau)
! 	!Mali
! 	! I' = (I - eta) Psi + Psi* eta with constant background opacities during succesive iterations
! 
! 		integer, intent(in) :: iray, id, icell
! 		real(kind=dp), dimension(NLTEspec%Nwaves), intent(in) :: chi, dtau
!    
! 		NLTEspec%Psi(:,iray,id) = (1d0 - exp(-dtau)) / chi
!    
! 		NLTEspec%etau(:,iray,id) = exp(-dtau)
! 
! 
! 	return
! 	end subroutine calc_psi_operator_m
! 	
! 	
! 	subroutine calc_psi_operator(id, icell, iray, chi, dtau, S)
! 	
! 	! I' = I * etau + S * Psi
! 
! 		integer, intent(in) :: iray, id, icell
! 		real(kind=dp), dimension(NLTEspec%Nwaves), intent(in) :: chi, dtau, S
!    
! 		NLTEspec%Psi(:,iray,id) = (1d0 - exp(-dtau)) 
!    
! 		NLTEspec%etau(:,iray,id) = exp(-dtau)
! 		
! 		NLTEspec%S(:,iray,id) = S(:)
! 
! 
! 	return
! 	end subroutine calc_psi_operator
