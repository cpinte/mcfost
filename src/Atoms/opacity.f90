module Opacity

	use atmos_type, only								: hydrogen, ntotal_atom, T, ne, NactiveAtoms, Natom, Npassiveatoms, ActiveAtoms, PassiveAtoms, Atoms, Elements, &
															icompute_atomRT, lmali_scheme, lhogerheijde_scheme, nhtot
	use atom_type
	use spectrum_type, only								: Itot, Icont, lambda, nlambda, Nlambda_cont, lambda_cont, dk, dk_min, dk_max, &
															chi0_bb, eta0_bb, chi_c, sca_c, eta_c, chi, eta, chi_c_nlte, eta_c_nlte, Jnu_cont, &
															rho_p, etaQUV_p, chiQUV_p
	use constant
	use constantes, only								: tiny_dp, huge_dp, AU_to_m
	use messages
	use broad, only										: Damping, line_damping
	use parametres
	use voigtfunctions, only							: Voigt
! 	use stark, only										: Stark_profile
	use profiles!, only									: Profile
	use getlambda, only									: hv, Nlambda_max_line, define_local_profile_grid
	use planck, only 									: bpnu
	use occupation_probability, only					: wocc_n, D_i
! 	use lte, only 										: phi_T

	use background_opacity, only						: Thomson, Hydrogen_ff, Hminus_bf, Hminus_bf_geltman, &
															Hminus_bf_geltman, Hminus_bf_Wishart, Hminus_ff_john, Hminus_ff, Hminus_ff_bell_berr, lte_bound_free, H_bf_Xsection
	use mcfost_env, only								: dp
	use input, only										: ds
	use molecular_emission, only						: v_proj


	implicit none
	real(kind=dp), parameter :: prec_pops = 1d-100, frac_ne_limit = 1d-10, frac_ntot_limit = 1d-10!1d-50!1d-15
	real(kind=dp) :: frac_limit_pops = frac_ntot_limit
    !for one ray
	real(kind=dp), allocatable :: eta_atoms(:,:,:), Uji_down(:,:,:,:), chi_up(:,:,:,:), chi_down(:,:,:,:)
	real(kind=dp), allocatable :: R_xcc(:,:,:,:)
 
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
!-> not allocated at the moment
! 							if (associated(Profile,local_profile_thomson)) then
! 								deallocate(atom%lines(kc)%eta)
! 							endif
          

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
   !To Do: Tex and Tion from populations instead of T(:). Allows to restard from previous calc
	subroutine alloc_atom_quantities
  
		type(AtomType), pointer :: atom
		integer :: n, k, kc, j, alloc_status, la, n_cells_d, icell, Nlam
		real(kind=dp) :: n_eff, u
		integer(kind=8) :: size_phi_all = 0, mem_alloc_local = 0

    
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
						allocate(atom%lines(kc)%Rij(nb_proc), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%Rij")
						atom%lines(kc)%Rij(:) = 0d0
						allocate(atom%lines(kc)%Rji(nb_proc), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error line%Rji")
						atom%lines(kc)%Rji(:) = 0d0
					endif
					
					mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%Tex) + sizeof(atom%lines(kc)%Rij)*2
          

          
			
					if (atom%lines(kc)%voigt) then
!-> do not keep adamp if profile is interpolated ! Profile are reevaluated					
							if (associated(profile, local_profile_v)) then
								allocate(atom%lines(kc)%a(n_cells), stat=alloc_status)
								if (alloc_status > 0) call ERROR("Allocation error line%adamp")
								atom%lines(kc)%a(:) = 0d0
								mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%a(:))
								if (allocated(atom%lines(kc)%u)) deallocate(atom%lines(kc)%u)
								!-> it is allocated with the definition of the line bounds, but not needed if we don't interpolate the profile
								
							elseif (associated(profile,local_profile_interp)) then
								allocate(atom%lines(kc)%phi(size(atom%lines(kc)%u), n_cells), stat=alloc_status)
								if (alloc_status > 0) call ERROR("Allocation error line%phi")
								atom%lines(kc)%phi(:,:) = 0d0
								size_phi_all = size_phi_all + sizeof(atom%lines(kc)%phi(:,:)) + sizeof(atom%lines(kc)%u(:))
								mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%phi(:,:)) + sizeof(atom%lines(kc)%u(:))
								
							elseif (associated(profile,local_profile_dk)) then
								!doesn't include the vel shift, just a local profile evaluated on the nlte grid.
								Nlam = atom%lines(kc)%Nlambda
								if (allocated(atom%lines(kc)%u)) deallocate(atom%lines(kc)%u)
								if (allocated(atom%lines(kc)%u)) call error("line%u should not be allocated here!")

								allocate(atom%lines(kc)%phi(Nlam, n_cells), stat=alloc_status)

								if (alloc_status > 0) call ERROR("Allocation error line%phi")
								atom%lines(kc)%phi(:,:) = 0d0
								size_phi_all = size_phi_all + sizeof(atom%lines(kc)%phi(:,:))
								mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%phi(:,:))
							
!-> not allocated anymore atm
!but in case, allocate only eta(n_cells) and creates a function to recompute aeff just like damping are possibly re evaluated
!if profile_v, we need to possibly re evaluate eta(n_cells; adamp).
! 							elseif (associated(Profile, local_profile_thomson)) then
! 								allocate(atom%lines(kc)%a(n_cells), stat=alloc_status)
! 								if (alloc_status > 0) call ERROR("Allocation error line%adamp")
! 								atom%lines(kc)%a(:) = 0d0
! 								mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%a(:))
! 							allocate(atom%lines(kc)%aeff(n_cells), stat=alloc_status)
! 							if (alloc_status > 0) call ERROR("Allocation error line%eff")
! 							atom%lines(kc)%aeff(:) = 0d0
! 							allocate(atom%lines(kc)%r(n_cells), stat=alloc_status)
! 							if (alloc_status > 0) call ERROR("Allocation error line%r")
! 							atom%lines(kc)%r(:) = 0d0
! 							allocate(atom%lines(kc)%r1(n_cells), stat=alloc_status)
! 							if (alloc_status > 0) call ERROR("Allocation error line%r1")
! 							atom%lines(kc)%r1(:) = 0d0
! 							
! 							mem_alloc_local = mem_alloc_local + 3*sizeof(atom%lines(kc)%aeff(n_cells))
! 
							endif
					else !Gaussian
						!local Gaussian profiles are also interpolated, no? if commented!
! 						if (associated(profile,local_profile_interp)) then
! 							allocate(atom%lines(kc)%phi(size(atom%lines(kc)%u), n_cells), stat=alloc_status)
! 							if (alloc_status > 0) call ERROR("Allocation error line%phi")
! 							atom%lines(kc)%phi(:,:) = 0d0
! 							size_phi_all = size_phi_all + sizeof(atom%lines(kc)%phi(:,:)) + sizeof(atom%lines(kc)%u(:))
! 							mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%phi(:,:)) + sizeof(atom%lines(kc)%u(:))
! 							!elseif
						if (associated(profile,local_profile_dk)) then
								!doesn't include the vel shift, just a local profile evaluated on the nlte grid.
							Nlam = atom%lines(kc)%Nlambda
							if (allocated(atom%lines(kc)%u)) deallocate(atom%lines(kc)%u)
							if (allocated(atom%lines(kc)%u)) call error("line%u should not be allocated here!")

							allocate(atom%lines(kc)%phi(Nlam, n_cells), stat=alloc_status)

							if (alloc_status > 0) call ERROR("Allocation error line%phi")
							atom%lines(kc)%phi(:,:) = 0d0
							size_phi_all = size_phi_all + sizeof(atom%lines(kc)%phi(:,:))
							mem_alloc_local = mem_alloc_local + sizeof(atom%lines(kc)%phi(:,:))
							
						!else local_profile_thomson or local_profile_voigt and gaussians are computed explictely
						endif
					endif

                  
				case ('ATOMIC_CONTINUUM')

          !Tex
          			if (atom%active) then
						allocate(atom%continua(kc)%Tex(n_cells), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%Tex")
						atom%continua(kc)%Tex(:) = T(:)
						mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kc)%Tex(:))

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
						mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kc)%Rij)*2

          !twohnu3_c2 for continuum waves (Not needed ?)
						allocate(atom%continua(kc)%twohnu3_c2(atom%continua(kc)%Nlambda),stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%twohnu3_c2")
						atom%continua(kc)%twohnu3_c2(:) = twohc / lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred)**3.
						mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kc)%twohnu3_c2(:))
          			endif
          			
          			
          !special for bound-free cross-sections which do not depend on cell
					if (atom%continua(kc)%Hydrogenic) then !Kramer's formula with quantum mechanical correction
						allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
						if (alloc_status > 0) call ERROR("Allocation error cont%alpha")

						atom%continua(kc)%alpha(:) = 0d0
        		!computes only where it is not zero
						atom%continua(kc)%alpha(:) = H_bf_Xsection(atom%continua(kc), lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))
						n_eff = atom%stage(atom%continua(kc)%j)*sqrt(atom%Rydberg/(atom%E(atom%continua(kc)%j)-atom%E(atom%continua(kc)%i)))
						mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kc)%alpha(:))
					else !interpolation of the read Cross-section
!       cont is not an alias but a copy of contiua(kc), so cont%alpha not deallocated
        				allocate(atom%continua(kc)%alpha(atom%continua(kc)%Nlambda), stat=alloc_status)
        				if (alloc_status > 0) call ERROR("Allocation error cont%alpha")
         
!          atom%continua(kc)%alpha(:) = linear_1D(size(atom%continua(kc)%lambda_file),atom%continua(kc)%lambda_file,&
!             atom%continua(kc)%alpha_file, atom%continua(kc)%Nlambda,NLTEspec%lambda(atom%continua(kc)%Nblue:atom%continua(kc)%Nred))
						call bezier3_interp(size(atom%continua(kc)%lambda_file), atom%continua(kc)%lambda_file, atom%continua(kc)%alpha_file, & !read values
     	atom%continua(kc)%Nlambda, lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), atom%continua(kc)%alpha) !interpolation grid
						mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kc)%alpha(:))
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
						mem_alloc_local = mem_alloc_local + sizeof(atom%continua(kc)%alpha_nu(:,:))

					endif
 
         !enddo
				case default
					call Error("Transition type unknown", atom%at(k)%trtype)
				end select
			enddo
			atom => NULL()
		enddo
		
		write(*,*) " Size of all local profiles for interp:", size_phi_all/1024./1024., ' MB'
		mem_alloc_tot = mem_alloc_tot + mem_alloc_local
		write(*,'("Total memory allocated in alloc_atom_quantities:"(1ES17.8E3)" MB")') mem_alloc_local / 1024./1024.
	

	RETURN
	end SUBROUTINE alloc_atom_quantities
  
    !to do: if verbose, store the warnings in an other file.
    !For large atoms and large grid, he number of warnings
    !can be overwhelming !
    !The warning are here for debugging and informations,
    !In general populations inversion are small and can be omitted in
    !the opacity.
	SUBROUTINE compute_atom_quantities(icell,verbose)
	!To DO, do not updated things that are constant in nlte_loop
		integer, intent(in) :: icell
		type(AtomType), pointer :: atom
		integer, optional :: verbose
		integer :: verbose_mode = 0
		integer :: n, k, kc, alloc_status, i, j, Nblue, Nred, la, icell_d
		real(kind=dp) :: vbroad, aL, cte, cte2, adamp, wi, wj, gij, neff, chi_ion
!     	real(kind=dp), dimension(:), allocatable :: test_phi

		if (present(verbose)) then
		
			if (verbose <= 3) verbose_mode = verbose
		
		endif

    
		do n=1,Natom
			atom => Atoms(n)%ptr_atom
			do k = 1, atom%Ntr   
     
				kc = atom%at(k)%ik 
        
				select case (atom%at(k)%trtype)
        
				case ('ATOMIC_LINE')
					vbroad = atom%vbroad(icell)
         
					i=atom%lines(kc)%i; j=atom%lines(kc)%j
					Nblue = atom%lines(kc)%Nblue; Nred = atom%lines(kc)%Nred
					wj = 1.0; wi = 1.0
					if (ldissolve) then
						if (atom%ID=="H") then! .or. atom%ID=="He") then
							wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1.0),hydrogen%n(1,icell))
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1.0),hydrogen%n(1,icell))
						endif
					endif

					if (atom%lines(kc)%voigt) then

						!call Damping(icell, atom, kc, adamp)
						adamp = line_damping(icell, atom%lines(kc))
						
!-> if profile is interpolated, adamp is not stored! the profile is
						if (allocated(atom%lines(kc)%a)) atom%lines(kc)%a(icell) = adamp

!-> if the profile is not interpolated, a is allocated and adamp is stored in it
						if (associated(profile,local_profile_interp)) then
							atom%lines(kc)%phi(:,icell) = Voigt(size(atom%lines(kc)%u), adamp, atom%lines(kc)%u(:)/vbroad) / (SQRTPI * vbroad)
						elseif (associated(profile,local_profile_dk)) then !evaluated on the nlte grid
							atom%lines(kc)%phi(:,icell) = Voigt(atom%lines(kc)%Nlambda, adamp, (lambda(Nblue:Nred)-atom%lines(kc)%lambda0)/atom%lines(kc)%lambda0 * clight/vbroad) / (SQRTPI * vbroad)								
						
!->these are not kept in memory anymore, only need to store eta(n_cells) but
!need also a function that call line_damping and recompute eta from aeff and ratio if necessary! 
! 						elseif (associated(Profile, local_profile_thomson)) then
! 							aL = line%a(icell) * vbroad !(m/s), adamp in doppler units
! 		
! 							aeff = (vbroad**5. + 2.69269*vbroad**4. * aL + 2.42843*vbroad**3. * aL**2. + &
! 									4.47163*vbroad**2. *aL**3. + 0.07842*vbroad*aL**4. + aL**5.)**(0.2)
! 							ratio = aL/aeff
! 							eta = 1.36603*ratio - 0.47719*(ratio*ratio) + 0.11116*(ratio*ratio*ratio)
! 			
						endif
					else !Gaussian
!-> gaussian profiles are interpolated in case of interpolation too.-> Not if commented!!! check profiles also
! 						if (associated(profile,local_profile_interp)) then
! 							atom%lines(kc)%phi(:,icell) = exp(-(atom%lines(kc)%u(:)/atom%vbroad(icell))**2) / (SQRTPI *vbroad)
						!elseif
						if (associated(profile,local_profile_dk)) then
							atom%lines(kc)%phi(:,icell) = exp(-( (lambda(Nblue:Nred)-atom%lines(kc)%lambda0)/atom%lines(kc)%lambda0 * clight/vbroad )**2 )/ (SQRTPI * vbroad)								
						endif
					endif !line voigt
					

!Include occupa, cannot be lower, since ni =  niwi and nj = njwj

					if (verbose_mode > 0) then
						if (verbose_mode == 2) then
							if ((wj/wi * atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij) <= 0 ) then
								call warning("Background line, population inversion!")
							endif
							write(*,"('at cell '(1I9), ' for atom '(1A2), ' line: '(1I3)' ->'(1I3) ) ") icell, atom%ID, i, j
							write(*,"(' lambda0='(1F12.5)' nm')") atom%lines(kc)%lambda0
							write(*,"('w(i)='(1ES20.7E3), '; w(j)='(1ES20.7E3))") wi, wj
							write(*,"('n(i)='(1ES20.7E3)' m^-3;',' n(j)='(1ES20.7E3)' m^-3;', ' gij='(1ES20.7E3))") atom%n(i,icell), atom%n(j,icell), atom%lines(kc)%gij
							write(*,"('nstar(i)='(1ES20.7E3)' m^-3;',' nstar(j)='(1ES20.7E3)' m^-3')") atom%nstar(i,icell), atom%nstar(j,icell)
!           						write(*,*) atom%n(i,icell), wi*atom%n(j, icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))!atom%continua(kc)%gij(la, icell)
          					write(*,"('(n(i)-n(j)gij)/n(i)='(1ES20.7E3)' %;',' (n(i)-n(j)gij)/n(j)gij='(1ES20.7E3)' %')") 100 * (wj/wi * atom%n(i,icell)-atom%n(j,icell)*atom%lines(kc)%gij)/atom%n(i,icell), 100 * (wj/wi * atom%n(i,icell)-atom%n(j,icell)*atom%lines(kc)%gij)/(atom%n(j,icell)*atom%lines(kc)%gij)
							write(*,"('n(i)/N='(1ES20.7E3)' %;', ' n(j)/N='(1ES20.7E3)' %')") 100 * atom%n(i,icell) / ntotal_atom(icell,atom), 100 * atom%n(j,icell) / ntotal_atom(icell,atom)
							write(*,"('nstar(i)/N='(1ES20.7E3)' %;', ' nstar(j)/N='(1ES20.7E3)' %')") 100 * atom%nstar(i,icell) / ntotal_atom(icell,atom), 100 * atom%nstar(j,icell) / ntotal_atom(icell,atom)
							write(*,"('n(i)/ne='(1ES20.7E3)' %;', ' n(i)/ne < prec ='(1L1))") 100 * atom%n(i,icell) / ne(icell), ( atom%n(i,icell) / ne(icell) < frac_ne_limit)
							write(*,"('n(j)/ne='(1ES20.7E3)' %;', ' n(j)/ne < prec ='(1L1))") 100 * atom%n(j,icell) / ne(icell), ( atom%n(j,icell) / ne(icell) < frac_ne_limit)
						elseif (verbose_mode == 1) then
							if ((wj/wi * atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij) <= 0 ) then
								call warning("Background line, population inversion!")
								write(*,"('at cell '(1I9), ' for atom '(1A2), ' line: '(1I3)' ->'(1I3) ) ") icell, atom%ID, i, j
								write(*,"(' lambda0='(1F12.5)' nm')") atom%lines(kc)%lambda0
								write(*,"('w(i)='(1ES20.7E3), '; w(j)='(1ES20.7E3))") wi, wj
								write(*,"('n(i)='(1ES20.7E3)' m^-3;',' n(j)='(1ES20.7E3)' m^-3;', ' gij='(1ES20.7E3))") atom%n(i,icell), atom%n(j,icell), atom%lines(kc)%gij
								write(*,"('nstar(i)='(1ES20.7E3)' m^-3;',' nstar(j)='(1ES20.7E3)' m^-3')") atom%nstar(i,icell), atom%nstar(j,icell)
!           						write(*,*) atom%n(i,icell), wi*atom%n(j, icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))!atom%continua(kc)%gij(la, icell)
          						write(*,"('(n(i)-n(j)gij)/n(i)='(1ES20.7E3)' %;',' (n(i)-n(j)gij)/n(j)gij='(1ES20.7E3)' %')") 100 * (wj/wi * atom%n(i,icell)-atom%n(j,icell)*atom%lines(kc)%gij)/atom%n(i,icell), 100 * (wj/wi * atom%n(i,icell)-atom%n(j,icell)*atom%lines(kc)%gij)/(atom%n(j,icell)*atom%lines(kc)%gij)
								write(*,"('n(i)/N='(1ES20.7E3)' %;', ' n(j)/N='(1ES20.7E3)' %')") 100 * atom%n(i,icell) / ntotal_atom(icell,atom), 100 * atom%n(j,icell) / ntotal_atom(icell,atom)
								write(*,"('nstar(i)/N='(1ES20.7E3)' %;', ' nstar(j)/N='(1ES20.7E3)' %')") 100 * atom%nstar(i,icell) / ntotal_atom(icell,atom), 100 * atom%nstar(j,icell) / ntotal_atom(icell,atom)
								write(*,"('n(i)/ne='(1ES20.7E3)' %;', ' n(i)/ne < prec ='(1L1))") 100 * atom%n(i,icell) / ne(icell), ( atom%n(i,icell) / ne(icell) < frac_ne_limit)
								write(*,"('n(j)/ne='(1ES20.7E3)' %;', ' n(j)/ne < prec ='(1L1))") 100 * atom%n(j,icell) / ne(icell), ( atom%n(j,icell) / ne(icell) < frac_ne_limit)
		
          					endif
          				elseif (verbose_mode == 3) then
							if ((wj/wi * atom%n(i,icell) - atom%n(j, icell)*atom%lines(kc)%gij) <= 0 ) then
								call warning("Background line, population inversion!")
								write(*,"(' => at cell '(1I9), ' for atom '(1A2), ' line: '(1I3)' ->'(1I3) ) ") icell, atom%ID, i, j
								write(*,"('  -- n(i)/N='(1ES20.9E3)' m^-3;',' n(j)/N='(1ES20.7E3)' m^-3;', ' gij='(1ES20.9E3))") atom%n(i,icell)/ntotal_atom(icell,atom), atom%n(j,icell)/ntotal_atom(icell,atom), atom%lines(kc)%gij
							endif
						endif
					endif
    
				case ("ATOMIC_CONTINUUM")
        
					i=atom%continua(kc)%i; j=atom%continua(kc)%j
					Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
					chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
					neff = atom%stage(j) * sqrt(atom%Rydberg / (atom%E(j) - atom%E(i)) )
					gij = atom%nstar(i,icell)/atom%nstar(j,icell)

					icell_d = 1
					wj = 1.0; wi = 1.0
					if (ldissolve) then
						if (atom%ID=="H") then! .or. atom%ID=="He") then
							wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(j)),hydrogen%n(1,icell))
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
!					if (atom%active) then
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
! 					endif
						
			
					if (verbose_mode > 0) then
						if (verbose_mode == 2) then
							!do la=1, atom%continua(kc)%Nlambda
							la = atom%continua(kc)%Nlambda
! 							if (atom%n(i,icell) - atom%n(j, icell)*atom%continua(kc)%gij(la, icell) <= 0 ) then
								if (atom%n(i,icell) - atom%n(j,icell) * gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)) <=0 ) then
									call Warning("Background continuum, population inversion !")
								endif
								write(*,"('at cell '(1I9), ' for atom '(1A2), ' cont: '(1I3)' ->'(1I3) ) ") icell, atom%ID, i, j
								write(*,"(' lambda0='(1F12.5)' nm;', ' lambda='(1F12.5)' nm')") atom%continua(kc)%lambda0, lambda_cont(Nblue+la-1)
								write(*,"('w(i)='(1ES20.7E3), '; w(j)='(1ES20.7E3))") wi, wj
								write(*,"('n(i)='(1ES20.7E3)' m^-3;',' n(j)='(1ES20.7E3)' m^-3;', ' gij='(1ES20.7E3))") atom%n(i,icell), atom%n(j,icell), gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))
								write(*,"('nstar(i)='(1ES20.7E3)' m^-3;',' nstar(j)='(1ES20.7E3)' m^-3')") atom%nstar(i,icell), atom%nstar(j,icell)
!           						write(*,*) atom%n(i,icell), wi*atom%n(j, icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))!atom%continua(kc)%gij(la, icell)
          						write(*,"('(n(i)-n(j)gij)/n(i)='(1ES20.7E3)' %')") 100 * (atom%n(i,icell)-wi*atom%n(j,icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)))/atom%n(i,icell)
								write(*,"('n(i)/N='(1ES20.7E3)' %;', ' n(j)/N='(1ES20.7E3)' %')") 100 * atom%n(i,icell) / ntotal_atom(icell,atom), 100 * atom%n(j,icell) / ntotal_atom(icell,atom)
								write(*,"('nstar(i)/N='(1ES20.7E3)' %;', ' nstar(j)/N='(1ES20.7E3)' %')") 100 * atom%nstar(i,icell) / ntotal_atom(icell,atom), 100 * atom%nstar(j,icell) / ntotal_atom(icell,atom)
								write(*,"('n(i)/ne='(1ES20.7E3)' %;', ' n(i)/ne < prec ='(1L1))") 100 * atom%n(i,icell) / ne(icell), ( atom%n(i,icell) / ne(icell) < frac_ne_limit)
								write(*,"('n(j)/ne='(1ES20.7E3)' %;', ' n(j)/ne < prec ='(1L1))") 100 * atom%n(j,icell) / ne(icell), ( atom%n(j,icell) / ne(icell) < frac_ne_limit)
								!exit
							!enddo
						elseif( verbose_mode==1 ) then
							do la=1, atom%continua(kc)%Nlambda
! 							if (atom%n(i,icell) - atom%n(j, icell)*atom%continua(kc)%gij(la, icell) <= 0 ) then
								if (atom%n(i,icell) - atom%n(j,icell) * gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)) <=0 ) then
									call Warning("Background continuum, population inversion !")

									write(*,"('at cell '(1I9), ' for atom '(1A2), ' cont: '(1I3)' ->'(1I3) ) ") icell, atom%ID, i, j
									write(*,"(' lambda0='(1F12.5)' nm;', ' lambda='(1F12.5)' nm')") atom%continua(kc)%lambda0, lambda_cont(Nblue+la-1)
									write(*,"('w(i)='(1ES20.7E3), '; w(j)='(1ES20.7E3))") wi, wj
									write(*,"('n(i)='(1ES20.7E3)' m^-3;',' n(j)='(1ES20.7E3)' m^-3;', ' gij='(1ES20.7E3))") atom%n(i,icell), atom%n(j,icell), gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))
									write(*,"('nstar(i)='(1ES20.7E3)' m^-3;',' nstar(j)='(1ES20.7E3)' m^-3')") atom%nstar(i,icell), atom%nstar(j,icell)
!           						write(*,*) atom%n(i,icell), wi*atom%n(j, icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))!atom%continua(kc)%gij(la, icell)
          							write(*,"('(n(i)-n(j)gij)/n(i)='(1ES20.7E3)' %')") 100 * (atom%n(i,icell)-wi*atom%n(j,icell)*gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)))/atom%n(i,icell)
									write(*,"('n(i)/N='(1ES20.7E3)' %;', ' n(j)/N='(1ES20.7E3)' %')") 100 * atom%n(i,icell) / ntotal_atom(icell,atom), 100 * atom%n(j,icell) / ntotal_atom(icell,atom)
									write(*,"('nstar(i)/N='(1ES20.7E3)' %;', ' nstar(j)/N='(1ES20.7E3)' %')") 100 * atom%nstar(i,icell) / ntotal_atom(icell,atom), 100 * atom%nstar(j,icell) / ntotal_atom(icell,atom)
									write(*,"('n(i)/ne='(1ES20.7E3)' %;', ' n(i)/ne < prec ='(1L1))") 100 * atom%n(i,icell) / ne(icell), ( atom%n(i,icell) / ne(icell) < frac_ne_limit)
									write(*,"('n(j)/ne='(1ES20.7E3)' %;', ' n(j)/ne < prec ='(1L1))") 100 * atom%n(j,icell) / ne(icell), ( atom%n(j,icell) / ne(icell) < frac_ne_limit)
									exit
								endif
							enddo
          				elseif (verbose_mode == 3) then
							do la=1, atom%continua(kc)%Nlambda
! 							if (atom%n(i,icell) - atom%n(j, icell)*atom%continua(kc)%gij(la, icell) <= 0 ) then
								if (atom%n(i,icell) - atom%n(j,icell) * gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell)) <=0 ) then
									call Warning("Background continuum, population inversion !")
									write(*,"(' => at cell '(1I9), ' for atom '(1A2), ' cont: '(1I3)' ->'(1I3) ) ") icell, atom%ID, i, j
									write(*,"('  -- n(i)/N='(1ES20.9E3)' m^-3;',' n(j)/N='(1ES20.7E3)' m^-3;', ' gij='(1ES20.9E3))") atom%n(i,icell)/ntotal_atom(icell,atom), atom%n(j,icell)/ntotal_atom(icell,atom), gij * exp(-hc_k/lambda_cont(Nblue+la-1)/T(icell))
									exit
								endif
							enddo
						endif
					endif

				case default
					call Error("Transition type unknown", atom%at(k)%trtype)
				end select
			enddo
			atom => NULL()
		enddo  
  
  
	RETURN
	end SUBROUTINE compute_atom_quantities
	
	!metal_bf typically
! 	subroutine opacity_atom_bf_lte(icell)
! 
! 		integer, intent(in)											:: icell
! 		integer														:: m, kr, kc, i, j, Nblue, Nred, la
! 		type (AtomType), pointer									:: atom
! 		real(kind=dp)												:: n_eff, wj, wi
! 		real(kind=dp)												:: Diss, chi_ion, gij
! 
! 
! 		do m=1,Npassiveatoms
! 
! 			atom => PassiveAtoms(m)%ptr_atom
! 
! 			do kc=atom%Ntr_line+1,atom%Ntr
! 				kr = atom%at(kc)%ik 
! 
! 				i = atom%continua(kr)%i
! 				j = atom%continua(kr)%j 
! 
! 				Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred
! 				wj = 1.0; wi = 1.0
! 
! 				if (ldissolve) then
! 					if (atom%ID=="H") then
! 						n_eff = real(i,kind=dp)
! 						wi = wocc_n(icell, n_eff, real(atom%stage(i)), real(atom%stage(j)),hydrogen%n(1,icell))
! 					else
! 						n_eff = atom%stage(j)*sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))
! 					endif
! 				endif
! 				
! 				chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
! 
! 				do la=1, atom%continua(kr)%nlambda
! 				
! 					Diss = D_i(icell, n_eff, real(atom%stage(i)), 1.0, lambda_cont(Nblue+la-1), atom%continua(kc)%lambda0, chi_ion)
! 					gij = atom%nstar(i,icell)/atom%nstar(j,icell) * exp(-hc_k/T(icell)/lambda_cont(Nblue+la-1))
! 
! 					
! 					if ( (atom%n(i,icell) - atom%n(j,icell) * gij) > 0.0) then
! 						
! 							chi_c(Nblue+la-1,icell) = chi_c(Nblue+la-1,icell) + &
! 							diss * atom%continua(kr)%alpha(la) * (atom%n(i,icell) - gij*atom%n(j,icell)) 
!        				
! 							eta_c(Nblue+la-1,icell) = eta_c(Nblue+la-1,icell) + &
! 							diss * atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * gij * atom%n(j,icell)
! 					else
! 					
! 							eta_c(Nblue+la-1,icell) = eta_c(Nblue+la-1,icell) + &
! 							diss * atom%continua(kr)%alpha(la) * atom%continua(kr)%twohnu3_c2(la) * gij * atom%n(j,icell)   
! 							
! 							chi_c(Nblue+la-1,icell) = chi_c(Nblue+la-1,icell) - &
! 							diss * atom%continua(kr)%alpha(la) * (atom%n(i,icell) - gij*atom%n(j,icell)) 
! 
! 					endif
! 				enddo				
! 
! 			end do ! loop over Ncont
! 			atom => NULL()
! 		end do !loop over metals
! 
! 	return
! 	end subroutine opacity_atom_bf_lte
	
	subroutine background_continua (icell)
		integer, intent(in) :: icell
		integer :: la, nat
		real(kind=dp), dimension(Nlambda_cont) :: chi, eta, Bp

		Bp = Bpnu(T(icell), lambda_cont)

		chi_c(:,icell) = Thomson(icell)!0.0_dp
		eta_c(:,icell) = 0.0_dp
! 		sca_c(:,icell) = 0.0_dp

! 		call Hminus_bf_wishart(icell, Nlambda_cont, lambda_cont, chi, eta) !->better with turbospec
		call Hminus_bf_geltman(icell, Nlambda_cont, lambda_cont, chi, eta) !->better with rh
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + eta(:)

		call Hminus_ff_bell_berr(icell, Nlambda_cont, lambda_cont, chi)
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + chi(:) * Bp(:)
!  
		call Hydrogen_ff(icell, Nlambda_cont, lambda_cont, chi)
		chi_c(:,icell) = chi_c(:,icell) + chi(:)
		eta_c(:,icell) = eta_c(:,icell) + chi(:) * Bp(:)
 
! 		!now atomic LTE bound-free
		if (NpassiveAtoms > 0) then
			call lte_bound_free(icell, Nlambda_cont, lambda_cont, chi, eta)
			chi_c(:,icell) = chi_c(:,icell) + chi(:)
			eta_c(:,icell) = eta_c(:,icell) + eta(:)
		endif

! 		sca_c(:,icell) = Thomson(icell)
! 		!Rayleigh scattering is included no more at the moment. Might be necessary for some temperature
! 		!below the Lyman limit
! 
! 		!Total opac once source functions are known
! 		chi_c(:,icell) = chi_c(:,icell) +  sca_c(:,icell)

	return
	end subroutine background_continua
	
	subroutine background_continua_lambda (icell, Nx, x, chiout, etaout)
	!Scattering coefficients at the moment recomputed in propagation since it is prop to Jnu * thomson only (no H scatt at the moment)
		integer, intent(in) :: icell, Nx
		integer :: la, nat
		real(kind=dp), dimension(Nx), intent(in) :: x
		real(kind=dp), dimension(Nx), intent(out) :: chiout, etaout
		real(kind=dp), dimension(Nx) :: chi, eta, Bp

		Bp = Bpnu(T(icell),x)

		chiout(:) = thomson(icell)
		etaout(:) = 0.0_dp

! 		!call Hminus_bf_wishart(icell, N, lambda, chi, eta) !->better with turbospec
		call Hminus_bf_geltman(icell,Nx, x, chi, eta) !->better with rh
		chiout(:) = chiout(:) + chi(:)
		etaout(:) = etaout(:) + eta(:)

		call Hminus_ff_bell_berr(icell, Nx, x, chi)
		chiout(:) = chiout(:) + chi(:)
		etaout(:) = etaout(:) + chi(:) * Bp(:)

		call Hydrogen_ff(icell, Nx, x, chi)
		chiout(:) = chiout(:) + chi(:)
		etaout(:) = etaout(:) + chi(:) * Bp(:)
 
		!now atomic LTE bound-free
		if (NpassiveAtoms > 0) then
			call lte_bound_free(icell, Nx, x, chi, eta)
			chiout(:) = chiout(:) + chi(:)
			etaout(:) = etaout(:) + eta(:)
		endif

	return
	end subroutine background_continua_lambda
	
	subroutine compute_background_continua()
	!$ use omp_lib
		integer :: icell, id, icell0
		icell0 = 1
		!$omp parallel &
		!$omp default(none) &
		!$omp private(icell,id,icell0) &
		!$omp shared(llimit_mem,Nlambda_cont, lambda_cont, chi_c, eta_c, Nlambda, lambda) & 
		!$omp shared(n_cells, icompute_atomRT, NactiveAtoms, chi0_bb, eta0_bb)
		!$omp do schedule(dynamic,1)
		do icell=1, n_cells
			!$ id = omp_get_thread_num() + 1
			if (icompute_atomRT(icell) > 0) then
				call compute_atom_quantities(icell,verbose=0)!1,2,3
				!!need for BackgroundContinua
    			!!and lines 
				!!!!call background_continua(icell)	
				call background_continua_lambda(icell, Nlambda_cont, lambda_cont, chi_c(:,icell), eta_c(:,icell))
				!!call background_continua_lambda(icell, Nlambda, lambda, chi0_bb(:,icell), eta0_bb(:,icell))
				if (NactiveAtoms > 0) then
					call NLTE_bound_free(icell)
				endif
				
				if (.not.llimit_mem) then
					call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
				endif
		
			endif
		end do
		!$omp end do
		!$omp end parallel

	return
	end subroutine compute_background_continua
	
	!not used anymore ?
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

  
! 		call Hminus_bf_wishart(icell, 1, l, chi, eta)
		call Hminus_bf_geltman(icell,1,l,chi,eta) !->better with rh
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
		real(kind=dp) :: wi, wj, chi_ion, Diss, nn, gij, ni_on_nj_star
		type (element), pointer :: elem
  
		chi_c_nlte(:,icell) = 0d0
		eta_c_nlte(:,icell) = 0d0

		atom_loop : do nact = 1, Nactiveatoms
			aatom => ActiveAtoms(nact)%ptr_atom
			elem => Elements(aatom%periodic_table)%ptr_elem
   
			tr_loop : do kr = aatom%Ntr_line+1,aatom%Ntr


        		kc = aatom%at(kr)%ik !relative index of a transition among continua or lines


				Nred = aatom%continua(kc)%Nred; Nblue = aatom%continua(kc)%Nblue    	
				i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
				nn = n_eff(aatom%Rydberg, aatom%E(j), aatom%E(i), aatom%stage(j))!real(atom%stage(j)) * sqrt(atom%Rydberg/(atom%E(j)-atom%E(i)))

				wj = 1.0; wi = 1.0
				if (ldissolve) then
					if (aatom%ID=="H") then
												!nn
						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(j)),hydrogen%n(1,icell))
					endif
! 				else
! 					if (lambda_cont(Nred) > aatom%continua(kc)%lambda0) then
! 						write(*,*)lambda_cont(Nred), lambda_cont(aatom%continua(kc)%Nlambda+Nblue-1),aatom%continua(kc)%lambda0
! 						write(*,*) i, j, Nred, locate(lambda_cont, aatom%continua(kc)%lambda0), locate(lambda_cont, lambda_cont(Nred))
! 						call error("lambda(Nred) can't be larger than lambda0!")
! 					endif
				endif

				!get the ionisation potential for the ion ot be use in the dissolve fraction
				chi_ion = elem%ionpot(aatom%stage(j))
! 				ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
				ni_on_nj_star = aatom%nstar(i,icell)/aatom%nstar(j,icell)
				
				gij = ni_on_nj_star * exp(-hc_k/T(icell)/aatom%continua(kc)%lambda0)

				if ((aatom%n(i,icell) - aatom%n(j,icell) * gij) <= 0.0_dp) then
					cycle tr_loop
				endif

				freq_loop : do la=1,aatom%continua(kc)%Nlambda

					!beware even without diss it appears that sometime lambda(Nred) = lambda0 + epsilon > lambda0
					Diss = D_i(icell, nn, real(aatom%stage(i)), 1.0, lambda_cont(Nblue+la-1), aatom%continua(kc)%lambda0, chi_ion)
					gij = ni_on_nj_star * exp(-hc_k/T(icell)/lambda_cont(Nblue+la-1))

! 					if ((aatom%n(i,icell) - aatom%n(j,icell) * gij) > 0.0_dp) then

						chi_c_nlte(Nblue+la-1,icell) = chi_c_nlte(Nblue+la-1,icell) + &
            			Diss * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - gij * aatom%n(j,icell))

						eta_c_nlte(Nblue+la-1,icell) = eta_c_nlte(Nblue+la-1,icell) + &
            			Diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * gij * aatom%n(j,icell)

            
! 					else !neg or null
! 						eta_c_nlte(Nblue+la-1,icell) = eta_c_nlte(Nblue+la-1,icell) + &
!             			Diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * gij * aatom%n(j,icell)
					!assuming small inversions
						!chi_c_nlte(Nblue+la-1,icell) = chi_c_nlte(Nblue+la-1,icell) - Diss * aatom%continua(kc)%alpha(la) * (aatom%n(i,icell) - gij * aatom%n(j,icell))
!  
! 					endif

				enddo freq_loop

			end do tr_loop
   
   			elem => null()
			aatom => NULL()
		end do atom_loop
	
 
	return
	end subroutine NLTE_bound_free
 
! 	subroutine compute_nlte_bound_free
! 	!$ use omp_lib
! 		integer :: icell, id
! 
! 
! 		!$omp parallel &
! 		!$omp default(none) &
! 		!$omp private(icell,id) &
! 		!$omp shared(n_cells,icompute_atomRT)
! 		!$omp do schedule(dynamic,1)
! 		do icell=1,n_cells
! 			!$ id = omp_get_thread_num() + 1
! 			if (icompute_atomRT(icell) > 0) then
! 				call NLTE_bound_free(icell)
! 			endif
! 		end do
! 		!$omp end do
! 		!$omp end parallel
!  
! 	return
! 	end subroutine compute_nlte_bound_free
	
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

		!!linear interpolation of both arrays
     	do i=i0,Nlambda_cont-1
        	do j=j0,Nlambda
           		if (lambda(j)>=lambda_cont(i) .and. lambda(j)<=lambda_cont(i+1)) then
            		u = (lambda(j) - lambda_cont(i)) / (lambda_cont(i+1) - lambda_cont(i))
              		chii(j) = (1.0_dp - u) * chi0(i) + u * chi0(i+1)
              		etai(j) = (1.0_dp - u) * eta0(i) + u * eta0(i+1)
           		endif
        	enddo
     	enddo
 

	return
	end subroutine interp_background_opacity
	
	subroutine interp_continuum_local(icell, chii, etai)
		integer, intent(in) :: icell
		integer :: i, j, i0,j0
		!it seems I'm having an interpolation error at the first point, so just avoinding it..
		real(kind=dp), dimension(Nlambda_cont) :: chi0, eta0
		real(kind=dp), dimension(Nlambda), intent(out) :: chii, etai
		real(kind=dp) :: u
		
		chii = 0.0_dp
		etai = 0.0_dp
     
     	chi0 = chi_c(:,icell)
     	eta0 = eta_c(:,icell)

     	
     	if (NactiveAtoms>0) then
     		chi0 = chi0 + chi_c_nlte(:,icell)
     		eta0 = eta0 + eta_c_nlte(:,icell)
     	endif
     	

     	j0=Nlambda+1
		do j=1, Nlambda
        	if (lambda(j) > lambda_cont(1)) then!>=
           		j0 = j
           		exit
        	endif
     	enddo
     	

     	i0 = 2
    	do j=j0,Nlambda
        loop_i : do i=i0, Nlambda_cont
        		if (lambda_cont(i) > lambda(j)) then
            		u = (lambda(j) - lambda_cont(i-1)) / (lambda_cont(i) - lambda_cont(i-1))
              		chii(j) = (1.0_dp - u) * chi0(i-1)  + u * chi0(i)
                  	etai(j) = (1.0_dp - u) * eta0(i-1)  + u * eta0(i)
              		i0 = i
              		exit loop_i
           		endif
        	enddo loop_i
    	 enddo
    	 
    	 
    	 !in case the last and first wavelength are the same
    	 if (lambda_cont(Nlambda_cont) == lambda(Nlambda)) then
    	 	chii(Nlambda) = chi0(Nlambda_cont)
    	 	etai(Nlambda) = eta0(Nlambda_cont)
    	 endif
    	 
    	 if (lambda_cont(1) == lambda(1)) then
    	 	chii(1) = chi0(1)
    	 	etai(1) = eta0(1)
    	 endif
 
!  	write(*,*) j0
!  	write(*,*) chii(Nlambda-3:Nlambda)
!  	write(*,*) etai(Nlambda-3:Nlambda)
!  	write(*,*) chi0(Nlambda_cont-3:Nlambda_cont)
!  	write(*,*) eta0(Nlambda_cont-3:Nlambda_cont)
!  	write(*,*) lambda(Nlambda-3:Nlambda)
!  	write(*,*) lambda_cont(Nlambda_cont-3:Nlambda_cont)
!  	stop

	return
	end subroutine interp_continuum_local

	subroutine opacity_atom_loc(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)

		integer, intent(in) :: id, icell, iray
		logical, intent(in) :: iterate
		real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
		integer :: nact, Nred, Nblue, kc, kr, i, j
		type(AtomType), pointer :: aatom
		integer :: dk0, dk1, Nlam
		real(kind=dp) :: wi, wj, adamp, dv_dl
		!if there is 40 points for the unperturbed profile for all lines, the size is 40. If there are velocity fields
		!we add two times the shift. This is however, much lower than Nlambda !!
		real(kind=dp), dimension(Nlambda_max_line+2*dk_max) :: phi0
! 		real(kind=dp), dimension(Nlambda) :: phi0

		!-> chi and eta contains already chi_c, chi_c_nlte and eta_c, eta_c_nlte

! 		dv_dl = gradv(icell, x,y,z,x1,y1,z1,u,v,w,l,dk0)
! 		if (dk0 >= 0) then
! 			dk0 = min(dk_max, 6 * dk0)
! 		else
! 			dk0 = max(dk_min, 6 * dk0)
! 		endif
! 		dk0 = dk_max
  
		atom_loop : do nact = 1, Natom!Nactiveatoms
			aatom => Atoms(nact)%ptr_atom!ActiveAtoms(nact)%ptr_atom

			tr_loop : do kr = 1,aatom%Ntr_line

				kc = aatom%at(kr)%ik !relative index of a transition among continua or lines

				Nred = aatom%lines(kc)%Nred; Nblue = aatom%lines(kc)%Nblue
				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j

							
				Nblue = Nblue + dk_min
				Nred = Nred + dk_max
				Nlam = Nred - Nblue + 1
				
 				wj = 1.0; wi = 1.0
				if (ldissolve) then
					if (aatom%ID=="H") then
												!nn
						wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1),hydrogen%n(1,icell)) !1 for H
						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1),hydrogen%n(1,icell))
					endif
				endif 

				if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0_dp) then
					cycle tr_loop
				endif

				phi0(1:Nlam) = profile(aatom%lines(kc),icell,iterate,Nlam,lambda(Nblue:Nred), x,y,z,x1,y1,z1,u,v,w,l)
	

! 				if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then


					chi(Nblue:Nred,id) = chi(Nblue:Nred,id) + &
						hc_fourPI * aatom%lines(kc)%Bij * phi0(1:Nlam) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))

					eta(Nblue:Nred,id)= eta(Nblue:Nred,id) + &
						hc_fourPI * aatom%lines(kc)%Aji * phi0(1:Nlam) * aatom%n(j,icell)

! 				else !neg or null
! 					eta(Nblue:Nred,id)= eta(Nblue:Nred,id) + &
! 						hc_fourPI * aatom%lines(kc)%Aji * aatom%n(j,icell) * phi0(1:Nlam)
! 				endif
				
				if (iterate) then  !only for active atoms ? if not hogerheijde only
					!if (aatom%active) &
					aatom%lines(kc)%phi_loc(:,iray,id) = phi0(1:Nlam)
				endif

    
			end do tr_loop

			aatom => NULL()

		end do atom_loop


	return
	end subroutine opacity_atom_loc
	
	!no level dissolution yet
	!fills also eta_atoms
	!check conditions on negative opac
	subroutine cross_coupling_terms(id, icell, iray)
		integer, intent(in) :: id, icell, iray
		integer :: nact, j, i, kr, kc, Nb, Nr, la, Nl, icell_d
		type (AtomType), pointer :: aatom
		real(kind=dp) :: gij, wi, wj, chicc, wl,  wphi, ni_on_nj_star
		
		!for one ray
		Uji_down(:,:,:,id) = 0.0_dp
		chi_down(:,:,:,id) = 0.0_dp
		chi_up(:,:,:,id)   = 0.0_dp
		
		eta_atoms(:,:,id) = 0.0_dp
	
		aatom_loop : do nact=1, Nactiveatoms
			aatom => ActiveAtoms(nact)%ptr_atom
			
				cont_loop : do kr = aatom%Ntr_line+1, aatom%Ntr
			
					kc = aatom%at(kr)%ik

					j = aatom%continua(kc)%j
					i = aatom%continua(kc)%i
					Nb = aatom%continua(kc)%Nb; Nr = aatom%continua(kc)%Nr
					Nl = Nr - Nb + 1

! 					ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
					ni_on_nj_star = aatom%nstar(i,icell)/aatom%nstar(j,icell)

				
					icell_d = 1
					if (ldissolve) then
						if (aatom%ID=="H") icell_d = icell
					endif
					
					gij = ni_on_nj_star * exp(-hc_k/T(icell)/aatom%continua(kc)%lambda0)			

					if (aatom%n(i,icell) - gij*aatom%n(j,icell) <= 0.0_dp) then
						cycle cont_loop
					endif			

		
					freq_loop : do la=1, Nl
						if (la==1) then
							wl = 0.5*(lambda(Nb+1)-lambda(Nb)) / lambda(Nb)
						elseif (la==Nl) then
							wl = 0.5*(lambda(Nb+la-1)-lambda(Nb+la-2)) / lambda(Nb+la-1)
						else
							wl = 0.5*(lambda(Nb+la)-lambda(Nb+la-2)) / lambda(Nb+la-1)
						endif

						gij = ni_on_nj_star * exp(-hc_k/T(icell)/lambda(Nb+la-1))			
										
 						!small inversions
						!chicc = wl * fourpi_h * aatom%continua(kc)%alpha_nu(la,icell_d) * abs(aatom%n(i,icell) - gij*aatom%n(j,icell))
						chicc = wl * fourpi_h * aatom%continua(kc)%alpha_nu(la,icell_d) * (aatom%n(i,icell) - gij*aatom%n(j,icell))
! 						if (chicc < 0.0) chicc = 0.0_dp !should not happend

						Uji_down(Nb+la-1,j,nact,id) = Uji_down(Nb+la-1,j,nact,id) + &
							aatom%continua(kc)%alpha_nu(la,icell_d) * (twohc/lambda(Nb+la-1)**3) * gij
							
						chi_down(Nb+la-1,j,nact,id) = chi_down(Nb+la-1,j,nact,id) + chicc
							
						chi_up(Nb+la-1,i,nact,id) = chi_up(Nb+la-1,i,nact,id) + chicc
!check consistency with nlte b-f and how negative opac is handled 
						!if (chicc > 0.0_dp) &
						eta_atoms(Nb+la-1,nact,id) = eta_atoms(Nb+la-1,nact,id) + &
							aatom%continua(kc)%alpha_nu(la,icell_d) * (twohc/lambda(Nb+la-1)**3)  * gij * aatom%n(j,icell)
					enddo freq_loop
			
				enddo cont_loop		
			
			!for each line eventually 
			wi = 1.0; wj = 1.0
			if (ldissolve) then
				if (aatom%ID=="H") then
												!nn
					wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1),hydrogen%n(1,icell))!1 for H
					wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1),hydrogen%n(1,icell))
				endif
			endif 

			line_loop : do kr = 1, aatom%Ntr_line
			
				kc = aatom%at(kr)%ik
				

				j = aatom%lines(kc)%j
				i = aatom%lines(kc)%i
				
				if (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell) <= 0.0_dp) then
					cycle line_loop
				endif				
				
				Nb = aatom%lines(kc)%Nblue; Nr = aatom%lines(kc)%Nred
			
				Nl = Nr-dk_min+dk_max-Nb+1
				wphi = 0.0
				freq2_loop : do la=1, Nl
					if (la==1) then
						wl = 0.5*1d3*hv
						!wl = 0.5*(lambda(Nb+1)-lambda(Nb)) * clight / aatom%lines(kc)%lambda0
					elseif (la==Nl) then
						wl = 0.5*1d3*hv
						!wl = 0.5*(lambda(Nr)-lambda(Nr-1)) * clight / aatom%lines(kc)%lambda0
					else
						wl = 1d3*hv
						!wl = 0.5*(lambda(Nb+la+1)-lambda(Nb+la-1)) * clight / aatom%lines(kc)%lambda0
					endif
						
					
					Uji_down(Nb+dk_min-1+la,j,nact,id) = Uji_down(Nb+dk_min-1+la,j,nact,id) + &
						hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi_loc(la,iray,id)
					
					!small inversions
! 					if (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell) >= 0.0_dp) then
					
					chi_down(Nb+dk_min-1+la,j,nact,id) = chi_down(Nb+dk_min-1+la,j,nact,id) + &
						wl * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(la,iray,id) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
					
					chi_up(Nb+dk_min-1+la,i,nact,id) = chi_up(Nb+dk_min-1+la,i,nact,id) + &
						wl * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(la,iray,id) * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))

						
! 					endif
										
					eta_atoms(Nb+dk_min-1+la,nact,id) = eta_atoms(Nb+dk_min-1+la,nact,id) + &
						hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi_loc(la,iray,id) * aatom%n(j,icell)

					
					wphi = wphi + wl * aatom%lines(kc)%phi_loc(la,iray,id)
				enddo freq2_loop

				chi_down(Nb+dk_min:Nr+dk_max,j,nact,id) = chi_down(Nb+dk_min:Nr+dk_max,j,nact,id) / wphi
				chi_up(Nb+dk_min:Nr+dk_max,i,nact,id) = chi_up(Nb+dk_min:Nr+dk_max,i,nact,id) / wphi

			enddo line_loop

					
		enddo aatom_loop
	
	return
	end subroutine cross_coupling_terms

	subroutine opacity_atom_zeeman_loc(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l, iterate)

		integer, intent(in) :: id, icell, iray
		logical, intent(in) :: iterate
		real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
		integer :: nact, Nred, Nblue, kc, kr, i, j, m
		type(AtomType), pointer :: aatom
		integer :: dk0, dk1, Nlam
		real(kind=dp) :: wi, wj, adamp, dv_dl, chil, etal
		real(kind=dp), dimension(Nlambda_max_line+2*dk_max) :: phi0
		real(kind=dp), dimension(Nlambda_max_line+2*dk_max,3) :: phiZ, psiZ

 		chiQUV_p(:,:,id) = 0.0_dp
 		etaQUV_p(:,:,id) = 0.0_dp
 		rho_p(:,:,id) = 0.0_dp
  
		atom_loop : do nact = 1, Natom
			aatom => Atoms(nact)%ptr_atom

			tr_loop : do kr = 1,aatom%Ntr_line

				kc = aatom%at(kr)%ik

				Nred = aatom%lines(kc)%Nred; Nblue = aatom%lines(kc)%Nblue
				i = aatom%lines(kc)%i;j = aatom%lines(kc)%j

							
				Nblue = Nblue + dk_min
				Nred = Nred + dk_max
				Nlam = Nred - Nblue + 1
				
 				wj = 1.0; wi = 1.0
				if (ldissolve) then
					if (aatom%ID=="H") then
												!nn
						wj = wocc_n(icell, real(j,kind=dp), real(aatom%stage(j)), real(aatom%stage(j)+1),hydrogen%n(1,icell)) !1 for H
						wi = wocc_n(icell, real(i,kind=dp), real(aatom%stage(i)), real(aatom%stage(i)+1),hydrogen%n(1,icell))
					endif
				endif 

				if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) <= 0.0_dp) then
					cycle tr_loop
				endif
				
				!fixed at the moment
				call local_profile_zv(aatom%lines(kc),icell,iterate,Nlam,lambda(Nblue:Nred),phi0,phiZ, psiZ, x,y,z,x1,y1,z1,u,v,w,l)

				etal = hc_fourPI * aatom%lines(kc)%Aji * aatom%n(j,icell)
				chil = hc_fourPI * aatom%lines(kc)%Bij * (aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))

! 				if ((aatom%n(i,icell)*wj/wi - aatom%n(j,icell)*aatom%lines(kc)%gij) > 0.0_dp) then

					
					chi(Nblue:Nred,id) = chi(Nblue:Nred,id) + chil * phi0(1:Nlam)
					eta(Nblue:Nred,id)= eta(Nblue:Nred,id) + etal * phi0(1:Nlam) 
						
					do m=1,3
						chiQUV_p(Nblue:Nred,m,id) = chiQUV_p(Nblue:Nred,m,id) + chil * phiz(1:Nlam,m)
						etaQUV_p(Nblue:Nred,m,id) = etaQUV_p(Nblue:Nred,m,id) + etal * phiz(1:Nlam,m)
						rho_p(Nblue:Nred,m,id) = rho_p(Nblue:Nred,m,id) + chil * psiz(1:Nlam,m)
					enddo

! 				else !neg or null
! 					eta(Nblue:Nred,id)= eta(Nblue:Nred,id) + etal * phi0(1:Nlam) 
! 					do m=1,3
! ! 						chiQUV_p(Nblue:Nred,m,id) = chiQUV_p(Nblue:Nred,m,id) + chil * phiz(1:Nlam,m)
! 						etaQUV_p(Nblue:Nred,m,id) = etaQUV_p(Nblue:Nred,m,id) + etal * phiz(1:Nlam,m)
! ! 						rho_p(Nblue:Nred,m,id) = rho_p(Nblue:Nred,m,id) + chil * psiz(1:Nlam,m)
! 					enddo						
! 				endif

    
			end do tr_loop

			aatom => NULL()

		end do atom_loop


	return
	end subroutine opacity_atom_zeeman_loc
	
end module Opacity

!  	subroutine cross_coupling_terms(id, icell, iray)
! 		integer, intent(in) :: id, iray, icell
! 		integer :: nact, j, i, kr, kc, Nb, Nr, la
! 		type (AtomType), pointer :: aatom
! 		real(kind=dp) :: gij, wi, wj, diss
! 		
! 		if (iray==1) then
! 			!initialize for first ray
! 			Uji_down(:,:,:,:,id) = 0.0_dp
! 			chi_down(:,:,:,:,id) = 0.0_dp
! 			chi_up(:,:,:,:,id)   = 0.0_dp
! 		endif
! 
! 		wi = 1.0; wj = 1.0; diss = 1.0
! 	
! 		aatom_loop : do nact=1, Nactiveatoms
! 			aatom => ActiveAtoms(nact)%ptr_atom
! 			
! 			
! 			if (iray==1) then !initialize
! 
! 				do kr = aatom%Ntr_line+1, aatom%Ntr
! 			
! 					kc = aatom%at(kr)%ik
! 
! 					j = aatom%continua(kc)%j
! 					i = aatom%continua(kc)%i
! 					Nb = aatom%continua(kc)%Nblue; Nr = aatom%continua(kc)%Nred
! 					
! 					do la=1, aatom%continua(kc)%Nlambda
! 						gij = aatom%nstar(i,icell)/aatom%nstar(j,icell) * exp(-hc_k/T(icell)/lambda_cont(Nb+la-1))							
! 	
! 						
! 						Uji_down(Nb+la-1,:,j,nact,id) = Uji_down(Nb+la-1,:,j,nact,id) + &
! 							diss * aatom%continua(kc)%alpha(la) * aatom%continua(kc)%twohnu3_c2(la) * gij
! 						chi_down(Nb+la-1,:,j,nact,id) = chi_down(Nb+la-1,:,j,nact,id) + &
! 							aatom%continua(kc)%alpha(la) * abs(aatom%n(i,icell) - gij*aatom%n(j,icell)) 
! 						chi_up(Nb+la-1,:,i,nact,id) = chi_up(Nb+la-1,:,i,nact,id) + &
! 							aatom%continua(kc)%alpha(la) * abs(aatom%n(i,icell) - gij*aatom%n(j,icell))
! 						
! 						
! 					enddo
! 			
! 				enddo			
! 			
! 			endif !over iray
! 
! 			do kr = 1, aatom%Ntr_line
! 			
! 				kc = aatom%at(kr)%ik
! 
! 				j = aatom%lines(kc)%j
! 				i = aatom%lines(kc)%i
! 				Nb = aatom%lines(kc)%Nblue; Nr = aatom%lines(kc)%Nred
! 						
! 				Uji_down(Nb+dk_min:Nr+dk_max,iray,j,nact,id) = Uji_down(Nb+dk_min:Nr+dk_max,iray,j,nact,id) + &
! 					hc_fourPI * aatom%lines(kc)%Aji * aatom%lines(kc)%phi_loc(:,iray,id)
! 				chi_down(Nb+dk_min:Nr+dk_max,iray,j,nact,id) = chi_down(Nb+dk_min:Nr+dk_max,iray,j,nact,id) + &
! 					hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,iray,id) * abs(aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
! 				chi_up(Nb+dk_min:Nr+dk_max,iray,i,nact,id) = chi_up(Nb+dk_min:Nr+dk_max,iray,i,nact,id) + &
! 					hc_fourPI * aatom%lines(kc)%Bij * aatom%lines(kc)%phi_loc(:,iray,id) * abs(aatom%n(i,icell)*wj/wi - aatom%lines(kc)%gij*aatom%n(j,icell))
! 			
! 			enddo
! 			
! ! 			if (iray > 1) cycle aatom_loop
! ! 			!only this part for the first ray otherwise go to the next atom
! 					
! 		enddo aatom_loop
! 	
! 	return
! 	end subroutine cross_coupling_terms
	
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
