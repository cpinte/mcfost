MODULE init_solution

	use statequil_atoms, only						: collision_matrix_atom, initGamma_atom, SEE_atom
 						     
	use atmos_type, only							: atmos
	use atom_type, only								: AtomType, AtomicLine, AtomicContinuum
	use spectrum_type, only							: NLTEspec, alloc_weights, dealloc_weights, alloc_Jnu
	use mcfost_env, only							: dp
	use accelerate
	use messages
	use parametres
	use math
	use opacity, only : chi_loc, eta_loc
	use constant, only : CLIGHT
	use input!, only : ds, nb_proc
 

	IMPLICIT NONE

	real(kind=dp), dimension(:,:,:), allocatable :: gpop_old, Tex_old
	real(kind=dp), dimension(:), allocatable :: flatpops
	logical, dimension(:), allocatable :: lcell_converged

 
 
  
	CONTAINS
 
  
 !the factor 0.5 is here because dw is computed as w(k+1) - w(k-1) / 2.
 !and int(f) = sum(dw_k*f_k)
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
		real(kind=dp), dimension(line%Nlambda) :: wlam
		integer :: la, Nblue, Nred, la_start, la_end, la0
		real :: factor = 0.5
		real(kind=dp) :: norm !beware this is not the result of the integral 
  						   ! just a factor to convert dv/c from dlambda/lambda
   !la0: index of wavelengths on the frequency grid (size Nwaves). 
   !la:  index of wavelengths on the lambda grid of the line (size Nlambda).
   !la0 = Nblue - 1 + la; line expands from Nblue to Nred on the frequency grid.
   !la=1 <=> la0=Nblue; la=Nlambda <=> la0 = Nred = Nblue - 1 + Nlambda
   !dlambda = (lambda(la0 + 1) - lambda(la0 - 1)) * 0.5 <=> mean value.

                   !should be line%lambda instead of lambda0 ??
		norm = CLIGHT / line%lambda0 !because we want dv
		Nblue = line%Nblue; Nred = line%Nred
		la_start = 1; la_end = line%Nlambda

		wlam(1) = factor * (NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) * norm
  
		wlam(line%Nlambda) = factor * (NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) * norm
  !write(*,*) 1, wlam(1)
		do la=2,line%Nlambda-1

			la0 = Nblue - 1 + la
			wlam(la) = factor*(NLTEspec%lambda(la0 + 1)-NLTEspec%lambda(la0 - 1)) * norm
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
	! a 1 / h is missing
	! --------------------------------------------------------- !
		type(AtomicContinuum), intent(in) :: cont
		real(kind=dp), dimension(cont%Nlambda) :: wlam
		integer :: la, Nblue, Nred, la_start, la_end , la0
		real :: factor = 0.5

		Nblue = cont%Nblue; Nred = cont%Nred
		la_start = 1; la_end = cont%Nlambda

		wlam(1) = factor * &
			(NLTEspec%lambda(Nblue+1)-NLTEspec%lambda(Nblue)) / NLTEspec%lambda(Nblue)
  	
		wlam(cont%Nlambda) = factor * & 
			(NLTEspec%lambda(Nred)-NLTEspec%lambda(Nred-1)) / NLTEspec%lambda(Nred)
  	
		do la=2,cont%Nlambda-1
  
			la0 = Nblue - 1 + la
			wlam(la) = factor * &
			(NLTEspec%lambda(la0+1)-NLTEspec%lambda(la0-1)) / NLTEspec%lambda(la0)
   	
		end do

	RETURN
	END FUNCTION cont_wlam

	SUBROUTINE Init_NLTE(sub_iterations_enabled)
	! ---------------------------------------------------- !
	! Initialize quantities for NLTE atomic line transfer
	! set up pointers to the right FillGamma_bb and _bf.
	! ---------------------------------------------------- !
		type (AtomType), pointer :: atom
		logical, intent(in), optional :: sub_iterations_enabled
		integer :: nact, Nmaxlevel, icell, kr, alloc_status, l, NmaxTr
		real(kind=dp) :: dM, dpop
   
		CALL alloc_weights() !and phi_ray
		write(*,*) " Setting up wavelength integration weights for all transitions..."
		do nact=1,atmos%Nactiveatoms
			do kr=1, atmos%ActiveAtoms(nact)%ptr_atom%Nline
   
				atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%w_lam(:) = line_wlam(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr))
 
			enddo
			do kr=1, atmos%ActiveAtoms(nact)%ptr_atom%Ncont
   
				atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)%w_lam(:) = cont_wlam(atmos%ActiveAtoms(nact)%ptr_atom%continua(kr))

   
			 enddo  
		enddo

		!Total LTE opacity and emissivity that doesn't change during sub iterations, used to compute S_loc
		allocate(chi_loc(NLTEspec%Nwaves, atmos%Nrays, NLTEspec%Nproc), eta_loc(NLTEspec%Nwaves, atmos%Nrays, NLTEspec%Nproc))
		if (allocated(ds)) deallocate(ds)
		allocate(ds(atmos%Nrays, NLTEspec%Nproc))
   

		if (lNg_acceleration) CALL Warning("Acceleration enabled", "Allocating space for Structure")
  ! if OLD_POPULATIONS, the init is done at the reading
  ! Or this subroutine should be invoked after readatom, because we need pops for ne

		NmaxLevel = 0
		NmaxTr = 0
		do nact=1,atmos%Nactiveatoms
			atom => atmos%ActiveAtoms(nact)%ptr_atom
			allocate(atom%Gamma(atom%Nlevel,atom%Nlevel,NLTEspec%NPROC))
			allocate(atom%C(atom%Nlevel,atom%Nlevel,NLTEspec%NPROC)) 
     
     !Now we can set it to .true. The new background pops or the new ne pops
     !will used the H%n
     !!Force it
     !!atom%initial_solution = "ZERO_RADIATION"
			atom%NLTEpops = .true.
			NmaxLevel = max(NmaxLevel, atom%Nlevel)
			NmaxTr = max(NmaxTr, atom%Ncont + atom%Nline)
			write(*,*) "Setting initial solution ", atom%initial_solution," for active atom ", atom%ID, atom%active
			SELECT CASE (atom%initial_solution)
				CASE ("LTE_POPULATIONS")
					atom%n(:,:) = atom%nstar(:,:)

				CASE ("ZERO_RADIATION")
					dpop = 0.0_dp
					write(*,*) "ZERO_RAD, deactiavted presently"
					stop
!         do icell=1,atmos%Nspace
!          if (atmos%icompute_atomRT(icell)>0) then
!            !atom level version of initGamma and FillGamma and updatePopulations
!            CALL collision_matrix_atom(1,icell,atom)
!            CALL initGamma_atom(1, atom) !for this cell or id
!            !no rays integration here
!            !!CALL FillGamma_bf_zero_radiation(1, icell, atom, 1)
!            !the next one is cell independent but we init Gamma at each cell
!            !even if the cell index appears, it is just here for consistance with the others FillGamma_bb_XXXX
!            !!CALL FillGamma_bb_zero_radiation(1, icell, atom, 1) !no parallèle yet
!            CALL FillGamma_atom_zero_radiation(1, icell, atom)
!            CALL SEE_atom(1, icell,atom, dM)
! !            write(*,*) "***"
! !            write(*,*) atmos%nHtot(icell), sum(atom%nstar(:,icell)), sum(atom%n(:,icell))
! !            write(*,*) "nstar", (atom%nstar(l,icell),l=1,atom%Nlevel)
! !            write(*,*) "n", (atom%n(l,icell),l=1,atom%Nlevel)
! !            write(*,*) "***"
! 
!            dpop = max(dM, dpop)
!          end if
!         enddo
!stop
			CASE ("OLD_POPULATIONS")
				CALL Warning("OLD_POPULATIONS handled in readatom.f90 A.T.M")
			CASE DEFAULT
				CALL Error("Initial solution unkown or not implemented", atom%initial_solution)
			END SELECT

     !CALL allocNetCoolingRates(atmos%ActiveAtoms(nact)%ptr_atom)
			if (lNg_acceleration) then !Nord >0 already tested
				Write(*,*) " Allocating space for Ng structure"
				CALL initNg(atmos%Nspace*atom%Nlevel, iNg_Ndelay, iNg_Norder , iNg_Nperiod, atom%Ngs)
			endif

     !!CALL writeAtomData(atmos%ActiveAtoms(nact)%ptr_atom) !to move elsewhere
     !deallocate(atom%C) !not used anymore if stored on RAM
			atom => NULL()
		enddo
   
		!!Cannot yet remove cell converged, otherwise I break the randomness 
		!!allocate(lcell_converged(atmos%Nspace),stat=alloc_status)
		!!if (alloc_status > 0) call error("Allocation error lcell_converged")
		!!cell_converged(:) = .false.

		write(*,*) " ** Number max of levels among all atoms:", Nmaxlevel
		write(*,*) " ** Number max of transitions among all atoms:", NmaxTr

! 		if (present(sub_iterations_enabled)) then
! ! 			if (sub_iterations_enabled) then
! ! 				Write(*,*) " Allocating space for sub-iterations populations"
! ! 				allocate(pop_old(atmos%NactiveAtoms, NmaxLevel, NLTEspec%NPROC)); pop_old = 0d0
! !       !allocate(pops(atmos%NactiveAtoms, NmaxLevel, NLTEspec%NPROC)); pops = 0d0
! ! 			endif
! 		endif
		
		!One will be removed, depends on the global criteria of convergence
		allocate(gpop_old(atmos%NactiveAtoms, Nmaxlevel,atmos%Nspace)); gpop_old = 0d0
		allocate(Tex_old(atmos%NactiveAtoms, NmaxTr,atmos%Nspace)); Tex_old = 0d0


	RETURN
	END SUBROUTINE init_nlte

	SUBROUTINE free_nlte_sol(sub_iterations_enabled)
		logical, intent(in) :: sub_iterations_enabled
		type (AtomType), pointer :: atom
		integer :: nact

		deallocate(gpop_old, Tex_old)
  
		CALL dealloc_weights()
		deallocate(chi_loc, eta_loc, ds)
		if (allocated(lcell_converged)) deallocate(lcell_converged)


		do nact=1,atmos%Nactiveatoms
			atom  => atmos%ActiveAtoms(nact)%ptr_atom
   !!!!CALL closeCollisionFile(atom)
			if (allocated(atmos%ActiveAtoms(nact)%ptr_atom%gamma)) & !otherwise we have never enter the loop
				deallocate(atmos%ActiveAtoms(nact)%ptr_atom%gamma)
   !!if (allocated(atmos%ActiveAtoms(nact)%ptr_atom%Ckij)) deallocate(atmos%ActiveAtoms(nact)%ptr_atom%Ckij)
			if (lNg_acceleration) CALL freeNg(atom%Ngs)
			atom => NULL()
		end do

	RETURN
	END SUBROUTINE free_nlte_sol

END MODULE init_solution
