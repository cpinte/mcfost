! ------------------------------------------------------------------- !
! ------------------------------------------------------------------- !
! Module that solves for electron density for given
! temperature grid, Hydrogen total populations
! elements and their abundance and their LTE or NLTE populations.
!
!
! If ne_intial_slution is "NEMODEL" then first guess is computed using 
! the value stored in ne as first guess.
! This allows to iterate the electron density through the NLTE scheme.
!
! See: Hubeny & Mihalas (2014)
!         "Theory of Stellar Atmospheres",
!                         from p. 588 to p. 593
!
!
! ------------------------------------------------------------------- !
! ------------------------------------------------------------------- !

MODULE solvene

 use atmos_type, only : Nelem, Elements, Atoms, Natom, nHtot, ne, T, icompute_atomRT, ntotal_atom
 use atom_type, only : Element, AtomType
 use math, only : interp1D
 use constant
 use lte
 !use accelerate, only : initNg, freeNg, NgAcceleration
 use messages, only : Error, Warning
 use input, only : nb_proc
 use mcfost_env, only : dp
 use constantes, only : tiny_dp, huge_dp
 use parametres, only : n_cells
 !$ use omp_lib

 IMPLICIT NONE

 real(kind=dp), parameter :: MAX_ELECTRON_ERROR=1d-6
 integer, parameter :: N_MAX_ELECTRON_ITERATIONS=50
 integer, parameter :: N_MAX_ELEMENT=26 !100 is the max
 real(kind=dp), parameter :: ne_min_limit = 1d-50 !if ne < ne_min_limt, ne = 0
											!tiny_dp
											!if this limit is to low, f*1/ne could create nan or infinity
 CONTAINS

 ! ----------------------------------------------------------------------!
 ! do not forget to use the good value of parition function and potential
 ! parition function read in logarithm and used is 10**(U)
 ! potential in cm-1, converted in the routines in J.
 ! ----------------------------------------------------------------------!

 SUBROUTINE ne_Hionisation0 (k, U0, U1, ne)
 ! ----------------------------------------------------------------------!
  ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
  ! of H.
  ! Njl = Nj1l * ne * phi_jl
  ! ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
  ! eq. 17.77 and 17.78
  ! ne(H) = (sqrt(N*phiH + 1)-1)/phiH
 ! ----------------------------------------------------------------------!

  integer, intent(in) :: k
  real(kind=dp), intent(in) :: U0, U1
  real(kind=dp) :: phiH
  real(kind=dp), intent(out) :: ne
  phiH = phi_jl(k, U0, U1, Elements(1)%ptr_elem%ionpot(1))

  !if no free e-, Nt = NH + NH+ with NH+=ne
  if (phiH>0) then
   ne = (sqrt(nHtot(k)*phiH*4. + 1)-1)/(2.*phiH) !without H minus
  else
   ne = 0d0
  endif

 RETURN
 END SUBROUTINE ne_Hionisation0

 SUBROUTINE ne_Hionisation (k, U0, U1, ne)
 ! ----------------------------------------------------------------------!
  ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
  ! of H.
  ! Njl = Nj1l * ne * phi_jl
  ! ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
  ! eq. 17.77 and 17.78
  ! ne(H) = (sqrt(N*phiH + 1)-1)/phiH
 ! ----------------------------------------------------------------------!

  integer, intent(in) :: k
  real(kind=dp), intent(in) :: U0, U1
  real(kind=dp) :: phiH
  real(kind=dp), intent(out) :: ne
  
  phiH = phi_jl(k, U0, U1, Elements(1)%ptr_elem%ionpot(1))

  !if no free e-, Nt = NH + NH+ with NH+=ne
  !ne = (sqrt(nHtot(k)*phiH*4. + 1)-1)/(2.*phiH) !without H minus
  !if free e-, Nt=Nhot + NH+ + ne
  if (phiH>0) then
   ne = (sqrt(nHtot(k)*phiH + 1)-1)/(phiH)
  else
   ne = 0d0
  endif

 RETURN
 END SUBROUTINE ne_Hionisation

 SUBROUTINE ne_Metal(k, U0, U1, chi, A, ne)
 ! ----------------------------------------------------------------------!
  ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
  ! of a single metal.
 ! ----------------------------------------------------------------------!

  integer, intent(in) :: k
  real(kind=dp), intent(in) :: U0, U1, chi, A
  real(kind=dp) :: phiM, alphaM
  real(kind=dp), intent(out) :: ne

  phiM = phi_jl(k, U0, U1, chi)
  alphaM = A ! relative to H, for instance 1.-6 etc
  
  if (phiM > 0) then
   ne = (sqrt(alphaM*nHtot(k)*phiM +0.25*(1+alphaM)**2) - 0.5*(1+alphaM) ) / phiM
  else
   ne = 0d0
  endif
  
 RETURN
 END SUBROUTINE ne_Metal


!  FUNCTION getPartitionFunctionk(elem, stage, k) result (Uk)
!  ! ----------------------------------------------------------------------!
!   ! Interpolate the partition function of Element elem in ionisation stage
!   ! stage at cell point k.
!  ! ----------------------------------------------------------------------!
! 
!   type(Element) :: elem
!   integer, intent(in) :: stage, k
!   real(kind=dp) :: Uk, part_func(Npf)
!   
!   part_func = elem%pf(stage,:)
!   Uk = Interp1D(Tpf,part_func,T(k))
!   !!Uk = (10.d0)**(Uk) !! 29/12/2019 changed log10 in log
!   Uk = exp(Uk)
!  RETURN
! END FUNCTION getPartitionFunctionk

subroutine calc_ionisation_frac(elem, k, ne, fjk, dfjk)

	! ------------------------------------------------------------------------!
	! fractional population f_j(ne,T)=N_j/N for element Elem
	! and its partial derivative with ne. If Elem is an element with
	! detailed model and if NLTE populations for it exist, their are used
	! instead of LTE.
	!
	! to test: If detailed_model is not applied only to active atoms
	! it means that some atoms will contribute their LTE values, therefore
	! they have to be updated. Is it okey or not ?
	! OR it better to finder a consistant solution with LTE atoms (elements or not)
	! ------------------------------------------------------------------------!

	real(kind=dp), intent(in) :: ne
	integer, intent(in) :: k
	type (Element), intent(in) :: Elem
	real(kind=dp), dimension(:), intent(inout) :: fjk, dfjk
	real(kind=dp) :: Uk, Ukp1, sum1, sum2
	logical :: detailed_model
	integer :: j, i

	detailed_model = .false.

	if (associated(elem%model)) then
		!otherwise, the first call (wihtout knowing ne at all)
		!would not work since atom%n = 0.0 (or atom%nstar)
		!.true. before the maxval
		!detailed_model = (maxval(elem%model%n)>0.0).and.(elem%model%active)
		!-> two much time to test n > 0 for large model
		detailed_model = (elem%model%NLTEpops .and. elem%model%active)
		!can add condition on active or not, but check bckgr opac and eval of lte pops after
	endif
  
	fjk(:) = 0d0
	dfjk(:) = 0d0

	if (detailed_model) then
		if (.not.elem%model%active) call warning("Beware passive atoms used!")
		fjk = 0d0
		dfjk = 0d0

		do i = 1, elem%model%Nlevel
			fjk(elem%model%stage(i)+1) = fjk(elem%model%stage(i)+1)+(elem%model%stage(i))*elem%model%n(i,k)
		end do  
		                                     
		fjk(:) = fjk(:)/(nHtot(k)*elem%model%Abund)
		
	else !no model use LTE
		fjk(1)=1.
		dfjk(1)=0.
		sum1 = 1.
		sum2 = 0.
		Uk = getPartitionFunctionk(elem,1,k)
		do j=2,Elem%Nstage
			Ukp1 = getPartitionFunctionk(elem,j,k)
			fjk(j) = Sahaeq(k,fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)

! 			if (ne>0) then
! 				dfjk(j) = -(j-1)*fjk(j)/ne
! 			else
! 				dfjk(j) = 0d0
! 			endif
			!-> fjk is zero if ne=0 in Sahaeq
			dfjk(j) = -(j-1)*fjk(j)/( ne + tiny_dp)
    
			sum1 = sum1 + fjk(j)
			sum2 = sum2 + dfjk(j)

			Uk = Ukp1
		end do
		
		fjk(:)=fjk(:)/sum1
		dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1

	endif

return
end subroutine calc_ionisation_frac
	
subroutine solve_ne_nlte_loc(k, ne_initial_solution, dne)
! for non-lte loop
!there is no iteration on the electron donors at LTE
	integer :: k
	real(kind=dp), intent(inout) :: dne
	integer, parameter :: Nmaxstage = 50
	character(len=20) :: ne_initial_solution
	real(kind=dp) :: ne_old, akj, ne_new
	real(kind=dp):: PhiHmin
	real(kind=dp), dimension(Nmaxstage) :: fjk, dfjk
	integer :: n, j, i
	type (Element), pointer :: elem

	ne_old = ne(k)
	ne_new = 0.0_dp
	
	do n=1, N_MAX_ELEMENT

		elem => Elements(n)%ptr_elem
		
		call calc_ionisation_frac(elem, k, ne_old, fjk, dfjk)

		if (n.eq.1)  then ! H minus for H
			PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
       ! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
			ne_new = ne_new + ne_old*fjk(1)*PhiHmin
		end if
		
		do j=2,elem%Nstage
			akj = elem%abund*(j-1) !(j=0 for neutrals, 1 for signly ionised etc..)
			ne_new = ne_new + akj * fjk(j)
		end do
		
		elem => NULL()
	end do !loop over elem

	ne_new = ne_new * nHtot(k)
	dne = dabs(1.0_dp - ne_old/ne(k))
	
	!update for this cell
	ne(k) = ne_new

return
end subroutine solve_ne_nlte_loc

!try without para ! for epsilon
subroutine solve_electron_density(ne_initial_solution, verbose, epsilon)
	! ----------------------------------------------------------------------!
	! Solve for electron density for a set of elements
	! stored in Elements. Elements up to N_MAX_ELEMENT are
	! used. If an element has also an atomic model, and if for
	! this element NLTE populations are known these pops are
	! used to compute the electron density. Otherwise, LTE is used.
	! 
	! When ne_initial_solution is set to HIONISA, uses
	! sole hydgrogen ionisation to estimate the initial ne density.
	! If set to NPROTON or NEMODEL, protons number or density
	! read from the model are used. Note that NPROTON supposes
	! that NLTE populations are present for hydrogen, since
	! nprot = hydrogen%n(Nlevel,:).
	!
	! Also returns epsilon, the change in electron density from the intial guess
	! to the converged values. THIS IS NOT the convergence threshold used inside
	! the electron loop.
	!
	!
	! epsilon is shared, thats' okay ?
	! ----------------------------------------------------------------------!
	character(len=20), intent(in) :: ne_initial_solution
	logical, intent(in) :: verbose
	real(kind=dp), intent(inout) :: epsilon !difference wrt the initial solution, not criterion of convergence (dne)
	real(kind=dp) :: error, ne_old, akj, sum, Uk, Ukp1, dne, ne0
	real(kind=dp):: ne_oldM, UkM, PhiHmin
	real(kind=dp), dimension(50) :: fjk, dfjk
	integer :: n, k, niter, j, ZM, id
	type (Element), pointer :: elem
	integer :: unconverged_cells(nb_proc)

	!!if (verbose) write(*,*) "Initial solution for electron loop:", ne_initial_solution

	unconverged_cells(:) = 0
	id = 1

	epsilon = 0.0

	!$omp parallel &
	!$omp default(none) &
	!$omp private(k,n,j,fjk,dfjk,ne_old,niter,error,sum,PhiHmin,Uk,Ukp1,ne_oldM) &
	!$omp private(dne, akj, id, ne0, elem) &
	!$omp shared(n_cells, Elements, ne_initial_solution,Hydrogen, ZM, unconverged_cells, Nelem) &
	!$omp shared(ne, T, icompute_atomRT, nHtot, epsilon)
	!$omp do
	do k=1,n_cells
		!$ id = omp_get_thread_num() + 1
		if (icompute_atomRT(k) <= 0) CYCLE !transparent or dark

		!Initial solution for this cell

		if (ne_initial_solution.eq."N_PROTON") then
			ne_old = Hydrogen%n(Hydrogen%Nlevel,k)!np(k)
		else if (ne_initial_solution.eq."NE_MODEL") then
			ne_old = ne(k)
		else !"H_IONISATION" or unkown
   
			!Initial solution ionisation of H
			elem => Elements(1)%ptr_elem
			Uk = getPartitionFunctionk(elem, 1, k)
			Ukp1 = getPartitionFunctionk(elem, 2, k)
			elem => NULL()

			if (Ukp1 /= 1.0_dp) then 
				CALL Warning("Partition function of H+ should be 1")
				write(*,*) Uk, Ukp1
				stop
			end if
			CALL ne_Hionisation (k, Uk, Ukp1, ne_old)

			!Metal
			Zm = 11
			elem => Elements(ZM)%ptr_elem
			Uk = getPartitionFunctionk(Elements(ZM)%ptr_elem, 1, k)
			Ukp1 = getPartitionFunctionk(Elements(ZM)%ptr_elem, 2, k)
			CALL ne_Metal(k, Uk, Ukp1, elements(ZM)%ptr_elem%ionpot(1),elements(ZM)%ptr_elem%Abund, ne_oldM)
			!write(*,*) "neMetal=",ne_oldM
			!if Abund << 1. and chiM << chiH then
			! ne (H+M) = ne(H) + ne(M)
			ne_old = ne_old + ne_oldM
    

			!to avoid having very low values, between say tiny_dp and 1e-100
			!just set to 0 if < 0: crude ?
			if (ne_oldM < ne_min_limit) ne_oldM = 0d0
			if (ne_old < ne_min_limit) ne_old = 0d0
			
		end if

		!Loop starts
		ne(k) = ne_old
		ne0 = ne_old
		niter=0
		do while (niter < N_MAX_ELECTRON_ITERATIONS)
			error = ne_old/nHtot(k)
			sum = 0.0

			do n=1,Nelem
				if (n > N_MAX_ELEMENT) exit
				elem => Elements(n)%ptr_elem

				call calc_ionisation_frac(elem, k, ne_old, fjk, dfjk)

				if (n.eq.1)  then ! H minus for H
					!2 = partition function of HI, should be replace by getPartitionFunctionk(elements(1)%ptr_elem, 1, icell)
					PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
					! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
					error = error + ne_old*fjk(1)*PhiHmin
					sum = sum-(fjk(1)+ne_old*dfjk(1))*PhiHmin
				end if
				!neutrals do not contribute right
				
				do j=2, elem%Nstage
					akj = elem%Abund*(j-1) !because j starts at 0 for neutrals, 1 for singly ionised etc
					error = error -akj*fjk(j)
					sum = sum + akj*dfjk(j)
				end do
				
				elem => NULL()
			end do !loop over elem

			ne(k) = ne_old - nHtot(k)*error / (1.-nHtot(k)*sum)
			!!dne = abs((ne(k)-ne_old)/(ne_old+tiny_dp))
			dne = abs ( 1.0_dp - ne_old / ne(k) )
			ne_old = ne(k)
			
			!can I do that in parallel ? epsilon is shared
			epsilon = max(epsilon, abs(1.0_dp - ne0 / ne(k)))
    
			if (is_nan_infinity(ne(k))>0) then
				write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, sum
				stop
			else if (ne(k) <= 0.0) then
				write(*,*) "ne <= 0"
				write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, sum	
				stop
			endif
			niter = niter + 1
			if (dne <= MAX_ELECTRON_ERROR) then
				if (ne(k) < ne_min_limit) ne(k) = 0d0 !tiny_dp or < tiny_dp
				exit
			else if (niter >= N_MAX_ELECTRON_ITERATIONS) then
				if (dne >= 1) then !shows the warning only if dne is actually greater than 1
					CALL Warning("Electron density has not converged for this cell")
					write(*,*) "icell=",k,"maxIter=",N_MAX_ELECTRON_ITERATIONS,"dne=",dne,"max(err)=", MAX_ELECTRON_ERROR, "ne=",ne(k), "T=",T(k)," nH=",nHtot(k)
					!set articially ne to some value ?
					ne(k) = ne_oldM !already tested if < ne_min_limit
					write(*,*) "Setting ne to ionisation of ", elements(ZM)%ptr_elem%ID, ne_oldM
				end if
				unconverged_cells(id) = unconverged_cells(id) + 1 !but count each cell for with dne > MAX_ELECTRON_ERROR
			end if !convergence test
    
		end do !while loop
	end do !loop over spatial points
	!$omp end do
	!$omp end parallel

	if (verbose) then
			write(*,*) " ------------------------------------------------ "
			write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
			write(*,'("   >>>  epsilon="(1ES17.8E3))') epsilon
			write(*,*) " ------------------------------------------------ "
	endif

  !that's because I use the sum as a variable so the function doesn't exist.
  do k=2, nb_proc
   unconverged_cells(1) = unconverged_cells(1) + unconverged_cells(k)
  end do
  if (unconverged_cells(1) > 0) write(*,*) "Found ", unconverged_cells(1)," unconverged cells (%):", &
  													100*real(unconverged_Cells(1))/real(n_cells)
  
return
end subroutine solve_electron_density

 SUBROUTINE getfjk (Elem, ne, k, fjk, dfjk)
 ! ----------------------------------------------------------------------!
  ! fractional population f_j(ne,T)=N_j/N for element Elem
  ! and its partial derivative with ne. If Elem is an element with
  ! detailed model and if NLTE populations for it exist, their are used
  ! instead of LTE.
 ! ----------------------------------------------------------------------!

  real(kind=dp), intent(in) :: ne
  integer, intent(in) :: k
  type (Element), intent(in) :: Elem
  real(kind=dp), dimension(:), intent(inout) :: fjk, dfjk
  real(kind=dp) :: Uk, Ukp1, sum1, sum2
  logical :: has_nlte_pops
  integer :: j, i

  has_nlte_pops = .false. !reset for all atom

  if (associated(elem%model)) then
  	!it also check that it is active ?
  	!because a model can be associated but it doesn't mean it has nlte pops !
  	if (elem%model%NLTEpops) then
  		has_nlte_pops = .true.
  	endif
  endif
  
  fjk(:) = 0d0
  dfjk(:) = 0d0

  !may be active without NLTEpops or passive with read NLTE pops
  if (has_nlte_pops) then
   fjk = 0d0
   dfjk = 0d0
   !For Nlevel, Nlevel stages
   !fj = Nj/Ntot
   !first, Nj = Sum_i nij
   do i = 1, elem%model%Nlevel
    fjk(elem%model%stage(i)+1) = fjk(elem%model%stage(i)+1)+(elem%model%stage(i))*elem%model%n(i,k)
   end do                                        !is there a + 1 here so that fjk(1)=1*atom%n(1)/Ntot
   !Divide by Ntotal and retrieve fj = Nj/N for all j
   fjk(:) = fjk(:)/(nHtot(k)*elem%model%Abund)
  else !not active or active but first iteration of the NLTEloop so that
  	!NLTEpos has been set to .false., whateveeer, use LTE
   fjk(1)=1.
   dfjk(1)=0.
   sum1 = 1.
   sum2 = 0.
   Uk = getPartitionFunctionk(elem,1,k)
   do j=2,Elem%Nstage !-->vectorize ?
    Ukp1 = getPartitionFunctionk(elem,j,k)
    ! fjk(j) = Nj / Sum_j Nj
    ! Nj = Nj-1/(phi*ne) Saha Equation
    ! fj = Nj/N = Nj/N0 / N/N0
    ! -> first computes Nj/N0 using Nj-1/N0 relative to N0
    ! -> sum up = 1 + N1/N0 + Nj/N0
    ! -> divide Nj/N0 by this sum and retrive Nj/N for all j
    !  --> Nj/N0 / (1+Nj/N0) = Nj/(N0*(1+Nj/N0)) = Nj / (N0 + Nj) = fj
!-> Check that again
    fjk(j) = Sahaeq(k,fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)
    !write(*,*) "j=",j," fjk=",fjk(j)
    !write(*,*) fjk(j-1), fjk(j), fjk(j)/(ne+tiny_dp)
    
    !why without this I have nan ???
    if (ne>0) then
    dfjk(j) = -(j-1)*fjk(j)/ne  !-> fjk(j) should be 0 if ne = 0
     			     			 !see Saheq.
     !dfjk(j) = -(j-1) * fjk(j) * min(huge_dp,1/ne)
    else
     dfjk(j) = 0d0
    endif
    
    sum1 = sum1 + fjk(j)
    sum2 = sum2 + dfjk(j)
    
    !write(*,*) "j=",j," dfjk=",dfjk(j), "fjk=", fjk(j), sum1, sum2


    Uk = Ukp1
   end do
   fjk(:)=fjk(:)/sum1
   dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1
   
   if (any_nan_infinity_vector(dfjk)>0) then
       write(*,*) "j=",j," dfjk=",dfjk(j), "fjk=", fjk(j), ne, T(k), sum1, sum2
   stop
   endif
  end if

 RETURN
 END SUBROUTINE getfjk

 SUBROUTINE Solve_Electron_Density_old(ne_initial_solution)
 ! ----------------------------------------------------------------------!
  ! Solve for electron density for a set of elements
  ! stored in Elements. Elements up to N_MAX_ELEMENT are
  ! used. If an element has also an atomic model, and if for
  ! this element NLTE populations are known these pops are
  ! used to compute the electron density. Otherwise, LTE is used.
  ! 
  ! When ne_initial_solution is set to HIONISA, uses
  ! sole hydgrogen ionisation to estimate the initial ne density.
  ! If set to NPROTON or NEMODEL, protons number or density
  ! read from the model are used. Note that NPROTON supposes
  ! that NLTE populations are present for hydrogen, since
  ! nprot = hydrogen%n(Nlevel,:).
  ! If keyword is not set, HIONISATION is used.
 ! ----------------------------------------------------------------------!
  character(len=20), optional :: ne_initial_solution
  character(len=20) :: initial
  real(kind=dp) :: error, ne_old, akj, sum, Uk, dne, Ukp1
  real(kind=dp):: ne_oldM, UkM, PhiHmin
!   real(kind=dp), dimension(n_cells) :: np
  real(kind=dp), dimension(:), allocatable :: fjk, dfjk
  integer :: Nmaxstage=0, n, k, niter, j, ZM, id, ninit, nfini
  integer :: unconverged_cells(nb_proc)

  if (.not. present(ne_initial_solution)) then
      initial="H_IONISATION"!use HIONISAtion
  else
    initial=ne_initial_solution
!     if ((initial .neq. "N_PROTON") .neq. (initial /= "NE_MODEL")) then
!      write(*,*) 'NE initial solution unkown, set to H_IONISATION'
!      initial = "H_IONISATION"
!     end if
  end if
	write(*,*) "Initial solution for electron loop:", initial

  unconverged_cells(:) = 0
  id = 1
  do n=1,Nelem
   Nmaxstage=max(Nmaxstage,Elements(n)%ptr_elem%Nstage)
  end do
  allocate(fjk(Nmaxstage))
  allocate(dfjk(Nmaxstage))

  !np is the number of protons, by default the last level
  !of Hydrogen.
  
  !Note that will raise errors if np is 0d0 because Hydrogen pops are not
  !known. np is known if:
  ! 1) NLTE populations from previous run are read
  ! 2) The routine is called in the NLTE loop meaning that all H levels are known

!   if (initial.eq."N_PROTON") &
!      np=Hydrogen%n(Hydrogen%Nlevel,:)
write(*,*) "Test Nelem",Nelem

  !$omp parallel &
  !$omp default(none) &
  !$omp private(k,n,j,fjk,dfjk,ne_old,niter,error,sum,PhiHmin,Uk,Ukp1,ne_oldM) &
  !$omp private(dne, akj, id, ninit, nfini) &
  !$omp shared(n_cells, Elements, initial,Hydrogen, ZM, unconverged_cells, Nelem) &
  !$omp shared(ne, T, icompute_atomRT, nHtot)
  !$omp do
  do k=1,n_cells
   !$ id = omp_get_thread_num() + 1
   if (icompute_atomRT(k) <= 0) CYCLE !transparent or dark


   !write(*,*) "The thread,", omp_get_thread_num() + 1," is doing the cell ", k
   if (initial.eq."N_PROTON") then
    ne_old = Hydrogen%n(Hydrogen%Nlevel,k)!np(k)
   else if (initial.eq."NE_MODEL") then
    ne_old = ne(k)
   else !"H_IONISATION" or unkown
   
    !Initial solution ionisation of H
    Uk = getPartitionFunctionk(Elements(1)%ptr_elem, 1, k)
    Ukp1 = getPartitionFunctionk(Elements(1)%ptr_elem, 2, k)

    if (Ukp1 /= 1d0) then 
     CALL Warning("Partition function of H+ should be 1")
     write(*,*) Uk, Ukp1
     stop
    end if
    CALL ne_Hionisation (k, Uk, Ukp1, ne_old)


    if (T(k) >= 20d3) then
      ZM = 2 !He
    !else if (T(k) <= 1000d0) then
    else
      ZM = 11 ! Na !trade off between higher A and lower ionpot
    !else
    !  ZM = 26 ! Fe
    end if    
    !add Metal
    Uk = getPartitionFunctionk(Elements(ZM)%ptr_elem, 1, k)
    Ukp1 = getPartitionFunctionk(Elements(ZM)%ptr_elem, 2, k)
    CALL ne_Metal(k, Uk, Ukp1, elements(ZM)%ptr_elem%ionpot(1),elements(ZM)%ptr_elem%Abund, ne_oldM)
    !write(*,*) "neMetal=",ne_oldM
    !if Abund << 1. and chiM << chiH then
    ! ne (H+M) = ne(H) + ne(M)
    ne_old = ne_old + ne_oldM
    
    !!!!Check here !!!
    !to avoid having very low values, between say tiny_dp and 1e-100
    !just set to 0 if < 0: crude ?
    if (ne_oldM < ne_min_limit) ne_oldM = 0d0
    if (ne_old < ne_min_limit) ne_old = 0d0
    !!!!!!!!!!!!!!!!!
   end if

   ne(k) = ne_old
   niter=0
   do while (niter < N_MAX_ELECTRON_ITERATIONS)
    error = ne_old/nHtot(k)
    sum = 0.

    !ninit = (1. * (id-1)) / nb_proc * Nelem + 1
    !nfini = (1. * id) / nb_proc * Nelem
    do n=1,Nelem
     if (n > N_MAX_ELEMENT) exit

     CALL getfjk(Elements(n)%ptr_elem,ne_old,k,fjk,dfjk)!

     if (n.eq.1)  then ! H minus for H
       !2 = partition function of HI, should be replace by getPartitionFunctionk(elements(1)%ptr_elem, 1, icell)
       PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
       ! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
       error = error + ne_old*fjk(1)*PhiHmin
       sum = sum-(fjk(1)+ne_old*dfjk(1))*PhiHmin
       !write(*,*) "phiHmin=",PhiHmin,error, sum
     end if
     !neutrals do not contribute right
     do j=2,elements(n)%ptr_elem%Nstage
      akj = elements(n)%ptr_elem%Abund*(j-1) !because j starts at 0 for neutrals, 1 for singly ionised etc
      error = error -akj*fjk(j)
      sum = sum + akj*dfjk(j)
      !write(*,*) n-1, j-1, akj, error, sum, dfjk(j)
     end do
    end do !loop over elem

    ne(k) = ne_old - nHtot(k)*error /&
          (1.-nHtot(k)*sum)
    dne = dabs((ne(k)-ne_old)/(ne_old+tiny_dp))
    ne_old = ne(k)
    
    if (is_nan_infinity(ne(k))>0) then
    write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), &
              "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, sum
    stop
    endif
    niter = niter + 1
    if (dne <= MAX_ELECTRON_ERROR) then
      !write(*,*) niter, "icell=",k," T=",T(k)," ne=",ne(k), " dne=", dne
      if (ne(k) < ne_min_limit) ne(k) = 0d0 !tiny_dp or < tiny_dp
     exit
    else if (niter >= N_MAX_ELECTRON_ITERATIONS) then
      if (dne >= 1) then !shows the warning only if dne is actually greater than 1
          CALL Warning("Electron density has not converged for this cell")
          write(*,*) "icell=",k,"maxIter=",N_MAX_ELECTRON_ITERATIONS,"dne=",dne,"max(err)=", MAX_ELECTRON_ERROR, &
          "ne=",ne(k), "T=",T(k)," nH=",nHtot(k)
          !set articially ne to some value ?
          ne(k) = ne_oldM !already tested if < ne_min_limit
          write(*,*) "Setting ne to ionisation of ", elements(ZM)%ptr_elem%ID, ne_oldM
!       else
!           write(*,*) "icell=",k,"maxIter=",N_MAX_ELECTRON_ITERATIONS,"dne=",dne,"max(err)=", MAX_ELECTRON_ERROR, &
!           "ne=",ne(k), "T=",T(k)," nH=",nHtot(k)
      end if
      unconverged_cells(id) = unconverged_cells(id) + 1 !but count each cell for with dne > MAX_ELECTRON_ERROR
    end if !convergence test
    
   end do !while loop
  end do !loop over spatial points
  !$omp end do
  !$omp end parallel

  deallocate(fjk, dfjk)
  
  write(*,*) "Maximum/minimum Electron density in the model (m^-3):"
  write(*,*) MAXVAL(ne),MINVAL(ne,mask=icompute_atomRT>0)
  
  !that's because I use the sum as a variable so the function doesn't exist.
  do k=2, nb_proc
   unconverged_cells(1) = unconverged_cells(1) + unconverged_cells(k)
  end do
  if (unconverged_cells(1) > 0) write(*,*) "Found ", unconverged_cells(1)," unconverged cells (%):", &
  													100*real(unconverged_Cells(1))/real(n_cells)
  
 RETURN
 END SUBROUTINE Solve_Electron_Density_old


END MODULE solvene
