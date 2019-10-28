MODULE init_solution

 use statequil_atoms, only : FillGamma_atom_zero_radiation, initGamma_atom, SEE_atom, &
 						     FillGamma_atom_mali, FillGamma_Atom, FillGamma_atom_hogerheijde, &
 						     FillGamma_atom_mu, FillGamma_atom_mali_mu, FillGamma_atom_hogerheijde_mu
 use atmos_type, only : atmos
 use atom_type, only : AtomType
 use spectrum_type, only : alloc_phi_lambda, NLTEspec, dealloc_phi_lambda, alloc_weights, dealloc_weights
 use mcfost_env, only : dp
 use collision, only : openCollisionFile, closeCollisionFile, keep_collision_lines
 use accelerate
 use messages
 use parametres
 use math
 
 IMPLICIT NONE
 
 real(kind=dp), dimension(:,:,:), allocatable :: gpop_old, pop_old
 real(kind=dp), dimension(:), allocatable :: flatpops
 
  
 CONTAINS
 
 SUBROUTINE Init_NLTE(sub_iterations_enabled)
 ! ---------------------------------------------------- !
  ! Initialize quantities for NLTE atomic line transfer
  ! set up pointers to the right FillGamma_bb and _bf.
 ! ---------------------------------------------------- !
  type (AtomType), pointer :: atom
  logical, intent(in), optional :: sub_iterations_enabled
  integer :: nact, Nmaxlevel, icell
   
   !alloc space for line profiles for each line: line%phi(line%Nlambda, Nrays, Nproc)
   CALL alloc_phi_lambda()
   CALL alloc_weights()
   
  if (lNg_acceleration) CALL Warning("Acceleration enabled", "Allocating space for Structure")
  ! if OLD_POPULATIONS, the init is done at the reading
  ! Or this subroutine should be invoked after readatom, because we need pops for ne
  
   NmaxLevel = 0
   do nact=1,atmos%Nactiveatoms
     atom => atmos%ActiveAtoms(nact)%ptr_atom
     allocate(atom%Gamma(atom%Nlevel,atom%Nlevel,NLTEspec%NPROC))
     !Now we can set it to .true. The new background pops or the new ne pops
     !will used the H%n
     !!Force it
     !!atom%initial_solution = "ZERO_RADIATION"
     atom%NLTEpops = .true.
     NmaxLevel = max(NmaxLevel, atom%Nlevel)
     write(*,*) "Setting initial solution ", atom%initial_solution," for active atom ", atom%ID, atom%active
     SELECT CASE (atom%initial_solution)
      CASE ("LTE_POPULATIONS")
        atom%n(:,:) = atom%nstar(:,:)
        
      CASE ("ZERO_RADIATION")
        do icell=1,atmos%Nspace
         if (atmos%icompute_atomRT(icell)>0) then
           !atom level version of initGamma and FillGamma and updatePopulations
           CALL initGamma_atom(1, icell, atom)
           !no rays integration here
           !!CALL FillGamma_bf_zero_radiation(1, icell, atom, 1)
           !the next one is cell independent but we init Gamma at each cell
           !even if the cell index appears, it is just here for consistance with the others FillGamma_bb_XXXX
           !!CALL FillGamma_bb_zero_radiation(1, icell, atom, 1) !no parallèle yet
           CALL FillGamma_atom_zero_radiation(1, icell, atom)
           CALL SEE_atom(1, icell,atom)
!            write(*,*) "n(zerR)", atmos%ActiveAtoms(1)%ptr_atom%n(:,icell)
!            write(*,*) "nstar", atmos%ActiveAtoms(1)%ptr_atom%nstar(:,icell)
         end if
        enddo
     !CASE ("OLD_POPULATIONS")
     !CASE ("ITERATIONS")
     !CASE ("SOBOLEV")
     CASE DEFAULT
      CALL Error("Initial solution unkown or not implemented", atom%initial_solution) 
     END SELECT 

     CALL openCollisionFile(atom) !closed at the end of the NLTE, it is not nice to do that
       								!but cheap in memory. Cause problem if parallel or
       								!run on a machine. Extra I/O overhead expected

     !CALL allocNetCoolingRates(atmos%ActiveAtoms(nact)%ptr_atom)
     if (lNg_acceleration) then !Nord >0 already tested
      Write(*,*) " Allocating space for Ng structure"
      CALL initNg(atmos%Nspace*atom%Nlevel, iNg_Ndelay, iNg_Norder , iNg_Nperiod, atom%Ngs)
     endif

	 CALL Keep_collision_lines(atom) !an array containing the lines in file read from atom%offset_col to END
	 !it avoids reading in the file, instead it reads an array (small)
	 CALL closeCollisionFile(atom)
    !!CALL writeAtomData(atmos%ActiveAtoms(nact)%ptr_atom) !to move elsewhere
  	 !deallocate(atom%C) !not used anymore if stored on RAM
     atom => NULL()
   enddo
   
   if (present(sub_iterations_enabled)) then
    if (sub_iterations_enabled) then
      Write(*,*) " Allocating space for sub-iterations populations"
      allocate(pop_old(atmos%NactiveAtoms, NmaxLevel, NLTEspec%NPROC)); pop_old = 0d0
    endif
   endif
   allocate(gpop_old(atmos%NactiveAtoms, Nmaxlevel,atmos%Nspace)); gpop_old = 0d0

  if (atmos%include_xcoupling) then
   write(*,*) "Using MALI"
    FillGamma_atom => FillGamma_atom_mali
    FillGamma_atom_mu => FillGamma_atom_mali_mu
  else
   write(*,*) "Using Hogerheijde"
    FillGamma_atom => FillGamma_atom_hogerheijde
    FillGamma_atom_mu => FillGamma_atom_hogerheijde_mu
  end if! 

   
 RETURN
 END SUBROUTINE init_nlte
 
 SUBROUTINE free_nlte_sol(sub_iterations_enabled)
  logical, intent(in) :: sub_iterations_enabled
  type (AtomType), pointer :: atom
  integer :: nact
  
!   FillGamma_bb => NULL()
!   FillGamma_bf => NULL() !not needed anymore at the end of the NLTE
 
  deallocate(gpop_old)
  if (sub_iterations_enabled) then
   !if (allocated(pop)) deallocate(pop)
   if (allocated(pop_old)) deallocate(pop_old)
  endif
  CALL dealloc_phi_lambda() !do not deallocate if kept in profile
  CALL dealloc_weights()
  
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