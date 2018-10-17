!
! This module solves for the radiative transfer equation by ray-tracing for 
! multi-level atoms.
!
! Outputs:
! - Flux (\lambda) [J.s^{-1}.m^{-2}.Hz^{-1}]
! - Irradiation map around a line [J.s^{-1}.m^{-2}.Hz^{-1}.sr^{-1}]
! - levels population n and nstar
! - Cooling rates PHIij = njRij - niRji
! - Electron density
!
! It uses all the NLTE modules in NLTE/ and it is called in mcfost.f90 similar to
! mol_transfer.f90
!
!
MODULE AtomicTransfer

 use metal, only : Background
 use spectrum_type
 use grid_type, only : atmos
 use readatom
 

 IMPLICIT NONE

 CONTAINS
 
 SUBROUTINE Atomic_transfer()
  integer :: atomunit = 1, nact, nat, integer

 
  ! Initialize the NLTE grid and space 
  !CALL initGrid(Nspace, T, ner, nHtotr, vturbr, Vchar)
  ! if the electron density is not provided by the model or simply want to
  ! recompute it
  if (atmos%calc_ne) CALL SolveElectronDensity(atmos%ne)

  !Read atomic models and allocate space for n, nstar, vbroad, ntotal, Rij, Rji
  ! on the whole grid space.
  CALL readAtomicModels(atomunit)

  NLTEspec%Nact = atmos%Nactiveatoms
  NLTEspec%atmos => atmos
  ! Initialize the wavelength grid depening on the number of transitions PASSIVE/ACTIVE
  CALL make_wavelength_grid(atmos%Natom, atmos%Atoms,NLTEspec%lambda)
  NLTEspec%Nwaves = size(NLTEspec%lambda)
  
  re_init = .false. !first allocation
  !when re_init is .false., table for opacities are allocated for all wavelengths.
  ! when it is .true., opacities are set to zero for next points.
  allocate(NLTEspec%I(NLTEspec%Nwaves))
  allocate(NLTEspec%Flux(NLTEspec%Nwaves))
  CALL initAS(re_init)

  !Compute LTE populations for all atoms, nstar. Open collision file for active atoms
  ! compute nHmin density from neutral hydrogen density (sum over neutral levels)
  Call setLTEcoefficients ()

  !initiate NLTE popuplations for ACTIVE atoms, depending on the choice of the solution
  ! CALL initialSol()
  ! for now initialSol() is replace by this if loop on active atoms
  if (atmos%Nactiveatoms.gt.0) then
    write(*,*) "solving for ", atmos%Nactiveatoms, " active atoms"
    do nact=1,atmos%Nactiveatoms
     write(*,*) "Setting initial solution for active atom ", atmos%ActiveAtoms(nact)%ID, &
      atmos%ActiveAtoms(nact)%active
     atmos%ActiveAtoms(nact)%n = 1d0 * atmos%ActiveAtoms(nact)%nstar
    end do
  end if !end replacing initSol()
  
  ! Main loop
  do icell=1, atmos%Nspace
  
   ! Compute background opacities for PASSIVE bound-bound and bound-free transitions
   ! at all wavelength points
   CALL Background(icell)


  ! Read collisional data and fill collisional matrix C(Nlevel**2) for each ACTIVE atoms.
  ! Initialize at C=0.0 for each cell points.
  ! the matrix is allocated for ACTIVE atoms only in setLTEcoefficients and the file open 
  ! before the transfer starts and closed at the end.
   do nact=1,atmos%Nactiveatoms 
    if (atmos%ActiveAtoms(nact)%active) &
     CALL CollisionRate(icell, atmos%ActiveAtoms(nact)) 
   end do
   
  ! if atmos.Niter > 0
  ! here NLTE loop
  ! integrate Rij, Rji over angle
  ! construct Gamma(Nlevel*Nlevel)
  ! Call StatEquil() ! Solve the SEE
  ! if (recompute_ne_density) CALL SolveElectronDensity(atmos%ne, "NEMODEL")
  ! if (solve_chemical_equilibrium) CALL ChemicalEq()
  ! end for NLTE loop
  ! CALL initAS(re_init=.true.)
  ! has populations of atoms change, background opacities change
  ! CALL Background(icell)
  ! end over Niter 
  
  !Note, this is note the line integrated J, but the J entering in the isotropic
  ! scattering emission.
  do iray=1, n_rayons
    tau -> NLTEspec%chi_c + NLTEspec%chi
    Sny = NLTEspec%eta_c + NLTEspec%eta + NLTEspec%sca%NLTEspecJ + other_eta
    NLTEspec%J = Iincident * tau + (Sny/NLTEspec%chi + NLTEspec%chi_c + other_chi) 
    NLTEspec%J = NLTEspec%J += Icenident *exp(tau) + Sny * (1. - exp(tau)))
  end do
   
  ! set opacities to 0.0 for next cell point.
   CALL initAS(re_init=.true.)
  end do !loop over cells


 ! Transfer ends, save data, free some space and leave
 !write to disk ...
 do nact=1,atmos%Nactiveatoms
  CALL closeCollisionFile(atmos%ActiveAtoms(nact))
 end do
 CALL freeAS() !deallocate opacity arrays
 CALL freeSpectrum()
 CALL freeGrid()  
  

 
 RETURN
 END SUBROUTINE

END MODULE AtomicTransfer