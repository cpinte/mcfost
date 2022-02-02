! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atomic (only) systems using the MALI method
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
module atom_transfer

   use parametres
   use io_atom, only : read_atomic_models, write_pops_atom
   use elecdensity, only : solve_ne, write_electron, read_electron
   use grid, only : lcalc_ne
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops
   use atom_type, only : atoms
   use gas_contopac, only : background_continua_lambda, alloc_rayleigh_xsections

   implicit none
                    
   integer :: omp_chunk_size
   real(kind=dp), allocatable :: Itot(:), Icont(:)

   contains



   subroutine atom_line_transfer()
   ! --------------------------------------------------------------------------- !
   ! This routine initialises the necessary quantities for atomic line transfer
   ! and calls the appropriate routines for LTE or NLTE transfer.
   ! --------------------------------------------------------------------------- !
      integer  :: ne_initial
      logical :: lelectron_read
      real(kind=dp) :: dne, chit(4),etat(4),xlam(4)

      omp_chunk_size = 1!max(nint( 0.01 * n_cells / nb_proc ),1)
      mem_alloc_tot = 0

      !read atomic models
      call read_atomic_Models()

      !is there an electron density file in the current directory ?
      call read_electron(lelectron_read)
      if (lelectron_read) then
         !yes, use that as initial solution
         ne_initial = 1
         lcalc_ne = lsolve_for_ne
      else
         !no, use hydrogen ionisation
         ne_initial = 0
         if (lcalc_ne) then
            write(*,*) "Solving for electron density"
            write(*,*) " -> initial solution from H+M ionisation."
         endif
      endif

      !no electron density in the model nor previous fits file calc_ne == True.
      !if a previous file is found (lelectron_read == True) and lsolve_for_ne, then
      !calc_ne is set to .true. Electron density is re-evaluated using the populations
      !in the fits files (H_Ionisation otherwise).
      if (lcalc_ne) then !eventually computed with the previous non-LTE pops !
         call Solve_ne(ne_initial, .true., dne)
         call write_Electron
      else
         if (lsolve_for_ne) then
            write(*,*) "Forcing calculation of electron density"
            write(*,*) " -> initial solution from model."
            call Solve_ne(ne_initial, .true., dne)
            call write_Electron
         endif
      endif

      call ltepops_atoms() 

      !call init_Spectrum(Nrayone,lam0=lam0,vacuum_to_air=lvacuum_to_air)
      !if (n_etoiles > 0) call init_stellar_disk
      !call alloc_atom_quantities
      !call compute_background_continua


      xlam(1) = 10.0
      xlam(2) = 300.0
      xlam(3) = 500.0
      xlam(4) = 3000.
      call alloc_rayleigh_xsections(4,xlam)
      call background_continua_lambda(1, 4, xlam, chit, etat)
      write(*,*) chit 
      write(*,*) etat

      do ne_initial=1, atoms(1)%p%Ntr
         write(*,*) "trans", atoms(1)%p%tab_trans(ne_initial)
      enddo


      return
   end subroutine atom_line_transfer

end module atom_transfer