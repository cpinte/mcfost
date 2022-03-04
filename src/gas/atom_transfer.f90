! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atomic (only) systems using the MALI method
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
module atom_transfer

   use parametres
   use io_atom, only : read_atomic_models, write_pops_atom
   use wavelengths_gas, only : make_wavelength_nlte, lambda
   use elecdensity, only : solve_ne, write_electron, read_electron
   use grid, only : lcalc_ne, move_to_grid
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops
   use atom_type, only : atoms
   use init_mcfost, only :  nb_proc
   use gas_contopac, only : background_continua_lambda
   use opacity_atom, only : alloc_atom_opac, Itot, psi, opacity_atom_bf_loc, opacity_atom_bb_loc
   use optical_depth, only : integ_ray_atom
   use utils, only : cross_product, gauss_legendre_quadrature

   implicit none
                    
   integer :: omp_chunk_size

   contains



   subroutine atom_line_transfer()
   ! --------------------------------------------------------------------------- !
   ! This routine initialises the necessary quantities for atomic line transfer
   ! and calls the appropriate routines for LTE or NLTE transfer.
   ! --------------------------------------------------------------------------- !
      integer  :: ne_initial
      logical :: lelectron_read
      real(kind=dp) :: dne

      integer :: la, j, icell0
      logical :: lintersect
      real(kind=dp) :: rr, u,v,w,u0,w0,v0,x0,y0,z0,x(3),y(3),uvw(3),cos_theta(1), weight_mu(1)
      real(kind=dp), allocatable :: chit(:),etat(:),xlam(:), Ftot(:,:)

      omp_chunk_size = max(nint( 0.01 * n_cells / nb_proc ),1)
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

      !every indexes are defined on the lambda frequency. 
      !So even when lambda is passed as an argument, it is expected
      !that the indexes correspond at that grid.
      call make_wavelength_nlte
      !allocate quantities in space and for this frequency grid
      call alloc_atom_opac(size(lambda), lambda)

      ! allocate(chit(size(lambda)),etat(size(lambda)))
      ! chit = 0
      ! etat = 0
      ! call background_continua_lambda(1, size(lambda), lambda, chit, etat)
      ! call opacity_atom_bf_loc(1,size(lambda),lambda,chit,etat)
      ! call opacity_atom_bb_loc(1,1,1,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0.0d0,0d0,0.d0,1.d0,0.0d0,0.0d0,.false.,size(lambda),lambda,chit,etat)

      u = 0.0_dp
      v = 0.0_dp
      w = 1.0_dp

      uvw = (/u,v,w/) !vector position
      x = (/1.0_dp,0.0_dp,0.0_dp/)
      y = -cross_product(x, uvw)

      ! Ray tracing : on se propage dans l'autre sens
      u0 = -u ; v0 = -v ; w0 = -w
      call gauss_legendre_quadrature(0.0_dp, 1.0_dp, size(cos_theta), cos_theta, weight_mu)
      allocate(Ftot(size(lambda),nb_proc))
      Ftot = 0.0
      do j=1,size(cos_theta)

         rr = Rmax * sqrt(1.0 - cos_theta(j)**2)

         x0 = 10.0*Rmax*u + rr * y(1)
         y0 = 10.0*Rmax*v + rr * y(2)
         z0 = 10.0*Rmax*w + rr * y(3)

         call move_to_grid(1,x0,y0,z0,u0,v0,w0,icell0,lintersect)
         if (lintersect) then
            call integ_ray_atom(1,icell0,x0,y0,z0,u0,v0,w0,1,.false.,size(lambda),lambda)
         endif
         Ftot(:,1) = Ftot(:,1) + Itot(:,1,1) * weight_mu(j)
      enddo
      open(1,file="sun.txt",status="unknown")
      do la=1,size(lambda)

         write(1,*) lambda(la), Ftot(la,1)
      enddo
      close(1)



      return
   end subroutine atom_line_transfer

end module atom_transfer