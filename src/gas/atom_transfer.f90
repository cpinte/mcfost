! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atomic (only) systems using the MALI method
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
module atom_transfer

   use parametres
   use input, only               : lkeplerian, linfall, limb_darkening, mu_limb_darkening, RT_line_method
   use constantes, only : nm_to_m, m_to_km, km_to_m, au_to_m, deg_to_rad, tiny_real, tiny_dp, pi, deux_pi, pc_to_au, sqrtpi, c_light
   use io_atom, only : read_atomic_models, write_pops_atom
   use wavelengths, only : n_lambda, tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda, n_lambda2, tab_lambda2
   use wavelengths_gas, only : make_wavelengths_nlte, tab_lambda_nm, tab_lambda_cont, n_lambda_cont, &
                                 deallocate_wavelengths_gasrt, make_wavelengths_raytracing, make_wavelengths_flux, Nlambda_max_line
   use elecdensity, only : solve_ne, write_electron, read_electron
   use grid, only : T, vturb,nHtot, nHmin, pos_em_cellule, lcalc_ne, move_to_grid, vfield3d, icompute_atomRT, &
        ne, Voronoi, r_grid, phi_grid, z_grid
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops, LTEpops_atom_loc, LTEpops_H_loc, nH_minus
   use atom_type, only : atoms, atomtype, n_atoms, nactiveatoms, activeAtoms, passiveAtoms, npassiveatoms, &
        hydrogen, helium, adjust_cswitch_atoms, &
                           maxval_cswitch_atoms, lcswitch_enabled, vbroad
   use init_mcfost, only :  nb_proc
   use gas_contopac, only : background_continua_lambda
   use opacity_atom, only : alloc_atom_opac, Itot, psi, dealloc_atom_opac, xcoupling, write_opacity_emissivity_bin, &
        lnon_lte_loop, vlabs, calc_contopac_loc, set_max_damping, deactivate_lines, activate_lines, activate_continua, deactivate_continua
   use see, only : ngpop, Neq_ng, ngpop, alloc_nlte_var, dealloc_nlte_var, frac_limit_pops, &
                  init_rates, update_populations, accumulate_radrates_mali, write_rates, init_radrates_atom
   use optical_depth, only : integ_ray_atom
   use utils, only : cross_product, gauss_legendre_quadrature, progress_bar, rotation_3d, Ng_accelerate, Accelerate, check_ng_pops
   use dust_ray_tracing, only    : RT_n_incl, RT_n_az, init_directions_ray_tracing,tab_u_RT, tab_v_RT, tab_w_RT, &
                                   tab_RT_az,tab_RT_incl
   use stars, only               : intersect_stars, max_Tshock, min_Tshock, max_Thp, min_Thp, max_Facc, min_Facc
   use output, only : allocate_atom_maps, flux_total, write_total_flux, write_atomic_maps
   use mcfost_env, only          : dp, time_tick, time_max
   use molecular_emission, only  : ds
   use messages, only : error, warning
   use voigts, only : Voigt
   use broad, only : line_damping

   use healpix_mod
   !$ use omp_lib
   use  naleat, only : seed, stream, gtype

   implicit none

#include "sprng_f.h"
   integer :: omp_chunk_size
   real(kind=dp), allocatable, dimension(:) :: tab_lambda_Sed

   contains

   subroutine io_write_convergence_maps(lmap, map)
      real(kind=dp), intent(in) :: map(n_cells)
      logical, intent(in) :: lmap(n_cells)
      integer :: unit, status

      write(*,*) " ** writing convergence maps..."
      unit = 100
      open(unit, file="dloc.b",form="unformatted",status='unknown',access="sequential",iostat=status)
      write(unit,iostat=status) lmap
      write(unit,iostat=status) map
      close(unit)

      return
   end subroutine io_write_convergence_maps
!TO DO
!      add Trad, Tion
!      checkpointing
   subroutine nlte_loop_mali()
   ! ------------------------------------------------------------------------------------ !
   ! Solve the set of statistical equilibrium equations (SEE) with the 
   ! Multi-level Accelerated Lambda Iterations method (Rybicki & Hummer 92, Apj 202 209).
   !
   ! By default, the SEE are solved twice.
   !  1) A first step with rays starting at the centre of each cell (healpix phase)
   !  2) Randomly distributed rays for random distribution of points of each cell.
   ! ------------------------------------------------------------------------------------ !
      integer :: etape, etape_start, etape_end, iray
      integer :: n_iter, id, i, alloc_status, n_rayons
      ! integer :: , iray_start, n_rayons_max
      integer :: nact, imax, icell_max, icell_max_2
      integer :: icell, ilevel, nb, nr, unconverged_cells
      integer, parameter :: maxIter = 300!150!, maxIter3 = 10
      !ray-by-ray integration of the SEE
      integer, parameter :: one_ray = 1!, n_rayons_start3 = 100
      logical :: lfixed_Rays, lconverged, lprevious_converged
      real :: rand, rand2, rand3, unconverged_fraction
      real(kind=dp) :: precision, vth
      real(kind=dp), dimension(:), allocatable :: xmu,wmu,xmux,xmuy,diff_loc, dne_loc
      real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, weight
      real(kind=dp) :: diff, diff_old, dne, dN, dN1, dNc
      real(kind=dp), allocatable :: dTM(:), dM(:), Tion_ref(:), Tex_ref(:)
      logical, allocatable :: lcell_converged(:)

      logical :: labs
      logical :: l_iterate, l_iterate_ne
      
      !Ng's acceleration
      integer, parameter :: Ng_Ndelay_init = 1 !minimal number of iterations before starting the cycle.
                                               !min is 1 (so that conv_speed can be evaluated)
      real(kind=dp), parameter :: conv_speed_limit_healpix = 5d-2 !1d-3
      real(kind=dp), parameter :: conv_speed_limit_mc = 1d1!1d-1
      real(kind=dp) :: conv_speed_limit
      integer :: iorder, i0_rest, n_iter_accel, iacc
      integer :: Ng_Ndelay, ng_index
      logical :: lng_turned_on, ng_rest, lconverging, accelerated 
      logical, parameter :: lextrapolate_electron = .false. !extrapolate electrons with populations
      real(kind=dp), dimension(:,:), allocatable :: ngtmp

      !convergence check
      type (AtomType), pointer :: at
      integer(kind=8) :: mem_alloc_local
      real(kind=dp) :: diff_cont, conv_speed, conv_acc
      real(kind=dp), allocatable :: Jnu(:,:)

      !timing and checkpointing
      ! NOTE: cpu time does not take multiprocessing (= nb_proc x the real exec time.)
      !-> overall the non-LTE loop
      real :: time_nlte, time_nlte_loop, time_nlte_cpu
      real :: cpu_time_begin, cpu_time_end, time_nlte_loop_cpu
      integer :: count_start, count_end, itime
      !-> for a single iteration
      integer :: cstart_iter, cend_iter
      real :: time_iteration, cpustart_iter, cpuend_iter, time_iter_avg
      integer :: ibar, n_cells_done, n_cells_remaining

      ! -------------------------------- INITIALIZATION -------------------------------- !
      write(*,*) '-------------------------- NON-LTE LOOP ------------------------------ '
      time_nlte_loop = 0!total time in the non-LTE loop for all steps
      time_nlte_loop_cpu = 0
      !non-LTE mode
      lnon_lte_loop = .true. !for substracting speed of reference cell.
      labs = .true.
      lfixed_rays = .true.
      id = 1

      !Use non-LTE pops for electrons and background opacities.
      do nact=1,NactiveAtoms
         if (.not.activeatoms(nact)%p%nltepops) activeatoms(nact)%p%nltepops = .true.
      enddo

      allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
      allocate(dTM(Nactiveatoms)); dTM=0d0 !keep tracks of Tex for all cells for each atom
      allocate(Tex_ref(Nactiveatoms)); Tex_ref=0d0 !keep tracks of max Tex for all cells for each line of each atom
      allocate(Tion_ref(Nactiveatoms)); Tion_ref=0d0 !keep tracks of max Tion for all cells for each cont of each atom
      diff_old = 1.0_dp
      dM(:) = 1.0_dp
      !init in case l_iterate_ne  is .false. (does not change a thing if l_iterate_ne)
      dne = 0.0_dp

      !-> negligible
      mem_alloc_local = 0
      mem_alloc_local = mem_alloc_local + sizeof(dM)+sizeof(dTm)+sizeof(Tex_ref)+sizeof(Tion_ref)


      !How many steps and which one
      etape_start = istep_start
      etape_end = istep_end
      ! lprecise_pop = .false.
      ! if (lprecise_pop) then
      ! !iray_start reset to 1 we recompute with twice has much rays but from the start
      !    etape_end = 3
      !    n_rayons_max =  n_rayons_start3 * (2**(maxIter3-1))
      ! endif
      if (allocated(stream)) deallocate(stream)
      allocate(stream(nb_proc),stat=alloc_status)
      if (alloc_status > 0) call error("Allocation error stream")
      if (allocated(ds)) deallocate(ds)
      allocate(ds(one_ray,nb_proc))
      allocate(vlabs(one_ray,nb_proc))
      allocate(lcell_converged(n_cells),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error lcell_converged")
      write(*,*) " size lcell_converged:", sizeof(lcell_converged) / 1024./1024./1024.," GB"
      allocate(diff_loc(n_cells), dne_loc(n_cells),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error diff/dne_loc")
      write(*,*) " size diff_loc:", 2*sizeof(diff_loc) / 1024./1024./1024.," GB"
      write(*,*) ""

      !-> negligible
      mem_alloc_local = mem_alloc_local + sizeof(ds) + sizeof(stream)
      ! allocate(jnu(n_lambda,n_cells))
      ! mem_alloc_local = mem_alloc_local + sizeof(jnu)
      call alloc_nlte_var(one_ray,mem=mem_alloc_local)

      ! --------------------------- OUTER LOOP ON STEP --------------------------- !

      step_loop : do etape=etape_start, etape_end

         write(*,*) ""
         !generate rays for a given step
         if (etape==1) then
            !use stepan if healpix_lorder < 2?
            n_rayons = healpix_npix(healpix_lorder)
            allocate(xmux(n_rayons),xmu(n_rayons),xmuy(n_rayons),wmu(n_rayons))
            wmu(:) = healpix_weight(healpix_lorder)
            write(*,'(" ****-> Using "(1I8)" pixels for healpix, resolution of "(1F12.3)" degrees")') n_rayons,  &
                 healpix_angular_resolution(healpix_lorder)
            call healpix_sphere(healpix_lorder,xmu,xmux,xmuy)
            if (etape_end > 1) then
               !use that etape as an initial solution for step 2
               precision = 1d-1
            else
               precision = dpops_max_error
            endif
            conv_speed_limit = conv_speed_limit_healpix
         else if (etape==2) then
            precision = dpops_max_error
            n_rayons = N_rayons_mc
            write(*,'(" ****-> Using "(1I4)" rays for Monte-Carlo step, ~resolution of "(1F12.3)" degrees")') n_rayons, &
                 360.0 * sqrt(pi/real(n_rayons))/pi
            allocate(wmu(n_rayons))
            wmu(:) = 1.0_dp / real(n_rayons,kind=dp)
            conv_speed_limit = conv_speed_limit_mc
         ! else if (etape==3) then
         ! !or solution with fixed rays that increase after each iteration??
         !    lfixed_rays = .false.
         !    n_rayons = n_rayons_start3 !start, same as step 2
         !    conv_speed_limit = conv_speed_limit_mc
         !    precision = min(1d-1,10.0*dpops_max_error)
         !    !only once for all iterations on this step
         !    stream = 0.0
         !    stream(:) = [(init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT),i=1,nb_proc)]
         else
            call error("(n_lte_loop_mali) etape unkown")
         end if
         write(*,*) ""
         call system_clock(count_start,count_rate=time_tick,count_max=time_max)
         call cpu_time(cpu_time_begin)

         !common variables for a given step
         lconverged = .false.
         lprevious_converged = lforce_lte
         lcell_converged(:) = .false.
         diff_loc(:) = 1d50
         dne_loc(:) = 0.0_dp
         lng_turned_on = .false.
         n_iter = 0
         conv_speed = 0.0
         conv_acc = 0.0
         diff_old = 0.0
         time_iter_avg = 0.0
         unconverged_fraction = 0.0
         unconverged_cells = 0

         !***********************************************************!
         ! *************** Main convergence loop ********************!
         do while (.not.lconverged)
            ! evaluate the time for an iteration independently of nb_proc
            ! time of a single iteration for this step
            call system_clock(cstart_iter,count_rate=time_tick,count_max=time_max)
            call cpu_time(cpustart_iter)

            if (lcswitch_enabled) then
               !deactivate
               if (maxval_cswitch_atoms()==1.0_dp) then
                  write(*,*) " ** cswitch off."
                  lcswitch_enabled = .false.
               endif
            endif

            n_iter = n_iter + 1
            !-> ng_index depends only on the value of n_iter for a given Neq_ng
            ! it is the index of solutions stored in ngpop. ng_index is 1 for the current solution.
            ! at which point we can accelerate. Therefore, ng_index if far from 1 at init.
            ng_index = Neq_ng - mod(n_iter-1,Neq_ng) !overwritten if lng_acceleration.
            ! NOTE :: the index 1 is hard-coded within the non-LTE loop as it stores the actual (running) value
            !           of the populations and electronic density.
            ! 
            !                    goes with maxIter
            write(*,'(" *** Iteration #"(1I4)"; step #"(1I1)"; threshold: "(1ES11.2E3)"; Nrays: "(1I5))') &
                     n_iter, etape, precision, n_rayons
            ! write(*,'("  -- min(diff_loc)="(1ES13.5E3))') minval(diff_loc,mask=(icompute_atomRT>0).and.(diff_loc >= 1d-2 * precision))
            ! if (n_iterate_ne > 0) &
            !    write(*,'("  -- min(dne_loc)="(1ES13.5E3))') minval(dne_loc,mask=(icompute_atomRT>0).and.(dne_loc >= 1d-2 * precision))
            ibar = 0
            n_cells_done = 0

            if (lfixed_rays) then
               stream = 0.0
               stream(:) = [(init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT),i=1,nb_proc)]
            ! else
            !    !update rays weight
            !    if (allocated(wmu)) deallocate(wmu)
            !    allocate(wmu(n_rayons));wmu(:) = 1.0_dp / real(n_rayons,kind=dp)
            end if

            !init here, to be able to stop/start electronic density iterations within MALI iterations
            l_iterate_ne = .false.
            if( n_iterate_ne > 0 ) then
               l_iterate_ne = ( mod(n_iter,n_iterate_ne)==0 ) .and. (n_iter>Ndelay_iterate_ne)
               ! if (lforce_lte) l_iterate_ne = .false.
            endif
            ! Jnu = 0.0_dp
            call progress_bar(0)
            !$omp parallel &
            !$omp default(none) &
            !$omp private(id,icell,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02,argmt)&
            !$omp private(l_iterate,weight,diff)&
            !$omp private(nact, at) & ! Acceleration of convergence
            !$omp shared(ne,ngpop,ng_index,Ng_Norder, accelerated, lng_turned_on, Jnu) & ! Ng's Acceleration of convergence
            !$omp shared(etape,lforce_lte,n_cells,voronoi,r_grid,z_grid,phi_grid,n_rayons,xmu,wmu,xmux,xmuy,n_cells_remaining) &
            !$omp shared(pos_em_cellule,labs,n_lambda,tab_lambda_nm, icompute_atomRT,lcell_converged,diff_loc,seed,nb_proc,gtype) &
            !$omp shared(stream,n_rayons_mc,lvoronoi,ibar,n_cells_done,l_iterate_ne,Itot,omp_chunk_size,precision,lcswitch_enabled)
            !$omp do schedule(static,omp_chunk_size)
            do icell=1, n_cells
               !$ id = omp_get_thread_num() + 1
               l_iterate = (icompute_atomRT(icell)>0)
               ! stream(id) = init_sprng(gtype, id-1,nb_proc,seed,SPRNG_DEFAULT)
               ! if( (diff_loc(icell) < 1d-2 * precision).and..not.lcswitch_enabled ) cycle

               if (l_iterate) then

                  !Init upward radiative rates to 0 and downward radiative rates to Aji or "Aji_cont"
                  !Init collisional rates
                  call init_rates(id,icell)

                  if (etape==1) then

                     !cell centre
                     if (lvoronoi) then
                        x0 = Voronoi(icell)%xyz(1)
                   	   y0 = Voronoi(icell)%xyz(2)
                   	   z0 = Voronoi(icell)%xyz(3)
                     else
                   	   x0 = r_grid(icell)*cos(phi_grid(icell))
                   	   y0 = r_grid(icell)*sin(phi_grid(icell))
                   	   z0 = z_grid(icell)
                     endif


                     do iray=1, n_rayons

                        w0 = xmu(iray)
                        u0 = xmux(iray)
                        v0 = xmuy(iray)

                        weight = wmu(iray)

                        !for one ray
                        if (.not.lforce_lte) then
                           !-> cannot compute radiative rates here if lforce_lte
                           call integ_ray_atom(id,icell,x0,y0,z0,u0,v0,w0,1,labs,n_lambda,tab_lambda_nm)
                           call xcoupling(id, icell,1)
                           call accumulate_radrates_mali(id, icell,1, weight)
                           ! Jnu(:,icell) = Jnu(:,icell) + weight * Itot(:,1,id)
                        endif
                     enddo


                  else
                   ! Position aleatoire dans la cellule
                     do iray=1,n_rayons

                        rand  = sprng(stream(id))
                        rand2 = sprng(stream(id))
                        rand3 = sprng(stream(id))

                        call pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
                        rand = sprng(stream(id))
                        W0 = 2.0_dp * rand - 1.0_dp !nz
                        W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
                        SRW02 = sqrt(W02)
                        rand = sprng(stream(id))
                        ARGMT = PI * (2.0_dp * rand - 1.0_dp)
                        U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
                        V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)

                        weight = wmu(iray)

                        !for one ray
                        if (.not.lforce_lte) then
                           call integ_ray_atom(id,icell,x0,y0,z0,u0,v0,w0,1,labs,n_lambda,tab_lambda_nm)
                           call xcoupling(id,icell,1)
                           call accumulate_radrates_mali(id, icell,1, weight)
                           ! Jnu(:,icell) = Jnu(:,icell) + weight * Itot(:,1,id)
                        endif
                     enddo !iray

                  end if !etape

                  !update populations of all atoms including electrons
                  call update_populations(id, icell, (l_iterate_ne.and.icompute_atomRT(icell)==1), diff)

                  ! ************************** NG's **************************!
                  ! accelerate locally, cell-by-cell for all levels only
                  ! ************************** END! **************************!

               end if !icompute_atomRT

               ! Progress bar
               !$omp atomic
               n_cells_done = n_cells_done + 1
               n_cells_remaining = size(pack(diff_loc, &
                                    mask=(diff_loc < 1d-2 * precision)))
               if (real(n_cells_done) > 0.02*ibar*n_cells) then
             	   call progress_bar(ibar)
             	   !$omp atomic
             	   ibar = ibar+1
               endif

            end do !icell
            !$omp end do
            !$omp end parallel
            call progress_bar(50)
            write(*,*) " " !for progress bar

            !***********************************************************!
            ! **********  Ng's acceleration administration *************!
            if (lng_acceleration) then
          	!be sure we are converging before extrapolating
               if (dpops_max_error > 1d-2) then
                  lconverging = (conv_speed < 0) .and. (-conv_speed < conv_speed_limit)
               else
                  lconverging = (diff_old < 5d-2)!; Ng_Nperiod = 0
               endif
               !or if the number of iterations is too large
               lconverging = lconverging .or. (n_iter > int(real(maxIter)/3.0)) !futur deprec of this one (non-local op)
          	   if ( (n_iter>max(Ng_Ndelay_init,1)).and.lconverging ) then
          		   if (.not.lng_turned_on) then
                     write(*,*) " +++ Activating Ng's acceleration +++ "
          			   lng_turned_on = .true.
          			   Ng_Ndelay = n_iter; iacc = 0
                     n_iter_accel = 0; i0_rest = 0; ng_rest = .false.
          		   endif
            	endif
               if (lng_turned_on) then
                  if (unconverged_fraction < 15.0) lng_turned_on = .false.
               endif
            endif
            !***********************************************************!
            ! ********************** GLOBAL NG's ***********************!
            ! Here minimize Ncells x Nlevels x Ng_Norder per atoms.
            accelerated = .false.!(maxval_cswitch_atoms()==1.0_dp)
            if ( (lNg_acceleration .and. lng_turned_on).and.(.not.lcswitch_enabled)&
                  .and.(.not.lprevious_converged) ) then
               iorder = n_iter - Ng_Ndelay
               if (ng_rest) then
                  write(*,'(" -> Acceleration relaxes... "(1I2)" / "(1I2))') iorder-i0_rest, Ng_Nperiod
                  if (iorder-i0_rest == Ng_Nperiod) ng_rest = .false.
               else
                  i0_rest = iorder
                  iacc = iacc + 1
                  ng_index = Neq_Ng - mod(iacc-1,Neq_ng)
                  write(*,'(" -> Accumulate solutions... "(1I2)" / "(1I2))') iacc, Neq_ng
                  !index 1 is the running one in the main loop for new values.
                  ngpop(:,:,:,ng_index) = ngpop(:,:,:,1)
                  accelerated = (iacc==Neq_ng)

                  if (accelerated) then
                     do nact=1,NactiveAtoms
                        at => ActiveAtoms(nact)%p
                        allocate(ngtmp(N_cells*at%Nlevel,Neq_ng))
                        !Flatten the array such that for each cell there are all levels
                        ngtmp(:,:) = reshape(ngpop(1:at%Nlevel,nact,:,:),(/n_cells*at%Nlevel,Neq_ng/))
                        !Flatten the array such that for each level there are all cells.
                        ! ngtmp(:,:) = reshape(& ! first transpose to (N_cells,Nlevel,Neq_ng)
                        !    reshape(ngpop(1:at%Nlevel,nact,:,:),shape=[n_cells,at%Nlevel,Neq_ng],order=[2,1,3]),&!then flatten
                        !    (/n_cells*at%Nlevel,Neq_ng/))

                        call Accelerate(n_cells*at%Nlevel,Ng_Norder,ngtmp) 

                        !reform in (Nlevel,N_cells,Neq_ng)
                        ngpop(1:at%Nlevel,nact,:,:) = reshape(ngtmp,(/at%Nlevel,n_cells,Neq_ng/))
                        ! ngpop(1:at%Nlevel,nact,:,:) = reshape(&!first reform
                        !    reshape(ngtmp(:,:),shape=[n_cells,at%Nlevel,Neq_Ng]),&!then tranpose back
                        !    shape=[at%Nlevel,n_cells,Neq_ng],order=[2,1,3])
                        deallocate(ngtmp)
                        !check negative populations and mass conservation for that atom
                        call check_ng_pops(at%Nlevel,n_cells,Neq_ng,ngpop(1:at%Nlevel,nact,:,:),at%Abund*nHtot(:))
                        at => NULL()
                     enddo
                     ! Accelerate electrons ? but still needs rest of SEE+ne loop.
                     if ( (lextrapolate_electron).and.(l_iterate_ne) ) then
                        write(*,*) " -- extrapolating electrons..."
                        call Accelerate(n_cells,Ng_Norder,ngpop(1,NactiveAtoms+1,:,:))
                     endif
                     n_iter_accel = n_iter_accel + 1
                     ng_rest = (Ng_Nperiod > 0); iacc = 0
                  endif
               endif
            endif
            !***********************************************************!

            !***********************************************************!
            ! *******************  update ne metals  *******************!
            if (l_iterate_ne) then
               id = 1
               dne = 0.0_dp !=max(dne_loc)
               dne_loc(:) = abs(1.0_dp - ne(:)/(1d-100 + ngpop(1,NactiveAtoms+1,:,1)))
               ! write(*,'("  OLD ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') &
               !    minval(ne,mask=(icompute_atomRT>0)), maxval(ne)
               !$omp parallel &
               !$omp default(none) &
               !$omp private(id,icell,l_iterate,nact,ilevel,vth,nb,nr)&
               !$omp shared(T,vturb,nHmin,icompute_atomRT,ne,n_cells,ngpop,ng_index,lng_acceleration)&
               !$omp shared(tab_lambda_nm,atoms,n_atoms,dne,PassiveAtoms,NactiveAtoms,NpassiveAtoms,Voigt)
               !$omp do schedule(dynamic,1)
               do icell=1,n_cells
                  !$ id = omp_get_thread_num() + 1
                  l_iterate = (icompute_atomRT(icell)==1)
                  if (l_iterate) then
                     ! dne = max(dne, abs(1.0_dp - ne(icell)/ne_new(icell)))
                     dne = max(dne, abs(1.0_dp - ne(icell)/ngpop(1,NactiveAtoms+1,icell,1)))
                     ! -> needed for lte pops !
                     ne(icell) = ngpop(1,NactiveAtoms+1,icell,1) !ne(icell) = ne_new(icell)
                     if (.not.lng_acceleration) &
                        ngpop(1,NactiveAtoms+1,icell,ng_index) = ne(icell)
                     call LTEpops_H_loc(icell)
                     nHmin(icell) = nH_minus(icell)
                     do nact = 2, n_atoms
                        call LTEpops_atom_loc(icell,Atoms(nact)%p,.false.)
                     enddo
                     !update profiles only for passive atoms. For active atoms we need new non-LTE pops.
                     !If LTE pops are not updated (so no ne_new) profiles of LTE elements are unchanged.
                     do nact = 1, NpassiveAtoms
                        do ilevel=1,PassiveAtoms(nact)%p%nline
                           if (.not.PassiveAtoms(nact)%p%lines(ilevel)%lcontrib) cycle
                           if (PassiveAtoms(nact)%p%lines(ilevel)%Voigt) then
                              nb = PassiveAtoms(nact)%p%lines(ilevel)%nb; nr = PassiveAtoms(nact)%p%lines(ilevel)%nr
                              PassiveAtoms(nact)%p%lines(ilevel)%a(icell) = line_damping(icell,PassiveAtoms(nact)%p%lines(ilevel))
                              !tmp because of vbroad!
                              vth = vbroad(T(icell),PassiveAtoms(nact)%p%weight, vturb(icell))
                              !-> beware, temporary array here. because line%v(:) is fixed! only vth changes
                              PassiveAtoms(nact)%p%lines(ilevel)%phi(:,icell) = &
                                   Voigt(PassiveAtoms(nact)%p%lines(ilevel)%Nlambda, PassiveAtoms(nact)%p%lines(ilevel)%a(icell), &
                                   PassiveAtoms(nact)%p%lines(ilevel)%v(:)/vth)&
                                    / (vth * sqrtpi)
                           endif
                        enddo
                     enddo !over passive atoms
                  end if !if l_iterate
               end do
               !$omp end do
               !$omp end parallel
               ! write(*,'("  NEW ne(min)="(1ES16.8E3)" m^-3 ;ne(max)="(1ES16.8E3)" m^-3")') &
               !    minval(ne,mask=(icompute_atomRT>0)), maxval(ne)
               ! write(*,*) ''
               ! if ((dne < 1d-2 * precision).and.(.not.lcswitch_enabled)) then
               !    !Or compare with 3 previous values of dne ? that should be below 1e-2 precision
               !    !Do we need to restart it eventually ?
               !    write(*,*) " *** stopping electronic density convergence at iteration ", n_iter
               !    n_iterate_ne = 0
               ! endif
            end if
            !***********************************************************!

            !***********************************************************!
            ! ****************  Convergence checking  ******************!
            id = 1
            dM(:) = 0.0 !all pops
            diff_cont = 0.0
            diff = 0.0
            !$omp parallel &
            !$omp default(none) &
            !$omp private(id,icell,l_iterate,dN1,dN,dNc,ilevel,nact,at,nb,nr,vth)&
            !$omp shared(ngpop,Neq_ng,ng_index,Activeatoms,lcell_converged,vturb,T,lng_acceleration,diff_loc)&!diff,diff_cont,dM)&
            !$omp shared(icompute_atomRT,n_cells,precision,NactiveAtoms,nhtot,llimit_mem,tab_lambda_nm,voigt,dne_loc) &
            !$omp reduction(max:dM,diff,diff_cont)
            !$omp do schedule(dynamic,1)
            cell_loop2 : do icell=1,n_cells
               !$ id = omp_get_thread_num() + 1
               l_iterate = (icompute_atomRT(icell)>0)

               if (l_iterate) then

                  !Local only
                  dN = 0.0 !for all levels of all atoms of this cell
                  dNc = 0.0!cont levels
                  do nact=1,NactiveAtoms
                     at => ActiveAtoms(nact)%p

                     do ilevel=1,at%Nlevel
                        if ( ngpop(ilevel,nact,icell,1) >= frac_limit_pops * at%Abund*nHtot(icell) ) then
                           dN1 = abs(1d0-at%n(ilevel,icell)/ngpop(ilevel,nact,icell,1))
                           dN = max(dN1, dN)
                           dM(nact) = max(dM(nact), dN1)
                        endif
                     end do !over ilevel
                     dNc = max(dNc, dN1)!the last dN1 is necessarily the last ion of each species

                     at => NULL()
                  end do !over atoms

                  !compare for all atoms and all cells
                  diff = max(diff, dN) ! pops
                  diff_cont = max(diff_cont,dNc)

                  lcell_converged(icell) = (dN < precision).and.(dne_loc(icell) < precision) !(diff < precision)
                  diff_loc(icell) = max(dN, dne_loc(icell))

                  !Re init for next iteration if any
                  do nact=1, NactiveAtoms
                     at => ActiveAtoms(nact)%p
                     !first update the populations with the new value
                     at%n(:,icell) = ngpop(1:at%Nlevel,nact,icell,1)
                     !Then, store Neq_ng previous iterations. index 1 is for the current one.
                     !if Ng_neq = 3 :
                     !  -> the oldest solutions stored is in ng_index = 3 - mod(1-1,3) = 3
                     !  -> the newest (current) one is in ng_index = 3 - mod(3-1,3) = 1
                     ! cannot accumulate here if global Ng's
                     ! !! if condition breaks if local Ng's !!
                     if (.not.lng_acceleration) &
                        ngpop(1:at%Nlevel,nact,icell,ng_index) = at%n(:,icell)

                     !Recompute damping and profiles once with have set the new non-LTE pops (and new ne) for next ieration.
                     !Only for Active Atoms here. PAssive Atoms are updated only if electronic density is iterated.
                     !Also need to change if profile interp is used or not! (a and phi)
                     do ilevel=1,at%nline
                        if (.not.at%lines(ilevel)%lcontrib) cycle
                        if (at%lines(ilevel)%Voigt) then
                           nb = at%lines(ilevel)%nb; nr = at%lines(ilevel)%nr
                           at%lines(ilevel)%a(icell) = line_damping(icell,at%lines(ilevel))
                           !tmp because of vbroad!
                           vth = vbroad(T(icell),at%weight, vturb(icell))
                           !-> beware, temporary array here. because line%v(:) is fixed! only vth changes
                           at%lines(ilevel)%phi(:,icell) = &
                                Voigt(at%lines(ilevel)%Nlambda, at%lines(ilevel)%a(icell), at%lines(ilevel)%v(:)/vth) &
                                / (vth * sqrtpi)
                        endif
                     enddo
                  end do
                  at => null()

                  !if electronic density is not updated, it is not necessary
                  !to compute the lte continous opacities.
                  !but the overhead should be negligible at this point.
                  if (.not.llimit_mem) call calc_contopac_loc(icell)

               end if !if l_iterate
            end do cell_loop2 !icell
            !$omp end do
            !$omp end parallel

            if (n_iter > 1) then
               conv_speed = (diff - diff_old) !<0 converging.
               conv_acc = conv_speed - conv_acc !>0 accelerating
            endif

            if (accelerated) then
               ! write(*,'("            "(1A13))') "(Accelerated)"
               write(*,'("            ("(1A11)" #"(1I4)")")') "Accelerated", n_iter_accel
            endif
            do nact=1,NactiveAtoms
               write(*,'("                  "(1A2))') ActiveAtoms(nact)%p%ID
               ! write(*,'("             Atom "(1A2))') ActiveAtoms(nact)%p%ID
               write(*,'("   >>> dpop="(1ES13.5E3))') dM(nact)
            enddo
            if (l_iterate_ne) then
               write(*,*) ""
               write(*,'("                  "(1A2))') "ne"
               write(*,'("   >>> dne="(1ES13.5E3))') dne
            endif
            write(*,*) ""
            write(*,'(" <<->> diff="(1ES13.5E3)," old="(1ES13.5E3))') diff, diff_old
            write(*,'("   ->> diffcont="(1ES13.5E3))') diff_cont
            write(*,'("   ->> speed="(1ES12.4E3)"; acc="(1ES12.4E3))') conv_speed, conv_acc
            unconverged_cells = size(pack(lcell_converged,mask=(.not.lcell_converged).and.(icompute_atomRT>0)))
            unconverged_fraction = 100.*real(unconverged_cells) / real(size(pack(icompute_atomRT,mask=icompute_atomRT>0)))
            write(*,"('   ->> Unconverged cells #'(1I6)'; fraction:'(1F6.2)'%')") unconverged_cells, unconverged_fraction
            write(*,*) "      -------------------------------------------------- "
            diff_old = diff
            conv_acc = conv_speed

            !force convergence if there are only few unconverged cells remaining
            if ((n_iter > min(maxIter/2,50)).and.(unconverged_fraction < 5.0)) then
               write(*,'("WARNING: there are less than "(1F6.2)" % of unconverged cells after "(1I4)" iterations")') &
                  unconverged_fraction, n_iter
               write(*,*) " -> forcing convergence"
               lconverged = .true.
            endif
            !
            if ((diff < precision).and.(.not.lcswitch_enabled))then!maxval_cswitch_atoms()==1.0_dp
            ! if ( (unconverged_fraction < 3.0).and.maxval_cswitch_atoms()==1.0_dp)then
            ! if ( ((unconverged_fraction < 3.0).and.maxval_cswitch_atoms()==1.0_dp).or.&
               !  ((diff < precision).and.maxval_cswitch_atoms()==1.0_dp) )then
               if (lprevious_converged) then
                  lconverged = .true.
               else
                  lprevious_converged = .true.
               endif
            else
               lprevious_converged = .false.
               if ((lcswitch_enabled).and.(maxval_cswitch_atoms() > 1.0_dp)) then
                  call adjust_cswitch_atoms()
               endif
               if (lfixed_rays) then
                  if (n_iter > maxIter) then
                     call warning("not enough iterations to converge !!")
                     lconverged = .true.
                  end if
               ! else
               !    !increase number of rays only if it converges already ?
               !    n_rayons = n_rayons * 2
               !       ! On continue en calculant 2 fois plus de rayons
               !       ! On les ajoute a l'ensemble de ceux calcules precedemment
               !       ! iray_start = iray_start + n_rayons
               !    if (n_iter >= maxIter3) then
               !       call warning("not enough rays to converge in step 3!!")
               !       lconverged = .true.
               !    endif
              endif
            end if
            !***********************************************************!

            !***********************************************************!
            ! ********** timing and checkpointing **********************!
            ! count time from the start of the step, at each iteration
            call system_clock(count_end,count_rate=time_tick,count_max=time_max)
            call cpu_time(cpu_time_end)
            if (count_end < count_start) then
               time_nlte=real(count_end + (1.0 * time_max)- count_start)/real(time_tick)
            else
               time_nlte=real(count_end - count_start)/real(time_tick)
            endif
            time_nlte_cpu = cpu_time_end - cpu_time_begin

            ! evaluate the time for an iteration independently of nb_proc
            ! time of a single iteration for this step
            call system_clock(cend_iter,count_rate=time_tick,count_max=time_max)
            call cpu_time(cpuend_iter)
            time_iteration = real(cend_iter - cstart_iter)/real(time_tick)
            write(*,'("  --> time iteration="(1F12.4)" min (cpu : "(1F8.4)")")') &
                  mod(time_iteration/60.0,60.0), mod((cpuend_iter-cpustart_iter)/60.0,60.0)
            time_iter_avg = time_iter_avg + time_iteration
            ! -> will be averaged with the number of iterations done for this step


            if (lsafe_stop) then
               if ((time_nlte + time_iteration >=  safe_stop_time)) then
                  lconverged = .true.
                  lprevious_converged = .true.
                  call warning("Time limit would be exceeded, leaving...")
                  write(*,*) " time limit:", mod(safe_stop_time/60.,60.) ," min"
                  write(*,*) " time etape:", mod(time_iter_avg/60.,60.), ' <time iter>=', &
                                                mod(time_iter_avg/real(n_iter)/60.,60.)," min"
                  write(*,*) " time etape (cpu):", mod(time_iter_avg * nb_proc/60.,60.), " min"
                  write(*,*) ' time =',mod(time_nlte/60.,60.), " min"
                  lexit_after_nonlte_loop = .true.
                  exit step_loop
               endif
            endif
            !***********************************************************!

         !***********************************************************!
         !***********************************************************!

         end do !while
         !combine all steps
         time_nlte_loop = time_nlte_loop + time_nlte
         time_nlte_loop_cpu = time_nlte_loop_cpu + time_nlte_cpu
         !-> for this step
         time_iter_avg = time_iter_avg / real(n_iter)

         write(*,'("<time iteration>="(1F8.4)" min; <time step>="(1F8.4)" min")') mod(time_iter_avg/60.,60.), &
                 mod(n_iter * time_iter_avg/60.,60.)
         write(*,'(" --> ~time step (cpu)="(1F8.4)" min")') mod(n_iter * time_iter_avg * nb_proc/60.,60.)

         if (allocated(xmu)) deallocate(xmux,xmuy,xmu)
         if (allocated(wmu)) deallocate(wmu)

         ! write(*,'("  <step #"(1I1)" done! :: threshold="(1ES17.8E3)">")') etape,  precision
         write(*,'("  <step #"(1I1)" done!>")') etape

      end do step_loop

      ! -------------------------------- CLEANING ------------------------------------------ !

      if (n_iterate_ne > 0) then
         write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
         write(*,'("nH/ne "(1ES13.5E3, 1ES13.5E3))') maxval(nHtot/ne,mask=nhtot*ne>0), minval(nHtot/ne,mask=nHtot*ne>0)
         call write_electron
      endif

      do nact=1,N_atoms
         call write_pops_atom(Atoms(nact)%p)
      end do


      if (loutput_rates) call write_rates()

      ! call io_write_convergence_maps(lcell_converged, diff_loc)
      ! open(100, file="jnu.b",form="unformatted",status='unknown',access="sequential")
      ! write(100) tab_lambda_nm
      ! write(100) Jnu
      ! deallocate(Jnu); close(150)
      call dealloc_nlte_var()
      deallocate(dM, dTM, Tex_ref, Tion_ref)
      deallocate(diff_loc, dne_loc)
      deallocate(stream, ds, vlabs, lcell_converged)

      ! --------------------------------    END    ------------------------------------------ !
      lnon_lte_loop = .false.
      if (mod(time_nlte_loop/60./60.,60.) > 1.0_dp) then
         write(*,*) ' (non-LTE loop) time =',mod(time_nlte_loop/60./60.,60.), " h"
      else
         write(*,*) ' (non-LTE loop) time =',mod(time_nlte_loop/60.,60.), " min"
      endif
      write(*,*) '-------------------------- ------------ ------------------------------ '

      return
   end subroutine nlte_loop_mali

   subroutine solve_for_nlte_pops()
      ! Wrapper to solve for the non-LTE populations.
      ! First solve the continuum radiative transfer, 
      ! then for the continuum + gas lines.
      ! In the continuum radiative transfer, bound-bound
      ! transitions are in LTE (only cij and cji).

      ! TO DO allocate cont grid only first to speed up

      ! add sobolev/escape here ?

      call deactivate_lines()
      call alloc_atom_opac(n_lambda, tab_lambda_nm)
      call nlte_loop_mali()
      call activate_lines()
      call dealloc_atom_opac()
      call alloc_atom_opac(n_lambda, tab_lambda_nm)
      Ndelay_iterate_ne = max(Ndelay_iterate_ne,3)
      call nlte_loop_mali()
      return
   end subroutine solve_for_nlte_pops

  subroutine compute_max_relative_velocity(dv)
  !
  ! TO DO: building, add maximum velocity difference between
  !   all cells in all directions.
  !   Compute mean velocity gradient for Sobolev
      use molecular_emission, only : v_proj
      ! use grid, only : cross_cell, test_exit_grid
      real(kind=dp), intent(out) :: dv
      real(kind=dp) :: v1,v2
      integer :: i1, i2, id

      ! integer, parameter :: n_rayons = 10000
      ! integer :: i, next_cell, previous_cell
      ! integer :: i_star, icell_star
      ! logical :: lintersect_stars, lcellule_non_vide
      ! real :: rand, rand2, rand3
      ! real(kind=dp) :: u, v, w, x0, y0, z0, x, y, z, x1, y1, z1
      ! real(kind=dp) :: w02, srw02, argmt,l, l_contrib, l_void_before

      ! if (lvoronoi) then
      !    dv = sqrt( maxval(Voronoi(:)%vxyz(1)**2+Voronoi(:)%vxyz(2)**2+Voronoi(:)%vxyz(3)**2) )
      ! else
      !    dv = sqrt( maxval(sum(vfield3d**2,dim=2)) )
      ! endif

      ! if (dv==0.0_dp) return

      ! if (allocated(stream)) deallocate(stream)
      ! allocate(stream(nb_proc))
      ! stream = 0.0
      ! do i=1,nb_proc
      !    stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
      ! end do

      ! id = 1
      ! dv = 0.0_dp
      ! !$omp parallel &
      ! !$omp default(none) &
      ! !$omp private(id,i1,i2,v1,v2,rand,rand2,rand3,u,v,w,x0,y0,z0,argmt,w02,srw02,x,y,z,i_star,icell_star)&
      ! !$omp private(i,x1,y1,z1,l,l_contrib,l_void_before,next_cell,previous_cell,lcellule_non_vide,lintersect_stars)&
      ! !$omp shared(dv, icompute_atomRT,vfield3d, n_Cells,stream,cross_cell,pos_em_cellule, test_exit_grid)
      ! !$omp do schedule(static,1)
      ! do i1=1,n_cells
      !    !$ id = omp_get_thread_num() + 1
      !    if (icompute_atomRT(i1) <= 0) cycle
      !    !compute velocity gradient between all cells in the ray path from i1
      !    do i=1, n_rayons
      !       rand  = sprng(stream(id))
      !       rand2 = sprng(stream(id))
      !       rand3 = sprng(stream(id))
      !       call pos_em_cellule(i1,rand,rand2,rand3,x,y,z)
      !       rand = sprng(stream(id))
      !       W = 2.0_dp * rand - 1.0_dp
      !       W02 =  1.0_dp - W*W
      !       SRW02 = sqrt(W02)
      !       rand = sprng(stream(id))
      !       ARGMT = PI * (2.0_dp * rand - 1.0_dp)
      !       U = SRW02 * cos(ARGMT)
      !       V = SRW02 * sin(ARGMT)

      !       call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)
      !       x1=x;y1=y;z1=z
      !       x0=x;y0=y;z0=z
      !       next_cell = i1
      !       v1 = v_proj(i1,x0,y0,z0,u,v,w)
      !       infinie : do ! Boucle infinie
      !          i2 = next_cell
      !          x0=x1 ; y0=y1 ; z0=z1
      !          lcellule_non_vide = (i2 <= n_cells)
      !          if (test_exit_grid(i2, x0, y0, z0)) exit infinie
      !          if (lintersect_stars) then
      !             if (i2 == icell_star) exit infinie
      !          end if
      !          if (i2 <= n_cells) then
      !             lcellule_non_vide = (icompute_atomRT(i2) > 0)
      !             if (icompute_atomRT(i2) < 0) exit infinie
      !          endif

      !          call cross_cell(x0,y0,z0, u,v,w, i2, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)

      !          !compute velocity gradient between cell i1 and crossed cells i2
      !          if (lcellule_non_vide) then
      !             v2 = v_proj(i2,x0,y0,z0,u,v,w)
      !             dv = max(dv,abs(v2-v1))
      !          endif

      !       enddo infinie
      !    enddo
      ! enddo
      ! !$omp end do
      ! !$omp end parallel
      ! deallocate(stream)
      ! write(*,'("maximum gradv="(1F12.3)" km/s")') dv*1d-3

      ! return

      if (lvoronoi) then
         dv = sqrt( maxval(Voronoi(:)%vxyz(1)**2+Voronoi(:)%vxyz(2)**2+Voronoi(:)%vxyz(3)**2,mask=icompute_atomRT>0) )
         if (dv == 0.0_dp) return

         ! dv = 0.0_dp
         ! !$omp parallel &
         ! !$omp default(none) &
         ! !$omp private(id,i1, i2,v1, v2)&
         ! !$omp shared(dv, icompute_atomRT, Voronoi, n_Cells)
         ! !$omp do schedule(dynamic,1)
         ! do i1=1,n_cells
         !    !$ id = omp_get_thread_num() + 1
         !    v1 = sqrt(sum(Voronoi(i1)%vxyz**2))
         !    do i2=1,n_cells
         !       v2 = sqrt(sum(Voronoi(i2)%vxyz**2))
         !       if ((icompute_atomRT(i1)>0).and.(icompute_atomRT(i2)>0)) then
         !          dv = max(dv,abs(v1-v2))
         !       endif
         !    enddo
         ! enddo
         ! !$omp end do
         ! !$omp end parallel
      else
         dv = sqrt( maxval(sum(vfield3d**2,dim=2),mask=icompute_atomRT>0) )
         if (dv==0.0_dp) return

      !    dv = 0.0_dp
      !    !$omp parallel &
      !    !$omp default(none) &
      !    !$omp private(id,i1,i2,v1,v2)&
      !    !$omp shared(dv, icompute_atomRT,vfield3d, n_Cells)
      !    !$omp do schedule(dynamic,1)
      !    do i1=1,n_cells
      !       !$ id = omp_get_thread_num() + 1
      !       v1 = sqrt(sum(vfield3d(i1,:)**2))
      !       do i2=1,n_cells
      !          v2 = sqrt(sum(vfield3d(i2,:)**2))
      !          if ((icompute_atomRT(i1)>0).and.(icompute_atomRT(i2)>0)) then
      !             dv = max(dv,abs(v1-v2))
      !          endif
      !       enddo
      !    enddo
      !    !$omp end do
      !    !$omp end parallel
      endif

      write(*,'("maximum gradv="(1F12.3)" km/s")') dv*1d-3

      return
   end subroutine compute_max_relative_velocity

   subroutine setup_image_grid()
   !to do lmono -limg
   !keep somewhere tab_lambda_sed = tab_lambda because tab_lambda is overwritten in non-LTE
      integer :: kr, nat

      call deallocate_wavelengths_gasrt(tab_lambda)
      call dealloc_atom_opac()
      call init_directions_ray_tracing()

      !defined in alloc_Atomic_maps but I need n_lambda so RT_line_method
      !before
      if (npix_x_save > 1) then
         RT_line_method = 2
         if (npix_y_save == 1) RT_line_method = 3
      else
         RT_line_method = 1
      endif

      if (RT_line_method==1) then
         if (lsed) then
            !tab lambda from file
            call make_wavelengths_flux(tab_lambda_sed,.true.)
         else
            !tab lambda similar to the non-LTE grid with lines edge from parameter file.
            call make_wavelengths_flux(tab_lambda_nm,.false.)
         endif
      else
         !create a wavelength grid around ray-traced lines
         !keep all transitions overlapping with those lines;
         call make_wavelengths_raytracing(tab_lambda_nm)
      endif
      n_lambda = size(tab_lambda_nm)
      tab_lambda = tab_lambda_nm * nm_to_m!micron

      call alloc_atom_opac(n_lambda, tab_lambda_nm)
      call allocate_atom_maps()
      if (laccretion_shock) then
         max_Tshock = 0.0; min_Tshock = 1d8
         max_Thp = 0.0; min_Thp = 1d8
         max_Facc = 0.0; min_Facc = 1d8
      endif

      !deallocated after non-LTE loops
      allocate(Itot(n_lambda, 1, nb_proc))

      return
   end subroutine setup_image_grid

   subroutine atom_line_transfer()
   ! --------------------------------------------------------------------------- !
   ! This routine initialises the necessary quantities for atomic line transfer
   ! and calls the appropriate routines for LTE or NLTE transfer.
   ! --------------------------------------------------------------------------- !
      integer  :: ne_initial, ibin, iaz, nat
      logical :: lelectron_read
      real(kind=dp) :: dne, v_char

      lnon_lte_loop = .false.
      omp_chunk_size = max(nint( 0.01 * n_cells / nb_proc ),1)
      mem_alloc_tot = 0
      if (lsed) then
         allocate(tab_lambda_sed(size(tab_lambda)))
         tab_lambda_sed = tab_lambda
      endif

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
      !used for the extension of Voigt profiles
      call set_max_damping()

      if( Nactiveatoms > 0) then

         call compute_max_relative_velocity(v_char)

         !every indexes are defined on the lambda frequency.
         !So even when lambda is passed as an argument, it is expected
         !that the indexes correspond at that grid.
         !Even if we split the transfer un group of wavelength that's probably best to keep indexes.
         !otherwise, we pass group wavelength and check contributing opacities for each group
         !but it is like indexes at the end?
         deallocate(tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda)
         call make_wavelengths_nlte(tab_lambda_nm,vmax_overlap=v_char)
         n_lambda = size(tab_lambda_nm)
         tab_lambda = tab_lambda_nm * m_to_km

         ! !allocate quantities in space and for this frequency grid
         call alloc_atom_opac(n_lambda, tab_lambda_nm)

         call nlte_loop_mali()
         ! call solve_for_nlte_pops
         call write_opacity_emissivity_bin(n_lambda,tab_lambda_nm)
         if (lexit_after_nonlte_loop) return

      end if !active atoms

      if (lmodel_1d) then
         call spectrum_1d()
         !deallocate and exit code
         return !from atomic transfer!
      endif


      !re alloc lambda here.
      call setup_image_grid()
      write(*,*) "Computing emerging flux..."
      do ibin=1,RT_n_incl
         do iaz=1,RT_n_az
            call emission_line_map(ibin,iaz)
         end do
      end do

      if (laccretion_shock) then
         !3 + lg(Facc) = lg(Facc) erg/cm2/s
         if (max_Facc>0) write(*,'("max(lgFacc)="(1F7.3)" W/m^2; min(lgFacc)="(1F7.3)" W/m^2")') & 
            log10(max_Facc), log10(min_Facc)
         if (max_Tshock>0) write(*,'("max(Tshock)="(1I8)" K; min(Tshock)="(1I8)" K")') &
            nint(max_Tshock), nint(min_Tshock)
         if (max_Thp>0) write(*,'("max(Thp)="(1I6)" K; min(Thp)="(1I6)" K")') &
            nint(max_Thp), nint(min_Thp)
      endif

      if (RT_line_method==1) then
         call write_total_flux()
      else
         do nat=1, n_atoms
            if (atoms(nat)%p%lline) call write_atomic_maps(atoms(nat)%p) !onyly with RT = 2
         enddo
      endif
      write(*,*) " ..done"

      write(*,*) "TO DO: add deallocation here"
      call dealloc_atom_opac()
      deallocate(Itot)

      return
   end subroutine atom_line_transfer


  subroutine intensite_pixel_atom(id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
   ! -------------------------------------------------------------- !
   ! stellar map ?
   ! TO DO line by line
   ! Magnetic fields
   ! -------------------------------------------------------------- !

      integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
      real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
      real(kind=dp), intent(in) :: pixelsize,u,v,w
      integer, parameter :: maxSubPixels = 32
      real(kind=dp) :: x0,y0,z0,u0,v0,w0
      real(kind=dp), dimension(N_lambda) :: Iold, I0
      real(kind=dp), dimension(3) :: sdx, sdy
      real(kind=dp):: npix2, diff, normF
      real(kind=dp), parameter :: precision = 1d-2
      integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter
      logical :: lintersect, labs
      integer :: kr, krr, nat
      type (AtomType), pointer :: atom

      labs = .false.
      u0 = -u ; v0 = -v ; w0 = -w

      iray = 1
      ri = 2*n_rad ; zj=1 ; phik=1

      ! le nbre de subpixel en x est 2^(iter-1)
      subpixels = 1
      iter = 1
      diff = 0.
      I0 = 0.0_dp

      infinie : do ! Boucle infinie tant que le pixel n'est pas converge
         npix2 =  real(subpixels)**2
         Iold = I0
         I0 = 0.0_dp
         ! Vecteurs definissant les sous-pixels
         sdx(:) = dx(:) / real(subpixels,kind=dp)
         sdy(:) = dy(:) / real(subpixels,kind=dp)

       !      !L'obs est en dehors de la grille
       !      ri = 2*n_rad ; zj=1 ; phik=1

         ! Boucle sur les sous-pixels qui calcule l'intensite au centre
         ! de chaque sous pixel
         do i = 1,subpixels
            do j = 1,subpixels
            ! Centre du sous-pixel
               x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
               y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
               z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)
               ! On se met au bord de la grille : propagation a l'envers
               call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
               if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
                  call integ_ray_atom(id,icell,x0,y0,z0,u0,v0,w0,iray,labs,n_lambda,tab_lambda_nm)

                  I0 = I0 + Itot(:,iray,id)

                 !else !Outside the grid, no radiation flux
               endif
            end do !j
         end do !i

         I0 = I0 / npix2

         if (iter < n_iter_min) then
            ! On itere par defaut
            subpixels = subpixels * 2
         else if (iter >= n_iter_max) then
            ! On arrete pour pas tourner dans le vide
            exit infinie
         else
            ! On fait le test sur a difference
            diff = maxval( abs(I0 - Iold) / (I0 + 1e-300_dp) )
            ! There is no iteration for Q, U, V, assuming that if I is converged, then Q, U, V also.
            ! Can be added and then use diff with I = (sqrt(I**2 + Q**2 + U**2 + V**2))
            if (diff > precision ) then
               ! On est pas converge
               subpixels = subpixels * 2
            else
             !write(*,*) "Pixel converged", ipix, jpix, i, j, iter, diff
               exit infinie
            end if
         end if ! iter
         iter = iter + 1
      end do infinie

      ! Prise en compte de la surface du pixel (en sr)

      ! Flux out of a pixel in W/m2/Hz/pix
      normF = ( pixelsize / (distance*pc_to_AU) )**2

      if (RT_line_method==1) then
         Flux_total(:,ibin,iaz,id) = Flux_total(:,ibin,iaz,id) + I0(:) * normF
      else
         do nat=1,N_atoms
            atom => atoms(nat)%p
            do kr=1,atom%nTrans_rayTracing
               krr = atom%ij_to_trans(atom%i_Trans_rayTracing(kr),atom%j_Trans_rayTracing(kr))
               atom%lines(krr)%map(ipix,jpix,:,ibin,iaz) = I0(atom%lines(krr)%nb:atom%lines(krr)%nr)*normF
            enddo
            atom => NULL()
         enddo
      endif


      return
   end subroutine intensite_pixel_atom

  subroutine emission_line_map(ibin,iaz)
   ! -------------------------------------------------------------- !
   !
   ! Line emission map in a given direction n(ibin,iaz),
   ! using ray-tracing.
   ! if only one pixel it gives the total Flux.
   ! See: emission_line_map in mol_transfer.f90
   !
   ! -------------------------------------------------------------- !
      integer, intent(in) :: ibin, iaz
      real(kind=dp) :: x0,y0,z0,l,u,v,w
      real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
      real(kind=dp), dimension(3,nb_proc) :: pixelcorner
      real(kind=dp):: taille_pix
      integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max
      integer, parameter :: n_rad_RT = 200, n_phi_RT = 150
      real(kind=dp), dimension(n_rad_RT) :: tab_r
      real(kind=dp):: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
      integer :: ri_RT, phi_RT
      logical :: lresolved

      real(kind=dp) :: z1, z2
      integer, dimension(:), allocatable :: tab_pix_healpix

      write(*,*) "Vector to observer =", real(tab_u_rt(ibin,iaz)),real(tab_v_rt(ibin,iaz)),real(tab_w_rt(ibin))
      write(*,*) "i=", real(tab_RT_incl(ibin)), "az=", real(tab_RT_az(iaz))

      u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin)
      uvw = (/u,v,w/) !vector position

      ! Definition des vecteurs de base du plan image dans le repere universel
      ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
      x = (/cos(tab_RT_az(iaz) * deg_to_rad),sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

      ! Vecteur x image avec PA
      if (abs(ang_disque) > tiny_real) then
         ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
         x_plan_image = rotation_3d(uvw, ang_disque, x)
      else
         x_plan_image = x
      endif

      ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
      y_plan_image = -cross_product(x_plan_image, uvw)

      write(*,*) "x-image =           ", real(x_plan_image(:))
      write(*,*) "y-image =           ", real(y_plan_image(:))

      ! position initiale hors modele (du cote de l'observateur)
      ! = centre de l'image
      l = 10.*Rmax  ! on se met loin ! in AU

      x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
      center(1) = x0 ; center(2) = y0 ; center(3) = z0
      center(:) = center(:) - image_offset_centre(:)

      ! Coin en bas gauche de l'image
      Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)

      if (RT_line_method==1) then !log pixels

         !no sub pixels here
         n_iter_min = 1
         n_iter_max = 1

         dx(:) = 0.0_dp
         dy(:) = 0.0_dp

         i = 1
         j = 1
         lresolved = .false.

         rmin_RT = 0.001_dp * Rmin!max(w*0.9_dp,0.05_dp) * Rmin
         rmax_RT = Rmax * 1.0_dp ! Rmax  * 2.0_dp

         tab_r(1) = rmin_RT
         fact_r = exp( (1.0_dp/(real(n_rad_RT,kind=dp) -1))*log(rmax_RT/rmin_RT) )
         do ri_RT = 2, n_rad_RT
            tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
         enddo

         fact_A = sqrt(pi * (fact_r - 1.0_dp/fact_r)  / n_phi_RT )


         ! Boucle sur les rayons d'echantillonnage
         !$omp parallel &
         !$omp default(none) &
         !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
         !$omp shared(tab_r,fact_A,x,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
         !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz,fact_r)
         id = 1 ! pour code sequentiel

         if (l_sym_ima) then
            cst_phi = pi  / real(n_phi_RT,kind=dp)
         else
            cst_phi = deux_pi  / real(n_phi_RT,kind=dp)
         endif

         !$omp do schedule(dynamic,1)
         do ri_RT=1, n_rad_RT
            !$ id = omp_get_thread_num() + 1
            r = tab_r(ri_RT)
            taille_pix =  fact_A * r ! racine carree de l'aire du pixel
            do phi_RT=1,n_phi_RT ! de 0 + eps à 2pi - eps (eps = pi/n_phi_RT)
               phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)

               pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image
               ! C'est le centre en fait car dx = dy = 0.
               call intensite_pixel_atom(id,ibin,iaz,n_iter_min,n_iter_max, i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
            end do !j
         end do !i
         !$omp end do
         !$omp end parallel

         !-> if star_map set dx and dy to
         !!taille_pix = (map_size/zoom)  ! en AU
         !!dx(:) = x_plan_image * taille_pix
         !!dy(:) = y_plan_image * taille_pix


      elseif (RT_line_method == 2) then !method 2
         ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
         taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
         lresolved = .true.

         dx(:) = x_plan_image * taille_pix
         dy(:) = y_plan_image * taille_pix

         !-> dust_transfer
         !!Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

         if (l_sym_ima) then
            npix_x_max = npix_x/2 + modulo(npix_x,2)
         else
            npix_x_max = npix_x
         endif

         !$omp parallel &
         !$omp default(none) &
         !$omp private(i,j,id) &
         !$omp shared(Icorner,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y) &
         !$omp shared(n_iter_min,n_iter_max,ibin,iaz)

         ! loop on pixels
         id = 1 ! pour code sequentiel
         n_iter_min = 1
         n_iter_max = 1
         !$omp do schedule(dynamic,1)
         do i = 1,npix_x_max
            !$ id = omp_get_thread_num() + 1
            do j = 1,npix_y
               pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
               call intensite_pixel_atom(id,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
            end do !j
         end do !i
         !$omp end do
         !$omp end parallel

     else !method 3, healpix array map !
         call error("healpix pixel map not yet!")
         ! lresolved = .true.
         ! zoom = 1.0 !here


         ! npix_x_max = npix_x - (2**healpix_lorder) * ((2**healpix_lorder) - 1) * 2
         ! allocate(tab_pix_healpix(npix_x_max))
         ! call healpix_ring_northern_pixels(healpix_lorder,tab_pix_healpix)

         ! taille_pix = 1.0 / real(npix_x_max) ! en AU
         ! !with healpix we already pass the centre of each pixel
         ! dx(:) = 0.0_dp
         ! dy(:) = 0.0_dp

         ! !$omp parallel &
         ! !$omp default(none) &
         ! !$omp private(i,id,z1,z2,j) &
         ! !$omp shared(ibin, iaz, pixelcorner,u,v,w,taille_pix)&
         ! !$omp shared(dx, dy, npix_x_max,l,healpix_lorder,tab_pix_healpix)

         ! ! loop on pixels
         ! id = 1 ! pour code sequentiel
         ! !$omp do schedule(dynamic,1)
         ! do i=1, npix_x_max !from healpix
         !    !$ id = omp_get_thread_num() + 1
         !    ! healpix pixel centre
         !    call healpix_ring_mu_and_phi(healpix_lorder, tab_pix_healpix(i), z1, z2)
         !    pixelcorner(:,id) = l*(/sqrt(1.0 - z1*z1)*cos(z2),sqrt(1.0-z1*z1)*sin(z2),z1/)
         !    call intensite_pixel_atom(id,ibin,iaz,1,1,i,1,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
         ! end do !i
         ! !$omp end do
         ! !$omp end parallel

         ! deallocate(tab_pix_healpix)

      end if

      return
   end subroutine emission_line_map

   subroutine spectrum_1d()
      !TO DO log rays
      !TO DO fits format

      !NOTE: cos_theta in that case is minimum at Rmax which is not the limb.
      ! the limb is at p = 1.0 (rr=Rstar). In plan geomtry, cos_theta is minimum at the limb.
      ! the plan-parralel cos_theta is np.sqrt(1.0 - p**2) for p < 1.0.
      integer :: la, j, icell0, id
      logical :: lintersect, labs
      integer, parameter :: Nimpact = 30
      real(kind=dp) :: rr, u,v,w,u0,w0,v0,x0,y0,z0,x(3),y(3),uvw(3), v_char
      real(kind=dp), allocatable :: cos_theta(:), weight_mu(:), p(:)
      real(kind=dp), allocatable ::I_1d(:,:)
      integer :: n_rays_done, ibar, alloc_status

      write(*,*) "**** computing CLV intensity for 1d model..."
      call compute_max_relative_velocity(v_char)
      if (v_char > 0.0) then
         write(*,*) "WARNING spectrum_1d():"
         write(*,*) " velocity fields are present, but Inu(p) is "
         write(*,*) " computed for a single radius of the stellar disc!!"
      endif

      allocate(cos_theta(Nimpact), weight_mu(Nimpact), p(Nimpact))
      !prepare wavelength grid and indexes for continua and lines
      !including max extension due to velocity shiftS.
      call deallocate_wavelengths_gasrt(tab_lambda)
      call dealloc_atom_opac()
      if (lsed) then
         call make_wavelengths_flux(tab_lambda_sed,.true.)
      else
         call make_wavelengths_flux(tab_lambda_nm,.false.)
      endif
      n_lambda = size(tab_lambda_nm)
      tab_lambda = tab_lambda_nm * m_to_km
      call alloc_atom_opac(n_lambda, tab_lambda_nm)
      !should not be allocated already
      allocate(Itot(N_lambda,Nimpact,nb_proc),stat=alloc_status); Itot = 0.0_dp
      if (alloc_status > 0) call error("spectrum_1d: cannot allocate Itot")
      allocate(I_1d(n_lambda,Nimpact))

      u = 0.0_dp
      v = 0.0_dp
      w = 1.0_dp

      id = 1
      n_rays_done = 0
      ibar = 0
      labs = .false.

      uvw = (/u,v,w/) !vector position
      x = (/1.0_dp,0.0_dp,0.0_dp/)
      y = -cross_product(x, uvw)

      ! Ray tracing : on se propage dans l'autre sens
      u0 = -u ; v0 = -v ; w0 = -w
      call gauss_legendre_quadrature(0.0_dp, 1.0_dp, Nimpact, cos_theta, weight_mu)
      I_1d = 0.0
      call progress_bar(0)
      !$omp parallel &
      !$omp default(none) &
      !$omp private(id,icell0,j) &
      !$omp private(x0,y0,z0,lintersect,rr) &
      !$omp shared(cos_theta,p,Itot,N_lambda,tab_lambda_nm,x,y,weight_mu,ibar,n_rays_done)&
      !$omp shared(u,v,w,u0,v0,w0,Rmax,Rmin,labs,move_to_grid)
      !$omp do schedule(dynamic,1)
      do j=1,Nimpact
         !$ id = omp_get_thread_num() + 1

         rr = Rmax * sqrt(1.0 - cos_theta(j)**2)
         p(j) = rr/Rmin

         x0 = 10.0*Rmax*u + rr * y(1)
         y0 = 10.0*Rmax*v + rr * y(2)
         z0 = 10.0*Rmax*w + rr * y(3)

         call move_to_grid(id,x0,y0,z0,u0,v0,w0,icell0,lintersect)
         if (lintersect) then
            call integ_ray_atom(id,icell0,x0,y0,z0,u0,v0,w0,j,labs,n_lambda,tab_lambda_nm)
         endif

         !$omp atomic
         n_rays_done = n_rays_done + 1
         if (real(n_rays_done) > 0.02*ibar*Nimpact) then
            call progress_bar(ibar)
            !$omp atomic
            ibar = ibar+1
         endif
      enddo
      !$omp end do
      !$omp end parallel
      call progress_bar(50)
      I_1d(:,:) = Itot(:,:,1)
      do j=2,nb_proc !direct sum on nb_proc raises a segmentation fault
         I_1d(:,:) = I_1d(:,:) + Itot(:,:,j)
      enddo

      open(1,file="spectrum_1d.txt",status="unknown")
      write(1,*) Nimpact, N_lambda!, atmos_1d%s
      write(1,'(*(1ES17.8E3))') (p(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (cos_theta(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (weight_mu(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (tab_lambda_nm(la), la=1,N_lambda)!1F12.4
      do la=1,N_lambda
            write(1,'(*(1ES17.8E3))') (I_1d(la,j), j=1,Nimpact)
      enddo
      close(1)

      if (loutput_rates) then
         call write_opacity_emissivity_bin(n_lambda, tab_lambda_nm)
      endif
      deallocate(cos_theta, weight_mu, p, I_1d, Itot, tab_lambda_nm, tab_lambda)

   return
   end subroutine spectrum_1d

end module atom_transfer
         ! if (lvoronoi) then
         !    v_char = sqrt( maxval(Voronoi(:)%vxyz(1)**2+Voronoi(:)%vxyz(2)**2+Voronoi(:)%vxyz(3)**2) )
         ! else
         !    v_char = sqrt( maxval(sum(vfield3d**2,dim=2)) )
         ! endif
         ! write(*,'("maximum gradv="(1F12.3)" km/s")') v_char*1d-3
