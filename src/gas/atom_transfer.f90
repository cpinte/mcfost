! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atomic (only) systems using the MALI method
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
module atom_transfer

   use parametres
   use input, only               : lkeplerian, linfall, limb_darkening, mu_limb_darkening, RT_line_method
   use constantes, only : nm_to_m, m_to_km, km_to_m, au_to_m, deg_to_rad, tiny_real, tiny_dp, pi, deux_pi, pc_to_au
   use io_atom, only : read_atomic_models, write_pops_atom
   use wavelengths, only : n_lambda, tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda, n_lambda2, tab_lambda2
   use wavelengths_gas, only : make_wavelengths_nlte, tab_lambda_nm, tab_lambda_cont, n_lambda_cont, &
                                 deallocate_wavelengths_gasrt, make_wavelengths_raytracing, make_wavelengths_flux, Nlambda_max_line
   use elecdensity, only : solve_ne, write_electron, read_electron
   use grid, only : nHtot, pos_em_cellule, lcalc_ne, move_to_grid, vfield3d, icompute_atomRT, ne, Voronoi, r_grid, phi_grid, z_grid
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops
   use atom_type, only : atoms, atomtype, n_atoms, nactiveatoms, activeAtoms, hydrogen, helium, adjust_cswitch_atoms, maxval_cswitch_atoms, lcswitch_enabled
   use init_mcfost, only :  nb_proc
   use gas_contopac, only : background_continua_lambda
   use opacity_atom, only : alloc_atom_opac, Itot, psi, dealloc_atom_opac, xcoupling
   use see, only : n_new, lcell_converged, ne_new, ngpop, alloc_nlte_var, dealloc_nlte_var, frac_limit_pops, &
                  rate_matrix, init_rates, update_populations, accumulate_radrates_mali
   use optical_depth, only : integ_ray_atom
   use utils, only : cross_product, gauss_legendre_quadrature, progress_bar, rotation_3d
   use dust_ray_tracing, only    : RT_n_incl, RT_n_az, init_directions_ray_tracing,tab_u_RT, tab_v_RT, tab_w_RT, &
                                   tab_RT_az,tab_RT_incl
   use stars, only               : intersect_stars, laccretion_shock, max_Tshock, min_Tshock
   use output, only : allocate_atom_maps, flux_total, write_total_flux, write_atomic_maps
   use mcfost_env, only          : dp, time_tick, time_max
   use molecular_emission, only  : ds
   use messages, only : error, warning

   use healpix_mod
   !$ use omp_lib

   implicit none
                    
   integer :: omp_chunk_size
   real(kind=dp), allocatable, dimension(:) :: tab_lambda_Sed

   contains

!TO DO 
!      add iterate ne
!      add Ng
!      add Trad, Tion
   subroutine nlte_loop_mali()
   ! ----------------------------------------------------------------------- !
   ! Descriptor
   ! ----------------------------------------------------------------------- !
      use  naleat, only : seed, stream, gtype
#include "sprng_f.h"
      integer :: etape, etape_start, etape_end, iray, n_rayons
      integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
      integer, dimension(nb_proc) :: max_n_iter_loc
      integer, parameter :: maxIter = 1000, maxIter_loc = 100
      logical :: lfixed_Rays, lconverged, lconverged_loc, lprevious_converged
      real :: rand, rand2, rand3, fac_etape
      real(kind=dp) :: precision
      integer, parameter :: n_rayons_max = 1 !ray-by-ray except for subiterations
                                             !but in that case only itot and phi_loc and ds needs
                                             !to have n_rayons_sub size.
      integer :: n_rayons_sub
      real(kind=dp), dimension(:), allocatable :: xmu,wmu,xmux,xmuy
      real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, weight
      real(kind=dp) :: diff, diff_old, dne, dN, dN1
      real(kind=dp), allocatable :: dTM(:), dM(:), Tion_ref(:), Tex_ref(:)

      logical :: labs, lsubiteration = .false.
      logical :: l_iterate, l_iterate_ne, lupdate_ne_other_nlte
      ! logical :: accelerated, ng_rest, lng_turned_on
      ! integer :: iorder, i0_rest, n_iter_accel, iacc

      integer :: nact, imax, icell_max, icell_max_2
      integer :: icell, ilevel, Ng_Ndelay
      logical :: lng_turned_on

      type (AtomType), pointer :: atom
      integer(kind=8) :: mem_alloc_local
      real(kind=dp) :: diff_cont, conv_speed, conv_acc

      real :: time_iteration, time_nlte
      integer :: count_start, count_end
      integer :: ibar, n_cells_done
      integer, parameter :: n_iter_counted = 1

      ! -------------------------------- INITIALIZATION -------------------------------- !
      write(*,*) '-------------------------- NON-LTE LOOP ------------------------------ '
      !non-LTE mode
      labs = .true.
      id = 1
      if (lforce_lte) lsubiteration = .false.

      l_iterate_ne = .false.
      lupdate_ne_other_nlte = (hydrogen%active.and.NactiveAtoms > 1).or.(.not.hydrogen%active)
      if (associated(helium)) then

         lupdate_ne_other_nlte = (helium%active.and.hydrogen%active.and.Nactiveatoms > 2).or.&
                                 (.not.hydrogen%active.and.helium.active.and.Nactiveatoms > 1).or.&
                                 (hydrogen%active.and..not.helium%active.and.NactiveAtoms > 1)
      endif
      !Use non-LTE pops for electrons and to write to file.
      do nact=1,NactiveAtoms
         if (.not.activeatoms(nact)%p%nltepops) activeatoms(nact)%p%nltepops = .true.
      enddo

      allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
      allocate(dTM(Nactiveatoms)); dM=0d0 !keep tracks of Tex for all cells for each atom
      allocate(Tex_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tex for all cells for each line of each atom
      allocate(Tion_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tion for all cells for each cont of each atom
      diff_old = 1.0_dp
      dM(:) = 1.0_dp

      !-> negligible
      mem_alloc_local = 0
      mem_alloc_local = mem_alloc_local + sizeof(dM)+sizeof(dTm)+sizeof(Tex_ref)+sizeof(Tion_ref)


      !How many steps and which one
      etape_start = istep_start
      n_rayons = healpix_npix(healpix_lorder)
      n_rayons_sub = 1 !for allocation of phi and Itot on the maximum number of rays!
      if (laccurate_integ) then
         etape_end = 2
         if (lsubiteration) n_rayons_sub = max(n_rayons,N_rayons_mc)
      else
         etape_end = 1
         if (lsubiteration) n_rayons_sub = n_rayons
         if (istep_start==2) then 
            etape_end = 2
            if (lsubiteration) n_rayons_sub = N_rayons_mc
         endif
      endif

      if (allocated(stream)) deallocate(stream)
      allocate(stream(nb_proc),stat=alloc_status)
      if (alloc_status > 0) call error("Allocation error stream")
      if (allocated(ds)) deallocate(ds)
      if (lsubiteration) then
         allocate(ds(n_rayons_sub,nb_proc))
      else
         allocate(ds(n_rayons_max,nb_proc))
      endif
      !-> negligible
      mem_alloc_local = mem_alloc_local + sizeof(ds) + sizeof(stream)
      call alloc_nlte_var(n_rayons_max,n_rayons_sub,lsubiteration)

    
      lfixed_rays = .true.
      iray_start = 1

      ! --------------------------- OUTER LOOP ON STEP --------------------------- !

      step_loop : do etape=etape_start, etape_end

         write(*,*) ""
         if (etape==1) then
            precision = dpops_max_error
            !use stepan if healpix_lorder < 2?
            allocate(xmux(n_rayons),xmu(n_rayons),xmuy(n_rayons),wmu(n_rayons))
            wmu(:) = healpix_weight(healpix_lorder)
            write(*,'(" ****-> Using "(1I8)" pixels for healpix, resolution of "(1F12.3)" degrees")') n_rayons,  healpix_angular_resolution(healpix_lorder)
            call healpix_sphere(healpix_lorder,xmu,xmux,xmuy)
         else if (etape==2) then
            lng_turned_on = .false.
            fac_etape = 0.1
            if (etape_start==1) then
               precision = min(1e-1,10.0*dpops_max_error)
            else
               precision = dpops_max_error
            endif
            n_rayons = N_rayons_mc
            write(*,'(" ****-> Using "(1I4)" rays for Monte-Carlo step, ~resolution of "(1F12.3)" degrees")') n_rayons, 360.0 * sqrt(pi/real(n_rayons))/pi
            allocate(wmu(n_rayons))
            wmu(:) = 1.0_dp / real(n_rayons,kind=dp)
         else
            call error("(n_lte_loop_mali) etape unkown")
         end if
         write(*,*) ""
         call system_clock(count_start,count_rate=time_tick,count_max=time_max)

         lconverged = .false.
         lprevious_converged =  lforce_lte
         lcell_converged(:) = .false.
         lng_turned_on = .false.
         n_iter = 0
         dne = 0.0_dp
         conv_speed = 0.0
         conv_acc = 0.0
         diff_old = 0.0
         time_iteration = 0.

         !***********************************************************!
         ! *************** Main convergence loop ********************!
         do while (.not.lconverged)

            n_iter = n_iter + 1
            write(*,'(" *** Iteration #"(1I3)"; step #"(1I1))') n_iter, etape
            ibar = 0
            n_cells_done = 0

            if (lfixed_rays) then
               stream = 0.0
               stream(:) = [(init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT),i=1,nb_proc)]
            end if

            max_n_iter_loc = 0

            if( n_iterate_ne > 0 ) then
               l_iterate_ne = ( mod(n_iter,n_iterate_ne)==0 ) .and. (n_iter>ndelay_iterate_ne)
            endif

            call progress_bar(0)
            !$omp parallel &
            !$omp default(none) &
            !$omp private(id,icell,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02,argmt)&
            !$omp private(l_iterate,weight,diff,lconverged_loc,n_iter_loc)&
            !$omp shared(etape,lforce_lte,n_cells,voronoi,r_grid,z_grid,phi_grid,n_rayons,xmu,wmu,xmux,xmuy) &
            !$omp shared(pos_em_cellule,labs,n_lambda,tab_lambda_nm, icompute_atomRT,max_n_iter_loc) &
            !$omp shared(stream,n_rayons_mc,lvoronoi,ibar,n_cells_done,l_iterate_ne)
            !$omp do schedule(static,omp_chunk_size)
            do icell=1, n_cells
               !$ id = omp_get_thread_num() + 1
               l_iterate = (icompute_atomRT(icell)>0)

               if (l_iterate) then

                  !also init atom%Gamma
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
                           call xcoupling(id, icell,1, .false.)
                           call accumulate_radrates_mali(id, icell,1, weight)
                        endif

                     enddo


                  elseif ( etape == 2) then
                   ! Position aleatoire dans la cellule
                     do iray=1,n_rayons

                        rand  = sprng(stream(id))
                        rand2 = sprng(stream(id))
                        rand3 = sprng(stream(id))

                        call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

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
                           call xcoupling(id, icell,1, .false.)
                           call accumulate_radrates_mali(id, icell,1, weight)
                        endif

                     enddo !iray
                  else
                     call error("etape inconnue")
                  end if !etape

                  call rate_matrix(id)
                  call update_populations(id, icell, (l_iterate_ne.and.icompute_atomRT(icell)==1), diff)

                  ! n_iter_loc = 0
                  ! lconverged_loc = .not.lsubiteration
                  ! npl(1:hydrogen%Nlevel,1,id) = hydrogen%n(:,icell)
                  ! !this is n_new that is to be updated here, not atom%n ! ("old pops")
                  ! do while (.not.lconverged_loc)
                  !    n_iter_loc = n_iter_loc + 1

                  !    hydrogen%n(:,icell) = n_new(1:hydrogen%Nlevel,1,icell)

                  !    !psi is updated.
                  !    !currenlty includes only the lines opacity!
                  !    !loop over rays of this step (not n_rayons_sub)
                  !    do iray=1,n_rayons
                  !       call xcoupling(id, icell, iray, .true.)
                  !       call accumulate_radrates_mali(id, icell, iray, weight)                        
                  !    enddo
                  !    call rate_matrix(id)
                  !    !electronic density fixed for the sub-iterations (because it is already an iterative loop)
                  !    !do the electron density iterations after the sub-iterations ??
                  !    call update_populations(id, icell,.false.,diff)

                  !    ! write(*,*) "subit:", id, n_iter_loc, diff
                  !    if (real(diff) < dpops_sub_max_error) lconverged_loc = .true.

                  !    if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
                  !    if (n_iter_loc == maxIter_loc) lconverged_loc = .true.
                  ! enddo
                  ! hydrogen%n(:,icell) = npl(1:hydrogen%Nlevel,1,id)

                  ! if ( (lNg_acceleration .and. lng_turned_on).and.&
                  !       (maxval_cswitch_atoms()==1.0_dp).and.(.not.lprevious_converged) ) then

                  ! endif

               end if !icompute_atomRT
             
               ! Progress bar
               !$omp atomic
               n_cells_done = n_cells_done + 1
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
            ! ************************** Ng's **************************!
 
            !***********************************************************!
            ! ******************* checkpointing ************************!

            !***********************************************************!
            ! *******************  update ne metals  *******************!

            !-> evaluate ne for other non-LTE atoms not included in the non-LTE ionisation scheme ?
            !with fixed non-LTE populations and non-LTE contributions
            if (l_iterate_ne) then
               !add in solve_ne_loc cell-by-cell
               !add in other atom in non_lte_charge_conservation only for the ionisation frac.
               !evaluate non-lte donors other than h + he
               !update lte opac and related quantities
               !damping and profiles for instance
               !re update nHmin(icell) = nH_minus(icell)

            end if

            !***********************************************************!
            ! ****************  Convergence checking  ******************!

            id = 1
            dM(:) = 0.0 !all pops
            diff_cont = 0
            diff = 0.0
            !$omp parallel &
            !$omp default(none) &
            !$omp private(id,icell,l_iterate,dN1,dN,ilevel,nact,atom)&
            !$omp shared(diff, diff_cont, dM,n_new,Activeatoms,lcell_converged)&
            !$omp shared(icompute_atomRT,n_cells,precision,NactiveAtoms,nhtot)
            !$omp do schedule(dynamic,1)
            cell_loop2 : do icell=1,n_cells
               !$ id = omp_get_thread_num() + 1
               l_iterate = (icompute_atomRT(icell)>0)

               if (l_iterate) then

                  !Local only
                  dN = 0.0 !for all levels of all atoms of this cell
                  do nact=1,NactiveAtoms
                     atom => ActiveAtoms(nact)%p

                     do ilevel=1,atom%Nlevel
                        if ( n_new(nact,ilevel,icell) >= frac_limit_pops * atom%Abund*nHtot(icell) ) then
                           dN1 = abs(1d0-atom%n(ilevel,icell)/n_new(nact,ilevel,icell))
                           dN = max(dN1, dN)
                           dM(nact) = max(dM(nact), dN1)
                        endif
                     end do !over ilevel
                     diff_cont = max(diff_cont, abs(1.0 - atom%n(atom%Nlevel,icell)/(1d-50 + n_new(nact,atom%Nlevel,icell) )))


                     atom => NULL()
                  end do !over atoms

                  !compare for all atoms and all cells
                  diff = max(diff, dN) ! pops
                  lcell_converged(icell) = (real(diff) < precision)!.and.(dne < precision)

                  !Re init for next iteration if any
                  do nact=1, NactiveAtoms
                     atom => ActiveAtoms(nact)%p
                     atom%n(:,icell) = n_new(nact,1:atom%Nlevel,icell)
                     atom => NULL()
                  end do

               end if !if l_iterate
            end do cell_loop2 !icell
            !$omp end do
            !$omp end parallel
  
            if (n_iter > 1) then
               conv_speed = (diff - diff_old) !<0 converging.
               conv_acc = conv_speed - conv_acc !>0 accelerating
            endif
          

            if (maxval(max_n_iter_loc)>1) write(*,'(" -> "(1I10)" sub-iterations")') maxval(max_n_iter_loc)
            do nact=1,NactiveAtoms
               write(*,'("             Atom "(1A2))') ActiveAtoms(nact)%p%ID
               write(*,'("   >>> dpop="(1ES17.8E3))') dM(nact)
            enddo
            if (dne /= 0.0_dp) write(*,'("   >>> dne="(1ES17.8E3))') dne
            write(*,*) ""
            write(*,'(" <<->> diff="(1ES17.8E3)," old="(1ES17.8E3))') diff, diff_old
            write(*,'("   ->> diffcont="(1ES17.8E3))') diff_cont
            write(*,'("   ->> speed="(1ES17.8E3)"; acc="(1ES17.8E3))') conv_speed, conv_acc
            write(*,"('   ->> Unconverged cells #'(1I6), ' fraction :'(1F12.3)' %')") &
               size(pack(lcell_converged,mask=(lcell_converged==.false.).and.(icompute_atomRT>0))), &
               100.*real(size(pack(lcell_converged,mask=(lcell_converged==.false.).and.(icompute_atomRT>0)))) / &
               real(size(pack(icompute_atomRT,mask=icompute_atomRT>0)))
            write(*,*) "      -------------------------------------------------- "
            diff_old = diff
            conv_acc = conv_speed
               

            if ((diff < precision).and.maxval_cswitch_atoms()==1.0_dp)then
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
               if (n_iter > maxIter) then
                  call warning("not enough iterations to converge !!")
                  lconverged = .true.
               end if
            end if

            !***********************************************************!
            ! **********  Ng's acceleration administration *************!

            if (lng_acceleration) then
          	!be sure we are converging before extrapolating
          	   if ( (n_iter>1).and.(conv_speed <= 0).and.(-conv_speed < 1d-3) ) then
          		   if (.not.lng_turned_on) then
          			   lng_turned_on = .true.
          			   Ng_Ndelay = n_iter
          		   endif
            	endif
            endif

            !***********************************************************!
            ! ********** timing and checkpointing **********************!

            call system_clock(count_end,count_rate=time_tick,count_max=time_max)
            if (count_end < count_start) then
               time_nlte=real(count_end + (1.0 * time_max)- count_start)/real(time_tick)
            else
               time_nlte=real(count_end - count_start)/real(time_tick)
            endif

            !average time of a single iteration for this step !
            if (n_iter <= n_iter_counted) then
               time_iteration = time_iteration + time_nlte / real(n_iter_counted)
            endif


            if (lsafe_stop) then
               if ((time_nlte + time_iteration >=  safe_stop_time).and.(n_iter >= n_iter_counted)) then
                  lconverged = .true.
                  lprevious_converged = .true.
                  call warning("Time limit would be exceeded, leaving...")
                  write(*,*) " time limit:", mod(safe_stop_time/60.,60.) ," min"
                  write(*,*) " ~<time> etape:", mod(n_iter * time_iteration/60.,60.), ' <time iter>=', &
                                                mod(time_iteration/60.,60.)," min"
                  write(*,*) " ~<time> etape (cpu):", mod(n_iter * time_iteration * nb_proc/60.,60.), " min"
                  write(*,*) ' time =',mod(time_nlte/60.,60.), " min"
                  lexit_after_nonlte_loop = .true.
                  exit step_loop
               endif
            endif

         !***********************************************************!
         !***********************************************************!

         end do !while
  
         if (n_iter >= n_iter_counted) then
            write(*,'("<time iteration>="(1F6.3)" min; <time step>="(1F6.3)" min")') mod(time_iteration/60.,60.), mod(n_iter * time_iteration/60.,60.)
            write(*,'(" --> ~time step (cpu)="(1F6.3)" min")') mod(n_iter * time_iteration * nb_proc/60.,60.)
         endif

         if (etape==1) deallocate(xmux,xmuy,xmu)
         if (allocated(wmu)) deallocate(wmu)

         write(*,'("  <step #"(1I1)" done! >>threshold="(1ES17.8E3))') etape,  precision

      end do step_loop

      ! -------------------------------- CLEANING ------------------------------------------ !

      if (n_iterate_ne > 0) then
         write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
         call write_electron
      endif

      do nact=1,N_atoms
         call write_pops_atom(Atoms(nact)%p)
      end do


      call dealloc_nlte_var()
      deallocate(dM, dTM, Tex_ref, Tion_ref)
      deallocate(stream,ds)

      ! --------------------------------    END    ------------------------------------------ !
      write(*,*) '-------------------------- ------------ ------------------------------ '

      return
   end subroutine nlte_loop_mali

   subroutine setup_image_grid()
      use wavelengths_gas, only : compute_line_bound
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
         if (npix_y == 1) RT_line_method = 3
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
         max_Tshock = 0.0
         min_Tshock = 1d8
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
      real(kind=dp) :: dne, v_char, max_vel_shift, v1, v2

      write(*,*) " *** BIG WARNING **** "
      write(*,*) " !!!!! CHECK FOR BILINEAR INTERP in utils.f90 ! !!!!"
      write(*,*) "check solve ne" 
      write(*,*) " ******************** "

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

      if( Nactiveatoms > 0) then

         v_char = sqrt( maxval(sum(vfield3d**2,dim=2)) )
         max_vel_shift = 0.0
         do ibin=1,n_cells
            v1 = sqrt(sum(vfield3d(ibin,:)**2))
            do iaz=1,n_cells
               v2 = sqrt(sum(vfield3d(iaz,:)**2))
               if ((icompute_atomRT(ibin)>0).and.(icompute_atomRT(iaz)>0)) then
                  max_vel_shift = max(max_vel_shift,abs(v1-v2))
               endif
            enddo
         enddo
         v_char = max_vel_shift
         write(*,*) "v_char", v_char * 1d-3
         !because for non-LTE we compute the motion of each cell wrt to the actual cell
         !so the maximum shift is what we have if all cells move at their max speed (v(icell))
         ! stop

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

         !allocate quantities in space and for this frequency grid
         call alloc_atom_opac(n_lambda, tab_lambda_nm)

         call nlte_loop_mali()
         if (lexit_after_nonlte_loop) return

      end if !active atoms


      if (lmodel_1d) then
         call spectrum_1d()
         !deallocate and exit code
         return !from atomic transfer!
      endif


      !re alloc lambda here.
      call setup_image_grid()
      write(*,*) "Computing emission flux map..."
      do ibin=1,RT_n_incl
         do iaz=1,RT_n_az
            call emission_line_map(ibin,iaz)
         end do
      end do
      if (RT_line_method==1) then
         call write_total_flux()
      else
         do nat=1, n_atoms
            if (atoms(nat)%p%lline) call write_atomic_maps(atoms(nat)%p) !onyly with RT = 2
         enddo
      endif
      write(*,*) " ..done"


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
      integer, parameter :: n_rad_RT = 150, n_phi_RT = 360
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
      integer, parameter :: Nimpact = 10
      real(kind=dp) :: rr, u,v,w,u0,w0,v0,x0,y0,z0,x(3),y(3),uvw(3)
      real(kind=dp), allocatable :: cos_theta(:), weight_mu(:), p(:)
      real(kind=dp), allocatable ::I_1d(:,:)
      integer :: n_rays_done, ibar, alloc_status

      write(*,*) "**** computing CLV intensity for 1d model..."

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
      !$omp shared(u,v,w,u0,v0,w0,Rmax,Rmin,labs)
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
      write(1,*) Nimpact, N_lambda
      write(1,'(*(1ES17.8E3))') (p(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (cos_theta(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (weight_mu(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (tab_lambda_nm(la), la=1,N_lambda)!1F12.4
      do la=1,N_lambda
            write(1,'(*(1ES17.8E3))') (I_1d(la,j), j=1,Nimpact)
      enddo
      close(1)

      deallocate(cos_theta, weight_mu, p, I_1d, Itot)

   return
   end subroutine spectrum_1d

end module atom_transfer