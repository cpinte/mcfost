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
   use grid, only : lcalc_ne, move_to_grid, vfield3d, icompute_atomRT, ne
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops
   use atom_type, only : atoms, atomtype, n_atoms, nactiveatoms, activeAtoms, hydrogen, helium
   use init_mcfost, only :  nb_proc
   use gas_contopac, only : background_continua_lambda
   use opacity_atom, only : alloc_atom_opac, Itot, psi, dealloc_atom_opac, xcoupling, uji_down, chi_up, chi_down, eta_atoms, phi_loc
   use optical_depth, only : integ_ray_atom
   use utils, only : cross_product, gauss_legendre_quadrature, progress_bar, rotation_3d
   use dust_ray_tracing, only    : RT_n_incl, RT_n_az, init_directions_ray_tracing,tab_u_RT, tab_v_RT, tab_w_RT, &
                                   tab_RT_az,tab_RT_incl
   use stars, only               : intersect_stars, laccretion_shock, max_Tshock, min_Tshock
   use output, only : allocate_atom_maps, flux_total, write_total_flux, write_atomic_maps
   use mcfost_env, only          : dp, time_begin, time_end, time_tick, time_max
   use molecular_emission, only  : ds
   use messages, only : error

   use healpix_mod
   !$ use omp_lib

   implicit none
                    
   integer :: omp_chunk_size
   real(kind=dp) :: dne
   real(kind=dp), allocatable, dimension(:) :: tab_lambda_Sed

   contains

!TO DO Ng
!Trad, Tion
   subroutine nlte_loop_mali()
   ! ----------------------------------------------------------------------- !
   ! Descriptor
   ! ----------------------------------------------------------------------- !
      use  naleat, only : seed, stream, gtype
#include "sprng_f.h"
      integer :: etape, etape_start, etape_end, iray, n_rayons
      integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
      integer, dimension(nb_proc) :: max_n_iter_loc
      logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged
      real :: rand, rand2, rand3, precision, fac_etape
      integer, parameter :: n_rayons_max = 1
      integer :: n_rayons_sub

      real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, weight
      real(kind=dp) :: diff, norme, dN, dN1
      real(kind=dp) :: dT, dN2, dN3, dN4, diff_old, dne
      real(kind=dp), allocatable :: dTM(:), dM(:), Tion_ref(:), Tex_ref(:)

      logical :: labs, lsubiteration = .false.
      logical :: l_iterate, l_iterate_ne, update_ne_other_nlte
      ! logical :: accelerated, ng_rest, lng_turned_on
      ! integer :: iorder, i0_rest, n_iter_accel, iacc

      integer :: nact, imax, icell_max, icell_max_2
      integer :: icell, ilevel, imu, iphi, Nmaxlevel, NMaxtrans, Nmaxline

      type (AtomType), pointer :: atom
      integer(kind=8) :: mem_alloc_local = 0
      real(kind=dp) :: diff_cont, conv_speed, conv_acc

      real(kind=dp), dimension(:), allocatable :: lcell_converged, ne_new, ngpop!check size

      integer :: ibar, n_cells_done
      integer, parameter :: n_iter_counted = 1

      !non-LTE mode
      labs = .true.
      id = 1

      NmaxLevel = 0
      Nmaxtrans = 0
      do i=1, NactiveAtoms
         NmaxLevel = max(Nmaxlevel, ActiveAtoms(i)%p%Nlevel)
         NmaxTrans = max(NmaxTrans, ActiveAtoms(i)%p%Ntr)
         Nmaxline = max(Nmaxline, ActiveAtoms(i)%p%Nline)
      enddo

      !time for individual steps + check the time from the checkpointing time if any
      !and if it is set lconverged = .true. and .lprevious converged == true
      call system_clock(time_begin,count_rate=time_tick,count_max=time_max)

      update_ne_other_nlte = (hydrogen%active.and.NactiveAtoms > 1).or.(.not.hydrogen%active)
      if (associated(helium)) then

         update_ne_other_nlte = (helium%active.and.hydrogen%active.and.Nactiveatoms > 2).or.&
                                 (.not.hydrogen%active.and.helium.active.and.Nactiveatoms > 1).or.&
                                 (hydrogen%active.and..not.helium%active.and.NactiveAtoms > 1)
      endif

      allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
      allocate(dTM(Nactiveatoms)); dM=0d0 !keep tracks of Tex for all cells for each atom
      allocate(Tex_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tex for all cells for each line of each atom
      allocate(Tion_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tion for all cells for each cont of each atom
      diff_old = 1.0_dp
      dM(:) = 1.0_dp

      mem_alloc_local = mem_alloc_local + sizeof(dM)+sizeof(dTm)+sizeof(Tex_ref)+sizeof(Tion_ref)

      ! if (lNg_acceleration) then
      !    n_iter_accel = 0
      !    i0_rest = 0
      !    ng_rest = .false.
      !    iacc = 0
      !    allocate(ngpop(n_cells * NmaxLevel, iNg_Norder+2, NactiveAtoms), stat=alloc_status)
      !    if (alloc_status > 0) then
      !       call error("Cannot allocate Ng table !")
      !    endif
      !    mem_alloc_local = mem_alloc_local + sizeof(ngpop)
      !    lng_turned_on = .false.
      ! endif

      if (allocated(stream)) deallocate(stream)
      allocate(stream(nb_proc),stat=alloc_status)
      if (alloc_status > 0) call error("Allocation error stream")
      if (allocated(ds)) deallocate(ds)
      allocate(ds(1,nb_proc))

      !How many steps and which one
      etape_start = istep_start
      n_rayons_sub = 100 !put healpix here
      if (laccurate_integ) then
         etape_end = 2
         n_rayons_sub = 100 !put mc here
      else
         etape_end = 1
         if (istep_start==2) then 
            etape_end = 2
            n_rayons_sub = 100 !put mc here
         endif
      endif

      !-> allocate space for non-LTE only quantities !
      allocate(ne_new(n_cells), stat=alloc_status)
      ne_new(:) = ne(:)
      if (alloc_status > 0) call error("Allocation error ne_new")
      write(*,*) " size ne_new:", sizeof(ne_new) / 1024./1024.," MB"
      ! allocate(n_new(NactiveAtoms,Nmaxlevel,n_cells),stat=alloc_status)
      ! if (alloc_status > 0) call error("Allocation error n_new")
      ! write(*,*) " size atom%n, n_new:", 2*sizeof(n_new) / 1024./1024.," MB"
      ! n_new(:,:,:) = 0.0_dp
      allocate(psi(n_lambda, n_rayons_max, nb_proc), stat=alloc_status); psi = 0.0_dp
      write(*,*) " size psi:", sizeof(psi) / 1024./1024.," MB"
      if (alloc_Status > 0) call error("Allocation error psi in nlte loop")
      allocate(eta_atoms(n_lambda,NactiveAtoms,nb_proc),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error eta_atoms")
      write(*,*) " size eta_atoms:", sizeof(eta_atoms) / 1024./1024.," MB"
      allocate(Uji_down(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error Uji_down")
      write(*,*) " size cross-coupling:", 3 * sizeof(Uji_down) / 1024./1024./1024.," GB"
      allocate(chi_down(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error chi_down")
      allocate(chi_up(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error chi_up")
      allocate(lcell_converged(n_cells),stat=alloc_status)
      if (alloc_Status > 0) call error("Allocation error lcell_converged")
      write(*,*) " size lcell_converged:", sizeof(lcell_converged) / 1024./1024./1024.," GB"
    
      if (lNg_acceleration) write(*,*) " size ngpop:", sizeof(ngpop)/1024./1024., " MB"

      if (lsubiteration) then
         allocate(Itot(n_lambda,n_rayons_sub,nb_proc),stat=alloc_status)
         if (alloc_Status > 0) call error("Allocation error Itot")
         allocate(phi_loc(Nlambda_max_line,Nmaxline,NactiveAtoms,n_rayons_sub,nb_proc),stat=alloc_status)
         if (alloc_Status > 0) call error("Allocation error phi_loc")
      else
         allocate(Itot(n_lambda,1,nb_proc),stat=alloc_status)
         if (alloc_Status > 0) call error("Allocation error Itot")
         allocate(phi_loc(Nlambda_max_line,Nmaxline,NactiveAtoms,1,nb_proc),stat=alloc_status)
         if (alloc_Status > 0) call error("Allocation error phi_loc")
      endif
      write(*,*) " size Itot:", sizeof(Itot) / 1024./1024./1024.," GB"
      write(*,*) " size phi_loc:", sizeof(phi_loc) / 1024./1024./1024.," GB"

      mem_alloc_local = mem_alloc_local + sizeof(psi) + sizeof(eta_atoms) + sizeof(Itot) + sizeof(phi_loc) + &
         sizeof(uji_down)+sizeof(chi_down)+sizeof(chi_up) + sizeof(ne_new) + sizeof(lcell_converged) !+ sizeof(n_new)
      write(*,'("Total memory allocated in NLTEloop:"(1F14.3)" GB")') mem_alloc_local / 1024./1024./1024.

      lfixed_rays = .true.
      n_rayons = 1 !always
      iray_start = 1

!       !outer loop on etapes
!       step_loop : do etape=etape_start, etape_end

!          if (etape==1) then
!             time_iteration = 0
!             lprevious_converged = .false.
!             lcell_converged(:) = .false.
!             precision = dpops_max_error

!          else if (etape==2) then
!             time_iteration = 0
!             lprevious_converged = .false.
!             lcell_converged(:) = .false.
!             fac_etape = 0.1
!             if (etape_start==1) then
!                precision = min(1e-1,10.0*dpops_max_error)
!             else
!                precision = dpops_max_error
!             endif

!             !iterate ng even in Monte Carlo ?!
!             ! if (lNg_acceleration) then
!             !    deallocate(ngpop)
!             !    if (allocated(ng_cur)) deallocate(ng_cur)
!             !    lNg_acceleration = .false.
!             ! endif
!          else
!             call ("(n_lte_loop_mali) etape unkown")
!          end if

!          lprevious_converged = lforce_lte
!          lnotfixed_rays = .not.lfixed_rays
!          lconverged = .false.
!          n_iter = 0
!          dne = 0.0_dp
!          conv_speed = 0.0
!          conv_acc = 0.0

!          !for this etape, iterate until convergence
!          do while (.not.lconverged)

!             n_iter = n_iter + 1
!             write(*,*) " -> Iteration #", n_iter, " Step #", etape
!             ibar = 0
!             n_cells_done = 0

!             if (lfixed_rays) then
!                stream = 0.0
!                stream(:) = (init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT), i=1,nb_proc)
!             end if

!             max_n_iter_loc = 0
!           !reset in case of step 3

!             if( n_iterate_ne > 0 ) then
!                l_iterate_ne = ( mod(n_iter,n_iterate_ne)==0 ) .and. (n_iter>ndelay_iterate_ne)
!             endif

!             call progress_bar(0)
!             !$omp parallel &
!             !$omp default(none) &
!             !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dM, dN, dN1,imu,iphi)&
!             !$omp shared()
!             !$omp do schedule(static,omp_chunk_size)
!             do icell=1, n_cells
!                !$ id = omp_get_thread_num() + 1
!                l_iterate = (icompute_atomRT(icell)>0)

!                if (l_iterate) then


!                 call fill_Collision_matrix(id, icell) !computes = C(i,j) before computing rates
!                 call initGamma(id) !init Gamma to C and init radiative rates
!                 !psi(:,iray,id) = 0.0_dp !no use full here


!                 if ( etape == 2) then !ray-by-ray, n_rayons fixed
!                    ! Position aleatoire dans la cellule
!                    do iray=iray_start, iray_start-1+n_rayons


!                       rand  = sprng(stream(id))
!                       rand2 = sprng(stream(id))
!                       rand3 = sprng(stream(id))

!                       call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

!                       ! Direction de propagation aleatoire
!                       rand = sprng(stream(id))
!                       W0 = 2.0_dp * rand - 1.0_dp !nz
!                       W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
!                       SRW02 = sqrt(W02)
!                       rand = sprng(stream(id))
!                       ARGMT = PI * (2.0_dp * rand - 1.0_dp)
!                       U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
!                       V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)

!                       call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, 1, labs)

!                       if (lelectron_scattering.and..not.lfixed_j) then
!                          Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,1,id) / n_rayons
!                       endif

!                       !for one ray
!                       if (.not.lforce_lte) then
!                          call cross_coupling_terms(id, icell, 1)
!                          call calc_rates_mali(id, icell, 1, 1.0_dp / real(n_rayons,kind=dp))
!                       endif

!                       if (loutput_Rates) then
!                          !need to take into account the fact that for MALI no quandities are store for all ray so Rij needs to be computed ray by ray
!                          call store_radiative_rates_mali(id, icell,(iray==1), 1.0_dp / real(n_rayons,kind=dp), &
!                            Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell))
!                       endif


!                    enddo !iray


!                 else if (etape==3) then	!accumulate rays, n_rayons not fixed

!                    do iray=iray_start, iray_start-1+n_rayons !for the new added rays


!                       rand  = sprng(stream(id))
!                       rand2 = sprng(stream(id))
!                       rand3 = sprng(stream(id))

!                       call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

!                       ! Direction de propagation aleatoire
!                       rand = sprng(stream(id))
!                       W0 = 2.0_dp * rand - 1.0_dp !nz
!                       W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
!                       SRW02 = sqrt(W02)
!                       rand = sprng(stream(id))
!                       ARGMT = PI * (2.0_dp * rand - 1.0_dp)
!                       U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
!                       V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)


!                       call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, iray, labs)

!                    enddo !iray

!                    !for all rays accumulated
!                    do iray=1, n_rayons !for all
!                       if (lelectron_scattering.and..not.lfixed_j) then
!                          Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,iray,id) / n_rayons
!                       endif

!                       !for one ray
!                       if (.not.lforce_lte) then
!                          call cross_coupling_terms(id, icell, iray)
!                          call calc_rates_mali(id, icell, iray, 1.0_dp / real(n_rayons,kind=dp))
!                       endif

!                       !might not work cause profile of iray==1 always used
!                       if (loutput_Rates) then
!                          call store_radiative_rates_mali(id, icell, (iray==1), 1.0_dp / &
!                               real(n_rayons,kind=dp), Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell))
!                       endif
!                    enddo

!                 elseif (etape==1) then !ray-by-ray, n_rayons fixed

!                    if (lvoronoi) then
!                    	x0 = Voronoi(icell)%xyz(1)
!                    	y0 = Voronoi(icell)%xyz(2)
!                    	z0 = Voronoi(icell)%xyz(3)
!                    else
!                    	x0 = r_grid(icell)*cos(phi_grid(icell))
!                    	y0 = r_grid(icell)*sin(phi_grid(icell))
!                    	z0 = z_grid(icell)
!                    endif

!                    do imu=1, size(xmu)

!                       w0 = xmu(imu)
!                       u0 = xmux(imu)
!                       v0 = xmuy(imu)

!                       weight = wmu(imu)

!                       call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, 1, labs)

!                       if (lelectron_scattering.and..not.lfixed_j) then
!                          Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,1,id) * weight
!                          !!psi_mean(:,icell) = psi_mean(:,icell) + chi_loc(:,1,id) * psi(:,1,id) * weight
!                       endif

!                       !for one ray
!                       if (.not.lforce_lte) then
!                          call cross_coupling_terms(id, icell, 1)
!                          call calc_rates_mali(id, icell, 1, weight)
!                       endif

!                       if (loutput_Rates) then
!                          !need to take into account the fact that for MALI no quandities are store for all ray so Rij needs to be computed ray by ray
!                          call store_radiative_rates_mali(id, icell, (imu==1), weight, Nmaxtr, &
!                               Rij_all(:,:,icell), Rji_all(:,:,icell))
!                       endif


!                    enddo !imu
!                 end if !etape

!                 call calc_rate_matrix(id, icell, lforce_lte)
! !                 call update_populations(id, icell, diff, .false., n_iter)
!                 call update_populations_and_electrons(id, icell, diff, .false., n_iter, &
!                      (l_iterate_ne.and.icompute_atomRT(icell)==1))!The condition on icompute_atomRT is here because if it is > 0 but /= 1, ne is not iterated!

!                 n_iter_loc = 0
!                 if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc

!                 if (loutput_Rates) then
!                    call store_rate_matrices(id,icell,Nmaxlevel,Gammaij_all(:,:,:,icell))
!                 endif

!              end if !icompute_atomRT
             
!              ! Progress bar
!              !$omp atomic
!              n_cells_done = n_cells_done + 1
!              if (real(n_cells_done) > 0.02*ibar*n_cells) then
!              	call progress_bar(ibar)
!              	!$omp atomic
!              	ibar = ibar+1
!              endif             
             
!           end do !icell
!           !$omp end do
!           !$omp end parallel
!           call progress_bar(50)
!           write(*,*) " " !for progress bar

!           !Ng acceleration
!           !should be para
!           accelerated = .false.
! !           if ( (lNg_acceleration .and. ((n_iter > iNg_Ndelay).and.(maxval(dM) < 1d-2))).and. &
! !                (maxval_cswitch_atoms()==1.0_dp).and.(.not.lprevious_converged) ) then
! !           if ( (lNg_acceleration .and. (n_iter > iNg_Ndelay)).and.(maxval_cswitch_atoms()==1.0_dp)&
! !           	.and.(.not.lprevious_converged) ) then
!           if ( (lNg_acceleration .and. lng_turned_on).and.(maxval_cswitch_atoms()==1.0_dp)&
!           	.and.(.not.lprevious_converged) ) then
!              iorder = n_iter - iNg_Ndelay
!              if (ng_rest) then
!                 write(*,'(" -> Acceleration relaxes... "(1I2)" /"(1I2))') iorder-i0_rest, iNg_Nperiod
!                 if (iorder-i0_rest == iNg_Nperiod) ng_rest = .false.
!              else
!                 !Still, all atoms run at the same speed
!                 i0_rest = iorder
!                 iacc = iacc + 1
!                 write(*,'(" -> Accumulate solutions... "(1I2)" /"(1I2))') iacc, iNg_Norder+2
!                 do nact=1,NactiveAtoms
!                    atom => ActiveAtoms(nact)%ptr_atom
!                    allocate(ng_cur(n_cells * atom%Nlevel)); ng_cur(:) = 0.0_dp
!                    ng_cur = flatten2(atom%Nlevel, n_cells,n_new(nact,1:atom%Nlevel,:))
!                    !has to be parallel in the future

!                    !for many atoms, increment niter only with the first one, as they go at the same speed.
!                    accelerated = ng_accelerate(iacc, ng_cur, n_cells * atom%Nlevel, iNg_Norder, &
!                         ngpop(1:atom%Nlevel*n_cells,:,nact), check_negative_pops=.false.)
                        
!                     !if check_negative_pops is present :
!                     ! if .true., negative pops raise an error.
!                     ! if .false. should I cancel the acceleration ?
!                     !if not present, there is no checking nor handling.

!                    if (accelerated) then
!                       !-> only reshape because we only print accelerated iteration for all atoms at once
!                       n_new(nact, 1:atom%Nlevel,:) = reform2(atom%Nlevel, n_cells, ng_cur)
!                    endif

!                    deallocate(ng_cur)
!                    atom => NULL()
!                 enddo
!                 !for all atoms should be accelerated at the same time
!                 if (accelerated) then
!                    n_iter_accel = n_iter_accel + 1 !True number of accelerated iter
!                    write(*,'("     ++> Ng iteration #"(1I4))') n_iter_accel
!                    ng_rest = (iNg_Nperiod > 0)!.true.
!                    iacc = 0
!                 endif

!              endif
!           endif

!           if ((mod(n_iter, checkpoint_period)==0).and.(lcheckpoint)) then
!              do nact=1, NactiveAtoms
!                 call write_pops_atom(ActiveAtoms(nact)%ptr_atom,iter=n_iter,step=etape)
!              enddo
!           endif

!           !-> evaluate ne for other non-LTE atoms not included in the non-LTE ionisation scheme ?
!           !with fixed non-LTE populations and non-LTE contributions
!           update_bckgr_opac = .false.
!           if (l_iterate_ne) then

!              !update ne from non-LTE ionisation
!              !-> no because dne need tu be zeroed for each new iteration of all ne.
!              ! 					dne = max(dne, maxval(abs(1.0 - ne(:)/(ne_new(:)+1d-50))))
!              dne = maxval(abs(1.0 - ne(:)/(ne_new(:)+1d-50)),mask=icompute_atomRT==1)
!              ! 					dne = 0
!              ! 					do icell=1,n_cells
!              ! 						if (icompute_atomRT(icell)) then
!              ! 							dne = max(dne, abs(1.0 - ne(icell)/ne_new(icell)))
!              ! 						endif
!              ! 					enddo
!              write(*,'("OLD ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') &
!                minval(ne,mask=(icompute_atomRT>0)), maxval(ne)
!              write(*,'("NEW ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') &
!                minval(ne_new,mask=(icompute_atomRT>0)), maxval(ne_new)

!              ne(:) = ne_new(:)
!              !-> For all elements use the old version
! !!! WARNING NOT TESTED YET 05 April 2021 !!!!
!              if ( update_ne_other_nlte ) then
!                 if (helium_is_active .or. hydrogen%active) then
!                    call warning("old ne recipe for non-LTE metals not tested yet!")
!                 endif
!                 call solve_electron_density(ne_start_sol, .true., dne)
!              endif

!              do nact=1, Natom
!                 if (Atoms(nact)%ptr_atom%ID=="H") then
!                    write(*,*) " -> updating LTE pops of hydrogen"
!                    call LTEpops_H
!                    !check pops actually are updated ??
!                 else
!                    write(*,*) " -> updating LTE pops of ",Atoms(nact)%ptr_atom%ID
!                    call LTEpops(Atoms(nact)%ptr_atom,.false.)
!                 endif
!              enddo
!              write(*,*) ''

!              update_bckgr_opac = .true.
!           end if

!           !Global convergence Tests
!           id = 1
!           !should be para
!           dM(:) = 0.0 !all pops
!           diff_cont = 0
!           !add a dpops cont ?
!           diff = 0.0
!           dJ = 0.0
!           lambda_max = 0.0
!           dT = 0.0
!           dTM(:) = 0.0 !all temperature (ion+line)
!           Tex_ref(:) = 0.0
!           Tion_ref(:) = 0.0
!           icell_max = 1
!           icell_max_2 = 1
!           !for all cells
!           dN2 = 0.0 !among all Tex
!           dN4 = 0.0 !among all Tion
!           cell_loop2 : do icell=1,n_cells

!              l_iterate = (icompute_atomRT(icell)>0)

!              if (l_iterate) then

!                 !Local only
!                 dN = 0.0 !for all levels of all atoms of this cell
!                 ! 						dN2 = 0.0 !among all Tex
!                 ! 						dN4 = 0.0 !among all Tion

!                 ! 						write(*,'(1I4, 1F14.4, 2ES17.8E3)') icell, T(icell), ne(icell), nHtot(icell)
!                 do nact=1,NactiveAtoms
!                    atom => ActiveAtoms(nact)%ptr_atom

!                    do ilevel=1,atom%Nlevel
!                       ! 								write(*,'(1I2, " fne="(1L1)," fntot="(1L1))') ilevel, n_new(nact,ilevel,icell) >= frac_ne_limit * ne(icell), n_new(nact,ilevel,icell) >= frac_limit_pops * ntotal_atom(icell, atom)
!                       ! 								write(*,'("--> n="(1ES17.8E3)," fne="(1ES17.8E3)," fntot="(1ES17.8E3))') n_new(nact,ilevel,icell), frac_ne_limit * ne(icell), frac_limit_pops * ntotal_atom(icell, atom)
!                       ! 								write(*,'("--> n/ne="(1ES17.8E3)," n/ntot="(1ES17.8E3))') n_new(nact,ilevel,icell) / ne(icell), n_new(nact,ilevel,icell) / ntotal_atom(icell, atom)
!                       ! 								write(*,*) " "
!                       ! 								if ( n_new(nact,ilevel,icell) >= frac_ne_limit * ne(icell) ) then
!                       if ( n_new(nact,ilevel,icell) >= frac_limit_pops * ntotal_atom(icell, atom) ) then
!                          dN1 = abs(1d0-atom%n(ilevel,icell)/n_new(nact,ilevel,icell))
!                          ! 									dn1 = abs(-2*atom%n(ilevel,icell) + n_new(nact,ilevel,icell) + gpop_old(nact,ilevel,icell)) /  n_new(nact,ilevel,icell)
!                          dN = max(dN1, dN)
!                          dM(nact) = max(dM(nact), dN1)
!                       endif
!                       !convergence_map(icell, ilevel, nact, etape) = dN1
!                    end do !over ilevel
!                    diff_cont = max(diff_cont, abs(1.0 - atom%n(atom%Nlevel,icell)/(1d-50 + n_new(nact,atom%Nlevel,icell) )))
!                    ! 							diff_cont = max(diff_cont, abs(-2*atom%n(atom%Nlevel,icell) + n_new(nact,atom%Nlevel,icell) + gpop_old(nact,atom%Nlevel,icell)) / n_new(nact,atom%Nlevel,icell))

!                    !I keep negative Temperatures for info.debug.
!                    !hence the /= 0.0. But, some work should be done in update_pops and compute_Tex
!                    do kr=1, atom%Nline

!                       if (atom%lines(kr)%Tex(icell) /= 0.0_dp) then
!                          dN3 = abs(1.0 - Tex_old(nact, kr, icell) /  atom%lines(kr)%Tex(icell) )
!                          dN2 = max(dN3,dN2)
!                          dTM(nact) = max(dTM(nact), dN3)
!                          if (dN3 >= dN2) then
!                             Tex_ref(nact) = atom%lines(kr)%Tex(icell)
!                             icell_max = icell
!                          endif
!                       endif

!                    enddo

!                    do kr=1, atom%Ncont

!                       if ( atom%continua(kr)%Tex(icell) /= 0.0_dp) then
!                          dN3 = abs(1.0 - Tex_old(nact, kr+atom%Nline, icell) /  atom%continua(kr)%Tex(icell) )
!                          dN4 = max(dN3,dN4)
!                          dTM(nact) = max(dTM(nact), dN3)
!                          if (dN3 >= dN4) then
!                             Tion_ref(nact) = atom%continua(kr)%Tex(icell)
!                             icell_max_2 = icell
!                          endif
!                       endif

!                    enddo

!                    atom => NULL()
!                 end do !over atoms

!                 !compare for all atoms and all cells
!                 diff = max(diff, dN) ! pops
!                 ! 						diff = max(diff, dN2) ! Texi, only lines


!                 !do not update if lfixed_J
!                 if (lelectron_scattering.and..not.lfixed_j) then
!                    do la=1, Nlambda_cont
!                       dN1 = abs( 1.0_dp - Jnu_cont(la,icell)/Jnew_cont(la,icell) )
!                       ! 								if (Jnu(la,icell) > 0.0_dp) then
!                       ! 									dN1 = abs( (Jnew(la,icell) - Jnu(la,icell))/Jnu(la,icell) )
!                       ! 								endif
!                       if (dN1 > dJ) then
!                          dJ = dN1
!                          lambda_max = lambda_cont(la)
!                       endif
!                       !updating
!                       Jnu_cont(la,icell) = Jnew_cont(la,icell)
!                    enddo

!                 endif

!                 lcell_converged(icell) = (real(diff) < precision)!.and.(dne < precision)

!                 !Re init for next iteration if any
!                 do nact=1, NactiveAtoms
!                    atom => ActiveAtoms(nact)%ptr_atom
!                    ! 							gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
!                    atom%n(:,icell) = n_new(nact,1:atom%Nlevel,icell)
!                    do kr=1,atom%Nline
!                       Tex_old(nact, kr, icell) = atom%lines(kr)%Tex(icell)
!                    enddo
!                    do kr=1, atom%Ncont
!                       Tex_old(nact, kr+Atom%Nline,icell) = atom%continua(kr)%Tex(icell)
!                    enddo
!                    !if (.not.lfix_backgrnd_opac) then
!                    ! !evalute LTE here ?? if so, need a local version of lte pops.
!                    ! !or store the populations nstar in a new array?
!                    ! recompute profiles or damping
!                    !
!                    !endif
!                    atom => NULL()
!                 end do

!                 !update only if ne has been iterated at this iteration
!                 !move elsewhere if electron density is iterated with SEE
!                 if (update_bckgr_opac) then
!                    !safe, used only if not lfixbackground
!                    nHmin(icell) = nH_minus(icell)
!                    ! 							call compute_atom_quantities(icell)
!                    !evaluate profiles anyway ?
!                    !if I recompute all background continua ?  and damping ? evaluate damping locally ?
!                    if (.not.lfix_backgrnd_opac) then!only ne changes
!                       !-> this evaluate profiles or damping but also continuum quantities not needed ???
!                       !-> do not check ni-njgij < 0 because a non converged state can lead to inversion of populations
!                       ! modify it to update only nlte quantities ?
!                       call compute_atom_quantities(icell)
!                       call background_continua_lambda(icell, Nlambda_cont, lambda_cont, chi_c(:,icell), eta_c(:,icell))
!                       !or just re evaluate ne in chi_c like chi_c -ne_old + ne_new
!                    endif
!                 endif

!                 call NLTE_bound_free(icell)

!                 !because of opac nlte changed
!                 if (.not.llimit_mem) then
!                    call interp_continuum_local(icell, chi0_bb(:,icell), eta0_bb(:,icell))
!                 endif
!                 !end init


!              end if !if l_iterate
!           end do cell_loop2 !icell
!           write(*,'("  ---> dnHII="(1ES17.8E3))') diff_cont  
!           !backward derivative
!           if (n_iter > 1) then
!           	conv_speed = (diff - diff_old) !<0 converging.
!           	conv_acc = conv_speed - conv_acc !>0 accelerating
!           endif
!           !!a more dynamic criterion should be use, that also depends on the atom.
!           if (ldamp_jacobi) then
!           	!conv_speed is negative if we converge!
!           	if (conv_speed <= 0) then
!           		if (-conv_speed < 1d-4) then
!           			omega_sor_atom(:) = 1.8
!           		else
!           			omega_sor_atom(:) = 1.0
!           		endif
!           	else
!           		if (conv_speed > 5) then
!           			omega_sor_atom(:) = 0.8
!           		else
!           			omega_sor_atom(:) = 1.0
!           		endif
!           	endif
!           endif
          
!           if (lng_acceleration) then
!           	!be sure we are converging before extrapolating
!           	if ( (n_iter>1).and.(conv_speed <= 0).and.(-conv_speed < 1d-3) ) then
!           		if (.not.lng_turned_on) then
!           			lng_turned_on = .true.
!           			iNg_Ndelay = n_iter
!           		endif
!           	endif
!           endif

!           if (maxval(max_n_iter_loc)>0) write(*,'(" -> "(1I10)" sub-iterations")') maxval(max_n_iter_loc)
!           write(*,'(" -> icell_max1 #"(2I6)," icell_max2 #"(2I6))') icell_max, icell_max_2
!           write(*,*) " ------------------------------------------------ "
!           do nact=1,NactiveAtoms
!              write(*,'("             Atom "(1A2))') ActiveAtoms(nact)%ptr_atom%ID
!              if (accelerated) then
!                 write(*,'("   >>> dpop="(1ES17.8E3)" (Accelerated)")') dM(nact)
!              else
!                 write(*,'("   >>> dpop="(1ES17.8E3))') dM(nact)
!              endif
!              write(*,'("   >>>   dT="(1ES17.8E3))') dTM(nact)
!              write(*,'("    --->   dT(line)="(1ES17.8E3), " dT(cont)="(1ES17.8E3))') dN2, dN4
!              write(*,'("    ->> Te(icell_max2)="(1F14.4)" K", " Tion="(1F14.4)" K")') T(icell_max_2), Tion_ref(nact)
!              write(*,'("    ->> Te(icell_max1)="(1F14.4)" K", " Texi="(1F14.4)" K")') T(icell_max), Tex_ref(nact)
!              write(*,*) " ------------------------------------------------ "
!           enddo
!           if (dne /= 0.0_dp) write(*,'("   >>> dne="(1ES17.8E3))') dne
!           if (dJ /= 0.0_dp)  write(*,'("   >>> dJ="(1ES14.5E3)" @"(1F14.4)" nm")') dJ, lambda_max !at the end of the loop over n_cells
!           write(*,'(" <<->> diff="(1ES17.8E3)," old="(1ES17.8E3))') diff, diff_old !at the end of the loop over n_cells
!           write(*,'("   ->> speed="(1ES17.8E3)"; acc="(1ES17.8E3))') conv_speed, conv_acc
!           write(*,"('Unconverged cells #'(1I6), ' fraction :'(1F12.3)' %')") &
!                size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))), &
!                100.*real(size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0)))) / &
!                real(size(pack(icompute_atomRT,mask=icompute_atomRT>0)))
!           write(*,*) " *************************************************************** "
!           diff_old = diff
!           conv_acc = conv_speed
               

!           ! 				write(*,*) " set lprevious_converged to true for test"
!           ! 				lprevious_converged = .true.

!           !(real(dne) < precision).and.
!           if ((real(diff) < precision).and.(maxval_cswitch_atoms() == 1.0_dp)) then
!              if (lprevious_converged) then
!                 lconverged = .true.
!              else
!                 lprevious_converged = .true.
!              endif
!           else !continue to iterate even if n_rayons max is reached ?
!              lprevious_converged = .false.
!              if ((cswitch_enabled).and.(maxval_cswitch_atoms() > 1.0)) then
!                 call adjust_cswitch_atoms
!              endif
!              if (.not.lfixed_rays) then !step 3 only
!                 n_rayons = n_rayons * 2
!                 write(*,*) ' -- Increasing number of rays :', n_rayons
!                 if (n_rayons > n_rayons_max) then
!                    if (n_iter >= maxIter) then
!                       write(*,*) "Warning : not enough rays to converge !!"
!                       lconverged = .true.
!                    end if
!                 end if
!              else !steps < 3
!                 if (n_iter > 1000) then
!                    write(*,*) "Warning : not enough iterations to converge !!"
!                    lconverged = .true.
!                 end if
!              end if
!           end if

!           !***********************************************************!
!           ! ********** timing and checkpointing **********************!

!           call system_clock(time_end,count_rate=time_tick,count_max=time_max)
!           if (time_end < time_begin) then
!              time_nlte=real(time_end + (1.0 * time_max)- time_begin)/real(time_tick)
!           else
!              time_nlte=real(time_end - time_begin)/real(time_tick)
!           endif

!           if (n_iter <= n_iter_counted) then
!              time_iteration = time_iteration + time_nlte / real(n_iter_counted)
!           endif


!           if (lsafe_stop) then

!              if ((time_nlte + time_iteration >=  safe_stop_time).and.(n_iter >= n_iter_counted)) then
!                 lconverged = .true.
!                 lprevious_converged = .true.
!                 call warning("Time limit would be exceeded, leaving...")
!                 write(*,*) " time limit:", mod(safe_stop_time/60.,60.) ," min"
!                 write(*,*) " ~<time> etape:", mod(n_iter * time_iteration/60.,60.), ' <time iter>=', &
!                      mod(time_iteration/60.,60.)," min"
!                 write(*,*) " ~<time> etape (cpu):", mod(n_iter * time_iteration * nb_proc/60.,60.), " min"
!                 write(*,*) ' time =',mod(time_nlte/60.,60.), " min"
!                 lexit_after_nonlte_loop = .true.
!                 !lsafe_stop only would leave the code even if the time is below the walltime
!                 exit step_loop
!              endif

!           endif
!           !***********************************************************!

!        end do !while
!        write(*,*) "step: ", etape, "Threshold: ", precision!dpops_max_error
!        !real values are possible, not needed and would increase the amount
!        !of lign of codes.
!        if (n_iter >= n_iter_counted) then
!           write(*,*) " ~<time> etape:", mod(n_iter * time_iteration/60.,60.), ' <time iter>=', mod(time_iteration/60.,60.)," min"
!           write(*,*) " ~<time> etape (cpu):", mod(n_iter * time_iteration * nb_proc/60.,60.), " min"
!        endif

!     end do step_loop

!     ! -------------------------------- CLEANING ------------------------------------------ !

!     if (lNg_acceleration) then
!        if (allocated(ngpop)) deallocate(ngpop)
!        if (allocated(ng_cur)) deallocate(ng_cur)
!     endif
    
!     if (ldamp_jacobi) then
!     	deallocate(omega_sor_atom)
!     endif


!     if (allocated(xyz_pos)) then
!        ! 			open(unit=20,file="xyz_pos.txt",status="unknown")
!        ! 			write(20,*) 1, 1
!        ! 			write(20,*) 1627, cell_map_i(1627), cell_map_j(1627), cell_map_k(1627)
!        ! 			do iray=1,1
!        ! 				write(20,'(*(1E20.7E3))') (xyz_pos(i,iray,1),i=1,3)
!        ! 			enddo
!        ! 			close(20)
!        ! 			open(unit=20,file="uvw_pos.txt",status="unknown")
!        ! 			write(20,*) 1, size(xmu), size(xmu)/8
!        ! 			do imu=1,size(xmu)
!        ! 				write(20,'(*(1E20.7E3))') (uvw_pos(i,imu,1),i=1,3)
!        ! 			enddo
!        ! 			close(20)
!     endif

! !     if (allocated(threeKminusJ)) then
! !        !only at id_ref until fits file
! !        write(*,*) " Writing anisotropy to ascii file..."
! ! 
! !        open(unit=20, file="anis_ascii.txt", status="unknown")
! !        write(20,*) n_cells, 1!Nlambda
! !        do icell=1, n_cells
! !           !do la=1, Nlambda
! !           do la=id_ref,id_ref
! !              write(20,'(1F12.5,5E20.7E3)') lambda(la), threeKminusJ(la,icell), 0.0, 0.0, 0.0, 0.0
! !           enddo
! !        enddo
! !        close(20)
! !        write(*,*) "done"
! !        deallocate(threeKminusJ)
! !     endif
!     if (allocated(psi_mean)) then
!        write(*,*) " Writing <Psi> to ascii file..."
!        open(unit=20, file="psim_ascii.txt", status="unknown")
!        write(20,*) n_cells, Nlambda
!        do icell=1, n_cells
!           do la=1, Nlambda
!              write(20,'(1F12.5,5E20.7E3)') lambda(la), psi_mean(la,icell), 0.0, 0.0, 0.0, 0.0
!           enddo
!        enddo
!        close(20)
!        write(*,*) "done"
!     endif

!     if (allocated(convergence_map)) then
!        do nact=1, NactiveAtoms
!           call write_convergence_map_atom(ActiveAtoms(nact)%ptr_atom, etape_end, &
!                convergence_map(:,1:ActiveAtoms(nact)%ptr_atom%Nlevel, nact, :))
!        enddo
!        if (l_iterate_ne) call write_convergence_map_electron(etape_end, convergence_map(:,1,NactiveAtoms+1,:))
!        deallocate(convergence_map)
!     endif

!     deallocate(dM, dTM, Tex_ref, Tion_ref)
!     if (allocated(Jnew_cont)) deallocate(Jnew_cont)
!     deallocate(psi, chi_up, chi_down, Uji_down, eta_atoms, n_new, ne_new)
!     deallocate(stream)
!     if (allocated(xmu)) then
!       deallocate(xmu,wmu,xmux,xmuy)
!     endif

!     !!close(unit=unit_invfile)

!     if (n_iterate_ne > 0) then
!        write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
!        call write_electron
!        call write_Hminus
!     endif



!       do m=1,N_atoms
!          call write_pops_atom(Atoms(m)%p)
!       end do

      !deallocate non-LTE variables and leave.
      !deallocate Itot. Will be allocated again for images.
      deallocate(ne_new,psi,eta_atoms,chi_up,chi_down,uji_down,lcell_converged,itot,phi_loc)!n_new)
      if (allocated(ngpop)) deallocate(ngpop)

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
      real(kind=dp) :: v_char, max_vel_shift, v1, v2

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

         stop
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
      real(kind=dp), parameter :: precision = 1.e-2
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