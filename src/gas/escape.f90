

module escape

    use parametres
    use grid
    use utils, only : progress_bar, Bpnu
    use molecular_emission, only : v_proj, ds
    use naleat, only  : seed, stream, gtype
    use mcfost_env, only  : time_begin, time_end, time_tick, time_max
    use stars, only : is_inshock, intersect_stars, star_rad, T_preshock
    use mem, only : emissivite_dust
    use dust_prop, only : kappa_abs_lte, kappa, kappa_factor
    use opacity_atom, only : contopac_atom_loc, vlabs, Itot, calc_contopac_loc, psi, cross_coupling_cont
    use atom_type, only : n_atoms, Atoms, ActiveAtoms, atomtype, NactiveAtoms, atomPointerArray, &
        lcswitch_enabled, vbroad, maxval_cswitch_atoms, NpassiveAtoms, PassiveAtoms, adjust_cswitch_atoms
    use see, only : alloc_nlte_var, neq_ng, ngpop, small_nlte_fraction, tab_Aji_cont, tab_Vij_cont, &
        dealloc_nlte_var, update_populations, init_colrates_atom
    use wavelengths, only : n_lambda
    use wavelengths_gas, only : tab_lambda_nm, tab_lambda_cont, n_lambda_cont, Nlambda_max_trans, Nlambda_max_line
    use lte, only : ltepops_atoms, LTEpops_atom_loc, LTEpops_H_loc, nH_minus
    use broad, only : line_damping
    use collision_atom, only : collision_rates_atom_loc, collision_rates_hydrogen_loc
    use voigts, only : voigt
    use gas_contopac, only : background_continua_lambda

   !$ use omp_lib

   implicit none

#include "sprng_f.h"

    type (atomPointerArray), dimension(:), allocatable :: sob_atoms
    integer :: Nsob_atoms
    integer, allocatable, dimension(:) :: sob_atom_to_activeatoms

    real(kind=dp), allocatable, dimension(:,:) :: d_to_star, wdi, domega_star, domega_core, domega_shock, dv_proj
    real(kind=dp), allocatable, dimension(:) :: mean_grad_v, mean_length_scale, Tchoc_average
    real(kind=dp), allocatable :: I0_line(:,:,:,:), Jext(:,:), Istar(:,:), Ishock(:,:)

    contains

    subroutine init_rates_escape(id,icell)
        integer, intent(in) :: id, icell
        integer :: n

        do n=1,NactiveAtoms
            call init_radrates_escape_atom(id,icell,ActiveAtoms(n)%p)

            !x ne included. Derivatives to ne not included.
            if (activeatoms(n)%p%id=='H') then
                ! call init_colrates_atom(id,ActiveAtoms(n)%p)
                call collision_rates_hydrogen_loc(id,icell)
            else
                call init_colrates_atom(id,ActiveAtoms(n)%p)
                call collision_rates_atom_loc(id,icell,ActiveAtoms(n)%p)
            endif
        enddo

        return
    end subroutine init_rates_escape

    subroutine init_radrates_escape_atom(id,icell,atom)
        integer, intent(in) :: id, icell
        type (atomtype), intent(inout) :: atom
        integer :: kr, i, j
        !here nlte_fact is .true. if lforce_lte, we don't enter in the Sobolev routine.

        do kr=1,atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            atom%continua(kr)%Rij(id) = 0.0_dp
            !updated value of ni and nj!
            !-> 0 if cont does not contribute to opac.
            ! atom%continua(kr)%Rji(id) = nlte_fact * tab_Aji_cont(kr,atom%activeindex,icell) * &
            !      atom%nstar(i,icell)/atom%nstar(j,icell)
            atom%continua(kr)%Rji(id) = tab_Aji_cont(kr,atom%activeindex,icell) * &
                 atom%ni_on_nj_star(i,icell)
            !check ratio ni_on_nj_star for the continua when there are multiple ionisation states (works for H, test for He)
        enddo

        do kr=1,atom%Nline
            atom%lines(kr)%Rij(id) = 0.0_dp
            atom%lines(kr)%Rji(id) = 0.0_dp !no init at Aji
        enddo

        return

    end subroutine init_radrates_escape_atom

    subroutine alloc_escape_variables()

        allocate(d_to_star(n_cells,n_etoiles),wdi(n_cells,n_etoiles),domega_core(n_cells,n_etoiles))
        !domega_core does not discriminates between the "shock" and the unperturbed stellar surface.
        !If the accretion shock is there, we split and uses 2 different quantities.
        if (laccretion_shock) then
            allocate(domega_star(n_cells,n_etoiles),domega_shock(n_cells,n_etoiles))
            allocate(Tchoc_average(n_etoiles))
        endif
        allocate(mean_grad_v(n_cells), mean_length_scale(n_cells))

        return
    end subroutine alloc_escape_variables
    subroutine dealloc_escape_variables()

        if (allocated(d_to_star)) deallocate(d_to_star,wdi,domega_core)
        if (allocated(domega_star)) deallocate(domega_star)
        if (laccretion_shock) then
            deallocate(Tchoc_average,domega_shock)
        endif
        if (allocated(mean_grad_v)) deallocate(mean_grad_v, mean_length_scale)

        return
    end subroutine dealloc_escape_variables

    subroutine mean_velocity_gradient()
    !TO DO: include dust in the dOmega and also background opacities
        integer :: icell, id, i, previous_cell
        real(kind=dp) :: x0,y0,z0,x1,y1,z1,u,v,w
        real(kind=dp) :: xa,xb,xc,xa1,xb1,xc1,l1,l2,l3
        integer :: next_cell, iray, icell_in
        integer, parameter :: n_rayons = 1000
        real :: rand, rand2, rand3
        real(kind=dp) :: W02,SRW02,ARGMT,v0,v1, r0, wei, F1, T1
        integer :: n_rays_shock(n_etoiles)
        real(kind=dp) :: l,l_contrib, l_void_before, Tchoc, rho_shock(n_etoiles)
        integer :: ibar, n_cells_done
        integer :: i_star, icell_star
        logical :: lintersect_stars, lintersect

        write(*,*) " *** computing mean velocity gradient for each cell.."
        !but also mean solid angle subtended by each cell to the different stars, including
        !shock regions if any. These are used to weight the contributions from the star in the escape prob.
        !SEE.
        !TO DO: occulation by dust should be included somehow.

        ibar = 0
        n_cells_done = 0

        mean_grad_v = tiny_dp!0.0 !un-signed v gradient of the projected velocity.
        mean_length_scale = 0.0
        domega_core = 0.0
        d_to_star = 0.0
        wdi = 0.0

        if (laccretion_shock) then
            domega_shock = 0.0; Tchoc_average = 0.0; domega_star = 0.0
            n_rays_shock(:) = 0
            rho_shock(:) = 0.0
        endif

        !uses only Monte Carlo to estimate these quantities
        id = 1
        if (.not.allocated(stream)) allocate(stream(nb_proc))
        stream = 0.0
        stream(:) = [(init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT),i=1,nb_proc)]

        call progress_bar(0)
        !$omp parallel &
        !$omp default(none) &
        !$omp private(id,icell,iray,rand,rand2,rand3,x0,y0,z0,x1,y1,z1,u,v,w) &
        !$omp private(wei,i_star,icell_star,lintersect_stars,v0,v1,r0)&
        !$omp private(l_contrib,l_void_before,l,W02,SRW02,ARGMT,previous_cell,next_cell) &
        !$omp private(l1,l2,l3,xa,xb,xc,xa1,xb1,xc1,icell_in,Tchoc,F1,T1,lintersect) &
        !$omp shared(Wdi,d_to_star, dOmega_core,etoile,Tchoc_average,rho_shock,nHtot)&
        !$omp shared(phi_grid,r_grid,z_grid,pos_em_cellule,ibar, n_cells_done,stream,n_cells)&
        !$omp shared (mean_grad_v,mean_length_scale,icompute_atomRT,n_etoiles)&
        !$omp shared(laccretion_shock,domega_shock,domega_star,n_rays_shock)
        !$omp do schedule(static,1)
        do icell=1, n_cells
            !$ id = omp_get_thread_num() + 1
            if (icompute_atomRT(icell) <= 0) cycle  !non-empty cells or dusty regions (with no atomic gas)
                                                    ! are not counted yet.
                                                    ! Simply change the condition to < 0 or set the regions to 1.

            wei = 1.0/real(n_rayons)
            rays_loop : do iray=1, n_rayons

                rand  = sprng(stream(id))
                rand2 = sprng(stream(id))
                rand3 = sprng(stream(id))


                call  pos_em_cellule(icell,rand,rand2,rand3,x0,y0,z0)
                do i_star = 1, n_etoiles
                    r0 = sqrt((x0-etoile(i_star)%x)**2 + (y0-etoile(i_star)%y)**2 + (z0-etoile(i_star)%z)**2)
                    !distance and solid-angle to the star assuming that all cells can see a star.
                    !it is an average distance from random points in the cell.
                    d_to_star(icell,i_star) = d_to_star(icell,i_star) + wei * r0
                    ! !it is like domega_core, but here I'm assuming all rays "would" hit the star.
                    !  wdi(icell,i_star) = wdi(icell,i_star) + wei * 0.5*(1.0 - sqrt(1.0 - (etoile(i_star)%r/r0)**2))
                enddo

                ! Direction de propagation aleatoire
                rand = sprng(stream(id))
                W = 2.0_dp * rand - 1.0_dp !nz
                W02 =  1.0_dp - W*W !1-mu**2 = sin(theta)**2
                SRW02 = sqrt(W02)
                rand = sprng(stream(id))
                ARGMT = PI * (2.0_dp * rand - 1.0_dp)
                U = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
                V = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)

                call intersect_stars(x0,y0,z0, u,v,w, lintersect_stars, i_star, icell_star)!will intersect

                !to do propagate along the ray to get dv_max    

                previous_cell = 0 ! unused, just for Voronoi
                call cross_cell(x0,y0,z0, u,v,w,icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)

                if (test_exit_grid(next_cell, x0, y0, z0)) cycle rays_loop

                v0 = v_proj(icell,x0,y0,z0,u,v,w)
                v1 = v_proj(icell,x1,y1,z1,u,v,w)

                if (lintersect_stars) then !"will interesct"
                    dOmega_core(icell,i_star) = dOmega_core(icell,i_star) + wei
                    icell_in = icell
                    !can I move directly at the stellar surface ?
                    ! call move_to_grid(id,x0,y0,z0,u,v,w,next_cell,lintersect)
                    if (laccretion_shock) then !will intersect the shock?
                        !propagate until we reach the stellar surface
                        xa1 = x0; xb1 = y0; xc1 = z0
                        xa = x0; xb = y0; xc = z0
                        !can I just move to the star ?
                        inf : do
                            if (next_cell==icell_star) then
                                !in shock ??
                                if (is_inshock(id, iray, i_star, icell_in, xa, xb, xc, Tchoc, F1, T1)) then 
                                    dOmega_shock(icell,i_star) = dOmega_shock(icell,i_star) + wei
                                    Tchoc_average(i_star) = Tchoc_average(i_star) + Tchoc * nHtot(icell_in)
                                    n_rays_shock(i_star) = n_rays_shock(i_star) + 1
                                    rho_shock(i_star) = rho_shock(i_star) + nHtot(icell_in)
                                else
                                    dOmega_star(icell,i_star) = domega_star(icell,i_star) + wei
                                endif
                                exit inf !because we touch the star
                            endif
                            icell_in = next_cell
                            if (test_exit_grid(icell_in, xa, xb, xc)) exit inf
                            if (icell_in <= n_cells) then
                                if (icompute_atomRT(icell_in) == -2) exit inf
                            endif
                            call cross_cell(xa,xb,xc,u,v,w,icell_in,previous_cell,xa1,xb1,xc1, next_cell,l1,l2,l3)
                            xa = xa1; xb = xb1; xc = xc1
                        enddo inf
                    endif!shock surface
                    !thus, dOmega_shock = f_shock * dOmega_core
                    !dOmega_core  *= (1.0 - f_shock)
                end if !touch star, computes dOmega_core

                !average for all rays
                mean_grad_v(icell) = mean_grad_v(icell) + wei * abs(v0-v1)/(l_contrib * AU_to_m)
                mean_length_scale(icell) = mean_length_scale(icell) + wei * l_contrib*AU_to_m
            enddo rays_loop !rays

            ! Progress bar
            !$omp atomic
            n_cells_done = n_cells_done + 1
            if (real(n_cells_done) > 0.02*ibar*n_cells) then
                call progress_bar(ibar)
                !$omp atomic
                ibar = ibar+1
            endif             
        enddo
        !$omp end do
        !$omp end parallel
        call progress_bar(50)
        if (laccretion_shock) then
            ! Tchoc_average = Tchoc_average / real(n_rays_shock)
            Tchoc_average = Tchoc_average / rho_shock
        endif
      
        write(*,'("max(<gradv>)="(1ES17.8E3)" s^-1; min(<gradv>)="(1ES17.8E3)" s^-1")') maxval(mean_grad_v), minval(mean_grad_v,icompute_atomRT>0)
        write(*,'("max(<l>)="(1ES17.8E3)" m; min(<l>)="(1ES17.8E3)" m")') maxval(mean_length_scale), minval(mean_length_scale,icompute_atomRT>0)
        do i_star=1, n_etoiles
            write(*,*) "star #", i_star
            !it is like domega_core, but here I'm assuming all rays "would" hit the star.
            !it uses the average distance to the star as reference point in the cell.
            where (d_to_star(:,i_star) > 0) wdi(:,i_star) = 0.5*(1.0 - sqrt(1.0 - (etoile(i_star)%r/d_to_star(:,i_star))**2))
            write(*,'("  -- max(<d>)="(1ES17.8E3)"; min(<d>)="(1ES17.8E3))') maxval(d_to_star(:,i_star))/etoile(i_star)%r, minval(d_to_star(:,i_star),icompute_atomRT>0)/etoile(i_star)%r
            write(*,'("  -- max(W)="(1ES17.8E3)"; min(W)="(1ES17.8E3))') maxval(Wdi(:,i_star)), minval(Wdi(:,i_star),icompute_atomRT>0)
            write(*,'("  -- max(dOmegac)="(1ES17.8E3)"; min(dOmegac)="(1ES17.8E3))') maxval(domega_core(:,i_star)), minval(domega_core(:,i_star),icompute_atomRT>0)
            if (laccretion_shock) then
                write(*,'("  -- max(dOmega_shock)="(1ES17.8E3)"; min(dOmega_shock)="(1ES17.8E3))') maxval(domega_shock(:,i_star)), minval(domega_shock(:,i_star),icompute_atomRT>0)
                write(*,'("  -- max(dOmega*)="(1ES17.8E3)"; min(dOmega*)="(1ES17.8E3))') maxval(domega_star(:,i_star)), minval(domega_star(:,i_star),icompute_atomRT>0)
                write(*,*) "  -- <Tshock> = ", Tchoc_average(i_star), ' K'
            endif
        enddo
        if (allocated(stream)) deallocate(stream)
        return
    end subroutine mean_velocity_gradient

    subroutine nlte_loop_sobolev()
    !Large velocity gradient / supersonic approximation
    !   - local Sobolev with no background continua for lines
    !   - optically thin excitation with no lines for continua
    !   - MC rays (no healpix)
        integer :: iray, n_iter, id, i, alloc_status
        integer, parameter :: n_rayons = 100
        integer :: nact, imax, icell_max, icell_max_2
        integer :: icell, ilevel, nb, nr, unconverged_cells
        integer, parameter :: maxIter = 80
        !ray-by-ray integration of the SEE
        integer, parameter :: one_ray = 1
        logical :: lfixed_Rays, lconverged, lprevious_converged
        real :: rand, rand2, rand3, unconverged_fraction
        real(kind=dp) :: precision, vth
        real(kind=dp), dimension(:), allocatable :: wmu, diff_loc, dne_loc
        real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, weight
        real(kind=dp) :: diff, diff_old, dne, dN, dN1, dNc
        real(kind=dp), allocatable :: dM(:)
        integer :: NmaxLine, ng_index
        logical, allocatable :: lcell_converged(:)

        logical :: labs
        logical :: l_iterate, l_iterate_ne
      
        !convergence check
        type (AtomType), pointer :: at
        integer(kind=8) :: mem_alloc_local
        real(kind=dp) :: diff_cont, conv_speed, conv_acc

        !timing and checkpointing
        ! NOTE: cpu time does not take multiprocessing (= nb_proc x the real exec time.)
        !-> overall the non-LTE loop
        real :: time_nlte, time_nlte_loop, time_nlte_cpu
        real :: cpu_time_begin, cpu_time_end, time_nlte_loop_cpu
        integer :: count_start, count_end, itime
        !-> for a single iteration
        integer :: cstart_iter, cend_iter
        real :: time_iteration, cpustart_iter, cpuend_iter, time_iter_avg
        integer :: ibar, n_cells_done

        ! -------------------------------- INITIALIZATION -------------------------------- !
        write(*,*) '---------------------- SOBOLEV NON-LTE LOOP -------------------------- '
        time_nlte_loop = 0!total time in the non-LTE loop for all steps
        time_nlte_loop_cpu = 0
        labs = .true.
        lfixed_rays = .true.
        id = 1

        call alloc_escape_variables
        call mean_velocity_gradient

!TO DO:
        !Only loop on atoms that have initial solution Sobolev / escape
        Nsob_atoms = 0
        do i=1, NactiveAtoms
            if (ActiveAtoms(i)%p%initial == 4)  Nsob_atoms =  Nsob_atoms + 1
        enddo
        allocate(sob_atoms(Nsob_atoms),sob_atom_to_activeatoms(Nsob_atoms))
        i = 0
        do nact=1, NactiveAtoms
            if (ActiveAtoms(nact)%p%initial == 4)  then
                i = i + 1
                sob_atoms(i)%p => ActiveAtoms(nact)%p
                sob_atom_to_activeatoms(i) = nact
            endif
        enddo

        !Use non-LTE pops for electrons and background opacities.
        NmaxLine = 0
        do nact=1,NactiveAtoms
            if (.not.activeatoms(nact)%p%nltepops) activeatoms(nact)%p%nltepops = .true.
            Nmaxline = max(NmaxLine,activeatoms(nact)%p%nline)
        enddo
!TO DO: use n_lambda_cont instead of full lambda.
        !Radiation field
        allocate(I0_line(Nlambda_max_line,NmaxLine,NactiveAtoms,nb_proc)); I0_line = 0.0_dp
        allocate(Jext(n_lambda, nb_proc)) !or n_lambda_cont if I only uses cont wavelengths
        allocate(Istar(n_lambda, n_etoiles))
        !see stars/star_rad()
        do i=1, n_etoiles
            Istar(:,i) = Bpnu(n_lambda, tab_lambda_nm, etoile(i)%T*1d0)
            !TO DO: FUV flux here. The coronal flux is added on top of the stellar photosphere radiation, excluding the shock region.
        enddo
        if (laccretion_shock) then
            allocate(Ishock(n_lambda, n_etoiles))
            do i=1, n_etoiles
                Ishock(:,i) = Bpnu(n_lambda, tab_lambda_nm, Tchoc_average(i))
                !par of the spectrum is due to the pre-shock region.
                if (T_preshock > 0.0 ) then
                    Ishock(:,i) = 0.0
                    where (tab_lambda_nm > 364.2096)
                        Ishock(:,i) = Ishock(:,i) + Bpnu(n_lambda,tab_lambda_nm,Tchoc_average(i))
                    elsewhere (tab_lambda_nm <= 364.2096) !-> can be very large between ~ 40 and 91 nm.
                        Ishock(:,i) = Ishock(:,i) + Bpnu(n_lambda,tab_lambda_nm,T_preshock)
                    endwhere
                endif
            enddo
        endif

        allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
        diff_old = 1.0_dp
        dM(:) = 1.0_dp
        !init in case l_iterate_ne  is .false. (does not change a thing if l_iterate_ne)
        dne = 0.0_dp

        !-> negligible
        mem_alloc_local = 0
        mem_alloc_local = mem_alloc_local + sizeof(dM)

        if (allocated(stream)) deallocate(stream)
        allocate(stream(nb_proc),stat=alloc_status)
        if (alloc_status > 0) call error("Allocation error stream")
        if (allocated(ds)) deallocate(ds)
        allocate(ds(one_ray,nb_proc))
        allocate(vlabs(one_ray,nb_proc))
        allocate(dv_proj(one_ray,nb_proc))
        allocate(lcell_converged(n_cells),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error lcell_converged")
        write(*,*) " size lcell_converged:", sizeof(lcell_converged) / 1024./1024./1024.," GB"
        allocate(diff_loc(n_cells), dne_loc(n_cells),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error diff/dne_loc")
        write(*,*) " size diff_loc:", 2*sizeof(diff_loc) / 1024./1024./1024.," GB"
        write(*,*) ""

        !-> negligible
        mem_alloc_local = mem_alloc_local +3*sizeof(ds) + sizeof(stream)

        call alloc_nlte_var(one_ray,mem=mem_alloc_local)

        ! --------------------------- OUTER LOOP ON STEP --------------------------- !
        precision = dpops_max_error

        write(*,*) ""

        allocate(wmu(n_rayons))
        wmu(:) = 1.0_dp / real(n_rayons,kind=dp)
        write(*,*) ""
        call system_clock(count_start,count_rate=time_tick,count_max=time_max)
        call cpu_time(cpu_time_begin)

        lconverged = .false.
        lprevious_converged = lforce_lte
        lcell_converged(:) = .false.
        diff_loc(:) = 1d50
        dne_loc(:) = 0.0_dp
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

!use less in Sobolev ?
            !TO DO:
            !in principle we don't need that if only atoms with Sobolev itterated
            if (lcswitch_enabled) then
               !deactivate
               if (maxval_cswitch_atoms()==1.0_dp) then
                  write(*,*) " ** cswitch off."
                  lcswitch_enabled = .false.
               endif
            endif

            n_iter = n_iter + 1
            ng_index = Neq_ng - mod(n_iter-1,Neq_ng)
            write(*,'(" *** Iteration #"(1I4)"; threshold: "(1ES11.2E3)"; Nrays: "(1I5))') &
                     n_iter, precision, n_rayons
            ibar = 0
            n_cells_done = 0

            l_iterate_ne = .false.
            if( n_iterate_ne > 0 ) then
               l_iterate_ne = ( mod(n_iter,n_iterate_ne)==0 ) .and. (n_iter>Ndelay_iterate_ne)
               ! if (lforce_lte) l_iterate_ne = .false.
            endif

            call progress_bar(0)
            !$omp parallel &
            !$omp default(none) &
            !$omp private(id,icell,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02,argmt)&
            !$omp private(l_iterate,weight,diff)&
            !$omp shared(lforce_lte,n_cells,voronoi,r_grid,z_grid,phi_grid,Jext,wmu) &
            !$omp shared(pos_em_cellule,labs,n_lambda,tab_lambda_nm, icompute_atomRT,lcell_converged,diff_loc,seed,nb_proc,gtype) &
            !$omp shared(stream,n_rayons_mc,lvoronoi,ibar,n_cells_done,l_iterate_ne,Itot,precision,lcswitch_enabled)
            !$omp do schedule(static,1)
            do icell=1, n_cells
                !$ id = omp_get_thread_num() + 1
                stream(id) = init_sprng(gtype, id-1,nb_proc,seed,SPRNG_DEFAULT)
                if (icompute_atomRT(icell)<=0) then
                    !$omp atomic
                    n_cells_done = n_cells_done + 1
                    if (real(n_cells_done) > 0.02*ibar*n_cells) then
                        call progress_bar(ibar)
             	        !$omp atomic
             	        ibar = ibar+1
                    endif
                    cycle
                endif
                ! if( (diff_loc(icell) < 5d-2 * precision).and..not.lcswitch_enabled ) cycle

                !Init upward radiative rates to 0 and downward radiative rates to 0 or "Aji_cont" for Sobolev.
                !Init collisional rates
                call init_rates_escape(id,icell)
                ! Jext(:,id) = 0.0
                ! ! Position aleatoire dans la cellule
                ! do iray=1,n_rayons

                !     rand  = sprng(stream(id))
                !     rand2 = sprng(stream(id))
                !     rand3 = sprng(stream(id))

                !     call pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                !     ! Direction de propagation aleatoire
                !     rand = sprng(stream(id))
                !     W0 = 2.0_dp * rand - 1.0_dp !nz
                !     W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
                !     SRW02 = sqrt(W02)
                !     rand = sprng(stream(id))
                !     ARGMT = PI * (2.0_dp * rand - 1.0_dp)
                !     U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
                !     V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)

                !     weight = wmu(iray)

                !     !TO DO: use only continuum wavelength in this one to speed up calc
                !     !and assumes no feed back of lines on continua and spectrally flat core radiation for lines.
                !     Itot(:,1,id) = I_sobolev_1ray(id,icell,x0,y0,z0,u0,v0,w0,1,n_lambda,tab_lambda_nm)
                !     Jext(:,id) = Jext(:,id) + weight * Itot(:,1,id)
                !     !using Itot iray too
                !     call partial_coupling_cont(id,icell,weight)
                !     ! call accumulate_radrates_sobolev_1ray(id,icell,1,weight)
                ! enddo !iray
                call radrates_sobolev_average(id, icell)

                !update populations of all atoms including electrons
                call update_populations(id, icell, (l_iterate_ne.and.icompute_atomRT(icell)==1), diff)

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
            ! *******************  update ne metals  *******************!
            if (l_iterate_ne) then
               id = 1
               dne = 0.0_dp !=max(dne_loc)
               dne_loc(:) = abs(1.0_dp - ne(:)/(1d-100 + ngpop(1,NactiveAtoms+1,:,1)))
               !$omp parallel &
               !$omp default(none) &
               !$omp private(id,icell,l_iterate,nact,ilevel,vth,nb,nr)&
               !$omp shared(T,vturb,nHmin,icompute_atomRT,ne,n_cells,ngpop,ng_index)&
               !$omp shared(tab_lambda_nm,atoms,n_atoms,dne,PassiveAtoms,NactiveAtoms,NpassiveAtoms,Voigt)
               !$omp do schedule(dynamic,1)
               do icell=1,n_cells
                  !$ id = omp_get_thread_num() + 1
                  l_iterate = (icompute_atomRT(icell)==1)
                  if (l_iterate) then
                     dne = max(dne, abs(1.0_dp - ne(icell)/ngpop(1,NactiveAtoms+1,icell,1)))
                     ! -> needed for lte pops !
                    ne(icell) = ngpop(1,NactiveAtoms+1,icell,1)
                    ngpop(1,NactiveAtoms+1,icell,ng_index) = ne(icell)
                     call LTEpops_H_loc(icell)
                     nHmin(icell) = nH_minus(icell)
                     do nact = 2, n_atoms
                        call LTEpops_atom_loc(icell,Atoms(nact)%p,.false.)
                     enddo
!use less in Sobolev ?
                     !update profiles only for passive atoms. For active atoms we need new non-LTE pops.
                     !If LTE pops are not updated (so no ne_new) profiles of LTE elements are unchanged.
                     do nact = 1, NpassiveAtoms
                        do ilevel=1,PassiveAtoms(nact)%p%nline
                           if (.not.PassiveAtoms(nact)%p%lines(ilevel)%lcontrib) cycle
                           if (PassiveAtoms(nact)%p%lines(ilevel)%Voigt) then
                              nb = PassiveAtoms(nact)%p%lines(ilevel)%nb; nr = PassiveAtoms(nact)%p%lines(ilevel)%nr
                              PassiveAtoms(nact)%p%lines(ilevel)%a(icell) = line_damping(icell,PassiveAtoms(nact)%p%lines(ilevel))
                              !-> do not allocate if thomson and humlicek profiles !
                              !tmp because of vbroad!
                              ! vth = vbroad(T(icell),PassiveAtoms(nact)%p%weight, vturb(icell))
                              ! !-> beware, temporary array here. because line%v(:) is fixed! only vth changes
                              ! PassiveAtoms(nact)%p%lines(ilevel)%phi(:,icell) = &
                              !      Voigt(PassiveAtoms(nact)%p%lines(ilevel)%Nlambda, PassiveAtoms(nact)%p%lines(ilevel)%a(icell), &
                              !      PassiveAtoms(nact)%p%lines(ilevel)%v(:)/vth)&
                              !       / (vth * sqrtpi)
                           endif
                        enddo
                     enddo !over passive atoms
                  end if !if l_iterate
               end do
               !$omp end do
               !$omp end parallel
               if ((dne < 1d-4 * precision).and.(.not.lcswitch_enabled)) then
                  !Do we need to restart it eventually ?
                  write(*,*) " *** dne", dne
                  write(*,*) " *** stopping electronic density convergence at iteration ", n_iter
                  n_iterate_ne = 0
               endif
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
            !$omp shared(ngpop,Neq_ng,ng_index,Activeatoms,lcell_converged,vturb,T,diff_loc)&!diff,diff_cont,dM)&
            !$omp shared(icompute_atomRT,n_cells,precision,NactiveAtoms,nhtot,voigt,dne_loc)&
            !$omp shared(limit_mem,n_lambda,tab_lambda_nm,n_lambda_cont,tab_lambda_cont) &
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
                        if ( ngpop(ilevel,nact,icell,1) >= small_nlte_fraction * at%Abund*nHtot(icell) ) then
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
                     ngpop(1:at%Nlevel,nact,icell,ng_index) = at%n(:,icell)
!use less in Sobolev ?
                     !Recompute damping and profiles once with have set the new non-LTE pops (and new ne) for next ieration.
                     !Only for Active Atoms here. PAssive Atoms are updated only if electronic density is iterated.
                     !Also need to change if profile interp is used or not! (a and phi)
                     do ilevel=1,at%nline
                        if (.not.at%lines(ilevel)%lcontrib) cycle
                        if (at%lines(ilevel)%Voigt) then
                           nb = at%lines(ilevel)%nb; nr = at%lines(ilevel)%nr
                           at%lines(ilevel)%a(icell) = line_damping(icell,at%lines(ilevel))
                           !-> do not allocate if thomson and humlicek profiles !
                           !tmp because of vbroad!
                           ! vth = vbroad(T(icell),at%weight, vturb(icell))
                           ! !-> beware, temporary array here. because line%v(:) is fixed! only vth changes
                           ! at%lines(ilevel)%phi(:,icell) = &
                           !      Voigt(at%lines(ilevel)%Nlambda, at%lines(ilevel)%a(icell), at%lines(ilevel)%v(:)/vth) &
                           !      / (vth * sqrtpi)
                        endif
                     enddo
                  end do
                  at => null()
!use less in Sobolev ?
                  !if electronic density is not updated, it is not necessary
                  !to compute the lte continous opacities.
                  !but the overhead should be negligible at this point.
                  select case (limit_mem)
                     case (0)
                        call calc_contopac_loc(icell,n_lambda,tab_lambda_nm)
                     case (1)
                        call calc_contopac_loc(icell,n_lambda_cont,tab_lambda_cont)
                     !case (2) evaluated on-the-fly.
                  end select

               end if !if l_iterate
            end do cell_loop2 !icell
            !$omp end do
            !$omp end parallel

            if (n_iter > 1) then
               conv_speed = (diff - diff_old) !<0 converging.
               conv_acc = conv_speed - conv_acc !>0 accelerating
            endif

            do nact=1,NactiveAtoms
               write(*,'("                  "(1A2))') ActiveAtoms(nact)%p%ID
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

            if ((diff < precision).and.(.not.lcswitch_enabled))then
               if (lprevious_converged) then
                  lconverged = .true.
               else
                  lprevious_converged = .true.
               endif
            else
               lprevious_converged = .false.
               !TO DO: remove cswitch in Sobolev ??
               if ((lcswitch_enabled).and.(maxval_cswitch_atoms() > 1.0_dp)) then
                  call adjust_cswitch_atoms()
               endif
               if (lfixed_rays) then
                  if (n_iter > maxIter) then
                     call warning("not enough iterations to converge !!")
                     lconverged = .true.
                  end if
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

            !force convergence if there are only few unconverged cells remaining
            ! if ((n_iter > maxIter/4).and.(unconverged_fraction < 5.0)) then
            if ( (n_iter > maxIter/4).and.(unconverged_fraction < 3.0).or.&
               ((unconverged_fraction < 3.0).and.(time_nlte + time_iteration >= 0.5*safe_stop_time)) ) then
               write(*,'("WARNING: there are less than "(1F6.2)" % of unconverged cells after "(1I4)" iterations")') &
                  unconverged_fraction, n_iter
               write(*,*) " -> forcing convergence"
               lconverged = .true.
            endif


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
                  exit
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

        if (allocated(wmu)) deallocate(wmu)
        call dealloc_escape_variables

      ! -------------------------------- CLEANING ------------------------------------------ !

      if (n_iterate_ne > 0) then
        write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
        !  call write_electron
      endif

    !   do nact=1,N_atoms
    !      call write_pops_atom(Atoms(nact)%p)
    !   end do


    !   if (loutput_rates) call write_rates()

      call dealloc_nlte_var()
      deallocate(dM, diff_loc, dne_loc, I0_line, Jext)
      deallocate(dv_proj, stream, ds, vlabs, lcell_converged, Istar)
      if (allocated(Ishock)) deallocate(Ishock)
      deallocate(sob_atoms,sob_atom_to_activeatoms)

      ! --------------------------------    END    ------------------------------------------ !
      if (mod(time_nlte_loop/60./60.,60.) > 1.0_dp) then
         write(*,*) ' (non-LTE loop) time =',mod(time_nlte_loop/60./60.,60.), " h"
      else
         write(*,*) ' (non-LTE loop) time =',mod(time_nlte_loop/60.,60.), " min"
      endif
      write(*,*) '-------------------------- ------------ ------------------------------ '

      return
    end subroutine nlte_loop_sobolev
!building
    subroutine partial_coupling_cont(id,icell,domega)
    use opacity_atom, only : eta_atoms, Uji_down, chi_up, chi_down
        integer, intent(in) :: id, icell
        real(kind=dp), intent(in) :: domega
        integer :: ns
        type (AtomType), pointer :: at

        Uji_down(:,:,:,id) = 0.0_dp
        chi_down(:,:,:,id) = 0.0_dp
        chi_up(:,:,:,id)   = 0.0_dp

        eta_atoms(:,:,id) = 0.0_dp

        do ns=1, Nsob_atoms
            at => sob_atoms(ns)%p
            call cross_coupling_cont(id,icell,at) !limit_mem == 0 because we will use the tab_lambda_cont in the future
        enddo

        !accumulate rates
        do ns=1, Nsob_atoms
            at => sob_atoms(ns)%p

        enddo
        at => null()
        return
    end subroutine partial_coupling_cont

    subroutine radrates_sobolev_average(id, icell)
        integer, intent(in) :: id, icell
        real, parameter :: fact_tau = 3.0
        real(kind=dp), parameter :: prec_vel = 1.0 / Rsun ! [s^-1]
        integer :: ns, nact, i, j, kc, kr, n0, nb, nr, Nl, l, i0
        real(kind=dp) :: tau0, beta, chi_ij, Icore, l0, dvds
        type(AtomType), pointer :: at
        real(kind=dp) :: ni_on_nj_star, tau_escape, vth, gij, chi0(1), eta0(1)
        real(kind=dp) :: jbar_down, jbar_up, ehnukt, anu, tau_max
        real(kind=dp), dimension(Nlambda_max_trans) :: Ieff, xl

        dvds = mean_grad_v(icell)

        Itot(:,1,id) = 0.0
        if (laccretion_shock) then
            do i=1, n_etoiles
                Itot(:,1,id) = Itot(:,1,id) + domega_shock(icell,i) * Ishock(:,i) + domega_star(icell,i)*Istar(:,i)
            enddo
        else
            do i=1, n_etoiles
                Itot(:,1,id) = Itot(:,1,id) + domega_core(icell,i) * Istar(:,i)
            enddo
        endif

        !-> no overlapping transitions and background continua
        ! at_loop : do nact = 1, Nactiveatoms
        !     at => ActiveAtoms(nact)%p
        at_loop : do ns = 1, Nsob_atoms
            nact = sob_atom_to_activeatoms(ns)
            at => sob_atoms(ns)%p

            vth = vbroad(T(icell),at%weight, vturb(icell))

            !-> purely local
            atr_loop : do kr = 1, at%nline

                i = at%lines(kr)%i
                j = at%lines(kr)%j
                n0 = (at%lines(kr)%Nr - at%lines(kr)%Nb + 1) / 2

                !pops inversions
                if (at%n(i,icell)-at%lines(kr)%gij*at%n(j,icell) <= 0.0) cycle atr_loop
         
                !assumes the underlying radiation is flat across the line.
                !or take n0 as a function of velocity
                Icore = Itot(n0,1,id)

   
                !could be stored in mem.
                chi_ij = hc_fourPI * at%lines(kr)%Bij * (at%n(i,icell)-at%lines(kr)%gij*at%n(j,icell))
                tau_escape = fact_tau * chi_ij * mean_length_scale(icell) / vth
                if (dvds > prec_vel) then
                    !s^-1
                    tau0 = fact_tau * chi_ij / dvds
                    if (tau0 > tau_escape) tau0 = tau_escape
                else !escape probability
                    tau0 = tau_escape
                endif
                beta = min((1.0 - exp(-tau0))/tau0,1.0_dp)

                at%lines(kr)%Rji(id) = beta * at%lines(kr)%Aji + beta * at%lines(kr)%Bji * Icore
                at%lines(kr)%Rij(id) = beta * at%lines(kr)%Bij * Icore


            enddo atr_loop

            !-> opt thin excitation for the continua
            !-> if the continuum is optically thick, let the rates to 0 (LTE approx.)
            ctr_loop : do kr = 1, at%ncont
                if (.not.at%continua(kr)%lcontrib) cycle ctr_loop

                i = at%continua(kr)%i; j = at%continua(kr)%j

                ni_on_nj_star = at%ni_on_nj_star(i,icell)
                gij = ni_on_nj_star * exp(-hc_k/T(icell)/at%continua(kr)%lambda0)
                if ((at%n(i,icell) - at%n(j,icell) * gij) <= 0.0_dp) cycle ctr_loop

                Nb = at%continua(kr)%Nb; Nr = at%continua(kr)%Nr
                Nl = Nr - Nb + 1
                i0 = Nb - 1

                !not only active continua but passive (but no dust here)
                !-> length x chi at lambda0 (beware if dissociation)
                ! call background_continua_lambda(icell, 1, (/at%continua(kr)%lambda0/), chi0(:), eta0(:))
                tau_max = mean_length_scale(icell) * tab_Vij_cont(Nl,kr,nact) * (at%n(i,icell) - at%n(j,icell) * gij)! + chi0(1))
                if (tau_max > 10.0) then
                    at%continua(kr)%Rij(id) = 0.0_dp; at%continua(kr)%Rji(id) = 0.0_dp
                    cycle ctr_loop !force LTE (only collisions)
                endif
                Ieff(1:Nl) = Itot(nb:nr,1,id) !from the core, opt thin

                Jbar_down = 0.0
                Jbar_up = 0.0

                !TO DO MALI precondition only in the trans
                ! Ieff(1:Nl) = Jext(nb:nr,id)
                ! if (tau_max > 10.0) then
                !     Ieff(1:Nl) = at%n(j,icell) * (twohc / tab_lambda_nm(nb:nr)**3) * ni_on_nj_star * exp(-hc_k/T(icell)/tab_lambda_nm(nb:nr))
                !     Ieff(1:Nl) = Ieff(1:Nl) / (at%n(i,icell) - ni_on_nj_star * exp(-hc_k/T(icell)/tab_lambda_nm(nb:nr)) * at%n(j,icell))
                ! endif
                ! xl(1:Nl) = tab_lambda_nm(Nb:nr)


                ! Jbar_up = 0.5 * sum ( &
                !     (tab_vij_cont(1:Nl-1,kr,nact)*Ieff(1:Nl-1)/xl(1:Nl-1) + tab_vij_cont(2:Nl,kr,nact)*Ieff(2:Nl)/xl(2:Nl)) * &
                !     (xl(2:Nl)-xl(1:Nl-1)) &
                ! )
                ! Jbar_down = 0.5 * sum ( &
                !     (tab_vij_cont(1:Nl-1,kr,nact)*Ieff(1:Nl-1)/xl(1:Nl-1) * exp(-hc_k/T(icell)/xl(1:Nl-1)) + &
                !     tab_vij_cont(2:Nl,kr,nact)*Ieff(2:Nl)/xl(2:Nl) * exp(-hc_k/T(icell)/xl(2:Nl))) * &
                !     (xl(2:Nl)-xl(1:Nl-1)) &
                ! )

                do l=1, Nl
                    if (l==1) then
                        xl(l) = 0.5*(tab_lambda_nm(Nb+1)-tab_lambda_nm(Nb)) / tab_lambda_nm(Nb)
                    elseif (l==n_lambda) then
                        xl(l) = 0.5*(tab_lambda_nm(Nr)-tab_lambda_nm(Nr-1)) / tab_lambda_nm(Nr)
                    else
                        xl(l) = 0.5*(tab_lambda_nm(i0+l+1)-tab_lambda_nm(i0+l-1)) / tab_lambda_nm(i0+l)
                    endif
                    anu = tab_Vij_cont(l,kr,nact)

                    Jbar_up = Jbar_up + anu*Ieff(l)*xl(l)

                    ehnukt = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l))

                    Jbar_down = jbar_down + anu*Ieff(l)*ehnukt*xl(l)

                enddo

                at%continua(kr)%Rij(id) = fourpi_h * Jbar_up
                !init at tab_Aji_cont(kr,nact,icell) <=> Aji
                at%continua(kr)%Rji(id) = at%continua(kr)%Rji(id) + fourpi_h * Jbar_down * ni_on_nj_star

            enddo ctr_loop

            at => NULL()

        end do at_loop

        return
    end subroutine radrates_sobolev_average

    function I_sobolev_1ray(id,icell_in,x,y,z,u,v,w,iray,N,lambda)
    ! ------------------------------------------------------------------------------- !
    ! ------------------------------------------------------------------------------- !
      integer, intent(in) :: id, icell_in, iray
      real(kind=dp), intent(in) :: u,v,w
      real(kind=dp), intent(in) :: x,y,z
      integer, intent(in) :: N
      real(kind=dp), dimension(N), intent(in) :: lambda
      real(kind=dp), dimension(N) :: I_sobolev_1ray, Icore(N)
      real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, Q, P(4)
      real(kind=dp), dimension(N) :: Snu, tau, dtau, chi, coronal_irrad
      integer :: kr, nat, nl
      integer, target :: icell
      integer, pointer :: p_icell
      integer :: nbr_cell, next_cell, previous_cell, icell_star, i_star, icell_prev
      logical :: lcellule_non_vide, lintersect_stars

      x1=x;y1=y;z1=z
      x0=x;y0=y;z0=z
      next_cell = icell_in
      nbr_cell = 0
      icell_prev = icell_in

      tau(:) = 0.0_dp

      I_sobolev_1ray(:) = 0.0_dp

      p_icell => icell_ref
      if (lvariable_dust) p_icell => icell

      ! Will the ray intersect a star
      call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)
      ! Boucle infinie sur les cellules (we go over the grid.)
      infinie : do ! Boucle infinie
      ! Indice de la cellule
         icell = next_cell
         x0=x1 ; y0=y1 ; z0=z1

         lcellule_non_vide = (icell <= n_cells)

         ! Test sortie ! "The ray has reach the end of the grid"
         if (test_exit_grid(icell, x0, y0, z0)) return

         if (lintersect_stars) then !"will interesct"
            if (icell == icell_star) then!"has intersected"
                Icore(:) = star_rad(id,iray,i_star,icell_prev,x0,y0,z0,u,v,w,N,lambda)
                I_sobolev_1ray = I_sobolev_1ray + exp(-tau(:)) * Icore
                do nat=1, NactiveAtoms
                    do kr=1, ActiveAtoms(nat)%p%nline
                        nl = ActiveAtoms(nat)%p%lines(kr)%Nr - ActiveAtoms(nat)%p%lines(kr)%Nb + 1
                        I0_line(1:nl,kr,nat,id) = Icore(ActiveAtoms(nat)%p%lines(kr)%Nb:ActiveAtoms(nat)%p%lines(kr)%Nr) 
                    enddo
                enddo
               return
            end if
         endif

         !Special handling of coronal irradiation from "above".
         !mainly for 1d stellar atmosphere
         if (lcellule_non_vide) then
            if (icompute_atomRT(icell) == -2) then
               !Does not return but cell is empty (lcellule_non_vide is .false.)
               coronal_irrad = linear_1D_sorted(atmos_1d%Ncorona,atmos_1d%x_coro(:), &
                                                   atmos_1d%I_coro(:,1),N,lambda)
               I_sobolev_1ray = I_sobolev_1ray + exp(-tau) * coronal_irrad
               lcellule_non_vide = .false.
            endif
         endif

         nbr_cell = nbr_cell + 1

         ! Calcul longeur de vol et profondeur optique dans la cellule
         previous_cell = 0 ! unused, just for Voronoi
         call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)

         !count opacity only if the cell is filled, else go to next cell
         if (lcellule_non_vide) then
            chi(:) = 1d-300; Snu(:) = 0.0_dp
            ! opacities in m^-1, l_contrib in au

            if (icompute_atomRT(icell) > 0) then
               !re-init chi
               call contopac_atom_loc(icell, N, lambda, chi, Snu)
            endif
            if (ldust_atom) then
               chi = chi + kappa_abs_LTE(p_icell,:) * kappa_factor(icell) * m_to_AU ! [m^-1]
               ! Snu = Snu + kappa_abs_LTE(p_icell,:) * kappa_factor(icell) * m_to_AU * Bpnu(N,lambda,T(icell)) ! [W m^-3 Hz^-1 sr^-1]
               Snu = Snu + emissivite_dust(:,icell) ! [W m^-3 Hz^-1 sr^-1]
            endif

            dtau(:) = l_contrib * chi(:) * AU_to_m !au * m^-1 * au_to_m

            if (nbr_cell==1) then
               ds(iray,id) = l_contrib * AU_to_m
               psi(:,1,id) = (1_dp - exp(-dtau) ) / chi
               dv_proj(iray,id) = v_proj(icell,x1,y1,z1,u,v,w) - v_proj(icell,x0,y0,z0,y,v,w)
            endif

            I_sobolev_1ray = I_sobolev_1ray + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu / chi
            tau(:) = tau(:) + dtau(:) !for next cell

         end if  ! lcellule_non_vide

         icell_prev = icell
         !duplicate with previous_cell, but this avoid problem with Voronoi grid here

        end do infinie

        return
    end function I_sobolev_1ray

    subroutine accumulate_radrates_sobolev_1ray(id, icell, iray, dOmega)
        integer, intent(in) :: id, icell, iray
        real(kind=dp), intent(in) :: dOmega
        real, parameter :: fact_tau = 3.0
        real(kind=dp), parameter :: prec_vel = 1.0 / Rsun ! [s^-1]
        integer :: ns, nact, i, j, kc, kr, n0, nb, nr, Nl, l, i0
        real(kind=dp) :: tau0, beta, chi_ij, Icore, l0, dvds, tau_max
        type(AtomType), pointer :: at
        real(kind=dp) :: ni_on_nj_star, tau_escape, vth, gij, jbar_down, jbar_up, ehnukt, anu
        real(kind=dp), dimension(Nlambda_max_trans) :: Ieff, xl

        dvds = abs(dv_proj(iray,id)) / ds(iray,id)

        !-> no overlapping transitions and background continua
        ! at_loop : do nact = 1, Nactiveatoms
        !     at => ActiveAtoms(nact)%p
        at_loop : do ns = 1, Nsob_atoms
            nact = sob_atom_to_activeatoms(ns)
            at => sob_atoms(ns)%p

            vth = vbroad(T(icell),at%weight, vturb(icell))

            !-> purely local
            atr_loop : do kr = 1, at%nline

                i = at%lines(kr)%i
                j = at%lines(kr)%j
                n0 = (at%lines(kr)%Nr - at%lines(kr)%Nb + 1) / 2

                !pops inversions
                if (at%n(i,icell)-at%lines(kr)%gij*at%n(j,icell) <= 0.0) cycle atr_loop
         
                !assumes the underlying radiation is flat across the line.
                !or take n0 as a function of velocity
                Icore = maxval(I0_line(:,kr,nact,id)) !shock + star weighted by the solid angle
   
                !could be stored in mem.
                chi_ij = hc_fourPI * at%lines(kr)%Bij * (at%n(i,icell)-at%lines(kr)%gij*at%n(j,icell))
                tau_escape = fact_tau * chi_ij * ds(iray,id) / vth
                if (dvds > prec_vel) then
                    !s^-1
                    tau0 = fact_tau * chi_ij / dvds
                    if (tau0 > tau_escape) tau0 = tau_escape
                else !escape probability
                    tau0 = tau_escape
                endif
                beta = min((1.0 - exp(-tau0))/tau0,1.0_dp)

                at%lines(kr)%Rji(id) = at%lines(kr)%Rji(id) + dOmega * ( beta * at%lines(kr)%Aji + beta * at%lines(kr)%Bji * Icore )
                at%lines(kr)%Rij(id) = at%lines(kr)%Rij(id) + dOmega * beta * at%lines(kr)%Bij * Icore


            enddo atr_loop

            !-> opt thin excitation for the continua
            ctr_loop : do kr = 1, at%ncont
                if (.not.at%continua(kr)%lcontrib) cycle ctr_loop

                i = at%continua(kr)%i; j = at%continua(kr)%j

                ni_on_nj_star = at%ni_on_nj_star(i,icell)
                gij = ni_on_nj_star * exp(-hc_k/T(icell)/at%continua(kr)%lambda0)
                if ((at%n(i,icell) - at%n(j,icell) * gij) <= 0.0_dp) cycle ctr_loop

                Nb = at%continua(kr)%Nb; Nr = at%continua(kr)%Nr
                Nl = Nr - Nb + 1
                i0 = Nb - 1

                !not only active continua but passive (but no dust here)
                ! call contopac_atom_loc(icell, N, lambda, chi0, Snu0)
                tau_max = mean_length_scale(icell) * tab_Vij_cont(Nl,kr,nact) * (at%n(i,icell) - at%n(j,icell) * gij)
                if (tau_max > 10.0) then
                    at%continua(kr)%Rij(id) = 0.0_dp; at%continua(kr)%Rji(id) = 0.0_dp
                    cycle ctr_loop !force LTE (only collisions)
                endif

                Jbar_down = 0.0
                Jbar_up = 0.0

!TO DO: MALI ?? for cont
                Ieff(1:Nl) = Itot(Nb:Nr,1,id)
                ! xl(1:Nl) = tab_lambda_nm(Nb:nr)


                ! Jbar_up = 0.5 * sum ( &
                !     (tab_vij_cont(1:Nl-1,kr,nact)*Ieff(1:Nl-1)/xl(1:Nl-1) + tab_vij_cont(2:Nl,kr,nact)*Ieff(2:Nl)/xl(2:Nl)) * &
                !     (xl(2:Nl)-xl(1:Nl-1)) &
                ! )
                ! Jbar_down = 0.5 * sum ( &
                !     (tab_vij_cont(1:Nl-1,kr,nact)*Ieff(1:Nl-1)/xl(1:Nl-1) * exp(-hc_k/T(icell)/xl(1:Nl-1)) + &
                !     tab_vij_cont(2:Nl,kr,nact)*Ieff(2:Nl)/xl(2:Nl) * exp(-hc_k/T(icell)/xl(2:Nl))) * &
                !     (xl(2:Nl)-xl(1:Nl-1)) &
                ! )

                do l=1, Nl
                    if (l==1) then
                        xl(l) = 0.5*(tab_lambda_nm(Nb+1)-tab_lambda_nm(Nb)) / tab_lambda_nm(Nb)
                    elseif (l==n_lambda) then
                        xl(l) = 0.5*(tab_lambda_nm(Nr)-tab_lambda_nm(Nr-1)) / tab_lambda_nm(Nr)
                    else
                        xl(l) = 0.5*(tab_lambda_nm(i0+l+1)-tab_lambda_nm(i0+l-1)) / tab_lambda_nm(i0+l)
                    endif
                    anu = tab_Vij_cont(l,kr,nact)

                    Jbar_up = Jbar_up + anu*Ieff(l)*xl(l)

                    ehnukt = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l))

                    Jbar_down = jbar_down + anu*Ieff(l)*ehnukt*xl(l)

                enddo

                at%continua(kr)%Rij(id) = at%continua(kr)%Rij(id) + dOmega*fourpi_h * Jbar_up
                at%continua(kr)%Rji(id) = at%continua(kr)%Rji(id) + ni_on_nj_star * dOmega * fourpi_h * Jbar_down
                ! at%continua(kr)%Rji(id) = at%continua(kr)%Rji(id) + ni_on_nj_star * dOmega * (&
                !     tab_Aji_cont(kr,nact,icell) + fourpi_h * Jbar_down)

            enddo ctr_loop

            at => NULL()

        end do at_loop

        return
    end subroutine accumulate_radrates_sobolev_1ray

end module escape