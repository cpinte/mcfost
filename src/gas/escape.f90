

module escape

    use parametres
    use grid
    use utils, only : progress_bar, Bpnu
    use molecular_emission, only : v_proj
    use naleat, only  : seed, stream, gtype
    use mcfost_env, only  : time_begin, time_end, time_tick, time_max
    use stars, only : is_inshock, intersect_stars

   !$ use omp_lib

   implicit none

#include "sprng_f.h"

    real(kind=dp), allocatable, dimension(:,:) :: d_to_star, wdi, domega_star, domega_core, domega_shock
    real(kind=dp), allocatable, dimension(:) :: mean_grad_v, mean_length_scale

    contains

    subroutine alloc_escape_variables()

        allocate(d_to_star(n_cells,n_etoiles),wdi(n_cells,n_etoiles),domega_core(n_cells,n_etoiles))
        !domega_core does not discriminates between the "shock" and the unperturbed stellar surface.
        !If the accretion shock is there, we split and uses 2 different quantities.
        if (laccretion_shock) then
            allocate(domega_star(n_cells,n_etoiles),domega_shock(n_cells,n_etoiles))
        endif
        allocate(mean_grad_v(n_cells), mean_length_scale(n_cells))

        return
    end subroutine alloc_escape_variables
    subroutine dealloc_escape_variables()

        if (allocated(d_to_star)) deallocate(d_to_star,wdi,domega_core)
        if (allocated(domega_star)) then
            deallocate(domega_star,domega_shock)
        endif
        if (allocated(mean_grad_v)) deallocate(mean_grad_v, mean_length_scale)

        return
    end subroutine dealloc_escape_variables

    subroutine mean_velocity_gradient()
        integer :: icell, id, i, previous_cell
        real(kind=dp) :: x0,y0,z0,x1,y1,z1,u,v,w
        real(kind=dp) :: xa,xb,xc,xa1,xb1,xc1,l1,l2,l3
        integer :: next_cell, iray, icell_in
        integer, parameter :: n_rayons = 100
        real :: rand, rand2, rand3
        real(kind=dp) :: W02,SRW02,ARGMT,v0,v1, r0, wei, F1, T1
        integer :: n_rays_shock(n_etoiles)
        real(kind=dp) :: l,l_contrib, l_void_before, Tchoc, Tchoc_average(n_etoiles)
        integer :: ibar, n_cells_done
        integer :: i_star, icell_star
        logical :: lintersect_stars, lintersect

        write(*,*) " *** computing mean velocity gradient for each cell.."

        ibar = 0
        n_cells_done = 0

        mean_grad_v = tiny_dp!0.0
        mean_length_scale = 0.0
        domega_core = 0.0
        d_to_star = 0.0
        wdi = 0.0

        if (laccretion_shock) then
            domega_shock = 0.0; Tchoc_average = 0.0; domega_star = 0.0
            n_rays_shock(:) = 0
        endif

        !uses only Monte Carlo to estimate these quantities
        id = 1
        stream = 0.0
        stream(:) = [(init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT),i=1,nb_proc)]

        call progress_bar(0)
        !$omp parallel &
        !$omp default(none) &
        !$omp private(id,icell,iray,rand,rand2,rand3,x0,y0,z0,x1,y1,z1,u,v,w) &
        !$omp private(wei,i_star,icell_star,lintersect_stars,v0,v1,r0)&
        !$omp private(l_contrib,l_void_before,l,W02,SRW02,ARGMT,previous_cell,next_cell) &
        !$omp private(l1,l2,l3,xa,xb,xc,xa1,xb1,xc1,icell_in,Tchoc,F1,T1,lintersect) &
        !$omp shared(Wdi,d_to_star, dOmega_core,etoile,Tchoc_average)&
        !$omp shared(phi_grid,r_grid,z_grid,pos_em_cellule,ibar, n_cells_done,stream,n_cells)&
        !$omp shared (mean_grad_v,mean_length_scale,icompute_atomRT)&
        !$omp shared(laccretion_shock,domega_shock,domega_star,n_rays_shock)
        !$omp do schedule(static,1)
        do icell=1, n_cells
            !$ id = omp_get_thread_num() + 1
            if (icompute_atomRT(icell) <= 0) cycle  !non-empty cells or dusty regions (with no atomic gas)
                                                    ! are not counted yet.
                                                    ! Simply change the condition to < 0 or set the regions to 1.

            wei = 1.0/real(n_rayons)
            do iray=1, n_rayons

                rand  = sprng(stream(id))
                rand2 = sprng(stream(id))
                rand3 = sprng(stream(id))


                call  pos_em_cellule(icell,rand,rand2,rand3,x0,y0,z0)
                r0 = sqrt(x0**2+y0**2+z0**2)
                !distance and solid-angle to the star assuming that all cells can see a star
                d_to_star(icell,:) = d_to_star(icell,:) + wei * ( r0 - etoile(:)%r )
                wdi(icell,:) = wdi(icell,:) + wei * 0.5*(1.0 - sqrt(1.0 - (etoile(:)%r/r0)**2))

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
                if (test_exit_grid(next_cell, x0, y0, z0)) then
                    v1 = 0.0
                else
                    v1 = v_proj(next_cell,x1,y1,z1,u,v,w)
                endif

                v0 = v_proj(icell,x0,y0,z0,u,v,w)

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
                                    Tchoc_average(i_star) = Tchoc_average(i_star) + wei * Tchoc
                                    n_rays_shock(i_star) = n_rays_shock(i_star) + 1
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
            enddo !rays

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
            Tchoc_average = Tchoc_average / real(n_rays_shock)
        endif
      
    write(*,'("max(<gradv>)="(1ES17.8E3)" s^-1; min(<gradv>)="(1ES17.8E3)" s^-1")') maxval(mean_grad_v), minval(mean_grad_v,icompute_atomRT>0)
    write(*,'("max(dOmegac)="(1ES17.8E3)"; min(dOmegac)="(1ES17.8E3))') maxval(domega_core), minval(domega_core(:,i_star),icompute_atomRT>0)
    if (laccretion_shock) then
        do i=1, i_star
            write(*,*) "star #", i_star
            write(*,'("max(dOmega_shock)="(1ES17.8E3)"; min(dOmega_shock)="(1ES17.8E3))') maxval(domega_shock(:,i_star)), minval(domega_shock(:,i_star),icompute_atomRT>0)
            write(*,'("max(dOmega*)="(1ES17.8E3)"; min(dOmega*)="(1ES17.8E3))') maxval(domega_star(:,i_star)), minval(domega_star(:,i_star),icompute_atomRT>0)
            write(*,*) "<Tshock> = ", Tchoc_average(i_star), ' K'
        enddo
    endif
    write(*,*) "max vshift=", maxval(mean_grad_v * mean_length_scale) * 1d-3
        return
    end subroutine mean_velocity_gradient

    subroutine escape_prob
    !Probabilistic solution of the non-LTE transfer.

    end subroutine escape_prob

    subroutine averaged_sobolev()
    ! A Sobolev method with averaged quantities
    ! dv/ds -> <dv/ds>

        return
    end subroutine averaged_sobolev

    subroutine sobolev()
    !Large velocity gradient / supersonic approximation
    !   - local Sobolev with no background continua for lines
    !   - optically thin excitation with no lines for continua


        return
    end subroutine sobolev


end module escape