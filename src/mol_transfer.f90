module mol_transfer

  use parameters
  use naleat
  use utils
  use molecular_emission
  !$ use omp_lib

  use input
  use benchmarks
  use output
  use dust_prop
  use scattering
  use optical_depth
  use ProDiMo, only: read_ProDiMo2mcfost
  use dust_ray_tracing, only: init_directions_ray_tracing
  use dust_transfer, only : compute_stars_map
  use stars
  use grid
  use mem
  use output

  implicit none

  contains

subroutine mol_line_transfer()

  implicit none

  integer :: imol, ibin, iaz

  if (lProDiMo2mcfost) ldust_mol = .true.

!   Default case for molecules and dust
!   optical_length_tot => dust_and_mol_optical_length_tot


  ! Free memory
  call deallocate_em_th_mol()
  lscatt_ray_tracing = .false. ! tmp : scatt ray-tracing has no sense yet for mol emssion
  call init_directions_ray_tracing() ! TODO : can be done later

  do imol=1,n_molecules
     if (lbenchmark_vanzadelhoff1) then
        call readmolecule_benchmark1()
     else if (lbenchmark_vanzadelhoff2) then
        call readmolecule_benchmark2()
     else
        call readmolecule(imol)
     endif

     call alloc_emission_mol(imol) ! does not use much memory

     ! Champ externe
     call init_tab_Cmb_mol()

     if (lDutrey94) then
        call init_GG_Tau_mol()
        call init_molecular_disk(imol)
     else if (lHH30mol) then
        call init_HH_30_mol()
        call init_molecular_disk(imol)
     else if (lbenchmark_vanzadelhoff1) then
        call init_benchmark_vanzadelhoff1()
     else if (lbenchmark_vanzadelhoff2) then
        call init_benchmark_vanzadelhoff2()
     else if (lbenchmark_water1) then
        call init_benchmark_water1()
     else if (lbenchmark_water2) then
        call init_benchmark_water2()
     else if (lbenchmark_water3) then
        call init_benchmark_water3()
     else ! Default case
        call init_molecular_disk(imol)
        ! TODO : missing default case for spherical geometry
     endif

     ! Freeze-out & photo-dissociation eventuels
     if (lfreeze_out) call freeze_out()
     if (lphoto_dissociation) call photo_dissociation()
     if (lwrite_abundance) call write_abundance(imol)

     if (lProDiMo2mcfost) call read_ProDiMo2mcfost(imol)

     ! Dust absorption and emissivity
     call init_dust_mol(imol)

     ! Recompute stellar fluxes at the new wavelengths
     call star_energy_distribution()

     call init_Doppler_profiles(imol)  ! does not depend on nTrans_tot

     ! Initial level populations
     if (.not.lProDiMo2mcfost) then
        ! populations in the optically thin limit
        if (ldouble_RT) call equilibre_othin_mol_pop2()

        ! Initial condition: LTE populations
        call equilibre_LTE_mol(imol)
     endif

     call opacite_mol(imol)
     call integ_tau_mol(imol)

     if (lwrite_mol_column_density) call write_mol_column_density(imol)

     ! Solve for nLTE level populations
     if (.not.lmol_LTE) then
        call NLTE_mol_line_transfer(imol)
     endif

     ! Create molecular emission map via ray-tracing
     if (mol(imol)%lline) then
        do ibin=1,RT_n_incl
           do iaz=1,RT_n_az
              call emission_line_map(imol,ibin,iaz)
              if (ltau_surface) call emission_line_tau_surface_map(imol,tau_surface, ibin,iaz)
              if (lflux_fraction_surface) call  emission_line_energy_fraction_surface_map(imol,flux_fraction, ibin,iaz)
           enddo
        enddo

        call ecriture_spectre(imol)
        if (ltau_surface .or. lflux_fraction_surface) call write_tau_surface(imol)

     endif ! lline

     call dealloc_emission_mol()

  enddo ! imol

  return

end subroutine mol_line_transfer

!********************************************************************

subroutine NLTE_mol_line_transfer(imol)
  ! C. Pinte
  ! 30/06/07

  ! TODO : finish adding packets in step 2 + verify usefulness with profiling
  ! The photon compute time must be negligible

  ! TODO : why does this fail at high optical depth in van Zadelhoff benchmark 1? :-(
  ! TODO : benchmark water3 is not understood ?????


  ! WARNING : this routine is not truly 3D

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: imol
  integer, parameter :: n_rayons_start = 100 ! increasing this reduces the time for step 2, which is the longest
  integer, parameter :: n_rayons_start2 = 100
  integer, parameter :: n_iter2_max = 10
  integer, parameter :: n_rayons_max = n_rayons_start2 * (2**(n_iter2_max-1))
  integer :: n_level_comp
  real, parameter :: precision_sub = 1.0e-3
  real, parameter :: precision = 1.0e-1 !1e-1

  integer :: etape, etape_start, etape_end, iray, n_rayons, icell
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, iv, n_speed, n_cells_done, ibar
  integer, dimension(nb_proc) :: max_n_iter_loc
  real, dimension(nLevels) :: tab_nLevel_old

  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged

  real :: rand, rand2, rand3, fac_etape, factor

  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, maxdiff, norm

  real(kind=dp), dimension(nLevels,nb_proc)  :: pop, pop_old

  logical :: labs

  integer, dimension(2) :: ispeed
  real(kind=dp), dimension(:,:), allocatable :: tab_speed

  integer, dimension(nTrans_tot) :: tab_Trans
  character(len=32) :: step

  labs = .true.

  id = 1

  n_speed = mol(imol)%n_speed_rt ! using the same value now
  n_level_comp = min(mol(imol)%iLevel_max,nLevels)

  do i=1, nTrans_tot
     tab_Trans(i) = i ! we use all the transitions here
  enddo

  if (n_level_comp < 2) call error("n_level_comp must be > 2")
  write(*,*) "NLTE line transfer on", n_level_comp, "levels"

  etape_start=1
  if (lprecise_pop) then
     etape_end = 3
  else
     etape_end = 2  ! 1
  endif

  allocate(tab_speed(-n_speed:n_speed,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_speed')
  tab_speed = 0.0_dp

  allocate(I0(-n_speed:n_speed,nTrans_tot,n_rayons_start,nb_proc), I0c(nTrans_tot,n_rayons_start,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error I0')
  I0 = 0.0_dp
  I0c = 0.0_dp

  if (ldouble_RT) then
     allocate(I02(-n_speed:n_speed,nTrans_tot,n_rayons_start,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error I0')
     I02 = 0.0_dp
  endif

  allocate(ds(n_rayons_start,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error ds')
  ds = 0.0_dp

  allocate(Doppler_P_x_freq(-n_speed:n_speed,n_rayons_start,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error Doppler_P_x_freq')
  Doppler_P_x_freq = 0.0_dp


  ! 2 steps: 1) fixed ray directions to eliminate random fluctuations
  ! and 2) random directions to ensure sufficient sampling
  do etape=etape_start, etape_end  ! Stopping at step 1 is OK for line profiles :-)

     if (etape==1) then
        lfixed_Rays = .true. ;  ispeed(1) = -n_speed ; ispeed(2) = n_speed
        n_rayons = 2 ! 1D RT: 1 ray going up, 1 going down
        iray_start=1
        fac_etape = 1.0
        lprevious_converged = .false.
     else if (etape==2) then
        lfixed_Rays = .true. ;  ispeed(1) = -n_speed ; ispeed(2) = n_speed
        n_rayons = n_rayons_start
        iray_start=1
        fac_etape = 1.0
        lprevious_converged = .false.
     else if (etape==3) then
        lfixed_Rays = .false.;  ispeed(1) = 1 ; ispeed(2) = 1
        n_rayons = n_rayons_start2
        fac_etape = 1.0
        lprevious_converged = .false.

        ! Switch to single-frequency mode
        deallocate(tab_speed,I0,I0c,ds,Doppler_P_x_freq)

        allocate(tab_speed(ispeed(1):ispeed(2),nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_speed')
        tab_speed = 0.0_dp

        allocate(I0(ispeed(1):ispeed(2),nTrans_tot,n_rayons_max,nb_proc), &
             I0c(nTrans_tot,n_rayons_max,nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error I0')
        I0 = 0.0_dp
        I0c = 0.0_dp

        allocate(ds(n_rayons_max,nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error ds')
        ds = 0.0_dp

        allocate(Doppler_P_x_freq(ispeed(1):ispeed(2),n_rayons_max,nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error( 'Allocation error Doppler_P_x_freq')
        Doppler_P_x_freq = 0.0_dp
     endif

     lnotfixed_Rays = .not.lfixed_Rays
     lconverged = .false.
     n_iter = 0

     do while (.not.lconverged)
        n_iter = n_iter + 1
        write(*,*) "--------------------------------------"
        write(*,*) "Step", etape, "Iteration", n_iter

        n_cells_done = 0
        ibar = 0
        max_n_iter_loc = 0

        if (lfixed_Rays) then
           ! Reset the random number generator
           stream = 0.0
           do i=1, nb_proc
              stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
           enddo
        endif

        ! Loop over cells
        call progress_bar(0)
        !$omp parallel &
        !$omp default(none) &
        !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02) &
        !$omp private(argmt,n_iter_loc,lconverged_loc,diff,norm,iv,icell,factor,tab_nLevel_old) &
        !$omp shared(imol,stream,n_rad,nz,n_az,n_rayons,iray_start,Doppler_P_x_freq,tab_nLevel,n_level_comp) &
        !$omp shared(deltaVmax,ispeed,r_grid,z_grid,phi_grid,lcompute_molRT,lkeplerian,n_cells) &
        !$omp shared(tab_speed,lfixed_Rays,lnotfixed_Rays,pop_old,pop,labs,n_speed,max_n_iter_loc,etape,pos_em_cell) &
        !$omp shared(nTrans_tot,tab_Trans,n_cells_done,ibar) &
        !$omp reduction(max:maxdiff)

        !$ id = omp_get_thread_num() + 1
        maxdiff=0

        !$omp do schedule(static)
        do icell=1, n_cells
           if (lcompute_molRT(icell)) then
              tab_nLevel_old(1:n_level_comp) = tab_nLevel(icell,1:n_level_comp)

              ! Uniform sampling of the line profile
              factor = deltaVmax(icell) / real(n_speed,kind=dp)
              if (lfixed_rays) then
                 do iv=-n_speed, n_speed
                    tab_speed(iv,id) =  factor * real(iv,kind=dp)
                 enddo ! iv
              endif

              ! Ray propagation
              do iray=iray_start, iray_start-1+n_rayons

                 if (etape==1) then
                    ! Position = cell centre
                    x0 = r_grid(icell) * cos(phi_grid(icell))
                    y0 = r_grid(icell) * sin(phi_grid(icell))
                    z0 = z_grid(icell)

                    if (lkeplerian) then
                       ! Vertical direction
                       if (iray==1) then
                          w0=1.0_dp
                       else
                          w0=-1.0_dp
                       endif
                       u0 = 0.0_dp
                       v0 = 0.0_dp
                    else
                       norm = sqrt(x0*x0 + y0*y0 + z0*z0)
                       if (iray==1) then
                          u0 = x0/norm
                          v0 = y0/norm
                          w0 = z0/norm
                       else
                          u0 = -x0/norm
                          v0 = -y0/norm
                          w0 = -z0/norm
                       endif
                    endif
                 else
                    ! Random position within the cell
                    rand  = sprng(stream(id))
                    rand2 = sprng(stream(id))
                    rand3 = sprng(stream(id))
                    call pos_em_cell(icell ,rand,rand2,rand3,x0,y0,z0)

                    ! Random propagation direction
                    rand = sprng(stream(id))
                    W0 = 2.0_dp * rand - 1.0_dp
                    W02 =  1.0_dp - W0*W0
                    SRW02 = sqrt(W02)
                    rand = sprng(stream(id))
                    ARGMT = PI * (2.0_dp * rand - 1.0_dp)
                    U0 = SRW02 * cos(ARGMT)
                    V0 = SRW02 * sin(ARGMT)
                 endif

                 ! Random sampling of the velocity field
                 if (lnotfixed_Rays) then
                    do iv=ispeed(1),ispeed(2)
                       rand = sprng(stream(id)) ; tab_speed(iv,id) =  2.0_dp * (rand - 0.5_dp) * deltaVmax(icell)
                    enddo
                 endif

                 ! Integration along the ray
                 call integ_ray_mol(id,imol,icell,x0,y0,z0,u0,v0,w0,iray,labs, ispeed,tab_speed(:,id), &
                      nTrans_tot, tab_Trans)
              enddo ! iray


              ! Solve the statistical equilibrium
              n_iter_loc = 0
              pop(:,id) = tab_nLevel(icell,:)
              lconverged_loc = .false.
              ! Loop to converge the local radiation field and populations
              ! with a fixed external field
              do while (.not.lconverged_loc)
                 n_iter_loc = n_iter_loc + 1

                 ! Save old local populations
                 pop_old(:,id) = pop(:,id)

                 ! Compute the radiation field
                 call J_mol_loc(id,icell,n_rayons,ispeed)  ! includes loops over transitions

                 call equilibre_rad_mol_loc(id,icell)
                 pop(:,id) = tab_nLevel(icell,:)

                 ! Local convergence criterion
                 diff = maxval( abs(pop(1:n_level_comp,id) - pop_old(1:n_level_comp,id)) &
                      / (pop_old(1:n_level_comp,id) + 10*tiny_real) )

                 if (diff < precision_sub) then
                    lconverged_loc = .true.
                 else
                    ! Not converged, recompute opacities and source functions
                    call opacite_mol_loc(icell,imol)
                 endif

              enddo ! while: local field convergence
              if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc

              diff = maxval( abs( tab_nLevel(icell,1:n_level_comp) - tab_nLevel_old(1:n_level_comp) ) / &
                   tab_nLevel_old(1:n_level_comp) + 10*tiny_real)

              if (diff > maxdiff) maxdiff = diff
           endif ! lcompute_molRT

           ! Progress bar
           !$omp atomic
           n_cells_done = n_cells_done + 1
           if (real(n_cells_done) > 0.02*ibar*n_cells) then
              call progress_bar(ibar)
              !$omp atomic
              ibar = ibar+1
           endif

        enddo ! icell
        !$omp end do
        !$omp end parallel
        call progress_bar(50)

        ! Global convergence criterion
        write(*,*) maxval(max_n_iter_loc), "sub-iterations"
        write(*,*) "Relative difference =", real(maxdiff)

        write(*,*) "Threshold =", precision*fac_etape

        if (maxdiff < precision*fac_etape) then
           if (lprevious_converged) then
              lconverged = .true.
           else
              lprevious_converged = .true.
           endif
        else
           lprevious_converged = .false.
           if (.not.lfixed_rays) then
              n_rayons = n_rayons * 2
             ! if (n_rayons > n_rayons_max) then
              if (n_iter >= n_iter2_max) then
              write(*,*) "Warning : not enough rays to converge !!"
                 lconverged = .true.
              endif

              ! Continue by computing twice as many rays
              ! and add them to the previously computed set
!              iray_start = iray_start + n_rayons

           endif
        endif

        !write(*,*) "STAT", minval(tab_nLevel(:,1:n_level_comp)), maxval(tab_nLevel(:,1:n_level_comp))
        call integ_tau_mol(imol)

     enddo ! while: global convergence

     write(step, "(A5, I1)") "_step", etape

     call ecriture_pops(imol,step)
     call ecriture_Tex(imol,step)
  enddo ! etape

  deallocate(ds, Doppler_P_x_freq, I0, I0c)

  return

end subroutine NLTE_mol_line_transfer

!***********************************************************

subroutine emission_line_map(imol,ibin,iaz)
  ! Create the molecular emission map
  ! (or the spectrum if there is only one pixel)
  ! via ray-tracing in a given direction
  ! C. Pinte
  ! 12/04/07

  implicit none

  integer, intent(in) :: imol, ibin, iaz

  real(kind=dp) :: x0,y0,z0,l, u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcorner
  real(kind=dp) :: taille_pix
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36  ! OK, also works with n_rad_RT = 1000
  integer, parameter :: n_ray_star = 1000
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp) :: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, nTrans_raytracing

  integer :: n_speed_rt, n_speed_center_rt, n_extraV_rt, lambda, iv
  real :: vmin_center_rt, vmax_center_rt, extra_deltaV_rt
  logical :: lresolved, l_sym_ima_mol

 ! Line-of-sight direction for ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  n_speed_rt = mol(imol)%n_speed_rt
  n_speed_center_rt = mol(imol)%n_speed_center_rt
  n_extraV_rt = mol(imol)%n_extraV_rt

  vmax_center_rt = mol(imol)%vmax_center_rt
  vmin_center_rt = mol(imol)%vmin_center_rt
  extra_deltaV_rt = mol(imol)%extra_deltaV_rt

  nTrans_raytracing = mol(imol)%nTrans_raytracing

  if ((ibin == 1).and.(iaz==1)) then
     allocate(I0(n_speed_rt,nTrans_raytracing,1,nb_proc),I0c(nTrans_raytracing,1,nb_proc))
     allocate(tab_speed_rt(n_speed_rt))

     ! line centre
     tab_speed_rt(1:n_speed_center_rt) = span(vmin_center_rt,vmax_center_rt,n_speed_center_rt)
     !! line wings
     !tab_speed_rt(n_speed_center_rt+1:n_speed_rt) = indgen(n_extraV_rt) * extra_deltaV_rt + vmax_center_rt
     !do i = -n_speed_rt, -n_speed_center_rt-1
     !   tab_speed_rt(i) = (i+n_speed_center_rt) * extra_deltaV_rt - vmax_center_rt
     !enddo

     if (lorigine) then
        allocate(origine_mol(n_speed_rt,nTrans_raytracing,n_cells,nb_proc))
        origine_mol = 0.0
     endif
  endif

  I0 = 0.0_dp
  I0c = 0.0_dp
  if (lorigine) origine_mol = 0.0

  ! Definition of the image plane basis vectors in the universal frame

  ! Image x-vector without PA: lies in the (x,y) plane and is orthogonal to uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad),sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Image x-vector with PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo: can simplify since the rotation axis is perpendicular to x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Image y-vector with PA: orthogonal to x_plan_image and uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! initial position outside the model (on the observer side)
  ! = image centre
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Bottom-left corner of the image
  Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)

  l_sym_ima_mol = mol(imol)%l_sym_ima

  if (RT_line_method == 1) then ! method 1: log sampling
     ! No sub-pixel since the pixels are not square
     n_iter_min = 1
     n_iter_max = 1

     ! dx and dy are only required for stellar map here
     dx(:) = x_plan_image * (map_size/zoom)  ! en au
     dy(:) = y_plan_image * (map_size/zoom)  ! en au

     i = 1
     j = 1
     lresolved = .false.

     rmin_RT = max(w*0.9_dp,0.05_dp) * Rmin
     rmax_RT = 2.0_dp * Rmax

     tab_r(1) = rmin_RT
     fact_r = exp( (1.0_dp/(real(n_rad_RT,kind=dp) -1))*log(rmax_RT/rmin_RT) )

     do ri_RT = 2, n_rad_RT
        tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
     enddo

     fact_A = sqrt(pi * (fact_r - 1.0_dp/fact_r)  / n_phi_RT )

     ! Loop over sampling rays
     !$omp parallel &
     !$omp default(none) &
     !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
     !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
     !$omp shared(n_iter_min,n_iter_max,l_sym_ima_mol,cst_phi,imol,ibin,iaz)
     id =1 ! for sequential code

     if (l_sym_ima_mol) then
        cst_phi = pi  / real(n_phi_RT,kind=dp)
     else
        cst_phi = two_pi  / real(n_phi_RT,kind=dp)
     endif

     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1

        r = tab_r(ri_RT)
        taille_pix =  fact_A * r ! square root of the pixel area

        do phi_RT=1,n_phi_RT ! from 0 to pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)

           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image ! This is actually the centre since dx = dy = 0.
           call intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

  else ! method 2: linear sampling with sub-pixels
     lresolved = .true.

     ! Vectors defining the pixels (dx,dy) in the universal frame
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en au
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     if (l_sym_ima_mol) then
        npix_x_max = npix_x/2 + modulo(npix_x,2)
     else
        npix_x_max = npix_x
     endif

     ! Loop over image pixels
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y) &
     !$omp shared(n_iter_min,n_iter_max,imol,ibin,iaz)

     id =1 ! for sequential code
     n_iter_min = 1 ! 3
     n_iter_max = 1 ! 6

     !$omp do schedule(dynamic,1)
     do i = 1,npix_x_max
        !$ id = omp_get_thread_num() + 1
        do j = 1,npix_y
           !write(*,*) i,j
           ! Bottom-left corner of the pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           call intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

  endif

  ! --------------------------
  ! Ajout flux star
  ! --------------------------
  do i = 1, mol(imol)%nTrans_raytracing
     lambda =  mol(imol)%index_trans_ray_tracing(i) ! == iTrans
     call compute_stars_map(lambda, ibin, iaz, u, v, w, taille_pix, dx, dy, lresolved)

     do iv = 1, n_speed_rt
        spectrum(:,:,iv,i,ibin,iaz) = spectrum(:,:,iv,i,ibin,iaz) + stars_map(:,:,1)
     enddo
     continu(:,:,i,ibin,iaz) = continu(:,:,i,ibin,iaz) + stars_map(:,:,1)
  enddo

  return

end subroutine emission_line_map

!***********************************************************

subroutine intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! Computes the intensity of a square pixel with arbitrary size, position and orientation
  ! using a ray-tracing method
  ! (u,v,w) pointe vers l'observateur
  ! Romberg integration to determine the number of sub-pixels
  ! necessaire
  ! Unite : W.m-2 : nu.F_nu
  ! C. Pinte
  ! 12/04/07

  implicit none

  integer, intent(in) :: ipix,jpix,id, imol, n_iter_min, n_iter_max, ibin, iaz
  real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
  real(kind=dp), intent(in) :: pixelsize,u,v,w
  real(kind=dp), dimension(:,:), allocatable :: IP, IP_old
  real(kind=dp), dimension(:), allocatable :: IPc

  integer, parameter :: maxSubPixels = 32

  real(kind=dp) :: x0,y0,z0,u0,v0,w0
  real(kind=dp), dimension(3) :: sdx, sdy
  real :: npix2, diff

  real, parameter :: precision = 1.e-2
  integer :: i, j, subpixels, iray, ri, zj, phik, icell, iTrans, iter, n_speed_rt, nTrans_raytracing

  logical :: lintersect, labs

  integer, dimension(2) :: ispeed

  n_speed_rt = mol(imol)%n_speed_rt
  nTrans_raytracing = mol(imol)%nTrans_raytracing

  allocate(IP(n_speed_rt,nTrans_raytracing), IP_old(n_speed_rt,nTrans_raytracing), IPc(nTrans_raytracing))

  ispeed(1) = 1 ; ispeed(2) = n_speed_rt


  labs = .false.

  ! Ray tracing: propagating in the opposite direction
  u0 = -u ; v0 = -v ; w0 = -w

  IP = 0.0_dp
  IPc = 0.0_dp

  ! the number of subpixels in x is 2^(iter-1)
  subpixels = 1
  iter = 1

  infinie : do ! Infinite loop until the pixel has converged
     npix2 =  real(subpixels)**2
     IP_old = IP
     IP = 0.0_dp
     IPc = 0.0_dp

     ! Vectors defining the sub-pixels
     sdx(:) = dx(:) / real(subpixels,kind=dp)
     sdy(:) = dy(:) / real(subpixels,kind=dp)

     iray = 1

     ! The observer is outside the grid
     ri = 2*n_rad ; zj=1 ; phik=1

     ! Loop over sub-pixels computing the intensity at the centre
     ! of each sub-pixel
     do i = 1,subpixels
        do j = 1,subpixels
           ! Centre of the sub-pixel
           x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
           y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
           z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)

           ! Start at the grid boundary: backward propagation
           call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

           if (lintersect) then ! Intersected the grid, meaning potential flux
              call integ_ray_mol(id,imol,icell,x0,y0,z0,u0,v0,w0,iray,labs,ispeed,tab_speed_rt, &
                   nTrans_raytracing, mol(imol)%index_trans_ray_tracing)
              ! Flux received in the pixel
              IP(:,:) = IP(:,:) +  I0(:,:,iray,id)
              IPc(:) = IPc(:) +  I0c(:,iray,id)
           else
              ! Only the CMB remains
              ! TODO : not sure if I want the CMB outside the Disk, ...
              ! IP(:,:) = IP(:,:) +  Cmb(ispeed,tab_speed)
           endif
        enddo !j
     enddo !i

     IP = IP / npix2
     IPc = IPc / npix2

     if (iter < n_iter_min) then
        ! Iterate by default
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! Stop to prevent infinite loops
        ! write(*,*) "Warning : converging pb in ray-tracing"
        ! write(*,*) " Pixel", ipix, jpix
        exit infinie
     else
        ! Test the difference
        diff = maxval( abs(IP - IP_old) / (IP + 1e-300_dp) )
        if (diff > precision ) then
           ! Not converged
           subpixels = subpixels * 2
        else
           ! Converged
           exit infinie
        endif
     endif ! iter

     iter = iter + 1

     ! TODO : Romberg integration
!!$     if(any(abs((I - oldintensite_pixel)) > precision * I)) then
!!$        oldintensite_pixel = I
!!$        ! Isn't there a better way to use the lower resolution computations ??
!!$
!!$        subpixels = subpixels * 2
!!$        !if(subpixels .gt. 15) write(*,*)"large",index
!!$     else
!!$        I = (real(2**(log(real(subpixels))/log(2.d0)))*I - Oldintensite_pixel) &
!!$             /(real(2**(log(real(subpixels))/log(2.d0)))-1.d0) ! Richardson Extrapolation
!!$        ! Ok but only uses the last 2 calculations: there must be a better way !!!
!!$
!!$     endif
  enddo infinie

  ! Take into account the pixel surface (in sr)
  IP = IP * (pixelsize / (distance*pc_to_AU) )**2
  IPc = IPc * (pixelsize / (distance*pc_to_AU) )**2

  ! and multiply by the frequency to get nu.F_nu
  ! Warning IP, IPc are smaller array (dimension mol(imol)%nTrans_raytracing)
  do i=1,mol(imol)%nTrans_raytracing
     iTrans = mol(imol)%index_trans_ray_tracing(i)
     IP(:,i) = IP(:,i) * transfreq(iTrans)
     IPc(i) = IPc(i) * transfreq(iTrans)
  enddo
  ! Units tested OK for CMB
  ! unconvolved line profile tested ok with torus

  if (RT_line_method==1) then ! Implicit summation over pixels
     spectrum(1,1,:,:,ibin,iaz) = spectrum(1,1,:,:,ibin,iaz) + IP(:,:)
     continu(1,1,:,ibin,iaz) = continu(1,1,:,ibin,iaz) + IPc(:)
  else
     spectrum(ipix,jpix,:,:,ibin,iaz) = IP(:,:)
     continu(ipix,jpix,:,ibin,iaz) = IPc(:)
  endif

  return

end subroutine intensite_pixel_mol

!***********************************************************

subroutine init_dust_mol(imol)
  ! computes the opacities and emissivities of the cells
  ! at the wavelengths of the molecular lines
  ! to take into account the radiative interaction between
  ! the 2 phases.
  ! C. Pinte
  ! 17/10/07
  ! TODO : handle the case where albedo (and thus scattering) is not negligible

  use mem, only : realloc_dust_mol, clean_mem_dust_mol

  implicit none

  integer, intent(in) :: imol

  integer :: iTrans
  integer, target :: icell
  integer, pointer :: p_icell
  real(kind=dp) :: freq!, Jnu
  real :: T, wl, kap

  real(kind=dp) :: cst_E

  real, parameter :: gas_dust = 100
  real, parameter :: delta_lambda = 0.025

  integer :: iTrans_min, iTrans_max

  iTrans_min = mol(imol)%iTrans_min
  iTrans_max = mol(imol)%iTrans_max

  cst_E=2.0*hp*c_light**2

  ! Reallocate the dust property arrays
  ! n_lambda =   mol(imol)%nTrans_raytracing ! dust opacities considered constant over the line profile
  n_lambda = nTrans_tot ! dust opacities considered constant over the line profile

  ! We are only interested in absorption properties: no need for Mueller matrices
  ! -> no polarisation, we use an HG function
  scattering_method=1 ; lscattering_method1 = .true.
  aniso_method = 2 ; lmethod_aniso1 = .false.

  lsepar_pola = .false.
  ltemp = .false.
  lmono = .true. ! equivalent to sed2 mode

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell1
  endif

  call realloc_dust_mol(imol)

  tab_lambda = 1e-30 ! to avoid error in stars.f90

  if (ldust_mol) then
     ! Wavelength array
     do iTrans=iTrans_min,iTrans_max
        tab_lambda(iTrans) = c_light/Transfreq(iTrans) * 1.0e6 ! in microns
        tab_lambda_sup(iTrans)= tab_lambda(iTrans)*delta_lambda
        tab_lambda_inf(iTrans)= tab_lambda(iTrans)/delta_lambda
        tab_delta_lambda(iTrans) = tab_lambda_sup(iTrans) - tab_lambda_inf(iTrans)
     enddo

     if (lbenchmark_water3) then ! power-law opacity
        write(*,*) "WARNING : hard-coded gas_dust =", gas_dust

        do iTrans=1,nTrans_tot
           wl = tab_lambda(iTrans)

           ! Opacity law (cm^2 per g of dust)
           if (wl > 250) then
              kap = 10. * (wl/250.)**(-2.0)
           else
              kap = 10. * (wl/250.)**(-1.3)
           endif

           ! Multiplication by density
           ! AU_to_cm**2 because we want kappa_abs_LTE in AU-1
           write(*,*) "TODO : the water benchmark 3 needs to be updated for cell pointer in opacity table"
           do icell=1,n_cells
              kappa_abs_LTE(icell,iTrans) =  kap * (gas_density(icell) * cm_to_m**3) * mu_mH / &
                   gas_dust / cm_to_AU
           enddo

        enddo ! iTrans

        ! No scattering
        kappa(:,:) = kappa_abs_LTE(:,:)

     else ! default case
        call init_optical_indices()

        ! Recompute the optical properties
        write(*,*) "Computing dust properties for", nTrans_tot, "wavelength"
        do iTrans=iTrans_min,iTrans_max
           call prop_grains(iTrans)
           call opacity(iTrans, iTrans, no_scatt=.true.)
        enddo
     endif


     ! Unit change: kappa in m-1 for RT in the lines !!!!!!!
     ! Will be converted back to AU-1 in opacite_mol_loc
     !kappa = kappa * m_to_AU
     !kappa_abs_eg = kappa_abs_eg * m_to_AU

     ! Computation of the dust emissivity
     do iTrans=iTrans_min,iTrans_max
        freq = Transfreq(iTrans)

        ! TODO : speed up this loop via Bnu_disk routine (does this take time ???)
        ! TODO : generalize for all grain types (should already exist no ???)
        ! TODO : this could also be used to do ray-tracing in the continuum 8-)
        do icell=1, n_cells
           !-- ! Interpolate radiation field across wavelength
           !-- if (lProDiMo2mcfost) then
           !--    Jnu = interp(m2p%Jnu(ri,zj,:), m2p%wavelengths, real(tab_lambda(iTrans)))
           !-- else
           !--    Jnu = 0.0 ! todo : to take scattering into account
           !-- endif

           T = Tdust(icell)
           ! On ne fait que du scattering isotropique dans les raies pour le moment ...
           emissivite_dust(icell,iTrans) = kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell) * Bnu(freq,T) ! + kappa_sca(iTrans,ri,zj,phik) * Jnu
        enddo ! icell
     enddo ! itrans

  else ! .not.ldust_mol
     kappa = 0.0
     emissivite_dust = 0.0
  endif

  ! Deallocation les tableaux dont on a pas besoin pour la suite
  call clean_mem_dust_mol()

  return

end subroutine init_dust_mol

!***********************************************************

subroutine init_molecular_disk(imol)
  ! definie les tableaux vfield, v_turb et tab_abundance
  ! et lcompute_molRT

  integer, intent(in) :: imol

  logical, save :: lfirst_time = .true.
  real(dp) :: factor
  integer :: icell

  ldust_mol  = .true.
  lkeplerian = .true.

  ! Turn image symmetry off if we read the velocity field from a file,
  ! because the velocity field may not be symmetric.
  if (lvelocity_file) l_sym_ima = .false.

  if (lfirst_time) then
     lfirst_time = .false.
     call init_Tgas()

     ! Velocity field in  m.s-1
     ! Warning : assume all stars are at the center of the disk
     if (.not.lVoronoi) then ! Velocities are defined from SPH files in Voronoi mode
        if (lcylindrical_rotation) then ! Midplane Keplerian velocity
           do icell=1, n_cells
              vfield(icell) = sqrt(Ggrav * sum(star%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
           enddo
        else ! dependance en z
           do icell=1, n_cells
              vfield(icell) = sqrt(Ggrav * sum(star%M) * Msun_to_kg * r_grid(icell)**2 / &
                   ((r_grid(icell)**2 + z_grid(icell)**2)**1.5 * AU_to_m) )
           enddo
        endif
     endif

     if (lvturb_in_cs) then
        factor = vitesse_turb**2
        do icell=1, n_cells
           v_turb2(icell) =  (kb*Tcin(icell) / (mu_mH * g_to_kg)) * factor  ! cs**2 * factor
        enddo
     else
        v_turb2(:) = vitesse_turb**2 ! constant vturb
     endif
  endif ! lfirst_time

  ! Abondance
  call init_abundance(imol)

  return

end subroutine init_molecular_disk

!***********************************************************

subroutine init_Tgas()

  use ML_prodimo, only : xgb_predict_Tgas

  if (lML) then
     write(*,*) "Predicting gas temperature"
     call xgb_predict_Tgas()
     write(*,*) "Max gas temperature=", maxval(Tcin)
  else
     ! Temperature gaz = dust
     if (lcorrect_Tgas) then
        write(*,*) "Correcting Tgas by", correct_Tgas
        Tcin(:) = Tdust(:)  * correct_Tgas
     else
        Tcin(:) = Tdust(:)
     endif
  endif

  return

end subroutine init_Tgas

!***********************************************************

subroutine init_abundance(imol)

  use ML_prodimo, only : xgb_predict_abundance

  integer, intent(in) :: imol

  integer :: icell

  if (lML) then
     write(*,*) "Predicting  molecular abundance"
     call xgb_predict_abundance(imol)
     write(*,*) "Min-Max abundance=", maxval(tab_abundance)
  else
     if (mol(imol)%lcst_abundance) then
        write(*,*) "Setting constant abundance"
        tab_abundance = mol(imol)%abundance
     else
        write(*,*) "Reading abundance from file"
        call read_abundance(imol)
     endif
  endif

  do icell=1, n_cells
     lcompute_molRT(icell) = (tab_abundance(icell) > tiny_real) .and. &
          (gas_density(icell) > tiny_real) .and. (Tcin(icell) > 1.)
  enddo

  return

end subroutine init_abundance

!***********************************************************

subroutine emission_line_tau_surface_map(imol,tau,ibin,iaz)

  real, intent(in) :: tau
  integer, intent(in) :: imol, ibin, iaz
  real(kind=dp) :: u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcenter

  integer :: i,j, id, icell, iTrans
  real(kind=dp) :: l, taille_pix, x0, y0, z0, u0, v0, w0
  logical :: lintersect, flag_sortie, lpacket_alive
  integer, dimension(4) :: ispeed

  ! Line-of-sight direction for ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition of the image plane basis vectors in the universal frame

  ! Image x-vector without PA: lies in the (x,y) plane and is orthogonal to uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Image x-vector with PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo: can simplify since the rotation axis is perpendicular to x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Image y-vector with PA: orthogonal to x_plan_image and uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! initial position outside the model (on the observer side)
  ! = image centre
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Vectors defining the pixels (dx,dy) in the universal frame
  taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
  dx(:) = x_plan_image * taille_pix
  dy(:) = y_plan_image * taille_pix

  ! Bottom-left corner of the image
  Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

  ! We only consider the 1st transition for now
  iTrans = mol(imol)%index_trans_ray_tracing(1)

  ! Array velocity
  !nTrans_raytracing = mol(imol)%nTrans_raytracing
  ispeed(1) = 1 ; ispeed(2) = mol(imol)%n_speed_rt

  ! Loop over image pixels
  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,id,icell,lintersect,x0,y0,z0,u0,v0,w0) &
  !$omp private(flag_sortie,lpacket_alive,pixelcenter) &
  !$omp shared(tau,Icorner,imol,iTrans,dx,dy,u,v,w,ispeed,tab_speed_rt) &
  !$omp shared(taille_pix,npix_x,npix_y,ibin,iaz,tau_surface_map,move_to_grid)
  id = 1 ! for sequential code

  tau_surface_map = 0.0_dp

  !$omp do schedule(dynamic,1)
  do i = 1, npix_x
     !$ id = omp_get_thread_num() + 1
     do j = 1,npix_y
        ! Bottom-left corner of the pixel
        pixelcenter(:,id) = Icorner(:) + (i-0.5_dp) * dx(:) + (j-0.5_dp) * dy(:)

        x0 = pixelcenter(1,id)
        y0 = pixelcenter(2,id)
        z0 = pixelcenter(3,id)

        ! Ray tracing: propagating in the opposite direction
        u0 = -u ; v0 = -v ; w0 = -w

        ! On se met au bord de la grid : propagation a l'envers
        call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

        if (lintersect) then ! On rencontre la grid, on a potentiellement du flux
           lpacket_alive = .true.
           call physical_length_mol(imol,iTrans,icell,x0,y0,z0,u0,v0,w0,ispeed,tab_speed_rt,tau,flag_sortie)
           if (flag_sortie) then ! We do not reach the surface tau=1
              tau_surface_map(i,j,ibin,iaz,:,id) = 0.0
           else
              tau_surface_map(i,j,ibin,iaz,1,id) = x0
              tau_surface_map(i,j,ibin,iaz,2,id) = y0
              tau_surface_map(i,j,ibin,iaz,3,id) = z0
           endif
        else ! We do not reach the disk
           tau_surface_map(i,j,ibin,iaz,:,id) = 0.0
        endif

     enddo !j
  enddo !i
  !$omp end do
  !$omp end parallel

  return

end subroutine emission_line_tau_surface_map

!***********************************************************

subroutine emission_line_energy_fraction_surface_map(imol,flux_fraction,ibin,iaz)

  real(kind=dp), intent(in) :: flux_fraction
  integer, intent(in) :: imol, ibin, iaz
  real(kind=dp) :: u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcenter

  integer :: i,j, id, icell, iTrans, iiTrans
  real(kind=dp) :: l, taille_pix, x0, y0, z0, u0, v0, w0, Flux, factor
  logical :: lintersect, flag_sortie, lpacket_alive
  integer, dimension(4) :: ispeed

  ! Line-of-sight direction for ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  ! Definition of the image plane basis vectors in the universal frame

  ! Image x-vector without PA: lies in the (x,y) plane and is orthogonal to uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad), sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Image x-vector with PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo: can simplify since the rotation axis is perpendicular to x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Image y-vector with PA: orthogonal to x_plan_image and uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! initial position outside the model (on the observer side)
  ! = image centre
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Vectors defining the pixels (dx,dy) in the universal frame
  taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en au
  dx(:) = x_plan_image * taille_pix
  dy(:) = y_plan_image * taille_pix

  ! Bottom-left corner of the image
  Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

  ! We only consider the 1st transition for now
  iTrans = 1
  iiTrans = mol(imol)%index_trans_ray_tracing(iTrans)

  ! Corrective factor : we want the flux before we take into account the pixel size
  factor = flux_fraction / (taille_pix / (distance*pc_to_AU))**2 /  transfreq(iiTrans)

  ! Array velocity
  !nTrans_raytracing = mol(imol)%nTrans_raytracing
  ispeed(1) = 1 ; ispeed(2) = mol(imol)%n_speed_rt

  ! Loop over image pixels
  !$omp parallel &
  !$omp default(none) &
  !$omp private(i,j,id,icell,lintersect,x0,y0,z0,u0,v0,w0) &
  !$omp private(flag_sortie,lpacket_alive,pixelcenter,Flux) &
  !$omp shared(flux_fraction,Icorner,imol,iTrans,iiTrans,dx,dy,u,v,w,ispeed,tab_speed_rt) &
  !$omp shared(taille_pix,npix_x,npix_y,ibin,iaz,tau_surface_map,move_to_grid,spectrum,factor)
  id = 1 ! for sequential code

  tau_surface_map = 0.0_dp

  !$omp do schedule(dynamic,1)
  do i = 1, npix_x
     !$ id = omp_get_thread_num() + 1
     do j = 1,npix_y
        ! Bottom-left corner of the pixel
        pixelcenter(:,id) = Icorner(:) + (i-0.5_dp) * dx(:) + (j-0.5_dp) * dy(:)

        x0 = pixelcenter(1,id)
        y0 = pixelcenter(2,id)
        z0 = pixelcenter(3,id)

        ! Ray tracing: propagating in the opposite direction
        u0 = -u ; v0 = -v ; w0 = -w

        ! On se met au bord de la grid : propagation a l'envers
        call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)

        ! Maximum flux in pixel
        Flux = maxval(spectrum(i,j,:,iTrans,ibin,iaz)) * factor

        if (lintersect) then ! On rencontre la grid, on a potentiellement du flux
           lpacket_alive = .true.
           call physical_length_mol_Flux(imol,iiTrans,icell,x0,y0,z0,u0,v0,w0,ispeed,tab_speed_rt,Flux,flag_sortie)
           if (flag_sortie) then ! We do not reach the surface tau=1
              tau_surface_map(i,j,ibin,iaz,:,id) = 0.0
           else
              tau_surface_map(i,j,ibin,iaz,1,id) = x0
              tau_surface_map(i,j,ibin,iaz,2,id) = y0
              tau_surface_map(i,j,ibin,iaz,3,id) = z0
           endif
        else ! We do not reach the disk
           tau_surface_map(i,j,ibin,iaz,:,id) = 0.0
        endif

     enddo !j
  enddo !i
  !$omp end do
  !$omp end parallel

  return

end subroutine emission_line_energy_fraction_surface_map

!***********************************************************

end module mol_transfer
