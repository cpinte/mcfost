module mol_transfer

  use parametres
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
  
  use accelerate
  use math, only : flatten, reform

  implicit none
  
  type (Ng) :: Ngmol

  contains

subroutine mol_line_transfer()

  implicit none

  integer :: imol, ibin, iaz

  if (lProDiMo2mcfost) ldust_mol = .true.
  
  !Default case for molecules and dust
  optical_length_tot => dust_and_mol_optical_length_tot

  ! Liberation memoire
  call deallocate_em_th_mol()
  lscatt_ray_tracing = .false. ! tmp : scatt ray-tracing has no sense yet for mol emssion
  call init_directions_ray_tracing() ! TODO : on peut le faire apres

  do imol=1,n_molecules
     if (lbenchmark_vanzadelhoff1) then
        call readmolecule_benchmark1()
     else if (lbenchmark_vanzadelhoff2) then
        call readmolecule_benchmark2()
     else
        call readmolecule(imol)
     endif

     call alloc_emission_mol(imol) ! ne prend pas bcp de memoire

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
     else ! Cas par defaut
        call init_molecular_disk(imol)
        ! TODO : il manque le cas defaut pour geometrie spherique
     endif

     ! Freeze-out & photo-dissociation eventuels
     if (lfreeze_out) call freeze_out()
     if (lphoto_dissociation) call photo_dissociation()

     if (lProDiMo2mcfost) call read_ProDiMo2mcfost(imol)

     ! Absorption et emissivite poussiere
     call init_dust_mol(imol)

     ! recalcul des flux stellaires aux nouvelles longeurs d'onde
     call repartition_energie_etoiles()

     call init_Doppler_profiles(imol)  ! ne depend pas de nTrans_tot

     ! Population des niveaux initiale
     if (.not.lProDiMo2mcfost) then
        ! population dans le cas limite optiquement mince
        if (ldouble_RT) call equilibre_othin_mol_pop2()

        ! Condition initiale : population a l'ETL
        call equilibre_LTE_mol()
     endif

     call opacite_mol(imol)
     call integ_tau_mol(imol)

     ! Resolution population des niveaux nLTE
     if (.not.lmol_LTE) then
        call NLTE_mol_line_transfer(imol)
        call ecriture_pops(imol)
        call ecriture_Tex(imol)
     endif

     !--- Creation carte emission moleculaire : ray-tracing
     if (mol(imol)%lline) then
        do ibin=1,RT_n_incl
           do iaz=1,RT_n_az
              call emission_line_map(imol,ibin,iaz)
           enddo
        enddo

        call ecriture_spectre(imol)
     endif ! lline

     call dealloc_emission_mol()

  enddo ! imol

  return

end subroutine mol_line_transfer

!********************************************************************

subroutine NLTE_mol_line_transfer(imol)
  ! C. Pinte
  ! 30/06/07

  ! TODO : finir l'ajout des paquets dans etape 2 + verifier utilite avec profiling
  ! Il faut que le tps de calcul des photons soit negligeable

  ! TODO : pourquoi ca merde a haute profondeur optique dans le benchmark 1 de van Zadelhoff ??? :-(
  ! TODO : je capte pas le benchmark water3 ?????


  ! WARNING : cette routine n'est pas vraiment 3D

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: imol
                                         
  integer, parameter :: n_rayons_start = 100 ! l'augmenter permet de reduire le tps de l'etape 2 qui est la plus longue
  integer, parameter :: n_rayons_start2 = 100
  integer, parameter :: n_iter2_max = 10
  integer, parameter :: n_rayons_max = n_rayons_start2 * (2**(n_iter2_max-1))
  integer :: n_level_comp
  real, parameter :: precision_sub = 1.0e-3
  real, parameter :: precision = 1.0e-1 !1e-1

  integer :: etape, etape_start, etape_end, iray, n_rayons, icell
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, iv, n_speed
  integer, dimension(nb_proc) :: max_n_iter_loc

  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged

  real :: rand, rand2, rand3, fac_etape, factor

  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, maxdiff, norme

  real(kind=dp), dimension(nLevels,nb_proc)  :: pop, pop_old

  logical :: labs

  integer, dimension(2) :: ispeed
  real(kind=dp), dimension(:,:), allocatable :: tab_speed

  integer, dimension(nTrans_tot) :: tab_Trans
  
  !Ng'acceleration
  logical :: accelerated, ng_rest
  integer :: iorder, i0_rest, n_iter_accel
  real(kind=dp), allocatable, dimension(:) :: flatpops, popsT(:,:)

  labs = .true.

  id = 1
  
  n_speed = mol(imol)%n_speed_rt ! j'utilise le meme maintenant
  n_level_comp = min(mol(imol)%iLevel_max,nLevels)
  
  if (lNg_acceleration) then
   write(*,*) " *********** "
    write(*,*) "Ng acceleration not tested for molecules yet"
   write(*,*) " *********** "
  endif
  if (lNg_acceleration) then !Nord >0 already tested
    Write(*,*) " Allocating space for Ng structure"
    CALL initNg(n_cells*n_level_comp,iNg_Ndelay, iNg_Norder, iNg_Nperiod, Ngmol)
    n_iter_accel = 0
    i0_rest = 0
    ng_rest = .false.
    allocate(flatpops(n_cells*n_level_comp))
    allocate(popsT(n_level_comp, n_cells)); popsT = 0d0
  endif

  do i=1, nTrans_tot
     tab_Trans(i) = 1 ! we use all the transitions here
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

  allocate(I0(-n_speed:n_speed,nTrans_tot,n_rayons_start,nb_proc), &
       I0c(nTrans_tot,n_rayons_start,nb_proc), stat=alloc_status)
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

  ! 2 etapes : 1) direction des rayons fixee pour eliminer les fluctuations aletoires
  ! et 2) direction aleatoire pour assurer un echantillonnage suffisant
  do etape=etape_start, etape_end  ! S'arreter a la 1 est OK pour les profils de raies :-)

     if (etape==1) then
        lfixed_Rays = .true. ;  ispeed(1) = -n_speed ; ispeed(2) = n_speed
        n_rayons = 2 ! TR 1D : 1 rayon qui monte, 1 rayon qui descend
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
        fac_etape = 1.
        lprevious_converged = .false.

        ! On passe en mode mono-frequence
        deallocate(tab_speed,I0,I0c,ds,Doppler_P_x_freq)

        allocate(tab_speed(ispeed(1):ispeed(2),nb_proc), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_speed')
        tab_speed = 0.0_dp

        allocate(I0(ispeed(1):ispeed(2),nTrans_tot,n_rayons_max,nb_proc), &
             I0c(nTrans_tot,n_rayons_start,nb_proc), stat=alloc_status)
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

        if (lfixed_Rays) then
           ! On remet le generateur aleatoire a zero
           stream = 0.0
           do i=1, nb_proc
              stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
           enddo
        endif


        ! Sauvegarde des populations globales
        tab_nLevel_old = tab_nLevel
        max_n_iter_loc = 0

        ! Boucle sur les cellules

        !$omp parallel &
        !$omp default(none) &
        !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02) &
        !$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme,iv,icell,factor) &
        !$omp shared(imol,stream,n_rad,nz,n_az,n_rayons,iray_start,Doppler_P_x_freq,tab_nLevel,n_level_comp) &
        !$omp shared(deltaVmax,ispeed,r_grid,z_grid,lcompute_molRT,lkeplerian,n_cells) &
        !$omp shared(tab_speed,lfixed_Rays,lnotfixed_Rays,pop_old,pop,labs,n_speed,max_n_iter_loc,etape,pos_em_cellule) &
        !$omp shared(nTrans_tot,tab_Trans)
        !$omp do schedule(static,1)
        do icell=1, n_cells
           !$ id = omp_get_thread_num() + 1

           ! Echantillonage uniforme du profil de raie
           factor = deltaVmax(icell) / real(n_speed,kind=dp)
           if (lfixed_rays) then
              do iv=-n_speed, n_speed
                 tab_speed(iv,id) =  factor * real(iv,kind=dp)
              enddo ! iv
           endif

           if (lcompute_molRT(icell)) then

              ! Propagation des rayons
              do iray=iray_start, iray_start-1+n_rayons

                 if (etape==1) then
                    ! Position = milieu de la cellule
                    x0 = r_grid(icell)
                    y0 = 0.0_dp
                    z0 = z_grid(icell)

                    if (lkeplerian) then
                       ! Direction verticale
                       if (iray==1) then
                          w0=1.0_dp
                       else
                          w0=-1.0_dp
                       endif
                       u0 = 0.0_dp
                       v0 = 0.0_dp
                    else
                       norme = sqrt(x0*x0 + y0*y0 + z0*z0)
                       if (iray==1) then
                          u0 = x0/norme
                          v0 = y0/norme
                          w0 = z0/norme
                       else
                          u0 = -x0/norme
                          v0 = -y0/norme
                          w0 = -z0/norme
                       endif
                    endif
                 else
                    ! Position aleatoire dans la cellule
                    rand  = sprng(stream(id))
                    rand2 = sprng(stream(id))
                    rand3 = sprng(stream(id))
                    call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                    ! Direction de propagation aleatoire
                    rand = sprng(stream(id))
                    W0 = 2.0_dp * rand - 1.0_dp
                    W02 =  1.0_dp - W0*W0
                    SRW02 = sqrt(W02)
                    rand = sprng(stream(id))
                    ARGMT = PI * (2.0_dp * rand - 1.0_dp)
                    U0 = SRW02 * cos(ARGMT)
                    V0 = SRW02 * sin(ARGMT)
                 endif

                 ! Echantillonnage aleatoire du champ de vitesse
                 if (lnotfixed_Rays) then
                    do iv=ispeed(1),ispeed(2)
                       rand = sprng(stream(id)) ; tab_speed(iv,id) =  2.0_dp * (rand - 0.5_dp) * deltaVmax(icell)
                    enddo
                 endif

                 ! Integration le long du rayon
                 call integ_ray_mol(id,imol,icell,x0,y0,z0,u0,v0,w0,iray,labs, ispeed,tab_speed(:,id), &
                      nTrans_tot, tab_Trans)
              enddo ! iray


              ! Resolution de l'equilibre statistique
              n_iter_loc = 0
              pop(:,id) = tab_nLevel(icell,:)
              lconverged_loc = .false.
              ! Boucle pour converger le champ local et les populations
              ! avec champ externe fixe
              do while (.not.lconverged_loc)
                 n_iter_loc = n_iter_loc + 1

                 ! Sauvegarde ancienne pop locale
                 pop_old(:,id) = pop(:,id)

                 ! Calcul du champ de radiation
                 call J_mol_loc(id,icell,n_rayons,ispeed)  ! inclus les boucles sur Transition

                 call equilibre_rad_mol_loc(id,icell)
                 pop(:,id) = tab_nLevel(icell,:)

                 ! Critere de convergence locale
                 diff = maxval( abs(pop(1:n_level_comp,id) - pop_old(1:n_level_comp,id)) &
                      / (pop_old(1:n_level_comp,id) + 1e-30) )

                 if (diff < precision_sub) then
                    lconverged_loc = .true.
                 else
                    ! On est pas converge, on recalcule les opacites et fonctions source
                    call opacite_mol_loc(icell,imol)
                 endif

              enddo ! while : convergence champ local
              if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
           endif ! lcompute_molRT

        enddo ! icell
        !$omp end do
        !$omp end parallel
        
     	!not parallel yet
     	accelerated=.false.
     	!to be compatible with atomictransfer routines
     	popsT = transpose(tab_nLevel(:,1:n_level_comp))
     	if (lNg_acceleration .and. (n_iter > iNg_Ndelay)) then
     		  iorder = n_iter - iNg_Ndelay !local number of iterations accumulated
     		  if (ng_rest) then 
     		    write(*,*) iorder-i0_rest, " Acceleration relaxes for ", iNg_Nperiod
     		    if (iorder - i0_rest == iNg_Nperiod) ng_rest = .false.
     		  else
     		    i0_rest = iorder
                 write(*,*) " Ng iorder=",iorder, " Ndelay=", iNg_Ndelay, " Nper=", iNg_Nperiod
     		     !In the flattened array, there are:
     		     ! ilvl=1
     		     !    icell=1->ncells
     		     !                   ilvl=2
     		     !                   icell=1->Ncells
     		     flatpops=flatten(n_level_comp,n_cells,popsT)
     		     accelerated = Acceleration(Ngmol, flatpops)
     		     if (accelerated) then 
     		      popsT = reform(n_level_comp,n_cells, flatpops)
                  n_iter_accel = n_iter_accel + 1 !True number of accelerated iter
                  ng_rest = .true.
     		     endif

             endif
     		endif
        tab_nLevel(:,1:n_level_comp) = transpose(popsT)
        ! Critere de convergence totale
        maxdiff = 0.0
        do icell = 1, n_cells
           if (lcompute_molRT(icell)) then
              diff = maxval( abs( tab_nLevel(icell,1:n_level_comp) - tab_nLevel_old(icell,1:n_level_comp) ) / &
                   tab_nLevel_old(icell,1:n_level_comp) + 1e-300_dp)
              if (diff > maxdiff) maxdiff = diff
           endif
        enddo ! icell

        write(*,*) maxval(max_n_iter_loc), "sub-iterations"
        if (accelerated) then
         write(*,*) "Relative difference =", real(maxdiff), " (Accelerated)"
        else
         write(*,*) "Relative difference =", real(maxdiff)
        endif
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

              ! On continue en calculant 2 fois plus de rayons
              ! On les ajoute a l'ensemble de ceux calcules precedemment
!              iray_start = iray_start + n_rayons

           endif
        endif

        write(*,*) "STAT", minval(tab_nLevel(:,1:n_level_comp)), maxval(tab_nLevel(:,1:n_level_comp))
        call integ_tau_mol(imol)

     enddo ! while : convergence totale
  enddo ! etape

  deallocate(ds, Doppler_P_x_freq, I0, I0c)
  
  if (lNg_acceleration) then 
   CALL freeNg(Ngmol)
   deallocate(flatpops,popsT)
  endif

  return

end subroutine NLTE_mol_line_transfer

!***********************************************************

subroutine emission_line_map(imol,ibin,iaz)
  ! Creation de la carte d'emission moleculaire
  ! (ou du spectre s'il n'y a qu'un seul pixel)
  ! par ray-tracing dans une direction donnee
  ! C. Pinte
  ! 12/04/07

  implicit none

  integer, intent(in) :: imol, ibin, iaz

  real(kind=dp) :: x0,y0,z0,l, u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcorner
  real(kind=dp) :: taille_pix
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36  ! OK, ca marche avec n_rad_RT = 1000
  integer, parameter :: n_ray_star = 1000
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp) :: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, nTrans_raytracing

  integer :: n_speed_rt, n_speed_center_rt, n_extraV_rt, lambda, iv
  real :: vmax_center_rt, extra_deltaV_rt
  logical :: lresolved

 ! Direction de visee pour le ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  n_speed_rt = mol(imol)%n_speed_rt
  n_speed_center_rt = mol(imol)%n_speed_center_rt
  n_extraV_rt = mol(imol)%n_extraV_rt

  vmax_center_rt = mol(imol)%vmax_center_rt
  extra_deltaV_rt = mol(imol)%extra_deltaV_rt

  nTrans_raytracing = mol(imol)%nTrans_raytracing

  if ((ibin == 1).and.(iaz==1)) then
     allocate(I0(-n_speed_rt:n_speed_rt,nTrans_raytracing,1,nb_proc), &
          I0c(nTrans_raytracing,1,nb_proc))

     allocate(tab_speed_rt(-n_speed_rt:n_speed_rt))

     ! centre de la raie
     tab_speed_rt(-n_speed_center_rt:n_speed_center_rt) = span(-vmax_center_rt,vmax_center_rt,2*n_speed_center_rt+1)
     ! ailes de la raie
     tab_speed_rt(n_speed_center_rt+1:n_speed_rt) = indgen(n_extraV_rt) * extra_deltaV_rt + vmax_center_rt
     do i = -n_speed_rt, -n_speed_center_rt-1
        tab_speed_rt(i) = (i+n_speed_center_rt) * extra_deltaV_rt - vmax_center_rt
     enddo

     if (lorigine) then
        allocate(origine_mol(-n_speed_rt:n_speed_rt,nTrans_raytracing,n_cells,nb_proc))
        origine_mol = 0.0
     endif
  endif

  I0 = 0.0_dp
  I0c = 0.0_dp
  if (lorigine) origine_mol = 0.0

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

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Coin en bas gauche de l'image
  Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)

  if (RT_line_method == 1) then ! method 1 : echantillonanage log
     ! Pas de sous-pixel car les pixels ne sont pas carres
     n_iter_min = 1
     n_iter_max = 1

     ! dx and dy are only required for stellar map here
     taille_pix = (map_size/zoom)  ! en AU
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

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

     ! Boucle sur les rayons d'echantillonnage
     !$omp parallel &
     !$omp default(none) &
     !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
     !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
     !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,imol,ibin,iaz)
     id =1 ! pour code sequentiel

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

        do phi_RT=1,n_phi_RT ! de 0 a pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)

           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image ! C'est le centre en fait car dx = dy = 0.
           call intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

  else ! method 2 : echantillonnage lineaire avec sous-pixels
     lresolved = .true.

     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     if (l_sym_ima) then
        npix_x_max = npix_x/2 + modulo(npix_x,2)
     else
        npix_x_max = npix_x
     endif

     ! Boucle sur les pixels de l'image
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y) &
     !$omp shared(n_iter_min,n_iter_max,imol,ibin,iaz)

     id =1 ! pour code sequentiel
     n_iter_min = 1 ! 3
     n_iter_max = 1 ! 6

     !$omp do schedule(dynamic,1)
     do i = 1,npix_x_max
        !$ id = omp_get_thread_num() + 1
        do j = 1,npix_y
           !write(*,*) i,j
           ! Coin en bas gauche du pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           call intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

  endif

  ! --------------------------
  ! Ajout flux etoile
  ! --------------------------
  do i = 1, mol(imol)%nTrans_raytracing
     lambda =  mol(imol)%indice_Trans_raytracing(i) ! == iTrans
     call compute_stars_map(lambda, ibin, iaz, u, v, w, taille_pix, dx, dy, lresolved)

     do iv =  -n_speed_rt, n_speed_rt
        spectre(:,:,iv,i,ibin,iaz) = spectre(:,:,iv,i,ibin,iaz) + stars_map(:,:,1)
     enddo
     continu(:,:,i,ibin,iaz) = continu(:,:,i,ibin,iaz) + stars_map(:,:,1)
  enddo

  return

end subroutine emission_line_map

!***********************************************************

subroutine intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! Calcule l'intensite d'un pixel carre de taille, position et orientation arbitaires
  ! par une methode de Ray-tracing
  ! (u,v,w) pointe vers l'observateur
  ! Integration par methode de Romberg pour determiner le nbre de sous-pixel
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

  allocate(IP(-n_speed_rt:n_speed_rt,nTrans_raytracing), IP_old(-n_speed_rt:n_speed_rt,nTrans_raytracing), IPc(nTrans_raytracing))

  ispeed(1) = -n_speed_rt ; ispeed(2) = n_speed_rt


  labs = .false.

  ! Ray tracing : on se propage dans l'autre sens
  u0 = -u ; v0 = -v ; w0 = -w

  IP = 0.0_dp
  IPc = 0.0_dp

  ! le nbre de subpixel en x est 2^(iter-1)
  subpixels = 1
  iter = 1

  infinie : do ! Boucle infinie tant que le pixel n'est pas converge
     npix2 =  real(subpixels)**2
     IP_old = IP
     IP = 0.0_dp
     IPc = 0.0_dp

     ! Vecteurs definissant les sous-pixels
     sdx(:) = dx(:) / real(subpixels,kind=dp)
     sdy(:) = dy(:) / real(subpixels,kind=dp)

     iray = 1

     ! L'obs est en dehors de la grille
     ri = 2*n_rad ; zj=1 ; phik=1

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
              call integ_ray_mol(id,imol,icell,x0,y0,z0,u0,v0,w0,iray,labs,ispeed,tab_speed_rt, &
                   nTrans_raytracing, mol(imol)%indice_Trans_raytracing)
              ! Flux recu dans le pixel
              IP(:,:) = IP(:,:) +  I0(:,:,iray,id)
              IPc(:) = IPc(:) +  I0c(:,iray,id)
           else
              ! Il n'y a que le Cmb
              ! TODO : je ne suis pas sur de vouloir le Cmb en dehors du disque, ...
              ! IP(:,:) = IP(:,:) +  Cmb(ispeed,tab_speed)
           endif
        enddo !j
     enddo !i

     IP = IP / npix2
     IPc = IPc / npix2

     if (iter < n_iter_min) then
        ! On itere par defaut
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! On arrete pour pas tourner dans le vide
        ! write(*,*) "Warning : converging pb in ray-tracing"
        ! write(*,*) " Pixel", ipix, jpix
        exit infinie
     else
        ! On fait le test sur a difference
        diff = maxval( abs(IP - IP_old) / (IP + 1e-300_dp) )
        if (diff > precision ) then
           ! On est pas converge
           subpixels = subpixels * 2
        else
           ! On est converge
           exit infinie
        endif
     endif ! iter

     iter = iter + 1

     ! TODO : Integration Romberg
!!$     if(any(abs((I - oldintensite_pixel)) > precision * I)) then
!!$        oldintensite_pixel = I
!!$        ! Il n'y a pas un truc mieux pour utiliser les calculs a plus faible resol ??
!!$
!!$        subpixels = subpixels * 2
!!$        !if(subpixels .gt. 15) write(*,*)"large",index
!!$     else
!!$        I = (real(2**(log(real(subpixels))/log(2.d0)))*I - Oldintensite_pixel) &
!!$             /(real(2**(log(real(subpixels))/log(2.d0)))-1.d0) ! Richardson Extrapolation
!!$        ! Ok mais n'utilise que les 2 derniers calculs : il doit y avoir mieux !!!
!!$
!!$     endif
  enddo infinie

  ! Prise en compte de la surface du pixel (en sr)
  IP = IP * (pixelsize / (distance*pc_to_AU) )**2
  IPc = IPc * (pixelsize / (distance*pc_to_AU) )**2

  ! et multiplication par la frequence pour avoir du nu.F_nu

  ! Warning IP, IPc are smaller array (dimension mol(imol)%nTrans_raytracin)
  do i=1,mol(imol)%nTrans_raytracing
     iTrans = mol(imol)%indice_Trans_raytracing(i)
     IP(:,i) = IP(:,i) * transfreq(iTrans)
     IPc(i) = IPc(i) * transfreq(iTrans)
  enddo
  ! Unite teste OK pour le Cmb
  ! profil de raie non convolue teste ok avec torus

  if (RT_line_method==1) then ! Sommation implicite sur les pixels
     spectre(1,1,:,:,ibin,iaz) = spectre(1,1,:,:,ibin,iaz) + IP(:,:)
     continu(1,1,:,ibin,iaz) = continu(1,1,:,ibin,iaz) + IPc(:)
  else
     spectre(ipix,jpix,:,:,ibin,iaz) = IP(:,:)
     continu(ipix,jpix,:,ibin,iaz) = IPc(:)
  endif

  return

end subroutine intensite_pixel_mol

!***********************************************************

subroutine init_dust_mol(imol)
  ! calcul les opacites et emissivites des cellules
  ! aux longueurs d'onde des raies moleculaires
  ! pour prendre en compte l'interaction radiative entre
  ! les 2 phases.
  ! C. Pinte
  ! 17/10/07
  ! TODO : gerer le cas ou l'albedo (et donc scattering) non negligeable

  use mem, only : realloc_dust_mol, clean_mem_dust_mol

  implicit none

  integer, intent(in) :: imol
  integer :: iTrans, p_lambda, icell
  real(kind=dp) :: freq!, Jnu
  real :: T, wl, kap

  real(kind=dp) :: cst_E

  real, parameter :: gas_dust = 100
  real, parameter :: delta_lambda = 0.025

  cst_E=2.0*hp*c_light**2

  ! Reallocation des tableaux de proprietes de poussiere
  ! n_lambda =   mol(imol)%nTrans_raytracing ! opacites dust considerees cst sur le profil de raie
  n_lambda = nTrans_tot ! opacites dust considerees cst sur le profil de raie

  ! On n'est interesse que par les prop d'abs : pas besoin des matrices de mueller
  ! -> pas de polarisation, on utilise une HG
  scattering_method=1 ; lscattering_method1 = .true. ; p_lambda = 1
  aniso_method = 2 ; lmethod_aniso1 = .false.

  lsepar_pola = .false.
  ltemp = .false.
  lmono = .true. ! equivalent au mode sed2

  call realloc_dust_mol()

  if (ldust_mol) then
     ! Tableau de longeur d'onde
     do iTrans=1,nTrans_tot
        tab_lambda(iTrans) = c_light/Transfreq(iTrans) * 1.0e6 ! en microns
        tab_lambda_sup(iTrans)= tab_lambda(iTrans)*delta_lambda
        tab_lambda_inf(iTrans)= tab_lambda(iTrans)/delta_lambda
        tab_delta_lambda(iTrans) = tab_lambda_sup(iTrans) - tab_lambda_inf(iTrans)
     enddo

     if (lbenchmark_water3) then ! opacite en loi de puissance
        write(*,*) "WARNING : hard-coded gas_dust =", gas_dust

        do iTrans=1,nTrans_tot
           wl = tab_lambda(iTrans)

           ! Loi d'opacite (cm^2 par g de poussiere)
           if (wl > 250) then
              kap = 10. * (wl/250.)**(-2.0)
           else
              kap = 10. * (wl/250.)**(-1.3)
           endif

           ! Multiplication par densite
           ! AU_to_cm**2 car on veut kappa_abs_LTE en AU-1
           do icell=1,n_cells
              kappa_abs_LTE(icell,iTrans) =  kap * (densite_gaz(icell) * cm_to_m**3) * masse_mol_gaz / &
                   gas_dust / cm_to_AU
           enddo

        enddo ! iTrans

        ! Pas de scattering
        kappa(:,:) = kappa_abs_LTE(:,:)

     else ! cas par defaut
        call init_indices_optiques()

        ! On recalcule les proprietes optiques
        write(*,*) "Computing dust properties for", nTrans_tot, "wavelength"
        do iTrans=1,nTrans_tot
           call prop_grains(iTrans)
           call opacite(iTrans, iTrans, no_scatt=.true.)
        enddo
     endif


     ! Changement d'unite : kappa en m-1 pour le TR dans les raies !!!!!!!
     ! Sera reconverti en AU-1 dans opacite_mol_loc
     !kappa = kappa * m_to_AU
     !kappa_abs_eg = kappa_abs_eg * m_to_AU

     ! calcul de l'emissivite de la poussiere
     do iTrans=1,nTrans_tot
        freq = Transfreq(iTrans)

        ! TODO : accelerer cette boucle via routine Bnu_disk (ca prend du tps ???)
        ! TODO : generaliser pour tous les types de grains (ca doit deja exister non ???)
        ! TODO : ca peut aussi servir pour faire du ray-tracing ds le continu 8-)
        do icell=1, n_cells
           !-- ! Interpolation champ de radiation en longeur d'onde
           !-- if (lProDiMo2mcfost) then
           !--    Jnu = interp(m2p%Jnu(ri,zj,:), m2p%wavelengths, real(tab_lambda(iTrans)))
           !-- else
           !--    Jnu = 0.0 ! todo : pour prendre en compte scattering
           !-- endif

           T = Tdust(icell)
           ! On ne fait que du scattering isotropique dans les raies pour le moment ...
           emissivite_dust(icell,iTrans) = kappa_abs_LTE(icell,iTrans) * Bnu(freq,T) ! + kappa_sca(iTrans,ri,zj,phik) * Jnu
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
  integer :: icell

  ldust_mol  = .true.
  lkeplerian = .true.

  if (lfirst_time) then
     lfirst_time = .false.
     call init_Tgas()

     ! Velocity field in  m.s-1
     ! Warning : assume all stars are at the center of the disk
     if (.not.lVoronoi) then ! Velocities are defined from SPH files in Voronoi mode
        if (lcylindrical_rotation) then ! Midplane Keplerian velocity
           do icell=1, n_cells
              vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
           enddo
        else ! dependance en z
           do icell=1, n_cells
              vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg * r_grid(icell)**2 / &
                   ((r_grid(icell)**2 + z_grid(icell)**2)**1.5 * AU_to_m) )
           enddo
        endif
     endif
     v_turb = vitesse_turb
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
     ! Temperature gaz = poussiere
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
          (densite_gaz(icell) > tiny_real) .and. (Tcin(icell) > 1.)
  enddo

  return

end subroutine init_abundance

!***********************************************************


end module mol_transfer
