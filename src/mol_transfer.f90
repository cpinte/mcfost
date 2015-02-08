module mol_transfer

  use parametres
  use molecular_emission
  use naleat
  use disk
  use resultats
  use utils
  use opacity
  use molecules
  !$ use omp_lib

  use input
  use benchmarks
  use output
  use molecules
  use dust
  use scattering
  use optical_depth
  use ProDiMo, only: read_ProDiMo2mcfost
  use dust_ray_tracing, only: init_directions_ray_tracing
  use dust_transfer, only : compute_stars_map
  use stars

  implicit none

  contains

subroutine mol_line_transfer()

  implicit none

  real(kind=db) :: u0, v0, w0
  integer :: iTrans, imol, ibin, iaz

  if (lProDiMo2mcfost) ldust_mol = .true.

  ! Liberation memoire
  call dealloc_em_th()

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

     nTrans = nTrans_tot
     allocate(indice_Trans(nTrans_tot))
     do iTrans=1,nTrans
        indice_Trans(iTrans) = iTrans
     enddo

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

     ! Freeze out eventuel
     if (lfreeze_out) call freeze_out()

     if (lProDiMo2mcfost) call read_ProDiMo2mcfost(imol)

     ! Absorption et emissivite poussiere
     call init_dust_mol(imol)

     ! recalcul des flux stellaires aux nouvelles longeurs d'onde
     call repartition_energie_etoiles()

     call init_Doppler_profiles(imol)  ! ne depend pas de nTrans

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
     if (.not.lmol_LTE) call NLTE_mol_line_transfer(imol)

     call ecriture_pops(imol)
     call ecriture_Tex(imol)

     !--- Creation carte emission moleculaire : ray-tracing
     if (mol(imol)%lline) then
        ! Selection des transitions
        if (mol(imol)%nTrans_raytracing /= nTrans_tot) then
           deallocate(indice_Trans)
           nTrans = mol(imol)%nTrans_raytracing
           allocate(indice_Trans(mol(imol)%nTrans_raytracing))
           indice_Trans(:) = mol(imol)%indice_Trans_raytracing(:)
        endif

        u0 = sin(angle_interet/180._db*pi) ;  v0=0.0_db ; w0 = sqrt(1.0_db - u0*u0)

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
  integer, parameter :: n_speed3 = 1
  integer, parameter :: n_rayons_max = n_rayons_start2 * (2**(n_iter2_max-1))
  integer :: n_level_comp
  real, parameter :: precision_sub = 1.0e-3
  real, parameter :: precision = 1.0e-1

  integer :: etape, etape_start, etape_end, ri, zj, phik, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, iv, n_speed
  integer, dimension(nb_proc) :: max_n_iter_loc

  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged

  real :: rand, rand2, rand3, fac_etape

  real(kind=db) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, maxdiff, norme

  real(kind=db), dimension(nLevels,nb_proc)  :: pop, pop_old

  logical :: labs

  integer, dimension(2) :: ispeed
  real(kind=db), dimension(:,:), allocatable :: tab_speed


  labs = .true.

  id = 1

  n_speed = mol(imol)%n_speed_rt ! j'utilise le meme maintenant
  n_level_comp = min(mol(imol)%iLevel_max,nLevels)

  if (n_level_comp < 2) then
     write(*,*) "ERROR : n_level_comp must be > 2"
     write(*,*) "Exiting."
     stop
  endif
  write(*,*) "NLTE line transfer on", n_level_comp, "levels"

  etape_start=1
  if (lprecise_pop) then
     etape_end = 3
  else
     etape_end = 1  ! 2
  endif

  allocate(tab_speed(-n_speed:n_speed,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_speed'
     stop
  endif
  tab_speed = 0.0_db

  allocate(I0(-n_speed:n_speed,nTrans,n_rayons_start,nb_proc), &
       I0c(nTrans,n_rayons_start,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error I0'
     stop
  endif
  I0 = 0.0_db
  I0c = 0.0_db

  if (ldouble_RT) then
     allocate(I02(-n_speed:n_speed,nTrans,n_rayons_start,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error I0'
        stop
     endif
     I02 = 0.0_db
  endif

  allocate(ds(n_rayons_start,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error ds'
     stop
  endif
  ds = 0.0_db

  allocate(Doppler_P_x_freq(-n_speed:n_speed,n_rayons_start,nb_proc), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error Doppler_P_x_freq'
     stop
  endif
  Doppler_P_x_freq = 0.0_db

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
        lfixed_Rays = .false.;  ispeed(1) = 1 ; ispeed(2) = n_speed3
        n_rayons = n_rayons_start2
        fac_etape = 1.
        lprevious_converged = .false.

        ! On passe en mode mono-frequence
        deallocate(tab_speed,I0,I0c,ds,Doppler_P_x_freq)

        allocate(tab_speed(ispeed(1):ispeed(2),nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_speed'
           stop
        endif
        tab_speed = 0.0_db

        allocate(I0(ispeed(1):ispeed(2),nTrans,n_rayons_max,nb_proc), &
             I0c(nTrans,n_rayons_start,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error I0'
           stop
        endif
        I0 = 0.0_db
        I0c = 0.0_db

        allocate(ds(n_rayons_max,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error ds'
           stop
        endif
        ds = 0.0_db

        allocate(Doppler_P_x_freq(ispeed(1):ispeed(2),n_rayons_max,nb_proc), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error Doppler_P_x_freq'
           stop
        endif
        Doppler_P_x_freq = 0.0_db
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
        !$omp private(id,ri,zj,phik,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02) &
        !$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme,iv) &
        !$omp shared(imol,stream,n_rad,nz,n_az,n_rayons,iray_start,Doppler_P_x_freq,tab_nLevel,n_level_comp) &
        !$omp shared(tab_deltaV,deltaVmax,ispeed,r_grid,z_grid,lcompute_molRT,lkeplerian) &
        !$omp shared(tab_speed,lfixed_Rays,lnotfixed_Rays,pop_old,pop,labs,n_speed,max_n_iter_loc,etape)
        !$omp do schedule(static,1)
        do ri=1, n_rad
           !$ id = omp_get_thread_num() + 1
           do zj=1, nz

              do phik=1, n_az
                 ! Echantillonage uniforme du profil de raie
                 if (lfixed_rays) then
                    tab_speed(:,id) = tab_deltaV(:,ri,zj,phik)
                 endif

                 if (lcompute_molRT(ri,zj,phik)) then

                    ! Propagation des rayons
                    do iray=iray_start, iray_start-1+n_rayons

                       if (etape==1) then
                          ! Position = milieu de la cellule
                          x0 = r_grid(ri,zj)
                          y0 = 0.0_db
                          z0 = z_grid(ri,zj)

                          if (lkeplerian) then
                             ! Direction verticale
                             if (iray==1) then
                                w0=1.0_db
                             else
                                w0=-1.0_db
                             endif
                             u0 = 0.0_db
                             v0 = 0.0_db
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
                          call  pos_em_cellule(ri,zj,phik,rand,rand2,rand3,x0,y0,z0)

                          ! Direction de propagation aleatoire
                          rand = sprng(stream(id))
                          W0 = 2.0_db * rand - 1.0_db
                          W02 =  1.0_db - W0*W0
                          SRW02 = sqrt(W02)
                          rand = sprng(stream(id))
                          ARGMT = PI * (2.0_db * rand - 1.0_db)
                          U0 = SRW02 * cos(ARGMT)
                          V0 = SRW02 * sin(ARGMT)
                       endif

                       ! Echantillonnage aleatoire du champ de vitesse
                       if (lnotfixed_Rays) then
                          do iv=ispeed(1),ispeed(2)
                             !tab_speed(1,id) = gauss_random(id) * deltaVmax(ri,zj)
                             rand = sprng(stream(id)) ; tab_speed(iv,id) =  2.0_db * (rand - 0.5_db) * deltaVmax(ri,zj,phik)
                          enddo
                       endif


                       ! Integration le long du rayon
                       call integ_ray_mol(id,ri,zj,phik,x0,y0,z0,u0,v0,w0,iray,labs,ispeed,tab_speed(:,id))

                    enddo ! iray


                    ! Resolution de l'equilibre statistique
                    n_iter_loc = 0
                    pop(:,id) = tab_nLevel(ri,zj,phik,:)
                    lconverged_loc = .false.
                    ! Boucle pour converger le champ local et les populations
                    ! avec champ externe fixe
                    do while (.not.lconverged_loc)
                       n_iter_loc = n_iter_loc + 1

                       ! Sauvegarde ancienne pop locale
                       pop_old(:,id) = pop(:,id)

                       ! Calcul du champ de radiation
                       call J_mol_loc(id,ri,zj,phik,n_rayons,ispeed)  ! inclus les boucles sur Transition

                       call equilibre_rad_mol_loc(id,ri,zj,phik)
                       pop(:,id) = tab_nLevel(ri,zj,phik,:)

                       ! Critere de convergence locale
                       diff = maxval( abs(pop(1:n_level_comp,id) - pop_old(1:n_level_comp,id)) &
                            / (pop_old(1:n_level_comp,id) + 1e-30) )

                       if (diff < precision_sub) then
                          lconverged_loc = .true.
                       else
                          ! On est pas converge, on recalcule les opacites et fonctions source
                          call opacite_mol_loc(ri,zj,phik,imol)
                       endif

                    enddo ! while : convergence champ local
                    if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc

                 endif ! lcompute_molRT

              enddo ! phik
           enddo !zj
        enddo !ri
        !$omp end do
        !$omp end parallel

        phik = 1 ! pas 3D
        ! Critere de convergence totale
        maxdiff = 0.0
        do ri=1,n_rad
           do zj=1,nz
              if (lcompute_molRT(ri,zj,phik)) then
                 diff = maxval( abs( tab_nLevel(ri,zj,phik,1:n_level_comp) - tab_nLevel_old(ri,zj,phik,1:n_level_comp) ) / &
                      tab_nLevel_old(ri,zj,phik,1:n_level_comp) + 1e-300_db)

             !    write(*,*) abs(tab_nLevel(ri,zj,1:n_level_comp) - tab_nLevel_old(ri,zj,1:n_level_comp)) / &
              !        tab_nLevel_old(ri,zj,1:n_level_comp)
               !  write(*,*) ri, zj, diff, tab_nLevel_old(ri,zj,1:n_level_comp)
                 if (diff > maxdiff) maxdiff = diff
              endif
           enddo
        enddo

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

              ! On continue en calculant 2 fois plus de rayons
              ! On les ajoute a l'ensemble de ceux calcules precedemment
!              iray_start = iray_start + n_rayons

           endif
        endif

        write(*,*) "STAT", minval(tab_nLevel(:,:,:,1:n_level_comp)), maxval(tab_nLevel(:,:,:,1:n_level_comp))
        call integ_tau_mol(imol)

     enddo ! while : convergence totale
  enddo ! etape

  deallocate(ds, Doppler_P_x_freq, I0, I0c)

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

  real(kind=db) :: x0,y0,z0,l, uv, u,v,w

  real(kind=db), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=db), dimension(3,nb_proc) :: pixelcorner
  real(kind=db) :: taille_pix
  integer :: i,j, id, igridx_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36  ! OK, ca marche avec n_rad_RT = 1000
  integer, parameter :: n_ray_star = 1000
  real(kind=db), dimension(n_rad_RT) :: tab_r
  real(kind=db) :: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, nTrans_raytracing

  integer :: n_speed_rt, n_speed_center_rt, n_extraV_rt, lambda, iv
  real :: vmax_center_rt, extra_deltaV_rt

 ! Direction de visee pour le ray-tracing
  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin) ;
  uvw = (/u,v,w/)

  n_speed_rt = mol(imol)%n_speed_rt
  n_speed_center_rt = mol(imol)%n_speed_center_rt
  n_extraV_rt = mol(imol)%n_extraV_rt

  vmax_center_rt = mol(imol)%vmax_center_rt
  extra_deltaV_rt = mol(imol)%extra_deltaV_rt

  nTrans_raytracing = mol(imol)%nTrans_raytracing

  if (ibin == 1) then
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
        allocate(origine_mol(-n_speed_rt:n_speed_rt,nTrans_raytracing,n_rad,nz,nb_proc))
        origine_mol = 0.0
     endif
  endif

  I0 = 0.0_db
  I0c = 0.0_db
  if (lorigine) origine_mol = 0.0

  ! Definition des vecteurs de base du plan image dans le repere universel

  ! Definition des vecteurs de base du plan image dans le repere universel

  ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/sin(tab_RT_az(iaz) * deg_to_rad),-cos(tab_RT_az(iaz) * deg_to_rad),0/)

  ! Vecteur x image avec PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
  y_plan_image = cross_product(x_plan_image, uvw)

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

     dx(:) = 0.0_db
     dy(:) = 0.0_db
     i = 1
     j = 1

     rmin_RT = max(w*0.9_db,0.05_db) * rmin
     rmax_RT = 2.0_db * Rmax

     tab_r(1) = rmin_RT
     fact_r = exp( (1.0_db/(real(n_rad_RT,kind=db) -1))*log(rmax_RT/rmin_RT) )

     do ri_RT = 2, n_rad_RT
        tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
     enddo

     fact_A = sqrt(pi * (fact_r - 1.0_db/fact_r)  / n_phi_RT )


     ! Boucle sur les rayons d'echantillonnage
     !$omp parallel &
     !$omp default(none) &
     !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
     !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
     !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,imol,ibin,iaz)
     id =1 ! pour code sequentiel

     if (l_sym_ima) then
        cst_phi = pi  / real(n_phi_RT,kind=db)
     else
        cst_phi = deux_pi  / real(n_phi_RT,kind=db)
     endif


     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1

        r = tab_r(ri_RT)
        taille_pix =  fact_A * r ! racine carree de l'aire du pixel

        do phi_RT=1,n_phi_RT ! de 0 a pi
           phi = cst_phi * (real(phi_RT,kind=db) -0.5_db)

           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image ! C'est le centre en fait car dx = dy = 0.
           call intensite_pixel_mol(id,imol,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel

  else ! method 2 : echantillonnage lineaire avec sous-pixels

     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = map_size / real(max(igridx,igridy),kind=db) ! en AU
     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     if (l_sym_ima) then
        igridx_max = igridx/2 + modulo(igridx,2)
     else
        igridx_max = igridx
     endif

     ! Boucle sur les pixels de l'image
     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,pixelcorner,dx,dy,u,v,w,taille_pix,igridx_max,igridy) &
     !$omp shared(n_iter_min,n_iter_max,imol,ibin,iaz)

     id =1 ! pour code sequentiel
     n_iter_min = 1 ! 3
     n_iter_max = 1 ! 6

     !$omp do schedule(dynamic,1)
     do i = 1,igridx_max
        !$ id = omp_get_thread_num() + 1
        do j = 1,igridy
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
  do lambda = 1, mol(imol)%nTrans_raytracing
     call compute_stars_map(lambda,ibin, u, v, w)

     do iv =  -n_speed_rt, n_speed_rt
        spectre(:,:,iv,lambda,ibin,iaz) = spectre(:,:,iv,lambda,ibin,iaz) + stars_map(:,:)
     enddo
     continu(:,:,lambda,ibin,iaz) = continu(:,:,lambda,ibin,iaz) + stars_map(:,:)
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
  real(kind=db), dimension(3), intent(in) :: pixelcorner,dx,dy
  real(kind=db), intent(in) :: pixelsize,u,v,w
  real(kind=db), dimension(:,:), allocatable :: IP, IP_old
  real(kind=db), dimension(:), allocatable :: IPc

  integer, parameter :: maxSubPixels = 32

  real(kind=db) :: x0,y0,z0,u0,v0,w0
  real(kind=db), dimension(3) :: sdx, sdy
  real :: npix2, diff

  real, parameter :: precision = 1.e-2
  integer :: i, j, subpixels, iray, ri, zj, phik, iTrans, iiTrans, iter, n_speed_rt, nTrans_raytracing

  logical :: lintersect, labs

  integer, dimension(2) :: ispeed

  n_speed_rt = mol(imol)%n_speed_rt
  nTrans_raytracing = mol(imol)%nTrans_raytracing

  allocate(IP(-n_speed_rt:n_speed_rt,nTrans_raytracing), IP_old(-n_speed_rt:n_speed_rt,nTrans_raytracing), IPc(nTrans_raytracing))

  ispeed(1) = -n_speed_rt ; ispeed(2) = n_speed_rt

  labs = .false.

  ! Ray tracing : on se propage dans l'autre sens
  u0 = -u ; v0 = -v ; w0 = -w

  IP = 0.0_db
  IPc = 0.0_db

  ! le nbre de subpixel en x est 2^(iter-1)
  subpixels = 1
  iter = 1

  infinie : do ! Boucle infinie tant que le pixel n'est pas converge
     npix2 =  real(subpixels)**2
     IP_old = IP
     IP = 0.0_db
     IPc = 0.0_db

     ! Vecteurs definissant les sous-pixels
     sdx(:) = dx(:) / real(subpixels,kind=db)
     sdy(:) = dy(:) / real(subpixels,kind=db)

     iray = 1

     ! L'obs est en dehors de la grille
     ri = 2*n_rad ; zj=1 ; phik=1

     ! Boucle sur les sous-pixels qui calcule l'intensite au centre
     ! de chaque sous pixel
     do i = 1,subpixels
        do j = 1,subpixels
           ! Centre du sous-pixel
           x0 = pixelcorner(1) + (i - 0.5_db) * sdx(1) + (j-0.5_db) * sdy(1)
           y0 = pixelcorner(2) + (i - 0.5_db) * sdx(2) + (j-0.5_db) * sdy(2)
           z0 = pixelcorner(3) + (i - 0.5_db) * sdx(3) + (j-0.5_db) * sdy(3)

           ! On se met au bord de la grille : propagation a l'envers
           call move_to_grid(x0,y0,z0,u0,v0,w0,ri,zj,phik,lintersect)

           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
              call integ_ray_mol(id,ri,zj,phik,x0,y0,z0,u0,v0,w0,iray,labs,ispeed,tab_speed_rt)
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
        diff = maxval( abs(IP - IP_old) / (IP + 1e-300_db) )
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
  do iTrans=1,nTrans
     iiTrans = indice_Trans(iTrans)
     IP(:,iTrans) = IP(:,iTrans) * transfreq(iiTrans)
     IPc(iTrans) = IPc(iTrans) * transfreq(iiTrans)
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

end module mol_transfer
