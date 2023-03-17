module optical_depth

  use parametres
  use dust_prop
  use constantes
  use molecular_emission
  use utils
  use dust_ray_tracing
  use grid
  use cylindrical_grid
  use radiation_field, only : save_radiation_field
  use density
  use stars, only : intersect_stars, star_rad
  use opacity_atom, only : opacity_atom_bb_loc, contopac_atom_loc, Itot, psi

  implicit none


  contains

subroutine physical_length(id,lambda,p_lambda,Stokes,icell,xio,yio,zio,u,v,w,flag_star,flag_direct_star,&
     extrin,ltot,flag_sortie,lpacket_alive)
  ! Integration par calcul de la position de l'interface entre cellules
  ! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
  ! C. Pinte
  ! 05/02/05

  implicit none

  integer, intent(in) :: id,lambda, p_lambda
  integer, intent(inout) :: icell
  real(kind=dp), dimension(4), intent(in) :: Stokes
  logical, intent(in) :: flag_star, flag_direct_star
  real(kind=dp), intent(inout) :: u,v,w
  real, intent(in) :: extrin
  real(kind=dp), intent(inout) :: xio,yio,zio
  real, intent(out) :: ltot
  logical, intent(out) :: flag_sortie
  logical, intent(inout) :: lpacket_alive

  real(kind=dp) :: x0, y0, z0, x1, y1, z1, x_old, y_old, z_old, extr
  real(kind=dp) :: l, tau, opacite, l_contrib, l_void_before
  integer :: icell_in, icell_old, next_cell, previous_cell, icell_star, i_star
  integer, target :: icell0

  logical :: lcellule_non_vide, lstop, lintersect_stars

  integer, pointer :: p_icell

  lstop = .false.
  flag_sortie = .false.

  x0=xio;y0=yio;z0=zio
  x1=xio;y1=yio;z1=zio

  extr=extrin
  icell_in = icell

  next_cell = icell
  icell0 = 0 ! to define previous_cell

  ltot = 0.0

  ! Calcule les angles de diffusion pour la direction de propagation donnee
  if ((.not.letape_th).and.lscatt_ray_tracing1) call angles_scatt_rt1(id,u,v,w)

  ! Will the packet intersect a star
  call intersect_stars(x0,y0,z0, u,v,w, lintersect_stars, i_star, icell_star)

  if (lvariable_dust) then
     p_icell => icell0
  else
     p_icell => icell_ref
  endif

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     icell_old = icell0
     x_old = x0 ; y_old = y0 ; z_old = z0

     x0=x1 ; y0=y1 ; z0=z1
     previous_cell = icell0
     icell0 = next_cell

     ! Test sortie
     if (test_exit_grid(icell0, x0, y0, z0)) then
        flag_sortie = .true.
        return
     endif
     if (lintersect_stars) then
        if (icell0 == icell_star) then
           lpacket_alive = .false.
           flag_sortie = .true.
           return
        endif
     endif

     ! Pour cas avec approximation de diffusion
     if (icell0 <= n_cells) then
        lcellule_non_vide=.true.
        opacite = kappa(p_icell,lambda) * kappa_factor(icell0)

        if (l_dark_zone(icell0)) then
           ! On renvoie le paquet dans l'autre sens
           u = -u ; v = -v ; w=-w
           ! et on le renvoie au point de depart
           icell = icell_old
           xio = x_old ; yio = y_old ; zio = z_old
           flag_sortie= .false.
           return
        endif
     else
        lcellule_non_vide=.false.
        opacite = 0.0_dp
     endif

     ! Calcul longeur de vol et profondeur optique dans la cellule
     call cross_cell(x0,y0,z0, u,v,w,  icell0, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     ! opacity wall
     !---if (ri0 == 1) then
     !---   ! Variation de hauteur du mur en cos(phi/2)
     !---   phi = atan2(y0,x0)
     !---   hh = h_wall * abs(cos(phi/2.))
     !---   hhm = -h_wall * abs(cos((phi+pi)/2.))
     !---
     !---   ! Ajout de l'opacite du mur le cas echeant
     !---   if ((z0 <= hh).and.(z0 >= hhm)) then
     !---      opacite = opacite + kappa_wall
     !---   endif
     !---endif

     tau = l_contrib * opacite ! opacite constante dans la cellule

     ! Comparaison integrale avec tau
     ! et ajustement longueur de vol eventuellement
     if(tau > extr) then ! On a fini d'integrer
        lstop = .true.
        l_contrib = l_contrib * (extr/tau) ! on rescale l_contrib pour que tau=extr et on ajoute la longeur de vol dans le vide
        l = l_void_before + l_contrib
        ltot=ltot+l
     else ! Il reste extr - tau a integrer dans la cellule suivante
        extr=extr-tau
        ltot=ltot+l
     endif

     ! Stockage des champs de radiation
     if (lcellule_non_vide) call save_radiation_field(id,lambda,p_lambda, icell0, Stokes, l_contrib, &
          x0,y0,z0, x1,y1,z1, u,v,w, flag_star, flag_direct_star)

     ! On a fini d'integrer : sortie de la routine
     if (lstop) then
        flag_sortie = .false.
        xio=x0+l*u
        yio=y0+l*v
        zio=z0+l*w

        icell = icell0

        ! TODO : here
        if (.not.lVoronoi) then
           if (l3D) then
              if (lcylindrical) call indice_cellule(xio,yio,zio, icell)
            ! following lines are useless --> icell0 is not returned
           !else
           !   if (lcylindrical) then
           !      call verif_cell_position_cyl(icell0, xio, yio, zio)
           !   else if (lspherical) then
           !      call verif_cell_position_sph(icell0, xio, yio, zio)
           !   endif
           endif
        endif ! todo : on ne fait rien dans la cas Voronoi ???

        return
     endif ! lstop

  enddo ! boucle infinie
  write(*,*) "BUG"
  return

end subroutine physical_length

!********************************************************************

subroutine integ_tau(lambda)

  implicit none

  integer, intent(in) :: lambda

  integer :: icell!, i

  real(kind=dp), dimension(4) :: Stokes
  ! angle de visee en deg
  real :: angle
  real(kind=dp) :: x0, y0, z0, u0, v0, w0
  real :: tau
  real(kind=dp) :: lmin, lmax

  angle=angle_interet

  x0=0.0 ; y0=0.0 ; z0=0.0
  Stokes = 0.0_dp ; Stokes(1) = 1.0_dp
  w0 = 0.0 ; u0 = 1.0 ; v0 = 0.0

  call indice_cellule(x0,y0,z0, icell)
  call optical_length_tot(1,lambda,Stokes,icell,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)

  !tau = 0.0
  !do i=1, n_rad
  !   tau=tau+kappa(cell_map(i,1,1),lambda)*(r_lim(i)-r_lim(i-1))
  !enddo
  write(*,*) 'Integ tau dans plan eq. = ', tau

  if (.not.lvariable_dust) then
     icell = icell_not_empty
     if (kappa(icell_ref,lambda) * kappa_factor(icell) > tiny_real) then
        write(*,*) " Column density (g/cm^2)   = ", real(tau*(masse(icell)/(volume(icell)*AU_to_cm**3))/ &
             (kappa(icell_ref,lambda) * kappa_factor(icell)/AU_to_cm))
     endif
  endif

  Stokes = 0.0_dp ; Stokes(1) = 1.0_dp
  w0 = cos((angle)*pi/180.)
  u0 = sqrt(1.0-w0*w0)
  v0 = 0.0

  call indice_cellule(x0,y0,z0, icell)
  call optical_length_tot(1,lambda,Stokes,icell,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)

  write(*,fmt='(" Integ tau (i =",f4.1," deg)   = ",E12.5)') angle, tau

  if (.not.lvariable_dust) then
     icell = icell_not_empty
     if (kappa(icell_ref,lambda) * kappa_factor(icell) > tiny_real) then
        write(*,*) " Column density (g/cm^2)   = ", real(tau*(masse(icell)/(volume(icell)*AU_to_cm**3))/ &
             (kappa(icell_ref,lambda) * kappa_factor(icell)/AU_to_cm))
     endif
  endif

  return

end subroutine integ_tau

!***********************************************************

subroutine optical_length_tot(id,lambda,Stokes,icell,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
! Integration par calcul de la position de l'interface entre cellules
! de l'opacite totale dans une direction donnée
! Grille a geometrie cylindrique
! C. Pinte
! 19/04/05

  implicit none

  integer, intent(in) :: id,lambda, icell
  real(kind=dp),dimension(4), intent(in) :: Stokes
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: xi,yi,zi
  real, intent(out) :: tau_tot_out
  real(kind=dp), intent(out) :: lmin,lmax


  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, ltot, tau, opacite, tau_tot, correct_plus, correct_moins, l_contrib, l_void_before
  integer :: previous_cell, next_cell
  integer, target :: icell0
  integer, pointer :: p_icell

  correct_plus = 1.0_dp + prec_grille
  correct_moins = 1.0_dp - prec_grille

  x1=xi;y1=yi;z1=zi

  tau_tot=0.0_dp

  lmin=0.0_dp
  ltot=0.0_dp

  next_cell = icell
  icell0 = 0 ! for previous_cell, just for Voronoi

  if (lvariable_dust) then
     p_icell => icell0
  else
     p_icell => icell_ref
  endif

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     previous_cell = icell0
     icell0 = next_cell
     x0=x1;y0=y1;z0=z1

     ! Test sortie
     if (test_exit_grid(icell0, x0, y0, z0)) then
        tau_tot_out=tau_tot
        lmax=ltot
        return
     endif

     if (icell0 <= n_cells) then
        opacite = kappa(p_icell,lambda) * kappa_factor(icell0)
     else
        opacite = 0.0_dp
     endif

     ! Calcul longeur de vol et profondeur optique dans la cellule
     call cross_cell(x0,y0,z0, u,v,w,  icell0, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     tau=l_contrib*opacite ! opacite constante dans la cellule

     tau_tot = tau_tot + tau
     ltot= ltot + l

     if (tau_tot < tiny_real) lmin=ltot

  enddo ! boucle infinie

  write(*,*) "BUG"
  return

end subroutine optical_length_tot

!************************************************************

subroutine compute_column(type, column, lambda)

  use density, only : densite_gaz
  !$ use omp_lib

  integer, intent(in) ::type ! 1 = column_density, 2 = optical_depth
  integer, intent(in), optional :: lambda

  integer, parameter :: n_directions = 4
  real, dimension(n_cells,n_directions), intent(out) :: column

  integer :: icell, next_cell, previous_cell, direction, icell0, p_icell

  real(kind=dp) :: x0,y0,z0, x1,y1,z1, norme, l, u,v,w, l_contrib, l_void_before, CD_units, factor, sum

  if (type==1) then
     CD_units = AU_to_m * masse_mol_gaz / (m_to_cm)**2 ! g/cm^-2 and AU_to_m factor as l_contrib is in AU
  else if (type==3) then
     CD_units = AU_to_m / (m_to_cm)**2 ! particle/cm^-2 and AU_to_m factor as l_contrib is in AU
  endif

  column(:,:) = 0.0
  do direction = 1, n_directions
     !$omp parallel default(none) &
     !$omp shared(densite_gaz,tab_abundance,lVoronoi,Voronoi,direction,column,r_grid,z_grid,phi_grid,n_cells,cross_cell) &
     !$omp shared(CD_units,kappa,kappa_factor,lambda,type,test_exit_grid,icell_ref,lvariable_dust) &
     !$omp private(icell,previous_cell,next_cell,icell0,p_icell,x0,y0,z0,x1,y1,z1,norme,u,v,w,l,l_contrib,l_void_before,factor,sum)
     p_icell = icell_ref
     !$omp do
     do icell=1,n_cells
        if (lVoronoi) then
           x1 = Voronoi(icell)%xyz(1)
           y1 = Voronoi(icell)%xyz(2)
           z1 = Voronoi(icell)%xyz(3)
        else
           x1 = r_grid(icell) * cos(phi_grid(icell))
           y1 = r_grid(icell) * sin(phi_grid(icell))
           z1 = z_grid(icell)
        endif

        if (direction == 1) then ! to star (assumed to in 0,0,0 for now + only 1 star)
           norme = 1./sqrt(x1*x1 + y1*y1 + z1*z1)
           u  = -x1 * norme ; v = -y1 * norme ; w = -z1 * norme
        else if (direction == 2) then ! vertical +z
           u = 0.0 ; v = 0.0 ; w = 1.0
        else if (direction == 3) then ! vertical -z
           u = 0.0 ; v = 0.0 ; w = -1.0
        else ! radial
           u = x1 ; v = y1 ; w = 0
           norme = 1./sqrt(u**2 + v**2)
           u = u * norme ; v = v* norme
        endif

        next_cell = icell
        icell0 = 0

        sum = 0.
        inf_loop: do
           if (lvariable_dust) p_icell = icell0
           previous_cell = icell0
           icell0 = next_cell
           x0 = x1 ; y0 = y1 ; z0 = z1

           ! Test sortie
           if (test_exit_grid(icell0, x0, y0, z0)) exit inf_loop

           call cross_cell(x0,y0,z0, u,v,w,  icell0, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

           if (icell0 <= n_cells) then
              if (type==1) then
                 factor = CD_units * densite_gaz(icell0) ! column density
              else if (type==2) then
                 factor = kappa(p_icell,lambda) * kappa_factor(icell0)! optical depth, kappa in AU^-1
              else
                 factor = CD_units * densite_gaz(icell0) * tab_abundance(icell0) ! molecular column density
              endif
              sum = sum + l_contrib * factor
           endif
        enddo inf_loop
        column(icell,direction) = sum
     enddo ! icell
     !$omp enddo
     !$omp end parallel
  end do ! direction

  return

end subroutine compute_column

!***********************************************************

subroutine integ_ray_mol(id,imol,icell_in,x,y,z,u,v,w,iray,labs, ispeed,tab_speed, nTrans, tab_Trans)
  ! Generalisation de la routine physical_length
  ! pour le cas du transfert dans les raies
  ! Propage un paquet depuis un point d'origine donne
  ! et integre l'equation du transfert radiatif
  ! La propagation doit etre a l'envers pour faire du
  ! ray tracing  !!
  !
  ! C. Pinte
  ! 12/04/07

  implicit none

  integer, intent(in) :: id, imol,icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs

  integer, dimension(2), intent(in) :: ispeed
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  integer, intent(in) :: nTrans
  integer, dimension(nTrans), intent(in) :: tab_Trans

  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: P, dtau, dtau2, Snu, opacite
  real(kind=dp), dimension(ispeed(1):ispeed(2),nTrans) :: tau, tau2
  real(kind=dp), dimension(nTrans) :: tau_c
  real(kind=dp) :: dtau_c, Snu_c, kappa_cont
  integer :: i, iTrans, nbr_cell, next_cell, previous_cell, icell_star, i_star
  integer, target :: icell

  real :: facteur_tau

  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

  integer, pointer :: p_icell

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau(:,:) = 0.0_dp
  I0(:,:,iray,id) = 0.0_dp

  tau_c(:) = 0.0_dp
  I0c(:,iray,id) = 0.0_dp

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif


  ! Will the ray intersect a star
  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! propagation dans la grille
  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     icell = next_cell
     x0=x1 ; y0=y1 ; z0=z1

     if (icell <= n_cells) then
        lcellule_non_vide=.true.
     else
        lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) then
        return
     endif
     if (lintersect_stars) then
        if (icell == icell_star) return
     endif

     nbr_cell = nbr_cell + 1

     ! Calcul longeur de vol et profondeur optique dans la cellule
     previous_cell = 0 ! unused, just for Voronoi
     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     if (lcellule_non_vide) then
        lsubtract_avg = ((nbr_cell == 1).and.labs)

        ! local line profile mutiplied by frequency
        P(:) = local_line_profile(icell,lsubtract_avg,x0,y0,z0,x1,y1,z1,u,v,w,l_void_before,l_contrib,ispeed,tab_speed)

        if ((nbr_cell == 1).and.labs) then
           ds(iray,id) = l_contrib
           Doppler_P_x_freq(:,iray,id) = P(:)
        endif

        ! surface superieure ou inf
        facteur_tau = 1.0
        if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
        if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

        do i=1,nTrans
           iTrans = tab_Trans(i) ! selecting the proper transition for ray-tracing

           kappa_cont = kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell)

           opacite(:) = kappa_mol_o_freq(icell,iTrans) * P(:) * facteur_tau + kappa_cont

           ! Epaisseur optique
           dtau(:) =  l_contrib * opacite(:)
           dtau_c = l_contrib * kappa_cont

           ! Fonction source
           Snu(:) = ( emissivite_mol_o_freq(icell,iTrans) * P(:) * facteur_tau &
                + emissivite_dust(icell,iTrans) ) / (opacite(:) + 1.0e-300_dp)
           Snu_c = emissivite_dust(icell,iTrans) / (kappa_cont  + 1.0e-300_dp)

           ! Warning I0, I0c (and origine_mol) are smaller arrays (dimension nTrans)
           I0(:,i,iray,id) = I0(:,i,iray,id) + &
                exp(-tau(:,i)) * (1.0_dp - exp(-dtau(:))) * Snu(:)
           I0c(i,iray,id) = I0c(i,iray,id) + &
                exp(-tau_c(i)) * (1.0_dp - exp(-dtau_c)) * Snu_c

           if (lorigine.and.(.not.labs)) then
              origine_mol(:,i,icell,id) = origine_mol(:,i,icell,id) + &
                   exp(-tau(:,i)) * (1.0_dp - exp(-dtau(:))) * Snu(:)
           endif

           ! Mise a jour profondeur optique pour cellule suivante
           ! Warning tau and  tau_c are smaller array (dimension nTrans)
           tau(:,i) = tau(:,i) + dtau(:)
           tau_c(i) = tau_c(i) + dtau_c
        enddo ! i

        if (ldouble_RT) then
           do i=1,nTrans
              iTrans = tab_Trans(i) ! selecting the proper transition for ray-tracing

              opacite(:) = kappa_mol_o_freq2(icell,iTrans) * P(:) + kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell)
              dtau(:) =  l_contrib * opacite(:)

              ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
              ! la profondeur optique jusqu'a la cellule
              Snu(:) = ( emissivite_mol_o_freq2(icell,iTrans) * P(:) + &
                   emissivite_dust(icell,iTrans) ) / (opacite(:) + 1.0e-30_dp)
              I02(:,iTrans,iray,id) = I02(:,iTrans,iray,id) + &
                   exp(-tau2(:,iTrans)) * (1.0_dp - exp(-dtau2(:))) * Snu(:)

              ! Mise a jour profondeur optique pour cellule suivante
              ! Warning tau2 is a smaller array (dimension nTrans)
              tau2(:,i) = tau2(:,i) + dtau2(:)
           enddo
        endif

     endif  ! lcellule_non_vide

  enddo infinie

  ! Ajout cmb, pondere par la profondeur optique totale
  !  tspeed(:) = tab_speed(:) + v_avg0
  !  I0(:,:,iray,id) = I0(:,:,iray,id) + Cmb(ispeed,tspeed) * exp(-tau(:,:))
  !  if (ldouble_RT) then
  !     I02(:,:,iray,id) = I02(:,:,iray,id) + Cmb(ispeed,tspeed) * exp(-tau2(:,:))
  !  endif

  do i=1,nTrans
     iTrans = tab_Trans(i)
     I0(:,i,iray,id) = I0(:,i,iray,id) + tab_Cmb_mol(iTrans) * exp(-tau(:,i))
  enddo

  if (ldouble_RT) then
     do i=1,nTrans
        iTrans = tab_Trans(i)
        I02(:,i,iray,id) = I02(:,i,iray,id) + tab_Cmb_mol(iTrans) * exp(-tau2(:,i))
     enddo
  endif

  return

end subroutine integ_ray_mol

!***********************************************************

subroutine physical_length_mol(imol,iTrans,icell_in,x,y,z,u,v,w, ispeed, tab_speed,tau_threshold,flag_sortie)
  ! Copmputes position where a given optical depth is reached
  ! This is simplified version of integ_ray_mol (also inspired by physical_length)
  ! only computes optical depth and stopes where given tau_max is reached

  integer, intent(in) :: imol, icell_in
  real(kind=dp), intent(inout) :: x,y,z
  real(kind=dp), intent(in) :: u,v,w
  real :: tau_threshold
  integer, dimension(2), intent(in) :: ispeed
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed
  logical, intent(out) :: flag_sortie

  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, ltot
  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: P, tau_mol, dtau_mol, opacite
  real(kind=dp) :: tau_max, tau_previous, facteur_tau

  integer :: i, iTrans, nbr_cell, next_cell, previous_cell, iv
  integer, target :: icell

  logical :: lcellule_non_vide, lstop
  logical, parameter :: lsubtract_avg = .false.

  integer, pointer :: p_icell

  x1=x;y1=y;z1=z
  next_cell = icell_in
  nbr_cell = 0

  ltot = 0.0_dp
  lstop = .false.
  flag_sortie = .false.
  tau_mol(:) = 0.0_dp

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     icell = next_cell
     x0=x1 ; y0=y1 ; z0=z1

     if (icell <= n_cells) then
        lcellule_non_vide=.true.
     else
        lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) then
         flag_sortie = .true.
        return
     endif

     nbr_cell = nbr_cell + 1

     ! Calcul longeur de vol et profondeur optique dans la cellule
     previous_cell = 0 ! unused, just for Voronoi
     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     if (lcellule_non_vide) then
        ! local line profile mutiplied by frequency
        P(:) = local_line_profile(icell,lsubtract_avg,x0,y0,z0,x1,y1,z1,u,v,w,l_void_before,l_contrib,ispeed,tab_speed)

        ! surface superieure ou inf
        facteur_tau = 1.0
        if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
        if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

        !do i=1,nTrans
        !iTrans = tab_Trans(i) ! selecting the proper transition for ray-tracing
        opacite(:) = kappa_mol_o_freq(icell,iTrans) * P(:) * facteur_tau &
             + kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell)

        ! Epaisseur optique
        dtau_mol(:) =  l_contrib * opacite(:)

        ! Mise a jour profondeur optique pour cellule suivante
        ! Warning tau and  tau_c are smaller array (dimension nTrans)
        tau_mol(:) = tau_mol(:) + dtau_mol(:)
        tau_max =  maxval(tau_mol(:))

        if (tau_max > tau_threshold) then
           lstop = .true.
           iv = maxloc(tau_mol(:),dim=1) + ispeed(1)-1

           tau_previous = tau_mol(iv) - dtau_mol(iv)

           ! rescaling l_contrib so that tau_max = tau_threshold
           l_contrib = l_contrib  * (tau_threshold-tau_previous)/(tau_max-tau_previous)
           l = l_void_before + l_contrib
           ltot=ltot+l
        else
           ltot=ltot+l
        endif
     else
        ltot=ltot+l
     endif  ! lcellule_non_vide


     ! On a fini d'integrer : sortie de la routine
     if (lstop) then
        flag_sortie = .false.
        ! we recompute the position
        x=x0+l*u
        y=y0+l*v
        z=z0+l*w

        if (.not.lVoronoi) then
           if (l3D) then
              if (lcylindrical) call indice_cellule(x,y,z, previous_cell)
           endif
        endif ! todo : on ne fait rien dans la cas Voronoi ???

        return
     endif ! lstop

  enddo infinie

  return

end subroutine physical_length_mol

!********************************************************************

subroutine physical_length_mol_Flux(imol,iTrans,icell_in,x,y,z,u,v,w, ispeed, tab_speed,Flux_threshold,flag_sortie)
  ! Copmputes position where a given optical depth is reached
  ! This is simplified version of integ_ray_mol (also inspired by physical_length)
  ! only computes optical depth and stopes where given tau_max is reached

  integer, intent(in) :: imol, icell_in
  real(kind=dp), intent(inout) :: x,y,z
  real(kind=dp), intent(in) :: u,v,w, Flux_threshold
  integer, dimension(2), intent(in) :: ispeed
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed
  logical, intent(out) :: flag_sortie

  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, ltot
  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: P, tau_mol, dtau_mol, opacite, I_mol, Snu, dI_mol
  real(kind=dp) :: I_max, I_previous

  integer :: i, iTrans, nbr_cell, next_cell, previous_cell
  integer, target :: icell

  logical :: lcellule_non_vide, lstop
  logical, parameter :: lsubtract_avg = .false.

  integer, pointer :: p_icell

  x1=x;y1=y;z1=z
  next_cell = icell_in
  nbr_cell = 0

  ltot = 0.0_dp
  lstop = .false.
  flag_sortie = .false.
  tau_mol(:) = 0.0_dp
  I_mol(:) = 0.0_dp

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     icell = next_cell
     x0=x1 ; y0=y1 ; z0=z1

     if (icell <= n_cells) then
        lcellule_non_vide=.true.
     else
        lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) then
         flag_sortie = .true.
        return
     endif

     nbr_cell = nbr_cell + 1

     ! Calcul longeur de vol et profondeur optique dans la cellule
     previous_cell = 0 ! unused, just for Voronoi
     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     if (lcellule_non_vide) then
        ! local line profile mutiplied by frequency
        P(:) = local_line_profile(icell,lsubtract_avg,x0,y0,z0,x1,y1,z1,u,v,w,l_void_before,l_contrib,ispeed,tab_speed)

        !do i=1,nTrans
        !iTrans = tab_Trans(i) ! selecting the proper transition for ray-tracing

        opacite(:) = kappa_mol_o_freq(icell,iTrans) * P(:) + kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell)

        ! Epaisseur optique
        dtau_mol(:) =  l_contrib * opacite(:)

        ! Fonction source
        Snu(:) = ( emissivite_mol_o_freq(icell,iTrans) * P(:) &
             + emissivite_dust(icell,iTrans) ) / (opacite(:) + 1.0e-300_dp)

        ! Specific intensity
        dI_mol(:) = exp(-tau_mol(:)) * (1.0_dp - exp(-dtau_mol(:))) * Snu(:)
        I_mol(:) = I_mol(:) + dI_mol(:)

        ! Mise a jour profondeur optique pour cellule suivante
        ! Warning tau and  tau_c are smaller array (dimension nTrans)
        tau_mol(:) = tau_mol(:) + dtau_mol(:)

        I_max =  maxval(I_mol(:))

        if (I_max > Flux_threshold) then
           lstop = .true.
           I_previous = maxval(I_mol(:)-dI_mol(:))

           ! rescaling l_contrib so that tau_max = tau_threshold
           if (I_max-I_previous > 0) l_contrib = l_contrib  * (Flux_threshold-I_previous)/(I_max-I_previous)
           l = l_void_before + l_contrib
           ltot=ltot+l
        else
           ltot=ltot+l
        endif
     else
        ltot=ltot+l
     endif  ! lcellule_non_vide


     ! On a fini d'integrer : sortie de la routine
     if (lstop) then
        flag_sortie = .false.
        ! we recompute the position
        x=x0+l*u
        y=y0+l*v
        z=z0+l*w

        if (.not.lVoronoi) then
           if (l3D) then
              if (lcylindrical) call indice_cellule(x,y,z, previous_cell)
           endif
        endif ! todo : on ne fait rien dans la cas Voronoi ???

        return
     endif ! lstop

  enddo infinie

  return

end subroutine physical_length_mol_Flux

!********************************************************************

function local_line_profile(icell,lsubtract_avg, x0,y0,z0,x1,y1,z1,u,v,w,l_void_before,l_contrib,ispeed,tab_speed)

  integer, intent(in) :: icell
  logical, intent(in) :: lsubtract_avg
  real(kind=dp), intent(in) :: x0,y0,z0,x1,y1,z1,u,v,w,l_void_before,l_contrib
  integer, dimension(2), intent(in) :: ispeed
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: local_line_profile

  integer, parameter :: n_vpoints_max = 200 ! pas super critique
  ! presque OK avec 2 pour la simu Herbig de Peter (2x plus vite)
  real(kind=dp), dimension(n_vpoints_max) :: vitesse
  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: tspeed
  real(kind=dp) ::  v0, v1, v_avg0, delta_vol_phi, xphi, yphi, zphi
  integer :: ivpoint, n_vpoints

  ! Differentiel de vitesse au travers de la cellule
  !dv = dv_proj(ri0,zj0,x0,y0,z0,x1,y1,z1,u,v,w)
  v0 = v_proj(icell,x0,y0,z0,u,v,w)

  if (lVoronoi) then ! Velocity is constant in cell
     n_vpoints = 1
     vitesse(1) = v0
  else ! Velocity is varying accross cell
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = abs(v1 - v0)

     ! Nbre de points d'integration en fct du differentiel de vitesse
     ! compare a la largeur de raie de la cellule de depart
     n_vpoints  = min(max(2,nint(dv/v_line(icell)*20.)),n_vpoints_max)

     ! Vitesse projete le long du trajet dans la cellule
     do ivpoint=2, n_vpoints-1
        delta_vol_phi = l_void_before + (real(ivpoint,kind=dp))/(real(n_vpoints,kind=dp)) * l_contrib
        xphi=x0+delta_vol_phi*u
        yphi=y0+delta_vol_phi*v
        zphi=z0+delta_vol_phi*w
        vitesse(ivpoint) = v_proj(icell,xphi,yphi,zphi,u,v,w)
     enddo
     vitesse(1) = v0
     vitesse(n_vpoints) = v1
  endif

  if (lsubtract_avg) then
     v_avg0 = 0.0_dp
     do ivpoint=1,n_vpoints
        v_avg0 = v_avg0 + vitesse(ivpoint)
     enddo
     v_avg0 = v_avg0 / real(n_vpoints,kind=dp)
  else
     v_avg0 = 0.0_dp
  endif

  ! Profil de raie local integre a multiplier par la frequence de la transition
  local_line_profile(:) = 0.0_dp
  do ivpoint=1,n_vpoints
     tspeed(:) = tab_speed(:) - (vitesse(ivpoint) - v_avg0)
     local_line_profile(:) = local_line_profile(:) + phiProf(icell,ispeed,tspeed)
  enddo
  local_line_profile(:) = local_line_profile(:)/n_vpoints

  return

end function local_line_profile

!***********************************************************

subroutine integ_tau_mol(imol)

  implicit none

  integer, intent(in) :: imol

  real ::  norme, norme1, vmax
  integer :: i, j, iTrans, n_speed, icell, it

  integer, dimension(2) :: ispeed
  real(kind=dp), dimension(0:0) :: tab_speed, P
  real(kind=dp) :: x0, y0, z0, u0,v0,w0

  real(kind=dp), dimension(0:0,mol(imol)%nTrans_raytracing) :: tau_mol
  real(kind=dp), dimension(mol(imol)%nTrans_raytracing) :: tau_dust

  ! Here, we are only looking at the line center (for the observer)
  ispeed(1) = 0 ; ispeed(2)  = 0 ; tab_speed(0) = 0.0_dp


  x0=0.0 ; y0=0.0 ; z0=0.0
  u0=1.0 ; v0=0.0 ; w0=0.0

  call indice_cellule(x0,y0,z0, icell)

  call optical_length_tot_mol(imol,icell,x0,y0,z0,u0,v0,w0, ispeed, tab_speed, &
       mol(imol)%nTrans_raytracing ,mol(imol)%indice_Trans_raytracing,tau_mol,tau_dust)

  do it=1, mol(imol)%nTrans_rayTracing
     iTrans = mol(imol)%indice_Trans_rayTracing(it)

     write(*,*) "-------------------------------"
     write(*,*) "Transition J=", j_qnb(itransUpper(iTrans)), "-", j_qnb(itransLower(iTrans))
     write(*,*) "tau_mol = ", real(tau_mol(0,it))
     write(*,*) "tau_dust=", real(tau_dust(it))

     ! Compute altitude at which tau=1 is reached
     if (.not.lVoronoi) then
        loop_r : do i=1,n_rad
           icell = cell_map(i,1,1)
           if (r_grid(icell) > 100.0) then
              norme=0.0
              loop_z : do j=nz, 1, -1
                 icell = cell_map(i,j,1)
                 P(:) = phiProf(icell,ispeed,tab_speed)
                 norme=norme+kappa_mol_o_freq(icell,iTrans)*(z_lim(i,j+1)-z_lim(i,j))*P(0)
                 if (norme > 1.0) then
                    write(*,*) "Vertical Tau_mol=1 (for r=100 au) at z=", real(z_grid(icell)), "au"
                    exit loop_z
                 endif
              enddo loop_z
              if (norme < 1.0) write(*,*) "Vertical Tau_mol=1 (for r=100 au) not reached, tau_max=", norme
              exit loop_r
           endif
        enddo loop_r
     endif

  enddo

  return

end subroutine integ_tau_mol

!********************************************************************

subroutine optical_length_tot_mol(imol,icell_in,x,y,z,u,v,w, ispeed, tab_speed, nTrans, tab_Trans,tau_mol,tau_c)
  ! Optical depth only : simplified version of integ_ray_mol

  integer, intent(in) :: imol, nTrans, icell_in
  real(kind=dp), intent(in) :: x,y,z,u,v,w
  real(kind=dp), dimension(:,:), intent(out) :: tau_mol
  real(kind=dp), dimension(nTrans), intent(out) :: tau_c

  integer, dimension(2), intent(in) :: ispeed
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed
  integer, dimension(nTrans), intent(in) :: tab_Trans

  !integer, dimension(nTrans), intent(in) :: tab_Trans

  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: P, dtau_mol, opacite

  real(kind=dp) :: dtau_c
  integer :: i, iTrans, nbr_cell, next_cell, previous_cell
  integer, target :: icell

  logical :: lcellule_non_vide

  logical, parameter :: lsubtract_avg = .false.

  integer, pointer :: p_icell

  x1=x;y1=y;z1=z
  next_cell = icell_in
  nbr_cell = 0

  tau_mol(:,:) = 0.0_dp
  tau_c(:) = 0.0_dp

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif

  !*** propagation dans la grille

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     icell = next_cell
     x0=x1 ; y0=y1 ; z0=z1

     if (icell <= n_cells) then
        lcellule_non_vide=.true.
     else
        lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) then
        return
     endif

     nbr_cell = nbr_cell + 1

     ! Calcul longeur de vol et profondeur optique dans la cellule
     previous_cell = 0 ! unused, just for Voronoi
     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     if (lcellule_non_vide) then
        ! local line profile mutiplied by frequency
        P(:) = local_line_profile(icell,lsubtract_avg,x0,y0,z0,x1,y1,z1,u,v,w,l_void_before,l_contrib,ispeed,tab_speed)

        do i=1,nTrans
           iTrans = tab_Trans(i) ! selecting the proper transition for ray-tracing

           opacite(:) = kappa_mol_o_freq(icell,iTrans) * P(:) + kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell)

           ! Epaisseur optique
           dtau_mol(:) =  l_contrib * opacite(:)
           dtau_c = l_contrib * kappa_abs_LTE(p_icell,iTrans) * kappa_factor(icell)

           ! Mise a jour profondeur optique pour cellule suivante
           ! Warning tau and  tau_c are smaller array (dimension nTrans)
           tau_mol(:,i) = tau_mol(:,i) + dtau_mol(:)
           tau_c(i) = tau_c(i) + dtau_c
        enddo ! i
     endif  ! lcellule_non_vide

  enddo infinie

  return

end subroutine optical_length_tot_mol

!********************************************************************
   subroutine integ_ray_atom(id,icell_in,x,y,z,u,v,w,iray,labs,N,lambda)
   ! ------------------------------------------------------------------------------- !
   ! TO DO: merge integ_ray_atom + integ_ray_line
   ! Zeeman
   ! scattering
   ! level dissolution
   ! dust
   ! ------------------------------------------------------------------------------- !
      integer, intent(in) :: id, icell_in, iray
      real(kind=dp), intent(in) :: u,v,w
      real(kind=dp), intent(in) :: x,y,z
      logical, intent(in) :: labs
      integer, intent(in) :: N
      real(kind=dp), dimension(N), intent(in) :: lambda
      real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, Q, P(4)
      real(kind=dp), dimension(N) :: Snu, tau, dtau, chi, coronal_irrad
      integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la, icell_prev
      logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

      x1=x;y1=y;z1=z
      x0=x;y0=y;z0=z
      next_cell = icell_in
      nbr_cell = 0
      icell_prev = icell_in

      tau(:) = 0.0_dp

      Itot(:,iray,id) = 0.0_dp

      ! Will the ray intersect a star
      call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)
      ! Boucle infinie sur les cellules (we go over the grid.)
      infinie : do ! Boucle infinie
      ! Indice de la cellule
         icell = next_cell
         x0=x1 ; y0=y1 ; z0=z1

         lcellule_non_vide = (icell <= n_cells)
         ! if (icell <= n_cells) then
         !    lcellule_non_vide=.true.
         ! else
         !    lcellule_non_vide=.false.
         ! endif

         ! Test sortie ! "The ray has reach the end of the grid"
         if (test_exit_grid(icell, x0, y0, z0)) return

         if (lintersect_stars) then !"will interesct"
            if (icell == icell_star) then!"has intersected"
               Itot(:,iray,id) = Itot(:,iray,id) + exp(-tau(:)) * &
                  star_rad(id,iray,i_star,icell_prev,x0,y0,z0,u,v,w,N,lambda)
               return
            end if
         endif
         !With the Voronoi grid, somme cells can have a negative index
         !therefore we need to test_exit_grid before using icompute_atom_rt
         if (icell <= n_cells) then
            lcellule_non_vide = (icompute_atomRT(icell) > 0)
            if (icompute_atomRT(icell) < 0) then
               if (icompute_atomRT(icell) == -1) then
                  !If the optically thick region (dark zone) has a temperature
                  !add a black body emission and leave.
                  if (T(icell) > 0.0_dp) Itot(:,iray,id) = Itot(:,iray,id) + &
                                    exp(-tau) * Bpnu(N,lambda,T(icell))
                  return
               else
                  !Does not return but cell is empty (lcellule_non_vide is .false.)
                  coronal_irrad = linear_1D_sorted(atmos_1d%Ncorona,atmos_1d%x_coro,atmos_1d%I_coro(:,1),N,lambda)
                  Itot(:,iray,id) = Itot(:,iray,id) + exp(-tau) * coronal_irrad
               endif
            endif
         endif

         nbr_cell = nbr_cell + 1

         ! Calcul longeur de vol et profondeur optique dans la cellule
         previous_cell = 0 ! unused, just for Voronoi
         call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)

         !count opacity only if the cell is filled, else go to next cell
         if (lcellule_non_vide) then
            lsubtract_avg = ((nbr_cell == 1).and.labs)
            ! opacities in m^-1, l_contrib in au


            call contopac_atom_loc(icell, N, lambda, chi, Snu)
            call opacity_atom_bb_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,&
               l_void_before,l_contrib,lsubtract_avg,N,lambda,chi,Snu)

            dtau(:) = l_contrib * chi(:) * AU_to_m !au * m^-1 * au_to_m

            if (lsubtract_avg) then
               !Lambda operator / chi_dag
               !force PSI to be ray-by-ray but not ds !
               !local, unaffected by vel.
               psi(:,1,id) = ( 1.0_dp - exp( -dtau(:) ) ) / chi
               ds(iray,id) = l_contrib * AU_to_m
            endif

            ! if (lorigine) then
            !    if (maxval(ori(:,icell,id))==0.0_dp) then
            !       ori(:,icell,id) = ori(:,icell,id) + eta(:,id) * exp(-tau(:))
            !       tet(:,icell,id) = tet(:,icell,id) + tau(:) * exp(-tau(:))
            !    endif
            ! endif

            Snu = Snu / chi

            Itot(:,iray,id) = Itot(:,iray,id) + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu
            tau(:) = tau(:) + dtau(:) !for next cell

         end if  ! lcellule_non_vide

         icell_prev = icell
         !duplicate with previous_cell, but this avoid problem with Voronoi grid here

      end do infinie

      return
   end subroutine integ_ray_atom

!********************************************************************

function integ_ray_dust(lambda,icell_in,x,y,z,u,v,w)
  ! Generalisation de la routine physical_length
  ! Propage un paquet depuis un point d'origine donne
  ! et integre l'equation du transfert radiatif
  ! La propagation doit etre a l'envers pour faire du
  ! ray tracing  !!
  !
  ! C. Pinte
  ! 23/01/08

  ! TODO : faire peter le phi ??? Ne sert que pour les champs de vitesse

  implicit none

  integer, intent(in) :: lambda, icell_in
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z

  real(kind=dp), dimension(N_type_flux) :: integ_ray_dust

  real(kind=dp) :: x0, y0, z0, x1, y1, z1, xm, ym, zm, l, l_contrib, l_void_before
  integer :: previous_cell, next_cell, icell_star, i_star
  integer, target :: icell

  real(kind=dp) :: tau, dtau

  logical :: lcellule_non_vide, lintersect_stars

  integer, pointer :: p_icell

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in

  tau = 0.0_dp
  integ_ray_dust(:) = 0.0_dp

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif

  ! Will the ray intersect a star
  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! propagation dans la grille
  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     icell=next_cell
     x0=x1 ; y0=y1 ; z0=z1

     if (icell <= n_cells) then
        lcellule_non_vide=.true.
     else
        lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) return
     if (lintersect_stars) then
        if (icell == icell_star) return
     endif

     ! Calcul longeur de vol et profondeur optique dans la cellule
     previous_cell = 0 ! unused, just for Voronoi
     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

     if (lcellule_non_vide) then
        ! Epaisseur optique de la cellule
        dtau =  l_contrib * kappa(p_icell,lambda) * kappa_factor(icell)

        ! Fct source au milieu du parcours dans la cellule
        xm = 0.5 * (x0 + x1)
        ym = 0.5 * (y0 + y1)
        zm = 0.5 * (z0 + z1)

        ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
        ! la profondeur optique jusqu'a la cellule
        integ_ray_dust(:) = integ_ray_dust(:) + &
             exp(-tau) * (1.0_dp - exp(-dtau)) * dust_source_fct(icell, xm,ym,zm)

        ! Mise a jour profondeur optique pour cellule suivante
        tau = tau + dtau

        ! Pas besoin d'integrer trop profond
        if (tau > tau_dark_zone_obs) return
     endif  ! lcellule_non_vide

  enddo infinie

  return

end function integ_ray_dust

!***********************************************************

subroutine define_dark_zone(lambda,p_lambda,tau_max,ldiff_approx)
! Definition l'etendue de la zone noire
! definie le tableau logique l_dark_zone
! et les rayons limites r_in_opacite pour le premier rayon
! C. Pinte
! 22/04/05

  implicit none

  integer, parameter :: nbre_angle = 11

  integer, intent(in) :: lambda, p_lambda
  real, intent(in) :: tau_max
  logical, intent(in) :: ldiff_approx
  integer :: i, j, pk, n, id, jj
  integer, target :: icell
  real(kind=dp) :: x0, y0, z0, u0, v0, w0
  real :: somme, angle, dvol1, phi, r0

  logical :: flag_direct_star = .false.
  logical :: flag_star = .false.
  logical :: flag_sortie, lpacket_alive

  real(kind=dp), dimension(4) :: Stokes

  integer, pointer :: p_icell

  lpacket_alive = .true.

  if (lvariable_dust) then
     p_icell => icell
  else
     p_icell => icell_ref
  endif

  do pk=1, n_az
     ri_in_dark_zone(pk)=n_rad
     ri_out_dark_zone(pk)=1
     ! étape 1 : radialement depuis le centre
     somme = 0.0
     do1 : do i=1,n_rad
        icell = cell_map(i,1,pk)
        somme=somme+kappa(p_icell,lambda)*kappa_factor(icell)*(r_lim(i)-r_lim(i-1))
        if (somme > tau_max) then
           ri_in_dark_zone(pk) = i
           exit do1
        endif
     enddo do1

     ! étape 2 : radialement depuis rout
     somme = 0.0
     do2 : do i=n_rad,1,-1
        icell = cell_map(i,1,pk)
        somme=somme+kappa(p_icell,lambda)*kappa_factor(icell)*(r_lim(i)-r_lim(i-1))
        if (somme > tau_max) then
           ri_out_dark_zone(pk) = i
           exit do2
        endif
     enddo do2
     if (ri_out_dark_zone(pk)==n_rad) ri_out_dark_zone(pk)=n_rad-1

     if (lcylindrical) then
        ! étape 3 : verticalement
        do i=ri_in_dark_zone(pk), ri_out_dark_zone(pk)
           somme = 0.0
           do3 : do j=nz, 1, -1
              icell = cell_map(i,j,pk)
              somme=somme+kappa(p_icell,lambda)*kappa_factor(icell)*(z_lim(i,j+1)-z_lim(i,j))
              if (somme > tau_max) then
                 zj_sup_dark_zone(i,pk) = j
                 exit do3
              endif
           enddo do3
        enddo

        ! étape 3.5 : verticalement dans autre sens
        if (l3D) then
           do i=ri_in_dark_zone(pk), ri_out_dark_zone(pk)
              somme = 0.0
              do3_5 : do j=-nz, -1
                 icell = cell_map(i,j,pk)
                 somme=somme+kappa(p_icell,lambda)*kappa_factor(icell)*(z_lim(i,abs(j)+1)-z_lim(i,abs(j)))
                 if (somme > tau_max) then
                    zj_inf_dark_zone(i,pk) = j
                    exit do3_5
                 endif
              enddo do3_5
           enddo
        endif
     else ! spherical
        zj_sup_dark_zone(:,pk) = nz
     endif

  enddo !pk


  l_is_dark_zone = .false.
  l_dark_zone(:) = .false.

  ! étape 4 : test sur tous les angles
  if (.not.l3D) then
     cell : do i=max(ri_in_dark_zone(1),2), ri_out_dark_zone(1)
        do j=zj_sup_dark_zone(i,1),1,-1
           icell = cell_map(i,j,1)
           do n=1,nbre_angle
              id=1
              ! position et direction vol
              angle= pi * real(n)/real(nbre_angle+1)! entre 0 et pi
              x0=r_grid(icell) !x0=1.00001*r_lim(i-1) ! cellule 1 traitee a part
              y0=0.0
              z0=z_grid(icell) !z0=0.99999*z_lim(i,j+1)
              u0=cos(angle)
              v0=0.0
              w0=sin(angle)
              Stokes(:) = 0.0_dp ; !Stokes(1) = 1.0_dp ; ! Pourquoi c'etait a 1 ?? ca fausse les chmps de radiation !!!
              call physical_length(id,lambda,p_lambda,Stokes,icell, x0,y0,z0,u0,v0,w0, &
                   flag_star,flag_direct_star,tau_max,dvol1,flag_sortie,lpacket_alive)
              if (.not.flag_sortie) then ! le photon ne sort pas
                 ! la cellule et celles en dessous sont dans la zone noire
                 do jj=1,j
                    icell = cell_map(i,jj,1)
                    l_dark_zone(icell) = .true.
                 enddo
                 l_is_dark_zone = .true.
                 ! on passe a la cellule suivante
                 cycle cell
              endif
           enddo
        enddo
     enddo cell
  else !3D
     do pk=1, n_az
        phi = 2*pi * (real(pk)-0.5)/real(n_az)
        cell_3D : do i=max(ri_in_dark_zone(pk),2), ri_out_dark_zone(pk)
           do j=zj_sup_dark_zone(i,pk),1,-1
              icell = cell_map(i,j,pk)
              do n=1,nbre_angle
                 id=1
                 ! position et direction vol
                 angle= pi * real(n)/real(nbre_angle+1)! entre 0 et pi
                 r0=r_grid(icell)!1.00001*r_lim(i-1) ! cellule 1 traitee a part
                 x0 = r0 *cos(phi)
                 y0 = r0 * sin(phi)
                 z0=z_grid(icell)!z0.99999*z_lim(i,j+1)
                 u0=cos(angle)
                 v0=0.0
                 w0=sin(angle)
                 Stokes(:) = 0.0_dp ; Stokes(1) = 1.0_dp ;
                 call physical_length(id,lambda,p_lambda,Stokes,icell,x0,y0,z0,u0,v0,w0, &
                      flag_star,flag_direct_star,tau_max,dvol1,flag_sortie,lpacket_alive)
                 if (.not.flag_sortie) then ! le photon ne sort pas
                    ! la cellule et celles en dessous sont dans la zone noire
                    do jj=1,j
                       icell = cell_map(i,jj,pk)
                       l_dark_zone(icell) = .true.
                    enddo
                    ! on passe a la cellule suivante
                    cycle cell_3D
                 endif
              enddo
           enddo
        enddo cell_3D

        cell_3D_2 : do i=max(ri_in_dark_zone(pk),2), ri_out_dark_zone(pk)
           do j=zj_inf_dark_zone(i,pk),-1
              icell = cell_map(i,j,pk)
              do n=1,nbre_angle
                 id=1
                 ! position et direction vol
                 angle= pi * real(n)/real(nbre_angle+1)! entre 0 et pi
                 r0=r_grid(icell)!1.00001*r_lim(i-1) ! cellule 1 traitee a part
                 x0 = r0 *cos(phi)
                 y0 = r0 * sin(phi)
                 z0=-z_grid(icell)!-0.99999*z_lim(i,abs(j)+1)
                 u0=cos(angle)
                 v0=0.0
                 w0=sin(angle)
                 Stokes(:) = 0.0_dp ; Stokes(1) = 1.0_dp ;
                 call physical_length(id,lambda,p_lambda,Stokes,icell,x0,y0,z0,u0,v0,w0, &
                      flag_star,flag_direct_star,tau_max,dvol1,flag_sortie,lpacket_alive)
                 if (.not.flag_sortie) then ! le photon ne sort pas
                    ! la cellule et celles en dessous sont dans la zone noire
                    do jj=1,-1
                       icell = cell_map(i,j,pk)
                       l_dark_zone(icell) = .true.
                    enddo
                    l_is_dark_zone=.true.
                    ! on passe a la cellule suivante
                    cycle cell_3D_2
                 endif
              enddo
           enddo
        enddo cell_3D_2
     enddo !pk
  endif !l3D

  do pk=1, n_az
     zj_sup_dark_zone(1:ri_in_dark_zone(pk)-1,pk) = zj_sup_dark_zone(ri_in_dark_zone(pk),pk)
     zj_sup_dark_zone(ri_out_dark_zone(pk)+1:n_rad,pk) = zj_sup_dark_zone(ri_out_dark_zone(pk),pk)
     if (l3D) then
        zj_inf_dark_zone(1:ri_in_dark_zone(pk)-1,pk) = zj_inf_dark_zone(ri_in_dark_zone(pk),pk)
        zj_inf_dark_zone(ri_out_dark_zone(pk)+1:n_rad,pk) = zj_inf_dark_zone(ri_out_dark_zone(pk),pk)
     endif
  enddo

  if ((ldiff_approx).and.(n_rad > 1)) then
     if (minval(ri_in_dark_zone(:))==1) call error("first cell is in diffusion approximation zone", &
          msg2="Increase spatial grid resolution")
  endif


  if (n_zones > 1) then
     do icell=1, n_cells
        if (sum(densite_pouss(:,icell)) < tiny_real) l_dark_zone(icell) = .false.
     enddo
  endif

  do i=1, n_regions
     do j=1,nz
        l_dark_zone(cell_map(regions(i)%iRmin,j,1)) = .false.
        l_dark_zone(cell_map(regions(i)%iRmax,j,1)) = .false.
     enddo
  enddo

  return

end subroutine define_dark_zone

!***********************************************************

subroutine no_dark_zone()
! Definie les variables quand il n'y a pas de zone noire
! C . Pinte
! 22/04/05

  implicit none

  l_dark_zone(:)=.false.

  return

end subroutine no_dark_zone

!***********************************************************

end module optical_depth
