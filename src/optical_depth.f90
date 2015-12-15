module optical_depth

  use wall
  use parametres
  use disk
  use opacity
  use constantes
  use em_th
  use molecular_emission
  use ray_tracing
  use grains, only : tab_lambda
  use utils
  use molecules

  use dust_ray_tracing
  use grid

  implicit none

  contains

subroutine length_deg2(id,lambda,Stokes,ri,zj,phik,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
! Routine chapeau aux sous-routines suivant la geometrie
! C. Pinte
! 03/06/07

  implicit none

  integer, intent(in) :: id,lambda
  integer, intent(inout) :: ri,zj,phik
  real(kind=db), dimension(4), intent(in) :: Stokes
  logical, intent(in) :: flag_star,flag_direct_star
  real(kind=db), intent(inout) :: u,v,w
  real, intent(in) :: extrin
  real(kind=db), intent(inout) :: xio,yio,zio
  real, intent(out) :: ltot
  logical, intent(out) :: flag_sortie

  !call test_convert()

  if (lcylindrical) then
     if (l3D) then
        call length_deg2_3D(id,lambda,Stokes,ri,zj,phik,xio,yio,zio,u,v,w,flag_star,extrin,ltot,flag_sortie)
     else
        if (lopacity_wall) then
           call length_deg2_opacity_wall(id,lambda,Stokes,ri,zj,xio,yio,zio,u,v,w,extrin,ltot,flag_sortie)
        else
           call length_deg2_cyl(id,lambda,Stokes,ri,zj,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
        endif
     endif
  else ! spherical
     call length_deg2_sph(id,lambda,Stokes,ri,zj,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
  endif

  return

end subroutine length_deg2

!********************************************************************

subroutine length_deg2_cyl(id,lambda,Stokes,ri,zj,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
! Integration par calcul de la position de l'interface entre cellules
! par eq deg2 en r et deg 1 en z
! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
! C. Pinte
! 05/02/05

  implicit none

  integer, intent(in) :: id,lambda
  integer, intent(inout) :: ri,zj
  real(kind=db), dimension(4), intent(in) :: Stokes
  logical, intent(in) :: flag_star, flag_direct_star
  real(kind=db), intent(inout) :: u,v,w
  real, intent(in) :: extrin
  real(kind=db), intent(inout) :: xio,yio,zio
  real, intent(out) :: ltot
  logical, intent(out) :: flag_sortie

  real(kind=db) :: x0, y0, z0, x1, y1, z1, x_old, y_old, z_old, extr, inv_w
  real(kind=db) :: l, tau
  real(kind=db) :: correct_moins, correct_plus
  real(kind=db) :: phi_vol, factor
  integer :: ri0, zj0, ri1, zj1, ri_old, zj_old
  integer :: ri_in, zj_in, tmp_k

  integer :: icell0, next_cell, previous_cell, icell

  logical :: lcellule_non_vide, lstop

!  flag_direct_star = .false. ! Ok c'est moins bon edge-on :-)

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille
  lstop = .false.

  x0=xio;y0=yio;z0=zio
  x1=xio;y1=yio;z1=zio
  extr=extrin
  ri_in = ri
  zj_in = zj
  ri0=ri
  zj0=zj
  ri1=ri
  zj1=zj

  ltot=0.0

  if (abs(w) > tiny_real) then
     inv_w=1.0_db/w
  else
     inv_w=sign(huge_db,w) ! huge_real avant
  endif

  ! pour ray-tracing
  phi_vol = atan2(v,u) + deux_pi ! deux_pi pour assurer diff avec phi_pos > 0

  ! Calcule les angles de diffusion pour la direction de propagation donnee
  if ((.not.letape_th).and.lscatt_ray_tracing1) call angles_scatt_rt1(id,u,v,w)

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri_old = ri0 ; zj_old = zj0
     x_old = x0 ; y_old = y0 ; z_old = z0

     ri0=ri1 ; zj0=zj1
     x0=x1 ; y0=y1 ; z0=z1

     icell0 = cell_map(ri0,zj0,1)

     ! Pour cas avec approximation de diffusion
     if (icell0 <= n_cells) then
        if (l_dark_zone(icell0)) then
           ! On revoie le paquet dans l'autre sens
           u = -u ; v = -v ; w=-w
           inv_w = -inv_w
           ! et on le renvoie au point de depart
           ri = ri_old ; zj = zj_old
           xio = x_old ; yio = y_old ; zio = z_old
           flag_sortie= .false.
           return
        endif
     endif

     lcellule_non_vide=.true.
     if (ri0==0) lcellule_non_vide=.false.

     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        flag_sortie = .true.
        return
     elseif (zj0>nz) then
        lcellule_non_vide=.false.

        ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           flag_sortie = .true.
           return
        endif
     endif ! Test sortie


     !call cylindrical2cell(ri0,zj0,1, cell) ! tmp : this routine should only know cell in the long term, not ri0, etc
     previous_cell = 0 ! unusued, just for Voronoi
     call cross_cylindrical_cell(lambda, x0,y0,z0, u,v,w,  icell0, previous_cell, x1,y1,z1, next_cell, l, tau)
     call cell2cylindrical(next_cell, ri1,zj1,tmp_k) ! tmp : this routine should only know cell in the long term


     !call cross_cylindrical_cell_tmp(lambda, x0,y0,z0, u,v,w, ri0, zj0,  x1,y1,z1, ri1, zj1, l, tau)

     ! Comparaison integrale avec tau
     ! et ajustement longueur de vol evntuellement
     if(tau > extr) then ! On a fini d'integrer
        lstop = .true.
        l = l * (extr/tau) ! on rescale l pour que tau=extr
        ltot=ltot+l
     else ! Il reste extr - tau a integrer dans la cellule suivante
        extr=extr-tau
        ltot=ltot+l
     endif

     ! Stockage des champs de radiation
     if (lcellule_non_vide) call save_radiation_field(id,lambda,icell0, Stokes, l,  &
          x0,y0,z0, x1,y1,z1, u,v,w, flag_star, flag_direct_star)

     ! On a fini d'integrer : sortie de la routine
     if (lstop) then
        flag_sortie = .false.
        xio=x0+l*u
        yio=y0+l*v
        zio=z0+l*w
        call indice_cellule(xio,yio,zio,ri,zj)

        ! Patch pour eviter BUG sur position radiale
        ! a cause de limite de precision
        if (ri==0) then
           factor = rmin/ sqrt(xio*xio+yio*yio) * correct_plus
           xio = xio * factor
           yio = yio * factor
           zio = zio * factor

           ! On verifie que c'est OK maintenant
           call indice_cellule(xio,yio,zio,ri,zj)
           if (ri==0) then
              write(*,*) "BUG integ_deg2_cyl"
              write(*,*) "Exiting"
              stop
           endif
        endif

        icell = cell_map(ri,zj,1)
        if (l_dark_zone(icell)) then ! Petit test de securite
           ! On resort le paquet
           if (zj < zj0) then
              zj = zj0
              zio = z_lim(ri0,zj0)*correct_plus
           endif
           if (ri < ri0) then
              ri = ri0
              xio = xio * correct_plus
              yio = yio * correct_plus
           else if (ri > ri0) then
              ri = ri0
              xio = xio * correct_moins
              yio = yio * correct_moins
           endif
        endif
        return
     endif ! lstop

  enddo ! boucle infinie
  write(*,*) "BUG"
  return

end subroutine length_deg2_cyl

!********************************************************************

subroutine cross_cylindrical_cell(lambda, x0,y0,z0, u,v,w,  cell, previous_cell, x1,y1,z1, next_cell, l, tau)
!subroutine cross_cylindrical_cell_tmp(lambda, x0,y0,z0, u,v,w, ri0, zj0,  x1,y1,z1, ri1, zj1, l, tau)

  integer, intent(in) :: lambda, cell, previous_cell
  real(kind=db), intent(in) :: x0,y0,z0
  real(kind=db), intent(in) :: u,v,w ! Todo : check that


  real(kind=db), intent(out) :: x1, y1, z1
  integer, intent(out) :: next_cell
  real(kind=db), intent(out) :: l, tau

  ! Variables to be sorted out
  integer :: ri0,zj0,k0

  real(kind=db) :: inv_a, a, b, c, s, rac, t, delta, inv_w, r_2
  real(kind=db) :: delta_vol, zlim, dotprod, opacite
  real(kind=db) :: correct_moins, correct_plus

  integer ::  delta_rad, delta_zj, ri1, zj1


  ! TODO: Can be calculated outside
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  a=u*u+v*v

  if (a > tiny_real) then
     inv_a=1.0_db/a
  else
     inv_a=huge_real
  endif

  if (abs(w) > tiny_real) then
     inv_w=1.0_db/w
  else
     inv_w=sign(huge_db,w) ! huge_real avant
  endif
  ! End : TODO : Can be calculated outside

  ! 3D cell indices
  call cell2cylindrical(cell, ri0,zj0,k0)

  ! Detection interface
  r_2=x0*x0+y0*y0
  b=(x0*u+y0*v)*inv_a

  if (ri0==0) then
     opacite=0.0_db
     ! Si on est avant le bord interne,  on passe forcement par rmin
     ! et on cherche forcement la racine positive (unique)
     c=(r_2-r_lim_2(0))*inv_a
     delta=b*b-c
     rac=sqrt(delta)
     s = (-b+rac) * correct_plus
     t=huge_real
     delta_rad=1
  else
     if (cell > n_cells) then
        opacite = 0.0_db
     else
        opacite=kappa(cell,lambda)
     endif
     ! 1) position interface radiale
     ! on avance ou recule en r ? -> produit scalaire
     dotprod=u*x0+v*y0  ! ~ b
     if (dotprod < 0.0_db) then
        ! on recule : on cherche rayon inférieur
        c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
        delta=b*b-c
        if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
           ! on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        else
           delta_rad=-1
        endif
     else
        ! on avance : on cherche le rayon supérieur
        c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
        delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
        delta_rad=1
     endif !dotprod
     rac=sqrt(delta)
     s=(-b-rac) * correct_plus
     if (s < 0.0_db) then
        s=(-b+rac) * correct_plus
     else if (s==0.0_db) then
        s=prec_grille
     endif


     ! 2) position interface verticale
     ! on monte ou on descend par rapport au plan équatorial ?
     dotprod=w*z0
     if (dotprod == 0.0_db) then
        t=1.0e10
     else
        if (dotprod > 0.0_db) then
           ! on monte
           if (zj0==nz+1) then
              delta_zj=0
              if (z0 > 0.0_db) then
                 zlim=1.0e10
              else
                 zlim=-1.0e10
              endif
           else
              if (z0 > 0.0) then
                 zlim=z_lim(ri0,zj0+1)*correct_plus
              else
                 zlim=-z_lim(ri0,zj0+1)*correct_plus
              endif
              delta_zj=1
           endif
        else
           ! on descend
           if (zj0==1) then
              ! on traverse le plan eq donc on va remonter
              ! et z va changer de signe
              delta_zj=1
              if (z0 > 0.0_db) then
                 zlim=-z_lim(ri0,2)*correct_moins
              else
                 zlim=z_lim(ri0,2)*correct_moins
              endif
           else !(zj0==1)
              ! on ne traverse pas z=0.
              if (z0 > 0.0_db) then
                 zlim=z_lim(ri0,zj0)*correct_moins
              else
                 zlim=-z_lim(ri0,zj0)*correct_moins
              endif
              delta_zj=-1
           endif !(zj0==1)
        endif ! monte ou descend
        t=(zlim-z0)*inv_w
        ! correct pb precision
        if (t < 0.0_db) t=prec_grille
     endif !dotprod=0.0
  endif ! ri0==0


  ! 3) interface en r ou z ?
  if (s < t) then ! r
     l=s
     delta_vol=s
     ! Position au bord de la cellule suivante
     x1=x0+delta_vol*u
     y1=y0+delta_vol*v
     z1=z0+delta_vol*w
     ri1=ri0+delta_rad
     if ((ri1<1).or.(ri1>n_rad)) then
        zj1=zj0
     else
        zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
        if (zj1>nz) zj1=nz+1
     endif
  else ! z
     l=t
     delta_vol=t
     ! Position au bord de la cellule suivante
     x1=x0+delta_vol*u
     y1=y0+delta_vol*v
     z1=z0+delta_vol*w
     ri1=ri0
     zj1=zj0+delta_zj
  endif

  ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
  if (z1 == 0.0_db) z1 = prec_grille

  ! Calcul longeur de vol et profondeur optique dans la cellule
  tau = l*opacite ! opacite constante dans la cellule

  !call cylindrical2cell(ri1,zj1,1, next_cell)
  next_cell =  cell_map(ri1,zj1,1)

  return

end subroutine cross_cylindrical_cell

!********************************************************************

subroutine cross_cylindrical_cell_tmp(lambda, x0,y0,z0, u,v,w, ri0, zj0,  x1,y1,z1, ri1, zj1, l, tau)

  integer, intent(in) :: lambda, ri0, zj0
  real(kind=db), intent(in) :: x0,y0,z0
  real(kind=db), intent(in) :: u,v,w ! Todo : check that
  real(kind=db), intent(out) :: x1,y1,z1 ! Todo : check that

  integer, intent(out) :: ri1, zj1
  real(kind=db), intent(out) :: l, tau

  real(kind=db) :: inv_a, a, b, c, s, rac, t, delta, inv_w, r_2
  real(kind=db) :: delta_vol, zlim, dotprod, opacite
  real(kind=db) :: correct_moins, correct_plus

  integer ::  delta_rad, delta_zj, icell0

  ! TODO: Can be calculated outside
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  a=u*u+v*v

  if (a > tiny_real) then
     inv_a=1.0_db/a
  else
     inv_a=huge_real
  endif

  if (abs(w) > tiny_real) then
     inv_w=1.0_db/w
  else
     inv_w=sign(huge_db,w) ! huge_real avant
  endif
  ! End : TODO : Can be calculated outside


  ! Detection interface
  r_2=x0*x0+y0*y0
  b=(x0*u+y0*v)*inv_a

  if (ri0==0) then
     opacite=0.0_db
     ! Si on est avant le bord interne,  on passe forcement par rmin
     ! et on cherche forcement la racine positive (unique)
     c=(r_2-r_lim_2(0))*inv_a
     delta=b*b-c
     rac=sqrt(delta)
     s = (-b+rac) * correct_plus
     t=huge_real
     delta_rad=1
  else
     icell0 = cell_map(ri0,zj0,1)
     if (icell0 > n_cells) then
        opacite = 0.0_db
     else
        opacite=kappa(icell0,lambda)
     endif
     ! 1) position interface radiale
     ! on avance ou recule en r ? -> produit scalaire
     dotprod=u*x0+v*y0  ! ~ b
     if (dotprod < 0.0_db) then
        ! on recule : on cherche rayon inférieur
        c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
        delta=b*b-c
        if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
           ! on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        else
           delta_rad=-1
        endif
     else
        ! on avance : on cherche le rayon supérieur
        c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
        delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
        delta_rad=1
     endif !dotprod
     rac=sqrt(delta)
     s=(-b-rac) * correct_plus
     if (s < 0.0_db) then
        s=(-b+rac) * correct_plus
     else if (s==0.0_db) then
        s=prec_grille
     endif


     ! 2) position interface verticale
     ! on monte ou on descend par rapport au plan équatorial ?
     dotprod=w*z0
     if (dotprod == 0.0_db) then
        t=1.0e10
     else
        if (dotprod > 0.0_db) then
           ! on monte
           if (zj0==nz+1) then
              delta_zj=0
              if (z0 > 0.0_db) then
                 zlim=1.0e10
              else
                 zlim=-1.0e10
              endif
           else
              if (z0 > 0.0) then
                 zlim=z_lim(ri0,zj0+1)*correct_plus
              else
                 zlim=-z_lim(ri0,zj0+1)*correct_plus
              endif
              delta_zj=1
           endif
        else
           ! on descend
           if (zj0==1) then
              ! on traverse le plan eq donc on va remonter
              ! et z va changer de signe
              delta_zj=1
              if (z0 > 0.0_db) then
                 zlim=-z_lim(ri0,2)*correct_moins
              else
                 zlim=z_lim(ri0,2)*correct_moins
              endif
           else !(zj0==1)
              ! on ne traverse pas z=0.
              if (z0 > 0.0_db) then
                 zlim=z_lim(ri0,zj0)*correct_moins
              else
                 zlim=-z_lim(ri0,zj0)*correct_moins
              endif
              delta_zj=-1
           endif !(zj0==1)
        endif ! monte ou descend
        t=(zlim-z0)*inv_w
        ! correct pb precision
        if (t < 0.0_db) t=prec_grille
     endif !dotprod=0.0
  endif ! ri0==0


  ! 3) interface en r ou z ?
  if (s < t) then ! r
     l=s
     delta_vol=s
     ! Position au bord de la cellule suivante
     x1=x0+delta_vol*u
     y1=y0+delta_vol*v
     z1=z0+delta_vol*w
     ri1=ri0+delta_rad
     if ((ri1<1).or.(ri1>n_rad)) then
        zj1=zj0
     else
        zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
        if (zj1>nz) zj1=nz+1
     endif
  else ! z
     l=t
     delta_vol=t
     ! Position au bord de la cellule suivante
     x1=x0+delta_vol*u
     y1=y0+delta_vol*v
     z1=z0+delta_vol*w
     ri1=ri0
     zj1=zj0+delta_zj
  endif

  ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
  if (z1 == 0.0_db) z1 = prec_grille

  ! Calcul longeur de vol et profondeur optique dans la cellule
  tau = l*opacite ! opacite constante dans la cellule

  return

end subroutine cross_cylindrical_cell_tmp

!*************************************************************************************

subroutine save_radiation_field(id,lambda,icell0, Stokes, l,  x0,y0,z0, x1,y1,z1, u,v, w, flag_star, flag_direct_star)

  integer, intent(in) :: id,lambda,icell0
  real(kind=db), dimension(4), intent(in) :: Stokes
  real(kind=db) :: l, x0,y0,z0, x1,y1,z1, u,v,w
  logical, intent(in) :: flag_star, flag_direct_star


  real(kind=db) :: xm,ym,zm, phi_pos, phi_vol
  integer :: psup, phi_I, theta_I, phi_k

  integer :: ri0, zj0, k0


 ! 3D cell indices : TMP
  call cell2cylindrical(icell0, ri0,zj0,k0)

  phi_vol = atan2(v,u) + deux_pi ! deux_pi pour assurer diff avec phi_pos > 0

  if (letape_th) then
     if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + kappa_abs_eg(icell0,lambda) &
          * l * Stokes(1)
     if (lxJ_abs) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
  else
     if (lProDiMo) then
        xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
        ! Pour statistique: nbre de paquet contribuant a intensite specifique
        xN_abs(icell0,lambda,id) = xN_abs(icell0,lambda,id) + 1.0
     endif ! lProDiMo

     if (lscatt_ray_tracing1) then
        xm = 0.5_db * (x0 + x1)
        ym = 0.5_db * (y0 + y1)
        zm = 0.5_db * (z0 + z1)

        phi_pos = atan2(ym,xm)
        phi_k = floor(  modulo(phi_pos, deux_pi) / deux_pi * n_az_rt ) + 1
        if (phi_k > n_az_rt) phi_k=n_az_rt

        if (zm > 0.0_db) then
           psup = 1
        else
           psup = 2
        endif

        if (lsepar_pola) then
           call calc_xI_scatt_pola(id,lambda,icell0,phi_k,psup,l,Stokes(:),flag_star) ! phik & psup=1 in 3D
        else
           ! ralentit d'un facteur 5 le calcul de SED
           ! facteur limitant
           call calc_xI_scatt(id,lambda,icell0,phi_k,psup,l,Stokes(1),flag_star)
        endif

     else if (lscatt_ray_tracing2) then
        if (flag_direct_star) then
           I_spec_star(icell0,id) = I_spec_star(icell0,id) + l * Stokes(1)
        else
           xm = 0.5_db * (x0 + x1)
           ym = 0.5_db * (y0 + y1)
           zm = 0.5_db * (z0 + z1)
           phi_pos = atan2(ym,xm)

           !  if (l_sym_ima) then
           !     delta_phi = modulo(phi_vol - phi_pos, deux_pi)
           !     if (delta_phi > pi) delta_phi = deux_pi - delta_phi
           !     phi_I =  nint( delta_phi  / pi * (n_phi_I -1) ) + 1
           !     if (phi_I > n_phi_I) phi_I = n_phi_I
           !  else
           phi_I =  floor(  modulo(phi_vol - phi_pos, deux_pi) / deux_pi * n_phi_I ) + 1
           if (phi_I > n_phi_I) phi_I = 1
           !  endif

           if (zm > 0.0_db) then
              theta_I = floor(0.5_db*( w + 1.0_db) * n_theta_I) + 1
           else
              theta_I = floor(0.5_db*(-w + 1.0_db) * n_theta_I) + 1
           endif
           if (theta_I > n_theta_I) theta_I = n_theta_I

           I_spec(1:n_Stokes,theta_I,phi_I,icell0,id) = I_spec(1:n_Stokes,theta_I,phi_I,icell0,id) + l * Stokes(1:n_Stokes)

           if (lsepar_contrib) then
              if (flag_star) then
                 I_spec(n_Stokes+2,theta_I,phi_I,icell0,id) = I_spec(n_Stokes+2,theta_I,phi_I,icell0,id) + l * Stokes(1)
              else
                 I_spec(n_Stokes+4,theta_I,phi_I,icell0,id) = I_spec(n_Stokes+4,theta_I,phi_I,icell0,id) + l * Stokes(1)
              endif
           endif ! lsepar_contrib

        endif ! flag_direct_star
     endif !lscatt_ray_tracing
  endif !letape_th

  return

end subroutine save_radiation_field

!*************************************************************************************


subroutine length_deg2_sph(id,lambda,Stokes,ri,thetaj,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
! Integration par calcul de la position de l'interface entre cellules
! par eq deg2 en r et deg 1 en z
! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
! C. Pinte
! 03/06/07

  implicit none

  integer, intent(in) :: id,lambda
  integer, intent(inout) :: ri,thetaj
  real(kind=db), dimension(4), intent(in) :: Stokes
  logical, intent(in) :: flag_star,flag_direct_star
 real(kind=db), intent(inout) :: u,v,w
  real, intent(in) :: extrin
  real(kind=db), intent(inout) :: xio,yio,zio
  real, intent(out) :: ltot
  logical, intent(out) :: flag_sortie

  real(kind=db) :: x0, y0, z0, x1, y1, z1, x_old, y_old, z_old, factor
  real(kind=db) :: b, c, s, rac, t, delta, r0_2, r0_cyl, r0_2_cyl
  real(kind=db) :: delta_vol, l, tau, extr, dotprod, opacite
  real(kind=db) :: correct_moins, correct_plus, uv, precision
  real(kind=db) :: phi_pos, phi_vol, delta_phi, xm, ym, zm
  integer :: ri0, thetaj0, ri1, thetaj1, delta_rad, delta_theta, nbr_cell, p_ri0, p_thetaj0, icell0
  integer :: theta_I, phi_I, thetaj_old, ri_old, icell
  logical :: lcellule_non_vide, lstop

  real(kind=db) :: a_theta, b_theta, c_theta, tan2, tan_angle_lim1, tan_angle_lim2, t1, t2
  real(kind=db) :: t1_1, t1_2, t2_1, t2_2, a_theta_m1

  logical :: detect_bug = .false.

!  detect_bug =  (abs(xio - 134.48909960993615_db ) < 1.0e-10)
  detect_bug=.false.

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille_sph
  correct_plus = 1.0_db + prec_grille_sph
  precision = 1.0e-20_db ! pour g95
  lstop = .false.

  x0=xio;y0=yio;z0=zio
  x1=xio;y1=yio;z1=zio
  extr=extrin
  ri0=ri
  thetaj0=thetaj
  ri1=ri
  thetaj1=thetaj

  ! Cas sans strat
  p_ri0=1
  p_thetaj0=1

  ! TODO : cas avec strat

  ltot=0.0
  nbr_cell = 0

  uv = sqrt(u*u + v*v)

  ! pour ray-tracing
  phi_vol = atan2(v,u) + deux_pi ! deux_pi pour assurer diff avec phi_pos > 0

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri_old = ri0 ; thetaj_old = thetaj0
     x_old = x0 ; y_old = y0 ; z_old = z0

     ri0=ri1;thetaj0=thetaj1
     x0=x1;y0=y1;z0=z1

     icell0 = cell_map(ri0,thetaj0,1)

     ! Pour cas avec approximation de diffusion
     if (l_dark_zone(icell0)) then
        ! On revoie le paquet dans l'autre sens
        u = -u ; v = -v ; w=-w
        ! et on le renvoie au point de depart
        ri = ri_old ; thetaj = thetaj_old
        xio = x_old ; yio = y_old ; zio = z_old
        return
     endif

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        flag_sortie = .true.
        return
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r0_2_cyl  = x0*x0+y0*y0
     r0_cyl = sqrt(r0_2_cyl)
     r0_2 = r0_2_cyl + z0*z0
     b   = (x0*u+y0*v+z0*w)

     if (ri0==0) then
        lcellule_non_vide=.false.
        opacite=0.0_db
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r0_2-r_lim_2(0)*correct_plus)
        delta=b*b-c
        rac=sqrt(delta)
        s = (-b+rac) * correct_plus
        t=huge_real
        delta_rad=1
     else
        if (icell0 > n_cells) then
           opacite = 0.0_db
        else
           opacite=kappa(icell0,lambda)
        endif

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod= b
        if (dotprod < 0.0_db) then
           ! on recule : on cherche rayon inférieur
           c=(r0_2-r_lim_2(ri0-1)*correct_moins)
           delta=b*b-c
           if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r0_2-r_lim_2(ri0)*correct_plus)
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r0_2-r_lim_2(ri0)*correct_plus)
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0_db) then
           s=-b+rac
        else if (s==0.0_db) then
           s=prec_grille
        endif

        ! 2) position interface inclinaison
        if (z0 >= 0.0_db) then ! on est dans le bon sens
           tan_angle_lim1 = tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = tan_theta_lim(thetaj0-1) * correct_moins
        else
           tan_angle_lim1 = - tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = - tan_theta_lim(thetaj0-1) * correct_moins
        endif ! z0

        ! Premiere limite theta
        tan2 = tan_angle_lim1 * tan_angle_lim1
        a_theta = w*w - tan2 * (u*u + v*v)
        a_theta_m1 = 1.0_db/a_theta
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t1 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t1_1 = (- b_theta - rac) * a_theta_m1
           t1_2 = (- b_theta + rac) * a_theta_m1

           if (t1_1 <= precision) then
              if (t1_2 <= precision) then ! les 2 sont <0
                 t1=1.0e30_db
              else   ! Seul t1_2 > 0
                 t1 = t1_2
              endif
           else
              if (t1_2 <= precision) then ! Seul t1_1 >0
                 t1 = t1_1
              else ! les 2 sont > 0
                 t1 = min(t1_1,t1_2)
              endif
           endif
        endif ! signe delta

        ! Deuxieme limite theta
        tan2 = tan_angle_lim2 * tan_angle_lim2
        a_theta = w*w - tan2 * (u*u + v*v)
        a_theta_m1 = 1.0_db/a_theta
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t2 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t2_1 = (- b_theta - rac) * a_theta_m1
           t2_2 = (- b_theta + rac) * a_theta_m1

           if (t2_1 <= precision) then
              if (t2_2 <= precision) then ! les 2 sont <0
                 t2=1.0e30_db
              else   ! Seul t2_2 > 0
                 t2 = t2_2
              endif
           else
              if (t2_2 <= precision) then ! Seul t2_1 >0
                 t2 = t2_1
              else ! les 2 sont > 0
                 t2 = min(t2_1,t2_2)
              endif
           endif
        endif ! signe delta

        ! Selection limite theta
        if (t1 < t2) then
           t=t1
           delta_theta = 1
           if (thetaj0 == nz) delta_theta = 0
        else
           t=t2
           delta_theta = -1
           if (thetaj0 == 1) delta_theta = 0
        endif

     endif ! ri0==0


     ! 3) interface en r ou theta ?
     if (s < t) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        thetaj1 = thetaj0

        ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
        if (ri0==0) then
           call indice_cellule_sph_theta(x1,y1,z1,thetaj1)
        endif

     else ! theta
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        thetaj1=thetaj0+delta_theta
     endif

     ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
     if (z1 == 0.0_db) z1 = prec_grille

     ! Calcul longeur de vol et profondeur optique dans la cellule
     tau=l*opacite ! opacite constante dans la cellule
     ! Comparaison integrale avec tau
     ! et ajustement longueur de vol eventuellement
     if(tau > extr) then ! On a fini d'integrer
        lstop = .true.
        l = l*extr/tau ! on rescale l pour que tau=extr
        ltot=ltot+l
     else ! Il reste extr - tau a integrer dans la cellule suivante
        extr=extr-tau
        ltot=ltot+l
     endif

     ! Stokage des champs de radiation
     if (lcellule_non_vide) then
        if (letape_th) then
           if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + kappa_abs_eg(icell0,lambda) * l * Stokes(1)
           if (lxJ_abs) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
        else if (lscatt_ray_tracing2) then
           if (flag_direct_star) then
              I_spec_star(icell0,id) = I_spec_star(icell0,id) + l * Stokes(1)
           else
              xm = 0.5_db * (x0 + x1)
              ym = 0.5_db * (y0 + y1)
              zm = 0.5_db * (z0 + z1)

              if (l_sym_ima) then
                 phi_pos = atan2(ym,xm)
                 delta_phi = modulo(phi_vol - phi_pos, deux_pi)
                 if (delta_phi > pi) delta_phi = abs(delta_phi - deux_pi)
                 phi_I =  nint( delta_phi  / pi * (n_phi_I ) ) + 1
                 if (phi_I > n_phi_I) phi_I = n_phi_I
              else
                 phi_pos = atan2(ym,xm)
                 phi_I =  nint(  modulo(phi_vol - phi_pos, deux_pi) / deux_pi * n_phi_I ) + 1
                 if (phi_I > n_phi_I) phi_I = 1
              endif

              if (zm > 0._db) then
                 theta_I = nint(0.5_db*( w + 1.0_db) * n_theta_I) + 1
              else
                 theta_I = nint(0.5_db*(-w + 1.0_db) * n_theta_I) + 1
              endif
              if (theta_I > n_theta_I) theta_I = n_theta_I

              I_spec(1:n_Stokes,theta_I,phi_I,icell0,id) = I_spec(1:n_Stokes,theta_I,phi_I,icell0,id) + l * Stokes(1:n_Stokes)

              if (lsepar_contrib) then
                 if (flag_star) then
                    I_spec(n_Stokes+2,theta_I,phi_I,icell0,id) = I_spec(n_Stokes+2,theta_I,phi_I,icell0,id) + l * Stokes(1)
                 else
                    I_spec(n_Stokes+4,theta_I,phi_I,icell0,id) = I_spec(n_Stokes+4,theta_I,phi_I,icell0,id) + l * Stokes(1)
                 endif
              endif ! lsepar_contrib

           endif

        endif

     endif ! lcellule_non_vide

     ! On a fini d'integrer
     if (lstop) then
        flag_sortie = .false.
        xio=x0+l*u
        yio=y0+l*v
        zio=z0+l*w
        call indice_cellule_sph(xio,yio,zio,ri,thetaj)

        ! Patch pour eviter BUG sur position radiale
        ! a cause de limite de precision
        if (ri==0) then
           factor = rmin/ sqrt(xio*xio+yio*yio+zio*zio) * correct_plus
           xio = xio * factor
           yio = yio * factor
           zio = zio * factor

           ! On verifie que c'est OK maintenant
           call indice_cellule_sph(xio,yio,zio,ri,thetaj)
           if (ri==0) then
              write(*,*) "BUG integ_deg2_SPH"
              write(*,*) "Exiting"
              stop
           endif
        endif

      !  if (lmono) then
      !     write(*,*) "xio", real(xio),real(yio),real(zio)
      !     write(*,*), "r", ri0, thetaj0, ri,thetaj
      !  endif
        icell = cell_map(ri,thetaj,1)
        if (l_dark_zone(icell)) then ! Petit test de securite
           ! On resort le paquet
           if (thetaj < thetaj0) then
              thetaj = thetaj0
              zio = z_lim(ri0,thetaj0)*correct_plus
           endif
           if (ri < ri0) then
              ri = ri0
              xio = xio * correct_plus
              yio = yio * correct_plus
           else if (ri > ri0) then
              ri = ri0
              xio = xio * correct_moins
              yio = yio * correct_moins
           endif
        endif
        return
     endif ! lstop


  enddo ! boucle infinie

  write(*,*) "BUG"
  return

end subroutine length_deg2_sph

!********************************************************************

subroutine length_deg2_3D(id,lambda,Stokes,ri,zj,phik,xio,yio,zio,u,v,w,flag_star,extrin,ltot,flag_sortie)
! Integration par calcul de la position de l'interface entre cellules
! par eq deg2 en r et deg 1 en z
! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
! C. Pinte
! 05/02/05

  implicit none

  integer, intent(in) :: id,lambda
  integer, intent(inout) :: ri,zj, phik
  real(kind=db), dimension(4), intent(in) :: Stokes
  real(kind=db), intent(in) :: u,v,w
  logical, intent(in) :: flag_star
  real, intent(in) :: extrin
  real(kind=db), intent(inout) :: xio,yio,zio
  real, intent(out) :: ltot
  logical, intent(out) :: flag_sortie

  real(kind=db) :: x0, y0, z0, x1, y1, z1, phi
  real(kind=db) :: inv_a, a, b, c, s, rac, t, t_phi, delta, inv_w, r_2, tan_angle_lim, den
  real(kind=db) :: delta_vol, l, tau, zlim, extr, dotprod, opacite
  real(kind=db) :: correct_plus, correct_moins
  real(kind=db) :: xm, ym, zm
  integer :: ri0, zj0, ri1, zj1, phik0, phik1, delta_rad, delta_zj, nbr_cell, delta_phi, phik0m1, icell0

  logical :: lcellule_non_vide

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  x1=xio;y1=yio;z1=zio
  extr=extrin
  ri0=ri
  zj0=zj
  phik0=phik
  ri1=ri
  zj1=zj
  phik1=phik

  ltot=0.0
  nbr_cell = 0

  a=u*u+v*v

  if (a > tiny_real) then
     inv_a=1.0/a
  else
     inv_a=huge_real
  endif

  if (abs(w) > tiny_real) then
     inv_w=1.0/w
  else
     inv_w=sign(huge_db,w)
  endif

  ! Calcule les angles de diffusion pour la direction de propagation donnee
  if ((.not.letape_th).and.lscatt_ray_tracing1) call angles_scatt_rt1(id,u,v,w)

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1 ; zj0=zj1 ; phik0=phik1
     x0=x1 ; y0=y1 ; z0=z1

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        flag_sortie = .true.
        return
     elseif (abs(zj0) > nz) then
        lcellule_non_vide=.false.
        ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           flag_sortie = .true.
           return
        endif
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r_2=x0*x0+y0*y0
     b=(x0*u+y0*v)*inv_a

     if (ri0==0) then
        lcellule_non_vide=.false.
        opacite=0.0
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r_2-r_lim_2(0)*correct_plus)*inv_a
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        t_phi= huge_real
        delta_rad=1
     else
        icell0 = cell_map(ri0,zj0,phik0)
        if (icell0 > n_cells) then
           opacite = 0.0_db
        else
           opacite=kappa(icell0,lambda)
        endif
        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod=u*x0+v*y0
        if (dotprod < 0.0) then
           ! on recule : on cherche rayon inférieur
           c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
           delta=b*b-c
           if (delta < 0.0) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0) then
           s=-b+rac
        else if (s==0.0) then
           s=prec_grille
        endif

        ! 2) position interface verticale
        ! on monte ou on descend par plan équatorial ?
        dotprod=w*z0
        if (dotprod == 0.0) then
           t=1.0e10
        else
           if (dotprod > 0.0) then
              ! on s'éloigne du midplane
              if (zj0==nz+1) then
                 delta_zj=0
                 zlim=1.0e10
              else if (zj0==-(nz+1)) then
                 delta_zj=0
                 zlim=-1.0e10
              else
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,zj0+1)*correct_plus
                    delta_zj=1
                 else
                    zlim=-z_lim(ri0,abs(zj0)+1)*correct_plus
                    delta_zj=-1
                 endif
              endif
           else
              ! on se rapproche du midplane
              if (z0 > 0.0) then
                 zlim=z_lim(ri0,abs(zj0))*correct_moins
                 delta_zj=-1
                 if (zj0==1) delta_zj=-2 ! pas d'indice 0
              else
                 zlim=-z_lim(ri0,abs(zj0))*correct_moins
                 delta_zj=1
                 if (zj0==-1) delta_zj=2 ! pas d'indice 0
              endif
           endif ! monte ou descend
           t=(zlim-z0)*inv_w
           ! correct pb precision
           if (t < 0.0) t=prec_grille
        endif !dotprod=0.0


        ! 3) position interface azimuthale
        dotprod =  x0*v - y0*u
        if (abs(dotprod) < 1.0e-10) then
           ! on ne franchit pas d'interface azimuthale
           t_phi = 1.0e30
        else
           ! Quelle cellule on va franchir
           if (dotprod > 0.0) then
              tan_angle_lim = tan_phi_lim(phik0)
              delta_phi=1
           else
              phik0m1=phik0-1
              if (phik0m1==0) phik0m1=N_az
              tan_angle_lim = tan_phi_lim(phik0m1)
              delta_phi=-1
           endif
           ! Longueur av interserction
           if (tan_angle_lim > 1.0d299) then
              t_phi = -x0/u
           else
              den= v-u*tan_angle_lim
              if (abs(den) > 1.0e-6) then
                 t_phi = -(y0-x0*tan_angle_lim)/den
              else
                 t_phi = 1.0e30
              endif
           endif
           if (t_phi < 0.0) t_phi = 1.0e30
        endif !dotprod = 0.0

     endif ! ri==0

     ! 4) interface en r ou z ou phi ?
     if ((s < t).and.(s < t_phi)) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        if ((ri1<1).or.(ri1>n_rad)) then
           zj1=zj0
        else
           zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
           if (zj1>nz) zj1=nz+1
           if (z1 < 0.0) zj1=-zj1
        endif

        phik1=phik0
        ! We need to find the azimuth when we enter the disc
        ! It can be different from the initial azimuth if the star is not centered
        ! so we need to compute it here
        if (ri0==0) then
           phi=modulo(atan2(y1,x1),2*real(pi,kind=db))
           phik1=floor(phi*un_sur_deux_pi*real(N_az))+1
           if (phik1==n_az+1) phik1=n_az
        endif

     else if (t < t_phi) then ! z
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1=zj0+delta_zj
        phik1=phik0
     else
        l=t_phi
        delta_vol=correct_plus*t_phi
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1= floor(abs(z1)/zmax(ri1)*nz) + 1
        if (zj1>nz) zj1=nz+1
        if (z1 < 0.0) zj1=-zj1
        phik1=phik0+delta_phi
        if (phik1 == 0) phik1=N_az
        if (phik1 == N_az+1) phik1=1
     endif

     ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
     if (z1 == 0.0_db) z1 = sign(prec_grille,w)

     ! Calcul longeur de vol et profondeur optique dans la cellule
     tau=l*opacite ! opacite constante dans la cellule
     ! Comparaison integrale avec tau
     if(tau > extr) then ! On a fini d'integrer
        l = l*extr/tau ! on rescale l pour que tau=extr
        ltot=ltot+l
        if (letape_th.and.lcellule_non_vide) then
           if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + kappa_abs_eg(icell0,lambda) * l * Stokes(1)
           if (lxJ_abs) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
        endif !l_abs
        flag_sortie = .false.
        xio=x0+l*u
        yio=y0+l*v
        zio=z0+l*w
        call indice_cellule_3D(xio,yio,zio,ri,zj,phik)
        return
     else ! Il reste extr - tau a integrer dans la cellule suivante
        extr=extr-tau
        ltot=ltot+l

        if (lcellule_non_vide) then
           if (letape_th) then
              if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + &
                   kappa_abs_eg(icell0,lambda) * l * Stokes(1)
              if (lxJ_abs) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
           else ! letape_th
              if (lscatt_ray_tracing1) then
                 xm = 0.5_db * (x0 + x1)
                 ym = 0.5_db * (y0 + y1)
                 zm = 0.5_db * (z0 + z1)

                 ! TODO : RT : clean the 0 in the indices
                 if (lsepar_pola) then
                    call calc_xI_scatt_pola(id,lambda,icell0,1,1,l,Stokes(:),flag_star) ! phik & psup=1 in 3D
                 else
                    call calc_xI_scatt(id,lambda,icell0,1,1,l,Stokes(1),flag_star) ! phik & psup=1 in 3D
                 endif
              endif

           endif ! letape_th
        endif ! lcellule_non_vide

     endif ! tau > extr

  enddo ! boucle infinie
  write(*,*) "BUG"
  return

end subroutine length_deg2_3D

!********************************************************

subroutine integ_tau(lambda)

  implicit none

  integer, intent(in) :: lambda

  real :: norme
  integer :: i, ri, zj, j, icell

  real(kind=db), dimension(4) :: Stokes
  ! angle de visee en deg
  real :: angle
  real(kind=db) :: x0, y0, z0, u0, v0, w0
  real :: tau
  real(kind=db) :: lmin, lmax

  angle=angle_interet

  norme=0.0
  do i=1, n_rad
     icell = cell_map(i,1,1)
     norme=norme+kappa(icell,lambda)*(r_lim(i)-r_lim(i-1))
  enddo

  norme = norme
  write(*,*) 'Integ tau dans plan eq. = ', norme
  ! 1.49597870691e13 car kappa est en AU**-1
  ! Ok si pas de sedimentation
  if (.not.lvariable_dust) then
     icell = cell_map(1,1,1)
     if (kappa(icell,lambda) > tiny_real) then
        write(*,*) " Column density (g/cm²)   = ", real(norme*(masse(icell)/(volume(1)*AU_to_cm**3))/ &
             (kappa(icell,lambda)/AU_to_cm))
     endif
  endif


  norme=0.0
  do j=1, nz
     icell = cell_map(1,j,1)
     norme=norme+kappa(icell,lambda)*(z_lim(1,j+1)-z_lim(1,j))
  enddo
  norme = norme * 2.

  write(*,*) 'Integ tau vert = ', norme
  ! 1.49597870691e13 car kappa est en AU**-1
  ! Ok si pas de sedimentation
  if (.not.lvariable_dust) then
     icell = cell_map(1,1,1)
     if (kappa(icell,lambda) > tiny_real) then
        write(*,*) " Column density (g/cm²)   = ", real(norme*(masse(icell)/(volume(1)*AU_to_cm**3))/ &
             (kappa(icell,lambda)/AU_to_cm))
     endif
  endif

  ri=0 ; zj=1 ; x0=0.0 ; y0=0.0 ; z0=0.0
  Stokes = 0.0_db ; Stokes(1) = 1.0_db
  w0 = cos((angle)*pi/180.)
  u0 = sqrt(1.0-w0*w0)
  v0 = 0.0

  call length_deg2_tot(1,lambda,Stokes,ri,zj,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)

  write(*,fmt='(" Integ tau (i =",f4.1," deg)   = ",E12.5)') angle, tau

  if (.not.lvariable_dust) then
     icell = cell_map(1,1,1)
     if (kappa(icell,lambda) > tiny_real) then
        write(*,*) " Column density (g/cm²)   = ", real(tau*(masse(icell)/(volume(1)*3.347929d39))/ &
             (kappa(icell,lambda)/1.49597870691e13))
     endif
  endif

  if (norme > tau_seuil) then
     write(*,*) "Tau greater than", tau_seuil, "Exiting"
     stop
  endif


  ! Verif
!!$  do i=1,30
!!$     tau=0.5*(tau_max+tau_min)
!!$     ri=0 ; zj=1 ; x0=0.0 ; y0=0.0 ; z0=0.0
!!$     call length_deg2(1,lambda,1,ri,zj,x0,y0,y0,u0,v0,w0,tau,dvol1,flag_sortie)
!!$     if (flag_sortie==4) then
!!$        tau_max=tau
!!$     else
!!$        tau_min=tau
!!$     endif
!!$  enddo
!!$  write(*,fmt='(" Integ tau (i =",f4.1," deg)   = ",F12.2)') angle, 0.5*(tau_max+tau_min)

  return

end subroutine integ_tau

!***********************************************************

subroutine length_deg2_tot(id,lambda,Stokes,ri,zj,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
! Integration par calcul de la position de l'interface entre cellules
! de l'opacite totale dans une direction donnée
! C. Pinte
! 19/06/07

  implicit none

  integer, intent(in) :: id,lambda, ri,zj
  real(kind=db), dimension(4), intent(in) :: Stokes
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: xi,yi,zi
  real, intent(out) :: tau_tot_out
  real(kind=db), intent(out) :: lmin,lmax

  if (lcylindrical) then
     call length_deg2_tot_cyl(id,lambda,Stokes,ri,zj,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
  else
     call length_deg2_tot_sph(id,lambda,Stokes,ri,zj,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
  endif

  return

end subroutine length_deg2_tot

!***********************************************************

subroutine length_deg2_tot_cyl(id,lambda,Stokes,ri,zj,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
! Integration par calcul de la position de l'interface entre cellules
! de l'opacite totale dans une direction donnée
! Grille a geometrie cylindrique
! C. Pinte
! 19/04/05

  implicit none

  integer, intent(in) :: id,lambda, ri,zj
  real(kind=db),dimension(4), intent(in) :: Stokes
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: xi,yi,zi
  real, intent(out) :: tau_tot_out
  real(kind=db), intent(out) :: lmin,lmax


  real(kind=db) :: x0, y0, z0, x1, y1, z1
  real(kind=db) :: inv_a, a, b, c, s, rac, t, delta, inv_w, r_2
  real(kind=db) :: delta_vol, l, ltot, tau, zlim, dotprod, opacite, tau_tot
  real(kind=db) :: correct_plus, correct_moins
  integer :: ri0, zj0, ri1, zj1, delta_rad, delta_zj, nbr_cell, icell0

  correct_plus = 1.0_db + prec_grille
  correct_moins = 1.0_db - prec_grille

  x1=xi;y1=yi;z1=zi
  ri0=ri
  zj0=zj
  ri1=ri
  zj1=zj

  tau_tot=0.0_db
  nbr_cell = 0

  a=u*u+v*v

  if (a > tiny(0.)) then
     inv_a=1.0/a
  else
     inv_a=huge_real
  endif

  if (abs(w) > tiny_real) then
     inv_w=1.0/w
  else
     inv_w=sign(huge_db,w)
  endif


  lmin=0.0_db
  ltot=0.0_db

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1;zj0=zj1
     x0=x1;y0=y1;z0=z1

     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        tau_tot_out=tau_tot
        lmax=ltot
        return
     elseif (zj0>nz) then
        ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           tau_tot_out=tau_tot
           lmax=ltot
           return
        endif
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r_2=x0*x0+y0*y0
     b=(x0*u+y0*v)*inv_a

     if (ri0==0) then
        opacite=0.0
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r_2-r_lim_2(0)*correct_plus)*inv_a
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        delta_rad=1
     else
        icell0 = cell_map(ri0,zj0,1)
        if (icell0 > n_cells) then
           opacite = 0.0_db
        else
           opacite=kappa(icell0,lambda)
        endif
        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod=u*x0+v*y0
        if (dotprod < 0.0) then
           ! on recule : on cherche le rayon inférieur
           c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
           delta=b*b-c
           if (delta < 0.0) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0) then ! test a faire car on croise le rayon inferieur 2 fois
           s=-b+rac
       else if (s==0.0) then
           s=prec_grille
        endif

        ! 2) position interface verticale
        ! on monte ou on descend par plan équatorial ?
        dotprod=w*z0
        if (abs(w) < tiny_real) then
           t=huge_real
        else
           if (dotprod > 0.0) then
              ! on monte
              if (zj0==nz+1) then
                 delta_zj=0
                 if (z0 > 0.0) then
                    zlim=1.0e10
                 else
                    zlim=-1.0e10
                 endif
              else
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,zj0+1)*correct_plus
                 else
                    zlim=-z_lim(ri0,zj0+1)*correct_plus
                 endif
                 delta_zj=1
              endif
           else
              ! on descend
              if (zj0==1) then
                 ! on traverse le plan eq donc on va remonter
                 ! et z va changer de signe
                 delta_zj=1
                 if (z0 > 0.0) then
                    zlim=-z_lim(ri0,2)*correct_moins
                 else
                    zlim=z_lim(ri0,2)*correct_moins
                 endif
              else !(zj0==1)
                 ! on ne traverse pas z=0.
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,zj0)*correct_moins
                 else
                    zlim=-z_lim(ri0,zj0)*correct_moins
                 endif
                 delta_zj=-1
              endif !(zj0==1)
           endif ! monte ou descend
           t=(zlim-z0)*inv_w
           ! correct pb precision
           if (t < 0.0) t=prec_grille
        endif !dotprod=0.0
     endif ! ri0==0


     ! 3) interface en r ou z ?
     if (s < t) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        if ((ri1<1).or.(ri1>n_rad)) then
           zj1=zj0
        else
           zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
           if (zj1>nz) zj1=nz+1
        endif
     else ! z
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1=zj0+delta_zj
     endif

     ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
     if (z1 == 0.0_db) z1 = prec_grille

     tau=l*opacite ! opacite constante dans la cellule

     tau_tot = tau_tot + tau
     ltot= ltot + l
     if (tau_tot < tiny_real) lmin=ltot

  enddo ! boucle infinie

  write(*,*) "BUG"
  return

end subroutine length_deg2_tot_cyl

!***********************************************************

subroutine length_deg2_tot_3D(id,lambda,Stokes,ri,zj,phik,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
! C. Pinte
! 14/07/14

  implicit none

  integer, intent(in) :: id,lambda, ri,zj, phik
  real(kind=db),dimension(4), intent(in) :: Stokes
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: xi,yi,zi
  real, intent(out) :: tau_tot_out
  real(kind=db), intent(out) :: lmin,lmax

  real(kind=db) :: x0, y0, z0, x1, y1, z1, phi
  real(kind=db) :: inv_a, a, b, c, s, rac, t, t_phi, delta, inv_w, r_2, tan_angle_lim, den
  real(kind=db) :: delta_vol, l, ltot, tau, zlim, dotprod, opacite, tau_tot
  real(kind=db) :: correct_plus, correct_moins
  integer :: ri0, zj0, ri1, zj1, phik0, phik1, delta_rad, delta_zj, nbr_cell, delta_phi, phik0m1, icell0

  logical :: lcellule_non_vide

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  x1=xi;y1=yi;z1=zi
  ri0=ri
  zj0=zj
  ri1=ri
  zj1=zj
  phik1=phik

  tau_tot=0.0_db
  nbr_cell = 0

  a=u*u+v*v

  if (a > tiny(0.)) then
     inv_a=1.0/a
  else
     inv_a=huge_real
  endif

  if (abs(w) > tiny_real) then
     inv_w=1.0/w
  else
     inv_w=sign(huge_db,w)
  endif


  lmin=0.0_db
  ltot=0.0_db

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1 ; zj0=zj1 ; phik0=phik1
     x0=x1 ; y0=y1 ; z0=z1

     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        tau_tot_out=tau_tot
        lmax=ltot
        return
     elseif (abs(zj0) > nz) then
         ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           tau_tot_out=tau_tot
           lmax=ltot
           return
        endif
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r_2=x0*x0+y0*y0
     b=(x0*u+y0*v)*inv_a

     if (ri0==0) then
        lcellule_non_vide=.false.
        opacite=0.0
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r_2-r_lim_2(0)*correct_plus)*inv_a
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        t_phi= huge_real
        delta_rad=1
     else
        icell0 = cell_map(ri0,zj0,phik0)
        if (icell0 > n_cells) then
           opacite = 0.0_db
        else
           opacite=kappa(icell0,lambda)
        endif
        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod=u*x0+v*y0
        if (dotprod < 0.0) then
           ! on recule : on cherche rayon inférieur
           c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
           delta=b*b-c
           if (delta < 0.0) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0) then
           s=-b+rac
        else if (s==0.0) then
           s=prec_grille
        endif

        ! 2) position interface verticale
        ! on monte ou on descend par plan équatorial ?
        dotprod=w*z0
        if (dotprod == 0.0) then
           t=1.0e10
        else
           if (dotprod > 0.0) then
              ! on s'éloigne du midplane
              if (zj0==nz+1) then
                 delta_zj=0
                 zlim=1.0e10
              else if (zj0==-(nz+1)) then
                 delta_zj=0
                 zlim=-1.0e10
              else
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,zj0+1)*correct_plus
                    delta_zj=1
                 else
                    zlim=-z_lim(ri0,abs(zj0)+1)*correct_plus
                    delta_zj=-1
                 endif
              endif
           else
              ! on se rapproche du midplane
              if (z0 > 0.0) then
                 zlim=z_lim(ri0,abs(zj0))*correct_moins
                 delta_zj=-1
                 if (zj0==1) delta_zj=-2 ! pas d'indice 0
              else
                 zlim=-z_lim(ri0,abs(zj0))*correct_moins
                 delta_zj=1
                 if (zj0==-1) delta_zj=2 ! pas d'indice 0
              endif
           endif ! monte ou descend
           t=(zlim-z0)*inv_w
           ! correct pb precision
           if (t < 0.0) t=prec_grille
        endif !dotprod=0.0


        ! 3) position interface azimuthale
        dotprod =  x0*v - y0*u
        if (abs(dotprod) < 1.0e-10) then
           ! on ne franchit pas d'interface azimuthale
           t_phi = 1.0e30
        else
           ! Quelle cellule on va franchir
           if (dotprod > 0.0) then
              tan_angle_lim = tan_phi_lim(phik0)
              delta_phi=1
           else
              phik0m1=phik0-1
              if (phik0m1==0) phik0m1=N_az
              tan_angle_lim = tan_phi_lim(phik0m1)
              delta_phi=-1
           endif
           ! Longueur av interserction
           if (tan_angle_lim > 1.0d299) then
              t_phi = -x0/u
           else
              den= v-u*tan_angle_lim
              if (abs(den) > 1.0e-6) then
                 t_phi = -(y0-x0*tan_angle_lim)/den
              else
                 t_phi = 1.0e30
              endif
           endif
           if (t_phi < 0.0) t_phi = 1.0e30
        endif !dotprod = 0.0

     endif ! ri==0

     ! 4) interface en r ou z ou phi ?
     if ((s < t).and.(s < t_phi)) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        if ((ri1<1).or.(ri1>n_rad)) then
           zj1=zj0
        else
           zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
           if (zj1>nz) zj1=nz+1
           if (z1 < 0.0) zj1=-zj1
        endif

        phik1=phik0
        ! We need to find the azimuth when we enter the disc
        ! It can be different from the initial azimuth if the star is not centered
        ! so we need to compute it here
        if (ri0==0) then
           phi=modulo(atan2(y1,x1),2*real(pi,kind=db))
           phik1=floor(phi/(2*pi)*real(N_az))+1
           if (phik1==n_az+1) phik1=n_az
        endif

     else if (t < t_phi) then ! z
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1=zj0+delta_zj
        phik1=phik0
     else
        l=t_phi
        delta_vol=correct_plus*t_phi
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1= floor(abs(z1)/zmax(ri1)*nz) + 1
        if (zj1>nz) zj1=nz+1
        if (z1 < 0.0) zj1=-zj1
        phik1=phik0+delta_phi
        if (phik1 == 0) phik1=N_az
        if (phik1 == N_az+1) phik1=1
     endif

     ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
     if (z1 == 0.0_db) z1 = sign(prec_grille,w)

     tau=l*opacite ! opacite constante dans la cellule

     tau_tot = tau_tot + tau
     ltot= ltot + l
     if (tau_tot < tiny_real) lmin=ltot

  enddo ! boucle infinie
  write(*,*) "BUG"
  return

end subroutine length_deg2_tot_3D

!********************************************************

subroutine length_deg2_tot_sph(id,lambda,Stokes,ri,thetaj,xi,yi,zi,u,v,w,tau_tot_out,lmin,lmax)
! Integration par calcul de la position de l'interface entre cellules
! de l'opacite totale dans une direction donnée
! Grille a geometrie spherique
! C. Pinte
! 19/06/07

  implicit none

  integer, intent(in) :: id, lambda, ri, thetaj
  real(kind=db), dimension(4), intent(in) :: Stokes
  real(kind=db), intent(in) :: xi, yi, zi, u,v,w
  real, intent(out) :: tau_tot_out
  real(kind=db), intent(out) :: lmin,lmax

  real(kind=db) :: x0, y0, z0, x1, y1, z1, ltot, tau_tot
  real(kind=db) :: b, c, s, rac, t, delta, r0_2, r0_cyl, r0_2_cyl
  real(kind=db) :: delta_vol, l, tau, dotprod, opacite
  real(kind=db) :: correct_moins, correct_plus, uv
  integer :: ri0, thetaj0, ri1, thetaj1, delta_rad, delta_theta, nbr_cell, icell0

  logical :: lcellule_non_vide

  real(kind=db) :: a_theta, b_theta, c_theta, tan2, tan_angle_lim1, tan_angle_lim2, t1, t2
  real(kind=db) :: t1_1, t1_2, t2_1, t2_2

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - 1.0e-10_db
  correct_plus = 1.0_db + 1.0e-10_db

  x1=xi;y1=yi;z1=zi
  ri1=ri
  thetaj1=thetaj

  tau_tot=0.0_db
  lmin=0.0_db
  ltot=0.0_db
  nbr_cell = 0

  uv = sqrt(u*u + v*v)

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1;thetaj0=thetaj1
     x0=x1;y0=y1;z0=z1

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        tau_tot_out=tau_tot
        lmax=ltot
        return
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r0_2_cyl  = x0*x0+y0*y0
     r0_cyl = sqrt(r0_2_cyl)
     r0_2 = r0_2_cyl + z0*z0
     b   = (x0*u+y0*v+z0*w)

     if (ri0==0) then
        lcellule_non_vide=.false.
        opacite=0.0_db
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r0_2-r_lim_2(0)*correct_plus)
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        delta_rad=1
     else
        icell0 = cell_map(ri0,thetaj0,1)
        if (icell0 > n_cells) then
           opacite = 0.0_db
        else
           opacite=kappa(icell0,lambda)
        endif

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod= b
        if (dotprod < 0.0_db) then
           ! on recule : on cherche rayon inférieur
           c=(r0_2-r_lim_2(ri0-1)*correct_moins)
           delta=b*b-c
           if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r0_2-r_lim_2(ri0)*correct_plus)
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r0_2-r_lim_2(ri0)*correct_plus)
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0_db) then
           s=-b+rac
        else if (s==0.0_db) then
           s=prec_grille
        endif


        ! 2) position interface inclinaison
        ! Meme methode que azimuth dans version cylindrique 3D
        if (z0 >= 0.0_db) then ! on est dans le bon sens
           tan_angle_lim1 = tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = tan_theta_lim(thetaj0-1) * correct_moins
        else
           tan_angle_lim1 = - tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = - tan_theta_lim(thetaj0-1) * correct_moins
        endif ! z0

        ! Premiere limite theta
        tan2 = tan_angle_lim1 * tan_angle_lim1
        a_theta = w*w - tan2 * (u*u + v*v)
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t1 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t1_1 = (- b_theta - rac)/a_theta
           t1_2 = (- b_theta + rac)/a_theta

           if (t1_1 <= 0.0_db) then
              if (t1_2 <= 0.0_db) then ! les 2 sont <0
                 t1=1.0e30_db
              else   ! Seul t1_2 > 0
                 t1 = t1_2
              endif
           else
              if (t1_2 <= 0.0_db) then ! Seul t1_1 >0
                 t1 = t1_1
              else ! les 2 sont > 0
                 t1 = min(t1_1,t1_2)
              endif
           endif
        endif ! signe delta

        ! Deuxieme limite theta
        tan2 = tan_angle_lim2 * tan_angle_lim2
        a_theta = w*w - tan2 * (u*u + v*v)
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t2 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t2_1 = (- b_theta - rac)/a_theta
           t2_2 = (- b_theta + rac)/a_theta

           if (t2_1 <= 0.0_db) then
              if (t2_2 <= 0.0_db) then ! les 2 sont <0
                 t2=1.0e30_db
              else   ! Seul t2_2 > 0
                 t2 = t2_2
              endif
           else
              if (t2_2 <= 0.0_db) then ! Seul t2_1 >0
                 t2 = t2_1
              else ! les 2 sont > 0
                 t2 = min(t2_1,t2_2)
              endif
           endif
         endif ! signe delta

        ! Selection limite theta
        if (t1 < t2) then
           t=t1
           delta_theta = 1
           if (thetaj0 == nz) delta_theta = 0
        else
           t=t2
           delta_theta = -1
           if (thetaj0 == 1) delta_theta = 0
        endif

     endif ! ri0==0

     ! 3) interface en r ou theta ?
     if (s < t) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        thetaj1 = thetaj0

        ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
        if (ri0==0) then
           call indice_cellule_sph_theta(x1,y1,z1,thetaj1)
        endif

     else ! theta
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        thetaj1=thetaj0+delta_theta
     endif

     ! Calcul longeur de vol et profondeur optique dans la cellule
     tau=l*opacite ! opacite constante dans la cellule

     tau_tot = tau_tot + tau
     ltot= ltot + l
     if (tau_tot < tiny_real) lmin=ltot

  enddo ! boucle infinie

  write(*,*) "BUG"
  return

end subroutine length_deg2_tot_sph

!***********************************************************

subroutine integ_ray_mol(id,ri,zj,phik,x,y,z,u,v,w,iray,labs,ispeed,tab_speed)
  ! Generalisation de la routine length_deg2
  ! pour le cas du transfert dans les raies
  ! Propage un paquet depuis un point d'origine donne
  !
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: id, ri,zj, phik, iray
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: x,y,z
  logical, intent(in) :: labs
  integer, dimension(2), intent(in) :: ispeed
  real(kind=db), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  if (lcylindrical) then
     call integ_ray_mol_cyl(id,ri,zj,phik,x,y,z,u,v,w,iray,labs,ispeed,tab_speed)
  else
     call integ_ray_mol_sph(id,ri,zj,phik,x,y,z,u,v,w,iray,labs,ispeed,tab_speed)
  endif

  return

end subroutine integ_ray_mol

!***********************************************************

subroutine integ_ray_mol_cyl(id,ri_in,zj_in,phik_in,x,y,z,u,v,w,iray,labs,ispeed,tab_speed)
  ! Generalisation de la routine length_deg2
  ! pour le cas du transfert dans les raies
  ! Propage un paquet depuis un point d'origine donne
  ! et integre l'equation du transfert radiatif
  ! La propagation doit etre a l'envers pour faire du
  ! ray tracing  !!
  !
  ! C. Pinte
  ! 12/04/07

  implicit none

  integer, intent(in) :: id, ri_in, zj_in, phik_in, iray
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: x,y,z
  logical, intent(in) :: labs
  integer, dimension(2), intent(in) :: ispeed
  real(kind=db), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  real(kind=db), dimension(ispeed(1):ispeed(2)) :: tspeed
  real(kind=db) :: x0, y0, z0, x1, y1, z1, xphi, yphi, zphi
  real(kind=db) :: inv_a, a, b, c, s, rac, t, t_phi, delta, inv_w, r_2, tan_angle_lim, den
  real(kind=db) :: delta_vol, l, zlim, dotprod, delta_vol_phi
  real(kind=db), dimension(ispeed(1):ispeed(2)) :: P, dtau, dtau2, Snu, opacite
  real(kind=db), dimension(ispeed(1):ispeed(2),nTrans) :: tau, tau2
  real(kind=db) :: dtau_c, Snu_c
  real(kind=db), dimension(nTrans) :: tau_c
  real(kind=db) :: correct_plus, correct_moins, v0, v1, v_avg0
  integer :: ri0, zj0, ri1, zj1, phik0, phik1, delta_rad, delta_zj, nbr_cell, delta_phi, phik0m1, icell0
  integer :: iTrans, ivpoint, iiTrans, n_vpoints

  real :: facteur_tau

  logical :: lcellule_non_vide

  integer, parameter :: n_vpoints_max = 200 ! pas super critique, presque OK avec 2 pour la simu Herbig de Peter (2x plus vite)
  real(kind=db), dimension(n_vpoints_max) :: vitesse

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  ri1=ri_in
  zj1=zj_in
  phik1=phik_in

  nbr_cell = 0

  a=u*u+v*v

  if (a > tiny_real) then
     inv_a=1.0/a
  else
     inv_a=huge_real
  endif


  if (abs(w) > tiny_real) then
     inv_w=1.0/w
  else
     inv_w=sign(huge_db,w)
  endif

  tau(:,:) = 0.0_db
  I0(:,:,iray,id) = 0.0_db
  v_avg0 = 0.0_db

  tau_c(:) = 0.0_db
  I0c(:,iray,id) = 0.0_db

  !*** propagation dans la grille

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule

     ri0=ri1 ; zj0=zj1 ; phik0=phik1
     x0=x1 ; y0=y1 ; z0=z1

     icell0 = cell_map(ri0,zj0,phik0)

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        exit infinie
     elseif (abs(zj0)>nz) then
        lcellule_non_vide=.false.
        ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           exit infinie
        endif
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r_2=x0*x0+y0*y0
     b=(x0*u+y0*v)*inv_a

     if (ri0==0) then
        lcellule_non_vide=.false.
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r_2-r_lim_2(0)*correct_plus)*inv_a
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        t_phi= huge_real
        delta_rad=1
     else

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod=u*x0+v*y0
        if (dotprod < 0.0) then
           ! on recule : on cherche rayon inférieur
           c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
           delta=b*b-c
           if (delta < 0.0) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0) then
           s=-b+rac
        else if (s==0.0) then
           s=prec_grille
        endif

        ! 2) position interface verticale
        ! on monte ou on descend par plan équatorial ?
        dotprod=w*z0
        if (dotprod == 0.0_db) then
           t=1.0e10
        else
           if (dotprod > 0.0_db) then ! on se rapproche de la surface
              if (l3D) then
                 if (zj0==nz+1) then
                    delta_zj=0
                    zlim=1.0e10
                 else if (zj0==-(nz+1)) then
                    delta_zj=0
                    zlim=-1.0e10
                 else
                    if (z0 > 0.0) then
                       zlim=z_lim(ri0,zj0+1)*correct_plus
                       delta_zj=1
                    else
                       zlim=-z_lim(ri0,abs(zj0)+1)*correct_plus
                       delta_zj=-1
                    endif
                 endif
              else ! 2D
                 if (zj0==nz+1) then
                    delta_zj=0
                    if (z0 > 0.0_db) then
                       zlim=1.0e10
                    else
                       zlim=-1.0e10
                    endif
                 else
                    if (z0 > 0.0) then
                       zlim=z_lim(ri0,zj0+1)*correct_plus
                    else
                       zlim=-z_lim(ri0,zj0+1)*correct_plus
                    endif
                    delta_zj=1
                 endif
              endif ! l3D
           else ! on se rappoche du midplane
              if (l3D) then
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,abs(zj0))*correct_moins
                    delta_zj=-1
                    if (zj0==1) delta_zj=-2 ! pas d'indice 0
                 else
                    zlim=-z_lim(ri0,abs(zj0))*correct_moins
                    delta_zj=1
                    if (zj0==-1) delta_zj=2 ! pas d'indice 0
                 endif
              else ! 2D
                 if (zj0==1) then
                    ! on traverse le plan eq donc on va remonter
                    ! et z va changer de signe
                    delta_zj=1
                    if (z0 > 0.0_db) then
                       zlim=-z_lim(ri0,2)*correct_moins
                    else
                       zlim=z_lim(ri0,2)*correct_moins
                    endif
                 else !(zj0==1)
                    ! on ne traverse pas z=0.
                    if (z0 > 0.0_db) then
                       zlim=z_lim(ri0,zj0)*correct_moins
                    else
                       zlim=-z_lim(ri0,zj0)*correct_moins
                    endif
                    delta_zj=-1
                 endif !(zj0==1)
              endif ! l3D
           endif ! monte ou descend
           t=(zlim-z0)*inv_w
           ! correct pb precision
           if (t < 0.0_db) t=prec_grille
        endif !dotprod=0.0

        ! 3) position interface azimuthale
        dotprod =  x0*v - y0*u


        if (abs(dotprod) < 1.0e-10) then
           ! on ne franchit pas d'interface azimuthale
           t_phi = 1.0e30
        else
           ! Quelle cellule on va franchir
           if (dotprod > 0.0) then
              tan_angle_lim = tan_phi_lim(phik0)
              delta_phi=1
           else
              phik0m1=phik0-1
              if (phik0m1==0) phik0m1=N_az
              tan_angle_lim = tan_phi_lim(phik0m1)
              delta_phi=-1
           endif
           ! Longueur av interserction
           if (tan_angle_lim > 1.0d299) then
              t_phi = -x0/u
           else
              den= v-u*tan_angle_lim
              if (abs(den) > 1.0e-6) then
                 t_phi = -(y0-x0*tan_angle_lim)/den
              else
                 t_phi = 1.0e30
              endif
           endif
           if (t_phi < 0.0) t_phi = 1.0e30
        endif !dotprod = 0.0

     endif ! ri==0

     ! 4) interface en r ou z ou phi ?
     if ((s < t).and.(s < t_phi)) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        if ((ri1<1).or.(ri1>n_rad)) then
           zj1=zj0
        else
           zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
           if (zj1>nz) zj1=nz+1
           if (l3D) then
              if (z1 < 0.0) zj1=-zj1
           endif
        endif
        phik1=phik0

        ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
        if (ri0==0) call indice_cellule_3D_phi(x1,y1,z1,phik1)
     else if (t < t_phi) then ! z
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1=zj0+delta_zj
        phik1=phik0
     else
        l=t_phi
        delta_vol=correct_plus*t_phi
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1= floor(abs(z1)/zmax(ri1)*nz) + 1
        if (zj1>nz) zj1=nz+1
        if (l3D) then
           if (z1 < 0.0) zj1=-zj1
        endif
        phik1=phik0+delta_phi
        if (phik1 == 0) phik1=N_az
        if (phik1 == N_az+1) phik1=1
     endif

     ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
     if (z1 == 0.0_db) then
        if (l3D) then
           z1 = sign(prec_grille,w)
        else ! 2D
           z1 = prec_grille
        endif ! l3D
     endif

     if (lcellule_non_vide) then
        ! Differentiel de vitesse au travers de la cellule
        !dv = dv_proj(ri0,zj0,x0,y0,z0,x1,y1,z1,u,v,w)
        v0 = v_proj(ri0,zj0,x0,y0,z0,u,v,w)
        v1 = v_proj(ri0,zj0,x1,y1,z1,u,v,w)
        dv = abs(v1 - v0)

        ! Nbre de points d'integration en fct du differentiel de vitesse
        ! compare a la largeur de raie de la cellule de depart
        n_vpoints  = min(max(2,nint(dv/v_line(icell0)*20.)),n_vpoints_max)

        ! Vitesse projete le long du trajet dans la cellule
        do ivpoint=2, n_vpoints-1
           delta_vol_phi = (real(ivpoint,kind=db))/(real(n_vpoints,kind=db)) * delta_vol
           xphi=x0+delta_vol_phi*u
           yphi=y0+delta_vol_phi*v
           zphi=z0+delta_vol_phi*w
           vitesse(ivpoint) = v_proj(ri0,zj0,xphi,yphi,zphi,u,v,w)
        enddo
        vitesse(1) = v0
        vitesse(n_vpoints) = v1

        if ((nbr_cell == 1).and.labs) then
           v_avg0 = 0.0_db
           do ivpoint=1,n_vpoints
              v_avg0 = v_avg0 + vitesse(ivpoint)
           enddo
           v_avg0 = v_avg0 / real(n_vpoints,kind=db)
        endif

        ! Profil de raie local integre a multiplie par la frequence de la transition
        P(:) = 0.0_db
        do ivpoint=1,n_vpoints
           tspeed(:) = tab_speed(:) - (vitesse(ivpoint) - v_avg0)
           P(:) = P(:) + phiProf(icell0,ispeed,tspeed)
        enddo
        P(:) = P(:)/n_vpoints

        if ((nbr_cell == 1).and.labs) then
           ds(iray,id) = delta_vol
            Doppler_P_x_freq(:,iray,id) = P(:)
        endif

        do iTrans=1,nTrans
           iiTrans = indice_Trans(iTrans)

           opacite(:) = kappa_mol_o_freq(icell0,iiTrans) * P(:) + kappa(icell0,iiTrans)

           ! Epaisseur optique
           dtau(:) =  l * opacite(:)
           dtau_c = l * kappa(icell0,iiTrans)

           ! Fonction source
           Snu(:) = ( emissivite_mol_o_freq(icell0,iiTrans) * P(:) &
                + emissivite_dust(icell0,iiTrans) ) / (opacite(:) + 1.0e-300_db)
           Snu_c = emissivite_dust(icell0,iiTrans) / (kappa(icell0,iiTrans) + 1.0e-300_db)

           ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
           ! la profondeur optique jusqu'a la cellule
           !---write(*,*) ri0, zj0
           !---write(*,*) "kappa", kappa_mol_o_freq(ri0,zj0,iiTrans), kappa(iiTrans,ri0,zj0,1)
           !---write(*,*) "eps", emissivite_mol_o_freq(ri0,zj0,iiTrans)
           !---
           !---write(*,*) minval(tau(:,iTrans)), maxval(tau(:,iTrans))
           !---write(*,*) minval(dtau(:)), maxval(dtau(:))
           !---write(*,*) minval(Snu(:)), maxval(Snu(:))
           I0(:,iTrans,iray,id) = I0(:,iTrans,iray,id) + &
                exp(-tau(:,iTrans)) * (1.0_db - exp(-dtau(:))) * Snu(:)
           I0c(iTrans,iray,id) = I0c(iTrans,iray,id) + &
                exp(-tau_c(iTrans)) * (1.0_db - exp(-dtau_c)) * Snu_c

           if (lorigine.and.(.not.labs)) then
              origine_mol(:,iiTrans,ri0,zj0,id) = origine_mol(:,iiTrans,ri0,zj0,id) + &
                   exp(-tau(:,iTrans)) * (1.0_db - exp(-dtau(:))) * Snu(:)
           endif

           ! surface superieure ou inf
           facteur_tau = 1.0
           if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
           if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

           ! Mise a jour profondeur optique pour cellule suivante
           tau(:,iTrans) = tau(:,iTrans) + dtau(:) * facteur_tau
           tau_c(iTrans) = tau_c(iTrans) + dtau_c
        enddo

        if (ldouble_RT) then
           do iTrans=1,nTrans
              iiTrans = indice_Trans(iTrans)
              opacite(:) = kappa_mol_o_freq2(icell0,iiTrans) * P(:) + kappa(icell0,iiTrans)
              dtau(:) =  l * opacite(:)

              ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
              ! la profondeur optique jusqu'a la cellule
              Snu(:) = ( emissivite_mol_o_freq2(icell0,iiTrans) * P(:) + &
                   emissivite_dust(icell0,iiTrans) ) / (opacite(:) + 1.0e-30_db)
              I02(:,iTrans,iray,id) = I02(:,iTrans,iray,id) + &
                   exp(-tau2(:,iTrans)) * (1.0_db - exp(-dtau2(:))) * Snu(:)

              ! Mise a jour profondeur optique pour cellule suivante
              tau2(:,iTrans) = tau2(:,iTrans) + dtau2(:)
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

  do iTrans=1,nTrans
     iiTrans = indice_Trans(iTrans)
     !I0(:,iTrans,iray,id) = I0(:,iTrans,iray,id) + tab_Cmb_mol(iiTrans) * exp(-tau(:,iTrans))
  enddo

  if (ldouble_RT) then
     do iTrans=1,nTrans
        iiTrans = indice_Trans(iTrans)
        I02(:,iTrans,iray,id) = I02(:,iTrans,iray,id) + tab_Cmb_mol(iiTrans) * exp(-tau2(:,iTrans))
     enddo
  endif

  return

end subroutine integ_ray_mol_cyl

!***********************************************************

subroutine integ_ray_mol_sph(id,ri_in,thetaj_in,phik_in,x,y,z,u,v,w,iray,labs,ispeed,tab_speed)
  ! Generalisation de la routine length_deg2
  ! pour le cas du transfert dans les raies
  ! Propage un paquet depuis un point d'origine donne
  ! et integre l'equation du transfert radiatif
  ! La propagation doit etre a l'envers pour faire du
  ! ray tracing  !!
  !
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: id, ri_in,thetaj_in, phik_in, iray
  real(kind=db), intent(in) :: u,v,w, x,y,z
  logical, intent(in) :: labs
  integer, dimension(2), intent(in) :: ispeed
  real(kind=db), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  real(kind=db), dimension(ispeed(1):ispeed(2)) :: tspeed
  real(kind=db) :: x0, y0, z0, x1, y1, z1, xphi, yphi, zphi
  real(kind=db) :: b, c, s, rac, t, t_phi, delta, tan_angle_lim, den
  real(kind=db) :: delta_vol, l, dotprod, delta_vol_phi
  real(kind=db),  dimension(ispeed(1):ispeed(2)) :: P
  real(kind=db), dimension(ispeed(1):ispeed(2),nTrans) :: tau, tau2
  real(kind=db), dimension(ispeed(1):ispeed(2)) :: dtau, dtau2, opacite, Snu
  real(kind=db) :: correct_plus, correct_moins, v0, v1, v_avg0, precision
  integer :: ri0, thetaj0, ri1, thetaj1, phik0, phik1, delta_rad, delta_theta, delta_phi, phik0m1, icell0
  integer :: iTrans, nbr_cell, ivpoint, iiTrans, n_vpoints

  real :: facteur_tau

  logical :: lcellule_non_vide

  real(kind=db) :: tan_angle_lim1, tan_angle_lim2, a_theta, b_theta, c_theta
  real(kind=db) :: t1, t1_1, t1_2, t2, t2_1, t2_2, uv, r0_cyl, r0_2_cyl, r0_2, tan2
  integer, parameter :: n_vpoints_max = 200
  real(kind=db), dimension(n_vpoints_max) :: vitesse

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - 1.0e-10_db
  correct_plus = 1.0_db + 1.0e-10_db
  precision = 1.0e-20_db ! pour g95

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  ri1=ri_in
  thetaj1=thetaj_in
  phik1=phik_in

  uv = sqrt(u*u + v*v)

  nbr_cell = 0

  tau(:,:) = 0.0_db
  I0(:,:,iray,id) = 0.0_db
  v_avg0 = 0.0_db

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1 ; thetaj0=thetaj1 ; phik0=phik1
     x0=x1 ; y0=y1 ; z0=z1

     icell0 = cell_map(ri0,thetaj0,phik0)

     nbr_cell = nbr_cell + 1

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        exit infinie
     endif ! Test sortie


     ! Detection interface
     r0_2_cyl  = x0*x0+y0*y0
     r0_cyl = sqrt(r0_2_cyl)
     r0_2 = r0_2_cyl + z0*z0
     b   = (x0*u+y0*v+z0*w)

     if (ri0==0) then
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        lcellule_non_vide = .false.
        c=(r0_2-r_lim_2(0)*correct_plus)
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        delta_rad=1
     else

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod= b
        if (dotprod < 0.0_db) then
           ! on recule : on cherche rayon inférieur
           c=(r0_2-r_lim_2(ri0-1)*correct_moins)
           delta=b*b-c
           if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r0_2-r_lim_2(ri0)*correct_plus)
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r0_2-r_lim_2(ri0)*correct_plus)
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0_db) then
           s=-b+rac
        else if (s==0.0_db) then
           s=prec_grille
        endif


        ! 2) position interface inclinaison
        ! Meme methode que azimuth dans version cylindrique 3D
        if (z0 >= 0.0_db) then ! on est dans le bon sens
           tan_angle_lim1 = tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = tan_theta_lim(thetaj0-1) * correct_moins
        else
           tan_angle_lim1 = - tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = - tan_theta_lim(thetaj0-1) * correct_moins
        endif ! z0

        ! Premiere limite theta
        tan2 = tan_angle_lim1 * tan_angle_lim1
        a_theta = w*w - tan2 * (u*u + v*v)
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t1 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t1_1 = (- b_theta - rac)/a_theta
           t1_2 = (- b_theta + rac)/a_theta

           if (t1_1 <= precision) then
              if (t1_2 <= precision) then ! les 2 sont <0
                 t1=1.0e30_db
              else   ! Seul t1_2 > 0
                 t1 = t1_2
              endif
           else
              if (t1_2 <= precision) then ! Seul t1_1 >0
                 t1 = t1_1
              else ! les 2 sont > 0
                 t1 = min(t1_1,t1_2)
              endif
           endif
        endif ! signe delta

        ! Deuxieme limite theta
        tan2 = tan_angle_lim2 * tan_angle_lim2
        a_theta = w*w - tan2 * (u*u + v*v)
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t2 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t2_1 = (- b_theta - rac)/a_theta
           t2_2 = (- b_theta + rac)/a_theta

           if (t2_1 <= precision) then
              if (t2_2 <= precision) then ! les 2 sont <0
                 t2=1.0e30_db
              else   ! Seul t2_2 > 0
                 t2 = t2_2
              endif
           else
              if (t2_2 <= 0.0_db) then ! Seul t2_1 >0
                 t2 = t2_1
              else ! les 2 sont > 0
                 t2 = min(t2_1,t2_2)
              endif
           endif
        endif ! signe delta

        ! Selection limite theta
        if (t1 < t2) then
           t=t1
           delta_theta = 1
           if (thetaj0 == nz) delta_theta = 0
        else
           t=t2
           delta_theta = -1
           if (thetaj0 == 1) delta_theta = 0
        endif


        ! 3) position interface azimuthale
        dotprod =  x0*v - y0*u
        if (abs(dotprod) < 1.0e-10) then
           ! on ne franchit pas d'interface azimuthale
           t_phi = 1.0e30
           delta_phi = 0
        else
           ! Quelle cellule on va franchir
           if (dotprod > 0.0) then
              tan_angle_lim = tan_phi_lim(phik0)
              delta_phi=1
           else
              phik0m1=phik0-1
              if (phik0m1==0) phik0m1=N_az
              tan_angle_lim = tan_phi_lim(phik0m1)
              delta_phi=-1
           endif
           ! Longueur av interserction
           if (tan_angle_lim > 1.0d299) then
              t_phi = -x0/u
              delta_phi=0
           else
              den= v-u*tan_angle_lim
              if (abs(den) > 1.0e-6) then
                 t_phi = -(y0-x0*tan_angle_lim)/den
              else
                 t_phi = 1.0e30
                 delta_phi = 0
              endif
           endif
           if (t_phi < 0.0) then
              t_phi = 1.0e30
              delta_phi = 0
           endif
        endif !dotprod = 0.0

     endif ! ri0==0


     ! 4) interface en r ou z ou phi ?
     if ((s < t).and.(s < t_phi)) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        thetaj1 = thetaj0
        phik1=phik0

        ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
        if (ri0==0) then
           call indice_cellule_sph_theta(x1,y1,z1,thetaj1)
        endif
     else if (t < t_phi) then ! theta
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        thetaj1=thetaj0+delta_theta
        phik1=phik0
     else
        l=t_phi
        delta_vol=correct_plus*t_phi
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        thetaj1 = thetaj0
        phik1=phik0+delta_phi
        if (phik1 == 0) phik1=N_az
        if (phik1 == N_az+1) phik1=1
     endif


     if (lcellule_non_vide) then
        ! Differentiel de vitesse au travers de la cellule
        !dv = dv_proj(ri0,zj0,x0,y0,z0,x1,y1,z1,u,v,w)

        v0 = v_proj(ri0,thetaj0,x0,y0,z0,u,v,w)
        v1 = v_proj(ri0,thetaj0,x1,y1,z1,u,v,w)
        dv = abs(v1 - v0)

        ! Nbre de points d'integration en fct du differentiel de vitesse
        ! compare a la largeur de raie de la cellule de depart
        n_vpoints  = min(max(2,nint(dv/v_line(icell0)*20.)),n_vpoints_max)
     !   n_vpoints = 2

!        write(*,*) ri0, real(dv), real(v0), real(v1), real(v_line(ri0,thetaj0)), nint(dv/v_line(ri0,thetaj0)*20.)
!        read(*,*)


        ! Vitesse projete le long du trajet dans la cellule
        do ivpoint=2, n_vpoints-1
           delta_vol_phi = (real(ivpoint,kind=db))/(real(n_vpoints,kind=db)) * delta_vol
           xphi=x0+delta_vol_phi*u
           yphi=y0+delta_vol_phi*v
           zphi=z0+delta_vol_phi*w
           vitesse(ivpoint) = v_proj(ri0,thetaj0,xphi,yphi,zphi,u,v,w)
        enddo
        vitesse(1) = v0
        vitesse(n_vpoints) = v1

     !   write(*,*) ri0
     !   write(*,*) vitesse(1:n_vpoints)
     !   read(*,*)


        if ((nbr_cell == 1).and.labs) then
           v_avg0 = 0.0_db
           do ivpoint=1,n_vpoints
              v_avg0 = v_avg0 + vitesse(ivpoint)
           enddo
           v_avg0 = v_avg0 / real(n_vpoints,kind=db)
        endif
       ! v_avg0 = 0.0_db


        ! Profil de raie local integre multiplie par la frequence de la transition
        P(:) = 0.0_db
        do ivpoint=1,n_vpoints
           tspeed(:) = tab_speed(:) - (vitesse(ivpoint) - v_avg0)
           P(:) = P(:) + phiProf(icell0,ispeed,tspeed)
        enddo
        P(:) = P(:)/n_vpoints

        if ((nbr_cell == 1).and.labs) then
           ds(iray,id) = delta_vol
            Doppler_P_x_freq(:,iray,id) = P(:)
        endif

        do iTrans=1,nTrans
           iiTrans = indice_Trans(iTrans)
           opacite(:) = kappa_mol_o_freq(icell0,iiTrans) * P(:) + kappa(icell0,iiTrans)
           dtau(:) =  l * opacite(:)

           ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
           ! la profondeur optique jusqu'a la cellule
           Snu(:) = ( emissivite_mol_o_freq(icell0,iiTrans) * P(:) + &
                emissivite_dust(icell0,iiTrans) ) / (opacite(:) + 1.0e-30_db)
           I0(:,iTrans,iray,id) = I0(:,iTrans,iray,id) + exp(-tau(:,iTrans)) * (1.0_db - exp(-dtau(:))) * Snu(:)

           ! surface superieure ou inf
           facteur_tau = 1.0
           if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
           if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

           ! Mise a jour profondeur optique pour cellule suivante
           tau(:,iTrans) = tau(:,iTrans) + dtau(:) * facteur_tau
        enddo

        if (ldouble_RT) then
           do iTrans=1,nTrans
              iiTrans = indice_Trans(iTrans)
              opacite(:) = kappa_mol_o_freq2(icell0,iiTrans) * P(:) + kappa(icell0,iiTrans)
              dtau(:) =  l * opacite(:)

              ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
              ! la profondeur optique jusqu'a la cellule
              Snu(:) = ( emissivite_mol_o_freq2(icell0,iiTrans) * P(:) + &
                   emissivite_dust(icell0,iiTrans) ) / (opacite(:) + 1.0e-30_db)
              I0(:,iTrans,iray,id) = I0(:,iTrans,iray,id) + exp(-tau2(:,iTrans)) * (1.0_db - exp(-dtau2(:))) * Snu(:)

              ! Mise a jour profondeur optique pour cellule suivante
              tau2(:,iTrans) = tau2(:,iTrans) + dtau2(:)
           enddo
        endif


     endif

  enddo infinie

  ! Ajout Cmb, pondere par la profondeur optique totale
!  tspeed(:) = tab_speed(:) + v_avg0
!  I0(:,:,iray,id) = I0(:,:,iray,id) + Cmb(ispeed,tspeed) * exp(-tau(:,:))
!  ! pas sur que v_avg0 soit ds le bon sens mais c'est necligeable par rapport a c
!  if (ldouble_RT) then
!     I02(:,:,iray,id) = I02(:,:,iray,id) + Cmb(ispeed,tspeed) * exp(-tau2(:,:))
!  endif

  do iTrans=1,nTrans
     iiTrans = indice_Trans(iTrans)
     I0(:,iTrans,iray,id) = I0(:,iTrans,iray,id) + tab_Cmb_mol(iiTrans) * exp(-tau(:,iTrans))
  enddo

  if (ldouble_RT) then
     do iTrans=1,nTrans
        iiTrans = indice_Trans(iTrans)
        I02(:,iTrans,iray,id) = I02(:,iTrans,iray,id) + tab_Cmb_mol(iiTrans) * exp(-tau2(:,iTrans))
     enddo
  endif

  return

end subroutine integ_ray_mol_sph

!***********************************************************

subroutine move_to_grid(x,y,z,u,v,w,ri,zj,phik,lintersect)
  ! Calcule la position au bord de la grille dans
  ! la direction donnee
  ! C. Pinte
  ! 01/08/07

  implicit none

  real(kind=db), intent(inout) :: x,y,z
  real(kind=db), intent(in) :: u,v,w
  integer, intent(out) :: ri, zj, phik
  logical, intent(out) :: lintersect

  if (lcylindrical) then
     call move_to_grid_cyl(x,y,z,u,v,w,ri,zj,phik,lintersect)
  else
     call move_to_grid_sph(x,y,z,u,v,w,ri,zj,lintersect)
     phik=1 ! pas 3D
  endif

  if (zj == nz) then ! indice_cellule force nz si z > zmax, il faut corriger ca
     if (ri < n_rad) zj = nz+1 ; ! suppose que zmaxmax est atteint dans la derniere cellule
  endif

  return

end subroutine move_to_grid

!***********************************************************

subroutine move_to_grid_cyl(x,y,z,u,v,w,ri,zj,phik,lintersect)
  ! Calcule la position au bord de la grille dans
  ! la direction donnee pour grille cylindrique
  ! C. Pinte
  ! 19/09/07

  implicit none

  real(kind=db), intent(inout) :: x,y,z
  real(kind=db), intent(in) :: u,v,w
  integer, intent(out) :: ri, zj, phik
  logical, intent(out) :: lintersect

  real(kind=db) :: x0, y0, z0, z1, a, inv_a, r_2, b, c, delta, rac, s1, s2, dotprod, t1, t2
  real(kind=db) :: zlim, zlim2, delta_vol, inv_w, correct_moins, correct_plus

  correct_moins = 1.0_db - 1.0e-10_db
  correct_plus = 1.0_db + 1.0e-10_db

  x0=x ; y0=y ; z0=z

  a=u*u+v*v
  if (a > tiny_real) then
     inv_a=1.0_db/a
  else
     inv_a=huge_real
  endif

  if (abs(w) > tiny_real) then
     inv_w=1.0_db/w
  else
     inv_w=sign(huge_db,w) ! huge_real avant
  endif

  ! Longueur de vol pour atteindre le rayon cylindrique rout
  r_2=x0*x0+y0*y0
  b=(x0*u+y0*v)*inv_a

  c=(r_2-r_lim_2(n_rad)*correct_moins)*inv_a
  delta=b*b-c
  if (delta < 0.0_db) then
     ! On ne rencontre pas le cylindre
     s1 = huge_real
     s2 = huge_real
  else
     ! On rencontre le cylindre
     rac=sqrt(delta)

     ! Les deux racines doivent etre positives sinon BUG !!
     ! s1 < s2
     s1=-b-rac
     s2=-b+rac

     ! TMP : BUG : ca plante si rayon vertical !!!
   !  if (s1 < 0.0) then
   !     write(*,*) "Bug dans ray tracing !!!", s1, s2
   !  endif
     ! END TMP
  endif

  ! longueur de vol pour atteindre zmax
  ! le rayon monte ou descend ?
  dotprod=w*z0

  if (abs(dotprod) < tiny_real) then ! rayon horizontal
     t1=huge_real
     t2=huge_real
  else
     if (z0 > 0.0_db) then
        zlim=zmaxmax*correct_moins
        zlim2=-zmaxmax*correct_moins
     else
        zlim=-zmaxmax*correct_moins
        zlim2=zmaxmax*correct_moins
     endif
     t1=(zlim-z0)*inv_w
     t2=(zlim2-z0)*inv_w
  endif !dotprod=0.0

  ! On ne rencontre ni le cylindre ni les plans
  if (t1 > 1e20) then
     if (s1 > 1e20) then
        lintersect = .false.
        return
     endif
  endif

  ! On rentre ?? et si oui, par le dessus ou par rout ??
  if (t1 > s1) then ! On rentre d'abord dans le cylindre
     if (t1 > s2) then ! on ressort du cylindre avant de croiser la tranche
        ! On se place au bord du cylindre
        delta_vol = s1
        z1 = z0 + delta_vol * w
        if (abs(z1) > zmaxmax) then ! On est constamment en dehors des 2 plans
           lintersect=.false.
           return
        else ! on est constamment entre les 2 plans
           lintersect = .true.
        endif
     else  ! on va croiser la surface dans le cylindre
        lintersect = .true.
        delta_vol = t1
     endif

  else  ! on croise d'abord le plan
     if (t2 < s1) then
        ! on ressort de la tranche avant de croiser le cylindre
        lintersect=.false.
        return
     else
        lintersect = .true.
        delta_vol = s1
     endif
  endif

  ! Position au bord de la grille
  x=x0+delta_vol*u!*correct_plus
  y=y0+delta_vol*v!*correct_plus
  z=z0+delta_vol*w!*correct_plus

  ! Determination de l'indice de la premiere cellule traversee
  ! pour initialiser la propagation
  if (l3D) then
     call indice_cellule_3D(x,y,z,ri,zj,phik)
  else
     call indice_cellule(x,y,z,ri,zj) ; phik=1
  endif

  return

end subroutine move_to_grid_cyl

!***********************************************************

subroutine move_to_grid_sph(x,y,z,u,v,w,ri,thetaj,lintersect)
  ! Calcule la position au bord de la grille dans
  ! la direction donnee pour grille spherique
  ! C. Pinte
  ! 01/08/07

  implicit none

  real(kind=db), intent(inout) :: x,y,z
  real(kind=db), intent(in) :: u,v,w
  integer, intent(out) :: ri, thetaj
  logical, intent(out) :: lintersect

  real(kind=db) :: x0, y0, z0, x1, y1, z1, r0_2, b, c, rac, delta, s1, delta_vol, correct_moins

  correct_moins = 1.0_db - 1.0e-10_db

  x0=x ; y0=y ; z0 =z


  r0_2 = x0*x0+y0*y0+z0*z0
  b   = (x0*u+y0*v+z0*w)
  c=(r0_2-r_lim_2(n_rad)*correct_moins)
  delta=b*b-c

  if (delta < 0.0_db) then
     ! On ne rencontre pas la sphere
     lintersect = .false.
     ri = n_rad+1
     thetaj = nz+1
     return
  endif
  lintersect = .true.
  rac=sqrt(delta)

  ! Les deux racines doivent etre positives sinon BUG !!
  ! s1 < s2
  s1=-b-rac
  !  s2=-b+rac

  ! TMP
  if (s1 < 0.0) then
     write(*,*) "Bug dans ray tracing !!!"
  endif
  ! END TMP

  delta_vol = s1

  ! Position au bord de la grille
  x1=x0+delta_vol*u
  y1=y0+delta_vol*v
  z1=z0+delta_vol*w

  ! Determination de l'indice de la premiere cellule traversee
  ! pour initialiser la propagation
  call indice_cellule_sph(x1,y1,z1,ri,thetaj)

  x=x1 ; y=y1 ; z=z1

  return


end subroutine move_to_grid_sph

!***********************************************************

subroutine integ_tau_mol(imol)

  implicit none

  integer, intent(in) :: imol

  real ::  norme, norme1, vmax, angle
  integer :: i, j, iTrans, n_speed, icell

  integer, dimension(2) :: ispeed
  real(kind=db), dimension(:), allocatable :: tab_speed, P


  n_speed = mol(imol)%n_speed_rt
  vmax = mol(imol)%vmax_center_rt

  allocate(tab_speed(-n_speed:n_speed), P(-n_speed:n_speed))

  ispeed(1) = -n_speed ; ispeed(2) = n_speed
  tab_speed(:) = span(-vmax,vmax,2*n_speed+1)

  angle=angle_interet

  iTrans = minval(mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing))

  norme=0.0
  norme1=0.0
  do i=1, n_rad
     icell = cell_map(i,1,1)
     P(:) = phiProf(icell,ispeed,tab_speed)
     norme=norme+kappa_mol_o_freq(icell,iTrans)*(r_lim(i)-r_lim(i-1))*P(0)
     norme1=norme1 + kappa(icell,1) * (r_lim(i)-r_lim(i-1))
  enddo
  write(*,*) "tau_mol = ", norme
  write(*,*) "tau_dust=", norme1

  loop_r : do i=1,n_rad
     icell = cell_map(i,1,1)
     if (r_grid(icell) > 100.0) then
        norme=0.0
        loop_z : do j=nz, 1, -1
           icell = cell_map(i,j,1)
           P(:) = phiProf(icell,ispeed,tab_speed)
           norme=norme+kappa_mol_o_freq(icell,1)*(z_lim(i,j+1)-z_lim(i,j))*p(0)
           if (norme > 1.0) then
              write(*,*) "Vertical Tau_mol=1 (for r=100AU) at z=", z_grid(icell), "AU"
              exit loop_z
           endif
        enddo loop_z
        exit loop_r
     endif
  enddo loop_r

  !read(*,*)

  return

end subroutine integ_tau_mol

!***********************************************************

subroutine length_deg2_opacity_wall(id,lambda,Stokes,ri,zj,xio,yio,zio,u,v,w,extrin,ltot,flag_sortie)
! Integration par calcul de la position de l'interface entre cellules
! par eq deg2 en r et deg 1 en z
! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
! C. Pinte
! 05/02/05

! ATTENTION : ce n'est qu'un mur analytique en opacite, pas une structure en densite
! donc, pas d'effet thermique, de dependance en lambda, etc

  implicit none

  integer, intent(in) :: id,lambda
  integer, intent(inout) :: ri,zj
  real(kind=db), dimension(4), intent(in) :: Stokes
  real(kind=db), intent(inout) :: u,v,w
  real, intent(in) :: extrin
  real(kind=db), intent(inout) :: xio,yio,zio
  real, intent(out) :: ltot
  logical, intent(out) :: flag_sortie

  real(kind=db) :: x0, y0, z0, x1, y1, z1
  real(kind=db) :: inv_a, a, b, c, s, rac, t, delta, inv_w, r_2
  real(kind=db) :: delta_vol, l, tau, zlim, extr, dotprod, opacite
  real(kind=db) :: correct_moins, correct_plus
  integer :: ri0, zj0, ri1, zj1, delta_rad, delta_zj, nbr_cell, icell0, icell

  logical :: lcellule_non_vide

  real(kind=db) :: hh, hhm, phi


  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  x1=xio;y1=yio;z1=zio
  extr=extrin
  ri0=ri
  zj0=zj
  ri1=ri
  zj1=zj


  ltot=0.0
  nbr_cell = 0

  a=u*u+v*v

  if (a > tiny_real) then
     inv_a=1.0_db/a
  else
     inv_a=huge_real
  endif


  if (abs(w) > tiny_real) then
     inv_w=1.0_db/w
  else
     inv_w=sign(huge_db,w) ! huge_real avant
  endif

  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1;zj0=zj1
     x0=x1;y0=y1;z0=z1

     icell0 = cell_map(ri0,zj0,1)

     ! Pour cas avec approximation de diffusion
     if (l_dark_zone(icell0)) then
        ! On revoie le paquet dans l'autre sens
        u = -u ; v = -v ; w=-w
        inv_w = -inv_w
        ! et on le renvoie au point de depart
        return
     endif

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        flag_sortie = .true.
        return
     elseif (zj0>nz) then
        lcellule_non_vide = .false.
        ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           flag_sortie = .true.
           return
        endif
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r_2=x0*x0+y0*y0
     b=(x0*u+y0*v)*inv_a

     if (ri0==0) then
        lcellule_non_vide = .false.
        opacite=0.0_db
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r_2-r_lim_2(0))*inv_a
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        delta_rad=1
     else
        if (icell0 > n_cells) then
           opacite = 0.0_db
        else
           opacite=kappa(icell0,lambda)
        endif

        if (ri0 == 1) then
           ! Variation de hauteur du mur en cos(phi/2)
           phi = atan2(y0,x0)
           hh = h_wall * abs(cos(phi/2.))
           hhm = -h_wall * abs(cos((phi+pi)/2.))

           ! Ajout de l'opacite du mur le cas echeant
           if ((z0 <= hh).and.(z0 >= hhm)) then
              opacite = opacite + kappa_wall
           endif
        endif

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod=u*x0+v*y0  ! ~ b
        if (dotprod < 0.0_db) then
           ! on recule : on cherche rayon inférieur
           c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
           delta=b*b-c
           if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0_db) then
           s=-b+rac
        else if (s==0.0_db) then
           s=prec_grille
        endif


        ! 2) position interface verticale
        ! on monte ou on descend par plan équatorial ?
        dotprod=w*z0
        if (dotprod == 0.0_db) then
           t=1.0e10
        else
           if (dotprod > 0.0_db) then
              ! on monte
              if (zj0==nz+1) then
                 delta_zj=0
                 if (z0 > 0.0_db) then
                    zlim=1.0e10
                 else
                    zlim=-1.0e10
                 endif
              else
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,zj0+1)*correct_plus
                 else
                    zlim=-z_lim(ri0,zj0+1)*correct_plus
                 endif
                 delta_zj=1
              endif
           else
              ! on descend
              if (zj0==1) then
                 ! on traverse le plan eq donc on va remonter
                 ! et z va changer de signe
                 delta_zj=1
                 if (z0 > 0.0_db) then
                    zlim=-z_lim(ri0,2)*correct_moins
                 else
                    zlim=z_lim(ri0,2)*correct_moins
                 endif
              else !(zj0==1)
                 ! on ne traverse pas z=0.
                 if (z0 > 0.0_db) then
                    zlim=z_lim(ri0,zj0)*correct_moins
                 else
                    zlim=-z_lim(ri0,zj0)*correct_moins
                 endif
                 delta_zj=-1
              endif !(zj0==1)
           endif ! monte ou descend
           t=(zlim-z0)*inv_w
           ! correct pb precision
           if (t < 0.0_db) t=prec_grille
        endif !dotprod=0.0
     endif ! ri0==0


     ! 3) interface en r ou z ?
     if (s < t) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        if ((ri1<1).or.(ri1>n_rad)) then
           zj1=zj0
        else
           zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
           if (zj1>nz) zj1=nz+1
        endif
     else ! z
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1=zj0+delta_zj
     endif



     ! Calcul longeur de vol et profondeur optique dans la cellule
     tau=l*opacite ! opacite constante dans la cellule
     ! Comparaison integrale avec tau
     if(tau > extr) then ! On a fini d'integrer
        l = l*extr/tau ! on rescale l pour que tau=extr
        ltot=ltot+l
        if (letape_th.and.lcellule_non_vide) then
           if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + kappa_abs_eg(icell0,lambda) * l * Stokes(1)
           if (lRE_nLTE.or.lnRE) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
        endif !l_abs
        flag_sortie = .false.
        xio=x0+l*u
        yio=y0+l*v
        zio=z0+l*w
        call indice_cellule(xio,yio,zio,ri,zj)

        icell = cell_map(ri,zj,1)
        if (l_dark_zone(icell)) then ! Petit test de securite
           ! On resort le paquet
           if (zj < zj0) then
              zj = zj0
              zio = z_lim(ri0,zj0)*correct_plus
           endif
           if (ri < ri0) then
              ri = ri0
              xio = xio * correct_plus
              yio = yio * correct_plus
           else if (ri > ri0) then
              ri = ri0
              xio = xio * correct_moins
              yio = yio * correct_moins
           endif
        endif
        return
     else ! Il reste extr - tau a integrer dans la cellule suivante
        extr=extr-tau
        ltot=ltot+l
        if (letape_th.and.lcellule_non_vide) then
           if (lRE_LTE) xKJ_abs(icell0,id) = xKJ_abs(icell0,id) + kappa_abs_eg(icell0,lambda) * l * Stokes(1)
           if (lRE_nLTE.or.lnRE) xJ_abs(icell0,lambda,id) = xJ_abs(icell0,lambda,id) + l * Stokes(1)
        endif !l_abs
     endif ! tau > extr

  enddo ! boucle infinie
  write(*,*) "BUG"
  return

end subroutine length_deg2_opacity_wall

!********************************************************************

function integ_ray_dust(lambda,ri,zj,phik,x,y,z,u,v,w)
  ! Generalisation de la routine length_deg2
  ! Modif depuis la routine moleculaire
  ! Propage un paquet depuis un point d'origine donne
  !
  ! C. Pinte
  ! 23/01/08

  implicit none

  integer, intent(in) :: lambda, ri,zj, phik
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: x,y,z

  real(kind=db), dimension(N_type_flux) :: integ_ray_dust

  if (lcylindrical) then
     integ_ray_dust(:) = integ_ray_dust_cyl(lambda,ri,zj,phik,x,y,z,u,v,w)
  else
     integ_ray_dust(:) = integ_ray_dust_sph(lambda,ri,zj,phik,x,y,z,u,v,w)
  endif

  return

end function integ_ray_dust

!***********************************************************

function integ_ray_dust_cyl(lambda,ri_in,zj_in,phik_in,x,y,z,u,v,w)
  ! Generalisation de la routine length_deg2
  ! Propage un paquet depuis un point d'origine donne
  ! et integre l'equation du transfert radiatif
  ! La propagation doit etre a l'envers pour faire du
  ! ray tracing  !!
  !
  ! C. Pinte
  ! 23/01/08

  ! TODO : faire peter le phi ??? Ne sert que pour les champs de vitesse

  implicit none

  integer, intent(in) :: lambda, ri_in, zj_in, phik_in
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: x,y,z

  real(kind=db), dimension(N_type_flux) :: integ_ray_dust_cyl

  real(kind=db) :: x0, y0, z0, x1, y1, z1, xm, ym, zm
  real(kind=db) :: inv_a, a, b, c, s, rac, t, t_phi, delta, inv_w, r_2, tan_angle_lim, den
  real(kind=db) :: delta_vol, l, zlim, dotprod
  real(kind=db) :: correct_plus, correct_moins
  integer :: ri0, zj0, ri1, zj1, phik0, phik1, delta_rad, delta_zj, nbr_cell, delta_phi, phik0m1

  real(kind=db) :: tau, dtau

  logical :: lcellule_non_vide

  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - prec_grille
  correct_plus = 1.0_db + prec_grille

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  ri1=ri_in
  zj1=zj_in
  phik1=phik_in

  nbr_cell = 0

  a=u*u+v*v

  if (a > tiny_real) then
     inv_a=1.0/a
  else
     inv_a=huge_real
  endif


  if (abs(w) > tiny_real) then
     inv_w=1.0/w
  else
     inv_w=sign(huge_db,w)
  endif

  tau = 0.0_db
  integ_ray_dust_cyl(:) = 0.0_db


  !*** propagation dans la grille

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1 ; zj0=zj1 ; phik0=phik1
     x0=x1 ; y0=y1 ; z0=z1

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        exit infinie
     elseif (abs(zj0)>nz) then
        lcellule_non_vide=.false.
        ! Test sortie vericale
        if (abs(z0) > zmaxmax) then
           exit infinie
        endif
     endif ! Test sortie

     nbr_cell = nbr_cell + 1

     ! Detection interface
     r_2=x0*x0+y0*y0
     b=(x0*u+y0*v)*inv_a

     if (ri0==0) then
        lcellule_non_vide=.false.
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        c=(r_2-r_lim_2(0)*correct_plus)*inv_a
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        t_phi= huge_real
        delta_rad=1
     else

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod=u*x0+v*y0
        if (dotprod < 0.0) then
           ! on recule : on cherche rayon inférieur
           c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
           delta=b*b-c
           if (delta < 0.0) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0) then
           s=-b+rac
        else if (s==0.0) then
           s=prec_grille
        endif

        ! 2) position interface verticale
        ! on monte ou on descend par plan équatorial ?
        dotprod=w*z0
        if (dotprod == 0.0_db) then
           t=1.0e10
        else
           if (dotprod > 0.0_db) then ! on se rapproche de la surface
              if (l3D) then
                 if (zj0==nz+1) then
                    delta_zj=0
                    zlim=1.0e10
                 else if (zj0==-(nz+1)) then
                    delta_zj=0
                    zlim=-1.0e10
                 else
                    if (z0 > 0.0) then
                       zlim=z_lim(ri0,zj0+1)*correct_plus
                       delta_zj=1
                    else
                       zlim=-z_lim(ri0,abs(zj0)+1)*correct_plus
                       delta_zj=-1
                    endif
                 endif
              else ! 2D
                 if (zj0==nz+1) then
                    delta_zj=0
                    if (z0 > 0.0_db) then
                       zlim=1.0e10
                    else
                       zlim=-1.0e10
                    endif
                 else
                    if (z0 > 0.0) then
                       zlim=z_lim(ri0,zj0+1)*correct_plus
                    else
                       zlim=-z_lim(ri0,abs(zj0)+1)*correct_plus
                    endif
                    delta_zj=1
                 endif
              endif ! l3D
           else ! on se rappoche du midplane
              if (l3D) then
                 if (z0 > 0.0) then
                    zlim=z_lim(ri0,abs(zj0))*correct_moins
                    delta_zj=-1
                    if (zj0==1) delta_zj=-2 ! pas d'indice 0
                 else
                    zlim=-z_lim(ri0,abs(zj0))*correct_moins
                    delta_zj=1
                    if (zj0==-1) delta_zj=2 ! pas d'indice 0
                 endif
              else ! 2D
                 if (zj0==1) then
                    ! on traverse le plan eq donc on va remonter
                    ! et z va changer de signe
                    delta_zj=1
                    if (z0 > 0.0_db) then
                       zlim=-z_lim(ri0,2)*correct_moins
                    else
                       zlim=z_lim(ri0,2)*correct_moins
                    endif
                 else !(zj0==1)
                    ! on ne traverse pas z=0.
                    if (z0 > 0.0_db) then
                       zlim=z_lim(ri0,zj0)*correct_moins
                    else
                       zlim=-z_lim(ri0,abs(zj0))*correct_moins
                    endif
                    delta_zj=-1
                 endif !(zj0==1)
              endif !l3D
           endif ! monte ou descend
           t=(zlim-z0)*inv_w
           ! correct pb precision
           if (t < 0.0_db) t=prec_grille
        endif !dotprod=0.0

        ! 3) position interface azimuthale
        dotprod =  x0*v - y0*u
        if (abs(dotprod) < 1.0e-10) then
           ! on ne franchit pas d'interface azimuthale
           t_phi = 1.0e30
        else
           ! Quelle cellule on va franchir
           if (dotprod > 0.0) then
              tan_angle_lim = tan_phi_lim(phik0)
              delta_phi=1
           else
              phik0m1=phik0-1
              if (phik0m1==0) phik0m1=N_az
              tan_angle_lim = tan_phi_lim(phik0m1)
              delta_phi=-1
           endif
           ! Longueur av interserction
           if (tan_angle_lim > 1.0d299) then
              t_phi = -x0/u
           else
              den= v-u*tan_angle_lim
              if (abs(den) > 1.0e-6) then
                 t_phi = -(y0-x0*tan_angle_lim)/den
              else
                 t_phi = 1.0e30
              endif
           endif
           if (t_phi < 0.0) t_phi = 1.0e30
        endif !dotprod = 0.0

     endif ! ri==0

     ! 4) interface en r ou z ou phi ?
     if ((s < t).and.(s < t_phi)) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        if ((ri1<1).or.(ri1>n_rad)) then
           zj1=zj0
        else
           zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
           if (zj1>nz) zj1=nz+1
           if (l3D) then
              if (z1 < 0.0) zj1=-zj1
           endif
        endif
        phik1=phik0

        ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
        if (ri0==0) call indice_cellule_3D_phi(x1,y1,z1,phik1)
     else if (t < t_phi) then ! z
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1=zj0+delta_zj
        phik1=phik0
     else
        l=t_phi
        delta_vol=correct_plus*t_phi
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        zj1= floor(abs(z1)/zmax(ri1)*nz) + 1
        if (zj1>nz) zj1=nz+1
        if (l3D) then
           if (z1 < 0.0) zj1=-zj1
        endif
        phik1=phik0+delta_phi
        if (phik1 == 0) phik1=N_az
        if (phik1 == N_az+1) phik1=1
     endif

     ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
     if (z1 == 0.0_db) then
        if (l3D) then
           z1 = sign(prec_grille,w)
        else ! 2D
           z1 = prec_grille
        endif ! l3D
     endif

     if (lcellule_non_vide) then
        ! Epaisseur optique de la cellule
        dtau =  l * kappa(cell_map(ri0,zj0,phik0),lambda)

        ! Fct source au milieu du parcours dans la cellule
        xm = 0.5 * (x0 + x1)
        ym = 0.5 * (y0 + y1)
        zm = 0.5 * (z0 + z1)

        ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
        ! la profondeur optique jusqu'a la cellule
        integ_ray_dust_cyl(:) = integ_ray_dust_cyl(:) + &
             exp(-tau) * (1.0_db - exp(-dtau)) * dust_source_fct(ri0,zj0,phik0,xm,ym,zm)

        ! Mise a jour profondeur optique pour cellule suivante
        tau = tau + dtau

        ! Pas besoin d'integrer trop profond
        if (tau > tau_dark_zone_obs) return
     endif  ! lcellule_non_vide

     ! TODO:stocker les ri0, zj0, l, xm, ym, zm pour chacun des rayons

  enddo infinie

  return

end function integ_ray_dust_cyl

!***********************************************************

function integ_ray_dust_sph(lambda,ri_in,thetaj_in,phik_in,x,y,z,u,v,w)
  ! Generalisation de la routine length_deg2
  ! Propage un paquet depuis un point d'origine donne
  ! et integre l'equation du transfert radiatif
  ! La propagation doit etre a l'envers pour faire du
  ! ray tracing  !!
  !
  ! C. Pinte
  ! 15/08/08

  implicit none

  integer, intent(in) :: lambda, ri_in, thetaj_in, phik_in
  real(kind=db), intent(in) :: u,v,w
  real(kind=db), intent(in) :: x,y,z

  real(kind=db), dimension(N_type_flux) :: integ_ray_dust_sph

  real(kind=db) :: x0, y0, z0, x1, y1, z1, xm, ym, zm, t1, t2
  real(kind=db) :: b, c, s, rac, t, t_phi, delta, tan_angle_lim, den
  real(kind=db) :: delta_vol, l, dotprod
  real(kind=db) :: correct_plus, correct_moins, precision
  real(kind=db) :: tan_angle_lim1, tan_angle_lim2, a_theta, b_theta, c_theta
  real(kind=db) :: t1_1, t1_2, t2_1, t2_2, uv, r0_cyl, r0_2_cyl, r0_2, tan2
  integer :: ri0, thetaj0, ri1, thetaj1, phik0, phik1, delta_rad, delta_theta, nbr_cell, delta_phi, phik0m1

  real(kind=db) :: tau, dtau

  logical :: lcellule_non_vide


  ! Petit delta pour franchir la limite de la cellule
  ! et ne pas etre pile-poil dessus
  correct_moins = 1.0_db - 1.0e-10_db
  correct_plus = 1.0_db + 1.0e-10_db
  precision = 1.0e-20_db ! pour g95

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  ri1=ri_in
  thetaj1=thetaj_in
  phik1=phik_in

  uv = sqrt(u*u + v*v)

  nbr_cell = 0

  tau = 0.0_db
  integ_ray_dust_sph(:) = 0.0_db

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     ri0=ri1 ; thetaj0=thetaj1 ; phik0=phik1
     x0=x1 ; y0=y1 ; z0=z1

     nbr_cell = nbr_cell + 1

     lcellule_non_vide=.true.
     ! Test sortie
     if (ri0>n_rad) then ! On est dans la derniere cellule
        ! Le photon sort du disque
        exit infinie
     endif ! Test sortie


     ! Detection interface
     r0_2_cyl  = x0*x0+y0*y0
     r0_cyl = sqrt(r0_2_cyl)
     r0_2 = r0_2_cyl + z0*z0
     b   = (x0*u+y0*v+z0*w)

     if (ri0==0) then
        ! Si on est avant le bord interne,  on passe forcement par rmin
        ! et on cherche forcement la racine positive (unique)
        lcellule_non_vide = .false.
        c=(r0_2-r_lim_2(0)*correct_plus)
        delta=b*b-c
        rac=sqrt(delta)
        s=-b+rac
        t=huge_real
        delta_rad=1
     else

        ! 1) position interface radiale
        ! on avance ou recule en r ? -> produit scalaire
        dotprod= b
        if (dotprod < 0.0_db) then
           ! on recule : on cherche rayon inférieur
           c=(r0_2-r_lim_2(ri0-1)*correct_moins)
           delta=b*b-c
           if (delta < 0.0_db) then ! on ne rencontre pas le rayon inférieur
              ! on cherche le rayon supérieur
              c=(r0_2-r_lim_2(ri0)*correct_plus)
              delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
              delta_rad=1
           else
              delta_rad=-1
           endif
        else
           ! on avance : on cherche le rayon supérieur
           c=(r0_2-r_lim_2(ri0)*correct_plus)
           delta=max(b*b-c,0.0_db) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
           delta_rad=1
        endif !dotprod
        rac=sqrt(delta)
        s=-b-rac
        if (s < 0.0_db) then
           s=-b+rac
        else if (s==0.0_db) then
           s=prec_grille
        endif


        ! 2) position interface inclinaison
        ! Meme methode que azimuth dans version cylindrique 3D
        if (z0 >= 0.0_db) then ! on est dans le bon sens
           tan_angle_lim1 = tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = tan_theta_lim(thetaj0-1) * correct_moins
        else
           tan_angle_lim1 = - tan_theta_lim(thetaj0) * correct_plus
           tan_angle_lim2 = - tan_theta_lim(thetaj0-1) * correct_moins
        endif ! z0

        ! Premiere limite theta
        tan2 = tan_angle_lim1 * tan_angle_lim1
        a_theta = w*w - tan2 * (u*u + v*v)
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t1 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t1_1 = (- b_theta - rac)/a_theta
           t1_2 = (- b_theta + rac)/a_theta

           if (t1_1 <= precision) then
              if (t1_2 <= precision) then ! les 2 sont <0
                 t1=1.0e30_db
              else   ! Seul t1_2 > 0
                 t1 = t1_2
              endif
           else
              if (t1_2 <= precision) then ! Seul t1_1 >0
                 t1 = t1_1
              else ! les 2 sont > 0
                 t1 = min(t1_1,t1_2)
              endif
           endif
        endif ! signe delta

        ! Deuxieme limite theta
        tan2 = tan_angle_lim2 * tan_angle_lim2
        a_theta = w*w - tan2 * (u*u + v*v)
        b_theta = w*z0 - tan2 * (x0*u + y0*v)
        c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

        delta = b_theta * b_theta - a_theta * c_theta
        if (delta < 0.0_db) then ! Pas de sol reelle
           t2 = 1.0e30_db
        else ! Au moins une sol reelle
           rac = sqrt(delta)
           t2_1 = (- b_theta - rac)/a_theta
           t2_2 = (- b_theta + rac)/a_theta

           if (t2_1 <= precision) then
              if (t2_2 <= precision) then ! les 2 sont <0
                 t2=1.0e30_db
              else   ! Seul t2_2 > 0
                 t2 = t2_2
              endif
           else
              if (t2_2 <= 0.0_db) then ! Seul t2_1 >0
                 t2 = t2_1
              else ! les 2 sont > 0
                 t2 = min(t2_1,t2_2)
              endif
           endif
        endif ! signe delta

        ! Selection limite theta
        if (t1 < t2) then
           t=t1
           delta_theta = 1
           if (thetaj0 == nz) delta_theta = 0
        else
           t=t2
           delta_theta = -1
           if (thetaj0 == 1) delta_theta = 0
        endif


        ! 3) position interface azimuthale
        dotprod =  x0*v - y0*u
        if (abs(dotprod) < 1.0e-10) then
           ! on ne franchit pas d'interface azimuthale
           t_phi = 1.0e30
           delta_phi = 0
        else
           ! Quelle cellule on va franchir
           if (dotprod > 0.0) then
              tan_angle_lim = tan_phi_lim(phik0)
              delta_phi=1
           else
              phik0m1=phik0-1
              if (phik0m1==0) phik0m1=N_az
              tan_angle_lim = tan_phi_lim(phik0m1)
              delta_phi=-1
           endif
           ! Longueur av interserction
           if (tan_angle_lim > 1.0d299) then
              t_phi = -x0/u
              delta_phi=0
           else
              den= v-u*tan_angle_lim
              if (abs(den) > 1.0e-6) then
                 t_phi = -(y0-x0*tan_angle_lim)/den
              else
                 t_phi = 1.0e30
                 delta_phi = 0
              endif
           endif
           if (t_phi < 0.0) then
              t_phi = 1.0e30
              delta_phi = 0
           endif
        endif !dotprod = 0.0

     endif ! ri0==0


     ! 4) interface en r ou z ou phi ?
     if ((s < t).and.(s < t_phi)) then ! r
        l=s
        delta_vol=s
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0+delta_rad
        thetaj1 = thetaj0
        phik1=phik0

        ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
        if (ri0==0) then
           call indice_cellule_sph_theta(x1,y1,z1,thetaj1)
        endif
     else if (t < t_phi) then ! theta
        l=t
        delta_vol=t
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        thetaj1=thetaj0+delta_theta
        phik1=phik0
     else
        l=t_phi
        delta_vol=correct_plus*t_phi
        ! Position au bord de la cellule suivante
        x1=x0+delta_vol*u
        y1=y0+delta_vol*v
        z1=z0+delta_vol*w
        ri1=ri0
        thetaj1 = thetaj0
        phik1=phik0+delta_phi
        if (phik1 == 0) phik1=N_az
        if (phik1 == N_az+1) phik1=1
     endif


     if (lcellule_non_vide) then
        ! Epaisseur optique de la cellule
        dtau =  l * kappa(cell_map(ri0,thetaj0,1),lambda)

        ! Fct source au milieu du parcours dans la cellule
        xm = 0.5 * (x0 + x1)
        ym = 0.5 * (y0 + y1)
        zm = 0.5 * (z0 + z1)

        ! Ajout emission en sortie de cellule (=debut car on va a l'envers) ponderee par
        ! la profondeur optique jusqu'a la cellule
        integ_ray_dust_sph(:) = integ_ray_dust_sph(:) + &
             exp(-tau) * (1.0_db - exp(-dtau)) * dust_source_fct(ri0,thetaj0,phik0,xm,ym,zm)

        ! Mise a jour profondeur optique pour cellule suivante
        tau = tau + dtau

        ! Pas besoin d'integrer trop profond
        if (tau > tau_dark_zone_obs) return
     endif

  enddo infinie

  return

end function integ_ray_dust_sph

!***********************************************************

subroutine angle_max(lambda)

  implicit none

  real, parameter :: tau_lim = 10.0

  integer, intent(in) :: lambda
  integer :: i, ri, zj

  real(kind=db) :: x0, y0, z0
  real(kind=db) :: u0, v0, w0
  real :: tau
  real(kind=db) ::  lmin, lmax, cos_max, cos_min
  real(kind=db), dimension(4) :: Stokes

  cos_max = sqrt(1.0-cos_max2)
  cos_min = 0.0

  v0= 0.0
  do i=1,20
     w0= 0.5*(cos_max + cos_min)
     u0=sqrt(1.0-w0*w0)
     ri = 0 ; zj=1
     x0=0.0 ; y0=0.0 ; z0=0.0
     Stokes = 0.0_db
     Stokes(1) = 1.0_db
     call length_deg2_tot(1,lambda,Stokes,ri,zj,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)
!     write(*,*) i, cos_min, w0, cos_max, tau
     if (tau > tau_lim) then
        cos_min = w0
     else
        cos_max = w0
     endif
  enddo

  w0_sup = sqrt(1.0-cos_max2)
  w0_inf = cos_min
 ! write(*,*)  w0_sup, w0_inf

  if (.not.lmono0) w0_inf=0.0

  return

end subroutine angle_max

!***********************************************************

subroutine define_dark_zone(lambda,tau_max,ldiff_approx)
! Definition l'etendue de la zone noire
! definie le tableau logique l_dark_zone
! et les rayons limites r_in_opacite pour le premier rayon
! C. Pinte
! 22/04/05

  implicit none

  integer, parameter :: nbre_angle = 11

  integer, intent(in) :: lambda
  real, intent(in) :: tau_max
  logical, intent(in) :: ldiff_approx
  integer :: i, j, pk, n, id, ri, zj, phik, icell, jj
  real(kind=db) :: x0, y0, z0, u0, v0, w0
  real :: somme, angle, dvol1, d_r, phi, r0

  logical :: flag_direct_star = .false.
  logical :: flag_star = .false.
  logical :: flag_sortie

  real(kind=db), dimension(4) :: Stokes

!  write(*,*) "defining DZ", tau_max

  do pk=1, n_az

     ri_in_dark_zone(pk)=n_rad
     ri_out_dark_zone(pk)=1
     ! étape 1 : radialement depuis le centre
     somme = 0.0
     do1 : do i=1,n_rad
        somme=somme+kappa(cell_map(i,1,pk),lambda)*(r_lim(i)-r_lim(i-1))
        if (somme > tau_max) then
           ri_in_dark_zone(pk) = i
           exit do1
        endif
     enddo do1

     ! étape 2 : radialement depuis rout
     somme = 0.0
     do2 : do i=n_rad,1,-1
        somme=somme+kappa(cell_map(i,1,pk),lambda)*(r_lim(i)-r_lim(i-1))
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
              somme=somme+kappa(cell_map(i,j,pk),lambda)*(z_lim(i,j+1)-z_lim(i,j))
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
                 somme=somme+kappa(cell_map(i,j,pk),lambda)*(z_lim(i,abs(j)+1)-z_lim(i,abs(j)))
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

  ! Cas premiere cellule
  r_in_opacite(:,:) = r_lim(1) ! bord externe de la premiere cellule
  r_in_opacite2(:,:) = r_lim_2(1)

  do pk=1, n_az
     if (ri_in_dark_zone(pk)==1) then
        do j=1, zj_sup_dark_zone(ri_in_dark_zone(pk),pk)
           d_r=tau_max/kappa(cell_map(1,j,pk),lambda)
           r_in_opacite(j,pk) = (rmin + d_r)
           r_in_opacite2(j,pk) = r_in_opacite(j,pk)**2
        enddo

        if (l3D) then
           do j=-1, zj_inf_dark_zone(ri_in_dark_zone(pk),pk), -1
              d_r=tau_max/kappa(cell_map(1,j,pk),lambda)
              r_in_opacite(j,pk) = (rmin + d_r)
              r_in_opacite2(j,pk) = r_in_opacite(j,pk)**2
           enddo
        endif

     endif
  enddo

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
              ! test
              ri=i
              zj=j
              phik=1
              Stokes(:) = 0.0_db ; !Stokes(1) = 1.0_db ; ! Pourquoi c'etait a 1 ?? ca fausse les chmps de radiation !!!
              call length_deg2(id,lambda,Stokes,ri,zj,phik,x0,y0,z0,u0,v0,w0,flag_star,flag_direct_star,tau_max,dvol1,flag_sortie)
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
                 ! test
                 ri=i
                 zj=j
                 phik=pk
                 Stokes(:) = 0.0_db ; Stokes(1) = 1.0_db ;
                 call length_deg2_3D(id,lambda,Stokes,ri,zj,phik,x0,y0,z0,u0,v0,w0,flag_star,tau_max,dvol1,flag_sortie)
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
                 ! test
                 ri=i
                 zj=j
                 phik=pk
                 Stokes(:) = 0.0_db ; Stokes(1) = 1.0_db ;
                 call length_deg2_3D(id,lambda,Stokes,ri,zj,phik,x0,y0,z0,u0,v0,w0,flag_star,tau_max,dvol1,flag_sortie)
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
     if (minval(ri_in_dark_zone(:))==1) then
        write(*,*) "WARNING : first cell is in diffusion approximation zone"
        write(*,*) "Increase spatial grid resolution"
        stop
     endif
  endif


  if (n_zones > 1) then
     do icell=1, n_cells
        if (sum(densite_pouss(icell,:)) < tiny_real) l_dark_zone(icell) = .false.
     enddo
  endif

  do i=1, n_regions
     do j=1,nz
        l_dark_zone(cell_map(regions(i)%iRmin,j,1)) = .false.
        l_dark_zone(cell_map(regions(i)%iRmax,j,1)) = .false.
     enddo
  enddo

!  write(*,*) l_dark_zone(:,nz,1)
!  write(*,*) "**"
 ! write(*,*) l_dark_zone(10,:,1)

  !  do i=1, n_rad
  !     do j=1,nz
  !        write(*,*) i, j, l_dark_zone(i,j), r_lim(1)
  !     enddo
  !  enddo
  !  write(*,*) sqrt(r_in_opacite2(1))
  !  read(*,*)

  return

end subroutine define_dark_zone

!***********************************************************

subroutine no_dark_zone()
! Definie les variables quand il n'y a pas de zone noire
! C . Pinte
! 22/04/05

  implicit none

  r_in_opacite(:,:) = r_lim(1)
  r_in_opacite2(:,:) = r_lim_2(1)

  l_dark_zone(:)=.false.

  return

end subroutine no_dark_zone

!***********************************************************

logical function test_dark_zone(ri,zj,phik,x,y)
! Test si on est dans la zone noire
! C. Pinte
! 22/04/05

  implicit none

  integer, intent(in) :: ri, zj, phik
  real(kind=db), intent(in) :: x, y

  integer :: icell

  if (ri==1) then
     if (x*x+y*y > r_in_opacite2(zj,phik)) then
        test_dark_zone = .true.
     else
        test_dark_zone = .false.
     endif
  else
     icell = cell_map(ri,zj,phik)
     test_dark_zone = l_dark_zone(icell)
  endif


  return

end function  test_dark_zone

!***********************************************************

subroutine define_proba_weight_emission(lambda)
  ! Augmente le poids des cellules pres de la surface
  ! par exp(-tau)
  ! Le poids est applique par weight_repartion_energie
  ! C. Pinte
  ! 19/11/08

  implicit none


  integer, intent(in) :: lambda

  real, dimension(n_rad,nz) :: tau_min
  real(kind=db), dimension(4) :: Stokes
  real(kind=db) :: x0, y0, z0, u0, v0, w0, angle, lmin, lmax
  real :: tau
  integer :: i, j, n, id, ri, zj, icell
  integer, parameter :: nbre_angle = 101

  tau_min(:,:) = 1.e30 ;

  do i=1,n_rad
     do j=1,nz
        do n=1,nbre_angle
           id=1
           ! position et direction vol
           angle= pi * real(n)/real(nbre_angle+1)! entre 0 et pi
           x0=1.00001*r_lim(i-1) ! cellule 1 traitee a part
           y0=0.0
           z0=0.99999*z_lim(i,j+1)
           u0=cos(angle)
           v0=0.0
           w0=sin(angle)

           ri=i ;  zj=j ; Stokes(:) = 0.0_db ;
           call length_deg2_tot(id,lambda,Stokes,ri,zj,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)
           if (tau < tau_min(i,j)) tau_min(i,j) = tau

           x0 = 0.99999*r_lim(i)
           call length_deg2_tot(id,lambda,Stokes,ri,zj,x0,y0,y0,u0,v0,w0,tau,lmin,lmax)
           if (tau < tau_min(i,j)) tau_min(i,j) = tau

        enddo
     enddo
  enddo

  do i=1,n_rad
     do j=1,nz
        icell = cell_map(i,j,1)
        weight_proba_emission(icell) =  exp(-tau_min(i,j))
     enddo
  enddo

  ! correct_E_emission sera normalise dans repartition energie
  correct_E_emission(:) = 1.0_db / weight_proba_emission(:)

  return

end subroutine define_proba_weight_emission

!***********************************************************

end module optical_depth
