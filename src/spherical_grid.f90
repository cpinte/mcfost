module spherical_grid

  use mem
  use cylindrical_grid, only: cell2cylindrical

  implicit none

contains

  subroutine indice_cellule_sph(xin,yin,zin,ri_out,thetaj_out)

  implicit none

  real(kind=db), intent(in) :: xin,yin,zin
  integer, intent(out) :: ri_out, thetaj_out

  real(kind=db) :: r2, r02, tan_theta
  integer :: ri, ri_min, ri_max, thetaj, thetaj_min, thetaj_max

  r02 = xin*xin+yin*yin
  r2 = r02+zin*zin

  ! Peut etre un bug en 0, 0 due a la correction sur grid_rmin dans define_grid3
  if (r2 < r_lim_2(0)) then
     ri_out=0
  else if (r2 > Rmax2) then
     ri_out=n_rad
  else
     ri_min=0
     ri_max=n_rad
     ri=(ri_min+ri_max)/2

     do while((ri_max-ri_min) > 1)
        if(r2 > r_lim_2(ri)) then
           ri_min=ri
        else
           ri_max=ri
        endif
        ri=(ri_min+ri_max)/2
     enddo
     ri_out=ri+1
  endif

  ! thetaj_out
  if (r02 > tiny_db) then
     tan_theta = abs(zin) / sqrt(r02)
  else
     tan_theta = 1.0e30
  endif

  thetaj_min = 0
  thetaj_max = nz
  thetaj=(thetaj_min+thetaj_max)/2

  do while((thetaj_max-thetaj_min) > 1)
     if(tan_theta > tan_theta_lim(thetaj)) then
        thetaj_min=thetaj
     else
        thetaj_max=thetaj
     endif
     thetaj=(thetaj_min+thetaj_max)/2
  enddo
  thetaj_out=thetaj+1

  return

end subroutine indice_cellule_sph

!******************************************************************************

subroutine indice_cellule_sph_theta(xin,yin,zin,thetaj_out)

  implicit none

  real(kind=db), intent(in) :: xin,yin,zin
  integer, intent(out) :: thetaj_out

  real(kind=db) :: r02, tan_theta
  integer :: thetaj, thetaj_min, thetaj_max

  r02 = xin*xin+yin*yin

  ! thetaj_out
  if (r02 > tiny_db) then
     tan_theta = abs(zin) / sqrt(r02)
  else
     tan_theta = 1.0e30
  endif

  thetaj_min = 0
  thetaj_max = nz
  thetaj=(thetaj_min+thetaj_max)/2

  do while((thetaj_max-thetaj_min) > 1)
     if(tan_theta > tan_theta_lim(thetaj)) then
        thetaj_min=thetaj
     else
        thetaj_max=thetaj
     endif
     thetaj=(thetaj_min+thetaj_max)/2
  enddo
  thetaj_out=thetaj+1

  return

end subroutine indice_cellule_sph_theta

!******************************************************************************

  subroutine cross_spherical_cell(lambda, x0,y0,z0, u,v,w,  cell, previous_cell, x1,y1,z1, next_cell, l, tau)

    integer, intent(in) :: lambda, cell, previous_cell
    real(kind=db), intent(in) :: x0,y0,z0
    real(kind=db), intent(in) :: u,v,w ! Todo : check that

    real(kind=db), intent(out) :: x1, y1, z1
    integer, intent(out) :: next_cell
    real(kind=db), intent(out) :: l, tau


    real(kind=db) :: correct_moins, correct_plus, uv, precision
    real(kind=db) :: b, c, s, rac, t, delta, r0_2, r0_cyl, r0_2_cyl
    real(kind=db) :: delta_vol, dotprod, opacite
    integer :: ri0, thetaj0, ri1, thetaj1, delta_rad, delta_theta, nbr_cell, p_ri0, p_thetaj0, icell0

    integer :: phik0

    real(kind=db) :: a_theta, b_theta, c_theta, tan2, tan_angle_lim1, tan_angle_lim2, t1, t2
    real(kind=db) :: t1_1, t1_2, t2_1, t2_2, a_theta_m1

    ! Todo : can be calculated outside

    ! Petit delta pour franchir la limite de la cellule
    ! et ne pas etre pile-poil dessus
    correct_moins = 1.0_db - prec_grille_sph
    correct_plus = 1.0_db + prec_grille_sph
    precision = 1.0e-20_db ! pour g95

    uv = sqrt(u*u + v*v)
    ! end     ! Todo : can be calculated outside

    ! 3D cell indices
    call cell2cylindrical(cell, ri0,thetaj0,phik0)

    ! Detection interface
    r0_2_cyl  = x0*x0+y0*y0
    r0_cyl = sqrt(r0_2_cyl)
    r0_2 = r0_2_cyl + z0*z0
    b   = (x0*u+y0*v+z0*w)

    if (ri0==0) then
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

    !call cylindrical2cell(ri1,zj1,1, next_cell)
    next_cell =  cell_map(ri1,thetaj1,1)

    return

  end subroutine cross_spherical_cell

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

end module spherical_grid
