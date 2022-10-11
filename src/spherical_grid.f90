module spherical_grid

  use constantes
  use parametres
  use cylindrical_grid, only: cell2cylindrical, cell_map, cell_map_i, cell_map_j, cell_map_k, lexit_cell, r_lim, r_lim_2, r_lim_3, &
       delta_z, dr2_grid, r_grid, z_grid, phi_grid, tab_region, z_lim, w_lim, theta_lim, tan_theta_lim, tan_phi_lim, &
       volume, l_dark_zone, zmax
  use messages

  implicit none

  public :: cross_spherical_cell, pos_em_cellule_sph, indice_cellule_sph, test_exit_grid_sph, &
       move_to_grid_sph


  private

  real(kind=dp), parameter :: prec_grille=1.0e-14_dp
  real(kind=dp), parameter :: prec_grille_sph=1.0e-7_dp


contains

  pure logical function test_exit_grid_sph(icell, x, y, z)

    integer, intent(in) :: icell
    real(kind=dp), intent(in) :: x,y,z

    if (icell <= n_cells) then
       test_exit_grid_sph = .false.
       return
    endif

    if (lexit_cell(icell)==1) then ! radial
       test_exit_grid_sph = .true.
    else
       test_exit_grid_sph = .false.
    endif

    return

    return

  end function test_exit_grid_sph

!******************************************************************************

  subroutine indice_cellule_sph(xin,yin,zin, icell)

  implicit none

  real(kind=dp), intent(in) :: xin,yin,zin
  integer, intent(out) :: icell

  real(kind=dp) :: r2, r02, tan_theta, phi
  integer :: ri, ri_min, ri_max, thetaj, thetaj_min, thetaj_max, ri_out, thetaj_out, phik_out

  r02 = xin*xin+yin*yin
  r2 = r02+zin*zin

  ! Peut etre un bug en 0, 0 due a la correction sur grid_rmin dans define_grid3
  if (r2 < r_lim_2(0)) then
     ri_out=0
     thetaj_out=1
     phik_out=1
  else if (r2 > Rmax2) then
     ri_out=n_rad+1
     thetaj_out=1
     phik_out=1
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

     ! thetaj_out
     if (r02 > tiny_dp) then
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

     if (l3D) then
        if (zin < 0.) thetaj_out = -thetaj_out

        if (zin /= 0.0) then
           phi=modulo(atan2(yin,xin),2*real(pi,kind=dp))
           phik_out=floor(phi/(2*pi)*real(N_az))+1
           if (phik_out==n_az+1) phik_out=n_az
        else
           phik_out=1
        endif
     else ! 2D
        phik_out = 1
     endif
  endif ! r < rmin or r > rmax

  icell = cell_map(ri_out,thetaj_out,phik_out)

  return

end subroutine indice_cellule_sph

!******************************************************************************

subroutine indice_cellule_sph_theta(xin,yin,zin,thetaj_out,phik_out)

  implicit none

  real(kind=dp), intent(in) :: xin,yin,zin
  integer, intent(out) :: thetaj_out, phik_out

  real(kind=dp) :: r02, tan_theta, phi
  integer :: thetaj, thetaj_min, thetaj_max

  r02 = xin*xin+yin*yin

  ! thetaj_out
  if (r02 > tiny_dp) then
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

  if (l3D) then
     if (zin < 0) thetaj_out = - thetaj_out

     if (zin /= 0.0) then
        phi=modulo(atan2(yin,xin),2*real(pi,kind=dp))
        phik_out=floor(phi/(2*pi)*real(N_az))+1
        if (phik_out==n_az+1) phik_out=n_az
     else
        phik_out=1
     endif
  else ! 2D
     phik_out=1
  endif

  return

end subroutine indice_cellule_sph_theta

!******************************************************************************


  subroutine cross_spherical_cell(x0,y0,z0, u,v,w,  cell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

    integer, intent(in) :: cell, previous_cell
    real(kind=dp), intent(in) :: x0,y0,z0
    real(kind=dp), intent(in) :: u,v,w ! Todo : check that

    real(kind=dp), intent(out) :: x1, y1, z1
    integer, intent(out) :: next_cell
    real(kind=dp), intent(out) :: l, l_contrib, l_void_before


    real(kind=dp) :: correct_moins, correct_plus, uv, precision
    real(kind=dp) :: b, c, s, rac, t, delta, r0_2, r0_cyl, r0_2_cyl
    real(kind=dp) :: delta_vol, dotprod, t_phi, tan_angle_lim, den
    integer :: ri0, thetaj0, ri1, thetaj1, delta_rad, delta_theta

    integer :: phik0, phik1, delta_phi, phik0m1

    real(kind=dp) :: a_theta, b_theta, c_theta, tan2, tan_angle_lim1, tan_angle_lim2, t1, t2
    real(kind=dp) :: t1_1, t1_2, t2_1, t2_2, a_theta_m1

    ! Todo : can be calculated outside

    ! Petit delta pour franchir la limite de la cellule
    ! et ne pas etre pile-poil dessus
    correct_moins = 1.0_dp - prec_grille_sph
    correct_plus = 1.0_dp + prec_grille_sph
    precision = 1.0e-20_dp ! pour g95

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
       ! Si on est avant le bord interne,  on passe forcement par rmin
       ! et on cherche forcement la racine positive (unique)
       c=(r0_2-r_lim_2(0)*correct_plus)
       delta=b*b-c
       rac=sqrt(delta)
       s = (-b+rac) * correct_plus
       t=huge_real
       delta_rad=1
       t_phi = huge_real
    else
       ! 1) position interface radiale
       ! on avance ou recule en r ? -> produit scalaire
       dotprod = b
       if (dotprod < 0.0_dp) then
          ! on recule : on cherche rayon inf�rieur
          c=(r0_2-r_lim_2(ri0-1)*correct_moins)
          delta=b*b-c
          if (delta < 0.0_dp) then ! on ne rencontre pas le rayon inf�rieur
             ! on cherche le rayon sup�rieur
             c=(r0_2-r_lim_2(ri0)*correct_plus)
             delta=max(b*b-c,0.0_dp) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
             delta_rad=1
          else
             delta_rad=-1
          endif
       else
          ! on avance : on cherche le rayon sup�rieur
          c=(r0_2-r_lim_2(ri0)*correct_plus)
          delta=max(b*b-c,0.0_dp) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
          delta_rad=1
       endif !dotprod
       rac=sqrt(delta)
       s=-b-rac
       if (s < 0.0_dp) then
          s=-b+rac
       else if (s==0.0_dp) then
          s=prec_grille
       endif

       ! 2) position interface inclinaison
       if (z0 >= 0.0_dp) then ! on est dans le bon sens
          tan_angle_lim1 = tan_theta_lim(abs(thetaj0)) * correct_plus
          tan_angle_lim2 = tan_theta_lim(abs(thetaj0)-1) * correct_moins
       else
          tan_angle_lim1 = - tan_theta_lim(abs(thetaj0)) * correct_plus
          tan_angle_lim2 = - tan_theta_lim(abs(thetaj0)-1) * correct_moins
       endif ! z0

       ! Premiere limite theta
       tan2 = tan_angle_lim1 * tan_angle_lim1
       a_theta = w*w - tan2 * (u*u + v*v)
       a_theta_m1 = 1.0_dp/a_theta
       b_theta = w*z0 - tan2 * (x0*u + y0*v)
       c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

       delta = b_theta * b_theta - a_theta * c_theta
       if (delta < 0.0_dp) then ! Pas de sol reelle
          t1 = 1.0e30_dp
       else ! Au moins une sol reelle
          rac = sqrt(delta)
          t1_1 = (- b_theta - rac) * a_theta_m1
          t1_2 = (- b_theta + rac) * a_theta_m1

          if (t1_1 <= precision) then
             if (t1_2 <= precision) then ! les 2 sont <0
                t1=1.0e30_dp
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
       a_theta_m1 = 1.0_dp/a_theta
       b_theta = w*z0 - tan2 * (x0*u + y0*v)
       c_theta = z0*z0 - tan2 * (x0*x0 + y0*y0)

       delta = b_theta * b_theta - a_theta * c_theta
       if (delta < 0.0_dp) then ! Pas de sol reelle
          t2 = 1.0e30_dp
       else ! Au moins une sol reelle
          rac = sqrt(delta)
          t2_1 = (- b_theta - rac) * a_theta_m1
          t2_2 = (- b_theta + rac) * a_theta_m1

          if (t2_1 <= precision) then
             if (t2_2 <= precision) then ! les 2 sont <0
                t2=1.0e30_dp
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
          if (abs(thetaj0) == nz) delta_theta = 0
       else
          t=t2
          delta_theta = -1
          if (abs(thetaj0) == 1) delta_theta = 0
       endif

       ! 3) position interface azimuthale
       if (l3D) then
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
             if (tan_angle_lim > 1.0d299) then ! inclus le cas 90 deg
                t_phi = -x0/u
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
       else ! l3D
          t_phi = huge_real
       endif
    endif ! ri0==0


    ! 4) interface en r ou theta ou phi ?
    if ((s < t).and.(s < t_phi)) then ! r
       l=s
       delta_vol=s
       ! Position au bord de la cellule suivante
       x1=x0+delta_vol*u
       y1=y0+delta_vol*v
       z1=z0+delta_vol*w
       ri1=ri0+delta_rad
       thetaj1 = thetaj0
       phik1 = phik0

       ! Calcul de l'indice theta quand on rentre dans la cellule ri=1
       ! Todo : calcul indice phi aussi ???
       if (ri0 == 0) then
          call indice_cellule_sph_theta(x1,y1,z1,thetaj1,phik1)
       endif

       ! Forcing indices to be 1 if ri1==0, as this defines the cell where the star is
       if (ri1 == 0) then
          thetaj1 = 1
          phik1 = 1
       endif

    else if (t < t_phi) then ! theta
       l=t
       delta_vol=t
       ! Position au bord de la cellule suivante
       x1=x0+delta_vol*u
       y1=y0+delta_vol*v
       z1=z0+delta_vol*w
       ri1=ri0
       thetaj1=abs(thetaj0)+delta_theta
       if (l3D) then
          if (z1<0) thetaj1 = -thetaj1
       endif
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

    ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
    if (z1 == 0.0_dp) z1 = prec_grille

    !call cylindrical2cell(ri1,zj1,1, next_cell)
    next_cell = cell_map(ri1,thetaj1,phik1)

    l_contrib = l
    l_void_before = 0.0_dp

    return

  end subroutine cross_spherical_cell

!***********************************************************

  subroutine verif_cell_position_sph(icell, x, y, z)

    real(kind=dp), intent(inout) :: x,y,z
    integer, intent(inout) :: icell

    integer :: ri, zj, ri0, zj0, tmp_k
    real(kind=dp) :: factor, correct_moins, correct_plus

    correct_moins = 1.0_dp - prec_grille_sph
    correct_plus = 1.0_dp + prec_grille_sph

    ! tmp :
    call cell2cylindrical(icell, ri0,zj0, tmp_k) ! converting current cell index

    ! locate current cell index
    call indice_cellule_sph(x,y,z, icell)
    ri = cell_map_i(icell)

    ! Patch pour eviter BUG sur position radiale
    ! a cause de limite de precision
    if (ri==0) then
       factor = rmin/ sqrt(x*x+y*y+z*z) * correct_plus
       x = x * factor
       y = y * factor
       z = z * factor

       ! On verifie que c'est OK maintenant
       call indice_cellule_sph(x,y,z, icell)
       ri = cell_map_i(icell)
       if (ri==0) then
          call error("BUG in verif_cell_position_sph")
          !write(*,*) x,y,z, sqrt(x*x+y*y+z*z)
       endif
    endif

    if (l_dark_zone(icell)) then ! Petit test de securite
       zj = cell_map_j(icell)
       ! On resort le paquet
       if (zj < zj0) then
          zj = zj0
          z = z_lim(ri0,zj0)*correct_plus
       endif
       if (ri < ri0) then
          ri = ri0
          x = x * correct_plus
          y = y * correct_plus
       else if (ri > ri0) then
          ri = ri0
          x = x * correct_moins
          y = y * correct_moins
       endif
    endif

    return

  end subroutine verif_cell_position_sph

!***********************************************************

  subroutine move_to_grid_sph(id, x,y,z,u,v,w, icell,lintersect)
    ! Calcule la position au bord de la grille dans
    ! la direction donnee pour grille spherique
    ! C. Pinte
    ! 01/08/07

    implicit none

    integer, intent(in) :: id
    real(kind=dp), intent(inout) :: x,y,z
    real(kind=dp), intent(in) :: u,v,w
    integer, intent(out) :: icell
    logical, intent(out) :: lintersect

    real(kind=dp) :: x0, y0, z0, x1, y1, z1, r0_2, b, c, rac, delta, s1, delta_vol, correct_moins

    correct_moins = 1.0_dp - 1.0e-10_dp

    x0=x ; y0=y ; z0=z

    r0_2 = x0*x0+y0*y0+z0*z0
    b   = (x0*u+y0*v+z0*w)
    c=(r0_2-r_lim_2(n_rad)*correct_moins)
    delta=b*b-c

    if (delta < 0.0_dp) then
       ! On ne rencontre pas la sphere
       lintersect = .false.
       icell = 0
       return
    endif
    lintersect = .true.
    rac=sqrt(delta)

    ! Les deux racines doivent etre positives sinon BUG !!
    ! s1 < s2
    s1=-b-rac
    !  s2=-b+rac

    delta_vol = s1

    ! Position au bord de la grille
    x1=x0+delta_vol*u
    y1=y0+delta_vol*v
    z1=z0+delta_vol*w

    ! Determination de l'indice de la premiere cellule traversee
    ! pour initialiser la propagation
    call indice_cellule_sph(x1,y1,z1, icell)
    x=x1 ; y=y1 ; z=z1

    return

  end subroutine move_to_grid_sph

!***********************************************************

  subroutine pos_em_cellule_sph(icell ,aleat1,aleat2,aleat3,x,y,z)
! Choisit la position d'emission uniformement
! dans la cellule (ri,thetaj)
! Geometrie spherique
! C. Pinte
! 08/06/07

  implicit none

  integer, intent(in) :: icell
  real, intent(in) :: aleat1, aleat2, aleat3
  real(kind=dp), intent(out) :: x,y,z

  real(kind=dp) :: r, theta, phi, r_cos_theta
  integer :: ri, thetaj, phik

  ri = cell_map_i(icell)
  thetaj = cell_map_j(icell)
  phik = cell_map_k(icell)

  ! Position radiale
  r=(r_lim_3(ri-1)+aleat1*(r_lim_3(ri)-r_lim_3(ri-1)))**un_tiers

  ! Position theta
  if (l3D) then
     theta = theta_lim(abs(thetaj)-1)+aleat2*(theta_lim(abs(thetaj))-theta_lim(abs(thetaj)-1))
     !if (thetaj < 0) theta = - theta
  else
     if (aleat2 > 0.5_dp) then
        theta = theta_lim(thetaj-1)+(2.0_dp*(aleat2-0.5_dp))*(theta_lim(thetaj)-theta_lim(thetaj-1))
     else
        theta = -(theta_lim(thetaj-1)+(2.0_dp*aleat2)*(theta_lim(thetaj)-theta_lim(thetaj-1)))
     endif
  endif
  ! BUG ??? : ca doit etre uniforme en w, non ??

  ! Position azimuthale
  phi = 2.0_dp*pi * (real(phik)-1.0_dp+aleat3)/real(n_az)

  ! x et y
  z=r*sin(theta)
  r_cos_theta = r*cos(theta)
  x=r_cos_theta*cos(phi)
  y=r_cos_theta*sin(phi)

!!$  ! Position theta
!!$  if (aleat2 > 0.5) then
!!$     w=w_lim(thetaj-1)+(2.0_dp*(aleat2-0.5_dp))*(w_lim(thetaj)-w_lim(thetaj-1))
!!$  else
!!$     w=-(w_lim(thetaj-1)+(2.0_dp*aleat2)*(w_lim(thetaj)-w_lim(thetaj-1)))
!!$  endif
!!$
!!$  ! Position azimuthale
!!$  phi = 2.0_dp*pi * (real(phik)-1.0_dp+aleat3)/real(n_az)
!!$
!!$  ! x et y
!!$  z=r*w
!!$  r_cos_theta = r*sqrt(1.0_dp-w*w)
!!$  x=r_cos_theta*cos(phi)
!!$  y=r_cos_theta*sin(phi)

!  z = z_grid(ri,thetaj)
!  r = r_grid(ri,thetaj)
!  x = r*cos(phi)
!  y = r*sin(phi)

 ! call indice_cellule_sph(x,y,z,ri_tmp,thetaj_tmp)
 ! if (ri /= ri_tmp) then
 !    write(*,*) "Bug ri", ri, ri_tmp, sqrt(x*x+y*y)
 !    read(*,*)
 ! else if (thetaj /= thetaj_tmp) then
 !    write(*,*) "Bug zj", w, thetaj, thetaj_tmp
 !    call indice_cellule_sph(x,y,-z,ri_tmp,thetaj_tmp)
 !    write(*,*) -z, thetaj_tmp
 !    read(*,*)
 ! endif


  return

end subroutine pos_em_cellule_sph

end module spherical_grid
