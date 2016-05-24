module cylindrical_grid

  use mem

  implicit none

contains

  subroutine build_cylindrical_cell_mapping() ! work also in sph

    integer :: i,j,k,icell, ntot, ntot2, alloc_status
    integer :: istart,iend,jstart,jend,kstart,kend, istart2,iend2,jstart2,jend2,kstart2,kend2

    icell_ref = 1

    istart = 1
    iend = n_rad

    jstart = j_start
    jend = nz

    kstart=1
    kend = n_az

    if (j_start < 0) then
       ntot = (iend - istart + 1) * (jend - jstart) * (kend - kstart + 1)
    else
       ntot = (iend - istart + 1) * (jend - jstart +1) * (kend - kstart + 1)
    endif

    if (ntot /= n_cells) then
       write(*,*) "ERROR in 'build_cylindrical_cell_mapping'"
       write(*,*) "The number of cells is not matching :"
       write(*,*) "ntot=", ntot, "should be", n_cells
       write(*,*) "Exiting."
       stop
    endif

    istart2 = 0
    iend2 = n_rad + 1

    jstart2 = min(1,j_start)-1
    jend2 = nz+1

    kstart2=1
    kend2 = n_az

    if (jstart2 < 0) then
       ntot2 = (iend2 - istart2 + 1) * (jend2 - jstart2) * (kend2 - kstart2 + 1)
    else
       ntot2 = (iend2 - istart2 + 1) * (jend2 - jstart2 +1) * (kend2 - kstart2 + 1)
    endif
    allocate(cell_map(istart2:iend2,jstart2:jend2,kstart2:kend2))
    allocate(cell_map_i(ntot2), cell_map_j(ntot2), cell_map_k(ntot2))

    ! Actual cells
    icell = 0
    do k=kstart, kend
       bz : do j=j_start, jend
          if (j==0) cycle bz
          do i=istart, iend

             icell = icell+1
             if (icell > ntot) then
                write(*,*) "ERROR : there is an issue in the cell mapping"
                write(*,*) "Exiting"
                stop
             endif

             cell_map_i(icell) = i
             cell_map_j(icell) = j
             cell_map_k(icell) = k

             cell_map(i,j,k) = icell
          enddo
       enddo bz
    enddo

    if (icell /= ntot) then
       write(*,*) "Something went wrong in the call mapping"
       write(*,*) "I am missing some real cells"
       write(*,*) icell, ntot
       write(*,*)
       stop
    endif


    ! Virtual cell indices for when the packets are just around the grid

    ! Can the packet exit from this cell : 0 -> no, 1 -> radially, 2 -> vertically
    allocate(lexit_cell(1:ntot2), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error lexit_cell'
       stop
    endif
    lexit_cell(:) = 0

    ! Cases j=0 and j=nz+1
    do k=kstart, kend
       do j = jstart2, jend2, jend2 - jstart2
          do i=istart2, iend2

             icell = icell+1
             if (icell > ntot2) then
                write(*,*) "ERROR : there is an issue in the cell mapping"
                write(*,*) "Exiting"
                stop
             endif

             if (j==jend2) lexit_cell(icell) = 2
             if (i==iend2) lexit_cell(icell) = 1

             cell_map_i(icell) = i
             cell_map_j(icell) = j
             cell_map_k(icell) = k

             cell_map(i,j,k) = icell
          enddo
       enddo
    enddo

    ! Cases i=0 and i=n_rad+1 (except j=0 and j=nz+1 done above)
    do k=kstart, kend
       bz2 : do j = jstart, jend
          if (j==0) cycle bz2
          do i=istart2,iend2, iend2-istart2

             icell = icell+1
             if (icell > ntot2) then
                write(*,*) "ERROR : there is an issue in the cell mapping"
                write(*,*) "Extra cells:", icell, ntot2
                write(*,*) i,j,k
                write(*,*) "Exiting"
                stop
             endif

             if (i==iend2) lexit_cell(icell) = 1

             cell_map_i(icell) = i
             cell_map_j(icell) = j
             cell_map_k(icell) = k

             cell_map(i,j,k) = icell
          enddo
       enddo bz2
    enddo

    if (icell /= ntot2) then
       write(*,*) "Something went wrong in the cell mapping"
       write(*,*) "I am missing some virtual cells"
       write(*,*) icell, ntot2
       write(*,*)
       stop
    endif

    !if (cell_map(1,1,1) /= 1) then
    !   write(*,*) "WARNING : mapping of cell (1,1,1) is not 1"
    !   write(*,*) "(1,1,1) --->", cell_map(1,1,1)
    !   write(*,*) "MCFOST might crash"
    !   !write(*,*) "Exiting"
    !   !stop
    !endif

    return

  end subroutine build_cylindrical_cell_mapping

  !******************************************************************************

  pure logical function exit_test_cylindrical(icell, x, y, z)

    integer, intent(in) :: icell
    real(kind=db), intent(in) :: x,y,z

    if (icell <= n_cells) then
       exit_test_cylindrical = .false.
       return
    endif

    if (lexit_cell(icell)==0) then
       exit_test_cylindrical = .false.
    else if (lexit_cell(icell)==1) then ! radial
       exit_test_cylindrical = .true.
    else ! 2 --> vertical
       if (abs(z) > zmaxmax) then
          exit_test_cylindrical = .true.
       else
          exit_test_cylindrical = .false.
       endif
    endif

    return

  end function exit_test_cylindrical

  !******************************************************************************


  subroutine test_convert()

    integer :: i, j, k, icell
    integer :: i2,j2,k2


    write(*,*)
    write(*,*) "TEST CONVERT"

    do k=1, n_az
       do j=1, nz+1
          do i=0, n_rad

             icell = cell_map(i,j,k)
             write(*,*) "convert", i,j,k, "-->", icell

             call cell2cylindrical(icell, i2,j2,k2)
             if (i>0) then
                if ((i/=i2).or.(j/=j2).or.(k2/=k)) then
                   write(*,*) "PB test convert"
                   write(*,*) i,j,k, "-->", icell
                   write(*,*) icell, "-->", i2,j2,k2
                   stop
                endif
             else
                if ((i/=i2)) then ! seul i est defini ds la cas 0
                   write(*,*) "PB test convert"
                   write(*,*) i,j,k, "-->", icell
                   write(*,*) icell, "-->", i2,j2,k2
                   stop
                endif
             endif
          enddo
       enddo
    enddo

    write(*,*) "DONE"
    stop
    return


  end subroutine test_convert

  !******************************************************************************
  !pure subroutine cylindrical2cell(i,j,k, icell)
  !
  !  integer, intent(in) :: i,j,k
  !  integer, intent(out) :: icell
  !
  !  icell = cell_map(i,j,k)
  !
  !  return
  !
  !end subroutine cylindrical2cell

  !******************************************************************************

  pure subroutine cell2cylindrical(icell, i,j,k) ! work also in sph

    integer, intent(in) :: icell
    integer, intent(out) :: i,j,k

    i = cell_map_i(icell)
    j = cell_map_j(icell)
    k = cell_map_k(icell)

    return

  end subroutine cell2cylindrical

  !******************************************************************************

  subroutine cylindrical2cell_old(i,j,k, icell)
    ! icell is between 1 and n_rad * (n_z+1) * n_az

    integer, intent(in) :: i,j,k
    integer, intent(out) :: icell

    if ((i==0).and.(j==0)) then
       icell = 0
    else if (j>nz+1) then
       icell = -i
    else
       icell = i + n_rad * ( j-1 + nz * (k-1))
    endif

    return

  end subroutine cylindrical2cell_old

  !******************************************************************************

  subroutine cell2cylindrical_old(icell, i,j,k)

    integer, intent(in) :: icell
    integer, intent(out) :: i,j,k

    integer :: ij ! indice combine i et j, ie : i + (j-1) * n_rad

    if (icell==0) then
       i=0
       j=0
       k=1
    else if (icell < 0) then
       i = -icell
       j = nz+2
       k = 1
    else
       k = (icell-1)/nrz + 1 ; if (k > n_az) k=n_az

       ij = icell - (k-1)*nrz
       j = (ij-1)/n_rad + 1 ; if (j > nz+1) j=nz+1

       i = ij - (j-1)*n_rad

       !write(*,*) "TEST ij", ij,  ij/n_rad
       !write(*,*) "i,j", i, j
    endif

    return

  end subroutine cell2cylindrical_old

  !******************************************************************************

  subroutine indice_cellule(xin,yin,zin,ri_out,zj_out)

    implicit none

    real(kind=db), intent(in) :: xin,yin,zin
    integer, intent(out) :: ri_out, zj_out

    real(kind=db) :: r2
    integer :: ri, ri_min, ri_max

    r2 = xin*xin+yin*yin


    if (r2 < r_lim_2(0)) then
       ri_out=0
       zj_out=1
       return
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

    zj_out = floor(min(real(abs(zin)/zmax(ri_out) * nz),max_int))+1
    if (zj_out > nz) then
       zj_out = nz + 1
    endif

    return

  end subroutine indice_cellule

  !******************************************************************************

  subroutine indice_cellule_3D(xin,yin,zin,ri_out,zj_out,phik_out)

    implicit none

    real(kind=db), intent(in) :: xin,yin,zin
    integer, intent(out) :: ri_out, zj_out, phik_out

    real(kind=db) :: r2, phi
    integer :: ri, ri_min, ri_max

    r2 = xin*xin+yin*yin

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

    if (ri_out > 0) then
       zj_out = floor(min(real(abs(zin)/zmax(ri_out) * nz),max_int))+1
    else
       zj_out = 0
    endif
    if (zj_out > nz) zj_out = nz
    if (zin < 0.0)  zj_out = -zj_out

    if (zin /= 0.0) then
       phi=modulo(atan2(yin,xin),2*real(pi,kind=db))
       phik_out=floor(phi/(2*pi)*real(N_az))+1
       if (phik_out==n_az+1) phik_out=n_az
    else
       phik_out=1
    endif

    return

  end subroutine indice_cellule_3D

  !******************************************************************************

  subroutine indice_cellule_3D_phi(xin,yin,zin,phik_out)

    implicit none

    real(kind=db), intent(in) :: xin,yin,zin
    integer, intent(out) :: phik_out

    real(kind=db) :: phi

    if (zin /= 0.0) then
       phi=modulo(atan2(yin,xin),2*real(pi,kind=db))
       phik_out=floor(phi/(2*pi)*real(N_az))+1
       if (phik_out==n_az+1) phik_out=n_az
    else
       phik_out=1
    endif

    return

  end subroutine indice_cellule_3D_phi

  !******************************************************************************

  subroutine cross_cylindrical_cell(lambda, x0,y0,z0, u,v,w,  cell, previous_cell, x1,y1,z1, next_cell, l)

    integer, intent(in) :: lambda, cell, previous_cell
    real(kind=db), intent(in) :: x0,y0,z0
    real(kind=db), intent(in) :: u,v,w ! Todo : check that

    real(kind=db), intent(out) :: x1, y1, z1
    integer, intent(out) :: next_cell
    real(kind=db), intent(out) :: l

    ! Variables to be sorted out
    integer :: ri0,zj0,k0, k0m1
    integer ::  delta_rad, delta_zj, delta_phi, ri1, zj1, k1

    real(kind=db) :: inv_a, a, b, c, s, rac, t, t_phi, delta, inv_w, r_2, den, tan_angle_lim
    real(kind=db) :: phi, delta_vol, zlim, dotprod
    real(kind=db) :: correct_moins, correct_plus


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
       ! Si on est avant le bord interne,  on passe forcement par rmin
       ! et on cherche forcement la racine positive (unique)
       c=(r_2-r_lim_2(0))*inv_a
       delta=b*b-c
       rac=sqrt(delta)
       s = (-b+rac) * correct_plus
       t=huge_real
       t_phi= huge_real
       delta_rad=1
    else
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
             ! on s'eloigne du midplane (ou on monte en 2D)
             if (abs(zj0)==nz+1) then
                delta_zj=0
                zlim=sign(1.0e10_db,z0)
             else
                zlim= sign(z_lim(ri0,abs(zj0)+1)*correct_plus, z0)  ! BUUG HERE TODO
                delta_zj=1
                if (l3D.and.(z0 < 0.0))  delta_zj=-1
             endif
          else
             !  on se rapproche du midplane (ou on descend en 2D)
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
             endif ! 3D
          endif ! monte ou descend
          t=(zlim-z0)*inv_w
          ! correct pb precision
          if (t < 0.0_db) t=prec_grille
       endif !dotprod=0.0


       ! 3) position interface azimuthale
       if (l3D) then
          dotprod =  x0*v - y0*u
          if (abs(dotprod) < 1.0e-10) then
             ! on ne franchit pas d'interface azimuthale
             t_phi = 1.0e30
          else
             ! Quelle cellule on va franchir
             if (dotprod > 0.0) then
                tan_angle_lim = tan_phi_lim(k0)
                delta_phi=1
             else
                k0m1=k0-1
                if (k0m1==0) k0m1=N_az
                tan_angle_lim = tan_phi_lim(k0m1)
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
       else ! l3D
          t_phi = huge_real
       endif
    endif ! ri0==0


    ! 4) interface en r ou z ?
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
          if (l3D.and.(z1 < 0.0)) zj1=-zj1
       endif

       k1=k0
       if (l3D) then
          ! We need to find the azimuth when we enter the disc
          ! It can be different from the initial azimuth if the star is not centered
          ! so we need to compute it here
          if (ri0==0) then
             phi=modulo(atan2(y1,x1),2*real(pi,kind=db))
             k1=floor(phi*un_sur_deux_pi*real(N_az))+1
             if (k1==n_az+1) k1=n_az
          endif
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
       k1=k0
    else ! phi --> only happens in 3D
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
       k1=k0+delta_phi
       if (k1 == 0) k1=N_az
       if (k1 == N_az+1) k1=1
    endif

    ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
    if (z1 == 0.0_db) then
       if (l3D) then
          z1 = sign(prec_grille,w)
       else
          z1 = prec_grille
       endif
    endif

    !call cylindrical2cell(ri1,zj1,1, next_cell)
    next_cell = cell_map(ri1,zj1,k1)

    return

  end subroutine cross_cylindrical_cell

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

end module cylindrical_grid
