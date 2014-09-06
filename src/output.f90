module output

  use parametres
  use resultats
  use disk
  use grains
  use em_th
  use constantes
  use opacity
  use molecular_emission
  use ray_tracing
  use utils

  implicit none

  contains

subroutine capteur(id,lambda,ri0,zj0,xin,yin,zin,uin,vin,win,stokin,flag_star,flag_scatt)

  implicit none

  real(kind=db), intent(in) ::  xin,yin,zin, uin,vin,win
  real(kind=db) ::  x1,y1,z1,u1,v1,w1
  integer, intent(in) :: id, lambda, ri0, zj0
  real(kind=db), dimension(4), intent(in)  :: stokin
  logical, intent(in) :: flag_star, flag_scatt
  real(kind=db), dimension(4)  :: stok

  integer :: capt, c_phi, imap1, jmap1, imap2, jmap2, i
  real(kind=db) :: xprim, yprim, zprim, ytmp, ztmp

  x1=xin ; y1=yin ; z1=zin
  u1=uin ; v1=vin ; w1=win
  stok=stokin


  !* ------------ SECTION CAPTEURS------------------
  !*
  !*
  !*     utilisation de la symetrie N-S, selon l'axe Z
  !*     changement de Z1 -> -Z1 si W1<0
  !*             et de W1 -> -W1 si W1<0
  !*
  !        if (W1 < 0.0) then
  !           W1 = -W1
  !           Z1 = -Z1
  !           STOK(3,1) = -STOK(3,1)
  !        endif

  !*     utilisation de la symetrie centrale
  if (w1 < 0.0) then
     if (l_sym_centrale) then
        x1 = -x1
        y1 = -y1
        z1 = -z1
        u1 = -u1
        v1 = -v1
        w1 = -w1
        stok(3) = -stok(3)
     else
        return
     endif
  endif


  !* Selection angle theta
  CAPT=int(( -1.0*W1 + 1.0 ) * N_thet) + 1
  if (CAPT == (N_thet+1)) then
     CAPT = N_thet
  endif

  if (lonly_capt_interet) then
     if ((capt > capt_sup).or.(capt < capt_inf)) return
  endif


  if (l_sym_axiale) then ! symetrie axiale du systeme
     ! Utilisation de la symetrie Est-Ouest
     ! (le centre du mur est dans la direction u)
     if (v1 < 0.0) then
        v1 = - v1
        y1 = - y1
        stok(3) = -stok(3)
     endif
     ! Selection angle phi
     ! (atan2 renvoie dans [-Pi,Pi] mais ici dans [0,Pi] car v1 > 0)
     if (w1==1.0) then
        c_phi=1
     else
        c_phi = int(atan2(v1,u1)/pi*N_phi) + 1
     endif
  else !l_sym_axiale
     ! Pas de sym : on regarde entre 0 et 2Pi
     if (w1==1.0) then
        c_phi=1
     else
        c_phi = int(modulo(atan2(U1,V1)+pi/2,2*real(pi,kind=db))/(2*pi)*N_phi) + 1
     endif
  endif !l_sym_axiale

  if (c_phi == (N_phi + 1)) then
     c_phi = N_phi
  elseif (c_phi==0) then
     c_phi=1
  endif

  ! Origine du paquet
  if (lorigine) then
     if (capt == capt_interet) then
        if (ri0 == 0) then
           star_origin(lambda,id) = star_origin(lambda,id) + stok(1)
        else
           disk_origin(lambda,ri0,zj0,id) = disk_origin(lambda, ri0,zj0, id) + stok(1)
        endif
     endif
  endif


  if (lmono0) then ! Creation carte
     !*****************************************************
     !*----DETERMINATION DE LA POSITION SUR LA CARTE
     !*----IL FAUT FAIRE UNE ROTATION DU POINT
     !*    (X1,Y1,Z1) POUR RAMENER LES COORDONNEES DANS
     !*    LE SYSTEME OU (U1,V1,W1)=(1,0,0)
     !*
     !*    utilisation de la symetrie gauche-droite pour
     !*    produire une demi-carte seulement, dans les
     !*    "yprim positifs"
     !*****************************************************

     call ROTATION(X1,Y1,Z1,u1,v1,w1,XPRIM,YPRIM,ZPRIM)

     !*     LA CARTE IGRIDX * IGRIDY EST PRODUITE DANS LE PLAN Y'Z'
     !*
     !*     utilisation de la symetrie gauche-droite, selon l'axe Y'
     !*     changement de yprim -> -yprim si yprim<0
     !*
     !        if (YPRIM < 0.0) then
     !           YPRIM = -YPRIM
     !           STOK(3,1) = -STOK(3,1)
     !        endif

     !        IMAP = int(((YPRIM / size_neb) * IGRIDY / 2.) + 0.5) + 1
     !        if (IMAP <= 0 ) IMAP = 1
     !        if (IMAP >= (IGRIDX + 1)) IMAP = IGRIDX

     ! Rotation eventuelle du disque
     ytmp = yprim
     ztmp = zprim
     yprim = ytmp * cos_disk + ztmp * sin_disk
     zprim = ztmp * cos_disk - ytmp * sin_disk

     IMAP1 = int((YPRIM*zoom + 0.5*map_size)*size_pix) + deltapix_x
     if (IMAP1 <= 0 ) return !cycle photon
     if (IMAP1 > IGRIDX)  return !cycle photon

     JMAP1 = int((ZPRIM*zoom + 0.5*map_size)*size_pix)  + deltapix_y
     if (JMAP1 <= 0) return !cycle photon
     if (JMAP1 > IGRIDY) return !cycle photon

     if (l_sym_ima) then
        ! 1/2 photon
        STOKEI(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)

        if (lsepar_pola) then
           STOKEQ(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEQ(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(2)
           STOKEU(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEU(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(3)
           STOKEV(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEV(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(4)
        endif

        if (lsepar_contrib) then
           if (flag_star) then ! photon étoile
              if (flag_scatt) then
                 STOKEI_star_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_star_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
              else
                 STOKEI_star(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_star(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
              endif
           else ! photon thermique
              if (flag_scatt) then
                 STOKEI_disk_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_disk_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
              else
                 STOKEI_disk(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_disk(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
              endif
           endif ! type de photon
        endif !lsepar

        ! 1/2 photon symetrique
        ytmp = - ytmp
        yprim = ytmp * cos_disk - ztmp * sin_disk
        zprim = ztmp * cos_disk + ytmp * sin_disk

        IMAP2 = int((YPRIM*zoom + 0.5*map_size)*size_pix)  + deltapix_x
        if (IMAP2 <= 0 ) return !cycle photon
        if (IMAP2 > IGRIDX)  return !cycle photon


        JMAP2 = int((ZPRIM*zoom + 0.5*map_size)*size_pix)  + deltapix_y
        if (JMAP2 <= 0) return !cycle photon
        if (JMAP2 > IGRIDY) return !cycle photon


        if ((IMAP1==IMAP2).and.(JMAP1==JMAP2)) then ! Pas de sym on est dans le meme pixel
           ! on rajoute la 2eme moitie du photon
           STOKEI(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)

           if (lsepar_pola) then
              STOKEQ(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEQ(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(2)
              STOKEU(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEU(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(3)
              STOKEV(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEV(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(4)
           endif

           if (lsepar_contrib) then
              if (flag_star) then ! photon étoile
                 if (flag_scatt) then
                    STOKEI_star_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                         STOKEI_star_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
                 else
                    STOKEI_star(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                         STOKEI_star(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
                 endif
              else ! photon thermique
                 if (flag_scatt) then
                    STOKEI_disk_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                         STOKEI_disk_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
                 else
                    STOKEI_disk(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                         STOKEI_disk(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
                 endif
              endif ! type de photon
           endif ! lsepar
        else ! symetrie
           ! on rajoute la 2eme moitie du photon dans le pix oppose avec prop symetrie du vecteur de Stokes
           STOKEI(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEI(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)

           if (lsepar_pola) then
              STOKEQ(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEQ(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(2)
              STOKEU(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEU(lambda,IMAP2,JMAP2,capt,c_phi,id) - 0.5 * STOK(3)
              STOKEV(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEV(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(4)
           endif

           if (lsepar_contrib) then
              if (flag_star) then ! photon étoile
                 if (flag_scatt) then
                    STOKEI_star_scat(lambda,IMAP2,JMAP2,capt,c_phi,id) = &
                         STOKEI_star_scat(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
                 else
                    STOKEI_star(lambda,IMAP2,JMAP2,capt,c_phi,id) = &
                         STOKEI_star(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
                 endif
              else ! photon thermique
                 if (flag_scatt) then
                    STOKEI_disk_scat(lambda,IMAP2,JMAP2,capt,c_phi,id) = &
                         STOKEI_disk_scat(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
                 else
                    STOKEI_disk(lambda,IMAP2,JMAP2,capt,c_phi,id) = &
                         STOKEI_disk(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
                 endif
              endif ! type de photon
           endif !lsepar
        endif
     else
        ! Pas de symetrie
        STOKEI(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)

        if (lsepar_pola) then
           STOKEQ(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEQ(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(2)
           STOKEU(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEU(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(3)
           STOKEV(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEV(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(4)
        endif
        if (lsepar_contrib) then
           if (flag_star) then ! photon étoile
              if (flag_scatt) then
                 STOKEI_star_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_star_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
              else
                 STOKEI_star(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_star(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
              endif
           else ! photon thermique
              if (flag_scatt) then
                 STOKEI_disk_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_disk_scat(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
              else
                 STOKEI_disk(lambda,IMAP1,JMAP1,capt,c_phi,id) = &
                      STOKEI_disk(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
              endif
           endif ! type de photon
        endif !lsepar
     endif

!!!!!!!!!!!!!!
     !Warning : marche que sans rotation !!!
     if (n_cartes > 1) then
        IMAP1 = int((YPRIM*zoom + 0.5*map_size)*size_pix2(1)) + deltapix_x2(1)
        if (IMAP1 <= 0 ) return !cycle photon
        if (IMAP1 > IGRIDX2(1))  return !cycle photon


        JMAP1 = int((ZPRIM*zoom + 0.5*map_size)*size_pix2(1))  + deltapix_y2(1)
        if (JMAP1 <= 0) return !cycle photon
        if (JMAP1 > IGRIDY2(1)) return !cycle photon

        if (l_sym_ima) then
           ! 1/2 photon
           STOKEI1(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI1(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)

           ! 1/2 photon symetrique
          ! ytmp = - ytmp
          ! yprim = ytmp * cos_disk - ztmp * sin_disk
          ! zprim = ztmp * cos_disk + ytmp * sin_disk

           IMAP2 = int((YPRIM*zoom + 0.5*map_size)*size_pix2(1))  + deltapix_x2(1)
           if (IMAP2 <= 0 ) return !cycle photon
           if (IMAP2 > IGRIDX2(1))  return !cycle photon


           JMAP2 = int((ZPRIM*zoom + 0.5*map_size)*size_pix2(1))  + deltapix_y2(1)
           if (JMAP2 <= 0) return !cycle photon
           if (JMAP2 > IGRIDY2(1)) return !cycle photon

           if ((IMAP1==IMAP2).and.(JMAP1==JMAP2)) then ! Pas de sym on est dans le meme pixel
              ! on rajoute la 2eme moitie du photon
              STOKEI1(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI1(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
           else ! symetrie
              ! on rajoute la 2eme moitie du photon dansle pix oppose avec prop symetrie du vecteur de Stokes
              STOKEI1(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEI1(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
           endif
        else
           ! Pas de symetrie
           STOKEI1(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI1(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
        endif
!!!
     endif
     if (n_cartes > 2) then
        IMAP1 = int((YPRIM*zoom + 0.5*map_size)*size_pix2(2)) + deltapix_x2(2)
        if (IMAP1 <= 0 ) return !cycle photon
        if (IMAP1 > IGRIDX2(2))  return !cycle photon


        JMAP1 = int((ZPRIM*zoom + 0.5*map_size)*size_pix2(2))  + deltapix_y2(2)
        if (JMAP1 <= 0) return !cycle photon
        if (JMAP1 > IGRIDY2(2)) return !cycle photon

        if (l_sym_ima) then
           ! 1/2 photon
           STOKEI2(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI2(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)

           ! 1/2 photon symetrique
          ! ytmp = - ytmp
          ! yprim = ytmp * cos_disk - ztmp * sin_disk
          ! zprim = ztmp * cos_disk + ytmp * sin_disk

           IMAP2 = int((YPRIM*zoom + 0.5*map_size)*size_pix2(2))  + deltapix_x2(2)
           if (IMAP2 <= 0 ) return !cycle photon
           if (IMAP2 > IGRIDX2(2))  return !cycle photon


           JMAP2 = int((ZPRIM*zoom + 0.5*map_size)*size_pix2(2))  + deltapix_y2(2)
           if (JMAP2 <= 0) return !cycle photon
           if (JMAP2 > IGRIDY2(2)) return !cycle photon

           if ((IMAP1==IMAP2).and.(JMAP1==JMAP2)) then ! Pas de sym on est dans le meme pixel
              ! on rajoute la 2eme moitie du photon
              STOKEI2(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI2(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
           else ! symetrie
              ! on rajoute la 2eme moitie du photon dansle pix oppose avec prop symetrie du vecteur de Stokes
              STOKEI2(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEI2(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
           endif
        else
           ! Pas de symetrie
           STOKEI2(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI2(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
        endif

     endif
     if (n_cartes > 3) then
        IMAP1 = int((YPRIM*zoom + 0.5*map_size)*size_pix2(3)) + deltapix_x2(3)
        if (IMAP1 <= 0 ) return !cycle photon
        if (IMAP1 > IGRIDX2(3))  return !cycle photon


        JMAP1 = int((ZPRIM*zoom + 0.5*map_size)*size_pix2(3))  + deltapix_y2(3)
        if (JMAP1 <= 0) return !cycle photon
        if (JMAP1 > IGRIDY2(3)) return !cycle photon

        if (l_sym_ima) then
           ! 1/2 photon
           STOKEI3(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI3(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)

           ! 1/2 photon symetrique
          ! ytmp = - ytmp
          ! yprim = ytmp * cos_disk - ztmp * sin_disk
          ! zprim = ztmp * cos_disk + ytmp * sin_disk

           IMAP2 = int((YPRIM*zoom + 0.5*map_size)*size_pix2(3))  + deltapix_x2(3)
           if (IMAP2 <= 0 ) return !cycle photon
           if (IMAP2 > IGRIDX2(3))  return !cycle photon


           JMAP2 = int((ZPRIM*zoom + 0.5*map_size)*size_pix2(3))  + deltapix_y2(3)
           if (JMAP2 <= 0) return !cycle photon
           if (JMAP2 > IGRIDY2(3)) return !cycle photon

           if ((IMAP1==IMAP2).and.(JMAP1==JMAP2)) then ! Pas de sym on est dans le meme pixel
              ! on rajoute la 2eme moitie du photon
              STOKEI3(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI3(lambda,IMAP1,JMAP1,capt,c_phi,id) + 0.5 * STOK(1)
           else ! symetrie
              ! on rajoute la 2eme moitie du photon dansle pix oppose avec prop symetrie du vecteur de Stokes
              STOKEI3(lambda,IMAP2,JMAP2,capt,c_phi,id) = STOKEI3(lambda,IMAP2,JMAP2,capt,c_phi,id) + 0.5 * STOK(1)
           endif
        else
           ! Pas de symetrie
           STOKEI3(lambda,IMAP1,JMAP1,capt,c_phi,id) = STOKEI3(lambda,IMAP1,JMAP1,capt,c_phi,id) + STOK(1)
        endif
     endif !n_cartes

!!!!!!!!!!!!!!
  else ! Creation sed
     sed(id,lambda,capt,c_phi) = sed(id,lambda,capt,c_phi) + stok(1)
     sed_q(id,lambda,capt,c_phi) = sed_q(id,lambda,capt,c_phi) + stok(2)
     sed_u(id,lambda,capt,c_phi) = sed_u(id,lambda,capt,c_phi) + stok(3)
     sed_v(id,lambda,capt,c_phi) = sed_v(id,lambda,capt,c_phi) + stok(4)
     n_phot_sed(id,lambda,capt,c_phi) = n_phot_sed(id,lambda,capt,c_phi) + 1.0_db
     n_phot_sed2(id,lambda,capt,c_phi) = n_phot_sed2(id,lambda,capt,c_phi) + 1.0_db
     if (flag_star) then ! photon étoile
        if (flag_scatt) then
           sed_star_scat(id,lambda,capt,c_phi) = sed_star_scat(id,lambda,capt,c_phi) + stok(1)
        else
           sed_star(id,lambda,capt,c_phi) = sed_star(id,lambda,capt,c_phi) + stok(1)
        endif
     else ! photon thermique
        if (flag_scatt) then
           sed_disk_scat(id,lambda,capt,c_phi) = sed_disk_scat(id,lambda,capt,c_phi) + stok(1)
        else
           sed_disk(id,lambda,capt,c_phi) = sed_disk(id,lambda,capt,c_phi) + stok(1)
        endif
     endif ! type de photon
  endif

  return

end subroutine capteur

!**********************************************************************

subroutine find_pixel(x,y,z,u,v,w, i, j, in_map)

  real(kind=db), intent(in) :: x,y,z,u,v,w
  integer, intent(out) :: i,j
  logical, intent(out) :: in_map

  real(kind=db) :: x2,y2,z2, y_map,z_map

  !*****************************************************
  !*----DETERMINATION DE LA POSITION SUR LA CARTE
  !*----IL FAUT FAIRE UNE ROTATION DU POINT
  !*    (X1,Y1,Z1) POUR RAMENER LES COORDONNEES DANS
  !*    LE SYSTEME OU (U1,V1,W1)=(1,0,0)
  !*****************************************************

  call rotation(x,y,z, u,v,w, x2,y2,z2)

  ! rotation eventuelle du disque
  y_map = y2 * cos_disk + z2 * sin_disk
  z_map = z2 * cos_disk - y2 * sin_disk

  i = int((y_map*zoom + 0.5*map_size)*size_pix) + deltapix_x
  j = int((z_map*zoom + 0.5*map_size)*size_pix) + deltapix_y

  if ((i<1).or.(i>igridx).or.(j<1).or.(j>igridy)) then
     in_map = .false.
  else
     in_map = .true.
  endif

  return

end subroutine find_pixel

!**********************************************************************

subroutine write_stokes_fits()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, alloc_status, id

  integer :: lambda=1
  integer :: image_type=1

  character(len = 512) :: filename
  logical :: simple, extend

  real, dimension(n_lambda) :: n_photons_envoyes
  real :: facteur, pixel_scale_x, pixel_scale_y

  real, dimension(:,:,:,:), allocatable :: tmp

  real :: o_star, frac_star, somme_disk
  real, dimension(n_rad, nz) :: o_disk

  if (lorigine) then

     o_star = sum(star_origin(1,:))
     o_disk(:,:) = sum(disk_origin(1,:,:,:), dim=3)
     somme_disk = sum(o_disk)
     frac_star = o_star / (o_star+somme_disk)
     if (somme_disk > tiny_real) o_disk=o_disk/somme_disk


     filename = trim(data_dir)//"/origin.fits.gz"
      !  Get an unused Logical Unit Number to use to open the FITS file.
     status=0
     call ftgiou (unit,status)

      !  Create the new empty FITS file.
     blocksize=1
     call ftinit(unit,trim(filename),blocksize,status)

      !  Initialize parameters about the FITS image
     simple=.true.
      ! le signe - signifie que l'on ecrit des reels dans le fits
     bitpix=-32
     extend=.true.

     naxis=2
     naxes(1)=n_rad
     naxes(2)=nz

     ! Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)

      ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,o_disk,status)

      !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

      !  Check for any error, and if so print out error messages
     if (status > 0) then
	call print_error(status)
     end if
  endif



  filename = trim(data_dir)//"/MC.fits.gz"

  stoke_io(:,:,:,:,:)=0.0
  do id=1, nb_proc
     stoke_io(:,:,:,:,1)=stoke_io(:,:,:,:,1)+stokei(lambda,:,:,:,:,id)
  enddo
  deallocate(stokei)

  if (lsepar_pola) then
     do id=1, nb_proc
        stoke_io(:,:,:,:,2)=stoke_io(:,:,:,:,2)+stokeq(lambda,:,:,:,:,id)
        stoke_io(:,:,:,:,3)=stoke_io(:,:,:,:,3)+stokeu(lambda,:,:,:,:,id)
        stoke_io(:,:,:,:,4)=stoke_io(:,:,:,:,4)+stokev(lambda,:,:,:,:,id)
     enddo
     deallocate(stokeq,stokeu,stokev)
  endif

  if (lsepar_contrib) then
     do id=1, nb_proc
        stoke_io(:,:,:,:,n_Stokes+1)=stoke_io(:,:,:,:,n_Stokes+1)+stokei_star(lambda,:,:,:,:,id)
        stoke_io(:,:,:,:,n_Stokes+2)=stoke_io(:,:,:,:,n_Stokes+2)+stokei_star_scat(lambda,:,:,:,:,id)
        stoke_io(:,:,:,:,n_Stokes+3)=stoke_io(:,:,:,:,n_Stokes+3)+stokei_disk(lambda,:,:,:,:,id)
        stoke_io(:,:,:,:,n_Stokes+4)=stoke_io(:,:,:,:,n_Stokes+4)+stokei_disk_scat(lambda,:,:,:,:,id)
     enddo
     deallocate(stokei_star,stokei_star_scat,stokei_disk,stokei_disk_scat)
  endif

  status=0


  lambda=1
  n_photons_envoyes(lambda) = real(sum(n_phot_envoyes(lambda,:)))
  E_totale(lambda) = E_totale(lambda)/n_photons_envoyes(lambda)
  facteur = E_totale(lambda) * tab_lambda(lambda) * 1.0e-6
  stoke_io = stoke_io * facteur

  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou ( unit, status )

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  naxis=5
  naxes(1)=igridx
  naxes(2)=igridy
  naxes(3)=N_incl
  naxes(4)=N_phi
  naxes(5)=N_type_flux

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)


  ! Write  optional keywords to the header
  !  wavelength
  call ftpkye(unit,'WAVE',tab_lambda(lambda),-3,'wavelength [microns]',status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-3,'RAD',status)
  call ftpkyj(unit,'CRPIX1',igridx/2+1,'',status)
  pixel_scale_x = -map_size / (igridx * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-3,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-3,'DEC',status)
  call ftpkyj(unit,'CRPIX2',igridy/2+1,'',status)
  pixel_scale_y = map_size / (igridy * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-3,'pixel scale y [deg]',status)
  call ftpkys(unit,'BUNIT',"W.m-2.pixel-1",' ',status)

  call ftpkys(unit,'FLUX_1',"I = total flux",' ',status)
  if (lsepar_pola) then
     call ftpkys(unit,'FLUX_2',"Q",' ',status)
     call ftpkys(unit,'FLUX_3',"U",' ',status)
     call ftpkys(unit,'FLUX_4',"V",' ',status)
  endif
  if (lsepar_contrib) then
     if (lsepar_pola) then
        call ftpkys(unit,'FLUX_5',"direct star light",' ',status)
        call ftpkys(unit,'FLUX_6',"scattered star light",' ',status)
        call ftpkys(unit,'FLUX_7',"direct thermal emission",' ',status)
        call ftpkys(unit,'FLUX_8',"scattered thermal emssion",' ',status)
     else
        call ftpkys(unit,'FLUX_2',"direct star light",' ',status)
        call ftpkys(unit,'FLUX_3',"scattered star light",' ',status)
        call ftpkys(unit,'FLUX_4',"direct thermal emission",' ',status)
        call ftpkys(unit,'FLUX_5',"scattered thermal emssion",' ',status)
     endif
  endif

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,stoke_io,status)

  ! Cartes de differentes resolutions
  if (n_cartes > 1) then
     ! create new hdu
     naxis=4
     naxes(1)=igridx2(1)
     naxes(2)=igridy2(1)

     call FTCRHD(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
     ! le e signifie real*4
     allocate(tmp(igridx2(1),igridy2(1),N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tmp 1'
        stop
     endif
     tmp = real(sum(stokei1(lambda,:,:,:,:,:),dim=5))
     call ftppre(unit,group,fpixel,nelements,tmp,status)
     deallocate(tmp)

  endif
  if (n_cartes > 2) then
     ! create new hdu
     call FTCRHD(unit, status)
     naxis=4
     naxes(1)=igridx2(2)
     naxes(2)=igridy2(2)


     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
     ! le e signifie real*4
     allocate(tmp(igridx2(2),igridy2(2),N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tmp 2'
        stop
     endif
     tmp = real(sum(stokei2(lambda,:,:,:,:,:),dim=5))
     call ftppre(unit,group,fpixel,nelements,tmp,status)
     deallocate(tmp)

  endif
  if (n_cartes > 3) then
     ! create new hdu
     call FTCRHD(unit, status)
     naxis=4
     naxes(1)=igridx2(3)
     naxes(2)=igridy2(3)


     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
     ! le e signifie real*4
     allocate(tmp(igridx2(3),igridy2(3),N_thet,N_phi), stat=alloc_status)
     if (alloc_status > 0) then
        write(*,*) 'Allocation error tmp 3'
        stop
     endif
     tmp = real(sum(stokei3(lambda,:,:,:,:,:),dim=5))
     call ftppre(unit,group,fpixel,nelements,tmp,status)
     deallocate(tmp)

  endif


  ! Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  ! Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return
end subroutine write_stokes_fits

!***********************************************************

subroutine ecriture_map_ray_tracing()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, alloc_status, id, xcenter, lambda, itype, ibin

  character(len = 512) :: filename
  logical :: simple, extend
  real :: pixel_scale_x, pixel_scale_y

  ! Allocation dynamique pour passer en stack
  real, dimension(:,:,:,:,:), allocatable :: image

  allocate(image(igridx, igridy, RT_n_ibin, 1, N_type_flux), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error RT image'
     stop
  endif
  image = 0.0 ;

  lambda=1
  filename = trim(data_dir)//"/RT.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit, status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  group=1
  fpixel=1
  extend=.true.
  bitpix=-32

  naxis=5
  naxes(1)=igridx
  naxes(2)=igridy
  naxes(3)= RT_n_ibin
  naxes(4)=1 ! N_phi
  naxes(5)=N_type_flux
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  !  wavelength
  call ftpkye(unit,'WAVE',tab_lambda(lambda),-3,'wavelength [microns]',status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-3,'RAD',status)
  call ftpkyj(unit,'CRPIX1',igridx/2+1,'',status)
  pixel_scale_x = -map_size / (igridx * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-3,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-3,'DEC',status)
  call ftpkyj(unit,'CRPIX2',igridy/2+1,'',status)
  pixel_scale_y = map_size / (igridy * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-3,'pixel scale y [deg]',status)
  call ftpkys(unit,'BUNIT',"W.m-2.pixel-1",' ',status)

  call ftpkys(unit,'FLUX_1',"I = total flux",' ',status)
  if (lsepar_pola) then
     call ftpkys(unit,'FLUX_2',"Q",' ',status)
     call ftpkys(unit,'FLUX_3',"U",' ',status)
     call ftpkys(unit,'FLUX_4',"V",' ',status)
  endif
  if (lsepar_contrib) then
     if (lsepar_pola) then
        call ftpkys(unit,'FLUX_5',"direct star light",' ',status)
        call ftpkys(unit,'FLUX_6',"scattered star light",' ',status)
        call ftpkys(unit,'FLUX_7',"direct thermal emission",' ',status)
        call ftpkys(unit,'FLUX_8',"scattered thermal emssion",' ',status)
     else
        call ftpkys(unit,'FLUX_2',"direct star light",' ',status)
        call ftpkys(unit,'FLUX_3',"scattered star light",' ',status)
        call ftpkys(unit,'FLUX_4',"direct thermal emission",' ',status)
        call ftpkys(unit,'FLUX_5',"scattered thermal emssion",' ',status)
     endif
  endif

  !----- Images
  ! Boucles car ca ne passe pas avec sum directement (ifort sur mac)
  do itype=1,N_type_flux
     do ibin=1,RT_n_ibin
        do j=1,igridy
           do i=1,igridx
              image(i,j,ibin,1,itype) = sum(Stokes_ray_tracing(lambda,i,j,ibin,itype,:))
           enddo !i
        enddo !j
     enddo !ibin
  enddo ! itype

  if (l_sym_ima) then
     xcenter = igridx/2 + modulo(igridx,2)
     do i=xcenter+1,igridx
        image(i,:,:,1,1) = image(igridx-i+1,:,:,1,1)

        if (lsepar_pola) then
           image(i,:,:,1,2) = image(igridx-i+1,:,:,1,2)
           image(i,:,:,1,3) = - image(igridx-i+1,:,:,1,3)
           image(i,:,:,1,4) = image(igridx-i+1,:,:,1,4)
        endif

        if (lsepar_contrib) then
           image(i,:,:,1,n_Stokes+1) = image(igridx-i+1,:,:,1,n_Stokes+1)
           image(i,:,:,1,n_Stokes+2) = image(igridx-i+1,:,:,1,n_Stokes+2)
           image(i,:,:,1,n_Stokes+3) = image(igridx-i+1,:,:,1,n_Stokes+3)
           image(i,:,:,1,n_Stokes+4) = image(igridx-i+1,:,:,1,n_Stokes+4)
        endif
     enddo
  endif ! l_sym_image

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,image,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  deallocate(image)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return

end subroutine ecriture_map_ray_tracing

!**********************************************************************

subroutine ecriture_sed_ray_tracing()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, alloc_status, id, lambda

  character(len = 512) :: filename
  logical :: simple, extend

  real, dimension(n_lambda2, RT_n_ibin, 1, N_type_flux) :: sed_rt

  lambda=1

  filename = trim(data_dir)//"/sed_rt.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit, status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  naxis=4
  naxes(1)=n_lambda2
  naxes(2)= RT_n_ibin
  naxes(3)=1 ! N_phi
  naxes(4)=N_type_flux

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  sed_rt = 0.0_db
  if (RT_sed_method == 1) then
     sed_rt(:,:,1,:) = sum(Stokes_ray_tracing(:,1,1,:,:,:),dim=4)
  else
     do i=1,igridx
        do j=1,igridy
           sed_rt(:,:,1,:) = sed_rt(:,:,1,:) +sum(Stokes_ray_tracing(:,i,j,:,:,:),dim=4)
        enddo
     enddo
  endif

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,sed_rt,status)

  ! Second HDU avec longueur d'onde
  call FTCRHD(unit, status)
  bitpix=-32
  naxis=1
  naxes(1)=n_lambda
  nelements=naxes(1)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,tab_lambda,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return
end subroutine ecriture_sed_ray_tracing

!**********************************************************************

subroutine write_origin()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, alloc_status, id

  integer :: lambda=1
  integer :: image_type=1

  character(len = 512) :: filename
  logical :: simple, extend

  real :: o_star, frac_star, somme_disk
  real, dimension(n_rad, nz) :: o_disk

  filename = trim(data_dir)//"/origine.fits.gz"

  ! Normalisation
  o_disk(:,:) = o_disk(:,:) / (sum(o_disk) + o_star)

  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou ( unit, status )

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  naxis=2
  naxes(1)=n_rad
  naxes(2)=nz

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,o_disk,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  return

end subroutine write_origin

!**********************************************************************

subroutine print_error(status)
  ! PRINT_ERROR prints out the FITSIO error messages to the user.

  integer status
  character ( len = 30 ) errtext
  character ( len = 80 ) errmessage

  !  Check if status is OK (no error); if so, simply return.
  if (status <= 0) then
     return
  end if

  !  Get the text string which describes the error
  call ftgerr(status,errtext)
  print *,'FITSIO Error Status =',status,': ',errtext

  !  Read and print out all the error messages on the FITSIO stack
  call ftgmsg(errmessage)
  do while (errmessage .ne. ' ')
     print *,errmessage
     call ftgmsg(errmessage)
  end do

  return
end subroutine print_error

!***********************************************************

subroutine calc_optical_depth_map(lambda)

  implicit none

  integer :: lambda

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,k, group,fpixel,nelements, alloc_status, id

  character(len = 512) :: filename
  logical :: simple, extend, lmilieu

  real, dimension(n_rad,nz,n_az,2) :: optical_depth_map

  lmilieu = .true. ! opacite au milieu ou a la "fin" de la cellule

  if (lmilieu) then
     ! Opacite radiale
     do k=1, n_az
        do j=1, nz
           i=1 ; optical_depth_map(i,j,k,1) = kappa(lambda,i,j,k)* 0.5 * (r_lim(i)-r_lim(i-1))
           do i=2, n_rad
              optical_depth_map(i,j,k,1) = optical_depth_map(i-1,j,k,1) + 0.5 * kappa(lambda,i-1,j,k)*(r_lim(i-1)-r_lim(i-2)) &
                   + 0.5 * kappa(lambda,i,j,k)*(r_lim(i)-r_lim(i-1))
           enddo
        enddo
     enddo

     ! Opacite verticale
     do i=1, n_rad
        do k=1, n_az
           j=nz ; optical_depth_map(i,j,k,2) = kappa(lambda,i,j,k)* 0.5 * (z_lim(i,j+1)-z_lim(i,j))
           do j=nz-1,1,-1
              optical_depth_map(i,j,k,2) = optical_depth_map(i,j+1,k,2) + 0.5 * kappa(lambda,i,j+1,k)*(z_lim(i,j+2)-z_lim(i,j+1)) &
                   + 0.5 * kappa(lambda,i,j,k)*(z_lim(i,j+1)-z_lim(i,j))
           enddo
        enddo
     enddo

  else
     ! Opacite radiale
     do k=1, n_az
        do j=1, nz
           i=1 ; optical_depth_map(i,j,k,1) = kappa(lambda,i,j,k)*(r_lim(i)-r_lim(i-1))
           do i=2, n_rad
              optical_depth_map(i,j,k,1) = optical_depth_map(i-1,j,k,1) +kappa(lambda,i,j,k)*(r_lim(i)-r_lim(i-1))
           enddo
        enddo
     enddo

     ! Opacite verticale
     do i=1, n_rad
        do k=1, n_az
           j=nz ; optical_depth_map(i,j,k,2) = kappa(lambda,i,j,k)*(z_lim(i,j+1)-z_lim(i,j))
           do j=nz-1,1,-1
              optical_depth_map(i,j,k,2) = optical_depth_map(i,j+1,k,2) + kappa(lambda,i,j,k)*(z_lim(i,j+1)-z_lim(i,j))
           enddo
        enddo
     enddo
  endif

  write(*,*) "Writing optical_depth_map.fits.gz for wl=", tab_lambda(lambda), "microns"
  filename = "optical_depth_map.fits.gz"

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  naxis=4
  naxes(1)=n_rad
  naxes(2)=nz
  naxes(3)=n_az
  naxes(4)=2

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,optical_depth_map,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)

  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  write(*,*) "Done"
  stop
  return

end subroutine calc_optical_depth_map

!***********************************************************

subroutine reemission_stats()

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, alloc_status, id

  character(len = 512) :: filename
  logical :: simple, extend

  real, dimension(n_rad,nz) :: N_reemission

  N_reemission = sum(nbre_reemission(:,:,1,:),dim=3)

  filename = trim(data_dir)//"/reemission_stats.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou ( unit, status )

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  naxis=2
  naxes(1)=n_rad
  naxes(2)=nz

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,N_reemission,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  return

end subroutine reemission_stats

!********************************************************************************

subroutine write_disk_struct()
! Ecrit les table de densite du gaz en g/cm^3
! de la poussiere en g/cm^3 et en particules
! + coordonnees r et z en AU
! + les tailles et masse des grains
! C. Pinte
! 3/06/06

  implicit none

  integer :: i, j, k

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id

  logical :: simple, extend
  character(len=512) :: filename

  real, dimension(n_rad,nz,n_az) :: dens
  real(kind=db), dimension(:,:,:,:), allocatable :: dust_dens
  real, dimension(n_rad) :: vol
  real(kind=db), dimension(n_rad,nz,2) :: grid


  allocate(dust_dens(n_rad,nz,n_az,n_grains_tot), stat = alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error dust density table for fits file'
     stop
  endif


  ! ********************************************************************************
  filename = "data_disk/gas_density.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=2
  naxes(1)=n_rad
  naxes(2)=nz
  nelements=naxes(1)*naxes(2)

  if (l3D) then
     naxis=3
     naxes(3)=n_az
     nelements=naxes(1)*naxes(2)*naxes(3)
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"g.cm^-3",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  dens =  densite_gaz(:,1:nz,:) * masse_mol_gaz / m3_to_cm3 ! nH2/m**3 --> g/cm**3

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,dens,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/dust_particule_density.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-64
  extend=.true.

  naxis=4
  naxes(1) = n_rad
  naxes(2) = nz
  naxes(3) = n_az
  naxes(4) = n_grains_tot

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"part.m^-3 [per grain size bin N(a).da]",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

  !  dens =  densite_pouss
  ! le d signifie real*8
  dust_dens = densite_pouss(:,1:nz,:,:) * m3_to_cm3
  call ftpprd(unit,group,fpixel,nelements,dust_dens,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/dust_mass_density.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=3
  naxes(1) = n_rad
  naxes(2) = nz
  naxes(3) = n_az

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"g.cm^-3",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)

  !  dens =  densite_pouss
  ! le d signifie real*8
  dust_dens = densite_pouss(:,1:nz,:,:)
  dens = 0.0
  do k=1,n_az
     do j=1,nz
        do i=1,n_rad
           dens(i,j,k) = sum(densite_pouss(i,j,k,:) * M_grain(:)) ! M_grain en g
        enddo !i
     enddo !j
  enddo !k
  call ftppre(unit,group,fpixel,nelements,dens,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/grain_sizes.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=1
  naxes(1:1) = shape(r_grain)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"microns",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,r_grain,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/grain_sizes_min.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=1
  naxes(1:1) = shape(r_grain)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"microns",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,r_grain_min,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/grain_sizes_max.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=1
  naxes(1:1) = shape(r_grain)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"microns",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,r_grain_max,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/grain_masses.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=1
  naxes(1:1) = shape(r_grain)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"g",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,M_grain,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = "data_disk/volume.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=1
  naxes(1)=n_rad

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"AU^3",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! conversion en simple precision
  vol = volume !* (AU_to_cm)**3

 ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,vol,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if


  ! ********************************************************************************
  filename = "data_disk/grid.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-64
  extend=.true.

  naxis=3
  naxes(1)=n_rad
  naxes(2)=nz
  naxes(3)=2

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"AU",' ',status)
  call ftpkys(unit,'DIM_1',"cylindrical radius",' ',status)
  call ftpkys(unit,'DIM_2',"elevation above midplane",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)

   do i=1, n_rad
     grid(i,:,1) = sqrt(r_lim(i) * r_lim(i-1))
     do j=1,nz
        grid(i,j,2) = (real(j)-0.5)*delta_z(i)
        !write(*,*) i, j, grid(i,j,1), grid(i,j,2)
     enddo
  enddo

  ! le d signifie real*8
  call ftpprd(unit,group,fpixel,nelements,grid,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  if (lstop_after_init) then
     write(*,*) "Exiting"
     stop
  endif

  return

end subroutine write_disk_struct

!********************************************************************

subroutine ecriture_J()
! Ecrit la table du champ de radiation
! C. Pinte
! 26/09/07

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(3) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id, lambda, ri, zj

  logical :: simple, extend
  character(len=512) :: filename

  real(kind=db), dimension(n_lambda,n_rad,nz) :: J
  real, dimension(n_rad,nz,n_lambda) :: Jio

  filename = trim(data_dir)//"/J.fits.gz"

  ! 1/4pi est inclus dans n_phot_l_tot
  J(:,:,:) = (sum(xJ_abs,dim=4) + J0(:,:,:,1)) * n_phot_L_tot !* 4* pi

  do ri=1, n_rad
     J(:,ri,:) = J(:,ri,:) / volume(ri)
  enddo

  ! xJ_abs est par bin de lambda donc Delta_lambda.F_lambda
  ! Jio en W.m-2 (lambda.F_lambda)
  ! teste OK par rapport a fct bb de yorick
  do lambda=1, n_lambda
     J(lambda,:,:) = J(lambda,:,:) * tab_lambda(lambda) / tab_delta_lambda(lambda)
     ! Inversion de l'ordre des dimensions + passage en simple precision
     do ri=1, n_rad
        do zj=1,nz
           Jio(ri,zj,lambda) = J(lambda,ri,zj)
        enddo
     enddo
  enddo

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=3
  naxes(1)=n_rad
  naxes(2)=nz
  naxes(3)=n_lambda

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,Jio,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return

end subroutine ecriture_J

!********************************************************************

subroutine ecriture_UV_field()
! Ecrit la table du champ de radiation
! C. Pinte
! 26/09/07

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(3) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id, lambda, ri, zj, l

  logical :: simple, extend
  character(len=512) :: filename

  integer, parameter :: n=200

  real(kind=db), dimension(n_lambda,n_rad,nz) :: J
  real(kind=db), dimension(n_rad,nz) :: G
  real(kind=db), dimension(n) :: wl, J_interp
  real(kind=db), dimension(n_lambda) :: lamb
  real(kind=db) :: delta_wl

  real, dimension(n_rad,nz) :: Gio


  filename = trim(data_dir)//"/UV_field.fits.gz"

  ! 1/4pi est inclus dans n_phot_l_tot
  J(:,:,:) = (sum(xJ_abs,dim=4) + J0(:,:,:,1)) * n_phot_L_tot !* 4* pi

  do ri=1, n_rad
     J(:,ri,:) = J(:,ri,:) / volume(ri)
  enddo

  ! xJ_abs est par bin de lambda donc Delta_lambda.F_lambda
  ! J en W.m-2.m-1 (F_lambda)
  do lambda=1, n_lambda
     J(lambda,:,:) = J(lambda,:,:) / (tab_delta_lambda(lambda) * 1.0e-6)
  enddo

  lamb = tab_lambda

  wl(:) = span(0.0912,0.24,n)
  delta_wl = (wl(n) - wl(1))/(n-1.) * 1e-6

  do ri=1,n_rad
     do zj=1,nz
        do l=1,n
           J_interp(l) = interp( J(:,ri,zj),lamb(:),wl(l))
        enddo

        ! Le Petit et al 2006 page 19-20
        ! integration trapeze
        G(ri,zj) = (sum(J_interp(:)) - 0.5 * (J_interp(1) + J_interp(n)) ) * delta_wl  &
             * 4*pi/c_light / (5.6e-14 * erg_to_J * m_to_cm**3)

        Gio(ri,zj) = G(ri,zj) ! Teste OK par comparaison avec yorick
     enddo
  enddo

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=2
  naxes(1)=n_rad
  naxes(2)=nz

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,Gio,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return

end subroutine ecriture_UV_field


!********************************************************************

subroutine ecriture_temperature(iTemperature)

  implicit none


  integer, intent(in) :: iTemperature

  integer :: n,k, i, j, lambda, l
  real :: flux, nbre_photons, facteur


  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id

  logical :: simple, extend
  character(len=512) :: filename

  integer, dimension(n_rad,nz,grain_nRE_start:grain_nRE_end) :: tmp


  if (lRE_LTE) then
     if (iTemperature == 2) then
        filename = trim(data_dir)//"/Temperature2.fits.gz"
     else
        filename = trim(data_dir)//"/Temperature.fits.gz"
     endif
     !  Get an unused Logical Unit Number to use to open the FITS file.
     status=0
     call ftgiou (unit,status)

     !  Create the new empty FITS file.
     blocksize=1
     call ftinit(unit,trim(filename),blocksize,status)

     !  Initialize parameters about the FITS image
     simple=.true.
     ! le signe - signifie que l'on ecrit des reels dans le fits
     bitpix=-32
     extend=.true.

     if (l3D) then
        naxis=3
        naxes(1)=n_rad
        naxes(2)=nz
        naxes(3)=n_az

        !  Write the required header keywords.
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        !call ftphps(unit,simple,bitpix,naxis,naxes,status)

        !  Write the array to the FITS file.
        group=1
        fpixel=1
        nelements=naxes(1)*naxes(2)*naxes(3)

        ! le e signifie real*4
        call ftppre(unit,group,fpixel,nelements,temperature(:,1:nz,:),status)
     else
        naxis=2
        naxes(1)=n_rad
        naxes(2)=nz

        !  Write the required header keywords.
        call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
        !call ftphps(unit,simple,bitpix,naxis,naxes,status)

        !  Write the array to the FITS file.
        group=1
        fpixel=1
        nelements=naxes(1)*naxes(2)

        ! le e signifie real*4
        call ftppre(unit,group,fpixel,nelements,temperature(:,1:nz,1),status)
     endif

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0) then
        call print_error(status)
     end if
  endif

  if (lRE_nLTE .and. iTemperature==1) then
     if (l3D) then
        write(*,*) "ERROR : 3D nLTE version not written yet"
        stop
     endif

     filename = trim(data_dir)//"/Temperature.fits.gz"
     !  Get an unused Logical Unit Number to use to open the FITS file.
     status=0
     call ftgiou (unit,status)

     !  Create the new empty FITS file.
     blocksize=1
     call ftinit(unit,trim(filename),blocksize,status)

     !  Initialize parameters about the FITS image
     simple=.true.
     ! le signe - signifie que l'on ecrit des reels dans le fits
     bitpix=-32
     extend=.true.

     naxis=3
     naxes(1)=n_rad
     naxes(2)=nz
     naxes(3)=n_grains_RE_nLTE

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,temperature_1grain,status)

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0) then
        call print_error(status)
     end if
  endif



  if (lnRE .and. iTemperature==1) then
     if (l3D) then
        write(*,*) "ERROR : 3D nRE version not written yet"
        stop
     endif

     filename = trim(data_dir)//"/Temperature_nRE.fits.gz"
     !  Get an unused Logical Unit Number to use to open the FITS file.
     status=0
     call ftgiou (unit,status)

     !  Create the new empty FITS file.
     blocksize=1
     call ftinit(unit,trim(filename),blocksize,status)

     !  Initialize parameters about the FITS image
     simple=.true.
     extend=.true.

     ! 1er HDU : Temperature ou Proba Temparature
     bitpix=32
     naxis=3
     naxes(1)=n_rad
     naxes(2)=nz
     naxes(3)=n_grains_nRE


     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)

     ! le j signifie integer
     do i=1, n_rad
        do j=1,nz
           do l=grain_nRE_start, grain_nRE_end
              if (l_RE(i,j,l)) then
                 tmp(i,j,l) = 1
              else
                 tmp(i,j,l) = 0
              endif
           enddo
        enddo
     enddo
     call ftpprj(unit,group,fpixel,nelements,tmp,status)

     ! 2eme HDU : Proba Temperature
     call FTCRHD(unit, status)
     bitpix=-32


     naxis=4
     naxes(1)=n_T
     naxes(2)=n_rad
     naxes(3)=nz
     naxes(4)=n_grains_nRE

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,Proba_Temperature,status)


     ! 3eme HDU : Temperature
     call FTCRHD(unit, status)
     bitpix=-32

     naxis=3
     naxes(1)=n_rad
     naxes(2)=nz
     naxes(3)=n_grains_nRE

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,temperature_1grain_nRE,status)

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0) then
        call print_error(status)
     end if
  endif

  return

end subroutine ecriture_temperature


!********************************************************************

subroutine ecriture_Tex(imol)

  implicit none

  integer, intent(in) :: imol

  integer :: i, j, iTrans, iUp, iLow, k
  real(kind=db) :: nUp, nLow, cst

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id

  logical :: simple, extend
  character(len=512) :: filename

  !real, dimension(n_rad,nz,nTrans_tot) :: Tex
  real, dimension(:,:,:), allocatable :: Tex

  allocate(Tex(n_rad,nz,nTrans_tot), stat = alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error Tex in ecriture_Tex'
     stop
  endif
  Tex = 0.0

  k=1 ! Cette version n'est pas en 3D

  do iTrans=1,nTrans_tot
     iUp = iTransUpper(iTrans)
     iLow = iTransLower(iTrans)
     cst = - hp * Transfreq(iTrans) / kb

     do j=1, nz
        do i=1, n_rad
           nUp = tab_nLevel(i,j,k,iUp)
           nLow =  tab_nLevel(i,j,k,iLow)
           if ((nUp > tiny_real) .and. (nLow > tiny_real) ) then
              Tex(i,j,iTrans) = cst / log(  (nUp * poids_stat_g(iLow))  / (nLow * poids_stat_g(iUp) ))
           endif
        enddo ! i
     enddo !j
  enddo !iTrans

  filename = trim(data_dir2(imol))//"/Tex.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=3
  naxes(1)=n_rad
  naxes(2)=nz
  naxes(3)=nTrans_tot

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,Tex,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  deallocate(Tex)

  return

end subroutine ecriture_Tex

!********************************************************************

subroutine taille_moyenne_grains()

  implicit none

  real(kind=db) :: somme
  integer ::  i, j, k, l
  real, dimension(n_rad,nz) :: a_moyen

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id

  logical :: simple, extend
  character(len=512) :: filename


  write(*,*) "Writing average_grain_size.fits.gz"

  a_moyen(:,:) = 0.

  do j=1,nz
     do i=1,n_rad
        somme=0.0
        do l=1, n_grains_tot
           a_moyen(i,j) = a_moyen(i,j) + densite_pouss(i,j,1,l) * r_grain(l)**2
           somme = somme + densite_pouss(i,j,1,l)
        enddo
        a_moyen(i,j) = a_moyen(i,j) / somme
     enddo
  enddo
  a_moyen = sqrt(a_moyen)

  filename = "average_grain_size.fits.gz"

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=2
  naxes(1)=n_rad
  naxes(2)=nz

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,a_moyen,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return

end subroutine taille_moyenne_grains

!******************************************************************************

subroutine ecriture_sed(ised)

  implicit none

  integer, intent(in) :: ised

  integer :: n,k, i, j, lambda
  real :: flux, flux_q, flux_u, flux_v, flux_star, flux_star_scat, flux_disk, flux_disk_scat, facteur, nbre_photons
  real, dimension(n_lambda2) :: n_photons_envoyes

  character(len=512) :: filename

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id

  logical :: simple, extend

  !! Ecriture SED
  if (ised ==1) then
     filename = trim(data_dir)//"/sed1.fits.gz"
  else
     filename = trim(data_dir)//"/sed2.fits.gz"
     if (lProDiMo) write(*,*) "WARNING: sed2.fits.gz will be noisy due to the use of ProDiMo mode"
  endif

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.
  group=1
  fpixel=1


  if (ised == 1) then
     naxis=3
     naxes(1)=n_lambda
     naxes(2)=N_thet
     naxes(3)=N_phi

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     nelements=naxes(1)*naxes(2)*naxes(3)

     do lambda=1,n_lambda
        facteur =  E_photon * tab_lambda(lambda)/tab_delta_lambda(lambda)
        sed1_io(lambda,:,:) = sum(sed(:,lambda,:,:),dim=1) * facteur
     enddo

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,sed1_io,status)

     L_bol1 = sum(sed) * E_photon

  else
     ! Energie totale emise a une distance emise egale au rayon stellaire
     ! on chosit cette distance pour calibrer le flux / pi*B(lambda)
     do lambda=1, n_lambda2
        n_photons_envoyes(lambda) = real(sum(n_phot_envoyes(lambda,:)))
        E_totale(lambda) = E_totale(lambda)/n_photons_envoyes(lambda)
     enddo

     naxis=4
     naxes(1)=n_lambda2
     naxes(2)=N_thet
     naxes(3)=N_phi
     naxes(4)=9

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

     do lambda=1,n_lambda2
        facteur = E_totale(lambda) * tab_lambda(lambda) * 1.0e-6
        sed2_io(lambda,:,:,1) = sum(sed(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,2) = sum(sed_q(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,3) = sum(sed_u(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,4) = sum(sed_v(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,5) = sum(sed_star(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,6) = sum(sed_star_scat(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,7) = sum(sed_disk(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,8) = sum(sed_disk_scat(:,lambda,:,:),dim=1) * facteur
        sed2_io(lambda,:,:,9) = sum(n_phot_sed(:,lambda,:,:),dim=1)
     enddo

     call ftppre(unit,group,fpixel,nelements,sed2_io,status)

     call ftpkys(unit,'FLUX_1',"I = total flux",' ',status)
     call ftpkys(unit,'FLUX_2',"Q",' ',status)
     call ftpkys(unit,'FLUX_3',"U",' ',status)
     call ftpkys(unit,'FLUX_4',"V",' ',status)
     call ftpkys(unit,'FLUX_5',"direct star light",' ',status)
     call ftpkys(unit,'FLUX_6',"scattered star light",' ',status)
     call ftpkys(unit,'FLUX_7',"direct thermal emission",' ',status)
     call ftpkys(unit,'FLUX_8',"scattered thermal emssion",' ',status)
     call ftpkys(unit,'FLUX_9',"number of packets",' ',status)

     L_bol2 = 0.0
     do lambda=1,n_lambda
        L_bol2 = L_bol2 + sum(sed(:,lambda,:,:))*tab_delta_lambda(lambda)*1.0e-6*E_totale(lambda)
     enddo

     if ((lsed_complete).and.(ltemp)) then
        write(*,*) "Flux conservation to within ",  (L_bol2-L_bol1)/L_bol1 * 100, "%  between the two methods"
     endif

  endif ! ised

  call ftpkys(unit,'BUNIT',"W.m-2",' ',status)

  ! Second HDU avec longueur d'onde
  call FTCRHD(unit, status)
  bitpix=-32
  naxis=1
  naxes(1)=n_lambda
  nelements=naxes(1)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,tab_lambda,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

!  write(*,*) "Nbre paquets recus", real(nbre_photons_loop)*real(nbre_photons_eq_th), sum(sed), sum(sed)/(real(nbre_photons_loop)*real(nbre_photons_eq_th))

  return


end subroutine ecriture_sed

!**********************************************************************

subroutine ecriture_pops(imol)

  implicit none

  integer, intent(in) :: imol

  character(len=512) :: filename
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(3) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id, iv, iTrans
  logical :: simple, extend

  real, dimension(n_rad,nz,1,nLevels) :: tab_nLevel_io ! pas en 3D

  filename = trim(data_dir2(imol))//'/populations.fits.gz'

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.

  naxis=3
  naxes(1)=n_rad
  naxes(2)=nz
  naxes(3)=nLevels

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)

  tab_nLevel_io = tab_nLevel
  call ftppre(unit,group,fpixel,nelements,tab_nLevel_io,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) call print_error(status)

  return

end subroutine ecriture_pops

!**********************************************************************

subroutine ecriture_spectre(imol)

  implicit none

  integer, intent(in) :: imol

  character(len=512) :: filename
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements, alloc_status, id, iv, iTrans, xcenter,i, iiTrans
  logical :: simple, extend

  real, dimension(:,:,:,:), allocatable ::  O ! nv, nTrans, n_rad, nz
  real, dimension(nTrans) :: freq

  integer :: n_speed_rt, nTrans_raytracing
  real :: pixel_scale_x, pixel_scale_y

  filename = trim(data_dir2(imol))//'/lines.fits.gz'

  n_speed_rt = mol(imol)%n_speed_rt
  nTrans_raytracing = mol(imol)%nTrans_raytracing

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  simple=.true.
  extend=.true.
  group=1
  fpixel=1

  !------------------------------------------------------------------------------
  ! Line map
  !------------------------------------------------------------------------------
  bitpix=-32
  naxis=5
  if (RT_line_method==1) then
     naxes(1)=1
     naxes(2)=1
  else
     naxes(1)=igridx
     naxes(2)=igridy
  endif
  naxes(3)=2*n_speed_rt+1
  naxes(4)=ntrans
  naxes(5)=RT_n_ibin
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)

  ! create new hdu
  !call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-3,'RAD',status)
  call ftpkyj(unit,'CRPIX1',igridx/2+1,'',status)
  pixel_scale_x = -map_size / (igridx * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-3,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-3,'DEC',status)
  call ftpkyj(unit,'CRPIX2',igridy/2+1,'',status)
  pixel_scale_y = map_size / (igridy * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-3,'pixel scale y [deg]',status)

!  call ftpkye(unit,'vmax_center',mol(imol)%vmax_center_output,-8,'m/s',status)

  if (l_sym_ima) then ! BUG : ca inverse aussi la pente du Cmb mais c'est pas tres grave
     if (RT_line_method==1) then
        ! On ajoute les 2 parties du spectres
        do iv = -n_speed_rt, -1
           spectre(1,1,iv,:,:) = spectre(1,1,iv,:,:) + spectre(1,1,-iv,:,:)
        enddo
        ! On symetrise
        do iv =1, n_speed_rt
           spectre(1,1,iv,:,:) = spectre(1,1,-iv,:,:)
        enddo
        spectre(1,1,0,:,:) = spectre(1,1,0,:,:) * 2.
        ! On divise par deux
        spectre = spectre * 0.5
     else
        xcenter = igridx/2 + modulo(igridx,2)
        if (lkeplerian) then ! profil de raie inverse des 2 cotes
           do i=xcenter+1,igridx
              do iv=-n_speed_rt,n_speed_rt
                 spectre(i,:,iv,:,:) = spectre(igridx-i+1,:,-iv,:,:)
              enddo
           enddo
        else ! infall : meme profil de raie des 2 cotes
           do i=xcenter+1,igridx
              spectre(i,:,:,:,:) = spectre(igridx-i+1,:,:,:,:)
           enddo
        endif
     endif ! lkeplerian
  endif ! l_sym_image

  !  Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,spectre,status)

  !------------------------------------------------------------------------------
  ! HDU 2 : Continuum map
  !------------------------------------------------------------------------------
  bitpix=-32
  naxis=4
  if (RT_line_method==1) then
     naxes(1)=1
     naxes(2)=1
  else
     naxes(1)=igridx
     naxes(2)=igridy
  endif
  naxes(3)=ntrans
  naxes(4)=RT_n_ibin
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

  ! create new hdu
  call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  if (l_sym_ima.and.(RT_line_method==2)) then
     xcenter = igridx/2 + modulo(igridx,2)
     do i=xcenter+1,igridx
        continu(i,:,:,:) = continu(igridx-i+1,:,:,:)
     enddo
  endif ! l_sym_image

  !  Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,continu,status)

  !------------------------------------------------------------------------------
  ! HDU 3 : Transition numbers
  !------------------------------------------------------------------------------
  bitpix=32
  naxis = 1
  naxes(1) = ntrans
  nelements = naxes(1)

  ! create new hdu
  call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  call ftpprj(unit,group,fpixel,nelements,indice_Trans,status)

  !------------------------------------------------------------------------------
  ! HDU 4 : Transition frequencies
  !------------------------------------------------------------------------------
  bitpix=-32
  naxis = 1
  naxes(1) = ntrans
  nelements = naxes(1)

  ! create new hdu
  call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  do i=1,ntrans
     iiTrans= indice_Trans(i)
     freq(i) = Transfreq(iiTrans)
  enddo

  !  Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,freq,status)

  !------------------------------------------------------------------------------
  ! HDU 5 : Velocities
  !------------------------------------------------------------------------------
  bitpix=-32
  naxis = 1
  naxes(1) = 2*n_speed_rt+1
  nelements = naxes(1)

  ! create new hdu
  call ftcrhd(unit, status)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,real(tab_speed_rt),status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)


  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  endif


  !------------------------------------------------------------------------------
  ! Origine
  !------------------------------------------------------------------------------
  if (lorigine) then
     allocate(O(-n_speed_rt:n_speed_rt,nTrans_raytracing,n_rad,nz))
     O = 0.0;
     do i=1, nb_proc
        O(:,:,:,:) =  O(:,:,:,:) +  origine_mol(:,:,:,:,i)
     enddo


     filename = trim(data_dir2(imol))//'/origine.fits.gz'

     !  Get an unused Logical Unit Number to use to open the FITS file.
     status=0
     call ftgiou (unit,status)

     !  Create the new empty FITS file.
     blocksize=1
     call ftinit(unit,trim(filename),blocksize,status)

     simple=.true.
     ! le signe - signifie que l'on ecrit des reels dans le fits
     bitpix=-32
     extend=.true.

     naxis=4
     naxes(1)=2*n_speed_rt+1
     naxes(2)=ntrans
     naxes(3)=n_rad
     naxes(4)=nz

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

     call ftppre(unit,group,fpixel,nelements,O,status)

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0) then
        call print_error(status)
     endif

  endif ! lorigine

  return

end subroutine ecriture_spectre

!**********************************************************************

subroutine cfitsWrite(filename,tab,dim)

  implicit none

  character(len=*), intent(in) :: filename
  real, dimension(*), intent(in) :: tab
  integer, dimension(:), intent(in) :: dim ! dim == shape(tab)

  integer :: status,unit,blocksize,bitpix,naxis
  integer :: group,fpixel,nelements, alloc_status
  logical :: simple, extend

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.


  !  Write the required header keywords.
  naxis=size(dim)
  call ftphpr(unit,simple,bitpix,naxis,dim,0,1,extend,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=product(dim)

  call ftppre(unit,group,fpixel,nelements,tab,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  endif

  return

end subroutine cfitsWrite

!**********************************************************************

end module output
