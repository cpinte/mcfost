module output

  use parametres
  use grains
  use temperature
  use radiation_field, only : J0, xJ_abs
  use thermal_emission, only : E_totale, nbre_reemission
  use constantes
  use dust_prop
  use molecular_emission
  use dust_ray_tracing, only : n_phot_envoyes, Stokes_ray_tracing, tau_surface, stars_map
  use utils
  use Voronoi_grid
  use grid
  use fits_utils, only : cfitsWrite, print_error
  use wavelengths
  use stars, only : E_stars
  use thermal_emission, only : L_packet_th
  use density

  implicit none
  save

  real(kind=dp), dimension(:,:,:,:,:,:), allocatable :: STOKEI, STOKEQ, STOKEU, STOKEV !id,lambda,x,y,n_thet,n_phi
  real(kind=dp), dimension(:,:,:,:,:,:), allocatable :: STOKEI_star, STOKEI_star_scat, STOKEI_disk, STOKEI_disk_scat
  real, dimension(:,:,:,:,:), allocatable :: stoke_io ! x,y, theta, phi, type

  real(kind=dp), dimension(:,:,:,:), allocatable :: sed, n_phot_sed, sed_u, sed_q, sed_v
  real(kind=dp), dimension(:,:,:,:), allocatable, target :: sed_star, sed_star_scat, sed_disk, sed_disk_scat!id,lambda,n_thet,n_phi,x,y
  real, dimension(:,:,:), allocatable :: sed1_io
  real, dimension(:,:,:,:), allocatable :: sed2_io
  real, dimension(:,:), allocatable :: wave2_io

  real(kind=dp), dimension(:,:,:), allocatable :: disk_origin
  real(kind=dp), dimension(:,:), allocatable :: star_origin

  ! Line transfer
  real, dimension(:,:,:,:,:,:), allocatable :: spectre ! speed,trans,thetai,phi,x,y
  real, dimension(:,:,:,:,:), allocatable :: continu ! trans,thetai,phi,x,y


  real(kind=dp) :: sin_disk, cos_disk


  contains

subroutine allocate_mc_images()

  integer :: alloc_status

  allocate(STOKE_io(npix_x,npix_y,capt_debut:capt_fin,N_phi,N_type_flux), stat=alloc_status)
  if (alloc_status > 0) call error( 'Allocation error STOKE_io')
  STOKE_io = 0.0

  allocate(STOKEI(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error STOKEI')
  STOKEI = 0.0

  if (lsepar_pola) then
     allocate(STOKEQ(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEQ')
     STOKEQ = 0.0

     allocate(STOKEU(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEU')
     STOKEU=0.0

     allocate(STOKEV(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEV')
     STOKEV = 0.0
  endif

  if (lsepar_contrib) then
     allocate(STOKEI_star(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEI_star')
     STOKEI_star = 0.0

     allocate(STOKEI_star_scat(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEI_star_scat')
     STOKEI_star_scat = 0.0

     allocate(STOKEI_disk(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEI_disk')
     STOKEI_disk = 0.0

     allocate(STOKEI_disk_scat(n_lambda,npix_x,npix_y,capt_debut:capt_fin,N_phi,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error STOKEI_disk_scat')
     STOKEI_disk_scat = 0.0
  endif !lsepar

  cos_disk = cos(ang_disque * deg_to_rad)
  sin_disk = sin(ang_disque * deg_to_rad)

  return

end subroutine allocate_mc_images

!**********************************************************************

subroutine allocate_sed()

  integer :: alloc_status

  allocate(sed(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed')
  sed = 0.0

  allocate(sed_q(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_q')
  sed_q = 0.0

  allocate(sed_u(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_u')
  sed_u = 0.0

  allocate(sed_v(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_v')
  sed_v = 0.0

  allocate(sed_star(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_star')
  sed_star = 0.0

  allocate(sed_star_scat(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_star_scat')
  sed_star_scat = 0.0

  allocate(sed_disk(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_disk')
  sed_disk = 0.0

  allocate(sed_disk_scat(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_disk_scat')
  sed_disk_scat = 0.0

  allocate(n_phot_sed(n_lambda,N_thet,N_phi,nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error n_phot_sed')
  n_phot_sed = 0.0

  allocate(sed1_io(n_lambda,N_thet,N_phi),sed2_io(n_lambda,N_thet,N_phi,9),wave2_io(n_lambda,2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error sed_io')
  sed1_io=0.0
  sed2_io=0.0
  wave2_io=0.0

  return

end subroutine allocate_sed

!**********************************************************************

subroutine deallocate_sed()

  deallocate(sed,sed_q,sed_u,sed_v)
  deallocate(sed_star,sed_star_scat,sed_disk,sed_disk_scat,n_phot_sed)
  deallocate(sed1_io,sed2_io,wave2_io)

  return

end subroutine deallocate_sed

!**********************************************************************

subroutine allocate_mol_maps(imol)

  integer, intent(in) :: imol

  integer :: alloc_status, n_speed_rt, nTrans_raytracing

  n_speed_rt = mol(imol)%n_speed_rt
  nTrans_raytracing = mol(imol)%nTrans_raytracing

 if (npix_x_save > 1) then
     RT_line_method = 2 ! creation d'une carte avec pixels carres
     npix_x = npix_x_save ; npix_y = npix_y_save ! we update the value after the SED calculation

     write(*,*) "WARNING : memory size if lots of pixels"
     allocate(spectre(npix_x,npix_y,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_incl,RT_n_az), &
          continu(npix_x,npix_y,nTrans_raytracing,RT_n_incl,RT_n_az), stars_map(npix_x,npix_y,1), stat=alloc_status)
  else
     RT_line_method = 1 ! utilisation de pixels circulaires
     allocate(spectre(1,1,-n_speed_rt:n_speed_rt,nTrans_raytracing,RT_n_incl,RT_n_az), &
          continu(1,1,nTrans_raytracing,RT_n_incl,RT_n_az), stars_map(1,1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) call error('Allocation error spectre')
  spectre=0.0
  continu=0.0
  stars_map=0.0

  return

end subroutine allocate_mol_maps

!**********************************************************************

subroutine deallocate_mol_maps()

  deallocate(spectre,continu)
  return

end subroutine deallocate_mol_maps

!**********************************************************************

subroutine allocate_origin()

  integer :: alloc_status

  if (allocated(disk_origin)) deallocate(disk_origin, star_origin)
  allocate(disk_origin(n_lambda2, n_cells, nb_proc), star_origin(n_lambda2, nb_proc), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error disk_origin')
  disk_origin = 0.0
  star_origin = 0.0

  return

end subroutine allocate_origin

!**********************************************************************

subroutine deallocate_origin()

  deallocate(disk_origin, star_origin)
  return

end subroutine deallocate_origin

!**********************************************************************

subroutine capteur(id,lambda,icell,xin,yin,zin,uin,vin,win,stokin,flag_star,flag_scatt, capt)

  implicit none

  real(kind=dp), intent(in) ::  xin,yin,zin, uin,vin,win
  integer, intent(in) :: id, lambda, icell
  real(kind=dp), dimension(4), intent(in)  :: stokin
  logical, intent(in) :: flag_star, flag_scatt
  integer, intent(out) :: capt

  real(kind=dp), dimension(4)  :: stok
  real(kind=dp) ::  x1,y1,z1,u1,v1,w1, xprim, yprim, zprim, ytmp, ztmp, size_pix
  integer :: c_phi, imap1, jmap1, imap2, jmap2, deltapix_x, deltapix_y, maxigrid


  maxigrid = max(npix_x, npix_y)
  if (npix_x == npix_y) then
     deltapix_x = 1
     deltapix_y = 1
  else if (npix_x > npix_y) then
     deltapix_x = 1
     deltapix_y = 1 - (npix_x/2) + (npix_y/2)
  else
     deltapix_x = 1 - (npix_y/2) + (npix_x/2)
     deltapix_y = 1
  endif
  size_pix=maxigrid/map_size


  x1=xin ; y1=yin ; z1=zin
  u1=uin ; v1=vin ; w1=win
  stok=stokin

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

  ! Origine du paquet
  if (lorigine) then
     if (capt == capt_interet) then
        if (flag_star) then
           star_origin(lambda,id) = star_origin(lambda,id) + stok(1)
        else
           disk_origin(lambda, icell,id) = disk_origin(lambda, icell, id) + stok(1)
        endif
     endif
  endif

  ! If we do not keep the MC map, we do not care about the pixels, etc
  if (lmono0 .and. .not.loutput_mc) return

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
        c_phi = int(modulo(atan2(U1,V1)+pi/2,2*real(pi,kind=dp))/(2*pi)*N_phi) + 1
     endif
  endif !l_sym_axiale

  if (c_phi == (N_phi + 1)) then
     c_phi = N_phi
  elseif (c_phi==0) then
     c_phi=1
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

     !*     LA CARTE NPIX_X * NPIX_Y EST PRODUITE DANS LE PLAN Y'Z'
     !*
     !*     utilisation de la symetrie gauche-droite, selon l'axe Y'
     !*     changement de yprim -> -yprim si yprim<0
     !*
     !        if (YPRIM < 0.0) then
     !           YPRIM = -YPRIM
     !           STOK(3,1) = -STOK(3,1)
     !        endif

     !        IMAP = int(((YPRIM / size_neb) * NPIX_Y / 2.) + 0.5) + 1
     !        if (IMAP <= 0 ) IMAP = 1
     !        if (IMAP >= (NPIX_X + 1)) IMAP = NPIX_X

     ! Rotation eventuelle du disque
     ytmp = yprim
     ztmp = zprim
     yprim = ytmp * cos_disk + ztmp * sin_disk
     zprim = ztmp * cos_disk - ytmp * sin_disk

     IMAP1 = int((YPRIM*zoom + 0.5*map_size)*size_pix) + deltapix_x
     if (IMAP1 <= 0 ) return !cycle photon
     if (IMAP1 > NPIX_X)  return !cycle photon

     JMAP1 = int((ZPRIM*zoom + 0.5*map_size)*size_pix)  + deltapix_y
     if (JMAP1 <= 0) return !cycle photon
     if (JMAP1 > NPIX_Y) return !cycle photon

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
        if (IMAP2 > NPIX_X)  return !cycle photon


        JMAP2 = int((ZPRIM*zoom + 0.5*map_size)*size_pix)  + deltapix_y
        if (JMAP2 <= 0) return !cycle photon
        if (JMAP2 > NPIX_Y) return !cycle photon


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

  else ! Creation sed
     sed(lambda,capt,c_phi,id) = sed(lambda,capt,c_phi,id) + stok(1)
     sed_q(lambda,capt,c_phi,id) = sed_q(lambda,capt,c_phi,id) + stok(2)
     sed_u(lambda,capt,c_phi,id) = sed_u(lambda,capt,c_phi,id) + stok(3)
     sed_v(lambda,capt,c_phi,id) = sed_v(lambda,capt,c_phi,id) + stok(4)
     n_phot_sed(lambda,capt,c_phi,id) = n_phot_sed(lambda,capt,c_phi,id) + 1.0_dp
     if (flag_star) then ! photon étoile
        if (flag_scatt) then
           sed_star_scat(lambda,capt,c_phi,id) = sed_star_scat(lambda,capt,c_phi,id) + stok(1)
        else
           sed_star(lambda,capt,c_phi,id) = sed_star(lambda,capt,c_phi,id) + stok(1)
        endif
     else ! photon thermique
        if (flag_scatt) then
           sed_disk_scat(lambda,capt,c_phi,id) = sed_disk_scat(lambda,capt,c_phi,id) + stok(1)
        else
           sed_disk(lambda,capt,c_phi,id) = sed_disk(lambda,capt,c_phi,id) + stok(1)
        endif
     endif ! type de photon
  endif

  return

end subroutine capteur

!**********************************************************************

subroutine write_stokes_fits()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements, id

  integer :: lambda=1

  character(len = 512) :: filename
  logical :: simple, extend

  real, dimension(n_lambda) :: n_photons_envoyes
  real :: facteur, pixel_scale_x, pixel_scale_y

  real :: o_star, frac_star, somme_disk
  real, dimension(n_cells) :: o_disk

  real(kind=dp) :: cos_disk_x2, sin_disk_x2

  cos_disk_x2 = cos(2.* ang_disque * deg_to_rad)
  sin_disk_x2 = sin(2.* ang_disque * deg_to_rad)

  if (lorigine) then

     o_star = sum(star_origin(1,:))
     o_disk(:) = sum(disk_origin(1,:,:), dim=2)
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

     naxis=1
     naxes(1)=n_cells

     ! Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)

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
        if (abs(ang_disque) > 0.) then ! Rotation Q and U
           stoke_io(:,:,:,:,2)=stoke_io(:,:,:,:,2) + stokeq(lambda,:,:,:,:,id) * cos_disk_x2 + &
                stokeu(lambda,:,:,:,:,id) * sin_disk_x2
           stoke_io(:,:,:,:,3)=stoke_io(:,:,:,:,3) - stokeq(lambda,:,:,:,:,id) * sin_disk_x2 + &
                stokeu(lambda,:,:,:,:,id) * cos_disk_x2
        else ! No need for rotation
           stoke_io(:,:,:,:,2)=stoke_io(:,:,:,:,2)+stokeq(lambda,:,:,:,:,id)
           stoke_io(:,:,:,:,3)=stoke_io(:,:,:,:,3)+stokeu(lambda,:,:,:,:,id)
        endif
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
  naxes(1)=npix_x
  naxes(2)=npix_y
  naxes(3)=N_incl
  naxes(4)=N_phi
  naxes(5)=N_type_flux

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)


  ! Write  optional keywords to the header
  !  wavelength
  call ftpkyd(unit,'WAVE',tab_lambda(lambda),-7,'wavelength [microns]',status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
  call ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
  pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
  call ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
  pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)
  call ftpkys(unit,'BUNIT',"W.m-2.pixel-1",'lambda.F_lambda',status)

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
  integer :: i,j,group,fpixel,nelements, alloc_status, xcenter, lambda, itype, ibin, iaz

  character(len = 512) :: filename
  logical :: simple, extend
  real :: pixel_scale_x, pixel_scale_y, W2m2_to_Jy, Q, U
  real(kind=dp) :: sin_disk, cos_disk, cos_disk_x2, sin_disk_x2



  ! Allocation dynamique pour passer en stack
  real, dimension(:,:,:,:,:), allocatable :: image
  real, dimension(:,:,:,:), allocatable :: image_casa

  if (lcasa) then
     allocate(image_casa(npix_x, npix_y, 1, 1), stat=alloc_status) ! 3eme axe : pola, 4eme axe : frequence
     if (alloc_status > 0) call error('Allocation error RT image_casa')
     image_casa = 0.0 ;
  else
     allocate(image(npix_x, npix_y, RT_n_incl, RT_n_az, N_type_flux), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error RT image')
     image = 0.0 ;
  endif

  cos_disk = cos(ang_disque * deg_to_rad)
  sin_disk = sin(ang_disque * deg_to_rad)
  cos_disk_x2 = cos(2.* ang_disque * deg_to_rad)
  sin_disk_x2 = sin(2.* ang_disque * deg_to_rad)

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

  if (lcasa) then
     naxis=4
     naxes(1)=npix_x
     naxes(2)=npix_y
     naxes(3)= 1 ! pola
     naxes(4)=1 ! freq
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
  else
     naxis=5
     naxes(1)=npix_x
     naxes(2)=npix_y
     naxes(3)= RT_n_incl
     naxes(4)= RT_n_az
     naxes(5)=N_type_flux
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  !  wavelength
  call ftpkyd(unit,'WAVE',tab_lambda(lambda),-7,'wavelength [microns]',status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
  call ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
  pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
  call ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
  pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)

  if (lcasa) then
     call ftpkys(unit,'BUNIT',"JY/PIXEL",'',status)
     call ftpkyd(unit,'RESTFREQ',c_light/(tab_lambda(lambda)*1e-6),-14,'Hz',status)
     call ftpkys(unit,'BTYPE',"Intensity",' ',status)

     ! 3eme axe
     call ftpkys(unit,'CTYPE3',"STOKES",' ',status)
     call ftpkye(unit,'CRVAL3',1.0,-7,'',status)
     call ftpkye(unit,'CDELT3',1.0,-7,'',status)
     call ftpkyj(unit,'CRPIX3',1,'',status)

     ! 4eme axe
     call ftpkys(unit,'CTYPE4',"FREQ",' ',status)
     call ftpkyd(unit,'CRVAL4',c_light/(tab_lambda(lambda)*1e-6),-14,'Hz',status)
     call ftpkye(unit,'CDELT4',2e9,-7,'Hz',status) ! 2GHz by default
     call ftpkyj(unit,'CRPIX4',0,'',status)
  else
     call ftpkys(unit,'BUNIT',"W.m-2.pixel-1",'lambda.F_lambda',status)
  endif

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
  if (lcasa) then
     itype=1 ; ibin=1 ; iaz = 1

     W2m2_to_Jy = 1e26 * (tab_lambda(lambda)*1e-6)/c_light;

     do j=1,npix_y
        do i=1,npix_x
           image_casa(i,j,1,1) = sum(Stokes_ray_tracing(lambda,i,j,ibin,iaz,itype,:)) * W2m2_to_Jy
        enddo !i
     enddo !j

     if (l_sym_ima) then
        xcenter = npix_x/2 + modulo(npix_x,2)
        do i=xcenter+1,npix_x
           image_casa(i,:,1,1) = image_casa(npix_x-i+1,:,1,1)
        enddo
     endif

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,image_casa,status)
  else ! mode non casa
     type_loop : do itype=1,N_type_flux
        do ibin=1,RT_n_incl
           do iaz=1,RT_n_az
              do j=1,npix_y
                 do i=1,npix_x
                    if (lsepar_pola) then
                       if (itype==2) then
                          Q = sum(Stokes_ray_tracing(lambda,i,j,ibin,iaz,2,:))
                          U = sum(Stokes_ray_tracing(lambda,i,j,ibin,iaz,3,:))
                          image(i,j,ibin,iaz,2) = Q * cos_disk_x2 + U * sin_disk_x2
                          image(i,j,ibin,iaz,3) = - Q * sin_disk_x2 + U * cos_disk_x2
                       else if (itype==3) then ! itype 3 already done together with itype 2
                          cycle type_loop
                       else ! Other images do not need rotation
                          image(i,j,ibin,iaz,itype) = sum(Stokes_ray_tracing(lambda,i,j,ibin,iaz,itype,:))
                       endif
                    else ! No need for rotation when there is no pola
                       image(i,j,ibin,iaz,itype) = sum(Stokes_ray_tracing(lambda,i,j,ibin,iaz,itype,:))
                    endif
                 enddo !i
              enddo !j
           enddo ! iaz
        enddo !ibin
     enddo type_loop ! itype

     if (l_sym_ima) then
        xcenter = npix_x/2 + modulo(npix_x,2)
        do i=xcenter+1,npix_x
           image(i,:,:,:,1) = image(npix_x-i+1,:,:,:,1)

           if (lsepar_pola) then
              image(i,:,:,:,2) = image(npix_x-i+1,:,:,:,2)
              image(i,:,:,:,3) = - image(npix_x-i+1,:,:,:,3)
              image(i,:,:,:,4) = image(npix_x-i+1,:,:,:,4)
           endif

           if (lsepar_contrib) then
              image(i,:,:,:,n_Stokes+1) = image(npix_x-i+1,:,:,:,n_Stokes+1)
              image(i,:,:,:,n_Stokes+2) = image(npix_x-i+1,:,:,:,n_Stokes+2)
              image(i,:,:,:,n_Stokes+3) = image(npix_x-i+1,:,:,:,n_Stokes+3)
              image(i,:,:,:,n_Stokes+4) = image(npix_x-i+1,:,:,:,n_Stokes+4)
           endif
        enddo
     endif ! l_sym_image

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,image,status)
  endif

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)
  if (lcasa) then
     deallocate(image_casa)
  else
     deallocate(image)
  endif

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return

end subroutine ecriture_map_ray_tracing

!**********************************************************************

subroutine write_tau_surface()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, alloc_status, xcenter, lambda, itype, ibin, iaz

  character(len = 512) :: filename
  logical :: simple, extend
  real :: pixel_scale_x, pixel_scale_y, W2m2_to_Jy, Q, U

  ! Allocation dynamique pour passer en stack
  real, dimension(:,:,:,:,:), allocatable :: image

  allocate(image(npix_x, npix_y, RT_n_incl, RT_n_az, 3), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tau_surface image')
  image = 0.0 ;

  lambda=1
  filename = trim(data_dir)//"/tau=1_surface.fits.gz"

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
  naxes(1)=npix_x
  naxes(2)=npix_y
  naxes(3)= RT_n_incl
  naxes(4)= RT_n_az
  naxes(5)= 3 !xyz
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  !  wavelength
  call ftpkyd(unit,'WAVE',tab_lambda(lambda),-7,'wavelength [microns]',status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
  call ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
  pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
  call ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
  pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)

  call ftpkys(unit,'BUNIT',"AU",'',status)

  call ftpkys(unit,'COORD_1',"x",' ',status)
  call ftpkys(unit,'COORD_2',"y",' ',status)
  call ftpkys(unit,'COORD_3',"z",' ',status)

  !----- Images
  ! Boucles car ca ne passe pas avec sum directement (ifort sur mac)
  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
        do j=1,npix_y
           do i=1,npix_x
              image(i,j,ibin,iaz,:) = sum(tau_surface(i,j,ibin,iaz,:,:), dim=2)
           enddo !i
        enddo !j
     enddo ! iaz
  enddo !ibin

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

end subroutine write_tau_surface

!**********************************************************************

subroutine ecriture_sed_ray_tracing()

  implicit none

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,group,fpixel,nelements, lambda

  character(len = 512) :: filename
  logical :: simple, extend

  real, dimension(n_lambda2, RT_n_incl, RT_n_az, N_type_flux) :: sed_rt

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
  naxes(2)= RT_n_incl
  naxes(3)= RT_n_az
  naxes(4)=N_type_flux

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  sed_rt = 0.0_dp
  if (RT_sed_method == 1) then
     sed_rt(:,:,:,:) = sum(Stokes_ray_tracing(:,1,1,:,:,:,:),dim=5)
  else
     do i=1,npix_x
        do j=1,npix_y
           sed_rt(:,:,:,:) = sed_rt(:,:,:,:) + sum(Stokes_ray_tracing(:,i,j,:,:,:,:),dim=5)
        enddo
     enddo
  endif

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,sed_rt,status)

  call ftpkys(unit,'BUNIT',"W.m-2",'lambda.F_lambda',status)

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
  call ftppre(unit,group,fpixel,nelements,real(tab_lambda,kind=sp),status)

  call ftpkys(unit,'BUNIT',"micron",' ',status)

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
  integer :: group,fpixel,nelements

  character(len = 512) :: filename
  logical :: simple, extend

  real :: o_star
  real, dimension(n_cells) :: o_disk

  filename = trim(data_dir)//"/origine.fits.gz"

  ! Normalisation
  o_star = sum(star_origin(1,:))
  o_disk(:) = o_disk(:) / (sum(o_disk) + o_star)

  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou ( unit, status )

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,trim(filename),blocksize,status)

  !  Initialize parameters about the FITS image
  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  naxis=1
  naxes(1)=n_cells

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,o_disk,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  return

end subroutine write_origin

!**********************************************************************

subroutine calc_optical_depth_map(lambda)

  implicit none

  integer :: lambda

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: i,j,k, group,fpixel,nelements

  character(len = 512) :: filename
  logical :: simple, extend, lmilieu

  real, dimension(n_rad,nz,n_az,2) :: optical_depth_map

  lmilieu = .true. ! opacite au milieu ou a la "fin" de la cellule

  if (lmilieu) then
     ! Opacite radiale
     do k=1, n_az
        do j=1, nz
           i=1 ; optical_depth_map(i,j,k,1) = kappa(cell_map(i,j,k),lambda)* 0.5 * (r_lim(i)-r_lim(i-1))
           do i=2, n_rad
              optical_depth_map(i,j,k,1) = optical_depth_map(i-1,j,k,1) + &
                   0.5 * kappa(cell_map(i-1,j,k),lambda)*(r_lim(i-1)-r_lim(i-2)) + &
                   0.5 * kappa(cell_map(i,j,k),lambda)*(r_lim(i)-r_lim(i-1))
           enddo
        enddo
     enddo

     ! Opacite verticale
     do i=1, n_rad
        do k=1, n_az
           j=nz ; optical_depth_map(i,j,k,2) = kappa(cell_map(i,j,k),lambda)* 0.5 * (z_lim(i,j+1)-z_lim(i,j))
           do j=nz-1,1,-1
              optical_depth_map(i,j,k,2) = optical_depth_map(i,j+1,k,2) + &
                   0.5 * kappa(cell_map(i,j+1,k),lambda)*(z_lim(i,j+2)-z_lim(i,j+1)) + &
                   0.5 * kappa(cell_map(i,j,k),lambda)*(z_lim(i,j+1)-z_lim(i,j))
           enddo
        enddo
     enddo

  else
     ! Opacite radiale
     do k=1, n_az
        do j=1, nz
           i=1 ; optical_depth_map(i,j,k,1) = kappa(cell_map(i,j,k),lambda)*(r_lim(i)-r_lim(i-1))
           do i=2, n_rad
              optical_depth_map(i,j,k,1) = optical_depth_map(i-1,j,k,1) +kappa(cell_map(i,j,k),lambda)*(r_lim(i)-r_lim(i-1))
           enddo
        enddo
     enddo

     ! Opacite verticale
     do i=1, n_rad
        do k=1, n_az
           j=nz ; optical_depth_map(i,j,k,2) = kappa(cell_map(i,j,k),lambda)*(z_lim(i,j+1)-z_lim(i,j))
           do j=nz-1,1,-1
              optical_depth_map(i,j,k,2) = optical_depth_map(i,j+1,k,2) + kappa(cell_map(i,j,k),lambda)*(z_lim(i,j+1)-z_lim(i,j))
           enddo
        enddo
     enddo
  endif

  write(*,*) "Writing optical_depth_map.fits.gz for wl=", real(tab_lambda(lambda),kind=sp), "microns"
  filename = "!optical_depth_map.fits.gz"

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
  call exit(0)

end subroutine calc_optical_depth_map

!***********************************************************

subroutine write_column_density()
  ! Only works if the star in in 0, 0, 0 at the moment

  integer, parameter :: n_directions = 4

  real, dimension(n_cells,n_directions) :: CD
  integer :: icell, icell0, next_cell, previous_cell
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, direction

  logical :: simple, extend, ltest
  character(len=512) :: filename

  real(kind=dp) :: x0,y0,z0, x1,y1,z1, norme, l, u,v,w, l_contrib, l_void_before

  CD(:,:) = 0.0
  do direction = 1, n_directions
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

        if (direction == 1) then ! to star : test needs to be updated for Voronoi
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

        ! Voronoi : next_cell > 1 (not a wall) and Voronoi(icell)%is_star == .false, (not a star)
        ! Cylindrical : next_cell <= n_cells
        ltest = .true.
        do while(ltest)
           previous_cell = icell0
           icell0 = next_cell
           x0 = x1 ; y0 = y1 ; z0 = z1
           call cross_cell(x0,y0,z0, u,v,w,  icell0, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)
           CD(icell,direction) = CD(icell,direction) + (l_contrib * AU_to_m) * densite_gaz(icell) * masse_mol_gaz

           if (lVoronoi) then
              ltest = (next_cell > 0)
              if (ltest) ltest = (.not.Voronoi(next_cell)%is_star)
           else
              ltest = (next_cell <= n_cells)
           endif
        enddo
     enddo ! icell
  end do ! direction
  CD(:,:) = CD(:,:) / (m_to_cm)**2 ! g/cm^-2

  filename = trim(root_dir)//"/data_disk/column_density.fits.gz"

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

  if (lVoronoi) then
     naxis=2
     naxes(1)=n_cells
     naxes(2)=n_directions
     nelements=naxes(1)*naxes(2)
  else
     if (l3D) then
        naxis=4
        naxes(1)=n_rad
        naxes(2)=2*nz
        naxes(3)=n_az
        naxes(4)=n_directions
        nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
     else
        naxis=3
        naxes(1)=n_rad
        naxes(2)=nz
        naxes(3)=n_directions
        nelements=naxes(1)*naxes(2)*naxes(3)
     endif
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  call ftpkys(unit,'BUNIT',"g.cm-2",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,CD,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  return

end subroutine write_column_density

!***********************************************************


subroutine reemission_stats()

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements

  character(len = 512) :: filename
  logical :: simple, extend

  real, dimension(n_cells) :: N_reemission

  N_reemission = sum(nbre_reemission(:,:),dim=2)

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
  naxis=1
  naxes(1)=n_cells

  extend=.true.

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,N_reemission,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  return

end subroutine reemission_stats

!********************************************************************************

subroutine write_disk_struct(lparticle_density)
! Ecrit les table de densite du gaz en g/cm^3
! de la poussiere en g/cm^3 et en particules
! + coordonnees r et z en AU
! + les tailles et masse des grains
! C. Pinte
! 3/06/06

  implicit none

  logical, intent(in) :: lparticle_density

  integer :: i, j, icell

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, alloc_status

  logical :: simple, extend
  character(len=512) :: filename

  real, dimension(:), allocatable :: dens, vol
  real(kind=dp), dimension(:,:), allocatable :: dust_dens
  real(kind=dp), dimension(:,:,:,:), allocatable :: grid

  write(*,*) "Writing disk structure files in data_disk ..."
  allocate(dens(n_cells), vol(n_cells), dust_dens(n_cells,n_grains_tot), stat = alloc_status)
  if (alloc_status > 0) call error('Allocation error density tables for fits file')

  filename = trim(root_dir)//"/data_disk/gas_density.fits.gz"

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

  if (lVoronoi) then
     naxis=1
     naxes(1) = n_cells
     nelements=naxes(1)
  else
     if (l3D) then
        naxis=3
        naxes(1)=n_rad
        naxes(2)=2*nz
        naxes(3)=n_az
        nelements=naxes(1)*naxes(2)*naxes(3)
     else
        naxis=2
        naxes(1)=n_rad
        naxes(2)=nz
        nelements=naxes(1)*naxes(2)
     endif
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"g.cm^-3",' ',status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  dens = 0.0

  dens(:) = densite_gaz(1:n_cells) * masse_mol_gaz / m3_to_cm3 ! nH2/m**3 --> g/cm**3

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
  if (lparticle_density) then
     filename = trim(root_dir)//"/data_disk/dust_particle_density.fits.gz"

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
     group=1
     fpixel=1

     if (lVoronoi) then
        naxis=2
        naxes(1) = n_cells
        naxes(2) = n_grains_tot
        nelements=naxes(1)*naxes(2)
     else
        if (l3D) then
           naxis=4
           naxes(1)=n_rad
           naxes(2)=2*nz
           naxes(3)=n_az
           naxes(4) = n_grains_tot
           nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           naxis=3
           naxes(1)=n_rad
           naxes(2)=nz
           naxes(3) = n_grains_tot
           nelements=naxes(1)*naxes(2)*naxes(3)
        endif
     endif

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     ! Write  optional keywords to the header
     call ftpkys(unit,'UNIT',"part.m^-3 [per grain size bin N(a).da]",' ',status)

     !  Write the array to the FITS file.
     !  dens =  densite_pouss
     ! le d signifie real*8
     dust_dens(:,:) = 0.0
     do icell=1,n_cells
        dust_dens(icell,:) = densite_pouss(:,icell) * m3_to_cm3  ! Todo : inverting dimensions is not a good idea
     enddo !icell
     call ftpprd(unit,group,fpixel,nelements,dust_dens,status)

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0)  call print_error(status)
  endif ! lparticle_density

  ! ********************************************************************************
  filename = trim(root_dir)//"/data_disk/dust_mass_density.fits.gz"

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

  if (lVoronoi) then
     naxis=1
     naxes(1) = n_cells
     nelements=naxes(1)
  else
     if (l3D) then
        naxis=3
        naxes(1)=n_rad
        naxes(2)=2*nz
        naxes(3)=n_az
        nelements=naxes(1)*naxes(2)*naxes(3)
     else
        naxis=2
        naxes(1)=n_rad
        naxes(2)=nz
        nelements=naxes(1)*naxes(2)
     endif
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"g.cm^-3",' ',status)

  !  Write the array to the FITS file.
  !  dens =  densite_pouss
  ! le d signifie real*8
  dens(:) = 0.0
  do icell=1,n_cells
     dens(icell) = sum(densite_pouss(:,icell) * M_grain(:)) ! M_grain en g
  enddo
  call ftppre(unit,group,fpixel,nelements,dens,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = trim(root_dir)//"/data_disk/grain_sizes.fits.gz"

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
  filename = trim(root_dir)//"/data_disk/grain_sizes_min.fits.gz"

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
  filename = trim(root_dir)//"/data_disk/grain_sizes_max.fits.gz"

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
  filename = trim(root_dir)//"/data_disk/grain_masses.fits.gz"

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
  filename = trim(root_dir)//"/data_disk/volume.fits.gz"

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
  group=1
  fpixel=1

  if (lVoronoi) then
     naxis=1
     naxes(1) = n_cells
     nelements=naxes(1)
  else
     if (l3D) then
        naxis=3
        naxes(1)=n_rad
        naxes(2)=2*nz
        naxes(3)=n_az
        nelements=naxes(1)*naxes(2)*naxes(3)
     else
        naxis=2
        naxes(1)=n_rad
        naxes(2)=nz
        nelements=naxes(1)*naxes(2)
     endif
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  ! Write  optional keywords to the header
  call ftpkys(unit,'UNIT',"AU^3",' ',status)

  !  Write the array to the FITS file.
  vol(:) = volume(:) ! conversion to single precision
  call ftppre(unit,group,fpixel,nelements,vol,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! ********************************************************************************
  filename = trim(root_dir)//"/data_disk/grid.fits.gz"

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
  group=1
  fpixel=1


  if (lVoronoi) then
     naxis=2
     naxes(1)=n_cells
     naxes(2)=3 !xyz
     nelements=naxes(1)*naxes(2)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     ! Write  optional keywords to the header
     call ftpkys(unit,'UNIT',"AU",' ',status)
     call ftpkys(unit,'DIM_1',"x",' ',status)
     call ftpkys(unit,'DIM_2',"y",' ',status)
     call ftpkys(unit,'DIM_3',"z",' ',status)

     allocate(grid(n_cells,3,1,1))
     do icell=1, n_cells
        grid(icell,:,1,1) = Voronoi(icell)%xyz(:)
     enddo

     ! le d signifie real*8
     call ftpprd(unit,group,fpixel,nelements,grid,status)
  else
     naxis=4
     naxes(1)=n_rad
     naxes(2)=nz
     naxes(3)=1
     naxes(4)=2

     if (l3D) then
        naxes(2)=2*nz+1
        naxes(3)=n_az
        naxes(4)=3

        allocate(grid(n_rad,2*nz+1,n_az,3)) ; grid = 0.0
        do i=1, n_rad
           grid(i,:,:,1) = sqrt(r_lim(i) * r_lim(i-1))
           do j=1,nz
              grid(i,nz+1+j,:,2) = (real(j)-0.5)*delta_z(i)
              grid(i,nz+1-j,:,2) = -(real(j)-0.5)*delta_z(i)
           enddo
        enddo

        do i=1, n_az
           grid(:,:,i,3) = (i-0.5)/n_az * deux_pi
        enddo
     else
        allocate(grid(n_rad,nz,1,2))

        do i=1, n_rad
           grid(i,:,1,1) = sqrt(r_lim(i) * r_lim(i-1))
           do j=1,nz
              grid(i,j,1,2) = (real(j)-0.5)*delta_z(i)
           enddo
        enddo
     endif
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     ! Write  optional keywords to the header
     call ftpkys(unit,'UNIT',"AU",' ',status)
     call ftpkys(unit,'DIM_1',"cylindrical radius",' ',status)
     call ftpkys(unit,'DIM_2',"elevation above midplane",' ',status)
     if (l3D) call ftpkys(unit,'DIM_3',"azimuth [rad]",' ',status)

     ! le d signifie real*8
     call ftpprd(unit,group,fpixel,nelements,grid,status)
  endif ! lVoronoi

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
  end if

  ! Wrting the column density
  !call write_column_density()

  if (lstop_after_init) then
     write(*,*) "Exiting"
     call exit(0)
  endif

  write(*,*) "Done"

  return

end subroutine write_disk_struct

!********************************************************************

subroutine ecriture_J(step)
! Ecrit la table du champ de radiation
! C. Pinte
! 26/09/07

  implicit none

  integer, intent(in) :: step

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, lambda, icell

  logical :: simple, extend
  character(len=512) :: filename

  real, dimension(n_cells,n_lambda) :: Jio
  real(kind=dp) :: n_photons_envoyes, energie_photon


  ! Step1
  ! xJ_abs est par bin de lambda donc Delta_lambda.F_lambda
  ! Jio en W.m-2 (lambda.F_lambda)
  ! 1/4pi est inclus dans n_phot_l_tot
  ! teste OK par rapport a fct bb de yorick
  if (step==1) then
     filename = trim(data_dir)//"/J_step1.fits.gz"
     do lambda=1, n_lambda
        do icell=1, n_cells
           Jio(icell,lambda) = sum(xJ_abs(icell,lambda,:) + J0(icell,lambda)) * L_packet_th / volume(icell) &
                * tab_lambda(lambda) / tab_delta_lambda(lambda)
        enddo
     enddo
  else
     ! Step 2
     filename = trim(data_dir)//"/J.fits.gz"
     do lambda=1, n_lambda2
        n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
        energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda)) / n_photons_envoyes * tab_lambda(lambda) * 1.0e-6  !lambda.F_lambda
        do icell=1, n_cells
           Jio(icell,lambda) = sum(xJ_abs(icell,lambda,:) + J0(icell,lambda)) * energie_photon/volume(icell)
        enddo
     enddo
  endif

  write(*,*) "Writing "//trim(filename)

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

  if (lVoronoi) then
     naxis=2
     naxes(1)=n_cells
     naxes(2)=n_lambda
     nelements=naxes(1)*naxes(2)
  else
     if (l3D) then
        naxis=4
        naxes(1)=n_rad
        naxes(2)=2*nz
        naxes(3)=n_az
        naxes(4)=n_lambda
        nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
     else
        naxis=3
        naxes(1)=n_rad
        naxes(2)=nz
        naxes(3)=n_lambda
        nelements=naxes(1)*naxes(2)*naxes(3)
     endif
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,Jio,status)

  call ftpkys(unit,'BUNIT',"W.m-2",'lambda.F_lambda',status)

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
  call ftppre(unit,group,fpixel,nelements,real(tab_lambda,kind=sp),status)

  call ftpkys(unit,'BUNIT',"micron",' ',status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) call print_error(status)

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
  integer :: group,fpixel,nelements, lambda, l, icell

  logical :: simple, extend
  character(len=512) :: filename

  integer, parameter :: n=200

  real(kind=dp) :: n_photons_envoyes, energie_photon
  real(kind=dp), dimension(n_lambda,n_cells) :: J
  real, dimension(n_cells) :: G
  real(kind=dp), dimension(n) :: wl, J_interp
  real(kind=dp), dimension(n_lambda) :: lamb
  real(kind=dp) :: delta_wl


  ! D'après van Dishoeck et al. (2008) le champs FUV de Draine est de 2.67e-3 erg cm-2 s-1
  ! Soit 2.67e-6 W / m2
  ! C'est entre 912 et 2000 Angströms (la borne supérieure varie un peu d'un papier à l'autre).

  filename = trim(data_dir)//"/UV_field.fits.gz"

  write(*,*) "Writing "//trim(filename)

  ! Step 2
  do lambda=1, n_lambda
     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda)) / n_photons_envoyes ! F_lambda here
     do icell=1, n_cells
        J(lambda,icell) = sum(xJ_abs(icell,lambda,:) + J0(icell,lambda)) * energie_photon/volume(icell)
     enddo
  enddo
  lamb(:) = tab_lambda(:) * 1e-6 ! en m

  ! van Dischoeck 2008
  wl(:) = span(0.0912,0.2,n) * 1e-6 ! en m
  delta_wl = (wl(n) - wl(1))/(n-1.) ! en m
  do icell=1,n_cells
     do l=1,n
        J_interp(l) = interp( J(:,icell),lamb(:),wl(l))
     enddo

     ! integration trapeze
     G(icell) = (sum(J_interp(:)) - 0.5 * (J_interp(1) + J_interp(n)) ) * delta_wl / 2.67e-6
  enddo ! icell

!---  ! Le petit et al
!---  wl(:) = span(0.0912,0.24,n) * 1e-6 ! en m
!---  delta_wl = (wl(n) - wl(1))/(n-1.) * 1e-6  ! must be in AA ??
!---
!---  do ri=1,n_rad
!---     do zj=j_start,nz
!---        do phik=1, n_az
!---           do l=1,n
!---              J_interp(l) = interp( J(:,ri,zj,phik),lamb(:),wl(l))
!---           enddo
!---
!---           ! Le Petit et al 2006 page 19-20
!---           ! integration trapeze
!---           G(ri,zj,phik) = (sum(J_interp(:)) - 0.5 * (J_interp(1) + J_interp(n)) ) * delta_wl  &
!---                * 4*pi/c_light / (5.6e-14 * erg_to_J * m_to_cm**3)
!---
!---           Gio(ri,zj,phik) = G(ri,zj,phik) ! Teste OK par comparaison avec yorick
!---        enddo
!---     enddo
!---  enddo

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
     naxes(2)=2*nz+1
     naxes(3)=n_az
     nelements=naxes(1)*naxes(2)*naxes(3)
  else
     naxis=2
     naxes(1)=n_rad
     naxes(2)=nz
     nelements=naxes(1)*naxes(2)
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1

  ! le e signifie real*4
  call ftppre(unit,group,fpixel,nelements,G,status)

  call ftpkys(unit,'BUNIT',"Habing",'[912-2000] AA, following formulae by van Dischoeck et al 2008',status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) call print_error(status)

  return

end subroutine ecriture_UV_field


!********************************************************************

subroutine ecriture_temperature(iTemperature)

  implicit none


  integer, intent(in) :: iTemperature

  integer :: l, icell
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements, alloc_status

  logical :: simple, extend
  character(len=512) :: filename

  integer, dimension(:,:), allocatable :: tmp


  if (lRE_LTE) then
     if (iTemperature == 2) then
        filename = "!"//trim(data_dir)//"/Temperature_DA.fits.gz" ! "!" to overwrite file if computing diffusion approx twice
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

     if (lVoronoi) then
        naxis=1
        naxes(1)=n_cells
        nelements=naxes(1)
     else
        if (l3D) then
           naxis=3
           naxes(1)=n_rad
           naxes(2)=2*nz
           naxes(3)=n_az
           nelements=naxes(1)*naxes(2)*naxes(3)
        else
           naxis=2
           naxes(1)=n_rad
           naxes(2)=nz
           nelements=naxes(1)*naxes(2)
        endif
     endif

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,Tdust(1:n_cells),status)

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0) then
        call print_error(status)
     end if
  endif

  if (lRE_nLTE) then
     if (iTemperature == 2) then
        filename = "!"//trim(data_dir)//"/Temperature_nLTE_DA.fits.gz" ! "!" to overwrite file if computing diffusion approx twice
     else
        filename = trim(data_dir)//"/Temperature_nLTE.fits.gz"
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

     if (lVoronoi) then
        naxis=2
        naxes(1)=n_grains_RE_nLTE
        naxes(2)=n_cells
        nelements=naxes(1)*naxes(2)
     else
        if (l3D) then
           naxis=4
           naxes(1)=n_grains_RE_nLTE
           naxes(2)=n_rad
           naxes(3)=2*nz
           naxes(4)=n_az
           nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           naxis=3
           naxes(1)=n_grains_RE_nLTE
           naxes(2)=n_rad
           naxes(3)=nz
           nelements=naxes(1)*naxes(2)*naxes(3)
        endif
     endif

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
     !call ftphps(unit,simple,bitpix,naxis,naxes,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,Tdust_1grain,status)

     !  Close the file and free the unit number.
     call ftclos(unit, status)
     call ftfiou(unit, status)

     !  Check for any error, and if so print out error messages
     if (status > 0) then
        call print_error(status)
     end if
  endif

  if (lnRE) then
     if (iTemperature == 2) then
        filename = "!"//trim(data_dir)//"/Temperature_nRE_DA.fits.gz" ! "!" to overwrite file if computing diffusion approx twice
     else
        filename = trim(data_dir)//"/Temperature_nRE.fits.gz"
     endif

     !  Get an unused Logical Unit Number to use to open the FITS file.
     status=0
     call ftgiou (unit,status)

     !  Create the new empty FITS file.
     blocksize=1
     call ftinit(unit,trim(filename),blocksize,status)

     !  Initialize parameters about the FITS image
     simple=.true.
     extend=.true.
     group=1
     fpixel=1

     !------------------------------------------------------------------------------
     ! 1er HDU : Temperature d'equilibre
     !------------------------------------------------------------------------------
     bitpix=-32

     if (lVoronoi) then
        naxis=2
        naxes(1)=n_grains_nRE
        naxes(2)=n_cells
        nelements=naxes(1)*naxes(2)
     else
        if (l3D) then
           naxis=4
           naxes(1)=n_grains_nRE
           naxes(2)=n_rad
           naxes(3)=2*nz
           naxes(4)=n_az
           nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           naxis=3
           naxes(1)=n_grains_nRE
           naxes(2)=n_rad
           naxes(3)=nz
           nelements=naxes(1)*naxes(2)*naxes(3)
        endif
     endif

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,Tdust_1grain_nRE,status)

     !------------------------------------------------------------------------------
     ! 2eme HDU : is the grain at equilibrium ?
     !------------------------------------------------------------------------------
     bitpix=32
     if (lVoronoi) then
        naxis=2
        naxes(1)=n_grains_nRE
        naxes(2)=n_cells
        nelements=naxes(1)*naxes(2)
     else
        if (l3D) then
           naxis=4
           naxes(1)=n_grains_nRE
           naxes(2)=n_rad
           naxes(3)=2*nz+1
           naxes(4)=n_az
           nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           naxis=3
           naxes(1)=n_grains_nRE
           naxes(2)=n_rad
           naxes(3)=nz
           nelements=naxes(1)*naxes(2)*naxes(3)
        endif
     endif
     allocate(tmp(grain_nRE_start:grain_nRE_end,n_cells), stat=alloc_status)

     if (alloc_status /= 0) call error("Allocation error in ecriture_temperature")

      ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     ! le j signifie integer
     do icell=1, n_cells
        do l=grain_nRE_start, grain_nRE_end
           if (l_RE(l,icell)) then
              tmp(l,icell) = 1
           else
              tmp(l,icell) = 0
           endif
        enddo
     enddo
     call ftpprj(unit,group,fpixel,nelements,tmp,status)

     !------------------------------------------------------------------------------
     ! 3eme HDU temperature table
     !------------------------------------------------------------------------------
     bitpix=-32
     naxis=1
     naxes(1)=n_T
     nelements=naxes(1)

     ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,tab_Temp,status)

     !------------------------------------------------------------------------------
     ! 4eme HDU : Proba Temperature
     !------------------------------------------------------------------------------
     bitpix=-32
     if (lVoronoi) then
        naxis=3
        naxes(1)=n_T
        naxes(2)=n_grains_nRE
        naxes(3)=n_cells
        nelements=naxes(1)*naxes(2)*naxes(3)
     else
        if (l3D) then
           naxis=5
           naxes(1)=n_T
           naxes(2)=n_grains_nRE
           naxes(3)=n_rad
           naxes(4)=2*nz
           naxes(5)=n_az
           nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
        else
           naxis=4
           naxes(1)=n_T
           naxes(2)=n_grains_nRE
           naxes(3)=n_rad
           naxes(4)=nz
           nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        endif
     endif

     ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,Proba_Tdust,status)

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

  integer :: iTrans, iUp, iLow, k, icell
  real(kind=dp) :: nUp, nLow, cst

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, alloc_status

  logical :: simple, extend
  character(len=512) :: filename

  real, dimension(:,:), allocatable :: Tex

  allocate(Tex(n_cells,nTrans_tot), stat = alloc_status)
  if (alloc_status > 0) call error('Allocation error Tex in ecriture_Tex')
  Tex = 0.0

  k=1 ! Cette version n'est pas en 3D

  do iTrans=1,nTrans_tot
     iUp = iTransUpper(iTrans)
     iLow = iTransLower(iTrans)
     cst = - hp * Transfreq(iTrans) / kb

     do icell=1, n_cells
        nUp = tab_nLevel(icell,iUp)
        nLow =  tab_nLevel(icell,iLow)
        if ((nUp > tiny_real) .and. (nLow > tiny_real) ) then
           Tex(icell,iTrans) = cst / log(  (nUp * poids_stat_g(iLow))  / (nLow * poids_stat_g(iUp) ))
        endif
     enddo ! icell
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

  naxis=2
  naxes(1)=n_cells
  naxes(2)=nTrans_tot

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)

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

  real(kind=dp) :: somme
  integer :: l, icell
  real, dimension(n_cells) :: a_moyen

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements

  logical :: simple, extend
  character(len=512) :: filename


  write(*,*) "Writing average_grain_size.fits.gz"

  a_moyen(:) = 0.

  do icell=1, n_cells
     somme=0.0
     do l=1, n_grains_tot
        a_moyen(icell) = a_moyen(icell) + densite_pouss(l,icell) * r_grain(l)**2
        somme = somme + densite_pouss(l,icell)
     enddo
     a_moyen(icell) = a_moyen(icell) / somme
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

  naxis=1
  naxes(1)=n_cells

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  !call ftphps(unit,simple,bitpix,naxis,naxes,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)

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

  integer :: lambda
  real :: facteur
  real, dimension(n_lambda2) :: n_photons_envoyes

  character(len=512) :: filename

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements

  logical :: simple, extend

  real, save :: L_bol1
  real :: L_bol2, E_photon1

  !! Ecriture SED
  if (ised ==1) then
     filename = trim(data_dir)//"/.sed_th.fits.gz"
  else
     filename = trim(data_dir)//"/sed_mc.fits.gz"
     if (lProDiMo) write(*,*) "WARNING: sed_mc.fits.gz will be noisy due to the use of ProDiMo mode"
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

  ! Energie d'un photon / Temps
  ! Ca donne lambda Flambda a une distance = R_etoile
  ! E_photon = sigma*T_etoile**4  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)) * real(N_thet)*real(N_phi)
  ! Ca donne lambda Flambda sur le detecteur
  if (l_sym_centrale) then
     !E_photon = L_tot  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)*(distance*pc_to_AU)**2) * real(N_thet)*real(N_phi)
     E_photon1 = L_packet_th * (real(N_thet)*real(N_phi)/quatre_pi) / (distance*pc_to_AU)**2
  else
     E_photon1 = L_packet_th * (real(2*N_thet)*real(N_phi)/quatre_pi) / (distance*pc_to_AU)**2
  endif

  if (ised == 1) then
     ! Ca donne lambda Flambda sur le detecteur

     naxis=3
     naxes(1)=n_lambda
     naxes(2)=N_thet
     naxes(3)=N_phi

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     nelements=naxes(1)*naxes(2)*naxes(3)

     do lambda=1,n_lambda
        facteur =  E_photon1 * tab_lambda(lambda)/tab_delta_lambda(lambda)
        sed1_io(lambda,:,:) = sum(sed(lambda,:,:,:),dim=3) * facteur
     enddo

     ! le e signifie real*4
     call ftppre(unit,group,fpixel,nelements,sed1_io,status)

     L_bol1 = sum(sed) * E_photon1
     !L_bol0 = L_etoile  / ((distance*pc_to_AU)**2) * real(N_thet)*real(N_phi) ! A comparer a L_bol1 OK par def
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
        sed2_io(lambda,:,:,1) = sum(sed(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,2) = sum(sed_q(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,3) = sum(sed_u(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,4) = sum(sed_v(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,5) = sum(sed_star(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,6) = sum(sed_star_scat(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,7) = sum(sed_disk(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,8) = sum(sed_disk_scat(lambda,:,:,:),dim=3) * facteur
        sed2_io(lambda,:,:,9) = sum(n_phot_sed(lambda,:,:,:),dim=3)
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
        L_bol2 = L_bol2 + sum(sed(lambda,:,:,:))*tab_delta_lambda(lambda)*1.0e-6*E_totale(lambda)
     enddo

     if ((lsed_complete).and.(ltemp)) then
        if (L_bol1 > 0.) then
           write(*,*) "Flux conservation to within ",  (L_bol2-L_bol1)/L_bol1 * 100, "%  between the two methods"
        else
           write(*,*) "Lbol1=", L_bol1
           write(*,*) "Lbol2=", L_bol2
        endif
     endif

  endif ! ised

  call ftpkys(unit,'BUNIT',"W.m-2",'lambda.F_lambda',status)

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
  call ftppre(unit,group,fpixel,nelements,real(tab_lambda,kind=sp),status)

  call ftpkys(unit,'BUNIT',"micron",' ',status)

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
  integer :: status,unit,blocksize,bitpix,naxis,icell
  integer, dimension(3) :: naxes
  integer :: group,fpixel,nelements
  logical :: simple, extend

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

  naxis=2
  naxes(1)=n_cells
  naxes(2)=nLevels

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=naxes(1)*naxes(2)

  call ftppre(unit,group,fpixel,nelements,tab_nLevel,status)

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
  integer, dimension(6) :: naxes
  integer :: group,fpixel,nelements, iv, xcenter,i
  logical :: simple, extend

  real, dimension(:,:,:), allocatable ::  O ! nv, nTrans, n_cells
  real, dimension(mol(imol)%nTrans_rayTracing) :: freq
  integer, dimension(mol(imol)%nTrans_rayTracing) :: indice_Trans

  real, dimension(:,:,:), allocatable :: spectre_casa

  real :: pixel_scale_x, pixel_scale_y, W2m2_to_Jy

  filename = trim(data_dir2(imol))//'/lines.fits.gz'

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

  if (lcasa) then
     naxis=3
     naxes(1)=npix_x
     naxes(2)=npix_y
     naxes(3)= 2*mol(imol)%n_speed_rt+1 ! velocity
     nelements=naxes(1)*naxes(2)*naxes(3)
     allocate(spectre_casa(naxes(1),naxes(2),naxes(3)))
  else
     naxis=6
     if (RT_line_method==1) then
        naxes(1)=1
        naxes(2)=1
     else
        naxes(1)=npix_x
        naxes(2)=npix_y
     endif
     naxes(3)=2*mol(imol)%n_speed_rt+1
     naxes(4)=mol(imol)%nTrans_rayTracing
     naxes(5)=RT_n_incl
     naxes(6)=RT_n_az
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
  endif

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! RAC, DEC, reference pixel & pixel scale en degres
  call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  call ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
  call ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
  pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  call ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)

  call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  call ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
  call ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
  pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
  call ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)

  call ftpkys(unit,'CTYPE3',"VELO-LSR",'[km/s]',status)
  call ftpkye(unit,'CRVAL3',0.,-7,'',status)
  call ftpkyj(unit,'CRPIX3',mol(imol)%n_speed_rt+1,'',status)
  call ftpkye(unit,'CDELT3',real(mol(imol)%vmax_center_rt/mol(imol)%n_speed_rt * m_to_km),-7,'delta_V [km/s]',status)

  if (lcasa) then
     call ftpkys(unit,'BUNIT',"JY/PIXEL",'',status)
     call ftpkyd(unit,'RESTFREQ',Transfreq(mol(imol)%indice_Trans_rayTracing(1)),-7,'',status)
     call ftpkys(unit,'BTYPE',"Intensity",' ',status)
  else
     call ftpkys(unit,'BUNIT',"W.m-2.pixel-1",'nu.F_nu',status)
  endif
  !  call ftpkye(unit,'vmax_center',mol(imol)%vmax_center_output,-8,'m/s',status)

  if (l_sym_ima) then ! BUG : ca inverse aussi la pente du Cmb mais c'est pas tres grave
     if (RT_line_method==1) then
        ! On ajoute les 2 parties du spectres
        do iv = -mol(imol)%n_speed_rt, -1
           spectre(1,1,iv,:,:,:) = spectre(1,1,iv,:,:,:) + spectre(1,1,-iv,:,:,:)
        enddo
        ! On symetrise
        do iv =1, mol(imol)%n_speed_rt
           spectre(1,1,iv,:,:,:) = spectre(1,1,-iv,:,:,:)
        enddo
        spectre(1,1,0,:,:,:) = spectre(1,1,0,:,:,:) * 2.
        ! On divise par deux
        spectre = spectre * 0.5
     else
        xcenter = npix_x/2 + modulo(npix_x,2)
        if (lkeplerian) then ! profil de raie inverse des 2 cotes
           do i=xcenter+1,npix_x
              do iv=-mol(imol)%n_speed_rt,mol(imol)%n_speed_rt
                 spectre(i,:,iv,:,:,:) = spectre(npix_x-i+1,:,-iv,:,:,:)
              enddo
           enddo
        else ! infall : meme profil de raie des 2 cotes
           do i=xcenter+1,npix_x
              spectre(i,:,:,:,:,:) = spectre(npix_x-i+1,:,:,:,:,:)
           enddo
        endif
     endif ! lkeplerian
  endif ! l_sym_image

  if (lcasa) then
     W2m2_to_Jy = 1e26 / Transfreq(mol(imol)%indice_Trans_rayTracing(1))

     spectre_casa(:,:,:) = spectre(:,:,:,1,1,1) * W2m2_to_Jy

     !  Write the array to the FITS file.
     call ftppre(unit,group,fpixel,nelements,spectre_casa,status)
     deallocate(spectre_casa)
  else
     !  Write the array to the FITS file.
     call ftppre(unit,group,fpixel,nelements,spectre,status)
  endif

  if (.not.lcasa) then
     !------------------------------------------------------------------------------
     ! HDU 2 : Continuum map
     !------------------------------------------------------------------------------
     bitpix=-32
     naxis=5
     if (RT_line_method==1) then
        naxes(1)=1
        naxes(2)=1
     else
        naxes(1)=npix_x
        naxes(2)=npix_y
     endif
     naxes(3)=mol(imol)%nTrans_rayTracing
     naxes(4)=RT_n_incl
     naxes(5)=RT_n_az
     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)

     ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     if (l_sym_ima.and.(RT_line_method==2)) then
        xcenter = npix_x/2 + modulo(npix_x,2)
        do i=xcenter+1,npix_x
           continu(i,:,:,:,:) = continu(npix_x-i+1,:,:,:,:)
        enddo
     endif ! l_sym_image

     !  Write the array to the FITS file.
     call ftppre(unit,group,fpixel,nelements,continu,status)

     !------------------------------------------------------------------------------
     ! HDU 3 : Transition numbers
     !------------------------------------------------------------------------------
     bitpix=32
     naxis = 1
     naxes(1) = mol(imol)%nTrans_rayTracing
     nelements = naxes(1)

     ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     do i=1, mol(imol)%nTrans_rayTracing
        indice_Trans(i) = mol(imol)%indice_Trans_rayTracing(i)
     enddo

     !  Write the array to the FITS file.
     call ftpprj(unit,group,fpixel,nelements,indice_Trans,status)

     !------------------------------------------------------------------------------
     ! HDU 4 : Transition frequencies
     !------------------------------------------------------------------------------
     bitpix=-32
     naxis = 1
     naxes(1) = mol(imol)%nTrans_rayTracing
     nelements = naxes(1)

     ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     do i=1,mol(imol)%nTrans_rayTracing
        freq(i) = Transfreq(mol(imol)%indice_Trans_rayTracing(i))
     enddo

     !  Write the array to the FITS file.
     call ftppre(unit,group,fpixel,nelements,freq,status)

     !------------------------------------------------------------------------------
     ! HDU 5 : Velocities
     !------------------------------------------------------------------------------
     bitpix=-32
     naxis = 1
     naxes(1) = 2*mol(imol)%n_speed_rt+1
     nelements = naxes(1)

     ! create new hdu
     call ftcrhd(unit, status)

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     call ftppre(unit,group,fpixel,nelements,real(tab_speed_rt),status)
  endif ! lcasa

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
     allocate(O(-mol(imol)%n_speed_rt:mol(imol)%n_speed_rt,mol(imol)%nTrans_rayTracing,n_cells))
     O = 0.0;
     do i=1, nb_proc
        O(:,:,:) =  O(:,:,:) +  origine_mol(:,:,:,i)
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
     naxes(1)=2*mol(imol)%n_speed_rt+1
     naxes(2)=mol(imol)%nTrans_rayTracing
     naxes(3)=n_cells

     !  Write the required header keywords.
     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

     !  Write the array to the FITS file.
     group=1
     fpixel=1
     nelements=naxes(1)*naxes(2)*naxes(3)

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

subroutine write_temperature_for_phantom(n_SPH)

  integer, intent(in) :: n_SPH
  integer :: i_SPH, icell
  real, dimension(n_SPH) :: T_SPH

  T_SPH = -1.0 ;

  do icell=1, n_cells
     i_SPH = Voronoi(icell)%id
     if (i_SPH > 0) T_SPH(i_SPH) = Tdust(icell)
  enddo

  open(1,file="T_for_phantom.tmp",status='replace',form='unformatted')
  write(1) T_SPH
  close(unit=1)

  return

end subroutine write_temperature_for_phantom

!**********************************************************************

end module output
