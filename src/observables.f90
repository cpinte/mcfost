module observables

  use disk
  use ki2
  use grains
  use resultats

  implicit none

  contains

subroutine crea_noyau()
  ! Cree le tableau contenant l'intensite + des 0 autour
  ! ainsi que le noyau de convolution
  ! C. pinte   18/11/04

  implicit none

  integer :: i,j, alloc_status
  real :: coeff, somme, fwhm

  ! FWHM en arc. sec.
  if (abs(tab_lambda(1) - 0.81) < 1.0e-3) then ! I
     fwhm=0.0737
  else if (abs(tab_lambda(1) - 1.6) < 1.0e-3) then ! H
     fwhm=0.0960
  else if (abs(tab_lambda(1) - 2.2) < 1.0e-3) then ! K
     fwhm=0.0719
  else if (abs(tab_lambda(1) - 3.8) < 1.0e-3) then ! L
     fwhm=0.0812
  else
     write(*,*) 'WARNING : fwhm not defined, using 0.1"'
     fwhm=0.1
  endif
  sigma=fwhm/(2.0*sqrt(2.0*log(2.0))) ! arc. sec.
  sigma=sigma*distance*igridx/(2*size_neb)  ! pixels

  ! Noyau de convolution  
  n_ker=4*ceiling(sigma)

  allocate(noyau(-n_ker:n_ker,-n_ker:n_ker), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  noyau=0.0

  coeff=-0.5/sigma**2
  somme=0.0
  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        noyau(i,j)=exp(coeff*(real(i)**2+real(j)**2))
        somme=somme+noyau(i,j)
     enddo !j
  enddo !i

  ! Normalisation
  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        noyau(i,j)=noyau(i,j)/somme
     enddo !j
  enddo !i

  ! Fichier avec image du disque + 0
  allocate(intensite(1-n_ker:igridx+n_ker,1-n_ker:igridy+n_ker), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  intensite=0.0

end subroutine crea_noyau

!********************************************************************

subroutine crea_noyau_2mass()
  ! Cree le tableau contenant l'intensite + des 0 autour
  ! ainsi que le noyau de convolution
  ! C. pinte   21/06/05

  implicit none

  integer :: i,j, alloc_status
  real :: coeff, somme, fwhm

  ! FWHM en arc. sec.
  if (abs(tab_lambda(1) - 0.81) < 1.0e-3) then ! I
     fwhm=0.40
  else if (abs(tab_lambda(1) - 2.2) < 1.0e-3) then ! K
     fwhm=0.42
  else if (abs(tab_lambda(1) - 3.8) < 1.0e-3) then ! L
     fwhm=0.42
  else
     write(*,*) 'WARNING : fwhm not defined, using 0.1"'
     fwhm=0.1
  endif
  sigma=fwhm/(2.0*sqrt(2.0*log(2.0))) ! arc. sec.
  sigma=sigma*distance*igridx/(2*size_neb)  ! pixels

  ! Noyau de convolution  
  n_ker=4*ceiling(sigma)

  allocate(noyau(-n_ker:n_ker,-n_ker:n_ker), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  noyau=0.0

  coeff=-0.5/sigma**2
  somme=0.0
  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        noyau(i,j)=exp(coeff*(real(i)**2+real(j)**2))
        somme=somme+noyau(i,j)
     enddo !j
  enddo !i

  ! Normalisation
  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        noyau(i,j)=noyau(i,j)/somme
     enddo !j
  enddo !i

  ! Fichier avec image du disque + 0
  allocate(intensite(1-n_ker:igridx+n_ker,1-n_ker:igridy+n_ker), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  intensite=0.0

end subroutine crea_noyau_2mass

!********************************************************************

subroutine crea_noyau_PDS144()
  ! Cree le tableau contenant l'intensite + des 0 autour
  ! ainsi que le noyau de convolution
  ! C. pinte   18/11/04

  implicit none

  integer :: i,j, ii, jj, alloc_status
  real :: coeff, somme, fwhm

  ! Noyau de convolution  
  n_ker=32

  allocate(noyau(-n_ker:n_ker,-n_ker:n_ker), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  noyau=0.0

  open(unit=1,file="psf_h.data",status="old")

  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        read(1,*) ii,jj, noyau(i,j)
        somme=somme+noyau(i,j)
     enddo !j
  enddo !i

  close(unit=1)

  ! Normalisation
  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        noyau(i,j)=noyau(i,j)/somme
     enddo !j
  enddo !i

  ! Fichier avec image du disque + 0
  allocate(intensite(1-n_ker:igridx+n_ker,1-n_ker:igridy+n_ker), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  intensite=0.0

  return

end subroutine crea_noyau_PDS144

!********************************************************************

subroutine crea_image(capt)
  ! Remplit le tableau image
  ! C. pinte   18/11/04

  implicit none

  integer, intent(in) :: capt
  integer :: i,j

  do i=1,igridx
     do j=1,igridy
        intensite(i,j)= sum(stokei(:,1,i,j,capt,1))
     enddo !j
  enddo !j

end subroutine crea_image

!********************************************************************

real function convol(x,y)
  ! Convolution en x,y
  ! C. pinte   18/11/04

  implicit none

  integer, intent(in) :: x,y

  integer :: i,j

  convol=0.0
  do i=-n_ker,n_ker
     do j=-n_ker,n_ker
        convol=convol+intensite(x+i,y+j)*noyau(i,j)
     enddo!j
  enddo !j


  return

end function convol

!********************************************************************

subroutine observables_HK_tau()
  ! Calcule des ki2 sur une des images donnée par le MC
  ! Convolution partielle sur les pixels intéressants.

  implicit none

  integer :: i,j,k, kmax, kmin,x0,y0,ysup,yinf,capt,ymin
  real :: ftot,sup,inf, pix_conv, max_sup, max_inf, conv_sec, dist
  real :: rap_flux, seuil, correct, delta_pix,fwhm_h,value,vmax,vmin, pix_conv_p1
  character(len=512) :: filename
  real, dimension(igridx,igridy) :: image_conv

  real :: Q_tot, U_tot, pola, max_1_2, rap_flux_pic, hwhm_v1, hwhm_v2, rap_flux_1_2, angle
  integer :: x_1_2, y_1_2

  x0 = int(igridx/2.0)+1
  y0 = int(igridy/2.0)+1
  conv_sec=(2*size_neb/igridy)/distance

  filename=trim(data_dir)//'/obs.txt'
100 format(I3,3X,10(ES10.3,3X))
  open(unit=10,file=filename, status='new')

  do capt=capt_inf,capt_sup

     call crea_image(capt)

     do i=1,igridx
        do j=1,igridy
           image_conv(i,j)=convol(i,j)
        enddo !j
     enddo !i

     ! Recherche max nebuleuse brillante (convoluee)
     max_sup=0.0
     do j = y0+1,igridy!y0+(igridy-1)/10
        pix_conv = image_conv(x0,j)
        if (pix_conv > max_sup) then
           max_sup = pix_conv
           ysup=j
        endif
     enddo  !j
     ! Recherche min entre nebuleuses (convoluee)  
     pix_conv_p1=max_sup
     min2 :do j=ysup-1, 1, -1
        pix_conv = image_conv(x0,j)
        if (pix_conv > pix_conv_p1) then
           ymin = j+1
           exit min2
        endif
        pix_conv_p1=pix_conv
     enddo min2
     ymin=ymin-1
     ! Recherche max nebuleuse faible (convoluee)
     max_inf=0.0
     do j=ymin,1,-1
        pix_conv = image_conv(x0,j)
        if (pix_conv > max_inf) then
           max_inf = pix_conv
           yinf=j
        endif
     enddo

     rap_flux_pic=max_sup/max_inf
     dist=(ysup-yinf)*conv_sec


     ! Rapport des Flux des 2 nébuleuses et flux total
     ! sur image non convoluee, limite fixee par bande sombre sur image convoluée
     inf=0.0
     sup=0.0
     ftot=0.0
     do i=1,igridx
        do j=1,igridy
           if (j>ymin) then 
              sup=sup+intensite(i,j)
           else if (j<ymin) then
              inf=inf+intensite(i,j)
           else
              ftot=ftot+intensite(i,j)
           endif
        enddo !j
     enddo !j
     ftot=ftot+sup+inf
     if (inf > 0.) then
        rap_flux=sup/inf
     else
        rap_flux = -1.e9
     endif


     ! FWHM horizontale
     seuil=0.5*max_sup

     fh : do i = x0,igridx
        if (image_conv(i,ysup) < seuil) then
           k=i
           exit fh
        endif
     enddo fh  !j

     correct=(seuil-image_conv(k,ysup))/(image_conv(k-1,ysup)-image_conv(k,ysup))
     delta_pix=real(k) - correct - real(x0)

     fwhm_h=2.0*delta_pix*conv_sec


     ! FWHM verticales
     fv1 : do j = ysup,igridy
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv1
        endif
     enddo fv1  !j

     correct=(seuil-image_conv(x0,k))/(image_conv(x0,k-1)-image_conv(x0,k))
     delta_pix=real(k) - real(ysup) - correct
     hwhm_v1=delta_pix*conv_sec

     fv2 : do j = ysup,1,-1
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv2
        endif
     enddo fv2  !j

     correct=(seuil-image_conv(x0,k))/(image_conv(x0,k+1)-image_conv(x0,k))
     delta_pix=real(ysup) - real(k) - correct
     hwhm_v2=delta_pix*conv_sec


     ! Recherche max nebuleuse brillante a mi rayon(convoluee)
     x_1_2=20
     max_1_2 = 0.0
     do j = y0+1,igridy
        pix_conv = image_conv(x_1_2,j)
        if (pix_conv > max_1_2) then
           max_1_2 = pix_conv
           y_1_2=j
        endif
     enddo  !j
     angle=57.29578*atan(real(y_1_2-ysup)/31.)
     rap_flux_1_2=max_1_2/max_sup


     ! Pola integrée
     Q_tot=0.0
     U_tot=0.0
     do i=1,igridx
        do j=1,igridy
           Q_tot=Q_tot+sum(STOKEQ(:,1,i,j,capt,1))
           U_tot=U_tot+sum(STOKEU(:,1,i,j,capt,1))
        enddo !j
     enddo !i
     pola=sqrt(Q_tot*Q_tot+U_tot*U_tot)/ftot

     !    write(*,*) ' '
     !    write(*,*) '************************************************************'
     !    write(*,*) "Capteur : ",capt
     !    write(*,*) 'Flux integres'
     !    write(*,*) '******************************'
     !    write(*,*) 'Flux total =', ftot
     !    write(*,*) 'sup/inf = ', rap_flux
     !    write(*,*) ' '
     !    write(*,*) 'Pic a pic'
     !    write(*,*) '******************************'
     !    write(*,*) 'sup/inf = ', max_sup/max_inf
     !    write(*,*) 'Distance (sec. arc) ', dist
     !    write(*,*) ' '
     !    write(*,*) 'Profil horizontal'
     !    write(*,*) '******************************'
     !    write(*,*) 'FWHM  (sec. arc) = ', fwhm 

     write(10,100) capt,ftot, rap_flux, rap_flux_pic, dist, fwhm_h, hwhm_v1, hwhm_v2, angle, rap_flux_1_2, pola

  enddo!capt

  close(unit=10)

  deallocate(noyau)
  deallocate(intensite)

  return

end subroutine observables_HK_tau

!********************************************************************

subroutine observables_PDS144()
  ! Calcule des ki2 sur une des images donnée par le MC
  ! Convolution partielle sur les pixels intéressants.

  implicit none

  integer :: i,j,k, kmax, kmin,x0,y0,ysup,yinf,capt,ymin
  real :: ftot,sup,inf, pix_conv, max_sup, max_inf, conv_sec, dist
  real :: rap_flux, seuil, correct, delta_pix,fwhm_h,value,vmax,vmin, pix_conv_p1
  real :: inf_conv, sup_conv, rap_flux_conv, ftot_conv, mini
  character(len=512) :: filename
  real, dimension(igridx,igridy) :: image_conv

  real :: Q_tot, U_tot, pola, max_1_2, rap_flux_pic, hwhm_v1, hwhm_v2, rap_flux_1_2, angle, w_01_inf, w_01_sup
  integer :: x_1_2, y_1_2

  x0 = int(igridx/2.0)+1
  y0 = int(igridy/2.0)+1
  conv_sec=(2*size_neb/igridy)/distance

  filename=trim(data_dir)//'/obs.txt'
100 format(I3,3X,13(ES10.3,3X))
  open(unit=10,file=filename, status='new')

  do capt=capt_inf,capt_sup

     call crea_image(capt)

     do i=1,igridx
        do j=1,igridy
           image_conv(i,j)=convol(i,j)
        enddo !j
     enddo !i

     ! Recherche max nebuleuse brillante (convoluee)
     max_sup=0.0
     do j = igridy, 1, -1
        pix_conv = image_conv(x0,j)
        if (pix_conv > max_sup) then
           max_sup = pix_conv
           ysup=j
        endif
     enddo  !j
     ! Recherche min entre nebuleuses (convoluee)  
     pix_conv_p1=max_sup
     min2 :do j=ysup-1, 1, -1
        pix_conv = image_conv(x0,j)
        if (pix_conv > pix_conv_p1) then
           ymin = j+1
           exit min2
        endif
        pix_conv_p1=pix_conv
     enddo min2
     ymin=ymin-1
     ! Recherche max nebuleuse faible (convoluee)
     max_inf=0.0
     do j=ymin,1,-1
        pix_conv = image_conv(x0,j)
        if (pix_conv > max_inf) then
           max_inf = pix_conv
           yinf=j
        endif
     enddo
     ! ReRecherche min entre nebuleuses (convoluee)  
     mini=max_sup
     do j=ysup,yinf,-1
        pix_conv = image_conv(x0,j)
        if (pix_conv < mini) then
           mini = pix_conv
           ymin=j
        endif
     enddo


     rap_flux_pic=max_sup/max_inf
     dist=(ysup-yinf)*conv_sec


     ! Rapport des Flux des 2 nébuleuses et flux total
     ! sur image non convoluee, limite fixee par bande sombre sur image convoluée
     inf=0.0
     sup=0.0
     ftot=0.0
     do i=1,igridx
        do j=1,igridy
           if (j>ymin) then 
              sup=sup+intensite(i,j)
           else if (j<ymin) then
              inf=inf+intensite(i,j)
           else
              ftot=ftot+intensite(i,j)
           endif
        enddo !j
     enddo !j
     ftot=ftot+sup+inf
     if (inf > tiny(0.0)) then
        rap_flux=sup/inf
     else
        rap_flux = 1.0e6
     endif

     ! Rapport des Flux des 2 nébuleuses et flux total
     ! sur image convoluee, limite fixee par bande sombre sur image convoluée
     inf_conv=0.0
     sup_conv=0.0
     ftot_conv=0.0
     do i=1,igridx
        do j=1,igridy
           if (j>ymin) then 
              sup_conv=sup_conv+image_conv(i,j)
           else if (j<ymin) then
              inf_conv=inf_conv+image_conv(i,j)
           else
              ftot_conv=ftot_conv+image_conv(i,j)
           endif
        enddo !j
     enddo !j
     ftot_conv=ftot_conv+sup_conv+inf_conv
     if (inf_conv > tiny(0.0)) then
        rap_flux_conv=sup_conv/inf_conv
     else
        rap_flux_conv = 1.0e6
     endif


     ! FWHM horizontale
     seuil=0.5*max_sup

     fh : do i = x0,igridx
        if (image_conv(i,ysup) < seuil) then
           k=i
           exit fh
        endif
     enddo fh  !j

     if ((image_conv(k-1,ysup)-image_conv(k,ysup)) /= 0.0) then
        correct=(seuil-image_conv(k,ysup))/(image_conv(k-1,ysup)-image_conv(k,ysup))
     else
        correct=0.0
     endif
     delta_pix=real(k) - correct - real(x0)

     fwhm_h=2.0*delta_pix*conv_sec


     ! FWHM verticales
     fv1 : do j = ysup,igridy
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv1
        endif
     enddo fv1  !j

     if ((image_conv(x0,k-1)-image_conv(x0,k)) /= 0.0) then
        correct=(seuil-image_conv(x0,k))/(image_conv(x0,k-1)-image_conv(x0,k))
     else
        correct=0.0
     endif
     delta_pix=real(k) - real(ysup) - correct
     hwhm_v1=delta_pix*conv_sec

     fv2 : do j = ysup,1,-1
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv2
        endif
     enddo fv2  !j

     if ((image_conv(x0,k+1)-image_conv(x0,k)) /= 0.0) then
        correct=(seuil-image_conv(x0,k))/(image_conv(x0,k+1)-image_conv(x0,k))
     else
        correct=0.0
     endif
     delta_pix=real(ysup) - real(k) - correct
     hwhm_v2=delta_pix*conv_sec


     ! Recherche max nebuleuse brillante a une certaine distance x_1_2 du bord (convoluee)
     x_1_2=49
     max_1_2 = 0.0
     do j = y0+1,igridy
        pix_conv = image_conv(x_1_2,j)
        if (pix_conv > max_1_2) then
           max_1_2 = pix_conv
           y_1_2=j
        endif
     enddo  !j
     angle=57.29578*atan(real(y_1_2-ysup)/(x0-x_1_2))
     rap_flux_1_2=max_1_2/max_sup


     ! Pola integrée
     Q_tot=0.0
     U_tot=0.0
     do i=1,igridx
        do j=1,igridy
           Q_tot=Q_tot+sum(STOKEQ(:,1,i,j,capt,1))
           U_tot=U_tot+sum(STOKEU(:,1,i,j,capt,1))
        enddo !j
     enddo !i
     pola=sqrt(Q_tot*Q_tot+U_tot*U_tot)/ftot

     ! Flux=10% du max sup
     seuil=0.1*max_sup

     fh2 : do i = x0,igridx
        if (image_conv(i,ysup) < seuil) then
           k=i
           exit fh2
        endif
     enddo fh2  !j

     correct=(seuil-image_conv(k,ysup))/(image_conv(k-1,ysup)-image_conv(k,ysup))
     delta_pix=real(k) - correct - real(x0)

     w_01_sup=delta_pix*conv_sec

     ! Flux=10% du max inf
     seuil=0.1*max_inf

     fh3 : do i = x0,igridx
        if (image_conv(i,yinf) < seuil) then
           k=i
           exit fh3
        endif
     enddo fh3  !j

     correct=(seuil-image_conv(k,yinf))/(image_conv(k-1,yinf)-image_conv(k,yinf))
     delta_pix=real(k) - correct - real(x0)

     w_01_inf=delta_pix*conv_sec


     !    write(*,*) ' '
     !    write(*,*) '************************************************************'
     !    write(*,*) "Capteur : ",capt
     !    write(*,*) 'Flux integres'
     !    write(*,*) '******************************'
     !    write(*,*) 'Flux total =', ftot
     !    write(*,*) 'sup/inf = ', rap_flux
     !    write(*,*) ' '
     !    write(*,*) 'Pic a pic'
     !    write(*,*) '******************************'
     !    write(*,*) 'sup/inf = ', max_sup/max_inf
     !    write(*,*) 'Distance (sec. arc) ', dist
     !    write(*,*) ' '
     !    write(*,*) 'Profil horizontal'
     !    write(*,*) '******************************'
     !    write(*,*) 'FWHM  (sec. arc) = ', fwhm 

     write(10,100) capt,ftot, rap_flux, rap_flux_conv, rap_flux_pic, dist, fwhm_h, hwhm_v1, hwhm_v2, &
          angle, rap_flux_1_2, w_01_sup,  w_01_inf, pola

  enddo!capt

  close(unit=10)

  deallocate(noyau)
  deallocate(intensite)

  return

end subroutine observables_PDS144

!********************************************************************

subroutine observables_2mass()
  ! Calcule des ki2 sur une des images donnée par le MC
  ! Convolution partielle sur les pixels intéressants.
  ! C . Pinte 21/06/05

  implicit none

  integer :: i,j,k, kmax, kmin,x0,y0,ysup,yinf,capt,ymin
  real :: ftot,sup,inf, pix_conv, max_sup, max_inf, conv_sec, dist
  real :: rap_flux, seuil, correct, delta_pix,fwhm_h,value,vmax,vmin, pix_conv_p1
  real :: inf_conv, sup_conv, rap_flux_conv, ftot_conv, mini
  character(len=512) :: filename
  real, dimension(igridx,igridy) :: image_conv

  real :: Q_tot, U_tot, pola, max_1_2, rap_flux_pic, hwhm_v1, hwhm_v2, rap_flux_1_2, angle, w_01_inf, w_01_sup
  real :: fwhm_h_inf, hwhm_v1_inf, hwhm_v2_inf
  integer :: x_1_2, y_1_2

  x0 = int(igridx/2.0)+1
  y0 = int(igridy/2.0)+1
  conv_sec=(2*size_neb/igridy)/distance

  filename=trim(data_dir)//'/obs.txt'
100 format(I3,3X,16(ES10.3,3X))
  open(unit=10,file=filename, status='new')

  do capt=capt_inf,capt_sup

     call crea_image(capt)

     do i=1,igridx
        do j=1,igridy
           image_conv(i,j)=convol(i,j)
        enddo !j
     enddo !i

     ! Recherche max nebuleuse brillante (convoluee)
     max_sup=0.0
     do j = igridy, 1, -1
        pix_conv = image_conv(x0,j)
        if (pix_conv > max_sup) then
           max_sup = pix_conv
           ysup=j
        endif
     enddo  !j
     ! Recherche min entre nebuleuses (convoluee)  
     pix_conv_p1=max_sup
     min2 :do j=ysup-1, 1, -1
        pix_conv = image_conv(x0,j)
        if (pix_conv > pix_conv_p1) then
           ymin = j+1
           exit min2
        endif
        pix_conv_p1=pix_conv
     enddo min2
     ymin=ymin-1
     ! Recherche max nebuleuse faible (convoluee)
     max_inf=0.0
     do j=ymin,1,-1
        pix_conv = image_conv(x0,j)
        if (pix_conv > max_inf) then
           max_inf = pix_conv
           yinf=j
        endif
     enddo
     ! ReRecherche min entre nebuleuses (convoluee)  
     mini=max_sup
     do j=ysup,yinf,-1
        pix_conv = image_conv(x0,j)
        if (pix_conv < mini) then
           mini = pix_conv
           ymin=j
        endif
     enddo


     rap_flux_pic=max_sup/max_inf
     dist=(ysup-yinf)*conv_sec


     ! Rapport des Flux des 2 nébuleuses et flux total
     ! sur image non convoluee, limite fixee par bande sombre sur image convoluée
     inf=0.0
     sup=0.0
     ftot=0.0
     do i=1,igridx
        do j=1,igridy
           if (j>ymin) then 
              sup=sup+intensite(i,j)
           else if (j<ymin) then
              inf=inf+intensite(i,j)
           else
              ftot=ftot+intensite(i,j)
           endif
        enddo !j
     enddo !j
     ftot=ftot+sup+inf
     if (inf > tiny(0.0)) then
        rap_flux=sup/inf
     else
        rap_flux = 1.0e6
     endif

     ! Rapport des Flux des 2 nébuleuses et flux total
     ! sur image convoluee, limite fixee par bande sombre sur image convoluée
     inf_conv=0.0
     sup_conv=0.0
     ftot_conv=0.0
     do i=1,igridx
        do j=1,igridy
           if (j>ymin) then 
              sup_conv=sup_conv+image_conv(i,j)
           else if (j<ymin) then
              inf_conv=inf_conv+image_conv(i,j)
           else
              ftot_conv=ftot_conv+image_conv(i,j)
           endif
        enddo !j
     enddo !j
     ftot_conv=ftot_conv+sup_conv+inf_conv
     if (inf_conv > tiny(0.0)) then
        rap_flux_conv=sup_conv/inf_conv
     else
        rap_flux_conv = 1.0e6
     endif


     ! FWHM horizontale
     seuil=0.5*max_sup

     fh : do i = x0,igridx
        if (image_conv(i,ysup) < seuil) then
           k=i
           exit fh
        endif
     enddo fh  !j

     if ((image_conv(k-1,ysup)-image_conv(k,ysup)) /= 0.0) then
        correct=(seuil-image_conv(k,ysup))/(image_conv(k-1,ysup)-image_conv(k,ysup))
     else
        correct=0.0
     endif
     delta_pix=real(k) - correct - real(x0)

     fwhm_h=2.0*delta_pix*conv_sec


     ! FWHM verticales
     fv1 : do j = ysup,igridy
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv1
        endif
     enddo fv1  !j

     if ((image_conv(x0,k-1)-image_conv(x0,k)) /= 0.0) then
        correct=(seuil-image_conv(x0,k))/(image_conv(x0,k-1)-image_conv(x0,k))
     else
        correct=0.0
     endif
     delta_pix=real(k) - real(ysup) - correct
     hwhm_v1=delta_pix*conv_sec

     fv2 : do j = ysup,1,-1
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv2
        endif
     enddo fv2  !j

     if ((image_conv(x0,k+1)-image_conv(x0,k)) /= 0.0) then
        correct=(seuil-image_conv(x0,k))/(image_conv(x0,k+1)-image_conv(x0,k))
     else
        correct=0.0
     endif
     delta_pix=real(ysup) - real(k) - correct
     hwhm_v2=delta_pix*conv_sec

     ! FWHM horizontale (neb inf)
     seuil=0.5*max_inf

     fh_inf : do i = x0,igridx
        if (image_conv(i,yinf) < seuil) then
           k=i
           exit fh_inf
        endif
     enddo fh_inf  !j

     if ((image_conv(k-1,yinf)-image_conv(k,yinf)) /= 0.0) then
        correct=(seuil-image_conv(k,yinf))/(image_conv(k-1,yinf)-image_conv(k,yinf))
     else
        correct=0.0
     endif
     delta_pix=real(k) - correct - real(x0)

     fwhm_h_inf=2.0*delta_pix*conv_sec


     ! FWHM verticales
     fv1_inf : do j = yinf,igridy
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv1_inf
        endif
     enddo fv1_inf  !j

     if ((image_conv(x0,k-1)-image_conv(x0,k)) /= 0.0) then
        correct=(seuil-image_conv(x0,k))/(image_conv(x0,k-1)-image_conv(x0,k))
     else
        correct=0.0
     endif
     delta_pix=real(k) - real(ysup) - correct
     hwhm_v1_inf=delta_pix*conv_sec

     fv2_inf : do j = yinf,1,-1
        if (image_conv(x0,j) < seuil) then
           k=j
           exit fv2_inf
        endif
     enddo fv2_inf  !j

     if ((image_conv(x0,k+1)-image_conv(x0,k)) /= 0.0) then
        correct=(seuil-image_conv(x0,k))/(image_conv(x0,k+1)-image_conv(x0,k))
     else
        correct=0.0
     endif
     delta_pix=real(ysup) - real(k) - correct
     hwhm_v2_inf=delta_pix*conv_sec


     ! Recherche max nebuleuse brillante a une certaine distance x_1_2 du bord (convoluee)
     x_1_2=20
     max_1_2 = 0.0
     do j = y0+1,igridy
        pix_conv = image_conv(x_1_2,j)
        if (pix_conv > max_1_2) then
           max_1_2 = pix_conv
           y_1_2=j
        endif
     enddo  !j
     angle=57.29578*atan(real(y_1_2-ysup)/(x0-x_1_2))
     rap_flux_1_2=max_1_2/max_sup


     ! Pola integrée
     Q_tot=0.0
     U_tot=0.0
     do i=1,igridx
        do j=1,igridy
           Q_tot=Q_tot+sum(STOKEQ(:,1,i,j,capt,1))
           U_tot=U_tot+sum(STOKEU(:,1,i,j,capt,1))
        enddo !j
     enddo !i
     pola=sqrt(Q_tot*Q_tot+U_tot*U_tot)/ftot

     ! Flux=10% du max sup
     seuil=0.1*max_sup

     fh2 : do i = x0,igridx
        if (image_conv(i,ysup) < seuil) then
           k=i
           exit fh2
        endif
     enddo fh2  !j

     correct=(seuil-image_conv(k,ysup))/(image_conv(k-1,ysup)-image_conv(k,ysup))
     delta_pix=real(k) - correct - real(x0)

     w_01_sup=delta_pix*conv_sec

     ! Flux=10% du max inf
     seuil=0.1*max_inf

     fh3 : do i = x0,igridx
        if (image_conv(i,yinf) < seuil) then
           k=i
           exit fh3
        endif
     enddo fh3  !j

     correct=(seuil-image_conv(k,yinf))/(image_conv(k-1,yinf)-image_conv(k,yinf))
     delta_pix=real(k) - correct - real(x0)

     w_01_inf=delta_pix*conv_sec


     !    write(*,*) ' '
     !    write(*,*) '************************************************************'
     !    write(*,*) "Capteur : ",capt
     !    write(*,*) 'Flux integres'
     !    write(*,*) '******************************'
     !    write(*,*) 'Flux total =', ftot
     !    write(*,*) 'sup/inf = ', rap_flux
     !    write(*,*) ' '
     !    write(*,*) 'Pic a pic'
     !    write(*,*) '******************************'
     !    write(*,*) 'sup/inf = ', max_sup/max_inf
     !    write(*,*) 'Distance (sec. arc) ', dist
     !    write(*,*) ' '
     !    write(*,*) 'Profil horizontal'
     !    write(*,*) '******************************'
     !    write(*,*) 'FWHM  (sec. arc) = ', fwhm 

     write(10,100) capt,ftot, rap_flux, rap_flux_conv, rap_flux_pic, dist, fwhm_h, hwhm_v1, hwhm_v2, &
          fwhm_h_inf, hwhm_v1_inf, hwhm_v2_inf, angle, rap_flux_1_2, w_01_sup,  w_01_inf, pola

  enddo!capt

  close(unit=10)

  deallocate(noyau)
  deallocate(intensite)

  return

end subroutine observables_2mass

!********************************************************************

end module observables
