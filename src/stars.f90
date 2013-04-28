module stars

  use parametres
  use disk
  use prop_star
  use grains
  use em_th
  use constantes

  use scattering
  use ProDiMo
  use optical_depth

  implicit none

  contains

subroutine select_etoile(lambda,aleat,n_star)
! Sélection d'étoile qui va émettre le photon
! C. Pinte
! 21/05/05

  implicit none

  integer, intent(in) :: lambda
  real, intent(in) :: aleat
  integer, intent(out) :: n_star
  integer :: k, kmin, kmax

  ! Dichotomie
  kmin=0
  kmax=n_etoiles
  k=(kmax-kmin)/2

  do while ((kmax-kmin) > 1)
     if (prob_E_star(lambda,k) < aleat) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   n_star=kmax

   return

end subroutine select_etoile

!**********************************************************************

subroutine em_sphere_uniforme(n_star,aleat1,aleat2,aleat3,aleat4,ri,zj,x,y,z,u,v,w,w2)
! Choisit la position d'emission uniformement
! sur la surface de l'etoile et la direction de vol
! suivant le cos de l'angle / normale
! C. Pinte
! 21/05/05

  implicit none

  integer, intent(in) :: n_star
  real, intent(in) :: aleat1, aleat2, aleat3, aleat4
  integer, intent(out) :: ri, zj
  real(kind=db), intent(out) :: x, y, z, u, v, w, w2

  real(kind=db) :: srw02, argmt, r_etoile, cospsi, phi

  ! Position de depart aleatoire sur une sphere de rayon 1
  z = 2.0_db * aleat1 - 1.0_db
  srw02 = sqrt(1.0-z*z)
  argmt = pi*(2.0_db*aleat2-1.0_db)
  x = srw02 * cos(argmt)
  y = srw02 * sin(argmt)


  ! Choix direction de vol : sphere uniforme
  cospsi = sqrt(aleat3) ! !TODO : c'est bizarre ca. aleat3 marche avec le benchmark : mais pas la meme chose pour la surface stellaire dans RT. sqrt(aleat3) ~ OK avec TORUS en mode MC
  phi = 2.0_db*pi*aleat4
  ! (x,y,z) definit la normale (ici, il a encore une norme=1)
  ! cospsi et phi sont definis / cette normale
  ! il faut faire une rotation
  call cdapres(cospsi, phi, x, y, z, u, v, w)

  w2=1.0_db-w*w

  ! Position de depart aleatoire sur une sphere de rayon r_etoile
  r_etoile = etoile(n_star)%r
  x = x * r_etoile
  y = y * r_etoile
  z = z * r_etoile

  ! Ajout position de l'étoile
  x=x+etoile(n_star)%x
  y=y+etoile(n_star)%y
  z=z+etoile(n_star)%z

  ri=0 ; zj=1

  return

end subroutine em_sphere_uniforme

!**********************************************************************

subroutine em_etoile_ponctuelle(n_star,aleat1,aleat2,ri,zj,x,y,z,u,v,w,w2)
! Emission isotrope
! C. Pinte
! 21/05/05

  implicit none

  integer, intent(in) :: n_star
  real, intent(in) :: aleat1, aleat2
  integer, intent(out) :: ri, zj
  real(kind=db), intent(out) :: x, y, z, u, v, w, w2

  real(kind=db) :: srw02, argmt

  ! Emission isotrope
  w = 2.0_db * aleat1 - 1.0_db
  w2 = 1.0_db-w*w
  srw02 = sqrt(w2)
  argmt = pi*(2.0_db*aleat2-1.0_db)
  u = srw02 * cos(argmt)
  v = srw02 * sin(argmt)

  ! Position de l'étoile
  x=etoile(n_star)%x
  y=etoile(n_star)%y
  z=etoile(n_star)%z

  ri=0 ; zj=1

  return

end subroutine em_etoile_ponctuelle

!**********************************************************************

subroutine loi_Wien(lambda)

  implicit none

  integer, intent(out) :: lambda
  real :: T
  real :: wl

!  T = maxval(etoile%T)

  wl = 2898./T

  loop : do lambda=1, n_lambda
     if (tab_lambda(lambda) > wl) exit loop
  enddo loop
  lambda = lambda - 1

end subroutine loi_Wien

!***********************************************************

subroutine repartition_energie_etoiles()
! Calcul l'energie totale emise par les etoiles
! et la proba d'emission de chaque etoile pout tous les lambda
! Pour toutes versions du code
! C. Pinte
! 21/05/05

! Modif : 26/07/05 (C. Pinte)
! Fusion avec create_b_body : la routine cree aussi spectre_etoile
! la proba d'emettre a lambda pour toutes les etoiles
! (Calculs redondants)

  implicit none

  integer :: lambda, i, n, l1, l2, lambda_mini, lambda_maxi, n_lambda_spectre, l
  real(kind=db) :: wl, cst_wl, delta_wl, surface, terme, terme0, spectre, spectre0
  real(kind=db) :: somme_spectre, somme_bb
  real ::  wl_inf, wl_sup, UV_ProDiMo, p, cst_UV_ProDiMo
  real(kind=db) :: fact_sup, fact_inf, cst_spectre_etoiles, L_save, correct, correct_step2
  real(kind=db) :: wl_spectre_max, wl_spectre_min
  real, dimension(:,:), allocatable :: spectre_tmp, tab_lambda_spectre, tab_spectre, tab_spectre0, tab_bb
  real(kind=db), dimension(:), allocatable :: log_spectre, log_spectre0, log_wl_spectre
  real(kind=db), dimension(n_lambda) :: log_lambda, spectre_etoiles0

  character(len=512) :: filename

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j, naxes1_ref
  real :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull

  real(kind=db) :: L_star_spectre, Cst0, L_UV, L_star_spectre0, correct_UV
  real(kind=db) :: wl_spectre_avg, wl_deviation


  ! pour obtenir un spectre normalise a 1 Rsun et 1pc
  ! Cst0 sert a remormaliser le corps noir pour l'aligner sur les spectres
  Cst0 =  2.0*hp*c_light**2 * 1e-6
  surface= pi*Rsun_to_AU**2 ! pour 1 rayon solaire
  Cst0 = Cst0 * surface / (pc_to_AU**2)

  ! pour le spectre de MCFOST
  cst_spectre_etoiles = 2.0*pi*hp*c_light**2 * (AU_to_m)**2 ! cst BB + r_etoile en AU

  if (letape_th) then
     fact_sup = exp(0.5_db/real(n_lambda,kind=db)*log(lambda_max/lambda_min))
     fact_inf = 1.0_db/fact_sup
  else
     fact_sup = 1.0001_db ;
     fact_inf = 0.9999_db ;
  endif

  ! Emission des etoiles
  ! Le corps noir est ici defini a une cst pres
  ! car on regarde en relatif par a ce qui est fait dans repartition_energie pour le disque
  ! En particulier, on a enleve un pi partout


  ! Luminosite etoile (le facteur 4Pi a ete simplifie TODO: avec quoi ???)
  !  L_etoile=r_etoile**2*sigma*T_etoile**4   ! tout en SI sauf R en AU
  L_etoile=sum((etoile(:)%r)**2*sigma*(etoile(:)%T)**4 )  ! tout en SI sauf R en AU

  do i=2, n_etoiles
     if (etoile(i)%lb_body .neqv. (etoile(1)%lb_body)) then
        write(*,*) "Error : all stars must be black bodies or"
        write(*,*) "all stars must not be black bodies : you cannot mix"
        stop
     endif
  enddo


  if (etoile(1)%lb_body) then ! les étoiles sont des corps noirs
     ! Creation d'un corps a haute resolution en F_lambda
     ! R = 1Rsun et distance = 1pc
     n_lambda_spectre = 1000

     do i=1, n_etoiles

        if (i==1) then
           allocate(tab_lambda_spectre(n_etoiles,n_lambda_spectre), &
                tab_spectre(n_etoiles,n_lambda_spectre), tab_spectre0(n_etoiles,n_lambda_spectre), &
                tab_bb(n_etoiles,n_lambda_spectre))
           tab_lambda_spectre = 0.0 ; tab_spectre = 0.0 ;  tab_spectre0 = 0.0 ;  tab_bb = 0.0
           allocate(log_spectre(n_lambda_spectre), log_spectre0(n_lambda_spectre), log_wl_spectre(n_lambda_spectre))
        endif
        tab_lambda_spectre(i,:) = 1.0_db * spanl(lambda_min, lambda_max, n_lambda_spectre)

        do l=1, n_lambda_spectre
           wl = tab_lambda_spectre(i,l) *1.e-6
           cst_wl=cst_th/(etoile(i)%T*wl)
           tab_spectre(i,l) = max(Cst0/ ( ((exp(min(cst_wl,700.)) -1.)+1.e-30) * (wl**5)), 1e-200_db) ;
        enddo ! l

     enddo !i

  else ! les étoiles ne sont pas des corps noirs
     ! On calcule 2 trucs en meme tps :
     ! - prob_E_star : proba a lambda fixe d'emission en fonction de l'etoile
     ! - spectre_etoile : proba d'emettre a lambda pour toutes les étoiles
     ! Lecture des spectres
     do i=1, n_etoiles
        ! --- Lecture du spectre stellaire
        ! --- Les flux sont normalises pour R=1 rayon solaire vu de 1 pc
        ! --- Unites : F_lambda : W.m-2 / micron
        !              lambda : micron
        !     ---> lambda x F_lambda : W.m-2

        filename=trim(star_dir)//trim(etoile(i)%spectre)

        status=0
        !  Get an unused Logical Unit Number to use to open the FITS file.
        call ftgiou(unit,status)

        readwrite=0
        call ftopen(unit,filename,readwrite,blocksize,status)
        if (status /= 0)then
           write(*,*) 'ERROR opening fits file '//trim(filename)
           write(*,*) "Exiting"
           stop
        end if

        !  determine the size of the image
        call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
        !  check that it found both NAXIS1 and NAXIS2 keywords
        if (nfound /= 2)then
           write(*,*) 'ERROR failed to read the NAXISn keywords in '//trim(filename)
           stop
        end if

        if (naxes(2) /= 3) then
           print *, "ERROR : "//trim(filename)//" does not have the right shape"
           stop
        endif

        ! check all spectra have same dimensions
        if (i==1) then
           naxes1_ref = naxes(1)
        else
           if (naxes1_ref /= naxes(1)) then
              print *, "Error : all stellar spactra must have the same shape"
              stop
           endif
        endif

        !  initialize variables.
        npixels=naxes(1)*naxes(2)
        group=1
        firstpix=1
        nullval=-999
        nbuffer=npixels

        n_lambda_spectre = naxes(1);

        ! Il faut que tous les spectres aient la meme dim sinon ca merde
        if (i==1) allocate(spectre_tmp(n_lambda_spectre,3))

        ! read_image
        call ftgpve(unit,group,firstpix,nbuffer,nullval,spectre_tmp,anynull,status)

        if (i==1) then
           allocate(tab_lambda_spectre(n_etoiles,n_lambda_spectre), &
                tab_spectre(n_etoiles,n_lambda_spectre), tab_spectre0(n_etoiles,n_lambda_spectre), &
           tab_bb(n_etoiles,n_lambda_spectre))
           tab_lambda_spectre = 0.0 ; tab_spectre = 0.0 ;  tab_spectre0 = 0.0 ;  tab_bb = 0.0

           allocate(log_spectre(n_lambda_spectre), log_spectre0(n_lambda_spectre), log_wl_spectre(n_lambda_spectre))
        endif

        tab_lambda_spectre(i,:) = spectre_tmp(:,1)
        tab_spectre(i,:) = spectre_tmp(:,2)
        tab_bb(i,:) = spectre_tmp(:,3)
        !tab_spectre(i,:) = tab_bb(i,:)

        call ftclos(unit, status)
        call ftfiou(unit, status)

        if (i==n_etoiles) deallocate(spectre_tmp)

     enddo ! n_etoiles

  endif ! bb

  tab_spectre0(:,:) = tab_spectre(:,:)


  ! Luminosite etoile integree sur le spectre
  L_star_spectre0 = 0.0
  do l = 2, n_lambda_spectre
     L_star_spectre0 = L_star_spectre0 + 0.5 * (tab_spectre(1,l) + tab_spectre(1,l-1)) &
          * (tab_lambda_spectre(1,l) - tab_lambda_spectre(1,l-1))
  enddo
  correct_step2 = (sigma*(etoile(1)%T)**4 * (Rsun_to_AU/pc_to_AU)**2) / L_star_spectre0


  ! Flux UV additionnel pour ProDiMo
  fUV_ProDiMo = etoile(1)%fUV
  p = etoile(1)%slope_UV

  wl_inf = 91.2e-9 ; ! en m
  wl_sup = 250e-9 ; ! en m

  ! wl doit etre en microns ici
  if (abs(p+1.0) > 1.0e-5) then
     cst_UV_ProDiMo =  fUV_ProDiMo * L_star_spectre0 * (p+1) / (wl_sup**(p+1) - wl_inf**(p+1)) / 1e6 !/ (1e6)**(p+1)
  else
     cst_UV_ProDiMo =  fUV_ProDiMo * L_star_spectre0 * log(wl_sup/wl_inf) / 1e6 !/ (1e6)**(p+1)
  endif


  ! On ajoute l'UV que avant le pic de Wien
  do l = 1, n_lambda_spectre
     if (tab_lambda_spectre(1,l)  < 2898./etoile(1)%T ) then
        wl = tab_lambda_spectre(1,l) * 1e-6
        UV_ProDiMo =  cst_UV_ProDiMo * wl**p
        if (UV_ProDiMo >  tab_spectre(1,l)) tab_spectre(1,l) = UV_ProDiMo
     endif
  enddo


  L_UV = 0.0
  do l = 2, n_lambda_spectre
     wl = tab_lambda_spectre(1,l) * 1e-6  ! en m
     if ((wl > wl_inf).and.(wl < wl_sup)) then
        L_UV = L_UV + 0.5 * (tab_spectre(1,l) + tab_spectre(1,l-1)) * (tab_lambda_spectre(1,l) - tab_lambda_spectre(1,l-1))
     endif
  enddo
  ! OK: L_UV/L_star_spectre0 est egal a fUV

  ! Luminosite etoile integree sur le spectre
  L_star_spectre = 0.0
  do l = 2, n_lambda_spectre
     L_star_spectre = L_star_spectre + 0.5 * (tab_spectre(1,l) + tab_spectre(1,l-1)) &
          * (tab_lambda_spectre(1,l) - tab_lambda_spectre(1,l-1))
  enddo

  !correct_UV = L_star_spectre / L_star_spectre0 ! Bug --> donne une petite difference: le calcul doit etre fait apres resampling
!  write(*,*) "Verif", real(L_star_spectre), real(sigma*(etoile(1)%T)**4 * (Rsun_to_AU/pc_to_AU)**2)


  !---------------------------------------------------------------------------
  ! On calcule 2 trucs en meme tps :
  ! - prob_E_star : proba a lambda fixe d'emission en fonction de l'etoile
  ! - spectre_etoile : proba d'emettre a lambda pour toutes les étoiles
  !---------------------------------------------------------------------------

  !---------------------------------------------------
  ! Integration du spectre dans les bandes de MCFOST
  !---------------------------------------------------
  wl_spectre_max = maxval(tab_lambda_spectre(1,:))
  wl_spectre_min = minval(tab_lambda_spectre(1,:))

  log_spectre = log(tab_spectre(1,:) + 1e-30)
  log_spectre0 = log(tab_spectre0(1,:) + 1e-30)
  log_wl_spectre = log(tab_lambda_spectre(1,:))
  log_lambda = log(tab_lambda(:))

  do lambda=1, n_lambda
     wl = tab_lambda(lambda)*1.e-6
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     ! delta_wl est la largeur du bin d'intégration
     spectre = 0.0 ;
     spectre0 = 0.0 ;
     prob_E_star(lambda,0) = 0.0

     wl_inf =  tab_lambda_inf(lambda)
     wl_sup =  tab_lambda_sup(lambda)

     do i=1,n_etoiles
        surface=4*pi*(etoile(i)%r**2)

        ! calcul de terme en binnant le spectre d'entree
        terme = 0.0 ; terme0 = 0.0 ; N = 0

        wl_spectre_avg = 0.0
        do l=1, n_lambda_spectre
           if ( (tab_lambda_spectre(i,l) > wl_inf).and.(tab_lambda_spectre(i,l) < wl_sup) ) then
              terme = terme + tab_spectre(i,l)
              terme0 = terme0 + tab_spectre0(i,l)
              wl_spectre_avg = wl_spectre_avg + tab_lambda_spectre(i,l)
              N = N + 1
           endif
        enddo ! l

        if (N>1) wl_spectre_avg = wl_spectre_avg / N
        wl_deviation = wl_spectre_avg / tab_lambda(lambda) ! Deviation between bin center and averaged wl

        ! Correction eventuelles
        if ((terme > tiny_db) .and. (N>3) .and. (abs(wl_deviation-1.0) < 3e-2)) then
           terme = terme / N * (surface / Cst0)
           terme0 = terme0 / N * (surface / Cst0)
        else ! on est en dehors du spectre fournit
           if (tab_lambda(lambda) < wl_spectre_min) then
              cst_wl=cst_th/(etoile(i)%T*wl)
              terme =  surface/((wl**5)*(exp(cst_wl)-1.0))  ! BB faute de mieux
              terme0 = terme
           else if (tab_lambda(lambda) > wl_spectre_max) then ! extrapolation loi de puissance -2 en lambda.F_lambda
              ! Le corps noir ne marche pas car le niveau est plus eleve a cause du line-blanketing
              ! Donc je fais une extrapolation en loi de puissance
              terme = (surface / Cst0) * tab_spectre(i,n_lambda_spectre) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
              terme0 = (surface / Cst0) * tab_spectre0(i,n_lambda_spectre) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
           else
              terme = (surface / Cst0) * exp(interp(log_spectre, log_wl_spectre, log_lambda(lambda)))
              terme0 = (surface / Cst0) * exp(interp(log_spectre0, log_wl_spectre, log_lambda(lambda)))
           endif
        endif ! Fin correction
        spectre = spectre + terme * delta_wl
        spectre0 = spectre0 + terme0 * delta_wl

        ! Pas de le delta_wl car il faut comparer a emission du disque
        prob_E_star(lambda,i) = prob_E_star(lambda,i-1) +  terme
     enddo ! i, etoiles
     spectre_etoiles(lambda) = spectre
     spectre_etoiles0(lambda) = spectre0

     ! Emission totale des etoiles
     ! Correction par le Teff dans le fichier de parametres
     E_stars(lambda) = prob_E_star(lambda,n_etoiles)  * correct_step2

     ! Normalisation a 1 de la proba d'emissition des etoiles
     prob_E_star(lambda,:) = prob_E_star(lambda,:)/prob_E_star(lambda,n_etoiles)
  enddo ! lambda

  correct_UV = sum(spectre_etoiles) / sum(spectre_etoiles0)


  ! Multiplication par rayon etoile et distance (en Rsun et pc)
  ! spectre_etoiles est F_lambda * dlambda
  spectre_etoiles(:) =  spectre_etoiles(:) * cst_spectre_etoiles


  if ( (lProDiMo).and.(.not.allocated(ProDiMo_star_HR)) ) then
     ! 1 seule etoile en mode ProDiMo
     ! ProDiMo_star_HR est du lambda.Flambda (idem spectre_etoiles mais avec tab_lambda au lieu de tab_delta_lambda)
     allocate(ProDiMo_star_HR(n_lambda_spectre,2))
     ProDiMo_star_HR(:,1) = tab_lambda_spectre(1,:)
     ProDiMo_star_HR(:,2) = tab_spectre(1,:) * (surface / Cst0) * cst_spectre_etoiles  * tab_lambda_spectre(1,:)
  endif


  !  TODO : L_etoile doit etre recalcule
  ! L_etoile fixe le flux dans sed1
  ! E_stars fixe le flux dans sed2
  L_etoile = sum((etoile(:)%r)**2*sigma*(etoile(:)%T)**4 ) * correct_UV  ! tout en SI sauf R en AU

  !write(*,*) "Verif", real(L_star_spectre), real(sigma*(etoile(1)%T)**4 * (Rsun_to_AU/pc_to_AU)**2)  * correct_UV

  ! Energie d'un photon / Temps
  ! Ca donne lambda Flambda a une distance = R_etoile
  ! E_photon = sigma*T_etoile**4  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)) * real(N_thet)*real(N_phi)
  ! Ca donne lambda Flambda sur le detecteur
  if (l_sym_centrale) then
     E_photon = L_etoile  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)*(distance*pc_to_AU)**2) * real(N_thet)*real(N_phi)
  else
     E_photon = L_etoile  / (real(nbre_photons_loop)*real(nbre_photons_eq_th)*(distance*pc_to_AU)**2) * real(2*N_thet)*real(N_phi)
  endif

  L_bol0 = L_etoile  / ((distance*pc_to_AU)**2) * real(N_thet)*real(N_phi) ! A comparer a L_bol1 OK

  nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)


  ! Proba cumulee
  if (n_lambda > 1) then
     ! Normalisation a 1 du spectre
     do lambda =1, n_lambda
        spectre_etoiles_cumul(lambda) = spectre_etoiles_cumul(lambda-1) + spectre_etoiles(lambda)
     enddo
     do lambda=1,n_lambda
        spectre_etoiles_cumul(lambda)=spectre_etoiles_cumul(lambda)/spectre_etoiles_cumul(n_lambda)
     enddo
  endif

  return

end subroutine repartition_energie_etoiles

!***********************************************************

subroutine emit_packet_ISM(id,ri,zj,x,y,z,u,v,w,stokes,lintersect)
! Choisit la position d'emission uniformement
! sur une sphere et la direction de vol
! suivant le cos de l'angle / normale
! C. Pinte
! 27/05/09

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: id
  integer, intent(out) :: ri, zj
  real(kind=db), intent(out) :: x, y, z, u, v, w
  real(kind=db), dimension(4), intent(out) :: stokes
  logical, intent(out) :: lintersect

  real :: aleat1, aleat2, aleat3, aleat4
  real(kind=db) :: srw02, argmt, r_etoile, cospsi, phi, l, w2

  ! Energie a 1
  stokes(:) = 0. ; stokes(1)  = 1.

  ! Position de depart aleatoire sur une sphere de rayon 1
  aleat1 = sprng(stream(id))
  aleat2 = sprng(stream(id))

  z = 2.0_db * aleat1 - 1.0_db
  srw02 = sqrt(1.0-z*z)
  argmt = pi*(2.0_db*aleat2-1.0_db)
  x = srw02 * cos(argmt)
  y = srw02 * sin(argmt)

  ! Choix direction de vol : sphere uniforme
  ! emission vers l'interieur
  aleat3 = sprng(stream(id))
  aleat4 = sprng(stream(id))

  cospsi = -sqrt(aleat3)
  phi = 2.0_db*pi*aleat4
  ! (x,y,z) definit la normale (ici, il a encore une norme=1)
  ! cospsi et phi sont definis / cette normale
  ! il faut faire une rotation
  call cdapres(cospsi, phi, x, y, z, u, v, w)

  w2=1.0_db-w*w

  ! Position de depart aleatoire sur une sphere de rayon r_etoile
  l = R_ISM * Rmax
  x = x * l
  y = y * l
  z = z * l

  call move_to_grid(x,y,z,u,v,w,ri,zj,lintersect)


  return

end subroutine emit_packet_ISM

!***********************************************************

end module stars
