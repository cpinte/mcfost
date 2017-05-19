module stars

  use parametres
  use prop_star
  use utils
  use constantes
  use ProDiMo
  use grid

  implicit none

  public :: spectre_etoiles

  public :: allocate_stellar_spectra, deallocate_stellar_spectra, em_sphere_uniforme, emit_packet_ism, &
       repartition_energie_ism, repartition_energie_etoiles, select_etoile, stars_cell_indices

  private


  contains

subroutine allocate_stellar_spectra()

  integer :: alloc_status

  allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), E_stars(n_lambda),  &
       E_disk(n_lambda), E_ISM(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error CDF_E_star'
     stop
  endif
  CDF_E_star = 0.0
  prob_E_star = 0.0
  E_stars = 0.0
  E_disk = 0.0
  E_ISM = 0.0

  allocate(spectre_etoiles_cumul(0:n_lambda),spectre_etoiles(n_lambda), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error spectre_etoile'
     stop
  endif
  spectre_etoiles_cumul = 0.0
  spectre_etoiles = 0.0

  return

end subroutine allocate_stellar_spectra

!**********************************************************************

subroutine deallocate_stellar_spectra()

  if (allocated(spectre_etoiles)) deallocate(spectre_etoiles,spectre_etoiles_cumul)
  if (allocated(CDF_E_star)) deallocate(CDF_E_star,prob_E_star,E_stars,E_disk,E_ISM)

  return

end subroutine deallocate_stellar_spectra

!**********************************************************************

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
     if (CDF_E_star(lambda,k) < aleat) then
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

subroutine em_sphere_uniforme(id, i_star,aleat1,aleat2,aleat3,aleat4, icell,x,y,z,u,v,w,w2,lintersect)
! Choisit la position d'emission uniformement
! sur la surface de l'etoile et la direction de vol
! suivant le cos de l'angle / normale
! C. Pinte
! 21/05/05

  implicit none

  integer, intent(in) :: id, i_star
  real, intent(in) :: aleat1, aleat2, aleat3, aleat4
  integer, intent(out) :: icell

  real(kind=dp), intent(out) :: x, y, z, u, v, w, w2
  logical, intent(out) :: lintersect

  real(kind=dp) :: srw02, argmt, r_etoile, cospsi, phi

  ! Position de depart aleatoire sur une sphere de rayon 1
  z = 2.0_dp * aleat1 - 1.0_dp
  srw02 = sqrt(1.0-z*z)
  argmt = pi*(2.0_dp*aleat2-1.0_dp)
  x = srw02 * cos(argmt)
  y = srw02 * sin(argmt)

  ! Choix direction de vol : sphere uniforme
  cospsi = sqrt(aleat3) ! !TODO : c'est bizarre ca. aleat3 marche avec le benchmark : mais pas la meme chose pour la surface stellaire dans RT. sqrt(aleat3) ~ OK avec TORUS en mode MC
  phi = 2.0_dp*pi*aleat4
  ! (x,y,z) definit la normale (ici, il a encore une norme=1)
  ! cospsi et phi sont definis / cette normale
  ! il faut faire une rotation
  call cdapres(cospsi, phi, x, y, z, u, v, w)

  w2=1.0_dp-w*w

  ! Position de depart aleatoire sur une sphere de rayon r_etoile
  r_etoile = etoile(i_star)%r
  x = x * r_etoile
  y = y * r_etoile
  z = z * r_etoile

  ! Ajout position de l'étoile
  x=x+etoile(i_star)%x
  y=y+etoile(i_star)%y
  z=z+etoile(i_star)%z

  if (lVoronoi) then
     icell = etoile(i_star)%icell
  else ! star can overlap several cells on a classical grid
     call indice_cellule(x,y,z, icell)
  endif

  if (etoile(i_star)%out_model) then
     call move_to_grid(id, x,y,z,u,v,w, icell,lintersect)
  else
     lintersect = .true.
  endif

  return

end subroutine em_sphere_uniforme

!**********************************************************************

!subroutine em_etoile_ponctuelle(n_star,aleat1,aleat2,ri,zj,phik,x,y,z,u,v,w,w2)
!! Emission isotrope
!! C. Pinte
!! 21/05/05
!
!  implicit none
!
!  integer, intent(in) :: n_star
!  real, intent(in) :: aleat1, aleat2
!  integer, intent(out) :: ri, zj, phik
!  real(kind=dp), intent(out) :: x, y, z, u, v, w, w2
!
!  real(kind=dp) :: srw02, argmt
!
!  ! Emission isotrope
!  w = 2.0_dp * aleat1 - 1.0_dp
!  w2 = 1.0_dp-w*w
!  srw02 = sqrt(w2)
!  argmt = pi*(2.0_dp*aleat2-1.0_dp)
!  u = srw02 * cos(argmt)
!  v = srw02 * sin(argmt)
!
!  ! Position de l'étoile
!  x=etoile(n_star)%x
!  y=etoile(n_star)%y
!  z=etoile(n_star)%z
!
!  ri=etoile(n_star)%ri
!  zj=etoile(n_star)%zj
!  phik=etoile(n_star)%phik
!
!  return
!
!end subroutine em_etoile_ponctuelle

!**********************************************************************

subroutine loi_Wien(lambda)

  implicit none

  integer, intent(out) :: lambda
  real :: T
  real :: wl

  T = maxval(etoile%T)

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
!
! renvoie :
! - spectre_etoiles
! - CDF_E_star, prob_E_star
! - E_stars
! - spectre_etoiles_cumul, spectre_etoiles
! - L_etoile

  implicit none

  real, dimension(n_lambda,n_etoiles) :: prob_E_star0
  real(kind=dp), dimension(n_lambda) :: log_lambda, spectre_etoiles0
  real(kind=dp), dimension(n_etoiles) ::  L_star0, correct_Luminosity

  real, dimension(:,:), allocatable :: spectre_tmp, tab_lambda_spectre, tab_spectre, tab_spectre0, tab_bb
  real(kind=dp), dimension(:), allocatable :: log_spectre, log_spectre0, log_wl_spectre
  character(len=512) :: filename, dir

  integer :: lambda, i, n, n_lambda_spectre, l, ios
  real(kind=dp) :: wl, cst_wl, delta_wl, surface, terme, terme0, spectre, spectre0, Cst0
  real ::  wl_inf, wl_sup, UV_ProDiMo, p, cst_UV_ProDiMo, correct_UV
  real(kind=dp) :: fact_sup, fact_inf, cst_spectre_etoiles

  real(kind=dp) :: wl_spectre_max, wl_spectre_min, wl_spectre_avg, wl_deviation

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, naxes1_ref
  real :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull

  if (n_etoiles < 1) return

  ! pour obtenir un spectre normalise a 1 Rsun et 1pc
  ! Cst0 sert a remormaliser le corps noir pour l'aligner sur les spectres
  Cst0 =  2.0*hp*c_light**2 * 1e-6
  surface= pi*Rsun_to_AU**2 ! pour 1 rayon solaire
  Cst0 = Cst0 * surface / (pc_to_AU**2)

  ! pour le spectre de MCFOST
  cst_spectre_etoiles = 2.0*pi*hp*c_light**2 * (AU_to_m)**2 ! cst BB + r_etoile en AU

  if (letape_th) then
     fact_sup = exp(0.5_dp/real(n_lambda,kind=dp)*log(lambda_max/lambda_min))
     fact_inf = 1.0_dp/fact_sup
  else
     fact_sup = 1.0001_dp ;
     fact_inf = 0.9999_dp ;
  endif

  ! Emission des etoiles
  ! Le corps noir est ici defini a une cst pres
  ! car on regarde en relatif par a ce qui est fait dans repartition_energie pour le disque
  ! En particulier, on a enleve un pi partout


  ! Luminosite etoile (le facteur 4Pi a ete simplifie avec le 1/4pi de l'estimateur de Lucy 1999)
  ! L_etoile=r_etoile**2*sigma*T_etoile**4   ! tout en SI sauf R en AU
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
        tab_lambda_spectre(i,:) = 1.0_dp * spanl(lambda_min, lambda_max, n_lambda_spectre)

        do l=1, n_lambda_spectre
           wl = tab_lambda_spectre(i,l) *1.e-6
           cst_wl=cst_th/(etoile(i)%T*wl)
           tab_spectre(i,l) = max(Cst0/ ( ((exp(min(cst_wl,700.)) -1.)+1.e-30) * (wl**5)), 1e-200_dp) ;
        enddo ! l

     enddo !i

  else ! les étoiles ne sont pas des corps noirs
     ! On calcule 2 trucs en meme tps :
     ! - CDF_E_star : proba cumulee a lambda fixe d'emission en fonction de l'etoile
     ! - spectre_etoile : proba d'emettre a lambda pour toutes les étoiles
     ! Lecture des spectres
     do i=1, n_etoiles
        ! --- Lecture du spectre stellaire
        ! --- Les flux sont normalises pour R=1 rayon solaire vu de 1 pc
        ! --- Unites : F_lambda : W.m-2 / micron
        !              lambda : micron
        !     ---> lambda x F_lambda : W.m-2

        filename=trim(etoile(i)%spectre)
        dir = in_dir(filename, star_dir,  status=ios)
        if (ios /=0) then
           write(*,*) "ERROR: star file cannot be found:",trim(filename)
           write(*,*) "Exiting"
           stop
        else
           filename = trim(dir)//trim(filename) ;
           write(*,*) "Reading "//trim(filename) ;
        endif

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
        call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
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

  ! Luminosite etoile integree sur le spectre
  L_star0(:) = 0.0
  do l = 2, n_lambda_spectre
     L_star0(:) = L_star0 + 0.5 * (tab_spectre(:,l) + tab_spectre(:,l-1)) &
             * (tab_lambda_spectre(:,l) - tab_lambda_spectre(:,l-1))
  enddo
  correct_Luminosity(:) = (sigma*(etoile(:)%T)**4 * (Rsun_to_AU/pc_to_AU)**2) / L_star0

  ! normalisation du spectre a la luminosite indique dans le fichier de para
  do i=1, n_etoiles
     tab_spectre(i,:) = tab_spectre(i,:) * correct_Luminosity(i)
  enddo

  ! Sauvegarde spectre stellaire avant d'ajouter UV
  tab_spectre0(:,:) = tab_spectre(:,:)

  ! Flux UV additionnel pour ProDiMo
  wl_inf = 91.2e-9 ; ! en m
  wl_sup = 250e-9 ; ! en m

  ! wl doit etre en microns ici
  do i=1, n_etoiles
     fUV_ProDiMo = etoile(i)%fUV
     p = etoile(i)%slope_UV

     if (abs(p+1.0) > 1.0e-5) then
        cst_UV_ProDiMo =  fUV_ProDiMo * L_star0(i) * (p+1) / (wl_sup**(p+1) - wl_inf**(p+1)) / 1e6 !/ (1e6)**(p+1)
     else
        cst_UV_ProDiMo =  fUV_ProDiMo * L_star0(i) * log(wl_sup/wl_inf) / 1e6 !/ (1e6)**(p+1)
     endif

     ! On ajoute l'UV que avant le pic de Wien
     do l = 1, n_lambda_spectre
        if (tab_lambda_spectre(i,l)  < 2898./etoile(i)%T ) then
           wl = tab_lambda_spectre(1,l) * 1e-6
           UV_ProDiMo =  cst_UV_ProDiMo * wl**p
           if (UV_ProDiMo >  tab_spectre(i,l)) tab_spectre(i,l) = UV_ProDiMo
        endif
     enddo

  enddo ! n_etoiles


  !---------------------------------------------------------------------------
  ! On calcule 2 trucs en meme tps :
  ! - CDF_E_star : proba cumulee a lambda fixe d'emission en fonction de l'etoile
  ! - spectre_etoile : proba d'emettre a lambda pour toutes les étoiles
  !---------------------------------------------------------------------------

  !---------------------------------------------------
  ! Integration du spectre dans les bandes de MCFOST
  !---------------------------------------------------
  spectre_etoiles(:) = 0.0
  spectre_etoiles0(:) = 0.0

  log_lambda = log(tab_lambda(:))
  do i=1,n_etoiles
     surface=4*pi*(etoile(i)%r**2)

     wl_spectre_max = maxval(tab_lambda_spectre(i,:))
     wl_spectre_min = minval(tab_lambda_spectre(i,:))

     log_spectre = log(tab_spectre(i,:) + 1e-30)
     log_spectre0 = log(tab_spectre0(i,:) + 1e-30)
     log_wl_spectre = log(tab_lambda_spectre(i,:))

     do lambda=1, n_lambda
        wl = tab_lambda(lambda)*1.e-6
        delta_wl=tab_delta_lambda(lambda)*1.e-6
        ! delta_wl est la largeur du bin d'intégration
        CDF_E_star(lambda,0) = 0.0

        wl_inf =  tab_lambda_inf(lambda)
        wl_sup =  tab_lambda_sup(lambda)

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
        if ((terme > tiny_dp) .and. (N>3) .and. (abs(wl_deviation-1.0) < 3e-2)) then
           terme = terme / N * (surface / Cst0)
           terme0 = terme0 / N * (surface / Cst0)
        else ! on est en dehors du spectre fournit
           if (tab_lambda(lambda) < wl_spectre_min) then
              cst_wl=cst_th/(etoile(i)%T*wl)
              if (cst_wl < 500.) then
                 terme =  surface/((wl**5)*(exp(cst_wl)-1.0))  ! BB faute de mieux
              else
                 terme = tiny_real
              endif
              terme0 = terme
           else if (tab_lambda(lambda) > wl_spectre_max) then ! extrapolation loi de puissance -2 en lambda.F_lambda
              ! Le corps noir ne marche pas car le niveau est plus eleve a cause du line-blanketing
              ! Donc je fais une extrapolation en loi de puissance
              terme = (surface / Cst0) * tab_spectre(i,n_lambda_spectre) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
              terme0 = (surface / Cst0) * tab_spectre0(i,n_lambda_spectre) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
           else
              !write(*,*) log_spectre
              terme = (surface / Cst0) * exp(interp(log_spectre, log_wl_spectre, log_lambda(lambda)))
              terme0 = (surface / Cst0) * exp(interp(log_spectre0, log_wl_spectre, log_lambda(lambda)))
           endif
        endif ! Fin correction

        ! Pas de le delta_wl car il faut comparer a emission du disque
        prob_E_star(lambda,i) = terme
        prob_E_star0(lambda,i) = terme0
      enddo ! lambda
  enddo ! etoiles


  ! On inverse les boucles pour faire les sommations sur les etoiles a lambda donne
  do lambda=1, n_lambda
     delta_wl=tab_delta_lambda(lambda)*1.e-6

     spectre = 0.0
     spectre0 = 0.0 ;
     CDF_E_star(lambda,0) = 0.0

     do i=1, n_etoiles
        terme = prob_E_star(lambda,i)
        terme0 = prob_E_star0(lambda,i)

        spectre = spectre + terme * delta_wl  ! Somme sur toutes les etoiles
        spectre0 = spectre0 + terme0 * delta_wl

        ! Pas de le delta_wl car il faut comparer a emission du disque
        CDF_E_star(lambda,i) = CDF_E_star(lambda,i-1) +  terme
     enddo ! i, etoiles

     spectre_etoiles(lambda) = spectre_etoiles(lambda) + spectre
     spectre_etoiles0(lambda) =  spectre_etoiles0(lambda) + spectre0

     ! Emission totale des etoiles
     ! Correction par le Teff dans le fichier de parametres
     E_stars(lambda) = CDF_E_star(lambda,n_etoiles)

     ! Normalisation a 1 de la proba d'emissition des etoiles
     if (CDF_E_star(lambda,n_etoiles) > 0.) then
        prob_E_star(lambda,:) = prob_E_star(lambda,:)/CDF_E_star(lambda,n_etoiles)
        CDF_E_star(lambda,:) = CDF_E_star(lambda,:)/CDF_E_star(lambda,n_etoiles)
     endif
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
     ProDiMo_star_HR(:,2) = tab_spectre(1,:) * (surface / Cst0) * cst_spectre_etoiles  * tab_lambda_spectre(1,:) * 1e-6
  endif

  !  TODO : L_etoile doit etre recalcule
  ! L_etoile fixe le flux dans sed1
  ! E_stars fixe le flux dans sed2
  L_etoile = sum((etoile(:)%r)**2*sigma*(etoile(:)%T)**4 ) * correct_UV  ! tout en SI sauf R en AU

  !write(*,*) "Verif", real(L_star_spectre), real(sigma*(etoile(1)%T)**4 * (Rsun_to_AU/pc_to_AU)**2)  * correct_UV

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

  ! Cellule d'origine dans laquelle est l'etoile
  if (.not.lVoronoi) then ! already done during tesselation for Voronoi grid
     do i=1, n_etoiles
        call indice_cellule(etoile(i)%x,etoile(i)%y,etoile(i)%z, etoile(i)%icell)
     enddo
  endif

  return

end subroutine repartition_energie_etoiles

!***********************************************************

subroutine repartition_energie_ISM(ISM_model)

  integer, intent(in) :: ISM_model ! 0 : no ISM radiation field, 1 : ProDiMo, 2 : Bate & Keto

  real :: wl, nu, ev_to_Hz, E, nu_p_MIR, Jnu, E_ev
  integer :: lambda, k

  real, dimension(5) :: wavelengths, power, W, T


  eV_to_Hz = electron_charge/hp
  nu_p_MIR = c_light/100.e-6

  if (ISM_model==0) then
     E_ISM(:) = 0.
  else if (ISM_model==1) then ! ProdiMo
     do lambda=1, n_lambda
        wl = tab_lambda(lambda) * 1e-6
        E_ISM(lambda) = (chi_ISM * 1.71 * Wdil * Blambda(wl,T_ISM_stars) + Blambda(wl,TCmb))
     enddo
  else if (ISM_model==2) then ! Bate & Keto
     ! Interstellar Radiation (ISR) field from Zucconi et al. (2001), following Black (1994)
     ! But adds Draine (1978) UV field (equation 11) which is similar to
     ! Black (1994) but doesn't seem to be included in Zucconi.
     ! Code adapted from sphNG routine
     wavelengths = (/0.4E-4,0.75E-4,1.0E-4,140.0E-4,1.06E-1/) * cm_to_m
     power = (/0.,0.,0.,1.65,0./)
     W = (/1e-14,1e-13,4e-13,2e-4,1/)
     T = (/7500.,4000.,3000.,23.3,2.728/)

     do lambda=1, n_lambda
        wl = tab_lambda(lambda) * 1e-6 ; nu = c_light/wl

        E = 0.
        ! Sum of black-bodies
        do k=1,5
           E = E + (wavelengths(k)/wl)**power(k) * W(k) * Blambda(wl,T(k))
        enddo

        ! Add mid-infrared which has a cut-off longer than 100 microns
        if (tab_lambda(lambda) < 100) then
           Jnu = 5.0E-7* (2.0*hp*nu_p_MIR**3/c_light**2) * (tab_lambda(lambda)/100.)**1.8

           E = E + c_light/wl**2 * Jnu ! Adding J_lambda
        endif

        ! Add Draine 1978 UV
        if ((nu > 5*ev_to_Hz).and.(nu < 13.6*eV_to_Hz)) then
           E_eV = nu/ev_to_Hz
           E = E + (1.658e6*E_eV - 2.152e5*E_eV**2 + 6.919e3*E_eV**3) * hp*E_eV
        endif

        E_ISM(lambda) = E
     enddo
  else
     write(*,*) "Unknown ISM model"
     write(*,*) "Exiting"
     stop
  endif

  ! Normalization for MCFOST
  E_ISM(:) = E_ISM(:) * (4.*(R_ISM*Rmax)**2) * 2.0/(hp *c_light**2) * 0.4

  return

end subroutine repartition_energie_ISM

!***********************************************************

subroutine emit_packet_ISM(id, icell,x,y,z,u,v,w,stokes,lintersect)
! Choisit la position d'emission uniformement
! sur une sphere et la direction de vol
! suivant le cos de l'angle / normale
! C. Pinte
! 27/05/09

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: id
  integer, intent(out) :: icell
  real(kind=dp), intent(out) :: x, y, z, u, v, w
  real(kind=dp), dimension(4), intent(out) :: stokes
  logical, intent(out) :: lintersect

  real :: aleat1, aleat2, aleat3, aleat4
  real(kind=dp) :: srw02, argmt, cospsi, phi, l, w2

  ! Energie a 1
  stokes(:) = 0. ; stokes(1)  = 1.

  ! Position de depart aleatoire sur une sphere de rayon 1
  aleat1 = sprng(stream(id))
  aleat2 = sprng(stream(id))

  z = 2.0_dp * aleat1 - 1.0_dp
  srw02 = sqrt(1.0-z*z)
  argmt = pi*(2.0_dp*aleat2-1.0_dp)
  x = srw02 * cos(argmt)
  y = srw02 * sin(argmt)

  ! Choix direction de vol : sphere uniforme
  ! emission vers l'interieur
  aleat3 = sprng(stream(id))
  aleat4 = sprng(stream(id))

  cospsi = -sqrt(aleat3)
  phi = 2.0_dp*pi*aleat4
  ! (x,y,z) definit la normale (ici, il a encore une norme=1)
  ! cospsi et phi sont definis / cette normale
  ! il faut faire une rotation
  call cdapres(cospsi, phi, x, y, z, u, v, w)

  w2=1.0_dp-w*w

  ! Position de depart aleatoire sur une sphere de rayon r_etoile
  l = R_ISM * Rmax
  x = x * l
  y = y * l
  z = z * l

  call move_to_grid(id, x,y,z,u,v,w, icell,lintersect)

  return

end subroutine emit_packet_ISM

!***********************************************************

subroutine stars_cell_indices()

  real(kind=dp) :: x, y, z
  integer :: i_star, icell

  do i_star=1, n_etoiles
     x = etoile(i_star)%x
     y = etoile(i_star)%y
     z = etoile(i_star)%z

     ! todo : l'etoile peut occuper plusieurs cellules
     call indice_cellule(x,y,z, icell)
     etoile(i_star)%icell = icell

     etoile(i_star)%out_model = test_exit_grid(icell, x,y,z)
  enddo

  return

end subroutine stars_cell_indices

end module stars
