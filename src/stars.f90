module stars

  use parametres
  use utils
  use constantes
  use messages
  use wavelengths
  use grid

  implicit none

  public :: spectre_etoiles, E_stars, ProDiMo_star_HR, R_ISM, E_ISM, prob_E_star

  public :: allocate_stellar_spectra, deallocate_stellar_spectra, em_sphere_uniforme, emit_packet_ism, &
       repartition_energie_ism, repartition_energie_etoiles, select_etoile, stars_cell_indices, find_spectra, &
       intersect_stars

  private

  real, dimension(:,:), allocatable :: CDF_E_star, prob_E_star
  real, dimension(:), allocatable :: E_stars !n_lambda
  real, dimension(:), allocatable :: spectre_etoiles_cumul, spectre_etoiles !(0:n_lambda)

  real, dimension(:), allocatable :: E_ISM
  real(kind=dp) :: R_ISM = 0._dp ! radius of the sphere from which the ISM radiation is emitted
  real(kind=dp), dimension(3) :: centre_ISM  ! centre of the ISM emitting sphere

  real, dimension(:,:), allocatable :: ProDiMo_star_HR

  contains

subroutine allocate_stellar_spectra(n_wl)

  integer, intent(in) :: n_wl
  integer :: alloc_status

  allocate(CDF_E_star(n_wl,0:n_etoiles), prob_E_star(n_wl,n_etoiles), E_stars(n_wl),  &
       E_ISM(n_wl), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error CDF_E_star')
  CDF_E_star = 0.0
  prob_E_star = 0.0
  E_stars = 0.0
  E_ISM = 0.0

  allocate(spectre_etoiles_cumul(0:n_wl),spectre_etoiles(n_wl), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error spectre_etoile')
  spectre_etoiles_cumul = 0.0
  spectre_etoiles = 0.0

  return

end subroutine allocate_stellar_spectra

!**********************************************************************

subroutine deallocate_stellar_spectra()

  if (allocated(spectre_etoiles)) deallocate(spectre_etoiles,spectre_etoiles_cumul)
  if (allocated(CDF_E_star)) deallocate(CDF_E_star,prob_E_star,E_stars,E_ISM)

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
  real(kind=dp), parameter :: precision = 1e-6_dp

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
  r_etoile = etoile(i_star)%r * (1._dp + precision)
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

  implicit none

  real, dimension(n_lambda,n_etoiles) :: prob_E_star0
  real(kind=dp), dimension(n_lambda) :: log_lambda, spectre_etoiles0
  real(kind=dp), dimension(n_etoiles) ::  L_star0, correct_Luminosity
  real(kind=dp), dimension(n_etoiles) :: Lacc, Tacc

  real, dimension(:,:), allocatable :: spectre_tmp, tab_lambda_spectre, tab_spectre, tab_spectre0, tab_bb
  real(kind=dp), dimension(:), allocatable :: log_spectre, log_spectre0, log_wl_spectre
  character(len=512) :: filename, dir

  integer, dimension(n_etoiles) :: n_lambda_spectre, unit
  integer :: lambda, i, n, l, ios, n_lambda_spectre_max, n_wl
  real(kind=dp) :: wl, cst_wl, delta_wl, surface, terme, terme0, spectre, spectre0, Cst0
  real ::  wl_inf, wl_sup, UV_ProDiMo, p, cst_UV_ProDiMo, correct_UV, fUV
  real(kind=dp) :: fact_sup, fact_inf, cst_spectre_etoiles

  real(kind=dp) :: wl_spectre_max, wl_spectre_min, wl_spectre_avg, wl_deviation

  integer :: status, readwrite, blocksize,nfound,group,firstpix,nbuffer,npixels, naxes1_ref
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

  do i=2, n_etoiles
     if (etoile(i)%lb_body .neqv. (etoile(1)%lb_body)) then
        call error("all stars must be black bodies or", &
             msg2="all stars must not be black bodies : you cannot mix")
     endif
  enddo
  if (etoile(1)%lb_body) etoile(:)%find_spectrum = .false.

  call find_spectra()

  if (etoile(1)%lb_body) then ! les étoiles sont des corps noirs
     ! Creation d'un corps a haute resolution en F_lambda
     ! R = 1Rsun et distance = 1pc
     n_lambda_spectre(:) = 1000

     do i=1, n_etoiles

        if (i==1) then
           allocate(tab_lambda_spectre(n_etoiles,n_lambda_spectre(1)), &
                tab_spectre(n_etoiles,n_lambda_spectre(1)), tab_spectre0(n_etoiles,n_lambda_spectre(1)), &
                tab_bb(n_etoiles,n_lambda_spectre(1)))
           tab_lambda_spectre = 0.0 ; tab_spectre = 0.0 ;  tab_spectre0 = 0.0 ;  tab_bb = 0.0
           allocate(log_spectre(n_lambda_spectre(1)), log_spectre0(n_lambda_spectre(1)), log_wl_spectre(n_lambda_spectre(1)))
        endif
        tab_lambda_spectre(i,:) = 1.0_dp * spanl(lambda_min, lambda_max, n_lambda_spectre(1))

        do l=1, n_lambda_spectre(1)
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
     n_lambda_spectre_max = 0
     do i=1, n_etoiles
        ! --- Lecture du spectre stellaire
        ! --- Les flux sont normalises pour R=1 rayon solaire vu de 1 pc
        ! --- Unites : F_lambda : W.m-2 / micron
        !              lambda : micron
        !     ---> lambda x F_lambda : W.m-2

        filename=trim(etoile(i)%spectre)
        dir = in_dir(filename, star_dir,  status=ios)
        if (ios /=0) then
           call error("star file cannot be found:",trim(filename))
        else
           filename = trim(dir)//trim(filename) ;
           write(*,*) "Reading "//trim(filename) ;
        endif

        status=0
        !  Get an unused Logical Unit Number to use to open the FITS file.
        call ftgiou(unit(i),status)

        readwrite=0
        call ftopen(unit(i),filename,readwrite,blocksize,status)
        if (status /= 0) call error("cannot open fits file "//trim(filename))

        !  determine the size of the image
        call ftgknj(unit(i),'NAXIS',1,10,naxes,nfound,status)
        !  check that it found both NAXIS1 and NAXIS2 keywords
        if (nfound /= 2) call error("failed to read the NAXISn keywords in "//trim(filename))

        if (naxes(2) /= 3) call error(trim(filename)//" does not have the right shape")

        ! We first read the length of the spectrum
        n_lambda_spectre(i) = naxes(1)
        if (naxes(1) > n_lambda_spectre_max) n_lambda_spectre_max = naxes(1)
     enddo

     ! We allocate the array to the maximum length
     allocate(tab_lambda_spectre(n_etoiles,n_lambda_spectre_max), &
          tab_spectre(n_etoiles,n_lambda_spectre_max), tab_spectre0(n_etoiles,n_lambda_spectre_max), &
          tab_bb(n_etoiles,n_lambda_spectre_max))
     tab_lambda_spectre = 0.0 ; tab_spectre = 0.0 ;  tab_spectre0 = 0.0 ;  tab_bb = 0.0

     allocate(log_spectre(n_lambda_spectre_max), log_spectre0(n_lambda_spectre_max), log_wl_spectre(n_lambda_spectre_max))


     do i=1, n_etoiles
        n_wl = n_lambda_spectre(i)
        ! We read again the dimensions
        call ftgknj(unit(i),'NAXIS',1,10,naxes,nfound,status)

        !  initialize variables.
        npixels=n_wl * 3
        group=1
        firstpix=1
        nullval=-999
        nbuffer=npixels

        ! read_image
        allocate(spectre_tmp(n_wl,3))
        call ftgpve(unit(i),group,firstpix,nbuffer,nullval,spectre_tmp,anynull,status)

        tab_lambda_spectre(i,1:n_wl) = spectre_tmp(:,1)
        tab_spectre(i,1:n_wl) = spectre_tmp(:,2)
        tab_bb(i,1:n_wl) = spectre_tmp(:,3)

        call ftclos(unit(i), status)
        call ftfiou(unit(i), status)
        deallocate(spectre_tmp)
     enddo ! n_etoiles

  endif ! bb

  ! Luminosite etoile integree sur le spectre
  L_star0(:) = 0.0
  do i=1, n_etoiles
     do l = 2, n_lambda_spectre(i)
        L_star0(i) = L_star0(i) + 0.5 * (tab_spectre(i,l) + tab_spectre(i,l-1)) &
             * (tab_lambda_spectre(i,l) - tab_lambda_spectre(i,l-1))
     enddo
  enddo
  correct_Luminosity(:) = (sigma*(etoile(:)%T)**4 * (Rsun_to_AU/pc_to_AU)**2) / L_star0(:)

  ! normalisation du spectre a la luminosite indique dans le fichier de para
  do i=1, n_etoiles
     tab_spectre(i,:) = tab_spectre(i,:) * correct_Luminosity(i)
  enddo

  ! Sauvegarde spectre stellaire avant d'ajouter UV
  tab_spectre0(:,:) = tab_spectre(:,:)


  !--------------------------------
  ! Flux UV additionnel pour ProDiMo
  !--------------------------------
  wl_inf = 91.2e-9 ; ! en m
  wl_sup = 250e-9 ; ! en m

  ! wl doit etre en microns ici
  do i=1, n_etoiles
     fUV = etoile(i)%fUV
     if (fUV > tiny_real) then
        p = etoile(i)%slope_UV

        if (abs(p+1.0) > 1.0e-5) then
           cst_UV_ProDiMo =  fUV * L_star0(i) * (p+1) / (wl_sup**(p+1) - wl_inf**(p+1)) / 1e6 !/ (1e6)**(p+1)
        else
           cst_UV_ProDiMo =  fUV * L_star0(i) * log(wl_sup/wl_inf) / 1e6 !/ (1e6)**(p+1)
        endif

        ! On ajoute l'UV que avant le pic de Wien
        do l = 1, n_lambda_spectre(i)
           if (tab_lambda_spectre(i,l)  < 2898./etoile(i)%T ) then
              wl = tab_lambda_spectre(1,l) * 1e-6
              UV_ProDiMo =  cst_UV_ProDiMo * wl**p
              if (UV_ProDiMo >  tab_spectre(i,l)) tab_spectre(i,l) = UV_ProDiMo
           endif
        enddo
     endif
  enddo ! n_etoiles

  !--------------------------------
  ! Calculate accretion spectrum
  !--------------------------------
  if (.not.lturn_off_Lacc) then
     if (maxval(etoile(:)%Mdot) > tiny_real) then
        ! Luminosity from the accretion [au^2 W / m^2]
        Lacc(:) = Ggrav/AU3_to_m3 &                         ! convert G to AU^3 / s^s / kg
             * etoile(:)%M*Msun_to_kg &                     ! M in kg
             * etoile(:)%Mdot*Msun_to_kg/year_to_s &        ! Mdot in kg / s
             / etoile(:)%r                                  ! R in AU
        ! Converting Lacc to Tacc
        Tacc(:) = (Lacc(:)/(quatre_pi * sigma * etoile(:)%r**2))**0.25

        write(*,*) "Accretion onto stars: "
        write(*,*) "Mdot=", etoile(:)%Mdot, "Msun/yr"
        write(*,*) "Tacc=", real(Tacc(:)), "K"

        ! We add a black-body to the stellar spectrum
        do i=1, n_etoiles
           if (Tacc(i) > tiny_real) then
              do l=1, n_lambda_spectre(i)
                 wl = tab_lambda_spectre(i,l) * 1.e-6
                 cst_wl=cst_th/(Tacc(i)*wl)
                 tab_spectre(i,l) = tab_spectre(i,l) +  max(Cst0/ ( ((exp(min(cst_wl,700.)) -1.)+1.e-30) * (wl**5)), 1e-200_dp) ;
              enddo ! l
           endif
        enddo
     endif
  else
     write(*,*) "Turning off accretion luminosity"
  endif

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
     n_wl = n_lambda_spectre(i)
     surface=4*pi*(etoile(i)%r**2)

     wl_spectre_max = maxval(tab_lambda_spectre(i,1:n_wl))
     wl_spectre_min = minval(tab_lambda_spectre(i,1:n_wl))

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
        do l=1, n_wl
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
        if ((terme > tiny_dp) .and. (N>3) .and. (abs(wl_deviation-1.0) < 0.1)) then
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
            terme = (surface / Cst0) * tab_spectre(i,n_wl) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
            terme0 = (surface / Cst0) * tab_spectre0(i,n_wl) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
          else
             !write(*,*) log_spectre
             terme = (surface / Cst0) * exp(interp(log_spectre(1:n_wl), log_wl_spectre(1:n_wl), log_lambda(lambda)))
             terme0 = (surface / Cst0) * exp(interp(log_spectre0(1:n_wl), log_wl_spectre(1:n_wl), log_lambda(lambda)))
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
     allocate(ProDiMo_star_HR(n_lambda_spectre_max,2))
     ProDiMo_star_HR(:,1) = tab_lambda_spectre(1,:)
     ProDiMo_star_HR(:,2) = tab_spectre(1,:) * (surface / Cst0) * cst_spectre_etoiles  * tab_lambda_spectre(1,:) * 1e-6
  endif

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

  ! Defining the ISM sphere
  if (lcylindrical) then
     R_ISM = 1.000001_dp * (sqrt(Rmax**2 + zmax(n_rad)**2))
     centre_ISM(:) = 0._dp
  else if (lspherical) then
     R_ISM = 1.000001_dp * Rmax
     centre_ISM(:) = 0._dp
  else if (lVoronoi) then
     ! Defining the ISM sphere
     R_ISM = 1.000001_dp * 0.5_dp * Rmax
     !centre_ISM(:) = 0.5_dp * (/limits(2)+limits(1), limits(4)+limits(3), limits(6)+limits(5)/)
     centre_ISM(:) = 0._dp
  endif

  eV_to_Hz = electron_charge/hp
  nu_p_MIR = c_light/100.e-6

  if (R_ISM < tiny_dp) call error("the ISM emitting sphere is not defined")

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
     W = (/1e-14,1e-13,4e-13,2e-4,1./)
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
     call error("Unknown ISM model")
  endif

  ! Normalization for MCFOST
  E_ISM(:) = E_ISM(:) * (4.*R_ISM**2) * 2.0/(hp *c_light**2) * 0.4

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
  l = R_ISM
  x = centre_ISM(1) + x * l
  y = centre_ISM(2) + y * l
  z = centre_ISM(3) + z * l

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

     ! todo : l'etoile peut occuper plusieurs cellules : Non,
     call indice_cellule(x,y,z, icell)
     etoile(i_star)%icell = icell

     etoile(i_star)%out_model = test_exit_grid(icell, x,y,z)
  enddo

  return

end subroutine stars_cell_indices

!***********************************************************

subroutine intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! lintersect_star is true is the ray/packet will intersect a star
  ! i_star returns the id of the star intersected
  ! icell_star returns the id of the cell where the star is

  ! This routine implies that a star is in a unique cell

  real(kind=dp), intent(in) :: x,y,z, u,v,w
  logical, intent(out) :: lintersect_stars
  integer, intent(out) :: i_star, icell_star

  real(kind=dp), dimension(3) :: r, k, delta_r
  real(kind=dp) :: b,c, delta, rac, s1, s2, d_to_star
  integer :: i


  r(1) = x ; r(2) = y ; r(3) = z
  k(1) = u ; k(2) = v ; k(3) = w

  d_to_star = huge(1.0_dp)

  i_star = 0
  star_loop : do i = 1, n_etoiles
     delta_r(:)  = r(:) - (/etoile(i)%x, etoile(i)%y, etoile(i)%z/)
     b = dot_product(delta_r,k)
     c = dot_product(delta_r,delta_r) - (etoile(i)%r)**2
     delta = b*b - c

     if (delta >= 0.) then ! the packet will encounter (or has encoutered) the star
        rac = sqrt(delta)
        s1 = -b - rac

        if (s1 < 0) then ! we already entered the star
           ! We can probably skip that test, s1 must be positive as we must be outside the star
           s2 = -b + rac
           if (s2 > 0) then ! for s2 < 0: we already exited the star
              ! We are still in the sphere and will exit it
              ! This means that we had a round-off error somewhere
              d_to_star = 0.0_dp
              i_star = i
           endif
        else ! We will enter in the star
           if (s1 < d_to_star) then
              d_to_star = s1
              i_star = i
           endif
        endif ! s1 < 0

     endif ! delta < 0

  enddo star_loop

  lintersect_stars = (i_star > 0)
  if (lintersect_stars) then
     icell_star = etoile(i_star)%icell
  else
     icell_star = 0
  end if

  return

end subroutine intersect_stars

!***********************************************************

subroutine find_spectra()
  ! Find an appropriate spectrum for all star based on the Teff, mass and radius (i.e. log(g))

  real :: Teff, r, M, logg, min_logg, max_logg
  integer :: delta_T, i_star

  real, parameter :: delta_logg = 0.5

  character(len=32) :: sTeff, slogg, type

  write(*,*) "Trying to find appropriate stellar spectra ..."

  do i_star = 1, n_etoiles
     if (etoile(i_star)%find_spectrum) then
        Teff = etoile(i_star)%T

        r = etoile(i_star)%r / Rsun_to_AU
        M = etoile(i_star)%M

        if (M < tiny_real) then
           call warning("Stellar mass is not set, forcing log(g) = 3.5")
           logg = 3.5 ! stellar mass is not defined in the parameter file, we fix logg
        else
           logg = logg_Sun + log10(M/r**2)
        endif

        if (Teff < 100) then
           call warning("Teff below 100K needs to be implemented")
           Teff = 100
           type = "cond"
           min_logg = 2.5
           max_logg = 6
           delta_T = 100
        else if (Teff < 1500) then
           type = "cond"
           min_logg = 2.5
           max_logg = 6
           delta_T = 100
        else if (Teff < 2700) then
           type = "dusty"
           if ((Teff < 2250).and.(Teff > 2050)) then
              min_logg = 4.0 ! Some models appear to be missing
           else
              min_logg = 3.5
           endif
           max_logg = 6
           delta_T = 100
        else if (Teff < 10000) then
           type = "NextGen"
           min_logg = 3.5
           max_logg = 5.5
           if (Teff <= 4000) then
              delta_T = 100
           else
              delta_T = 200
           endif
        else if (Teff <= 35000) then
           type="Kurucz"

           if (Teff < 11000) then
              min_logg = 2.0
           else if (Teff < 19000) then
              min_logg = 2.5
           else if (Teff < 27000) then
              min_logg = 3.0
           else if (Teff < 32000) then
              min_logg = 3.5
           else
              min_logg = 4.0
           endif
           max_logg = 5.0

           if (Teff <= 12000) then
              delta_T = 500
           else
              delta_T = 1000
           endif

        else
           call error("Teff above 35000K needs to be implemented")
        endif

        ! Rounding off at the nearest point in the grid of stellar atmospheres
        Teff = nint(Teff/delta_T) * delta_T
        logg = nint(logg/delta_logg) * delta_logg
        logg = min(max(logg,min_logg), max_logg)

        if (Teff < 1000) then
           write(sTeff, "(I3)") int(Teff)
        else if (Teff < 10000) then
           write(sTeff, "(I4)") int(Teff)
        else
           write(sTeff, "(I5)") int(Teff)
        endif
        write(slogg, "(F3.1)") logg

        if (type=="Kurucz") then
           etoile(i_star)%spectre = "Kurucz"//trim(sTeff)//"-"//trim(slogg)//".fits.gz"
        else
           etoile(i_star)%spectre = "lte"//trim(sTeff)//"-"//trim(slogg)//"."//trim(type)//".fits.gz"
        endif

        write(*,*) "Star #", i_star, " --> ", trim(etoile(i_star)%spectre)
     else ! We do not update the spectrum
        if (etoile(i_star)%lb_body) then
           write(*,*) "Star #", i_star, " --> BB at T=", etoile(i_star)%T, "K"
        else
           write(*,*) "Star #", i_star, " --> ", trim(etoile(i_star)%spectre), " (forced)"
        endif
     endif
  enddo

  write(*,*) "Done"

  return

end subroutine find_spectra

end module stars
