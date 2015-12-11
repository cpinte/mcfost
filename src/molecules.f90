module molecules

  use parametres
  use em_th
  use molecular_emission
  use opacity
  use prop_star
  use grid
  use dust
  use mem

  use ProDiMo, only : m2p ! to get the continuum radiation field

  implicit none

  contains

real function tau_collision(temperature, iPart, iTrans)
  ! Interpolation lineaire des taux de collisions avec dichotomie
  ! C. Pinte
  ! 25/06/07

  implicit none

  real, intent(in) :: temperature
  integer, intent(in) :: iPart, iTrans

  real :: frac
  integer :: k, kmin, kmax

  if (nCollTemps(iPart) == 1) then ! Pas d'interpolation
     tau_collision = collRates(iPart, iTrans, 1)
  else if (temperature < collTemps(iPart,1)) then ! extrapolation cubique
     tau_collision = collRates(iPart, iTrans, 1) * sqrt(temperature/collTemps(iPart,1))
  else ! interpolation lineaire
     ! dichotomie
     kmin=1
     kmax=nCollTemps(iPart)

     k=(kmax+kmin)/2

     do while ((kmax - kmin) > 1)
        if (collTemps(iPart,k) < temperature) then
           kmin = k
        else
           kmax = k
        endif
        k = (kmin + kmax)/2
     enddo   ! while
     k=kmax


     frac = (temperature - collTemps(iPart,k-1)) / (collTemps(iPart,k) - collTemps(iPart, k-1))
     tau_collision = collRates(iPart, iTrans, k-1) + frac * ( collRates(iPart, iTrans, k) -  collRates(iPart, iTrans, k-1))
  endif


  return

end function tau_collision

!***********************************************************

subroutine init_GG_Tau_mol()

  implicit none

  integer :: icell

  ldust_mol = .true.

  do icell=1,n_cells
     Temperature(icell) = 30.0 * (r_grid(icell)/100.)**(-0.5)
     Tcin(icell) = 30.0 * (r_grid(icell)/100.)**(-0.5)
  enddo

  icell = cell_map(1,1,1)
  write(*,*) "Density @ 100 AU", real(densite_gaz(icell) / 100.**3 * (sqrt(r_grid(icell)**2 + z_grid(icell)) / 100.)**2.75)

  return

end subroutine init_GG_Tau_mol

!***********************************************************

subroutine init_HH_30_mol()

  implicit none

  integer :: icell

  ldust_mol = .true.

  do icell=1,n_cells
     Temperature(icell) = 12.0 * (r_grid(icell)/100.)**(-0.55)
     vfield(icell) = 2.0 * (r_grid(icell)/100.)**(-0.55)
  enddo

  v_turb = 230.

  return

end subroutine init_HH_30_mol

!***********************************************************

subroutine init_benchmark_vanZadelhoff1()

  implicit none

  ldust_mol = .false.

  v_turb = 150._db !m.s-1
  Temperature = 20.
  Tcin = 20._db
  vfield(:) = 0.0

  masse_mol = 1.0

  linfall = .true.
  lkeplerian = .false.

  tab_abundance = abundance

  write(*,*) "Density", real(densite_gaz(cell_map(1,1,1)) * &
       (r_grid(cell_map(1,1,1))**2 + z_grid(cell_map(1,1,1))**2)/rmin**2)

  return

end subroutine init_benchmark_vanZadelhoff1

!***********************************************************

subroutine init_benchmark_vanzadelhoff2()
  ! Lecture de la structure du modele de van Zadelhoff 2a/b
  ! et ajustement sur la grille de mcfost
  ! C. Pinte
  ! 13/07/07

  implicit none

  integer, parameter :: n_lines = 50

  integer :: i, j, ri, l, icell, zj, k
  real, dimension(n_lines) :: tmp_r, tmp_nH2, tmp_T, tmp_v, tmp_vturb, log_tmp_r, log_tmp_nH2

  real :: junk, rayon, log_rayon, frac

  ldust_mol = .false.

  ! Lecture du fichier modele
  open(unit=1,file="model_1.d",status="old")
  do i=1,7
     read(1,*)
  enddo
  do i=1,n_lines
     j = n_lines - i + 1
     read(1,*) tmp_r(j), tmp_nH2(j), junk, tmp_T(j), tmp_v(j), tmp_vturb(j)
  enddo
  close(unit=1)

  ! Conversion en AU
  tmp_r(:) = tmp_r(:) * cm_to_AU

  ! Interpolation sur la grille de mcfost
  ! Distance en log
  ! densite en log
  ! T, V et vturb en lineaire cf Fig 2 van Zadelhoff 2002
  log_tmp_r(:) = log(tmp_r(:))
  log_tmp_nH2(:) = log(tmp_nH2(:))

  do ri=1, n_rad
     ! Recherche rayon dans def bench
     icell = cell_map(ri,1,1)
     rayon = sqrt(r_grid(icell)**2+z_grid(icell)**2)
     log_rayon = log(rayon)

     l=2
     ! rayon est entre tmp_r(l-1) et tmp_r(l)
     search : do i=l, n_lines
        if (tmp_r(i) >= rayon) then
           l = i
           exit search
        endif
     enddo search
     if (l > n_lines) l = n_lines

     frac = (log_rayon - log_tmp_r(l-1)) / (log_tmp_r(l) - log_tmp_r(l-1))

     do zj=1, nz
        k=1
        icell = cell_map(ri,zj,k)
        densite_gaz(icell) = exp( log_tmp_nH2(l-1) + frac * (log_tmp_nH2(l) - log_tmp_nH2(l-1)) )
        Temperature(icell) =  tmp_T(l-1) + frac * (tmp_T(l) - tmp_T(l-1))
        Tcin(icell) = tmp_T(l-1) + frac * (tmp_T(l) - tmp_T(l-1))
        vfield(icell) = tmp_v(l-1) + frac * (tmp_v(l) - tmp_v(l-1))
        v_turb(icell) = tmp_vturb(l-1) + frac * (tmp_vturb(l) - tmp_vturb(l-1))
     enddo
  enddo

  ! Conversion vitesses en m.s-1
  vfield = vfield * 1.0e3
  v_turb = v_turb * 1.0e3

  ! Conversion part.m-3
  densite_gaz(:) = densite_gaz(:) / (cm_to_m)**3

  linfall = .true.
  lkeplerian = .false.

  tab_abundance = abundance

  return

end subroutine init_benchmark_vanzadelhoff2

!***********************************************************

subroutine init_benchmark_water1()
  ! Initialisation de la structure du modele d'eau 1
  ! C. Pinte
  ! 15/10/07

  implicit none

  ldust_mol = .false.

  densite_gaz = 1.e4 / (cm_to_m)**3 ! part.m-3
  Tcin = 40.
  vfield = 0.0
  v_turb = 0.0

  linfall = .true.
  lkeplerian = .false.

  tab_abundance = abundance


  ! Pas de Cmb
  tab_Cmb_mol = 0.0

  return

end subroutine init_benchmark_water1

!***********************************************************

subroutine init_benchmark_water2()
  ! Initialisation de la structure du modele d'eau 2
  ! C. Pinte
  ! 15/10/07

  implicit none

  integer :: icell

  ldust_mol = .false.

  densite_gaz = 1.e4 / (cm_to_m)**3 ! part.m-3
  Tcin = 40.
  v_turb = 0.0

  do icell=1, n_cells
     vfield(icell) = 1e5 * sqrt(r_grid(icell)**2 + z_grid(icell)**2) * AU_to_pc
  enddo

  linfall = .true.
  lkeplerian = .false.

  tab_abundance = abundance

  ! Pas de Cmb
  tab_Cmb_mol = 0.0

  return

end subroutine init_benchmark_water2

!***********************************************************

subroutine init_benchmark_water3()
  ! Initialisation de la structure du modele d'eau 2
  ! C. Pinte
  ! 15/10/07

  implicit none

  integer, parameter :: n_lines = 100

  integer :: i, j, ri, l, zj, k, icell
  real, dimension(n_lines) :: tmp_r, tmp_nH2, tmp_Tkin, tmp_T, tmp_v, tmp_vturb
  real, dimension(n_lines) :: log_tmp_r, log_tmp_nH2, log_tmp_T, log_tmp_Tkin, log_tmp_v

  real :: rayon, log_rayon, frac

  ldust_mol = .true.

  ! Lecture du fichier modele
  open(unit=1,file="mc_100.d",status="old")
  read(1,*)
  !  radius [cm]  n(H2) [cm^-3] Tkin [K]    Tdust [K]   Vrad [km/s] FWHM [km/s]
  do i=1,n_lines
     j = n_lines - i + 1
     read(1,*) tmp_r(j), tmp_nH2(j), tmp_Tkin(j), tmp_T(j), tmp_v(j), tmp_vturb(j)
  enddo
  close(unit=1)

  ! Conversion en AU
  tmp_r(:) = tmp_r(:) * cm_to_AU

  ! Interpolation sur la grille de mcfost
  ! Distance en log
  ! densite en log
  ! Tkin, T, V et vturb
  log_tmp_r(:) = log(tmp_r(:))
  log_tmp_nH2(:) = log(tmp_nH2(:))
  log_tmp_T(:) = log(tmp_T(:))
  log_tmp_Tkin(:) = log(tmp_Tkin(:))
  log_tmp_v(:) = log(tmp_v(:) + 1.0e-30)


  do ri=1, n_rad
     ! Recherche rayon dans def bench
     icell = cell_map(ri,1,1)
     rayon = sqrt(r_grid(icell)**2+z_grid(icell)**2)
     log_rayon = log(rayon)


     if (rayon < 2.0) then
        do zj=1, nz
           k=1
           icell = cell_map(ri,zj,k)
           densite_gaz(icell) = tmp_nH2(1)
           Temperature(icell) = tmp_T(1)
           Tcin(icell) =  tmp_Tkin(1)
        enddo

     else
        l=2
        ! rayon est entre tmp_r(l-1) et tmp_r(l)
        search : do i=l, n_lines
           if (tmp_r(i) >= rayon) then
              l = i
              exit search
           endif
        enddo search
        if (l > n_lines) l = n_lines

        frac = (log_rayon - log_tmp_r(l-1)) / (log_tmp_r(l) - log_tmp_r(l-1))
        do zj=1, nz
           k=1
           icell = cell_map(ri,zj,k)
           densite_gaz(icell) = exp( log_tmp_nH2(l-1) + frac * (log_tmp_nH2(l) - log_tmp_nH2(l-1)) )
           Temperature(icell) = exp( log_tmp_T(l-1) + frac * (log_tmp_T(l) - log_tmp_T(l-1)) )
           Tcin(icell) = exp( log_tmp_Tkin(l-1) + frac * (log_tmp_Tkin(l) - log_tmp_Tkin(l-1)) )

           if (rayon < 5.95) then
              vfield(icell) = 0.0
              v_turb(icell) = 3.
           else
              vfield(icell) =  exp( log_tmp_v(l-1) + frac * (log_tmp_v(l) - log_tmp_v(l-1)) )
              v_turb(icell) = 1.
           endif
        enddo
     endif

  enddo

  ! Conversion FWHM ---> vitesse
  v_turb = v_turb / (2.*sqrt(log(2.)))

  ! Conversion vitesses en m.s-1
  vfield = vfield * 1.0e3
  v_turb = v_turb * 1.0e3

  ! Conversion part.m-3
  densite_gaz(:) = densite_gaz(:) / (cm_to_m)**3

  linfall = .true.
  lkeplerian = .false.

  tab_abundance = abundance

  return

end subroutine init_benchmark_water3

!***********************************************************

subroutine init_molecular_disk(imol)
  ! definie les tableaux vfield, v_turb et tab_abundance
  ! et lcompute_molRT

  implicit none

  integer, intent(in) :: imol
  integer :: icell

  ldust_mol  = .true.
  lkeplerian = .true.
  linfall = .false.

  ! Temperature gaz = poussiere
  if (lcorrect_Tgas) then
     write(*,*) "Correcting Tgas by", correct_Tgas
     Tcin(:) = Temperature(:)  * correct_Tgas
  else
     Tcin(:) = Temperature(:)
  endif

  ! En m.s-1
  do icell=1, n_cells
     vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
  enddo

  v_turb = vitesse_turb

  ! Abondance
  if (mol(imol)%lcst_abundance) then
     write(*,*) "Setting constant abundance"
     tab_abundance = mol(imol)%abundance
  else
     call read_abundance(imol)
  endif

  do icell=1, n_cells
     lcompute_molRT(icell) = (tab_abundance(icell) > tiny_real) .and. &
          (densite_gaz(icell) > tiny_real) .and. (Tcin(icell) > 1.)
  enddo

  return

end subroutine init_molecular_disk

!***********************************************************

subroutine init_Doppler_profiles(imol)
  ! Creation des constantes pour les profiles Doppler dans chaque cellule
  ! a vitesse systematique nulle + tableau de vitesses
  ! C. Pinte
  ! 18/03/07

  implicit none

  integer, intent(in) :: imol

  real(kind=db) :: sigma2, sigma2_m1, vmax
  integer :: icell, iv, n_speed

  n_speed = mol(imol)%n_speed_rt

  do icell=1, n_cells
     ! Utilisation de la temperature LTE de la poussiere comme temperature cinetique
     ! WARNING : c'est pas un sigma mais un delta, cf Cours de Boisse p47
     ! Consistent avec benchmark
     sigma2 =  2.0_db * (kb*Tcin(icell) / (masse_mol* g_to_kg)) + v_turb(icell)**2
     v_line(icell) = sqrt(sigma2)

     !  write(*,*) "FWHM", sqrt(sigma2 * log(2.)) * 2.  ! Teste OK bench water 1
     sigma2_m1 = 1.0_db / sigma2
     sigma2_phiProf_m1(icell) = sigma2_m1
     ! phi(nu) et non pas phi(v) donc facteur c_light et il manque 1/f0
     ! ATTENTION : il ne faut pas oublier de diviser par la freq apres
     norme_phiProf_m1(icell) = c_light / sqrt(pi * sigma2)

     ! Echantillonage du profil de vitesse dans la cellule
     ! 2.15 correspond a l'enfroit ou le profil de la raie faut 1/100 de
     ! sa valeur au centre : exp(-2.15^2) = 0.01
     vmax = sqrt(sigma2)
     tab_dnu_o_freq(icell) = largeur_profile * vmax / (real(n_speed))
     do iv=-n_speed, n_speed
        tab_deltaV(iv,icell) = largeur_profile * real(iv,kind=db)/real(n_speed,kind=db) * vmax
     enddo ! iv
     deltaVmax(icell) = largeur_profile * vmax !* 2.0_db  ! facteur 2 pour tirage aleatoire
  enddo !icell

  return

end subroutine init_Doppler_profiles

!***********************************************************

function phiProf(icell,ispeed,tab_speed)
  ! renvoie le profil de raie local a une cellule (ri, zj)
  ! sur le tableau de vitesse tab_speed
  ! Il faut diviser par la frequence de la transition pour que la normalisation soit OK !!!
  ! C. Pinte
  ! 13/07/07

  implicit none

  integer, dimension(2), intent(in) :: ispeed
  integer, intent(in) :: icell
  real(kind=db), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  real(kind=db), dimension(ispeed(1):ispeed(2)) :: phiProf

  real(kind=db) :: norme_m1, sigma2_m1

  norme_m1 = norme_phiProf_m1(icell) ! ATTENTION : il manque la frequence ici !!!!
  sigma2_m1 = sigma2_phiProf_m1(icell)

  phiProf(:) =  norme_m1 * exp(- sigma2_m1 * tab_speed(:)**2)

  return

end function phiProf

!***********************************************************

! --- function Cmb(ispeed,tab_speed)
! --- ! Loi de Planck aux differentes frequences de la
! --- ! molecule etudiee
! --- ! Bnu en SI : W.m-2.Hz-1.sr-1
! --- ! C. Pinte
! --- ! 17/07/7
! ---
! ---   implicit none
! ---
! ---   integer, dimension(2), intent(in) :: ispeed
! ---   real(kind=db), dimension(ispeed(1):ispeed(2)),intent(in) :: tab_speed
! ---
! ---   real(kind=db), dimension(ispeed(1):ispeed(2),nTrans) :: Cmb
! ---
! ---
! ---   integer :: iTrans, iv, iiTrans
! ---   real(kind=db) :: cst, cst_nu, nu
! ---
! ---   cst = 2.0_db*hp/c_light**2
! ---
! ---   ! On ne calcule le Cmb qu'aux frequences specifiees
! ---   do iTrans=1,nTrans
! ---      iiTrans = indice_Trans(iTrans)
! ---      do iv=ispeed(1),ispeed(2)
! ---         nu = Transfreq(iiTrans) * (1.0_db + tab_speed(iv)/c_light)
! ---
! ---         cst_nu = (hp * nu) / (kb * T_min)
! ---         if (cst_nu > 100._db) then
! ---            Cmb(iv,iTrans) = 0.0_db
! ---         else
! ---            Cmb(iv,iTrans) = cst * nu**3 / (exp(cst_nu)-1.0_db)
! ---         endif
! ---      enddo ! iv
! ---   enddo ! iTrans
! ---
! ---   return
! ---
! --- end function Cmb

!***********************************************************

subroutine init_tab_Cmb_mol()
  ! Loi de Planck aux differentes frequences de la
  ! molecule etudiee
  ! Bnu en SI : W.m-2.Hz-1.sr-1
  ! Sans prendre en compte le champ de vitesse << c
  ! C. Pinte
  ! 17/07/7

  implicit none

  integer :: iTrans
  real(kind=db) :: cst, cst_nu, nu

  cst = 2.0_db*hp/c_light**2

  do iTrans=1,nTrans
     nu = Transfreq(iTrans)

     cst_nu = (hp * nu) / (kb * T_Cmb)
     if (cst_nu > 100._db) then
        tab_Cmb_mol(iTrans) = 0.0_db
     else
        tab_Cmb_mol(iTrans) = cst * nu**3 / (exp(cst_nu)-1.0_db)
     endif
  enddo

  return

end subroutine init_tab_Cmb_mol

!***********************************************************

subroutine opacite_mol(imol)
  ! Calcule la fonction source dans chaque cellule
  ! etant donne les populations des niveaux
  ! C. Pinte
  ! 20/03/07

  implicit none

  integer, intent(in) :: imol
  integer :: i,j, k

  do i=1,n_rad
     bz : do j=j_start,nz
        if (j==0) cycle bz
        do k=1, n_az
           call opacite_mol_loc(i,j,k,imol)
        enddo ! k
     enddo bz ! j
  enddo ! i

  return

end subroutine opacite_mol

!***********************************************************

subroutine opacite_mol_loc(ri,zj,phik,imol)
  ! Calcule la fonction source dans chaque cellule
  ! etant donne les populations des niveaux
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: ri, zj, phik, imol

  integer :: iTrans, iiTrans
  real(kind=db) :: nu, nl, kap, eps

  logical, save :: lmaser = .false.

  character(len=128) :: filename

  integer :: icell

  filename = trim(data_dir2(imol))//"/maser_map.fits.gz"

  icell = cell_map(ri,zj,phik)

  do iTrans=1,nTrans
     iiTrans = indice_Trans(iTrans) ! Pas fondamental ici mais bon ...
     nu = tab_nLevel(icell,iTransUpper(iiTrans))
     nl = tab_nLevel(icell,iTransLower(iiTrans))


     ! Opacite et emissivite raie
     kap = (nl*fBlu(iiTrans) - nu*fBul(iiTrans))
     eps =  nu*fAul(iiTrans)

     if (kap < 0.) then
        lmaser = .true.
        ! inversion value (inversion population is > 1 )
        maser_map(icell,iTrans) = (nu * poids_stat_g(iTransLower(iiTrans))) / &
             (poids_stat_g(iTransUpper(iiTrans)) * nl)
        kap = 0.
     endif

     ! longueur de vol en AU, a multiplier par le profil de raie
     kappa_mol_o_freq(icell,iiTrans) = kap / Transfreq(iiTrans) * AU_to_m
     emissivite_mol_o_freq(icell,iiTrans) = eps /  Transfreq(iiTrans) * AU_to_m
  enddo

!  if ( (lmaser) .and. (ri==n_rad) .and. (zj==nz) ) then
!      write(*,*) "*************************************************"
!      write(*,*) "WARNING : There are some inversion populations"
!      write(*,*) "   --> forcing kappa = 0."
!      write(*,*) "Inversion values written to :"
!      write(*,*) trim(filename)
!      write(*,*) "Max. inversion value =", maxval(maser_map)
!      write(*,*) "*************************************************"
!      call cfitsWrite(trim(filename),maser_map,shape(maser_map))
!  endif


  if (ldouble_RT) then
     do iTrans=1,nTrans
        iiTrans = indice_Trans(iTrans) ! Pas fondamental ici mais bon ...
        nu = tab_nLevel2(icell,iTransUpper(iiTrans))
        nl = tab_nLevel2(icell,iTransLower(iiTrans))

        ! Opacite et emissivite raie
        kap = (nl*fBlu(iiTrans) - nu*fBul(iiTrans))
        eps =  nu*fAul(iiTrans)

        ! longueur de vol en AU, a multiplier par la profil de raie
        kappa_mol_o_freq2(icell,iiTrans) = kap / Transfreq(iiTrans) * AU_to_m
        emissivite_mol_o_freq2(icell,iiTrans) = eps /  Transfreq(iiTrans) * AU_to_m
     enddo
  endif ! ldouble_RT

  return

end subroutine opacite_mol_loc

!***********************************************************

subroutine init_dust_mol(imol)
  ! calcul les opacites et emissivites des cellules
  ! aux longueurs d'onde des raies moleculaires
  ! pour prendre en compte l'interaction radiative entre
  ! les 2 phases.
  ! C. Pinte
  ! 17/10/07
  ! TODO : gerer le cas ou l'albedo (et donc scattering) non negligeable

  implicit none

  integer, intent(in) :: imol
  integer :: iTrans, ri, zj, phik, p_lambda, icell
  real(kind=db) :: f, Jnu
  real :: T, wl, kap

  real(kind=db) :: cst_E

  real, parameter :: gas_dust = 100
  real, parameter :: delta_lambda = 0.025

  phik=1

  cst_E=2.0*hp*c_light**2

  p_n_lambda = n_lambda

  ! Reallocation des tableaux de proprietes de poussiere
 ! n_lambda =   mol(imol)%nTrans_raytracing ! opacites dust considerees cst sur le profil de raie
  n_lambda =   nTrans ! opacites dust considerees cst sur le profil de raie
  call realloc_dust_mol()

  if (ldust_mol) then
     ! Tableau de longeur d'onde
     do iTrans=1,nTrans
        tab_lambda(iTrans) = c_light/Transfreq(iTrans) * 1.0e6 ! en microns
        tab_lambda_sup(iTrans)= tab_lambda(iTrans)*delta_lambda
        tab_lambda_sup(iTrans)= tab_lambda(iTrans)/delta_lambda
        tab_delta_lambda(iTrans) = tab_lambda_sup(iTrans) - tab_lambda_inf(iTrans)
     enddo

     if (lbenchmark_water3) then ! opacite en loi de puissance
        write(*,*) "WARNING : hard-coded gas_dust =", gas_dust

        do iTrans=1,mol(imol)%nTrans_raytracing
           wl = tab_lambda(iTrans)

           ! Loi d'opacite (cm^2 par g de poussiere)
           if (wl > 250) then
              kap = 10. * (wl/250.)**(-2.0)
           else
              kap = 10. * (wl/250.)**(-1.3)
           endif

           ! Multiplication par densite
           ! AU_to_cm**2 car on veut kappa_abs_eg en AU-1
           do icell=1,n_cells
              kappa_abs_eg(icell,iTrans) =  kap * (densite_gaz(icell) * cm_to_m**3) * masse_mol_gaz / &
                   gas_dust / cm_to_AU
           enddo

        enddo ! iTrans

        ! Pas de scattering
        kappa(:,:) = kappa_abs_eg(:,:)

     else ! cas par defaut
        call init_indices_optiques()

        ! On n'est interesse que par les prop d'abs : pas besoin des matrices de mueller
        scattering_method=1 ; lscattering_method1 = .true. ; p_lambda = 1
        aniso_method = 2 ; lmethod_aniso1 = .false.

        lsepar_pola = .false.
        ltemp = .false.
        lmono = .true. ! equivalent au mode sed2

        ! On recalcule les proprietes optiques
        write(*,*) "Computing dust properties for", nTrans, "wavelength"
        do iTrans=1,nTrans
           call prop_grains(iTrans, p_lambda)
           call opacite(iTrans)
        enddo
     endif


     ! Changement d'unite : kappa en m-1 pour le TR dans les raies !!!!!!!
     ! Sera reconverti en AU-1 dans opacite_mol_loc
     !kappa = kappa * m_to_AU
     !kappa_abs_eg = kappa_abs_eg * m_to_AU

     ! calcul de l'emissivite de la poussiere
     do iTrans=1,nTrans
        f = Transfreq(iTrans)

        ! TODO : accelerer cette boucle via routine Bnu_disk (ca prend du tps ???)
        ! TODO : generaliser pour tous les types de grains (ca doit deja exister non ???)
        ! TODO : ca peut aussi servir pour faire du ray-tracing ds le continu 8-)
        do ri=1,n_rad
           bz : do zj=j_start,nz
              if (zj==0) cycle bz
              do phik=1, n_az
                 icell = cell_map(ri,zj,phik)
                 ! Interpolation champ de radiation en longeur d'onde
                 if (lProDiMo2mcfost) then
                    Jnu = interp(m2p%Jnu(ri,zj,:), m2p%wavelengths, real(tab_lambda(iTrans)))
                 else
                    Jnu = 0.0 ! todo : pour prendre en compte scattering
                 endif

                 T = Temperature(icell)
                 ! On ne fait que du scattering isotropique dans les raies pour le moment ...
                 emissivite_dust(icell,iTrans) = kappa_abs_eg(icell,iTrans) * Bnu(f,T) ! + kappa_sca(iTrans,ri,zj,phik) * Jnu
              enddo ! phik
           enddo bz ! zj
        enddo ! ri
     enddo ! itrans

  else ! .not.ldust_mol
     kappa = 0.0
     emissivite_dust = 0.0
  endif

  ! Deallocation les tableaux dont on a pas besoin pour la suite
  call clean_mem_dust_mol()

  return

end subroutine init_dust_mol

!***********************************************************

subroutine equilibre_LTE_mol()
  ! Calcul les niveaux d'une molecule dans le cas LTE
  ! Pour initialisation
  ! Calcule au passage le nombre total de mol dans chaque cellule
  ! reutilise par equilibre_rad_mol
  ! C. Pinte
  ! 18/03/07

  implicit none

  integer :: l, icell

  !$omp parallel &
  !$omp default(none) &
  !$omp private(l,icell) &
  !$omp shared(n_cells,nLevels,tab_nLevel,poids_stat_g,Transfreq,Tcin,densite_gaz,tab_abundance)
  !$omp do
  do icell=1, n_cells
     tab_nLevel(icell,1) = 1.0
     do l=2, nLevels
        ! Utilisation de la temperature de la poussiere comme temperature LTE
        tab_nLevel(icell,l) = tab_nLevel(icell,l-1) * poids_stat_g(l)/poids_stat_g(l-1) * &
             exp(- hp * Transfreq(l-1)/ (kb*Tcin(icell)))
     enddo
     ! Teste OK : (Nu*Bul) / (Nl*Blu) = exp(-hnu/kT)
     ! write(*,*) "Verif", i, j, tab_nLevel(i,j,l) * Bul(l-1) / (tab_nLevel(i,j,l-1) * Blu(l-1)) ,  exp(- hp * Transfreq(l-1)/ (kb*Temperature(i,j,1)))
     ! read(*,*)

     ! Normalisation
     tab_nLevel(icell,:) = densite_gaz(icell) * tab_abundance(icell) * tab_nLevel(icell,:)  / sum(tab_nLevel(icell,:))
  enddo!icell
  !$omp end do
  !$omp  end parallel

  ! write(*,*) real( (sum(masse) * g_to_kg * gas_dust / masse_mol_gaz) / (4.*pi/3. * (rout * AU_to_cm)**3 ) )
  ! write(*,*) (sum(masse) * g_to_kg * gas_dust / masse_mol_gaz) * abundance * fAul(1) * hp * transfreq(1)
  ! stop

  return

end subroutine equilibre_LTE_mol

!********************************************************************

subroutine equilibre_rad_mol_loc(id,ri,zj,phik)
  ! Calcul les populations sur les differents niveaux a l'aide de
  ! l'equilibre radiatif et collisionnel
  ! C. Pinte
  ! 01/07/07
  ! 11/10/07 : modif pour traiter simultanement 2 champs de radiation

  implicit none

  integer, intent(in) :: id, ri, zj, phik

  real :: Temp, Tdust
  real(kind=db), dimension(nLevels,nLevels) :: A, C
  real(kind=db), dimension(nLevels) :: B
  real(kind=db), dimension(nLevels) :: cTot
  real(kind=db) :: boltzFac, JJmol
  integer :: i, j, eq, n_eq
  integer :: itrans, l, k, iPart, icell
  real(kind=db) :: collEx, colldeEx

  if (ldouble_RT) then
     n_eq = 2
  else
     n_eq = 1
  endif

  ! Matrice d'excitations/desexcitations collisionnelles
  icell = cell_map(ri,zj,phik)
  Temp = Tcin(icell)
  nH2 = densite_gaz(icell) * cm_to_m**3

  ! Pour test
  Tdust = Temperature(icell)

  C = 0._db
  do iPart = 1, nCollPart
     do iTrans = 1, nCollTrans(iPart)
        k = iCollUpper(iPart, iTrans)
        l = iCollLower(iPart, iTrans)

        boltzFac =  exp(-abs(Level_energy(k)-Level_energy(l)) * ev_to_J / (kb*Temp))
        colldeEx = tau_collision(Temp, iPart, iTrans) * nH2

        collEx = colldeEx * boltzFac * poids_stat_g(k)/poids_stat_g(l)

        C(l, k) = C(l, k) + collEx
        C(k, l) = C(k, l) + colldeEx
     enddo ! iTrans
  enddo ! iPart

  ! Test
  !C = 0. ! OK : je retouve le niveau d'exitation par le Cmb en zone externe

  cTot = 0.
  do k = 1, nLevels
     do l = 1, nLevels
        cTot(k) = cTot(k) + C(k,l)
     enddo
  enddo


  ! Calcul des populations de niveaux en eq avec 1 ou 2 champs de radiation
  do eq=1,n_eq

     ! Matrice de transitions radiatives
     B = 0.0_db
     A = 0.0_db
     do iTrans = 1, nTrans ! toutes les transition ici : pas de iiTrans

        k = iTransUpper(iTrans)
        l = iTransLower(iTrans)

        if (eq==1) then
           JJmol = Jmol(iTrans,id)
        else
           JJmol = Jmol(iTrans,id)
        endif

        !write(*,*) "Jmol", iTrans, JJmol, Bnu(transfreq(iTrans),350.)

        !JJmol = Bnu(transfreq(iTrans),Tdust) !  OK, redonne Tex = 350.

        A(k,k) = A(k,k) + Bul(iTrans) * JJmol  + Aul(iTrans)
        A(l,l) = A(l,l) + Blu(iTrans) * JJmol
        A(k,l) = A(k,l) - Blu(iTrans) * JJmol
        A(l,k) = A(l,k) - Bul(iTrans) * JJmol  - Aul(iTrans)
     enddo

     ! Test
     !A=0.0 ! Ok : je retombe sur le LTE

     ! Matrice totale
     do i = 1, nLevels
        A(i,i) = A(i,i) + cTot(i)
        do j = 1, nLevels
           if (i /= j) then
              A(i,j) = A(i,j) - C(j, i)
           endif
        enddo
     enddo

     ! On remplace la derniere equations (qui est une CL des autres)
     ! par l'equation de conservation
     A(nLevels,:) = 1.0_db ! la somme de toutes les pops
     B(nLevels) = 1.0_db   ! = 1

     ! Resolution systeme matriciel par methode Gauss-Jordan
     call GaussSlv(A, B, nLevels)

     icell = cell_map(ri,zj,phik)
     if (eq==1) then
        tab_nLevel(icell,:) = densite_gaz(icell) * tab_abundance(icell) * B(:)
     else
        tab_nLevel2(icell,:) = densite_gaz(icell) * tab_abundance(icell) * B(:)
     endif

  enddo !n_eq

  return

end subroutine equilibre_rad_mol_loc

!***********************************************************

subroutine equilibre_othin_mol()
  ! Calcul les populations dans le cas optiquement mince
  ! equilibre avec le champ de radiation externe
  ! C. Pinte
  ! 11/10/07

  implicit none

  integer :: ri, zj, phik, id

  Jmol(:,:) = 0.0_db

  id =1 ! TODO : parallelisation

  ! Equilibre avec Cmb pour toutes les cellules
  do ri=1, n_rad
     do zj=1, nz
        do phik=1,n_az
           ! Le champ de radiation est egal au Cmb
           Jmol(:,id) = tab_Cmb_mol(:)

           ! Equilibre
           call equilibre_rad_mol_loc(id,ri,zj,phik)

        enddo !phik
     enddo !zj
  enddo !ri

  return

end subroutine equilibre_othin_mol

!***********************************************************

subroutine equilibre_othin_mol_pop2()
  ! Calcul les populations dans le cas optiquement mince
  ! equilibre avec le champ de radiation externe et sauvegarde
  ! dans la seconde population
  ! C. Pinte
  ! 15/10/07

  implicit none

  real(kind=db) :: tab_nLevel_tmp(n_cells,nLevels)  ! pas 3D
  logical :: ldouble_RT_tmp

  Jmol(:,:) = 0.0_db

  ! Par securite : sauvegarde population 1
  tab_nLevel_tmp(:,:) =  tab_nLevel(:,:)
  ldouble_RT_tmp = ldouble_RT
  ldouble_RT = .false.

  call equilibre_othin_mol()

  ! Initialisation de la population 2
  tab_nLevel2(:,:) = tab_nLevel(:,:)

  ! Restauration population 1
  tab_nLevel(:,:) =  tab_nLevel_tmp(:,:)
  ldouble_RT = ldouble_RT_tmp

  return

end subroutine equilibre_othin_mol_pop2

!***********************************************************

subroutine J_mol_loc(id,ri,zj,phik,n_rayons,ispeed)
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: id,ri, zj, phik, n_rayons
  integer, dimension(2), intent(in) :: ispeed

  real(kind=db), dimension(ispeed(1):ispeed(2)) :: etau, P, opacite, Snu
  real(kind=db) :: somme, J

  integer :: iTrans, iray, icell

  icell = cell_map(ri,zj,phik)

  Jmol(:,id) = 0.0_db

  do iTrans=1, nTrans ! Toutes les transitions ici : pas de iiTrans
     somme = 0.0_db
     J = 0.0_db
     do iray=1, n_rayons
        P(:) =  Doppler_P_x_freq(:,iray,id)
        opacite(:) = kappa_mol_o_freq(icell,iTrans) * P(:) + kappa(icell,iTrans)
        etau(:) = exp(-ds(iray,id) * opacite(:)) ! exp(-tau)

        Snu(:) = ( emissivite_mol_o_freq(icell,iTrans) * P(:) + &
             emissivite_dust(iTrans,icell) ) / (opacite(:) + 1.0e-30_db)

        J = J + sum( (I0(:,iTrans,iray,id) * etau(:) + Snu(:) * (1.0_db - etau(:))) * P(:))
        somme = somme + sum(P(:))
     enddo ! iray
     Jmol(iTrans,id) =  J / somme

  enddo ! iTrans

  ! Normalisation par le nombre de rayons utilises
!  Jmol(:,id) = Jmol(:,id) / n_rayons

  if (ldouble_RT) then
     Jmol2(:,id) = 0.0_db

     do iray=1, n_rayons
        do iTrans=1, nTrans ! Toutes les transitions ici : pas de iiTrans
           P(:) =  Doppler_P_x_freq(:,iray,id)
           opacite(:) = kappa_mol_o_freq2(icell,iTrans) * P(:) + kappa(icell,iTrans)
           etau(:) = exp(-ds(iray,id) * opacite(:)) ! exp(-tau)

           Snu(:) = ( emissivite_mol_o_freq2(icell,iTrans) * P(:) + emissivite_dust(icell,iTrans) ) &
                / (opacite(:) + 1.0e-30_db)
           J = sum( (I0(:,iTrans,iray,id) * etau(:) + Snu(:) * (1.0_db - etau(:))) &
                * P(:)) / sum(P(:))

           Jmol2(iTrans,id) = Jmol2(iTrans,id) + J
        enddo ! iTrans
     enddo ! iray

     ! Normalisation par le nombre de rayons utilises
     Jmol2(:,id) = Jmol2(:,id) / n_rayons
  endif !ldouble_RT

  return

end subroutine J_mol_loc

!***********************************************************

function v_proj(ri,zj,x,y,z,u,v,w) !
  ! Vitesse projete en 1 point d'une cellule
  ! C. Pinte
  ! 13/07/07

  implicit none

  real(kind=db) :: v_proj
  integer, intent(in) :: ri, zj
  real(kind=db), intent(in) :: x,y,z,u,v,w

  real(kind=db) :: vitesse, vx, vy, vz, norme, r
  integer :: icell

  icell = cell_map(ri,zj,1)
  vitesse = vfield(icell)

  if (linfall) then
     r = sqrt(x*x+y*y+z*z)
     !  if (lbenchmark_water2)  vitesse = -r  * 1e5 * AU_to_pc    ! TMP pour bench water2 : ca change rien !!!????
     if (r > tiny_db) then
        norme = 1.0_db/r
        vx = x * norme * vitesse
        vy = y * norme * vitesse
        vz = z * norme * vitesse
        v_proj = vx * u + vy * v + vz * w
     else
        v_proj = 0.0_db
     endif
  else if (lkeplerian) then
     r = sqrt(x*x+y*y)
     if (r > tiny_db) then
        norme = 1.0_db/r
        vx = -y * norme * vitesse
        vy = x * norme * vitesse
        v_proj = vx * u + vy * v
   else
        v_proj = 0.0_db
     endif
  else
     write(*,*) "Error : velocity field not defined"
     stop
  endif

  return

end function v_proj

!***********************************************************

real(kind=db) function dv_proj(ri,zj,x0,y0,z0,x1,y1,z1,u,v,w) !
  ! Differentiel de vitesse projete entre 2 points
  ! au sein d'une cellule
  ! C. Pinte
  ! 13/07/07

  implicit none

  integer, intent(in) :: ri, zj
  real(kind=db), intent(in) :: x0,y0,z0,x1,y1,z1,u,v,w

  real(kind=db) :: vitesse, vx0, vy0, vz0, vx1, vy1, vz1, norme
  integer :: icell

  icell = cell_map(ri,zj,1)
  vitesse = vfield(icell)

  if (linfall) then
     ! Champ de vitesse au point 0
     norme = 1.0_db/sqrt(x0*x0+y0*y0+z0*z0)
     vx0 = x0 * norme * vitesse
     vy0 = y0 * norme * vitesse
     vz0 = z0 * norme * vitesse

     ! Champ de vitesse au point 1
     norme = 1.0_db/sqrt(x1*x1+y1*y1+z1*z1)
     vx1 = x1 * norme * vitesse
     vy1 = y1 * norme * vitesse
     vz1 = z1 * norme * vitesse

     dv_proj = (vx1 - vx0) * u + (vy1 - vy0) * v + (vz1 - vz0) * w
  else if (lkeplerian) then
      ! Champ de vitesse au point 0
     norme = 1.0_db/sqrt(x0*x0+y0*y0)
     vx0 = -y0 * norme * vitesse
     vy0 = x0 * norme * vitesse

     ! Champ de vitesse au point 1
     norme = 1.0_db/sqrt(x1*x1+y1*y1)
     vx1 = -y1 * norme * vitesse
     vy1 = x1 * norme * vitesse

     dv_proj = (vx1 - vx0) * u + (vy1 - vy0) * v
  endif

  return

end function dv_proj

!***********************************************************

subroutine freeze_out()
  ! supprime une molecule si T < T_freeze_out
  ! C. Pinte
  ! 18/01/08

  implicit none

  integer :: i, j, k, icell

  write (*,*) "Freezing out of molecules"

  do k=1, n_az
     bz : do j=j_start, nz
        if (j==0) cycle bz
        do i=1, n_rad
           icell = cell_map(i,j,k)
           if (Temperature(icell) < T_freeze_out)  tab_abundance(icell) = 0.
        enddo
     enddo bz
  enddo

  return

end subroutine freeze_out

!***********************************************************

end module molecules
