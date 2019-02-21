module molecular_emission

  use parametres
  use temperature
  use constantes
  use grid
  use density
  use dust_prop

  implicit none
  save

  logical :: ldouble_RT

  real, dimension(:), allocatable :: Level_energy
  integer, dimension(:), allocatable ::  j_qnb
  ! g est dp car les calculs utilisant g sont en dp
  real(kind=dp), dimension(:), allocatable :: poids_stat_g
  integer :: nTrans_tot

  real(kind=dp), dimension(:), allocatable :: Aul, Blu, Bul, fAul, fBlu, fBul, transfreq
  integer, dimension(:), allocatable :: itransUpper, itransLower
  integer :: nCollPart
  character(len=512), dimension(:), allocatable :: collBetween
  integer, dimension(:), allocatable :: nCollTrans, nCollTemps
  real, dimension(:,:), allocatable :: collTemps

  integer, dimension(:,:), pointer :: iCollUpper, iCollLower
  real, dimension(:,:,:), pointer :: collRates

  real, dimension(:), allocatable :: Tcin ! Temperature cinetique
  real :: correct_Tgas
  logical :: lcorrect_Tgas

  real :: nH2, masse_mol
  ! masse_mol_gaz sert uniquement pour convertir masse disque en desnite de particule
  real(kind=dp), dimension(:,:), allocatable :: kappa_mol_o_freq, kappa_mol_o_freq2 ! n_cells, nTrans
  real(kind=dp), dimension(:,:), allocatable :: emissivite_mol_o_freq,  emissivite_mol_o_freq2 ! n_cells, nTrans
  real, dimension(:,:), allocatable :: tab_nLevel, tab_nLevel2, tab_nLevel_old ! n_cells, nLevels

  real, dimension(:), allocatable :: v_turb, v_line ! n_cells

  real ::  vitesse_turb, dv, dnu
  integer, parameter :: n_largeur_Doppler = 15
  real(kind=dp), dimension(:), allocatable :: tab_v ! n_speed

  real(kind=dp), dimension(:,:), allocatable :: ds
  real(kind=dp), dimension(:,:,:,:), allocatable :: I0, I02 ! nSpeed,nTrans,iray,ncpus
  real(kind=dp), dimension(:,:,:), allocatable :: I0c ! Intensite dans le continu: nTrans,iray,ncpus
  real(kind=dp), dimension(:,:,:), allocatable :: Doppler_P_x_freq

  real(kind=dp), dimension(:,:), allocatable :: Jmol, Jmol2 ! nTrans, n_cpu
  real(kind=dp), dimension(:), allocatable :: tab_Cmb_mol ! nTrans

  logical :: linfall, lkeplerian, lcylindrical_rotation
  real :: chi_infall

  real(kind=dp), dimension(:), allocatable :: deltaVmax ! n_cells
  real(kind=dp), dimension(:), allocatable :: tab_dnu_o_freq ! n_cells
  real(kind=dp), dimension(:), allocatable :: norme_phiProf_m1, sigma2_phiProf_m1 ! n_cells

  real, dimension(:), allocatable :: tab_abundance ! n_cells
  logical, dimension(:), allocatable :: lcompute_molRT ! n_cells

  logical ::  lfreeze_out, lphoto_dissociation, lphoto_desorption
  real :: T_freeze_out, freeze_out_depletion

  real(kind=dp), dimension(:,:,:,:), allocatable ::  origine_mol ! nv, nTrans, n_cells, nb_proc

  integer :: RT_line_method, n_molecules

  type molecule
     integer :: n_speed_rt, n_speed_center_rt, n_extraV_rt, nTrans_raytracing, iLevel_max
     real :: vmax_center_rt, extra_deltaV_rt, abundance
     logical :: lcst_abundance, lline
     character(len=512) :: abundance_file, filename
     character(len=32) :: name
     integer, dimension(100) :: indice_Trans_rayTracing
  end type molecule

  type(molecule), dimension(:), allocatable :: mol

  real(kind=dp), dimension(:), allocatable :: tab_speed_rt

  real, dimension(:,:), allocatable :: maser_map ! n_cells, n_trans

  real(kind=dp), dimension(:,:), allocatable :: emissivite_dust ! emissivite en SI (pour mol)

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

subroutine init_Doppler_profiles(imol)
  ! Creation des constantes pour les profiles Doppler dans chaque cellule
  ! a vitesse systematique nulle + tableau de vitesses
  ! C. Pinte
  ! 18/03/07

  implicit none

  integer, intent(in) :: imol

  real(kind=dp) :: sigma2, sigma2_m1, vmax
  integer :: icell, iv, n_speed

  n_speed = mol(imol)%n_speed_rt

  do icell=1, n_cells
     ! Utilisation de la temperature LTE de la poussiere comme temperature cinetique
     ! WARNING : c'est pas un sigma mais un delta, cf Cours de Boisse p47
     ! Consistent avec benchmark
     sigma2 =  2.0_dp * (kb*Tcin(icell) / (masse_mol * g_to_kg)) + v_turb(icell)**2
     v_line(icell) = sqrt(sigma2)

     !  write(*,*) "FWHM", sqrt(sigma2 * log(2.)) * 2.  ! Teste OK bench water 1
     sigma2_m1 = 1.0_dp / sigma2
     sigma2_phiProf_m1(icell) = sigma2_m1
     ! phi(nu) et non pas phi(v) donc facteur c_light et il manque 1/f0
     ! ATTENTION : il ne faut pas oublier de diviser par la freq apres
     norme_phiProf_m1(icell) = c_light / sqrt(pi * sigma2)

     ! Echantillonage du profil de vitesse dans la cellule
     ! 2.15 correspond a l'enfroit ou le profil de la raie faut 1/100 de
     ! sa valeur au centre : exp(-2.15^2) = 0.01
     vmax = sqrt(sigma2)
     tab_dnu_o_freq(icell) = largeur_profile * vmax / (real(n_speed))
     deltaVmax(icell) = largeur_profile * vmax !* 2.0_dp  ! facteur 2 pour tirage aleatoire
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
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: phiProf

  real(kind=dp) :: norme_m1, sigma2_m1

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
! ---   real(kind=dp), dimension(ispeed(1):ispeed(2)),intent(in) :: tab_speed
! ---
! ---   real(kind=dp), dimension(ispeed(1):ispeed(2),nTrans) :: Cmb
! ---
! ---
! ---   integer :: i, iv, iTrans
! ---   real(kind=dp) :: cst, cst_nu, nu
! ---
! ---   cst = 2.0_dp*hp/c_light**2
! ---
! ---   ! On ne calcule le Cmb qu'aux frequences specifiees
! ---   do i=1,nTrans
! ---      iTrans = indice_Trans(i)
! ---      do iv=ispeed(1),ispeed(2)
! ---         nu = Transfreq(iTrans) * (1.0_dp + tab_speed(iv)/c_light)
! ---
! ---         cst_nu = (hp * nu) / (kb * T_min)
! ---         if (cst_nu > 100._dp) then
! ---            Cmb(iv,i) = 0.0_dp
! ---         else
! ---            Cmb(iv,i) = cst * nu**3 / (exp(cst_nu)-1.0_dp)
! ---         endif
! ---      enddo ! iv
! ---   enddo ! i
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
  real(kind=dp) :: cst, cst_nu, nu

  cst = 2.0_dp*hp/c_light**2

  do iTrans=1,nTrans_tot
     nu = Transfreq(iTrans)

     cst_nu = (hp * nu) / (kb * T_Cmb)
     if (cst_nu > 100._dp) then
        tab_Cmb_mol(iTrans) = 0.0_dp
     else
        tab_Cmb_mol(iTrans) = cst * nu**3 / (exp(cst_nu)-1.0_dp)
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
  integer :: icell

  do icell=1,n_cells
     call opacite_mol_loc(icell,imol)
  enddo ! icell

  return

end subroutine opacite_mol

!***********************************************************

subroutine opacite_mol_loc(icell,imol)
  ! Calcule la fonction source dans chaque cellule
  ! etant donne les populations des niveaux
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: icell, imol

  integer :: iTrans
  real(kind=dp) :: nu, nl, kap, eps

  logical, save :: lmaser = .false.

  character(len=128) :: filename

  filename = trim(data_dir2(imol))//"/maser_map.fits.gz"

  do iTrans=1,nTrans_tot
     nu = tab_nLevel(icell,iTransUpper(iTrans))
     nl = tab_nLevel(icell,iTransLower(iTrans))

     ! Opacite et emissivite raie
     kap = (nl*fBlu(iTrans) - nu*fBul(iTrans))
     eps =  nu*fAul(iTrans)

     if (kap < 0.) then
        lmaser = .true.
        ! inversion value (inversion population is > 1 )
        maser_map(icell,iTrans) = (nu * poids_stat_g(iTransLower(iTrans))) / &
             (poids_stat_g(iTransUpper(iTrans)) * nl)
        kap = 0.
     endif

     ! longueur de vol en AU, a multiplier par le profil de raie
     kappa_mol_o_freq(icell,iTrans) = kap / Transfreq(iTrans) * AU_to_m
     emissivite_mol_o_freq(icell,iTrans) = eps /  Transfreq(iTrans) * AU_to_m
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
     do iTrans=1,nTrans_tot
        nu = tab_nLevel2(icell,iTransUpper(iTrans))
        nl = tab_nLevel2(icell,iTransLower(iTrans))

        ! Opacite et emissivite raie
        kap = (nl*fBlu(iTrans) - nu*fBul(iTrans))
        eps =  nu*fAul(iTrans)

        ! longueur de vol en AU, a multiplier par la profil de raie
        kappa_mol_o_freq2(icell,iTrans) = kap / Transfreq(iTrans) * AU_to_m
        emissivite_mol_o_freq2(icell,iTrans) = eps /  Transfreq(iTrans) * AU_to_m
     enddo
  endif ! ldouble_RT

  return

end subroutine opacite_mol_loc

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
     ! write(*,*) "Verif", i, j, tab_nLevel(i,j,l) * Bul(l-1) / (tab_nLevel(i,j,l-1) * Blu(l-1)) ,  exp(- hp * Transfreq(l-1)/ (kb*Tdust(i,j,1)))
     ! read(*,*)

     ! Normalisation
     tab_nLevel(icell,:) = densite_gaz(icell) * tab_abundance(icell) * tab_nLevel(icell,:)  / sum(tab_nLevel(icell,:))
  enddo!icell
  !$omp end do
  !$omp  end parallel

  ! write(*,*) real( (sum(masse) * g_to_kg * gas_dust / masse_mol_gaz) / (4.*pi/3. * (rout * AU_to_cm)**3 ) )
  ! write(*,*) (sum(masse) * g_to_kg * gas_dust / masse_mol_gaz) * abundance * fAul(1) * hp * transfreq(1)

  return

end subroutine equilibre_LTE_mol

!********************************************************************

subroutine equilibre_rad_mol_loc(id,icell)
  ! Calcul les populations sur les differents niveaux a l'aide de
  ! l'equilibre radiatif et collisionnel
  ! C. Pinte
  ! 01/07/07
  ! 11/10/07 : modif pour traiter simultanement 2 champs de radiation

  implicit none

  integer, intent(in) :: id, icell

  real :: Temp
  real(kind=dp), dimension(nLevels,nLevels) :: A, C
  real(kind=dp), dimension(nLevels) :: B
  real(kind=dp), dimension(nLevels) :: cTot
  real(kind=dp) :: boltzFac, JJmol
  integer :: i, j, eq, n_eq
  integer :: itrans, l, k, iPart
  real(kind=dp) :: collEx, colldeEx

  if (ldouble_RT) then
     n_eq = 2
  else
     n_eq = 1
  endif

  ! Matrice d'excitations/desexcitations collisionnelles
  Temp = Tcin(icell)
  nH2 = densite_gaz(icell) * cm_to_m**3

  C = 0._dp
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
     B = 0.0_dp
     A = 0.0_dp
     do iTrans = 1, nTrans_tot
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
     A(nLevels,:) = 1.0_dp ! la somme de toutes les pops
     B(nLevels) = 1.0_dp   ! = 1

     ! Resolution systeme matriciel par methode Gauss-Jordan
     call GaussSlv(A, B, nLevels)

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

  integer :: icell, id

  Jmol(:,:) = 0.0_dp

  id = 1 ! TODO : parallelisation

  ! Equilibre avec Cmb pour toutes les cellules
  do icell=1, n_cells
     ! Le champ de radiation est egal au Cmb
     Jmol(:,id) = tab_Cmb_mol(:)

     ! Equilibre
     call equilibre_rad_mol_loc(id,icell)
  enddo !icell

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

  real(kind=dp) :: tab_nLevel_tmp(n_cells,nLevels)  ! pas 3D
  logical :: ldouble_RT_tmp

  Jmol(:,:) = 0.0_dp

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

subroutine J_mol_loc(id,icell,n_rayons,ispeed)
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: id, icell, n_rayons
  integer, dimension(2), intent(in) :: ispeed

  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: etau, P, opacite, Snu
  real(kind=dp) :: somme, J

  integer :: iTrans, iray

  Jmol(:,id) = 0.0_dp

  do iTrans=1, nTrans_tot
     somme = 0.0_dp
     J = 0.0_dp
     do iray=1, n_rayons
        P(:) =  Doppler_P_x_freq(:,iray,id)
        opacite(:) = kappa_mol_o_freq(icell,iTrans) * P(:) + kappa(icell,iTrans)
        etau(:) = exp(-ds(iray,id) * opacite(:)) ! exp(-tau)

        Snu(:) = ( emissivite_mol_o_freq(icell,iTrans) * P(:) + &
             emissivite_dust(iTrans,icell) ) / (opacite(:) + 1.0e-30_dp)

        J = J + sum( (I0(:,iTrans,iray,id) * etau(:) + Snu(:) * (1.0_dp - etau(:))) * P(:))
        somme = somme + sum(P(:))
     enddo ! iray
     Jmol(iTrans,id) =  J / somme

  enddo ! iTrans

  ! Normalisation par le nombre de rayons utilises
!  Jmol(:,id) = Jmol(:,id) / n_rayons

  if (ldouble_RT) then
     Jmol2(:,id) = 0.0_dp

     do iray=1, n_rayons
        do iTrans=1, nTrans_tot
           P(:) =  Doppler_P_x_freq(:,iray,id)
           opacite(:) = kappa_mol_o_freq2(icell,iTrans) * P(:) + kappa(icell,iTrans)
           etau(:) = exp(-ds(iray,id) * opacite(:)) ! exp(-tau)

           Snu(:) = ( emissivite_mol_o_freq2(icell,iTrans) * P(:) + emissivite_dust(icell,iTrans) ) &
                / (opacite(:) + 1.0e-30_dp)
           J = sum( (I0(:,iTrans,iray,id) * etau(:) + Snu(:) * (1.0_dp - etau(:))) &
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

function v_proj(icell,x,y,z,u,v,w) !
  ! Vitesse projete en 1 point d'une cellule
  ! C. Pinte
  ! 13/07/07

  implicit none

  real(kind=dp) :: v_proj
  integer, intent(in) :: icell
  real(kind=dp), intent(in) :: x,y,z,u,v,w

  real(kind=dp) :: vitesse, vx, vy, vz, norme, r

  if (lVoronoi) then
     vx = Voronoi(icell)%vxyz(1)
     vy = Voronoi(icell)%vxyz(2)
     vz = Voronoi(icell)%vxyz(3)

     v_proj = vx * u + vy * v + vz * w
  else
     if (ldensity_file) then
        vx = vfield_x(icell) ; vy = vfield_y(icell) ; vz = vfield_z(icell)
        v_proj = vx * u + vy * v + vz * w
     else ! Using analytical velocity field
        vitesse = vfield(icell)

        if (lkeplerian) then
           r = sqrt(x*x+y*y)
           if (r > tiny_dp) then
              norme = 1.0_dp/r
              vx = -y * norme * vitesse
              vy = x * norme * vitesse
              vz = 0.
              if (linfall) then ! Adding extra velocity
                 r = sqrt(x*x+y*y+z*z) ; norme = 1.0_dp/r
                 vx = vx - chi_infall * x * norme * vitesse
                 vy = vy - chi_infall * y * norme * vitesse
                 vz = vz - chi_infall * z * norme * vitesse
                 v_proj = vx * u + vy * v + vz * w
              else
                 v_proj = vx * u + vy * v
              endif
           else
              v_proj = 0.0_dp
           endif
        else if (linfall) then
           r = sqrt(x*x+y*y+z*z)
           !  if (lbenchmark_water2)  vitesse = -r  * 1e5 * AU_to_pc    ! TMP pour bench water2 : ca change rien !!!????
           if (r > tiny_dp) then
              norme = 1.0_dp/r
              vx = x * norme * vitesse
              vy = y * norme * vitesse
              vz = z * norme * vitesse
              v_proj = vx * u + vy * v + vz * w
           else
              v_proj = 0.0_dp
           endif
        else
           call error("velocity field not defined")
        endif
     endif ! ldensity_file
  endif

  return

end function v_proj

!***********************************************************

real(kind=dp) function dv_proj(icell,x0,y0,z0,x1,y1,z1,u,v,w) !
  ! Differentiel de vitesse projete entre 2 points
  ! au sein d'une cellule
  ! C. Pinte
  ! 13/07/07

  implicit none

  integer, intent(in) :: icell
  real(kind=dp), intent(in) :: x0,y0,z0,x1,y1,z1,u,v,w

  real(kind=dp) :: vitesse, vx0, vy0, vz0, vx1, vy1, vz1, norme

  vitesse = vfield(icell)

  if (linfall) then
     ! Champ de vitesse au point 0
     norme = 1.0_dp/sqrt(x0*x0+y0*y0+z0*z0)
     vx0 = x0 * norme * vitesse
     vy0 = y0 * norme * vitesse
     vz0 = z0 * norme * vitesse

     ! Champ de vitesse au point 1
     norme = 1.0_dp/sqrt(x1*x1+y1*y1+z1*z1)
     vx1 = x1 * norme * vitesse
     vy1 = y1 * norme * vitesse
     vz1 = z1 * norme * vitesse

     dv_proj = (vx1 - vx0) * u + (vy1 - vy0) * v + (vz1 - vz0) * w
  else if (lkeplerian) then
      ! Champ de vitesse au point 0
     norme = 1.0_dp/sqrt(x0*x0+y0*y0)
     vx0 = -y0 * norme * vitesse
     vy0 = x0 * norme * vitesse

     ! Champ de vitesse au point 1
     norme = 1.0_dp/sqrt(x1*x1+y1*y1)
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
  ! 6/11/16 : add a depletion factor

  implicit none

  real, parameter :: threshold_CD = 0.8 * 1.59e21  / (cm_to_m**2) ! m^-2
  ! 1e4x photodissociation from Qi et al 2011, ajsuted by hand for IM Lupi

  integer :: icell
  real(kind=dp) :: CD
  logical :: ldeplete

  write (*,*) "Freezing-out of molecules"

  do icell=1,n_cells
     ldeplete = .false.
     if (Tdust(icell) < T_freeze_out)  then
        if (lphoto_desorption) then
           CD = compute_vertical_CD(icell)
           if (CD < threshold_CD) then
              ldeplete = .false. ! photo-desorption
           else
              ldeplete = .true. ! freeze-out
           endif
        else
           ldeplete = .true. ! freeze-out
        endif
        if (ldeplete) tab_abundance(icell) = tab_abundance(icell) * freeze_out_depletion
     endif
  enddo

  return

end subroutine freeze_out

!***********************************************************

subroutine photo_dissociation()
  ! supprime une molecule si cd < value
  ! C. Pinte
  ! 6/11/16

  implicit none

  integer :: icell
  real(kind=dp) :: CD

  real, parameter :: threshold_CD = 0.8 * 1.59e21 * 0.65  / (cm_to_m**2) ! m^-2 ! Value from Qi et al 2011
  ! It makes sense only for constant dust --> needs to be updated
  real, parameter :: photo_dissocation_depletion = 1.e-6

  write (*,*) "Photo-dissociating molecules"!, threshold_CD

  do icell=1, n_cells
     CD = compute_vertical_CD(icell)
     if (CD < threshold_CD) then
        tab_abundance(icell) = tab_abundance(icell) * photo_dissocation_depletion
     endif
  enddo ! icell

  return

end subroutine photo_dissociation

!***********************************************************

function compute_vertical_CD(icell) result(CD)

  integer, intent(in) :: icell
  integer :: icell0, next_cell, previous_cell
  real(kind=dp) :: CD, x0, y0, z0, x1, y1, z1, u,v,w, l, l_contrib, l_void_before

  if (lVoronoi) then
     x1 = Voronoi(icell)%xyz(1)
     y1 = Voronoi(icell)%xyz(2)
     z1 = Voronoi(icell)%xyz(3)
  else
     x1 = r_grid(icell) * cos(phi_grid(icell))
     y1 = r_grid(icell) * sin(phi_grid(icell))
     z1 = z_grid(icell)
  endif

  u = 0.0 ; v = 0.0 ;
  if (z1 >= 0) then
     w = 1.0
  else
     w = -1.0
  endif

  next_cell = icell
  icell0 = 0
  CD = 0.0
  ! The test to check if we exit the model volume depends on the grid:
  !  - next_cell > n_cells for cylindrical and spherical grids
  !  - next_cell < 0 for Voronoi mesh (ie next_cell is a wall)
  do while((next_cell > 0).and.(next_cell <= n_cells))
     previous_cell = icell0
     icell0 = next_cell
     x0 = x1 ; y0 = y1 ; z0 = z1
     call cross_cell(x0,y0,z0, u,v,w,  icell0, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)
     CD = CD + (l_contrib * AU_to_m) * densite_gaz(icell) ! part.m^-2
  enddo

  return

end function compute_vertical_CD

!***********************************************************

end module molecular_emission
