module molecular_emission

  use parameters
  use temperature
  use constants
  use grid
  use density
  use dust_prop

  implicit none
  save

  logical :: ldouble_RT, lrovib

  real, dimension(:), allocatable :: Level_energy
  character(len=12), dimension(:), allocatable ::  j_qnb, v_qnb

  ! g is dp because calculations using g are in double precision
  real(kind=dp), dimension(:), allocatable :: stat_weight_g
  integer :: nTrans_tot

  real(kind=dp), dimension(:), allocatable :: Aul, Blu, Bul, fAul, fBlu, fBul, transfreq
  integer, dimension(:), allocatable :: itransUpper, itransLower
  integer :: nCollPart
  character(len=512), dimension(:), allocatable :: collBetween
  integer, dimension(:), allocatable :: nCollTrans, nCollTemps
  real, dimension(:,:), allocatable :: collTemps

  integer, dimension(:,:), pointer :: iCollUpper, iCollLower
  real, dimension(:,:,:), pointer :: collRates

  real, dimension(:), allocatable :: Tcin ! Kinetic temperature
  real :: correct_Tgas
  logical :: lcorrect_Tgas

  real :: nH2, masse_mol ! masse_mol is the mass of the current molecule
  real(kind=dp), dimension(:,:), allocatable :: kappa_mol_o_freq, kappa_mol_o_freq2 ! n_cells, nTrans
  real(kind=dp), dimension(:,:), allocatable :: emissivite_mol_o_freq,  emissivite_mol_o_freq2 ! n_cells, nTrans
  real, dimension(:,:), allocatable :: tab_nLevel, tab_nLevel2 ! n_cells, nLevels

  real, dimension(:), allocatable :: v_turb2, dv_line ! n_cells, v_turb2 is (v_turb in m/s)**2

  logical :: lvturb_in_cs
  real ::  vitesse_turb, dv, dnu, v_syst
  character(len=8) :: v_turb_unit
  integer, parameter :: n_doppler_width = 15
  real(kind=dp), dimension(:), allocatable :: tab_v ! n_speed

  real(kind=dp), dimension(:,:), allocatable :: ds!, gradv !local velocity difference between two cells in a specific direction
  real(kind=dp), dimension(:,:,:,:), allocatable :: I0, I02 ! nSpeed,nTrans,iray,ncpus
  real(kind=dp), dimension(:,:,:), allocatable :: I0c ! Continuum intensity: nTrans,iray,ncpus
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
  real :: T_freeze_out, freeze_out_depletion, photodissociation_factor

  real(kind=dp), dimension(:,:,:,:), allocatable ::  origine_mol ! nv, nTrans, n_cells, nb_proc

  integer :: RT_line_method, n_molecules

  type molecule
     integer :: n_speed_rt, n_speed_center_rt, n_extraV_rt, nTrans_raytracing, iLevel_max
     real :: vmin_center_rt, vmax_center_rt, extra_deltaV_rt, abundance, molecularWeight
     logical :: lcst_abundance, lline, l_sym_ima
     character(len=512) :: abundance_file, filename
     character(len=32) :: name
     integer, dimension(100) :: index_trans_ray_tracing

     integer ::  iTrans_min, iTrans_max, level_min, level_max ! transitions and levels used in practice
  end type molecule

  type(molecule), dimension(:), allocatable :: mol

  real(kind=dp), dimension(:), allocatable :: tab_speed_rt

  !real, dimension(:,:), allocatable :: maser_map ! n_cells, n_trans

  real(kind=dp), dimension(:,:), allocatable :: emissivite_dust ! emissivity in SI units (for molecules)

  contains

real function tau_collision(temperature, iPart, iTrans)
  ! Linear interpolation of collision rates with bisection
  ! C. Pinte
  ! 25/06/07

  implicit none

  real, intent(in) :: temperature
  integer, intent(in) :: iPart, iTrans

  real :: frac
  integer :: k, kmin, kmax

  if (nCollTemps(iPart) == 1) then ! Pas d'interpolation
     tau_collision = collRates(iPart, iTrans, 1)
  else if (temperature < collTemps(iPart,1)) then ! cubic extrapolation
     tau_collision = collRates(iPart, iTrans, 1) * sqrt(temperature/collTemps(iPart,1))
  else ! linear interpolation
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
  ! Create constants for Doppler profiles in each cell
  ! at zero systematic velocity + velocity array
  ! C. Pinte
  ! 18/03/07

  implicit none

  integer, intent(in) :: imol

  real(kind=dp) :: sigma2, sigma2_m1
  integer :: icell, n_speed

  n_speed = mol(imol)%n_speed_rt

  do icell=1, n_cells
     ! Using the LTE dust temperature as the kinetic temperature
     ! WARNING : this is a delta, not a sigma; see Boisse lecture notes p47
     ! Consistent avec benchmark
     sigma2 =  2.0_dp * (kb*Tcin(icell) / (mol(imol)%molecularWeight * mH  * g_to_kg)) + v_turb2(icell)
     dv_line(icell) = sqrt(sigma2)

     !  write(*,*) "FWHM", sqrt(sigma2 * log(2.)) * 2.  ! Teste OK bench water 1
     sigma2_m1 = 1.0_dp / sigma2
     sigma2_phiProf_m1(icell) = sigma2_m1
     ! phi(nu) not phi(v), hence factor c_light and missing 1/f0
     ! WARNING: must divide by the frequency afterwards
     norme_phiProf_m1(icell) = c_light / sqrt(pi * sigma2)

     ! Sampling of the velocity profile in the cell
     ! 2.15 corresponds to where the line profile reaches 1/100 of
     ! sa valeur au centre : exp(-2.15^2) = 0.01
     tab_dnu_o_freq(icell) = profile_width * dv_line(icell) / (real(n_speed))
     deltaVmax(icell) = profile_width * dv_line(icell) !* 2.0_dp  ! factor 2 for random sampling
  enddo !icell

  return

end subroutine init_Doppler_profiles

!***********************************************************

function phiProf(icell,ispeed,tab_speed)
  ! Returns the local line profile at a cell (ri, zj)
  ! on the velocity array tab_speed
  ! Must divide by the transition frequency for correct normalisation!!!
  ! C. Pinte
  ! 13/07/07

  implicit none

  integer, dimension(2), intent(in) :: ispeed
  integer, intent(in) :: icell
  real(kind=dp), dimension(ispeed(1):ispeed(2)), intent(in) :: tab_speed

  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: phiProf

  real(kind=dp) :: norme_m1, sigma2_m1

  norme_m1 = norme_phiProf_m1(icell) ! WARNING: frequency is missing here!!!!
  sigma2_m1 = sigma2_phiProf_m1(icell)

  phiProf(:) =  norme_m1 * exp(- sigma2_m1 * tab_speed(:)**2)

  return

end function phiProf

!***********************************************************

! --- function Cmb(ispeed,tab_speed)
! --- ! Planck's law at the different frequencies of the
! --- ! studied molecule
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
! ---   ! We only compute the CMB at the specified frequencies
! ---   do i=1,nTrans
! ---      iTrans = index_trans(i)
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
  ! Planck's law at the different frequencies of the
  ! studied molecule
  ! Bnu in SI : W.m-2.Hz-1.sr-1
  ! Without taking into account the velocity field << c
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
  ! Computes the source function in each cell
  ! given the level populations
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
  ! Computes the source function in each cell
  ! given the level populations
  ! C. Pinte
  ! 01/07/07

  implicit none

  integer, intent(in) :: icell, imol

  integer :: iTrans
  real(kind=dp) :: nu, nl, kap, eps

  do iTrans=mol(imol)%iTrans_min,mol(imol)%iTrans_max
     nu = tab_nLevel(icell,iTransUpper(iTrans))
     nl = tab_nLevel(icell,iTransLower(iTrans))

     ! opacity and line emissivity
     kap = (nl*fBlu(iTrans) - nu*fBul(iTrans))
     eps =  nu*fAul(iTrans)

     if (kap < 0.) then
        kap = 0.
        !lmaser = .true.
        ! inversion value (inversion population is > 1 )
        !maser_map(icell,iTrans) = (nu * stat_weight_g(iTransLower(iTrans))) / &
        !     (stat_weight_g(iTransUpper(iTrans)) * nl)
     endif

     ! flight length in AU, to be multiplied by the line profile
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

        ! opacity and line emissivity
        kap = (nl*fBlu(iTrans) - nu*fBul(iTrans))
        eps =  nu*fAul(iTrans)

        ! flight length in AU, to be multiplied by the line profile
        kappa_mol_o_freq2(icell,iTrans) = kap / Transfreq(iTrans) * AU_to_m
        emissivite_mol_o_freq2(icell,iTrans) = eps /  Transfreq(iTrans) * AU_to_m
     enddo
  endif ! ldouble_RT

  return

end subroutine opacite_mol_loc

!***********************************************************

subroutine equilibre_LTE_mol(imol)
  ! Computes the molecule levels in the LTE case
  ! For initialization
  ! Computes the total number of molecules in each cell simultaneously
  ! reused by equilibre_rad_mol
  ! C. Pinte
  ! 18/03/07

  implicit none

  integer, intent(in) :: imol

  integer :: l, icell, lmin, lmax

  real, dimension(nLevels) :: pop_levels ! local population levels to a cell
  real(kind=dp) :: norm

  lmin=mol(imol)%level_min
  lmax=mol(imol)%level_max

  !$omp parallel &
  !$omp default(none) &
  !$omp private(l,icell,pop_levels,norm) &
  !$omp shared(lmin,lmax,n_cells,nLevels,tab_nLevel,stat_weight_g,Transfreq,Tcin,gas_density,tab_abundance)
  !$omp do
  do icell=1, n_cells
     pop_levels(1) = 1.0
     do l=2, nLevels
        ! Use of the dust temperature as LTE temperature
        pop_levels(l) = pop_levels(l-1) * stat_weight_g(l)/stat_weight_g(l-1) * &
             exp(- hp * Transfreq(l-1)/ (kb*Tcin(icell)))
     enddo
     ! Test OK : (Nu*Bul) / (Nl*Blu) = exp(-hnu/kT)
     ! write(*,*) "Verif", i, j, tab_nLevel(i,j,l) * Bul(l-1) / (tab_nLevel(i,j,l-1) * Blu(l-1)) ,  exp(- hp * Transfreq(l-1)/ (kb*Tdust(i,j,1)))
     ! read(*,*)

     ! Normalisation
     norm = (gas_density(icell) * tab_abundance(icell)  / sum(pop_levels(:)))
     pop_levels(:) = pop_levels(:) * norm

     ! Saving a fraction of the array for the transition we consider
     do l=lmin,lmax
        tab_nLevel(icell,l) = pop_levels(l)
     enddo
  enddo!icell
  !$omp end do
  !$omp  end parallel

  ! write(*,*) real( (sum(dust_mass) * g_to_kg * gas_dust / mu_mH) / (4.*pi/3. * (rout * AU_to_cm)**3 ) )
  ! write(*,*) (sum(dust_mass) * g_to_kg * gas_dust / mu_mH) * abundance * fAul(1) * hp * transfreq(1)

  return

end subroutine equilibre_LTE_mol

!********************************************************************

subroutine equilibre_rad_mol_loc(id,icell)
  ! Computes the populations on the different levels using
  ! radiative and collisional equilibrium
  ! C. Pinte  16/2/06
  !
  !
  ! 11/10/07 : modification to handle 2 radiation fields simultaneously

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

  ! Collisional excitation/de-excitation matrix
  Temp = Tcin(icell)
  nH2 = gas_density(icell) * cm_to_m**3

  C = 0._dp
  do iPart = 1, nCollPart
     do iTrans = 1, nCollTrans(iPart)
        k = iCollUpper(iPart, iTrans)
        l = iCollLower(iPart, iTrans)

        boltzFac =  exp(-abs(Level_energy(k)-Level_energy(l)) * ev_to_J / (kb*Temp))
        colldeEx = tau_collision(Temp, iPart, iTrans) * nH2

        collEx = colldeEx * boltzFac * stat_weight_g(k)/stat_weight_g(l)

        C(l, k) = C(l, k) + collEx
        C(k, l) = C(k, l) + colldeEx
     enddo ! iTrans
  enddo ! iPart

  ! Test
  !C = 0. ! OK : I retrieve the excitation level by the CMB in the outer zone

  cTot = 0.
  do k = 1, nLevels
     do l = 1, nLevels
        cTot(k) = cTot(k) + C(k,l)
     enddo
  enddo


  ! Computes the level populations in equilibrium with 1 or 2 radiation fields
  do eq=1,n_eq

     ! Radiative transition matrix
     B = 0.0_dp
     A = 0.0_dp
     do iTrans = 1, nTrans_tot
        k = iTransUpper(iTrans)
        l = iTransLower(iTrans)

        if (eq==1) then
           JJmol = Jmol(iTrans,id)
        else
           JJmol = Jmol2(iTrans,id)
        endif

        !write(*,*) "Jmol", iTrans, JJmol, Bnu(transfreq(iTrans),350.)

        !JJmol = Bnu(transfreq(iTrans),Tdust) !  OK, redonne Tex = 350.

        A(k,k) = A(k,k) + Bul(iTrans) * JJmol  + Aul(iTrans)
        A(l,l) = A(l,l) + Blu(iTrans) * JJmol
        A(k,l) = A(k,l) - Blu(iTrans) * JJmol
        A(l,k) = A(l,k) - Bul(iTrans) * JJmol  - Aul(iTrans)
     enddo

     ! Test
     !A=0.0 ! Ok : I revert to LTE

     ! Matrice totale
     do i = 1, nLevels
        A(i,i) = A(i,i) + cTot(i)
        do j = 1, nLevels
           if (i /= j) then
              A(i,j) = A(i,j) - C(j, i)
           endif
        enddo
     enddo

     ! We replace the last equation (which is a LC of the others)
     ! by the conservation equation
     A(nLevels,:) = 1.0_dp ! the total sum of all populations
     B(nLevels) = 1.0_dp   ! = 1

     ! Resolution systeme matriciel par methode Gauss-Jordan
     call GaussSlv(A, B, nLevels)

     if (eq==1) then
        tab_nLevel(icell,:) = gas_density(icell) * tab_abundance(icell) * B(:)
     else
        tab_nLevel2(icell,:) = gas_density(icell) * tab_abundance(icell) * B(:)
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

  ! Equilibrium with the CMB for all cells
  do icell=1, n_cells
     ! The radiation field is equal to the CMB
     Jmol(:,id) = tab_Cmb_mol(:)

     ! Equilibrium
     call equilibre_rad_mol_loc(id,icell)
  enddo !icell

  return

end subroutine equilibre_othin_mol

!***********************************************************

subroutine equilibre_othin_mol_pop2()
  ! Computes the populations in the optically thin case
  ! equilibrium with the external radiation field and stored
  ! in the second population
  ! C. Pinte
  ! 15/10/07

  implicit none

  real(kind=dp) :: tab_nLevel_tmp(n_cells,nLevels)  ! not 3D
  logical :: ldouble_RT_tmp

  Jmol(:,:) = 0.0_dp

  ! For safety: safeguard population 1
  tab_nLevel_tmp(:,:) =  tab_nLevel(:,:)
  ldouble_RT_tmp = ldouble_RT
  ldouble_RT = .false.

  call equilibre_othin_mol()

  ! Initialisation of population 2
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

  real(kind=dp), dimension(ispeed(1):ispeed(2)) :: etau, P, opacity, Snu
  real(kind=dp) :: total_sum, J
  integer :: iTrans, iray

  Jmol(:,id) = 0.0_dp

  do iTrans=1, nTrans_tot
     total_sum = 0.0_dp
     J = 0.0_dp
     do iray=1, n_rayons
        P(:) =  Doppler_P_x_freq(:,iray,id)
        opacity(:) = kappa_mol_o_freq(icell,iTrans) * P(:) + kappa(icell,iTrans)
        etau(:) = exp(-ds(iray,id) * opacity(:)) ! exp(-tau)

        Snu(:) = ( emissivite_mol_o_freq(icell,iTrans) * P(:) + &
             emissivite_dust(icell,iTrans) ) / (opacity(:) + 1.0e-30_dp)

        J = J + sum( (I0(:,iTrans,iray,id) * etau(:) + Snu(:) * (1.0_dp - etau(:))) * P(:))
        total_sum = total_sum + sum(P(:))
     enddo ! iray
     Jmol(iTrans,id) =  J / total_sum

  enddo ! iTrans

  ! Normalization par le nombre de radiuss utilisesby the number of used radii
  !Jmol(:,id) = Jmol(:,id) / n_rayons

  if (ldouble_RT) then
     Jmol2(:,id) = 0.0_dp

     do iray=1, n_rayons
        do iTrans=1, nTrans_tot
           P(:) =  Doppler_P_x_freq(:,iray,id)
           opacity(:) = kappa_mol_o_freq2(icell,iTrans) * P(:) + kappa(icell,iTrans)
           etau(:) = exp(-ds(iray,id) * opacity(:)) ! exp(-tau)

           Snu(:) = ( emissivite_mol_o_freq2(icell,iTrans) * P(:) + emissivite_dust(icell,iTrans) ) &
                / (opacity(:) + 1.0e-30_dp)
           J = sum( (I0(:,iTrans,iray,id) * etau(:) + Snu(:) * (1.0_dp - etau(:))) &
                * P(:)) / sum(P(:))

           Jmol2(iTrans,id) = Jmol2(iTrans,id) + J
        enddo ! iTrans
     enddo ! iray

     ! Normalisation by the number of rays used
     Jmol2(:,id) = Jmol2(:,id) / n_rayons
  endif !ldouble_RT

  return

end subroutine J_mol_loc

!***********************************************************

function v_proj(icell,x,y,z,u,v,w) !
  ! projected velocity at 1 point of a cell
  ! C. Pinte
  ! 13/07/07

  implicit none

  real(kind=dp) :: v_proj
  integer, intent(in) :: icell
  real(kind=dp), intent(in) :: x,y,z,u,v,w

  real(kind=dp) :: velocity, vx, vy, vz, v_r, v_phi, v_theta, v_rcyl, norm, r, phi, rcyl, rcyl2, r2
  real(kind=dp) :: cos_phi, sin_phi, cos_theta, sin_theta

  if (lVoronoi) then
     vx = Voronoi(icell)%vxyz(1)
     vy = Voronoi(icell)%vxyz(2)
     vz = Voronoi(icell)%vxyz(3)

     v_proj = vx * u + vy * v + vz * w
  else
    !also work with pluto's models, model ascci and model 1d!
     if (lvelocity_file) then
        if (vfield_coord == 1) then ! cartesian velocity field
           vx = vfield3d(icell,1) ; vy = vfield3d(icell,2) ; vz = vfield3d(icell,3)
        else if (vfield_coord == 2) then
           ! Convert the velocity field from cylindrical to Cartesian coordinates
           v_r = vfield3d(icell,1) ; v_phi = vfield3d(icell,2) ;  vz = vfield3d(icell,3)
           phi = atan2(y, x)
           vx = cos(phi) * v_r - sin(phi) * v_phi
           vy = sin(phi) * v_r + cos(phi) * v_phi
           if ((.not.l3d).and.(z < 0)) vz = -vz
        else ! vfield == 3 --> spherical
           ! Convert the velocity field from spherical to Cartesian coordinates
           v_r = vfield3d(icell,1) ; v_phi = vfield3d(icell,2) ;  v_theta = vfield3d(icell,3)
           if ((.not.l3d).and.(z < 0)) v_theta = -v_theta

           rcyl2 = x*x + y*y
           r2 = rcyl2 + z*z
           rcyl = sqrt(rcyl2)
           r = sqrt(r2)

           cos_theta = rcyl/r
           sin_theta = z/r

           cos_phi = x/rcyl
           sin_phi = y/rcyl
           ! write(*,*) cos_phi, cos(atan2(y, x)) ! OK

           vz = -v_theta * cos_theta + v_r * sin_theta
           v_rcyl = v_theta * sin_theta + v_r * cos_theta

           vx = v_rcyl * cos_phi - v_phi * sin_phi
           vy = v_rcyl * sin_phi + v_phi * cos_phi
        endif

        v_proj = vx * u + vy * v + vz * w

     else ! Using analytical velocity field
        velocity = vfield(icell)

        if (lkeplerian) then
           r = sqrt(x*x+y*y)
           if (r > tiny_dp) then
              norm = 1.0_dp/r
              vx = -y * norm * velocity
              vy = x * norm * velocity
              vz = 0.
              if (linfall) then ! Adding extra velocity
                 r = sqrt(x*x+y*y+z*z) ; norm = 1.0_dp/r
                 vx = vx - chi_infall * x * norm * velocity
                 vy = vy - chi_infall * y * norm * velocity
                 vz = vz - chi_infall * z * norm * velocity
                 v_proj = vx * u + vy * v + vz * w
              else
                 v_proj = vx * u + vy * v
              endif
           else
              v_proj = 0.0_dp
           endif
        else if (linfall) then
           r = sqrt(x*x+y*y+z*z)
           !  if (lbenchmark_water2)  velocity = -r  * 1e5 * AU_to_pc    ! TMP for bench water2 : it changes nothing !!!????
           if (r > tiny_dp) then
              norm = 1.0_dp/r
              vx = x * norm * velocity
              vy = y * norm * velocity
              vz = z * norm * velocity
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

real(kind=dp) function dv_proj(icell,x0,y0,z0,x1,y1,z1,u,v,w)
  ! Projected velocity differential between 2 points
  ! au sein d'une cell
  ! C. Pinte
  ! 13/07/07

  implicit none

  integer, intent(in) :: icell
  real(kind=dp), intent(in) :: x0,y0,z0,x1,y1,z1,u,v,w

  real(kind=dp) :: velocity, vx0, vy0, vz0, vx1, vy1, vz1, norm

  velocity = vfield(icell)

  if (linfall) then
     ! Champ de velocity au point 0
     norm = 1.0_dp/sqrt(x0*x0+y0*y0+z0*z0)
     vx0 = x0 * norm * velocity
     vy0 = y0 * norm * velocity
     vz0 = z0 * norm * velocity

     ! Champ de velocity au point 1
     norm = 1.0_dp/sqrt(x1*x1+y1*y1+z1*z1)
     vx1 = x1 * norm * velocity
     vy1 = y1 * norm * velocity
     vz1 = z1 * norm * velocity

     dv_proj = (vx1 - vx0) * u + (vy1 - vy0) * v + (vz1 - vz0) * w
  else if (lkeplerian) then
      ! Champ de velocity au point 0
     norm = 1.0_dp/sqrt(x0*x0+y0*y0)
     vx0 = -y0 * norm * velocity
     vy0 = x0 * norm * velocity

     ! Champ de velocity au point 1
     norm = 1.0_dp/sqrt(x1*x1+y1*y1)
     vx1 = -y1 * norm * velocity
     vy1 = x1 * norm * velocity

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

  !$omp parallel default(none) &
  !$omp private(icell,CD,ldeplete) &
  !$omp shared(tab_abundance,freeze_out_depletion,T_freeze_out,lphoto_desorption,n_cells,Tdust)
  !$omp do
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
  !$omp end do
  !$omp end parallel

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

  write (*,*) "Photo-dissociating molecules ..."!, threshold_CD

  !$omp parallel default(none) &
  !$omp private(icell,CD) &
  !$omp shared(tab_abundance,n_cells, photodissociation_factor)
  !$omp do
  do icell=1, n_cells
     CD = compute_vertical_CD(icell)
     if (CD < threshold_CD * photodissociation_factor) then
        tab_abundance(icell) = tab_abundance(icell) * photo_dissocation_depletion
     endif
  enddo ! icell
  !$omp end do
  !$omp end parallel

  write(*,*) "Done"

  return

end subroutine photo_dissociation


!***********************************************************

subroutine write_abundance(imol)

  use fits_utils, only : cfitsWrite

  integer, intent(in) :: imol
  character(len=512) :: filename

  filename = trim(data_dir2(imol))//"/abundance.fits"
  if (.not.lVoronoi) then
     if (l3D) then
        call cfitsWrite(trim(filename),tab_abundance,[n_rad,2*nz,n_az])
     else
        call cfitsWrite(trim(filename),tab_abundance,[n_rad,nz])
     endif
  else
     call cfitsWrite(trim(filename),tab_abundance,[n_cells])
  endif

  return

end subroutine write_abundance

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
     CD = CD + (l_contrib * AU_to_m) * gas_density(icell0) ! part.m^-2
  enddo

  return

end function compute_vertical_CD

!***********************************************************

end module molecular_emission
