module benchmarks

  use parametres
  use mcfost_env
  use constantes
  use grains
  use messages
  use molecular_emission
  use wavelengths

  implicit none

  contains

subroutine init_Pascucci_benchmark()

  character(len=32) :: filename = "Pascucci_optSi.dat"

  write(*,*) "------------------------------------------------"
  write(*,*) "| Setting up the Pascucci et al 2004 benchmark |"
  write(*,*) "------------------------------------------------"

  if (n_grains_tot /= 1) call error("Benchmark : n_grains_tot must be 1")

  write(*,*) "! Forcing dust grain density to 3.6 g.cm-3      |"
  dust_pop(1)%component_rho1g(1) = 3.6
  dust_pop(1)%rho1g_avg = 3.6

  write(*,*) "! Forcing isotropic scattering                  |"
  write(*,*) "-------------------------------------------------"
  lisotropic=.true.

  return

end subroutine init_Pascucci_benchmark

!*************************************************************

!-- subroutine read_Pascucci_cross_sections(lambda, Cext, Csca)
!-- ! Lecture table section efficace
!-- ! Benchmark Pascucci et al. 2004
!-- ! bande V : ligne 15
!--
!--   use utils, only : in_dir
!--
!--   implicit none
!--
!--   integer, intent(in) :: lambda
!--   real, intent(out) :: Cext, Csca
!--
!--   integer :: l, ios
!--   real, dimension(61), save :: Cext_Pascucci, Csca_Pascucci
!--
!--   character(len=512) :: filename, dir
!--
!--   if (lambda == 1) then
!--      filename = dust_pop(1)%indices(1)
!--      dir = in_dir(filename, dust_dir,  status=ios)
!--      if (ios /=0) call error("dust file cannot be found:"//trim(filename))
!--      filename = trim(dir)//trim(filename) ;
!--
!--      open(unit=1,file=filename,status="old")
!--      ! Skipping header
!--      do l=1,11
!--         read(1,*)
!--      enddo
!--
!--      do l=1,n_lambda
!--         read(1,*) tab_lambda(l), Csca_Pascucci(l), Cext_Pascucci(l)
!--      enddo
!--      close(unit=1)
!--   endif
!--
!--   ! Convert to mum**2
!--   Cext = Cext_Pascucci(lambda) * (m_to_mum)**2
!--   Csca = Csca_Pascucci(lambda) * (m_to_mum)**2
!--
!--   return
!--
!-- end subroutine read_Pascucci_cross_sections

!*************************************************************

subroutine readMolecule_benchmark1()

   implicit none

  character(len=80) :: junk
  real(kind=dp) :: Delta_E, Kul

  mol(1)%filename = "mol_benchmark.dat"

  open(unit=1, file=mol(1)%filename, status="old")

  read(1,*) junk

  nLevels = 2
  nTrans_tot  = 1

  allocate(Level_energy(nLevels),poids_stat_g(nLevels),j_qnb(nLevels))
  allocate(Aul(1:nTrans_tot),fAul(nTrans_tot))
  allocate(Bul(1:nTrans_tot),fBul(nTrans_tot))
  allocate(Blu(1:nTrans_tot),fBlu(nTrans_tot))
  allocate(transfreq(1:nTrans_tot))
  allocate(itransUpper(1:nTrans_tot))
  allocate(itransLower(1:nTrans_tot))

  poids_stat_g(1) = 1.0_dp

  read(1,*) Delta_E, poids_stat_g(2), Aul(1), Kul

  Level_energy(1) = 0.0_dp
  Level_energy(2) = Delta_E  / 8065.541  ! per cm to ev


  transfreq(1) = c_light * (Delta_E * 100_dp) ! en Hz
  itransUpper(1) = 2
  itransLower(1) = 1

  ! Transformation Aul -> Bul
  Bul(1) = Aul(1) * (c_light**2)/(2.d0*hp*(transfreq(1))**3)
  ! Transformation Bul -> Blu
  Blu(1) = Bul(1) * poids_stat_g(2)/poids_stat_g(1)

  fAul(:) = Aul(:) * hp * transfreq(:)/(4*pi)
  fBul(:) = Bul(:) * hp * transfreq(:)/(4*pi)
  fBlu(:) = Blu(:) * hp * transfreq(:)/(4*pi)


  nCollPart = 1

  allocate(nCollTrans(1:nCollPart))
  allocate(nCollTemps(1:nCollPart))
  allocate(collTemps(1:nCollPart, 1:20))
  allocate(collBetween(1:nCollPart))

  nCollTrans(1) = 1
  nCollTemps(1) = 1

  collTemps = 1.0_dp

  allocate(collRates(1:nCollPart, 1:nCollTrans(1), 1:nCollTemps(1))) ! TODO : passage par des pointeurs, c'est crade
  allocate(iCollUpper(1:nCollPart, 1:nCollTrans(1)))
  allocate(iCollLower(1:nCollPart, 1:nCollTrans(1)))
  collRates(1,1,1) = Kul
  iCollUpper(1,1) = 2
  iCollLower(1,1) = 1

  close(unit=1)

  write(*,*) "Molecular file read successfully : Benchmark 1"

  return

end subroutine readMolecule_benchmark1

!*************************************************************

subroutine readMolecule_benchmark2()
  ! Lecture du fichier moleculaire pour le benchmark 2
  ! de van Zadelhoff 2002
  ! adapte de Torus
  ! C. Pinte
  ! 26/06/07

  implicit none

  real(kind=dp), parameter :: ergToEv = 6.24145d11
  real(kind=dp), parameter :: hCgs = 6.626205d-27

  character(len=80) :: junk
  integer :: i, iPart

  real :: molecularWeight

  open(1, file=mol(1)%filename, status="old", form="formatted")

  read(1,'(a)') mol(1)%name

  read(1,*) molecularWeight
  masse_mol = masseH * molecularWeight

  read(1,*) nLevels, nTrans_tot

  allocate(Level_energy(1:nLevels))
  allocate(poids_stat_g(1:nLevels))
  allocate(j_qnb(1:nLevels))

  allocate(Aul(nTrans_tot),fAul(nTrans_tot))
  allocate(Bul(nTrans_tot),fBul(nTrans_tot))
  allocate(Blu(nTrans_tot),fBlu(nTrans_tot))
  allocate(transfreq(nTrans_tot))
  allocate(itransUpper(nTrans_tot))
  allocate(itransLower(nTrans_tot))

  read(1,'(7f11.7)') Level_energy(1:nLevels)

  Level_energy = Level_energy / 8065.541  ! per cm to ev

  read(1,*) poids_stat_g(1:nLevels)
!  j_qnb(1:nLevels) = (poids_stat_g(1:nLevels) - 1.)/2.

  read(1,*) itransUpper(1:nTrans_tot)
  read(1,*) itransLower(1:nTrans_tot)

  read(1,*) Aul(1:nTrans_tot)

  do i = 1, nTrans_tot
     transFreq(i) = ((Level_energy(iTransUpper(i)) - Level_energy(iTransLower(i)))/ergToEV)/hcgs

     ! Transformation Aul -> Bul
     Bul(i) = Aul(i) * (c_light**2)/(2.d0*hp*(transfreq(i))**3)
     ! Transformation Bul -> Blu
     Blu(i) = Bul(i) * poids_stat_g(iTransUpper(i))/poids_stat_g(iTransLower(i))
  enddo

  fAul(:) = Aul(:) * hp * transfreq(:)/(4*pi)
  fBul(:) = Bul(:) * hp * transfreq(:)/(4*pi)
  fBlu(:) = Blu(:) * hp * transfreq(:)/(4*pi)

  read(1,*) junk

  nCollPart = 1 ! only one collision partner considered in benchmark
  iPart = 1 ! follows from ^

  allocate(nCollTrans(1:nCollPart))
  allocate(nCollTemps(1:nCollPart))

  allocate(collTemps(1:nCollPart, 1:4))

  read(1,*) nCollTrans(1), nCollTemps(1), collTemps(1,1:4)

  allocate(iCollUpper(1:nCollPart, 1:nCollTrans(iPart)))
  allocate(iCollLower(1:nCollPart, 1:nCollTrans(iPart)))

  allocate(collRates(1:nCollPart, &
       1:nCollTrans(iPart), 1:nCollTemps(ipart)))

  read(1,*) iCollUpper(1, 1:nCollTrans(1))
  read(1,*) iCollLower(1, 1:nCollTrans(1))
  do i = 1, nCollTemps(1)
     read(1,*) collRates(iPart,  1:nCollTrans(ipart), i)
  enddo

  close(1)

  return

end subroutine readMolecule_benchmark2

!*************************************************************

subroutine init_GG_Tau_mol()

  implicit none

  integer :: icell

  ldust_mol = .true.

  do icell=1,n_cells
     Tdust(icell) = 30.0 * (r_grid(icell)/100.)**(-0.5)
     Tcin(icell) = 30.0 * (r_grid(icell)/100.)**(-0.5)
  enddo

  icell = icell_ref
  write(*,*) "Density @ 100 AU", real(densite_gaz(icell) / 100.**3 * (sqrt(r_grid(icell)**2 + z_grid(icell)) / 100.)**2.75)

  return

end subroutine init_GG_Tau_mol

!***********************************************************

subroutine init_HH_30_mol()

  implicit none

  integer :: icell

  ldust_mol = .true.

  do icell=1,n_cells
     Tdust(icell) = 12.0 * (r_grid(icell)/100.)**(-0.55)
     vfield(icell) = 2.0 * (r_grid(icell)/100.)**(-0.55)
  enddo

  v_turb = 230.

  return

end subroutine init_HH_30_mol

!***********************************************************

subroutine init_benchmark_vanZadelhoff1()

  implicit none

  ldust_mol = .false.

  v_turb = 150._dp !m.s-1
  Tdust = 20.
  Tcin = 20._dp
  vfield(:) = 0.0

  masse_mol = 1.0

  linfall = .true.
  lkeplerian = .false.

  !tab_abundance = abundance

  write(*,*) "Density", real(densite_gaz(icell_ref) * &
       (r_grid(icell_ref)**2 + z_grid(icell_ref)**2)/rmin**2)

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
        Tdust(icell) =  tmp_T(l-1) + frac * (tmp_T(l) - tmp_T(l-1))
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

  !tab_abundance = abundance

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

  !tab_abundance = abundance

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

  !tab_abundance = abundance

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
           Tdust(icell) = tmp_T(1)
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
           Tdust(icell) = exp( log_tmp_T(l-1) + frac * (log_tmp_T(l) - log_tmp_T(l-1)) )
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

  !tab_abundance = abundance

  return

end subroutine init_benchmark_water3

!***********************************************************

end module benchmarks
