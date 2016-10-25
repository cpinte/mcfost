module benchmarks

  use parametres
  use constantes
  use grains
  use opacity
  use em_th
  use molecular_emission

  implicit none

  contains

subroutine lect_section_eff()
! Lecture table section efficace
! Benchmark Pascucci et al. 2004
! bande V : ligne 15

  implicit none

  integer :: lambda, l, icell
  real, dimension(61) :: qext, qsca
  real, dimension(0:61) :: tab_lambda_lim

  open(unit=1,file="optSi.dat",status="old")

  if (n_lambda /= 61) then
     write(*,*) "Benchmark error"
     write(*,*) "n_lambda must be 61"
     stop
  endif

  if (n_grains_tot /= 1) then
     write(*,*) "Benchmark error"
     write(*,*) "n_grains_tot must be 1"
     stop
  endif

  do lambda=1,n_lambda
     read(1,*) tab_lambda(lambda), qsca(lambda), qext(lambda)
  enddo

  close(unit=1)

  ! Propriétés optiques des cellules
  do lambda=1,n_lambda
     tab_albedo_pos(:,lambda)=qsca(lambda)/qext(lambda)
     tab_s11_pos(:,:,lambda)=1.0
     tab_s12_o_s11_pos(:,:,lambda)=0.0
     tab_s33_o_s11_pos(:,:,lambda)=0.0
     tab_s34_o_s11_pos(:,:,lambda)=0.0
     do icell=1, n_cells
        ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
        ! a² m² -> 1e4 cm²                    \    /\ Cas particulier benchmark !!
        ! density en cm-3                      > reste facteur 1.49595e17
        ! longueur de vol en AU = 1.5e13 cm   /
        kappa(icell,lambda)=densite_pouss(1,icell) * qext(lambda) * 1.49595e17
        kappa_abs_LTE(icell,lambda)=densite_pouss(1,icell)*(qext(lambda)-qsca(lambda))*1.49595e17
     enddo !icell
  enddo !lambda

  tab_g_pos=0.0
  proba_abs_RE = 1.0


  ! Largeur des bins en longueur d'onde
  ! tab_lambda_lim : limite sup du bin
  do l=1,n_lambda-1
     tab_lambda_lim(l)=sqrt(tab_lambda(l)*tab_lambda(l+1))
  enddo

  tab_lambda_lim(0)=tab_lambda(1)/(tab_lambda_lim(1)/tab_lambda(1))
  tab_lambda_lim(n_lambda)=tab_lambda(n_lambda)*(tab_lambda(n_lambda)/tab_lambda_lim(n_lambda-1))

  do l=1,n_lambda
     tab_delta_lambda(l)=tab_lambda_lim(l)-tab_lambda_lim(l-1)
  enddo

  return

end subroutine lect_section_eff

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

end module benchmarks
