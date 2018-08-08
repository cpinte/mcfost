module benchmarks

  use parametres
  use mcfost_env
  use constantes
  use grains
  use opacity
  use em_th
  use molecular_emission
  use messages
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

subroutine read_Pascucci_cross_sections(lambda, Cext, Csca)
! Lecture table section efficace
! Benchmark Pascucci et al. 2004
! bande V : ligne 15

  use utils, only : in_dir

  implicit none

  integer, intent(in) :: lambda
  real, intent(out) :: Cext, Csca

  integer :: l, ios
  real, dimension(61), save :: Cext_Pascucci, Csca_Pascucci

  character(len=512) :: filename, dir

  if (lambda == 1) then
     filename = dust_pop(1)%indices(1)
     dir = in_dir(filename, dust_dir,  status=ios)
     if (ios /=0) call error("dust file cannot be found:"//trim(filename))
     filename = trim(dir)//trim(filename) ;

     open(unit=1,file=filename,status="old")
     ! Skipping header
     do l=1,11
        read(1,*)
     enddo

     do l=1,n_lambda
        read(1,*) tab_lambda(l), Csca_Pascucci(l), Cext_Pascucci(l)
     enddo
     close(unit=1)
  endif

  ! Convert to mum**2
  Cext = Cext_Pascucci(lambda) * (m_to_mum)**2
  Csca = Csca_Pascucci(lambda) * (m_to_mum)**2

  return

end subroutine read_Pascucci_cross_sections

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
