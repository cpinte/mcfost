PROGRAM test

use uplow
use grid_type
use atom_type
use spectrum_type
use constant
use math
use readatom
use getline
use barklem
use solvene
use collision
use lte
use voigtfunctions, only : VoigtArmstrong, VoigtHumlicek, Voigt
use metal
use hydrogen_opacities
use rayleigh_scattering
use thomson_scattering
use broad
use profilefunctions, only : Profile
use planck
implicit none

integer :: Nspace, k, Ncell, icell, nat, nact, n
logical :: re_init
type (ActiveSetType) :: as
double precision, allocatable, dimension(:) :: T, ner, nHtotr, vturbr, velfield, &
  magfield, chil, etal
double precision, allocatable, dimension(:) :: scatt, chi, eta, tmp1, tmp2, tmp3, Bpnu
double precision :: Bchar, Vchar
double precision, dimension(2) :: x1, x2, y
double precision, dimension(2) :: x3, x4
double precision, dimension(2,2) :: yy

Vchar = 1d0
Bchar = 0d0
!
!  x1(1)=-5
!  x2(1)=1
!  x1(2)=-5
!  x2(2)=1
! x3(1)=-4.5
! x4(1)=-3.1
! ! x3(2)=-4.5
! ! x4(2)=-4.5
!  yy(1,2)=0.2
!  yy(2,1)=1.1
!  yy(2,2)=1.
!  yy(1,1)=0.8
! ! !test interp2D
!  write(*,*) interp2D(x1,x2, yy, x3(1), x4(1)), interp2Darr(2, x1,2, x2,yy,1, x3,1, x4)
!  stop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!done in reading the model
!Do not forget to flat all quantities if needed


Nspace = 2 !flattened size of arrays
Ncell = Nspace

allocate(T(Nspace))
allocate(ner(Nspace))
allocate(nHtotr(Nspace))
allocate(vturbr(Nspace))
allocate(NLTEspec%I(1))
allocate(NLTEspec%Flux(1))

do k=1,Nspace
 T(k) = 5000.
 ner(k) = 0.
 nHtotr(k) = 0.
 vturbr(k) = 0.
end do
T(1) = 4500. !dep 28 from bottom in FAL_C 82-28
T(2) = 9400. !depth 82 FAL_C from comparison
ner(2)=3.831726E+15 !/cm3
ner(1)=2.441017E+11
nHtotr(2)=1.2887d+17+1.7545d+12 +3.8349d+11 +3.0146d+11+ 3.2285d+11 +3.7897d+15
nHtotr(1)= 2.0529d15+3.1103d04 + 7.1209d02 +2.0204d+02 +1.4078d+02 +4.4674d+09!/cm3
vturbr(1)=8.043894d-01 !km/s
vturbr(2)=1.806787d+00
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


CALL initGrid(Nspace, T, ner, nHtotr, vturbr, Vchar)

!actually electron density for all the grid
if (atmos%calc_ne) CALL SolveElectronDensity(atmos%ne)

!reading atomic models and allocate variable for all the grid: n, nstar and Rij, Rji
CALL readAtomicModels(1) !1 is a unit for reading freed out of the subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! init wavelength grid
CALL make_wavelength_grid(atmos%Natom, atmos%Atoms,&
     NLTEspec%lambda)
NLTEspec%Nwaves = size(NLTEspec%lambda)
NLTEspec%Nact = atmos%Nactiveatoms
re_init = .false. !first allocation

CALL initAS(re_init)
NLTEspec%atmos => atmos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!Setting LTEpops for the all grid and allocate space for C(Nl*Nl) for active atoms
! and repopen atomic file for active atoms for reading collision
Call setLTEcoefficients ()

!initiate NLTE pops
! CALL initialSol()
write(*,*) atmos%Nactiveatoms, " active atoms"
do nact=1,atmos%Nactiveatoms
  write(*,*) "Setting initial solution for active atom ", atmos%ActiveAtoms(nact)%ID, &
      atmos%ActiveAtoms(nact)%active
 atmos%ActiveAtoms(nact)%n = 0.97 * atmos%ActiveAtoms(nact)%nstar
end do

open(unit=12, file="../../test_suite/test_NLTE/opacities2.dat",status='unknown')

do icell=1, atmos%Nspace

! for all passive atoms and all wavelength at this cell point
CALL Background(icell)

!!!For testing!!
do n=1,NLTEspec%Nwaves
write(12,'(4E)') NLTEspec%lambda(n), NLTEspec%ActiveSet%sca_c(n), &
  NLTEspec%ActiveSet%chi_c(n), NLTEspec%ActiveSet%eta_c(n)
end do
!!!!!!!!


!read collision data and fill collisional matrix, initizlize to C(:)=0 at each cell
! the matrix is allocated for ACTIVE atoms only in setLTEcoefficients and the file open
! before the transfer starts, and closed at the end
do nact=1,atmos%Nactiveatoms   !ActiveAtoms(nact) is an alias for Atoms(n) if Atoms(n)%active
                        !so operation on ActiveAtoms should fill also Atoms(n)
if (atmos%ActiveAtoms(nact)%active) then
 CALL CollisionRate(icell, atmos%ActiveAtoms(nact))
end if
end do

CALL initAS(re_init=.true.)  !when re_init=.true. only set to 0 the different opacities
                         ! this prevents adding extra opac from cell to cell
end do !loop over cells
close(12)


!end transfer, freeing atoms, and model
! close atomic file kept open for collision data
do nact=1,atmos%Nactiveatoms
 CALL closeCollisionFile(atmos%ActiveAtoms(nact))
end do
CALL freeAS() !completely deallocate, unlike initAs(re_init=.true.)
CALL freeSpectrum()
CALL freeGrid()

END PROGRAM test
