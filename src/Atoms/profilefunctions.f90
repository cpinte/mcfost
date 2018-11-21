MODULE ProfileFunctions

 use atom_type, only : AtomicLine, AtomType, ZeemanType
 use atmos_type, only : atmos, Hydrogen, Helium
 use voigtfunctions, only : Voigt
 use broad, only : Damping
 use zeeman, only : ZeemanMultiplet
 use constant

 IMPLICIT NONE

 ! temporaire, la grille de vitesse de MCFOST
 ! projetée sur une cellule sera utilisée ici
 integer, parameter :: NVEL=1

 double precision, dimension(NVEL) :: velocarr
 data velocarr /1000./ !m/s

 CONTAINS

 SUBROUTINE Profile(icell, atom, kr)
  ! There are significant differences with the RH code here
  ! because, in RH, angle_dep routines use short char -like
  ! solvers and the transfer equation is solved from top to bottom
  ! and bottom to top so that the profiles are saved for mu>0 and mu<0.
  ! In MC method we compute the (local) profile in each cell, which have
  ! a velocity field known. Further we do not need to save the profile
  ! has a function Nspace*Nrays, or rather we only compute the profile
  ! for the Nspace size of the cell, at all the relevent frequency points.
  type (AtomType), intent(inout) :: atom
  type (AtomicLine) :: line
  type (ZeemanType) :: zm
  integer, intent(in) :: kr, icell
  logical :: angle_dep = .true. !most of time, atmosphere
                                ! is moving or eventually magnetized
  double precision :: H, F, vbroad(atmos%Nspace), adamp, Larmor
  integer :: Nlamu, k, la

  line = atom%lines(kr)

  Nlamu = line%Nlambda

  ! Initialize the ratio between emission and absorption coefficient
  ! if PRD is included
  if (line%PRD) then
  allocate(line%rho_prd(Nlamu,atmos%Nspace))
  do k=1,atmos%Nspace
   do la=1,Nlamu
    line%rho_prd(la,k) = 1.
   end do
  end do
  end if

  if (.not.atmos%moving .or. .not.atmos%magnetized) angle_dep=.false.
  ! even in the (rare) case of angle independent (angle_dep=.false.)
  ! their will be no big differences in the line profile, mainly because
  ! short char and MC methods do not act in the same way.

  vbroad = atom%vbroad
  CALL Damping(icell, atom, kr, adamp)

  !allocate(line%wphi(atmos%Nspace))

  if (line%polarizable) then
   Larmor = (Q_ELECTRON / (4.*PI*M_ELECTRON)) * (line%lambda0*NM_TO_M)
   ! fill Zeeman structure
   CALL ZeemanMultiplet(zm) !empty
   ! write(*,*) "Number of Zeeman component for line j->i
  end if

 RETURN
 END SUBROUTINE Profile

END MODULE ProfileFunctions
