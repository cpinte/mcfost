! Routines for line broadening, including van der Waals and
! Stark broadening by electronic collisions.
!
! Adapted from RH H. Uitenbroek


MODULE broad

 use constant
 use atom_type
 use atmos_type, only : atmos, Hydrogen, Helium
 use math, only : dpow, SQ, CUBE

 IMPLICIT NONE

 real(8), parameter :: AVERAGE_ATOMIC_WEIGHT = 28.0

 CONTAINS


 SUBROUTINE VanderWaals(icell, atom, kr, GvdW)
 ! Compute van der Waals broadening in Lindholm theory with
 !Unsold's approximation for the interaction coefficient C_6.
 !
 ! See: Traving 1960, "Uber die Theorie der Druckverbreiterung
 ! von Spektrallinien", p 91-97
 !
 !  -- Mihalas 1978, p. 282ff, and Table 9-1.
 !
 ! --->
 !  -- Hubeny & Mihalas 2014, chap 8, sect. 8.3
 ! and table 8.1
 ! --->
 !
 !Gamma = 8.08 * vrel^3/5 * C_6^2/5 * atmos.H.ntotal_neutral
 !
 !
 ! Or with parametrization using Smirnov-Roueff potential.
 !
 ! See: DeRidder & van Rensbergen 1976, A&A Suppl. 23, 147-165
 !
 ! Gamma = alpha * T^{beta} * atmos.H.ntotal_neutral
 !
 !
 ! Or with the Anstee, Barklem & O'Mara formalism.
 !
 ! See: Anstee & O'Mara 1995, MNRAS 276, 859-866
 ! Barklem & O'Mara 1998, MNRAS 300, 863-871
 !
 ! Gamma =
 ! (4/pi)^(alpha/2) * &
 !     Gam((4-alpha)/2) v_0 sigma(v_0)(vmean/v_0)^(1-alpha)
  double precision, intent(out) :: GvdW
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (AtomicLine) :: line
  type (Element) :: Helium
  integer, intent(in) :: kr
  integer :: i, j, ic, Z, k
  double precision :: vrel35_H, vrel35_He, fourPIeps0, deltaR
  double precision :: cross, gammaCorrH, gammaCorrHe, C625

  Helium = atmos%Elements(2) ! Helium basic informations
  line = atom%lines(kr)

  ! here, starts a 1 because it is j saved in line%j, which
  ! is corrected by +1 wrt the j read in the atomic file
  j = line%j
  i = line%i

  if (line%vdWaals.eq."UNSOLD" .or. line%vdWaals.eq."BARKLEM") then
   fourPIeps0 = 4.*PI*EPSILON_0
   vrel35_He = dpow(8.*KBOLTZMANN/(PI*AMU*atom%weight) * &
         (1.+atom%weight/Helium%weight),3d-1) !0.3

   Z = atom%stage(j)+1 !remember, stage starts at 0 like in C
   !--> voir barklem.f90, Z is not an index but a physical value
   ! independent of the language (??)

   ic = j+1 !ic is an index, the law of +1 wrt C index applies
   do while (atom%stage(ic).lt.atom%stage(j)+1)
    ic = ic+1
   end do
   !write(*,*) Z, ic
   deltaR = SQ(E_RYDBERG/(atom%E(ic)-atom%E(j))) - &
            SQ(E_RYDBERG/(atom%E(ic)-atom%E(i)))
   C625 = dpow(2.5*(SQ(Q_ELECTRON)/fourPIeps0)*&
           (ABARH/fourPIeps0) * 2.*PI*&
           SQ(Z*RBOHR)/HPLANCK * deltaR,4d-1) !0.4
  end if

  SELECT CASE (line%vdWaals)
  CASE ("UNSOLD")
   ! Relative velocity of radiator and perturber with Maxwellian
   ! velocity distributions
   vrel35_H = dpow(8.*KBOLTZMANN/(PI*AMU*atom%weight) *&
         (1.+atom%weight/Hydrogen%weight),3d-1)
   cross = 8.08 * (line%cvdWaals(1)*vrel35_H+line%cvdWaals(3)*&
           Helium%abund*vrel35_He)*C625
   !write(*,*) Helium%abund, C625
    GvdW = cross * dpow(atmos%T(icell), 3d-1)
    !write(*,*) GvdW(k)
  CASE ("BARKLEM")
   write(*,*) "Warning-> vdWaals broadening following BARKLEM ",&
     "not tested yet!"
   ! UNSOLD for Helium
   cross = 8.08 * line%cvdWaals(3)*Helium%abund*vrel35_He*C625
    GvdW = line%cvdWaals(1) * dpow(atmos%T(icell),real(1.-&
          line%cvdWaals(2),kind=8)/2.) + cross*dpow(atmos%T(icell),3d-1)
  CASE ("RIDDER_RENSBERGEN")
   write(*,*) "Warning-> vdWaals broadening following ",&
     "RIDDER_RENSBERGEN not tested yet!"
  !alpha = 1.0E-8 * cvdW[0]  (Hydrogen broadening)
  !      = 1.0E-9 * cvdW[1]  (Helium broadening)
   gammaCorrH = 1d-8 * CUBE(CM_TO_M) * &
     dpow(1.+Hydrogen%weight/atom%weight,&
         real(line%cvdWaals(2),kind=8))
   gammaCorrHe = 1d-9 * CUBE(CM_TO_M) * &
     dpow(1.+Helium%weight/atom%weight, &
         real(line%cvdWaals(4),kind=8))
    GvdW = gammaCorrH*line%cvdWaals(1) * &
      dpow(atmos%T(icell),real(line%cvdWaals(2),kind=8)) + &
      gammaCorrHe*line%cvdWaals(3) * &
      dpow(atmos%T(icell),real(line%cvdWaals(4),kind=8)) *&
           Helium%abund
  CASE DEFAULT
   write(*,*) "Method for van der Waals broadening unknown"
   write(*,*) "exiting..."
   stop
  END SELECT

  ! Multiply with the Hydrogen ground level population
   !write(*,*) "GvdW=", GvdW(k)
   GvdW = GvdW * Hydrogen%n(1,icell)

 RETURN
 END SUBROUTINE VanderWaals

 SUBROUTINE Stark(icell, atom, kr, GStark)
 ! Quadratic Stark broadening by electrons and singly charged
 ! ions.
 !
 ! Gamma = 11.37 * vrel^1/3 * C_4^2/3 * (ne + nion)
 !
 ! Use estimate for C_4 from Traving.
 !
 ! See: Traving 1960, "Uber die Theorie der Druckverbreiterung
 !      von Spektrallinien", p 93
 !
 ! -- Mihalas 1978, p. 282ff, and Table 9-1.
 !
 ! -- David F. Gray, Observation and Analysis of Stellar
 !    Photospheres (1992), 2nd ed., p. 216, eq. 11.33
 !
 !
 ! if (line%cStark .lt.0) then Gamma = abs(line%cStark) * ne
 !
  double precision, intent(out) :: GStark
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (AtomicLine) :: line
  integer, intent(in) :: kr
  integer :: k, ic, Z
  double precision :: C4, C, Cm, melONamu, cStark23, cStark
  double precision :: neff_u, neff_l, vrel, ERYDBERG

  melONamu = M_ELECTRON / AMU
  line = atom%lines(kr)

  if (line%cStark .lt. 0.) then
   cStark = dabs(real(line%cStark,kind=8))
    GStark = cStark * atmos%ne(icell)
  else
   ! Constants for relative velocity. Assume that nion = ne
   ! and that the average atomic weight of ionic pertubers is
   ! given by AVERAGE_ATOMIC_WEIGHT

   C = 8. * KBOLTZMANN / (PI*AMU*atom%weight)
   Cm = dpow(1.+atom%weight/melONamu, 1.6666667d-1) + &
      dpow(1.+atom%weight/AVERAGE_ATOMIC_WEIGHT,1.6666667d-1)

   ! Find core charge Z and effective quantumer numbers neff_u
   ! and neff_l for upper and lower level
   Z = atom%stage(line%i) + 1 !physical, not index
   ic = line%i + 1
   do while ((atom%stage(ic).lt.atom%stage(line%i)+1).and.&
        ic.lt.atom%Nlevel)
    ic = ic + 1
   end do
   
   if (atom%stage(ic).eq.atom%stage(line%i)) then
    write(*,*) "(Broad.Stark) Cannot find overlying continuum for level ",&
      line%i," of atom ",atom%ID," exiting..."
    write(*,*) atom%stage(ic), atom%stage(line%i)
    stop
   end if
   ERYDBERG = E_RYDBERG / (1.+M_ELECTRON / (atom%weight*AMU))
   neff_l = Z*dsqrt(ERYDBERG/(atom%E(ic)-atom%E(line%i)))
   neff_u = Z*dsqrt(ERYDBERG/(atom%E(ic)-atom%E(line%j)))
   C4 = (SQ(Q_ELECTRON)/(4.*PI*EPSILON_0))*RBOHR *   &
     (2.*PI*SQ(RBOHR)/HPLANCK) / (18.*Z*Z*Z*Z) * &
     (SQ(neff_u*(5.0*SQ(neff_u) + 1.0)) - SQ(neff_l* &
     (5.0*SQ(neff_l) + 1.0)))
   cStark23 = 11.37 * dpow(line%cStark * C4,6.6666667d-1)

    vrel = dpow(C*atmos%T(icell),1.6666667d-1) * Cm
    GStark = cStark23 * vrel * atmos%ne(icell)
  end if

 RETURN
 END SUBROUTINE Stark

 SUBROUTINE StarkLinear(icell, atom, kr, GStark)
 ! Linear Stark broadening by electrons for hydrogen lines.
 !
 ! See: K. Sutton (1978), JQSRT 20, 333-343
 !
 ! GStark = a_1 * &
 !       [0.60 * (n_u^2 - n_l^2) * (N_e)^(2/3) * CM_TO_M^2]
 !
  double precision, intent(out) :: GStark
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (AtomicLine) :: line
  integer, intent(in) :: kr
  double precision :: n_upper, n_lower
  integer :: k
  double precision :: a1, C

  line = atom%lines(kr)

  ! or n_lower = dsqrt(atom%g(i)/2.) !with g the statistical weight
  n_lower = dsqrt(atom%g(line%i)/2.)
  n_upper = dsqrt(atom%g(line%j)/2.) !works because stark only for Hydrogen.
!   if (.not.getPrincipal(atom%label(line%i),n_lower)) then
!    write(*,*) "Cannot find principal quantum number for label ", &
!     atom%label(line%i), " exiting..."
!    stop
!   end if
!   if (.not.getPrincipal(atom%label(line%j),n_upper)) then
!    write(*,*) "Cannot find principal quantum number for label ", &
!    atom%label(line%j), " exiting..."
!    stop
!   end if

  if ((n_upper-n_lower).eq.1) then
   a1 = 0.642
  else
   a1 = 1.
  end if
  C = a1 * 0.6 * (n_upper*n_upper-n_lower*n_lower) * SQ(CM_TO_M)

   GStark = C*dpow(atmos%ne(icell), 6.6666667d-1)

 RETURN
 END SUBROUTINE StarkLinear

 SUBROUTINE Damping(icell, atom, kr, adamp)
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (AtomicLine) :: line
  integer, intent(in) :: kr
  integer :: k
  double precision :: cDop
  double precision, intent(out) :: adamp
  double precision :: Qelast

  Qelast = 0.0

  line=atom%lines(kr)
  !write(*,*) "Computing damping for line ",line%j,"-",line%i

  cDop = (NM_TO_M*line%lambda0) / (4.*PI)
  ! van der Waals broadening
  if ((line%cvdWaals(1).gt.0).or.(line%cvdWaals(3).gt.0)) then
   CALL VanderWaals(icell, atom, kr, adamp)
    Qelast = Qelast + adamp
    !write(*,*) "vdW(k)=",adamp(k), Qelast(k)
  end if
  ! Quadratic Stark broadening
  if (line%cStark.ne.0.) then
   CALL Stark(icell, atom, kr, adamp)
    Qelast = Qelast + adamp
    !write(*,*) "QStark(k)=",adamp(k), Qelast(k)
  end if
  ! Linear Stark broadening only for Hydrogen
  if (atom%ID.eq."H ") then
   CALL StarkLinear(icell, atom,kr, adamp)
    Qelast = Qelast + adamp
    !write(*,*) "LStark(k)=",adamp(k), Qelast(k)
  end if
  ! Store Qelast, total rate of elastic collisions, if PFR only
  if (line%PFR) then
   write(*,*) "PFR not handled yet, avoiding"
   ! line%Qelast = Qelast
  end if

   adamp = (line%Grad + Qelast)*cDop / atom%vbroad(icell)
   !write(*,*) adamp(k), Qelast(k), cDop, atom%vbroad(k)

 RETURN
 END SUBROUTINE Damping

END MODULE broad
