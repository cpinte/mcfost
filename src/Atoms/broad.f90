
MODULE broad

 use constant
 use atom_type
 use atmos_type, only : atmos, Hydrogen, VBROAD_atom !, Helium !alias to ptr_atom of Helium, which might not exist
 
 use messages, only : error, warning

 IMPLICIT NONE

 real(8), parameter :: AVERAGE_ATOMIC_WEIGHT = 28.0

 CONTAINS

 SUBROUTINE Broad_Kurosawa(icell, atom, kr, adamp)
 ! For testing, implementing the broadening for H alpha used by some people
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: kr, icell
  double precision, intent(out) :: adamp
  double precision :: Crad, Cvdw, Cstark, G
  
  !AA
  Crad = 6.5d-4; CvdW = 4.4d-4; Cstark = 1.17d-3
  !in AA
  G = Crad + Cvdw * (1d-6 * sum(atom%n(1:atom%Nlevel-1,icell))/1d16) * (atmos%T(icell)/5d3)**(0.3) + &
      Cstark * (atmos%ne(icell)/1d12)**(2./3.) * 1d-4 !conversion from 12cm^-3 ^ 2/3 to 12m^-3 ^ 2/3
  G = G * 1d-10 !in m
  !conversion in s^-1 -> G(m) = G(s^-1) * lambda**2 / c : G(s^-1)/nu = G(m)/lambda -> G(m)*nu/lambda
  adamp  = G * CLIGHT / atom%lines(kr)%lambda0**2 * 1d18!s^-1
  !will be divided by 4PI later
 RETURN
 END SUBROUTINE Broad_Kurosawa

 !Building
 SUBROUTINE RadiativeDamping(icell, atom, kr, Grad)
 ! -------------------------------------------------------------------------- !
  ! Radiative damping of a line transition j->i,
  ! assuming the radiation field is given by a Planck function, and psi=phi
  !
  ! Grad = Sum_l<m Aml + Bml * int2(psiml*I*dv) + Sum_n>m mn int2(phimn*Idv)
 ! -------------------------------------------------------------------------- !
  integer, intent(in) :: icell, kr
  type (AtomType), intent(in) :: atom
  real(kind=dp), intent(out) :: Grad
  integer :: kp,l,n
  real(kind=dp) :: hc_lakT, gamma_j, gamma_i
  type (AtomicLine) :: other_line, line
  !Remove radiation field to compare with RH
  
  hc_lakT = hc_k / atmos%T(icell) !nm factor
  line = atom%lines(kr)
  Grad = 0d0
  gamma_j = 0d0; gamma_i = 0d0
  
  do kp=1,atom%Nline
    other_line = atom%lines(kp)
    l = other_line%i; n = other_line%j
    if (n==line%j) then !our upper level is also the upper level of another transition
     gamma_j = gamma_j + other_line%Aji/(1d0 - dexp(-hc_lakT/other_line%lambda0))
    elseif (l==line%j) then !our upper level is a lower level of another transition
     gamma_j = gamma_j + other_line%Bji/other_line%Bij * other_line%Aji / (dexp(hc_lakT/other_line%lambda0)-1d0)
!     elseif (l==line%i) then !our lower level is zn upper level
!      gamma_i = gamma_i + ohter
    endif
  
  enddo
  
  Grad = gamma_j + gamma_i

 RETURN
 END SUBROUTINE RadiativeDamping

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
  type (Element), pointer :: Helium
  integer, intent(in) :: kr
  integer :: i, j, ic, Z, k
  double precision :: vrel35_H, vrel35_He, fourPIeps0, deltaR
  double precision :: cross, gammaCorrH, gammaCorrHe, C625

  Helium => atmos%Elements(2)%ptr_elem ! Helium basic informations
  line = atom%lines(kr)

  ! here, starts a 1 because it is j saved in line%j, which
  ! is corrected by +1 wrt the j read in the atomic file
  j = line%j
  i = line%i

  !could be nice to keep vrel35 for each atom ? as it depends only on the weight

  fourPIeps0 = 4.*PI*EPSILON_0
   !write(*,*) atom%weight, Helium%weight, Hydrogen%weight, atmos%Elements(2)%ptr_elem%weight
  vrel35_He = (8.*KBOLTZMANN/(PI*AMU*atom%weight) * (1.+atom%weight/Helium%weight))**0.3

  Z = atom%stage(j)+1
  ic = find_Continuum(line%atom, j)

  deltaR = (E_RYDBERG/(atom%E(ic)-atom%E(j)))**2. - (E_RYDBERG/(atom%E(ic)-atom%E(i)))**2.
  C625 = (2.5*((Q_ELECTRON)**2./fourPIeps0)*(ABARH/fourPIeps0) * 2.*PI*(Z*RBOHR)**2./HPLANCK * deltaR)**(4d-1)


  SELECT CASE (line%vdWaals)
  CASE ("UNSOLD")
   ! Relative velocity of radiator and perturber with Maxwellian
   ! velocity distributions
   vrel35_H = (8.*KBOLTZMANN/(PI*AMU*atom%weight) * (1.+atom%weight/Hydrogen%weight))**(3d-1)
   cross = 8.08 * (line%cvdWaals(1)*vrel35_H+line%cvdWaals(3)*Helium%abund*vrel35_He)*C625
   !write(*,*) Helium%abund, C625
   GvdW = cross * atmos%T(icell)**(3d-1)
  CASE ("BARKLEM")
   write(*,*) "Warning-> vdWaals broadening following BARKLEM ",&
     "not tested yet!"
   ! UNSOLD for Helium
   cross = 8.08 * line%cvdWaals(3)*Helium%abund*vrel35_He*C625
    GvdW = line%cvdWaals(1) * atmos%T(icell)**(real(1.-&
          line%cvdWaals(2),kind=8)/2.) + cross*atmos%T(icell)**(3d-1)
  CASE DEFAULT
   write(*,*) "Method for van der Waals broadening unknown"
   write(*,*) "exiting..."
   stop
  END SELECT

   !write(*,*) "GvdW=", GvdW(k)
   GvdW = GvdW * atmos%nHtot(icell)
   !GvdW = GvdW * Hydrogen%n(1,icell) !ground level pops
   !GvdW = GvdW * sum(Hydrogen%n(1:Hydrogen%Nlevel-1,icell)) !total nH_I pops

 RETURN
 END SUBROUTINE VanderWaals

 SUBROUTINE Stark(icell, atom, kr, GStark)
 ! Quadratic Stark broadening by electrons and singly charged
 ! ions.
 !
 ! Gamma4 = 11.37 * vrel^1/3 * C_4^2/3 * (ne + nion)
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
   cStark = dabs(real(line%cStark,kind=dp))
    GStark = cStark * atmos%ne(icell)
  else
   ! Constants for relative velocity. Assume that nion = ne
   ! and that the average atomic weight of ionic pertubers is
   ! given by AVERAGE_ATOMIC_WEIGHT

   C = 8. * KBOLTZMANN / (PI*AMU*atom%weight)
   Cm = (1.+atom%weight/melONamu)**(1.6666667d-1) + &
      (1.+atom%weight/AVERAGE_ATOMIC_WEIGHT)**(1.6666667d-1)

   ! Find core charge Z and effective quantumer numbers neff_u
   ! and neff_l for upper and lower level
   Z = atom%stage(line%i) + 1 !physical, not index
   ic = find_continuum(line%atom, line%i)
   
   ERYDBERG = E_RYDBERG / (1.+M_ELECTRON / (atom%weight*AMU)) !corection
   neff_l = Z*dsqrt(ERYDBERG/(atom%E(ic)-atom%E(line%i)))
   neff_u = Z*dsqrt(ERYDBERG/(atom%E(ic)-atom%E(line%j)))
   
    C4 = ((Q_ELECTRON)**2./(4.*PI*EPSILON_0))*RBOHR *   &
     (2.*PI*RBOHR**2./HPLANCK) / (18.*Z*Z*Z*Z) * &
     ((neff_u*(5.0*(neff_u)**2. + 1.0))**2. - (neff_l*(5.0*(neff_l)**2. + 1.0))**2.)
     
    cStark23 = 11.37 * (line%cStark * C4)**(6.6666667d-1)

    vrel = (C*atmos%T(icell))**(1.6666667d-1) * Cm
    GStark = cStark23 * vrel * atmos%ne(icell)
  end if

 RETURN
 END SUBROUTINE Stark

 SUBROUTINE StarkLinear(icell, atom, kr, GStark)
 ! Linear Stark broadening by electrons for hydrogen lines.
 !
 ! See: K. Sutton (1978), JQSRT 20, 333-343 (formula for z = nu)
 !
 ! GStark = a_1 * &
 !       [0.60 * (n_u^2 - n_l^2) * (N_e)^(2/3) * CM_TO_M^2] !s^-1, with N_e here in cm^-3
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

  if ((n_upper-n_lower).eq.1) then !first line of the serie, dn=1
   a1 = 0.642
  else
   a1 = 1.
  end if
  C = a1 * 0.6 * (n_upper*n_upper-n_lower*n_lower) * (CM_TO_M)**2.

   GStark = C*atmos%ne(icell)**(2d0/3d0)

 RETURN
 END SUBROUTINE StarkLinear

 SUBROUTINE Damping(icell, atom, kr, adamp)
 ! I need to change something here, only a line from Atom%lines is needed since,
 !line%atom point to atom
 ! Beware here, because if I remove some transitions for images, but I need to do the calculations of
 !broadening by summing over transitions..
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (AtomicLine) :: line
  integer, intent(in) :: kr
  integer :: k
  double precision :: cDop
  double precision, intent(out) :: adamp
  double precision :: Qelast, Gj

  Qelast = 0.0

  line=atom%lines(kr)
  !write(*,*) "Computing damping for line ",line%j,"-",line%i

!   if (atom%g(line%j)==18.and.atom%g(line%i)==8) then
!   CALL RadiativeDamping(icell, atom, kr, Gj)
!   write(*,*) "Grad compare:", line%Grad, Gj
!   
!   stop
!   endif

  cDop = (NM_TO_M*line%lambda0) / (4.*PI)
  ! van der Waals broadening
  if ((line%cvdWaals(1).gt.0).or.(line%cvdWaals(3).gt.0)) then
   CALL VanderWaals(icell, atom, kr, adamp)
    Qelast = Qelast + adamp
    !write(*,*) " Van der Waals:", adamp
    !write(*,*),"vdW", adamp* cDop/VBROAD_atom(icell,atom)
  end if
  ! Quadratic Stark broadening
  if (line%cStark > 0.) then
   CALL Stark(icell, atom, kr, adamp)
    Qelast = Qelast + adamp
    !write(*,*) " Stark qu:", adamp
    !write(*,*),"Stark2", adamp* cDop/VBROAD_atom(icell,atom)
  end if
  ! Linear Stark broadening only for Hydrogen
  if (atom%ID.eq."H ") then
   CALL StarkLinear(icell, atom,kr, adamp)
    Qelast = Qelast + adamp
    !write(*,*) " Stark linear:", adamp
    !write(*,*),"Stark1", adamp* cDop/VBROAD_atom(icell,atom)
  end if
  ! Store Qelast, total rate of elastic collisions, if PFR only
  if (line%PFR) then
   write(*,*) "PFR not handled yet, avoiding"
   ! line%Qelast = Qelast
  end if
  !write(*,*) "Grad=", line%Grad * cDop/VBROAD_atom(icell,atom), cDop/VBROAD_atom(icell,atom)
  adamp = (line%Grad + Qelast)*cDop / VBROAD_atom(icell,atom)


 RETURN
 END SUBROUTINE Damping

END MODULE broad
