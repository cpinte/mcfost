
MODULE broad

 use constant
 use atom_type
 use atmos_type, only : T, ne, elements, Hydrogen, Helium !alias to ptr_atom of Helium, which might not exist
 
 use messages, only : error, warning
 use mcfost_env, only : dp

 IMPLICIT NONE

 real, parameter :: AVERAGE_ATOMIC_WEIGHT = 28.0

 CONTAINS

 SUBROUTINE Broad_Kurosawa(icell, atom, kr, adamp)
 ! For testing, implementing the broadening for H alpha used by some people
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: kr, icell
  real(kind=dp), intent(out) :: adamp
  real(kind=dp) :: Crad, Cvdw, Cstark, G
  
  !AA
  Crad = 6.5d-4; CvdW = 4.4d-4; Cstark = 1.17d-3
  !in AA
  G = Crad + Cvdw * (1d-6 * sum(atom%n(1:atom%Nlevel-1,icell))/1d16) * (T(icell)/5d3)**(0.3) + &
      Cstark * (ne(icell)/1d12)**(2./3.) * 1d-4 !conversion from 12cm^-3 ^ 2/3 to 12m^-3 ^ 2/3
  G = G * 1d-10 !in m
  !conversion in s^-1 -> G(m) = G(s^-1) * lambda**2 / c : G(s^-1)/nu = G(m)/lambda -> G(m)*nu/lambda
  adamp  = G * CLIGHT / atom%lines(kr)%lambda0**2 * 1d18!s^-1
  !will be divided by 4PI later
 RETURN
 END SUBROUTINE Broad_Kurosawa
 
 
 SUBROUTINE Radiative_Damping(atom, line, gamma_j, gamma_i)
 ! -------------------------------------------------------------------------- !
  ! Radiative damping of a line transition j->i
  ! Already included in model atom.
  ! check on H
  ! need a check on other atoms, like Na and Ca.
  ! Do not use if the model atom is small, because the damping would be underestimated
  ! even by an order of magnitude
 ! -------------------------------------------------------------------------- !
  type (AtomType), intent(in) :: atom
  type (AtomicLine), intent(in) :: line
  real(kind=dp) :: Grad
  integer :: kp,l,n
  real(kind=dp), intent(out) :: gamma_j, gamma_i
  
  Grad = 0d0
  gamma_j = 0d0
  gamma_i = 0.0_dp
  
  !it also adds line%Aji for the current line if l  = line%i == the same line
!   write(*,*) "Radiative damping for line ", line%j, line%i
  
  do kp=1,atom%Nline !sum over all lines

    l = atom%lines(kp)%i; n = atom%lines(kp)%j
    !write(*,*) "------"
    !write(*,*) "l=",l, "u=",n
    !write(*,*) line%i, line%j
    !write(*,*) "------"
    if (n==line%j) then !our upper level is also the upper level of another transition
!      write(*,*) "j is upper level of ", line%j, n, l
     gamma_j = gamma_j + atom%lines(kp)%Aji
    endif
    
    if (n==line%i) then !our lower level is also the upper level of other transitions
!      write(*,*) "i is upper level of ", line%i, n, l
     gamma_i = gamma_i + atom%lines(kp)%Aji
    endif
  
  enddo
  
  Grad = gamma_j + gamma_i
!   write(*,*) "Grad=", Grad/1d8, gamma_j/1d8, gamma_i/1d8

 RETURN
 END SUBROUTINE Radiative_Damping
 
 SUBROUTINE VanderWaals_new(icell, atom, kr, GvdW)
 ! Compute van der Waals broadening in Lindholm theory with
 !Unsold's approximation for the interaction coefficient C_6.
 !
 !Gamma = 8.08 * vrel^3/5 * C_6^2/5 * nH(1)
 !
 ! vrel is the maxwellian velocity of the perturber
 !
  real(kind=dp), intent(out) :: GvdW
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (Element), pointer :: Helium
  integer, intent(in) :: kr
  integer :: i, j, ic, Z
  real(kind=dp) :: vrel35_H, vrel35_He, deltaR
  real(kind=dp) :: C625, mu, mub, nH, nHe, alpha
  
  alpha = ABARH / fourpi / EPSILON_0!polarizability of H in m^3

  nH = hydrogen%n(1,icell)
  if (associated(helium)) then
   nHe = helium%n(1,icell)
   !if helium is a perturber
   mub = atom%weight * helium%weight / (atom%weight + helium%weight)
   vrel35_He =(vtherm/mub * T(icell))**(0.3)
  else
   nHe = 0.0
   vrel35_He = 0.0
  endif
    
  j = atom%lines(kr)%j
  i = atom%lines(kr)%i

  Z = atom%stage(j)+1
  ic = find_continuum(atom, j)

  !reduce mass for mean velocity
  mu = atom%weight * hydrogen%weight / (atom%weight + hydrogen%weight)
  
  !square
  !in total there is a 2 square factor
  deltaR = RBOHR**2 * 2.5 / Z / Z * (n_eff(atom%Rydberg, atom%E(ic), atom%E(j), Z)**4 - &
  									n_eff(atom%Rydberg, atom%E(ic), atom%E(i), Z)**4)
  									
  C625 = (2.0*pi/fourpi/EPSILON_0)**(0.4) * (Q_ELECTRON*Q_ELECTRON/HPLANCK * alpha * deltaR)**(0.4) 
  vrel35_H = (vtherm/mu * T(icell))**(0.3) !sqvel**0.3 = vel**(3./5.)

  !enhancement 
  !*atom%lines(kr)%cvdWaals(1), * atom%lines(kr)%cvdWaals(3)
  Gvdw = 8.08 * (nH * vrel35_H + nHe * vrel35_He) * C625 
  
  !GvdW = 10**(6.33 + 0.4 * log10(n_eff(atom%Rydberg, atom%E(ic), atom%E(j), Z)**4 - &
  !			n_eff(atom%Rydberg, atom%E(ic), atom%E(i), Z)**4) -0.7*log10(T(icell)) + &
  !			log10(nHtot(icell)*KBOLTZMANN*T(icell)*0.1))


 RETURN
 END SUBROUTINE VanderWaals_new


 SUBROUTINE VanderWaals(icell, atom, kr, GvdW)
 ! Compute van der Waals broadening in Lindholm theory with
 !Unsold's approximation for the interaction coefficient C_6.
 !
 !Gamma = 8.08 * vrel^3/5 * C_6^2/5 * nH(1)
 !
 !
 ! Or with the Anstee, Barklem & O'Mara formalism (not yet)
 !
 ! See: Anstee & O'Mara 1995, MNRAS 276, 859-866
 ! Barklem & O'Mara 1998, MNRAS 300, 863-871
 !
 ! Gamma =
 ! (4/pi)^(alpha/2) * &
 !     Gam((4-alpha)/2) v_0 sigma(v_0)(vmean/v_0)^(1-alpha)
 !
  real(kind=dp), intent(out) :: GvdW
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  type (Element), pointer :: Helium_elem
  integer, intent(in) :: kr
  integer :: i, j, ic, Z, k
  real(kind=dp) :: vrel35_H, vrel35_He, fourPIeps0, deltaR
  real(kind=dp) :: cross, gammaCorrH, gammaCorrHe, C625


  Helium_elem => Elements(2)%ptr_elem ! Helium basic informations


  j = atom%lines(kr)%j
  i = atom%lines(kr)%i

  !could be nice to keep vrel35 for each atom ? as it depends only on the weight

  fourPIeps0 = 4.*PI*EPSILON_0
  vrel35_He = (8.*KBOLTZMANN/(PI*AMU*atom%weight) * (1.+atom%weight/Helium_elem%weight))**0.3

  Z = atom%stage(j)+1
  ic = find_continuum(atom%lines(kr)%atom, j)

  deltaR = (atom%Rydberg/(atom%E(ic)-atom%E(j)))**2. - (atom%Rydberg/(atom%E(ic)-atom%E(i)))**2.
  C625 = (2.5*((Q_ELECTRON)**2./fourPIeps0)*(ABARH/fourPIeps0) * 2.*PI*(Z*RBOHR)**2./HPLANCK * deltaR)**(4d-1)


  SELECT CASE (atom%lines(kr)%vdWaals)
  CASE ("UNSOLD")
   ! Relative velocity of radiator and perturber with Maxwellian
   ! velocity distributions
   vrel35_H = (8.*KBOLTZMANN/(PI*AMU*atom%weight) * (1.+atom%weight/Hydrogen%weight))**(3d-1)
   cross = 8.08 * (atom%lines(kr)%cvdWaals(1)*vrel35_H+atom%lines(kr)%cvdWaals(3)*Helium_elem%abund*vrel35_He)*C625
   GvdW = cross * T(icell)**(3d-1)
   
  CASE ("BARKLEM")
   write(*,*) "Warning-> vdWaals broadening following BARKLEM ",&
     "not tested yet!"
   ! UNSOLD for Helium
   cross = 8.08 * atom%lines(kr)%cvdWaals(3)*Helium_elem%abund*vrel35_He*C625
    GvdW = atom%lines(kr)%cvdWaals(1) * T(icell)**(real(1.-&
          atom%lines(kr)%cvdWaals(2),kind=8)/2.) + cross*T(icell)**(3d-1)
  CASE DEFAULT
   write(*,*) "Method for van der Waals broadening unknown"
   write(*,*) "exiting..."
   stop
  END SELECT


   GvdW = GvdW * Hydrogen%n(1,icell) !total nH_I in the ground state

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
  real(kind=dp), intent(out) :: GStark
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: kr
  integer :: k, ic, Z
  real(kind=dp) :: C4, C, Cm, melONamu, cStark23, cStark
  real(kind=dp) :: neff_u, neff_l, vrel, ERYDBERG

  melONamu = M_ELECTRON / AMU
  
  if ( atom%lines(kr)%cStark < 0.) then
   cStark = dabs(real( atom%lines(kr)%cStark,kind=dp))
    GStark = cStark * ne(icell)
  else
   ! Constants for relative velocity. Assume that nion = ne
   ! and that the average atomic weight of ionic pertubers is
   ! given by AVERAGE_ATOMIC_WEIGHT

   C = 8. * KBOLTZMANN / (PI*AMU*atom%weight)
   Cm = (1.+atom%weight/melONamu)**(1.6666667d-1) + &
      (1.+atom%weight/AVERAGE_ATOMIC_WEIGHT)**(1.6666667d-1)

   ! Find core charge Z and effective quantumer numbers neff_u
   ! and neff_l for upper and lower level
   Z = atom%stage( atom%lines(kr)%i) + 1 !physical, not index
   ic = find_continuum(atom, atom%lines(kr)%i)
   
   ERYDBERG = E_RYDBERG / (1.+M_ELECTRON / (atom%weight*AMU)) !corection
   neff_l = Z*sqrt(ERYDBERG/(atom%E(ic)-atom%E( atom%lines(kr)%i)))
   neff_u = Z*sqrt(ERYDBERG/(atom%E(ic)-atom%E( atom%lines(kr)%j)))
   
    C4 = ((Q_ELECTRON)**2./(4.*PI*EPSILON_0))*RBOHR *   &
     (2.*PI*RBOHR**2./HPLANCK) / (18.*Z*Z*Z*Z) * &
     ((neff_u*(5.0*(neff_u)**2. + 1.0))**2. - (neff_l*(5.0*(neff_l)**2. + 1.0))**2.)
     
    cStark23 = 11.37 * ( atom%lines(kr)%cStark * C4)**(6.6666667d-1)

    vrel = (C*T(icell))**(1.6666667d-1) * Cm
    GStark = cStark23 * vrel * ne(icell)
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
 ! This is highly inaccurate for lines with n > 2
 !
  real(kind=dp), intent(out) :: GStark
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: kr
  real(kind=dp) :: n_upper, n_lower, Gstark2
  real :: Z, K, dn, a1, gs
  real(kind=dp) :: C, nelectric
  
  Z = 1.0 + atom%stage(atom%lines(kr)%i)
  
  !include proton ?
  nelectric = 1d-6 * (ne(icell) + 0.0 * hydrogen%n(hydrogen%Nlevel,icell)) !cm^-3

  n_lower = sqrt(atom%g(atom%lines(kr)%i)/2.)
  n_upper = sqrt(atom%g(atom%lines(kr)%j)/2.) !works because stark linear only for Hydrogen.
  									 !at the moment. otherwise use n_eff, wich works for H
  									 
  dn = real(n_upper - n_lower)

  if (dn == 1.0) then !first line of the serie, dn=1
   a1 = 0.642
  else
   a1 = 1.0
  end if
  
  !Sutton recommends to multiply dz with gs=0.425 if to be used with a dispersion profile
  !but since it is very inaccurate and we miss broadening 
  gs = 1.0 !

  GStark = gs * a1 * 0.6 * (n_upper*n_upper-n_lower*n_lower) * ( nelectric**(2./3.) )



 RETURN
 END SUBROUTINE StarkLinear
 
 !building
 subroutine resonance_broadening(icell, atom, kr, adamp)
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: kr
  real(kind=dp), intent(out) :: adamp
  real(kind=dp) :: gr, nu1j, a1, a2, f1j, np
  integer :: ic, i, j
  
  i = atom%lines(kr)%i
  j = atom%lines(kr)%j
  
  a1 = 0.0
  a2 = 0.0
  
  np = 1.0
  ic = find_continuum(atom, i)
  
  if (i==1) then !transition between ground state and excited levels
   gr = atom%g(atom%lines(kr)%i)/atom%g(atom%lines(kr)%j)
   nu1j = M_TO_NM * CLIGHT / atom%lines(kr)%lambda0 !ground state is i
   a1 = sqrt(gr) * atom%lines(kr)%fosc / nu1j
  else !transition between two exited states
  	   !only one is resonance boradened
   nu1j = (atom%E(i) - atom%E(1)) / HPLANCK
   gr = atom%g(1) / atom%g(i)
   np = n_eff(atom%rydberg, atom%E(ic),atom%E(i), atom%stage(ic))
   f1j = line_oscillator_strength(1.0_dp, np)
   a1 = f1j * sqrt(gr)/nu1j
!    nu1j = (atom%E(atom%lines(kr)%j) - atom%E(1)) / HPLANCK
!    gr = atom%g(1) / atom%g(atom%lines(kr)%j)
!    np = n_eff(atom%rydberg, atom%E(ic),atom%E(atom%lines(kr)%j), atom%stage(atom%lines(kr)%j)+1)
!    f1j = line_oscillator_strength(1.0_dp, np)   
!    a2 = f1j * sqrt(gr)/nu1j
  endif

  adamp = 5.48 * pi * Q_ELECTRON**2  * max(0.0, max(a1,a2)) * atom%n(1,icell) / M_ELECTRON!/ fourpi
 
 return
 end subroutine resonance_broadening

 SUBROUTINE Damping(icell, atom, kr, adamp)

  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: kr
  integer :: k
  real(kind=dp) :: cDop
  real(kind=dp), intent(out) :: adamp
  real(kind=dp) :: Qelast, Gj

  Qelast = 0.0 !used for PFR
  !stored everything that is not radiative

  !conversion factor from gamma in sr / s^-1 to Gamma in m/s
  !and then to adamp as Gamma / vrboad
  cDop = (NM_TO_M*atom%lines(kr)%lambda0) / (4.*PI)
  

  ! van der Waals broadening
  !interaction with neutral H and neutral He
  !Not for H and He ?
  if ((atom%lines(kr)%cvdWaals(1) > 0).or.(atom%lines(kr)%cvdWaals(3) > 0)) then
   CALL VanderWaals(icell, atom, kr, adamp)
   !CALL VanderWaals_new(icell, atom, kr, adamp)
    Qelast = Qelast + adamp
  end if

  ! Quadratic Stark broadening
  !Interaction with charged particles
  !-> for H cStark is 0
  if (atom%lines(kr)%cStark /= 0.0) then
   CALL Stark(icell, atom, kr, adamp)
    Qelast = Qelast + adamp
  end if

  ! Linear Stark broadening only for Hydrogen
  !Approximate treatment at the moment
  if (atom%ID == "H") then
   CALL StarkLinear(icell, atom,kr, adamp)
    Qelast = Qelast + adamp
!    call resonance_broadening(icell, atom, kr, adamp)
!     Qelast = Qelast + adamp
  end if
  
  adamp = (atom%lines(kr)%Grad + Qelast)*cDop / atom%vbroad(icell)
!   write(*,*) atom%ID, atom%lines(kr)%j, atom%lines(kr)%i, atom%vbroad(icell)/1e3
!   write(*,*) T(icell), "a=", adamp, "Grad", atom%lines(kr)%Grad


 RETURN
 END SUBROUTINE Damping

END MODULE broad
