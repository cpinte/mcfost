module broad

   use constantes
   use atom_type
   use elements_type
   use grid, only : T, ne, vturb
   use mcfost_env, only : dp

   implicit none

   real, parameter :: AVERAGE_ATOMIC_WEIGHT = 28.0

   contains


   function line_damping(icell, line)
      !main function to compute line broadening
      !for each transition individually
      integer, intent(in) :: icell
      type (AtomicLine), intent(in) :: line
      integer :: k
      real(kind=dp) :: radfreq_to_vel
      real(kind=dp) :: adamp, vdw, vth
      real(kind=dp) :: line_damping
      real(kind=dp) :: Qelast, Gj

      line_damping = line%Grad
      vth = vbroad(T(icell), line%atom%weight, vturb(icell))

      !conversion factor from gamma in sr / s^-1 to Gamma in m/s
      !and then to adamp as Gamma / vrboad
      radfreq_to_vel = (NM_TO_M*line%lambda0) / (4.0*pi)


      ! van der Waals broadening
      !interaction with neutral H and neutral He

      if ((line%cvdWaals(1) > 0).or.(line%cvdWaals(3) > 0)) then
         call VanderWaals_line(icell, line, adamp)
         line_damping = line_damping + adamp
         vdw = adamp !store for checking
      end if

      ! Quadratic Stark broadening
      !Interaction with charged particles
      !-> for H cStark is 0
      if (line%cStark /= 0.0) then
         call Stark_line(icell, line, adamp)
         line_damping = line_damping + adamp
      end if

      ! Linear Stark broadening only for Hydrogen
      !Approximate treatment at the moment
      if (line%atom%ID == "H") then
         call StarkLinear_line(icell, line, adamp)
         line_damping = line_damping + adamp
      end if
      !unit less
      line_damping = line_damping * radfreq_to_vel / vth

      return
   end function line_damping

   subroutine VanderWaals_line(icell, line, GvdW)
      ! Compute van der Waals broadening in Lindholm theory with
      ! Unsold's approximation for the interaction coefficient C_6.
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
      type (AtomicLine), intent(in) :: line
      type (Element), pointer :: Helium_elem
      integer :: i, j, ic, Z, k
      real(kind=dp) :: vrel35_H, vrel35_He, fourPIeps0, deltaR
      real(kind=dp) :: cross, gammaCorrH, gammaCorrHe, C625

      j = line%j
      i = line%i


      fourPIeps0 = 4.0*pi*EPSILON_0
      vrel35_He = (8.*KB/(PI*amu_kg*line%atom%weight) * (1.+line%atom%weight/Elems(2)%weight))**0.3

      Z = line%atom%stage(j)+1
      ic = find_continuum(line%atom, j)

      deltaR = (line%atom%Rydberg/(line%atom%E(ic)-line%atom%E(j)))**2. - (line%atom%Rydberg/(line%atom%E(ic)-line%atom%E(i)))**2.
      C625 = (2.5*((electron_charge)**2./fourPIeps0)*(ABARH/fourPIeps0) * 2.*PI*(Z*RBOHR)**2./HP * deltaR)**(0.4)


      select case (line%vdWaals)
         case ("UNSOLD")
            ! Relative velocity of radiator and perturber with Maxwellian
            ! velocity distributions
            vrel35_H = (8.*KB/(PI*amu_kg*line%atom%weight) * (1.+line%atom%weight/Hydrogen%weight))**(0.3)
            cross = 8.08 * (line%cvdWaals(1)*vrel35_H+line%cvdWaals(3)*Elems(2)%abund*vrel35_He)*C625
            GvdW = cross * T(icell)**(3d-1)
         case ("BARKLEM")
            cross = 8.08 * line%cvdWaals(3)*Helium_elem%abund*vrel35_He*C625 * T(icell)**(0.3)  !Unsold
            GvdW = line%cvdWaals(1) * T(icell)**(0.5 - 0.5*line%cvdWaals(2)) + cross
         case default
            write(*,*) "Method for van der Waals broadening unknown", line%vdwaals
            write(*,*) "exiting..."
            stop
      end select

      GvdW = GvdW * Hydrogen%n(1,icell)!total nH_I in the ground state

      return
   end subroutine VanderWaals_line

   subroutine Stark_line(icell,line, GStark)
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
      type (AtomicLine), intent(in) :: line
      integer :: k, ic, Z
      real(kind=dp) :: C4, C, Cm, melONamu, cStark23, cStark
      real(kind=dp) :: neff_u, neff_l, vrel, ERYDBERG

      melONamu = Mel / amu_kg

      if ( line%cStark < 0.) then
         cStark = abs(real( line%cStark,kind=dp))
         GStark = cStark * ne(icell)
      else
         ! Constants for relative velocity. Assume that nion = ne
         ! and that the average atomic weight of ionic pertubers is
         ! given by AVERAGE_ATOMIC_WEIGHT

         C = 8. * KB / (PI*amu_kg*line%atom%weight)
         Cm = (1.+line%atom%weight/melONamu)**(1.6666667d-1) + &
            (1.+line%atom%weight/AVERAGE_ATOMIC_WEIGHT)**(1.6666667d-1)

         ! Find core charge Z and effective quantumer numbers neff_u
         ! and neff_l for upper and lower level
         Z = line%atom%stage( line%i) + 1 !physical, not index
         ic = find_continuum(line%atom, line%i)

         ERYDBERG = E_RYDBERG / (1.+Mel / (line%atom%weight*amu_kg)) !corection
         neff_l = Z*sqrt(ERYDBERG/(line%atom%E(ic)-line%atom%E( line%i)))
         neff_u = Z*sqrt(ERYDBERG/(line%atom%E(ic)-line%atom%E( line%j)))

         C4 = ((electron_charge)**2./(4.*PI*EPSILON_0))*RBOHR *   &
            (2.*PI*RBOHR**2./HP) / (18.*Z*Z*Z*Z) * &
            ((neff_u*(5.0*(neff_u)**2. + 1.0))**2. - (neff_l*(5.0*(neff_l)**2. + 1.0))**2.)

         cStark23 = 11.37 * ( line%cStark * C4)**(6.6666667d-1)

         vrel = (C*T(icell))**(1.6666667d-1) * Cm
         GStark = cStark23 * vrel * ne(icell)
      end if

      return
   end subroutine Stark_line

   subroutine StarkLinear_line(icell, line, GStark)
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
      type (AtomicLine), intent(in) :: line
      real(kind=dp) :: n_upper, n_lower, Gstark2
      real :: Z, K, dn, a1, gs
      real(kind=dp) :: C, nelectric

      Z = 1.0 + line%atom%stage(line%i)

      !include proton ?
      nelectric = 1d-6 * (ne(icell) + 0.0 * hydrogen%n(hydrogen%Nlevel,icell)) !cm^-3

      n_lower = sqrt(line%atom%g(line%i)/2.)
      n_upper = sqrt(line%atom%g(line%j)/2.) !works because stark linear only for Hydrogen.
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

      return
   end subroutine StarkLinear_line

   ! subroutine Radiative_Damping(atom, line, gamma_j, gamma_i)
   !    ! -------------------------------------------------------------------------- !
   !    ! Radiative damping of a line transition j->i
   !    ! Already included in model atom.
   !    ! check on H
   !    ! need a check on other atoms, like Na and Ca.
   !    ! Do not use if the model atom is small, because the damping would be underestimated
   !    ! even by an order of magnitude.
   !    ! Stimulated emission is neglected, the lower lovel is infinitely sharp.
   !    ! -------------------------------------------------------------------------- !
   !    type (AtomType), intent(in) :: atom
   !    type (AtomicLine), intent(in) :: line
   !    real(kind=dp) :: Grad
   !    integer :: kp,l,n
   !    real(kind=dp), intent(out) :: gamma_j, gamma_i

   !    Grad = 0d0
   !    gamma_j = 0d0
   !    gamma_i = 0.0_dp

   !    do kp=1,atom%Nline !sum over all lines

   !       l = atom%lines(kp)%i; n = atom%lines(kp)%j
   !       if (n==line%j) then !our upper level is also the upper level of another transition
   !          gamma_j = gamma_j + atom%lines(kp)%Aji
   !       endif

   !    enddo

   !    Grad = gamma_j + gamma_i
   !    return
   ! end subroutine Radiative_Damping

  ! !building
   ! subroutine resonance_broadening(icell, atom, kr, adamp)
   !    integer, intent(in) :: icell
   !    type (AtomType), intent(in) :: atom
   !    integer, intent(in) :: kr
   !    real(kind=dp), intent(out) :: adamp
   !    real(kind=dp) :: gr, nu1j, a1, a2, f1j, np
   !    integer :: ic, i, j

   !    i = atom%lines(kr)%i
   !    j = atom%lines(kr)%j

   !    a1 = 0.0
   !    a2 = 0.0

   !    np = 1.0
   !    ic = find_continuum(atom, i)

   !    if (i==1) then !transition between ground state and excited levels
   !       gr = atom%g(atom%lines(kr)%i)/atom%g(atom%lines(kr)%j)
   !       nu1j = M_TO_NM * CLIGHT / atom%lines(kr)%lambda0 !ground state is i
   !       a1 = sqrt(gr) * atom%lines(kr)%fosc / nu1j
   !    else !transition between two exited states
   !       !only one is resonance boradened
   !       nu1j = (atom%E(i) - atom%E(1)) / HPLANCK
   !       gr = atom%g(1) / atom%g(i)
   !       np = n_eff(atom%rydberg, atom%E(ic),atom%E(i), atom%stage(ic))
   !       f1j = line_oscillator_strength(1.0_dp, np)
   !       a1 = f1j * sqrt(gr)/nu1j
   !     !    nu1j = (atom%E(atom%lines(kr)%j) - atom%E(1)) / HPLANCK
   !     !    gr = atom%g(1) / atom%g(atom%lines(kr)%j)
   !     !    np = n_eff(atom%rydberg, atom%E(ic),atom%E(atom%lines(kr)%j), atom%stage(atom%lines(kr)%j)+1)
   !     !    f1j = line_oscillator_strength(1.0_dp, np)
   !     !    a2 = f1j * sqrt(gr)/nu1j
   !    endif

   !    adamp = 5.48 * pi * Q_ELECTRON**2  * max(0.0, max(a1,a2)) * atom%n(1,icell) / M_ELECTRON!/ fourpi

   !    return
   ! end subroutine resonance_broadening

end module broad