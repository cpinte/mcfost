! Hydrogen, bound-free, free-free and continuum opacities
! - H b-f and f-f --> H b-f is now computed in metal.f90 with the other atoms
! - H^- b-f and f-f
!  free-free emissivity is obtained by
! -> eta_ff = chi_ff * Bplanck
!
!
! Opacities in m^-1 (chi)
! Emissivities in J/s/m3/Hz/sr (eta)
MODULE hydrogen_opacities

 use atom_type, only : AtomicContinuum, find_continuum, AtomType
 use atmos_type, only : atmos, Hydrogen, Helium
 use constant
 use spectrum_type, only : NLTEspec
 use math, only : bezier3_interp, interp2Darr
 use lte, only : nH_minus
 
 use constantes, only : tiny_dp
 use mcfost_env, only : dp

 IMPLICIT NONE

 integer, parameter :: NBF=34, NJOHN=6, NFF=17, NTHETA=16

 CONTAINS

 ELEMENTAL FUNCTION Gaunt_bf(u, n_eff)
  ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
  ! See also Menzel & Pekeris 1935
  real(kind=dp), intent(in) :: n_eff,  u ! = n_Eff**2 * eps = hnu/Z/Z/E_RYDBERG - 1
  real(kind=dp) :: Gaunt_bf
  Gaunt_bf = 1d0 + 0.1728 * (n_eff**(-2./3.)) * (u+1d0)**(-2./3.) * (u-1d0) &
             - 0.0496*(n_eff**(-4./3.)) * (u+1d0)**(-4./3.) * (u*u + 4./3. * u + 1d0)
             
  if (Gaunt_bf <= 0d0 .or. Gaunt_bf > 2d0) then
   Gaunt_bf = 1d0
  end if

 RETURN
 END FUNCTION Gaunt_bf
  

 Elemental FUNCTION Gaunt_ff(lambda, Z, T)
 ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
 !
 ! Note: There is a problem with this expansion at higher temperatures
 ! (T > 3.0E4 and longer wavelengths (lambda > 2000 nm). Set to
 ! 1.0 when the value goes below 1.0
  real(kind=dp),intent(in) :: lambda, T
  real(kind=dp) :: x, x3, y, Gaunt_ff
  integer, intent(in) :: Z

  x = ((HPLANCK * CLIGHT)/(lambda * NM_TO_M)) / &
       (E_RYDBERG * (Z)**(2d0))
  x3 = (x**(3.3333333d-1))
  y  = (2.0 * lambda * NM_TO_M * KBOLTZMANN*T) / &
       (HPLANCK*CLIGHT)

  gaunt_ff = 1.0 + 0.1728*x3 * (1.0 + y) - &
        0.0496*(x3*x3) * (1.0 + (1.0 + y)*0.33333333*y)

  if (gaunt_ff <= 0d0 .or. gaunt_ff > 2d0) gaunt_ff = 1d0

 RETURN
 END FUNCTION Gaunt_ff

 ELEMENTAL FUNCTION  H_bf_Xsection(cont, lambda) result(alpha) !_lambda
  Type (AtomicContinuum), intent(in) :: cont
  real(kind=dp), intent(in) :: lambda
  real(kind=dp) :: n_eff, g_bf, u, Z, u0, g_bf0, alpha
  
  
   Z = real(cont%atom%stage(cont%i) + 1,kind=dp)
   n_eff = Z*dsqrt(E_RYDBERG / (cont%atom%E(cont%j) - cont%atom%E(cont%i))) 
   

    u = n_eff**2 * HPLANCK*CLIGHT / (NM_TO_M * lambda) / Z*Z / E_RYDBERG - 1
    u0 = n_eff*n_eff * HPLANCK*CLIGHT / (NM_TO_M * cont%lambda0) / Z / Z / E_RYDBERG - 1.
    
    g_bf = Gaunt_bf(u, n_eff)
    g_bf0 = Gaunt_bf(u0, n_eff)
    
    !There is a factor n_eff/Z**2 absorbed in alpha0 beware
    !alpha = cont%alpha0 * (lambda/cont%lambda0)**3  * g_bf / g_bf0 
 	alpha = 1d-4* 2.815d29 * (Z**4) * g_bf /n_eff**5 * (NM_TO_M*lambda/CLIGHT)**3

 RETURN
 END FUNCTION H_bf_Xsection!_lambda
 
 SUBROUTINE  test_bf_xs(cont, lambda)
  !test the computation of hydrogenic bound-free cross-sections
  Type (AtomicContinuum), intent(in) :: cont
  real(kind=dp), intent(in) :: lambda
  real(kind=dp) :: n_eff, g_bf, u, Z, u0, g_bf0, alpha
  
  
   Z = real(cont%atom%stage(cont%i) + 1,kind=dp)
!    if (cont%atom%ID == "H ") then
!       n_eff = dsqrt(Hydrogen%g(cont%i)/2.)  !only for Hydrogen !
!       !n_eff = Z*dsqrt(E_RYDBERG / (cont%atom%E(cont%j) - cont%atom%E(cont%i))) 
!    else
!      !obtained_n = getPrincipal(metal%label(continuum%i), n_eff)
!      !if (.not.obtained_n) &
!      n_eff = Z*dsqrt(E_RYDBERG / (cont%atom%E(cont%j) - cont%atom%E(cont%i))) 
!    end if
   n_eff = Z*dsqrt(E_RYDBERG / (cont%atom%E(cont%j) - cont%atom%E(cont%i))) 
   

    u = n_eff**2 * HPLANCK*CLIGHT / (NM_TO_M * lambda) / Z*Z / E_RYDBERG - 1
    u0 = n_eff*n_eff * HPLANCK*CLIGHT / (NM_TO_M * cont%lambda0) / Z / Z / E_RYDBERG - 1.
    
    g_bf = Gaunt_bf(u, n_eff)
    g_bf0 = Gaunt_bf(u0, n_eff)
    
    alpha = cont%alpha0 * (lambda/cont%lambda0)**3  * g_bf / g_bf0
    !alpha0 = 1d-4*2.815d29 * Z**4 * g_bf0 / n**5
 	!alpha =1d-4* 2.815d29 * (Z**4) * g_bf /n_eff**5 * (NM_TO_M*lambda/CLIGHT)**3
write(*,*) lambda, cont%lambda0, g_bf, g_bf0, n_eff
write(*,*) alpha, 1d-4* 2.815d29 * (Z**4) * g_bf /n_eff**5 * (NM_TO_M*lambda/CLIGHT)**3
 RETURN
 END SUBROUTINE test_bf_xs
 
 
 FUNCTION  H_ff_Xsection(Z, T, lambda) result(alpha) !ELEMENTAL 
  real(kind=dp), intent(in) :: lambda, T
  real(kind=dp) :: g_ff, nu3, alpha, K0
  integer, intent(in) :: Z
  
   !K0 = (Q_ELECTRON**2)/(4.0*PI*EPSILON_0) / dsqrt(M_ELECTRON)
   !K0 = (K0**3) * 4./3. * dsqrt(2*pi/3./KBOLTZMANN) / HPLANCK / CLIGHT
   !sigma0_H_ff = K0
  
   nu3 = (NM_TO_M * lambda / CLIGHT)**3 ! = 1 / nu**3
   
   g_ff = Gaunt_ff(lambda, Z, T)

   alpha =  sigma0_H_ff * real(Z) * real(Z) * nu3 * g_ff / dsqrt(T)

   !write(*,*) "alphaff", alpha, sigma0_H_ff, g_ff, T, nu3

 RETURN
 END FUNCTION H_ff_Xsection



 SUBROUTINE Hydrogen_ff(icell, chi)
 ! Hubeny & Mihalas eq. 7.100 (from cgs to SI)
 ! takes place at LTE because it is collisional
  integer, intent(in) :: icell
  integer :: la
  real(kind=dp), dimension(NLTEspec%Nwaves), intent(out) :: chi
  real(kind=dp) :: stim, np, arg_exp, exp_val, C0, alpha

  np = Hydrogen%n(Hydrogen%Nlevel,icell) !nH+=Nion H
  !should be nstar instead of %n ?
  
  chi(:) = 0d0

  if (atmos%ne(icell) == 0d0) return
  

  do la=1,NLTEspec%Nwaves
   stim = 1. - dexp(-hc_k/NLTEspec%lambda(la)/atmos%T(icell))
   
   ! = alpha0 /nu**3 / sqrt(T) = m^5
   !I now consider ne as part of the cross-sections to be in the same units as
   !the bound-free cross-sections
   alpha = H_ff_Xsection(1, atmos%T(icell), NLTEspec%lambda(la)) * atmos%ne(icell)
   
   chi(la) =  alpha * np * stim
   !chi(la) = alpha * 1.38985d21 * 1.37901d21 * stim / atmos%ne(icell)
   !write(*,*) NLTEspec%lambda(la), chi(la), atmos%ne(icell), np, stim

  enddo


 RETURN
 END SUBROUTINE Hydrogen_ff
 
!building
 SUBROUTINE atom_ff_transitions(atom, icell, chi)
 ! Hubeny & Mihalas eq. 7.100 (from cgs to SI)
 ! takes place at LTE because it is collisional
 !should work with Hydrogen
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  real(kind=dp), dimension(NLTEspec%Nwaves), intent(out) :: chi
  real(kind=dp) :: stim, nion, arg_exp, exp_val
  integer :: ic, Z, la

 chi(:) = 0d0
 if (atmos%ne(icell)==0d0) return

 ic = atom%Nlevel 
 Z = atom%stage(ic)

 nion = atom%n(ic, icell)

 do la=1,NLTEspec%Nwaves

    stim = 1.- dexp(-hc_k/NLTEspec%lambda(la)/atmos%T(icell))
   
    chi(la) = H_ff_Xsection(Z, atmos%T(icell), NLTEspec%lambda(la)) * atmos%ne(icell) * nion * stim

 enddo
  

 RETURN
 END SUBROUTINE atom_ff_transitions
 
 SUBROUTINE atom_ff_transitions_old(atom, icell, chi)
 ! Hubeny & Mihalas eq. 7.100 (from cgs to SI)
 ! takes place at LTE because it is collisional
 !should work with Hydrogen
  integer, intent(in) :: icell
  type (AtomType), intent(in) :: atom
  real(kind=dp), dimension(NLTEspec%Nwaves), intent(out) :: chi
  real(kind=dp) :: stim, nion, arg_exp, exp_val
  integer :: ic, Z, la, kr, kc

 chi(:) = 0d0
 if (atmos%ne(icell)==0d0) return

 
  do kr=atom%Ntr_line+1,atom%Ntr
   kc = atom%at(kr)%ik
   ic = atom%continua(kc)%j

   Z = atom%stage(ic)! = cont%atom%stage(cont%i) + 1  

   nion = atom%n(ic, icell)

   do la=1,NLTEspec%Nwaves

    stim = 1.- dexp(-hc_k/NLTEspec%lambda(la)/atmos%T(icell))
   
    chi(la) =  chi(la) + H_ff_Xsection(Z, atmos%T(icell), NLTEspec%lambda(la)) &
     * atmos%ne(icell) * nion * stim

   enddo
  
  enddo !over continua

 RETURN
 END SUBROUTINE atom_ff_transitions_old
 
!  SUBROUTINE Hminus_bf_ff(icell, chi, eta)
!    integer, intent(in) :: icell
!    real(kind=dp) :: lambda, theta, Pe_cgs, fact
!    real(kind=dp), dimension(NLTEspec%Nwaves), intent(inout) :: chi, eta, alpha
!    
!    chi(:) = 0.0_dp
!    eta(:) = 0.0_dp
!    
!    !1e-18 cm2/H-
!    alpha(:) = 1.99654 - 1.18267d-5*(NLTEspec%lambda(:)*10) + &
!               2.64243e-6d-6 * (NLTEspec%lambda * 10)**2 + &
!               -4.40524d-10 * (NLTEspec%lambda*10)**3 + &
!               3.23992d-14 * (NLTEspec%lambda*10)**4 + &
!               -1.39568d-18 * (NLTEspec%lambda*10)**5 + &
!               2.78701d-23 * (NLTEspec%lambda*10)**6
!    !cm2/HI -> m^-1   
!    fact = 1d-4 * sum(Hydrogen%n(1:hydrogen%Nlevel-1,icell))
!      
!    theta = 5040./atmos%T(icell)
!    Pe_cgs = 1d7 * KBOLTZMANN * atmos%ne(icell) * 1d-6 * atmos%T(icell)  
!    !b-f   
!    chi(:) = fact * 4.158d-10 * alpha(:) * Pe_cgs * theta**(2.5) * 10**(0.754*theta)
!    
!    eta(:) = alpha(:) * (twohc/NLTEspec%lambda**3) * 
!  
!  RETURN
!  END SUBROUTINE Hminus_bf_ff
 
!fit
 SUBROUTINE Hminus_bf(icell, chi, eta)
 !-----------------------------------------------------------------
 ! Calculates the negative hydrogen (H-) bound-free continuum absorption coefficient per
 ! hydrogen atom
 !  John 1988 A&A 193, 189
 !		lambda: wavelength greater than 125 nm
 !-----------------------------------------------------------------
   integer, intent(in) :: icell
   real(kind=dp) :: lambda, lambda0, K, alpha, sigma, flambda, Cn(6), diff, stm,funit, cte, Keta
   integer :: la,n
   real(kind=dp) :: arg_exp, nHmin
   real(kind=dp), dimension(NLTEspec%Nwaves), intent(inout) :: chi, eta
   
   chi(:) = 0.0_dp
   eta(:) = 0.0_dp
   
   !1dyne/cm2 = 1e-1 Pa
   !1dyne = 1e-1 Pa cm2
   !cm4/dyne = cm4 / (1e-1 Pa cm2)
   !cm4/dyne = 10 * cm2 / Pa
   !cm4/dyne = 10 * (1e-2)**2 m2/Pa = 1e-3 * m2/Pa
   !m2/Pa = cm4/dyne  * 1e3
   funit = 1d-3 !cm4/dyn -> m2/Pa
   
   Cn(:) = (/152.519d0,49.534d0,-118.858d0,92.536d0,-34.194d0,4.982d0/)

   alpha = hc_k / MICRON_TO_NM ! hc_k = hc/k/nm_to_m

   !alpha = hc_k / micron_to_nm * nm_to_m
   cte = 1d-18 * 0.75 * atmos%T(icell)**(-2.5) * KBOLTZMANN * atmos%T(icell) * atmos%ne(icell)
   lambda0 = 1.6419 !micron, photo detachement threshold
   K = funit * cte * dexp(alpha/lambda0/atmos%T(icell)) * Hydrogen%n(1,icell)
   Keta = 1d-18 * 1d-4 * atmos%ne(icell) * Hydrogen%n(1,icell)

  nHmin = nH_minus(icell)

   do la=1, NLTEspec%Nwaves
   
    lambda = NLTEspec%lambda(la) / MICRON_TO_NM !nm->micron
    if (NLTEspec%lambda(la) > 125. .and. lambda < lambda0) then

     diff = (1d0/lambda - 1d0/lambda0)
    
     stm = 1. - dexp(-alpha/lambda/atmos%T(icell))
        
     flambda = Cn(1) + Cn(2) * diff**(0.5) + Cn(3) * diff + Cn(4) * diff**(1.5) + &
              Cn(5)*diff**(2.) + Cn(6) * diff**(2.5)
    
     sigma = lambda**3d0 * diff**(1.5) * flambda
     !chi(la) = K * stm * sigma
     chi(la) = 1d-18 * 1d-4 * sigma * nHmin * stm !nHmin is ni and nj hydrogen%n(1)
     !eta(la) = Keta * twohc/NLTEspec%lambda(la)**3 * sigma * (1.-stm)
     eta(la) = twohc/NLTEspec%lambda(la)**3 * 1d-18 * 1d-4 * sigma * nHmin * (1-stm)
    endif   !remember: at LTE eta = nj * nistar/njstar exp(-de/kt) = nistar * exp(-de/kT) 
   enddo



 RETURN
 END SUBROUTINE Hminus_bf

 
!fit to the best data
 SUBROUTINE Hminus_ff(icell, chi)
!-----------------------------------------------------------------
! Calculates the negative hydrogen (H-) free-free continuum absorption coefficient per
! hydrogen atom (cm^2)
!  REF: John 1989 A&A 193, 189
!  INPUT:
!		T : temperature in K (1400 < T < 100080)
!		wavelength greater than 188 nm
!-----------------------------------------------------------------
   integer, intent(in) :: icell
   !tab. 3a, lambda > 0.3645micron
   real(kind=dp), dimension(6) :: An, Bn, Cn, Dn, En, Fn, com1
   !tab 3b, lambda > 0.1823 micron and < 0.3645 micron, size of 4 isntead 5, because the last two rows are 0
   real(kind=dp), dimension(4) :: An2, Bn2, Cn2, Dn2, En2, Fn2, com2
   real(kind=dp) :: funit, K, kappa
   integer :: la, n
   real(kind=dp), dimension(NLTEspec%Nwaves), intent(inout) :: chi
   real(kind=dp) :: lambda, theta
	
   chi(:) = 0.0_dp
	
	An = (/0.d0,2483.346d0,-3449.889d0,2200.04d0,-696.271d0,88.283d0/)
	Bn = (/0.d0,285.827d0,-1158.382d0,2427.719d0,-1841.4d0,444.517d0/)
    Cn = (/0.d0,-2054.291d0,8746.523d0,-13651.105d0,8624.97d0,-1863.864d0/)
    Dn = (/0.d0,2827.776d0,-11485.632d0,16755.524d0,-10051.53d0,2095.288d0/)
    En = (/0.d0,-1341.537d0,5303.609d0,-7510.494d0,4400.067d0,-901.788d0/)
    Fn = (/0.d0,208.952d0,-812.939d0,1132.738d0,-655.02d0,132.985d0/)
    An2 = (/518.1021d0,473.2636d0,-482.2089d0,115.5291d0/)
    Bn2 = (/-734.8666d0,1443.4137d0,-737.1616d0,169.6374d0/)
    Cn2 = (/1021.1775d0,-1977.3395d0,1096.8827d0,-245.649d0/)
    Dn2 = (/-479.0721d0,922.3575d0,-521.1341d0,114.243d0/)
    En2 = (/93.1373d0,-178.9275d0,101.7963d0,-21.9972d0/)
    Fn2 = (/-6.4285d0,12.36d0,-7.0571d0,1.5097d0/)
    
    funit = 1d-3 !cm4/dyne to m2/Pa
    K = 1d-29 * funit * atmos%ne(icell) * KBOLTZMANN * atmos%T(icell) * Hydrogen%n(1,icell)
    
    if (atmos%ne(icell)==0d0) return
      
    do la=1, NLTEspec%Nwaves
      lambda = NLTEspec%lambda(la) / MICRON_TO_NM
      if (lambda < 0.1823) cycle
      
      theta = 5040d0 / atmos%T(icell) 
	  if (lambda < 0.3645) then
	    kappa = 0d0
	    do n=1,4
	     kappa = kappa + theta**((n+1)/2d0) * (lambda**2 * An2(n) + Bn2(n) + Cn2(n)/lambda + &
	     	Dn2(n)/lambda**2 + En2(n)/lambda**3 + Fn2(n)/lambda**4)
	    enddo
	  else
	    kappa = 0D0
	    do n=1,6
	     kappa = kappa + theta**((n+1)/2d0) * (lambda**2 * An(n) + Bn(n) + Cn(n)/lambda + &
	     	Dn(n)/lambda**2 + En(n)/lambda**3 + Fn(n)/lambda**4)
	    enddo
	  endif
	  chi(la) = kappa * K
	enddo
 
 RETURN
 END SUBROUTINE Hminus_ff

!Interpolation to the best data
!  SUBROUTINE Hminus_ff_longwavelength(icell, chi)
!  ! H- free-free opacity for wavelength beyond 9113.0 nm
!  ! see: T. L. John (1988), A&A 193, 189-192 (see table 3a).
!  ! His results are based on calculations by K. L. Bell and
!  ! K. A. Berrington (1987), J. Phys. B 20, 801-806.
!   integer :: k, n
!   integer, intent(in) :: icell
!   real(kind=dp), intent(inout), dimension(NLTEspec%Nwaves) :: chi
!   real(kind=dp), dimension(NJOHN) :: AJ, BJ, CJ, DJ, EJ, FJ
!   real(kind=dp), dimension(NJOHN,NLTEspec%Nwaves) :: Clam
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: lambda_mu, lambda_inv
!   real(kind=dp) :: sqrt_theta, theta_n, CK
! 
!   data AJ / 0.000,  2483.346, -3449.889,  2200.040, -696.271, 88.283   /
!   data BJ / 0.000,   285.827, -1158.382,  2427.719,-1841.400, 444.517  /
!   data CJ / 0.000, -2054.291,  8746.523,-13651.105,8624.970, -1863.864 /
!   data DJ / 0.000,2827.776,-11485.632, 16755.524,-10051.530, 2095.288  /
!   data EJ / 0.000, -1341.537,  5303.609, -7510.494,4400.067,  -901.788 /
!   data FJ / 0.000,   208.952,  -812.939,  1132.738, -655.020, 132.985  /
! 
!   chi = 0d0
! 
!   CK= (KBOLTZMANN * THETA0 * 1.0E-32);
! 
!   lambda_mu = NLTEspec%lambda / MICRON_TO_NM
!   lambda_inv = 1. / lambda_mu
! 
!   do n=1,NJOHN
!    Clam(n,:) = (lambda_mu)**2 * AJ(n) + BJ(n) + lambda_inv * &
!     (CJ(n) + lambda_inv*(DJ(n) + lambda_inv*(EJ(n) + &
!       lambda_inv*FJ(n))))
!   end do
! 
!   theta_n = 1.
!   sqrt_theta = dsqrt(THETA0/atmos%T(icell))
!   do n=2,NJOHN
!     theta_n = theta_n * sqrt_theta
!     chi = chi + theta_n * Clam(n,:)
!   end do
!   chi= chi* Hydrogen%n(1,icell) * (atmos%ne(icell)*CK)
! 
!  RETURN
!  END SUBROUTINE Hminus_ff_longwavelength
! 
!  SUBROUTINE Hminus_ff(icell, chi)
!  ! from RH
!  ! Hminus free-free coefficients in 1d-29 m^5/J
!  ! see Stilley and Callaway 1970 ApJ 160
!   integer, intent(in) :: icell
!   real(kind=dp), intent(out), dimension(NLTEspec%Nwaves) :: chi
!   integer :: k, index, index2
!   real(kind=dp) :: lambdaFF(NFF), thetaFF(NTHETA)
!   real(kind=dp), dimension(NFF*NTHETA) :: kappaFF_flat
!   real(kind=dp) :: theta(1), pe, kappa(1, NLTEspec%Nwaves)
!   real(kind=dp), dimension(NTHETA,NFF) :: kappaFF
!   data lambdaFF / 0.0, 303.8, 455.6, 506.3, 569.5, 650.9, &
!                   759.4, 911.3, 1013.0, 1139.0, 1302.0,   &
!                   1519.0, 1823.0, 2278.0, 3038.0, 4556.0, &
!                   9113.0 /
! 
!   !theta = 5040. K / T
!   data thetaFF / 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,&
!                  1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 /
! 
!   data kappaFF_flat / 0.00d0, 0.00d0, 0.00d0, 0.00d0,&! 0nm
!                  0.00d0, 0.00d0, 0.00d0, 0.00d0,&
!                  0.00d0, 0.00d0, 0.00d0, 0.00d0,&
!                  0.00d0, 0.00d0, 0.00d0, 0.00d0,&
! 3.44d-2, 4.18d-2, 4.91d-2, 5.65d-2, 6.39d-2,&!303.8nm
! 7.13d-2,7.87d-2, 8.62d-2, 9.36d-2, 1.01d-1, &
! 1.08d-1, 1.16d-1, 1.23d-1, 1.30d-1, 1.38d-1, 1.45d-1,&
! 7.80d-2, 9.41d-2, 1.10d-1, 1.25d-1, 1.40d-1,&!455.6nm
! 1.56d-1,1.71d-1, 1.86d-1, 2.01d-1, 2.16d-1, &
! 2.31d-1, 2.45d-1,2.60d-1, 2.75d-1, 2.89d-1, 3.03d-1,&
! 9.59d-2, 1.16d-1, 1.35d-1, 1.53d-1, 1.72d-1,&!506.3 nm
! 1.90d-1,2.08d-1, 2.25d-1, 2.43d-1, 2.61d-1, &
! 2.78d-1, 2.96d-1,3.13d-1, 3.30d-1, 3.47d-1, 3.64d-1,&
! 1.21d-1, 1.45d-1, 1.69d-1, 1.92d-1, 2.14d-1,&!569.5 nm
! 2.36d-1,2.58d-1, 2.80d-1, 3.01d-1, 3.22d-1, &
! 3.43d-1, 3.64d-1,3.85d-1, 4.06d-1, 4.26d-1, 4.46d-1,&
! 1.56d-1, 1.88d-1, 2.18d-1, 2.47d-1, 2.76d-1,&!650.9 nm
! 3.03d-1,3.31d-1, 3.57d-1, 3.84d-1, 4.10d-1, &
! 4.36d-1, 4.62d-1,4.87d-1, 5.12d-1, 5.37d-1, 5.62d-1,&
! 2.10d-1, 2.53d-1, 2.93d-1, 3.32d-1, 3.69d-1,&!759.4 nm
! 4.06d-1, 4.41d-1, 4.75d-1, 5.09d-1, 5.43d-1, &
! 5.76d-1, 6.08d-1, 6.40d-1, 6.72d-1, 7.03d-1, 7.34d-1,&
! 2.98d-1, 3.59d-1, 4.16d-1, 4.70d-1, 5.22d-1,&!911.3 nm
! 5.73d-1,6.21d-1, 6.68d-1, 7.15d-1, 7.60d-1, 8.04d-1,&
! 8.47d-1,8.90d-1, 9.32d-1, 9.73d-1, 1.01d0,&
! 3.65d-1, 4.39d-1, 5.09d-1, 5.75d-1, 6.39d-1,&!1013 nm
! 7.00d-1, 7.58d-1, 8.15d-1, 8.71d-1, 9.25d-1, &
!     9.77d-1, 1.03d0, 1.08d0, 1.13d0, 1.18d0, 1.23d0,&
!     4.58d-1, 5.50d-1, 6.37d-1, 7.21d-1, 8.00d-1,&!1139 nm
!     8.76d-1,9.49d-1, 1.02d0, 1.09d0, 1.15d0, 1.22d0,&
!     1.28d0,1.34d0, 1.40d0, 1.46d0, 1.52d0,&
!     5.92d-1, 7.11d-1, 8.24d-1, 9.31d-1, 1.03d0,&!1302 nm
!     1.13d0, 1.23d0, 1.32d0, 1.40d0, 1.49d0,&
!     1.57d0, 1.65d0, 1.73d0, 1.80d0, 1.88d0, 1.95d0,&
!     7.98d-1, 9.58d-1, 1.11d0, 1.25d0, 1.39d0,&!1519 nm
!     1.52d0, 1.65d0, 1.77d0, 1.89d0, 2.00d0, &
!     2.11d0, 2.21d0, 2.32d0, 2.42d0, 2.51d0, 2.61d0,&
!     1.14d0, 1.36d0, 1.58d0, 1.78d0, 1.98d0,&!1823 nm
!     2.17d0, 2.34d0, 2.52d0, 2.68d0, 2.84d0,&
!     3.00d0, 3.15d0, 3.29d0, 3.43d0, 3.57d0, 3.70d0,&
!     1.77d0, 2.11d0, 2.44d0, 2.75d0, 3.05d0,&!2278 nm
!     3.34d0, 3.62d0, 3.89d0, 4.14d0, 4.39d0,&
!     4.63d0, 4.86d0, 5.08d0, 5.30d0, 5.51d0, &
!     5.71d0,3.10d0, 3.71d0, 4.29d0, 4.84d0, 5.37d0,&!3038 nm
!     5.87d0, 6.36d0, 6.83d0, 7.28d0, 7.72d0, &
!     8.14d0, 8.55d0,8.95d0, 9.33d0, 9.71d0, 1.01d1,&
!     6.92d0, 8.27d0, 9.56d0, 1.08d1, 1.19d1,&!4556 nm
!     1.31d1,1.42d1, 1.52d1, 1.62d1, 1.72d1, 1.82d1,&
!     1.91d1,2.00d1, 2.09d1, 2.17d1, 2.25d1,&
!     2.75d1, 3.29d1, 3.80d1, 4.28d1, 4.75d1,&!9113 nm
!     5.19d1,5.62d1, 6.04d1, 6.45d1, 6.84d1, 7.23d1,&
!     7.60d1,7.97d1, 8.32d1, 8.67d1, 9.01d1 /
! 
!   chi = 0d0
! 
!   pe = atmos%ne(icell) * KBOLTZMANN * atmos%T(icell)
!   ! 2-dimensionalize kappaFF_flat for interp2D
!   kappaFF = RESHAPE(kappaFF_flat,(/NTHETA, NFF/))
!   !!write(*,*) "Check reshape + interp2D"
!   !do index=1,NFF
!   ! do index2=1,NTHETA
!   !  kappaFF(index,index2) = &
!   !     kappaFF_flat(NTHETA*(index-1)+index2)
!   ! end do
!   !end do
!   
!   ! for long wavelengths
!   ! Can do more efficiently !! that computing for all wavelength
!   ! as it was all long wavelength and then computing
!   ! only where it is below for the other wavelengths
!   if ((MAXVAL(NLTEspec%lambda) >= lambdaFF(NFF)) .or. &
!        (MINVAL(NLTEspec%lambda) >= lambdaFF(NFF))) then !if at least one
!    CALL Hminus_ff_longwavelength(icell, chi)
!    !!RETURN
!   end if
! 
!    theta(1:1) = THETA0 /  atmos%T(icell)
!    ! interpolate kappaFF at theta and lambda using
!    ! 2x 1D cubic Bezier splines
!    kappa = interp2Darr(NTHETA, thetaFF,NFF,lambdaFF,kappaFF,&
!                   1,theta,NLTEspec%Nwaves, NLTEspec%lambda)
!    where(NLTEspec%lambda < lambdaFF(NFF))
!    chi = (Hydrogen%n(1,icell)*1d-29) * pe * kappa(1,:)
!    end where
! 
! 
!  RETURN
!  END SUBROUTINE Hminus_ff
!  SUBROUTINE Hminus_bf(icell, chi, eta)
!  ! H minus bound-free coefficients from RH
!  ! in units of 1e-21 m^2
!  ! See: Geltman 1962, ApJ 136, 935-945 and
!  ! Mihalas Stellar atmospheres
!   logical :: res
!   integer :: icell
!   real(kind=dp), intent(out), dimension(NLTEspec%Nwaves) :: &
!    chi, eta
!   real(kind=dp), dimension(NBF) :: lambdaBF, alphaBF
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: hc_kla, stimEmis, twohnu3_c2, alpha
! 
!   data lambdaBF / 0.0, 50.0, 100.0, 150.0, 200.0, 250.0,  &
!                  300.0, 350.0, 400.0, 450.0, 500.0, 550.0,&
!                  600.0, 650.0, 700.0, 750.0, 800.0, 850.0,&
!                  900.0, 950.0, 1000.0, 1050.0, 1100.0,    &
!                  1150.0, 1200.0, 1250.0, 1300.0, 1350.0,  &
!                  1400.0, 1450.0, 1500.0, 1550.0, 1600.0,  &
!                  1641.9 /
! 
!   data alphaBF / 0.0,  0.15, 0.33, 0.57, 0.85, 1.17, 1.52,&
!                  1.89, 2.23, 2.55, 2.84, 3.11, 3.35, 3.56,&
!                  3.71, 3.83, 3.92, 3.95, 3.93, 3.85, 3.73,&
!                  3.58, 3.38, 3.14, 2.85, 2.54, 2.20, 1.83,&
!                  1.46, 1.06, 0.71, 0.40, 0.17, 0.0 /
! 
!   chi = 0d0
!   eta = 0d0
! 
! 
!   ! interpolate cross-section at lambda
!   !!alpha = 1e-21 * interp1D(lambdaBF,alphaBF,lambda) !m^2
!   CALL bezier3_interp(NBF, lambdaBF, alphaBF, NLTEspec%Nwaves, NLTEspec%lambda, alpha)
!   alpha = alpha * 1d-21
!   
!   ! hnu/kB
!   hc_kla = (HPLANCK*CLIGHT) / (KBOLTZMANN*NLTEspec%lambda*NM_TO_M)
!   twohnu3_c2 = (2.*HPLANCK*CLIGHT) / (NM_TO_M*NLTEspec%lambda)**3
!   stimEmis = dexp(-hc_kla/atmos%T(icell))
!   
!   !non zero only in this region
!   where((NLTEspec%lambda > lambdaBF(1)).and.(NLTEspec%lambda < lambdaBF(NBF)))
!    ! note that nhmin is deduced from hydrogen nTotal = A*nHtot
!    ! which is the total number of hydrogen in neutral and
!    ! ionised form (H and H+).
!    ! nHtot is the total number of Hydrogen present in
!    ! all form = HI, HII, H-, H2, CH etc etc
!    chi = nH_minus(icell) * (1.-stimEmis) * alpha
!    eta = nH_minus(icell) * twohnu3_c2 * stimEmis * alpha
!   end where
! 
!  RETURN
!  END SUBROUTINE Hminus_bf
!  

END MODULE hydrogen_opacities
