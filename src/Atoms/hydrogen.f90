! Hydrogen, bound-free, free-free and continuum opacities
! - H b-f and f-f
! - H^- b-f and f-f
! - H2^- f-f  --> Not yet
! - H2^+ f-f  --> Not yet
! - Rayleigh scattering by H2 --> Not yet
! - Rayleigh scattering by H
!
!
! Note:
! to include molecular hydrogen opacities,
! chemical equilibrium as well as molecules-related
! routines have to be created. TBD
!
! Note 2:
!  free-free emissivity is obtained by
! -> eta_ff = chi_ff * Bplanck
!
!
! Opacities in m^-1 (chi)
! Emissivities in J/s/m3/Hz/sr (eta)
! Note:
!  Opacitiy and Emissivity arrays are allocated outside the
!  module.
! Source function S = eta/chi in J/s/m3/Hz/sr / m^-1=W/m2/Hz/sr
!
! Note 2: chi and eta given for the all wavelength grid at the cell icell
MODULE hydrogen_opacities

 use atom_type, only : AtomicContinuum
 use atmos_type, only : atmos, Hydrogen, Helium
 use constant
 use spectrum_type, only : NLTEspec
 use math, only : bezier3_interp, interp2Darr

 IMPLICIT NONE

 integer, parameter :: NBF=34, NJOHN=6, NFF=17, NTHETA=16

 CONTAINS

 FUNCTION Gaunt_bf(Nl,u, n_eff) result(GII)
  ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
  ! See also Menzel & Pekeris 1935
  integer :: Nl
  double precision :: n_eff, GII(Nl)
  double precision :: u(Nl), n23

  n23 = n_eff**(-6.6666666666667d-1) !n^-2/3

  GII = 1.+0.1728 * n23 * ((u+1)**(-6.6666666666667d-1)) * &
      (u-1) - 0.0496*n23*n23 * (u*u + 4./3. * u + 1) !+...

  if (MINVAL(GII) < 0d0) then
   !write(*,*) "Warning, Gaunt factor is less than 0"
   !write(*,*) Nl, minval(u), maxval(u), n_eff, minval(GII)
   GII = dabs(GII)
  end if

  if (MAXVAL(GII) > 2d0) then
    write(*,*) "Bound-free Gaunt factor gII = ", MAXVAL(GII)
  end if

 RETURN
 END FUNCTION Gaunt_bf

 FUNCTION Gaunt_ff(lambda, charge, T) result(GIII)
 ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
 !
 ! Note: There is a problem with this expansion at higher temperatures
 ! (T > 3.0E4 and longer wavelengths (lambda > 2000 nm). Set to
 ! 1.0 when the value goes below 1.0
  double precision, intent(in) :: lambda(NLTEspec%Nwaves)
  double precision, dimension(NLTEspec%Nwaves) :: x, x3, y, GIII
  double precision :: charge, T

  x = ((HPLANCK * CLIGHT)/(lambda * NM_TO_M)) / &
       (E_RYDBERG * (charge)**(2d0))
  x3 = (x**(3.3333333d-1))
  y  = (2.0 * lambda * NM_TO_M * KBOLTZMANN*T) / &
       (HPLANCK*CLIGHT)

  gIII = 1.0 + 0.1728*x3 * (1.0 + y) - &
        0.0496*(x3*x3) * (1.0 + (1.0 + y)*0.33333333*y)

  where (GIII <= 1d0) gIII = 1.

  if (MAXVAL(GIII) > 2d0) write(*,*) "free-free Gaunt factor gIII = ", &
      MAXVAL(gIII)

 RETURN
 END FUNCTION Gaunt_ff


 SUBROUTINE Hydrogen_bf (icell, chi, eta)
  ! Computes LTE hydrogen bound-free opacity
  ! See Hubeny & Mihalas 2014, chap. 7
  double precision, dimension(NLTEspec%Nwaves) :: uu
  integer, intent(in) :: icell
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: &
      chi, eta
  type (AtomicContinuum) :: continuum
  integer :: i, kr, k, nc
!  integer, dimension(:), allocatable :: iLam
  double precision :: lambdaEdge, sigma(NLTEspec%Nwaves), sigma0, g_bf(NLTEspec%Nwaves), &
    twohnu3_c2(NLTEspec%Nwaves), twohc, gijk(NLTEspec%Nwaves), hc_k, hc_kla(NLTEspec%Nwaves), &
    expla(NLTEspec%Nwaves), n_eff,npstar, sigma02, sigma2(NLTEspec%Nwaves), np


  ! initialize for this cell point
   chi = 0.
   eta = 0.
   uu = 0d0
   g_bf = 0d0
   sigma = 0d0

  twohc = (2.*HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
  hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
  sigma0 = (32.)/(PI*3.*dsqrt(3d0)) * EPSILON_0 * &
          (HPLANCK**(3d0)) / (CLIGHT * &
          (M_ELECTRON*Q_ELECTRON)**(2d0))



  ! the facotre 1/NM_TO_M per lambda is in hc_k and twohc
  hc_kla = hc_k / NLTEspec%lambda
  twohnu3_c2 = twohc / (NLTEspec%lambda)**3 !J/s/m2/Hz/sr
  ! because Planck'law (nu) = c*U(nu) / 4pi = m/s/sr * &
  !     unit(spectral energy density)

         ! = 7.907d-22 ! m^2
  !! Or use this
  !sigma02 = (1.)/(PI*48.*sqrt(3d0)) * M_ELECTRON * &
  !     dpow(Q_ELECTRON,10d0) / (CLIGHT * dpow(HPLANCK,6d0) * &
  !     dpow(EPSILON_0,5d0)) !2.8154d25 ! m^2 * Hz^3


  ! LTE number of protons (H+)
  npstar = Hydrogen%nstar(Hydrogen%Nlevel,icell)
  do kr=1,Hydrogen%Ncont
   continuum = Hydrogen%continua(kr)
   lambdaEdge = continuum%lambda0 !ionisation frequency (min)
   i = continuum%i
    ! evaluate effective principal quantum number ->
    ! n_eff = dsqrt(E_RYDBERG / &
    !   (Hydrogen%E(continuum%j)-Hydrogen%E(i))) !Z=1

    ! Principal quantum number n from statistical weight
    ! of the continuum level
    n_eff = dsqrt(Hydrogen%g(i)/2.)  !only for Hydrogen !
    np = Hydrogen%n(Hydrogen%Nlevel,icell)

!    if ((Hydrogen%n(i,icell) <= 0.).or.(npstar <= 0.)) CYCLE
    ! -> prevents dividing by zero
    if ((Hydrogen%n(i,icell) <= 0.).or.(npstar <= 0.)) then
       write(*,*) "(Hydrogen_bf) Warning at icell=", icell," T(K)=", atmos%T(icell)
       if (npstar <= 0) then
          write(*,*) "np density <= 0"
          write(*,*) "skipping this level"
         CYCLE
       else
          write(*,*) "Hydrogen%n(i) density <= 0 for i=", i
          write(*,*) "skipping this level"
         CYCLE
       end if
     end if


!     allocate(iLam(continuum%Nlambda)) !not used yet
!     iLam = (/ (nc, nc=continuum%Nblue, continuum%Nblue+continuum%Nlambda-1) /)

   ! lambda has to be lower than the edge but greater than
   ! lambda0, see Hubeny & Mihalas chap. 7 photoionisation
   ! and Rutten hydrogen bound-free parts.
   ! lambda0 or lambdaEdge is the maximum wavelength (
   ! minimal frequency) required for the photon to be unbound
   ! i.e., it is the wavelength (or frequency) of photoionisation.
   ! further, if lambda is lower than lambda min (or freq max)
   ! defined in the atomic file, it is not defined.
   ! Below the Edge (or beyond is frequency), the photoioniz
   ! cross-section alpha falls off in lambda^3 (or nu^-3).


    ! u = n**2 * h * nu / (Z**2 * R) - 1
!     uu(ilam) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*NLTEspec%lambda(ilam)) / &
!       ((Hydrogen%stage(i)+1)* (Hydrogen%stage(i)+1)) / E_RYDBERG - 1.
    uu(continuum%Nblue:continuum%Nred) = &
      n_eff*n_eff*HPLANCK*CLIGHT/ & 
       (NM_TO_M*NLTEspec%lambda(continuum%Nblue:continuum%Nred)) / &
      ((Hydrogen%stage(i)+1)*(Hydrogen%stage(i)+1)) / E_RYDBERG - 1.
!    g_bf(ilam) = Gaunt_bf(continuum%Nlambda, uu(ilam), n_eff)
    g_bf(continuum%Nblue:continuum%Nred) = &
     Gaunt_bf(continuum%Nlambda, uu(continuum%Nblue:continuum%Nred), n_eff)
    ! Z scaled law
!     sigma(ilam) = &
!      sigma0 * g_bf(ilam) * (NLTEspec%lambda(ilam)/lambdaEdge)**3 * n_eff / (Hydrogen%stage(i)+1)**2 !m^2
    sigma(continuum%Nblue:continuum%Nred) = &
     sigma0 * g_bf(continuum%Nblue:continuum%Nred) * &
       (NLTEspec%lambda(continuum%Nblue:continuum%Nred)/lambdaEdge)**3 * &
       n_eff / (Hydrogen%stage(i)+1)**2 !m^2


    !! or use this (with sigma02)
    !sigma2 = &
    ! sigma02 * g_bf * &
    ! dpow(real(Hydrogen%stage(i)+1,kind=8),4d0) / &
    !   (dpow(n_eff,5d0) * &
    !     CUBEarr(NLTEspec%lambda*NM_TO_M/(HPLANCK*CLIGHT)))


    ! now computes emissivity and extinction
!     expla(iLam) = dexp(-hc_kla(iLam)/atmos%T(icell))
    expla(continuum%Nblue:continuum%Nred) = dexp(-hc_kla(continuum%Nblue:continuum%Nred)/atmos%T(icell))
!     gijk(iLam) = Hydrogen%nstar(i,icell)/npstar * expla(iLam)
    gijk(continuum%Nblue:continuum%Nred) = Hydrogen%nstar(i,icell)/npstar * &
      expla(continuum%Nblue:continuum%Nred)
     ! at LTE only for chi
     ! see Hubeny & Mihalas eq. 14.16 to 14.18
     ! if LTE Hydrogen%n points to %nstar
!    chi(iLam) = chi(iLam) + sigma(iLam) * (1.-expla(iLam)) * Hydrogen%n(i,icell)
    chi(continuum%Nblue:continuum%Nred) = chi(continuum%Nblue:continuum%Nred) + &
      sigma(continuum%Nblue:continuum%Nred) * (1.-expla(continuum%Nblue:continuum%Nred)) * Hydrogen%n(i,icell)
!   eta(iLam) = eta(iLam) + twohnu3_c2(iLam) * gijk(iLam) * sigma(iLam) * np       
    eta(continuum%Nblue:continuum%Nred) = eta(continuum%Nblue:continuum%Nred) + &
      twohnu3_c2(continuum%Nblue:continuum%Nred) * gijk(continuum%Nblue:continuum%Nred) * &
        sigma(continuum%Nblue:continuum%Nred) * np

!    deallocate(iLam)
  end do

 RETURN
 END SUBROUTINE Hydrogen_bf

 SUBROUTINE Hydrogen_ff(icell, chi)
 ! Hubeny & Mihalas eq. 7.100 (from cgs to SI)
 ! takes place at LTE because it is collisional
  integer, intent(in) :: icell
  double precision, dimension(NLTEspec%Nwaves), intent(out) :: chi
  double precision :: gff(NLTEspec%Nwaves), hc_kla(NLTEspec%Nwaves), &
     stim(NLTEspec%Nwaves), nu3(NLTEspec%Nwaves), np, sigma0, charge, e0cgs, KBcgs

  charge = 1d0
  ! CGS units constants
  e0cgs = 4.80320427/(1d10) !Fr (ESU, Gaussian)
  KBcgs = 1.3806504/(1d16) !erg/K
  sigma0 = dsqrt(32.*PI) / (3.*dsqrt(3d0)) * &
    ((e0cgs)**(6d0)) / (1d2*CLIGHT*1d7 * HPLANCK*dsqrt(KBcgs*&
      (1d3*M_ELECTRON)**(3d0))) !cm^5 K^1/2 Hz^3
  ! = 3.6923284d8 ! cm^5 K^1/2 Hz^3
  nu3 = (NLTEspec%lambda*NM_TO_M/CLIGHT)**3 !inverse of nu3=1/nu3
  hc_kla = (HPLANCK*CLIGHT) / (KBOLTZMANN*NM_TO_M*NLTEspec%lambda)

  np = Hydrogen%n(Hydrogen%Nlevel,icell) !nH+=Nion H
  stim = 1.-dexp(-hc_kla/atmos%T(icell))

  gff = Gaunt_ff(NLTEspec%lambda, charge, atmos%T(icell))

   !write(*,*) gff
   !1e-10 = cm5 ->m5
  chi = 1d-10 * gff * charge*charge*sigma0 / &
    dsqrt(atmos%T(icell)) * nu3 * atmos%ne(icell) * np * stim
   !write(*,*) chi(k), 1d-10 * sigma0, nu3, atmos%ne(k), np(k), stim, gff

 RETURN
 END SUBROUTINE Hydrogen_ff

 SUBROUTINE Hminus_bf(icell, chi, eta)
 ! H minus bound-free coefficients from RH
 ! in units of 1e-21 m^2
 ! See: Geltman 1962, ApJ 136, 935-945 and
 ! Mihalas Stellar atmospheres
  logical :: res
  integer :: icell
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: &
   chi, eta
  double precision, dimension(NBF) :: lambdaBF, alphaBF
  double precision, dimension(NLTEspec%Nwaves) :: hc_kla, stimEmis, twohnu3_c2, alpha

  data lambdaBF / 0.0, 50.0, 100.0, 150.0, 200.0, 250.0,  &
                 300.0, 350.0, 400.0, 450.0, 500.0, 550.0,&
                 600.0, 650.0, 700.0, 750.0, 800.0, 850.0,&
                 900.0, 950.0, 1000.0, 1050.0, 1100.0,    &
                 1150.0, 1200.0, 1250.0, 1300.0, 1350.0,  &
                 1400.0, 1450.0, 1500.0, 1550.0, 1600.0,  &
                 1641.9 /

  data alphaBF / 0.0,  0.15, 0.33, 0.57, 0.85, 1.17, 1.52,&
                 1.89, 2.23, 2.55, 2.84, 3.11, 3.35, 3.56,&
                 3.71, 3.83, 3.92, 3.95, 3.93, 3.85, 3.73,&
                 3.58, 3.38, 3.14, 2.85, 2.54, 2.20, 1.83,&
                 1.46, 1.06, 0.71, 0.40, 0.17, 0.0 /

  chi = 0d0
  eta = 0d0


  ! interpolate cross-section at lambda
  !!alpha = 1e-21 * interp1D(lambdaBF,alphaBF,lambda) !m^2
  CALL bezier3_interp(NBF, lambdaBF, alphaBF, NLTEspec%Nwaves, NLTEspec%lambda, alpha)
  alpha = alpha * 1d-21
  
  ! hnu/kB
  hc_kla = (HPLANCK*CLIGHT) / (KBOLTZMANN*NLTEspec%lambda*NM_TO_M)
  twohnu3_c2 = (2.*HPLANCK*CLIGHT) / (NM_TO_M*NLTEspec%lambda)**3
  stimEmis = dexp(-hc_kla/atmos%T(icell))
  
  !non zero only in this region
  where((NLTEspec%lambda > lambdaBF(1)).and.(NLTEspec%lambda < lambdaBF(NBF)))
   ! note that nhmin is deduced from hydrogen%ntotal
   ! which is the total number of hydrogen in neutral and
   ! ionised form (H and H+).
   ! nHtot is the total number of Hydrogen present in
   ! all form = HI, HII, H-, H2, CH etc etc
   chi = atmos%nHmin(icell) * (1.-stimEmis) * alpha
   eta = atmos%nHmin(icell) * twohnu3_c2 * stimEmis * alpha
  end where

 RETURN
 END SUBROUTINE Hminus_bf

 SUBROUTINE Hminus_ff_longwavelength(icell, chi)
 ! H- free-free opacity for wavelength beyond 9113.0 nm
 ! see: T. L. John (1988), A&A 193, 189-192 (see table 3a).
 ! His results are based on calculations by K. L. Bell and
 ! K. A. Berrington (1987), J. Phys. B 20, 801-806.
  integer :: k, n
  integer, intent(in) :: icell
  double precision, intent(inout), dimension(NLTEspec%Nwaves) :: chi
  double precision, dimension(NJOHN) :: AJ, BJ, CJ, DJ, EJ, FJ
  double precision, dimension(NJOHN,NLTEspec%Nwaves) :: Clam
  double precision, dimension(NLTEspec%Nwaves) :: lambda_mu, lambda_inv
  double precision :: sqrt_theta, theta_n, CK

  data AJ / 0.000,  2483.346, -3449.889,  2200.040, -696.271, 88.283   /
  data BJ / 0.000,   285.827, -1158.382,  2427.719,-1841.400, 444.517  /
  data CJ / 0.000, -2054.291,  8746.523,-13651.105,8624.970, -1863.864 /
  data DJ / 0.000,2827.776,-11485.632, 16755.524,-10051.530, 2095.288  /
  data EJ / 0.000, -1341.537,  5303.609, -7510.494,4400.067,  -901.788 /
  data FJ / 0.000,   208.952,  -812.939,  1132.738, -655.020, 132.985  /

  chi = 0d0

  CK= (KBOLTZMANN * THETA0 * 1.0E-32);

  lambda_mu = NLTEspec%lambda / MICRON_TO_NM
  lambda_inv = 1. / lambda_mu

  do n=1,NJOHN
   Clam(n,:) = (lambda_mu)**2 * AJ(n) + BJ(n) + lambda_inv * &
    (CJ(n) + lambda_inv*(DJ(n) + lambda_inv*(EJ(n) + &
      lambda_inv*FJ(n))))
  end do

  theta_n = 1.
  sqrt_theta = dsqrt(THETA0/atmos%T(icell))
  do n=2,NJOHN
    theta_n = theta_n * sqrt_theta
    chi = chi + theta_n * Clam(n,:)
  end do
  chi= chi* Hydrogen%n(1,icell) * (atmos%ne(icell)*CK)

 RETURN
 END SUBROUTINE Hminus_ff_longwavelength

 SUBROUTINE Hminus_ff(icell, chi)
 ! from RH
 ! Hminus free-free coefficients in 1d-29 m^5/J
 ! see Stilley and Callaway 1970 ApJ 160
  integer, intent(in) :: icell
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: chi
  integer :: k, index, index2
  double precision :: lambdaFF(NFF), thetaFF(NTHETA)
  double precision, dimension(NFF*NTHETA) :: kappaFF_flat
  double precision :: theta(1), pe, kappa(1, NLTEspec%Nwaves)
  double precision, dimension(NTHETA,NFF) :: kappaFF
  data lambdaFF / 0.0, 303.8, 455.6, 506.3, 569.5, 650.9, &
                  759.4, 911.3, 1013.0, 1139.0, 1302.0,   &
                  1519.0, 1823.0, 2278.0, 3038.0, 4556.0, &
                  9113.0 /

  !theta = 5040. K / T
  data thetaFF / 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,&
                 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 /

  data kappaFF_flat / 0.00d0, 0.00d0, 0.00d0, 0.00d0,&! 0nm
                 0.00d0, 0.00d0, 0.00d0, 0.00d0,&
                 0.00d0, 0.00d0, 0.00d0, 0.00d0,&
                 0.00d0, 0.00d0, 0.00d0, 0.00d0,&
3.44d-2, 4.18d-2, 4.91d-2, 5.65d-2, 6.39d-2,&!303.8nm
7.13d-2,7.87d-2, 8.62d-2, 9.36d-2, 1.01d-1, &
1.08d-1, 1.16d-1, 1.23d-1, 1.30d-1, 1.38d-1, 1.45d-1,&
7.80d-2, 9.41d-2, 1.10d-1, 1.25d-1, 1.40d-1,&!455.6nm
1.56d-1,1.71d-1, 1.86d-1, 2.01d-1, 2.16d-1, &
2.31d-1, 2.45d-1,2.60d-1, 2.75d-1, 2.89d-1, 3.03d-1,&
9.59d-2, 1.16d-1, 1.35d-1, 1.53d-1, 1.72d-1,&!506.3 nm
1.90d-1,2.08d-1, 2.25d-1, 2.43d-1, 2.61d-1, &
2.78d-1, 2.96d-1,3.13d-1, 3.30d-1, 3.47d-1, 3.64d-1,&
1.21d-1, 1.45d-1, 1.69d-1, 1.92d-1, 2.14d-1,&!569.5 nm
2.36d-1,2.58d-1, 2.80d-1, 3.01d-1, 3.22d-1, &
3.43d-1, 3.64d-1,3.85d-1, 4.06d-1, 4.26d-1, 4.46d-1,&
1.56d-1, 1.88d-1, 2.18d-1, 2.47d-1, 2.76d-1,&!650.9 nm
3.03d-1,3.31d-1, 3.57d-1, 3.84d-1, 4.10d-1, &
4.36d-1, 4.62d-1,4.87d-1, 5.12d-1, 5.37d-1, 5.62d-1,&
2.10d-1, 2.53d-1, 2.93d-1, 3.32d-1, 3.69d-1,&!759.4 nm
4.06d-1, 4.41d-1, 4.75d-1, 5.09d-1, 5.43d-1, &
5.76d-1, 6.08d-1, 6.40d-1, 6.72d-1, 7.03d-1, 7.34d-1,&
2.98d-1, 3.59d-1, 4.16d-1, 4.70d-1, 5.22d-1,&!911.3 nm
5.73d-1,6.21d-1, 6.68d-1, 7.15d-1, 7.60d-1, 8.04d-1,&
8.47d-1,8.90d-1, 9.32d-1, 9.73d-1, 1.01d0,&
3.65d-1, 4.39d-1, 5.09d-1, 5.75d-1, 6.39d-1,&!1013 nm
7.00d-1, 7.58d-1, 8.15d-1, 8.71d-1, 9.25d-1, &
    9.77d-1, 1.03d0, 1.08d0, 1.13d0, 1.18d0, 1.23d0,&
    4.58d-1, 5.50d-1, 6.37d-1, 7.21d-1, 8.00d-1,&!1139 nm
    8.76d-1,9.49d-1, 1.02d0, 1.09d0, 1.15d0, 1.22d0,&
    1.28d0,1.34d0, 1.40d0, 1.46d0, 1.52d0,&
    5.92d-1, 7.11d-1, 8.24d-1, 9.31d-1, 1.03d0,&!1302 nm
    1.13d0, 1.23d0, 1.32d0, 1.40d0, 1.49d0,&
    1.57d0, 1.65d0, 1.73d0, 1.80d0, 1.88d0, 1.95d0,&
    7.98d-1, 9.58d-1, 1.11d0, 1.25d0, 1.39d0,&!1519 nm
    1.52d0, 1.65d0, 1.77d0, 1.89d0, 2.00d0, &
    2.11d0, 2.21d0, 2.32d0, 2.42d0, 2.51d0, 2.61d0,&
    1.14d0, 1.36d0, 1.58d0, 1.78d0, 1.98d0,&!1823 nm
    2.17d0, 2.34d0, 2.52d0, 2.68d0, 2.84d0,&
    3.00d0, 3.15d0, 3.29d0, 3.43d0, 3.57d0, 3.70d0,&
    1.77d0, 2.11d0, 2.44d0, 2.75d0, 3.05d0,&!2278 nm
    3.34d0, 3.62d0, 3.89d0, 4.14d0, 4.39d0,&
    4.63d0, 4.86d0, 5.08d0, 5.30d0, 5.51d0, &
    5.71d0,3.10d0, 3.71d0, 4.29d0, 4.84d0, 5.37d0,&!3038 nm
    5.87d0, 6.36d0, 6.83d0, 7.28d0, 7.72d0, &
    8.14d0, 8.55d0,8.95d0, 9.33d0, 9.71d0, 1.01d1,&
    6.92d0, 8.27d0, 9.56d0, 1.08d1, 1.19d1,&!4556 nm
    1.31d1,1.42d1, 1.52d1, 1.62d1, 1.72d1, 1.82d1,&
    1.91d1,2.00d1, 2.09d1, 2.17d1, 2.25d1,&
    2.75d1, 3.29d1, 3.80d1, 4.28d1, 4.75d1,&!9113 nm
    5.19d1,5.62d1, 6.04d1, 6.45d1, 6.84d1, 7.23d1,&
    7.60d1,7.97d1, 8.32d1, 8.67d1, 9.01d1 /

  chi = 0d0

  pe = atmos%ne(icell) * KBOLTZMANN * atmos%T(icell)
  ! 2-dimensionalize kappaFF_flat for interp2D
  kappaFF = RESHAPE(kappaFF_flat,(/NTHETA, NFF/))
  !!write(*,*) "Check reshape + interp2D"
  !do index=1,NFF
  ! do index2=1,NTHETA
  !  kappaFF(index,index2) = &
  !     kappaFF_flat(NTHETA*(index-1)+index2)
  ! end do
  !end do
  
  ! for long wavelengths
  ! Can do more efficiently !! that computing for all wavelength
  ! as it was all long wavelength and then computing
  ! only where it is below for the other wavelengths
  if ((MAXVAL(NLTEspec%lambda) >= lambdaFF(NFF)) .or. &
       (MINVAL(NLTEspec%lambda) >= lambdaFF(NFF))) then !if at least one
   CALL Hminus_ff_longwavelength(icell, chi)
   !!RETURN
  end if

   theta(1:1) = THETA0 /  atmos%T(icell)
   ! interpolate kappaFF at theta and lambda using
   ! 2x 1D cubic Bezier splines
   kappa = interp2Darr(NTHETA, thetaFF,NFF,lambdaFF,kappaFF,&
                  1,theta,NLTEspec%Nwaves, NLTEspec%lambda)
   where(NLTEspec%lambda < lambdaFF(NFF))
   chi = (Hydrogen%n(1,icell)*1d-29) * pe * kappa(1,:)
   end where


 RETURN
 END SUBROUTINE Hminus_ff

 ! Now Molecular Hydrogen continuum opacities

END MODULE hydrogen_opacities
