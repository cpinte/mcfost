! Hydrogen, bound-free, free-free and continuum opacities
! - H b-f and f-f --> H b-f is now computed in metal.f90 with the other atoms
! - H^- b-f and f-f
!  free-free emissivity is obtained by
! -> eta_ff = chi_ff * Bplanck
!
!
! Opacities in m^-1 (chi)
! Emissivities in W/m3/Hz/sr (eta)
module gas_contopac


   use atom_type, only : AtomicContinuum, find_continuum, AtomType, hydrogen, helium, &
      PassiveAtoms, Npassiveatoms
   use elements_type, only : elems
   use grid, only : T, ne, nHmin, nHtot
   use constantes
   use utils, only : locate, interp, bilinear, linear_1D_sorted, Bpnu
   use occupation_probability, only : D_i, wocc_n
   use parametres, only : ldissolve, n_cells

   implicit none
	
   real, parameter :: lambda_limit_HI_rayleigh = 91.1763062831680!102.6 !nm
   real, parameter :: lambda_limit_HeI_rayleigh = 140.0 !
   real(kind=dp), dimension(:), allocatable :: Hray_lambda, HeIray_lambda
   real(kind=dp), dimension(:), allocatable :: alpha_geltman, alpha_wishart
   integer, parameter :: N_geltman = 34, N_wishart = 63
   real(kind=dp), dimension(N_wishart) :: lambdai_wishart, alphai_wishart
   real(kind=dp), dimension(N_geltman) :: lambdai_geltman, alphai_geltman
   real(kind=dp), parameter :: lambda_base = 500.0_dp
   real(kind=dp), dimension(:), allocatable :: exphckT !exp(-hc_k/T)

   data lambdai_geltman / 0.0, 50.0, 100.0, 150.0, 200.0, 250.0,  &
   300.0, 350.0, 400.0, 450.0, 500.0, 550.0,&
   600.0, 650.0, 700.0, 750.0, 800.0, 850.0,&
   900.0, 950.0, 1000.0, 1050.0, 1100.0,    &
   1150.0, 1200.0, 1250.0, 1300.0, 1350.0,  &
   1400.0, 1450.0, 1500.0, 1550.0, 1600.0,  &
   1641.9 /

   !in 1e-21
   data alphai_geltman / 0.0,  0.15, 0.33, 0.57, 0.85, 1.17, 1.52,&
   1.89, 2.23, 2.55, 2.84, 3.11, 3.35, 3.56,&
   3.71, 3.83, 3.92, 3.95, 3.93, 3.85, 3.73,&
   3.58, 3.38, 3.14, 2.85, 2.54, 2.20, 1.83,&
   1.46, 1.06, 0.71, 0.40, 0.17, 0.0 /

    !36
    !    data lambdai_wishart / 0.00,1250.00 ,  1750.00  , 2250.70  , 2750.81  , 3250.94,&
    !    3751.07  , 4251.20  , 4751.33   ,5251.46 ,  5751.59 ,  6251.73, &
    !    6751.86  , 7252.00  , 7752.13  , 8252.27  , 8752.40   ,9252.54,&
    !    9752.67 ,  10252.81  ,10752.95 , 11253.08 , 11753.22 , 12253.35,&
    !   12753.49  ,13253.62 , 13753.76 , 14253.90 , 14754.03  ,15254.17,&
    !   15504.24 , 15754.30  ,16004.37 , 16104.40 , 16204.43 , 16304.45 /

    !1d-18 cm^2 * (1d-2)**2 for m^2 = 1d-22
    !   data alphai_wishart /0.0  ,     5.431 ,    7.918    , 11.08  ,   14.46   ,  17.92, &
    !      21.35 ,    24.65   ,  27.77  ,   30.62   ,  33.17   ,  35.37,&
    !      37.17  ,   38.54  ,   39.48  ,   39.95   ,  39.95 ,    39.48,&
    !      38.53  ,   37.13 ,   35.28  ,   33.01 ,    30.34     ,27.33,&
    !      24.02  ,   20.46   ,  16.74   ,  12.95  ,   9.211   ,  5.677,&
    !      4.052 ,    2.575 ,    1.302   ,  .8697    , .4974   ,  .1989 /

   data lambdai_wishart / 250, 1500, 1750, 2000, 2250, 2500, 2750, 3000, 3250, 3500, 3750,      &
   4000, 4250, 4500, 4750, 5000, 5250, 5500, 5750, 6000, 6250, 6500,     &
   6750, 7000, 7250, 7500, 7750, 8000, 8250, 8500, 8750, 9000, 9250,     &
   9500, 9750, 10000, 10250, 10500, 10750, 11000, 11250, 11500, 11750,   &
   12000, 12250, 12500, 12750, 13000, 13250, 13500, 13750, 14000,        &
   14250, 14500, 14750, 15000, 15250, 15500, 15750, 16000, 16100, 16200, &
   16300 /
   data alphai_wishart / 5.431, 6.512, 7.918, 9.453, 11.08,12.75,14.46,16.19,17.92,19.65,21.35,   &
   23.02,24.65,26.24, 27.77, 29.23,30.62,31.94,33.17,34.32, 35.37,36.32,    &
   37.17,37.91,38.54,39.07,39.48,39.77,39.95, 40.01, 39.95, 39.77, 39.48,   &
   39.06,38.53,37.89,37.13,36.25,35.28,34.19,33.01,31.72,30.34,28.87,       &
   27.33,25.71,24.02,22.26,20.46,18.62,16.74,14.85,12.95,11.07,9.211,7.407, &
   5.677,4.052,2.575,1.302,0.8697,0.4974, 0.1989 /

   contains

   subroutine background_continua_lambda(icell, Nx, x, chiout, etaout)
      integer, intent(in) :: icell, Nx
      integer :: la, nat
      real(kind=dp), dimension(Nx), intent(in) :: x
      real(kind=dp), dimension(Nx), intent(out) :: chiout, etaout
      real(kind=dp), dimension(Nx) :: chi, eta, Bp
  
      Bp(:) = Bpnu(x, T(icell))
  
      !HI Rayleigh must be initialised first
      chiout(:) = thomson(icell) + HI_rayleigh(icell)
      etaout(:) = 0.0_dp
      
      if (associated(helium)) then
         chiout(:) = chiout(:) + HeI_rayleigh(icell)
      endif

      call Hydrogen_ff(icell, Nx, x, chi)
      chiout(:) = chiout(:) + chi(:)
      etaout(:) = etaout(:) + chi(:) * Bp(:)
  
  !     call Hminus_bf_wishart(icell, N, lambda, chi, eta) !->better with turbospec
      call Hminus_bf_geltman(icell,Nx, x, chi, eta) !->better with rh
      chiout(:) = chiout(:) + chi(:)
      etaout(:) = etaout(:) + eta(:)
  
      call Hminus_ff_bell_berr(icell, Nx, x, chi)
      chiout(:) = chiout(:) + chi(:)
      etaout(:) = etaout(:) + chi(:) * Bp(:)

      !now atomic LTE bound-free
      !elsewhere even if LTE

      return
   end subroutine background_continua_lambda

   function Thomson(icell)
    ! ------------------------------------------------------------- !
    ! Thomson scattering cross-section in non relativistic limit.
    ! (i.e., wavelength independent) x ne density.
    ! unit m^-1
    ! ------------------------------------------------------------- !
      integer, intent(in) 			:: icell
      real(kind=dp)				   :: Thomson

      Thomson = ne(icell) * sigma_e

      return
   end function Thomson

   function HI_rayleigh(icell)
      integer, intent(in) :: icell 
      real(kind=dp) :: HI_rayleigh(size(Hray_lambda))

      HI_rayleigh = sigma_e * Hray_lambda(:) * Hydrogen%n(1,icell)

      return
   end function HI_rayleigh

   function HeI_rayleigh(icell)
      integer, intent(in) :: icell
      real(kind=dp) :: HeI_rayleigh(size(HeIray_lambda))

      !if helium has no neutral stage HeI_rayleigh is 0.

      HeI_rayleigh = sigma_e *  HeIray_lambda(:) * Helium%n(1,icell)

      return
   end function HeI_rayleigh

   subroutine alloc_gas_contopac(N, lambda)
      integer, intent(in) :: N
      real(kind=dp), intent(in) :: Lambda(N)

      allocate(Hray_lambda(N))

      where(lambda > lambda_limit_HI_rayleigh)
         Hray_lambda = (1d0 + (156.6d0/lambda)**2.d0 + &
          (148.d0/lambda)**4d0)*(96.6d0/lambda)**4d0
      elsewhere
         Hray_lambda = 0.0
      end where  

      if (associated(helium)) then
         allocate(HeIray_lambda(N))
         where(lambda > lambda_limit_HeI_rayleigh)
            HeIray_lambda = 4.0 * (1d0 + (66.9d0/lambda)**2.d0 + &
              (64.1d0/lambda)**4d0)*(37.9d0/lambda)**4d0
         elsewhere
            HeIray_lambda = 0.0
         end where 
      endif

      !cheap !
      !allocate cross-sections and compute then on lambda grid.
      allocate(alpha_geltman(N), alpha_wishart(N))
      alpha_geltman = 0.0
      alpha_wishart = 0.0
      alpha_wishart = 1d-22 * linear_1D_sorted(N_wishart, lambdai_wishart, alphai_wishart, N, 10*lambda)
      alpha_geltman = 1d-21 * linear_1D_sorted(N_geltman, lambdai_geltman, alphai_geltman, N, lambda) !1e-17 cm^2 to m^2


      !for obvious reason a lambda base is used to avoid exp(-hc_k/T) to goes with zero at low T.
      allocate(exphckT(n_cells))
      exphckT(:) = exp(-hc_k/T/lambda_base)!exp(-hnu/kT) = exphckT**(lambda_base/lambda(nm)
    
      return 
   end subroutine alloc_gas_contopac

   subroutine dealloc_gas_contopac

      deallocate(Hray_lambda)
      if (allocated(HeIray_lambda)) deallocate(HeIray_lambda)

      !deallocate or not because depends on cells so unchanged if we change lambda grid...
      deallocate(exphckT)
      !anyway should not cost anything!

      deallocate(alpha_wishart, alpha_geltman)

      return
   end subroutine dealloc_gas_contopac
  
   elemental function Gaunt_bf(u, n_eff)
      ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
      ! See also Menzel & Pekeris 1935
      real(kind=dp), intent(in) :: n_eff,  u ! = n_Eff**2 * eps = hnu/Z/Z/E_RYDBERG - 1
      real(kind=dp) :: Gaunt_bf
      Gaunt_bf = 1d0 + 0.1728 * (n_eff**(-2./3.)) * (u+1d0)**(-2./3.) * (u-1d0) &
         - 0.0496*(n_eff**(-4./3.)) * (u+1d0)**(-4./3.) * (u*u + 4./3. * u + 1d0)

    !   if (Gaunt_bf <= 0d0 .or. Gaunt_bf > 2d0) then
    !    Gaunt_bf = 1d0
    !   end if

      ! need a proper extrapolation for when it is negative
      if (Gaunt_bf < 0.0) Gaunt_bf = 0.0
      if (Gaunt_bf > 2.0) Gaunt_bf = 1.0

      return
   end function Gaunt_bf


   Elemental function Gaunt_ff(lambda, Z, T)
   ! M. J. Seaton (1960), Rep. Prog. Phys. 23, 313
   !
   ! Note: There is a problem with this expansion at higher temperatures
   ! (T > 3.0E4 and longer wavelengths (lambda > 2000 nm). Set to
   ! 1.0 when the value goes below 1.0
      real(kind=dp),intent(in) :: lambda, T
      real(kind=dp) :: x, x3, y, Gaunt_ff
      integer, intent(in) :: Z

      x = (HC/(lambda * NM_TO_M)) / (E_RYDBERG * (Z)**(2d0))
      x3 = (x**(3.3333333d-1))
      y  = (2.0 * lambda * NM_TO_M * KB*T) / HC

      gaunt_ff = 1.0 + 0.1728*x3 * (1.0 + y) - &
         0.0496*(x3*x3) * (1.0 + (1.0 + y)*0.33333333*y)

      if (gaunt_ff <= 0d0 .or. gaunt_ff > 2d0) gaunt_ff = 1d0

      return
   end function Gaunt_ff

   elemental function  H_bf_Xsection(cont, lambda) result(alpha)
      Type (AtomicContinuum), intent(in) :: cont
      real(kind=dp), intent(in) :: lambda
      real(kind=dp) :: n_eff, g_bf, u, Z, u0, g_bf0, alpha, u1


      Z = real(cont%atom%stage(cont%i) + 1,kind=dp)

      if (cont%atom%ID=='H') then
         n_eff = real(cont%i,kind=dp)
      else
         n_eff = Z*sqrt(cont%atom%Rydberg / (cont%atom%E(cont%j) - cont%atom%E(cont%i)))
      endif

      u = n_eff**2 * HC / (NM_TO_M * lambda) / Z*Z / E_RYDBERG - 1
      ! u0 = n_eff*n_eff * HC / (NM_TO_M * cont%lambda0) / Z / Z / E_RYDBERG - 1.

      g_bf = Gaunt_bf(u, n_eff)
   !  g_bf0 = Gaunt_bf(u0, n_eff)

    !     if (lambda > cont%lambda0) then !linear  extrapolation of g_bf
    !       u1 = n_eff**2 * HPLANCK*CLIGHT / (NM_TO_M * 0.8 * cont%lambda0 ) / Z*Z / E_RYDBERG - 1
    !       g_bf = g_bf0 + (u - u1) / (u0 - u1) * (g_bf0 -  Gaunt_bf(u1, n_eff))
    !     endif

    !There is a factor n_eff/Z**2 absorbed in alpha0 beware
    !alpha = cont%alpha0 * (lambda/cont%lambda0)**3  * g_bf / g_bf0
    !alpha = n_eff/Z**2 * cont%alpha0 * (lambda/cont%lambda0)**3  * g_bf / g_bf0
    !1d-4* 2.815d29
      alpha = 2.815d25 * (Z**4) * g_bf / n_eff**5 * (NM_TO_M*lambda/C_LIGHT)**3

      return
   end function H_bf_Xsection
  

   elemental function  H_ff_Xsection(Z, T, lambda) result(alpha)
      real(kind=dp), intent(in) :: lambda, T
      real(kind=dp) :: g_ff, nu3, alpha, K0
      integer, intent(in) :: Z

    !    K0 = (Q_ELECTRON**2)/(4.0*PI*EPSILON_0) / sqrt(M_ELECTRON)
    !    K0 = (K0**3) * 4./3. * sqrt(2*pi/3./KBOLTZMANN) / HPLANCK / CLIGHT
    !sigma0_H_ff = K0

      nu3 = (NM_TO_M * lambda / C_LIGHT)**3 ! = 1 / nu**3

      g_ff = Gaunt_ff(lambda, Z, T)

      alpha =  sigma0_H_ff * real(Z) * real(Z) * nu3 * g_ff / sqrt(T)

      return
   end function H_ff_Xsection



   !To do, add contributions from dissolve states to dissolve states
   subroutine Hydrogen_ff(icell, N, lambda, chi)
   ! Hubeny & Mihalas eq. 7.100 (from cgs to SI)
   ! takes place at LTE because it is collisional
      integer, intent(in) :: icell, N
      real(kind=dp), intent(in) :: lambda(N)
      real(kind=dp), intent(out) :: chi(N)
      integer :: la
      real(kind=dp) :: stim, np, arg_exp, exp_val, C0, alpha

      np = Hydrogen%n(Hydrogen%Nlevel,icell) !nH+=Nion H
      !should be nstar instead of %n ?

      if (ne(icell) == 0.0) then
         chi = 0.0
         return
      endif

      ! do la=1,N
      !    ! stim = 1. - exp(-hc_k/lambda(la)/T(icell))
      !    stim = 1.0 - exphckT(icell)**(lambda_base/lambda(la))

      !    ! = alpha0 /nu**3 / sqrt(T) = m^5
      !    !I now consider ne as part of the cross-sections to be in the same units as
      !    !the bound-free cross-sections
      !    alpha = H_ff_Xsection(1, T(icell), lambda(la)) * ne(icell)


      !    chi(la) =  alpha * np * stim
      ! enddo
      chi(:) = np * H_ff_Xsection(1, T(icell), lambda(:)) * ne(icell) * &
         (1.0 - exphckT(icell)**(lambda_base/lambda(:)) )
    return
   end subroutine Hydrogen_ff

  
   subroutine atom_ff

      return
   end subroutine atom_ff

   subroutine Hminus_bf_john(icell, N, lambda, chi, eta)
   !-----------------------------------------------------------------
   ! Calculates the negative hydrogen (H-) bound-free continuum
   ! absorption coefficient per
   ! hydrogen atom in the ground state from
   !  John 1988 A&A 193, 189
   ! Includes stimulated emission
   !-----------------------------------------------------------------
      integer, intent(in) :: icell, N
      real(kind=dp), intent(in), dimension(N) :: lambda
      real(kind=dp) :: lam, lambda0, alpha, sigma, flambda, Cn(6)
      real(kind=dp) :: diff, stm, funit, cte, pe
      integer :: la
      real(kind=dp) :: arg_exp, nH
      real(kind=dp), dimension(N), intent(out) :: chi, eta

      chi(:) = 0.0_dp
      eta(:) = 0.0_dp

      !1dyne/cm2 = 1e-1 Pa
      !1dyne = 1e-1 Pa cm2
      !cm4/dyne = cm4 / (1e-1 Pa cm2)
      !cm4/dyne = 10 * cm2 / Pa
      !cm4/dyne = 10 * (1e-2)**2 m2/Pa = 1e-3 * m2/Pa
      !m2/Pa = cm4/dyne  * 1e3
      funit = 1d-3 !m2/Pa -> cm2/dyn

      pe = ne(icell) * KB * T(icell) !nkT in Pa

      Cn(:) = (/152.519d0,49.534d0,-118.858d0,92.536d0,-34.194d0,4.982d0/)

      !check here km_to_m = micron_to_nm
      alpha = hc_k / km_to_m ! hc_k = hc/k/nm_to_m
      lambda0 = 1.6419 !micron, photo detachement threshold

      nH = hydrogen%n(1,icell)

    !alpha = hc_k / micron_to_nm * nm_to_m
      cte = 0.75d-18 * T(icell)**(-2.5) * exp(alpha/lambda0/T(icell)) * pe * funit * nH


      do la=1, N

         lam = lambda(la) / km_to_m !nm->micron
         !if (lambda > 0.125 .and. lambda < lambda0) then
         if (lam <= lambda0) then

            diff = (1d0/lam - 1d0/lambda0)

            stm = 1. - exp(-alpha/lam/T(icell))

            flambda = Cn(1) + Cn(2) * diff**(0.5) + Cn(3) * diff + Cn(4) * diff**(1.5) + &
               Cn(5)*diff**(2.) + Cn(6) * diff**(2.5)

            sigma = lam**3d0 * diff**(1.5) * flambda !cm2
            chi(la) = cte * stm * sigma! m^-1
            !exp(-hnu/kt) * 2hnu3/c2 = Bp * (1.-exp(-hnu/kt))
            eta(la) = cte * (1.-stm) * sigma * twohc/lambda(la)**3

         endif

      enddo



      return
   end subroutine Hminus_bf_john

   subroutine Hminus_bf_geltman(icell, N, lambda, chi, eta)
    !-----------------------------------------------------------------
    ! Calculates the negative hydrogen (H-) bound-free continuum
    ! absorption coefficient per
    ! H minus atom from Geltman 1962, ApJ 136, 935-945
    ! Stimulated emission included
    !-----------------------------------------------------------------
    integer, intent(in) :: icell, N
    real(kind=dp), dimension(N), intent(in) :: lambda
    integer :: la
    real(kind=dp), intent(out), dimension(N) :: chi, eta
    real(kind=dp) :: lam, stm, twohnu3_c2

    

    do la=1, N
       lam = lambda(la)
       !do not test negativity of lambda
       if (lam >= lambdai_geltman(N_geltman)) then
         chi(la:N) = 0.0
         eta(la:N) = 0.0
         exit
       endif

       stm = exphckT(icell)**(lambda_base/lam) !exp(-hc_k/T(icell)/lam)
       twohnu3_c2 = twohc / lam**3.

       chi(la) = nHmin(icell) * (1.-stm) * alpha_geltman(la)
       eta(la) = nHmin(icell) * twohnu3_c2 * stm * alpha_geltman(la)

    enddo


      return
   end subroutine Hminus_bf_geltman


   subroutine Hminus_bf_Wishart(icell, N, lambda, chi, eta)
    !The one use in Turbospectrum, number 1 in opacity
    !-----------------------------------------------------------------
    ! Calculates the negative hydrogen (H-) bound-free continuum
    ! absorption coefficient per
    ! H minus atom from Wishart A.W., 1979, MNRAS 187, 59P
    !-----------------------------------------------------------------
    integer, intent(in) :: icell, N
    integer :: la
    real(kind=dp), dimension(N), intent(in)	:: lambda
    real(kind=dp), dimension(N), intent(out) :: chi, eta
    real(kind=dp) :: lam, stm, chi_extr(1), eta_extr(1), lambda_extr(1)

    chi(:) = 0.0_dp
    eta(:) = 0.0_dp

    freq_loop : do la=1, N
       lam = lambda(la) * 10. !AA
       !stm = 0.0
       stm = exphckT(icell)**(lambda_base/lambda(la)) !exp(-hc_k / T(icell) / lambda(la))

       !if (lam < minval(lambdai_wishart)) then
       if (lam < lambdai_wishart(1)) then
          !cyle
          !other formula for very low frequencies
          lambda_extr(1) = lambda(la)
          call Hminus_bf_geltman(icell, 1, lambda_extr, chi_extr, eta_extr)
          chi(la) = chi_extr(1)
          eta(la) = eta_extr(1)
          !cycle
          !else if (lam > maxval(lambdai)) then
       elseif (lam > lambdai_wishart(N_wishart)) then
          exit freq_loop
       else

         chi(la) = alpha_wishart(la) * (1.0 - stm) * nHmin(icell)
         eta(la) = alpha_wishart(la) * twohc/lambda(la)**3  * stm * nHmin(icell)

       endif
    enddo freq_loop


      return
   end subroutine Hminus_bf_Wishart

   elemental function Hminus_ff_john_lam(icell, lambda) result(chi)
   !-----------------------------------------------------------------
   ! Calculates the negative hydrogen (H-) free-free continuum
   ! absorption coefficient per
   ! hydrogen atom in the ground state from John 1988
   !-----------------------------------------------------------------
      integer, intent(in) :: icell
      real(kind=dp), intent(in) :: lambda
      !tab. 3a, lambda > 0.3645micron
      real(kind=dp), dimension(6) :: An, Bn, Cn, Dn, En, Fn
      !tab 3b, lambda > 0.1823 micron and < 0.3645 micron, size of 4 isntead 5, because the last two rows are 0
      real(kind=dp), dimension(4) :: An2, Bn2, Cn2, Dn2, En2, Fn2
      real(kind=dp) :: funit, K, kappa
      integer :: la, n
      real(kind=dp) :: chi
      real(kind=dp) :: lam, theta, nH

      chi = 0.0
      theta = 5040d0 / T(icell)
      lam = lambda / km_to_m

      if (theta < 0.5 .or. theta > 3.6) return
      if (lam <= 0.1823) return

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
      nH = hydrogen%n(1,icell)

      K = 1d-29 * funit * ne(icell) * KB * T(icell) * nH


      kappa = 0.0_dp
      if (lam < 0.3645) then
         do n=1,4
            kappa = kappa + theta**((n+1)/2d0) * (lam**2 * An2(n) + Bn2(n) + Cn2(n)/lam + &
               Dn2(n)/lam**2 + En2(n)/lam**3 + Fn2(n)/lam**4)
         enddo
      else! if (lambda > 0.3645) then
         do n=1,6
            kappa = kappa + theta**((n+1)/2d0) * (lam**2 * An(n) + Bn(n) + Cn(n)/lam + &
               Dn(n)/lam**2 + En(n)/lam**3 + Fn(n)/lam**4)
         enddo
      endif
      chi = kappa * K


      return
   end function Hminus_ff_John_lam

   subroutine Hminus_ff_john(icell, N, lambda, chi)
   !-----------------------------------------------------------------
   ! Wrapper around Hminus_ff_john
   !-----------------------------------------------------------------
      integer, intent(in) :: icell, N
      real(kind=dp), dimension(N), intent(in) :: lambda
      integer :: la
      real(kind=dp), dimension(N), intent(out) :: chi
      real(kind=dp) :: lam, theta!, nH

      chi(:) = 0.0_dp

      theta = 5040d0 / T(icell)

      if (theta < 0.5 .or. theta > 3.6) return

      chi(:) = Hminus_ff_john_lam(icell, lambda)

      return
   end subroutine Hminus_ff_john

   subroutine Hminus_ff_bell_berr(icell, N, lambda, chi)
   !-----------------------------------------------------------------
   ! Calculates the negative hydrogen (H-) free-free continuum
   ! absorption coefficient per
   ! hydrogen atom in the ground state from
   ! Bell & Berrington 1986, J. Phys. B 20, 801
   !
   ! at long wavelength Hminus_ff_john is used.
   !
   !-----------------------------------------------------------------
    !stm is already included

      integer, intent(in) :: icell, N
      real(kind=dp), dimension(N), intent(in) :: lambda
      integer :: la
      real(kind=dp), dimension(N), intent(out) :: chi
      real(kind=dp), dimension(:) :: lambdai(23), thetai(11)
      real(kind=dp), dimension(23,11) :: alphai
      real :: inter
      real(kind=dp) :: lam, stm, sigma, theta, pe, nH

      chi(:) = 0.0_dp
      theta = 5040. / T(icell)

      if (theta < 0.5 .or. theta > 3.6) then
         !return
         if (theta < 0.5) theta = 0.5_dp
         if (theta > 3.6) theta = 3.6_dp
      endif

      !AA, index 12 (11 in C) is 9113.0 AA
      data lambdai  /0.00, 1823.00, 2278.0, 2604.0, 3038.0, 3645.0,			&
         4557.0, 5063.0, 5696.0, 6510.0, 7595.0, 9113.0, 					&
         10126.0, 11392.0, 13019.0, 15189.0, 18227.0, 22784.0,				&
         30378.0, 45567.0, 91134.0, 113918.0, 151890.0						/

      data thetai /0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.8, 3.6		/

      !1d-26 cm^4/dyn = 1d-26 * 1d-3 m^2 / Pa = 1d-29 m^2/Pa
      !t = 0.5
      data alphai(:,1) /0.0, 1.78e-2, 2.28e-2, 2.77e-2, 3.64e-2, 5.20e-2,		&
         7.91e-2, 9.65e-2, 1.21e-1, 1.54e-1, 2.08e-1, 2.93e-1,			&
         3.58e-1, 4.48e-1, 5.79e-1, 7.81e-1, 1.11, 1.73,					&
         3.04, 6.79, 2.70e1, 4.23e1, 7.51e1								/
      !t = 0.6
      data alphai(:,2) /0.0, 2.22e-2, 2.80e-2, 3.42e-2, 4.47e-2, 6.33e-2,		&
         9.59e-2, 1.17e-1, 1.46e-1, 1.88e-1, 2.50e-1, 3.54e-1, 4.32e-1,	&
         5.39e-1, 6.99e-1, 9.40e-1, 1.34, 2.08, 3.65, 8.16, 3.24e1,		&
         5.06e1, 9.0e1													/
      !t = 0.8
      data alphai(:,3) /0.0, 3.08e-2, 3.88e-2, 4.76e-2, 6.16e-2, 8.59e-2,		&
         1.29e-1, 1.57e-1, 1.95e-1, 2.49e-1, 3.32e-1, 4.68e-1, 5.72e-1,	&
         7.11e-1, 9.24e-1, 1.24, 1.77, 2.74, 4.80, 1.07e1, 4.26e1, 		&
         6.64e1, 1.18e2													/
      !t = 1.0
      data alphai(:,4) /0.0, 4.02e-2, 4.99e-2, 6.15e-2, 7.89e-2, 1.08e-1,		&
         1.61e-1, 1.95e-1, 2.41e-1, 3.09e-1, 4.09e-1, 5.76e-1, 7.02e-1,		&
         8.71e-1, 1.13, 1.52, 2.17, 3.37, 5.86, 1.31e1, 5.19e1, 8.08e1,		&
         1.44e2																/
      !t = 1.2
      data alphai(:,5) /0.0, 4.98e-2, 6.14e-2, 7.60e-2, 9.66e-2, 1.31e-1,		&
         1.94e-1, 2.34e-1, 2.88e-1, 3.67e-1, 4.84e-1, 6.77e-1, 8.25e-1, 	&
         1.02, 1.33, 1.78, 2.53, 3.90, 6.86, 1.53e1, 6.07e1, 9.45e1,		&
         1.68e2															/
      !t = 1.4
      data alphai(:,6) /0.0, 5.96e-2, 7.32e-2, 9.08e-2, 1.14e-1, 1.54e-1, 	&
         2.27e-1, 2.72e-1, 3.34e-1, 4.24e-1, 5.57e-1, 7.77e-1, 9.43e-1, 	&
         1.16, 1.51, 2.02, 2.87, 4.50, 7.79, 1.74e1, 6.89e1, 1.07e2, 	&
         1.91e2															/
      !t = 1.6
      data alphai(:,7) /0.0, 6.95e-2, 8.51e-2, 1.05e-1, 1.32e-1, 1.78e-1, 	&
         2.60e-1, 3.11e-1, 3.81e-1, 4.82e-1, 6.30e-1, 8.74e-1, 1.06,		&
         1.29, 1.69, 2.26, 3.20, 5.01, 8.67, 1.94e1, 7.68e1, 1.20e2, 	&
         2.12e2															/

      !t = 1.8
      data alphai(:,8) /0.0, 7.95e-2, 9.72e-2, 1.21e-1, 1.50e-1, 2.01e-1, 	&
         2.93e-1,3.51e-1, 4.28e-1, 5.39e-1, 7.02e-1,9.69e-1,1.17, 1.43, 	&
         1.86,2.48, 3.51,5.50, 9.50, 2.12e1, 8.42e1, 1.31e2, 2.34e2		/

      !t = 2.0
      data alphai(:,9) /0.0, 8.96e-2, 1.1e-1, 1.36e-1, 1.69e-1, 2.25e-1, 		&
         3.27e-1, 3.9e-1, 4.75e-1, 5.97e-1, 7.74e-1, 1.06, 1.28, 1.57, 	&
         2.02, 2.69, 3.8, 5.95, 1.03e1, 2.30e1, 9.14e1, 1.42e2, 2.53e2	/

      !t = 2.8
      data alphai(:,10) /0.0, 1.31e-1, 1.60e-1, 1.99e-1, 2.43e-1, 3.21e-1, 	&
         4.63e-1, 5.49e-1, 6.67e-1, 8.30e-1, 1.06, 1.45, 1.73, 2.09, 	&
         2.67, 3.52, 4.92, 7.59, 1.32e1, 2.95e1, 1.17e2, 1.83e2, 		&
         3.25e2															/

      !t = 3.6
      data alphai(:,11) /0.0, 1.72e-1, 2.11e-1, 2.62e-1, 3.18e-1, 4.18e-1,	&
         6.02e-1, 7.11e-1, 8.61e-1, 1.07, 1.36, 1.83, 2.17, 2.60, 3.31, 	&
         4.31, 5.97, 9.06, 1.56e1, 3.50e1, 1.40e2, 2.19e2, 3.88e2		/


      pe = ne(icell) * KB * T(icell)
      nH = hydrogen%n(1,icell)

      do la=1, N
         lam = lambda(la) * 10.
         if (lam > lambdai(23)) then
            chi(la) = Hminus_ff_john_lam(icell,lambda(la))
         else
            stm = exphckT(icell)**(lambda_base/lambda(la)) !exp(-hc_k/T(icell)/lambda(la))
            sigma = 1d-29 * bilinear(23,lambdai,11,thetai,alphai,lam,theta) !m^2/Pa
            chi(la) = sigma * pe * nH!m^-1
         endif
      enddo

      return
   end subroutine Hminus_ff_bell_berr


end module gas_contopac