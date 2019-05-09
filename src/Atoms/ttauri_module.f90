MODULE TTauri_module

  use parametres
  use atmos_type, only : atmos
  use Voronoi_grid, only : Voronoi
  use utils, only : interp_dp
  use math, only : locate
  use constantes
  use stars
  use grid, only : r_grid, z_grid

  IMPLICIT NONE

  CONTAINS  
  
  SUBROUTINE TTauri_Temperature(rmi, rmo, Macc)
  ! ---------------------------------------------------------------- !
   ! Computes temperature of accretion columns of T Tauri stars
   ! given the density, using the prescription of Hartmann 94
   ! and the cooling rates for Hartmann 82
   !
   ! TO DO: Numerically solve for the temperature, iterating
   ! between NLTE solution (Cooling rates) and Temperature->ne->npops
  ! ------------------------------------------------------------------ !
   integer 	:: icell, icell0
   double precision, intent(in) :: rmi, rmo, Macc
   									!with the one in magneto_accretion_model() if used
   double precision :: Q0, nH0, L, TL(8), Lambda(8) !K and erg/cm3/s
   double precision :: Rstar, Mstar, Mdot, Lr, thetai, thetao, Tring, rho_to_nH, Sr
   double precision, parameter ::  year_to_sec = 3.154d7, r0=1d0
   double precision, allocatable, dimension(:) :: r
   
   !T = K
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 /
   !Lambda = Q/nH**2, Q in J/m9/s
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 /

   allocate(r(atmos%Nspace))
   if (lVoronoi) then
    r(:) = dsqrt(Voronoi(:)%xyz(1)**2 + &
    	Voronoi(:)%xyz(2)**2+Voronoi(:)%xyz(3)**2 )
   else
    r = dsqrt(r_grid**2 + z_grid**2)
   end if

   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot = Macc * Msun_to_kg / year_to_sec !kg/s
   Lr = Ggrav * Mstar * Mdot / Rstar * (1. - 2./(rmi + rmo))

   thetai = asin(dsqrt(1d0/rmi)) !rmi and rmo in Rstar
   thetao = asin(dsqrt(1d0/rmo))

   Tring = Lr / (4d0 * PI * Rstar**2 * sigma * abs(cos(thetai)-cos(thetao)))
   Tring = Tring**0.25

   Sr = abs(cos(thetai)-cos(thetao)) !4pi*Rstar**2 * abs(c1-c2) / 4pi Rstar**2
   rho_to_nH = 1d3/masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3
   
   !Arbitray definition, could be varying also
   nH0 = rho_to_nH * (Mdot * Rstar) /  (4d0*PI*(r0/rmi - r0/rmo)) * &
                        (Rstar * r0)**( real(-5./2.) ) / dsqrt(2d0 * Ggrav * Mstar) * &
                        dsqrt(4d0-3d0*(r0/rmi)) / dsqrt(1d0-(r0/rmi))
   !!Q0 in J/m9/s
   !!Q0 = nH0**2 * 10**(interp1D(TL, Lambda, log10(Tring)))*0.1 !0.1 = erg/cm3 to J/m3
   Q0 = nH0**2 * 10**(interp_dp(Lambda, TL, log10(Tring)))*0.1

   
   write(*,*) "Ring T", Tring, "K"
   write(*,*) "Surface", 100*Sr, "%"
   write(*,*) "Luminosity", Lr/Lsun, "Lsun"   
   write(*,*) "nH0", nH0,"m^-3"
   
   atmos%T = 0d0

   do icell=1, atmos%Nspace
    if ((atmos%nHtot(icell)>0).and.(r(icell)/etoile(1)%r>=r0)) then !Iam assuming that it is .true. only
    									 !along field lines
     L = 10 * Q0*(r0*etoile(1)%r/r(icell))**3 / atmos%nHtot(icell)**2!erg/cm3/s
     atmos%T(icell) = 10**(interp_dp(TL, Lambda, log10(L)))
     !write(*,*) icell, L, atmos%T(icell), atmos%nHtot(icell), nH0, r
    end if
   end do
   
   deallocate(r)
!stop
  RETURN
  END SUBROUTINE TTauri_Temperature
  
  SUBROUTINE intersect_spots(i_star,u,v,w,x,y,z,ispot,lintersect)
  ! ------------------------------------------------------------ !
   ! Will a ray/packet hit a spot ?
   ! suppose that no spot overlap.
   
  ! ------------------------------------------------------------ !
   real(kind=dp), intent(in) :: x,y,z, u,v,w
   logical, intent(out) :: lintersect
   integer, intent(in) :: i_star
   integer, intent(out) :: ispot

   real(kind=dp), dimension(3) :: r, k, delta_r
   real(kind=dp) :: b,c, delta, rac, s1, s2
   integer :: i
   
   r(1) = x ; r(2) = y ; r(3) = z
   k(1) = u ; k(2) = v ; k(3) = w

   ispot = 0
   lintersect = .false.
   spot_loop : do i = 1, etoile(i_star)%Nspot
!      delta_r(:)  = r(:) - (/etoile(i_star)%Starspots(i)%xs, \
!                         etoile(i_star)%Starspots(i)%ys, etoile(i_star)%Starspots(i)%zs/)
!      b = dot_product(delta_r,k)
!      c = dot_product(delta_r,delta_r) - (etoile(i_star)%Starspots(i)%rs)**2
!      delta = b*b - c
     if (delta>=0) then
      ispot = i
      lintersect = .true.
     end if
   end do spot_loop
     
  RETURN
  END SUBROUTINE intersect_spots
  
  
END MODULE TTauri_module