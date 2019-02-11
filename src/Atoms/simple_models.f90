MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos, define_atomRT_domain
 
 ! MCFOST's modules
 use input
 use parametres
 use grid
 use density
 use constantes
 use fits_utils, only : print_error
 
 IMPLICIT NONE

 CONTAINS
 
  SUBROUTINE rotateZ(Vx, Vy, Vz, angle)
   !rotation around z-axis = rotation_3d((/0,0,1/), -angle, (/Vx, Vy, Vz/))
   double precision, intent(inout) :: Vx, Vy, Vz
   double precision, intent(in)    :: angle
   double precision				   :: x, y, z
   
   x = Vx
   y = Vy
   z = Vz
   Vz = z
   Vx = cos(angle) * x + sin(angle) * y
   Vy = -sin(angle) * x + cos(angle) * y
  
  RETURN
  END SUBROUTINE rotateZ
 
  SUBROUTINE magneto_accretion_model()
  ! ----------------------------------------------------------- !
   ! Magnetospheric accretion model
   ! Field line equation:
   !	r(theta, Rm) = Rm * sin2(theta), with
   ! sintheta = R/sqrt(z**2+R**2) and r/sin2(theta) = cte along
   ! a line.
   !
   ! No rotation of the star yet
  ! ----------------------------------------------------------- !
   use math, only : interp1D
   use utils, only : interp_dp
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Tmax = 8d3, days_to_sec = 86400d0, Prot = 8. !days
   double precision, parameter :: rmi=2.2d0, rmo=3.0d0, Tshk=0d0, Macc = 1d-7
   double precision, parameter :: year_to_sec = 3.154d7, r0 = 1d0, deltaz = 0.2!Rstar
   double precision ::  OmegasK, Rstar, Mstar, thetao, thetai, Lr, Tring, Sr, Q0, nH0
   double precision :: vp, y, rcyl, z, r, phi, Mdot, sinTheta, Rm, L, proj_phix, proj_phiy
   double precision :: Omega, Vphi, NormV(n_cells), TL(8), Lambda(8) !K and erg/cm3/s
   
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 / !log10 ?
   !Lambda = Q/nH**2
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 / !log10? 
   
   !if (.not.allocated(Vfiedl)) allocate(Vfield(n_cells))
   Omega = 2.*PI / (Prot * days_to_sec) ! Angular velocity in rad/s

   CALL init_atomic_atmos(n_cells)
   atmos%vturb = 5d3 !m/s
   !vx = 0d0; vy=0d0; vz=0d0

   !Define some useful quantities
   !write(*,*) "Mstar ", etoile(1)%M, " Rstar ", etoile(1)%r*AU_to_Rsun
   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot = Macc * Msun_to_kg / year_to_sec !kg/s
   Lr = Ggrav * Mstar * Mdot / Rstar * (1. - 2./(rmi + rmo))
   !sqrt(2*GMstar/rstar)
   OmegasK = dsqrt(Mstar * Ggrav * 2d0 / Rstar)
   !maximum and minimum angle on the stellar surface
   !all field lines are between theta0 and thetai at the stellar surface (r=Rstar)
   thetai = asin(dsqrt(1d0/rmi)) !rmi and rmo in Rstar
   thetao = asin(dsqrt(1d0/rmo))
   !write(*,*) "Thetamax/Thetamin", thetai*rad_to_deg, thetao*rad_to_deg
   if (Tshk <= 0d0) then
    Tring = Lr / (4d0 * PI * Rstar**2 * sigma * abs(cos(thetai)-cos(thetao)))
    Tring = Tring**0.25
   else
    Tring = Tshk
   end if
   Sr = abs(cos(thetai)-cos(thetao)) !4pi*Rstar**2 * abs(c1-c2) / 4pi Rstar**2
   write(*,*) "Ring T ", Tring, "K"
   write(*,*) "Surface ", 100*Sr, "%"
   write(*,*) "Luminosity", Lr/Lsun, "Lsun"
   write(*,*) "Mean molecular weight", atmos%avgWeight
   !now nH0 and Q0 are computed for each field lines assuming that Tring is the same
   !for all.
!    nH0 = 1d3/masseH/atmos%avgWeight * (Mdot * Rstar) /  (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
!                        (Rstar * r0)**( real(-5./2.) ) / dsqrt(2d0 * Ggrav * Mstar) * &
!                        dsqrt(4d0-3d0*(r0/rmi)) / dsqrt(1d0-(r0/rmi))
!    if (dabs(1d0 - (r0/rmi)) <= 1e-15) write(*,*) "Error nH0",nH0
!    !!Q0 in J/m9/s
!    !!Q0 = nH0**2 * 10**(interp1D(TL, Lambda, log10(Tring)))*0.1 !0.1 = erg/cm3 to J/m3
!    Q0 = nH0**2 * 10**(interp_dp(Lambda, TL, log10(Tring)))*0.1

    do i=1, n_rad
     do j=j_start,nz !j_start = -nz in 3D
      do k=1, n_az
       if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
       else
        icell = cell_map(i,j,k)
        rcyl = r_grid(icell)
        z = z_grid(icell)/z_scaling_env
       end if
       phi = phi_grid(icell)
       r = dsqrt(z**2 + rcyl**2)
       Vphi = 235d0! Omega * (r*AU_to_m) !m/s
       !from trigonometric relations we get y using rcyln 
       !and the radius at r(rcyln, zn)
       sinTheta = rcyl/r
       ! and cosTheta = z/r = sqrt(1.-y) = sqrt(1.-sinTheta**2)
       proj_phix = -sin(phi) !cos(phi)*sinTheta
       proj_phiy =  cos(phi) !sin(phi)*sinTheta
       !from definition of y
       y = sinTheta**2!min(real(sinTheta**2), 0.9999)
       !here the trick is to handle y = 1 at z=0 because rcyl=r and 1/(1-y) = inf.
       !r is a field line radial distance (from the centre of the model to the centre
       !of the cell.
       !if r**3/rcyl**2 or Rm=r/sint**2 is between rmi and rmo it means that the
       !cell icell (r, theta) belongs to a field line.
       Rm = r**3 / rcyl**2 / etoile(1)%r !in Rstar, same as rmi,o
       if (Rm>=rmi .and. Rm<=rmo) then !r>etoile(1)%r
          nH0 = 1d3/masseH/atmos%avgWeight * (Mdot * Rstar) /  (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (Rstar * r0)**( real(-5./2.) ) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*(r0/Rm)) / dsqrt(1d0-(r0/Rm))
          Q0 = nH0**2 * 10**(interp_dp(Lambda, TL, log10(Tring)))*0.1
          !write(*,*) icell, Rm, nH0, Q0, log10(Q0/nH0**2)
          !Interface funnels/disk
          if (dabs(z*AU_to_m/Rstar)<=deltaZ) then!(z==0d0) then
           !at z=0, r = rcyl = Rm here so 1/r - 1/Rm = 0d0 and y=1: midplane
           vp = 0d0
           !atmos%Vxyz is initialized to 0 everywhere
           if (.not.lstatic) then
            !V . xhat = (0rhat,0thetahat,Omega*r phihat) dot (..rhat,..thetahat,-sin(phi)phihat)
            !V . yhat = (0rhat,0thetahat,Omega*r phihat) dot (..rhat,..thetahat,cos(phi)phihat)
            !V .zhat = (0rhat,0thetahat,Omega*r phihat) dot (..rhat,..thetahat,0phihat)
            atmos%Vxyz(icell,1) = Vphi * proj_phix
            atmos%Vxyz(icell,2) = Vphi * proj_phiy
            !here it is better to have keplerian rotation than stellar rotation?
           end if
           !Density midplane of the disk
           !using power law of the dust-gas disk (defined in another zone?)
           atmos%nHtot(icell) = 5e18!
           !Temperature still given by the same law as in accretion funnels yet.
           !But probably different at this interface funnels/disk
          else !accretion funnels
          !Density
           atmos%nHtot(icell) = ( Mdot * Rstar )/ (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (AU_to_m * r)**( real(-2.5) ) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*y) / dsqrt(1d0-y) * &
                       1d3/masseH /atmos%avgWeight  !kg/m3 ->/m3
           vp = OmegasK * dsqrt(etoile(1)%r/r - 1./Rm)
           vp = -vp / dsqrt(4d0-3d0*y)
           if (.not.lstatic) then
             atmos%Vxyz(icell,1) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) * cos(phi)
             atmos%Vxyz(icell,2) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) * sin(phi)
             atmos%Vxyz(icell,3) = vp * (2d0 - 3d0 * y) !z
             ! or theta > PI/2d0
             if (z < 0) atmos%Vxyz(icell,3) = -atmos%Vxyz(icell,3)
             atmos%Vxyz(icell,1) = atmos%Vxyz(icell,1) + Vphi * proj_phix
             atmos%Vxyz(icell,2) = atmos%Vxyz(icell,2) + Vphi * proj_phiy
           end if
          end if
          L = 10 * Q0*(r0*etoile(1)%r/r)**3 / atmos%nHtot(icell)**2!erg/cm3/s
          !atmos%T(icell) = 10**(interp1D(Lambda, TL, log10(L)))
          atmos%T(icell) = 10**(interp_dp(TL, Lambda, log10(L)))
          !write(*,*) log10(L), Rm, rcyl/etoile(1)%r, atmos%T(icell), nH0, Q0
       end if
      end do
     end do
    end do
    !The disk could be included as new zone or by extending the density to r>rmo
    !but it has keplerian rotation and not infall but mcfost cannot treat two
    !kind of velocities right?
   
   !cannot use both simultaneously
   lkeplerian = .false.
   linfall = .false.
   if (maxval(atmos%vturb) > 0) atmos%v_char = minval(atmos%vturb,mask=atmos%vturb>0)
   if (.not.lstatic) then
    NormV = dsqrt(atmos%Vxyz(:,1)**2+atmos%Vxyz(:,2)**2+atmos%Vxyz(:,3)**2)
    atmos%v_char = atmos%v_char + minval(NormV,mask=NormV>0)
   end if
   write(*,*) "maxVfield (km/s)", maxval(NormV)/1d3, &
              " minVfield", minval(NormV,mask=normV>0)/1d3
   write(*,*) "atmos%v_char (km/s) =", atmos%v_char/1d3
   CALL define_atomRT_domain()

  RETURN
  END SUBROUTINE magneto_accretion_model
  
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid
   ! Values width idk are from a Solar model atmosphere by 
   ! Fontenla et al. namely the FAL C model.
   ! idk = 0, top of the atmosphere, idk = 81 (max) bottom.
  ! ----------------------------------------------------------- !
   CALL init_atomic_atmos(n_cells)

   !idk = 10
!    atmos%nHtot =  2.27414200581936d16
!    atmos%T = 45420d0
!    atmos%ne = 2.523785d16

   !idk = 75
!   atmos%T=7590d0
!   atmos%ne = 4.446967d20
!   atmos%nHtot = 1.259814d23

    !idk = 81   
    atmos%T=9400d0
    atmos%ne = 3.831726d21
    atmos%nHtot = 1.326625d23

    !idk = 0 
!     atmos%T = 100000d0
!     atmos%ne = 1.251891d16
!     atmos%nHtot = 1.045714d16

   CALL define_atomRT_domain()

   lstatic = .true. !force to be static for this case
   !tmos%vturb = 9.506225d3 !m/s !idk=10
   !atmos%vturb = 1.696164d3 !idk = 75
   atmos%vturb = 1.806787d3 !idk=81
   !atmos%vturb = 10.680960d3 !idk=0
   atmos%v_char = minval(atmos%vturb,mask=atmos%vturb>0)

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE prop_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to molecular/dust density and temperature.
  ! ----------------------------------------------------------- !
   double precision :: Vkepmax, Vkepmin
   integer :: icell
   
   CALL init_atomic_atmos(n_cells)

   atmos%nHtot = 1d14 * densite_gaz/MAXVAL(densite_gaz) + 1d12
   atmos%T = 3000d0+Tdust
   atmos%ne = 0d0

   lkeplerian = .true.
   atmos%vturb = 1d-5
   atmos%v_char = minval(atmos%vturb,mask=atmos%vturb>0)
   if (.not.lstatic) then
!     Vkepmax = dsqrt(Ggrav * sum(etoile%M) * Msun_to_kg * Rmin**2 / &
!                 ((Rmin**2 + minval(z_grid)**2)**1.5 * AU_to_m) )
    write(*,*) "Implement Vfield here"
!     if (lcylindrical_rotation) then ! Midplane Keplerian velocity
!       do icell=1, n_cells
!         Vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
!       end do
!     else ! dependance en z
!        do icell=1, n_cells
!          Vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg * r_grid(icell)**2 / &
!                 ((r_grid(icell)**2 + z_grid(icell)**2)**1.5 * AU_to_m) )
!        end do
!     end if
    !Vkepmin = minval(Vfield, mask=atmos%V>0)
    atmos%v_char = atmos%v_char! + Vkepmin
   end if
   CALL define_atomRT_domain() !deallocate atmos%V if static

  RETURN
  END SUBROUTINE prop_law_model

 END MODULE simple_models