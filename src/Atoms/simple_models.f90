MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos
 
 ! MCFOST's modules
 use input
 use parametres
 use grid
 use density
 use constantes
 use fits_utils, only : print_error
 
 IMPLICIT NONE

 CONTAINS
 
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
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Tmax = 8d3, Omegas = 0d0 !Stellar angular velocity
   double precision, parameter :: rmi=2.2d0, rmo=3.0d0, Tshk=0d0, Macc = 1d-7
   double precision, parameter :: year_to_sec = 3.154d7
   double precision, dimension(n_cells) :: T, rho, vx, vy, vz, V, ne
   double precision ::  OmegasK, Rstar, Mstar, thetao, thetai, Lr, Tring, Sr
   double precision :: vp, y, cosTheta, rcyl, z, r, phi, Mdot, sinTheta, Rm, theta
   
   !Init parameters
   T = 0d0
   rho = 0d0
   V = 0d0
   ne = 0d0
   vx = 0d0; vy=0d0; vz=0d0
   !
   
   !Define some useful quantities
   !write(*,*) "Mstar ", etoile(1)%M, " Rstar ", etoile(1)%r*AU_to_Rsun
   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot = Macc * Msun_to_kg / year_to_sec !kg/s
   Lr = Ggrav * Mstar * Mdot / Rstar * (1. - 2./(rmi + rmo))
   !sqrt(2*GMstar/rstar)
   OmegasK = dsqrt(Mstar * Ggrav * 2d0 / Rstar)
   !write(*,*) "Vmax (km/s)", OmegasK/1d3
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
   write(*,*) "Ring T (K) ", Tring, 100*Sr, Lr/Lsun
   
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

       r = sqrt(z**2 + rcyl**2)
       !from trigonometric relations we get y using rcyln 
       !and the radius at r(rcyln, zn)
       sinTheta = rcyl/r
       cosTheta = z/r
       theta = asin(sinTheta)
       
       !All lines are between the spheres of radii 1d0 and rmo
       !if (rn < 1d0 .or. rn > rmo*unit) CYCLE
       !from definition of y
       y = sinTheta**2
       !rn is a field line radial distance (from the centre of the model to the center
       !of the cell.
       !if rn**3/rcyln**2 or Rm=rn/sint**2 is between rmi and rmo it means that the
       !cell icell (rn, theta) belongs to a field line.
       Rm = r**3 / rcyl**2 / etoile(1)%r
       !if (Rm < rmi .or. Rm >rmo .or. rn < 1) CYCLE
       if (Rm>=rmi .and. Rm<=rmo .and. r >= etoile(1)%r) then
          if (y>1) write(*,*) "Warning, y > 1", y, icell
          !Density
          rho(icell) = Mdot * Rstar / (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                      (AU_to_m * r)**(-5./2.) / dsqrt(2d0 * Ggrav * Mstar) * &
                      dsqrt((4d0-3d0*y)/(1d0-min(y,0.99))) * 1d3/masseH  !kg/m3 ->/m3
          !Velocities
          vp = OmegasK * dsqrt(max(etoile(1)%r/r - 1./Rm,0d0))
          vp = -vp * cosTheta * 3d0 / dsqrt(4d0-3d0*y)
          vx(icell) = vp * sinTheta !R
          vy(icell) = 0.0 !stellar rotation velocity
          vz(icell) = vp * cosTheta !z
          if (theta > PI/2.) vz(icell) = -vz(icell)
          V(icell) = dsqrt(vx(icell)**2 + vy(icell)**2 + vz(icell)**2)
          write(*,*) theta*rad_to_deg, vx(icell), V(icell)*sinTheta
       end if
      end do
     end do
    end do
    do icell=1,n_cells
     !Temperature
     if (rho(icell) > 0) &
      T(icell) = Tmax !max(5000.d0,Tmax-2500.*(rho(icell) - minval(rho,mask=rho>0))/&
        !(maxval(rho) - minval(rho,mask=rho>0)))
    end do

   CALL init_atomic_atmos(n_cells, T, ne, rho)
   allocate(atmos%Vmap(n_cells), atmos%vx(n_cells), atmos%vy(n_cells),atmos%vz(n_cells))
   atmos%Vmap = V
   atmos%vz = vz
   atmos%vy = vy
   atmos%vx = vx
   
   if (lstatic) atmos%moving = .false.
   lkeplerian = .false.
   atmos%v_char = maxval(atmos%vturb)
   if (atmos%moving) then
    atmos%v_char = atmos%v_char + maxval(abs(V))
    write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
   else
    write(*,*) " >-< Model is static"
   end if
   
   write(*,*) log(1e-6*maxval(rho)), log(1e-6*minval(rho,mask=rho>0))
   write(*,*) maxval(V)/1d3, minval(V,mask=V>0)/1d3
   
  RETURN
  END SUBROUTINE magneto_accretion_model
  
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid
   ! Values width idk are from a Solar model atmosphere by 
   ! Fontenla et al. namely the FAL C model.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne
   double precision :: Vkepmax
   !idk = 10
!    nHtot =  2.27414200581936d16
!    T = 45420d0
!    ne = 2.523785d16

   !idk = 75
   T=7590d0
   ne = 4.446967d20
   nHtot = 1.259814d23

    !idk = 81   
!     T=9400d0
!     ne = 3.831726d21
!     nHtot = 1.326625d23

    !idk = 0 
!     T = 100000d0
!     ne = 1.251891d16
!     nHtot = 1.045714d16
    
   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving = .false. !force to be static for this case
   !atmos%vturb = 0d0
   !tmos%vturb = 9.506225d3 !m/s !idk=10
   atmos%vturb = 1.696164d3 !idk = 75
   !atmos%vturb = 1.806787d3 !idk=81
   !atmos%vturb = 10.680960d3 !idk=0
   atmos%v_char = maxval(atmos%vturb)
   write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE prop_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to dust density and temperature.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne
   double precision :: Vkepmax
   
   nHtot = 1d14! * densite_gaz/MAXVAL(densite_gaz) + 1d12
   T = 3000d0!+Tdust
   ne = 0d0
   
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   if (lstatic) atmos%moving = .false.
   lkeplerian = .true.
   atmos%v_char = maxval(atmos%vturb)
   if (atmos%moving) then
    Vkepmax = dsqrt(Ggrav * sum(etoile%M) * Msun_to_kg * Rmin**2 / &
                ((Rmin**2 + minval(z_grid)**2)**1.5 * AU_to_m) )

    atmos%v_char = atmos%v_char + Vkepmax
    write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
   else
    write(*,*) " >-< Model is static"
   end if

  RETURN
  END SUBROUTINE prop_law_model

 END MODULE simple_models