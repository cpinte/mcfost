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
   use math, only : locate
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Tmax = 8d3, Omegas = 0d0 !Stellar angular velocity
   double precision, parameter :: rmi=2.2d0, rmo=3.0d0, Tshk=0d0, Macc = 1d-7
   double precision, parameter :: Rmmax = 3.5d0, Rmmin = 0d0, year_to_sec = 3.154d7, Zmmax=2.
   double precision, dimension(n_cells) :: T, rho, vx, vy, vz, V, ne
   double precision ::  OmegasK, Rstar, Mstar, thetao, thetai, Lr, Tring, unit=0d0
   double precision :: vp, y, cosTheta, rcyl, z, r, Mdot, sinTheta, rcyln, zn, rn, Rm
   double precision :: TL(8), Lambda(8) !K and erg/cm3/s
   
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 / !log
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 / !log
   !Lambda = Prad/nH**2
   ! 1e-7 * 1e6 = ergcm^-3 to Jm^3
   
   !Init parameters
   T = 0d0
   rho = 0d0
   V = 0d0
   ne = 0d0
   !
   
   !Define some useful quantities
   !write(*,*) "Mstar ", etoile(1)%M, " Rstar ", etoile(1)%r*AU_to_Rsun
   Rstar = etoile(1)%r * AU_to_Rsun * Rsun !m
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
   write(*,*) "Ring T (K) ", Tring!, Lr/Lsun

   main : do izone=1, n_zones
    !dz = disk_zone(izone) ! for now the parameters are hardcoded
    rcyl_loop : do i=1, n_rad
     z_loop : do j=j_start,nz
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


       if (unit <= 1d0) then
        zn = (z - minval(z_grid/z_scaling_env)) / &
         (maxval(z_grid/z_scaling_env)-minval(z_grid/z_scaling_env)) * Zmmax
        rcyln = (rcyl-Rmin) * (Rmmax - Rmmin)/ (Rmax - Rmin) + Rmmin
        unit = 1d0
       else
        zn = z/unit
        rcyln = rcyl/unit !in Rstar, assuming Rstar=unit
       end if
       rn = sqrt(zn**2 + rcyln**2)
       
       !All lines are between the spheres of radii 1d0 and rmo
       !if (rn < 1d0 .or. rn > rmo*unit) CYCLE
       
       !from trigonometric relations we get y using rcyln 
       !and the radius at r(rcyln, zn)
       cosTheta = zn/rn
       sinTheta = rcyln/rn
       !from definition of y
       y = sinTheta**2
       
       !rn is a field line radial distance (from the centre of the model to the center
       !of the cell.
       !if rn**3/rcyln**2 or Rm=rn/sint**2 is between rmi and rmo it means that the
       !cell icell (rn, theta) belongs to a field line.
       Rm = rn**3 / rcyln**2
       if (Rm < rmi .or. Rm >rmo) CYCLE
       !if (Rm>=rmi .and. Rm<=rmo) then
          if (y>1) write(*,*) "Warning, y > 1", y, icell
          if (y==1d0) write(*,*) "Warning, y == 1", y, icell
          if (y==1d0) y=y-1d-16 !avoiding /1-1
          !Density
          rho(icell) = Mdot / (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                      (Rstar*rn)**(-2.5d0) / dsqrt(2d0 * Ggrav * Mstar) * &
                      dsqrt((4d0-3d0*y)/(1d0-y)) * 1d3/masseH  !kg/m3 ->/m3
          !Velocities
          vp = OmegasK * dsqrt(1./rn - 1./Rm) 
          vx(icell) = -vp * 3d0*dsqrt(y)*cosTheta/dsqrt(4d0-3d0*y)
          vy(icell) = 0.0 !stellar rotation velocity
          vz(icell) = -vp * 3d0*cosTheta**2/dsqrt(4d0-3d0*y)
          V(icell) = dsqrt(vx(icell)**2 + vy(icell)**2 + vz(icell)**2)
          !Temperature
          T(icell) = Tmax * (1.5 / rn)**3
       !end if
      end do !k
     end do z_loop!j
    end do rcyl_loop !i 
   end do main !over n_zones 
   
   rho = rho/1.28876812421607  !move in init_atomic_atmos
   
   write(*,*)Rmin, Rmax, minval(z_grid/z_scaling_env), maxval(z_grid/z_scaling_env)
   CALL writeVfield(V)
   CALL writeTfield(T)
   CALL init_atomic_atmos(n_cells, T, ne, rho)
   
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
   atmos%velocity_law = 0 !keplerian = 0
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
   atmos%velocity_law = 0 !keplerian
   atmos%v_char = maxval(atmos%vturb)
   if (atmos%moving) then
    Vkepmax = dsqrt(Ggrav * sum(etoile%M) * Msun_to_kg * Rmin**2 / &
                ((Rmin**2 + minval(z_grid)**2)**1.5 * AU_to_m) )

    atmos%v_char = atmos%v_char + Vkepmax
    write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
    if (atmos%velocity_law == 0) then
      if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
      atmos%Vmap = 0d0
    end if
   else
    write(*,*) " >-< Model is static"
   end if

  RETURN
  END SUBROUTINE prop_law_model
  
 SUBROUTINE writeVfield(Vmap)
  double precision, intent(in), dimension(n_cells) :: Vmap
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,trim("Vfield.fits.gz"),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 !??
  fpixel = 1
  extend = .false. !??
  bitpix = -64  

  if (lVoronoi) then
   naxis = 1
   naxes(1) = n_cells
   nelements = naxes(1)
  else
   if (l3D) then
    naxis = 3
    naxes(1) = n_rad
    naxes(2) = 2*nz
    naxes(3) = n_az
    nelements = naxes(1) * naxes(2) * naxes(3) ! != n_cells ? should be
   else
    naxis = 2
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2) ! should be n_cells also
   end if
  end if

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "m/s", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,Vmap,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writeVfield
 
 SUBROUTINE writeTfield(Tfield)
  double precision, intent(in), dimension(n_cells) :: Tfield
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,trim("Tfield.fits.gz"),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 !??
  fpixel = 1
  extend = .false. !??
  bitpix = -64  

  if (lVoronoi) then
   naxis = 1
   naxes(1) = n_cells
   nelements = naxes(1)
  else
   if (l3D) then
    naxis = 3
    naxes(1) = n_rad
    naxes(2) = 2*nz
    naxes(3) = n_az
    nelements = naxes(1) * naxes(2) * naxes(3) ! != n_cells ? should be
   else
    naxis = 2
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2) ! should be n_cells also
   end if
  end if

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "K", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,Tfield,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writeTfield

 END MODULE simple_models