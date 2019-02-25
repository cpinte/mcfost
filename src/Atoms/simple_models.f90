MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos, define_atomRT_domain, &
 						free_atomic_atmos
 
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
  ! ----------------------------------------------------------- !
   use math, only : interp1D
   use utils, only : interp_dp
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Tmax = 8d3, days_to_sec = 86400d0, Prot = 8. !days
   double precision, parameter :: rmi=2.2d0, rmo=3.0d0, Tshk=7d3, Macc = 1d-7
   double precision, parameter :: year_to_sec = 3.154d7, r0 = 1d0, deltaz = 0.1!Rstar
   double precision ::  OmegasK, Rstar, Mstar, thetao, thetai, Lr, Tring, Sr, Q0, nH0
   double precision :: vp, y, rcyl, z, r, phi, Mdot, sinTheta, Rm, L
   double precision :: Omega, Vphi, NormV(n_cells), TL(8), Lambda(8) !K and erg/cm3/s
   double precision :: rho_to_nH
   
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 / !log10 ?
   !Lambda = Q/nH**2
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 / !log10? 
   
   Omega = 2.*PI / (Prot * days_to_sec) ! Angular velocity in rad/s

   linfall = .false.
   lkeplerian = .false.
   lmagnetoaccr = .false.
   CALL init_atomic_atmos()
   atmos%vturb = 0d3 !m/s
   !lvoronoi = .true.!has to be done after init with lvoronoi=.false.
   !otherwise atmos%Vxyz is not allocated
   
   rho_to_nH = 1d3/masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3

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
       Vphi = Omega * (r*AU_to_m) !m/s
       if (.not.lstatic.and.lmagnetoaccr) atmos%Vxyz(icell,3) = Vphi !phihat
       !from trigonometric relations we get y using rcyln 
       !and the radius at r(rcyln, zn)
       sinTheta = rcyl/r
       ! and cosTheta = z/r = sqrt(1.-y) = sqrt(1.-sinTheta**2)
       !from definition of y
       y = sinTheta**2!min(real(sinTheta**2), 0.9999)
       !here the trick is to handle y = 1 at z=0 because rcyl=r and 1/(1-y) = inf.
       !r is a field line radial distance (from the centre of the model to the centre
       !of the cell.
       !if r**3/rcyl**2 or Rm=r/sint**2 is between rmi and rmo it means that the
       !cell icell (r, theta) belongs to a field line.
       Rm = r**3 / rcyl**2 / etoile(1)%r !in Rstar, same as rmi,o
       if (Rm>=rmi .and. Rm<=rmo) then !r>etoile(1)%r  
          nH0 = rho_to_nH * (Mdot * Rstar) /  (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (Rstar * r0)**( real(-5./2.) ) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*(r0/Rm)) / dsqrt(1d0-(r0/Rm))
          Q0 = nH0**2 * 10**(interp_dp(Lambda, TL, log10(Tring)))*0.1
          !write(*,*) icell, Rm, nH0, Q0, log10(Q0/nH0**2)
          !Interface funnels/disk
          if (dabs(z*AU_to_m/Rstar)<=deltaZ) then!(z==0d0) then
           !at z=0, r = rcyl = Rm here so 1/r - 1/Rm = 0d0 and y=1: midplane
           vp = 0d0
           !atmos%Vxyz is initialized to 0 everywhere
           if (.not.lstatic.and..not.lmagnetoaccr) then!only project if we are not
           													!using analytical velocity law
           													!hence if atmos%magnetoaccr = .false.
            !V . xhat = (0rhat,0thetahat,Omega*r phihat) dot (..rhat,..thetahat,-sin(phi)phihat)
            !V . yhat = (0rhat,0thetahat,Omega*r phihat) dot (..rhat,..thetahat,cos(phi)phihat)
            !V .zhat = (0rhat,0thetahat,Omega*r phihat) dot (..rhat,..thetahat,0phihat)
            atmos%Vxyz(icell,1) = Vphi * -sin(phi) !xhat
            atmos%Vxyz(icell,2) = Vphi *  cos(phi) !yhat
            !here it is better to have keplerian rotation than stellar rotation?
           end if
           !Density midplane of the disk
           !using power law of the dust-gas disk (defined in another zone?)
           atmos%nHtot(icell) = 1d19!
           !Temperature still given by the same law as in accretion funnels yet.
           !But probably different at this interface funnels/disk
          else !accretion funnels
          !Density
           atmos%nHtot(icell) = ( Mdot * Rstar )/ (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (AU_to_m * r)**( real(-2.5) ) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*y) / dsqrt(1d0-y) * &
                       rho_to_nH  !kg/m3 ->/m3
           vp = OmegasK * dsqrt(etoile(1)%r/r - 1./Rm)
           vp = -vp / dsqrt(4d0-3d0*y)
           if (.not.lstatic) then
             atmos%Vxyz(icell,1) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) !Rhat
             atmos%Vxyz(icell,2) = vp * (2d0 - 3d0 * y) !zhat
             atmos%Vxyz(icell,3) = Vphi !phihat
             if (.not.lmagnetoaccr) then !projection
              atmos%Vxyz(icell,3) = atmos%Vxyz(icell,2) !zhat         
              atmos%Vxyz(icell,1) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) * cos(phi) &
              						-Vphi*sin(phi)!xhat
              atmos%Vxyz(icell,2) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) * sin(phi) &
              						+Vphi*cos(phi)!yhat
             if (z < 0) atmos%Vxyz(icell,3) = -atmos%Vxyz(icell,3)
             end if
             ! or theta > PI/2d0
             !Because z, is the second axis if Vxyz=(Vr, Vz, Vphi) and the third if the
             !velocity is projected on the cartesian coordinates
             !!the minus sign in case of lmagnetoaccr is now taken in v_proj
             !!if (z<0 .and.lmagnetoaccr) atmos%Vxyz(icell,2) = -atmos%Vxyz(icell,2)

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

   !! inclued in atom%vbroad so we do not count it here
   !!if (maxval(atmos%vturb) > 0) atmos%v_char = minval(atmos%vturb,mask=atmos%vturb>0)
   if (.not.lstatic) then
    NormV = dsqrt(atmos%Vxyz(:,1)**2+atmos%Vxyz(:,2)**2+atmos%Vxyz(:,3)**2)
    atmos%v_char = atmos%v_char + minval(NormV,mask=NormV>0)
   end if
   write(*,*) "maxVfield (km/s)", maxval(NormV)/1d3, &
              " minVfield", minval(NormV,mask=normV>0)/1d3
   write(*,*) "atmos%v_char (km/s) =", atmos%v_char/1d3
   CALL define_atomRT_domain()
   
   if (.not.lmagnetoaccr) then !for now
    CALL writeAtmos_ascii()
    !leave here, we enter here only to write the model for Voronoi tesselation
    CALL free_atomic_atmos()
    !densities will be fill again during the Voronoi tesselation
    !leave here with an empty atomic atmos that will be filled during the tesselation
    !using this model
    write(*,*) "Stop after writing for the moment"
    stop
   end if
   
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
   CALL init_atomic_atmos()!(n_cells)

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
   atmos%v_char = 0d0!would be only atom%vbroad which contains vturb already

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE spherical_shells_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with moving spherical shells
  ! ----------------------------------------------------------- !
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Tmax = 3d3, days_to_sec = 86400d0, Prot = 8. !days
   double precision, parameter :: Macc = 1d-7, Vexp = 200d3
   double precision, parameter :: year_to_sec = 3.154d7, r0 = 1d0
   double precision :: Rstar, Mstar, normV(n_Cells), Omega
   double precision :: rcyl, z, r, phi, Mdot
   


   linfall = .true.
   lkeplerian = .false.
   lmagnetoaccr = .false.
   CALL init_atomic_atmos()
   atmos%vturb = 0d3 !m/s
   
   if (linfall.or.lkeplerian) then 
   	if (.not.allocated(Vfield)) allocate(Vfield(n_cells))
   	deallocate(atmos%Vxyz)
   end if

   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot = Macc * Msun_to_kg / year_to_sec !kg/s
   
   Omega = 2.*PI / (Prot * days_to_sec) ! Angular velocity in rad/s
   
   atmos%nHtot = 1d14
   atmos%T = Tmax

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
       
       if (.not.lstatic) then
       !!->expansion
        if (.not.linfall.and..not.lmagnetoaccr) then !expansion not analytic
         atmos%Vxyz(icell,1) = Vexp * z/r * cos(phi)
         atmos%Vxyz(icell, 2) = Vexp * z/r * sin(phi)
         atmos%Vxyz(icell,3) = Vexp * rcyl/r
        else if (linfall.and..not.lmagnetoaccr) then !expansion analytic
         Vfield(icell) = Vexp
        else if (.not.linfall.and.lmagnetoaccr) then !rotation analytic
         atmos%Vxyz(icell,1) = 0d0
         atmos%Vxyz(icell,2) = 0d0 
         atmos%Vxyz(icell,3) = Omega * (r*AU_to_m)
        else if (lkeplerian) then !rotation analytic also
         Vfield(icell) = Omega * (r*AU_to_m)
        end if
       end if !static
      end do
     end do
    end do

   if (.not.lstatic) then
    if (linfall.or.lkeplerian) then
     NormV = Vfield
    else
     NormV = dsqrt(atmos%Vxyz(:,1)**2+atmos%Vxyz(:,2)**2+atmos%Vxyz(:,3)**2)
    end if
    atmos%v_char = atmos%v_char + minval(NormV,mask=NormV>0)
   end if
   write(*,*) "maxVfield (km/s)", maxval(NormV)/1d3, &
              " minVfield", minval(NormV,mask=normV>0)/1d3
   write(*,*) "atmos%v_char (km/s) =", atmos%v_char/1d3
   CALL define_atomRT_domain()

  RETURN
  END SUBROUTINE spherical_shells_model
  
  SUBROUTINE prop_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to molecular/dust density and temperature.
  ! ----------------------------------------------------------- !
   double precision :: Vkepmax, Vkepmin
   integer :: icell
   
   CALL init_atomic_atmos()!(n_cells)

   atmos%nHtot = 1d14 * densite_gaz/MAXVAL(densite_gaz) + 1d12
   atmos%T = 3000d0+Tdust
   atmos%ne = 0d0

   lkeplerian = .true.
   atmos%vturb = 1d0

   if (.not.lstatic) then
   !we need vfield
    if (.not.allocated(Vfield)) allocate(Vfield(n_cells))
    deallocate(atmos%Vxyz) !supposing lkep or linf, Vfield is used and veloc is varying across cells
    Vkepmax = dsqrt(Ggrav * sum(etoile%M) * Msun_to_kg * Rmin**2 / &
                 ((Rmin**2 + minval(z_grid)**2)**1.5 * AU_to_m) )
    if (lcylindrical_rotation) then ! Midplane Keplerian velocity
      do icell=1, n_cells
        Vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
      end do
    else ! dependance en z
       do icell=1, n_cells
         Vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg * r_grid(icell)**2 / &
                ((r_grid(icell)**2 + z_grid(icell)**2)**1.5 * AU_to_m) )
       end do
    end if
    Vkepmin = minval(Vfield, mask=Vfield>0)
    atmos%v_char = atmos%v_char + Vkepmin
   end if
   CALL define_atomRT_domain()

  RETURN
  END SUBROUTINE prop_law_model
  
  SUBROUTINE TTauri_Temperature()
  ! ---------------------------------------------------------------- !
   ! Computes temperature of accretion columns of T Tauri stars
   ! given the density, using the prescription of Hartmann 94
   ! and the cooling rates for Hartmann 82
   !
   ! TO DO: Numerically solve for the temperature, iterating
   ! between NLTE solution (Cooling rates) and Temperature->ne->npops
  ! ------------------------------------------------------------------ !
   use Voronoi_grid, only : Voronoi
   use utils, only : interp_dp
   integer 	:: icell, icell0
   double precision, parameter :: T0 = 7d3 !Defined the Temperature, must be consitant
   									!with the one in magneto_accretion_model() if used
   double precision :: Q0, nH0, L, TL(8), Lambda(8) !K and erg/cm3/s
   double precision :: r0 = 1d0, r
   
   !T = K
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 /
   !Lambda = Q/nH**2, Q in J/m9/s
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 /

   if (.not.lVoronoi) then
    write(*,*) "Only defined for external model read"
    stop
   end if


   icell0 = etoile(1)%icell - 1
   !We need to define a Q0... but we do not know the position along a field line
   !and cell, contruary to the case where the model is defined
   do icell=1, atmos%Nspace
    r = dsqrt(Voronoi(icell)%xyz(1)**2 +Voronoi(icell)%xyz(2)**2+Voronoi(icell)%xyz(3)**2 )
    Q0 = atmos%nHtot(icell0)**2 * 10**(interp_dp(Lambda, TL, log10(T0)))*0.1
    L = 10 * Q0*(r0*etoile(1)%r/r)**3 / atmos%nHtot(icell)**2!erg/cm3/s
    atmos%T(icell) = 10**(interp_dp(TL, Lambda, log10(L)))
    write(*,*) icell0, Q0, L, atmos%T(icell)
   end do

  RETURN
  END SUBROUTINE TTauri_Temperature
  
  SUBROUTINE writeAtmos_ascii()
  ! ------------------------------------- !
   ! Write GridType atmos to ascii file.
   ! ultimately this model will be mapped
   ! onto a Voronoi mesh.
  ! ------------------------------------- !
   character(len=10)	:: filename="model.s"
   integer :: icell, i, j, k
   double precision, dimension(n_cells) :: XX, YY, ZZ
   double precision :: r, rcyl, phi, z, sinTheta, rho_to_nH, fact
   
   rho_to_nH = 1d3/masseH /atmos%avgWeight
   
   if (n_az <= 1) then
    write(*,*) "Error, in this case, lmagnetoaccr is false, and", &
    	' n_az > 1'
    stop
   end if
   !X, Y, Z cartesian coordinates
   do i=1, n_rad
    do j=j_start,nz
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
      fact = 1d0
      if (z<0) fact = -1
      sinTheta = fact*dsqrt(1.-(z/r)**2)

   
        !XX(icell) = rcyl*cos(phi)*AU_to_m!r*cos(phi)*sinTheta * AU_to_m
        !YY(icell) = rcyl*sin(phi)*AU_to_m!r*sin(phi)*sinTheta * AU_to_m
        XX(icell) = r*cos(phi)*sinTheta * AU_to_m
        YY(icell) = r*sin(phi)*sinTheta * AU_to_m   
        ZZ(icell) = z
     end do
    end do
   end do
   
   open(unit=1,file=trim(filename),status="replace")
   atmos%nHtot = atmos%nHtot / rho_to_nH !using total density isntead of nHtot for writing
   do icell=1, atmos%Nspace
    write(1,"(7E)") XX(icell), YY(icell), ZZ(icell), atmos%nHtot(icell), &
    	 atmos%Vxyz(icell,1), atmos%Vxyz(icell,2), atmos%Vxyz(icell,3)
   end do
   
   close(unit=1)

  RETURN
  END SUBROUTINE writeAtmos_ascii

 END MODULE simple_models