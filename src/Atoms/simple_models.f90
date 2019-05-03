MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos, define_atomRT_domain, &
 						free_atomic_atmos, init_magnetic_field
 
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
   use TTauri_module, only : TTauri_temperature
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Tmax = 0d3, days_to_sec = 86400d0, Prot = 8. !days
   double precision, parameter :: rmi=2.2d0, rmo=3.0d0, Tshk=0d4, Macc = 1d-7
   double precision, parameter :: year_to_sec = 3.154d7, r0 = 1d0, deltaz = 0.2!Rstar
   double precision ::  OmegasK, Rstar, Mstar, thetao, thetai, Lr, Tring, Sr, Q0, nH0
   double precision :: vp, y, rcyl, z, r, phi, Mdot, sinTheta, Rm, L
   double precision :: Omega, Vphi, TL(8), Lambda(8), rho_to_nH !K and erg/cm3/s
   double precision :: Bp, BR, Bz
   logical :: lwrite_model_ascii = .false.
   !wind 
   double precision :: Mdotwind=1d-9, Twind=8d3, beta=5d-1, vinf=400d0, vinit10d0
   
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 / !log10 ?
   !Lambda = Q/nH**2
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 / !log10? 
   
   Omega = 2.*PI / (Prot * days_to_sec) ! Angular velocity in rad/s

   linfall = .false.
   lkeplerian = .false.
   lmagnetoaccr = .not.(lwrite_model_ascii)
   CALL init_atomic_atmos()
   !For now we do not necessarily need to compute B if no polarization at all
   atmos%magnetized = (lmagnetic_field) .and. (PRT_SOLUTION /= "NO_STOKES")
   if (atmos%magnetized) CALL init_magnetic_field()
   
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
   !1G = 1e-4 T
   Bp = 1d-4 * 4.2*1d2 * (rmi/2.2)*(2*Mstar*kg_to_Msun)**(0.25) * &
   	    (Macc * 1d8)**(0.5) * (Rstar/(2*Rsun))**(-3.)
   write(*,*) "Equatorial magnetic field", Bp*1d4, 'G'
  
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
       y = min(y,0.99)
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
           atmos%nHtot(icell) = 1d23!
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
           ! Now magnetic field
            BR = Bp * (etoile(1)%r / r)**(3.0) * 3*dsqrt(y)*dsqrt(1.-y**2)
            Bz = Bp * (etoile(1)%r / r)**(3.0) * (3*(1.-y)-1)
            if (atmos%magnetized) then
             atmos%Bxyz(icell,1) = BR
             atmos%Bxyz(icell,2) = Bz
             atmos%Bxyz(icell,3) = 0d0
              if (.not.lmagnetoaccr) then !projection on the cell
               atmos%Bxyz(icell,1) = BR*cos(phi)
               atmos%Bxyz(icell,2) = BR*sin(phi)
               atmos%Bxyz(icell,3) = Bz
               if (z<0) atmos%Bxyz(icell,3) = -atmos%Bxyz(icell,3)
              end if
            end if
           end if
          L = 10 * Q0*(r0*etoile(1)%r/r)**3 / atmos%nHtot(icell)**2!erg/cm3/s
          !atmos%T(icell) = 10**(interp1D(Lambda, TL, log10(L)))
          atmos%T(icell) = 10**(interp_dp(TL, Lambda, log10(L)))
          !write(*,*) log10(L), Rm, rcyl/etoile(1)%r, atmos%T(icell), nH0, Q0
          
         !! Stellar wind 
!        else if (r/Rstar>=3d0 .and. y <= sin(50d0*PI/180d0)**2) then !wind
!            if (.not.lstatic) then
!              vp = (vinit + (vinf - vinit) * (1d0 - 3d0*Rstar/r)**beta) * 1d3
!              atmos%Vxyz(icell,1) = vp * dsqrt(y) !Rhat
!              atmos%Vxyz(icell,2) = vp * dsqrt(1.-y) !zhat
!              atmos%Vxyz(icell,3) = 0d0 !phihat
!              if (.not.lmagnetoaccr) then !projection
!               atmos%Vxyz(icell,3) = vp * dsqrt(1.-y) !zhat         
!               atmos%Vxyz(icell,1) = vp * 3d0 * dsqrt(y) * cos(phi) 
!               atmos%Vxyz(icell,2) = vp * 3d0 * dsqrt(y) * sin(phi) 
!              if (z < 0) atmos%Vxyz(icell,3) = -atmos%Vxyz(icell,3)
!              end if
!            end if
!           atmos%nHtot(icell) = Mdotwind * Msun_to_kg / year_to_sec / (4d0*PI*vp*(1.-cos(50d0*PI/180d0)**2))) / &
!            (r*AU_to_m)**2 * rho_to_nH
       end if
      end do
     end do
    end do
    !The disk could be included as new zone or by extending the density to r>rmo
    !but it has keplerian rotation and not infall but mcfost cannot treat two
    !kind of velocities right?
    
    if (Tmax > 0d0) atmos%T(:) = Tmax * atmos%T/maxval(atmos%T)
    !atmos%T(:) = minval(atmos%nHtot)/atmos%nHtot * Tmax

   !! inclued in atom%vbroad so we do not count it here
   !!if (maxval(atmos%vturb) > 0) atmos%v_char = minval(atmos%vturb,mask=atmos%vturb>0)
   if (.not.lstatic) then
    atmos%v_char = atmos%v_char + &
    minval(dsqrt(sum(atmos%Vxyz**2,dim=2)),&
    	   dim=1,mask=sum(atmos%Vxyz**2,dim=2)>0)
   end if
   
   if (atmos%magnetized) then
    atmos%B_char = maxval(sum(atmos%Bxyz(:,:)**2,dim=2))
    write(*,*)  "Typical Magnetic field modulus (G)", atmos%B_char * 1d4
   end if

!    atmos%T = 0d0
!    CALL TTauri_Temperature(rmi, rmo, Macc)
   CALL define_atomRT_domain()
   
   if (lwrite_model_ascii.and..not.lmagnetoaccr) then !for now
    CALL writeAtmos_ascii()
    !leave here, we enter here only to write the model for Voronoi tesselation
    CALL free_atomic_atmos()
    !densities will be fill again during the Voronoi tesselation
    !leave here with an empty atomic atmos that will be filled during the tesselation
    !using this model
     write(*,*) "Stop after writing for the moment"
     stop
   end if
   
   if (.not.lstatic) then
    write(*,*) "Maximum/minimum velocities in the model (km/s):"
    write(*,*) maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1)*1d-3,  &
     		  minval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1,&
     		  mask=sum(atmos%Vxyz**2,dim=2)>0)*1d-3
   end if
   write(*,*) "Typical velocity in the model (km/s):"
   write(*,*) atmos%v_char/1d3

   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%lcompute_atomRT==.true.)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%lcompute_atomRT==.true.)  
   
  RETURN
  END SUBROUTINE magneto_accretion_model
  
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid.
   ! Values width idk are from a Solar model atmosphere by 
   ! Fontenla et al. namely the FAL C model.
   ! idk = 0, top of the atmosphere, idk = 81 (max) bottom.
  ! ----------------------------------------------------------------- !
   INTEGER :: i, j, k, icell
   double precision :: r, rcyl, z

   CALL init_atomic_atmos()
   lstatic = .true. !force to be static for this case
   atmos%magnetized = .false.
   atmos%calc_ne = .false.

! 
!     do i=1, n_rad
!      do j=j_start,nz !j_start = -nz in 3D
!       do k=1, n_az
!        if (j==0) then !midplane
!         icell = cell_map(i,1,k)
!         rcyl = r_grid(icell) !AU
!         z = 0.0_dp
!        else
!         icell = cell_map(i,j,k)
!         rcyl = r_grid(icell)
!         z = z_grid(icell)/z_scaling_env
!        end if
!        r = dsqrt(z**2 + rcyl**2)
!        
!        if (r<=Rmax/2) then
! 
!          atmos%nHtot(icell) =  2.27414200581936d16
!          atmos%T(icell) = 45420d0
!          atmos%ne(icell) = 2.523785d16
!         atmos%vturb(icell) = 9.506225d3 !m/s
!         
!        end if
!       end do
!      end do
!     end do

   !idk = 10
    atmos%nHtot =  2.27414200581936d16
    atmos%T = 45420d0
    atmos%ne = 2.523785d16
    atmos%vturb = 9.506225d3 !m/s
   
   !idk = 75
!   atmos%T=7590d0
!   atmos%ne = 4.446967d20
!   atmos%nHtot = 1.259814d23

    !idk = 81   
!     atmos%T=9400d0
!     atmos%ne = 3.831726d21
!     atmos%nHtot = 1.326625d23
!     atmos%vturb = 1.806787d3 !idk=81

    !idk = 0 
!      atmos%T = 100000d0
!      atmos%ne = 1.251891d16
!      atmos%nHtot = 1.045714d16 
!      atmos%vturb = 10.680960d3 !idk=0
   CALL define_atomRT_domain()
   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%lcompute_atomRT==.true.)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%lcompute_atomRT==.true.) 
   write(*,*) "Maximum/minimum Electron density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%ne),MINVAL(atmos%ne,mask=atmos%lcompute_atomRT==.true.)

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE spherical_shells_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with moving spherical shells
  ! ----------------------------------------------------------- !
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Mdot = 1d-5, Vexp = 200d3, beta = 5d-1, nH0 = 1d14
   double precision, parameter :: year_to_sec = 3.154d7, r0 = 1d0, v0 = 0d0
   double precision :: Rstar, Mstar, Vr, rhor, rho_to_nH
   double precision :: rcyl, z, r, phi, Mdot_si
   


   linfall = .true.
   lkeplerian = .false.
   lmagnetoaccr = .false.
   lstatic = .false.
   CALL init_atomic_atmos()
   atmos%vturb = 1d3 !m/s
   rho_to_nH = 1d3/masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3

   
   if (linfall) then 
   	if (.not.allocated(Vfield)) allocate(Vfield(n_cells))
   	deallocate(atmos%Vxyz)
   end if

   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot_si = Mdot * Msun_to_kg / year_to_sec !kg/s
   
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
       
       Vr = v0 + (Vexp - v0) * (1d0 - (etoile(1)%r*r0)/r)**beta 
       rhor = Mdot_si / 4d0 / PI / ((r*AU_to_m)**2 * Vr) !kg/m3
       atmos%nHtot(icell) = rhor * rho_to_nH
       atmos%T(icell) = (r0*etoile(1)%r/r)**2 * 3d3

       !!->expansion
        if (.not.linfall) then !expansion not analytic
         atmos%Vxyz(icell,1) = Vr * z/r * cos(phi)  
         atmos%Vxyz(icell,2) = Vr * z/r * sin(phi)
         atmos%Vxyz(icell,3) = Vr * rcyl/r
        else!expansion analytic
         Vfield(icell) = Vr
        end if

      end do
     end do
    end do


    if (.not.linfall) then
    	atmos%v_char = atmos%v_char + &
    	minval(dsqrt(sum(atmos%Vxyz**2,dim=2)),&
    	   dim=1,mask=sum(atmos%Vxyz**2,dim=2)>0)
    else
        atmos%v_char = atmos%v_char + minval(Vfield)
    end if


   CALL define_atomRT_domain()
   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%lcompute_atomRT==.true.)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%lcompute_atomRT==.true.) 

  RETURN
  END SUBROUTINE spherical_shells_model
  
  SUBROUTINE rectangular_slab_model()
  
  
  RETURN
  END SUBROUTINE rectangular_slab_model
  
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
   double precision :: dx, dy, dz, smooth_scale
   
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
       ZZ(icell) = z*AU_to_m
     end do
    end do
   end do
   
   if (lstatic) then 
    allocate(atmos%Vxyz(n_cells,3))
    atmos%Vxyz = 0d0
   end if
   
   open(unit=1,file=trim(filename),status="replace")
   atmos%nHtot = atmos%nHtot / rho_to_nH !using total density instead of nHtot for writing
   !write Masses per cell instead of densities. Voronoi cells Volume is different than
   !grid cell right ? 
   write(*,*) "Maximum/minimum density in the model (kg/m^3):"
   write(*,*) maxval(atmos%nHtot), minval(atmos%nHtot,mask=atmos%lcompute_atomRT==.true.)
   atmos%nHtot = atmos%nHtot * volume * AU3_to_m3  * kg_to_Msun!Msun
   do icell=1, atmos%Nspace
     dx = XX(icell) - etoile(1)%x*AU_to_m
     dy = YY(icell) - etoile(1)%y*AU_to_m
     dz = ZZ(icell) - etoile(1)%z*AU_to_m
     smooth_scale = 5. * Rmax * AU_to_m
     if (min(dx,dy,dz) < 2*etoile(1)%r*AU_to_m) then
      if ((dx**2+dy**2+dz**2) < (2*etoile(1)%r*AU_to_m)**2) then
         smooth_scale = dsqrt(dx**2+dy**2+dz**2)/3. * AU_to_m
         !=etoile(1)%r*AU_to_m!etoile(1)%r/3. * AU_to_m!m
      end if
     end if
        write(1,"(8E)") XX(icell), YY(icell), ZZ(icell), atmos%nHtot(icell), &
    	 atmos%Vxyz(icell,1), atmos%Vxyz(icell,2), atmos%Vxyz(icell,3), smooth_scale
   end do
   write(*,*) "Maximum/minimum mass in the model (Msun):"
   write(*,*) maxval(atmos%nHtot), minval(atmos%nHtot,mask=atmos%lcompute_atomRT==.true.)
   write(*,*) "Maximum/minimum mass in the model (Kg):"
   write(*,*) maxval(atmos%nHtot)*Msun_to_kg, &
   				minval(atmos%nHtot,mask=atmos%lcompute_atomRT==.true.)*Msun_to_kg
   write(*,*) "Maximum/minimum velocities in the model (km/s):"
   write(*,*) maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1)*1d-3,  &
     		  minval(dsqrt(sum(atmos%Vxyz**2,dim=2)), dim=1,&
     		  mask=sum(atmos%Vxyz**2,dim=2)>0)*1d-3
   close(unit=1)

  RETURN
  END SUBROUTINE writeAtmos_ascii

 END MODULE simple_models