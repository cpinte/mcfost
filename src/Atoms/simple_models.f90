MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos, define_atomRT_domain, &
 						free_atomic_atmos, init_magnetic_field, writeAtmos_ascii, &
 						writeVfield, write_atmos_domain
 
 ! MCFOST's modules
 use input
 use parametres
 use grid
 use density
 use constantes
 use writeAtom, only : writeTemperature, writeHydrogendensity, writeElectron
 use solvene, only : SolveElectronDensity
 use math, only : interp1D, locate
 use utils, only : interp_dp
 
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
   use TTauri_module, only : gett
   integer :: n_zones = 1, izone, i, j, k, icell, southern_hemp, idmax, niter
   integer, parameter :: NiterMax = 0!50
   double precision, parameter :: Tmax = 7.5d3, days_to_sec = 86400d0, Prot = 6.53 !days
   double precision, parameter :: rmi=2.2, rmo=3.0, Tshk=8d3, Macc = 1d-9, Tiso=0d3
   double precision, parameter :: yroot = (15d0 - dsqrt(33d0)) / 12d0 !value of y for which T is maximum (max(-nH**2 * r**3))
   double precision, parameter :: year_to_sec = 3.154d7, y_lim_z0 = 0.99!99
   double precision ::  OmegasK, Rstar, Mstar, thetao, thetai, Lr, Tring, Sr, Q0, nH0, fp
   double precision :: vp, y, rcyl, z, r, phi, Theta, Mdot, sinTheta, Rm, L, r0
   double precision :: Omega, Vphi, TL(8), Lambda(8), rho_to_nH !K and erg/cm3/s
   double precision :: Bp, BR, Bz, tc, phic, Tmax_old, Den, Mdotfid = 0d0
   logical, dimension(:), allocatable :: dark_zone
   logical, parameter :: lwrite_model_ascii = .false., accretion_spots = .true., include_dark_zone=.false.
   logical :: is_not_dark
   double precision :: zd_max = 0., dT, ne1
   !double precision, dimension(:), allocatable :: Told
      
   data TL / 3.70, 3.80, 3.90, 4.00, 4.20, 4.60, 4.90, 5.40 / !log10 ?
   !Lambda = Q/nH**2
   data Lambda / -28.3, -26., -24.5, -23.6, -22.6, -21.8,-21.2, -21.2 / !Q/nH**2
   
   Omega = 2.*PI / (Prot * days_to_sec) ! Angular velocity in rad/s

   linfall = .false.
   lkeplerian = .false.
   lmagnetoaccr = .not.(lwrite_model_ascii)
   lmagnetoaccr = .not.(lstatic)
   !!lmagnetoaccr = .false. !test
   CALL init_atomic_atmos()
   !For now we do not necessarily need to compute B if no polarization at all
   atmos%magnetized = (lmagnetic_field) .and. (PRT_SOLUTION /= "NO_STOKES")
   if (atmos%magnetized) CALL init_magnetic_field()
   
   rho_to_nH = 1d3/masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3
   fp = 1d0
   
   r0 = rmo * yroot
   
   is_not_dark = .true.
   if (include_dark_zone) then 
    allocate(dark_zone(n_cells)); dark_zone(:) = .false.
   end if

   !Define some useful quantities
   !write(*,*) "Mstar ", etoile(1)%M, " Rstar ", etoile(1)%r*AU_to_Rsun
   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot = Macc * Msun_to_kg / year_to_sec !kg/s
   Mdotfid = Mdotfid * Msun_to_kg / year_to_sec !kg/s
   Lr = Ggrav * Mstar * Mdot / Rstar * (1. - 2./(rmi + rmo))
   !sqrt(2*GMstar/rstar)
   OmegasK = dsqrt(Mstar * Ggrav * 2d0 / Rstar)
   !maximum and minimum angle on the stellar surface
   !all field lines are between theta0 and thetai at the stellar surface (r=Rstar)
   thetai = asin(dsqrt(1d0/rmi)) !rmi and rmo in Rstar
   thetao = asin(dsqrt(1d0/rmo))
   write(*,*) "Angles at stellar surface (deg)", thetao*rad_to_deg, thetai*rad_to_deg
   Tring = Lr / (4d0 * PI * Rstar**2 * sigma * abs(cos(thetai)-cos(thetao)))
   Tring = Tring**0.25
   
        !2 is for two rings, 2Pi is dphi from 0 to 2pi / total sphere surface
   Sr = abs(cos(thetai)-cos(thetao)) !2 * 2PI * Rstar**2 * abs(c1-c2) / 4pi Rstar**2
   write(*,*) "Ring T ", Tring, "K", " Shock T ", Tshk, ' K'
   write(*,*) "Surface ", 100*Sr, "%" !couverte par les deux anneaux sur toute la sphere
   write(*,*) "Luminosity", Lr/Lsun, "Lsun"
   write(*,*) "Mean molecular weight", atmos%avgWeight
   !1G = 1e-4 T
   Bp = 1d-4 * 4.2*1d2 * (rmi/2.2)*(2*Mstar*kg_to_Msun)**(0.25) * &
   	    (Macc * 1d8)**(0.5) * (Rstar/(2*Rsun))**(-3.)
   write(*,*) "Equatorial magnetic field", Bp*1d4, 'G'
 

   if (Tshk > 0) then 
    Tring = Tshk
    write(*,*) " Changing value for Tring to ", Tshk
   end if
 
!    write(*,*) Tring, Mstar, Rstar
!    write(*,*) masseH * 1d-3, Msun_to_kg, Ggrav
!    write(*,*) sigma, Mdot
 !     Surface on the apparent radius of one ring
!     etoile(1)%SurfB(k)%Sp = (etoile(1)%SurfB(k)%limits(4)-etoile(1)%SurfB(k)%limits(3)) *&
!     			abs(etoile(1)%SurfB(k)%limits(1)**2-etoile(1)%SurfB(k)%limits(2)**2) / (2d0 * PI)
!     etoile(1)%SurfB(k)%dOmega=dsqrt(1d0-etoile(1)%SurfB(k)%Sp)
!     write(*,*) " Surface on the disk of the spot ", k, etoile(1)%SurfB(k)%Sp, etoile(1)%SurfB(k)%dOmega, &
!     (etoile(1)%SurfB(k)%limits(4)-etoile(1)%SurfB(k)%limits(3)) *&
!     			abs(etoile(1)%SurfB(k)%limits(1)-etoile(1)%SurfB(k)%limits(2)) / (4d0 * PI)  
   !Add ring
   if (accretion_spots) then
   etoile(1)%Nr = 2
   	allocate(etoile(1)%SurfB(etoile(1)%Nr))
   	southern_hemp = 1
   	do k=1,etoile(1)%Nr!should be 1 ring for testing
    	tc = 0d0 * PI/180 !center of vector position
    	phic = 0d0 !vector position, pointing to the center of the spot
    	if (k==2) southern_hemp = -1 !for the ring in the southern hemisphere
    	etoile(1)%SurfB(k)%T = Tring
    !center of the spot
   	 etoile(1)%SurfB(k)%r(1) = cos(phic)*sin(tc)
   	 etoile(1)%SurfB(k)%r(2) = sin(phic)*sin(tc)
   	 etoile(1)%SurfB(k)%r(3) = cos(tc)*southern_hemp
    	etoile(1)%SurfB(k)%muo = cos(thetao); etoile(1)%SurfB(k)%mui = cos(thetai)
    	etoile(1)%SurfB(k)%phio = 2*PI; etoile(1)%SurfB(k)%phii = 0d0
  	 end do
   end if !accretion_spots

    do i=1, n_rad
     do j=j_start,nz !j_start = -nz in 3D
      do k=1, n_az
       if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
        if (include_dark_zone) dark_zone(icell) = (rcyl>=rmi*etoile(1)%r)
        sinTheta = 1d0
        y = 1d0
        rM = rcyl/etoile(1)%r
        r = rcyl
       else
        icell = cell_map(i,j,k)
        rcyl = r_grid(icell)
        z = z_grid(icell)/z_scaling_env
        r = dsqrt(z**2 + rcyl**2)
        sinTheta = rcyl/r
        y = sinTheta**2
        Rm = r**3 / rcyl**2 / etoile(1)%r !in Rstar, same as rmi,o
        if ((zd_max > 0 .and. abs(abs(z/etoile(1)%r)-zd_max) <= zd_max) .and.rcyl/etoile(1)%r >= rmi) then
         if (include_dark_zone) dark_zone(icell) = .true.
        end if
       end if
       phi = phi_grid(icell)
       

       Vphi = Omega * (r*AU_to_m) * dsqrt(y) !m/s *sinTheta, the pole doesn't rotate

       if (.not.lstatic.and.lmagnetoaccr) atmos%Vxyz(icell,3) = Vphi !phihat       

       Den = 1d0 / dsqrt(1d0 - min(y,y_lim_z0))
       y = min(y, y_lim_z0)

       !always true if not dark zone.
       if (include_dark_zone) is_not_dark = .not.dark_zone(icell)

       !if not dark, possibly matter. But if Rm is out of range, it is transparent
       if ((Rm>=rmi .and. Rm<=rmo).and.is_not_dark) then
          if (Mdotfid > 0) then
!            nH0 = rho_to_nH * (Mdotfid * Rstar) /  (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
!                        (Rstar * r0)**(-5./2.) / dsqrt(2d0 * Ggrav * Mstar) * &
!                        dsqrt(4d0-3d0*(r0/Rm)) / dsqrt(1d0-(r0/Rm))
           nH0 = rho_to_nH * (Mdotfid * Rstar) /  (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (Rstar * r0)**(-5./2.) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*yroot) / dsqrt(1d0-yroot)
          else
            nH0 = rho_to_nH * (Mdot * Rstar) /  (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (Rstar * r0)**(-5./2.) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*yroot) / dsqrt(1d0-yroot)         
		  end if
          Q0 = 10**(interp_dp(Lambda, TL, log10(Tmax)))*0.1 / fp * nH0**2 !/ dsqrt(4d0-3d0*yroot)

          
          
!           write(*,*) icell, "Rm=", Rm, "y=",y," r = ", r/etoile(1)%r
!           write(*,*) " logQ0=", log10(Q0), "rho0=", nH0/rho_to_nH, "den=", den
!           write(*,*) " nH0(log(cm^-3)) = ", log10(1d-6 * nH0)

           vp = OmegasK * dsqrt(etoile(1)%r/r - 1./Rm)
           vp = -vp / dsqrt(4d0-3d0*y)
           
           if (.not.lstatic) then
             atmos%Vxyz(icell,1) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) !Rhat
             atmos%Vxyz(icell,2) = vp * (2d0 - 3d0 * y) !zhat
             atmos%Vxyz(icell,3) = Vphi !phihat
             !if ( z < 0. ) atmos%Vxyz(icell,2) = -atmos%Vxyz(icell,2) !done in v_proj
             if (.not.lmagnetoaccr) then !projection
              atmos%Vxyz(icell,3) = atmos%Vxyz(icell,2) !zhat         
              atmos%Vxyz(icell,1) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) * cos(phi) &
              						-Vphi*sin(phi)!xhat
              atmos%Vxyz(icell,2) = vp * 3d0 * dsqrt(y) * dsqrt(1.-y) * sin(phi) &
              						+Vphi*cos(phi)!yhat
             if ( z < 0 ) atmos%Vxyz(icell,3) = -atmos%Vxyz(icell,3) !not done in v_proj
             end if
             ! or theta > PI/2d0
             !Because z, is the second axis if Vxyz=(Vr, Vz, Vphi) and the third if the
             !velocity is projected on the cartesian coordinates
             !!the minus sign in case of lmagnetoaccr is now taken in v_proj
             !!if (z<0 .and.lmagnetoaccr) atmos%Vxyz(icell,2) = -atmos%Vxyz(icell,2)

           end if !not static
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
            
          !Den = 1d0 / dsqrt(1d0 - min(y,0.99))
          atmos%nHtot(icell) = ( Mdot * Rstar )/ (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (AU_to_m * r)**(-5./2.) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*y) * Den * rho_to_nH  !kg/m3 ->/m3
                       
          if (Mdotfid > 0) then
          !Using T(nH) computed with nH fiducial
		   atmos%T(icell) = ( Mdotfid * Rstar )/ (4d0*PI*(1d0/rmi - 1d0/rmo)) * &
                       (AU_to_m * r)**(-5./2.) / dsqrt(2d0 * Ggrav * Mstar) * &
                       dsqrt(4d0-3d0*y) * Den * rho_to_nH
           L = fp * 10 * (r0*etoile(1)%r/r)**3 / atmos%T(icell)**2 * Q0 * dsqrt(4d0-3d0*y)!erg/cm3/s
          else
           L = fp * 10 * (r0*etoile(1)%r/r)**3 / atmos%nHtot(icell)**2 * Q0 !* dsqrt(4d0-3d0*y)!erg/cm3/s
          end if
          atmos%T(icell) = 10**(interp_dp(TL, Lambda, log10(L)))
!           write(*,*) "log(L)=", log10(L)
!           write(*,*) "log(Q)=",log10((r0*etoile(1)%r/r)**3*Q0)," rho=", atmos%nHtot(icell)/rho_to_nH, " T=", atmos%T(icell)
!           write(*,*) "nH(log(cm^-3)) = ", log10(1d-6*atmos%nHtot(icell))
!  		  write(*,*) "***************"
       end if !rM <= rmo and >= rmi
      end do
     end do
    end do

   !! vturb if any inclued in VBROAD_atom(icell,atom) so we do not count it here
   if (.not.lstatic) then
    atmos%v_char = atmos%v_char + maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)),&
    	   dim=1,mask=sum(atmos%Vxyz**2,dim=2)>0)
   end if
   
   if (atmos%magnetized) then
    atmos%B_char = maxval(sum(atmos%Bxyz(:,:)**2,dim=2))
    write(*,*)  "Typical Magnetic field modulus (G)", atmos%B_char * 1d4
   end if

  CALL define_atomRT_domain() !can be defined earlier so that we avoid create a dark_zone array?	
  							  !for instance, cell by cell init of icompute_atomRT, when density
  							  !is computed

!  write(*,*) "Density and T set to zero for tetsting"
!  atmos%icompute_atomRT(:) = 0

   if (include_dark_zone) then !pass a mask to define_atomRT_domain() for that case?
    where(dark_zone)
     atmos%icompute_atomRT = -1
    end where
   end if

   CALL write_atmos_domain()

   if (Tiso > 0d0) then
    where(atmos%icompute_atomRT > 0)
     atmos%T = Tiso
    end where
   end if

   CALL SolveElectronDensity()
   CALL WriteElectron(.false.)
   atmos%calc_ne = .false.
   CALL writeTemperature()
   CALL writeHydrogenDensity()
   if (.not.lstatic) CALL writeVfield ()


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
   
!    write(*,*) "Typical microturbulence in the model (km/s):"
!    write(*,*) sum(atmos%vturb)/sum(atmos%icompute_atomRT,mask=atmos%icompute_atomRT>0) /1d3

   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%icompute_atomRT>0)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%icompute_atomRT>0)  

   if (include_dark_zone) deallocate(dark_zone)
  RETURN
  END SUBROUTINE magneto_accretion_model
  

  
  SUBROUTINE spherical_shells_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with moving spherical shells
  ! ----------------------------------------------------------- !
   integer :: n_zones = 1, izone, i, j, k, icell
   double precision, parameter :: Lstar = 3.5e4, i_incl = 45d0 !deg
   double precision, parameter :: Mdot = 3.16d-7, Vinf_p = 2000d3, gamma = 0.86, Vinf_eq = 200d3
   double precision, parameter :: Vrot_eq = 326d3, Vsini = 230d3, Mdisk = 4.43d-9 !Msun
   double precision, parameter :: year_to_sec = 3.154d7, r0 = 1d0, v0_p = 0.11d3, v0_eq = 0d0
   double precision :: Rstar, Mstar, Vr, rhor, rho_to_nH, phi_p, phi_t, v0_t, Vinf_t, bigGamma, Vbreak
   double precision :: rcyl, z, r, phi, Mdot_si, rho0, O, rp, theta, sTheta, C1, m1, m2, chi, vphi
   
   linfall = .false.
   lkeplerian = .false.
   lmagnetoaccr = .true. !not really magneto accretion, but we have rotation + expansion
   						 !but, cannot be factorised in lkep + linfall with infall vel prop to keplerian
   !Remmeber, Vxyz is either (Vx, Vy, Vz) or (VR, vz, vphi); equatorial, z and azimutal 
   lstatic = .false.
   lwind_rotation = .not.lmagnetoaccr
   CALL init_atomic_atmos()
   rho_to_nH = 1d3/masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3


   Rstar = etoile(1)%r * AU_to_m !AU_to_Rsun * Rsun !m
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Mdot_si = Mdot * Msun_to_kg / year_to_sec !kg/s
   !!not used here
   !nH0 = Mdot_si / 4d0 / PI / ((r0*etoile(1)%r*AU_to_m)**2 * v0) *rho_to_nH

   m1 = 6d0
   C1 = 30d0
   m2 = 0.3

   bigGamma = 2.5e-5 * Lstar/etoile(1)%M !grad/g, M in Msun, L in Lsun
   vbreak = (Ggrav*Mstar*(1.-bigGamma)/Rstar)**0.5
   !vbreak = (Ggrav*Mstar/Rstar)**0.5 !* dsqrt(2.)
   chi = Vsini / sin(i_incl * pi / 180d0) / vbreak
   !mass flux to pole
   phi_p = 1.7d-9 * Msun_to_kg / year_to_sec / (etoile(1)%r * AU_to_m)**2 !Msun/yr/sr -> kg/s/sr/m2
   rho0 = 5.12d-11 * 1d3 !kg/m3 !!density, theta independent
   write(*,*) "Vbreak", vbreak/1d3, " chi(%)=", chi*100, " Gamma=", bigGamma
   write(*,*) "rho0 (kg/m3) = ", rho0, " phi(0) = ", phi_p
   all_loop : do i=1, n_rad
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
       theta = asin(rcyl/r)
       !theta = acos(z/r)
       sTheta = rcyl/r
       O = (1d0 - (etoile(1)%r*r0)/r)**gamma

       rp = (r/(etoile(1)%r*r0))
       
!        write(*,*) icell, "R=",rcyl, " z=", z, " r=", r, rp, O, " theta=", theta*180/pi, sTheta
       
       if (sTheta==0d0) then !pole
         Vr = v0_p + (Vinf_p - v0_p) * O
         rhor = phi_p / rp**2 / Vr
         vphi = 0d0
!           if (Vr > Vinf_p) then
!            write(*,*) Vr, v0_p, Vinf_p, O, rhor
!            stop
!           end if
       else
         phi_t = phi_p * (1+(C1-1)*sTheta**m1)
         Vinf_t = Vinf_p + (Vinf_eq - Vinf_p)*sTheta**m2
         v0_t = phi_t / rho0
         Vr = V0_t + (Vinf_t - V0_t) * O 
         vphi = sTheta * rp**-0.5 * chi * vbreak! * (Ggrav*Mstar*(1.-bigGamma)/Rstar)**0.5 
         rhor = phi_t / rp**2 / Vr

!           if (Vr > Vinf_p) then
!            write(*,*) Vr, phi_t, Vinf_t, v0_t, Vinf_eq, O
!            stop
!           end if
       end if

       atmos%T(icell) = etoile(1)%T * rp**-0.5
       atmos%nHtot(icell) = rho_to_nH * rhor
       !Project Ve_r onto e_R and e_z, to be used with v_proj()
       atmos%Vxyz(icell,1) = Vr * sTheta !R !R=rsin(t); z=rcos(t)
       atmos%Vxyz(icell,2) = Vr * dsqrt(1-sTheta**2)!z
! 	   atmos%Vxyz(icell,1) = Vr; atmos%Vxyz(icell,2) = 0d0
!        atmos%Vxyz(icell,3) = vphi!phi
!        write(*,*) "Vr=", Vr/1d3, " vphi=", vphi/1d3, rhor, " T=", atmos%T(icell), " nH=", atmos%nHtot(icell)

      end do
     end do
    end do all_loop


    atmos%v_char = atmos%v_char + maxval(dsqrt(sum(atmos%Vxyz**2,dim=2)),&
    	   dim=1,mask=sum(atmos%Vxyz**2,dim=2)>0)



   CALL define_atomRT_domain()
   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%icompute_atomRT>0)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%icompute_atomRT>0)
   write(*,*) "Typical velocity in the model (km.s^-1):"
   write(*,*) atmos%v_char/1d3
   
   CALL writeTemperature()
   CALL writeHydrogenDensity()

  RETURN
  END SUBROUTINE spherical_shells_model
  
 
  FUNCTION feqv(Rgrid_min, Rgrid_max, R1, R2, rr, inverse) result (rrp)
  !convert from mcfost grid to stellar grid
  !1.00000000000740     
  ! 49.7370802905699     

   double precision :: R1, R2, rr, rrp, Rgrid_min, Rgrid_max
   logical, optional :: inverse
   logical :: get_rs
   
   get_rs = .false.
   
   if (present(inverse)) then
    get_rs = inverse
   end if
   
   
   if (get_rs) then
     rrp = (rr - R1) * (Rgrid_max-Rgrid_min) / (R2-R1) + Rgrid_min 
    return
   end if
   
   rrp = (rr-Rgrid_min) / (Rgrid_max-Rgrid_min) * (R2-R1) + R1

  RETURN
  END FUNCTION feqv
  
  SUBROUTINE FALC_MODEL()
  !Only for testing, set atmos%XX = FALC_xx
  ! it sets an EXISTING (important !!) model to the FAL C model
   integer, parameter :: NDEP = 82
   double precision, dimension(82) :: T, ne, nHtot, nHmin
   double precision :: nH6(6,82) 
  
   data ne /1.251891d16, 1.304293d16, 1.366348d16, 1.467464d16,&
       1.603707d16, 1.694766d16, 1.811689d16, 1.969550d16,&
       2.192829d16, 2.344633d16, 2.523785d16, 2.748732d16,&
       3.043595d16, 3.394113d16, 3.804151d16, 4.291177d16,&
       4.860836d16, 5.337548d16, 5.700496d16, 6.137226d16,&
       6.520094d16, 6.747161d16, 6.918466d16, 6.991762d16,&
       6.967801d16, 6.848673d16, 6.599262d16, 6.265075d16,&
       5.801540d16, 5.841143d16, 6.141942d16, 6.651616d16,&
       7.168647d16, 7.561758d16, 7.846423d16, 8.095357d16,&
       8.445100d16, 8.838377d16, 9.185673d16, 9.612962d16,&
       1.031897d17, 1.137706d17, 1.220873d17, 1.280117d17,&
       1.340907d17, 1.308237d17, 1.229867d17, 1.111659d17,&
       9.742137d16, 8.363294d16, 7.544027d16, 9.114873d16,&
       1.307890d17, 1.818187d17, 2.441017d17, 3.294358d17,&
       4.648930d17, 7.167248d17, 1.099721d18, 1.680047d18,&
       2.557346d18, 3.880555d18, 4.803401d18, 5.984860d18,&
       7.584832d18, 9.816085d18, 1.329338d19, 1.950519d19,&
       2.749041d19, 4.028244d19, 5.456858d19, 7.641340d19,&
       1.099800d20, 1.719226d20, 2.789487d20, 4.446967d20,&
       6.869667d20, 1.041290d21, 1.531806d21, 2.194603d21,&
       2.952398d21, 3.831726d21 / 
       
  data T /100000.,  95600.,  90820.,  83890.,  75930.,  71340.,  66150.,&
        60170.,  53280.,  49390.,  45420.,  41180.,  36590.,  32150.,&
        27970.,  24060.,  20420.,  17930.,  16280.,  14520.,  13080.,&
        12190.,  11440.,  10850.,  10340.,   9983.,   9735.,   9587.,&
         9458.,   9358.,   9228.,   8988.,   8635.,   8273.,   7970.,&
         7780.,   7600.,   7410.,   7220.,   7080.,   6910.,   6740.,&
         6570.,   6370.,   6180.,   5950.,   5760.,   5570.,   5380.,&
         5160.,   4900.,   4680.,   4560.,   4520.,   4500.,   4510.,&
         4540.,   4610.,   4690.,   4780.,   4880.,   4990.,   5060.,&
         5150.,   5270.,   5410.,   5580.,   5790.,   5980.,   6180.,&
         6340.,   6520.,   6720.,   6980.,   7280.,   7590.,   7900.,&
         8220.,   8540.,   8860.,   9140.,   9400. /
  
  data nHtot /1.04571358d16, 1.09401812d16, 1.15182218d16, 1.24723003d16,&
       1.37824624d16, 1.46666937d16, 1.58073849d16, 1.73522158d16,&
       1.95244380d16, 2.09965070d16, 2.27414201d16, 2.49529201d16,&
       2.78867302d16, 3.14878205d16, 3.58761009d16, 4.12490015d16,&
       4.80215027d16, 5.41514040d16, 5.91966052d16, 6.57624071d16,&
       7.24630094d16, 7.75050116d16, 8.24780142d16, 8.72100175d16,&
       9.21940221d16, 9.65480272d16, 1.00765034d17, 1.04571040d17,&
       1.09119049d17, 1.12308055d17, 1.16667062d17, 1.28173075d17,&
       1.50485091d17, 1.83065106d17, 2.26075118d17, 2.71927128d17,&
       3.45149141d17, 4.76695155d17, 7.11571167d17, 1.04788518d18,&
       1.73738021d18, 3.04612025d18, 5.24351028d18, 9.77814031d18,&
       1.85672803d19, 3.66282003d19, 6.15684103d19, 9.80419662d19,&
       1.48043185d20, 2.29178909d20, 3.65754559d20, 6.20960707d20,&
       1.00590598d21, 1.47210490d21, 2.05290447d21, 2.84430487d21,&
       4.10450602d21, 6.41610951d21, 9.93341605d21, 1.52160292d22,&
       2.30240577d22, 3.43691248d22, 4.17101990d22, 5.01483483d22,&
       5.97136909d22, 7.04214461d22, 8.21612134d22, 9.45184976d22,&
       1.01434202d23, 1.08385618d23, 1.12638909d23, 1.16619660d23,&
       1.20282003d23, 1.22812486d23, 1.24607714d23, 1.25981369d23,&
       1.27401250d23, 1.28402892d23, 1.29180598d23, 1.30020496d23,&
       1.31275567d23, 1.32662462d23/
  

   if (n_cells >= NDEP) then
    atmos%ne(1:NDEP) = ne
    atmos%T(1:NDEP) = T
    atmos%nHtot(1:NDEP) = nHtot
   else
    CALL Error("Cannot associate FAL_C to MCFOST, increase number of depth points")
   end if
  
  RETURN
  END SUBROUTINE FALC_MODEL
  
  
  !building
  SUBROUTINE spherical_star()
  ! ----------------------------------------------------------- !
   ! star with MCFOST
  ! ----------------------------------------------------------- !
   integer :: n_zones = 1, izone, i, j, k, icell, m, NDEP_lim
   integer, parameter :: Pmax = 4, NDEP = 82
   double precision, parameter :: r0 = 1d0, v0 = 1d0, Prot = 8.
   double precision :: Vr, rhor, rho_to_nH, R0s, Mstar, Rstar
   double precision :: rcyl, z, r, phi, rs, omega, Vrot, cm, Rsphere_mcfost(n_cells), Rsmin, Rsmax
   double precision, dimension(NDEP) :: m_mod, r_mod, rho_mod, T_mod, ne_mod, V_mod, f_mod, mass_mod, cm2, S_mod
   double precision, dimension(Pmax) :: p0, p1, p2, p3
   
   Rsphere_mcfost(:) = dsqrt(r_grid**2 + z_grid**2)
   Rsmin = minval(Rsphere_mcfost); Rsmax = maxval(Rsphere_mcfost)
   NDEP_lim = NDEP
   
   Data m_mod / -4.58623126, -4.5857112 , -4.58516185, -4.58438554, -4.58353816,&
       -4.58306938, -4.58256495, -4.58200811, -4.58140527, -4.5810861 ,&
       -4.58076716, -4.58042942, -4.58006653, -4.57968072, -4.57927198,&
       -4.5788341 , -4.57827672, -4.5777787 , -4.57733569, -4.57670219,&
       -4.57592287, -4.57519269, -4.5743236 , -4.57320574, -4.57154354,&
       -4.56932934, -4.56601185, -4.56169625, -4.55427337, -4.52520214,&
       -4.47294526, -4.36681641, -4.23030561, -4.09593027, -3.96939025,&
       -3.86606284, -3.74136425, -3.58681894, -3.41002803, -3.2474327 ,&
       -3.04481118, -2.82707945, -2.62103882, -2.38742203, -2.1459382 ,&
       -1.88455528, -1.67860661, -1.49403941, -1.33345229, -1.16537024,&
       -0.98799858, -0.78195113, -0.58617383, -0.42570423, -0.2838909 ,&
       -0.1417349 ,  0.01997558,  0.21965262,  0.41650311,  0.60999289,&
        0.79959933,  0.98478766,  1.07597976,  1.16545643,  1.25280141,&
        1.33769133,  1.41984591,  1.49890025,  1.54458787,  1.5887034 ,&
        1.61721505,  1.64493425,  1.67183161,  1.69780313,  1.72273204,&
        1.74660052,  1.76947003,  1.7913973 ,  1.81241214,  1.83257963,&
        1.85200451,  1.87078221/
    
   !
   data r_mod / 1.00355514, 1.00355306, 1.00355097, 1.00354819, 1.00354541,&
       1.00354399, 1.00354256, 1.00354111, 1.0035397 , 1.00353901,&
       1.00353838, 1.00353777, 1.00353717, 1.00353661, 1.00353608,&
       1.00353559, 1.00353504, 1.00353462, 1.00353428, 1.00353383,&
       1.00353334, 1.00353291, 1.00353243, 1.00353185, 1.00353103,&
       1.00352999, 1.00352849, 1.0035266 , 1.00352342, 1.00351091,&
       1.00348702, 1.00343246, 1.0033508 , 1.00325907, 1.00316402,&
       1.00308104, 1.00297584, 1.00284057, 1.0026836 , 1.00253963,&
       1.00236651, 1.00219052, 1.00203416, 1.00187142, 1.00171724,&
       1.001564  , 1.00144842, 1.00134871, 1.00126514, 1.00118178,&
       1.00109871, 1.00100801, 1.00092508, 1.00085801, 1.00079907,&
       1.00074009, 1.00067293, 1.00058963, 1.00050623, 1.00042273,&
       1.00033912, 1.00025537, 1.00021285, 1.00017032, 1.00012776,&
       1.00008519, 1.0000426 , 1.        , 0.9999744 , 0.99994879,&
       0.99993172, 0.99991465, 0.99989757, 0.9998805 , 0.99986343,&
       0.99984635, 0.99982928, 0.99981221, 0.99979514, 0.99977807,&
       0.999761  , 0.99974392 /
   
   data T_mod / 100000.,  95600.,  90820.,  83890.,  75930.,  71340.,  66150.,&
        60170.,  53280.,  49390.,  45420.,  41180.,  36590.,  32150.,&
        27970.,  24060.,  20420.,  17930.,  16280.,  14520.,  13080.,&
        12190.,  11440.,  10850.,  10340.,   9983.,   9735.,   9587.,&
         9458.,   9358.,   9228.,   8988.,   8635.,   8273.,   7970.,&
         7780.,   7600.,   7410.,   7220.,   7080.,   6910.,   6740.,&
         6570.,   6370.,   6180.,   5950.,   5760.,   5570.,   5380.,&
         5160.,   4900.,   4680.,   4560.,   4520.,   4500.,   4510.,&
         4540.,   4610.,   4690.,   4780.,   4880.,   4990.,   5060.,&
         5150.,   5270.,   5410.,   5580.,   5790.,   5980.,   6180.,&
         6340.,   6520.,   6720.,   6980.,   7280.,   7590.,   7900.,&
         8220.,   8540.,   8860.,   9140.,   9400. /
   
   data rho_mod / 2.08373931d-11, 2.17999328d-11, 2.29517644d-11, 2.48529073d-11,&
       2.74635995d-11, 2.92255613d-11, 3.14985577d-11, 3.45768623d-11,&
       3.89053370d-11, 4.18386528d-11, 4.53156507d-11, 4.97223924d-11,&
       5.55684439d-11, 6.27441500d-11, 7.14884493d-11, 8.21947503d-11,&
       9.56899629d-11, 1.07904700d-10, 1.17958011d-10, 1.31041345d-10,&
       1.44393288d-10, 1.54440225d-10, 1.64349670d-10, 1.73778888d-10,&
       1.83710256d-10, 1.92386256d-10, 2.00789267d-10, 2.08373299d-10,&
       2.17435881d-10, 2.23790448d-10, 2.32476416d-10, 2.55403853d-10,&
       2.99863851d-10, 3.64784360d-10, 4.50488186d-10, 5.41855113d-10,&
       6.87760828d-10, 9.49885761d-10, 1.41791103d-09, 2.08806656d-09,&
       3.46198761d-09, 6.06984615d-09, 1.04484715d-08, 1.94843941d-08,&
       3.69980585d-08, 7.29871189d-08, 1.22684184d-07, 1.95363151d-07,&
       2.94997992d-07, 4.56672950d-07, 7.28820178d-07, 1.23735626d-06,&
       2.00441679d-06, 2.93338724d-06, 4.09071647d-06, 5.66769910d-06,&
       8.17883670d-06, 1.27850493d-05, 1.97938040d-05, 3.03201937d-05,&
       4.58788481d-05, 6.84855761d-05, 8.31137546d-05, 9.99280179d-05,&
       1.18988381d-04, 1.40325171d-04, 1.63718397d-04, 1.88342118d-04,&
       2.02122684d-04, 2.15974412d-04, 2.24449724d-04, 2.32381960d-04,&
       2.39679721d-04, 2.44722083d-04, 2.48299342d-04, 2.51036554d-04,&
       2.53865877d-04, 2.55861796d-04, 2.57411491d-04, 2.59085112d-04,&
       2.61586027d-04, 2.64349622d-04 /
   
   data ne_mod / 1.251891d16, 1.304293d16, 1.366348d16, 1.467464d16,&
       1.603707d16, 1.694766d16, 1.811689d16, 1.969550d16,&
       2.192829d16, 2.344633d16, 2.523785d16, 2.748732d16,&
       3.043595d16, 3.394113d16, 3.804151d16, 4.291177d16,&
       4.860836d16, 5.337548d16, 5.700496d16, 6.137226d16,&
       6.520094d16, 6.747161d16, 6.918466d16, 6.991762d16,&
       6.967801d16, 6.848673d16, 6.599262d16, 6.265075d16,&
       5.801540d16, 5.841143d16, 6.141942d16, 6.651616d16,&
       7.168647d16, 7.561758d16, 7.846423d16, 8.095357d16,&
       8.445100d16, 8.838377d16, 9.185673d16, 9.612962d16,&
       1.031897d17, 1.137706d17, 1.220873d17, 1.280117d17,&
       1.340907d17, 1.308237d17, 1.229867d17, 1.111659d17,&
       9.742137d16, 8.363294d16, 7.544027d16, 9.114873d16,&
       1.307890d17, 1.818187d17, 2.441017d17, 3.294358d17,&
       4.648930d17, 7.167248d17, 1.099721d18, 1.680047d18,&
       2.557346d18, 3.880555d18, 4.803401d18, 5.984860d18,&
       7.584832d18, 9.816085d18, 1.329338d19, 1.950519d19,&
       2.749041d19, 4.028244d19, 5.456858d19, 7.641340d19,&
       1.099800d20, 1.719226d20, 2.789487d20, 4.446967d20,&
       6.869667d20, 1.041290d21, 1.531806d21, 2.194603d21,&
       2.952398d21, 3.831726d21 /     

  !put in decreasing order, because first point is the inner point
   m_mod(NDEP:1:-1) = m_mod(1:NDEP:1)
   r_mod(NDEP:1:-1) = r_mod(1:NDEP:1)
   T_mod(NDEP:1:-1) = T_mod(1:NDEP:1)
   ne_mod(NDEP:1:-1) = ne_mod(1:NDEP:1)
   rho_mod(NDEP:1:-1) = rho_mod(1:NDEP:1)

   linfall = .false. 
   lkeplerian = .false.
   lmagnetoaccr = .false. 
   lstatic = .true.
   CALL init_atomic_atmos()
   rho_to_nH = 1d3/masseH /atmos%avgWeight !density kg/m3 -> nHtot m^-3
   
   f_mod(:) = ne_mod(:)/(rho_mod(:)*rho_to_nH)
   !V_mod(1) = (abs(r_mod(2)-r_mod(1))/2d0)**3 * Rsun**3 !m**3
   !V_mod(1) = 4./3. * pi * (abs(r_mod(2)-r_mod(1)))**3 * Rsun**3 !m**3
   V_mod(1) = 4./3. * pi * r_mod(1)**3 * Rsun**3
   !V_mod(1) = (4./3 * pi * Rsun**3) * abs(r_mod(2)-r_mod(1))**3 * 2d0**-3
   
   S_mod(1) = 4.*pi * r_mod(1)**2 * Rsun**2
   k = 1
   mass_mod(1) = rho_mod(1) * V_mod(1)
   write(*,*) k, r_mod(k), mass_mod(k), V_mod(k)
   do k=2, NDEP
    !V_mod(k) = 4./3. * pi * (abs(r_mod(k-1)-r_mod(k))/2d0)**3 * Rsun**3
    !V_mod(k) = 4./3. * pi * (abs(r_mod(k)-r_mod(k-1)))**3 * Rsun**3 !m**3
    !V_mod(k) = (4./3 * pi * Rsun**3) * abs(r_mod(k)-r_mod(k-1))**3 * 2d0**-3
    V_mod(k) = r_mod(k)**3 * 4./3. * pi * Rsun**3
    !V_mod(k) = (4./3 * pi * Rsun**3) * abs(r_mod(k)-r_mod(k-1))**3
    S_mod(k) = 4. * pi * r_mod(k)**2 * Rsun**2
    mass_mod(k) = rho_mod(k) * V_mod(k)
    write(*,*)k, r_mod(k), mass_mod(k), V_mod(k)
   end do
   
   cm2 = m_mod * S_mod !kg


   R0s = etoile(1)%r!inner boundary in au
   Mstar = etoile(1)%M * Msun_to_kg !kg
   Rstar = feqv(Rsmin, Rsmax, minval(r_mod), maxval(r_mod), 1d0, .true.) !in au
   write(*,*) "Rstar in mcfost grid (au) = ", Rstar, (Rsmax - Rstar)/Rstar, 1., (maxval(r_mod)-1.)
   
   if (lstatic) then
    if (allocated(atmos%Vxyz)) deallocate(atmos%Vxyz)
   end if
   
   all_loop : do i=1, n_rad
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

       rs = feqv(Rsmin, Rsmax, minval(r_mod), maxval(r_mod), r)!in Rstar
       !cm = interp_dp(m_mod, r_mod, rs)
       cm = interp_dp(cm2, r_mod, rs)
       !then interpolate
  
       !! UNIFORM LAW
!    CALL Warning(&
!    	"Beware, I do not take into account vturb anymore. Small discrepancies can arise")
!       if (r <= 0.8*Rsmax) then !remove if for uniform squared
!        atmos%nHtot(icell) = 1d15
!        atmos%T(icell) = 3000d0
!        atmos%ne(icell) = 1d12
!       end if
       !!
! 
       !used a different scale for Rstar to au
       atmos%nHtot(icell) = interp_dp(mass_mod, cm2, cm) * rho_to_nH / (Volume(icell) * Au3_to_m3)
       atmos%T(icell) = interp_dp(T_mod, cm2, cm)
       !atmos%ne(icell) = interp_dp(ne_mod, m_mod, cm)
       atmos%ne(icell) = interp_dp(f_mod, cm2, cm) * atmos%nHtot(icell)
       write(*,*) r, rs, cm, interp_dp(mass_mod, m_mod, cm)

      end do
     end do
    end do all_loop
    

   CALL define_atomRT_domain()

   CALL writeTemperature()
   CALL writeHydrogenDensity()
   CALL write_atmos_domain()

   atmos%calc_ne = .false.
   
   if (lkeplerian) atmos%v_char = maxval(Vfield)

   write(*,*) "Maximum/minimum Temperature in the model (K):"
   write(*,*) MAXVAL(atmos%T), MINVAL(atmos%T,mask=atmos%icompute_atomRT>0)
   write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHtot), MINVAL(atmos%nHtot,mask=atmos%icompute_atomRT>0) 

  RETURN
  END SUBROUTINE spherical_star

 END MODULE simple_models