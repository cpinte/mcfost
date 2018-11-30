MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos
 use input
 use parametres
 use grid
 use density
 use constantes
 
 
 IMPLICIT NONE

 CONTAINS
 
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid
   ! Values are from a Solar model atmosphere by Fontenla et al.
   ! Namely the FAL C model.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne

!    !idk = 10
   nHtot =  2.27414200581936d16
   T = 45420d0
   ne = 2.523785d16

!    !idk = 75
!    T=7590d0
!    ne = 4.446967d20
!    nHtot = 1.259814d23

!    !idk = 81   
!     T=9400d0
!     ne = 3.831726d21
!     nHtot = 1.326625d23

!    !idk = 0 
!     T = 100000d0
!     ne = 1.251891d16
!     nHtot = 1.045714d16
    
   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving = .true.
   atmos%vturb = 9.506225d3 !m/s !idk=10
   !atmos%vturb = 1.696164d3 !idk = 75
   !atmos%vturb = 1.806787d3 !idk=81
   !atmos%vturb = 10.680960d3 !idk=0
   if (atmos%moving .and. .not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
   atmos%Vmap = 0d0
   atmos%velocity_law = 0 !keplerian

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE prop_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to dust density and temperature.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne

   !nHtot = 1d23 * densite_gaz/MAXVAL(densite_gaz)
   nHtot =  1d9 * densite_gaz * masse_mol_gaz / m3_to_cm3 / masseH
   T = Tdust * 100d0!100d0, depends on the stellar flux
   ne = 1d-2 * nHtot
   
!   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving=.true.
   if (atmos%moving .and. .not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
   atmos%Vmap = 0d0
   atmos%velocity_law = 0 !keplerian
   
  RETURN
  END SUBROUTINE prop_law_model
  
  SUBROUTINE radial_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to r^-n
  ! ----------------------------------------------------------- !
  integer :: izone, i, j, k, icell, n_zones=1 !Only one component firstly
  double precision, dimension(n_cells) :: Vr, Tr,vturbr, nHtotr, ner
  type(disk_zone_type) :: dz ! to define the properties of the model
  double precision :: r, nH0, ne0, T0, vt, rcyl, z, v0
  
  nH0 = 1d16
  ne0 = 1d-2 * nH0
  T0 = 1d4
  vt = 0d0
  v0 = 40d3
  
  Tr = 1d0
  nHtotr = 1d0
  Vr = 0d0
  vturbr = 0d0
  ner = 1d0

  do izone=1, n_zones
   !dz = disk_zone(izon) ! for now the parameters are hardcoded
    do i=1, n_rad
     !do j=j_start,nz
     j = 0
      do k=1, n_az
       if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
        Tr(icell) = T0 * 5d0/rcyl
        nHtotr(icell) = nH0 * 5d0/rcyl
        ner(icell) = 1d-2 * nHtotr(icell)
        Vr(icell) = v0 * dsqrt(1.-5d0/rcyl)**0.5 !m/s
        vturbr(icell) = vt
!        else
!         icell = cell_map(i,j,k)
!         rcyl = r_grid(icell)
!         z = z_grid(icell)/z_scaling_env
       endif
!        r = sqrt(rcyl**2 + z**2)
!        if ((r <= 300).and.(r>=5)) then
!         Tr(icell) = T0 * 5d0/r
!         nHtotr(icell) = nH0 * 5d0/r
!         ner(icell) = 1d-2 * nHtotr(icell)
!         Vr(icell) = v0 * 5d0/r !m/s
!         vturbr(icell) = vt
!        end if
      end do !k
     !end do !j
    end do !i
  end do !over n_zones
 

  CALL init_atomic_atmos(n_cells, Tr, ner, nHtotr)
  atmos%moving=.true.
  if (atmos%moving .and. .not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
  atmos%Vmap = 0d0
  atmos%Vmap = Vr !and set keyword lkeplerian=.false. and linfall=.true.
  				  ! otherwise vinfall/vkep = cte = expected.
  atmos%velocity_law = -1
  atmos%v_char = 0.5*MAXVAL(atmos%Vmap)
  RETURN
  END SUBROUTINE radial_model

  SUBROUTINE define_circumstellar_envelope()
  ! ----------------------------------------------------------------------------------------- !
   ! Generate a model to test line RT modules.
   ! This test case corresponds to the following MARCS (cgs units) model:
   ! s3500_g+0.5_m1.0_t05_st_z+0.00_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod
   ! Which corresponds to a BB for Teff=3500K,
   !   R = 92.79 Rsun, logg = 0.5 (g in g/cm2). Note however that the model expand up to
   !   R * (1+dr) with dr = (r0-R)/R about 0.07
   !
   ! Computes a spherically symmetric circumstellar envelope
   ! from mass conservation law and radiative equilibrium for the temperature (T \propto r^-2)
   ! and constant mass-loss.
   !
   ! V(r) = v0 + (vinf-v0) * (1.- r0/r) ** beta
   ! T(r) = T0 * sqrt(r0/r)
   ! v0 = Mdot / (4*PI * r0**2 * rho0)
   ! rho(r) = Mdot / (4*PI * r**2 * V(r))
   !
   ! Bonus:  Add the dust contribution -> 
   !           modelling circumstellar environnement of evolved stars
   !
  ! ----------------------------------------------------------------------------------------- ! 
  integer :: izone, i, j, k, icell, n_zones=1 !Only one component firstly
  double precision, dimension(n_cells) :: Venv, Tenv, rhoenv, neenv, vturb, nHtot
  type(disk_zone_type) :: dz ! to define the properties of the model
  
  ! Radius in Rstar, units SI
  double precision :: rhorho, r, z, vel_char, rcyl, RstarSun = 92.79, Rstar_to_AU
  double precision :: r0, v0, T0, rho0, rmax, Rstar, Mdot, yr_2_sec, Msun, PI, beta, vinf
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Fixed initial parameters
  neenv = 0d0
  rhoenv = 0d0
  Venv = 0d0
  Tenv = 0d0
  nHtot = 0d0
  beta = 5d-1 ! CAK original-law
  vinf = 40d0 * 1000 !m/s
  Msun = 2.0d30 !kg
  yr_2_sec = 3.154d7 !s/yr
  Rstar=6.4953d10!m
  PI = 314.0d-1
  r0 = -1 * (-4.580d9) + Rstar! m ! r = -depth + Rstar
  T0 = 2362.0
  rho0 = 7.107d-6!kg/m3
  Mdot = 1d-4 * Msun / yr_2_sec ! Msun/yr -> kg/s
  v0 = Mdot / (4*PI * r0**2 * rho0) !m/s
  rmax = 10d0 !Rstar
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  r0 = r0 / Rstar!Rstar
  Rstar_to_AU = Rstar / AU_to_m

  do izone=1, n_zones
   !dz = disk_zone(izon) ! for now the parameters are hardcoded
    do i=1, n_rad
     do j=j_start,nz
      do k=1, n_az
       if (j==0) then !pourquoi ici
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) * AU_TO_RSUN / RstarSun !Rstar
        z = 0.0_dp
       else
        icell = cell_map(i,j,k)
        rcyl = r_grid(icell) * AU_TO_RSUN / RstarSun
        z = z_grid(icell) * AU_TO_RSUN / RstarSun
       endif
        r = sqrt(rcyl**2 + z**2) !Rstar
        if ((r <= rmax).and.(r>=r0).and.(j/=0)) then
         !if (j/=0) then 
           Venv(icell) = v0 + (vinf-v0) * (1d0 - r0/r) ** beta !m/s
           Tenv(icell) = T0 * dsqrt(r0/r)
           rhoenv(icell) = Mdot / (4*pi * (Rstar*r)**2 * Venv(icell))
           nHtot(icell) = 1000 * rhoenv(icell) / (masseH) !m^-3
!             write(*,*) "icell=", icell, "rcyl=",rcyl*Rstar_to_AU, &
!                                      "z=",z*Rstar_to_AU,       &
!                                      ' r=',r*Rstar_to_AU,      &
!                                      " rmax=",rmax*Rstar_to_AU,&
!                                      " r0=", r0*Rstar_to_AU,   & 
!                                      rhoenv(icell), nHtot(icell), Tenv(icell), Venv(icell)
         !end if
        !else
           !Venv(icell) = 1d-30
           !Tenv(icell) = 1d-30
           !rhoenv(icell) = 1d-30
           !nHtot(icell) = 1d-30
        end if
!         write(*,*) "icell=", icell, "rcyl=",rcyl*Rstar_to_AU, &
!                                     "z=",z*Rstar_to_AU,       &
!                                     ' r=',r*Rstar_to_AU,      &
!                                     " rmax=",rmax*Rstar_to_AU,&
!                                     " r0=", r0*Rstar_to_AU,   & 
!                                     rhoenv(icell), nHtot(icell), Tenv(icell), Venv(icell)
      end do !k
     end do !j
    end do !i
  end do
 
 ! Initialise the atmosphere for line RT
 neenv = 1d-4 * nHtot ! can be recomputed
 CALL init_atomic_atmos(n_cells, Tenv, neenv, nHtot)
 ! if commented, velocity from file or from keplerian is used
  if (atmos%moving .and. .not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
  atmos%Vmap = 0d0
  atmos%Vmap = Venv
  atmos%velocity_law = -1
 
 RETURN
 END SUBROUTINE define_circumstellar_envelope
 

 
!  SUBROUTINE ContributionFunction()
!   integer :: icell, iray, id, NrecStokes, nbr_cell, previous_cell, next_cell,icellinf
!   double precision :: x0, y0, z0, norme, u0, v0, w0, x1, y1, z1
!   logical :: labs, lcellule_non_vide, lsubtract_avg, lonly_top, lonly_bottom
!   double precision :: ksi(NLTEspec%Nwaves,atmos%Nspace,atmos%Nrays,1), l
!                                                            !Replace 1 by 4 if pol.
!   double precision :: l_contrib, l_void_before, Ic(1,NLTEspec%Nwaves,atmos%Nrays)
!   double precision, dimension(NLTEspec%Nwaves) :: dtau, dtau_c, tau, tau_c, facteur_tau, Snu, Snu_c
!   ! for writing                                                        
!   integer :: status,unit,blocksize,bitpix,naxis
!   integer, dimension(7) :: naxes
!   integer :: group,fpixel,nelements, i, xcenter
!   logical :: simple, extend
!   character(len=512) :: CNTRB_FILE
!   CNTRB_FILE = "CNTRB.fits.gz"
! 
!   ! ------------------------------------------ !  
!   write(*,*) "Computing contribution function(s)..."
!   
!   NrecStokes = 1 !up to now
!   Ic = 0d0
!   dtau = 0d0
!   dtau_c = 0d0
!   tau = 0d0
!   tau_c = 0d0
!   nbr_cell = 0
!   Snu = 0d0
!   Snu_c = 0d0
!   id = 1 !sequentiel
!   do icell=1, n_cells
!    ! Propagation des rayons
!    do iray=1, atmos%Nrays
!     ! Position = milieu de la cellule
!     x0 = r_grid(icell)
!     y0 = 0.0_dp
!     z0 = z_grid(icell)
! 
!     norme = sqrt(x0*x0 + y0*y0 + z0*z0)
!     if (iray==1) then
!      u0 = x0/norme
!      v0 = y0/norme
!      w0 = z0/norme
!     else
!      u0 = -x0/norme
!      v0 = -y0/norme
!      w0 = -z0/norme
!     end if
!     ! Integration le long du rayon
!     ! Boucle infinie sur les cellules
!     next_cell = icell
!     infinie : do ! Boucle infinie
!      ! Indice de la cellule
!      icellinf = next_cell
!      x0=x1 ; y0=y1 ; z0=z1
!      if (icellinf <= n_cells) then
!         lcellule_non_vide=.true.
!      else
!         lcellule_non_vide=.false.
!      endif
!      
!      ! Test sortie
!      if (test_exit_grid(icellinf, x0, y0, z0)) then
!         exit infinie
!      endif
! 
!      nbr_cell = nbr_cell + 1
! 
!      ! Calcul longeur de vol et profondeur optique dans la cellule
!      previous_cell = 0 ! unused, just for Voronoi
!      call cross_cell(x0,y0,z0, u0,v0,w0,  icellinf, previous_cell, x1,y1,z1, next_cell, &
!                      l, l_contrib, l_void_before)
! 
!      if (lcellule_non_vide) then
!      lsubtract_avg = ((nbr_cell == 1).and.labs)
! 
!       ! opacities in m^-1
!       l_contrib = l_contrib * AU_to_m
! 
! 
!       CALL Background(1,icellinf, x0, y0, z0, u0, v0, w0)
!       dtau(:) =  l_contrib * (NLTEspec%ActiveSet%chi_c(:)+NLTEspec%ActiveSet%chi(:)) !scattering + thermal
!       dtau_c(:) = l_contrib * NLTEspec%ActiveSet%chi_c_bf(:)
! 
!       Snu = (NLTEspec%ActiveSet%eta_c + & 
!                   NLTEspec%ActiveSet%eta) / &
!                  (NLTEspec%ActiveSet%chi_c + NLTEspec%ActiveSet%chi)
!       ! continuum source function
!       Snu_c = (NLTEspec%ActiveSet%eta_c_bf) / NLTEspec%ActiveSet%chi_c_bf
! 
!       Ic(id,:,iray) = Ic(id,:,iray)*dexp(-dtau_c) + Snu_c * dexp(-dtau_c) * dtau_c
! 
!       facteur_tau = 1d0
!       if (lonly_top    .and. z0 < 0.) facteur_tau = 0d0
!       if (lonly_bottom .and. z0 > 0.) facteur_tau = 0d0
! 
!       tau = tau + dtau * facteur_tau
!       tau_c = tau_c + dtau_c
! 
!         ! set opacities to 0.0 for next cell point.
!         CALL initAS(re_init=.true.)
!      end if  ! lcellule_non_vide
!     end do infinie
!   !---------------------------------------------!    
!     ksi(:,icell,iray,1) = (Ic(id, :, iray) - Snu)*dexp(-tau)
!     
!    end do ! iray
!   end do ! icell
!   
!   ! ------------------------------------------ !
!   ! Now write them to disk
!   status=0
!   CALL ftgiou (unit,status)
! 
!   !  Create the new empty FITS file.
!   blocksize=1
!   CALL ftinit(unit,trim(CNTRB_FILE),blocksize,status)
! 
!   simple=.true.
!   extend=.false.
!   if (NrecStokes > 1) extend = .true.
!   group=1
!   fpixel=1  
!   
!   bitpix=-64
!   naxes(1)=NLTEspec%Nwaves
!   if (lVoronoi) then
!    if (NrecStokes > 1) then 
!     naxis = 4
!     naxes(2)=n_cells
!     naxes(3)=atmos%Nrays
!     naxes(4)=NrecStokes !if only Q NrecStokes = 2 for IQ etc
!     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
!    else
!     naxis = 3 !only ksi_I
!     naxes(2)=n_cells
!     naxes(3)=atmos%Nrays
!     nelements=naxes(1)*naxes(2)*naxes(3)
!    end if
!   else
!    if (l3D) then
!     if (NrecStokes > 1) then 
!      naxis = 6
!      naxes(2) = n_rad
!      naxes(3) = 2*nz
!      naxes(4) = n_az
!      naxes(5)=atmos%Nrays
!      naxes(6)=NrecStokes !if only Q NrecStokes = 2 for IQ etc
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
!     else
!      naxis = 5
!      naxes(2) = n_rad
!      naxes(3) = 2*nz
!      naxes(4) = n_az
!      naxes(5)=atmos%Nrays
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
!     end if
!   else !not l3D
!     if (NrecStokes > 1) then 
!      naxis = 5
!      naxes(2) = n_rad
!      naxes(3) = nz
!      naxes(4)=atmos%Nrays
!      naxes(5)=NrecStokes !if only Q NrecStokes = 2 for IQ etc
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
!     else
!      naxis = 4
!      naxes(2) = n_rad
!      naxes(3) = nz
!      naxes(4)=atmos%Nrays
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
!     end if
!   end if !not l3D
!  end if !end if not Voronoi
! 
! 
!   !  Write the required header keywords.
!   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!   CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.sr-1",'ksi',status)
! 
!   CALL ftpprd(unit,group,fpixel,nelements,ksi(:,:,:,1),status)
! !   if (NrecStokes > 1) then
! !    !write other ksi
! !   end if
! 
!   !  Close the file and free the unit number.
!   CALL ftclos(unit, status)
!   CALL ftfiou(unit, status)
! 
!   !  Check for any error, and if so print out error messages
!   if (status > 0) then
!      CALL print_error(status)
!   endif 
!  RETURN
!  END SUBROUTINE ContributionFunction 
 
 END MODULE simple_models