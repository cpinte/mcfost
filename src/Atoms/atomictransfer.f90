! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atoms, using the MALI scheme.
!
! Outputs:
! - Flux (\lambda) [J.s^{-1}.m^{-2}.Hz^{-1}]
! - Irradiation map around a line [J.s^{-1}.m^{-2}.Hz^{-1}.pix^{-1}] !sr^{-1}
!
!
! Note: SI units, velocity in m/s, density in kg/m3, radius in m
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

MODULE AtomicTransfer

 use metal, only                        : Background, BackgroundContinua, BackgroundLines, &
 										  storeBackground
 use opacity
 use Profiles
 use spectrum_type
 use atmos_type
 use readatom
 use lte
 use constant, only 					: MICRON_TO_NM
 use collision
 use solvene
 use statequil_atoms
 use writeatom
 use simple_models, only 				: magneto_accretion_model, magneto_accretion_diskwind_model, &
 										 FALC_MODEL, spherical_star, feqv
 !$ use omp_lib

 !MCFOST's original modules
 use input
 use parametres
 use grid
 use optical_depth, only				: atom_optical_length_tot, optical_length_tot
 use dust_prop
 use dust_transfer, only 				: compute_stars_map
 use dust_ray_tracing, only 			: init_directions_ray_tracing ,            &
                              			  tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az, &
                              			  tab_RT_incl, stars_map, kappa, stars_map_cont
 use stars
 use wavelengths
 use density

 IMPLICIT NONE
 
 logical :: lmcfost_star !TMP
 
 PROCEDURE(INTEG_RAY_LINE_I), pointer :: INTEG_RAY_LINE => NULL()
 
 !Try to give a good idea to the user on how much RAM memory will be used
 integer :: Total_size_1 = 0, Total_size_2 = 0

 CONTAINS
 
 SUBROUTINE INTEG_RAY_LINE_I(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray
  ! crossing different cells.
  ! --> Atomic Lines case.
  !
  ! voir integ_ray_mol from mcfost
  ! All necessary quantities are initialised before this routine, and everything
  ! call, update, rewrite atoms% and spectrum%
  ! if atmos%nHtot or atmos%T is 0, the cell is empty
  !
  ! id = processor information, iray = index of the running ray
  ! what is labs? lsubtract_avg?
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  double precision, intent(in) :: u,v,w
  double precision, intent(in) :: x,y,z
  logical, intent(in) :: labs !used in NLTE but why?
  double precision :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  double precision, dimension(NLTEspec%Nwaves) :: Snu, Snu_c, etau, Istar
  double precision, dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  double precision :: facteur_tau !used only in molecular line to have emission for one
                                  !  half of the disk only. Note used in AL-RT.
                                  !  this is passed through the lonly_top or lonly_bottom.
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp !go from surface down to the star
  tau_c = 0.0_dp

  ! Reset, because when we compute the flux map
  ! with emission_line_map, we only use one ray, iray=1
  ! and the same space is used for emergent I.
  ! Therefore it is needed to reset it.
  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0
  Istar(:) = 0d0
  ! -------------------------------------------------------------- !
  !*** propagation dans la grille ***!
  ! -------------------------------------------------------------- !
  ! Will the ray intersect a star
  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)
  !Computes the stellar radiation for this direction
  !cannot be here if limbD depends on x0,y0,z0 than x,y,z
!   if (lintersect_stars) CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x, y, z, u,v,w,Istar)

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
    ! Indice de la cellule
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     !lcellule_non_vide=.true.
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif

    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN
    
    nbr_cell = nbr_cell + 1
!     if (abs(abs(z0)-0.05*etoile(1)%R)<=0.05*etoile(1)%R) RETURN
    if (lintersect_stars) then
      if (icell == icell_star) then !this is equivalent to compute_stars_map()
       !if we start at icell, computes the radiation from the star to icell.
       !in particular, the cell icell can be the cell for which we solve the SEE/or
       !a cell at the border of the grid for an image
       CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,Istar)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + Istar*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + Istar*dexp(-tau_c)
       RETURN
      end if
    endif

    !nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
    previous_cell = 0 ! unused, just for Voronoi

    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, &
                     l, l_contrib, l_void_before)

    ! ************* TEST ************* !
     if (lmcfost_star) then
       if (nbr_cell == 1) CALL Warning("Rescaling path along cell!!")
       l = feqv(0.99974392*1d0,1.00355514*1d0,l); l_contrib = feqv(0.99974392*1d0,1.00355514*1d0,l_contrib)
     end if
    ! ******************************** !

!     if (.not.atmos%lcompute_atomRT(icell)) lcellule_non_vide = .false. !chi and chi_c = 0d0, cell is transparent  
!	   this makes the code to break down  --> the atomRT is defined on n_cells
!		but the infinite loop runs over "virtual" cells, which can be be 
!		out of atomRT boundary

    !count opacity only if the cell is filled, else go to next cell
    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m !l_contrib in AU
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m !we enter at the cell, icell_in
     														  !and we need to compute the length path
     
     !! Compute background opacities for PASSIVE bound-bound and bound-free transitions
     !! at all wavelength points including vector fields in the bound-bound transitions
     eval_operator = (labs .and. (nbr_cell == 1)) !labs if false for images
     											  !so no pb if Nact>0 and we use a different grid
     CALL initAtomOpac(id,eval_operator) !set opac to 0 for this cell and thread id
     CALL NLTEopacity(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
     if (eval_operator) then
	   !in practice only if MALI
	   CALL initCrossCoupling(id)
       CALL FillCrossCoupling_terms(id, icell)
       CALL add_to_psi_operator(id, icell, iray, ds(iray,id))
     end if
     !never enter NLTEopacity if no activeatoms
     if (lstore_opac) then !not updated during NLTE loop, just recomputed using initial pops
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      dtau(:)   = l_contrib * (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1))
      dtau_c(:) = l_contrib * NLTEspec%AtomOpac%Kc(icell,:,1)
      Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
               NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / &
                 (NLTEspec%AtomOpac%Kc(icell,:,1) + &
                 NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%chi(:,id) + 1d-300)
                 ! + &NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%J(:,id)
                 ! + &NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%Jc(:,id)
      Snu_c = (NLTEspec%AtomOpac%jc(icell,:)) / &
            (NLTEspec%AtomOpac%Kc(icell,:,1) + 1d-300) !it is missing NLTE cont
     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l) !x,y,z,u,v,w,x1,y1,z1
                                !define the projection of the vector field (velocity, B...)
                                !at each spatial location.

      ! Epaisseur optique
      ! chi_p contains both thermal and continuum scattering extinction
      dtau(:)   =  l_contrib * (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id))
      dtau_c(:) = l_contrib * NLTEspec%AtomOpac%chi_c(:,id)

      ! Source function
      ! No dust yet
      ! J and Jc are the mean radiation field for total and continuum intensities
      ! it multiplies the continuum scattering coefficient for isotropic (unpolarised)
      ! scattering. chi, eta are opacity and emissivity for ACTIVE lines.
      Snu = (NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / &
                 (NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%chi(:,id) + 1d-300)
      ! + &NLTEspec%AtomOpac%sca_c(:,id) * NLTEspec%J(:,id)
      ! + &NLTEspec%AtomOpac%sca_c(:,id) * NLTEspec%Jc(:,id)
      
      ! continuum source function
      Snu_c = (NLTEspec%AtomOpac%eta_c(:,id)) / (NLTEspec%AtomOpac%chi_c(:,id) + 1d-300)
    end if
    if (minval(Snu) < 0) then
      write(*,*) "eta (max/min)", maxval(NLTEspec%AtomOpac%eta(:,id)), minval(NLTEspec%AtomOpac%eta(:,id))
      write(*,*) "chi (max/min)", maxval(NLTEspec%AtomOpac%chi(:,id)), minval(NLTEspec%AtomOpac%chi(:,id))
      call Warning("Snu negative")   
    end if
    !In ray-traced map (which is used this subroutine) we integrate from the observer
    ! to the source(s), but we actually get the outgoing intensity/flux. Direct
    !intgreation of the RTE from the observer to the source result in having the
    ! inward flux and intensity. Integration from the source to the observer is done
    ! when computing the mean intensity: Sum_ray I*exp(-dtau)+S*(1-exp(-dtau))
    ! NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dtau*Snu*dexp(-tau)
    ! NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dtau_c*Snu_c*dexp(-tau_c)
    NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + &
    						 dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
    NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + &
                             dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c

!      !!surface superieure ou inf, not used with AL-RT
     facteur_tau = 1.0
!      if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
!      if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

     ! Mise a jour profondeur optique pour cellule suivante
     ! dtau = chi * ds
     tau = tau + dtau * facteur_tau
     tau_c = tau_c + dtau_c

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_I
  
 SUBROUTINE INTEG_RAY_LINE_Z(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  double precision, intent(in) :: u,v,w
  double precision, intent(in) :: x,y,z
  logical, intent(in) :: labs
  double precision :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  double precision, dimension(NLTEspec%Nwaves) :: Snu, Snu_c, Istar
  double precision, dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chiI
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  double precision :: facteur_tau
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp
  tau_c = 0.0_dp


  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0
  NLTEspec%StokesV(:,iray,id) = 0d0
  NLTEspec%StokesQ(:,iray,id) = 0d0
  NLTEspec%StokesU(:,iray,id) = 0d0


  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  infinie : do 
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     lcellule_non_vide = (atmos%icompute_atomRT(icell)>0)
    else
     lcellule_non_vide=.false.
    endif

    
    if (test_exit_grid(icell, x0, y0, z0)) RETURN
    if (lintersect_stars) then
      if (icell == icell_star) then
       CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,Istar)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + Istar*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + Istar*dexp(-tau_c)
       RETURN
       endif
    endif

    nbr_cell = nbr_cell + 1

    previous_cell = 0

    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, &
                     l, l_contrib, l_void_before)


   
    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs)
     l_contrib = l_contrib * AU_to_m
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m

     CALL initAtomOpac(id,(labs.and.(nbr_cell==1)))
     CALL NLTEopacity(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l, (labs.and.(nbr_cell==1)))

     if (lstore_opac) then
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chiI = (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp)
      dtau(:)   = l_contrib * chiI
      dtau_c(:) = l_contrib * NLTEspec%AtomOpac%Kc(icell,:,1)
      Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
               NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / chiI

      Snu_c = (NLTEspec%AtomOpac%jc(icell,:)) / &
            (NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp)
     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chiI = (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id) + tiny_dp)
      dtau(:)   =  l_contrib * chiI
      dtau_c(:) = l_contrib * NLTEspec%AtomOpac%chi_c(:,id)


      Snu = (NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / chiI


      Snu_c = (NLTEspec%AtomOpac%eta_c(:,id)) / (NLTEspec%AtomOpac%chi_c(:,id) + tiny_dp)
    end if
    !continuum not affected by polarisation yet
     NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + &
                             exp(-tau_c) * (1.0_dp - exp(-dtau_c)) * Snu_c
                             
    !Correct line source fonction from polarisation
    !explicit product of Seff = S - (K/chiI - 1) * I
    !Is there a particular initialization to do for polarisation ?
    !because it will be zero at the first place
     Snu(:) = Snu(:) + &
     	-NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI *  NLTEspec%StokesQ(:,iray,id) +&
     	-NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI *  NLTEspec%StokesU(:,iray,id) +&
     	-NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI *  NLTEspec%StokesV(:,iray,id)

     NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu
                             
     NLTEspec%StokesQ(:,iray,id) = NLTEspec%StokesQ(:,iray,id) + &
                             exp(-tau) * (1.0_dp - exp(-dtau)) * (&
			NLTEspec%AtomOpac%etaQUV_p(:,1,id) /chiI - &
			NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI * NLTEspec%I(:,iray,id) - &
			NLTEspec%AtomOpac%rho_p(:,3,id)/chiI * NLTEspec%StokesU(:,iray,id) + &
			NLTEspec%AtomOpac%rho_p(:,2,id)/chiI * NLTEspec%StokesV(:,iray, id) &
     ) !end Snu_Q
     NLTEspec%StokesU(:,iray,id) = NLTEspec%StokesU(:,iray,id) + &
                             exp(-tau) * (1.0_dp - exp(-dtau)) * (&
			NLTEspec%AtomOpac%etaQUV_p(:,2,id) /chiI - &
			NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI * NLTEspec%I(:,iray,id) + &
			NLTEspec%AtomOpac%rho_p(:,3,id)/chiI * NLTEspec%StokesQ(:,iray,id) - &
			NLTEspec%AtomOpac%rho_p(:,1,id)/chiI * NLTEspec%StokesV(:,iray, id) &
     ) !end Snu_U
     NLTEspec%StokesV(:,iray,id) = NLTEspec%StokesV(:,iray,id) + &
                             exp(-tau) * (1.0_dp - exp(-dtau)) * (&
			NLTEspec%AtomOpac%etaQUV_p(:,3,id) /chiI - &
			NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI * NLTEspec%I(:,iray,id) - &
			NLTEspec%AtomOpac%rho_p(:,2,id)/chiI * NLTEspec%StokesQ(:,iray,id) + &
			NLTEspec%AtomOpac%rho_p(:,1,id)/chiI * NLTEspec%StokesU(:,iray, id) &
     ) !end Snu_V

     facteur_tau = 1.0

     tau = tau + dtau * facteur_tau
     tau_c = tau_c + dtau_c

    end if
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_Z
  
  SUBROUTINE FLUX_PIXEL_LINE(&
         id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! -------------------------------------------------------------- !
   ! Computes the flux emerging out of a pixel.
   ! see: mol_transfer.f90/intensite_pixel_mol()
  ! -------------------------------------------------------------- !

   integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
   double precision, dimension(3), intent(in) :: pixelcorner,dx,dy
   double precision, intent(in) :: pixelsize,u,v,w
   integer, parameter :: maxSubPixels = 32
   double precision :: x0,y0,z0,u0,v0,w0
   double precision, dimension(NLTEspec%Nwaves) :: Iold, nu, I0, I0c
   double precision, dimension(3,NLTEspec%Nwaves) :: QUV
   double precision, dimension(3) :: sdx, sdy
   double precision:: npix2, diff
   double precision, parameter :: precision = 1.e-2
   integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter
   logical :: lintersect, labs

   labs = .false.
   ! reset local Fluxes
   !I0c = 0d0
   !I0 = 0d0
   !if (atmos%magnetized .and. PRT_SOLUTION == "FULL_STOKES") QUV(:,:) = 0d0

   ! Ray tracing : on se propage dans l'autre sens
   u0 = -u ; v0 = -v ; w0 = -w

   ! le nbre de subpixel en x est 2^(iter-1)
   subpixels = 1
   iter = 1

   infinie : do ! Boucle infinie tant que le pixel n'est pas converge
     npix2 =  real(subpixels)**2
     Iold = I0
     I0 = 0d0
     I0c = 0d0
     if (atmos%magnetized.and. PRT_SOLUTION == "FULL_STOKES") QUV(:,:) = 0d0
     ! Vecteurs definissant les sous-pixels
     sdx(:) = dx(:) / real(subpixels,kind=dp)
     sdy(:) = dy(:) / real(subpixels,kind=dp)

     iray = 1 ! because the direction is fixed and we compute the flux emerging
               ! from a pixel, by computing the Intensity in this pixel

     ! L'obs est en dehors de la grille
     ri = 2*n_rad ; zj=1 ; phik=1

     ! Boucle sur les sous-pixels qui calcule l'intensite au centre
     ! de chaque sous pixel
     do i = 1,subpixels
        do j = 1,subpixels
           ! Centre du sous-pixel
           x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
           y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
           z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)
           ! On se met au bord de la grille : propagation a l'envers
           !write(*,*) x0, y0, z0
           CALL move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
           !write(*,*) i, j, lintersect, labs, n_cells, icell
             CALL INTEG_RAY_LINE(id, icell, x0,y0,z0,u0,v0,w0,iray,labs)
             I0 = I0 + NLTEspec%I(:,iray,id)
             I0c = I0c + NLTEspec%Ic(:,iray,id)
             if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
             	QUV(3,:) = QUV(3,:) + NLTEspec%STokesV(:,iray,id)
             	QUV(1,:) = QUV(1,:) + NLTEspec%STokesQ(:,iray,id)
             	QUV(2,:) = QUV(2,:) + NLTEspec%STokesU(:,iray,id)
             end if
           !else !Outside the grid, no radiation flux
           endif
        end do !j
     end do !i

     I0 = I0 / npix2
     I0c = I0c / npix2
     if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") QUV(:,:) = QUV(:,:) / npix2

     if (iter < n_iter_min) then
        ! On itere par defaut
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! On arrete pour pas tourner dans le vide
        ! write(*,*) "Warning : converging pb in ray-tracing"
        ! write(*,*) " Pixel", ipix, jpix
        exit infinie
     else
        ! On fait le test sur a difference
        diff = maxval( abs(I0 - Iold) / (I0 + 1e-300_dp) )
        if (diff > precision ) then
           ! On est pas converge
           subpixels = subpixels * 2
        else
           exit infinie
        end if
     end if ! iter
     iter = iter + 1
   end do infinie

  !Prise en compte de la surface du pixel (en sr)

  !Because unit(Bnu)=unit(I) = W/m2/Hz/sr = unit(Blambda*lambda**2/c)
  !if we want I in W/m2 --> I*nu
  nu = 1d0 !c_light / NLTEspec%lambda * 1d9 !to get W/m2 instead of W/m2/Hz !in Hz
  ! Flux out of a pixel in W/m2/Hz
  I0 = nu * I0 * (pixelsize / (distance*pc_to_AU) )**2
  I0c = nu * I0c * (pixelsize / (distance*pc_to_AU) )**2
  if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then 
   QUV(1,:) = QUV(1,:) * nu * (pixelsize / (distance*pc_to_AU) )**2
   QUV(2,:) = QUV(2,:) * nu * (pixelsize / (distance*pc_to_AU) )**2
   QUV(3,:) = QUV(3,:) * nu * (pixelsize / (distance*pc_to_AU) )**2
  end if

  ! adding to the total flux map.
  if (RT_line_method==1) then
    NLTEspec%Flux(:,1,1,ibin,iaz) = NLTEspec%Flux(:,1,1,ibin,iaz) + I0
    NLTEspec%Fluxc(:,1,1,ibin,iaz) = NLTEspec%Fluxc(:,1,1,ibin,iaz) + I0c
    if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then 
     NLTEspec%F_QUV(1,:,1,1,ibin,iaz) = QUV(1,:) !U
     NLTEspec%F_QUV(2,:,1,1,ibin,iaz) = QUV(2,:) !Q
     NLTEspec%F_QUV(3,:,1,1,ibin,iaz) = QUV(3,:) !V
    end if
  else
    NLTEspec%Flux(:,ipix,jpix,ibin,iaz) = I0
    NLTEspec%Fluxc(:,ipix,jpix,ibin,iaz) = I0c
    if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then 
     NLTEspec%F_QUV(1,:,ipix,jpix,ibin,iaz) = QUV(1,:) !U
     NLTEspec%F_QUV(2,:,ipix,jpix,ibin,iaz) = QUV(2,:) !Q
     NLTEspec%F_QUV(3,:,ipix,jpix,ibin,iaz) = QUV(3,:) !V
    end if
  end if

  RETURN
  END SUBROUTINE FLUX_PIXEL_LINE

 SUBROUTINE EMISSION_LINE_MAP(ibin,iaz)
 ! -------------------------------------------------------------- !
  ! Line emission map in a given direction n(ibin,iaz),
  ! using ray-tracing.
  ! if only one pixel it gives the total Flux.
  ! See: emission_line_map in mol_transfer.f90
 ! -------------------------------------------------------------- !
  integer, intent(in) :: ibin, iaz !define the direction in which the map is computed
  double precision :: x0,y0,z0,l,u,v,w

  double precision, dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  double precision, dimension(3,nb_proc) :: pixelcorner
  double precision:: taille_pix, nu
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36
  integer, parameter :: n_ray_star = 1000
  double precision, dimension(n_rad_RT) :: tab_r
  double precision:: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda, lM, lR, lB, ll
  logical :: lresolved

  write(*,*) "incl (deg) = ", tab_RT_incl(ibin), "azimuth (deg) = ", tab_RT_az(iaz)

  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin)
  uvw = (/u,v,w/) !vector position

  ! Definition des vecteurs de base du plan image dans le repere universel
  ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad),sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Vecteur x image avec PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin ! in AU

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Coin en bas gauche de l'image
  Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)

  if (RT_line_method==1) then !log pixels
    n_iter_min = 1
    n_iter_max = 1

    ! dx and dy are only required for stellar map here
    taille_pix = (map_size/zoom)  ! en AU
    dx(:) = x_plan_image * taille_pix
    dy(:) = y_plan_image * taille_pix

    i = 1
    j = 1
    lresolved = .false.

    rmin_RT = max(w*0.9_dp,0.05_dp) * Rmin
    rmax_RT = 2.0_dp * Rmax

    tab_r(1) = rmin_RT
    fact_r = exp( (1.0_dp/(real(n_rad_RT,kind=dp) -1))*log(rmax_RT/rmin_RT) )

    do ri_RT = 2, n_rad_RT
      tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
    enddo

    fact_A = sqrt(pi * (fact_r - 1.0_dp/fact_r)  / n_phi_RT )

    ! Boucle sur les rayons d'echantillonnage
    !$omp parallel &
    !$omp default(none) &
    !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
    !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
    !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz)
    id =1 ! pour code sequentiel

    if (l_sym_ima) then
      cst_phi = pi  / real(n_phi_RT,kind=dp)
    else
      cst_phi = deux_pi  / real(n_phi_RT,kind=dp)
    endif

     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1

        r = tab_r(ri_RT)
        taille_pix =  fact_A * r ! racine carree de l'aire du pixel

        do phi_RT=1,n_phi_RT ! de 0 a pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)
           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image
            ! C'est le centre en fait car dx = dy = 0.
           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, &
                      i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        end do !j
     end do !i
     !$omp end do
     !$omp end parallel
  else !method 2
     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
     lresolved = .true.

     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     if (l_sym_ima) then
        npix_x_max = npix_x/2 + modulo(npix_x,2)
     else
        npix_x_max = npix_x
     endif

     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y) &
     !$omp shared(n_iter_min,n_iter_max,ibin,iaz)

     ! loop on pixels
     id = 1 ! pour code sequentiel
     n_iter_min = 1 ! 3
     n_iter_max = 1 ! 6
     !$omp do schedule(dynamic,1)
     do i = 1,npix_x_max
        !$ id = omp_get_thread_num() + 1
        do j = 1,npix_y
           ! Coin en bas gauche du pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, &
                      i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        end do !j
     end do !i
     !$omp end do
     !$omp end parallel
  end if

!    write(*,*) " -> Adding stellar flux" !Stellar flux by MCFOST is in W/m2 = Blambda*lambda
! ! 									 !Divide by nu to obtain W/m2/Hz (Bnu = Blambda*c/nu**2)
! !   !! This is slow in my implementation actually
!   do lambda = 1, NLTEspec%Nwaves
!    nu = c_light / NLTEspec%lambda(lambda) * 1d9 !if NLTEspec%Flux in W/m2 set nu = 1d0 Hz
!                                              !else it means that in FLUX_PIXEL_LINE, nu
!                                              !is 1d0 (to have flux in W/m2/Hz)
!    CALL compute_stars_map(lambda, u, v, w, taille_pix, dx, dy, lresolved)
!    NLTEspec%Flux(lambda,:,:,ibin,iaz) =  NLTEspec%Flux(lambda,:,:,ibin,iaz) +  &
!                                          stars_map(:,:,1) / nu
!    NLTEspec%Fluxc(lambda,:,:,ibin,iaz) = NLTEspec%Fluxc(lambda,:,:,ibin,iaz) + &
!                                          stars_map_cont(:,:,1) / nu
!   end do

 RETURN
 END SUBROUTINE EMISSION_LINE_MAP

 SUBROUTINE Atomic_transfer()
 ! --------------------------------------------------------------------------- !
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
 ! --------------------------------------------------------------------------- !
  integer :: atomunit = 1, nact
  integer :: icell
  integer :: ibin, iaz, Nrayone = 1
  real(kind=dp) :: dumm
  character(len=20) :: ne_start_sol = "H_IONISATION"
  character(len=20)  :: newPRT_SOLUTION = "FULL_STOKES"
  logical :: lwrite_waves = .false.
   
!   
!   integer :: i, j
!   double precision, dimension(:), allocatable :: test_col
!   double precision, dimension(:,:), allocatable :: test_col2, Cji

 ! -------------------------------INITIALIZE AL-RT ------------------------------------ !
  Profile => IProfile
  INTEG_RAY_LINE => INTEG_RAY_LINE_I
  optical_length_tot => atom_optical_length_tot !not used with the new version of Istar

  if (.not.associated(optical_length_tot)) then
   write(*,*) "pointer optical_length_tot not associated"
   stop
  end if

  if (npix_x_save > 1) then
   RT_line_method = 2 ! creation d'une carte avec pixels carres
   npix_x = npix_x_save ; npix_y = npix_y_save
  else
   RT_line_method = 1 !pixels circulaires
  end if

  if (PRT_SOLUTION == "FULL_STOKES") then
   CALL Warning(" Full Stokes solution not allowed. Stokes polarization not handled in SEE yet.")
   PRT_SOLUTION = "FIELD_FREE"
  end if

  if (atmos%magnetized .and. PRT_SOLUTION /= "NO_STOKES") then
   if (PRT_SOLUTION == "FIELD_FREE") newPRT_SOLUTION = "FULL_STOKES"
   CALL adjustStokesMode(PRT_SOLUTION)
  end if

  !read in dust_transfer.f90 in case of -pluto.
  !mainly because, RT atomic line and dust RT are not decoupled
  !move elsewhere, in the model reading/definition?
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
!! ----------------------- Read Model ---------------------- !!
  ! -> only for not Voronoi
  !later use the same density and T as defined in dust_transfer.f90
  !apply a correction for atomic line if needed.
  !if not flag atom in parafile, never enter this subroutine
  lmcfost_star = .false.
  if (.not.lpluto_file) then 
   !test
   !CALL spherical_star(); lmcfost_star = .true.
   CALL magneto_accretion_model()  
  end if
!! --------------------------------------------------------- !!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------- READATOM and INITIZALIZE POPS ------------------------------ !

  !Read atomic models and allocate space for n, nstar
  ! on the whole grid space.
  CALL readAtomicModels(atomunit)
  if (atmos%NactiveAtoms > 0) then 
   atmos%Nrays = 100 !maximum number of rays allowed
  end if
  
!   test Johnson
!   atmos%T(1) = 2000.
!   allocate(test_col(atmos%Atoms(1)%ptr_atom%Nlevel),&
!   test_col2(atmos%Atoms(1)%ptr_atom%Nlevel,atmos%Atoms(1)%ptr_atom%Nlevel),&
!   Cji(atmos%Atoms(1)%ptr_atom%Nlevel,atmos%Atoms(1)%ptr_atom%Nlevel))
!   test_col = 0d0
!   test_col2 = 0d0
!   CALL Johnson_CI(1, test_col)
!   CALL Johnson_CE(1, test_col2)
!   Cji = Collision_Hydrogen(1)
!   do j=1,Hydrogen%Nlevel
!    write(*,*) j, (Cji(j,i), i=1,Hydrogen%Nlevel)
!   end do
!   stop
  ! if the electron density is not provided by the model or simply want to
  ! recompute it
  ! if atmos%ne given, but one want to recompute ne
  ! else atmos%ne not given, its computation is mandatory
  ! if NLTE pops are read in readAtomicModels, H%n(Nlevel,:) and atom%n
  !can be used for the electron density calculation.
  if (.not.atmos%calc_ne) atmos%calc_ne = lsolve_for_ne
  if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
  if (atmos%calc_ne) then
   if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
   if (.not.lsolve_for_ne) write(*,*) "Solving for electron density"
   write(*,*) " Starting solution : ", ne_start_sol
   if ((ne_start_sol == "NE_MODEL") .and. (atmos%calc_ne)) then
    write(*,*) "WARNING, ne from model is presumably 0 (or not given)."
    write(*,*) "Changing ne starting solution NE_MODEL in H_IONISATION"
    ne_start_sol = "H_IONISATION"
   end if
   CALL SolveElectronDensity(ne_start_sol)
   CALL writeElectron()
  end if

! atmos%Nspace = 82
! CALL FALC_MODEL()
! atmos%icompute_atomRT = 0
! atmos%icompute_atomRT(1:82) = 1 
  CALL setLTEcoefficients () !write pops at the end because we probably have NLTE pops also
  !set Hydrogen%n = Hydrogen%nstar for the background opacities calculations.
! stop
 ! ------------------------------------------------------------------------------------ !
 ! ------------- INITIALIZE WAVELNGTH GRID AND BACKGROUND OPAC ------------------------ !
 ! ------------------------------------------------------------------------------------ !
  if (ltab_wavelength_image) NLTEspec%write_wavelength_grid = .true.
  !otherwise not necessary to write them, because they are in flux.fits
  CALL initSpectrum(vacuum_to_air=lvacuum_to_air)
  if ((atmos%NactiveAtoms > 0) .or.&
      (atmos%NactiveAtoms==0 .and. .not.ltab_wavelength_image)) then
    if (lstore_opac) CALL storeBackground()
  end if
 ! ------------------------------------------------------------------------------------ !
 ! --------------------- ADJUST MCFOST FOR STELLAR MAP AND VELOCITY ------------------- !
  ! ----- ALLOCATE SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ------!
  CALL reallocate_mcfost_vars()
  ! --- END ALLOCATING SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT --!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- NLTE LOOP -------------------------------------- !
  !The BIG PART IS HERE
  if (atmos%Nactiveatoms > 0) then 
   CALL NLTEloop()
  if (lstore_opac) then
   !Recompute Background opac if ne and Hydrogen%n have changed
   ! if ne all storeBackground + LTEpops before
   ! if not, only nH min bound-free
!    if (lstore_opac) then 
!     if (.not.ltab_wavelength_image)  CALL storeBackground()
!    end if  
  open(unit=12, file="Snu_nlte.dat",status="unknown")
  icell = 1 !select cell
  do while (atmos%icompute_atomRT(icell) /= 1) !search the next non dark non transparent cell
   if (icell==atmos%Nspace) then
    icell = 1
    exit
   end if
   icell = icell + 1
  end do
  CALL initAtomOpac(1,.false.)
  CALL BackgroundLines(1,icell, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, dumm)
  CALL NLTEopacity(1, icell, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, dumm, .false.)
  do nact=1, NLTEspec%Nwaves
! LTE
!   write(12,"(5E)") NLTEspec%lambda(nact), NLTEspec%AtomOpac%eta_p(nact,1)+NLTEspec%AtomOpac%jc(icell,nact),&
!    NLTEspec%AtomOpac%chi_p(nact,1)+NLTEspec%AtomOpac%Kc(icell,nact,1), &
! NLTEspec%AtomOpac%eta(nact,1), NLTEspec%AtomOpac%chi(nact, 1)
! NLTE
  write(12,"(5E)") NLTEspec%lambda(nact), NLTEspec%AtomOpac%eta(nact,1), NLTEspec%AtomOpac%chi(nact, 1), &
    NLTEspec%AtomOpac%eta_p(nact,1)+NLTEspec%AtomOpac%jc(icell,nact), NLTEspec%AtomOpac%chi_p(nact,1)+NLTEspec%AtomOpac%Kc(icell,nact,1)
  end do
  close(12) 
  !stop
  end if !test Sf
  end if !NLTE loop
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------- WRITE CONVERGED POPULATIONS ------------------------------ !
 ! ------------------------------------------------------------------------------------ !
  do icell=1,atmos%NactiveAtoms
   CALL writePops(atmos%Atoms(icell)%ptr_atom)
  end do
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- MAKE IMAGES ------------------------------------ !
  if (ltab_wavelength_image) then
   !Check smarter to deallocate/reallocated NLTE wavelength arrays
   CALL initSpectrumImage() !deallocate waves arrays/ define a new grid
   !shorter than the grid for NLTE / reallocate waves arrays
   !write(*,*) maxval(Hydrogen%continua(1)%alpha), maxval(atmos%Atoms(1)%ptr_atom%continua(1)%alpha)
   !write(*,*) loc(Hydrogen)==loc(Atmos%Atoms(1)%ptr_atom)

   if (lstore_opac) CALL storeBackground() !recompute background opac
   ! TO DO: add NLTE continua and LTE/NLTE lines if possible

   CALL reallocate_mcfost_vars()
  end if
!   write(*,*) atmos%atoms(1)%ptr_atom%lines(1)%i, atmos%atoms(1)%ptr_atom%lines(1)%j, maxval(atmos%atoms(1)%ptr_atom%continua(1)%alpha)
!   write(*,*) atmos%passiveatoms(1)%ptr_atom%lines(1)%i, atmos%passiveatoms(1)%ptr_atom%lines(1)%j, maxval(atmos%passiveatoms(1)%ptr_atom%continua(1)%alpha)
!   write(*,*) Hydrogen%lines(1)%i, Hydrogen%lines(1)%j, maxval(Hydrogen%continua(1)%alpha)

  if (atmos%Nrays /= Nrayone) CALL reallocate_rays_arrays(Nrayone)
  !Take into account the wavelength grid for images 
  !Except ds(Nray,Nproc) and Polarized arrays which are: 1) deallocated if PRT_SOL
  !is FIELD_FREE (or not allocated if no magnetic field), reallocated if FULL_STOKES.
  !And, atmos%Nrays = Nrayone if they are different.
  !This in order to reduce the memory of arrays in the map calculation
  
  !Add also a test, if the old solution is the same we do not beed that
  !To do
  if (atmos%magnetized .and. PRT_SOLUTION /= "NO_STOKES") then
   if (PRT_SOLUTION == "FIELD_FREE") PRT_SOLUTION = newPRT_SOLUTION
    CALL adjustStokesMode(PRT_SOLUTION)
  end if

  write(*,*) "Computing emission flux map..."
  !Use converged NLTEOpac
  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
       CALL EMISSION_LINE_MAP(ibin,iaz)
     end do
  end do
  CALL WRITE_FLUX()
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------- CLEANING ------------------------------------------ !
 !close file after NLTE loop
!Temporary: because we kept in memory, so file is closed earlier
!  do nact=1,atmos%Nactiveatoms
!   CALL closeCollisionFile(atmos%ActiveAtoms(nact)%ptr_atom) !if opened
!  end do
 !CALL WRITEATOM() !keep C in memory for that ?
 CALL freeSpectrum() !deallocate spectral variables
 CALL free_atomic_atmos()
 NULLIFY(optical_length_tot, Profile)

 RETURN
 END SUBROUTINE
 ! ------------------------------------------------------------------------------------ !

 SUBROUTINE NLTEloop() !for all active atoms
 ! -------------------------------------------------------- !
  ! CASE : MALI
  ! CASE : HOGEREIJDE
  !      -> should be simpler (in reading and implementing)
  !			than MALI
  !		 -> less accurate than MALI
  !		 -> faster than MALI
  !		 -> should have a lower convergence speed ?
 ! -------------------------------------------------------- !
 
#include "sprng_f.h"

  integer, parameter :: n_rayons_start = 100 ! l'augmenter permet de reduire le tps de l'etape 2 qui est la plus longue
  integer, parameter :: n_rayons_start2 = 100
  integer, parameter :: n_iter2_max = 3
  integer :: n_rayons_max = 0!n_rayons_start2 * (2**(n_iter2_max-1))
  integer :: n_level_comp
  real, parameter :: precision_sub = 1.0e-3 !1e-2
  real, parameter :: precision = 1.0e-1
  integer :: etape, etape_start, etape_end, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  integer, dimension(nb_proc) :: max_n_iter_loc

  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged

  real :: rand, rand2, rand3, fac_etape

  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dN, dN1

  logical :: labs, disable_subit, Ng_acceleration = .false., iterate_ne = .true.
  integer :: atomunit = 1, nact
  integer :: icell
  integer :: Nlevel_total = 0, NmaxLevel, ilevel, max_sub_iter
  character(len=20) :: ne_start_sol = "NE_MODEL" !iterate from the starting guess
  double precision, allocatable, dimension(:,:,:) :: pop_old, pop
  double precision, allocatable, dimension(:,:,:) :: gpop_old
  real(kind=dp), dimension(3, atmos%Nrays, nb_proc) :: xyz0, uvw0
  type (AtomType), pointer :: atom

  write(*,*) "   -> Solving for kinetic equations for ", atmos%Nactiveatoms, " atoms"


  if (allocated(ds)) deallocate(ds)
  allocate(ds(atmos%Nrays,NLTEspec%NPROC))
  ds = 0d0 !meters
  
  n_rayons_max = atmos%Nrays
  labs = .true. !to have ds at cell icell
  id = 1
  etape_start = 1
  etape_end = 1
  disable_subit = .true. !set to true to avoid subiterations over the emissivity
  max_sub_iter = 25
  atmos%nLTE_methode="HOGEREIJDE" !force Hogereijde, MALI not okay yet
 ! ----------------------------  INITIAL POPS------------------------------------------ !
  ! CALL initialSol()
  ! if OLD_POPULATIONS, the init is done at the reading
  ! for now initialSol() is replaced by this if loop on active atoms
  NmaxLevel = 0
  do nact=1,atmos%Nactiveatoms
     atom => atmos%ActiveAtoms(nact)%ptr_atom
     !Now we can set it to .true. The new background pops or the new ne pops
     !will used the H%n
     atom%NLTEpops = .true.
     Nlevel_total = Nlevel_total + atom%Nlevel**2
     NmaxLevel = max(NmaxLevel, atom%Nlevel)
     write(*,*) "Setting initial solution for active atom ", atom%ID, atom%active
     atom%n(:,:) = 1d0*atom%nstar(:,:)

!!Check collisionRate new before completely removing atom%Ckij
!        allocate(atom%Ckij(atmos%Nspace,atom%Nlevel*atom%Nlevel))
!        !open collision file
!        atom%Ckij = 0d0
     CALL openCollisionFile(atom) !closed at the end of the NLTE, it is not nice to do that
       								!but cheap in memory. Cause problem if parallel or
       								!run on a machine. Extra I/O overhead expected

     !CALL allocNetCoolingRates(atmos%ActiveAtoms(nact)%ptr_atom)
     !!Allocate Ng structure of all levels of all atoms, updated at each cell
     !!-> check allocation in statequil
     !!if (Ng_acceleration) CALL initNg(atom%Nlevel, 2, 3, 6,atom%n(:,1), atom%Ngs)
     allocate(atom%Gamma(atom%Nlevel, atom%Nlevel,NLTEspec%NPROC))

	 CALL Keep_collision_lines(atom) !an array containing the lines in file read from atom%offset_col to END
	 !it avoids reading in the file, instead it reads an array (small)
	 CALL closeCollisionFile(atom)
!      do icell=1,atmos%Nspace
!         if (atmos%lcompute_atomRT(icell)) CALL CollisionRate(icell, atom) !open and allocated in LTE.f90
!       !try keeping in memory until better collision routine !
!   	 end do	
!   	 write(*,*) "Fill collision rates for that atom.."
!   	 CALL closeCollisionFile(atom) !if opened
    !!CALL writeAtomData(atmos%ActiveAtoms(nact)%ptr_atom) !to move elsewhere
  	 !deallocate(atom%C) !not used anymore if stored on RAM
     atom => NULL()
  end do
  
  ! Temporary keep collision on RAM, BEWARE IT CAN BE LARGE, need a better collision routine
  ! which stores the required data to compute on the fly the Cij, instead of reading it
  ! each cell
!    if (real(Nlevel_total*n_cells)/(1024**3) < 1.) then
!     write(*,*) "Keeping", real(Nlevel_total*n_cells)/(1024**2), " MB of memory", &
!    	 " for Collisional matrix."
!    else
!     write(*,*) "Keeping", real(Nlevel_total*n_cells)/(1024**3), " GB of memory", &
!    	 " for Collisional matrix."
!    end if
  !end replacing initSol()
 ! ------------------------------------------------------------------------------------ !
  allocate(pop_old(atmos%NactiveAtoms, NmaxLevel, NLTEspec%NPROC)); pop_old = 0d0
  allocate(pop(atmos%NactiveAtoms, NmaxLevel, NLTEspec%NPROC)); pop = 0d0
  allocate(gpop_old(atmos%NactiveAtoms, Nmaxlevel,n_cells)); gpop_old = 0d0
  
  CALL alloc_wlambda() !only for lines actually
          
  SELECT CASE (atmos%nLTE_methode)
   CASE ("MALI")
   
   		lfixed_rays = .true.
  		n_rayons = n_rayons_max
  		iray_start = 1
  		lprevious_converged = .false.
  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0
        write(*,*) " Implementing ..."
        stop
 
    CASE ("HOGEREIJDE")
     do etape=etape_start, etape_end

      if (etape==1) then !two rays
        lfixed_rays=.true.
        n_rayons = 2
        iray_start = 1
        fac_etape = 1.
        lprevious_converged = .false.
      else if (etape==2) then
  		lfixed_rays = .true.
  		n_rayons = n_rayons_start
  		iray_start = 1
  		lprevious_converged = .false.
  		fac_etape = 1. !1d-2
  	  else
  	    CALL ERROR("etape unkown")
  	  end if
  	  
  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0 

        do while (.not.lconverged)
        	n_iter = n_iter + 1    
            write(*,*) " -> Iteration #", n_iter, " Step #", etape
            	
  			if (lfixed_rays) then
    			stream = 0.0
    			do i=1,nb_proc
     				stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    			end do
 			 end if

            max_n_iter_loc = 0
            
 			!$omp parallel &
            !$omp default(none) &
            !$omp private(icell, id, atom) &
            !$omp shared(atmos, gpop_old,nact)
            !$omp do schedule(static,1)
            do icell=1, atmos%Nspace
   				!$ id = omp_get_thread_num() + 1
   				if (atmos%icompute_atomRT(icell)>0) then
            	do nact=1, atmos%NactiveAtoms
            	    atom => atmos%ActiveAtoms(nact)%ptr_atom
             		gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
             		atom => NULL()
            	end do
            	end if
            end do
        	!$omp end do
        	!$omp end parallel

        	if (iterate_ne .and. n_iter > 1)  then 
        	 write(*,*) "  --> old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
        	 CALL SolveElectronDensity(ne_start_sol)
        	end if
        	

 			!$omp parallel &
            !$omp default(none) &
            !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02) &
            !$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, atom) &
            !$omp shared(xyz0, uvw0, lkeplerian,n_iter,nact) &
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,max_sub_iter) &
            !$omp shared(atmos, n_cells, pop_old, pop, ds,disable_subit, dN, gpop_old) &
            !$omp shared(NLTEspec, lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape)
            !$omp do schedule(static,1)
  			do icell=1, n_cells  !!atom is shared or private???
   				!$ id = omp_get_thread_num() + 1
   				if (atmos%icompute_atomRT(icell)>0) then !nor transparent nor dark
  			        CALL initGamma(id,icell)
     				do iray=iray_start, iray_start-1+n_rayons
      					if (etape==1) then
                        ! Position = milieu de la cellule
                         x0 = r_grid(icell)
                         y0 = 0.0_dp
                         z0 = z_grid(icell)
                         if (lkeplerian) then
                       ! Direction verticale
                          if (iray==1) then
                           w0=1.0_dp
                          else
                           w0=-1.0_dp
                          endif
                          u0 = 0.0_dp
                          v0 = 0.0_dp
                          else !not keplerian, spherical, infall even lmagneto_accr
                          norme = sqrt(x0*x0 + y0*y0 + z0*z0)
                          if (iray==1) then
                           u0 = x0/norme
                           v0 = y0/norme
                           w0 = z0/norme
                          else
                           u0 = -x0/norme
                           v0 = -y0/norme
                           w0 = -z0/norme
                          endif
                         endif !lkeplerian      					 
      					else !etape 2, 3
                   	    ! Position aleatoire dans la cellule
                         rand  = sprng(stream(id))
                         rand2 = sprng(stream(id))
                         rand3 = sprng(stream(id))
                        
                         CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)
                        
                        ! Direction de propagation aleatoire
                         rand = sprng(stream(id))
                         W0 = 2.0_dp * rand - 1.0_dp
                         W02 =  1.0_dp - W0*W0
                         SRW02 = sqrt(W02)
                         rand = sprng(stream(id))
                         ARGMT = PI * (2.0_dp * rand - 1.0_dp)
                         U0 = SRW02 * cos(ARGMT)
                         V0 = SRW02 * sin(ARGMT)  
						end if !etape    					
						xyz0(1,iray,id) = x0; xyz0(2,iray,id) = y0; xyz0(3,iray,id) = z0
						uvw0(1,iray,id) = U0; uvw0(2,iray,id) = V0; uvw0(3,iray,id) = W0

! if (iray==1 .or. iray==n_rayons) &
! write(*,*) iray, icell, "id", id,"rays ", xyz0(:,iray,id), uvw0(:,iray,id)

      					CALL INTEG_RAY_LINE(id, icell, x0, y0, z0, u0, v0, w0, iray, labs)
                        !+add cross_coupling for this cell in this direction in Gamma
                        !then for next ray they are re init
 						CALL fillGamma_Hogereijde(id, icell, iray, n_rayons)
      				end do !iray
    				 if (check_nan_infinity_dp(atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,id))==1) then
    				  write(*,*) "Gamma is nan"
    				  write(*,*) atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,id)
    				  stop
    				 else if(check_nan_infinity_dp(atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,id))==2)  then
    				  write(*,*) "Gamma is inf"
    				  write(*,*) atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,id)
    				  stop
    				 end if
      				!!CALL Gamma_LTE(id,icell) !G(j,i) = C(j,i) + ...
!       				NLTEspec%J(:,id) = 0d0; NLTEspec%Jc(:,id) = 0d0
!       	            do iray=1,n_rayons
!       	             NLTEspec%J(:,id)  NLTEspec%J(:,id) + NLTEspec%I(:,iray,id)/n_rayons
!       	             NLTEspec%Jc(:,id)  NLTEspec%Jc(:,id) + NLTEspec%Ic(:,iray,id)/n_rayons
!       	            end do		
    				!lconverged_loc = disable_subit !to disable sub-it
    		     	n_iter_loc = 0
    				if (disable_subit) then
! atom => atmos%ActiveAtoms(1)%ptr_atom
!     				 write(*,*) "Before update"
!     				 write(*,*) maxval(atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,id)), &
!     				 maxval(atmos%ActiveAtoms(1)%ptr_atom%n(:,icell))
!!!!$omp critical
    				 CALL updatePopulations(id, icell)
!!!$omp end critical
!     				 CALL SEE_atom(id, icell, atom)
!     				 write(*,*) maxval(atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,id)), &
!     				 maxval(atmos%ActiveAtoms(1)%ptr_atom%n(:,icell))
!     				 write(*,*) "after update"
    				 ! stop
    				 lconverged_loc = .true.
!     if (n_iter == 2) then
!     	write(*,*) id, icell, atom%ID, loc(atom)
!     	write(*,*) atom%n(:,icell)
!     	write(*,*) gpop_old(1,:,icell)
!    end if
! atom=>NULL()
    				else
    				 lconverged_loc = .false.
    				!save pops for all active atoms
    				!pop(NactAtom, Nmaxlvel, Nthreads)
    				 do nact=1,atmos%NactiveAtoms
    				  atom => atmos%ActiveAtoms(nact)%ptr_atom
    				  pop(nact,1:atom%Nlevel,id) = atom%n(:,icell)
    				  atom => NULL()
    				 end do
    				end if
    				
     				!!sub iteration on the local emissivity, keeping Idag fixed
     				!!only iterate on cells which are not converged
     				do while (.not.lconverged_loc)
     				!write(*,*) "Starting subit for cell ", icell
       					n_iter_loc = n_iter_loc + 1
       					
                        pop_old(:,:,id) = pop(:,:,id)
						
						!Solve SEE for all atoms
						CALL updatePopulations(id, icell)
						
    				    do nact=1,atmos%NactiveAtoms
    				     atom => atmos%ActiveAtoms(nact)%ptr_atom
    				     pop(nact,1:atom%Nlevel,id) = atom%n(:,icell)
    				     atom => NULL()
    				    end do
						
						diff = 0.0 !keeps track of the maximum dpops(dN) among each atom.
								   !for this cell
     					do nact=1,atmos%NactiveAtoms
     					    atom => atmos%ActiveAtoms(nact)%ptr_atom
     						do ilevel=1,atom%Nlevel
     				    		dN = abs((pop_old(nact,ilevel,id)-pop(nact,ilevel,id))/&
     				    			(pop_old(nact,ilevel,id)+1d-300))
     				    		diff = max(diff, dN)
!      				    		if (diff > 1) then !(diff == dN)
!      							 write(*,*) "icell#",icell," level#", ilevel, " sub-it->#", &
!      								n_iter_loc, atom%ID, " dpops = ", dN, diff
!       							end if
     						end do
     						atom => NULL()
     					end do

       					if (diff < precision_sub) then
       						lconverged_loc = .true.
       						!write(*,*) " sub it converged! icell#", icell," maxdpops=", diff
       					else
        					!recompute opacity of this cell., but I need angles and pos...
       						!NLTEspec%I not changed
  			        		CALL initGamma(id,icell)
 							do iray=iray_start, iray_start-1+n_rayons
      							!I unchanged
!     if (iray==1 .or. iray==n_rayons) &
!       	write(*,*) iray, icell, "id", id,"rays sub-it", xyz0(:,iray,id), uvw0(:,iray,id)
							    CALL init_local_field_atom(id, icell, iray, &
							         xyz0(1,iray,id), xyz0(2,iray,id), xyz0(3,iray,id), &
							         uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id))
                               !! +add Xcoupling in Gamma for this cell/ray
 						       CALL fillGamma_Hogereijde(id, icell, iray, n_rayons)
      						end do !iray
      						!!CALL Gamma_LTE(id,icell)
       					end if
       					if (n_iter_loc >= max_sub_iter) then 
       					  if (diff>1) write(*,*) " sub-it not converged after", n_iter_loc, &
       					  	" iterations; diff=", diff
       					  lconverged_loc = .true.
       					end if
     				end do
     	            if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
     			end if !icompute_atomRT
     		end do !icell
        	!$omp end do
        	!$omp end parallel
!if (n_iter == 2) stop
        	!Global convergence criterion
        	!I cannot iterate on unconverged cells because the radiation coming for each cell
        	!depends on the other cell.
        	diff = 0d0
  			cell_loop2 : do icell=1, atmos%Nspace
  				if (atmos%icompute_atomRT(icell)>0) then
						dN = 0d0 !for all levels of all atoms of this cell
     					do nact=1,atmos%NactiveAtoms
     					    atom => atmos%ActiveAtoms(nact)%ptr_atom
     						do ilevel=1,atom%Nlevel
     				    		dN1 = abs((gpop_old(nact,ilevel,icell)-atom%n(ilevel,icell))/&
     				    			(gpop_old(nact,ilevel,icell)+1d-300))
     				    		dN = max(dN1, dN)
 !     						if (dN > 1) & !(dN == dN1)
!      							write(*,*) icell, " **->", atom%ID, " ilevel", ilevel, " dpops = ", dN1
     						end do
     						atom => NULL()
     					end do
     					diff = max(diff, dN) !for all cells
!      					if (diff>1) write(*,*) icell, " dpops(atoms) = ", dN, " dpops(icell) =", diff
     			end if
     		end do cell_loop2

         	!if (maxval(max_n_iter_loc)> max_sub_iter) &
         		write(*,*) maxval(max_n_iter_loc), "sub-iterations"
         	write(*,*) "Relative difference =", diff !at the end of the loop over n_cells
        	write(*,*) "Threshold =", precision*fac_etape

        	if (diff < precision*fac_etape) then
           		if (lprevious_converged) then
            	  lconverged = .true.
           		else
            	  lprevious_converged = .true.
          	    endif
        	else
           		lprevious_converged = .false.
           		if (.not.lfixed_rays) then
              		n_rayons = n_rayons * 2
              		write(*,*) ' -- Increasing number of rays'
             		if (n_rayons > n_rayons_max) then
              			if (n_iter >= n_iter2_max) then
             		 		write(*,*) "Warning : not enough rays to converge !!"
                 			lconverged = .true.
              			end if
              		end if

              ! On continue en calculant 2 fois plus de rayons
              ! On les ajoute a l'ensemble de ceux calcules precedemment
!              iray_start = iray_start + n_rayons

          	   end if
        	end if
        	
	    end do !while
	  end do !over etapes

  CASE DEFAULT
   CALL ERROR("Methode for SEE unknown", atmos%nLTE_methode)
  END SELECT
  
  if (iterate_ne) then
   write(*,*) "  --> old max/min nHmin", maxval(atmos%nHmin), minval(atmos%nHmin,mask=atmos%nHmin>0)
   CALL calc_nHmin()
   do nact=1,atmos%NpassiveAtoms
    atom => atmos%PassiveAtoms(nact)%ptr_atom
    CALL LTEpops(atom,.true.)
    atom => Null()
   end do
   CALL writeElectron(.true.) !the .true. means append, to compare with initial solution.
  end if
 ! -------------------------------- CLEANING ------------------------------------------ !
  ! Remove NLTE quantities not useful now
  
  deallocate(pop_old)
  if (allocated(pop)) deallocate(pop)
  do nact=1,atmos%Nactiveatoms
   atom  => atmos%ActiveAtoms(nact)%ptr_atom
   !!!!CALL closeCollisionFile(atom) 
   if (allocated(atmos%ActiveAtoms(nact)%ptr_atom%gamma)) & !otherwise we have never enter the loop
     deallocate(atmos%ActiveAtoms(nact)%ptr_atom%gamma)
   if (allocated(atmos%ActiveAtoms(nact)%ptr_atom%Ckij)) deallocate(atmos%ActiveAtoms(nact)%ptr_atom%Ckij)
   if (Ng_acceleration) CALL freeNg(atmos%ActiveAtoms(nact)%ptr_atom%Ngs)
   atom => NULL()
  end do
  deallocate(ds)
 ! ------------------------------------------------------------------------------------ !
 RETURN
 END SUBROUTINE NLTEloop
 ! ------------------------------------------------------------------------------------ !
 
 SUBROUTINE adjustStokesMode(Solution)
 !
 !
   character(len=*), intent(in) :: Solution
   
!     if (.not.atmos%magnetized) then
!      Profile => IProfile
!      RETURN
!     end if
      
    if (SOLUTION=="FIELD_FREE") then
      !if (associated(Profile)) Profile => NULL()
      !Profile => Iprofile !already init to that
      if (allocated(NLTEspec%AtomOpac%rho_p)) deallocate(NLTEspec%AtomOpac%rho_p)
      if (allocated(NLTEspec%AtomOpac%etaQUV_p)) deallocate(NLTEspec%AtomOpac%etaQUV_p)
      if (allocated(NLTEspec%AtomOpac%chiQUV_p)) deallocate(NLTEspec%AtomOpac%chiQUV_p)
      if (allocated(NLTEspec%F_QUV)) deallocate(NLTEspec%F_QUV)
      if (allocated(NLTEspec%StokesV)) deallocate(NLTEspec%StokesV, NLTEspec%StokesQ, &
      												NLTEspec%StokesU)
      write(*,*) " Using FIELD_FREE solution for the SEE!"
      CALL Warning("  Set PRT_SOLUTION to FULL_STOKES for images")
      
    else if (SOLUTION=="POLARIZATION_FREE") then
     CALL ERROR("POLARIZATION_FREE solution not implemented yet")
     !NLTE not implemented in integ_ray_line_z and Zprofile always compute the full
     !propagation dispersion matrix, but only the first row is needed
     if (associated(Profile)) Profile => NULL()
     Profile => Zprofile
     if (associated(INTEG_RAY_LINE)) INTEG_RAY_LINE => NULL()
     INTEG_RAY_LINE => INTEG_RAY_LINE_Z
     
     !beware, only I! polarization is neglected hence not allocated
      if (.not.allocated(NLTEspec%AtomOpac%rho_p)) then !same for all, dangerous.
      	allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%etaQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%chiQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))

        write(*,*) " Using POLARIZATION_FREE solution for the SEE!"
     end if
							
    else if (SOLUTION=="FULL_STOKES") then !Solution unchanged here
      if (associated(Profile)) Profile => NULL()
      Profile => Zprofile
      if (associated(INTEG_RAY_LINE)) INTEG_RAY_LINE => NULL()
      INTEG_RAY_LINE => INTEG_RAY_LINE_Z
      !allocate space for Zeeman polarisation
      !only LTE for now
       if (.not.allocated(NLTEspec%StokesQ)) & !all allocated, but dangerous ?
            allocate(NLTEspec%StokesQ(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC), & 
             NLTEspec%StokesU(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC), &
             NLTEspec%StokesV(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC))!, &
             !!NLTEspec%S_QUV(3,NLTEspec%Nwaves))
      NLTEspec%StokesQ=0d0
      NLTEspec%StokesU=0d0
      NLTEspec%StokesV=0d0
      if (.not.allocated(NLTEspec%F_QUV)) &
      	allocate(NLTEspec%F_QUV(3,NLTEspec%Nwaves,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ))
      NLTEspec%F_QUV = 0d0
      if (.not.allocated(NLTEspec%AtomOpac%rho_p)) then !same for all, dangerous.
      	allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%etaQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%chiQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))

        write(*,*) " Using FULL_STOKES solution for the SEE!"
      end if
    else
     CALL ERROR("Error in adjust StokesMode with solution=",Solution)
    end if
 RETURN
 END SUBROUTINE adjustStokesMode
 ! ------------------------------------------------------------------------------------ !

 SUBROUTINE reallocate_mcfost_vars()
  !--> should move them to init_atomic_atmos ? or elsewhere
  !need to be deallocated at the end of molecule RT or its okey ?`
  integer :: la
  CALL init_directions_ray_tracing()
  if (.not.allocated(stars_map)) allocate(stars_map(npix_x,npix_y,3))
  if (.not.allocated(stars_map_cont)) allocate(stars_map_cont(npix_x,npix_y,3))
  n_lambda = NLTEspec%Nwaves
  if (allocated(tab_lambda)) deallocate(tab_lambda)
  allocate(tab_lambda(n_lambda))
  if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
  allocate(tab_delta_lambda(n_lambda))
  tab_lambda = NLTEspec%lambda / MICRON_TO_NM
  tab_delta_lambda(1) = 0d0
  do la=2,NLTEspec%Nwaves
   tab_delta_lambda(la) = tab_lambda(la) - tab_lambda(la-1)
  end do
  if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
  allocate(tab_lambda_inf(n_lambda))
  if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
  allocate(tab_lambda_sup(n_lambda))
  tab_lambda_inf = tab_lambda
  tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
  ! computes stellar flux at the new wavelength points
  CALL deallocate_stellar_spectra()
  ! probably do not deallocate, 'cause maybe dust is in here too.
  if (allocated(kappa)) deallocate(kappa) !do not use it anymore
  !used for star map ray-tracing.
  CALL allocate_stellar_spectra(n_lambda)
  CALL repartition_energie_etoiles()
  !If Vfield is used in atomic line transfer, with lkep or linfall
  !atmos%Vxyz is not used and lmagnetoaccr is supposed to be zero
  !and Vfield allocated and filled in the model definition if not previously allocated
  !nor filled.
  if (allocated(Vfield).and.lkeplerian.or.linfall) then
    write(*,*) "Vfield max/min", maxval(Vfield), minval(Vfield)
  else if (allocated(Vfield).and.lmagnetoaccr) then
   write(*,*) "deallocating Vfield for atomic lines when lmagnetoaccr = .true."
   deallocate(Vfield)
  end if

 RETURN
 END SUBROUTINE reallocate_mcfost_vars
 
 SUBROUTINE reallocate_mcfost_wavelength_arrays()
  !--> should move them to init_atomic_atmos ? or elsewhere
  !need to be deallocated at the end of molecule RT or its okey ?`
  integer :: la
  n_lambda = NLTEspec%Nwaves
  if (allocated(tab_lambda)) deallocate(tab_lambda)
  allocate(tab_lambda(n_lambda))
  if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
  allocate(tab_delta_lambda(n_lambda))
  tab_lambda = NLTEspec%lambda / MICRON_TO_NM
  tab_delta_lambda(1) = 0d0
  do la=2,NLTEspec%Nwaves
   tab_delta_lambda(la) = tab_lambda(la) - tab_lambda(la-1)
  end do
  !unlikely in NLTE because of the different atomic grids
  !but if an image at only one wavelength --> tab_Delta_lambda=0
  if (size(tab_lambda)==1) tab_delta_lambda = tab_delta_lambda + tiny_dp
  if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
  allocate(tab_lambda_inf(n_lambda))
  if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
  allocate(tab_lambda_sup(n_lambda))
  tab_lambda_inf = tab_lambda
  tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
  ! computes stellar flux at the new wavelength points
  CALL deallocate_stellar_spectra()
  CALL allocate_stellar_spectra(n_lambda)
  CALL repartition_energie_etoiles()

 RETURN
 END SUBROUTINE reallocate_mcfost_wavelength_arrays
 
 SUBROUTINE calc_stellar_radiation(N,i_star,x,y,z,u,v,w,Istar)
 ! ---------------------------------------------------------------!
  ! Compute the stellar radiation field and in Istar:
  ! Radiation emerging from the star (at the surface), corrected
  ! by limb darkening.
  !
  ! tab_lambda has to be defined and repartition_energie_etoiles()
 ! -------------------------------------------------------------- !
  use input, only : limb_darkening, mu_limb_darkening
  use Planck, only : uLD, Bplanck
  integer, intent(in) :: N, i_star
  real(kind=dp), dimension(N), intent(out) :: Istar
  real(kind=dp), intent(in) :: u, v, w, x, y, z
  real(kind=dp) :: energie(N), gamma(N)
  real(kind=dp) :: mu, ulimb, LimbDarkening, surface, HC
  integer :: ns
  logical :: lintersect_spot
  
   HC = HPLANCK * CLIGHT / KBOLTZMANN / NM_TO_M
   gamma(:) = 1d0
   !cos(theta) = dot(r,n)/module(r)/module(n)
   mu = abs(x*u + y*v + z*w)/dsqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
   if (real(mu)>1d0) then !to avoid perecision error
    write(*,*) "mu=",mu, x, y, z, u, v, w
    CALL Error(" mu limb > 1!")
   end if

   !Re-norm E_stars(:)*Prob_E_Star(:,i_star) at the stellar surface
   surface = (etoile(i_star)%r * AU_to_Rsun)**2
   
   !1) Get the energy radiated by the star at the stellar surface
   !I need unit of I which is the unit of Bnu = W/m2/Hz/sr
   !Here, E_stars in W/m2. But I is in W/m2/Hz/sr unit of Bnu.
   !E_stars propto Blambda*lambda; remembering, Bnu = Blambda*c/nu**2 = Blambda*lambda**2/c
   ! we get E_stars (in W/m2/Hz) = E_stars / lambda * lambda**2 / c
!    energie(:) = Prob_E_Star(:,i_star) * E_stars(:) * tab_lambda(:) * 1.0e-6 / surface &
!              * 1.35e-12 * (tab_lambda(:) * 1d-6) / C_LIGHT
   !write(*,*) maxval(energie)
   !write(*,*) maxval(E_stars* (tab_lambda(:) * 1.0e-6)/CLIGHT)
   CALL Bplanck(etoile(i_star)%T*1d0, energie) !it is not factorised for test cheks, but can be computed outside loop
   !write(*,*) maxval(energie)
   !stop
   
   !2) Correct with the contrast gamma of a hotter/cooler region if any
   CALL intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect_spot)
   if (lintersect_spot) then
     gamma(:) = (dexp(HC/NLTEspec%lambda(:)/real(etoile(i_star)%T,kind=dp))-1)/&
     			(dexp(HC/NLTEspec%lambda(:)/etoile(i_star)%SurfB(ns)%T)-1)
     !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
   end if

   !3) Apply Limb darkening   
   if (llimb_darkening) then
     LimbDarkening = Interp1D(real(limb_darkening,kind=dp), real(mu_limb_darkening,kind=dp), mu)
     !LimbDarkening = interp(limb_darkening, mu_limb_darkening, mu)
     !pol not included yet
     stop
   else
     !write(*,*) maxval(uLD(real(etoile(i_star)%T,kind=dp))), minval(uLD(real(etoile(i_star)%T,kind=dp)))
     ulimb = 0.0 ! could use BB slope
     LimbDarkening = 1d0 - ulimb*(1d0-mu)
   end if
   Istar(:) = energie(:) * LimbDarkening * gamma(:)
 
 RETURN
 END SUBROUTINE calc_stellar_radiation

END MODULE AtomicTransfer
!         	  icell = 1 !select cell
!   do while (atmos%icompute_atomRT(icell) /= 1) !search the next non dark non transparent cell
!    if (icell==atmos%Nspace) then
!     icell = 1
!     exit
!    end if
!    icell = icell + 1
!   end do
!             atmos%T(icell) = 5000d0
!             Hydrogen%nstar(:,icell)=1d0
!             CALL initGamma(1,icell)
!             write(*,*) Hydrogen%Gamma(:,:,1)
!             Hydrogen%Gamma(:,:,1) = Collision_Hydrogen(icell)
!             write(*,*) atmos%ActiveAtoms(1)%ptr_atom%Gamma(:,:,1)
!         	stop
!         	
!      		atom => atmos%Atoms(1)%ptr_atom
!  			!$omp parallel &
!             !$omp default(none) &
!             !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02) &
!             !$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell) &
!             !$omp shared(xyz0, uvw0, lkeplerian, atom) &
!             !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,max_sub_iter) &
!             !$omp shared(atmos, n_cells, pop_old, pop, ds,disable_subit, dN, gpop_old) &
!             !$omp shared(NLTEspec, lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape)
!             !$omp do schedule(static,1)
!   			do icell=1, n_cells
!      				!$ id = omp_get_thread_num() + 1
!    				if (atmos%icompute_atomRT(icell)>0) then
!    				!!!!!$omp critical
!    				!!!!!atom%Gamma(:,:,id) = CollisionRate(icell, atom)
!    				CALL initGamma(id,icell)
!    				write(*,*) icell, atom%Gamma(:,:,id)
!    				!!!!$omp end critical 
!    				end if
!    			end do
!         	!$omp end do
!         	!$omp end parallel
! 	stop