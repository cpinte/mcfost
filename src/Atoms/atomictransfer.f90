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
 use broad, only 						: Damping, RadiativeDamping !for tests
 use spectrum_type
 use atmos_type
 use readatom
 use lte
 use constant, only 					: MICRON_TO_NM
 use collision, only					: CollisionRate !future deprecation
 use impact
 use solvene
 use statequil_atoms
 use init_solution, only 			: Init_NLTE, free_NLTE_sol, gpop_old, pop_old, flatpops
 use accelerate
 use voigtfunctions
 use writeatom
 use math
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
 use mcfost_env, only : dp
 use constantes, only : tiny_dp, huge_dp

 IMPLICIT NONE

 PROCEDURE(INTEG_RAY_LINE_I), pointer :: INTEG_RAY_LINE => NULL()
 real(kind=dp), dimension(:,:), allocatable :: QUV
 real(kind=dp), allocatable :: S_contrib(:,:,:), S_contrib2(:,:), chil(:), Sl(:)

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
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs !used in NLTE but why?
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, chiI, Istar
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, dtau_ds
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp !go from surface down to the star
  tau_c = 0.0_dp
  chiI = 0d0

  ! Reset, because when we compute the flux map
  ! with emission_line_map, we only use one ray, iray=1
  ! and the same space is used for emergent I.
  ! Therefore it is needed to reset it.
  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0

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
     !atmos%nHtot(icell)>0
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif

    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN

    if (lintersect_stars) then
      if (icell == icell_star) then !this is equivalent to compute_stars_map()
       !if we start at icell, computes the radiation from the star to icell.
       !in particular, the cell icell can be the cell for which we solve the SEE/or
       !a cell at the border of the grid for an image
       Istar(:) = 0d0
       CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,Istar)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + Istar*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + Istar*dexp(-tau_c)
       RETURN
      end if
    endif

    nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
    previous_cell = 0 ! unused, just for Voronoi
    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)

!     if (.not.atmos%lcompute_atomRT(icell)) lcellule_non_vide = .false. !chi and chi_c = 0d0, cell is transparent
!	   this makes the code to break down  --> the atomRT is defined on n_cells
!		but the infinite loop runs over "virtual" cells, which can be be
!		out of atomRT boundary

    !count opacity only if the cell is filled, else go to next cell
    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m !l_contrib in AU
     if ((nbr_cell == 1).and.labs) ds(iray,id) = l * AU_to_m
     !! Compute background opacities for PASSIVE bound-bound and bound-free transitions
     !! at all wavelength points including vector fields in the bound-bound transitions
     ! evaluate Psi operator ?
     eval_operator = (labs .and. (nbr_cell == 1)) !labs if false for images
     											  !so no pb if Nact>0 and we use a different grid
     CALL initAtomOpac(id) !set opac to 0 for this cell and thread id
     if (eval_operator) CALL init_psi_operator(id, iray) !set xcoupling also and eta !

     CALL NLTEopacity(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
     !CALL NLTEopacity_o(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
     !never enter NLTEopacity if no activeatoms
     if (lstore_opac) then !not updated during NLTE loop, just recomputed using initial pops
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chiI(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp
      dtau(:)   = l_contrib * chiI(:)
      dtau_c(:) = l_contrib * (NLTEspec%AtomOpac%Kc(icell,:,1)+ NLTEspec%AtomOpac%chic_nlte(:,id))
      Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
               NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / chiI(:)

      Snu_c = (NLTEspec%AtomOpac%jc(icell,:) + NLTEspec%AtomOpac%etac_nlte(:,id)) / &
            (NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp + NLTEspec%AtomOpac%chic_nlte(:,id))
      if (atmos%coherent_scattering) then
        Snu = Snu + NLTEspec%J(:,icell)*NLTEspec%AtomOpac%Kc(icell,:,2) / chiI(:)
        Snu_c = NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%Jc(:,icell) / (NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp + NLTEspec%AtomOpac%chic_nlte(:,id))
      endif
     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l) !x,y,z,u,v,w,x1,y1,z1
                                !define the projection of the vector field (velocity, B...)
                                !at each spatial location.

      ! Epaisseur optique
      ! chi_p contains both thermal and continuum scattering extinction
      chiI(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id) + tiny_dp
      dtau(:)   =  l_contrib * chiI(:)
      dtau_c(:) = l_contrib * (NLTEspec%AtomOpac%chi_c(:,id) + NLTEspec%AtomOpac%chic_nlte(:,id))

      ! Source function
      ! No dust yet
      ! J and Jc are the mean radiation field for total and continuum intensities
      ! it multiplies the continuum scattering coefficient for isotropic (unpolarised)
      ! scattering. chi, eta are opacity and emissivity for ACTIVE lines.
      Snu = (NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / chiI(:)

      ! continuum source function
      Snu_c = (NLTEspec%AtomOpac%eta_c(:,id) + NLTEspec%AtomOpac%etac_nlte(:,id)) / &
      			 (NLTEspec%AtomOpac%chi_c(:,id) + tiny_dp + NLTEspec%AtomOpac%chic_nlte(:,id))

      if (atmos%coherent_scattering) then
        Snu = Snu + NLTEspec%J(:,icell)*NLTEspec%AtomOpac%sca_c(:,id) / chiI(:)
        Snu_c = NLTEspec%AtomOpac%sca_c(:,id) * NLTEspec%Jc(:,icell) / (NLTEspec%AtomOpac%chi_c(:,id) + tiny_dp + NLTEspec%AtomOpac%chic_nlte(:,id))
      endif

    end if

    !In ray-traced map (which is used this subroutine) we integrate from the observer
    ! to the source(s), but we actually get the outgoing intensity/flux. Direct
    !intgreation of the RTE from the observer to the source result in having the
    ! inward flux and intensity. Integration from the source to the observer is done
    ! when computing the mean intensity: Sum_ray I*exp(-dtau)+S*(1-exp(-dtau))
    ! NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dtau*Snu*dexp(-tau)
    ! NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dtau_c*Snu_c*dexp(-tau_c)
    NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
    NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c


     !here because %chip, chi have to be allocated
     if (eval_operator) then
       dtau_ds = ds(iray,id) * chiI(:)
       CALL calc_psi_operator(id, icell, iray, chiI, dtau_ds)
     end if

     ! Mise a jour profondeur optique pour cellule suivante
     ! dtau = chi * ds
     !tau is updated after Snu, because at icell, we weight Snu with the previous optical depth
     !computed up to icell. For icell+1 we add tau of this icell.
     tau = tau + dtau
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
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, Istar
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chiI
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  real(kind=dp) :: facteur_tau
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
  Istar(:) = 0d0

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
       Istar(:) = 0d0
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

     CALL initAtomOpac(id)
     CALL NLTEopacity(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)!(labs.and.(nbr_cell==1)))

     if (lstore_opac) then
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chiI = (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp)
      dtau(:)   = l_contrib * chiI
      dtau_c(:) = l_contrib * (NLTEspec%AtomOpac%Kc(icell,:,1)+ NLTEspec%AtomOpac%chic_nlte(:,id))
      Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
               NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / chiI

      Snu_c = (NLTEspec%AtomOpac%jc(icell,:) + NLTEspec%AtomOpac%etac_nlte(:,id)) / &
            (NLTEspec%AtomOpac%Kc(icell,:,1) + tiny_dp)
     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chiI = (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id) + tiny_dp)
      dtau(:)   =  l_contrib * chiI
      dtau_c(:) = l_contrib * (NLTEspec%AtomOpac%chi_c(:,id) + NLTEspec%AtomOpac%chic_nlte(:,id))


      Snu = (NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / chiI


      Snu_c = (NLTEspec%AtomOpac%eta_c(:,id) + NLTEspec%AtomOpac%etac_nlte(:,id)) / &
              (NLTEspec%AtomOpac%chi_c(:,id) + tiny_dp + NLTEspec%AtomOpac%chic_nlte(:,id))
    end if
    !continuum not affected by polarisation yet
     NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + exp(-tau_c) * (1.0_dp - exp(-dtau_c)) * Snu_c

    !Correct line source fonction from polarisation
    !explicit product of Seff = S - (K/chiI - 1) * I
    !Is there a particular initialization to do for polarisation ?
    !because it will be zero at the first place
     Snu(:) = Snu(:) -NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI *  NLTEspec%StokesQ(:,iray,id) - &
          NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI *  NLTEspec%StokesU(:,iray,id) - &
          NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI *  NLTEspec%StokesV(:,iray,id)

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
   real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
   real(kind=dp), intent(in) :: pixelsize,u,v,w
   integer, parameter :: maxSubPixels = 32
   real(kind=dp) :: x0,y0,z0,u0,v0,w0
   real(kind=dp), dimension(NLTEspec%Nwaves) :: Iold, I0, I0c !, nu
   real(kind=dp), dimension(3) :: sdx, sdy
   real(kind=dp):: npix2, diff, normF, R0
   real(kind=dp), parameter :: precision = 1.e-2
   integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter
   logical :: lintersect, labs

   labs = .false.

   ! Ray tracing : on se propage dans l'autre sens
   u0 = -u ; v0 = -v ; w0 = -w

   ! le nbre de subpixel en x est 2^(iter-1)
   subpixels = 1
   iter = 1
   diff = 0.

   infinie : do ! Boucle infinie tant que le pixel n'est pas converge
     npix2 =  real(subpixels)**2
     Iold = I0
     I0 = 0d0
     I0c = 0d0
     if (atmos%magnetized.and. PRT_SOLUTION == "FULL_STOKES") QUV(:,:) = 0d0 !move outside
     if (lcontrib_function) S_contrib2(:,:) = 0d0

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

             I0 = I0 + NLTEspec%I(:,iray,id) !/ R0**2
             I0c = I0c + NLTEspec%Ic(:,iray,id)

             if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
             	QUV(3,:) = QUV(3,:) + NLTEspec%STokesV(:,iray,id)
             	QUV(1,:) = QUV(1,:) + NLTEspec%STokesQ(:,iray,id)
             	QUV(2,:) = QUV(2,:) + NLTEspec%STokesU(:,iray,id)
             end if
             if (lcontrib_function) &
                          S_contrib2(:,:) = S_contrib2(:,:) + S_contrib(:,:,id)

           !else !Outside the grid, no radiation flux
           endif
        end do !j
     end do !i

     !I0 = I0 / npix2
     !I0c = I0c / npix2
     !below now

     if (iter < n_iter_min) then
        ! On itere par defaut
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! On arrete pour pas tourner dans le vide
         !write(*,*) "Warning : converging pb in ray-tracing"
         !write(*,*) " Pixel", ipix, jpix
        if (diff > 1) write(*,*) 'pixel not converged:', iter, subpixels, diff
        exit infinie
     else
        ! On fait le test sur a difference
        diff = maxval( abs(I0 - Iold) / (I0 + 1e-300_dp) )
        ! There is no iteration for Q, U, V, assuming that if I is converged, then Q, U, V also.
        ! Can be added and then use diff as max(diff, diffQ, diffU, diffV)
        !write(*,*) 'iter pixel:', ipix, jpix, i, j, iter, subpixels, diff
        if (diff > precision ) then
           ! On est pas converge
           subpixels = subpixels * 2
        else
           !write(*,*) "Pixel converged", ipix, jpix, i, j, iter, diff
           exit infinie
        end if
     end if ! iter
     iter = iter + 1
   end do infinie

  !Prise en compte de la surface du pixel (en sr)

  !Because unit(Bnu)=unit(I) = W/m2/Hz/sr = unit(Blambda*lambda**2/c)
  !if we want I in W/m2 --> I*nu
  !nu = 1d0 !c_light / NLTEspec%lambda * 1d9 !to get W/m2 instead of W/m2/Hz !in Hz
  ! --> Deprecated, convert in post processing or before adding disk, dust, molecular flux to
  ! atomic flux or reversed case.

  ! Flux out of a pixel in W/m2/Hz
  normF = (pixelsize / (distance*pc_to_AU) )**2  / npix2 !divide by number of subpixels
  !normF = (pixelsize / R0)**2  / npix2
!   I0 = nu * I0 * (pixelsize / (distance*pc_to_AU) )**2
!   I0c = nu * I0c * (pixelsize / (distance*pc_to_AU) )**2

  ! adding to the total flux map.
  if (RT_line_method==1) then
    NLTEspec%Flux(:,1,1,ibin,iaz) = NLTEspec%Flux(:,1,1,ibin,iaz) + I0 * normF
    NLTEspec%Fluxc(:,1,1,ibin,iaz) = NLTEspec%Fluxc(:,1,1,ibin,iaz) + I0c * normF
    if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
     NLTEspec%F_QUV(1,:,1,1,ibin,iaz) = NLTEspec%F_QUV(1,:,1,1,ibin,iaz)+QUV(1,:) * normF !U
     NLTEspec%F_QUV(2,:,1,1,ibin,iaz) = NLTEspec%F_QUV(2,:,1,1,ibin,iaz)+QUV(2,:) * normF !Q
     NLTEspec%F_QUV(3,:,1,1,ibin,iaz) = NLTEspec%F_QUV(3,:,1,1,ibin,iaz)+QUV(3,:) * normF !V
    end if
  else
    NLTEspec%Flux(:,ipix,jpix,ibin,iaz) = I0 * normF
    NLTEspec%Fluxc(:,ipix,jpix,ibin,iaz) = I0c * normF
    if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
     NLTEspec%F_QUV(1,:,ipix,jpix,ibin,iaz) = QUV(1,:) * normF !U
     NLTEspec%F_QUV(2,:,ipix,jpix,ibin,iaz) = QUV(2,:) * normF !Q
     NLTEspec%F_QUV(3,:,ipix,jpix,ibin,iaz) = QUV(3,:) * normF !V
    end if
  end if
  if (lcontrib_function) NLTEspec%Ksi(:,:,ibin,iaz) = &
                         NLTEspec%Ksi(:,:,ibin,iaz) + S_contrib2(:,:) * (pixelsize*au_to_m)**2 / npix2 !In W/Hz


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
  real(kind=dp) :: x0,y0,z0,l,u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcorner
  real(kind=dp):: taille_pix, nu
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36
  integer, parameter :: n_ray_star = 1000
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp):: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda, lM, lR, lB, ll
  logical :: lresolved

  write(*,*) "Vector to observer =", real(tab_u_rt(ibin,iaz)),real(tab_v_rt(ibin,iaz)),real(tab_w_rt(ibin))
  write(*,*) "i=", real(tab_RT_incl(ibin)), "az=", real(tab_RT_az(iaz))

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
     n_iter_min = 1 !1 !3
     n_iter_max = 1 !1 !6
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

 RETURN
 END SUBROUTINE EMISSION_LINE_MAP

 SUBROUTINE Atomic_transfer()
 ! --------------------------------------------------------------------------- !
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
 ! --------------------------------------------------------------------------- !
  integer :: atomunit = 1, nact
  integer :: icell
  integer :: ibin, iaz
  integer, parameter :: Nrayone = 1
  real(kind=dp) :: dumm
  character(len=20) :: ne_start_sol = "H_IONISATION"
  character(len=20)  :: newPRT_SOLUTION = "FULL_STOKES"

!
!   integer :: i, j
!   real(kind=dp), dimension(:), allocatable :: test_col
!   real(kind=dp), dimension(:,:), allocatable :: test_col2, Cji
!     type (Atomicline) :: line, line1
!     type (AtomicContinuum) :: cont
!  real(kind=dp) :: M(3,4), B(12), C(12), M2(10,1000), B2(10000)
!   M(1,1) = 11; M(1,3) = 13; M(1,2) = 12; M(1,4) = 14
!   M(2,1) = 21; M(2,2) = 22; M(2,3) = 23; M(2,4) = 24
!   M(3,1) = 31; M(3,2) = 32; M(3,3) = 33; M(3,4) = 34
!   !M(Nlevel, Ncells)
!   B = flatten(3,4, M)
!   C = flatten2(3,4, M)
! !  M2(:,:) = 0d0
! !  B2 = flatten(10,1000,M2)
! !  B2 = flatten2(10,1000,M2)
!   write(*,*)"B=", B ! M(1:), M(2,:).. = M11, M12, M13 ..M21, M22 ..
! write(*,*) "C=", C  !M(:,1) M(:,2) ... = M11, 21, M31 .. M12, M22, M32
! write(*,*) reform2(3,4, C) - M
! stop
 ! -------------------------------INITIALIZE AL-RT ------------------------------------ !
  Profile => IProfile
  INTEG_RAY_LINE => INTEG_RAY_LINE_I
  !only one available yet, I need one unpolarised, faster and more accurate.
  Voigt => VoigtHumlicek

  !---> Futur deprecation
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
  if (.not.lpluto_file) then
   !CALL spherical_shells_model()
   !CALL spherical_star()
   !CALL magneto_accretion_model()
   if (lmodel_ascii) then
    CALL readAtmos_ascii(density_file)
    CALL writeHydrogenDensity()
    CALL writeTemperature()
    !could write also T, V ect and grid domain icompute_atomRT, but present in the model.s
    !I write the density here just to check the cell mapping.
   end if
  end if
!! --------------------------------------------------------- !!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------- READATOM and INITIZALIZE POPS ------------------------------ !

  !Read atomic models and allocate space for n, nstar
  ! on the whole grid space.
  CALL readAtomicModels(atomunit)
  if (atmos%NactiveAtoms > 0) atmos%Nrays = 200

 ! ------------------------------------------------------------------------------------ !
!   test impacts
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
 ! ------------------------------------------------------------------------------------ !

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

  CALL setLTEcoefficients () !write pops at the end because we probably have NLTE pops also

 ! ------------------------------------------------------------------------------------ !
 ! ------------- INITIALIZE WAVELNGTH GRID AND BACKGROUND OPAC ------------------------ !
 ! ------------------------------------------------------------------------------------ !
  atmos%coherent_scattering=lcoherent_scattering
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
  CALL reallocate_mcfost_vars() !assumes more than 1 wavelength otherwise delta_wl is 0!
  ! --- END ALLOCATING SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT --!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- NLTE LOOP -------------------------------------- !
  !The BIG PART IS HERE
if (atmos%Nactiveatoms > 0) then
write(*,*) " "
write(*,*) " "

write(*,*) " !!!!!!OPENMP DEACTIVATED FOR NLTE LOOP ATM!!!!!!!!!"

write(*,*) " "
write(*,*) " "
  !if (atmos%Nactiveatoms > 0) &
     CALL NLTEloop()
endif
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------- WRITE CONVERGED POPULATIONS ------------------------------ !
 ! ------------------------------------------------------------------------------------ !
  do icell=1,atmos%NactiveAtoms
   !write final atom%n (count==0) + lte pops (ilte_only = .false.)
   CALL writePops(atmos%Atoms(icell)%ptr_atom, 0)
  end do
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- MAKE IMAGES ------------------------------------ !

   !should add an other flag to force to avoid continua lines
   !Define a wavelength grid for image with only lines, if not input wavelength
!    if (.not.(ltab_wavelength_image)) then !or for all passive lines also
!     ltab_wavelength_image = .true.
!     !!CALL write_wavelengths_table_NLTE_lines(NLTEspec%lambda) !!after grid creation actually
!     tab_wavelength_image = "line_waves.s"
!    endif

  if (ltab_wavelength_image) then
   !Check smarter to deallocate/reallocated NLTE wavelength arrays
   CALL initSpectrumImage() !deallocate waves arrays/ define a new grid
   !shorter than the grid for NLTE / reallocate waves arrays
   !write(*,*) maxval(Hydrogen%continua(1)%alpha), maxval(atmos%Atoms(1)%ptr_atom%continua(1)%alpha)
   !write(*,*) loc(Hydrogen)==loc(Atmos%Atoms(1)%ptr_atom)

   if (lstore_opac) CALL storeBackground() !recompute background opac
   ! TO DO: add NLTE continua and LTE/NLTE lines if possible
   !!CALL reallocate_mcfost_vars() !wavelength only?
   CALL reallocate_mcfost_wavelength_arrays()
  end if !grid for image
!check
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
    allocate(QUV(3, NLTEspec%Nwaves))
    QUV = 0d0
  end if

  CALL alloc_flux_image()
  write(*,*) "Computing emission flux map..."

  !Use converged NLTEOpac
  if (lcontrib_function) then
   write(*,*) "   -> including contribution functions."
   allocate(S_contrib(NLTEspec%Nwaves,atmos%Nspace,NLTEspec%NPROC))
   allocate(chil(NLTEspec%Nwaves), Sl(NLTEspec%Nwaves),S_contrib2(NLTEspec%Nwaves,n_cells))
   S_contrib(:,:,:) = 0d0; S_contrib2(:,:) = 0d0
   INTEG_RAY_LINE => NULL()
   INTEG_RAY_LINE => INTEG_RAY_LINE_I_CNTRB !same as INTEG_RAY_LINE_I, without nlte (but %n is known)
  end if

  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
       CALL EMISSION_LINE_MAP(ibin,iaz)
     end do
  end do
  CALL WRITE_FLUX()
  if (allocated(QUV)) deallocate(QUV)
  if (lcontrib_function) then
   CALL WRITE_CNTRB_FUNC_pix()
   deallocate(S_contrib, chil, Sl, S_contrib2)
  end if
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------- CLEANING ------------------------------------------ !
 !close file after NLTE loop
!Temporary: because we kept in memory, so file is closed earlier
!  do nact=1,atmos%Nactiveatoms
!   CALL closeCollisionFile(atmos%ActiveAtoms(nact)%ptr_atom) !if opened
!  end do
 !CALL WRITEATOM() !keep C in memory for that ?
 if (atmos%coherent_scattering) CALL freeJ()
 CALL freeSpectrum() !deallocate spectral variables
 CALL free_atomic_atmos()
 NULLIFY(optical_length_tot, Profile, Voigt, INTEG_RAY_LINE)

 RETURN
 END SUBROUTINE
 ! ------------------------------------------------------------------------------------ !

 SUBROUTINE NLTEloop() !for all active atoms
 ! -------------------------------------------------------- !
  ! Descriptor here
 ! -------------------------------------------------------- !
#include "sprng_f.h"

  integer, parameter :: n_rayons_start = 50
  integer, parameter :: n_rayons_start2 = 10
  integer :: n_rayons_max != n_rayons_start2 * (2**(n_iter2_max-1)
  !! it is atmos%Nrays in atom transfer
  !!                   True absolute change. fact_etape not used anymore
  !!real, parameter :: precision_sub = 1e-3 !1e-4
  !!real, parameter :: precision = 1.0e-4 !1e-4
  integer :: etape, etape_start, etape_end, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  integer, dimension(nb_proc) :: max_n_iter_loc

  integer :: la, imu, ncells_filled

  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged

  real :: rand, rand2, rand3!, fac_etape

  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dN, dN1, dJ
  real(kind=dp), allocatable :: dM(:), Jold(:,:), taub(:,:,:)
                                   !futur global flag
  logical :: labs, disable_subit, iterate_ne = .false., accelerated, ng_rest
  integer :: atomunit = 1, nact, maxIter
  integer :: icell, iorder, i0_rest, n_iter_accel
  integer :: Nlevel_total = 0, NmaxLevel, ilevel, max_sub_iter
  character(len=20) :: ne_start_sol = "H_IONISATION"!"NE_MODEL"
  real(kind=dp), dimension(3, atmos%Nrays, nb_proc) :: xyz0, uvw0
  !real(kind=dp), dimension(:), allocatable :: xmu, wmu
  !integer :: to_obs0
  type (AtomType), pointer :: atom

  write(*,*) "   -> Solving for kinetic equations for ", atmos%Nactiveatoms, " atoms"

  if (allocated(ds)) deallocate(ds)
  allocate(ds(atmos%Nrays,NLTEspec%NPROC))
  ds = 0d0 !meters

  if (atmos%include_xcoupling) then
   allocate(taub(NLTEspec%Nwaves, atmos%Nrays, NLTEspec%NPROC))
   taub(:,:,:) = 0d0
  endif

  !move to initSol
  allocate(dM(atmos%Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
  if (atmos%coherent_scattering) then
   allocate(Jold(NLTEspec%Nwaves, atmos%Nspace))!move to initsol
   Jold = 0d0
  endif

  n_rayons_max = atmos%Nrays
  xyz0(:,:,:) = 0d0
  uvw0(:,:,:) = 0d0

  labs = .true. !to have ds at cell icell
  id = 1
  etape_start = 2
  etape_end = 2
  if (etape_start==0) then
   write(*,*) "Warn: etape 0 not accurate"
  else if (etape_start==1 .or. etape_end==1) then
   write(*,*) " Etape 1 with new angle quad not implemented, setting to 2"
   if (etape_start==1) etape_start=2
   if (etape_end==1) etape_end=2
  endif

  disable_subit = atmos%include_xcoupling!set to true to avoid subiterations over the emissivity

  if (lNg_acceleration) then
   n_iter_accel = 0
   i0_rest = 0
   ng_rest = .false.
  endif

  iterate_ne = (n_iterate_ne>0)

  if (lforce_lte) disable_subit = .true.

  max_sub_iter = 10000 !to continue until convergence, reduce to force a number, or increase for infinity
  maxIter = 100
 ! ----------------------------  INITIAL POPS------------------------------------------ !
   CALL Init_NLTE(sub_iterations_enabled=.not.disable_subit)
 !  ------------------------------------------------------------------------------------ !

! ncells_filled = 0
! do icell=1, atmos%Nspace
!  if (atmos%icompute_atomRT(icell) >0) ncells_filled = ncells_filled + 1
! end do
! open(16,file="testI3", status='old')
! write(16, *) etape_end-etape_start+1, ncells_filled, NLTEspec%Nwaves, n_rayons_start

     do etape=etape_start, etape_end

      !precision = fac_etape* precision, not anymore
      if (etape==0) then !two rays, not accurate
        lfixed_rays=.true.
        n_rayons = 2
        iray_start = 1
        !fac_etape = 1e-1
        lprevious_converged = .false.

      !building
      else if (etape==1) then ! Try new angular quadrature
      !N random positions + N fixed directions for each
      !adding to obs with negative direction vector
        !!CALL Gauleg(0d0, 1d0, xmu, wmu,n_rayons_start)
        !!
        !!to_obs0 = -1
        lfixed_rays = .true.
        n_rayons = min(n_rayons_max,n_rayons_start)
  		iray_start = 1
  		lprevious_converged = .false.

      else if (etape==2) then ! a random direction associated with a random position
  		lfixed_rays = .true.
  		n_rayons = min(n_rayons_max,n_rayons_start2)
  		iray_start = 1
  		lprevious_converged = .false.

  	  else
  	    CALL ERROR("etape unkown")
  	  end if

  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0

        do while (.not.lconverged)

        	n_iter = n_iter + 1
        	if (n_iter > maxIter) exit !change step
            write(*,*) " -> Iteration #", n_iter, " Step #", etape

  			if (lfixed_rays) then !beware if lfixed_rays = .false. in etape 2
  								  ! A.T.M fixed rays is assumed
    			stream = 0.0
    			do i=1,nb_proc
     				stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    			end do
 			 end if

            max_n_iter_loc = 0

 			!!$omp parallel &
            !!$omp default(none) &
            !!$omp private(icell, id, atom,nact) &
            !!$omp shared(gpop_old,atmos)
            !!$omp do schedule(static,1)
            do icell=1, atmos%Nspace
   				!!$ id = omp_get_thread_num() + 1
   				if (atmos%icompute_atomRT(icell)>0) then
            	do nact=1, atmos%NactiveAtoms
            	    atom => atmos%ActiveAtoms(nact)%ptr_atom
             		gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
             		atom => NULL()
            	end do
            	end if
            end do
        	!!$omp end do
        	!!$omp end parallel


 			!!$omp parallel &
            !!$omp default(none) &
            !!$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, imu)&
            !!$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, nact, atom) &
            !!$omp shared(atmos,NLTEspec, taub,dpops_sub_max_err) &
            !!$omp shared(xyz0, uvw0, lkeplerian,n_iter) & !before nact was shared
            !!$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,max_sub_iter) &
            !!$omp shared(n_cells, pop_old, ds,disable_subit, dN, dN1,gpop_old,lforce_lte) & !pop
            !!$omp shared(lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape)
            !!$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !!$ id = omp_get_thread_num() + 1

   				!Should consider to add a different positive flag for cells in NLTE
   				!Because, if we keep low T, low nH and low ne region, the transfer
   				!could become very unstable, no ?
   				if (atmos%icompute_atomRT(icell)>0) then !nor transparent nor dark

					!comute directions, store them and compute line weights = ray integral
					!then do the ray propagation and build gamma ray by ray
  			        CALL initGamma(id,icell)

  			        !compute directions + postitions of rays for angular integration
  			        !set up also the norm of the integral over continua and lines
  			        !it has a loop on the rays, which is fast and allows to built the Rate matrix rays by ray and avoid some loop on rays
                    CALL compute_directions(etape, id, icell, iray_start, n_rayons, stream(id), xyz0(:,:,id), uvw0(:,:,id))

					do iray=iray_start, iray_start-1+n_rayons
					    !!$omp critical
						CALL INTEG_RAY_LINE(id, icell, xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
											uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, labs)
						!!$omp end critical
						if (atmos%include_xcoupling) then
							CALL integrate_tau_bound(id,icell,xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
								uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, 1, taub)

						    NLTEspec%Psi(:,iray,id) = NLTEspec%Psi(:,iray,id) * dexp(-taub(:,iray,id))
						endif
                        CALL FillGamma_mu(id, icell, iray, n_rayons, lforce_lte) !ray by ray since wphi is known
      			    enddo !iray

      				if (atmos%coherent_scattering) CALL calc_J_coherent(id, icell, n_rayons)


! write(16, *) etape, n_iter, icell
! do la=1,NLTEspec%Nwaves
!    write(16,'(52E)') NLTEspec%lambda(la),(NLTEspec%I(la, imu, id), imu=1,n_rayons), sum(NLTEspec%I(la, :, id))/n_rayons
!    !make sure J is the same
!    !write(*,*) sum(NLTEspec%I(la, :, id))/n_rayons, NLTEspec%J(la, icell)
! end do

    		     	n_iter_loc = 0
    				if (disable_subit) then
    				 CALL updatePopulations(id, icell)
    				 lconverged_loc = .true.
    				else
    				 lconverged_loc = .false.
    				end if


     				!!sub iteration on the local emissivity, keeping Idag fixed
     				!!only iterate on cells which are not converged
     				do while (.not.lconverged_loc)
     				!write(*,*) "Starting subit for cell ", icell
       					n_iter_loc = n_iter_loc + 1

            			do nact=1, atmos%NactiveAtoms
            	    		atom => atmos%ActiveAtoms(nact)%ptr_atom
             				pop_old(nact, 1:atom%Nlevel,id) = atom%n(:,icell)
             				atom => NULL()
            			end do

						!Solve SEE for all atoms
						CALL updatePopulations(id, icell)

						diff = 0.0 !keeps track of the maximum dpops(dN) among each atom.
								   !for this cell
     					do nact=1,atmos%NactiveAtoms
     					    atom => atmos%ActiveAtoms(nact)%ptr_atom
     					    dN1 = 0.0 !for one atom
     						do ilevel=1,atom%Nlevel
     				    		dN = dabs(1d0 - pop_old(nact,ilevel,id)/(tiny_dp + atom%n(ilevel,icell)))
     				    		diff = max(diff, dN)
     				    		dN1 = max(dN1, dN)
     						end do
     						!if (dN1 >= 1) &
     						!	write(*,*) id, " --> subit",n_iter_loc,atom%ID, " dpops = ", dN1
     						atom => NULL()
     					end do

       					if (diff < dpops_sub_max_error) then!precision_sub
       						lconverged_loc = .true.
     					    !write(*,*) id, n_iter_loc, "dpops(sub) = ", diff
       					else
        					!recompute opacity of this cell., but I need angles and pos...
       						!NLTEspec%I not changed
  			        		CALL initGamma(id,icell)
 							do iray=iray_start, iray_start-1+n_rayons
      							!I unchanged
							    CALL init_local_field_atom(id, icell, iray, &
							         xyz0(1,iray,id), xyz0(2,iray,id), xyz0(3,iray,id), &
							         uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id))
									 !if (atmos%include_xcoupling) then
										!CALL integrate_tau_bound(id,icell,xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
								        !      uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, 1, taub)

						   			    !NLTEspec%Psi(:,iray,id) = NLTEspec%Psi(:,iray,id) * dexp(-taub(:,iray,id))
									!endif
							    CALL FillGamma_mu(id, icell, iray, n_rayons, lforce_lte)
      						enddo !iray

       					end if
       					if (n_iter_loc >= max_sub_iter) then
       					  if (diff>1) write(*,*) id, " sub-it not converged after", n_iter_loc, &
       					  	" iterations; diff=", diff
       					  lconverged_loc = .true.
       					  !write(*,*) id, n_iter_loc, "dpops(sub) = ", diff

       					end if
     				end do !local sub iteration
     	            if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
     			end if !icompute_atomRT
     		end do !icell
        	!!$omp end do
        	!!$omp end parallel

        	!Global convergence criterion
        	!I cannot iterate on unconverged cells because the radiation coming for each cell
        	!depends on the other cell.
        	!can be parallel
        	diff = 0d0
        	dM(:) = 0d0
        	dJ = 0.0

     		!not parallel yet
     		accelerated=.false.
     		if (lNg_acceleration .and. (n_iter > iNg_Ndelay)) then
     		  iorder = n_iter - iNg_Ndelay !local number of iterations accumulated
     		  if (ng_rest) then
     		    write(*,*) "    --> Acceleration relaxes...", iorder-i0_rest
     		    if (iorder - i0_rest == iNg_Nperiod) ng_rest = .false.
     		  else
     		    i0_rest = iorder
                do nact=1,atmos%NactiveAtoms
     			 atom => atmos%ActiveAtoms(nact)%ptr_atom
     		     allocate(flatpops(atom%Ngs%N))
     		     !flatpops=flatten(atom%Nlevel,atmos%Nspace,atom%n)
     		     flatpops=flatten2(atom%Nlevel,atmos%Nspace,atom%n)
     		     accelerated = Acceleration(atom%Ngs, flatpops)
     		     if (accelerated) then
     		      !atom%n(:,:) = reform(atom%Nlevel, atmos%Nspace, flatpops)
     		      atom%n(:,:) = reform2(atom%Nlevel, atmos%Nspace, flatpops)
                  n_iter_accel = n_iter_accel + 1 !True number of accelerated iter
                  write(*,*) "    ++> ", atom%id, "accelerated iteration #", n_iter_accel
                  ng_rest = .true.
     		     endif
     		     deallocate(flatpops)

     		     atom => NULL()
     		    enddo
             endif
     		endif

  			cell_loop2 : do icell=1, atmos%Nspace
  				if (atmos%icompute_atomRT(icell)>0) then
						dN = 0d0 !for all levels of all atoms of this cell
     					do nact=1,atmos%NactiveAtoms
     					    atom => atmos%ActiveAtoms(nact)%ptr_atom
     						!do ilevel=1,atom%Nlevel
     						ilevel = 1


     				    		dN1 = dabs(1d0-gpop_old(nact,ilevel,icell)/(atom%n(ilevel,icell)+tiny_dp))
     				    		dN = max(dN1, dN) !compare with all atoms and levels
     				    						  !for this cell

!      						    if (dN >= 1) then
!      				    		   write(*,*) atom%ID, "Pops:",icell, ilevel," dN =",dN
!      				    		   write(*,*) "nold(i) =",gpop_old(nact,ilevel,icell), "n(i) =",atom%n(ilevel,icell)
!      				    		 end if

     				    		 dM(nact) = max(dM(nact), dN1) !compare for one atom
     				    		 							   !for all cells
     						!end do
     						atom => NULL()
     					end do
     					diff = max(diff, dN) !compare for all atoms and all cells
     					if (atmos%coherent_scattering) &
     						dJ = max(dJ,dabs(1d0 - maxval(Jold(:,icell)/(tiny_dp + NLTEspec%J(:,icell)))))
     			end if
     		end do cell_loop2
     		if (dJ /= 0) write(*,*) " dJ = ", dJ
     		if (atmos%coherent_scattering) Jold(:,:) = NLTEspec%J(:,:)


         	!if (maxval(max_n_iter_loc)> max_sub_iter) &
         	if (.not.disable_subit)	write(*,*) maxval(max_n_iter_loc), "sub-iterations"
         	do nact=1,atmos%NactiveAtoms
         	 if (accelerated) then
         	  write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dM = ", dM(nact)," (Accelerated)"
         	 else
         	  write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dM = ", dM(nact)
         	 endif
         	enddo
         	write(*,*) "dpops =", diff !at the end of the loop over n_cells
        	!Use dJ as a criterion for convergence ?
        	!if (atmos%coherent_scattering) diff = dJ

        	lconverged = (real(diff) < dpops_max_error)

!         	if (real(diff) < dpops_max_err) then !precision
!            		if (lprevious_converged) then
!             	  lconverged = .true.
!            		else
!             	  lprevious_converged = .true.
!           	    endif
!         	else
!            		lprevious_converged = .false.
!            		if (.not.lfixed_rays) then
!               		n_rayons = n_rayons * 2
!               		write(*,*) ' -- Increasing number of rays'
!              		if (n_rayons > n_rayons_max) then
!               			if (n_iter >= maxIter) then
!              		 		write(*,*) "Warning : not enough rays to converge !!"
!                  			lconverged = .true.
!               			end if
!               		end if
!
!               ! On continue en calculant 2 fois plus de rayons
!               ! On les ajoute a l'ensemble de ceux calcules precedemment
! !              iray_start = iray_start + n_rayons
!
!           	   end if
!         	end if
        									!if ==1 for Ne_period=3: 1 ok, 2, 3, 4 ok, 5, 6, 7 ok etc
        									! if 0: 1, 2, 3ok, 4, 5, 6ok ect
        	!Only if specified
        	!if (disable_subit) then
        		if (iterate_ne .and. (mod(n_iter,n_iterate_ne)==0))  then
        		write(*,*) " Solve ne Global iteration:", n_iter
        	 	write(*,*) " --> old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
        	 	CALL SolveElectronDensity(ne_start_sol)
        	 !Recompute LTE pops used in continua radiative rates
        	 !for Activeatoms only ?
        	 	do nact=1,atmos%NactiveAtoms !if over Active only, should recompute for passive also
        	 							  !but this preserve the constant background opac
        	 							  !as passiveatoms%n = passiveatoms%nstar with ne init.
               	atom => atmos%ActiveAtoms(nact)%ptr_atom
             !do nact=1,atmos%NAtom
               !atom => atmos%Atoms(nact)%ptr_atom
               	CALL LTEpops(atom,.true.)
               	atom => Null()
             	end do
        		end if
        	!endif

!         	do nact=1, atmos%NActiveAtoms
!         	 CALL writePops(atmos%Activeatoms(nact)%ptr_atom, 1)
!         	enddo
	    end do !while
        write(*,*) "Threshold =", dpops_max_error
	  end do !over etapes

!close(16)

  !Force to compute a new value of electron density after convergence
  !and new lte populations for all atoms
  if (n_iterate_ne < 0) then
   write(*,*) "END LOOP: old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
   CALL SolveElectronDensity(ne_start_sol)
   !Recompute for all atoms the LTE pops
   write(*,*) "   recompute LTE pops for all atoms"
   do nact=1,atmos%NAtom
      atom => atmos%Atoms(nact)%ptr_atom
      CALL LTEpops(atom,.true.)
      atom => Null()
   end do
  else if (iterate_ne) then
   write(*,*) "END LOOP: recompute LTE pops of passive atoms"
   !Recompute for all passive atoms the LTE pops
   !because the lte pops of active atoms only is updated with iterate_ne
   do nact=1,atmos%NpassiveAtoms
      atom => atmos%Atoms(nact)%ptr_atom
      CALL LTEpops(atom,.true.)
      atom => Null()
   end do
  end if

  if (iterate_ne .or. (n_iterate_ne < 0)) &
  	CALL writeElectron(.true.) !the .true. means append, to compare with initial solution.
 ! -------------------------------- CLEANING ------------------------------------------ !
  !to move inside free_nlte_sol
  deallocate(dM)
  if (allocated(Jold)) deallocate(Jold)
  CALL free_nlte_sol(disable_subit)
  deallocate(ds)
  if (atmos%include_xcoupling) deallocate(taub)

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

     !if (associated(Voigt)) Voigt => NULL
     !Voigt => VoigtHumlicek !Already the algortihm used

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
     !if (associated(Voigt)) Voigt => NULL
     !Voigt => VoigtHumlicek !Already the algortihm used
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
  tab_delta_lambda(1) = 0d0! + tiny_dp !! I assume that Nwaves > 1 here
  !!if (NLTEspec%Nwaves==1) tab_delta_lambda(1) = tab_delta_lambda(1) + tiny_dp

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
  !!if I do not put a tiny_dp here, if only one wavelength delta_wl is 0
  !which creates an overflow in stars.f90, because spectre and spectre0 are 0
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

 SUBROUTINE compute_directions(etape, id, icell, rstart, n_rayons, streami, xyz, uvw)
#include "sprng_f.h"

 	integer, intent(in) :: etape, id, icell, rstart, n_rayons
 	SPRNG_POINTER, intent(in) :: streami
 	real(kind=dp), dimension(3,n_rayons), intent(out) :: xyz, uvw
 	integer :: iray
 	real(kind=dp) :: x0, y0, z0, u0, w0, v0, norme
 	real :: rand, rand2, rand3, W02, ARGMT, srW02

    do iray=rstart, rstart-1+n_rayons

    	if (etape==0) then
        	! Position = milieu de la cellule
        	x0 = r_grid(icell)
        	y0 = 0.0_dp
        	z0 = z_grid(icell)
        	if (lkeplerian) then
                       ! Direction verticale "z"
        		if (iray==1) then
            		w0=1.0_dp !nz
       			else
            		w0=-1.0_dp
        		endif
            	u0 = 0.0_dp !nx
            	v0 = 0.0_dp !ny
        	else
           		norme = sqrt(x0*x0 + y0*y0 + z0*z0)
       			if (iray==1) then
            		u0 = x0/norme !sin(theta)sin(phi) = nx
            		v0 = y0/norme !ny
            		w0 = z0/norme !nz = mu = cos(theta)
        		else
            		u0 = -x0/norme !-1  backward
            		v0 = -y0/norme
           			w0 = -z0/norme
           		endif
        	endif !lkeplerian

         else if (etape==1) then
            write(*,*) "step 1 not implemented yet"
            stop
            !random position + fixed directions

            !w0 = xmu(iray) * real(to_obs)
            !u0 = dsqrt(1.-xmu(iray)**2) * cos(pmu(iray)) * to_obs
            !v0 = dsqrt(1.-xmu(iray)**2) * sin(pmu(iray)) * to_obs

      	 else !etape 2
                   	    ! Position aleatoire dans la cellule
         	rand  = sprng(streami)
            rand2 = sprng(streami)
            rand3 = sprng(streami)

            CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
            rand = sprng(streami)
            W0 = 2.0_dp * rand - 1.0_dp !nz
            W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
            SRW02 = sqrt(W02)
            rand = sprng(streami)
            ARGMT = PI * (2.0_dp * rand - 1.0_dp)
            U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
            V0 = SRW02 * sin(ARGMT) !ny = sin(theta) * sin(phi)
		 end if !etape

		 xyz(1,iray) = x0; xyz(2,iray) = y0; xyz(3,iray) = z0
		 uvw(1,iray) = U0; uvw(2,iray) = V0; uvw(3,iray) = W0
		 CALL compute_integral_weight(id, icell, iray, n_rayons, x0, y0, z0, u0, v0, w0)
 	enddo
 RETURN
 END SUBROUTINE compute_directions

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
     gamma(:) = (dexp(hc_k/NLTEspec%lambda(:)/real(etoile(i_star)%T,kind=dp))-1)/&
     			(dexp(hc_k/NLTEspec%lambda(:)/etoile(i_star)%SurfB(ns)%T)-1)
     !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
   end if

   !3) Apply Limb darkening
   if (llimb_darkening) then
     !!LimbDarkening = Interp1D(real(mu_limb_darkening,kind=dp),real(limb_darkening,kind=dp),mu)
     !LimbDarkening = interp_dp(limb_darkening, mu_limb_darkening, mu)
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

 SUBROUTINE INTEG_RAY_LINE_I_CNTRB(id,icell_in,x,y,z,u,v,w,iray,labs)
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
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs !used in NLTE but why?
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, chir
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, Istar
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chi_I
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star!, idref
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp !go from surface down to the star
  tau_c = 0.0_dp


  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0
  S_contrib(:,:,id) = 0d0 !global


  ! -------------------------------------------------------------- !
  !*** propagation dans la grille ***!
  ! -------------------------------------------------------------- !
  ! Will the ray intersect a star
  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
    ! Indice de la cellule
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif

    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN


    if (lintersect_stars) then
      if (icell == icell_star) then
       Istar(:) = 0d0
       CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,Istar)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + Istar*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + Istar*dexp(-tau_c)
       RETURN
      end if
    endif

    nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
    previous_cell = 0 ! unused, just for Voronoi

    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)

    !count opacity only if the cell is filled, else go to next cell
    if (lcellule_non_vide) then

     lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m !l_contrib in AU
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m

    ! eval_operator = (labs .and. (nbr_cell == 1))

     CALL initAtomOpac(id)
     !if (eval_operator) CALL init_psi_operator(id, iray)

     CALL NLTEopacity(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
     !never enter NLTEopacity if no activeatoms
     if (lstore_opac) then !not updated during NLTE loop, just recomputed using initial pops
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chi_I(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1)

      dtau(:)   = l_contrib * chi_I(:)
      dtau_c(:) = l_contrib * (NLTEspec%AtomOpac%Kc(icell,:,1) + NLTEspec%AtomOpac%chic_nlte(:,id))
      Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
               NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / &
                 (NLTEspec%AtomOpac%Kc(icell,:,1) + &
                 NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%chi(:,id) + 1d-300)

      Snu_c = (NLTEspec%AtomOpac%jc(icell,:) + NLTEspec%AtomOpac%etac_nlte(:,id)) / &
            (NLTEspec%AtomOpac%Kc(icell,:,1) + 1d-300 + NLTEspec%AtomOpac%chic_nlte(:,id))


      chil(:) = NLTEspec%AtomOpac%chi(:,id) -  NLTEspec%AtomOpac%chic_nlte(:,id) + &
        			   NLTEspec%AtomOpac%chi_p(:,id) + 1d-300
      Sl(:) = NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id) - &
                 NLTEspec%AtomOpac%etac_nlte(:,id)

      !chir = NLTEspec%AtomOpac%Kc(icell,idref,1) + NLTEspec%AtomOpac%chic_nlte(idref,id)

     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chi_I(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)
      ! Epaisseur optique
      ! chi_p contains both thermal and continuum scattering extinction
      dtau(:)   =  l_contrib * chi_I(:)
      dtau_c(:) = l_contrib * (NLTEspec%AtomOpac%chi_c(:,id) + NLTEspec%AtomOpac%chic_nlte(:,id))

      Snu = (NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / &
                 (NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%chi(:,id) + 1d-300)

      ! continuum source function
      Snu_c = (NLTEspec%AtomOpac%eta_c(:,id) + NLTEspec%AtomOpac%etac_nlte(:,id)) / &
      			 (NLTEspec%AtomOpac%chi_c(:,id) + 1d-300 + NLTEspec%AtomOpac%chic_nlte(:,id))

      chil(:) = NLTEspec%AtomOpac%chi_p(:,id)-NLTEspec%AtomOpac%chi_c(:,id) + &
        			NLTEspec%AtomOpac%chi(:,id) -  NLTEspec%AtomOpac%chic_nlte(:,id) + 1d-300
        !actually eta here
      Sl(:) = NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id) - &
                 NLTEspec%AtomOpac%eta_c(:,id) - NLTEspec%AtomOpac%etac_nlte(:,id) !) / chil(:)

       !cont only
      !chir = NLTEspec%AtomOpac%chi_c(idref,id) + NLTEspec%AtomOpac%chic_nlte(idref,id)
      !chir = chi_I(idref)
    end if


    NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
    NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c


     !!
     !! In general I m not sure of what is the attenuation coefficient here
     !! tau here = tau at the edge of the cell, before entering
     !! tau after = tau(k+1) = tau(k) + dtau(k) == tau at the edge of the cell before leaving
     !! which one should be used, as CF, is the integrand of the RTE, thus tau (or tau+dtau).


     S_contrib(:,icell,id) = (chil(:) * NLTEspec%Ic(:,iray,id) - Sl(:)) * dexp(-tau-dtau) / (tiny_dp + chi_I)

     tau = tau + dtau
     tau_c = tau_c + dtau_c

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_I_CNTRB

 SUBROUTINE integrate_tau_bound(id,icell_in,x,y,z,u,v,w,iray,to_obs, tau)
 ! ------------------------------------------------------------------------------- !
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray, to_obs
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(NLTEspec%Nwaves) :: dtau, chi_I
  real(kind=dp), intent(out), dimension(NLTEspec%NWaves,atmos%Nrays, NLTEspec%Nproc) :: tau
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau(:,iray,id) = 0d0
  !!write(*,*) " vector direction times ", to_obs, u,v,w

  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
    ! Indice de la cellule
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif

    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN


    if (lintersect_stars) then
      if (icell == icell_star) RETURN
    endif

    nbr_cell = nbr_cell + 1

    previous_cell = 0

    call cross_cell(x0,y0,z0, u,v,w, icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
!     if (l_contrib /= l) then
!       write(*,*) "l_contrib, l"
! 	  write(*,*) l_contrib/ etoile(1)%r, l / etoile(1)%r
! 	endif
    if (lcellule_non_vide) then


     l_contrib = l_contrib * AU_to_m !l_contrib in AU

     CALL initAtomOpac(id)

     CALL NLTEopacity(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)
     if (lstore_opac) then
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chi_I(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1)

      dtau(:)   = l_contrib * chi_I(:)

     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      chi_I(:) = NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)
      ! Epaisseur optique
      ! chi_p contains both thermal and continuum scattering extinction
      dtau(:)   =  l_contrib * chi_I(:)

    end if

	 !do not take into account dtau of the cell if we integrate tau in the same direction of propagation
	 !but in opposite direction yes.
	 if (to_obs == 1) then !going in the same direction as the ray
	   !do not take into account the path in the cell for this case as we want tau at the edge
       if (nbr_cell > 1) tau(:,iray,id) = tau(:,iray,id) + dtau
     else if (to_obs==-1) then !opposite direction, count the travel path in the cell
       tau(:,iray,id) = tau(:,iray,id) + dtau
     else
      write(*,*) to_obs
      CALL ERROR("unknown value for to_obs")
     endif

    end if
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE integrate_tau_bound

END MODULE AtomicTransfer
