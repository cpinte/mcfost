! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atoms
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

MODULE AtomicTransfer

	use metal, only				: metal_bb_new, compute_opacities, alloc_atom_quantities, dealloc_Atom_quantities
	use opacity
	use Planck, only			: bplanck
	use Profiles, only			: Iprofile, Iprofile_thomson, Iprofile_cmf_to_obs, Zprofile, write_profiles_ascii
	use spectrum_type
	use atmos_type
	use readatom
	use lte
	use constant, only			: MICRON_TO_NM
	use collision, only			: CollisionRate !future deprecation
	use impact
	use solvene
	use statequil_atoms
	use init_solution, only		: Init_NLTE, free_NLTE_sol, gpop_old, Tex_old, flatpops, lcell_converged
	use accelerate
	use voigtfunctions
	use writeatom, only			: writePops, writeelectron, writehydrogendensity
	use write_opacity!, only 				: write_Jnu, write_taur, write_atom_xsections_bf_ff
	use math
	!$ use omp_lib

 !MCFOST's original modules
	use input, only				: lkeplerian, linfall, RT_line_method
	use parametres
	use grid
	use dust_transfer, only		: compute_stars_map
	use dust_ray_tracing, only	: init_directions_ray_tracing ,            &
                              			  tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az, &
                              			  tab_RT_incl, stars_map, kappa, stars_map_cont
	use stars
	use wavelengths
	use mcfost_env, only		: dp
	use constantes, only		: tiny_dp, huge_dp

	IMPLICIT NONE
 
	real(kind=dp), parameter :: tiny_chi = 1d-50
 !Pointer to formal solver
	PROCEDURE(INTEG_RAY_LINE_I), pointer :: INTEG_RAY_LINE => NULL()
 !Temporary variable for Zeeman calculations
	real(kind=dp), dimension(:,:), allocatable :: QUV
 !Temporary variables for Contribution functions
	real(kind=dp), allocatable :: S_contrib(:,:,:), S_contrib2(:,:), mean_formation_depth

 CONTAINS

  
	SUBROUTINE INTEG_RAY_LINE_I(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	! This routine performs integration of the transfer equation along a ray
	! crossing different cells.
	! --> Atomic Lines case.
	! ------------------------------------------------------------------------------- !

		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs !used in NLTE but why?
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
		real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, chiI, LimbD, chiIc
		real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau
		real(kind=dp) :: etau, edtau, etau_c, edtau_c
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp
		tau_c(:) = 0.0_dp
		chiI(:) = 0.0_dp
		chiIc(:) = 0.0_dp


		NLTEspec%I(:,iray,id) = 0d0
		NLTEspec%Ic(:,iray,id) = 0d0


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
    !write(*,*) "Boucle infinie, icell=", icell

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
       !AT the moment the stellar Flux is only a BB
					call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)
       !if (maxval(limbD == 0.0_dp)) return
       !CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
					NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
       				NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau_c)
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
				l_contrib = l_contrib * AU_to_m !l_contrib in m
				if ((nbr_cell == 1).and.labs) ds(iray,id) = l * AU_to_m

				eval_operator = (labs .and. (nbr_cell == 1)) !labs if false for images
     											  !so no pb if Nact>0 and we use a different grid
				CALL initAtomOpac(id)
				if (atmos%NpassiveAtoms>0) &
					CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)

				chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell) + tiny_chi 
				chiI(:)  = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)

				if (atmos%NactiveAtoms>0) then 
					!(line+continuum) LTE opacities, not updated during the local sub iterations.
					if (eval_operator) then 
						chi_loc(:,iray,id) = chiI(:)
						eta_loc(:,iray,id) = NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) + &
						NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell)
						!Do not updated the scattering emissivity for local iterations ?? better ??
						
						!Because Kc_nlte is the total nlte b-f, stored on the whole grid
						!if (iray==1) call calc_etac_atom_loc(id, icell)
						!-> needed only if MALI ? since Stot is computed with Psi or locally 
					endif
					
					CALL initAtomOpac_nlte(id)
      
					CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
					
					!NLTE bound-free computed before propagation on the whole grid.
					!Or Before each iterations
					chiIc(:) = chiIc(:) + NLTEspec%AtomOpac%Kc_nlte(:,icell)

					chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id) +  NLTEspec%AtomOpac%Kc_nlte(:,icell)
				endif !active atoms

				dtau(:)   = l_contrib * chiI(:)
				dtau_c(:) = l_contrib * chiIc(:)
            
				Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
				Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
      
      			!Always include electron scattering if NLTE loop.
      			!But ATM, uses %Jc even for lines.
				if (atmos%Nactiveatoms>0) then 
					Snu = Snu + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell) ) / chiI(:)
					Snu_c = Snu_c + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + &
					NLTEspec%AtomOpac%jc_nlte(:,icell) )/ chiIc(:)

				else if ((atmos%electron_scattering).and.(atmos%NactiveAtoms==0)) then
					Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
					Snu_c = Snu_c + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) / chiIc(:)
				endif


! 				NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
! 				NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c
				do la=1, NLTEspec%Nwaves
! 					if (tau(la) < ..) then
! 						etau = 1 + tau(la) ? 
! 					...
				
					NLTEspec%I(la,iray,id) = NLTEspec%I(la,iray,id) + dexp(-tau(la)) * (1.0_dp - dexp(-dtau(la))) * Snu(la)
					NLTEspec%Ic(la,iray,id) = NLTEspec%Ic(la,iray,id) + dexp(-tau_c(la)) * (1.0_dp - dexp(-dtau_c(la))) * Snu_c(la) 
				
				enddo
     
				if (eval_operator) then
					!!dtau_ds = ds(iray,id) * chiI(:)
					!!CALL calc_psi_operator(id, icell, iray, chiI,ds(iray,id) * chiI(:),Snu)
					NLTEspec%S(:,iray,id) = Snu(:)
					NLTEspec%chi(:,iray,id) = chiI(:)
				end if
     
				tau = tau + dtau
				tau_c = tau_c + dtau_c

			end if  ! lcellule_non_vide
		end do infinie

	RETURN
	END SUBROUTINE INTEG_RAY_LINE_I

 !building and ACTIVE opacities
 SUBROUTINE INTEG_RAY_LINE_Z(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! This routine is "LTE", although NLTE pops can be used
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, LimbD
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chiI, chiIc
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
       call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)
       !CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau_c)
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
     CALL initAtomOpac_zeeman(id)

      if (atmos%NpassiveAtoms>0) CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)

      chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell) + tiny_chi 
      if (atmos%NactiveAtoms > 0) chiIc(:) = chiIc(:) + NLTEspec%AtomOpac%Kc_nlte(:,icell)
      chiI(:) = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)

      if (atmos%NactiveAtoms>0) then 
       CALL initAtomOpac_nlte(id)

       CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)

       chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id)
      endif
                    
      dtau(:)   = l_contrib * chiI(:)
      dtau_c(:) = l_contrib * chiIc(:)
      
      Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
      Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
      
      			!Always include electron scattering if NLTE loop.
      			!But ATM, uses %Jc even for lines.
				if (atmos%Nactiveatoms>0) then 
					Snu = Snu + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell) ) / chiI(:)
					Snu_c = Snu_c + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + &
					NLTEspec%AtomOpac%jc_nlte(:,icell) )/ chiIc(:)
				else if ((atmos%electron_scattering).and.(atmos%NactiveAtoms==0)) then
					Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
					Snu_c = Snu_c + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) / chiIc(:)
				endif

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
           CALL move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
             CALL INTEG_RAY_LINE(id, icell, x0,y0,z0,u0,v0,w0,iray,labs)

             I0 = I0 + NLTEspec%I(:,iray,id) / npix2
             I0c = I0c + NLTEspec%Ic(:,iray,id) / npix2
             
             if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
             	QUV(3,:) = QUV(3,:) + NLTEspec%STokesV(:,iray,id)
             	QUV(1,:) = QUV(1,:) + NLTEspec%STokesQ(:,iray,id)
             	QUV(2,:) = QUV(2,:) + NLTEspec%STokesU(:,iray,id)
             end if
             !loop over all directions for each cell-frequency
             if (lcontrib_function) S_contrib2(:,:) = S_contrib2(:,:) + S_contrib(:,:,id)

           !else !Outside the grid, no radiation flux
           endif
        end do !j
     end do !i

     I0 = I0 / npix2
     I0c = I0c / npix2

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

  ! Flux out of a pixel in W/m2/Hz/pix
  normF = (pixelsize / (distance*pc_to_AU) )**2

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
                         NLTEspec%Ksi(:,:,ibin,iaz) + S_contrib2(:,:)


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

  integer, parameter :: n_rad_RT = 600, n_phi_RT = 100 !(100, 36)
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
    rmax_RT = Rmax !* 2.0_dp

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
    !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz,etoile)
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
           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
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
     !do i=npix_x_max/2+1, npix_x_max/2+1
        !$ id = omp_get_thread_num() + 1

        do j = 1,npix_y
        !do j=npix_y/2+1,npix_y/2+1
        !write(*,*) "ipix, jpix", i, j
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
  integer :: icell, m
  integer :: ibin, iaz
  integer, parameter :: Nrayone = 1
  character(len=20) :: ne_start_sol = "H_IONISATION"
  character(len=20)  :: newPRT_SOLUTION = "FULL_STOKES"
  
 ! -------------------------------INITIALIZE AL-RT ------------------------------------ !
  !only one available yet, I need one unpolarised, faster and more accurate.
  Voigt => VoigtHumlicek
  !Voigt => dirac_line
  !Profile => Iprofile_cmf_to_obs
  Profile => IProfile
  !Profile => Iprofile_thomson
  INTEG_RAY_LINE => INTEG_RAY_LINE_I



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

 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
        !! ----------------------- Read Model ---------------------- !!

  if (.not.lpluto_file) then
   if (lmodel_ascii) then
    CALL readAtmos_ascii(density_file)
   end if
  end if
        !! --------------------------------------------------------- !!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------- READATOM and INITIZALIZE POPS ------------------------------ !

  CALL readAtomicModels(atomunit)
  
  if (atmos%NactiveAtoms > 0) then 
   atmos%Nrays = Nrays_atom_transfer!100!0
  else
   atmos%Nrays = Nrayone
   if (lelectron_scattering) atmos%Nrays = 2000 !here before I is allocated and fixed !
  endif

  !compute first guess of electron density ??
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
 ! ---------- INITIALIZE WAVELNGTH, BACKGROUND OPAC AND ATOMC QUANTITIES -------------- !
 ! ------------------------------------------------------------------------------------ !
  atmos%electron_scattering=lelectron_scattering !but also H and He, futur deprec of atmos
  if (ltab_wavelength_image) NLTEspec%write_wavelength_grid = .true.
  !call set_up_sub_wavelengths_grid()
  CALL initSpectrum(lam0=500d0,vacuum_to_air=lvacuum_to_air)
  
 ! ----------------------------  INITIAL POPS------------------------------------------ !
  if (atmos%NactiveAtoms > 0) then
  		!Allows to have atom%n init which is needed in compute_atom_quantities() (compute_opacities)
		CALL Init_NLTE(sub_iterations_enabled=.true.)
		!free at the end of NLTEloop()
  endif
 !  ------------------------------------------------------------------------------------ !

  if ((atmos%NactiveAtoms>0) .or. .not.(ltab_wavelength_image)) then
   if (n_etoiles > 0) CALL init_stellar_disk !for all wavelengths, all stars at disk centre
   write(*,*) " Computing background opacities..."
   CALL alloc_atom_quantities !call compute_atom_quantities(icell) for each cell
   CALL compute_opacities
   write(*,*) " ..done"
  endif
 ! ------------------------------------------------------------------------------------ !
 ! --------------------- ADJUST MCFOST FOR STELLAR MAP AND VELOCITY ------------------- !
  ! ----- ALLOCATE SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ------!
  CALL reallocate_mcfost_vars() !assumes more than 1 wavelength otherwise delta_wl is 0!
  ! --- END ALLOCATING SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT --!
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------------- Jny ----------------------------------------- !
   !IF NLTE loop we should also be able to read or compute an estimate of Jcont by the way...
   ! -> TBD
   if (atmos%NactiveAtoms == 0 .and. lelectron_scattering) then
!    call alloc_jnu(.false.) -> already allocated in alloc_spectrum
!    if (.not.lread_jnu_atom) then
	if (lread_jnu_atom) then
		call read_jnu_ascii
!      call read_jnu() !I now write Jnu in ascii file to be interpolated if lread_jnu_atom
!	   the interpolated version of Jnu is written at the end in write_jnu()
    else
     call iterate_Jnu()
     call write_Jnu
     if (lstop_after_jnu) then
      write(*,*) " Jnu calculation done." 
      stop
     endif
    endif
   endif
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- NLTE LOOP -------------------------------------- !
  if (atmos%Nactiveatoms > 0) then
     CALL NLTEloop()
  endif
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------- WRITE CONVERGED POPULATIONS ------------------------------ !
 ! ------------------------------------------------------------------------------------ !
  !! alternatively, could be invoked earlier, if the system send a message to stop the code
  !! before convergence, so that we restart with these pops.
  !! Or define a function to store pops every N iter.
  do icell=1,atmos%Natom
   CALL writePops(atmos%Atoms(icell)%ptr_atom)
  end do
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- MAKE IMAGES ------------------------------------ !
 ! ------------------------------------------------------------------------------------ !

   !should add an other flag to force to avoid continua lines
   !Define a wavelength grid for image with only lines, if not input wavelength
!    if (.not.(ltab_wavelength_image)) then !or for all passive lines also
!     ltab_wavelength_image = .true.
!     !!CALL write_wavelengths_table_NLTE_lines(NLTEspec%lambda) !!after grid creation actually
!     tab_wavelength_image = "line_waves.s"
!    endif
  
  !not needed, but avoid to reset it, if only lte, it is already integ_ray_line_i
  if (atmos%NactiveAtoms > 0) then
   INTEG_RAY_LINE => NULL()
   INTEG_RAY_LINE => INTEG_RAY_LINE_I
  endif

  if (ltab_wavelength_image) then
   NLTEspec%atmos%Nrays = Nrayone

   CALL initSpectrumImage() !deallocate waves arrays/ define a new grid also resample Jnu if any
   if (n_etoiles > 0) CALL init_stellar_disk 
   write(*,*) " Computing continuum opacities for image..."
   if (atmos%NactiveAtoms >0) CALL dealloc_atom_quantities
   CALL alloc_atom_quantities
   CALL compute_opacities !recompute background opac
   if (atmos%NactiveAtoms > 0) then
    CALL compute_nlte_bound_free
    !! can be added to Kc and jc to save memory for image
    !NLTEspec%AtomOpac%Kc(:,:,1) = NLTEspec%AtomOpac%Kc(:,:,1) + NLTEspec%AtomOpac%Kc_nlte(:,:)
    !NLTEspec%AtomOpac%jc(:,:,1) = NLTEspec%AtomOpac%jc(:,:,1) + NLTEspec%AtomOpac%jc_nlte(:,:)
	!deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)
   endif
   write(*,*) " ..done"
   !add NLTE contiuna
   !!CALL reallocate_mcfost_vars() !wavelength only?
   CALL reallocate_mcfost_wavelength_arrays()
  else
    !same wave grid, Jnu the same if any
    !if (atmos%NactiveAtoms > 0) then
     !NLTEspec%AtomOpac%Kc(:,:,1) = NLTEspec%AtomOpac%Kc(:,:,1) + NLTEspec%AtomOpac%Kc_nlte(:,:)
     !NLTEspec%AtomOpac%jc(:,:,1) = NLTEspec%AtomOpac%jc(:,:,1) + NLTEspec%AtomOpac%jc_nlte(:,:)
	 !deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)
	!endif
    CALL reallocate_rays_arrays(Nrayone) !rellocate rays array, but keep the same wavelength grid
  end if

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
   write(*,*) "   >>> including contribution function Ksi..."
   allocate(S_contrib(NLTEspec%Nwaves,n_cells,NLTEspec%NPROC))
   allocate(S_contrib2(NLTEspec%Nwaves,n_cells))
   !allocate(mean_formation_depth(npix_x, npix_y))
   S_contrib(:,:,:) = 0d0; S_contrib2(:,:) = 0d0
   INTEG_RAY_LINE => NULL()
   INTEG_RAY_LINE => INTEG_RAY_LINE_I_CNTRB
  end if

  !Actual flux calculation
  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
       CALL EMISSION_LINE_MAP(ibin,iaz)
     end do
  end do
  write(*,*) " ..done"
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------- WRITE RESULTS -------------------------------------- !
 ! ------------------------------------------------------------------------------------ !
  write(*,*) "Writing result to file..."
  call WRITE_FLUX_ASCII() !for testing
  if (3*size(NLTEspec%Flux)/(1024.**3) <= 6.) call WRITE_FLUX
  deallocate(NLTEspec%Flux, NLTEspec%Fluxc) !free here to release a bit of memory
  if (lcontrib_function) then
   CALL WRITE_CNTRB_FUNC_pix()
   deallocate(S_contrib, S_contrib2, NLTEspec%Ksi)
 end if
    
  !call write_contrib_lambda_ascii(654.0_dp,0.0_dp,0.0_dp,1.0_dp)
  !call write_contrib_lambda_ascii(1083.0_dp,0.0_dp,0.0_dp,1.0_dp)
  call write_contrib_lambda_ascii(500.0_dp,0.0_dp,0.0_dp,1.0_dp)
!   call write_taur(500._dp,0._dp,0._dp,1.0_dp)
!-> building
!  call write_contrib_lambda(654.0_dp,0.0_dp,0.0_dp,1.0_dp)
!   call write_contrib_lambda(1080._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp)

  !write it after because sca_c (and maybe %Kc, %jc, %Kc_nlte, %Jc_nlte) is modified !!
  do icell=1, atmos%Natom
    CALL write_cont_opac_ascii(atmos%Atoms(icell)%ptr_atom)
    !!CALL write_cont_opac(atmos%Atoms(icell)%ptr_atom) !-> building
    !!CALL write_atom_xsections_bf_ff(atmos%Atoms(icell)%ptr_atom)
    !!CALL write_profiles_ascii(52, atmos%Atoms(icell)%ptr_atom, 1)
  enddo

  if (atmos%electron_scattering) then
   !Jnu is written to ascii file if read
   if ((atmos%NactiveAtoms > 0).or.(atmos%NactiveAtoms==0.and.lread_jnu_atom)) call write_Jnu
   CALL dealloc_Jnu()
  endif

 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------- CLEANING ------------------------------------------ !
  if (allocated(QUV)) deallocate(QUV)
 !close file after NLTE loop
!Temporary: because we kept in memory, so file is closed earlier
!  do nact=1,atmos%Nactiveatoms
!   CALL closeCollisionFile(atmos%ActiveAtoms(nact)%ptr_atom) !if opened
!  end do
 !CALL WRITEATOM() !keep C in memory for that ?
 CALL freeSpectrum() !deallocate spectral variables
 CALL free_atomic_atmos()
 NULLIFY(Profile, Voigt, INTEG_RAY_LINE)
 !CALL dealloc_profile_line_profiles

 RETURN
 END SUBROUTINE
 ! ------------------------------------------------------------------------------------ !

	SUBROUTINE NLTEloop() !for all active atoms
	! -------------------------------------------------------- !
	! Descriptor here
	! -------------------------------------------------------- !
#include "sprng_f.h"


		integer, parameter :: n_rayons_start = 1000 
		integer, parameter :: n_rayons_start2 = 1000
		integer :: n_rayons_max
		integer :: etape, etape_start, etape_end, iray, n_rayons
		integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, la, kr
		integer, dimension(nb_proc) :: max_n_iter_loc
		logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged
		real :: rand, rand2, rand3
		real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt
		real(kind=dp) :: diff, norme, dT_max, dN, dN1, dJ, lambda_max
		real(kind=dp) :: dT, dN2, dN3, diff_old
		real(kind=dp), allocatable :: dM(:), Jold(:,:), dTM(:), Tion_ref(:), Tex_ref(:)
                                 !futur global flag
		logical :: labs, iterate_ne = .false., accelerated, ng_rest
		logical :: verbose, l_unconverged
		integer :: nact, maxIter, imax, icell_max
		integer :: icell, iorder, i0_rest, n_iter_accel
		integer :: Nlevel_total = 0, NmaxLevel, ilevel, max_sub_iter
		character(len=20) :: ne_start_sol = "NE_MODEL"
! 		real(kind=dp), dimension(3, atmos%Nrays, nb_proc) :: xyz0, uvw0
		type (AtomType), pointer :: atom


		write(*,*) "   -> Solving for kinetic equations for ", atmos%Nactiveatoms, " atoms"
		write(*,*) " Max error : ", dpops_max_error, dpops_sub_max_error
		

open(unit=unit_invfile, file=invpop_file, status="unknown")
write(unit_invfile,*) n_cells

  
		verbose = .true.

		ds = 0.0_dp !meters
		chi_loc = 0.0_dp

  !move to initSol
		allocate(dM(atmos%Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
		allocate(dTM(atmos%Nactiveatoms)); dM=0d0 !keep tracks of Tex for all cells for each atom
		allocate(Jold(NLTEspec%Nwaves, atmos%Nspace))!move to initsol
		allocate(Tex_ref(atmos%Nactiveatoms)); dM=0d0 !keep tracks of max Tex for all cells for each line of each atom
		allocate(Tion_ref(atmos%Nactiveatoms)); dM=0d0 !keep tracks of max Tion for all cells for each cont of each atom
		diff_old = 0.0_dp

		n_rayons_max = atmos%Nrays
! 		xyz0(:,:,:) = 0d0
! 		uvw0(:,:,:) = 0d0

		labs = .true. !to have ds at cell icell = eval_operator
		id = 1
		etape_start = 2
		etape_end = 2
		if (etape_start==1) then
			write(*,*) "Warn: etape 1 not accurate"
		endif

		accelerated  = .false.
!   if (lNg_acceleration) then
!    CALL error("+Ng no tested yet")
!    n_iter_accel = 0
!    i0_rest = 0
!    ng_rest = .false.
!   endif

		iterate_ne = (n_iterate_ne>0)
		if (iterate_ne) then
			write(*,*) " before iterate ne I need te recompute gij for continua !!"
		endif


		max_sub_iter = 1000
		maxIter = 1000
		
!-> move after LTE pops
! 	----------------------------  INITIAL POPS------------------------------------------ !
! 		CALL Init_NLTE(sub_iterations_enabled=.true.)
! 	 ------------------------------------------------------------------------------------ !

		do etape=etape_start, etape_end

			if (etape==1) then !two rays, not accurate
				lfixed_rays=.true.
				n_rayons = 2
				iray_start = 1
				lprevious_converged = .false.

			else if (etape==2) then 
				lfixed_rays = .true.
				n_rayons = n_rayons_max!min(n_rayons_max,n_rayons_start)
				iray_start = 1
				lprevious_converged = .false.

			else
				CALL ERROR("etape unkown")
			end if
  	  
			write(*,*)  "-> Using ", n_rayons, ' rays for step', etape

			lnotfixed_rays = .not.lfixed_rays
			lconverged = .false.
			n_iter = 0

			do while (.not.lconverged)

				n_iter = n_iter + 1
				if (n_iter > maxIter) exit !change step
				write(*,*) " -> Iteration #", n_iter, " Step #", etape

				if (lfixed_rays) then
					stream = 0.0
					do i=1,nb_proc
						stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
					end do		
				end if

				max_n_iter_loc = 0
				dT_max = 0.
            
            
				!Initialize some quantities that change from one iteration to another,
				!but which are constant for sub iterations

				!$omp parallel &
				!$omp default(none) &
				!$omp private(icell, id, atom,nact) &
				!$omp shared(gpop_old,Tex_old, atmos) &
				!$omp shared (Jold, NLTEspec)
				!$omp do schedule(static,1)
				do icell=1, atmos%Nspace
					!$ id = omp_get_thread_num() + 1
					if (atmos%icompute_atomRT(icell)>0) then
						!For the propagation we need this opacity at all cells
						!It doesn't change with rays. Will ne updated cell by cell for local iterations
						!call compute_atom_quantities(icell) ? update, gij with new nstar and ne
						!update profiles and adamp
						call NLTE_bound_free(icell) !takes the new populations, after sub it eventually 
						do nact=1, atmos%NactiveAtoms
							atom => atmos%ActiveAtoms(nact)%ptr_atom
							gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
							do kr=1,atom%Nline
								Tex_old(nact, kr, icell) = atom%lines(kr)%Tex(icell)
							enddo
							do kr=1, atom%Ncont
								Tex_old(nact, kr+Atom%Nline,icell) = atom%continua(kr)%Tex(icell)
							enddo
							atom => NULL()
             		!!Recompute Profiles if needed here
						end do
					end if
					!-> Would be replaced later, when I'll know what to do for electron scattering
					Jold(:,icell) = NLTEspec%Jc(:,icell)
				end do
				!$omp end do
				!$omp end parallel
        	
				!$omp parallel &
				!$omp default(none) &
				!$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dM, dN, dN1,dT_max)&
				!$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, nact, atom, l_unconverged) &
				!$omp shared(atmos,NLTEspec, dpops_sub_max_error, verbose,lkeplerian,n_iter) & !!xyz0, uvw0 & !before nact was shared
				!$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,max_sub_iter,lcell_converged) &
				!$omp shared(n_cells, gpop_old,lforce_lte, INTEG_RAY_LINE) &
				!$omp shared(lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape,pos_em_cellule)
				!$omp do schedule(static,1)
				do icell=1, n_cells
					!$ id = omp_get_thread_num() + 1
					if (atmos%icompute_atomRT(icell)>0) then

						CALL fill_Collision_matrix(id, icell) !Gamma(i,j) = C(i,j) before computing rates

                    	!Psi, exp(-dtau) and eta_line(lambda,ray,proc) and eta_cont(lambda, proc)
						!!CALL init_psi_operator(id)
                    
						do iray=iray_start, iray_start-1+n_rayons
							if (etape==1) then
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

							else !etape 2
                   	    ! Position aleatoire dans la cellule
								rand  = sprng(stream(id))
								rand2 = sprng(stream(id))
								rand3 = sprng(stream(id))

								CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
								rand = sprng(stream(id))
								W0 = 2.0_dp * rand - 1.0_dp !nz
								W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
								SRW02 = sqrt(W02)
								rand = sprng(stream(id))
								ARGMT = PI * (2.0_dp * rand - 1.0_dp)
								U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
								V0 = SRW02 * sin(ARGMT) !ny = sin(theta) * sin(phi)
							end if !etape

						!keep them for subiterations or keep phi(:,iray,id)
						!Because the shift is what matters now, I don't need them ?
		
! 							xyz0(1,iray,id) = x0; xyz0(2,iray,id) = y0; xyz0(3,iray,id) = z0
! 							uvw0(1,iray,id) = U0; uvw0(2,iray,id) = V0; uvw0(3,iray,id) = W0					
! 							
! 							CALL INTEG_RAY_LINE(id, icell, xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
! 											uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, labs)

							CALL INTEG_RAY_LINE(id, icell, x0,y0, z0, u0, v0, w0, iray, labs)
											
						enddo !iray
						!always in NLTEloop, although %J is not used ATM.
						!Done here, otherwise, J is change in integration
						!Pb: What if the transfer converges before J ?
						!if force_lte, J is not updated because there is only one iteration.
						CALL calc_Jnu(id, icell, n_rayons)
						!Is it this one for continua radiative rates ?

						if (lforce_lte) then
							call initGamma(id) !because we do not enter the local sub iterations
							!so Gamma need to be initialized to the collision matrix.
							call update_populations(id, icell, diff, verbose,n_iter)
							lconverged_loc = .true.
						else
							lconverged_loc = .false.
							!Will start at LTE, since radiative rates are added in the sub iterations.
						endif

						n_iter_loc = 0
     				!!only iterate on cells which are not converged
						do while (.not.lconverged_loc)
							n_iter_loc = n_iter_loc + 1
                        
						!Computes Radiative rates
						!if Jnu is used instead of J accelerated for continua, can be spread in two.
							do nact=1, atmos%Nactiveatoms
								atom => atmos%Activeatoms(nact)%ptr_atom
								call initGamma_atom(id, atom)
								call init_rates_atom(id, atom) !init to 0 before accumulate
															   !except line%Rji = line%Aji
								do iray=iray_start, iray_start-1+n_rayons
									CALL calc_rates_atom(id, icell, iray, atom, n_rayons)
								enddo
								atom => NULL()
							enddo

							do nact=1, atmos%Nactiveatoms
								atom => atmos%Activeatoms(nact)%ptr_atom
                         !fill rate matrix with radiative rates
								CALL rate_matrix_atom(id, icell, atom, lforce_lte)
								atom => NULL()
							enddo
                        
                        !Solve populations and remove delta(l,l')
							if (mod(n_iter_loc,max(1,max_sub_iter/10))==0) then
								write(*,*) " -- Global iteration #", n_iter, " step #", etape, "eps: ", dpops_sub_max_error
							endif
							CALL update_Populations(id, icell, diff,(mod(n_iter_loc,max(1,max_sub_iter/10))==0),n_iter_loc)
							if (mod(n_iter_loc,max(1,max_sub_iter/10))==0) write(*,*) " --------------------------------------------------------------- "

                        !for this icell (id) and local iteration
							if (diff < dpops_sub_max_error) then!precision_sub
								lconverged_loc = .true.
     					    !write(*,*) "converged", id, icell, "dpops(sub) = ", diff
							else
								!!CALL init_psi_operator(id) !init also NLTEspec%S

								!etac is now updated in calc_eta_atom_loc if iray==1.
  			        			!and chi = chi(LTE) + chiline_nlte (new) + chic_nlte (new)
  			        			!At the start of the next iteration, NLTE_bound_free is recomputed
  			        			!and saved for the whole grid.
  			        			!chic_nlte is stored for each id, to be added for each rays, even
  			        			!if computed at iray==1. It is reset each time iray ==1 (for an id/icell)
  			        		
								do iray=iray_start, iray_start-1+n_rayons
      							!I unchanged
      							
      							
      							!also adds chi_loc, which is takes LTE opacities, that do not change
      							!Also compute atom%etac and chic if iray==1
      							!These are then added to atom%eta and chitot for each ray
									!CALL calc_eta_atom_loc(id,icell,iray, (iray==1))
									
									CALL calc_total_source_loc(id,icell,iray, (iray==1))
      						    

								enddo !iray

							end if
							if (n_iter_loc >= max_sub_iter) then 

								lconverged_loc = .true.

							end if

						end do !local sub iteration
						if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
     	            
					end if !icompute_atomRT
				end do !icell
				!$omp end do
				!$omp end parallel
        	
     		!Global convergence Tests

            !should be para
				dM(:) = 0d0
				diff = 0d0
				dJ = 0.0
				dT = 0.0
				dTM(:) = 0.0
				Tex_ref(:) = 0.0
				Tion_ref(:) = 0.0
				cell_loop2 : do icell=1, atmos%Nspace
					if (atmos%icompute_atomRT(icell)>0) then
						dN = 0d0 !for all levels of all atoms of this cell
						dN2 = 0.0
						do nact=1,atmos%NactiveAtoms
							atom => atmos%ActiveAtoms(nact)%ptr_atom
     					         					    
							do ilevel=1,atom%Nlevel
     						
								if (atom%n(ilevel,icell) >= prec_pops) then
									dN1 = dabs(1d0-gpop_old(nact,ilevel,icell)/atom%n(ilevel,icell))
     				    		 !dabs(atom%n(l,icell) - ndag(l))/atom%n(l,icell)
									dN = max(dN1, dN) !compare with all atoms and levels
     				    						  !for this cell
									dM(nact) = max(dM(nact), dN1) !compare for one atom
     				    		 							   !for all cells
								endif
							end do !over ilevel
							do kr=1, atom%Nline
							
								dN3 = dabs(1.0 - Tex_old(nact, kr, icell) / atom%lines(kr)%Tex(icell))
								dN2 = max(dN3, dN2)
								!dTM(nact) = max(dTM(nact), dN3)
								if (dN3 > dTM(nact)) then
									dTM(nact) = dN3
									Tex_ref(nact) = atom%lines(kr)%Tex(icell)
									icell_max = icell
								endif
							enddo
							do kr=1, atom%Ncont
								dN3 = dabs(1.0 - Tex_old(nact, kr+atom%Nline, icell) / atom%continua(kr)%Tex(icell))
								dN2 = max(dN3, dN2)
								!dTM(nact) = max(dTM(nact), dN3)
								if (dN3 > dTM(nact)) then
									dTM(nact) = dN3
									Tion_ref(nact) = atom%continua(kr)%Tex(icell)
									icell_max = icell
								endif
							enddo
							atom => NULL()
						end do !over atoms
						!compare for all atoms and all cells
						!diff = max(diff, dN) ! pops
						diff = max(diff, dN2) ! Tex
     					
						dN1 = dabs(1d0 - maxval(Jold(:,icell)/(tiny_dp + NLTEspec%Jc(:,icell))))
						if (dN1 > dJ) then
							dJ = dN1
							imax = locate(dabs(1d0 - Jold(:,icell)/(tiny_dp + NLTEspec%Jc(:,icell))), dJ)
							lambda_max = NLTEspec%lambda(imax)
						endif
						!lcell_converged(icell) = (real(diff) < dpops_max_error)
						!diff = max(diff, dJ) ! ? to have bother converged
     					
					end if
				end do cell_loop2 !icell
     		
				write(*,*) " --------------------------------------------------------------- "
				do nact=1,atmos%NactiveAtoms
					write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dT = ", dTM(nact)
					if (Tion_ref(nact) /= 0.0) write(*,*) "       <::>  Tion (K) = ", Tion_ref(nact),  " Te (K) = ", atmos%T(icell_max)
					if (Tex_ref(nact) /= 0.0) write(*,*) "       <::>  Tex (K) = ", Tex_ref(nact),  " Te (K) = ", atmos%T(icell_max)
					if (accelerated) then
						write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dM = ", dM(nact)," (Accelerated)"
					else
						write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dM = ", dM(nact)
					endif
				enddo
				write(*,*) " -> diff =", diff, " old = ", diff_old !at the end of the loop over n_cells
				write(*,*) " -> dJ = ", dJ, " @ ", lambda_max, " nm"
				write(*,*) maxval(max_n_iter_loc), "sub-iterations"
				write(*,*) " --------------------------------------------------------------- "
				diff_old = diff

				!Not used if the next is not commented out
				lconverged = (real(diff) < dpops_max_error) !global convergence for all iterations

				if (real(diff) < dpops_max_error) then !precision
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
              				if (n_iter >= maxIter) then
             		 			write(*,*) "Warning : not enough rays to converge !!"
                 			lconverged = .true.
              				end if
              			end if

              ! On continue en calculant 2 fois plus de rayons
              ! On les ajoute a l'ensemble de ceux calcules precedemment
!              iray_start = iray_start + n_rayons

          	   		end if
        		end if
        									
			!Not optimal, but fast relatively
			!must be parra
				if (iterate_ne .and. (mod(n_iter,n_iterate_ne)==0))  then
					write(*,*) " Solve ne Global iteration:", n_iter
					write(*,*) " --> old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
					do icell=1, n_cells
						if (atmos%icompute_atomRT(icell) > 0) then
							NLTEspec%AtomOpac%sca_c(:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) - sigma_e * atmos%ne(icell)
							do nact=1,atmos%NactiveAtoms
								atom => atmos%ActiveAtoms(nact)%ptr_atom
								do kr=1, atom%Ncont
									atom%continua(kr)%gij(:,icell) = atom%continua(kr)%gij(:,icell) / atmos%ne(icell)
								enddo
								atom => Null()
							end do
						endif
					enddo
					CALL SolveElectronDensity(ne_start_sol)
					do icell=1, n_cells
						if (atmos%icompute_atomRT(icell) > 0) then
							NLTEspec%AtomOpac%sca_c(:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) + sigma_e * atmos%ne(icell)

							do nact=1,atmos%NactiveAtoms
								atom => atmos%ActiveAtoms(nact)%ptr_atom
								do kr=1, atom%Ncont
									atom%continua(kr)%gij(:,icell) = atom%continua(kr)%gij(:,icell) * atmos%ne(icell)
								enddo
								atom => Null()
							end do
						endif
					enddo
				end if


			end do !while
			write(*,*) "step: ", etape, "Threshold: ", dpops_max_error

		end do !over etapes


		!I recompute for all even for active atoms, for consistency. 
		!Because, only gij for active continua and ne are updated, not %nstar
		if (n_iterate_ne < 0) then
			write(*,*) "END LOOP: old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
			CALL SolveElectronDensity(ne_start_sol)
			!Recompute for all atoms the LTE pops
			write(*,*) "   recompute LTE pops for all atoms"
			do nact=1,atmos%NAtom
				atom => atmos%Atoms(nact)%ptr_atom
				if (atom%ID=="H") then
					call LTEpops_H()
				else
					CALL LTEpops(atom,.true.)
				endif
				atom => Null()
			end do
		endif

		if (iterate_ne .or. (n_iterate_ne < 0)) then
			CALL writeElectron(.true.) !the .true. means append, to compare with initial solution.
			!Recompute Background  + NLTE bound free.
			!Check it is not done after
			write(*,*) "Check the consistency of NLTE b-f and background opac with iterate_ne"
		endif
! -------------------------------- CLEANING ------------------------------------------ !
		!to move inside free_nlte_sol
		deallocate(dM, dTM, Tex_ref, Tion_ref)
		if (allocated(Jold)) deallocate(Jold)
		CALL free_nlte_sol(.false.)
		
close(unit=unit_invfile)
! ------------------------------------------------------------------------------------ !
	RETURN
	END SUBROUTINE NLTEloop
! ------------------------------------------------------------------------------------ !

 !building
 SUBROUTINE adjustStokesMode(Solution)
 !
 !
   character(len=*), intent(in) :: Solution


write(*,*) "Note polarized profile not allocated yet, check that in profile and line%phiZ, line%psi"
write(*,*) "Magnetic profile not ready yet"
stop

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
!should be done more properly and I should use functions from mcfostfost for stellar flux
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
 
 !futur move to initial_solution
 SUBROUTINE init_stellar_disk
  integer :: i_star!, lam
  
   write(*,*) " Computing Istar(mu=1) for each star..."
   !move elsewhere if we read a file, but normally should well fit with MCFOST routines
   do i_star=1, n_etoiles
    if (etoile(i_star)%T <= 1e-6) then
     call warning("Setting stellar radiation to 0: T* < 1e-6 K")
     NLTEspec%Istar(:,i_star) = 0._dp
    else
     CALL Bplanck(etoile(i_star)%T*1d0, NLTEspec%Istar(:,i_star))
    endif
   enddo
   write(*,*) " ..done"
 
 RETURN
 END SUBROUTINE init_stellar_disk

 
 SUBROUTINE calc_stellar_surface_brightness(N,lambda,i_star,x,y,z,u,v,w,gamma)
 ! ---------------------------------------------------------------!
  ! Compute the stellar radiation field surface brightness.
  ! Istar = B(x,y,z;mu) * Stellar spectrum or B * BlackBody
  ! return gamma, the brightness. For uniform disk, gamma is 1.
  !
  ! For a point on the stellar disk, returns gamma * limbdarkenig
 ! -------------------------------------------------------------- !
  use input, only : limb_darkening, mu_limb_darkening
  use Planck, only : uLD, Bplanck
  integer, intent(in) :: N, i_star
  real(kind=dp), dimension(N), intent(in) :: lambda
  real(kind=dp), dimension(N), intent(out) :: gamma
  real(kind=dp), intent(in) :: u, v, w, x, y, z
  real(kind=dp) :: energie(N)
  real(kind=dp) :: mu, ulimb, LimbDarkening, surface, HC
  integer :: ns
  logical :: lintersect_spot

   if (etoile(i_star)%T <= 1e-6) then
    gamma(:) = 0.0_dp
    return
   else
    gamma(:) = 1d0
   endif
   
   !cos(theta) = dot(r,n)/module(r)/module(n)
   mu = abs(x*u + y*v + z*w)/dsqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
   if (real(mu)>1d0) then !to avoid perecision error
    write(*,*) "mu=",mu, x, y, z, u, v, w
    CALL Error(" mu limb > 1!")
   end if
   
   !1) Compute stellar flux from mcfost 
   ! .... 
   
   !2) Correct with the contrast gamma of a hotter/cooler region if any
   CALL intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect_spot)
   if (lintersect_spot) then
     gamma(:) = (dexp(hc_k/lambda/real(etoile(i_star)%T,kind=dp))-1)/&
     			(dexp(hc_k/lambda/etoile(i_star)%SurfB(ns)%T)-1)
     !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
   end if

   !3) Apply Limb darkening
   if (llimb_darkening) then
     CALL ERROR("option for reading limb darkening not implemented")
   else
     !write(*,*) maxval(uLD(real(etoile(i_star)%T,kind=dp))), minval(uLD(real(etoile(i_star)%T,kind=dp)))
     ulimb = 0.0 ! could use BB slope
     LimbDarkening = 1d0 - ulimb*(1d0-mu)
   end if
   !Istar(:) = energie(:) * LimbDarkening * gamma(:)
   gamma(:) = LimbDarkening * gamma(:)

 RETURN
 END SUBROUTINE calc_stellar_surface_brightness

   SUBROUTINE INTEG_RAY_JNU(id,icell_in,x,y,z,u,v,w,iray,labs, kappa_tot, Snu, Istar, psi, Ic)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray to
  ! compute coherent Jnu in the continuum
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(in), dimension(:,:) :: kappa_tot, Snu
  real(kind=dp), intent(in), dimension(:) :: Istar
  real(kind=dp), intent(out) :: Ic(:,:), psi(:,:)
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
!   real(kind=dp) :: chil(size(Ic(:,1)), etal(size(Ic(:,1))
  real(kind=dp), dimension(size(Ic(:,1))) :: LimbD
  real(kind=dp), dimension(size(Ic(:,1))) :: tau_c!, tau
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  !tau(:) = 0.0_dp !go from surface down to the star
  tau_c(:) = 0.0_dp

  Ic(:,id) = 0.0
  psi(:,id) = 0.0
  
  eval_operator = .false.
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
     !lcellule_non_vide=.true.
     !atmos%nHtot(icell)>0
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif
    
    !if (minval(tau_c) > 50.) return
    
    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN

    if (lintersect_stars) then
      if (icell == icell_star) then
       !call calc_stellar_surface_brightness(size(Ic(:,1)),lambda,i_star,x0,y0,z0,u,v,w,LimbD)
       Ic(:,id) =  Ic(:,id) + Istar(:) * dexp(-tau_c)
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
     l_contrib = l_contrib * AU_to_m !l_contrib in m

      Ic(:,id) = Ic(:,id) + dexp(-tau_c) * (1.0_dp - dexp(-l_contrib * kappa_tot(:,icell))) * Snu(:,icell)

!    	  if (any_nan_infinity_vector(dexp(-tau_c) * (1.0_dp - dexp(-l_contrib * kappa_tot(:,icell))) * Snu(:,icell)) > 0) then
!    	   write(*,*) " Error inconsistant values found in I(icell)"
!    	   write(*,*) "icell=", icell, 	l_contrib, atmos%t(icell), atmos%nhtot(icell), atmos%ne(icell)
!    	   write(*,*) "tau=",tau_c
!    	   write(*,*) "K=",kappa_tot(:,icell)
!    	   write(*,*) "Snu=",Snu(:,icell)
! 	   write(*,*) "l*kappa=", l_contrib * kappa_tot(:,icell)
! 	   stop
!    	  endif
     
     if ((nbr_cell == 1).and.labs) then 
      ds(iray,id) = l * AU_to_m
      psi(:,id) = (1d0 - exp(-ds(iray,id)*kappa_tot(:,icell)))
     endif


     !tau = tau + dtau
     tau_c = tau_c + l_contrib * kappa_tot(:,icell)

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_JNU
  
 SUBROUTINE Iterate_Jnu()
 ! -------------------------------------------------------- !
  ! Compute the mean radiation field at all cells
  ! evaluated on a small grid and then interpolated on the
  ! wavelength grid
 ! -------------------------------------------------------- !
#include "sprng_f.h"

  integer, parameter :: n_rayons_start = 2000, maxIter = 1000
  !integer, parameter :: Nlambda = 350
  integer :: n_rayons_max, Nlambda
  real, parameter :: precision = 1e-2
  real, parameter :: lambda_min = 5., lambda_max0 = 100000
  integer :: etape, etape_start, etape_end, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lprevious_converged, write_convergence_file

  real :: rand, rand2, rand3, a1, a0, a2, a3!, fac_etape
  real(kind=dp) :: dSource, dN
  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dJ, diffs
  real(kind=dp), allocatable :: Jold(:,:), Jflat(:), lambda(:), Jnew(:,:), lambda_1(:,:), Jnew_l(:,:), Jold_l(:,:)
  real(kind=dp), allocatable :: Snew(:,:), Kappa_tot(:,:), beta(:,:), Istar(:), Ic(:,:)
  real(kind=dp), allocatable :: lambda_star(:,:), Sth(:,:), Sold(:,:), Sline(:,:)
                                   !futur global flag
  logical :: labs, l_unconverged, accelerated, ng_rest, ng_acc
  integer :: la, icell, imax, icell_max, iorder,i0_rest, icell_max_s, imax_s
  real(kind=dp) :: lambda_max
  type(Ng) :: NgJ
  
	Nlambda = size(NLTEspec%lambda_cont)
  write(*,*) "   --> Lambda iterating Jnu with Nlambda ", Nlambda
  write_convergence_file = .false.
  
  if (allocated(ds)) deallocate(ds)
  allocate(ds(atmos%Nrays, NLTEspec%NPROC))
  ds = 0.0_dp !meters
  
  !this would allow to compute easily Jnu with lines or I need a version
  !of Profile which use as input the wavelength grid
  !call initSpectrum_jnu(Nlambda, lambda_min, lambda_max0)

  !only one star
  allocate(Istar(Nlambda), Ic(Nlambda, NLTEspec%NPROC), lambda_star(Nlambda, NLTEspec%NPROC), lambda_1(Nlambda, NLTEspec%NPROC))
  Ic = 0.0_dp; Istar = 0.0_dp; lambda_star = 0.0_dp; lambda_1 = 0.0

  allocate(Jold(Nlambda, atmos%Nspace), Jnew(Nlambda, atmos%Nspace))
  Jold = 0d0; Jnew(:,:) = 0.0_dp
  allocate(Sth(Nlambda, atmos%Nspace), Snew(Nlambda, atmos%Nspace), Sold(Nlambda, atmos%Nspace))
  Sth = 0.; Sold = 0.0; Snew = 0.0
  allocate(Kappa_tot(Nlambda, atmos%Nspace)); Kappa_tot = 0.
  allocate(beta(Nlambda, atmos%Nspace)); beta = 0.
    
!   sampling lambda for Jnu, at the end interpolated on the NLTE grid
!   allocate(lambda(Nlambda))
!   a1 = real(NLTEspec%lambda(1)) !max(lambda_min, real(NLTEspec%lambda(1)))
!   a2 = min(lambda_max0, real(NLTEspec%lambda(NLTEspec%Nwaves)))
!   a0 = 368.
!   a3 = 91.
!   evaluate opacities on the exact same grid, better than interpolate
!   lambda(1:20) = span(a1, a3-0.1, 20)
!   lambda(21:50) = span(a3, a3+2, 30)
!   lambda(51:51+249) = span(a3+2+0.1, a0, 250)
!   lambda(301:350) = spanl(a0+10, a2, 50)
!   lambda(1:20) = span(a1, a3-0.1, 20)
!   lambda(21:20+150) = span(a3, a3+2, 150)
!   lambda(171:170+300) = span(a3+2+0.1, a0, 300)
!   lambda(471:470+530) = spanl(a0+10, a2, 530)

	lambda = NLTEspec%lambda_cont
  
  Istar(:) = Bpnu (real(etoile(1)%T,kind=dp), lambda)

  write(*,*) "  -> interpolating contopac on Jnu grid for each frequency.."
  write(*,*) "       lambda (min, max in nm):", minval(lambda), maxval(lambda)
  do icell=1, atmos%Nspace
   if (atmos%icompute_atomRT(icell) > 0) then
    
    !LTE at start ?
    Jold(:,icell) = 0.0!bpnu(atmos%T(icell), lambda)
    
   					
   	CALL bezier2_interp(NLTEspec%Nwaves,NLTEspec%lambda,&
   					NLTEspec%AtomOpac%Kc(:,icell),Nlambda, lambda, kappa_tot(:,icell))
   	if (any_nan_infinity_vector(kappa_tot(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in kappa_tot after interpolation.."
   	 write(*,*) "icell=", icell, "kappa=",kappa_tot(:,icell)
   	 write(*,*) "Kc=", NLTEspec%AtomOpac%Kc(:,icell)
   	endif
   					
   					
   	CALL bezier2_interp(NLTEspec%Nwaves,NLTEspec%lambda,&
   					NLTEspec%AtomOpac%sca_c(:,icell),Nlambda, lambda, beta(:,icell))
   	beta(:,icell) = beta(:,icell) / kappa_tot(:,icell)
   	
   	if (any_nan_infinity_vector(beta(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in beta after interpolation.."
   	 write(*,*) "icell=", icell, "beta=",beta(:,icell)
   	 write(*,*) "sigma=", NLTEspec%AtomOpac%sca_c(:,icell)
   	endif   					

   	CALL bezier2_interp(NLTEspec%Nwaves,NLTEspec%lambda,&
   					NLTEspec%AtomOpac%jc(:,icell),Nlambda, lambda, Sth(:,icell))
   					
   	if (any_nan_infinity_vector(Sth(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in Sth after interpolation.."
   	 write(*,*) "icell=", icell, "eta_th=",Sth(:,icell)
   	 write(*,*) "eta=", NLTEspec%AtomOpac%jc(:,icell)
   	endif   
   	
   	Sth(:,icell) = Sth(:,icell) / ( kappa_tot(:,icell) * (1.-beta(:,icell)) + 1e-30)

   	if (any_nan_infinity_vector(Sth(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in Sth after interpolation.."
   	 write(*,*) "icell=", icell, "Sth=",Sth(:,icell)
   	 write(*,*) "kappa_abs=", ( kappa_tot(:,icell) * (1.-beta(:,icell)) + 1e-30)
   	endif     	
   					
    Sold(:,icell) = (1.-beta(:,icell))*Sth(:,icell) + beta(:,icell) * Jold(:,icell)
    
   	if (any_nan_infinity_vector(Sold(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in Sold after interpolation.."
   	 write(*,*) "icell=", icell, "Sold=",Sold(:,icell)
   	 write(*,*) "Jold=", Jold(:,icell)
   	endif  

   endif
  enddo
  write(*,*) "Sold (max,min)", maxval(Sold), minval(Sold,mask=Sold>0)
  write(*,*) "Sth (max,min)", maxval(STh), minval(Sth,mask=Sold>0)
  write(*,*) "Jold (max,min)", maxval(Jold), minval(Jold,mask=Sold>0)
  write(*,*) "beta (max,min)", maxval(beta), minval(beta,mask=Sold>0)
  write(*,*) "kappa (max,min)", maxval(kappa_tot), minval(kappa_tot,mask=Sold>0)
  write(*,*) "  ->..done"


  if (.not.allocated(lcell_converged)) allocate(lcell_converged(n_cells))
  lcell_converged(:) = .false.

  n_rayons_max = atmos%Nrays
  
  ng_acc = .false. !temp, because don't know if it works
  accelerated  = .false.
  if (ng_acc) then
   allocate(Jflat(atmos%Nspace*Nlambda))
   CALL initNg(atmos%Nspace*Nlambda, 6, 2, 3, NgJ)
  endif


  labs = .true.
  id = 1
  etape_start = 2
  etape_end = 2

if (write_convergence_file ) then
 open(unit=20, file="Jnu_convergence.s", status="unknown")
 write(20,*) maxIter, Nlambda, n_cells
endif
 
     do etape=etape_start, etape_end

      if (etape==1) then 
        lfixed_Rays = .true.
        n_rayons = 2
        iray_start=1
        lprevious_converged = .false.
		write(*,*) " Using ", n_rayons, " rays for Jnu."
      else if (etape==2) then 
  		lfixed_rays = .true.
  		n_rayons = n_rayons_max!min(n_rayons_max,n_rayons_start)
  		iray_start = 1
  		lprevious_converged = .false.
		write(*,*) " Using ", n_rayons, " rays for Jnu."
  	  else
  	    CALL ERROR("etape unkown")
  	  end if
  	  
if (write_convergence_file ) write(20,*) etape, n_rayons

  	  
  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0
		 				

        do while (.not.lconverged)

        	n_iter = n_iter + 1
        	if (n_iter > maxIter)  exit
            write(*,*) " -> Iteration #", n_iter, " Step #", etape

  			if (lfixed_rays) then
    			stream = 0.0
    			do i=1,nb_proc
     				stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    			end do
 			 end if

if (write_convergence_file ) write(20,*) " -> Iteration #", n_iter
            imax = 1
            imax_s = 1
            
 			!$omp parallel &
            !$omp default(none) &
            !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02)&
            !$omp private(argmt,norme, icell, l_unconverged) &
            !$omp shared(atmos,NLTEspec,lambda_star, Snew, Sold, Sth, lambda_1, Istar) &
            !$omp shared(lkeplerian,n_iter) &
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,lcell_converged) &
            !$omp shared(n_cells,ds, Jold, Jnew, beta, kappa_tot, Ic) &
            !$omp shared(lfixed_Rays,lnotfixed_Rays,labs,etape,pos_em_cellule)
            !$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !$ id = omp_get_thread_num() + 1
   				!!l_unconverged = (atmos%icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
   				!!if (l_unconverged) then 
      			l_unconverged = (.not.lcell_converged(icell))
   				if (atmos%icompute_atomRT(icell)>0) then
   				   Jnew(:,icell) = 0. !here otherwise if all cell converged -> Jnew = 0 and we never enter here
           		   Snew(:,icell) = 0. 
           		   lambda_1(:,id) = 0.        		   
           		   
					do iray=iray_start, iray_start-1+n_rayons
					
    					if (etape==1) then
        					! Position = milieu de la cellule
        					x0 = r_grid(icell)
        					y0 = 0.0_dp
        					z0 = z_grid(icell)
        					if (lkeplerian) then
        						write(*,*) "lkeplerian should be false here"
        						stop
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

      	 				else !etape 2
                   	    ! Position aleatoire dans la cellule
         					rand  = sprng(stream(id))
            				rand2 = sprng(stream(id))
            				rand3 = sprng(stream(id))

            				CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
            				rand = sprng(stream(id))
           					W0 = 2.0_dp * rand - 1.0_dp !nz
            				W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
            				SRW02 = sqrt(W02)
            				rand = sprng(stream(id))
            				ARGMT = PI * (2.0_dp * rand - 1.0_dp)
            				U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
            				V0 = SRW02 * sin(ARGMT) !ny = sin(theta) * sin(phi)
		 				end if !etape

						CALL INTEG_RAY_JNU(id, icell, x0, y0, z0, u0, v0, w0, iray, labs, &
											kappa_tot, Sold, Istar, lambda_star, Ic)
						!LI
		                Jnew(:,icell) = Jnew(:,icell) + Ic(:,id)/n_rayons
! 		                lambda_1(:,id) = lambda_1(:,id) + lambda_star(:,id) / n_rayons
           
		                !hogerheijde
! 		                Jnew(:,icell) = Jnew(:,icell) + ( Ic(:,id)*dexp(-ds(iray,id)*kappa_tot(:,icell)) +&
!  		                 (1.-dexp(-ds(iray,id)*kappa_tot(:,icell)))*Sold(:,icell) )/n_rayons	
!  		                 					
      			    enddo !iray

      		   endif !icompute_AtomRT
      		   
        	   Snew(:,icell) = Sth(:,icell) * (1.-beta(:,icell)) + Jnew(:,icell) * beta(:,icell)
        	   !!ALI, but lambda is not exactly the same as for hogerheijde
!         	   Snew(:,icell) = (Snew(:,icell)- beta(:,icell)*lambda_1(:,id) * Sold(:,icell)) / (1.-beta(:,icell)*lambda_1(:,id))
        	   
     		end do !icell
        	!$omp end do
        	!$omp end parallel

            !should be para
        	diff = 0d0
        	dSource = 0.0_dp
  			cell_loop2 : do icell=1, atmos%Nspace
  				if (atmos%icompute_atomRT(icell)>0) then
  					
  						dN = dabs(1.0 - maxval(Sold(:,icell) / Snew(:,icell)))
  						dSource = max(dSource, dN)

						dJ = 0.0_dp
						do la=1, Nlambda
						 if (Jnew(la, icell) > 0) then 
						  dJ = max(dJ,dabs(1.-Jold(la,icell)/Jnew(la,icell)))
						  imax = locate(dabs(1.-Jold(:,icell)/Jnew(:,icell)),dabs(1.-Jold(la,icell)/Jnew(la,icell)))
						 endif 
						enddo
						if (mod(icell,10)==0) write(*,*) icell, " ::> dJ(icell)", real(dJ)!," Jmax/min:", maxval(Jnew(:,icell)), minval(Jnew(:,icell))
if (write_convergence_file ) write(20,*) icell, " ::> dJ(icell)", real(dJ), " Jmax/min:", maxval(Jnew(:,icell)), minval(Jnew(:,icell))
						if (dJ > diff) then
						  diff = dJ
						  icell_max = icell
						  !write(*,*) icell_max, " ::> dJ(icell_max)", real(diff), maxval(Jnew(:,icell_max)), minval(Jnew(:,icell_max))
						endif
     					lcell_converged(icell) = (real(dJ) < precision)		
     			end if
     			Jold(:,icell) = Jnew(:,icell)
     			Sold(:,icell) = Snew(:,icell)
     		end do cell_loop2 !icell
     		
     		! include cells only that have converge for all frequencies, or for each cell / frequencies?
!      		Jold(:,:) = Jnew(:,:)
!      		Sold(:,:) = Snew(:,:)

     		write(*,*) " >>> dS = ", dSource
if (write_convergence_file ) write(20,*) " >>> dS = ", dSource
         	if (accelerated) then
         	  write(*,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff," (Accelerated)"
if (write_convergence_file ) write(20,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff," (Accelerated)"
         	else
         	  write(*,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff
if (write_convergence_file ) write(20,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff
         	endif
         	write(*,*) "   >>> ", "Jmax/min (icell_max):", maxval(Jnew(:,icell_max)), minval(Jnew(:,icell_max)), "T=",atmos%T(icell_max), "ne=", atmos%ne(icell_max)
         	write(*,*) "   >>> ", "Jmax/min (all):", maxval(Jnew(:,:)), minval(Jnew(:,:))
if (write_convergence_file ) write(20,*) "   >>> ", "Jmax/min (icell_max):", maxval(Jnew(:,icell_max)), minval(Jnew(:,icell_max))
if (write_convergence_file ) write(20,*) "   >>> ", "Jmax/min (all):", maxval(Jnew(:,:)), minval(Jnew(:,:))
	        write(*,*) " <-> # unconverged cells : ", size(pack(lcell_converged,mask=lcell_converged.eqv..false.)), &
	          100.*real(size(pack(lcell_converged,mask=lcell_converged.eqv..false.)))/real(n_cells), "%"
if (write_convergence_file ) write(20,*) " <-> # unconverged cells : ", size(pack(lcell_converged,mask=lcell_converged.eqv..false.)), &
	          100.*real(size(pack(lcell_converged,mask=lcell_converged.eqv..false.)))/real(n_cells), "%"
!-> gfortran does not like to compare logical with ==
! 	        write(*,*) " <-> # unconverged cells : ", size(pack(lcell_converged,mask=lcell_converged==.false.)), &
! 	          100.*real(size(pack(lcell_converged,mask=lcell_converged==.false.)))/real(n_cells), "%"
! if (write_convergence_file ) write(20,*) " <-> # unconverged cells : ", size(pack(lcell_converged,mask=lcell_converged==.false.)), &
! 	          100.*real(size(pack(lcell_converged,mask=lcell_converged==.false.)))/real(n_cells), "%"

        	!!lconverged = (real(diff) < precision)

        	
        	if (real(diff) < precision) then
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
              			if (n_iter >= maxIter) then
             		 		write(*,*) "Warning : not enough rays to converge !!"
                 			lconverged = .true.
              			end if
              		end if

          	   end if
        	end if

	    end do !while
        write(*,*) etape, "Threshold =", precision
	  end do !over etapes
if (write_convergence_file ) close(20)
	  
  if (.not.lstop_after_jnu) then
  
!
!      call freeSpectrum
!      !redfine spectrum grid here
!      if (ltab_wavelength_image) NLTEspec%write_wavelength_grid = .true.
!      CALL initSpectrum(vacuum_to_air=lvacuum_to_air)
!      if ((atmos%NactiveAtoms>0) .or. .not.(ltab_wavelength_image)) then
!       if (n_etoiles > 0) CALL init_stellar_disk !for all wavelengths, all stars at disk centre
!       write(*,*) " recomputing background opacities after Jnu..."
!       CALL alloc_atom_quantities
!       CALL compute_opacities
!       write(*,*) " ..done"
!      endif
!      CALL reallocate_mcfost_vars()

   write(*,*) " -> Interpolating Jnu on image grid.."
   do icell=1, atmos%Nspace
  	CALL bezier2_interp(Nlambda, lambda, Jnew(:,icell), NLTEspec%Nwaves, NLTEspec%lambda, NLTEspec%Jc(:,icell))
   enddo
   write(*,*) "Jmax/min after interpolation on the image grid"
   write(*,*) maxval(NLTEspec%Jc(:,:)), minval(NLTEspec%Jc(:,:))
   write(*,*) " -> ..done"
  endif
  
  
  open(unit=20, file="Jnu_no_interp.s", status="unknown")
  write(20,*) n_cells, Nlambda
  do icell=1, n_cells
  do la=1, Nlambda
    write(20,'(1F12.5,5E20.7E3)') lambda(la), Jnew(la,icell), Sth(la,icell), beta(la,icell), kappa_tot(la,icell), Sold(la,icell)
   enddo
  enddo
  close(20)


  if (ng_acc) CALL freeNg(NgJ)
  deallocate(Jold, Sth, kappa_tot,  beta, Jnew, Istar, Ic, lambda, Sold, lambda_star, Snew)
  if (allocated(Jflat)) deallocate(jflat)

 ! ------------------------------------------------------------------------------------ !
 RETURN
 END SUBROUTINE Iterate_Jnu
 

	SUBROUTINE INTEG_RAY_LINE_I_CNTRB(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	! Computes the contribution functions along to the Intensities using converged
	! populations
	! ------------------------------------------------------------------------------- !
		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs !used in NLTE but why?
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
		real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, chiI, LimbD, chiIc
		real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chil, etal
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars
		integer :: la


		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp
		tau_c(:) = 0.0_dp
		chiI(:) = 0.0_dp
		chiIc(:) = 0.0_dp


		NLTEspec%I(:,iray,id) = 0d0
		NLTEspec%Ic(:,iray,id) = 0d0

		!S_contrib(:,icell,id) = 0.0_dp

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

					call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)

					NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
       				NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau_c)
       				RETURN
      			end if
   			 endif

			nbr_cell = nbr_cell + 1


			previous_cell = 0 ! unused, just for Voronoi
			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)



			if (lcellule_non_vide) then
				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
     ! opacities in m^-1
				l_contrib = l_contrib * AU_to_m !l_contrib in m
				if ((nbr_cell == 1).and.labs) ds(iray,id) = l * AU_to_m


				CALL initAtomOpac(id)
				if (atmos%NpassiveAtoms>0) &
					CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
					
				chil(:) = NLTEspec%AtomOpac%chi_p(:,id)
				etal(:) = NLTEspec%AtomOpac%eta_p(:,id)

				chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell) + tiny_chi 
				chiI(:)  = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)

				if (atmos%NactiveAtoms>0) then 
					
					CALL initAtomOpac_nlte(id)
      
					CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)
					
					chil(:) = chil(:) + NLTEspec%ATomOpac%chi(:,id)
					etal(:) = etal(:) + NLTEspec%AtomOpac%eta(:,id)
					
					chiIc(:) = chiIc(:) + NLTEspec%AtomOpac%Kc_nlte(:,icell)

					chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id) +  NLTEspec%AtomOpac%Kc_nlte(:,icell)
				endif !active atoms

				dtau(:)   = l_contrib * chiI(:)
				dtau_c(:) = l_contrib * chiIc(:)
            
				Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
				Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
      
      			!Always include electron scattering if NLTE loop.
      			!But ATM, uses %Jc even for lines.
				if (atmos%Nactiveatoms>0) then 
					Snu = Snu + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell) ) / chiI(:)
					Snu_c = Snu_c + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + &
					NLTEspec%AtomOpac%jc_nlte(:,icell) )/ chiIc(:)

				else if ((atmos%electron_scattering).and.(atmos%NactiveAtoms==0)) then
					Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
					Snu_c = Snu_c + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) / chiIc(:)
				endif


! 				NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
! 				NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c
				do la=1, NLTEspec%Nwaves
				
					!!Seems that chil * Ic << etal / chiI even for absorption lines	
					!S_contrib(la,icell,id) =  chil(la) * ( NLTEspec%Ic(la,iray,id) - etal(la) / chiI(la) ) * dexp(-tau(la))
										
					NLTEspec%I(la,iray,id) = NLTEspec%I(la,iray,id) + dexp(-tau(la)) * (1.0_dp - dexp(-dtau(la))) * Snu(la)
					NLTEspec%Ic(la,iray,id) = NLTEspec%Ic(la,iray,id) + dexp(-tau_c(la)) * (1.0_dp - dexp(-dtau_c(la))) * Snu_c(la) 
					
					S_contrib(la,icell,id) =  (etal(la) / chiI(la)) * dexp(-tau(la))
					!S_contrib(la,icell,id) =  etal(la) * dexp(-tau(la))			
			
					!->Flow chart
					!S_contrib(la,icell,id) = ( etal(la)/chiI(la) ) * (1.0_dp - dexp(-dtau(la))) * dexp(-tau(la))

				enddo
     

     
				tau = tau + dtau
				tau_c = tau_c + dtau_c

			end if  ! lcellule_non_vide
		end do infinie
	RETURN
	END SUBROUTINE INTEG_RAY_LINE_I_CNTRB
  

END MODULE AtomicTransfer
