! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atoms
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

module atom_transfer

	use atom_type, only			: AtomType
	use opacity, only			: dealloc_atom_quantities, alloc_atom_quantities, compute_atom_quantities, compute_background_continua, &
									interp_background_opacity, opacity_atom_loc, interp_contopac, compute_nlte_bound_free, &
									nlte_bound_free, background_continua
	use background_opacity, only: Thomson
	use Planck, only			: bpnu
	use Profiles, only			: Iprofile, Iprofile_thomson, Iprofile_cmf_to_obs, Zprofile
	use spectrum_type, only     : dk, sca_c, chi, eta, chi_c, eta_c, eta_c_nlte, chi_c_nlte, eta0_bb, chi0_bb, lambda, Nlambda, lambda_cont, Nlambda_cont, Itot, Icont, Istar_tot, Istar_cont, &
									Stokes_Q, Stokes_U, Stokes_V, Flux, Fluxc, F_QUV, Ksi, init_spectrum, init_spectrum_image, dealloc_spectrum, Jnu_cont, eta_es, alloc_flux_image, &
									dealloc_jnu, write_image, write_flux_ascii, reallocate_rays_arrays
	use atmos_type, only		: icompute_atomRT, lmagnetized, ds, Nactiveatoms, Atoms, calc_ne, Natom, ne, T, &
									readatmos_ascii, dealloc_atomic_atmos, ActiveAtoms, nHmin, hydrogen, helium
	use readatom, only			: readAtomicModels
	use lte, only				: set_LTE_populations, nH_minus
	use constant, only			: MICRON_TO_NM, hc_k
	use solvene, only			: solve_electron_density
	use getlambda, only			: hv
	use voigtfunctions, only	: Voigt, VoigtHumlicek
	use profiles, only			: Iprofile, Profile
	use broad, only 			: damping
	use io_atomic_pops, only	: write_pops, write_electron, write_hydrogen, write_Hminus
	use io_opacity, only 		: write_Jnu, write_taur, write_contrib_lambda_ascii, read_jnu_ascii, Jnu_File_ascii, &
									write_collision_matrix_atom, write_collision_matrix_atom_ascii, &
									write_radiative_rates_atom, write_rate_matrix_atom, write_cont_opac_ascii
	use math
	!$ use omp_lib
	use molecular_emission, only: v_proj
	use input, only				: lkeplerian, linfall, RT_line_method, nb_proc, RT_line_method, limb_darkening, mu_limb_darkening
	use parametres, only		: Rmax, Rmin, map_size, zoom, n_cells, prt_solution, lcontrib_function, lelectron_scattering, n_rad, nz, n_az, distance, ang_disque, &
									l_sym_ima, etoile, npix_x, npix_y, npix_x_save, npix_y_save, lpluto_file, lmodel_ascii, density_file, lsolve_for_ne, ltab_wavelength_image, &
									lvacuum_to_air, n_etoiles, lread_jnu_atom, lstop_after_jnu, llimb_darkening, dpops_max_error, laccurate_integ, NRAYS_ATOM_TRANSFER, &
									DPOPS_SUB_MAX_ERROR, n_iterate_ne,lforce_lte, loutput_rates

	use grid, only				: test_exit_grid, cross_cell, pos_em_cellule, move_to_grid
	use dust_transfer, only		: compute_stars_map
	use dust_ray_tracing, only	: RT_n_incl, RT_n_az, init_directions_ray_tracing,tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az,tab_RT_incl, stars_map, kappa
	use stars, only				: intersect_spots, intersect_stars
	!use wavelengths, only		:
	use mcfost_env, only		: dp
	use constantes, only		: tiny_dp, huge_dp, au_to_m, pc_to_au, deg_to_rad, tiny_real, pi, deux_pi, rad_to_deg
	use utils, only				: rotation_3d, cross_product
	use naleat, only 			: seed, stream, gtype
	use cylindrical_grid, only	: r_grid, z_grid
	use messages, only 			: error, warning
	
	use statequil_atoms, only   : invpop_file, profiles_file, unit_invfile, unit_profiles, prec_pops, calc_bb_rates, calc_bf_rates, calc_rate_matrix, update_populations, fill_collision_matrix, &
									init_bb_rates_atom, initgamma, initgamma_atom , init_rates_atom, store_radiative_rates, calc_rates, store_rate_matrices
	use collision, only			: CollisionRate !future deprecation
	use impact, only			: collision_hydrogen

	implicit none

 !OpenMp debug
 	integer, dimension(:), allocatable :: npix_x_id, npix_y_id, cell_id
 !RT procedure
	procedure(integ_ray_line_i), pointer :: integ_ray_line => NULL()
 !Temporary variable for Zeeman calculations
	real(kind=dp), dimension(:,:), allocatable :: QUV
 !Temporary variables for Contribution functions
	real(kind=dp), allocatable :: S_contrib(:,:,:), S_contrib2(:,:), mean_formation_depth
 !NLTE
	logical, allocatable :: lcell_converged(:)
	real(kind=dp), dimension(:,:,:), allocatable :: gpop_old, Tex_old, local_pops
	real(kind=dp), allocatable :: Gammaij_all(:,:,:,:), Rij_all(:,:,:), Rji_all(:,:,:) !need a flag to output it
	integer :: NmaxLevel, NmaxTr

 	contains

  	subroutine integ_ray_line_i(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	!
	! ------------------------------------------------------------------------------- !

		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
		real(kind=dp), dimension(Nlambda) :: Snu, LD, tau, dtau
		real(kind=dp), dimension(Nlambda_cont) :: Snu_c, LDc, dtau_c, tau_c
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp
		tau_c(:) = 0.0_dp

		Itot(:,iray,id) = 0d0
		Icont(:,id) = 0d0
		

!!write(*,*) "id=",id
!!write(*,*) "****"
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
				lcellule_non_vide = (icompute_atomRT(icell) > 0)
				if (icompute_atomRT(icell) < 0) return !-1 if dark
			else
				lcellule_non_vide=.false.
			endif
    
    ! Test sortie ! "The ray has reach the end of the grid"

			if (test_exit_grid(icell, x0, y0, z0)) return

			if (lintersect_stars) then
				if (icell == icell_star) then
					!can be done better
					call calc_stellar_surface_brightness(Nlambda_cont,lambda_cont,i_star, x0, y0, z0, u,v,w,LDc)
       				Icont(:,id) =  Icont(:,id) + LDc(:) * Istar_cont(:,i_star)*exp(-tau_c(:))
					call calc_stellar_surface_brightness(Nlambda,lambda,i_star, x0, y0, z0, u,v,w,LD)
					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * LD(:)
       				return
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

				if ((nbr_cell == 1).and.labs) then
					ds(iray,id) = l * AU_to_m
				endif
! 				dk(iray,id) = nint(1e-3/hv * v_proj(icell,(x0+x1)*0.5,(y0+y1)*0.5,(z0+z1)*0.5,u,v,w))
! 				dk(iray,id) = nint(1e-3/hv * 0.5 * (v_proj(icell,x0,y0,z0,u,v,w)+v_proj(icell,x1,y1,z1,u,v,w)))

     									
				!total bound-bound + bound-free + background opacities lte + nlte
				chi(:,id) = chi0_bb(:,icell)
				eta(:,id) = eta0_bb(:,icell)
				!includes a loop over all bound-bound, passive and active
				call opacity_atom_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,( (nbr_cell==1).and.labs ) ) 

				dtau(:)   = l_contrib * chi(:,id)
				            
				if (lelectron_scattering) then
					Snu = ( eta(:,id) + eta_es(:,icell) ) / ( chi(:,id) + tiny_dp )
					Snu_c = eta_c(:,icell) + sca_c(:,icell) * Jnu_cont(:,icell)
				else
					Snu = eta(:,id) / ( chi(:,id) + tiny_dp )
					Snu_c = eta_c(:,icell)
				endif
				
				if (Nactiveatoms > 0) then
					Snu_c = ( Snu_c + eta_c_nlte(:,icell) ) / ( chi_c(:,icell) + chi_c_nlte(:,icell) + tiny_dp)
					dtau_c(:) = l_contrib * (chi_c(:,icell) + chi_c_nlte(:,icell))
				else
					Snu_c = Snu_c / (chi_c(:,icell) + tiny_dp)
					dtau_c(:) = l_contrib * chi_c(:,icell)
				endif


				Itot(:,iray,id) = Itot(:,iray,id) + exp(-tau(:)) * (1.0_dp - exp(-dtau(:))) * Snu(:)
				Icont(:,id) = Icont(:,id) + exp(-tau_c(:)) * (1.0_dp - exp(-dtau_c(:))) * Snu_c(:)

				tau = tau + dtau
				tau_c = tau_c + dtau_c

			end if  ! lcellule_non_vide
		end do infinie

	return
	end subroutine integ_ray_line_i


 !rebuilding
 subroutine integ_ray_line_z(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! This routine is "LTE", although NLTE pops can be used
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
		real(kind=dp), dimension(Nlambda) :: Snu, LD, tau, dtau
		real(kind=dp), dimension(Nlambda_cont) :: Snu_c, LDc, dtau_c, tau_c
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  real(kind=dp) :: facteur_tau
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp
  tau_c = 0.0_dp


  Itot(:,iray,id) = 0d0
  Icont(:,id) = 0d0
  Stokes_V(:,id) = 0d0
  Stokes_Q(:,id) = 0d0
  Stokes_U(:,id) = 0d0

  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  infinie : do
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     lcellule_non_vide = (icompute_atomRT(icell)>0)
    else
     lcellule_non_vide=.false.
    endif


    if (test_exit_grid(icell, x0, y0, z0)) return
    
			if (lintersect_stars) then
				if (icell == icell_star) then
					!can be done better
					call calc_stellar_surface_brightness(Nlambda_cont,lambda_cont,i_star, x0, y0, z0, u,v,w,LDc)
       				Icont(:,id) =  Icont(:,id) + LDc(:) * Istar_cont(:,i_star)*exp(-tau_c(:))
					call calc_stellar_surface_brightness(Nlambda,lambda,i_star, x0, y0, z0, u,v,w,LD)
					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * LD(:)
       				return
      			end if
   			 endif

    nbr_cell = nbr_cell + 1

    previous_cell = 0

    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)



    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs)
     l_contrib = l_contrib * AU_to_m
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m
	 dk(iray,id) = nint(1e-3/hv * v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w))

    !.....

    !Correct line source fonction from polarisation
    !explicit product of Seff = S - (K/chiI - 1) * I
    !Is there a particular initialization to do for polarisation ?
    !because it will be zero at the first place
!      Snu(:) = Snu(:) -NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI *  NLTEspec%Stokes_Q(:,iray,id) - &
!           NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI *  NLTEspec%Stokes_U(:,iray,id) - &
!           NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI *  NLTEspec%Stokes_V(:,iray,id)
! 
!      NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu
! 
!      NLTEspec%Stokes_Q(:,iray,id) = NLTEspec%Stokes_Q(:,iray,id) + &
!                              exp(-tau) * (1.0_dp - exp(-dtau)) * (&
! 			NLTEspec%AtomOpac%etaQUV_p(:,1,id) /chiI - &
! 			NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI * NLTEspec%I(:,iray,id) - &
! 			NLTEspec%AtomOpac%rho_p(:,3,id)/chiI * NLTEspec%Stokes_U(:,iray,id) + &
! 			NLTEspec%AtomOpac%rho_p(:,2,id)/chiI * NLTEspec%Stokes_V(:,iray, id) &
!      ) !end Snu_Q
!      NLTEspec%Stokes_U(:,iray,id) = NLTEspec%Stokes_U(:,iray,id) + &
!                              exp(-tau) * (1.0_dp - exp(-dtau)) * (&
! 			NLTEspec%AtomOpac%etaQUV_p(:,2,id) /chiI - &
! 			NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI * NLTEspec%I(:,iray,id) + &
! 			NLTEspec%AtomOpac%rho_p(:,3,id)/chiI * NLTEspec%Stokes_Q(:,iray,id) - &
! 			NLTEspec%AtomOpac%rho_p(:,1,id)/chiI * NLTEspec%Stokes_V(:,iray, id) &
!      ) !end Snu_U
!      NLTEspec%Stokes_V(:,iray,id) = NLTEspec%Stokes_V(:,iray,id) + &
!                              exp(-tau) * (1.0_dp - exp(-dtau)) * (&
! 			NLTEspec%AtomOpac%etaQUV_p(:,3,id) /chiI - &
! 			NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI * NLTEspec%I(:,iray,id) - &
! 			NLTEspec%AtomOpac%rho_p(:,2,id)/chiI * NLTEspec%Stokes_Q(:,iray,id) + &
! 			NLTEspec%AtomOpac%rho_p(:,1,id)/chiI * NLTEspec%Stokes_U(:,iray, id) &
!      ) !end Snu_V

!      facteur_tau = 1.0
! 
!      tau = tau + dtau * facteur_tau
!      tau_c = tau_c + dtau_c

    end if
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  return
  end subroutine integ_ray_line_z

  subroutine flux_pixel_line(id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! -------------------------------------------------------------- !
   ! Computes the flux emerging out of a pixel.
   ! see: mol_transfer.f90/intensite_pixel_mol()
  ! -------------------------------------------------------------- !

   integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
   real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
   real(kind=dp), intent(in) :: pixelsize,u,v,w
   integer, parameter :: maxSubPixels = 32
   real(kind=dp) :: x0,y0,z0,u0,v0,w0
   real(kind=dp), dimension(Nlambda) :: Iold, I0
   real(kind=dp), dimension(Nlambda_cont) :: I0c
   !real(kind=dp), dimension(NLTEspec%Nwaves_cont) :: I0c !-> removed at the moment
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
     if (lmagnetized.and. PRT_SOLUTION == "FULL_STOKES") QUV(:,:) = 0d0 !move outside
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
           call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
             call integ_ray_line(id, icell, x0,y0,z0,u0,v0,w0,iray,labs)

             I0 = I0 + Itot(:,iray,id) / npix2
             I0c = I0c + Icont(:,id) / npix2
             
             if (lmagnetized.and.PRT_SOLUTION == "FULL_STOKES") then
             	QUV(3,:) = QUV(3,:) + Stokes_V(:,id)
             	QUV(1,:) = QUV(1,:) + Stokes_Q(:,id)
             	QUV(2,:) = QUV(2,:) + Stokes_U(:,id)
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
        !if (diff > 1) write(*,*) 'pixel not converged:', iter, subpixels, diff
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
    Flux(:,1,1,ibin,iaz) = Flux(:,1,1,ibin,iaz) + I0 * normF
    Fluxc(:,1,1,ibin,iaz) = Fluxc(:,1,1,ibin,iaz) + I0c * normF
    if (lmagnetized.and.PRT_SOLUTION == "FULL_STOKES") then
     F_QUV(1,:,1,1,ibin,iaz) = F_QUV(1,:,1,1,ibin,iaz)+QUV(1,:) * normF !U
     F_QUV(2,:,1,1,ibin,iaz) = F_QUV(2,:,1,1,ibin,iaz)+QUV(2,:) * normF !Q
     F_QUV(3,:,1,1,ibin,iaz) = F_QUV(3,:,1,1,ibin,iaz)+QUV(3,:) * normF !V
    end if
  else
    Flux(:,ipix,jpix,ibin,iaz) = I0 * normF
    Fluxc(:,ipix,jpix,ibin,iaz) = I0c * normF
    if (lmagnetized.and.PRT_SOLUTION == "FULL_STOKES") then
     F_QUV(1,:,ipix,jpix,ibin,iaz) = QUV(1,:) * normF !U
     F_QUV(2,:,ipix,jpix,ibin,iaz) = QUV(2,:) * normF !Q
     F_QUV(3,:,ipix,jpix,ibin,iaz) = QUV(3,:) * normF !V
    end if
  end if
  if (lcontrib_function) Ksi(:,:,ibin,iaz) = Ksi(:,:,ibin,iaz) + S_contrib2(:,:)


  return
  end subroutine flux_pixel_line

 subroutine emission_line_map(ibin,iaz)
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

  integer, parameter :: n_rad_RT = 1000, n_phi_RT = 1500 !(100, 36)
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp):: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda!, lM, lR, lB, ll
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
    rmax_RT = Rmax * 2.0_dp

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
    id = 1 ! pour code sequentiel

    if (l_sym_ima) then
      cst_phi = pi  / real(n_phi_RT,kind=dp)
    else
      cst_phi = deux_pi  / real(n_phi_RT,kind=dp)
    endif

     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1
        r = tab_r(ri_RT)
		!!write(*,*) "r = ", r / etoile(1)%r
        taille_pix =  fact_A * r ! racine carree de l'aire du pixel

        do phi_RT=1,n_phi_RT ! de 0 a pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)
           !!write(*,*) " phi = ", phi * rad_to_deg
           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image
            ! C'est le centre en fait car dx = dy = 0.
           call flux_pixel_line(id,ibin,iaz,n_iter_min,n_iter_max, i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
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
!!write(*,*) "*****"
!!write(*,*) "i=",i, "id=",id
        do j = 1,npix_y

!!write(*,*) "j=",j
        !do j=npix_y/2+1,npix_y/2+1
        !write(*,*) "ipix, jpix", i, j
           ! Coin en bas gauche du pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           call flux_pixel_line(id,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        end do !j
     end do !i
     !$omp end do
     !$omp end parallel
  end if

 return
 end subroutine emission_line_map

 subroutine atom_line_transfer()
 ! --------------------------------------------------------------------------- !
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
 ! --------------------------------------------------------------------------- !
  integer :: atomunit = 1, nact
  integer :: icell, m
  integer :: ibin, iaz
  integer :: n_rayons_start
  !!integer, parameter :: n_rayons_start = 100
  integer, parameter :: n_rayons_start2 = 100, maxIter = 11
  integer :: n_rayons_max = n_rayons_start2 * (2**(maxIter-1))
  integer, parameter :: Nrayone = 1
  character(len=20) :: ne_start_sol = "H_IONISATION"
  character(len=20)  :: newPRT_SOLUTION = "FULL_STOKES"
  logical :: lread_jnu_ascii = .false.
  type (AtomType), pointer :: atom
  integer :: alloc_status
  
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


!   if (PRT_SOLUTION == "FULL_STOKES") then
!    call Warning(" Full Stokes solution not allowed. Stokes polarization not handled in SEE yet.")
!    PRT_SOLUTION = "FIELD_FREE"
!   end if

!   if (lmagnetized .and. PRT_SOLUTION /= "NO_STOKES") then
!    if (PRT_SOLUTION == "FIELD_FREE") newPRT_SOLUTION = "FULL_STOKES"
!    call adjustStokesMode(PRT_SOLUTION)
!   end if

 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
        !! ----------------------- Read Model ---------------------- !!

  if (.not.lpluto_file) then
   if (lmodel_ascii) then
    call readAtmos_ascii(density_file)
   end if
  end if
        !! --------------------------------------------------------- !!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------- READATOM and INITIZALIZE POPS ------------------------------ !

  call readAtomicModels(atomunit)

  !compute first guess of electron density ??
  if (.not.calc_ne) calc_ne = lsolve_for_ne
  if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
  if (calc_ne) then
   if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
   if (.not.lsolve_for_ne) write(*,*) "Solving for electron density"
   write(*,*) " Starting solution : ", ne_start_sol
   if ((ne_start_sol == "NE_MODEL") .and. (calc_ne)) then
    write(*,*) "WARNING, ne from model is presumably 0 (or not given)."
    write(*,*) "Changing ne starting solution NE_MODEL in H_IONISATION"
    ne_start_sol = "H_IONISATION"
   end if
   call Solve_Electron_Density(ne_start_sol)
   call write_Electron()
  end if

  	call set_LTE_populations !write pops at the end because we probably have NLTE pops also

	if (NactiveAtoms > 0) then

		n_rayons_start = Nrays_atom_transfer

		if (.not.laccurate_integ) n_rayons_max = n_rayons_start
		
		call init_Spectrum(n_rayons_max,lam0=500d0,vacuum_to_air=lvacuum_to_air)
		if (n_etoiles > 0) call init_stellar_disk
		write(*,*) " Computing background opacities..."
		call alloc_atom_quantities !call compute_atom_quantities(icell) for each cell
		call compute_background_continua
		write(*,*) " ..done"
!-> done during nlte iterations even the interpolation as chi_c_nlte changes
! 		write(*,*) " init bound-free opacities..."
! 		call compute_nlte_bound_free
! 		write(*,*) " ..done"
! 		call interp_contopac

		!use the same rays as nlteloop
		if (lelectron_scattering) call iterate_Jnu
		
		if (allocated(ds)) deallocate(ds)
		allocate(ds(n_Rayons_max, nb_proc), stat=alloc_status)
		if (alloc_status > 0) call error ("Allocation error ds")
		!allocated(dk ....)
		
		NmaxLevel = 0
		NmaxTr = 0
		do nact=1,Nactiveatoms
			atom => ActiveAtoms(nact)%ptr_atom
			allocate(atom%Gamma(atom%Nlevel,atom%Nlevel,nb_proc),stat=alloc_status)
			atom%Gamma(:,:,:) = 0.0_dp
			if (alloc_status > 0) call error("Allocation error atom%Gamma")
			
			allocate(atom%C(atom%Nlevel,atom%Nlevel,nb_proc),stat=alloc_status)!,n_cells)) 
			atom%C(:,:,:) = 0.0_dp
			if (alloc_status > 0) then
				write(*,*) nact, atom%ID
				call error("Allocation error atom%C")
			endif

			atom%NLTEpops = .true.
			NmaxLevel = max(NmaxLevel, atom%Nlevel)
			NmaxTr = max(NmaxTr, atom%Ncont + atom%Nline)

			atom => NULL()
		enddo
		!remove some arrays in depth to save memory, like vbroad per atom not needed, it is fast to compute locally
   
		allocate(lcell_converged(n_cells),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error lcell_converged")
		lcell_converged = .false.
		allocate(gpop_old(NactiveAtoms,Nmaxlevel,n_cells),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error gpop_old")
		gpop_old = 0.0_dp
		allocate(Tex_old(NactiveAtoms,NmaxTr,n_cells),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error Tex_old")
		Tex_old = 0.0_dp
		!sub iter needed ?
		!allocate(pops(NactiveAtoms, NmaxLevel, nb_proc)); pops = 0d0
		
		!remove some arrays in depth to save memory, like vbroad per atom not needed, it is fast to compute locally
		if (loutput_rates) then
			allocate(Rij_all(NactiveAtoms,NmaxTr,n_cells),Rji_all(NactiveAtoms,NmaxTr,n_cells),stat=alloc_status)
			if (alloc_status > 0) call error("Allocation error Rij/Rji _all")
			Rij_all(:,:,:) = 0.0_dp; Rji_all(:,:,:) = 0.0_dp
			allocate(Gammaij_all(NactiveAtoms,Nmaxlevel,Nmaxlevel,n_cells),stat=alloc_status)
			if (alloc_status > 0) call error("Allocation error Gammij_all")
			Gammaij_all(:,:,:,:) = 0.0_dp
		endif

		write(*,*) " ** Number max of levels among all atoms:", Nmaxlevel
		write(*,*) " ** Number max of transitions among all atoms:", NmaxTr


		if (NactiveAtoms > 1) then
			write(*,*) "   -> Solving for kinetic equations for ", Nactiveatoms, " atoms"
		else
			write(*,*) "   -> Solving for kinetic equations for ", Nactiveatoms, " atom"
		endif
		write(*,*) " Max error : ", dpops_max_error, dpops_sub_max_error
				
		call NLTEloop(n_rayons_max, n_rayons_start, n_rayons_start2, maxIter,.true.)

		
		deallocate(lcell_converged, gpop_old, Tex_old)
		if (allocated(dk)) deallocate(dk)
		deallocate(ds)
		
		!dealloc some quantity not needed anymore
		!and write some output
		do nact=1,Natom
			atom => Atoms(nact)%ptr_atom
			if (atom%active) then
				if (allocated(atom%Gamma)) deallocate(atom%Gamma)
				if (allocated(atom%C)) deallocate(atom%C)
				if (loutput_rates) then
					call write_radiative_rates_atom(atom, Rij_all(nact,1:atom%Ntr,:), Rji_all(nact,1:atom%Ntr,:))
					call write_rate_matrix_atom(atom, Gammaij_all(nact,1:atom%Nlevel,1:atom%Nlevel,:))
				endif
			endif
			if (associated(hydrogen,atom)) call write_collision_matrix_atom(hydrogen)
			do m=1, atom%Nline
				deallocate(atom%lines(m)%phi_loc)
			enddo
			atom => NULL()
		enddo
		if (loutput_rates) deallocate(Rij_all,Rji_all,Gammaij_all)
		
		call write_Jnu
		!recompute some opacities ?
		if (lelectron_scattering) call iterate_Jnu
		!-> write_Jnu and iterate_Jnu ??


		!->> Case of ltab or not and add chi_c_nlte and eta_c_nlte
		!then if ltab_wavelength_image ...
		if (ltab_wavelength_image) then
!  I need to set up grid image without having defined the full image in the first place !
!    	call initSpectrumImage() !deallocate waves arrays/ define a new grid also resample Jnu if any
!    if (n_etoiles > 0) call init_stellar_disk 
!    write(*,*) " Computing continuum opacities for image..."
!    if (atmos%NactiveAtoms >0) call dealloc_atom_quantities
!    call alloc_atom_quantities
!    call compute_opacities !recompute background opac
!    if (atmos%NactiveAtoms > 0) then
!     call compute_nlte_bound_free
!     !! can be added to Kc and jc to save memory for image
!     !NLTEspec%AtomOpac%Kc(:,:,1) = NLTEspec%AtomOpac%Kc(:,:,1) + NLTEspec%AtomOpac%Kc_nlte(:,:)
!     !NLTEspec%AtomOpac%jc(:,:,1) = NLTEspec%AtomOpac%jc(:,:,1) + NLTEspec%AtomOpac%jc_nlte(:,:)
! 	!deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)
!    endif
!    write(*,*) " ..done"
!    !add NLTE contiuna
!    !!call reallocate_mcfost_vars() !wavelength only?
!    call reallocate_mcfost_wavelength_arrays()
		else
			!just remove the Nrays dependencies
			call reallocate_rays_arrays(nrayOne)
		endif

	else !no nlte or using old populations so no chi_c_nlte no eta_c_nlte and atom is passive with old_populations
  
  		if (ltab_wavelength_image) then !atomic lines transfer with either lte or nlte populations
  									!using user defined grid
  									
  			!and electron scattering here ?
  	
  		else
			call init_Spectrum(Nrayone,lam0=500d0,vacuum_to_air=lvacuum_to_air)
			if (n_etoiles > 0) call init_stellar_disk
			write(*,*) " Computing background opacities..."
			call alloc_atom_quantities !call compute_atom_quantities(icell) for each cell
			call compute_background_continua
			write(*,*) " ..done" 				
			if (lelectron_scattering) then
			
!    call alloc_jnu(.false.) -> already allocated in alloc_spectrum
!    if (.not.lread_jnu_atom) then
				if (lread_jnu_atom) then
					call read_jnu_ascii
					lread_jnu_ascii = .true.
! 		call read_jnu() !I now write Jnu in ascii file to be interpolated if lread_jnu_atom
!	   the interpolated version of Jnu is written at the end in write_jnu()
				else
					call iterate_Jnu()
					call write_Jnu
					if (lstop_after_jnu) then
						write(*,*) " Jnu calculation done." 
						stop
					endif
				endif !lelectron scattering	
			
			endif !electron scatt
  		endif !atomic lines transfer with either lte or nlte populations
  		  !using the same grid as the nlte loop
    	call interp_contopac
		
	endif !NLTE

	!for all so after nlte case (debug)
	do m=1,Natom
		if (Atoms(m)%ptr_atom%initial_solution /= "OLD_POPULATIONS") then
			call write_pops(Atoms(m)%ptr_atom)
		endif
	end do

! 	if (lmagnetized .and. PRT_SOLUTION /= "NO_STOKES") then
! 		if (PRT_SOLUTION == "FIELD_FREE") PRT_SOLUTION = newPRT_SOLUTION
! 		call adjustStokesMode(PRT_SOLUTION)
! 		allocate(QUV(3, Nlambda))
! 		QUV = 0d0
! 	end if

	call alloc_flux_image()
	write(*,*) "Computing emission flux map..."

! 	if (lcontrib_function) then
! 		write(*,*) "   >>> including contribution function Ksi..."
! 		allocate(S_contrib(Nlambda,n_cells,nb_proc))
! 		allocate(S_contrib2(Nlambda,n_cells))
! 		!allocate(mean_formation_depth(npix_x, npix_y))
! 		S_contrib(:,:,:) = 0d0; S_contrib2(:,:) = 0d0
! 		INTEG_RAY_LINE => NULL()
! 		INTEG_RAY_LINE => INTEG_RAY_LINE_I_CNTRB
! 	end if

  !Actual flux calculation
	call init_directions_ray_tracing()
	do ibin=1,RT_n_incl
		do iaz=1,RT_n_az
			call emission_line_map(ibin,iaz)
		end do
	end do
	write(*,*) " ..done"

	write(*,*) "Writing result to file..."
	call write_flux_ascii
	
	if (3*size(Flux)/(1024.**3) <= 6.) call write_image
	deallocate(Flux, Fluxc) !free here to release a bit of memory
! 	if (lcontrib_function) then
! 		call WRITE_CNTRB_FUNC_pix()
! 		deallocate(S_contrib, S_contrib2, Ksi)
! 	end if
    
    !ascii files can be large !
    if (loutput_rates) then
		call write_cont_opac_ascii(hydrogen)
 		!chi0_bb and eta0_bb not modified
		call write_contrib_lambda_ascii(654.0_dp,0.0_dp,0.0_dp,1.0_dp,.true.)
		
		call write_taur(500._dp,0._dp,0._dp,1.0_dp)
	endif

	if (lelectron_scattering.and.(NactiveAtoms==0)) then !otherwise if active, written even if not included in emissivity
		!Jnu is written to ascii file if read
		if (lread_jnu_atom) then
			if (lread_jnu_ascii) call write_Jnu
		endif
		call dealloc_Jnu()
	endif
	
! 	if (allocated(npix_x_id)) then
! 		open (unit=3,file="npixxid.txt",status="unknown")
! 		write(3,*) size(npix_x_id), nb_proc
! 		do ibin=1,size(npix_x_id)
! 			write(3,*) ibin, npix_x_id(ibin)
! 		enddo
! 		write(3,*) size(npix_y_id)
! 		do ibin=1,size(npix_y_id)
! 			write(3,*) ibin, npix_y_id(ibin)
! 		enddo
! 		close(3)
! 		deallocate(npix_x_id, npix_y_id)
! 	endif

 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------- CLEANING ------------------------------------------ !
	if (allocated(QUV)) deallocate(QUV)

	call dealloc_Spectrum !deallocate spectral variables
	call dealloc_atomic_atmos
	NULLIFY(Profile, Voigt, INTEG_RAY_LINE)
	!call dealloc_profile_line_profiles

	return
	end subroutine atom_line_transfer
! ------------------------------------------------------------------------------------ !
	subroutine NLTEloop(n_rayons_max,n_rayons_1,n_rayons_2,maxIter,verbose)
	! -------------------------------------------------------- !
	! Descriptor here
	! -------------------------------------------------------- !
#include "sprng_f.h"
		integer, intent(in) :: n_rayons_1, n_rayons_2, maxIter, n_rayons_max
		logical, intent(in) :: verbose
		integer, parameter :: max_sub_iter = 10000, max_global_iter = 1000
		integer :: etape, etape_start, etape_end, iray, n_rayons
		integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, la, kr
		integer, dimension(nb_proc) :: max_n_iter_loc
		logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged
		real :: rand, rand2, rand3, precision, precision_sub
		real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt
		real(kind=dp) :: diff, norme, dT_max, dN, dN1, dJ, lambda_max
		real(kind=dp) :: dT, dN2, dN3, diff_old
		real(kind=dp), allocatable :: dTM(:), dM(:), TM(:), Tion_ref(:), Tex_ref(:)
		real(kind=dp), allocatable :: Jnew(:,:), Jold(:,:), Jnew_cont(:,:)
		logical :: labs, iterate_ne = .false.
		logical :: l_iterate
		integer :: nact, imax, icell_max, icell_max_2
		integer :: icell, ilevel
		character(len=20) :: ne_start_sol = "NE_MODEL"
! 		real(kind=dp), dimension(3, atmos%Nrays, nb_proc) :: xyz0, uvw0
		type (AtomType), pointer :: atom

		open(unit=unit_invfile, file=trim(invpop_file), status="unknown")
		write(unit_invfile,*) n_cells
		open(unit=unit_profiles, file=trim(profiles_file), status="unknown")	
  
		allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
		allocate(dTM(Nactiveatoms)); dM=0d0 !keep tracks of Tex for all cells for each atom
		allocate(Jnew(Nlambda, n_cells)); Jnew = 0.0
		allocate(Jold(Nlambda, n_cells)); Jold = 0.0
		allocate(Jnew_cont(Nlambda_cont, n_cells)); Jnew_cont = 0.0
		allocate(Tex_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tex for all cells for each line of each atom
		allocate(Tion_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tion for all cells for each cont of each atom
		diff_old = 0.0_dp

		!!n_rayons_max = atmos%Nrays
! 		xyz0(:,:,:) = 0d0
! 		uvw0(:,:,:) = 0d0

		labs = .true. !to have ds at cell icell = eval_operator
		id = 1
		etape_start = 2
		etape_end = 2
		if (laccurate_integ) then
			write(*,*) " Using step 3 with ", n_rayons_max, " nrays max"
			etape_end = 3
		else
			etape_end = 2
		endif
		if (etape_start==1) then
			call error("step 1 not implemented")
		endif


		iterate_ne = (n_iterate_ne>0)
		if (iterate_ne) then
			write(*,*) " before iterate ne I need te recompute gij for continua !!"
		endif

		do etape=etape_start, etape_end

			if (etape==1) then
				stop

			else if (etape==2) then 
				lfixed_rays = .true.
				n_rayons = n_rayons_1
				iray_start = 1
				lprevious_converged = .false.
				lcell_converged(:) = .false.
				precision = dpops_max_error
				precision_sub = dpops_sub_max_error
			else if (etape==3) then
		  		lfixed_rays = .false.
  				n_rayons = n_rayons_2
  				iray_start = 1
  				lprevious_converged = .false.
				lcell_converged(:) = .false.
				precision = 1e-1 !dpops_max_error
				precision_sub = dpops_sub_max_error
			else
				call ERROR("etape unkown")
			end if
  	  
			write(*,*)  "-> Using ", n_rayons, ' rays for step', etape

			lnotfixed_rays = .not.lfixed_rays
			lconverged = .false.
			n_iter = 0

			do while (.not.lconverged)

				n_iter = n_iter + 1
				!if (n_iter > max_global_iter) exit !change step !maxIter
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
				!$omp shared(gpop_old,Tex_old, lcell_converged, chi0_bb, eta0_bb, lnotfixed_Rays, lfixed_rays) &
				!$omp shared (n_cells,icompute_atomRT, NactiveAtoms, l_iterate, lelectron_scattering, eta_es, Jnu_cont) &
				!$omp shared(ActiveAtoms)
				!$omp do schedule(static,1)
				do icell=1, n_cells
					!$ id = omp_get_thread_num() + 1
! 					if (lfixed_rays) then
   						l_iterate = (icompute_atomRT(icell)>0)
! 					else
!    						l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! 					endif							
   					if (l_iterate) then
						do nact=1, NactiveAtoms
							atom => ActiveAtoms(nact)%ptr_atom
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
						!but also update din sub it ? so we repeat calculation here
						call NLTE_bound_free(icell)
						!chi_c, sca_c and eta_c unchanged except if ne has changed during iteration and updated
						!at the end, only nlte changed 
						call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
					end if
				end do
				!$omp end do
				!$omp end parallel
        	
				!$omp parallel &
				!$omp default(none) &
				!$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dM, dN, dN1,dT_max)&
				!$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, nact, atom, l_iterate) &
				!$omp shared(icompute_atomRT, dpops_sub_max_error, verbose,lkeplerian,n_iter,precision_sub) &
				!$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,lcell_converged,loutput_rates) &
				!$omp shared(n_cells, gpop_old,lforce_lte, integ_ray_line, Itot, Icont, Jnu_cont, eta_es) &
				!$omp shared(Jnew, Jnew_cont, Jold, lelectron_scattering,chi0_bb, etA0_bb) &
				!$omp shared(nHmin, chi_c, chi_c_nlte, eta_c, eta_c_nlte, ds, Rij_all, Rji_all, Nmaxtr, Gammaij_all, Nmaxlevel) &
				!$omp shared(lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape,pos_em_cellule,Nactiveatoms)
				!$omp do schedule(static,1)
				do icell=1, n_cells
					!$ id = omp_get_thread_num() + 1
! 					if (lfixed_rays) then
   						l_iterate = (icompute_atomRT(icell)>0)
! 					else
!    						l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! 					endif				
   					if (l_iterate) then
   						Jnew(:,icell) = 0.0
   						Jnew_cont(:,icell) = 0.0
					!!if (atmos%icompute_atomRT(icell)>0) then

						call fill_Collision_matrix(id, icell) !computes = C(i,j) before computing rates
						call initGamma(id) !init Gamma to C and init radiative rates
                    
						do iray=iray_start, iray_start-1+n_rayons
						
							if ( (etape==2).or.(etape==3) ) then
                   	    ! Position aleatoire dans la cellule
								rand  = sprng(stream(id))
								rand2 = sprng(stream(id))
								rand3 = sprng(stream(id))

								call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

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

							call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, iray, labs)
							Jnew(:,icell) = Jnew(:,icell) + Itot(:,iray,id) / n_rayons
							Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,id) / n_rayons

						enddo !iray
						!!if (lelectron_scattering) then
						!!always allocated if nlte ?, need to change that but remember we need
						!! in the nlte loop, to write the full mean intensity to fits file.
							eta_es(:,icell) = Jnew(:,icell) * thomson(icell)
						!!endif
		
						if (loutput_Rates) &
							call store_radiative_rates(id, icell, n_rayons, Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell), Jnew(:,icell))

						if (lforce_lte) then
							call update_populations(id, icell, diff, .false., n_iter)
							lconverged_loc = .true.
						else
							lconverged_loc = .false.
							do iray=1, n_rayons
								call calc_rates(id, icell, iray, n_rayons)
							enddo
							call calc_rate_matrix(id, icell, lforce_lte)
! 							call update_populations(id, icell, diff, .false., n_iter)
! 							lconverged_loc = .true.
						endif

						n_iter_loc = 0
     				!!only iterate on cells which are not converged
						do while (.not.lconverged_loc)
							n_iter_loc = n_iter_loc + 1

! 							if ((mod(n_iter_loc,max(1,max_sub_iter/10))==0)) then
! 								write(*,*) " -- Global iteration #", n_iter, " step #", etape, "eps: ", dpops_sub_max_error
! 							endif!(mod(n_iter_loc,max(1,max_sub_iter/10))==0)
							call update_Populations(id, icell, diff,.false.,n_iter_loc)
! 							if ((mod(n_iter_loc,max(1,max_sub_iter/10))==0)) write(*,*) " --------------------------------------------------------------- "

							!recompute with new populations
							!bound-bound are computed on the fly during the propagation or bb rates
							!so they always use the new populations
							call NLTE_bound_free(icell)
! 							nHmin(icell) = nH_minus(icell)
! 							call background_continua(icell)
							call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
							
                        !for this icell (id) and local iteration
							if (diff < precision_sub) then
								lconverged_loc = .true.
     					    !write(*,*) "converged", id, icell, "dpops(sub) = ", diff
							else

								call initGamma(id)
								!local value of the continuum is computed at lambda0 to evaluate bb rates
								do iray=1,n_rayons
									call calc_rates(id, icell, iray, n_rayons)
								enddo
								
								call calc_rate_matrix(id, icell, lforce_lte)

							end if
! 							if (n_iter_loc >= max_sub_iter) then 
! 
! 								lconverged_loc = .true.
! 								write(*,*) " step #", etape, "eps: ", precision_sub, " niter_loc:", n_iter_loc
! 								write(*,*) "id=",id, " icell=",icell, " diff=", diff
! 
! 							end if

						end do !local sub iteration
						if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
						if (loutput_Rates) &
							call store_rate_matrices(id,icell,Nmaxlevel,Gammaij_all(:,:,:,icell))
	
					end if !icompute_atomRT
				end do !icell
				!$omp end do
				!$omp end parallel

     		!Global convergence Tests

            !should be para
				dM(:) = 0.0
				diff = 0.0
				dJ = 0.0
				dT = 0.0
				dTM(:) = 0.0
				Tex_ref(:) = 0.0
				Tion_ref(:) = 0.0
   				icell_max = 1
   				icell_max_2 = 1
				cell_loop2 : do icell=1,n_cells
! 					if (lfixed_rays) then
   						l_iterate = (icompute_atomRT(icell)>0)
! 					else
!    						l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! 					endif				
   					if (l_iterate) then
						dN = 0.0 !for all levels of all atoms of this cell
						dN2 = 0.0
						do nact=1,NactiveAtoms
							atom => ActiveAtoms(nact)%ptr_atom
     					         					    
							do ilevel=1,atom%Nlevel
							
								!if (atom%n(ilevel,icell) > 0) then
									!dN1 = abs(1d0-gpop_old(nact,ilevel,icell)/atom%n(ilevel,icell))
     				    		dN1 = abs(atom%n(ilevel,icell)-gpop_old(nact,ilevel,icell))/(gpop_old(nact,ilevel,icell)+tiny_dp)
								dN = max(dN1, dN) !compare with all atoms and levels
     				    						  !for this cell
								dM(nact) = max(dM(nact), dN1) !compare for one atom
     				    		 							   !for all cells
								!endif
								!!write(*,*) "pops", dN1,gpop_old(nact,ilevel,icell),atom%n(ilevel,icell)
							end do !over ilevel
							do kr=1, atom%Nline
							
								!dN3 = abs(1.0 - Tex_old(nact, kr, icell) / ( atom%lines(kr)%Tex(icell) + tiny_dp ) )
								dN3 = abs(atom%lines(kr)%Tex(icell) - Tex_old(nact, kr, icell)) / ( Tex_old(nact, kr, icell) + tiny_dp ) 
								dN2 = max(dN3, dN2)
								!!write(*,*) "line", icell, T(icell), Tex_old(nact, kr, icell),atom%lines(kr)%Tex(icell), dN2, dN3
								!dTM(nact) = max(dTM(nact), dN3)
								if (dN3 >= dTM(nact)) then
									dTM(nact) = dN3
									Tex_ref(nact) = atom%lines(kr)%Tex(icell)
									icell_max = icell
								endif
							enddo
							
							do kr=1, atom%Ncont
								!dN3 = abs(1.0 - Tex_old(nact, kr+atom%Nline, icell) / ( atom%continua(kr)%Tex(icell) + tiny_dp ))
								dN3 = abs(atom%continua(kr)%Tex(icell) - Tex_old(nact, kr+atom%Nline, icell)) / ( Tex_old(nact, kr+atom%Nline, icell) + tiny_dp )
								dN2 = max(dN3, dN2)
								!!write(*,*) "cont", kr, " icell=",icell, " T=", T(icell), " Told=",Tex_old(nact, kr+atom%Nline, icell)," Tex=",atom%continua(kr)%Tex(icell), dN2, dN3
								!dTM(nact) = max(dTM(nact), dN3)
								if (dN3 >= dTM(nact)) then
									dTM(nact) = dN3
									Tion_ref(nact) = atom%continua(kr)%Tex(icell)
									icell_max_2 = icell
								endif
							enddo
							atom => NULL()
						end do !over atoms
						!compare for all atoms and all cells
						diff = max(diff, dN) ! pops
						!diff = max(diff, dN2) ! Tex
						
						dN1 = abs(1d0 - maxval(Jold(:,icell)/(tiny_dp + Jnew(:,icell))))
						!dN1 = maxval(abs(Jold(:,icell)-Jnew(:,icell))/(1d-30 + Jold(:,icell)))
						if (dN1 > dJ) then
							dJ = dN1
							imax = locate(abs(1d0 - Jold(:,icell)/(tiny_dp + Jnew(:,icell))), dJ)
							!imax = locate(abs(Jold(:,icell)-Jnew(:,icell))/(tiny_dp + Jold(:,icell)),dJ)
							lambda_max = lambda(imax)
						endif
						!replace old value by new
						Jold(:,icell) = Jnew(:,icell)
						Jnu_cont(:,icell) = Jnew_cont(:,icell)
! 						if (lelectron_scattering) then
						!diff = max(diff,dJ)
! 						endif
						lcell_converged(icell) = (real(diff) < precision) !(real(diff) < dpops_max_error)
					end if
				end do cell_loop2 !icell

				write(*,'(" -> "(1I5)" sub-iterations")') maxval(max_n_iter_loc)
				write(*,'(" -> icell_max1 #"(1I5)," icell_max2 #"(1I5))') icell_max, icell_max_2
				write(*,'(" -> dJ="(1ES14.5E3)" @"(1F14.4)" nm")') dJ, lambda_max !at the end of the loop over n_cells
				write(*,*) " ------------------------------------------------ "
				do nact=1,NactiveAtoms
					write(*,'("             Atom "(1A2))') ActiveAtoms(nact)%ptr_atom%ID
					write(*,'("   >>> dpop="(1ES14.5E3))') dM(nact)
					write(*,'("   >>>   dT="(1ES14.5E3))') dTM(nact)
					write(*,'("   >>> Te(icell_max2)="(1F14.4)" K", " Tion="(1F14.4)" K")') T(icell_max_2), Tion_ref(nact)
					write(*,'("   >>> Te(icell_max1)="(1F14.4)" K", " Texi="(1F14.4)" K")') T(icell_max), Tex_ref(nact)
					!if (Tion_ref(nact) /= 0.0)
					!if (Tex_ref(nact) /= 0.0)
					write(*,*) " ------------------------------------------------ "
				enddo
				write(*,'(" <<->> diff="(1ES14.5E3)," old="(1ES14.5E3))') diff, diff_old !at the end of the loop over n_cells
	        	write(*,"('Unconverged cells #'(1I5), ' fraction :'(1F12.3)' %')") size(pack(lcell_converged,mask=lcell_converged.eqv..false.)), 100.*real(size(pack(lcell_converged,mask=lcell_converged.eqv..false.)))/real(n_cells)
				write(*,*) " *************************************************************** "
				diff_old = diff
				!Not used if the next is not commented out
! 				lconverged = (real(diff) < dpops_max_error) !global convergence for all iterations


				if (real(diff) < precision) then
           			if (lprevious_converged) then
            	  		lconverged = .true.
           			else
            	  		lprevious_converged = .true.
          	    	endif
        		else !continue to iterate even if n_rayons max is reached ?
           			lprevious_converged = .false.
           			if (.not.lfixed_rays) then
              			n_rayons = n_rayons * 2
              			write(*,*) ' -- Increasing number of rays :', n_rayons
             			if (n_rayons > n_rayons_max) then
              				if (n_iter >= maxIter) then
             		 			write(*,*) "Warning : not enough rays to converge !!"
                 				lconverged = .true.
              				end if
              			end if
          	   		end if
        		end if
        									
			!Not optimal, but fast relatively
			!must be parra
				if (iterate_ne .and. (mod(n_iter,n_iterate_ne)==0))  then
				write(*,*) "Recompute lte opacities"
! 					write(*,*) " Solve ne Global iteration:", n_iter
! 					write(*,*) " --> old max/min ne", maxval(ne), minval(ne,mask=ne>0)
! 					do icell=1, n_cells
! 						if (icompute_atomRT(icell) > 0) then
! 							sca_c(:,icell) = sca_c(:,icell) - sigma_e * ne(icell)
! 							do nact=1,NactiveAtoms
! 								atom => ActiveAtoms(nact)%ptr_atom
! 								do kr=1, atom%Ncont
! 									atom%continua(kr)%gij(:,icell) = atom%continua(kr)%gij(:,icell) / ne(icell)
! 								enddo
! 								atom => Null()
! 							end do
! 						endif
! 					enddo
! 					call Solve_Electron_Density(ne_start_sol)
! 					do icell=1, n_cells
! 						if (icompute_atomRT(icell) > 0) then
! 							sca_c(:,icell) =sca_c(:,icell) + sigma_e * atmos%ne(icell)
! 
! 							do nact=1,NactiveAtoms
! 								atom => ActiveAtoms(nact)%ptr_atom
! 								do kr=1, atom%Ncont
! 									atom%continua(kr)%gij(:,icell) = atom%continua(kr)%gij(:,icell) * ne(icell)
! 								enddo
! 								atom => Null()
! 							end do
							!avoid to call this one because otherwise might be some problem with ni-nj_gij as it uses nlte pops
							!and we do not need to updated only either adamp or phi, and eventually gij or sigma(:,icell)
							!call compute_atom_quantities(icell) ? update, gij with new nstar and ne
							!also update adamp if needed ect
							!call background_continua(icell)
! 						endif
! 					enddo
	
				end if


			end do !while
			write(*,*) "step: ", etape, "Threshold: ", precision!dpops_max_error

		end do !over etapes
		
! 				do icell=1, n_cells
!    					if (icompute_atomRT(icell)>0) then
! 						call NLTE_bound_free(icell)
! 						call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
! 					end if
! 				end do


! 		do icell=1,n_cells
! 			if (icompute_atomRT(icell) > 0) then
! 				nHmin(icell) = nH_minus(icell)
! 				do nact=1,Natom
! 					do kr=1, Atoms(nact)%ptr_atom%Nline
! 						call Damping(icell, Atoms(nact)%ptr_atom, kr, dN)
! 						if (allocated(Atoms(nact)%ptr_atom%lines(kr)%a)) Atoms(nact)%ptr_atom%lines(kr)%a(icell) = dN
! 					enddo
! 				enddo
! 				call background_continua (icell)
! 			endif
! 		enddo

		!I recompute for all even for active atoms, for consistency. 
		!Because, only gij for active continua and ne are updated, not %nstar
! 		if (n_iterate_ne < 0) then
! 			write(*,*) "end LOOP: old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
! 			call SolveElectronDensity(ne_start_sol)
! 			!Recompute for all atoms the LTE pops
! 			write(*,*) "   recompute LTE pops for all atoms"
! 			do nact=1,atmos%NAtom
! 				atom => atmos%Atoms(nact)%ptr_atom
! 				if (atom%ID=="H") then
! 					call LTEpops_H()
! 				else
! 					call LTEpops(atom,.true.)
! 				endif
! 				atom => Null()
! 			end do
! 		endif

! 		if (iterate_ne .or. (n_iterate_ne < 0)) then
! 			call writeElectron(.true.) !the .true. means append, to compare with initial solution.
! 			!Recompute Background  + NLTE bound free.
! 			!Check it is not done after
! 			write(*,*) "Check the consistency of NLTE b-f and background opac with iterate_ne"
! 		endif
! -------------------------------- CLEANING ------------------------------------------ !

		!Always written in nlte
		!if (lelectron_scattering) then
  		open(unit=20, file=Jnu_File_ascii, status="unknown")
  		write(20,*) n_cells, Nlambda_cont
  		do icell=1, n_cells
  			do la=1, Nlambda_cont
    			write(20,'(1F12.5,5E20.7E3)') lambda_cont(la), Jnu_cont(la,icell), 0.0, sca_c(la,icell)/(chi_c(la,icell)+chi_c_nlte(la,icell)), 0.0, 0.0
   			enddo
  		enddo
  		close(20)
  		!endif

		deallocate(dM, dTM, Tex_ref, Tion_ref)
		if (allocated(Jnew)) deallocate(Jnew)
		if (allocated(Jnew_cont)) deallocate(Jnew_cont)
		if (allocated(Jold)) deallocate(Jold)
	
		close(unit=unit_invfile)
		close(unit=unit_profiles)
	return
	end subroutine NLTEloop
! ------------------------------------------------------------------------------------ !

 !building
!  subroutine adjustStokesMode(Solution)
!  !
!  !
!    character(len=*), intent(in) :: Solution
! 
! 
! write(*,*) "Note polarized profile not allocated yet, check that in profile and line%phiZ, line%psi"
! write(*,*) "Magnetic profile not ready yet"
! stop
! 
!     if (SOLUTION=="FIELD_FREE") then
!       !if (associated(Profile)) Profile => NULL()
!       !Profile => Iprofile !already init to that
!       if (allocated(NLTEspec%AtomOpac%rho_p)) deallocate(NLTEspec%AtomOpac%rho_p)
!       if (allocated(NLTEspec%AtomOpac%etaQUV_p)) deallocate(NLTEspec%AtomOpac%etaQUV_p)
!       if (allocated(NLTEspec%AtomOpac%chiQUV_p)) deallocate(NLTEspec%AtomOpac%chiQUV_p)
!       if (allocated(NLTEspec%F_QUV)) deallocate(NLTEspec%F_QUV)
!       if (allocated(NLTEspec%StokesV)) deallocate(NLTEspec%StokesV, NLTEspec%StokesQ, &
!       												NLTEspec%StokesU)
!       write(*,*) " Using FIELD_FREE solution for the SEE!"
!       call Warning("  Set PRT_SOLUTION to FULL_STOKES for images")
! 
!     else if (SOLUTION=="POLARIZATION_FREE") then
!      call ERROR("POLARIZATION_FREE solution not implemented yet")
!      !NLTE not implemented in integ_ray_line_z and Zprofile always compute the full
!      !propagation dispersion matrix, but only the first row is needed
!      if (associated(Profile)) Profile => NULL()
!      Profile => Zprofile
!      if (associated(INTEG_RAY_LINE)) INTEG_RAY_LINE => NULL()
!      INTEG_RAY_LINE => INTEG_RAY_LINE_Z
! 
!      !if (associated(Voigt)) Voigt => NULL
!      !Voigt => VoigtHumlicek !Already the algortihm used
! 
!      !beware, only I! polarization is neglected hence not allocated
!       if (.not.allocated(NLTEspec%AtomOpac%rho_p)) then !same for all, dangerous.
!       	allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
!         allocate(NLTEspec%AtomOpac%etaQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
!         allocate(NLTEspec%AtomOpac%chiQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
! 
!         write(*,*) " Using POLARIZATION_FREE solution for the SEE!"
!      end if
! 
!     else if (SOLUTION=="FULL_STOKES") then !Solution unchanged here
!       if (associated(Profile)) Profile => NULL()
!       Profile => Zprofile
!       if (associated(INTEG_RAY_LINE)) INTEG_RAY_LINE => NULL()
!       INTEG_RAY_LINE => INTEG_RAY_LINE_Z
!      !if (associated(Voigt)) Voigt => NULL
!      !Voigt => VoigtHumlicek !Already the algortihm used
!       !allocate space for Zeeman polarisation
!       !only LTE for now
!        if (.not.allocated(NLTEspec%StokesQ)) & !all allocated, but dangerous ?
!             allocate(NLTEspec%StokesQ(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC), &
!              NLTEspec%StokesU(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC), &
!              NLTEspec%StokesV(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC))!, &
!              !!NLTEspec%S_QUV(3,NLTEspec%Nwaves))
!       NLTEspec%StokesQ=0d0
!       NLTEspec%StokesU=0d0
!       NLTEspec%StokesV=0d0
!       if (.not.allocated(NLTEspec%F_QUV)) &
!       	allocate(NLTEspec%F_QUV(3,NLTEspec%Nwaves,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ))
!       NLTEspec%F_QUV = 0d0
!       if (.not.allocated(NLTEspec%AtomOpac%rho_p)) then !same for all, dangerous.
!       	allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
!         allocate(NLTEspec%AtomOpac%etaQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
!         allocate(NLTEspec%AtomOpac%chiQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
! 
!         write(*,*) " Using FULL_STOKES solution for the SEE!"
!       end if
!     else
!      call ERROR("Error in adjust StokesMode with solution=",Solution)
!     end if
!  return
!  end subroutine adjustStokesMode
 ! ------------------------------------------------------------------------------------ !
!should be done more properly and I should use functions from mcfostfost for stellar flux
!  subroutine reallocate_mcfost_vars()
!   !--> should move them to init_atomic_atmos ? or elsewhere
!   !need to be deallocated at the end of molecule RT or its okey ?`
!   integer :: la
!   call init_directions_ray_tracing()
!   if (.not.allocated(stars_map)) allocate(stars_map(npix_x,npix_y,3))
!   if (.not.allocated(stars_map_cont)) allocate(stars_map_cont(npix_x,npix_y,3))
!   n_lambda = NLTEspec%Nwaves
!   if (allocated(tab_lambda)) deallocate(tab_lambda)
!   allocate(tab_lambda(n_lambda))
!   if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
!   allocate(tab_delta_lambda(n_lambda))
!   tab_lambda = NLTEspec%lambda / MICRON_TO_NM
!   tab_delta_lambda(1) = 0d0! + tiny_dp !! I assume that Nwaves > 1 here
!   !!if (NLTEspec%Nwaves==1) tab_delta_lambda(1) = tab_delta_lambda(1) + tiny_dp
! 
!   do la=2,NLTEspec%Nwaves
!    tab_delta_lambda(la) = tab_lambda(la) - tab_lambda(la-1)
!   end do
! 
!   if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
!   allocate(tab_lambda_inf(n_lambda))
!   if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
!   allocate(tab_lambda_sup(n_lambda))
!   tab_lambda_inf = tab_lambda
!   tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
!   ! computes stellar flux at the new wavelength points
!   call deallocate_stellar_spectra()
!   ! probably do not deallocate, 'cause maybe dust is in here too.
!   if (allocated(kappa)) deallocate(kappa) !do not use it anymore
!   !used for star map ray-tracing.
!   call allocate_stellar_spectra(n_lambda)
!   call repartition_energie_etoiles()
!   !If Vfield is used in atomic line transfer, with lkep or linfall
!   !atmos%Vxyz is not used and lmagnetoaccr is supposed to be zero
!   !and Vfield allocated and filled in the model definition if not previously allocated
!   !nor filled.
!   if (allocated(Vfield).and.lkeplerian.or.linfall) then
!     write(*,*) "Vfield max/min", maxval(Vfield), minval(Vfield)
!   else if (allocated(Vfield).and.lmagnetoaccr) then
!    write(*,*) "deallocating Vfield for atomic lines when lmagnetoaccr = .true."
!    deallocate(Vfield)
!   end if
! 
!  return
!  end subroutine reallocate_mcfost_vars

!  subroutine reallocate_mcfost_wavelength_arrays()
!   !--> should move them to init_atomic_atmos ? or elsewhere
!   !need to be deallocated at the end of molecule RT or its okey ?`
!   integer :: la
!   n_lambda = NLTEspec%Nwaves
!   if (allocated(tab_lambda)) deallocate(tab_lambda)
!   allocate(tab_lambda(n_lambda))
!   if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
!   allocate(tab_delta_lambda(n_lambda))
!   tab_lambda = NLTEspec%lambda / MICRON_TO_NM
!   tab_delta_lambda(1) = 0d0
! 
!   do la=2,NLTEspec%Nwaves
!    tab_delta_lambda(la) = tab_lambda(la) - tab_lambda(la-1)
!   end do
!   !unlikely in NLTE because of the different atomic grids
!   !but if an image at only one wavelength --> tab_Delta_lambda=0
!   !!if I do not put a tiny_dp here, if only one wavelength delta_wl is 0
!   !which creates an overflow in stars.f90, because spectre and spectre0 are 0
!   if (size(tab_lambda)==1) tab_delta_lambda = tab_delta_lambda + tiny_dp
! 
!   if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
!   allocate(tab_lambda_inf(n_lambda))
!   if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
!   allocate(tab_lambda_sup(n_lambda))
!   tab_lambda_inf = tab_lambda
!   tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
!   ! computes stellar flux at the new wavelength points
!   call deallocate_stellar_spectra()
!   call allocate_stellar_spectra(n_lambda)
!   call repartition_energie_etoiles()
! 
!  return
!  end subroutine reallocate_mcfost_wavelength_arrays
 
 !futur move to initial_solution
 subroutine init_stellar_disk
  integer :: i_star!, lam

   write(*,*) " Computing Istar(mu=1) for each star..."
   !move elsewhere if we read a file, but normally should well fit with MCFOST routines
   do i_star=1, n_etoiles
    if (etoile(i_star)%T <= 1e-6) then
     call warning("Setting stellar radiation to 0: T* < 1e-6 K")
     Istar_tot(:,i_star) = 0.0_dp
     Istar_cont(:,i_star) = 0.0_dp
    else
     write(*,"( 'T(i_star='(1I1)') = '(1F14.7)' K' )") i_star, etoile(i_star)%T
     Istar_tot(:,i_star) = Bpnu(etoile(i_star)%T*1d0, lambda)
     Istar_cont(:,i_star) = Bpnu(etoile(i_star)%T*1d0, lambda_cont)
    endif
   enddo
   write(*,*) " ..done"
 
 return
 end subroutine init_stellar_disk

 
 subroutine calc_stellar_surface_brightness(N,lambda,i_star,x,y,z,u,v,w,gamma)
 ! ---------------------------------------------------------------!
  ! Compute the stellar radiation field surface brightness.
  ! Istar = B(x,y,z;mu) * Stellar spectrum or B * BlackBody
  ! return gamma, the brightness. For uniform disk, gamma is 1.
  !
  ! For a point on the stellar disk, returns gamma * limbdarkenig
 ! -------------------------------------------------------------- !
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
   mu = abs(x*u + y*v + z*w)/sqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
   if (real(mu)>1d0) then !to avoid perecision error
    write(*,*) "mu=",mu, x, y, z, u, v, w
    call Error(" mu limb > 1!")
   end if
   
   !1) Compute stellar flux from mcfost 
   ! .... 
   
   !2) Correct with the contrast gamma of a hotter/cooler region if any
   call intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect_spot)
   if (lintersect_spot) then
     gamma(:) = (exp(hc_k/lambda/real(etoile(i_star)%T,kind=dp))-1)/&
     			(exp(hc_k/lambda/etoile(i_star)%SurfB(ns)%T)-1)
     !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
   end if

   !3) Apply Limb darkening
   if (llimb_darkening) then
     call ERROR("option for reading limb darkening not implemented")
   else
     !write(*,*) maxval(uLD(real(etoile(i_star)%T,kind=dp))), minval(uLD(real(etoile(i_star)%T,kind=dp)))
     ulimb = 0.0 ! could use BB slope
     LimbDarkening = 1d0 - ulimb*(1d0-mu)
   end if
   !Istar(:) = energie(:) * LimbDarkening * gamma(:)
   gamma(:) = LimbDarkening * gamma(:)

 return
 end subroutine calc_stellar_surface_brightness

   subroutine INTEG_RAY_JNU(id,icell_in,x,y,z,u,v,w,iray,labs, kappa_tot, Snu, Istar, Ic, tau_tot)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray to
  ! compute coherent Jnu in the continuum
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(in), dimension(:,:) :: kappa_tot, Snu
  real(kind=dp), intent(in), dimension(:) :: Istar
  real(kind=dp), intent(out) :: Ic(:,:), tau_tot(:,:)
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
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
  tau_tot(:,id) = 0.0
  
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
     lcellule_non_vide = (icompute_atomRT(icell) > 0)
     if (icompute_atomRT(icell) < 0) return !-1 if dark
    else
     lcellule_non_vide=.false.
    endif
    
    !if (minval(tau_c) > 50.) return
    
    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) return

    if (lintersect_stars) then
      if (icell == icell_star) then
       !call calc_stellar_surface_brightness(size(Ic(:,1)),lambda,i_star,x0,y0,z0,u,v,w,LimbD)
       Ic(:,id) =  Ic(:,id) + Istar(:) * exp(-tau_c)
       return
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

      Ic(:,id) = Ic(:,id) + exp(-tau_c) * (1.0_dp - exp(-l_contrib * kappa_tot(:,icell))) * Snu(:,icell)

     if ((nbr_cell == 1).and.labs) then 
      ds(iray,id) = l * AU_to_m
     endif

     !tau = tau + dtau
     tau_c = tau_c + l_contrib * kappa_tot(:,icell)
     if ((nbr_cell) > 1) tau_tot(:,id) = tau_tot(:,id) + l_contrib * kappa_tot(:,icell)

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  return
  end subroutine INTEG_RAY_JNU
  
 subroutine Iterate_Jnu()
 ! -------------------------------------------------------- !
  ! Compute the mean radiation field at all cells
  ! evaluated on a small grid and then interpolated on the
  ! wavelength grid
 ! -------------------------------------------------------- !
#include "sprng_f.h"

  integer, parameter :: maxIter = 11, n_rayons_start2 = 1000
  integer :: n_rayons_start
  integer, parameter :: n_rayons_max = n_rayons_start2 * (2**(maxIter-1))
  real :: precision = 1e-1!, parameter
  real, parameter :: lambda_min = 5., lambda_max0 = 100000
  integer :: etape, etape_start, etape_end, iray, n_rayons, n_rayons_old
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lprevious_converged, write_convergence_file

  real :: rand, rand2, rand3, a1, a0, a2, a3!, fac_etape
  real(kind=dp) :: dSource, dN
  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dJ, diffs
  real(kind=dp), allocatable :: Jold(:,:), Jnew(:,:), tau_tot(:,:), Jnew_l(:,:), Jold_l(:,:)
  real(kind=dp), allocatable :: Snew(:,:), Kappa_tot(:,:), beta(:,:), Istar(:), Ic(:,:)
  real(kind=dp), allocatable :: lambda_star(:,:), Sth(:,:), Sold(:,:), Sline(:,:), delta_J(:)
                                   
  logical :: labs, l_iterate
  integer :: la, icell, imax, icell_max, icell_max_s, imax_s
  integer :: imu, Nmu, Nphi, iphi
  real(kind=dp) :: lambda_max, phia, mua, dphi, dmu
  
  n_rayons_start = Nrays_atom_transfer

  write(*,*) "   --> Lambda iterating Jnu with Nlambda ", Nlambda_cont
  precision = dpops_max_error

  write(*,*) "   precision in J is ", precision
  write(*,*) "   n_rayons_max is ", n_rayons_max
  write_convergence_file = .false.
  
  if (n_rayons_max <= 0) call error("Error n_rayons_max pb !")
  if (allocated(ds)) deallocate(ds)
  allocate(ds(n_rayons_max, nb_proc))
  ds = 0.0_dp !meters

  !only one star
  allocate(Istar(Nlambda_cont), Ic(Nlambda_cont, nb_proc), lambda_star(Nlambda_cont,nb_proc), delta_J(Nlambda_cont))
  Ic = 0.0_dp; Istar = 0.0_dp; lambda_star = 0.0_dp; delta_J = 0.0
  allocate(tau_tot(Nlambda_cont, nb_proc)); tau_tot = 0.0

  allocate(Jold(Nlambda_cont, n_cells))!, Jnew(Nlambda_cont, n_cells))
  Jold = 0d0!; Jnew(:,:) = 0.0_dp
  allocate(Sth(Nlambda_cont, n_cells), Snew(Nlambda_cont, n_cells), Sold(Nlambda_cont, n_cells))
  Sth = 0.; Sold = 0.0; Snew = 0.0
  allocate(Kappa_tot(Nlambda_cont, n_cells)); Kappa_tot = 0.
  allocate(beta(Nlambda_cont, n_cells)); beta = 0.
    
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


  
  Istar(:) = Istar_cont(:,1)
  Jold = Jnu_cont
  
  write(*,*) "  -> interpolating contopac on Jnu grid for each frequency.."
  write(*,*) "       lambda (min, max in nm):", minval(lambda_cont), maxval(lambda_cont)
  do icell=1, n_cells
   if (icompute_atomRT(icell) > 0) then
    
   					
   	kappa_tot(:,icell) = chi_c(:,icell)


   	if (any_nan_infinity_vector(kappa_tot(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in kappa_tot after interpolation.."
   	 write(*,*) "icell=", icell, "kappa=",kappa_tot(:,icell)
   	 write(*,*) "Kc=", chi_c(:,icell)
   	endif
   					
   					
   	beta(:,icell) = sca_c(:,icell) / ( kappa_tot(:,icell) + tiny_dp)

   	
   	if (any_nan_infinity_vector(beta(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in beta after interpolation.."
   	 write(*,*) "icell=", icell, "beta=",beta(:,icell)
   	 write(*,*) "sigma=", sca_c(:,icell)
   	endif   					

	Sth(:,icell) = eta_c(:,icell) / (tiny_dp + chi_c(:,icell)-sca_c(:,icell))
!    	Sth(:,icell) = Sth(:,icell) / ( tiny_dp + kappa_tot(:,icell) * (1.-beta(:,icell)))

   	if (any_nan_infinity_vector(Sth(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in Sth after interpolation.."
   	 write(*,*) "icell=", icell, "Sth=",Sth(:,icell)
   	 write(*,*) "kappa_abs=", ( kappa_tot(:,icell) * (1.-beta(:,icell)))
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

  labs = .true.
  id = 1
  etape_start = 2
  if (laccurate_integ) then
  	write(*,*) " Using step 3"
  	etape_end = 3
  else
  	etape_end = 2
  endif

if (write_convergence_file ) then
 open(unit=20, file="Jnu_convergence.s", status="unknown")
 write(20,*) "maxIter:", maxIter, " Nlam:", Nlambda, " Ncells:", n_cells
endif
 
     do etape=etape_start, etape_end

      if (etape==1) then 
	    call error("not step 1!")
      else if (etape==2) then 
  		lfixed_rays = .true.
  		n_rayons = n_rayons_start
  		iray_start = 1
  		lprevious_converged = .false.
		write(*,*) "Step 2: Using ", n_rayons, " rays for Jnu."
		lcell_converged(:) = .false.
  	  else if (etape==3) then 
  		write(*,*) "Step 3: max_rayons =  ", n_rayons_max, " n_rayons_start = ", n_rayons_start2
  		lfixed_rays = .false.
  		n_rayons = n_rayons_start2
  		iray_start = 1
  		lprevious_converged = .false.
		lcell_converged(:) = .false.
	  else 
	  	call error("step unknown")
  	  end if
  	  
if (write_convergence_file ) write(20,*) "step:", etape, " nrays:", n_rayons

  	  
  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0
		 				

        do while (.not.lconverged)

        	n_iter = n_iter + 1
        	!if (n_iter > 100 * maxIter)  exit
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
            !$omp private(argmt,norme, icell, l_iterate, delta_J) &
            !$omp shared(lambda_cont, lambda_star, Snew, Sold, Sth, Istar, tau_tot) &
            !$omp shared(lkeplerian,n_iter) &
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,lcell_converged) &
            !$omp shared(n_cells,ds, Jold, Jnu_cont, beta, kappa_tot, Ic,icompute_atomRT) &
            !$omp shared(lfixed_Rays,lnotfixed_Rays,labs,etape,pos_em_cellule)
            !$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !$ id = omp_get_thread_num() + 1
				if (lnotfixed_rays) then
   					l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
				else
   					l_iterate = (icompute_atomRT(icell)>0)
				endif				
   				if (l_iterate) then
   				   Jnu_cont(:,icell) = 0.
           		   Snew(:,icell) = 0. 
           		   lambda_star(:,id) = 0.0_dp                  		   		   
           		   
					do iray=iray_start, iray_start-1+n_rayons
					
    					if (etape==1) then
							stop
      	 				else! if (etape==2) then
                   	    ! Position aleatoire dans la cellule
         					rand  = sprng(stream(id))
            				rand2 = sprng(stream(id))
            				rand3 = sprng(stream(id))

            				call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
            				rand = sprng(stream(id))
           					W0 = 2.0_dp * rand - 1.0_dp
            				W02 =  1.0_dp - W0*W0
            				SRW02 = sqrt(W02)
            				rand = sprng(stream(id))
            				ARGMT = PI * (2.0_dp * rand - 1.0_dp)
            				U0 = SRW02 * cos(ARGMT)
            				V0 = SRW02 * sin(ARGMT)
            				
							call INTEG_RAY_JNU(id, icell, x0, y0, z0, u0, v0, w0, iray, labs, kappa_tot, Sold, Istar, Ic, tau_tot)
												
							!LI
! 		                	Jnu_cont(:,icell) = Jnu_cont(:,icell) + Ic(:,id)/n_rayons
		                	!ALI
		                	!lambda_star(:,id) = lambda_star(:,id) + exp(-tau_tot(:,id)) * (1d0 - exp(-ds(iray,id)*kappa_tot(:,icell))) / n_rayons
           					
		                	!hogerheijde
		                	Jnu_cont(:,icell) = Jnu_cont(:,icell) + ( Ic(:,id)*exp(-ds(iray,id)*kappa_tot(:,icell)) + &
 		                	(1.0_dp - exp(-ds(iray,id)*kappa_tot(:,icell)))*Sold(:,icell) ) / n_rayons

		 				end if !etape
  		                 					
      			    enddo !iray

      		   endif !icompute_AtomRT

			   !delta_J(:) = beta(:,icell) * (Jnu_cont(:,icell) - Jold(:,icell)) / (1.0 - beta(:,icell) * lambda_star(:,id))
			   !Snew(:,icell) = Sold(:,icell) + delta_J(:)
        	   Snew(:,icell) = Sth(:,icell) * (1.0_dp - beta(:,icell)) + Jnu_cont(:,icell) * beta(:,icell)

	   
     		end do !icell
        	!$omp end do
        	!$omp end parallel

            !should be para
        	diff = 0d0
        	dSource = 0.0_dp
  			cell_loop2 : do icell=1, n_cells
				if (lnotfixed_rays) then
   					l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
				else
   					l_iterate = (icompute_atomRT(icell)>0)
				endif
  				if (l_iterate) then
  					
  						dN = maxval(abs(Snew(:,icell) - Sold(:,icell))/Snew(:,icell))
  						dSource = max(dSource, dN)

						dJ = 0.0_dp
						do la=1, Nlambda_cont
						 if (Jnu_cont(la, icell) > 0) then 
						  dJ = max(dJ,abs(1.-Jold(la,icell)/Jnu_cont(la,icell)))
						  imax = locate(abs(1.-Jold(:,icell)/Jnu_cont(:,icell)),abs(1.-Jold(la,icell)/Jnu_cont(la,icell)))
						 endif 
						enddo
						!if (mod(icell,10)==0)
						write(*,'((1I5)" ::> dJ="(1ES14.5E3), " Jmax="(1ES14.5E3), " Jmin="(1ES14.5E3), " beta="(1ES14.5E3))') icell, real(dJ), maxval(Jnu_cont(:,icell)), minval(Jnu_cont(:,icell)), maxval(beta(:,icell))
if (write_convergence_file ) write(20,'((1I5)" ::> dJ="(1ES14.5E3), " Jmax="(1ES14.5E3), " Jmin="(1ES14.5E3), " beta="(1ES14.5E3))') icell, real(dJ), maxval(Jnu_cont(:,icell)), minval(Jnu_cont(:,icell)), maxval(beta(:,icell))
						if (dJ > diff) then
						  diff = dJ
						  icell_max = icell
						endif
     					lcell_converged(icell) = (real(dJ) < precision)		
     			end if
!      			Jold(:,icell) = Jnu_cont(:,icell)
!      			Sold(:,icell) = Snew(:,icell)
     		end do cell_loop2 !icell
     		Sold(:,:) = Snew(:,:)
     		Jold(:,:) = Jnu_cont(:,:)

     		write(*,'(" >>> dS = "(1ES14.5E3)," T(icell_max)="(1F14.5)" K", " ne(icell_max)="(1ES14.5E2))') dSource, T(icell_max), ne(icell_max)
if (write_convergence_file ) write(20,'(" >>> dS = "(1ES14.5E3)," T(icell_max)="(1F14.5)" K", " ne(icell_max)="(1ES14.5E2))') dSource, T(icell_max), ne(icell_max)
			write(*,'("  -- beta ="(1ES14.5E3))') beta(imax,icell_max)
if (write_convergence_file ) write(20,'("  -- beta="(1ES14.5E3))') beta(imax,icell_max)
         	write(*,'(" >>> icell_max "(1I5)," lambda="(1F14.7) " nm"," dJ="(1ES14.5E3))') icell_max, lambda(imax), diff
if (write_convergence_file ) write(20,'(" >>> icell_max "(I1)," lambda="(1F14.7)" nm"," dJ="(1ES14.5E3))') icell_max, lambda(imax), diff

         	
         	write(*,'(" @icell_max : Jmax="(1ES14.5E3)," Jmin="(1ES14.5E3) )') maxval(Jnu_cont(:,icell_max)), minval(Jnu_cont(:,icell_max))
         	write(*,'(" global     : Jmax="(1ES14.5E3)," Jmin="(1ES14.5E3) )') maxval(Jnu_cont(:,:)), minval(Jnu_cont(:,:))
if (write_convergence_file ) write(20,'(" @icell_max : Jmax="(1ES14.5E3)," Jmin="(1ES14.5E3) )') maxval(Jnu_cont(:,icell_max)), minval(Jnu_cont(:,icell_max))
if (write_convergence_file ) write(20,'(" global     : Jmax="(1ES14.5E3)," Jmin="(1ES14.5E3) )') maxval(Jnu_cont(:,:)), minval(Jnu_cont(:,:))
	        write(*,"('# unconverged cells :'(1I5), '('(1F12.3)' %)')") size(pack(lcell_converged,mask=lcell_converged.eqv..false.)), &
	          100.*real(size(pack(lcell_converged,mask=lcell_converged.eqv..false.)))/real(n_cells)
if (write_convergence_file ) write(20,"('# unconverged cells : '(1I5), '('(1F12.3)' %)')") size(pack(lcell_converged,mask=lcell_converged.eqv..false.)), &
	          100.*real(size(pack(lcell_converged,mask=lcell_converged.eqv..false.)))/real(n_cells)


! 			lconverged = (real(diff) < precision)        	
        	if (real(diff) < precision) then
           		if (lprevious_converged) then
            	  lconverged = .true.
           		else
            	  lprevious_converged = .true.
          	    endif
        	else
           		lprevious_converged = .false.
           		if (.not.lfixed_rays) then
           			n_rayons_old = n_rayons
              		n_rayons = n_rayons * 2
              		write(*,*) ' -- Increasing number of rays :', n_rayons
              		!!write(*,"('   '(1I10)' -> '(1I10))") n_rayons_old, n_rayons
              		if (n_iter >= maxIter) then
             		 	write(*,*) "Warning : not enough rays to converge !!"
                 		lconverged = .true.
              		end if

          	   end if
        	end if

	    end do !while
        write(*,*) etape, "Threshold =", precision
	  end do !over etapes
if (write_convergence_file ) close(20)
	
  
  if (.not.lstop_after_jnu) then
    do icell=1, n_cells
     if (icompute_atomRT(icell)>0) then
  	  call bezier2_interp(Nlambda_cont, lambda_cont, sca_c(:,icell) * Jnu_cont(:,icell), Nlambda, lambda, eta_es(:,icell))
  	 else
  	  eta_es(:,icell) = 0.0_dp
  	 endif
    enddo
    write(*,*) "Maxval eta_es", maxval(eta_es), minval(eta_es)
   
  endif
 
  
  
  open(unit=20, file="Jnu_no_interp.s", status="unknown")
  write(20,*) n_cells, Nlambda_cont
  do icell=1, n_cells
  	do la=1, Nlambda_cont !if you change the format, check read_Jnu_ascii()
  		if (lambda_cont(la) < 1d5) then
    		write(20,'(1F12.5,5E20.7E3)') lambda_cont(la), Jnu_cont(la,icell), Sth(la,icell), beta(la,icell), kappa_tot(la,icell), Sold(la,icell)
    	else
    		write(20,'(1F15.5,5E20.7E3)') lambda_cont(la), Jnu_cont(la,icell), Sth(la,icell), beta(la,icell), kappa_tot(la,icell), Sold(la,icell)
    	endif
   	enddo
  enddo
  close(20)


  deallocate(ds, lcell_converged)
  deallocate(Jold, Sth, kappa_tot, beta, Istar, Ic, Sold, lambda_star, Snew, delta_J, tau_tot)

 ! ------------------------------------------------------------------------------------ !
 return
 end subroutine Iterate_Jnu
 

! 	subroutine INTEG_RAY_LINE_I_CNTRB(id,icell_in,x,y,z,u,v,w,iray,labs)
! 	! ------------------------------------------------------------------------------- !
! 	! Computes the contribution functions along to the Intensities using converged
! 	! populations
! 	! ------------------------------------------------------------------------------- !
! 		integer, intent(in) :: id, icell_in, iray
! 		real(kind=dp), intent(in) :: u,v,w
! 		real(kind=dp), intent(in) :: x,y,z
! 		logical, intent(in) :: labs !used in NLTE but why?
! 		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
! 		real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, chiI, LimbD, chiIc
! 		real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chil, etal
! 		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
! 		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars
! 		integer :: la
! 
! 
! 		x1=x;y1=y;z1=z
! 		x0=x;y0=y;z0=z
! 		next_cell = icell_in
! 		nbr_cell = 0
! 
! 		tau(:) = 0.0_dp
! 		tau_c(:) = 0.0_dp
! 		chiI(:) = 0.0_dp
! 		chiIc(:) = 0.0_dp
! 
! 
! 		NLTEspec%I(:,iray,id) = 0d0
! 		NLTEspec%Ic(:,iray,id) = 0d0
! 
! 		!S_contrib(:,icell,id) = 0.0_dp
! 
!   ! -------------------------------------------------------------- !
!   !*** propagation dans la grille ***!
!   ! -------------------------------------------------------------- !
!   ! Will the ray intersect a star
! 		call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)
! 
!   ! Boucle infinie sur les cellules
! 		infinie : do ! Boucle infinie
!     ! Indice de la cellule
! 			icell = next_cell
! 			x0=x1 ; y0=y1 ; z0=z1
! 
! 
! 			if (icell <= n_cells) then
! 				lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
! 				if (atmos%icompute_atomRT(icell) < 0) return !-1 if dark
! 			else
! 				lcellule_non_vide=.false.
! 			endif
!     
!     ! Test sortie ! "The ray has reach the end of the grid"
! 
! 			if (test_exit_grid(icell, x0, y0, z0)) return
! 
! 			if (lintersect_stars) then
! 				if (icell == icell_star) then
! 
! 					call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)
! 
! 					NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*exp(-tau)
!        				NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*exp(-tau_c)
!        				return
!       			end if
!    			 endif
! 
! 			nbr_cell = nbr_cell + 1
! 
! 
! 			previous_cell = 0 ! unused, just for Voronoi
! 			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
! 
! 
! 
! 			if (lcellule_non_vide) then
! 				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
!      ! opacities in m^-1
! 				l_contrib = l_contrib * AU_to_m !l_contrib in m
! 				if ((nbr_cell == 1).and.labs) ds(iray,id) = l * AU_to_m
! 
! 
! 				call initAtomOpac(id)
! 				if (atmos%NpassiveAtoms>0) &
! 					call metal_bb(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l)
! 					
! 				chil(:) = NLTEspec%AtomOpac%chi_p(:,id)
! 				etal(:) = NLTEspec%AtomOpac%eta_p(:,id)
! 
! 				chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell) + tiny_dp 
! 				chiI(:)  = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)
! 
! 				if (atmos%NactiveAtoms>0) then 
! 					
! 					call initAtomOpac_nlte(id)
!       
! 					call NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l)
! 					
! 					chil(:) = chil(:) + NLTEspec%ATomOpac%chi(:,id)
! 					etal(:) = etal(:) + NLTEspec%AtomOpac%eta(:,id)
! 					
! 					chiIc(:) = chiIc(:) + NLTEspec%AtomOpac%Kc_nlte(:,icell)
! 
! 					chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id) +  NLTEspec%AtomOpac%Kc_nlte(:,icell)
! 				endif !active atoms
! 
! 				dtau(:)   = l_contrib * chiI(:)
! 				dtau_c(:) = l_contrib * chiIc(:)
!             
! 				Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
! 				Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
!       
!       			!Always include electron scattering if NLTE loop.
!       			!But ATM, uses %Jc even for lines.
! 				if (atmos%Nactiveatoms>0) then 
! 					Snu = Snu + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell) ) / chiI(:)
! 					Snu_c = Snu_c + ( NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) + &
! 					NLTEspec%AtomOpac%jc_nlte(:,icell) )/ chiIc(:)
! 
! 				else if ((atmos%electron_scattering).and.(atmos%NactiveAtoms==0)) then
! 					Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
! 					Snu_c = Snu_c + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) / chiIc(:)
! 				endif
! 
! 
! ! 				NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu
! ! 				NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + exp(-tau_c) * (1.0_dp - exp(-dtau_c)) * Snu_c
! 
! 				do la=1, NLTEspec%Nwaves
! 				
! 					!!Seems that chil * Ic << etal / chiI even for absorption lines	
! 					!S_contrib(la,icell,id) =  chil(la) * ( NLTEspec%Ic(la,iray,id) - etal(la) / chiI(la) ) * exp(-tau(la))
! 										
! 					NLTEspec%I(la,iray,id) = NLTEspec%I(la,iray,id) + exp(-tau(la)) * (1.0_dp - exp(-dtau(la))) * Snu(la)
! 					NLTEspec%Ic(la,iray,id) = NLTEspec%Ic(la,iray,id) + exp(-tau_c(la)) * (1.0_dp - exp(-dtau_c(la))) * Snu_c(la) 
! 					
! 					S_contrib(la,icell,id) =  (etal(la) / chiI(la)) * exp(-tau(la))
! 					!S_contrib(la,icell,id) =  etal(la) * exp(-tau(la))			
! 			
! 					!->Flow chart
! 					!S_contrib(la,icell,id) = ( etal(la)/chiI(la) ) * (1.0_dp - exp(-dtau(la))) * exp(-tau(la))
! 
! 				enddo
!      
! 
!      
! 				tau = tau + dtau
! 				tau_c = tau_c + dtau_c
! 
! 			end if  ! lcellule_non_vide
! 		end do infinie
! 	return
! 	end subroutine INTEG_RAY_LINE_I_CNTRB
  

end module atom_transfer
