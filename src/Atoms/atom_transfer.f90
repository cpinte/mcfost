! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atoms
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

module atom_transfer

	use atom_type, only			: AtomType
	use opacity, only			: dealloc_atom_quantities, alloc_atom_quantities, compute_atom_quantities, compute_background_continua, &
									interp_background_opacity, opacity_atom_loc, interp_contopac, interp_continuum_local, & !compute_nlte_bound_free, &
									nlte_bound_free, background_continua, Uji_down, chi_up, chi_down, eta_atoms, cross_coupling_terms, background_continua_lambda, opacity_atom_zeeman_loc, &
									prec_pops, frac_limit_pops
	use background_opacity, only: Thomson
	use Planck, only			: bpnu
	use spectrum_type, only     : dk, dk_max, dk_min, sca_c, chi, eta, chi_c, eta_c, eta_c_nlte, chi_c_nlte, eta0_bb, chi0_bb, lambda, Nlambda, lambda_cont, Nlambda_cont, Itot, Icont, Istar_tot, Istar_cont, &
									Stokes_Q, Stokes_U, Stokes_V, Flux, Fluxc, F_QUV, rho_p, etaQUV_p, chiQUV_p, cntrb_ray, init_spectrum, init_spectrum_image, dealloc_spectrum, Jnu_cont, eta_es, alloc_flux_image, allocate_stokes_quantities, &
									dealloc_jnu, reallocate_rays_arrays, write_contribution_functions_ray, write_flux, write_1D_arr_ascii, cntrb, flow_chart, write_lambda_cell_array, cf_file, flow_file
									
	use atmos_type, only		: nHtot, icompute_atomRT, lmagnetized, ds, Nactiveatoms, Atoms, calc_ne, Natom, ne, T, vr, vphi, v_z, vtheta, wght_per_H, &
									readatmos_ascii, dealloc_atomic_atmos, ActiveAtoms, nHmin, hydrogen, helium, lmali_scheme, lhogerheijde_scheme, &
									compute_angular_integration_weights, wmu, xmu, xmux, xmuy, v_char, angular_quadrature, Taccretion, laccretion_shock, ntotal_atom, &
									Ncentre, frac_pos, compute_anglequad_centres, largest_integer_smaller_than_x, healpix_sphere, healpix_npix, healpix_weight, healpix_mu_and_phi
									
	use readatom, only			: readAtomicModels, cswitch_enabled, maxval_cswitch_atoms, adjust_cswitch_atoms
	use lte, only				: set_LTE_populations, nH_minus, ltepops, ltepops_h
	use constant, only			: MICRON_TO_NM, hc_k, sigma_e
	use solvene, only			: solve_electron_density_old, solve_electron_density
	use getlambda, only			: hv
	use voigtfunctions, only	: Voigt, VoigtHumlicek, dirac_line
	use profiles, only			: profile, local_profile_v, local_profile_thomson, local_profile_interp, local_profile_dk
	use io_atomic_pops, only	: write_pops_atom, write_electron, read_electron, write_hydrogen, write_Hminus, write_convergence_map_atom, write_convergence_map_electron, prepare_check_pointing
	use io_opacity, only 		: write_Jnu, write_taur, write_contrib_lambda_ascii, read_jnu_ascii, Jnu_File_ascii, read_Jnu, &
									write_collision_matrix_atom, write_collision_matrix_atom_ascii, &
									write_radiative_rates_atom, write_rate_matrix_atom, write_cont_opac_ascii
	use math
	!$ use omp_lib
	use molecular_emission, only: v_proj
	use input, only				: lkeplerian, linfall, RT_line_method, nb_proc, RT_line_method, limb_darkening, mu_limb_darkening
	use parametres, only		: Rmax, Rmin, map_size, zoom, n_cells, lcontrib_function, lorigin_atom, lelectron_scattering, n_rad, nz, n_az, distance, ang_disque, &
									l_sym_ima, etoile, npix_x, npix_y, npix_x_save, npix_y_save, lpluto_file, lmodel_ascii, density_file, lsolve_for_ne, ltab_wavelength_image, &
									lvacuum_to_air, n_etoiles, lread_jnu_atom, lstop_after_jnu, llimb_darkening, dpops_max_error, laccurate_integ, NRAYS_ATOM_TRANSFER, &
									DPOPS_SUB_MAX_ERROR, n_iterate_ne,lforce_lte, loutput_rates, ing_norder, ing_nperiod, ing_ndelay, lng_acceleration, mem_alloc_tot, &
									ndelay_iterate_ne, llimit_mem, lfix_backgrnd_opac, lsafe_stop, safe_stop_time, checkpoint_period, lcheckpoint, istep_start, lno_iterate_ne_mc, &
									healpix_lorder, healpix_lmin, healpix_lmax

	use grid, only				: test_exit_grid, cross_cell, pos_em_cellule, move_to_grid
	use dust_transfer, only		: compute_stars_map
	use dust_ray_tracing, only	: RT_n_incl, RT_n_az, init_directions_ray_tracing,tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az,tab_RT_incl, stars_map, kappa
	use stars, only				: intersect_spots, intersect_stars
	!use wavelengths, only		:
	use mcfost_env, only		: dp, time_begin, time_end, time_tick, time_max
	use constantes, only		: tiny_dp, huge_dp, au_to_m, pc_to_au, deg_to_rad, tiny_real, pi, deux_pi, rad_to_deg, masseH, sigma
	use utils, only				: rotation_3d, cross_product
	use naleat, only 			: seed, stream, gtype
	use cylindrical_grid, only	: volume, r_grid, z_grid, phi_grid, cell_map_i, cell_map_j, cell_map_k, area
	use spherical_grid, only	: solid_angle_cell_sph
	use messages, only 			: error, warning
	
	use statequil_atoms, only   : invpop_file, profiles_file, unit_invfile, unit_profiles, calc_bb_rates, calc_bf_rates, calc_rate_matrix, update_populations, fill_collision_matrix, &
									init_bb_rates_atom, initgamma, initgamma_atom , init_rates_atom, store_radiative_rates,store_radiative_rates_mali,calc_rates, store_rate_matrices, &
									psi, calc_rates_mali, n_new, radiation_free_pops_atom, omega_sor_atom
	use collision, only			: CollisionRate !future deprecation
	use impact, only			: collision_hydrogen

	implicit none

 !Look-up table of exp(-tau) if interpolation is used. ATM too slow
!  	integer, parameter :: N_etau_table = 100
!  	real, parameter :: table_tau_min = 0.0, table_tau_max = 650.0
!  	real(kind=dp), dimension(N_etau_table) :: table_etau, table_xtau
 !OpenMp and debug
 	real(kind=dp), allocatable :: convergence_map(:,:,:,:) !say, N_cells, Nlevel, Natom, Nstep
 	integer, dimension(:), allocatable :: npix_x_id, npix_y_id, cell_id
 	real(kind=dp), dimension(:,:,:), allocatable :: xyz_pos, uvw_pos, chi_loc
 	real(kind=dp), dimension(:,:), allocatable :: threeKminusJ, Stest_tot, psi_mean
 	!!real(kind=dp), allocatable :: I0_rf(:), Iray_rf(:,:)
 !RT procedure
	procedure(integ_ray_line_i), pointer :: integ_ray_line => NULL()
 !Temporary variable for Zeeman calculations
	real(kind=dp), dimension(:,:), allocatable :: QUV
 !NLTE
 	real(kind=dp) :: dne
	logical :: ljacobi_sor = .false., lfixed_J = .false. !Option to compute J assuming LTE and keep this value for NLTE transfer. eta_es = ne * Jcont in that case
! 	logical :: lNg_acceleration = .true.
! 	integer :: iNg_Norder=2, iNg_ndelay=5, iNg_Nperiod=5
	real(kind=dp), allocatable :: ng_cur(:)
	logical, allocatable :: lcell_converged(:)
	real(kind=dp), allocatable :: gpop_old(:,:,:), Tex_old(:,:,:), ngpop(:,:,:), ne_old(:)
	real(kind=dp), allocatable :: Gammaij_all(:,:,:,:), Rij_all(:,:,:), Rji_all(:,:,:) !need a flag to output it
	integer :: NmaxLevel, NmaxTr
	!Check-pointing and stopping and timing
	logical :: lexit_after_nonlte_loop = .false.
	real :: time_iteration, time_nlte


 	contains


  	subroutine integ_ray_line_i(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	!
	! ------------------------------------------------------------------------------- !

		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, edt, et
		real(kind=dp), dimension(Nlambda) :: Snu, tau, dtau!, LD
		real(kind=dp), dimension(Nlambda_cont) :: Snu_c, dtau_c, tau_c!, LDc
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la, icell_prev
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp
		tau_c(:) = 0.0_dp

		Itot(:,iray,id) = 0d0
		Icont(:,iray,id) = 0d0
				

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
    		icell_prev = icell !duplicate with previous_cell, but this avoid problem with Voronoi grid here
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
					!It is completely possible to merge Icont and Itot in the say function
					!as the operation are similar
					!!call calc_stellar_surface_brightness(Nlambda_cont,lambda_cont,i_star, icell_prev,x0, y0, z0, u,v,w,LDc)
       				!!Icont(:,iray,id) =  Icont(:,iray,id) + LDc(:) * Istar_cont(:,i_star)*exp(-tau_c(:))
					!!call calc_stellar_surface_brightness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w,LD)
					!!Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * LD(:)
       				Icont(:,iray,id) =  Icont(:,iray,id) + exp(-tau_c) * Istar_cont(:,i_star) * local_stellar_brigthness(Nlambda_cont,lambda_cont,i_star, icell_prev,x0, y0, z0, u,v,w)
					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * local_stellar_brigthness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w)
       				return
      			end if
   			 endif

			nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
			previous_cell = 0 ! unused, just for Voronoi
			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
			
    !count opacity only if the cell is filled, else go to next cell
			if (lcellule_non_vide) then
				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet same as iterate?
     ! opacities in m^-1

				l_contrib = l_contrib * AU_to_m !l_contrib in m
			
				!total bound-bound + bound-free + background opacities lte + nlte
				if (llimit_mem) then
					call interp_continuum_local(icell, chi(:,id), eta(:,id))
				else
					chi(:,id) = chi0_bb(:,icell)
					eta(:,id) = eta0_bb(:,icell)				
				endif

				!includes a loop over all bound-bound, passive and active
				call opacity_atom_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,( (nbr_cell==1).and.labs ) ) 
				
				dtau(:) = l_contrib * chi(:,id)
				
				if ((nbr_cell == 1).and.labs) then
					if (lmali_scheme) then
						if (minval(chi(:,id)) == 0) write(*,*) lambda(locate(chi(:,id), minval(chi(:,id)))), T(icell), ne(icell), nhtot(icell)
						if (any_nan_infinity_vector(chi(:,id)) > 0) then
							write(*,*) chi(:,id)
							call error( "inf or nan in chi")
						endif
						psi(:,iray,id) = ( 1.0_dp - exp( -l*AU_to_m*chi(:,id) ) ) / chi(:,id)
						if (allocated(chi_loc)) chi_loc(:,iray,id)  = chi(:,id)
					else
						ds(iray,id) = l*au_to_m
					endif
				endif
				            
				if (lelectron_scattering) then
					Snu = ( eta(:,id) + eta_es(:,icell) ) / ( chi(:,id) + tiny_dp )
					Snu_c = eta_c(:,icell) + thomson(icell) * Jnu_cont(:,icell)
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

				!Itot(:,iray,id) = Itot(:,iray,id) + exp(-tau(:)) * (1.0_dp - exp(-dtau(:))) * Snu(:)
!-> faster written that way ?
!! 				Itot(:,iray,id) = Itot(:,iray,id) + ( exp(-tau(:)) - exp(-tau(:) - dtau(:))) * Snu(:)
!
! 				Icont(:,iray,id) = Icont(:,iray,id) + exp(-tau_c(:)) * (1.0_dp - exp(-dtau_c(:))) * Snu_c(:)

! 				tau = tau + dtau
! 				tau_c = tau_c + dtau_c

!-> avoid multiple loop over frequencies and is faster for large Nlambda
				do la=1,Nlambda
					Itot(la,iray,id) = Itot(la,iray,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * Snu(la)
					tau(la) = tau(la) + dtau(la) !for next cell
				enddo

				do la=1, Nlambda_cont
					Icont(la,iray,id) = Icont(la,iray,id) + ( exp(-tau_c(la)) - exp(-(tau_c(la) + dtau_c(la))) ) * Snu_c(la)
					tau_c(la) = tau_c(la) + dtau_c(la)
				enddo
				
! 				if (minval(tau) > 1000.) then
! 					write(*,*) "taumin > 100 stopping propagation"
! 					return
! 				endif


			end if  ! lcellule_non_vide
		end do infinie

	return
	end subroutine integ_ray_line_i
	
  	subroutine integ_ray_line_flow(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	! For contribution functions and origin emission
	! ------------------------------------------------------------------------------- !

		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, edt, et
		real(kind=dp), dimension(Nlambda) :: Snu, tau, dtau, chi_line, tau_line, Sline!, LD
		real(kind=dp), dimension(Nlambda_cont) :: Snu_c, dtau_c, tau_c!, LDc
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la, icell_prev
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp
		tau_c(:) = 0.0_dp

		Itot(:,iray,id) = 0d0
		Icont(:,iray,id) = 0d0
		
		tau_line = 0.0_dp
		chi_line = 0.0_dp
		Sline = 0.0_dp
				

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
    		icell_prev = icell
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
! 					call calc_stellar_surface_brightness(Nlambda_cont,lambda_cont,i_star, icell_prev, x0, y0, z0, u,v,w,LDc)
!        			Icont(:,iray,id) =  Icont(:,iray,id) + LDc(:) * Istar_cont(:,i_star)*exp(-tau_c(:))
! 					call calc_stellar_surface_brightness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w,LD)
! 					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * LD(:)
       				Icont(:,iray,id) =  Icont(:,iray,id) + exp(-tau_c) * Istar_cont(:,i_star) * local_stellar_brigthness(Nlambda_cont,lambda_cont,i_star, icell_prev,x0, y0, z0, u,v,w)
					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * local_stellar_brigthness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w)
       				return
      			end if
   			 endif

			nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
			previous_cell = 0 ! unused, just for Voronoi
			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
			
    !count opacity only if the cell is filled, else go to next cell
			if (lcellule_non_vide) then
				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet same as iterate?
     ! opacities in m^-1

				l_contrib = l_contrib * AU_to_m !l_contrib in m
			
				!total bound-bound + bound-free + background opacities lte + nlte
! 				!rename it chi0_bckgr etc
				if (llimit_mem) then
					call interp_continuum_local(icell, chi(:,id), eta(:,id))
				else
					chi(:,id) = chi0_bb(:,icell)
					eta(:,id) = eta0_bb(:,icell)				
				endif

				!includes a loop over all bound-bound, passive and active
				call opacity_atom_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,( (nbr_cell==1).and.labs ) ) 
				
				Sline = eta(:,id) - eta0_bb(:,icell)
				chi_line = chi(:,id) - chi0_bb(:,icell)
				
				dtau(:) = l_contrib * chi(:,id)
				            
				if (lelectron_scattering) then
					Snu = ( eta(:,id) + eta_es(:,icell) ) / ( chi(:,id) + tiny_dp )
					Snu_c = eta_c(:,icell) + thomson(icell) * Jnu_cont(:,icell)
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
				
				if (lorigin_atom .and. lcontrib_function) then
					do la=1,Nlambda
						if (x0>=0 .and. y0>=0 .and. z0>=0) then
							flow_chart(la,icell) = flow_chart(la,icell) + ( exp(-tau_line(la)) - exp(-(tau_line(la)+l_contrib * chi_line(la))) ) * Sline(la) / (tiny_dp + chi_line(la))!Snu(la)
							if (tau_line(la) > 0) cntrb(la,icell) = cntrb(la,icell) + E2(tau(la)) * eta(la,id) !tau(la) * Snu(la) in logtau
						endif
						Itot(la,iray,id) = Itot(la,iray,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * Snu(la)
						tau(la) = tau(la) + dtau(la) !for next cell
						tau_line(la) = tau_line(la) + l_contrib * chi_line(la)
					enddo
				else if (lorigin_atom) then
					do la=1,Nlambda
						if (x0>=0 .and. y0>=0 .and. z0>=0) then
							flow_chart(la,icell) = flow_chart(la,icell) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * Snu(la)
						endif
						Itot(la,iray,id) = Itot(la,iray,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * Snu(la)
						tau(la) = tau(la) + dtau(la) !for next cell
					enddo
				else if (lcontrib_function) then
					do la=1,Nlambda
						if (x0>=0 .and. y0>=0 .and. z0>=0) then
							if (tau(la) > 0) cntrb(la,icell) = cntrb(la,icell) + E2(tau(la)) * eta(la,id)
						endif
						Itot(la,iray,id) = Itot(la,iray,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * Snu(la)
						tau(la) = tau(la) + dtau(la) !for next cell
					enddo
				endif
				

				!cont rad only
				do la=1, Nlambda_cont
					Icont(la,iray,id) = Icont(la,iray,id) + ( exp(-tau_c(la)) - exp(-(tau_c(la) + dtau_c(la))) ) * Snu_c(la)
					tau_c(la) = tau_c(la) + dtau_c(la)
				enddo


			end if  ! lcellule_non_vide
		end do infinie

	return
	end subroutine integ_ray_line_flow

  	subroutine integ_ray_line_z(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	!
	! ------------------------------------------------------------------------------- !

		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, edt, et
		real(kind=dp), dimension(Nlambda) :: Snu, tau, dtau, S_Q, S_U, S_V!, LD
		real(kind=dp), dimension(Nlambda_cont) :: Snu_c, dtau_c, tau_c!, LDc
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la, icell_prev
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp
		tau_c(:) = 0.0_dp

		Itot(:,iray,id) = 0d0
		Icont(:,iray,id) = 0d0
		
 		Stokes_V(:,id) = 0d0
  		Stokes_Q(:,id) = 0d0
  		Stokes_U(:,id) = 0d0
				

  ! -------------------------------------------------------------- !
  !*** propagation dans la grille ***!
  ! -------------------------------------------------------------- !
  ! Will the ray intersect a star
		call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! Boucle infinie sur les cellules
		infinie : do ! Boucle infinie
    ! Indice de la cellule
    		icell_prev = icell
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
! 					call calc_stellar_surface_brightness(Nlambda_cont,lambda_cont,i_star, icell_prev, x0, y0, z0, u,v,w,LDc)
!        				Icont(:,iray,id) =  Icont(:,iray,id) + LDc(:) * Istar_cont(:,i_star)*exp(-tau_c(:))
! 					call calc_stellar_surface_brightness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w,LD)
! 					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * LD(:)
       				Icont(:,iray,id) =  Icont(:,iray,id) + exp(-tau_c) * Istar_cont(:,i_star) * local_stellar_brigthness(Nlambda_cont,lambda_cont,i_star, icell_prev,x0, y0, z0, u,v,w)
					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * local_stellar_brigthness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w)
       				return
      			end if
   			 endif

			nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
			previous_cell = 0 ! unused, just for Voronoi
			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
			
    !count opacity only if the cell is filled, else go to next cell
			if (lcellule_non_vide) then
				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet same as iterate?
     ! opacities in m^-1

				l_contrib = l_contrib * AU_to_m !l_contrib in m
			
				!total bound-bound + bound-free + background opacities lte + nlte
				if (llimit_mem) then
					call interp_continuum_local(icell, chi(:,id), eta(:,id))
				else
					chi(:,id) = chi0_bb(:,icell)
					eta(:,id) = eta0_bb(:,icell)				
				endif

				!includes a loop over all bound-bound, passive and active
				call opacity_atom_zeeman_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,( (nbr_cell==1).and.labs ) ) 
				
				dtau(:) = l_contrib * chi(:,id)

				            
				if (lelectron_scattering) then
					Snu = ( eta(:,id) + eta_es(:,icell) ) / ( chi(:,id) + tiny_dp )
					Snu_c = eta_c(:,icell) + thomson(icell) * Jnu_cont(:,icell)
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
				
    			!explicit product of Seff = S - (K/chiI - 1) * I
    			!Is there a particular initialization to do for polarisation ?
    			!because it will be zero at the first place
     			Snu(:) = Snu(:) -chiQUV_p(:,1,id) / chi(:,id) * Stokes_Q(:,id) - &
          					chiQUV_p(:,2,id) / chi(:,id) *  Stokes_U(:,id) - &
          					chiQUV_p(:,3,id) / chi(:,id) *  Stokes_V(:,id)


     			S_Q(:) = etaQUV_p(:,1,id) /chi(:,id) - &
							chiQUV_p(:,1,id) / chi(:,id) * Itot(:,iray,id) - &
							rho_p(:,3,id)/chi(:,id) * Stokes_U(:,id) + &
							rho_p(:,2,id)/chi(:,id) * Stokes_V(:,id)

     			S_U(:) = etaQUV_p(:,2,id) /chi(:,id) - &
							chiQUV_p(:,2,id) / chi(:,id) * Itot(:,iray,id) + &
							rho_p(:,3,id)/chi(:,id) * Stokes_Q(:,id) - &
							rho_p(:,1,id)/chi(:,id) * Stokes_V(:,id)

     			S_V(:) = etaQUV_p(:,3,id) /chi(:,id) - &
     						chiQUV_p(:,3,id) / chi(:,id) * Itot(:,iray,id) - &
							rho_p(:,2,id)/chi(:,id) * Stokes_Q(:,id) + &
							rho_p(:,1,id)/chi(:,id) * Stokes_U(:,id)


				do la=1,Nlambda
					Itot(la,iray,id) = Itot(la,iray,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * Snu(la)
					Stokes_Q(la,id) = Stokes_Q(la,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * S_Q(la)
					Stokes_U(la,id) = Stokes_U(la,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * S_U(la)
					Stokes_V(la,id) = Stokes_V(la,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau(la))) ) * S_V(la)

					tau(la) = tau(la) + dtau(la) !for next cell
				enddo

				do la=1, Nlambda_cont
					Icont(la,iray,id) = Icont(la,iray,id) + ( exp(-tau_c(la)) - exp(-(tau_c(la) + dtau_c(la))) ) * Snu_c(la)
					tau_c(la) = tau_c(la) + dtau_c(la)
				enddo


			end if  ! lcellule_non_vide
		end do infinie

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
   real(kind=dp), dimension(3) :: sdx, sdy
   real(kind=dp):: npix2, diff, normF, R0
   real(kind=dp), parameter :: precision = 1.e-2
   integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter
   logical :: lintersect, labs
   integer :: kr, nr, nb, nat
   type (AtomType), pointer :: atom

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
     if (lmagnetized) QUV(:,:) = 0d0 !move outside
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

             I0 = I0 + Itot(:,iray,id)
             I0c = I0c + Icont(:,iray,id)
                          
             if (lmagnetized) then
             	QUV(:,3) = QUV(:,3) + Stokes_V(:,id) / npix2
             	QUV(:,2) = QUV(:,2) + Stokes_Q(:,id) / npix2
             	QUV(:,1) = QUV(:,1) + Stokes_U(:,id) / npix2
             end if

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

  ! Prise en compte de la surface du pixel (en sr)

  ! Flux out of a pixel in W/m2/Hz/pix
  normF = ( pixelsize / (distance*pc_to_AU) )**2

  ! adding to the total flux map.
  !storing by proc is necessary here
  Flux(:,ibin,iaz,id) = Flux(:,ibin,iaz,id) + I0(:) * normF
  Fluxc(:,ibin,iaz,id) = Fluxc(:,ibin,iaz,id) + I0c(:) * normF
  
  
  if (lmagnetized) then
  	F_QUV(:,ibin,iaz,1) = F_QUV(:,ibin,iaz,1) + normF * QUV(:,1) / npix2
  	F_QUV(:,ibin,iaz,2) = F_QUV(:,ibin,iaz,2) + normF * QUV(:,2) / npix2
  	F_QUV(:,ibin,iaz,3) = F_QUV(:,ibin,iaz,3) + normF * QUV(:,3) / npix2
  end if
  
  
  !Flux map for lines
  !adding a map for the continuum point too ?
  if (RT_line_method==1) then
  !summation over pixels
  	do nat=1,Natom
  		atom => atoms(nat)%ptr_atom
  		do kr=1,atom%Nline
  		
  			if (atom%lines(kr)%write_flux_map) then
  				nr = atom%lines(kr)%Nred + dk_max
  				nb = atom%lines(kr)%Nblue + dk_min
  			
  				atom%lines(kr)%map(:,1,1,ibin,iaz) = &
  					atom%lines(kr)%map(:,1,1,ibin,iaz) + I0(nb:nr)*normF
  				
  			endif
  		
  		enddo
  		atom => NULL()
  	enddo
  else!2D maps
  	do nat=1,Natom
  		atom => atoms(nat)%ptr_atom
  		do kr=1,atom%Nline
  		
  			if (atom%lines(kr)%write_flux_map) then
  				nr = atom%lines(kr)%Nred + dk_max
  				nb = atom%lines(kr)%Nblue + dk_min
  			
  				atom%lines(kr)%map(:,ipix,jpix,ibin,iaz) = I0(nb:nr)*normF
  				
  			endif
  		
  		enddo
  		atom => NULL()
  	enddo
  
  endif


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

  integer, parameter :: n_rad_RT = 151, n_phi_RT = 100 !(100, 36)
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp):: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda
  logical :: lresolved
  real(kind=dp) :: aphi(n_phi_RT), ap(n_rad_RT)

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
  !-> at the moment it is like rt_line_method==2 but with explicit integration over the map
!     n_iter_min = 1
!     n_iter_max = 1
!     
!     ! dx and dy are only required for stellar map here
!     taille_pix = (map_size/zoom)  ! en AU
!     dx(:) = x_plan_image * taille_pix
!     dy(:) = y_plan_image * taille_pix
! 
!     i = 1
!     j = 1
!     lresolved = .false.
! 
!     rmin_RT = max(w*0.9_dp,0.05_dp) * Rmin
!     rmax_RT = Rmax * 2.0_dp
! 
!     tab_r(1) = rmin_RT
!     fact_r = exp( (1.0_dp/(real(n_rad_RT,kind=dp) -1))*log(rmax_RT/rmin_RT) )
!     do ri_RT = 2, n_rad_RT
!       tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
!     enddo
! 
!     fact_A = sqrt(pi * (fact_r - 1.0_dp/fact_r)  / n_phi_RT )
! 
! 
!     ! Boucle sur les rayons d'echantillonnage
!     !$omp parallel &
!     !$omp default(none) &
!     !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
!     !$omp shared(tab_r,fact_A,x,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
!     !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz,fact_r)
!     id = 1 ! pour code sequentiel
! 
!     if (l_sym_ima) then
!       cst_phi = pi  / real(n_phi_RT,kind=dp)
!     else
!       cst_phi = deux_pi  / real(n_phi_RT,kind=dp)
!     endif
! 
!      !$omp do schedule(dynamic,1)
!      do ri_RT=1, n_rad_RT
!         !$ id = omp_get_thread_num() + 1
!         r = tab_r(ri_RT)
! 		
!         taille_pix =  fact_A * r ! racine carree de l'aire du pixel
! 
!         do phi_RT=1,n_phi_RT ! de 0 + eps à 2pi - eps (eps = pi/n_phi_RT)
!            phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)
! 			
!            pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image
!             ! C'est le centre en fait car dx = dy = 0.
!            call flux_pixel_line(id,ibin,iaz,n_iter_min,n_iter_max, i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
!         end do !j
!      end do !i
!      !$omp end do
!      !$omp end parallel

     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = (map_size/zoom) / real(n_rad_RT,kind=dp) ! en AU
     lresolved = .false.

     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

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
     do i = 1,n_rad_RT
        !$ id = omp_get_thread_num() + 1
        do j = 1,n_rad_RT
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           call flux_pixel_line(id,ibin,iaz,n_iter_min,n_iter_max,1,1,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
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
  
  !recombine for all proc
  flux(:,ibin,iaz,1) = sum(flux(:,ibin,iaz,:),dim=2)
  fluxc(:,ibin,iaz,1) = sum(fluxc(:,ibin,iaz,:),dim=2)

 return
 end subroutine emission_line_map

 subroutine atom_line_transfer()
 ! --------------------------------------------------------------------------- !
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
 ! --------------------------------------------------------------------------- !
  integer :: atomunit = 1, nact
  integer :: icell, m
  integer :: ibin, iaz, kr
  integer :: n_rayons_start
  !!integer, parameter :: n_rayons_start = 100
  integer, parameter :: n_rayons_start2 = 1000, maxIter = 11
  integer :: n_rayons_max = n_rayons_start2 * (2**(maxIter-1))
  integer, parameter :: Nrayone = 1
  character(len=20)  :: ne_initial
  logical :: lread_jnu_ascii = .false., lelectron_read
  type (AtomType), pointer :: atom
  integer :: alloc_status
  
  !init at 0
  mem_alloc_tot = 0

  lhogerheijde_scheme = .false. !not ready yet 
  lmali_scheme = (.not.lhogerheijde_scheme)

 ! -------------------------------INITIALIZE AL-RT ------------------------------------ !
  !only one available yet, I need one unpolarised, faster and more accurate.
  Voigt => VoigtHumlicek
  profile => local_profile_interp 
  integ_ray_line => integ_ray_line_i

  if (npix_x_save > 1) then
   RT_line_method = 2 ! creation d'une carte avec pixels carres
   npix_x = npix_x_save ; npix_y = npix_y_save
  else
   RT_line_method = 1 !pixels circulaires
  end if


 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
        !! ----------------------- Read Model ---------------------- !!

  if (.not.lpluto_file) then
   if (lmodel_ascii) then
    call readAtmos_ascii(density_file)
   end if
  end if
  
  if (lmagnetized) then
  	!or need to allocate line%adamp after or interpolate the dispersion profile
  	write(*,*) "Cannot yet accomodate local_profile_interp with polarised calc"
  	profile => local_profile_v
  endif 

        !! --------------------------------------------------------- !!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------- READATOM and INITIZALIZE POPS ------------------------------ !
  call readAtomicModels(atomunit)!if old_pops not compatible with fixed background
  								 !because of the electron density evaluated in non-LTE
  								 !if background fixed, they will be computed with the current estimate of
  								 !electron density and non-LTE pops !

  !not in the model but fits file exists from previous run ?
  call read_electron(lelectron_read)
  if (lelectron_read) then
  	ne_initial = "NE_MODEL"
  	calc_ne = lsolve_for_ne !only if we force solving it. Otherwise it is read
  	!from old values and eventually iterated during non-LTE loop. This way,
  	!the iteration starts exactly where it stopped.
  	!Forcing evaluation of ne (lsolve_for_ne == True) might smooth the solution.
  else
    ne_initial = "H_IONISATION"
	write(*,*) "Solving for electron density from H+M ionisation"
  endif	
 
  !no electron density in the model nor previous fits file calc_ne == True.
  !if a previous file is found (lelectron_read == True) and lsolve_for_ne, then
  !calc_ne is set to .true. Electron density is re-evaluated using the populations
  !in the fits files (H_Ionisation otherwise).
  if (calc_ne) then !eventually computed with the previous non-LTE pops ! 
	call Solve_Electron_Density(ne_initial, .true., dne) 
! 	call Solve_Electron_Density_old(ne_initial)
	call write_Electron
  else
  	if (lsolve_for_ne) then
  		ne_initial = "NE_MODEL"
		write(*,*) "Forcing calculation of electron density"	
		write(*,*) " solution from model!"
		call Solve_Electron_Density(ne_initial, .true., dne) 
		call write_Electron 
  	endif
  endif

  	call set_LTE_populations !write pops at the end because we probably have NLTE pops also

	if (NactiveAtoms > 0) then


		n_rayons_start = Nrays_atom_transfer
	
		!I do not handle step 3 at the moment as it can be two long
		!so laccurate_integ is used for step 2!
! 		if (.not.laccurate_integ) then
! 			if (lmali_scheme) then
! 				n_rayons_max = 1
! 			else
! 				n_rayons_max = n_rayons_start
! 			endif
! 		endif
		if (lmali_scheme) then
			n_rayons_max = 1
		else
			n_rayons_max = n_rayons_start
		endif
		
		call init_Spectrum(n_rayons_max,lam0=500d0,vacuum_to_air=lvacuum_to_air)


		if (n_etoiles > 0) call init_stellar_disk
		write(*,*) " Computing background opacities..."!+ init nlte cont opac + interpolation for lines
		call alloc_atom_quantities !call compute_atom_quantities(icell) for each cell

		call compute_background_continua

		!use the same rays as nlteloop
		if (lelectron_scattering) then
			!otherwise, J = 0 at init for nlte loop
			if (lfixed_J) call iterate_Jnu
		endif
		
		if (allocated(ds)) deallocate(ds)
! 		if (lmali_scheme) then !not use full here
! 			allocate(ds(1, nb_proc), stat=alloc_status)
! 		else
		if (lhogerheijde_scheme) then
			allocate(ds(n_rayons_max, nb_proc), stat=alloc_status)
			if (alloc_status > 0) call error ("Allocation error ds")
		endif
		!allocated(dk ....)
		
		NmaxLevel = 0
		NmaxTr = 0
		do nact=1,Nactiveatoms
			atom => ActiveAtoms(nact)%ptr_atom
			atom%NLTEpops = .true.
			!now we have computed background, lte pops and electron density we can set it to true
			!to use nlte populations in electron density
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
			
			!calc_delta_Tex_atom
			
			if (atom%initial_solution=="ZERO_RADIATION") then
				write(*,*) "-> Initial solution at SEE with I = 0 for atom ", atom%ID
				!factorize in a new subroutine for all cells ?
				do icell=1,n_cells
					if (icompute_atomRT(icell) > 0) then
						call radiation_free_pops_atom(1, icell, atom, .false.)
					endif
				enddo
    			do m=1,atom%Nlevel
    				write(*,"('Level #'(3I1))") m
    				write(*,'("  -- min(n)="(1ES20.7E3)" m^-3; max(n)="(1ES20.7E3)" m^-3")') , minval(atom%n(m,:),mask=(icompute_atomRT>0)), maxval(atom%n(m,:))
    				write(*,'("  -- min(nstar)="(1ES20.7E3)" m^-3; max(nstar)="(1ES20.7E3)" m^-3")')  minval(atom%nstar(m,:),mask=(icompute_atomRT>0)), maxval(atom%nstar(m,:))
				enddo
			endif
			mem_alloc_tot = mem_alloc_tot + sizeof(atom%C) + sizeof(atom%Gamma)
			atom => NULL()
		enddo
		!remove some arrays in depth to save memory, like vbroad per atom not needed, it is fast to compute locally
   
		allocate(lcell_converged(n_cells),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error lcell_converged")
		lcell_converged = .false.
		mem_alloc_tot = mem_alloc_tot + sizeof(lcell_converged)

		allocate(Tex_old(NactiveAtoms,NmaxTr,n_cells),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error Tex_old")
		Tex_old = 0.0_dp
		mem_alloc_tot = mem_alloc_tot + sizeof(Tex_old)
		
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
			write(*,*) " Size of radiative rates and rate matrix:", (2*sizeof(Rij_all)+sizeof(Gammaij_all))/1024./1024., " MB"
		endif

		write(*,*) " ** Number max of levels among all atoms:", Nmaxlevel
		write(*,*) " ** Number max of transitions among all atoms:", NmaxTr


		if (NactiveAtoms > 1) then
			write(*,*) "   -> Solving for kinetic equations for ", Nactiveatoms, " atoms"
		else
			write(*,*) "   -> Solving for kinetic equations for ", Nactiveatoms, " atom"
		endif
		write(*,*) " Max error : ", dpops_max_error, dpops_sub_max_error


		!!allocate(psi_mean(nlambda, n_cells)); psi_mean = 0.0
		!!allocate(chi_loc(nlambda, n_rayons_max, nb_proc)); chi_loc = 0.0

		write(*,'("->Total memory allocated before NLTEloop:"(1ES17.8E3)" MB")') mem_alloc_tot/1024./1024.
		
		if (lcheckpoint) call prepare_check_pointing()

		call NLTEloop_mali(n_rayons_max, n_rayons_start, n_rayons_start2, maxIter)
		
		!reevaluate here only if fixed during the loop but ne was iterated
		!otherwise constant because 1) ne constant or 2), ne evaluated just after
		if ((lfix_backgrnd_opac).and.(n_iterate_ne > 0)) then
			call compute_background_continua
		endif
		
		!evaluate electron density once SEE is solved
		if (n_iterate_ne == 0) then
			write(*,*) "Evaluate Electron density from converged populations..."
			do icell=1,n_cells
				if ((icompute_atomRT(icell) > 0).and.(lelectron_scattering)) eta_es(:,icell) = eta_es(:,icell) / thomson(icell)
			enddo
			ne_initial = "NE_MODEL"
			call solve_electron_density(ne_initial, .true., dne)
			write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)

			do icell=1,n_cells
				if (icompute_atomRT(icell) > 0) then
					if (lelectron_scattering) eta_es(:,icell) = eta_es(:,icell) * thomson(icell)
					nHmin(icell) = nH_minus(icell)
				endif
			enddo
			do nact=1, Natom!see calc_ionisation_frac, for why we updated also passive atoms
				if (Atoms(nact)%ptr_atom%ID=="H") then
					write(*,*) " -> updating LTE pops of hydrogen with non-LTE ne"
					call LTEpops_H
				else
					write(*,*) " -> updating LTE pops of ",Atoms(nact)%ptr_atom%ID, " with non-LTE ne"
					call LTEpops(Atoms(nact)%ptr_atom,.false.)
				endif	
			enddo
			
			call compute_background_continua
							
			call write_electron
			call write_Hminus
			write(*,*) "...done."
		endif


		deallocate(lcell_converged, Tex_old)
		if (allocated(gpop_old)) deallocate(gpop_old)
		if (allocated(dk)) deallocate(dk)
		if (allocated(ds)) deallocate(ds)
		
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
			!!if (associated(hydrogen,atom)) call write_collision_matrix_atom(hydrogen)
			call write_collision_matrix_atom(atom)
			do m=1, atom%Nline
				deallocate(atom%lines(m)%phi_loc)
			enddo
			atom => NULL()
		enddo
		if (loutput_rates) deallocate(Rij_all,Rji_all,Gammaij_all)
		
		if (lelectron_scattering) call write_Jnu

		
		!-> I do not do that for nlte calc otherwise, Eta_es is thomson * Jnu_cont
		!if (lelectron_scattering) call iterate_Jnu
		!-> write_Jnu and iterate_Jnu ??


		!->> Case of ltab or not and add chi_c_nlte and eta_c_nlte
		!then if ltab_wavelength_image ...
		if (ltab_wavelength_image) then
			call error("tab wavelength not implemented")
		else
			!just remove the Nrays dependencies
			call reallocate_rays_arrays(nrayOne)
		endif

	else !no nlte or using old populations so no chi_c_nlte no eta_c_nlte and atom is passive with old_populations
  
  		if (ltab_wavelength_image) then !atomic lines transfer with either lte or nlte populations
  									!using user defined grid
  						
  			call error("tab wavelength not implemented")			
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
					!call read_jnu_ascii
					!lread_jnu_ascii = .true.
					call read_jnu
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
    	!call interp_contopac !-> done in compute_background_continua because now there is 
    							! a seperation between eta_es and Jcont
		
	endif !NLTE

	!for all so after nlte case (debug)
	do m=1,Natom
		call write_pops_atom(Atoms(m)%ptr_atom)
	end do
	

	if (lexit_after_nonlte_loop) then
		return
		!free all data or because we leave after not needed
	endif

	if (lmagnetized) then
		integ_ray_line => integ_ray_line_z
		lorigin_atom = .false.
		lcontrib_function = .false.
 		allocate(QUV(Nlambda,3))
		call allocate_stokes_quantities
	endif

	call alloc_flux_image()
	write(*,*) "Computing emission flux map..."
	write(*,'("Total memory allocated for image calculation:"(1ES17.8E3)" GB")') mem_alloc_tot / 1024./1024./1024.

	if (lorigin_atom .or. lcontrib_function) then
		integ_ray_line => integ_ray_line_flow 
	endif

	call init_directions_ray_tracing()
	do ibin=1,RT_n_incl
		do iaz=1,RT_n_az
			call emission_line_map(ibin,iaz)
		end do
	end do
	write(*,*) " ..done"
	
	if (lmagnetized) then
		deallocate(QUV)
	endif


	if (allocated(flow_chart)) then
		call write_lambda_cell_array(flow_chart, flow_file)
		deallocate(flow_chart)
	endif	
	if (allocated(cntrb)) then
		call write_lambda_cell_array(cntrb, cf_file)
		deallocate(cntrb)
	endif
	
	call write_flux
	deallocate(flux,fluxc)

	!one ray
! 	if (lcontrib_function) then
! 		call compute_contribution_functions
! 	end if
	
    
    !ascii files can be large !
    if (loutput_rates) then
		call write_cont_opac_ascii(hydrogen)
 		!chi0_bb and eta0_bb not modified
 		do kr=1,hydrogen%Nline
 			if (hydrogen%lines(kr)%write_flux_map) then
				call write_contrib_lambda_ascii(hydrogen%lines(kr)%lambda0,0.0_dp,0.0_dp,1.0_dp,.true.)
			endif
		enddo
		call write_taur(500._dp,0._dp,0._dp,1.0_dp)
	endif

	if (lelectron_scattering.and.(NactiveAtoms==0)) then !otherwise if active, written even if not included in emissivity
		!Jnu is written to ascii file if read
		if (lread_jnu_atom) then
			if (lread_jnu_ascii) call write_Jnu
		endif
		call dealloc_Jnu()
	endif

	if (allocated(psi_mean)) deallocate(psi_mean, chi_loc)
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

	subroutine NLTEloop_mali(n_rayons_max,n_rayons_1,n_rayons_2,maxIter)
	! ----------------------------------------------------------------------- !
	! removed all dependency on nrayons since everything is ray by ray built
	! ----------------------------------------------------------------------- !
#include "sprng_f.h"
		integer, intent(in) :: n_rayons_1, n_rayons_2, maxIter, n_rayons_max
		integer :: etape, etape_start, etape_end, iray, n_rayons, iray_p
		integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, la, kr
		integer, dimension(nb_proc) :: max_n_iter_loc
		logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged
		real :: rand, rand2, rand3, precision, fac_etape!, precision_sub
		real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, weight
		real(kind=dp) :: diff, norme, dN, dN1, dJ, lambda_max
		real(kind=dp) :: dT, dN2, dN3, dN4, diff_old
		real(kind=dp), allocatable :: dTM(:), dM(:), Tion_ref(:), Tex_ref(:)
		real(kind=dp), allocatable :: Jnew(:,:), Jnew_cont(:,:), Jnu_loc(:,:)
		!healpix
		integer :: l_order
		real, parameter :: critical_ratio = 4 !12,4, dOmega*Npix/4pi < critical_ratio
		real(kind=dp) :: domega, healpix_phi
		!
! 		real(kind=dp), allocatable :: err_pop(:,:)
		logical :: labs, iterate_ne = .false.
		logical :: l_iterate
		logical :: accelerated, ng_rest, evaluate_background, lmean_intensity = .false.!,lapply_sor_correction
		integer :: iorder, i0_rest, n_iter_accel, iacc!, iter_sor
		integer :: nact, imax, icell_max, icell_max_2
		integer :: icell, ilevel, imu, iphi, id_ref !for anisotropy
		character(len=20) :: ne_start_sol = "NE_MODEL"
		type (AtomType), pointer :: atom
		integer(kind=8) :: mem_alloc_local = 0
		
		write(*,*) " USING MALI METHOD FOR NLTE LOOP"
		
		!only for step 1 at the moment
		allocate(threeKminusJ(nlambda, n_cells))
		!or (nlambda, nb_proc) if only used locally as a threshold
		allocate(Jnu_loc(nlambda,nb_proc))
		id_ref = locate(lambda, 300.0_dp)
		!write anisotropy at ref wavelength only otherwise too big file (in ascii, n_lambda * n_cells * 8 bits)
		!if used, we only need it locally or at one wavelenegth
		
		!time for individual steps + check the time from the checkpointing time if any
		!and if it is set lconverged = .true. and .lprevious converged == true
 		call system_clock(time_begin,count_rate=time_tick,count_max=time_max)
		
		!we are missing some information about dJ in the case lmean_intensity is .false.
		!how to correct that without much overhead and without allocating Jnew, Jold (or eta_es new/old)
		if (lelectron_scattering) lmean_intensity = .true.
		
		!If we iterate electron scattering with the SEE. if == 0, ONCE SEE is set.
		iterate_ne = (n_iterate_ne>0)

		!!open(unit=unit_invfile, file=trim(invpop_file), status="unknown")
		!!write(unit_invfile,*) n_cells
		

! 		if (ljacobi_sor) then
! 			iter_sor = 0
! 			lapply_sor_correction = .false.
! 			allocate(omega_sor_atom(NactiveAtoms)); omega_sor_atom(:) = 1.0_dp
! 			allocate(err_pop(Nactiveatoms,3)); err_pop(:,:) = 0.0_dp
! 		endif

! 		if (laccurate_integ) then
! 			if (iterate_ne) then
! 				allocate(convergence_map(n_cells,NmaxLevel,1+NactiveAtoms,2))
! 			else
! 				allocate(convergence_map(n_cells,NmaxLevel,NactiveAtoms,2))			
! 			endif
! 		else
! 			if (iterate_ne) then
! 				allocate(convergence_map(n_cells,NmaxLevel,1+NactiveAtoms,1))
! 			else
! 				allocate(convergence_map(n_cells,NmaxLevel,NactiveAtoms,1))			
! 			endif
! 		endif
! 		convergence_map = -1.0_dp
  
		allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
		allocate(dTM(Nactiveatoms)); dM=0d0 !keep tracks of Tex for all cells for each atom
		if (lmean_intensity) then 
			allocate(Jnew(Nlambda, n_cells)); Jnew = 0.0
! 			allocate(Jold(Nlambda, n_cells)); Jold = 0.0
			allocate(Jnew_cont(Nlambda_cont, n_cells)); Jnew_cont = 0.0
			write(*,*) "size Jtot:", 2*sizeof(Jnew)/1024./1024.," MB" !+eta_es
			write(*,*) "size Jcont:", 2*sizeof(Jnu_cont)/1024./1024.," MB"
			mem_alloc_local = mem_alloc_local + 2*(sizeof(Jnew)+sizeof(Jnew_cont))
		endif
		
		allocate(Tex_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tex for all cells for each line of each atom
		allocate(Tion_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tion for all cells for each cont of each atom
		diff_old = 1.0_dp
		
		mem_alloc_local = mem_alloc_local + sizeof(dM)+sizeof(dTm)+sizeof(Tex_ref)+sizeof(Tion_ref)
		
		if (lNg_acceleration) then 
			n_iter_accel = 0
			i0_rest = 0
			ng_rest = .false.
			iacc = 0
			allocate(ngpop(n_cells * NmaxLevel, iNg_Norder+2, NactiveAtoms), stat=alloc_status)
			if (alloc_status > 0) then
				call error("Cannot allocate Ng table !")
			endif
			write(*,*) " Size ngpop:", sizeof(ngpop)/1024./1024., " MB"
			mem_alloc_local = mem_alloc_local + sizeof(ngpop)
		endif
		
		deallocate(stream)
		allocate(stream(nb_proc),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error stream")

		labs = .true. !to have ds at cell icell = eval_operator
		id = 1
		!In case we restart a calculation in the step 2
		etape_start = istep_start!1
		if (laccurate_integ) then
			etape_end = 2!3, step 3 not implemented yet
		else
			etape_end = 1
			if (istep_start==2) etape_end = 2
		endif
		
		if (angular_quadrature == "HEALpix_adapt") then
			etape_start = 4
			etape_end = 4
		endif
		
		write(*,*) "---------------------------------------"
		write(*,*) " step start ",etape_start, istep_start
		write(*,*) " step end", etape_end
		write(*,*) "---------------------------------------"

		! ds is not used for this scheme at  the moment. Psi is used instead
		!since there is no local subit needing recomputing psi
		!! ds(:,:)
		allocate(n_new(NactiveAtoms,Nmaxlevel,n_cells),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error n_new")
		write(*,*) " size atom%n, n_new:", 2*sizeof(n_new) / 1024./1024.," MB"
		n_new(:,:,:) = 0.0_dp
		allocate(psi(Nlambda, n_rayons_max, nb_proc), stat=alloc_status); psi = 0.0_dp
		write(*,*) " size psi:", sizeof(psi) / 1024./1024.," MB"
		if (alloc_Status > 0) call error("Allocation error psi in nlte loop")
		allocate(eta_atoms(Nlambda,NactiveAtoms,nb_proc),stat=alloc_status)
		if (alloc_Status > 0) call error("Allocation error eta_atoms")
		write(*,*) " size eta_atoms:", sizeof(eta_atoms) / 1024./1024.," MB"
		allocate(Uji_down(Nlambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
		if (alloc_Status > 0) call error("Allocation error Uji_down")
		write(*,*) " size cross-coupling:", 3 * sizeof(Uji_down) / 1024./1024.," MB"
		allocate(chi_down(Nlambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
		if (alloc_Status > 0) call error("Allocation error chi_down")
		allocate(chi_up(Nlambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
		if (alloc_Status > 0) call error("Allocation error chi_up")
		
		mem_alloc_local = mem_alloc_local + sizeof(n_new) + sizeof(psi) + sizeof(eta_atoms) + &
			sizeof(uji_down)+sizeof(chi_down)+sizeof(chi_up)
			
		mem_alloc_tot = mem_alloc_tot + mem_alloc_local
		write(*,'("Total memory allocated in NLTEloop:"(1ES17.8E3)" MB")') mem_alloc_local / 1024./1024.
		write(*,'("Total memory allocated up to now:"(1ES17.8E3)" GB")') mem_alloc_tot / 1024./1024./1024.
				

		step_loop : do etape=etape_start, etape_end
		
			if (etape==1) then
				time_iteration = 0
				
				!Only one ray until a proper scheme for choosing positions inside the cell exist
      			call compute_angular_integration_weights()
      			call compute_anglequad_centres()

  				lfixed_rays = .true.
  				n_rayons = 1 !always
    			iray_start = 1

  				write(*,*) " Using step 1 with ", size(xmu), " rays"
  				if (Ncentre + 1 > 0) then
  					write(*,'("   --> using "(1I3)" ray centres.")') Ncentre + 1
  		  			write(*,'("   --> a total of "(1I5)" elements.")')( Ncentre + 1 ) * size(xmu) 
!   		  		if (Ncentre > 0) iray_start = 2 !exclude the centre of the cell ? 
  				endif

  				lprevious_converged = .false.
				lcell_converged(:) = .false.
				precision = dpops_max_error
  				 		
			else if (etape==2) then 		
				time_iteration = 0
				
				!or no ? because we resample the radiation field
! 				if (etape_start == 1) then
! 					Ndelay_iterate_ne = 0
! 				endif
				
				!-> no iteration in MC but after the solution ??
				if (iterate_ne) then
					if (lno_iterate_ne_mc) then
						iterate_ne = .false.
						n_iterate_ne = 0	
					endif
				endif
			
				write(*,*) " Using step 2 with ", n_rayons_1, " rays"
				lfixed_rays = .true.
				n_rayons = n_rayons_1
				iray_start = 1
				lprevious_converged = .false.
				lcell_converged(:) = .false.
				fac_etape = 0.1
				if (etape_start==1) then
					precision = fac_etape * 1e-1!1e-1!1e-3, 1e-2, fac_etape * 1.0 / sqrt(real(n_rayons))
				else
					precision = dpops_max_error				
				endif
				write(*,*) " threshold:", precision
				
				if (lNg_acceleration) then
					deallocate(ngpop)
					if (allocated(ng_cur)) deallocate(ng_cur)
					lNg_acceleration = .false.
				endif
				
			else if (etape==3) then
				time_iteration = 0
				write(*,*) " Using step 3 with ", n_rayons_max, " nrays max"
		  		lfixed_rays = .false.
  				n_rayons = n_rayons_2
  				lprevious_converged = .false.
				lcell_converged(:) = .false.
				precision = 0.1 !dpops_max_error
				
			else if (etape==4) then
				time_iteration = 0
				
      			call compute_angular_integration_weights()
  				lfixed_rays = .true.
  				n_rayons = 1
    			iray_start = 1

    			write(*,*) " Using step 1 with adapatative healpix"
    			l_order = healpix_lorder !init
    			write(*,*) " -> l_order max for cells:", max(min(largest_integer_smaller_than_x(log(etoile(1)%r**2 * critical_ratio*pi/3.0/minval(area))/log(4.0)),healpix_lmax),healpix_lmin)

  				lprevious_converged = .false.
				lcell_converged(:) = .false.
				precision = dpops_max_error
			else
				call ERROR("etape unkown")
			end if
  	  
			lnotfixed_rays = .not.lfixed_rays
			lconverged = .false.
			n_iter = 0

			do while (.not.lconverged)

				n_iter = n_iter + 1
				write(*,*) " -> Iteration #", n_iter, " Step #", etape
				!!write(unit_invfile,*) "************************************************"
				!!write(unit_invfile,*) "step ", etape, ' iter ', n_iter
				!!write(unit_invfile,*) "************************************************"

				!in step 3 we use the same seed (same sets of random numbers) to add rays to
				!the previously computed.
				if (lfixed_rays) then
					stream = 0.0
					do i=1,nb_proc
						stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
					end do
				end if

				max_n_iter_loc = 0
				!reset in case of step 3

        	
				!$omp parallel &
				!$omp default(none) &
				!$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dM, dN, dN1,iray_p,imu,iphi)&
				!$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, nact, atom, l_iterate, weight, l_order, domega, healpix_phi) &
				!$omp shared(icompute_atomRT, dpops_sub_max_error,lkeplerian,lforce_lte,n_iter, threeKminusJ, psi_mean, psi, chi_loc, Jnu_loc,healpix_lorder,healpix_lmin, healpix_lmax) &
				!$omp shared(stream,n_rayons,iray_start, area, r_grid, z_grid, phi_grid, lcell_converged,loutput_rates, Nlambda_cont, Nlambda, lambda_cont) &
				!$omp shared(n_cells, gpop_old,integ_ray_line, Itot, Icont, Jnu_cont, eta_es, xmu, xmux, xmuy,wmu,etoile,id_ref, frac_pos,Ncentre,xyz_pos,uvw_pos) &
				!$omp shared(Jnew, Jnew_cont, lelectron_scattering,chi0_bb, etA0_bb, T,eta_atoms, lmean_intensity) &
				!$omp shared(nHmin, chi_c, chi_c_nlte, eta_c, eta_c_nlte, ds, Rij_all, Rji_all, Nmaxtr, Gammaij_all, Nmaxlevel) &
				!$omp shared(lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape,pos_em_cellule,Nactiveatoms,lambda)
				!$omp do schedule(static,1)
				do icell=1, n_cells
					!$ id = omp_get_thread_num() + 1
   					l_iterate = (icompute_atomRT(icell)>0)
			
   					if (l_iterate) then
   						if (lmean_intensity) then
   							Jnew(:,icell) = 0.0
   							Jnew_cont(:,icell) = 0.0
   							!!psi_mean(:,icell) = 0.0
   						endif
   						threeKminusJ(:,icell) = 0.0
   						Jnu_loc(:,id) = 0.0

						call fill_Collision_matrix(id, icell) !computes = C(i,j) before computing rates
						call initGamma(id) !init Gamma to C and init radiative rates
						!psi(:,iray,id) = 0.0_dp !no use full here


						if ( etape == 2) then !ray-by-ray, n_rayons fixed
                   	    ! Position aleatoire dans la cellule
							do iray=iray_start, iray_start-1+n_rayons


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
								V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)
								

								call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, 1, labs)			

								if (lmean_intensity) then
									Jnew(:,icell) = Jnew(:,icell) + Itot(:,1,id) / n_rayons
									Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,1,id) / n_rayons
								endif

								!for one ray
								if (.not.lforce_lte) then
									call cross_coupling_terms(id, icell, 1)
									call calc_rates_mali(id, icell, 1, 1.0_dp / real(n_rayons,kind=dp))
								endif	
								
								if (loutput_Rates) then
									!need to take into account the fact that for MALI no quandities are store for all ray so Rij needs to be computed ray by ray
									call store_radiative_rates_mali(id, icell,(iray==1), 1.0_dp / real(n_rayons,kind=dp), Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell))
								endif	


							enddo !iray	
		
						
						else if (etape==3) then	!accumulate rays, n_rayons not fixed
						
							do iray=iray_start, iray_start-1+n_rayons !for the new added rays


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
								V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)
								

								call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, iray, labs)			

							enddo !iray		
														
							!for all rays accumulated
							do iray=1, n_rayons !for all
								if (lmean_intensity) then
									Jnew(:,icell) = Jnew(:,icell) + Itot(:,iray,id) / n_rayons
									Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,iray,id) / n_rayons
								endif
								
								!for one ray
								if (.not.lforce_lte) then
									call cross_coupling_terms(id, icell, iray)
									call calc_rates_mali(id, icell, iray, 1.0_dp / real(n_rayons,kind=dp))
								endif	
								
								!might not work cause profile of iray==1 alwaus used
								if (loutput_Rates) then
									call store_radiative_rates_mali(id, icell, (iray==1), 1.0_dp / real(n_rayons,kind=dp), Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell))
								endif	
							enddo			
								
						elseif (etape==1) then !ray-by-ray, n_rayons fixed

  		         			do iray = iray_start, Ncentre+1
  		         			
								if (iray==1) then
									x0 = r_grid(icell) * cos(phi_grid(icell))
									y0 = r_grid(icell) * sin(phi_grid(icell))
									z0 = z_grid(icell)
								else
									call  pos_em_cellule(icell,frac_pos(iray-1,1),frac_pos(iray-1,2),frac_pos(iray-1,3),x0,y0,z0)
								endif
							
  		         				do imu=1, size(xmu)
  		         					w0 = xmu(imu)
									u0 = xmux(imu); v0 = xmuy(imu)

									weight = wmu(imu) / real(Ncentre + 2 - iray_start)
									
									call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, 1, labs)			

									if (lmean_intensity) then
										Jnew(:,icell) = Jnew(:,icell) + Itot(:,1,id) * weight
										Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,1,id) * weight
										!!psi_mean(:,icell) = psi_mean(:,icell) + chi_loc(:,1,id) * psi(:,1,id) * weight
									endif
									
									threeKminusJ(:,icell) = threeKminusJ(:,icell) +  (3.0 * (u0*x0+v0*y0+w0*z0)**2/(x0**2+y0**2+z0**2) - 1.0) * Itot(:,1,id) * weight
									Jnu_loc(:,id) = Jnu_loc(:,id) + Itot(:,1,id) * weight

								!for one ray
									if (.not.lforce_lte) then
										call cross_coupling_terms(id, icell, 1)
										call calc_rates_mali(id, icell, 1, weight)
									endif	
								
									if (loutput_Rates) then
									!need to take into account the fact that for MALI no quandities are store for all ray so Rij needs to be computed ray by ray
										call store_radiative_rates_mali(id, icell, (iray==1 .and. imu==1), weight, Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell))
									endif	


      			   				enddo !imu	
      			   			enddo !pos / iray
      			   			
      			   			threeKminusJ(:,icell) = 0.5 * threeKminusJ(:,icell) / Jnu_loc(:,id)
      			   		
      			   		elseif (etape==4) then
      			   		
  		         			
							x0 = r_grid(icell) * cos(phi_grid(icell))
							y0 = r_grid(icell) * sin(phi_grid(icell))
							z0 = z_grid(icell)

								
							!domega = solid_angle_cell_sph(icell) * (x0*x0+y0*y0+z0*z0)/etoile(1)%r**2
							domega = area(icell) / etoile(1)%r**2
							!dOmega_cell / dOmega_healpix; dOmega_healpix = SQ(angular_resolution in rad) = pi/real(3*4**l_order)
							!dOmega = area_choc / (x0**2+y0**2+z0**2) !area_choc == sum(area(cell touching the star which acreates))
							l_order = max(min(largest_integer_smaller_than_x(log(critical_ratio*pi/3.0/domega)/log(4.0)),healpix_lmax),healpix_lmin)
							!if (n_iter==1) &
							!write(*,*) icell, " lorder=", l_order, healpix_npix(l_order), sqrt(x0*x0+y0*y0+z0*z0)/etoile(1)%r

							weight = healpix_weight(l_order)
							
  		         			do imu=1, healpix_npix(l_order)
  		         				call healpix_mu_and_phi(l_order,imu,w0,healpix_phi)
								u0 = sqrt(1.0 - w0*w0)*cos(healpix_phi)
								v0 = sqrt(1.0 - w0*w0)*sin(healpix_phi)

									
								call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, 1, labs)			

								if (lmean_intensity) then
									Jnew(:,icell) = Jnew(:,icell) + Itot(:,1,id) * weight
									Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,1,id) * weight
								endif
									
								threeKminusJ(:,icell) = threeKminusJ(:,icell) +  (3.0 * (u0*x0+v0*y0+w0*z0)**2/(x0**2+y0**2+z0**2) - 1.0) * Itot(:,1,id) * weight
								Jnu_loc(:,id) = Jnu_loc(:,id) + Itot(:,1,id) * weight

								!for one ray
								if (.not.lforce_lte) then
									call cross_coupling_terms(id, icell, 1)
									call calc_rates_mali(id, icell, 1, weight)
								endif	
								
								if (loutput_Rates) then
									!need to take into account the fact that for MALI no quandities are store for all ray so Rij needs to be computed ray by ray
									call store_radiative_rates_mali(id, icell, (iray==1 .and. imu==1), weight, Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell))
								endif	

      			   			enddo !imu	
      			   			
      			   			!if kept add for all steps
      			   			threeKminusJ(:,icell) = 0.5 * threeKminusJ(:,icell) / Jnu_loc(:,id)
!       			   			write(*,*) id, icell, " Anis = ", 100*abs(threeKminusJ(id_ref,icell)), " dOmega=", solid_angle_cell_sph(icell) * 3*(4**healpix_lorder)/pi      			   			
			
						end if !etape
						
						call calc_rate_matrix(id, icell, lforce_lte)
						call update_populations(id, icell, diff, .false., n_iter)


						n_iter_loc = 0
						if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
						
						if (loutput_Rates) then
							call store_rate_matrices(id,icell,Nmaxlevel,Gammaij_all(:,:,:,icell))
							!!call store_radiative_rates(id, icell, n_rayons, Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell), Jnew(:,icell), .false.)
						endif
	
					end if !icompute_atomRT
				end do !icell
				!$omp end do
				!$omp end parallel
				
				!Ng acceleration
   				accelerated = .false.
   				!do not accelerate if lfprevious_converged ??.and.((cswitch==1.0_dp).and..not.lprevious_converged))
   				if ( (lNg_acceleration .and. (n_iter > iNg_Ndelay)).and.(maxval_cswitch_atoms()==1.0_dp) ) then
   					!remove also the cswtich iteration from the calculation or accumulate immediately after cswitch == 1?
   					iorder = n_iter - iNg_Ndelay
   					if (ng_rest) then
						write(*,'(" -> Acceleration relaxes... "(1I2)" /"(1I2))') iorder-i0_rest, iNg_Nperiod
   						if (iorder-i0_rest == iNg_Nperiod) ng_rest = .false.
   					else
   						!Still, all atoms run at the same speed
   						i0_rest = iorder
            	        iacc = iacc + 1
            	        write(*,'(" -> Accumulate solutions... "(1I2)" /"(1I2))') iacc, iNg_Norder+2
   						do nact=1,NactiveAtoms
   							atom => ActiveAtoms(nact)%ptr_atom
   							allocate(ng_cur(n_cells * atom%Nlevel)); ng_cur(:) = 0.0_dp
   							!flatten2 works better with 2D arrays (tested on 1D and 2D models)
   							ng_cur = flatten2(atom%Nlevel, n_cells,n_new(nact,1:atom%Nlevel,:))
   								!has to be parallel in the future
   								
   							!for many atoms, increment niter only with the first one, as they go at the same speed.
   							accelerated = ng_accelerate(iacc, ng_cur, n_cells * atom%Nlevel, iNg_Norder, ngpop(1:atom%Nlevel*n_cells,:,nact), check_negative_pops=.false.)
   							
   							if (accelerated) then
   								!handle negative pops by simpling cancel Ng's iteration ?
   								!Extra overheads to do that. Might be better to raise an error since it can be
   								!a bug that needs proper correction.
!    							if (minval(ng_cur)<0) then
!    								accelerated = .false.
!    								write(*,*) "   --> avoiding acceleration"
!    							else
!                   				n_iter_accel = n_iter_accel + 1 !True number of accelerated iter
!             					write(*,'("     ++> Ngs iteration #"(1I4))') n_iter_accel

								!-> only reshape because we only print accelerated iteration for all atoms at once
                  					n_new(nact, 1:atom%Nlevel,:) = reform2(atom%Nlevel, n_cells, ng_cur)
!                   			endif
! 
! 								!after acceleration we rest for Nperiod iterations.
! 								!But if iNg_Nperiod == 0 we accumulate directly the next solutions.
!                   				ng_rest = (iNg_Nperiod > 0)!.true.
!                   				iacc = 0	
   							endif
   							
   							deallocate(ng_cur)
   							atom => NULL()
   						enddo
   						!for all atoms should be accelerated at the same time
   						if (accelerated) then
                  			n_iter_accel = n_iter_accel + 1 !True number of accelerated iter
            				write(*,'("     ++> Ng iteration #"(1I4))') n_iter_accel
                  			ng_rest = (iNg_Nperiod > 0)!.true.
                  			iacc = 0	
   						endif
   						
   					endif
   				endif
   				
   				if ((mod(n_iter, checkpoint_period)==0).and.(lcheckpoint)) then
   					do nact=1, NactiveAtoms
						call write_pops_atom(ActiveAtoms(nact)%ptr_atom,iter=n_iter,step=etape)
   					enddo
   				endif
   				
   				
			!evaluate electron density with the ionisation fraction
			!I have to evaluate background continua if I update all atoms lte pops ?
			!still update chi_c if electron scatt ?
				evaluate_background = .false.
				dne = 0.0_dp											
				!Updating for non Ng iterations ? add a condition on the dM < 1e-2 ? (better ndelay)
				if ((iterate_ne .and. (mod(n_iter,n_iterate_ne)==0)).and.(n_iter>ndelay_iterate_ne)) then!.and..not.accelerated)  then
					!do not evaluate if accelerated. Do not evaluate if cswitch ?? Or do not care?
					
					!chic(:,:) = chi(:,:) - ne(:)
					!allocate(ne_old(n_cells))
					!!ne_old(:) = ne(:)
					!write(*,'("ne old ="(1ES17.8E3)"("(1ES17.8E3)") m^-3")') maxval(ne), minval(ne,mask=ne>0)
					call solve_electron_density(ne_start_sol, .true., dne)
					!call Solve_Electron_Density_old(ne_start_sol)
					!dne = maxval(abs(1.0_dp - ne_old(:)/(ne(:)+tiny_dp)))
				!might be correct to do nstar_new = nstar_old * ne_new/ne_old ?
! -> If i update here, for the first iteration nold is not equal to the lte populations used for
! the first solution, so the first test on the convergence is wrong. But,
! after the first iteration it starts to  be consistent
					do nact=1, Natom!see calc_ionisation_frac, for why we updated also passive atoms
						!if evaluate background update passive atoms ?
						if (Atoms(nact)%ptr_atom%ID=="H") then
							write(*,*) " -> updating LTE pops of hydrogen"
							call LTEpops_H
							!check pops actually are updated ??
						else
							write(*,*) " -> updating LTE pops of ",Atoms(nact)%ptr_atom%ID 
							call LTEpops(Atoms(nact)%ptr_atom,.false.)
						endif	
					enddo
					write(*,*) ''
					!chi_c(:,:) = chi_c(:,:) + ne(:)
					!deallocate(ne_old)
					evaluate_background = (.not.lfix_backgrnd_opac) !.true.
					!convergence_map(:,1,NactiveAtoms+1,:) = dne
				end if

     		!Global convergence Tests
     			id = 1
            !should be para
				dM(:) = 0.0 !all pops
				diff = 0.0
				dJ = 0.0
				lambda_max = 0.0
				dT = 0.0
				dTM(:) = 0.0 !all temperature (ion+line)
				Tex_ref(:) = 0.0
				Tion_ref(:) = 0.0
   				icell_max = 1
   				icell_max_2 = 1
   				!for all cells
				dN2 = 0.0 !among all Tex
				dN4 = 0.0 !among all Tion
				cell_loop2 : do icell=1,n_cells
				
   					l_iterate = (icompute_atomRT(icell)>0)


   					if (l_iterate) then
   					
   						!Local only 
						dN = 0.0 !for all levels of all atoms of this cell
! 						dN2 = 0.0 !among all Tex
! 						dN4 = 0.0 !among all Tion
							
						do nact=1,NactiveAtoms
							atom => ActiveAtoms(nact)%ptr_atom
     					         					    
							do ilevel=1,atom%Nlevel
! 								if ( n_new(nact,ilevel,icell) >= frac_ne_limit * ne(icell) ) then
								if ( n_new(nact,ilevel,icell) >= frac_limit_pops * ntotal_atom(icell, atom) ) then
									dN1 = abs(1d0-atom%n(ilevel,icell)/n_new(nact,ilevel,icell))
									dN = max(dN1, dN)
									dM(nact) = max(dM(nact), dN1)
								endif
									!convergence_map(icell, ilevel, nact, etape) = dN1
							end do !over ilevel

							!I keep negative Temperatures for info.debug.
							!hence the /= 0.0. But, some work should be done in update_pops and compute_Tex
							do kr=1, atom%Nline
							
								if (atom%lines(kr)%Tex(icell) /= 0.0_dp) then 
									dN3 = abs(1.0 - Tex_old(nact, kr, icell) /  atom%lines(kr)%Tex(icell) )!( atom%lines(kr)%Tex(icell) + tiny_dp ) )
									dN2 = max(dN3,dN2)
									dTM(nact) = max(dTM(nact), dN3)
									if (dN3 >= dN2) then
										Tex_ref(nact) = atom%lines(kr)%Tex(icell)
										icell_max = icell
									endif
								endif
								
							enddo

							do kr=1, atom%Ncont
							
								if ( atom%continua(kr)%Tex(icell) /= 0.0_dp) then
									dN3 = abs(1.0 - Tex_old(nact, kr+atom%Nline, icell) /  atom%continua(kr)%Tex(icell) )
									dN4 = max(dN3,dN4)
									dTM(nact) = max(dTM(nact), dN3)
									if (dN3 >= dN4) then
										Tion_ref(nact) = atom%continua(kr)%Tex(icell)
										icell_max_2 = icell
									endif
								endif
																	
							enddo

							atom => NULL()
						end do !over atoms
						
						!compare for all atoms and all cells
						diff = max(diff, dN) ! pops

						
						!do not update if lfixed_J
						if ((lmean_intensity).and..not.(lfixed_J)) then
							Jnu_cont(:,icell) = Jnew_cont(:,icell)
							do la=1, Nlambda
								dN1 = abs( 1.0_dp - eta_es(la,icell)/(thomson(icell)*Jnew(la,icell)) )
								if (dN1 > dJ) then
									dJ = dN1
									lambda_max = lambda(la)
								endif
							enddo
						
							eta_es(:,icell) = Jnew(:,icell) * thomson(icell)
						endif
												
						lcell_converged(icell) = (real(diff) < precision) !(real(diff) < dpops_max_error)
						
						!Re init for next iteration if any
						do nact=1, NactiveAtoms
							atom => ActiveAtoms(nact)%ptr_atom
							!gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
							atom%n(:,icell) = n_new(nact,1:atom%Nlevel,icell)
							do kr=1,atom%Nline
								Tex_old(nact, kr, icell) = atom%lines(kr)%Tex(icell)
							enddo
							do kr=1, atom%Ncont
								Tex_old(nact, kr+Atom%Nline,icell) = atom%continua(kr)%Tex(icell)
							enddo
							!if (evaluate_background) then
							! !evalute LTE here ?? if so, need a local version of lte pops.
							! !or store the populations nstar in a new array?
							! recompute profiles or damping
							!
							!endif
							atom => NULL()
						end do
						
						!if I recompute all background continua ?  and damping ? evaluate damping locally ?
						if (evaluate_background) then
							nHmin(icell) = nH_minus(icell)
							!-> this evaluate profiles or damping but also continuum quantities not needed ???
							!-> do not check ni-njgij < 0 because a non converged state can lead to inversion of populations
							! modify it to update only nlte quantities ?
							call compute_atom_quantities(icell)
							call background_continua_lambda(icell, Nlambda_cont, lambda_cont, chi_c(:,icell), eta_c(:,icell))
							!or just re evaluate ne in chi_c like chi_c -ne_old + ne_new
						endif
						
						call NLTE_bound_free(icell)

						!because of opac nlte changed
						if (.not.llimit_mem) then
							call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
						endif
						!end init

						
					end if !if l_iterate
				end do cell_loop2 !icell
				!! test on the solar case
				!! use dN2 when step == 2 ?
				!! change how Tex are initialised (presenly with T, use pops instead)
   				!!diff = dN2 !dT(line)
   				
! 				if ((ljacobi_sor).and..not.(lapply_sor_correction)) then
! 
! 					nact = 1
! 					iter_sor = n_iter - 6 !3 previous dn/n are stored, and used for the next iterations
! 					if ((iter_sor > 0).and.(iter_sor < 4)) then
! 						err_pop(:,iter_sor) = dM(:)
! 						write(*,*) n_iter, iter_sor, err_pop(nact,iter_sor)
! 					endif
! 					
! 					!!if (iter_sor == 4) write(*,*) "alpha=", abs (err_pop(nact,1) - 2.0*err_pop(nact,2) + err_pop(nact,3) )
! 
! 					if ( (mod(n_iter,10) == 0).and.(err_pop(nact,3)/err_pop(nact,2) <= 1.0) ) then
! 						omega_sor_atom(nact) = 2.0_dp / ( 1.0 + sqrt(1.0 - err_pop(nact,3)/err_pop(nact,2) ) )
! 						lapply_sor_correction = .true.
! 						write(*,*) n_iter, "omega_sor=", omega_sor_atom(nact), err_pop(nact,:)
! 					endif
! 				endif

				if (maxval(max_n_iter_loc)>0) write(*,'(" -> "(1I10)" sub-iterations")') maxval(max_n_iter_loc)
				write(*,'(" -> icell_max1 #"(1I6)," icell_max2 #"(1I6))') icell_max, icell_max_2
				if (dJ /= 0.0_dp) write(*,'(" -> dJ="(1ES14.5E3)" @"(1F14.4)" nm")') dJ, lambda_max !at the end of the loop over n_cells
				write(*,*) " ------------------------------------------------ "
				do nact=1,NactiveAtoms
					write(*,'("             Atom "(1A2))') ActiveAtoms(nact)%ptr_atom%ID
					if (accelerated) then
						write(*,'("   >>> dpop="(1ES17.8E3)" (Accelerated)")') dM(nact)
					else
						write(*,'("   >>> dpop="(1ES17.8E3))') dM(nact)
					endif
					write(*,'("   >>>   dT="(1ES17.8E3))') dTM(nact)
					write(*,'("    --->   dT(line)="(1ES17.8E3), " dT(cont)="(1ES17.8E3))') dN2, dN4
					write(*,'("    ->> Te(icell_max2)="(1F14.4)" K", " Tion="(1ES17.8E3)" K")') T(icell_max_2), Tion_ref(nact)
					write(*,'("    ->> Te(icell_max1)="(1F14.4)" K", " Texi="(1ES17.8E3)" K")') T(icell_max), Tex_ref(nact)
					write(*,*) " ------------------------------------------------ "
				enddo
				if (dne /= 0.0_dp) write(*,'("   >>> dne="(1ES17.8E3))') dne

				write(*,'(" <<->> diff="(1ES17.8E3)," old="(1ES17.8E3))') diff, diff_old !at the end of the loop over n_cells
	        	write(*,"('Unconverged cells #'(1I5), ' fraction :'(1F12.3)' %')") size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))), 100.*real(size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))))/real(size(pack(icompute_atomRT,mask=icompute_atomRT>0)))
				write(*,*) " *************************************************************** "
				diff_old = diff
				
				!Not used if the next is not commented out
				
				if ((real(diff) < precision).and.(maxval_cswitch_atoms() == 1.0_dp)) then
           			if (lprevious_converged) then
            	  		lconverged = .true.
           			else
            	  		lprevious_converged = .true.
          	    	endif
        		else !continue to iterate even if n_rayons max is reached ?
           			lprevious_converged = .false.
           			if ((cswitch_enabled).and.(maxval_cswitch_atoms() > 1.0)) then
           				call adjust_cswitch_atoms
					endif
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
        		
				!***********************************************************!
				! ********** timing and checkpointing **********************!
				
				call system_clock(time_end,count_rate=time_tick,count_max=time_max)	
  				if (time_end < time_begin) then
     				time_nlte=real(time_end + (1.0 * time_max)- time_begin)/real(time_tick)
  				else
     				time_nlte=real(time_end - time_begin)/real(time_tick)
  				endif
	
        		if (n_iter <= 4) then
        			time_iteration = time_iteration + time_nlte  * 0.25
        			!if the problem converge in less than 4 iterations pb
        		endif

        		
        		if (lsafe_stop) then

					if ((time_nlte + time_iteration >=  safe_stop_time).and.(n_iter >= 4)) then
						lconverged = .true.
						lprevious_converged = .true.
						call warning("Time limit would be exceeded, leaving...")
						write(*,*) " time limit:", mod(safe_stop_time/60.,60.) ," min"
						write(*,*) " ~<time> etape:", mod(n_iter * time_iteration/60.,60.), ' <time iter>=', mod(time_iteration/60.,60.)," min"
  						write(*,*) " ~<time> etape (cpu):", mod(n_iter * time_iteration * nb_proc/60.,60.), " min"
  						write(*,*) ' time =',mod(time_nlte/60.,60.), " min"
  						lexit_after_nonlte_loop = .true.
  						!lsafe_stop only would leave the code even if the time is below the walltime
						exit step_loop
					endif

				endif
				!***********************************************************!

			end do !while
			write(*,*) "step: ", etape, "Threshold: ", precision!dpops_max_error
			!real values are possible, not needed and would increase the amount
			!of lign of codes.
			if (n_iter >= 4) then
				write(*,*) " ~<time> etape:", mod(n_iter * time_iteration/60.,60.), ' <time iter>=', mod(time_iteration/60.,60.)," min"
  				write(*,*) " ~<time> etape (cpu):", mod(n_iter * time_iteration * nb_proc/60.,60.), " min"
  			endif

		end do step_loop

! -------------------------------- CLEANING ------------------------------------------ !

		if (lNg_acceleration) then
			if (allocated(ngpop)) deallocate(ngpop)
			if (allocated(ng_cur)) deallocate(ng_cur)
		endif

		!if (lelectron_scattering) then
! 		write(*,*) " Writing Jnu cont to ascii file..."
!   		open(unit=20, file=Jnu_File_ascii, status="unknown")
!   		write(20,*) n_cells, Nlambda_cont
!   		do icell=1, n_cells
!   			do la=1, Nlambda_cont
!     			write(20,'(1F12.5,5E20.7E3)') lambda_cont(la), Jnu_cont(la,icell), 0.0, thomson(icell)/(chi_c(la,icell)+chi_c_nlte(la,icell)), 0.0, 0.0
!    			enddo
!   		enddo
!   		close(20)
!   		write(*,*) "done"
  		!endif
		if (allocated(xyz_pos)) then
			open(unit=20,file="xyz_pos.txt",status="unknown")
			write(20,*) 1, Ncentre+1
			write(20,*) 1627, cell_map_i(1627), cell_map_j(1627), cell_map_k(1627)
			do iray=1,Ncentre+1
				write(20,'(*(1E20.7E3))') (xyz_pos(i,iray,1),i=1,3)
			enddo
! 			write(20,'(*(1E20.7E3))') (frac_pos(iray-1,i),i=1,3)
			do iray=2,Ncentre+1
				write(20,'(*(1E20.7E3))') (frac_pos(iray-1,i),i=1,3)
			enddo
			close(20)
			open(unit=20,file="uvw_pos.txt",status="unknown")
			write(20,*) 1, size(xmu), size(xmu)/8
			do imu=1,size(xmu)
				write(20,'(*(1E20.7E3))') (uvw_pos(i,imu,1),i=1,3)
			enddo
			close(20)
		endif
  		
		if (allocated(threeKminusJ)) then
			!only at id_ref until fits file
			write(*,*) " Writing anisotropy to ascii file..."

  			open(unit=20, file="anis_ascii.txt", status="unknown")
  			write(20,*) n_cells, 1!Nlambda
  			do icell=1, n_cells
  				!do la=1, Nlambda
  				do la=id_ref,id_ref
    				write(20,'(1F12.5,5E20.7E3)') lambda(la), threeKminusJ(la,icell), 0.0, 0.0, 0.0, 0.0
   				enddo
  			enddo
  			close(20)
  			write(*,*) "done" 
  			deallocate(threeKminusJ) 			
		endif
		if (allocated(psi_mean)) then
			write(*,*) " Writing <Psi> to ascii file..."	
  			open(unit=20, file="psim_ascii.txt", status="unknown")
  			write(20,*) n_cells, Nlambda
  			do icell=1, n_cells
  				do la=1, Nlambda
    				write(20,'(1F12.5,5E20.7E3)') lambda(la), psi_mean(la,icell), 0.0, 0.0, 0.0, 0.0
   				enddo
  			enddo
  			close(20)
   			write(*,*) "done"  			 			
		endif
		
		if (allocated(convergence_map)) then
			do nact=1, NactiveAtoms
				call write_convergence_map_atom(ActiveAtoms(nact)%ptr_atom, etape_end, convergence_map(:,1:ActiveAtoms(nact)%ptr_atom%Nlevel, nact, :))
			enddo
			if (iterate_ne) call write_convergence_map_electron(etape_end, convergence_map(:,1,NactiveAtoms+1,:))
			deallocate(convergence_map)
		endif
		
		deallocate(dM, dTM, Tex_ref, Tion_ref)
		if (allocated(Jnew)) deallocate(Jnew)
		if (allocated(Jnew_cont)) deallocate(Jnew_cont)
		deallocate(psi, chi_up, chi_down, Uji_down, eta_atoms, n_new)
		deallocate(stream)
	
		!!close(unit=unit_invfile)
		
		if (n_iterate_ne > 0) then
			write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
			call write_electron
			call write_Hminus
		endif

	return
	end subroutine NLTEloop_mali
	

	subroutine init_stellar_disk
		integer :: i_star!, lam
		!read stellar radiation here or compute from mcfsot

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

	function local_stellar_brigthness(N,lambda,i_star,icell_prev,x,y,z,u,v,w)
	! ---------------------------------------------------------------!
	!
	! -------------------------------------------------------------- !
  		integer, intent(in) :: N, i_star, icell_prev
  		real(kind=dp), dimension(N), intent(in) :: lambda
  		real(kind=dp), intent(in) :: u, v, w, x, y, z
  		real(kind=dp), dimension(N) :: local_stellar_brigthness

  		real(kind=dp) :: Tchoc
  		real(kind=dp) :: mu, ulimb, LimbDarkening
  		integer :: ns,la
  		logical :: lintersect = .false.

   		if (etoile(i_star)%T <= 1e-6) then !even with spots
    		local_stellar_brigthness(:) = 0.0_dp !look at init_stellar_disk
    		return !no radiation from the starr
   		endif
   
  		local_stellar_brigthness = 1.0_dp
  		!add X-rays, UV flux etc..
   
		!if (ladd_xrays) then
		!	....
		!endif
   
		!cos(theta) = dot(r,n)/module(r)/module(n)
		if (llimb_darkening) then
     		call ERROR("option for reading limb darkening not implemented")
			mu = abs(x*u + y*v + z*w)/sqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
			if (real(mu)>1d0) then !to avoid perecision error
				write(*,*) "mu=",mu, x, y, z, u, v, w
				call Error(" mu limb > 1!")
			end if
		else
			LimbDarkening = 1.0_dp
		end if

   
   		!better to store a map(1:n_cells) with all Tchoc
   		!and if map(icell_prev) > 0 -> choc
		if ((laccretion_shock).and.(icell_prev<=n_cells)) then
			if (icompute_atomRT(icell_prev)) then
				if (vr(icell_prev) < 0.0_dp) then
					if (Taccretion>0) then
						Tchoc = Taccretion
					else!need a condition to use vtheta or vphi. Or an array that contains for each cell Tshock or zero (only for cell close to the star)
						Tchoc = (1d-3 * masseH * wght_per_H * nHtot(icell_prev)*abs(vr(icell_prev))/sigma * &
							(0.5 * (vr(icell_prev)**2+v_z(icell_prev)**2+vphi(icell_prev)**2)))**0.25
					endif
					lintersect = (Tchoc > etoile(i_star)%T)
				endif
			endif
			if (lintersect) then
				local_stellar_brigthness(:) = local_stellar_brigthness(:) + &
				(exp(hc_k/max(lambda,10.0)/etoile(i_star)%T)-1)/(exp(hc_k/max(lambda,10.0)/Tchoc)-1)
			endif

   		endif


		local_stellar_brigthness(:) = LimbDarkening * local_stellar_brigthness(:)

	return
	end function local_stellar_brigthness
 
!  subroutine calc_stellar_surface_brightness(N,lambda,i_star,icell_prev,x,y,z,u,v,w,gamma)
!  ! ---------------------------------------------------------------!
!   ! Compute the stellar radiation field surface brightness.
!   ! Istar = B(x,y,z;mu) * Stellar spectrum or B * BlackBody
!   ! return gamma, the brightness. For uniform disk, gamma is 1.
!   !
!   ! If there is a shock spot at the surface, gamma returned is :
!   ! gamma = 1 + ratio, such that the radiation from the star Istar is
!   ! Istar = I(photosphere) + Ishock, with Ishock = I(photosphere) * ratio.
!   ! (previously, Istar was Ishock, now it is the sum of the two contrib)
!  ! -------------------------------------------------------------- !
!   integer, intent(in) :: N, i_star, icell_prev
!   real(kind=dp), dimension(N), intent(in) :: lambda
!   real(kind=dp), dimension(N), intent(out) :: gamma
!   real(kind=dp), intent(in) :: u, v, w, x, y, z
!   real(kind=dp) :: energie(N), Tchoc
!   real(kind=dp) :: mu, ulimb, LimbDarkening, surface, HC
!   integer :: ns,la
!   logical :: lintersect = .false.
! 
!    if (etoile(i_star)%T <= 1e-6) then !even with spots
!     gamma(:) = 0.0_dp
!     return !no radiation from the starr
!    endif
!    
!    gamma(:) = 1.0_dp !init
!    					 !if no spots (and no other sources like X rays, UV flux)
!    					 !it is the outvalue
!    
!    !if (ladd_xrays) then
!     !such that Ixray = Iphot * gamma and Istar = Iphot + Ixray = Iphot * (1 + gamma)
! !     gamma(:) = gamma(:) + (exp(hc_k/lambda/etoile(i_star)%T)-1)/(exp(hc_k/lambda/1d6)-1)
! !     where (lambda <= 50.)
! !         gamma(:) = gamma(:) + (exp(hc_k/lambda/etoile(i_star)%T)-1)/(exp(hc_k/lambda/1d6)-1)
! !     end where
!    !endif
!    
!    !cos(theta) = dot(r,n)/module(r)/module(n)
!    mu = abs(x*u + y*v + z*w)/sqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
!    if (real(mu)>1d0) then !to avoid perecision error
!     write(*,*) "mu=",mu, x, y, z, u, v, w
!     call Error(" mu limb > 1!")
!    end if
!    
!    
!   ! Correct with the contrast gamma of a hotter/cooler region if any
! !    call intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect)
! !    !avoid error with infinity for small lambda
! !    if (lintersect) then
! !    		!Means that Ispot = Bp(Tspot) = gamma * Iphot  = Ispot
! !    		!gamma is initialized to one here.
! !    		!The +1 (i.e., the gamma = gamma + ...) means that Istar = Iphot + Ispot = Iphot * (1 + gamma)
! ! 		gamma(:) = gamma(:) + (exp(hc_k/max(lambda,10.0)/etoile(i_star)%T)-1)/(exp(hc_k/max(lambda,10.0)/etoile(i_star)%SurfB(ns)%T)-1)
! !      !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
! !      	!Lambda independent spots, Ts = 2*Tphot means Fspot = 2 * Fphot
! !  		!gamma = gamma + (etoile(i_star)%SurfB(ns)%T - etoile(i_star)%T) / etoile(i_star)%T
! !    end if
!    
!    if ((laccretion_shock).and.(icell_prev<=n_cells)) then
!    	if (icompute_atomRT(icell_prev)) then
!    		if (vr(icell_prev) < 0.0_dp) then
! !    		write(*,*) "Accretion E (K):", (1d-3 * masseH * wght_per_H * nHtot(icell_prev)*abs(vr(icell_prev))/sigma * (0.5 * (vr(icell_prev)**2+v_z(icell_prev)**2+vphi(icell_prev)**2)))**0.25
! !    			lintersect = .true.
!    			if (Taccretion>0) then
!    				Tchoc = Taccretion
!    			else!need a condition to use vtheta or vphi. Or an array that contains for each cell Tshock or zero (only for cell close to the star)
!    				Tchoc = (1d-3 * masseH * wght_per_H * nHtot(icell_prev)*abs(vr(icell_prev))/sigma * (0.5 * (vr(icell_prev)**2+v_z(icell_prev)**2+vphi(icell_prev)**2)))**0.25
!    			endif
! !    			if (Taccretion>0) then
! !    				Tchoc = Taccretion
! !    			else!need a condition to use vtheta or vphi. Or an array that contains for each cell Tshock or zero (only for cell close to the star)
! !    				Tchoc = (1d-3 * masseH * wght_per_H * nHtot(icell_prev)*abs(vr(icell_prev))/sigma * (0.5 * (vr(icell_prev)**2+v_z(icell_prev)**2+vphi(icell_prev)**2)))**0.25
! !    			endif
!    			lintersect = (Tchoc > etoile(i_star)%T)
!    		endif
! 	endif !rho > 0
! 	if (lintersect) then
! ! 		write(*,*) "intersect spot"
! 		gamma(:) = gamma(:) + (exp(hc_k/max(lambda,10.0)/etoile(i_star)%T)-1)/(exp(hc_k/max(lambda,10.0)/Tchoc)-1)
! 	endif
! 
!    endif!cells <= n_cells
! 
! 
!    !Apply Limb darkening
!    if (llimb_darkening) then
!      call ERROR("option for reading limb darkening not implemented")
!    else
!      !ulimb = 0.6
!      LimbDarkening = 1.0_dp!0.4 + 0.6 * mu!1.0_dp
!    end if
!    !Istar(:) = energie(:) * LimbDarkening * gamma(:)
!    gamma(:) = LimbDarkening * gamma(:)
! 
!  return
!  end subroutine calc_stellar_surface_brightness

   subroutine INTEG_RAY_JNU(id,icell_in,x,y,z,u,v,w,iray,labs, kappa_tot, Snu, Istar, Ic)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray to
  ! compute coherent Jnu in the continuum
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(in), dimension(:,:) :: kappa_tot, Snu
  real(kind=dp), intent(in), dimension(:) :: Istar
  real(kind=dp), intent(out) :: Ic(:,:)
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
  previous_cell = icell_in

  !tau(:) = 0.0_dp !go from surface down to the star
  tau_c(:) = 0.0_dp

  Ic(:,id) = 0.0
  
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
       !call calc_stellar_surface_brightness(size(Ic(:,1)),lambda,i_star,previous_cell,x0,y0,z0,u,v,w,LimbD)
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
  integer :: n_rayons_start, n_rayons_max
  integer, parameter :: n_rayons_max2 = n_rayons_start2 * (2**(maxIter-1))
  real :: precision = 1e-1!, parameter
  real, parameter :: lambda_min = 5., lambda_max0 = 100000
  integer :: etape, etape_start, etape_end, iray, n_rayons, n_rayons_old
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lprevious_converged, write_convergence_file

  real :: rand, rand2, rand3, smu, a1, a0, a2, a3
  real, allocatable, dimension(:,:) :: randz
  real(kind=dp) :: dSource, dN
  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dJ, diffs
  real(kind=dp), allocatable :: Jold(:,:), Jnew(:,:), tau_tot(:,:), Jnew_l(:,:), Jold_l(:,:)
  real(kind=dp), allocatable :: Snew(:,:), Kappa_tot(:,:), beta(:,:), Istar(:), Ic(:,:,:), J20_cont(:,:)
  real(kind=dp), allocatable :: lambda_star(:,:), Sth(:,:), Sold(:,:), Sline(:,:), delta_J(:)
                                   
  logical :: labs, l_iterate
  integer :: la, icell, imax, icell_max, icell_max_s, imax_s
  integer :: imu, iphi
  real(kind=dp) :: lambda_max, weight
  
  integer :: iorder, i0_rest, n_iter_accel, iacc
  logical :: ng_rest, accelerated
  real(kind=dp), dimension(:), allocatable :: ng_cur
  real(kind=dp), dimension(:,:), allocatable :: ngJ
  
  	lNg_acceleration = .false.!no acceleration for continuum
  	
	if (lNg_acceleration) then 
		n_iter_accel = 0
		i0_rest = 0
		ng_rest = .false.
		iacc = 0										!+1 or +2
		allocate(ngJ(n_cells * Nlambda_cont, iNg_Norder+2), stat=alloc_status)
		if (alloc_status > 0) then
			call error("Cannot allocate Ng table !")
		endif
	endif
   		
  
  n_rayons_start = Nrays_atom_transfer
  n_rayons_max = n_rayons_start

  write(*,*) "   --> Lambda iterating Jnu with Nlambda ", Nlambda_cont
  precision = dpops_max_error

  write(*,*) "   precision in J is ", precision
  write(*,*) "   n_rayons_max is ", n_rayons_max
  write_convergence_file = .false.
  
  if (allocated(ds)) deallocate(ds)
  if  (laccurate_integ) n_rayons_max = n_rayons_max2
  if (n_rayons_max <= 0) call error("Error n_rayons_max pb !")

  !only one star
  allocate(Istar(Nlambda_cont), delta_J(Nlambda_cont), lambda_star(Nlambda_cont, nb_proc))
  Istar = 0.0_dp; delta_J = 0.0; lambda_star = 0.0_dp

  allocate(Jold(Nlambda_cont, n_cells))!, Jnew(Nlambda_cont, n_cells))
  Jold = 0d0!; Jnew(:,:) = 0.0_dp
  allocate(Sth(Nlambda_cont, n_cells), Snew(Nlambda_cont, n_cells), Sold(Nlambda_cont, n_cells))
  Sth = 0.; Sold = 0.0; Snew = 0.0
  allocate(Kappa_tot(Nlambda_cont, n_cells)); Kappa_tot = 0.
  allocate(beta(Nlambda_cont, n_cells)); beta = 0.

  if (allocated(threeKminusJ)) deallocate(threekMinusJ)
  allocate(threeKminusJ(Nlambda,n_cells),J20_cont(Nlambda_cont, n_cells))


  
  Istar(:) = Istar_cont(:,1)
  Jold = Jnu_cont
  
  !write(*,*) "  -> interpolating contopac on Jnu grid for each frequency.."
  !write(*,*) "       lambda (min, max in nm):", minval(lambda_cont), maxval(lambda_cont)
  do icell=1, n_cells
   if (icompute_atomRT(icell) > 0) then
    
   					
   	kappa_tot(:,icell) = chi_c(:,icell)


   	if (any_nan_infinity_vector(kappa_tot(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in kappa_tot after interpolation.."
   	 write(*,*) "icell=", icell, "kappa=",kappa_tot(:,icell)
   	 write(*,*) "Kc=", chi_c(:,icell)
   	endif
   					
   					
   	beta(:,icell) = thomson(icell) / ( kappa_tot(:,icell) + tiny_dp)

   	
   	if (any_nan_infinity_vector(beta(:,icell)) > 0) then
   	 write(*,*) " Error inconsistant values found in beta after interpolation.."
   	 write(*,*) "icell=", icell, "beta=",beta(:,icell)
   	 write(*,*) "sigma=", thomson(icell)
   	endif   					

   	Sth(:,icell) = eta_c(:,icell) / ( tiny_dp + kappa_tot(:,icell) * (1.-beta(:,icell)))

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
  etape_start = 1
  if (laccurate_integ) then
  	write(*,*) " Using step 3"
  	etape_end = 3
  	allocate(ds(n_rayons_max, nb_proc))
  	ds = 0.0_dp !meters
  else
  	etape_end = 1
  	allocate(ds(1, nb_proc))
  	ds = 0.0_dp !meters
  endif


	if (write_convergence_file ) then
		open(unit=20, file="Jnu_convergence.s", status="unknown")
		write(20,*) "maxIter:", maxIter, " Nlam:", Nlambda, " Ncells:", n_cells
	endif

     do etape=etape_start, etape_end

      if (etape==1) then 
      	!!call error("no step 1 implemented")
      	call compute_angular_integration_weights()
  		lfixed_rays = .true.
  		n_rayons = 1 + (n_rayons_start-1)
  		iray_start = 1
  		lprevious_converged = .false.
		write(*,*) "Step 1: Using ", size(xmu), " rays for Jnu and ", n_rayons, " positions."
		lcell_converged(:) = .false.
  		!only one star
  		allocate(Ic(Nlambda_cont, 1, nb_proc))
  		Ic = 0.0_dp
    	!allocate(xyz_pos(3,n_rayons,2),uvw_pos(3,size(xmu),1))
		if (n_rayons > 1) then
  			stream(1) = init_sprng(gtype,0,1,seed,SPRNG_DEFAULT)
  			allocate(randz(n_rayons-1,3))
  			do iray=1,n_rayons-1
  				do i=1,3
  					randz(iray,i)  = sprng(stream(1))
  				enddo
  				write(*,*) "pos=",iray, (randz(iray,i),i=1,3)
  			enddo
  		endif
      else if (etape==2) then 
  		lfixed_rays = .true.
  		n_rayons = n_rayons_start
  		iray_start = 1
  		lprevious_converged = .false.
		write(*,*) "Step 2: Using ", n_rayons, " rays for Jnu."
		lcell_converged(:) = .false.
		if(allocated(Ic)) deallocate(Ic)
  		allocate(Ic(Nlambda_cont, 1, nb_proc))
  		Ic = 0.0_dp
  	  else if (etape==3) then 
  		write(*,*) "Step 3: max_rayons =  ", n_rayons_max2, " n_rayons_start = ", n_rayons_start2
  		lfixed_rays = .false.
  		n_rayons = n_rayons_start2
  		lprevious_converged = .false.
		lcell_converged(:) = .false.
		precision = 0.1
		if(allocated(Ic)) deallocate(Ic)
  		allocate(Ic(Nlambda_cont, n_rayons_max, nb_proc))
  		Ic = 0.0_dp
	  else 
	  	call error("step unknown")
  	  end if
  	  
		if (write_convergence_file ) write(20,*) "step:", etape, " nrays:", n_rayons

  	  
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

			if (write_convergence_file ) write(20,*) " -> Iteration #", n_iter
            imax = 1
            imax_s = 1
            
 			!$omp parallel &
            !$omp default(none) &
            !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02,iphi,imu) &
            !$omp private(argmt,norme, icell, delta_J, l_iterate,weight) &
            !$omp shared(lambda_cont, lambda_star, Snew, Sold, Sth, Istar, xmu,wmu, xmux, xmuy,randz) &
            !$omp shared(lkeplerian,n_iter,gtype, nb_proc,seed,xyz_pos,uvw_pos,etoile) &
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid, phi_grid, lcell_converged) &
            !$omp shared(n_cells,ds, Jold, Jnu_cont, beta, kappa_tot, Ic,icompute_atomRT,J20_cont) &
            !$omp shared(lfixed_Rays,lnotfixed_Rays,labs,etape,pos_em_cellule)
            !$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !$ id = omp_get_thread_num() + 1

! 				if (lfixed_rays) then
! 					l_iterate = (icompute_atomRT(icell)>0)
! 				else
! 					l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! 				endif
!    				if (l_iterate) then
   				if (icompute_atomRT(icell)>0) then
   				   Jnu_cont(:,icell) = 0.0_dp
           		   Snew(:,icell) = 0.0_dp
           		   lambda_star(:,id) = 0.0_dp 
           		   J20_cont(:,icell) = 0.0_dp
                 		   		            		
           		
      	 			if (etape==2) then

						!ray by ray integration, use iray==1 for quantity
						do iray=iray_start, iray_start-1+n_rayons
					
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

            				
							call INTEG_RAY_JNU(id, icell, x0, y0, z0, u0, v0, w0, 1, labs, kappa_tot, Sold, Istar, Ic(:,1,:))
							
							!LI
		                	Jnu_cont(:,icell) = Jnu_cont(:,icell) + Ic(:,1,id) / n_rayons
		                	!ALI
		                	lambda_star(:,id) = lambda_star(:,id) + (1d0 - exp(-ds(1,id)*kappa_tot(:,icell))) / n_rayons
          		          	J20_cont(:,icell) = J20_cont(:,icell) +  (3.0 * (u0*x0+v0*y0+w0*z0)**2/(x0**2+y0**2+z0**2) - 1.0) * Ic(:,1,id) / n_rayons
      					
      			   		enddo !iray	
		            	
		            else if (etape==3) then

						do iray=iray_start, iray_start-1+n_rayons
					
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

            				
							call INTEG_RAY_JNU(id, icell, x0, y0, z0, u0, v0, w0, iray, labs, kappa_tot, Sold, Istar, Ic(:,iray,:))
                					
      			   		enddo !iray	

		 			    !because we accumulate rays in step 3
						do iray=1, n_rayons
						
							!LI
		                	Jnu_cont(:,icell) = Jnu_cont(:,icell) + Ic(:,iray,id) / n_rayons
		                	!ALI
		                	lambda_star(:,id) = lambda_star(:,id) + (1d0 - exp(-ds(iray,id)*kappa_tot(:,icell))) / n_rayons
		                	J20_cont(:,icell) = J20_cont(:,icell) +  (3.0 * (u0*x0+v0*y0+w0*z0)**2/(x0**2+y0**2+z0**2) - 1.0) * Ic(:,iray,id) / n_rayons

		            	enddo

      			   	else if (etape==1) then
  		         
  		         		do iray = 1,n_rayons
							if (iray==1) then
								x0 = r_grid(icell)*cos(phi_grid(icell))
								y0 = r_grid(icell)*sin(phi_grid(icell))
								z0 = z_grid(icell)
							else
! 						!-> the same for all cells and proc
								call  pos_em_cellule(icell ,randz(iray-1,1),randz(iray-1,2),randz(iray-1,3),x0,y0,z0)
							endif
! 							if (n_iter==1) then
! 								if (icell==1) then
! 									xyz_pos(:,iray,1) = (/x0/etoile(1)%r,y0/etoile(1)%r,z0/etoile(1)%r/)
! 								else if (icell==n_cells) then
! 									xyz_pos(:,iray,2) = (/x0/etoile(1)%r,y0/etoile(1)%r,z0/etoile(1)%r/)
! 								endif
! 							endif


  		         			do imu=1, size(xmu)
  		         				w0 = xmu(imu)
								u0 = xmux(imu); v0 = xmuy(imu)
								weight = wmu(imu) / n_rayons
							
								call INTEG_RAY_JNU(id, icell, x0, y0, z0, u0, v0, w0, 1, labs, kappa_tot, Sold, Istar, Ic(:,1,:))
		
							!LI
		                		Jnu_cont(:,icell) = Jnu_cont(:,icell) + Ic(:,1,id) * weight
		                  		J20_cont(:,icell) = J20_cont(:,icell) + (3.0 * (u0*x0+v0*y0+w0*z0)**2/(x0**2+y0**2+z0**2) - 1.0) * Ic(:,1,id) * weight

		                	!ALI
		                		lambda_star(:,id) = lambda_star(:,id) + (1d0 - exp(-ds(1,id)*kappa_tot(:,icell))) * weight
	
! 								if (n_iter==1 .and. icell==1) uvw_pos(:,imu,1) = (/u0,v0,w0/)

      			   			enddo !imu	
      			   		enddo !pos / iray		         

		            else
						call error("step unknown!")
		 			end if !etape		            

			   		!ALI
			   		delta_J(:) = Jnu_cont(:,icell) - Lambda_star(:,id) * Sold(:,icell)
			   		Snew(:,icell) = (beta(:,icell)*delta_J(:) + (1.0-beta(:,icell))*Sth(:,icell))/(1.0-beta(:,icell)*Lambda_star(:,id))
			   		!LI
        	   		!Snew(:,icell) = Sth(:,icell) * (1.0_dp - beta(:,icell)) + Jnu_cont(:,icell) * beta(:,icell)

      		   endif !icompute_AtomRT

     		end do !icell
        	!$omp end do
        	!$omp end parallel
        	
				!Ng acceleration, testing
	accelerated = .false.
	if ((lNg_acceleration) .and. (n_iter > iNg_Ndelay)) then
		iorder = n_iter - iNg_Ndelay
		if (ng_rest) then
			write(*,'(" -> Acceleration relaxes... "(1I2)" /"(1I2))') iorder-i0_rest, iNg_Nperiod
			if (write_convergence_file ) write(20,'(" -> Acceleration relaxes "(1I2)" /"(1I2))') iorder-i0_rest, iNg_Nperiod
			if (iorder-i0_rest == iNg_Nperiod) ng_rest = .false.
		else
			i0_rest = iorder
			iacc = iacc + 1
            write(*,'(" -> Accumulate solutions... "(1I2)" /"(1I2))') iacc, iNg_Norder+2
			if (write_convergence_file ) write(20,'(" -> Acceleration relaxes "(1I2)" /"(1I2))') iacc, iNg_Norder+2	
			allocate(ng_cur(n_cells * Nlambda_cont))
   							!or flatten2, reform2 ?
   			ng_cur = flatten2(Nlambda_cont, n_cells,Jnu_cont)
   								!has to be parallel in the future
   			accelerated = ng_accelerate(iacc, ng_cur, n_cells * Nlambda_cont, iNg_Norder, ngJ(:,:))
   			if (accelerated) then
            	n_iter_accel = n_iter_accel + 1 !True number of accelerated iter
            	write(*,'("     ++> accelerated iteration #"(1I4))') n_iter_accel
				if (write_convergence_file ) write(20,'("++> accelerated iteration #"(1I4))') n_iter_accel
                ng_rest = .true.
                Jnu_cont(:,:) = reform2(Nlambda_cont, n_cells, ng_cur)
                iacc = 0 
   			endif
   			deallocate(ng_cur)
   		endif
   	endif

            !should be para
        	diff = 0d0
        	dSource = 0.0_dp
  			cell_loop2 : do icell=1, n_cells

! 				if (lfixed_rays) then
! 					l_iterate = (icompute_atomRT(icell)>0)
! 				else
! 					l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! 				endif
!    				if (l_iterate) then
     			if (icompute_atomRT(icell)>0) then
  					
  						dN = maxval(abs(Snew(:,icell) - Sold(:,icell))/Snew(:,icell))
  						dSource = max(dSource, dN)
						
! 						dN = 0.0_dp
						dJ = 0.0_dp
						do la=1, Nlambda_cont
! 							if (Snew(la,icell) > 0) then
! 								dN = max(dN, abs(1.0 - Sold(la,icell) / Snew(la,icell)))
! 								imax = locate(abs(1.0 - Sold(:,icell) / Snew(:,icell)), dN)
! 							endif
						
						 	if (Jnu_cont(la, icell) > 0) then 
						  		dJ = max(dJ,abs(1.-Jold(la,icell)/Jnu_cont(la,icell)))
						  		imax = locate(abs(1.-Jold(:,icell)/Jnu_cont(:,icell)),abs(1.-Jold(la,icell)/Jnu_cont(la,icell)))
							 endif 
						enddo
						
						!if (mod(icell,10)==0)
! 						write(*,'((1I5)" ::> dJ="(1ES14.5E3), " Jmax="(1ES14.5E3), " Jmin="(1ES14.5E3), " beta="(1ES14.5E3))') icell, real(dJ), maxval(Jnu_cont(:,icell)), minval(Jnu_cont(:,icell)), maxval(beta(:,icell))
						if (write_convergence_file ) write(20,'((1I5)" ::> dJ="(1ES14.5E3), " Jmax="(1ES14.5E3), " Jmin="(1ES14.5E3), " beta="(1ES14.5E3))') icell, real(dJ), maxval(Jnu_cont(:,icell)), minval(Jnu_cont(:,icell)), maxval(beta(:,icell))
						if (dJ > diff) then
						  diff = dJ
						  icell_max = icell
						endif
     					lcell_converged(icell) = (real(dJ) < precision)		
! 						if (dN > diff) then
! 						  diff = dN
! 						  icell_max = icell
! 						  dSource = diff
! 						endif
!      					lcell_converged(icell) = (real(dN) < precision)	

     			end if

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
	        write(*,"('# unconverged cells :'(1I5), '('(1F12.3)' %)')") size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))), &
	          100.*real(size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))))/real(size(pack(icompute_atomRT,mask=icompute_atomRT>0)))
			if (write_convergence_file ) write(20,"('# unconverged cells : '(1I5), '('(1F12.3)' %)')") size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))), &
	          100.*real(size(pack(lcell_converged,mask=(lcell_converged.eqv..false.).and.(icompute_atomRT>0))))/real(size(pack(icompute_atomRT,mask=icompute_atomRT>0)))


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

	if (lNg_acceleration) then
		deallocate(ngJ)
		if (allocated(ng_cur)) deallocate(ng_cur)
	endif	
  
  if (.not.lstop_after_jnu) then
    do icell=1, n_cells
     if (icompute_atomRT(icell)>0) then
  	  call bezier2_interp(Nlambda_cont, lambda_cont, Jnu_cont(:,icell), Nlambda, lambda, eta_es(:,icell))
  	  eta_es(:,icell) = eta_es(:,icell) * Thomson(icell) 
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

		if (allocated(threeKminusJ)) then
	
  			open(unit=20, file="anis_ascii.txt", status="unknown")
  			write(20,*) n_cells, Nlambda
  			do icell=1, n_cells
    	  		call bezier2_interp(Nlambda_cont, lambda_cont, J20_cont(:,icell)/Jnu_cont(:,icell), Nlambda, lambda, threeKminusJ(:,icell))
  				do la=1, Nlambda
    				write(20,'(1F12.5,5E20.7E3)') lambda(la), 0.5*threeKminusJ(la,icell), 0.0, 0.0, 0.0, 0.0
   				enddo
  			enddo
  			close(20)
  			
  			deallocate(threeKminusJ)
  			
		endif

!   if (allocated(xyz_pos)) then
!   	open(unit=20,file="xyz_pos.txt",status="unknown")
!   	!first cell then last cell
!   	write(20,*) 2, n_rayons
!   	do iray=1,n_rayons
!   		write(20,*) (xyz_pos(i,iray,1),i=1,3)
!   	enddo
!   	do iray=1,n_rayons
!   		write(20,*) (xyz_pos(i,iray,2),i=1,3)
!   	enddo
!   	close(20)
!   	open(unit=20,file="uvw_pos.txt",status="unknown")
!   	!only one cell to test step1
!   	write(20,*) 1, size(xmu), size(xmu)/8
!   	do imu=1,size(xmu)
!   		write(20,*) (uvw_pos(i,imu,1),i=1,3)
!   	enddo
!   	close(20)
!   endif
  
  
  if (allocated(xmu)) deallocate(xmu,wmu,xmux,xmuy) !in case we call nlte_loop after

  if (allocated(randz)) deallocate(randz)
  if (allocated(J20_cont)) deallocate(J20_cont)
  deallocate(ds, lcell_converged)
  deallocate(Jold, Sth, kappa_tot, beta, Istar, Ic, Sold, lambda_star, Snew, delta_J)

 ! ------------------------------------------------------------------------------------ !
 return
 end subroutine Iterate_Jnu   

	subroutine compute_contribution_functions!_uvw(u,v,w)
	!for one ray at the centre of the image
	!unpolarised case at the moment.
   !not para yet
		real(kind=dp) :: u, v, w !, intent(in)
		integer :: ibin, iaz
		integer :: id, icell0, icell, iray
		real(kind=dp) :: x0,y0,z0,u0,v0,w0
		logical :: lintersect, labs
		
		write(*,*) " allocating space for contribution function for one ray"
		allocate(cntrb_ray(Nlambda, n_cells))
		cntrb_ray = 0.0_dp
		
		ibin = 1; iaz = 1
  		u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin)
	
		id = 1
		labs = .false.
		iray = 1
		write(*,*) "   Computing contribution function at centre of the image..."

		! Ray tracing : on se propage dans l'autre sens
		u0 = -u ; v0 = -v ; w0 = -w
		x0 = 10*Rmax*u; y0 = 10*Rmax*v; z0 = 10*Rmax*w


		call move_to_grid(id,x0,y0,z0,u0,v0,w0,icell0,lintersect)


		if (lintersect) then
			call integ_ray_cntrb(id,icell0,x0,y0,z0,u0,v0,w0,iray,labs)
		else
			write(*,*) "Nothing in this direction, cntrb_ray = 0"
			return
		endif 
		
		call write_contribution_functions_ray
		deallocate(cntrb_ray)  
   
  return
  end subroutine compute_contribution_functions!_uvw
  
  	subroutine integ_ray_cntrb(id,icell_in,x,y,z,u,v,w,iray,labs)
	! ------------------------------------------------------------------------------- !
	!
	! ------------------------------------------------------------------------------- !

		integer, intent(in) :: id, icell_in, iray
		real(kind=dp), intent(in) :: u,v,w
		real(kind=dp), intent(in) :: x,y,z
		logical, intent(in) :: labs
		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
		real(kind=dp), dimension(Nlambda) :: tau, dtau!,etal, chil
		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

		x1=x;y1=y;z1=z
		x0=x;y0=y;z0=z
		next_cell = icell_in
		nbr_cell = 0

		tau(:) = 0.0_dp


		call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)


		infinie : do

			icell = next_cell
			x0=x1 ; y0=y1 ; z0=z1

			if (icell <= n_cells) then
				lcellule_non_vide = (icompute_atomRT(icell) > 0)
				if (icompute_atomRT(icell) < 0) return !-1 if dark
			else
				lcellule_non_vide=.false.
			endif
    

			if (test_exit_grid(icell, x0, y0, z0)) return

			if (lintersect_stars) then
				if (icell == icell_star) return
   			 endif

			nbr_cell = nbr_cell + 1

			previous_cell = 0 
			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
			
			if (lcellule_non_vide) then
				lsubtract_avg = ((nbr_cell == 1).and.labs)

				l_contrib = l_contrib * AU_to_m

				!Total source fonction or set eta to 0 to have only line emissivity
				chi(:,id) = chi0_bb(:,icell)
				eta(:,id) = eta0_bb(:,icell)

				!includes a loop over all bound-bound, passive and active
				call opacity_atom_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,.false.) 
							
				dtau(:) = l_contrib * chi(:,id)
				            
				if (lelectron_scattering) then
					eta(:,id) = eta(:,id) + eta_es(:,icell)
				endif

				!Use etal only ?  here I plot dI/dr integrated over directions
				!S_contrib(:,icell,id) =  eta(:,id) * exp(-tau(:)) !or S * -dtau/dr = S * dtau/l_contrib
				do la=1, Nlambda
					if (tau(la) > 0) cntrb_ray(la,icell) = eta(la,id) * E2(tau(la))!exp(-tau(:))
					!->Flow chart
					!cntrb_ray(la,icell) = ( eta(la,id)/chi(la,id) ) * (1.0_dp - exp(-dtau(la))) * exp(-tau(la))
				enddo
				

				tau = tau + dtau

			end if
		end do infinie

	return
	end subroutine integ_ray_Cntrb

	!building
	subroutine NLTEloop(n_rayons_max,n_rayons_1,n_rayons_2,maxIter,verbose)
	! -------------------------------------------------------- !
	! Descriptor here
	! -------------------------------------------------------- !
#include "sprng_f.h"
		integer, intent(in) :: n_rayons_1, n_rayons_2, maxIter, n_rayons_max
		logical, intent(in) :: verbose
! 		integer, parameter :: max_sub_iter = 10000, max_global_iter = 1000
! 		integer :: etape, etape_start, etape_end, iray, n_rayons, iray_p
! 		integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, la, kr
! 		integer, dimension(nb_proc) :: max_n_iter_loc
! 		logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged
! 		real :: rand, rand2, rand3, precision, precision_sub
! 		real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, smu
! 		real(kind=dp) :: diff, norme, dN, dN1, dJ, lambda_max
! 		real(kind=dp) :: dT, dN2, dN3, diff_old
! 		real(kind=dp), allocatable :: dTM(:), dM(:), Tion_ref(:), Tex_ref(:)
! 		real(kind=dp), allocatable :: Jnew(:,:), Jold(:,:), Jnew_cont(:,:)
! 		logical :: labs, iterate_ne = .false.
! 		logical :: l_iterate
! 		integer :: nact, imax, icell_max, icell_max_2
! 		integer :: icell, ilevel, imu, iphi
! 		character(len=20) :: ne_start_sol = "NE_MODEL"
! 		type (AtomType), pointer :: atom
! 
! 		open(unit=unit_invfile, file=trim(invpop_file), status="unknown")
! 		write(unit_invfile,*) n_cells
! 		open(unit=unit_profiles, file=trim(profiles_file), status="unknown")
!
! 		allocate(gpop_old(NactiveAtoms,Nmaxlevel,n_cells),stat=alloc_status)
! 		if (alloc_status > 0) call error("Allocation error gpop_old")
! 		gpop_old(:,:,:) = 0.0_dp
! 		do nact=1,NactiveAtoms
! 			atom => ActiveAtoms(nact)%ptr_atom
! 			gpop_old(nact,1:atom%Nlevel,:) = atom%n(:,:)
! 			atom => NULL()
! 		enddo	
!   
! 		allocate(n_new(NactiveAtoms,Nmaxlevel,n_cells),stat=alloc_status)
! 		if (alloc_status > 0) call error("Allocation error n_new")
! 		n_new(:,:,:) = 0.0_dp
! 		allocate(dM(Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
! 		allocate(dTM(Nactiveatoms)); dM=0d0 !keep tracks of Tex for all cells for each atom
! 		allocate(Jnew(Nlambda, n_cells)); Jnew = 0.0
! 		allocate(Jold(Nlambda, n_cells)); Jold = 0.0
! 		allocate(Jnew_cont(Nlambda_cont, n_cells)); Jnew_cont = 0.0
! 		allocate(Tex_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tex for all cells for each line of each atom
! 		allocate(Tion_ref(Nactiveatoms)); dM=0d0 !keep tracks of max Tion for all cells for each cont of each atom
! 		diff_old = 0.0_dp
! 
! 
! 		labs = .true. !to have ds at cell icell = eval_operator
! 		id = 1
! 		etape_start = 2
! 		if (laccurate_integ) then
! 			write(*,*) " Using step 3 with ", n_rayons_max, " nrays max"
! 			etape_end = 3
! 		else
! 			etape_end = 2
! 		endif
!  
! 
! 		iterate_ne = (n_iterate_ne>0)
! 		if (iterate_ne) then
! 			write(*,*) " before iterate ne I need te recompute gij for continua !!"
! 		endif
! 
! 		do etape=etape_start, etape_end
! 
! 			if (etape==1) then
! 				write(*,*) " step 1 not implemented"
! 				stop
! 			else if (etape==2) then 
! 				lfixed_rays = .true.
! 				n_rayons = n_rayons_1
! 				iray_start = 1
! 				lprevious_converged = .false.
! 				lcell_converged(:) = .false.
! 				precision = dpops_max_error
! 				precision_sub = dpops_sub_max_error
! 			else if (etape==3) then
! 		  		lfixed_rays = .false.
!   				n_rayons = n_rayons_2
!   				lprevious_converged = .false.
! 				lcell_converged(:) = .false.
! 				precision = 1e-1 !dpops_max_error
! 				precision_sub = dpops_sub_max_error
! 			else
! 				call ERROR("etape unkown")
! 			end if
!   	  
! 			write(*,*)  "-> Using ", n_rayons, ' rays for step', etape
! 
! 			lnotfixed_rays = .not.lfixed_rays
! 			lconverged = .false.
! 			n_iter = 0
! 
! 			do while (.not.lconverged)
! 
! 				n_iter = n_iter + 1
! 				!if (n_iter > max_global_iter) exit !change step !maxIter
! 				write(*,*) " -> Iteration #", n_iter, " Step #", etape
! 
! 				if (lfixed_rays) then
! 					stream = 0.0
! 					do i=1,nb_proc
! 						stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
! 					end do		
! 				end if
! 
! 				max_n_iter_loc = 0
! 
!         	
! 				!$omp parallel &
! 				!$omp default(none) &
! 				!$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dM, dN, dN1,iray_p,imu,iphi, smu)&
! 				!$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, nact, atom, l_iterate) &
! 				!$omp shared(icompute_atomRT, dpops_sub_max_error, verbose,lkeplerian,n_iter,precision_sub,threeKminusJ,stest_tot) &
! 				!$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,lcell_converged,loutput_rates,xyz_pos,uvw_pos, gtype, seed, nb_proc) &
! 				!$omp shared(n_cells, gpop_old,lforce_lte, integ_ray_line, Itot, Icont, Jnu_cont, eta_es) &
! 				!$omp shared(Jnew, Jnew_cont, Jold, lelectron_scattering,chi0_bb, etA0_bb, activeatoms,n_new) &
! 				!$omp shared(nHmin, chi_c, chi_c_nlte, eta_c, eta_c_nlte, ds, Rij_all, Rji_all, Nmaxtr, Gammaij_all, Nmaxlevel) &
! 				!$omp shared(lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape,pos_em_cellule,Nactiveatoms)
! 				!$omp do schedule(static,1)
! 				do icell=1, n_cells
! 					!$ id = omp_get_thread_num() + 1
! ! 					if (lfixed_rays) then
!    						l_iterate = (icompute_atomRT(icell)>0)
! ! 					else
! !    						l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! ! 					endif				
!    					if (l_iterate) then
!    						Jnew(:,icell) = 0.0
!    						Jnew_cont(:,icell) = 0.0
!    						threeKminusJ(:,icell) = 0.0
! 					!!if (atmos%icompute_atomRT(icell)>0) then
! 
! 						call fill_Collision_matrix(id, icell) !computes = C(i,j) before computing rates
! 						call initGamma(id) !init Gamma to C and init radiative rates
!                     
! 						do iray=iray_start, iray_start-1+n_rayons
! 						
! 							if ( (etape==2).or.(etape==3) ) then
!                    	    ! Position aleatoire dans la cellule
! 
! 								rand  = sprng(stream(id))
! 								rand2 = sprng(stream(id))
! 								rand3 = sprng(stream(id))
!  
! 								call  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)
! 
!                         ! Direction de propagation aleatoire
! 								rand = sprng(stream(id))
! 								W0 = 2.0_dp * rand - 1.0_dp !nz
! 								W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
! 								SRW02 = sqrt(W02)
! 								rand = sprng(stream(id))
! 								ARGMT = PI * (2.0_dp * rand - 1.0_dp)
! 								U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
! 								V0 = SRW02 * sin(ARGMT) !ny = sin(theta)*sin(phi)
! 								
! 								call integ_ray_line(id, icell, x0, y0, z0, u0, v0, w0, iray, labs)
! 								Jnew(:,icell) = Jnew(:,icell) + Itot(:,iray,id)
! 								Jnew_cont(:,icell) = Jnew_cont(:,icell) + Icont(:,id)
! 								threeKminusJ(:,icell) = threeKminusJ(:,icell) + (3.0 * (u0*x0+v0*y0+w0*z0)**2/(x0**2+y0**2+z0**2) - 1.0) * Itot(:,iray,id)
! 
! 							elseif (etape==1) then	
! 								stop
! 							end if !etape
! 
! 						enddo !iray
! 						Jnew(:,icell) = Jnew(:,icell) / n_rayons
! 						Jnew_cont(:,icell) = Jnew_cont(:,icell) / n_rayons
! 						threeKminusJ(:,icell) = threeKminusJ(:,icell) / n_rayons
! 						!!if (lelectron_scattering) then
! 						!!always allocated if nlte ?, need to change that but remember we need
! 						!! in the nlte loop, to write the full mean intensity to fits file.
! 							eta_es(:,icell) = Jnew(:,icell) * thomson(icell)
! 						!!endif
! 		
! 						if (loutput_Rates) &
! 							call store_radiative_rates(id, icell, n_rayons, Nmaxtr, Rij_all(:,:,icell), Rji_all(:,:,icell), Jnew(:,icell),.false.)
! 
! 						if (lforce_lte) then
! 							call update_populations(id, icell, diff, .false., n_iter)
! 							lconverged_loc = .true.
! 						else
! 							lconverged_loc = .false.
! 							do iray=1,n_rayons !-> including n_new
! 								call calc_rates(id, icell, iray, n_rayons)
! 							enddo
! 							call calc_rate_matrix(id, icell, lforce_lte)
! 						endif
! 
! 						n_iter_loc = 0
!      				!!only iterate on cells which are not converged
! 						do while (.not.lconverged_loc)
! 							n_iter_loc = n_iter_loc + 1
! 
! 							if ((mod(n_iter_loc,max(1,1000))==0)) then
! 								write(*,*) " -- Global iteration #", n_iter, " step #", etape, "eps: ", dpops_sub_max_error
! 							endif!(mod(n_iter_loc,max(1,max_sub_iter/10))==0)
! 							call update_Populations(id, icell, diff,(mod(n_iter_loc,max(1,1000))==0),n_iter_loc)
! 							if ((mod(n_iter_loc,max(1,1000))==0)) write(*,*) " --------------------------------------------------------------- "
! 
! 
!                         !for this icell (id) and local iteration
! 							if (diff < precision_sub) then
! 								lconverged_loc = .true.
!      					    !write(*,*) "converged", id, icell, "dpops(sub) = ", diff
! 							else
! !need to include cont nlte
! 								call initGamma(id)
! ! 								call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
! 								do iray=1,n_rayons
! 									call calc_rates(id, icell, iray, n_rayons)
! 								enddo
! 								
! 								call calc_rate_matrix(id, icell, lforce_lte)
! 
! 							end if
! 
! 						end do !local sub iteration
! 						if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
! 						if (loutput_Rates) &
! 							call store_rate_matrices(id,icell,Nmaxlevel,Gammaij_all(:,:,:,icell))
! 					end if !icompute_atomRT
! 				end do !icell
! 				!$omp end do
! 				!$omp end parallel
! 
!      		!Global convergence Tests
! 
!             !should be para
! 				dM(:) = 0.0
! 				diff = 0.0
! 				dJ = 0.0
! 				dT = 0.0
! 				dTM(:) = 0.0
! 				Tex_ref(:) = 0.0
! 				Tion_ref(:) = 0.0
!    				icell_max = 1
!    				icell_max_2 = 1
! 				cell_loop2 : do icell=1,n_cells
! ! 					if (lfixed_rays) then
!    						l_iterate = (icompute_atomRT(icell)>0)
! ! 					else
! !    						l_iterate = (icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
! ! 					endif				
!    					if (l_iterate) then
! 						dN = 0.0 !for all levels of all atoms of this cell
! 						dN2 = 0.0
! 						do nact=1,NactiveAtoms
! 							atom => ActiveAtoms(nact)%ptr_atom
! 
! 							do ilevel=1,atom%Nlevel
! 								!if (atom%n(ilevel,icell) > 0) then
! 								dN1 = abs(1d0-atom%n(ilevel,icell)/n_new(nact,ilevel,icell))
! ! 								dN1 = abs(1d0-gpop_old(nact,ilevel,icell)/atom%n(ilevel,icell))
!      				    		!dN1 = abs(atom%n(ilevel,icell)-gpop_old(nact,ilevel,icell))/(gpop_old(nact,ilevel,icell)+tiny_dp)
! 								dN = max(dN1, dN) !compare with all atoms and levels
!      				    						  !for this cell
! 								dM(nact) = max(dM(nact), dN1) !compare for one atom
!      				    		 							   !for all cells
! 								!endif
! 								!!write(*,*) "pops", dN1,gpop_old(nact,ilevel,icell),atom%n(ilevel,icell)
! 							end do !over ilevel
! 							do kr=1, atom%Nline
! 							
! 								dN3 = abs(1.0 - Tex_old(nact, kr, icell) / ( atom%lines(kr)%Tex(icell) + tiny_dp ) )
! 								!dN3 = abs(atom%lines(kr)%Tex(icell) - Tex_old(nact, kr, icell)) / ( Tex_old(nact, kr, icell) + tiny_dp ) 
! 								dN2 = max(dN3, dN2)
! 								!!write(*,*) "line", icell, T(icell), Tex_old(nact, kr, icell),atom%lines(kr)%Tex(icell), dN2, dN3
! 								!dTM(nact) = max(dTM(nact), dN3)
! 								if (dN3 >= dTM(nact)) then
! 									dTM(nact) = dN3
! 									Tex_ref(nact) = atom%lines(kr)%Tex(icell)
! 									icell_max = icell
! 								endif
! 							enddo
! 							
! 							do kr=1, atom%Ncont
! 								dN3 = abs(1.0 - Tex_old(nact, kr+atom%Nline, icell) / ( atom%continua(kr)%Tex(icell) + tiny_dp ))
! 								!dN3 = abs(atom%continua(kr)%Tex(icell) - Tex_old(nact, kr+atom%Nline, icell)) / ( Tex_old(nact, kr+atom%Nline, icell) + tiny_dp )
! 								dN2 = max(dN3, dN2)
! 								!!write(*,*) "cont", kr, " icell=",icell, " T=", T(icell), " Told=",Tex_old(nact, kr+atom%Nline, icell)," Tex=",atom%continua(kr)%Tex(icell), dN2, dN3
! 								!dTM(nact) = max(dTM(nact), dN3)
! 								if (dN3 >= dTM(nact)) then
! 									dTM(nact) = dN3
! 									Tion_ref(nact) = atom%continua(kr)%Tex(icell)
! 									icell_max_2 = icell
! 								endif
! 							enddo
! 							atom => NULL()
! 						end do !over atoms
! 						!compare for all atoms and all cells
! 						diff = max(diff, dN) ! pops
! ! 						diff = max(diff, dN2) ! Tex
! 						
! 						dN1 = abs(1d0 - maxval(Jold(:,icell)/(tiny_dp + Jnew(:,icell))))
! 						!dN1 = maxval(abs(Jold(:,icell)-Jnew(:,icell))/(1d-30 + Jold(:,icell)))
! 						if (dN1 > dJ) then
! 							dJ = dN1
! 							imax = locate(abs(1d0 - Jold(:,icell)/(tiny_dp + Jnew(:,icell))), dJ)
! 							!imax = locate(abs(Jold(:,icell)-Jnew(:,icell))/(tiny_dp + Jold(:,icell)),dJ)
! 							lambda_max = lambda(imax)
! 						endif
! 						!replace old value by new
! 						Jold(:,icell) = Jnew(:,icell)
! 						Jnu_cont(:,icell) = Jnew_cont(:,icell)
! 
! 						lcell_converged(icell) = (real(diff) < precision) !(real(diff) < dpops_max_error)
! 						
! 						!Re init for next iteration if any
! 						do nact=1, NactiveAtoms
! 							atom => ActiveAtoms(nact)%ptr_atom
! 							!!gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
! 							atom%n(:,icell) = n_new(nact,1:atom%Nlevel,icell)
! 							do kr=1,atom%Nline
! 								Tex_old(nact, kr, icell) = atom%lines(kr)%Tex(icell)
! 							enddo
! 							do kr=1, atom%Ncont
! 								Tex_old(nact, kr+Atom%Nline,icell) = atom%continua(kr)%Tex(icell)
! 							enddo
! 							atom => NULL()
!              				!!Recompute Profiles if needed here
! 						end do
! 						call NLTE_bound_free(icell)
! 						!because of opac nlte changed
! 						call interp_background_opacity(icell, chi0_bb(:,icell), eta0_bb(:,icell))
! 						!end init
! 						
! 					end if
! 				end do cell_loop2 !icell
! 
! 				write(*,'(" -> "(1I10)" sub-iterations")') maxval(max_n_iter_loc)
! 				write(*,'(" -> icell_max1 #"(1I6)," icell_max2 #"(1I6))') icell_max, icell_max_2
! 				write(*,'(" -> dJ="(1ES14.5E3)" @"(1F14.4)" nm")') dJ, lambda_max !at the end of the loop over n_cells
! 				write(*,*) " ------------------------------------------------ "
! 				do nact=1,NactiveAtoms
! 					write(*,'("             Atom "(1A2))') ActiveAtoms(nact)%ptr_atom%ID
! 					write(*,'("   >>> dpop="(1ES17.8E3))') dM(nact)
! 					write(*,'("   >>>   dT="(1ES17.8E3))') dTM(nact)
! 					write(*,'("   >>> Te(icell_max2)="(1F14.4)" K", " Tion="(1F14.4)" K")') T(icell_max_2), Tion_ref(nact)
! 					write(*,'("   >>> Te(icell_max1)="(1F14.4)" K", " Texi="(1F14.4)" K")') T(icell_max), Tex_ref(nact)
! 					!if (Tion_ref(nact) /= 0.0)
! 					!if (Tex_ref(nact) /= 0.0)
! 					write(*,*) " ------------------------------------------------ "
! 				enddo
! 				write(*,'(" <<->> diff="(1ES17.8E3)," old="(1ES17.8E3))') diff, diff_old !at the end of the loop over n_cells
! 	        	write(*,"('Unconverged cells #'(1I5), ' fraction :'(1F12.3)' %')") size(pack(lcell_converged,mask=lcell_converged.eqv..false.)), 100.*real(size(pack(lcell_converged,mask=lcell_converged.eqv..false.)))/real(n_cells)
! 				write(*,*) " *************************************************************** "
! 				diff_old = diff
! 				!Not used if the next is not commented out
! ! 				lconverged = (real(diff) < dpops_max_error) !global convergence for all iterations
! 
! 				if (real(diff) < precision) then
!            			if (lprevious_converged) then
!             	  		lconverged = .true.
!            			else
!             	  		lprevious_converged = .true.
!           	    	endif
!         		else !continue to iterate even if n_rayons max is reached ?
!            			lprevious_converged = .false.
!            			if (.not.lfixed_rays) then
!               			n_rayons = n_rayons * 2
!               			write(*,*) ' -- Increasing number of rays :', n_rayons
!              			if (n_rayons > n_rayons_max) then
!               				if (n_iter >= maxIter) then
!              		 			write(*,*) "Warning : not enough rays to converge !!"
!                  				lconverged = .true.
!               				end if
!               			end if
!           	   		end if
!         		end if
!         									
! 			!Not optimal, but fast relatively
! 			!must be parra
! 				if (iterate_ne .and. (mod(n_iter,n_iterate_ne)==0))  then
! 				write(*,*) "Recompute lte opacities"
! 				write(*,*) " not yet"
! 					stop
! 	
! 				end if
! 
! 
! 			end do !while
! 			write(*,*) "step: ", etape, "Threshold: ", precision!dpops_max_error
! 
! 		end do !over etapes
! 
! ! -------------------------------- CLEANING ------------------------------------------ !
! 
! 		!Always written in nlte
! 		!if (lelectron_scattering) then
!   		open(unit=20, file=Jnu_File_ascii, status="unknown")
!   		write(20,*) n_cells, Nlambda_cont
!   		do icell=1, n_cells
!   			do la=1, Nlambda_cont
!     			write(20,'(1F12.5,5E20.7E3)') lambda_cont(la), Jnu_cont(la,icell), 0.0, thomson(icell)/(chi_c(la,icell)+chi_c_nlte(la,icell)), 0.0, 0.0
!    			enddo
!   		enddo
!   		close(20)
!   		!endif
!   	
! 		if (allocated(threeKminusJ)) then
! 	
!   			open(unit=20, file="anis_ascii.txt", status="unknown")
!   			write(20,*) n_cells, Nlambda
!   			do icell=1, n_cells
!   				do la=1, Nlambda
!     				write(20,'(1F12.5,5E20.7E3)') lambda(la), 0.5*threeKminusJ(la,icell)/Jnew(la,icell), 0.0, 0.0, 0.0, 0.0
!    				enddo
!   			enddo
!   			close(20)
!   			
! 		endif
! 
! 		deallocate(dM, dTM, Tex_ref, Tion_ref)
! 		if (allocated(Jnew)) deallocate(Jnew)
! 		if (allocated(Jnew_cont)) deallocate(Jnew_cont)
! 		if (allocated(Jold)) deallocate(Jold)
! 		if (allocated(n_new)) deallocate(n_new)
! 	
! 		close(unit=unit_invfile)
! 		close(unit=unit_profiles)
	return
	end subroutine NLTEloop
! ------------------------------------------------------------------------------------ !  
  
end module atom_transfer
!   	subroutine integrate_i_icell(id,icell_in,x,y,z,u,v,w,iray,labs, I, Ic)
! 	! ------------------------------------------------------------------------------- !
! 	!
! 	! ------------------------------------------------------------------------------- !
! 
! 		integer, intent(in) :: id, icell_in, iray
! 		real(kind=dp), intent(in) :: u,v,w
! 		real(kind=dp), intent(in) :: x,y,z
! 		logical, intent(in) :: labs
! 		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
! 		real(kind=dp), dimension(Nlambda, n_cells), intent(out) :: I
! 		real(kind=dp), dimension(Nlambda_cont, n_cells), intent(out) :: Ic
! 		real(kind=dp), dimension(Nlambda) :: Snu, tau, dtau
! 		real(kind=dp), dimension(Nlambda_cont) :: Snu_c, dtau_c, tau_c
! 		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la, icell_p
! 		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars
! 
! 		x1=x;y1=y;z1=z
! 		x0=x;y0=y;z0=z
! 		next_cell = icell_in
! 		nbr_cell = 0
! 
! 		tau(:) = 0.0_dp
! 		tau_c(:) = 0.0_dp
! 
! 		I(:,:) = 0d0
! 		Ic(:,:) = 0d0
! 		icell_p = icell_in
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
!     !write(*,*) "Boucle infinie, icell=", icell
! 
! 			if (icell <= n_cells) then
! 				lcellule_non_vide = (icompute_atomRT(icell) > 0)
! 				if (icompute_atomRT(icell) < 0) return !-1 if dark
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
! 					!can be done better
!       				Icont(:,iray,id) =  Icont(:,iray,id) + exp(-tau_c) * Istar_cont(:,i_star) * local_stellar_brigthness(Nlambda_cont,lambda_cont,i_star, icell_prev,x0, y0, z0, u,v,w)
!					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * local_stellar_brigthness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w)
!        				return
!       			end if
!    			 endif
! 
! 			nbr_cell = nbr_cell + 1
! 
!     ! Calcul longeur de vol et profondeur optique dans la cellule
! 			previous_cell = 0 ! unused, just for Voronoi
! 			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
! 			
!     !count opacity only if the cell is filled, else go to next cell
! 			if (lcellule_non_vide) then
! 				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
!      ! opacities in m^-1
! 
! 				l_contrib = l_contrib * AU_to_m !l_contrib in m
! 
! 				if ((nbr_cell == 1).and.labs) then
! 					ds(iray,id) = l * AU_to_m
! 				endif
!      									
! 				!total bound-bound + bound-free + background opacities lte + nlte
! 				chi(:,id) = chi0_bb(:,icell)
! 				eta(:,id) = eta0_bb(:,icell)
! 				!includes a loop over all bound-bound, passive and active
! 				call opacity_atom_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,( (nbr_cell==1).and.labs ) ) 
! 
! 				dtau(:)   = l_contrib * chi(:,id)
! 				            
! 				if (lelectron_scattering) then
! 					Snu = ( eta(:,id) + eta_es(:,icell) ) / ( chi(:,id) + tiny_dp )
! 					Snu_c = eta_c(:,icell) + thomson(icell) * Jnu_cont(:,icell)
! 				else
! 					Snu = eta(:,id) / ( chi(:,id) + tiny_dp )
! 					Snu_c = eta_c(:,icell)
! 				endif
! 				
! 				if (Nactiveatoms > 0) then
! 					Snu_c = ( Snu_c + eta_c_nlte(:,icell) ) / ( chi_c(:,icell) + chi_c_nlte(:,icell) + tiny_dp)
! 					dtau_c(:) = l_contrib * (chi_c(:,icell) + chi_c_nlte(:,icell))
! 				else
! 					Snu_c = Snu_c / (chi_c(:,icell) + tiny_dp)
! 					dtau_c(:) = l_contrib * chi_c(:,icell)
! 				endif
! 
! 
! 				I(:,icell) = I(:,icell_p) + exp(-tau(:)) * (1.0_dp - exp(-dtau(:))) * Snu(:)
! 				Ic(:,icell) = Ic(:,icell_p) + exp(-tau_c(:)) * (1.0_dp - exp(-dtau_c(:))) * Snu_c(:)
! 				icell_p = icell
! 
! 				tau = tau + dtau
! 				tau_c = tau_c + dtau_c
! 
! 			end if  ! lcellule_non_vide
! 		end do infinie
! 
! 	return
! 	end subroutine integrate_i_icell
!  
!   subroutine compute_Imu()
!    use utils, only: Gauss_Legendre_quadrature, span
!    use constantes, only : Au_to_Rsun
!    
!    integer, parameter :: Nmu = 500
!    real(kind=dp) :: u, v, w, uvw(3), x(3), y(3)
!    integer :: id, icell0, icell, i, j, iray, alloc_status
!    real(kind=dp), dimension(:), allocatable :: weight_mu, cos_theta, p
!    integer, parameter :: maxSubPixels = 32
!    real(kind=dp) :: x0,y0,z0,u0,v0,w0, r0, r1, rr, phi
!    real(kind=dp), dimension(:,:,:), allocatable :: Imu, Imuc
!    real(kind=dp):: normF
!    logical :: lintersect, labs
! 
!    allocate(weight_mu(Nmu), cos_theta(Nmu), p(Nmu), stat=alloc_status)
!    if (alloc_status > 0) call error(" Allocation error cos_theta(Nmu)")
! !    allocate(Imu(Nlambda,Nmu,n_cells), Imuc(Nlambda_cont,Nmu,n_cells), stat=alloc_status)
! !    if (alloc_status > 0) call error(" Allocation error Imu, Imuc")
!    allocate(Imu(Nlambda,Nmu,1), Imuc(Nlambda_cont,Nmu,1), stat=alloc_status)
!    if (alloc_status > 0) call error(" Allocation error Imu, Imuc")
! 
! 
!    Imu(:,:,:) = 0.0_dp
!    Imuc(:,:,:) = 0.0_dp
!    labs = .false.
!    id = 1
!    iray = 1
!    phi = pi
! 
!    !look parallel to z
!    u = 0.0_dp
!    v = -1.745d-22!0.0_dp
!    w = 1.0_dp 
! 
!    uvw = (/u,v,w/) !vector position
!    x = (/1.0_dp,0.0_dp,0.0_dp/)
!    y = -cross_product(x, uvw)
! 
!    ! Ray tracing : on se propage dans l'autre sens
!    u0 = -u ; v0 = -v ; w0 = -w
!    
! 
!    !r0 = 1d-5
!    !r1 = Rmax * 1.005
!    !q = span(real(r0),real(r1), Nmu)
!    
!    call gauss_legendre_quadrature(0.0_dp, 1.0_dp, Nmu, cos_theta, weight_mu)
! 
!   !!!$omp parallel &
!   !!!$omp default(none) &
!   !!!$omp private(j,i,id,u0,v0,w0,lintersect,icell0) &
!   !!!$omp shared(Imu, Imuc, u,v,w,Rmax,iray,labs,Itot, Icont)
!    do j=1,Nmu
!    !!!$ id = omp_get_thread_num() + 1
! 
!     rr = Rmax * sqrt(1.0 - cos_theta(j)**2)
!     !!write(*,*) "rr=", rr, " q=", q(j)
!    
!     x0 = 10.0*Rmax*u + rr * sin(phi) * x(1) + rr * cos(phi) * y(1)
!     y0 = 10.0*Rmax*v + rr * sin(phi) * x(2) + rr * cos(phi) * y(2)
!     z0 = 10.0*Rmax*w + rr * sin(phi) * x(3) + rr * cos(phi) * y(3)
!     
!     !!write(*,*) "x'=",rr * sin(phi) * x(1) + rr * cos(phi) * y(1)
!     !!write(*,*) "y'=",rr * sin(phi) * x(2) + rr * cos(phi) * y(2)
!     !!write(*,*) "z'=",rr * sin(phi) * x(3) + rr * cos(phi) * y(3)
! 
! 
!     call move_to_grid(id,x0,y0,z0,u0,v0,w0,icell0,lintersect)
!     
!     !cos_theta(j)  = abs(x0*u + y0*v + z0*w) / sqrt(z0*z0 + x0*x0 + y0*y0)
!     p(j) = rr*AU_to_rsun
!     !!write(*,*) "mu=", cos_theta(j), abs(x0*u + y0*v + z0*w) / sqrt(z0*z0 + x0*x0 + y0*y0)
! 
! 	!!write(*,*) "x0=",x0, " y0=",y0, " z0=",z0
! 
! 	if (cos_theta(j) < 0.05) then
! 		if( (icell0 > n_cells).or.(icell0 < 1) ) then
! 			call warning("for this cos(theta) the ray is put in the star! check prec grid")
! 			write(*,*) "cos(theta)=",cos_theta(j), " icell=", icell0, " intersect?:", lintersect
! 			cycle !because of prec grid ??
! 		endif
! 	endif
! 
!      if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
!      	call integ_ray_line(id, icell0, x0,y0,z0,u0,v0,w0,iray,labs)
!         Imu(:,j,1) = Itot(:,iray,id)
!         Imuc(:,j,1) = Icont(:,iray,id)
! ! 		call integrate_i_icell(id, icell0, x0,y0,z0,u0,v0,w0,iray,labs, Imu(:,j,:), Imuc(:,j,:))
! 
!      endif
!    enddo
!    !!!$omp end do
!    !!!$omp end parallel
! 
!    open(unit=14,file="Imu.s",status="unknown")
!    open(unit=15,file="Imuc.s",status="unknown")
!    write(14,*) Nlambda, Nmu, 1
!    write(15,*) Nlambda_cont, Nmu, 1
!    write(14,'(*(E20.7E3))') (p(j), j=1,Nmu)
!    write(15,'(*(E20.7E3))') (p(j), j=1,Nmu)
!    write(14,'(*(E20.7E3))') (cos_theta(j), j=1,Nmu)
!    write(15,'(*(E20.7E3))') (cos_theta(j), j=1,Nmu)
!    do i=1,Nlambda
!     write(14,'(1F12.5, *(E20.7E3))') lambda(i), (Imu(i,j,1),  j=1,Nmu)
!    enddo
!    do i=1,Nlambda_cont
!     write(15,'(1F12.5, *(E20.7E3))') lambda_cont(i), (Imuc(i,j,1),  j=1,Nmu)
!    enddo
!    close(14)
!    close(15)
! 
! !    open(unit=14,file="Imu.s",status="unknown")
! !    open(unit=15,file="Imuc.s",status="unknown")
! !    write(14,*) Nlambda, Nmu, n_cells
! !    write(15,*) Nlambda_cont, Nmu, n_cells
! !    write(14,'(*(E20.7E3))') (p(j), j=1,Nmu)
! !    write(15,'(*(E20.7E3))') (p(j), j=1,Nmu)
! !    write(14,'(*(E20.7E3))') (cos_theta(j), j=1,Nmu)
! !    write(15,'(*(E20.7E3))') (cos_theta(j), j=1,Nmu)
! !    do icell=1, n_cells
! !    	do i=1,Nlambda
! !     	write(14,'(1F12.5, *(E20.7E3))') lambda(i), (Imu(i,j,icell),  j=1,Nmu)
! !   	enddo
! !    enddo
! !    do icell=1, n_cells
! !    	do i=1,Nlambda_cont
! !     	write(15,'(1F12.5, *(E20.7E3))') lambda_cont(i), (Imuc(i,j,icell),  j=1,Nmu)
! !    	enddo
! !    enddo
! !    close(14)
! !    close(15)
!    
!    deallocate(weight_mu, cos_theta, p, Imu, Imuc)
!    
!   return
!   end subroutine compute_Imu


!   	subroutine integ_ray_line_i_new(id,icell_in,x,y,z,u,v,w,iray,labs)
! 	! ------------------------------------------------------------------------------- !
! 	!
! 	! ------------------------------------------------------------------------------- !
! 
! 		integer, intent(in) :: id, icell_in, iray
! 		real(kind=dp), intent(in) :: u,v,w
! 		real(kind=dp), intent(in) :: x,y,z
! 		logical, intent(in) :: labs
! 		real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, edt, et
! 		real(kind=dp), dimension(Nlambda) :: tau
! 		real(kind=dp) :: Snu, dtau, Snu_c, dtau_c
! 		real(kind=dp), dimension(Nlambda_cont) :: tau_c
! 		integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
! 		logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars
! 
! 		x1=x;y1=y;z1=z
! 		x0=x;y0=y;z0=z
! 		next_cell = icell_in
! 		nbr_cell = 0
! 
! 		tau(:) = 0.0_dp
! 		tau_c(:) = 0.0_dp
! 
! 		Itot(:,iray,id) = 0d0
! 		Icont(:,iray,id) = 0d0
! 		
! 
! !!write(*,*) "id=",id
! !!write(*,*) "****"
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
!     !write(*,*) "Boucle infinie, icell=", icell
! 
! 			if (icell <= n_cells) then
! 				lcellule_non_vide = (icompute_atomRT(icell) > 0)
! 				if (icompute_atomRT(icell) < 0) return !-1 if dark
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
!       				Icont(:,iray,id) =  Icont(:,iray,id) + exp(-tau_c) * Istar_cont(:,i_star) * local_stellar_brigthness(Nlambda_cont,lambda_cont,i_star, icell_prev,x0, y0, z0, u,v,w)
!					Itot(:,iray,id) =  Itot(:,iray,id) + exp(-tau) * Istar_tot(:,i_star) * local_stellar_brigthness(Nlambda,lambda,i_star, icell_prev, x0, y0, z0, u,v,w)
!        				return
!       			end if
!    			 endif
! 
! 			nbr_cell = nbr_cell + 1
! 
!     ! Calcul longeur de vol et profondeur optique dans la cellule
! 			previous_cell = 0 ! unused, just for Voronoi
! 			call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
! 			
!     !count opacity only if the cell is filled, else go to next cell
! 			if (lcellule_non_vide) then
! 				lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet same as iterate?
!      ! opacities in m^-1
! 
! 				l_contrib = l_contrib * AU_to_m !l_contrib in m
!      									
! 				!total bound-bound + bound-free + background opacities lte + nlte
! 				chi(:,id) = chi0_bb(:,icell)
! 				eta(:,id) = eta0_bb(:,icell)
! 				!rename it chi0_bckgr etc
! 
! 				!includes a loop over all bound-bound, passive and active
! 				call opacity_atom_loc(id,icell,iray,x0,y0,z0,x1,y1,z1,u,v,w,l,( (nbr_cell==1).and.labs ) ) 
! 				
! 				do la=1, Nlambda
! 					dtau = l_contrib * chi(la,id)
! 				
! 					if ((nbr_cell == 1).and.labs) then
! 						if (lmali_scheme) then
! 							psi(la,iray,id) = ( 1.0_dp - exp( -l*AU_to_m*chi(la,id) ) ) / chi(la,id)
! 						else
! 							ds(iray,id) = l*au_to_m
! 						endif
! 					endif
! 				
! 
! 					if (lelectron_scattering) then
! 						Snu = ( eta(la,id) + eta_es(la,icell) ) / ( chi(la,id) + tiny_dp )
! 					else
! 						Snu = eta(la,id) / ( chi(la,id) + tiny_dp )
! 					endif
! 			
! 			
! 					Itot(la,iray,id) = Itot(la,iray,id) + ( exp(-tau(la)) - exp(-(tau(la)+dtau)) ) * Snu
! 					tau(la) = tau(la) + dtau !for next cell
! 				enddo
! 
! 				do la=1, Nlambda_cont
! 					if (lelectron_scattering) then
! 						Snu_c = eta_c(la,icell) + thomson(icell) * Jnu_cont(la,icell)
! 					else
! 						Snu_c = eta_c(la,icell)
! 					endif
! 				
! 					if (Nactiveatoms > 0) then
! 						Snu_c = ( Snu_c + eta_c_nlte(la,icell) ) / ( chi_c(la,icell) + chi_c_nlte(la,icell) + tiny_dp)
! 						dtau_c = l_contrib * (chi_c(la,icell) + chi_c_nlte(la,icell))
! 					else
! 						Snu_c = Snu_c / (chi_c(la,icell) + tiny_dp)
! 						dtau_c = l_contrib * chi_c(la,icell)
! 					endif
! 					
! 					Icont(la,iray,id) = Icont(la,iray,id) + ( exp(-tau_c(la)) - exp(-(tau_c(la) + dtau_c)) ) * Snu_c
! 					tau_c(la) = tau_c(la) + dtau_c
! 				enddo
! 
! 				
! 			end if  ! lcellule_non_vide
! 		end do infinie
! 
! 	return
! 	end subroutine integ_ray_line_i_new