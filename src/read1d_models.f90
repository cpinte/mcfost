module read1d_models
!
! Interface with stellar atmosphere models such as MARCS, Kurucz, CMFGEN and MULTI.
! The input file is in a format which is common to all of these codes (that have very different formatting).
! Still, the format is easily obtained from a small editing of the original model files (from MARCS for instance).
!
! TO DO:
!	direct interface with these different codes (and format)
!
	use parametres
	use messages
	use mcfost_env
	use constantes
    use elements_type
    use grid, only : cell_map, vfield3d, alloc_atomrt_grid, nHtot, ne, v_char, lmagnetized, vturb, T, icompute_atomRT, &
         lcalc_ne, check_for_zero_electronic_density

	implicit none

	!.false. mainly for debug
	logical, parameter :: lcell_centered = .true.

	contains
	
	
	subroutine read_model_1d(mod_file)
		character(len=*), intent(in) :: mod_file
		integer :: i, nrad_0
		real(kind=dp) :: Fnorm, dnu, Rmax_c
		
		grid_type = 2
		n_rad_in = 1
		n_az = 1
		nz = 1
		theta_max = 0.5 * pi !unused

		lvelocity_file = .true.
		vfield_coord = 3
		lmagnetized = .false.
		lcalc_ne = .false.
		
		open(unit=1,file=mod_file, status="old")
		read(1,*) atmos_1d%rstar
		read(1,*) atmos_1d%nr
		allocate(atmos_1d%r(atmos_1d%nr),atmos_1d%T(atmos_1d%nr),atmos_1d%ne(atmos_1d%nr),atmos_1d%rho(atmos_1d%nr))
		allocate(atmos_1d%iz(atmos_1d%nr),atmos_1d%vt(atmos_1d%nr),atmos_1d%v(atmos_1d%nr,3))
		! write(*,*) "Rstar = ", Rstar_mod
		! write(*,*) "Nradii = ", Nr
		do i=1, atmos_1d%nr
			read(1,*) atmos_1d%r(i), atmos_1d%T(i), atmos_1d%rho(i), atmos_1d%ne(i), &
				atmos_1d%vt(i), atmos_1d%v(i,1), atmos_1d%v(i,2), atmos_1d%v(i,3), atmos_1d%iz(i)
			! write(*,*) atmos_1d%r(i), atmos_1d%T(i), atmos_1d%rho(i), atmos_1d%ne(i), &
				! atmos_1d%vt(i), atmos_1d%v(i,1), atmos_1d%v(i,2), atmos_1d%v(i,3), atmos_1d%iz(i)
		enddo
		!if corona will be zero at some point
		! write(*,*) "T_limits (read):", atmos_1d%T(1), atmos_1d%T(atmos_1d%nr)
		! write(*,*) "rho_limits (read):", atmos_1d%rho(1), atmos_1d%rho(atmos_1d%nr)
		! write(*,*) "ne_limits (read):", atmos_1d%ne(1), atmos_1d%ne(atmos_1d%nr)
		! write(*,*) "vturb_limits (read):", atmos_1d%vt(1), atmos_1d%vt(atmos_1d%nr)
		! write(*,*) "vr_limits (read):", atmos_1d%v(1,1), atmos_1d%v(atmos_1d%nr,1)
		! write(*,*) "vtheta_limits (read):", atmos_1d%v(1,2), atmos_1d%v(atmos_1d%nr,2)
		! write(*,*) "vphi_limits (read):", atmos_1d%v(1,3), atmos_1d%v(atmos_1d%nr,3)


		call check_for_coronal_illumination()
		if (atmos_1d%lcoronal_illumination) then
			write(*,*) " -> Reading (isotropic) coronal illumination from 1d model!"
			read(1,*) atmos_1d%Ncorona, atmos_1d%E_corona
			allocate(atmos_1d%x_coro(atmos_1d%Ncorona), atmos_1d%I_coro(atmos_1d%Ncorona,1))!1 ray at the moment
			Fnorm = 0
			do i=1, atmos_1d%Ncorona
				read(1,*) atmos_1d%x_coro(i), atmos_1d%I_coro(i,1)
				if (i>1) then
					dnu = abs(1d9 * c_light / atmos_1d%x_coro(i) - 1d9 * c_light / atmos_1d%x_coro(i-1))
					Fnorm = Fnorm + 0.5 * pi * (atmos_1d%I_coro(i-1,1) + atmos_1d%I_coro(i,1)) * dnu
				endif
			enddo
			if (atmos_1d%E_corona > 0) then!log10(F_SI) + 3 = log10(F_cgs)
				atmos_1d%I_coro(:,1) = atmos_1d%I_coro(:,1) * atmos_1d%E_corona / Fnorm
				write(*,'("lgFcorona="(1F12.4)" W/m^2; lgFcorona(old)="(1F12.4)" W/m^2")') &
					log10(atmos_1d%E_corona), log10(Fnorm)
				write(*,'("max(I)="(1ES14.5E3)" W/m^2/Hz/sr; min(I)="(ES14.5E3)" W/m^2/Hz/sr")') &
					maxval(atmos_1d%I_coro), minval(atmos_1d%I_coro)
			else
				write(*,'("lgFcorona="(1F12.4)" W/m^2")') log10(Fnorm)
			endif
		endif
		close(unit=1)

		n_rad = atmos_1d%nr - 1
		n_cells = n_rad
		n_etoiles = 1 !force		

		Rmin = minval(atmos_1d%r) * atmos_1d%rstar * m_to_au
		Rmax = maxval(atmos_1d%r,mask=(atmos_1d%iz>0)) * atmos_1d%rstar* m_to_au
		Rmax_c = maxval(atmos_1d%r) * atmos_1d%rstar * m_to_au !corona extent
		!it is not the star but the inner boundary of the model !
		!Because the star is the model (in general !) !
		etoile(1)%r = Rmin
		etoile(1)%x = 0.0_dp
		etoile(1)%y = 0.0_dp
		etoile(1)%z = 0.0_dp

		!convert to au
		atmos_1d%r = atmos_1d%r * atmos_1d%rstar * m_to_au

		allocate(atmos_1d%m(atmos_1d%nr))
		nrad_0 = atmos_1d%nr
		if (atmos_1d%lcoronal_illumination) nrad_0 = atmos_1d%nr - 1
		atmos_1d%m(nrad_0) = 0.5 * (atmos_1d%rho(nrad_0-1)+atmos_1d%rho(nrad_0)) * &
									au_to_m * (atmos_1d%r(nrad_0) - atmos_1d%r(nrad_0-1))
		do i=nrad_0-1,1,-1
			atmos_1d%m(i) = atmos_1d%m(i+1) + 0.5 * (atmos_1d%rho(i)+atmos_1d%rho(i+1)) * &
							au_to_m * (atmos_1d%r(i+1) - atmos_1d%r(i))
		enddo

		disk_zone(1)%rin  = Rmin
		disk_zone(1)%edge = 0.0
		disk_zone(1)%rmin = disk_zone(1)%rin
		disk_zone(1)%rout = Rmax_c
		disk_zone(1)%rmax = disk_zone(1)%rout
		atmos_1d%s = (atmos_1d%rstar / Rmax)**2

		write(*,*) "WARNING distance and map size set to : "
		distance = Rmax * au_to_pc !pc
		map_size = 2.005 * Rmax !au
		write(*,*) distance, ' pc', map_size * 1e3, 'mau'	

		return
	end subroutine read_model_1d

	subroutine setup_model1d_to_mcfost()
		real(kind=dp) :: rho_to_nH
		integer :: i, icell, nrad_0, icell0
		real(kind=dp), dimension(:), allocatable :: cm

		call alloc_atomrt_grid()
		call read_abundance
		rho_to_nH = 1d3 / masseH / wght_per_H

		!- MCFOST is a cell centered radiative transfer code while the 1D spherically symmetric codes
		! are based on a set of points corresponding to the nodes of multi-dimensional models.
		!Therefore, there is one cell less than the number of points in the model.
		!If lcell_centered, the points in the file are interpreted at cell edges, and the value
		!of thermodynamics quantities for each cell is computed from 1/2 * (q(k) + q(k-1)).
		!This is more accurate, preserve the column mass, but makes it harder to benchmark with other 1D
		!codes, as the point values does not correspond to the node values, hence the emission/absorption coefficients
		!are different (while they should be equivalent for the same values of rho, T, ne, except for bugs).
		!- If not lcell_centered, the outermost cell (limit = Rmax) will be given the outermost point values (rho(nr), etc.).
		!The inner most point will be skipped (limit at Rmin). This option makes it easier to benchmark opacities and emissivities
		!locally with other codes. USE that option for DEBUG (if the emerging flux/intensity are too much different from the reference).
		!NOTE: because the last point is remove, its contribution to continuum and line is removed. While lines generally do not form
		!at the bottom of the photosphere, it can slighly modify the continuum emission. Playing with the inner boundary irradiation
		!allows to properly balance that effect.
		if (lcell_centered) then
			etoile(1)%T = real(	atmos_1d%T(1) ) ! + something else here 
			icell = n_cells + 1
			nrad_0 = n_rad
			icell0 = n_cells
			if (atmos_1d%lcoronal_illumination) then
				nrad_0 = n_rad - 1 !avoid the coronal point
				icell = icell - 1
				icell0 = n_cells - 1
				icompute_atomRT(n_cells) = atmos_1d%iz(n_rad+1)
			endif
			do i=nrad_0, 1, -1
				icell = icell - 1
				icompute_atomRT(icell) = atmos_1d%iz(i)
				T(icell) = 0.5 * (atmos_1d%T(i) + atmos_1d%T(i+1))
				nHtot(icell) = 0.5 * (atmos_1d%rho(i) + atmos_1d%rho(i+1)) * rho_to_nH
				ne(icell) = 0.5 * (atmos_1d%ne(i) + atmos_1d%ne(i+1))
				vfield3d(icell,1) = 0.5 * ( atmos_1d%v(i,1) + atmos_1d%v(i+1,1) )! vr
				vfield3d(icell,2) = 0.5 * ( atmos_1d%v(i,3) + atmos_1d%v(i+1,3) ) ! vphi
				vfield3d(icell,3) = 0.5 * ( atmos_1d%v(i,2) + atmos_1d%v(i+1,2) ) ! vtheta
				vturb(icell) = 0.5 * (atmos_1d%vt(i) + atmos_1d%vt(i+1))
				! write(*,*) i, icell, T(icell), nHtot(icell), ne(icell), icompute_atomRT(icell)
			enddo
		else !for debug better, direct association points-to-cells
			nrad_0 = n_rad
			icell0 = n_cells
			if (atmos_1d%lcoronal_illumination) then
				nrad_0 = n_rad - 1 !avoid the coronal point
				icell0 = n_cells - 1
				icompute_atomRT(n_cells) = atmos_1d%iz(n_rad+1)
			endif
			etoile(1)%T = real(	atmos_1d%T(1) ) ! + irrad ? 
			icompute_atomRT(:) = atmos_1d%iz(2:n_cells+1)
			T(:) = atmos_1d%T(2:n_cells+1)
			nHtot(:) = atmos_1d%rho(2:n_cells+1) * rho_to_nH
			ne(:) = atmos_1d%ne(2:n_cells+1)
			vfield3d(:,1) = atmos_1d%v(2:n_cells+1,1) ! vr
			vfield3d(:,2) = atmos_1d%v(2:n_cells+1,3) ! vphi
			vfield3d(:,3) = atmos_1d%v(2:n_cells+1,2) ! vtheta
			vturb(:) = atmos_1d%vt(2:n_cells+1)
		endif
! TO DO: renorm to the column mass of the input points

		write(*,*) "Tbot (1d model)", etoile(1)%T, ' K'

		!column mass
		allocate(cm(n_cells))
		cm(icell0) = au_to_m * (atmos_1d%r(icell0+1) - atmos_1d%r(icell0)) * nHtot(icell0) / rho_to_nH
		write(*,*) "icell=", icell0, " T=", T(icell0), " nH=", nHtot(icell0), " ne=", ne(icell0), &
				" iz=", icompute_atomRT(icell0)
		do icell=icell0-1,1,-1
			write(*,*) "icell=", icell, " T=", T(icell), " nH=", nHtot(icell), " ne=", ne(icell), &
				" iz=", icompute_atomRT(icell)
			cm(icell) = cm(icell+1) +  au_to_m * (atmos_1d%r(icell+1) - atmos_1d%r(icell)) * nHtot(icell) / rho_to_nH
			write(*,*) " cm = ", cm(icell), " model:", atmos_1d%m(icell)
		enddo
		write(*,*) "cm in MCFOST:", maxval(cm), " cm in model: ", maxval(atmos_1d%m)


		call check_for_zero_electronic_density()
		call print_info_model()

		!atmos_1d%r !deallocated elsewhere
		! deallocate(atmos_1d%T, atmos_1d%ne, atmos_1d%rho, atmos_1d%v, atmos_1d%vt, atmos_1d%iz)
		deallocate(cm)

		return
	end subroutine setup_model1d_to_mcfost

	subroutine print_info_model
		real(kind=dp) :: v_char

		v_char = sqrt( maxval(sum(vfield3d**2,dim=2)) )

		write(*,*) "Maximum/minimum velocities in the model (km/s):"
		write(*,*) " Vfield(1) = ", &
			1e-3 * maxval(abs(vfield3d(:,1))), 1d-3*minval(abs(vfield3d(:,1)),mask=icompute_atomRT>0)
		write(*,*) " Vfield(2) = ",  &
			1d-3 * maxval(abs(vfield3d(:,2))), 1d-3*minval(abs(vfield3d(:,2)),mask=icompute_atomRT>0)
		write(*,*) " Vfield(3) = ",  &
			1d-3 * maxval(abs(vfield3d(:,3))), 1d-3*minval(abs(vfield3d(:,3)),mask=icompute_atomRT>0)


		write(*,*) "Typical line extent due to V fields (km/s):"
		write(*,*) v_char/1d3

		write(*,*) "Maximum/minimum turbulent velocity (km/s):"
		write(*,*) maxval(vturb)/1d3, minval(vturb, mask=icompute_atomRT>0)/1d3

		write(*,*) "Maximum/minimum Temperature in the model (K):"
		write(*,*) real(maxval(T)), real(minval(T,mask=icompute_atomRT>0))
		write(*,*) "Maximum/minimum Hydrogen total density in the model (m^-3):"
		write(*,*) real(maxval(nHtot)), real(minval(nHtot,mask=icompute_atomRT>0))
		if (.not.lcalc_ne) then
			write(*,*) "Maximum/minimum ne density in the model (m^-3):"
			write(*,*) real(maxval(ne)), real(minval(ne,mask=icompute_atomRT>0))
		endif

		if (maxval(icompute_atomRT)<=0) then
			call warning(" There is no gas density zone in the model!")
			return
		endif
        write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT>0)), " density zones"
        write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT==0)), " transparent zones"
        write(*,*) "Read ", size(pack(icompute_atomRT,mask=icompute_atomRT<0)), " dark zones"
		write(*,'("-- Solving RTE for "(1F6.2)" % of cells")') &
			100.0*real(size(pack(icompute_atomRT,mask=icompute_atomRT>0))) / real(n_cells)

		return
	end subroutine print_info_model

	subroutine check_for_coronal_illumination()
		integer :: i, Ncorona_points
		atmos_1d%lcoronal_illumination = .false.

		i_loop : do i=1,size(atmos_1d%iz(:))
			if (atmos_1d%iz(i) == -2) then

				Ncorona_points = size(pack(atmos_1d%iz,mask=(atmos_1d%iz==-2)))
				write(*,*) "  *** Found ", Ncorona_points," Coronal illumination zones in 1d model! *** "
				atmos_1d%lcoronal_illumination = .true.
				exit i_loop

			endif
		enddo i_loop

		return
	end subroutine check_for_coronal_illumination
	
end module read1d_models