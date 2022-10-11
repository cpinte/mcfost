module read1d_models
!!!
! Read 1d atmospheric models in ascii format.
! WARNING: the number of cells is then npoints - 1 where npoints
! is the number of radius in the atmospheric model.
! TO DO: add a ghost point to prevent that ???
!
! In the future, atmos_type disappear which will simplify the dependency of the module
! and brings clarity.
!
!!!
	use parametres
	use messages
	use mcfost_env
	use constantes

	implicit none

	real(kind=dp), dimension(:), allocatable :: tab_r_mod1d, tab_T_mod1, tab_rho_mod1, tab_ne_mod1
	real(kind=dp), allocatable :: tab_v_mod1(:,:), tab_vt_mod1(:), Icorona(:,:), xcorona(:)
	integer, dimension(:), allocatable :: tab_zone_mod1
	logical :: lcoronal_illumination !from above

	contains
	
	
	subroutine read_grid_1d(mod_file)
		!first loop on the model to read the list of radii to generate the grid
		character(len=*), intent(in) :: mod_file
		real(kind=dp) ::  Rstar_mod, Fcorona
		integer :: Nr, i, Nlambda_corona
		real(kind=dp) :: Fnorm, dnu, Rmax_c
		
		grid_type = 2
		n_rad_in = 1
		n_az = 1
		nz = 1

		lvelocity_file = .true.
		vfield_coord = 3
		
		open(unit=1,file=mod_file, status="old")
		read(1,*) Rstar_mod
		read(1,*) Nr
		allocate(tab_r_mod1d(Nr),tab_T_mod1(Nr),tab_ne_mod1(Nr),tab_rho_mod1(Nr))
		allocate(tab_zone_mod1(Nr),tab_vt_mod1(Nr),tab_v_mod1(3,Nr))
		! write(*,*) "Rstar = ", Rstar_mod
		! write(*,*) "Nradii = ", Nr
		do i=1, Nr
			read(1,*) tab_r_mod1d(i), tab_T_mod1(i), tab_rho_mod1(i), tab_ne_mod1(i), &
				tab_vt_mod1(i), tab_v_mod1(1,i), tab_v_mod1(2,i), tab_v_mod1(3,i), tab_zone_mod1(i)
			! write(*,*) i, tab_r_mod1d(i), tab_T_mod1(i), tab_rho_mod1(i), tab_ne_mod1(i), &
			! 	tab_vt_mod1(i), tab_v_mod1(1,i), tab_v_mod1(2,i), tab_v_mod1(3,i), tab_zone_mod1(i)
		enddo
		call check_for_coronal_illumination()
		if (lcoronal_illumination) then
			write(*,*) " -> Reading (isotropic) coronal illumination from 1d model!"
			read(1,*) Nlambda_corona, Fcorona
			allocate(xcorona(Nlambda_corona), Icorona(Nlambda_corona,1))!1 ray at the moment
			Fnorm = 0
			do i=1, Nlambda_corona
				read(1,*) xcorona(i), Icorona(i,1)
				! write(*,*) xcorona(i), Icorona(i,1)
				if (i>1) then
					dnu = abs(1d9 * c_light / xcorona(i) - 1d9 * c_light / xcorona(i-1))
					Fnorm = Fnorm + 0.5 * pi * (Icorona(i-1,1) + Icorona(i,1)) * dnu
				endif
			enddo
			if (Fcorona > 0) then!log10(F_SI) + 3 = log10(F_cgs)
				Icorona(:,1) = Icorona(:,1) * Fcorona / Fnorm
				write(*,'("lgFcorona="(1F12.4)" W/m^2; lgFcorona(old)="(1F12.4)" W/m^2")') log10(Fcorona), log10(Fnorm)
				write(*,'("max(I)="(1ES14.5E3)" W/m^2/Hz/sr; min(I)="(ES14.5E3)" W/m^2/Hz/sr")') maxval(Icorona), minval(Icorona)
			else
				write(*,'("lgFcorona="(1F12.4)" W/m^2")') log10(Fnorm)
			endif
		endif
		close(unit=1)

		n_rad = Nr - 1
		n_cells = n_rad
		
		n_etoiles = 1 !force
		!it is not the star but the inner boundary of the model !
		! Rmin = minval(tab_r_mod1d) * Rstar_mod * m_to_au
		Rmin = minval(tab_r_mod1d) * Rstar_mod * m_to_au
		Rmax = maxval(tab_r_mod1d,mask=(tab_zone_mod1>0)) * Rstar_mod * m_to_au
		Rmax_c = maxval(tab_r_mod1d) * Rstar_mod * m_to_au
		etoile(1)%r = Rmin
		etoile(1)%x = 0.0_dp
		etoile(1)%y = 0.0_dp
		etoile(1)%z = 0.0_dp
		!Temperature read from file
		!etoile(1)%T = real(tab_T_mod1(1)) or 2 if there is a ghost point
		!other elements not useful in 1d mode !
		tab_r_mod1d = tab_r_mod1d * Rstar_mod * m_to_au

		disk_zone(1)%rin  = Rmin
		disk_zone(1)%edge=0.0
		disk_zone(1)%rmin = disk_zone(1)%rin
		disk_zone(1)%rout = Rmax_c!Rmax
		disk_zone(1)%rmax = disk_zone(1)%rout

		write(*,*) "WARNING distance and map size set to : "
		distance = Rmax * au_to_pc !pc
		map_size = 2.005 * Rmax !au
		write(*,*) distance, ' pc', map_size * 1e3, 'mau'		
	
		return
	end subroutine read_grid_1d

	subroutine check_for_coronal_illumination()
		integer :: i, Ncorona_points
		lcoronal_illumination = .false.

		i_loop : do i=1,size(tab_zone_mod1)
			if (tab_zone_mod1(i) == -2) then

				Ncorona_points = size(pack(tab_zone_mod1,mask=(tab_zone_mod1==-2)))
            	write(*,*) "  *** Found ", Ncorona_points," Coronal illumination zones in 1d model! *** "
            	lcoronal_illumination = .true.
            	exit i_loop

			endif
   		enddo i_loop

	return
	end subroutine check_for_coronal_illumination
	
end module read1d_models