module read1d_models
!!!
! Read 1d atmospheric models in ascii format.
! WARNING: the number of cells is then npoints - 1 where npoints
! is the number of radius in the atmospheric model.
!
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
	real(kind=dp), allocatable :: tab_v_mod1(:,:), tab_vt_mod1(:)
	integer, dimension(:), allocatable :: tab_zone_mod1

	contains
	
	
	subroutine read_grid_1d()
		!first loop on the model to read the list of radii to generate the grid
		real(kind=dp) ::  Rstar_mod
		integer :: Nr, i
		
		grid_type = 2
		n_rad_in = 1
		n_az = 1
		nz = 1

		! lmagnetoaccr = .false.
		! lspherical_velocity = .true.
		! lvoronoi = .false.
		! lmagnetized = .false.
		! calc_ne = .false.	
		
		open(unit=1,file=density_file, status="old")
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
		close(unit=1)
		n_rad = Nr - 1
		n_cells = n_rad
		
		n_etoiles = 1 !force
		!it is not the star but the inner boundary of the model !
		Rmin = minval(tab_r_mod1d) * Rstar_mod * m_to_au
		Rmax = maxval(tab_r_mod1d) * Rstar_mod * m_to_au
		etoile(1)%r = Rmin
		etoile(1)%x = 0.0_dp
		etoile(1)%y = 0.0_dp
		etoile(1)%z = 0.0_dp
		!other elements not useful in 1d mode !
		tab_r_mod1d = tab_r_mod1d * Rstar_mod * m_to_au

		disk_zone(1)%rin  = Rmin
		disk_zone(1)%edge=0.0
		disk_zone(1)%rmin = disk_zone(1)%rin
		disk_zone(1)%rout = Rmax
		disk_zone(1)%rmax = Rmax

		write(*,*) "WARNING distance and map size set to : "
		distance = Rmax * au_to_pc !pc
		map_size = 2.005 * Rmax !au
		write(*,*) distance, ' pc', map_size * 1e3, 'mau'
	
		return
	end subroutine read_grid_1d
	
end module read1d_models