!
! Read model generated on a structured spherical mesh (e.g. with numpy.meshgrid).!
!
module read_spherical_grid

    use constantes
    use elements_type
    use grid, only : cell_map, vfield3d, alloc_atomrt_grid, nHtot, ne, v_char, lmagnetized, vturb, T, icompute_atomRT, &
         lcalc_ne, check_for_zero_electronic_density
    use parametres
    use messages
    use utils
    use cylindrical_grid
    use density
    use dust_prop, only : dust_pop
    use stars, only : T_hp, T_preshock
    use read1d_models, only : print_info_model

    implicit none

    integer :: header_pos

    contains

    subroutine read_spherical_grid_parameters(filename)
    ! ----------------------------------------------------- !
    ! Read a model on a spherical grid (r, theta, phi)
    ! First call read the necessary informations to
    ! generate the grid.
    !
    ! The limits of the grid (not the cell centres !) are stored
    ! in pluto%x1, pluto%x2, pluto%x3.
    ! ----------------------------------------------------- !
        character(len=*), intent(in) :: filename
        integer :: ios, Nsize, acc, i, ipos
        real :: dphi,theta_0, theta_1

        lvelocity_file = .true.
        vfield_coord = 3 ! spherical
        grid_type = 2
        n_rad_in = 1
        n_zones = 1

        open(unit=1, file=trim(filename), status="old",access="stream",form='unformatted')
        !read size along each direction + cell limits (size + 1) 
        read(1,iostat=ios) pluto%nx1
        allocate(pluto%x1(pluto%nx1+1))
        read(1,iostat=ios) pluto%x1(:)
        pluto%x1_min = minval(pluto%x1); pluto%x1_max = maxval(pluto%x1)
        ! write(*,*) "r_limits [Rstar(1)] (read)=", pluto%x1(:) / etoile(1)%r
        write(*,*) "r_limits [Rstar(1)] (read)=", pluto%x1(1) / etoile(1)%r, pluto%x1(pluto%nx1+1) / etoile(1)%r
        read(1,iostat=ios) pluto%nx2
        allocate(pluto%x2(pluto%nx2+1))
        read(1,iostat=ios) pluto%x2(:)
        pluto%x2_min = minval(pluto%x2); pluto%x2_max = maxval(pluto%x2)
        ! write(*,*) "sin(theta)_limits (read)=", pluto%x2(:)
        write(*,*) "sin(theta)_limits (read)=", pluto%x2(1), pluto%x2(pluto%nx2+1)
        !correct prec error ? 
        ! if (pluto%x2(1)-0.5*pi < 1e-6) pluto%x2(1) = 0.5*pi
        theta_0 = 180.0 * real(asin(real(pluto%x2(1),kind=dp))) / pi
        theta_1 = 180.0 * real(asin(real(pluto%x2(pluto%nx2+1),kind=dp))) / pi
        ! write(*,*) "theta_limits [°] (read)=",180.0 * real(asin(real(pluto%x2(:),kind=dp))) / pi
        write(*,*) "theta_limits [°] (read)=", theta_0, theta_1
        read(1,iostat=ios) pluto%nx3
        !special, only 1 azimuth if not 3D (' no limits ')
        Nsize = pluto%nx3
        if (pluto%nx3 > 1) then
            Nsize = Nsize + 1
        endif
        allocate(pluto%x3(Nsize))
        read(1,iostat=ios) pluto%x3(:)
        pluto%x3_min = minval(pluto%x3); pluto%x3_max = maxval(pluto%x3)
        ! write(*,*) "phi_limits [rad] (read)=", pluto%x3(:)
        write(*,*) "phi_limits [°] (read)=", 180.0 * pluto%x3(1) / pi, 180.0 * pluto%x3(size(pluto%x3)) / pi

        read(1, iostat=ios) acc
        read(1, iostat=ios) T_hp
        read(1, iostat=ios) T_preshock
        inquire(1,pos=header_pos) !position at which to start reading data.
        close(unit=1)

        if (pluto%x2(1) < pluto%x2(pluto%nx2)) then
            call error("(spherical input grid) theta(1) must be the largest value (pi or pi/2)")
        endif

        !either 2.5d (min(sin(theta))=-1) or 3d (Nphi > 1)
        l3d = (minval(pluto%x2) < 0.0).or.(pluto%nx3 > 1)
        write(*,*) "3d mode ? ", l3d

        laccretion_shock = (acc == 1)
        if (T_hp == 0.0_dp) T_hp = -1.0_dp

        ! Updating mcfost parameters
        n_rad = pluto%nx1
        n_az = pluto%nx3
        nz = pluto%nx2
        theta_max = pluto%x2_max

        disk_zone(1)%rin  = pluto%x1_min
        disk_zone(1)%edge = 0.0
        disk_zone(1)%rmin = disk_zone(1)%rin

        disk_zone(1)%rout = pluto%x1_max
        ! if (maxval(pluto%x1) == pluto%x1_max) then
        !     pluto%x1_max = pluto%x1_max + pluto%x1(pluto%nx1) - pluto%x1(pluto%nx1-1)
        ! endif
        disk_zone(1)%rmax = disk_zone(1)%rout

        !change nz to two hemispheres and check that the grid of phi is linear.
        if (l3d) then
            nz = pluto%nx2/2
            if (mod(pluto%nx2,2)>0) call warning("odd nx2")
            dphi = pluto%x3(2) - pluto%x3(1)
            do i=2, size(pluto%x3)
                if ( abs((pluto%x3(i) - pluto%x3(i-1)) - dphi) > 1e-6 ) then
                    write(*,*) dphi, (pluto%x3(i) - pluto%x3(i-1)),abs((pluto%x3(i) - pluto%x3(i-1)) - dphi)
                    call error("(spherical input grid) Non-linear phi space is not allowed!")
                endif
                dphi = pluto%x3(i) - pluto%x3(i-1)
            enddo
        endif

        Nsize = pluto%nx2 * pluto%nx3 * pluto%nx1
        write(*,*) "Nsize=", Nsize, " nx1=", pluto%nx1, " nx2=", pluto%nx2," nx3=", pluto%nx3
        write(*,*) "n_rad=", n_rad, "nz=", nz, "n_az=", n_az
        write(*,*) "rin=", real(disk_zone(1)%rin), "[au]; ", "rout=", real(disk_zone(1)%rout), "[au]"

        return
    endsubroutine read_spherical_grid_parameters

    subroutine read_spherical_model(filename)
    ! ----------------------------------------------------- !
    ! read spherical grid data defined at cell centres.
    ! ----------------------------------------------------- !
    ! use grains, only : M_Grain
        character(len=*), intent(in) :: filename
        integer :: ios, i, Nsize
        
        integer :: j, k, jj, icell
        integer, allocatable :: dz(:,:,:)
        real, allocatable :: vtmp(:,:,:,:)
        real(kind=dp), allocatable :: rho(:,:,:), rho_dust(:,:,:)
        real(kind=dp), allocatable :: T_tmp(:,:,:), ne_tmp(:,:,:), vt_tmp(:,:,:)
        real(kind=dp) :: mass

        call alloc_atomrt_grid
        call read_abundance !can be move in atom_transfer, but then rho must be changed in nHtot

        open(unit=1, file=trim(filename), status="old",access="stream",form='unformatted')
        !skip header and read data
        read (1,iostat=ios,pos=header_pos)

        Nsize =  pluto%nx1*pluto%nx2*pluto%nx3
        write(*,*) "n_cells=", n_cells, Nsize

        ! --> explicit cell mapping of the 3d arrays
        allocate(rho(pluto%nx1,pluto%nx2,pluto%nx3), ne_tmp(pluto%nx1,pluto%nx2,pluto%nx3), &
                    T_tmp(pluto%nx1,pluto%nx2,pluto%nx3), vt_tmp(pluto%nx1,pluto%nx2,pluto%nx3), &
                    dz(pluto%nx1,pluto%nx2,pluto%nx3),vtmp(pluto%nx1,pluto%nx2,pluto%nx3,3))
        read(1, iostat=ios) T_tmp(:,:,:)
        read(1, iostat=ios) rho(:,:,:)
        read(1, iostat=ios) ne_tmp(:,:,:)
        read(1, iostat=ios) vtmp(:,:,:,:)
        read(1, iostat=ios) vt_tmp(:,:,:)
        read(1, iostat=ios) dz(:,:,:)
        read(1, iostat=ios) disk_zone(1)%gas_to_dust
        write(*,*) "Gas/Dust from model:", real(disk_zone(1)%gas_to_dust)
        !read total dust density
        allocate(rho_dust(pluto%nx1,pluto%nx2,pluto%nx3))
        read(1, iostat=ios) rho_dust(:,:,:)
        close(unit=1)

        densite_pouss = 0.0_dp ! init just in case.
        masse_gaz = 0.0_dp
        disk_zone(1)%diskmass = 0.0
        do i=1, n_rad
            j = 0
            bz : do jj=j_start+1,nz-1 ! 1 extra empty cell in theta on each side
            if (jj==0) cycle bz
                j = j + 1
                do k=1, n_az
                    icell = cell_map(i,jj,k)
                    T(icell) = T_tmp(i,j,k)
                    nHtot(icell) = rho(i,j,k) * 1d3 / masseH / wght_per_H ! [H/m^3]
                    icompute_atomRT(icell) = dz(i,j,k)
                    vfield3d(icell,:) = vtmp(i,j,k,:)

                    !-> wrapper for dust RT.
                    !-> taking into account proper weights assuming only molecular gas in dusty regions
                    if (rho_dust(i,j,k) > 0.0) then ! dusty region
                        densite_gaz(icell) = rho(i,j,k) * 1d3 / masse_mol_gaz !total molecular gas density in H2/m^3
                        densite_pouss(:,icell) = rho_dust(i,j,k) * 1d3 / masse_mol_gaz ! [m^-3]
                        disk_zone(1)%diskmass = disk_zone(1)%diskmass + rho_dust(i,j,k) * volume(icell)
                    else !No dust.
                        densite_gaz(icell) = nHtot(icell) * wght_per_H !total atomic gas density in [m^-3]
                    endif
                    masse_gaz(icell) = rho(i,j,k) * volume(icell)
                enddo ! phi
            enddo bz ! theta
        enddo ! r
        masse_gaz = masse_gaz * AU3_to_m3 * 1d3 ! [g]
        deallocate(rho,ne_tmp,T_tmp,vt_tmp,dz,vtmp)
        !total dust mass
        disk_zone(1)%diskmass = disk_zone(1)%diskmass * AU3_to_m3 * kg_to_Msun
        if (allocated(rho_dust)) deallocate(rho_dust)

        !dust part see read_pluto.f90
        ! ********************************** !
        mass = sum(masse_gaz) * g_to_Msun
        ! mass = sum(densite_gaz * volume) * AU3_to_m3 * masseH * g_to_Msun
        if (mass <= 0.0) call error('Gas mass is 0')
        ! masse_gaz(:) = masseH * densite_gaz(:) * volume(:) * AU3_to_m3 ! [g]

        !--> no normalisation of density. The dust and gas densities are provided in the model.

        write(*,*) 'Total  gas mass in model:', real(sum(masse_gaz) * g_to_Msun),' Msun'
        ! write(*,*) 'Total  dust mass in model:', real(disk_zone(1)%diskmass),' Msun'
        if (disk_zone(1)%diskmass > 0.0_dp) then
            dust_pop(:)%masse = disk_zone(1)%diskmass
            !here the densite_pouss gets normalised such that sum(densite_pouss) = 1 [units less]
            call normalize_dust_density()
            !the units of densite_pouss comes from M_grain in g/cm^3, normalize to give the dust mass if summed.
            !the density of grains is M_grain x densite_pouss
            if (lemission_atom) ldust_atom = .true.
        endif
        ! ********************************** !

        call check_for_zero_electronic_density()
        call print_info_model()

        return
    endsubroutine read_spherical_model

end module read_spherical_grid