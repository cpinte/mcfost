!Interface between MCFOST and an arbitrary set of points

module arb2mcfost

    use constantes
    use grid, only : cell_map, vfield3d, alloc_atomrt_grid, nHtot, ne, v_char, lmagnetized, vturb, T, icompute_atomRT, &
         lcalc_ne, check_for_zero_electronic_density
    use parametres
    use messages
    use utils
    use sph2mcfost, only : SPH_to_Voronoi, Hydro_to_Voronoi_atomic
    use sort, only : find_kth_smallest_inplace
    use elements_type
    use stars, only : T_hp, T_preshock

    implicit none

    contains



      ! subroutine SPH_to_Voronoi(n_SPH, ndusttypes, particle_id, x,y,z,h, vx,vy,vz, T_gas, massgas,massdust,rho,rhodust,&
      !      SPH_grainsizes, SPH_limits, check_previous_tesselation, is_ghost, mask)
      !
      !   ! ************************************************************************************ !
      !   ! n_sph : number of points in the input model
      !   ! ndusttypes : number of dust species in the input model
      !   ! particle_id : index of the particle / cell centre in the input model
      !   ! x,y,z : coordinates of the particle / cell
      !   ! h : smoothing length to cut elongated cells
      !   ! vx,vy,vz : vector components of the velocity fields (cartesian) in the input model
      !   ! massgas : total mass of the gas
      !   ! massdust : total mass of the dust
      !   ! rho : gas density
      !   ! rhodust : dust density
      !   ! sph_grainsizes : size of grains for sph models
      !   ! sph_limits : limit of the input model box
      !   ! check_previous_tesselation :
      !   ! mask :
      !   ! ************************************************************************************ !
      !   use Voronoi_grid
      !   use density, only : densite_gaz, masse_gaz, densite_pouss, masse
      !   use grains, only : n_grains_tot, M_grain
      !   use disk_physics, only : compute_othin_sublimation_radius
      !   use mem
      !
      !   integer, intent(in) :: n_SPH, ndusttypes
      !   real(dp), dimension(n_SPH), intent(inout) :: x,y,z,h,massgas!,rho, !move rho to allocatable, assuming not always allocated
      !   real(dp), dimension(:), allocatable, intent(inout) :: rho
      !   real(dp), dimension(:), allocatable, intent(inout) :: vx,vy,vz ! dimension n_SPH or 0
      !   real(dp), dimension(:), allocatable, intent(in) :: T_gas
      !   integer, dimension(n_SPH), intent(in) :: particle_id
      !   real(dp), dimension(:,:), allocatable, intent(inout) :: rhodust, massdust ! ndusttypes,n_SPH
      !   real(dp), dimension(:), allocatable, intent(in) :: SPH_grainsizes ! ndusttypes
      !   real(dp), dimension(6), intent(in) :: SPH_limits
      !   logical, intent(in) :: check_previous_tesselation
      !   logical, dimension(:), allocatable, intent(in), optional :: mask
      !
      !   integer, dimension(:), allocatable, intent(out) :: is_ghost





    subroutine setup_arb_to_mcfost(x, y, z, h, vx, vy, vz, mass_gas, particle_id)
        integer                                  :: n_points,k,icell
        real(kind=dp), allocatable, dimension(:), intent(inout) :: x,y,z,h,vx,vy,vz,mass_gas
        integer,  allocatable, dimension(:), intent(inout)      :: particle_id
        real(kind=dp), allocatable, dimension(:) :: mass_ne_on_massgas, T_tmp, vt_tmp, dz
        real(kind=dp), allocatable, dimension(:) :: rho, hydro_grainsizes
        real(kind=dp), allocatable, dimension(:,:) :: rhodust, massdust, dust_moments
        integer,       allocatable, dimension(:) :: is_ghost
        ! logical, allocatable, dimension(:) :: mask ! size == np, not n_SPH, index is original SPH id
        real(dp) :: mass_per_H

        integer, parameter                       :: Nheader = 2 !Add more, for ascii file
        integer                                  :: syst_status, acspot, alloc_status
        character(len=512)                       :: inputline, FormatLine, cmd
        logical                                  :: check_previous_tesselation, ldust_moments
        ! integer,  allocatable, dimension(:)      :: particle_id
        real(kind=dp), dimension(6)              :: hydro_limits
        integer                                  :: ndusttypes !, voroindex, N_fixed_ne = 0
        real, parameter                          :: limit_factor = 1.005!, Lextent = 1.01
        integer, dimension(:), allocatable :: mask


        !There will be an error if lphantom_file is true. Because density_files and density_file
        !stores pluto's model name. But also the filename from phantom.. So to date, the two
        !codes cannot be merged.
        !This is to be able to use read previous tesselation


        lmagnetized = .false.
        ldust_moments = .false.
        ! lfix_star = .true.
        lphantom_file = .false.
        larg_voronoi = .true.
        lignore_dust = .true.
        lrandomize_voronoi = .false.
        check_previous_tesselation = (.not.lrandomize_voronoi)
        n_points = size(x)

        ! !needed for Voronoi
        ! if (allocated(density_files)) deallocate(density_files)
        ! allocate(density_files(1)); density_files(1) = density_file

        if (lignore_dust) then
           ndusttypes = 0
           if (allocated(rhodust)) deallocate(rhodust,massdust)
        else
           call error("Dust not handled yet for arb models!")
        endif

           ! N_points = N_points - Nheader


        hydro_limits(:) = 0

        k = 1
        hydro_limits(1) = find_kth_smallest_inplace(k,real(x))*limit_factor
        hydro_limits(3) = find_kth_smallest_inplace(k,real(y))*limit_factor
        hydro_limits(5) = find_kth_smallest_inplace(k,real(z))*limit_factor

        k = n_points
        hydro_limits(2) = find_kth_smallest_inplace(k,real(x))*limit_factor
        hydro_limits(4) = find_kth_smallest_inplace(k,real(y))*limit_factor
        hydro_limits(6) = find_kth_smallest_inplace(k,real(z))*limit_factor

        !also work with grid-based code
        !massdust, rhodust, hydro_grainsizes not allocated if ndusttypes = 0 !
        call sph_to_voronoi(n_points, ndusttypes, particle_id, x, y, z, h, vx, vy, vz, &
             T_tmp, mass_gas, massdust, rho, rhodust, hydro_grainsizes, hydro_limits, check_previous_tesselation, is_ghost, &
             ldust_moments, dust_moments, mass_per_H, mask)
       ! subroutine SPH_to_Voronoi(n_SPH, ndusttypes, particle_id, x,y,z,h, vx,vy,vz, T_gas, massgas,massdust,rho,rhodust,&
       !      SPH_grainsizes, SPH_limits, check_previous_tesselation, is_ghost, ldust_moments, dust_moments, mass_per_H, mask)

        ! -> correction for small density applied on mass_gas directly inside


        ! call hydro_to_Voronoi_atomic(n_points,T_tmp,vt_tmp,mass_gas,mass_ne_on_massgas,dz)
        ! 	call empty_cells

        !deallocating temporary variables from input file.
        ! deallocate(h,vx,vy,vz,mass_gas, mass_ne_on_massgas, x,y,z,T_tmp, vt_tmp, dz)
        deallocate(h,vx,vy,vz,mass_gas, x,y,z)

        return
      end subroutine setup_arb_to_mcfost

end module arb2mcfost
