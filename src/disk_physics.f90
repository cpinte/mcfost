module disk_physics

  use grains
  use mcfost_env
  use dust_prop
  use density, only : dust_density_o_n_grains, dust_mass, find_non_empty_cell, &
       gas_density, gas_mass, gas_density_midplane, redistribute_dust_hydrostatic
  use constants
  use stars, only : star_spectrum, prob_E_star
  use messages
  use wavelengths
  use temperature
  use cylindrical_grid
  use parameters

  implicit none

  character(len=32) :: sublimationFile = "sublimation_radius.txt"
  real(kind=dp), parameter :: sub_factor = 1.6_dp

contains


subroutine compute_othin_sublimation_radius()
  ! In the optically thin case, depends only on the temperature (and spectrum) of the star

  real(kind=dp) :: E_dust, E_etoile, coeff_exp, cst_wl, sublimation_radius
  real(kind=dp) :: cst, wl, delta_wl, star_flux
  integer :: lambda, icell, i

  E_dust = 0.0
  cst=thermal_const/dust_pop(1)%T_sub

  icell = icell1

  do lambda=1, n_lambda
     ! wavelength in metres
     wl = tab_lambda(lambda) * 1.0e-6_dp
     delta_wl = tab_delta_lambda(lambda) * 1.0e-6_dp
     cst_wl=cst/wl
     if (cst_wl < 500.0_dp) then
        coeff_exp=exp(cst_wl)
        E_dust = E_dust + 4.0 * kappa_abs_LTE(icell,lambda)/((wl**5)*(coeff_exp-1.0)) *delta_wl
     endif
  enddo
  E_dust = E_dust * 2.0*pi*hp*c_light**2

  if (E_dust < tiny_real) then
     call error("Sublimation radius : opacity is not defined yet", &
          msg2="Maybe the parameter file is old ?")
  endif

  ! Emission etoiles: per-star flux via prob_E_star * star_spectrum
  do i=1, n_stars
     E_etoile = 0.0
     do lambda=1, n_lambda
        star_flux = prob_E_star(lambda,i) * star_spectrum(lambda)
        E_etoile = E_etoile + kappa_abs_LTE(icell,lambda) * star_flux / ( 4*pi * AU_to_m**2)
     enddo
     star(i)%othin_sublimation_radius = sqrt(E_etoile/E_dust)
     write(*,'(a,i0,a,f6.3,a)') "Optically thin sublimation radius (star #", i, ") = ", &
          real(star(i)%othin_sublimation_radius), " AU"
  enddo

  open(unit=1,file=trim(data_dir)//"/"//trim(sublimationFile),status="replace")
  write(1,*) star(:)%othin_sublimation_radius
  close(1)

  if (lVoronoi) then
     call sublimate_dust_Voronoi_radius(sub_factor)
  else
     sublimation_radius = real(star(1)%othin_sublimation_radius) * sub_factor
     call set_sublimation_radius(sublimation_radius)
  endif

  return

end subroutine compute_othin_sublimation_radius

!***********************************************************

subroutine set_sublimation_radius(sublimation_radius)

  real(kind=dp), intent(in) :: sublimation_radius

  integer :: i

  do i=1,n_zones
     !write(*,*) "zone", i,sublimation_radius, disk_zone(i)%rmin, sublimation_radius > disk_zone(i)%rmin
     if (sublimation_radius < disk_zone(i)%rmin) then
        !write(*,*) "SUBLIMATING DUST IN ZONE", i
        disk_zone(i)%rmin = sublimation_radius
        disk_zone(i)%rin = disk_zone(i)%rmin + 5* disk_zone(1)%edge
        disk_zone(i)%edge = 0.0
     endif
  enddo !i
  rmin = minval(disk_zone(:)%rmin)
  write(*,*) "New minimum radius = ", rmin, "AU"

  do i = 1, n_regions
     regions(i)%Rmin = max(regions(i)%Rmin,sublimation_radius)
  enddo !i

end subroutine set_sublimation_radius

!***********************************************************

subroutine read_sublimation_radius()

  real(kind=dp) :: sublimation_radius
  integer :: ios

  write(*,*) "Reading sublimation file : ./data_th/"//trim(sublimationfile)

  open(unit=1,file=trim(root_dir)//"/"//trim(seed_dir)//"/data_th/"//trim(sublimationfile), status="old", iostat=ios)
  if (ios /= 0) then
     call error("read_sublimation_radius: cannot open "//trim(root_dir)//"/"//trim(seed_dir)//"/data_th/"//trim(sublimationfile))
  endif
  read(1,*,iostat=ios) sublimation_radius
  close(unit=1)
  if (ios /= 0) then
     call error("read_sublimation_radius: error reading sublimation radius (file may be from an older run or have wrong format)")
  endif

  if (lVoronoi) then
     call sublimate_dust_Voronoi_radius(sub_factor)
  else
     sublimation_radius = real(star(1)%othin_sublimation_radius) * sub_factor
     call set_sublimation_radius(sublimation_radius)
  endif

  return

end subroutine read_sublimation_radius

!**********************************************************************

subroutine sublimate_dust_Voronoi_radius(rsub_factor)
  ! Remove dust inside the sublimation radius of each star on a Voronoi grid.
  ! Radius for star i is rsub_factor * star(i)%othin_sublimation_radius.

  use Voronoi_grid, only : Voronoi

  real(kind=dp), intent(in) :: rsub_factor

  integer :: icell, istar, n_delete_star, n_delete_total
  real(kind=dp) :: dx, dy, dz, d2, r2, sublimation_radius
  logical, dimension(:), allocatable :: lcleared

  allocate(lcleared(n_cells))
  lcleared = .false.
  n_delete_total = 0

  do istar=1, n_stars
     sublimation_radius = real(star(istar)%othin_sublimation_radius, kind=dp) * rsub_factor
     r2 = sublimation_radius**2
     n_delete_star = 0

     write(*,'(a,i0,a,f6.3,a,f6.3,a,f3.1,a)') "Star #", istar, ": sublimation radius = ", &
          real(sublimation_radius), " AU (= ", real(star(istar)%othin_sublimation_radius), &
          " x ", real(rsub_factor), ")"

     do icell=1, n_cells
        ! Skip star sink cells (no dust)
        if (Voronoi(icell)%is_star) cycle

        dx = Voronoi(icell)%xyz(1) - star(istar)%x
        if (abs(dx) > sublimation_radius) cycle
        dy = Voronoi(icell)%xyz(2) - star(istar)%y
        if (abs(dy) > sublimation_radius) cycle
        dz = Voronoi(icell)%xyz(3) - star(istar)%z
        if (abs(dz) > sublimation_radius) cycle

        d2 = dx**2 + dy**2 + dz**2
        if (d2 < r2) then
           dust_density_o_n_grains(:,icell) = 0.0
           dust_mass(icell) = 0.0
           n_delete_star = n_delete_star + 1
           if (.not.lcleared(icell)) then
              lcleared(icell) = .true.
              n_delete_total = n_delete_total + 1
           endif
        endif
     enddo

     write(*,'(a,i0,a)') "  -> removing dust in ", n_delete_star, " cells"
  enddo

  write(*,'(a,i0)') "Total cells cleared inside sublimation radii: ", n_delete_total
  deallocate(lcleared)
  call find_non_empty_cell()

  return

end subroutine sublimate_dust_Voronoi_radius

!**********************************************************************

subroutine sublimate_dust()
  ! Remove grains whose temperature exceeds the sublimation temperature
  ! C. Pinte 07/08/12
  ! D. Price 25/06/26: updated to work with all grid types (esp. Voronoi)

  integer :: icell, k, ipop, p_k
  real(kind=dp) :: mass
  logical :: lchanged

  write(*,*) "Sublimating dust"

  do icell=1, n_cells
     lchanged = .false.
     do k=1, n_grains_tot
        ipop = grain(k)%pop

        if (.not.dust_pop(ipop)%is_PAH) then
           if (Tdust(icell) > dust_pop(ipop)%T_sub) then
              ! (n_grains_tot,n_cells) if lvariable_dust, else (n_zones,n_cells)
              p_k = merge(k, grain(k)%zone, lvariable_dust)
              dust_density_o_n_grains(p_k,icell) = 0.0
              lchanged = .true.
           endif
        endif
     enddo

     ! Keep dust_mass consistent with density (needed for kappa_factor when .not.lvariable_dust)
     if (lchanged) then
        if (lvariable_dust) then
           dust_mass(icell) = 0.0
           do k=1, n_grains_tot
              dust_mass(icell) = dust_mass(icell) + dust_density_o_n_grains(k,icell) * n_grains(k) &
                   * M_grain(k) * volume(icell) * AU3_to_cm3
           enddo
        else
           ! Zone-indexed density: if wiped, mass is zero (SPH and analytic single-zone)
           if (maxval(dust_density_o_n_grains(:,icell)) <= 0.0) then
              dust_mass(icell) = 0.0
           endif
        endif
     endif
  enddo

  call find_non_empty_cell()
  mass = real(sum(dust_mass) * g_to_Msun)
  write(*,*) 'New total dust mass in model :', mass,' Msun'

end subroutine sublimate_dust


!**********************************************************

subroutine equilibre_hydrostatique()
  ! Compute hydrostatic equilibrium for the gas, then re-derive the dust
  ! vertical distribution consistently with the new gas structure
  ! (see redistribute_dust_hydrostatic in density.f90).
  !
  ! Equation 2.4.3 of the thesis (page 38, 52 of the PDF)
  ! Valid for a geometrically thin, non-self-gravitating perfect gas disk.
  !
  ! Restrictions :
  !  - cylindrical grid only, and symmetric about the midplane (j_start=1) :
  !    guarded against below (Voronoi/spherical/3D grids are rejected)
  !  - if several disk zones overlap radially, the solve works on the
  !    combined gas density, not zone by zone
  !
  ! C. Pinte
  ! 25/09/07
  ! Fixed the vertical step bug, made it actually update the density
  ! (conserving the surface density of each column), and coupled the
  ! dust settling to the new structure : C. Pinte, 08/2026

  real(kind=dp), dimension(nz) :: rho, ln_rho, Tsmooth
  real(kind=dp) :: dz_local, dTdz, fac1, fac2, M_stars, M_mol, total_sum, cst
  real(kind=dp), save, allocatable, dimension(:,:) :: Sigma_gas_target
  integer :: i,j, k, icell, alloc_status

  real(kind=dp), parameter :: hydro_relax = 0.3_dp ! under-relaxation factor for the gas density update
  real(kind=dp), parameter :: min_decline_fraction = 0.3_dp ! floor on the vertical decline rate, as a fraction of the gravity-only rate

  if (lVoronoi.or.lspherical.or.l3D) then
     call error("Hydrostatic equilibrium is only implemented for a cylindrical grid", &
          msg2="symmetric about the midplane (not a Voronoi, spherical or 3D grid)")
  endif

  M_stars = sum(star(:)%M) * Msun_to_kg
  M_mol = mu_mH * g_to_kg

  ! cst converts G*M*mu_mH/kb (SI) into a rate per AU, matching dTdz which is
  ! computed as delta-T per AU (cell_height is in AU) : z_m/r_m**3 introduces
  ! one power of 1/AU_to_m**2 (SI/AU**2), and expressing the whole rate per AU
  ! rather than per metre reintroduces one power of AU_to_m, leaving a single
  ! power of AU_to_m overall.
  ! NB: the previous AU_to_m**2 here made fac1 (the stellar tidal/gravity
  ! term) about 1e11x too small -- effectively inert -- so the density
  ! profile was driven almost entirely by noisy dTdz, with nothing to force
  ! it back down at large z. This routine was dormant (never called) until
  ! it was wired into dust_transfer.f90's hydrostatic iteration in 08/2026,
  ! which is presumably why this had never been caught.
  cst = Ggrav * M_stars * M_mol / (kb * AU_to_m)

  if (.not.allocated(Sigma_gas_target)) then
     ! Snapshot the surface density of each column the first time we are
     ! called : every hydrostatic iteration afterwards redistributes this
     ! same mass vertically, it never changes the radial (Sigma) structure.
     allocate(Sigma_gas_target(n_rad,n_az), stat=alloc_status)
     if (alloc_status > 0) call error("Allocation error Sigma_gas_target")
     do k=1, n_az
        do i=1, n_rad
           total_sum = 0.0_dp
           do j=1, nz
              total_sum = total_sum + gas_density(cell_map(i,j,k)) * cell_height(i,j)
           enddo
           Sigma_gas_target(i,k) = total_sum
        enddo
     enddo
  endif

  do k=1, n_az
     do i=1, n_rad
        ! Smooth Tdust vertically before differentiating it : dTdz feeds an
        ! exponential, cumulative integral, so raw MC noise in Tdust (cells
        ! high in the disk atmosphere receive very few photons) gets
        ! amplified into spurious density structure -- in the worst case a
        ! noise-driven upward blip that the monotonic gas_shape floor below
        ! then freezes in place for the rest of the column, leaving an
        ! unphysically puffy, optically-thin "shelf" exposed to direct
        ! starlight. A 5-cell running mean suppresses this without
        ! affecting genuine, well-resolved vertical structure.
        do j=1, nz
           icell = cell_map(i,j,k)
           Tsmooth(j) = sum(Tdust(cell_map(i,max(1,j-2):min(nz,j+2),k))) / real(min(nz,j+2)-max(1,j-2)+1,dp)
        enddo

        ln_rho(1) = 0.0_dp
        rho(1) = 1.0_dp

        do j = 2, nz
           dz_local = cell_height(i,j)
           icell = cell_map(i,j,k)
           dTdz = (Tsmooth(j)-Tsmooth(j-1)) / dz_local
           fac1 = cst * z_grid(icell)/ (r_grid(icell)**3)
           fac2 = -1.0_dp * (dTdz + fac1) / Tsmooth(j)
           ! A bound disk always settles under its own gravity, so fac2 must
           ! stay negative even where a noisy dTdz partially offsets fac1 :
           ! floor the decline rate at a fraction of the noise-free,
           ! gravity-only (isothermal) rate rather than letting it approach
           ! zero or flip sign, which would freeze or invert the column.
           fac2 = min(fac2, -min_decline_fraction * fac1 / Tsmooth(j))
           ln_rho(j) = ln_rho(j-1) + fac2 * dz_local
           rho(j) = exp(ln_rho(j))
        enddo !j

        ! Renormalisation : conserve the surface density of the column
        total_sum = 0.0_dp
        do j=1, nz
           total_sum = total_sum + rho(j) * cell_height(i,j)
        enddo

        if (total_sum > tiny_dp) then
           ! Under-relaxation : Tdust is only known from a noisy MC estimate,
           ! and dTdz feeds into an exponential, cumulative integral, so a
           ! small temperature fluctuation (worst near the inner edge, where
           ! the stellar gravity term fac1 is largest) can otherwise produce
           ! a wildly wrong density jump on a single iteration. Blending with
           ! the previous profile damps this without biasing the converged
           ! solution : both profiles integrate to the same Sigma_gas_target,
           ! so the blend does too.
           do j = 1, nz
              icell = cell_map(i,j,k)
              gas_density(icell) = (1.0_dp-hydro_relax) * gas_density(icell) + &
                   hydro_relax * (rho(j) * Sigma_gas_target(i,k) / total_sum)
           enddo
           gas_density_midplane(i,k) = gas_density(cell_map(i,1,k))
        endif

     enddo !i
  enddo !k

  do icell=1, n_cells
     gas_mass(icell) = gas_density(icell) * volume(icell) * mu_mH * AU3_to_m3
  enddo

  call redistribute_dust_hydrostatic()

  return

end subroutine equilibre_hydrostatique

!**********************************************************

end module disk_physics
