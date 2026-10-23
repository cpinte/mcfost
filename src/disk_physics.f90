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
  use sort, only : find_kth_smallest

  implicit none

  character(len=32) :: sublimationFile = "sublimation_radius.txt"
  ! The innermost zone is extended inwards to sub_factor x the optically thin
  ! sublimation radius if it starts further out; the sublimation front itself is
  ! found by the iterations (sublimation mask)
  real(kind=dp), parameter :: sub_factor = 1.0_dp

  ! Surface density snapshot taken on the first hydrostatic iteration,
  ! freed by free_hydrostatic once the iterations are complete
  real(kind=dp), allocatable, dimension(:,:), private :: Sigma_gas_target

  ! Hydrostatic density convergence : cells below this fraction of their
  ! column midplane density are excluded from the residual
  real(kind=dp), parameter :: hydro_residual_floor = 1.0e-3_dp
  real(kind=dp), parameter :: hydro_density_precision = 5.0e-2_dp

  ! Dust sublimation, iterated with the radiative transfer in the physics loop.
  ! dust_density_pristine is the dust density before sublimation (it follows the
  ! hydrostatic redistribution); lsublimated has the first dimension of
  ! dust_density_o_n_grains : one flag per grain size (lvariable_dust) or per zone.
  real(kind=dp), allocatable, dimension(:,:), private :: dust_density_pristine
  logical, allocatable, dimension(:,:) :: lsublimated
  integer(kind=1), allocatable, dimension(:,:), private :: n_sublimation_flips
  integer, allocatable, dimension(:), private :: pop_of_pk ! population of each first index of dust_density_o_n_grains
  character(len=32), parameter :: sublimationMaskFile = "sublimation_mask"

  ! A trace of the pristine dust is left in sublimated cells : negligible for the
  ! transport and the observables, but the cell still records the radiation field,
  ! giving the equilibrium ("test particle") temperature of its dust, which is
  ! independent of the dust density. Sublimated dust can then re-condense.
  real(kind=dp), parameter :: sublimation_trace = 1.0e-20_dp
  ! Hysteresis band around T_sub : sublimation if T > T_sub (1+delta), re-condensation
  ! if T < T_sub (1-delta), with delta from the Monte Carlo noise on T (as for the
  ! hydrostatic convergence test), between sublimation_precision and sublimation_band_max
  real, parameter :: sublimation_precision = 0.02
  real, parameter :: sublimation_n_sigma = 3.0
  real, parameter :: sublimation_band_max = 0.3
  real, parameter :: sublimation_min_packets = 10.0
  ! a grain that changed state this many times keeps its current state (no oscillation)
  integer, parameter :: sublimation_max_flips = 4
  ! The sublimation front is converged when fewer (grain or zone, cell) pairs than this
  ! fraction of the sublimated ones would change state : the remaining changes are
  ! Monte Carlo noise near T_sub (in particular with one temperature per grain size)
  real, parameter :: sublimation_change_tolerance = 1.0e-3
  ! floor on the gas density in the power-law sublimation temperature [g/cm3]
  real, parameter :: sublimation_rho_floor = 1.0e-30
  ! return values of sublimation_change
  integer, parameter :: sublimation_no_change = 0, sublimation_sublimates = 1, &
       sublimation_recondenses = -1, sublimation_frozen = 2

contains


subroutine compute_othin_sublimation_radius()
  ! In the optically thin case, depends only on the temperature (and spectrum) of the star

  real(kind=dp) :: E_dust, E_etoile, coeff_exp, cst_wl, sublimation_radius
  real(kind=dp) :: cst, wl, delta_wl, star_flux
  real :: T_sub
  integer :: lambda, icell, i

  ! Highest sublimation temperature of the first population in the model (density-
  ! dependent laws) : no dust can survive inside the corresponding radius
  T_sub = 0.0
  do icell=1, n_cells
     T_sub = max(T_sub, population_sublimation_temperature(1, real(gas_density(icell) * mu_mH / m3_to_cm3)))
  enddo
  write(*,'(a,f7.1,a)') " Sublimation temperature for the optically thin sublimation radius = ", T_sub, " K"

  E_dust = 0.0
  cst=thermal_const/T_sub

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

  ! data_dir is not created when the run stops after the initialisation (eg -disk_struct)
  if (.not.lstop_after_init) then
     open(unit=1,file=trim(data_dir)//"/"//trim(sublimationFile),status="replace")
     write(1,*) star(:)%othin_sublimation_radius
     close(1)
  endif

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
  real(kind=dp) :: rmin_zones

  ! Only the innermost zone(s) are extended inwards to the sublimation radius :
  ! extending the other zones would spread their dust inside the inner zones
  ! (and make the grid regions overlap)
  rmin_zones = minval(disk_zone(:)%rmin)
  do i=1,n_zones
     if ((disk_zone(i)%rmin <= rmin_zones) .and. (sublimation_radius < disk_zone(i)%rmin)) then
        !write(*,*) "SUBLIMATING DUST IN ZONE", i
        disk_zone(i)%rmin = sublimation_radius
        disk_zone(i)%rin = disk_zone(i)%rmin + 5* disk_zone(i)%edge
        disk_zone(i)%edge = 0.0
     endif
  enddo !i
  ! The radial grid is built from the regions (cylindrical_grid), with the first
  ! cell starting at rmin : the regions must follow the zones, as in
  ! define_physical_zones, otherwise the first radial cell spans from the new rmin
  ! to the old inner radius of the region.
  do i = 1, n_regions
     regions(i)%Rmin = minval(disk_zone(:)%rmin, mask=(disk_zone(:)%region == i))
  enddo !i
  if (z_scaling_env < 1.) regions(:)%Rmin = regions(:)%Rmin * z_scaling_env
  rmin = minval(regions(:)%Rmin)
  write(*,*) "New minimum radius = ", rmin, "AU"

end subroutine set_sublimation_radius

!***********************************************************

subroutine read_sublimation_radius()

  real(kind=dp) :: sublimation_radius
  integer :: ios
  character(len=512) :: filename

  filename = trim(root_dir)//"/"//trim(seed_dir)//"/data_th/"//trim(sublimationfile)
  write(*,*) "Reading sublimation file : "//trim(filename)

  open(unit=1,file=trim(filename), status="old", iostat=ios)
  if (ios /= 0) then
     call error("read_sublimation_radius: cannot open "//trim(filename))
  endif
  ! one optically thin sublimation radius per star, as written by compute_othin_sublimation_radius
  read(1,*,iostat=ios) star(:)%othin_sublimation_radius
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

subroutine equilibre_hydrostatique(density_residual, lapply)
  ! Compute hydrostatic equilibrium for the gas, then re-derive the dust
  ! vertical distribution consistently with the new gas structure
  ! (see redistribute_dust_hydrostatic in density.f90).
  !
  ! density_residual returns the 95th percentile of |ln(rho_target/rho)| over
  ! the cells denser than hydro_residual_floor x their column midplane
  ! density, ie how far the current density is from the hydrostatic profile
  ! set by the current Tdust. With lapply=.false. the density is left
  ! untouched and only the residual is computed (convergence test).
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

  real(kind=dp), intent(out), optional :: density_residual
  logical, intent(in), optional :: lapply

  real(kind=dp), dimension(nz) :: rho, ln_rho, Tsmooth
  real(kind=dp) :: dz_local, dTdz, fac1, fac2, M_stars, M_mol, total_sum, cst, rho_target, rho_mid
  real(sp), dimension(:), allocatable :: residual
  integer :: i,j, k, icell, alloc_status, n_res
  logical :: lupdate

  real(kind=dp), parameter :: hydro_relax = 0.3_dp ! under-relaxation factor for the gas density update
  real(kind=dp), parameter :: min_decline_fraction = 0.3_dp ! floor on the vertical decline rate, as a fraction of the gravity-only rate

  lupdate = .true.
  if (present(lapply)) lupdate = lapply

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

  allocate(residual(n_cells), stat=alloc_status)
  if (alloc_status > 0) call error("Allocation error residual")
  n_res = 0

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
           icell = cell_map(i,j,k)
           dz_local = z_grid(icell) - z_grid(cell_map(i,j-1,k)) ! Distance between cell centers (exact for both constant and variable dz)
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
           ! Distance to the hydrostatic profile, restricted to the cells
           ! that carry the column's mass and optical depth : the tenuous
           ! upper atmosphere is too noisy to converge cell by cell.
           rho_mid = gas_density(cell_map(i,1,k))
           do j = 1, nz
              icell = cell_map(i,j,k)
              rho_target = rho(j) * Sigma_gas_target(i,k) / total_sum
              if ((gas_density(icell) > hydro_residual_floor * rho_mid) .and. (rho_target > 0.0_dp)) then
                 n_res = n_res + 1
                 residual(n_res) = real(abs(log(rho_target / gas_density(icell))), sp)
              endif
           enddo

           ! Under-relaxation : Tdust is only known from a noisy MC estimate,
           ! and dTdz feeds into an exponential, cumulative integral, so a
           ! small temperature fluctuation (worst near the inner edge, where
           ! the stellar gravity term fac1 is largest) can otherwise produce
           ! a wildly wrong density jump on a single iteration. Blending with
           ! the previous profile damps this without biasing the converged
           ! solution : both profiles integrate to the same Sigma_gas_target,
           ! so the blend does too.
           if (lupdate) then
              do j = 1, nz
                 icell = cell_map(i,j,k)
                 gas_density(icell) = (1.0_dp-hydro_relax) * gas_density(icell) + &
                      hydro_relax * (rho(j) * Sigma_gas_target(i,k) / total_sum)
              enddo
              gas_density_midplane(i,k) = gas_density(cell_map(i,1,k))
           endif
        endif

     enddo !i
  enddo !k

  if (present(density_residual)) then
     if (n_res > 0) then
        density_residual = find_kth_smallest(max(1,nint(0.95*n_res)), residual(1:n_res))
     else
        density_residual = 0.0_dp
     endif
  endif
  deallocate(residual)

  if (.not.lupdate) return

  do icell=1, n_cells
     gas_mass(icell) = gas_density(icell) * volume(icell) * mu_mH * AU3_to_m3
  enddo

  call redistribute_dust_hydrostatic()

  call print_hydrostatic_scale_height()

  return

end subroutine equilibre_hydrostatique

!**********************************************************

subroutine free_hydrostatic()
  ! Releases the hydrostatic state so that a subsequent run in the same
  ! process takes a fresh surface density snapshot

  if (allocated(Sigma_gas_target)) deallocate(Sigma_gas_target)

  return

end subroutine free_hydrostatic


!**********************************************************

subroutine compute_column_scale_heights(i, H68, Hrms, Hmid, Hth, r_col)
  ! Computes vertical scale heights for radial column i:
  !   H68  : height enclosing 68.27% of the mass column
  !   Hrms : root-mean-square height sqrt(<z^2>)
  !   Hmid : isothermal sound speed scale height cs(Tmid)/OmegaK
  !   Hth  : sound speed scale height from column mass-weighted temperature cs(<T>)/OmegaK
  ! C. Pinte, 08/2026

  integer, intent(in) :: i
  real(kind=dp), intent(out) :: H68, Hrms, Hmid, Hth, r_col

  real(kind=dp) :: total_mass, target_mass, cum_mass, prev_mass
  real(kind=dp) :: sum_z2, sum_rho_T, z_cell, rho_val, T_val, T_mid
  real(kind=dp) :: cst, M_stars, M_mol
  integer :: j, icell

  r_col = r_grid(cell_map(i,1,1))
  total_mass = 0.0_dp
  sum_z2 = 0.0_dp
  sum_rho_T = 0.0_dp

  do j = 1, nz
     icell = cell_map(i,j,1)
     rho_val = gas_density(icell)
     T_val = Tdust(icell)
     total_mass = total_mass + rho_val * cell_height(i,j)
     z_cell = z_grid(icell)
     sum_z2 = sum_z2 + rho_val * cell_height(i,j) * z_cell**2
     sum_rho_T = sum_rho_T + rho_val * cell_height(i,j) * T_val
  enddo

  if (total_mass <= 0.0_dp) then
     H68 = 0.0_dp
     Hrms = 0.0_dp
     Hmid = 0.0_dp
     Hth = 0.0_dp
     return
  endif

  ! 1. H_68: height enclosing 68.27% of the mass column
  target_mass = 0.682689492_dp * total_mass
  cum_mass = 0.0_dp
  H68 = z_lim(i,nz+1)
  do j = 1, nz
     icell = cell_map(i,j,1)
     rho_val = gas_density(icell)
     prev_mass = cum_mass
     cum_mass = cum_mass + rho_val * cell_height(i,j)
     if (cum_mass >= target_mass) then
        if (rho_val > 0.0_dp) then
           H68 = z_lim(i,j) + (target_mass - prev_mass) / rho_val
        else
           H68 = z_lim(i,j)
        endif
        exit
     endif
  enddo

  ! 2. H_rms: second moment
  Hrms = sqrt(sum_z2 / total_mass)

  ! 3. Physical constants for thermal scale heights
  M_stars = sum(star(:)%M) * Msun_to_kg
  M_mol = mu_mH * g_to_kg
  cst = Ggrav * M_stars * M_mol / (kb * AU_to_m)

  ! 4. H_mid: isothermal scale height from midplane temperature
  T_mid = Tdust(cell_map(i,1,1))
  if (cst > 0.0_dp .and. T_mid > 0.0_dp) then
     Hmid = sqrt(T_mid * r_col**3 / cst)
  else
     Hmid = 0.0_dp
  endif

  ! 5. H_th: sound-speed scale height from mass-weighted temperature
  if (cst > 0.0_dp .and. sum_rho_T > 0.0_dp) then
     Hth = sqrt((sum_rho_T / total_mass) * r_col**3 / cst)
  else
     Hth = 0.0_dp
  endif
end subroutine compute_column_scale_heights


!**********************************************************

subroutine get_interpolated_scale_heights(R_target, H68, Hrms, Hmid, Hth)
  ! Radial interpolation of scale heights to an arbitrary radius R_target
  ! C. Pinte, 08/2026

  real(kind=dp), intent(in) :: R_target
  real(kind=dp), intent(out) :: H68, Hrms, Hmid, Hth

  real(kind=dp) :: r1, r2, w, r_dummy
  real(kind=dp) :: H68_1, Hrms_1, Hmid_1, Hth_1
  real(kind=dp) :: H68_2, Hrms_2, Hmid_2, Hth_2
  integer :: i

  if (R_target <= r_grid(cell_map(1,1,1))) then
     call compute_column_scale_heights(1, H68, Hrms, Hmid, Hth, r_dummy)
     return
  else if (R_target >= r_grid(cell_map(n_rad,1,1))) then
     call compute_column_scale_heights(n_rad, H68, Hrms, Hmid, Hth, r_dummy)
     return
  endif

  ! Find radial index i such that r(i) <= R_target <= r(i+1)
  do i = 1, n_rad - 1
     r1 = r_grid(cell_map(i,1,1))
     r2 = r_grid(cell_map(i+1,1,1))
     if (R_target >= r1 .and. R_target <= r2) exit
  enddo

  call compute_column_scale_heights(i, H68_1, Hrms_1, Hmid_1, Hth_1, r1)
  call compute_column_scale_heights(i+1, H68_2, Hrms_2, Hmid_2, Hth_2, r2)

  if (r2 > r1) then
     w = (log(R_target) - log(r1)) / (log(r2) - log(r1))
  else
     w = 0.5_dp
  endif

  H68 = exp((1.0_dp - w) * log(max(1.0e-30_dp, H68_1)) + w * log(max(1.0e-30_dp, H68_2)))
  Hrms = exp((1.0_dp - w) * log(max(1.0e-30_dp, Hrms_1)) + w * log(max(1.0e-30_dp, Hrms_2)))
  Hmid = exp((1.0_dp - w) * log(max(1.0e-30_dp, Hmid_1)) + w * log(max(1.0e-30_dp, Hmid_2)))
  Hth = exp((1.0_dp - w) * log(max(1.0e-30_dp, Hth_1)) + w * log(max(1.0e-30_dp, Hth_2)))
end subroutine get_interpolated_scale_heights


!**********************************************************

subroutine print_hydrostatic_scale_height()
  ! Prints hydrostatic scale heights at Rref and the local flaring index
  ! C. Pinte, 08/2026

  real(kind=dp) :: Rref, H0_init, R1, R2
  real(kind=dp) :: H68, Hrms, Hmid, Hth
  real(kind=dp) :: H68_R1, Hrms_R1, Hmid_R1, Hth_R1
  real(kind=dp) :: H68_R2, Hrms_R2, Hmid_R2, Hth_R2
  real(kind=dp) :: beta_local, diff_pct
  real(kind=dp) :: r_min_grid, r_max_grid

  if (lVoronoi .or. lspherical .or. l3D) return
  if (.not. allocated(disk_zone)) return
  if (n_zones < 1) return

  Rref = disk_zone(1)%Rref
  H0_init = disk_zone(1)%sclht

  if (Rref <= 0.0_dp) Rref = 100.0_dp

  ! 1. Scale heights at Rref
  call get_interpolated_scale_heights(Rref, H68, Hrms, Hmid, Hth)

  ! 2. Local flaring index between ~ 0.9 Rref and 1.1 Rref
  r_min_grid = r_grid(cell_map(1,1,1))
  r_max_grid = r_grid(cell_map(n_rad,1,1))
  R1 = max(r_min_grid, 0.9_dp * Rref)
  R2 = min(r_max_grid, 1.1_dp * Rref)

  if (R2 > R1) then
     call get_interpolated_scale_heights(R1, H68_R1, Hrms_R1, Hmid_R1, Hth_R1)
     call get_interpolated_scale_heights(R2, H68_R2, Hrms_R2, Hmid_R2, Hth_R2)
     if (H68_R1 > 0.0_dp .and. H68_R2 > 0.0_dp) then
        beta_local = (log(H68_R2) - log(H68_R1)) / (log(R2) - log(R1))
     else
        beta_local = 0.0_dp
     endif
  else
     beta_local = 0.0_dp
  endif

  ! 3. Difference percentage between H_rms and H_thermal
  if (Hrms > 0.0_dp) then
     diff_pct = abs(Hrms - Hth) / Hrms * 100.0_dp
  else
     diff_pct = 0.0_dp
  endif

  ! 4. Display 2-line summary requested by user
  write(*,'("  Scale height @ Rref (", f6.1, " AU): H_68 = ", f6.2, " AU, H_mid = ", f6.2, &
       & " AU (initial H0 = ", f6.2, " AU), flaring beta = ", f5.2)') &
       Rref, H68, Hmid, H0_init, beta_local

  write(*,'("  Hydrostatic check: H_rms = ", f6.2, " AU, H_thermal(mean_T) = ", f6.2, &
       & " AU (diff = ", f5.2, "%)")') &
       Hrms, Hth, diff_pct

end subroutine print_hydrostatic_scale_height

!**********************************************************

subroutine init_sublimation()
  ! Takes the current dust density as the pristine (unsublimated) density

  integer :: k

  call free_sublimation()

  allocate(dust_density_pristine(size(dust_density_o_n_grains,1),n_cells), &
       lsublimated(size(dust_density_o_n_grains,1),n_cells), &
       n_sublimation_flips(size(dust_density_o_n_grains,1),n_cells), &
       pop_of_pk(size(dust_density_o_n_grains,1)))
  dust_density_pristine(:,:) = dust_density_o_n_grains(:,1:n_cells)
  lsublimated(:,:) = .false.
  n_sublimation_flips(:,:) = 0_1

  ! first index of dust_density_o_n_grains : grain size (lvariable_dust) or zone.
  ! In zone mode there is a single population per zone (lvariable_dust is set otherwise).
  pop_of_pk(:) = 0
  do k=1, n_grains_tot
     pop_of_pk(merge(k, grain(k)%zone, lvariable_dust)) = grain(k)%pop
  enddo

  return

end subroutine init_sublimation

!**********************************************************

subroutine set_pristine_dust_density()
  ! The dust density has been recomputed (hydrostatic redistribution) :
  ! it becomes the new pristine density, before the sublimation mask is applied

  dust_density_pristine(:,:) = dust_density_o_n_grains(:,1:n_cells)

  return

end subroutine set_pristine_dust_density

!**********************************************************

subroutine update_sublimation(n_new, n_recondensed, n_frozen, lapply)
  ! Compares the temperature of each (grain size or zone, cell) with the sublimation
  ! temperature of its population, with a hysteresis band set by the Monte Carlo noise :
  !  - not sublimated and T > T_sub (1+delta) : sublimates
  !  - sublimated and T < T_sub (1-delta) : re-condenses (T is then the test-particle
  !    temperature of the trace dust left in the cell)
  ! Grains that already changed state sublimation_max_flips times are frozen.
  ! n_new / n_recondensed return the number of changes, n_frozen the number of frozen
  ! grains that would otherwise change. With lapply=.false., the mask is not modified.
  ! LTE grains : one temperature per cell; nLTE grains and nRE grains at equilibrium :
  ! one temperature per grain size. PAHs are not sublimated.

  use radiation_field, only : xN_abs

  integer, intent(out) :: n_new, n_recondensed, n_frozen
  logical, intent(in) :: lapply

  integer :: icell, p_k, pop, change
  real :: delta, n_packets, T_grain, rho_gas

  n_new = 0 ; n_recondensed = 0 ; n_frozen = 0
  do icell=1, n_cells
     n_packets = 0.0
     if (allocated(xN_abs)) n_packets = sum(xN_abs(icell,1,:))
     delta = sublimation_band(n_packets)
     rho_gas = real(gas_density(icell) * mu_mH / m3_to_cm3) ! g/cm3

     do p_k=1, size(lsublimated,1)
        pop = pop_of_pk(p_k)
        if (pop == 0) cycle
        if (dust_pop(pop)%is_PAH) cycle
        if (dust_density_pristine(p_k,icell) <= 0.0_dp) cycle

        ! Temperature of the grain : one per cell for LTE grains, one per grain size for
        ! nLTE grains and for nRE grains at equilibrium (p_k is the grain index, as
        ! lvariable_dust is set for nLTE / nRE grains). Stochastically heated grains have
        ! no single temperature : their state is not changed.
        T_grain = Tdust(icell)
        if (lvariable_dust) then
           select case(dust_pop(pop)%methode_chauffage)
           case(2)
              T_grain = Tdust_1grain(p_k,icell)
           case(3)
              if (.not.l_RE(p_k,icell)) cycle
              T_grain = Tdust_1grain_nRE(p_k,icell)
           end select
        endif

        change = sublimation_change(lsublimated(p_k,icell), T_grain, population_sublimation_temperature(pop, rho_gas), delta, &
             int(n_sublimation_flips(p_k,icell)))
        select case(change)
        case(sublimation_sublimates)
           n_new = n_new + 1
        case(sublimation_recondenses)
           n_recondensed = n_recondensed + 1
        case(sublimation_frozen)
           n_frozen = n_frozen + 1
           cycle
        case default
           cycle
        end select

        if (lapply) then
           lsublimated(p_k,icell) = .not.lsublimated(p_k,icell)
           n_sublimation_flips(p_k,icell) = n_sublimation_flips(p_k,icell) + 1_1
        endif
     enddo
  enddo

  return

end subroutine update_sublimation

!**********************************************************

pure real function sublimation_temperature(law, T_sub_constant, T0, rho0, gamma, rho_gas)
  ! Sublimation temperature of a dust component [K] for a gas density rho_gas [g/cm3] :
  !  - sublimation_law_constant  : T_sub_constant (from the optical indices file)
  !  - sublimation_law_power_law : T0 (rho_gas / rho0)^gamma (Pollack et al. 1994,
  !    Isella & Natta 2005), with rho_gas floored at sublimation_rho_floor

  integer, intent(in) :: law
  real, intent(in) :: T_sub_constant, T0, rho0, gamma, rho_gas

  select case(law)
  case(sublimation_law_power_law)
     sublimation_temperature = T0 * (max(rho_gas, sublimation_rho_floor) / rho0)**gamma
  case default
     sublimation_temperature = T_sub_constant
  end select

end function sublimation_temperature

!**********************************************************

real function population_sublimation_temperature(pop, rho_gas)
  ! Sublimation temperature of a population : the lowest one of its components
  ! (as for the constant T_sub in dust_prop)

  integer, intent(in) :: pop
  real, intent(in) :: rho_gas

  integer :: k

  if (dust_pop(pop)%is_opacity_file) then
     population_sublimation_temperature = dust_pop(pop)%T_sub
     return
  endif

  population_sublimation_temperature = huge(1.0)
  do k=1, dust_pop(pop)%n_components
     population_sublimation_temperature = min(population_sublimation_temperature, &
          sublimation_temperature(dust_pop(pop)%component_sub_law(k), dust_pop(pop)%component_T_sub(k), &
          dust_pop(pop)%component_sub_T0(k), dust_pop(pop)%component_sub_rho0(k), &
          dust_pop(pop)%component_sub_gamma(k), rho_gas))
  enddo

end function population_sublimation_temperature

!**********************************************************

pure real function sublimation_band(n_packets)
  ! Relative half-width delta of the hysteresis band around T_sub, from the Monte
  ! Carlo noise on the temperature : T ~ E^(1/4), relative noise on E ~ 1/sqrt(N),
  ! so sigma_T/T ~ 0.25/sqrt(N). Bounded by sublimation_precision and
  ! sublimation_band_max; sublimation_band_max when there are too few packets.

  real, intent(in) :: n_packets

  if (n_packets >= sublimation_min_packets) then
     sublimation_band = min(max(sublimation_precision, sublimation_n_sigma * 0.25 / sqrt(n_packets)), &
          sublimation_band_max)
  else
     sublimation_band = sublimation_band_max
  endif

end function sublimation_band

!**********************************************************

pure integer function sublimation_change(lsublimated_in, T, T_sub, delta, n_flips)
  ! State change of a grain (or zone) in a cell, given its temperature :
  !  sublimation_sublimates  : not sublimated and T > T_sub (1+delta)
  !  sublimation_recondenses : sublimated and T < T_sub (1-delta)
  !  sublimation_frozen      : would change, but already changed sublimation_max_flips times
  !  sublimation_no_change   : otherwise

  logical, intent(in) :: lsublimated_in
  real, intent(in) :: T, T_sub, delta
  integer, intent(in) :: n_flips

  logical :: lchange

  if (lsublimated_in) then
     lchange = (T < T_sub * (1.0 - delta))
  else
     lchange = (T > T_sub * (1.0 + delta))
  endif

  if (.not.lchange) then
     sublimation_change = sublimation_no_change
  else if (n_flips >= sublimation_max_flips) then
     sublimation_change = sublimation_frozen
  else if (lsublimated_in) then
     sublimation_change = sublimation_recondenses
  else
     sublimation_change = sublimation_sublimates
  endif

end function sublimation_change

!**********************************************************

subroutine apply_sublimation()
  ! dust density = pristine density where the grains have not sublimated, and a trace
  ! (sublimation_trace) of it where they have; then recomputes the dust mass of each
  ! cell and the reference (first non-empty) cell

  integer :: icell, k, pop
  real(kind=dp) :: mass

  where (lsublimated)
     dust_density_o_n_grains(:,1:n_cells) = dust_density_pristine * sublimation_trace
  elsewhere
     dust_density_o_n_grains(:,1:n_cells) = dust_density_pristine
  end where

  do icell=1, n_cells
     dust_mass(icell) = 0.0_dp
     if (lvariable_dust) then
        do k=1, n_grains_tot
           dust_mass(icell) = dust_mass(icell) + dust_density_o_n_grains(k,icell) * n_grains(k) * M_grain(k)
        enddo
     else
        ! zone mode : dust_density_o_n_grains(izone,:) is the mass density of the zone,
        ! split between its populations by frac_mass (as in normalize_dust_density)
        do pop=1, n_pop
           dust_mass(icell) = dust_mass(icell) + dust_density_o_n_grains(dust_pop(pop)%zone,icell) * dust_pop(pop)%frac_mass
        enddo
     endif
     dust_mass(icell) = dust_mass(icell) * volume(icell) * AU3_to_cm3
  enddo

  call find_non_empty_cell()
  mass = sum(dust_mass) * g_to_Msun
  write(*,'(" Sublimation : ", i0, " sublimated (grain or zone, cell) pairs, dust mass = ", es10.3, " Msun")') &
       count(lsublimated), mass
  if (all(lsublimated .or. (dust_density_pristine <= 0.0_dp))) call error("Sublimation : all the dust has sublimated")

  return

end subroutine apply_sublimation

!**********************************************************

subroutine write_sublimation_mask(iter)
  ! Saves the sublimation mask (1 = sublimated) in data_th, for restart and
  ! for image/SED runs. iter > 0 : per-iteration file, iter = 0 : final mask.

  use fits_utils, only : cfitsWrite

  integer, intent(in) :: iter

  character(len=512) :: filename
  character(len=4) :: iter_str

  if (iter > 0) then
     write(iter_str,'(i4.4)') iter
     filename = "!"//trim(data_dir)//"/"//trim(sublimationMaskFile)//"_iter"//iter_str//".fits.gz"
  else
     filename = "!"//trim(data_dir)//"/"//trim(sublimationMaskFile)//".fits.gz"
  endif
  call cfitsWrite(trim(filename), merge(1.0, 0.0, lsublimated), shape(lsublimated))

  return

end subroutine write_sublimation_mask

!**********************************************************

subroutine read_sublimation_mask()
  ! Reads the final sublimation mask written by the temperature run

  integer :: status, unit, blocksize, readwrite, nfound, group, firstpix, npixels
  integer, dimension(2) :: naxes
  real :: nullval
  logical :: anynull
  real, dimension(:), allocatable :: mask
  character(len=512) :: filename

  if (.not.allocated(lsublimated)) call error("read_sublimation_mask: init_sublimation must be called first")

  filename = trim(root_dir)//"/"//trim(seed_dir)//"/data_th/"//trim(sublimationMaskFile)//".fits.gz"
  write(*,*) "Reading sublimation mask : "//trim(filename)

  status = 0
  call ftgiou(unit,status)
  readwrite = 0
  call ftopen(unit,trim(filename),readwrite,blocksize,status)
  if (status /= 0) call error("cannot open "//trim(filename), &
       msg2="the temperature run must be done with this version of mcfost")
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  if ((status /= 0).or.(nfound /= 2)) call error("cannot read the dimensions of "//trim(filename))
  if ((naxes(1) /= size(lsublimated,1)).or.(naxes(2) /= size(lsublimated,2))) then
     call error("sublimation mask does not match the current grid and dust populations")
  endif
  npixels = naxes(1)*naxes(2)
  allocate(mask(npixels))
  group = 1 ; firstpix = 1 ; nullval = -999
  call ftgpve(unit,group,firstpix,npixels,nullval,mask,anynull,status)
  if (status /= 0) call error("cannot read "//trim(filename))
  call ftclos(unit,status)
  call ftfiou(unit,status)

  lsublimated = reshape(mask > 0.5, shape(lsublimated))
  deallocate(mask)

  return

end subroutine read_sublimation_mask

!**********************************************************

subroutine free_sublimation()

  if (allocated(dust_density_pristine)) deallocate(dust_density_pristine)
  if (allocated(lsublimated)) deallocate(lsublimated)
  if (allocated(n_sublimation_flips)) deallocate(n_sublimation_flips)
  if (allocated(pop_of_pk)) deallocate(pop_of_pk)

  return

end subroutine free_sublimation

!**********************************************************

end module disk_physics
