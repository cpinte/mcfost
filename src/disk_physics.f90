module disk_physics

  use grains
  use mcfost_env
  use dust_prop
  use density, only : dust_density_o_n_grains, dust_mass, find_non_empty_cell
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
  real :: cst, wl, delta_wl, star_flux
  integer :: lambda, icell, i

  E_dust = 0.0
  cst=thermal_const/dust_pop(1)%T_sub

  icell = icell1

  do lambda=1, n_lambda
     ! wavelength in metres
     wl = tab_lambda(lambda)*1.e-6
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     cst_wl=cst/wl
     if (cst_wl < 500.0) then
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
  ! Compute hydrostatic equilibrium for each radius
  ! Equation 2.4.3 of the thesis (page 38, 52 of the PDF)
  ! Valid for a geometrically thin, non-self-gravitating perfect gas disk
  !
  ! C. Pinte
  ! 25/09/07

  real, dimension(nz) :: rho, ln_rho
  real :: dz, dz_m1, dTdz, fac1, fac2, M_stars, M_mol, total_sum, cst
  integer :: i,j, k, icell, icell_m1

  real, parameter :: gas_dust = 100

  M_stars = sum(star(:)%M) * Msun_to_kg
  M_mol = mu_mH * g_to_kg

  cst = Ggrav * M_stars * M_mol / (kb * AU_to_m**2)

  do k=1, n_az
     do i=1, n_rad
        ln_rho(1) = 0.
        rho(1) = 1.
        dz_m1 = 1.0/dz
        total_sum = rho(1)
        do j = 2, nz
           dz = cell_height(i,j)
           icell = cell_map(i,j,k)
           icell_m1 = cell_map(i,j-1,k)
           dTdz = (Tdust(icell)-Tdust(icell_m1)) * dz_m1
           fac1 = cst * z_grid(icell)/ (r_grid(icell)**3)
           fac2 = -1.0 * (dTdz + fac1) / Tdust(icell)
           ln_rho(j) = ln_rho(j-1) + fac2 * dz
           rho(j) = exp(ln_rho(j))
           total_sum = total_sum + rho(j)
        enddo !j

        ! Renormalisation
        do j = 1, nz
          ! icell = cell_map(i,j,k)
          ! fac = gas_dust * masse_rayon(i,k) / (volume(icell) * total_sum) ! TODO : densite est en particule, non ???
          ! gas_density(icell) =  rho(j) * fac
        enddo

     enddo !i
  enddo !k

  return

end subroutine equilibre_hydrostatique

!**********************************************************

end module disk_physics
