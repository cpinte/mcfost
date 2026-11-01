program test_mrw
  ! Unit test: Modified Random Walk (MRW) diffusion kernel (MRW.f90) and
  ! nearest-wall Euclidean distance routines on cylindrical, spherical, and
  ! 3D Voronoi grids (cylindrical_grid.f90, spherical_grid.f90, Voronoi.f90).
  !
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use mcfost_env, only : dp, nb_proc
  use constants, only : pi
  use parameters
  use naleat, only : init_random_streams
  use MRW, only : n, zeta, y_MRW, gamma_MRW, cumulative_zeta_from_y, &
       initialize_cumulative_zeta, interp_y_from_zeta, make_MRW_step, &
       sample_lambertian_direction
  use cylindrical_grid, only : build_cylindrical_cell_mapping, define_cylindrical_grid, &
       cross_cylindrical_cell, index_cell_cyl, distance_to_closest_wall_cyl, &
       cell2cylindrical, cell_map, cell_map_i, cell_map_j, cell_map_k, lexit_cell, &
       r_lim, r_lim_2, r_lim_3, cell_height, r_grid, z_grid, phi_grid, region_map, &
       z_lim, w_lim, theta_lim, tan_theta_lim, cos_theta_lim, tan_phi_lim, &
       cos_phi_lim, sin_phi_lim, volume, zmax
  use spherical_grid, only : cross_spherical_cell, index_cell_sph, distance_to_closest_wall_sph
  use Voronoi_grid, only : Voronoi, Voronoi_xyz, Voronoi_neighbour_xyz, neighbours_list, &
       wall, n_walls, voro, init_Platonic_Solid, PS, cross_Voronoi_cell, &
       distance_to_closest_wall_Voronoi
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use, intrinsic :: iso_c_binding, only : c_int, c_double, c_bool

  implicit none

  integer, parameter :: n_dirs = 1000
  real(dp), dimension(3, n_dirs) :: fib_dirs

  integer :: n_fail, n_test

  n_fail = 0
  n_test = 0

  call init_fibonacci_sphere(fib_dirs)
  call init_rng()

  call test_mrw_kernel()
  call test_cylindrical_walls()
  call test_spherical_walls()
  call test_voronoi_walls(50.0_dp)
  call test_voronoi_walls(0.05_dp)

  write(*,'(a,i0,a,i0)') "test_mrw: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  subroutine check(name, cond, msg)
    character(len=*), intent(in) :: name
    logical, intent(in) :: cond
    character(len=*), intent(in), optional :: msg

    n_test = n_test + 1
    if (.not. cond) then
       n_fail = n_fail + 1
       if (present(msg)) then
          write(*,'(a,a,a,a)') "FAIL ", name, ": ", msg
       else
          write(*,'(a,a)') "FAIL ", name
       endif
    endif
  end subroutine check

  subroutine init_rng()
    nb_proc = 1
    call init_random_streams(1)
  end subroutine init_rng

  subroutine init_fibonacci_sphere(dirs)
    ! Quasi-uniform Fibonacci lattice on S^2 plus the 6 cardinal axis directions
    real(dp), dimension(3, n_dirs), intent(out) :: dirs
    real(dp), parameter :: golden_angle = pi * (3.0_dp - sqrt(5.0_dp))
    integer :: i, n_fib
    real(dp) :: z_i, r_i, phi_i

    dirs(:,1) = [ 1.0_dp,  0.0_dp,  0.0_dp]
    dirs(:,2) = [-1.0_dp,  0.0_dp,  0.0_dp]
    dirs(:,3) = [ 0.0_dp,  1.0_dp,  0.0_dp]
    dirs(:,4) = [ 0.0_dp, -1.0_dp,  0.0_dp]
    dirs(:,5) = [ 0.0_dp,  0.0_dp,  1.0_dp]
    dirs(:,6) = [ 0.0_dp,  0.0_dp, -1.0_dp]

    n_fib = n_dirs - 6
    do i = 1, n_fib
       z_i = 1.0_dp - 2.0_dp * (real(i, dp) - 0.5_dp) / real(n_fib, dp)
       r_i = sqrt(max(0.0_dp, 1.0_dp - z_i * z_i))
       phi_i = golden_angle * real(i, dp)
       dirs(:, i + 6) = [r_i * cos(phi_i), r_i * sin(phi_i), z_i]
    enddo
  end subroutine init_fibonacci_sphere

  !======================================================================
  ! 1. MRW analytical diffusion kernel & step mechanics
  !======================================================================
  subroutine test_mrw_kernel()
    integer, parameter :: n_mc = 80000, n_quad = 200000, n_chain = 4000
    integer :: i
    real(dp) :: z_direct, z_jacobi, term, max_inv_err, xi, y_val
    real(dp) :: quad_mean, exact_mean_lny
    real(dp) :: x0, y0, z0, x, y, z, u, v, w, ct, R0, chi
    real(dp) :: dx, dy, dz, step_r, norm_u, mu
    real(dp) :: sum_ct, sum_mu, sum_mu2, sum_dx, sum_dy, sum_dz
    real(dp) :: R_macro, d_rem, r_cur, sum_r2_stop
    real(dp) :: max_r_err, max_u_err, min_mu
    logical :: all_finite
    character(len=128) :: msg

    ! Continuity between direct series and Jacobi theta_4 identity at y = 0.5
    z_direct = 0.0_dp
    do i = 1, 20
       term = 0.5_dp**(i * i)
       if (mod(i, 2) == 0) then
          z_direct = z_direct - term
       else
          z_direct = z_direct + term
       endif
    enddo
    z_direct = 2.0_dp * z_direct
    z_jacobi = cumulative_zeta_from_y(0.5_dp + 1.0e-15_dp)
    write(msg,'(a,es11.4)') "err = ", abs(z_direct - z_jacobi)
    call check("MRW zeta continuity at y=0.5", abs(z_direct - z_jacobi) < 1.0e-13_dp, trim(msg))

    call initialize_cumulative_zeta()

    max_inv_err = 0.0_dp
    do i = 1, n
       max_inv_err = max(max_inv_err, abs(cumulative_zeta_from_y(y_MRW(i)) - zeta(i)))
    enddo
    write(msg,'(a,es11.4)') "max |zeta(y_MRW) - zeta| = ", max_inv_err
    call check("MRW inverse zeta table accuracy", max_inv_err < 1.0e-13_dp, trim(msg))

    ! Endpoint FPE safety (zeta = 0, zeta = 1, zeta > 1)
    y_val = interp_y_from_zeta(0.0_dp)
    call check("MRW zeta=0 strictly positive & finite -log(y)", &
         (y_val > 0.0_dp) .and. ieee_is_finite(-log(y_val)))
    y_val = interp_y_from_zeta(1.0_dp)
    call check("MRW zeta=1 finite -log(y)", &
         (y_val > 0.0_dp) .and. (y_val <= 1.0_dp) .and. ieee_is_finite(-log(y_val)))

    ! Midpoint quadrature of <-ln(y)> against analytical Dirichlet eta(2)*2 = pi^2 / 6
    quad_mean = 0.0_dp
    do i = 1, n_quad
       xi = (real(i, dp) - 0.5_dp) / real(n_quad, dp)
       quad_mean = quad_mean - log(interp_y_from_zeta(xi))
    enddo
    quad_mean = quad_mean / real(n_quad, dp)
    exact_mean_lny = (pi * pi) / 6.0_dp
    write(msg,'(a,es12.5,a,es12.5)') "got ", quad_mean, " expected ", exact_mean_lny
    call check("MRW quadrature <-ln(y)> = pi^2/6", &
         abs(quad_mean - exact_mean_lny) < 5.0e-5_dp * exact_mean_lny, trim(msg))

    ! North & south pole Lambertian direction checks (including wn = -1.0)
    min_mu = 1.0_dp
    max_u_err = 0.0_dp
    do i = 1, 1000
       call sample_lambertian_direction(1, 0.0_dp, 0.0_dp, -1.0_dp, u, v, w)
       max_u_err = max(max_u_err, abs(sqrt(u*u + v*v + w*w) - 1.0_dp))
       min_mu = min(min_mu, -w)
       call sample_lambertian_direction(1, 0.0_dp, 0.0_dp, 1.0_dp, u, v, w)
       max_u_err = max(max_u_err, abs(sqrt(u*u + v*v + w*w) - 1.0_dp))
       min_mu = min(min_mu, w)
    enddo
    call check("MRW Lambertian poles unit norm & outward", &
         (max_u_err < 1.0e-14_dp) .and. (min_mu >= 0.0_dp))

    ! Monte Carlo verification of make_MRW_step (catches A11 y-coordinate overwrite,
    ! verifies <ct> = 0.5 * R0^2 * chi, isotropic sphere displacement, Lambertian moments)
    x0 = 1.5_dp ; y0 = -7.25_dp ; z0 = 3.1_dp
    R0 = 2.4_dp ; chi = 100.0_dp

    max_r_err = 0.0_dp ; max_u_err = 0.0_dp ; min_mu = 1.0_dp
    sum_ct = 0.0_dp ; sum_mu = 0.0_dp ; sum_mu2 = 0.0_dp
    sum_dx = 0.0_dp ; sum_dy = 0.0_dp ; sum_dz = 0.0_dp
    all_finite = .true.

    do i = 1, n_mc
       x = x0 ; y = y0 ; z = z0
       call make_MRW_step(1, x, y, z, u, v, w, R0, chi, ct)
       if (.not. (ieee_is_finite(x) .and. ieee_is_finite(y) .and. ieee_is_finite(z) .and. &
            ieee_is_finite(u) .and. ieee_is_finite(v) .and. ieee_is_finite(w) .and. &
            ieee_is_finite(ct) .and. (ct >= 0.0_dp))) then
          all_finite = .false.
       endif
       dx = x - x0 ; dy = y - y0 ; dz = z - z0
       step_r = sqrt(dx*dx + dy*dy + dz*dz)
       norm_u = sqrt(u*u + v*v + w*w)
       mu = (dx*u + dy*v + dz*w) / R0

       max_r_err = max(max_r_err, abs(step_r - R0))
       max_u_err = max(max_u_err, abs(norm_u - 1.0_dp))
       min_mu = min(min_mu, mu)

       sum_ct = sum_ct + ct
       sum_mu = sum_mu + mu
       sum_mu2 = sum_mu2 + mu*mu
       sum_dx = sum_dx + dx ; sum_dy = sum_dy + dy ; sum_dz = sum_dz + dz
    enddo

    call check("MRW step all outputs finite & ct >= 0", all_finite)
    write(msg,'(a,es11.4)') "max ||dr| - R0| = ", max_r_err
    call check("MRW step exact sphere displacement |dr| = R0 (A11)", max_r_err < 1.0e-13_dp, trim(msg))
    call check("MRW step direction unit norm & outward mu >= 0", &
         (max_u_err < 1.0e-13_dp) .and. (min_mu >= -1.0e-14_dp))

    write(msg,'(a,es12.5,a,es12.5)') "got ", sum_ct / n_mc, " expected ", 0.5_dp * R0 * R0 * chi
    call check("MRW step <ct> = 0.5 * R0^2 * chi", &
         abs(sum_ct / n_mc - 0.5_dp * R0 * R0 * chi) < 0.015_dp * (0.5_dp * R0 * R0 * chi), trim(msg))
    call check("MRW step isotropic <dr> ~ 0", &
         max(abs(sum_dx), abs(sum_dy), abs(sum_dz)) / (n_mc * R0) < 0.015_dp)
    call check("MRW step Lambertian <mu>=2/3 and <mu^2>=1/2", &
         (abs(sum_mu / n_mc - 2.0_dp/3.0_dp) < 0.01_dp) .and. &
         (abs(sum_mu2 / n_mc - 0.5_dp) < 0.01_dp))

    ! Multi-step chained adaptive MRW walk inside a macro-sphere of radius R_macro:
    ! Using R0 = 0.5 * (R_macro - |r|) forces ~15-30 chained steps per walk;
    ! by Dynkin's formula, E[sum ct_k] = 0.5 * E[r_stop^2] * chi.
    R_macro = 1.0_dp
    chi = 5.0e3_dp
    sum_ct = 0.0_dp
    sum_r2_stop = 0.0_dp
    do i = 1, n_chain
       x = 0.0_dp ; y = 0.0_dp ; z = 0.0_dp
       do
          r_cur = sqrt(x*x + y*y + z*z)
          d_rem = R_macro - r_cur
          if (d_rem * chi <= gamma_MRW) exit
          R0 = 0.5_dp * d_rem
          call make_MRW_step(1, x, y, z, u, v, w, R0, chi, ct)
          sum_ct = sum_ct + ct
       enddo
       sum_r2_stop = sum_r2_stop + (x*x + y*y + z*z)
    enddo
    write(msg,'(a,es12.5,a,es12.5)') "got ", sum_ct / n_chain, " expected ", 0.5_dp * (sum_r2_stop / n_chain) * chi
    call check("MRW chained multi-step <sum ct> = 0.5 * <r_stop^2> * chi", &
         abs(sum_ct / n_chain - 0.5_dp * (sum_r2_stop / n_chain) * chi) < &
         0.025_dp * (0.5_dp * (sum_r2_stop / n_chain) * chi), trim(msg))
  end subroutine test_mrw_kernel

  !======================================================================
  ! Grid setup/teardown helpers for cylindrical & spherical tests
  !======================================================================
  subroutine free_structured_grid()
    if (allocated(cell_map)) deallocate(cell_map)
    if (allocated(cell_map_i)) deallocate(cell_map_i)
    if (allocated(cell_map_j)) deallocate(cell_map_j)
    if (allocated(cell_map_k)) deallocate(cell_map_k)
    if (allocated(lexit_cell)) deallocate(lexit_cell)
    if (allocated(r_lim)) deallocate(r_lim)
    if (allocated(r_lim_2)) deallocate(r_lim_2)
    if (allocated(r_lim_3)) deallocate(r_lim_3)
    if (allocated(cell_height)) deallocate(cell_height)
    if (allocated(r_grid)) deallocate(r_grid)
    if (allocated(z_grid)) deallocate(z_grid)
    if (allocated(phi_grid)) deallocate(phi_grid)
    if (allocated(region_map)) deallocate(region_map)
    if (allocated(z_lim)) deallocate(z_lim)
    if (allocated(w_lim)) deallocate(w_lim)
    if (allocated(theta_lim)) deallocate(theta_lim)
    if (allocated(tan_theta_lim)) deallocate(tan_theta_lim)
    if (allocated(cos_theta_lim)) deallocate(cos_theta_lim)
    if (allocated(tan_phi_lim)) deallocate(tan_phi_lim)
    if (allocated(cos_phi_lim)) deallocate(cos_phi_lim)
    if (allocated(sin_phi_lim)) deallocate(sin_phi_lim)
    if (allocated(volume)) deallocate(volume)
    if (allocated(zmax)) deallocate(zmax)
    if (allocated(disk_zone)) deallocate(disk_zone)
    if (allocated(regions)) deallocate(regions)
  end subroutine free_structured_grid

  subroutine init_structured_grid(is_cyl, is_3d, nr, n_vert, naz, reg_theta)
    logical, intent(in) :: is_cyl, is_3d, reg_theta
    integer, intent(in) :: nr, n_vert, naz

    call free_structured_grid()

    lcylindrical = is_cyl
    lspherical = .not. is_cyl
    lVoronoi = .false.
    l3D = is_3d
    llinear_rgrid = .false.
    lregular_theta = reg_theta
    theta_max = 0.5 * real(pi)
    lidefix = .false.
    lsphere_model = .false.
    lmodel_1d = .false.
    lfargo3d = .false.
    lSeb_Charnoz = .false.

    n_rad = nr
    n_rad_in = 1
    nz = n_vert
    n_az = naz
    if (l3D) then
       j_start = -nz
       n_cells = 2 * n_rad * nz * n_az
    else
       j_start = 1
       n_cells = n_rad * nz
    endif

    n_zones = 1
    n_regions = 1
    allocate(disk_zone(1), regions(1))
    allocate(regions(1)%zones(1))
    regions(1)%n_zones = 1
    regions(1)%zones(1) = 1
    Rmin = 1.0_dp
    Rmax = 20.0_dp
    regions(1)%Rmin = Rmin
    regions(1)%Rmax = Rmax

    disk_zone(1)%region = 1
    disk_zone(1)%Rin = Rmin
    disk_zone(1)%Rmin = Rmin
    disk_zone(1)%Rout = Rmax
    disk_zone(1)%Rmax = Rmax
    disk_zone(1)%Rref = 5.0_dp
    disk_zone(1)%sclht = 0.8_dp
    disk_zone(1)%exp_beta = 1.15_dp
    disk_zone(1)%surf = -1.0_dp
    cutoff = 3.0

    call build_cylindrical_cell_mapping()
    call define_cylindrical_grid()
  end subroutine init_structured_grid

  !======================================================================
  ! 2. Cylindrical nearest-wall distance vs ray-sweep & sphere containment
  !======================================================================
  subroutine test_cylindrical_walls()
    call verify_cylindrical_config(.false., 4, 5, 1, "cyl 2D")
    call verify_cylindrical_config(.true.,  3, 4, 1, "cyl 3D n_az=1")
    call verify_cylindrical_config(.true.,  3, 4, 2, "cyl 3D n_az=2")
    call verify_cylindrical_config(.true.,  3, 4, 4, "cyl 3D n_az=4 (phi=pi/2,3pi/2)")
    call verify_cylindrical_config(.true.,  3, 4, 5, "cyl 3D n_az=5")
    call free_structured_grid()
  end subroutine test_cylindrical_walls

  subroutine verify_cylindrical_config(is_3d, nr, n_vert, naz, label)
    logical, intent(in) :: is_3d
    integer, intent(in) :: nr, n_vert, naz
    character(len=*), intent(in) :: label

    integer :: icell, ri, zj, ja, k, idir, c_chk, next_c, hop, n_pts
    real(dp) :: f1, f2, f3, r, z, phi, dphi, x, y, d_wall, d_ray_min
    real(dp) :: xr, yr, zr, x1, y1, z1, l_step, l_contrib, l_void, l_tot
    real(dp) :: u, v, w
    logical :: ok_contain, ok_upper, ok_tight, ok_bounds
    character(len=128) :: msg

    call init_structured_grid(.true., is_3d, nr, n_vert, naz, .false.)

    ! Virtual/invalid cell indices must safely return 0
    ok_bounds = (distance_to_closest_wall_cyl(0, 5.0_dp, 0.0_dp, 0.1_dp) == 0.0_dp) .and. &
                (distance_to_closest_wall_cyl(n_cells + 1, 5.0_dp, 0.0_dp, 0.1_dp) == 0.0_dp)
    call check(label // " invalid/virtual cell -> 0", ok_bounds)

    ok_contain = .true.
    ok_upper   = .true.
    ok_tight   = .true.
    n_pts = 0

    dphi = 2.0_dp * pi / real(n_az, dp)

    do icell = 1, n_cells
       call cell2cylindrical(icell, ri, zj, k)
       ja = abs(zj)

       ! Sample an asymmetric interior point in every cell (and near z=0 for 2D j=1)
       f1 = 0.37_dp
       f2 = 0.29_dp
       f3 = 0.23_dp
       r = r_lim(ri-1) + f1 * (r_lim(ri) - r_lim(ri-1))
       if (l3D) then
          z = z_lim(ri, ja) + f2 * (z_lim(ri, ja+1) - z_lim(ri, ja))
          if (zj < 0) z = -z
       else
          if (ja == 1) then
             ! Explicitly test very close to z = 0 in 2D midplane cell
             z = 0.05_dp * z_lim(ri, 2)
          else
             z = -(z_lim(ri, ja) + f2 * (z_lim(ri, ja+1) - z_lim(ri, ja)))
          endif
       endif
       phi = (real(k - 1, dp) + f3) * dphi
       x = r * cos(phi)
       y = r * sin(phi)

       call index_cell_cyl(x, y, z, c_chk)
       if (c_chk /= icell) cycle

       d_wall = distance_to_closest_wall_cyl(icell, x, y, z)
       if ((.not. ieee_is_finite(d_wall)) .or. (d_wall <= 0.0_dp)) then
          ok_upper = .false.
          exit
       endif

       d_ray_min = huge(1.0_dp)
       do idir = 1, n_dirs
          u = fib_dirs(1, idir) ; v = fib_dirs(2, idir) ; w = fib_dirs(3, idir)

          ! 1) Sphere of radius 0.999 * d_wall must stay strictly inside icell
          call index_cell_cyl(x + 0.999_dp*d_wall*u, y + 0.999_dp*d_wall*v, z + 0.999_dp*d_wall*w, c_chk)
          if (c_chk /= icell) ok_contain = .false.

          ! 2) Ray distance to exit icell via cross_cylindrical_cell
          xr = x ; yr = y ; zr = z
          l_tot = 0.0_dp
          do hop = 1, 4
             call cross_cylindrical_cell(xr, yr, zr, u, v, w, icell, 0, &
                  x1, y1, z1, next_c, l_step, l_contrib, l_void)
             l_tot = l_tot + l_step
             if (next_c /= icell) exit
             xr = x1 ; yr = y1 ; zr = z1
          enddo
          d_ray_min = min(d_ray_min, l_tot)
       enddo

       if (d_wall > d_ray_min * (1.0_dp + 1.0e-4_dp)) then
          write(msg,'(a,i0,a,es12.5,a,es12.5)') "cell ", icell, " d_wall=", d_wall, " > d_ray=", d_ray_min
          ok_upper = .false.
          exit
       endif
       if (d_wall < d_ray_min * 0.96_dp) then
          write(msg,'(a,i0,a,es12.5,a,es12.5)') "cell ", icell, " d_wall=", d_wall, " << d_ray=", d_ray_min
          ok_tight = .false.
          exit
       endif
       n_pts = n_pts + 1
    enddo

    call check(label // " sphere contained in cell", ok_contain .and. (n_pts == n_cells))
    call check(label // " d_wall <= min ray exit distance", ok_upper, trim(msg))
    call check(label // " d_wall tight against min ray distance", ok_tight, trim(msg))

    if (is_3d .and. (naz == 4)) then
       ! Wrong-hemisphere & outside-sector guards + chained MRW walk inside cell 1
       icell = cell_map(2, 2, 1)
       r = 0.5_dp * (r_lim(1) + r_lim(2))
       z = 0.5_dp * (z_lim(2, 2) + z_lim(2, 3))
       phi = 0.25_dp * pi
       x = r * cos(phi) ; y = r * sin(phi)
       call check(label // " wrong z hemisphere -> 0", &
            distance_to_closest_wall_cyl(icell, x, y, -z) == 0.0_dp)
       call check(label // " outside phi sector -> 0", &
            distance_to_closest_wall_cyl(icell, -x, y, z) == 0.0_dp)

       ok_contain = .true.
       do idir = 1, 200
          xr = x ; yr = y ; zr = z
          do hop = 1, 500
             d_wall = distance_to_closest_wall_cyl(icell, xr, yr, zr)
             if (d_wall * 1.0e4_dp <= gamma_MRW) exit
             call make_MRW_step(1, xr, yr, zr, u, v, w, (1.0_dp - 1.0e-6_dp) * d_wall, 1.0e4_dp, l_step)
             call index_cell_cyl(xr, yr, zr, c_chk)
             if (c_chk /= icell) then
                ok_contain = .false.
                exit
             endif
          enddo
       enddo
       call check(label // " 200 chained MRW walks stay inside cell", ok_contain)
    endif
  end subroutine verify_cylindrical_config

  !======================================================================
  ! 3. Spherical nearest-wall distance vs ray-sweep & sphere containment
  !======================================================================
  subroutine test_spherical_walls()
    ! Note: nz = 6 > n_az = 1 and n_az = 4 explicitly tests the old OOB bug
    ! where cos_phi_lim(1:n_az) was indexed by thetaj0 (up to nz, or negative in 3D).
    call verify_spherical_config(.false., 3, 6, 1, .false., "sph 2D cos-grid (nz=6>n_az=1)")
    call verify_spherical_config(.false., 3, 5, 1, .true.,  "sph 2D regular-theta")
    call verify_spherical_config(.true.,  3, 6, 1, .false., "sph 3D n_az=1 (z<0 & z>0)")
    call verify_spherical_config(.true.,  3, 6, 4, .false., "sph 3D n_az=4 (nz=6>n_az, phi=pi/2)")
    call free_structured_grid()
  end subroutine test_spherical_walls

  subroutine verify_spherical_config(is_3d, nr, n_vert, naz, reg_theta, label)
    logical, intent(in) :: is_3d, reg_theta
    integer, intent(in) :: nr, n_vert, naz
    character(len=*), intent(in) :: label

    integer :: icell, ri, tj, ja, k, idir, c_chk, next_c, hop, n_pts
    real(dp) :: f1, f2, f3, r, th, phi, dphi, rcyl, x, y, z, d_wall, d_ray_min
    real(dp) :: xr, yr, zr, x1, y1, z1, l_step, l_contrib, l_void, l_tot
    real(dp) :: u, v, w
    logical :: ok_contain, ok_upper, ok_tight, ok_bounds
    character(len=128) :: msg

    call init_structured_grid(.false., is_3d, nr, n_vert, naz, reg_theta)

    ok_bounds = (distance_to_closest_wall_sph(0, 5.0_dp, 0.0_dp, 0.1_dp) == 0.0_dp) .and. &
                (distance_to_closest_wall_sph(n_cells + 1, 5.0_dp, 0.0_dp, 0.1_dp) == 0.0_dp)
    call check(label // " invalid/virtual cell -> 0", ok_bounds)

    ok_contain = .true.
    ok_upper   = .true.
    ok_tight   = .true.
    n_pts = 0
    msg = ""

    dphi = 2.0_dp * pi / real(n_az, dp)

    do icell = 1, n_cells
       call cell2cylindrical(icell, ri, tj, k)
       ja = abs(tj)

       f1 = 0.41_dp
       f2 = 0.33_dp
       f3 = 0.27_dp
       r = r_lim(ri-1) + f1 * (r_lim(ri) - r_lim(ri-1))
       if (l3D) then
          th = theta_lim(ja-1) + f2 * (theta_lim(ja) - theta_lim(ja-1))
          if (tj < 0) th = -th
       else
          if (ja == 1) then
             th = 0.05_dp * theta_lim(1)
          else
             th = -(theta_lim(ja-1) + f2 * (theta_lim(ja) - theta_lim(ja-1)))
          endif
       endif
       phi = (real(k - 1, dp) + f3) * dphi
       z = r * sin(th)
       rcyl = r * cos(th)
       x = rcyl * cos(phi)
       y = rcyl * sin(phi)

       call index_cell_sph(x, y, z, c_chk)
       if (c_chk /= icell) cycle

       d_wall = distance_to_closest_wall_sph(icell, x, y, z)
       if ((.not. ieee_is_finite(d_wall)) .or. (d_wall <= 0.0_dp)) then
          write(msg,'(a,i0,a,es12.5)') "cell ", icell, " non-positive d_wall=", d_wall
          ok_upper = .false.
          exit
       endif

       d_ray_min = huge(1.0_dp)
       do idir = 1, n_dirs
          u = fib_dirs(1, idir) ; v = fib_dirs(2, idir) ; w = fib_dirs(3, idir)

          ! 1) Sphere of radius 0.999 * d_wall must stay strictly inside icell
          call index_cell_sph(x + 0.999_dp*d_wall*u, y + 0.999_dp*d_wall*v, z + 0.999_dp*d_wall*w, c_chk)
          if (c_chk /= icell) ok_contain = .false.

          ! 2) Ray distance to exit icell via cross_spherical_cell
          xr = x ; yr = y ; zr = z
          l_tot = 0.0_dp
          do hop = 1, 4
             call cross_spherical_cell(xr, yr, zr, u, v, w, icell, 0, &
                  x1, y1, z1, next_c, l_step, l_contrib, l_void)
             l_tot = l_tot + l_step
             if (next_c /= icell) exit
             xr = x1 ; yr = y1 ; zr = z1
          enddo
          d_ray_min = min(d_ray_min, l_tot)
       enddo

       if (d_wall > d_ray_min * (1.0_dp + 1.0e-4_dp)) then
          write(msg,'(a,i0,a,es12.5,a,es12.5)') "cell ", icell, " d_wall=", d_wall, " > d_ray=", d_ray_min
          ok_upper = .false.
          exit
       endif
       if (d_wall < d_ray_min * 0.96_dp) then
          write(msg,'(a,i0,a,es12.5,a,es12.5)') "cell ", icell, " d_wall=", d_wall, " << d_ray=", d_ray_min
          ok_tight = .false.
          exit
       endif
       n_pts = n_pts + 1
    enddo

    call check(label // " sphere contained in cell", ok_contain .and. (n_pts == n_cells))
    call check(label // " d_wall <= min ray exit distance", ok_upper, trim(msg))
    call check(label // " d_wall tight against min ray distance", ok_tight, trim(msg))

    if (is_3d .and. (naz == 4)) then
       icell = cell_map(2, -3, 3)
       r = 0.5_dp * (r_lim(1) + r_lim(2))
       th = -0.5_dp * (theta_lim(2) + theta_lim(3))
       phi = 1.25_dp * pi
       z = r * sin(th) ; rcyl = r * cos(th)
       x = rcyl * cos(phi) ; y = rcyl * sin(phi)
       call check(label // " wrong z hemisphere -> 0", &
            distance_to_closest_wall_sph(icell, x, y, -z) == 0.0_dp)
       call check(label // " outside phi sector -> 0", &
            distance_to_closest_wall_sph(icell, abs(x), abs(y), z) == 0.0_dp)

       ok_contain = .true.
       do idir = 1, 200
          xr = x ; yr = y ; zr = z
          do hop = 1, 500
             d_wall = distance_to_closest_wall_sph(icell, xr, yr, zr)
             if (d_wall * 1.0e4_dp <= gamma_MRW) exit
             call make_MRW_step(1, xr, yr, zr, u, v, w, (1.0_dp - 1.0e-6_dp) * d_wall, 1.0e4_dp, l_step)
             call index_cell_sph(xr, yr, zr, c_chk)
             if (c_chk /= icell) then
                ok_contain = .false.
                exit
             endif
          enddo
       enddo
       call check(label // " 200 chained MRW walks stay inside cell", ok_contain)
    endif
  end subroutine verify_spherical_config

  !======================================================================
  ! 4. 3D Voronoi tessellation nearest-wall distance (at arbitrary AU scale)
  !======================================================================
  subroutine test_voronoi_walls(L_scale)
    real(dp), intent(in) :: L_scale

    integer, parameter :: np = 35, max_nb = 60
    real(c_double), dimension(6) :: limits
    real(c_double), dimension(np) :: xv, yv, zv, hv, vol_c
    integer(c_int), dimension(1) :: n_nb_cpu, stars_cell_c
    real(c_double), dimension(1) :: stars_rad_c
    integer(c_int), dimension(np) :: first_nb, last_nb
    integer(c_int), dimension(np * max_nb) :: nb_list_tmp
    logical(c_bool), dimension(np) :: cut_c, star_nb_c
    integer(c_int) :: n_in, ierr

    integer :: i, j, k, idx, icell, idir, next_c, n_interior, best_cell
    real(dp) :: x, y, z, u, v, w, x1, y1, z1, l_step, l_contrib, l_void
    real(dp) :: d_wall, d_ray_min, d_star_expected
    real(dp), dimension(3) :: p_out, n_face
    logical :: has_wall_nb, ok_contain, ok_upper, ok_tight, ok_wall_zero
    character(len=64) :: label
    character(len=128) :: msg

    write(label,'(a,es8.1)') "Voronoi 3D L=", L_scale

    call free_voronoi_grid()

    lVoronoi = .true.
    lcylindrical = .false.
    lspherical = .false.
    n_cells = np

    ! Build a jittered 3x3x3 interior lattice + 7 inner points + 1 central star at index np
    idx = 0
    do k = -1, 1
       do j = -1, 1
          do i = -1, 1
             if ((i == 0) .and. (j == 0) .and. (k == 0)) cycle
             idx = idx + 1
             xv(idx) = L_scale * (0.32_dp * i + 0.03_dp * mod(idx, 5) - 0.06_dp)
             yv(idx) = L_scale * (0.32_dp * j + 0.025_dp * mod(idx, 7) - 0.07_dp)
             zv(idx) = L_scale * (0.32_dp * k + 0.02_dp * mod(idx, 3) - 0.02_dp)
          enddo
       enddo
    enddo
    ! Add extra inner shell points around origin so central cells are purely interior (no outer walls)
    xv(27:34) = L_scale * [ 0.14_dp, -0.15_dp,  0.02_dp, -0.01_dp,  0.03_dp, -0.02_dp,  0.11_dp, -0.12_dp]
    yv(27:34) = L_scale * [ 0.02_dp, -0.01_dp,  0.16_dp, -0.13_dp, -0.02_dp,  0.03_dp,  0.10_dp, -0.11_dp]
    zv(27:34) = L_scale * [-0.01_dp,  0.02_dp,  0.01_dp, -0.02_dp,  0.15_dp, -0.14_dp,  0.09_dp, -0.10_dp]
    ! Star particle at origin (cell np)
    xv(np) = 0.0_dp ; yv(np) = 0.0_dp ; zv(np) = 0.0_dp
    hv(:) = 10.0_dp * L_scale

    limits = [-L_scale, L_scale, -L_scale, L_scale, -L_scale, L_scale]

    call init_Platonic_Solid(12, 3.0_dp)

    n_stars = 1
    allocate(star(1))
    star(1)%x = 0.0_dp ; star(1)%y = 0.0_dp ; star(1)%z = 0.0_dp
    star(1)%r = 0.03_dp * L_scale
    star(1)%icell = np
    stars_cell_c(1) = np - 1 ! C 0-based index for voro_C
    stars_rad_c(1) = star(1)%r

    call voro(int(np, c_int), int(max_nb, c_int), limits, xv, yv, zv, hv, &
         1.0e10_c_double, int(PS%n_faces, c_int), PS%vectors, PS%cutting_distance_o_h, &
         0_c_int, int(np, c_int), 0_c_int, 1_c_int, int(np, c_int), &
         n_in, vol_c, first_nb, last_nb, n_nb_cpu, nb_list_tmp, cut_c, &
         1_c_int, stars_cell_c, stars_rad_c, star_nb_c, ierr)

    call check(trim(label) // " voro++ tessellation succeeded", (ierr == 0) .and. (n_in == np))

    n_walls = 6
    allocate(wall(6))
    wall(:)%type = "plane"
    wall(1)%x1 = -1.0 ; wall(1)%x2 =  0.0 ; wall(1)%x3 =  0.0 ; wall(1)%x4 = real(limits(1))
    wall(2)%x1 =  1.0 ; wall(2)%x2 =  0.0 ; wall(2)%x3 =  0.0 ; wall(2)%x4 = real(limits(2))
    wall(3)%x1 =  0.0 ; wall(3)%x2 = -1.0 ; wall(3)%x3 =  0.0 ; wall(3)%x4 = real(limits(3))
    wall(4)%x1 =  0.0 ; wall(4)%x2 =  1.0 ; wall(4)%x3 =  0.0 ; wall(4)%x4 = real(limits(4))
    wall(5)%x1 =  0.0 ; wall(5)%x2 =  0.0 ; wall(5)%x3 = -1.0 ; wall(5)%x4 = real(limits(5))
    wall(6)%x1 =  0.0 ; wall(6)%x2 =  0.0 ; wall(6)%x3 =  1.0 ; wall(6)%x4 = real(limits(6))

    allocate(Voronoi(np), Voronoi_xyz(3, np), Voronoi_neighbour_xyz(3, 0, np))
    allocate(neighbours_list(n_nb_cpu(1)))
    neighbours_list(:) = nb_list_tmp(1:n_nb_cpu(1))

    do i = 1, np
       Voronoi(i)%xyz = [xv(i), yv(i), zv(i)]
       Voronoi_xyz(:, i) = real(Voronoi(i)%xyz)
       Voronoi(i)%h = hv(i)
       Voronoi(i)%first_neighbour = first_nb(i) + 1
       Voronoi(i)%last_neighbour = last_nb(i) + 1
       Voronoi(i)%was_cut = logical(cut_c(i))
       Voronoi(i)%is_star = (i == np)
       Voronoi(i)%is_star_neighbour = .false.
       Voronoi(i)%exist = .true.
    enddo

    ! Star cell must return 0
    call check(trim(label) // " is_star -> 0", &
         distance_to_closest_wall_Voronoi(np, 0.0_dp, 0.0_dp, 0.0_dp) == 0.0_dp)

    ok_contain   = .true.
    ok_upper     = .true.
    ok_tight     = .true.
    ok_wall_zero = .true.
    n_interior   = 0
    msg = ""

    do icell = 1, np - 1
       has_wall_nb = any(neighbours_list(Voronoi(icell)%first_neighbour:Voronoi(icell)%last_neighbour) < 0)
       x = Voronoi(icell)%xyz(1) + 0.01_dp * L_scale
       y = Voronoi(icell)%xyz(2) - 0.008_dp * L_scale
       z = Voronoi(icell)%xyz(3) + 0.006_dp * L_scale

       if (has_wall_nb) then
          if (distance_to_closest_wall_Voronoi(icell, x, y, z) /= 0.0_dp) ok_wall_zero = .false.
          cycle
       endif

       ! Verify point belongs to icell
       best_cell = closest_voronoi_generator(x, y, z, np)
       if (best_cell /= icell) then
          x = Voronoi(icell)%xyz(1)
          y = Voronoi(icell)%xyz(2)
          z = Voronoi(icell)%xyz(3)
       endif

       d_wall = distance_to_closest_wall_Voronoi(icell, x, y, z)
       if ((.not. ieee_is_finite(d_wall)) .or. (d_wall <= 0.0_dp)) then
          ok_upper = .false.
          exit
       endif

       d_ray_min = huge(1.0_dp)
       do idir = 1, n_dirs
          u = fib_dirs(1, idir) ; v = fib_dirs(2, idir) ; w = fib_dirs(3, idir)

          ! 1) Sphere of radius 0.999 * d_wall must stay strictly closest to generator icell
          best_cell = closest_voronoi_generator(x + 0.999_dp*d_wall*u, &
               y + 0.999_dp*d_wall*v, z + 0.999_dp*d_wall*w, np)
          if (best_cell /= icell) ok_contain = .false.

          ! 2) Ray distance via cross_Voronoi_cell
          call cross_Voronoi_cell(x, y, z, u, v, w, icell, 0, &
               x1, y1, z1, next_c, l_step, l_contrib, l_void)
          d_ray_min = min(d_ray_min, l_step)
       enddo

       if (d_wall > d_ray_min * (1.0_dp + 1.0e-3_dp)) then
          write(msg,'(a,i0,a,es12.5,a,es12.5)') "cell ", icell, " d_wall=", d_wall, " > d_ray=", d_ray_min
          ok_upper = .false.
          exit
       endif
       if (d_wall < d_ray_min * 0.96_dp) then
          write(msg,'(a,i0,a,es12.5,a,es12.5)') "cell ", icell, " d_wall=", d_wall, " << d_ray=", d_ray_min
          ok_tight = .false.
          exit
       endif

       ! Point slightly outside a face must return 0 (regression for s_tmp < 0 -> huge bug)
       j = neighbours_list(Voronoi(icell)%first_neighbour)
       n_face = Voronoi(j)%xyz - Voronoi(icell)%xyz
       p_out = 0.5_dp * (Voronoi(j)%xyz + Voronoi(icell)%xyz) + 1.0e-4_dp * n_face
       if (distance_to_closest_wall_Voronoi(icell, p_out(1), p_out(2), p_out(3)) /= 0.0_dp) then
          ok_upper = .false.
          msg = "point outside Voronoi face did not return 0"
          exit
       endif

       ! Star-neighbour bounding check
       Voronoi(icell)%is_star_neighbour = .true.
       star(1)%r = norm2([x, y, z]) - 0.25_dp * d_wall
       d_star_expected = 0.25_dp * d_wall
       if (abs(distance_to_closest_wall_Voronoi(icell, x, y, z) - d_star_expected) > 1.0e-12_dp * L_scale) then
          ok_upper = .false.
          msg = "is_star_neighbour distance not bounded by stellar surface"
          exit
       endif
       Voronoi(icell)%is_star_neighbour = .false.
       star(1)%r = 0.03_dp * L_scale

       ! was_cut cell must return 0
       Voronoi(icell)%was_cut = .true.
       if (distance_to_closest_wall_Voronoi(icell, x, y, z) /= 0.0_dp) then
          ok_upper = .false.
          msg = "was_cut cell did not return 0"
          exit
       endif
       Voronoi(icell)%was_cut = .false.

       n_interior = n_interior + 1
    enddo

    call check(trim(label) // " boundary-wall cells -> 0", ok_wall_zero)
    call check(trim(label) // " interior cells tested & sphere contained", &
         ok_contain .and. (n_interior >= 6))
    call check(trim(label) // " d_wall <= min ray exit & edge guards", ok_upper, trim(msg))
    call check(trim(label) // " d_wall tight against min ray distance (AU units)", ok_tight, trim(msg))

    call free_voronoi_grid()
  end subroutine test_voronoi_walls

  pure integer function closest_voronoi_generator(x, y, z, np) result(best_i)
    real(dp), intent(in) :: x, y, z
    integer, intent(in) :: np
    integer :: i
    real(dp) :: d2, d2_min

    d2_min = huge(1.0_dp)
    best_i = 1
    do i = 1, np
       d2 = (Voronoi(i)%xyz(1) - x)**2 + (Voronoi(i)%xyz(2) - y)**2 + (Voronoi(i)%xyz(3) - z)**2
       if (d2 < d2_min) then
          d2_min = d2
          best_i = i
       endif
    enddo
  end function closest_voronoi_generator

  subroutine free_voronoi_grid()
    if (allocated(Voronoi)) deallocate(Voronoi)
    if (allocated(Voronoi_xyz)) deallocate(Voronoi_xyz)
    if (allocated(Voronoi_neighbour_xyz)) deallocate(Voronoi_neighbour_xyz)
    if (allocated(neighbours_list)) deallocate(neighbours_list)
    if (allocated(wall)) deallocate(wall)
    if (allocated(star)) deallocate(star)
  end subroutine free_voronoi_grid

end program test_mrw
