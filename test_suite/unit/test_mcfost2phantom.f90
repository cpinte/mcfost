program test_mcfost2phantom
  ! Unit test: mcfost called as a library by phantom (mcfost2phantom), as during a
  ! live radiation-hydrodynamics run: one initialisation, then several temperature
  ! calculations on particle sets that change between calls (more particles, a new
  ! sink). Each call must give a dust temperature to the gas particles it keeps, and
  ! to the particles on top of another one (ghosts). The stellar parameters given by the
  ! sinks (Teff, luminosity) must be used at each call, and the mass-based ones
  ! (isochrones) when the sinks give none.
  ! The particles are a Keplerian disc around a 1 Msun sink, in phantom units. As in
  ! phantom, the maximum number of particles is given at initialisation and the sink
  ! array has the maximum number of sinks: mcfost allocates once.
  ! With the argument "overflow", max_particles is not given and the second call has
  ! more particles than the first: mcfost must stop with an error (test_unit.py).
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py, from
  ! test_suite/ (parameter file unit/mcfost2phantom.para). Exits with status 1 on failure.

  use mcfost_env, only : dp
  use mcfost2phantom, only : init_mcfost_phantom, run_mcfost_phantom
  use parameters, only : n_stars, star

  implicit none

  real(dp), parameter :: au_cm = 1.495978707e13_dp, Msun_g = 1.988409870698051e33_dp, G_cgs = 6.6743e-8_dp
  real(dp), parameter :: pi = 3.141592653589793_dp
  real(dp), parameter :: udist = au_cm, umass = Msun_g
  real(dp), parameter :: hfact = 1.2_dp
  integer, parameter :: nsinkproperties = 18, ntypes = 1, ndusttypes = 0
  integer, parameter :: maxp = 8000, maxptmass = 10

  real(dp) :: utime, SPH_limits(6)
  integer :: ierr, n_fail
  character(len=16) :: mode

  n_fail = 0
  utime = sqrt(udist**3/(G_cgs*umass))
  mode = ""
  if (command_argument_count() > 0) call get_command_argument(1, mode)

  if (mode == "overflow") then
     call init_mcfost_phantom("unit/mcfost2phantom.para", ndusttypes, .false., "", SPH_limits, ierr)
     call run_disc(5000, 1, "overflow call 1: 5000 particles")
     call run_disc(8000, 1, "overflow call 2: 8000 particles") ! must stop with an error
     write(*,*) "FAIL: no error with more particles than allocated"
     call exit(0) ! test_unit.py expects a non-zero exit status
  endif

  call init_mcfost_phantom("unit/mcfost2phantom.para", ndusttypes, .false., "", SPH_limits, ierr, &
       max_particles=maxp)
  call check(ierr == 0, "init_mcfost_phantom")
  if (ierr /= 0) call finish()

  ! first call, then more particles and a new sink, as when phantom injects
  ! particles or forms a planet. Calls 1 and 2: Teff and luminosity of the star given by
  ! phantom; call 3: stellar parameters from the sink mass (isochrones)
  call run_disc(5000, 1, "call 1: 5000 particles, 1 sink", Teff=4000.0_dp, ghosts=.true.)
  call check(n_stars == 1 .and. abs(star(1)%T - 4000.) < 1., "call 1: Teff from phantom")
  call run_disc(8000, 2, "call 2: 8000 particles, 2 sinks", Teff=5000.0_dp, ghosts=.true.)
  call check(n_stars == 2, "call 2: new sink used")
  call check(abs(star(1)%T - 5000.) < 1., "call 2: Teff from phantom updated")
  call run_disc(3000, 1, "call 3: 3000 particles, 1 sink")
  call check(n_stars == 1 .and. abs(star(1)%T - 5000.) > 1., "call 3: Teff from the sink mass")
  write(*,*) "Teff of the star in call 3 (isochrones):", star(1)%T

  call finish()

contains

  subroutine run_disc(np, nptmass, label, Teff, ghosts)
    integer, intent(in) :: np, nptmass
    character(len=*), intent(in) :: label
    real(dp), intent(in), optional :: Teff ! of the central sink, with a luminosity of 1 Lsun
    logical, intent(in), optional :: ghosts ! particles 1-10 accreted (h < 0), 12 on 11, 14 and 15 on 13

    real(dp), parameter :: Lsun_W = 3.828e26_dp

    real(dp), parameter :: Mdisc = 1.0e-2_dp, rin = 10.0_dp, rout = 100.0_dp, h_r = 0.05_dp
    real(dp), allocatable :: xyzh(:,:), vxyzu(:,:), xyzmh_ptmass(:,:), vxyz_ptmass(:,:), dustfrac(:,:), T_gas(:)
    integer(kind=1), allocatable :: iphase(:), apr_level(:)
    real, allocatable :: Tphantom(:), n_packets(:)
    real(dp) :: grainsize(0), graindens(0), dudt(0), massoftype(ntypes), mu_gas
    integer :: npoftype(ntypes), i, n_T, ierr_run
    real(dp) :: r, phi, z, H, sigma, rho, m, vk, u(3)

    allocate(xyzh(4,np), vxyzu(4,np), xyzmh_ptmass(nsinkproperties,maxptmass), vxyz_ptmass(3,maxptmass), &
         dustfrac(ndusttypes,np), T_gas(np), iphase(np), apr_level(np), Tphantom(np), n_packets(np))

    ! surface density ~ 1/r between rin and rout (uniform in r), Gaussian vertical profile
    m = Mdisc/np
    do i=1, np
       call random3(i, np, u)
       r = rin + (rout - rin)*u(1)
       phi = 2*pi*u(2)
       H = h_r*r
       z = H*sqrt(-2*log(max(u(3),1.0e-6_dp)))*cos(2*pi*u(1)*7)
       sigma = Mdisc/(2*pi*(rout - rin)*r)
       rho = sigma/(sqrt(2*pi)*H)*exp(-0.5_dp*(z/H)**2)
       vk = sqrt(1.0_dp/r) ! G = M = 1 in code units
       xyzh(:,i) = [r*cos(phi), r*sin(phi), z, hfact*(m/rho)**(1.0_dp/3)]
       vxyzu(:,i) = [-vk*sin(phi), vk*cos(phi), 0.0_dp, 1.0e-4_dp]
    enddo
    if (present(ghosts)) then
       if (ghosts) then
          ! accreted particles (h < 0) are skipped by mcfost: the mcfost and phantom
          ! indices of the following particles differ
          xyzh(4,1:10) = -xyzh(4,1:10)
          xyzh(:,12) = xyzh(:,11) ; xyzh(:,14) = xyzh(:,13) ; xyzh(:,15) = xyzh(:,13)
       endif
    endif
    iphase = 1_1 ; apr_level = 1_1 ; T_gas = 20.0_dp
    npoftype = np ; massoftype = m

    ! sinks: the star at the centre, then a planet at 50 au. When the star has a Teff from
    ! phantom, the planet is a 100 K black body with its accretion radius (0.1 au) as
    ! radius: a large star in the Voronoi mesh (star neighbour flags, star fields not set
    ! by the reader such as fUV)
    xyzmh_ptmass = 0.0_dp ; vxyz_ptmass = 0.0_dp
    xyzmh_ptmass(4,1) = 1.0_dp ; xyzmh_ptmass(5:6,1) = 0.5_dp
    if (present(Teff)) then ! luminosity in code units (erg/s -> W), Teff in K
       xyzmh_ptmass(12,1) = Lsun_W / (umass*udist**2/utime**3 * 1.0e-7_dp)
       xyzmh_ptmass(13,1) = Teff
    endif
    if (nptmass > 1) then
       xyzmh_ptmass(1,2) = 50.0_dp ; xyzmh_ptmass(4,2) = 1.0e-3_dp ; xyzmh_ptmass(5:6,2) = 0.1_dp
       vxyz_ptmass(2,2) = sqrt(1.0_dp/50.0_dp)
    endif

    call run_mcfost_phantom(np, nptmass, ntypes, ndusttypes, 1, npoftype, xyzh, vxyzu, iphase, grainsize, &
         graindens, dustfrac, massoftype, xyzmh_ptmass, vxyz_ptmass, hfact, umass, utime, udist, 0, dudt, &
         .false., SPH_limits, Tphantom, n_packets, mu_gas, ierr_run, write_T_files=.false., ISM=0, &
         T_gas=T_gas, apr_level=apr_level, use_apr=.false.)

    ! every gas particle kept by mcfost (all but a few percent at the outer edge) has a
    ! dust temperature; the others are flagged with -1. Particles at the 1 K floor are
    ! midplane cells reached by no packet (no diffusion approximation in the library,
    ! 1e5 packets here): a few percent, all of them if the packets carry no energy
    n_T = count(Tphantom > 0.)
    call check(ierr_run == 0, label//": ierr")
    call check(n_T > 0.95*np, label//": particles with a temperature")
    call check(count(Tphantom > 0. .and. Tphantom <= 1.) < 0.25*np, label//": particles heated")
    call check(all(Tphantom > 0. .or. Tphantom == -1.), label//": temperature or -1 flag")
    call check(maxval(Tphantom) < 2000. .and. minval(Tphantom, mask=Tphantom > 0.) >= 1., &
         label//": temperature range")
    write(*,*) label, ": ", n_T, " particles with T in", minval(Tphantom, mask=Tphantom > 0.), maxval(Tphantom), &
         ", at the 1 K floor:", count(Tphantom > 0. .and. Tphantom <= 1.)

    ! ghosts get the temperature of the particle they are on top of
    if (present(ghosts)) then
       if (ghosts) then
          call check(all(Tphantom(1:10) == -1.), label//": accreted particles flagged")
          call check(Tphantom(11) > 0. .and. Tphantom(12) == Tphantom(11), label//": ghost of a pair")
          call check(Tphantom(13) > 0. .and. Tphantom(14) == Tphantom(13) .and. Tphantom(15) == Tphantom(13), &
               label//": ghosts of a triplet")
          write(*,*) label, ": T of particles 9-15:", Tphantom(9:15)
       endif
    endif

  end subroutine run_disc

  subroutine random3(i, n, u)
    ! deterministic pseudo-random numbers in ]0,1[ (no dependence on the RNG of the compiler)
    integer, intent(in) :: i, n
    real(dp), intent(out) :: u(3)
    real(dp), parameter :: alpha(3) = [0.6180339887498949_dp, 0.7548776662466927_dp, 0.5698402909980532_dp]
    integer :: k
    do k=1,3
       u(k) = modulo((i + 0.5_dp)*alpha(k) + k*real(n,dp)*1.0e-7_dp, 1.0_dp)
    enddo
  end subroutine random3

  subroutine check(ok, label)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: label
    if (.not.ok) then
       write(*,*) "FAIL: ", label
       n_fail = n_fail + 1
    endif
  end subroutine check

  subroutine finish()
    if (n_fail > 0) then
       write(*,*) n_fail, " failure(s)"
       call exit(1)
    endif
    write(*,*) "test_mcfost2phantom: all checks passed"
    call exit(0)
  end subroutine finish

end program test_mcfost2phantom
