module stars

  use parameters
  use utils
  use constants
  use messages
  use wavelengths
  use grid

  implicit none

  public :: star_spectrum, E_stars, ProDiMo_star_HR, R_ISM, E_ISM, prob_E_star

  public :: allocate_stellar_spectra, deallocate_stellar_spectra, emit_packet_uniform_sphere, emit_packet_ism, &
       ism_energy_distribution, star_energy_distribution, select_star, stars_cell_indices, find_spectra, &
       intersect_stars, distance_to_star, compute_stellar_parameters
  !-> to move in parameters ?
  public :: star_rad, is_inshock, laccretion_shock, max_Tshock, min_Tshock, min_Thp, max_Thp, T_hp, max_Facc, min_Facc, T_preshock

  private

  real, dimension(:,:), allocatable :: CDF_E_star, prob_E_star
  real, dimension(:), allocatable :: E_stars !n_lambda
  real, dimension(:), allocatable :: star_spectrum_cumul, star_spectrum !(0:n_lambda)

  real, dimension(:), allocatable :: E_ISM
  real(kind=dp) :: R_ISM = 0._dp ! radius of the sphere from which the ISM radiation is emitted
  real(kind=dp), dimension(3) :: centre_ISM  ! centre of the ISM emitting sphere

  real, dimension(:,:), allocatable :: ProDiMo_star_HR

  !onto the star(s)
  real(kind=dp) :: T_hp, max_Thp = 0.0, min_Thp = 1d8 !photosphere heated.
  real(kind=dp) :: max_Tshock = 0.0, min_Tshock = 1d8 !soft X-rays emission from the shock
  real(kind=dp) :: max_Facc = 0.0, min_Facc = 1d8 !Accretion flux in W/m2
  real(kind=dp) :: T_preshock !temperature of the opt-thin pre-shock region

  contains

subroutine allocate_stellar_spectra(n_wl)

  integer, intent(in) :: n_wl
  integer :: alloc_status

  allocate(CDF_E_star(n_wl,0:n_stars), prob_E_star(n_wl,n_stars), E_stars(n_wl),  &
       E_ISM(n_wl), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error CDF_E_star')
  CDF_E_star = 0.0
  prob_E_star = 0.0
  E_stars = 0.0
  E_ISM = 0.0

  allocate(star_spectrum_cumul(0:n_wl),star_spectrum(n_wl), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error spectre_etoile')
  star_spectrum_cumul = 0.0
  star_spectrum = 0.0

  return

end subroutine allocate_stellar_spectra

!**********************************************************************

subroutine deallocate_stellar_spectra()

  if (allocated(star_spectrum)) deallocate(star_spectrum,star_spectrum_cumul)
  if (allocated(CDF_E_star)) deallocate(CDF_E_star,prob_E_star,E_stars,E_ISM)

  return

end subroutine deallocate_stellar_spectra

!**********************************************************************

subroutine select_star(lambda,rand,n_star)
! Selection of the star that will emit the photon
! C. Pinte
! 21/05/05

  implicit none

  integer, intent(in) :: lambda
  real, intent(in) :: rand
  integer, intent(out) :: n_star
  integer :: k, kmin, kmax

  ! Dichotomy (Binary search)
  kmin=0
  kmax=n_stars
  k=(kmax-kmin)/2

  do while ((kmax-kmin) > 1)
     if (CDF_E_star(lambda,k) < rand) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   n_star=kmax

   return

end subroutine select_star

!**********************************************************************

subroutine emit_packet_uniform_sphere(id, i_star,rand1,rand2,rand3,aleat4, icell,x,y,z,u,v,w,w2,lintersect)
! Chooses the emission position uniformly
! on the surface of the star and the flight direction
! according to the cosine of the angle with the normal
! C. Pinte
! 21/05/05

  implicit none

  integer, intent(in) :: id, i_star
  real, intent(in) :: rand1, rand2, rand3, aleat4
  integer, intent(out) :: icell

  real(kind=dp), intent(out) :: x, y, z, u, v, w, w2

  real(kind=dp) :: srw02, argmt, r_star, cospsi, phi
  real(kind=dp), parameter :: precision = 1e-6_dp
  logical, intent(out) :: lintersect

  ! Random starting position on a sphere of radius 1
  z = 2.0_dp * rand1 - 1.0_dp
  srw02 = sqrt(1.0-z*z)
  argmt = pi*(2.0_dp*rand2-1.0_dp)
  x = srw02 * cos(argmt)
  y = srw02 * sin(argmt)

  ! Choice of flight direction: uniform sphere
  cospsi = sqrt(rand3) ! !TODO : This is strange. rand3 works with the benchmark: but not the same thing for the stellar surface in RT. sqrt(rand3) ~ OK with TORUS in MC mode
  phi = 2.0_dp*pi*aleat4
  ! (x,y,z) defines the normal (here, it still has norm=1)
  ! cospsi and phi are defined relative to this normal
  ! a rotation is required
  call cdapres(cospsi, phi, x, y, z, u, v, w)

  w2=1.0_dp-w*w

  ! Random starting position on a sphere of radius r_star
  r_star = star(i_star)%r * (1._dp + precision)
  x = x * r_star
  y = y * r_star
  z = z * r_star

  ! Adding position of the star
  x=x+star(i_star)%x
  y=y+star(i_star)%y
  z=z+star(i_star)%z

  if (lVoronoi) then
     icell = star(i_star)%icell
  else ! star can overlap several cells on a classical grid
     call index_cell(x,y,z, icell)
  endif

  if (star(i_star)%out_model) then
     call move_to_grid(id, x,y,z,u,v,w, icell,lintersect)
  else
     lintersect = .true.
  endif

  return

end subroutine emit_packet_uniform_sphere

!**********************************************************************

!subroutine em_star_ponctuelle(n_star,rand1,rand2,ri,zj,phik,x,y,z,u,v,w,w2)
!! Isotropic emission
!! C. Pinte
!! 21/05/05
!
!  implicit none
!
!  integer, intent(in) :: n_star
!  real, intent(in) :: rand1, rand2
!  integer, intent(out) :: ri, zj, phik
!  real(kind=dp), intent(out) :: x, y, z, u, v, w, w2
!
!  real(kind=dp) :: srw02, argmt
!
!  ! Isotropic emission
!  w = 2.0_dp * rand1 - 1.0_dp
!  w2 = 1.0_dp-w*w
!  srw02 = sqrt(w2)
!  argmt = pi*(2.0_dp*rand2-1.0_dp)
!  u = srw02 * cos(argmt)
!  v = srw02 * sin(argmt)
!
!  ! Star position
!  x=star(n_star)%x
!  y=star(n_star)%y
!  z=star(n_star)%z
!
!  ri=star(n_star)%ri
!  zj=star(n_star)%zj
!  phik=star(n_star)%phik
!
!  return
!
!end subroutine em_star_ponctuelle

!**********************************************************************

subroutine wien_law(lambda)

  implicit none

  integer, intent(out) :: lambda
  real :: T
  real :: wl

  T = maxval(star%T)

  wl = 2898./T

  loop : do lambda=1, n_lambda
     if (tab_lambda(lambda) > wl) exit loop
  enddo loop
  lambda = lambda - 1

end subroutine wien_law

!***********************************************************

subroutine star_energy_distribution()
! Calculates the total energy emitted by the stars
! and the emission probability of each star for all lambdas
! For all code versions
! C. Pinte
! 21/05/05

! Modif : 26/07/05 (C. Pinte)
! Merged with create_b_body: the routine also creates the star spectrum,
! the probability of emitting at lambda for all stars
! (Redundant calculations)
!
! returns:
! - star_spectrum
! - CDF_E_star, prob_E_star
! - E_stars
! - star_spectrum_cumul, star_spectrum

  implicit none

  real, dimension(n_lambda,n_stars) :: prob_E_star0
  real(kind=dp), dimension(n_lambda) :: log_lambda, star_spectrum0
  real(kind=dp), dimension(n_stars) ::  L_star0, correct_Luminosity
  real(kind=dp), dimension(n_stars) :: Lacc, Tacc

  real, dimension(:,:), allocatable :: spectre_tmp, tab_lambda_spectre, tab_spectre, tab_spectre0, tab_bb
  real(kind=dp), dimension(:), allocatable :: log_spectre, log_spectre0, log_wl_spectre
  character(len=512) :: filename, dir

  integer, dimension(n_stars) :: n_lambda_spectre, unit
  integer :: lambda, i, n, l, ios, n_lambda_spectre_max, n_wl
  real(kind=dp) :: wl, cst_wl, delta_wl, surface, terme, terme0, spectrum, spectre0, Cst0
  real ::  wl_inf, wl_sup, UV_ProDiMo, p, cst_UV_ProDiMo, correct_UV, fUV
  real(kind=dp) :: fact_sup, fact_inf, cst_spectre_etoiles

  real(kind=dp) :: wl_spectre_max, wl_spectre_min, wl_spectre_avg, wl_deviation

  integer :: status, readwrite, blocksize,nfound,group,firstpix,nbuffer,npixels
  real :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull

  if (n_stars < 1) return

  ! to obtain a spectrum normalized to 1 Rsun and 1 pc
  ! Cst0 is used to renormalize the black body to align it with the spectra
  Cst0 =  2.0*hp*c_light**2 * 1e-6
  surface= pi*Rsun_to_AU**2 ! for 1 solar radius
  Cst0 = Cst0 * surface / (pc_to_AU**2)

  ! for the MCFOST spectrum
  cst_spectre_etoiles = 2.0*pi*hp*c_light**2 * (AU_to_m)**2 ! BB constant + r_star in AU

  if (letape_th) then
     fact_sup = exp(0.5_dp/real(n_lambda,kind=dp)*log(lambda_max/lambda_min))
     fact_inf = 1.0_dp/fact_sup
  else
     fact_sup = 1.0001_dp ;
     fact_inf = 0.9999_dp ;
  endif

  ! Star emission
  ! The black body is defined here up to a constant,
  ! because we look at it relative to what is done in the energy distribution for the Disk.
  ! In particular, we removed a pi everywhere.

  do i=2, n_stars
     if (star(i)%lb_body .neqv. (star(1)%lb_body)) then
        call error("all stars must be black bodies or", &
             msg2="all stars must not be black bodies : you cannot mix")
     endif
  enddo
  if (star(1)%lb_body) star(:)%find_spectrum = .false.

  call find_spectra()

  if (star(1)%lb_body) then ! the stars are black bodies
     ! Creation of a high-resolution body in F_lambda
     ! R = 1 Rsun and distance = 1 pc
     n_lambda_spectre(:) = 1000

     do i=1, n_stars

        if (i==1) then
           allocate(tab_lambda_spectre(n_stars,n_lambda_spectre(1)), &
                tab_spectre(n_stars,n_lambda_spectre(1)), tab_spectre0(n_stars,n_lambda_spectre(1)), &
                tab_bb(n_stars,n_lambda_spectre(1)))
           tab_lambda_spectre = 0.0 ; tab_spectre = 0.0 ;  tab_spectre0 = 0.0 ;  tab_bb = 0.0
           allocate(log_spectre(n_lambda_spectre(1)), log_spectre0(n_lambda_spectre(1)), log_wl_spectre(n_lambda_spectre(1)))
        endif
        tab_lambda_spectre(i,:) = 1.0_dp * spanl(lambda_min, lambda_max, n_lambda_spectre(1))

        do l=1, n_lambda_spectre(1)
           wl = tab_lambda_spectre(i,l) *1.e-6
           cst_wl=thermal_const/(star(i)%T*wl)
           tab_spectre(i,l) = max(Cst0/ ( ((exp(min(cst_wl,700.)) -1.)+1.e-30) * (wl**5)), 1e-200_dp) ;
        enddo ! l

     enddo !i

  else ! the stars are not black bodies
     ! We calculate two things at the same time:
     ! - CDF_E_star: cumulative probability at fixed lambda for emission as a function of the star
     ! - star_spectrum: probability of emitting at lambda for all stars
     ! Reading the spectra
     n_lambda_spectre_max = 0
     do i=1, n_stars
        ! --- Reading the stellar spectrum
        ! --- Fluxes are normalized for R=1 solar radius seen from 1 pc
        ! --- units: F_lambda: W.m-2 / micron
        !             lambda: micron
        !     ---> lambda x F_lambda: W.m-2

        filename=trim(star(i)%spectrum)
        dir = in_dir(filename, star_dir,  status=ios)
        if (ios /=0) then
           call error("star file cannot be found:",trim(filename))
        else
           filename = trim(dir)//trim(filename) ;
           write(*,*) "Reading "//trim(filename) ;
        endif

        status=0
        !  Get an unused Logical Unit Number to use to open the FITS file.
        call ftgiou(unit(i),status)

        readwrite=0
        call ftopen(unit(i),filename,readwrite,blocksize,status)
        if (status /= 0) call error("cannot open fits file "//trim(filename))

        !  determine the size of the image
        call ftgknj(unit(i),'NAXIS',1,2,naxes,nfound,status)
        if (status /= 0) call error("reading axes of "//trim(filename))
        !  check that it found both NAXIS1 and NAXIS2 keywords
        if (nfound /= 2) call error("failed to read the NAXISn keywords in "//trim(filename))

        if (naxes(2) /= 3) then
           write(*,*) "NAXIS2 =", naxes(2)
           call error(trim(filename)//" does not have the right shape")
        endif

        ! We first read the length of the spectrum
        n_lambda_spectre(i) = naxes(1)
        if (naxes(1) > n_lambda_spectre_max) n_lambda_spectre_max = naxes(1)
     enddo

     ! We allocate the array to the maximum length
     allocate(tab_lambda_spectre(n_stars,n_lambda_spectre_max), &
          tab_spectre(n_stars,n_lambda_spectre_max), tab_spectre0(n_stars,n_lambda_spectre_max), &
          tab_bb(n_stars,n_lambda_spectre_max))
     tab_lambda_spectre = 0.0 ; tab_spectre = 0.0 ;  tab_spectre0 = 0.0 ;  tab_bb = 0.0

     allocate(log_spectre(n_lambda_spectre_max), log_spectre0(n_lambda_spectre_max), log_wl_spectre(n_lambda_spectre_max))


     do i=1, n_stars
        n_wl = n_lambda_spectre(i)
        ! We read again the dimensions
        call ftgknj(unit(i),'NAXIS',1,10,naxes,nfound,status)

        !  initialize variables.
        npixels=n_wl * 3
        group=1
        firstpix=1
        nullval=-999
        nbuffer=npixels

        ! read_image
        allocate(spectre_tmp(n_wl,3))
        call ftgpve(unit(i),group,firstpix,nbuffer,nullval,spectre_tmp,anynull,status)

        tab_lambda_spectre(i,1:n_wl) = spectre_tmp(:,1)
        tab_spectre(i,1:n_wl) = spectre_tmp(:,2)
        tab_bb(i,1:n_wl) = spectre_tmp(:,3)

        call ftclos(unit(i), status)
        call ftfiou(unit(i), status)
        deallocate(spectre_tmp)
     enddo ! n_stars

  endif ! bb

  ! Star luminosity integrated over the spectrum
  L_star0(:) = 0.0
  do i=1, n_stars
     do l = 2, n_lambda_spectre(i)
        L_star0(i) = L_star0(i) + 0.5 * (tab_spectre(i,l) + tab_spectre(i,l-1)) &
             * (tab_lambda_spectre(i,l) - tab_lambda_spectre(i,l-1))
     enddo
  enddo
  correct_Luminosity(:) = (sigma*(star(:)%T)**4 * (Rsun_to_AU/pc_to_AU)**2) / L_star0(:)

  ! Normalization of the spectrum to the luminosity indicated in the parameter file
  do i=1, n_stars
     tab_spectre(i,:) = tab_spectre(i,:) * correct_Luminosity(i)
  enddo

  ! Saving stellar spectrum before adding UV
  tab_spectre0(:,:) = tab_spectre(:,:)


  !--------------------------------
  ! Additional UV flux for ProDiMo
  !--------------------------------
  wl_inf = 91.2e-9 ; ! in m
  wl_sup = 250e-9 ; ! en m

  ! wl must be in microns here
  do i=1, n_stars
     fUV = star(i)%fUV
     if (fUV > tiny_real) then
        p = star(i)%slope_UV

        if (abs(p+1.0) > 1.0e-5) then
           cst_UV_ProDiMo =  fUV * L_star0(i) * (p+1) / (wl_sup**(p+1) - wl_inf**(p+1)) / 1e6 !/ (1e6)**(p+1)
        else
           cst_UV_ProDiMo =  fUV * L_star0(i) * log(wl_sup/wl_inf) / 1e6 !/ (1e6)**(p+1)
        endif

        ! We add UV only before the Wien peak
        do l = 1, n_lambda_spectre(i)
           if (tab_lambda_spectre(i,l)  < 2898./star(i)%T ) then
              wl = tab_lambda_spectre(1,l) * 1e-6
              UV_ProDiMo =  cst_UV_ProDiMo * wl**p
              if (UV_ProDiMo >  tab_spectre(i,l)) tab_spectre(i,l) = UV_ProDiMo
           endif
        enddo
     endif
  enddo ! n_stars

  !--------------------------------
  ! Calculate accretion spectrum
  !--------------------------------
  if (.not.lturn_off_Lacc) then
     if (maxval(star(:)%Mdot) > tiny_real) then
        ! Luminosity from the accretion [au^2 W / m^2]
        Lacc(:) = Ggrav/AU3_to_m3 &                         ! convert G to AU^3 / s^s / kg
             * star(:)%M*Msun_to_kg &                     ! M in kg
             * star(:)%Mdot*Msun_to_kg/year_to_s &        ! Mdot in kg / s
             / star(:)%r                                  ! R in AU
        ! Converting Lacc to Tacc
        Tacc(:) = (Lacc(:)/(four_pi * sigma * star(:)%r**2))**0.25

        write(*,*) "Accretion onto stars: "
        write(*,*) "Mdot=", star(:)%Mdot, "Msun/yr"
        write(*,*) "Tacc=", real(Tacc(:)), "K"

        ! We add a black-body to the stellar spectrum
        do i=1, n_stars
           if (Tacc(i) > tiny_real) then
              do l=1, n_lambda_spectre(i)
                 wl = tab_lambda_spectre(i,l) * 1.e-6
                 cst_wl=thermal_const/(Tacc(i)*wl)
                 tab_spectre(i,l) = tab_spectre(i,l) +  max(Cst0/ ( ((exp(min(cst_wl,700.)) -1.)+1.e-30) * (wl**5)), 1e-200_dp) ;
              enddo ! l
           endif
        enddo
     endif
  else
     write(*,*) "Turning off accretion luminosity"
  endif

  !---------------------------------------------------------------------------
  ! We calculate two things at the same time:
  ! - CDF_E_star: cumulative probability at fixed lambda for emission as a function of the star
  ! - star_spectrum: probability of emitting at lambda for all stars
  !---------------------------------------------------------------------------

  !---------------------------------------------------
  ! Integration of the spectrum in the MCFOST bands
  !---------------------------------------------------
  star_spectrum(:) = 0.0
  star_spectrum0(:) = 0.0

  log_lambda = log(tab_lambda(:))
  do i=1,n_stars
     n_wl = n_lambda_spectre(i)
     surface=4*pi*(star(i)%r**2)

     wl_spectre_max = maxval(tab_lambda_spectre(i,1:n_wl))
     wl_spectre_min = minval(tab_lambda_spectre(i,1:n_wl))

     log_spectre = log(tab_spectre(i,:) + 1e-30)
     log_spectre0 = log(tab_spectre0(i,:) + 1e-30)
     log_wl_spectre = log(tab_lambda_spectre(i,:))

     do lambda=1, n_lambda

        wl = tab_lambda(lambda)*1.e-6
        delta_wl=tab_delta_lambda(lambda)*1.e-6
        ! delta_wl is the width of the integration bin
        CDF_E_star(lambda,0) = 0.0

        wl_inf =  tab_lambda_inf(lambda)
        wl_sup =  tab_lambda_sup(lambda)

        ! Calculation of term by binning the input spectrum
        terme = 0.0 ; terme0 = 0.0 ; N = 0

        wl_spectre_avg = 0.0
        do l=1, n_wl
           if ( (tab_lambda_spectre(i,l) > wl_inf).and.(tab_lambda_spectre(i,l) < wl_sup) ) then
              terme = terme + tab_spectre(i,l)
              terme0 = terme0 + tab_spectre0(i,l)
              wl_spectre_avg = wl_spectre_avg + tab_lambda_spectre(i,l)
              N = N + 1
           endif
        enddo ! l

        if (N>1) wl_spectre_avg = wl_spectre_avg / N
        wl_deviation = wl_spectre_avg / tab_lambda(lambda) ! Deviation between bin center and averaged wl

        ! Potential corrections
        if ((terme > tiny_dp) .and. (N>3) .and. (abs(wl_deviation-1.0) < 0.1)) then
           terme = terme / N * (surface / Cst0)
           terme0 = terme0 / N * (surface / Cst0)
        else ! outside the provided spectrum
          if (tab_lambda(lambda) < wl_spectre_min) then
             cst_wl=thermal_const/(star(i)%T*wl)
             if (cst_wl < 500.) then
                terme =  surface/((wl**5)*(exp(cst_wl)-1.0))  ! BB for lack of a better option
             else
                terme = tiny_real
             endif
             terme0 = terme
          else if (tab_lambda(lambda) > wl_spectre_max) then ! extrapolation loi de puissance -2 en lambda.F_lambda
            ! A blackbody does not work as the level is higher due to line blanketing
            ! So a power-law extrapolation is used instead
            terme = (surface / Cst0) * tab_spectre(i,n_wl) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
            terme0 = (surface / Cst0) * tab_spectre0(i,n_wl) * (tab_lambda(lambda) / wl_spectre_max)**(-4)
          else
             !write(*,*) log_spectre
             terme = (surface / Cst0) * exp(interp(log_spectre(1:n_wl), log_wl_spectre(1:n_wl), log_lambda(lambda)))
             terme0 = (surface / Cst0) * exp(interp(log_spectre0(1:n_wl), log_wl_spectre(1:n_wl), log_lambda(lambda)))
          endif
        endif ! Fin correction

        ! No delta_wl here as we must compare with disk emission
        prob_E_star(lambda,i) = terme
        prob_E_star0(lambda,i) = terme0
      enddo ! lambda
  enddo ! etoiles


  ! We swap the loops to compute the summations over the stars at a given lambda
  do lambda=1, n_lambda
     delta_wl=tab_delta_lambda(lambda)*1.e-6

     spectrum = 0.0
     spectre0 = 0.0 ;
     CDF_E_star(lambda,0) = 0.0

     do i=1, n_stars
        terme = prob_E_star(lambda,i)
        terme0 = prob_E_star0(lambda,i)

        spectrum = spectrum + terme * delta_wl  ! total_sum over all stars
        spectre0 = spectre0 + terme0 * delta_wl

        ! No delta_wl here as we must compare with disk emission
        CDF_E_star(lambda,i) = CDF_E_star(lambda,i-1) +  terme
     enddo ! i, etoiles

     star_spectrum(lambda) = star_spectrum(lambda) + spectrum
     star_spectrum0(lambda) =  star_spectrum0(lambda) + spectre0

     ! Total star emission
     ! Correction by Teff in the parameter file
     E_stars(lambda) = CDF_E_star(lambda,n_stars)

     ! Normalization to 1 of the stars' emission probability
     if (CDF_E_star(lambda,n_stars) > 0.) then
        prob_E_star(lambda,:) = prob_E_star(lambda,:)/CDF_E_star(lambda,n_stars)
        CDF_E_star(lambda,:) = CDF_E_star(lambda,:)/CDF_E_star(lambda,n_stars)
     endif
  enddo ! lambda

  correct_UV = sum(star_spectrum) / sum(star_spectrum0)

  ! Multiplication by star radius and distance (in Rsun and pc)
  ! star_spectrum is F_lambda * dlambda
  star_spectrum(:) =  star_spectrum(:) * cst_spectre_etoiles

  if ( (lProDiMo).and.(.not.allocated(ProDiMo_star_HR)) ) then
     ! Only 1 star in ProDiMo mode
     ! ProDiMo_star_HR is lambda * F_lambda (same as star_spectrum but with tab_lambda at instead of tab_delta_lambda)
     allocate(ProDiMo_star_HR(n_lambda_spectre_max,2))
     ProDiMo_star_HR(: ,1) = tab_lambda_spectre(1,:)
     ProDiMo_star_HR(:,2) = tab_spectre(1,:) * (surface / Cst0) * cst_spectre_etoiles  * tab_lambda_spectre(1,:) * 1e-6
  endif

  ! Cumulative probability
  if (n_lambda > 1) then
     ! Normalization to 1 of the spectrum
     do lambda =1, n_lambda
        star_spectrum_cumul(lambda) = star_spectrum_cumul(lambda-1) + star_spectrum(lambda)
     enddo
     do lambda=1,n_lambda
        star_spectrum_cumul(lambda)=star_spectrum_cumul(lambda)/star_spectrum_cumul(n_lambda)
     enddo
  endif

  ! origin cell where the star is
  if (.not.lVoronoi) then ! already done during tessellation for Voronoi grid
     do i=1, n_stars
        call index_cell(star(i)%x,star(i)%y,star(i)%z, star(i)%icell)
     enddo
  endif

  return

end subroutine star_energy_distribution

!***********************************************************

subroutine ism_energy_distribution(ISM_model)

  integer, intent(in) :: ISM_model ! 0 : no ISM radiation field, 1 : ProDiMo, 2 : Bate & Keto

  real :: wl, nu, ev_to_Hz, E, nu_p_MIR, Jnu, E_ev
  integer :: lambda, k

  real, dimension(5) :: wavelengths, power, W, T

  ! Defining the ISM sphere
  if (lcylindrical) then
     R_ISM = 1.000001_dp * (sqrt(Rmax**2 + zmax(n_rad)**2))
     centre_ISM(:) = 0._dp
  else if (lspherical) then
     R_ISM = 1.000001_dp * Rmax
     centre_ISM(:) = 0._dp
  else if (lVoronoi) then
     ! Defining the ISM sphere
     R_ISM = 1.000001_dp * 0.5_dp * Rmax
     !centre_ISM(:) = 0.5_dp * (/limits(2)+limits(1), limits(4)+limits(3), limits(6)+limits(5)/)
     centre_ISM(:) = 0._dp
  endif

  eV_to_Hz = electron_charge/hp
  nu_p_MIR = c_light/100.e-6

  if (R_ISM < tiny_dp) call error("the ISM emitting sphere is not defined")

  if (ISM_model==0) then
     E_ISM(:) = 0.
  else if (ISM_model==1) then ! ProdiMo
     do lambda=1, n_lambda
        wl = tab_lambda(lambda) * 1e-6
        E_ISM(lambda) = (chi_ISM * 1.71 * Wdil * Blambda(wl,T_ISM_stars) + Blambda(wl,TCmb))
     enddo
  else if (ISM_model==2) then ! Bate & Keto
     ! Interstellar Radiation (ISR) field from Zucconi et al. (2001), following Black (1994)
     ! But adds Draine (1978) UV field (equation 11) which is similar to
     ! Black (1994) but doesn't seem to be included in Zucconi.
     ! Code adapted from sphNG routine
     wavelengths = (/0.4E-4,0.75E-4,1.0E-4,140.0E-4,1.06E-1/) * cm_to_m
     power = (/0.,0.,0.,1.65,0./)
     W = (/1e-14,1e-13,4e-13,2e-4,1./)
     T = (/7500.,4000.,3000.,23.3,2.728/)

     do lambda=1, n_lambda
        wl = tab_lambda(lambda) * 1e-6 ; nu = c_light/wl

        E = 0.
        ! Sum of black-bodies
        do k=1,5
           E = E + (wavelengths(k)/wl)**power(k) * W(k) * Blambda(wl,T(k))
        enddo

        ! Add mid-infrared which has a cut-off longer than 100 microns
        if (tab_lambda(lambda) < 100) then
           Jnu = 5.0E-7* (2.0*hp*nu_p_MIR**3/c_light**2) * (tab_lambda(lambda)/100.)**1.8

           E = E + c_light/wl**2 * Jnu ! Adding J_lambda
        endif

        ! Add Draine 1978 UV
        if ((nu > 5*ev_to_Hz).and.(nu < 13.6*eV_to_Hz)) then
           E_eV = nu/ev_to_Hz
           E = E + (1.658e6*E_eV - 2.152e5*E_eV**2 + 6.919e3*E_eV**3) * hp*E_eV
        endif

        E_ISM(lambda) = E
     enddo
  else
     call error("Unknown ISM model")
  endif

  ! Normalization for MCFOST
  E_ISM(:) = E_ISM(:) * (4.*R_ISM**2) * 2.0/(hp *c_light**2) * 0.4

  return

end subroutine ism_energy_distribution

!***********************************************************

subroutine emit_packet_ISM(id, icell,x,y,z,u,v,w,stokes,lintersect)
! Chooses the emission position uniformly
! on a sphere and the flight direction
! according to the cosine of the angle with the normal
! C. Pinte
! 27/05/09

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: id
  integer, intent(out) :: icell
  real(kind=dp), intent(out) :: x, y, z, u, v, w
  real(kind=dp), dimension(4), intent(out) :: stokes
  logical, intent(out) :: lintersect

  real :: rand1, rand2, rand3, aleat4
  real(kind=dp) :: srw02, argmt, cospsi, phi, l, w2

  ! Energy at 1
  stokes(:) = 0. ; stokes(1)  = 1.

  ! Random starting position on a sphere of radius 1
  rand1 = sprng(stream(id))
  rand2 = sprng(stream(id))

  z = 2.0_dp * rand1 - 1.0_dp
  srw02 = sqrt(1.0-z*z)
  argmt = pi*(2.0_dp*rand2-1.0_dp)
  x = srw02 * cos(argmt)
  y = srw02 * sin(argmt)

  ! Choice of flight direction: uniform sphere
  ! emission towards the interior
  rand3 = sprng(stream(id))
  aleat4 = sprng(stream(id))

  cospsi = -sqrt(rand3)
  phi = 2.0_dp*pi*aleat4
  ! (x,y,z) defines the normal (here, it still has norm=1)
  ! cospsi and phi are defined relative to this normal
  ! a rotation is required
  call cdapres(cospsi, phi, x, y, z, u, v, w)

  w2=1.0_dp-w*w

  ! Random starting position on a sphere of radius r_star
  l = R_ISM
  x = centre_ISM(1) + x * l
  y = centre_ISM(2) + y * l
  z = centre_ISM(3) + z * l

  call move_to_grid(id, x,y,z,u,v,w, icell,lintersect)

  return

end subroutine emit_packet_ISM

!***********************************************************

subroutine stars_cell_indices()

  real(kind=dp) :: x, y, z
  integer :: i_star, icell

  do i_star=1, n_stars
     x = star(i_star)%x
     y = star(i_star)%y
     z = star(i_star)%z

     ! todo: the star can occupy several cells: No,
     call index_cell(x,y,z, icell)
     star(i_star)%icell = icell

     star(i_star)%out_model = test_exit_grid(icell, x,y,z)
  enddo

  return

end subroutine stars_cell_indices

!***********************************************************

subroutine intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! lintersect_star is true is the ray/packet will intersect a star
  ! i_star returns the id of the star intersected
  ! icell_star returns the id of the cell where the star is

  ! This routine implies that a star is in a unique cell

  real(kind=dp), intent(in) :: x,y,z, u,v,w
  logical, intent(out) :: lintersect_stars
  integer, intent(out) :: i_star, icell_star

  real(kind=dp), dimension(3) :: r, k, delta_r
  real(kind=dp) :: b,c, delta, rac, s1, s2, d_to_star
  integer :: i

!   d_to_star = distance_to_stars(x,y,z,u,v,w,i_star)
!   lintersect_stars = (i_star > 0)
!   if (lintersect_stars) then
!      icell_star = star(i_star)%icell
!   else
!      icell_star = 0
!   end if
!
!   return


  r(1) = x ; r(2) = y ; r(3) = z
  k(1) = u ; k(2) = v ; k(3) = w

  d_to_star = huge(1.0_dp)

  i_star = 0
  star_loop : do i = 1, n_stars
     delta_r(:)  = r(:) - (/star(i)%x, star(i)%y, star(i)%z/)
     b = dot_product(delta_r,k)
     c = dot_product(delta_r,delta_r) - (star(i)%r)**2
     delta = b*b - c

     if (delta >= 0.) then ! the packet will encounter (or has encoutered) the star
        rac = sqrt(delta)
        s1 = -b - rac

        if (s1 < 0) then ! we already entered the star
           ! We can probably skip that test, s1 must be positive as we must be outside the star
           s2 = -b + rac
           if (s2 > 0) then ! for s2 < 0: we already exited the star
              ! We are still in the sphere and will exit it
              ! This means that we had a round-off error somewhere
              d_to_star = 0.0_dp
              i_star = i
           endif
        else ! We will enter in the star
           if (s1 < d_to_star) then
              d_to_star = s1
              i_star = i
           endif
        endif ! s1 < 0
     endif ! delta < 0
  enddo star_loop

  lintersect_stars = (i_star > 0)
  if (lintersect_stars) then
     icell_star = star(i_star)%icell
  else
     icell_star = 0
  end if

  return

end subroutine intersect_stars

!***********************************************************

   function star_rad(id,iray,i_star,icell0,x,y,z,u,v,w,N,lambda)
   ! ---------------------------------------------------------------!
   ! routine to manage the radiation of the stellar boundary
   !
   ! TO  DO:
   ! - add limb darkening
   ! - add reading spectrum direclty (either flux or full CLV)
   ! -------------------------------------------------------------- !
      integer, intent(in) :: i_star, icell0, id, iray, N
      real(kind=dp), intent(in) :: u, v, w, x, y, z, lambda(N)
      real(kind=dp) :: Tshock, Thp, Facc
      real(kind=dp) :: star_rad(N)
      real(kind=dp) :: Icorona

      if (star(i_star)%T <= 1e-6) then !even with spots
         star_rad(:) = 0.0_dp
         return !no radiation from the star
      endif

      !Stellar "surface" contribution
      star_rad(:) = Bpnu(N,lambda,star(i_star)%T*1d0)

      if (is_inshock(id, iray, i_star, icell0, x, y, z, Thp, Tshock, Facc)) then
         if (T_preshock > 0.0 ) then
            star_rad(:) = 0.0_dp !set to zero before accumulation
            where (lambda > 364.2096)
               star_rad(:) = star_rad(:) + Bpnu(N,lambda,Thp)
            elsewhere (lambda <= 364.2096) !-> can be very large between ~ 40 and 91 nm.
            ! elsewhere (lambda > 91.176 .and. lambda <= 364.2096)
               star_rad(:) = star_rad(:) + Bpnu(N,lambda,T_preshock)
            endwhere
         else
            star_rad(:) = Bpnu(N,lambda,Thp) !no contribution from the star
         endif
         return
      endif

      ! add coronal illumination at the "unresolved" stellar surface
      ! At the moment, the EUV radiation is constant across the wl range.
      ! Icorona in W/m2/Hz/sr. It is normalised such that
      ! int (dOmega int( dnu Icorona ) ) = F_EUV in W/m2
      ! **** At the moment, F_EUV is stored in star%fuv (W/m2) **** !
      if ( star(i_star)%fuv > 0.0) then
         Icorona = star(i_star)%fuv * 1d-9 * 911.76 / (pi*c_light*81.176)
         where (lambda >= 10 .and. lambda <= 91.176)
         ! Icorona = star(i_star)%fuv * 1d-9 * 504.0 / (pi*c_light*40.4)
         ! where (lambda >= 10 .and. lambda <= 50.4)
            star_rad(:) = star_rad(:) + Icorona
         endwhere
      endif

      return
   end function star_rad

  function is_inshock(id, iray, i_star, icell0, x, y, z, Thp, Tshock, Facc)
  !
  ! Computes the temperature of the heated photosphere and of the
  ! pre-shock region that will radiate away respectively 3/4 and 1/4 of the
  ! black body radiation at Thp and Tshock respectively (with Tshock >> Thp).
  !
  ! Thp ~ (Facc/sigma)**0.25 = (0.5 * rho * vs**3 / sigma)**3
  ! Ts ~ 3/16 * mu * amu / kb * vs**2
  !
   use grid, only : voronoi
   use constants, only : sigma, kb
   use elements_type, only : wght_per_H
   ! use density, only : gas_density
   logical :: is_inshock
   integer :: i_star, icell0, id, iray
   real(kind=dp), intent(out) :: Thp, Tshock, Facc
   real(kind=dp) :: x, y, z, rho
   real(kind=dp) :: Tloc, vaccr, vmod2, rr, sign_z

   is_inshock = .false.
   if (.not.laccretion_shock) return

   if (icell0<=n_cells) then
   !TO DO: gas_density(icell0) instead of nHtot
      rho = nHtot(icell0) * wght_per_H
      if (rho > 0.0) then ! even if icompute_atomRT(icell0) /= 0
         rr = sqrt( x*x + y*y + z*z)
         ! Get vaccr : the accretion velocity above the shock.
         if (lvoronoi) then !always 3d
            vaccr = Voronoi(icell0)%vxyz(1)*x/rr + Voronoi(icell0)%vxyz(2)*y/rr + Voronoi(icell0)%vxyz(3) * z/rr
            vmod2 = sum( Voronoi(icell0)%vxyz(:)**2 )
         else
         	if (vfield_coord==1) then
               if (l3D) then !needed here if not 2.5d
                  sign_z = 1.0_dp
               else
                  sign_z = sign(1.0_dp, z)
               endif
         		vaccr = vfield3d(icell0,1) * x/rr + vfield3d(icell0,2) * y/rr + vfield3d(icell0,3) * z/rr * sign_z
            elseif (vfield_coord==2) then
               if (l3D) then !needed here if not 2.5d
                  sign_z = 1.0_dp
               else
                  sign_z = sign(1.0_dp, z)
               endif
               vaccr = vfield3d(icell0,1) * sqrt(1.0 - (z/rr)**2) + sign_z * vfield3d(icell0,3) * z/rr
            else !spherical vector here
               vaccr = vfield3d(icell0,1) !always negative for accretion
            endif
            vmod2 = sum(vfield3d(icell0,:)**2)
         endif


         if (vaccr < 0.0_dp) then
            !Facc = 1/2 rho vs^3
            Facc = 0.5 * (1d-3 * mH * rho) * abs(vaccr)**3
            Tloc = ( 0.75 * Facc / sigma )**0.25
            ! is_inshock = (Tloc > 0.5 * star(i_star)%T)
            is_inshock = (T_hp > 1.0_dp * star(i_star)%T)
            Thp = T_hp
            if (T_hp<=0.0) then
               !depends on the local value
               is_inshock = (abs(T_hp) * Tloc > 1.0_dp*star(i_star)%T)
               Thp = abs(T_hp) * Tloc
            endif
            !assuming mu is 0.5
            Tshock = 0.5 * (3.0/16.0) * (1d-3 * mH) / kb * vaccr**2
            max_Thp = max(max_Thp, Thp); min_Thp = min(min_Thp, Thp)
            max_Tshock = max(max_Tshock, Tshock); min_Tshock = min(min_Tshock, Tshock)
            max_Facc = max(max_Facc,Facc); min_Facc = min(min_Facc, Facc)
         endif

      endif !icompute_atomRT
   endif !laccretion_shock

   return
  end function is_inshock

!***********************************************************

subroutine find_spectra()
  ! Find an appropriate spectrum for all star based on the Teff, mass and radius (i.e. log(g))

  real :: Teff, r, M, logg, min_logg, max_logg
  integer :: delta_T, i_star

  real, parameter :: delta_logg = 0.5

  character(len=32) :: sTeff, slogg, type

  write(*,*) "Trying to find appropriate stellar spectra ..."

  do i_star = 1, n_stars
     if (star(i_star)%find_spectrum) then
        Teff = star(i_star)%T

        r = star(i_star)%r / Rsun_to_AU
        M = star(i_star)%M

        if (M < tiny_real) then
           call warning("Stellar mass is not set, forcing log(g) = 3.5")
           logg = 3.5 ! stellar mass is not defined in the parameter file, we fix logg
        else
           logg = logg_Sun + log10(M/r**2)
        endif

        if (Teff < 100) then
           call warning("Teff below 100K needs to be implemented")
           Teff = 100
           type = "cond"
           min_logg = 2.5
           max_logg = 6
           delta_T = 100
        else if (Teff < 1500) then
           type = "cond"
           min_logg = 2.5
           max_logg = 6
           delta_T = 100
        else if (Teff < 2700) then
           type = "dusty"
           if ((Teff < 2250).and.(Teff > 2050)) then
              min_logg = 4.0 ! Some models appear to be missing
           else
              min_logg = 3.5
           endif
           max_logg = 6
           delta_T = 100
        else if (Teff < 10000) then
           type = "NextGen"
           min_logg = 3.5
           max_logg = 5.5
           if (Teff <= 4000) then
              delta_T = 100
           else
              delta_T = 200
           endif
        else if (Teff <= 35000) then
           type="Kurucz"

           if (Teff < 11000) then
              min_logg = 2.0
           else if (Teff < 19000) then
              min_logg = 2.5
           else if (Teff < 27000) then
              min_logg = 3.0
           else if (Teff < 32000) then
              min_logg = 3.5
           else
              min_logg = 4.0
           endif
           max_logg = 5.0

           if (Teff <= 12000) then
              delta_T = 500
           else
              delta_T = 1000
           endif

        else
           call error("Teff above 35000K needs to be implemented")
        endif

        ! Rounding off at the nearest point in the grid of stellar atmospheres
        Teff = nint(Teff/delta_T) * delta_T
        logg = nint(logg/delta_logg) * delta_logg
        logg = min(max(logg,min_logg), max_logg)

        if (Teff < 1000) then
           write(sTeff, "(I3)") int(Teff)
        else if (Teff < 10000) then
           write(sTeff, "(I4)") int(Teff)
        else
           write(sTeff, "(I5)") int(Teff)
        endif
        write(slogg, "(F3.1)") logg

        if (type=="Kurucz") then
           star(i_star)%spectrum = "Kurucz"//trim(sTeff)//"-"//trim(slogg)//".fits.gz"
        else
           star(i_star)%spectrum = "lte"//trim(sTeff)//"-"//trim(slogg)//"."//trim(type)//".fits.gz"
        endif

        write(*,*) "Star #", i_star, " --> ", trim(star(i_star)%spectrum)
     else ! We do not update the spectrum
        if (star(i_star)%lb_body) then
           write(*,*) "Star #", i_star, " --> BB at T=", star(i_star)%T, "K"
        else
           write(*,*) "Star #", i_star, " --> ", trim(star(i_star)%spectrum), " (forced)"
        endif
     endif
  enddo

  write(*,*) "Done"

  return

end subroutine find_spectra

  !*********************************************************

  subroutine compute_stellar_parameters()

    integer :: i

    character(len=512) :: isochrone_file, filename
    character(len=100) :: line_buffer
    character(len=1)   :: s_age

    character(len=2) :: SpT
    real :: L, R, T, M, maxM, logg, minM_Allard, maxM_Allard, minM_Siess, maxM_Siess
    real(kind=dp) :: Gcm_to_Rsun
    integer :: age, k

    logical :: lread_Siess, lread_Allard
    integer, parameter :: nSpT_Siess = 29
    real, dimension(nSpT_Siess) :: logR_Siess, logTeff_Siess, logM_Siess
    integer, parameter :: nSpT_Allard = 50
    real, dimension(nSpT_Allard) :: logR_Allard, logTeff_Allard, logM_Allard

    if (n_stars < 1) return

    minM_Siess = 0.1
    maxM_siess = 7.0

    minM_Allard = 0.0005
    maxM_Allard = minM_Siess

    ! Which models do we need to read ?
    lread_Siess = .False. ; lread_Allard = .False.
    do i=1, n_stars
       M = star(i)%M
       if (M > minM_Siess)  then
          lread_Siess = .True.
       else
          lread_Allard = .True.
       endif
    enddo

    if (lread_Siess) then
       ! Siess models
       isochrone_file = "Siess/isochrone_"//trim(system_age)//".txt"
       write(*,*) "Reading isochrone file: "//trim(isochrone_file)
       filename = trim(mcfost_utils)//"/Isochrones/"//trim(isochrone_file)

       open(unit=1,file=filename,status="old")
       do i=1,3
          read(1,*) line_buffer
       enddo
       minM_Siess = 1.e30 ; maxM = 0 ;
       do i=1, nSpT_Siess
          read(1,*) SpT, L, r, T, M
          logR_Siess(i) = log(r) ; logTeff_Siess(i) = log(T) ; logM_Siess(i) = log(M)
       enddo
       close(unit=1)
    endif

    if (lread_Allard) then
       ! Allard models if mass < 0.1Msun
       isochrone_file = "Allard/model.AMES-Cond-2000.M-0.0.2MASS.AB"
       write(*,*) "Reading isochrone file: "//trim(isochrone_file)
       filename = trim(mcfost_utils)//"/Isochrones/"//trim(isochrone_file)

       s_age = system_age(1:1)
       read(s_age,*) age ! age is an int with age in Myr

       open(unit=1,file=filename,status="old")
       Gcm_to_Rsun = 1e9 * cm_to_m/Rsun
       ! Skipping age block
       do k=1,age-1
          ! header
          do i=1,4
             read(1,*) line_buffer
          enddo
          ! data
          line_loop : do i=1,nSpT_Allard
             read(1,'(A)') line_buffer
             if (line_buffer(1:1) == "-") exit line_loop
          enddo line_loop
       enddo

       ! header
       do i=1,4
          read(1,*) line_buffer
       enddo
       ! data
       k=0
       line_loop2 : do i=1,nSpT_Allard
          read(1,'(A)') line_buffer
          if (line_buffer(1:1) == "-") exit line_loop2
          k = k+1
          read(line_buffer,*) M, T, L, logg, R
          logR_Allard(i) = log(r * Gcm_to_Rsun) ; logTeff_Allard(i) = log(T) ; logM_Allard(i) = log(M)
       enddo line_loop2
       close(unit=1)
    endif

    ! interpolate L and T, the functions are smoother
    write(*,*) ""
    write(*,*) "New stellar parameters:"
    do i=1, n_stars
       if (lturn_off_planets .and. i>1) then
          write(*,*) " "
          write(*,*) "*** WARNING : turning off emission fron sink particle"
          write(*,*) "*** object #", i, "M=", star(i)%M, "Msun"
          write(*,*) "*** The object will not radiate"
          star(i)%T = 3.
          star(i)%r = 1e-4
       else if (star(i)%M < minM_Allard) then
          write(*,*) " "
          write(*,*) "*** WARNING : stellar object mass is below isochrone range"
          write(*,*) "*** object #", i, "M=", star(i)%M, "Msun"
          write(*,*) "*** The object will not radiate"
          star(i)%T = 3.
          star(i)%r = 0.01
       else if (star(i)%M < maxM_Allard) then
          star(i)%T = exp(interp(logTeff_Allard(1:k), logM_Allard(1:k), log(star(i)%M)))
          star(i)%r = exp(interp(logR_Allard(1:k), logM_Allard(1:k), log(star(i)%M)))
       else ! using Siess' models
          if (star(i)%M > maxM_Siess) then
             write(*,*) " "
             write(*,*) "*** WARNING : stellar object mass is above in isochrone range"
             write(*,*) "*** object #", i, "M=", star(i)%M, "Msun"
             write(*,*) "*** Stellar properties are extrapolated"
          endif
          star(i)%T = exp(interp(logTeff_Siess, logM_Siess, log(star(i)%M)))
          star(i)%r = exp(interp(logR_Siess, logM_Siess, log(star(i)%M)))
       endif

       ! No fUV and no stellar spectrum for the moment
       star(i)%fUV = 0.0 ; star(i)%slope_UV = 0.0 ;
       star(i)%lb_body = .false. ! Not a bb by default
       star(i)%spectrum = "None"

       write(*,*) "Star #",i,"  Teff=", star(i)%T, "K, r=", real(star(i)%r), "Rsun"
    enddo

    ! Passage radius en AU
    star(:)%r = star(:)%r * Rsun_to_AU

    ! We force again a black-body if needed
    if (lstar_bb) star(:)%lb_body = .true.

    return

  end subroutine compute_stellar_parameters

end module stars
