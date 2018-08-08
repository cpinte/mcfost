module ProDiMo

  use parametres
  use opacity
  use em_th
  use constantes
  use molecular_emission
  use utils, only: get_NH, Blambda, Bnu
  use ray_tracing, only:  RT_imin, RT_imax, RT_n_incl, lRT_i_centered
  use radiation_field, only : xN_abs, xJ_abs
  use stars, only : spectre_etoiles, ProDiMo_star_HR, E_stars, R_ISM
  use read_params
  use sha
  use utils, only : appel_syst
  use messages

  implicit none

  save

  character(len=512), parameter :: mcfost2ProDiMo_file = "forProDiMo.fits.gz"
  integer :: mcfost2ProDiMo_version, ProDiMo2mcfost_version

  ! directory with *.in files and output dir
  character(len=512) :: ProDiMo_input_dir, data_ProDiMo

  ! Pour champ UV
  !real, parameter ::  slope_UV_ProDiMo = 2.2 - 2 ! Fnu -> F_lambda
  real :: fUV_ProDiMo, slope_UV_ProDiMo

  ! Variables a passer dans le Parameter.in
  real, dimension(:), allocatable :: ProDiMo_fPAH, ProDiMo_dust_gas, ProDiMo_Mdisk  ! n_regions
  integer :: ProDiMo_PAH_NC, ProDiMo_PAH_NH
  logical :: lPAH, ProDiMo_other_PAH
  character(len=10) :: sProDiMo_fPAH

  ! Grille de longeurs d'onde
  character(len=32) :: ProDiMo_tab_wavelength = "ProDiMo_UV3_9.lambda"

  real(kind=dp), dimension(:,:,:), allocatable :: J_ProDiMo, N_ProDiMo
  real(kind=dp), dimension(:,:), allocatable :: n_phot_envoyes_ISM

  ! ProDiMo2mcfost
  integer,parameter :: MC_NSP=5      ! number of species
  integer,parameter :: MC_LEVMAX=10  ! maximum number of levels
  integer, parameter :: NXX=1, NZZ=1 ! c'est pas dans le fichier de Peter

  type lpop_prodimo
     character(len=10) :: name(MC_NSP)       ! species name
     integer :: isp(MC_NSP)                  ! species indices
     integer :: lmax(MC_NSP)                 ! maximum level recorded
     real :: amass(MC_NSP)                   ! mass of atom/molecule
     real :: deltaV(MC_NSP,NXX,NZZ)          ! turbulent+thermal line widths
     real :: pops(MC_NSP,MC_LEVMAX,NXX,NZZ)  ! level populations
  end type lpop_prodimo

  type(lpop_prodimo) :: lpops

  ! Lecture du fichier for MCFOST.fits.gz
  integer, parameter :: nLevel_CII = 2
  integer, parameter :: nLevel_OI = 3
  integer, parameter :: nLevel_CO = 33  !52
  !integer, parameter :: nLevel_oH2O = 30 !10
  !integer, parameter :: nLevel_pH2O = 30 !9
  integer :: nLevel_oH2O, nLevel_pH2O

  real, dimension(:,:), allocatable :: nCII, dvCII, nOI, dvOI,  noH2O, dvoH2O, npH2O, dvpH2O, nCO, dvCO
  real, dimension(:,:,:), allocatable :: pop_CII
  real, dimension(:,:,:), allocatable :: pop_OI
  real, dimension(:,:,:), allocatable :: pop_CO
  real, dimension(:,:,:), allocatable :: pop_oH2O
  real, dimension(:,:,:), allocatable :: pop_pH2O


  type mcfost2ProDiMo_model

     ! parameters
     real*4  :: diskmass, Teff, Rstar, Mstar, fUV, pUV, rin, rout
     real*4  :: h0, rref, alpha, beta, amin, amax, aexp, settle, asettle

     ! Spatial grid (Unit: AU)
     real*8, allocatable,dimension(:,:) :: r_grid, z_grid

     ! dust material density
     real*4, allocatable,dimension(:) :: rho_grain

     ! dust temperature structure (Unit: K)
     real*4, allocatable,dimension(:,:) :: Tdust

     ! Wavelength bins (Unit: microns)
     real*4, allocatable,dimension(:) :: wavelengths

     ! Fraquency bins (Unit: Hz)
     real*4, allocatable,dimension(:) :: nu

     ! Stellar spectrum (lambda.I_lambda) (Units : W/m^2/sr)
     ! int I_lambda x dlambda = sigma/pi Teff^4 in case without extra UV excess
     real*4, allocatable,dimension(:) :: lamIlamStar

     ! ISM spectrum (lambda.F_lambda) (Units : W/m^2/sr)
     real*4, allocatable,dimension(:) :: lamIlamISM

     ! Radiation field
     ! lambda*J_lambda (Unit: W/m^2/sr)
     real*4, allocatable,dimension(:,:,:) :: lamJlam,lamJlamStar,lamJlamISM
     real*4, allocatable,dimension(:,:,:) :: Jnu !(Unit : W/m^2/sr/Hz)

     ! Statistics of the radiation field (number of packets)
     real*4, allocatable,dimension(:,:,:) :: nJ, nJ_Star, nJ_ISM

     ! Gas mass density (Unit: g cm^-3)
     real*4, allocatable,dimension(:,:) :: density

     ! 0th, 1st, 2nd and 3rd moment of the grain size distribution
     ! 0th moment: unit = part.m-3
     ! 1th, 2nd and 3rd moment: units = micron, micron^2, micron^3
     real*4, allocatable,dimension(:,:,:) :: grain_size

     ! Dust opacities (Units : AU^-1)
     real*4, allocatable,dimension(:,:,:) :: kappa_ext, kappa_abs

     ! local dust/gas ratio
     real*4, allocatable,dimension(:,:) :: dust_to_gas

     character(len=8) :: version
     integer*4 :: mcfost2prodimo
  end type mcfost2ProDiMo_model

  type(mcfost2ProDiMo_model) :: m2p


contains

  subroutine setup_ProDiMo()
    ! - distance 100pc
    ! - 1 pop de grains loi de puissance
    !
    ! Must be run after routines taille_grain & define_physical_zones()

    ! TODO : linit_gaz

    real :: eps_PAH, mPAH, masse_PAH, norme, a
    integer :: i, pop, k, NC, NH, alloc_status, ir, iz, NC_0
    real, dimension(:), allocatable :: fPAH

    ! Maximum 4 zones dans ProDiMo
    if (n_regions > 4) call error("ProDiMo cannot deal with more than 4 zones.")

    ! Directories
    if (.not.lprodimo_input_dir) then
       ProDiMo_input_dir = trim(mcfost_utils)//"/forProDiMo"
       write(*,*) "Using "//trim(ProDiMo_input_dir)//" for ProDiMo input"
    endif
    data_ProDiMo = trim(root_dir)//"/"//trim(seed_dir)//"/data_ProDiMo"

    if (lforce_ProDiMo_PAH) data_ProDiMo = trim(data_ProDiMO)//"_fPAH="//sProDiMo_fPAH

    ! Limite 10x plus haut pour avoir temperature plus propre pour ProDiMo
    !tau_dark_zone_eq_th = 15000.

    fUV_ProDiMo = etoile(1)%fUV
    slope_UV_ProDiMo = etoile(1)%slope_UV

    allocate(xN_abs(n_cells,n_lambda2,nb_proc),  stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error xN_abs')

    allocate(J_ProDiMo(n_lambda2,n_rad,nz),  stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error J_ProDiMo')
    J_ProDiMo = 0.0

    allocate(N_ProDiMo(n_lambda2,n_rad,nz),  stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error J_ProDiMo')
    N_ProDiMo = 0.0

    allocate(n_phot_envoyes_ISM(n_lambda2,nb_proc),  stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_phot_envoyes_ISM')
    n_phot_envoyes_ISM = 0.0


    ! Calcul des fPAH pour ProDiMo
    lPAH = .false.
    test_PAH : do i=1, n_pop
       if (dust_pop(i)%is_PAH) then
          lPAH = .true.
          exit test_PAH
       endif
    enddo test_PAH

    allocate(fPAH(n_zones))
    NC_0 = 0
    fPAH = 0.0 ;
    do i=1, n_zones
       do pop=1, n_pop
          if (dust_pop(pop)%zone == i) then
             if (dust_pop(pop)%is_PAH) then
                masse_PAH = 0.0 ;
                norme = 0.0 ;
                do  k=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                   a = r_grain(k)
                   NC = nint((a*1.e3)**3*468.)   ! number of Carbon atoms (Draine & Li 2001 Eq 8+)
                   NH = get_NH(NC)               ! number of Hydrogen atoms
                   mPAH = 12 * NC + NH    ! masse du PAH

                   if (NC_0 > 0) then
                      if (NC /= NC_0) call error("there can be only 1 type of PAH for ProDiMo")
                   else
                      NC_0 = NC
                   endif

                   masse_PAH = masse_PAH + mPAH * nbre_grains(k)
                   norme = norme + nbre_grains(k)
                enddo ! k
                masse_PAH = masse_PAH / norme  ! masse moyenne des PAHs en mH
                eps_PAH = dust_pop(pop)%frac_mass /disk_zone(i)%gas_to_dust  ! fraction en masse (de gaz) des PAHS

                ! 1.209274 = (mH2*nH2 + mHe * mHe) / (nH2 + nHe) avec nHe/nH2 = 10^-1.125
                ! abondance en nombre par rapport à H-nuclei + correction pour NC
                fPAH(i) = fPAH(i) + (1.209274/masse_PAH) * eps_PAH/3e-7 * (NC/50.)
                !write(*,*) i, fPAH(i), real(dust_pop(pop)%frac_mass * disk_zone(i)%diskmass), real(dust_pop(pop)%frac_mass),  real(disk_zone(i)%diskmass)
             endif  ! PAH
          endif ! pop dans la zone
       enddo ! pop

       fPAH(i) = max(fPAH(i),1e-9)
    enddo ! i


    ProDiMo_other_PAH = .false.
    if (NC_0 > 0) then
       ProDiMo_other_PAH = .true.
       ProDiMo_PAH_NC = NC
       ProDiMo_PAH_NH = NH
    endif

    allocate(ProDiMo_fPAH(n_regions), ProDiMo_dust_gas(n_regions), ProDiMo_Mdisk(n_regions))
    do ir=1, n_regions
       iz = regions(ir)%zones(1)

       ProDiMo_fPAH(ir) = fPAH(iz)
       ProDiMo_dust_gas(ir) = 1.0/disk_zone(iz)%gas_to_dust
       ProDiMo_Mdisk(ir) = disk_zone(iz)%diskmass * disk_zone(iz)%gas_to_dust

       do i=2,regions(ir)%n_zones
          iz = regions(ir)%zones(i)
          if ( abs(ProDiMo_fPAH(ir)-fPAH(iz)) > 1e-3*ProDiMo_fPAH(ir)) then
             do k=1, n_zones
                write(*,*) "zone ", k,  "region", disk_zone(k)%region, "fPAH =", fPAH(k)
             enddo ! k
             call error("fPAH must be contant within a region for ProDiMo")
          endif

          if ( abs(ProDiMo_dust_gas(ir)-1.0/disk_zone(iz)%gas_to_dust) > 1e-3*ProDiMo_dust_gas(ir)) then
             do k=1, n_zones
                write(*,*) "zone ", k, "region", disk_zone(k)%region, "gas_dust =", disk_zone(k)%gas_to_dust
             enddo ! k
             call error("gas_to_dust must be contant within a region for ProDiMo")
          endif

          ProDiMo_Mdisk(ir) = ProDiMo_Mdisk(ir) + disk_zone(iz)%diskmass * disk_zone(iz)%gas_to_dust
       enddo ! i

    enddo ! ir

    call create_input_ProDiMo()

    return

  end subroutine setup_ProDiMo


  !********************************************************************

  subroutine save_J_prodimo(lambda)
    ! C. Pinte
    ! 19/06/09
    ! sauvegarde le champ de radiation pour ProDiMo
    ! avant de calculer le champ ISM

    integer, intent(in) :: lambda
    integer :: ri, zj, phik, icell
    real(kind=dp) :: n_photons_envoyes, energie_photon, facteur

    ! Step2
    n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
    energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda)) / n_photons_envoyes &
         * tab_lambda(lambda) * 1.0e-6  !lambda.F_lambda  ! ICI

    do ri=1, n_rad
       do zj=j_start,nz
          if (zj==0) cycle
          phik=1
          icell = cell_map(ri,zj,phik)
          facteur = energie_photon / volume(icell)
          J_prodimo(lambda,ri,zj) =  J_prodimo(lambda,ri,zj) +  facteur * sum(xJ_abs(icell,lambda,:))
          N_ProDiMo(lambda,ri,zj) =  sum(xN_abs(icell,lambda,:))
       enddo
    enddo

    xJ_abs(:,lambda,:) = 0.0
    xN_abs(:,lambda,:) = 0.0

    return

  end subroutine save_J_prodimo

  !********************************************************************

  subroutine allocate_m2p(n_rad,nz,n_lambda)
    ! Routine similaire a celle de ProDiMo dans readMCFOST
    ! pour le moment, ne sert que pour recuperer le champ de radiation
    ! calcule par le run initial de mcfost
    ! C. Pinte
    ! 29/08/2012

    integer, intent(in) :: n_lambda, n_rad, nz

    integer, dimension(8) :: istat
    integer :: i

    istat = 0
    allocate(m2p%r_grid(n_rad,nz), &
           & m2p%z_grid(n_rad,nz),m2p%Tdust(n_rad,nz), &
           & m2p%density(n_rad,nz),m2p%dust_to_gas(n_rad,nz),stat=istat(2))
    !allocate(m2p%rho_grain(n_rad),stat=istat(3))
    allocate(m2p%wavelengths(n_lambda),m2p%nu(n_lambda),m2p%lamIlamStar(n_lambda), &
           & m2p%lamIlamISM(n_lambda),stat=istat(4))
    allocate(m2p%lamJlam(n_rad,nz,n_lambda),m2p%lamJlamStar(n_rad,nz,n_lambda), &
           & m2p%lamJlamISM(n_rad,nz,n_lambda),m2p%Jnu(n_rad,nz,n_lambda),stat=istat(5))
    allocate(m2p%nJ(n_rad,nz,n_lambda), m2p%nJ_Star(n_rad,nz,n_lambda), &
           & m2p%nJ_ISM(n_rad,nz,n_lambda), stat=istat(6))
    !allocate(m2p%kappa_ext(n_rad,nz,n_lambda),m2p%kappa_abs(n_rad,nz,n_lambda), &
    !       & stat=istat(7))
    !allocate(m2p%grain_size(n_rad,nz,0:3),stat=istat(8))
    do i=1,8
       if(istat(i) /= 0) then
          write(*,*) istat(i)," at ",i
          call error("memory allocation problem in allocate_m2p")
       endif
    enddo

    return

  end subroutine allocate_m2p

  !********************************************************************

  subroutine mcfost2ProDiMo()
    ! C. Pinte
    ! 18/02/09
    ! - cell position
    ! - dust temp
    ! - mean intensity Jnu
    !  ---> nbre de bins entiers entre 91.2 et 205nm
    ! - Qabs_nu et Qext_nu pour chaque cellule
    ! - gas density pour un rapport gas sur poussiere fixe ---> pour
    ! verification seleument
    ! - zero, first et second moment of the grain size distribution
    !  N = int f(a) da   <a^i> = 1/N * int int f(a) a^i da
    ! - separation of constribution

    use fits_utils, only : print_error

    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(5) :: naxes
    integer :: group,fpixel,nelements, alloc_status, lambda, ri, zj, l, icell, i, iRegion, k
    integer :: iPAH_start, iPAH_end, n_grains_PAH
    real (kind=dp) :: n_photons_envoyes, energie_photon, facteur, N
    real :: wl, norme, Ttmp

    logical :: simple, extend, lPAH_nRE
    character(len=512) :: filename
    character :: s

    real(kind=dp), dimension(n_rad,nz,2) :: grid
    real, dimension(n_rad,nz) :: dens
    real, dimension(n_rad,nz,0:3) :: N_grains
    real, dimension(n_lambda) :: spectre
    integer, dimension(n_rad) :: which_region

    character(len=512) :: para

    ! Allocation dynamique pour eviter d'utiliser la memeoire stack
    real(kind=dp), dimension(:,:,:), allocatable :: J_ProDiMo_ISM    ! n_lambda,n_rad,nz
    real, dimension(:,:,:), allocatable :: J_io ! n_rad,nz,n_lambda, joue aussi le role de N_io
    real, dimension(:,:,:,:), allocatable :: opacite ! (n_rad,nz,2,n_lambda)

    integer, dimension(:,:,:), allocatable :: is_eq
    real, dimension(:,:,:), allocatable :: TPAH_eq
    logical, dimension(n_grains_tot) :: mask_not_PAH
    real, dimension(:,:,:,:), allocatable :: P_TPAH

    call get_command_argument(1,para)

    lPAH_nRE = .false.
    test_PAH_nRE : do i=1, n_pop
       if (dust_pop(i)%is_PAH.and.(dust_pop(i)%methode_chauffage==3)) then
          lPAH_nRE = .true.
          exit test_PAH_nRE
       endif
    enddo test_PAH_nRE

    if (lPAH_nRE) then
       iPAH_start = grain_nRE_start
       iPAH_end = grain_nRE_end
       n_grains_PAH=n_grains_nRE
    else
       iPAH_start = grain_RE_nLTE_start
       iPAH_end = grain_RE_nLTE_end
       n_grains_PAH=n_grains_RE_nLTE
    endif

    mask_not_PAH(:) = .not.grain(:)%is_PAH

    !allocate(is_eq(n_rad,nz,iPAH_start:iPAH_end), stat=alloc_status)
    allocate(is_eq(n_rad,nz,1), TPAH_eq(n_rad,nz,1), P_TPAH(n_T,n_rad,nz,1), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error is_eq forProDiMo.fits.gz')
    is_eq = 0 ; TPAH_eq = 0.0 ; P_TPAH = 0.0

    allocate(opacite(n_rad,nz,2,n_lambda), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error opacite forProDiMo.fits.gz')
    opacite = 0

    allocate(J_io(n_rad,nz,n_lambda), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error J_io forProDiMo.fits.gz')
    J_io = 0

    allocate(J_ProDiMo_ISM(n_lambda,n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error J_ProDiMo_ISM forProDiMo.fits.gz')
    J_ProDiMo_ISM = 0

    filename = trim(data_ProDiMo)//"/"//trim(mcfost2ProDiMo_file)
    ! Initialize parameters about the FITS image
    simple=.true.
    extend=.true.
    group=1
    fpixel=1

    !  Get an unused Logical Unit Number to use to open the FITS file.
    status=0
    call ftgiou(unit,status)

    !  Create the new empty FITS file.
    blocksize=1
    call ftinit(unit,trim(filename),blocksize,status)

    !------------------------------------------------------------------------------
    ! HDU 1 : Grille
    !------------------------------------------------------------------------------
    bitpix=-64
    naxis=3
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=2
    nelements=naxes(1)*naxes(2)*naxes(3)

    ! Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Write optional keywords to the header
    call ftpkys(unit,'mcfost',mcfost_release,'',status)
    call ftpkys(unit,'mcfost_id',sha_id,'',status)
    call ftpkyj(unit,'mcfost2prodimo',mcfost2ProDiMo_version,'',status)
    call ftpkys(unit,'mcfost_model_name',trim(para),'',status)

    if (mcfost2ProDiMo_version >=3) then
       call ftpkyj(unit,'n_zones',n_zones,' ',status)
       call ftpkyj(unit,'n_regions',n_regions,' ',status)
    endif

    if (mcfost2ProDiMo_version >=5) then
       if (lPAH) then
          call ftpkyj(unit,'PAH_present',1,' ',status)
       else
          call ftpkyj(unit,'PAH_present',0,' ',status)
       endif

       if (lPAH_nRE) then
          call ftpkyj(unit,'PAH_force_eq',0,' ',status)
       else
          call ftpkyj(unit,'PAH_force_eq',1,' ',status)
       endif
    endif

    call ftpkye(unit,'Teff',etoile(1)%T,-8,'[K]',status)
    call ftpkye(unit,'Rstar',real(etoile(1)%r*AU_to_Rsun),-8,'[Rsun]',status)
    call ftpkye(unit,'Mstar',real(etoile(1)%M),-8,'[Msun]',status)
    call ftpkye(unit,'fUV',fUV_ProDiMo,-3,'',status)
    call ftpkye(unit,'slope_UV',slope_UV_ProDiMo+2.0,-3,'',status) ! on remet le 2
    call ftpkye(unit,'distance',real(distance),-8,'[pc]',status)

    if (mcfost2ProDiMo_version >=3) then
       call ftpkye(unit,'disk_dust_mass_tot',real(diskmass),-8,'[Msun]',status) ! Version >=3

       ! Adding regions information
       do i=1,n_regions
          write(s,'(i1)') i
          call ftpkye(unit,'Rmin_region_'//s,real(regions(i)%Rmin),-8,'[AU]',status)
          call ftpkye(unit,'Rmax_region_'//s,real(regions(i)%Rmax),-8,'[AU]',status)
       enddo
    endif ! mcfost2ProDiMo_version 3

    ! 1st zone or only zone if mcfost2ProDiMo_version <=2
    call ftpkye(unit,'disk_dust_mass',real(disk_zone(1)%diskmass),-8,'[Msun]',status)
    call ftpkye(unit,'Rin',real(disk_zone(1)%rmin),-8,'[AU]',status)
    call ftpkye(unit,'Rout',real(disk_zone(1)%rout),-8,'[AU]',status)
    call ftpkye(unit,'Rref',real(disk_zone(1)%rref),-8,'[AU]',status)
    call ftpkye(unit,'H0',real(disk_zone(1)%sclht),-8,'[AU]',status)
    call ftpkye(unit,'edge',real(disk_zone(1)%edge),-8,'[AU]',status)
    call ftpkye(unit,'beta',real(disk_zone(1)%exp_beta),-8,'',status)
    call ftpkye(unit,'alpha',real(disk_zone(1)%surf),-8,'',status)

    if (mcfost2ProDiMo_version >=3) then
       ! Adding zones information
       do i=2, n_zones
          write(s,'(i1)') i
          call ftpkye(unit,'disk_dust_mass_'//s,real(disk_zone(i)%diskmass),-8,'[Msun]',status)
          call ftpkye(unit,'Rin_'//s,real(disk_zone(i)%Rmin),-8,'[AU]',status)
          call ftpkye(unit,'Rout_'//s,real(disk_zone(i)%Rmax),-8,'[AU]',status)
          call ftpkye(unit,'Rref_'//s,real(disk_zone(i)%Rref),-8,'[AU]',status)
          call ftpkye(unit,'H0_'//s,real(disk_zone(i)%sclht),-8,'[AU]',status)
          call ftpkye(unit,'edge_'//s,real(disk_zone(i)%edge),-8,'[AU]',status)
          call ftpkye(unit,'beta_'//s,real(disk_zone(i)%exp_beta),-8,'',status)
          call ftpkye(unit,'alpha_'//s,real(disk_zone(i)%surf),-8,'',status)
       enddo
    endif ! mcfost2ProDiMo_version 3

    call ftpkye(unit,'amin',dust_pop(1)%amin,-8,'[micron]',status)
    call ftpkye(unit,'amax',dust_pop(1)%amax,-8,'[micron]',status)
    call ftpkye(unit,'aexp',dust_pop(1)%aexp,-8,'slope of grain size distribution',status)
    call ftpkye(unit,'strat',exp_strat,-8,'stratification exponent',status)
    call ftpkye(unit,'a_settle',a_strat,-8,'[micron]',status)
    call ftpkye(unit,'rho_grain',dust_pop(1)%rho1g_avg,-8,'[g.cm^-3] 1st pop',status)
    call ftpkys(unit,'optical_indices',trim(dust_pop(1)%indices(1)),'1st component 1st pop',status)

    call ftpkyj(unit,'n_grains',dust_pop(1)%n_grains,' ',status)
    call ftpkyj(unit,'n_rad',n_rad,' ',status)
    call ftpkyj(unit,'nz',nz,' ',status)
    call ftpkyj(unit,'n_rad_in',n_rad_in,' ',status)

    !  Write the array to the FITS file.
    do ri=1, n_rad
       do zj=1,nz
          icell = cell_map(ri,zj,1)
          grid(ri,:,1) = r_grid(icell)
          grid(ri,zj,2) = z_grid(icell)
       enddo
    enddo

    call ftpprd(unit,group,fpixel,nelements,grid,status)

    !------------------------------------------------------------------------------
    ! HDU 2: Temperature
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=2
    naxes(1)=n_rad
    naxes(2)=nz
    nelements=naxes(1)*naxes(2)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    norme = sum(r_grain(:)**2 * nbre_grains(:))
    !if (lRE_nLTE) then
    !   do ri=1, n_rad
    !      do zj=1, nz
    !         Temperature(ri,zj,1) = sum( Temperature_1grain(ri,zj,:) * r_grain(:)**2 * nbre_grains(:)) / norme
    !      enddo !j
    !   enddo !i
    !
    !   write(*,*) "************************************************************"
    !   write(*,*) "WARNING: nLTE mode in MCFOST"
    !   write(*,*) "Crude averaging of the temperature distribution for ProDiMo"
    !   write(*,*) "T = \int T(a) * n(a) * a^2 da / \int n(a) * a^2 da"
    !   write(*,*) "Max. averaged temperature = ", maxval(Temperature(:,:,1))
    !   write(*,*) "************************************************************"
    !
    !endif

    !  Write the array to the FITS file.
    call ftppre(unit,group,fpixel,nelements,Temperature,status)

    !------------------------------------------------------------------------------
    ! HDU 3 : Longueurs d'onde
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=1
    naxes(1)=n_lambda
    nelements=naxes(1)

    ! create new hdu
    call ftcrhd(unit, status)

    ! Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Write the array to the FITS file.
    call ftppre(unit,group,fpixel,nelements,real(tab_lambda,kind=sp),status)

    !------------------------------------------------------------------------------
    ! HDU 4 : Spectre stellaire
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=1
    naxes(1)=n_lambda
    nelements=naxes(1)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Verif Spectre stellaire
    !write(*,*) sum(spectre_etoiles) / (sigma * etoile(1)%T**4 * 4 * pi * (etoile(1)%r * AU_to_m)**2 )

    ! Conversion en lambda.F_lamda
    ! Division par 4 * pi * (etoile(1)%r * AU_to_m)**2 pour passer de luminosite a flux
    ! Division par pi pour passer en intensite
    spectre(:) = spectre_etoiles(:) * tab_lambda(:) / tab_delta_lambda(:) &
         / (4 * pi * (etoile(1)%r * AU_to_m)**2) / pi

    !  Write the array to the FITS file.
    call ftppre(unit,group,fpixel,nelements,spectre,status)

    !------------------------------------------------------------------------------
    ! HDU 5 : Spectre ISM (input)
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=1
    naxes(1)=n_lambda
    nelements=naxes(1)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    do lambda=1,n_lambda
       wl = tab_lambda(lambda) * 1e-6
       spectre(lambda) = (chi_ISM * 1.71 * Wdil * Blambda(wl,T_ISM_stars) + Blambda(wl,TCmb)) * wl
    enddo

    !  Write the array to the FITS file.
    call ftppre(unit,group,fpixel,nelements,spectre,status)

    !------------------------------------------------------------------------------
    ! HDU 6 : Champ de radiation en W.m-2 (lambda.F_lambda)
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=3
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=n_lambda
    nelements=naxes(1)*naxes(2)*naxes(3)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Inversion de l'ordre des dimensions + passage en simple precision
    do ri=1, n_rad
       do zj=1,nz
          J_io(ri,zj,:) = J_ProDiMo(:,ri,zj)
       enddo
    enddo

    call ftppre(unit,group,fpixel,nelements,J_io,status)

    !------------------------------------------------------------------------------
    ! HDU 7 : Statistique du champ de radiation (nombre de paquet)
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=3
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=n_lambda
    nelements=naxes(1)*naxes(2)*naxes(3)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Inversion de l'ordre des dimensions
    do zj=1,nz
       do ri=1, n_rad
          J_io(ri,zj,:) = N_ProDiMo(:,ri,zj)
       enddo
    enddo

    call ftppre(unit,group,fpixel,nelements,J_io,status)

    !------------------------------------------------------------------------------
    ! HDU 8 : Champ de radiation ISM en W.m-2 (lambda.F_lambda)
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=3
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=n_lambda
    nelements=naxes(1)*naxes(2)*naxes(3)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)


    do lambda=1, n_lambda
       n_photons_envoyes = sum(n_phot_envoyes_ISM(lambda,:))

       wl = tab_lambda(lambda) * 1e-6
       energie_photon = (chi_ISM * 1.71 * Wdil * Blambda(wl,T_ISM_stars) + Blambda(wl,TCmb)) * wl & !lambda.F_lambda
           * (4.*pi*(R_ISM*Rmax)**2) / n_photons_envoyes / pi  ! ici

       do ri=1, n_rad
          do zj=1,nz
             icell = cell_map(ri,zj,1)
             facteur = energie_photon / volume(icell)
             J_prodimo_ISM(lambda,ri,zj) =  J_prodimo_ISM(lambda,ri,zj) +  facteur * sum(xJ_abs(icell,lambda,:))
          enddo
       enddo
    enddo ! lambda



    ! Inversion de l'ordre des dimensions + passage en simple precision
    do ri=1, n_rad
       do zj=1,nz
          J_io(ri,zj,:) = J_ProDiMo_ISM(:,ri,zj)
       enddo
    enddo

    call ftppre(unit,group,fpixel,nelements,J_io,status)

    !------------------------------------------------------------------------------
    ! HDU 9 : Statistique du champ de radiation ISM (nombre de paquet)
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=3
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=n_lambda
    nelements=naxes(1)*naxes(2)*naxes(3)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! Inversion de l'ordre des dimensions et somaation
    do zj=1,nz
       do ri=1, n_rad
          icell = cell_map(ri,zj,1)
          J_io(ri,zj,:) = sum(xN_abs(icell,:,:), dim=2)
       enddo
    enddo

    call ftppre(unit,group,fpixel,nelements,J_io,status)


    !------------------------------------------------------------------------------
    ! HDU 10 : Densite de gaz [g.cm^-3]
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=2
    naxes(1)=n_rad
    naxes(2)=nz
    nelements=naxes(1)*naxes(2)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    !  Write the array to the FITS file.
    do ri=1, n_rad
       do zj=1,nz
          icell = cell_map(ri,zj,1)
          dens(ri,zj) =  densite_gaz(icell) * masse_mol_gaz / m3_to_cm3 ! g.cm^-3
       enddo
    enddo

    call ftppre(unit,group,fpixel,nelements,dens,status)

    !------------------------------------------------------------------------------
    ! HDU 11 : Opacites
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=4
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=2
    naxes(4)=n_lambda
    nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
    ! C_ext = a² microns² -> 1e-8 cm²             \
    ! density en cm-3                      > reste facteur 149595.0
    ! longueur de vol en AU = 1.5e13 cm   /
    facteur = AU_to_cm * mum_to_cm**2
    opacite = 0.0
    do zj=1,nz
       do ri=1,n_rad
          icell = cell_map(ri,zj,1)
          do lambda=1,n_lambda
             do l= grain_RE_LTE_start, grain_RE_LTE_end
                opacite(ri,zj,1,lambda) = opacite(ri,zj,1,lambda) + C_ext(l,lambda) * densite_pouss(l,icell) * facteur
                opacite(ri,zj,2,lambda) = opacite(ri,zj,2,lambda) + C_abs(l,lambda) * densite_pouss(l,icell) * facteur
             enddo ! l
          enddo ! lambda
       enddo ! ri
    enddo !zj

    call ftppre(unit,group,fpixel,nelements,opacite,status)

    !------------------------------------------------------------------------------
    ! HDU 12 : Moments de la distribution en tailles des grains
    !------------------------------------------------------------------------------
    bitpix=-32
    naxis=3
    naxes(1)=n_rad
    naxes(2)=nz
    naxes(3)=4
    nelements=naxes(1)*naxes(2)*naxes(3)

    ! create new hdu
    call ftcrhd(unit, status)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

    do zj=1,nz
       do ri=1,n_rad
          ! Nbre total de grain : le da est deja dans densite_pouss
          icell = cell_map(ri,zj,1)
          N = sum(densite_pouss(:,icell),mask=mask_not_PAH)
          N_grains(ri,zj,0) = N
          if (N > 0) then
             N_grains(ri,zj,1) = sum(densite_pouss(:,icell) * r_grain(:),mask=mask_not_PAH) / N
             N_grains(ri,zj,2) = sum(densite_pouss(:,icell) * r_grain(:)**2,mask=mask_not_PAH) / N
             N_grains(ri,zj,3) = sum(densite_pouss(:,icell) * r_grain(:)**3,mask=mask_not_PAH) / N
          else
             N_grains(ri,zj,1) = 0.0
             N_grains(ri,zj,2) = 0.0
             N_grains(ri,zj,3) = 0.0
          endif
       enddo
    enddo

    ! part.cm^-3 --> part.m^-3
    N_grains(:,:,0) = N_grains(:,:,0) /  (cm_to_m**3)
    call ftppre(unit,group,fpixel,nelements,N_grains,status)


    if (mcfost2ProDiMo_version >=3) then
       !------------------------------------------------------------------------------
       ! HDU 13 : Region index
       !------------------------------------------------------------------------------
       bitpix=32
       naxis=1
       naxes(1)=n_rad
       nelements=naxes(1)

       ! create new hdu
       call ftcrhd(unit, status)

       ! Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       which_region = 0
       do iRegion=1,n_regions
          do i=1, n_rad
             icell = cell_map(i,1,1)
             if ( (r_grid(icell) > regions(iRegion)%Rmin).and.(r_grid(icell) < regions(iRegion)%Rmax) ) which_region(i) = iRegion
          enddo !i
       enddo ! iRegion
       call ftpprj(unit,group,fpixel,nelements,which_region,status)

    endif ! mcfost2ProDiMo_version >=3

    if (mcfost2ProDiMo_version >=4) then
       !------------------------------------------------------------------------------
       ! HDU 14 : Spectre stellaire HR
       !------------------------------------------------------------------------------
       bitpix=-32
       naxis=2
       naxes(1:2)=shape(ProDiMo_star_HR)
       nelements=naxes(1)*naxes(2)

       ! create new hdu
       call ftcrhd(unit, status)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       ! Conversion en lambda.F_lamda
       ! Division par 4 * pi * (etoile(1)%r * AU_to_m)**2 pour passer de luminosite a flux
       ! Division par pi pour passer en intensite
       ProDiMo_star_HR(:,2) = ProDiMo_star_HR(:,2) &
            / (4 * pi * (etoile(1)%r * AU_to_m)**2) / pi

       !  Write the array to the FITS file.
       call ftppre(unit,group,fpixel,nelements,ProDiMo_star_HR,status)
    endif

    if ((mcfost2ProDiMo_version >=5).and.lPAH) then

       !------------------------------------------------------------------------------
       ! HDU 15 : PAH density
       !------------------------------------------------------------------------------
       bitpix=-32
       naxis=2
       naxes(1)=n_rad
       naxes(2)=nz
       nelements=naxes(1)*naxes(2)

       ! create new hdu
       call ftcrhd(unit, status)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       ! Write  optional keywords to the header
       call ftpkys(unit,'UNIT',"g.cm^-3",' ',status)

       !  Write the array to the FITS file.
       k=1 ! azimuth
       do ri=1, n_rad
          do zj=1,nz
             icell = cell_map(ri,zj,k)
             dens(ri,zj) = sum(densite_pouss(iPAH_start:iPAH_end, icell) * M_grain(iPAH_start:iPAH_end)) ! M_grain en g
          enddo
       enddo
       call ftppre(unit,group,fpixel,nelements,dens,status)

       !------------------------------------------------------------------------------
       ! HDU 16 : total PAH opacity
       !------------------------------------------------------------------------------
       bitpix=-32
       naxis=4
       naxes(1)=n_rad
       naxes(2)=nz
       naxes(3)=2
       naxes(4)=n_lambda
       nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

       ! create new hdu
       call ftcrhd(unit, status)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       opacite = 0.0
       do zj=1,nz
          do ri=1,n_rad
             icell = cell_map(ri,zj,1)
             do lambda=1,n_lambda
                do l= iPAH_start, iPAH_end
                   opacite(ri,zj,1,lambda) = opacite(ri,zj,1,lambda) + C_ext(l,lambda) * densite_pouss(l,icell)
                   opacite(ri,zj,2,lambda) = opacite(ri,zj,2,lambda) + C_abs(l,lambda) * densite_pouss(l,icell)
                enddo ! l
             enddo ! lambda
          enddo ! ri
       enddo !zj
       call ftppre(unit,group,fpixel,nelements,opacite,status)

       !------------------------------------------------------------------------------
       ! HDU 17 : PAH Teq
       !------------------------------------------------------------------------------
       bitpix=-32
       naxis=3
       naxes(1)=n_rad
       naxes(2)=nz
       naxes(3)= 1 !n_grains_PAH
       nelements=naxes(1)*naxes(2)*naxes(3)

       ! create new hdu
       call ftcrhd(unit, status)

       ! Average temperature at the moment, as ProDiMo can only handle 1 PAH so far
       TPAH_eq = 0.0
       do zj=1,nz
          do ri=1,n_rad
             icell = cell_map(ri,zj,1)
             norme = 0.0
             do l= iPAH_start, iPAH_end
                if (lPAH_nRE) then
                   Ttmp = temperature_1grain_nRE(l,icell)
                else
                   Ttmp = temperature_1grain(l,icell)
                endif
                TPAH_eq(ri,zj,1) = TPAH_eq(ri,zj,1) + Ttmp**4 * densite_pouss(l,icell)
                norme = norme + densite_pouss(l,icell)
             enddo ! l
             TPAH_eq(ri,zj,1) = (TPAH_eq(ri,zj,1)/norme)**0.25
          enddo ! ri
       enddo !zj

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       call ftppre(unit,group,fpixel,nelements,TPAH_eq,status)

       !------------------------------------------------------------------------------
       ! HDU 18 : is PAH at equilibrium
       !------------------------------------------------------------------------------
       bitpix=32
       naxis=3
       naxes(1)=n_rad
       naxes(2)=nz
       naxes(3)=1 !n_grains_PAH
       nelements=naxes(1)*naxes(2)*naxes(3)

       ! create new hdu
       call ftcrhd(unit, status)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

       !  Write the array to the FITS file.
       group=1
       fpixel=1
       nelements=naxes(1)*naxes(2)*naxes(3)

       ! le j signifie integer
       if (lPAH_nRE) then
          do ri=1, n_rad
             do zj=1,nz
                icell = cell_map(ri,zj,1)
                do l=grain_nRE_start, grain_nRE_end
                   if (l_RE(l,icell)) then
                      is_eq(ri,zj,1) = 1
                   else
                      is_eq(ri,zj,1) = 0
                   endif
                enddo
             enddo
          enddo
       else
          is_eq(:,:,1) = 1
       endif
       call ftpprj(unit,group,fpixel,nelements,is_eq,status)

       if (lPAH_nRE) then
          !------------------------------------------------------------------------------
          ! HDU 19 : temperature table
          !------------------------------------------------------------------------------
          bitpix=-32
          naxis=1
          naxes(1)=n_T
          nelements=naxes(1)

          ! create new hdu
          call ftcrhd(unit, status)

          !  Write the required header keywords.
          call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

          ! le e signifie real*4
          call ftppre(unit,group,fpixel,nelements,tab_Temp,status)

          !------------------------------------------------------------------------------
          ! HDU 20 : PAH temperature probability density
          !------------------------------------------------------------------------------
          bitpix=-32
          naxis=4
          naxes(1)=n_T
          naxes(2)=n_rad
          naxes(3)=nz
          naxes(4)=1 !n_grains_nRE
          nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

          ! create new hdu
          call ftcrhd(unit, status)

          !  Write the required header keywords.
          call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

          ! Trick to moved pack all the grain on index 1 using grain index = tab_region(ri)
          do ri=1, n_rad
             do zj=1, nz
                icell = cell_map(ri,zj,1)
                if (tab_region(i) > 0) then
                   P_TPAH(:,ri,zj,1) = Proba_Temperature(:,tab_region(ri),icell)
                else
                   P_TPAH(:,ri,zj,1) = 0.0
                endif
             enddo ! j
          enddo ! j

           ! le e signifie real*4
          call ftppre(unit,group,fpixel,nelements,P_TPAH,status)
       endif !lPAH_nRE

    endif

    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

    !  Check for any error, and if so print out error messages
    if (status > 0) then
       call print_error(status)
       call exit(1)
    endif

    return

  end subroutine mcfost2ProDiMo

  !********************************************************************

  subroutine create_input_ProDiMo()
    ! Create the *.in files for ProDiMo
    ! C. Pinte 17/03/2012

    character(len=512) :: cmd

    character(len=128) :: line
    character(len=10) :: fmt

    integer :: status, i, iProDiMo, syst_status, ipop
    logical :: unchanged

    character :: s

    write (*,*) 'Creating directory '//trim(data_ProDiMo)
    ! Copy *.in files ()
    cmd = 'mkdir -p '//trim(data_ProDiMo)//" ; "// &
         "cp "//trim(ProDiMo_input_dir)//"/Elements.in "&
         //trim(ProDiMo_input_dir)//"/LineTransferList.in "&
         //trim(ProDiMo_input_dir)//"/Reactions.in "&
         //trim(ProDiMo_input_dir)//"/Species.in "//trim(data_ProDiMo)
    call appel_syst(cmd,syst_status)


    ! Copy and modify Parameter.in
    open(unit=1,file=trim(ProDiMo_input_dir)//"/Parameter.in",status='old')
    open(unit=2,file=trim(data_ProDiMo)//"/Parameter.in",status='replace')

    status = 0
    read_loop : do
       read(1,'(A128)',iostat=status) line
       if (status /= 0) exit read_loop
       unchanged = .true.

       ! Rin, Rout, Mdisk, dust_to_gas ratio and fPAH
       if (INDEX(line,"! Rout") > 0) then
          write(2,*)  n_regions, " ! NZONES : number of MCFOST regions"
          ! 1st ProDiMo zone is the last region (no number)
          if (n_regions==1) then
             write(2,'(A)') "----- single ProDiMo-zone ----- : set by MCFOST"
          else
             write(2,'(A)') "------ outer ProDiMo-zone ----- : MCFOST region"
          endif
          write(2,*) real(regions(n_regions)%Rmin)," ! Rin  : set by MCFOST"
          write(2,*) real(regions(n_regions)%Rmax)," ! Rout : set by MCFOST"
          write(2,*) ProDiMo_dust_gas(n_regions),  " ! dust_to_gas : set by MCFOST"
          write(2,*) ProDiMo_fPAH(n_regions),      " ! fPAH : set by MCFOST"
          write(2,*) ProDiMo_Mdisk(n_regions),     " ! Mdisk : set by MCFOST"

          do i=1, n_regions-1
             iProDiMo = i+1
             write(s,'(i1)') iProDiMo
             write(2,'(A)') "------ ProDiMo-zone #"//s//"----- : set by MCFOST"
             write(2,*) real(regions(i)%Rmin)," ! R"//s//"in  : weird ProDiMo order, set by MCFOST"
             write(2,*) real(regions(i)%Rmax)," ! R"//s//"out : set by MCFOST"
             write(2,*) ProDiMo_dust_gas(i), " ! d"//s//"ust_to_gas : set by MCFOST"
             write(2,*) ProDiMo_fPAH(i),     " ! f"//s//"PAH : set by MCFOST"
             write(2,*) ProDiMo_Mdisk(i),    " ! M"//s//"disk : set by MCFOST"
          enddo

          write(2,*) ProDiMo_other_PAH, " ! other_PAH : set by MCFOST"
          if (ProDiMo_other_PAH) then
             write(2,*) ProDiMo_PAH_NC, " ! PAH_NC : set by MCFOST"
             write(2,*) ProDiMo_PAH_NH, " ! PAH_NH : set by MCFOST"
          endif

          unchanged = .false.
       endif
       ! On zappe Rin, dust_to_gas, fPAH et Mdisk
       if ( (INDEX(line,"! Rin") > 0).or.(INDEX(line,"! dust_to_gas") > 0).or. &
          (INDEX(line,"! fPAH") > 0).or.(INDEX(line,"! Mdisk") > 0) ) then
          unchanged = .false.
       endif

       if (INDEX(line,"! Rphoto_bandint") > 0) then
          if ((etoile(1)%T > 6000.).or.(etoile(1)%fUV > 0.05)) then
             write(2,*) .true.,  " ! Rphoto_bandint : set by MCFOST"
          else
             write(2,*) .false., " ! Rphoto_bandint : set by MCFOST"
          endif
          unchanged = .false.
       endif

       if (INDEX(line,"! PAH_in_RT") > 0) then
          lPAH = .false.
          do ipop=1, n_pop
             if (dust_pop(ipop)%is_PAH) then
                lPAH = .true.
             endif
          enddo
          write(2,*) lPAH, " ! PAH_in_RT : set by MCFOST"
          unchanged = .false.
       endif


       if (INDEX(line,"! v_turb") > 0) then
          write(2,*) real(vitesse_turb * m_to_km), " ! v_turb [km/s] : set by MCFOST"
          unchanged = .false.
       endif

       if (INDEX(line,"! incl") > 0) then
          write(2,*) RT_imin, " ! incl : set by MCFOST (1st MCFOST inclination)" ! tab_RT_incl n'est pas necessairement alloue
          unchanged = .false.
       endif

       if (INDEX(line,"! dist") > 0) then
          write(2,*) distance, " ! dist : set by MCFOST"
          unchanged = .false.
       endif

       if (unchanged) then
          if (len(trim(line)) > 0) then
             write(fmt,'("(A",I3.3,")")') len(trim(line))
             write(2,fmt) trim(line)
          else
             write(2,*) ""
          endif
       endif
    enddo read_loop

    if (etoile(1)%lb_body) then
       write(2,*) ".true.  ! PlanckSpec : set by MCFOST"
    else
       write(2,*) ".false.  ! PlanckSpec : set by MCFOST"
    endif

    close(1)
    close(2)

    return

  end subroutine create_input_ProDiMo

  !********************************************************************

  subroutine read_mcfost2ProDiMo(para)
    ! Relit le fichier for ProDiMo.fits.gz cree par mcfost pour ProDiMo
    ! afin de redemarrer dessus pour le transfert dans les raies
    ! Lit les parametres et la structure en temperature
    ! A terme le champ de radiation pour calculer le contribution de
    ! lumiere diffusee dans les raies
    !
    ! C.Pinte
    ! 12/07/09

    character(len=*) :: para

    integer :: status, readwrite, unit, blocksize, nfound, group, firstpix, npixels, hdutype
    integer :: imol, syst_status, n_files, n_rad_m2p, nz_m2p, n_lambda_m2p, i, j, l
    real :: nullval, Tdust
    logical :: anynull

    integer, dimension(4) :: naxes

    character(len=30) :: errtext
    character (len=80) :: errmessage
    character(len=512) :: filename, cmd

    !**************************************************************
    ! 1) Lecture des parametres dans le fichier .para de data_th
    !**************************************************************

    ! Copie temporaire et lecture du fichier de parametres
    ! car je ne connais pas le nom du fichier .par dans data_th
    n_files = 0
    cmd = "ls data_th/*.par* | wc -l > n_files.tmp" ; call appel_syst(cmd, syst_status)
    open(unit=1, file="n_files.tmp", status='old')
    read(1,*) n_files
    close(unit=1,status="delete")

    if (n_files > 1) call error("There are more than 1 parameter file in data_th")
    if (n_files < 1) call error("There are less than 1 parameter file in data_th")

    cmd = "cp data_th/*.par* data_th/forMCFOST.par" ; call appel_syst(cmd, syst_status)
    call read_para("data_th/forMCFOST.par")
    cmd = "rm -rf data_th/forMCFOST.par" ; call appel_syst(cmd, syst_status)

    ! Parametres par defaut
    ltemp = .false.
    lsed = .false.
    lsed_complete = .false.
    lemission_mol = .true.

    ! on ne force plus depuis version 2 de l'interface
    !n_lambda = 39
    !lambda_min = 0.0912
    !lambda_max = 3410.85

    deallocate(mol) ! car on va le remplacer

    open(unit=1, file=para, status='old')

    ! Inclinaisons
    read(1,*)
    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered

    read(1,*)
    read(1,*)
    !Molecules
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_center_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%n_extraV_rt, mol(imol)%extra_deltaV_rt
       mol(imol)%extra_deltaV_rt = mol(imol)%extra_deltaV_rt * km_to_m ! Conversion en m.s-1
       mol(imol)%n_speed_rt = mol(imol)%n_speed_center_rt + mol(imol)%n_extraV_rt
       mol(imol)%lcst_abundance = .true.
       mol(imol)%lline = .true.
       read(1,*) mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       !mol(imol)%n_speed = 1 ! inutilise si on ne calcule pas le NLTE


       mol(imol)%abundance = 1e-6
       mol(imol)%iLevel_max = maxval(mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)) + 2
    enddo
    vitesse_turb = 0.0
    largeur_profile = 15.

    !lpop = .false. ; lprecise_pop = .false. ; lmol_LTE = .true. ! lmol_LTE force l'utilisation des pop de ProDiMo

    if (lpop) then
       write(*,*) "Calculating level population"
    endif
    if (lmol_LTE) then
       write(*,*) "Using ProDiMo levels"
    endif

    para_version = mcfost_version
    close(unit=1)


    !******************************************************************
    ! 2) Lecture du champ de radiation continu dans forProDiMo.fits.gz
    !******************************************************************
    filename = "data_ProDiMo/"//trim(mcfost2ProDiMo_file)

    write(*,*) "Reading "//trim(filename)

    status=0
    group=1
    firstpix=1
    nullval=-999

    !  Get an unused Logical Unit Number to use to open the FITS file.
    call ftgiou(unit,status)

    ! Open file
    readwrite=0 ; blocksize = 1
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status /= 0) call error(trim(mcfost2ProDiMo_file)//" file needed")

    !------------------------------------------------------------------------
    ! HDU 6 : Read dimensions : radiation field(n_rad,nz,n_lambda)
    !------------------------------------------------------------------------
    call ftmahd(unit,6,hdutype,status)
    call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
    if (nfound /= 3) call error('failed to read the NAXISn keyworrds of HDU 6')
    n_rad_m2p    = naxes(1)
    nz_m2p       = naxes(2)
    n_lambda_m2p = naxes(3)

    if ( (n_rad_m2p /= n_rad) .or. (nz_m2p /= nz) ) then
       write(*,*) "forProDiMo.fits.gz: n_rad=", n_rad_m2p, "nz=", nz_m2p
       call error("Spatial grid if forProDiMo.fits.gz must be the same as the mcfost's one")
    endif

    call allocate_m2p(n_rad_m2p,nz_m2p,n_lambda_m2p)

    !------------------------------------------------------------------------
    ! HDU 6: Radiation field caused by stellar photons
    !------------------------------------------------------------------------
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,m2p%lamJlamStar,anynull,status)

    !------------------------------------------------------------------------
    ! HDU 2: Temperature
    !------------------------------------------------------------------------
    !  move to hdu 2
    call ftmahd(unit,2,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
    if (nfound /= 2) call error("failed to read the NAXISn keywords of HDU 2")
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) call error("HDU 2 does not have the right dimensions.")
    npixels=naxes(1)*naxes(2)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,m2p%Tdust,anynull,status)


    !------------------------------------------------------------------------
    ! HDU 3: Wavelengths
    !------------------------------------------------------------------------
    !  move to hdu 3
    call ftmahd(unit,3,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
    if (nfound /= 1) call error("failed to read the NAXISn keywords of HDU 3")
    if (naxes(1) /= n_lambda_m2p) call error("HDU 3 does not have the right dimensions.")
    npixels=naxes(1)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,m2p%wavelengths,anynull,status)
    m2p%nu(:) = c_light * m_to_mum / m2p%wavelengths(:) ! Hz

    !--------------------------------------------------------------------------
    ! HDU 7 : Statistic of the radiation fied
    !--------------------------------------------------------------------------
    !  move to hdu 7
    call ftmahd(unit,7,hdutype,status)


    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
    if (nfound /= 3) call error("failed to read the NAXISn keywords of HDU 7")
    if ((naxes(1) /= n_rad_m2p).or.(naxes(2) /= nz_m2p).or.(naxes(3) /= n_lambda_m2p))  call error("HDU 7 does not have the right dimensions.")
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,m2p%nJ_Star,anynull,status)

    !--------------------------------------------------------------------------
    ! HDU 8 : Radiation field caused by ISM photons
    !--------------------------------------------------------------------------
    !  move to hdu 8
    call ftmahd(unit,8,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
    if (nfound /= 3) call error("failed to read the NAXISn keywords of HDU 8")
    if ((naxes(1) /= n_rad_m2p).or.(naxes(2) /= nz_m2p).or.(naxes(3) /= n_lambda_m2p)) call error("HDU 8 does not have the right dimensions.")
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,m2p%lamJlamISM,anynull,status)

    !--------------------------------------------------------------------------
    ! HDU 9 : Statistic of the ISM radiation field
    !--------------------------------------------------------------------------
    !  move to hdu 9
    call ftmahd(unit,9,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
    if (nfound /= 3) call error("failed to read the NAXISn keywords of HDU 9")
    if ((naxes(1) /= n_rad_m2p).or.(naxes(2) /= nz_m2p).or.(naxes(3) /= n_lambda_m2p)) call error("HDU 9 does not have the right dimensions.")
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,m2p%nJ_ISM,anynull,status)

    ! close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

    !  Check for any error, and if so print out error messages
    !  Get the text string which describes the error
    if (status > 0) then
       call ftgerr(status,errtext)
       write(*,*) "ERROR : "//trim(filename)
       write(*,*) 'FITSIO Error Status =',status,': ',errtext

       !  Read and print out all the error messages on the FITSIO stack
       call ftgmsg(errmessage)
       do while (errmessage .ne. ' ')
          write(*,*) errmessage
          call ftgmsg(errmessage)
       end do
    endif


    !---------------------------------------------------------------------
    ! Suppress Jnu for low-number statistics : same criteria as in ProDiMo
    !---------------------------------------------------------------------
    m2p%nJ = m2p%nJ_Star + m2p%nJ_ISM
    do i=1,n_rad
       do j=1,nz
          Tdust = m2p%Tdust(i,j)
          do l=1,n_lambda_m2p
             if (m2p%nJ(i,j,l) < 10) then
                m2p%lamJlamStar(i,j,l) = m2p%nu(l) * Bnu(m2p%nu(l)*1.0_dp,Tdust)
                m2p%lamJlamISM(i,j,l) = 0.0
          endif
        enddo  !l
      enddo ! j
    enddo ! i
    m2p%lamJlam = m2p%lamJlamStar + m2p%lamJlamISM

    do i=1,n_rad
       do j=1,nz
          m2p%Jnu(i,j,:) = m2p%lamJlam(i,j,:) / m2p%nu(:)
       enddo !j
    enddo !i

    return

  end subroutine read_mcfost2ProDiMo

  !********************************************************************

  subroutine read_ProDiMo2mcfost(imol)
    ! C. Pinte  17/09/09
    ! Lit les resultats de ProDiMo et met a jour les valeurs dans mcfost
    ! -->  Tgas, abondances a terme
    ! population et vitesses pour le moment

    integer, intent(in) :: imol

    integer :: fits_status, readwrite, unit, blocksize, nfound, group, firstpix, npixels, hdutype,alloc_status
    real :: nullval
    logical :: anynull
    integer, dimension(4) :: naxes
    character(len=30) errtext
    character (len=80) errmessage, comment
    character(len=512) :: filename

    character(len=40) :: used_sha_id
    character(len=8) :: used_mcfost_version
    integer :: used_mcfost2ProDiMo_version


    logical, save :: l_first_time = .true.

    real(kind=dp), dimension(MC_NSP,n_rad,nz) :: TMP
    real(kind=dp), dimension(n_rad,nz) :: Tgas ! Pas utilise pour le moment, pour futurs calculs NLTE
    real, dimension(n_rad,nz,2) :: grid ! Seulement pout test
    real, dimension(n_rad,nz) :: sum_pops
    real, dimension(:,:,:), allocatable :: MCpops

    logical :: lCII, lOI, lCO, loH2O, lpH2O

    real :: sigma2, sigma2_m1, r_sph, factor_CO
    integer :: i,j, ri, zj, n_speed_rt, l, keyword_status, icell

    real, parameter :: ratio_12_13_CO = 70.
    real, parameter :: ratio_12_18_CO = 500.

    n_speed_rt = mol(imol)%n_speed_rt

    lCII = .false. ; lOI = .false. ; lCO = .false. ; loH2O = .false. ;  lpH2O = .false. ;
    if (mol(imol)%name=="C+") lCII = .true.
    if (mol(imol)%name=="O") lOI = .true.
    if (mol(imol)%name=="CO") then
       lCO = .true.
       factor_CO = 1.
    endif
    if (mol(imol)%name=="13C16O") then
       lCO = .true.
       factor_CO = 1./ ratio_12_13_CO
    endif
    if (mol(imol)%name=="C18O") then
       lCO = .true.
       factor_CO = 1./ ratio_12_18_CO
    endif
    if (mol(imol)%name=="o-H2O") loH2O = .true.
    if (mol(imol)%name=="p-H2O") lpH2O = .true.

    ldust_mol = .true.

    if (l_first_time) then
       l_first_time = .false.

       filename = "data_ProDiMo/forMCFOST.fits.gz"
       write(*,*)
       write(*,*) "Reading ProDiMo calculations"
       write(*,*) trim(filename)

       fits_status = 0
       !  Get an unused Logical Unit Number to use to open the FITS file.
       call ftgiou(unit,fits_status)

       ! Open file
       readwrite=0
       call ftopen(unit,filename,readwrite,blocksize,fits_status)
       if (fits_status > 0) call error(trim(filename)//" file needed")

       group=1
       firstpix=1
       nullval=-999

       !---------------------------------------------------------
       ! HDU 1 : MCgrid, to make a test inside mcfost
       ! Format is the same as for the mcfost2ProDiMo interface
       !---------------------------------------------------------
       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,10,naxes,nfound,fits_status)
       if (nfound /= 3) call error("failed to read the NAXISn keywords of HDU 1")
       if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= 2)) call error("HDU 1 does not have the right dimensions.")
       npixels=naxes(1)*naxes(2)*naxes(3)

       ! reading keywords
       used_mcfost_version = 'unknown'
       comment = ' '
       used_sha_id = 'unknown'
       used_mcfost2ProDiMo_version = 0
       ProDiMo2mcfost_version = 0

       ! Read model info
       keyword_status = 0
       call ftgkys(unit,'mcfost',used_mcfost_version,comment,keyword_status)
       if (keyword_status /= 0) then
          write(*,*) "Fits keyword error: mcfost"
          keyword_status = 0
       endif
       call ftgkys(unit,'mcfost_id',used_sha_id,comment,keyword_status)
       if (keyword_status /= 0) then
          write(*,*) "Fits keyword error: mcfost_id"
          keyword_status = 0
       endif
       call ftgkyj(unit,'mcfost2prodimo',used_mcfost2ProDiMo_version,comment,keyword_status)
       if (keyword_status /= 0) then
          write(*,*) "Fits keyword error: mcfost2prodimo"
          keyword_status = 0
       endif
       call ftgkyj(unit,'prodimo2mcfost',ProDiMo2mcfost_version,comment,keyword_status)
       if (keyword_status /= 0) then
          write(*,*) "Fits keyword error: prodimo2mcfost"
          keyword_status = 0
       endif

       if (used_mcfost_version /= mcfost_release) then
          write(*,*) "************************************************"
          write(*,*) "WARNING: MCFOST has been updated"
          write(*,*) "- initial model was computed with mcfost "//trim(used_mcfost_version)
          write(*,*) "id="//used_sha_id
          write(*,*) "- this is mcfost "//trim(mcfost_release)
          write(*,*) "id="//sha_id
          write(*,*) "Problems can occur"
          write(*,*) "************************************************"
       else
          if (used_sha_id /= sha_id) then
             write(*,*) "************************************************"
             write(*,*) "WARNING: MCFOST has been updated"
             write(*,*) "- initial model was computed with mcfost:"
             write(*,*) "id="//used_sha_id
             write(*,*) "- this is mcfost:"
             write(*,*) "id="//sha_id
             write(*,*) "Problems can occur"
             write(*,*) "************************************************"
          endif
       endif

       !---------------------------------------------------------
       ! Memory allocation
       !---------------------------------------------------------
       write(*,*) "ProDiMo2mcfost version:", ProDiMo2mcfost_version
       if (ProDiMo2mcfost_version <= 2) then
            nLevel_oH2O = 10
            nLevel_pH2O = 9
       else
          nLevel_oH2O = 30
          nLevel_pH2O = 30
       endif
       lpops%lmax(1) = nlevel_CII
       lpops%lmax(2) = nlevel_OI
       lpops%lmax(3) = nlevel_CO
       lpops%lmax(4) = nlevel_oH2O
       lpops%lmax(5) = nlevel_pH2O

       allocate(nCII(n_rad,nz), dvCII(n_rad,nz), nOI(n_rad,nz), dvOI(n_rad,nz),  &
            noH2O(n_rad,nz), dvoH2O(n_rad,nz), npH2O(n_rad,nz), dvpH2O(n_rad,nz), &
            nCO(n_rad,nz), dvCO(n_rad,nz), stat=alloc_status )
       if (alloc_status > 0) call error('Allocation error nCII')
       nCII = 0.0 ; dvCII = 0.0 ; nOI = 0.0 ; dvOI = 0.0 ; nCO=0.0 ;  dvCO = 0.0
       noH2O = 0.0 ;  dvoH2O = 0.0 ;  npH2O = 0.0 ;  dvpH2O = 0.0

       allocate(pop_CII(nLevel_CII, n_rad, nz), stat=alloc_status )
       allocate(pop_OI(nLevel_OI, n_rad, nz), stat=alloc_status )
       allocate(pop_CO(nLevel_CO, n_rad, nz), stat=alloc_status )
       allocate(pop_oH2O(nLevel_oH2O, n_rad, nz), stat=alloc_status )
       allocate(pop_pH2O(nLevel_pH2O,n_rad, nz) , stat=alloc_status )
       if (alloc_status > 0) call error('Allocation error pop_CII')
       pop_CII = 0.0 ; pop_OI = 0.0 ; pop_oH2O = 0.0 ; pop_pH2O = 0.0 ; pop_CO=0.0 ;

       ! read_image : HDU 1 : MCgrid, to make a test inside mcfost
       call ftgpve(unit,group,firstpix,npixels,nullval,grid,anynull,fits_status)
       !   write(*,*) "Status1 = ", status

       ! Verification grille
       do ri=1,n_rad
          do zj=1,nz
             icell = cell_map(ri,zj,1)
             if (abs(r_grid(icell) - grid(ri,zj,1)) > 1e-5 * r_grid(icell)) then
                write(*,*) ri, "MCFOST=", r_grid(icell), "ProDiMo=", grid(ri,zj,1)
                call error("Non matching radius")
             endif
             if (abs(z_grid(icell) - grid(ri,zj,2)) > 1e-5 * z_grid(icell)) then
                write(*,*) ri, zj, "MCFOST=", z_grid(icell), "ProDiMo=", grid(ri,zj,2)
                call error("Non matching altitude")
             endif
          enddo
       enddo

       !---------------------------------------------------------
       ! HDU 2 : gas temperature [K]  64 bits
       !---------------------------------------------------------
       !  move to next hdu
       call ftmrhd(unit,1,hdutype,fits_status)

       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,10,naxes,nfound,fits_status)
       if (nfound /= 2) call error("failed to read the NAXISn keywords of HDU 2")
       if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) call error("HDU 2 does not have the right dimensions.")
       npixels=naxes(1)*naxes(2)

       ! read_image
       call ftgpvd(unit,group,firstpix,npixels,nullval,Tgas,anynull,fits_status)

       !---------------------------------------------------------
       ! HDU 3 : molecular particle densities [1/cm3]  64 bits
       !---------------------------------------------------------
       !  move to next hdu
       call ftmrhd(unit,1,hdutype,fits_status)

       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,10,naxes,nfound,fits_status)
       if (nfound /= 3)  call error("failed to read the NAXISn keywords of HDU 3")
       if ((naxes(1) /= MC_NSP).or.(naxes(2) /= n_rad).or.(naxes(3) /= nz)) call error("HDU 3 does not have the right dimensions.")
       npixels=naxes(1)*naxes(2)*naxes(3)

       ! read_image
       call ftgpvd(unit,group,firstpix,npixels,nullval,TMP,anynull,fits_status)

       ! Conversion en part.m^-3
       nCII = TMP(1,:,:) * m_to_cm**3
       nOI = TMP(2,:,:) * m_to_cm**3
       nCO = TMP(3,:,:) * m_to_cm**3
       noH2O = TMP(4,:,:) * m_to_cm**3
       npH2O = TMP(5,:,:) * m_to_cm**3

       !------------------------------------------------------------------------------
       ! HDU 4 : line broadening parameter [km/s]   64bits
       !------------------------------------------------------------------------------
       !  move to next hdu
       call ftmrhd(unit,1,hdutype,fits_status)

       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,10,naxes,nfound,fits_status)
       if (nfound /= 3) call error("failed to read the NAXISn keywords of HDU 4")
       if ((naxes(1) /= MC_NSP).or.(naxes(2) /= n_rad).or.(naxes(3) /= nz)) call error("HDU 4 does not have the right dimensions.")
       npixels=naxes(1)*naxes(2)*naxes(3)

       ! read_image
       call ftgpvd(unit,group,firstpix,npixels,nullval,TMP,anynull,fits_status)

       ! Conversion en m.s-1
       dvCII = TMP(1,:,:) * km_to_m
       dvOI = TMP(2,:,:) * km_to_m
       dvCO = TMP(3,:,:) * km_to_m
       dvoH2O = TMP(4,:,:) * km_to_m
       dvpH2O = TMP(5,:,:) * km_to_m

       !------------------------------------------------------------------------------
       ! HDU 5 ... 4+MC_NSP : relative population numbers
       !------------------------------------------------------------------------------
       do i=1,MC_NSP
          !  move to next hdu
          call ftmrhd(unit,1,hdutype,fits_status)

          ! Check dimensions
          call ftgknj(unit,'NAXIS',1,10,naxes,nfound,fits_status)
          if (nfound /= 3) call error("failed to read the NAXISn keywords of HDU 5+")
          if ((naxes(1) /= lpops%lmax(i)).or.(naxes(2) /= n_rad).or.(naxes(3) /= nz)) call error("HDU 5+ does not have the right dimensions.")
          npixels=naxes(1)*naxes(2)*naxes(3)

          allocate(MCpops(lpops%lmax(i),n_rad,nz), stat=alloc_status)
          if (alloc_status > 0) call error('Allocation error MCpops')
          MCpops = 0.0

          ! read_image
          call ftgpve(unit,group,firstpix,npixels,nullval,MCpops,anynull,fits_status)

          if (ProDiMo2mcfost_version > 1) then
             ! Les pops sont definies en relatif par rapport au niveau precedent
             sum_pops(:,:) = MCpops(1,:,:)
             do l=2, lpops%lmax(i)
                MCpops(l,:,:) = MCpops(l,:,:) * MCpops(l-1,:,:)
                sum_pops(:,:) = sum_pops(:,:) + MCpops(l,:,:)
             enddo !l

             ! On renormalise pour avoir une somme des populations a 1
             do l=1, lpops%lmax(i)
                MCpops(l,:,:) = MCpops(l,:,:) / sum_pops(:,:)
             enddo !l
          endif ! interface_version

          if (i==1) pop_CII = MCpops
          if (i==2) pop_OI = MCpops
          if (i==3) pop_CO = MCpops
          if (i==4) pop_oH2O = MCpops
          if (i==5) pop_pH2O = MCpops
          deallocate(MCpops)

       enddo !i


       !  Close the file and free the unit number.
       if (fits_status > 0) write(*,*) "Status1=", fits_status ! renvoie non 0 sauf si on l'imprime (????)
       call ftclos(unit, fits_status)
       if (fits_status > 0) write(*,*) "Status2=", fits_status ! renvoie non 0 sauf si on l'imprime (????)
       call ftfiou(unit, fits_status)
       !call ftfiou(-1, fits_status) ! deallocate toutes les unites fitsio, semble resoudre le probleme
       if (fits_status > 0) write(*,*) "Status3=", fits_status ! renvoie non 0 sauf si on l'imprime (????)

       !  Check for any error, and if so print out error messages
       !  Get the text string which describes the error
       if (fits_status > 0) then
          call ftgerr(fits_status,errtext)
          write(*,*) 'FITSIO Error Status =',fits_status,': ',errtext

          !  Read and print out all the error messages on the FITSIO stack
          call ftgmsg(errmessage)
          do while (errmessage .ne. ' ')
             write(*,*) errmessage
             call ftgmsg(errmessage)
          end do
       endif

       ! Setting the kinetic temperature to Prodimo's Tgas
       do i=1, n_rad
          do j=1, nz
             icell =cell_map(i,j,1)
             Tcin(icell) = Tgas(i,j)
          enddo
       enddo
    endif ! l_first_time

    ! Niveaux et populations
    write(*,*) "Setting ProDiMo abundances, population levels and Tgas"
    tab_abundance(:) = 0.0
    tab_nLevel(:,:) = 0.0
    do i=1, n_rad
       do j=1, nz
          icell =cell_map(i,j,1)
          if (lCII) then
             tab_nLevel(icell,1:nLevel_CII) = pop_CII(:,i,j) * nCII(i,j)
             tab_abundance(icell) = nCII(i,j)
          endif
          if (lOI) then
             tab_nLevel(icell,1:nLevel_OI) = pop_OI(:,i,j) * nOI(i,j)
             tab_abundance(icell) = nOI(i,j)
          endif
          if (lCO) then
             tab_nLevel(icell,1:nLevel_CO) = pop_CO(:,i,j) * nCO(i,j)  * factor_CO
             tab_abundance(icell) = nCO(i,j)  * factor_CO
          endif
          if (loH2O) then
             tab_nLevel(icell,1:nLevel_oH2O) = pop_oH2O(:,i,j) * noH2O(i,j)
             tab_abundance(icell) = noH2O(i,j)
          endif
          if (lpH2O) then
             tab_nLevel(icell,1:nLevel_pH2O) = pop_pH2O(:,i,j) * npH2O(i,j)
             tab_abundance(icell) = npH2O(i,j)
          endif
       enddo
    enddo
    do icell=1, n_cells
       tab_abundance(icell) = tab_abundance(icell) / densite_gaz(icell) ! conversion nbre en abondance
    enddo
    write(*,*) "Max =", maxval(tab_abundance), "min =", minval(tab_abundance)



!    iTrans = 2 ;
!    iUp = iTransUpper(iTrans)
!       iLow = iTransLower(iTrans)
!       cst = - hp * Transfreq(iTrans) / kb
!
!       write(*,*) "iTrans", iTrans, iLow, iUp
!       write(*,*), "g", poids_stat_g(iLow), poids_stat_g(iUp)
!       do j=1, nz
!          do i=1, n_rad
!             nUp = tab_nLevel(i,j,iUp)
!             nLow =  tab_nLevel(i,j,iLow)
!             if ((nUp > tiny_real) .and. (nLow > tiny_real) ) then
!                Tex = cst / log(  (nUp * poids_stat_g(iLow))  / (nLow * poids_stat_g(iUp) ))
!
!                write(*,*) "Cell", i, j, nUp, nLow
!                write(*,*) "ratio", nUp/nLow, poids_stat_g(iUp)/poids_stat_g(iLow) * exp(- hp * Transfreq(iTrans)/ (kb*Tgas(i,j)))
!                write(*,*) Tex, Tgas(i,j)
!                read(*,*)
!             endif
!          enddo ! i
!       enddo !j
!    !enddo !iTrans
!
!    write(*,*) "test done"
!    stop
!
!    Tcin(:,:,1) = Tgas
!    do i=1, n_rad
!       do j=1, nz
!          tab_nLevel(i,j,1) = 1.0
!          do l=2, nLevels
!             tab_nLevel(i,j,l) = tab_nLevel(i,j,l-1) * poids_stat_g(l)/poids_stat_g(l-1) * &
!                  exp(- hp * Transfreq(l-1)/ (kb*Tcin(i,j,1)))
!          enddo
!
!          Somme = sum(tab_nLevel(i,l,:))
!          if (Somme < 1e-30) then
!             write(*,*) "ERROR"
!             stop
!          endif
!          tab_nLevel(i,j,:) = tab_nLevel(i,j,:) * nCO(i,j) / Somme
!       enddo
!    enddo


    ! Vitesse keplerienne
    do i=1, n_rad
       do j=1, nz
          icell = cell_map(i,j,1)
          r_sph = sqrt(r_grid(icell)**2 + z_grid(icell)**2)
          ! vfield(i,j) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(i,j) * AU_to_m) ) ! Midplane Keplerian velocity
          vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg  * (r_grid(icell) * AU_to_m)**2 /  (r_sph * AU_to_m)**3 )
       enddo
    enddo

    ! Vitesse Doppler
    do i=1, n_rad
       do j=1, nz
          if (lCII) sigma2 =  dvCII(i,j)**2 !/ (2 * log(2.))
          if (lOI) sigma2 =  dvOI(i,j)**2
          if (lCO) sigma2 =  dvCO(i,j)**2
          if (loH2O) sigma2 =  dvoH2O(i,j)**2
          if (lpH2O) sigma2 =  dvpH2O(i,j)**2
          icell = cell_map(i,j,1)
          v_line(icell) = sqrt(sigma2)

          !  write(*,*) "FWHM", sqrt(sigma2 * log(2.)) * 2.  ! Teste OK bench water 1
          if (sigma2 <=0.) call error("ProDiMo data, dv = 0")

          sigma2_m1 = 1.0_dp / sigma2
          sigma2_phiProf_m1(icell) = sigma2_m1
          ! phi(nu) et non pas phi(v) donc facteur c_light et il manque 1/f0
          ! ATTENTION : il ne faut pas oublier de diviser par la freq apres
          norme_phiProf_m1(icell) = c_light / sqrt(pi * sigma2)
       enddo !i
    enddo !j

    write(*,*) "Done"
    write(*,*) " "

    return

  end subroutine read_ProDiMo2mcfost

  !********************************************************************

end module ProDiMo
