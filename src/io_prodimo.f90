module ProDiMo

  use parametres
  use opacity
  use em_th
  use constantes
  use molecular_emission
  use output, only: print_error
  use disk
  use prop_star
  use resultats
  use utils
  use ray_tracing, only:  RT_imin, RT_imax, RT_n_ibin, lRT_i_centered 
  use read_params
  use sha

  implicit none

  save

  character(len=512), parameter :: mcfost2ProDiMo_file = "forProDiMo.fits.gz"
  integer :: mcfost2ProDiMo_version, ProDiMo2mcfost_version

  ! directory with *.in files and output dir
  character(len=512) :: ProDiMo_input_dir, data_ProDiMo  

  ! Pour champ UV
  !real, parameter ::  slope_UV_ProDiMo = 2.2 - 2 ! Fnu -> F_lambda
  real :: fUV_ProDiMo, slope_UV_ProDiMo, ProDiMo_fPAH
  ! fPAH = (2.2/mPAH) * dust_to_gas * eps_MCFOST/3e-7 
  ! + PAH_NC, PAH_NH + distance + inclinaison
  character(len=10) :: sProDiMo_fPAH
  
  ! Grille de longeurs d'onde
  character(len=32) :: ProDiMo_tab_wavelength = "ProDiMo.lambda"

  ! Pour champ ISM
  real, parameter :: Wdil =  9.85357e-17 
  real, parameter :: TCmb = 2.73 
  real, parameter :: T_ISM_stars = 20000.
  real, parameter :: chi_ISM = 1.0
  real(kind=db), dimension(:,:,:), allocatable :: J_ProDiMo, N_ProDiMo

  real(kind=db), dimension(:,:), allocatable :: n_phot_envoyes_ISM
        

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
  integer, parameter :: nLevel_oH2O = 10
  integer, parameter :: nLevel_pH2O = 9

  real, dimension(:,:), allocatable :: nCII, dvCII, nOI, dvOI,  noH2O, dvoH2O, npH2O, dvpH2O, nCO, dvCO
  real, dimension(:,:,:), allocatable :: pop_CII
  real, dimension(:,:,:), allocatable :: pop_OI
  real, dimension(:,:,:), allocatable :: pop_CO
  real, dimension(:,:,:), allocatable :: pop_oH2O
  real, dimension(:,:,:), allocatable :: pop_pH2O


contains

  subroutine setup_ProDiMo()
    ! - distance 100pc
    ! - 1 pop de grains loi de puissance

    ! lambda : 13 bins entre 0.0912 et 3410.85 microns
    ! donne les 11 bins de ProDiMo

    ! TODO : linit_gaz

    real :: wl_min
    integer :: alloc_status

    ! Directories
    if (.not.lprodimo_input_dir) then
       ProDiMo_input_dir = trim(mcfost_utils)//"/forProDiMo"
       write(*,*) "Using "//trim(ProDiMo_input_dir)//" for ProDiMo input"
    endif
    data_ProDiMo = trim(root_dir)//"/"//trim(seed_dir)//"/data_ProDiMo"
    
    if (lforce_ProDiMo_PAH) data_ProDiMo = trim(data_ProDiMO)//"_fPAH="//sProDiMo_fPAH
    
    ! Limite 10x plus haut pour avoir temperature plus propre pour ProDiMo
    tau_dark_zone_eq_th = 15000. 

    ! TODO : verif fichier parametres : pour relecture avec Yorick
    ! N_etoiles = 1
    if (mcfost2ProDiMo_version == 1) then
       write(*,*) "***************************"
       write(*,*) "* Modelling for ProDiMo   *"
       write(*,*) "* Forcing wavelength grid *"
       write(*,*) "***************************"
       n_lambda = 39
       lambda_min = 0.0912
       lambda_max = 3410.85

       ! on veut changer les lambda pour le step 2
       lsed_complete = .true.
    endif

    fUV_ProDiMo = etoile(1)%fUV
    slope_UV_ProDiMo = etoile(1)%slope_UV

    allocate(xN_abs(n_lambda2,n_rad,nz,nb_proc),  stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error xN_abs'
       stop
    endif
    xN_abs = 0.0

    allocate(J_ProDiMo(n_lambda2,n_rad,nz),  stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error J_ProDiMo'
       stop
    endif
    J_ProDiMo = 0.0

    allocate(N_ProDiMo(n_lambda2,n_rad,nz),  stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error J_ProDiMo'
       stop
    endif
    N_ProDiMo = 0.0

    allocate(n_phot_envoyes_ISM(n_lambda2,nb_proc),  stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_phot_envoyes_ISM'
       stop
    endif
    n_phot_envoyes_ISM = 0.0

    return

  end subroutine setup_ProDiMo


  !********************************************************************
  
  subroutine save_J_prodimo(lambda)
    ! C. Pinte
    ! 19/06/09
    ! sauvegarde le champ de radiation pour ProDiMo
    ! avant de calculer le champ ISM

    integer, intent(in) :: lambda
    integer :: ri, zj
    real (kind=db) :: n_photons_envoyes, energie_photon, facteur

    ! Step1
    ! 1/4pi est inclus dans n_phot_l_tot
    !    J_ProDiMo(:,:,:) = (sum(xJ_abs,dim=4) + J0(:,:,:,1)) * n_phot_L_tot !* 4* pi 
    !
    !    do ri=1, n_rad
    !       J_ProDiMo(:,ri,:) = J_ProDiMo(:,ri,:) / volume(ri)
    !    enddo
    !
    !    ! xJ_abs est par bin de lambda donc Delta_lambda.F_lambda
    !    ! J en W.m-2 (lambda.F_lambda)
    !    ! teste OK par rapport a fct bb de yorick
    !    do lambda=1, n_lambda
    !       J_ProDiMo(lambda,:,:) = J_ProDiMo(lambda,:,:) * tab_lambda(lambda) / tab_delta_lambda(lambda)
    !    enddo

    ! Step2
    n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
    energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda)) / n_photons_envoyes &
         * tab_lambda(lambda) * 1.0e-6  !lambda.F_lambda

    do ri=1, n_rad
       do zj=1,nz
          facteur = energie_photon / volume(ri) 
          J_prodimo(lambda,ri,zj) =  J_prodimo(lambda,ri,zj) +  facteur * sum(xJ_abs(lambda,ri,zj,:))
          N_ProDiMo(lambda,ri,zj) =  sum(xN_abs(lambda,ri,zj,:))
       enddo
    enddo
        
    xJ_abs(lambda,:,:,:) = 0.0
    xN_abs(lambda,:,:,:) = 0.0

    return

  end subroutine save_J_prodimo

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

    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(5) :: naxes
    integer :: group,fpixel,nelements, alloc_status, id, lambda, ri, zj, l, i
    real (kind=db) :: n_photons_envoyes, energie_photon, facteur, N
    real :: wl, norme

    logical :: simple, extend
    character(len=512) :: filename

    real(kind=db), dimension(n_rad,nz,2) :: grid
    real, dimension(n_rad,nz) :: dens
    real, dimension(n_rad,nz,0:3) :: N_grains
    real, dimension(n_lambda) :: spectre

    ! Allocation dynamique pour eviter d'utiliser la memeoire stack
    real(kind=db), dimension(:,:,:), allocatable :: J_ProDiMo_ISM    ! n_lambda,n_rad,nz)
    real, dimension(:,:,:), allocatable :: J_io ! n_rad,nz,n_lambda, joue aussi le role de N_io
    real, dimension(:,:,:,:), allocatable :: opacite ! (n_rad,nz,2,n_lambda)


    call create_input_ProDiMo()

    allocate(opacite(n_rad,nz,2,n_lambda), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error opacite forProDiMo.fits.gz'
       stop
    endif
    opacite = 0
   
    allocate(J_io(n_rad,nz,n_lambda), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error J_io forProDiMo.fits.gz'
       stop
    endif
    J_io = 0    

    allocate(J_ProDiMo_ISM(n_lambda,n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then
       write(*,*) 'Allocation error J_ProDiMo_ISM forProDiMo.fits.gz'
       stop
    endif
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

    call ftpkye(unit,'Teff',etoile(1)%T,-8,'[K]',status)
    call ftpkye(unit,'Rstar',real(etoile(1)%r*AU_to_Rsun),-8,'[Rsun]',status)
    call ftpkye(unit,'Mstar',real(etoile(1)%M),-8,'[Msun]',status)
    call ftpkye(unit,'fUV',fUV_ProDiMo,-3,'',status)
    call ftpkye(unit,'slope_UV',slope_UV_ProDiMo+2.0,-3,'',status) ! on remet le 2 
    call ftpkye(unit,'distance',real(distance),-8,'[pc]',status)

    call ftpkye(unit,'disk_dust_mass',real(diskmass),-8,'[Msun]',status)
    call ftpkye(unit,'Rin',real(disk_zone(1)%rin),-8,'[AU]',status)
    call ftpkye(unit,'Rout',real(disk_zone(1)%rout),-8,'[AU]',status)
    call ftpkye(unit,'Rref',real(disk_zone(1)%rref),-8,'[AU]',status)
    call ftpkye(unit,'H0',real(disk_zone(1)%sclht),-8,'[AU]',status)
    call ftpkye(unit,'edge',real(disk_zone(1)%edge),-8,'[AU]',status)
    call ftpkye(unit,'beta',real(disk_zone(1)%exp_beta),-8,'',status)
    call ftpkye(unit,'alpha',real(disk_zone(1)%surf),-8,'',status)

    call ftpkye(unit,'amin',dust_pop(1)%amin,-8,'[micron]',status)
    call ftpkye(unit,'amax',dust_pop(1)%amax,-8,'[micron]',status)
    call ftpkye(unit,'aexp',dust_pop(1)%aexp,-8,'slope of grain size distribution',status)
    call ftpkye(unit,'strat',exp_strat,-8,'stratification exponent',status)
    call ftpkye(unit,'a_settle',a_strat,-8,'[micron]',status)
    call ftpkye(unit,'rho_grain',dust_pop(1)%component_rho1g(1),-8,'[g.cm^-3]',status)
    call ftpkys(unit,'optical_indices',trim(dust_pop(1)%indices(1)),'',status)

    call ftpkyj(unit,'n_grains',dust_pop(1)%n_grains,' ',status)
    call ftpkyj(unit,'n_rad',n_rad,' ',status)
    call ftpkyj(unit,'nz',nz,' ',status)
    call ftpkyj(unit,'n_rad_in',n_rad_in,' ',status)

    !  Write the array to the FITS file.
    do ri=1, n_rad
       grid(ri,:,1) = r_grid(ri,:)
       do zj=1,nz
          grid(ri,zj,2) = z_grid(ri,zj)
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
    if (lRE_nLTE) then
       do ri=1, n_rad
          do zj=1, nz
             Temperature(ri,zj,1) = sum( Temperature_1grain(ri,zj,:) * r_grain(:)**2 * nbre_grains(:)) / norme
          enddo !j
       enddo !i
       
       write(*,*) "************************************************************"
       write(*,*) "WARNING: nLTE mode in MCFOST"
       write(*,*) "Crude averaging of the temperature distribution for ProDiMo"
       write(*,*) "T = \int T(a) * n(a) * a^2 da / \int n(a) * a^2 da"
       write(*,*) "Max. averaged temperature = ", maxval(Temperature(:,:,1))
       write(*,*) "************************************************************"

    endif


    !  Write the array to the FITS file.
    call ftppre(unit,group,fpixel,nelements,Temperature(:,:,1),status)

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
    call ftppre(unit,group,fpixel,nelements,tab_lambda,status)

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
           * (4.*pi*(R_ISM*size_neb)**2) / n_photons_envoyes / pi
              
       do ri=1, n_rad
          do zj=1,nz
             facteur = energie_photon / volume(ri) 
             J_prodimo_ISM(lambda,ri,zj) =  J_prodimo_ISM(lambda,ri,zj) +  facteur * sum(xJ_abs(lambda,ri,zj,:))
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
          J_io(ri,zj,:) = sum(xN_abs(:,ri,zj,:), dim=2) 
       enddo
    enddo

    call ftppre(unit,group,fpixel,nelements,J_io,status)


    !------------------------------------------------------------------------------
    ! HDU 10 : Densite de gaz pour un rapport de masse de 100 / poussiere
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
          dens(ri,zj) =  densite_gaz(ri,zj,1) * masse_mol_gaz / m_to_cm**3 ! g.cm^-3
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

    do zj=1,nz
       do ri=1,n_rad
          opacite(ri,zj,1,:) = kappa(:,ri,zj,1)
          if (lRE_LTE) then 
             opacite(ri,zj,2,:) = kappa_abs_eg(:,ri,zj,1)
          else if (lRE_nLTE) then ! Patch, c'est pas super propre
             do lambda=1,n_lambda 
                do l=grain_RE_nLTE_start,grain_RE_nLTE_end
                   opacite(ri,zj,2,lambda) = opacite(ri,zj,2,lambda) + q_abs(lambda,l) * densite_pouss(ri,zj,1,l)

                enddo ! k
             enddo ! lambda
          endif
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
          N = sum(densite_pouss(ri,zj,1,:))
          N_grains(ri,zj,0) = N
          if (N > 0) then
             N_grains(ri,zj,1) = sum(densite_pouss(ri,zj,1,:) * r_grain(:)) / N
             N_grains(ri,zj,2) = sum(densite_pouss(ri,zj,1,:) * r_grain(:)**2) / N
             N_grains(ri,zj,3) = sum(densite_pouss(ri,zj,1,:) * r_grain(:)**3) / N
          else 
             N_grains(ri,zj,1) = 0.0
             N_grains(ri,zj,2) = 0.0
             N_grains(ri,zj,3) = 0.0
          endif
       enddo
    enddo

    ! part.cm^-3 --> part.m^-3
    N_grains(:,:,0) = N_grains(:,:,0) /  (cm_to_m**3)

!    write(*,*)  "TEST_0", sum(N_grains(1,:,0))
!    write(*,*)  "TEST_3", sum(N_grains(1,:,3) * N_grains(1,:,0))

    call ftppre(unit,group,fpixel,nelements,N_grains,status)
    
    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

    !  Check for any error, and if so print out error messages
    if (status > 0) then
       call print_error(status)
       stop
    end if

    return

  end subroutine mcfost2ProDiMo

  !********************************************************************

  subroutine create_input_ProDiMo()
    ! Create the *.in files for ProDiMo
    ! C. Pinte 17/03/2012

    character(len=512) :: cmd
    character(len=128) :: line
    character(len=10) :: fmt
    integer :: j, status, syst_status
    logical :: unchanged

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
       if (INDEX(line,"! dust_to_gas") > 0) then
          write(2,*) 1.0/gas_dust," ! dust_to_gas : set by MCFOST"
          unchanged = .false.
       endif

       if (INDEX(line,"! fPAH") > 0) then
          write(2,*) ProDiMo_fPAH," ! fPAH : set by MCFOST"
          unchanged = .false.
       endif

       if (INDEX(line,"! Rout") > 0) then
          write(2,*) Rout," ! Rout : set by MCFOST"
          unchanged = .false.
       endif
       
       if (unchanged) then
          write(fmt,'("(A",I3.3,")")') len(trim(line)) 
          write(2,fmt) trim(line)
       endif
    enddo read_loop

    close(1)
    close(2)

    return

  end subroutine create_input_ProDiMo

  !********************************************************************

  subroutine read_mcfost2ProDiMo()
    ! Relit le fichier for ProDiMo.fits.gz cree par mcfost pour ProDiMo
    ! afin de redemarrer dessus pour le transfert dans les raies
    ! Lit les parametres et la structure en temperature
    ! A terme le champ de radiation pour calculer le contribution de 
    ! lumiere diffusee dans les raies
    !
    ! C.Pinte
    ! 12/07/09

    integer :: status, readwrite, unit, blocksize, nfound, group, firstpix, nbuffer, npixels, hdunum, hdutype, imol, syst_status, n_files
    real :: nullval, buffer
    logical :: anynull

    integer, dimension(4) :: naxes

    character(len=30) :: errtext
    character (len=80) :: errmessage, comment
    character(len=512) :: filename, mol_para, cmd


    real, dimension(:,:), allocatable :: temperature_ProDiMo


    !*********************************************************
    ! Lecture des parametres dans le fichier .para de data_th
    !*********************************************************

    ! Sauvegarde du nom du fichier de parametres pour les raies
    mol_para = para
    
    ! Copie temporaire et lecture du fichier de parametres
    ! car je ne connais pas le nom du fichier .par dans data_th
    n_files = 0
    cmd = "ls data_th/*.par* | wc -l > n_files.tmp" ; call appel_syst(cmd, syst_status)
    open(unit=1, file="n_files.tmp", status='old')
    read(1,*) n_files
    close(unit=1,status="delete")
    
    if (n_files > 1) then
       write(*,*) "There are more than 1 parameter file in data_th"
       write(*,*) "Exiting"
       stop
    else if (n_files < 1) then
       write(*,*) "There are less than 1 parameter file in data_th"
       write(*,*) "Exiting"
       stop
    endif
    
    cmd = "cp data_th/*.par* data_th/forMCFOST.par" ; call appel_syst(cmd, syst_status)
    para = "data_th/forMCFOST.par" ; call read_para()
    cmd = "rm -rf data_th/forMCFOST.par" ; call appel_syst(cmd, syst_status)

    ! Restaure nom fichier parametres
    para = mol_para

    ! Parametres par defaut
    ltemp = .false.
    lsed = .false.
    lsed_complete = .false.
    lemission_mol = .true.
    lsetup_gas=.true.
    lfits = .false.

    ! on ne force plus depuis version 2 de l'interface
    !n_lambda = 39 
    !lambda_min = 0.0912
    !lambda_max = 3410.85

    deallocate(mol) ! car on va le remplacer

    
!---    !**********************************************
!---    ! Lecture des parametres dans le fichier fits
!---    !**********************************************
!---    
!---    filename = "data_th/"//trim(mcfost2ProDiMo_file)
!---
!---
!---    ! Parametres par defaut
!---    lcheckpoint=.false.
!---    ltemp = .false.
!---    lsed = .false.
!---    lemission_mol = .true.
!---    lsetup_gas=.true.
!---    lfits = .false.
!---    n_lambda = 39
!---    lambda_min = 0.0912
!---    lambda_max = 3410.85
!---    grid_type=1
!---    n_az =1
!---    angle_interet = 45. ! TODO
!---    lcavity = .false.
!---    lRE_LTE = .true.
!---    lRE_nLTE = .false.
!---    lnRE = .false.
!---
!---    l_sym_ima=.true.
!---    l_sym_centrale=.true.
!---    l_sym_axiale=.true.
!---
!---    n_etoiles = 1
!---    n_zones = 1
!---    n_pop = 1
!---    allocate(etoile(n_etoiles))
!---    allocate(disk_zone(n_zones))
!---    allocate(dust_pop(n_pop))
!---    dust_pop(1)%zone = 1
!---    disk_zone(1)%geometry = 1
!---
!---    !  Get an unused Logical Unit Number to use to open the FITS file.
!---    status=0
!---    call ftgiou (unit,status)
!---
!---    write(*,*) "Reading "//trim(filename)
!---
!---    ! Open file
!---    readwrite=0
!---    call ftopen(unit,filename,readwrite,blocksize,status)
!---    if (status /= 0) then ! le fichier temperature n'existe pas
!---       write(*,*) "ERROR : "//trim(mcfost2ProDiMo_file)//" file needed"
!---       stop
!---    endif
!---
!---    !------------------------------------------------------------------------------
!---    ! Spatial grid : pas besoin de relire
!---    !------------------------------------------------------------------------------
!---    ! hdu 1
!---
!---    ! Keywords
!---    ! FTGKY[EDJKLS](unit,keyword, > keyval,comment,status)
!---    call ftgkye(unit,'Teff',etoile(1)%T,comment,status)
!---    call ftgkye(unit,'Rstar',buffer,comment,status)
!---    etoile(1)%R = buffer/AU_to_Rsun
!---    call ftgkye(unit,'Mstar',buffer,comment,status)
!---    etoile(1)%M = buffer
!---    call ftgkye(unit,'fUV',fUV_ProDiMo,comment,status)
!---    call ftgkye(unit,'distance',buffer,comment,status)
!---    distance = buffer
!---
!---    call ftgkye(unit,'disk_dust_mass',buffer,comment,status)
!---    disk_zone(1)%diskmass = buffer    
!---    diskmass = buffer
!---    call ftgkye(unit,'Rin',buffer,comment,status)
!---    disk_zone(1)%rin = buffer
!---    call ftgkye(unit,'Rout',buffer,comment,status)
!---    disk_zone(1)%rout = buffer
!---    call ftgkye(unit,'Rref',buffer,comment,status)
!---    disk_zone(1)%rref = buffer
!---    call ftgkye(unit,'H0',buffer,comment,status)
!---    disk_zone(1)%sclht = buffer
!---    call ftgkye(unit,'edge',buffer,comment,status)
!---    disk_zone(1)%edge = buffer
!---    call ftgkye(unit,'beta',buffer,comment,status)
!---    disk_zone(1)%exp_beta = buffer
!---    call ftgkye(unit,'alpha',buffer,comment,status)
!---    disk_zone(1)%surf = buffer
!---
!---    rmin = disk_zone(1)%rin - 5 * disk_zone(1)%edge
!---    grid_rmin = rmin
!---    rout = disk_zone(1)%rout
!---    
!---    call ftgkye(unit,'amin',dust_pop(1)%amin,comment,status)
!---    call ftgkye(unit,'amax',dust_pop(1)%amax,comment,status)
!---    call ftgkye(unit,'aexp',dust_pop(1)%aexp,comment,status)
!---    call ftgkye(unit,'strat',exp_strat,comment,status)
!---    call ftgkye(unit,'rho_grain',dust_pop(1)%rho1g,comment,status)
!---    call ftgkys(unit,'optical_indices',dust_pop(1)%indices,comment,status)
!---
!---    call ftgkyj(unit,'n_grains',dust_pop(1)%n_grains,comment,status)
!---    call ftgkyj(unit,'n_rad',n_rad,comment,status)
!---    call ftgkyj(unit,'nz',nz,comment,status)
!---    call ftgkyj(unit,'n_rad_in',n_rad_in,comment,status)
!---
!---    dust_pop(1)%ind_debut = 1
!---    dust_pop(1)%ind_fin = dust_pop(1)%n_grains
!---    n_grains_RE_LTE = dust_pop(1)%n_grains
!---    n_grains_RE_nLTE = 0
!---    n_grains_nRE = 0
!---    n_grains_tot = dust_pop(1)%n_grains
!---    grain_RE_LTE_start = 1
!---    grain_RE_LTE_end = dust_pop(1)%n_grains
!---    dust_pop(1)%masse = disk_zone(1)%diskmass
!---
!---    !  Close the file and free the unit number.
!---    call ftclos(unit, status)
!---    call ftfiou(unit, status)
!---
!---    !  Check for any error, and if so print out error messages
!---    !  Get the text string which describes the error
!---    if (status > 0) then
!---       call ftgerr(status,errtext)
!---       write(*,*) "ERROR : "//trim(filename)
!---       print *,'FITSIO Error Status =',status,': ',errtext
!---
!---       !  Read and print out all the error messages on the FITSIO stack
!---       call ftgmsg(errmessage)
!---       do while (errmessage .ne. ' ')
!---          print *,errmessage
!---          call ftgmsg(errmessage)
!---       end do
!---    endif

    open(unit=1, file=para, status='old')
    
    ! Inclinaisons
    read(1,*)     
    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered 

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
       mol(imol)%n_speed = 1 ! inutilise si on ne calcule pas le NLTE
    

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

    return

  end subroutine read_mcfost2ProDiMo

  !********************************************************************
  
  subroutine read_ProDiMo2mcfost(imol)
    ! C. Pinte  17/09/09
    ! Lit les resultats de ProDiMo et met a jour les valeurs dans mcfost
    ! -->  Tgas, abondances a terme
    ! population et vitesses pour le moment
    
    integer, intent(in) :: imol

    integer :: fits_status, readwrite, unit, blocksize, nfound, group, firstpix, nbuffer, npixels, hdunum, hdutype,alloc_status
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

    real(kind=db), dimension(MC_NSP,n_rad,nz) :: TMP
    real(kind=db), dimension(n_rad,nz) :: Tgas ! Pas utilise pour le moment, pour futurs calculs NLTE 
    real, dimension(n_rad,nz,2) :: grid ! Seulement pout test
    real, dimension(n_rad,nz) :: sum_pops
    real, dimension(:,:,:), allocatable :: MCpops

    logical :: lCII, lOI, lCO, loH2O, lpH2O 

    real :: sigma2, vmax, sigma2_m1, r_cyl
    real(kind=db) :: S, dz, dV, Somme
    integer :: i,j, ri, zj, n1, n2, iv, n_speed_rt, n_speed, l, keyword_status


    ! TMP
    integer :: iTrans, iUp, iLow
    real(kind=db) :: cst, nUp, nLow
    real :: Tex

    n_speed_rt = mol(imol)%n_speed_rt

    lCII = .false. ; lOI = .false. ; lCO = .false. ; loH2O = .false. ;  lpH2O = .false. ; 
    if (mol(imol)%name=="C+") lCII = .true.
    if (mol(imol)%name=="O") lOI = .true.
    if (mol(imol)%name=="CO") lCO = .true.
    if (mol(imol)%name=="o-H2O") loH2O = .true.
    if (mol(imol)%name=="p-H2O") lpH2O = .true.

    lpops%lmax(1) = nlevel_CII
    lpops%lmax(2) = nlevel_OI
    lpops%lmax(3) = nlevel_CO
    lpops%lmax(4) = nlevel_oH2O
    lpops%lmax(5) = nlevel_pH2O

    ldust_mol = .true.


    if (l_first_time) then
       l_first_time = .false.

       allocate(nCII(n_rad,nz), dvCII(n_rad,nz), nOI(n_rad,nz), dvOI(n_rad,nz),  &
            noH2O(n_rad,nz), dvoH2O(n_rad,nz), npH2O(n_rad,nz), dvpH2O(n_rad,nz), &
            nCO(n_rad,nz), dvCO(n_rad,nz), stat=alloc_status )
       if (alloc_status > 0) then
          write(*,*) 'Allocation error nCII'
          stop
       endif
       nCII = 0.0 ; dvCII = 0.0 ; nOI = 0.0 ; dvOI = 0.0 ; nCO=0.0 ;  dvCO = 0.0
       noH2O = 0.0 ;  dvoH2O = 0.0 ;  npH2O = 0.0 ;  dvpH2O = 0.0 
       

       allocate(pop_CII(nLevel_CII, n_rad, nz), stat=alloc_status )
       allocate(pop_OI(nLevel_OI, n_rad, nz), stat=alloc_status )
       allocate(pop_CO(nLevel_CO, n_rad, nz), stat=alloc_status )
       allocate(pop_oH2O(nLevel_oH2O, n_rad, nz), stat=alloc_status )
       allocate(pop_pH2O(nLevel_pH2O,n_rad, nz) , stat=alloc_status )
       if (alloc_status > 0) then
          write(*,*) 'Allocation error pop_CII'
          stop
       endif
       pop_CII = 0.0 ; pop_OI = 0.0 ; pop_oH2O = 0.0 ; pop_pH2O = 0.0 ; pop_CO=0.0 ; 


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
       if (fits_status > 0) then ! le fichier temperature n'existe pas
          write(*,*) "Status0=", fits_status
          write(*,*) "ERROR : "//trim(filename)//" file needed"
          stop
       endif

       group=1
       firstpix=1
       nullval=-999

       !---------------------------------------------------------
       ! HDU 1 : MCgrid, to make a test inside mcfost
       ! Format is the same as for the mcfost2ProDiMo interface
       !---------------------------------------------------------
       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,3,naxes,nfound,fits_status)
       if (nfound /= 3) then
          write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
          write(*,*) 'of HDU 1 file. Exiting.'
          stop
       endif
       if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= 2)) then
          write(*,*) "Error : HDU 1 does not have the right dimensions. Exiting."
          write(*,*) naxes(1:3)
          stop
       endif
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


       
       ! read_image
       call ftgpve(unit,group,firstpix,npixels,nullval,grid,anynull,fits_status)
       !   write(*,*) "Status1 = ", status
    
       ! Verification grille
       do ri=1,n_rad
          do zj=1,nz
             if (abs(r_grid(ri,zj) - grid(ri,zj,1)) > 1e-5 * r_grid(ri,zj)) then 
                write(*,*) "ERROR R. Exiting."
                write(*,*) ri, "MCFOST=", r_grid(ri,zj), "ProDiMo=", grid(ri,zj,1)
                !stop
             endif
             
             if (abs(z_grid(ri,zj) - grid(ri,zj,2)) > 1e-5 * z_grid(ri,zj)) then 
                write(*,*) "ERROR Z. Exiting."
                write(*,*) ri, zj, z_grid(ri,zj), grid(ri,zj,2)
                !stop
             endif
          enddo
       enddo

       !---------------------------------------------------------
       ! HDU 2 : gas temperature [K]  64 bits
       !---------------------------------------------------------
       !  move to next hdu
       call ftmrhd(unit,1,hdutype,fits_status)    
       
       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,2,naxes,nfound,fits_status)
       if (nfound /= 2) then
          write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
          write(*,*) 'of HDU 2 file. Exiting.'
          stop
       endif
       if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
          write(*,*) "Error : HDU 2 does not have the"
          write(*,*) "right dimensions. Exiting."
          stop
       endif
       npixels=naxes(1)*naxes(2)

       ! read_image
       call ftgpvd(unit,group,firstpix,npixels,nullval,Tgas,anynull,fits_status)

       !---------------------------------------------------------
       ! HDU 3 : molecular particle densities [1/cm3]  64 bits
       !---------------------------------------------------------
       !  move to next hdu
       call ftmrhd(unit,1,hdutype,fits_status)    
       
       ! Check dimensions
       call ftgknj(unit,'NAXIS',1,3,naxes,nfound,fits_status)
       if (nfound /= 3) then
          write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
          write(*,*) 'of HDU 3 file. Exiting.'
          stop
       endif
       if ((naxes(1) /= MC_NSP).or.(naxes(2) /= n_rad).or.(naxes(3) /= nz)) then
          write(*,*) "Error : HDU 3 does not have the"
          write(*,*) "right dimensions. Exiting."
          write(*,*) naxes(1), n_lambda
          stop
       endif
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
       call ftgknj(unit,'NAXIS',1,3,naxes,nfound,fits_status)
       if (nfound /= 3) then
          write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
          write(*,*) 'of HDU 4 file. Exiting.'
          stop
       endif
       if ((naxes(1) /= MC_NSP).or.(naxes(2) /= n_rad).or.(naxes(3) /= nz)) then
          write(*,*) "Error : HDU 4 does not have the"
          write(*,*) "right dimensions. Exiting."
          stop
       endif
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
          call ftgknj(unit,'NAXIS',1,3,naxes,nfound,fits_status)
          if (nfound /= 3) then
             write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
             write(*,*) 'of HDU ',i+4,' file. Exiting.'
             stop
          endif
          if ((naxes(1) /= lpops%lmax(i)).or.(naxes(2) /= n_rad).or.(naxes(3) /= nz)) then
             write(*,*) "Error : HDU ",i+4," does not have the"
             write(*,*) "right dimensions. Exiting."
             write(*,*) naxes(1:3)
             write(*,*) lpops%lmax(i)
             stop
          endif
          npixels=naxes(1)*naxes(2)*naxes(3)

          allocate(MCpops(lpops%lmax(i),n_rad,nz), stat=alloc_status)
          if (alloc_status > 0) then
             write(*,*) 'Allocation error MCpops'
             stop
          endif
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
          print *,'FITSIO Error Status =',fits_status,': ',errtext

          !  Read and print out all the error messages on the FITSIO stack
          call ftgmsg(errmessage)
          do while (errmessage .ne. ' ')
             print *,errmessage
             call ftgmsg(errmessage)
          end do
       endif

    endif ! l_first_time

    ! Niveaux et populations
    write(*,*) "Setting ProDiMo abundances, population levels and Tgas"
    tab_abundance(:,:) = 0.0
    tab_nLevel(:,:,:) = 0.0
    do i=1, n_rad
       do j=1, nz    
          if (lCII) then
             tab_nLevel(i,j,1:nLevel_CII) = pop_CII(:,i,j) * nCII(i,j)
             tab_abundance(i,j) = nCII(i,j)
          endif
          if (lOI) then
             tab_nLevel(i,j,1:nLevel_OI) = pop_OI(:,i,j) * nOI(i,j)
             tab_abundance(i,j) = nOI(i,j)
          endif
          if (lCO) then
             tab_nLevel(i,j,1:nLevel_CO) = pop_CO(:,i,j) * nCO(i,j)
             tab_abundance(i,j) = nCO(i,j)
          endif
          if (loH2O) then
             tab_nLevel(i,j,1:nLevel_oH2O) = pop_oH2O(:,i,j) * noH2O(i,j) 
             tab_abundance(i,j) = noH2O(i,j)
          endif
          if (lpH2O) then
             tab_nLevel(i,j,1:nLevel_pH2O) = pop_pH2O(:,i,j) * npH2O(i,j)
             tab_abundance(i,j) = npH2O(i,j)
          endif
       enddo
    enddo
    tab_abundance = tab_abundance / densite_gaz(:,:,1) ! conversion nbre en abondance
    write(*,*) "Max =", maxval(tab_abundance), "min =", minval(tab_abundance)
    
    Tcin(:,:,1) = Tgas



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
          r_cyl = sqrt(r_grid(i,j)**2 + z_grid(i,j)**2)
!          vfield(i,j) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(i,j) * AU_to_m) )
          vfield(i,j) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg  * (r_grid(i,j) * AU_to_m)**2 /  (r_cyl * AU_to_m)**3 )
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
          v_line(i,j) = sqrt(sigma2)
             
          !  write(*,*) "FWHM", sqrt(sigma2 * log(2.)) * 2.  ! Teste OK bench water 1
          if (sigma2 <=0.) then
             write(*,*) "ERROR in ProDiMo data, dv = 0"
             stop
          endif

          sigma2_m1 = 1.0_db / sigma2
          sigma2_phiProf_m1(i,j) = sigma2_m1
          ! phi(nu) et non pas phi(v) donc facteur c_light et il manque 1/f0
          ! ATTENTION : il ne faut pas oublier de diviser par la freq apres
          norme_phiProf_m1(i,j) = c_light / sqrt(pi * sigma2) 
          
!----          ! Echantillonage du profil de vitesse dans la cellule
!----          ! 2.15 correspond a l'enfroit ou le profil de la raie faut 1/100 de 
!----          ! sa valeur au centre : exp(-2.15^2) = 0.01
!----          vmax = sqrt(sigma2)
!----          tab_dnu_o_freq(i,j) = largeur_profile * vmax / (real(n_speed) )
!----          do iv=-n_speed, n_speed
!----             tab_deltaV(iv,i,j) = largeur_profile * real(iv,kind=db)/real(n_speed,kind=db) * vmax
!----          enddo ! iv
!----          deltaVmax(i,j) = largeur_profile * vmax !* 2.0_db  ! facteur 2 pour tirage aleatoire
       enddo !i
    enddo !j

    write(*,*) "Done"
    write(*,*) " "
     
    return
    
  end subroutine read_ProDiMo2mcfost

  !********************************************************************

end module ProDiMo
