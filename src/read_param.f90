module read_params

  use parametres
  use grains
  use molecular_emission
  use wavelengths
  use constantes
  use sha
  use messages

  implicit none

contains

  subroutine read_para(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old', iostat=ios)
    if (ios/=0) call error("cannot open "//trim(para))

    read(1,*) para_version

    ! Version 2.21 and 3.0 are the same
    if (abs(para_version - 2.21) < 1.e-4) para_version = 3.0

    ! Petit test pour le ray-tracing
    if (lscatt_ray_tracing .and. (abs(para_version) < 2.0901)) then
       call error("Parameter version >= 2.10 required for ray-tracing.")
    endif

    ! Petit test pour cavite
    if (lcavity .and. (abs(para_version) < 3.0)) then
       call error("Parameter version >= 3.0 required for option -lcavity.", &
            msg2="Use cavity section is parameter file for earlier version.")
    endif

    correct_Rsub = 1.0_dp
    lmigration = .false.
    lhydrostatic = .false.
    lread_Misselt=.false.
    lread_DustEM=.false.

    if (abs(para_version - 3.0) > 1.e-4) then
       write(*,*) "Wrong version of the parameter file."
       if (abs(para_version-2.20) < 1.e-4) then
          write(*,*) "Trying to read 2.20 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para220(para)
          return
       else if (abs(para_version-2.19) < 1.e-4) then
          write(*,*) "Trying to read 2.19 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para219(para)
          return
    else if (abs(para_version-2.18) < 1.e-4) then
          write(*,*) "Trying to read 2.18 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para218(para)
          return
       else if (abs(para_version-2.17) < 1.e-4) then
          write(*,*) "Trying to read 2.17 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para217(para)
          return
       else if (abs(para_version-2.16) < 1.e-4) then
          write(*,*) "Trying to read 2.16 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para216(para)
          return
       else if (abs(para_version-2.15) < 1.e-4) then
          write(*,*) "Trying to read 2.15 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para215(para)
          return
       else if (abs(para_version-2.14) < 1.e-4) then
          write(*,*) "Trying to read 2.14 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para214(para)
          return
       else if (abs(para_version-2.13) < 1.e-4) then
          write(*,*) "Trying to read 2.13 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para213(para)
          return
       else if (abs(para_version-2.12) < 1.e-4) then
          write(*,*) "Trying to read 2.12 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para212(para)
          return
       else if (abs(para_version-2.11) < 1.e-4) then ! for the DENT grid
          write(*,*) "Trying to read 2.11 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para211(para)
          return
       else
          close(unit=1)
          call error("Unsupported version of the parameter file")
       endif
    endif

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max

    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, map_size
    zoom = 1.0

    ! default values for MC inclinations and azimuth bins
    N_thet = 10 ; N_phi = 1

    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    read(1,*) RT_az_min, RT_az_max, RT_n_az
    if ((RT_n_incl < 1).or.(RT_n_az < 1)) call error("The number of inclination and azimuth must be >= 1")
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) settling_type, exp_strat, a_strat
    if (settling_type == 0) then
       lvariable_dust = .false.
    else
       lvariable_dust = .true.
       if (exp_strat < 0.) then
          exp_strat = -exp_strat
          write(*,*) "Setting exp_strat > 0"
       endif
    endif
    read(1,*) lmigration
    read(1,*,IOSTAT=status) ldust_sublimation , correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lhydrostatic
    read(1,*) lchauff_int, alpha
    T_min= 1.0 ; T_max=3000. ; n_T=100
    if (lchange_Tmax_PAH) T_max = Tmax_PAH

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       if ((disk_zone(j)%geometry == 4).and.(settling_type > 0)) then
          write(*,*) "WARNING : debris disk, setting settling to 0"
          settling_type=0
          lvariable_dust = .false.
       endif

       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref, disk_zone(j)%vert_exponent
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%edge, disk_zone(j)%Rout, disk_zone(j)%Rc
       disk_zone(j)%Rmax = disk_zone(j)%Rout
       if ((disk_zone(j)%geometry == 2).and.(disk_zone(j)%Rout < tiny_real)) then
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc ! tappered-edge
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf, disk_zone(j)%moins_gamma_exp
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%type, dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass, dust_pop_tmp(n_pop)%dhs_maxf
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          if ((dust_pop_tmp(n_pop)%lcoating) .and. ((dust_pop_tmp(n_pop)%type=="DHS").or. &
               (dust_pop_tmp(n_pop)%type=="dhs")) ) then
             call error("cannot use DHS and coating for the same dust grains")
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) &
                  / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j

          ! Checking which type of opacity file it is
          dust_pop_tmp(n_pop)%is_PAH = .false.
          dust_pop_tmp(n_pop)%is_opacity_file = .false.
          dust_pop_tmp(n_pop)%is_Misselt_opacity_file = .false.
          dust_pop_tmp(n_pop)%is_DustEM_opacity_file = .false.
          if (dust_pop_tmp(n_pop)%indices(1)(1:3) == "PAH") then
             dust_pop_tmp(n_pop)%is_PAH = .true. ; dust_pop_tmp(n_pop)%is_opacity_file = .true.
          else  if (dust_pop_tmp(n_pop)%indices(1)(1:7) == "Misselt") then
             dust_pop_tmp(n_pop)%is_opacity_file = .true.
             dust_pop_tmp(n_pop)%is_Misselt_opacity_file = .true.
             lread_Misselt = .true.
             if (dust_pop_tmp(n_pop)%indices(1)(9:11) == "PAH") then
                dust_pop_tmp(n_pop)%is_PAH = .true.
             endif
          else  if (dust_pop_tmp(n_pop)%indices(1)(1:6) == "DustEM") then
             dust_pop_tmp(n_pop)%is_opacity_file = .true.
             dust_pop_tmp(n_pop)%is_DustEM_opacity_file = .true.
             lread_DustEM = .true.
             if (dust_pop_tmp(n_pop)%indices(1)(8:10) == "PAH") then
                dust_pop_tmp(n_pop)%is_PAH = .true.
             endif
             ! Removing DustEM prefix
             dust_pop_tmp(n_pop)%indices(1) = trim(dust_pop_tmp(n_pop)%indices(1)(8:512))
          endif
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    !dust_pop%is_PAH = .false.
    !dust_pop%is_opacity_file = .false.
    !dust_pop%is_Misselt_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 3000.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             if (dust_pop(ind_pop)%indices(1)(1:7) == "Misselt") then
                dust_pop(ind_pop)%is_opacity_file = .true.
                dust_pop(ind_pop)%is_Misselt_opacity_file = .true.
                lread_Misselt = .true.
                if (dust_pop(ind_pop)%indices(1)(9:11) == "PAH") then
                   dust_pop(ind_pop)%is_PAH = .true.
                   T_max = 3000.
                   if (lchange_Tmax_PAH) T_max = Tmax_PAH
                endif
             endif

             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU
       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para

  !**********************************************************************

  subroutine read_para220(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max

    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, map_size
    zoom = 1.0
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    read(1,*) RT_az_min, RT_az_max, RT_n_az
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) settling_type, exp_strat, a_strat
    if (settling_type == 0) then
       lvariable_dust = .false.
    else
       lvariable_dust = .true.
       if (exp_strat < 0.) then
          exp_strat = -exp_strat
          write(*,*) "Setting exp_strat > 0"
       endif
    endif
    read(1,*) lmigration
    read(1,*,IOSTAT=status) ldust_sublimation , correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lhydrostatic
    read(1,*) lchauff_int, alpha
    T_min= 1.0 ; T_max=3000. ; n_T=100
    if (lchange_Tmax_PAH) T_max = Tmax_PAH

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       if ((disk_zone(j)%geometry == 4).and.(settling_type > 0)) then
          write(*,*) "WARNING : debris disk, setting settling to 0"
          settling_type=0
          lvariable_dust = .false.
       endif

       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref, disk_zone(j)%vert_exponent
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%edge, disk_zone(j)%Rout, disk_zone(j)%Rc
       disk_zone(j)%Rmax = disk_zone(j)%Rout
       if ((disk_zone(j)%geometry == 2).and.(disk_zone(j)%Rout < tiny_real)) then
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc ! tappered-edge
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf, disk_zone(j)%moins_gamma_exp
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%type, dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass, dust_pop_tmp(n_pop)%dhs_maxf
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          if ((dust_pop_tmp(n_pop)%lcoating) .and. ((dust_pop_tmp(n_pop)%type=="DHS").or. &
               (dust_pop_tmp(n_pop)%type=="dhs")) ) then
             call error("cannot use DHS and coating for the same dust grains")
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) &
                  / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j

          ! Checking which type of opacity file it is
          dust_pop_tmp(n_pop)%is_PAH = .false.
          dust_pop_tmp(n_pop)%is_opacity_file = .false.
          dust_pop_tmp(n_pop)%is_Misselt_opacity_file = .false.
          if (dust_pop_tmp(n_pop)%indices(1)(1:3) == "PAH") then
             dust_pop_tmp(n_pop)%is_PAH = .true. ; dust_pop_tmp(n_pop)%is_opacity_file = .true.
          else  if (dust_pop_tmp(n_pop)%indices(1)(1:7) == "Misselt") then
             dust_pop_tmp(n_pop)%is_opacity_file = .true.
             dust_pop_tmp(n_pop)%is_Misselt_opacity_file = .true.
             lread_Misselt = .true.
             if (dust_pop_tmp(n_pop)%indices(1)(9:11) == "PAH") then
                dust_pop_tmp(n_pop)%is_PAH = .true.
             endif
          endif

       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    !dust_pop%is_PAH = .false.
    !dust_pop%is_opacity_file = .false.
    !dust_pop%is_Misselt_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 3000.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             if (dust_pop(ind_pop)%indices(1)(1:7) == "Misselt") then
                dust_pop(ind_pop)%is_opacity_file = .true.
                dust_pop(ind_pop)%is_Misselt_opacity_file = .true.
                lread_Misselt = .true.
                if (dust_pop(ind_pop)%indices(1)(9:11) == "PAH") then
                   dust_pop(ind_pop)%is_PAH = .true.
                   T_max = 3000.
                   if (lchange_Tmax_PAH) T_max = Tmax_PAH
                endif
             endif

             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU
       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para220

  !**********************************************************************

  subroutine read_para219(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, map_size
    zoom = 1.0
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) settling_type, exp_strat, a_strat
    if (settling_type == 0) then
       lvariable_dust = .false.
    else
       lvariable_dust = .true.
       if (exp_strat < 0.) then
          exp_strat = -exp_strat
          write(*,*) "Setting exp_strat > 0"
       endif
    endif
    read(1,*) lmigration
    read(1,*,IOSTAT=status) ldust_sublimation , correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lhydrostatic
    read(1,*) lchauff_int, alpha
    T_min= 1.0 ; T_max=3000. ; n_T=100
    if (lchange_Tmax_PAH) T_max = Tmax_PAH

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       if ((disk_zone(j)%geometry == 4).and.(settling_type > 0)) then
          write(*,*) "WARNING : debris disk, setting settling to 0"
          settling_type=0
          lvariable_dust = .false.
       endif

       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref, disk_zone(j)%vert_exponent
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%edge, disk_zone(j)%Rout, disk_zone(j)%Rc
       disk_zone(j)%Rmax = disk_zone(j)%Rout
       if ((disk_zone(j)%geometry == 2).and.(disk_zone(j)%Rout < tiny_real)) then
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc ! tappered-edge
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf, disk_zone(j)%moins_gamma_exp
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%type, dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass, dust_pop_tmp(n_pop)%dhs_maxf
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          if ((dust_pop_tmp(n_pop)%lcoating) .and. ((dust_pop_tmp(n_pop)%type=="DHS").or. &
               (dust_pop_tmp(n_pop)%type=="dhs")) ) then
             call error("cannot use DHS and coating for the same dust grains")
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) &
                  / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.
    dust_pop%is_Misselt_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 3000.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             if (dust_pop(ind_pop)%indices(1)(1:7) == "Misselt") then
                dust_pop(ind_pop)%is_opacity_file = .true.
                dust_pop(ind_pop)%is_Misselt_opacity_file = .true.
                lread_Misselt = .true.
                if (dust_pop(ind_pop)%indices(1)(9:11) == "PAH") then
                   dust_pop(ind_pop)%is_PAH = .true.
                   T_max = 3000.
                   if (lchange_Tmax_PAH) T_max = Tmax_PAH
                endif
             endif

             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) call error("you are using a 2-zone disk parameter file", &
         msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para219

  !**********************************************************************

  subroutine read_para218(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, map_size
    zoom = 1.0
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) settling_type, exp_strat, a_strat
    if (settling_type == 0) then
       lvariable_dust = .false.
    else
       lvariable_dust = .true.
       if (exp_strat < 0.) then
          exp_strat = -exp_strat
          write(*,*) "Setting exp_strat > 0"
       endif
    endif
    read(1,*) lmigration
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lhydrostatic
    read(1,*) lchauff_int, alpha
    T_min=1. ; T_max=1500. ; n_T=100

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref
       disk_zone(j)%vert_exponent = 2.
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%edge, disk_zone(j)%Rout, disk_zone(j)%Rc
       disk_zone(j)%Rmax = disk_zone(j)%Rout
       if ((disk_zone(j)%geometry == 2).and.(disk_zone(j)%Rout < tiny_real)) then
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc ! tappered-edge
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf, disk_zone(j)%moins_gamma_exp
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%type, dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass, dust_pop_tmp(n_pop)%dhs_maxf
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          if ((dust_pop_tmp(n_pop)%lcoating) .and. ((dust_pop_tmp(n_pop)%type=="DHS").or. &
               (dust_pop_tmp(n_pop)%type=="dhs")) ) then
             call error("cannot use DHS and coating for the same dust grains")
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) &
                  / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para218

  !**********************************************************************

  subroutine read_para217(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, map_size, zoom
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) settling_type, exp_strat, a_strat
    if (settling_type == 0) then
       lvariable_dust = .false.
    else
       lvariable_dust = .true.
       if (exp_strat < 0.) then
          exp_strat = -exp_strat
          write(*,*) "Setting exp_strat > 0"
       endif
    endif
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    T_min=1. ; T_max=1500. ; n_T=100

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%Rout, disk_zone(j)%edge
       if (disk_zone(j)%geometry == 2) then ! tappered-edge
          disk_zone(j)%Rc =  disk_zone(j)%Rout
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc
       else
          disk_zone(j)%Rmax = disk_zone(j)%Rout
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf, disk_zone(j)%moins_gamma_exp
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%type, dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass, dust_pop_tmp(n_pop)%dhs_maxf
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          if ((dust_pop_tmp(n_pop)%lcoating) .and. ((dust_pop_tmp(n_pop)%type=="DHS").or. &
               (dust_pop_tmp(n_pop)%type=="dhs")) ) then
             call error("cannot use DHS and coating for the same dust grains")
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) &
                  / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para217

  !**********************************************************************

  subroutine read_para216(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_dp
    dust_pop_tmp(:)%dhs_maxf = 0.9

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, map_size, zoom
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) settling_type, exp_strat, a_strat
    if (settling_type == 0) then
       lvariable_dust = .false.
    else
       lvariable_dust = .true.
       if (exp_strat < 0.) then
          exp_strat = -exp_strat
          write(*,*) "Setting exp_strat > 0"
       endif
    endif
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    T_min=1. ; T_max=1500. ; n_T=100

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%Rout, disk_zone(j)%edge
       if (disk_zone(j)%geometry == 2) then ! tappered-edge
          disk_zone(j)%Rc =  disk_zone(j)%Rout
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc
       else
          disk_zone(j)%Rmax = disk_zone(j)%Rout
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf
       disk_zone(j)%moins_gamma_exp = disk_zone(j)%surf
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%type, dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          if ((dust_pop_tmp(n_pop)%lcoating) .and. ((dust_pop_tmp(n_pop)%type=="DHS").or. &
               (dust_pop_tmp(n_pop)%type=="dhs")) ) then
             call error("cannot use DHS and coating for the same dust grains")
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) &
                  / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.


    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para216

  !**********************************************************************

  subroutine read_para215(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_dp
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) npix_x, npix_y, map_size, zoom
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) lvariable_dust, exp_strat, a_strat
    if (exp_strat < 0.) then
       exp_strat = -exp_strat
       write(*,*) "Setting exp_strat > 0"
    endif
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    T_min=1. ; T_max=1500. ; n_T=100

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref
       read(1,*) disk_zone(j)%Rin, disk_zone(j)%Rout, disk_zone(j)%edge
       if (disk_zone(j)%geometry == 2) then ! tappered-edge
          disk_zone(j)%Rc =  disk_zone(j)%Rout
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc
       else
          disk_zone(j)%Rmax = disk_zone(j)%Rout
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf
       disk_zone(j)%moins_gamma_exp = disk_zone(j)%surf
    enddo ! n_zones

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return


end subroutine read_para215

!**********************************************************************

 subroutine read_para214(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: size_neb_tmp, somme, V_somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_dp
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)
    tau_seuil  = 1.0e31
    wl_seuil = 0.81


    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*) npix_x, npix_y, zoom
    read(1,*,iostat=ios) N_thet, N_phi
    capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif

    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) gas_dust
    read(1,*) lvariable_dust, exp_strat, a_strat
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    T_min=1. ; T_max=1500. ; n_T=100

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry <=2) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
       if (disk_zone(j)%geometry == 2) then ! tappered-edge
          disk_zone(j)%Rc =  disk_zone(j)%Rout
          disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc
       else
          disk_zone(j)%Rmax = disk_zone(j)%Rout
       endif
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf
       disk_zone(j)%moins_gamma_exp = disk_zone(j)%surf

       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) call error("different values for size_neb")
       endif
    enddo ! n_zones

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    diskmass = sum(disk_zone(:)%diskmass)

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)


    return

  end subroutine read_para214


  !**********************************************************************

 subroutine read_para213(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: size_neb_tmp, somme, V_somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_dp
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)
    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*,iostat=ios) N_thet, N_phi, npix_x, npix_y, zoom
    read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif
    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) gas_dust
    read(1,*) lvariable_dust, exp_strat, a_strat
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry ==1) is_there_disk = .true.
       if (disk_zone(j)%geometry == 2) disk_zone(j)%geometry = 3 ! update avec tappered-edge
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf

       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) call error("different values for size_neb")
       endif
    enddo ! n_zones

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) call error("r_min < 0.0")

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          !read(1,*) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*,iostat=ios) dust_pop_tmp(n_pop)%n_components, dust_pop_tmp(n_pop)%mixing_rule, &
               dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          if ( (dust_pop_tmp(n_pop)%n_components > 1).and.(dust_pop_tmp(n_pop)%mixing_rule == 2) ) then
             dust_pop_tmp(n_pop)%lcoating = .true.
          else
             dust_pop_tmp(n_pop)%lcoating = .false.
          endif
          V_somme = 0.0
          do k=1, dust_pop_tmp(n_pop)%n_components
             read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
             V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
          enddo
          if (V_somme < tiny_real) then
             write(*,*) "ERROR: population #", n_pop,  ": sum of volume fraction is 0"
             write(*,*) "Exiting"
             call exit(1)
          endif
          ! renormalisation des fraction en volume
          do k=1, dust_pop_tmp(n_pop)%n_components
             dust_pop_tmp(n_pop)%component_volume_fraction(k) = dust_pop_tmp(n_pop)%component_volume_fraction(k) / V_somme
          enddo
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo
    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para213

  !**********************************************************************

  subroutine read_para212(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, alloc_status, ios, ind_pop, imol, status
    real(kind=dp) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes
    character(len=100) :: line_buffer

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_dp
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)
    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*,iostat=ios) N_thet, N_phi, npix_x, npix_y, zoom
    read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif
    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    RT_az_min = 0.0 ; RT_az_max = 0.0 ; RT_n_az = 1 ;
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale

    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) gas_dust
    read(1,*) lvariable_dust, exp_strat, a_strat
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lvariable_dust is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry ==1) is_there_disk = .true.
       if (disk_zone(j)%geometry == 2) disk_zone(j)%geometry = 3 ! update avec tappered-edge
       if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf

       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) call error("different values for size_neb")
       endif
    enddo ! n_zones

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) call error("r_min < 0.0")

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo

    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ; dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * km_to_m ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * km_to_m ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error('you are using a 2-zone disk parameter file', &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*) line_buffer
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para212

  !**********************************************************************

  subroutine read_para211(para)

    implicit none

    character(len=*), intent(in) :: para

    integer :: i, j, alloc_status, ios, tmpint, ind_pop, status, imol

    real(kind=dp) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    character(len=100) :: line_buffer

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_dp
    dust_pop_tmp(:)%type = "Mie"
    dust_pop_tmp(:)%mixing_rule=1

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*) line_buffer
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    nbre_photons_tot=real(nbre_photons_loop)*real(nbre_photons_eq_th)
    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*) line_buffer
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*) line_buffer
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif

    ! ----
    ! Maps
    ! ----
    read(1,*) line_buffer
    read(1,*,iostat=ios) N_thet, N_phi, npix_x, npix_y, zoom
    read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet

    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    if (lonly_capt_interet) then
       N_incl = capt_sup - capt_inf + 1
       capt_debut = capt_inf
       capt_fin = capt_sup
    else
       N_incl = N_thet
       capt_debut = 1
       capt_fin = N_thet
    endif
    read(1,*) RT_imin, RT_imax, RT_n_incl, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque

    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*) line_buffer
    read(1,*) scattering_method0 ; scattering_method = scattering_method0
    read(1,*) aniso_method

    ! ----------
    ! Symmetries
    ! ----------
    read(1,*) line_buffer
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale


    ! ----------------------
    ! Dust global properties
    ! ----------------------
    read(1,*) line_buffer
    read(1,*) gas_dust
    read(1,*) lvariable_dust, exp_strat, a_strat
    read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
    if (status/=0) correct_Rsub = 1.0
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T

    ! ---------------
    ! Number of zones
    ! ---------------
    read(1,*) line_buffer
    read(1,*) n_zones
    if (n_zones > 1) then
       lvariable_dust=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lstrat is set to true and exp_strat to 0."
    endif
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error disk parameters')
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    read(1,*) line_buffer
    do j=1,n_zones
       read(1,*) disk_zone(j)%geometry
       if (disk_zone(j)%geometry ==1) is_there_disk = .true.
       if ((disk_zone(j)%geometry == 2).and.(grid_type == 1)) then
          write(*,*) "WARNING : you are using an envelope density structure"
          write(*,*) "          with a cylindrical grid !!!!"
       endif
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
           disk_zone(j)%Rmax = disk_zone(j)%Rout
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf

       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) call error("different values for size_neb")
       endif
    enddo ! n_zones

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    rmax = maxval(disk_zone(:)%rout)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) call error("r_min < 0.0")

    ! ------
    ! Cavity
    ! ------
    read(1,*) line_buffer
    read(1,*) lcavity
    read(1,*) cavity%sclht, cavity%rref
    read(1,*) cavity%exp_beta

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*) line_buffer
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
          read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
          if (dust_pop_tmp(n_pop)%methode_chauffage == 1) lRE_LTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 2) lRE_nLTE=.true.
          if (dust_pop_tmp(n_pop)%methode_chauffage == 3) lnRE=.true.
          somme = somme + dust_pop_tmp(n_pop)%frac_mass
          dust_pop_tmp(n_pop)%zone = j

          dust_pop_tmp(n_pop)%n_components = 1
          dust_pop_tmp(n_pop)%component_volume_fraction(1) =  1.0
          dust_pop_tmp(n_pop)%is_PAH = .false.
          dust_pop_tmp(n_pop)%is_opacity_file = .false.
          dust_pop_tmp(n_pop)%is_Misselt_opacity_file = .false.
          dust_pop_tmp(n_pop)%is_DustEM_opacity_file = .false.
       enddo

       ! renormalisation des fraction en masse
       do i=1,n_especes(j)
          dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass = dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass/somme
          dust_pop_tmp(n_pop-n_especes(j)+i)%masse =  dust_pop_tmp(n_pop-n_especes(j)+i)%frac_mass * disk_zone(j)%diskmass
       enddo

    enddo !n_zones

    if (lRE_LTE.and.lRE_nLTE) call error("cannot mix grains in LTE and nLTE")

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error n_pop tmp')
    dust_pop%is_PAH = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 1) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains

          endif
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 2) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
             grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end
    if (lnRE) then
       do i=1, n_pop
          if (dust_pop_tmp(i)%methode_chauffage == 3) then
             ind_pop=ind_pop+1
             dust_pop(ind_pop) = dust_pop_tmp(i)
             if (dust_pop(ind_pop)%indices(1)(1:3) == "PAH") then
                dust_pop(ind_pop)%is_PAH = .true. ;dust_pop(ind_pop)%is_opacity_file = .true.
                T_max = 2500.
                if (lchange_Tmax_PAH) T_max = Tmax_PAH
             endif
             dust_pop(ind_pop)%ind_debut = grain_nRE_end + 1
             dust_pop(ind_pop)%ind_fin = grain_nRE_end + dust_pop(ind_pop)%n_grains
             grain_nRE_end = grain_nRE_end +  dust_pop(ind_pop)%n_grains
          endif
       enddo
    endif

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    ! ---------------------
    ! Molecular RT settings
    ! ---------------------
    read(1,*) line_buffer
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) vitesse_turb
    vitesse_turb = vitesse_turb * 1.e3 ! Conversion en m.s-1
    read(1,*) n_molecules
    allocate(mol(n_molecules))
    do imol=1,n_molecules
       read(1,*) mol(imol)%filename, mol(imol)%iLevel_max
       read(1,*) mol(imol)%vmax_center_rt, mol(imol)%n_speed_rt
       mol(imol)%vmax_center_rt = mol(imol)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
       read(1,*) mol(imol)%lcst_abundance, mol(imol)%abundance, mol(imol)%abundance_file
       read(1,*) mol(imol)%lline, mol(imol)%nTrans_raytracing
       read(1,*) mol(imol)%indice_Trans_rayTracing(1:mol(imol)%nTrans_raytracing)
       mol(imol)%n_speed_center_rt = mol(imol)%n_speed_rt
       mol(imol)%n_extraV_rt = 0 ; mol(imol)%extra_deltaV_rt = 0.0
    enddo

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*) line_buffer
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       call error("you are using a 2-zone disk parameter file", &
            msg2='You must use the [-2zone] option to calculate a 2-zone disk')
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) call error('Allocation error etoile')

    if (n_etoiles > 1) then
       write(*,*) "Multiple illuminating stars! Cancelling all image symmetries"
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    endif

    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body, etoile(i)%fUV
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU
    enddo
    etoile(:)%vx = 0.0 ; etoile(:)%vy = 0.0 ; etoile(:)%vz = 0.0
    etoile(:)%Mdot = 0.0

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para211

  !**********************************************************************

end module read_params
