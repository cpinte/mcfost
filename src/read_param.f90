module read_params

  use parametres
  use grains
  use disk
  use prop_star
  use em_th
  use molecular_emission
  use ray_tracing
  use input
  use sha

  implicit none

contains

  subroutine read_para()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old', iostat=ios)
    if (ios/=0) then
       write(*,*) "ERROR : cannot open "//trim(para)
       write(*,*) "Exiting"
       stop
    endif

    read(1,*) para_version

    ! Petit test pour le ray-tracing
    if (lscatt_ray_tracing .and. (abs(para_version) < 2.0901)) then
       write(*,*) "Parameter version >= 2.10 required for ray-tracing."
       write(*,*) "Exiting."
       stop
    endif

    if (para_version < 0.) then
       if (ldebris) then
          write(*,*) "You cannot use the -debris option to feed a FITS file density table"
          write(*,*) "Exiting!"
          stop
       else
          write(*,*) "You are running MCFOST with a user-specified density structure"
          write(*,*) ""
          lfits=.true.
          para_version=-para_version
       endif
    endif

    correct_Rsub = 1.0_db
    lmigration = .false.
    lhydrostatic = .false.
    lread_misselt=.false.

    if (abs(para_version - 2.19) > 1.e-4) then
       write(*,*) "Wrong version of the parameter file."
      if (abs(para_version-2.18) < 1.e-4) then
          write(*,*) "Trying to read 2.17 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para218()
          return
       else if (abs(para_version-2.17) < 1.e-4) then
          write(*,*) "Trying to read 2.17 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para217()
          return
       else if (abs(para_version-2.16) < 1.e-4) then
          write(*,*) "Trying to read 2.16 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para216()
          return
       else if (abs(para_version-2.15) < 1.e-4) then
          write(*,*) "Trying to read 2.15 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para215()
          return
       else if (abs(para_version-2.14) < 1.e-4) then
          write(*,*) "Trying to read 2.14 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para214()
          return
       else if (abs(para_version-2.13) < 1.e-4) then
          write(*,*) "Trying to read 2.13 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para213()
          return
       else if (abs(para_version-2.12) < 1.e-4) then
          write(*,*) "Trying to read 2.12 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para212()
          return
       else if (abs(para_version-2.11) < 1.e-4) then
          write(*,*) "Trying to read 2.11 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para211()
          return
       else if (abs(para_version-2.10) < 1.e-4) then
          write(*,*) "Trying to read 2.10 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para210()
          return
       else if (abs(para_version-2.09) < 1.e-4) then
          write(*,*) "Trying to read 2.09 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para209()
          return
       else if (abs(para_version-2.08) < 1.e-4) then
          write(*,*) "Trying to read 2.08 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para208()
          return
       else if (abs(para_version-2.07) < 1.e-4) then
          write(*,*) "Trying to read 2.07 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para207()
          return
       else if (abs(para_version-2.06) < 1.e-4) then
          write(*,*) "Trying to read 2.06 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para206()
          return
       else if (abs(para_version-2.05) < 1.e-4) then
          write(*,*) "Trying to read 2.05 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para205()
          return
       else if (abs(para_version-2.04) < 1.e-4) then
          write(*,*) "Trying to read 2.04 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para204()
          return
       else if (abs(para_version-2.03) < 1.e-4) then
          write(*,*) "Trying to read 2.03 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para203()
          return
       else if (abs(para_version-2.02) < 1.e-4) then
          write(*,*) "Trying to read 2.02 parameter file."
          write(*,*) "Pbs can appear. Parameter file should be updated !!!"
          close(unit=1)
          call read_para202()
          return
       else
          close(unit=1)
          write(*,*) "Unsupported version of the parameter file :", para_version
          write(*,*) "Exiting."
          stop
       endif
    endif

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) igridx, igridy, map_size
    zoom = 1.0
    read(1,*,iostat=ios) N_thet, N_phi, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo

    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    endif  ! lfits
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

    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size ! compatibilite avec size_neb
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) settling_type, exp_strat, a_strat
       if (settling_type == 0) then
          lstrat = .false.
       else
          lstrat = .true.
       endif
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*) lmigration
       read(1,*,IOSTAT=status) ldust_sublimation , correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lhydrostatic
       read(1,*) lchauff_int, alpha
       T_min= 1.0 ; T_max=3000. ; n_T=300
       if (lchange_Tmax_PAH) T_max = Tmax_PAH
    endif  ! lfits

    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
          disk_zone(j)%gas_to_dust = 100.
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
             lstrat = .false.
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
    endif ! lfits

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    Rmin = minval(disk_zone(:)%Rmin)
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (Rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
                write(*,*) "ERROR: cannot use DHS and coating for the same dust garins"
                write(*,*) "Exiting"
                stop
             endif
             V_somme = 0.0
             do k=1, dust_pop_tmp(n_pop)%n_components
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
                V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
             enddo
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
    endif ! lfits

!    if (lRE_LTE.and.lRE_nLTE) then
!       write(*,*) "Error : cannot mix grains in LTE and nLTE"
!       write(*,*) " Is it usefull anyway ???"
!       write(*,*) "Exiting"
!       stop
!    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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

    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU


       if ( (abs(etoile(i)%x) < tiny_real).and.(abs(etoile(i)%x) < tiny_real).and.(abs(etoile(i)%x) < tiny_real) ) then
          if (etoile(i)%r > Rmin) then
             write(*,*) "ERROR : inner disk radius is smaller than stellar radius"
             write(*,*) "Exiting"
             stop
          endif
       endif

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para

  !**********************************************************************

  subroutine read_para218()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) igridx, igridy, map_size
    zoom = 1.0
    read(1,*,iostat=ios) N_thet, N_phi, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo

    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    endif  ! lfits
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

    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size ! compatibilite avec size_neb
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) settling_type, exp_strat, a_strat
       if (settling_type == 0) then
          lstrat = .false.
       else
          lstrat = .true.
       endif
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*) lmigration
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lhydrostatic
       read(1,*) lchauff_int, alpha
       T_min=1. ; T_max=1500. ; n_T=100
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
          disk_zone(j)%gas_to_dust = 100.
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
       do j=1,n_zones
          read(1,*) disk_zone(j)%geometry
          if (disk_zone(j)%geometry <=2) is_there_disk = .true.
          if ((disk_zone(j)%geometry == 3).and.(grid_type == 1)) then
             write(*,*) "WARNING : you are using an envelope density structure"
             write(*,*) "          with a cylindrical grid !!!!"
          endif
          read(1,*) disk_zone(j)%diskmass, disk_zone(j)%gas_to_dust
          read(1,*) disk_zone(j)%sclht, disk_zone(j)%Rref
          read(1,*) disk_zone(j)%Rin, disk_zone(j)%edge, disk_zone(j)%Rout, disk_zone(j)%Rc
          disk_zone(j)%Rmax = disk_zone(j)%Rout
          if ((disk_zone(j)%geometry == 2).and.(disk_zone(j)%Rout < tiny_real)) then
             disk_zone(j)%Rmax = 8 * disk_zone(j)%Rc ! tappered-edge
          endif
          read(1,*) disk_zone(j)%exp_beta
          read(1,*) disk_zone(j)%surf, disk_zone(j)%moins_gamma_exp
       enddo ! n_zones
    endif ! lfits

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    Rmin = minval(disk_zone(:)%Rmin)
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (Rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
                write(*,*) "ERROR: cannot use DHS and coating for the same dust garins"
                write(*,*) "Exiting"
                stop
             endif
             V_somme = 0.0
             do k=1, dust_pop_tmp(n_pop)%n_components
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
                V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
             enddo
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para218

  !**********************************************************************

  subroutine read_para217()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) igridx, igridy, map_size, zoom
    read(1,*,iostat=ios) N_thet, N_phi, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo

    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    endif  ! lfits
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

    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size ! compatibilite avec size_neb
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) settling_type, exp_strat, a_strat
       if (settling_type == 0) then
          lstrat = .false.
       else
          lstrat = .true.
       endif
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       T_min=1. ; T_max=1500. ; n_T=100
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
          disk_zone(j)%gas_to_dust = 100.
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
    endif ! lfits

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    Rmin = minval(disk_zone(:)%Rmin)
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (Rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
                write(*,*) "ERROR: cannot use DHS and coating for the same dust garins"
                write(*,*) "Exiting"
                stop
             endif
             V_somme = 0.0
             do k=1, dust_pop_tmp(n_pop)%n_components
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
                V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
             enddo
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para217

  !**********************************************************************

  subroutine read_para216()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_db
    dust_pop_tmp(:)%dhs_maxf = 0.9

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) igridx, igridy, map_size, zoom
    read(1,*,iostat=ios) N_thet, N_phi, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo

    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    endif  ! lfits
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

    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size ! compatibilite avec size_neb
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) settling_type, exp_strat, a_strat
       if (settling_type == 0) then
          lstrat = .false.
       else
          lstrat = .true.
       endif
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       T_min=1. ; T_max=1500. ; n_T=100
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
          disk_zone(j)%gas_to_dust = 100.
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
    endif ! lfits

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    Rmin = minval(disk_zone(:)%Rmin)
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (Rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
                write(*,*) "ERROR: cannot use DHS and coating for the same dust garins"
                write(*,*) "Exiting"
                stop
             endif
             V_somme = 0.0
             do k=1, dust_pop_tmp(n_pop)%n_components
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(k), dust_pop_tmp(n_pop)%component_volume_fraction(k)
                V_somme = V_somme + dust_pop_tmp(n_pop)%component_volume_fraction(k)
             enddo
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para216

  !**********************************************************************

  subroutine read_para215()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: somme, V_somme

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    real :: fnbre_photons_eq_th, fnbre_photons_lambda, fnbre_photons_image

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_db
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) fnbre_photons_eq_th ;
    read(1,*) fnbre_photons_lambda ;
    read(1,*) fnbre_photons_image
    nbre_photons_loop = 128 ;
    nbre_photons_eq_th = max(fnbre_photons_eq_th / nbre_photons_loop,1.)
    nbre_photons_lambda = max(fnbre_photons_lambda / nbre_photons_loop,1.)
    nbre_photons_image = max(fnbre_photons_image / nbre_photons_loop,1.)

    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) igridx, igridy, map_size, zoom
    read(1,*,iostat=ios) N_thet, N_phi, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo

    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    endif  ! lfits
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

    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size ! compatibilite avec size_neb
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) lstrat, exp_strat, a_strat
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       T_min=1. ; T_max=1500. ; n_T=100
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
          disk_zone(j)%gas_to_dust = 100.
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
    endif ! lfits

    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    Rmin = minval(disk_zone(:)%Rmin)
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (Rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return


end subroutine read_para215


!**********************************************************************

 subroutine read_para214()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: size_neb_tmp, somme, V_somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_db
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    tau_seuil  = 1.0e31
    wl_seuil = 0.81


    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*) igridx, igridy, zoom
    read(1,*,iostat=ios) N_thet, N_phi, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo

    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       capt_interet= 1     ; delta_capt=1 ; angle_interet=75. ; lonly_capt_interet=.false.
    endif  ! lfits
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

    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       gas_dust=100.
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) gas_dust
       read(1,*) lstrat, exp_strat, a_strat
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       T_min=1. ; T_max=1500. ; n_T=100
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
             if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
                write(*,*) "Error : different values for size_neb"
                write(*,*) "Exiting"
                stop
             endif
          endif
       enddo ! n_zones
    endif ! lfits

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    Rmin = minval(disk_zone(:)%Rmin)
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)


    return

  end subroutine read_para214


  !**********************************************************************

 subroutine read_para213()

    implicit none

    integer :: i, j, k, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: size_neb_tmp, somme, V_somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_db
    dust_pop_tmp(:)%type = "Mie"

    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet
    endif  ! lfits
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
    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2* map_size
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       gas_dust=100.
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) gas_dust
       read(1,*) lstrat, exp_strat, a_strat
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       read(1,*) T_min, T_max, n_T
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
             if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
                write(*,*) "Error : different values for size_neb"
                write(*,*) "Exiting"
                stop
             endif
          endif
       enddo ! n_zones
    endif ! lfits

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else ! lfits
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para213

  !**********************************************************************

  subroutine read_para212()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop, imol, status
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    correct_Rsub = 1.0_db
    dust_pop_tmp(:)%type = "Mie"


    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif


    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet
    endif  ! lfits
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
    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       gas_dust=100.
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) gas_dust
       read(1,*) lstrat, exp_strat, a_strat
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       read(1,*) T_min, T_max, n_T
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
             if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
                write(*,*) "Error : different values for size_neb"
                write(*,*) "Exiting"
                stop
             endif
          endif
       enddo ! n_zones
    endif ! lfits

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices, dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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


    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU

       read(1,*) etoile(i)%fUV, etoile(i)%slope_UV
       etoile(i)%slope_UV = etoile(i)%slope_UV - 2.0  ! Fnu -> F_lambda
    enddo

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para212

  !**********************************************************************

  subroutine read_para211()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop, status, imol
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    tau_seuil  = 1.0e31
    wl_seuil = 0.81

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif

    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet
    endif  ! lfits
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
    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       gas_dust=100.
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) gas_dust
       read(1,*) lstrat, exp_strat, a_strat
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*,IOSTAT=status) ldust_sublimation, correct_Rsub
       if (status/=0) correct_Rsub = 1.0
       read(1,*) lchauff_int, alpha
       read(1,*) T_min, T_max, n_T
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
             if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
                write(*,*) "Error : different values for size_neb"
                write(*,*) "Exiting"
                stop
             endif
          endif
       enddo ! n_zones
    endif ! lfits

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif


    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       read(1,*)
       read(1,*)
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

    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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

    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para211

  !**********************************************************************


  subroutine read_para210()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    ! -------------------------
    ! Number of photon packages
    ! -------------------------
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image
    read(1,*) !lcheckpoint,  checkpoint_time

    ! ----------
    ! Wavelength
    ! ----------
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

   if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif

    if (lfits) then
       tau_seuil=1.0e16; wl_seuil=0.81
    else
       read(1,*) tau_seuil, wl_seuil
    endif  ! lfits
    ! -------------------------------
    ! Grid geometry / input FITS file
    ! -------------------------------
    read(1,*)
    read(1,*)
    if (.not.lfits) then
       read(1,*) grid_type
       read(1,*) n_rad, nz, n_az, n_rad_in
    else
       read(1,*) struct_fits_file
       call read_struct_fits_file()
    endif ! lfits

    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    ! ----
    ! Maps
    ! ----
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    maxigrid = max(igridx, igridy)

    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    if (lfits) then
       capt_interet=N_thet ; delta_capt=1 ; angle_interet=90. ; lonly_capt_interet=.false.
    else
       read(1,*) capt_interet, delta_capt, angle_interet, lonly_capt_interet
    endif  ! lfits
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
    read(1,*) RT_imin, RT_imax, RT_n_ibin, lRT_i_centered
    read(1,*) distance
    read(1,*) ang_disque
    if (lfits) then
       read(1,*) map_size
       map_size = 2*map_size
    endif  ! lfits
    ! -----------------
    ! Scattering method
    ! -----------------
    read(1,*)
    read(1,*)
    if (lfits) then
       scattering_method=0  ! Use "auto" mode to store dust properties
    else
       read(1,*) scattering_method
    endif ! lfits
    read(1,*) aniso_method
    ! ----------
    ! Symmetries
    ! ----------
    if (lfits) then
       l_sym_ima=.false.
       l_sym_centrale=.false.
       l_sym_axiale=.false.
    else
       read(1,*)
       read(1,*)
       read(1,*) l_sym_ima
       read(1,*) l_sym_centrale
       read(1,*) l_sym_axiale
    endif ! lfits
    ! ----------------------
    ! Dust global properties
    ! ----------------------
    if (lfits) then
       gas_dust=100.
       lstrat=.true. ; exp_strat=0.0 ; a_strat=1.0
       ldust_sublimation=.false.
       lchauff_int=.false. ; alpha=0.0
       T_min=1. ; T_max=1500. ; n_T=100
    else
       read(1,*)
       read(1,*)
       read(1,*) gas_dust
       read(1,*) lstrat, exp_strat, a_strat
       if (ldebris) then
          lstrat=.true.
       endif
       read(1,*) ldust_sublimation
       read(1,*) lchauff_int, alpha
       read(1,*) T_min, T_max, n_T
    endif  ! lfits
    ! ---------------
    ! Number of zones
    ! ---------------
    if (lfits) then
       n_zones=1
    else
       read(1,*)
       read(1,*)
       read(1,*) n_zones
       if (n_zones > 1) then
          lstrat=.true. ; exp_strat=0.
          write(*,*) "You are using a n-zone parameter file"
          write(*,*) "lstrat is set to true and exp_strat to 0."
       endif
    endif ! lfits
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! -----------------
    ! Density structure
    ! -----------------
    if (lfits) then
       do j=1,n_zones
          disk_zone(j)%geometry=1
          is_there_disk = .true.
          disk_zone(j)%diskmass=1.e-5
          disk_zone(j)%sclht=struct_file_zmax/cutoff
          disk_zone(j)%rref=struct_file_rref
          disk_zone(j)%rin=struct_file_rin ; disk_zone(j)%rout=struct_file_rout ; disk_zone(j)%edge=0.0
          disk_zone(j)%exp_beta=struct_file_beta
          disk_zone(j)%surf=0.0
       enddo ! n_zones
    else
       read(1,*)
       read(1,*)
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
             if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
                write(*,*) "Error : different values for size_neb"
                write(*,*) "Exiting"
                stop
             endif
          endif
       enddo ! n_zones
    endif ! lfits

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)

    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do

    ! ------
    ! Cavity
    ! ------
    if (lfits) then
       lcavity=.false.
       cavity%sclht=15. ; cavity%rref=50.
       cavity%exp_beta=1.5
    else
       read(1,*)
       read(1,*)
       read(1,*) lcavity
       read(1,*) cavity%sclht, cavity%rref
       read(1,*) cavity%exp_beta
    endif  !lfits

    ! ----------------
    ! Grain properties
    ! ----------------
    read(1,*)
    read(1,*)
    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.
    n_pop=0
    if (lfits) then
       do j=1, n_zones
          read(1,*) n_especes(j)
          if (n_especes(j) /= struct_file_nspecies) then
             write(*,*) "ERROR! Number of species in parameter file does not match structure of input FITS file"
             write(*,*) "Exiting."
             stop
          else
             somme=0.0
             do i=1, n_especes(j)
                n_pop = n_pop+1
                dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
                read(1,*,iostat=ios) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%porosity, dust_pop_tmp(n_pop)%frac_mass
                if (ios/=0) then
                   write(*,*) 'Error reading file: Incorrect number of lines in parameter file!'
                   write(*,*) 'Check the coherence of the number of species'
                   write(*,*) 'Exiting'
                   stop
                endif

                read(1,*) dust_pop_tmp(n_pop)%methode_chauffage
                read(1,*) dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amin=struct_file_amin*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%amax=struct_file_amax*dust_pop_tmp(n_pop)%sblow
                dust_pop_tmp(n_pop)%aexp=0.0;  dust_pop_tmp(n_pop)%n_grains=struct_file_n_grains
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
          endif
       enddo !n_zones
    else
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
    endif ! lfits

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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
    if (.not.lfits) then
       n_molecules = 1 ;
       allocate(mol(1))
       read(1,*)
       read(1,*)
       read(1,*) mol(1)%vmax_center_rt, vitesse_turb, mol(1)%n_speed_rt
       mol(1)%vmax_center_rt = mol(1)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
       vitesse_turb = vitesse_turb * 1.e3
       read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
       read(1,*) mol(1)%filename, mol(1)%iLevel_max
       read(1,*) mol(1)%lcst_abundance, mol(1)%abundance, mol(1)%abundance_file
       read(1,*) mol(1)%lline, mol(1)%nTrans_raytracing
       read(1,*) mol(1)%indice_Trans_rayTracing(1:mol(1)%nTrans_raytracing)
       mol(1)%n_speed_center_rt = mol(1)%n_speed_rt
       mol(1)%n_extraV_rt = 0 ; mol(1)%extra_deltaV_rt = 0.0
    endif ! lfits

    ! ---------------
    ! Star properties
    ! ---------------
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0

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
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * Rsun_to_AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV = 2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)


    close(unit=1)

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para210

  !**********************************************************************

  subroutine read_para209()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image

    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif

    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
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
    read(1,*) distance
    read(1,*) ang_disque
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) gas_dust
    read(1,*) lstrat, exp_strat, a_strat
    if (ldebris) then
       lstrat=.true.
    endif
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)
    read(1,*) n_zones
    if (n_zones > 1) then
       lstrat=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lstrat is set to true and exp_strat to 0."
    endif

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    read(1,*)
    read(1,*)

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
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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

    read(1,*)
    read(1,*)

    n_molecules = 1
    allocate(mol(n_molecules))
    read(1,*) mol(1)%vmax_center_rt, vitesse_turb, mol(1)%n_speed_rt
    mol(1)%vmax_center_rt = mol(1)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
    vitesse_turb = vitesse_turb * 1.e3 ! Conversion en m.s-1
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) mol(1)%filename, mol(1)%iLevel_max
    read(1,*) mol(1)%lcst_abundance, mol(1)%abundance, mol(1)%abundance_file
    read(1,*) mol(1)%lline, mol(1)%nTrans_raytracing
    read(1,*) mol(1)%indice_Trans_rayTracing(1:mol(1)%nTrans_raytracing)
    mol(1)%n_speed_center_rt = mol(1)%n_speed_rt
    mol(1)%n_extraV_rt = 0 ; mol(1)%extra_deltaV_rt = 0.0
    read(1,*)
    read(1,*)

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0



    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para209

  !**********************************************************************

  subroutine read_para208()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image

    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif

    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
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
    read(1,*) distance
    read(1,*) ang_disque
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) gas_dust
    read(1,*) lstrat, exp_strat, a_strat
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)
    read(1,*) n_zones
    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    read(1,*)
    read(1,*)

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
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    rmax = maxval(disk_zone(:)%rout)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.


    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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

    read(1,*)
    read(1,*)

    n_molecules = 1
    allocate(mol(n_molecules))
    read(1,*) mol(1)%vmax_center_rt, vitesse_turb, mol(1)%n_speed_rt
    mol(1)%vmax_center_rt = mol(1)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
    vitesse_turb = vitesse_turb * 1.e3
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) mol(1)%filename, mol(1)%iLevel_max
    read(1,*) mol(1)%lcst_abundance, mol(1)%abundance, mol(1)%abundance_file
    read(1,*) mol(1)%lline, mol(1)%nTrans_raytracing
    read(1,*) mol(1)%indice_Trans_rayTracing(1:mol(1)%nTrans_raytracing)
    mol(1)%n_speed_center_rt = mol(1)%n_speed_rt
    mol(1)%n_extraV_rt = 0 ; mol(1)%extra_deltaV_rt = 0.0
    read(1,*)
    read(1,*)

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0



    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para208

  !**********************************************************************

  subroutine read_para207()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) para_version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image

    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete = .false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif

    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
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
    read(1,*) distance
    ang_disque = 0.
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) gas_dust
    read(1,*) lstrat, exp_strat
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)
    read(1,*) n_zones
    if (n_zones > 1) then
       lstrat=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lstrat is set to true and exp_strat to 0."
    endif


    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    read(1,*)
    read(1,*)

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
          map_size=size_neb_tmp
       else
          if (abs(map_size-2.0*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.


    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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

    read(1,*)
    read(1,*)

    n_molecules = 1
    allocate(mol(1))
    read(1,*) mol(1)%vmax_center_rt, vitesse_turb, mol(1)%n_speed_rt
    mol(1)%vmax_center_rt = mol(1)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
    vitesse_turb = vitesse_turb * 1.e3
    read(1,*) lpop, lprecise_pop, lmol_LTE, largeur_profile
    read(1,*) mol(1)%filename, mol(1)%iLevel_max
    read(1,*) mol(1)%lcst_abundance, mol(1)%abundance, mol(1)%abundance_file
    read(1,*) mol(1)%lline, mol(1)%nTrans_raytracing
    read(1,*) mol(1)%indice_Trans_rayTracing(1:mol(1)%nTrans_raytracing)
    mol(1)%n_speed_center_rt = mol(1)%n_speed_rt
    mol(1)%n_extraV_rt = 0 ; mol(1)%extra_deltaV_rt = 0.0
    read(1,*)
    read(1,*)

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0



    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    a_strat = minval(dust_pop(:)%amin)

    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    return

  end subroutine read_para207

  !**********************************************************************

  subroutine read_para206()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real :: version
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image

    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete=.false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_pola) then
       n_Stokes = 4
       if (lsepar_contrib) then
          N_type_flux = 8
       else
          N_type_flux = 4
       endif
    else
       n_Stokes = 1
       if (lsepar_contrib) then
          N_type_flux = 5
       else
          N_type_flux = 1
       endif
    endif

    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
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
    read(1,*) distance
    ang_disque = 0.
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) gas_dust
    read(1,*) lstrat, exp_strat
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)
    read(1,*) n_zones
    if (n_zones > 1) then
       lstrat=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lstrat is set to true and exp_strat to 0."
    endif

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    read(1,*)
    read(1,*)

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
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.


    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .true.

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

    read(1,*)
    read(1,*)

    n_molecules = 1
    allocate(mol(1))
    read(1,*) mol(1)%vmax_center_rt, vitesse_turb, mol(1)%n_speed_rt
    mol(1)%vmax_center_rt = mol(1)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
    vitesse_turb = vitesse_turb * 1.e3
    read(1,*) lpop, lprecise_pop, mol(1)%lline, lmol_LTE, largeur_profile
    read(1,*) mol(1)%filename, mol(1)%iLevel_max
    read(1,*) mol(1)%lcst_abundance, mol(1)%abundance, mol(1)%abundance_file
    mol(1)%nTrans_raytracing = -1
    mol(1)%n_speed_center_rt = mol(1)%n_speed_rt
    mol(1)%n_extraV_rt = 0 ; mol(1)%extra_deltaV_rt = 0.0

    read(1,*)
    read(1,*)

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0



    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    a_strat = minval(dust_pop(:)%amin)


    return

  end subroutine read_para206

  !**********************************************************************

  subroutine read_para205()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real :: version
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image ; read(1,*) nbre_photons_spectre

    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete=.false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib, lsepar_pola

    if (lsepar_contrib) then
       if (lsepar_pola) then
          N_type_flux = 8
       else
          N_type_flux = 5
       endif
    else
       if (lsepar_pola) then
          N_type_flux = 4
       else
          N_type_flux = 1
       endif
    endif

    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) grid_type
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
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
    read(1,*) distance
    ang_disque =0.
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) lstrat, exp_strat
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)
    read(1,*) n_zones
    if (n_zones > 1) then
       lstrat=.true. ; exp_strat=0.
       write(*,*) "You are using a n-zone parameter file"
       write(*,*) "lstrat is set to true and exp_strat to 0."
    endif

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    read(1,*)
    read(1,*)

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
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.


    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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

    read(1,*)
    read(1,*)

    n_molecules = 1
    allocate(mol(1))
    read(1,*) mol(1)%vmax_center_rt, vitesse_turb, mol(1)%n_speed_rt
    mol(1)%vmax_center_rt = mol(1)%vmax_center_rt * 1.e3 ! Conversion en m.s-1
    vitesse_turb = vitesse_turb * 1.e3
    read(1,*) lpop, lprecise_pop, mol(1)%lline, lmol_LTE, largeur_profile
    read(1,*) mol(1)%filename, mol(1)%iLevel_max
    read(1,*) gas_dust, mol(1)%abundance
    mol(1)%nTrans_raytracing=-1
    mol(1)%lcst_abundance = .true.
    mol(1)%n_speed_center_rt = mol(1)%n_speed_rt
    mol(1)%n_extraV_rt = 0 ; mol(1)%extra_deltaV_rt = 0.0

    read(1,*)
    read(1,*)

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0



    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    a_strat = minval(dust_pop(:)%amin)

    return

  end subroutine read_para205

  !**********************************************************************

  subroutine read_para204()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real :: version
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    logical :: ljunk

    dust_pop_tmp(:)%type = "Mie"

    grid_type = 1
    n_zones = 1
    lsepar_pola = .true.
    lonly_capt_interet = .false.

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')
    read(1,*) version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;
    read(1,*) nbre_photons_image ; read(1,*) nbre_photons_spectre
    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete=.false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib
    if (lsepar_contrib) then
       N_type_flux = 8
    else
       N_type_flux = 5
    endif
    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    N_incl = N_thet ; capt_debut=1 ; capt_fin=N_thet
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    read(1,*) capt_interet, delta_capt, angle_interet
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    read(1,*) distance
    ang_disque = 0.
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) lstrat, exp_strat
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)

    do j=1,n_zones
       disk_zone(j)%geometry = 1
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf
       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.


    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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

    read(1,*)
    read(1,*)

    n_molecules = 1
    allocate(mol(1))
    read(1,*) lmol_LTE, ljunk, ljunk, mol(1)%n_speed_rt
    read(1,*) mol(1)%filename
    read(1,*) gas_dust, mol(1)%abundance

    read(1,*)
    read(1,*)

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0



    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    a_strat = minval(dust_pop(:)%amin)

    return

  end subroutine read_para204

  !**********************************************************************

  subroutine read_para203()

    implicit none

    integer :: i, j, alloc_status, ios, ind_pop
    real :: version
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    grid_type = 1
    gas_dust=100
    lsepar_pola = .true.
    lonly_capt_interet = .false. ;

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0

    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*) version
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;  read(1,*) nbre_photons_image
    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda, lambda_min, lambda_max
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete=.false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib
    if (lsepar_contrib) then
       N_type_flux = 8
    else
       N_type_flux = 5
    endif
    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    read(1,*) capt_interet, delta_capt, angle_interet
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    N_incl = N_thet ; capt_debut=1 ; capt_fin=N_thet
    read(1,*) distance
    ang_disque = 0.
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) lstrat, exp_strat
    read(1,*) ldust_sublimation
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)

    do j=1,n_zones
       disk_zone(j)%geometry = 1
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf
       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)

    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    lRE_LTE=.false. ; lRE_nLTE=.false. ; lnRE=.false.


    n_pop=0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme=0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
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

    if (lRE_LTE.and.lRE_nLTE) then
       write(*,*) "Error : cannot mix grains in LTE and nLTE"
       write(*,*) " Is it usefull anyway ???"
       write(*,*) "Exiting"
       stop
    endif

    ! variables triees
    allocate(dust_pop(n_pop), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
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

    read(1,*)
    read(1,*)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0


    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)


    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf


    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    a_strat = minval(dust_pop(:)%amin)

    return

  end subroutine read_para203

  !**********************************************************************

  subroutine read_para202()

    implicit none

    integer :: i, alloc_status, ios, j, ind_pop
    real(kind=db) :: size_neb_tmp, somme
    real :: gas_dust

    type(dust_pop_type), dimension(100) :: dust_pop_tmp
    integer, dimension(100) :: n_especes

    dust_pop_tmp(:)%type = "Mie"

    grid_type = 1
    n_zones = 1
    gas_dust=100
    lsepar_pola = .true.
    lonly_capt_interet = .false. ;
    lambda_min = 0.1 ; lambda_max = 3000.0
    lnRE=.false.

    ! Allocation des variables pour disque a une zone
    allocate(disk_zone(n_zones), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error disk parameters'
       stop
    endif
    disk_zone%exp_beta=0.0; disk_zone%surf=0.0; disk_zone%sclht=0.0; disk_zone%diskmass=0.0; disk_zone%rref=0.0
    disk_zone%rin=0.0 ; disk_zone%rout=0.0 ; disk_zone%edge=0.0


    ! Lecture du fichier de parametres
    open(unit=1, file=para, status='old')

    read(1,*)
    read(1,*)
    read(1,*)
    read(1,*) nbre_photons_loop ;  read(1,*) nbre_photons_eq_th ; read(1,*) nbre_photons_lambda ;  read(1,*) nbre_photons_image
    read(1,*) !lcheckpoint,  checkpoint_time
    read(1,*)
    read(1,*)
    read(1,*) n_lambda
    lmono0 = (n_lambda==1) ; lmono = lmono0
    read(1,*) ltemp, lsed, lsed_complete
    if (.not.lsed) lsed_complete=.false.
    read(1,*) tab_wavelength
    read(1,*) l_em_disk_image
    read(1,*) lsepar_contrib
    if (lsepar_contrib) then
       N_type_flux = 8
    else
       N_type_flux = 5
    endif
    read(1,*) tau_seuil, wl_seuil
    read(1,*)
    read(1,*)
    read(1,*) n_rad, nz, n_az, n_rad_in
    if ((.not.l3D).and.(n_az > 1)) then
       write(*,*) "WARNING : n_az > 1 in 2D configuration, forcing n_az=1"
       n_az=1
    endif
    n_cell_max = nz *n_rad
    read(1,*)
    read(1,*)
    read(1,*,iostat=ios) N_thet, N_phi, igridx, igridy, zoom, n_cartes
    if (ios/=0) then
       n_cartes=1
    endif
    if (n_cartes > n_cartes_max) then
       write(*,*) "Erreur : n_cartes >", n_cartes_max
       stop
    endif
    allocate(igridx2(n_cartes-1), igridy2(n_cartes-1), maxigrid2(n_cartes-1))
    maxigrid = max(igridx, igridy)
    do i=1, n_cartes-1
       read(1,*) igridx2(i), igridy2(i)
       maxigrid2(i) = max(igridx2(i), igridy2(i))
    enddo
    read(1,*) capt_interet, delta_capt, angle_interet
    capt_inf=max(1,capt_interet-delta_capt)
    capt_sup=min(N_thet,capt_interet+delta_capt)
    N_incl = N_thet ; capt_debut=1 ; capt_fin=N_thet
    read(1,*) distance
    ang_disque = 0.
    read(1,*)
    read(1,*)
    read(1,*) scattering_method
    read(1,*) aniso_method
    read(1,*)
    read(1,*)
    read(1,*) l_sym_ima
    read(1,*) l_sym_centrale
    read(1,*) l_sym_axiale
    read(1,*)
    read(1,*)
    read(1,*) lstrat, exp_strat
    read(1,*) ldust_sublimation
    read(1,*) lRE_LTE
    lRE_nLTE = .not.lRE_LTE
    read(1,*) lchauff_int, alpha
    read(1,*) T_min, T_max, n_T
    read(1,*)
    read(1,*)

    do j=1,n_zones
       disk_zone(j)%geometry = 1
       read(1,*) disk_zone(j)%diskmass
       read(1,*) disk_zone(j)%sclht, disk_zone(j)%rref
       read(1,*) disk_zone(j)%rin, disk_zone(j)%rout, size_neb_tmp, disk_zone(j)%edge
       read(1,*) disk_zone(j)%exp_beta
       read(1,*) disk_zone(j)%surf
       if (j==1) then
          map_size=2*size_neb_tmp
       else
          if (abs(map_size-2*size_neb_tmp) > 1.e-6*map_size) then
             write(*,*) "Error : different values for size_neb"
             write(*,*) "Exiting"
             stop
          endif
       endif
    enddo

    disk_zone(:)%gas_to_dust = gas_dust
    disk_zone(:)%rmin = disk_zone(:)%rin - 5*disk_zone(:)%edge
    rmin = minval(disk_zone(:)%rmin)
    disk_zone(:)%Rmax = disk_zone(:)%rout
    Rmax = maxval(disk_zone(:)%Rmax)
    diskmass = sum(disk_zone(:)%diskmass)


    if (rmin < 0.0) then
       write(*,*) "Error : r_min < 0.0"
       stop
    endif

    if (igridx == igridy) then
       deltapix_x = 1
       deltapix_y = 1
    else if (igridx > igridy) then
       deltapix_x = 1
       deltapix_y = 1 - (igridx/2) + (igridy/2)
    else
       deltapix_x = 1 - (igridy/2) + (igridx/2)
       deltapix_y = 1
    endif
    size_pix=maxigrid/(map_size)
    allocate(deltapix_x2(n_cartes-1),deltapix_y2(n_cartes-1),size_pix2(n_cartes-1))
    do i=1,n_cartes-1
       if (igridx2(i) == igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1
       else if (igridx2(i) > igridy2(i)) then
          deltapix_x2(i) = 1
          deltapix_y2(i) = 1 - (igridx2(i)/2) + (igridy2(i)/2)
       else
          deltapix_x2(i) = 1 - (igridy2(i)/2) + (igridx2(i)/2)
          deltapix_y2(i) = 1
       endif
       size_pix2(i)=maxigrid2(i)/(map_size)
    end do
    read(1,*)
    read(1,*)

    n_pop = 0
    do j=1, n_zones
       read(1,*) n_especes(j)
       somme = 0.0
       do i=1, n_especes(j)
          n_pop = n_pop+1
          dust_pop_tmp(n_pop)%n_components = 1 ; dust_pop_tmp(n_pop)%component_volume_fraction(1) = 1.0
          read(1,*) dust_pop_tmp(n_pop)%indices(1), dust_pop_tmp(n_pop)%component_rho1g(1), dust_pop_tmp(n_pop)%frac_mass
          read(1,*) dust_pop_tmp(n_pop)%amin, dust_pop_tmp(n_pop)%amax, dust_pop_tmp(n_pop)%aexp, dust_pop_tmp(n_pop)%n_grains
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
    if (alloc_status > 0) then
       write(*,*) 'Allocation error n_pop tmp'
       stop
    endif
    dust_pop%is_PAH = .false.
    dust_pop%is_opacity_file = .false.

    ! Classement des populations de grains : LTE puis nLTE puis nRE
    ind_pop = 0
    grain_RE_LTE_start = 1
    grain_RE_LTE_end = 0
    if (lRE_LTE) then
       do i=1, n_pop
          ind_pop=ind_pop+1
          dust_pop(ind_pop) = dust_pop_tmp(i)
          dust_pop(ind_pop)%ind_debut = grain_RE_LTE_end + 1
          dust_pop(ind_pop)%ind_fin = grain_RE_LTE_end + dust_pop(ind_pop)%n_grains
          dust_pop(ind_pop)%methode_chauffage = 1
          grain_RE_LTE_end = grain_RE_LTE_end +  dust_pop(ind_pop)%n_grains
       enddo
    endif

    grain_RE_nLTE_start = grain_RE_LTE_end + 1
    grain_RE_nLTE_end =  grain_RE_LTE_end
    if (lRE_nLTE) then
       do i=1, n_pop
          ind_pop=ind_pop+1
          dust_pop(ind_pop) = dust_pop_tmp(i)
          dust_pop(ind_pop)%ind_debut = grain_RE_nLTE_end + 1
          dust_pop(ind_pop)%ind_fin = grain_RE_nLTE_end + dust_pop(ind_pop)%n_grains
          dust_pop(ind_pop)%methode_chauffage = 2
          grain_RE_nLTE_end = grain_RE_nLTE_end +  dust_pop(ind_pop)%n_grains
       enddo
    endif

    grain_nRE_start = grain_RE_nLTE_end + 1
    grain_nRE_end =  grain_RE_nLTE_end

    n_grains_RE_LTE = grain_RE_LTE_end  - grain_RE_LTE_start + 1
    n_grains_RE_nLTE = grain_RE_nLTE_end  - grain_RE_nLTE_start + 1
    n_grains_nRE = grain_nRE_end  - grain_nRE_start + 1

    n_grains_tot=sum(dust_pop(:)%n_grains)

    read(1,*)
    read(1,*)

    if (ldust_sublimation) then
       if (n_lambda==1) then
          write(*,*) "error : sub radius"
       endif
    endif

    read(1,*,iostat=ios) n_etoiles
    if (ios/=0) then
       write(*,*) 'Error reading file: you are using a 2-zone disk parameter file'
       write(*,*) 'You must use the [-2zone] option to calculate a 2-zone disk'
       !   write(*,*) ' '
       write(*,*) 'Exiting'
       stop
    endif
    allocate(etoile(n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error etoile'
       stop
    endif

    allocate(CDF_E_star(n_lambda,0:n_etoiles), prob_E_star(n_lambda,n_etoiles), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error CDF_E_star'
       stop
    endif
    CDF_E_star = 0.0 ; prob_E_star = 0.0


    do i=1,n_etoiles
       read(1,*) etoile(i)%T, etoile(i)%r, etoile(i)%M, etoile(i)%x, etoile(i)%y, etoile(i)%z, etoile(i)%lb_body
       if (.not.etoile(i)%lb_body) then
          read(1,*) etoile(i)%spectre
       else
          read(1,*)
       endif
       ! Passage rayon en AU
       etoile(i)%r = etoile(i)%r * 0.00466666666 ! 1 rayon solaire en AU
    enddo

    etoile(:)%fUV = 0.0
    etoile(:)%slope_UV =  2.2 - 2.0  ! Fnu -> F_lambda, par defaut pour version < 2.12 (y compris DENT)

    close(unit=1)

    !  write(*,*) amin ; write(*,*) amax ; write(*,*) aexp
    !  write(*,*) exp_strat ; write(*,*) diskmass ; write(*,*) sclht
    !  write(*,*) rout ; write(*,*) rin ; write(*,*) exp_beta ; write(*,*) surf

    nbre_photons_lim = nbre_photons_lim*real(N_thet)*real(N_phi)*real(nbre_photons_lambda)

    a_strat = minval(dust_pop(:)%amin)

    return

  end subroutine read_para202

  !**********************************************************************

end module read_params
