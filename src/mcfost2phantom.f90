module mcfost2phantom

  implicit none

contains

  subroutine init_mcfost_phantom(mcfost_para_filename, ndusttypes, use_SPH_limits_file, SPH_limits_file, SPH_limits, ierr, &
       keep_particles, fix_star)

    use parametres
    use init_mcfost, only : set_default_variables, get_mcfost_utils_dir
    use read_params, only : read_para
    use dust_prop, only : build_grain_size_distribution, init_indices_optiques, prop_grains
    use grid, only : define_physical_zones, order_zones
    use wavelengths, only : init_lambda, n_lambda
    use optical_depth, only : no_dark_zone
    use mem, only : alloc_dust_prop
    use SPH2mcfost, only : read_SPH_limits_file
    use messages, only : error

    character(len=*), intent(in) :: mcfost_para_filename, SPH_limits_file
    integer, intent(in) :: ndusttypes
    logical, intent(in) :: use_SPH_limits_file
    real(dp), dimension(6), intent(out) :: SPH_limits
    integer, intent(out) :: ierr
    real, intent(in), optional :: keep_particles
    logical, intent(in), optional :: fix_star

    integer, target :: lambda, lambda0
    integer, pointer, save :: p_lambda

    write(*,*)
    write(*,*) "------------------------------"
    write(*,*) "Initializing MCFOST library"
    write(*,*) "------------------------------"
    write(*,*)

    ! Global logical variables
    call set_default_variables()

    if (present(keep_particles)) then
       SPH_keep_particles = keep_particles
       write(*,*) "WARNING: updating SPH_keep_particles to" , SPH_keep_particles
    endif

    if (present(fix_star)) then
       lfix_star = .true.
       write(*,*) "WARNING: using mcfost parameters for the stars"
    endif

    ! Model limits
    if (use_SPH_limits_file) then
       call read_SPH_limits_file(SPH_limits_file, SPH_limits)
    else
       SPH_limits = 0.
    endif

    write(*,*) "WARNING: internal heating is turned off in mcfost"
    lno_internal_energy = .true.

    ! Looking for the mcfost utils directory
    call get_mcfost_utils_dir()

    ! parameter file
    call read_para(mcfost_para_filename)

    ! Setting option for the mcfost2phantom interface
    ltemp = .true. ; lsed = .false. ; lsed_complete = .false.
    lVoronoi = .true. ; l3D = .true.
    lsepar_pola = .false.
    lvariable_dust = (ndusttypes > 0)

    ! We do not use a file with limits yet
    llimits_file = .false.

    ! A few modes are not available via the interface
    if (lnRE) then
       call error("Non-equilibrium grains are not yet implemented in libmcfost",ierr=ierr)
       if (ierr /= 0) return
    endif

    if (n_zones > 1) then
       call error("mcfost2phantom only works with a 1zone parameter file",ierr=ierr)
       if (ierr /= 0) return
    endif

    ! Building the wavelength & basic dust properties grid
    call init_lambda()
    call init_indices_optiques()

    ! Building the dust grain population
    call build_grain_size_distribution()

    ! Building the model volume and corresponding grid
    call order_zones()
    call define_physical_zones()

    if (lscattering_method1) then
       lambda = 1
       p_lambda => lambda
    else
       if (p_n_lambda_pos == n_lambda) then
          lambda = 1
          p_lambda => lambda
       else
          lambda0 = 1
          p_lambda => lambda0
       endif
    endif

    ! Dust properties
    write(*,'(a30, $)') "Computing dust properties ..."
    call alloc_dust_prop()
    do lambda=1,n_lambda
       call prop_grains(lambda)
    enddo !n
    write(*,*) "Done"

    lonly_LTE = .false.
    lonly_nLTE = .false.
    if (lRE_LTE .and. .not.lRE_nLTE .and. .not. lnRE) lonly_LTE = .true.
    if (lRE_nLTE .and. .not.lRE_LTE .and. .not. lnRE) lonly_nLTE = .true.

    if (aniso_method==1) then
       lmethod_aniso1 = .true.
    else
       lmethod_aniso1 = .false.
    endif

    ierr = 0
    return

  end subroutine init_mcfost_phantom

  !*************************************************************************

  subroutine run_mcfost_phantom(np,nptmass,ntypes,ndusttypes,dustfluidtype,&
    npoftype,xyzh,vxyzu,iphase,grainsize,graindens,dustfrac,massoftype,&
    xyzmh_ptmass,hfact,umass,utime,udist,ndudt,dudt,compute_Frad,SPH_limits,&
    Tphantom,Frad,mu_gas,ierr,write_T_files,ISM,T_to_u)

    use parametres
    use constantes, only : mu
    use read_phantom
    use stars, only : E_ISM
    use radiation_field, only : reset_radiation_field
    use thermal_emission, only : select_wl_em, repartition_energie, init_reemission, &
         temp_finale, temp_finale_nlte, repartition_wl_em, set_min_temperature, reset_temperature, E_abs_nRE
    use mem, only : alloc_dynamique, deallocate_densities
    use naleat, only : seed, stream, gtype
    use SPH2mcfost, only : SPH_to_Voronoi, compute_stellar_parameters
    use Voronoi_grid, only : Voronoi, deallocate_Voronoi
    use dust_transfer, only : emit_packet, propagate_packet
    use utils, only : progress_bar
    use stars, only : repartition_energie_etoiles, repartition_energie_ISM
    use grid,only : setup_grid
    use scattering, only : setup_scattering
    use optical_depth, only : no_dark_zone, integ_tau, opacite
    use wavelengths, only : n_lambda, tab_lambda
    use Temperature, only : Tdust
    use mcfost_env
    !$ use omp_lib

#include "sprng_f.h"

    integer, intent(in) :: np, nptmass, ntypes,ndusttypes,dustfluidtype
    real(dp), dimension(4,np), intent(in) :: xyzh,vxyzu
    integer(kind=1), dimension(np), intent(in) :: iphase
    real(dp), dimension(ndusttypes,np), intent(in) :: dustfrac
    real(dp), dimension(ndusttypes), intent(in) :: grainsize, graindens
    real(dp), dimension(ntypes), intent(in) :: massoftype
    real(dp), intent(in) :: hfact, umass, utime, udist, T_to_u
    real(dp), dimension(:,:), intent(in) :: xyzmh_ptmass
    integer, dimension(ntypes), intent(in) :: npoftype

    integer, parameter :: n_files = 1 ! the library only works on 1 set of phantom particles
    integer(kind=1), dimension(np) :: ifiles

    logical, intent(in), optional :: write_T_files

    logical, intent(in) :: compute_Frad ! does mcfost need to compute the radiation pressure
    real(dp), dimension(6), intent(in) :: SPH_limits ! not used yet, llimits_file is set to false

    integer, intent(in) :: ndudt
    real(dp), dimension(ndudt), intent(in) :: dudt

    integer, intent(in) :: ISM ! ISM heating: 0 -> no ISM radiation field, 1 -> ProDiMo, 2 -> Bate & Keto

    real(sp), dimension(np), intent(out) :: Tphantom ! mcfost stores Tdust as real, not dp
    real(sp), dimension(3,ndusttypes,np), intent(out) :: Frad
    real(dp), intent(out) :: mu_gas
    integer, intent(out) :: ierr

    real, parameter :: Tmin = 1.

    real(dp), dimension(:), allocatable :: x_SPH,y_SPH,z_SPH,h_SPH,rhogas, massgas, vx_SPH,vy_SPH,vz_SPH, SPH_grainsizes
    integer, dimension(:), allocatable :: particle_id
    real(dp), dimension(:,:), allocatable :: rhodust, massdust
    real, dimension(:), allocatable :: extra_heating
    real(dp), dimension(n_files,ntypes) :: massoftype2

    real(kind=dp), dimension(4) :: Stokes
    real(kind=dp) :: nnfot2
    real(kind=dp) :: x,y,z, u,v,w
    real :: rand, time, cpu_time_begin, cpu_time_end
    integer :: n_SPH, icell, nbre_phot2, ibar, id, nnfot1_cumul, i_SPH, i, lambda_seuil
    integer :: itime !, time_begin, time_end, time_tick, time_max
    logical :: lpacket_alive, lintersect, laffichage, flag_star, flag_scatt, flag_ISM
    integer, target :: lambda, lambda0
    integer, pointer, save :: p_lambda

    logical, save :: lfirst_time = .true.

    integer :: i_Phantom

    ! We use the phantom_2_mcfost interface with 1 file
    ifiles(:) = 1 ; massoftype2(1,:) = massoftype(:)

    write(*,*)
    write(*,*) "------------------------------"
    write(*,*) "Running MCFOST via the library"
    write(*,*) "------------------------------"
    write(*,*)

    ! debut de l'execution
    call system_clock(time_begin,count_rate=time_tick,count_max=time_max)
    call cpu_time(cpu_time_begin)

    ierr = 0
    mu_gas = mu ! Molecular weight
    Frad = 0.

    call phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,n_files,dustfluidtype,xyzh,&
         vxyzu,iphase,grainsize,dustfrac(1:ndusttypes,np),massoftype2(1,1:ntypes),xyzmh_ptmass,hfact,&
         umass,utime,udist,graindens,ndudt,dudt,ifiles,&
         n_SPH,x_SPH,y_SPH,z_SPH,h_SPH,vx_SPH,vy_SPH,vz_SPH,particle_id,&
         SPH_grainsizes,massgas,massdust,rhogas,rhodust,extra_heating,T_to_u)

    if (.not.lfix_star) call compute_stellar_parameters()

    ! Performing the Voronoi tesselation & defining density arrays
    call SPH_to_Voronoi(n_SPH, ndusttypes, x_SPH,y_SPH,z_SPH,h_SPH,vx_SPH,vy_SPH,vz_SPH, &
         massgas,massdust,rhogas,rhodust,SPH_grainsizes, SPH_limits, .false.)

    call setup_grid()
    call setup_scattering()

    ! We allocate the total number of SPH cells as the number of Voronoi cells mays vary
    if (lfirst_time) then
       call alloc_dynamique(n_cells_max = n_SPH + n_etoiles)
       lfirst_time = .false.
    endif
    call no_dark_zone()
    lapprox_diffusion=.false.

    ! init random number generator
    stream = 0.0
    do i=1, nb_proc
       stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    enddo

    if (lscattering_method1) then
       lambda = 1
       p_lambda => lambda
    else
       if (p_n_lambda_pos == n_lambda) then
          lambda = 1
          p_lambda => lambda
       else
          lambda0 = 1
          p_lambda => lambda0
       endif
    endif

    ! ToDo : needs to be made parallel
    do lambda=1,n_lambda
       call opacite(lambda, p_lambda)
    enddo !n

    call repartition_energie_etoiles()
    call repartition_energie_ISM(ISM)

    ! ToDo : needs to be made parallel
    call init_reemission(lextra_heating,extra_heating)

    !$omp parallel default(none) private(lambda) shared(n_lambda)
    !$omp do schedule(static,1)
    do lambda=1, n_lambda
       call repartition_energie(lambda)
    enddo
    !$omp end do
    !$omp end parallel

    call repartition_wl_em()

    letape_th = .true.
    laffichage=.true.
    nbre_phot2 = nbre_photons_eq_th

    test_tau : do lambda=1,n_lambda
       if (tab_lambda(lambda) > wl_seuil) then
          lambda_seuil=lambda
          exit test_tau
       endif
    enddo test_tau
    write(*,*) "lambda =", real(tab_lambda(lambda_seuil))
    call integ_tau(lambda_seuil)

    write(*,*) "Computing temperature structure ..."
    ! Making the MC run
    if (laffichage) call progress_bar(0)
    !$omp parallel &
    !$omp default(none) &
    !$omp firstprivate(lambda,p_lambda) &
    !$omp private(id,icell,lpacket_alive,lintersect,nnfot2,rand) &
    !$omp private(x,y,z,u,v,w,Stokes,flag_star,flag_ISM,flag_scatt) &
    !$omp shared(nbre_photons_loop) &
    !$omp shared(nbre_phot2,nb_proc) &
    !$omp shared(stream,laffichage,nbre_photons_lambda, nnfot1_cumul,ibar) &
    !$omp reduction(+:E_abs_nRE)
    E_abs_nRE = 0.0

    id = 1 ! Pour code sequentiel
    !$ id = omp_get_thread_num() + 1
    ibar=1 ;  nnfot1_cumul = 0

    !$omp do schedule(dynamic,1)
    do nnfot1=1,nbre_photons_loop
       nnfot2 = 0.0_dp
       photon : do while ((nnfot2 < nbre_phot2))
          nnfot2=nnfot2+1.0_dp

          ! Choix longueur d'onde
          rand = sprng(stream(id))
          call select_wl_em(rand,lambda)

          ! Emission du paquet
          call emit_packet(id,lambda, icell,x,y,z,u,v,w,stokes,flag_star,flag_ISM,lintersect)
          Stokes = 0.0_dp ; Stokes(1) = 1.0_dp
          lpacket_alive = .true.

          ! Propagation du packet
          if (lintersect) call propagate_packet(id,lambda,p_lambda,icell,x,y,z,u,v,w,stokes, &
               flag_star,flag_ISM,flag_scatt,lpacket_alive)
       enddo photon !nnfot2

       ! Progress bar
       !$omp atomic
       nnfot1_cumul = nnfot1_cumul+1
       if (laffichage) then
          if (real(nnfot1_cumul) > 0.02*ibar * real(nbre_photons_loop)) then
             call progress_bar(ibar)
             !$omp atomic
             ibar = ibar+1
          endif
       endif
    enddo !nnfot1
    !$omp end do
    !$omp end parallel
    if (laffichage) call progress_bar(50)

    if (lRE_LTE) then
       call Temp_finale()
    end if
    if (lRE_nLTE) then
       call Temp_finale_nLTE()
    endif
    if (lnRE) then
       ! TBD
    endif

    call set_min_Temperature(Tmin)
    Tphantom = -1.0 ;

    if (present(write_T_files)) then
       if (write_T_files) call write_mcfost2phantom_temperature()
    endif

    do icell=1, n_cells
       i_SPH = Voronoi(icell)%id
       if (i_SPH > 0) then
          i_Phantom = particle_id(i_SPH)
          Tphantom(i_Phantom) = Tdust(icell)
       endif
    enddo

    call deallocate_Voronoi()
    call deallocate_densities()

    ! reset energy and temperature arrays
    call reset_radiation_field()
    call reset_temperature()

    ! Temps d'execution
    call system_clock(time_end)
    if (time_end < time_begin) then
       time=(time_end + (1.0 * time_max)- time_begin)/real(time_tick)
    else
       time=(time_end - time_begin)/real(time_tick)
    endif
    if (time > 60) then
       itime = real(time)
       write (*,'(" Processing complete in ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
    else
       itime = int(time)
       write (*,'(" Processing complete in ", F5.2, "s")')  time
    endif
    call cpu_time(cpu_time_end)
    time = cpu_time_end - cpu_time_begin
    if (time > 60) then
       itime = int(time)
       write (*,'(" CPU time used          ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
    else
       write (*,'(" CPU time used          ", F5.2, "s")')  time
    endif

    ! Verifications et codes d'erreur
    if (maxval(Tphantom) < 0.) then
       write(*,*) "***********************************************"
       write(*,*) "ERROR : PB setting T_DUST", n_cells
       write(*,*)  minval(Tphantom), maxval(Tphantom)
       write(*,*) "***********************************************"

       ierr = 1
       return
    endif

    write(*,*)
    write(*,*) "------------------------------"
    write(*,*) "End of MCFOST run"
    write(*,*) "------------------------------"
    write(*,*)

    return

  end subroutine run_mcfost_phantom

!*************************************************************************

  subroutine write_mcfost2phantom_temperature()

    use mcfost_env, only : data_dir
    use utils, only : appel_syst
    use output, only : ecriture_temperature, write_disk_struct
    use messages, only : error

    integer, save :: n_call = -1 ! to match phantom's dump numbers
    character(len=1) :: s
    integer :: syst_status

    n_call = n_call + 1
    if (n_call > 9) call error("STOPPING to write T files at dump #9")

    syst_status = 0

    write(s,'(i1)') n_call ; data_dir = "mcfost_"//s
    write(*,*)
    write(*,*) "**************************************************"
    write(*,*) "Saving Temperature structure in "//trim(data_dir)
    write(*,*) "**************************************************"
    write(*,*)

    call appel_syst("rm -rf "//trim(data_dir)//" ; mkdir -p "//trim(data_dir), syst_status)
    call ecriture_temperature(1)

    call appel_syst("rm -rf data_disk ; mkdir -p data_disk", syst_status)
    call write_disk_struct(.false.)
    call appel_syst("mv data_disk/grid.fits.gz "//trim(data_dir), syst_status)

    return

  end subroutine write_mcfost2phantom_temperature

end module mcfost2phantom
