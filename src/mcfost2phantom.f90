module mcfost2phantom

!  use parametres

  implicit none

contains

  subroutine init_mcfost_phantom(mcfost_para_filename, ierr) !,  np, nptmass, ntypes, ndusttypes, npoftype)

    use parametres
    use init_mcfost, only : set_default_variables, get_mcfost_utils_dir
    use read_params, only : read_para
    use disk, only : n_zones
    use dust_transfer, only : transfert_poussiere

    character(len=*), intent(in) :: mcfost_para_filename
    integer, intent(out) :: ierr

    integer, target :: lambda, lambda0
    integer, pointer, save :: p_lambda

    ! Global logical variables
    call set_default_variables()

    ! Looking for the mcfost utils directory
    call get_mcfost_utils_dir()

    ! parameter file
    call read_para(mcfost_para_filename)

    if (lnRE) then
       write(*,*) "Non-equilibrium grains are not yet implemented in libmcfost"
       write(*,*) "Exiting"
       stop
    endif

    ! Setting option for the mcfost2phantom interface
    ltemp = .true. ; lsed = .false. ; lsed_complete = .false.
    lVoronoi = .true. ; l3D = .true.

    ! We do not use a file with limits yet
    llimits_file = .false.

    if (n_zones > 1) then
       write(*,*) "ERROR: mcfost2phantom only works with a 1zone parameter file"
       write(*,*) "Exiting"
       ierr = 1
       return
    endif

    ! Building the wavelength & basic dust properties grid
    call init_lambda()
    call init_indices_optiques()

    ! Building the dust grain population
    call build_grain_size_distribution()

    ! Building the model volume and corresponding grid
    call order_zones()
    call define_physical_zones()

    ! Dust properties
    write(*,'(a30, $)') "Computing dust properties ..."
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

    do lambda=1,n_lambda
       call prop_grains(lambda)
       call opacite(lambda, p_lambda)
    enddo !n
    write(*,*) "Done"

    call no_dark_zone()
    lapprox_diffusion=.false.

    ierr = 0
    return

  end subroutine init_mcfost_phantom

  !*************************************************************************

  subroutine run_mcfost_phantom(np,nptmass,ntypes,ndusttypes,dustfluidtype,npoftype,xyzh,iphase,grainsize,&
       dustfrac, massoftype,xyzmh_ptmass,hfact,umass,utime,udist,graindens, &
       compute_Frad,SPH_limits_file, & ! options
       Tdust,Frad,mu_gas,ierr) ! intent(out)

    use parametres
    use constantes, only : mu
    use read_phantom
    use prop_star, only : n_etoiles
    use em_th, only : temperature, E_abs_nRE, frac_E_stars
    use thermal_emission, only : reset_radiation_field
    use mem, only : alloc_dynamique
    use naleat, only : seed, stream, gtype
    use SPH2mcfost, only : SPH_to_Voronoi, compute_stellar_parameters
    use Voronoi_grid, only : Voronoi
    !$ use omp_lib

#include "sprng_f.h"

    integer, intent(in) :: np, nptmass, ntypes,ndusttypes,dustfluidtype
    real(db), dimension(4,np), intent(in) :: xyzh
    integer(kind=1), dimension(np), intent(in) :: iphase
    real(db), dimension(ndusttypes,np), intent(in) :: dustfrac
    real(db), dimension(ndusttypes), intent(in) :: grainsize
    real(db), dimension(ntypes), intent(in) :: massoftype
    real(db), intent(in) :: hfact, umass, utime, udist, graindens
    real(db), dimension(:,:), intent(in) :: xyzmh_ptmass
    integer, dimension(ntypes), intent(in) :: npoftype

    logical, intent(in) :: compute_Frad ! does mcfost need to compute the radiation pressure
    character(len=512), intent(in) :: SPH_limits_file ! not used yet, llimits_file is set to false


    real, dimension(np), intent(out) :: Tdust ! mcfost stores Tdust as real, not db
    real, dimension(3,ndusttypes,np), intent(out) :: Frad
    real(db), intent(out) :: mu_gas
    integer, intent(out) :: ierr


    real(db), dimension(:), allocatable :: XX,YY,ZZ,rhogas, massgas
    real(db), dimension(:,:), allocatable :: rhodust, massdust

    real(kind=db), dimension(4) :: Stokes
    real(kind=db) :: nnfot2
    real(kind=db) :: x,y,z, u,v,w
    real :: rand, time, cpu_time_begin, cpu_time_end
    integer :: ncells, n_SPH, icell, nbre_phot2, ibar, id, nnfot1_cumul, i_SPH, i
    integer :: itime !, time_begin, time_end, time_tick, time_max
    logical :: lpacket_alive, lintersect, laffichage, flag_star, flag_scatt, flag_ISM
    integer, target :: lambda, lambda0
    integer, pointer, save :: p_lambda

    logical, save :: lfirst_time


    ! debut de l'execution
    call system_clock(time_begin,count_rate=time_tick,count_max=time_max)
    call cpu_time(cpu_time_begin)

    ierr = 0
    mu_gas = mu ! Molecular weight
    Frad = 0.

    call phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,dustfluidtype,xyzh,iphase,grainsize,dustfrac,&
         massoftype(1:ntypes),xyzmh_ptmass,hfact,umass,utime,udist,graindens,XX,YY,ZZ,massgas,massdust,rhogas,rhodust,n_SPH)
    if (ncells <= 0) then
       ierr = 1
       return
    endif

    call compute_stellar_parameters()

    ! Performing the Voronoi tesselation & defining density arrays
    call SPH_to_Voronoi(n_SPH, ndusttypes, XX,YY,ZZ,massgas,massdust,rhogas,rhodust,grainsize, SPH_limits_file)

    ! Allocation dynamique
    ! We allocate the total number of SPH cells as the number of Voronoi cells mays vary
    if (lfirst_time) call alloc_dynamique(n_cells_max= n_SPH + n_etoiles)

    ! init random number generator
    stream = 0.0
    do i=1, nb_proc
       stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    enddo

    frac_E_stars=1.0 ! tous les photons partent de l'etoile
    call repartition_energie_etoiles()

    call init_reemission()
    call chauffage_interne()

    !$omp parallel default(none) private(lambda) shared(n_lambda)
    !$omp do schedule(static,1)
    do lambda=1, n_lambda
       call repartition_energie(lambda)
    enddo
    !$omp end do
    !$omp end parallel

    letape_th = .true.
    nbre_phot2 = nbre_photons_eq_th

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

    lambda=1
    laffichage=.true.
    nbre_phot2 = nbre_photons_lambda

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
       nnfot2 = 0.0_db
       photon : do while ((nnfot2 < nbre_phot2))
          nnfot2=nnfot2+1.0_db

          ! Choix longueur d'onde
          rand = sprng(stream(id))
          call select_wl_em(rand,lambda)

          ! Emission du paquet
          call emit_packet(id,lambda, icell,x,y,z,u,v,w,stokes,flag_star,flag_ISM,lintersect)
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

    Tdust = -1.0 ;

    do icell=1, n_cells
       i_SPH = Voronoi(icell)%id
       if (i_SPH > 0) Tdust(i_SPH) = Temperature(icell)
    enddo


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

    return

  end subroutine run_mcfost_phantom

!*************************************************************************

end module mcfost2phantom
