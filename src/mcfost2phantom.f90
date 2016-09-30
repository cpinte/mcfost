module mcfost2phantom

  use parametres
  use constantes
  use utils
  use Voronoi_grid
  use em_th, only : temperature

  implicit none

contains

  subroutine init_mcfost_phantom(mcfost_para_filename, ierr) !,  np, nptmass, ntypes, ndusttypes, npoftype)

    use init_mcfost, only : set_default_variables, get_mcfost_utils_dir
    use read_params, only : read_para
    use disk, only : n_zones
    use dust_transfer, only : transfert_poussiere


    ! This routine should be in mcfost2phantom
    character(len=*), intent(in) :: mcfost_para_filename
    integer, intent(out) :: ierr

    ! Global logical variables
    call set_default_variables()

    ! Looking for the mcfost utils directory
    call get_mcfost_utils_dir()

    ! parameter file
    call read_para(mcfost_para_filename)

    ! Setting option for the mcfost2phantom interface
    ltemp = .true. ; lsed = .false. ; lsed_complete = .false.
    lVoronoi = .true. ; l3D = .true.

    if (n_zones > 1) then
       write(*,*) "ERROR: mcfost2phantom only works with a 1zone parameter file"
       write(*,*) "Exiting"
       ierr = 1
       return
    endif

    ! dust properties

    ! making the MC run
    call transfert_poussiere()

    ierr = 0
    return

  end subroutine init_mcfost_phantom

  !*************************************************************************

!--  subroutine run_mcfost_phantom(np,nptmass,ntypes,ndusttypes,npoftype,xyzh,iphase,grainsize,&
!--       dustfrac, massoftype,xyzmh_ptmass,hfact,umass,utime,udist,graindens,compute_Frad, Tdust, Frad,mu_gas,ierr)
!--
!--
!--    ! This routine should be in mcfost2phantom
!--    use io_phantom
!--
!--    integer, intent(in) :: np, nptmass, ntypes,ndusttypes
!--    real(db), dimension(4,np), intent(in) :: xyzh
!--    integer(kind=1), dimension(np), intent(in) :: iphase
!--    real(db), dimension(ndusttypes,np), intent(in) :: dustfrac
!--    real(db), dimension(ndusttypes), intent(in) :: grainsize
!--    real(db), dimension(ntypes), intent(in) :: massoftype
!--    real(db), intent(in) :: hfact, umass, utime, udist, graindens
!--    real(db), dimension(:,:), intent(in) :: xyzmh_ptmass
!--    integer, dimension(ntypes), intent(in) :: npoftype
!--    logical, intent(in) :: compute_Frad ! does mcfost need to compute the radiation pressure
!--
!--    real, dimension(np), intent(out) :: Tdust ! mcfost stores Tdust as real, not db
!--    real, dimension(3,ndusttypes,np), intent(out) :: Frad
!--    real(db), intent(out) :: mu_gas
!--    integer, intent(out) :: ierr
!--
!--
!--    real(db), dimension(:), allocatable :: x,y,z,rhogas
!--    real(db), dimension(:,:), allocatable :: rhodust
!--    integer :: ncells
!--
!--    ierr = 0
!--    mu_gas = mu ! Molecular weight
!--    Frad = 0.
!--
!--    call phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,xyzh,iphase,grainsize,dustfrac,&
!--         massoftype(1:ntypes),xyzmh_ptmass,hfact,umass,utime,udist,graindens,x,y,z,rhogas,rhodust,ncells)
!--
!--    if (ncells <= 0) then
!--       ierr = 1
!--       return
!--    endif
!--
!--    !call setup_mcfost_Voronoi_grid()
!--
!--    Tdust = 2.73
!--
!--    return
!--
!--  end subroutine run_mcfost_phantom

!*************************************************************************

end module mcfost2phantom
