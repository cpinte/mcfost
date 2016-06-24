!--  subroutine init_mcfost_phantom(mcfost_filename, ierr)
!--
!--    ! This routine should be in mcfost2phantom
!--    character(len=*), intent(in) :: mcfost_filename
!--    integer, intent(out) :: ierr
!--
!--    ! options
!--    write(*,*) "INIT MCFOST"
!--
!--    ! parameter file
!--
!--    ! dust properties
!--
!--    ierr = 0
!--    return
!--
!--  end subroutine init_mcfost_phantom
!--
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
