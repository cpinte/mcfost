module phantom2mcfost

  use parametres
  use constantes

  implicit none

contains


  subroutine init_mcfost_phantom(mcfost_filename, ierr)


    character(len=*), intent(in) :: mcfost_filename
    integer, intent(out) :: ierr

    ! options
    write(*,*) "INIT MCFOST"

    ! parameter file

    ! dust properties

    ierr = 0
    return

  end subroutine init_mcfost_phantom

  !*************************************************************************

  subroutine run_mcfost_phantom(np,xyzh, iphase, dustfrac, ntypes, massoftype, hfact, npoftype, Tdust, mu_gas, ierr)

    use io_phantom

    integer, intent(in) :: np, ntypes
    real(db), dimension(4,np), intent(in) :: xyzh
    integer, dimension(np), intent(in) :: iphase
    real(sl), dimension(np), intent(in) :: dustfrac
    real(db), dimension(ntypes), intent(in) ::massoftype
    real(db), intent(in) :: hfact
    integer, dimension(ntypes), intent(in) :: npoftype

    real(db), dimension(np), intent(out) :: Tdust
    integer, intent(out) :: ierr
    real(db), intent(out) :: mu_gas

    real(db), dimension(:), allocatable :: x,y,z,rhogas,rhodust
    integer :: ncells

    mu_gas = mu

    ierr = 0

    call phantom_2_mcfost(np,xyzh,iphase,dustfrac,ntypes,massoftype(1:ntypes),hfact, x,y,z,rhogas,rhodust,ncells)
    if (ncells <= 0) then
       ierr = 1
       return
    endif

    call setup_mcfost_Voronoi_grid()




  end subroutine run_mcfost_phantom

!*************************************************************************

  subroutine setup_mcfost_Voronoi_grid()



    write(*,*) "MCFOST setup Voronoi"

  end subroutine setup_mcfost_Voronoi_grid

end module phantom2mcfost
