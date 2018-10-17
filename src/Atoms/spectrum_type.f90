MODULE spectrum_type

  use atom_type
  use grid_type, only : GridType

  IMPLICIT NONE


   integer, parameter :: N_MAX_OVERLAP=100
   real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
   real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352

  ! One activeset (active transitions of each atom) at each wavelength
  TYPE ActiveSetType
   integer, allocatable, dimension(:,:) :: lower_levels, upper_levels
   double precision, allocatable, dimension(:) :: chi, eta, chip, eta_c, chi_c, sca_c, &
      eta_Q, eta_U
               ! in case of Zeeman polarisation, chi, eta are 4*Nspace long instead of Nspace
   type (AtomicLine), allocatable, dimension(:,:) :: art
   type (AtomicContinuum), allocatable, dimension(:,:) :: crt
   !type (MolecularLine) :: mrt
  END TYPE ActiveSetType

  TYPE Spectrum
   type  (GridType), pointer :: atmos
   logical :: vacuum_to_air, updateJ
   integer :: Nwaves, Nact
   double precision, allocatable, dimension(:) :: lambda
   real, allocatable, dimension(:) :: J, I, StokesQ, StokesU, StokesV, J20, Flux
   type (ActiveSetType) :: ActiveSet !one active set for all wavelength at a cell point
   character:: Jfile, J20file
  END TYPE Spectrum

  type (Spectrum) :: NLTEspec

CONTAINS


  SUBROUTINE freeSpectrum()
   deallocate(NLTEspec%lambda)
   NULLIFY(NLTEspec%atmos)
   deallocate(NLTEspec%J, NLTEspec%I, NLTEspec%Flux)

  END SUBROUTINE freeSpectrum

  SUBROUTINE initAS(re_init)
    logical, intent(in) :: re_init
    integer :: Nsize
    
    Nsize = NLTEspec%Nwaves
    !if mag field, Nsize=4*Nsize
     ! Nsize2 = 3*Nsize
   ! first time allocate space on the wavelength grid
   if (.not.re_init) then
    if (NLTEspec%Nact.gt.0) then
     allocate(NLTEspec%Activeset%chi(Nsize))
    ! if mag field allocate(NLTEspec%Activesets%chip(NLTEspec%Nwaves))
     allocate(NLTEspec%Activeset%eta(Nsize))
     NLTEspec%Activeset%chi = 0.
     NLTEspec%Activeset%eta = 0.
    end if 
    allocate(NLTEspec%Activeset%eta_c(Nsize))
    allocate(NLTEspec%Activeset%chi_c(Nsize))
    allocate(NLTEspec%Activeset%sca_c(Nsize))
    allocate(NLTEspec%J(NLTEspec%Nwaves))
   ! if line pol
    !allocate(NLTEspec%Activesets%eta_Q(NLTEspec%Nwaves))
     !allocate(NLTEspec%Activesets%eta_U(NLTEspec%Nwaves))
     ! ...
    NLTEspec%Activeset%chi_c = 0.
    NLTEspec%Activeset%eta_c = 0.
    NLTEspec%Activeset%sca_c = 0.
    NLTEspec%J = 0.0
    !!as% etc etc = 0
   else
    !reset for next cell points
    if (NLTEspec%Nact.gt.0) then
     NLTEspec%Activeset%chi = 0.
     NLTEspec%Activeset%eta = 0.
    end if
    NLTEspec%Activeset%chi_c = 0.
    NLTEspec%Activeset%eta_c = 0.
    NLTEspec%Activeset%sca_c = 0.
    NLTEspec%J = 0.0
   end if
   
  RETURN
  END SUBROUTINE initAS

  SUBROUTINE freeAS()
    type (ActiveSetType) :: as

    !! etc etc ...
    if (NLTEspec%Nact.gt.0) then
     deallocate(NLTEspec%Activeset%chi)
    ! if mag field allocate(NLTEspec%Activesets%chip(NLTEspec%Nwaves))
     deallocate(NLTEspec%Activeset%eta)
    end if 
    deallocate(NLTEspec%Activeset%eta_c)
    deallocate(NLTEspec%Activeset%chi_c)
    deallocate(NLTEspec%Activeset%sca_c)
    
   RETURN
  END SUBROUTINE freeAS

  FUNCTION vaccum2air(Nlambda, lambda_vac) result(lambda_air)
   !wavelength in nm
   integer, intent(in) :: Nlambda
   real, dimension(:), intent(in) :: lambda_vac
   real, dimension(Nlambda) :: lambda_air
   integer :: la
   real :: sqwave, reduction

   do la=0,Nlambda
    if (lambda_vac(la).ge.VACUUM_TO_AIR_LIMIT) then
     sqwave = 1./(lambda_vac(la)**2)
     reduction = 1. + 2.735182e-4 + &
            (1.314182 + 2.76249e+4 * sqwave) * sqwave
     lambda_air(la) = lambda_vac(la) / reduction
    else
     lambda_air(la) = lambda_vac(la)
    end if
   end do

  RETURN
  END FUNCTION vaccum2air

  FUNCTION air2vaccum(Nlambda, lambda_air) result(lambda_vac)
  !wavelength in nm
  integer, intent(in) :: Nlambda
  real, dimension(:), intent(in) :: lambda_air
  real, dimension(Nlambda) :: lambda_vac
  integer :: la
  real :: sqwave, increase

  do la=0,Nlambda
   if (lambda_air(la).ge.AIR_TO_VACUUM_LIMIT) then
    sqwave = (1.0e+7 / lambda_air(la))**2
    increase = 1.0000834213E+00 + &
            2.406030E+06/(1.30E+10 - sqwave) + &
            1.5997E+04/(3.89E+09 - sqwave)
    lambda_vac(la) = lambda_air(la) * increase
   else
    lambda_vac(la) = lambda_air(la)
   end if
  end do

  RETURN
  END FUNCTION air2vaccum


END MODULE spectrum_type

