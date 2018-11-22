MODULE spectrum_type

  use atom_type
  use atmos_type, only : GridType
  use getlambda, only : make_wavelength_grid
  
  !MCFOST's original modules
  use fits_utils, only : print_error

  IMPLICIT NONE


   real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
   real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352
   character(len=*), parameter :: WAVES_FILE="wavelength_grid.fits.gz"
   ! not that FLUX_FILE is 1D only if one pixel, otherwise it is a map.
   ! F(x,y,0,lambda) in several directions.
   character(len=*), parameter :: FLUX_FILE="flux.fits.gz"
   ! I(x,y,0,lambda,imu), not used yet
   !character(len=*), parameter :: SPEC_FILE="spectrum.fits.gz"

  ! at all wavelength at a time
  TYPE ActiveSetType
   integer, allocatable, dimension(:,:) :: lower_levels, upper_levels
   double precision, allocatable, dimension(:,:) :: chi, eta, chip, eta_c, chi_c, sca_c, &
      eta_Q, eta_U, eta_c_bf, chi_c_bf ! *_bf is used for pure continuum opacities.
               ! Remember, chi_c, eta_c contain all opacities that are LTE (including lines of PASSIVE atoms)
               ! in case of Zeeman polarisation, chi, eta are 4*Nspace long instead of Nspace
   type (AtomicLine), allocatable, dimension(:,:) :: art
   type (AtomicContinuum), allocatable, dimension(:,:) :: crt
   double precision, allocatable, dimension(:) :: Psi
  END TYPE ActiveSetType

  TYPE Spectrum
   type  (GridType), pointer :: atmos
   logical :: vacuum_to_air=.false., updateJ, write_wavelength_grid=.false.
   integer :: Nwaves, Nact, NPROC
   double precision :: wavelength_ref=0.d0 !nm optionnal
   double precision, allocatable, dimension(:) :: lambda
   !nproc, nlambda, ncells
   double precision, allocatable, dimension(:,:,:) :: I, StokesQ, StokesU, StokesV, Ic
   !nproc, nlambda
   double precision, allocatable, dimension(:,:) :: J, J20, Jc
   double precision, allocatable, dimension(:,:,:,:,:) :: Flux, Fluxc
   ! Flux is a map of Nlambda, xpix, ypix, nincl, nazimuth
   type (ActiveSetType) :: ActiveSet !one active set for all wavelength at a cell point
   character:: Jfile, J20file
  END TYPE Spectrum

  type (Spectrum) :: NLTEspec

CONTAINS

  SUBROUTINE initSpectrum(NPROC, lam0, vacuum_to_air, write_wavelength)
   integer, intent(in) :: NPROC
   double precision, optional :: lam0
   logical, optional :: vacuum_to_air, write_wavelength
   
   NLTEspec%Nact = NLTEspec%atmos%Nactiveatoms
   NLTEspec%NPROC = NPROC
   
   if (present(lam0)) NLTEspec%wavelength_ref = lam0
   if (present(write_wavelength)) NLTEspec%write_wavelength_grid = write_wavelength
   if (present(vacuum_to_air)) NLTEspec%vacuum_to_air = vacuum_to_air
  
   ! Initialize the wavelength grid depending on the number of transitions PASSIVE/ACTIVE
   CALL make_wavelength_grid(NLTEspec%atmos%Natom, NLTEspec%atmos%Atoms, & 
                        NLTEspec%lambda, NLTEspec%wavelength_ref)
   NLTEspec%Nwaves = size(NLTEspec%lambda)
   CALL writeWavelength()  
  
  RETURN
  END SUBROUTINE initSpectrum

  SUBROUTINE allocSpectrum(NPIX_X, NPIX_Y, N_INCL, N_AZIMUTH)
   integer, intent(in) :: NPIX_X, NPIX_Y, N_INCL, N_AZIMUTH
   
   if (allocated(NLTEspec%I)) then
    write(*,*) "Error I already allocated"
    stop
   end if
   
   allocate(NLTEspec%I(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS))
   allocate(NLTEspec%Ic(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS))

   if (NLTEspec%atmos%Magnetized) then
    allocate(NLTEspec%StokesQ(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS), & 
             NLTEspec%StokesU(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS), &
             NLTEspec%StokesV(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS))
   end if
   allocate(NLTEspec%J(NLTEspec%NPROC, NLTEspec%Nwaves))
   allocate(NLTEspec%Jc(NLTEspec%NPROC, NLTEspec%Nwaves))
   !! allocate(NLTEspec%J20(NLTEspec%NPROC, NLTEspec%Nwaves))
   
   allocate(NLTEspec%Flux(NLTEspec%Nwaves,NPIX_X, NPIX_Y,N_INCL,N_AZIMUTH))
   allocate(NLTEspec%Fluxc(NLTEspec%Nwaves,NPIX_X, NPIX_Y,N_INCL,N_AZIMUTH))
   
   ! --> initialised in integ_ray_line  
   !   NLTEspec%I = 0.0
   !   NLTEspec%Ic = 0.0
   NLTEspec%Flux = 0.0
   NLTEspec%Fluxc = 0.0
   ! J and Jc are initialised in initAS()
  
  RETURN
  END SUBROUTINE allocSpectrum

  SUBROUTINE freeSpectrum()
   deallocate(NLTEspec%lambda)
   deallocate(NLTEspec%J, NLTEspec%I, NLTEspec%Flux)
   deallocate(NLTEspec%Jc, NLTEspec%Ic, NLTEspec%Fluxc)
   if (NLTEspec%atmos%Magnetized) deallocate(NLTEspec%StokesQ, NLTEspec%StokesU, NLTEspec%StokesV)
   !! deallocate(NLTEspec%J20)
   NULLIFY(NLTEspec%atmos)
  RETURN
  END SUBROUTINE freeSpectrum

  SUBROUTINE initAS(id, re_init)
    ! allocate arrays for re_init == FALSE
    ! if re_init is .true. only set opac arrays to 0 for next cell point.
    ! when re_init is FALSE, opacities for all threads are initialised
    ! so that initAS becomes independent on id
    logical, intent(in) :: re_init
    integer, intent(in) :: id
    integer :: Nproc, Nsize !, Nsize2
    
    Nsize = NLTEspec%Nwaves
    Nproc = NLTEspec%NPROC
    !if mag field, Nsize=4*Nsize
     ! Nsize2 = 3*Nsize
   ! first time allocate space on the wavelength grid
   if (.not.re_init) then !set to 0d0 for all threads
    !if (NLTEspec%Nact.gt.0) then
     allocate(NLTEspec%Activeset%chi(Nproc,Nsize))
    ! if mag field allocate(NLTEspec%Activesets%chip(Nproc,Nsize))
     allocate(NLTEspec%Activeset%eta(Nproc,Nsize))
     NLTEspec%Activeset%chi = 0.
     NLTEspec%Activeset%eta = 0.
    !end if 
    allocate(NLTEspec%Activeset%eta_c(Nproc,Nsize))
    allocate(NLTEspec%Activeset%chi_c(Nproc,Nsize))
    allocate(NLTEspec%Activeset%sca_c(Nproc,Nsize))
    allocate(NLTEspec%Activeset%eta_c_bf(Nproc,Nsize))
    allocate(NLTEspec%Activeset%chi_c_bf(Nproc,Nsize))

   ! if line pol
    !allocate(NLTEspec%Activesets%eta_Q(Nproc,Nsize))
     !allocate(NLTEspec%Activesets%eta_U(Nproc,Nsize))
     ! ...
    NLTEspec%Activeset%chi_c = 0.
    NLTEspec%Activeset%eta_c = 0.
    NLTEspec%Activeset%sca_c = 0.
    NLTEspec%Activeset%chi_c_bf = 0.
    NLTEspec%Activeset%eta_c_bf = 0.
    NLTEspec%J = 0.0
    NLTEspec%Jc = 0.0

   else
    !reset for next cell points but only for a specific thread id
    !if (NLTEspec%Nact.gt.0) then
     NLTEspec%Activeset%chi(id,:) = 0.
     NLTEspec%Activeset%eta(id,:) = 0.
    !end if
    NLTEspec%Activeset%chi_c(id,:) = 0.
    NLTEspec%Activeset%eta_c(id,:) = 0.
    NLTEspec%Activeset%sca_c(id,:) = 0.
    NLTEspec%Activeset%chi_c_bf(id,:) = 0.
    NLTEspec%Activeset%eta_c_bf(id,:) = 0.
    NLTEspec%J(id,:) = 0.0
    NLTEspec%Jc(id,:) = 0.0
   end if
   
  RETURN
  END SUBROUTINE initAS

  SUBROUTINE freeAS()
    type (ActiveSetType) :: as

    !! etc etc ...
    if (NLTEspec%Nact.gt.0) then
     deallocate(NLTEspec%Activeset%chi)
     deallocate(NLTEspec%Activeset%eta)
    end if 
    deallocate(NLTEspec%Activeset%eta_c)
    deallocate(NLTEspec%Activeset%chi_c)
    deallocate(NLTEspec%Activeset%sca_c)
    deallocate(NLTEspec%Activeset%chi_c_bf)
    deallocate(NLTEspec%Activeset%eta_c_bf)
    
   RETURN
  END SUBROUTINE freeAS
  
  SUBROUTINE writeWavelength()
   ! Write wavelength grid build with each transition of each atom
   integer :: unit, EOF = 0, blocksize, naxes(1), naxis,group, bitpix, fpixel
   logical :: extend, simple
   character(len=6) :: comment="VACUUM"
   double precision :: lambda_vac(NLTEspec%Nwaves)
   
   if (.not.allocated(NLTEspec%lambda).or.&
       .not.NLTEspec%write_wavelength_grid) RETURN !
   CALL ftgiou(unit,EOF)
   blocksize=1
   CALL ftinit(unit,trim(WAVES_FILE),blocksize,EOF)
   simple = .true.
   group = 1
   fpixel = 1
   extend = .false. 
   bitpix = -64   
   naxis = 1
   naxes(1) = NLTEspec%Nwaves
   
   if (NLTEspec%vacuum_to_air) then
     comment="AIR"
     lambda_vac = NLTEspec%lambda
     NLTEspec%lambda = vacuum2air(NLTEspec%Nwaves, lambda_vac)
   end if 
   
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
   CALL ftpkys(unit, "UNIT", "nm", comment, EOF)
   CALL ftpprd(unit,group,fpixel,naxes(1),NLTEspec%lambda,EOF)
   CALL ftclos(unit, EOF)
   CALL ftfiou(unit, EOF)
   
   if (EOF > 0) CALL print_error(EOF)       
 
  RETURN
  END SUBROUTINE writeWavelength

  FUNCTION vacuum2air(Nlambda, lambda_vac) result(lambda_air)
   !wavelength in nm
   integer, intent(in) :: Nlambda
   double precision, dimension(:), intent(in) :: lambda_vac
   double precision, dimension(Nlambda) :: lambda_air
   double precision, dimension(Nlambda) :: sqwave, reduction
!   integer :: la
!   double precision :: sqwave, reduction

!    do la=0,Nlambda
!     if (lambda_vac(la).ge.VACUUM_TO_AIR_LIMIT) then
!      sqwave = 1./(lambda_vac(la)**2)
!      reduction = 1. + 2.735182e-4 + &
!             (1.314182 + 2.76249e+4 * sqwave) * sqwave
!      lambda_air(la) = lambda_vac(la) / reduction
!     else
!      lambda_air(la) = lambda_vac(la)
!     end if
!    end do

   where (lambda_vac.ge.VACUUM_TO_AIR_LIMIT) 
     sqwave = 1./(lambda_vac**2)
     reduction = 1. + 2.735182e-4 + &
            (1.314182 + 2.76249e+4 * sqwave) * sqwave
     lambda_air = lambda_vac / reduction
   else where(lambda_vac.lt.VACUUM_TO_AIR_LIMIT)
     lambda_air = lambda_vac
   end where


  RETURN
  END FUNCTION vacuum2air

  FUNCTION air2vacuum(Nlambda, lambda_air) result(lambda_vac)
  !wavelength in nm
  integer, intent(in) :: Nlambda
  double precision, dimension(:), intent(in) :: lambda_air
  double precision, dimension(Nlambda) :: lambda_vac
  double precision, dimension(Nlambda) :: sqwave, increase
!   integer :: la
!   double precision :: sqwave, increase

!   do la=0,Nlambda
!    if (lambda_air(la).ge.AIR_TO_VACUUM_LIMIT) then
!     sqwave = (1.0e+7 / lambda_air(la))**2
!     increase = 1.0000834213E+00 + &
!             2.406030E+06/(1.30E+10 - sqwave) + &
!             1.5997E+04/(3.89E+09 - sqwave)
!     lambda_vac(la) = lambda_air(la) * increase
!    else
!     lambda_vac(la) = lambda_air(la)
!    end if
!   end do

   where (lambda_air.ge.AIR_TO_VACUUM_LIMIT)
    sqwave = (1.0e+7 / lambda_air)**2
    increase = 1.0000834213E+00 + &
            2.406030E+06/(1.30E+10 - sqwave) + &
            1.5997E+04/(3.89E+09 - sqwave)
    lambda_vac = lambda_air * increase
   else where(lambda_air.lt.AIR_TO_VACUUM_LIMIT)
    lambda_vac = lambda_air
   end where

  RETURN
  END FUNCTION air2vacuum


END MODULE spectrum_type

