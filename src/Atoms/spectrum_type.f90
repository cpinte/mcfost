MODULE spectrum_type

  use atom_type
  use atmos_type, only : GridType
  use getlambda, only  : make_wavelength_grid, Nred_array, Nmid_array, Nblue_array, &
  					    adjust_wavelength_grid
  
  !MCFOST's original modules
  use fits_utils, only : print_error
  use parametres
  use input

  IMPLICIT NONE


   real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
   real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352
   character(len=*), parameter :: WAVES_FILE="wavelength_grid.fits.gz"
   ! not that FLUX_FILE is 1D only if one pixel, otherwise it is a map.
   ! F(x,y,0,lambda) in several directions.
   character(len=*), parameter :: FLUX_FILE="flux.fits.gz"
   ! I(x,y,0,lambda,imu), not used yet
   !character(len=*), parameter :: SPEC_FILE="spectrum.fits.gz"
   character(len=*), parameter :: OPAC_CONTRIB_FILE="opacities.fits.gz"

  !Like RH, should vanish in the futur
!   TYPE ActiveSetType
!    integer, allocatable, dimension(:,:) :: lower_levels, upper_levels
!    type (AtomicLine), allocatable, dimension(:,:) :: art
!    type (AtomicContinuum), allocatable, dimension(:,:) :: crt
!    double precision, allocatable, dimension(:,:) :: Psi !proc and wavelength
!   END TYPE ActiveSetType
  
  TYPE AtomicOpacity
   !active opacities
   double precision, allocatable, dimension(:,:)   :: chi, eta, rho
   !passive opacities
   double precision, allocatable, dimension(:,:)   :: rho_p, eta_p, chi_p
   double precision, allocatable, dimension(:,:)   :: eta_c, chi_c, sca_c
   double precision, allocatable, dimension(:,:)   :: jc
   double precision, allocatable, dimension(:,:,:) :: Kc
   !type (ActiveSetType) :: ActiveSet
  END TYPE AtomicOpacity

  TYPE Spectrum
   type  (GridType), pointer :: atmos
   logical :: vacuum_to_air=.false., updateJ, write_wavelength_grid=.false.
   integer :: Nwaves, Nact, Npass, Ntrans, NPROC=1
   double precision :: wavelength_ref=0.d0 !nm optionnal
   double precision, allocatable, dimension(:) :: lambda
   !nproc, nlambda, nrays
   double precision, allocatable, dimension(:,:,:) :: I, StokesQ, StokesU, StokesV, Ic
   !nproc, nlambda
   double precision, allocatable, dimension(:,:) :: J, J20, Jc
   !Nlambda, xpix, ypix, Nincl, Naz
   double precision, allocatable, dimension(:,:,:,:,:) :: Flux, Fluxc
   !Contribution function
   double precision, allocatable, dimension(:,:,:) :: Ksi 
   ! Flux is a map of Nlambda, xpix, ypix, nincl, nazimuth
   type (AtomicOpacity) :: AtomOpac
   character:: Jfile, J20file
  END TYPE Spectrum

  type (Spectrum) :: NLTEspec

  CONTAINS

  SUBROUTINE initSpectrum(lam0, vacuum_to_air, write_wavelength)
   ! ------------------------------------------- !
    ! Allocate and store the wavelength grid for
    ! NLTE atomic line transfer.
   ! ------------------------------------------- !

   !integer, intent(in) :: NPROC
   double precision, optional :: lam0
   logical, optional :: vacuum_to_air, write_wavelength
      
   NLTEspec%Nact = NLTEspec%atmos%Nactiveatoms
   NLTEspec%Npass = NLTEspec%atmos%Npassiveatoms

   NLTEspec%NPROC = nb_proc!NPROC
   
   if (present(lam0)) NLTEspec%wavelength_ref = lam0
   if (present(write_wavelength)) NLTEspec%write_wavelength_grid = write_wavelength
   if (present(vacuum_to_air)) NLTEspec%vacuum_to_air = vacuum_to_air
  
   ! Initialize the wavelength grid depending on the number of transitions PASSIVE/ACTIVE
   ! and allocate NLTEspec%lambda
   ! Set also Nblue and Nred for each transitions: extension of the transtion
   ! in absolute indexes on the whole grid.
   CALL make_wavelength_grid(NLTEspec%atmos%Natom, NLTEspec%atmos%Atoms, & 
                        NLTEspec%lambda, NLTEspec%Ntrans, NLTEspec%wavelength_ref)
   NLTEspec%Nwaves = size(NLTEspec%lambda)
   CALL writeWavelength()
  
  RETURN
  END SUBROUTINE initSpectrum
  
  SUBROUTINE initSpectrumImage(lam0, dL, resol)
  ! -------------------------------------------------------------------- !
   ! Allocate a special wavelength grid for emission_line map.
   !This allows to solve for the RT equation only for a specific line.
   ! a LINEAR wavelength grid spanning from lam0-dL to lam0+dL with a
   ! constant resolution in km/s given by resol is built.
  ! -------------------------------------------------------------------- !
   double precision, intent(in) :: lam0, resol, dL
      
   NLTEspec%Nact = NLTEspec%atmos%Nactiveatoms
   NLTEspec%Npass = NLTEspec%atmos%Npassiveatoms
   NLTEspec%NPROC = nb_proc
   
   !reallocate wavelength arrays and indexes
   CALL adjust_wavelength_grid()
   
   NLTEspec%Nwaves = size(NLTEspec%lambda)
  
  RETURN
  END SUBROUTINE initSpectrumImage

  SUBROUTINE allocSpectrum()!NPIX_X, NPIX_Y, N_INCL, N_AZIMUTH)
   !integer, intent(in) :: NPIX_X, NPIX_Y, N_INCL, N_AZIMUTH
   
   integer :: Nsize
   
   Nsize = NLTEspec%Nwaves
   ! if pol, Nsize = 3*Nsize for rho and 4*Nsize for eta, chi
   
   if (allocated(NLTEspec%I)) then
    write(*,*) "Error I already allocated"
    stop
   end if
   
   allocate(NLTEspec%I(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS))
   allocate(NLTEspec%Ic(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS))
   NLTEspec%I = 0d0
   NLTEspec%Ic = 0d0

   if (NLTEspec%atmos%Magnetized) then
    allocate(NLTEspec%StokesQ(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS), & 
             NLTEspec%StokesU(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS), &
             NLTEspec%StokesV(NLTEspec%NPROC, NLTEspec%NWAVES, NLTEspec%atmos%NRAYS))
    NLTEspec%StokesQ=0d0
    NLTEspec%StokesU=0d0
    NLTEspec%StokesV=0d0
   end if
   allocate(NLTEspec%J(NLTEspec%NPROC, NLTEspec%Nwaves))
   allocate(NLTEspec%Jc(NLTEspec%NPROC, NLTEspec%Nwaves))
   !! allocate(NLTEspec%J20(NLTEspec%NPROC, NLTEspec%Nwaves))
   NLTEspec%J = 0.0
   NLTEspec%Jc = 0.0
      
   allocate(NLTEspec%Flux(NLTEspec%Nwaves,NPIX_X, NPIX_Y,RT_N_INCL,RT_N_AZ))
   allocate(NLTEspec%Fluxc(NLTEspec%Nwaves,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ))
   !create pol flux map or add to Flux
   
   NLTEspec%Flux = 0.0
   NLTEspec%Fluxc = 0.0
   
   !Now opacities
   if (lstore_opac) then !keep continuum LTE opacities in memory
     !sca_c = Kc(:,:,2), chi_c = Kc(:,:,1), eta_c = jc
     allocate(NLTEspec%AtomOpac%Kc(NLTEspec%atmos%Nspace,NLTEspec%Nwaves,2), &
       NLTEspec%AtomOpac%jc(NLTEspec%atmos%Nspace,NLTEspec%Nwaves))
     NLTEspec%AtomOpac%Kc = 0d0
     NLTEspec%AtomOpac%jc = 0d0
   else
    allocate(NLTEspec%AtomOpac%eta_c(NLTEspec%NPROC,Nsize))
    allocate(NLTEspec%AtomOpac%chi_c(NLTEspec%NPROC,Nsize))
    allocate(NLTEspec%AtomOpac%sca_c(NLTEspec%NPROC,NLTEspec%Nwaves))
    NLTEspec%AtomOpac%chi_c = 0.
    NLTEspec%AtomOpac%eta_c = 0.
    NLTEspec%AtomOpac%sca_c = 0.
   end if
   
   allocate(NLTEspec%AtomOpac%chi(NLTEspec%NPROC,Nsize))
   allocate(NLTEspec%AtomOpac%eta(NLTEspec%NPROC,Nsize))
   ! if pol allocate AtomOpac%rho
   NLTEspec%AtomOpac%chi = 0.
   NLTEspec%AtomOpac%eta = 0.

   allocate(NLTEspec%AtomOpac%eta_p(NLTEspec%NPROC,Nsize))
   allocate(NLTEspec%AtomOpac%chi_p(NLTEspec%NPROC,Nsize))
   !allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%NPROC,Nsize))
   ! if pol same as atice opac case

   NLTEspec%AtomOpac%chi_p = 0.
   NLTEspec%AtomOpac%eta_p = 0.
   
   !Contribution function
   if (lcontrib_function) then
    allocate(NLTEspec%Ksi(NLTEspec%Ntrans, NLTEspec%atmos%Nspace, NLTEspec%Nwaves))
    NLTEspec%Ksi = 0d0
   end if

  RETURN
  END SUBROUTINE allocSpectrum

  SUBROUTINE freeSpectrum()   
   deallocate(NLTEspec%lambda)
   deallocate(NLTEspec%J, NLTEspec%I, NLTEspec%Flux)
   deallocate(NLTEspec%Jc, NLTEspec%Ic, NLTEspec%Fluxc)
   if (NLTEspec%atmos%Magnetized) deallocate(NLTEspec%StokesQ, NLTEspec%StokesU, NLTEspec%StokesV)
   !! deallocate(NLTEspec%J20)
   
   !active
   deallocate(NLTEspec%AtomOpac%chi)
   deallocate(NLTEspec%AtomOpac%eta)

   !passive
   deallocate(NLTEspec%AtomOpac%eta_p) !contains only passive lines if store_opac
   deallocate(NLTEspec%AtomOpac%chi_p)
   !deallocate(NLTEspec%AtomOpac%rho_p)
   if (lstore_opac) then !keep continuum LTE opacities in memory
     deallocate(NLTEspec%AtomOpac%Kc,  NLTEspec%AtomOpac%jc)
   else !they are not allocated if we store background continua on ram
    deallocate(NLTEspec%AtomOpac%eta_c)
    deallocate(NLTEspec%AtomOpac%chi_c)
    deallocate(NLTEspec%AtomOpac%sca_c)
   end if
   NULLIFY(NLTEspec%atmos)
   
   deallocate(Nblue_array, Nmid_array, Nred_array)
  RETURN
  END SUBROUTINE freeSpectrum

  SUBROUTINE initAtomOpac(id)
    ! set opacities to 0d0 for thread id
    integer, intent(in) :: id
    
    if (id <= 0) then
     write(*,*) "(initAtomOpac) thread id has to be >= 1!"
     stop
    end if

    NLTEspec%AtomOpac%chi(id,:) = 0d0
    NLTEspec%AtomOpac%eta(id,:) = 0d0
    if (.not.lstore_opac) then
      NLTEspec%AtomOpac%chi_c(id,:) = 0d0
      NLTEspec%AtomOpac%eta_c(id,:) = 0d0
      NLTEspec%AtomOpac%sca_c(id,:) = 0d0
    end if !else thay are not allocated
    NLTEspec%AtomOpac%chi_p(id,:) = 0d0
    NLTEspec%AtomOpac%eta_p(id,:) = 0d0
    NLTEspec%J(id,:) = 0d0
    NLTEspec%Jc(id,:) = 0d0
  
  RETURN
  END SUBROUTINE initAtomOpac
  

 SUBROUTINE WRITE_FLUX()
 ! -------------------------------------------------- !
  ! Write the spectral Flux map on disk.
  ! FLUX map:
  ! NLTEspec%Flux total and NLTEspec%Flux continuum
 ! --------------------------------------------------- !
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements, i, xcenter
  integer :: la, Nred, Nblue, kr, kc, m, Nmid
  logical :: simple, extend
  character(len=6) :: comment="VACUUM"
  double precision :: lambda_vac(NLTEspec%Nwaves)
  real :: pixel_scale_x, pixel_scale_y 
  
  write(*,*) "Writing Flux-map"
  write(*,*) "npix_x = ", npix_x, " npix_y = ", npix_y, ' RT method:', RT_line_method
  write(*,*) "Wavelength points:", NLTEspec%Nwaves
  
  if (l_sym_ima) then 
   write(*,*) " image is symmetric."
  else
   write(*,*) " image is not symmetric."
  end if
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
   status=0
   CALL ftgiou (unit,status)

   !  Create the new empty FITS file.
   blocksize=1
   CALL ftinit(unit,trim(FLUX_FILE),blocksize,status)

   simple=.true.
   extend=.true.
   group=1
   fpixel=1

   bitpix=-64
   naxis=5
   naxes(1)=NLTEspec%Nwaves!1!1 if only one wavelength

   if (RT_line_method==1) then
     naxes(2)=1
     naxes(3)=1
   else
     naxes(2)=npix_x
     naxes(3)=npix_y
   endif
   naxes(4)=RT_n_incl
   naxes(5)=RT_n_az
   nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
  ! write(*,*) (naxes(i), i=1,naxis)

  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

   !!RAC, DEC, reference pixel & pixel scale en degres
  CALL ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  CALL ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
  CALL ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
  pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  CALL ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)
 
  CALL ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  CALL ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
  CALL ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
  pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
  CALL ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)

  CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'F_nu',status)

  if (l_sym_ima) then 
   write(*,*) "Check symmetries with RT methods"
   if (RT_line_method==1) then ! what should I do in my case ?
    write(*,*) 'RT1 symmetry?'
    
   else if (RT_line_method == 2) then !/= 1, ??
    xcenter = npix_x/2 + modulo(npix_x,2) 
    if (.not.lstatic) then !due to velocity the lines could be asym even if the model is sym
    					   ! for instance keplerian profile is not the same for the two parts
    					   ! of the final image
     if (lkeplerian) then !line profile reversed
      do i=xcenter+1,npix_x
!        do m=1,NLTEspec%atmos%Natom
!         !!lines are sensitive to symmetry
!         do kr=1,NLTEspec%atmos%Atoms(m)%Nline
!          Nred = NLTEspec%atmos%Atoms(m)%lines(kr)%Nred
!          Nblue = NLTEspec%atmos%Atoms(m)%lines(kr)%Nblue
!          Nmid = NLTEspec%atmos%Atoms(m)%lines(kr)%Nmid
!          !Keplerian rotation inverses profile slope
!          NLTEspec%Flux(Nblue:Nmid-1,i,:,:,:) = NLTEspec%Flux(Nmid+1:Nred,npix_x-i+1,:,:,:)
!         end do 
!         !!continua are not sensitive to symmetry right ?
!         do kc=1,NLTEspec%atmos%Atoms(m)%Ncont
!          Nred = NLTEspec%atmos%Atoms(m)%lines(kr)%Nred
!          Nblue = NLTEspec%atmos%Atoms(m)%lines(kr)%Nblue
!          NLTEspec%Flux(:,i,:,:,:) = NLTEspec%Flux(:,npix_x-i+1,:,:,:)
!         end do
!        end do !atom  
       write(*,*) "RT2 symmetry ? "
       do m=1,NLTEspec%Ntrans
        Nred = Nred_array(m)
        Nblue = Nblue_array(m)
        Nmid = Nmid_array(m)
        NLTEspec%Flux(Nblue:Nmid-1,i,:,:,:) = NLTEspec%Flux(Nmid+1:Nred,npix_x-i+1,:,:,:)
       end do  
      end do !pix
     else ! not keplerian (infall or expansion)
      do i=xcenter+1,npix_x
       NLTEspec%Flux(:,i,:,:,:) = NLTEspec%Flux(:,npix_x-i+1,:,:,:)
      end do
     end if !lkeplerian  
    else !static, not moving but image still has symmetries
     do i=xcenter+1,npix_x !case of symmetric lines or continua
    					  !copy the left part of the image to the right
      NLTEspec%Flux(:,i,:,:,:) = NLTEspec%Flux(:,npix_x-i+1,:,:,:)
     end do 
    end if !if not static or if moving
   else !RT line method /= 0, 2
    write(*,*) "RT_line_method", RT_line_method, " not valid"
    stop
   endif !RT_line_method
  end if ! l_sym_image

  !  Write the array to the FITS file.
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%Flux,status)
  
  ! create new hdu for continuum
  CALL ftcrhd(unit, status)

  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  CALL ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
  CALL ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
  CALL ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
  pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
  CALL ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)
 
  CALL ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
  CALL ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
  CALL ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
  pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
  CALL ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)
  CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'F_nu',status)
  
  !the continuum should not be affected by velocity shift
  if (l_sym_ima.and.(RT_line_method == 2)) then
   xcenter = npix_x/2 + modulo(npix_x,2)
   do i=xcenter+1,npix_x
    NLTEspec%Fluxc(:,i,:,:,:) = NLTEspec%Fluxc(:,npix_x-i+1,:,:,:)
   end do
  end if ! l_sym_image  
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%Fluxc,status)
  
  
  ! create new hdu for wavelength grid
  CALL ftcrhd(unit, status)
   if (NLTEspec%vacuum_to_air) then
     comment="AIR"
     lambda_vac = NLTEspec%lambda
     NLTEspec%lambda = vacuum2air(NLTEspec%Nwaves, lambda_vac)
   end if 
   
   naxis = 1
   naxes(1) = NLTEspec%Nwaves
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
   CALL ftpkys(unit, "UNIT", "nm", comment, status)
   CALL ftpprd(unit,group,fpixel,naxes(1),NLTEspec%lambda,status)

  !  Close the file and free the unit number.
  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     CALL print_error(status)
  endif

 RETURN
 END SUBROUTINE WRITE_FLUX
  
  SUBROUTINE writeWavelength()
  ! --------------------------------------------------------------- !
   ! Write wavelength grid build with each transition of each atom
  ! --------------------------------------------------------------- !
   integer :: unit, EOF = 0, blocksize, naxes(1), naxis,group, bitpix, fpixel
   logical :: extend, simple
   character(len=6) :: comment="VACUUM"
   double precision :: lambda_vac(NLTEspec%Nwaves)
   
   if (.not.allocated(NLTEspec%lambda).or.&
       .not.NLTEspec%write_wavelength_grid) RETURN !
   write(*,*) " -> Writing wavelengths to ", WAVES_FILE       
   
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
  
  SUBROUTINE writeContrib()
  ! --------------------------------------------------------------- !
   ! Write the different contributions to opacities to fits file
  ! --------------------------------------------------------------- !
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(4) :: naxes
  integer :: group,fpixel,nelements, i, xcenter
  integer :: la, Nred, Nblue, kr, kc, m
  logical :: simple, extend
  character(len=6) :: comment=""
  
  write(*,*) "Writing opacities"
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
   status=0
   CALL ftgiou (unit,status)

   !  Create the new empty FITS file.
   blocksize=1
   CALL ftinit(unit,trim(OPAC_CONTRIB_FILE),blocksize,status)

   simple=.true.
   extend=.true.
   group=1
   fpixel=1

   bitpix=-64
   naxis=4
   naxes(1)=NLTEspec%Nwaves
   naxes(2)=NLTEspec%atmos%Nspace
   naxes(3)=RT_n_incl
   naxes(4)=RT_n_az
   nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
  
  RETURN
  END SUBROUTINE

END MODULE spectrum_type

