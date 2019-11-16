! ---------------------------------------------------------------------- !
! This module writes to separate files the contributions of each atom.
! ---------------------------------------------------------------------- !
MODULE write_opacity

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type
 use math, only : locate

 !MCFOST's original modules
 use fits_utils, only : print_error
 use grid
 use utils, only : appel_syst
 use messages
 use mcfost_env, only : dp

 IMPLICIT NONE

 character(len=15), parameter :: Jnu_File = "Jnu.fits.gz"
 character(len=15), parameter :: taur_File = "taur.fits.gz"



 CONTAINS


 SUBROUTINE write_Jnu
 ! ------------------------------------ !
 ! write the mean radiation field
 ! ------------------------------------ !
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements

  if (.not.atmos%electron_scattering) return
  !get unique unit number
  CALL ftgiou(unit,EOF)
  blocksize=1
  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .true.
  bitpix = -64

  if (lVoronoi) then
   	naxis = 2
   	naxes(1) = NLTEspec%Nwaves
   	naxes(2) = atmos%Nspace ! equivalent n_cells
   	nelements = naxes(1)*naxes(2)
  else
   	if (l3D) then
    	naxis = 4
    	naxes(1) = NLTEspec%Nwaves
    	naxes(2) = n_rad
   	 	naxes(3) = 2*nz
    	naxes(4) = n_az
    	nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
   	else
    	naxis = 3
    	naxes(1) = NLTEspec%Nwaves
    	naxes(2) = n_rad
    	naxes(3) = nz
    	nelements = naxes(1) * naxes(2) * naxes(3)
   	end if
  end if

  CALL ftinit(unit,trim(Jnu_file),blocksize,EOF)
  !  Initialize parameters about the FITS image
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "W.m-2.Hz^-1.sr^-1", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,atmos%nHmin,EOF)


  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)

  if (EOF > 0) CALL print_error(EOF)

 RETURN
 END SUBROUTINE write_Jnu
 
 SUBROUTINE write_taur(lambda)
 ! ----------------------------------------------------------- !
 ! write the continuum optical depth at reference
 ! wavelength
 ! ----------------------------------------------------------- !
  real(kind=dp), intent(in) :: lambda
  integer :: iref
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  real(kind=dp), dimension(atmos%Nspace) :: taur
  
  taur(:) = 0.0_dp
  iref = locate(NLTEspec%lambda, lambda)
  
  taur(:) = NLTEspec%AtomOpac%Kc(iref,:,1)
  if (atmos%NactiveAtoms > 0) then
   taur(:) = taur(:) + NLTEspec%AtomOpac%Kc_nlte(iref,:)
  endif

  !get unique unit number
  CALL ftgiou(unit,EOF)
  blocksize=1
  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .true.
  bitpix = -64

  if (lVoronoi) then
   	naxis = 1
   	naxes(1) = atmos%Nspace ! equivalent n_cells
   	nelements = naxes(1)
  else
   	if (l3D) then
    	naxis = 3
    	naxes(1) = n_rad
   	 	naxes(2) = 2*nz
    	naxes(3) = n_az
    	nelements = naxes(1) * naxes(2) * naxes(3)
   	else
    	naxis = 2
    	naxes(1) = n_rad
    	naxes(2) = nz
    	nelements = naxes(1) * naxes(2)
   	end if
  end if
  

  CALL ftinit(unit,trim(taur_file),blocksize,EOF)
  !  Initialize parameters about the FITS image
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", " ", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,taur,EOF)


  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)

  if (EOF > 0) CALL print_error(EOF)
  
 RETURN
 END SUBROUTINE write_taur

END MODULE write_opacity
