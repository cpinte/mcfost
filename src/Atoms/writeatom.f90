! ---------------------------------------------------------------------- !
! This module writes to separate files the contribution of each atom.
! - n
! - eta and chi (only for cheking)
! -> Only if the atom is ACTIVE:
! - Rij and Rji
! - rho = psi/phi (only for cheking)
! 
! Atomic data are written to atom%dataFile
! Electronic density is also written in a separate file
! ---------------------------------------------------------------------- !
MODULE writeatom

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type
 
 !MCFOST's original modules
 use fits_utils, only : print_error
 use grid

 IMPLICIT NONE 
 
 character(len=15), parameter :: neFile = "ne.fits.gz"
 character(len=15), parameter :: nHFile = "nHtot.fits.gz"

 
 CONTAINS
 
 SUBROUTINE writeElectron()
 ! ------------------------------------ !
 ! write Electron density if computed 
 ! by the code.
 ! The electron density, ne, is in m^-3
 ! ------------------------------------ !
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,trim(neFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 !??
  fpixel = 1
  extend = .false. !??
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
    nelements = naxes(1) * naxes(2) * naxes(3) ! != n_cells ? should be
   else
    naxis = 2
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2) ! should be n_cells also
   end if
  end if
  ! write(*,*) 'yoyo2', n_rad, nz, n_cells, atmos%Nspace
  !  Write the required header keywords.

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,atmos%ne,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writeElectron
 
 SUBROUTINE writeHydrogenDensity()
 ! ------------------------------------ !
 ! write nHtot density for testing in m^-3
 ! ------------------------------------ !
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,trim(nHFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 !??
  fpixel = 1
  extend = .false. !??
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
    nelements = naxes(1) * naxes(2) * naxes(3) ! != n_cells ? should be
   else
    naxis = 2
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2) ! should be n_cells also
   end if
  end if
  ! write(*,*) 'yoyo2', n_rad, nz, n_cells, atmos%Nspace
  !  Write the required header keywords.

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,atmos%nHtot,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writeHydrogenDensity
 
 !!!!!!!!!! BUILDING !!!!!!!!!!!!!!!
 SUBROUTINE writeAtomData(atom)
 ! ------------------------------------ !
 ! write AtomType atom to atom%dataFile
 ! ------------------------------------ !
  type (AtomType), intent(in) :: atom
  integer :: unit, EOF = 0, blocksize, naxes(6), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements, Nlevel2
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,trim(atom%dataFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true.
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
    nelements = naxes(1) * naxes(2) * naxes(3) ! != n_cells ? should be
   else
    naxis = 2
    naxes(1) = n_rad
    naxes(2) = nz
    nelements = naxes(1) * naxes(2) ! should be n_cells also
   end if
  end if
  
  !temporary
  Nlevel2 = atom%Nlevel*atom%Nlevel
  ! First write collision
  !add 1 axis for ij
  if (lVoronoi) then
   naxes(2) = Nlevel2
  else
   if (l3D) then
    naxes(4) = Nlevel2
   else
    naxes(3) = Nlevel2
   end if
  end if
  CALL ftphpr(unit,simple,bitpix,naxis+1,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements*Nlevel2,atom%Ckij,EOF)
  !write Collision
  ! create new hdu
  !CALL ftcrhd(unit, status)
  !  Write the required header keywords.
  !CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)  
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)  
 
 RETURN
 END SUBROUTINE writeAtomData 
 
 !!!!!!!!!! BUILDING !!!!!!!!!!!!!!!
 SUBROUTINE readElectron(fileNe)
 ! ------------------------------------ !
 ! read Electron density from file
 ! ------------------------------------ !
  character(len=15), intent(in) :: fileNe
  integer :: EOF = 0, unit, blocksize, hdutype, anynull
  integer :: naxis, naxes(4)

  !get unique unit number
  CALL ftgiou(unit,EOF)
  ! open fits file in readmode'
  CALL ftopen(unit, TRIM(fileNe), 0, blocksize, EOF)
  
  CALL FTMAHD(unit,1,hdutype,EOF) !only one, the electron density
  if (lVoronoi) then
   CALL ftgkyj(unit, "NAXIS1", naxes(1), " ", EOF)
   if (naxes(1).ne.(atmos%Nspace)) then
      write(*,*) "Cannot read electron density (1)"
      RETURN
   end if
   RETURN !not implemented
  else
   if (l3D) then
    CALL ftgkyj(unit, "NAXIS1", naxes(1), " ", EOF)
    CALL ftgkyj(unit, "NAXIS2", naxes(2), " ", EOF)
    CALL ftgkyj(unit, "NAXIS3", naxes(3), " ", EOF)
    if ((naxes(1)*naxes(2)*naxes(3)).ne.(atmos%Nspace)) then
      write(*,*) "Cannot read electron density (2)"
      RETURN
    end if
    RETURN ! Not implemented
   else
    CALL ftgkyj(unit, "NAXIS1", naxes(1), " ", EOF)
    CALL ftgkyj(unit, "NAXIS2", naxes(2), " ", EOF)
    if ((naxes(1)*naxes(2)).ne.(atmos%Nspace)) then
      write(*,*) "Cannot read electron density (3)"
      RETURN
    end if
    CALL FTG2Dd(unit,1,-999,shape(atmos%ne),naxes(1),naxes(2),atmos%ne,anynull,EOF)
   end if
  end if 
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 RETURN
 END SUBROUTINE readElectron
 
 
 !!!!!!!!!! BUILDING !!!!!!!!!!!!!!!
 SUBROUTINE readAtomData(atom)
 ! ------------------------------------ !
 ! reads and fills atomic data from
 ! atomic file atom%dataFile.
 ! 
 ! Note that is that case, atom is filled
 ! with quantum data and with the name of
 ! the file to read populations from.
 ! ------------------------------------ !
  type (AtomType), intent(inout) :: atom
 
 RETURN
 END SUBROUTINE readAtomData

END MODULE writeatom