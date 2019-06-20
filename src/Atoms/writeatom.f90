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
MODULE writeatom !Futur write.... 

 use atmos_type, only : atmos
 use atom_type
 use spectrum_type
 
 !MCFOST's original modules
 use fits_utils, only : print_error
 use grid
 use utils, only : appel_syst
 use messages

 IMPLICIT NONE 
 
 !Compressed files may only be opened with read-only permission
 ! So if ne is iterated we re open ne.fits to write the nlte pops, so It cannot be
 ! compressed.
 !
 !
 character(len=15), parameter :: neFile = "ne.fits"
 character(len=15), parameter :: nHFile = "nHtot.fits.gz"
 character(len=50), parameter :: TemperatureFile = "Temperature.fits.gz"

 
 CONTAINS
 
 SUBROUTINE writeIteration()
 
 
 RETURN
 END SUBROUTINE writeIteration
 
 SUBROUTINE make_checkpoint()
 
 RETURN
 END SUBROUTINE make_checkpoint
 
 SUBROUTINE writeElectron(write_ne_nlte)
 ! ------------------------------------ !
 ! write Electron density if computed 
 ! by the code.
 ! The electron density, ne, is in m^-3
 ! ------------------------------------ !
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple, already_exists = .false.
  logical, optional :: write_ne_nlte
  integer :: nelements, nfirst
  
  if (present(write_ne_nlte)) already_exists = (write_ne_nlte==.true.) !append the file
  !get unique unit number
  CALL ftgiou(unit,EOF)
  blocksize=1
  simple = .true. !Standard fits
  group = 1 !??
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
  ! write(*,*) 'yoyo2', n_rad, nz, n_cells, atmos%Nspace
  !  Write the required header keywords.
  
  if (already_exists) then !exits, expect to write nlte, read lte ne, append the file and write
  	CALL ftopen(unit, TRIM(neFile), 1, blocksize, EOF) !1 stends for read and write
  	!Compressed files may only be opened with read-only permission, so neFile is not .gz
   	CALL ftcrhd(unit, EOF)
   	CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  	CALL ftpkys(unit, "UNIT", "m^-3", '(NLTE)', EOF)
  	CALL ftpprd(unit,group,fpixel,nelements,atmos%ne,EOF)
  else !not exist, expected not nlte, create and write
  	CALL ftinit(unit,trim(neFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
  	CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  	CALL ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
  	CALL ftpprd(unit,group,fpixel,nelements,atmos%ne,EOF)
  end if !alreayd_exists
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writeElectron
 
 SUBROUTINE writeHydrogenDensity()
 ! ------------------------------------ !
 ! write nHtot density m^-3
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
 
 SUBROUTINE writeTemperature()
 ! ------------------------------------ !
 ! T in K
 ! ------------------------------------ !
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL ftinit(unit,trim(TemperatureFile),blocksize,EOF)
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
  CALL ftpkys(unit, "UNIT", "K", ' ', EOF)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,atmos%T,EOF)
  
  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writeTemperature
 
 SUBROUTINE writePops(atom)
 ! -------------------------------------------- !
 ! write Atom populations.
 ! First, NLTE populations, then LTE populations.
 ! -------------------------------------------- !
  type (AtomType), intent(in) :: atom
  integer :: unit, EOF = 0, blocksize, naxes(5), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements, syst_status
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  blocksize=1
  CALL appel_syst('mkdir -p '//trim(atom%ID),syst_status) !create a dir for populations
  CALL ftinit(unit,trim(atom%ID)//"/"//trim(atom%dataFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 
  fpixel = 1
  extend = .true.
  bitpix = -64  

  naxes(1) = atom%Nlevel

  if (lVoronoi) then
   naxis = 2
   naxes(2) = atmos%Nspace ! equivalent n_cells
   nelements = naxes(1)*naxes(2)
  else
   if (l3D) then
    naxis = 4
    naxes(2) = n_rad
    naxes(3) = 2*nz
    naxes(4) = n_az
    nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
   else
    naxis = 3
    naxes(2) = n_rad
    naxes(3) = nz
    nelements = naxes(1) * naxes(2) * naxes(3)
   end if
  end if
  
  

  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", "n (m^-3)", ' ', EOF)
  !write NLTE pops
  CALL ftpprd(unit,group,fpixel,nelements,atom%n,EOF)
  !now LTE
  CALL ftcrhd(unit, EOF)
  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  CALL ftpkys(unit, "UNIT", "nstar (m^-3)", ' ', EOF)

  CALL ftpprd(unit,group,fpixel,nelements,atom%nstar,EOF)

  
  CALL ftclos(unit, EOF) !close
  CALL ftfiou(unit, EOF) !free
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE writePops
 
 SUBROUTINE readPops(atom)
 ! -------------------------------------------- !
 ! read Atom populations.
 ! First, NLTE populations, then LTE populations.
 ! -------------------------------------------- !
  type (AtomType), intent(in) :: atom
  integer :: unit, EOF = 0, blocksize, naxis,group, bitpix, fpixel
  logical :: extend, simple, anynull
  integer :: nelements, naxis2(4), Nl, naxis_found, hdutype
  character(len=256) :: some_comments
  
  !get unique unit number
  CALL ftgiou(unit,EOF)

  CALL ftopen(unit, trim(atom%ID)//"/"//TRIM(atom%dataFile), 0, blocksize, EOF)
  CALL FTMAHD(unit,1,hdutype,EOF) !move to HDU, atom%n is first

  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 
  fpixel = 1
  extend = .true.
  bitpix = -64  

  if (lVoronoi) then
   naxis = 2
   CALL ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, EOF)
   CALL ftgkyj(unit, "NAXIS1", Nl, some_comments, EOF)
   if (Nl /= atom%Nlevel) CALL Error(" Nlevel read from file different from actual model!")
   nelements = Nl
   CALL ftgkyj(unit, "NAXIS2", naxis_found, some_comments, EOF)
   nelements = nelements * naxis_found
   if (nelements /= atom%Nlevel * atmos%Nspace) &
    CALL Error (" Model read does not match simulation box")
  else
   if (l3D) then
    naxis = 4
    CALL ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, EOF)
    !if (naxis2 /= naxis) CALL Error(" Naxis read from file different from actual model!")
    CALL ftgkyj(unit, "NAXIS1", Nl, some_comments, EOF)
    if (Nl /= atom%Nlevel) CALL Error(" Nlevel read from file different from actual model!")
    nelements = Nl
    CALL ftgkyj(unit, "NAXIS2", naxis_found, some_comments, EOF)
    nelements = nelements * naxis_found
    CALL ftgkyj(unit, "NAXIS3", naxis_found, some_comments, EOF)
    nelements = nelements * naxis_found
    CALL ftgkyj(unit, "NAXIS4", naxis_found, some_comments, EOF)
    nelements = nelements * naxis_found
    if (nelements /= atom%Nlevel * atmos%Nspace) &
    CALL Error (" Model read does not match simulation box")
   else
    naxis = 3
    CALL ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, EOF)
    CALL ftgkyj(unit, "NAXIS1", Nl, some_comments, EOF)
    if (Nl /= atom%Nlevel) CALL Error(" Nlevel read from file different from actual model!")
    nelements = Nl
    CALL ftgkyj(unit, "NAXIS2", naxis_found, some_comments, EOF)
    nelements = nelements * naxis_found
    CALL ftgkyj(unit, "NAXIS3", naxis_found, some_comments, EOF)
    nelements = nelements * naxis_found
    if (nelements /= atom%Nlevel * atmos%Nspace) &
    CALL Error (" Model read does not match simulation box")
   end if
  end if

  CALL FTG2Dd(unit,1,-999,shape(atom%n),atom%Nlevel,atmos%Nspace,atom%N,anynull,EOF)

  CALL ftclos(unit, EOF) !close
  CALL ftfiou(unit, EOF) !free
  
  if (EOF > 0) CALL print_error(EOF)
 
 RETURN
 END SUBROUTINE readPops
 
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