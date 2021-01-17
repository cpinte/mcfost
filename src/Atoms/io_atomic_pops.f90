! ---------------------------------------------------------------------- !
! This module writes to separate files the populations of each atom.
! ---------------------------------------------------------------------- !
MODULE io_atomic_pops

 use atmos_type, only : ne, nHmin, nHtot, hydrogen, ntotal_atom
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
 character(len=15), parameter :: neFile = "ne.fits.gz"
 character(len=15), parameter :: nHminFile = "nHmin.fits.gz"
 character(len=15), parameter :: nHFile = "nHtot.fits.gz"

 !integer, parameter :: Nmax_kept_iter = 1 !keep the Nmax_kept_iter latest iterations.
 !the first is the oldest one


 CONTAINS
 
subroutine write_Hminus
	! ------------------------------------ !
	! write H- density if computed
	! by the code.
	! The electron density, ne, is in m^-3
	! ------------------------------------ !
	integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
	logical :: extend, simple
	integer :: nelements, nfirst,sys_status
	character(len=512) :: cmd
	
	
		!check if data file already exist, can be the case if initial solutions is OLD_POPULATIONS
	cmd = "ls "//trim(nHminFile)
	call appel_syst(cmd, sys_status)
	if (sys_status == 0) then !means the file exist
	
		cmd = "mv "//trim(nHminFile)//" "//"nHmin_old.fits.gz"
		call appel_syst(cmd, sys_status)
		if (sys_status /= 0) then 
			call error("Error in copying old ne!")
		endif
		
	endif

	!get unique unit number
	call ftgiou(unit,EOF)
	blocksize=1
	simple = .true. !Standard fits
	group = 1 !??
	fpixel = 1
	extend = .true.
	bitpix = -64

	if (lVoronoi) then
   		naxis = 1
   		naxes(1) = n_cells
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


	call ftinit(unit,trim(nHminFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
	call ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
	call ftpprd(unit,group,fpixel,nelements,nHmin,EOF)

	call ftclos(unit, EOF)
	call ftfiou(unit, EOF)

	if (EOF > 0) call print_error(EOF)

return
end subroutine write_Hminus


subroutine write_Electron
	! ------------------------------------ !
	! write Electron density if computed
	! by the code.
	! The electron density, ne, is in m^-3
	! ------------------------------------ !
	integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
	logical :: extend, simple
	integer :: nelements, nfirst,sys_status
	character(len=512) :: cmd
	
	
		!check if data file already exist, can be the case if initial solutions is OLD_POPULATIONS
	cmd = "ls "//trim(nefile)
	call appel_syst(cmd, sys_status)
	if (sys_status == 0) then !means the file exist
	
		cmd = "mv "//trim(nefile)//" "//"ne_old.fits.gz"
		call appel_syst(cmd, sys_status)
		if (sys_status /= 0) then 
			call error("Error in copying old ne!")
		endif
		
	endif

	!get unique unit number
	call ftgiou(unit,EOF)
	blocksize=1
	simple = .true. !Standard fits
	group = 1 !??
	fpixel = 1
	extend = .true.
	bitpix = -64

	if (lVoronoi) then
   		naxis = 1
   		naxes(1) = n_cells
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


	call ftinit(unit,trim(neFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
	call ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
	call ftpprd(unit,group,fpixel,nelements,ne,EOF)

	call ftclos(unit, EOF)
	call ftfiou(unit, EOF)

	if (EOF > 0) call print_error(EOF)

return
end subroutine write_Electron

subroutine read_electron(lelectron_read)
	! ------------------------------------ !
	! read Electron density from file
	! ------------------------------------ !
	logical, intent(out) :: lelectron_read
	integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
	logical :: extend, simple, anynull
	integer :: nelements, naxis2(4), sys_status, naxis_found, hdutype
	character(len=512) :: cmd, some_comments
			
	!check if data file already exist, otherwise lead and return .false.
	cmd = "ls "//trim(nefile)
	call appel_syst(cmd, sys_status)
	if (sys_status == 0) then !means the file exist
		write(*,*) " Reading electron density from previous calc."
		lelectron_read = .true.
	else
		write(*,*) " found no electron density file. Evaluating from read populations."
		lelectron_read = .false.
		return
	endif
  
	status = 0
	call ftgiou(unit,status)
	
	call ftopen(unit, trim(nefile), 0, blocksize, status)
	if (status > 0) then
		write(*,*) "Cannot open electron file! ", nefile
		call print_error(status)
	endif


	simple = .true. !Standard fits
	group = 1
	fpixel = 1
	extend = .true.
	bitpix = -64
  

	call ftmahd(unit,1,hdutype,status)
	if (status > 0) then
		write(*,*) "Cannot read ne! "
		call print_error(status)
		stop
	endif 

	if (lVoronoi) then
		naxis = 1
		call error("read_electron does not handled Voronoi grid yet!")
	else

		if (l3D) then
			naxis = 3
		else
			naxis = 2
		end if
	
		!Number of axis
		call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
		if (status > 0) then
			write(*,*) "error reading number of axis (naxis)"
			call print_error(status)
			stop
		endif
		
		!R
		call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
		if (status > 0) then
			write(*,*) "error reading nrad from file (naxis1)"
			call print_error(status)
			stop
		endif
		nelements = naxis_found
		
		!Z
		call ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
		if (status > 0) then
			write(*,*) "error reading nz from file (naxis2)"
			call print_error(status)
			stop
		endif
		nelements = nelements * naxis_found
		
		
		!phi axis ?
		if (l3D) then
			call ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
			if (status > 0) then
				write(*,*) "error reading naz from file (naxis3)"
				call print_error(status)
				stop
			endif
			nelements = nelements * naxis_found
		endif
		
		if (nelements /= n_cells) then
			write(*,*) " read_electron does not do interpolation yet!"
			call Error (" Model read does not match simulation box")
		endif

		call FTG2Dd(unit,group,-999,shape(ne),n_cells,1,ne,anynull,status)
		!call ftgpvd(unit,group,1,n_cells,-999,ne,anynull,status)
		
		if (status > 0) then
			write(*,*) "error reading ne density"
			call print_error(status)
			stop
		endif

	endif !lvoronoi

	call ftclos(unit, status) !close
	if (status > 0) then
		write(*,*) "error cannot close file in ", trim(nefile)
		call print_error(status)
		stop
	endif
	
	call ftfiou(unit, status) !free
	if (status > 0) then
		write(*,*) "error cannot free file unit!"
		call print_error(status)
! 		stop
	endif	
	
	!call FTG2Dd(unit,1,-999,shape(ne),naxes(1),naxes(2),ne,anynull,EOF)

    write(*,'("  -- min(ne)="(1ES20.7E3)" m^-3; max(ne)="(1ES20.7E3)" m^-3")') , minval(ne,mask=(icompute_atomRT>0)), maxval(ne)

return
end subroutine read_electron

subroutine write_convergence_map_electron(nstep, conv_map)
	! ------------------------------------ !
	! write convergence map for an atom atom
	! ------------------------------------ !
	integer, intent(in) :: nstep
	real(kind=dp), intent(in) :: conv_map(n_cells,nstep)
	integer :: unit, EOF = 0, blocksize, naxes(6), naxis,group, bitpix, fpixel
	logical :: extend, simple
	integer :: nelements, nfirst,sys_status
	character(len=512) :: cmd
	
	!get unique unit number
	call ftgiou(unit,EOF)
	blocksize=1
	simple = .true. !Standard fits
	group = 1 !??
	fpixel = 1
	extend = .true.
	bitpix = -64

	if (lVoronoi) then
   		naxis = 2
   		naxes(1) = n_cells
   		naxes(2) = nstep
   		nelements = naxes(1) * naxes(2)
	else
		if (l3D) then
			naxis = 4
			naxes(1) = n_rad
			naxes(2) = 2*nz
			naxes(3) = n_az
   			naxes(4) = nstep
			nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
		else
			naxis = 3
			naxes(1) = n_rad
			naxes(2) = nz
   			naxes(3) = nstep
			nelements = naxes(1) * naxes(2) * naxes(3)
		end if
	end if


	call ftinit(unit,"ne_convmap.fits.gz",blocksize,EOF)
  !  Initialize parameters about the FITS image
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  !write data
	call ftpprd(unit,group,fpixel,nelements,conv_map,EOF)

	call ftclos(unit, EOF)
	call ftfiou(unit, EOF)

	if (EOF > 0) call print_error(EOF)

return
end subroutine write_convergence_map_electron

subroutine write_convergence_map_atom(atom, nstep, conv_map)
	! ------------------------------------ !
	! write convergence map for an atom atom
	! ------------------------------------ !
	type (AtomType), intent(in) :: atom
	integer, intent(in) :: nstep
	real(kind=dp), intent(in) :: conv_map(n_cells,atom%Nlevel,nstep)
	integer :: unit, EOF = 0, blocksize, naxes(6), naxis,group, bitpix, fpixel
	logical :: extend, simple
	integer :: nelements, nfirst,sys_status
	character(len=512) :: cmd
	

	!get unique unit number
	call ftgiou(unit,EOF)
	blocksize=1
	simple = .true. !Standard fits
	group = 1 !??
	fpixel = 1
	extend = .true.
	bitpix = -64

	if (lVoronoi) then
   		naxis = 3
   		naxes(1) = n_cells
   		naxes(2) = atom%Nlevel
   		naxes(3) = nstep
   		nelements = naxes(1) * naxes(2) * naxes(3)
	else
		if (l3D) then
			naxis = 5
			naxes(1) = n_rad
			naxes(2) = 2*nz
			naxes(3) = n_az
   			naxes(4) = atom%Nlevel
   			naxes(5) = nstep
			nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4) * naxes(5)
		else
			naxis = 4
			naxes(1) = n_rad
			naxes(2) = nz
   			naxes(3) = atom%Nlevel
   			naxes(4) = nstep
			nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
		end if
	end if


	call ftinit(unit,trim(atom%ID)//"_convmap.fits.gz",blocksize,EOF)
  !  Initialize parameters about the FITS image
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  !write data
	call ftpprd(unit,group,fpixel,nelements,conv_map,EOF)

	call ftclos(unit, EOF)
	call ftfiou(unit, EOF)

	if (EOF > 0) call print_error(EOF)

return
end subroutine write_convergence_map_atom
 

 subroutine write_Hydrogen()
 ! ------------------------------------ !
 ! write nHtot density m^-3
 ! ------------------------------------ !
  integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements

  !get unique unit number
  call ftgiou(unit,EOF)

  blocksize=1
  call ftinit(unit,trim(nHFile),blocksize,EOF)
  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1 !??
  fpixel = 1
  extend = .false. !??
  bitpix = -64

  if (lVoronoi) then
   naxis = 1
   naxes(1) = n_cells ! equivalent n_cells
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
  ! write(*,*) 'yoyo2', n_rad, nz, n_cells, n_cells
  !  Write the required header keywords.

  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
  ! Additional optional keywords
  call ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
  !write data
  call ftpprd(unit,group,fpixel,nelements,nHtot,EOF)

  call ftclos(unit, EOF)
  call ftfiou(unit, EOF)

  if (EOF > 0) call print_error(EOF)

 return
 end subroutine write_Hydrogen

!  subroutine write_Hminus()
!  ! ------------------------------------ !
!  ! write H- density.
!  ! ------------------------------------ !
!   integer :: unit, blocksize, naxes(5), naxis,group, bitpix, fpixel
!   logical :: extend, simple
!   integer :: nelements, hdutype, status
!   character(len=20) :: popsF
! 
!   status = 0
!   !get unique unit number
!   call ftgiou(unit,status)
!   blocksize=1
!   simple = .true. !Standard fits
!   group = 1 !??
!   fpixel = 1
!   extend = .true.
!   bitpix = -64
!   
!   call ftinit(unit, trim(root_dir)//"/"//TRIM(nHminFile), blocksize, status)
!   if (status > 0) then
!       write(*,*) "Cannot create fits file ", nHminFile
!       call print_error(status)
!   endif
! 
!   if (lVoronoi) then
!    naxis = 1
!    naxes(1) = n_cells
!    nelements = naxes(1)
!   else
!    if (l3D) then
!     naxis = 3
!     naxes(1) = n_rad
!     naxes(2) = 2*nz
!     naxes(3) = n_az
!     nelements = naxes(1) * naxes(2) * naxes(3)
!    else
!     naxis = 2
!     naxes(1) = n_rad
!     naxes(2) = nz
!     nelements = naxes(1) * naxes(2)
!    end if
!   end if
!   
! 
! 
!   
!     call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!     if (status > 0) then
!       write(*,*) "Error writing LTE H- to file"
!       call print_error(status)
!     endif
! 
!     ! Additional optional keywords
!     call ftpkys(unit, "UNIT", "m^-3", "(LTE)", status)
!     if (status > 0) then
!       write(*,*) "Error writing LTE H- to file (2)"
!       call print_error(status)
!     endif
!     !write data
!     call ftpprd(unit,group,fpixel,nelements,nHmin,status)
!     if (status > 0) then
!       write(*,*) "Error writing LTE H- to file (3)"
!       call print_error(status)
!     endif
! 
! 
!   call ftclos(unit,status) !close
!   call ftfiou(unit,status) !free
! 
!   if (status > 0) call print_error(status)
! 
!  return
!  end subroutine write_Hminus



!  subroutine writePops(atom, count)
!  ! ----------------------------------------------------------------- !
!  ! write Atom populations. The file for writing exists already
!  ! First, NLTE populations, then LTE populations.
!  ! It is not a compressed file (.fits only)
!  ! so that populations after ith iteration can be
!  ! append to the file.
!  ! Keep in memory the three latest iterations
!  ! and overwrite them. At the end write the final
!  ! NLTE pops and then the LTE pops.
!  !
!  ! Contains:
!  !   1st oldest iter among the three last iter: count = 3
!  !   2nd 									  : count = 2
!  !   3rd youngest iter among the three last iter : count = 1
!  !      count = 0
!  !   atom%n converged (can be equal to atom%n(iter==3)
!  !   atom%nstar
!  !
!  !  count is decreasing
!  ! ----------------------------------------------------------------- !
!   type (AtomType), intent(inout) :: atom
!   integer, intent(in) :: count
!   integer :: unit, EOF = 0, blocksize, naxes(5), naxis,group, bitpix, fpixel
!   logical :: extend, simple, lte_only
!   integer :: nelements, syst_status, hdutype
!   character(len=20) :: popsF
! 
!   lte_only = .not.atom%active
!   if (count < 0) lte_only = .true. !force it in case it is active
! 
!   if (count > Nmax_kept_iter) return
! 
!   !get unique unit number
!   call ftgiou(unit,EOF)
!   call ftopen(unit, TRIM(atom%dataFile), 1, blocksize, EOF) !1 stends for read and write
! 
!   blocksize=1
! 
!   !  Initialize parameters about the FITS image
!   simple = .true. !Standard fits
!   group = 1
!   fpixel = 1
!   extend = .true.
!   bitpix = -64
! 
!   naxes(1) = atom%Nlevel
! 
!   if (lVoronoi) then
!    naxis = 2
!    naxes(2) = n_cells ! equivalent n_cells
!    nelements = naxes(1)*naxes(2)
!   else
!    if (l3D) then
!     naxis = 4
!     naxes(2) = n_rad
!     naxes(3) = 2*nz
!     naxes(4) = n_az
!     nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
!    else
!     naxis = 3
!     naxes(2) = n_rad
!     naxes(3) = nz
!     nelements = naxes(1) * naxes(2) * naxes(3)
!    end if
!   end if
! 
!   if (count>0) then !iterations
!     if (count <= Nmax_kept_iter) then
!     call FTMAHD(unit,count,hdutype,EOF)
!   	call ftpprd(unit,group,fpixel,nelements,atom%n,EOF)
!   	endif
!   else if (count<=0) then !final solutions + LTE
!    if (lte_only) then
!     call FTMAHD(unit,1,hdutype,EOF) !only one if lte_only
!   	call ftpprd(unit,group,fpixel,nelements,atom%nstar,EOF)
!    else !at least two
!     call FTMAHD(unit,Nmax_kept_iter+1,hdutype,EOF)
!   	call ftpprd(unit,group,fpixel,nelements,atom%n,EOF)
!     call FTMAHD(unit,Nmax_kept_iter+2,hdutype,EOF)
!   	call ftpprd(unit,group,fpixel,nelements,atom%nstar,EOF)
!    end if
!   end if
! 
!   call ftclos(unit, EOF) !close
!   call ftfiou(unit, EOF) !free
! 
!   if (EOF > 0) call print_error(EOF)
! 
!  return
!  end subroutine writePops
!-> debug needed above

!building
!  subroutine read_departure_bfactor(atom)
!  ! -------------------------------------------- !
!  ! -------------------------------------------- !
!   type (AtomType), intent(in) :: atom
!   integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
!   logical :: extend, simple, anynull
!   integer :: nelements, naxis2(4), Nl, naxis_found, hdutype, l, icell
!   character(len=256) :: some_comments
!   
!   
!   !get unique unit number
!   status = 0
!   call ftgiou(unit,status)
! 
!   call ftopen(unit, TRIM(atom%dataFile), 0, blocksize, status)
!   if (status > 0) then
!       write(*,*) "Read departure fitsError (1) "
!       call print_error(status)
!   endif
! 
!   !  Initialize parameters about the FITS image
!   simple = .true. !Standard fits
!   group = 1
!   fpixel = 1
!   extend = .true.
!   bitpix = -64
! 
!   if (lVoronoi) then
!    call error ("read_departure: Voronoi ot supported yet")
!   else
! ! we always read 3D if not Voronoi ? We write 1D models as (1, 1, nr, nl)
! !											  2D models as (1, nz, nr, nl)
! !											  3D models as (naz, nz, nr, nl)
! !    if (l3D) then
!     naxis = 4
!     call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
!     if (status > 0) then
!       write(*,*) "Read departure fitsError (2) "
!       call print_error(status)
!     endif
!     !if (naxis2 /= naxis) call Error(" Naxis read from file different from actual model!")
!     call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
!     if (status > 0) then
!       write(*,*) "Read departure fitsError (3) "
!       call print_error(status)
!     endif
! !     if (naxis_found /= n_az) then 
! !       call Error(" (3D) wrong naz read from file different from actual model!")
! !     endif
!     nelements = naxis_found
!     call ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
!     if (status > 0) then
!       write(*,*) "Read departures fitsError (4) "
!       call print_error(status)
!     endif
! !     if (naxis_found /= nz) then 
! !       call Error(" (3D) wrong nz read from file different from actual model!")
! !     endif
!     nelements = nelements * naxis_found
!     call ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
!     if (status > 0) then
!       write(*,*) "Read departure fitsError (5) "
!       call print_error(status)
!     endif
! !     if (naxis_found /= n_rad) then 
! !       call Error(" (3D) wrong nrad read from file different from actual model!")
! !     endif
!     nelements = nelements * naxis_found
!     call ftgkyj(unit, "NAXIS4", naxis_found, some_comments, status)
!     if (status > 0) then
!       write(*,*) "Read departure fitsError (6) "
!       call print_error(status)
!     endif
! !     if (naxis_found /= atom%Nlevel) then 
! !       call Error(" (3D) wrong Nlevel read from file different from actual model!")
! !     endif
!     Nl = atom%Nlevel
!     nelements = nelements * naxis_found
!     if (nelements /= atom%Nlevel * n_cells) call Error (" Model read does not match simulation box")
!   end if
! 
!   call FTG2Dd(unit,1,-999,shape(atom%b),atom%Nlevel,n_cells,atom%b,anynull,status)
! 
!   do icell=1,n_cells
!    write(*,*) '--------------------'
!    write(*,*) "b=", (atom%b(l, icell), l=1,atom%Nlevel)
!    write(*,*) '--------------------'
!   enddo
! 
!   if (status > 0) then
!       write(*,*) "Read departure fitsError (7) "
!       call print_error(status)
!   endif
! 
!   call ftclos(unit, status) !close
!   if (status > 0) then
!       write(*,*) "Read departure fitsError (8) "
!       call print_error(status)
!   endif
!   call ftfiou(unit, status) !free
!   if (status > 0) then
!       write(*,*) "Read departure fitsError (9) "
!       call print_error(status)
!   endif
! 
!  return
!  end subroutine  read_departure_bfactor

subroutine prepare_check_pointing()
	integer :: sys_status
	character(len=512) :: cmd
	
	cmd = "mkdir iterations"
	call appel_syst(cmd, sys_status)
	if (sys_status /= 0) then 
		call Warning("Iterations folder already exists!!")
	endif

return
end subroutine
 
subroutine write_pops_atom(atom,iter,step)
	! ----------------------------------------------------------------- !
	! write Atom populations
	! First, NLTE populations, then LTE populations.
	! ----------------------------------------------------------------- !
	type (AtomType), intent(in) :: atom
	integer,  optional :: iter, step
	integer :: unit, blocksize, naxes(5), naxis,group, bitpix, fpixel
	logical :: extend, simple, lte_only
	integer :: nelements, hdutype, status, k, sys_status
	character(len=512) :: cmd, popsF
	character(len=10000) :: step_c, iter_c

	!lte_only = .not.atom%active
	lte_only = (.not.atom%NLTEpops) !can be true if atom ACTIVE or INITIAL_SOLUTION==OLD_POPULATIONS
	!because maybe NLTE pops are read but atom is considered as PASSIVE (i.e., no nlte iterations)
  
	status = 0
	!get unique unit number
	call ftgiou(unit,status)
	blocksize=1
	simple = .true. !Standard fits
	group = 1 !??
	fpixel = 1
	extend = .true.
	bitpix = -64
	
	popsF = trim(atom%ID)//".fits.gz"
	
	
	if (present(iter)) then
	!suppose that the file does not exist
		write(iter_c, '(i0)') iter
		if (.not.present(step)) step=0
		write(step_c, '(i0)') step
		!folder iterations must exist
		popsF = "iterations/"//trim(atom%ID)//'_iter'//trim(iter_c)//'_step'//trim(step_c)//'.fits.gz'
		
	else
	!check if data file already exist, can be the case if initial solutions is OLD_POPULATIONS
		cmd = "ls "//popsF!trim(atom%ID)//".fits.gz"
		call appel_syst(cmd, sys_status)
		if (sys_status == 0) then !means the file exist
	
			cmd = "mv "//trim(atom%ID)//".fits.gz"//" "//trim(atom%ID)//"_oldpop.fits.gz"
			call appel_syst(cmd, sys_status)
			if (sys_status /= 0) then 
				call error("Error in copying old pops!")
			endif
		
		endif
	endif
  
	call ftinit(unit, trim(root_dir)//"/"//TRIM(popsF), blocksize, status)
	if (status > 0) then
      write(*,*) "Cannot create fits file ", popsF
      call print_error(status)
	endif
	naxes(1) = atom%Nlevel

	if (lVoronoi) then
		naxis = 2
		naxes(2) = n_cells 
		nelements = naxes(1)*naxes(2)
	else
		if (l3D) then
			naxis = 4
			naxes(2) = n_rad
			naxes(3) = 2*nz
			naxes(4) = n_az
			nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
		else
			naxis = 3!4
			naxes(2) = n_rad
			naxes(3) = nz
! 			naxes(4) = 1 !write naz which is one
			nelements = naxes(1) * naxes(2) * naxes(3)
		end if
	end if


	if (lte_only) then
  
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		if (status > 0) then
			write(*,*) "Error writing LTE pops to file"
			call print_error(status)
		endif

    ! Additional optional keywords
		call ftpkys(unit, "UNIT", "m^-3", "(LTE)", status)
		if (status > 0) then
			write(*,*) "Error writing LTE pops to file (2)"
			call print_error(status)
		endif
    !write data
		call ftpprd(unit,group,fpixel,nelements,atom%nstar,status)
		if (status > 0) then
			write(*,*) "Error writing LTE pops to file (3)"
			call print_error(status)
		endif
	else !NLTE + LTE
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		if (status > 0) then
			write(*,*) "Error writing nLTE pops to file"
			call print_error(status)
		endif
    ! Additional optional keywords
		call ftpkys(unit, "UNIT", "m^-3", "(NLTE) ", status)
		if (status > 0) then
			write(*,*) "Error writing nLTE pops to file (2)"
			call print_error(status)
		endif

		call ftpprd(unit,group,fpixel,nelements,atom%n,status)
		if (status > 0) then
			write(*,*) "Error writing nLTE pops to file (3)"
			call print_error(status)
		endif
    
		call ftcrhd(unit, status)
		if (status > 0) then
			write(*,*) "Error writing nLTE pops to file (4)"
			call print_error(status)
		endif
  
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		if (status > 0) then
			write(*,*) "Error writing LTE pops to file (5)"
			call print_error(status)
		endif
    ! Additional optional keywords
		call ftpkys(unit, "UNIT", "m^-3", "(LTE)",status)
		if (status > 0) then
			write(*,*) "Error writing LTE pops to file (6)"
			call print_error(status)
		endif
    !write data
		call ftpprd(unit,group,fpixel,nelements,atom%nstar,status)
		if (status > 0) then
			write(*,*) "Error writing LTE pops to file (7)"
			call print_error(status)
		endif
  	
	end if  


	call ftclos(unit,status) !close
	call ftfiou(unit,status) !free

	if (status > 0) call print_error(status)

return
end subroutine write_pops_atom
 
subroutine read_pops_atom(atom)
	! ---------------------------------------------------------- !
	! read Atom populations.
	! First, NLTE populations, then LTE populations.
	!
	! To DO:
	! If the read populations are on a different grid than the
	! actual grid, they are interpolated. 
	! However: there is not extrapolation at the moment.
	! ---------------------------------------------------------- !
	type (AtomType), intent(inout) :: atom
	integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
	logical :: extend, simple, anynull, show_warning = .true.
	integer :: nelements, naxis2(4), Nl, naxis_found, hdutype, l, icell
	character(len=256) :: some_comments, popsF
  
	status = 0
	call ftgiou(unit,status)
	
	popsF = TRIM(atom%ID)//".fits.gz"

	call ftopen(unit, trim(popsF), 0, blocksize, status)
	if (status > 0) then
		write(*,*) "Cannot open file with populations! ", atom%id
		call print_error(status)
	endif


	simple = .true. !Standard fits
	group = 1
	fpixel = 1
	extend = .true.
	bitpix = -64
  
	!NLTE
	call ftmahd(unit,1,hdutype,status)
	if (status > 0) then
		write(*,*) "Cannot read nlte populations ! ", atom%id
		call print_error(status)
		stop
	endif 

	if (lVoronoi) then
		naxis = 2
		call error("read_pops_atom does not handled Voronoi grid yet!")
	else

		if (l3D) then
			naxis = 4
		else
			!no naz
			naxis = 3!4
		end if
	
		!Number of axis
		call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
		!naxis2 = (Nlevel,nrad, naz, nphi) = (Naxis1,Naxis2,Naxis3,Naxis4)
		!naxis4 is 0 if not l3D
		if (status > 0) then
			write(*,*) "error reading number of axis (naxis)"
			call print_error(status)
			stop
		endif
		
		!Nlevel
		call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
		if (status > 0) then
			write(*,*) "error reading Nlevel from file (naxis1)"
			call print_error(status)
			stop
		endif
		
		if (naxis_found /= atom%Nlevel) then
			if (naxis_found /= naxis2(1)) then
				write(*,*) "Nlevel read does not match atom !", atom%Nlevel
				stop
			endif
		endif
		nelements = naxis_found
		
		!R
		call ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
		if (status > 0) then
			write(*,*) "error reading nrad from file (naxis2)"
			call print_error(status)
			stop
		endif
! 		if (naxis_found /= nrad) then
! 			write(*,*) "nrad read does not match grid !", naxis_found
! 			stop
! 		endif
		nelements = nelements * naxis_found
		
		!z
		call ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
		if (status > 0) then
			write(*,*) "error reading nz from file (naxis3)"
			call print_error(status)
			stop
		endif
! 		if (naxis_found /= nz) then
! 			write(*,*) "nz read does not match grid !", naxis_found
! 			stop
! 		endif
		nelements = nelements * naxis_found
		
		!phi axis ?
		if (l3D) then
			call ftgkyj(unit, "NAXIS4", naxis_found, some_comments, status)
			if (status > 0) then
				write(*,*) "error reading naz from file (naxis4)"
				call print_error(status)
				stop
			endif
! 			if (naxis_found /= n_az) then
! 				write(*,*) "n_az read does not match grid !", naxis_found
! 				stop
! 			endif
			nelements = nelements * naxis_found
		endif
		
		if (nelements /= atom%Nlevel * n_cells) then
			write(*,*) " read_pops_atom does not do interpolation yet!"
			call Error (" Model read does not match simulation box")
		endif

		!READ NLTE
		call FTG2Dd(unit,1,-999,shape(atom%n),atom%Nlevel,n_cells,atom%n,anynull,status)
		if (status > 0) then
			write(*,*) "error reading non-LTE populations ", atom%id
			call print_error(status)
			stop
		endif


    	!now LTE
		call FTMAHD(unit,2,hdutype,status)
		if (status > 0) then
			write(*,*) "error opening LTE hdu ", atom%id
			call print_error(status)
			stop
		endif
    
		atom%nstar(:,:) = 0.0_dp
		call FTG2Dd(unit,1,-999,shape(atom%nstar),atom%Nlevel,n_cells,atom%nstar,anynull,status)
		if (status > 0) then
			write(*,*) "error reading LTE populations ", atom%id
			call print_error(status)
			stop
		endif

    
	endif !lvoronoi

	call ftclos(unit, status) !close
	if (status > 0) then
		write(*,*) "error cannot close file in ", trim(popsF)
		call print_error(status)
		stop
	endif
	
	call ftfiou(unit, status) !free
	if (status > 0) then
		write(*,*) "error cannot free file unit!"
		call print_error(status)
! 		stop
	endif

!     write(*,*) " min/max pops for each level:"
    do l=1,atom%Nlevel
    	!-> correct edge effect due to interpolation?
    	!pre-smoothing the data could help or a better interpolation.
    	!This steps does not cost anything but it requires a careful knowledge about the data
    	!and how this affects the flux.
    	!Alternatively those region could be adressed in another way (interpolation, LTE values, non-LTE solution...)
		do icell=1,n_cells
			if (icompute_atomRT(icell)) then!if on the exact grid (computed) is not empty but the interpolations are, empty the grid point.
				if ( (atom%n(l,icell) <= 0.0_dp) .or. (atom%nstar(l,icell) <= 0.0_dp)) then
					icompute_atomRT(icell) = 0
					if (show_warning) then
						call warning(" ++ Accommodating grid to match interpolated data...")
						show_warning = .false.
					endif
				endif
			endif
		enddo
    	write(*,"('Level #'(1I3))") l
    	write(*,'("  -- min(n)="(1ES20.7E3)" m^-3; max(n)="(1ES20.7E3)" m^-3")') , minval(atom%n(l,:),mask=(icompute_atomRT>0)), maxval(atom%n(l,:))
    	write(*,'("  -- min(nstar)="(1ES20.7E3)" m^-3; max(nstar)="(1ES20.7E3)" m^-3")')  minval(atom%nstar(l,:),mask=(icompute_atomRT>0)), maxval(atom%nstar(l,:))
	enddo


return
end subroutine read_pops_atom

end MODULE io_atomic_pops
