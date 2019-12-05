! ---------------------------------------------------------------------- !
! This module writes to separate files the contributions of each atom.
! ---------------------------------------------------------------------- !
MODULE write_opacity

 use atmos_type, only : atmos, Hydrogen, Helium
 use atom_type
 use spectrum_type
 use math, only : locate, Linear_1D, bezier2_interp

 use constant
 use Planck
 use hydrogen_opacities
 use rayleigh_scattering

 !MCFOST's original modules
 use fits_utils, only : print_error, cfitsWrite
 use grid
 use utils, only : appel_syst
 use messages
 use mcfost_env, only : dp

 IMPLICIT NONE

 character(len=15), parameter :: Jnu_File = "Jnu.fits.gz", Jnu_File_ascii="Jnu_no_interp.s"
 character(len=15), parameter :: taur_File = "taur.fits.gz"



 CONTAINS
 
 SUBROUTINE read_Jnu_ascii()
 ! -------------------------------------------- !
 ! read Jnu from previous calculations
 ! This is the Jnu computed in Iterate_Jnu() if 
 ! we never enter in the NLTE loop.
 ! It is then interpolated on the whole grid.
 ! The bounds of lambda_jnu should contain
 ! the actual wavelength grid.
 ! The only requierment is that Ncell_j == n_cells
 ! I am not doing cells interpolations nor I planning
 ! to do it.
 ! -------------------------------------------- !
  integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
  logical :: extend, simple, anynull
  integer :: nelements, naxis2(4), Nl, naxis_found, hdutype, l, icell, ncells_j, nlam_j
  character(len=256) :: some_comments
  real(kind=dp), dimension(:), allocatable :: Jtmp, xtmp
  
  !get unique unit number
  status = 0
  CALL ftgiou(unit,status)
  
  write(*,*) " Reading Jnu from a previous calculation..."

  open(unit,file=Jnu_file_ascii,status="old")
  read(unit,*) ncells_j, nlam_j
  
  allocate(Jtmp(nlam_j), xtmp(nlam_j))
  Jtmp = 0.
  
  if (ncells_j /= n_cells) then
   write(*,*) "reading ", ncells_j," cells in J file,  but actual grid has ", n_cells, " cells" 
  endif

   write(*,*) "    -> interpolating on RT grid..."

  do icell=1, n_cells
    !unfortunatelly reading xtmp n_cells time
  	do l=1, nlam_j
     read(unit, *) xtmp(l), Jtmp(l)
    enddo
    CALL bezier2_interp(Nlam_j, xtmp, Jtmp, NLTEspec%Nwaves, NLTEspec%lambda, NLTEspec%Jc(:,icell))
    xtmp = 0.; Jtmp = 0.
  enddo
   write(*,*) "    -> ..done"

      
  close(unit)
    
  CALL ftfiou(unit, status) !free
  if (status > 0) then
      write(*,*) "Read_Jnu_ascii free unit error "
      CALL print_error(status)
  endif
  
  write(*,*) " ..done"

  
  deallocate(Jtmp, xtmp)
  
 RETURN
 END SUBROUTINE read_Jnu_ascii

 SUBROUTINE read_Jnu()
 ! -------------------------------------------- !
 ! read Jnu from previous calculations
 ! ATM assumes:
 !    1) that only Jc is written
 !    2) Jnu is on the exact same wavelength grid
 !       and spatial grid as NLTEspec%I
 ! -------------------------------------------- !
  integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
  logical :: extend, simple, anynull
  integer :: nelements, naxis2(4), Nl, naxis_found, hdutype, l, icell
  character(len=256) :: some_comments
  real(kind=dp), dimension(:), allocatable :: Jtmp
  
  !get unique unit number
  status = 0
  CALL ftgiou(unit,status)

  CALL ftopen(unit, TRIM(Jnu_File), 0, blocksize, status)
  if (status > 0) then
      write(*,*) "Read_Jnu fitsError (1) "
      CALL print_error(status)
  endif

  !  Initialize parameters about the FITS image
  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .true.
  bitpix = -64
  

!   CALL FTMAHD(unit,1,hdutype,status)
!    if (status > 0) then
!       write(*,*) "Readpops cannot move to first HDU "
!       CALL print_error(status)
!   endif 


  if (lVoronoi) then
   CALL error("read_Jnu: Voronoi not supported yet")
  else
   if (l3D) then
    naxis = 4
    CALL ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (2) "
      CALL print_error(status)
    endif !the first one is the number of axis in the fits right ?
    !!nelements = naxis_found
    
    
    
    CALL ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (3) "
      CALL print_error(status)
    endif
    Nl = naxis_found

    CALL ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (4) "
      CALL print_error(status)
    endif
    nelements = naxis_found
    CALL ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Readpops fitsError (5) "
      CALL print_error(status)
    endif
    nelements = nelements * naxis_found
    CALL ftgkyj(unit, "NAXIS4", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (6) "
      CALL print_error(status)
    endif
    nelements = nelements * Nl *  naxis_found
   else
   
    naxis = 3
    CALL ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (2) "
      CALL print_error(status)
    endif
    !nelements = naxis_found
    
    
    CALL ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (3) "
      CALL print_error(status)
    endif
    Nl = naxis_found
    
    
    CALL ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (3) "
      CALL print_error(status)
    endif
    nelements = naxis_found
    CALL ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
    if (status > 0) then
      write(*,*) "Read_Jnu fitsError (4) "
      CALL print_error(status)
    endif
    nelements = nelements * naxis_found * Nl
   endif !l3D
  end if
!      allocate(Jtmp(Nl, n_cells),stat=status)
!     if (status > 0)
!      Call error("Cannot allocate Jtmp read_Jnu")
!     endif

	if (Nl /= NLTEspec%Nwaves) then
	 CALL Warning("Jnu size does not match actual grid")
	endif

    CALL FTG2Dd(unit,1,-999,shape(NLTEspec%Jc),Nl,n_Cells,NLTEspec%Jc,anynull,status)
    if (status > 0) then
      write(*,*) "Read_Jnu cannot read Jnu "
      CALL print_error(status)
    endif


  CALL ftclos(unit, status) !close
  if (status > 0) then
      write(*,*) "Read_Jnu fitsError (8) "
      CALL print_error(status)
  endif
  
  CALL ftfiou(unit, status) !free
  if (status > 0) then
      write(*,*) "Read_Jnu fitsError (9) "
      CALL print_error(status)
  endif

 RETURN
 END SUBROUTINE read_Jnu

 SUBROUTINE write_Jnu
 ! ------------------------------------ !
 ! write the mean radiation field
 ! ------------------------------------ !
  integer :: unit, status = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements

  if (.not.atmos%electron_scattering) return
  !get unique unit number
  CALL ftgiou(unit,status)
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

  CALL ftinit(unit,trim(Jnu_file),blocksize,status)
  
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  ! Additional optional keywords
  CALL ftpkys(unit, "Jc(nu)", "W.m-2.Hz^-1.sr^-1", ' ', status)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%Jc,status)
  
  if (allocated(NLTEspec%J)) then
   CALL ftcrhd(unit, status)

   !  Initialize parameters about the FITS image
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
   ! Additional optional keywords
   CALL ftpkys(unit, "J(nu)", "W.m-2.Hz^-1.sr^-1", ' ', status)
   !write data
   CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%J,status)
  endif


  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  if (status > 0) CALL print_error(status)

 RETURN
 END SUBROUTINE write_Jnu
 
 SUBROUTINE write_taur(lambda)
 ! ----------------------------------------------------------- !
 ! write the continuum optical depth at reference
 ! wavelength
 ! ----------------------------------------------------------- !
  real(kind=dp), intent(in) :: lambda
  integer :: iref
  integer :: unit, status = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
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
  CALL ftgiou(unit,status)
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
  

  CALL ftinit(unit,trim(taur_file),blocksize,status)
  !  Initialize parameters about the FITS image
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  ! Additional optional keywords
  CALL ftpkys(unit, "UNIT", " ", ' ', status)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,taur,status)


  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  if (status > 0) CALL print_error(status)
  
 RETURN
 END SUBROUTINE write_taur
 
!works only if lstore_opac, which will be default mode and option deprecated
 SUBROUTINE write_cont_opac(atom)
 ! ----------------------------------------------------------- !
 ! write separate total continuum opacities for atom
 ! ----------------------------------------------------------- !
  type(AtomType), intent(in) :: atom
  integer :: unit, status
  integer :: naxis, nelements
  integer(8) :: blocksize, naxes(4), group, bitpix, fpixel
  logical :: extend, simple
  integer :: kc, kr, i, j, Nblue, Nred, Z, icell, iend, istart
  real(kind=dp), dimension(:,:), allocatable :: chi_bf, eta_bf, sca
  real(kind=dp), dimension(:,:), allocatable ::chi_ff, eta_ff
  real(kind=dp), dimension(:), allocatable :: Bpnu, chi, eta
  character(len=100) :: filename
  
  allocate(Bpnu(NLTEspec%Nwaves),chi(NLTEspec%Nwaves),eta(NLTEspec%Nwaves))
  allocate(chi_bf(NLTEspec%Nwaves,atmos%Nspace),eta_bf(NLTEspec%Nwaves,atmos%Nspace),sca(NLTEspec%Nwaves,atmos%Nspace))

  
  !get unique unit number
  call ftgiou(unit, status)
  if (status > 0) then 
   write(*,*) " Cannot get iou"
   CALL print_error(status)
   unit = 1
  endif

  

  filename = trim(atom%ID)//'_contopac.fits.gz'

!   endif
  blocksize=1
  CALL ftinit(unit,trim(root_dir)//"/"//trim(filename),blocksize,status)
  if (status > 0) then
   write(*,*) " Cannot write contopac file"
   CALL print_error(status)
  endif

  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .true.
  bitpix = -64
  
  if (lVoronoi) then
   	naxis = 1
   	iend = 2
   	naxes(2) = atmos%Nspace
   	nelements = naxes(2)
  else
   	if (l3D) then
    	naxis = 3
    	iend = 4
    	naxes(2) = n_rad
   	 	naxes(3) = 2*nz
    	naxes(4) = n_az
    	nelements = naxes(2) * naxes(3) * naxes(4)
   	else
    	naxis = 2
    	iend = 3
    	naxes(2) = n_rad
    	naxes(3) = nz
    	nelements = naxes(2) * naxes(3)
   	end if
  end if
  
  !always first
  naxis = naxis + 1
  naxes(1) = NLTEspec%Nwaves
  nelements = nelements * NLTEspec%Nwaves
  
  

!Thomson scat opac

!init fits also
!   write them at the end as they are atom dependent

  sca(:,:) = 0.0_dp  
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (1) "
   CALL print_error(status)
  endif
  CALL ftpkys(unit, "sca_thomson", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (2) "
   CALL print_error(status)
  endif
  do icell=1,n_cells
   sca(:,icell) = atmos%ne(icell) * sigma_e
  enddo
  CALL ftpprd(unit,group,fpixel,nelements,sca,status)  
  if (status > 0) then
   write(*,*) " contopac (3) "
   CALL print_error(status)
  endif! 
!   
  NLTEspec%AtomOpac%Kc(:,:,:) = 0.0_dp
  NLTEspec%AtomOpac%jc(:,:) = 0.0_dp

  
  if (atom%ID=="H") then !only H free now (and H-)
   allocate(chi_ff(NLTEspec%Nwaves, atmos%Nspace),&
            eta_ff(NLTEspec%Nwaves, atmos%Nspace),stat=status)
   if (status > 0) call error("allocation error in write_cont_opac")
  endif
 
!H Rayleigh or He Rayleigh 
 sca = 0.0_dp
!  First scattering opacities
  if (atom%ID=="H") then
   do icell=1,n_cells
       CALL HI_Rayleigh(1, icell)
       sca(:,icell) = NLTEspec%AtomOpac%Kc(:,icell,1) / atom%n(1,icell) /sigma_e!ground state pops, sigT*sigH
       NLTEspec%AtomOpac%Kc(:,icell,1) = 0.0_dp
   enddo
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (4) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (5) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "sca_H", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (6) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,sca,status)
  if (status > 0) then
   write(*,*) " contopac (7) "
   CALL print_error(status)
  endif

  else if ((atom%ID=="He".and.associated(Helium))) then
   do icell=1,n_cells
       CALL HeI_Rayleigh(1, icell)
       sca(:,icell) = NLTEspec%AtomOpac%Kc(:,icell,1) / atom%n(1,icell) / sigma_e
       NLTEspec%AtomOpac%Kc(:,icell,1) = 0.0_dp
   enddo
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (8) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (9) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "sca_He", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (10) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,sca,status) 
  if (status > 0) then
   write(*,*) " contopac (11) "
   CALL print_error(status)
  endif 
  endif


  chi = 0.0_dp
  eta = 0.0_dp
  chi_bf = 0.0_dp
  eta_bf = 0.0_dp

!Hminus bf and ff
  
  !write Hminus bf and ff if He
  if (atom%ID=='H') then
   do icell=1,n_cells
    CALL Hminus_bf(icell, chi, eta)
    chi_bf(:,icell) = chi
    eta_bf(:,icell) = eta
   enddo
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (12) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (13) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "H- chi_bf", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (14) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,chi_bf,status)
  if (status > 0) then
   write(*,*) " contopac (15) "
   CALL print_error(status)
  endif
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (16) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (17) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "H- eta_bf", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (18) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,eta_bf,status)
  if (status > 0) then
   write(*,*) " contopac (19) "
   CALL print_error(status)
  endif
   chi_bf = 0.0_dp
   eta_bf = 0.0_dp
   chi_ff = 0.0_dp
   eta_ff = 0.0_dp
   do icell=1, n_cells
    CALL Bplanck(atmos%T(icell), Bpnu)
    CALL Hminus_ff(icell, chi)
    chi_ff(:,icell) = chi
    eta_ff(:,icell) = chi * Bpnu
   enddo


   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (20) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (21) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "H- chi_ff", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (22) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,chi_ff,status)
  if (status > 0) then
   write(*,*) " contopac (23) "
   CALL print_error(status)
  endif
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (24) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (25) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "H- eta_ff", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (26) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,eta_ff,status)
  if (status > 0) then
   write(*,*) " contopac (27) "
   CALL print_error(status)
  endif
   
   chi_ff = 0.0_dp
   eta_ff = 0.0_dp
   do icell=1, n_cells
    CALL Bplanck(atmos%T(icell), Bpnu)
    CALL Hydrogen_ff(icell, chi)
    chi_ff(:,icell) = chi
    eta_ff(:,icell) = chi * Bpnu
   enddo

   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (28) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (29) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "H chi_ff", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (30) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,chi_ff,status)
  if (status > 0) then
   write(*,*) " contopac (31) "
   CALL print_error(status)
  endif
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (32) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (33) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, "H eta_ff", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (34) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,eta_ff,status)
  if (status > 0) then
   write(*,*) " contopac (35) "
   CALL print_error(status)
  endif
   
   deallocate(chi_ff, eta_ff)
  else
   write(*,*) " metal f-f not written yet"
  endif

!now bound-free

  chi_bf = 0.0_dp
  eta_bf = 0.0_dp
  
  do icell=1, n_cells
    chi = 0.0_dp
    eta = 0.0_dp
    do kc=atom%Ntr_line+1,atom%Ntr !sum for this cell opac of all cont
     kr = atom%at(kc)%ik 
     i = atom%continua(kr)%i
     j = atom%continua(kr)%j 
     Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred

   
     chi(Nblue:Nred) = chi(Nblue:Nred) + atom%continua(kr)%alpha(:) * &
      (atom%n(i,icell)-atom%continua(kr)%gij(:,icell)*atom%n(j,icell))
      
     eta(Nblue:Nred) = eta(Nblue:Nred) + atom%continua(kr)%alpha(:) * &
      atom%continua(kr)%twohnu3_c2(:) * atom%continua(kr)%gij(:,icell) * atom%n(j,icell)


    end do ! loop over Ncont
    chi_bf(:,icell) = chi
    eta_bf(:,icell) = eta
  enddo
  
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (36) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (37) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, atom%ID//" chi_bf", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (38) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,chi_bf,status)
  if (status > 0) then
   write(*,*) " contopac (39) "
   CALL print_error(status)
  endif
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (40) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (41) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit, atom%ID//" eta_bf", "(m^-1)", ' ', status)
  if (status > 0) then
   write(*,*) " contopac (42) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,eta_bf,status)
  if (status > 0) then
   write(*,*) " contopac (43) "
   CALL print_error(status)
  endif
  

 if (allocated(chi_ff)) deallocate(chi_ff, eta_ff)
 deallocate(Bpnu,chi,eta,chi_bf, eta_bf, sca)
 
 CALL ftclos(unit, status)
  if (status > 0) then
   write(*,*) " contopac (45) "
   CALL print_error(status)
  endif
 CALL ftfiou(unit, status)
  if (status > 0) then
   write(*,*) " contopac (46) "
   CALL print_error(status)
  endif  

 RETURN
 END SUBROUTINE write_cont_opac
 
 SUBROUTINE write_atom_xsections_bf_ff(atom)
 ! ----------------------------------------------------------- !
 ! write Hydrogenic bf and ff cross-sections per atom
 ! ----------------------------------------------------------- !
  type(AtomType), intent(in) :: atom
  integer :: unit, status = 0, blocksize, naxes(5), naxis,group, bitpix, fpixel
  logical :: extend, simple
  integer :: nelements, kc, kr, i, j, Nblue, Nred, Z
  real(kind=dp), dimension(NLTEspec%Nwaves) :: alpha_ff, alpha_bf
  character(len=50) :: filename
  
  alpha_ff(:) = 0.0_dp
  alpha_bf(:) = 0.0_dp

    do kc=atom%Ntr_line+1,atom%Ntr
     kr = atom%at(kc)%ik 
     i = atom%continua(kr)%i
     j = atom%continua(kr)%j 
     Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred

   
     Z = atom%stage(j) !1 if H
 
    !per gaunt_ff which is about 1 anyway; alpha_ff (m^5 H^3 K^1/2 * ne(m^-3) in m^2 Hz^3 K^1/2 ) given per unit electron density
    !and per nu**3					so alpha_ff in m^-2 K^1/2
    alpha_ff(Nblue:Nred) = alpha_ff(Nblue:Nred) + sigma0_H_ff * real(Z*Z) * (NLTEspec%lambda(Nblue:Nred)*NM_TO_M/CLIGHT)**3
       				
    !in m2
    alpha_bf(Nblue:Nred) = alpha_bf(Nblue:Nred) + atom%continua(kr)%alpha(:)


    end do ! loop over Ncont

  !get unique unit number
  CALL ftgiou(unit,status)
  if (atom%ID=="H") then
   filename = trim(atom%ID)//"_anu.fits.gz"
  else
   filename = atom%ID//'_anu.fits.gz'
  endif
  blocksize=1
  simple = .true. !Standard fits
  group = 1
  fpixel = 1
  extend = .true.
  bitpix = -64


  naxis = 1
  naxes(1) = NLTEspec%Nwaves
  nelements = naxes(1)
  

  CALL ftinit(unit,trim(filename),blocksize,status)
  !  Initialize parameters about the FITS image
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  ! Additional optional keywords
  CALL ftpkys(unit, "alpha_bf ", "(m^2)", ' ', status)
  !write data
  CALL ftpprd(unit,group,fpixel,nelements,alpha_bf,status)
  
  CALL ftcrhd(unit, status)
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  CALL ftpkys(unit, "alpha_ff", "(m^2)", ' ', status)
  CALL ftpprd(unit,group,fpixel,nelements,alpha_ff,status)


  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  if (status > 0) CALL print_error(status)
  
 RETURN
 END SUBROUTINE write_atom_xsections_bf_ff

END MODULE write_opacity
