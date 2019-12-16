! ---------------------------------------------------------------------- !
! This module writes to separate files the contributions of each atom.
! ---------------------------------------------------------------------- !
MODULE write_opacity

 use atmos_type, only : atmos, Hydrogen, Helium
 use atom_type
 use spectrum_type
 use math, only : locate, Linear_1D, bezier2_interp, bezier3_interp
 use metal, only : metal_bb_new
 use opacity, only : NLTE_bound_bound
 use constant
 use Planck
 use hydrogen_opacities
 use rayleigh_scattering

 !MCFOST's original modules
 use stars, only : intersect_stars
 use fits_utils, only : print_error, cfitsWrite
 use grid
 use utils, only : appel_syst, interp_dp
 use messages
 use mcfost_env, only : dp

 IMPLICIT NONE

 character(len=15), parameter :: Jnu_File = "Jnu.fits.gz", Jnu_File_ascii="Jnu_no_interp.s"
 character(len=15), parameter :: taur_File = "taur.fits.gz"
 character(len=15), parameter :: contrib_file = "contrib.fits.gz"



 CONTAINS
 
 SUBROUTINE integ_cont_optical_depth(id,icell_in,xi,yi,zi,u,v,w,chic,tauc)
  ! Given chic, the continuum opacity, 
  !integrate the continuum optical depth in the direction (u,v,w).
  !the continuum opacity is an input
  integer, intent(in) :: id, icell_in
  real(kind=dp), intent(in) :: xi,yi,zi
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: chic(:,:)
  real(kind=dp), intent(out) :: tauc(:,:)
  integer :: icell, previous_cell, next_cell, icell_star, i_star, nbr_cell, icell_before
  real(kind=dp) :: x1, y1, z1, x0, y0, z0, l, ltot
  real(kind=dp) :: l_contrib, l_void_before, lmin, lmax
  logical :: lcellule_non_vide, lintersect_stars
  
  


  x1=xi;y1=yi;z1=zi

  next_cell = icell_in
  icell_before = icell_in
  nbr_cell = 0
  
  call intersect_stars(xi,yi,zi, u,v,w, lintersect_stars, i_star, icell_star)


  ! Boucle infinie sur les cellules
  do ! Boucle infinie
     ! Indice de la cellule
     icell = next_cell
     x0=x1;y0=y1;z0=z1

     if (icell <= n_cells) then
      lcellule_non_vide = (atmos%icompute_atomRT(icell)>0)
     else
      lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) then
		return
     endif
     
     nbr_cell = nbr_cell + 1
     
    if (lintersect_stars) then
      if (icell == icell_star) return
    endif

	previous_cell = 0
    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)
    
    if (lcellule_non_vide) then
    
     	tauc(:,icell) = tauc(:,icell_before) + l_contrib * chic(:,icell) * au_to_m
		icell_before = icell
     
    endif

  enddo ! boucle infinie

  write(*,*) "BUG"
  return
 
 RETURN
 END SUBROUTINE integ_cont_optical_depth
 
 SUBROUTINE integ_optical_depth_bb(id,iray,icell_in,xi,yi,zi,u,v,w,tau0,chibb,etabb)
  !Integrate total (line + continuum) optical depth in the direction (u,v,w) 
  !and keep also chibb and etabb the total opacity and total emissivity.
  !Opacities are computed on the fly.
  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: xi,yi,zi
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(out) :: tau0(:,:), chibb(:,:), etabb(:,:)
  integer :: icell, previous_cell, next_cell, icell_star, i_star, nbr_cell, icell_before, itest
  real(kind=dp) :: x1, y1, z1, x0, y0, z0, l, ltot
  real(kind=dp) :: l_contrib, l_void_before, lmin, lmax
  logical :: lcellule_non_vide, lintersect_stars

  x1=xi;y1=yi;z1=zi
  lmax=0.
  ltot=0.

  next_cell = icell_in
  icell_before = icell_in
  nbr_cell = 0
  itest = locate(NLTEspec%lambda, 654.0_dp)

  call intersect_stars(xi,yi,zi, u,v,w, lintersect_stars, i_star, icell_star)

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
     ! Indice de la cellule
     icell = next_cell
     x0=x1;y0=y1;z0=z1
     
     if (icell <= n_cells) then
      lcellule_non_vide = (atmos%icompute_atomRT(icell)>0)
     else
      lcellule_non_vide=.false.
     endif

     ! Test sortie
     if (test_exit_grid(icell, x0, y0, z0)) then
        return
     endif
     
     nbr_cell = nbr_cell + 1
          
    if (lintersect_stars) then
      if (icell == icell_star) return
    endif


    previous_cell = 0
    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

    if (lcellule_non_vide) then
    
        chibb(:,icell) = NLTEspec%AtomOpac%Kc(:,icell)
        etabb(:,icell) = NLTEspec%AtomOpac%jc(:,icell)
        if (atmos%electron_scattering) etabb(:,icell) = etabb(:,icell) + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell)
        if (atmos%NactiveAtoms>0) then 
         chibb(:,icell) =  chibb(:,icell) + NLTEspec%AtomOpac%Kc_nlte(:,icell)
         etabb(:,icell) = etabb(:,icell) + NLTEspec%AtomOpac%jc_nlte(:,icell)
        endif
        
		if (atmos%NpassiveAtoms>0) then
    		call initAtomOpac(id)
    		call Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
    		chibb(:,icell) = chibb(:,icell) + NLTEspec%AtomOpac%chi_p(:,id)
    		etabb(:,icell) = etabb(:,icell) + NLTEspec%AtomOpac%eta_p(:,id)
    	endif
 		if (atmos%NactiveAtoms>0) then 
			call initAtomOpac_nlte(id)
       		call NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)
    		chibb(:,icell) = chibb(:,icell) + NLTEspec%AtomOpac%chi(:,id)
    		etabb(:,icell) = etabb(:,icell) + NLTEspec%AtomOpac%eta(:,id)
		endif
		
     	tau0(:,icell) = tau0(:,icell_before) + l_contrib * chibb(:,icell) * AU_to_m
     	lmax = lmax + l_contrib
     	ltot = ltot + l_contrib * chibb(itest,icell) * AU_to_m
	    !write(*,*) "icell=",icell, "lentgh=",lmax, "tau=",tau0(itest,icell), "taubis=", ltot
		icell_before = icell

    endif

  enddo infinie

  write(*,*) "BUG"


 RETURN
 END SUBROUTINE integ_optical_depth_bb
 
!building
 SUBROUTINE write_contrib_lambda_ascii(lam,u,v,w)
 ! ----------------------------------------------------- !
  ! - tau(lam)
  ! - chi_c (lam)
  ! - sigma_scatt (lam)
  ! - chi_l (lam)
  ! - eta_l (lam)
  ! - eta_c (lam)
  ! - eta_scatt (lam)
  ! - tauR (lam)
  ! - chiR (lam)
  !
  ! TO DO: loop over wavelengths to write contrib at
  ! specific wavelength for each.
 ! ---------------------------------------------------- !
  real(kind=dp), intent(in) :: lam
  real(kind=dp), intent(in) :: u,v,w
  integer :: ilam, ila, i, j, kr, kc
  integer :: icell, icell0, id, iray
  real(kind=dp) :: u0, v0, w0
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, lmin, lmax, l_contrib, l_void_before
  logical :: lintersect
  real(kind=dp), dimension(:), allocatable :: taul, chiR, chic, sigmac, etac, etasc
  real(kind=dp), dimension(:), allocatable :: tauc, tauR, chil, etal
  real(kind=dp), dimension(:,:), allocatable :: chibb, etabb, tau
  character(len=50) :: filename
  
  write(filename, "(1I)") nint(lam)
  filename = "contrib_at"//trim(adjustl(filename))//"_nm.s"
  
  ilam = locate(NLTEspec%lambda, lam)
  write(*,*) " Write opacities at ", NLTEspec%lambda(ilam), lam, "nm...", filename
  id = 1
  iray = 1 !not used if eval_operator = .false.
  !centre of each cell
  
  x0 = 10*Rmax*u; y0 = 10*Rmax*v; z0 = 10*Rmax*w
  u0 = -u; v0=-v; w0=-w

  call move_to_grid(id,x0,y0,z0,u0,v0,w0,icell0,lintersect)




  if (.not.lintersect) then
   write(*,*) "Cannot write taur in this direction there is nothing.."
   return
  endif
  
  allocate(taul(atmos%Nspace), chil(atmos%Nspace), etal(atmos%Nspace))
  taul = 0.0_dp; chil = 0.0_dp; etal = 0.0_dp
  allocate(chic(atmos%Nspace), etac(atmos%Nspace))
  chic = 0.0_dp; etac = 0.0_dp
  allocate(sigmac(atmos%Nspace), etasc(atmos%Nspace))
  sigmac = 0.0_dp; etasc = 0.0_dp

  do icell=1, n_cells !should be zero if icompute_atomRT < 1
  	chic(icell) = NLTEspec%AtomOpac%Kc(ilam,icell)!interp_dp(NLTEspec%AtomOpac%Kc(:,icell), NLTEspec%lambda, lam)
	sigmac(icell) = NLTEspec%AtomOpac%sca_c(ilam,icell)!interp_dp(NLTEspec%AtomOpac%sca_c(:,icell), NLTEspec%lambda, lam)
  	etac(icell) = NLTEspec%AtomOpac%jc(ilam,icell)!interp_dp(NLTEspec%AtomOpac%jc(:,icell),NLTEspec%lambda, lam)
  	
   if (atmos%electron_scattering) then 
  	etasc(icell) = NLTEspec%AtomOpac%sca_c(ilam,icell) * NLTEspec%Jc(ilam,icell)!interp_dp(NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell), NLTEspec%lambda, lam)
  	etac(icell) = etac(icell) + etasc(icell)
   endif
   if (atmos%NactiveAtoms>0) then
  	etac(icell) = etac(icell) + NLTEspec%AtomOpac%jc_nlte(ilam,icell)!interp_dp(NLTEspec%AtomOpac%jc_nlte(:,icell), NLTEspec%lambda, lam)
  	chic(icell) = chic(icell) + NLTEspec%AtomOpac%Kc_nlte(ilam,icell)!interp_dp(NLTEspec%AtomOpac%Kc_nlte(:,icell), NLTEspec%lambda, lam)
   endif
  		
  enddo !over cells
  
  if (lintersect) then
    allocate(chibb(NLTEspec%Nwaves, n_cells), etabb(NLTEspec%NWaves, n_cells), tau(NLTEspec%NWaves, n_cells))
    chibb = 0.; etabb = 0.; tau=0.
    call integ_optical_depth_bb(id,iray,icell0,x0,y0,z0,u0,v0,w0,tau,chibb, etabb)
    write(*,*) "Tau(max/min)", maxval(tau(ilam,:)), minval(tau(ilam,:))

    do icell=1, n_cells
     taul(icell) = tau(ilam,icell)!interp_dp(tau(:,icell), NLTEspec%lambda, lam)
     chil(icell) = chibb(ilam,icell)!interp_dp(chibb(:,icell), NLTEspec%lambda, lam)
     etal(icell) = etabb(ilam,icell)!interp_dp(etabb(:,icell), NLTEspec%lambda, lam)
     !write(*,*) icell," taul=", taul(icell), " etal=", etal(icell), " chil=", chil(icell)
    enddo
    deallocate(chibb, etabb, tau)
  endif
  
  !write results
  open(unit=1, file=filename, status="unknown")
  write(1,*) lam, n_cells
  write(1,*) " header"
  
  do icell=1, n_cells !how to read integer then float ?
   write(1,"(8E)") real(icell), taul(icell), chil(icell), etal(icell), etac(icell), chic(icell), etasc(icell), sigmac(icell)
  enddo
  
  close(1)
    
  deallocate(taul, chil, etal)
  deallocate(chic, etac)
  deallocate(sigmac, etasc)
 
 RETURN
 END SUBROUTINE write_contrib_lambda_ascii
 
!building
 SUBROUTINE write_taur(lambda,u,v,w)
 ! ----------------------------------------------------------- !
 ! write the continuum optical depth at reference
 ! wavelength in a specific direction
 ! the wavelength-point is found by locate()
 ! ----------------------------------------------------------- !
  real(kind=dp), intent(in) :: lambda
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp) ::  lmin, lmax, x0, y0, z0, u0, w0, v0
  integer :: iref, icell, id, icell0
  integer :: unit, status = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
  logical :: extend, simple, lintersect
  integer :: nelements
  real(kind=dp), dimension(:), allocatable   :: tauc
  real(kind=dp), dimension(:,:), allocatable :: chic, taur
  
  id = 1
  iref = locate(NLTEspec%lambda, lambda)
  write(*,*) "   Writing optical depth at ", NLTEspec%lambda(iref), lambda, " nm..."
  
  !presently I only compute it for one lambda, the the integration subroutine is more general
  !and accept chic(:,:)
  allocate(taur(1,n_cells), chic(1,n_cells))
  taur = 0.0
  
  do icell=1, n_cells !should be zero if icompute_atomRT < 1
  
  	chic(1,icell) = NLTEspec%AtomOpac%Kc(iref,icell)
  	
    if (atmos%NactiveAtoms>0) then
  	 chic(1,icell) = chic(1,icell) + NLTEspec%AtomOpac%Kc_nlte(iref,icell)
    endif
  		
  enddo !over cells

  x0 = 10*Rmax*u; y0 = 10*Rmax*v; z0 = 10*Rmax*w
  u0 = -u; v0=-v; w0=-w
  call move_to_grid(id,x0,y0,z0,u0,v0,w0,icell0,lintersect)


  if (lintersect) then
   !Actually, the continuum optical depth at all frequencies could be computed
   call integ_cont_optical_depth(id,icell0,x0,y0,z0,u0,v0,w0,chic,taur)
   write(*,*) " Maxval/minval of taur", maxval(taur), minval(taur)
  else
   write(*,*) " Writing taur, there is nothing in this direction.."
  endif
  
  deallocate(chic)
  allocate(tauc(n_cells))
  tauc(:) = taur(1,:)
  deallocate(taur)
  
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
  CALL ftpprd(unit,group,fpixel,nelements,tauc,status)
  

  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  if (status > 0) CALL print_error(status)
  
  deallocate(tauc)
  
  write(*,*) "   ..done"
  
 RETURN
 END SUBROUTINE write_taur
 
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
  
  write(*,*) " Reading Jnu from Jnu_file_ascii..."

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
    CALL bezier3_interp(Nlam_j, xtmp, Jtmp, NLTEspec%Nwaves, NLTEspec%lambda, NLTEspec%Jc(:,icell))
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
 
!building
 SUBROUTINE write_cont_opac_ascii(atom)
 ! ----------------------------------------------------------- !
 ! write separate total continuum opacities for atom
 ! ----------------------------------------------------------- !
  type(AtomType), intent(in) :: atom
  integer :: unit, status
  integer :: naxis, nelements
  integer(8) :: blocksize, naxes(10), group, bitpix, fpixel
  logical :: extend, simple
  integer :: kc, kr, i, j, Nblue, Nred, Z, icell, iend, istart, la, Nwrite
  real(kind=dp), dimension(:), allocatable :: Bpnu, chi, eta
  real(kind=dp), dimension(:,:,:), allocatable :: big_tab
  character(len=100) :: filename
  
  allocate(Bpnu(NLTEspec%Nwaves),chi(NLTEspec%Nwaves),eta(NLTEspec%Nwaves))
  
  !get unique unit number
!   call ftgiou(unit, status)
!   if (status > 0) then 
!    write(*,*) " Cannot get iou"
!    CALL print_error(status)
!    unit = 1
!   endif

  unit = 1
  

  filename = trim(atom%ID)//'_contopac.s'
  
  !by default 4
  Nwrite = 4 !b-b and b-f
  if (atom%ID=="H") then
                 !thomson + rayleigh + H_ff, H-ff H-bf
   Nwrite = Nwrite + 2 + 2 + 2 + 2
  else if (atom%ID=="He") then
   Nwrite = Nwrite + 1 !rayleigh HE 
  endif
  allocate(big_tab(Nwrite, NLTEspec%Nwaves, n_cells))
  big_tab = 0.0
  !free free for other atom not included yet

  open(unit, file=filename, status="unknown")
  
  write(unit,*) n_cells, NLTEspec%Nwaves, Nwrite
  
  !write wavelengths
  do la=1, NLTEspec%Nwaves
   write(unit,*) NLTEspec%lambda(la)
  enddo


  NLTEspec%AtomOpac%sca_c(:,:) = 0.0_dp ! set it to zero for routines called after this one
  
 istart = 0
!H Rayleigh or He Rayleigh 
!  First scattering opacities only for H and He presently
  if (atom%ID=="H") then
   NLTEspec%AtomOpac%sca_c(:,:) = 0.0_dp
   do icell=1,n_cells
       !CALL HI_Rayleigh(1, icell)
       CALL Rayleigh(1, icell, hydrogen)
       big_tab(1,:,icell) = atmos%ne(icell) * sigma_e
       big_tab(2,:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) / atom%n(1,icell) /sigma_e!ground state pops, sigT*sigH
   enddo

  else if ((atom%ID=="He".and.associated(Helium))) then
   istart = 1
   NLTEspec%AtomOpac%sca_c(:,:) = 0.0_dp !set to zero, we only want He
   do icell=1,n_cells
       CALL HeI_Rayleigh(1, icell) !for helium starts at 1
       big_tab(1,:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) / atom%n(1,icell) /sigma_e!ground state pops, sigT*sigH
   enddo
 endif


  chi = 0.0_dp
  eta = 0.0_dp

!Hminus bf and ff
  
  !write Hminus bf and ff if He
  if (atom%ID=='H') then
   do icell=1,n_cells
    CALL Hminus_bf(icell, chi, eta)
    big_tab(3,:,icell) = chi
    big_tab(4,:,icell) = eta
    
    chi = 0.0
    eta = 0.
    CALL Bplanck(atmos%T(icell), Bpnu)
    CALL Hminus_ff(icell, chi)
    big_tab(5,:,icell) = chi
    big_tab(6,:,icell) = chi * Bpnu
   
    chi = 0.0_dp
    eta = 0.0_dp
!     CALL Bplanck(atmos%T(icell), Bpnu)
    CALL Hydrogen_ff(icell, chi)
    big_tab(7,:,icell) = chi
    big_tab(8,:,icell) = chi * Bpnu
   enddo
   
   deallocate(Bpnu)
   istart = 8 + istart ! + 0 if H
  else
   write(*,*) " metal f-f not written yet"
   istart = 1 + istart ! = 2 if He because of Rayleigh
  endif
  

!now bound-free
  
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
    big_tab(istart+1,:,icell) = chi
    big_tab(istart+2,:,icell) = eta
  enddo
 
 !now bound bound in the zero-velocity
 !Assume that profile are already computed
 if (.not.allocated(atom%lines(1)%phi)) then
  call error("Cannot write bound bound, phi not allocated")
 endif 
  
  do icell=1, n_cells
    chi = 0.0_dp
    eta = 0.0_dp
    do kc=1,atom%Ntr_line
     kr = atom%at(kc)%ik 
     i = atom%lines(kr)%i
     j = atom%lines(kr)%j 
     Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred

   
     chi(Nblue:Nred) = chi(Nblue:Nred) + hc_fourpi * (atom%lines(kr)%Bij * atom%n(i,icell) - &
     						atom%lines(kr)%Bji * atom%n(j,icell)) * atom%lines(kr)%phi(:,icell)
     						
     eta(Nblue:Nred) = eta(Nblue:Nred) + hc_fourpi * atom%lines(kr)%Aji * atom%n(j,icell) * atom%lines(kr)%phi(:,icell)


    end do ! loop over Ncont
    big_tab(istart+3,:,icell) = chi
    big_tab(istart+4,:,icell) = eta
  enddo
  
  !assumes no free free other than for H, and no rayleigh other than H, He
  if (atom%ID=='H') then
   write(unit, *) "header"
  else if (atom%ID=="He") then
   write(unit, *) "header"
  else
   write(unit, *) "header"
  endif
  
  !can be better with format line
  do icell=1, n_cells
   write(unit, *) icell
   do la=1, NLTEspec%Nwaves
  
    if (Nwrite == 12) then
    write(unit, "(12E)") (big_tab(i,la,icell), i=1, 12)
    
    else if (Nwrite==5) then
    write(unit, "(5E)") (big_tab(i,la,icell), i=1, 5)
    else
    write(unit, "(4E)") (big_tab(i,la,icell), i=1, 4)
    endif
   enddo
  enddo
 
 deallocate(chi, eta, big_tab)

 
  close(unit)
!   CALL ftfiou(unit, status)
!   if (status > 0) then
!    write(*,*) " contopac (53) "
!    CALL print_error(status)
!   endif  

 RETURN
 END SUBROUTINE write_cont_opac_ascii

!building
!works only if lstore_opac, which will be default mode and option deprecated
 SUBROUTINE write_cont_opac(atom)
 ! ----------------------------------------------------------- !
 ! write separate total continuum opacities for atom
 ! ----------------------------------------------------------- !
  type(AtomType), intent(in) :: atom
  integer :: unit, status
  integer :: naxis, nelements
  integer(8) :: blocksize, naxes(10), group, bitpix, fpixel
  logical :: extend, simple
  integer :: kc, kr, i, j, Nblue, Nred, Z, icell, iend, istart
  real(kind=dp), dimension(:,:), allocatable :: chi_bf, eta_bf, sca
  real(kind=dp), dimension(:,:), allocatable ::chi_ff, eta_ff,chi_bb, eta_bb
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
  
  !write(*,*) "naxis contopac=", naxis, nelements
  !write(*,*) "naxes=", naxes(1:naxis)

!Thomson scat opac

!init fits also
!   write them at the end as they are atom dependent

  sca(:,:) = 0.0_dp  
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (1) "
   CALL print_error(status)
  endif
  !CALL ftpkys(unit, "sca_thomson", "(m^-1)", '', status)
  CALL ftpkys(unit,'BUNIT',"m-1",'Thomson',status)
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
  
!-> set to zero for routines called after this one
!   
!   NLTEspec%AtomOpac%Kc(:,:) = 0.0_dp
!   NLTEspec%AtomOpac%jc(:,:) = 0.0_dp
  NLTEspec%AtomOpac%sca_c(:,:) = 0.0_dp

  
 
!H Rayleigh or He Rayleigh 
 sca = 0.0_dp
!  First scattering opacities
  if (atom%ID=="H") then
   do icell=1,n_cells
       CALL HI_Rayleigh(1, icell)
       sca(:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) / atom%n(1,icell) /sigma_e!ground state pops, sigT*sigH
       NLTEspec%AtomOpac%sca_c(:,icell) = 0.0_dp
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
  CALL ftpkys(unit,'BUNIT',"m-1",'Rayleigh H',status)
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
       sca(:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) / atom%n(1,icell) / sigma_e
       NLTEspec%AtomOpac%sca_c(:,icell) = 0.0_dp
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
  CALL ftpkys(unit,'BUNIT',"m-1",'Rayleigh He',status)
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
  
  deallocate(sca)


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
   !CALL ftpkys(unit, "H- chi_bf", "(m^-1)", '', status)
   CALL ftpkys(unit,'BUNIT',"m-1",'H- chi_bf',status)
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
   CALL ftpkys(unit,'BUNIT',"W/m3/Hz/sr",'H- eta_bf',status)
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
   
   !if (atom%ID=="H") then !only H free now (and H-)
   allocate(chi_ff(NLTEspec%Nwaves, atmos%Nspace),&
            eta_ff(NLTEspec%Nwaves, atmos%Nspace),stat=status)
   if (status > 0) call error("allocation error in write_cont_opac")
  !endif  
   
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
   CALL ftpkys(unit,'BUNIT',"m-1",'H- chi_ff',status)
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
   CALL ftpkys(unit,'BUNIT',"W/m3/Hz/sr",'H- eta_ff',status)
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
   
   deallocate(Bpnu)


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
   CALL ftpkys(unit,'BUNIT',"m-1",'H chi_ff',status)
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
   CALL ftpkys(unit,'BUNIT',"W/m3/Hz/sr",'H eta_ff',status)
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
  
  if (allocated(chi_ff)) deallocate(chi_ff, eta_ff)

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
   CALL ftpkys(unit,'BUNIT',"m-1",atom%ID//' eta_ff',status)
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
   CALL ftpkys(unit,'BUNIT',"W/m3/Hz/sr",atom%ID//' eta_bf',status)
  if (status > 0) then
   write(*,*) " contopac (42) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,eta_bf,status)
  if (status > 0) then
   write(*,*) " contopac (43) "
   CALL print_error(status)
  endif
  
 deallocate(chi_bf, eta_bf)
 
 !now bound bound in the zero-velocity
 !Assume that profile are already computed
 if (.not.allocated(atom%lines(1)%phi)) then
  call error("Cannot write bound bound, phi not allocated")
 endif
 allocate(chi_bb(NLTEspec%Nwaves,atmos%Nspace),eta_bb(NLTEspec%Nwaves,atmos%Nspace))
 
  chi_bb = 0.0_dp
  eta_bb = 0.0_dp
  
  do icell=1, n_cells
    chi = 0.0_dp
    eta = 0.0_dp
    do kc=1,atom%Ntr_line
     kr = atom%at(kc)%ik 
     i = atom%lines(kr)%i
     j = atom%lines(kr)%j 
     Nblue = atom%lines(kr)%Nblue; Nred = atom%lines(kr)%Nred

   
     chi(Nblue:Nred) = chi(Nblue:Nred) + hc_fourpi * (atom%lines(kr)%Bij * atom%n(i,icell) - &
     						atom%lines(kr)%Bji * atom%n(j,icell)) * atom%lines(kr)%phi(:,icell)
     						
     eta(Nblue:Nred) = eta(Nblue:Nred) + hc_fourpi * atom%lines(kr)%Aji * atom%n(j,icell) * atom%lines(kr)%phi(:,icell)


    end do ! loop over Ncont
    chi_bb(:,icell) = chi
    eta_bb(:,icell) = eta
  enddo
  
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (44) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (45) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit,'BUNIT',"m-1",atom%ID//' chi_bb',status)
  if (status > 0) then
   write(*,*) " contopac (46) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,chi_bb,status)
  if (status > 0) then
   write(*,*) " contopac (47) "
   CALL print_error(status)
  endif
   CALL ftcrhd(unit, status)
  if (status > 0) then
   write(*,*) " contopac (48) "
   CALL print_error(status)
  endif
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  if (status > 0) then
   write(*,*) " contopac (49) "
   CALL print_error(status)
  endif
   CALL ftpkys(unit,'BUNIT',"W/m3/Hz/sr",atom%ID//' eta_bb',status)
  if (status > 0) then
   write(*,*) " contopac (50) "
   CALL print_error(status)
  endif
   CALL ftpprd(unit,group,fpixel,nelements,eta_bb,status)
  if (status > 0) then
   write(*,*) " contopac (51) "
   CALL print_error(status)
  endif
 
 deallocate(chi, eta)
 deallocate(chi_bb, eta_bb)

 
 CALL ftclos(unit, status)
  if (status > 0) then
   write(*,*) " contopac (52) "
   CALL print_error(status)
  endif
 CALL ftfiou(unit, status)
  if (status > 0) then
   write(*,*) " contopac (53) "
   CALL print_error(status)
  endif  

 RETURN
 END SUBROUTINE write_cont_opac
 
 SUBROUTINE write_atom_xsections_bf_ff(atom)
 ! ----------------------------------------------------------- !
 ! write Hydrogenic bf and ff cross-sections per atom
 ! the ff cross-sections is written for all cells
 ! because the cross section is a function of ne and T (and the
 ! absorption of ne, T and nion)
 ! ----------------------------------------------------------- !
  type(AtomType), intent(in) :: atom
  integer :: unit, status = 0, blocksize, naxes(5), naxis,group, bitpix, fpixel, icell
  logical :: extend, simple
  integer :: nelements, kc, kr, i, j, Nblue, Nred, Z, la
  real(kind=dp), allocatable :: alpha_ff(:,:), alpha_bf(:)
  character(len=50) :: filename
  
  allocate(alpha_ff(NLTEspec%Nwaves, n_cells))
  allocate(alpha_bf(NLTEspec%Nwaves))
  alpha_ff(:,:) = 0.0_dp
  alpha_bf(:) = 0.0_dp

    do kc=atom%Ntr_line+1,atom%Ntr
     kr = atom%at(kc)%ik 
     i = atom%continua(kr)%i
     j = atom%continua(kr)%j 
     Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred

   
     Z = atom%stage(j) !1 if H
       				
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
  
  if (lVoronoi) then
   	naxis = 2
   	naxes(1) = NLTEspec%Nwaves
   	naxes(2) = atmos%Nspace
   	nelements = naxes(2)*naxes(1)
  else
   	if (l3D) then
    	naxis = 4
    	naxes(1) =  NLTEspec%Nwaves
    	naxes(2) = n_rad
   	 	naxes(3) = 2*nz
    	naxes(4) = n_az
    	nelements = naxes(2) * naxes(3) * naxes(4) * naxes(1)
   	else
    	naxis = 3
    	naxes(1) =  NLTEspec%Nwaves
    	naxes(2) = n_rad
    	naxes(3) = nz
    	nelements = naxes(2) * naxes(3) * naxes(1)
   	end if
  end if
  
  do icell=1, n_cells
    do kc=atom%Ntr_line+1,atom%Ntr
     kr = atom%at(kc)%ik 
     i = atom%continua(kr)%i
     j = atom%continua(kr)%j 
     Nblue = atom%continua(kr)%Nblue; Nred = atom%continua(kr)%Nred

   
     Z = atom%stage(j) !1 if H
       				
    !in m2
            !per gaunt_ff which is about 1 anyway; alpha_ff (m^5 H^3 K^1/2 * ne(m^-3) in m^2 Hz^3 K^1/2 ) given per unit electron density
    !and per nu**3					so alpha_ff in m^-2 K^1/2
     do la=1, atom%continua(kr)%Nlambda
      alpha_ff(Nblue+la-1,icell) = alpha_ff(Nblue+la-1,icell) + atmos%ne(icell) * H_ff_Xsection(Z, atmos%T(icell), NLTEspec%lambda(Nblue+la-1))
	 enddo



    end do ! loop over Ncont
 enddo
  
  CALL ftcrhd(unit, status)
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  CALL ftpkys(unit, "alpha_ff", "(m^2)", ' ', status)
  CALL ftpprd(unit,group,fpixel,nelements,alpha_ff,status)


  deallocate(alpha_bf, alpha_ff)
  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  if (status > 0) CALL print_error(status)
  
 RETURN
 END SUBROUTINE write_atom_xsections_bf_ff

END MODULE write_opacity
