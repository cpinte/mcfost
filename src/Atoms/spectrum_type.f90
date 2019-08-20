MODULE spectrum_type

  use atom_type
  use atmos_type, only : GridType, atmos
  use getlambda, only  : make_wavelength_grid, adjust_wavelength_grid, Read_wavelengths_table!,&
  ! 						Nred_array, Nmid_array, Nblue_array
  
  !MCFOST's original modules
  use fits_utils, only : print_error
  use parametres
  use input

  IMPLICIT NONE

   real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
   real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352
   character(len=*), parameter :: WAVES_FILE="nLTE_waves_grid.fits.gz"
   ! not that FLUX_FILE is 1D only if one pixel, otherwise it is a map.
   ! F(x,y,0,lambda) in several directions.
   character(len=*), parameter :: FLUX_FILE="flux.fits.gz"
   ! I(x,y,0,lambda,imu), not used yet
   !character(len=*), parameter :: SPEC_FILE="spectrum.fits.gz"
   character(len=*), parameter :: OPAC_CONTRIB_FILE="opacities.fits.gz"
   logical :: lcompute_continuum !if false, do not allocate cont only opacities (save memory)
   								 !if lcontrib_function, set it to true, before computing contrib function
  
  !! Store S, Sc, Slte, Sclte, and jc and Kc, to save memory ???
  TYPE AtomicOpacity
   !active opacities
   double precision, allocatable, dimension(:,:)   :: chi, eta
   double precision, dimension(:,:), allocatable :: chic_nlte, etac_nlte
   ! NLTE magneto-optical elements and dichroism are stored in the background _p arrays.
   ! Mainly because we do not use them in SEE, even if etaQUV can be added to the total
   ! emissivity. But in case, etaQUV has to be atom dependent, so now we can store the LTE
   ! and NLTE is the same variable
   !passive opacities
   double precision, allocatable, dimension(:,:)   :: eta_p, chi_p
   double precision, allocatable, dimension(:,:)   :: eta_c, chi_c, sca_c
   !separate polarised opacities. Or regroud in one array with Nsize, or pointers ?
   double precision, allocatable, dimension(:,:,:)   :: rho_p, chiQUV_p, etaQUV_p
   double precision, allocatable, dimension(:,:)   :: jc
   double precision, allocatable, dimension(:,:,:) :: Kc
   !!logical, dimension(:), allocatable :: initialized
   									     !set to .true. for each cell, when iray=1.
   									     !.false. otherwise.
   									     !%initialized(id) = (iray == 1)
  END TYPE AtomicOpacity

  TYPE Spectrum
   !n_proc should be the last index
   type  (GridType), pointer :: atmos
   logical :: vacuum_to_air=.false., updateJ, write_wavelength_grid=.false.
   integer :: Nwaves, Nact, Npass, Ntrans, NPROC=1
   double precision :: wavelength_ref=0.d0 !nm optionnal
   double precision, allocatable, dimension(:) :: lambda
   double precision, allocatable, dimension(:) :: Istar !not polarised
   !nlambda, nrays, nproc
   double precision, allocatable, dimension(:,:,:) :: I, StokesQ, StokesU, StokesV, Ic
   !nlambda, nproc
   double precision, allocatable, dimension(:,:) :: J, Jc, J20
   !Nlambda, xpix, ypix, Nincl, Naz
   double precision, allocatable, dimension(:,:,:,:,:) :: Flux, Fluxc
   double precision, allocatable, dimension(:,:,:,:,:,:) :: F_QUV
   !!double precision, allocatable, dimension(:,:) :: S_QUV
   !Contribution function
   ! N_cells, Nlambda ...
   double precision, allocatable, dimension(:,:,:,:) :: Ksi 
   real(kind=dp), allocatable, dimension(:,:,:) :: scale_ref
   ! Flux is a map of Nlambda, xpix, ypix, nincl, nazimuth
   double precision, allocatable, dimension(:,:,:) :: Psi, dtau !for cell icell in direction iray, thread id
   !size of Psi could change during the devlopment
   type (AtomicOpacity) :: AtomOpac
   character:: Jfile, J20file
  END TYPE Spectrum

  type (Spectrum) :: NLTEspec

  CONTAINS

  SUBROUTINE initSpectrum(lam0, vacuum_to_air)
   ! ------------------------------------------- !
    ! Allocate and store the wavelength grid for
    ! NLTE atomic line transfer.
   ! ------------------------------------------- !

   !integer, intent(in) :: NPROC
   integer :: kr, nat
   double precision, optional :: lam0
   logical, optional :: vacuum_to_air
   logical :: alloc_nlte_vars, alloc_image_vars
    
   NLTEspec%atmos => atmos
   alloc_nlte_vars = NLTEspec%atmos%Nactiveatoms > 0
   alloc_image_vars = .not.ltab_wavelength_image

   
   NLTEspec%Nact = NLTEspec%atmos%Nactiveatoms
   NLTEspec%Npass = NLTEspec%atmos%Npassiveatoms
   NLTEspec%NPROC = nb_proc!NPROC
   
   if (present(lam0)) NLTEspec%wavelength_ref = lam0
   if (present(vacuum_to_air)) NLTEspec%vacuum_to_air = vacuum_to_air
  
   ! Initialize the wavelength grid depending on the number of transitions PASSIVE/ACTIVE
   ! and allocate NLTEspec%lambda
   ! Set also Nblue and Nred for each transitions: extension of the transtion
   ! in absolute indexes on the whole grid.
   CALL make_wavelength_grid(NLTEspec%atmos%Natom, NLTEspec%atmos%Atoms, & 
                        NLTEspec%lambda, NLTEspec%Ntrans, NLTEspec%wavelength_ref)
   NLTEspec%Nwaves = size(NLTEspec%lambda)
   
   !determine maximum wavelength per transition
!    NLTEspec%Nwaves_trans = 0
!    do nat = 1, NLTEspec%atmos%Natom
!     do kr=1, NLTEspec%atmos%Atoms(nat)%ptr_atom%Ncont
!      NLTEspec%Nwaves_trans = max(NLTEspec%Nwaves_trans, &
!      	NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%Nlambda)
!     end do
!     do kr=1, NLTEspec%atmos%Atoms(nat)%ptr_atom%Nline
!      NLTEspec%Nwaves_trans = max(NLTEspec%Nwaves_trans, &
!      	NLTEspec%atmos%Atoms(nat)%ptr_atom%lines(kr)%Nlambda)    
!     end do
!    end do
!    write(*,*)"  Maximum number of wavelengths per transition is", NLTEspec%Nwaves_trans
   
   CALL writeWavelength()
   CALL allocSpectrum(alloc_nlte_vars,& !.true. => alloc atom%eta, %chi, %Xcoupling if NLTE
  							 !.false. => if NLTE, skip this allocation
  							 !we probably do an image and we do not need these values.
       alloc_image_vars) !alloc flux only if we do not read an ascii waves list.
   !should consider to do that only if Nlte is on, otherwise we go directly to image
  RETURN
  END SUBROUTINE initSpectrum
  
  !building
  SUBROUTINE initSpectrumImage()
  ! -------------------------------------------------------------------- !
   ! Allocate a special wavelength grid for emission_line map.
   ! This allows to solve for the RT equation only for a specific line.
   ! This fasten LTE calculation in case of 1 line.
  ! -------------------------------------------------------------------- !
   double precision, dimension(NLTEspec%Nwaves) :: old_grid
   integer, dimension(:), allocatable :: Nlam_R
   integer :: nat, kr, Ntrans_new, kp
   
   old_grid = NLTEspec%lambda
   write(*,*) " -> Redefining a wavelength grid for image.."
   CALL freeSpectrum() !first free waves arrays
   NLTEspec%atmos => atmos

   !Reset the differents Nblue, Nred for atomic transitions
   !and recompute photoionisation cross-section
   CALL Read_wavelengths_table(NLTEspec%lambda, Nlam_R)
   NLTEspec%Nwaves = size(NLTEspec%lambda)
   !Works for Active and Passive atoms alike.
   !-> sort also the grid
   !-> make sure to search for transitions in individual regions
   ! and not in the total grid ! Otherwise, the transitions are kept if they fall
   ! between lambda(min) and lambda(max)
   CALL adjust_wavelength_grid(old_grid, NLTEspec%lambda, Nlam_R, NLTEspec%atmos%Atoms)
   deallocate(Nlam_R) !not used anymore
   
   !Reallocate J and Jc keeping their value if coherent_scattering
   !NB: the spatial grid is the same
   !if (NLTEspec%atmos%coherent_scattering) then
   ! Jold = J; Joldc = Jc
   ! dealloc(J, Jc)
   ! allocate(J, Jc)
   ! J = PACK(Jold, mask=NLTEspec%mlamba==old_grid) !may not work directly
   !end if
   
   Ntrans_new = 0
   write(*,*) " -> Using ", NLTEspec%Nwaves," wavelengths for image and spectrum."
!I do not removed transitions for the moment
!so I count only transitions falling in the new wavelength. Thi summation is implicit
!if transitions are removed.
   do nat=1,NLTEspec%atmos%Natom
     write(*,*) "  --> Atom ", NLTEspec%atmos%Atoms(nat)%ptr_atom%ID!, &
!      NLTEspec%atmos%Atoms(nat)%ptr_atom%Nline, "b-b", &
!      NLTEspec%atmos%Atoms(nat)%ptr_atom%Ncont, "b-f"
     kp = 0
     do kr=1,NLTEspec%atmos%Atoms(nat)%ptr_atom%Nline
      if (NLTEspec%atmos%Atoms(nat)%ptr_atom%lines(kr)%Nblue >= 1 .and. &
      		NLTEspec%atmos%Atoms(nat)%ptr_atom%lines(kr)%Nred <= NLTEspec%Nwaves) then
       kp = kp + 1
       write(*,*) "    b-b #", kp, NLTEspec%atmos%Atoms(nat)%ptr_atom%lines(kr)%lambda0, &
        "nm"
       Ntrans_new = Ntrans_new + 1    
      endif 
     end do
    write(*,*) "   --> ", kp, "b-b transitions"
     kp = 0
     do kr=1,NLTEspec%atmos%Atoms(nat)%ptr_atom%Ncont
      if (NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%Nblue >= 1 .and. &
      		NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%Nred <= NLTEspec%Nwaves) then
       kp = kp + 1
       write(*,*) "    b-f #", kp, NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%lambdamin, &
        "nm", NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%lambda0, "nm"  
       Ntrans_new = Ntrans_new + 1  
      endif
     end do
    write(*,*) "   --> ", kp, "b-f transitions"
   end do
   write(*,*) "  ** Total number of transitions for image:", Ntrans_new
   NLTEspec%Ntrans = Ntrans_new
   !reallocate wavelength arrays, except polarisation ? which are in adjustStokesMode
   CALL allocSpectrum(.false.,.true.) !do not realloc NLTE atom%chi;%eta;%Xcoupling
   							   !even if NLTEpops are present
   							   !NLTE eval_opartor condition never reached for images.
   !!CALL writeWavelength() !not useful because they are written with the flux and
   						    ! here the wavelength comes from a file..
  RETURN
  END SUBROUTINE initSpectrumImage
  
  SUBROUTINE reallocate_rays_arrays(newNray)
   integer, intent(in) :: newNray
  
   NLTEspec%atmos%Nrays = newNray
   
   if ((NLTEspec%atmos%Nrays) /= newNray .or. (newNray /=1 )) &
   	 CALL Warning("  Beware, check the number of rays for the ray-traced map!")
   
   deallocate(NLTEspec%I, NLTEspec%Ic)
   !except polarization which are (de)allocated in adjustStokesMode
   if (NLTEspec%Nact > 0 .and.allocated(NLTEspec%PSI)) then 
    deallocate(NLTEspec%Psi)
    if (allocated(NLTEspec%dtau)) deallocate(NLTEspec%dtau)
   	allocate(NLTEspec%Psi(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
   	allocate(NLTEspec%dtau(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))  
   end if 
    
   !Could be also LTE opac if line are kept in memory ?
   allocate(NLTEspec%I(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
   allocate(NLTEspec%Ic(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
   NLTEspec%I = 0d0; NLTEspec%Ic = 0d0
  
  RETURN
  END SUBROUTINE reallocate_rays_arrays

  SUBROUTINE allocSpectrum(alloc_atom_nlte, alloc_image)!NPIX_X, NPIX_Y, N_INCL, N_AZIMUTH)
   !integer, intent(in) :: NPIX_X, NPIX_Y, N_INCL, N_AZIMUTH
   !Polarized quantities allocated in adjustStokesMode
   integer :: nat, k, Nlambda_max
   type (AtomicContinuum) :: cont
   type (AtomicLine) 	  :: line
   logical, intent(in)    :: alloc_atom_nlte, alloc_image
   
   if (allocated(NLTEspec%I)) then
    write(*,*) "Error I already allocated"
    stop
   end if
   
   !stellar radiation; without LD
   !not included yet. Should be the stellar spectrum or
   !should be a BB at a ref T, so that we only need to multiply
   !by a correction depending on the real T of the star i_star
   !if (allocated(NLTEspec%Istar)) allocate(NLTEspec%Istar(NLTEspec%NWAVES))
   
   allocate(NLTEspec%I(NLTEspec%NWAVES, NLTEspec%atmos%NRAYS,NLTEspec%NPROC))
   allocate(NLTEspec%Ic(NLTEspec%NWAVES, NLTEspec%atmos%NRAYS,NLTEspec%NPROC))
   NLTEspec%I = 0d0
   NLTEspec%Ic = 0d0

!    if (NLTEspec%atmos%coherent_scattering) then
!     allocate(NLTEspec%J(NLTEspec%Nwaves,NLTEspec%atmos%Nspace))!%NPROC))
!     allocate(NLTEspec%Jc(NLTEspec%Nwaves,NLTEspec%atmos%Nspace))
!    !Just to try
!    !CALL Warning("(allocSpectrum()) J20 allocated")
!    !allocate(NLTEspec%J20(NLTEspec%Nwaves,NLTEspec%NPROC)); NLTEspec%J20 = 0d0
!     NLTEspec%J = 0.0
!     NLTEspec%Jc = 0.0
!    end if
      
   if (alloc_image) then
    allocate(NLTEspec%Flux(NLTEspec%Nwaves,NPIX_X, NPIX_Y,RT_N_INCL,RT_N_AZ))
    allocate(NLTEspec%Fluxc(NLTEspec%Nwaves,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ))

    NLTEspec%Flux = 0.0
    NLTEspec%Fluxc = 0.0
   end if
   
   !Now opacities
   if (lstore_opac) then !keep continuum LTE opacities in memory
     !sca_c = Kc(:,:,2), chi_c = Kc(:,:,1), eta_c = jc
     allocate(NLTEspec%AtomOpac%Kc(NLTEspec%atmos%Nspace,NLTEspec%Nwaves,2), &
       NLTEspec%AtomOpac%jc(NLTEspec%atmos%Nspace,NLTEspec%Nwaves))
     NLTEspec%AtomOpac%Kc = 0d0
     NLTEspec%AtomOpac%jc = 0d0
   else
    allocate(NLTEspec%AtomOpac%eta_c(NLTEspec%Nwaves ,NLTEspec%NPROC))
    allocate(NLTEspec%AtomOpac%chi_c(NLTEspec%Nwaves ,NLTEspec%NPROC))
    allocate(NLTEspec%AtomOpac%sca_c(NLTEspec%Nwaves,NLTEspec%NPROC))
    NLTEspec%AtomOpac%chi_c = 0.
    NLTEspec%AtomOpac%eta_c = 0.
    NLTEspec%AtomOpac%sca_c = 0.
   end if
   
   allocate(NLTEspec%AtomOpac%chi(NLTEspec%Nwaves ,NLTEspec%NPROC))
   allocate(NLTEspec%AtomOpac%eta(NLTEspec%Nwaves ,NLTEspec%NPROC))
   ! if pol allocate AtomOpac%rho
   NLTEspec%AtomOpac%chi = 0.
   NLTEspec%AtomOpac%eta = 0.
   allocate(NLTEspec%AtomOpac%chic_nlte(NLTEspec%Nwaves, NLTEspec%NPROC),&
      NLTEspec%AtomOpac%etac_nlte(NLTEspec%Nwaves,NLTEspec%NPROC))
   NLTEspec%AtomOpac%chic_nlte = 0d0; NLTEspec%AtomOpac%etac_nlte = 0d0

   allocate(NLTEspec%AtomOpac%eta_p(NLTEspec%Nwaves ,NLTEspec%NPROC))
   allocate(NLTEspec%AtomOpac%chi_p(NLTEspec%Nwaves ,NLTEspec%NPROC))

   NLTEspec%AtomOpac%chi_p = 0.
   NLTEspec%AtomOpac%eta_p = 0.
   
   !do not try to realloc after non-LTE for image.
   !Fursther, with labs=.false. for images, we do not enter in eval_operator condition
   !in NLTEOpacity()
   if (alloc_atom_nlte) then !NLTE loop activated
   
    if (lxcoupling) NLTEspec%atmos%include_xcoupling = .true.
     
    allocate(NLTEspec%Psi(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
    allocate(NLTEspec%dtau(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))   

    do nat=1,NLTEspec%atmos%Nactiveatoms
     allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%eta(NLTEspec%Nwaves,NLTEspec%atmos%Nrays,NLTEspec%NPROC))

     if (NLTEspec%atmos%include_xcoupling) then 
      do k=1, NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Ncont
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%continua(k)%U&
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%continua(k)%Nlambda, NLTEspec%NPROC))
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%continua(k)%chi&
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%continua(k)%Nlambda, NLTEspec%NPROC))
      end do
      do k=1, NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Nline
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%lines(k)%U&
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%lines(k)%Nlambda, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%lines(k)%chi&
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%lines(k)%Nlambda, NLTEspec%atmos%Nrays, NLTEspec%NPROC))      
      end do
     end if

    end do
   end if

   !Contribution function
   if (lcontrib_function .and. ltab_wavelength_image) then !the tab is here to prevent allocating a large array
		!hence ksi should be computed for small waves intervals.
	 !because otherwise, it is allocated with the alloc spectrum before initSpectrumImage()
     if (alloc_image) then 
      allocate(NLTEspec%Ksi(NLTEspec%atmos%Nspace, NLTEspec%Nwaves,RT_N_INCL,RT_N_AZ))  
      allocate(NLTEspec%scale_ref(NLTEspec%atmos%Nspace, RT_N_INCL, RT_N_AZ))
      NLTEspec%Ksi(:,:,:,:) = 0d0
      NLTEspec%scale_ref(:,:,:) = 0d0
     end if
   else if (lcontrib_function .and. .not.ltab_wavelength_image) then
    CALL Warning(" Contribution function not taken into account. Use a wavelength table.")
    lcontrib_function = .false.
   end if

  RETURN
  END SUBROUTINE allocSpectrum

  SUBROUTINE freeSpectrum() 
  
   if (allocated(NLTEspec%Istar)) deallocate(NLTEspec%Istar)
   deallocate(NLTEspec%lambda)
   deallocate(NLTEspec%Ic, NLTEspec%I)
    if (allocated(NLTEspec%Flux)) deallocate(NLTEspec%Flux, NLTEspec%Fluxc)
    
   !if wavelength image do not remove the scattering part ?
!    if (.not.ltab_wavelength_image .and. NLTEspec%atmos%coherent_scattering) then
!     deallocate(NLTEspec%Jc, NLTEspec%J)
! 
!     if (allocated(NLTEspec%J20)) deallocate(NLTEspec%J20)
!    end if

   if (NLTEspec%atmos%Magnetized) then 
    !check allocation due to field_free sol
    if (allocated(NLTEspec%StokesQ)) & !same for all dangerous
    	deallocate(NLTEspec%StokesQ, NLTEspec%StokesU, NLTEspec%StokesV, NLTEspec%F_QUV)
    if (allocated(NLTEspec%AtomOpac%rho_p)) deallocate(NLTEspec%AtomOpac%rho_p)
    if (allocated(NLTEspec%AtomOpac%etaQUV_p)) deallocate(NLTEspec%AtomOpac%etaQUV_p)
    if (allocated(NLTEspec%AtomOpac%chiQUV_p)) deallocate(NLTEspec%AtomOpac%chiQUV_p)
   end if
   
   !active
   deallocate(NLTEspec%AtomOpac%chi)
   deallocate(NLTEspec%AtomOpac%eta)
   deallocate(NLTEspec%AtomOpac%etac_nlte, NLTEspec%AtomOpac%chic_nlte)
   if (allocated(NLTEspec%Psi)) then
    deallocate(NLTEspec%Psi)
    if (allocated(NLTEspec%dtau)) deallocate(NLTEspec%dtau)!, NLTEspec%AtomOpac%initialized)
   end if

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
!    deallocate(Nblue_array, Nmid_array, Nred_array)

   if (allocated(NLTEspec%Ksi)) deallocate(NLTEspec%ksi)
   if (allocated(NLTEspec%scale_ref)) deallocate(NLTEspec%scale_ref)

  RETURN
  END SUBROUTINE freeSpectrum

  SUBROUTINE initAtomOpac(id)!, eval_operator)
    ! set opacities to 0d0 for thread id
    integer, intent(in) :: id
    !logical, intent(in) :: eval_operator !: evaluate operator psi
    
    if (id <= 0) then
     write(*,*) "(initAtomOpac) thread id has to be >= 1!"
     stop
    end if
    
    NLTEspec%AtomOpac%chi(:,id) = 0d0
    NLTEspec%AtomOpac%eta(:,id) = 0d0
    NLTEspec%AtomOpac%chic_nlte(:,id) = 0d0 !ray indep
    NLTEspec%AtomOpac%etac_nlte(:,id) = 0d0
    if (.not.lstore_opac) then
      NLTEspec%AtomOpac%chi_c(:,id) = 0d0
      NLTEspec%AtomOpac%eta_c(:,id) = 0d0
      NLTEspec%AtomOpac%sca_c(:,id) = 0d0
    end if !else thay are not allocated
    NLTEspec%AtomOpac%chi_p(:,id) = 0d0
    NLTEspec%AtomOpac%eta_p(:,id) = 0d0
    
    !Currently LTE or NLTE Zeeman opac are not stored on memory. They change with 
    !direction. BUT the star is assumed to not emit polarised photons (from ZeemanEffect)
    !So we do not take into account this opac in Metal_lambda and futur NLTEOpacity_lambda
    if (NLTEspec%atmos%magnetized .and. PRT_SOLUTION == "FULL_STOKES") then
    !check allocation because even if magnetized, due to the FIELD_FREE solution or WF
    !they might be not allocated !if (allocated(NLTEspec%AtomOpac%rho_p)) 
     NLTEspec%AtomOpac%rho_p(:,:,id) = 0d0
     NLTEspec%AtomOpac%etaQUV_p(:,:,id) = 0d0
     NLTEspec%AtomOpac%chiQUV_p(:,:,id) = 0d0
     !both NLTE and LTE actually.
     !If one want to add them in SEE, it has to be atom (atom%eta) dependent for emissivity.
     !adding the off diagonal elements in SEE results in solving for the whole
     !density matrix, WEEEEEEELLLLL beyond our purpose.
    end if
    
  RETURN
  END SUBROUTINE initAtomOpac
  
  SUBROUTINE init_J_coherent()!(id)
   !integer, intent(in) :: id
  
    NLTEspec%J(:,:) = 0d0
    NLTEspec%Jc(:,:) = 0d0
!     if (allocated(NLTEspec%J20)) NLTEspec%J20(:,id) = 0d0
  
  RETURN
  END SUBROUTINE init_J_coherent
  
  SUBROUTINE init_psi_operator(id, iray)
    integer, intent(in) :: iray, id
    
   	NLTEspec%Psi(:,iray,id) = 0d0
   	NLTEspec%dtau(:,iray,id) = 0d0 !always allocated
  
  RETURN
  END SUBROUTINE init_psi_operator
  
!   SUBROUTINE init_nlte_eta(id)
!    integer :: n
!   
!    do n=1, NLTEspec%atmos%NActiveatoms
!     NLTEspec%atmos%ActiveAtoms(n)%ptr_atom%eta(:,:,id) = 0d0   
!    end do
!   
!   RETURN
!   END SUBROUTINE init_nlte_eta
  
  SUBROUTINE alloc_phi_lambda()
  ! --------------------------------------------------- !
   ! line%phi(line%Nlambda, atmos%Nrays, Nb_proc)
  ! --------------------------------------------------- !
   use atmos_type, only : atmos
   integer :: kr, nact
   
   do nact=1,atmos%Nactiveatoms
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
       		
       !!if I keep phi, I should condiser using it everywhere in Profile() !! ??
       allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%&
       		phi(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%Nlambda,atmos%Nrays,NLTEspec%Nproc))
    end do
   end do
 
  RETURN
  END SUBROUTINE alloc_phi_lambda
  
  SUBROUTINE dealloc_phi_lambda()
  ! --------------------------------------------------- !
   !
  ! --------------------------------------------------- !
   use atmos_type, only : atmos
   integer :: kr, nact
   
   do nact=1,atmos%Nactiveatoms
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
       		
       !!if I keep phi, I should condiser using it everywhere in Profile() !! ??
       deallocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%phi)
    end do
   end do
 
  RETURN
  END SUBROUTINE dealloc_phi_lambda


 SUBROUTINE WRITE_FLUX()
 ! -------------------------------------------------- !
  ! Write the spectral Flux map on disk.
  ! FLUX map:
  ! NLTEspec%Flux total and NLTEspec%Flux continuum
 ! --------------------------------------------------- !
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(6) :: naxes
  integer :: group,fpixel,nelements, i, xcenter
  integer :: la, Nred, Nblue, kr, kc, m, Nmid
  logical :: simple, extend
  character(len=6) :: comment="VACUUM"
  double precision :: lambda_vac(NLTEspec%Nwaves)
  real :: pixel_scale_x, pixel_scale_y 
  
  write(*,*) "Writing Flux-map"
  write(*,*) "npix_x = ", npix_x, " npix_y = ", npix_y, ' RT method:', RT_line_method
  write(*,*) "Wavelength points:", NLTEspec%Nwaves
  
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

  !Method1 not sure of what it does
  !error if Keplerian or infall or any ?
  if ((lkeplerian .or. linfall .or. lmagnetoaccr).and.(l_sym_ima)) &
  	write(*,*) "Warning, image symmetry might be wrong."
  if (l_sym_ima.and.RT_line_method == 2) then 
   xcenter = npix_x/2 + modulo(npix_x,2)
   do i=xcenter+1,npix_x
     NLTEspec%Flux(:,i,:,:,:) = NLTEspec%Flux(:,npix_x-i+1,:,:,:)
   end do
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
  
  if (l_sym_ima.and.(RT_line_method == 2)) then
   xcenter = npix_x/2 + modulo(npix_x,2)
   do i=xcenter+1,npix_x
    NLTEspec%Fluxc(:,i,:,:,:) = NLTEspec%Fluxc(:,npix_x-i+1,:,:,:)
   end do
  end if ! l_sym_image  
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%Fluxc,status)
  ! write polarized flux if any. Atmosphere magnetic does not necessarily
  								!means we compute polarization
  if ((NLTEspec%atmos%magnetized) .and. (PRT_SOLUTION /= "NO_STOKES") &
               .and. (RT_line_method /= 1)) then
   write(*,*) " -> Writing polarization"
   CALL ftcrhd(unit, status)
   naxis = 6
   naxes(1) = 3 !Q, U, V
   naxes(2)=NLTEspec%Nwaves
   naxes(3)=npix_x
   naxes(4)=npix_y
   naxes(5)=RT_n_incl
   naxes(6)=RT_n_az
   nelements = naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5) * naxes(6)
   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
   CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'Polarised Flux (Q, U, V)',status)
   CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%F_QUV,status)
  end if
  
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
   !
   ! To.Do: 
   ! Write for each transition of each atom the part of the grid
   ! associated to this transition.
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
  
 !building 
 SUBROUTINE WRITE_CNTRB_FUNC_PIX()
 ! -------------------------------------------------- !
  ! Write contribution function to disk.
 ! --------------------------------------------------- !
  integer :: status,unit,blocksize,bitpix,naxis, naxis2
  integer, dimension(7) :: naxes, naxes2
  integer :: group,fpixel,nelements, nelements2
  logical :: simple, extend
  character(len=6) :: comment=""
  real :: pixel_scale_x, pixel_scale_y 
  
   write(*,*)" -> writing contribution function.."
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
   status=0
   CALL ftgiou (unit,status)

   !  Create the new empty FITS file.
   blocksize=1
   CALL ftinit(unit,"cntrb.fits.gz",blocksize,status)

   simple=.true.
   extend=.true.
   group=1
   fpixel=1

   bitpix=-64
   
  if (lVoronoi) then
   naxis = 4
   naxes(1) = NLTEspec%atmos%Nspace ! equivalent n_cells
   naxes(2) = NLTEspec%Nwaves
   naxes(3) = RT_n_incl
   naxes(4) = RT_n_az
   nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
   naxis2 = 3
   naxes2(1) = NLTEspec%atmos%Nspace ! equivalent n_cells
   naxes2(2) = RT_n_incl
   naxes2(3) = RT_n_az
   nelements2 = naxes2(1) * naxes2(2) * naxes2(3)
  else
   if (l3D) then
    naxis = 6
    naxes(1) = n_rad
    naxes(2) = 2*nz
    naxes(3) = n_az
    naxes(4) = NLTEspec%Nwaves
    naxes(5) = RT_n_incl
    naxes(6) = RT_n_az
    nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4) * naxes(5) * naxes(6)
    naxis2 = 5
    naxes2(1) = n_rad
    naxes2(2) = 2*nz
    naxes2(3) = n_az
    naxes2(4) = RT_n_incl
    naxes2(5) = RT_n_az
    nelements2 = naxes2(1) * naxes2(2) * naxes2(3) * naxes2(4) * naxes2(5)
   else
    naxis = 5
    naxes(1) = n_rad
    naxes(2) = nz
    naxes(3) = NLTEspec%Nwaves
    naxes(4) = RT_n_incl
    naxes(5) = RT_n_az
    nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4) * naxes(5) 
    naxis2 = 4
    naxes2(1) = n_rad
    naxes2(2) = nz
    naxes2(3) = RT_n_incl
    naxes2(4) = RT_n_az
    nelements2 = naxes2(1) * naxes2(2) * naxes2(3) * naxes2(4)
   end if
  end if

  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  CALL ftpkys(unit,'UNIT',"W.m-2.Hz-1.pixel-1",'Ksi',status)


  !  Write the array to the FITS file.
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%ksi,status)
  
  !open new hdu for reference scale
  CALL ftcrhd(unit, status)
  
  CALL ftphpr(unit,simple,bitpix,naxis2,naxes2,0,1,extend,status)
  CALL ftpkys(unit, "Optical depth", "", "", status)
  CALL ftpprd(unit,group,fpixel,nelements2,NLTEspec%scale_ref,status)

  !  Close the file and free the unit number.
  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     CALL print_error(status)
  endif

 RETURN
 END SUBROUTINE WRITE_CNTRB_FUNC_PIX

END MODULE spectrum_type

