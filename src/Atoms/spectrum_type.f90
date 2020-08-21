module spectrum_type

	use atom_type, only : AtomicLine, AtomicContinuum, AtomType
	use atmos_type, only : helium, hydrogen, NactiveAtoms, Natom, Atoms, v_char, icompute_atomRT
	use getlambda, only  : hv, adjust_wavelength_grid, Read_wavelengths_table, make_wavelength_grid_new!, make_wavelength_grid, 
	use math, only : linear_1D
	use fits_utils, only : print_error
	use parametres, only : n_cells, lelectron_scattering, n_etoiles, npix_x, npix_y, rt_n_incl, rt_n_az, &
							lcontrib_function, n_rad, n_az, nz, map_size, distance, zoom, lmagnetoaccr, lmagnetic_field, prt_solution, &
							lvacuum_to_air, ltab_wavelength_image, lvoronoi, l3D
	use input, only : nb_proc, RT_line_method,lkeplerian, linfall, l_sym_ima
	use constantes, only : arcsec_to_deg
	use mcfost_env, only : dp
	use messages, only : error, warning
	!use hdf5

	implicit none

	real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
	real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352
	character(len=*), parameter :: WAVES_FILE="atom_transfer_waves_grid.fits.gz"
   ! not that FLUX_FILE is 1D only if one pixel, otherwise it is a map.
   ! F(x,y,0,lambda) in several directions.
	character(len=*), parameter :: FLUX_FILE="flux.fits.gz", CF_FILE="cntrb.fits.gz"
	character(len=*), parameter :: FLUX_FILE_H5="flux.h5", CF_FILE_H5="cntrb.h5"
	
	!shift in index of line profiles, function of iray and nb_proc
	integer, dimension(:,:), allocatable :: dk
	integer :: dk_max, dk_min

	integer :: Nlambda, Ntrans, Nlambda_cont
	real(kind=dp) :: wavelength_ref = 0d0
	real(kind=dp), dimension(:,:), allocatable :: chi_c, sca_c, eta_c, chi, eta, chi_c_nlte, eta_c_nlte, chi0_bb, eta0_bb
	real(kind=dp), dimension(:,:,:), allocatable :: Icont
	real(kind=dp), dimension(:), allocatable :: lambda, lambda_cont
	real(kind=dp), dimension(:,:), allocatable :: Jnu, Jnu_cont, eta_es
	real(kind=dp), dimension(:,:), allocatable :: Istar_tot, Istar_cont
	real(kind=dp), dimension(:,:,:), allocatable :: Itot
	real(kind=dp), allocatable, dimension(:,:,:,:,:) :: Flux !incl,az, npix,npix, lambda, Fluxc
	real(kind=dp), allocatable, dimension(:,:,:,:,:) :: Fluxc
														!If stored for each trans in a file, store Fluxc only at lambda0   
	real(kind=dp), dimension(:,:,:), allocatable :: rho_p, chiQUV_p, etaQUV_p
	real(kind=dp), allocatable, dimension(:,:) :: Stokes_Q, Stokes_U, Stokes_V
	real(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: F_QUV
	real(kind=dp), allocatable, dimension(:,:,:,:) :: Ksi 
   ! Flux is a map of Nlambda, xpix, ypix, nincl, nazimuth
   
!! 	real(kind=dp), allocatable, dimension(:,:,:) :: Psi, Stot, chitot
	
!   TYPE AtomicOpacity
!    !active opacities
!    real(kind=dp), allocatable, dimension(:,:)   :: chi, eta!, chic_nlte, etac_nlte
!    ! NLTE magneto-optical elements and dichroism are stored in the background _p arrays.
!    ! Mainly because we do not use them in SEE, even if etaQUV can be added to the total
!    ! emissivity. But in case, etaQUV has to be atom dependent, so now we can store the LTE
!    ! and NLTE is the same variable
!    !passive opacities
!    real(kind=dp), allocatable, dimension(:,:)   :: eta_p, chi_p
! !    real(kind=dp), allocatable, dimension(:,:)   :: eta_c, chi_c!, sca_c
!    real(kind=dp), allocatable, dimension(:,:,:)   :: rho_p, chiQUV_p, etaQUV_p
!    real(kind=dp), allocatable, dimension(:,:)   :: jc, jc_nlte, Kc_nlte
!    real(kind=dp), allocatable, dimension(:,:) :: Kc, sca_c
! 
!   end TYPE AtomicOpacity
! 
!   TYPE Spectrum
!    !n_proc should be the last index
!    type  (GridType), pointer :: atmos
!    logical :: vacuum_to_air=.false., write_wavelength_grid=.false.
!    integer :: Nwaves, Nact, Npass, Ntrans, NPROC=1, Nwaves_cont
!    real(kind=dp) :: wavelength_ref=0.d0 !nm optionnal
!    real(kind=dp), allocatable, dimension(:) :: lambda, lambda_cont
!    !nlambda, nrays, nproc
!    real(kind=dp), allocatable, dimension(:,:,:) :: I, Ic, Stokes_Q, Stokes_U, Stokes_V
!    real(kind=dp), allocatable, dimension(:,:) :: Istar
!    !nlambda, nproc
!    real(kind=dp), allocatable, dimension(:,:) :: J, Jc, J20
!    !Nlambda, xpix, ypix, Nincl, Naz
!    real(kind=dp), allocatable, dimension(:,:,:,:,:) :: Flux, Fluxc
!    real(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: F_QUV
!    !!real(kind=dp), allocatable, dimension(:,:) :: S_QUV
!    !Contribution function
!    !Nlambda,N_INCL, N_AZ, NCELLS
!    real(kind=dp), allocatable, dimension(:,:,:,:) :: Ksi 
!    ! Flux is a map of Nlambda, xpix, ypix, nincl, nazimuth
!    real(kind=dp), allocatable, dimension(:,:,:) :: Psi, etau, S, chi
!    real(kind=dp), allocatable, dimension(:,:) :: Jext, Jloc
!    !size of Psi could change during the devlopment
!    type (AtomicOpacity) :: AtomOpac
!    character:: Jfile, J20file
!   end TYPE Spectrum
! 
!   type (Spectrum) :: NLTEspec

	contains
  
 
	function getlambda_limit(atom)
	!reddest wavelength of the shortest b-b transition in the ground state
	real(kind=dp) :: getlambda_limit
	type(AtomType) :: atom
	integer :: kr
  
	getlambda_limit = 1d6
	do kr=1, atom%Nline

		if (atom%lines(kr)%i==1) then !transition to the ground state
			getlambda_limit = min(getlambda_limit, lambda(atom%lines(kr)%Nred))
		endif
  
	enddo
  
	!avoid problem
	getlambda_limit = max(lambda(1), getlambda_limit)
  
	return
	end function getlambda_limit

	subroutine init_Spectrum(Nray, lam0, vacuum_to_air)
	! ------------------------------------------- !
    ! Allocate and store the wavelength grid for
    ! NLTE atomic line transfer.
	! ------------------------------------------- !
		integer, intent(in) :: Nray
		integer :: kr, nat
		real(kind=dp), optional :: lam0
		logical, optional :: vacuum_to_air
		logical :: alloc_nlte_vars, alloc_image_vars
    
		alloc_nlte_vars = Nactiveatoms > 0

   
		if (present(lam0)) wavelength_ref = lam0
  
   ! Initialize the wavelength grid depending on the number of transitions PASSIVE/ACTIVE
   ! and allocate NLTEspec%lambda
   ! Set also Nblue and Nred for each transitions: extension of the transtion
   ! in absolute indexes on the whole grid.
!    call make_wavelength_grid(NLTEspec%atmos%Natom, NLTEspec%atmos%Atoms, & 
!                         NLTEspec%lambda, NLTEspec%Ntrans, NLTEspec%wavelength_ref)
                        
! 		call make_wavelength_grid_new(atmos%Natom, atmos%Atoms, wavelength_ref, &
!    										atmos%v_char, lambda, Ntrans, lambda_cont)
		call make_wavelength_grid_new(wavelength_ref, v_char, lambda, Ntrans, lambda_cont)

		!for each line
		dk_max = sign(1.0_dp, v_char) * int( 1e-3 * abs(v_char) / hv + 0.5 ) !nint( (1e-3 * v_char) / hv)
		dk_min = -dk_max
		write(*,*) "Maximum shift in index:", dk_max, (1e-3 * v_char + hv) / hv - 1.0
		if (dk_max * hv > v_char*1e-3 + hv) then
			call warning("Beware, maximum shift might be beyond the grid !")
		endif
! 		write(*,*) dk_max * hv, v_char*1e-3, v_char*1e-3 + hv
! 		stop
		Nlambda_cont = size(lambda_cont)
		Nlambda = size(lambda)
   !Futur deprecation, rayleigh scattering will undergo a revamp, and some informations
   !will be hardcoded. Atm, Rayleigh scattering is deactivated 
   !only for H and He atm, not changed for image even if we remove lines
		hydrogen%scatt_limit = getlambda_limit(hydrogen)
		if (associated(helium)) helium%scatt_limit = getlambda_limit(helium)
   
		call writeWavelength()
		call alloc_spectrum(alloc_nlte_vars, Nray)
		

	return
	end subroutine init_Spectrum
  

	subroutine init_Spectrum_Image()
	! -------------------------------------------------------------------- !
	! Allocate a special wavelength grid for emission_line map.
	! This allows to solve for the RT equation only for a specific line.
	! This fasten LTE calculation in case of 1 line.
	! -------------------------------------------------------------------- !
	real(kind=dp), dimension(Nlambda) :: old_grid
	
	integer, parameter :: Nray = 1 !only one for image !!
	
	integer, dimension(:), allocatable :: Nlam_R
	integer :: nat, kr, Ntrans_new, kp, kc
	!for Jnu
	integer :: icell, alloc_status, Nwaves_old
	real(kind=dp), dimension(:,:), allocatable :: Jnuo, Jnuoc

call error("initSpectrumImage not modified!!")
   
	old_grid = lambda
	write(*,*) " -> Redefining a wavelength grid for image.."
	call dealloc_spectrum() !first free waves arrays

   !Reset the differents Nblue, Nred for atomic transitions
   !and recompute photoionisation cross-section
	call Read_wavelengths_table(lambda, Nlam_R)
	Nlambda = size(lambda)
   !Works for Active and Passive atoms alike.
   !-> sort also the grid
   !-> make sure to search for transitions in individual regions
   ! and not in the total grid ! Otherwise, the transitions are kept if they fall
   ! between lambda(min) and lambda(max)
	call adjust_wavelength_grid(old_grid, lambda, Nlam_R)
	deallocate(Nlam_R) !not used anymore
	write(*,*) "check dkmin/max"
	dk_max = 1
	dk_min = 1
   
   !Reallocate J and Jc keeping their value if electron_scattering
   !NB: the spatial grid is the same
   !if (NLTEspec%atmos%electron_scattering) then
   ! Jold = J; Joldc = Jc
   ! dealloc(J, Jc)
   ! allocate(J, Jc)
   ! J = PACK(Jold, mask=NLTEspec%mlamba==old_grid) !may not work directly
   !end if
   
	Ntrans_new = 0
	write(*,*) " -> Using ", Nlambda," wavelengths for image and spectrum."
!I do not removed transitions for the moment
!so I count only transitions falling in the new wavelength. Thi summation is implicit
!if transitions are removed.
	do nat=1,Natom
		write(*,*) "  --> Atom ", Atoms(nat)%ptr_atom%ID!, &
!      NLTEspec%atmos%Atoms(nat)%ptr_atom%Nline, "b-b", &
!      NLTEspec%atmos%Atoms(nat)%ptr_atom%Ncont, "b-f"
		kp = 0
     
     !do sum over transitions
		do kr=1,atoms(nat)%ptr_Atom%Ntr
			kc = atoms(nat)%ptr_Atom%at(kr)%ik
			select case (atoms(nat)%ptr_Atom%at(kr)%trtype)
      
			case ("ATOMIC_LINE")
				write(*,*) "    b-b #", kc, Atoms(nat)%ptr_atom%lines(kc)%lambda0,"nm"

			case ("ATOMIC_CONTINUUM")
				write(*,*) "    b-f #", kc, Atoms(nat)%ptr_atom%continua(kc)%lambdamin, &
					"nm", Atoms(nat)%ptr_atom%continua(kc)%lambda0, "nm"       
			case default
				call error ("transition type unkown",atoms(nat)%ptr_Atom%at(kr)%trtype)
			end select
      
		end do
		Ntrans_new = Ntrans_new + atoms(nat)%ptr_atom%Ntr
	end do

	Ntrans = Ntrans_new

	write(*,*) "  ** Total number of transitions for image:", Ntrans

	call alloc_spectrum(.false., Nray)


	if (lelectron_scattering) then
		write(*,*)  " -> Resample mean radiation field for isotropic scattering..." 

		Nwaves_old = size(old_grid)
		allocate(Jnuo(Nwaves_old, n_cells))
		Jnuo(:,:) = Jnu(:,:)
		allocate(Jnuoc(Nwaves_old, n_cells))
		Jnuoc(:,:) = Jnu_cont(:,:)
   
		deallocate(Jnu, Jnu_cont)
   
		!new freq grid
		allocate(Jnu(Nlambda, n_cells),stat=alloc_status)
		if (alloc_status>0) call error("Allocation error, resample J")
		allocate(Jnu_cont(Nlambda_cont, n_cells),stat=alloc_status)
		if (alloc_status>0) call error("Allocation error, resample Jc")

		Jnu(:,:) = 0d0
		Jnu_cont(:,:) = 0d0
		!use better interpolation ??
		do icell=1, n_cells
			if (icompute_atomRT(icell)>0) then
				Jnu(:,icell) = Linear_1D(Nwaves_old, old_grid,Jnuo(:,icell),Nlambda,lambda)
				Jnu_cont(:,icell) = Linear_1D(Nwaves_old, old_grid,Jnuoc(:,icell),Nlambda_cont,lambda_cont)
			endif
		enddo  
	endif !resample Jnu
	write(*,*)  " ...Done!" 

   
	return
	end subroutine init_Spectrum_Image

  
	subroutine reallocate_rays_arrays(newNray)
		integer, intent(in) :: newNray
		integer :: kr, nat, alloc_status, oldRay
  
   		oldRay = size(Itot(1,:,1))
   		
		if (newNray /=1 ) call Warning("  Beware, check the number of rays for the ray-traced map!")

   
		deallocate(Itot, Icont)
    
		allocate(Itot(Nlambda, newNray, nb_proc))
		allocate(Icont(Nlambda_cont, newNray, nb_proc))
		Itot = 0.0_dp; Icont = 0.0_dp
  
	return
	end subroutine reallocate_rays_arrays

	subroutine alloc_Spectrum(alloc_atom_nlte, Nray)
		!Polarized quantities allocated in adjustStokes_Mode
		integer, intent(in) :: Nray
		integer :: nat, k, Nlambda_max, alloc_status, istar
		type (AtomType), pointer :: atom
		logical, intent(in)    :: alloc_atom_nlte
   
		if (allocated(Itot)) then
			write(*,*) "Error I already allocated"
			stop
		end if

		if (n_etoiles > 0) then
			allocate(Istar_tot(Nlambda,n_etoiles))
			Istar_tot(:,:) = 0.0_dp
			allocate(Istar_cont(Nlambda_cont,n_etoiles))
			Istar_cont(:,:) = 0.0_dp
		endif

		allocate(Itot(Nlambda, Nray, nb_proc),stat=alloc_status)
		if (alloc_status > 0) call error("Allocation error Itot")
		allocate(Icont(Nlambda_cont, Nray, nb_proc))
		Itot = 0.0_dp
		Icont = 0.0_dp
	
		!allocate(dk(Nray, nb_proc))
		!if (alloc_status > 0) call error("Allocation error dk")

		allocate(chi_c(Nlambda_cont,n_cells), eta_c(Nlambda_cont,n_cells), stat=alloc_status)
		!-> At the moment not needed because only Thomson scattering included
		!allocate(sca_c(Nlambda_cont,n_cells), stat=alloc_status)

		if (alloc_status > 0) then
			write(*,*) " mem = ", real(3 * n_cells * Nlambda_cont), " GB"
			call error("Allocation error, continuum opacities")
		endif
		chi_c = 0.0_dp
		eta_c = 0.0_dp
		if (allocated(sca_c)) sca_c = 0.0_dp


		allocate(eta(Nlambda ,nb_proc))
		allocate(chi(Nlambda ,nb_proc))

		eta(:,:) = 0.0_dp
		chi(:,:) = 0.0_dp
		
		!interpolated total continuum opacities on the lambda grid to be used with lines
		allocate(eta0_bb(Nlambda , n_cells))
		allocate(chi0_bb(Nlambda , n_cells))
		
		eta0_bb = 0.0_dp
		chi0_bb = 0.0_dp
		
		!otherwise allocated bellow for nlte.
		if ((lelectron_scattering).and.(.not.alloc_atom_nlte)) call alloc_jnu


		if (alloc_atom_nlte) then !NLTE loop activated
			allocate(chi_c_nlte(Nlambda_cont,n_cells),stat=alloc_status)
			if (alloc_status >0) call error("Allocation error chi_c_nlte")
			allocate(eta_c_nlte(Nlambda_cont,n_cells),stat=alloc_status)
			if (alloc_status >0) call error("Allocation error eta_c_nlte")
			chi_c_nlte = 0.0_dp
			eta_c_nlte = 0.0_dp
			
			call alloc_jnu
                			
			! .....  add other NLTE opac
			do nat=1, Natom
				atom => Atoms(nat)%ptr_atom
				do k=1, atom%Nline
					allocate(atom%lines(k)%phi_loc(atom%lines(k)%Nred-dk_min+dk_max+1-atom%lines(k)%Nblue,Nray,nb_proc),stat=alloc_status)
					if (alloc_status > 0) call error("Allocation error line%phi_loc for nlte loop")
				enddo
			enddo

		endif

	return
	end subroutine alloc_Spectrum
	
	subroutine alloc_flux_image
		integer :: alloc_status
		real :: mem_alloc, mem_cont
   
   		!Continuum image and continuum wavelengths not written atm but still computed for debug
   		mem_cont = real(npix_x*npix_y)*real(rt_n_incl*rt_n_az)*real(Nlambda_cont)/1024./1024. + real(Nlambda_cont)/1024./1024.
		mem_alloc = real(npix_x*npix_y)*real(rt_n_incl*rt_n_az)*real(Nlambda)/1024./1024. + real(Nlambda)/1024./1024.
		!Flux  + Flux cont + lines grid
		if (mem_alloc + mem_cont > 1d3) then !in MB
			write(*,*) " allocating ",( mem_alloc + mem_cont )/1024., " GB for flux arrays.."
		else
			write(*,*) " allocating ", mem_alloc + mem_cont, " MB for flux arrays.."  
		endif
		write(*,*) "  -> ", mem_alloc, " MB for fits file"

		allocate(Flux(Nlambda,NPIX_X, NPIX_Y,RT_N_INCL,RT_N_AZ), stat=alloc_status)
		if (alloc_Status > 0) call ERROR ("Cannot allocate Flux")
		allocate(Fluxc(Nlambda_cont,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ), stat=alloc_status)
		if (alloc_Status > 0) call ERROR ("Cannot allocate Flux continuum")

		Flux = 0.0_dp
		Fluxc = 0.0_dp
    
		!Contribution function
   
		if (lcontrib_function) then

			mem_alloc = real(n_cells,kind=dp)/1024. * real(Nlambda,kind=dp)/1024. * &
						real(RT_N_INCL*RT_N_AZ,kind=dp) !in MB
	 
			if (mem_alloc > 1d3) then
				write(*,*) " allocating ", mem_alloc, " GB for contribution function.."
			else
				write(*,*) " allocating ", mem_alloc, " MB for contribution function.."
			endif 
      
			if (mem_alloc >= 2.1d3) then !2.1 GB
				call Warning(" To large Ksi array. Use a wavelength table instead..")
				lcontrib_function = .false.
			else
      
				allocate(Ksi(Nlambda,n_cells,RT_N_INCL,RT_N_AZ),stat=alloc_status)
				if (alloc_status > 0) then
					call ERROR('Cannot allocate ksi')
					lcontrib_function = .false.
				else
					Ksi(:,:,:,:) = 0d0
				endif

			end if

		end if

   
	return
	end subroutine alloc_flux_image

	subroutine dealloc_spectrum() 
  
		deallocate(lambda)
		if (allocated(lambda_cont)) deallocate(lambda_cont)

		deallocate(Itot, Icont, Istar_tot, Istar_cont)
   	!can be deallocated before to save memory 
		if (allocated(Flux)) deallocate(Flux)
		if (allocated(Fluxc)) deallocate(Fluxc)
		
		if (allocated(dk)) deallocate(dk)

    !check allocation due to field_free sol
		if (allocated(Stokes_Q)) then !same for all dangerous
    		deallocate(Stokes_Q, Stokes_U, Stokes_V, F_QUV)
			deallocate(rho_p)
			deallocate(etaQUV_p)
			deallocate(chiQUV_p)
		endif
   
! 		if (allocated(psi)) deallocate(Psi)
! 		if (allocated(stot)) deallocate(Stot)
! 		if (allocated(chitot)) deallocate(chitot)

		deallocate(chi_c,  eta_c)
		if (allocated(sca_c)) deallocate(sca_c)
		deallocate(chi, eta)
		deallocate(chi0_bb, eta0_bb)

   !elsewhere
   !if (NLTEspec%Nact > 0) deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)


	!Can be deallocated before to save memory
		if (allocated(ksi)) deallocate(ksi)


	return
	end subroutine dealloc_Spectrum

! 	subroutine initAtomOpac(id)
! 		integer, intent(in) :: id
!     
! 
! 		NLTEspec%AtomOpac%chi_p(:,id) = 0d0
! 		NLTEspec%AtomOpac%eta_p(:,id) = 0d0
! 
!     
! 	return
! 	end subroutine initAtomOpac
  
!   subroutine initAtomOpac_nlte(id)!, eval_operator)
!     ! set opacities to 0d0 for thread id
!     integer, intent(in) :: id
!     !logical, intent(in) :: eval_operator !: evaluate operator psi
!     
!     !need to be sure that id is > 0
! !     if (id <= 0) then
! !      write(*,*) "(initAtomOpac) thread id has to be >= 1!"
! !      stop
! !     end if
! 
!     NLTEspec%AtomOpac%chi(:,id) = 0d0
!     NLTEspec%AtomOpac%eta(:,id) = 0d0
!     
! !     NLTEspec%AtomOpac%chic_nlte(:,id) = 0d0 
! !     NLTEspec%AtomOpac%etac_nlte(:,id) = 0d0
! 
! 
!     
!   return
!   end subroutine initAtomOpac_nlte
  
!   subroutine initAtomOpac_zeeman(id)!, eval_operator)
!     ! set opacities to 0d0 for thread id
!     integer, intent(in) :: id
!     !logical, intent(in) :: eval_operator !: evaluate operator psi
!     
!     !need to be sure that id is > 0
! !     if (id <= 0) then
! !      write(*,*) "(initAtomOpac) thread id has to be >= 1!"
! !      stop
! !     end if
! 
! 
!     
!     !Currently LTE or NLTE Zeeman opac are not stored on memory. They change with 
!     !direction. BUT the star is assumed to not emit polarised photons (from ZeemanEffect)
!     !So we do not take into account this opac in Metal_lambda and futur NLTEOpacity_lambda
!     !if (NLTEspec%atmos%magnetized .and. PRT_SOLUTION == "FULL_Stokes_") then
!     !check allocation because even if magnetized, due to the FIELD_FREE solution or WF
!     !they might be not allocated !if (allocated(NLTEspec%AtomOpac%rho_p)) 
!     NLTEspec%AtomOpac%rho_p(:,:,id) = 0d0
!     NLTEspec%AtomOpac%etaQUV_p(:,:,id) = 0d0
!     NLTEspec%AtomOpac%chiQUV_p(:,:,id) = 0d0
!      !both NLTE and LTE actually.
!      !If one want to add them in SEE, it has to be atom (atom%eta) dependent for emissivity.
!      !adding the off diagonal elements in SEE results in solving for the whole
!      !density matrix, WEEEEEEELLLLL beyond our purpose.
!     !end if
!     
!   return
!   end subroutine initAtomOpac_zeeman
  
	subroutine alloc_Jnu()
  
		!allocate(Jnu(Nlambda,n_cells))
		!Jnu = 0.0
		allocate(eta_es(Nlambda,n_cells))
		eta_es = 0.0

		allocate(Jnu_cont(Nlambda_cont,n_cells))
		Jnu_cont = 0.0

   
	return
	end subroutine alloc_Jnu
  
	subroutine dealloc_Jnu()
    
		if (allocated(Jnu)) deallocate(Jnu)
		if (allocated(Jnu_cont)) deallocate(Jnu_cont)
		if (allocated(eta_es)) deallocate(eta_es)

	return 
	end subroutine dealloc_Jnu
  
! 	subroutine init_psi_operator_m(id)!, iray)
! 		integer, intent(in) :: id!, iray
! 		integer :: nact
!     
! 		NLTEspec%Psi(:,:,id) = 0d0
! 		NLTEspec%etau(:,:,id) = 0d0
!    	
! 		do nact=1,NLTEspec%atmos%NactiveAtoms
! !works only for MALI
! 			NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%eta(:,:,id) = 0d0
! 			NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%etac(:,id) = 0d0
! 			NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%chic(:,id) = 0d0
! 		enddo
!   
! 	return
! 	end subroutine init_psi_operator_m

! 	subroutine init_psi_operator(id)
! 		integer, intent(in) :: id
! 		integer :: nact
!     
! 		NLTEspec%Psi(:,:,id) = 0d0
! 		NLTEspec%etau(:,:,id) = 0d0
! 		NLTEspec%S(:,:,id) = 0d0
! 
! !Only if for sub iterations, I compute the local nlte cont for iray==1 only.
! !Or if it is not too costly I recompute them for each ray, even if ray indep..   	
! ! 		do nact=1,NLTEspec%atmos%NactiveAtoms
! ! 			NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%etac(:,id) = 0d0
! ! 			NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%chic(:,id) = 0d0
! ! 		enddo
!   
! 	return
! 	end subroutine init_psi_operator
	
	
! 	subroutine init_Iext_line(iray,id)
! 		integer, intent(in) :: iray, id
! 		integer :: nact, kr, Nb, Nr
! 		real(kind=dp) :: dk
! 		
! 		do nact=1, NLTEspec%atmos%NActiveatoms
! 		
! 			do kr=1, NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%Nline
! 				dk = NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%dk(iray,id)
! 
! 				Nb = NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%Nblue
! 				Nr = NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%Nred
! 				
! 				NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%I_ij(:,iray,id) = &
! 					NLTEspec%I(Nb+dk:Nr+dk,iray,id) !!* exp(-)
! 			
! 			enddo
! 		
! 		enddo
! 	
! 	return
! 	end subroutine init_Iext_line

  
!   subroutine alloc_weights()
!   ! --------------------------------------------------- !
!    ! 
!   ! --------------------------------------------------- !
!    use atmos_type, only : atmos
!    integer :: kr, nact
!    type(AtomicLine) :: line
!    type(AtomicContinuum) :: cont
!    
!    do nact=1,atmos%Nactiveatoms
!     do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
!        line = atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)
!        allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%w_lam(line%Nlambda))
! !        allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%phi_ray(line%Nlambda,NLTEspec%atmos%Nrays, nb_proc))
! 
!     end do
!     do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Ncont
!        cont = atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)
!        allocate(atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)%w_lam(cont%Nlambda))
!     end do
!    end do
!  
!   return
!   end subroutine alloc_weights
  
!   subroutine dealloc_weights()
!   ! --------------------------------------------------- !
!    ! 
!   ! --------------------------------------------------- !
!    use atmos_type, only : atmos
!    integer :: kr, nact
!    
!    do nact=1,atmos%Nactiveatoms
!     do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
!        if (allocated(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%w_lam)) &
!        	deallocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%w_lam)
! !        deallocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%phi_ray)
!     end do
!     do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Ncont
!        deallocate(atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)%w_lam)
!     end do
!    end do
!  
!   return
!   end subroutine dealloc_weights



	subroutine write_image()
	! -------------------------------------------------- !
	! Write the spectral Flux map on disk.
	! FLUX map:
	! NLTEspec%Flux total and NLTEspec%Flux continuum
	! --------------------------------------------------- !
  !!use input
	integer :: status,unit,blocksize,bitpix,naxis
	integer, dimension(6) :: naxes
	integer :: group,fpixel,nelements, i, xcenter
	integer :: la, Nred, Nblue, kr, kc, m, Nmid
	logical :: simple, extend
	character(len=6) :: comment="VACUUM"
	real(kind=dp) :: lambda_vac(Nlambda), Fnu, lambdac_vac(Nlambda_cont)
	real :: pixel_scale_x, pixel_scale_y 
  
	write(*,*) "Writing Flux-map"
	write(*,*) "npix_x = ", npix_x, " npix_y = ", npix_y, ' RT method:', RT_line_method
	write(*,*) "Wavelength points:", Nlambda
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
	status=0
	call ftgiou (unit,status)

   !  Create the new empty FITS file.
	blocksize=1
	call ftinit(unit,trim(FLUX_FILE),blocksize,status)

	simple=.true.
	extend=.true.
	group=1
	fpixel=1

	bitpix=-64
	naxis=5
	naxes(1)=Nlambda!1!1 if only one wavelength

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
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	if (status > 0) then
		call print_error(status)
	endif

   !!RAC, DEC, reference pixel & pixel scale en degres
	call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
	call ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
	call ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
	pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
	call ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)
 
	call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
	call ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
	call ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
	pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
	call ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)

	call ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'F_nu',status)

	if ((lkeplerian .or. linfall .or. lmagnetoaccr).and.(l_sym_ima)) &
		write(*,*) "Warning, image symmetry might be wrong."
		if (l_sym_ima.and.RT_line_method == 2) then 
			xcenter = npix_x/2 + modulo(npix_x,2)
			do i=xcenter+1,npix_x
				Flux(:,i,:,:,:) = Flux(:,npix_x-i+1,:,:,:)
			end do
	end if ! l_sym_image

	!  Write the array to the FITS file.
	call ftpprd(unit,group,fpixel,nelements,Flux,status)
	if (status > 0) then
		call print_error(status)
	endif

!-> Continuum map not written ATM.
  ! create new hdu for continuum
!   call ftcrhd(unit, status)
!   if (status > 0) then
!      call print_error(status)
!   endif
! 
!   !naxis(1) = NLTEspec%Nwaves_cont
!   !nelements = naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
! 
!   !  Write the required header keywords.
!   call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!   call ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
!   call ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
!   call ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
!   pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
!   call ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)
!  
!   call ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
!   call ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
!   call ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
!   pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
!   call ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)
!   call ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'F_nu',status)
!   
!   if (l_sym_ima.and.(RT_line_method == 2)) then
!    xcenter = npix_x/2 + modulo(npix_x,2)
!    do i=xcenter+1,npix_x
!     NLTEspec%Fluxc(:,i,:,:,:) = NLTEspec%Fluxc(:,npix_x-i+1,:,:,:)
!    end do
!   end if ! l_sym_image
!     
!
!	write continuum wavelengths
!
!
!   call ftpprd(unit,group,fpixel,nelements,NLTEspec%Fluxc,status)
!   if (status > 0) then
!      call print_error(status)
!   endif
!   
  
  ! write polarized flux if any. Atmosphere magnetic does not necessarily
  								!means we compute polarization
	if ((lmagnetic_field) .and. (PRT_SOLUTION /= "NO_STOKES") .and. (RT_line_method == 2)) then
		write(*,*) " -> Writing polarization"
		call ftcrhd(unit, status)
		if (status > 0) then
			call print_error(status)
		endif
		naxis = 6
		naxes(1) = 3 !Q, U, V
		naxes(2)=Nlambda
		naxes(3)=npix_x
		naxes(4)=npix_y
		naxes(5)=RT_n_incl
		naxes(6)=RT_n_az
		nelements = naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5) * naxes(6)
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		call ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'Polarised Flux (Q, U, V)',status)
		call ftpprd(unit,group,fpixel,nelements,F_QUV,status)
		if (status > 0) then
			call print_error(status)
		endif
	end if
  
  ! create new hdu for wavelength grid
	call ftcrhd(unit, status)
  
	if (status > 0) then
		call print_error(status)
	endif
  
	if (lvacuum_to_air) then
		comment="AIR"
		lambda_vac = lambda
     	lambda = vacuum2air(Nlambda, lambda)
	end if 
   
	naxis = 1
	naxes(1) = Nlambda
	write(*,*) " (debug) writing lambda to image.."
	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
	call ftpkys(unit, "UNIT", "nm", comment, status)
	call ftpprd(unit,group,fpixel,Nlambda,lambda,status)
	if (status > 0) then
		call print_error(status)
	endif

  !  Close the file and free the unit number.
	call ftclos(unit, status)
	if (status > 0) then
		write(*,*) "error at closing"
		call print_error(status)
	endif
	call ftfiou(unit, status)
	if (status > 0) then
		write(*,*) "error at free unit"
		call print_error(status)
	endif

  !  Check for any error, and if so print out error messages
	if (status > 0) then
		call print_error(status)
	endif
	write(*,*) " (debug) done."


	return
	end subroutine write_image
 
	subroutine write_flux_ascii
	! -------------------------------------------------- !
	! written only for the first inclination / azimuth
	! --------------------------------------------------- !
  		integer :: status,unit
		integer :: la, j, i
		real(kind=dp), dimension(:), allocatable :: Fnu, Fnuc
  
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
		status=0
		unit = 10
		allocate(Fnu(Nlambda), Fnuc(Nlambda_cont))
		Fnu = 0.
		Fnuc = 0.
		if (RT_line_method==2) then
			do i=1, npix_x
				do j=1, npix_y
					Fnu = Fnu + Flux(:,i,j,1,1)
				enddo
			enddo
			do i=1, npix_x
				do j=1, npix_y
					Fnuc = Fnuc + Fluxc(:,i,j,1,1)
				enddo
			enddo
		else if (RT_line_method==1) then
		
			Fnu(:) = Flux(:,1,1,1,1)
			Fnuc(:) = Fluxc(:,1,1,1,1)
		else
			call error("RT_method unknown, write_flux_ascii")
		endif
   	
		open(unit,file="flux.s", status='unknown', iostat=status)
		do la=1, Nlambda
			write(unit, *) lambda(la), Fnu(la)
		enddo
		close(unit)
	
		open(unit,file="fluxc.s", status='unknown', iostat=status)
		do la=1, Nlambda_cont
			write(unit, *) lambda_cont(la), fnuc(la)
		enddo
		close(unit)

   
		deallocate(Fnu, Fnuc)

	return
	end subroutine write_flux_ascii
  
	subroutine writeWavelength()
	! --------------------------------------------------------------- !
	! Write wavelength grid build with each transition of each atom
	! --------------------------------------------------------------- !
		integer :: unit, EOF = 0, blocksize, naxes(1), naxis,group, bitpix, fpixel
		logical :: extend, simple
		character(len=6) :: comment="VACUUM"
		real(kind=dp) :: lambda_vac(Nlambda)
   
		if (.not.allocated(lambda).or..not.ltab_wavelength_image) return !
		write(*,*) " -> Writing wavelengths to ", WAVES_FILE       
   
		call ftgiou(unit,EOF)
		blocksize=1
		call ftinit(unit,trim(WAVES_FILE),blocksize,EOF)
		simple = .true.
		group = 1
		fpixel = 1
		extend = .false. 
		bitpix = -64   
		naxis = 1
		naxes(1) = Nlambda
   
		if (lvacuum_to_air) then
			comment="AIR"
			lambda_vac = lambda
			lambda = vacuum2air(Nlambda, lambda)
		end if 
   
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
		call ftpkys(unit, "UNIT", "nm", comment, EOF)
		call ftpprd(unit,group,fpixel,naxes(1),lambda,EOF)
		call ftclos(unit, EOF)
		call ftfiou(unit, EOF)
   
		if (EOF > 0) call print_error(EOF)       
 
	return
	end subroutine writeWavelength

	function vacuum2air(Nlambda, lambda_vac) result(lambda_air)
	!wavelength in nm
		integer, intent(in) :: Nlambda
		real(kind=dp), dimension(:), intent(in) :: lambda_vac
		real(kind=dp), dimension(Nlambda) :: lambda_air
		real(kind=dp), dimension(Nlambda) :: sqwave, reduction

		where (lambda_vac.ge.VACUUM_TO_AIR_LIMIT) 
     		sqwave = 1./(lambda_vac**2)
     		reduction = 1. + 2.735182e-4 + &
            (1.314182 + 2.76249e+4 * sqwave) * sqwave
     		lambda_air = lambda_vac / reduction
		else where(lambda_vac.lt.VACUUM_TO_AIR_LIMIT)
     		lambda_air = lambda_vac
		end where


	return
	end function vacuum2air

	function air2vacuum(Nlambda, lambda_air) result(lambda_vac)
	!wavelength in nm
		integer, intent(in) :: Nlambda
		real(kind=dp), dimension(:), intent(in) :: lambda_air
		real(kind=dp), dimension(Nlambda) :: lambda_vac
		real(kind=dp), dimension(Nlambda) :: sqwave, increase

		where (lambda_air.ge.AIR_TO_VACUUM_LIMIT)
			sqwave = (1.0e+7 / lambda_air)**2
			increase = 1.0000834213E+00 + &
            2.406030E+06/(1.30E+10 - sqwave) + &
            1.5997E+04/(3.89E+09 - sqwave)
			lambda_vac = lambda_air * increase
		else where(lambda_air.lt.AIR_TO_VACUUM_LIMIT)
			lambda_vac = lambda_air
		end where

	return
	end function air2vacuum
  
	subroutine WRITE_CNTRB_FUNC_PIX()
	! -------------------------------------------------- !
	! Write contribution function to disk.
	! --------------------------------------------------- !
		integer :: status,unit,blocksize,bitpix,naxis, naxis2
		integer, dimension(8) :: naxes, naxes2
		integer :: group,fpixel,nelements, nelements2
		logical :: simple, extend
		character(len=6) :: comment=""
		real :: pixel_scale_x, pixel_scale_y
  
		write(*,*)" -> writing contribution function..."
  
		!  Get an unused Logical Unit Number to use to open the FITS file.
		status=0
		call ftgiou (unit,status)

   !  Create the new empty FITS file.
		blocksize=1
		call ftinit(unit,TRIM(CF_FILE),blocksize,status)

		simple=.true.
		extend=.true.
		group=1
		fpixel=1

		bitpix=-64

		if (lVoronoi) then   
			naxis = 4
			naxes(1) = Nlambda
			naxes(2) = n_cells
			naxes(3) = RT_n_incl
			naxes(4) = RT_n_az
			nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
		else
			if (l3D) then
				naxis = 6
				naxes(1) = Nlambda
				naxes(2) = n_rad
				naxes(3) = 2*nz
				naxes(4) = n_az
				naxes(5) = RT_n_incl
				naxes(6) = RT_n_az
				nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4) * naxes(5) * naxes(6)
			else
				naxis = 5
				naxes(1) = Nlambda
				naxes(2) = n_rad
				naxes(3) = nz
				naxes(4) = RT_n_incl
				naxes(5) = RT_n_az
				nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4) * naxes(5)
			end if
		end if

  !  Write the required header keywords.
		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		call ftpkys(unit,'UNIT',"W.m-2.Hz-1.pixel-1",'Ksi',status)
  !  Write line CF to fits
		write(*,*) "Max,min abs(Ksi)=",maxval(dabs(Ksi)), minval(dabs(Ksi),mask=dabs(Ksi)>0)
		call ftpprd(unit,group,fpixel,nelements,ksi,status)


  !  Close the file and free the unit number.
		call ftclos(unit, status)
		call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
		if (status > 0) then
			call print_error(status)
		endif

	return
	end subroutine WRITE_CNTRB_FUNC_PIX

end module spectrum_type

