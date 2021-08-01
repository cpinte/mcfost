module spectrum_type

  use atom_type, only : AtomicLine, AtomicContinuum, AtomType
  use atmos_type, only : helium, hydrogen, NactiveAtoms, Natom, Atoms, v_char, icompute_atomRT, lmagnetized, laccretion_shock
  use getlambda, only  : hv, adjust_wavelength_grid, Read_wavelengths_table, make_wavelength_grid, make_wavelength_grid_new
  use fits_utils, only : print_error
  use parametres, only : n_cells, lelectron_scattering, n_etoiles, npix_x, npix_y, rt_n_incl, rt_n_az, &
       lorigin_atom, n_rad, n_az, nz, map_size, distance, zoom, lmagnetoaccr, lzeeman_polarisation, &
       lvacuum_to_air, ltab_wavelength_image, lvoronoi, l3D, mem_alloc_tot, llimit_mem
  use input, only : nb_proc, RT_line_method,lkeplerian, linfall, l_sym_ima
  use constantes, only : arcsec_to_deg
  use constant, only : clight
  use mcfost_env, only : dp
  use messages, only : error, warning
  !use hdf5

  implicit none

  real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
  real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352
  character(len=*), parameter :: WAVES_FILE="atom_transfer_waves_grid.fits.gz"
  ! not that FLUX_FILE is 1D only if one pixel, otherwise it is a map.
  ! F(x,y,0,lambda) in several directions.
  character(len=*), parameter :: FLUX_FILE="flux.fits.gz"
  character(len=*), parameter :: ORIGINC_FILE="origin_atom_cont.out", ORIGIN_FILE="origin_atom.out"
  character(len=*), parameter :: ORIGINC_FILE_FITS="origin_atom_cont.fits.gz", ORIGIN_FILE_FITS="origin_atom.fits.gz"

  !shift in index of line profiles, function of iray and nb_proc
  integer, dimension(:,:), allocatable :: dk
  integer :: dk_max, dk_min

  integer :: Nlambda, Ntrans, Nlambda_cont
  real(kind=dp) :: wavelength_ref = 0.0
  real(kind=dp), dimension(:,:), allocatable :: chi_c, sca_c, eta_c, chi, eta, chi_c_nlte, eta_c_nlte, chi0_bb, eta0_bb
  real(kind=dp), dimension(:,:,:), allocatable :: Icont
  real(kind=dp), dimension(:), allocatable :: lambda, lambda_cont
  real(kind=dp), dimension(:,:), allocatable :: Jnu_cont, Jnu
  real(kind=dp), dimension(:,:), allocatable :: Istar_tot, Istar_cont, Ishock, Istar_loc
  real(kind=dp), dimension(:,:,:), allocatable :: Itot
  real(kind=dp), allocatable, dimension(:,:,:,:) :: Flux !incl, az, lambda!, id ?
  real(kind=dp), allocatable, dimension(:,:,:,:) :: Fluxc
  real(kind=dp), allocatable, dimension(:,:,:,:) :: Flux_star, Flux_acc !total accretion flux (including lines if present)

  real(kind=dp), dimension(:,:,:), allocatable :: rho_p, chiQUV_p, etaQUV_p
  real(kind=dp), allocatable, dimension(:,:) :: Stokes_Q, Stokes_U, Stokes_V
  real(kind=dp), allocatable, dimension(:,:,:,:,:) :: F_QUV
  real(kind=dp), allocatable, dimension(:,:) :: origin_atom, originc_atom
  
  logical :: limage_at_lam0 = .false.
  real(kind=dp), allocatable, dimension(:,:,:,:) :: image_map

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


    if (present(lam0)) then
    	wavelength_ref = lam0
    	if (RT_line_method==2 .and. .not.lmagnetized) then
    		limage_at_lam0 = .true.
    	endif
    endif


    !for each line
    dk_max = int( sign(1.0_dp, v_char) * ( 1e-3 * abs(v_char) / hv + 0.5 ) )
    write(*,*) "Maximum shift in index:", dk_max, (1e-3 * v_char + hv) / hv
    if (1d3 * abs(dk_max) * hv / clight > 1.0) then
    	call error("Doppler shift larger than c!")
    endif

!     call make_wavelength_grid(wavelength_ref, v_char, lambda, Ntrans, lambda_cont)
    call make_wavelength_grid_new(wavelength_ref, dk_max, lambda, Ntrans, lambda_cont)
    dk_min = -dk_max

    mem_alloc_tot = mem_alloc_tot + sizeof(lambda) + sizeof(lambda_cont)

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
    integer :: nat, k, Nlambda_max, alloc_status, istar, size_phi_loc_tot
    type (AtomType), pointer :: atom
    logical, intent(in)    :: alloc_atom_nlte
    integer(kind=8) :: mem_alloc_local = 0

    if (allocated(Itot)) then
       write(*,*) "Error I already allocated"
       stop
    end if

    if (n_etoiles > 0) then
       allocate(Istar_tot(Nlambda,n_etoiles))
       Istar_tot(:,:) = 0.0_dp
       allocate(Istar_cont(Nlambda_cont,n_etoiles))
       Istar_cont(:,:) = 0.0_dp
       allocate(Istar_loc(Nlambda,nb_proc),stat=alloc_status)
       if (alloc_status > 0) call error("Allocation error Istar_loc")
       write(*,*) " -> using:", sizeof(Istar_loc)/1024./1024.," MB for local stellar radiation."
       mem_alloc_local = mem_alloc_local + sizeof(Istar_tot) + sizeof(Istar_cont) + sizeof(Istar_loc)
       Istar_loc = 0.0_dp
    endif


    allocate(Itot(Nlambda, Nray, nb_proc),stat=alloc_status)
    if (alloc_status > 0) call error("Allocation error Itot")
    allocate(Icont(Nlambda_cont, Nray, nb_proc))
    Itot = 0.0_dp
    Icont = 0.0_dp

    mem_alloc_local = mem_alloc_local + sizeof(Itot)+sizeof(Icont)


    if (laccretion_shock) then
       allocate(Ishock(Nlambda,nb_proc), stat=alloc_status)
       if (alloc_Status > 0) call ERROR ("Cannot allocate Ishock")
       write(*,*) " -> using:", sizeof(Ishock)/1024./1024.," MB for accretion intensity!"
       mem_alloc_local = mem_alloc_local + sizeof(Ishock)
       Ishock = 0.0_dp
    endif

    !allocate(dk(Nray, nb_proc))
    !if (alloc_status > 0) call error("Allocation error dk")

    allocate(chi_c(Nlambda_cont,n_cells), eta_c(Nlambda_cont,n_cells), stat=alloc_status)
    !-> At the moment not needed because only Thomson scattering included
    !allocate(sca_c(Nlambda_cont,n_cells), stat=alloc_status)
    ! 		write(*,*) " size chi_c/eta_c ", 2 * n_cells * Nlambda_cont/1024./1024./1024., " GB"
    ! 		write(*,*) " size chi_c/eta_c ", (size(chi_c)+size(eta_c))/1024./1024./1024., " GB"

    if (alloc_status > 0) then
       call error("Allocation error, continuum opacities")
    endif
    chi_c = 0.0_dp
    eta_c = 0.0_dp
    if (allocated(sca_c)) then
       sca_c = 0.0_dp
       mem_alloc_local = mem_alloc_local + sizeof(chi_c)
       write(*,*) " -> size contopac:", 3 * sizeof(chi_c) /1024./1024./1024.," GB"
    else
       write(*,*) " -> size contopac:", 2 * sizeof(chi_c) /1024./1024./1024.," GB"
    endif

    mem_alloc_local = mem_alloc_local + sizeof(chi_c) * 2

    allocate(eta(Nlambda ,nb_proc))
    allocate(chi(Nlambda ,nb_proc))

    eta(:,:) = 0.0_dp
    chi(:,:) = 0.0_dp

    !interpolated total continuum opacities on the lambda grid to be used with lines
    if (.not.llimit_mem) then
       allocate(eta0_bb(Nlambda , n_cells))
       allocate(chi0_bb(Nlambda , n_cells))
       eta0_bb = 0.0_dp
       chi0_bb = 0.0_dp
       write(*,*) " -> size contopac (line grid):", 2*sizeof(chi0_bb)/1024./1024./1024.," GB"
       mem_alloc_local = mem_alloc_local + sizeof(eta0_bb) * 2
    else
       write(*,*) " -> Reducing memory usage by interpolating continuous opacity!"
    endif

    !otherwise allocated bellow for nlte.
    if ((lelectron_scattering).and.(.not.alloc_atom_nlte)) then
       call alloc_jnu
       mem_alloc_local = mem_alloc_local +  sizeof(Jnu_cont) + sizeof(Jnu)
       write(*,*) " -> size Jnu(lambda):", sizeof(Jnu)/1024./1024.," MB"
       write(*,*) " -> size Jnu_cont:", sizeof(Jnu_cont)/1024./1024./1024.," GB"
    endif


    if (alloc_atom_nlte) then !NLTE loop activated
       allocate(chi_c_nlte(Nlambda_cont,n_cells),stat=alloc_status)
       if (alloc_status >0) call error("Allocation error chi_c_nlte")
       allocate(eta_c_nlte(Nlambda_cont,n_cells),stat=alloc_status)
       if (alloc_status >0) call error("Allocation error eta_c_nlte")
       chi_c_nlte = 0.0_dp
       eta_c_nlte = 0.0_dp
       mem_alloc_local = mem_alloc_local +  2 * sizeof(chi_c_nlte)

       if (lelectron_scattering) then
          call alloc_jnu
          mem_alloc_local = mem_alloc_local +  sizeof(Jnu_cont) + sizeof(Jnu)
          write(*,*) " -> size Jnu(lambda):", sizeof(Jnu)/1024./1024.," MB"
          write(*,*) " -> size Jnu_cont:", sizeof(Jnu_cont)/1024./1024./1024.," GB"
       endif

       ! .....  add other NLTE opac
       size_phi_loc_tot = 0
       do nat=1, Natom
          atom => Atoms(nat)%ptr_atom
          do k=1, atom%Nline
             allocate(atom%lines(k)%phi_loc(atom%lines(k)%Nred-dk_min+dk_max+1-atom%lines(k)%Nblue,Nray,nb_proc),stat=alloc_status)
             if (alloc_status > 0) call error("Allocation error line%phi_loc for nlte loop")
             size_phi_loc_tot = size_phi_loc_tot + sizeof(atom%lines(k)%phi_loc)
          enddo
       enddo
       mem_alloc_local = mem_alloc_local + size_phi_loc_tot

    endif

    mem_alloc_tot = mem_alloc_tot + mem_alloc_local
    write(*,'("Total memory allocated in alloc_spectrum:"(1F14.3)" GB")') mem_alloc_local / 1024./1024./1024.

    return
  end subroutine alloc_Spectrum


  subroutine alloc_flux_image
    !Store total flux and flux maps for selected lines
    integer :: alloc_status, kr, nat
    real :: mem_alloc, mem_flux, mem_cont, mem_for_file
    integer(kind=8) :: mem_alloc_local = 0,  Ntrans_tot
    real, dimension(:), allocatable :: mem_per_file
    integer :: Nlam, Ntrans

    allocate(mem_per_file(Natom))
    mem_per_file = 0.0

    Ntrans_tot = 0
    do nat=1, Natom

       Ntrans = 0
       Nlam = 0

       do kr=1,atoms(nat)%ptr_atom%Nline

          if (atoms(nat)%ptr_atom%lines(kr)%write_flux_map) then

             Ntrans = Ntrans + 1 !count number of transition for that atom
             Nlam = Nlam + atoms(nat)%ptr_atom%lines(kr)%Nred-dk_min+dk_max-atoms(nat)%ptr_atom%lines(kr)%Nblue+1


          endif


       enddo
       if (Ntrans == 0) cycle !go to next atom
       Ntrans_tot = Ntrans_tot + Ntrans
       !mean memory needed for that atom, Nlam+1 for the continuum, +Nlam for the grid
       mem_per_file(nat) = real(8 * (Nlam + 1) + 8 * Nlam) * real(rt_n_incl*rt_n_az) &
            * real(npix_x)/1024. * real(npix_y)/1024. / real(Ntrans)
    enddo
    mem_per_file = mem_per_file!total for each atom, average of all transitions.
    mem_for_file = sum(mem_per_file) / real(Natom) !average per atom in MB

    !memory allocated for each file, so for all transitions of each atom (the one we write)
    if (mem_for_file/1024 > 2.5) then
       call warning("Size of fits file to store flux map for each line might be large")
    else
       write(*,*) " Will write in average ", mem_for_file, " MB for each line for each atom!"
    endif

    !total flux and total cont  +  their wavelength grid
    mem_cont = Nlambda_cont + Nlambda_cont*rt_n_incl*rt_n_az*nb_proc
    mem_flux = Nlambda + Nlambda*rt_n_incl*rt_n_az*nb_proc
    if (lmagnetized) mem_flux = mem_flux + 3 * Nlambda*rt_n_incl*rt_n_az


    mem_flux = 8 * mem_flux / 1024. / 1024.
    mem_cont = 8 * mem_cont / 1024. / 1024.
    mem_alloc = mem_for_file * real(Natom) !for all transitions of all atoms !, the one we keep in memory

    !Flux  + Flux cont + lines grid
    if (mem_flux + mem_cont > 1d3) then !in MB
       write(*,*) " allocating ",( mem_flux + mem_cont )/1024., " GB for total flux arrays.."
    else
       write(*,*) " allocating ", mem_flux + mem_cont, " MB for total flux arrays.."
    endif

    !otherwise there is not flux map stored
    if (Ntrans_tot > 0) then
       !For all lines of all atoms!
       if (mem_alloc > 1d3) then
          write(*,*) " allocating ",mem_alloc / 1024., " GB of mem for flux map for selected lines!"
       else
          write(*,*) " allocating ",mem_alloc, " MB of mem for flux map for selected lines!"
       endif

    endif

    write(*,*) "  -> ", mem_alloc+mem_flux+mem_cont, " MB in total"

    !remove the pixel dimension
    allocate(Flux(Nlambda,RT_N_INCL,RT_N_AZ,nb_proc), stat=alloc_status)
    ! 		allocate(Flux(Nlambda,RT_N_INCL,RT_N_AZ), stat=alloc_status)
    if (alloc_Status > 0) call ERROR ("Cannot allocate Flux")
    allocate(Fluxc(Nlambda_cont,RT_N_INCL,RT_N_AZ,nb_proc), stat=alloc_status)
    ! 		allocate(Fluxc(Nlambda_cont,RT_N_INCL,RT_N_AZ), stat=alloc_status)
    if (alloc_Status > 0) call ERROR ("Cannot allocate Flux continuum")
    Flux = 0.0_dp
    Fluxc = 0.0_dp

    mem_alloc_local = mem_alloc_local + sizeof(Flux)+sizeof(fluxc)
    if (lmagnetized) mem_alloc_local = mem_alloc_local + 3*sizeof(flux)/nb_proc

	if (n_etoiles > 0) then
    	allocate(flux_star(Nlambda,RT_N_INCL,RT_N_AZ,nb_proc), stat=alloc_status)
    	if (alloc_Status > 0) call ERROR ("Cannot allocate Flux_star")
    	mem_alloc_local = mem_alloc_local + sizeof(Flux_star)
    	Flux_star = 0.0_dp
    endif

    if (laccretion_shock) then
       allocate(Flux_acc(Nlambda,RT_N_INCL,RT_N_AZ,nb_proc), stat=alloc_status)
       if (alloc_Status > 0) call ERROR ("Cannot allocate Flux_acc")
       write(*,*) " -> using ", sizeof(flux_acc)/1024./1024.," MB for accretion flux!"
       mem_alloc_local = mem_alloc_local + sizeof(flux_acc)
       Flux_acc = 0.0_dp
    endif

    !now for the lines
    do nat=1, Natom

       do kr=1,atoms(nat)%ptr_atom%Nline

          if (atoms(nat)%ptr_atom%lines(kr)%write_flux_map) then

             Nlam = (atoms(nat)%ptr_atom%lines(kr)%Nred-dk_min+dk_max- &
                  atoms(nat)%ptr_atom%lines(kr)%Nblue+1)

             allocate(atoms(nat)%ptr_atom%lines(kr)%map(Nlam, npix_x, npix_y, rt_n_incl, rt_n_az),stat=alloc_status)
             if (alloc_status > 0) then
                write(*,*) atoms(nat)%ptr_atom%ID,atoms(nat)%ptr_atom%lines(kr)%j, atoms(nat)%ptr_atom%lines(kr)%i
                write(*,*) atoms(nat)%ptr_atom%g(atoms(nat)%ptr_atom%lines(kr)%j), &
                     atoms(nat)%ptr_atom%g(atoms(nat)%ptr_atom%lines(kr)%i)
                call error("Cannot allocate map for this line !")
             endif
             allocate(atoms(nat)%ptr_atom%lines(kr)%mapc(npix_x, npix_y, rt_n_incl, rt_n_az),stat=alloc_status)
             if (alloc_status > 0) then
                write(*,*) atoms(nat)%ptr_atom%ID,atoms(nat)%ptr_atom%lines(kr)%j, atoms(nat)%ptr_atom%lines(kr)%i
                write(*,*) atoms(nat)%ptr_atom%g(atoms(nat)%ptr_atom%lines(kr)%j), &
                     atoms(nat)%ptr_atom%g(atoms(nat)%ptr_atom%lines(kr)%i)
                call error("Cannot allocate mapc for this line !")
             endif
             atoms(nat)%ptr_atom%lines(kr)%map = 0.0_dp
             atoms(nat)%ptr_atom%lines(kr)%mapc = 0.0_dp
             mem_alloc_local = mem_alloc_local + sizeof(atoms(nat)%ptr_atom%lines(kr)%map) + &
                  sizeof(atoms(nat)%ptr_atom%lines(kr)%mapc)
          endif


       enddo

    enddo

    !Contribution functions (for one ray it is allocated elsewhere)
    !Future: contribution function for selected lines only !
    ! 		if (lcontrib_function) then
    !
    ! 			mem_alloc = real(8 * n_cells * Nlambda) / real(1024*1024)
    !
    ! 			if (mem_alloc > 1000.0) then
    ! 				write(*,*) " allocating ", mem_alloc/real(1024), " GB for contribution function.."
    ! 			else
    ! 				write(*,*) " allocating ", mem_alloc, " MB for contribution function.."
    ! 			endif
    !
    ! 			if (mem_alloc >= 2.1d3) then !2.1 GB
    ! 				call Warning(" To large cntrb array. Use a wavelength table instead..")
    ! 				lcontrib_function = .false.
    ! 			else
    !
    ! 				allocate(cntrb(Nlambda,n_cells),stat=alloc_status)
    ! 				mem_alloc_local = mem_alloc_local + sizeof(cntrb)
    ! 				if (alloc_status > 0) then
    ! 					call ERROR('Cannot allocate cntrb_ray')
    ! 					lcontrib_function = .false.
    ! 				else
    ! 					cntrb(:,:) = 0.0_dp
    ! 				endif
    !
    ! 			end if
    !
    ! 		end if

    if (lorigin_atom) then
       mem_alloc = 8 * n_cells * Nlambda / 1024./ 1024. + 8 * n_cells * Nlambda_cont / 1024./1024.
       if (mem_alloc > 1d3) then
          write(*,*) " allocating ", mem_alloc/1024., " GB for local emission origin.."
       else
          write(*,*) " allocating ", mem_alloc, " MB for local emission origin.."
       endif

       if (mem_alloc >= 2.1d3) then !2.1 GB
          call Warning(" To large origin_atom array. Use a wavelength table instead..")
          !lorigin_atom = .false.
       else

          allocate(origin_atom(Nlambda,n_cells),stat=alloc_status)
          allocate(originc_atom(Nlambda_cont,n_cells),stat=alloc_status)

          mem_alloc_local = mem_alloc_local + sizeof(origin_atom) + sizeof(originc_atom)
          if (alloc_status > 0) then
             call ERROR('Cannot allocate origin_atom')
             !lorigin_atom = .false.
          else
             origin_atom = 0.0_dp
             originc_atom = 0.0_dp
          endif

       end if
    endif

	if ( limage_at_lam0 ) then

    	allocate(image_map(npix_x, npix_y, rt_n_incl, rt_n_az),stat=alloc_status)
        if (alloc_status > 0) then
            call error("Cannot allocate image_map !")
        endif	
        image_map = 0.0_dp
        mem_alloc_local = mem_alloc_local + sizeof(image_map)
	
	endif

    mem_alloc_tot = mem_alloc_tot + mem_alloc_local
    write(*,'("Total memory allocated in alloc_flux_image:"(1F14.3)" GB")') mem_alloc_local / 1024./1024./1024.


    deallocate(mem_per_file)

    return
  end subroutine alloc_flux_image

  subroutine allocate_stokes_quantities
    ! only available for flux calculations
    ! this routine is simplified there is only on solution, full_stokes!

    allocate(Stokes_Q(Nlambda, nb_proc), Stokes_U(Nlambda, nb_proc), Stokes_V(Nlambda, nb_proc))
    Stokes_Q = 0.0_dp
    Stokes_U = 0.0_dp
    Stokes_V = 0.0_dp

    allocate(F_QUV(Nlambda,RT_N_INCL,RT_N_AZ,3,nb_proc)); F_QUV = 0.0_dp
    allocate(rho_p(Nlambda, 3, nb_proc)); rho_p = 0.0_dp
    allocate(etaQUV_p(Nlambda, 3, nb_proc)); etaQUV_p = 0.0_dp
    allocate(chiQUV_p(Nlambda, 3, nb_proc)); chiQUV_p = 0.0_dp

    return
  end subroutine allocate_stokes_quantities

  ! 	subroutine fill_map(ibin,iaz,ipix,jpix,method, I0)
  ! 		real(kind=dp), intent(in) :: I0(Nlambda)!It's actually I0 * normF
  ! 		integer, intent(in) :: ibin, iaz, ipix, jpix, method
  ! 		integer :: nat, kr, nr, nb
  ! 		type (AtomType), pointer :: atom
  !
  ! 		if (method==1) then
  ! 		!summation over pixels
  ! 			do nat=1,Natom
  ! 				atom => atoms(nat)%ptr_atom
  ! 				do kr=1,atom%Nline
  !
  ! 					if (atom%lines(kr)%write_flux_map) then
  ! 						nr = atom%lines(kr)%Nred + dk_max
  ! 						nb = atom%lines(kr)%Nblue + dk_min
  !
  !   						atom%lines(kr)%map(:,1,1,ibin,iaz) = &
  !   						atom%lines(kr)%map(:,1,1,ibin,iaz) + I0(nb:nr)
  !
  !   					endif
  !
  !   				enddo
  !   				atom => NULL()
  !   			enddo
  ! 		else!2D maps
  ! 			do nat=1,Natom
  !   				atom => atoms(nat)%ptr_atom
  !   				do kr=1,atom%Nline
  !
  !   					if (atom%lines(kr)%write_flux_map) then
  !   						nr = atom%lines(kr)%Nred + dk_max
  !   						nb = atom%lines(kr)%Nblue + dk_min
  !
  !   						atom%lines(kr)%map(:,ipix,jpix,ibin,iaz) = I0(nb:nr)
  !
  !   					endif
  !
  !   				enddo
  !   				atom => NULL()
  !   			enddo
  !
  !   		endif
  !
  ! 	return
  ! 	end subroutine fill_map

  subroutine dealloc_spectrum()
  
  
    call dealloc_jnu

    deallocate(lambda)
    if (allocated(lambda_cont)) deallocate(lambda_cont)

    deallocate(Itot, Icont)
    if (allocated(Istar_tot)) deallocate(Istar_tot, Istar_cont, Istar_loc, flux_star)
    !can be deallocated before to save memory
    if (allocated(Flux)) deallocate(Flux)
    if (allocated(Fluxc)) deallocate(Fluxc)
    if (allocated(flux_acc)) deallocate(flux_acc,Ishock)

    if (allocated(dk)) deallocate(dk)

    !check allocation due to field_free sol
    if (allocated(Stokes_Q)) then !same for all dangerous
       deallocate(Stokes_Q, Stokes_U, Stokes_V, F_QUV)
       deallocate(rho_p)
       deallocate(etaQUV_p)
       deallocate(chiQUV_p)
    endif


    deallocate(chi_c,  eta_c)
    if (allocated(sca_c)) deallocate(sca_c)
    deallocate(chi, eta)

    if (allocated(chi0_bb)) then
       deallocate(chi0_bb, eta0_bb)
    endif


    if (allocated(origin_atom)) deallocate(origin_atom, originc_atom)

	if (limage_at_lam0) deallocate(image_map)

    return
  end subroutine dealloc_Spectrum


  subroutine alloc_Jnu()


	!Flat and locally interpolated !
    allocate(Jnu(Nlambda,1))
    Jnu = 0.0

	!At the moment stored on the whole grid !
    allocate(Jnu_cont(Nlambda_cont,n_cells))
    Jnu_cont = 0.0


    return
  end subroutine alloc_Jnu

  subroutine dealloc_Jnu()

    if (allocated(Jnu_cont)) deallocate(Jnu_cont)
    if (allocated(Jnu)) deallocate(Jnu)

    return
  end subroutine dealloc_Jnu


  subroutine write_flux(only_flux)
    !
    !write flux total and flux map for lines
    !
    logical, optional :: only_flux
    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(6) :: naxes
    integer :: group,fpixel,nelements
    integer :: la, Nred, Nblue, kr, m, n, i, j
    logical :: simple, extend
    character(len=6) :: comment="VACUUM"
    real(kind=dp) :: lambda_vac(Nlambda), Fnu, lambdac_vac(Nlambda_cont)
    real :: pixel_scale_x, pixel_scale_y
    type (AtomType), pointer :: atom
    logical, dimension(Natom) :: write_map
    real(kind=dp), allocatable :: g_i(:), g_j(:), maps(:,:,:,:,:,:), nu_trans(:), lam_trans(:), waves(:,:)
    integer, allocatable :: i_level(:), j_level(:)
    integer :: Ntrans, Nlam

    write(*,*) "Writing Flux"

    blocksize=1
    simple=.true.
    extend=.false.
    group=1
    fpixel=1
    bitpix=-64

    !  Get an unused Logical Unit Number to use to open the FITS file.
    status=0
    call ftgiou (unit,status)

    !  Create the new empty FITS file.

    call ftinit(unit,trim(FLUX_FILE),blocksize,status)

    naxis = 3
    naxes(1) = Nlambda
    naxes(2) = RT_n_incl
    naxes(3) = RT_n_az

    nelements = naxes(1)*naxes(2)*naxes(3)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    if (status > 0) then
       call print_error(status)
    endif
    call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm F_{\nu}}$',status)


    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,sum(Flux(:,:,:,:),dim=4),status)
    ! 	call ftpprd(unit,group,fpixel,nelements,Flux(:,:,:),status)
    if (status > 0) then
       call print_error(status)
    endif


    ! create new hdu for wavelength grid
    call ftcrhd(unit, status)

    if (status > 0) then
       call print_error(status)
    endif


    if (lvacuum_to_air) then
       comment="AIR"
       lambda_vac = lambda
       lambda = vacuum2air(Nlambda, lambda_vac)
    end if
    naxis = 1
    naxes(1) = Nlambda
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    call ftpkys(unit, "UNIT", "nm", comment, status)
    call ftpprd(unit,group,fpixel,Nlambda,lambda,status)
    if (status > 0) then
       call print_error(status)
    endif

    !now continuum flux
    call ftcrhd(unit, status)

    if (status > 0) then
       call print_error(status)
    endif

    naxis = 3
    naxes(1) = Nlambda_cont
    naxes(2) = RT_n_incl
    naxes(3) = RT_n_az

    nelements = naxes(1)*naxes(2)*naxes(3)

    !  Write the required header keywords.
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    if (status > 0) then
       call print_error(status)
    endif
    call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm F_{\nu}^{cont}}$',status)


    !  Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,sum(Fluxc(:,:,:,:),dim=4),status)
    ! 	call ftpprd(unit,group,fpixel,nelements,Fluxc(:,:,:),status)
    if (status > 0) then
       call print_error(status)
    endif


    ! create new hdu for cont wavelength grid
    call ftcrhd(unit, status)

    if (status > 0) then
       call print_error(status)
    endif


    if (lvacuum_to_air) then
       comment="AIR"
       lambdac_vac = lambda_cont
       lambda_cont = vacuum2air(Nlambda_cont, lambdac_vac)
    end if

    naxis = 1
    naxes(1) = Nlambda_cont
    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    call ftpkys(unit, "UNIT", "nm", comment, status)
    call ftpprd(unit,group,fpixel,Nlambda_cont,lambda_cont,status)
    if (status > 0) then
       call print_error(status)
    endif

    !Stellar flux
    if (n_etoiles > 0) then
    	call ftcrhd(unit, status)
    	if (status > 0) then
       		call print_error(status)
    	endif
    	naxis = 3
    	naxes(1) = Nlambda
    	naxes(2) = RT_n_incl
    	naxes(3) = RT_n_az

    	nelements = naxes(1)*naxes(2)*naxes(3)

    	!  Write the required header keywords.
    	call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
		if (status > 0) then
       		call print_error(status)
    	endif
    	call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm Fstar_{\nu}}$',status)
	

    	!  Write the array to the FITS file.
    	if (maxval(sum(Flux_star(:,:,:,:),dim=4)) <= 0.0) then
    		call warning("Bug with stellar flux !")
    		write(*,*) " Flux_star_max = ", maxval(sum(Flux_star(:,:,:,:),dim=4))
    	endif
    	call ftpprd(unit,group,fpixel,nelements,sum(Flux_star(:,:,:,:),dim=4),status)
    	if (status > 0) then
       		call print_error(status)
    	endif
    endif

    if (lmagnetized) then

       call ftcrhd(unit, status)
       if (status > 0) then
          call print_error(status)
       endif

       naxis = 4
       naxes(1) = Nlambda
       naxes(2) = RT_n_incl
       naxes(3) = RT_n_az
       naxes(4) = 3 !Q,U,V

       nelements = naxes(1)*naxes(2)*naxes(3)*naxes(4)

       !  Write the required header keywords.
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) then
          call print_error(status)
       endif
       call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm X_{\nu}}$',status)

       !  Write the array to the FITS file.
       ! 		call ftpprd(unit,group,fpixel,nelements,F_QUV,status)
       call ftpprd(unit,group,fpixel,nelements,sum(F_QUV,dim=5),status)
       if (status > 0) then
          call print_error(status)
       endif

    endif

    if (laccretion_shock) then
       call ftcrhd(unit, status)
       if (status > 0) then
          call print_error(status)
       endif
       write(*,*) " Writing accretion flux to fits file"
       naxis = 3
       naxes(1) = Nlambda
       naxes(2) = RT_n_incl
       naxes(3) = RT_n_az
       nelements = naxes(1)*naxes(2)*naxes(3)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) then
          call print_error(status)
       endif
       call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm Facc_{\nu}}$',status)
       if (maxval(sum(flux_acc(:,:,:,:),dim=4)) <= 0.0) then
       	call warning("Bug with accretion flux!")
       	write(*,*)  "   -> max acc flux :", maxval(sum(flux_acc(:,:,:,:),dim=4))
       endif
       call ftpprd(unit,group,fpixel,nelements,sum(Flux_acc(:,:,:,:),dim=4),status)
       if (status > 0) then
          call print_error(status)
       endif
    endif

    ! Close the file and free the unit number.
    call ftclos(unit, status)
    if (status > 0) then
       call print_error(status)
    endif

    call ftfiou(unit, status)
    if (status > 0) then
       call print_error(status)
    endif

    if (present(only_flux)) then
       if (only_flux) return
    endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    write(*,*) "You should use the new map!"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    blocksize=1
    simple=.true.
    extend=.true.
    group=1
    fpixel=1
    bitpix=-64

    !do we need to write line map ? for what atoms
    write_map(:) = .false.
    atm_loop : do n=1, Natom
       atom => atoms(n)%ptr_atom

       line_loop : do kr=1,atom%nline
          if (atom%lines(kr)%write_flux_map) then

             write_map(n) = .true.
             exit line_loop

          endif
       enddo line_loop

       atom => NULL()
    enddo atm_loop

    !Now write the individual files for each atom lines map
    !Continuum map not written yet
    !write the resolution of grid
    do n=1, Natom
       atom => atoms(n)%ptr_atom
       if (.not.write_map(n)) cycle

       !count transition and max frequency for this atom
       Nlam = 0
       Ntrans = 0
       do kr=1, atom%Nline
          if (atom%lines(kr)%write_flux_map) then
             Ntrans = Ntrans + 1
             Nlam = max(Nlam, size(atom%lines(kr)%map(:,1,1,1,1)))
          endif
       enddo
       !allocate temporary variables for storage
       allocate(i_level(Ntrans), j_level(Ntrans), g_i(Ntrans), g_j(Ntrans), maps(Nlam,npix_x,npix_y,rt_n_incl,rt_n_az,Ntrans), &
            nu_trans(Ntrans), lam_trans(Ntrans), waves(Nlam,Ntrans))

       waves = 0.0_dp
       maps = 0.0_dp
       j = 0
       do kr=1, atom%Nline
          if (atom%lines(kr)%write_flux_map) then
             j = j + 1
             i_level(j) = atom%lines(kr)%i
             j_level(j) = atom%lines(kr)%j
             g_i(j) = atom%g(i_level(j))
             g_j(j) = atom%g(j_level(j))
             lam_trans(j) = atom%lines(kr)%lambda0
             nu_trans(j) = 1d9 * clight / lam_trans(j)
             !includes the interaction zone due to velocity fields.
             Nred = atom%lines(kr)%Nred+dk_max
             Nblue = atom%lines(kr)%Nblue+dk_min
             waves(1:Nred-Nblue+1,j) = lambda(Nblue:Nred)
             maps(1:Nred-Nblue+1,:,:,:,:,j) = atom%lines(kr)%map(:,:,:,:,:)
          endif
       enddo

       !Get an unused Logical Unit Number to use to open the FITS file.
       status=0
       call ftgiou (unit,status)

       call ftinit(unit,trim(atom%ID)//"_lines.fits.gz",blocksize,status)

       naxis=6
       naxes(6)=Ntrans
       naxes(2)=npix_x
       naxes(3)=npix_y
       naxes(4)=RT_n_incl
       naxes(5)=RT_n_az

       naxes(1)=Nlam
       nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)

       !if we enter here there is necessarily at least one map!
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
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
       call ftpkye(unit,'dv',hv,-7,'km.s^-1',status)

       !Write the array to the FITS file.
       call ftpprd(unit,group,fpixel,nelements,maps,status)
       if (status > 0) call print_error(status)

       !create new hdu for cont wavelength grid
       naxis=2
       naxes(2)=Ntrans
       naxes(1)=Nlam
       nelements = naxes(1)*naxes(2)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpkys(unit, "UNIT", "nm", comment, status)
       call ftpprd(unit,group,fpixel,nelements,waves,status)


       !create new hdu for nu and lam
       naxis=1
       naxes(1)=Ntrans
       nelements = naxes(1)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpprd(unit,group,fpixel,nelements,nu_trans,status)
       naxis=1
       naxes(1)=Ntrans
       nelements = naxes(1)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpprd(unit,group,fpixel,nelements,lam_trans,status)

       !create new hdu for i and j
       naxis=1
       naxes(1)=Ntrans
       nelements = naxes(1)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,8,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpprj(unit,group,fpixel,nelements,i_level,status)
       naxis=1
       naxes(1)=Ntrans
       nelements = naxes(1)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,8,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpprj(unit,group,fpixel,nelements,j_level,status)

       !create new hdu for gi and gj
       naxis=1
       naxes(1)=Ntrans
       nelements = naxes(1)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpprd(unit,group,fpixel,nelements,g_i,status)
       naxis=1
       naxes(1)=Ntrans
       nelements = naxes(1)
       call ftcrhd(unit, status)
       if (status > 0) call print_error(status)
       call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
       if (status > 0) call print_error(status)
       call ftpprd(unit,group,fpixel,nelements,g_j,status)

       deallocate(i_level, j_level, g_i, g_j, maps, nu_trans, lam_trans, waves)

       atom => NULL()
       !  Close the file and free the unit number.
       call ftclos(unit, status)
       if (status > 0) then
          call print_error(status)
       endif
       call ftfiou(unit, status)
       if (status > 0) then
          call print_error(status)
       endif
    enddo ! over atoms


    return
  end subroutine write_flux

  !lambda air-> vacuum to do.
  subroutine write_atomic_maps()
    !
    !write 2d flux maps for atomic lines in individual files
    !
    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(5) :: naxes
    integer :: group,fpixel,nelements
    integer :: kr,n, j
    logical :: simple, extend
    character(len=10) :: transition
    character(len=6) :: comment="VACUUM"
    real(kind=dp), dimension(:), allocatable :: lambda_vac(:), Fnu, lambdac_vac(:)
    real :: pixel_scale_x, pixel_scale_y
    type (AtomType), pointer :: atom
    logical, dimension(Natom) :: write_map

    write(*,*) "Writing line Flux maps"

    if (lvacuum_to_air) then
       comment="AIR"
       write(*,*) " -> lines' wavelengths converted to air in the call of write_flux() !"
    endif

    blocksize=1
    simple=.true.
    extend=.false.
    group=1
    fpixel=1
    bitpix=-64

    naxis=5
    naxes(2)=npix_x
    naxes(3)=npix_y
    naxes(4)=RT_n_incl
    naxes(5)=RT_n_az
    nelements=naxes(2)*naxes(3)*naxes(4)*naxes(5)

    !do we need to write line map ? for what atoms
    write_map(:) = .false.
    atm_loop : do n=1, Natom
       atom => atoms(n)%ptr_atom

       line_loop : do kr=1,atom%nline
          if (atom%lines(kr)%write_flux_map) then

             write_map(n) = .true.
             exit line_loop

          endif
       enddo line_loop

       atom => NULL()
    enddo atm_loop

    !Now write the individual files for each atom lines map
    !Continuum map not written yet
    !write the resolution of grid
    do n=1, Natom
       atom => atoms(n)%ptr_atom
       if (.not.write_map(n)) cycle

       j = 0
       do kr=1, atom%Nline
          if (atom%lines(kr)%write_flux_map) then
             j = j + 1
             !Get an unused Logical Unit Number to use to open the FITS file.
             status=0
             call ftgiou (unit,status)
             write(transition, '(1I10)') j !line j that we follow among all lines.
             call ftinit(unit,trim(atom%ID)//"_line_"//trim(adjustl(transition))//".fits.gz",blocksize,status)

             !-> includes the wavelengths due to velocity
             naxes(1) = atom%lines(kr)%Nred+dk_max - atom%lines(kr)%Nblue-dk_min +1

             !if we enter here there is necessarily at least one map!
             call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
             if (status > 0) call print_error(status)
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

             call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm F_{\nu}}$',status)
             call ftpkye(unit,'dv',hv,-7,'${\rm km\,s^{-1}}$',status)
             !add dk_max and dk_min on Nblue and Nred ?
             call ftpkyj(unit,'dvshift',dk_max,'pixels',status)
             !of the local profile with v = 0
             call ftpkyj(unit,'ib',atom%lines(kr)%Nblue,'index blue',status)
             call ftpkyj(unit,'ir',atom%lines(kr)%Nred,'index red',status)
             call ftpkyj(unit,'i',atom%lines(kr)%i,'',status)
             call ftpkyj(unit,'j',atom%lines(kr)%j,'',status)
             call ftpkyd(unit,'gi',atom%g(atom%lines(kr)%i),-7,'',status)
             call ftpkyd(unit,'gj',atom%g(atom%lines(kr)%j),-7,'',status)
             call ftpkyd(unit,'lambda0',atom%lines(kr)%lambda0,-7,'nm',status)
             call ftpkyd(unit,'nu0',1d-6*clight / atom%lines(kr)%lambda0,-7,'10^15 Hz',status)

             !Write the array to the FITS file.
             call ftpprd(unit,group,fpixel,naxes(1)*nelements,atom%lines(kr)%map(:,:,:,:,:),status)
             if (status > 0) call print_error(status)

             !create new hdu for lambda grid
             !-> convert to air ?
             call ftcrhd(unit, status)
             if (status > 0) call print_error(status)
             call ftphpr(unit,simple,bitpix,1,naxes,0,1,extend,status)
             if (status > 0) call print_error(status)
             call ftpkys(unit, "UNIT", "nm", comment, status)
             call ftpprd(unit,group,fpixel,naxes(1),lambda(atom%lines(kr)%Nblue+dk_min:atom%lines(kr)%Nred+dk_max),status)

             !-> now continuum image at central wavelength

             call ftcrhd(unit, status)
             if (status > 0) call print_error(status)
             call ftphpr(unit,simple,bitpix,4,naxes(2:naxis),0,1,extend,status)
             if (status > 0) call print_error(status)
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

             call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm F_{\nu}}$',status)
             call ftpkyj(unit,'index', atom%lines(kr)%Nmid,'index',status)
             call ftpkyj(unit,'i',atom%lines(kr)%i,'',status)
             call ftpkyj(unit,'j',atom%lines(kr)%j,'',status)
             call ftpkyd(unit,'gi',atom%g(atom%lines(kr)%i),-7,'',status)
             call ftpkyd(unit,'gj',atom%g(atom%lines(kr)%j),-7,'',status)
             call ftpkyd(unit,'lambda0',atom%lines(kr)%lambda0,-7,'nm',status)
             call ftpkyd(unit,'nu0',1d-6*clight / atom%lines(kr)%lambda0,-7,'10^15 Hz',status)

             !Write the array to the FITS file.
             call ftpprd(unit,group,fpixel,nelements,atom%lines(kr)%mapc(:,:,:,:),status)
             if (status > 0) call print_error(status)

             !  Close the file and free the unit number.
             call ftclos(unit, status)
             if (status > 0) then
                call print_error(status)
             endif
             call ftfiou(unit, status)
             if (status > 0) then
                call print_error(status)
             endif
          endif !line has a map
       enddo !lines
       atom => NULL()
    enddo !atom


	if (limage_at_lam0) call write_image_map()

    return
  end subroutine write_atomic_maps
  
  subroutine write_image_map()
    integer :: status,unit,blocksize,bitpix,naxis
    integer, dimension(4) :: naxes
    integer :: group,fpixel,nelements
    logical :: simple, extend
    character(len=6) :: comment="VACUUM"
    real :: pixel_scale_x, pixel_scale_y
    character(len=10) :: wave_write


    blocksize=1
    simple=.true.
    extend=.false.
    group=1
    fpixel=1
    bitpix=-64

    naxis=4
    naxes(1)=npix_x
    naxes(2)=npix_y
    naxes(3)=RT_n_incl
    naxes(4)=RT_n_az
    nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)


    status=0
    call ftgiou (unit,status)
    write(wave_write, '(1I10)') nint(wavelength_ref)
    call ftinit(unit,"image_"//trim(adjustl(wave_write))//".fits.gz",blocksize,status)


    call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
    if (status > 0) call print_error(status)
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

    call ftpkys(unit,'BUNIT',"${\rm W \, m^{-2} \, Hz^{-1} \, pixel^{-1}}$",'${\rm F_{\nu}}$',status)
    call ftpkyd(unit,'lambda0',wavelength_ref,-7,'nm',status)

    !Write the array to the FITS file.
    call ftpprd(unit,group,fpixel,nelements,image_map,status)
    if (status > 0) call print_error(status)


    !  Close the file and free the unit number.
    call ftclos(unit, status)
    if (status > 0) then
        call print_error(status)
    endif
    call ftfiou(unit, status)
    if (status > 0) then
        call print_error(status)
    endif

    return
  end subroutine write_image_map

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


end module spectrum_type

! 	subroutine write_contribution_functions_ray()
! 	! -------------------------------------------------- !
! 	! Write contribution function to disk. + tau_one
! 	! --------------------------------------------------- !
! 		integer :: status,unit,blocksize,bitpix,naxis, naxis2
! 		integer, dimension(8) :: naxes, naxes2
! 		integer :: group,fpixel,nelements,la, icell
! 		logical :: simple, extend
! 		character(len=6) :: comment=""
! 		real :: pixel_scale_x, pixel_scale_y
!
! 		write(*,*)" -> writing contribution function for single ray ..."
!
! 		!  Get an unused Logical Unit Number to use to open the FITS file.
! 		status=0
! 		call ftgiou (unit,status)
!
!    !  Create the new empty FITS file.
! 		blocksize=1
! 		call ftinit(unit,"cntrb_ray.fits.gz",blocksize,status)
!
! 		simple=.true.
! 		extend=.true.
! 		group=1
! 		fpixel=1
!
! 		bitpix=-64
!
! 		if (lVoronoi) then
! 			naxis = 2
! 			naxes(1) = Nlambda
! 			naxes(2) = n_cells
! 			nelements = naxes(1) * naxes(2)
! 		else
! 			if (l3D) then
! 				naxis = 4
! 				naxes(1) = Nlambda
! 				naxes(2) = n_rad
! 				naxes(3) = 2*nz
! 				naxes(4) = n_az
! 				nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
! 			else
! 				naxis = 3
! 				naxes(1) = Nlambda
! 				naxes(2) = n_rad
! 				naxes(3) = nz
! 				nelements = naxes(1) * naxes(2) * naxes(3)
! 			end if
! 		end if
!
!   !  Write the required header keywords.
! 		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
! 		call ftpkys(unit,'UNIT',"W.m-2.Hz-1.sr-1.m-1",'cntrb I',status)
!   !  Write line CF to fits
! 		write(*,*) "Max,min abs(cntrb)=",maxval(abs(cntrb_ray)), minval(abs(cntrb_ray),mask=abs(cntrb_ray)>0)
! 		write(*,*) size(cntrb_ray), sizeof(cntrb_ray), nelements, shape(cntrb_ray)
! 		call ftpprd(unit,group,fpixel,nelements,cntrb_ray,status)
!
!
!   !  Close the file and free the unit number.
! 		call ftclos(unit, status)
! 		call ftfiou(unit, status)
!
!   !  Check for any error, and if so print out error messages
! 		if (status > 0) then
! 			call print_error(status)
! 		endif
!
! 		blocksize=1
! 		call ftinit(unit,"tauone_ray.fits.gz",blocksize,status)
!
! 		simple=.true.
! 		extend=.true.
! 		group=1
! 		fpixel=1
!
! 		bitpix=-64
!
! 		if (lVoronoi) then
! 			naxis = 2
! 			naxes(1) = Nlambda
! 			naxes(2) = n_cells
! 			nelements = naxes(1) * naxes(2)
! 		else
! 			if (l3D) then
! 				naxis = 4
! 				naxes(1) = Nlambda
! 				naxes(2) = n_rad
! 				naxes(3) = 2*nz
! 				naxes(4) = n_az
! 				nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
! 			else
! 				naxis = 3
! 				naxes(1) = Nlambda
! 				naxes(2) = n_rad
! 				naxes(3) = nz
! 				nelements = naxes(1) * naxes(2) * naxes(3)
! 			end if
! 		end if
!
!   !  Write the required header keywords.
! 		call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
! 		call ftpkys(unit,'UNIT'," ",'<tau>',status)
!   !  Write line CF to fits
! 		call ftpprd(unit,group,fpixel,nelements,tau_one_ray,status)
!
!
!   !  Close the file and free the unit number.
! 		call ftclos(unit, status)
! 		call ftfiou(unit, status)
!
!   !  Check for any error, and if so print out error messages
! 		if (status > 0) then
! 			call print_error(status)
! 		endif
!
! 	return
! 	end subroutine write_contribution_functions_ray
