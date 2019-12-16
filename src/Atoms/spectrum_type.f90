MODULE spectrum_type

  use atom_type, only : AtomicLine, AtomicContinuum, AtomType
  use atmos_type, only : GridType, atmos, helium, hydrogen
  use getlambda, only  : make_wavelength_grid, adjust_wavelength_grid, Read_wavelengths_table!,&
  ! 						Nred_array, Nmid_array, Nblue_array  
  use math, only : linear_1D
  !MCFOST's original modules
  use fits_utils, only : print_error
  use parametres
  use input!, only : nb_proc
  use mcfost_env, only : dp
  !use hdf5

  IMPLICIT NONE

   real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
   real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352
   character(len=*), parameter :: WAVES_FILE="atom_transfer_waves_grid.fits.gz"
   ! not that FLUX_FILE is 1D only if one pixel, otherwise it is a map.
   ! F(x,y,0,lambda) in several directions.
   character(len=*), parameter :: FLUX_FILE="flux.fits.gz", CF_FILE="cntrb.fits.gz"
   character(len=*), parameter :: FLUX_FILE_H5="flux.h5", CF_FILE_H5="cntrb.h5"
   

  !! Store S, Sc, Slte, Sclte, and jc and Kc, to save memory ???
  TYPE AtomicOpacity
   !active opacities
   real(kind=dp), allocatable, dimension(:,:)   :: chi, eta!, chic_nlte, etac_nlte
   ! NLTE magneto-optical elements and dichroism are stored in the background _p arrays.
   ! Mainly because we do not use them in SEE, even if etaQUV can be added to the total
   ! emissivity. But in case, etaQUV has to be atom dependent, so now we can store the LTE
   ! and NLTE is the same variable
   !passive opacities
   real(kind=dp), allocatable, dimension(:,:)   :: eta_p, chi_p
   real(kind=dp), allocatable, dimension(:,:)   :: eta_c, chi_c!, sca_c
   real(kind=dp), allocatable, dimension(:,:,:)   :: rho_p, chiQUV_p, etaQUV_p
   real(kind=dp), allocatable, dimension(:,:)   :: jc, jc_nlte, Kc_nlte
   !real(kind=dp), allocatable, dimension(:,:,:) :: Kc
   real(kind=dp), allocatable, dimension(:,:) :: Kc, sca_c
   !!logical, dimension(:), allocatable :: initialized
   									     !set to .true. for each cell, when iray=1.
   									     !.false. otherwise.
   									     !%initialized(id) = (iray == 1)
  END TYPE AtomicOpacity

  TYPE Spectrum
   !n_proc should be the last index
   type  (GridType), pointer :: atmos
   logical :: vacuum_to_air=.false., write_wavelength_grid=.false.
   integer :: Nwaves, Nact, Npass, Ntrans, NPROC=1
   real(kind=dp) :: wavelength_ref=0.d0 !nm optionnal
   real(kind=dp), allocatable, dimension(:) :: lambda
   !nlambda, nrays, nproc
   real(kind=dp), allocatable, dimension(:,:,:) :: I, StokesQ, StokesU, StokesV, Ic
   real(kind=dp), allocatable, dimension(:,:) :: Istar
   !nlambda, nproc
   real(kind=dp), allocatable, dimension(:,:) :: J, Jc, J20
   !Nlambda, xpix, ypix, Nincl, Naz
   real(kind=dp), allocatable, dimension(:,:,:,:,:) :: Flux, Fluxc
   real(kind=dp), allocatable, dimension(:,:,:,:,:,:) :: F_QUV
   !!real(kind=dp), allocatable, dimension(:,:) :: S_QUV
   !Contribution function
   !Nlambda,N_INCL, N_AZ, NCELLS
   real(kind=dp), allocatable, dimension(:,:,:,:) :: Ksi 
   ! Flux is a map of Nlambda, xpix, ypix, nincl, nazimuth
   real(kind=dp), allocatable, dimension(:,:,:) :: Psi, etau, tau !for cell icell in direction iray, thread id
   !size of Psi could change during the devlopment
   type (AtomicOpacity) :: AtomOpac
   character:: Jfile, J20file
  END TYPE Spectrum

  type (Spectrum) :: NLTEspec

  CONTAINS
  
 
 FUNCTION getlambda_limit(atom)
 !reddest wavelength of the sortest b-b transition in the ground state
  real(kind=dp) :: getlambda_limit
  type(AtomType) :: atom
  integer :: kr
  
  do kr=1, atom%Nline

   if (atom%lines(kr)%i==1) then !transition to the ground state
    getlambda_limit = min(1d6, NLTEspec%lambda(atom%lines(kr)%Nred))
   endif
  
  enddo
  
  !avoid problem
  getlambda_limit = max(NLTEspec%lambda(1), getlambda_limit)
  
 RETURN
 END FUNCTION getlambda_limit

  SUBROUTINE initSpectrum(lam0, vacuum_to_air)
   ! ------------------------------------------- !
    ! Allocate and store the wavelength grid for
    ! NLTE atomic line transfer.
   ! ------------------------------------------- !

   !integer, intent(in) :: NPROC
   integer :: kr, nat
   real(kind=dp), optional :: lam0
   logical, optional :: vacuum_to_air
   logical :: alloc_nlte_vars, alloc_image_vars
    
   NLTEspec%atmos => atmos
   alloc_nlte_vars = NLTEspec%atmos%Nactiveatoms > 0

   
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
   
   !only for H and He atm, not changed for image even if we remove lines
   hydrogen%scatt_limit = getlambda_limit(hydrogen)
   if (associated(helium)) helium%scatt_limit = getlambda_limit(helium)
   
   CALL writeWavelength()
   CALL allocSpectrum(alloc_nlte_vars)!.true. => alloc atom%eta, %chi, %Xcoupling if NLTE
  							 !.false. => if NLTE, skip this allocation
  							 !we probably do an image and we do not need these values.
  RETURN
  END SUBROUTINE initSpectrum
  
  !building
  SUBROUTINE initSpectrumImage()
  ! -------------------------------------------------------------------- !
   ! Allocate a special wavelength grid for emission_line map.
   ! This allows to solve for the RT equation only for a specific line.
   ! This fasten LTE calculation in case of 1 line.
  ! -------------------------------------------------------------------- !
   real(kind=dp), dimension(NLTEspec%Nwaves) :: old_grid
   integer, dimension(:), allocatable :: Nlam_R
   integer :: nat, kr, Ntrans_new, kp, kc
   !for Jnu
   integer :: icell, alloc_status, Nwaves_old
   real(kind=dp), dimension(:,:), allocatable :: Jnu, Jnuc

   
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
   
   !Reallocate J and Jc keeping their value if electron_scattering
   !NB: the spatial grid is the same
   !if (NLTEspec%atmos%electron_scattering) then
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
     
     !do sum over transitions
     do kr=1,NLTEspec%atmos%atoms(nat)%ptr_Atom%Ntr
      kc = NLTEspec%atmos%atoms(nat)%ptr_Atom%at(kr)%ik
      SELECT CASE (NLTEspec%atmos%atoms(nat)%ptr_Atom%at(kr)%trtype)
      
      CASE ("ATOMIC_LINE")
        write(*,*) "    b-b #", kc, NLTEspec%atmos%Atoms(nat)%ptr_atom%lines(kc)%lambda0,"nm"

      CASE ("ATOMIC_CONTINUUM")
        write(*,*) "    b-f #", kc, NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kc)%lambdamin, &
        "nm", NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kc)%lambda0, "nm"       
      CASE DEFAULT
       CALL error ("transition type unkown",NLTEspec%atmos%atoms(nat)%ptr_Atom%at(kr)%trtype)
      END SELECT
      
     end do
     Ntrans_new = Ntrans_new + NLTEspec%atmos%atoms(nat)%ptr_atom%Ntr
   end do
     
!      do kr=1,NLTEspec%atmos%Atoms(nat)%ptr_atom%Nline
!       if (NLTEspec%atmos%Atoms(nat)%ptr_atom%at(kr)%lcontrib_to_opac) then
!        kp = kp + 1
!        write(*,*) "    b-b #", kp, NLTEspec%atmos%Atoms(nat)%ptr_atom%lines(kr)%lambda0, &
!         "nm"
!        Ntrans_new = Ntrans_new + 1    
!       endif 
!      end do
!     write(*,*) "   --> ", kp, "b-b transitions"
!      kp = 0
!      do kr=1,NLTEspec%atmos%Atoms(nat)%ptr_atom%Ncont
!       if (NLTEspec%atmos%Atoms(nat)%ptr_atom%at(kr+NLTEspec%atmos%Atoms(nat)%ptr_atom%Nline)%lcontrib_to_opac) then
!        kp = kp + 1
!        write(*,*) "    b-f #", kp, NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%lambdamin, &
!         "nm", NLTEspec%atmos%Atoms(nat)%ptr_atom%continua(kr)%lambda0, "nm"  
!        Ntrans_new = Ntrans_new + 1  
!       endif
!      end do
!     write(*,*) "   --> ", kp, "b-f transitions"
!    end do
!    write(*,*) "  ** Total number of transitions for image:", Ntrans_new
   NLTEspec%Ntrans = Ntrans_new

   write(*,*) "  ** Total number of transitions for image:", NLTEspec%Ntrans
   !write(*,*) NLTEspec%atmos%Atoms(1)%ptr_atom%Nline, NLTEspec%atmos%Atoms(1)%ptr_atom%Ntr_line

   !reallocate wavelength arrays, except polarisation ? which are in adjustStokesMode
   CALL allocSpectrum(.false.) !do not realloc NLTE atom%chi;%eta;%Xcoupling
   							   !even if NLTEpops are present
   							   !NLTE eval_opartor condition never reached for images.
   							   
   !resample Jnu if any
!-> can be para if needed, but only once calc is needed
!   if (NLTEspec%Nact > 0 .and. NLTEspec%atmos%electron_scattering) then

   if (NLTEspec%atmos%electron_scattering) then
    write(*,*)  " -> Resample mean radiation field for isotropic scattering..." 

    Nwaves_old = size(old_grid)
    allocate(Jnu(Nwaves_old, NLTEspec%atmos%Nspace))
    Jnu(:,:) = NLTEspec%J(:,:)
    allocate(Jnuc(Nwaves_old, NLTEspec%atmos%Nspace))
    Jnuc(:,:) = NLTEspec%Jc(:,:)
   
    deallocate(NLTEspec%J, NLTEspec%Jc)
    !add  J20, J20c here
   
    !new freq grid
    allocate(NLTEspec%J(NLTEspec%Nwaves, atmos%Nspace),stat=alloc_status)
    if (alloc_status>0) call error("Allocation error, resample J")
    allocate(NLTEspec%Jc(NLTEspec%Nwaves, atmos%Nspace),stat=alloc_status)
    if (alloc_status>0) call error("Allocation error, resample Jc")

    NLTEspec%J(:,:) = 0d0
    NLTEspec%Jc(:,:) = 0d0
    !use better interpolation ??
    do icell=1, atmos%Nspace
     if (NLTEspec%atmos%icompute_atomRT(icell)>0) then
       NLTEspec%J(:,icell) = Linear_1D(Nwaves_old, old_grid,Jnu,NLTEspec%Nwaves,NLTEspec%lambda)
       NLTEspec%Jc(:,icell) = Linear_1D(Nwaves_old, old_grid,Jnuc,NLTEspec%Nwaves,NLTEspec%lambda)
     endif
    enddo  
   endif !resample Jnu
    write(*,*)  " ...Done!" 

   
  RETURN
  END SUBROUTINE initSpectrumImage
  
  !building
  SUBROUTINE initSpectrum_jnu(Nlambda,l0,l1)
  ! -------------------------------------------------------------------- !
   ! Allocate a special wavelength grid for Jnu calculation
  ! -------------------------------------------------------------------- !
   real(kind=dp), dimension(NLTEspec%Nwaves) :: old_grid
   integer, dimension(1) :: Nlam_R
   integer, intent(in) :: Nlambda
   real, intent(in) :: l0, l1
   integer :: la
   
   old_grid = NLTEspec%lambda
   write(*,*) " -> Defining a wavelength grid for Jnu.."
   CALL freeSpectrum() !first free waves arrays
   NLTEspec%atmos => atmos

   allocate(NLTEspec%lambda(Nlambda))
   NLTEspec%lambda(1) = l0
   do la=2, Nlambda
   	NLTEspec%lambda(la) = NLTEspec%lambda(la-1) + (l1-l0)/real(Nlambda-1)
   enddo
   
   Nlam_R(1) = Nlambda!Only one region with Nlambda points
   NLTEspec%Nwaves = size(NLTEspec%lambda)

   CALL adjust_wavelength_grid(old_grid, NLTEspec%lambda, Nlam_R, NLTEspec%atmos%Atoms)
 
   CALL allocSpectrum(.false.)
   							
   write(*,*) " -> done.."

  RETURN
  END SUBROUTINE initSpectrum_jnu
  
  SUBROUTINE reallocate_rays_arrays(newNray)
   integer, intent(in) :: newNray
   integer :: kr, nat, alloc_status
  
   NLTEspec%atmos%Nrays = newNray
   
   if ((NLTEspec%atmos%Nrays) /= newNray .or. (newNray /=1 )) &
   	 CALL Warning("  Beware, check the number of rays for the ray-traced map!")
   	 
   !reallocate phi_ray for NLTEloop only
   
!--> do not reallocate because it is used only during NLTEloop.
!    if we are here it means we are doing an image

!    do nat=1,atmos%Natom
!     do kr=1,atmos%Atoms(nat)%ptr_atom%Nline
!      deallocate(atmos%Atoms(nat)%ptr_atom%lines(kr)%phi_ray)
!      !only one ray and one proc need for images or LTE spectrum
!      allocate(atmos%Atoms(nat)%ptr_atom%lines(kr)%phi_ray(&
!      atmos%Atoms(nat)%ptr_atom%lines(kr)%Nlambda, newNray, NLTEspec%NPROC),stat=alloc_status)
!      if (alloc_status > 0) CALL ERROR("Allocation error line%phi_ray(:,newNray,id)")
!     enddo
!    enddo
   
   deallocate(NLTEspec%I, NLTEspec%Ic)
   !except polarization which are (de)allocated in adjustStokesMode
   if (NLTEspec%Nact > 0 .and.allocated(NLTEspec%PSI)) then 
    deallocate(NLTEspec%Psi)
    if (allocated(NLTEspec%etau)) deallocate(NLTEspec%etau)
    if (allocated(NLTEspec%tau)) deallocate(NLTEspec%tau)
   	allocate(NLTEspec%Psi(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
   	allocate(NLTEspec%tau(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))  
   end if 
    
   !Could be also LTE opac if line are kept in memory ?
   allocate(NLTEspec%I(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
   allocate(NLTEspec%Ic(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
   NLTEspec%I = 0d0; NLTEspec%Ic = 0d0
  
  RETURN
  END SUBROUTINE reallocate_rays_arrays

  SUBROUTINE allocSpectrum(alloc_atom_nlte)
   !Polarized quantities allocated in adjustStokesMode
   integer :: nat, k, Nlambda_max, alloc_status, istar
   type (AtomicContinuum) :: cont
   type (AtomicLine) 	  :: line
   logical, intent(in)    :: alloc_atom_nlte
   
   if (allocated(NLTEspec%I)) then
    write(*,*) "Error I already allocated"
    stop
   end if
   

   allocate(NLTEspec%Istar(NLTEspec%NWAVES,n_etoiles)); NLTEspec%Istar(:,:) = 0d0
   
   allocate(NLTEspec%I(NLTEspec%NWAVES, NLTEspec%atmos%NRAYS,NLTEspec%NPROC))
   allocate(NLTEspec%Ic(NLTEspec%NWAVES, NLTEspec%atmos%NRAYS,NLTEspec%NPROC))
   NLTEspec%I = 0d0
   NLTEspec%Ic = 0d0

   
   !Now opacities
   if (lstore_opac) then !keep continuum LTE opacities in memory
     !sca_c = Kc(:,:,2), chi_c = Kc(:,:,1), eta_c = jc
     allocate(NLTEspec%AtomOpac%Kc(NLTEspec%Nwaves,NLTEspec%atmos%Nspace), &
       NLTEspec%AtomOpac%jc(NLTEspec%Nwaves,NLTEspec%atmos%Nspace),NLTEspec%AtomOpac%sca_c(NLTEspec%Nwaves,NLTEspec%atmos%Nspace))
     NLTEspec%AtomOpac%Kc = 0d0
     NLTEspec%AtomOpac%jc = 0d0
     NLTEspec%AtomOpac%sca_c = 0.0_dp
   else
    allocate(NLTEspec%AtomOpac%eta_c(NLTEspec%Nwaves ,NLTEspec%NPROC))
    allocate(NLTEspec%AtomOpac%chi_c(NLTEspec%Nwaves ,NLTEspec%NPROC))
    !allocate(NLTEspec%AtomOpac%sca_c(NLTEspec%Nwaves,NLTEspec%NPROC))
    NLTEspec%AtomOpac%chi_c = 0.
    NLTEspec%AtomOpac%eta_c = 0.
    !NLTEspec%AtomOpac%sca_c = 0.
   end if

   !allocate(NLTEspec%AtomOpac%chic_nlte(NLTEspec%Nwaves, NLTEspec%NPROC),&
   !   NLTEspec%AtomOpac%etac_nlte(NLTEspec%Nwaves,NLTEspec%NPROC))
   !NLTEspec%AtomOpac%chic_nlte = 0.; NLTEspec%AtomOpac%etac_nlte = 0.

   allocate(NLTEspec%AtomOpac%eta_p(NLTEspec%Nwaves ,NLTEspec%NPROC))
   allocate(NLTEspec%AtomOpac%chi_p(NLTEspec%Nwaves ,NLTEspec%NPROC))

   NLTEspec%AtomOpac%chi_p = 0.
   NLTEspec%AtomOpac%eta_p = 0.
   
      
   !do not try to realloc after non-LTE for image.
   !Fursther, with labs=.false. for images, we do not enter in eval_operator condition
   !in NLTEOpacity()
   if (alloc_atom_nlte) then !NLTE loop activated
     allocate(NLTEspec%AtomOpac%chi(NLTEspec%Nwaves ,NLTEspec%NPROC))
     allocate(NLTEspec%AtomOpac%eta(NLTEspec%Nwaves ,NLTEspec%NPROC))
     NLTEspec%AtomOpac%chi = 0.
     NLTEspec%AtomOpac%eta = 0.
     allocate(NLTEspec%AtomOpac%Kc_nlte(NLTEspec%Nwaves,atmos%Nspace),stat=alloc_status)
     if (alloc_status >0) call error("Allocation error Kc_nlte")
     allocate(NLTEspec%AtomOpac%jc_nlte(NLTEspec%Nwaves,atmos%Nspace),stat=alloc_status)
     if (alloc_status >0) call error("Allocation error jc_nlte")
     NLTEspec%AtomOpac%Kc_nlte = 0.
     NLTEspec%AtomOpac%jc_nlte = 0.
          
    !if (NLTEspec%atmos%electron_scattering) CALL alloc_Jnu()
   
    !do not allocate if lxcoupling but no NLTE effects
    if (lxcoupling) NLTEspec%atmos%include_xcoupling = .true.
     
    allocate(NLTEspec%Psi(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC))
    
    if (.not.lxcoupling) then !not for MALI
       allocate(NLTEspec%tau(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC)) 
       allocate(NLTEspec%etau(NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%NPROC)) 
    endif  

    do nat=1,NLTEspec%atmos%Nactiveatoms
     !these are local to a cell, used to fill the Gamma matrix at a specific point: one cell per atom
     allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%eta(NLTEspec%Nwaves,NLTEspec%atmos%Nrays,NLTEspec%NPROC))
     allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%etac(NLTEspec%Nwaves,NLTEspec%NPROC),stat=alloc_status)
     if (alloc_status >0) call error("Allocation error atom%etac")!atmos%Nspace

     !Now the waves and angle integraed X coupling terms for each atom and for all transitions
     !(Sum over all active transitions for each transition) is kept.
     if (NLTEspec%atmos%include_xcoupling) then 
!        allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Xc&
!        (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Nlevel,NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Nlevel,&
!        NLTEspec%Nproc))
       
       !for testing
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Uji_down& ! U j->down
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Nlevel,NLTEspec%Nwaves, NLTEspec%atmos%Nrays,NLTEspec%Nproc))
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%chi_down& !chi j->down
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Nlevel,NLTEspec%Nwaves, NLTEspec%atmos%Nrays,NLTEspec%Nproc))
       allocate(NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%chi_up& ! i->up
       (NLTEspec%atmos%ActiveAtoms(nat)%ptr_atom%Nlevel,NLTEspec%Nwaves, NLTEspec%atmos%Nrays, NLTEspec%Nproc))
      !loop over atom%at if %U and %chi are stored for each transitions
     end if

    end do
  endif

  RETURN
  END SUBROUTINE allocSpectrum
  
  SUBROUTINE alloc_flux_image
   !allocated at the end only, to save memory
   !alloc also contribution function at the moment
   integer :: alloc_status
   real :: mem_alloc
   
   mem_alloc = real(npix_x*npix_y*rt_n_incl*rt_n_az*NLTEspec%Nwaves)/1024./1024.
   if (2*mem_alloc > 1e3) then
    write(*,*) " allocating ", 2 * mem_alloc/1024., " GB for flux arrays.."
   else
    write(*,*) " allocating ", 2 * mem_alloc, " MB for flux arrays.."  
   endif

    allocate(NLTEspec%Flux(NLTEspec%Nwaves,NPIX_X, NPIX_Y,RT_N_INCL,RT_N_AZ), stat=alloc_status)
    if (alloc_Status > 0) CALL ERROR ("Cannot allocate Flux")
    allocate(NLTEspec%Fluxc(NLTEspec%Nwaves,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ), stat=alloc_status)
    if (alloc_Status > 0) CALL ERROR ("Cannot allocate Flux continuum")

    NLTEspec%Flux = 0.0
    NLTEspec%Fluxc = 0.0
    
   !Contribution function
   
   if (lcontrib_function) then
    write(*,*) " CF not working ATM ... "
    lcontrib_function = .false.
    RETURN
   endif
   
   if (lcontrib_function .and. ltab_wavelength_image) then !the tab is here to prevent allocating a large array
		!hence ksi should be computed for small waves intervals.
	 !because otherwise, it is allocated with the alloc spectrum before initSpectrumImage()
!       write(*,*) "Trying to allocate",real(NLTEspec%atmos%Nspace)/1024. * real(NLTEspec%Nwaves)/1024. * &
!       	real(RT_N_INCL*RT_N_AZ)*real(npix_x*npix_y)/1024., " GB for ksi"  
!       allocate(NLTEspec%Ksi(NLTEspec%atmos%Nspace,NLTEspec%Nwaves,npix_x, npix_y, RT_N_INCL,RT_N_AZ),stat=alloc_status)
      write(*,*) "Trying to allocate",real(NLTEspec%atmos%Nspace)/1024. * real(NLTEspec%Nwaves)/1024. * &
      	real(RT_N_INCL*RT_N_AZ)/1024., " GB for ksi"  
      allocate(NLTEspec%Ksi(NLTEspec%atmos%Nspace,NLTEspec%Nwaves,RT_N_INCL,RT_N_AZ),stat=alloc_status)
      if (alloc_status > 0) then
       call ERROR('Cannot allocate ksi')
       lcontrib_function = .false.
      end if

      NLTEspec%Ksi(:,:,:,:) = 0d0
   else if (lcontrib_function .and. .not.ltab_wavelength_image) then
    CALL Warning(" Contribution function not taken into account. Use a wavelength table.")
    lcontrib_function = .false.
   end if

   
  RETURN
  END SUBROUTINE alloc_flux_image

  SUBROUTINE freeSpectrum() 
  
   deallocate(NLTEspec%lambda)
   deallocate(NLTEspec%Ic, NLTEspec%I, NLTEspec%Istar)
    if (allocated(NLTEspec%Flux)) deallocate(NLTEspec%Flux, NLTEspec%Fluxc)

   if (NLTEspec%atmos%Magnetized) then 
    !check allocation due to field_free sol
    if (allocated(NLTEspec%StokesQ)) & !same for all dangerous
    	deallocate(NLTEspec%StokesQ, NLTEspec%StokesU, NLTEspec%StokesV, NLTEspec%F_QUV)
    if (allocated(NLTEspec%AtomOpac%rho_p)) deallocate(NLTEspec%AtomOpac%rho_p)
    if (allocated(NLTEspec%AtomOpac%etaQUV_p)) deallocate(NLTEspec%AtomOpac%etaQUV_p)
    if (allocated(NLTEspec%AtomOpac%chiQUV_p)) deallocate(NLTEspec%AtomOpac%chiQUV_p)
   end if
   
   !active
   if (allocated(NLTEspec%AtomOpac%chi)) deallocate(NLTEspec%AtomOpac%chi,NLTEspec%AtomOpac%eta)
   !deallocate(NLTEspec%AtomOpac%etac_nlte, NLTEspec%AtomOpac%chic_nlte)
   if (allocated(NLTEspec%Psi)) then
    deallocate(NLTEspec%Psi)
    if (allocated(NLTEspec%tau)) deallocate(NLTEspec%tau)!, NLTEspec%AtomOpac%initialized)
   end if

   !passive
   deallocate(NLTEspec%AtomOpac%eta_p) !contains only passive lines if store_opac
   deallocate(NLTEspec%AtomOpac%chi_p)
   !deallocate(NLTEspec%AtomOpac%rho_p)
   if (lstore_opac) then !keep continuum LTE opacities in memory
     deallocate(NLTEspec%AtomOpac%Kc,  NLTEspec%AtomOpac%jc, NLTEspec%AtomOpac%sca_c)
   else !they are not allocated if we store background continua on ram
    deallocate(NLTEspec%AtomOpac%eta_c)
    deallocate(NLTEspec%AtomOpac%chi_c)
    !deallocate(NLTEspec%AtomOpac%sca_c)
   end if
   !elsewhere
   !if (NLTEspec%Nact > 0) deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)

   NULLIFY(NLTEspec%atmos)

   if (allocated(NLTEspec%Ksi)) deallocate(NLTEspec%ksi)


  RETURN
  END SUBROUTINE freeSpectrum

  SUBROUTINE initAtomOpac(id)!, eval_operator)
    ! set opacities to 0d0 for thread id
    integer, intent(in) :: id
    !logical, intent(in) :: eval_operator !: evaluate operator psi
    
    !need to be sure that id is > 0
!     if (id <= 0) then
!      write(*,*) "(initAtomOpac) thread id has to be >= 1!"
!      stop
!     end if

    !Should consider to allocate only if Nact now if I test
!     if (NLTEspec%Nact > 0) then !otherwise always 0
!      NLTEspec%AtomOpac%chi(:,id) = 0d0
!      NLTEspec%AtomOpac%eta(:,id) = 0d0
!      NLTEspec%AtomOpac%chic_nlte(:,id) = 0d0 !ray indep
!      NLTEspec%AtomOpac%etac_nlte(:,id) = 0d0
!     endif
    if (.not.lstore_opac) then !future remove
      NLTEspec%AtomOpac%chi_c(:,id) = 0d0
      NLTEspec%AtomOpac%eta_c(:,id) = 0d0
      !NLTEspec%AtomOpac%sca_c(:,id) = 0d0
    end if !else thay are not allocated
    NLTEspec%AtomOpac%chi_p(:,id) = 0d0
    NLTEspec%AtomOpac%eta_p(:,id) = 0d0

    
  RETURN
  END SUBROUTINE initAtomOpac
  
  SUBROUTINE initAtomOpac_nlte(id)!, eval_operator)
    ! set opacities to 0d0 for thread id
    integer, intent(in) :: id
    !logical, intent(in) :: eval_operator !: evaluate operator psi
    
    !need to be sure that id is > 0
!     if (id <= 0) then
!      write(*,*) "(initAtomOpac) thread id has to be >= 1!"
!      stop
!     end if

    NLTEspec%AtomOpac%chi(:,id) = 0d0
    NLTEspec%AtomOpac%eta(:,id) = 0d0
    
!     NLTEspec%AtomOpac%chic_nlte(:,id) = 0d0 
!     NLTEspec%AtomOpac%etac_nlte(:,id) = 0d0


    
  RETURN
  END SUBROUTINE initAtomOpac_nlte
  
  SUBROUTINE initAtomOpac_zeeman(id)!, eval_operator)
    ! set opacities to 0d0 for thread id
    integer, intent(in) :: id
    !logical, intent(in) :: eval_operator !: evaluate operator psi
    
    !need to be sure that id is > 0
!     if (id <= 0) then
!      write(*,*) "(initAtomOpac) thread id has to be >= 1!"
!      stop
!     end if


    
    !Currently LTE or NLTE Zeeman opac are not stored on memory. They change with 
    !direction. BUT the star is assumed to not emit polarised photons (from ZeemanEffect)
    !So we do not take into account this opac in Metal_lambda and futur NLTEOpacity_lambda
    !if (NLTEspec%atmos%magnetized .and. PRT_SOLUTION == "FULL_STOKES") then
    !check allocation because even if magnetized, due to the FIELD_FREE solution or WF
    !they might be not allocated !if (allocated(NLTEspec%AtomOpac%rho_p)) 
    NLTEspec%AtomOpac%rho_p(:,:,id) = 0d0
    NLTEspec%AtomOpac%etaQUV_p(:,:,id) = 0d0
    NLTEspec%AtomOpac%chiQUV_p(:,:,id) = 0d0
     !both NLTE and LTE actually.
     !If one want to add them in SEE, it has to be atom (atom%eta) dependent for emissivity.
     !adding the off diagonal elements in SEE results in solving for the whole
     !density matrix, WEEEEEEELLLLL beyond our purpose.
    !end if
    
  RETURN
  END SUBROUTINE initAtomOpac_zeeman
  
  SUBROUTINE alloc_Jnu()
   !only allocate Jc for the moment and use it even for lines
  
     !allocate(NLTEspec%J(NLTEspec%Nwaves,NLTEspec%atmos%Nspace))!%NPROC))
     allocate(NLTEspec%Jc(NLTEspec%Nwaves,NLTEspec%atmos%Nspace))
   !Just to try
   !CALL Warning("(allocSpectrum()) J20 allocated")
   !allocate(NLTEspec%J20(NLTEspec%Nwaves,NLTEspec%NPROC)); NLTEspec%J20 = 0d0
   !  NLTEspec%J = 0.0
     NLTEspec%Jc = 0.0

   
  RETURN
  END SUBROUTINE alloc_Jnu
  
  SUBROUTINE dealloc_Jnu()
    
    if (allocated(NLTEspec%J)) deallocate(NLTEspec%J)
    if (allocated(NLTEspec%Jc)) deallocate(NLTEspec%Jc)
    if (allocated(NLTEspec%J20)) deallocate(NLTEspec%J20)

  RETURN 
  END SUBROUTINE dealloc_Jnu
  
  SUBROUTINE init_psi_operator(id, iray)
    integer, intent(in) :: iray, id
    integer :: nact
    
   	NLTEspec%Psi(:,iray,id) = 0d0
   	if (.not.lxcoupling) then
   		NLTEspec%tau(:,iray,id) = 0d0
   		NLTEspec%etau(:,iray,id) = 0d0 !keep only exp(-dtau), dtau itself not needed,
   									   !avoid te recompute the expo each time
   	endif
   	
   	do nact=1,NLTEspec%atmos%NactiveAtoms
   		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%eta(:,iray,id) = 0d0
   		if (lxcoupling) then
    		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%Uji_down(:,:,iray,id) = 0d0 
    		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%chi_down(:,:,iray,id) = 0d0  	 
    		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%chi_up(:,:,iray,id) = 0d0	 
   		endif
   	enddo
  
  RETURN
  END SUBROUTINE init_psi_operator
  
  SUBROUTINE init_psi_operator_new(id)
    !for one cell
    integer, intent(in) :: id
    integer :: nact
    
   	NLTEspec%Psi(:,:,id) = 0d0
   	if (.not.lxcoupling) then
   	 NLTEspec%tau(:,:,id) = 0d0
   	 NLTEspec%etau(:,:,id) = 0d0
   	endif
   	
   	do nact=1,NLTEspec%atmos%NactiveAtoms
   		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%eta(:,:,id) = 0d0
   		if (lxcoupling) then
    		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%Uji_down(:,:,:,id) = 0d0 
    		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%chi_down(:,:,:,id) = 0d0  	 
    		NLTEspec%atmos%ActiveAtoms(nact)%ptr_atom%chi_up(:,:,:,id) = 0d0	 
   		endif
   	enddo
  
  RETURN
  END SUBROUTINE init_psi_operator_new

  
  SUBROUTINE alloc_weights()
  ! --------------------------------------------------- !
   ! 
  ! --------------------------------------------------- !
   use atmos_type, only : atmos
   integer :: kr, nact
   type(AtomicLine) :: line
   type(AtomicContinuum) :: cont
   
   do nact=1,atmos%Nactiveatoms
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
       line = atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)
       !allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%wphi(NLTEspec%Nproc))
       allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%w_lam(line%Nlambda))
       allocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%phi_ray(line%Nlambda,NLTEspec%atmos%Nrays, NLTEspec%NPROC))

    end do
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Ncont
       cont = atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)
       allocate(atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)%w_lam(cont%Nlambda))
       !allocate(atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)%wmu(NLTEspec%Nproc))
    end do
   end do
 
  RETURN
  END SUBROUTINE alloc_weights
  
  SUBROUTINE dealloc_weights()
  ! --------------------------------------------------- !
   ! 
  ! --------------------------------------------------- !
   use atmos_type, only : atmos
   integer :: kr, nact
   
   do nact=1,atmos%Nactiveatoms
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
       deallocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%w_lam)!wphi)
       deallocate(atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%phi_ray)
    end do
    do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Ncont
       deallocate(atmos%ActiveAtoms(nact)%ptr_atom%continua(kr)%w_lam)!%wmu)
    end do
   end do
 
  RETURN
  END SUBROUTINE dealloc_weights



 SUBROUTINE WRITE_FLUX()
 ! -------------------------------------------------- !
  ! Write the spectral Flux map on disk.
  ! FLUX map:
  ! NLTEspec%Flux total and NLTEspec%Flux continuum
 ! --------------------------------------------------- !
  use input
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(6) :: naxes
  integer :: group,fpixel,nelements, i, xcenter
  integer :: la, Nred, Nblue, kr, kc, m, Nmid
  logical :: simple, extend
  character(len=6) :: comment="VACUUM"
  real(kind=dp) :: lambda_vac(NLTEspec%Nwaves), Fnu
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
  if (status > 0) then
     CALL print_error(status)
  endif

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
  if (status > 0) then
     CALL print_error(status)
  endif

  ! create new hdu for continuum
  CALL ftcrhd(unit, status)
  if (status > 0) then
     CALL print_error(status)
  endif

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
  if (status > 0) then
     CALL print_error(status)
  endif
  ! write polarized flux if any. Atmosphere magnetic does not necessarily
  								!means we compute polarization
  if ((NLTEspec%atmos%magnetized) .and. (PRT_SOLUTION /= "NO_STOKES") &
               .and. (RT_line_method /= 1)) then
   write(*,*) " -> Writing polarization"
   CALL ftcrhd(unit, status)
   if (status > 0) then
     CALL print_error(status)
   endif
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
   if (status > 0) then
     CALL print_error(status)
   endif
  end if
  
  ! create new hdu for wavelength grid
  CALL ftcrhd(unit, status)
  
  if (status > 0) then
     CALL print_error(status)
  endif
  
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
   if (status > 0) then
     CALL print_error(status)
   endif

  !  Close the file and free the unit number.
  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     CALL print_error(status)
  endif

 RETURN
 END SUBROUTINE WRITE_FLUX
 
 SUBROUTINE WRITE_FLUX_ASCII()
 ! -------------------------------------------------- !
  ! written only for the first inclination / azimuth
 ! --------------------------------------------------- !
  use input
  integer :: status,unit
  integer :: la, j, i
  real(kind=dp), dimension(:), allocatable :: Fnu, Fnuc
  
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
   status=0
   unit = 10
   allocate(Fnu(NLTEspec%Nwaves), Fnuc(NLTEspec%Nwaves))
   Fnu = 0.
   Fnuc = 0.
   do i=1, npix_x
    do j=1, npix_y
     Fnu = Fnu + NLTEspec%Flux(:,i,j,1,1)
     Fnuc = Fnuc + NLTEspec%Fluxc(:,i,j,1,1)
    enddo
   enddo
   
   open(unit,file="flux.s", status='unknown', iostat=status)
    do la=1, NLTEspec%Nwaves
      write(unit, *) NLTEspec%lambda(la), Fnu(la), fnuc(la)
    enddo
   close(unit)
   
   deallocate(Fnu, Fnuc)

 RETURN
 END SUBROUTINE WRITE_FLUX_ASCII
  
  SUBROUTINE writeWavelength()
  ! --------------------------------------------------------------- !
   ! Write wavelength grid build with each transition of each atom
  ! --------------------------------------------------------------- !
  use input
   integer :: unit, EOF = 0, blocksize, naxes(1), naxis,group, bitpix, fpixel
   logical :: extend, simple
   character(len=6) :: comment="VACUUM"
   real(kind=dp) :: lambda_vac(NLTEspec%Nwaves)
   
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


  RETURN
  END FUNCTION vacuum2air

  FUNCTION air2vacuum(Nlambda, lambda_air) result(lambda_vac)
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

  RETURN
  END FUNCTION air2vacuum
  
 !building , move to write_opacity.f90
 SUBROUTINE WRITE_CNTRB_FUNC_PIX()
 ! -------------------------------------------------- !
  ! Write contribution function to disk.
 ! --------------------------------------------------- !
 use input
  integer :: status,unit,blocksize,bitpix,naxis, naxis2
  integer, dimension(8) :: naxes, naxes2
  integer :: group,fpixel,nelements, nelements2
  logical :: simple, extend
  character(len=6) :: comment=""
  real :: pixel_scale_x, pixel_scale_y
  
   write(*,*)" -> writing contribution function..."
  
   !  Get an unused Logical Unit Number to use to open the FITS file.
   status=0
   CALL ftgiou (unit,status)

   !  Create the new empty FITS file.
   blocksize=1
   CALL ftinit(unit,TRIM(CF_FILE),blocksize,status)

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
   else
    naxis = 5
    naxes(1) = n_rad
    naxes(2) = nz
    naxes(3) = NLTEspec%Nwaves
    naxes(4) = RT_n_incl
    naxes(5) = RT_n_az
    nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4) * naxes(5)
   end if
  end if

  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  CALL ftpkys(unit,'UNIT',"W.m-2.Hz-1.pixel-1",'Ksi',status)
  !  Write line CF to fits
  write(*,*) "Max,min abs(Ksi)=",maxval(dabs(NLTEspec%Ksi)), minval(dabs(NLTEspec%Ksi),mask=dabs(NLTEspec%Ksi)>0)
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%ksi,status)


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

