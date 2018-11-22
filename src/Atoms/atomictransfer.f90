!
! This module solves for the radiative transfer equation by ray-tracing for 
! multi-level atoms, using the MALI scheme.
!
! Outputs:
! - Flux (\lambda) [J.s^{-1}.m^{-2}.Hz^{-1}]
! - Irradiation map around a line [J.s^{-1}.m^{-2}.Hz^{-1}.pix^{-1}] !sr^{-1}
!   --> not exactly, it is multiplied by the pixel solid angle seen from the Earth.
! TBD - levels population n
! TBD - Cooling rates PHIij = njRij - niRji
! TBD - Contribution function around a line
! - Electron density
! TBD - Fixed radiative rates
! TBD - Full NLTE line transfer
! TBD - PFR on the "atoms'frame - observer's frame" approach
!
! It uses all the NLTE modules in NLTE/ and it is called in mcfost.f90 similar to
! mol_transfer.f90
!
! Note: SI units, velocity in m/s, density in kg/m3, radius in m
MODULE AtomicTransfer

 use metal, only : Background
 use spectrum_type
 use atmos_type
 use readatom
 use lte
 use collision
 use solvene
 use writeatom
 use readatmos, only : readatmos_1D !for testing
 
 !$ use omp_lib
 
 !MCFOST's original modules
 use input
 use parametres
 use grid
 use density
 use dust_prop
 use dust_transfer, only : compute_stars_map
 use dust_ray_tracing, only : init_directions_ray_tracing , & 
                              tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az, tab_RT_incl, & 
                              stars_map, kappa
 use stars
 use wavelengths
 
 IMPLICIT NONE


 CONTAINS

 
 SUBROUTINE INTEG_RAY_LINE(id,icell_in,x,y,z,u,v,w,iray,labs)
  ! This routine performs integration of the transfer equation along a ray
  ! crossing different cells.
  ! --> Atomic Lines case.
  !
  ! voir integ_ray_mol from mcfost
  ! All necessary quantities are initialised before this routine, and everything
  ! call, update, rewrite atoms% and spectrum%
  ! if atmos%nHtot or atmos%T is 0, the cell is empty
  !
  ! id = processor information, iray = index of the running ray
  integer, intent(in) :: id, icell_in, iray
  double precision, intent(in) :: u,v,w
  double precision, intent(in) :: x,y,z
  logical, intent(in) :: labs


  logical :: re_init ! to reset opac to 0 for next cell point
  double precision :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  double precision, dimension(NLTEspec%Nwaves) :: dtau, dtau2, Snu
  double precision, dimension(NLTEspec%Nwaves) :: tau, tau2
  double precision, dimension(NLTEspec%Nwaves) :: tau_c
  double precision, dimension(NLTEspec%Nwaves) :: dtau_c, Snu_c
  integer :: nbr_cell, icell, next_cell, previous_cell
  double precision :: facteur_tau
  logical :: lcellule_non_vide

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp !go from surface down to the star
  tau_c = 0.0_dp

  ! Reset, because when we compute the flux map
  ! with emission_line_map, we only use one ray, iray=1
  ! and the same space is used for emergent I.
  ! Therefore it is needed to reset it.
  NLTEspec%I(id,:,iray) = 0d0
  NLTEspec%Ic(id,:,iray) = 0d0

  !*** propagation dans la grille ***!
  !---------------------------------------------!

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
    ! Indice de la cellule
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     lcellule_non_vide=.true.
    else
     lcellule_non_vide=.false.
    endif
       
    ! Test sortie
    if (test_exit_grid(icell, x0, y0, z0)) then 
     RETURN
    end if
    
    nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
    previous_cell = 0 ! unused, just for Voronoi
    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, &
                     l, l_contrib, l_void_before)
                     
    if (.not.atmos%lcompute_atomRT(icell)) & 
         lcellule_non_vide = .false. !chi and chi_c = 0d0, cell is transparent                     

    if (lcellule_non_vide) then
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m ! fabs(s(k+1)-s(k)) along a ray ??
     ! a quoi il sert lui ?
        
!bug here ?
!         atmos_parcel%x = x0
!         atmos_parcel%y = y0
!         atmos_parcel%z = z0
!         atmos_parcel%u = u
!         atmos_parcel%v = v
!         atmos_parcel%w = w

     CALL initAS(id, re_init=.true.) !set opac to 0 for this cell and thread id
     !! Compute background opacities for PASSIVE bound-bound and bound-free transitions
     !! at all wavelength points including vector fields in the bound-bound transitions
     CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w) !x,y,z,u,v,w,x1,y1,z1 
                                !define the projection of the vector field (velocity, B...)
                                !at each spatial location.
     ! Epaisseur optique
     ! dtau = -chi*ds --> ds < 0 and ds in meter
     dtau(:) =  l_contrib * (NLTEspec%ActiveSet%chi_c(id,:)+NLTEspec%ActiveSet%chi(id,:)) !scattering + thermal
     dtau_c(:) = l_contrib * NLTEspec%ActiveSet%chi_c_bf(id,:)

     ! Source function
     ! No dust yet
     ! J and Jc are the mean radiation field for total and continuum intensities
     ! it multiplies the continuum scattering coefficient for isotropic (unpolarised)
     ! scattering.
     Snu = (NLTEspec%ActiveSet%eta_c(id,:) + NLTEspec%ActiveSet%eta(id,:) + & !NLTE lines 
                  NLTEspec%ActiveSet%sca_c(id,:) * NLTEspec%J(id,:)) / & !+(sca_c * (3mu2-1)*J20)
                 (NLTEspec%ActiveSet%chi_c(id,:) + NLTEspec%ActiveSet%chi(id,:)) ! LTE+NLTE
     ! continuum source function
     Snu_c = (NLTEspec%ActiveSet%eta_c_bf(id,:) + & 
            NLTEspec%ActiveSet%sca_c(id,:) * NLTEspec%Jc(id,:)) / NLTEspec%ActiveSet%chi_c_bf(id,:)

     NLTEspec%I(id,:,iray) = NLTEspec%I(id,:,iray) + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu
     NLTEspec%Ic(id,:,iray) = NLTEspec%Ic(id,:,iray) + exp(-tau_c) * (1.0_dp - exp(-dtau_c)) * Snu_c
!     NLTEspec%I(id,:,iray) = NLTEspec%I(id,:,iray)*exp(-dtau) + Snu * exp(-tau) * dtau
!     NLTEspec%Ic(id,:,iray) = NLTEspec%Ic(id,:,iray)*exp(-dtau_c) + Snu_c * exp(-tau_c) * dtau_c

     ! surface superieure ou inf
     facteur_tau = 1.0
     if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
     if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

     ! Mise a jour profondeur optique pour cellule suivante
     tau = tau + dtau * facteur_tau
     tau_c = tau_c + dtau_c

     ! Define PSI Operators here

!      !set opacities to 0.0 for next cell point, for the thread id.
!      CALL initAS(id, re_init=.true.)
    end if  ! lcellule_non_vide
  end do infinie
  !---------------------------------------------!
  
  RETURN
  END SUBROUTINE INTEG_RAY_LINE
 
  SUBROUTINE FLUX_PIXEL_LINE(&
         id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! see: mol_transfer.f90/intensite_pixel_mol()

   integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
   double precision, dimension(3), intent(in) :: pixelcorner,dx,dy
   double precision, intent(in) :: pixelsize,u,v,w
   integer, parameter :: maxSubPixels = 32
   double precision :: x0,y0,z0,u0,v0,w0
   double precision, dimension(NLTEspec%Nwaves) :: Iold, nuHz, I0, I0c
   double precision, dimension(3) :: sdx, sdy
   double precision:: npix2, diff
   double precision, parameter :: precision = 1.e-2
   integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter

   logical :: lintersect, labs

   labs = .false.
   ! reset local Fluxes
   I0c = 0d0
   I0 = 0d0

   ! Ray tracing : on se propage dans l'autre sens
   u0 = -u ; v0 = -v ; w0 = -w

   Iold = 0d0

   ! le nbre de subpixel en x est 2^(iter-1)
   subpixels = 1
   iter = 1

   infinie : do ! Boucle infinie tant que le pixel n'est pas converge
     npix2 =  real(subpixels)**2
     Iold = I0
     I0 = 0d0
     I0c = 0d0
     ! Vecteurs definissant les sous-pixels
     sdx(:) = dx(:) / real(subpixels,kind=dp)
     sdy(:) = dy(:) / real(subpixels,kind=dp)

     iray = 1 ! because the direction is fixed and we compute the flux emerging
               ! from a pixel, by computing the Intensity in this pixel

     ! L'obs est en dehors de la grille
     ri = 2*n_rad ; zj=1 ; phik=1

     ! Boucle sur les sous-pixels qui calcule l'intensite au centre
     ! de chaque sous pixel
     do i = 1,subpixels
        do j = 1,subpixels
           ! Centre du sous-pixel
           x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
           y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
           z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)
           ! On se met au bord de la grille : propagation a l'envers
           CALL move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
!            write(*,*) i, j, lintersect, labs, icell, x0, y0, z0, u0, v0, w0
!            stop
           
           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
             CALL INTEG_RAY_LINE(id, icell, x0,y0,z0,u0,v0,w0,iray,labs)
             I0 = I0 + NLTEspec%I(id,:,iray)
             I0c = I0c + NLTEspec%Ic(id,:,iray)
           endif
        end do !j
     end do !i

     I0 = I0 / npix2
     I0c = I0c / npix2

     if (iter < n_iter_min) then
        ! On itere par defaut
        subpixels = subpixels * 2
     else if (iter >= n_iter_max) then
        ! On arrete pour pas tourner dans le vide
        ! write(*,*) "Warning : converging pb in ray-tracing"
        ! write(*,*) " Pixel", ipix, jpix
        exit infinie
     else
        ! On fait le test sur a difference
        diff = maxval( abs(I0 - Iold) / (I0 + 1e-300_dp) )
        if (diff > precision ) then
           ! On est pas converge
           subpixels = subpixels * 2
        else
           exit infinie
        end if
     end if ! iter
     iter = iter + 1
   end do infinie
  !Prise en compte de la surface du pixel (en sr)
  !correction en D**2, en AU, 'cause pixelsize is in AU
  !distance = Rmax  ! (AU) je me mets au bord du modèle.
  nuHz = 1d0 !c_light / NLTEspec%lambda * 1d9 !s^-1 !to get W/m2 instead of W/m2/Hz
  ! Flux out of a pixel in W/m2/Hz
  I0 = nuHz * I0 * (pixelsize / (distance*pc_to_AU) )**2
  I0c = nuHz * I0c * (pixelsize / (distance*pc_to_AU) )**2
  
  ! adding to the total flux map.
  NLTEspec%Flux(:,ipix,jpix,ibin,iaz) = NLTEspec%Flux(:,ipix,jpix,ibin,iaz) + I0
  NLTEspec%Fluxc(:,ipix,jpix,ibin,iaz) = NLTEspec%Fluxc(:,ipix,jpix,ibin,iaz) + I0c

  RETURN
  END SUBROUTINE FLUX_PIXEL_LINE
  
 SUBROUTINE EMISSION_LINE_MAP(ibin,iaz)
  ! Line emission map in a given direction n(ibin,iaz).
  ! if only one pixel it gives the Flux.
  ! See: emission_line_map in mol_transfer.f90
  
  integer, intent(in) :: ibin, iaz !define the direction in which the map is computed
  double precision :: x0,y0,z0,l,u,v,w

  double precision, dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  double precision, dimension(3,nb_proc) :: pixelcorner
  double precision:: taille_pix, nuHz
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36
  integer, parameter :: n_ray_star = 1000
  double precision, dimension(n_rad_RT) :: tab_r
  double precision:: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda
  logical :: lresolved = .false.
  
npix_x = 101; npix_y = 101
  write(*,*) "incl (deg) = ", tab_RT_incl(ibin), "azimuth (deg) = ", tab_RT_az(iaz)

  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin)
  uvw = (/u,v,w/) !vector position

  ! Definition des vecteurs de base du plan image dans le repere universel
  ! Vecteur x image sans PA : il est dans le plan (x,y) et orthogonal a uvw
  x = (/cos(tab_RT_az(iaz) * deg_to_rad),sin(tab_RT_az(iaz) * deg_to_rad),0._dp/)

  ! Vecteur x image avec PA
  if (abs(ang_disque) > tiny_real) then
     ! Todo : on peut faire plus simple car axe rotation perpendiculaire a x
     x_plan_image = rotation_3d(uvw, ang_disque, x)
  else
     x_plan_image = x
  endif

  ! Vecteur y image avec PA : orthogonal a x_plan_image et uvw
  y_plan_image = -cross_product(x_plan_image, uvw)

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin ! in AU

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Coin en bas gauche de l'image
  Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)
  ! presently only squared pixels

     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
     lresolved = .true.
     
     !write(*,*) taille_pix, map_size, zoom, npix_x, npix_y, RT_n_incl, RT_n_az

     dx(:) = x_plan_image * taille_pix
     dy(:) = y_plan_image * taille_pix

     if (l_sym_ima) then
        npix_x_max = npix_x/2 + modulo(npix_x,2)
     else
        npix_x_max = npix_x
     endif

     !$omp parallel &
     !$omp default(none) &
     !$omp private(i,j,id) &
     !$omp shared(Icorner,pixelcorner,dx,dy,u,v,w,taille_pix,npix_x_max,npix_y) &
     !$omp shared(n_iter_min,n_iter_max,ibin,iaz)

     ! loop on pixels
     id = 1 ! pour code sequentiel
     n_iter_min = 1 ! 3
     n_iter_max = 1 ! 6
     
     !$omp do schedule(dynamic,1)
     do i = 1,npix_x_max
        !$ id = omp_get_thread_num() + 1
        do j = 1,npix_y
           !write(*,*) i,j
           ! Coin en bas gauche du pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)

           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, &
                      i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        enddo !j
     enddo !i
     !$omp end do
     !$omp end parallel
     

  ! --------------------------
  ! Ajout Flux etoile
  ! --------------------------
  write(*,*) " --> adding stellar flux map..."
  do i = 1, NLTEspec%Nwaves
   nuHz = c_light / NLTEspec%lambda(i) * 1d9 !if NLTEspec%Flux in W/m2 set nuHz = 1d0
                                             !else it means that in FLUX_PIXEL_LINE, nuHz
                                             !is 1d0 (to have flux in W/m2/Hz)
   CALL compute_stars_map(i, u, v, w, taille_pix, dx, dy, lresolved)
!    write(*,*) "Adding the star", NLTEspec%lambda(i), tab_lambda(i)*1000, & 
!       " maxFstar = ", MAXVAL(stars_map(:,:,1)), MINVAL(stars_map(:,:,1))
   NLTEspec%Flux(i,:,:,ibin,iaz) = NLTEspec%Flux(i,:,:,ibin,iaz) + stars_map(:,:,1) / nuHz
   NLTEspec%Fluxc(i,:,:,ibin,iaz) = NLTEspec%Fluxc(i,:,:,ibin,iaz) + stars_map(:,:,1) / nuHz
  end do

 RETURN
 END SUBROUTINE EMISSION_LINE_MAP
 
 ! NOTE: should inverse the order of frequencies and depth because in general
 !       n_cells >> n_lambda, in "real" cases.
 
 SUBROUTINE Atomic_transfer()
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
#include "sprng_f.h"

  integer :: atomunit = 1, nact, nat, la !atoms and wavelength
  integer :: icell !spatial variables
  integer :: ibin, iaz
  logical :: re_init, labs, lkeplerian
  
  !testing vars to deleted
  double precision :: Ttmp(n_cells), nHtot(n_cells), v_turb(n_cells), vchar, netmp(n_cells)
  !!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the following quantities are input parameters !!
  integer :: NiterMax = 20, Nrays = 1! Number of rays for angular integration and to compute Inu(mu)
  integer :: IterLimit
  logical :: SOLVE_FOR_NE = .false. !for calculation of electron density even if atmos%calc_ne
                                   ! is .false.
  character(len=7) :: NE0 = "NEMODEL"
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
npix_x = 101; npix_y = 101

  atmos%Nrays = Nrays
  if (atmos%Nrays == 0) then
   write(*,*) "Nrays should at least be 1!"
   stop
  end if
! --> move elsewhere (for instance before NLTE loop)
!   if (atmos%Nrays==1) then
!    write(*,*) "Solving for", atmos%Nrays,' direction.'
!   else
!    write(*,*) "Solving for", atmos%Nrays,' directions.'
!   end if
  
!! -------------------------------------------------------- !!
  
   nHtot = 2d21 * densite_gaz/MAXVAL(densite_gaz)
   !nHtot =  1d6 * densite_gaz * masse_mol_gaz / m3_to_cm3 / masseH
   Ttmp = Tdust * 100d0 !100d0, depends on the stellar flux
   netmp = 1d-2 * nHtot
!    nHtot = 1.7d15
!    Ttmp = 9.400000d03
!    netmp = 3.831726d15
  
  !!!! Lyman core
!   nHtot = 2.2770d15
!   Ttmp = 2.406000d4
!   netmp = 4.291177d16
  !!!! Lyman wing
!   nHtot = 8.8304d15
!   Ttmp = 1.452000d4
!   netmp = 6.137226d16
!  
  ! more or less the same role as init_molecular_disk
  CALL init_atomic_atmos(n_cells, Ttmp, netmp, nHtot)
 
  ! OR READ FROM MODEL (to move elsewhere) 
  !suppose the model is in utils/Atmos/
  CALL readatmos_1D("Atmos/FALC_mcfost.fits.gz")
  
  write(*,*) "maxTgas = ", MAXVAL(atmos%T), " minTgas = ", MINVAL(atmos%T)
  write(*,*) "maxnH = ", MAXVAL(atmos%nHtot), " minnH = ", MINVAL(atmos%nHtot)
  write(*,*) "maxNE = ", MAXVAL(atmos%ne), " minNE = ", MINVAL(atmos%ne)

!! -------------------------------------------------------- !!


  ! if the electron density is not provided by the model or simply want to
  ! recompute it
  ! if atmos%ne given, but one want to recompute ne
  ! else atmos%ne not given, its computation is mandatory
  if (.not.atmos%calc_ne) atmos%calc_ne = SOLVE_FOR_NE
  if (SOLVE_FOR_NE) write(*,*) "(Force) Solving for electron density"
  if (atmos%calc_ne) CALL SolveElectronDensity(atmos%ne,NE0)
  CALL writeElectron() !will be moved elsewhere
  ! do it in the reading process
  CALL writeHydrogenDensity()



  !Read atomic models and allocate space for n, nstar, vbroad, ntotal, Rij, Rji
  ! on the whole grid space.
  ! The following routines have to be invoked in the right order !
  CALL readAtomicModels(atomunit)
  NLTEspec%atmos => atmos
  CALL initSpectrum(nb_proc, 500d0, .false., .true.)
  re_init = .false. !first allocation it is done by setting re_init = .false.
  !when re_init is .false., table for opacities are allocated for all wavelengths and
  !all threads, and set to 0d0.
  ! when it is .true., opacities are set to zero for next points, for a specific thread.
  CALL allocSpectrum(npix_x, npix_y, RT_n_incl, RT_n_az)
  CALL initAS(0, re_init) !zero because when re_init=.false. it is independent of the
  ! threads. 0 ensures that an error occurs if the allocation is thread-dependent.
  !Compute LTE populations for all atoms, nstar. Open collision file for active atoms
  ! compute nHmin density from neutral hydrogen density (sum over neutral levels)
  Call setLTEcoefficients ()

  ! ----- ALLOCATE SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ------!
  !--> should move them to init_atomic_atmos ? or elsewhere
  !need to be deallocated at the end of molecule RT or its okey ?
  if (allocated(ds)) deallocate(ds)
  allocate(ds(atmos%Nrays,NLTEspec%NPROC))
  ds = 0d0 !meters
  CALL init_directions_ray_tracing()
  if (.not.allocated(stars_map)) allocate(stars_map(npix_x,npix_y,3))
  stars_map = 0
  n_lambda = NLTEspec%Nwaves
  if (allocated(tab_lambda)) deallocate(tab_lambda)
  allocate(tab_lambda(n_lambda))
  if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
  allocate(tab_delta_lambda(n_lambda))
  tab_lambda = NLTEspec%lambda * 1d-3 !nm to micron
  tab_delta_lambda(1) = 0d0
  do la=2,NLTEspec%Nwaves
   tab_delta_lambda(la) = tab_delta_lambda(la) - tab_delta_lambda(la-1) 
  end do
  if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
  allocate(tab_lambda_inf(n_lambda))
  if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
  allocate(tab_lambda_sup(n_lambda))
  tab_lambda_inf = tab_lambda
  tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
  ! computes stellar flux at the new wavelength points
  CALL deallocate_stellar_spectra()
  if (allocated(kappa)) deallocate(kappa)
  allocate(kappa(n_cells,n_lambda))
  kappa = 0.0 !Important to init !!
  !kappa will be computed on the fly in  optical_length_tot()
  !used for star map ray-tracing.
  CALL allocate_stellar_spectra(n_lambda)
  CALL repartition_energie_etoiles()
  ! Velocity field in  m.s-1
  if (.not.allocated(Vfield)) allocate(Vfield(n_cells))
  Vfield=atmos%Vmap !0 presently
  lkeplerian = .true.
  ! Warning : assume all stars are at the center of the disk
  if (.not.lVoronoi) then ! Velocities are defined from SPH files in Voronoi mode
     if (lcylindrical_rotation) then ! Midplane Keplerian velocity
        do icell=1, n_cells
           vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
        enddo
     else ! dependance en z
        do icell=1, n_cells
           vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg * r_grid(icell)**2 / &
                ((r_grid(icell)**2 + z_grid(icell)**2)**1.5 * AU_to_m) )
        enddo
     endif
  endif
  ! --- END ALLOCATING SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ----!

  
!! -------------------------------------------------------- !!
!! For NLTE do not forget ->
  !initiate NLTE popuplations for ACTIVE atoms, depending on the choice of the solution
  ! CALL initialSol()
  ! for now initialSol() is replaced by this if loop on active atoms
  if (atmos%Nactiveatoms.gt.0) then
    write(*,*) "solving for ", atmos%Nactiveatoms, " active atoms"
    do nact=1,atmos%Nactiveatoms
     write(*,*) "Setting initial solution for active atom ", atmos%ActiveAtoms(nact)%ID, &
      atmos%ActiveAtoms(nact)%active
     atmos%ActiveAtoms(nact)%n = 1d0 * atmos%ActiveAtoms(nact)%nstar
    end do
  end if !end replacing initSol()
!   ! Read collisional data and fill collisional matrix C(Nlevel**2) for each ACTIVE atoms.
!   ! Initialize at C=0.0 for each cell points.
!   ! the matrix is allocated for ACTIVE atoms only in setLTEcoefficients and the file open 
!   ! before the transfer starts and closed at the end.
!   do nact=1,atmos%Nactiveatoms 
!     if (atmos%ActiveAtoms(nact)%active) &
!      CALL CollisionRate(icell, atmos%ActiveAtoms(nact)) 
!   end do
! 
!! -------------------------------------------------------- !!


!  CALL ContributionFunction()

  write(*,*) "Computing emission flux map..."
  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
       CALL EMISSION_LINE_MAP(ibin,iaz)
     end do
  end do
  CALL WRITE_FLUX()

 ! Transfer ends, save data, free some space and leave
 do nact=1,atmos%Nactiveatoms
  CALL closeCollisionFile(atmos%ActiveAtoms(nact)) !if opened
 end do
 CALL freeAS() !deallocate opacity arrays
 CALL freeSpectrum() !deallocate spectral variables
 CALL free_atomic_atmos()  
 deallocate(ds)

 RETURN
 END SUBROUTINE

 SUBROUTINE WRITE_FLUX()

  character(len=512) :: filename !in case

  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  integer :: group,fpixel,nelements, i, xcenter
  logical :: simple, extend
  !integer :: a,b,c,d=1,e=1, idL
  real :: pixel_scale_x, pixel_scale_y 
npix_x = 101; npix_y = 101
  write(*,*) "Writing Flux-map"
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

   !------------------------------------------------------------------------------
   ! FLUX map
   !------------------------------------------------------------------------------
   bitpix=-64
   naxis=5
   !only squared pixels for now
!      if (RT_line_method==1) then
!         naxes(2)=1
!         naxes(3)=1
!      else
   naxes(1)=NLTEspec%Nwaves!1!1 if only one wavelength
   naxes(2)=npix_x
   naxes(3)=npix_y
!     endif
   naxes(4)=RT_n_incl
   naxes(5)=RT_n_az
   nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
  ! write(*,*) (naxes(i), i=1,naxis)

  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

!   !!RAC, DEC, reference pixel & pixel scale en degres
!   CALL ftpkys(unit,'CTYPE1',"RA---TAN",' ',status)
!   CALL ftpkye(unit,'CRVAL1',0.,-7,'RAD',status)
!   CALL ftpkyj(unit,'CRPIX1',npix_x/2+1,'',status)
!   pixel_scale_x = -map_size / (npix_x * distance * zoom) * arcsec_to_deg ! astronomy oriented (negative)
!   CALL ftpkye(unit,'CDELT1',pixel_scale_x,-7,'pixel scale x [deg]',status)
! 
!   CALL ftpkys(unit,'CTYPE2',"DEC--TAN",' ',status)
!   CALL ftpkye(unit,'CRVAL2',0.,-7,'DEC',status)
!   CALL ftpkyj(unit,'CRPIX2',npix_y/2+1,'',status)
!   pixel_scale_y = map_size / (npix_y * distance * zoom) * arcsec_to_deg
!   CALL ftpkye(unit,'CDELT2',pixel_scale_y,-7,'pixel scale y [deg]',status)

  CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'F_nu',status)
  
  if (l_sym_ima) then 
     !if (RT_line_method==1) then
     !squared pixels still
     !else
        xcenter = npix_x/2 + modulo(npix_x,2)
        do i=xcenter+1,npix_x
          NLTEspec%Flux(:,i,:,:,:) = NLTEspec%Flux(:,npix_x-i+1,:,:,:)
          NLTEspec%Fluxc(:,i,:,:,:) = NLTEspec%Fluxc(:,npix_x-i+1,:,:,:)
        end do
     !endif 
  endif ! l_sym_image

  !  Write the array to the FITS file.
!   do a=1,NLTEspec%Nwaves
!   do b=1,npix_x
!   do c=1,npix_y
!      if (NLTEspec%Flux(a,b,c,d,e)-1 == NLTEspec%Flux(a,b,c,d,e)) then 
!       write(*,*) "Infinite"
!       exit
!      end if
!      if (NLTEspec%Flux(a,b,c,d,e) /= NLTEspec%Flux(a,b,c,d,e)) then 
!       write(*,*) "Nan"
!       exit
!      end if
!   end do
!   end do
!   end do
  !idL = locate(NLTEspec%lambda,121.582d0)!121.568d0)
  !write(*,*) idL, NLTEspec%lambda(idL)
  !write(*,*) "max(II)=",MAXVAL(II), " min(II)",MINVAL(II)
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%Flux,status)
  
  ! create new hdu for continuum
  CALL ftcrhd(unit, status)

  !  Write the required header keywords.
  CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
  CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.pixel-1",'F_nu',status)
  CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%Fluxc,status)


  !  Close the file and free the unit number.
  CALL ftclos(unit, status)
  CALL ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     CALL print_error(status)
  endif

 RETURN
 END SUBROUTINE WRITE_FLUX
 
!  SUBROUTINE ContributionFunction()
!   integer :: icell, iray, id, NrecStokes, nbr_cell, previous_cell, next_cell,icellinf
!   double precision :: x0, y0, z0, norme, u0, v0, w0, x1, y1, z1
!   logical :: labs, lcellule_non_vide, lsubtract_avg, lonly_top, lonly_bottom
!   double precision :: ksi(NLTEspec%Nwaves,atmos%Nspace,atmos%Nrays,1), l
!                                                            !Replace 1 by 4 if pol.
!   double precision :: l_contrib, l_void_before, Ic(1,NLTEspec%Nwaves,atmos%Nrays)
!   double precision, dimension(NLTEspec%Nwaves) :: dtau, dtau_c, tau, tau_c, facteur_tau, Snu, Snu_c
!   ! for writing                                                        
!   integer :: status,unit,blocksize,bitpix,naxis
!   integer, dimension(7) :: naxes
!   integer :: group,fpixel,nelements, i, xcenter
!   logical :: simple, extend
!   character(len=512) :: CNTRB_FILE
!   CNTRB_FILE = "CNTRB.fits.gz"
! 
!   ! ------------------------------------------ !  
!   write(*,*) "Computing contribution function(s)..."
!   
!   NrecStokes = 1 !up to now
!   Ic = 0d0
!   dtau = 0d0
!   dtau_c = 0d0
!   tau = 0d0
!   tau_c = 0d0
!   nbr_cell = 0
!   Snu = 0d0
!   Snu_c = 0d0
!   id = 1 !sequentiel
!   do icell=1, n_cells
!    ! Propagation des rayons
!    do iray=1, atmos%Nrays
!     ! Position = milieu de la cellule
!     x0 = r_grid(icell)
!     y0 = 0.0_dp
!     z0 = z_grid(icell)
! 
!     norme = sqrt(x0*x0 + y0*y0 + z0*z0)
!     if (iray==1) then
!      u0 = x0/norme
!      v0 = y0/norme
!      w0 = z0/norme
!     else
!      u0 = -x0/norme
!      v0 = -y0/norme
!      w0 = -z0/norme
!     end if
!     ! Integration le long du rayon
!     ! Boucle infinie sur les cellules
!     next_cell = icell
!     infinie : do ! Boucle infinie
!      ! Indice de la cellule
!      icellinf = next_cell
!      x0=x1 ; y0=y1 ; z0=z1
!      if (icellinf <= n_cells) then
!         lcellule_non_vide=.true.
!      else
!         lcellule_non_vide=.false.
!      endif
!      
!      ! Test sortie
!      if (test_exit_grid(icellinf, x0, y0, z0)) then
!         exit infinie
!      endif
! 
!      nbr_cell = nbr_cell + 1
! 
!      ! Calcul longeur de vol et profondeur optique dans la cellule
!      previous_cell = 0 ! unused, just for Voronoi
!      call cross_cell(x0,y0,z0, u0,v0,w0,  icellinf, previous_cell, x1,y1,z1, next_cell, &
!                      l, l_contrib, l_void_before)
! 
!      if (lcellule_non_vide) then
!      lsubtract_avg = ((nbr_cell == 1).and.labs)
! 
!       ! opacities in m^-1
!       l_contrib = l_contrib * AU_to_m
! 
! 
!       CALL Background(1,icellinf, x0, y0, z0, u0, v0, w0)
!       dtau(:) =  l_contrib * (NLTEspec%ActiveSet%chi_c(:)+NLTEspec%ActiveSet%chi(:)) !scattering + thermal
!       dtau_c(:) = l_contrib * NLTEspec%ActiveSet%chi_c_bf(:)
! 
!       Snu = (NLTEspec%ActiveSet%eta_c + & 
!                   NLTEspec%ActiveSet%eta) / &
!                  (NLTEspec%ActiveSet%chi_c + NLTEspec%ActiveSet%chi)
!       ! continuum source function
!       Snu_c = (NLTEspec%ActiveSet%eta_c_bf) / NLTEspec%ActiveSet%chi_c_bf
! 
!       Ic(id,:,iray) = Ic(id,:,iray)*dexp(-dtau_c) + Snu_c * dexp(-dtau_c) * dtau_c
! 
!       facteur_tau = 1d0
!       if (lonly_top    .and. z0 < 0.) facteur_tau = 0d0
!       if (lonly_bottom .and. z0 > 0.) facteur_tau = 0d0
! 
!       tau = tau + dtau * facteur_tau
!       tau_c = tau_c + dtau_c
! 
!         ! set opacities to 0.0 for next cell point.
!         CALL initAS(re_init=.true.)
!      end if  ! lcellule_non_vide
!     end do infinie
!   !---------------------------------------------!    
!     ksi(:,icell,iray,1) = (Ic(id, :, iray) - Snu)*dexp(-tau)
!     
!    end do ! iray
!   end do ! icell
!   
!   ! ------------------------------------------ !
!   ! Now write them to disk
!   status=0
!   CALL ftgiou (unit,status)
! 
!   !  Create the new empty FITS file.
!   blocksize=1
!   CALL ftinit(unit,trim(CNTRB_FILE),blocksize,status)
! 
!   simple=.true.
!   extend=.false.
!   if (NrecStokes > 1) extend = .true.
!   group=1
!   fpixel=1  
!   
!   bitpix=-64
!   naxes(1)=NLTEspec%Nwaves
!   if (lVoronoi) then
!    if (NrecStokes > 1) then 
!     naxis = 4
!     naxes(2)=n_cells
!     naxes(3)=atmos%Nrays
!     naxes(4)=NrecStokes !if only Q NrecStokes = 2 for IQ etc
!     nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
!    else
!     naxis = 3 !only ksi_I
!     naxes(2)=n_cells
!     naxes(3)=atmos%Nrays
!     nelements=naxes(1)*naxes(2)*naxes(3)
!    end if
!   else
!    if (l3D) then
!     if (NrecStokes > 1) then 
!      naxis = 6
!      naxes(2) = n_rad
!      naxes(3) = 2*nz
!      naxes(4) = n_az
!      naxes(5)=atmos%Nrays
!      naxes(6)=NrecStokes !if only Q NrecStokes = 2 for IQ etc
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)*naxes(6)
!     else
!      naxis = 5
!      naxes(2) = n_rad
!      naxes(3) = 2*nz
!      naxes(4) = n_az
!      naxes(5)=atmos%Nrays
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
!     end if
!   else !not l3D
!     if (NrecStokes > 1) then 
!      naxis = 5
!      naxes(2) = n_rad
!      naxes(3) = nz
!      naxes(4)=atmos%Nrays
!      naxes(5)=NrecStokes !if only Q NrecStokes = 2 for IQ etc
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
!     else
!      naxis = 4
!      naxes(2) = n_rad
!      naxes(3) = nz
!      naxes(4)=atmos%Nrays
!      nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)
!     end if
!   end if !not l3D
!  end if !end if not Voronoi
! 
! 
!   !  Write the required header keywords.
!   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!   CALL ftpkys(unit,'BUNIT',"W.m-2.Hz-1.sr-1",'ksi',status)
! 
!   CALL ftpprd(unit,group,fpixel,nelements,ksi(:,:,:,1),status)
! !   if (NrecStokes > 1) then
! !    !write other ksi
! !   end if
! 
!   !  Close the file and free the unit number.
!   CALL ftclos(unit, status)
!   CALL ftfiou(unit, status)
! 
!   !  Check for any error, and if so print out error messages
!   if (status > 0) then
!      CALL print_error(status)
!   endif 
!  RETURN
!  END SUBROUTINE ContributionFunction

END MODULE AtomicTransfer