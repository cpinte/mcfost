! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
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
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

MODULE AtomicTransfer

 use metal, only                        : Background, BackgroundContinua, BackgroundLines
 use spectrum_type
 use atmos_type
 use readatom
 use lte
 use constant, only 					: MICRON_TO_NM
 use collision
 use solvene
 use writeatom
 use readatmos, only				    : readatmos_1D, readPLUTO
 use simple_models, only 				: spherically_symmetric_shells_model, &
 										  prop_law_model, uniform_law_model
 
 !$ use omp_lib
 
 !MCFOST's original modules
 use input
 use parametres
 use grid
 use optical_depth
 use dust_prop
 use dust_transfer, only 				: compute_stars_map
 use dust_ray_tracing, only 			: init_directions_ray_tracing ,            & 
                              			  tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az, &
                              			  tab_RT_incl, stars_map, kappa
 use stars
 use wavelengths
 use density
 
 IMPLICIT NONE

 CONTAINS

 
 SUBROUTINE INTEG_RAY_LINE(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
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
  ! what is labs? lsubtract_avg?
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  double precision, intent(in) :: u,v,w
  double precision, intent(in) :: x,y,z
  logical, intent(in) :: labs !used in NLTE but why?
  double precision :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  double precision, dimension(NLTEspec%Nwaves) :: Snu, Snu_c
  double precision, dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau
  integer :: nbr_cell, icell, next_cell, previous_cell
  double precision :: facteur_tau !used only in molecular line to have emission for one
                                  !  half of the disk only. Note used in AL-RT.
                                  !  this is passed through the lonly_top or lonly_bottom.
  logical :: lcellule_non_vide, lsubtract_avg

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

  ! -------------------------------------------------------------- !
  !*** propagation dans la grille ***!
  ! -------------------------------------------------------------- !


  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
    ! Indice de la cellule
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then ! .and. atmos%lcompute_atomRT(icell)
     lcellule_non_vide=.true.
    else
     lcellule_non_vide=.false.
    endif
       
    ! Test sortie ! "The ray has reach the end of the grid"
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

    if (lcellule_non_vide) then !count opacity only if the cell is filled, else go to next cell
     lsubtract_avg = ((nbr_cell == 1).and.labs)
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m !l_contrib in AU
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m

     CALL initAtomOpac(id) !set opac to 0 for this cell and thread id
     !! Compute background opacities for PASSIVE bound-bound and bound-free transitions
     !! at all wavelength points including vector fields in the bound-bound transitions
     if (NLTEspec%AtomOpac%store_opac) then
      CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
      dtau(:) =  l_contrib * (NLTEspec%AtomOpac%chi_p(id,:)+NLTEspec%AtomOpac%chi(id,:)+&
                    NLTEspec%AtomOpac%Kc(icell,:,1))
      dtau_c(:) = l_contrib * NLTEspec%AtomOpac%Kc(icell,:,1)


      Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
               NLTEspec%AtomOpac%eta_p(id,:) + NLTEspec%AtomOpac%eta(id,:) + &
                  NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%J(id,:)) / & 
                 (NLTEspec%AtomOpac%Kc(icell,:,1) + &
                 NLTEspec%AtomOpac%chi_p(id,:) + NLTEspec%AtomOpac%chi(id,:))

      Snu_c = (NLTEspec%AtomOpac%jc(icell,:) + & 
            NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%Jc(id,:)) / NLTEspec%AtomOpac%Kc(icell,:,1)
     else
      CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l) !x,y,z,u,v,w,x1,y1,z1 
                                !define the projection of the vector field (velocity, B...)
                                !at each spatial location.
      ! Epaisseur optique
      ! chi_p contains both thermal and continuum scattering extinction
      dtau(:) =  l_contrib * (NLTEspec%AtomOpac%chi_p(id,:)+NLTEspec%AtomOpac%chi(id,:))
      dtau_c(:) = l_contrib * NLTEspec%AtomOpac%chi_c(id,:)

      ! Source function
      ! No dust yet
      ! J and Jc are the mean radiation field for total and continuum intensities
      ! it multiplies the continuum scattering coefficient for isotropic (unpolarised)
      ! scattering. chi, eta are opacity and emissivity for ACTIVE lines.
      Snu = (NLTEspec%AtomOpac%eta_p(id,:) + NLTEspec%AtomOpac%eta(id,:) + &
                  NLTEspec%AtomOpac%sca_c(id,:) * NLTEspec%J(id,:)) / & 
                 (NLTEspec%AtomOpac%chi_p(id,:) + NLTEspec%AtomOpac%chi(id,:))
      ! continuum source function
      Snu_c = (NLTEspec%AtomOpac%eta_c(id,:) + & 
            NLTEspec%AtomOpac%sca_c(id,:) * NLTEspec%Jc(id,:)) / NLTEspec%AtomOpac%chi_c(id,:)
    end if
    !In ray-traced map (which is used this subroutine) we integrate from the observer
    ! to the source(s), but we actually get the outgoing intensity/flux. Direct
    !intgreation of the RTE from the observer to the source result in having the
    ! inward flux and intensity. Integration from the source to the observer is done
    ! when computing the mean intensity: Sum_ray I*exp(-dtau)+S*(1-exp(-dtau))
    NLTEspec%I(id,:,iray) = NLTEspec%I(id,:,iray) + & 
                             exp(-tau) * (1.0_dp - exp(-dtau)) * Snu
    NLTEspec%Ic(id,:,iray) = NLTEspec%Ic(id,:,iray) + &
                             exp(-tau_c) * (1.0_dp - exp(-dtau_c)) * Snu_c
!     NLTEspec%I(id,:,iray) = NLTEspec%I(id,:,iray) + dtau*Snu*dexp(-tau)
!     NLTEspec%Ic(id,:,iray) = NLTEspec%Ic(id,:,iray) + dtau_c*Snu_c*dexp(-tau_c)
 
    if (MINVAL(NLTEspec%I(id,:,iray))<0) then
     write(*,*) "Negative Opac !"
     if (NLTEspec%AtomOpac%store_opac) then
      write(*,*) id, icell, MINVAL(Snu), MINVAL(NLTEspec%AtomOpac%Kc(icell,:,1)) , MINVAL(NLTEspec%AtomOpac%Kc(icell,:,2)) , &
       MINVAL(NLTEspec%AtomOpac%jc(icell,:))
      stop
     else
      write(*,*) id, icell, MINVAL(Snu),&
      MINVAL(NLTEspec%AtomOpac%eta_p(id,:)),MINVAL(NLTEspec%AtomOpac%chi_p(id,:)),&
      MINVAL(NLTEspec%AtomOpac%sca_c(id,:)), MINVAL(NLTEspec%AtomOpac%eta_c(id,:)), MINVAL(NLTEspec%AtomOpac%chi_c(id,:))
      stop
     end if
    end if
!      !!surface superieure ou inf, not used with AL-RT
     facteur_tau = 1.0
!      if (lonly_top    .and. z0 < 0.) facteur_tau = 0.0
!      if (lonly_bottom .and. z0 > 0.) facteur_tau = 0.0

     ! Mise a jour profondeur optique pour cellule suivante
     ! dtau = chi * ds
     tau = tau + dtau * facteur_tau
     tau_c = tau_c + dtau_c

     ! Define PSI Operators here  
    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !

  
  RETURN
  END SUBROUTINE INTEG_RAY_LINE
 
  SUBROUTINE FLUX_PIXEL_LINE(&
         id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! -------------------------------------------------------------- !
   ! Computes the flux emerging out of a pixel.
   ! see: mol_transfer.f90/intensite_pixel_mol()
  ! -------------------------------------------------------------- !

   integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
   double precision, dimension(3), intent(in) :: pixelcorner,dx,dy
   double precision, intent(in) :: pixelsize,u,v,w
   integer, parameter :: maxSubPixels = 32
   double precision :: x0,y0,z0,u0,v0,w0
   double precision, dimension(NLTEspec%Nwaves) :: Iold, nu, I0, I0c
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
           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
             CALL INTEG_RAY_LINE(id, icell, x0,y0,z0,u0,v0,w0,iray,labs)
             I0 = I0 + NLTEspec%I(id,:,iray)
             I0c = I0c + NLTEspec%Ic(id,:,iray)
           !else !Outside the grid, no radiation flux
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

  nu = 1d0 !c_light / NLTEspec%lambda * 1d9 !to get W/m2 instead of W/m2/Hz !in Hz
  ! Flux out of a pixel in W/m2/Hz
  I0 = nu * I0 * (pixelsize / (distance*pc_to_AU) )**2
  I0c = nu * I0c * (pixelsize / (distance*pc_to_AU) )**2
  
  ! adding to the total flux map.
  if (RT_line_method==1) then
    NLTEspec%Flux(:,1,1,ibin,iaz) = NLTEspec%Flux(:,1,1,ibin,iaz) + I0
    NLTEspec%Fluxc(:,1,1,ibin,iaz) = NLTEspec%Fluxc(:,1,1,ibin,iaz) + I0c
  else
    NLTEspec%Flux(:,ipix,jpix,ibin,iaz) = I0
    NLTEspec%Fluxc(:,ipix,jpix,ibin,iaz) = I0c  
  end if

  RETURN
  END SUBROUTINE FLUX_PIXEL_LINE
  
 SUBROUTINE EMISSION_LINE_MAP(ibin,iaz)
 ! -------------------------------------------------------------- !
  ! Line emission map in a given direction n(ibin,iaz),
  ! using ray-tracing.
  ! if only one pixel it gives the total Flux.
  ! See: emission_line_map in mol_transfer.f90
 ! -------------------------------------------------------------- !
 
  integer, intent(in) :: ibin, iaz !define the direction in which the map is computed
  double precision :: x0,y0,z0,l,u,v,w

  double precision, dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  double precision, dimension(3,nb_proc) :: pixelcorner
  double precision:: taille_pix, nu
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 100, n_phi_RT = 36
  integer, parameter :: n_ray_star = 100 !1000
  double precision, dimension(n_rad_RT) :: tab_r
  double precision:: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda
  logical :: lresolved
  
  integer :: kc, kr, nat, Nred, Nblue, la
  logical :: lbottom_irradiated = .false. !to switch off the star in the map creation
  
npix_x = 101; npix_y = 101
  write(*,*) "incl (deg) = ", tab_RT_incl(ibin), "azimuth (deg) = ", tab_RT_az(iaz)

  u = tab_u_RT(ibin,iaz) ;  v = tab_v_RT(ibin,iaz) ;  w = tab_w_RT(ibin)
  uvw = (/u,v,w/) !vector position
  write(*,*) "u=", u, "v=",v, "w=",w

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

  if (RT_line_method==1) then !log pixels
    write(*,*) " WARNING: RT_line_method==1 not working correctly with AL-RT"
    n_iter_min = 1
    n_iter_max = 1
    
    dx(:) = 0.0_dp
    dy(:) = 0.0_dp
    i = 1
    j = 1
    lresolved = .false.

    rmin_RT = max(w*0.9_dp,0.05_dp) * Rmin
    rmax_RT = 2.0_dp * Rmax

    tab_r(1) = rmin_RT
    fact_r = exp( (1.0_dp/(real(n_rad_RT,kind=dp) -1))*log(rmax_RT/rmin_RT) )

    do ri_RT = 2, n_rad_RT
      tab_r(ri_RT) = tab_r(ri_RT-1) * fact_r
    enddo

    fact_A = sqrt(pi * (fact_r - 1.0_dp/fact_r)  / n_phi_RT )

    ! Boucle sur les rayons d'echantillonnage
    !$omp parallel &
    !$omp default(none) &
    !$omp private(ri_RT,id,r,taille_pix,phi_RT,phi,pixelcorner) &
    !$omp shared(tab_r,fact_A,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
    !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz)
    id =1 ! pour code sequentiel

    if (l_sym_ima) then
      cst_phi = pi  / real(n_phi_RT,kind=dp)
    else
      cst_phi = deux_pi  / real(n_phi_RT,kind=dp)
    endif

     !$omp do schedule(dynamic,1)
     do ri_RT=1, n_rad_RT
        !$ id = omp_get_thread_num() + 1

        r = tab_r(ri_RT)
        taille_pix =  fact_A * r ! racine carree de l'aire du pixel

        do phi_RT=1,n_phi_RT ! de 0 a pi
           phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)
           pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image
            ! C'est le centre en fait car dx = dy = 0.
           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, &
                      i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        end do !j
     end do !i
     !$omp end do
     !$omp end parallel  
  else !method 2  
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
     
  end if
 ! adding the stellar flux
 if (lbottom_irradiated) then
  write(*,*) " --> adding stellar flux map..."  
  if (etoile(1)%lb_body) then !all stars are Blackbodies ?
   !assumes a constant stellar flux across a transition.
   !computes flux star only at lambda0 and set Flux(Nred:Nblue) = Flux(lambda0)
   do nat=1,atmos%Natom
    do kc=1,atmos%Atoms(nat)%Ncont
     lambda = atmos%Atoms(nat)%continua(kc)%Nmid
     Nred = atmos%Atoms(nat)%continua(kc)%Nred
     Nblue = atmos%Atoms(nat)%continua(kc)%Nblue
     nu = c_light / NLTEspec%lambda(lambda) * 1d9 
     CALL compute_stars_map(lambda, u, v, w, taille_pix, dx, dy, lresolved)
      do la=Nblue,Nred
      NLTEspec%Flux(la,:,:,ibin,iaz) = NLTEspec%Flux(la,:,:,ibin,iaz) + stars_map(:,:,1) / nu
      NLTEspec%Fluxc(la,:,:,ibin,iaz) = NLTEspec%Fluxc(la,:,:,ibin,iaz) + stars_map(:,:,1) / nu
      end do
    end do
    do kr=1,atmos%Atoms(nat)%Nline
     lambda = atmos%Atoms(nat)%lines(kr)%Nmid
     Nred = atmos%Atoms(nat)%lines(kr)%Nred
     Nblue = atmos%Atoms(nat)%lines(kr)%Nblue
     nu = c_light / NLTEspec%lambda(lambda) * 1d9 
     CALL compute_stars_map(lambda, u, v, w, taille_pix, dx, dy, lresolved)
     !array = array + scalar
      do la=Nblue,Nred
      NLTEspec%Flux(la,:,:,ibin,iaz) = NLTEspec%Flux(la,:,:,ibin,iaz) + stars_map(:,:,1) / nu
      NLTEspec%Fluxc(la,:,:,ibin,iaz) = NLTEspec%Fluxc(la,:,:,ibin,iaz) + stars_map(:,:,1) / nu
      end do
    end do
   end do !over atoms
  else !a star is not a BB, computes for all wavelength
   do lambda = 1, NLTEspec%Nwaves
    nu = c_light / NLTEspec%lambda(lambda) * 1d9 !if NLTEspec%Flux in W/m2 set nu = 1d0 Hz
                                             !else it means that in FLUX_PIXEL_LINE, nu
                                             !is 1d0 (to have flux in W/m2/Hz)
    CALL compute_stars_map(lambda, u, v, w, taille_pix, dx, dy, lresolved)
!      write(*,*) "Stellar flux at lambda = ", NLTEspec%lambda(lambda), &
!       stars_map(npix_x/2+1,npix_y/2+1,1)/nu, maxval(NLTEspec%Flux(lambda,:,:,ibin,iaz)), &
!       maxval(stars_map(:,:,1))/nu

     NLTEspec%Flux(lambda,:,:,ibin,iaz) = NLTEspec%Flux(lambda,:,:,ibin,iaz) +  & 
                                         stars_map(:,:,1) / nu
     NLTEspec%Fluxc(lambda,:,:,ibin,iaz) = NLTEspec%Fluxc(lambda,:,:,ibin,iaz) + &
                                         stars_map(:,:,1) / nu

   end do  
  end if  
 end if

 RETURN
 END SUBROUTINE EMISSION_LINE_MAP
 
 ! NOTE: should inverse the order of frequencies and depth because in general
 !       n_cells >> n_lambda, in "real" cases.
 ! npix_x, xpix_y, RT_line_method
 ! some of shared quantities by the code that i don't know were they are !!
 ! lkeplerian, linfall, lvoro where to declare them, lcylindrical ects
 
 SUBROUTINE Atomic_transfer()
 ! --------------------------------------------------------------------------- !
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
 ! --------------------------------------------------------------------------- !

#include "sprng_f.h"

  integer :: atomunit = 1, nact,npass,  nat, la !atoms and wavelength
  integer :: icell !spatial variables
  integer :: ibin, iaz!, RT_line_method
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! the following quantities are input parameters !!
  integer :: NiterMax = 20, Nrays = 1! Number of rays for angular integration and to compute Inu(mu)
  integer :: IterLimit
  logical :: lsolve_for_ne = .false. !for calculation of electron density even if atmos%calc_ne
                                     ! is .false.
  character(len=7) :: NE0 = "HIONIS" !set this parameter according to the temperature.
  logical :: lstore_opac = .true.
  logical :: lstatic !to change the default of atmos%moving (.true.) in .false.
  					 ! should be considered only for testing purposes
  					 !or special geometry with static atmospheres
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  optical_length_tot => atom_optical_length_tot
  !because for now lemission_atom is not a case of readparameters
  write(*,*) npix_x, npix_y
if ((npix_x /= 101).or.(npix_y /= 101)) write(*,*) 'BEWARE: npix_x read is different from what it should be..'
npix_x = 101; npix_y = 101

!! ---------------- ---------------------------------------- !!
  
  atmos%Nrays = Nrays !remember, it is needed to allocate I, Ic
  if (atmos%Nrays == 0) then
   write(*,*) "Nrays should at least be 1!"
   stop
  end if

!! ----------------------- Read Model ---------------------- !!
  !CALL uniform_law_model()
  !CALL prop_law_model()
  CALL spherically_symmetric_shells_model()

  ! OR READ FROM MODEL (to move elsewhere) 
  !suppose the model is in utils/Atmos/
  !CALL readatmos_1D("Atmos/FALC_mcfost.fits.gz")
  
  write(*,*) "maxTgas = ", MAXVAL(atmos%T), " minTgas = ", MINVAL(atmos%T)
  write(*,*) "maxnH = ", MAXVAL(atmos%nHtot), " minnH = ", MINVAL(atmos%nHtot)  
!! --------------------------------------------------------- !!

  !Read atomic models and allocate space for n, nstar, vbroad, ntotal, Rij, Rji
  ! on the whole grid space.
  ! The following routines have to be invoked in the right order !
  CALL readAtomicModels(atomunit)
  !we need atmos%v_char here for computing wavelength grid but it is known far after,
  !except if we read velocity from a model. But we do not allow for a depth dependent
  !frequency grid.
  !If velocity is keplerian, vchar is the maximum value of Keplerian velocity

  ! if the electron density is not provided by the model or simply want to
  ! recompute it
  ! if atmos%ne given, but one want to recompute ne
  ! else atmos%ne not given, its computation is mandatory
  ! if NLTE pops are read in readAtomicModels, H%n(Nlevel,:) and atom%n
  !can be used for the electron density calculation.
  if (.not.atmos%calc_ne) atmos%calc_ne = lsolve_for_ne
  if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
  if (atmos%calc_ne) CALL SolveElectronDensity(atmos%ne,NE0)
  CALL writeElectron() !will be moved elsewhere
  ! do it in the reading process
  CALL writeHydrogenDensity()  

  Call setLTEcoefficients () !write pops at the end because we probably have NLTE pops also

  write(*,*) "maxNE = ", MAXVAL(atmos%ne), " minNE = ", MINVAL(atmos%ne)
  write(*,*) "maxnHmin = ", MAXVAL(atmos%nHmin), " minnHmin = ", MINVAL(atmos%nHmin)

!! --------------------------------------------------------- !!
  NLTEspec%atmos => atmos !this one is important because atmos_Type : atmos is not used.
  CALL initSpectrum(nb_proc, 500d0)  !optional vacuum2air and writewavelength
  NLTEspec%AtomOpac%store_opac = lstore_opac
  CALL allocSpectrum(npix_x, npix_y, RT_n_incl, RT_n_az)
  if (NLTEspec%AtomOpac%store_opac) then !only Background lines and active transitions
                                         ! chi and eta, are computed on the fly.
                                         ! Background continua (=ETL) are kept in memory.
   !$omp parallel &
   !$omp default(none) &
   !$omp private(icell) &
   !$omp shared(atmos) 
   !$omp do
   do icell=1,atmos%Nspace
    CALL BackgroundContinua(icell)
   end do
   !$omp end do
   !$omp end parallel
  end if

  ! ----- ALLOCATE SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ------!
  CALL adjusting_mcfost()
  ! --- END ALLOCATING SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ----!
!! --------------------- NLTE--------------------------------- !!
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
!     if (atmos%ActiveAtoms(nact)%active) then
!      CALL CollisionRate(icell, atmos%ActiveAtoms(nact)) 
!      CALL initGamma(icell, atmos%ActiveAtoms(nact)) !set Gamma to C for each active atom
!       !note that when updating populations, if ne is kept fixed (and T and nHtot etc)
!       !atom%C is fixed, therefore we only use initGamma. If they chane, call CollisionRate() again
!      end if
!   end do
!   if (atmos%Nrays==1) then !placed before computing J etc
!    write(*,*) "Solving for", atmos%Nrays,' direction.'
!   else
!    write(*,*) "Solving for", atmos%Nrays,' directions.'
!   end if
! 
!! -------------------------------------------------------- !!


  write(*,*) "Voronoi?",lVoronoi, " keplerian?",lkeplerian, &
    " cylindrical?", lcylindrical_rotation, " moving?", atmos%moving

  write(*,*) "Computing emission flux map..."
  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
       CALL EMISSION_LINE_MAP(ibin,iaz)
     end do
  end do
  CALL WRITE_FLUX()
  
!   open(10*atomunit, file="Snu_0.dat",action="write",status="unknown")
!   CALL initAtomOpac(1)
!   CALL Background(1, 1, 0d0,0d0, 0d0, 0d0, 0d0,0d0, 0d0,0d0,0d0)
!   do la=1,NLTEspec%Nwaves
!    ! write(*,*) NLTEspec%lambda(la), (NLTEspec%AtomOpac%eta_p(1,la)) / (NLTEspec%AtomOpac%chi_p(1,la)) 
!    write(10*atomunit,*) NLTEspec%lambda(la), (NLTEspec%AtomOpac%eta_p(1,la)) / & 
!     (NLTEspec%AtomOpac%chi_p(1,la)) 
!   end do
!   close(10*atomunit)

 ! Transfer ends, save data, free some space and leave
 do nact=1,atmos%Nactiveatoms
  CALL closeCollisionFile(atmos%ActiveAtoms(nact)) !if opened
 end do
 CALL freeSpectrum() !deallocate spectral variables
 CALL free_atomic_atmos()  
 deallocate(ds)
 NULLIFY(optical_length_tot)

 RETURN
 END SUBROUTINE
 
 SUBROUTINE adjusting_mcfost()
  !--> should move them to init_atomic_atmos ? or elsewhere
  !need to be deallocated at the end of molecule RT or its okey ?`
  integer :: icell, la
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
  tab_lambda = NLTEspec%lambda / MICRON_TO_NM
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
  if (allocated(kappa)) deallocate(kappa) !do not use it anymore
!   allocate(kappa(n_cells,n_lambda))
!   kappa = 0.0 !Important to init !!
!   !kappa will be computed on the fly in  optical_length_tot()
  !used for star map ray-tracing.
  CALL allocate_stellar_spectra(n_lambda)
  CALL repartition_energie_etoiles()
  ! Velocity field in  m.s-1
  if (.not.allocated(Vfield)) allocate(Vfield(n_cells))
  lkeplerian = .false.
  linfall = .false.
  Vfield=0d0! by default or something if velocity_law != 0
  if (atmos%moving) then
   if (atmos%velocity_law == 0) then !default is 0 for keplerian vel
    lkeplerian = .true.
    ! Warning : assume all stars are at the center of the disk
    if (.not.lVoronoi) then ! Velocities are defined from SPH files in Voronoi mode
     if (lcylindrical_rotation) then ! Midplane Keplerian velocity
        do icell=1, n_cells
           vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg /  (r_grid(icell) * AU_to_m) )
        end do
     else ! dependance en z
        do icell=1, n_cells
           vfield(icell) = sqrt(Ggrav * sum(etoile%M) * Msun_to_kg * r_grid(icell)**2 / &
                ((r_grid(icell)**2 + z_grid(icell)**2)**1.5 * AU_to_m) )
        end do
     end if
    end if
   else if (atmos%velocity_law == -1) then 
    linfall=.true.
    Vfield = atmos%Vmap
!      write(*,*) linfall, atmos%velocity_law, MAXVAL(Vfield)
!      stop
   end if
   if (allocated(atmos%Vmap)) deallocate(atmos%Vmap)   !free Vmap here because we do not use it	
  else
   deallocate(Vfield)						! allocated only if moving
   !vchar cannot be defined here because wavelength grid is already set up..
  end if 

 RETURN
 END SUBROUTINE adjusting_mcfost

END MODULE AtomicTransfer