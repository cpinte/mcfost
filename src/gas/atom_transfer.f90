! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atomic (only) systems using the MALI method
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
module atom_transfer

   use parametres
   use input, only               : lkeplerian, linfall, RT_line_method, RT_line_method, &
                                 limb_darkening, mu_limb_darkening
   use constantes, only : nm_to_m, m_to_km, au_to_m, deg_to_rad, tiny_real, tiny_dp, pi, deux_pi, pc_to_au
   use io_atom, only : read_atomic_models, write_pops_atom
   use wavelengths, only : n_lambda, tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda
   use wavelengths_gas, only : make_wavelength_nlte, tab_lambda_nm, tab_lambda_cont, n_lambda_cont
   use elecdensity, only : solve_ne, write_electron, read_electron
   use grid, only : lcalc_ne, move_to_grid, vfield3d, icompute_atomRT
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops
   use atom_type, only : atoms, atomtype, n_atoms
   use init_mcfost, only :  nb_proc
   use gas_contopac, only : background_continua_lambda
   use opacity_atom, only : alloc_atom_opac, Itot, psi, opacity_atom_bf_loc, opacity_atom_bb_loc
   use optical_depth, only : integ_ray_atom
   use utils, only : cross_product, gauss_legendre_quadrature, progress_bar, rotation_3d
   use dust_ray_tracing, only    : RT_n_incl, RT_n_az, init_directions_ray_tracing,tab_u_RT, tab_v_RT, tab_w_RT, &
                                   tab_RT_az,tab_RT_incl
   use stars, only               : intersect_stars, laccretion_shock, max_Tshock, min_Tshock
   use output, only : allocate_atom_maps, flux_total, write_total_flux, write_atomic_maps

   !$ use omp_lib

   implicit none
                    
   integer :: omp_chunk_size
   real(kind=dp) :: dne

   contains

   subroutine nlte_loop_mali()

      return
   end subroutine nlte_loop_mali

   subroutine setup_image_grid()

      return
   end subroutine setup_image_grid

   subroutine atom_line_transfer()
   ! --------------------------------------------------------------------------- !
   ! This routine initialises the necessary quantities for atomic line transfer
   ! and calls the appropriate routines for LTE or NLTE transfer.
   ! --------------------------------------------------------------------------- !
      integer  :: ne_initial, ibin, iaz, nat
      logical :: lelectron_read
      real(kind=dp) :: v_char, max_vel_shift, v1, v2


      write(*,*) " *** BIG WARNING **** "
      write(*,*) " !!!!! CHECK FOR BILINEAR INTERP in utils.f90 ! !!!!"
      write(*,*) "check solve ne" 
      write(*,*) " ******************** "

      omp_chunk_size = max(nint( 0.01 * n_cells / nb_proc ),1)
      mem_alloc_tot = 0
      v_char = sqrt( maxval(sum(vfield3d**2,dim=2)) )

      !read atomic models
      call read_atomic_Models()

      !is there an electron density file in the current directory ?
      call read_electron(lelectron_read)
      if (lelectron_read) then
         !yes, use that as initial solution
         ne_initial = 1
         lcalc_ne = lsolve_for_ne
      else
         !no, use hydrogen ionisation
         ne_initial = 0
         if (lcalc_ne) then
            write(*,*) "Solving for electron density"
            write(*,*) " -> initial solution from H+M ionisation."
         endif
      endif

      !no electron density in the model nor previous fits file calc_ne == True.
      !if a previous file is found (lelectron_read == True) and lsolve_for_ne, then
      !calc_ne is set to .true. Electron density is re-evaluated using the populations
      !in the fits files (H_Ionisation otherwise).
      if (lcalc_ne) then !eventually computed with the previous non-LTE pops !
         call Solve_ne(ne_initial, .true., dne)
         call write_Electron
      else
         if (lsolve_for_ne) then
            write(*,*) "Forcing calculation of electron density"
            write(*,*) " -> initial solution from model."
            call Solve_ne(ne_initial, .true., dne)
            call write_Electron
         endif
      endif

      call ltepops_atoms() 
      max_vel_shift = 0.0
      do ibin=1,n_cells
         v1 = sqrt(sum(vfield3d(ibin,:)**2))
         do iaz=1,n_cells
            v2 = sqrt(sum(vfield3d(iaz,:)**2))
            if ((icompute_atomRT(ibin)>0).and.(icompute_atomRT(iaz)>0)) then
               max_vel_shift = max(max_vel_shift,abs(v1-v2))
            endif
         enddo
      enddo
      write(*,*) max_vel_shift*1d-3, maxval(sqrt(sum(vfield3d(:,:)**2,dim=2)))*1d-3
      v_char = max_vel_shift
      !because for non-LTE we compute the motion of each cell wrt to the actual cell
      !so the maximum shift is what we have if all cells move at their max speed (v(icell))
      ! stop

      !every indexes are defined on the lambda frequency. 
      !So even when lambda is passed as an argument, it is expected
      !that the indexes correspond at that grid.
      !Even if we split the transfer un group of wavelength that's probably best to keep indexes.
      !otherwise, we pass group wavelength and check contributing opacities for each group
      !but it is like indexes at the end?
      deallocate(tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda)
      write(*,*) "v_char", v_char * 1d-3
      call make_wavelength_nlte(tab_lambda_nm,vmax_overlap=v_char)
      n_lambda = size(tab_lambda_nm)
      tab_lambda = tab_lambda_nm * nm_to_m * m_to_km !micron

      !allocate quantities in space and for this frequency grid
      call alloc_atom_opac(n_lambda, tab_lambda_nm)
      !-> still allocate cont on tab_lambda_cont 
      !-> use interpolation to go from cont grid to total grid
      !except for images ?

      !for image keep only the transitions in the file for ray_tracing
      !recompute a special wavelength grid with new n_lambda and tab_lambda
      !and solve the transfer for these wavelengths covering only the transitions
      !for ray tracing.
      !keep a copy of tab_lambda from the file if we want to compute a sed
      !check that the index on the lambda for map is consistent witht eh value of lambda for a 
      !line because in worst, line%nr and line%nb can be 1, but very unlikely meaning the line is on the grid.

      if (lmodel_1d) then
         !in 1d we don't care about the images, we compute cylindrically symmetric
         !intensity (still with velocity).
         call spectrum_1d()
         !deallocate and exit code
         return !from atomic transfer!
      endif


      !re alloc lambda here. Eventually use the tab_lambda in the file
      !for SED flux
      !add an limage mode
      !call setup_image_grid()!either for RT_line_method==1 or >1
      call allocate_atom_maps()
      write(*,*) "Computing emission flux map..."
      if (laccretion_shock) then
         max_Tshock = 0.0
         min_Tshock = 1d8
      endif
      call init_directions_ray_tracing()
      do ibin=1,RT_n_incl
         do iaz=1,RT_n_az
            call emission_line_map(ibin,iaz)
         end do
      end do
      do nat=1, n_atoms
         if (atoms(nat)%p%lline) call write_atomic_maps(atoms(nat)%p) !onyly with RT = 2
      enddo
      call write_total_flux() !tmp only with RT = 1
      write(*,*) " ..done"


      return
   end subroutine atom_line_transfer


  subroutine intensite_pixel_atom(id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
   ! -------------------------------------------------------------- !
   ! stellar map ?
   ! TO DO line by line
   ! Magnetic fields
   ! -------------------------------------------------------------- !

      integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
      real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
      real(kind=dp), intent(in) :: pixelsize,u,v,w
      integer, parameter :: maxSubPixels = 32
      real(kind=dp) :: x0,y0,z0,u0,v0,w0
      real(kind=dp), dimension(N_lambda) :: Iold, I0
      real(kind=dp), dimension(3) :: sdx, sdy
      real(kind=dp):: npix2, diff, normF
      real(kind=dp), parameter :: precision = 1.e-2
      integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter
      logical :: lintersect, labs
      integer :: kr, krr, nat
      type (AtomType), pointer :: atom

      labs = .false.
      u0 = -u ; v0 = -v ; w0 = -w

      iray = 1 
      ri = 2*n_rad ; zj=1 ; phik=1

      ! le nbre de subpixel en x est 2^(iter-1)
      subpixels = 1
      iter = 1
      diff = 0.
      I0 = 0.0_dp

      infinie : do ! Boucle infinie tant que le pixel n'est pas converge
         npix2 =  real(subpixels)**2
         Iold = I0
         I0 = 0.0_dp
         ! Vecteurs definissant les sous-pixels
         sdx(:) = dx(:) / real(subpixels,kind=dp)
         sdy(:) = dy(:) / real(subpixels,kind=dp)

       !      !L'obs est en dehors de la grille
       !      ri = 2*n_rad ; zj=1 ; phik=1

         ! Boucle sur les sous-pixels qui calcule l'intensite au centre
         ! de chaque sous pixel
         do i = 1,subpixels
            do j = 1,subpixels
            ! Centre du sous-pixel
               x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
               y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
               z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)
               ! On se met au bord de la grille : propagation a l'envers
               call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
               if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
                  call integ_ray_atom(id,icell,x0,y0,z0,u0,v0,w0,iray,labs,n_lambda,tab_lambda_nm)

                  I0 = I0 + Itot(:,iray,id)

                 !else !Outside the grid, no radiation flux
               endif
            end do !j
         end do !i

         I0 = I0 / npix2

         if (iter < n_iter_min) then
            ! On itere par defaut
            subpixels = subpixels * 2
         else if (iter >= n_iter_max) then
            ! On arrete pour pas tourner dans le vide
            exit infinie
         else
            ! On fait le test sur a difference
            diff = maxval( abs(I0 - Iold) / (I0 + 1e-300_dp) )
            ! There is no iteration for Q, U, V, assuming that if I is converged, then Q, U, V also.
            ! Can be added and then use diff with I = (sqrt(I**2 + Q**2 + U**2 + V**2))
            if (diff > precision ) then
               ! On est pas converge
               subpixels = subpixels * 2
            else
             !write(*,*) "Pixel converged", ipix, jpix, i, j, iter, diff
               exit infinie
            end if
         end if ! iter
         iter = iter + 1
      end do infinie

      ! Prise en compte de la surface du pixel (en sr)

      ! Flux out of a pixel in W/m2/Hz/pix
      normF = ( pixelsize / (distance*pc_to_AU) )**2

      !-> temporary becaus eI haven't redefined the grid for images!
      ! if (RT_line_method==1) then
         Flux_total(:,ibin,iaz,id) = Flux_total(:,ibin,iaz,id) + I0(:) * normF
      ! endif


      if (RT_line_method==2) then
         do nat=1,N_atoms
            atom => atoms(nat)%p
            do kr=1,atom%nTrans_rayTracing
               krr = atom%ij_to_trans(atom%i_Trans_rayTracing(kr),atom%j_Trans_rayTracing(kr))
               atom%lines(krr)%map(ipix,jpix,:,ibin,iaz) = I0(atom%lines(krr)%nb:atom%lines(krr)%nr)*normF
            enddo
            atom => NULL()
         enddo
      endif   


      return 
   end subroutine intensite_pixel_atom 

  subroutine emission_line_map(ibin,iaz)
   ! -------------------------------------------------------------- !
   !
   ! Line emission map in a given direction n(ibin,iaz),
   ! using ray-tracing.
   ! if only one pixel it gives the total Flux.
   ! See: emission_line_map in mol_transfer.f90
   !
   ! -------------------------------------------------------------- !
      integer, intent(in) :: ibin, iaz
      real(kind=dp) :: x0,y0,z0,l,u,v,w
      real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
      real(kind=dp), dimension(3,nb_proc) :: pixelcorner
      real(kind=dp):: taille_pix
      integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max
      integer, parameter :: n_rad_RT = 150, n_phi_RT = 360
      real(kind=dp), dimension(n_rad_RT) :: tab_r
      real(kind=dp):: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
      integer :: ri_RT, phi_RT
      logical :: lresolved

      write(*,*) "Vector to observer =", real(tab_u_rt(ibin,iaz)),real(tab_v_rt(ibin,iaz)),real(tab_w_rt(ibin))
      write(*,*) "i=", real(tab_RT_incl(ibin)), "az=", real(tab_RT_az(iaz))

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

      write(*,*) "x-image =           ", real(x_plan_image(:))
      write(*,*) "y-image =           ", real(y_plan_image(:))

      ! position initiale hors modele (du cote de l'observateur)
      ! = centre de l'image
      l = 10.*Rmax  ! on se met loin ! in AU

      x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
      center(1) = x0 ; center(2) = y0 ; center(3) = z0

      ! Coin en bas gauche de l'image
      Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)

      if (RT_line_method==1) then !log pixels

         !no sub pixels here
         n_iter_min = 1
         n_iter_max = 1

         dx(:) = 0.0_dp
         dy(:) = 0.0_dp

         i = 1
         j = 1
         lresolved = .false.

         rmin_RT = 0.001_dp * Rmin!max(w*0.9_dp,0.05_dp) * Rmin
         rmax_RT = Rmax * 1.0_dp ! Rmax  * 2.0_dp

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
         !$omp shared(tab_r,fact_A,x,x_plan_image,y_plan_image,center,dx,dy,u,v,w,i,j) &
         !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz,fact_r)
         id = 1 ! pour code sequentiel

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
            do phi_RT=1,n_phi_RT ! de 0 + eps à 2pi - eps (eps = pi/n_phi_RT)
               phi = cst_phi * (real(phi_RT,kind=dp) -0.5_dp)

               pixelcorner(:,id) = center(:) + r * sin(phi) * x_plan_image + r * cos(phi) * y_plan_image
               ! C'est le centre en fait car dx = dy = 0.
               call intensite_pixel_atom(id,ibin,iaz,n_iter_min,n_iter_max, i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
            end do !j
         end do !i
         !$omp end do
         !$omp end parallel

         !-> if star_map set dx and dy to
         !!taille_pix = (map_size/zoom)  ! en AU
         !!dx(:) = x_plan_image * taille_pix
         !!dy(:) = y_plan_image * taille_pix


      else !method 2
         ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
         taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
         lresolved = .true.

         dx(:) = x_plan_image * taille_pix
         dy(:) = y_plan_image * taille_pix

         !-> dust_transfer
         !!Icorner(:) = center(:) - ( 0.5 * npix_x * dx(:) +  0.5 * npix_y * dy(:))

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
         n_iter_min = 1 !1 !3
         n_iter_max = 1 !1 !6
         !$omp do schedule(dynamic,1)
         do i = 1,npix_x_max
            !$ id = omp_get_thread_num() + 1
            do j = 1,npix_y
               pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
               call intensite_pixel_atom(id,ibin,iaz,n_iter_min,n_iter_max,i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
            end do !j
         end do !i
         !$omp end do
         !$omp end parallel
      end if

      return
   end subroutine emission_line_map   

   subroutine spectrum_1d()
      !TO DO log rays
      !TO DO fits format

      !NOTE: cos_theta in that case is minimum at Rmax which is not the limb.
      ! the limb is at p = 1.0 (rr=Rstar). In plan geomtry, cos_theta is minimum at the limb.
      ! the plan-parralel cos_theta is np.sqrt(1.0 - p**2) for p < 1.0.
      integer :: la, j, icell0, id
      logical :: lintersect, labs
      integer, parameter :: Nimpact = 10
      real(kind=dp) :: rr, u,v,w,u0,w0,v0,x0,y0,z0,x(3),y(3),uvw(3)
      real(kind=dp), allocatable :: cos_theta(:), weight_mu(:), p(:)
      real(kind=dp), allocatable ::I_1d(:,:)
      integer :: n_rays_done, ibar

      write(*,*) "**** computing CLV intensity for 1d model..."

      allocate(cos_theta(Nimpact), weight_mu(Nimpact), p(Nimpact))
      allocate(I_1d(n_lambda,Nimpact))

      deallocate(Itot)
      allocate(Itot(N_lambda,Nimpact,nb_proc)); Itot = 0.0_dp

      u = 0.0_dp
      v = 0.0_dp
      w = 1.0_dp

      id = 1
      n_rays_done = 0
      ibar = 0
      labs = .false.

      uvw = (/u,v,w/) !vector position
      x = (/1.0_dp,0.0_dp,0.0_dp/)
      y = -cross_product(x, uvw)

      ! Ray tracing : on se propage dans l'autre sens
      u0 = -u ; v0 = -v ; w0 = -w
      call gauss_legendre_quadrature(0.0_dp, 1.0_dp, Nimpact, cos_theta, weight_mu)
      I_1d = 0.0
      call progress_bar(0)
      !$omp parallel &
      !$omp default(none) &
      !$omp private(id,icell0,j) &
      !$omp private(x0,y0,z0,lintersect,rr) &
      !$omp shared(cos_theta,p,Itot,N_lambda,tab_lambda_nm,x,y,weight_mu,ibar,n_rays_done)&
      !$omp shared(u,v,w,u0,v0,w0,Rmax,Rmin,labs)
      !$omp do schedule(dynamic,1)
      do j=1,Nimpact
         !$ id = omp_get_thread_num() + 1

         rr = Rmax * sqrt(1.0 - cos_theta(j)**2)
         p(j) = rr/Rmin

         x0 = 10.0*Rmax*u + rr * y(1)
         y0 = 10.0*Rmax*v + rr * y(2)
         z0 = 10.0*Rmax*w + rr * y(3)
         
         call move_to_grid(id,x0,y0,z0,u0,v0,w0,icell0,lintersect)
         if (lintersect) then
            call integ_ray_atom(id,icell0,x0,y0,z0,u0,v0,w0,j,labs,n_lambda,tab_lambda_nm)
         endif

         !$omp atomic
         n_rays_done = n_rays_done + 1
         if (real(n_rays_done) > 0.02*ibar*Nimpact) then
            call progress_bar(ibar)
            !$omp atomic
            ibar = ibar+1
         endif  
      enddo
      !$omp end do
      !$omp end parallel
      call progress_bar(50)
      I_1d(:,:) = Itot(:,:,1)
      do j=2,nb_proc !direct sum on nb_proc raises a segmentation fault
         I_1d(:,:) = I_1d(:,:) + Itot(:,:,j)
      enddo

      open(1,file="spectrum_1d.txt",status="unknown")
      write(1,*) Nimpact, N_lambda
      write(1,'(*(1ES17.8E3))') (p(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (cos_theta(j), j=1,Nimpact)
      write(1,'(*(1ES17.8E3))') (weight_mu(j), j=1,Nimpact)
      write(1,'(*(1F12.4))') (tab_lambda_nm(la), la=1,N_lambda)
      do la=1,N_lambda
            write(1,'(*(1ES17.8E3))') (I_1d(la,j), j=1,Nimpact)
      enddo
      close(1)

      deallocate(cos_theta, weight_mu, p, I_1d, Itot)

   return
   end subroutine spectrum_1d

end module atom_transfer