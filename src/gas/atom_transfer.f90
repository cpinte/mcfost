! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atomic (only) systems using the MALI method
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
module atom_transfer

   use parametres
   use constantes, only : km_to_m, m_to_km
   use io_atom, only : read_atomic_models, write_pops_atom
   use wavelengths, only : n_lambda, tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda
   use wavelengths_gas, only : make_wavelength_nlte
   use elecdensity, only : solve_ne, write_electron, read_electron
   use grid, only : lcalc_ne, move_to_grid
   use lte, only : ltepops_atoms, ltepops_atoms_1, print_pops
   use atom_type, only : atoms
   use init_mcfost, only :  nb_proc
   use gas_contopac, only : background_continua_lambda
   use opacity_atom, only : alloc_atom_opac, Itot, psi, opacity_atom_bf_loc, opacity_atom_bb_loc
   use optical_depth, only : integ_ray_atom
   use utils, only : cross_product, gauss_legendre_quadrature

   implicit none
                    
   integer :: omp_chunk_size

   contains



   subroutine atom_line_transfer()
   ! --------------------------------------------------------------------------- !
   ! This routine initialises the necessary quantities for atomic line transfer
   ! and calls the appropriate routines for LTE or NLTE transfer.
   ! --------------------------------------------------------------------------- !
      integer  :: ne_initial
      logical :: lelectron_read
      real(kind=dp) :: dne

      integer :: la, j, icell0
      logical :: lintersect
      real(kind=dp) :: rr, u,v,w,u0,w0,v0,x0,y0,z0,x(3),y(3),uvw(3),cos_theta(1), weight_mu(1)
      real(kind=dp), allocatable :: chit(:),etat(:),xlam(:), Ftot(:,:)

      omp_chunk_size = max(nint( 0.01 * n_cells / nb_proc ),1)
      mem_alloc_tot = 0

      !read atomic models
      call read_atomic_Models()
      !every indexes are defined on the lambda frequency. 
      !So even when lambda is passed as an argument, it is expected
      !that the indexes correspond at that grid.
      deallocate(tab_lambda, tab_lambda_inf, tab_lambda_sup, tab_delta_lambda)
      call make_wavelength_nlte(tab_lambda)
      n_lambda = size(tab_lambda)
      !-> the lines above can be mose elsewhere to that everything is defined in that grid
      !except the dust opacities that can be interpolated locally

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

      !allocate quantities in space and for this frequency grid
      call alloc_atom_opac(n_lambda, km_to_m*tab_lambda)

      !for image keep only the transitions in the file for ray_tracing
      !recompute a special wavelength grid with new n_lambda and tab_lambda
      !and solve the transfer for these wavelengths covering only the transitions
      !for ray tracing.
      !keep a copy of tab_lambda from the file if we want to compute a sed
      !check that the index on the lambda for map is consistent witht eh value of lambda for a 
      !line because in worst, line%nr and line%nb can be 1, but very unlikely meaning the line is on the grid.

      u = 0.0_dp
      v = 0.0_dp
      w = 1.0_dp

      uvw = (/u,v,w/) !vector position
      x = (/1.0_dp,0.0_dp,0.0_dp/)
      y = -cross_product(x, uvw)

      ! Ray tracing : on se propage dans l'autre sens
      u0 = -u ; v0 = -v ; w0 = -w
      call gauss_legendre_quadrature(0.0_dp, 1.0_dp, size(cos_theta), cos_theta, weight_mu)
      allocate(Ftot(n_lambda,nb_proc))
      Ftot = 0.0
      do j=1,size(cos_theta)

         rr = Rmax * sqrt(1.0 - cos_theta(j)**2)

         x0 = 10.0*Rmax*u + rr * y(1)
         y0 = 10.0*Rmax*v + rr * y(2)
         z0 = 10.0*Rmax*w + rr * y(3)

         call move_to_grid(1,x0,y0,z0,u0,v0,w0,icell0,lintersect)
         if (lintersect) then
            call integ_ray_atom(1,icell0,x0,y0,z0,u0,v0,w0,1,.false.,n_lambda,km_to_m*tab_lambda)
         endif
         Ftot(:,1) = Ftot(:,1) + Itot(:,1,1) * weight_mu(j)
      enddo
      open(1,file="sun.txt",status="unknown")
      do la=1,n_lambda

         write(1,*) tab_lambda(la), Ftot(la,1)
      enddo
      close(1)



      return
   end subroutine atom_line_transfer

   ! subroutine intensite_pixel_atom(id,imol,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
   !    ! Calcule l'intensite d'un pixel carre de taille, position et orientation arbitaires
   !    ! par une methode de Ray-tracing
   !    ! (u,v,w) pointe vers l'observateur
   !    ! Integration par methode de Romberg pour determiner le nbre de sous-pixel
   !    ! necessaire
   !    ! Unite : W.m-2 : nu.F_nu
   !    ! C. Pinte
   !    ! 12/04/07
    
   !    implicit none
    
   !    integer, intent(in) :: ipix,jpix,id, imol, n_iter_min, n_iter_max, ibin, iaz
   !    real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
   !    real(kind=dp), intent(in) :: pixelsize,u,v,w
   !    real(kind=dp), dimension(:,:), allocatable :: IP, IP_old
   !    real(kind=dp), dimension(:), allocatable :: IPc
    
   !    integer, parameter :: maxSubPixels = 32
    
   !    real(kind=dp) :: x0,y0,z0,u0,v0,w0
   !    real(kind=dp), dimension(3) :: sdx, sdy
   !    real :: npix2, diff
    
   !    real, parameter :: precision = 1.e-2
   !    integer :: i, j, subpixels, iray, ri, zj, phik, icell, iTrans, iter, n_speed_rt, nTrans_raytracing
    
   !    logical :: lintersect, labs
    
   !    integer, dimension(2) :: ispeed
    
   !    n_speed_rt = mol(imol)%n_speed_rt
   !    nTrans_raytracing = mol(imol)%nTrans_raytracing
    
   !    allocate(IP(-n_speed_rt:n_speed_rt,nTrans_raytracing), IP_old(-n_speed_rt:n_speed_rt,nTrans_raytracing), IPc(nTrans_raytracing))
    
   !    ispeed(1) = -n_speed_rt ; ispeed(2) = n_speed_rt
    
    
   !    labs = .false.
    
   !    ! Ray tracing : on se propage dans l'autre sens
   !    u0 = -u ; v0 = -v ; w0 = -w
    
   !    IP = 0.0_dp
   !    IPc = 0.0_dp
    
   !    ! le nbre de subpixel en x est 2^(iter-1)
   !    subpixels = 1
   !    iter = 1
    
   !    infinie : do ! Boucle infinie tant que le pixel n'est pas converge
   !       npix2 =  real(subpixels)**2
   !       IP_old = IP
   !       IP = 0.0_dp
   !       IPc = 0.0_dp
    
   !       ! Vecteurs definissant les sous-pixels
   !       sdx(:) = dx(:) / real(subpixels,kind=dp)
   !       sdy(:) = dy(:) / real(subpixels,kind=dp)
    
   !       iray = 1
    
   !       ! L'obs est en dehors de la grille
   !       ri = 2*n_rad ; zj=1 ; phik=1
    
   !       ! Boucle sur les sous-pixels qui calcule l'intensite au centre
   !       ! de chaque sous pixel
   !       do i = 1,subpixels
   !          do j = 1,subpixels
   !             ! Centre du sous-pixel
   !             x0 = pixelcorner(1) + (i - 0.5_dp) * sdx(1) + (j-0.5_dp) * sdy(1)
   !             y0 = pixelcorner(2) + (i - 0.5_dp) * sdx(2) + (j-0.5_dp) * sdy(2)
   !             z0 = pixelcorner(3) + (i - 0.5_dp) * sdx(3) + (j-0.5_dp) * sdy(3)
    
   !             ! On se met au bord de la grille : propagation a l'envers
   !             call move_to_grid(id, x0,y0,z0,u0,v0,w0, icell,lintersect)
    
   !             if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
   !                call integ_ray_mol(id,imol,icell,x0,y0,z0,u0,v0,w0,iray,labs,ispeed,tab_speed_rt, &
   !                     nTrans_raytracing, mol(imol)%indice_Trans_raytracing)
   !                ! Flux recu dans le pixel
   !                IP(:,:) = IP(:,:) +  I0(:,:,iray,id)
   !                IPc(:) = IPc(:) +  I0c(:,iray,id)
   !             else
   !                ! Il n'y a que le Cmb
   !                ! TODO : je ne suis pas sur de vouloir le Cmb en dehors du disque, ...
   !                ! IP(:,:) = IP(:,:) +  Cmb(ispeed,tab_speed)
   !             endif
   !          enddo !j
   !       enddo !i
    
   !       IP = IP / npix2
   !       IPc = IPc / npix2
    
   !       if (iter < n_iter_min) then
   !          ! On itere par defaut
   !          subpixels = subpixels * 2
   !       else if (iter >= n_iter_max) then
   !          ! On arrete pour pas tourner dans le vide
   !          ! write(*,*) "Warning : converging pb in ray-tracing"
   !          ! write(*,*) " Pixel", ipix, jpix
   !          exit infinie
   !       else
   !          ! On fait le test sur a difference
   !          diff = maxval( abs(IP - IP_old) / (IP + 1e-300_dp) )
   !          if (diff > precision ) then
   !             ! On est pas converge
   !             subpixels = subpixels * 2
   !          else
   !             ! On est converge
   !             exit infinie
   !          endif
   !       endif ! iter
    
   !       iter = iter + 1
    
   !       ! TODO : Integration Romberg
   !  !!$     if(any(abs((I - oldintensite_pixel)) > precision * I)) then
   !  !!$        oldintensite_pixel = I
   !  !!$        ! Il n'y a pas un truc mieux pour utiliser les calculs a plus faible resol ??
   !  !!$
   !  !!$        subpixels = subpixels * 2
   !  !!$        !if(subpixels .gt. 15) write(*,*)"large",index
   !  !!$     else
   !  !!$        I = (real(2**(log(real(subpixels))/log(2.d0)))*I - Oldintensite_pixel) &
   !  !!$             /(real(2**(log(real(subpixels))/log(2.d0)))-1.d0) ! Richardson Extrapolation
   !  !!$        ! Ok mais n'utilise que les 2 derniers calculs : il doit y avoir mieux !!!
   !  !!$
   !  !!$     endif
   !    enddo infinie
    
   !    ! Prise en compte de la surface du pixel (en sr)
   !    IP = IP * (pixelsize / (distance*pc_to_AU) )**2
   !    IPc = IPc * (pixelsize / (distance*pc_to_AU) )**2
    
   !    ! et multiplication par la frequence pour avoir du nu.F_nu
    
   !    ! Warning IP, IPc are smaller array (dimension mol(imol)%nTrans_raytracin)
   !    do i=1,mol(imol)%nTrans_raytracing
   !       iTrans = mol(imol)%indice_Trans_raytracing(i)
   !       IP(:,i) = IP(:,i) * transfreq(iTrans)
   !       IPc(i) = IPc(i) * transfreq(iTrans)
   !    enddo
   !    ! Unite teste OK pour le Cmb
   !    ! profil de raie non convolue teste ok avec torus
    
   !    if (RT_line_method==1) then ! Sommation implicite sur les pixels
   !       spectre(1,1,:,:,ibin,iaz) = spectre(1,1,:,:,ibin,iaz) + IP(:,:)
   !       continu(1,1,:,ibin,iaz) = continu(1,1,:,ibin,iaz) + IPc(:)
   !    else
   !       spectre(ipix,jpix,:,:,ibin,iaz) = IP(:,:)
   !       continu(ipix,jpix,:,ibin,iaz) = IPc(:)
   !    endif
    
   !    return
    
   ! end subroutine intensite_pixel_atom
    

end module atom_transfer