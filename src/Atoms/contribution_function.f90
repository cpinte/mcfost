 ! 
!  !building
!  SUBROUTINE compute_contribution_function()
!   ! ---------------------------------------------------- !   
!   ! Following Albrow & Cottrell 1995
!   ! Ksi(logtau0) = log(10)*tau0 / chi0 *
!   !   int_0^1 (chil * (Ic - Sl) * exp(-tau) ), dlogtau0
!   !
!   ! here expressed has the "line flux depression" CF.
!   ! See also Amarsi 2015
!   ! ---------------------------------------------------- ! 
! #include "sprng_f.h"
!    integer :: icell, az, ibin, iray, n_rayons, irayprop, id
!    integer :: i
!    integer :: ipix, jpix
!    logical :: labs
!    real(kind=dp) :: tauR, chiR, w0, v0, u0, x0, y0, z0, argmt, srw02, w02
!    real ::  rand, rand2, rand3, stream(nb_proc)
!    real(kind=dp), dimension(NLTEspec%NWaves) :: chil
!    real(kind=dp), dimension(NLTEspec%Nwaves,nb_proc) :: S_contrib
!    type (AtomType), pointer :: atom
!    
!    labs = .false.
!    n_rayons = 10
!    irayprop = 1
!    id = 1
!    !note that the propagation is done for one ray, we di directly the sum. Otherwise, the
!    !sizes of I and Ic have to be redifined, because we enter here after image ATM.
! 
!    stream = 0.0
!    do i=1,nb_proc
!      stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
!    enddo
!    
!    !$omp parallel &
!    !$omp default(none) &
!    !$omp private(id,iray,ipix, jpix, icell,atom,rand,rand2,rand3) &
!    !$omp private(w0, u0, v0, w02, srW02, argmt,x0,y0,z0, chil) &
!    !$omp shared(atmos, n_cells, NLTEspec, labs, n_rayons, stream, S_contrib, irayprop)
!    !$omp do schedule(static,1)   
!    do icell=1, n_cells
!    	!$ id = omp_get_thread_num() + 1
!    	if (atmos%icompute_atomRT(icell)>0) then !nor transparent nor dark
!      do iray=1, n_rayons
!        rand  = sprng(stream(id))
!        rand2 = sprng(stream(id))
!        rand3 = sprng(stream(id))
!                         
!        CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)
!                         
!        ! Direction de propagation aleatoire
!        rand = sprng(stream(id))
!        W0 = 2.0_dp * rand - 1.0_dp
!        W02 =  1.0_dp - W0*W0
!        SRW02 = sqrt(W02)
!        rand = sprng(stream(id))
!        ARGMT = PI * (2.0_dp * rand - 1.0_dp)
!        U0 = SRW02 * cos(ARGMT)
!        V0 = SRW02 * sin(ARGMT)  
! 
!       !in iray direction
!       S_contrib(:,id) = 0d0
!       
!       !iray is here to loop, but I do not store Ic(iray) and I(iray) nor Ksi(iray)
!       CALL propagate_ray_contrib(id,icell,x0,y0,z0,u0,v0,w0,irayprop,labs, S_contrib)
!       !for the moment not NLTE because the lack of NLTE cont separated from lines
!       !					line+cont			-          cont   = line
!       chil(:) = (NLTEspec%AtomOpac%chi_p(:,id) - NLTEspec%AtomOpac%chi_c(:,id))/NLTEspec%AtomOpac%chi_p(:,id)
!       
!       NLTEspec%ksi(icell,:) = NLTEspec%ksi(icell,:) + chil(:)*(NLTEspec%Ic(:,irayprop,id)-S_Contrib(:,id))/n_rayons
!       
!      enddo !iray
!     endif !not dark and nor transparent
!    enddo !icell
!    !$omp end do
!    !$omp end parallel
!    
!    CALL WRITE_CNTRB()
!    
!  RETURN
!  END SUBROUTINE compute_contribution_function
!  
!  SUBROUTINE propagate_ray_contrib(id,icell_in,x,y,z,u,v,w,iray,labs, S_contrib)
!  ! ------------------------------------------------------------------------------- !
!   !
!  ! ------------------------------------------------------------------------------- !
! 
!   integer, intent(in) :: id, icell_in, iray
!   real(kind=dp), intent(in) :: u,v,w
!   real(kind=dp), intent(in) :: x,y,z
!   real(kind=dp), intent(out) :: S_contrib(:,:)
!   logical, intent(in) :: labs !used in NLTE but why?
!   real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, etau, Istar
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau
!   integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
!   logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator
! 
!   x1=x;y1=y;z1=z
!   x0=x;y0=y;z0=z
!   next_cell = icell_in
!   nbr_cell = 0
! 
!   tau = 0.0_dp !go from surface down to the star
!   tau_c = 0.0_dp
! 
!   NLTEspec%I(:,iray,id) = 0d0
!   NLTEspec%Ic(:,iray,id) = 0d0
!   Istar(:) = 0d0
!   
!   ! -------------------------------------------------------------- !
!   !*** propagation dans la grille ***!
!   ! -------------------------------------------------------------- !
!   ! Will the ray intersect a star
!   call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)
! 
!   ! Boucle infinie sur les cellules
!   infinie : do ! Boucle infinie
!     ! Indice de la cellule
!     icell = next_cell
!     x0=x1 ; y0=y1 ; z0=z1
! 
!     if (icell <= n_cells) then
!      lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
!      if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
!     else
!      lcellule_non_vide=.false.
!     endif
!     
!     ! Test sortie ! "The ray has reach the end of the grid"
!     if (test_exit_grid(icell, x0, y0, z0)) RETURN
!     
! 
!     if (lintersect_stars) then
!       if (icell == icell_star) then
!        CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,Istar)
!        NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + Istar*dexp(-tau)
!        NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + Istar*dexp(-tau_c)
!        RETURN
!       end if
!     endif
! 
!     nbr_cell = nbr_cell + 1
! 
!     ! Calcul longeur de vol et profondeur optique dans la cellule
!     previous_cell = 0 ! unused, just for Voronoi
! 
!     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
! 
!     if (lcellule_non_vide) then
! 
!      ! opacities in m^-1
!      l_contrib = l_contrib * AU_to_m !l_contrib in AU
!      
!      !if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m
!      eval_operator = .false. !labs if false for cntrb function
!      
!      CALL initAtomOpac(id,eval_operator) !set opac to 0 for this cell and thread id
!      CALL NLTEopacity(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
!      !never enter NLTEopacity if no activeatoms
!      if (lstore_opac) then !not updated during NLTE loop, just recomputed using initial pops
!       CALL BackgroundLines(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
!       dtau(:)   = l_contrib * (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id)+&
!                     NLTEspec%AtomOpac%Kc(icell,:,1))
!       dtau_c(:) = l_contrib * NLTEspec%AtomOpac%Kc(icell,:,1)
!       Snu = (NLTEspec%AtomOpac%jc(icell,:) + &
!                NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / &
!                  (NLTEspec%AtomOpac%Kc(icell,:,1) + &
!                  NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%chi(:,id) + 1d-300)
!                  ! + &NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%J(:,id)
!                  ! + &NLTEspec%AtomOpac%Kc(icell,:,2) * NLTEspec%Jc(:,id)
!       Snu_c = (NLTEspec%AtomOpac%jc(icell,:)) / &
!             (NLTEspec%AtomOpac%Kc(icell,:,1) + 1d-300) !it is missing NLTE cont
!      else
!       CALL Background(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l) !x,y,z,u,v,w,x1,y1,z1
!                                 !define the projection of the vector field (velocity, B...)
!                                 !at each spatial location.
! 
!       dtau(:)   =  l_contrib * (NLTEspec%AtomOpac%chi_p(:,id)+NLTEspec%AtomOpac%chi(:,id))
!       dtau_c(:) = l_contrib * NLTEspec%AtomOpac%chi_c(:,id)
! 
!       Snu = (NLTEspec%AtomOpac%eta_p(:,id) + NLTEspec%AtomOpac%eta(:,id)) / &
!                  (NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%chi(:,id) + 1d-300)
!       ! + &NLTEspec%AtomOpac%sca_c(:,id) * NLTEspec%J(:,id)
!       ! + &NLTEspec%AtomOpac%sca_c(:,id) * NLTEspec%Jc(:,id)
!       
!       ! continuum source function
!       Snu_c = (NLTEspec%AtomOpac%eta_c(:,id)) / (NLTEspec%AtomOpac%chi_c(:,id) + 1d-300)
!     end if
!     if (minval(Snu) < 0. .or. maxval(Snu) < 0.) then
!       write(*,*) "eta (max/min)", maxval(NLTEspec%AtomOpac%eta(:,id)), minval(NLTEspec%AtomOpac%eta(:,id))
!       write(*,*) "chi (max/min)", maxval(NLTEspec%AtomOpac%chi(:,id)), minval(NLTEspec%AtomOpac%chi(:,id))
!       call Warning("Snu negative")   
!     end if
! 
!     NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + &
!     						 dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
!     NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + &
!                              dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c
!           
! 
!      tau = tau + dtau
!      tau_c = tau_c + dtau_c
!      
!      !at icell
!      S_contrib(:,id) = Snu(:) * dexp(-dtau) !*exp(-tau)
!      if (nbr_cell==1) RETURN !get the intensity at this cell
! 
!     end if  ! lcellule_non_vide
!   end do infinie
!   ! -------------------------------------------------------------- !
!   ! -------------------------------------------------------------- !
!   RETURN
!   END SUBROUTINE propagate_ray_contrib

 
!  SUBROUTINE WRITE_CNTRB_FUNC()
!  ! -------------------------------------------------- !
!   ! Write contribution function to disk.
!  ! --------------------------------------------------- !
!   integer :: status,unit,blocksize,bitpix,naxis
!   integer, dimension(7) :: naxes
!   integer :: group,fpixel,nelements
!   logical :: simple, extend
!   character(len=6) :: comment=""
!   real :: pixel_scale_x, pixel_scale_y 
!   
!    !  Get an unused Logical Unit Number to use to open the FITS file.
!    status=0
!    CALL ftgiou (unit,status)
! 
!    !  Create the new empty FITS file.
!    blocksize=1
!    CALL ftinit(unit,"cntrb.fits.gz",blocksize,status)
! 
!    simple=.true.
!    extend=.true.
!    group=1
!    fpixel=1
! 
!    bitpix=-64
!    
!   if (lVoronoi) then
!    naxis = 2
!    naxes(1) = NLTEspec%atmos%Nspace ! equivalent n_cells
!    naxes(2) = NLTEspec%Nwaves
!    nelements = naxes(1) * naxes(2)
!   else
!    if (l3D) then
!     naxis = 4
!     naxes(1) = n_rad
!     naxes(2) = 2*nz
!     naxes(3) = n_az
!     naxes(4) = NLTEspec%Nwaves
!     nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
!    else
!     naxis = 3
!     naxes(1) = n_rad
!     naxes(2) = nz
!     naxes(3) = NLTEspec%Nwaves
!     nelements = naxes(1) * naxes(2) * naxes(3) 
!    end if
!   end if
! 
!   !  Write the required header keywords.
!   CALL ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
!   CALL ftpkys(unit,'UNIT',"W.m-2.Hz-1.pixel-1",'Ksi',status)
! 
! 
!   !  Write the array to the FITS file.
!   CALL ftpprd(unit,group,fpixel,nelements,NLTEspec%ksi,status)
! 
!   !  Close the file and free the unit number.
!   CALL ftclos(unit, status)
!   CALL ftfiou(unit, status)
! 
!   !  Check for any error, and if so print out error messages
!   if (status > 0) then
!      CALL print_error(status)
!   endif
! 
!  RETURN
!  END SUBROUTINE WRITE_CNTRB_FUNC  