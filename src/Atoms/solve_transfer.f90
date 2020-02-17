

! A couple of forsaken subroutines
! Just here as a reminder

! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atoms.
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

! SUBROUTINE of AtomicTransfer line integ_ray_line will be move here
! generator for directions also ect, in atomictransfer, only the NLTEloop and the flux
!will remain

MODULE solve_transfer

 use metal, only                        : metal_bb_new, compute_opacities, alloc_atom_quantities, dealloc_Atom_quantities
 use opacity
 use spectrum_type
 use atmos_type
 use constant, only 					: MICRON_TO_NM
 use impact
 use solvene
 use math
 !$ use omp_lib
 use parametres
 use grid
 use stars
 use mcfost_env, only : dp
 use constantes, only : tiny_dp, huge_dp
 use init_solution, only : lcell_converged
 use accelerate

 IMPLICIT NONE
 
 private
 public :: Iterate_Jnu

 
 real(kind=dp), parameter :: tiny_chi = 1d-50

 CONTAINS
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
!  SUBROUTINE integrate_tau_bound(id,icell_in,x,y,z,u,v,w,iray,to_obs, tau)
!  ! ------------------------------------------------------------------------------- !
!  ! ------------------------------------------------------------------------------- !
! 
!   integer, intent(in) :: id, icell_in, iray, to_obs
!   real(kind=dp), intent(in) :: u,v,w
!   real(kind=dp), intent(in) :: x,y,z
!   real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
!   real(kind=dp), dimension(NLTEspec%Nwaves) :: dtau, chiI
!   real(kind=dp), intent(out), dimension(NLTEspec%NWaves,atmos%Nrays, NLTEspec%Nproc) :: tau
!   integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
!   logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars
! 
!   x1=x;y1=y;z1=z
!   x0=x;y0=y;z0=z
!   next_cell = icell_in
!   nbr_cell = 0
! 
!   tau(:,iray,id) = 0d0
!   !!write(*,*) " vector direction times ", to_obs, u,v,w
! 
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
!       if (icell == icell_star) RETURN
!     endif
! 
!     nbr_cell = nbr_cell + 1
! 
!     previous_cell = 0
! 
!     call cross_cell(x0,y0,z0, u,v,w, icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
! !     if (l_contrib /= l) then
! !       write(*,*) "l_contrib, l"
! ! 	  write(*,*) l_contrib/ etoile(1)%r, l / etoile(1)%r
! ! 	endif
!     if (lcellule_non_vide) then
! 
! 
!      l_contrib = l_contrib * AU_to_m !l_contrib in AU
! 
!      
!      CALL initAtomOpac(id) !set opac to 0 for this cell and thread id
!      
!       if (atmos%NpassiveAtoms>0) CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)
!       
!       chiI(:) = NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%Kc(:,icell,1) + tiny_chi
! 
!       !add only if non zero,; this should be fater to do than adding a large 0 array
!       if (atmos%NactiveAtoms>0) then 
!        CALL initAtomOpac_nlte(id)
!        CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)
!        chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id) + NLTEspec%AtomOpac%Kc_nlte(:,icell)
!       endif
!                     
!       dtau(:)   = l_contrib * chiI(:)
! 
! 
! 	 !do not take into account dtau of the cell if we integrate tau in the same direction of propagation
! 	 !but in opposite direction yes.
! 	 if (to_obs == 1) then !going in the same direction as the ray
! 	   !do not take into account the path in the cell for this case as we want tau at the edge
!        if (nbr_cell > 1) tau(:,iray,id) = tau(:,iray,id) + dtau
!      else if (to_obs==-1) then !opposite direction, count the travel path in the cell
!        tau(:,iray,id) = tau(:,iray,id) + dtau
!      else
!       write(*,*) to_obs
!       CALL ERROR("unknown value for to_obs")
!      endif
! 
!     end if
!   end do infinie
!   ! -------------------------------------------------------------- !
!   ! -------------------------------------------------------------- !
!   RETURN
!   END SUBROUTINE integrate_tau_bound

! 						CALL INTEG_RAY_JNU_ALI(id, icell, xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
! 											uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, labs, &
! 											kappa_tot, Sold, Istar, lambda_star(:,id), Ic(:,id))
! 											
! 						J_FS(:,icell) = J_FS(:,icell) + Ic(:,id) / n_rayons						
! 						lambda_star_j(:,id) = lambda_star_j(:,id) + lambda_star(:,id) / n_rayons
						
!       			    S_FS(:,id) = Sth(:,icell) * (1.-beta(:,icell)) + beta(:,icell) * J_FS(:,icell)
!       			    Snew(:,icell) = (S_FS(:,id) - beta(:,icell) * lambda_star_j(:,id) * Sold(:,icell)) / (1.-beta(:,icell)*lambda_star_j(:,id))
!       			    Jnew(:,icell) = J_FS(:,id)
!    SUBROUTINE INTEG_RAY_JNU_ALI(id,icell_in,x,y,z,u,v,w,iray,labs, kappa_tot, Sny, Istar, psi, Ic)
!  ! ------------------------------------------------------------------------------- !
!   ! This routine performs integration of the transfer equation along a ray to
!   ! compute coherent Jnu in the continuum
!  ! ------------------------------------------------------------------------------- !
! 
!   integer, intent(in) :: id, icell_in, iray
!   real(kind=dp), intent(in) :: u,v,w
!   real(kind=dp), intent(in) :: x,y,z
!   real(kind=dp), intent(in), dimension(:,:) :: kappa_tot, Sny
!   real(kind=dp), intent(in), dimension(:) :: Istar
!   real(kind=dp), intent(out) :: Ic(:), psi(:)
!   logical, intent(in) :: labs
!   real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
!   real(kind=dp), dimension(size(Ic)) :: LimbD!, Snu_c
!   real(kind=dp), dimension(size(Ic)) :: tau_c!, dtau_c, chiIc
!   integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
!   logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator
! 
!   x1=x;y1=y;z1=z
!   x0=x;y0=y;z0=z
!   next_cell = icell_in
!   nbr_cell = 0
! 
!   !tau(:) = 0.0_dp !go from surface down to the star
!   tau_c(:) = 0.0_dp
! 
!   Ic(:) = 0.0
!   psi(:) = 0.0
!   
!   eval_operator = .false.
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
!      !lcellule_non_vide=.true.
!      !atmos%nHtot(icell)>0
!      lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
!      if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
!     else
!      lcellule_non_vide=.false.
!     endif
!     
!     if (minval(tau_c) > 20.) return
!     
!     ! Test sortie ! "The ray has reach the end of the grid"
!     if (test_exit_grid(icell, x0, y0, z0)) RETURN
! 
!     if (lintersect_stars) then
!       if (icell == icell_star) then
!        !CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
!        Ic(:) =  Ic(:) + Istar(:) * dexp(-tau_c)
!        RETURN
!       end if
!     endif
! 
!     nbr_cell = nbr_cell + 1
! 
!     ! Calcul longeur de vol et profondeur optique dans la cellule
!     previous_cell = 0 ! unused, just for Voronoi
!     call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
! 
! 
!     !count opacity only if the cell is filled, else go to next cell
!     if (lcellule_non_vide) then
!      lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
!      ! opacities in m^-1
!      l_contrib = l_contrib * AU_to_m !l_contrib in m    
! 
!       Ic(:) = Ic(:) + dexp(-tau_c) * (1.0_dp - dexp(-l_contrib * kappa_tot(:,icell))) * Sny(:,icell)
! 
!      
!      if ((nbr_cell == 1).and.labs) then 
!       ds(iray,id) = l * AU_to_m
!       psi(:) = (1d0 - exp(-ds(iray,id)*kappa_tot(:,icell))) * dexp(-tau_c)!in this direction
!      endif
! 
! 
!      !tau = tau + dtau
!      tau_c = tau_c + l_contrib * kappa_tot(:,icell)
! 
!     end if  ! lcellule_non_vide
!   end do infinie
!   ! -------------------------------------------------------------- !
!   ! -------------------------------------------------------------- !
!   RETURN
!   END SUBROUTINE INTEG_RAY_JNU_ALI
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- Jnu -------------------------------------------- !
 ! ------------------------------------------------------------------------------------ !
 ! done here because Nrays should be > 1 and it it sets to 1 for images
!   if (lelectron_scattering .and. atmos%NactiveAtoms==0) then 
!    write(*,*) "  -- Computing Jnu"
!    CALL alloc_Jnu() !temporary, in NLTE allocated in initSol(), not in alloc_spectrum anymore
!    allocate(lcell_converged(atmos%Nspace))
!    lcell_converged(:) = .false.
!    !init to Planck function
!    do icell=1,n_cells
!     if (atmos%icompute_atomRT(icell)>0) then
!      CALL Bplanck(atmos%T(icell), NLTEspec%J(:,icell))
!      CALL Bplanck(atmos%T(icell), NLTEspec%Jc(:,icell))
!     endif
!    enddo
!    CALL Iterate_Jnu()
!   endif

  SUBROUTINE INTEG_RAY_LINE_Jnu(id,icell_in,x,y,z,u,v,w,iray,labs, chi1, chi2, eta1, eta2)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(out), dimension(NLTEspec%Nwaves) :: chi1, chi2, eta1, eta2
  logical, intent(in) :: labs !used in NLTE but why?
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, etau, edtau, etauc, edtauc
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, chiI, LimbD, chiIc
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, dtau_ds
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau(:) = 0.0_dp !go from surface down to the star
  tau_c(:) = 0.0_dp
  chiI(:) = 0.0_dp !not needed ?
  chiIc(:) = 0.0_dp

  chi1(:) = 0.0_dp
  chi2(:) = 0.0_dp
  
  eta1(:) = 0.0_dp
  eta2(:) = 0.0_dp

  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0
  
  eval_operator = .false.
  ! -------------------------------------------------------------- !
  !*** propagation dans la grille ***!
  ! -------------------------------------------------------------- !
  ! Will the ray intersect a star
  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  ! Boucle infinie sur les cellules
  infinie : do ! Boucle infinie
    ! Indice de la cellule
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     !lcellule_non_vide=.true.
     !atmos%nHtot(icell)>0
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif
    
    
    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN

    if (lintersect_stars) then
      if (icell == icell_star) then
       CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau_c)
       RETURN
      end if
    endif

    nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
    previous_cell = 0 ! unused, just for Voronoi
    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)


    !count opacity only if the cell is filled, else go to next cell
    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m !l_contrib in m

     CALL initAtomOpac(id)
     
      if (atmos%NpassiveAtoms>0) CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)

      chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell,1) + tiny_chi 
      chiI(:)  = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)
                    
      dtau(:)   = l_contrib * chiI(:)
      dtau_c(:) = l_contrib * chiIc(:)
            
      Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
      Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
      
      Snu = Snu + NLTEspec%J(:,icell) * NLTEspec%AtomOpac%Kc(:,icell,2) / chiI(:)
      Snu_c = Snu_c + NLTEspec%AtomOpac%Kc(:,icell,2) * NLTEspec%Jc(:,icell) / chiIc(:)

      NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
      NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c
     
     if ((nbr_cell == 1).and.labs) then 
      ds(iray,id) = l * AU_to_m
      chi1(:) = chiI
      chi2(:) = chiIc
      eta1(:) = Snu
      eta2(:) = Snu_c
     endif


     tau = tau + dtau
     tau_c = tau_c + dtau_c

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_Jnu


 SUBROUTINE Iterate_Jnu()
 ! -------------------------------------------------------- !
  ! Compute the mean radiation field at all cells
 ! -------------------------------------------------------- !
#include "sprng_f.h"

  integer, parameter :: n_rayons_start = 100
  integer :: n_rayons_max
  real, parameter :: precision = 1e-2
  integer :: etape, etape_start, etape_end, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lprevious_converged

  real :: rand, rand2, rand3!, fac_etape

  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dN, dJ, dJc, dNc, diffc
  real(kind=dp), allocatable :: Jold(:,:), Jcold(:,:), Jflat(:), chi(:,:,:), chic(:,:,:), Snu(:,:,:), Snu_c(:,:,:)
                                   !futur global flag
  logical :: labs, l_unconverged, accelerated, ng_rest, ng_acc
  integer :: maxIter
  integer :: la, icell, imax, icell_max, iorder,i0_rest
  real(kind=dp) :: lambda_max
  real(kind=dp), dimension(3, atmos%Nrays, nb_proc) :: xyz0, uvw0
  type(Ng) :: NgJ, NgJc
  
  if (allocated(ds)) deallocate(ds)
  allocate(ds(atmos%Nrays, NLTEspec%NPROC))
  ds = 0.0_dp !meters

  !move to initSol
  allocate(Jold(NLTEspec%Nwaves, atmos%Nspace))!move to initsol
  Jold = 0d0
  !only in "LTE" case
  allocate(Jcold(NLTEspec%Nwaves, atmos%Nspace))
  Jcold = 0d0
  
  allocate(chi(NLTEspec%Nwaves, atmos%nrays, NLTEspec%Nproc), chic(NLTEspec%Nwaves, atmos%nrays, NLTEspec%Nproc))
  allocate(Snu(NLTEspec%Nwaves, atmos%nrays, NLTEspec%Nproc), Snu_c(NLTEspec%Nwaves, atmos%nrays, NLTEspec%Nproc))

  n_rayons_max = atmos%Nrays
  xyz0(:,:,:) = 0d0
  uvw0(:,:,:) = 0d0
  
  ng_acc = .false. !temp, because don't know if it works
  if (ng_acc) then
   allocate(Jflat(atmos%Nspace*NLTEspec%Nwaves))
   CALL initNg(atmos%Nspace*NLTEspec%Nwaves, 6, 2, 3, NgJ)
  endif


  labs = .true.
  id = 1
  etape_start = 2
  etape_end = 2
  if (etape_start==0) then
   write(*,*) "Warn: etape 0 not accurate"
  else if (etape_start==1 .or. etape_end==1) then
   write(*,*) " Etape 1 with new angle quad not implemented, setting to 2"
   if (etape_start==1) etape_start=2
   if (etape_end==1) etape_end=2
  endif
  
  maxIter = 20
 
     do etape=etape_start, etape_end

      if (etape==0) then !two rays, not accurate
  	    CALL ERROR("etape 0 not allowed for Jnu")

      else if (etape==1) then 
  	    CALL ERROR("etape 1 not allowed for Jnu")

      else if (etape==2) then 
  		lfixed_rays = .true.
  		n_rayons = min(n_rayons_max,n_rayons_start)
  		iray_start = 1
  		lprevious_converged = .false.
		write(*,*) " Using ", n_rayons, " rays for Jnu."
  	  else
  	    CALL ERROR("etape unkown")
  	  end if
  	  
  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0

        do while (.not.lconverged)

        	n_iter = n_iter + 1
        	if (n_iter > maxIter) exit
            write(*,*) " -> Iteration #", n_iter, " Step #", etape

  			if (lfixed_rays) then
    			stream = 0.0
    			do i=1,nb_proc
     				stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    			end do
 			 end if

            diff = 0.0_dp
            
            Jold(:,:) = NLTEspec%J(:,:)
            Jcold(:,:) = NLTEspec%Jc(:,:)

 			!$omp parallel &
            !$omp default(none) &
            !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dN, dJ, dNc, dJc, diffc)&
            !$omp private(argmt,diff,norme, icell, l_unconverged) &
            !$omp shared(atmos,NLTEspec) &
            !$omp shared(xyz0, uvw0, lkeplerian,n_iter) &
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,lcell_converged) &
            !$omp shared(n_cells,chi,chic,snu,snu_c,ds) &
            !$omp shared(lfixed_Rays,lnotfixed_Rays,labs,etape)
            !$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !$ id = omp_get_thread_num() + 1
   				l_unconverged = (atmos%icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
   				if (l_unconverged) then

                    CALL compute_directions(etape, id, icell, iray_start, n_rayons, stream(id), xyz0(:,:,id), uvw0(:,:,id))
                    
					do iray=iray_start, iray_start-1+n_rayons

						CALL INTEG_RAY_LINE_Jnu(id, icell, xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
											uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, labs, &
											chi(:,iray,id), chic(:,iray,id), Snu(:,iray,id), Snu_c(:,iray,id))
		
      			    enddo !iray
      			    !apply ALI Hogerheijde here	
!       			    NLTEspec%J(:,icell) = 0.
!       			    NLTEspec%Jc(:,icell) = 0.
!       			    do iray=1, n_rayons
!       			     NLTEspec%J(:,icell) = NLTEspec%J(:,icell) + ( NLTEspec%I(:,iray,id)*dexp(-ds(iray,id)*chi(:,iray,id)) + &
!       			         (1.-dexp(-ds(iray,id)*chi(:,iray,id)))*Snu(:,iray,id) )/n_rayons
!       			         
!       			     NLTEspec%Jc(:,icell) = NLTEspec%Jc(:,icell) + ( NLTEspec%Ic(:,iray,id)*dexp(-ds(iray,id)*chic(:,iray,id)) + &
!       			         (1.-dexp(-ds(iray,id)*chic(:,iray,id)))*Snu_c(:,iray,id) )/n_rayons
!       			    enddo	
					CALL calc_Jnu(id, icell, n_rayons)
      		   endif
     		end do !icell
        	!$omp end do
        	!$omp end parallel
     		accelerated=.false.
     	if (ng_acc) then
     		if (n_iter > 6) then
     		  iorder = n_iter - 6 !local number of iterations accumulated
     		  if (ng_rest) then
     		    write(*,*) "    --> Acceleration relaxes...", iorder-i0_rest
     		    if (iorder - i0_rest == 3) ng_rest = .false.
     		  else
     		     i0_rest = iorder
     		     Jflat=flatten2(NLTEspec%Nwaves,atmos%Nspace,NLTEspec%J)
     		     accelerated = Acceleration(NgJ, Jflat)
     		     if (accelerated) then
     		      NLTEspec%J = reform2(NLTEspec%Nwaves,atmos%Nspace,Jflat)
                  ng_rest = .true.
     		     endif
             endif
     		endif
        endif
            !should be para
        	diff = 0d0
        	diffc = 0.
        	icell_max = 0
  			cell_loop2 : do icell=1, atmos%Nspace
  				if (atmos%icompute_atomRT(icell)>0) then

						dJ = 0.0_dp
						dJc = 0.0_dp
						do la=1, NLTEspec%Nwaves
						 if (NLTEspec%J(la, icell) > 0) then 
						  dJ = max(dJ,dabs(1.-Jold(la,icell)/NLTEspec%J(la,icell)))
						  imax = locate(dabs(1.-Jold(:,icell)/NLTEspec%J(:,icell)),dabs(1.-Jold(la,icell)/NLTEspec%J(la,icell)))
						 endif 
						 if (NLTEspec%Jc(la, icell) > 0) dJc = max(dJc,dabs(1.-Jcold(la,icell)/NLTEspec%Jc(la,icell))) 
						enddo
						if (dJ > diff) then
						  diff = dJ
						  icell_max = icell
						endif
     					diffc = max(diffc,dJc)
     
     					lcell_converged(icell) = (real(dJ) < precision)		
     			end if
     		end do cell_loop2 !icell
     		
         	if (accelerated) then
         	  write(*,*) "   >>> ", icell_max, NLTEspec%lambda(imax)," dJ = ", diff, " dJc=", diffc," (Accelerated)"
         	else
         	  write(*,*) "   >>> ", icell_max, NLTEspec%lambda(imax)," dJ = ", diff, " dJc=", diffc
         	endif
        	lconverged = (real(diff) < precision)

	    end do !while
        write(*,*) etape, "Threshold =", dpops_max_error
	  end do !over etapes

  CALL freeNg(NgJ)
  deallocate(Jold, Jcold, chi, chic, Snu, Snu_c)

 ! ------------------------------------------------------------------------------------ !
 RETURN
 END SUBROUTINE Iterate_Jnu
 ! ------------------------------------------------------------------------------------ !

 SUBROUTINE compute_directions(etape, id, icell, rstart, n_rayons, streami, xyz, uvw)
#include "sprng_f.h"

 	integer, intent(in) :: etape, id, icell, rstart, n_rayons
 	SPRNG_POINTER, intent(in) :: streami
 	real(kind=dp), dimension(3,n_rayons), intent(out) :: xyz, uvw
 	integer :: iray
 	real(kind=dp) :: x0, y0, z0, u0, w0, v0, norme
 	real :: rand, rand2, rand3, W02, ARGMT, srW02

    do iray=rstart, rstart-1+n_rayons

    	if (etape==0) then
        	! Position = milieu de la cellule
        	x0 = r_grid(icell)
        	y0 = 0.0_dp
        	z0 = z_grid(icell)
        	if (lkeplerian) then
                       ! Direction verticale "z"
        		if (iray==1) then
            		w0=1.0_dp !nz
       			else
            		w0=-1.0_dp
        		endif
            	u0 = 0.0_dp !nx
            	v0 = 0.0_dp !ny
        	else
           		norme = sqrt(x0*x0 + y0*y0 + z0*z0)
       			if (iray==1) then
            		u0 = x0/norme !sin(theta)sin(phi) = nx
            		v0 = y0/norme !ny
            		w0 = z0/norme !nz = mu = cos(theta)
        		else
            		u0 = -x0/norme !-1  backward
            		v0 = -y0/norme
           			w0 = -z0/norme
           		endif
        	endif !lkeplerian

         else if (etape==1) then
            write(*,*) "step 1 not implemented yet"
            stop
            !random position + fixed directions

            !w0 = xmu(iray) * real(to_obs)
            !u0 = dsqrt(1.-xmu(iray)**2) * cos(pmu(iray)) * to_obs
            !v0 = dsqrt(1.-xmu(iray)**2) * sin(pmu(iray)) * to_obs

      	 else !etape 2
                   	    ! Position aleatoire dans la cellule
         	rand  = sprng(streami)
            rand2 = sprng(streami)
            rand3 = sprng(streami)

            CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
            rand = sprng(streami)
            W0 = 2.0_dp * rand - 1.0_dp !nz
            W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
            SRW02 = sqrt(W02)
            rand = sprng(streami)
            ARGMT = PI * (2.0_dp * rand - 1.0_dp)
            U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
            V0 = SRW02 * sin(ARGMT) !ny = sin(theta) * sin(phi)
		 end if !etape

		 xyz(1,iray) = x0; xyz(2,iray) = y0; xyz(3,iray) = z0
		 uvw(1,iray) = U0; uvw(2,iray) = V0; uvw(3,iray) = W0
		 !for the moment disable, integrale frequency weights are known, and for angular it is simply 1/N
		 !but probably the integral of line profile is not exactly 1 that's why we can recompute the normalisation
		 !CALL compute_integral_weight(id, icell, iray, n_rayons, x0, y0, z0, u0, v0, w0)
 	enddo
 RETURN
 END SUBROUTINE compute_directions
 
 SUBROUTINE calc_stellar_radiation(N,i_star,x,y,z,u,v,w,gamma)
 ! ---------------------------------------------------------------!
  ! Compute the stellar radiation field and in Istar:
  ! Radiation emerging from the star (at the surface), corrected
  ! by limb darkening.
  !
  ! tab_lambda has to be defined and repartition_energie_etoiles()
 ! -------------------------------------------------------------- !
  use input, only : limb_darkening, mu_limb_darkening
  use Planck, only : uLD, Bplanck
  integer, intent(in) :: N, i_star
  real(kind=dp), dimension(N), intent(out) :: gamma
  real(kind=dp), intent(in) :: u, v, w, x, y, z
  real(kind=dp) :: energie(N)
  real(kind=dp) :: mu, ulimb, LimbDarkening, surface, HC
  integer :: ns
  logical :: lintersect_spot

   gamma(:) = 1d0
   !cos(theta) = dot(r,n)/module(r)/module(n)
   mu = abs(x*u + y*v + z*w)/dsqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
   if (real(mu)>1d0) then !to avoid perecision error
    write(*,*) "mu=",mu, x, y, z, u, v, w
    CALL Error(" mu limb > 1!")
   end if

   !if not using mcfost Energie, Istar is computed after the wavelength grid is done
   !at the centre of the disk for each star, so that it is need to apply only a correction factor
   !for LD and spots
   
   !Re-norm E_stars(:)*Prob_E_Star(:,i_star) at the stellar surface
   !surface = (etoile(i_star)%r * AU_to_Rsun)**2

   !1) Get the energy radiated by the star at the stellar surface
   !I need unit of I which is the unit of Bnu = W/m2/Hz/sr
   !Here, E_stars in W/m2. But I is in W/m2/Hz/sr unit of Bnu.
   !E_stars propto Blambda*lambda; remembering, Bnu = Blambda*c/nu**2 = Blambda*lambda**2/c
   ! we get E_stars (in W/m2/Hz) = E_stars / lambda * lambda**2 / c
!    energie(:) = Prob_E_Star(:,i_star) * E_stars(:) * tab_lambda(:) * 1.0e-6 / surface &
!              * 1.35e-12 * (tab_lambda(:) * 1d-6) / C_LIGHT
   !write(*,*) maxval(energie)
   !write(*,*) maxval(E_stars* (tab_lambda(:) * 1.0e-6)/CLIGHT)
   !!CALL Bplanck(etoile(i_star)%T*1d0, energie) !it is not factorised for test cheks, but can be computed outside loop
   !write(*,*) maxval(energie)
   !stop

   !2) Correct with the contrast gamma of a hotter/cooler region if any
   CALL intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect_spot)
   if (lintersect_spot) then
     gamma(:) = (dexp(hc_k/NLTEspec%lambda(:)/real(etoile(i_star)%T,kind=dp))-1)/&
     			(dexp(hc_k/NLTEspec%lambda(:)/etoile(i_star)%SurfB(ns)%T)-1)
     !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
   end if

   !3) Apply Limb darkening
   if (llimb_darkening) then
     !!LimbDarkening = Interp1D(real(mu_limb_darkening,kind=dp),real(limb_darkening,kind=dp),mu)
     !LimbDarkening = interp_dp(limb_darkening, mu_limb_darkening, mu)
     !pol not included yet
     stop
   else
     !write(*,*) maxval(uLD(real(etoile(i_star)%T,kind=dp))), minval(uLD(real(etoile(i_star)%T,kind=dp)))
     ulimb = 0.0 ! could use BB slope
     LimbDarkening = 1d0 - ulimb*(1d0-mu)
   end if
   !Istar(:) = energie(:) * LimbDarkening * gamma(:)
   gamma(:) = LimbDarkening * gamma(:)

 RETURN
 END SUBROUTINE calc_stellar_radiation

END MODULE solve_transfer