! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !
! This module solves for the radiative transfer equation by ray-tracing for
! multi-level atoms
! ----------------------------------------------------------------------------------- !
! ----------------------------------------------------------------------------------- !

MODULE AtomicTransfer

 use metal, only                        : metal_bb_new, compute_opacities, alloc_atom_quantities, dealloc_Atom_quantities
 use opacity
 use Planck, only 						: bplanck
 use Profiles, only 					: Iprofile, Iprofile_thomson, Iprofile_cmf_to_obs, Zprofile
 !use broad, only 						: Damping, RadiativeDamping !for tests
 use spectrum_type
 use atmos_type
 use readatom
 use lte
 use constant, only 					: MICRON_TO_NM
 use collision, only					: CollisionRate !future deprecation
 use impact
 use solvene
 use statequil_atoms
 use init_solution, only 			: Init_NLTE, free_NLTE_sol, gpop_old, pop_old, flatpops, lcell_converged
 use accelerate
 use voigtfunctions
 use writeatom, only 					: writePops, writeelectron, writehydrogendensity
 use write_opacity!, only 				: write_Jnu, write_taur, write_atom_xsections_bf_ff
 use math
 !$ use omp_lib

 !MCFOST's original modules
 use input, only 						: lkeplerian, linfall, RT_line_method
 use parametres
 use grid
 !use dust_prop
 use dust_transfer, only 				: compute_stars_map
 use dust_ray_tracing, only 			: init_directions_ray_tracing ,            &
                              			  tab_u_RT, tab_v_RT, tab_w_RT, tab_RT_az, &
                              			  tab_RT_incl, stars_map, kappa, stars_map_cont
 use stars
 use wavelengths
 !use density
 use mcfost_env, only : dp
 use constantes, only : tiny_dp, huge_dp

 IMPLICIT NONE
 
 real(kind=dp), parameter :: tiny_chi = 1d-50
 !Pointer to formal solver
 PROCEDURE(INTEG_RAY_LINE_I), pointer :: INTEG_RAY_LINE => NULL()
 !Temporary variable for Zeeman calculations
 real(kind=dp), dimension(:,:), allocatable :: QUV
 !Temporary variables for Contribution functions
 real(kind=dp), allocatable :: S_contrib(:,:,:), S_contrib2(:,:), chil(:), Sl(:)

 CONTAINS

 SUBROUTINE INTEG_RAY_LINE_I_NLTE(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray
  ! co continuum
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before, etau, edtau
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, chiI, LimbD
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, dtau, dtau_ds
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau(:) = 0.0_dp
  chiI(:) = 0.0_dp
  NLTEspec%I(:,iray,id) = 0d0

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
     lcellule_non_vide = (atmos%icompute_atomRT(icell) > 0)
     if (atmos%icompute_atomRT(icell) < 0) RETURN !-1 if dark
    else
     lcellule_non_vide=.false.
    endif
    
    
    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN

    if (lintersect_stars) then
      if (icell == icell_star) then
       call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)
       !CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
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
     if ((nbr_cell == 1).and.labs) ds(iray,id) = l * AU_to_m

     eval_operator = (labs .and. (nbr_cell == 1)) 
     
     !init bound-bound eta and chi
     CALL initAtomOpac(id)
     !compute LTE bound-bound opacity
      if (atmos%NpassiveAtoms>0) CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)

      chiI(:)  = NLTEspec%AtomOpac%chi_p(:,id) + NLTEspec%AtomOpac%Kc(:,icell) + tiny_chi
      Snu(:) = 0d0 !for this cell

      if (atmos%NactiveAtoms>0) then 
       !init NLTE bound-bound
       CALL initAtomOpac_nlte(id)
       CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)
       if (eval_operator) then
        chi_loc(:,iray,id) = chiI(:) !for this cell, for computing eta_Atom_loc = chi_loc+chi_nlte
       endif

       !add NLTE bound-bound and bound-free opacities
       chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id) +  NLTEspec%AtomOpac%Kc_nlte(:,icell) 
       
       !Source function
       Snu = Snu + (NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell)) / chiI(:)
      
        if (atmos%electron_scattering) then !coherent scattering approximation
         Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
        endif
      endif !active atoms
                    
      dtau(:)   = l_contrib * chiI(:)

      !add LTE part of the source function
      Snu = Snu + ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
      
      !compute intensity at this point of the propagation
      
      NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
     
     !Compute (1-exp(-dtau*chi))
     if (eval_operator) then
       dtau_ds = ds(iray,id) * chiI(:) !could be avoid id ds = l_contrib, it is simply dtau
       CALL calc_psi_operator(id, icell, iray, chiI, dtau_ds)
     end if

     tau = tau + dtau

    end if 
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_I_NLTE
  
 SUBROUTINE INTEG_RAY_LINE_I(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray
  ! crossing different cells.
  ! --> Atomic Lines case.
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
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

  tau(:) = 0.0_dp
  tau_c(:) = 0.0_dp
  chiI(:) = 0.0_dp
  chiIc(:) = 0.0_dp


  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0

    ! write(*,*) "Before propagating", "x=",x,"y=",y,"z=",z,"x1=",x1,"y1=",y1,"z1=",z1

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
    !write(*,*) "Boucle infinie, icell=", icell

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
       !AT the moment the stellar Flux is only a BB
       call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)
       !CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau_c)
       RETURN
      end if
    endif

    nbr_cell = nbr_cell + 1

    ! Calcul longeur de vol et profondeur optique dans la cellule
    previous_cell = 0 ! unused, just for Voronoi
     !write(*,*) "in integ_ray_line before X-cell:", icell, "x=",x0,"y=",y0,"z=",z0,"x1=",x1,"y1=",y1,"z1=",z1
    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell,l, l_contrib, l_void_before)
     !write(*,*) "in integ_ray_line after X-cell:", icell, "x=",x0,"y=",y0,"z=",z0,"x1=",x1,"y1=",y1,"z1=",z1


    !count opacity only if the cell is filled, else go to next cell
    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs) !not used yet
     ! opacities in m^-1
     l_contrib = l_contrib * AU_to_m !l_contrib in m
     if ((nbr_cell == 1).and.labs) ds(iray,id) = l * AU_to_m

     eval_operator = (labs .and. (nbr_cell == 1)) !labs if false for images
     											  !so no pb if Nact>0 and we use a different grid
     CALL initAtomOpac(id) !in principle, if lstore_opac is the only mode, this init could be in 
     						!the condition below.
 	!write(*,*) "Taille cellule", r_lim(icell)
      if (atmos%NpassiveAtoms>0) CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)

      chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell) + tiny_chi 
      chiI(:)  = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)

      if (atmos%NactiveAtoms>0) then 
       CALL initAtomOpac_nlte(id)
      
       CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, eval_operator)

       chiIc(:) = chiIc(:) + NLTEspec%AtomOpac%Kc_nlte(:,icell)

       chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id) +  NLTEspec%AtomOpac%Kc_nlte(:,icell)
      endif !active atoms

      dtau(:)   = l_contrib * chiI(:)
      dtau_c(:) = l_contrib * chiIc(:)
            
      Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
      Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
      
      if (atmos%Nactiveatoms>0) then 
        Snu = Snu + ( NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell) ) / chiI(:)
        Snu_c = Snu_c + NLTEspec%AtomOpac%jc_nlte(:,icell) / chiIc(:)
      endif
      if (atmos%electron_scattering) then
         Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
         Snu_c = Snu_c + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) / chiIc(:)
      endif

      NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + dexp(-tau) * (1.0_dp - dexp(-dtau)) * Snu
      NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + dexp(-tau_c) * (1.0_dp - dexp(-dtau_c)) * Snu_c
     
     if (eval_operator) then
       dtau_ds = ds(iray,id) * chiI(:)
       CALL calc_psi_operator(id, icell, iray, chiI, dtau_ds)
     end if
     
     tau = tau + dtau
     tau_c = tau_c + dtau_c

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_I

 SUBROUTINE INTEG_RAY_LINE_Z(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! This routine is "LTE", although NLTE pops can be used
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, LimbD
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, chiI, chiIc
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  real(kind=dp) :: facteur_tau
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  tau = 0.0_dp
  tau_c = 0.0_dp


  NLTEspec%I(:,iray,id) = 0d0
  NLTEspec%Ic(:,iray,id) = 0d0
  NLTEspec%StokesV(:,iray,id) = 0d0
  NLTEspec%StokesQ(:,iray,id) = 0d0
  NLTEspec%StokesU(:,iray,id) = 0d0

  call intersect_stars(x,y,z, u,v,w, lintersect_stars, i_star, icell_star)

  infinie : do
    icell = next_cell
    x0=x1 ; y0=y1 ; z0=z1

    if (icell <= n_cells) then
     lcellule_non_vide = (atmos%icompute_atomRT(icell)>0)
    else
     lcellule_non_vide=.false.
    endif


    if (test_exit_grid(icell, x0, y0, z0)) RETURN
    if (lintersect_stars) then
      if (icell == icell_star) then
       call calc_stellar_surface_brightness(NLTEspec%Nwaves,NLTEspec%lambda,i_star, x0, y0, z0, u,v,w,LimbD)
       !CALL calc_stellar_radiation(NLTEspec%Nwaves,i_star, x0, y0, z0, u,v,w,LimbD)
       NLTEspec%I(:,iray, id) =  NLTEspec%I(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau)
       NLTEspec%Ic(:,iray, id) =  NLTEspec%Ic(:,iray, id) + LimbD(:) * NLTEspec%Istar(:,i_star)*dexp(-tau_c)
       RETURN
       endif
    endif

    nbr_cell = nbr_cell + 1

    previous_cell = 0

    call cross_cell(x0,y0,z0, u,v,w,  icell, previous_cell, x1,y1,z1, next_cell, &
                     l, l_contrib, l_void_before)



    if (lcellule_non_vide) then
     lsubtract_avg = ((nbr_cell == 1).and.labs)
     l_contrib = l_contrib * AU_to_m
     if ((nbr_cell == 1).and.labs)  ds(iray,id) = l * AU_to_m

     CALL initAtomOpac(id)
     CALL initAtomOpac_zeeman(id)

      if (atmos%NpassiveAtoms>0) CALL Metal_bb_new(id, icell, x0, y0, z0, x1, y1, z1, u, v, w, l)

      chiIc(:) = NLTEspec%AtomOpac%Kc(:,icell) + tiny_chi 
      if (atmos%NactiveAtoms > 0) chiIc(:) = chiIc(:) + NLTEspec%AtomOpac%Kc_nlte(:,icell)
      chiI(:) = NLTEspec%AtomOpac%chi_p(:,id) + chiIc(:)

      if (atmos%NactiveAtoms>0) then 
       CALL initAtomOpac_nlte(id)

       CALL NLTE_bound_bound(id, icell, iray, x0, y0, z0, x1, y1, z1, u, v, w, l, .false.)

       chiI(:) = chiI(:) + NLTEspec%AtomOpac%chi(:,id)
      endif
                    
      dtau(:)   = l_contrib * chiI(:)
      dtau_c(:) = l_contrib * chiIc(:)
      
      Snu = ( NLTEspec%AtomOpac%jc(:,icell) + NLTEspec%AtomOpac%eta_p(:,id) ) / chiI(:)
      Snu_c = NLTEspec%AtomOpac%jc(:,icell) / chiIc(:)
      
      !add only if non zero,; this should be fater to do than adding a large 0 array
      if (atmos%Nactiveatoms>0) then 
        Snu = Snu + (NLTEspec%AtomOpac%eta(:,id) + NLTEspec%AtomOpac%jc_nlte(:,icell)) / chiI(:)
        Snu_c = Snu_c + NLTEspec%AtomOpac%jc_nlte(:,icell) / chiIc(:)
      endif
       if (atmos%electron_scattering) then !if pure LTE J computed before emission_line_map.
         Snu = Snu + NLTEspec%Jc(:,icell) * NLTEspec%AtomOpac%sca_c(:,icell) / chiI(:)
         Snu_c = Snu_c + NLTEspec%AtomOpac%sca_c(:,icell) * NLTEspec%Jc(:,icell) / chiIc(:)
       endif

    !continuum not affected by polarisation yet
     NLTEspec%Ic(:,iray,id) = NLTEspec%Ic(:,iray,id) + exp(-tau_c) * (1.0_dp - exp(-dtau_c)) * Snu_c

    !Correct line source fonction from polarisation
    !explicit product of Seff = S - (K/chiI - 1) * I
    !Is there a particular initialization to do for polarisation ?
    !because it will be zero at the first place
     Snu(:) = Snu(:) -NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI *  NLTEspec%StokesQ(:,iray,id) - &
          NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI *  NLTEspec%StokesU(:,iray,id) - &
          NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI *  NLTEspec%StokesV(:,iray,id)

     NLTEspec%I(:,iray,id) = NLTEspec%I(:,iray,id) + exp(-tau) * (1.0_dp - exp(-dtau)) * Snu

     NLTEspec%StokesQ(:,iray,id) = NLTEspec%StokesQ(:,iray,id) + &
                             exp(-tau) * (1.0_dp - exp(-dtau)) * (&
			NLTEspec%AtomOpac%etaQUV_p(:,1,id) /chiI - &
			NLTEspec%AtomOpac%chiQUV_p(:,1,id) / chiI * NLTEspec%I(:,iray,id) - &
			NLTEspec%AtomOpac%rho_p(:,3,id)/chiI * NLTEspec%StokesU(:,iray,id) + &
			NLTEspec%AtomOpac%rho_p(:,2,id)/chiI * NLTEspec%StokesV(:,iray, id) &
     ) !end Snu_Q
     NLTEspec%StokesU(:,iray,id) = NLTEspec%StokesU(:,iray,id) + &
                             exp(-tau) * (1.0_dp - exp(-dtau)) * (&
			NLTEspec%AtomOpac%etaQUV_p(:,2,id) /chiI - &
			NLTEspec%AtomOpac%chiQUV_p(:,2,id) / chiI * NLTEspec%I(:,iray,id) + &
			NLTEspec%AtomOpac%rho_p(:,3,id)/chiI * NLTEspec%StokesQ(:,iray,id) - &
			NLTEspec%AtomOpac%rho_p(:,1,id)/chiI * NLTEspec%StokesV(:,iray, id) &
     ) !end Snu_U
     NLTEspec%StokesV(:,iray,id) = NLTEspec%StokesV(:,iray,id) + &
                             exp(-tau) * (1.0_dp - exp(-dtau)) * (&
			NLTEspec%AtomOpac%etaQUV_p(:,3,id) /chiI - &
			NLTEspec%AtomOpac%chiQUV_p(:,3,id) / chiI * NLTEspec%I(:,iray,id) - &
			NLTEspec%AtomOpac%rho_p(:,2,id)/chiI * NLTEspec%StokesQ(:,iray,id) + &
			NLTEspec%AtomOpac%rho_p(:,1,id)/chiI * NLTEspec%StokesU(:,iray, id) &
     ) !end Snu_V

     facteur_tau = 1.0

     tau = tau + dtau * facteur_tau
     tau_c = tau_c + dtau_c

    end if
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_LINE_Z

  SUBROUTINE FLUX_PIXEL_LINE(&
         id,ibin,iaz,n_iter_min,n_iter_max,ipix,jpix,pixelcorner,pixelsize,dx,dy,u,v,w)
  ! -------------------------------------------------------------- !
   ! Computes the flux emerging out of a pixel.
   ! see: mol_transfer.f90/intensite_pixel_mol()
  ! -------------------------------------------------------------- !

   integer, intent(in) :: ipix,jpix,id, n_iter_min, n_iter_max, ibin, iaz
   real(kind=dp), dimension(3), intent(in) :: pixelcorner,dx,dy
   real(kind=dp), intent(in) :: pixelsize,u,v,w
   integer, parameter :: maxSubPixels = 32
   real(kind=dp) :: x0,y0,z0,u0,v0,w0
   real(kind=dp), dimension(NLTEspec%Nwaves) :: Iold, I0, I0c !, nu
   real(kind=dp), dimension(3) :: sdx, sdy
   real(kind=dp):: npix2, diff, normF, R0
   real(kind=dp), parameter :: precision = 1.e-2
   integer :: i, j, subpixels, iray, ri, zj, phik, icell, iter
   logical :: lintersect, labs

   labs = .false.
   ! Ray tracing : on se propage dans l'autre sens
   u0 = -u ; v0 = -v ; w0 = -w

   ! le nbre de subpixel en x est 2^(iter-1)
   subpixels = 1
   iter = 1
   diff = 0.
   
   infinie : do ! Boucle infinie tant que le pixel n'est pas converge
     npix2 =  real(subpixels)**2
     Iold = I0
     I0 = 0d0
     I0c = 0d0
     if (atmos%magnetized.and. PRT_SOLUTION == "FULL_STOKES") QUV(:,:) = 0d0 !move outside
     !if (lcontrib_function) S_contrib2(:,:) = 0d0

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
           if (lintersect) then ! On rencontre la grille, on a potentiellement du flux
             CALL INTEG_RAY_LINE(id, icell, x0,y0,z0,u0,v0,w0,iray,labs)

             I0 = I0 + NLTEspec%I(:,iray,id) / npix2
             I0c = I0c + NLTEspec%Ic(:,iray,id) / npix2
             
             if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
             	QUV(3,:) = QUV(3,:) + NLTEspec%STokesV(:,iray,id)
             	QUV(1,:) = QUV(1,:) + NLTEspec%STokesQ(:,iray,id)
             	QUV(2,:) = QUV(2,:) + NLTEspec%STokesU(:,iray,id)
             end if
             !if (lcontrib_function) S_contrib2(:,:) = S_contrib2(:,:) + S_contrib(:,:,id)

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
         !write(*,*) "Warning : converging pb in ray-tracing"
         !write(*,*) " Pixel", ipix, jpix
        if (diff > 1) write(*,*) 'pixel not converged:', iter, subpixels, diff
        exit infinie
     else
        ! On fait le test sur a difference
        diff = maxval( abs(I0 - Iold) / (I0 + 1e-300_dp) )
        ! There is no iteration for Q, U, V, assuming that if I is converged, then Q, U, V also.
        ! Can be added and then use diff as max(diff, diffQ, diffU, diffV)
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

  !Prise en compte de la surface du pixel (en sr)

  ! Flux out of a pixel in W/m2/Hz/pix
  normF = (pixelsize / (distance*pc_to_AU) )**2

  ! adding to the total flux map.
  if (RT_line_method==1) then
    NLTEspec%Flux(:,1,1,ibin,iaz) = NLTEspec%Flux(:,1,1,ibin,iaz) + I0 * normF
    NLTEspec%Fluxc(:,1,1,ibin,iaz) = NLTEspec%Fluxc(:,1,1,ibin,iaz) + I0c * normF
    if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
     NLTEspec%F_QUV(1,:,1,1,ibin,iaz) = NLTEspec%F_QUV(1,:,1,1,ibin,iaz)+QUV(1,:) * normF !U
     NLTEspec%F_QUV(2,:,1,1,ibin,iaz) = NLTEspec%F_QUV(2,:,1,1,ibin,iaz)+QUV(2,:) * normF !Q
     NLTEspec%F_QUV(3,:,1,1,ibin,iaz) = NLTEspec%F_QUV(3,:,1,1,ibin,iaz)+QUV(3,:) * normF !V
    end if
  else
    NLTEspec%Flux(:,ipix,jpix,ibin,iaz) = I0 * normF
    NLTEspec%Fluxc(:,ipix,jpix,ibin,iaz) = I0c * normF
    if (atmos%magnetized.and.PRT_SOLUTION == "FULL_STOKES") then
     NLTEspec%F_QUV(1,:,ipix,jpix,ibin,iaz) = QUV(1,:) * normF !U
     NLTEspec%F_QUV(2,:,ipix,jpix,ibin,iaz) = QUV(2,:) * normF !Q
     NLTEspec%F_QUV(3,:,ipix,jpix,ibin,iaz) = QUV(3,:) * normF !V
    end if
  end if
  !if (lcontrib_function) NLTEspec%Ksi(:,:,ibin,iaz) = &
  !                       NLTEspec%Ksi(:,:,ibin,iaz) + S_contrib2(:,:) * (pixelsize*au_to_m)**2 !In W/Hz


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
  real(kind=dp) :: x0,y0,z0,l,u,v,w

  real(kind=dp), dimension(3) :: uvw, x_plan_image, x, y_plan_image, center, dx, dy, Icorner
  real(kind=dp), dimension(3,nb_proc) :: pixelcorner
  real(kind=dp):: taille_pix, nu
  integer :: i,j, id, npix_x_max, n_iter_min, n_iter_max

  integer, parameter :: n_rad_RT = 600, n_phi_RT = 100 !(100, 36)
  real(kind=dp), dimension(n_rad_RT) :: tab_r
  real(kind=dp):: rmin_RT, rmax_RT, fact_r, r, phi, fact_A, cst_phi
  integer :: ri_RT, phi_RT, lambda, lM, lR, lB, ll
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

  ! position initiale hors modele (du cote de l'observateur)
  ! = centre de l'image
  l = 10.*Rmax  ! on se met loin ! in AU

  x0 = u * l  ;  y0 = v * l  ;  z0 = w * l
  center(1) = x0 ; center(2) = y0 ; center(3) = z0

  ! Coin en bas gauche de l'image
  Icorner(:) = center(:) - 0.5 * map_size * (x_plan_image + y_plan_image)
  
  if (RT_line_method==1) then !log pixels
    n_iter_min = 1
    n_iter_max = 1

    ! dx and dy are only required for stellar map here
    taille_pix = (map_size/zoom)  ! en AU
    dx(:) = x_plan_image * taille_pix
    dy(:) = y_plan_image * taille_pix

    i = 1
    j = 1
    lresolved = .false.

    rmin_RT = max(w*0.9_dp,0.05_dp) * Rmin
    rmax_RT = Rmax !* 2.0_dp

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
    !$omp shared(n_iter_min,n_iter_max,l_sym_ima,cst_phi,ibin,iaz,etoile)
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
           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        end do !j
     end do !i
     !$omp end do
     !$omp end parallel
  else !method 2
     ! Vecteurs definissant les pixels (dx,dy) dans le repere universel
     taille_pix = (map_size/zoom) / real(max(npix_x,npix_y),kind=dp) ! en AU
     lresolved = .true.

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
     n_iter_min = 1 !1 !3
     n_iter_max = 1 !1 !6
     !$omp do schedule(dynamic,1)
     do i = 1,npix_x_max
     !do i=npix_x_max/2+1, npix_x_max/2+1
        !$ id = omp_get_thread_num() + 1

        do j = 1,npix_y
        !do j=npix_y/2+1,npix_y/2+1
        !write(*,*) "ipix, jpix", i, j
           ! Coin en bas gauche du pixel
           pixelcorner(:,id) = Icorner(:) + (i-1) * dx(:) + (j-1) * dy(:)
           CALL FLUX_PIXEL_LINE(id,ibin,iaz,n_iter_min,n_iter_max, &
                      i,j,pixelcorner(:,id),taille_pix,dx,dy,u,v,w)
        end do !j
     end do !i
     !$omp end do
     !$omp end parallel
  end if

 RETURN
 END SUBROUTINE EMISSION_LINE_MAP

 SUBROUTINE Atomic_transfer()
 ! --------------------------------------------------------------------------- !
  ! This routine initialises the necessary quantities for atomic line transfer
  ! and calls the appropriate routines for LTE or NLTE transfer.
 ! --------------------------------------------------------------------------- !
  integer :: atomunit = 1, nact
  integer :: icell, m
  integer :: ibin, iaz
  integer, parameter :: Nrayone = 1
  character(len=20) :: ne_start_sol = "H_IONISATION"
  character(len=20)  :: newPRT_SOLUTION = "FULL_STOKES"
  
 ! -------------------------------INITIALIZE AL-RT ------------------------------------ !
  !only one available yet, I need one unpolarised, faster and more accurate.
  Voigt => VoigtHumlicek
  !Voigt => dirac_line
  !Profile => Iprofile_cmf_to_obs
  Profile => IProfile
  !Profile => Iprofile_thomson

  lstore_opac = .true. !futur deprecation always true but ATM needed

  if (npix_x_save > 1) then
   RT_line_method = 2 ! creation d'une carte avec pixels carres
   npix_x = npix_x_save ; npix_y = npix_y_save
  else
   RT_line_method = 1 !pixels circulaires
  end if


  if (PRT_SOLUTION == "FULL_STOKES") then
   CALL Warning(" Full Stokes solution not allowed. Stokes polarization not handled in SEE yet.")
   PRT_SOLUTION = "FIELD_FREE"
  end if

  if (atmos%magnetized .and. PRT_SOLUTION /= "NO_STOKES") then
   if (PRT_SOLUTION == "FIELD_FREE") newPRT_SOLUTION = "FULL_STOKES"
   CALL adjustStokesMode(PRT_SOLUTION)
  end if

 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
        !! ----------------------- Read Model ---------------------- !!

  if (.not.lpluto_file) then
   if (lmodel_ascii) then
    CALL readAtmos_ascii(density_file)
    !CALL writeHydrogenDensity()
    !CALL writeTemperature()
   end if
  end if
        !! --------------------------------------------------------- !!
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------- READATOM and INITIZALIZE POPS ------------------------------ !

  CALL readAtomicModels(atomunit)
  
  if (atmos%NactiveAtoms > 0) then 
   atmos%Nrays = 200
   !save time by not computing continuum flux, just needed for image
   !INTEG_RAY_LINE => INTEG_RAY_LINE_I_NLTE
   INTEG_RAY_LINE => INTEG_RAY_LINE_I
  else
   atmos%Nrays = Nrayone
   if (lelectron_scattering) atmos%Nrays = 600 !here before I is allocated
   INTEG_RAY_LINE => INTEG_RAY_LINE_I
  endif

  !compute first guess of electron density ??
  if (.not.atmos%calc_ne) atmos%calc_ne = lsolve_for_ne
  if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
  if (atmos%calc_ne) then
   if (lsolve_for_ne) write(*,*) "(Force) Solving for electron density"
   if (.not.lsolve_for_ne) write(*,*) "Solving for electron density"
   write(*,*) " Starting solution : ", ne_start_sol
   if ((ne_start_sol == "NE_MODEL") .and. (atmos%calc_ne)) then
    write(*,*) "WARNING, ne from model is presumably 0 (or not given)."
    write(*,*) "Changing ne starting solution NE_MODEL in H_IONISATION"
    ne_start_sol = "H_IONISATION"
   end if
   CALL SolveElectronDensity(ne_start_sol)
   CALL writeElectron()
  end if

  CALL setLTEcoefficients () !write pops at the end because we probably have NLTE pops also

 ! ------------------------------------------------------------------------------------ !
 ! ---------- INITIALIZE WAVELNGTH, BACKGROUND OPAC AND ATOMC QUANTITIES -------------- !
 ! ------------------------------------------------------------------------------------ !
  atmos%electron_scattering=lelectron_scattering !but also H and He, futur deprec of atmos
  if (ltab_wavelength_image) NLTEspec%write_wavelength_grid = .true.
  CALL initSpectrum(vacuum_to_air=lvacuum_to_air)

  if ((atmos%NactiveAtoms>0) .or. .not.(ltab_wavelength_image)) then
   if (n_etoiles > 0) CALL init_stellar_disk !for all wavelengths, all stars at disk centre
   write(*,*) " Computing background opacities..."
   CALL alloc_atom_quantities
   CALL compute_opacities
   write(*,*) " ..done"
  endif
 ! ------------------------------------------------------------------------------------ !
 ! --------------------- ADJUST MCFOST FOR STELLAR MAP AND VELOCITY ------------------- !
  ! ----- ALLOCATE SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT ------!
  CALL reallocate_mcfost_vars() !assumes more than 1 wavelength otherwise delta_wl is 0!
  ! --- END ALLOCATING SOME MCFOST'S INTRINSIC VARIABLES NEEDED FOR AL-RT --!
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------------- Jny ----------------------------------------- !
   !IF NLTE loop we should also be able to read or compute an estimate of Jcont by the way...
   ! -> TBD
   if (atmos%NactiveAtoms == 0 .and. lelectron_scattering) then
    call alloc_jnu()
!    if (.not.lread_jnu_atom) then
	if (lread_jnu_atom) then
		call read_jnu_ascii
!      call read_jnu() !I now write Jnu in ascii file to be interpolated if lread_jnu_atom
!	   the interpolated version of Jnu is written at the end in write_jnu()
    else
     call iterate_Jnu()
     call write_Jnu
     if (lstop_after_jnu) then
      write(*,*) " Jnu calculation done." 
      stop
     endif
    endif
   endif
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- NLTE LOOP -------------------------------------- !
  if (atmos%Nactiveatoms > 0) then
     CALL NLTEloop()
     !free not useful quantities
     do icell=1,atmos%Nactiveatoms
      deallocate(atmos%ActiveAtoms(icell)%ptr_atom%etac)
      !Collision etc
     enddo
  endif
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------- WRITE CONVERGED POPULATIONS ------------------------------ !
 ! ------------------------------------------------------------------------------------ !
  !! alternatively, could be invoked earlier, if the system send a message to stop the code
  !! before convergence, so that we restart with these pops.
  !! Or define a function to store pops every N iter.
  !!do icell=1,atmos%Natom
   !CALL writePops(atmos%Atoms(icell)%ptr_atom)
  !!end do
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! ----------------------------------- MAKE IMAGES ------------------------------------ !
 ! ------------------------------------------------------------------------------------ !

   !should add an other flag to force to avoid continua lines
   !Define a wavelength grid for image with only lines, if not input wavelength
!    if (.not.(ltab_wavelength_image)) then !or for all passive lines also
!     ltab_wavelength_image = .true.
!     !!CALL write_wavelengths_table_NLTE_lines(NLTEspec%lambda) !!after grid creation actually
!     tab_wavelength_image = "line_waves.s"
!    endif
  
  !not needed, but avoid to reset it, if only lte, it is already integ_ray_line_i
  if (atmos%NactiveAtoms > 0) then
   INTEG_RAY_LINE => NULL()
   INTEG_RAY_LINE => INTEG_RAY_LINE_I
  endif

  if (ltab_wavelength_image) then
   NLTEspec%atmos%Nrays = Nrayone

   CALL initSpectrumImage() !deallocate waves arrays/ define a new grid also resample Jnu if any
   if (n_etoiles > 0) CALL init_stellar_disk 
   write(*,*) " Computing continuum opacities for image..."
   if (atmos%NactiveAtoms >0) CALL dealloc_atom_quantities
   CALL alloc_atom_quantities
   CALL compute_opacities !recompute background opac
   if (atmos%NactiveAtoms > 0) then
    CALL compute_nlte_bound_free
    !! can be added to Kc and jc to save memory for image
    !NLTEspec%AtomOpac%Kc(:,:,1) = NLTEspec%AtomOpac%Kc(:,:,1) + NLTEspec%AtomOpac%Kc_nlte(:,:)
    !NLTEspec%AtomOpac%jc(:,:,1) = NLTEspec%AtomOpac%jc(:,:,1) + NLTEspec%AtomOpac%jc_nlte(:,:)
	!deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)
   endif
   write(*,*) " ..done"
   !add NLTE contiuna
   !!CALL reallocate_mcfost_vars() !wavelength only?
   CALL reallocate_mcfost_wavelength_arrays()
  else
    !same wave grid, Jnu the same if any
    !if (atmos%NactiveAtoms > 0) then
     !NLTEspec%AtomOpac%Kc(:,:,1) = NLTEspec%AtomOpac%Kc(:,:,1) + NLTEspec%AtomOpac%Kc_nlte(:,:)
     !NLTEspec%AtomOpac%jc(:,:,1) = NLTEspec%AtomOpac%jc(:,:,1) + NLTEspec%AtomOpac%jc_nlte(:,:)
	 !deallocate(NLTEspec%AtomOpac%Kc_nlte, NLTEspec%AtomOpac%jc_nlte)
	!endif
    CALL reallocate_rays_arrays(Nrayone) !rellocate rays array, but keep the same wavelength grid
  end if

  if (atmos%magnetized .and. PRT_SOLUTION /= "NO_STOKES") then
   if (PRT_SOLUTION == "FIELD_FREE") PRT_SOLUTION = newPRT_SOLUTION
    CALL adjustStokesMode(PRT_SOLUTION)
    allocate(QUV(3, NLTEspec%Nwaves))
    QUV = 0d0
  end if

  CALL alloc_flux_image()
  write(*,*) "Computing emission flux map..."

  !Use converged NLTEOpac
!   if (lcontrib_function) then
!    write(*,*) "   -> including contribution functions."
!    allocate(S_contrib(NLTEspec%Nwaves,atmos%Nspace,NLTEspec%NPROC))
!    allocate(chil(NLTEspec%Nwaves), Sl(NLTEspec%Nwaves),S_contrib2(NLTEspec%Nwaves,n_cells))
!    S_contrib(:,:,:) = 0d0; S_contrib2(:,:) = 0d0
!    INTEG_RAY_LINE => NULL()
!    INTEG_RAY_LINE => INTEG_RAY_LINE_I_CNTRB !same as INTEG_RAY_LINE_I, without nlte (but %n is known)
!   end if

  !Actual flux calculation
  do ibin=1,RT_n_incl
     do iaz=1,RT_n_az
       CALL EMISSION_LINE_MAP(ibin,iaz)
     end do
  end do
  write(*,*) " ..done"
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------- WRITE RESULTS -------------------------------------- !
 ! ------------------------------------------------------------------------------------ !
  write(*,*) "Writing result to file..."
  call WRITE_FLUX_ASCII() !for testing
  if (2*size(NLTEspec%Flux)/(1024.**3) <= 6.) call WRITE_FLUX
  deallocate(NLTEspec%Flux, NLTEspec%Fluxc) !free here to release a bit of memory
    
	call write_contrib_lambda_ascii(654.0_dp,0.0_dp,0.0_dp,1.0_dp)
	call write_contrib_lambda_ascii(1083.0_dp,0.0_dp,0.0_dp,1.0_dp)
	call write_taur(500._dp,0._dp,0._dp,1.0_dp)
!-> building
!  call write_contrib_lambda(654.0_dp,0.0_dp,0.0_dp,1.0_dp)
!   call write_contrib_lambda(1080._dp,0._dp,0._dp,0._dp,0._dp,0._dp,0._dp)

  !write it after because sca_c (and maybe %Kc, %jc, %Kc_nlte, %Jc_nlte) is modified !!
  do icell=1, atmos%Natom
    CALL write_cont_opac_ascii(atmos%Atoms(icell)%ptr_atom)
    !!CALL write_cont_opac(atmos%Atoms(icell)%ptr_atom) !-> building
    CALL write_atom_xsections_bf_ff(atmos%Atoms(icell)%ptr_atom)
  enddo

  if (atmos%electron_scattering) then
   !Jnu is written to ascii file if read
   if ((atmos%NactiveAtoms > 0).or.(atmos%NactiveAtoms==0.and.lread_jnu_atom)) call write_Jnu
   CALL dealloc_Jnu()
  endif
!   if (lcontrib_function) then
!    CALL WRITE_CNTRB_FUNC_pix()
!    deallocate(S_contrib, chil, Sl, S_contrib2)
!  end if
 ! ------------------------------------------------------------------------------------ !
 ! ------------------------------------------------------------------------------------ !
 ! -------------------------------- CLEANING ------------------------------------------ !
  if (allocated(QUV)) deallocate(QUV)
 !close file after NLTE loop
!Temporary: because we kept in memory, so file is closed earlier
!  do nact=1,atmos%Nactiveatoms
!   CALL closeCollisionFile(atmos%ActiveAtoms(nact)%ptr_atom) !if opened
!  end do
 !CALL WRITEATOM() !keep C in memory for that ?
 CALL freeSpectrum() !deallocate spectral variables
 CALL free_atomic_atmos()
 NULLIFY(Profile, Voigt, INTEG_RAY_LINE)
 !CALL dealloc_profile_line_profiles

 RETURN
 END SUBROUTINE
 ! ------------------------------------------------------------------------------------ !

 SUBROUTINE NLTEloop() !for all active atoms
 ! -------------------------------------------------------- !
  ! Descriptor here
 ! -------------------------------------------------------- !
#include "sprng_f.h"


  integer, parameter :: n_rayons_start = 100 
  integer, parameter :: n_rayons_start2 = 100
  integer :: n_rayons_max
  integer :: etape, etape_start, etape_end, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status, la
  integer, dimension(nb_proc) :: max_n_iter_loc
  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lconverged_loc, lprevious_converged
  real :: rand, rand2, rand3
  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt
  real(kind=dp) :: diff, norme, dT_max, dN, dN1, dJ, lambda_max
  real(kind=dp), allocatable :: dM(:), Jold(:,:), taub(:,:,:)
                                   !futur global flag
  logical :: labs, disable_subit, iterate_ne = .false., accelerated, ng_rest, l_unconverged
  integer :: nact, maxIter, imax
  integer :: icell, iorder, i0_rest, n_iter_accel
  integer :: Nlevel_total = 0, NmaxLevel, ilevel, max_sub_iter
  character(len=20) :: ne_start_sol = "NE_MODEL"
  real(kind=dp), dimension(3, atmos%Nrays, nb_proc) :: xyz0, uvw0
  type (AtomType), pointer :: atom


  write(*,*) "   -> Solving for kinetic equations for ", atmos%Nactiveatoms, " atoms"
  write(*,*) " Max error : ", dpops_max_error, dpops_sub_max_error
  
  if (atmos%include_xcoupling) CALL ERROR("MALI SCHEME NOT IMPLEMENTED YET")

  ds = 0.0_dp !meters
  chi_loc = 0.0_dp

  !move to initSol
  allocate(dM(atmos%Nactiveatoms)); dM=0d0 !keep tracks of dpops for all cells for each atom
  if (atmos%electron_scattering) then
   allocate(Jold(NLTEspec%Nwaves, atmos%Nspace))!move to initsol
   if (lread_jnu_atom) CALL read_Jnu() !initial solution from previous run or fix it ??
  endif


  n_rayons_max = atmos%Nrays
  xyz0(:,:,:) = 0d0
  uvw0(:,:,:) = 0d0

  labs = .true. !to have ds at cell icell = eval_operator
  id = 1
  etape_start = 2
  etape_end = 2
  if (etape_start==1) then
   write(*,*) "Warn: etape 1 not accurate"
  endif

  disable_subit = .false.

  accelerated  = .false.
  if (lNg_acceleration) then
   CALL error("+Ng no tested yet")
   n_iter_accel = 0
   i0_rest = 0
   ng_rest = .false.
  endif

  iterate_ne = (n_iterate_ne>0)
  if (iterate_ne) then
   write(*,*) " before iterate ne I need te recompute gij for continua !!"
  endif

  if (lforce_lte) disable_subit = .true.

  max_sub_iter = 1000
  maxIter = 1000
 ! ----------------------------  INITIAL POPS------------------------------------------ !
   CALL Init_NLTE(sub_iterations_enabled=.not.disable_subit)
 !  ------------------------------------------------------------------------------------ !

     do etape=etape_start, etape_end

      if (etape==1) then !two rays, not accurate
        lfixed_rays=.true.
        n_rayons = 2
        iray_start = 1
        lprevious_converged = .false.

      else if (etape==2) then 
  		lfixed_rays = .true.
  		n_rayons = min(n_rayons_max,n_rayons_start2)
  		iray_start = 1
  		lprevious_converged = .false.

  	  else
  	    CALL ERROR("etape unkown")
  	  end if
  	  
  	  write(*,*)  "-> Using ", n_rayons, ' rays for step', etape

  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0

        do while (.not.lconverged)

        	n_iter = n_iter + 1
        	if (n_iter > maxIter) exit !change step
            write(*,*) " -> Iteration #", n_iter, " Step #", etape

  			if (lfixed_rays) then
    			stream = 0.0
    			do i=1,nb_proc
     				stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    			end do
    			
    			!compute Profile(:,iray,id) for each line ??, it doesn't change for sub-iter
 			 end if

            max_n_iter_loc = 0
            dT_max = 0.
            

 			!$omp parallel &
            !$omp default(none) &
            !$omp private(icell, id, atom,nact) &
            !$omp shared(gpop_old,atmos)
            !$omp do schedule(static,1)
            do icell=1, atmos%Nspace
   				!$ id = omp_get_thread_num() + 1
   				if (atmos%icompute_atomRT(icell)>0) then
            	do nact=1, atmos%NactiveAtoms
            	    atom => atmos%ActiveAtoms(nact)%ptr_atom
             		gpop_old(nact, 1:atom%Nlevel,icell) = atom%n(:,icell)
             		atom => NULL()
            	end do
            	end if
            end do
        	!$omp end do
        	!$omp end parallel
        	
        	!0 or another start value
        	if (atmos%electron_scattering) Jold = NLTEspec%Jc

 			!$omp parallel &
            !$omp default(none) &
            !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02, la, dM, dN, dN1,dT_max)&
            !$omp private(argmt,n_iter_loc,lconverged_loc,diff,norme, icell, nact, atom, l_unconverged) &
            !$omp shared(atmos,NLTEspec, taub,dpops_sub_max_error) &
            !$omp shared(xyz0, uvw0, lkeplerian,n_iter) & !before nact was shared
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,max_sub_iter,lcell_converged) &
            !$omp shared(n_cells, pop_old, disable_subit, gpop_old,lforce_lte) & !dN, dN1,!pop
            !$omp shared(lfixed_Rays,lnotfixed_Rays,labs,max_n_iter_loc, etape)
            !$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !$ id = omp_get_thread_num() + 1
   				!l_unconverged = (atmos%icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
   				!if (l_unconverged) then
   				if (atmos%icompute_atomRT(icell)>0) then

					CALL fill_Collision_matrix(id, icell)
	                CALL NLTE_bound_free(id, icell, .true.) !also fill etac for the (this )iterate cell and set to 0 for this cell
					!for all atoms, set atom%Gamma to atom%C. atom%C is not updated until
					!next cell point. Gamma is initialized at each sub-iterations if any
  			        CALL initGamma(id)
                    
                    CALL init_psi_operator_new(id)
                    
					do iray=iray_start, iray_start-1+n_rayons
    					if (etape==1) then
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

      	 				else !etape 2
                   	    ! Position aleatoire dans la cellule
         					rand  = sprng(stream(id))
            				rand2 = sprng(stream(id))
            				rand3 = sprng(stream(id))

            				CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
            				rand = sprng(stream(id))
           					 W0 = 2.0_dp * rand - 1.0_dp !nz
            					W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
            				SRW02 = sqrt(W02)
            				rand = sprng(stream(id))
            				ARGMT = PI * (2.0_dp * rand - 1.0_dp)
            				U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
            				V0 = SRW02 * sin(ARGMT) !ny = sin(theta) * sin(phi)
		 				end if !etape

						!keep them for subiterations or keep phi(:,iray,id)
		 				xyz0(1,iray,id) = x0; xyz0(2,iray,id) = y0; xyz0(3,iray,id) = z0
		 				uvw0(1,iray,id) = U0; uvw0(2,iray,id) = V0; uvw0(3,iray,id) = W0					
					
					
					
						CALL INTEG_RAY_LINE(id, icell, xyz0(1,iray,id),xyz0(2,iray,id), xyz0(3,iray,id), &
											uvw0(1,iray,id), uvw0(2,iray,id), uvw0(3,iray,id), iray, labs)
											
                        !CALL FillGamma_mu(id, icell, iray, n_rayons, lforce_lte) !ray by ray since wphi is known
                        CALL FillGamma_part1(id, icell, iray, n_rayons, lforce_lte) !lambda+ray integ line, compute Jnu cont
      			    enddo !iray
      			    CALL FillGamma_part2(id, icell, lforce_lte)
      				if (atmos%electron_scattering) CALL calc_Jnu(id, icell, n_rayons)

    		     	n_iter_loc = 0
    				if (disable_subit) then
    				 !if lforce_lte, FillGamma_mu is returned after entered
    				 !and atom%Gamma = atom%C. delta term removed in SEE
    				 CALL updatePopulations(id, icell, dN,.false.)
    				 dT_max = max(dT_max, dN)
    				 lconverged_loc = .true.
    				else
    				 lconverged_loc = .false.
    				end if


     				!!only iterate on cells which are not converged
     				do while (.not.lconverged_loc)
       					n_iter_loc = n_iter_loc + 1
                        
                        !obtain for all lines and continua of atoms, the line%Jbar and cont%jnu
                        !with Hogerheijde method
                        do nact=1, atmos%Nactiveatoms
                        	atom => atmos%Activeatoms(nact)%ptr_atom
                        	do iray=iray_start, iray_start-1+n_rayons
                         		CALL calc_J_atom(id, icell, iray, atom, n_rayons)
                        	enddo
                        	atom => NULL()
                        enddo

                        do nact=1, atmos%Nactiveatoms
                         atom => atmos%Activeatoms(nact)%ptr_atom
                         CALL rate_matrix_atom(id, icell, atom, lforce_lte)
                         atom => NULL()
                        enddo
                        
						CALL updatePopulations(id, icell, diff,.false.)

                        !for this icell (id) and local iteration
       					if (diff < dpops_sub_max_error) then!precision_sub
       						lconverged_loc = .true.
     					    !write(*,*) "converged", id, icell, "dpops(sub) = ", diff
       					else
       						!for this cell (or these threads) set atom%Gamma to atom%C, and add Aji to Gamma(j,i)
       						!for each atom. atom%C is computed for this cell at the start of the loop
  			        		CALL initGamma(id)
  			        		CALL init_psi_operator_new(id)
  			        		!re-update, and compute etac(:,id), Vji(:,id) to be use in FillGamma_mu
  			        		CALL NLTE_bound_free(id, icell,.true.)
  			        		
  			        		
 							do iray=iray_start, iray_start-1+n_rayons
      							!I unchanged
      							
      							
      							!does not recompute profile phi(:,iray,id) and chi_loc(:,iray,id)
      							!wich is computed in integ_ray_line() above
      						    CALL init_eta_atom_loc(id,icell,iray)
      						    !iterate always true in there

      						enddo !iray

       					end if
       					if (n_iter_loc >= max_sub_iter) then 

       					  lconverged_loc = .true.

       					end if

     				end do !local sub iteration
     	            if (n_iter_loc > max_n_iter_loc(id)) max_n_iter_loc(id) = n_iter_loc
     	            
     			end if !icompute_atomRT
     		end do !icell
        	!$omp end do
        	!$omp end parallel
        	
     		

            !should be para
            dM(:) = 0d0
        	diff = 0d0
        	dJ = 0.0
  			cell_loop2 : do icell=1, atmos%Nspace
  				if (atmos%icompute_atomRT(icell)>0) then
						dN = 0d0 !for all levels of all atoms of this cell
     					do nact=1,atmos%NactiveAtoms
     					    atom => atmos%ActiveAtoms(nact)%ptr_atom
     					    
     					    !CALL calc_delta_Tex_atom(icell, atom, dN1)
     					    !dN = max(dN, dN1) !all atoms for this cell
     					    
     						do ilevel=1,atom%Nlevel
     						
								if (atom%n(ilevel,icell) >= prec_pops) then
     				    		 dN1 = dabs(1d0-gpop_old(nact,ilevel,icell)/atom%n(ilevel,icell))
     				    		 !dabs(atom%n(l,icell) - ndag(l))/atom%n(l,icell)
     				    		 dN = max(dN1, dN) !compare with all atoms and levels
     				    						  !for this cell
     				    		 dM(nact) = max(dM(nact), dN1) !compare for one atom
     				    		 							   !for all cells
     				    		endif
     						end do !over ilevel
     						atom => NULL()
     					end do !over atoms
     					diff = max(diff, dN) !compare for all atoms and all cells
     					
     					if (atmos%electron_scattering) then !local dJ
     					    dN1 = dabs(1d0 - maxval(Jold(:,icell)/(tiny_dp + NLTEspec%Jc(:,icell))))
     					    if (dN1 > dJ) then
     						   dJ = dN1
     						   imax = locate(dabs(1d0 - Jold(:,icell)/(tiny_dp + NLTEspec%Jc(:,icell))), dJ)
     						   lambda_max = NLTEspec%lambda(imax)
     						 endif
     					     lcell_converged(icell) = (real(dJ) < dpops_max_error)

     					else
     					     lcell_converged(icell) = (real(diff) < dpops_max_error)
     					     !lcell_converged(icell) = (real(dT_max) < dpops_max_error)								
     					endif
     					
     			end if
     		end do cell_loop2 !icell
     		

         	if (.not.disable_subit)	write(*,*) maxval(max_n_iter_loc), "sub-iterations"
         	if (atmos%electron_scattering) then 
         	 write(*,*) " -> dJ = ", dJ, " @ ", lambda_max, " nm"
         	endif
         	do nact=1,atmos%NactiveAtoms
         	 if (accelerated) then
         	  write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dM = ", dM(nact)," (Accelerated)"
         	 else
         	  write(*,*) "   >>> ", atmos%ActiveAtoms(nact)%ptr_atom%ID, " dM = ", dM(nact)
         	 endif
         	enddo
         	write(*,*) "dpops =", diff !at the end of the loop over n_cells

			if (atmos%electron_scattering) diff = dJ

        	lconverged = (real(diff) < dpops_max_error) !global convergence for all iterations

        	if (real(diff) < dpops_max_error) then !precision
           		if (lprevious_converged) then
            	  lconverged = .true.
           		else
            	  lprevious_converged = .true.
          	    endif
        	else
           		lprevious_converged = .false.
           		if (.not.lfixed_rays) then
              		n_rayons = n_rayons * 2
              		write(*,*) ' -- Increasing number of rays'
             		if (n_rayons > n_rayons_max) then
              			if (n_iter >= maxIter) then
             		 		write(*,*) "Warning : not enough rays to converge !!"
                 			lconverged = .true.
              			end if
              		end if

              ! On continue en calculant 2 fois plus de rayons
              ! On les ajoute a l'ensemble de ceux calcules precedemment
!              iray_start = iray_start + n_rayons

          	   end if
        	end if
        									
        	!Only if specified
        	!if (disable_subit) then
        		if (iterate_ne .and. (mod(n_iter,n_iterate_ne)==0))  then
        			write(*,*) " Solve ne Global iteration:", n_iter
        	 		write(*,*) " --> old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
        	 		CALL SolveElectronDensity(ne_start_sol)
        	 !Recompute LTE pops used in continua radiative rates
        	 !for Activeatoms only ?
        	 		do nact=1,atmos%NactiveAtoms !if over Active only, should recompute for passive also
        	 							  !but this preserve the constant background opac
        	 							  !as passiveatoms%n = passiveatoms%nstar with ne init.
               			atom => atmos%ActiveAtoms(nact)%ptr_atom
               			CALL LTEpops(atom,.true.)
               			atom => Null()
             		end do
        		end if
        	!endif

	    end do !while
        write(*,*) etape, "Threshold =", dpops_max_error
	  end do !over etapes

  !LTE opacities NOT updated if not ltab_Wavelength_grid
  !Force to compute a new value of electron density after convergence
  !and new lte populations for all atoms
  if (n_iterate_ne < 0) then
   write(*,*) "END LOOP: old max/min ne", maxval(atmos%ne), minval(atmos%ne,mask=atmos%ne>0)
   CALL SolveElectronDensity(ne_start_sol)
   !Recompute for all atoms the LTE pops
   write(*,*) "   recompute LTE pops for all atoms"
   do nact=1,atmos%NAtom
      atom => atmos%Atoms(nact)%ptr_atom
      CALL LTEpops(atom,.true.)
      atom => Null()
   end do
  else if (iterate_ne) then
   write(*,*) "END LOOP: recompute LTE pops of passive atoms"
   !Recompute for all passive atoms the LTE pops
   !because the lte pops of active atoms only is updated with iterate_ne
   do nact=1,atmos%NpassiveAtoms
      atom => atmos%Atoms(nact)%ptr_atom
      CALL LTEpops(atom,.true.)
      atom => Null()
   end do
  end if

  if (iterate_ne .or. (n_iterate_ne < 0)) CALL writeElectron(.true.) !the .true. means append, to compare with initial solution.
 ! -------------------------------- CLEANING ------------------------------------------ !
  !to move inside free_nlte_sol
  deallocate(dM)
  if (allocated(Jold)) deallocate(Jold)
  CALL free_nlte_sol(disable_subit)
  !!if (atmos%include_xcoupling) deallocate(taub)

 ! ------------------------------------------------------------------------------------ !
 RETURN
 END SUBROUTINE NLTEloop
 ! ------------------------------------------------------------------------------------ !

 SUBROUTINE adjustStokesMode(Solution)
 !
 !
   character(len=*), intent(in) :: Solution


write(*,*) "Note polarized profile not allocated yet, check that in profile and line%phiZ, line%psi"
write(*,*) "Magnetic profile not ready yet"
stop

    if (SOLUTION=="FIELD_FREE") then
      !if (associated(Profile)) Profile => NULL()
      !Profile => Iprofile !already init to that
      if (allocated(NLTEspec%AtomOpac%rho_p)) deallocate(NLTEspec%AtomOpac%rho_p)
      if (allocated(NLTEspec%AtomOpac%etaQUV_p)) deallocate(NLTEspec%AtomOpac%etaQUV_p)
      if (allocated(NLTEspec%AtomOpac%chiQUV_p)) deallocate(NLTEspec%AtomOpac%chiQUV_p)
      if (allocated(NLTEspec%F_QUV)) deallocate(NLTEspec%F_QUV)
      if (allocated(NLTEspec%StokesV)) deallocate(NLTEspec%StokesV, NLTEspec%StokesQ, &
      												NLTEspec%StokesU)
      write(*,*) " Using FIELD_FREE solution for the SEE!"
      CALL Warning("  Set PRT_SOLUTION to FULL_STOKES for images")

    else if (SOLUTION=="POLARIZATION_FREE") then
     CALL ERROR("POLARIZATION_FREE solution not implemented yet")
     !NLTE not implemented in integ_ray_line_z and Zprofile always compute the full
     !propagation dispersion matrix, but only the first row is needed
     if (associated(Profile)) Profile => NULL()
     Profile => Zprofile
     if (associated(INTEG_RAY_LINE)) INTEG_RAY_LINE => NULL()
     INTEG_RAY_LINE => INTEG_RAY_LINE_Z

     !if (associated(Voigt)) Voigt => NULL
     !Voigt => VoigtHumlicek !Already the algortihm used

     !beware, only I! polarization is neglected hence not allocated
      if (.not.allocated(NLTEspec%AtomOpac%rho_p)) then !same for all, dangerous.
      	allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%etaQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%chiQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))

        write(*,*) " Using POLARIZATION_FREE solution for the SEE!"
     end if

    else if (SOLUTION=="FULL_STOKES") then !Solution unchanged here
      if (associated(Profile)) Profile => NULL()
      Profile => Zprofile
      if (associated(INTEG_RAY_LINE)) INTEG_RAY_LINE => NULL()
      INTEG_RAY_LINE => INTEG_RAY_LINE_Z
     !if (associated(Voigt)) Voigt => NULL
     !Voigt => VoigtHumlicek !Already the algortihm used
      !allocate space for Zeeman polarisation
      !only LTE for now
       if (.not.allocated(NLTEspec%StokesQ)) & !all allocated, but dangerous ?
            allocate(NLTEspec%StokesQ(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC), &
             NLTEspec%StokesU(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC), &
             NLTEspec%StokesV(NLTEspec%NWAVES, atmos%NRAYS,NLTEspec%NPROC))!, &
             !!NLTEspec%S_QUV(3,NLTEspec%Nwaves))
      NLTEspec%StokesQ=0d0
      NLTEspec%StokesU=0d0
      NLTEspec%StokesV=0d0
      if (.not.allocated(NLTEspec%F_QUV)) &
      	allocate(NLTEspec%F_QUV(3,NLTEspec%Nwaves,NPIX_X,NPIX_Y,RT_N_INCL,RT_N_AZ))
      NLTEspec%F_QUV = 0d0
      if (.not.allocated(NLTEspec%AtomOpac%rho_p)) then !same for all, dangerous.
      	allocate(NLTEspec%AtomOpac%rho_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%etaQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))
        allocate(NLTEspec%AtomOpac%chiQUV_p(NLTEspec%Nwaves, 3, NLTEspec%NPROC))

        write(*,*) " Using FULL_STOKES solution for the SEE!"
      end if
    else
     CALL ERROR("Error in adjust StokesMode with solution=",Solution)
    end if
 RETURN
 END SUBROUTINE adjustStokesMode
 ! ------------------------------------------------------------------------------------ !
!should be done more properly and I should use functions from mcfostfost for stellar flux
 SUBROUTINE reallocate_mcfost_vars()
  !--> should move them to init_atomic_atmos ? or elsewhere
  !need to be deallocated at the end of molecule RT or its okey ?`
  integer :: la
  CALL init_directions_ray_tracing()
  if (.not.allocated(stars_map)) allocate(stars_map(npix_x,npix_y,3))
  if (.not.allocated(stars_map_cont)) allocate(stars_map_cont(npix_x,npix_y,3))
  n_lambda = NLTEspec%Nwaves
  if (allocated(tab_lambda)) deallocate(tab_lambda)
  allocate(tab_lambda(n_lambda))
  if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
  allocate(tab_delta_lambda(n_lambda))
  tab_lambda = NLTEspec%lambda / MICRON_TO_NM
  tab_delta_lambda(1) = 0d0! + tiny_dp !! I assume that Nwaves > 1 here
  !!if (NLTEspec%Nwaves==1) tab_delta_lambda(1) = tab_delta_lambda(1) + tiny_dp

  do la=2,NLTEspec%Nwaves
   tab_delta_lambda(la) = tab_lambda(la) - tab_lambda(la-1)
  end do

  if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
  allocate(tab_lambda_inf(n_lambda))
  if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
  allocate(tab_lambda_sup(n_lambda))
  tab_lambda_inf = tab_lambda
  tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
  ! computes stellar flux at the new wavelength points
  CALL deallocate_stellar_spectra()
  ! probably do not deallocate, 'cause maybe dust is in here too.
  if (allocated(kappa)) deallocate(kappa) !do not use it anymore
  !used for star map ray-tracing.
  CALL allocate_stellar_spectra(n_lambda)
  CALL repartition_energie_etoiles()
  !If Vfield is used in atomic line transfer, with lkep or linfall
  !atmos%Vxyz is not used and lmagnetoaccr is supposed to be zero
  !and Vfield allocated and filled in the model definition if not previously allocated
  !nor filled.
  if (allocated(Vfield).and.lkeplerian.or.linfall) then
    write(*,*) "Vfield max/min", maxval(Vfield), minval(Vfield)
  else if (allocated(Vfield).and.lmagnetoaccr) then
   write(*,*) "deallocating Vfield for atomic lines when lmagnetoaccr = .true."
   deallocate(Vfield)
  end if

 RETURN
 END SUBROUTINE reallocate_mcfost_vars

 SUBROUTINE reallocate_mcfost_wavelength_arrays()
  !--> should move them to init_atomic_atmos ? or elsewhere
  !need to be deallocated at the end of molecule RT or its okey ?`
  integer :: la
  n_lambda = NLTEspec%Nwaves
  if (allocated(tab_lambda)) deallocate(tab_lambda)
  allocate(tab_lambda(n_lambda))
  if (allocated(tab_delta_lambda)) deallocate(tab_delta_lambda)
  allocate(tab_delta_lambda(n_lambda))
  tab_lambda = NLTEspec%lambda / MICRON_TO_NM
  tab_delta_lambda(1) = 0d0

  do la=2,NLTEspec%Nwaves
   tab_delta_lambda(la) = tab_lambda(la) - tab_lambda(la-1)
  end do
  !unlikely in NLTE because of the different atomic grids
  !but if an image at only one wavelength --> tab_Delta_lambda=0
  !!if I do not put a tiny_dp here, if only one wavelength delta_wl is 0
  !which creates an overflow in stars.f90, because spectre and spectre0 are 0
  if (size(tab_lambda)==1) tab_delta_lambda = tab_delta_lambda + tiny_dp

  if (allocated(tab_lambda_inf)) deallocate(tab_lambda_inf)
  allocate(tab_lambda_inf(n_lambda))
  if (allocated(tab_lambda_sup)) deallocate(tab_lambda_sup)
  allocate(tab_lambda_sup(n_lambda))
  tab_lambda_inf = tab_lambda
  tab_lambda_sup = tab_lambda_inf + tab_delta_lambda
  ! computes stellar flux at the new wavelength points
  CALL deallocate_stellar_spectra()
  CALL allocate_stellar_spectra(n_lambda)
  CALL repartition_energie_etoiles()

 RETURN
 END SUBROUTINE reallocate_mcfost_wavelength_arrays
 
 !futur move to initial_solution
 SUBROUTINE init_stellar_disk
  integer :: i_star!, lam
  
   write(*,*) " Computing Istar(mu=1) for each star..."
   !move elsewhere if we read a file, but normally should well fit with MCFOST routines
   do i_star=1, n_etoiles
    CALL Bplanck(etoile(i_star)%T*1d0, NLTEspec%Istar(:,i_star))
   enddo
   write(*,*) " ..done"
 
 RETURN
 END SUBROUTINE init_stellar_disk
 
!  futur move to stars.f90 and merging with mcfost generator of stellar radiation
!  SUBROUTINE calc_stellar_radiation(N,i_star,x,y,z,u,v,w,gamma)
!  ---------------------------------------------------------------!
!   Compute the stellar radiation field and in Istar:
!   Radiation emerging from the star (at the surface), corrected
!   by limb darkening.
!   
!   tab_lambda has to be defined and repartition_energie_etoiles()
!  -------------------------------------------------------------- !
!   use input, only : limb_darkening, mu_limb_darkening
!   use Planck, only : uLD, Bplanck
!   integer, intent(in) :: N, i_star
!   real(kind=dp), dimension(N), intent(out) :: gamma
!   real(kind=dp), intent(in) :: u, v, w, x, y, z
!   real(kind=dp) :: energie(N)
!   real(kind=dp) :: mu, ulimb, LimbDarkening, surface, HC
!   integer :: ns
!   logical :: lintersect_spot
! 
!    gamma(:) = 1d0
!    cos(theta) = dot(r,n)/module(r)/module(n)
!    mu = abs(x*u + y*v + z*w)/dsqrt(x**2+y**2+z**2) !n=(u,v,w) is normalised
!    if (real(mu)>1d0) then !to avoid perecision error
!     write(*,*) "mu=",mu, x, y, z, u, v, w
!     CALL Error(" mu limb > 1!")
!    end if
!    
!    1) Compute stellar flux from mcfost 
! 
!    2) Correct with the contrast gamma of a hotter/cooler region if any
!    CALL intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect_spot)
!    if (lintersect_spot) then
!      gamma(:) = (dexp(hc_k/NLTEspec%lambda(:)/real(etoile(i_star)%T,kind=dp))-1)/&
!      			(dexp(hc_k/NLTEspec%lambda(:)/etoile(i_star)%SurfB(ns)%T)-1)
!      so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
!    end if
! 
!    3) Apply Limb darkening
!    if (llimb_darkening) then
!      CALL ERROR("option for reading limb darkening not implemented")
!    else
!      write(*,*) maxval(uLD(real(etoile(i_star)%T,kind=dp))), minval(uLD(real(etoile(i_star)%T,kind=dp)))
!      ulimb = 0.0 ! could use BB slope
!      LimbDarkening = 1d0 - ulimb*(1d0-mu)
!    end if
!    Istar(:) = energie(:) * LimbDarkening * gamma(:)
!    gamma(:) = LimbDarkening * gamma(:)
! 
!  RETURN
!  END SUBROUTINE calc_stellar_radiation
 
 SUBROUTINE calc_stellar_surface_brightness(N,lambda,i_star,x,y,z,u,v,w,gamma)
 ! ---------------------------------------------------------------!
  ! Compute the stellar radiation field surface brightness.
  ! Istar = B(x,y,z;mu) * Stellar spectrum or B * BlackBody
  ! return gamma, the brightness. For uniform disk, gamma is 1.
  !
  ! For a point on the stellar disk, returns gamma * limbdarkenig
 ! -------------------------------------------------------------- !
  use input, only : limb_darkening, mu_limb_darkening
  use Planck, only : uLD, Bplanck
  integer, intent(in) :: N, i_star
  real(kind=dp), dimension(N), intent(in) :: lambda
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
   
   !1) Compute stellar flux from mcfost 

   !2) Correct with the contrast gamma of a hotter/cooler region if any
   CALL intersect_spots(i_star,u,v,w,x,y,z, ns,lintersect_spot)
   if (lintersect_spot) then
     gamma(:) = (dexp(hc_k/lambda/real(etoile(i_star)%T,kind=dp))-1)/&
     			(dexp(hc_k/lambda/etoile(i_star)%SurfB(ns)%T)-1)
     !so that I_spot = Bplanck(Tspot) = Bp(Tstar) * gamma = Bp(Tstar)*B_spot/B_star
   end if

   !3) Apply Limb darkening
   if (llimb_darkening) then
     CALL ERROR("option for reading limb darkening not implemented")
   else
     !write(*,*) maxval(uLD(real(etoile(i_star)%T,kind=dp))), minval(uLD(real(etoile(i_star)%T,kind=dp)))
     ulimb = 0.0 ! could use BB slope
     LimbDarkening = 1d0 - ulimb*(1d0-mu)
   end if
   !Istar(:) = energie(:) * LimbDarkening * gamma(:)
   gamma(:) = LimbDarkening * gamma(:)

 RETURN
 END SUBROUTINE calc_stellar_surface_brightness

   SUBROUTINE INTEG_RAY_JNU(id,icell_in,x,y,z,u,v,w,iray,labs, kappa_tot, Snu, Istar, psi, Ic)
 ! ------------------------------------------------------------------------------- !
  ! This routine performs integration of the transfer equation along a ray to
  ! compute coherent Jnu in the continuum
 ! ------------------------------------------------------------------------------- !

  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  real(kind=dp), intent(in), dimension(:,:) :: kappa_tot, Snu
  real(kind=dp), intent(in), dimension(:) :: Istar
  real(kind=dp), intent(out) :: Ic(:,:), psi(:,:)
  logical, intent(in) :: labs
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
!   real(kind=dp) :: chil(size(Ic(:,1)), etal(size(Ic(:,1))
  real(kind=dp), dimension(size(Ic(:,1))) :: LimbD
  real(kind=dp), dimension(size(Ic(:,1))) :: tau_c!, tau
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star, la
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator

  x1=x;y1=y;z1=z
  x0=x;y0=y;z0=z
  next_cell = icell_in
  nbr_cell = 0

  !tau(:) = 0.0_dp !go from surface down to the star
  tau_c(:) = 0.0_dp

  Ic(:,id) = 0.0
  psi(:,id) = 0.0
  
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
    
    !if (minval(tau_c) > 50.) return
    
    ! Test sortie ! "The ray has reach the end of the grid"
    if (test_exit_grid(icell, x0, y0, z0)) RETURN

    if (lintersect_stars) then
      if (icell == icell_star) then
       !call calc_stellar_surface_brightness(size(Ic(:,1)),lambda,i_star,x0,y0,z0,u,v,w,LimbD)
       Ic(:,id) =  Ic(:,id) + Istar(:) * dexp(-tau_c)
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

      Ic(:,id) = Ic(:,id) + dexp(-tau_c) * (1.0_dp - dexp(-l_contrib * kappa_tot(:,icell))) * Snu(:,icell)

     
     if ((nbr_cell == 1).and.labs) then 
      ds(iray,id) = l * AU_to_m
      psi(:,id) = (1d0 - exp(-ds(iray,id)*kappa_tot(:,icell)))
     endif


     !tau = tau + dtau
     tau_c = tau_c + l_contrib * kappa_tot(:,icell)

    end if  ! lcellule_non_vide
  end do infinie
  ! -------------------------------------------------------------- !
  ! -------------------------------------------------------------- !
  RETURN
  END SUBROUTINE INTEG_RAY_JNU
  
 SUBROUTINE Iterate_Jnu()
 ! -------------------------------------------------------- !
  ! Compute the mean radiation field at all cells
  ! evaluated on a small grid and then interpolated on the
  ! wavelength grid
 ! -------------------------------------------------------- !
#include "sprng_f.h"

  integer, parameter :: n_rayons_start = 200, Nlambda = 350, maxIter = 1000
  integer :: n_rayons_max
  real, parameter :: precision = 1e-2
  real, parameter :: lambda_min = 5., lambda_max0 = 100000
  integer :: etape, etape_start, etape_end, iray, n_rayons
  integer :: n_iter, n_iter_loc, id, i, iray_start, alloc_status
  logical :: lfixed_Rays, lnotfixed_Rays, lconverged, lprevious_converged, write_convergence_file

  real :: rand, rand2, rand3, a1, a0, a2, a3!, fac_etape

  real(kind=dp) :: x0, y0, z0, u0, v0, w0, w02, srw02, argmt, diff, norme, dN, dJ, diffs, dSo
  real(kind=dp), allocatable :: Jold(:,:), Jflat(:), lambda(:), Jnew(:,:), lambda_1(:,:), Jnew_l(:,:), Jold_l(:,:)
  real(kind=dp), allocatable :: Snew(:,:), Kappa_tot(:,:), beta(:,:), Istar(:), Ic(:,:)
  real(kind=dp), allocatable :: lambda_star(:,:), Sth(:,:), Sold(:,:), Sline(:,:)
                                   !futur global flag
  logical :: labs, l_unconverged, accelerated, ng_rest, ng_acc
  integer :: la, icell, imax, icell_max, iorder,i0_rest, icell_max_s, imax_s
  real(kind=dp) :: lambda_max
  type(Ng) :: NgJ
  
  write(*,*) "   --> Lambda iterating Jnu with Nlambda ", Nlambda
  write_convergence_file = .false.
  
  if (allocated(ds)) deallocate(ds)
  allocate(ds(atmos%Nrays, NLTEspec%NPROC))
  ds = 0.0_dp !meters
  
  !this would allow to compute easily Jnu with lines or I need a version
  !of Profile which use as input the wavelength grid
  !call initSpectrum_jnu(Nlambda, lambda_min, lambda_max0)

  !only one star
  allocate(Istar(Nlambda), Ic(Nlambda, NLTEspec%NPROC), lambda_star(Nlambda, NLTEspec%NPROC), lambda_1(Nlambda, NLTEspec%NPROC))
  Ic = 0.0_dp; Istar = 0.0_dp; lambda_star = 0.0_dp; lambda_1 = 0.0

  allocate(Jold(Nlambda, atmos%Nspace), Jnew(Nlambda, atmos%Nspace))
  Jold = 0d0; Jnew(:,:) = 0.0_dp
  allocate(Sth(Nlambda, atmos%Nspace), Snew(Nlambda, atmos%Nspace), Sold(Nlambda, atmos%Nspace))
  Sth = 0.; Sold = 0.0; Snew = 0.0
  allocate(Kappa_tot(Nlambda, atmos%Nspace)); Kappa_tot = 0.
  allocate(beta(Nlambda, atmos%Nspace)); beta = 0.
    
  !sampling lambda for Jnu, at the end interpolated on the NLTE grid
  allocate(lambda(Nlambda))
  a1 = real(NLTEspec%lambda(1)) !max(lambda_min, real(NLTEspec%lambda(1)))
  a2 = min(lambda_max0, real(NLTEspec%lambda(NLTEspec%Nwaves)))
  a0 = 368.
  a3 = 91.
  !evaluate opacities on the exact same grid, better than interpolate
  lambda(1:20) = span(a1, a3-0.1, 20)
  lambda(21:50) = span(a3, a3+2, 30)
  lambda(51:51+249) = span(a3+2+0.1, a0, 250)
  lambda(301:350) = spanl(a0+10, a2, 50)
!   lambda(1:20) = span(a1, a3-0.1, 20)
!   lambda(21:20+150) = span(a3, a3+2, 150)
!   lambda(171:170+300) = span(a3+2+0.1, a0, 300)
!   lambda(471:470+530) = spanl(a0+10, a2, 530)

  
  Istar(:) = Bpnu (real(etoile(1)%T,kind=dp), lambda)

  write(*,*) "  -> interpolating contopac on Jnu grid for each frequency.."
  write(*,*) "       lambda (min, max in nm):", minval(lambda), maxval(lambda)
  do icell=1, atmos%Nspace
   if (atmos%icompute_atomRT(icell) /= 0.) then
    
    !LTE at start ?
    !Jold(:,icell) = bpnu(atmos%T(icell), lambda)
    
   					
   	CALL bezier2_interp(NLTEspec%Nwaves,NLTEspec%lambda,&
   					NLTEspec%AtomOpac%Kc(:,icell),Nlambda, lambda, kappa_tot(:,icell))
   					
   					
   	CALL bezier2_interp(NLTEspec%Nwaves,NLTEspec%lambda,&
   					NLTEspec%AtomOpac%sca_c(:,icell),Nlambda, lambda, beta(:,icell))
   	beta(:,icell) = beta(:,icell) / kappa_tot(:,icell)
   					

   	CALL bezier2_interp(NLTEspec%Nwaves,NLTEspec%lambda,&
   					NLTEspec%AtomOpac%jc(:,icell),Nlambda, lambda, Sth(:,icell))
   	Sth(:,icell) = Sth(:,icell) / ( kappa_tot(:,icell) * (1.-beta(:,icell)))
   	
   					
    Sold(:,icell) = (1.-beta(:,icell))*Sth(:,icell) + beta(:,icell) * Jold(:,icell) 

   endif
  enddo
  write(*,*) "Sold (max,min)", maxval(Sold), minval(Sold)
  write(*,*) "Sth (max,min)", maxval(STh), minval(Sth)
  write(*,*) "Jold (max,min)", maxval(Jold), minval(Jold)
  write(*,*) "beta (max,min)", maxval(beta), minval(beta)
  write(*,*) "kappa (max,min)", maxval(kappa_tot), minval(kappa_tot)
  write(*,*) "  ->..done"


  if (.not.allocated(lcell_converged)) allocate(lcell_converged(n_cells))
  lcell_converged(:) = .false.

  n_rayons_max = atmos%Nrays
  
  ng_acc = .false. !temp, because don't know if it works
  accelerated  = .false.
  if (ng_acc) then
   allocate(Jflat(atmos%Nspace*Nlambda))
   CALL initNg(atmos%Nspace*Nlambda, 6, 2, 3, NgJ)
  endif


  labs = .true.
  id = 1
  etape_start = 2
  etape_end = 2

if (write_convergence_file ) then
 open(unit=20, file="Jnu_convergence.s", status="unknown")
 write(20,*) maxIter, Nlambda, n_cells
endif
 
     do etape=etape_start, etape_end

      if (etape==1) then 
        lfixed_Rays = .true.
        n_rayons = 2
        iray_start=1
        lprevious_converged = .false.
		write(*,*) " Using ", n_rayons, " rays for Jnu."
      else if (etape==2) then 
  		lfixed_rays = .true.
  		n_rayons = min(n_rayons_max,n_rayons_start)
  		iray_start = 1
  		lprevious_converged = .false.
		write(*,*) " Using ", n_rayons, " rays for Jnu."
  	  else
  	    CALL ERROR("etape unkown")
  	  end if
  	  
if (write_convergence_file ) write(20,*) etape, n_rayons

  	  
  		lnotfixed_rays = .not.lfixed_rays
  		lconverged = .false.
  		n_iter = 0
		 				

        do while (.not.lconverged)

        	n_iter = n_iter + 1
        	if (n_iter > maxIter)  exit
            write(*,*) " -> Iteration #", n_iter, " Step #", etape

  			if (lfixed_rays) then
    			stream = 0.0
    			do i=1,nb_proc
     				stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    			end do
 			 end if

if (write_convergence_file ) write(20,*) " -> Iteration #", n_iter
            imax = 1
            imax_s = 1
            
 			!$omp parallel &
            !$omp default(none) &
            !$omp private(id,iray,rand,rand2,rand3,x0,y0,z0,u0,v0,w0,w02,srw02)&
            !$omp private(argmt,norme, icell, l_unconverged) &
            !$omp shared(atmos,NLTEspec,lambda_star, Snew, Sold, Sth, lambda_1, Istar) &
            !$omp shared(lkeplerian,n_iter) &
            !$omp shared(stream,n_rayons,iray_start, r_grid, z_grid,lcell_converged) &
            !$omp shared(n_cells,ds, Jold, Jnew, beta, kappa_tot, Ic) &
            !$omp shared(lfixed_Rays,lnotfixed_Rays,labs,etape)
            !$omp do schedule(static,1)
  			do icell=1, n_cells
   			    !$ id = omp_get_thread_num() + 1
   				!!l_unconverged = (atmos%icompute_atomRT(icell)>0).and.(.not.lcell_converged(icell))
   				!!if (l_unconverged) then 
      			l_unconverged = (.not.lcell_converged(icell))
   				if (atmos%icompute_atomRT(icell)>0) then
   				   Jnew(:,icell) = 0. !here otherwise if all cell converged -> Jnew = 0 and we never enter here
           		   Snew(:,icell) = 0. 
           		   lambda_1(:,id) = 0.        		   
           		   
					do iray=iray_start, iray_start-1+n_rayons
					
    					if (etape==1) then
        					! Position = milieu de la cellule
        					x0 = r_grid(icell)
        					y0 = 0.0_dp
        					z0 = z_grid(icell)
        					if (lkeplerian) then
        						write(*,*) "lkeplerian should be false here"
        						stop
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

      	 				else !etape 2
                   	    ! Position aleatoire dans la cellule
         					rand  = sprng(stream(id))
            				rand2 = sprng(stream(id))
            				rand3 = sprng(stream(id))

            				CALL  pos_em_cellule(icell ,rand,rand2,rand3,x0,y0,z0)

                        ! Direction de propagation aleatoire
            				rand = sprng(stream(id))
           					W0 = 2.0_dp * rand - 1.0_dp !nz
            				W02 =  1.0_dp - W0*W0 !1-mu**2 = sin(theta)**2
            				SRW02 = sqrt(W02)
            				rand = sprng(stream(id))
            				ARGMT = PI * (2.0_dp * rand - 1.0_dp)
            				U0 = SRW02 * cos(ARGMT) !nx = sin(theta)*cos(phi)
            				V0 = SRW02 * sin(ARGMT) !ny = sin(theta) * sin(phi)
		 				end if !etape

						CALL INTEG_RAY_JNU(id, icell, x0, y0, z0, u0, v0, w0, iray, labs, &
											kappa_tot, Sold, Istar, lambda_star, Ic)
						!LI
		                Jnew(:,icell) = Jnew(:,icell) + Ic(:,id)/n_rayons
! 		                lambda_1(:,id) = lambda_1(:,id) + lambda_star(:,id) / n_rayons
           
		                !hogerheijde
! 		                Jnew(:,icell) = Jnew(:,icell) + ( Ic(:,id)*dexp(-ds(iray,id)*kappa_tot(:,icell)) +&
!  		                 (1.-dexp(-ds(iray,id)*kappa_tot(:,icell)))*Sold(:,icell) )/n_rayons	
!  		                 					
      			    enddo !iray

      		   endif !icompute_AtomRT
      		   
        	   Snew(:,icell) = Sth(:,icell) * (1.-beta(:,icell)) + Jnew(:,icell) * beta(:,icell)
        	   !!ALI, but lambda is not exactly the same as for hogerheijde
!         	   Snew(:,icell) = (Snew(:,icell)- beta(:,icell)*lambda_1(:,id) * Sold(:,icell)) / (1.-beta(:,icell)*lambda_1(:,id))
        	   
     		end do !icell
        	!$omp end do
        	!$omp end parallel

            !should be para
        	diff = 0d0
  			cell_loop2 : do icell=1, atmos%Nspace
  				if (atmos%icompute_atomRT(icell)>0) then

						dJ = 0.0_dp
						do la=1, Nlambda
						 if (Jnew(la, icell) > 0) then 
						  dJ = max(dJ,dabs(1.-Jold(la,icell)/Jnew(la,icell)))
						  imax = locate(dabs(1.-Jold(:,icell)/Jnew(:,icell)),dabs(1.-Jold(la,icell)/Jnew(la,icell)))
						 endif 
						enddo
						if (mod(icell,10)==0) write(*,*) icell, " ::> dJ(icell)", real(dJ)!," Jmax/min:", maxval(Jnew(:,icell)), minval(Jnew(:,icell))
if (write_convergence_file ) write(20,*) icell, " ::> dJ(icell)", real(dJ), " Jmax/min:", maxval(Jnew(:,icell)), minval(Jnew(:,icell))
						if (dJ > diff) then
						  diff = dJ
						  icell_max = icell
						  !write(*,*) icell_max, " ::> dJ(icell_max)", real(diff), maxval(Jnew(:,icell_max)), minval(Jnew(:,icell_max))
						endif
     					lcell_converged(icell) = (real(dJ) < precision)		
     			end if
     			Jold(:,icell) = Jnew(:,icell)
     			Sold(:,icell) = Snew(:,icell)
     		end do cell_loop2 !icell
     		
     		! include cells only that have converge for all frequencies, or for each cell / frequencies?
!      		Jold(:,:) = Jnew(:,:)
!      		Sold(:,:) = Snew(:,:)

     		
         	if (accelerated) then
         	  write(*,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff," (Accelerated)"
if (write_convergence_file ) write(20,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff," (Accelerated)"
         	else
         	  write(*,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff
if (write_convergence_file ) write(20,*) "   >>> ", icell_max, lambda(imax)," dJ = ", diff
         	endif
         	write(*,*) "   >>> ", "Jmax/min (icell_max):", maxval(Jnew(:,icell_max)), minval(Jnew(:,icell_max))
         	write(*,*) "   >>> ", "Jmax/min (all):", maxval(Jnew(:,:)), minval(Jnew(:,:))
if (write_convergence_file ) write(20,*) "   >>> ", "Jmax/min (icell_max):", maxval(Jnew(:,icell_max)), minval(Jnew(:,icell_max))
if (write_convergence_file ) write(20,*) "   >>> ", "Jmax/min (all):", maxval(Jnew(:,:)), minval(Jnew(:,:))
	        write(*,*) " <-> # unconverged cells : ", size(pack(lcell_converged,mask=lcell_converged==.false.)), &
	          100.*real(size(pack(lcell_converged,mask=lcell_converged==.false.)))/real(n_cells), "%"
if (write_convergence_file ) write(20,*) " <-> # unconverged cells : ", size(pack(lcell_converged,mask=lcell_converged==.false.)), &
	          100.*real(size(pack(lcell_converged,mask=lcell_converged==.false.)))/real(n_cells), "%"

        	!!lconverged = (real(diff) < precision)

        	
        	if (real(diff) < precision) then
           		if (lprevious_converged) then
            	  lconverged = .true.
           		else
            	  lprevious_converged = .true.
          	    endif
        	else
           		lprevious_converged = .false.
           		if (.not.lfixed_rays) then
              		n_rayons = n_rayons * 2
              		write(*,*) ' -- Increasing number of rays'
             		if (n_rayons > n_rayons_max) then
              			if (n_iter >= maxIter) then
             		 		write(*,*) "Warning : not enough rays to converge !!"
                 			lconverged = .true.
              			end if
              		end if

          	   end if
        	end if

	    end do !while
        write(*,*) etape, "Threshold =", precision
	  end do !over etapes
if (write_convergence_file ) close(20)
	  
  if (.not.lstop_after_jnu) then
  
!
!      call freeSpectrum
!      !redfine spectrum grid here
!      if (ltab_wavelength_image) NLTEspec%write_wavelength_grid = .true.
!      CALL initSpectrum(vacuum_to_air=lvacuum_to_air)
!      if ((atmos%NactiveAtoms>0) .or. .not.(ltab_wavelength_image)) then
!       if (n_etoiles > 0) CALL init_stellar_disk !for all wavelengths, all stars at disk centre
!       write(*,*) " recomputing background opacities after Jnu..."
!       CALL alloc_atom_quantities
!       CALL compute_opacities
!       write(*,*) " ..done"
!      endif
!      CALL reallocate_mcfost_vars()

   write(*,*) " -> Interpolating Jnu on image grid.."
   do icell=1, atmos%Nspace
  	CALL bezier2_interp(Nlambda, lambda, Jnew(:,icell), NLTEspec%Nwaves, NLTEspec%lambda, NLTEspec%Jc(:,icell))
   enddo
   write(*,*) "Jmax/min after interpolation on the image grid"
   write(*,*) maxval(NLTEspec%Jc(:,:)), minval(NLTEspec%Jc(:,:))
   write(*,*) " -> ..done"
  endif
  
  
  open(unit=20, file="Jnu_no_interp.s", status="unknown")
  write(20,*) n_cells, Nlambda
  do icell=1, n_cells
  do la=1, Nlambda
    write(20,'(6E)') lambda(la), Jnew(la,icell), Sth(la,icell), beta(la,icell), kappa_tot(la,icell), Sold(la,icell)
   enddo
  enddo
  close(20)


  if (ng_acc) CALL freeNg(NgJ)
  deallocate(Jold, Sth, kappa_tot,  beta, Jnew, Istar, Ic, lambda, Sold, lambda_star, Snew)
  if (allocated(Jflat)) deallocate(jflat)

 ! ------------------------------------------------------------------------------------ !
 RETURN
 END SUBROUTINE Iterate_Jnu
 
 SUBROUTINE INTEG_RAY_LINE_I_CNTRB(id,icell_in,x,y,z,u,v,w,iray,labs)
 ! ------------------------------------------------------------------------------- !
  ! Not Working yet
 ! ------------------------------------------------------------------------------- !
  integer, intent(in) :: id, icell_in, iray
  real(kind=dp), intent(in) :: u,v,w
  real(kind=dp), intent(in) :: x,y,z
  logical, intent(in) :: labs !used in NLTE but why?
  real(kind=dp) :: x0, y0, z0, x1, y1, z1, l, l_contrib, l_void_before
  real(kind=dp), dimension(NLTEspec%Nwaves) :: Snu, Snu_c, chiI, LimbD, chiIc
  real(kind=dp), dimension(NLTEspec%Nwaves) :: tau, tau_c, dtau_c, dtau, dtau_ds
  integer :: nbr_cell, icell, next_cell, previous_cell, icell_star, i_star
  logical :: lcellule_non_vide, lsubtract_avg, lintersect_stars, eval_operator


  RETURN
  END SUBROUTINE INTEG_RAY_LINE_I_CNTRB
  

END MODULE AtomicTransfer

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