! This module computes the metals bound-free opacitiy according
! to Kramer's law if Hydrogenic or using detailed photoionisation
! cross-section read in the atomic model.
!
! Note that Hydrogen bound-free is computed in the hydrogen.f90 module
!
! The module also computes the bound-bound opacity of metals (including
! Hydrogen this time) treated as PASSIVE atoms (either LTE populations or
! NLTE pops for a previous run) and sums all PASSIVE opacities in NLTEspec%chi and %eta
!
! chi in m^-1, eta in J/m3/s/Hz/sr
MODULE metal

 use atmos_type, only                : atmos, Hydrogen, Helium, B_project, vbroad_atom
 use constant
 use math, only                      : bezier3_interp, cent_deriv
 use atom_type
 use spectrum_type, only			 : NLTEspec, initAtomOpac
 use hydrogen_opacities
 use voigtfunctions, only 			 : Voigt
 use Profiles, only					 : Profile!, Profile_lambda
 use broad, only 					 : Damping
 use thomson_scattering
 use Rayleigh_scattering
 use Planck

 ! MCFOST's original
 use mcfost_env, only : dp
 use molecular_emission, only		 : v_proj
 use parametres
 use input
 use constantes, only				 : tiny_dp, huge_dp

 IMPLICIT NONE

 CONTAINS
 
 FUNCTION bound_free_Xsection(cont) result(alpha)
  type (AtomicContinuum) :: cont
  !real(kind=dp) :: n_eff, g_bf(cont%Nlambda), uu(cont%Nlambda), Z, uu0(1), g_bf0(1)
  real(kind=dp) :: alpha(cont%Nlambda) !note that, size(cont%alpha) /= cont%Nlambda
  									 !one represent the wavelength tabulated and the other one
  									 !the number of points on the total grid between Nblue and Nred
  									 
!   Z = real(cont%atom%stage(cont%i) + 1,kind=dp)
   if (cont%Hydrogenic) then !Kramer's formula with quantum mechanical correction
     alpha = H_bf_Xsection(cont)
!      if (cont%atom%ID == "H ") then
!       n_eff = dsqrt(Hydrogen%g(cont%i)/2.)  !only for Hydrogen !
!      else
!      !obtained_n = getPrincipal(metal%label(continuum%i), n_eff)
!      !if (.not.obtained_n) &
!         n_eff = Z*dsqrt(E_RYDBERG / (cont%atom%E(cont%j) - cont%atom%E(cont%i))) 
!      end if
!      sigma0_H
!      uu(:) = n_eff*n_eff*HPLANCK*CLIGHT / (NM_TO_M*NLTEspec%lambda(cont%Nblue:cont%Nred)) &
!      	/ (Z*Z) / E_RYDBERG - 1.
!      uu0(1) = n_eff*n_eff * HPLANCK*CLIGHT / (NM_TO_M * cont%lambda0) / Z / Z / E_RYDBERG - 1.
!        
!      g_bf0(:) = Gaunt_bf(1, uu0, n_eff)
!     
!      g_bf(:) = Gaunt_bf(cont%Nlambda, uu(:), n_eff)
! 
!      alpha = &
!         cont%alpha0 * g_bf(:) * (NLTEspec%lambda(cont%Nblue:cont%Nred)/cont%lambda0)**3  / g_bf0(1)!m^2
   else !interpolation of the read Cross-section
    !alpha = interp_dp(cont%alpha, cont%lambda, NLTEspec%lambda(cont%Nblue:cont%Nred))
    CALL bezier3_interp(size(cont%alpha), cont%lambda, cont%alpha, & !read values
     	cont%Nlambda, NLTEspec%lambda(cont%Nblue:cont%Nred), alpha) !interpolation grid

   endif
 
 RETURN
 END FUNCTION bound_free_Xsection

 SUBROUTINE Metal_bf(id, icell)
 !cross-section in cm2 per particle is given by Kramers’ formula
  !with n the principal quantum number of the level i from
 !which the atom or ion is ionized, Z the ion charge, ν in Hz and gbf the dimensionless
 !Gaunt factor, a quantummechanical correction factor of order unity.
 !The Kramers cross-section decays ∼ ν−3 above the threshold (“edge”) frequency ν0,
 !being zero below it because the threshold energy is the required minimum. Think the
  ! inverse in terms of wavelengths
  integer, intent(in)							            :: icell, id
  logical 										            :: obtained_n
  integer                                                   :: m, kr, i, j, Z, nc, Nblue, Nred
  type (AtomType)                                           :: metal
  type (AtomicContinuum)                                    :: continuum
  double precision                                          :: lambdaEdge
  double precision, dimension(:), allocatable               :: twohnu3_c2, gijk, Vij
   !obtained_n = .false. !true if the routine to read principal quantum number is fine
   
   if (atmos%Npassiveatoms==1 .and. atmos%PassiveAtoms(1)%ptr_atom%ID == "H ") RETURN
   ! else if H is not passive (active) it will not be the first.
   ! If H is passive it is always the first one (And if active the first one).


  ! Go throught all bound-free transitions of each PASSIVE
  ! metal and add the opacity and emissivity if lambda
  ! is lower (greater) than the wavelength threshold lambdaEdge
  ! (the frequency threshold) and if greater (lower) than
  ! wavelength min (frequency max). See Hydrogen b-f for more
  ! informations, and Hubeny & Mihalas chap. 7

  do m=1,atmos%Npassiveatoms
  ! run over all passive atoms
   metal = atmos%PassiveAtoms(m)%ptr_atom!atmos%Atoms(m)
   if (metal%ID == "H ") CYCLE !H cont is treated in Hydrogen_bf()

    do kr=1,metal%Ncont
     continuum = metal%continua(kr)
     i = continuum%i
     j = continuum%j !+1 wrt C indexing
     Nblue = continuum%Nblue; Nred = continuum%Nred
     if (.not.continuum%lcontrib_to_opac) CYCLE !avoid continua not defined on the grid

     !if (Nred == -99 .and. Nblue == -99) CYCLE !avoid continua not defined on the grid
     
     allocate(twohnu3_c2(continuum%Nlambda), gijk(continuum%Nlambda), Vij(continuum%Nlambda))
     !hc_kla = hc_k/NLTEspec%lambda(Nblue:Nred) !factor 1/NM_TO_M in hc_k
     twohnu3_c2 = twohc / NLTEspec%lambda(Nblue:Nred)**3
     Vij(:) = bound_free_Xsection(continuum)

     lambdaEdge = continuum%lambda0! or ionisation wavelength or wavelength
               ! associated to the minimal frquency needed
               ! to unbound an electron

    ! -> prevents dividing by zero
     ! even if lcompute_atomRT(icell) it is still possible to not have a continuum transition
     ! between the level i and j, but not for the others.
!	  if (metal%nstar(j,icell) <= tiny_dp .or. metal%nstar(i,icell)<=tiny_dp) CYCLE
    if (metal%nstar(j,icell) < tiny_dp) then
        write(*,*) "(Metal_bf) Warning at icell=", icell," T(K)=", atmos%T(icell)
        write(*,*) metal%ID,"%n(j) density <= tiny dp for j=", j, metal%n(j,icell), i, metal%n(i,icell)
        write(*,*) "skipping this level"
        write(*,*) "nstar=", metal%nstar(:,icell)
       CYCLE
    end if

     
     gijk(:) = metal%nstar(i,icell)/metal%nstar(j,icell) *  &
     			dexp(-hc_k/NLTEspec%lambda(Nblue:Nred)/atmos%T(icell))

!      write(*,*)
!       write(*,*) maxval(metal%n(i,icell)*(1.-expla(Nblue:Nred))), &
!       	maxval(metal%n(i,icell)-gijk(Nblue:Nred)*metal%n(i,icell))
!      stop
     if (lstore_opac) then !we don't care about proc id id
!       NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) = NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) + &
!        				continuum%alpha(Nblue:Nred) * &
!        				(1.-expla(Nblue:Nred))*metal%n(i,icell)
      NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) = NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) + &
       				Vij(:) * (metal%n(i,icell)-gijk(:)*metal%n(i,icell))
       				
      NLTEspec%AtomOpac%jc(icell,Nblue:Nred) = NLTEspec%AtomOpac%jc(icell,Nblue:Nred) + &
       				twohnu3_c2(:) * gijk(:) * Vij(:)*metal%n(j,icell)
     else !proc id is important
!       NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
!        				continuum%alpha(Nblue:Nred) * &
!        				(1.-expla(Nblue:Nred))*metal%n(i,icell) !-->ETL only
      NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
       				Vij(:) * (metal%n(i,icell)-gijk(:)*metal%n(j,icell))

      NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
      				twohnu3_c2(:) * gijk(:) * Vij(:)*metal%n(j,icell)
     end if

     deallocate(gijk, twohnu3_c2, Vij)
    end do ! loop over Ncont
  end do !loop over metals

 RETURN
 END SUBROUTINE Metal_bf
 
  SUBROUTINE Metal_bf_old(id, icell)
 !cross-section in cm2 per particle is given by Kramers’ formula
  !with n the principal quantum number of the level i from
 !which the atom or ion is ionized, Z the ion charge, ν in Hz and gbf the dimensionless
 !Gaunt factor, a quantummechanical correction factor of order unity.
 !The Kramers cross-section decays ∼ ν−3 above the threshold (“edge”) frequency ν0,
 !being zero below it because the threshold energy is the required minimum. Think the
  ! inverse in terms of wavelengths
  integer, intent(in)							            :: icell, id
  logical 										            :: obtained_n
  integer                                                   :: m, kr, i, j, Z, nc, Nblue, Nred
  type (AtomType)                                           :: metal
  type (AtomicContinuum)                                    :: continuum
  double precision                                          :: lambdaEdge, twohc, hc_k!, &
       														   !n_eff
  !double precision, dimension(1) 							:: gbf_0, uu0
  double precision, dimension(NLTEspec%Nwaves)              :: twohnu3_c2, gijk, hc_kla,&
                                                               expla!, alpha_la, uu, gbf
   !obtained_n = .false. !true if the routine to read principal quantum number is fine

   twohc = (2. * HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
   hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)
   
   if (atmos%Npassiveatoms==1 .and. atmos%PassiveAtoms(1)%ptr_atom%ID == "H ") RETURN
   ! else if H is not passive (active) it will not be the first.
   ! If H is passive it is always the first one (And if active the first one).

   hc_kla = hc_k/NLTEspec%lambda !factor 1/NM_TO_M in hc_k
   twohnu3_c2 = twohc / NLTEspec%lambda**3
   expla = dexp(-hc_kla/atmos%T(icell))

  ! Go throught all bound-free transitions of each PASSIVE
  ! metal and add the opacity and emissivity if lambda
  ! is lower (greater) than the wavelength threshold lambdaEdge
  ! (the frequency threshold) and if greater (lower) than
  ! wavelength min (frequency max). See Hydrogen b-f for more
  ! informations, and Hubeny & Mihalas chap. 7


  do m=1,atmos%Npassiveatoms
  ! run over all passive atoms
   metal = atmos%PassiveAtoms(m)%ptr_atom!atmos%Atoms(m)
   if (metal%ID == "H ") CYCLE !H cont is treated in Hydrogen_bf()

    do kr=1,metal%Ncont
     continuum = metal%continua(kr)
     i = continuum%i
     j = continuum%j !+1 wrt C indexing
     Nblue = continuum%Nblue; Nred = continuum%Nred
     if (.not.continuum%lcontrib_to_opac) CYCLE !avoid continua not defined on the grid

     lambdaEdge = continuum%lambda0! or ionisation wavelength or wavelength
               ! associated to the minimal frquency needed
               ! to unbound an electron

    ! -> prevents dividing by zero
     ! even if lcompute_atomRT(icell) it is still possible to not have a continuum transition
     ! between the level i and j, but not for the others.
!	  if (metal%nstar(j,icell) <= tiny_dp .or. metal%nstar(i,icell)<=tiny_dp) CYCLE
    if (metal%nstar(j,icell) <= tiny_dp) then
       if (metal%nstar(j,icell) <= 0) then
        write(*,*) "(Metal_bf) Warning at icell=", icell," T(K)=", atmos%T(icell)
        write(*,*) metal%ID,"%n(j) density <= tiny dp for j=", j, metal%n(j,icell)
        write(*,*) "skipping this level"
       end if
       CYCLE
    end if

!      Z = metal%stage(i) + 1
!      !! Only for Hydrogen n_eff = dsqrt(metal%g(i)/2.)
!      !obtained_n = getPrincipal(metal%label(continuum%i), n_eff)
! 
!      if (.not.obtained_n) &
!         n_eff = Z*dsqrt(E_RYDBERG / (metal%E(continuum%j) - metal%E(continuum%i)))
! 
!      ! for this continuum of this line
!      alpha_la = 0d0
!      uu = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*Z*Z * E_RYDBERG * NLTEspec%lambda)
!      uu0 = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*Z*Z * E_RYDBERG * lambdaEdge)
!      if (.not.continuum%Hydrogenic) then
! 
!        CALL bezier3_interp(continuum%Nlambda,&
!                            continuum%lambda, &
!                            continuum%alpha,  & !Now the interpolation grid
!              continuum%Nlambda, &
!              NLTEspec%lambda(continuum%Nblue:continuum%Nred), &
!              alpha_la(continuum%Nblue:continuum%Nred)) !end
! 
!      else
!      
!         gbf_0 = Gaunt_bf(1, uu0, n_eff)
!          !continuum%alpha0 is the result of the photoionisation
!          !cross-section at lambdaEdge=continuum%lambda0
!          ! = HPLANCK*CLIGHT/(Ej - Ei)
!          ! therefore we scale wrt the lambdaEdge
!          ! alpha_la(lambda=lambdaEdge)=alpha0 (containing
!          ! already the gaunt factor !! and factor 1/Z^2)
!            
!         gbf(continuum%Nblue:continuum%Nred) = Gaunt_bf(continuum%Nlambda, &
!           uu(continuum%Nblue:continuum%Nred), n_eff)
! 
!         alpha_la(continuum%Nblue:continuum%Nred) = gbf(continuum%Nblue:continuum%Nred) * &
!            continuum%alpha0*((NLTEspec%lambda(continuum%Nblue:continuum%Nred)/lambdaEdge)**3)*&
!            n_eff/gbf_0(1)
! 
!      end if !continuum type
     
     gijk(Nblue:Nred) = metal%nstar(i,icell)/metal%nstar(j,icell) * expla(Nblue:Nred)
!      write(*,*)
!       write(*,*) maxval(metal%n(i,icell)*(1.-expla(Nblue:Nred))), &
!       	maxval(metal%n(i,icell)-gijk(Nblue:Nred)*metal%n(i,icell))
!      stop
     if (lstore_opac) then !we don't care about proc id id
!       NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) = NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) + &
!        				continuum%alpha(Nblue:Nred) * &
!        				(1.-expla(Nblue:Nred))*metal%n(i,icell)
      NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) = NLTEspec%AtomOpac%Kc(icell,Nblue:Nred,1) + &
       				continuum%alpha(Nblue:Nred) * &
       				(metal%n(i,icell)-gijk(Nblue:Nred)*metal%n(i,icell))
       				
      NLTEspec%AtomOpac%jc(icell,Nblue:Nred) = NLTEspec%AtomOpac%jc(icell,Nblue:Nred) + &
       				twohnu3_c2(Nblue:Nred) * gijk(Nblue:Nred) * &
         			continuum%alpha(Nblue:Nred)*metal%n(j,icell)
     else !proc id is important
!       NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
!        				continuum%alpha(Nblue:Nred) * &
!        				(1.-expla(Nblue:Nred))*metal%n(i,icell) !-->ETL only
      NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
       				continuum%alpha(Nblue:Nred) * &
       				(metal%n(i,icell)-gijk(Nblue:Nred)*metal%n(j,icell))

      NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
      				twohnu3_c2(Nblue:Nred) * gijk(Nblue:Nred) * &
         			continuum%alpha(Nblue:Nred)*metal%n(j,icell)
     end if

    end do ! loop over Ncont
  end do !loop over metals

 RETURN
 END SUBROUTINE Metal_bf_old

 SUBROUTINE Metal_bb (id, icell,x,y,z,x1,y1,z1,u,v,w,l)
  ! Computes the emissivity and extinction of passive lines.
  ! i.e., Atoms with detailed atomic structures read but
  ! not treated in NLTE.
  ! Because damping is wavelength independent and depend only on
  ! the grid (cell) points, here, if line%damping_initialized
  ! do not CALL Damping()
  ! the x,y,z and u,v,w quantities are used to compute the projected velocities at the
  ! cell point we are computing the opacities.
  ! Chip is only computed in Stokes transfer and contains the magneto-optical elements.
  integer 													:: kr, m, i, j, nk
  integer, intent(in) 							            :: icell, id
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  double precision, dimension(:), allocatable				:: Vij
  double precision 											:: twohnu3_c2, gij
  integer													:: Nred, Nblue
  type (AtomicLine)										    :: line
  type (AtomType)											:: atom
  double precision, dimension(:,:), allocatable 			:: psiZ, phiZ
  double precision, allocatable, dimension(:) 				:: phi(:)


  do m=1,atmos%Npassiveatoms
   atom = atmos%PassiveAtoms(m)%ptr_atom
    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i; j = line%j
     Nred = line%Nred; Nblue = line%Nblue
     if (.not.line%lcontrib_to_opac) CYCLE !avoid lines not defined on the grid
     allocate(Vij(line%Nlambda))

!     if ((atom%n(j,icell) <=0).or.(atom%n(i,icell) <=0)) CYCLE !"no contrib to opac"
     ! -> prevents dividing by zero
     ! even if lcompute_atomRT(icell) it is still possible to not have a transition
     ! between the levels i and j, but not for the others.

     if ((atom%n(j,icell) <tiny_dp).or.(atom%n(i,icell) <tiny_dp)) then !no transition
       !!but show the message only if pops is negative
      !if ((atom%n(j,icell) < 0 ).or.(atom%n(i,icell) < 0)) then
        write(*,*) "(Metal_bb) Warning at icell=", icell," T(K)=", atmos%T(icell)
        write(*,*) atom%ID," density <= tiny dp ", i, j, line%lambda0, atom%n(i,icell), atom%n(j,icell)
        write(*,*) "skipping this level"
        write(*,*) "nstar=", atom%nstar(:,icell)
        write(*,*) "n = ", atom%n(:,icell)
      !!end if
      CYCLE
     end if

     !allocate(Vij(line%Nlambda)); Vij = 0d0
     gij = line%Bji / line%Bij ! = gi/gj
     twohnu3_c2 = line%Aji / line%Bji

     allocate(phi(line%Nlambda))!; phi = 0d0
     if (PRT_SOLUTION == "FULL_STOKES") allocate(phiZ(3,line%Nlambda), psiZ(3,line%Nlambda))
     !write(*,*) allocated(phiZ), allocated(psiZ), line%polarizable, PRT_SOLUTION
     if (line%voigt) CALL Damping(icell, atom, kr, line%adamp)
     CALL Profile (line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phiZ, psiZ)

     !Sum up all contributions for this line with the other
     Vij(:) = hc_4PI * line%Bij * phi(:)!already normalized / (SQRTPI * VBROAD_atom(icell,atom))
      
     NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + Vij(:) * (atom%n(i,icell)-gij*atom%n(j,icell))

     NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + twohnu3_c2 * gij * Vij(:) * atom%n(j,icell)

!          write(*,*) "check"
!          write(*,*) gij,  atom%g(i)/atom%g(j)
!          write(*,*) atom%g(i)/atom%g(j)*atom%n(j,icell)/atom%n(i,icell), exp(-hc/KBOLTZMANN/atmos%T(icell)/NM_TO_M/line%lambda0)
!          stop

     if (line%polarizable .and. PRT_SOLUTION == "FULL_STOKES") then
       do nk = 1, 3
         !magneto-optical
         NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%rho_p(Nblue:Nred,nk,id) + &
           hc_4PI * line%Bij * (atom%n(i,icell)-gij*atom%n(j,icell)) * psiZ(nk,:)
         !dichroism
         NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%chiQUV_p(Nblue:Nred,nk,id) + &
           hc_4PI * line%Bij * (atom%n(i,icell)-gij*atom%n(j,icell)) * psiZ(nk,:)
         !emissivity
         NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) = NLTEspec%AtomOpac%etaQUV_p(Nblue:Nred,nk,id) + &
          twohnu3_c2 * gij * hc_4PI * line%Bij * atom%n(j,icell) * phiZ(nk,:)
       end do 
     end if
     
     deallocate(Vij,phi)
     !deallocate(phi)
     if (PRT_SOLUTION == "FULL_STOKES") deallocate(psiZ, phiZ)
    end do !end loop on lines for this atom
  end do !end loop over Natom

 RETURN
 END SUBROUTINE Metal_bb
 SUBROUTINE Metal_bb_old (id, icell,x,y,z,x1,y1,z1,u,v,w,l)
  ! Computes the emissivity and extinction of passive lines.
  ! i.e., Atoms with detailed atomic structures read but
  ! not treated in NLTE.
  ! Because damping is wavelength dependent and depend only on
  ! the grid (cell) points, here, if line%damping_initialized
  ! do not CALL Damping()
  ! the x,y,z and u,v,w quantities are used to compute the projected velocities at the
  ! cell point we are computing the opacities.
  ! Chip is only computed in Stokes transfer and contains the magneto-optical elements.
  integer 													:: kr, m, i, j
  integer, intent(in) 							            :: icell, id
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  double precision, dimension(NLTEspec%Nwaves)              :: phi, vvoigt, phip, &
 															   Vij, vv
  double precision 											:: twohnu3_c2, hc, fourPI, &
      														   hc_4PI, gij
  integer, parameter										:: NvspaceMax = 101
  double precision, dimension(NvspaceMax)					:: omegav
  integer													:: Nvspace, nv, Nred, Nblue
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv
  type (AtomicLine)										    :: line
  type (AtomType)											:: atom

  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI

  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   v0 = v_proj(icell,x,y,z,u,v,w)
   omegav(1) = v0
  end if

  do m=1,atmos%Npassiveatoms
   atom = atmos%PassiveAtoms(m)%ptr_atom
   
    !velocity projected along a path between one border of the cell to the other
    if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0) 
     !write(*,*) v0/1000, v1/1000, dv/1000
     Nvspace = max(2,nint(20*dv/VBROAD_atom(icell,atom)))
     Nvspace = min(Nvspace,NvspaceMax) !Ensure that we never have an allocation error
     omegav(Nvspace) = v1
     !write(*,*) Nvspace, Nvspacemax
     do nv=2,Nvspace-1
      delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
      xphi=x+delta_vol_phi*u
      yphi=y+delta_vol_phi*v
      zphi=z+delta_vol_phi*w
      omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
      !write(*,*) icell, omegav(nv)/1d3
     end do 
    end if

    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i; j = line%j
     Nred = line%Nred; Nblue = line%Nblue
     if (Nred == -99 .and. Nblue == -99) CYCLE !avoid lines not defined on the grid

!     if ((atom%n(j,icell) <=0).or.(atom%n(i,icell) <=0)) CYCLE !"no contrib to opac"
     ! -> prevents dividing by zero
     ! even if lcompute_atomRT(icell) it is still possible to not have a transition
     ! between the levels i and j, but not for the others.
     if ((atom%n(j,icell) <=tiny_dp).or.(atom%n(i,icell) <=tiny_dp)) then !no transition
       write(*,*) "(Metal_bb) Warning at icell=", icell," T(K)=", atmos%T(icell)
       write(*,*) atom%ID," density <= tiny dp ", i, j, line%lambda0, atom%n(i,icell), atom%n(j,icell)
       write(*,*) "skipping this level"
      CYCLE
     end if

     phi = 0d0
     ! line dependent only
     vv(Nblue:Nred) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * VBROAD_atom(icell,atom))

     gij = line%Bji / line%Bij
     twohnu3_c2 = line%Aji / line%Bji
     
     !if (line%voigt) CALL Damping(icell, atom, kr, line%adamp)
     !CALL Iprofile (line, icell,x,y,z,x1,y1,z1,u,v,w,l, phi, phip)
     
     if (line%voigt) then
      !some work to do here if line%damping_initialized = .true.==kept on the whole grid.
      CALL Damping(icell, atom, kr, line%adamp)
       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
        vvoigt(Nblue:Nred) = vv(Nblue:Nred) - &
                                       omegav(nv) / VBROAD_atom(icell,atom)
        !loop over components here

        phi(Nblue:Nred) = phi(Nblue:Nred) + &
            Voigt(line%Nlambda, line%adamp,vvoigt(Nblue:Nred), &
                  phip, VoigtMethod) / Nvspace !1 if no vel or Voronoi

      end do
     else !Gaussian !only for checking
      do nv=1, Nvspace
        vvoigt(Nblue:Nred) = vv(Nblue:Nred) - omegav(nv) / VBROAD_atom(icell,atom)
       phi(Nblue:Nred) = phi(Nblue:Nred) + dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace

      end do
     end if !line%voigt

     !Sum up all contributions for this line with the other

     Vij(Nblue:Nred) = &
      hc_4PI * line%Bij * phi(Nblue:Nred) / (SQRTPI * VBROAD_atom(icell,atom))
      
     NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
       		Vij(Nblue:Nred) * (atom%n(i,icell)-gij*atom%n(j,icell))
       		
     NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
       		twohnu3_c2 * gij * Vij(Nblue:Nred) * atom%n(j,icell)

    end do !end loop on lines for this atom
  end do !end loop over Natom

 RETURN
 END SUBROUTINE Metal_bb_old
 
 SUBROUTINE MetalZeeman_bb (id, icell,x,y,z,x1,y1,z1,u,v,w,l)
  integer 													:: kr, m, i, j, NrecStokes
  integer, intent(in) 							            :: icell, id
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  double precision, dimension(NLTEspec%Nwaves)              :: phi, vvoigt, phiPol, phip, &
 															   Vij, vv
  double precision 											:: twohnu3_c2, hc, fourPI, &
      														   hc_4PI, gij
  integer, parameter										:: NvspaceMax = 101, NbspaceMax=15
  double precision, dimension(NvspaceMax)					:: omegav
  double precision, dimension(NbspaceMax)					:: omegaB, gamma, chi
  integer													:: Nvspace, nv, Nred, Nblue, nc, &
  															   Nbspace, nb, Nzc
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv, dlamB, b0, b1,g1,c1,dB
  type (AtomicLine)										    :: line
  type (AtomType)											:: atom

  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI
  omegaB = 0d0

  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   v0 = v_proj(icell,x,y,z,u,v,w)
   omegav(1) = v0
  end if


  b0 = B_project(icell,x,y,z,u,v,w,g1,c1)
  omegaB(1) = b0
  gamma(1) = g1; chi(1)=c1

  do m=1,atmos%Npassiveatoms
   atom = atmos%PassiveAtoms(m)%ptr_atom
   
    if (.not.lstatic .and. .not.lVoronoi) then
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0) 
     Nvspace = max(2,nint(20*dv/VBROAD_atom(icell,atom)))
     Nvspace = min(Nvspace,NvspaceMax)
     omegav(Nvspace) = v1
     do nv=2,Nvspace-1
      delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
      xphi=x+delta_vol_phi*u
      yphi=y+delta_vol_phi*v
      zphi=z+delta_vol_phi*w
      omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
     end do 
    end if
    if (.not.lvoronoi) then
      b1 = B_project(icell,x1,y1,z1,u,v,w,g1,c1)
      Nbspace = NbspaceMax
      omegaB(Nbspace) = b1
      gamma(Nbspace) = g1; chi(Nbspace)=c1
      do nv=2,Nbspace-1
       delta_vol_phi = (real(nv,kind=dp))/(real(Nbspace,kind=dp)) * l
       xphi=x+delta_vol_phi*u
       yphi=y+delta_vol_phi*v
       zphi=z+delta_vol_phi*w
       omegaB(nv) = B_project(icell,xphi,yphi,zphi,u,v,w,g1,c1)
       gamma(nv) = g1; chi(nv)=c1
      end do      
    end if
    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i; j = line%j
     Nred = line%Nred; Nblue = line%Nblue
     Nzc = 0
     if (line%polarizable) Nzc = line%zm%Ncomponent
     if (.not.line%voigt) then
      CALL Warning("Skipping line because only Voigt profile for Zeeman calculation!")
      CYCLE
     end if

     if ((atom%n(j,icell) <=tiny_dp).or.(atom%n(i,icell) <=tiny_dp)) then !no transition
       write(*,*) "(Metal_bb) Warning at icell=", icell," T(K)=", atmos%T(icell)
       write(*,*) atom%ID," density <= tiny dp ", i, j, line%lambda0, atom%n(i,icell), atom%n(j,icell)
       write(*,*) "skipping this level"
      CYCLE
     end if


     phi = 0d0
     phip = 0d0
     phiPol = 0d0
     ! line dependent only
     vv(Nblue:Nred) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * VBROAD_atom(icell,atom))

     gij = line%Bji / line%Bij
     twohnu3_c2 = line%Aji / line%Bji
     
      !some work to do here if line%damping_initialized = .true.==kept on the whole grid.
      CALL Damping(icell, atom, kr, line%adamp)

       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
        vvoigt(Nblue:Nred) = vv(Nblue:Nred) - &
                                       omegav(nv) / VBROAD_atom(icell,atom)
         do nb=1,Nbspace
          if (Nzc == 0) & !line not polarizable or weak field
             phi(Nblue:Nred) = phi(Nblue:Nred) + &
          					Voigt(line%Nlambda, line%adamp,vvoigt(Nblue:Nred), &
                  			phip, VoigtMethod) / Nvspace / Nbspace
          do nc=1,Nzc
             vvoigt(Nblue:Nred) = vvoigt(Nblue:Nred) - omegaB(nb) / VBROAD_atom(icell,atom) * &
                                  LARMOR/CLIGHT * (line%lambda0 * NM_TO_M)**2 * line%zm%shift(nc)
             phi(Nblue:Nred) = phi(Nblue:Nred) + &
          					Voigt(line%Nlambda, line%adamp,vvoigt(Nblue:Nred), &
                  			phip, VoigtMethod) / Nvspace / Nbspace!1 if no vel or Voronoi
             phiPol(Nblue:Nred) = phiPol(Nblue:Nred) + phip(Nblue:Nred)/Nvspace/Nbspace
             !do something here
             CALL Error("Full profile not implemented")
          end do !components 
        end do !magnetic field     
       end do !velocity

     !Sum up all contributions for this line with the other

     Vij(Nblue:Nred) = &
      hc_4PI * line%Bij * phi(Nblue:Nred) / (SQRTPI * VBROAD_atom(icell,atom))
      
     NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%chi_p(Nblue:Nred,id) + &
       		Vij(Nblue:Nred) * (atom%n(i,icell)-gij*atom%n(j,icell))
       		
     NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) = &
     		NLTEspec%AtomOpac%eta_p(Nblue:Nred,id) + &
       		twohnu3_c2 * gij * Vij(Nblue:Nred) * atom%n(j,icell)
     
       
    end do !end loop on lines for this atom
  end do !end loop over Natom

 RETURN
 END SUBROUTINE MetalZeeman_bb

 !--> I should include the zeeman lines in phi, because it plays a role in opacity
 !and tau
 SUBROUTINE Metal_bb_lambda (id, la, icell,x,y,z,x1,y1,z1,u,v,w,l)
  integer 													:: kr, m, i, j, NrecStokes
  integer, intent(in) 							            :: icell, la, id
  double precision, intent(in) 					            :: x,y,z,u,v,w,& 
                                				               x1,y1,z1,l
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  double precision, dimension(1)                            :: phi, vvoigt, phiPol, phip, &
 															   Vij, vv
  double precision 											:: twohnu3_c2, hc, fourPI, &
      														   hc_4PI, gij
  integer, parameter										:: NvspaceMax = 101
  double precision, dimension(NvspaceMax)					:: omegav
  integer													:: Nvspace, nv
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv
  type (AtomicLine)										    :: line
  type (AtomType)											:: atom

  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI

  !check that NLTEspec%AtomOpac%eta_p(id,:) = 0d0 here

  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
    v0 = v_proj(icell,x,y,z,u,v,w)
    omegav(1) = v0
  end if


  !do m=1,atmos%Natom ! go over all atoms
  do m=1,atmos%Npassiveatoms
   atom = atmos%PassiveAtoms(m)%ptr_atom!atmos%Atoms(m)
   !if (atom%active) CYCLE ! go to next passive atom
   
    !velocity projected along a path between one border of the cell to the other
    !if (.not.lstatic .and. .not.atmos%Voronoi .and.lmagnetoaccr) then ! velocity is varying across the cell
    if (.not.lstatic .and. .not.lVoronoi) then
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0) 
     Nvspace = max(2,nint(20*dv/VBROAD_atom(icell,atom)))
     Nvspace = min(Nvspace,NvspaceMax) !Ensure that we never have an allocation error
     omegav(Nvspace) = v1
     do nv=2,Nvspace-1
      delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
      xphi=x+delta_vol_phi*u
      yphi=y+delta_vol_phi*v
      zphi=z+delta_vol_phi*w
      omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
     end do 
    end if

    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i
     j = line%j
     
     !if (line%Nred == -99 .and. line%Nblue == -99) CYCLE
     
     !wavelength does not fall inside line domaine? OK cylcle
     if ((NLTEspec%lambda(la) < NLTEspec%lambda(line%Nblue)).or.&
        (NLTEspec%lambda(la) > NLTEspec%lambda(line%Nred))) then
!        write(*,*) NLTEspec%lambda(la), " beyond line edge:", j, i
!        write(*,*) NLTEspec%lambda(line%Nblue), NLTEspec%lambda(line%Nred)
!        write(*,*) "Skipping"
       CYCLE
     end if


! Checked in LTEpops and NLTE loop
   if ((atom%n(j,icell) <= tiny_dp).or.(atom%n(i,icell) <= tiny_dp)) then !no transition
       write(*,*) "(Metal_bb_lambda) Warning at icell=", icell," T(K)=", atmos%T(icell)
       write(*,*) atom%ID," density <= tiny dp ", i, j, line%lambda0, atom%n(i,icell), atom%n(j,icell)
       write(*,*) "skipping this level"
      CYCLE
    end if

     phi = 0d0
     phip = 0d0
     phiPol = 0d0
     ! line dependent only
     vv(1) = (NLTEspec%lambda(la)-line%lambda0) * &
           CLIGHT / (line%lambda0 * VBROAD_atom(icell,atom))

     gij = line%Bji / line%Bij
     twohnu3_c2 = line%Aji / line%Bji
     
     if (line%voigt) then
      !some work to do here if line%damping_initialized = .true.==kept on the whole grid.
      CALL Damping(icell, atom, kr, line%adamp)
       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace 
        vvoigt(1) = vv(1) -  omegav(nv) / VBROAD_atom(icell,atom)

        phi = phi + Voigt(1, line%adamp,vvoigt,phip, VoigtMethod) / Nvspace 
      end do
     else !Gaussian
      do nv=1, Nvspace
        vvoigt(1) = vv(1) - omegav(nv) / VBROAD_atom(icell,atom)
       phi(1) = phi(1) + dexp(-(vvoigt(1))**2) / Nvspace
      end do
     end if !line%voigt
    

     Vij(1) = hc_4PI * line%Bij * phi(1) / (SQRTPI * VBROAD_atom(icell,atom))
     NLTEspec%AtomOpac%chi_p(la,id) = NLTEspec%AtomOpac%chi_p(la,id) +&
     								  Vij(1) * (atom%n(i,icell)-gij*atom%n(j,icell))
     NLTEspec%AtomOpac%eta_p(la,id) = NLTEspec%AtomOpac%eta_p(la,id) +&
     								  twohnu3_c2 * gij * Vij(1) * atom%n(j,icell)

    end do !end loop on lines for this atom
  end do !end loop over Natom

 RETURN
 END SUBROUTINE Metal_bb_lambda
 
 SUBROUTINE storeBackground()
 !$ use omp_lib
  integer :: icell, id
  if (.not.lstore_opac) RETURN

   if (real(3*n_cells*NLTEspec%Nwaves)/(1024**3) < 1.) then
    write(*,*) "Keeping", real(3*n_cells*NLTEspec%Nwaves)/(1024**2), " MB of memory", &
   	 " for Background continuum opacities."
   else
    write(*,*) "Keeping", real(3*n_cells*NLTEspec%Nwaves)/(1024**3), " GB of memory", &
   	 " for Background continuum opacities."
   end if
   !$omp parallel &
   !$omp default(none) &
   !$omp private(icell,id) &
   !$omp shared(atmos)
   !$omp do schedule(dynamic,1)
   do icell=1,atmos%Nspace
    !$ id = omp_get_thread_num() + 1
    CALL BackgroundContinua(icell)
   end do
   !$omp end do
   !$omp end parallel

 RETURN
 END SUBROUTINE storeBackground
 
 SUBROUTINE BackgroundContinua (icell)
  integer, intent(in) :: icell
  double precision, dimension(NLTEspec%Nwaves) :: chi, eta, Bpnu!,sca

  if (atmos%icompute_atomRT(icell)<1) RETURN
  
   CALL Bplanck(atmos%T(icell), Bpnu)

   chi = 0d0
   eta = 0d0

  NLTEspec%AtomOpac%Kc(icell,:,1) = Thomson(icell)

   CALL Rayleigh(1, icell, Hydrogen)
   if (associated(Helium)) CALL Rayleigh(1, icell, Helium)

   NLTEspec%AtomOpac%Kc(icell,:,2) = NLTEspec%AtomOpac%Kc(icell,:,1)

   CALL Hydrogen_ff(icell, chi)
   NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + chi
   NLTEspec%AtomOpac%jc(icell,:) = NLTEspec%AtomOpac%jc(icell,:) + chi * Bpnu

   
   CALL Hminus_bf(icell, chi,eta)
   NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + chi
   NLTEspec%AtomOpac%jc(icell,:) = NLTEspec%AtomOpac%jc(icell,:) + eta


   CALL Hminus_ff(icell, chi)
   NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + chi
   NLTEspec%AtomOpac%jc(icell,:) = NLTEspec%AtomOpac%jc(icell,:) + chi * Bpnu

   if (.not.Hydrogen%active) then !passive bound-free
    CALL Hydrogen_bf(icell, chi, eta)
    NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + chi
    NLTEspec%AtomOpac%jc(icell,:) = NLTEspec%AtomOpac%jc(icell,:) + eta
   end if

    if (atmos%Npassiveatoms == 0) RETURN !no passive bound-bound and bound-free
    CALL Metal_bf(1, icell) !here we do not need id, because lstore_atom_opac=.true.


 RETURN
 END SUBROUTINE BackgroundContinua
 
 SUBROUTINE BackgroundLines(id,icell,x,y,z,x1,y1,z1,u,v,w,l)
  integer, intent(in) :: icell, id
  double precision, intent(in) :: x, y, z, u, v, w, &
                                  x1, y1, z1, l

  if (atmos%icompute_atomRT(icell)<1) RETURN 

   if (atmos%Npassiveatoms == 0) RETURN
   CALL Metal_bb(id, icell, x, y, z, x1, y1, z1, u, v, w, l)


 RETURN
 END SUBROUTINE BackgroundLines
 
 SUBROUTINE BackgroundLines_lambda(la,id,icell,x,y,z,x1,y1,z1,u,v,w,l)
  !Only for stellar calculations so I assume no Zeeman polarisation from the star
  integer, intent(in) :: icell, id, la
  double precision, intent(in) :: x, y, z, u, v, w, &
                                  x1, y1, z1, l

  if (atmos%icompute_atomRT(icell)<1) RETURN 


   if (atmos%Npassiveatoms == 0) RETURN
   CALL Metal_bb_lambda(id, la, icell, x, y, z, x1, y1, z1, u, v, w, l)

 RETURN
 END SUBROUTINE BackgroundLines_lambda

 SUBROUTINE Background(id,icell,x,y,z,x1,y1,z1,u,v,w,l)
  integer, intent(in) :: icell, id
  double precision, intent(in) :: x, y, z, u, v, w, &
                                  x1, y1, z1, l!only relevant for b-b when vector fields are present
  double precision, dimension(NLTEspec%Nwaves) :: chi, eta, Bpnu, chip!, sca

  if (atmos%icompute_atomRT(icell)<1) RETURN !nH <= tiny_nH or T <= tiny_T == empty cell
  ! all opac are zero, return.

!   if ((atmos%nHtot(icell)==0d0).or.(atmos%T(icell)==0d0)) &
!     RETURN ! stoping for this cell,
            ! it is free of (significant) gas
            ! so no emission/absorption. Coefficients set to 0d0 for all wavelengths
  !Do not forget that it is still possible however, that lcompute_atomRT is TRUE,
  !but that the temperature is too low to have all the levels non zero.
  !The compute_atomRT ensures that at least one level has non-zero populations,
  !to avoid division by zero.

   CALL Bplanck(atmos%T(icell), Bpnu)

   chi = 0d0
   eta = 0d0
!   sca = 0d0
   chip = 0d0


!    if (Rayleigh(icell, Hydrogen, sca)) NLTEspec%AtomOpac%sca_c(id,:) = sca
!    NLTEspec%AtomOpac%sca_c(id,:) =  NLTEspec%AtomOpac%sca_c(id,:) + Thomson(icell)
!    if (associated(Helium)) then
!     if (Rayleigh(icell, Helium, sca)) NLTEspec%AtomOpac%sca_c(id,:) = &
!           NLTEspec%AtomOpac%sca_c(id,:) + sca
!    end if
   !CALL HRayleigh(icell,Hydrogen, sca)
   !NLTEspec%AtomOpac%sca_c(id,:) = NLTEspec%AtomOpac%sca_c(id,:) + sca
   NLTEspec%AtomOpac%sca_c(:,id) = Thomson(icell) !init then update
   CALL Rayleigh(id, icell, Hydrogen)
   if (associated(Helium)) CALL Rayleigh(id, icell, Helium)


   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%sca_c(:,id)

   CALL Hydrogen_ff(icell, chi)
   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
   NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + chi * Bpnu

   CALL Hminus_bf(icell, chi, eta)
   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
   NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + eta

   CALL Hminus_ff(icell, chi)
   NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
   NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + chi * Bpnu

   if (.not.Hydrogen%active) then !passive bound-free !do not enter if active !!!
    CALL Hydrogen_bf(icell, chi, eta)
    NLTEspec%AtomOpac%chi_p(:,id) = NLTEspec%AtomOpac%chi_p(:,id) + chi
    NLTEspec%AtomOpac%eta_p(:,id) = NLTEspec%AtomOpac%eta_p(:,id) + eta
   end if
   
   if (atmos%Npassiveatoms == 0) RETURN !no passive bound-bound and bound-free
   ! we avoid H passive in metal_bf if H is passive
   						!, chi, eta)
   						
   !--> at this point, eta_p and chi_p are 0 because of initAtomOpac(id), therefore
   !after metal_bf they only points to continuum bound-free.
    CALL Metal_bf(id, icell) !Return if Npassive=1 and PassiveAtoms(1)==" H"
!     NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
!     NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + eta

   !keep pure continuum opacities now
   if (MINVAL(NLTEspec%AtomOpac%eta_p(:,id))<0 .or. &
    MINVAL(NLTEspec%AtomOpac%chi_p(:,id)) < 0) then
    write(*,*) "err, negative opac"
    stop
   end if
   NLTEspec%AtomOpac%eta_c(:,id) = NLTEspec%AtomOpac%eta_p(:,id)
   NLTEspec%AtomOpac%chi_c(:,id) = NLTEspec%AtomOpac%chi_p(:,id)

   ! we already RETURNs if no passive transitions (H included)
   CALL Metal_bb(id, icell, x, y, z, x1, y1, z1, u, v, w, l)


 RETURN
 END SUBROUTINE Background

END MODULE metal
