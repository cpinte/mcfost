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

 use atmos_type, only : atmos, Hydrogen, Helium, atmos_parcel
 use constant
 use math, only : bezier3_interp
 use atom_type
 use spectrum_type, only : NLTEspec, initAtomOpac
 use hydrogen_opacities
 use voigtfunctions, only : Voigt
 use broad, only : Damping
 use thomson_scattering
 use Rayleigh_scattering
 use Planck

 ! MCFOST's original
 use molecular_emission, only : v_proj

 IMPLICIT NONE


 CONTAINS

 FUNCTION Metal_bf(icell, chi, eta) result(res)
 !cross-section in cm2 per particle is given by Kramers’ formula
  !with n the principal quantum number of the level i from
 !which the atom or ion is ionized, Z the ion charge, ν in Hz and gbf the dimensionless
 !Gaunt factor, a quantummechanical correction factor of order unity.
 !The Kramers cross-section decays ∼ ν−3 above the threshold (“edge”) frequency ν0,
 !being zero below it because the threshold energy is the required minimum. Think the
  ! inverse in terms of wavelengths
  integer :: icell
  logical :: res, obtained_n
  integer :: m, kr, i, j, Z, nc
! integer, dimension(:), allocatable :: iLam
  type (AtomType) :: metal
  type (AtomicContinuum) :: continuum
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: eta, chi
  double precision :: lambdaEdge, twohc, hc_k, n_eff
  double precision, dimension(NLTEspec%Nwaves) :: gbf_0, uu, gbf
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2, gijk, hc_kla, expla, alpha_la

  res = .true. !always now, return condition tested outside
  obtained_n = .false. !true if the routine to read principal quantum number is fine

  twohc = (2. * HPLANCK * CLIGHT) / (NM_TO_M)**(3d0)
  hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)

   chi = 0.
   eta = 0.

   hc_kla = hc_k/NLTEspec%lambda !factor 1/NM_TO_M in hc_k
   twohnu3_c2 = twohc / NLTEspec%lambda**3
   expla = dexp(-hc_kla/atmos%T(icell))

  ! Go throught all bound-free transitions of each PASSIVE
  ! metal and add the opacity and emissivity if lambda
  ! is lower (greater) than the wavelength threshold lambdaEdge
  ! (the frequency threshold) and if greater (lower) than
  ! wavelength min (frequency max). See Hydrogen b-f for more
  ! informations, and Hubeny & Mihalas chap. 7

  ! m=2 because we avoid Hydrogen =)
  do m=2,atmos%Natom
  ! run over all passive atoms
   metal = atmos%Atoms(m)
   if (metal%active) CYCLE ! go to next passive metal
   !if loc(%n)=loc(%nstar) pure ETL, or if associated(%n,%nstar)
   !else NLTEpops read from previous run for actual passive metal
   ! I use atom%n(level,icell) !if NLTE pops exist for this atom, atom%n != atom%nstar
   ! else atom%n = atom%nstar.
   !write(*,*) loc(metal%n), loc(metal%nstar), &
   !                  associated(metal%n,metal%nstar)


    uu = 0d0
    do kr=1,metal%Ncont
     continuum = metal%continua(kr)
     i = continuum%i
     j = continuum%j !+1 wrt C indexing
     lambdaEdge = continuum%lambda0! or ionisation wavelength or wavelength
               ! associated to the minimal frquency needed
               ! to unbound an electron

    ! -> prevents dividing by zero
     ! even if lcompute_atomRT(icell) it is still possible to not have a continuum transition
     ! between the level i and j, but not for the others.
    if (metal%nstar(j,icell) <=0) then
!        write(*,*) "(Metal_bf) Warning at icell=", icell," T(K)=", atmos%T(icell)
!        write(*,*) metal%ID,"%n(j) density <= 0 for j=", j
!        write(*,*) "skipping this level"
       CYCLE
    end if

!      allocate(iLam(continuum%Nlambda))
!      iLam = (/ (nc, nc=continuum%Nblue, continuum%Nblue+continuum%Nlambda-1) /) !from Nblue to Nblue + Nlambda

     Z = metal%stage(i) + 1
     !! Only for Hydrogen n_eff = dsqrt(metal%g(i)/2.)
     !obtained_n = getPrincipal(metal%label(continuum%i), n_eff)

     if (.not.obtained_n) &
        n_eff = Z*dsqrt(E_RYDBERG / (metal%E(continuum%j) - metal%E(continuum%i)))

     ! for this continuum of this line
     alpha_la = 0d0

     if (.not.continuum%Hydrogenic) then
!        CALL bezier3_interp(continuum%Nlambda,continuum%lambda,continuum%alpha, &
!              continuum%Nlambda,  NLTEspec%lambda(ilam), alpha_la(ilam)
       CALL bezier3_interp(continuum%Nlambda,continuum%lambda,continuum%alpha, &
             continuum%Nlambda, &
             NLTEspec%lambda(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1), &
             alpha_la(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1))

     else

        uu(1:1) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*lambdaEdge) /(Z*Z * E_RYDBERG)
        gbf_0(1:1) = Gaunt_bf(1, uu, n_eff)
         !continuum%alpha0 is the result of the photoionisation
         !cross-section at lambdaEdge=continuum%lambda0
         ! = HPLANCK*CLIGHT/(Ej - Ei)
         ! therefore we scale wrt the lambdaEdge
         ! alpha_la(lambda=lambdaEdge)=alpha0 (containing
         ! already the gaunt factor !! and factor 1/Z^2)

!         uu(iLam) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*NLTEspec%lambda(iLam)) /(Z*Z * E_RYDBERG)
!         gbf(iLam) = Gaunt_bf(continuum%Nlambda, uu, n_eff)
!         alpha_la(iLam) = gbf(iLam) * &
!                 continuum%alpha0*((NLTEspec%lambda(iLam)/lambdaEdge)**3)*n_eff/gbf_0(1)
        uu(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) = &
         n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*NLTEspec%lambda(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1)) /(Z*Z * E_RYDBERG)
        gbf(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) = &
         Gaunt_bf(continuum%Nlambda, uu, n_eff)
        alpha_la(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) = &
         gbf(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) * &
                continuum%alpha0*((NLTEspec%lambda(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1)/lambdaEdge)**3)*n_eff/gbf_0(1)

     end if !continuum type
!      gijk(iLam) = metal%nstar(i,icell)/metal%nstar(j,icell) * expla(iLam)
!      chi(iLam) = chi(iLam) + alpha_la(iLam) * (1.-expla(iLam))*metal%n(i,icell)
!      eta(iLam) = eta(iLam) + twohnu3_c2(iLam) * gijk(iLam) * alpha_la(iLam)*metal%n(j,icell)
     gijk(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) = &
      metal%nstar(i,icell)/metal%nstar(j,icell) * expla(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1)
     chi(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) = &
      chi(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) + &
       alpha_la(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) * &
        (1.-expla(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1))*metal%n(i,icell)
     eta(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) = &
      eta(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) + &
       twohnu3_c2(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) * &
        gijk(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1) * &
         alpha_la(continuum%Nblue:continuum%Nblue+continuum%Nlambda-1)*metal%n(j,icell)
!     deallocate(iLam)
    end do ! loop over Ncont
  end do !loop over metals

 RETURN

 END FUNCTION Metal_bf


 FUNCTION Metal_bb (icell,x,y,z,x1,y1,z1,u,v,w,chi, eta, chip) result(res)
  ! Computes the emissivity and extinction of passive lines.
  ! i.e., Atoms with detailed atomic structures read but
  ! not treated in NLTE.
  ! Because damping is wavelength dependent and depend only on
  ! the grid (cell) points, here, if line%damping_initialized
  ! do not CALL Damping()
  ! the x,y,z and u,v,w quantities are used to compute the projected velocities at the
  ! cell point we are computing the opacities.
  ! Chip is only computed in Stokes transfer and contains the magneto-optical elements.
  logical :: res
  integer :: kr, m, i, j, NrecStokes, nc
  integer, intent(in) :: icell
  double precision, intent(in) :: x,y,z,u,v,w,& !positions and angles used to project
                                  x1,y1,z1      ! velocity field and magnetic field
  character(len=20) :: VoigtMethod = "HUMLICEK"
!  integer, dimension(:), allocatable :: iLam
  double precision, dimension(NLTEspec%Nwaves) :: phi, vvoigt, phiPol, phip, Vij
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: chi, eta, chip
  double precision :: twohnu3_c2, hc, fourPI, hc_4PI, gij, v0, v1, dv, dlambda
  type (AtomicLine) :: line
  type (AtomType) :: atom

  res = .true.
  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI
  
  if (NLTEspec%precomputed_voigt) VoigtMethod = "LookupTable"

  !NrecStokes = 1 in unpolarised transfer and 4 in Stokes tranfer (I,Q,U,V)
  chi = 0.  !NrecStokes per cell at eachwavelength = Nsize=NrecStokes*Nwavelength
  eta = 0.  !NrecStokes per cell at each wavelength: Nspec sized in unpolarised
  chip = 0. !NrecStokes-1 per cell at each wavelength = Nsize=(NrecStokes-1)*Nwavelength


  ! v_proj in m/s at point icell
  v0 = 0d0
!   if (atmos%moving) then !note so easy beware, check local_line_prof in optical_depth.f90
!    v0 = v_proj(icell,x,y,z,u,v,w)
!    if (lVoronoi) then ! constant velocity accross cell
!     
!     else !velocity is varying accross cell
!      v1 = v_proj(icell, x1,y1,z1,u,v,w)
!      dv = dabs(v1-v0)
!      !velocity between v0 and v1:
!      
!   end if !atmos is moving
  
  if (atmos%magnetized) then
   write(*,*) "Passive bound-bound Zeeman transitions", &
     " not implemented yet!"
   !! Magneto-optical elements (f, and r, LL04)
  end if

  do m=1,atmos%Natom ! go over all atoms
   atom = atmos%Atoms(m)
   if (atom%active) CYCLE ! go to next passive atom
    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i
     j = line%j
     
!      dlambda = line%qwing * line%lambda0 * atmos%v_char/CLIGHT
!      if ((dabs(line%lambda0 - NLTEspec%lambda(line%Nblue)) > dlambda) .or. &
!         (dabs(line%lambda0 - NLTEspec%lambda(line%Nblue+line%Nlambda-1)) > dlambda)) CYCLE

     ! now at each cell we compute it !!
     ! intialized is .true. if damping are stored on the whole grid
     line%damping_initialized=.false. !always now, because we compute damping at each icell


     if ((atom%n(j,icell) <=0).or.(atom%n(i,icell) <=0)) CYCLE !"no contrib to opac"
     ! -> prevents dividing by zero
     ! even if lcompute_atomRT(icell) it is still possible to not have a transition
     ! between the levels i and j, but not for the others.
!    if ((atom%n(j,icell) <=0).or.(atom%n(i,icell) <=0)) then !no transition
!        write(*,*) "(Metal_bb) Warning at icell=", icell," T(K)=", atmos%T(icell)
!        write(*,*) atom%ID," density <= 0 ", i, j, line%lambda0, atom%n(i,icell), atom%n(j,icell)
!        write(*,*) "skipping this level"
!       CYCLE
!     end if

     ! using directly line%Nblue:line%Nblue+line%Nlambda-1 instead of iLam
     ! could save time, because for each line we allocate/deallocate iLam.
     ! just for a more readable code, which is not negligible at all.
!      allocate(iLam(line%Nlambda))
!      iLam = (/ (nc, nc=line%Nblue, line%Nblue+line%Nlambda-1) /)

     phi = 0d0
     phip = 0d0
     phiPol = 0d0

     gij = line%Bji / line%Bij
     twohnu3_c2 = line%Aji / line%Bji
    !convert to doppler units
    
    ! init for this line of this atom
!      vvoigt(iLam) = (NLTEspec%lambda(iLam)-line%lambda0) * &
!          CLIGHT / (line%lambda0 * atom%vbroad(icell))
     vvoigt(line%Nblue:line%Nblue+line%Nlambda-1) = &
      (NLTEspec%lambda(line%Nblue:line%Nblue+line%Nlambda-1)-line%lambda0) * &
         CLIGHT / (line%lambda0 * atom%vbroad(icell))
     if (line%voigt) then
      CALL Damping(icell, atom, kr, line%adamp)
      
!        vvoigt(iLam) = vvoigt(iLam) - v0 / atom%vbroad(icell) !add velocity field
!        phi(iLam) = Voigt(line%Nlambda, line%adamp, vvoigt(iLam), &
!                 phip, VoigtMethod)
       vvoigt(line%Nblue:line%Nblue+line%Nlambda-1) = &
        vvoigt(line%Nblue:line%Nblue+line%Nlambda-1) - v0 / atom%vbroad(icell) !add velocity field
       phi(line%Nblue:line%Nblue+line%Nlambda-1) = Voigt(line%Nlambda, line%adamp, &
          vvoigt(line%Nblue:line%Nblue+line%Nlambda-1), phip, VoigtMethod)
                !phip is on the whole grid, you take indexes after. Otherwise the function
                !Voigt will return an error
!       phiPol(iLam) = phip(iLam)
!       phiPol(line%Nblue:line%Nblue+line%Nlambda-1) = phip(line%Nblue:line%Nblue+line%Nlambda-1)
     else !Gaussian
!        vvoigt(iLam) = vvoigt(iLam) - v0 / atom%vbroad(icell)
       vvoigt(line%Nblue:line%Nblue+line%Nlambda-1) = &
        vvoigt(line%Nblue:line%Nblue+line%Nlambda-1) - v0 / atom%vbroad(icell)
!        phi(iLam) = dexp(-(vvoigt(iLam))**2)
       phi(line%Nblue:line%Nblue+line%Nlambda-1) = dexp(-(vvoigt(line%Nblue:line%Nblue+line%Nlambda-1))**2)
     end if
     !Sum up all contributions for this line with the other
!      Vij(iLam) = hc_4PI * line%Bij * phi(iLam) / (SQRTPI * atom%vbroad(icell))
     Vij(line%Nblue:line%Nblue+line%Nlambda-1) = &
      hc_4PI * line%Bij * phi(line%Nblue:line%Nblue+line%Nlambda-1) / (SQRTPI * atom%vbroad(icell))
!      chi(iLam) = chi(iLam) + Vij(iLam) * (atom%n(i,icell)-gij*atom%n(j,icell))
     chi(line%Nblue:line%Nblue+line%Nlambda-1) = &
      chi(line%Nblue:line%Nblue+line%Nlambda-1) + &
       Vij(line%Nblue:line%Nblue+line%Nlambda-1) * (atom%n(i,icell)-gij*atom%n(j,icell))
!      eta(iLam) = eta(iLam) + twohnu3_c2 * gij * Vij(iLam) * atom%n(j,icell)
     eta(line%Nblue:line%Nblue+line%Nlambda-1) = &
      eta(line%Nblue:line%Nblue+line%Nlambda-1) + &
       twohnu3_c2 * gij * Vij(line%Nblue:line%Nblue+line%Nlambda-1) * atom%n(j,icell)
       
!      write(*,*) atom%ID, line%j,"->",line%i, NLTEspec%lambda(line%NBlue+(line%Nlambda-1)/2) - line%lambda0, &
!       " nspect=", line%NBlue+(line%Nlambda-1)/2 - 1
!      write(*,*) Vij(line%NBlue+(line%Nlambda-1)/2), &
!         phi(line%Nblue+(line%Nlambda-1)/2), gij, line%adamp
       
     !dealloc indexes for next line
!     deallocate(iLam)
    end do !end loop on lines for this atom
  end do !end loop over Natom

 RETURN

 END FUNCTION Metal_bb

 SUBROUTINE Background(id,icell,x,y,z,x1,y1,z1,u,v,w)
  integer, intent(in) :: icell, id
  double precision, intent(in) :: x, y, z, u, v, w, &
                                  x1, y1, z1!only relevant for b-b when vector fields are present
  double precision, dimension(NLTEspec%Nwaves) :: chi, eta, sca, Bpnu, chip
  integer :: Npassive

  if (.not.atmos%lcompute_atomRT(icell)) RETURN !nH <= tiny_nH or T <= tiny_T == empty cell
  ! all opac are zero, return.

!   if ((atmos%nHtot(icell)==0d0).or.(atmos%T(icell)==0d0)) &
!     RETURN ! stoping for this cell,
            ! it is free of (significant) gas
            ! so no emission/absorption. Coefficients set to 0d0 for all wavelengths
  !Do not forget that it is still possible however, that lcompute_atomRT is TRUE,
  !but that the temperature is too low to have all the levels non zero.
  !The compute_atomRT ensures that at least one level has non-zero populations,
  !to avoid division by zero.

  CALL Bplanck(atmos%T(icell), NLTEspec%Nwaves, NLTEspec%lambda, Bpnu)

   chi = 0.
   eta = 0.
   sca = 0.
   chip = 0.


   if (Rayleigh(icell, Hydrogen, sca)) NLTEspec%AtomOpac%sca_c(id,:) = sca
   NLTEspec%AtomOpac%sca_c(id,:) =  NLTEspec%AtomOpac%sca_c(id,:) + Thomson(icell)
   if (associated(Helium)) then
    if (Rayleigh(icell, Helium, sca)) NLTEspec%AtomOpac%sca_c(id,:) = &
          NLTEspec%AtomOpac%sca_c(id,:) + sca
   end if

   NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%sca_c(id,:)

   CALL Hydrogen_ff(icell, chi)
   NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
   NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + chi * Bpnu

   if (Hminus_bf(icell, chi,eta)) then
    NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
    NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + eta
   end if

   CALL Hminus_ff(icell, chi)
    NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
    NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + chi * Bpnu

   if (.not.Hydrogen%active) then !passive bound-free
    if (Hydrogen_bf(icell, chi, eta)) then
     NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
     NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + eta
    end if
   end if

   Npassive = atmos%Natom - atmos%Nactiveatoms
   if (Npassive == 0) RETURN !no passive bound-bound and bound-free
   ! remove H if passive because it is computed in Hydrogen_bf()
   if (.not.Hydrogen%active) Npassive = Npassive - 1
   if (Npassive>0) then !other metal passive ?
    if (Metal_bf(icell, chi, eta)) then
     NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
     NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + eta
    end if
   end if


   !keep pure continuum opacities now
   NLTEspec%AtomOpac%eta_c(id,:) = NLTEspec%AtomOpac%eta_p(id,:)
   NLTEspec%AtomOpac%chi_c(id,:) = NLTEspec%AtomOpac%chi_p(id,:)

   Npassive = atmos%Natom - atmos%Nactiveatoms !including Hydrogen
   if (Npassive > 0) then !go into it even if only H is passive
    if (Metal_bb(icell, x, y, z, x1, y1, z1, u, v, w, chi, eta, chip)) then
     NLTEspec%AtomOpac%chi_p(id,:) = NLTEspec%AtomOpac%chi_p(id,:) + chi
     NLTEspec%AtomOpac%eta_p(id,:) = NLTEspec%AtomOpac%eta_p(id,:) + eta
    end if
   end if

 RETURN
 END SUBROUTINE Background

END MODULE metal
