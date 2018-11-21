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
 use math, only : CUBE, SQ,  bezier3_interp, interp1D, interp1Darr, SQarr, CUBEarr
 use atom_type
 use spectrum_type, only : NLTEspec, ActiveSetType, initAS
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
 FUNCTION back_Metal_bf(icell, chi, eta) result(res)
 !cross-section in cm2 per particle is given by Kramers’ formula
  !with n the principal quantum number of the level i from
 !which the atom or ion is ionized, Z the ion charge, ν in Hz and gbf the dimensionless
 !Gaunt factor, a quantummechanical correction factor of order unity.
 !The Kramers cross-section decays ∼ ν−3 above the threshold (“edge”) frequency ν0,
 !being zero below it because the threshold energy is the required minimum. Think the
  ! inverse in terms of wavelengths
  integer :: icell
  logical :: res, hunt, obtained_n
  integer :: m, k, kr, i, j, Z, la, Npassive
  type (AtomType) :: metal
  type (AtomicContinuum) :: continuum
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: eta, chi
  double precision :: lambdaEdge, twohc, hc_k, n_eff
  double precision, allocatable, dimension(:) :: n
  double precision, dimension(1) :: gbf_0, uu, gbf
  double precision :: twohnu3_c2, gijk, alpha_la, hc_kla, expla
  res = .false.

  twohc = (2. * HPLANCK * CLIGHT) / CUBE (NM_TO_M)
  hc_k = (HPLANCK * CLIGHT) / (KBOLTZMANN * NM_TO_M)

   chi = 0.
   eta = 0.

   Npassive = atmos%Natom - atmos%Nactiveatoms
   ! remove H if passive because it is computed in Hydrogen_bf()
   if (.not.Hydrogen%active) Npassive = Npassive - 1
   !write(*,*) "Npassive metals : ", Npassive
   if (Npassive <= 0) RETURN ! no passive metals

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
   allocate(n(metal%Nlevel))
   if (.not.metal%active) then

   !if loc(%n)=loc(%nstar) pure ETL, or if associated(%n,%nstar)
   !else NLTEpops read from previous run for actual passive metal
    n = metal%n(:,icell) !nstar if pure LTE, %n if NLTEpops.eq.true.
   !write(*,*) loc(metal%n), loc(metal%nstar), &
   !                  associated(metal%n,metal%nstar)


    uu = 0.
    do kr=1,metal%Ncont
     continuum = metal%continua(kr)
     i = continuum%i
     j = continuum%j !+1 wrt C indexing
     lambdaEdge = continuum%lambda0! or ionisation wavelength or wavelength
               ! associated to the minimal frquency needed
               ! to unbound an electron

     if (metal%nstar(j,icell).le.0.) then
       write(*,*) "Error, negative or null LTE pops of atom ",metal%ID," exiting..."
       stop
     end if

     !internal wavelength loop
     !cannot use easily where statement which is normally dedicated to this condition
     do la=1,NLTEspec%Nwaves
     !!where((NLTEspec%lambda(la).le.lambdaEdge).and. &
     !!     (NLTEspec%lambda(la).ge.continuum%lambda(1))) -> but inside 2 do loops
     !! and an if statement and contais if statements
      if ((NLTEspec%lambda(la).le.lambdaEdge).and. &
          (NLTEspec%lambda(la).ge.continuum%lambda(1))) then
       hc_kla = hc_k/NLTEspec%lambda(la) !factor 1/NM_TO_M in hc_k
       twohnu3_c2 = twohc / CUBE(NLTEspec%lambda(la))
       expla = dexp(-hc_kla/atmos%T(icell))

        if (continuum%Hydrogenic) then
         Z = metal%stage(i) + 1
         !! Only for Hydrogen n_eff = dsqrt(metal%g(i)/2.)
         obtained_n = getPrincipal(metal%label(continuum%i), n_eff)
          if (.not.obtained_n) &
            n_eff = Z*dsqrt(E_RYDBERG / (metal%E(continuum%j) - metal%E(continuum%i)))
          uu(1) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*lambdaEdge) &
               /(Z*Z * E_RYDBERG)
          gbf_0 = Gaunt_bf(1, uu, n_eff)
         !continuum%alpha0 is the result of the photoionisation
         !cross-section at lambdaEdge=continuum%lambda0
         ! = HPLANCK*CLIGHT/(Ej - Ei)
         ! therefore we scale wrt the lambdaEdge
         ! alpha_la(lambda=lambdaEdge)=alpha0 (containing
         ! already the gaunt factor !! and factor 1/Z^2)
         uu(1) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*NLTEspec%lambda(la)) &
               /(Z*Z * E_RYDBERG)
         gbf = Gaunt_bf(1, uu, n_eff)
         alpha_la = gbf(1) * &
                continuum%alpha0*CUBE(NLTEspec%lambda(la)/lambdaEdge)*n_eff/gbf_0(1)
        else ! use interpolation on data
         alpha_la = interp1D(continuum%lambda,continuum%alpha,NLTEspec%lambda(la))
        end if ! over continuum type


        gijk = metal%nstar(i,icell)/metal%nstar(j,icell) * expla
        chi(la) = chi(la) + alpha_la * (1.-expla)*n(i)
        eta(la) = eta(la) + twohnu3_c2 * gijk * alpha_la*n(j)
        !if (NLTEspec%lambda(la).ge.405) then
        ! write(*,*) kr, lambdaEdge, metal%ID, NLTEspec%lambda(la), alpha_la,
        !   icell, gijk, chi(la), eta(la)
        !end if
      end if ! loop if lambda is in the good range
      ! else chi(la)=eta(la)=0.0
     end do ! over wavelength la
    end do ! loop over Ncont

   end if !loop if not active
   deallocate(n)
  end do !loop over metals

  if ((MAXVAL(chi).gt.0.).and.(MAXVAL(eta).gt.0)) then
   res = .true. !should alway be true if at least one passive atom!
   RETURN
  else
   res = .false.
   write(*,*) "Error in Metal_bf()" !error if false but there are passive atoms
   RETURN
  end if
 !RETURN
 END FUNCTION  back_Metal_bf
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
  integer, dimension(:), allocatable :: iLam
  type (AtomType) :: metal
  type (AtomicContinuum) :: continuum
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: eta, chi
  double precision :: lambdaEdge, twohc, hc_k, n_eff
  double precision, dimension(NLTEspec%Nwaves) :: gbf_0, uu, gbf
  double precision, dimension(NLTEspec%Nwaves) :: twohnu3_c2, gijk, hc_kla, expla, alpha_la

  res = .true.
  obtained_n = .false.

  twohc = (2. * HPLANCK * CLIGHT) / CUBE (NM_TO_M)
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

     allocate(iLam(continuum%Nlambda))
     iLam = (/ (nc, nc=continuum%Nblue, continuum%Nblue+continuum%Nlambda-1) /) !from Nblue to Nblue + Nlambda

     Z = metal%stage(i) + 1
     !! Only for Hydrogen n_eff = dsqrt(metal%g(i)/2.)
     !obtained_n = getPrincipal(metal%label(continuum%i), n_eff)

     if (.not.obtained_n) &
        n_eff = Z*dsqrt(E_RYDBERG / (metal%E(continuum%j) - metal%E(continuum%i)))

     ! for this continuum of this line
     alpha_la = 0d0

     if (.not.continuum%Hydrogenic) then

!      if (.not.continuum%Hydrogenic) then
!       where((NLTEspec%lambda.le.lambdaEdge).and. &
!             (NLTEspec%lambda.ge.continuum%lambda(1)))
!         alpha_la = interp1Darr(continuum%lambda, continuum%alpha, NLTEspec%lambda)
!       end where
       CALL bezier3_interp(continuum%Nlambda,continuum%lambda,continuum%alpha, &
             continuum%Nlambda, NLTEspec%lambda(iLam), alpha_la(iLam))
!
!         where((NLTEspec%lambda.gt.lambdaEdge).and.&
!               (NLTEspec%lambda.lt.continuum%lambda(1)))
!           alpha_la = 0d0
!         end where
     else
!       where((NLTEspec%lambda.le.lambdaEdge).and. &
!             (NLTEspec%lambda.ge.continuum%lambda(1)))

!         hc_kla = hc_k/NLTEspec%lambda !factor 1/NM_TO_M in hc_k
!         twohnu3_c2 = twohc / NLTEspec%lambda**3
!         expla = dexp(-hc_kla/atmos%T(icell))

        uu(1:1) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*lambdaEdge) /(Z*Z * E_RYDBERG)
        gbf_0(1:1) = Gaunt_bf(1, uu, n_eff)
         !continuum%alpha0 is the result of the photoionisation
         !cross-section at lambdaEdge=continuum%lambda0
         ! = HPLANCK*CLIGHT/(Ej - Ei)
         ! therefore we scale wrt the lambdaEdge
         ! alpha_la(lambda=lambdaEdge)=alpha0 (containing
         ! already the gaunt factor !! and factor 1/Z^2)
!         uu = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*NLTEspec%lambda) /(Z*Z * E_RYDBERG)
!         gbf = Gaunt_bf(NLTEspec%Nwaves, uu, n_eff)
!         alpha_la = gbf * &
!                 continuum%alpha0*((NLTEspec%lambda/lambdaEdge)**3)*n_eff/gbf_0(1)
        uu(iLam) = n_eff*n_eff*HPLANCK*CLIGHT/(NM_TO_M*NLTEspec%lambda(iLam)) /(Z*Z * E_RYDBERG)
        gbf(iLam) = Gaunt_bf(continuum%Nlambda, uu, n_eff)
        alpha_la(iLam) = gbf(iLam) * &
                continuum%alpha0*((NLTEspec%lambda(iLam)/lambdaEdge)**3)*n_eff/gbf_0(1)
!       end where
     end if !continuum type
!      gijk = metal%nstar(i,icell)/metal%nstar(j,icell) * expla
!      chi = chi + alpha_la * (1.-expla)*metal%n(i,icell)
!      eta = eta + twohnu3_c2 * gijk * alpha_la*metal%n(j,icell)
     gijk(iLam) = metal%nstar(i,icell)/metal%nstar(j,icell) * expla(iLam)
     chi(iLam) = chi(iLam) + alpha_la(iLam) * (1.-expla(iLam))*metal%n(i,icell)
     eta(iLam) = eta(iLam) + twohnu3_c2(iLam) * gijk(iLam) * alpha_la(iLam)*metal%n(j,icell)
     !write(*,*) MAXVAL(chi), MAXVAL(eta), metal%ID
     deallocate(iLam)
    end do ! loop over Ncont
  end do !loop over metals

! Condition tested outside  so should be at least one passive atom if we enter
!   if ((MAXVAL(chi).gt.0.).and.(MAXVAL(eta).gt.0) &
!       .and.(atmos%Natom - atmos%Nactiveatoms)) then
!    res = .true. !should alway be true if at least one passive atom!
!    RETURN
!   else
!    res = .false.
!    write(*,*) "Error in Metal_bf()" !error if false but there are passive atoms
!    RETURN
!   end if
 RETURN

 END FUNCTION Metal_bf

 FUNCTION back_Metal_bb (icell,x,y,z,u,v,w,chi, eta, chip) result(res)
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
  integer :: k, kr, l, m, nc, i, j, la, NrecStokes
  integer, intent(in) :: icell
  double precision, intent(in) :: x,y,z,u,v,w !positions and angles used to project
                                              ! velocity field and magnetic field
  double precision :: dlambda, phi(1), vvoigt(1), phiPol(1)
  ! in scalar version of Voigt, phi is of size 1
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: chi, eta, chip
  double precision :: twohnu3_c2, hc, fourPI, hc_4PI, gij, Vij
  double precision, allocatable :: n(:)
  type (AtomicLine) :: line
  type (AtomType) :: atom

  res = .false.
  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI


  !NrecStokes = 1 in unpolarised transfer and 4 in Stokes tranfer (I,Q,U,V)
  chi = 0.  !NrecStokes per cell at eachwavelength = Nsize=NrecStokes*Nwavelength
  eta = 0.  !NrecStokes per cell at each wavelength: Nspec sized in unpolarised
  chip = 0. !NrecStokes-1 per cell at each wavelength = Nsize=NrecStokes*Nwavelength


  if (atmos%magnetized) then
   write(*,*) "Passive bound-bound Zeeman transitions", &
     " not implemented yet!"

   !! Magneto-optical elements (f, and r, LL04)
  end if

  do m=1,atmos%Natom ! go over all atoms
   atom = atmos%Atoms(m)
   allocate(n(atom%Nlevel))
   if (.not.atom%active) then ! but only passive atoms
     n = atom%n(:,icell) !if NLTE pops exist for this atom %n != %nstar

    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i
     j = line%j
     dlambda = line%lambda0 * line%qwing * &
           atmos%v_char/CLIGHT
     ! calculate line's damping parameter at this icell
     ! depends only on icell, not lam
     line%damping_initialized=.true.
     if (line%Voigt) CALL Damping(icell, atom, kr, line%adamp)

     !dlambda is the line domain. Extending up to qwing*l0/c
     !corrected by the typical (micro-turbulent) veloctity.

     ! only count this line, if lambda falls in the
     ! line domain
    gij = line%Bji / line%Bij
    twohnu3_c2 = line%Aji / line%Bji
    do la=1,NLTEspec%Nwaves
     if (dabs(NLTEspec%lambda(la)-line%lambda0).le.dlambda) then
        res = .true. !set the flag to true for this wavelength

!         !calculate line's damping parameter if not initialized
!         if ((line%Voigt).and.&
!            (.not.line%damping_initialized)) then
!            CALL Damping(icell, atom, kr, line%adamp)
!            !write(*,*) icell, line%adamp, line%damping_initialized
!            !line%damping_initialized=.true.
!            !again, adamp is only depend on the grid point
!            !not of the lambda at the opacity is calculated
!            !so as we loop over all line of an atom, we do not
!            !need to recompute it every time.
!         end if
      !Now absorption and emission coefficients
       do nc=1,line%Ncomponent !in case of multi-component line
         vvoigt(1) = (NLTEspec%lambda(la)-line%lambda0-line%c_shift(nc))
         !convert to doppler velocity
         vvoigt(1) = vvoigt(1) * CLIGHT / (line%lambda0 * atom%vbroad(icell))
         ! the position and the angle are used for projection
         vvoigt(1) = vvoigt(1) - 0d0 / atom%vbroad(icell)!v_proj(icell,x,y,z,u,v,w)
         if (line%Voigt) then
          ! A trick is used in Voigt() to avoid v=-v
          ! out of the function VoigtArmstrong if v<0.
          ! But v is unimportant because it is a temp var.
          !write(*,* ) "vbroad=",atom%vbroad(icell),"a=",line%adamp
          phi = Voigt(1, line%adamp, vvoigt, &
               phiPol, 'ARMSTRONG ') * line%c_fraction(nc)
          ! if magnetized and line%polarizable etc ...
         else
          phi(1) = dexp(-SQ(vvoigt(1))) * line%c_fraction(nc)
         end if
         Vij = hc_4PI * line%Bij * phi(1) / (SQRTPI * atom%vbroad(icell))
         chi(la) = chi(la) + Vij * (n(i)-gij*n(j))
         eta(la) = eta(la) + twohnu3_c2 * gij * Vij * n(j)
       end do
      end if !lambda is in the line region
     line%damping_initialized=.false. !for next cell point for this atom
     end do !over wavelength
    end do !end loop on lines for this atom
    deallocate(n)
   end if !atom is passive
  end do !end loop over Natom

  if ((MAXVAL(chi).gt.0.).and.(MAXVAL(eta).gt.0)) then
   res = .true. !should alway be true , because at least always H!
   RETURN
  else
   res = .false.
   write(*,*) "Error in Metal_bb()"
   RETURN
  end if
 END FUNCTION back_Metal_bb

 FUNCTION Metal_bb (icell,x,y,z,u,v,w,chi, eta, chip) result(res)
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
  integer :: kr, m, nc, i, j, NrecStokes
  integer, intent(in) :: icell
  double precision, intent(in) :: x,y,z,u,v,w !positions and angles used to project
                                              ! velocity field and magnetic field
  double precision :: dlambda
  integer, dimension(:), allocatable :: iLam
  double precision, dimension(NLTEspec%Nwaves) :: phi, vvoigt, phiPol, phip, Vij
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: chi, eta, chip
  double precision :: twohnu3_c2, hc, fourPI, hc_4PI, gij, v0
  type (AtomicLine) :: line
  type (AtomType) :: atom

  res = .true.
  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI


  !NrecStokes = 1 in unpolarised transfer and 4 in Stokes tranfer (I,Q,U,V)
  chi = 0.  !NrecStokes per cell at eachwavelength = Nsize=NrecStokes*Nwavelength
  eta = 0.  !NrecStokes per cell at each wavelength: Nspec sized in unpolarised
  chip = 0. !NrecStokes-1 per cell at each wavelength = Nsize=(NrecStokes-1)*Nwavelength


  ! v_proj in m/s at point icell
  v0 = 0d0
  !if (atmos%moving) v0 = v_proj(icell,x,y,z,u,v,w)

  if (atmos%magnetized) then
   write(*,*) "Passive bound-bound Zeeman transitions", &
     " not implemented yet!"
   !! Magneto-optical elements (f, and r, LL04)
  end if

  do m=1,atmos%Natom ! go over all atoms
   atom = atmos%Atoms(m)
   if (atom%active) CYCLE ! go to next passive atom
    ! I use atom%n(level,icell) !if NLTE pops exist for this atom, atom%n != atom%nstar
    ! else atom%n = atom%nstar.

    do kr=1,atom%Nline ! for this atom go over all transitions
                       ! bound-bound
     line = atom%lines(kr)
     i = line%i
     j = line%j
     dlambda = line%lambda0 * line%qwing * &
           atmos%v_char/CLIGHT

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
     allocate(iLam(line%Nlambda))
     iLam = (/ (nc, nc=line%Nblue, line%Nblue+line%Nlambda-1) /)
     !from Nblue to Nblue + Nlambda !try with line%Nblue:line%Nblue+line%Nlambda-1
                                      !instead of iLam if not working
     phi = 0d0
     phip = 0d0
     phiPol = 0d0

    ! now at each cell we compute it !!
     line%damping_initialized=.false.
     if (line%Voigt) CALL Damping(icell, atom, kr, line%adamp)

     !dlambda is the line domain. Extending up to qwing*l0/c
     !corrected by the typical (micro-turbulent) veloctity.

     ! only count this line, if lambda falls in the
     ! line domain
     gij = line%Bji / line%Bij
     twohnu3_c2 = line%Aji / line%Bji

    !convert to doppler units
    ! init for this line of this atom
     vvoigt(iLam) = (NLTEspec%lambda(iLam)-line%lambda0) * &
         CLIGHT / (line%lambda0 * atom%vbroad(icell))

     if (line%voigt) then
      do nc=1,line%Ncomponent !add multicomponent if any, else Ncomponent=1
       vvoigt(iLam) = vvoigt(iLam) - line%c_shift(nc) * &
             CLIGHT / (line%lambda0 * atom%vbroad(icell)) - v0 / atom%vbroad(icell) !add velocity field
       phi(iLam) = phi(iLam) + Voigt(line%Nlambda, line%adamp, vvoigt(iLam), &
                phip, 'ARMSTRONG ') * line%c_fraction(nc)
                !phip is on the whole grid, you take indexes after. Otherwise the function
                !Voigt will return an error
       phiPol(iLam) = phiPol(iLam) + phip(iLam)
      end do
     else !Gaussian
      do nc=1,line%Ncomponent
       vvoigt(iLam) = vvoigt(iLam) - line%c_shift(nc) * &
             CLIGHT / (line%lambda0 * atom%vbroad(icell)) - v0 / atom%vbroad(icell)
       phi(iLam) = phi(iLam) + dexp(-(vvoigt(iLam))**2) * line%c_fraction(nc)
      end do
     end if


!     if (line%voigt) then
!      do nc=1,line%Ncomponent !in case of multi-component line
!       where(dabs(NLTEspec%lambda-line%lambda0).le.dlambda)
!        phi = phi + Voigt(NLTEspec%Nwaves, line%adamp, vvoigt - &
!           line%c_shift(nc) * CLIGHT / (line%lambda0 * atom%vbroad(icell)), &
!                phip, 'ARMSTRONG ') * line%c_fraction(nc)
!        phiPol = phiPol + phip
!           ! if magnetized and line%polarizable etc ...
!       end where
!      end do
!     else !Gaussian
!      do nc=1,line%Ncomponent !
!       where(dabs(NLTEspec%lambda-line%lambda0).le.dlambda)
!        phi = phi + exp(-(vvoigt-line%c_shift(nc)* &
!                CLIGHT / (line%lambda0 * atom%vbroad(icell)))**2) * line%c_fraction(nc)
!        !phiPol = phiPol + phip !force weak field in that case
!       end where
!      end do
!     end if !type of profile

     !Sum up all contributions for this line with the other
     Vij(iLam) = hc_4PI * line%Bij * phi(iLam) / (SQRTPI * atom%vbroad(icell))
     chi(iLam) = chi(iLam) + Vij(iLam) * (atom%n(i,icell)-gij*atom%n(j,icell))
     eta(iLam) = eta(iLam) + twohnu3_c2 * gij * Vij(iLam) * atom%n(j,icell)
!      Vij = hc_4PI * line%Bij * phi / (SQRTPI * atom%vbroad(icell))
!      chi = chi + Vij * (atom%n(i,icell)-gij*atom%n(j,icell))
!      eta = eta + twohnu3_c2 * gij * Vij * atom%n(j,icell)
     !dealloc indexes for next line
     deallocate(iLam)
    end do !end loop on lines for this atom
  end do !end loop over Natom


! Condition tested outside    so should alxays be a passive atom if we enter
!   if ((MAXVAL(chi).gt.0.).and.(MAXVAL(eta).gt.0) &
!                .and.(atmos%Natom - atmos%Nactiveatoms)) then
!    res = .true. !should alway be true if passive atoms
!    RETURN
!   else
!    res = .false.
!    write(*,*) "Error in Metal_bb()"
!    RETURN
!   end if
 RETURN

 END FUNCTION Metal_bb

 SUBROUTINE Background(icell,x,y,z,u,v,w)
  integer, intent(in) :: icell
  double precision, intent(in) :: &
              x, y, z, u, v, w!only relevant for b-b when vector fields are present
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


   if (Rayleigh(icell, Hydrogen, sca)) NLTEspec%ActiveSet%sca_c = sca
   NLTEspec%ActiveSet%sca_c =  NLTEspec%ActiveSet%sca_c + Thomson(icell)
   if (associated(Helium)) then
    if (Rayleigh(icell, Helium, sca)) NLTEspec%ActiveSet%sca_c = &
          NLTEspec%ActiveSet%sca_c + sca
   end if

   NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%sca_c

   CALL Hydrogen_ff(icell, chi)
   NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
   NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + chi * Bpnu

   if (Hminus_bf(icell, chi,eta)) then
    NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
    NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
   end if

   CALL Hminus_ff(icell, chi)
    NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
    NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + chi * Bpnu

   if (.not.Hydrogen%active) then !passive bound-free
    if (Hydrogen_bf(icell, chi, eta)) then
     NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
     NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
    end if
   end if

   Npassive = atmos%Natom - atmos%Nactiveatoms
   if (Npassive == 0) RETURN !no passive bound-bound and bound-free
   ! remove H if passive because it is computed in Hydrogen_bf()
   if (.not.Hydrogen%active) Npassive = Npassive - 1
   if (Npassive>0) then !other metal passive ?
    if (Metal_bf(icell, chi, eta)) then
     NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
     NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
    end if
   end if


   !keep pure continuum opacities now
   NLTEspec%ActiveSet%eta_c_bf = NLTEspec%ActiveSet%eta_c
   NLTEspec%ActiveSet%chi_c_bf = NLTEspec%ActiveSet%chi_c

   Npassive = atmos%Natom - atmos%Nactiveatoms !including Hydrogen
   if (Npassive > 0) then !go into it even if only H is passive
    if (Metal_bb(icell, x, y, z, u, v, w, chi, eta, chip)) then
     NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
     NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
    end if
   end if


 RETURN
 END SUBROUTINE Background

END MODULE metal
