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

 use grid_type, only : atmos, Hydrogen, Helium
 use constant
 use math, only : CUBE, SQ,  bezier3_interp, CUBEarr, SQarr, interp1D
 use atom_type
 use spectrum_type, only : NLTEspec, ActiveSetType, initAS
 use hydrogen_opacities
 use voigtfunctions, only : Voigt
 use broad, only : Damping
 use thomson_scattering
 use Rayleigh_scattering
 use Planck

 IMPLICIT NONE

 integer, parameter :: N_MAX_OVERLAP=10 !transitions that overlap
                           ! at a given wavelength

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
  logical :: res, hunt, obtained_n
  integer :: m, k, kr, i, j, Z, la
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
       write(*,*) "Error, negative LTE pops of atom ",metal%ID," exiting..."
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
   write(*,*) "Error in Metal_bf()"
   write(*,*) "Ignore the previous message if no passive atoms yet"
   RETURN
  end if
 !RETURN
 END FUNCTION

 FUNCTION Metal_bb (icell, chi, eta, chip) result(res)
  ! Computes the emissivity and extinction of passive lines.
  ! i.e., Atoms with detailed atomic structures read but
  ! not treated in NLTE.
  ! Because damping is wavelength dependent and depend only on
  ! the grid (cell) points, here, if line%damping_initialized
  ! do not CALL Damping()
  logical :: res
  integer :: k, kr, l, m, nc, i, j, Nll, la
  ! Nll = number of lines at this wavelength
  integer, intent(in) :: icell
  double precision :: dlambda, phi, v, phiPol
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: chi, eta, chip
  double precision :: twohnu3_c2, hc, fourPI, hc_4PI, gij, Vij
  double precision, allocatable :: n(:)
  type (AtomicLine) :: line
  type (AtomType) :: atom

  res = .false.
  hc = HPLANCK * CLIGHT
  fourPI = 4.*PI
  hc_4PI = hc/fourPI
  Nll = 0
 
  chi = 0.
  eta = 0.
  chip = 0.


  if (atmos%magnetized) then
   write(*,*) "Passive bound-bound Zeeman transitions", &
     " not implemented yet!"
   !do k=atmos%Nspace,4*atmos%Nspace
    !chi(k) = 0.
    !eta(k) = 0.
   !end do
   !! Magneto-optical elements (f, and r, LL04)
   !do k=1,atmos%Nspace*3
    !chip(k) = 0.
   !end do
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
           atmos%vmicro_char/CLIGHT

     !dlambda is the line domain. Extending up to qwing*l0/c
     !corrected by the typical micro turbulent veloctity.
     !actually, vmicro_char could be the largest velocity

     ! only count this line, if lambda falls in the
     ! line domain
     ! count how many transitions at a wavelenght point
     ! if (Nll.ge.N_MAX_OVERLAP) then
     !     write(*,*) "Too many overlapping transitions"
     !     stop
     ! end if
    gij = line%Bji / line%Bij
    twohnu3_c2 = line%Aji / line%Bji
    do la=1,NLTEspec%Nwaves 
     if (dabs(NLTEspec%lambda(la)-line%lambda0).le.dlambda) then
        res = .true. !set the flag to true for this wavelength

        ! calculate line's damping parameter if not initialized
        if ((line%Voigt).and.&
           (.not.line%damping_initialized)) then
           CALL Damping(icell, atom, kr, line%adamp)
           !write(*,*) icell, line%adamp, line%damping_initialized
           line%damping_initialized=.true.
           ! again, adamp is only depend on the grid point
           ! not of the lambda at the opacity is calculated
           ! so as we loop over all line of an atom, we do not
           ! need to recompute it every time.
        end if
      !Now absorption and emission coefficients
       do nc=1,line%Ncomponent !in case of multi-component line
         v = (NLTEspec%lambda(la)-line%lambda0-line%c_shift(nc))
         !convert to doppler velocity
         v = v * CLIGHT / (line%lambda0 * atom%vbroad(icell))
         v = v !-vproject(icell) / atom%vbroad(icell)
         if (line%Voigt) then
          ! A trick is used in Voigt() to avoid v=-v
          ! out of the function VoigtArmstrong if v<0.
          ! But v is unimportant because it is a temp var.
          !write(*,* ) "vbroad=",atom%vbroad(k),"a=",line%adamp(k)
          phi = Voigt(line%adamp, v, &
               phiPol, 'ARMSTRONG ') * line%c_fraction(nc)
          ! if magnetized and line%polarizable etc ...
         else
          phi = dexp(-SQ(v))
         end if
         Vij = hc_4PI * line%Bij * phi / (SQRTPI * atom%vbroad(icell))
         chi(la) = chi(la) + Vij * (n(i)-gij*n(j))
         eta(la) = eta(la) + twohnu3_c2 * gij * Vij * n(j)
       end do
      end if !lambda is in the line region
     end do !over wavelength
    end do !end loop on lines for this atom
    deallocate(n)
   end if !atom is passive
  end do !end loop over Natom
  
  if ((MAXVAL(chi).gt.0.).and.(MAXVAL(eta).gt.0)) then
   res = .true. !should alway be true !
   RETURN
  else
   res = .false.
   write(*,*) "Error in Metal_bb()"
   RETURN
  end if
 END FUNCTION Metal_bb
 
 SUBROUTINE Background(icell)
  integer, intent(in) :: icell
  double precision, dimension(NLTEspec%Nwaves) :: chi, eta, sca, Bpnu, chip
  
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
   
   if (Hydrogen_bf(icell, chi, eta)) then
    NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
    NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
   end if

   if (Metal_bf(icell, chi, eta)) then
    NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi 
    NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
   end if 
   
   if (Metal_bb(icell, chi, eta, chip)) then
    NLTEspec%ActiveSet%chi_c = NLTEspec%ActiveSet%chi_c + chi
    NLTEspec%ActiveSet%eta_c = NLTEspec%ActiveSet%eta_c + eta
   end if
   

 
 RETURN
 END SUBROUTINE Background

END MODULE metal
