MODULE rayleigh_scattering

 use atom_type, only       : AtomType
 use atmos_type, only      : atmos, Hydrogen, Helium
 use spectrum_type, only   : NLTEspec
 use constant
 use parametres

 IMPLICIT NONE

 double precision, parameter :: LONG_RAYLEIGH_WAVE=1.d6 !nm

 CONTAINS
 
 SUBROUTINE HI_Rayleigh(id, icell)
 ! ------------------------------------------------------------- !
  ! Rayleigh scattering on neutral Hydrogen.
  ! Baschek & Scholz 1982

 ! ------------------------------------------------------------- !
  !type (AtomType), intent(in)                               :: atom
  integer, intent(in)                                       :: icell, id
  double precision 											:: lambda_limit!, sigma_e
  double precision, dimension(NLTEspec%Nwaves) 				:: scatt
  
  !if (atom%ID /= "H") RETURN
  
  lambda_limit = 121.6d0
  scatt = 0d0

  where(NLTEspec%lambda >= lambda_limit)
   scatt = (1d0 + (1566d0/NLTEspec%lambda)**2.d0 + &
   			(1484.d0/NLTEspec%lambda)**4d0)*(966.d0/NLTEspec%lambda)**4d0
  end where

  if (lstore_opac) then
   NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + &
   									 scatt * sigma_e * Hydrogen%n(icell,1) !m^-1
  else
   NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + &
   									 scatt * sigma_e * Hydrogen%n(icell,1)
  end if
   
 RETURN
 END SUBROUTINE HI_Rayleigh
 
 SUBROUTINE HeI_Rayleigh(id, icell)!, atom)
 ! ------------------------------------------------------------- !
  ! Rayleigh scattering on neutral Helium.
  ! Baschek & Scholz 1982

 ! ------------------------------------------------------------- !
  !type (AtomType), intent(in)                               :: atom
  integer, intent(in)                                       :: icell, id
  double precision 											:: lambda_limit!, sigma_e
  double precision, dimension(NLTEspec%Nwaves) 				:: scatt
  
  !if (atom%ID /= "He") RETURN
  if (.not.associated(Helium)) RETURN
  write(*,*) "Computing HE I rayleigh"
  
  lambda_limit = 121.6d0
  scatt = 0d0


  where(NLTEspec%lambda >= lambda_limit)
   scatt = 4d0 * (1d0 + (669d0/NLTEspec%lambda)**2.d0 + &
   			(641d0/NLTEspec%lambda)**4d0)*(379d0/NLTEspec%lambda)**4d0
  end where

  if (lstore_opac) then
   NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + &
   									 scatt * sigma_e * Helium%n(icell,1) !m^-1
  else
   NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + &
   									 scatt * sigma_e * Helium%n(icell,1)
  end if
   
 RETURN
 END SUBROUTINE HeI_Rayleigh
 
 SUBROUTINE Rayleigh(id, icell, atom)
 ! ------------------------------------------------------------- !
  ! Rayleigh scattering by transitions from the ground state of
  ! neutral atoms.
  ! Sums scattering crosssections of all bound-bound
  ! transitions from the groundstate of the atom with lamda_red
  ! (the longest wavelength covered by the transition) less than wavelength 
  ! lambda.
  !
  ! See:
  ! -- Hubeny & Mihalas 2014 p. 155, eq. 6.44-6.45
  ! SigmaR \propto sigmaT * sum(f/(lambda/lambda0 - 1))**2
  !
  !
  ! Return scattering cross-section for all wavelengths
 ! ------------------------------------------------------------- !
  use math, only            								 : dpow
  type (AtomType), intent(in)                               :: atom
  integer, intent(in)                                       :: icell, id
  logical                                                   :: res
  integer                                                   :: kr, k, la
  double precision                                          :: lambda_red, lambda_limit!, &
   															   !sigma_e
  double precision, dimension(NLTEspec%Nwaves) 				:: fomega, lambda2, scatt
  
  lambda_limit = LONG_RAYLEIGH_WAVE
!   sigma_e = 8.0*PI/3.0 * &
!          dpow(Q_ELECTRON/(dsqrt(4.0*PI*EPSILON_0) *&
!          (dsqrt(M_ELECTRON)*CLIGHT)), 4.d0)
  res = .false. !init, if .false. do not consider Rayleigh
                !of this atom in opac module.
                
  scatt = 0d0
  if (atom%ID=="He") write(*,*) "Rayleigh for Helium" !a test to remove later

  if (atom%stage(1) /= 0) then
   write(*,*) "Lowest level of atom is not a neutral stage"
   write(*,*) "Not computing rayleigh for this atom", atom%ID
   res = .false.
   RETURN
  end if
  
 !   find lambda_red, the longest wavelength covered by the
 !   transition less than lambda for this atom.
   do kr=1,atom%Nline
    !if (atom%lines(kr)%Nred==-99 .and. atom%lines(kr)%Nblue==-99) CYCLE
    ! BEWARE: a line expands from Nblue to Nred but the reddest wavelength
    ! is line%lambda(Nlambda) = NLTEspec%lambda(Nred) by construction.
    !In getlambda.f90, a new grid is constructed and lambda(Nlambda) might not
    !exist anymore..
    if (atom%lines(kr)%i == 1) then
     lambda_red = NLTEspec%lambda(atom%lines(kr)%Nred) !reddest wavelength including
    														!velocity shifts.
     lambda_limit = min(lambda_limit, lambda_red)
    end if
   end do
  fomega = 0.0
  lambda2 = 0d0

 do kr=1,atom%Nline
  !if (atom%lines(kr)%Nred==-99 .and. atom%lines(kr)%Nblue==-99) CYCLE
  lambda_red =  NLTEspec%lambda(atom%lines(kr)%Nred)
   if (atom%lines(kr)%i == 1) then
    where((NLTEspec%lambda > lambda_limit).and.(NLTEspec%lambda > lambda_red))
     lambda2 = 1./((NLTEspec%lambda/atom%lines(kr)%lambda0)**2 -1.)
     fomega = fomega + (lambda2)**2 *atom%lines(kr)%fosc
    end where
   end if 
 end do

  !at worst we add zero and res=.false.
  scatt = sigma_e * fomega * atom%n(1,icell) !m^-1 = m^2 * m^-3
  
  if ((MAXVAL(scatt) > 0))  res = .true.
  
  if (res .and. lstore_opac) then
     NLTEspec%AtomOpac%Kc(icell,:,1) = NLTEspec%AtomOpac%Kc(icell,:,1) + scatt
  else if (res .and. .not.lstore_opac) then
     NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + scatt
  end if

 RETURN
 END SUBROUTINE Rayleigh
!  
!  FUNCTION Rayleigh(icell, atom, scatt) result(res)
!  ! ------------------------------------------------------------- !
!   ! Rayleigh scattering by transitions from the ground state of
!   ! neutral atoms.
!   ! Sums scattering crosssections of all bound-bound
!   ! transitions from the groundstate of the atom with lamda_red
!   ! (the longest wavelength covered by the transition) less than wavelength 
!   ! lambda.
!   !
!   ! See:
!   ! -- Hubeny & Mihalas 2014 p. 155, eq. 6.44-6.45
!   ! SigmaR \propto sigmaT * sum(f/(lambda/lambda0 - 1))**2
!   !
!   !
!   ! Return scattering cross-section for all wavelengths
!  ! ------------------------------------------------------------- !
!   use math, only            								 : dpow
!   type (AtomType), intent(in)                               :: atom
!   integer, intent(in)                                       :: icell
!   logical                                                   :: res
!   double precision, dimension(NLTEspec%Nwaves), intent(out) :: scatt
!   integer                                                   :: kr, k, la
!   double precision                                          :: lambda_red, lambda_limit, &
!    															   sigma_e
!   double precision, dimension(NLTEspec%Nwaves) 				:: fomega, lambda2
!   
!   lambda_limit = LONG_RAYLEIGH_WAVE
!   sigma_e = 8.0*PI/3.0 * &
!          dpow(Q_ELECTRON/(dsqrt(4.0*PI*EPSILON_0) *&
!          (dsqrt(M_ELECTRON)*CLIGHT)), 4.d0)
!   res = .false. !init, if .false. do not consider Rayleigh
!                 !of this atom in opac module.
!                 
!   scatt = 0.
! 
!   if (atom%stage(1) /= 0) then
!    write(*,*) "Lowest level of atom is not a neutral stage"
!    write(*,*) "Not computing rayleigh for this atom", atom%ID
!    res = .false.
!    RETURN
!   end if
!   
!  !   find lambda_red, the longest wavelength covered by the
!  !   transition less than lambda for this atom.
!    do kr=1,atom%Nline
!     ! BEWARE: a line expands from Nblue to Nred but the reddest wavelength
!     ! is line%lambda(Nlambda) = NLTEspec%lambda(Nred) by construction.
!     !In getlambda.f90, a new grid is constructed and lambda(Nlambda) might not
!     !exist anymore..
!     if (atom%lines(kr)%i == 1) then
!      lambda_red = NLTEspec%lambda(atom%lines(kr)%Nred) !reddest wavelength including
!     														!velocity shifts.
!      lambda_limit = min(lambda_limit, lambda_red)
!     end if
!    end do
!   fomega = 0.0
!   lambda2 = 0d0
! 
!  do kr=1,atom%Nline
!   lambda_red =  NLTEspec%lambda(atom%lines(kr)%Nred)
!    if (atom%lines(kr)%i == 1) then
!     where((NLTEspec%lambda > lambda_limit).and.(NLTEspec%lambda > lambda_red))
!      lambda2 = 1./((NLTEspec%lambda/atom%lines(kr)%lambda0)**2 -1.)
!      fomega = fomega + (lambda2)**2 *atom%lines(kr)%fosc
!     end where
!    end if 
!  end do
! 
!   scatt = sigma_e * fomega * atom%n(1,icell) !m^-1 = m^2 * m^-3
!   
!   if ((MAXVAL(scatt).gt.0)) then
!     res = .true.
!     RETURN
!   end if
! 
!  RETURN
!  END FUNCTION Rayleigh
END MODULE rayleigh_scattering
