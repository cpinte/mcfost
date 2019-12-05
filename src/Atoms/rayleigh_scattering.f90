MODULE rayleigh_scattering

 use atom_type, only       : AtomType
 use atmos_type, only      : atmos, Hydrogen, Helium
 use spectrum_type, only   : NLTEspec
 use constant
 use parametres
 
 use mcfost_env, only : dp

 IMPLICIT NONE

 real(kind=dp), parameter :: LONG_RAYLEIGH_WAVE=1.d6 !nm

 CONTAINS
 
 SUBROUTINE HI_Rayleigh(id, icell)
 ! ------------------------------------------------------------- !
  ! Rayleigh scattering on neutral Hydrogen.
  ! Baschek & Scholz 1982

 ! ------------------------------------------------------------- !
  integer, intent(in)                                       :: icell, id
  real(kind=dp) 											:: lambda_limit!, sigma_e
  real(kind=dp), dimension(NLTEspec%Nwaves) 				:: scatt
  integer :: k
    
  lambda_limit = 102.6!121.6d0
  scatt = 0d0

  where(NLTEspec%lambda > lambda_limit)
   scatt = (1d0 + (156.6d0/NLTEspec%lambda)**2.d0 + &
   			(148.d0/NLTEspec%lambda)**4d0)*(96.6d0/NLTEspec%lambda)**4d0
  end where

  if (lstore_opac) then
   NLTEspec%AtomOpac%Kc(:,icell,2) = NLTEspec%AtomOpac%Kc(:,icell,2) + &
   									 scatt * sigma_e * Hydrogen%n(1,icell)
  else
   NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + &
   									 scatt * sigma_e * Hydrogen%n(1,icell)
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
  real(kind=dp) 											:: lambda_limit!, sigma_e
  real(kind=dp), dimension(NLTEspec%Nwaves) 				:: scatt
  integer													:: Neutr_index, l
 
  if (.not.associated(Helium)) RETURN
  
  lambda_limit = helium%scatt_limit 
  scatt = 0d0
  
  l = 1
  do while (Helium%stage(l)==0)
   l = l+1
  enddo
  Neutr_index = l - 1


  where(NLTEspec%lambda > lambda_limit)
   scatt = 4d0 * (1d0 + (66.9d0/NLTEspec%lambda)**2.d0 + &
   			(64.1d0/NLTEspec%lambda)**4d0)*(37.9d0/NLTEspec%lambda)**4d0
  end where

  if (lstore_opac) then
   NLTEspec%AtomOpac%Kc(:,icell,2) = NLTEspec%AtomOpac%Kc(:,icell,2) + &
   									 scatt * sigma_e * Helium%n(1,icell)!sum(Helium%n(1:Neutr_index,icell)) !m^-1

  else
   NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + &
   									 scatt * sigma_e * Helium%n(1,icell)
  end if
   
 RETURN
 END SUBROUTINE HeI_Rayleigh
 
 !Not working properly if we remove some lines
 !because some would have indexes of -99 meaning they are removed
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
  real(kind=dp)                                          :: lambda_red, lambda_limit!, &
   															   !sigma_e
  real(kind=dp), dimension(NLTEspec%Nwaves) 				:: fomega, lambda2, scatt
  
  lambda_limit = LONG_RAYLEIGH_WAVE
!   sigma_e = 8.0*PI/3.0 * &
!          dpow(Q_ELECTRON/(dsqrt(4.0*PI*EPSILON_0) *&
!          (dsqrt(M_ELECTRON)*CLIGHT)), 4.d0)
  res = .false. !init, if .false. do not consider Rayleigh
                !of this atom in opac module.
                
  scatt = 0d0
  !if (atom%ID=="He") write(*,*) "Rayleigh for Helium" !a test to remove later

  if (atom%stage(1) /= 0) then
   write(*,*) "Lowest level of atom is not a neutral stage"
   write(*,*) "Not computing rayleigh for this atom", atom%ID
   res = .false.
   RETURN
  end if
  
  !!!!
  !!!! CANNOT REMOVE TRANSITIONS EVEN IF THEY DO NO CONTRIBUTE TO OPAC. BECAUSE WE NEED TO SUM OVER ALL OF THEM
  !!!! TO COMPUTE THE CROSS-SECTIONS
  !!!!
  
 !   find lambda_red, the longest wavelength covered by the
 !   transition less than lambda for this atom.
   do kr=1,atom%Nline
    !!if (.not.atom%at(kr)%lcontrib_to_opac) CYCLE
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
  !!if (.not.atom%at(kr)%lcontrib_to_opac) CYCLE
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
     NLTEspec%AtomOpac%Kc(:,icell,2) = NLTEspec%AtomOpac%Kc(:,icell,2) + scatt
  else if (res .and. .not.lstore_opac) then
     NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + scatt
  end if

 RETURN
 END SUBROUTINE Rayleigh

END MODULE rayleigh_scattering
