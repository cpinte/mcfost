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
    
  lambda_limit = hydrogen%scatt_limit!102.6
  scatt = 0d0

  where(NLTEspec%lambda > lambda_limit)
   scatt = (1d0 + (156.6d0/NLTEspec%lambda)**2.d0 + &
   			(148.d0/NLTEspec%lambda)**4d0)*(96.6d0/NLTEspec%lambda)**4d0
  end where


   NLTEspec%AtomOpac%sca_c(:,icell) =  NLTEspec%AtomOpac%sca_c(:,icell) + &
   									 scatt * sigma_e * sum(Hydrogen%n(1:Hydrogen%Nlevel-1,icell))

   
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
  
  lambda_limit = helium%scatt_limit!50.0_dp
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

!   if (lstore_opac) then
   NLTEspec%AtomOpac%sca_c(:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) + &
   									 scatt * sigma_e * sum(Helium%n(1:Neutr_index,icell)) !m^-1
! 
!   else
!    NLTEspec%AtomOpac%sca_c(:,id) = NLTEspec%AtomOpac%sca_c(:,id) + &
!    									 scatt * sigma_e * Helium%n(1,icell)
!   end if
!    
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
  integer                                                   :: kr, k, la, l, Neutr_index
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
  l = 1
  do while (atom%stage(l)==0)
   l = l+1
  enddo
  Neutr_index = l - 1
  !write(*,*) atom%ID, " neutral index=", Neutr_index
  
  scatt = sigma_e * fomega * sum(atom%n(1:Neutr_index,icell))!atom%n(1,icell) !m^-1 = m^2 * m^-3
    
  NLTEspec%AtomOpac%sca_c(:,icell) = NLTEspec%AtomOpac%sca_c(:,icell) + scatt

  res = .true.

 RETURN
 END SUBROUTINE Rayleigh

END MODULE rayleigh_scattering
! 
! !
! ! Subroutine to compute the opacity due to Rayleigh scattering at
! ! ND depth points. This routine has good accuracy (5%) redward of the
! ! Lyman lines. Within 1000 km/s of Lyman alpha, cros-section is continuous
! ! and constant.
! !
! ! References:
! !         Gavrila (Phys Rev., 1967, 163, p147)i
! !         Hee-Won Lee and Hee Il Kim, 2004, MNRAS, 347, 802
! !         Nussbaumer, H., Schmid, H. M., Vogel, M., 1989, A&A, 211, L27 
! !
! 	SUBROUTINE RAYLEIGH_SCAT(RKI,HYD,AHYD,EDGE_HYD,NHYD,FREQ,ND)
! 	IMPLICIT NONE
! !
! ! Altered 11-Apr-2005 : Major bug fix.
! ! Created 17-Mar-2004
! !
! 	INTEGER NHYD
! 	INTEGER ND
! !
! 	REAL*8 AHYD(NHYD,NHYD)		!No longer used
! 	REAL*8 HYD(NHYD,ND)
! 	REAL*8 EDGE_HYD(NHYD)		!Ground state value used only.
! 	REAL*8 RKI(ND)
! 	REAL*8 FREQ
! !
! ! We include the oscillator strength in the file. That way we don't
! ! have to worry about whether the l states are split.
! !
! 	REAL*8, SAVE ::  FOSC(20)
! 	DATA FOSC(2:20)/4.162D-01,7.910D-02,2.899D-02,1.394D-02,7.799D-03,
! 	1         4.814D-03,3.183D-03,2.216D-03,1.605D-03,1.201D-03,
! 	1         9.214D-04,7.227D-04,5.774D-04,4.686D-04,3.856D-04,
! 	1         3.211D-04,2.702D-04,2.296D-04,1.967D-04/
! !
! 	REAL*8 EDGE_FREQ
! 	REAL*8 T1,T2
! 	REAL*8 DELF
! 	INTEGER I,IRES
! !
! ! Use the results of Gavrila (Phys Rev., 1967, 163, p147) below
! ! the Lyman limit. Below 0.647 we use the fitting formula due to
! ! Ferland (Cloudy). Blueward of Lyman edge did a simple fit to
! ! results of Gavrila.
! !
! 	IF(FREQ .GT. EDGE_HYD(1))THEN
! 	  T1=FREQ/EDGE_HYD(1)
! 	  RKI(1:ND)=RKI(1:ND)+6.65D-15*HYD(1,1:ND)*(1.0D0+1.66D0/SQRT(T1))
! 	  RETURN
! 	ELSE IF(FREQ .LT. 0.647D0*EDGE_HYD(1))THEN
! 	  T1=(FREQ/EDGE_HYD(1))**2
! 	  RKI(1:ND)=RKI(1:ND)+6.65D-15*HYD(1,1:ND)*
! 	1               T1*T1*(1.26D0+5.068D0*T1+708.3D0*(T1**5))
! 	  RETURN
! 	END IF
! !
! ! Now sum up over all the resonances.
! ! See: Nussbaumer, H., Schmid, H. M., Vogel, M., 1989, A&A, 211, L27 
! !
! 	T1=0.0D0
! 	IRES=1
! 	DELF=100
! 	DO I=2,20
! 	  EDGE_FREQ=EDGE_HYD(1)*(1.0D0-1.0D0/I/I)
! 	  T2=(EDGE_FREQ/FREQ)**2-1.0D0
! 	  IF(DELF .GT. ABS(T2))THEN
! 	    DELF=ABS(T2)
! 	    IRES=I
! 	  END IF
! 	  T2=MAX(ABS(T2),0.0001D0)		!Prevent overflow
! 	  IF(EDGE_FREQ .LT. FREQ)T2=-T2
! 	  T1=T1+FOSC(I)/T2
! 	END DO
! !
! ! The 0.35 gives a match to fitting formula above, and prevents resonances
! ! from going to zero.
! !
! 	T1=T1*T1+0.35D0
! 	IF(IRES .GT. 4)IRES=4
! 	T1=MIN(1.0D4*FOSC(IRES)*FOSC(IRES),T1)
! 	RKI(1:ND)=RKI(1:ND)+6.65D-15*T1*HYD(1,1:ND)
! !
! 	RETURN
! 	END
