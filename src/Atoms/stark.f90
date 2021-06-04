!building
MODULE stark

  use atom_type, only : AtomicContinuum, find_continuum, AtomType, AtomicLine
  use atmos_type, only : Hydrogen, Helium, ne, T
  use constant
  use spectrum_type, only : lambda
  use math, only : bezier3_interp, interp2Darr, interp1D, interp2D, convolve

  use constantes, only : tiny_dp
  use mcfost_env, only : dp

  IMPLICIT NONE




  !Griem 1960, 132, 883
  real, dimension(10), private :: gamma, beta
  real, dimension(10,10), private :: T_betgam ! T(beta, gamma)

  data gamma /0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0/
  data beta /0.0, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 15.0, 20.0/

  data T_betgam(1,:) /0.0004, 0.0, 0.992e-2, 0.936e-2, 0.846e-2, 0.591e-2, 0.292e-2, 0.860e-3, 0.370e-3, 0.996e-4/
  data T_betgam(2,:) /0.963e-2, 0.958e-2,0.943e-2,0.912e-2,0.828e-2,0.579e-2,0.313e-2, 0.975e-3, 0.455e-3, 0.159e-3/
  data T_betgam(3,:) /0.912e-2, 0.911e-2, 0.905e-2, 0.862e-2, 0.794e-2, 0.588e-2, 0.356e-2, 0.05e-3, 0.580e-3, 0.252e-3/
  data T_betgam(4,:) /0.874e-2,0.872e-2,0.864e-2,0.826e-2,0.768e-2, 0.593e-2, 0.385e-2, 0.092e-2, 0.666e-3, 0.354e-3/
  data T_betgam(5,:) /0.837e-2, 0.834e-2, 0.826e-2, 0.794e-2, 0.742e-2, 0.588e-2, 0.402e-2, 0.123e-2, 0.729e-3, 0.434e-3/
  data T_betgam(6,:) /0.768e-2, 0.766e-2, 0.759e-2, 0.733e-2, 0.692e-2, 0.570e-2, 0.415e-2, 0.175e-2, 0.822e-3, 0.550e-3/
  data T_betgam(7,:) /0.706e-2, 0.704e-2, 0.699e-2, 0.678e-2, 0.645e-2, 0.545e-2, 0.417e-2, 0.208e-2, 0.888e-3, 0.633e-3/
  data T_betgam(8,:) /0.651e-2, 0.649e-2, 0.646e-2, 0.628e-2, 0.600e-2, 0.518e-2, 0.410e-2, 0.228e-2, 0.938e-3, 0.698e-3/
  data T_betgam(9,:) /0.556e-2, 0.554e-2, 0.552e-2, 0.540e-2, 0.520e-2, 0.461e-2, 0.382e-2, 0.242e-2, 0.001e-2, 0.786e-3/
  data T_betgam(10,:) /0.440e-2, 0.439e-2, 0.438e-2, 0.430e-2, 0.417e-2, 0.379e-2, 0.324e-2, 0.228e-2, 0.043e-2, 0.864e-3/


CONTAINS
  !building
  subroutine Stark_profile(icell, line)
    !Convolve the Stark profile for Lymann and Balmer lines with the line%phi
    !Stark broadening not included in line%phi for these series
    type (AtomicLine), intent(inout) :: line
    integer, intent(in) :: icell
    ! 		integer :: t
    ! 		character(len=3) :: itxt, jtxt
    ! 		real(kind=dp) :: S(line%Nlambda), test_prof(line%Nlambda)
    !
    ! 		!line%lstark_profile = .false.
    !
    ! return
    ! 		if ((ne(icell) < 1d20) .or. (ne(icell) > 1e23)) then
    ! 			return
    ! 		else
    ! 			if (line%i==1) then
    ! 				CALL Stark_lyman(icell, line, S)
    ! 				!line%lstark_profile = .true.
    ! 			else if (line%i==2) then
    ! 				CALL Stark_balmer(icell, line, S)
    ! 				!line%lstark_profile = .true.
    ! 			else
    ! 				!just leave for other series
    ! 				return
    ! 			endif
    ! 		endif
    !
    !
    !   !do convolution now
    !   !It is S(nu) convolved with phi(v) ??
    ! !   line%phi(:,icell) = convolve(lambda(line%Nblue:line%Nred), S(:), line%phi(:,icell))
    !
    !   test_prof(:) = convolve(lambda(line%Nblue:line%Nred), S(:), line%phi(:,icell))
    !
    !   !find last icell to show to remove !!!!!
    ! !   do t=n_cells, 1, -1
    ! !    if (icompute_atomRT(t)>0) exit
    ! !   enddo
    ! !
    ! !  if (icell == t) then
    ! !  write(itxt,'(I3)') line%i
    ! !  write(jtxt,'(I3)') line%j
    ! !  write(*,*) itxt, jtxt
    ! !  open(80, file=trim(itxt)//trim("->")//trim(jtxt)//"_starkprof.txt", status='unknown')
    ! !  do t=1,line%Nlambda
    ! !   write(80,"(4E)") lambda(line%Nblue+t-1), test_prof(t), S(t), line%phi(t,icell)
    ! !  enddo
    ! !  close(80)
    ! !  endif
    !
    return
  end subroutine Stark_profile

  subroutine Stark_Lyman(icell,line, Snu)
    integer, intent(in) :: icell
    type (AtomicLine), intent(in) :: line
    real(kind=dp) :: delta_nus, S0, dnu, ne0, T0, logne, logT
    real(kind=dp), intent(out) :: Snu(line%Nlambda)
    real :: A, alp, B, c, bb
    integer :: np, n, la

    n = 1 !LYmann

    T0 = T(icell)
    ne0 = ne(icell)*1d-6 !cm^-3
    logne = log10(ne0)
    logT = log10(T0)

    np = int(sqrt(line%atom%g(line%j)/2.))

    alp = 2.5

    if (np == 2) then !A is in log
       A = -0.997
       alp = 2.44
    else if (np == 3) then
       A = -0.571
    else if (np == 4) then
       A = -0.572
    else !np > 4
       A = -0.499
    endif

    delta_nus = 10**(A) * ( ne0 )**(0.688) * real(np)**(2.257)

    if (np == 2) then !A is in log
       B = -0.23
    else if (np == 3) then
       B = -0.37
    else if (np == 4) then
       B = -0.26
    else !np > 4
       B = -0.31
    endif

    S0 = 10**(B) * delta_nus**(-1.012)

    if (np == 2) then !A is in log
       bb = 1.23 + 1.62e-2 * logne + 3.51e-2 * logT
    else if (np == 3) then
       bb = -2.44 + 3.98e-1 * logne - 8.79e-3 * logT
    else if (np == 4) then
       bb = 1.09 + 2.66e-2 * logne + 1.72e-1 * logT
    else !np > 4
       bb = -0.91 + 1.37*log10(real(np)) + 1.64e-1 * logne + 4.93e-2 * logT
    endif

    c = alp/bb

    do la=1, line%Nlambda
       dnu = abs( M_TO_NM * (lambda(line%Nblue+la-1) - line%lambda0) * CLIGHT / line%lambda0**2 ) / delta_nus
       Snu(la) = S0 / (1.0 + dnu**b )**c
    enddo

    return
  end subroutine stark_lyman

  subroutine Stark_balmer(icell,line, Snu)
    integer, intent(in) :: icell
    type (AtomicLine), intent(in) :: line
    real(kind=dp) :: delta_nus, S0, dnu, ne0, T0, logne, logT
    real(kind=dp), intent(out) :: Snu(line%Nlambda)
    real :: A, alp, B, c, bb
    integer :: np, n, la

    n = 2 !Balmer

    T0 = T(icell)
    ne0 = ne(icell)*1d-6 !cm^-3
    logne = log10(ne0)
    logT = log10(T0)

    np = int(sqrt(line%atom%g(line%j)/2.))

    alp = 2.5

    if (np == 3) then !A is in log10
       A = -1.16
       alp = 2.45
    else if (np == 4) then
       A = -0.747
       alp = 2.47
    else if (np == 5) then
       A = -0.714
    else !np > 5
       A = -0.656
    endif

    delta_nus = 10**(A) * ( ne0 )**(0.69) * real(np)**(2.377)

    if (np == 3) then !A is in log
       B = -0.21
    else if (np == 4) then
       B = -0.39
    else if (np == 5) then
       B = -0.26
    else !np > 5
       B = -0.31
    endif

    S0 = 10**(B) * delta_nus**(-1.012)

    if (np == 3) then !A is in log
       bb = 1.53 - 1.52e-2 * logne + 6.71e-2 * logT
    else if (np == 4) then
       bb = -0.12 + 2.30e-1 * logne + 1.80e-1 * logT
    else if (np == 5) then
       bb = 1.61 + 0.90e-2 * logne + 1.98e-1 * logT
    else !np > 5
       bb = -1.56 + 1.63 * log(real(np)) + 1.90e-1 * logne + 4.25e-2 * logT
    endif

    c = alp/bb

    do la=1, line%Nlambda
       dnu = abs( M_TO_NM * (lambda(line%Nblue+la-1) - line%lambda0) * CLIGHT / line%lambda0**2) / delta_nus
       Snu(la) = S0 / (1.0 + dnu**bb )**c
    enddo

    return
  end subroutine stark_balmer

  !needless to say that j > i

  FUNCTION Kij (i,j,Z)
    !eq. 33 Griem 1960, 132, 883
    integer, intent(in) :: i, j, Z
    real :: Kij, Z5

    Z5 = real(Z*Z*Z*Z*Z)

    !(i*j)**4 / (j**2 - i**2)
    Kij = 5.5e-5 * real(i*i*i*i * j*j*j*j) / real(j*j - i*i) / Z5

    RETURN
  END FUNCTION Kij

  FUNCTION electron_damping (i,j,Z,T,ne)
    !eq. 29 Griem 1960, 132, 883
    !with few corrections
    !ne is in cm^-3 in the formula
    real(kind=dp) :: electron_damping
    real(kind=dp), intent(in) :: ne, T
    real(kind=dp) :: A
    integer :: i,j,Z

    !   A = real(j*j*j*j*j + i*i*i*i*i) / real(j*j - i*i)
    !
    !   electron_damping = ( 5.6e-6 * ne**(1./3.) / real(Z) / sqrt(T) ) * &
    !     ( dlog10( 4d6 * T * Z / real(j*j) / sqrt(ne) ) - 0.125 ) * A

    A = real(j*j*j*j - 3*i*j*j) / real(j*j - i*i)

    electron_damping = ( 3.78e-5 * (1e-6*ne)**(1./3.) / real(Z) / sqrt(T) ) * &
         dlog10( 4d6 * T * Z / real(j*j) / sqrt(1e-6*ne) )  * A


    RETURN
  END FUNCTION electron_damping

  FUNCTION F0 (ne)
    ! Holstmark normal field strength
    real(kind=dp) F0
    real(kind=dp), intent(in) :: ne
    !! Remember
    !real :: gamma = 2.61
    !real(kind=dp) :: Cp = Q_ELECTRON
    !1 esu = 0.1 / clight Coulomb
    ! q cgs = Q_ELECTRON * 10 * CLIGHT

    !with electron charge in cgs units (esu = cm^3/2 g^2/2 s^-1)
    !and electron densities in cm-3, F0 = 2.61 * qel * ne**(2./3.)
    !with 2.61 * q = 1.25e-9. = 2.61 * (Q_ELECTRON * 10.0 * CLIGHT)
    F0 = 1e-4 * 1.25e-9 * ne**(2./3.) !cgs, field strength unit
    !the 1e-4 is here because, ne**2./3 (ne in m^-3) has to be converted in (cm^-3)**2./3.

    RETURN
  END FUNCTION F0

END MODULE stark
