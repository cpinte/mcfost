module occupation_probability

  use constantes
  use atom_type, only : atomic_orbital_sqradius

   implicit none

   contains


  function keq9dot70(i)
    real(kind=dp), intent(in) :: i
    real(kind=dp) :: keq9dot70

    if (i <= 3.) then
       keq9dot70 = 1.0_dp
    else !16/3 * i * (i+1)**-2
       keq9dot70 = 5.3333333333 * i / (i+1.0)**2.
    endif

    RETURN
   end function keq9dot70

   function wocc_n(t, ne, n, Zr, Zp, nH1, nHe1)
    !neglecting neutral at the moment: need nZ, and nL for the neutral radiator considered
    !with principal quantum number n (n_eff if hydrogenic like)

    !Assuming ionized perturbers are predominently singly ionisze Zp = 1 * ne
    !For neutral contributions, assumes that nH and nHe are in the ground state
    !if T>> only ions contribute, if T<< only neutrals and H is mostly in its ground state

    real(kind=dp), intent(in) :: t, ne
    real, intent(in)    :: Zr ! radiator charge
    real, intent(in)    :: Zp ! perturber (ions) charge
    real(kind=dp), intent(in) :: n !principal quantum number
    real(kind=dp), intent(in) :: nH1 !ground state neutral H
    real(kind=dp), optional :: nhe1
    real(kind=dp) :: wocc_n
    real(kind=dp) :: w_neutr, r1, rp1, rp2, npop1, npop2, a0fourpi_three
    real(kind=dp) :: w_ion, f, x, a , betac
    real, parameter ::  c1 = 0.1402, c2 = 0.1285
    integer :: nl, nz

    a0fourpi_three = (4./3.) * pi * RBOHR*RBOHR*RBOHR

    !n=1, l=0, Z=1 ground state of H I
    rp1 = sqrt(atomic_orbital_sqradius(1.0_dp, 0, 1)) ! about 1.74 a0
    !n=1, l=0, Z=1 ? ground state of He I
    rp2 = sqrt(atomic_orbital_sqradius(1.0_dp, 0, 2)) !about 0.87 a0

    nl = 0. !need to be extracted from label
    nZ = int(Zr) + 1!?
    !or
    !nZ = 2 for He, 1 for H etc
    r1  = sqrt(atomic_orbital_sqradius(n, nl, nZ))

    npop1 = nH1
    npop2 = 0.0_dp
    if (present(nhe1)) then
       !!write(*,*) "Check occupation probab with helium"
       !need to check though
       npop2 = nHe1
    endif

    !init at 1 because wocc is the product of the two probabilities
    w_neutr = exp( -a0fourpi_three * (npop1*(r1+rp1)**3 + npop2*(r1+rp2)**3))
    !write(*,*) a0fourpi_three, r1, rp1, rp2
    !write(*,*) "neutr=",w_neutr, npop1*1d-6,ne(icell)*1d-6,npop2*1d-6, n, T(icell)
    !!w_neutr = 1.0

    !eq. 4.71 Hummer & Mihalas 1988
    !Chapter 9, eq. 9.71, Hubeny & Mihalas, Stellar Atmospheres
    !see also Hubeny & Hummer & Lanz 1994

    !1d4 = ne*-2/3. with ne in m^-3 to cm^-3
    betac  = 1d4 * 8.3d14 * ne**(-2./3.) * Zp*Zp*Zp * keq9dot70(n) / n / n / n / n !unitless it is a ratio

    !1d-1 convert ne m^-3 to ne cm^-3
    a = 1d-1 * 0.09 * ne**(1./6.) / sqrt(T)
    x = (1.0_dp + a)**(3.15)
    f = c1 * (x + 4.0 * Zr * a*a*a)*betac*betac*betac / (1.0 + c2*x*sqrt(betac*betac*betac))

    w_ion = f / (1.0_dp + f)

    wocc_n = w_neutr * w_ion

    return
   end function wocc_n

   function D_i(t, ne, i, Zr, Zp, lambda, lambda0, chi0, nH1)
    !Dissolve fraction
    !Däppen, Anderson & Mihalas Apj 319, 195, 1987
    !chi0 is the ionisation potential of the ion on J
    real(kind=dp), intent(in) :: t, ne, nH1
    real(kind=dp) :: D_i
    real(kind=dp), intent(in) :: i, lambda,  lambda0, chi0!!, dEi !E(i+1) - E(i)
    real, intent(in) :: Zr, Zp
    real(kind=dp) :: m, w, hnu!!, taper, redcut
    real(kind=dp) :: w1, wi
    real :: Zsq

    Zsq = (Zr+1)*(Zr+1) ! 1 for H I

    hnu = HP * C_LIGHT  /  (NM_TO_M * lambda)
    m = 1.0/i/i - hnu / chi0 / Zsq !I'm not sure here, in Dam it is Z**2 * IonH
    !doesn't change for Hydrogen

    w1 = wocc_n(t, ne, real(i,kind=dp), Zr, Zp, nH1)

    if (lambda <= lambda0) then
       D_i = 1.0_dp
    else
       if (m > 0) then
          m = 1.0/sqrt(m)
          wi = wocc_n(t, ne, m, Zr, Zp, nH1)
          D_i = 1. - wi/w1
       else
          D_i = 1.0_dp
       endif
    endif

    return
   end function D_i

  !Building
!   FUNCTION D_i_b(icell, m, Zr, Zp, lambda, lambda0)
!     !Dissolve fraction
!     !Similar to CMFGEN, SYNSPEC ?
!     integer, intent(in) :: icell
!     real(kind=dp) :: D_i_b
!     real(kind=dp), intent(in) :: m, lambda,  lambda0
!     real, intent(in) :: Zr, Zp
!     real(kind=dp) :: w!!, taper, redcut

!     ! ? in Hubeny Mihalas
!     !m = 1.0/i/i - nu/nu0


!     if (lambda <= lambda0) then
!        D_i_b = 1.0_dp
!     else
!        if (m > 2*(Zr+1)) then
!           D_i_b = 1. - wocc_n(icell, m, Zr, Zp,hydrogen%n(1,icell))
!        else
!           D_i_b = 0.0_dp
!        endif
!     endif

!     RETURN
!   END FUNCTION D_i_b

  !still building but they suggest to use the DAM version anyway
  !   FUNCTION D_i_hhl(icell, kc, Zr, Zp, atom)
  !    !Hubeny, Hummer & Lanz A&A, 282, 151, 1994
  !    !dissolve fraction
  !    integer, intent(in) :: icell, kc
  !    type(AtomType), intent(in) :: atom
  !    real(kind=dp) :: D_i_hhl(atom%continua(kc)%Nlambda) !?
  !    real, intent(in) :: Zr, Zp
  !    integer :: la, Nl, Nb, Nr, kr, idl, lal
  !    real(kind=dp) :: w, j, lambda, l0, l1
  !
  !    D_i_hhl(:) = 0.0_dp
  !    Nb = atom%continua(kc)%Nblue
  !    Nr = atom%continua(kc)%Nred
  !    Nl = atom%continua(kc)%Nlambda
  !
  !    do la=1, Nl
  !     idl = la+Nb-1
  !     lambda = nltespec%lambda(idl)
  !
  !     if (lambda <= atom%continua(kc)%lambda0) then
  !
  !      D_i_hhl(la) = 1.0_dp !no extrapolation
  !
  !     else
  !      do kr=1, atom%Nline
  !       l0 = nltespec%lambda(atom%lines(kr)%Nblue)
  !       l1 = nltespec%lambda(atom%lines(kr)%Nred)
  !       if (l0 >= lambda .and. lambda <= l1) then
  !
  !         lal = idl + 1 - nltespec%lambda(atom%lines(kr)%Nblue)
  !         w = wocc_n(icell, real(atom%lines(kr)%j,kind=dp), Zr, Zp)
  !         D_i_hhl(la) = D_i_hhl(la) + (1.-w) * atom%lines(kr)%fosc * atom%lines(kr)%phi(lal,icell)
  !
  !       endif
  !
  !      enddo
  !      D_i_hhl(la) = D_i_hhl(la) / atom%continua(kc)%alpha(la)
  !
  !     endif
  !
  !    enddo
  !
  !   RETURN
  !   END FUNCTION D_i_hhl


END MODULE occupation_probability
