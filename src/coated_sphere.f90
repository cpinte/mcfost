module coated_sphere

  use parametres
  use constantes
  use grains
  use utils, only : gauleg

  implicit none

  complex(kind=dp), parameter :: czero = (0.0_dp,0.0_dp)
  complex(kind=dp), parameter :: ci = (0.0_dp,1.0_dp)

  integer, save ::   maxmsg = 0, nummsg = 100
  logical, save :: msglim = .false.

contains

  subroutine mueller_coated_sphere(lambda,taille_grain,wl,amu1,amu2,amu1_coat,amu2_coat,qext,qsca,gsca)
    !***************************************************************
    ! calcule les elements de la matrice de diffusion a partir de
    ! la sous-routine dmilay (coated grains)
    !
    !        calcule aussi "g" = le parametre d'asymetrie
    !
    ! C. Pinte 19 Mars 2008
    ! G. Duchene 22 Avril 2011
    !****************************************************************

    implicit none

    integer, intent(in) :: lambda, taille_grain
    real, intent(in) :: amu1, amu2, amu1_coat, amu2_coat
    real, intent(in) :: wl
    real, intent(out) :: qext, qsca, gsca

    integer :: j, nang

    complex, dimension(nang_scatt+1) :: S1,S2

    real :: rcore, rshell, wvno, gqsc
    real :: vi1, vi2, norme, somme_prob, qbs, x, theta, dtheta
    complex :: refrel, refrel_coat
    real, dimension(0:nang_scatt) ::  S11,S12,S33,S34


    refrel = cmplx(amu1,amu2)
    refrel_coat = cmplx(amu1_coat,amu2_coat)
    refrel = conjg(refrel)  ! to match convetion in dmilay (negative img part)
    refrel_coat = conjg(refrel_coat)  ! to match convention in dmilay (negative img part)

    if (modulo(nang_scatt,2)==1) then
       write(*,*) "ERROR : nang_scatt must be an EVEN number"
       write(*,*) "Exiting"
       stop
    endif

    ! Si fonction de HG, on ne calcule pas la fonction de phase
    if (aniso_method==2) then
       nang=1
    else
       nang= (nang_scatt+1) / 2 + 1
    endif

    rcore=r_core(taille_grain)
    rshell=r_grain(taille_grain)

    wvno= 2.0 * pi / wl
    x = rshell * wvno

    call dMiLay(rcore,rshell,wvno,refrel_coat,refrel,nang, qext,qsca,qbs,gqsc,s1,s2)
    gsca = gqsc / qsca ! dmilay return gsca * qsca

    if (lforce_HG) gsca = forced_g


    ! Passage des valeurs dans les tableaux de mcfost
    if (aniso_method==1) then

       !  QABS=QEXT-QSCA
       ! Calcul des elements de la matrice de diffusion
       ! indices decales de 1 par rapport a bhmie
       do J=0,nang_scatt
          vi1 = cabs(S2(J+1))*cabs(S2(J+1))
          vi2 = cabs(S1(J+1))*cabs(S1(J+1))
          s11(j) = 0.5*(vi1 + vi2)
          s12(j) = 0.5*(vi1 - vi2)
          s33(j)=real(S2(J+1)*conjg(S1(J+1)))
          s34(j)=aimag(S2(J+1)*conjg(S1(J+1)))
       enddo !j

       ! Integration S11 pour tirer angle
       if (scattering_method==1) then
          prob_s11(lambda,taille_grain,0)=0.0
          dtheta = pi/real(nang_scatt)
          do j=2,nang_scatt ! probabilite de diffusion jusqu'a l'angle j, on saute j=0 car sin(theta) = 0
             theta = real(j)*dtheta
             prob_s11(lambda,taille_grain,j)=prob_s11(lambda,taille_grain,j-1)+s11(j)*sin(theta)*dtheta
          enddo

          ! s11 est calculee telle que la normalisation soit: 0.5*x**2*qsca
          ! il y a un soucis numerique quand x >> 1 car la resolution en angle n'est pas suffisante
          ! On rate le pic de diffraction (en particulier entre 0 et 1)
          somme_prob = 0.5*x**2*qsca
          prob_s11(lambda,taille_grain,1:nang_scatt) = prob_s11(lambda,taille_grain,1:nang_scatt) + &
               somme_prob - prob_s11(lambda,taille_grain,nang_scatt)

          ! Normalisation de la proba cumulee a 1
          prob_s11(lambda,taille_grain,:)=prob_s11(lambda,taille_grain,:)/somme_prob
       endif ! scattering_method==1

       do j=0,nang_scatt
          if (scattering_method==1) then ! Matrice de Mueller par grain
             ! Normalisation pour diffusion selon fonction de phase (tab_s11=1.0 sert dans stokes)
             norme=s11(j)
          else ! Sinon normalisation a Qsca
             ! La normalisation par default : 0.5*x**2*Qsca --> correction par 0.5*x**2
             norme = 0.5 * x**2
          endif

          tab_s11(j,taille_grain,lambda) = s11(j) / norme
          tab_s12(j,taille_grain,lambda) = s12(j) / norme
          tab_s33(j,taille_grain,lambda) = s33(j) / norme
          tab_s34(j,taille_grain,lambda) = s34(j) / norme
       enddo

    endif ! aniso_method ==1

    return

  end subroutine mueller_coated_sphere

  ! **********************************************************************

  subroutine mueller_DHS(lambda,taille_grain,wl,amu1,amu2,qext,qsca,gsca)
    !***************************************************************
    ! Adapte de la routine q_dhs de Michiel Min
    !
    ! C. Pinte 30 janvier 2013
    !****************************************************************

    implicit none

    integer, intent(in) :: lambda, taille_grain
    real, intent(in) :: amu1, amu2
    real, intent(in) :: wl
    real, intent(out) :: qext, qsca, gsca

    integer :: i,j, nang, ipop

    real :: qext_HS, qsca_HS, qbs_HS, gqsc_HS, gsca_HS ! pour 1 HS
    complex, dimension(nang_scatt+1) :: S1,S2, s1_HS, s2_HS

    real :: a, rcore, rshell, wvno, factor, cext, csca
    real :: vi1, vi2, norme, somme_prob, theta, dtheta, x
    complex :: refrel, refrel_coat
    real, dimension(0:nang_scatt) :: S11,S12,S33,S34

    integer :: N_vf
    real(kind=dp), dimension(:), allocatable :: f, wf

    N_vf = min(max(20., deux_pi * r_grain(taille_grain)/tab_lambda(lambda)), 100.)
    allocate(f(N_vf), wf(N_vf))

    refrel = cmplx(1.0,0.0) ! vide
    refrel_coat = cmplx(amu1,amu2)
    refrel_coat = conjg(refrel_coat)  ! to match convention in dmilay (negative img part)

    if (modulo(nang_scatt,2)==1) then
       write(*,*) "ERROR : nang_scatt must be an EVEN number"
       write(*,*) "Exiting"
       stop
    endif

    ! Si fonction de HG, on ne calcule pas la fonction de phase
    if (aniso_method==2) then
       nang = 1
    else
       nang = (nang_scatt+1) / 2 + 1
    endif

    a = r_grain(taille_grain)

    wvno= 2.0 * pi / wl
    x = wvno * a

    ipop = grain(taille_grain)%pop

    ! Calcul des poids pour integration de Gauss-Legendre
    call GauLeg(0.0_dp,real(dust_pop(ipop)%dhs_maxf,kind=dp),f,wf,N_vf) ; wf = wf/sum(wf)
    ! todo : a ne faire que pour 1 taille de grain et 1 lambda, sauf si n_vf change

    cext=0 ; csca=0 ; gsca=0 ; s1=0 ; s2=0
    do i=1,N_vf
       rshell = a/((1.-f(i))**un_tiers)
       rcore = rshell * f(i)**un_tiers

       call dMilay(rcore,rshell,wvno,refrel_coat,refrel,nang, qext_HS,qsca_HS,qbs_HS,gqsc_HS,s1_HS,s2_HS)
       gsca_HS = gqsc_HS / qsca_HS ! dmilay return gsca * qsca

       if(qext_HS < 0_dp) qext_HS = 0_dp
       if(qsca_HS < 0_dp) qsca_HS = 0_dp

       factor = pi*rshell**2 * wf(i)
       cext = cext + factor * qext_HS  ! pas sur que ce soit bon, verifier normalization qext, (rshell/lambda)**2 doit deja y etre
       csca = csca + factor * qsca_HS
       gsca = gsca + factor * qsca_HS * gsca_HS

       s1 = s1 + s1_HS * wf(i)
       s2 = s2 + s2_HS * wf(i)
    enddo
    if (csca > 0.) gsca = gsca/csca

    if (lforce_HG) gsca = forced_g

    factor = pi*a**2 ! va etre remultiplie par pi*a**2 a la sortie dans dust.f90
    qext = cext/factor
    qsca = csca/factor

    ! Passage des valeurs dans les tableaux de mcfost
    if (aniso_method==1) then

       ! QABS=QEXT-QSCA
       ! Calcul des elements de la matrice de diffusion
       ! indices decales de 1 par rapport a bhmie
       do J=0,nang_scatt
          vi1 = cabs(S2(J+1))*cabs(S2(J+1))
          vi2 = cabs(S1(J+1))*cabs(S1(J+1))
          s11(j) = 0.5*(vi1 + vi2)
          s12(j) = 0.5*(vi1 - vi2)
          s33(j)=real(S2(J+1)*conjg(S1(J+1)))
          s34(j)=aimag(S2(J+1)*conjg(S1(J+1)))
       enddo !j

       prob_s11(lambda,taille_grain,0)=0.0
       dtheta = pi/real(nang_scatt)
       do j=2,nang_scatt ! probabilite de diffusion jusqu'a l'angle j, on saute j=0 car sin(theta) = 0
          theta = real(j)*dtheta
          prob_s11(lambda,taille_grain,j)=prob_s11(lambda,taille_grain,j-1)+s11(j)*sin(theta)*dtheta
       enddo

       ! s11 est calculee telle que la normalisation soit: 0.5*x**2*qsca
       ! il y a un soucis numerique quand x >> 1 car la resolution en angle n'est pas suffisante
       ! On rate le pic de diffraction (en particulier entre 0 et 1)
       somme_prob = 0.5*x**2*qsca
       prob_s11(lambda,taille_grain,1:nang_scatt) = prob_s11(lambda,taille_grain,1:nang_scatt) + &
            somme_prob - prob_s11(lambda,taille_grain,nang_scatt)

       ! Normalisation de la proba cumulee a 1
       prob_s11(lambda,taille_grain,:)=prob_s11(lambda,taille_grain,:)/somme_prob

       do j=0,nang_scatt
          if (scattering_method==1) then ! Matrice de Mueller par grain
             ! Normalisation pour diffusion selon fonction de phase (tab_s11=1.0 sert dans stokes)
             norme=s11(j)
          else ! Sinon normalisation a Qsca
             ! La normalisation par default : 0.5*x**2*Qsca --> correction par 0.5*x**2
             norme = 0.5 * x**2
          endif

          tab_s11(j,taille_grain,lambda) = s11(j) / norme
          tab_s12(j,taille_grain,lambda) = s12(j) / norme
          tab_s33(j,taille_grain,lambda) = s33(j) / norme
          tab_s34(j,taille_grain,lambda) = s34(j) / norme
       enddo

    endif ! aniso_method ==1

    deallocate(f,wf)
    return

  end subroutine mueller_DHS

  ! **********************************************************************

  subroutine dmilay(rcore, rshell, wvno, rindsh, rindco, &
       numang, qext, qsca, qbs, gqsc, s1, s2)
    ! *******************************************************************
    !              DOUBLE PRECISION version of MieLay
    ! *******************************************************************
    ! It is recommended to use this DOUBLE PRECISION version for IEEE
    ! arithmetic (32-bit floating point) computers, just to be safe.
    ! If computer time is critical, back-to-back tests with the single
    ! precision version should be done for ranges of radii and refractive
    ! index relevant to your particular problem, before adopting the
    ! single precision version.  This version is also recommended for
    ! cases of large size parameter (bigger than 10 or so) and/or large
    ! imaginary refractive index (bigger than 1 or so) and also whenever
    ! overflows or strange behavior are encountered in running the
    ! single precision version.  Sometimes the bigger exponent range in
    ! DOUBLE PRECISION is as important as the added precision.
    !
    ! This version is designed to be interchangeable with the single
    ! precision version:  all the floating-point input arguments are
    ! still single precision.  Only the name of the routine has been
    ! changed to prevent confusion (and it is strongly urged not to
    ! change it to MieLay for the same reason).
    !
    !    This subroutine computes electromagnetic scattering by a
    !    stratified sphere (a particle with a spherical core surrounded
    !    by a spherical shell).  The surrounding medium is assumed to
    !    have refractive index unity.  The formulas, manipulated to avoid
    !    the ill-conditioning that plagued earlier formulations, were
    !    published in:
    !
    !        Toon, O. and T. Ackerman, Applied Optics 20, 3657 (1981)
    !
    !    The latest version of this program is available by anonymous ftp
    !    from climate.gsfc.nasa.gov in directory pub/wiscombe.
    !
    !    The program was based on the famous homogeneous sphere
    !    program of Dave, published in:
    !
    !       Dave, J.V., "Subroutines for Computing the Parameters of the
    !          Electromagnetic Radiation Scattered by a Sphere",
    !          IBM Scientific Center, Palo Alto, California,
    !          Report No. 320 - 3236, May 1968
    !
    !       Dave, J.V., Applied Optics 8, 155 (1969)
    !
    !    This was done because the formulas are identical in structure to
    !    those for the homogeneous sphere, except that the coefficients of
    !    the Mie series (commonly denoted by little-a-sub-n and
    !    little-b-sub-n) are much more complicated.
    !
    !    The address of the first author is:
    !
    !         Dr. O. B. Toon (toon@sky.arc.nasa.gov)
    !         NASA Ames Research Center
    !         M.S. 245-3
    !         Moffett Field, CA (USA)
    !
    !    The program was explicitly designed to avoid the ill-conditioning
    !    which was latent in the standard analytic formulation of the
    !    Mie core-shell problem.  Oddly, this ill-conditioning had been
    !    exposed and eliminated in the late 1960s for the homogeneous sphere
    !    problem, but went unrecognized for the core-shell problem, leading
    !    to many incorrect results being published prior to 1981.  In
    !    particular, previous calculations were generally wrong for the
    !    case of a thin carbon shell.
    !
    !    After a number of years of experience, this program seems to have
    !    only two limitations:  a slow degradation as size increases,
    !    and a rapid degradation as size decreases due to lack of explicit
    !    handling of the ill-conditioning in that limit. It has been used
    !    successfully for cases with large imaginary refractive index (both
    !    in core and shell) and even with real refractive index less than
    !    unity.
    !
    !    For too-large particles, internal array sizes will be inadequate,
    !    but this generates an error message and an abort.
    !
    !    It is highly recommended to use the DOUBLE PRECISION version of
    !    this program, called 'DMiLay', on machines with 32-bit floating-
    !    point arithmetic (e.g. all IEEE arithmetic machines).  'DMiLay'
    !    is also available on the network.  'MieLay' may be adequate but
    !    this should be tested by running it back to back with 'DMiLay'.
    !
    !        Note that in 32-bit arithmetic it was impossible to run
    !        'MieLay' above shell radius = 3 (with WVNO=1, so shell size
    !        parameter = 3 also) due to overflow, whereas 'DMiLay' can be
    !        run well above this limit due to a larger range of exponents.
    !
    !    The original version of this program defaulted to a homogeneous
    !    sphere case when the core was less than 10**(-6) the size of the
    !    shell.  It could also have done so when core and shell radii were
    !    equal although it did not.  But this option was dropped since
    !    it could better be tested for in the calling program;  if a
    !    homogeneous sphere case is detected, it is far better to call one
    !    of the programs designed explicitly for that case.
    !
    !    NOTE:  This program requires input scattering angles between
    !           zero and 90 degrees.  Then it does calculations for
    !           those angles, plus all their supplements.  Thus, to get,
    !           e.g., 170 degrees, you must use 10 degrees.
    !
    !    The program was modified and further documented by W. Wiscombe
    !    (wiscombe@climate.gsfc.nasa.gov; NASA Goddard, Code 913,
    !    Greenbelt, MD 20771), including:
    !    ** complex refractive indices of shell and core submitted as 2
    !          complex arguments rather 4 real arguments
    !    ** scattering angles submitted as cosines rather than angles
    !          themselves (anticipating important usage where cosines remain
    !          constant while program is called repeatedly for different
    !          shell and/or core sizes in order to integrate over size)
    !    ** returning scattering matrix elements M1, M2, D21, S21 separately
    !          rather than bundled into awkward data structure ELTRMX
    !    ** defining new arrays S1 and S2, the complex scattering
    !          amplitudes, for which ELTRMX had formerly been used as
    !          temporary storage (allowing easy modification to return S1
    !          and S2 through argument list rather than M1, M2, D21, S21)
    !    ** use of internal work arrays, with guards against blowing their
    !          dimensions, rather than submitting these arrays as arguments;
    !          allows easy conversion to ALLOCATABLE arrays in Fortran-90
    !    ** elimination of all GO TOs but one
    !    ** elimination of dangerous EQUIVALENCE statement (it was used to
    !          equate variables, not to share storage)
    !    ** more mnemonic variable names
    !    ** error-checking input arguments
    !    ** making data structures more optimal for vectorization
    !          (particularly PI, TAU, and the replacements for ELTRMX);
    !          mainly so that innermost loops are over the first
    !          dimension of all arrays
    !    ** polishing and declaration standardization using NAG Fortran
    !          Tool nag_decs
    !    ** creation of a DOUBLE PRECISION version using NAG Fortran
    !          Tool nag_apt
    !    ** certification by 'flint' (Fortran 'lint', for C folks)
    !
    !    Suggestions for the future:
    !    ** much of this program reflects the belief, true in the late
    !         1960s but not any longer, that complex arithmetic should
    !         be broken into real and imaginary parts for efficiency,
    !         in spite of massive loss of clarity;  the program would
    !         simplify considerably if done in complex variables
    !         (now only partially true)
    !    ** improve treatment of angles greater than 90 degrees
    !         (an awkward and confusing hangover from the old Dave code)
    !    ** phrase program in terms of size parameters and eliminate
    !         input argument WVNO
    !    ** develop special-case formulas for Rcore-->0 for any Rshell,
    !         and for Rshell-->0 (comparing single and double precision
    !         versions showed larger and larger differences as size
    !         parameter fell below about 1.E-2, esp. for gQsc); the
    !         layered sphere formulae are ill-conditioned in this limit
    !         just as the homogeneous sphere formulae are
    !
    !
    !    I N P U T   A R G U M E N T S
    !
    !    (Definition:  size parameter = sphere circumference / wavelength )
    !
    !      Rshell      radius of shell
    !
    !      Rcore       radius of core
    !
    !      WVNO        2*pi / wavelength
    !
    !      RindSh      COMPLEX refractive index of shell (negative
    !                     imaginary part)
    !
    !      RindCo      COMPLEX refractive index of core (negative
    !                     imaginary part)
    !
    !      MU          array of cosines of scattering angles (angles between
    !                     directions of incident and scattered radiation).
    !                     For angles between 90 and 180 degrees, use the
    !                     supplement (180-angle) of the angle instead, so
    !                     0.le.MU.le.1 (see comments below on M1,M2,21,D21)
    !
    !      NumAng      Number of scattering angles for which computations
    !                     are required; should not exceed MaxAng
    !                     (NOTE:  NumAng=0 will suppress the calculation
    !                      of the scattering matrix quantitities  M1, M2,
    !                      S21, D21 and save a lot of computer time)
    !
    !      MaxAng      First dimension of M1,M2,21,D21 in calling program
    !
    !
    !
    !    O U T P U T   A R G U M E N T S
    !
    !      (Definitions for these arguments can be found in H.C. van de
    !       Hulst, Light Scattering By Small Particles, Dover Press, New
    !       York, 1981 (reprint of 1957 edition); abbreviated VDH below)
    !
    !      QEXT     Efficiency factor for extinction (VDH Sec 9.32)
    !               (same as corresponding quantity in MIEV)
    !
    !      Qsca     Efficiency factor for scattering (VDH Sec 9.32)
    !               (same as corresponding quantity in MIEV)
    !
    !      GQSC     average(cosine theta) * Qsca (VDH Sec 9.32)
    !                  (<cos theta> is usually denoted by g, hence
    !                   the name of the variable)
    !               (same as corresponding quantity in MIEV)
    !
    !      QBS      Backscatter cross section.
    !               ( Re(SBACK)**2 + Im(SBACK)**2 ) / (Pi*XSHELL**2)
    !               where the corresponding quantity from MIEV is
    !
    !               SBACK = 0.5*sum(n=1 to inf)((-1)**(n+1)(2n+1)(an-bn))
    !
    !               and an,bn are ACOE,BCOE below.
    !
    !      M1(j,k)  Element M1 of scattering matrix F' (VDH Sec 5.14);
    !                  M1(j,1) refers to angle with cosine MU(j);
    !                  M1(j,2) refers to supplement of that angle.
    !               (Be sure to type REAL in calling program.)
    !
    !      M2(j,k)  Element M2 of scattering matrix F' (VDH Sec 5.14);
    !                  M2(j,1) refers to angle with cosine MU(j);
    !                  M2(j,2) refers to supplement of that angle.
    !               (Be sure to type REAL in calling program.)
    !
    !     S21(j,k)  Element S21 of scattering matrix F' (VDH Sec 5.14);
    !                  S21(j,1) refers to angle with cosine MU(j);
    !                  S21(j,2) refers to supplement of that angle.
    !
    !     D21(j,k)  Element D21 of scattering matrix F' (VDH Sec 5.14);
    !                  D21(j,1) refers to angle with cosine MU(j);
    !                  D21(j,2) refers to supplement of that angle.
    !
    !
    !    L O C A L   V A R I A B L E S
    !
    !      ACOE     (COMPLEX) Mie coeff. little-a-sub-n
    !      BCOE     (COMPLEX) Mie coeff. little-b-sub-n
    !      ACOEM1   (COMPLEX) Mie coeff. little-a-sub-(n-1)
    !      BCOEM1   (COMPLEX) Mie coeff. little-b-sub-(n-1)
    !
    !      LL       dimension of ACAP; in original program was 7000; for
    !               conserving memory this should not be much bigger than
    !                  1.1*Abs(RindSh) * x + 1
    !               BUT it must ALWAYS exceed 150
    !
    !      TA(1,2)  real, imaginary parts of WFN(1)
    !      TA(3,4)  real, imaginary parts of WFN(2)
    !
    !      S1,S2   complex scattering amplitudes; these are processed to
    !               produce elements of the real scattering matrix
    !               (they would need to be multiplied by a factor to
    !               give the true scattering amplitudes, VDH Sec 9.31)
    !
    !      TOLER   tolerance for cosines of angles when slightly below 0
    !               or slightly above 1 (a common occurrence)
    !
    !    NOTE:  the definitions of U(i) in this program are not the same as
    !           the u-sub-i defined in the Toon/Ackerman paper.  The
    !           correspondence is:
    !             usub1 = u(1)    usub2 = u(5)
    !             usub3 = u(7)    usub4 = dumsq
    !             usub5 = u(2)    usub6 = u(3)
    !             usub7 = u(6)    usub8 = u(4)
    !             ratio of spherical Bessel to spherical Hankel func = u(8)
    !
    !    The Bessel function ratio A is always computed by downward
    !    recurrence.
    !
    ! **********************************************************************
    ! Modif : C. Pinte
    ! 19/03/08 : conversion en fortran 90, suppresion goto et continue,
    ! ajout des intent, suppresion data, passage en module et interfacage
    ! mcfost
    ! 21/04/11 : skipping the computation of M1, M2, S21 and D21. Adjusting
    ! the angle values to match those used in mueller2/BHMIE.
    ! 13/10/13 : bug fix : avoiding computing the case 90 degrees twice
    !
    ! **********************************************************************

    ! .. Parameters ..
    !    integer, parameter ::   mxang = 100
    !  integer, parameter :: ll = 500000

    real(kind=dp), parameter :: zero = 0.0_dp
    real(kind=dp), parameter :: one = 1.0_dp
    real(kind=dp), parameter :: two = 2.0_dp
    real(kind=dp), parameter :: toler = 1.0e-6_dp

    ! .. Scalar Arguments ..
    integer, intent(in) ::   numang
    real, intent(in) :: rcore, rshell, wvno
    real, intent(out) ::   gqsc, qbs, qext, qsca
    complex, intent(in) ::   rindco, rindsh
    ! ..
    ! .. Array Arguments ..
    !    real, dimension(numang,2), intent(out) :: D21, M1 ,M2, S21
    complex, dimension(nang_scatt+1), intent(out) :: S1, S2
    ! ..
    ! .. Local Scalars ..
    logical ::   inperr
    integer ::   j, m, n, nmx1, nmx2, nn, jj, ll, alloc_status

    real(kind=dp) :: aa, aim, am1im, am1re, are, bb, bim, bm1im, bm1re, bre, cc, cosx1, cosx4
    real(kind=dp) :: dd, denom, dgqsc, dqext, dqsca, e2y1, ey1, ey1my4, ey1py4, ey4
    real(kind=dp) :: rmm, rx, sinx1, sinx4, x1, x4, xcore, xshell, y1, y4

    complex(kind=dp) :: ac, acoe, acoem1, bc, bcoe, bcoem1, dh1, dh2, dh4, dummy, dumsq
    complex(kind=dp) :: k1, k2, k3, p24h21, p24h24, rrfx, sback, wm1

    ! .. Local Arrays ..
    real(kind=dp), dimension(numang) ::   AMU
    real(kind=dp) :: THETA
    real(kind=dp), dimension(numang,3) ::  tau, pi_tab
    real(kind=dp), dimension(numang) :: SI2THT
    real(kind=dp), dimension(5) :: T
    real(kind=dp), dimension(5) :: TA
    real(kind=dp) :: DANG, PII

    complex(kind=dp), dimension(8) :: U
    complex(kind=dp), dimension(2) :: wfn
    complex(kind=dp), dimension(4) :: z

    complex(kind=dp), dimension(:), allocatable :: acap
    complex(kind=dp), dimension(:,:), allocatable :: W

    !    write(*,*) wvno,rcore,rshell,numang,rindco,rindsh

    ! ==============================
    ! Bits copied over from BHMIE
    !
    !*** Obtain pi:
    PII=4.*atan(1.D0)
    DANG=0.
    if (numang > 1) then
       DANG=.5*PII/real(numang-1,kind=dp)
    endif
    do J=1,numang
       THETA=(real(J,kind=dp)-1.0)*dang
       AMU(J)=cos(THETA)
    end do

    NN=2*numang-1
    do J=1,NN
       S1(J)=(0._dp,0._dp)
       S2(J)=(0._dp,0._dp)
    end do

    xshell = rshell*wvno
    xcore  = rcore*wvno
    T(1) = xshell*abs(rindsh)
    NMX1   = max(2.2*T(1),150.)
    NMX2   = NMX1 / 1.1

    ! Dynamical allocation
    ll = nmx1 + 1
    allocate(acap(ll), W(3,ll), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error in Dmilay'
       stop
    endif
    acap = czero ; W = czero ;

    ! Nothing changed below this line (except skipping computation of M1, M2, S21 and D21)
    ! ==============================

    ! ** Check input arguments for gross errors
    inperr = .false.
    if (wvno <= 0.0) inperr = wrtbad('WVNO')
    if (rshell <= 0.0) inperr = wrtbad('Rshell')
    if (rcore <= 0.0 .or. rcore > rshell)  inperr = wrtbad('Rcore')
    if (real(rindsh) <= 0.0 .or. aimag(rindsh) > 0.0) inperr = wrtbad('rindsh')
    if (real(rindco) <= 0.0 .or. aimag(rindco) > 0.0) inperr = wrtbad('RindCo')
    if (numang < 0) inperr = wrtbad('NumAng')
!    if (numang > mxang) inperr = wrtdim('MxAng', numang)
!    if (numang > maxang) inperr = wrtdim('MaxAng', numang)
!    if (nmx1 + 1 > ll) inperr = wrtdim('LL', nmx1 + 1)
    do  j=1, numang
       if (amu(j) < - toler .or. amu(j) > 1.0+toler)   inperr = wrtbad('MU')
    enddo
    if (inperr) call errmsg('MIELAY--Input argument errors.  Aborting...', .true.)

    K1    = RINDCO*WVNO
    K2   = RINDSH*WVNO
    K3   = DCMPLX(WVNO)
    Z(1) = RINDSH*XSHELL
    Z(2) = XSHELL
    Z(3) = RINDCO*XCORE
    Z(4) = RINDSH*XCORE
    X1   =  real(Z(1),kind=dp)
    Y1   = DIMAG(Z(1))
    X4   =  real(Z(4),kind=dp)
    Y4   = DIMAG(Z(4))
    RX   = ONE / XSHELL

    ! ** Down-recurrence for A function
    acap(NMX1 + 1) = czero
    do m = 1, 3
       w(m, nmx1 + 1) = czero
    enddo

    rrfx  = one / (rindsh*xshell)
    do nn = nmx1, 1, -1
       acap(nn) = ((nn+1)*rrfx) - one / (((nn+1)*rrfx) + acap(nn+1))
       do m = 1, 3
          w(m,nn) = ((nn + 1) / z(m+1)) - one / (((nn+1)/z(m+1)) + w(m,nn+1))
       enddo !m
    enddo !nn

    do  j=1, numang
       si2tht(j) = one - amu(j)**2
       pi_tab(j,1) = zero
       pi_tab(j,2) = one
       tau(j,1) = zero
       tau(j,2) = amu(j)
    enddo

    ! ** Initialization of homogeneous sphere
    T(1) = COS(XSHELL)
    T(2) = SIN(XSHELL)
    WM1      = DCMPLX(T(1), - T(2))
    WFN(1) = DCMPLX(T(2), T(1))
    TA(1) = T(2)
    TA(2) = T(1)
    WFN(2) = RX*WFN(1) - WM1
    TA(3) =  real(WFN(2),kind=dp)
    TA(4) = DIMAG(WFN(2))

    ! ** Initialization procedure for stratified sphere
    N      = 1
    SINX1  = SIN(X1)
    SINX4  = SIN(X4)
    COSX1  = COS(X1)
    COSX4  = COS(X4)
    EY1    = EXP(Y1)
    E2Y1   = EY1**2
    EY4    = EXP(Y4)
    EY1MY4 = EXP(Y1 - Y4)
    EY1PY4 = EY1*EY4
    AA     = SINX4*(EY1PY4 + EY1MY4)
    BB     = COSX4*(EY1PY4 - EY1MY4)
    CC     = SINX1*(E2Y1 + ONE)
    DD     = COSX1*(E2Y1 - ONE)
    DENOM  = ONE + E2Y1*(4.0D0*SINX1**2 - TWO + E2Y1)
    DUMMY  = DCMPLX((AA*CC + BB*DD) / DENOM, &
         (BB*CC - AA*DD) / DENOM)
    DUMMY  = DUMMY*(ACAP(N) + N / Z(1)) / (W(3, N) + N / Z(4))
    DUMSQ  = DUMMY**2

    P24H24 = 0.5D0 + DCMPLX(SINX4**2 - 0.5D0, COSX4*SINX4)*EY4**2
    P24H21 = 0.5D0*DCMPLX(SINX1*SINX4 - COSX1*COSX4, &
         SINX1*COSX4 + COSX1*SINX4)*EY1PY4 &
         + 0.5D0*DCMPLX(SINX1*SINX4 + COSX1*COSX4, &
         - SINX1*COSX4 + COSX1*SINX4)*EY1MY4
    DH1    = Z(1) / (ONE + CI*Z(1)) - ONE / Z(1)
    DH2    = Z(2) / (ONE + CI*Z(2)) - ONE / Z(2)
    DH4    = Z(4) / (ONE + CI*Z(4)) - ONE / Z(4)
    P24H24 = P24H24 / ((DH4 + N/Z(4))*(W(3, N) + N/Z(4)))
    P24H21 = P24H21 / ((DH1 + N/Z(1))*(W(3, N) + N/Z(4)))

    U(1) = K3*ACAP(N) - K2*W(1, N)
    U(2) = K3*ACAP(N) - K2*DH2
    U(3) = K2*ACAP(N) - K3*W(1, N)
    U(4) = K2*ACAP(N) - K3*DH2
    U(5) = K1*W(3, N) - K2*W(2, N)
    U(6) = K2*W(3, N) - K1*W(2, N)
    U(7) = - CI*(DUMMY*P24H21 - P24H24)
    U(8) = TA(3) / WFN(2)

    ACOE  = U(8)*(U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5)) / &
         (U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5))

    BCOE  = U(8)*(U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6)) / &
         (U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6))

    ACOEM1 = ACOE
    BCOEM1 = BCOE
    ARE    =  real(ACOE,kind=dp)
    AIM    = DIMAG(ACOE)
    BRE    =  real(BCOE,kind=dp)
    BIM    = DIMAG(BCOE)

    DQEXT  = 3.D0*(ARE + BRE)
    DQSCA  = 3.D0*(ARE**2 + AIM**2 + BRE**2 + BIM**2)
    DGQSC  = ZERO
    SBACK  = 3.D0*(ACOE - BCOE)
    RMM    = ONE

    AC  = 1.5D0*ACOE
    BC  = 1.5D0*BCOE
    do j = 1,NUMANG
       S1(J) = AC*PI_tab(J,2) + BC*TAU(J,2)
       S2(J) = BC*PI_tab(J,2) + AC*TAU(J,2)
    enddo
    do j = 1,NUMANG-1
       S1(2*numang-J) = AC*PI_tab(J,2) - BC*TAU(J,2)
       S2(2*numang-J) = BC*PI_tab(J,2) - AC*TAU(J,2)
    enddo

    ! ***************** Start of Mie summing loop ******************

    N = 2

    infinie : do
       ! ** Recurrences for functions little-pi,
       ! little-tau of Mie theory
       T(1) = 2*N - 1
       T(2) = N - 1

       do j=1, NUMANG
          PI_tab(J,3) = (T(1)*PI_tab(J,2)*AMU(J) - N*PI_tab(J,1)) / T(2)
          TAU(J,3) = AMU(J)*(PI_tab(J,3) - PI_tab(J,1)) - T(1)*SI2THT(J)*PI_tab(J,2) + TAU(J, 1)
       enddo

       ! ** Here set up homogeneous sphere
       WM1    = WFN(1)
       WFN(1) = WFN(2)
       WFN(2) = T(1)*RX*WFN(1) - WM1
       TA(1) =  real(WFN(1),kind=dp)
       TA(2) = DIMAG(WFN(1))
       TA(3) =  real(WFN(2),kind=dp)
       TA(4) = DIMAG(WFN(2))

       ! ** Here set up stratified sphere
       DH1    = - N / Z(1) + ONE / (N / Z(1) - DH1)
       DH2    = - N / Z(2) + ONE / (N / Z(2) - DH2)
       DH4    = - N / Z(4) + ONE / (N / Z(4) - DH4)
       P24H24 = P24H24 / ((DH4 + N/Z(4))*(W(3, N) + N/Z(4)))
       P24H21 = P24H21 / ((DH1 + N/Z(1))*(W(3, N) + N/Z(4)))
       DUMMY  = DUMMY*(ACAP(N) + N / Z(1)) / (W(3, N) + N / Z(4))
       DUMSQ  = DUMMY**2

       U(1) = K3*ACAP(N) - K2*W(1, N)
       U(2) = K3*ACAP(N) - K2*DH2
       U(3) = K2*ACAP(N) - K3*W(1, N)
       U(4) = K2*ACAP(N) - K3*DH2
       U(5) = K1*W(3, N) - K2*W(2, N)
       U(6) = K2*W(3, N) - K1*W(2, N)
       U(7) = - CI*(DUMMY*P24H21 - P24H24)
       U(8) = TA(3) / WFN(2)

       ACOE  = U(8)*(U(1)*U(5)*U(7) + K1*U(1) - DUMSQ*K3*U(5)) / (U(2)*U(5)*U(7) + K1*U(2) - DUMSQ*K3*U(5))

       BCOE  = U(8)*(U(3)*U(6)*U(7) + K2*U(3) - DUMSQ*K2*U(6)) / (U(4)*U(6)*U(7) + K2*U(4) - DUMSQ*K2*U(6))
       ARE  =  real(ACOE,kind=dp)
       AIM  = DIMAG(ACOE)
       BRE  =  real(BCOE,kind=dp)
       BIM  = DIMAG(BCOE)

       ! ** Increment sums for efficiency factors
       AM1RE  =  real(ACOEM1,kind=dp)
       AM1IM  = DIMAG(ACOEM1)
       BM1RE  =  real(BCOEM1,kind=dp)
       BM1IM  = DIMAG(BCOEM1)
       T(4) = (2*N - ONE) / (N*(N - ONE))
       T(2) = (N - ONE)*(N + ONE) / N
       DGQSC  = DGQSC + T(2)*(AM1RE*ARE + AM1IM*AIM + BM1RE*BRE + BM1IM*BIM) + T(4)*(AM1RE*BM1RE + AM1IM*BM1IM)

       T(3)  = 2*N + 1
       DQEXT   = DQEXT + T(3)*(ARE + BRE)
       T(4)  = ARE**2 + AIM**2 + BRE**2 + BIM**2
       DQSCA   = DQSCA + T(3)*T(4)
       RMM     = - RMM
       SBACK  = SBACK + T(3) * RMM *(ACOE - BCOE)

       T(2) = N*(N + 1)
       T(1) = T(3) / T(2)

       AC  = T(1)*ACOE
       BC  = T(1)*BCOE
       do j=1, NUMANG
          S1(J) = S1(J) + AC*PI_tab(J,3) + BC*TAU(J,3)
          S2(J) = S2(J) + BC*PI_tab(J,3) + AC*TAU(J,3)
       enddo

       ! ** Scattering matrix elements for
       ! supplements of 0-90 degree scattering
       ! angles submitted by user
       if(mod(N,2) == 0) then
          do j= 1, NUMANG-1
             JJ = 2*numang - j
             S1(JJ) = S1(JJ) - AC*PI_tab(J,3) + BC*TAU(J,3)
             S2(JJ) = S2(JJ) - BC*PI_tab(J,3) + AC*TAU(J,3)
          enddo
       else
          do j= 1, NUMANG-1
             JJ = 2*numang - j
             S1(JJ) = S1(JJ) + AC*PI_tab(J,3) - BC*TAU(J,3)
             S2(JJ) = S2(JJ) + BC*PI_tab(J,3) - AC*TAU(J,3)
          enddo
       endif

       ! ** Test for convergence of sums
       if (T(4) >= 1.0e-14_dp) then
          N  = N + 1
          if (N > NMX2) then
             write(*,*) 2*pii/wvno,rcore,rshell
             CALL ERRMSG('MIELAY--Dimensions for W,ACAP not enough. Suggest'// &
                  ' get detailed output, modify routine', .true.)
          endif
          do j = 1, NUMANG
             PI_tab(J,1) = PI_tab(J,2)
             PI_tab(J,2) = PI_tab(J,3)
             TAU(J,1) = TAU(J,2)
             TAU(J,2) = TAU(J,3)
          enddo

          ACOEM1 = ACOE
          BCOEM1 = BCOE

          ! On boucle

       else ! On sort de la boucle
          exit infinie
       endif

    enddo infinie

    ! ***************** End of summing loop ******************

    ! ** Transform complex scattering amplitudes
    ! into elements of real scattering matrix

!    do j = 1, NUMANG
!       do K = 1, 2
!          M1(J,K) = DBLE(S1(J,K))**2 + DIMAG(S1(J,K))**2
!          M2(J,K) = DBLE(S2(J,K))**2 + DIMAG(S2(J,K))**2
!          S21(J,K) = DBLE( S1(J,K))*DBLE( S2(J,K)) + DIMAG(S1(J,K))*DIMAG(S2(J,K))
!          D21(J,K) = DIMAG(S1(J,K))*DBLE(S2(J,K)) - DIMAG(S2(J,K))*DBLE(S1(J,K))
!       enddo
!    enddo


    T(1) = TWO*RX**2
    QEXT   = T(1)*DQEXT
    QSCA   = T(1)*DQSCA
    GQSC   = TWO*T(1)*DGQSC
    SBACK  = 0.5*SBACK
    QBS    = (real(SBACK,kind=dp)**2 + DIMAG(SBACK)**2) / (pi*XSHELL**2)

    deallocate(acap,W)

    return

  end subroutine dmilay

  ! **********************************************************************

  subroutine errmsg(messag, fatal)
    ! Print out a warning or error message;  abort if error
    ! after making symbolic dump (machine-specific)

    ! Provenance:  the 3 error-handlers ErrMsg, WrtBad, WrtDim are
    ! borrowed from MIEV, the Wiscombe Mie program
    ! .. Scalar Arguments ..

    CHARACTER, intent(in) :: MESSAG*(*)
    LOGICAL, intent(in) ::   FATAL

    if (fatal) then
       write(*, '(//,2A,//)') ' ****** ERROR *****  ', MESSAG
       stop
    endif

    NUMMSG = NUMMSG + 1

    if (msglim) return

    if (nummsg <= maxmsg) then
       write(*, '(/,2A,/)') ' ****** WARNING *****  ', MESSAG
    else
       write(*, '(//,A,//)') ' ****** TOO MANY WARNING MESSAGES --  ' // &
            'They will no longer be printed *******'
       msglim = .true.
    endif

  end subroutine errmsg

  ! **********************************************************************

  logical function wrtbad(varnam)
    ! Write names of erroneous variables and return 'TRUE'
    ! INPUT :   VarNam = Name of erroneous variable to be written
    ! (CHARACTER, any length)

    ! .. Scalar Arguments ..
    CHARACTER, intent(in) :: VARNAM*(*)

    WRTBAD = .TRUE.
    NUMMSG = NUMMSG + 1
    WRITE(*, '(3A)') ' ****  Input variable  ', VARNAM, '  in error  ****'

    IF(NUMMSG == MAXMSG) CALL ERRMSG('Too many input errors.  Aborting...', .TRUE.)

  end function wrtbad

  ! **********************************************************************

  logical function wrtdim(dimnam, minval)
    ! Write name of too-small symbolic dimension and
    ! the value it should be increased to;  return 'TRUE'

    ! INPUT :  DimNam = Name of symbolic dimension which is too small
    ! (CHARACTER, any length)
    ! Minval = Value to which that dimension should be
    ! increased (at least)

    ! .. Scalar Arguments ..
    CHARACTER, intent(in) :: DIMNAM*(*)
    INTEGER, intent(in) ::   MINVAL

    write(*, '(3A,I7)') ' ****  Symbolic dimension  ', DIMNAM, '  should be increased to at least ', MINVAL
    wrtdim = .true.

  end function wrtdim

end module coated_sphere
