module scattering

  use parameters
  use constants
  use wavelengths
  use grains
  use utils
  use read_opacity

  implicit none

  procedure(get_Mueller_matrix_per_grain), pointer :: get_Mueller_matrix => null()

  contains

subroutine setup_scattering()

  ! parameterization of scattering method according to grid size
  ! 1 : per dust grain
  ! 2 : per cell
  call select_scattering_method(p_n_cells)

  lMueller_pos_multi = .false.
  if (lmono) then
     p_n_lambda_pos = 1
  else
     if (scattering_method==1) then
        p_n_lambda_pos = 1
     else
        p_n_lambda_pos = n_lambda
        lMueller_pos_multi = .true.
     endif
  endif

end subroutine setup_scattering

!***************************************************

subroutine select_scattering_method(p_n_cells)

  integer, intent(in) :: p_n_cells

  real :: mem_size

  if (scattering_method0 == 0) then
     if (.not.lmono) then
        mem_size = (1.0*p_n_cells) * (nang_scatt+1) * n_lambda * 4 / 1024**3
        if (mem_size > max_mem) then
           scattering_method = 1
        else
           scattering_method = 2
        endif
     else
        if (lscatt_ray_tracing) then
           scattering_method = 2 ! it needs to be 2 for ray-tracing
        else
           ! ??? TODO + + TODO en realloc lscatt_ray_tracing = .false, en mode ML 3D ???
           scattering_method = 2
        endif
     endif
  endif

  write(*,fmt='(" Using scattering method ",i1)') scattering_method
  lscattering_method1 = (scattering_method==1)

end subroutine select_scattering_method

!***************************************************

subroutine bhmie(x,refrel,nang,s1,s2,qext,qsca,qback,gsca)

  implicit none

! Declare parameters:
! Note: important that MXNANG be consistent with dimension of S1 and S2
!       in calling routine!

  ! Useless since the switch to dynamic allocation
  !  integer, parameter :: NMXX=20000000
  ! defaut : NMXX=20000
  ! Can use NMXX=2000000: works up to a=10cm in band B
  ! but you shouldn.t be in a hurry
  integer, parameter :: dp = selected_real_kind(p=13,r=200)

! Arguments:
  integer, intent(in) :: nang
  real, intent(in) :: x
  complex, intent(in) :: refrel
  real, intent(out) :: gsca,qback,qext,qsca
  complex, dimension(2*nang-1), intent(out) :: s1,s2

! Local variables:
  integer :: J,JJ,N,NSTOP,NMX,NN
  real (kind =dp) :: CHI,CHI0,CHI1,DANG,DX,EN,FN,P,PII,PSI,PSI0,PSI1,THETA,XSTOP,YMOD
  real (kind =dp), dimension(NANG) :: AMU,PI,PI0,PI1,TAU
!  real (kind =dp) :: AMU_0,PI_0,PI0_0,PI1_0,TAU_0
  complex (kind=dp) :: AN,AN1,BN,BN1,DREFRL,XI,XI1,Y
!  complex (kind=dp), dimension(nmxx):: D
  complex (kind=dp), dimension(:), allocatable :: D
  integer :: alloc_status


!***********************************************************************
! Subroutine BHMIE is the Bohren-Huffman Mie scattering subroutine
!    to calculate scattering and absorption by a homogenous isotropic
!    sphere.
! Given:
!    X = 2*pi*a/lambda
!    REFREL = (complex refr. index of sphere)/(real index of medium)
!    NANG = number of angles between 0 and 90 degrees
!           (will calculate 2*NANG-1 directions from 0 to 180 deg.)
!           if called with NANG<2, will set NANG=2 and will compute
!           scattering for theta=0,90,180.
! Returns:
!    S1(1 - 2*NANG-1) = -i*f_22 (incid. E perp. to scatt. plane,
!                                scatt. E perp. to scatt. plane)
!    S2(1 - 2*NANG-1) = -i*f_11 (incid. E parr. to scatt. plane,
!                                scatt. E parr. to scatt. plane)
!    QEXT = C_ext/pi*a**2 = efficiency factor for extinction
!    QSCA = C_sca/pi*a**2 = efficiency factor for scattering
!    QBACK = (dC_sca/domega)/pi*a**2
!          = backscattering efficiency [NB: this is (1/4*pi) smaller
!            than the "radar backscattering efficiency"; see Bohren &
!            Huffman 1983 pp. 120-123]
!    GSCA = <cos(theta)> for scattering
!
! Original program taken from Bohren and Huffman (1983), Appendix A
! Modified by B.T.Draine, Princeton Univ. Obs., 90/10/26
! in order to compute <cos(theta)>
! 91/05/07 (BTD): Modified to allow NANG=1
! 91/08/15 (BTD): Corrected error (failure to initialize P)
! 91/08/15 (BTD): Modified to enhance vectorizability.
! 91/08/15 (BTD): Modified to make NANG=2 if called with NANG=1
! 91/08/15 (BTD): Changed definition of QBACK.
! 92/01/08 (BTD): Converted to full double precision and double complex
!                 eliminated 2 unneed lines of code
!                 eliminated redundant variables (e.g. APSI,APSI0)
!                 renamed RN -> EN = double precision N
!                 Note that DOUBLE COMPLEX and DCMPLX are not part
!                 of f77 standard, so this version may not be fully
!                 portable.  In event that portable version is
!                 needed, use src/bhmie_f77.f
! 93/06/01 (BTD): Changed AMAX1 to generic function MAX
! 04/03/04 (CP): switched to fortran 90
! 13/10/04 (CP): switched to half-integer angles (middle of the bin) to perform s11(theta) integration
! Angle 0 is computed separately as it is needed for Qext
! Similarly, 180 would be needed for Qback but is not used, so it is not computed
! 16/12/04 (CP): Remplacement imag par aimag (Pas forcement ok avec tous les
! compilateurs) mais standard f95 et 2003
! 24/03/05 (CP): dynamic allocation of D. Allows placing the array in the data zone
! and allocating it with exactly the correct number of terms.
! 07/02/08 (CP): switch back to integer angles with 0 and 180 explicitly calculated
! necessary for ray-tracing integration: switched to 2*nang + 1 angles
!***********************************************************************


!*** Safety checks
!  if(NANG > MXNANG) then
!     write(*,*)'***Error: NANG > MXNANG in bhmie'
!     stop
!  endif
!  if(NANG < 2) then
!     write(*,*)'***Error: NANG must be >= 2'
!     stop
!  endif
!*** Obtain pi:
  PII=4.*atan(1.D0)
  DX=X
  DREFRL=REFREL
  Y=X*DREFRL
  YMOD=abs(Y)
!
!*** Series expansion terminated after NSTOP terms
!    Logarithmic derivatives calculated from NMX on down

  XSTOP=X+4.*X**0.3333+2.
  NMX=max(XSTOP,YMOD)+15

  ! Allocation dynamique
  allocate(D(nmx), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error BHMIE')
  D = 0


! BTD experiment 91/1/15: add one more term to series and compare result
!      NMX=AMAX1(XSTOP,YMOD)+16
! test: compute 7001 wavelengths between .0001 and 1000 micron
! for a=1.0micron SiC grain.  When NMX increased by 1, only a single
! computed number changed (out of 4*7001) and it only changed by 1/8387
! conclusion: we are indeed retaining enough terms in series!
  NSTOP=XSTOP

!*** Require NANG.GE.1 in order to calculate scattering intensities
  DANG=0.
  if (NANG > 1) then
     DANG=.5*PII/real(NANG-1,kind=dp)
  endif
  do J=1,NANG
     THETA=(real(J,kind=dp)-1.0)*DANG
     AMU(J)=cos(THETA)
  end do

  do J=1,NANG
     PI0(J)=0.
     PI1(J)=1.
  end do

  NN=2*NANG-1
  do J=1,NN
     S1(J)=(0._dp,0._dp)
     S2(J)=(0._dp,0._dp)
  end do
!
!*** Logarithmic derivative D(J) calculated by downward recurrence
!    beginning with initial value (0.,0.) at J=NMX
!
  D(NMX)=(0.,0.)
  NN=NMX-1
  do N=1,NN
     EN=NMX-N+1
     D(NMX-N)=(EN/Y)-(1./(D(NMX-N+1)+EN/Y))
  end do
!
!*** Riccati-Bessel functions with real argument X
!    calculated by upward recurrence
!
  PSI0=cos(DX)
  PSI1=sin(DX)
  CHI0=-sin(DX)
  CHI1=cos(DX)
  XI1=CMPLX(PSI1,-CHI1,dp)
  QSCA=0.E0
  GSCA=0.E0
  P=-1.
  do N=1,NSTOP
     EN=N
     FN=(2.E0*EN+1.)/(EN*(EN+1.))
! for given N, PSI  = psi_n        CHI  = chi_n
!              PSI1 = psi_{n-1}    CHI1 = chi_{n-1}
!              PSI0 = psi_{n-2}    CHI0 = chi_{n-2}
! Calculate psi_n and chi_n
     PSI=(2.E0*EN-1.)*PSI1/DX-PSI0
     CHI=(2.E0*EN-1.)*CHI1/DX-CHI0
     XI=CMPLX(PSI,-CHI,dp)

!*** Store previous values of AN and BN for use
!    in computation of g=<cos(theta)>
     if(N > 1)then
        AN1=AN
        BN1=BN
     endif
!
!*** Compute AN and BN:
     AN=(D(N)/DREFRL+EN/DX)*PSI-PSI1
     AN=AN/((D(N)/DREFRL+EN/DX)*XI-XI1)
     BN=(DREFRL*D(N)+EN/DX)*PSI-PSI1
     BN=BN/((DREFRL*D(N)+EN/DX)*XI-XI1)

!*** Augment sums for Qsca and g=<cos(theta)>
     QSCA=QSCA+(2.*EN+1.)*(abs(AN)**2+abs(BN)**2)
     GSCA=GSCA+((2.*EN+1.)/(EN*(EN+1.)))*(real(AN)*real(BN)+aimag(AN)*aimag(BN))
     if(N > 1)then
        GSCA=GSCA+((EN-1.)*(EN+1.)/EN)*(real(AN1)*real(AN)+aimag(AN1)*aimag(AN) &
             +real(BN1)*real(BN)+aimag(BN1)*aimag(BN))
     endif
!
!*** Now calculate scattering intensity pattern
!    First do angles from 0 to 90
     do  J=1,NANG
        PI(J)=PI1(J)
        TAU(J)=EN*AMU(J)*PI(J)-(EN+1.)*PI0(J)
        S1(J)=S1(J)+FN*(AN*PI(J)+BN*TAU(J))
        S2(J)=S2(J)+FN*(AN*TAU(J)+BN*PI(J))
     enddo
!
!*** Now do angles greater than 90 using PI and TAU from
!    angles less than 90.
!    P=1 for N=1,3,...; P=-1 for N=2,4,...
     P=-P
     do J=1,NANG-1
        JJ=2*NANG-J
        S1(JJ)=S1(JJ)+FN*P*(AN*PI(J)-BN*TAU(J))
        S2(JJ)=S2(JJ)+FN*P*(BN*PI(J)-AN*TAU(J))
     enddo
     PSI0=PSI1
     PSI1=PSI
     CHI0=CHI1
     CHI1=CHI
     XI1=CMPLX(PSI1,-CHI1,dp)
!
!*** Compute pi_n for next value of n
!    For each angle J, compute pi_n+1
!    from PI = pi_n , PI0 = pi_n-1
     do J=1,NANG
        PI1(J)=((2.*EN+1.)*AMU(J)*PI(J)-(EN+1.)*PI0(J))/EN
        PI0(J)=PI(J)
     enddo
  end do
!
!*** Have summed sufficient terms.
!    Now compute QSCA,QEXT,QBACK,and GSCA
  GSCA=2.*GSCA/QSCA
  QSCA=(2./(DX*DX))*QSCA
  QEXT=(4./(DX*DX))*real(S1(1))
  QBACK=(abs(S1(2*NANG-1))/DX)**2/PII

  deallocate(D)
  return

end subroutine BHMIE

!***************************************************

subroutine mueller_Mie(lambda,igrain,x,amu1,amu2, qext,qsca,gsca)
  !***************************************************************
  ! calculate les elements de la matrice de diffusion a partir de
  ! la sous-routine bhmie (grains spheriques)
  !
  ! C. Pinte Fevrier 2004
  !****************************************************************

  implicit none
  integer, intent(in) :: lambda, igrain
  real, intent(in) :: amu1, amu2
  real, intent(in) :: x ! 2*pi*a/wl
  real, intent(out) :: qext, qsca, gsca

  integer :: j, jj, nang

  complex, dimension(nang_scatt+1) :: S1,S2

  real :: vi1, vi2, qback, factor
  complex :: s, refrel
  real, dimension(0:nang_scatt) ::  s11,s12,s22,s33,s34,s44

  refrel = cmplx(amu1,amu2)

  if (modulo(nang_scatt,2) == 1) call error("nang_scatt must be an EVEN number")

  ! If HG function, the phase function is not calculated
  if (aniso_method==2) then
     nang=1
  else
     nang= (nang_scatt+1) / 2 + 1
  endif

  call bhmie(x,refrel,nang, s1,s2,qext,qsca,qback,gsca)

  if (lforce_HG) gsca = forced_g
  if (lisotropic) gsca = 0.0

  ! The Normalization of bhmie for s11 is 0.5*x**2*Qsca
  ! --> We correct by 0.5*x**2 to have a Normalization to Qsca Qsa
  factor = 1 / (0.5 * x**2)

  ! Transfer of values to mcfost arrays
  if (aniso_method==1) then
     ! Calculation of scattering matrix elements
     ! indices offset by 1 compared to bhmie
     do j=0,nang_scatt
        jj=j+1
        vi1 = cabs(S2(jj))*cabs(S2(jj))
        vi2 = cabs(S1(jj))*cabs(S1(jj))
        s11(j) = 0.5*(vi1 + vi2) * factor
        s12(j) = 0.5*(vi1 - vi2) * factor
        s = S2(jj)*conjg(S1(jj))
        s33(j)=real(s) * factor
        s34(j)=aimag(s) * factor
     enddo !j

     s22(:) = s11(:)
     s44(:) = s33(:)

     call normalise_Mueller_matrix(lambda,igrain, s11,s12,s22,s33,s34,s44, normalization=qsca)
  endif

  return

end subroutine mueller_Mie

!***************************************************

subroutine Mueller_input(lambda, k_abs,k_sca,gsca)

  integer, intent(in) :: lambda
  real, intent(out) :: k_abs, k_sca, gsca

  character(len=128) :: string, filename

  integer :: EoF, alloc_status,iformat,nlam,iang,nang, l
  logical :: lscat
  real(dp), dimension(:), allocatable :: angles,wavel,kabs,ksca,g
  real(dp), dimension(:,:), allocatable :: f11,f12,f22,f33,f34,f44
  real(dp) :: wl, frac, frac_m1

  real, dimension(0:nang_scatt) ::  s11,s12,s22,s33,s34,s44

  gsca = 0

  alloc_status = 0
  EOF=0
  open(unit=1,file=filename,status='old')

  ! Skip lines starting with "#"
  do while(EoF==0)
     read(1,*,iostat=EoF) string
     if (string(1:1) /= '#') exit
  enddo

  ! Read the format number.
  read(string,*) iformat
  if(iformat == 3) then
     lscat = .false.
  else if(iformat == 1) then
     lscat = .true.
  endif

  ! Read the number of wavelengths
  read(1,*) nlam

  ! Read the number of scattering angles
  if (lscat) read(1,*) nang

  ! Allocate arrays
  alloc_status = 0
  allocate(wavel(nlam),kabs(nlam),ksca(nlam),g(nlam), stat = alloc_status)

  ! Read wavelengths, opacities and g
  do l=1,nlam
     read(1,*) wavel(l), kabs(l), ksca(l), g(l)
  enddo

  wl = tab_lambda(lambda)
  if (wl < minval(wavel)) call error("Mueller matrices: wavelength is smaller than tabulated")
  if (wl > maxval(wavel)) call error("Mueller matrices: wavelength is larger than tabulated")

  if (wavel(2) >  wavel(1)) then ! increasing order
     do l=1,nlam
        if  (wavel(l) > wl) then
           frac = (wavel(l) - wl) / (wavel(l) - wavel(l-1))
           exit
        endif
     enddo
  else ! decreasing order
     do l=1,nlam
        if  (wavel(l) < wl) then
           frac = (wl - wavel(l)) / (wavel(l-1) - wavel(l))
           exit
        endif
     enddo
  endif
  if ((frac > 1) .or. (frac < 0)) call error("Mueller matrices: interpolation error")
  frac_m1 = 1.0_dp - frac

  if (lscat) then
     allocate(angles(nang), f11(nlam,nang),f12(nlam,nang),f22(nlam,nang),&
          f33(nlam,nang),f34(nlam,nang),f44(nlam,nang), stat = alloc_status)

     ! Read the angular grid
     do iang=1,nang
        read(1,*) angles(iang)
     enddo

     ! Read the scattering matrix.
     do l=1,nlam
        do iang=1,nang
           read(1,*) f11(l,iang),f12(l,iang), f22(l,iang), f33(l,iang), f34(l,iang), f44(l,iang)
        enddo
     enddo
  endif
  close(1)

  ! Interpolate at correct wavelength
  ! todo : find l
  call error("Mueller input: interpolation was never implemented")

  k_abs = kabs(l-1) * frac + kabs(l) * frac_m1
  k_sca = ksca(l-1) * frac + ksca(l) * frac_m1
  gsca = ksca(l-1) * frac + ksca(l) * frac_m1

  s11(:) = f11(l-1,:) * frac + f11(l,:) * frac_m1
  s12(:) = f11(l-1,:) * frac + f12(l,:) * frac_m1
  s22(:) = f11(l-1,:) * frac + f22(l,:) * frac_m1
  s33(:) = f11(l-1,:) * frac + f33(l,:) * frac_m1
  s34(:) = f11(l-1,:) * frac + f34(l,:) * frac_m1
  s44(:) = f11(l-1,:) * frac + f44(l,:) * frac_m1

  ! Deallocate arrays
  deallocate(wavel,kabs,ksca,g)
  if (lscat) then
     deallocate(angles)
     deallocate(f11,f12,f22,f33,f34,f44)
  endif

  return

end subroutine Mueller_input

!***************************************************

subroutine normalise_Mueller_matrix(lambda,igrain, s11,s12,s22,s33,s34,s44, normalization)
  ! read in the matrix elements, calculate the cumulative s11, niormnaklise and fill up tab_s11, etc

  integer, intent(in) :: lambda,igrain
  real, dimension(0:nang_scatt), intent(in) ::  s11,s12,s22,s33,s34,s44
  real, intent(in), optional :: normalization
  ! "normalization" is the integral of S11(theta) sin(theta) dtheta if known
  ! The missing "flux" from the numerical integration is placed between 0 and 1 degree (or first scattering angle)
  ! as we assume it is mostly diffraction that has been missed by the sampling

  integer :: j
  real(kind=dp) :: norm, theta, dtheta

  ! S11 integration for angle selection
  if (scattering_method==1) then
     prob_s11(lambda,igrain,0)=0.0
     dtheta = pi/real(nang_scatt)

     do j=2,nang_scatt ! scattering probability up to angle j, j=0 is skipped because sin(theta) = 0
        theta = real(j)*dtheta
        prob_s11(lambda,igrain,j)=prob_s11(lambda,igrain,j-1)+s11(j)*sin(theta)*dtheta
     enddo

     ! there is a numerical issue when x >> 1 because the angular resolution is not sufficient
     ! The diffraction peak is missed (particularly between 0 and 1)
     if (present(normalization)) then
        if (normalization > prob_s11(lambda,igrain,nang_scatt)) then
           prob_s11(lambda,igrain,1:nang_scatt) = prob_s11(lambda,igrain,1:nang_scatt) + &
                normalization - prob_s11(lambda,igrain,nang_scatt)
        else
           call error("normalise_Mueller_matrix: exact normalization is smaller than numerical integration")
        endif
     endif

     ! Normalization of the cumulative probability to 1
     prob_s11(lambda,igrain,:)=prob_s11(lambda,igrain,:)/prob_s11(lambda,igrain,nang_scatt)
  endif ! scattering_method==1

  do j=0,nang_scatt
     if (scattering_method==1) then ! Mueller matrix per grain
        ! Normalization for scattering according to phase function (tab_s11=1.0 used in Stokes)
        norm = 1./s11(j)
     else ! Otherwise normalization to Qsca
        norm = 1.
     endif

     tab_s11(j,igrain,lambda) = s11(j) * norm
     tab_s12(j,igrain,lambda) = s12(j) * norm
     tab_s22(j,igrain,lambda) = s22(j) * norm
     tab_s33(j,igrain,lambda) = s33(j) * norm
     tab_s34(j,igrain,lambda) = s34(j) * norm
     tab_s44(j,igrain,lambda) = s44(j) * norm
  enddo

  return

end subroutine normalise_Mueller_matrix

!***************************************************

subroutine overwrite_s12(Pmax)

  real, intent(in) :: Pmax

  real :: dtheta, theta
  integer :: j

  dtheta = pi/real(nang_scatt)
  do j=0,nang_scatt
     theta = real(j)*dtheta
     tab_s12(j,:,:) = - Pmax * (1-(cos(theta))**2)
  enddo

  return

end subroutine overwrite_s12

!***************************************************

subroutine mueller_GMM(lambda,igrain, qext,qsca,gsca)
!***************************************************************
! calculate les elements de la matrice de diffusion a partir du
! code gmm01TrA (clusters de sphres)
!     Aggrgats
!
!        ALSO calculate "G" = THE ASYMMETRY PARAMETER
!
! C. Pinte
! 04/07/2005
!****************************************************************

  implicit none

  integer, intent(in) :: lambda, igrain
  real, intent(out) :: qext, qsca

  integer, parameter :: nang2 = 2*nang_scatt+1

  integer :: j, i

  real :: gsca, norm, somme_sin, somme_prob, somme2
  real :: cext,cabs,csca,cbak,cpr,assym
  real :: cextv,cabsv,cscav,cbakv,cprv
  real :: cexts,cabss,cscas,cbaks,cprs

  real, dimension(nang2) :: dang
  real, dimension(4,4,nang2) :: mue

  real, dimension(4,4,2*nang_scatt) :: mueller

  integer :: idscmt=1

  character(len=128) :: string


  ! Correct for the new definition of nang_scatt
  ! Also correct the normalization of s11 (do like mueller_Mie)
  ! int S11 sin(theta) dtheta = Qsca, use exact normalization, not numerical
  call error("mueller_gmm needs to be updated")

  if(n_grains_tot > 1) call error("You must choose n_grains_tot=1")
  if (scattering_method /= 1) call error("You must choose scattering_method 1")

  ! Lecture du fichier de rsultats de gmm01TrA : 'gmm01TrA.out'
  open(unit=12,file=mueller_aggregate_file,status='old')
  read(12,*)
  read(12,*)
  if(idscmt < 0) then
     read(12,*) string
     read(12,*)
  endif
  read(12,*)
  read(12,*) string
  read(12,*) cext,cabs,csca,cbak,cpr,assym
  cpr=cext-cpr
  read(12,*) string
  read(12,*) cextv,cabsv,cscav,cbakv,cprv,assym
  cprv=cextv-cprv
  read(12,*) string
  read(12,*) cexts,cabss,cscas,cbaks,cprs,assym
  cprs=cexts-cprs
  if(idscmt > 0) then
     read(12,*)
     read(12,*)
     read(12,*) string
     do i=1,nang2
        !read(12,'(f6.1,e13.5,f8.4,4e13.5)') dang(i),inat(i),pol(i),i11(i),i21(i),i12(i),i22(i)
        read(12,*)
     enddo
     read(12,*)
     read(12,*)  !'Scattering matrix (4X4 for each scattering angle):'
     read(12,*) string
     do i=1,nang2
        read(12,'(f7.1,4e16.7)') dang(i),mue(1,1,i),mue(1,2,i),mue(1,3,i),mue(1,4,i)
        read(12,'(7x,4e16.7)')   mue(2,1,i),mue(2,2,i),mue(2,3,i),mue(2,4,i)
        read(12,'(7x,4e16.7)')   mue(3,1,i),mue(3,2,i),mue(3,3,i),mue(3,4,i)
        read(12,'(7x,4e16.7)')   mue(4,1,i),mue(4,2,i),mue(4,3,i),mue(4,4,i)
     enddo
  endif
  close(12)
  ! End of reading
  close(unit=1)


!  QABS=QEXT-QSCA
! Calculation of scattering matrix elements
! Calcul angle central du bin par interpolation linaire
  do J=1,2*NANG_scatt
     mueller(:,:,j) = 0.5*(mue(:,:,j)+mue(:,:,j+1))
  enddo !j

  ! S11 integration for angle selection
  prob_s11(lambda,igrain,0)=0.0
  somme_sin= 0.0
  somme2 = 0.0

  do j=1,2*nang_scatt
     prob_s11(lambda,igrain,j)=prob_s11(lambda,igrain,j-1)+&
          mueller(1,1,j)*sin((real(j)-0.5)/180.*pi)*pi/(2*nang_scatt)
     somme_sin = somme_sin + sin((real(j)-0.5)/180.*pi)*pi/(2*nang_scatt)
!     somme2=somme2+s12(j)*sin((real(j)-0.5)/180.*pi)*pi/(2*nang)
! Somme2 sert juste pour faire des plots
  enddo

  ! normalization
  somme_prob=prob_s11(lambda,igrain,2*nang_scatt) ! = (0.5*x**2*qsca)
  ! Soit int_0^\pi (i1(t)+i2(t)) sin(t) = x**2*qsca
  do j=1,2*nang_scatt
     prob_s11(lambda,igrain,j)=prob_s11(lambda,igrain,j)/somme_prob
  enddo


  do J=1,2*NANG_scatt
!     ! Normalization pour diffusion isotrope et E_sca(theta)
!     if (j == 1)  then
!        norm = somme_prob/somme_sin
!     endif

! Normalization REMOVED FOR TAB_POS CALCULATIONS (MATRICES DE MUELLER
! PAR cell)
! TO BE REINSTATED FOR MUELLER MATRICES PER GRAIN

     if (scattering_method==1) then
        ! Normalization for scattering according to phase function (tab_s11=1.0 used in Stokes)
        norm=mueller(1,1,j) !* qext/q sca
        mueller(:,:,j) = mueller(:,:,j) / norm
     endif

     tab_s11(j,igrain,lambda) = mueller(1,1,j)
     tab_s12(j,igrain,lambda) = mueller(1,2,j)
     tab_s22(j,igrain,lambda) = mueller(2,2,j)
     tab_s33(j,igrain,lambda) = mueller(3,3,j)
     tab_s34(j,igrain,lambda) = mueller(3,4,j)
     tab_s44(j,igrain,lambda) = mueller(4,4,j)
  enddo

  gsca = assym
  qext=cext/(pi*R_sph_same_M)
  qsca=csca/(pi*R_sph_same_M)

  return

end subroutine mueller_GMM

!***************************************************

subroutine Fresnel_input(lambda,igrain, qext,qsca,gsca)

!***************************************************************
! Routine driv de Mueller_GMM.
! calculate les elements de la matrice de diffusion a partir d'un
! fichier ascii contenant divers informations :
! La matrice de Mueller  chaque angle de diffraction de 0  180,
! Qsca, Qext et le factor d'assymtrie
! Here is an example file:
!        Qext        Qsca        <cos(theta)>
!     0.34402E+01  0.82702E+00     0.73069E+00
!
!
!                          Mueller Scattering Matrix
!    0.0   0.7171184E+04   0.1608979E+03  -0.1360425E-13 0.1008702E-12
!          0.1608979E+03   0.7171184E+04   0.1085675E-13 -0.1342563E-13
!         -0.1281373E-13  -0.1100342E-13   0.7169003E+04 0.7345735E+02
!          0.1013335E-12   0.1580635E-13  -0.7345735E+02 0.7169003E+04
!    1.0   0.6996034E+04   0.1585395E+03  -0.1087324E-12 -0.2014537E-12
!          0.1585395E+03   0.6996034E+04   0.7065835E-13 -0.1650172E-12
!.......
!
!
! F. Malaval
! 20/04/2023
!****************************************************************

  implicit none

  integer, intent(in) :: lambda, igrain
  real, intent(out) :: qext, qsca, gsca

  integer :: j, i

  real :: norm, somme_prob, theta, dtheta
  real :: ext,sca,assym

  real, dimension(0:nang_scatt) :: dang

  real, dimension(4,4,0:nang_scatt) :: mueller

  character(len=128) :: string

  if (modulo(nang_scatt,2) == 1) call error("nang_scatt must be an EVEN number")

  if (aniso_method == 2) call error ("You shoudn't use a hg function option when putting Mueller matrix in input")

  ! Lecture du fichier d'entre :
  open(unit=12,file=mueller_file,status='old')
  read(12,*) string
  read(12,*) ext,sca,assym
  read(12,*)
  read(12,*)  !'Matrice de mueller (4X4 coefficients pour chaque angle):'
  read(12,*) string
  do i=0,nang_scatt
     read(12,*) dang(i),mueller(1,1,i),mueller(1,2,i),mueller(1,3,i),mueller(1,4,i)
     read(12,*)         mueller(2,1,i),mueller(2,2,i),mueller(2,3,i),mueller(2,4,i)
     read(12,*)         mueller(3,1,i),mueller(3,2,i),mueller(3,3,i),mueller(3,4,i)
     read(12,*)         mueller(4,1,i),mueller(4,2,i),mueller(4,3,i),mueller(4,4,i)
  enddo
  close(12)
  ! End of reading
  close(unit=1)

  if (scattering_method == 1) then
     ! Integration de S11 pour tirer angle
     ! Initial condition: sin(0)=0
     prob_s11(lambda,igrain,0)=0.0
     dtheta = pi/real(nang_scatt)

     do j=1,nang_scatt
        theta = real(j)*dtheta
        prob_s11(lambda, igrain,j)=prob_s11(lambda, igrain,j-1)+&
             mueller(1,1,j)*sin(theta)*dtheta
     enddo

     ! normalization
     somme_prob = prob_s11(lambda, igrain,nang_scatt)
     prob_s11(lambda, igrain,:)=prob_s11(lambda, igrain,:)/somme_prob
  else  ! non-normalized values will be needed
     do i=1, 4
        do j=1, 4
           if (i*j>1) mueller(i,j,:) = mueller(i,j,:)*mueller(1,1,:)
        enddo
     enddo
     ! now the integral of s11 is calculated
     somme_prob=0.0
     dtheta = pi/real(nang_scatt)
     do j=1,nang_scatt
        theta = real(j)*dtheta
        somme_prob = somme_prob + mueller(1,1,j)*sin(theta)*dtheta
     enddo
  endif ! scatt_meth

  do J=0,NANG_scatt
     if (scattering_method==1) then ! Mueller matrix per grain
        ! Normalization for scattering according to phase function (tab_s11=1.0 used in Stokes)
        ! only S11 is normalized because the data are already normalized
        norm=mueller(1,1,j)
        if (norm > tiny_real) then
           mueller(1,1,j) = mueller(1,1,j) / norm
        else
           write (*,*) "at angle", real(j)*pi/real(nang_scatt)
           call error ("s11=0.0")
        endif
     else ! Otherwise normalization to Qsca. somme_prob is used.
        norm = somme_prob/sca
        if (norm > tiny_real) then
           mueller(:,:,j) = mueller(:,:,j)/norm
        else
           call error ("Error : s11 is normalised to 0")
        endif
     endif

     tab_s11(j,igrain,lambda) = mueller(1,1,j)
     tab_s12(j,igrain,lambda) = mueller(1,2,j)
     tab_s22(j,igrain,lambda) = mueller(2,2,j)
     tab_s33(j,igrain,lambda) = mueller(3,3,j)
     tab_s34(j,igrain,lambda) = mueller(3,4,j)
     tab_s44(j,igrain,lambda) = mueller(4,4,j)
  enddo
  qext = ext
  qsca = sca
  gsca = assym
  if (lforce_HG)  gsca = forced_g
  if (lisotropic)  gsca = 0.0

  return

end subroutine Fresnel_input

!***************************************************

 subroutine Fresnel_input_size(lambda,igrain, qext,qsca,gsca)

!***************************************************************
! Routine derived from Mueller_GMM.
! calculate the scattering matrix elements from an
! ascii file per grain size containing various information:
! The Mueller matrix at each diffraction angle from 0 to 180
! Qsca, Qext, the asymmetry factor
! Here is an example file:
!        Qext        Qsca        <cos(theta)>
!     0.34402E+01  0.82702E+00     0.73069E+00
!
!
!                          Mueller Scattering Matrix
!    0.0   0.7171184E+04   0.1608979E+03  -0.1360425E-13 0.1008702E-12
!          0.1608979E+03   0.7171184E+04   0.1085675E-13 -0.1342563E-13
!         -0.1281373E-13  -0.1100342E-13   0.7169003E+04 0.7345735E+02
!          0.1013335E-12   0.1580635E-13  -0.7345735E+02 0.7169003E+04
!    1.0   0.6996034E+04   0.1585395E+03  -0.1087324E-12 -0.2014537E-12
!          0.1585395E+03   0.6996034E+04   0.7065835E-13 -0.1650172E-12
!.......
!
! The address of these files (sorted in increasing order according to grain size)
! is contained in a text file.
!
! F. Malaval
! 20/04/2023
!****************************************************************

  implicit none

  integer, intent(in) :: lambda, igrain
  real, intent(out) :: qext, qsca, gsca

  integer :: j, i

  real :: norm, somme_prob, theta, dtheta
  real :: ext,sca,assym, size

  real, dimension(0:nang_scatt) :: dang

  real, dimension(4,4,0:nang_scatt) :: mueller

  character(len=128) :: string, path

  if (modulo(nang_scatt,2) == 1) call error("nang_scatt must be an EVEN number")

  if (aniso_method == 2) call error ("You shouldn't use a hg function option when putting Mueller matrix in input")

  open(unit=13, file = mueller_file, status = 'old')
  ! In this case, mueller_file contains the paths to the files for each matrix,
  ! sorted in increasing order according to grain size.
  do j=1, igrain
     read(13, *) size, path
  enddo
  close(unit=13)
  if (abs(size-r_grain(igrain))>r_grain(igrain)*0.00001) then
     write(*,*) "error, grain size in file is", size, "while expecting", r_grain(igrain)
     call error("Grain sizes do not match")
  endif

  ! Reading the input file:
  open(unit=12,file=path,status='old')
  read(12,*) string
  read(12,*) ext,sca,assym
  read(12,*)
  read(12,*)  !'Mueller matrix (4X4 coefficients for each angle):'
  read(12,*) string
  do i=0,nang_scatt
     read(12,*) dang(i),mueller(1,1,i),mueller(1,2,i),mueller(1,3,i),mueller(1,4,i)
     read(12,*)   mueller(2,1,i),mueller(2,2,i),mueller(2,3,i),mueller(2,4,i)
     read(12,*)   mueller(3,1,i),mueller(3,2,i),mueller(3,3,i),mueller(3,4,i)
     read(12,*)   mueller(4,1,i),mueller(4,2,i),mueller(4,3,i),mueller(4,4,i)
  enddo
  close(12)
  ! End of reading

  close(unit=1)

  if (scattering_method == 1) then
     ! S11 integration for angle selection
     ! Initial condition: sin(0)=0
     prob_s11(lambda,igrain,0)=0.0
     dtheta = pi/real(nang_scatt)

     do j=1,nang_scatt
        theta = real(j)*dtheta
        prob_s11(lambda, igrain,j)=prob_s11(lambda, igrain,j-1)+&
             mueller(1,1,j)*sin(theta)*dtheta
     enddo

     ! Normalization
     somme_prob= prob_s11(lambda, igrain,nang_scatt)
     prob_s11(lambda, igrain,:)=prob_s11(lambda, igrain,:)/somme_prob
  else  ! non-normalized values will be needed
     do i=1, 4
        do j=1, 4
           if (i*j>1) mueller(i,j,:) = mueller(i,j,:)*mueller(1,1,:)
        enddo
     enddo
     ! now the integral of s11 is calculated
     somme_prob=0.0
     dtheta = pi/real(nang_scatt)
     do j=1,nang_scatt
        theta = real(j)*dtheta
        somme_prob = somme_prob + mueller(1,1,j)*sin(theta)*dtheta
     enddo
  endif ! scatt_meth


  do J=0,NANG_scatt
     if (scattering_method==1) then ! Mueller matrix per grain
        ! Normalization for scattering according to phase function (tab_s11=1.0 used in Stokes)
        ! only S11 is normalized because the data are already normalized
        norm=mueller(1,1,j)
        if (norm > tiny_real) then
           mueller(1,1,j) = mueller(1,1,j) / norm
        else
           write (*,*) "at angle", real(j)*pi/real(nang_scatt)
           call error ("s11=0.0")
        endif
     else ! Otherwise normalization to Qsca. somme_prob is used.
        norm = somme_prob/sca
        if (norm > tiny_real) then
           mueller(:,:,j) = mueller(:,:,j)/norm
        else
           call error ("Error : s11 is normalised to 0")
        endif
     endif

     tab_s11(j,igrain,lambda) = mueller(1,1,j)
     tab_s12(j,igrain,lambda) = mueller(1,2,j)
     tab_s22(j,igrain,lambda) = mueller(2,2,j)
     tab_s33(j,igrain,lambda) = mueller(3,3,j)
     tab_s34(j,igrain,lambda) = mueller(3,4,j)
     tab_s44(j,igrain,lambda) = mueller(4,4,j)
  enddo
  qext = ext
  qsca = sca
  gsca = assym
  if (lforce_HG)  gsca = forced_g
  if (lisotropic)  gsca = 0.0

  return

end subroutine Fresnel_input_size

!***************************************************

subroutine mueller_opacity_file(lambda,igrain, qext,qsca,gsca)
  ! bi-linear interpolation (in log-log) of cross sections
  ! for grains after reading the opacity file
  ! Particularly for Draine's PAHs
  ! Assumes a HG for the phase function and zero polarizability !!
  ! C. Pinte
  ! 31/01/07

  implicit none

  integer, intent(in) :: igrain, lambda
  real, intent(out) :: qext,qsca,gsca

  real :: frac_a, frac_a_m1, frac_lambda, fact1, fact2, fact3, fact4
  integer :: i, j, pop, N

  real, dimension(0:nang_scatt) ::  s11,s12,s22,s33,s34,s44
  real :: norm, somme_prob, log_a, log_wavel, wl_min, wl_max, theta, dtheta

  log_a=log(r_grain(igrain))
  log_wavel = log(tab_lambda(lambda))

  pop=grain(igrain)%pop

  if ((igrain == dust_pop(pop)%ind_debut).and.(lambda==1)) call read_opacity_file(pop)

  ! Increasing order for grain sizes
  if (r_grain(igrain) < exp(op_file_log_r_grain(1,pop))) then
     if (lambda==1) then
        write(*,*) "WARNING: index=",igrain, "grain size=",r_grain(igrain)
        write(*,*) "Minimum grain size in opacity file is",  exp(op_file_log_r_grain(1,pop))
        write(*,*) "Smaller grains are assumed to have the same opacity"
     endif
    j = 2
    frac_a = 0.0 ; frac_a_m1 = 1.0
  else if (r_grain(igrain) > exp(op_file_log_r_grain(op_file_na(pop),pop))) then
     if (lambda==1) then
        write(*,*) "WARNING: index=",igrain, "grain size=",r_grain(igrain)
        write(*,*) "Maximum grain size in opacity file is",  exp(op_file_log_r_grain(op_file_na,pop))
        write(*,*) "Larger grains are assumed to have the same opacity"
     endif
     j = op_file_na(pop)
     frac_a = 1.0 ; frac_a_m1 = 0.
  else
     ! Search in grain size
     ! Increasing array
     do j=2,op_file_na(pop)
        if (op_file_log_r_grain(j,pop) > log_a) exit
     enddo
     frac_a = (log_a-op_file_log_r_grain(j-1,pop))/(op_file_log_r_grain(j,pop)-op_file_log_r_grain(j-1,pop))
     frac_a_m1 = 1 - frac_a
  endif

  ! Wavelength averaging
  wl_min = tab_lambda_inf(lambda)
  wl_max = tab_lambda_sup(lambda)

  qext = 0.0 ; qsca = 0.0 ; gsca = 0.0 ; norm = 0 ; N=0
  do i=1,op_file_n_lambda(pop)
     if ((op_file_lambda(i,pop) > wl_min).and.(op_file_lambda(i,pop) < wl_max)) then
        N = N+1
        norm = norm + op_file_delta_lambda(i,pop)
        qext = qext + (frac_a * op_file_Qext(i,j,pop) +  frac_a_m1 * op_file_Qext(i,j-1,pop) ) * op_file_delta_lambda(i,pop)
        qsca = qsca + (frac_a * op_file_Qsca(i,j,pop) +  frac_a_m1 * op_file_Qsca(i,j-1,pop) ) * op_file_delta_lambda(i,pop)
        gsca = gsca + (frac_a * op_file_g(i,j,pop) +  frac_a_m1 * op_file_g(i,j-1,pop) ) * op_file_delta_lambda(i,pop)
     endif
  enddo !i
  !write(*,*) lambda, N, norm

  if (norm > 0) then ! averaging over selected points
     qext = qext / norm
     qsca = qsca / norm
     gsca = gsca / norm
  else ! cannot average, doing log interpolation
     ! Wavelength search
     if (op_file_lambda(2,pop) > op_file_lambda(1,pop)) then  ! Increasing order
        do i=2,op_file_n_lambda(pop)
           if (log(op_file_lambda(i,pop)) > log_wavel) exit
        enddo
        !log(op_file_lambda(i)) > log_wavel >  log(op_file_lambda(i-1))
        frac_lambda = (log_wavel-log(op_file_lambda(i-1,pop)))/(log(op_file_lambda(i,pop))-log(op_file_lambda(i-1,pop)))

        fact1 = frac_a_m1 * (1.-frac_lambda)
        fact2 = frac_a_m1 * frac_lambda
        fact3 = frac_a * (1.-frac_lambda)
        fact4 = frac_a * frac_lambda

     else ! Decreasing order
        do i=2,op_file_n_lambda(pop)-1
           if (log(op_file_lambda(i,pop)) < log_wavel) exit
        enddo
        !log(op_file_lambda(i-1)) > log_wavel >  log(op_file_lambda(i))

        frac_lambda = (log_wavel-log(op_file_lambda(i,pop)))/(log(op_file_lambda(i-1,pop))-log(op_file_lambda(i,pop)))

        fact2 = (1.-frac_a)  * (1.-frac_lambda)
        fact1 = (1.-frac_a) * frac_lambda
        fact4 = frac_a * (1.-frac_lambda)
        fact3 = frac_a * frac_lambda
     endif

     qext = exp(log(op_file_Qext(i-1,j-1,pop)) * fact1 &
          + log(op_file_Qext(i,j-1,pop)) * fact2 &
          + log(op_file_Qext(i-1,j,pop)) *  fact3 &
          + log(op_file_Qext(i,j,pop)) * fact4)

     qsca = exp(log(max(op_file_Qsca(i-1,j-1,pop),tiny_real)) * fact1 &
          + log(max(op_file_Qsca(i,j-1,pop),tiny_real)) * fact2 &
          + log(max(op_file_Qsca(i-1,j,pop),tiny_real)) * fact3 &
          + log(max(op_file_Qsca(i,j,pop),tiny_real)) * fact4)

     gsca = op_file_g(i-1,j-1,pop) * fact1 &
          + op_file_g(i,j-1,pop) * fact2 &
          + op_file_g(i-1,j,pop) * fact3 &
          + op_file_g(i,j,pop) * fact4
  endif

  !if ((igrain == dust_pop(pop)%ind_fin).and.(lambda==n_lambda)) call free_mem_opacity_file()

  !! Mueller matrices
  if (aniso_method==1) then ! we have to compute a phase function from g

     ! HG with the g interpolated in the table
     do j=0,nang_scatt
        s11(j)=((1-gsca**2)/(2.0))*(1+gsca**2-2*gsca*cos((real(j))/real(nang_scatt)*pi))**(-1.5)
     enddo

     ! Zero polarizability
     s12=0.0 ; s22 = 0.0 ; s33 = 0.0 ; s34 = 0.0 ; s44 = 0.0

     if (scattering_method==1) then
        prob_s11(lambda,igrain,0)=0.0
        dtheta = pi/real(nang_scatt)
        do j=2,nang_scatt ! scattering probability up to angle j, j=0 is skipped because sin(theta) = 0
           theta = real(j)*dtheta
           prob_s11(lambda,igrain,j)=prob_s11(lambda,igrain,j-1)+s11(j)*sin(theta)*dtheta
        enddo

        ! s11 is calculated such that the normalization is: 1.0 (HG definition)
        ! there is a numerical issue when x >> 1 because the angular resolution is not sufficient
        ! The diffraction peak is missed (particularly between 0 and 1)
        somme_prob = 1.0
        prob_s11(lambda,igrain,1:nang_scatt) = prob_s11(lambda,igrain,1:nang_scatt) + &
             somme_prob - prob_s11(lambda,igrain,nang_scatt)

        ! Normalization of the cumulative probability to 1
        prob_s11(lambda,igrain,:)=prob_s11(lambda,igrain,:)/somme_prob
     endif ! scattering_method==1

     do j=0,nang_scatt
        if (scattering_method==1) then ! Mueller matrix per grain
           ! Normalization for scattering according to phase function (tab_s11=1.0 used in Stokes)
           norm=s11(j)
        else ! Otherwise normalization to Qsca
           ! Default normalization is 1 for HG
           if (qsca > 1e-35) then
              norm = 1./qsca
           else ! we don't care there won't be any scattering
              norm = huge_real* 1e-4
           endif
        endif

        tab_s11(j,igrain,lambda) = s11(j) / norm
        tab_s12(j,igrain,lambda) = s12(j) / norm
        tab_s22(j,igrain,lambda) = s22(j) / norm
        tab_s33(j,igrain,lambda) = s33(j) / norm
        tab_s34(j,igrain,lambda) = s34(j) / norm
        tab_s44(j,igrain,lambda) = s44(j) / norm
     enddo

  endif ! aniso_method ==1

  return

end subroutine mueller_opacity_file

!**********************************************************************

subroutine update_Stokes(S,u0,v0,w0,u1,v1,w1, M)
  ! Astronomical convention used: position angle
  ! Calculated counter-clockwise from celestial north
  !
  ! Francois Menard, Montreal, 15 February 1989
  !
  ! C. Pinte
  !  - Modification 22/12/03 (C. Pinte): index l for scattering grain size
  !  - Energy normalization: the photon leaves with the energy it entered with
  !
  ! Scatterer coordinates
  !  - x is towards the observer
  !  - y and z form a right-handed basis (in the correct direction)
  !   - y to the right and z upwards
  !
  ! Calculation of the omega angle between the scattering plane
  ! and the projected celestial north (equatorial coordinates)
  !
  ! In the nebula's coordinate system
  !    vector 1 (from point "0" to point "1") = (u0,v0,w0)
  !    vector 2 (from point "1" to point "2") = (u1,v1,w1)
  !
  !  Transform so that v2prime = (1,0,0)
  !  The observer is then at +x in the new system

  implicit none

  real(kind=dp), intent(in) ::  u0,v0,w0,u1,v1,w1
  real(kind=dp), dimension(4,4), intent(in) :: M
  real(kind=dp), dimension(4), intent(inout) :: S

  real :: sinw, cosw, omega, theta, costhet,  xnyp
  real(kind=dp) :: v1pi, v1pj, v1pk, S1_0
  real(kind=dp), dimension(4,4) :: ROP, RPO
  real(kind=dp), dimension(4) :: C, D


  ! Transformation for v1prime
  call rotation(u0,v0,w0,u1,v1,w1,v1pi,v1pj,v1pk)

  ! Calculation of angles for the rotation
  ! The normal yprime is the cross product of v1prime x v2prime
  !
  !  yprimei = 0.0
  !  yprimej = v1pk
  !  yprimek = -v1pj
  xnyp = sqrt(v1pk*v1pk + v1pj*v1pj)
  if (xnyp < 1e-10) then
     xnyp = 0.0
     costhet = 1.0
  else
     costhet = -1.0*v1pj / xnyp
  endif

  ! Calculation of the angle between the normal and the z axis (theta)
  theta = acos(costhet)
  if (theta >= pi) theta = 0.0

  ! the scattering plane is at +/- 90deg from the normal
  theta = theta + half_pi

  ! In the rotation matrices the angle is omega = 2 * theta
  omega = 2.0 * theta

  ! next if because arccos only goes from 0 a pi only
  ! The +/- to distinguish the rotation sense
  if (v1pk < 0.0) omega = -1.0 * omega

  ! Calculation of rotation matrix elements
  !
  ! RPO = rotation from the point to the original system
  ! ROP = rotation from the original system to the point
  ! (brings the z-axis into the scattering plane)
  cosw = cos(omega)
  sinw = sin(omega)

  if (abs(cosw) < 1e-06) cosw = 0.0
  if (abs(sinw) < 1e-06) sinw = 0.0

  RPO(:,:) = 0.0
  ROP(:,:) = 0.0

  RPO(1,1) = 1.0
  ROP(1,1) = 1.0
  RPO(2,2) = cosw
  ROP(2,2) = cosw


  RPO(2,3) = sinw
  ROP(3,2) = sinw
  RPO(3,2) = -1.0 * sinw
  ROP(2,3) = -1.0 * sinw

  RPO(3,3) = cosw
  ROP(3,3) = cosw
  RPO(4,4) = 1.0
  ROP(4,4) = 1.0

  S1_0 = S(1)
  ! Stoke final = RPO * M * ROP * Stoke initial
  C=matmul(ROP,S)
  D=matmul(M,C)
  S=matmul(RPO,D)

  ! Energy normalization: the photon leaves with the energy it entered with
  ! I outgoing = I incoming if scattered according to s11 (tab_s11=1.0 normalized in mueller2)
  ! I outgoing = I incoming * s11 if uniform scattering
  if (S(1) > tiny_real) S(:) = S(:) * M(1,1)*S1_0/S(1)

  return

end subroutine update_Stokes

!**********************************************************************

subroutine get_Mueller_matrix_per_grain(lambda,itheta,frac,igrain, M)

  integer, intent(in) :: lambda,itheta,igrain
  real, intent(in) :: frac
  real(kind=dp), dimension(4,4), intent(out) :: M

  real :: frac_m1

  frac_m1 = 1.0 - frac

  M(:,:) = 0.0_dp
  M(1,1) = tab_s11(itheta,igrain,lambda) * frac + tab_s11(itheta-1,igrain,lambda) * frac_m1
  M(2,2) = tab_s22(itheta,igrain,lambda) * frac + tab_s22(itheta-1,igrain,lambda) * frac_m1
  M(1,2) = tab_s12(itheta,igrain,lambda) * frac + tab_s12(itheta-1,igrain,lambda) * frac_m1
  M(2,1) = M(1,2)
  M(3,3) = tab_s33(itheta,igrain,lambda) * frac + tab_s33(itheta-1,igrain,lambda) * frac_m1
  M(4,4) = tab_s44(itheta,igrain,lambda) * frac + tab_s44(itheta-1,igrain,lambda) * frac_m1
  M(3,4) = -tab_s34(itheta,igrain,lambda)* frac - tab_s34(itheta-1,igrain,lambda) * frac_m1
  M(4,3) = -M(3,4)

  return

end subroutine get_Mueller_matrix_per_grain

!**********************************************************************

subroutine get_Mueller_matrix_per_cell(lambda,itheta,frac,icell, M)

  integer, intent(in) :: lambda,itheta,icell
  real, intent(in) :: frac
  real(kind=dp), dimension(4,4), intent(out) :: M

  real :: frac_m1

  frac_m1 = 1.0 - frac

  M(:,:) = 0.0_dp
  M(1,1) = 1.0 ! Mueller matrix is normalized to 1.0 as we select the scattering angle
  M(2,2) = tab_s22_o_s11_pos(itheta,icell,lambda) * frac +  tab_s22_o_s11_pos(itheta-1,icell,lambda) * frac_m1
  M(1,2) = tab_s12_o_s11_pos(itheta,icell,lambda) * frac +  tab_s12_o_s11_pos(itheta-1,icell,lambda) * frac_m1
  M(2,1) = M(1,2)
  M(3,3) = tab_s33_o_s11_pos(itheta,icell,lambda) * frac +  tab_s33_o_s11_pos(itheta-1,icell,lambda) * frac_m1
  M(4,4) = tab_s44_o_s11_pos(itheta,icell,lambda) * frac +  tab_s44_o_s11_pos(itheta-1,icell,lambda) * frac_m1
  M(3,4) = -tab_s34_o_s11_pos(itheta,icell,lambda)* frac -  tab_s34_o_s11_pos(itheta-1,icell,lambda) * frac_m1
  M(4,3) = -M(3,4)

  return

end subroutine get_Mueller_matrix_per_cell

!**********************************************************************

subroutine hg(g, rand, itheta, cospsi)
  ! Drawing of the scattering angle via the Henyey-Greenstein function
  ! itheta is the index i of the angle from 1 to 180 corresponding to i-0.5 degrees
  ! rand is a random number between 0 and 1
  !
  ! C. Pinte    23/10/2004

  implicit none

  real, intent(in) :: g, rand
  integer, intent(out) :: itheta
  real(kind=dp), intent(out) :: cospsi

  real(kind=dp) :: g1, g2, rand_dp

  rand_dp = min(real(rand,kind=dp), 1.0_dp-1e-6_dp)

  if (abs(g) > tiny_real) then
     g1 = g ! dp
     g2 = g1*g1
     cospsi = (1.0_dp + g2 - ((1.0_dp - g2) / (1.0_dp - g1 + 2.0_dp*g1*rand_dp))**2) / (2.0_dp * g1)
  else ! g=0 --> isotropic scattering
     cospsi=2.0_dp*rand_dp-1.0_dp
  endif

  itheta = floor(acos(cospsi)*180.0_dp/pi)+1
  if (itheta > nang_scatt) itheta = nang_scatt

  return
end subroutine hg

!***********************************************************

subroutine angle_diff_theta(lambda, igrain, rand, rand2, itheta, cospsi)
  ! Calculation of the cosine of the scattering angle
  ! from the integral of the pretabulated s11
  ! itheta is the index i of the angle from 1 to 180 corresponding to i-0.5 degree
  ! cospsi is drawn uniformly in the 1 degree bin around angle i-0.5 degrees
  ! itheta is used for pretabulated values
  ! cospsi is used for the flight direction
  ! C. Pinte 23/10/04

  implicit none

  integer, intent(in) :: lambda, igrain
  real, intent(in) :: rand, rand2
  integer, intent(out) :: itheta
  real(kind=dp), intent(out) :: cospsi

  integer :: k, kmin, kmax

  kmin=0
  kmax=nang_scatt
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_s11(lambda,igrain,k) < rand) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   k=kmax

   itheta=k
   !cospsi=cos((real(k)-0.5)*pi/180.)

   ! Random draw of the scattering angle between angle k and angle k-1
   ! uniform scattering (linear in cos)
   cospsi=cos((real(k)-1.0)*pi/real(nang_scatt)) + &
        rand2*(cos((real(k))*pi/real(nang_scatt))-cos((real(k)-1.0)*pi/real(nang_scatt)))

   return

end subroutine angle_diff_theta

!**********************************************************************

subroutine angle_diff_theta_pos(lambda, icell, rand, rand2, itheta, cospsi)
  ! Calculation of the cosine of the scattering angle
  ! from the integral of s11 pretabulated per cell
  ! itheta is the index i of the angle from 1 to 180 corresponding to i-0.5 degrees
  ! cospsi is drawn uniformly in the 1 degree bin around angle i-0.5 degrees
  ! itheta is used for pretabulated values
  ! cospsi is used for the flight direction
  ! C. Pinte 9/01/05

  implicit none

  integer, intent(in) :: lambda,icell
  real, intent(in) :: rand, rand2
  integer, intent(out) :: itheta
  real(kind=dp), intent(out) :: cospsi

  integer :: k, kmin, kmax

  kmin=0
  kmax=nang_scatt
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_s11_pos(k,icell,lambda) < rand) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   k=kmax

   itheta=k
   !cospsi=cos((real(k)-0.5)*pi/180.)

   ! Random draw of the scattering angle between angle k and angle k-1
   ! uniform scattering (linear in cos)
   cospsi=cos((real(k,kind=dp)-1.0_dp)*pi/real(nang_scatt,kind=dp)) + &
        rand2*(cos((real(k,kind=dp))*pi/real(nang_scatt,kind=dp))-cos((real(k,kind=dp)-1.0_dp)*pi/real(nang_scatt,kind=dp)))

   return

end subroutine angle_diff_theta_pos

!**********************************************************************

subroutine funcd(x,fval,fderiv,ppp,A)
  ! calculate the scattering angle distribution function phi
  ! and its derivative to determine the zero by the Newton method
  ! C. Pinte    23/10/2004

  implicit none

  real, intent(in) :: x,ppp,A
  real, intent(out) :: fval,fderiv

  ! Are you sure about the sign here?
  fval=2*x-ppp*sin(2*x)-A
  fderiv=2-ppp*2*cos(2*x)

  return

end subroutine funcd

!**********************************************************************

subroutine angle_diff_phi(lambda,igrain, I, Q, U, itheta, frac, rand, phi)
  ! Draw of the photon's scattering angle phi
  ! Uniform for an unpolarized wave
  ! Following the Mie phase function for a polarized photon
  ! C. Pinte    23/10/2004
  ! Added case where Mueller matrices are given as input
  ! 20/04/2023

  implicit none

  integer, intent(in) :: lambda,igrain, itheta
  real, intent(in) :: frac, I, Q, U, rand
  real(kind=dp), intent(out) :: phi

  real :: p, pp, ppp, phi1, phi2, frac_m1

  real(kind=dp) :: Q_dp, U_dp, Ip

  frac_m1 = 1.0 - frac

!  write(*,*) 'in',l, itheta, I, Q, U, rand

  ! Polarized flux and polarization rate
  Q_dp=Q;U_dp=U
  Ip=sqrt(Q_dp*Q_dp+U_dp*U_dp)
  p=Ip/I

  ! polarizability
  pp= (tab_s12(itheta,igrain,lambda) * frac + tab_s12(itheta-1,igrain,lambda) * frac_m1) &
       / (tab_s11(itheta,igrain,lambda) * frac + tab_s11(itheta-1,igrain,lambda) * frac_m1)

  ppp=p*pp
!  write(*,*) p,pp,ppp

  if (abs(ppp) > 1.e-3) then
     ! Measure the angle of the polarization plane relative to celestial North
     phi1=0.5*acos(Q_dp/Ip) ! This is where we need the dp
     if (U < 0.0) then
        phi1= -phi1
     endif
     ! draw of angle between the scattering plane and the polarization plane
     ! according to A = 2*phi - ppp*sin(2*phi)   A=4*pi*rand
     phi2=rtsafe(funcd,0.,2*real(pi),1.e-4,ppp,4*real(pi)*rand)
!     write(*,*) 'test', phi1 , Q/Ip
     phi=phi1+phi2
     if (phi > pi) phi = phi -2*pi
     if (phi < -pi) phi = phi +2*pi
  else ! Uniform selection
     phi= pi * (2._dp*rand -1.0_dp)
  endif
!  write(*,*) phi/pi

  return

end subroutine angle_diff_phi

!**********************************************************************

real function rtsafe(funcd,x1,x2,xacc,ppp,A)
  ! Find the zero of a function
  ! Here, the scattering angle in phi for a polarized photon
  ! following the distribution function given by funcd
  ! C. Pinte    23/10/2004

  implicit none
  real, intent(in) :: x1,x2,xacc,ppp,A

  integer, parameter :: i4b = selected_int_kind(9)
  integer(i4b), parameter :: maxit=100
  integer(i4b) :: j
  real :: df,dx,dxold,f,fh,fl,temp,xh,xl

  interface
     subroutine funcd(x,fval,fderiv,ppp,A)
       implicit none
       real, intent(in) :: x,ppp,A
       real, intent(out) :: fval,fderiv
     end subroutine funcd
  end interface

  call funcd(x1,fl,df,ppp,A)
  call funcd(x2,fh,df,ppp,A)
  if ((fl > 0.0 .and. fh > 0.0) .or. (fl < 0.0 .and. fh < 0.0)) call error('root must be bracketed in rtsafe')

  if (fl == 0.0) then
     rtsafe=x1
     return
  else if (fh == 0.0) then
     rtsafe=x2
     return
  else if (fl < 0.0) then
     xl=x1
     xh=x2
  else
     xh=x1
     xl=x2
  end if
  rtsafe=0.5*(x1+x2)
  dxold=abs(x2-x1)
  dx=dxold
  call funcd(rtsafe,f,df,ppp,A)
  do j=1,MAXIT
     if (((rtsafe-xh)*df-f)*((rtsafe-xl)*df-f) >= 0.0 .or. &
          abs(2.0*f) > abs(dxold*df) ) then
        dxold=dx
        dx=0.5*(xh-xl)
        rtsafe=xl+dx
        if (xl == rtsafe) return
     else
        dxold=dx
        dx=f/df
        temp=rtsafe
        rtsafe=rtsafe-dx
        if (temp == rtsafe) return
     end if
     if (abs(dx) < xacc) return
     call funcd(rtsafe,f,df,ppp,A)
     if (f < 0.0) then
        xl=rtsafe
     else
        xh=rtsafe
     end if
  end do

  call error('rtsafe: exceeded maximum iterations')

end function rtsafe

!**********************************************************************

subroutine radius_aggregate()

  implicit none

  integer :: i, alloc_status
  real :: wavelength
  real, dimension(:), allocatable :: x, y, z, r, eps1, eps2


  open(unit=1,file=trim(aggregate_file), status='old')
  read(1,*) wavelength
  if (abs(wavelength - tab_lambda(1)) < 1.0e-5) call error("wavelength does not correspond to wavelength of the Mueller matrix")

  read(1,*) n_grains_tot
  allocate(x(n_grains_tot), y(n_grains_tot), z(n_grains_tot), r(n_grains_tot), eps1(n_grains_tot), &
       eps2(n_grains_tot),stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error in radius_aggregate')

  do i=1, n_grains_tot
     read(1,*) x(i), y(i), z(i), r(i), eps1(i), eps2(i)
  enddo
  close(unit=1)

  R_sph_same_M = (sum(r(:)**3))**(1.0/3.0)
  return

end subroutine radius_aggregate

!***********************************************************

end module scattering
