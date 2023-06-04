module scattering

  use parametres
  use constantes
  use wavelengths
  use grains
  use utils
  use read_opacity

  implicit none

  procedure(get_Mueller_matrix_per_grain), pointer :: get_Mueller_matrix => null()

  contains

subroutine setup_scattering()

  ! parametrage methode de diffusion en fonction de la taille de la grille
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

  ! Inutile depuis le passage en allocation dynamique
  !  integer, parameter :: NMXX=20000000
  ! defaut : NMXX=20000
  ! On peut passer a  NMXX=2000000 : ca marche jusque a=10cm en bande B
  ! mais faut pas etre presse
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
! 04/03/04 (CP): passage en fortran 90
! 13/10/04 (CP): passage a des angles demi-entier (milieu du bin) pour pouvoir faire l'integration de s11(theta)
! l'angle 0° est calcule a part car on en a besoin pour Qext
! de meme, on aurait besoin de 180° pour Qback mais ca ne sert pas donc je ne calcule pas
! 16/12/04 (CP): Remplacement imag par aimag (Pas forcement ok avec tous les
! compilateurs) mais standard f95 et 2003
! 24/03/05 (CP): allocation dynamique de D. Permet de placer le tableau dans la zone data
! et de l'allouer avec juste le bon nombre de terme.
! 07/02/08 (CP): repassage en angles entiers avec 0 et 180 explicitement calcules
! necessaire pour integration du ray_tracing : passage a 2*nang + 1 angles
!***********************************************************************


!*** Safety checks
!  if(NANG > MXNANG) then
!     write(*,*)'***Error: NANG > MXNANG in bhmie'
!     stop
!  endif
!  if(NANG < 2) then
!     write(*,*)'***Error: NANG doit etre >= 2'
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

subroutine mueller_Mie(lambda,taille_grain,x,amu1,amu2, qext,qsca,gsca)
  !***************************************************************
  ! Calcule les elements de la matrice de diffusion a partir de
  ! la sous-routine bhmie (grains spheriques)
  !
  ! C. Pinte Fevrier 2004
  !****************************************************************

  implicit none
  integer, intent(in) :: lambda, taille_grain
  real, intent(in) :: amu1, amu2
  real, intent(in) :: x ! 2*pi*a/wl
  real, intent(out) :: qext, qsca, gsca

  integer :: j, jj, nang

  complex, dimension(nang_scatt+1) :: S1,S2

  real :: vi1, vi2, qback, factor
  complex :: s, refrel
  real, dimension(0:nang_scatt) ::  s11,s12,s33,s34

  refrel = cmplx(amu1,amu2)

  if (modulo(nang_scatt,2) == 1) call error("nang_scatt must be an EVEN number")

  ! Si fonction de HG, on ne calcule pas la fonction de phase
  if (aniso_method==2) then
     nang=1
  else
     nang= (nang_scatt+1) / 2 + 1
  endif

  call bhmie(x,refrel,nang, s1,s2,qext,qsca,qback,gsca)

  if (lforce_HG) gsca = forced_g
  if (lisotropic) gsca = 0.0

  ! The normalisation of bhmie for s11 is 0.5*x**2*Qsca
  ! --> We correct by 0.5*x**2 to have a normalisation to Qsca Qsa
  factor = 1 / (0.5 * x**2)

  ! Passage des valeurs dans les tableaux de mcfost
  if (aniso_method==1) then
     ! Calcul des elements de la matrice de diffusion
     ! indices decales de 1 par rapport a bhmie
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

     call normalise_Mueller_matrix(lambda,taille_grain, s11,s12,s11,s33,s34,s33, qsca)
  endif

  return

end subroutine mueller_Mie

!***************************************************

subroutine Mueller_input(lambda, k_abs,k_sca,gsca)

  integer, intent(in) :: lambda
  real, intent(out) :: k_abs, k_sca, gsca

  character(len=128) :: string, filename

  integer :: EoF, alloc_status,iformat,nlam,nmu,ilam,iline,iang,nang, l
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
  do ilam=1,nlam
     read(1,*) wavel(l), kabs(l), ksca(l), g(l)
  enddo

  wl = tab_lambda(wl)
  if (wl < minval(wavel)) call error("Mueller matrices: wavelength is smaller than tabulated")
  if (wl > maxval(wavel)) call error("Mueller matrices: wavelength is larger than tabulated")

  if (wavel(2) >  wavel(1)) then ! increasing order
     do ilam=1,nlam
        if  (wavel(l) > wl) then
           frac = (wavel(l) - wl) / (wavel(l) - wavel(l-1))
           exit
        endif
     enddo
  else ! decreasing order
     do ilam=1,nlam
        if  (wavel(l) < wl) then
           frac = (wl - wavel(l)) / (wavel(l-1) - wavel(l))
           exit
        endif
     enddo
  endif
  if ((frac > 1) .or. (frac < 0)) call error("Mueller matrices: interpolation error")

  if (lscat) then
     allocate(angles(nang), f11(nlam,nang),f12(nlam,nang),f22(nlam,nang),&
          f33(nlam,nang),f34(nlam,nang),f44(nlam,nang), stat = alloc_status)

     ! Read the angular grid
     do iang=1,nang
        read(1,*) angles(iang)
     enddo

     ! Read the scattering matrix.
     do ilam=1,nlam
        do iang=1,nang
           read(1,*) f11(ilam,iang),f12(ilam,iang), f22(ilam,iang), f33(ilam,iang), f34(ilam,iang), f44(ilam,iang)
        enddo
     enddo
  endif
  close(1)

  ! Interpolate at correct wavelength
  ! todo : find ilam

  k_abs = kabs(ilam-1) * frac + kabs(ilam) * frac_m1
  k_sca = ksca(ilam-1) * frac + ksca(ilam) * frac_m1
  gsca = ksca(ilam-1) * frac + ksca(ilam) * frac_m1

  s11(:) = f11(ilam-1,:) * frac + f11(ilam,:) * frac_m1
  s12(:) = f11(ilam-1,:) * frac + f12(ilam,:) * frac_m1
  s22(:) = f11(ilam-1,:) * frac + f22(ilam,:) * frac_m1
  s33(:) = f11(ilam-1,:) * frac + f33(ilam,:) * frac_m1
  s34(:) = f11(ilam-1,:) * frac + f34(ilam,:) * frac_m1
  s44(:) = f11(ilam-1,:) * frac + f44(ilam,:) * frac_m1

  ! Deallocate arrays
  deallocate(wavel,kabs,ksca,g)
  if (lscat) then
     deallocate(angles)
     deallocate(f11,f12,f22,f33,f34,f44)
  endif

  return

end subroutine Mueller_input

!***************************************************

subroutine normalise_Mueller_matrix(lambda,taille_grain, s11,s12,s22,s33,s34,s44, qsca)

  integer, intent(in) :: lambda,taille_grain
  real, dimension(0:nang_scatt), intent(in) ::  s11,s12,s22,s33,s34,s44
  real, intent(in) :: qsca

  integer :: j
  real :: norm, somme_prob, theta, dtheta

  ! Integration S11 pour tirer angle
  if (scattering_method==1) then
     prob_s11(lambda,taille_grain,0)=0.0
     dtheta = pi/real(nang_scatt)
     do j=2,nang_scatt ! probabilite de diffusion jusqu'a l'angle j, on saute j=0 car sin(theta) = 0
        theta = real(j)*dtheta
        prob_s11(lambda,taille_grain,j)=prob_s11(lambda,taille_grain,j-1)+s11(j)*sin(theta)*dtheta
     enddo

     ! il y a un soucis numerique quand x >> 1 car la resolution en angle n'est pas suffisante
     ! On rate le pic de diffraction (en particulier entre 0 et 1)
     somme_prob = qsca
     prob_s11(lambda,taille_grain,1:nang_scatt) = prob_s11(lambda,taille_grain,1:nang_scatt) + &
          somme_prob - prob_s11(lambda,taille_grain,nang_scatt)

     ! Normalisation de la proba cumulee a 1
     prob_s11(lambda,taille_grain,:)=prob_s11(lambda,taille_grain,:)/somme_prob
  endif ! scattering_method==1

  do j=0,nang_scatt
     if (scattering_method==1) then ! Matrice de Mueller par grain
        ! Normalisation pour diffusion selon fonction de phase (tab_s11=1.0 sert dans stokes)
        norm = 1./s11(j)
     else ! Sinon normalisation a Qsca
        norm = 1.
     endif

     tab_s11(j,taille_grain,lambda) = s11(j) * norm
     tab_s12(j,taille_grain,lambda) = s12(j) * norm
     tab_s22(j,taille_grain,lambda) = s22(j) * norm
     tab_s33(j,taille_grain,lambda) = s33(j) * norm
     tab_s34(j,taille_grain,lambda) = s34(j) * norm
     tab_s44(j,taille_grain,lambda) = s44(j) * norm
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

subroutine mueller_GMM(lambda,taille_grain, qext,qsca,gsca)
!***************************************************************
! calcule les elements de la matrice de diffusion a partir du
! code gmm01TrA (clusters de sphères)
!     Aggrégats
!
!        CALCULE AUSSI "G" = LE PARAMETRE D'ASYMETRIE
!
! C. Pinte
! 04/07/2005
!****************************************************************

  implicit none

  integer, intent(in) :: lambda, taille_grain
  real, intent(out) :: qext, qsca

  integer, parameter :: nang2 = 2*nang_scatt+1

  integer :: j, i

  real :: gsca, norme, somme_sin, somme_prob, somme2
  real :: cext,cabs,csca,cbak,cpr,assym
  real :: cextv,cabsv,cscav,cbakv,cprv
  real :: cexts,cabss,cscas,cbaks,cprs

  real, dimension(nang2) :: dang
  real, dimension(4,4,nang2) :: mue

  real, dimension(4,4,2*nang_scatt) :: mueller

  integer :: idscmt=1

  character(len=128) :: string


  ! Il faut coriger pour la nouvelle definition de nang_scatt
  ! Il faut aussi corriger la normalisation de s11 (faire comme mueller_Mie)
  ! int S11 sin(theta) dtheta = Qsca, utilise la normalisation exacte, pas numerique
  call error("mueller_gmm needs to be updated")

  if(n_grains_tot > 1) call error("You must choose n_grains_tot=1")
  if (scattering_method /= 1) call error("You must choose scattering_method 1")

  ! Lecture du fichier de résultats de gmm01TrA : 'gmm01TrA.out'
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
  ! Fin lecture
  close(unit=1)


!  QABS=QEXT-QSCA
! Calcul des elements de la matrice de diffusion
! Calcul angle central du bin par interpolation linéaire
  do J=1,2*NANG_scatt
     mueller(:,:,j) = 0.5*(mue(:,:,j)+mue(:,:,j+1))
  enddo !j

  ! Integration S11 pour tirer angle
  prob_s11(lambda,taille_grain,0)=0.0
  somme_sin= 0.0
  somme2 = 0.0

  do j=1,2*nang_scatt
     prob_s11(lambda,taille_grain,j)=prob_s11(lambda,taille_grain,j-1)+&
          mueller(1,1,j)*sin((real(j)-0.5)/180.*pi)*pi/(2*nang_scatt)
     somme_sin = somme_sin + sin((real(j)-0.5)/180.*pi)*pi/(2*nang_scatt)
!     somme2=somme2+s12(j)*sin((real(j)-0.5)/180.*pi)*pi/(2*nang)
! Somme2 sert juste pour faire des plots
  enddo

  ! Normalisation
  somme_prob=prob_s11(lambda,taille_grain,2*nang_scatt) ! = (0.5*x**2*qsca)
  ! Soit int_0^\pi (i1(t)+i2(t)) sin(t) = x**2*qsca
  do j=1,2*nang_scatt
     prob_s11(lambda,taille_grain,j)=prob_s11(lambda,taille_grain,j)/somme_prob
  enddo


  do J=1,2*NANG_scatt
!     ! Normalisation pour diffusion isotrope et E_sca(theta)
!     if (j == 1)  then
!        norme = somme_prob/somme_sin
!     endif

! NORMALISATION ENLEVEE POUR LES CALCULS DES TAB_POS (MATRICES DE MUELLER
! PAR CELLULE)
! A REMETTRE POUR MATRICES DE MUELLER PAR GRAINS

     if (scattering_method==1) then
        ! Normalisation pour diffusion selon fonction de phase (tab_s11=1.0 sert dans stokes)
        norme=mueller(1,1,j) !* qext/q sca
        mueller(:,:,j) = mueller(:,:,j) / norme
     endif

     tab_s11(j,taille_grain,lambda) = mueller(1,1,j)
     tab_s12(j,taille_grain,lambda) = mueller(1,2,j)
     tab_s22(j,taille_grain,lambda) = mueller(2,2,j)
     tab_s33(j,taille_grain,lambda) = mueller(3,3,j)
     tab_s34(j,taille_grain,lambda) = mueller(3,4,j)
     tab_s44(j,taille_grain,lambda) = mueller(4,4,j)
  enddo

  gsca = assym
  qext=cext/(pi*R_sph_same_M)
  qsca=csca/(pi*R_sph_same_M)

  return

end subroutine mueller_GMM

!***************************************************

subroutine Fresnel_input(lambda,taille_grain, qext,qsca,gsca)

!***************************************************************
! Routine dérivé de Mueller_GMM.
! Calcule les elements de la matrice de diffusion a partir d'un
! fichier ascii contenant divers informations :
! La matrice de Mueller à chaque angle de diffraction de 0 à 180°,
! Qsca, Qext et le facteur d'assymétrie
! Voici un exemple de fichier :
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

  integer, intent(in) :: lambda, taille_grain
  real, intent(out) :: qext, qsca, gsca

  integer :: j, i

  real :: norme, somme_prob, theta, dtheta
  real :: ext,sca,assym

  real, dimension(0:nang_scatt) :: dang

  real, dimension(4,4,0:nang_scatt) :: mueller

  character(len=128) :: string

  if (modulo(nang_scatt,2) == 1) call error("nang_scatt must be an EVEN number")

  if (aniso_method == 2) call error ("You shoudn't use a hg function option when putting Mueller matrix in input")

  ! Lecture du fichier d'entrée :
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
  ! Fin lecture
  close(unit=1)

  if (scattering_method == 1) then
     ! Integration de S11 pour tirer angle
     ! Condition initiale : sin(0)=0
     prob_s11(lambda,taille_grain,0)=0.0
     dtheta = pi/real(nang_scatt)

     do j=1,nang_scatt
        theta = real(j)*dtheta
        prob_s11(lambda, taille_grain,j)=prob_s11(lambda, taille_grain,j-1)+&
             mueller(1,1,j)*sin(theta)*dtheta
     enddo

	  ! Normalisation
     somme_prob = prob_s11(lambda, taille_grain,nang_scatt)
     prob_s11(lambda, taille_grain,:)=prob_s11(lambda, taille_grain,:)/somme_prob
  else  !on va avoir besoin des valeurs non normalisés
     do i=1, 4
        do j=1, 4
           if (i*j>1) mueller(i,j,:) = mueller(i,j,:)*mueller(1,1,:)
        enddo
     enddo
     !on calcul maintenant l'intégrale de s11
     somme_prob=0.0
     dtheta = pi/real(nang_scatt)
     do j=1,nang_scatt
        theta = real(j)*dtheta
        somme_prob = somme_prob + mueller(1,1,j)*sin(theta)*dtheta
     enddo
  endif !scatt_meth

  do J=0,NANG_scatt
     if (scattering_method==1) then ! Matrice de Mueller par grain
        ! Normalisation pour diffusion selon fonction de phase (tab_s11=1.0 sert dans stokes)
        ! on ne normalise que S11 car les données sont déjà normalisées
        norme=mueller(1,1,j)
        if (norme > tiny_real) then
           mueller(1,1,j) = mueller(1,1,j) / norme
        else
           write (*,*) "at angle", real(j)*pi/real(nang_scatt)
           call error ("s11=0.0")
        endif
     else ! Sinon normalisation a Qsca. On utilise somme_prob.
        norme = somme_prob/sca
        if (norme > tiny_real) then
           mueller(:,:,j) = mueller(:,:,j)/norme
        else
           call error ("Error : s11 is normalised to 0")
        endif
     endif

     tab_s11(j,taille_grain,lambda) = mueller(1,1,j)
     tab_s12(j,taille_grain,lambda) = mueller(1,2,j)
     tab_s22(j,taille_grain,lambda) = mueller(2,2,j)
     tab_s33(j,taille_grain,lambda) = mueller(3,3,j)
     tab_s34(j,taille_grain,lambda) = mueller(3,4,j)
     tab_s44(j,taille_grain,lambda) = mueller(4,4,j)
  enddo
  qext = ext
  qsca = sca
  gsca = assym
  if (lforce_HG)  gsca = forced_g
  if (lisotropic)  gsca = 0.0

  return

end subroutine Fresnel_input

!***************************************************

 subroutine Fresnel_input_size(lambda,taille_grain, qext,qsca,gsca)

!***************************************************************
! Routine dérivé de Mueller_GMM.
! Calcule les elements de la matrice de diffusion a partir d'un
! fichier ascii par taille de grain contenant divers informations :
! La matrice de Mueller à chaque angle de diffraction de 0 à 180°
! Qsca, Qext, le facteur d'assymétrie
! Voici un exemple de fichier :
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
! L'adresse de ces fichiers est contenu dans un fichier texte
! trié dans l'ordre croissant suivant la taille de grain.
!
! F. Malaval
! 20/04/2023
!****************************************************************

  implicit none

  integer, intent(in) :: lambda, taille_grain
  real, intent(out) :: qext, qsca, gsca

  integer :: j, i

  real :: norme, somme_prob, theta, dtheta
  real :: ext,sca,assym, size

  real, dimension(0:nang_scatt) :: dang

  real, dimension(4,4,0:nang_scatt) :: mueller

  character(len=128) :: string, path

  if (modulo(nang_scatt,2) == 1) call error("nang_scatt must be an EVEN number")

  if (aniso_method == 2) call error ("You shoudn't use a hg function option when putting Mueller matrix in input")

  open(unit=13, file = mueller_file, status = 'old')
  !Dans ce cas, mueller_file contient les chemins vers les fichiers de chaque matrice,
  !triés dans l'ordre croissant suivant la taille de grain.
  do j=1, taille_grain
     read(13, *) size, path
  enddo
  close(unit=13)
  if (abs(size-r_grain(taille_grain))>r_grain(taille_grain)*0.00001) then
     write(*,*) "error, grain size in file is", size, "while expecting", r_grain(taille_grain)
     call error("Grain sizes do not match")
  endif

  ! Lecture du fichier d'entrée :
  open(unit=12,file=path,status='old')
  read(12,*) string
  read(12,*) ext,sca,assym
  read(12,*)
  read(12,*)  !'Matrice de mueller (4X4 coefficients pour chaque angle):'
  read(12,*) string
  do i=0,nang_scatt
     read(12,*) dang(i),mueller(1,1,i),mueller(1,2,i),mueller(1,3,i),mueller(1,4,i)
     read(12,*)   mueller(2,1,i),mueller(2,2,i),mueller(2,3,i),mueller(2,4,i)
     read(12,*)   mueller(3,1,i),mueller(3,2,i),mueller(3,3,i),mueller(3,4,i)
     read(12,*)   mueller(4,1,i),mueller(4,2,i),mueller(4,3,i),mueller(4,4,i)
  enddo
  close(12)
  ! Fin lecture

  close(unit=1)

  if (scattering_method == 1) then
     ! Integration S11 pour tirer angle
     ! Condition initiale : sin(0)=0
     prob_s11(lambda,taille_grain,0)=0.0
     dtheta = pi/real(nang_scatt)

     do j=1,nang_scatt
        theta = real(j)*dtheta
        prob_s11(lambda, taille_grain,j)=prob_s11(lambda, taille_grain,j-1)+&
             mueller(1,1,j)*sin(theta)*dtheta
     enddo

     ! Normalisation
     somme_prob= prob_s11(lambda, taille_grain,nang_scatt)
     prob_s11(lambda, taille_grain,:)=prob_s11(lambda, taille_grain,:)/somme_prob
  else  !on va avoir besoin des valeurs non normalisés
     do i=1, 4
        do j=1, 4
           if (i*j>1) mueller(i,j,:) = mueller(i,j,:)*mueller(1,1,:)
        enddo
     enddo
     !on calcul maintenant l'intégrale de s11
     somme_prob=0.0
     dtheta = pi/real(nang_scatt)
     do j=1,nang_scatt
        theta = real(j)*dtheta
        somme_prob = somme_prob + mueller(1,1,j)*sin(theta)*dtheta
     enddo
  endif !scatt_meth


  do J=0,NANG_scatt
     if (scattering_method==1) then ! Matrice de Mueller par grain
        ! Normalisation pour diffusion selon fonction de phase (tab_s11=1.0 sert dans stokes)
        ! on ne normalise que S11 car les données sont déjà normalisées
        norme=mueller(1,1,j)
        if (norme > tiny_real) then
           mueller(1,1,j) = mueller(1,1,j) / norme
        else
           write (*,*) "at angle", real(j)*pi/real(nang_scatt)
           call error ("s11=0.0")
        endif
     else ! Sinon normalisation a Qsca. On utilise somme_prob.
        norme = somme_prob/sca
        if (norme > tiny_real) then
           mueller(:,:,j) = mueller(:,:,j)/norme
        else
           call error ("Error : s11 is normalised to 0")
        endif
     endif

     tab_s11(j,taille_grain,lambda) = mueller(1,1,j)
     tab_s12(j,taille_grain,lambda) = mueller(1,2,j)
     tab_s22(j,taille_grain,lambda) = mueller(2,2,j)
     tab_s33(j,taille_grain,lambda) = mueller(3,3,j)
     tab_s34(j,taille_grain,lambda) = mueller(3,4,j)
     tab_s44(j,taille_grain,lambda) = mueller(4,4,j)
  enddo
  qext = ext
  qsca = sca
  gsca = assym
  if (lforce_HG)  gsca = forced_g
  if (lisotropic)  gsca = 0.0

  return

end subroutine Fresnel_input_size

!***************************************************

subroutine mueller_opacity_file(lambda,taille_grain, qext,qsca,gsca)
  ! interpolation bi-lineaire (en log-log) des sections efficaces
  ! pour grains apres lecture du fichier d'opacite
  ! En particulier pour les PAHs de Draine
  ! Suppose une HG pour la fonction de phase et une polarisabilite nulle !!
  ! C. Pinte
  ! 31/01/07

  implicit none

  integer, intent(in) :: taille_grain, lambda
  real, intent(out) :: qext,qsca,gsca

  real :: frac_a, frac_a_m1, frac_lambda, fact1, fact2, fact3, fact4
  integer :: i, j, pop, N

  real, dimension(0:nang_scatt) ::  S11,S12,S33,S34
  real :: norme, somme_prob, log_a, log_wavel, wl_min, wl_max, theta, dtheta

  log_a=log(r_grain(taille_grain))
  log_wavel = log(tab_lambda(lambda))

  pop=grain(taille_grain)%pop

  if ((taille_grain == dust_pop(pop)%ind_debut).and.(lambda==1)) call read_opacity_file(pop)

  ! Ordre croissant pour les tailles de grains
  if (r_grain(taille_grain) < exp(op_file_log_r_grain(1,pop))) then
     if (lambda==1) then
        write(*,*) "WARNING: index=",taille_grain, "grain size=",r_grain(taille_grain)
        write(*,*) "Minimum grain size in opacity file is",  exp(op_file_log_r_grain(1,pop))
        write(*,*) "Smaller grains are assumed to have the same opacity"
     endif
    j = 2
    frac_a = 0.0 ; frac_a_m1 = 1.0
  else if (r_grain(taille_grain) > exp(op_file_log_r_grain(op_file_na(pop),pop))) then
     if (lambda==1) then
        write(*,*) "WARNING: index=",taille_grain, "grain size=",r_grain(taille_grain)
        write(*,*) "Maximum grain size in opacity file is",  exp(op_file_log_r_grain(op_file_na,pop))
        write(*,*) "Larger grains are assumed to have the same opacity"
     endif
     j = op_file_na(pop)
     frac_a = 1.0 ; frac_a_m1 = 0.
  else
     ! Recherche en taille de grain
     ! tableau croissant
     do j=2,op_file_na(pop)
        if (op_file_log_r_grain(j,pop) > log_a) exit
     enddo
     frac_a = (log_a-op_file_log_r_grain(j-1,pop))/(op_file_log_r_grain(j,pop)-op_file_log_r_grain(j-1,pop))
     frac_a_m1 = 1 - frac_a
  endif

  ! Moyennage en longueur d'onde
  wl_min = tab_lambda_inf(lambda)
  wl_max = tab_lambda_sup(lambda)

  qext = 0.0 ; qsca = 0.0 ; gsca = 0.0 ; norme = 0 ; N=0
  do i=1,op_file_n_lambda(pop)
     if ((op_file_lambda(i,pop) > wl_min).and.(op_file_lambda(i,pop) < wl_max)) then
        N = N+1
        norme = norme + op_file_delta_lambda(i,pop)
        qext = qext + (frac_a * op_file_Qext(i,j,pop) +  frac_a_m1 * op_file_Qext(i,j-1,pop) ) * op_file_delta_lambda(i,pop)
        qsca = qsca + (frac_a * op_file_Qsca(i,j,pop) +  frac_a_m1 * op_file_Qsca(i,j-1,pop) ) * op_file_delta_lambda(i,pop)
        gsca = gsca + (frac_a * op_file_g(i,j,pop) +  frac_a_m1 * op_file_g(i,j-1,pop) ) * op_file_delta_lambda(i,pop)
     endif
  enddo !i
  !write(*,*) lambda, N, norme

  if (norme > 0) then ! on fait la moyenne sur les points selectionnes
     qext = qext / norme
     qsca = qsca / norme
     gsca = gsca / norme
  else ! on peut pas moyenner, on fait une interpolation en log
     ! Recherche en longueur d'onde
     if (op_file_lambda(2,pop) > op_file_lambda(1,pop)) then  ! Ordre croisant
        do i=2,op_file_n_lambda(pop)
           if (log(op_file_lambda(i,pop)) > log_wavel) exit
        enddo
        !log(op_file_lambda(i)) > log_wavel >  log(op_file_lambda(i-1))
        frac_lambda = (log_wavel-log(op_file_lambda(i-1,pop)))/(log(op_file_lambda(i,pop))-log(op_file_lambda(i-1,pop)))

        fact1 = frac_a_m1 * (1.-frac_lambda)
        fact2 = frac_a_m1 * frac_lambda
        fact3 = frac_a * (1.-frac_lambda)
        fact4 = frac_a * frac_lambda

     else ! Ordre decroisant
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

  !if ((taille_grain == dust_pop(pop)%ind_fin).and.(lambda==n_lambda)) call free_mem_opacity_file()

  !! Matrices de mueller
  if (aniso_method==1) then ! we have to compute a phase function from g

     ! HG avec le g interpole dans la table
     do j=0,nang_scatt
        s11(j)=((1-gsca**2)/(2.0))*(1+gsca**2-2*gsca*cos((real(j))/real(nang_scatt)*pi))**(-1.5)
     enddo

     ! Polarisabilite nulle
     s12=0.0 ; s33 = 0.0 ; s34 = 0.0

     if (scattering_method==1) then
        prob_s11(lambda,taille_grain,0)=0.0
        dtheta = pi/real(nang_scatt)
        do j=2,nang_scatt ! probabilite de diffusion jusqu'a l'angle j, on saute j=0 car sin(theta) = 0
           theta = real(j)*dtheta
           prob_s11(lambda,taille_grain,j)=prob_s11(lambda,taille_grain,j-1)+s11(j)*sin(theta)*dtheta
        enddo

        ! s11 est calculee telle que la normalisation soit: 1.0 (def de la HG)
        ! il y a un soucis numerique quand x >> 1 car la resolution en angle n'est pas suffisante
        ! On rate le pic de diffraction (en particulier entre 0 et 1)
        somme_prob = 1.0
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
           ! La normalisation par default est 1 pour la HG
           if (qsca > 1e-35) then
              norme = 1./qsca
           else ! we don't care there won't be any scattering
              norme = huge_real* 1e-4
           endif
        endif

        tab_s11(j,taille_grain,lambda) = s11(j) / norme
        tab_s12(j,taille_grain,lambda) = s12(j) / norme
        tab_s33(j,taille_grain,lambda) = s33(j) / norme
        tab_s34(j,taille_grain,lambda) = s34(j) / norme
     enddo

  endif ! aniso_method ==1

  return

end subroutine mueller_opacity_file

!**********************************************************************

subroutine update_Stokes(S,u0,v0,w0,u1,v1,w1, M)
  ! Convention astronomique utilisee: angle de position
  ! calcule antihoraire a partir du nord celeste
  !
  ! Francois Menard, Montreal, 15 Fevrier 1989
  !
  ! C. Pinte
  !  - Modif 22/12/03 (C. Pinte) : indice l de taille du grain diffuseur
  !  - Normalisation de l'energie : le photon repart avec l'energie avec laquelle il est entré
  !
  ! Coordonnees des diffuseurs
  !  - x est vers l'observateur
  !  - y et z forme une base droite(dans le bon sens)
  !   - y vers la droite et z vers le haut
  !
  ! Calcul de l'angle omega entre le plan de diffusion
  ! et le nord celeste projete(coord. equatoriale)
  !
  ! Dans le systeme de coord de la nebuleuse
  !    vecteur 1 (du point "0" au point "1") = (u0,v0,w0)
  !    vecteur 2 (du point "1" au point "2") = (u1,v1,w1)
  !
  !  On transforme pour que v2prime = (1,0,0)
  !  L'observateur est alors a +x dans le nouveau systeme

  implicit none

  real(kind=dp), intent(in) ::  u0,v0,w0,u1,v1,w1
  real(kind=dp), dimension(4,4), intent(in) :: M
  real(kind=dp), dimension(4), intent(inout) :: S

  real :: sinw, cosw, omega, theta, costhet,  xnyp, stok_I0, norme, frac_m1
  real(kind=dp) :: v1pi, v1pj, v1pk, S1_0
  real(kind=dp), dimension(4,4) :: ROP, RPO
  real(kind=dp), dimension(4) :: C, D


  ! Transformation pour v1prime
  call rotation(u0,v0,w0,u1,v1,w1,v1pi,v1pj,v1pk)

  !  Calcul des angles pour la rotation
  !  La normale yprime est le produit vectoriel de v1prime x v2prime
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

  ! calcul de l'angle entre la normale et l'axe z (theta)
  theta = acos(costhet)
  if (theta >= pi) theta = 0.0

  ! le plan de diffusion est a +ou- 90deg de la normale
  theta = theta + pi_sur_deux

  ! Dans les matrices de rotation l'angle est omega = 2 * theta
  omega = 2.0 * theta

  !  prochain if car l'arccos va de 0 a pi seulement
  !   Le +/- pour faire la difference dans le sens de rotation
  if (v1pk < 0.0) omega = -1.0 * omega

  ! Calcul des elements des matrices de rotation
  !
  ! RPO = rotation du point vers le systeme original
  ! ROP = rotation du systeme original vers le point
  ! (amene l'axe z dans le plan de diffusion)
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

  ! Normalisation de l'energie : le photon repart avec l'energie avec laquelle il est entré
  ! I sortant = I entrant si diff selon s11  (tab_s11=1.0 normalisé dans mueller2)
  ! I sortant = I entrant * s11 si diff uniforme
  if (S(1) > tiny_real) S(:) = S(:) * M(1,1)*S1_0/S(1)

  return

end subroutine update_Stokes

!**********************************************************************

subroutine get_Mueller_matrix_per_grain(lambda,itheta,frac,taille_grain, M)

  integer, intent(in) :: lambda,itheta,taille_grain
  real, intent(in) :: frac
  real(kind=dp), dimension(4,4), intent(out) :: M

  real :: frac_m1

  frac_m1 = 1.0 - frac

  M(:,:) = 0.0_dp
  M(1,1) = tab_s11(itheta,taille_grain,lambda) * frac + tab_s11(itheta-1,taille_grain,lambda) * frac_m1
  M(2,2) = tab_s22(itheta,taille_grain,lambda) * frac + tab_s22(itheta-1,taille_grain,lambda) * frac_m1
  M(1,2) = tab_s12(itheta,taille_grain,lambda) * frac + tab_s12(itheta-1,taille_grain,lambda) * frac_m1
  M(2,1) = M(1,2)
  M(3,3) = tab_s33(itheta,taille_grain,lambda) * frac + tab_s33(itheta-1,taille_grain,lambda) * frac_m1
  M(4,4) = tab_s44(itheta,taille_grain,lambda) * frac + tab_s44(itheta-1,taille_grain,lambda) * frac_m1
  M(3,4) = -tab_s34(itheta,taille_grain,lambda)* frac - tab_s34(itheta-1,taille_grain,lambda) * frac_m1
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

  M(:,:) = 0.0
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

integer function seuil_n_dif(lambda)

  implicit none

  integer, intent(in) :: lambda
  integer :: n
  real :: albedo
  real, parameter :: seuil = 1.0e-4
  integer, parameter :: seuil_n = 15


  albedo = maxval(tab_albedo_pos(:,lambda))

  n = floor(log(seuil)/log(albedo))+1

  if (n < seuil_n) then
     if (albedo**seuil_n > tiny_real) then
        n=seuil_n
     else
        n = floor(log(tiny_real)/log(albedo))
     endif
  endif

  seuil_n_dif=n

  return

end function seuil_n_dif

!***********************************************************

subroutine isotrope(aleat1,aleat2,u,v,w)
! Choix direction de vol isotrope
! C. Pinte
! 24/05/05

  implicit none

  real, intent(in) :: aleat1, aleat2
  real(kind=dp), intent(out) :: u,v,w

  real(kind=dp) :: SRW02, ARGMT, w02

  w = 2.0_dp*aleat1-1.0_dp
  W02 =  1.0_dp - w*w
  SRW02 = sqrt(W02)
  ARGMT = PI * ( 2.0_dp * aleat2 - 1.0_dp )
  u = SRW02 * cos(ARGMT)
  v = SRW02 * sin(ARGMT)

  return

end subroutine isotrope

!********************************************************************

subroutine hg(g, aleat, itheta, cospsi)
!********************************************************
!* CALCUL DU COSINUS DE L'ANGLE DE DIFFUSION, COS(PSI)
!*
!*     1-PARAMETER HENYEY-GREENSTEIN PHASE FUNCTION
!*
!* FRANCOIS MENARD, 23 NOV 1988, UDEM
! Modif 22/12/03 (C. Pinte) : indice l de taille du grain diffuseur
! - gestion du cas isotrope proprement
! - passage en double necessaire
!********************************************************

  implicit none

  ! Le calcul de cospsi se fait en double precision mais on
  ! renvoie cospsi en simple precision
  ! Passage cospsi en dp lors passage physical_length en dp
  real, intent(in) :: g
  real, intent(in) :: aleat
  real(kind=dp), intent(out) :: cospsi
  integer, intent(out) :: itheta
  real (kind=dp) :: rand, g1, g2

  rand = min(real(aleat,kind=dp), 1.0_dp-1e-6_dp)

  if (abs(g) > tiny_real) then
     g1 = g ! dp
     g2 = g1*g1
     cospsi = (1.0_dp + g2 - ((1.0_dp - g2) / (1.0_dp - g1 + 2.0_dp*g1*rand))**2) / (2.0_dp * g1)
  else ! g=0 --> diffusion isotrope
     cospsi=2.0_dp*rand-1.0_dp
  endif

  itheta = floor(acos(cospsi)*180.0_dp/pi)+1
  if (itheta > nang_scatt) itheta = nang_scatt

  return
end subroutine hg

!***********************************************************

subroutine angle_diff_theta(lambda, taille_grain, aleat, aleat2, itheta, cospsi)
! Calcul du cosinus de l'angle de diffusion
! a partir de l'integrale des s11 pretabulées
! itheta est l'indice i de l'angle de 1 a 180 correspondant a i-0.5°
! cospsi est tire uniformément dans le bin de 1° autour de l'angle i-0.5°
! itheta est utilisé pour les valeurs prétabulées
! cospsi est utilisé pour la direction de vol
! C. Pinte 23/10/04

  implicit none

  integer, intent(in) :: lambda, taille_grain
  real, intent(in) :: aleat, aleat2
  integer, intent(out) :: itheta
  real(kind=dp), intent(out) :: cospsi

  integer :: k, kmin, kmax

  kmin=0
  kmax=nang_scatt
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_s11(lambda,taille_grain,k) < aleat) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   k=kmax

   itheta=k
   !cospsi=cos((real(k)-0.5)*pi/180.)

   ! Tirage aleatoire de l'angle de diffusion autour entre l'angle k et l'angle k-1
   ! diffusion uniforme (lineaire en cos)
   cospsi=cos((real(k)-1.0)*pi/real(nang_scatt)) + &
        aleat2*(cos((real(k))*pi/real(nang_scatt))-cos((real(k)-1.0)*pi/real(nang_scatt)))

   return

end subroutine angle_diff_theta

!**********************************************************************

subroutine angle_diff_theta_pos(lambda, icell, aleat, aleat2, itheta, cospsi)
! Calcul du cosinus de l'angle de diffusion
! a partir de l'integrale des s11 pretabulee par cellule
! itheta est l'indice i de l'angle de 1 a 180 correspondant Ã  i-0.5Â°
! cospsi est tire uniformÃ©ment dans le bin de 1Â° autour de l'angle i-0.5Â°
! itheta est utilise pour les valeurs pretabulee
! cospsi est utilise pour la direction de vol
! C. Pinte 9/01/05

  implicit none

  integer, intent(in) :: lambda,icell
  real, intent(in) :: aleat, aleat2
  integer, intent(out) :: itheta
  real(kind=dp), intent(out) :: cospsi

  integer :: k, kmin, kmax

  kmin=0
  kmax=nang_scatt
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_s11_pos(k,icell,lambda) < aleat) then
        kmin = k
     else
        kmax = k
     endif
     k = (kmin + kmax)/2
   enddo   ! while
   k=kmax

   itheta=k
   !cospsi=cos((real(k)-0.5)*pi/180.)

   ! Tirage aleatoire de l'angle de diffusion entre l'angle k et l'angle k-1
   ! diffusion uniforme (lineaire en cos)
   cospsi=cos((real(k,kind=dp)-1.0_dp)*pi/real(nang_scatt,kind=dp)) + &
        aleat2*(cos((real(k,kind=dp))*pi/real(nang_scatt,kind=dp))-cos((real(k,kind=dp)-1.0_dp)*pi/real(nang_scatt,kind=dp)))

   return

end subroutine angle_diff_theta_pos

!**********************************************************************

subroutine funcd(x,fval,fderiv,ppp,A)
! Calcule la fonction de repartition de l'angle de diffusion phi
! ainsi que sa dérivée pour déterminer le zero par la méthode des
! tangentes (Newton)
! C. Pinte    23/10/2004

  implicit none

  real, intent(in) :: x,ppp,A
  real, intent(out) :: fval,fderiv

  ! T'es sur que c'est le bon signe la ?????????
  fval=2*x-ppp*sin(2*x)-A
  fderiv=2-ppp*2*cos(2*x)

  return

end subroutine funcd

!**********************************************************************

subroutine angle_diff_phi(lambda,taille_grain, I, Q, U, itheta, frac, aleat, phi)
! Tirage de l'angle de diffusion phi du photon
! Uniforme pour une onde non polarisee
! Suivant fct de phase de Mie pour un photon polarise
! C. Pinte    23/10/2004
! Ajout du cas ou les matrices de Mueller sont donnees en entrees
! 20/04/2023

  implicit none

  integer, intent(in) :: lambda,taille_grain, itheta
  real, intent(in) :: frac, I, Q, U, aleat
  real(kind=dp), intent(out) :: phi

  real :: p, pp, ppp, phi1, phi2, frac_m1

  real(kind=dp) :: Q_dp, U_dp, Ip

  frac_m1 = 1.0 - frac

!  write(*,*) 'in',l, itheta, I, Q, U, aleat

  ! Flux polarisé et taux de pola
  Q_dp=Q;U_dp=U
  Ip=sqrt(Q_dp*Q_dp+U_dp*U_dp)
  p=Ip/I

  ! polarisabilite
  pp= (tab_s12(itheta,taille_grain,lambda) * frac + tab_s12(itheta-1,taille_grain,lambda) * frac_m1) &
       / (tab_s11(itheta,taille_grain,lambda) * frac + tab_s11(itheta-1,taille_grain,lambda) * frac_m1)

  ppp=p*pp
!  write(*,*) p,pp,ppp

  if (abs(ppp) > 1.e-3) then
     ! Mesure de l'angle du plan de pola par rapport au Nord céleste
     phi1=0.5*acos(Q_dp/Ip) ! C'est ici qu'on a besoin du dp
     if (U < 0.0) then
        phi1= -phi1
     endif
     ! Tirage de l'angle entre le plan de diffusion et le plan de pola
     ! selon A = 2*phi - ppp*sin(2*phi)   A=4*pi*aleat
     phi2=rtsafe(funcd,0.,2*real(pi),1.e-4,ppp,4*real(pi)*aleat)
!     write(*,*) 'test', phi1 , Q/Ip
     phi=phi1+phi2
     if (phi > pi) phi = phi -2*pi
     if (phi < -pi) phi = phi +2*pi
  else ! Tirage uniforme
     phi= pi * (2._dp*aleat -1.0_dp)
  endif
!  write(*,*) phi/pi

  return

end subroutine angle_diff_phi

!**********************************************************************

real function rtsafe(funcd,x1,x2,xacc,ppp,A)
! Trouve le zéro d'une fonction
! Ici, l'angle de diffusion en phi pour un photon polarisé
! suivant la fonction de répartion donnée par funcd
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
  if (abs(wavelength - tab_lambda(1)) < 1.0e-5) call error("walength does correspond to wavelength of the Mueller matrix")

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
