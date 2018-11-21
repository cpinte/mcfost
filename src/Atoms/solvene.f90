! ---------------------------------------------------------
! Module that solves for electron density for given
! temperature grid, Hydrogen total populations
! elements and their abundance and their LTE or NLTE populations.
!
!
! If ne_fromscratch is true then first guess is computed using pure Hydrogen
! ionisation (or +Helium)
! else the value in ne is used as first guess.
! This allows to iterate the electron density through the NLTE scheme
!
! See: Hubeny & Mihalas (2014)
!         "Theory of Stellar Atmospheres",
!                         from p. 588 to p. 593
!
!
! ---------------------------------------------------------


MODULE solvene

 use atmos_type, only : atmos, Nelem !in atmos%Elements
 use atom_type, only : Element, AtomType
 use math, only : interp1D
 use constant
 use lte
 !use accelerate, only : initNg, freeNg, NgAcceleration

 IMPLICIT NONE

 real(8), parameter :: MAX_ELECTRON_ERROR=1e-5
 integer, parameter :: N_MAX_ELECTRON_ITERATIONS=30
 integer, parameter :: N_MAX_ELEMENT=26 !100 is the max

 CONTAINS

 ! ---------------------------------------------------------
 ! do not forget to use the good value of parition function and potential
 ! parition function read in logarithm and used is 10**(U)
 ! potential in cm-1, converted in the routines in J.
 ! ---------------------------------------------------------


 SUBROUTINE ne_Hionisation (k, U0, U1, ne)
 ! ---------------------------------------------------------
 !Application of eq. 4.35 of Hubeny & Mihalas to ionisation
 !of H.
 !Njl = Nj1l * ne * phi_jl
 !ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
 !eq. 17.77 and 17.78
 !ne(H) = (sqrt(N*phiH + 1)-1)/phiH
 ! ---------------------------------------------------------

  integer, intent(in) :: k
  real(8), intent(in) :: U0, U1
  real(8) :: phiH
  real(8), intent(out) :: ne
  phiH = phi_jl(k, U0, U1, atmos%Elements(1)%ionpot(1))

  !ne = (sqrt(atmos%nHtot(k)*phiH*4. + 1)-1)/(2.*phiH)
  ne = (sqrt(atmos%nHtot(k)*phiH + 1)-1)/(phiH)

 RETURN
 END SUBROUTINE ne_Hionisation

 SUBROUTINE ne_Metal(k, U0, U1, chi, A, ne)
 ! ---------------------------------------------------------
 !Application of eq. 4.35 of Hubeny & Mihalas to ionisation
 !of a single metal.
 ! ---------------------------------------------------------

  integer, intent(in) :: k
  real(8), intent(in) :: U0, U1, chi, A
  real(8) :: phiM, alphaM
  real(8), intent(out) :: ne

  phiM = phi_jl(k, U0, U1, chi)
  alphaM = A ! relative to H, for instance 1.-6 etc
  ne = (sqrt(alphaM*atmos%nHtot(k)*phiM +0.25*(1+alphaM)**2)&
     - 0.5*(1+alphaM) ) / phiM
 RETURN
 END SUBROUTINE ne_Metal


 FUNCTION getPartitionFunctionk(elem, stage, k) result (Uk)
  type(Element) :: elem
  integer, intent(in) :: stage, k
  real(8) :: Uk

  Uk = Interp1D(atmos%Tpf,elem%pf(stage,:),atmos%T(k))
       !do not forget that Uk is base 10 logarithm !!
       ! note that in RH, natural (base e) logarithm
       ! is used instead
  Uk = (10.d0)**(Uk)

 RETURN
END FUNCTION getPartitionFunctionk



 SUBROUTINE getfjk (Elem, ne, k, fjk, dfjk)
 ! ---------------------------------------------------------
 ! fractional population f_j(ne,T)=N_j/N for element Elem
 ! and its partial derivative with ne
 ! ---------------------------------------------------------

  real(8), intent(in) :: ne
  integer, intent(in) :: k
  type (Element), intent(in) :: Elem
  type (AtomType) :: atom
  real(8), dimension(:), intent(inout) :: fjk, dfjk
  real(8) :: Uk, Ukp1, sum1, sum2
  logical :: is_active=.false.
  integer :: nll, j, i

  ! check if the element as an atomic model and it is active
  do nll=1,atmos%Nactiveatoms
   if (Elem%ID.eq.atmos%Atoms(nll)%ID .and. &
       atmos%Atoms(nll)%active) then
     is_active=.true.
     write(*,*) "Atom ",Elem%ID,atmos%Atoms(nll)%ID," is active"
     exit
   end if
  end do

  if (is_active) then
   atom = Elem%model
   do j=1,Elem%Nstage
    fjk(j) = 0.
    dfjk(j) = 0.
   end do
   do i=1,atom%Nlevel
    fjk(atom%stage(i))=fjk(atom%stage(i))+atom%stage(i)*atom%n(i,k)
   end do
   do j=1,Elem%Nstage
    fjk(j) = fjk(j)/atom%ntotal(k)
   end do
  else !not active, whateveeer, use LTE
   fjk(1)=1.
   dfjk(1)=0.
   sum1 = 1.
   sum2 = 0.
   Uk = getPartitionFunctionk(elem,1,k)
   do j=2,Elem%Nstage
    Ukp1 = getPartitionFunctionk(elem,j,k)
    ! fjk(j) = fjk(j-1)/(phi*ne) Saha equation
    ! Nj = Nj-1/(phi*ne)
    ! fj = Nj/N = Nj-1/N / (phi*ne)
    ! See LTEpops_elem in LTE.f90 for further details
    fjk(j) = Sahaeq(k,fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)
    !write(*,*) "j=",j," fjk=",fjk(j)
    dfjk(j) = -(j-1)*fjk(j)/ne
    !write(*,*) "j=",j," dfjk=",dfjk(j)
    sum1 = sum1 + fjk(j)
    sum2 = sum2 + dfjk(j)
    Uk = Ukp1
   end do
   do j=1,elem%Nstage
    fjk(j)=fjk(j)/sum1
    dfjk(j)=(dfjk(j)-fjk(j)*sum2)/sum1
    !write(*,*) j, fjk(j), dfjk(j)
   end do
  end if

 RETURN
 END SUBROUTINE getfjk

 SUBROUTINE SolveElectronDensity(ne, ne_initial_solution)
  ! ---------------------------------------------------------
  ! Solve for electron density for a set of elements
  ! stored in atmos%Elements. Elements up to N_MAX_ELEMENT are
  ! used. If an element has also a atomic model, and if for
  ! this element NLTE populations are known these pops are
  ! used to compute the electron density. Otherwise, LTE is used.
  ! when ne_initial_solution is set to HIONISA, use
  ! sole hydgrogen ionisation to estimate the initial ne density.
  ! is set to NPROTON or NEMODEL, protons number or density
  ! read from the model is used. Note that NPROTON suppose
  ! that NLTE populations are present for hydrogen, since
  ! nprot = hydrogen%n(Nlevel,:).
  ! If keyword is not set, HIONISATION is used.
  ! ---------------------------------------------------------

  real(8), dimension(atmos%Nspace), intent(inout) :: ne
  character(len=7), optional :: ne_initial_solution
  character(len=7) :: initial
  real(8) :: error, ne_old, akj, sum, Uk, dne, Ukp1
  real(8) :: ne_oldM, UkM, PhiHmin
  real(8), dimension(atmos%Nspace) :: np
  real(8), dimension(:), allocatable :: fjk, dfjk
  integer :: Nmaxstage=0, n, k, niter, j

  if (.not. present(ne_initial_solution)) then
      initial="HIONISA"!use HIONISAtion
  else
    initial=ne_initial_solution
  end if

  do n=1,Nelem
   Nmaxstage=max(Nmaxstage,atmos%Elements(n)%Nstage)
  end do
  allocate(fjk(Nmaxstage))
  allocate(dfjk(Nmaxstage))

  !np is the number of protons, by default the last level
  !of Hydrogen

  if (initial.eq."NPROTON") &
     np=Hydrogen%n(Hydrogen%Nlevel,:)

  do k=1,atmos%Nspace
   if (.not.atmos%lcompute_atomRT(k)) CYCLE
   !if ((atmos%nHtot(k)==0d0).or.(atmos%T(k)==0d0)) CYCLE ! go to next depth point

   if (initial.eq."NPROTON") then
    ne_old = np(k)
   else if (initial.eq."NEMODEL") then
    ne_old = atmos%ne(k)
   else
    !Initial solution ionisation of H
    Uk = getPartitionFunctionk(atmos%Elements(1), 1, k)
    Ukp1 = 1d0 !getPartitionFunctionk(atmos%Elements(1), 2, k)
    CALL ne_Hionisation (k, Uk, Ukp1, ne_old)
    Uk = getPartitionFunctionk(atmos%Elements(26), 1, k)
    Ukp1 = getPartitionFunctionk(atmos%Elements(26), 2, k)
    CALL ne_Metal(k, Uk, Ukp1, atmos%elements(26)%ionpot(1), &
         atmos%elements(26)%Abund, ne_oldM)
    !write(*,*) "neMetal=",ne_oldM
    !if Abund << 1. and chiM << chiH then
    ! ne (H+M) = ne(H) + ne(M)
    ne_old = ne_old + ne_oldM
   end if
   !write(*,*) "k=",k," ne_old=",ne_old, &
   !       " ne_mod=",atmos%ne(k)

   ne(k) = ne_old
   niter=0
   do while (niter.lt.N_MAX_ELECTRON_ITERATIONS)
    error = ne_old/atmos%nHtot(k)
    sum = 0.

    do n=1,Nelem
     if (n.gt.N_MAX_ELEMENT) exit

     CALL getfjk(atmos%Elements(n),ne_old,k,fjk,dfjk)

     if (n.eq.1)  then ! H minus
       PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
       ! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
       error = error + ne_old*fjk(1)*PhiHmin
       sum = sum-(fjk(1)+ne_old*dfjk(1))*PhiHmin
       !write(*,*) "phiHmin=",PhiHmin,error, sum
     end if
     do j=2,atmos%elements(n)%Nstage
      akj = atmos%elements(n)%Abund*(j-1) !because j starts at 1
      error = error -akj*fjk(j)
      sum = sum + akj*dfjk(j)
      !write(*,*) n-1, j-1, akj, error, sum
     end do
    end do !loop over elem
    ne(k) = ne_old - atmos%nHtot(k)*error /&
          (1.-atmos%nHtot(k)*sum)
    dne = dabs((ne(k)-ne_old)/ne_old)
    ne_old = ne(k)


    if (dne.le.MAX_ELECTRON_ERROR) then
      !write(*,*) "icell=",k," T=",atmos%T(k)," ne=",ne(k)
     exit
    !else
      !write(*,*) "icell=",k," T=",atmos%T(k)," nH=",atmos%nHtot(k), &
      !         "dne = ",dne, " ne=",ne(k)
    end if
    niter = niter + 1
   end do !while loop
  end do !loop over spatial points

  write(*,*) "maxium electron density (m^-3) =", MAXVAL(atmos%ne), &
            " minimum electron density (m^-3) =", MINVAL(atmos%ne)
  deallocate(fjk, dfjk)
 RETURN
 END SUBROUTINE SolveElectronDensity

END MODULE solvene
