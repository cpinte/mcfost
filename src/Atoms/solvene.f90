! ------------------------------------------------------------------- !
! ------------------------------------------------------------------- !
! Module that solves for electron density for given
! temperature grid, Hydrogen total populations
! elements and their abundance and their LTE or NLTE populations.
!
!
! If ne_intial_slution is "NEMODEL" then first guess is computed using 
! the value stored in atmos%ne as first guess.
! This allows to iterate the electron density through the NLTE scheme.
!
! See: Hubeny & Mihalas (2014)
!         "Theory of Stellar Atmospheres",
!                         from p. 588 to p. 593
!
!
! ------------------------------------------------------------------- !
! ------------------------------------------------------------------- !

MODULE solvene

 use atmos_type, only : atmos, Nelem !in atmos%Elements
 use atom_type, only : Element, AtomType
 use math, only : interp1D
 use constant
 use lte
 !use accelerate, only : initNg, freeNg, NgAcceleration
 
 !$ use omp_lib

 IMPLICIT NONE

 double precision, parameter :: MAX_ELECTRON_ERROR=1e-5
 integer, parameter :: N_MAX_ELECTRON_ITERATIONS=30
 integer, parameter :: N_MAX_ELEMENT=26 !100 is the max

 CONTAINS

 ! ----------------------------------------------------------------------!
 ! do not forget to use the good value of parition function and potential
 ! parition function read in logarithm and used is 10**(U)
 ! potential in cm-1, converted in the routines in J.
 ! ----------------------------------------------------------------------!


 SUBROUTINE ne_Hionisation (k, U0, U1, ne)
 ! ----------------------------------------------------------------------!
  ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
  ! of H.
  ! Njl = Nj1l * ne * phi_jl
  ! ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
  ! eq. 17.77 and 17.78
  ! ne(H) = (sqrt(N*phiH + 1)-1)/phiH
 ! ----------------------------------------------------------------------!

  integer, intent(in) :: k
  double precision, intent(in) :: U0, U1
  double precision :: phiH
  double precision, intent(out) :: ne
  phiH = phi_jl(k, U0, U1, atmos%Elements(1)%ionpot(1))

  !ne = (sqrt(atmos%nHtot(k)*phiH*4. + 1)-1)/(2.*phiH)
  ne = (sqrt(atmos%nHtot(k)*phiH + 1)-1)/(phiH)

 RETURN
 END SUBROUTINE ne_Hionisation

 SUBROUTINE ne_Metal(k, U0, U1, chi, A, ne)
 ! ----------------------------------------------------------------------!
  ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
  ! of a single metal.
 ! ----------------------------------------------------------------------!

  integer, intent(in) :: k
  double precision, intent(in) :: U0, U1, chi, A
  double precision :: phiM, alphaM
  double precision, intent(out) :: ne

  phiM = phi_jl(k, U0, U1, chi)
  alphaM = A ! relative to H, for instance 1.-6 etc
  ne = (sqrt(alphaM*atmos%nHtot(k)*phiM +0.25*(1+alphaM)**2)&
     - 0.5*(1+alphaM) ) / phiM
 RETURN
 END SUBROUTINE ne_Metal


 FUNCTION getPartitionFunctionk(elem, stage, k) result (Uk)
 ! ----------------------------------------------------------------------!
  ! Interpolate the partition function of Element elem in ionisation stage
  ! stage at cell point k.
 ! ----------------------------------------------------------------------!

  type(Element) :: elem
  integer, intent(in) :: stage, k
  double precision :: Uk

  Uk = Interp1D(atmos%Tpf,elem%pf(stage,:),atmos%T(k))
       !do not forget that Uk is base 10 logarithm !!
       ! note that in RH, natural (base e) logarithm
       ! is used instead
  Uk = (10.d0)**(Uk)

 RETURN
END FUNCTION getPartitionFunctionk



 SUBROUTINE getfjk (Elem, ne, k, fjk, dfjk)
 ! ----------------------------------------------------------------------!
  ! fractional population f_j(ne,T)=N_j/N for element Elem
  ! and its partial derivative with ne. If Elem is an element with
  ! detailed model and if NLTE populations for it exist, their are used
  ! instead of LTE.
 ! ----------------------------------------------------------------------!

  double precision, intent(in) :: ne
  integer, intent(in) :: k
  type (Element), intent(in) :: Elem
  type (AtomType) :: atom
  double precision, dimension(:), intent(inout) :: fjk, dfjk
  double precision :: Uk, Ukp1, sum1, sum2
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
   fjk = 0d0
   dfjk = 0d0

   !For Nlevel, Nlevel stages
   fjk(atom%stage(:))=fjk(atom%stage(:))+atom%stage(:)*atom%n(:,k)
   !For Nstage
   fjk(:) = fjk(:)/atom%ntotal(k)

  else !not active, whateveeer, use LTE
   fjk(1)=1.
   dfjk(1)=0.
   sum1 = 1.
   sum2 = 0.
   Uk = getPartitionFunctionk(elem,1,k)
   do j=2,Elem%Nstage !-->vectorize ?
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
   fjk(:)=fjk(:)/sum1
   dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1
  end if

 RETURN
 END SUBROUTINE getfjk

 SUBROUTINE SolveElectronDensity(ne_initial_solution)
 ! ----------------------------------------------------------------------!
  ! Solve for electron density for a set of elements
  ! stored in atmos%Elements. Elements up to N_MAX_ELEMENT are
  ! used. If an element has also an atomic model, and if for
  ! this element NLTE populations are known these pops are
  ! used to compute the electron density. Otherwise, LTE is used.
  ! 
  ! When ne_initial_solution is set to HIONISA, uses
  ! sole hydgrogen ionisation to estimate the initial ne density.
  ! If set to NPROTON or NEMODEL, protons number or density
  ! read from the model are used. Note that NPROTON supposes
  ! that NLTE populations are present for hydrogen, since
  ! nprot = hydrogen%n(Nlevel,:).
  ! If keyword is not set, HIONISATION is used.
 ! ----------------------------------------------------------------------!
  character(len=20), optional :: ne_initial_solution
  character(len=20) :: initial
  double precision :: error, ne_old, akj, sum, Uk, dne, Ukp1
  double precision :: ne_oldM, UkM, PhiHmin
!   double precision, dimension(atmos%Nspace) :: np
  double precision, dimension(:), allocatable :: fjk, dfjk
  integer :: Nmaxstage=0, n, k, niter, j

  if (.not. present(ne_initial_solution)) then
      initial="H_IONISATION"!use HIONISAtion
  else
    initial=ne_initial_solution
!     if ((initial .neq. "N_PROTON") .neq. (initial /= "NE_MODEL")) then
!      write(*,*) 'NE initial solution unkown, set to H_IONISATION'
!      initial = "H_IONISATION"
!     end if
  end if

  do n=1,Nelem
   Nmaxstage=max(Nmaxstage,atmos%Elements(n)%Nstage)
  end do
  allocate(fjk(Nmaxstage))
  allocate(dfjk(Nmaxstage))

  !np is the number of protons, by default the last level
  !of Hydrogen.
  
  !Note that will raise errors if np is 0d0 because Hydrogen pops are not
  !known. np is known if:
  ! 1) NLTE populations from previous run are read
  ! 2) The routine is called in the NLTE loop meaning that all H levels are known

!   if (initial.eq."N_PROTON") &
!      np=Hydrogen%n(Hydrogen%Nlevel,:)

  !$omp parallel &
  !$omp default(none) &
  !$omp private(k,n,j,fjk,dfjk,ne_old,niter,error,sum,PhiHmin,Uk,Ukp1,ne_oldM) &
  !$omp private(dne, akj) &
  !$omp shared(atmos, initial,Hydrogen)
  !$omp do
  do k=1,atmos%Nspace
   if (.not.atmos%lcompute_atomRT(k)) CYCLE

   !write(*,*) "The thread,", omp_get_thread_num() + 1," is doing the cell ", k
   if (initial.eq."N_PROTON") then
    ne_old = Hydrogen%n(Hydrogen%Nlevel,k)!np(k)
   else if (initial.eq."NE_MODEL") then
    ne_old = atmos%ne(k)
   else !"H_IONISATION" or unkown
   
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

   atmos%ne(k) = ne_old
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
    atmos%ne(k) = ne_old - atmos%nHtot(k)*error /&
          (1.-atmos%nHtot(k)*sum)
    dne = dabs((atmos%ne(k)-ne_old)/ne_old)
    ne_old = atmos%ne(k)


    if (dne.le.MAX_ELECTRON_ERROR) then
      !write(*,*) "icell=",k," T=",atmos%T(k)," ne=",atmos%ne(k)
     exit
    !else
      !write(*,*) "icell=",k," T=",atmos%T(k)," nH=",atmos%nHtot(k), &
      !         "dne = ",dne, " ne=",atmos%ne(k)
    end if
    niter = niter + 1
   end do !while loop
  end do !loop over spatial points
  !$omp end do
  !$omp end parallel

  deallocate(fjk, dfjk)
  
  write(*,*) "Maximum/minimum Electron density in the model (m^-3):"
  write(*,*) MAXVAL(atmos%ne),MINVAL(atmos%ne,mask=atmos%ne>0)
 RETURN
 END SUBROUTINE SolveElectronDensity

END MODULE solvene
