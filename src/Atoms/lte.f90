! -------------------------------------------------------------- !
! -------------------------------------------------------------- !
   ! Various routines to solve for LTE populations
! -------------------------------------------------------------- !
! -------------------------------------------------------------- !


MODULE lte

 use atom_type, only : AtomType, element
 use atmos_type, only : atmos, Hydrogen, Helium, ntotal_atom
 use constant
 use math
 use writeatom, only : writeHydrogenMinusDensity

 use constantes, only : tiny_dp, huge_dp

 !$ use omp_lib

 IMPLICIT NONE

 real(kind=dp), parameter :: phi_min_limit = 1d-50 !tiny_dp

 CONTAINS

 FUNCTION phi_jl(k, Ujl, Uj1l, ionpot) result(phi)
  ! -------------------------------------------------------------- !
  ! Hubeny & Mihalas (2014)
  ! "Theory of Stellar Atmospheres", p.  94, eq. 4.35
  !
  !
  ! ionisation state j of element l ! (j1=j+1)
  ! ionpot = ionpot of ion j of element l (ionpot0l(H) = 13.6 eV)
  ! ionpot in J
  !
  ! Njl = Nj1l * ne * phi_jl
  ! --> SahaEq: Nj1l = Njl / (phi_jl * ne)
  ! --> NH- = NH * ne * phi_jl
  ! NH = NH+ * ne * phi_jl
  ! -------------------------------------------------------------- !
  integer :: k
  real(kind=dp) :: ionpot, C1, phi
  real(kind=dp) :: Ujl, Uj1l, expo
  C1 = 0.5*(HPLANCK**2 / (2.*PI*M_ELECTRON*KBOLTZMANN))**1.5
  !!ionpot = ionpot * 100.*HPLANCK*CLIGHT !cm1->J
  !! -> avoiding dividing by big numbers causing overflow.
  !!maximum argument is 600, exp(600) = 3.77e260
  expo = dexp(min(600d0,ionpot/(KBOLTZMANN*atmos%T(k))))
  if (ionpot/(KBOLTZMANN*atmos%T(k)) >= 600d0) expo = huge_dp

  !if dexp(300) it means phi is "infinity", exp(300) == 2e130 so that's enough
  phi = C1 * (Ujl / Uj1l) * expo / (atmos%T(k)**1.5 + tiny_dp)
  if (phi < phi_min_limit) phi = 0d0 !tiny_dp ! or phi = 0d0 should be more correct ?
  								   ! but test to not divide by 0

  if (is_nan_infinity(phi)>0) then
   write(*,*) "error, phi_jl", k, phi
   stop
  endif

  RETURN
  END FUNCTION phi_jl

  FUNCTION BoltzmannEq9dot1(k, Ei, gi, UI) result(ni_NI)
  ! -------------------------------------------------------------- !
   ! Hubeny & Mihalas (2014)
   ! "Theory of Stellar Atmospheres" eq. 9.1
   !
   ! ni_NI = (gi/UI)*exp(-Ei/kT)
   ! OR
   ! nipjl/nijl = gipjl/gijl * exp(-(Eip-Ei)/kT)
   !
   ! Ei = level energy measured from ground state in J
   !
   ! UI = partition function of stage I (for instance the one of HI)
   !
   ! gi = statistical weight of level i of atom in stage I
   !
   ! NI represents the total number density of stage I
   ! ni population of level i belonging to NI
  ! -------------------------------------------------------------- !

   real(kind=dp) :: Ei, gi, Ui, ni_NI, kT
   integer :: k

   !Ei = Ei * 100.*HPLANCK*CLIGHT
   kT = KBOLTZMANN * atmos%T(k)
   ni_NI = (gi/UI) * exp(-Ei/kT)
  RETURN
  END FUNCTION BoltzmannEq9dot1

  FUNCTION BoltzmannEq4dot20b(k, Ei, gi, gi1) result(ni1_ni)
  ! -------------------------------------------------------------- !
   ! Hubeny & Mihalas (2014)
   ! "Theory of Stellar Atmospheres" eq. 4.20b
   !
   ! nipjl/nijl = gipjl/gijl * exp(-(Eip-Ei)/kT)
   !
   ! Ei = level energy of level i measured from ground state in J
   ! Ei = Eij - E0j
   !
   ! gi = statistical weight of level i of atom in stage I
   !
   ! ni population of level i belonging to NI
  ! -------------------------------------------------------------- !

  real(kind=dp) :: Ei,ni1_ni, kT
  real(kind=dp) ::  gi, gi1
  integer :: k

  kT = KBOLTZMANN * atmos%T(k)
  ni1_ni = (gi1/gi) * exp(-Ei/kT)
  RETURN
  END FUNCTION BoltzmannEq4dot20b

  FUNCTION SahaEq(k, NI, UI1, UI, chi, ne) result(NI1)
  ! -------------------------------------------------------------- !
  ! Hubeny & Mihalas (2014)
  ! "Theory of Stellar Atmospheres" eq. 9.2
  ! (here, the inverse, i.e., NI+1/NI)
  !
  ! NI1 = NI * 2.*(h**2/(2pimkT))^-1.5 * &
  !                   UI1_UI * np.exp(-chiI/kT) / ne
  ! NI1 = NI/(ne*phijl)
  ! where NI, NI1 are the total number of particles of an element
  ! in the stage I, I+1
  ! -------------------------------------------------------------- !
   integer :: k
   real(kind=dp) :: NI, UI1, UI, chi, ne, phi, NI1
   phi = phi_jl(k, UI, UI1, chi) ! chi in J

    !phi should be between phi_min_limit and dexp(600)
    !and ne between ne_limit and nemax, so never nan nor infinity
    !further in General, a larg phi implies a ne close to 0, so the product remains small

   if ((ne > 0).and.(is_nan_infinity(ne)==0)) then !is this correct or if phi/ne nan or infinty ?
    NI1 = NI/(phi*ne)
   else !all in ground state, phi->infinity if T->0 and phi->0 if T->infinty
   		!AND IF	ne -> 0, T likely goes to 0, all in neutral states
   		!meaning if ne->0 phi goes to infinity
    NI1 = 0d0
   endif

  RETURN
  END FUNCTION SahaEq

  SUBROUTINE calc_nHmin()
   integer :: k, j
   real(kind=dp) :: phiHmin
   write(*,*) " BEWARE, nHMIN NOT USED anymore. WAITING FOR CHEMECHIL SOLVER"
   return
   !$omp parallel &
   !$omp default(none) &
   !$omp private(k, PhiHmin,j) &
   !$omp shared(atmos,Hydrogen)
   !$omp do
   do k=1, atmos%Nspace
    if (atmos%icompute_atomRT(k) /= 1) CYCLE

    PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
     != 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
    !nHmin is proportial to the total number of neutral Hydrogen
    ! which for hydrogen is given by the sum of the indivdual
    !levels except the last one which is the protons (or H+)
    !number density. Remember:
    ! -> Njl = Nj1l * ne * phi_jl with j the ionisation state
    !j = -1 for Hminus (l=element=H)
    atmos%nHmin = atmos%ne(k) * PhiHmin * Hydrogen%n(1,k)
    !atmos%nHmin(k) = sum(Hydrogen%n(1:Hydrogen%Nlevel-1,k)) * atmos%ne(k)*PhiHmin!faster?
    !! Estimation
   end do !over depth points
   !$omp end do
   !$omp  end parallel
   write(*,*) "Maximum/minimum H minus density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHmin), MINVAL(atmos%nHmin,mask=atmos%icompute_atomRT > 0)
  RETURN
  END SUBROUTINE calc_nHmin

 SUBROUTINE LTEpops(atom, debeye)
  ! -------------------------------------------------------------- !
   ! Computes LTE populations of each level of the atom.
   ! Distributes populations among each level (ni) according
   ! to the total population (NI) of each stage (I),
   ! for which a level i belongs (Boltzmann's equation)
  ! -------------------------------------------------------------- !
  type (Atomtype), intent(inout) :: atom
  logical, intent(in) :: debeye
  real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
  integer :: Z, dZ, k, i, m
  integer, allocatable, dimension(:) :: nDebeye

  !if (Debeye) write(*,*) "Debeye Shielding activated"
  !if (.not.Debeye) write(*,*) "Debeye Shielding deactivated"

  ! Debeye shielding activated:
  ! It lowers the ionisation potential

  ! dE_ion = -Z*(qel**2 / 4PI*EPS0) / D in J
  ! D = sqrt(EPS0/2qel**2) * sqrt(kT/ne) in m

  ! Hubeny & Mihalas (2014)
  ! "Theory of Stellar Atmospheres" p. 244-246

   if (Debeye) then
    c2 = dsqrt(8.0*PI/KBOLTZMANN) * &
      (SQ(Q_ELECTRON)/(4.0*PI*EPSILON_0))**1.5
   !write(*,*) "c2=",c2,  atom%Abund * atmos%nHtot(1),  atom%Abund * atmos%nHtot(2)
    allocate(nDebeye(atom%Nlevel))
    nDebeye(:)=0
    do i=2,atom%Nlevel !reject ground level
     Z = atom%stage(i) !0 for neutrals, 1 for singly ionised
     do m=1,atom%stage(i)-atom%stage(1) !relative to ground level
      nDebeye(i) = nDebeye(i) + Z
      Z = Z + 1
     end do
    end do
   end if

   !$omp parallel &
   !$omp default(none) &
   !$omp private(k, phik, dEion,i,dZ,dE,m,sum,phiHmin) &
   !$omp shared(atmos,Hydrogen,c2,atom, Debeye,nDebeye)
   !$omp do
   do k=1,atmos%Nspace
    if (atmos%icompute_atomRT(k) /= 1) CYCLE
    if (Debeye) dEion = c2*sqrt(atmos%ne(k)/atmos%T(k))
    sum = 1.
    phik = atmos%ne(k)*phi_jl(k,1.d0,1.d0,0.d0)
    !a constant of the temperature and electron density

    do i=2, atom%Nlevel
     dE = atom%E(i) - atom%E(1) ! relative to ground level
     dZ = atom%stage(i) - atom%stage(1) !stage is 0 for neutrals

!    write(*,*) "dZ = ", dZ," nDebeye*dEion=",nDebeye(i)*dEion*&
!         JOULE_TO_EV, " dEion=",dEion*JOULE_TO_EV

     if (Debeye) dE = dE - nDebeye(i)*dEion !J

     ! --------- Boltzmann equation ------------------------------------------- !
     ! relative to ni=n1, the populations
     ! of the ground state in stage stage(1).

     atom%nstar(i,k)=BoltzmannEq4dot20b(k, dE, atom%g(1), atom%g(i))
!     write(*,*) "Atom=",atom%ID, " A=", atom%Abund
     !!relative to n1j(1)
!     write(*,*) k-1, "level=", i-1, " stage=", atom%stage(i), &
!                "n(i,k)=",atom%nstar(i,k)

     ! ---------- Saha equation ------------------------------------------------ !
     ! n1j+1 = n1j * exp(-ionpotj->j+1)/(ne*phi(T)) !or phik here
     ! we can express all pops nij in specific ion stage j
     ! wrt the groud state and ionisation stage of the ground state
     ! namely n1j=1. Remembering that E0j+1-E0j = chiIj->j+1 = ionpot

     ! for this level i in stage stage(i),
     ! we have dZ +1 ion stages from the
     ! level 1 in stage stage(1).
     ! the population of the ground state in stage (i)+dZ
     ! is given by nstar(i,k)/(phik)**dZ.
     ! because expressing dE has Eij - E00 will involved succesive division by phik
     ! dZ + 1 times -> 0 times if same ion stage, +1 times if next ion, +2 times if
     ! third ion ect ..

     ! is equivalent to nstar(i,k)*(phik)**(-dZ)
     ! if dZ = 0, same ion, no division, only Boltzmann eq.
     ! if dZ = 1, dZ+1=2 -> another (next) ion, one div by phik = Saha
     ! if dZ = 2, dZ+1=3 -> the next next ion, two div by phik ...
     do m=2,dZ+1
      if (atmos%ne(k)>0) then
       atom%nstar(i,k)=atom%nstar(i,k)/phik !if phik -> means large n(i) ??
      else
       atom%nstar(i,k) = 0d0
      endif
     end do

     ! by doing so we avoid the easiest case of using
     ! partition function.
     sum = sum+atom%nstar(i,k) !compute total pop = Sum_jSum_i nij
     						   !with Sum_i nij = Nj
    end do
    !Incorpore Hminus
!     if (atom%ID=="H") then
!      phiHmin = atmos%ne(k)*phi_jl(k, 1d0, 2d0, E_ION_HMIN)
!      atmos%nHmin(k) = 1d0 * phiHmin!ground state of HI wrt to nH00
!      do i = 2, atom%Nlevel-1
!       atmos%nHmin(k) = atmos%nHmin(k) + atom%nstar(i,k) * phiHmin !wrt nH00
!      end do
!      sum = sum + atmos%nHmin(k)
!     end if

     ! now we have, 1 + n2/n1 + ni/n1 = sum over all levels
     ! and each stage for each level,
     ! wrt the first population of the first level
     ! further, n1 + n2+n3+nij sum over all stage of each level
     ! gives ntotal.
     ! Therefore n1 = atom%ntotal/sum
!     write(*,*) "-------------------"
!     write(*,*) "Atom=",atom%ID, " A=", atom%Abund
!     write(*,*) "ntot", atom%ntotal(k), " nHtot=",atmos%nHtot(k)
    atom%nstar(1,k) = atom%Abund*atmos%nHtot(k)/sum

    !test positivity, can be 0
    if (atom%nstar(1,k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning too small ground state population ", atom%ID, "n0=", atom%nstar(1,k)
       write(*,*) "cell=",k, atom%ID, "dark?=",atmos%icompute_atomRT(k), "T=",atmos%T(k), "nH=",atmos%nHtot(k), "ne=",atmos%ne(k)
        !atom%nstar(1,k) = tiny_dp
       write(*,*) " ************************************* "
       stop !only if we test >= 0
    end if
    !write(*,*) "depth index=",k, "level=", 1, " n(1,k)=",atom%nstar(1,k)
    do i=2,atom%Nlevel !debug
      atom%nstar(i,k) = atom%nstar(i,k)*atom%nstar(1,k)
      !to avoid very low populations in calculations
      !if (atom%nstar(i,k) < 1d-30) atom%nstar(i,k) = 1d-30
      !write(*,*) "depth index=",k, "level=", i, " n(i,k)=",atom%nstar(i,k)

      !! check positivity. Can be 0 now, see above and phik and phi_jl
      if (atom%nstar(i,k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning population of atom ", atom%ID, "lvl=", i, "nstar=",atom%nstar(i,k), " lower than", &
        " tiny_dp."! Replacing by tiny_dp"
       write(*,*) "cell=",k, atom%ID, "dark?=",atmos%icompute_atomRT(k), "T=",atmos%T(k), "nH=",atmos%nHtot(k), "ne=",atmos%ne(k), &
       			  " n0=", atom%nstar(1,k)
        !atom%nstar(i,k) = tiny_dp
       write(*,*) " ************************************* "
      stop !only if we test >= 0
      end if
    end do

    if (maxval(atom%nstar(:,k)) >= huge_dp) then
    write(*,*) " ************************************* "
     write(*,*) "ERROR, populations of atom larger than huge_dp"
     write(*,*) "cell=",k, atom%ID, "dark?=",atmos%icompute_atomRT(k), "T=",atmos%T(k), "nH=",atmos%nHtot(k), "ne=",atmos%ne(k)
     write(*,*) "nstar=",atom%nstar(:,k)
     write(*,*) " ************************************* "
     stop
    end if

    !faster ?
!    atom%nstar(2:atom%Nlevel,k) = atom%nstar(2:atom%Nlevel,k) * atom%nstar(1,k)
	!!Incorpore Hminus
!     if (atom%ID=="H") atmos%nHmin(k) = atmos%nHmin(k) * atom%nstar(1,k)

   end do !over depth points
   !$omp end do
   !$omp  end parallel
!    write(*,*) "-------------------"

  if (allocated(nDebeye)) deallocate(nDebeye)
 RETURN
 END SUBROUTINE LTEpops


 SUBROUTINE setLTEcoefficients ()
  ! -------------------------------------------------------------- !
   ! For each atom, fills LTE populations and
   ! fill collisional matrix if it is active.
  ! -------------------------------------------------------------- !
  integer n, k
  logical :: debeye=.true.
    

  write(*,*) "Setting LTE populations..."
  do n=1,atmos%Natom

   ! fill lte populations for each atom, active or not
   CALL LTEpops(atmos%Atoms(n)%ptr_atom,debeye) !it is parralel 
   
  end do

  !!This is in case H is NLTE but has not converged pops yet, or not read from file
  !!because we need this in opacity routines for H b-f, f-f, H- b-f, H- f-f ect
  !!%NLTEpops is true only if %n /= 0 or after NLTE loop, or if pops are read from file.
  !! otherwise it is FALSE, even at the begening of NLTE.

  if (Hydrogen%active .and. .not.Hydrogen%NLTEpops) Hydrogen%n = Hydrogen%nstar 
                                            !for background Hydrogen and Hminus                                            
  !CALL calc_nHmin()

 RETURN
 END SUBROUTINE setLTEcoefficients


END MODULE lte
