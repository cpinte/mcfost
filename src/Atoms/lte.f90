! -------------------------------------------------------------- !
! -------------------------------------------------------------- !
   ! Various routines to solve for LTE populations
! -------------------------------------------------------------- !
! -------------------------------------------------------------- !


MODULE lte

 use atom_type, only : AtomType, element
 use atmos_type, only : atmos, Hydrogen, Helium
 use constant
 use math
 use collision, only : CollisionRate, openCollisionFile
 
 use constantes, only : tiny_dp, huge_dp
 
 !$ use omp_lib

 IMPLICIT NONE

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
  double precision :: ionpot, C1, phi
  double precision :: Ujl, Uj1l, arg_exp
  C1 = 0.5*(HPLANCK**2 / (2.*PI*M_ELECTRON*KBOLTZMANN))**1.5
  !!ionpot = ionpot * 100.*HPLANCK*CLIGHT !cm1->J
  !! -> avoiding dividing by big numbers causing overflow.
  !!maximum argument is 600, exp(600) = 3.77e260
  arg_exp = min(600d0,ionpot/(KBOLTZMANN*atmos%T(k)))
  phi = C1 * (Ujl / Uj1l) * dexp(arg_exp) / (atmos%T(k)**1.5)
  if (phi < tiny_dp) phi = tiny_dp !
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

   double precision :: Ei, gi, Ui, ni_NI, kT
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

  double precision :: Ei,ni1_ni, kT
  double precision ::  gi, gi1
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
   double precision :: NI, UI1, UI, chi, ne, phi, NI1
   phi = phi_jl(k, UI, UI1, chi) ! chi in J
   NI1 = NI/(phi*ne)

  RETURN
  END FUNCTION SahaEq

  SUBROUTINE calc_nHmin()
   integer :: k, j
   double precision :: phiHmin
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
    atmos%nHmin(k) = sum(Hydrogen%n(1:Hydrogen%Nlevel-1,k)) * atmos%ne(k)*PhiHmin!faster?
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
  double precision :: dEion, dE, sum, c2, phik, phiHmin
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
      atom%nstar(i,k)=atom%nstar(i,k)/phik
     end do

     ! by doing so we avoid the easiest case of using
     ! partition function.
     sum = sum+atom%nstar(i,k) !compute total pop = Sum_jSum_i nij
     						   !with Sum_i nij = Nj
    end do
    !Incorpore Hminus
    if (atom%ID=="H") then
     phiHmin = atmos%ne(k)*phi_jl(k, 1d0, 2d0, E_ION_HMIN)
     atmos%nHmin(k) = 1d0 * phiHmin!ground state of HI wrt to nH00
     do i = 2, atom%Nlevel-1
      atmos%nHmin(k) = atmos%nHmin(k) + atom%nstar(i,k) * phiHmin !wrt nH00
     end do
     sum = sum + atmos%nHmin(k)
    end if
    
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
    !write(*,*) "depth index=",k, "level=", 1, " n(1,k)=",atom%nstar(1,k)
    do i=2,atom%Nlevel !debug
      atom%nstar(i,k) = atom%nstar(i,k)*atom%nstar(1,k)
      !to avoid very low populations in calculations
      if (atom%nstar(i,k) < 1d-30) atom%nstar(i,k) = 1d-30
      !write(*,*) "depth index=",k, "level=", i, " n(i,k)=",atom%nstar(i,k)
    end do
    !faster ? 
!    atom%nstar(2:atom%Nlevel,k) = atom%nstar(2:atom%Nlevel,k) * atom%nstar(1,k)
    if (atom%ID=="H") atmos%nHmin(k) = atmos%nHmin(k) * atom%nstar(1,k)

    if (MAXVAL(atom%nstar(:,k)) <= tiny_dp) then !at the cell
     write(*,*) "cell=",k, atom%ID, atmos%icompute_atomRT(k), atmos%T(k), atmos%nHtot(k)
     write(*,*) atom%nstar(:,k)
     write(*,*) "Error, populations negative or lower than tiny dp for this atom at this cell point"
     write(*,*) "stop!"
     stop !beware if low T, np -> 0, check to not divide by 0 density
    end if
   end do !over depth points
   !$omp end do
   !$omp  end parallel
!    write(*,*) "-------------------"

!    if (MAXVAL(atom%nstar) <= tiny_dp .and. atmos%lcompute_atomRT(k)) then !Total
!      write(*,*) "cell=",k, atom%ID,atmos%lcompute_atomRT(k)
!      write(*,*) "Error, populations negative or lower than tiny dp for this atom on the grid!"
!      write(*,*) "Exciting..."
!      stop !beware if low T, np -> 0, check to not divide by 0 density
!    end if

  if (allocated(nDebeye)) deallocate(nDebeye)
 RETURN
 END SUBROUTINE LTEpops

!  SUBROUTINE LTEpops_element(elem)
!   -------------------------------------------------------------- !
!    Obtains the TOTAL LTE populations of each stage (NI)
!    for type(Element) elem. These populations are used
!    to compute the population of each level (ni)
!    using Boltzmann's equation, gi, gj and Ei, Ej.
!    this version is for testing,
!    electron density calculations
!    and LTE background 2-levels atoms
! 
!    At this step, we do not know which two levels atom
!    transititions are used, so only total populations
!    of each ion is computed. These populations
!    can esealy give the population of each level with
!    Botlzmann's equation
!   -------------------------------------------------------------- !
! 
!   type (Element), intent(inout) :: elem
!   integer :: j, k
!   double precision, dimension(atmos%Nspace) :: Uk, Ukp1, sum
!   double precision, dimension(atmos%Npf) :: pfel
!   double precision :: phik
! 
!   pfel = elem%pf(1,:)
! 
!   parition functions in logarithm
!   CALL bezier3_interp(atmos%Npf, atmos%Tpf, pfel, &
!                      atmos%Nspace, atmos%T, Uk)
! 
!   init
!   nj(1,:) is set to 1 -> nj(j>1,:) normalised to nj(1,:)
!   do k=1,atmos%Nspace
!    if (.not.atmos%lcompute_atomRT(k)) CYCLE
!    sum(k) = 1.
!    elem%n(1,k) = 1.
!   end do
!   if (MAXVAL(elem%n(1,:)) == 0d0) RETURN !no element elem here.
! 
!   do j=2,elem%Nstage !next ion stage
!    pfel = elem%pf(j,:)
!    CALL  bezier3_interp(atmos%Npf, atmos%Tpf, pfel, &
!        atmos%Nspace, atmos%T, Ukp1)
!    do k=1,atmos%Nspace
! 
!     phik = phi_jl(k, 10**(Uk(k)), 10**(Ukp1(k)),&
!                   elem%ionpot(j-1))
!     Saha equation here elem%n(j,k)=&
!     Sahaeq(k,elem%n(j-1,k),Ukp1,Uk,elem%ionpot(j-1),atmos%ne(k))
!     elem%n(j,k) = elem%n(j-1,k) / (phik*atmos%ne(k))
!     sum(k) = sum(k) + elem%n(j,k)
!     1.+n2/n1 + n3/n1 +...+nj/n1+...+nNstage/n1
!     Need new equations for the first stage!:
!     1 + n2/n1 + nj/n1 = sum (eq. 1)
!     n2 + n1 + nj = nElemTot = A * nHtot (eq. 2)
!     (eq. 2)/(eq. 1) ->
!     (n1 + n2 +.. +nj)/(1 + n2/n1 + ..+nj/n1) =
!     = nElemTot/sum -> remebering, 1=n1/n1!
!     -> n1 (n2 + n1 ...) / (n1 + n2 ...) = nElemTot/sum
!     n1 = nElemTot/sum = A * nHtot / sum
! 
!     set Uk=Ukp1 for next ion stage by swapping values
!     Uk(k)=Ukp1(k)
!    end do
!   end do
! 
!   because of the CYCLE statement, here sum(k) is at least 1
!   do k=1,atmos%Nspace
!    elem%n(1,k) = elem%abund*atmos%nHtot(k)/sum(k)
!   end do
!   do j=2,elem%Nstage
!    do k=1,atmos%Nspace
!     elem%n(j,k) = elem%n(j,k)*elem%n(1,k)
!    denormalise to n1
!    end do
!   end do
! 
!   RETURN
!  END SUBROUTINE LTEpops_element


 SUBROUTINE setLTEcoefficients ()
  ! -------------------------------------------------------------- !
   ! For each atom, fills LTE populations and
   ! fill collisional matrix if it is active.
  ! -------------------------------------------------------------- !
  integer n, k, j
  logical :: debeye=.true.
  double precision :: PhiHmin

  write(*,*) "Setting LTE populations..."
  do n=1,atmos%Natom

   ! fill lte populations for each atom, active or not
   CALL LTEpops(atmos%Atoms(n)%ptr_atom,debeye) !it is parralel 
!      if (atmos%Atoms(n)%ptr_atom%active) then
!        allocate(atmos%Atoms(n)%ptr_atom%C(atmos%Atoms(n)%ptr_atom%Nlevel*atmos%Atoms(n)%ptr_atom%Nlevel))
!        !Temporary
!        atmos%Atoms(n)%ptr_atom%C = 0d0
!        allocate(atmos%Atoms(n)%ptr_atom%Ckij(&
!        atmos%Nspace,atmos%Atoms(n)%ptr_atom%Nlevel*atmos%Atoms(n)%ptr_atom%Nlevel))
!        !open collision file
!        atmos%Atoms(n)%ptr_atom%Ckij = 0d0
!        CALL openCollisionFile(atmos%Atoms(n)%ptr_atom)
!      end if
  end do

  ! if (.not.chemEquil) then
  !atmos%nHmin = 0d0 !init
  
  !!This is in case H is NLTE but has not converged pops yet, or not read from file
  if (Hydrogen%active .and. .not.Hydrogen%NLTEpops) Hydrogen%n = Hydrogen%nstar 
                                            !for background Hydrogen and Hminus
  ! CALL calc_nHmin()
  !now Hminus in computed with Saha equation during LTEpops
  !should be updated with NLTE pops at the end of the NLTE loop                                         

  !!!$omp parallel &
  !!!$omp default(none) &
  !!!$omp private(k, PhiHmin,j) &
  !!!$omp shared(atmos,Hydrogen)
  !!!$omp do
!   do k=1, atmos%Nspace
!    if (.not.atmos%lcompute_atomRT(k)) CYCLE
! 
!    PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
   ! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
   ! nHmin is proportial to the total number of neutral Hydrogen
   ! which for hydrogen is given by the sum of the indivdual
   ! levels except the last one which is the protons (or H+)
   ! number density. Remember:
   ! -> Njl = Nj1l * ne * phi_jl with j the ionisation state
   ! j = -1 for Hminus (l=element=H)
!     do j=1,Hydrogen%Nlevel - 1
!      atmos%nHmin(k) = atmos%nHmin(k) + Hydrogen%n(j,k)
!     end do
!    atmos%nHmin(k) = atmos%nHmin(k) * atmos%ne(k) * Phihmin
!     atmos%nHmin(k) = sum(Hydrogen%n(1:Hydrogen%Nlevel-1,k)) * atmos%ne(k)*PhiHmin!faster?
   !! comme %n pointe vers %nstar si pureETL (car pointeurs)
   !! il n y pas besoin de preciser if %NLTEpops
   !!write(*,*) "k, H%n(1,k), Atoms(1)%n(1,k), Atoms(1)%nstar(1,k)",&
   !!   ", nHmin(k), np(k)"
!    write(*,*) "icell=",k, Hydrogen%n(1,k), atmos%Atoms(1)%ptr_atom%n(1,k), &
!    atmos%Atoms(1)%ptr_atom%nstar(1,k), atmos%nHmin(k), Hydrogen%n(Hydrogen%Nlevel,k),&
!    atmos%PAssiveAtoms(1)%ptr_atom%n(1,k)
!   end do !over depth points
  !!!$omp end do
  !!!!$omp  end parallel
  ! end if !if chemical equilibrium
   write(*,*) "Maximum/minimum H minus density in the model (m^-3):"
   write(*,*) MAXVAL(atmos%nHmin), MINVAL(atmos%nHmin,mask=atmos%icompute_atomRT>0)
! stop
 RETURN
 END SUBROUTINE setLTEcoefficients

END MODULE lte
