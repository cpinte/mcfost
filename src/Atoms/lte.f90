! -------------------------------------------------------------- !
! -------------------------------------------------------------- !
   ! Various routines to solve for LTE populations
! -------------------------------------------------------------- !
! -------------------------------------------------------------- !


MODULE lte

 use atom_type, only : AtomType, element, atomic_orbital_sqradius
 use atmos_type, only : Npf, Tpf, ne, T, Hydrogen, Helium, ntotal_atom, Elements, icompute_atomRT, nHmin, nHtot, wght_per_H, avgweight, &
 						Atoms, Natom
 use constant
 use math
 use occupation_probability, only : wocc_n
 use io_atomic_pops, only : write_Hminus

 use constantes, only : tiny_dp, huge_dp
 use parametres, only : ldissolve, n_cells, loutput_rates
 !$ use omp_lib

 IMPLICIT NONE

 real(kind=dp), parameter :: phi_min_limit = 1d-100!tiny_dp!1d-100 !1d-50, tiny_dp
 !!real(kind=dp), parameter :: CI = 0.5*(HPLANCK**2 / (2.*PI*M_ELECTRON*KBOLTZMANN))**1.5


 CONTAINS
 
 FUNCTION getPartitionFunctionk(elem, stage, k) result (Uk)
 ! ----------------------------------------------------------------------!
  ! Interpolate the partition function of Element elem in ionisation stage
  ! stage at cell point k.
 ! ----------------------------------------------------------------------!

  type(Element) :: elem
  integer, intent(in) :: stage, k
  real(kind=dp) :: Uk, part_func(Npf), Uka(1)
  
  part_func = elem%pf(stage,:)
  Uk = exp( Interp1D(Tpf,part_func,T(k)) )
  !Uk = (10.d0)**(Uk) !29/12/2019 -> part_func is ln(U)

!->faster
  !out of bound the function return 0 not the inner (outer) bound.
!   Uka(:) = linear_1D_sorted(Npf, Tpf, part_func, 1, T(k))
!   Uk = exp(Uka(1))
!   if( T(k) < Tpf(1) ) Uk = exp(part_func(1))
!   if (T(k) > Tpf(Npf)) Uk = exp(part_func(Npf))

 RETURN
END FUNCTION getPartitionFunctionk

 FUNCTION get_logPartitionFunctionk(elem, stage, k) result (Uk)
 ! ----------------------------------------------------------------------!
  ! Interpolate the log partition function of Element elem in ionisation stage
  ! stage at cell point k.
 ! ----------------------------------------------------------------------!

  type(Element) :: elem
  integer, intent(in) :: stage, k
  real(kind=dp) :: Uk, part_func(Npf), Uka(1)
  
  part_func = elem%pf(stage,:)
  Uk = Interp1D(Tpf,part_func,T(k)) !log(Uk)
  !out of bound the function return 0 not the inner (outer) bound.
!   Uka(:) = linear_1D_sorted(Npf, Tpf, part_func, 1, T(k))
!   Uk = Uka(1)
!   if( T(k) < Tpf(1) ) Uk = part_func(1)
!   if (T(k) > Tpf(Npf)) Uk = part_func(Npf)
       
 RETURN
END FUNCTION get_logPartitionFunctionk

! function phi_T(k, gi_on_gj, dE)!Ui_on_Uj, 
! 	real(kind=dp) :: phi_T
! 	integer, intent(in) :: k
! 	real(kind=dp), intent(in) ::dE, gi_on_gj!, Ui_on_Uj
! 
! 	phi_T = gi_on_gj * CI * T(k)**(-1.5_dp) * exp(dE / ( KBOLTZMANN*T(k) ))!* Ui_on_Uj * 
! 
! return
! end function phi_T


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
  !C1 = (HPLANCK**2 / (2.*PI*M_ELECTRON*KBOLTZMANN))**1.5

  !!ionpot = ionpot * 100.*HPLANCK*CLIGHT !cm1->J
  !! -> avoiding dividing by big numbers causing overflow.
  !!maximum argument is 600, exp(600) = 3.77e260
  expo = exp(min(600d0,ionpot/(KBOLTZMANN*T(k))))
  if (ionpot/(KBOLTZMANN*T(k)) >= 600d0) expo = huge_dp

  !if exp(300) it means phi is "infinity", exp(300) == 2e130 so that's enough
  phi = C1 * (Ujl / Uj1l) * expo / (T(k)**1.5 + tiny_dp)
  if (phi < phi_min_limit) phi = 0d0 !tiny_dp ! or phi = 0d0 should be more correct ?
  								   ! but test to not divide by 0

  if (is_nan_infinity(phi)>0) then
   write(*,*) "error, phi_jl", k, phi, T(k)
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
   kT = KBOLTZMANN * T(k)
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

  kT = KBOLTZMANN * T(k)
  ni1_ni = (gi1/gi) * exp(-Ei/kT)
  RETURN
  END FUNCTION BoltzmannEq4dot20b

  !Inverse Saha
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

    !phi should be between phi_min_limit and exp(600)
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

  

  FUNCTION nH_minus(icell)
   !Saha Ionisation equation applied to H-
   integer :: icell
   real(kind=dp) :: nH_minus, nH, UH, UHm, n00
   
   UH = 2.0
   !UH = getPartitionFunctionk(elements(1)%ptr_elem, 1, icell)
   !Mostly 2
   UHm = 1.0
   !nH = nHtot(icell)
   nH = hydrogen%n(1,icell)!sum(Hydrogen%n(1:hydrogen%Nlevel-1,icell))
   nH_minus =  ne(icell) * phi_jl(icell, UHm, UH, E_ION_HMIN) * nH
      
  RETURN
  END FUNCTION nH_minus
  
!-> Check dissolve too
 SUBROUTINE LTEpops_H_and_Hmin()
  ! -------------------------------------------------------------- !
   ! Take H- in the equation
   ! computed wrt the ground state of H-
  ! -------------------------------------------------------------- !
  logical :: locupa_prob
  real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
  real(kind=dp) :: n_eff, wocc, chi0, wocc0, E00, E
  integer :: Z, dZ, k, i, m


   !Take the occupation probability only once the population of each level is known ? 
   wocc0 = 1.0 !H- not hydrogenic
   
   E00 = 1.0 * 3e-11 * EV ! Joules

   !$omp parallel &
   !$omp default(none) &
   !$omp private(k, phik, dEion,i,dZ,dE,m,sum,phiHmin,wocc, n_eff, chi0, E) &
   !$omp shared(n_cells, icompute_atomRT, Elements, Hydrogen,c2,locupa_prob,wocc0, E00) &
   !$omp shared(ne, nHmin, nHtot, T)
   !$omp do
   do k=1,n_cells

!     if (icompute_atomRT(k) /= 1) CYCLE
    if (icompute_atomRT(k) <= 0) cycle !transparent or dark


    sum = 1.0
    phik = ne(k)*phi_jl(k,1.d0,1.d0,0.d0)
    !a constant of the temperature and electron density

    do i=1, hydrogen%Nlevel
    
     E = hydrogen%E(i)
     
	 !if (i==hydrogen%Nlevel) E = E - E00 * (ne(k)/T(k))**(0.5)

	 !wrt ground level of Hmin which has E = 0 eV, and shift all other levels to E_ION_HMIN
	 !so that E(1) - E(H-) = 0.7 - 0. = ionisation of H- and E(Nlevel) = 13.6eV + 0.7eV from the
	 !ground state of Hminus.
	 !dE = E' - E(H-)_groundstate, with E' = E + EH- (because taken by default wrt to E(1))
     dE = E + E_ION_HMIN !- 2.0 * E00 * (ne(k)/T(k))**(0.5) 
     dZ = hydrogen%stage(i) + 1

     chi0 = Elements(hydrogen%periodic_table)%ptr_elem%ionpot(1+hydrogen%stage(i))


     ! --------- Boltzmann equation ------------------------------------------- !

     hydrogen%nstar(i,k)=BoltzmannEq4dot20b(k, dE, 1.0_dp, hydrogen%g(i)) 

     ! ---------- Saha equation ------------------------------------------------ !

     do m=1,dZ
      if (ne(k)>0) then
       hydrogen%nstar(i,k)=hydrogen%nstar(i,k)/phik
      else
       hydrogen%nstar(i,k) = 0d0
      endif
     end do

 
     sum = sum+hydrogen%nstar(i,k)
    end do

    nHmin(k) = hydrogen%Abund*nHtot(k)/sum * wocc0
  
    !test positivity, can be 0
    if (nHmin(k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning too small H- pop", hydrogen%ID, "nH-", nHmin(k)
       write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k)
       write(*,*) " ************************************* "
       stop
    end if
    do i=1,hydrogen%Nlevel !debug
      hydrogen%nstar(i,k) = hydrogen%nstar(i,k)*nHmin(k)

      if (hydrogen%nstar(i,k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning population of hydrogen ", hydrogen%ID, "lvl=", i, "nstar=",hydrogen%nstar(i,k), " lower than", &
        " tiny_dp."
       write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k), &
       			  " n0=", hydrogen%nstar(1,k)
       write(*,*) " ************************************* "
      stop
      end if
    end do

    if (maxval(hydrogen%nstar(:,k)) >= huge_dp) then
    write(*,*) " ************************************* "
     write(*,*) "ERROR, populations of hydrogen larger than huge_dp"
     write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k)
     write(*,*) "nstar=",hydrogen%nstar(:,k)
     write(*,*) " ************************************* "
     stop
    end if

   end do !over depth points
   !$omp end do
   !$omp  end parallel


  if (ldissolve) then
   if (loutput_rates) call write_occupation_file(52, hydrogen, 1)
   do k=1, n_cells
    if (icompute_atomRT(k)>0) then
      do i=hydrogen%Nlevel-1,1,-1  !only for bound-levels ?
       wocc = wocc_n(k, real(i,kind=dp), real(hydrogen%stage(i)),real(hydrogen%stage(i)+1),hydrogen%nstar(1,k))
		hydrogen%nstar(i,k) = hydrogen%nstar(i,k) * wocc
      enddo
    endif
   
   enddo

  endif !if locupa_prob

 !!CALL write_ltepops_file(52, hydrogen)

 RETURN
 END SUBROUTINE LTEpops_H_and_Hmin

 subroutine determine_Hminus_density()
 ! Compute the H- number density by solving a chemical 
 ! equilibrium between H-, HI and HII
 
 
 return
 end subroutine determine_Hminus_density
        
 SUBROUTINE LTEpops_H()
  ! -------------------------------------------------------------- !
   ! computed wrt the ground state of H I
   ! ########Check occupation prob.#########
   ! Beware cannot determine ground state of nH to be use in w
   ! if computed before lte pops without w are known.
   !
   ! Diss Should be okay here.
   ! In turbospectrum nHI and nHII are known by chemical equi.
   ! before computing n(i) with occupation proba, which should be
   ! similar to mine since n(i) (so implicitely nHI) are computed
   ! with no dissolution and then dissolution is taken for bound
   ! levels.
   !
  ! -------------------------------------------------------------- !
  logical :: locupa_prob
  real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
  real(kind=dp) :: n_eff, wocc, chi0, wocc0, E00, E, Egs
  integer :: Z, dZ, k, i, m
  logical :: print_diff
  real(kind=dp) :: max_nstar(hydrogen%Nlevel), min_nstar(hydrogen%Nlevel)
   
   E00 = 1.0 * 3e-11 * EV ! Joules
   Egs = hydrogen%E(1)
   
   print_diff = (maxval(hydrogen%nstar) > 0.0_dp)
   if (print_diff) then
   	do i=1,hydrogen%Nlevel
   		max_nstar(i) = maxval(hydrogen%nstar(i,:))
   		min_nstar(i) = minval(hydrogen%nstar(i,:),mask=icompute_atomrt > 0)
   	enddo
   endif

   !$omp parallel &
   !$omp default(none) &
   !$omp private(k, phik, dEion,i,dZ,dE,m,sum,phiHmin,wocc, n_eff, chi0, E) &
   !$omp shared(n_cells, icompute_atomRT, Elements, Hydrogen,c2,locupa_prob,wocc0, E00) &
   !$omp shared(ne, nHmin, nHtot, T, Egs, ldissolve, print_diff)
   !$omp do
   do k=1,n_cells

!     if (icompute_atomRT(k) /= 1) CYCLE
    if (icompute_atomRT(k) <= 0) cycle !transparent or dark

    sum = 1.0
    phik = ne(k)*phi_jl(k,1.d0,1.d0,0.d0)
    !a constant of the temperature and electron density

    do i=2, hydrogen%Nlevel
    
     E = hydrogen%E(i)

     dE = E - Egs
     dZ = hydrogen%stage(i) - hydrogen%stage(1)

     chi0 = Elements(hydrogen%periodic_table)%ptr_elem%ionpot(1+hydrogen%stage(i))

     ! --------- Boltzmann equation ------------------------------------------- !

     hydrogen%nstar(i,k)=BoltzmannEq4dot20b(k, dE, hydrogen%g(1), hydrogen%g(i))

     ! ---------- Saha equation ------------------------------------------------ !

     do m=1,dZ
      if (ne(k)>0) then
       hydrogen%nstar(i,k)=hydrogen%nstar(i,k)/phik
      else
       hydrogen%nstar(i,k) = 0d0
      endif
     end do

 
     sum = sum+hydrogen%nstar(i,k)
    end do

	hydrogen%nstar(1,k) = hydrogen%Abund*nHtot(k)/sum

    !test positivity, can be 0
    if ((hydrogen%nstar(1,k) < 0)) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning too small gs pop", hydrogen%ID, hydrogen%nstar(i,k)
       write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k)
       write(*,*) " ************************************* "
       stop
    end if
    
    do i=2,hydrogen%Nlevel !debug
      hydrogen%nstar(i,k) = hydrogen%nstar(i,k)*hydrogen%nstar(1,k)

      if (hydrogen%nstar(i,k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning population of hydrogen ", hydrogen%ID, "lvl=", i, "nstar=",hydrogen%nstar(i,k), " lower than", &
        " tiny_dp."
       write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k), &
       			  " n0=", hydrogen%nstar(1,k)
       write(*,*) " ************************************* "
      stop
      end if
    end do

    if (maxval(hydrogen%nstar(:,k)) >= huge_dp) then
    write(*,*) " ************************************* "
     write(*,*) "ERROR, populations of hydrogen larger than huge_dp"
     write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k)
     write(*,*) "nstar=",hydrogen%nstar(:,k)
     write(*,*) " ************************************* "
     stop
    end if
    
!     write(*,*) " cell ", k
!     write(*,*) T(k), nHtot(k), ne(K)
!     write(*,*) hydrogen%nstar(:,k)

   end do !over depth points
   !$omp end do
   !$omp  end parallel


  if (ldissolve) then
   if (loutput_rates) call write_occupation_file(52, hydrogen, 1)
   do k=1, n_cells
    !!sum = 0.0
    if (icompute_atomRT(k)>0) then
      do i=2, hydrogen%Nlevel-1 !only for bound-levels ?
      
      	wocc = wocc_n(k, real(i,kind=dp), real(hydrogen%stage(i)),real(hydrogen%stage(i)+1), hydrogen%nstar(1,k))
      	
      	!added to the continuum
		!!sum = sum + hydrogen%nstar(i,k) * (1.0 - wocc)
		
		!remains b-b
		hydrogen%nstar(i,k) = hydrogen%nstar(i,k) * wocc
				
      enddo
   	  wocc = wocc_n(k, real(1,kind=dp), real(hydrogen%stage(1)),real(hydrogen%stage(1)+1), hydrogen%nstar(1,k))
      hydrogen%nstar(1,k) = hydrogen%nstar(1,k) * wocc


      !!hydrogen%nstar(hydrogen%Nlevel,k) = hydrogen%nstar(hydrogen%Nlevel,k) + sum
    endif
   
   enddo

  endif !if locupa_prob
  
  if (print_diff) then
  	write(*,*) " Old Max/min nstar:", max_nstar(hydrogen%Nlevel), min_nstar(hydrogen%Nlevel)
  	write(*,*) ' New max/min nstar', maxval(hydrogen%nstar(hydrogen%Nlevel,:)), minval(hydrogen%nstar(hydrogen%nlevel,:),mask=icompute_atomrt>0)
  endif

 RETURN
 END SUBROUTINE LTEpops_H

 !new debye to do
 SUBROUTINE LTEpops(atom, debye)
  ! -------------------------------------------------------------- !
   ! Computes LTE populations of each level of the atom.
   ! Distributes populations among each level (ni) according
   ! to the total population (NI) of each stage (I),
   ! for which a level i belongs (Boltzmann's equation)
  ! -------------------------------------------------------------- !
  type (Atomtype), intent(inout) :: atom
  logical, intent(in) :: debye
  logical :: locupa_prob
  real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
  real(kind=dp) :: n_eff, wocc, chi0, wocc0
  integer :: Z, dZ, k, i, m
  integer, allocatable, dimension(:) :: ndebye

  ! debye shielding activated:
  ! It lowers the ionisation potential taking into account the charge of all levels

  ! dE_ion = -Z*(qel**2 / 4PI*EPS0) / D in J
  ! D = sqrt(EPS0/2qel**2) * sqrt(kT/ne) in m

  ! Hubeny & Mihalas (2014)
  ! "Theory of Stellar Atmospheres" p. 244-246
  ! eq. 8.86 and 8.87
 
   !Not yet
   locupa_prob = .false. !Only for He
   !Not yet. For H, use the LTEpops_H()
   !only for the bound levels of each ionisation stage ? The last level beeing the last continuum ?
   if (debye .and. locupa_prob) locupa_prob = .false.
    
   if (debye) then
    c2 = sqrt(8.0*PI/KBOLTZMANN) * &
      (SQ(Q_ELECTRON)/(4.0*PI*EPSILON_0))**1.5
   !write(*,*) "c2=",c2,  atom%Abund * nHtot(1),  atom%Abund * nHtot(2)
    allocate(ndebye(atom%Nlevel))
    ndebye(:)=0
    !reduce the ionisation potential of ion i+1 (wrt to the ground state) by nDebye(i) * (ne/T)**1/2. 
    do i=2, atom%Nlevel
     if (atom%stage(i) - atom%stage(1) > 0) nDebye(i) = atom%stage(i) + 1
    enddo
!     do i=2,atom%Nlevel !reject ground level
!      Z = atom%stage(i) !0 for neutrals, 1 for singly ionised
!      !if stage(i) - stage(1) = 0, no ndebye
!      do m=1,atom%stage(i)-atom%stage(1) !relative to ground level
!       ndebye(i) = ndebye(i) + Z
!       Z = Z + 1
!      end do
!     end do
   end if

   !$omp parallel &
   !$omp default(none) &
   !$omp private(k, phik, dEion,i,dZ,dE,m,sum,phiHmin) &
   !$omp shared(n_cells, icompute_atomRT, Elements, Hydrogen,c2,atom, debye,ndebye) &
   !$omp shared(ne, nHmin, nHtot, T)
   !$omp do
   do k=1,n_cells

!     if (icompute_atomRT(k) /= 1) CYCLE
    if (icompute_atomRT(k) <= 0) cycle !transparent or dark

    if (debye) dEion = c2*sqrt(ne(k)/T(k))

    
    sum = 1.
    phik = ne(k)*phi_jl(k,1.d0,1.d0,0.d0)
    !a constant of the temperature and electron density

    do i=2, atom%Nlevel

     dE = atom%E(i) - atom%E(1) ! relative to ground level which has 0 energy
     dZ = atom%stage(i) - atom%stage(1) !stage is 0 for neutrals



!    write(*,*) "dZ = ", dZ," ndebye*dEion=",ndebye(i)*dEion*&
!         JOULE_TO_EV, " dEion=",dEion*JOULE_TO_EV

     if (debye) dE = dE - ndebye(i)*dEion !J

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
     do m=2,dZ+1 !actually m=1,dz works since dz starts at 0 not 1 for neutrals
      if (ne(k)>0) then
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

     ! now we have, 1 + n2/n1 + ni/n1 = sum over all levels
     ! and each stage for each level,
     ! wrt the first population of the first level
     ! further, n1 + n2+n3+nij sum over all stage of each level
     ! gives ntotal.
     ! Therefore n1 = atom%ntotal/sum
!     write(*,*) "-------------------"
!     write(*,*) "Atom=",atom%ID, " A=", atom%Abund
!     write(*,*) "ntot", ntotal_atom(k,atom), " nHtot=",nHtot(k)
    atom%nstar(1,k) = atom%Abund*nHtot(k)/sum
  
    !test positivity, can be 0
    if (atom%nstar(1,k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning too small ground state population ", atom%ID, "n0=", atom%nstar(1,k)
       write(*,*) "cell=",k, atom%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k)
        !atom%nstar(1,k) = tiny_dp
       write(*,*) " ************************************* "
       stop !only if we test >= 0
    end if
!     write(*,*) "depth index=",k, "level=", 1, " n(1,k)=",atom%nstar(1,k)
    do i=2,atom%Nlevel !debug
      atom%nstar(i,k) = atom%nstar(i,k)*atom%nstar(1,k)
      !to avoid very low populations in calculations
      !if (atom%nstar(i,k) < 1d-30) atom%nstar(i,k) = 1d-30
!       write(*,*) "depth index=",k, "level=", i, " n(i,k)=",atom%nstar(i,k)

      !! check positivity. Can be 0 now, see above and phik and phi_jl
      if (atom%nstar(i,k) < 0) then !<= tiny_dp) then
       write(*,*) " ************************************* "
       write(*,*) "Warning population of atom ", atom%ID, "lvl=", i, "nstar=",atom%nstar(i,k), " lower than", &
        " tiny_dp."! Replacing by tiny_dp"
       write(*,*) "cell=",k, atom%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k), &
       			  " n0=", atom%nstar(1,k)
        !atom%nstar(i,k) = tiny_dp
       write(*,*) " ************************************* "
      stop !only if we test >= 0
      end if
    end do

    if (maxval(atom%nstar(:,k)) >= huge_dp) then
    write(*,*) " ************************************* "
     write(*,*) "ERROR, populations of atom larger than huge_dp"
     write(*,*) "cell=",k, atom%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k)
     write(*,*) "nstar=",atom%nstar(:,k)
     write(*,*) " ************************************* "
     stop
    end if

    if (atom%ID=="H") then
    	nHmin(k) = nH_minus(k)
    endif

   end do !over depth points
   !$omp end do
   !$omp  end parallel
!    write(*,*) "-------------------"

  !! Write for all grid points, no parallel
  if (locupa_prob) then
    if (loutput_rates) CALL write_occupation_file(52, atom, 1)
    do k=1, n_cells
    	if (icompute_atomRT(k)>0) then
    		do i=1, atom%Nlevel-1
    			n_eff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(atom%Nlevel) - atom%E(i)))
    			wocc = wocc_n(k, n_eff, real(atom%stage(i)),real(atom%stage(i)+1), hydrogen%nstar(1,k))
    			atom%nstar(i,k) = atom%nstar(i,k) * wocc
    		enddo
    	endif
    enddo
  endif !if locupa_prob
  !CALL write_ltepops_file(52, atom, 1)

  if (allocated(ndebye)) deallocate(ndebye)
 RETURN
 END SUBROUTINE LTEpops


 SUBROUTINE set_LTE_populations ()
  ! -------------------------------------------------------------- !
   ! For each atom, fills LTE populations and
   ! fill collisional matrix if it is active.
  ! -------------------------------------------------------------- !
  integer n, k
  logical :: debye=.false.
     
  do n=1,Natom
   ! fill lte populations for each atom, active or not, only if not read from file
   if (Atoms(n)%ptr_atom%set_ltepops) then
     if ( Atoms(n)%ptr_atom%ID == "H") then 
      allocate(nHmin(n_cells))
      

      write(*,*) "Setting LTE populations for hydrogen"
      
      !Use a routine that includes H minus also ? use at is it ? use chemequil to determine nH- ??      
      CALL LTEpops_H()

   		do k=1,n_cells
   			if (icompute_atomRT(k) > 0) nHmin(k) = nH_minus(k)
   		enddo
   		call write_Hminus()
   		!nHmin(k) = 1d10

      
     else !other atoms
   
      if (debye) then
     	 write(*,*) "Setting LTE populations for ", Atoms(n)%ptr_atom%ID, " atom (shielded)"
      else
     	 write(*,*) "Setting LTE populations for ", Atoms(n)%ptr_atom%ID," atom "
      endif
      CALL LTEpops(Atoms(n)%ptr_atom,debye) !it is parralel 
     endif
     
    !Write only if set_ltepops i.e., if not read from file.

   else !no set lte pops
   
    if (n==1) then
    	write(*,*) " -> Setting H minus density from read Hydrogen populations..."
		allocate(nHmin(n_cells))
   		do k=1,n_cells
   			if (icompute_atomRT(k) > 0) nHmin(k) = nH_minus(k)
   		enddo
   		call write_Hminus()
   	endif
   endif !set_ltepops

   !If populations not read from file (that is with NLTEpops = .fals.) we need to set n=nstar
   !as the first solution.
   !If the atom initial solution is ZERO_RADIATION field, background opacities are evaluated with n=nstar
   !and then for the loop with n = SEE(I=0)
	if (Atoms(n)%ptr_atom%active .and. .not.Atoms(n)%ptr_atom%NLTEpops) then
		if (Atoms(n)%ptr_atom%initial_solution=='OLD_POPULATIONS') then
			write(*,*) "here in LTE, set n=nstar for bckgr cont only for LTE and ZERO rad initial sol"
			stop
		endif
		!if (Atoms(n)%ptr_atom%initial_solution=="LTE_POPULATIONS") 
		Atoms(n)%ptr_atom%n(:,:) = Atoms(n)%ptr_atom%nstar(:,:)
	endif
	
	if (loutput_rates) call write_ltepops_file(52, Atoms(n)%ptr_atom)

   
  enddo

  !!This is in case H is NLTE but has not converged pops yet, or not read from file
  !!because we need this in opacity routines for H b-f, f-f, H- b-f, H- f-f ect
  !!%NLTEpops is true only if %n /= 0 or after NLTE loop, or if pops are read from file.
  !! otherwise it is FALSE, even at the begening of NLTE.
  !!if (Hydrogen%active .and. .not.Hydrogen%NLTEpops) Hydrogen%n = Hydrogen%nstar 
                                            !for background Hydrogen and Hminus                                            

 RETURN
 END SUBROUTINE set_LTE_populations

!building -> need to change output format
 SUBROUTINE write_occupation_file(unit, atom, delta_k)
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: unit
  integer, intent(in), optional :: delta_k
  integer :: dk, k, i
  real(kind=dp) :: wocc
  
  if (present(delta_k)) then
   dk = delta_k
   if (dk <= 0) then
    write(*,*) " delta_k should be > 0, setting to 1"
    dk = 1
   else if (dk >= n_cells) then
    write(*,*) "Delta_k should be < n_cells"
    dk = 1
   endif
  else
   dk = 1
  endif
  
  open(unit, file=trim(atom%ID)//"_wi.txt", status="unknown")
  do k=1, n_cells, dk
    if (icompute_atomRT(k)>0) then
      write(unit,*) k
      write(unit,*) "T=",real(T(k)), "ne(/cc)=", real(ne(k))*1e-6
      write(unit,*) "   i      w(i)     nstar(i;wi=1.0) (/cc)   nstar(i) (/cc)    g(i)"
      do i=1, atom%Nlevel-1
       wocc = wocc_n(k, real(i,kind=dp), real(atom%stage(i)),real(atom%stage(i)+1),hydrogen%nstar(1,k))
       !if (wocc < 0.95) then
        write(unit,"(1I3, 3E14.7, 1F14.7)") i, wocc, atom%nstar(i,k)/wocc * 1d-6, atom%nstar(i,k)*1d-6, atom%g(i)
       !endif
      enddo
      write(unit,"(1I3, 1E14.7, 1E14.7, 1E14.7, 1F14.7)") atom%Nlevel, wocc, atom%nstar(atom%Nlevel,k)/wocc * 1d-6, atom%nstar(atom%Nlevel,k)*1d-6, atom%g(atom%Nlevel)
    endif
   
  enddo
  close(unit)  

 RETURN
 END SUBROUTINE write_occupation_file
 
 !building
 subroutine print_pops(atom, icell, unit)
 	integer, intent(in) :: icell
 	integer, intent(in), optional :: unit
 	type (AtomType), intent(in) :: atom
 	integer :: i
 	 	
 	do i=1,atom%Nlevel
 	
 	 	if (present(unit)) then
 			write(unit,'("level #"(1I3)," n="(1E20.7E3)" m^-3"," nstar="(1E20.7E3)" m^-3")') i, atom%n(i,icell), atom%nstar(i,icell)
 		else
  			write(*,'("level #"(1I3)," n="(1E20.7E3)" m^-3"," nstar="(1E20.7E3)" m^-3")') i, atom%n(i,icell), atom%nstar(i,icell)
 		endif
 	
 	enddo
 
 return
 end subroutine print_pops
 
 SUBROUTINE write_ltepops_file(unit, atom, delta_k)
 !Depth index reversed
  type (AtomType), intent(in) :: atom
  integer, intent(in) :: unit
  integer, intent(in), optional :: delta_k
  integer :: dk, k, i, j, Nstage, Nj, ip
  real(kind=dp) :: wocc
  
  if (present(delta_k)) then
   dk = delta_k
   if (dk <= 0) then
    write(*,*) " delta_k should be > 0, setting to 1"
    dk = 1
   else if (dk >= n_cells) then
    write(*,*) "Delta_k should be < n_cells"
    dk = 1
   endif
  else
   dk = 1
  endif
  
  Nstage = atom%stage(atom%Nlevel)+1
  Nj = Nstage
  if (atom%ID=="H") Nj = Nstage + 1 !for H-
  
  open(unit, file=trim(atom%ID)//"_ltepops.txt", status="unknown")
  write(unit,*) "atom ", atom%ID, ' with ', Nj, " stages"
  write(unit,*) "xmh=",AMU*1d3, "xmy=",wght_per_H,avgweight
  write(unit,*) "xmh*xmy=",AMU*1d3*avgweight
  do k=n_cells, 1, -dk
   if (icompute_atomRT(k)>0) then
      write(unit,*) n_cells - k  + 1, "T=",T(k), "ne=", ne(k) * 1d-6
      write(unit,*) "nTotal=", atom%Abund * nHtot(k)*1d-6
      if (atom%ID=="H") then
       write(unit,*) "nH-=", nHmin(k) * 1d-6, " Debye/Z (eV)=", 1.0 * 3e-11 * (ne(k)/T(k))**0.5
      endif
      ip = 1
     stage_loop : do j=1, Nstage
        write(unit,*) j, "Q=", getPartitionFunctionk(elements(atom%periodic_table)%ptr_elem, j, k)
        level_loop : do i=ip, atom%Nlevel
         if (atom%stage(i)+1 /= j) then
          ip = i
          cycle stage_loop 
         endif
         if (atom%n(i,k) /= atom%nstar(i,k)) then
         	write(unit,*) 'ilevel=',i, 'nstar=', atom%nstar(i, k)*1d-6, 'n=', atom%n(i, k)*1d-6
         else
         	write(unit,*) 'ilevel=',i, 'nstar=', atom%nstar(i, k)*1d-6
         endif
        enddo level_loop
     enddo stage_loop
   endif
  enddo !space loop
  close(unit)  
 
 RETURN
 END SUBROUTINE write_ltepops_file

END MODULE lte

!->> building
  !building
  !to do add, negative ions, like H-
  !with debye shielding for them
!   SUBROUTINE LTEpops_new(k, atom)
!    use atom_type, only : find_ground_state_ion
!    integer, intent(in) :: k
!    type (AtomType), intent(inout) :: atom
!    real(kind=dp) :: beta, dE, C0, Uj, Uj1, Ej0
!    real(kind=dp) :: n_eff, wi
!    integer :: i, Nstage, stage, j
!    !actually, I need an array of Nstage which is lower than Nlevel
!    !but to avoid allocation for each k I use Nlevel.
!    real(kind=dp), dimension(atom%Nlevel) :: NI !array of ground state pop of each ion
!    integer, dimension(atom%Nlevel) :: j0
!    type (Element), pointer :: elem
!    logical :: methode1
!    
!    methode1 = .false. !tmp
! 
!    elem => Elements(atom%periodic_table)%ptr_elem
! 
!    beta = 1.0_dp / (T(k) * KBOLTZMANN)
!    C0 = 0.5 * (HPLANCK*HPLANCK / (2.0*pi*M_ELECTRON) * beta)**1.5
!    !(M_ELECTRON/beta * 2.0*pi/HPLANCK**2.)**(1.5)
!    !The values Nstage and j0 could be stored in atom.
!       
!    !Nstage = atom%stage(atom%Nlevel) + 1 -> yes, but in general, we could have for instance
!    ! a model atom for ion II and III ground level (like MgII model atom) and Nstage
!    !is 2, not 3 !
!    Nstage = atom%stage(atom%Nlevel) - atom%stage(1) + 1
!    !maxval(atom%stage) - minval(atom%stage) + 1
! 
!    NI(1:Nstage) = 0.0_dp !beyon Nstage not used. Nstage << Nlevel by definition
!    j0(1:Nstage) = 0
!    !find the index for each stage (0 and 1 for H) of the ground state of each stage
!    !among all levels: for H 20 it is 1 (stage==0), and 20 (stage==1)
!    j0(1:Nstage) = find_ground_state_ion(atom%Nlevel, atom%stage, Nstage)
!    !write(*,*) Nstage, (j0(j), j=1,Nstage)
! 
!    NI(Nstage) = 1.0_dp
!    
!    !Saha equation to get the population of the state-0 each ionisation stage
!    Uj1 = get_logPartitionFunctionk(elem, Nstage, k)
!    do j=Nstage-1, 1, -1
!     !!write(*,*) Nstage, j, NI(j+1), exp(Uj1)
!     !j should be atom%stage + 1 because it is an index in an array, j=1 = ionstage 0 in the array
!     Uj = get_logPartitionFunctionk(elem, j, k)
!     !Uj/Uj1 * exp(elem%ionpot(j) * beta)
!     !elem%ionpot(j) wrt the ground state of ion (j)
!     NI(j) = NI(j+1) * C0 * ne(k) * exp(elem%ionpot(j) * beta + Uj - Uj1)
!     Uj1 = Uj
!     !!write(*,*) NI(j), exp(Uj), elem%ionpot(j)/EV
!    enddo
! 
!    if (methode1) then
!     NI(Nstage) = ntotal_atom(k, atom) / sum(NI(1:Nstage))
!     NI(1:Nstage-1) = NI(1:Nstage-1) * NI(Nstage)
!    endif
! 
! 
!    !Now boltzmann equation for each level of each ion
!    atom%nstar(:,k) = 0.0_dp
!    atom%nstar(atom%Nlevel,k) = NI(Nstage) !last level not treated in detail i.e., = state-0
!    do i=1, atom%Nlevel-1
!     if (atom%ID == "H") then
!      n_eff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(j0(Nstage)) - atom%E(i)))
!      wi = wocc_n(k, sqrt(atom%g(i)/2.), real(atom%stage(i)), 1.0) !perturber = H+ ?
!      !write(*,*) sqrt(atom%g(i)/2.), "n_eff=", n_eff, "wi=", wi
!     else
!      wi = 1.0
!     endif
!     j = j0(atom%stage(i)+1)
!     Ej0 = atom%E(j)! ground state energy of each ion
!     Uj = getPartitionFunctionk(elem, j, k)
!     
!     if (methode1) then
!       atom%nstar(i,k) = wi * BoltzmannEq9dot1(k, atom%E(i)-Ej0, atom%g(i), Uj) * NI(atom%stage(i)+1)    
! 
!     else
!      atom%nstar(i,k) = wi * C0 * ne(k) * atom%g(i)/atom%g(j+1) * &
!                     exp(beta*(elem%ionpot(j)-atom%E(i)-Ej0)) * NI(j+1)
!     endif
!    enddo
!    
! 
!    !Abolute populations
!    !!write(*,*) atom%ID, " populations of cell ", k
!    !!write(*,*) (atom%nstar(i,k), i=1,atom%Nlevel)
!    if (.not.methode1) then
!     atom%nstar(atom%Nlevel,k) = ntotal_atom(k, atom) / sum(atom%nstar(:,k))
!     atom%nstar(1:atom%Nlevel-1,k) = atom%nstar(1:atom%Nlevel-1,k) * atom%nstar(atom%Nlevel,k)
!    endif
!    !!write(*,*) " total population of each level (m^-3)"
!    !!write(*,*) (atom%nstar(i,k), i=1, atom%Nlevel)
! 
!   RETURN
!   END SUBROUTINE LTEpops_new