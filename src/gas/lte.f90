! -------------------------------------------------------------- !
! -------------------------------------------------------------- !
! Various routines to solve for LTE populations
! -------------------------------------------------------------- !
! -------------------------------------------------------------- !
module lte

   use atom_type
   use elements_type
   use grid, only : ne, T, icompute_atomRT, nHmin, nHtot
   use constantes
   use occupation_probability, only : wocc_n
   use parametres, only : ldissolve, n_cells, loutput_rates
   !$ use omp_lib

   implicit none

   integer, dimension(101) :: ndebye

   contains


   function BoltzmannEq9dot1(temp, Ei, gi, UI)
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

      real(kind=dp), intent(in) :: Ei, gi, Ui, temp
      real(kind=dp) :: kT
      real(kind=dp) :: BoltzmannEq9dot1

      kT = KB * temp
      BoltzmannEq9dot1 = (gi/UI) * exp(-Ei/kT) !ni/NI
      return
   end function BoltzmannEq9dot1

   function BoltzmannEq4dot20b(temp, Ei, gi, gi1)
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

      real(kind=dp) :: BoltzmannEq4dot20b
      real(kind=dp) :: kT
      real(kind=dp), intent(in) ::  gi, gi1, Ei, temp
      integer :: k

      kT = KB * temp
      BoltzmannEq4dot20b = (gi1/gi) * exp(-Ei/kT) !ni1/ni
      return

   end function BoltzmannEq4dot20b



   function nH_minus(icell)
   !Saha Ionisation equation applied to H-
      integer :: icell
      real(kind=dp) :: nH_minus, nH, UH, UHm, n00

      UH = 2.0
      !UH = get_pf(elements(1), 1, T(icell))
      !Mostly 2
      UHm = 1.0
      !nH = nHtot(icell)
      nH = hydrogen%n(1,icell)!sum(Hydrogen%n(1:hydrogen%Nlevel-1,icell))
      nH_minus =  ne(icell) * phi_jl(T(icell), UHm, UH, E_ION_HMIN) * nH

      return
   end function nH_minus

   subroutine determine_Hminus_density()
    ! Compute the H- number density by solving a chemical
    ! equilibrium between H-, HI and HII


      return
   end subroutine determine_Hminus_density

   subroutine LTEpops_H_loc (k)
   ! -------------------------------------------------------------- !
   ! TO DO: introduce nH-
   ! -------------------------------------------------------------- !
      integer, intent(in) :: k
      logical :: locupa_prob
      real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
      real(kind=dp) :: n_eff, wocc, chi0, wocc0, E00, E, Egs
      integer :: Z, dZ, i, m

      E00 = 1.0 * 3e-11 * electron_charge ! Joules
      Egs = hydrogen%E(1)


      sum = 1.0
      phik = ne(k)*phi_jl(T(k),1.d0,1.d0,0.d0)
      !a constant of the temperature and electron density

      do i=2, hydrogen%Nlevel

         E = hydrogen%E(i)

         dE = E - Egs
         dZ = hydrogen%stage(i) - hydrogen%stage(1)

         chi0 = Elems(hydrogen%periodic_table)%ionpot(1+hydrogen%stage(i))

         ! --------- Boltzmann equation ------------------------------------------- !

         hydrogen%nstar(i,k)=BoltzmannEq4dot20b(T(k), dE, hydrogen%g(1), hydrogen%g(i))

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


      if (ldissolve) then
         do i=2, hydrogen%Nlevel-1 !only for bound-levels ?
            wocc = wocc_n(T(k), ne(k), real(i,kind=dp), real(hydrogen%stage(i)),real(hydrogen%stage(i)+1), hydrogen%nstar(1,k))

            hydrogen%nstar(i,k) = hydrogen%nstar(i,k) * wocc
         enddo
         wocc = wocc_n(T(k), ne(k), real(1,kind=dp), real(hydrogen%stage(1)),real(hydrogen%stage(1)+1), hydrogen%nstar(1,k))
         hydrogen%nstar(1,k) = hydrogen%nstar(1,k) * wocc
      endif


      return
   end subroutine LTEpops_H_loc

   subroutine LTEpops_atom_loc(k,atom,debye)
   ! -------------------------------------------------------------- !
   ! Computes LTE populations of each level of the atom.
   ! Distributes populations among each level (ni) according
   ! to the total population (NI) of each stage (I),
   ! for which a level i belongs (Boltzmann's equation)
   ! -------------------------------------------------------------- !
      integer, intent(in) :: k
      type (Atomtype), intent(inout) :: atom
      logical, intent(in) :: debye
      logical :: locupa_prob, print_diff
      real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
      real(kind=dp) :: n_eff, wocc, chi0, wocc0
      integer :: Z, dZ, i, m

      ! debye shielding activated:
      ! It lowers the ionisation potential taking into account the charge of all levels

      ! dE_ion = -Z*(qel**2 / 4PI*EPS0) / D in J
      ! D = sqrt(EPS0/2qel**2) * sqrt(kT/ne) in m

      ! Hubeny & Mihalas (2014)
      ! "Theory of Stellar Atmospheres" p. 244-246
      ! eq. 8.86 and 8.87

      locupa_prob = .false.
      if (debye .and. locupa_prob) locupa_prob = .false.

      if (debye) then
         c2 = sqrt(8.0*PI/KB) * (electron_charge**2/(4.0*PI*EPSILON_0))**1.5
         !write(*,*) "c2=",c2,  atom%Abund * nHtot(1),  atom%Abund * nHtot(2)
         ndebye(1:atom%Nlevel)=0
         !reduce the ionisation potential of ion i+1 (wrt to the ground state) by nDebye(i) * (ne/T)**1/2.
         do i=2, atom%Nlevel
            if (atom%stage(i) - atom%stage(1) > 0) nDebye(i) = atom%stage(i) + 1
         enddo
      end if

      if (debye) dEion = c2*sqrt(ne(k)/T(k))

      sum = 1.0
      phik = ne(k)*phi_jl(T(k),1.d0,1.d0,0.d0)
      !a constant of the temperature and electron density

      do i=2, atom%Nlevel

         dE = atom%E(i) - atom%E(1) ! relative to ground level which has 0 energy
         dZ = atom%stage(i) - atom%stage(1) !stage is 0 for neutrals


         if (debye) dE = dE - ndebye(i)*dEion !J

          ! --------- Boltzmann equation ------------------------------------------- !
          ! relative to ni=n1, the populations
          ! of the ground state in stage stage(1).

         atom%nstar(i,k)=BoltzmannEq4dot20b(T(k), dE, atom%g(1), atom%g(i))

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
      do i=2,atom%Nlevel !debug
         atom%nstar(i,k) = atom%nstar(i,k)*atom%nstar(1,k)
         if (atom%nstar(i,k) < 0) then !<= tiny_dp) then
            write(*,*) " ************************************* "
            write(*,*) "Warning population of atom ", atom%ID, "lvl=", i, "nstar=",atom%nstar(i,k), " lower than", &
               " tiny_dp."! Replacing by tiny_dp"
            write(*,*) "cell=",k, atom%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k), &
               " n0=", atom%nstar(1,k)
            write(*,*) " ************************************* "
            stop
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

      !! Write for all grid points, no parallel
      if (locupa_prob) then
         do i=1, atom%Nlevel-1
            n_eff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(atom%Nlevel) - atom%E(i)))
            wocc = wocc_n(T(k), ne(k), n_eff, real(atom%stage(i)),real(atom%stage(i)+1), hydrogen%nstar(1,k))
            atom%nstar(i,k) = atom%nstar(i,k) * wocc
         enddo
      endif

      return
   end subroutine LTEpops_atom_loc

   subroutine ltepops_atoms
      ! -------------------------------------------------------------- !
      ! Computes LTE populations for each atom
      ! TO DO: print max/min for all atoms
      ! -------------------------------------------------------------- !
         integer n, k
         logical :: debye=.false.

         !$omp parallel &
         !$omp default(none) &
         !$omp private(k,n) &
         !$omp shared(n_cells, icompute_atomRT, T) &
         !$omp shared(nHmin, n_atoms, atoms, debye, hydrogen)
         !$omp do schedule(dynamic)
         do k=1, n_cells

            if (icompute_atomRT(k) <= 0) cycle !transparent or dark

            n = 1
            if (hydrogen%set_ltepops) then
               call LTEpops_H_loc(k)

               if (.not.hydrogen%active)  then
                  nHmin(k) = nH_minus(k)
               endif
            else !.not. set_ltepops, populations read from file, just compute H- if n==1 (H)
               nHmin(k) = nH_minus(k)
            endif !set_ltepops

            if (hydrogen%active .and. .not.hydrogen%NLTEpops) then
               if (hydrogen%initial==1) then
                  write(*,*) "here in LTE, set n=nstar for bckgr cont only for LTE and ZERO rad initial sol"
                  stop
               endif
               hydrogen%n(:,k) = hydrogen%nstar(:,k)
               nHmin(k) = nH_minus(k)
            endif

            !hydrogen always first and present
            do n=2, n_atoms
               if (Atoms(n)%p%set_ltepops) then
                  call LTEpops_atom_loc(k,Atoms(n)%p,debye)
               endif !set_ltepops
               if (Atoms(n)%p%active .and. .not.Atoms(n)%p%NLTEpops) then
                  if (Atoms(n)%p%initial==1) then
                     write(*,*) "here in LTE, set n=nstar for bckgr cont only for LTE and ZERO rad initial sol"
                     stop
                  endif
                  Atoms(n)%p%n(:,k) = Atoms(n)%p%nstar(:,k)
               endif
            enddo

         enddo
         !$omp end do
         !$omp  end parallel

         !call write_Hminus()
      return
   end subroutine ltepops_atoms

   subroutine LTEpops_H()
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
      ! TO DO: introduce nH-
      ! -------------------------------------------------------------- !
         logical :: locupa_prob
         real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
         real(kind=dp) :: n_eff, wocc, chi0, wocc0, E00, E, Egs
         integer :: Z, dZ, k, i, m
         logical :: print_diff
         real(kind=dp) :: max_nstar(hydrogen%Nlevel), min_nstar(hydrogen%Nlevel)

         E00 = 1.0 * 3e-11 * electron_charge ! Joules
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
         !$omp shared(n_cells, icompute_atomRT, Elems, Hydrogen,c2,locupa_prob,wocc0, E00) &
         !$omp shared(ne, nHmin, nHtot, T, Egs, ldissolve, print_diff)
         !$omp do schedule(dynamic)
         do k=1,n_cells

            if (icompute_atomRT(k) <= 0) cycle !transparent or dark

            sum = 1.0
            phik = ne(k)*phi_jl(T(k),1.d0,1.d0,0.d0)
            !a constant of the temperature and electron density

            do i=2, hydrogen%Nlevel

               E = hydrogen%E(i)

               dE = E - Egs
               dZ = hydrogen%stage(i) - hydrogen%stage(1)

               chi0 = Elems(hydrogen%periodic_table)%ionpot(1+hydrogen%stage(i))

             ! --------- Boltzmann equation ------------------------------------------- !

               hydrogen%nstar(i,k)=BoltzmannEq4dot20b(T(k), dE, hydrogen%g(1), hydrogen%g(i))

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
                  write(*,*) "Warning population of hydrogen ", hydrogen%ID, "lvl=", i, "nstar=",hydrogen%nstar(i,k), &
                       " lower than tiny_dp."
                  write(*,*) "cell=",k, hydrogen%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), &
                       "ne=",ne(k), " n0=", hydrogen%nstar(1,k)
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
               !!sum = 0.0
               if (icompute_atomRT(k)>0) then
                  do i=2, hydrogen%Nlevel-1 !only for bound-levels ?
                     wocc = wocc_n(T(k), ne(k), real(i,kind=dp), real(hydrogen%stage(i)),real(hydrogen%stage(i)+1), &
                          hydrogen%nstar(1,k))

                     ! added to the continuum
                     !!sum = sum + hydrogen%nstar(i,k) * (1.0 - wocc)

                     !remains b-b
                     hydrogen%nstar(i,k) = hydrogen%nstar(i,k) * wocc

                  enddo
                  wocc = wocc_n(T(k), ne(k), real(1,kind=dp), real(hydrogen%stage(1)),real(hydrogen%stage(1)+1), &
                       hydrogen%nstar(1,k))
                  hydrogen%nstar(1,k) = hydrogen%nstar(1,k) * wocc
                  !hydrogen%nstar(hydrogen%Nlevel,k) = hydrogen%nstar(hydrogen%Nlevel,k) + sum
               endif

            enddo

         endif !if locupa_prob

         if (print_diff) then
            write(*,'("max(rho)="(1ES20.7E3)" m^-3"," min(rho)="(1ES20.7E3)" m^-3")'), &
               maxval(hydrogen%Abund*nHtot), minval(hydrogen%Abund*nHtot,mask=nHtot > 0)
            write(*,'("Old max(nstar)="(1ES20.7E3)" m^-3"," min(nstar)="(1ES20.7E3)" m^-3")'), &
                maxval(max_nstar(:)), minval(min_nstar(:))
            write(*,'("New max(nstar)="(1ES20.7E3)" m^-3"," min(nstar)="(1ES20.7E3)" m^-3")'), &
               maxval(maxval(hydrogen%nstar(:,:),dim=2)),&
               minval(minval(hydrogen%nstar(:,:),dim=2,mask=hydrogen%nstar(:,:)>0))
         endif

         return
   end subroutine LTEpops_H

   subroutine LTEpops_atom(atom, debye)
   ! -------------------------------------------------------------- !
   ! Computes LTE populations of each level of the atom.
   ! Distributes populations among each level (ni) according
   ! to the total population (NI) of each stage (I),
   ! for which a level i belongs (Boltzmann's equation)
   ! -------------------------------------------------------------- !
      type (Atomtype), intent(inout) :: atom
      logical, intent(in) :: debye
      logical :: locupa_prob, print_diff
      real(kind=dp) :: dEion, dE, sum, c2, phik, phiHmin
      real(kind=dp) :: n_eff, wocc, chi0, wocc0
      integer :: Z, dZ, k, i, m
      real(kind=dp), dimension(atom%Nlevel) :: max_nstar, min_nstar !old_values

      ! debye shielding activated:
      ! It lowers the ionisation potential taking into account the charge of all levels

      ! dE_ion = -Z*(qel**2 / 4PI*EPS0) / D in J
      ! D = sqrt(EPS0/2qel**2) * sqrt(kT/ne) in m

      ! Hubeny & Mihalas (2014)
      ! "Theory of Stellar Atmospheres" p. 244-246
      ! eq. 8.86 and 8.87

      print_diff = (maxval(atom%nstar) > 0.0_dp)
      if (print_diff) then
         do i=1,atom%Nlevel
            max_nstar(i) = maxval(atom%nstar(i,:))
            min_nstar(i) = minval(atom%nstar(i,:),mask=icompute_atomrt > 0)
         enddo
      endif

      locupa_prob = .false. !Only for He
      !Not yet. For H, use the LTEpops_H()
      !only for the bound levels of each ionisation stage ? The last level beeing the last continuum ?
      if (debye .and. locupa_prob) locupa_prob = .false.

      if (debye) then
         c2 = sqrt(8.0*PI/KB) * (electron_charge**2/(4.0*PI*EPSILON_0))**1.5
         !write(*,*) "c2=",c2,  atom%Abund * nHtot(1),  atom%Abund * nHtot(2)
         ndebye(1:atom%Nlevel)=0
         !reduce the ionisation potential of ion i+1 (wrt to the ground state) by nDebye(i) * (ne/T)**1/2.
         do i=2, atom%Nlevel
            if (atom%stage(i) - atom%stage(1) > 0) nDebye(i) = atom%stage(i) + 1
         enddo
      end if

      !$omp parallel &
      !$omp default(none) &
      !$omp private(k, phik, dEion,i,dZ,dE,m,sum,phiHmin) &
      !$omp shared(n_cells, icompute_atomRT, Elems, Hydrogen,c2,atom, debye,ndebye) &
      !$omp shared(ne, nHmin, nHtot, T)
      !$omp do schedule(dynamic)
      do k=1,n_cells

         if (icompute_atomRT(k) <= 0) cycle !transparent or dark

         if (debye) dEion = c2*sqrt(ne(k)/T(k))

         sum = 1.0
         phik = ne(k)*phi_jl(T(k),1.d0,1.d0,0.d0)
         !a constant of the temperature and electron density

         do i=2, atom%Nlevel

            dE = atom%E(i) - atom%E(1) ! relative to ground level which has 0 energy
            dZ = atom%stage(i) - atom%stage(1) !stage is 0 for neutrals


            if (debye) dE = dE - ndebye(i)*dEion !J

          ! --------- Boltzmann equation ------------------------------------------- !
          ! relative to ni=n1, the populations
          ! of the ground state in stage stage(1).

            atom%nstar(i,k)=BoltzmannEq4dot20b(T(k), dE, atom%g(1), atom%g(i))

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
         do i=2,atom%Nlevel !debug
            atom%nstar(i,k) = atom%nstar(i,k)*atom%nstar(1,k)
            if (atom%nstar(i,k) < 0) then !<= tiny_dp) then
               write(*,*) " ************************************* "
               write(*,*) "Warning population of atom ", atom%ID, "lvl=", i, "nstar=",atom%nstar(i,k), " lower than", &
                  " tiny_dp."! Replacing by tiny_dp"
               write(*,*) "cell=",k, atom%ID, "dark?=",icompute_atomRT(k), "T=",T(k), "nH=",nHtot(k), "ne=",ne(k), &
                  " n0=", atom%nstar(1,k)
               write(*,*) " ************************************* "
               stop
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
                  wocc = wocc_n(T(k), ne(k), n_eff, real(atom%stage(i)),real(atom%stage(i)+1), hydrogen%nstar(1,k))
                  atom%nstar(i,k) = atom%nstar(i,k) * wocc
               enddo
            endif
         enddo
      endif !if locupa_prob

      if (print_diff) then
         write(*,'("max(rho)="(1ES20.7E3)" m^-3"," min(rho)="(1ES20.7E3)" m^-3")'), &
            maxval(atom%Abund*nHtot), minval(atom%Abund*nHtot,mask=nHtot > 0)
         write(*,'("Old max(nstar)="(1ES20.7E3)" m^-3"," min(nstar)="(1ES20.7E3)" m^-3")'), &
            maxval(max_nstar(:)), minval(min_nstar(:))
         write(*,'("New max(nstar)="(1ES20.7E3)" m^-3"," min(nstar)="(1ES20.7E3)" m^-3")'), &
            maxval(maxval(atom%nstar(:,:),dim=2)),&
            minval(minval(atom%nstar(:,:),dim=2,mask=atom%nstar(:,:)>0))
      endif

      return
   end subroutine LTEpops_atom

   subroutine ltepops_atoms_1 ()
   ! -------------------------------------------------------------- !
   ! Computes LTE populations for each atom
   ! -------------------------------------------------------------- !
      integer n, k
      logical :: debye=.false.

      do n=1,N_atoms
      ! fill lte populations for each atom, active or not, only if not read from file
         if (Atoms(n)%p%set_ltepops) then
            if ( Atoms(n)%p%ID == "H") then
               write(*,*) "Setting LTE populations for hydrogen"

               call LTEpops_H()

               if (.not.hydrogen%active)  then
                  write(*,*) " -> setting H- density (H is passive)"
                  do k=1,n_cells
                     if (icompute_atomRT(k) > 0) nHmin(k) = nH_minus(k)
                  enddo
                  !nHmin(k) = 1d10
               endif
            else !other atoms
               if (debye) then
                  write(*,*) "Setting LTE populations for ", Atoms(n)%p%ID, " atom (shielded)"
               else
                  write(*,*) "Setting LTE populations for ", Atoms(n)%p%ID," atom "
               endif
               call LTEpops_atom(Atoms(n)%p,debye)
          endif

          !Write only if set_ltepops i.e., if not read from file.

         else !.not. set_ltepops, populations read from file, just compute H- if n==1 (H)

            if (n==1) then
               write(*,*) " -> Setting H- density from read Hydrogen populations..."
               do k=1,n_cells
                  if (icompute_atomRT(k) > 0) nHmin(k) = nH_minus(k)
               enddo
            endif

         endif !set_ltepops

         !If populations not read from file (NLTEpops = .fals.) we need to set n=nstar
         !as the first solution.
         !If the atom initial solution is ZERO_RADIATION field, background opacities are evaluated with n=nstar
         !and then for the loop with n = SEE(I=0)
         if (Atoms(n)%p%active .and. .not.Atoms(n)%p%NLTEpops) then
            if (Atoms(n)%p%initial==1) then
               write(*,*) "here in LTE, set n=nstar for bckgr cont only for LTE and ZERO rad initial sol"
               stop
            endif
            Atoms(n)%p%n(:,:) = Atoms(n)%p%nstar(:,:)
            if (n==1) then
               !Here need to compute H- from initial solution of n
               write(*,*) " -> evaluating H- density with LTE H pops (H is active)"
               do k=1,n_cells
                  if (icompute_atomRT(k) > 0) nHmin(k) = nH_minus(k)
               enddo
            endif
         endif

         if (loutput_rates) call write_ltepops_file(52, Atoms(n)%p)


      enddo
      !call write_Hminus()

      return
   end subroutine ltepops_atoms_1

   !building -> need to change output format
   subroutine write_occupation_file(unit, atom, delta_k)
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
               wocc = wocc_n(T(k), ne(k), real(i,kind=dp), real(atom%stage(i)),real(atom%stage(i)+1),hydrogen%nstar(1,k))
             !if (wocc < 0.95) then
               write(unit,"(1I3, 3E14.7, 1F14.7)") i, wocc, atom%nstar(i,k)/wocc * 1d-6, atom%nstar(i,k)*1d-6, atom%g(i)
             !endif
            enddo
            write(unit,"(1I3, 1E14.7, 1E14.7, 1E14.7, 1F14.7)") atom%Nlevel, wocc, atom%nstar(atom%Nlevel,k)/wocc * 1d-6, &
               atom%nstar(atom%Nlevel,k)*1d-6, atom%g(atom%Nlevel)
         endif

      enddo
      close(unit)

      return
   end subroutine write_occupation_file

   subroutine print_pops(atom, icell, unit)
      integer, intent(in) :: icell
      integer, intent(in), optional :: unit
      type (AtomType), intent(in) :: atom
      integer :: i

      if (present(unit)) then
         write(unit,*) " ---------------------------------------------------- "
      else
         write(*,*) " ---------------------------------------------------- "
      endif
      do i=1,atom%Nlevel

         if (present(unit)) then
            write(unit,'("level #"(1I3)," n="(1ES20.7E3)" m^-3"," nstar="(1ES20.7E3)" m^-3")') i, atom%n(i,icell), &
                 atom%nstar(i,icell)
         else
            write(*,'("level #"(1I3)," n="(1ES20.7E3)" m^-3"," nstar="(1ES20.7E3)" m^-3")') i, atom%n(i,icell), &
                 atom%nstar(i,icell)
         endif

      enddo
      if (present(unit)) then
         write(unit,*) " ---------------------------------------------------- "
      else
         write(*,*) " ---------------------------------------------------- "
      endif

      return
   end subroutine print_pops

  subroutine write_ltepops_file(unit, atom, delta_k)
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
   write(unit,*) "xmh=",amu, "xmy=",wght_per_H,avgweight
   write(unit,*) "xmh*xmy=",amu*avgweight
   do k=n_cells, 1, -dk
      if (icompute_atomRT(k)>0) then
         write(unit,*) n_cells - k  + 1, "T=",T(k), "ne=", ne(k) * 1d-6
         write(unit,*) "nTotal=", atom%Abund * nHtot(k)*1d-6
         if (atom%ID=="H") then
            write(unit,*) "nH-=", nHmin(k) * 1d-6, " Debye/Z (eV)=", 1.0 * 3e-11 * (ne(k)/T(k))**0.5
         endif
         ip = 1
         stage_loop : do j=1, Nstage
            write(unit,*) j, "Q=", get_pf(elems(atom%periodic_table), j, T(k))
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

   return
   end subroutine write_ltepops_file

end module lte
