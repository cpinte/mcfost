module see

    use atom_type
    use elements_type
    use grid
    use parametres
    use elecdensity, only  : ionisation_frac_lte
    use gas_contopac, only        : H_bf_Xsection
    use wavelengths, only         :  n_lambda
    use wavelengths_gas, only     : Nlambda_max_line, Nlambda_max_trans, Nlambda_max_cont, n_lambda_cont, tab_lambda_cont, tab_lambda_nm
    use utils, only               : gaussslv, is_nan_infinity_vector, linear_1D_sorted, is_nan_infinity_matrix
    use opacity_atom, only : phi_loc, psi, chi_up, chi_down, uji_down, Itot, eta_atoms, chi_tot
    use messages, only : warning, error
    use collision_atom, only : collision_rates_atom_loc, collision_rates_hydrogen_loc

    implicit none

    real(kind=dp), parameter :: frac_limit_pops = 1d-50
    real(kind=dp), allocatable :: n_new(:,:,:), lcell_converged(:), ne_new(:), ngpop(:,:,:,:)
    real(kind=dp), allocatable :: tab_Uji_cont(:,:,:), tab_Vij_cont(:,:,:)



    contains

    subroutine alloc_nlte_var(n_rayons_max, n_rayons_sub, lsubiteration)
        integer, intent(in) :: n_rayons_max, n_rayons_sub
        integer :: Nmaxlevel,Nmaxline,Nmaxcont
        logical, intent(in) :: lsubiteration
        integer :: alloc_status, mem_alloc_local
        integer :: kr, n, icell, l, i, j
        real(kind=dp) :: anu1, anu2, e1, e2, b1, b2
        type(AtomType), pointer :: atom

        mem_alloc_local = 0
        NmaxLevel = 0
        NmaxLine = 0
        Nmaxcont = 0
        !neglect size of rates
        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            Nmaxlevel = max(Nmaxlevel,atom%Nlevel)
            Nmaxcont = max(NmaxCont,atom%Ncont)
            NmaxLine = max(NmaxLine,atom%Nline)
            allocate(atom%Gamma(atom%Nlevel,atom%Nlevel,nb_proc))
            do kr=1, atom%Nline
                allocate(atom%lines(kr)%Rij(nb_proc),atom%lines(kr)%Rji(nb_proc))
                allocate(atom%lines(kr)%Cij(nb_proc),atom%lines(kr)%Cji(nb_proc))
            enddo
            do kr=1, atom%Ncont
                allocate(atom%continua(kr)%Rij(nb_proc),atom%continua(kr)%Rji(nb_proc))
                allocate(atom%continua(kr)%Cij(nb_proc),atom%continua(kr)%Cji(nb_proc))
            enddo
            atom => null()
        enddo

        !-> allocate space for non-LTE only quantities !
        allocate(ne_new(n_cells), stat=alloc_status)
        ne_new(:) = ne(:)
        if (alloc_status > 0) call error("Allocation error ne_new")
        write(*,*) " size ne_new:", sizeof(ne_new) / 1024./1024.," MB"
        allocate(n_new(NactiveAtoms,Nmaxlevel,n_cells),stat=alloc_status)
        if (alloc_status > 0) call error("Allocation error n_new")
        !only count n_new has atom%n is global in the code and "normal"
        write(*,*) " size n_new:", sizeof(n_new) / 1024./1024.," MB"
        n_new(:,:,:) = 0.0_dp
        allocate(psi(n_lambda, n_rayons_max, nb_proc), stat=alloc_status); psi = 0.0_dp
        write(*,*) " size psi:", sizeof(psi) / 1024./1024.," MB"
        if (alloc_Status > 0) call error("Allocation error psi in nlte loop")
        allocate(eta_atoms(n_lambda,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error eta_atoms")
        write(*,*) " size eta_atoms:", sizeof(eta_atoms) / 1024./1024.," MB"
        allocate(Uji_down(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error Uji_down")
        write(*,*) " size cross-coupling:", 3 * sizeof(Uji_down) / 1024./1024./1024.," GB"
        allocate(chi_down(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error chi_down")
        allocate(chi_up(n_lambda,Nmaxlevel,NactiveAtoms,nb_proc),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error chi_up")
        allocate(lcell_converged(n_cells),stat=alloc_status)
        if (alloc_Status > 0) call error("Allocation error lcell_converged")
        write(*,*) " size lcell_converged:", sizeof(lcell_converged) / 1024./1024./1024.," GB"
    
        if (lNg_acceleration) write(*,*) " size ngpop:", sizeof(ngpop)/1024./1024., " MB"

        if (lsubiteration) then
            allocate(Itot(n_lambda,n_rayons_sub,nb_proc),stat=alloc_status)
            if (alloc_Status > 0) call error("Allocation error Itot")
            allocate(phi_loc(Nlambda_max_line,Nmaxline,NactiveAtoms,n_rayons_sub,nb_proc),stat=alloc_status)
            if (alloc_Status > 0) call error("Allocation error phi_loc")
            allocate(chi_tot(n_lambda))
        else
            allocate(Itot(n_lambda,n_rayons_max,nb_proc),stat=alloc_status)
            if (alloc_Status > 0) call error("Allocation error Itot")
            allocate(phi_loc(Nlambda_max_line,Nmaxline,NactiveAtoms,1,nb_proc),stat=alloc_status)
            if (alloc_Status > 0) call error("Allocation error phi_loc")
        endif
        write(*,*) " size Itot:", sizeof(Itot) / 1024./1024./1024.," GB"
        write(*,*) " size phi_loc:", sizeof(phi_loc) / 1024./1024./1024.," GB"

        mem_alloc_local = mem_alloc_local + sizeof(psi) + sizeof(eta_atoms) + sizeof(Itot) + sizeof(phi_loc) + &
            sizeof(uji_down)+sizeof(chi_down)+sizeof(chi_up) + sizeof(ne_new) + sizeof(lcell_converged) + sizeof(n_new)


        allocate(tab_Uji_cont(NmaxCont,NactiveAtoms,n_cells))
        allocate(tab_Vij_cont(Nlambda_max_cont,NmaxCont,NactiveAtoms))
        write(*,*) " size tab_Uji_cont:", sizeof(tab_Uji_cont) / 1024./1024./1024.," GB"
        write(*,*) " size tab_Vij_cont:", sizeof(tab_Vij_cont) / 1024./1024./1024.," GB"
        !integrate in frequency Uji which is fourpi/h * anu * 2hnu3/c2 * exp(-hnu/kT)
        !for each cont get the cross-sections
        mem_alloc_local = mem_alloc_local + sizeof(tab_Uji_cont) + sizeof(tab_Vij_cont)
        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            do kr=1,atom%Ncont
                i = atom%Continua(kr)%i
                j = atom%continua(kr)%j
                if (atom%continua(kr)%hydrogenic) then
                    tab_Vij_cont(1:atom%continua(kr)%Nlambda,kr,n) = H_bf_Xsection(atom%continua(kr), tab_lambda_nm(atom%continua(kr)%Nb:atom%continua(kr)%Nr))
                else
                    tab_Vij_cont(1:atom%continua(kr)%Nlambda,kr,n) = linear_1D_sorted(size(atom%continua(kr)%alpha_file),&
                        atom%continua(kr)%lambda_file,atom%continua(kr)%alpha_file,atom%continua(kr)%Nlambda,tab_lambda_nm(atom%continua(kr)%Nb:atom%continua(kr)%Nr))
                endif
                !loop for all cells here can be long
                !could be para
                do icell=1, n_cells
                    if (icompute_atomRT(icell) > 0) then
                        tab_Uji_cont(kr,n,icell) = 0.0_dp !frequncy and angle integrated!
                        do l=2,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1
                            anu1 = tab_Vij_cont(l,kr,n)
                            anu2 = tab_Vij_cont(l-1,kr,n)

                            e1 = exp(-hc_k/T(icell)/tab_lambda_nm(atom%continua(kr)%Nb+l-1))
                            e2 = exp(-hc_k/T(icell)/tab_lambda_nm(atom%continua(kr)%Nb+l-2))

                            b1 = twohc/tab_lambda_nm(atom%continua(kr)%Nb+l-1)**4
                            b2 = twohc/tab_lambda_nm(atom%continua(kr)%Nb+l-2)**4

                            tab_Uji_cont(kr,n,icell) = tab_Uji_cont(kr,n,icell) + fourpi_h * 0.5*(anu1*e1*b1 + anu2*e2*b2) * &
                                (tab_lambda_nm(atom%continua(kr)%Nb+l-1) - tab_lambda_nm(atom%continua(kr)%Nb+l-2))
  
                        enddo
                    endif
                enddo 
            enddo
        enddo
        atom => null()


        write(*,'("Total memory allocated in NLTEloop:"(1F14.3)" GB")') mem_alloc_local / 1024./1024./1024.

        return 
    end subroutine alloc_nlte_var

    subroutine dealloc_nlte_var()
        integer :: n, kr
        type (AtomType), pointer :: atom

        if (allocated(ngpop)) deallocate(ngpop)
        deallocate(ne_new,psi,eta_atoms,chi_up,chi_down,uji_down,lcell_converged)
        deallocate(itot,phi_loc,n_new)    

        deallocate(tab_Uji_cont, tab_Vij_cont) 

        do n=1, NactiveAtoms
            atom => ActiveAtoms(n)%p
            deallocate(atom%Gamma)
            do kr=1, atom%Nline
                deallocate(atom%lines(kr)%Rij,atom%lines(kr)%Rji)
                deallocate(atom%lines(kr)%Cij,atom%lines(kr)%Cji)
            enddo
            do kr=1, atom%Ncont
                deallocate(atom%continua(kr)%Rij,atom%continua(kr)%Rji)
                deallocate(atom%continua(kr)%Cij,atom%continua(kr)%Cji)
            enddo
            atom => null()
        enddo  

        if (allocated(chi_tot)) deallocate(chi_tot)
 
        return 
    end subroutine dealloc_nlte_var


    subroutine see_atom(id,icell,atom,dM)
    ! --------------------------------------------------------------------!
    ! For atom atom solves for the Statistical Equilibrium Equations (SEE)
    !
    !
    ! For numerical stability, the row with the largest populations is
    ! removed (i.e., it is replaced by the mass conservation law).
    !
    ! see Hubeny & Mihalas Stellar atmospheres, p. 448-450
    !	eqs 14.6 to 14.8c
    !
    !
    ! --------------------------------------------------------------------!

        integer, intent(in) :: icell, id
        type(AtomType), intent(inout) :: atom
        integer :: lp, imaxpop, l, Nsmall_pops, Nneg_or_null_pops
        integer :: level_index(atom%Nlevel)
        real(kind=dp), intent(out) :: dM
        real(kind=dp), dimension(atom%Nlevel) :: ndag, delta
        real(kind=dp) :: ntotal, Gamma_dag(atom%Nlevel,atom%Nlevel) !debug
        real(kind=dp) :: dn_n, wjac

    !if molecules (or H-, ntotal must be different from A*nHtot)
        ntotal = atom%Abund*nHtot(icell) !in atomic form here. Could be stored on mem.

        ndag(:) = atom%n(:,icell) / ntotal
        imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
    !imaxpop = atom%Nlevel

        !see Hubeny & Mihalas 2014 eq. 14.8a
        !calc the term in delta(l,l') (diagonal elements of rate matrix)
        !delta(l,l') = Gamma(l,l) = sum_l' R_ll' + C_ll' = sum_l' -Gamma(l',l)
        ! do l = 1, atom%Nlevel
        !     atom%Gamma(l,l,id) = 0.0_dp
        !     !Gamma(j,i) = -(Rij + Cij); Gamma(i,j) = -(Rji + Cji)
        !     !Gamma(i,i) = Rij + Cij = -sum(Gamma(j,i))
        !     atom%Gamma(l,l,id) = sum(-atom%Gamma(:,l,id)) !positive
        ! end do
        !check that diagonal is > 0 and off-diagonals < 0
        !call check_elements(id, icell, atom%Nlevel, atom%Gamma(:,:,id), "in see_atom: "//atom%ID)

        atom%n(:,icell) = 0d0
        atom%n(imaxpop,icell) = 1.0_dp

        Gamma_dag = atom%Gamma(:,:,id) !with diagonal
        atom%Gamma(imaxpop,:,id) = 1.0_dp

    !solve for residual
        wjac = 1.0_dp
        delta = atom%n(:,icell) - matmul(atom%Gamma(:,:,id), ndag)

        call GaussSlv(atom%Gamma(:,:,id), delta(:), atom%Nlevel)
   	    atom%n(:,icell) = ndag(:) + wjac * delta(:)
        ! call gaussslv(atom%Gamma(:,:,id), atom%n(:,icell), atom%Nlevel)

        if ((maxval(atom%n(:,icell)) < 0.0)) then
            write(*,*) atom%ID, id, icell
            write(*,'("nstar: "*(ES14.5E3))') (ndag(l),l=1,atom%Nlevel) !*ntotal
            write(*,'("n: "*(ES14.5E3))') (atom%n(l,icell),l=1,atom%Nlevel) !*ntotal
            write(*,'("b: "*(ES14.5E3))') (atom%n(l,icell)/ndag(l),l=1,atom%Nlevel) !*ntotal
            do l=1, atom%Nlevel
                write(*, '(1I3, *(ES14.5E3))') l, (atom%Gamma(l,lp,id), lp=1, atom%Nlevel)
            enddo
            call warning("All pops are negative after SEE!")
            !or error ?
        endif

        if ((is_nan_infinity_vector(atom%n(:,icell))>0)) then
            write(*,*) atom%ID
            write(*,*) "(SEE) BUG pops", " id=",id, " icell=",icell
            write(*,'("ilevel: "*(1I4))') (l, l=1, atom%Nlevel)
            write(*,'("n: "*(ES14.5E3))') (atom%n(l,icell),l=1,atom%Nlevel)
            write(*,'("ndag: "*(ES14.5E3))') (ndag(l),l=1,atom%Nlevel)
            write(*,*) "Gamma:"
            write(*,'(*(I14))') (l, l=1, atom%Nlevel)
            do l=1, atom%Nlevel
                write(*, '(1I3, *(ES14.5E3))') l, (Gamma_dag(l,lp), lp=1, atom%Nlevel)
                write(*, '(1I3, "out:", *(ES14.5E3))') l, (atom%Gamma(l,lp,id), lp=1, atom%Nlevel)
            enddo
            write(*,*) "Rates"
            do l=1, atom%Nline
                write(*,*) " line ", atom%lines(l)%i, atom%lines(l)%j
                write(*,*) "-> Rij=",atom%lines(l)%Rij(id)," Rji=",atom%lines(l)%Rji(id)
                write(*,*) "-> Cij=",atom%lines(l)%Cij(id)," Cji=",atom%lines(l)%Cji(id)
            enddo
            do l=1, atom%Ncont
                write(*,*) " cont ", atom%continua(l)%i, atom%continua(l)%j
                write(*,*) "-> Rij=",atom%continua(l)%Rij(id)," Rji=",atom%continua(l)%Rji(id)
                write(*,*) "-> Cij=",atom%continua(l)%Cij(id)," Cji=",atom%continua(l)%Cji(id)
            enddo
            call error("nan or infinity found in n after SEE!")
        end if

        dM = 0.0_dp
        ndag = ndag * ntotal

        !Handle negative pops and very small populations
        Nsmall_pops = 0
        Nneg_or_null_pops = 0
        lp = 1
        do l=1,atom%Nlevel
            atom%n(l,icell) = atom%n(l,icell) * ntotal

            if (atom%n(l,icell) < frac_limit_pops * ntotal) then
                Nsmall_pops = Nsmall_pops + 1
                level_index(lp) = l
                lp = lp + 1
                if (atom%n(l,icell) <= 0.0_dp) then
                    atom%n(l,icell) = abs(atom%n(l,icell))
                    Nneg_or_null_pops = Nneg_or_null_pops + 1
                endif
            else
                dn_n = (1.0_dp - ndag(l))/atom%n(l,icell)
                dM = max(dM, abs(1.0_dp - ndag(l)/atom%n(l,icell)))
                ! dM = max(dM, abs(atom%n(l,icell)-ndag(l))/ndag(l))
            endif

        enddo
    ! 		if (Nsmall_pops >= 0.75*atom%Nlevel) then
    ! 			write(*,*) "*************** Warning SEE ***************"
    ! 			write(*,'("id #"(1I2)", cell #"(1I6)," ne "(ES14.5E3)" m^-3, nTot "(ES14.5E3)" m^-3")') id, icell, ne(icell), ntotal
    ! 			write(*,'(" for atom "(1A2)" with "(1I3)" levels")') atom%ID, atom%Nlevel
    ! 			write(*,'("--> Found "(1I3)" impurity levels "(1F12.5)" %")') Nsmall_pops, 100.0*real(Nsmall_pops)/real(atom%Nlevel)
    ! 			write(*,'(" --> with "(1I3)" below "(ES14.5E3)" m^-3")') Nneg_or_null_pops, prec_pops * ntotal
    ! 			write(*,'("ilevel: "*(1I4))') (level_index(lp), lp=1, Nsmall_pops)
    ! 			write(*,'("n: "*(ES14.5E3)" m^-3")') (atom%n(level_index(lp),icell), lp=1, Nsmall_pops)
    ! 			write(*,'("n/nTot: "*(ES14.5E3)" m^-3")') (atom%n(level_index(lp),icell)/ntotal, lp=1, Nsmall_pops)
    ! 			write(*,'("n/ne: "*(ES14.5E3)" m^-3")') (atom%n(level_index(lp),icell)/ne(icell), lp=1, Nsmall_pops)
    ! 			write(*,*) "*************** *********** ***************"
    ! 		endif

        if (allocated(n_new)) then
            n_new(atom%activeindex,1:atom%Nlevel,icell) = atom%n(:,icell)
            atom%n(:,icell) = ndag
        endif
        return
    end subroutine see_atom


    subroutine rate_matrix(id)
        integer, intent(in) :: id
        integer :: nact

        do nact=1, NactiveAtoms
            call rate_matrix_atom(id, ActiveAtoms(nact)%p)
        enddo

        return
    end subroutine rate_matrix

    subroutine rate_matrix_atom(id, atom)
    !see Hubeny & Mihalas 2014 eq. 14.8a to 14.8c for the rate matrix elements.
    !Rij and Rji are 0 if lte.
    !atom%cswitch is 1.0 by default
    !TO DO:
    ! occupation probability
        integer, intent(in) :: id
        type (AtomType), intent(inout) :: atom
        integer :: kr, l, i, j
        real(kind=dp) :: nlte_fact
        nlte_fact = 1.0_dp
        if (lforce_lte) nlte_fact = 0.0_dp
        !only Rji are initialised to Aji or Uji for cont and lines!
        !even if lforce_lte.

        atom%Gamma(:,:,id) = 0.0_dp

        do kr=1, atom%Nline

            i = atom%lines(kr)%i; j = atom%lines(kr)%j
            ! write(*,*) i, j
            ! write(*,*) "Rij=",atom%lines(kr)%Rij(id),"Cij=",atom%lines(kr)%Cij(id)
            ! write(*,*)"Rji=", atom%lines(kr)%Rji(id),"Cji=",atom%lines(kr)%Cji(id)

            atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - (atom%lines(kr)%Rij(id)+atom%cswitch*atom%lines(kr)%Cij(id))
            atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - (nlte_fact*atom%lines(kr)%Rji(id)+atom%cswitch*atom%lines(kr)%Cji(id))

        enddo

        do kr=1, atom%Ncont

            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            ! write(*,*) i, j
            ! write(*,*) "Rij=",atom%continua(kr)%Rij(id),"Cij=",atom%continua(kr)%Cij(id)
            ! write(*,*) "Rji=",atom%continua(kr)%Rji(id),"Cji=",atom%continua(kr)%Cji(id)

            atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - (atom%continua(kr)%Rij(id)+atom%cswitch*atom%continua(kr)%Cij(id))
            atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - (nlte_fact*atom%continua(kr)%Rji(id)+atom%cswitch*atom%continua(kr)%Cji(id))

        enddo


        do l = 1, atom%Nlevel
            atom%Gamma(l,l,id) = 0.0_dp
            !Gamma(j,i) = -(Rij + Cij); Gamma(i,j) = -(Rji + Cji)
            !Gamma(i,i) = Rij + Cij = -sum(Gamma(j,i))
            atom%Gamma(l,l,id) = sum(-atom%Gamma(:,l,id)) !positive
        end do


        return
    end subroutine rate_matrix_atom

    subroutine init_rates(id,icell)
        integer, intent(in) :: id, icell
        integer :: n

        do n=1,NactiveAtoms
            call init_rates_atom(id,icell,ActiveAtoms(n)%p)

            !x ne included. Derivatives to ne not included.
            if (activeatoms(n)%p%id=='H') then
                call collision_rates_hydrogen_loc(id,icell)
            else 
                call collision_rates_atom_loc(id,icell,ActiveAtoms(n)%p)
            endif
        enddo

        return 
    end subroutine init_rates

    subroutine init_rates_atom(id,icell,atom)
        integer, intent(in) :: id, icell
        type (atomtype), intent(inout) :: atom 
        integer :: kr, i, j

        do kr=1,atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            atom%continua(kr)%Rij(id) = 0.0_dp
            !updated value of ni and nj!
            atom%continua(kr)%Rji(id) = tab_Uji_cont(kr,atom%activeindex,icell) * atom%nstar(i,icell)/atom%nstar(j,icell)
            atom%continua(kr)%Cij(id) = 0.0_dp
            atom%continua(kr)%Cji(id) = 0.0_dp
        enddo

        do kr=1,atom%Nline
            atom%lines(kr)%Rij(id) = 0.0_dp
            atom%lines(kr)%Rji(id) = atom%lines(kr)%Aji
            atom%lines(kr)%Cij(id) = 0.0_dp
            atom%lines(kr)%Cji(id) = 0.0_dp
        enddo

        return

    end subroutine init_rates_atom


!to do level dissolution
    subroutine accumulate_radrates_mali(id, icell, iray, domega)
    ! --------------------------------------------------------- !
    ! Accumulate radiative rates for the MALI method.
    !  Rates are integrated in frequency and accumulated 
    !  ray-by-ray.
    !
    ! BEWARE: Psi independent of iray (so iray can be different of 1 during sub-iter)
    !
    ! 
    ! --------------------------------------------------------- !
        integer, intent(in) :: id, icell, iray
        real(kind=dp), intent(in) :: dOmega
        real(kind=dp), dimension(Nlambda_max_trans) :: Ieff
        real(kind=dp), dimension(Nlambda_max_line) :: phi0
        type(AtomType), pointer :: atom, atom2
        integer :: kr, i, j, Nl, Nr, Nb, ip, jp, Nrp, Nbp
        integer :: i0, l, nact, krr
        real(kind=dp) :: jbar_up, jbar_down, xcc_down, xcc_up
        real(kind=dp) :: wl, wphi, anu, anu1, ni_on_nj_star
        real(kind=dp) :: ehnukt, ehnukt1

        atom_loop : do nact = 1, Nactiveatoms
            atom => ActiveAtoms(nact)%p

            line_loop : do kr=1, atom%Nline

                Nr = atom%lines(kr)%Nr;Nb = atom%lines(kr)%Nb
                i = atom%lines(kr)%i
                j = atom%lines(kr)%j

                Nl = Nr - Nb + 1
                i0 = Nb - 1

                phi0(1:Nl) = phi_loc(1:Nl,kr,nact,iray,id)

                Jbar_up = 0.0
                xcc_down = 0.0
                wphi = 0.0

                !(Psi - Psi^\ast) eta
                Ieff(1:Nl) = Itot(Nb:Nr,iray,id) - Psi(Nb:Nr,1,id) * eta_atoms(Nb:Nr,nact,id)

                do l=2, nl
                    wl = c_light * (tab_lambda_nm(i0+l) - tab_lambda_nm(i0+l-1))/atom%lines(kr)%lambda0
                    Jbar_up = Jbar_up + 0.5 * (Ieff(l)*phi0(l)+Ieff(l-1)*phi0(l-1)) * wl
                    wphi = wphi + 0.5*(phi0(l)+phi0(l-1)) * wl
                    !xcc_down = xcc_down + 0.5 * ...
                enddo  
                !the weights of this integral should change
                xcc_down = sum(chi_up(Nb:Nr,i,nact,id)*psi(Nb:nR,1,id)*Uji_down(Nb:Nr,j,nact,id))
                
                if (wphi <= 0.0) then
                    call error("critical! wphi in accumulate_rates_mali!")
                endif
                Jbar_up = Jbar_up / wphi

                jbar_down = jbar_up

                !init at Aji
                atom%lines(kr)%Rji(id) = atom%lines(kr)%Rji(id) + dOmega * Jbar_down * atom%lines(kr)%Bji
                atom%lines(kr)%Rij(id) = atom%lines(kr)%Rij(id) + dOmega * Jbar_up * atom%lines(kr)%Bij

                !cross-coupling terms
                xcc_up = 0.0_dp!
                do krr=1, atom%Nline
                    ip = atom%lines(krr)%i
                    jp = atom%lines(krr)%j
                    Nbp = atom%lines(krr)%Nb
                    Nrp = atom%lines(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo
                do krr=1, atom%Ncont
                    ip = atom%continua(krr)%i
                    jp = atom%continua(krr)%j
                    Nbp = atom%continua(krr)%Nb
                    Nrp = atom%continua(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo

                atom%lines(kr)%Rji(id) = atom%lines(kr)%Rji(id) - xcc_down * dOmega
                atom%lines(kr)%Rij(id) = atom%lines(kr)%Rij(id) + xcc_up * dOmega

            end do line_loop

            do kr = 1, atom%Ncont

                i = atom%continua(kr)%i; j = atom%continua(kr)%j

          !ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
                ni_on_nj_star = atom%nstar(i,icell)/atom%nstar(j,icell)

                Nb = atom%continua(kr)%Nb; Nr = atom%continua(kr)%Nr
                Nl = Nr-Nb + 1
                i0 = Nb - 1

                Jbar_down = 0.0
                Jbar_up = 0.0
                xcc_down = 0.0

                Ieff(1:Nl) = Itot(Nb:Nr,iray,id) - Psi(Nb:Nr,1,id) * eta_atoms(Nb:Nr,nact,id)

                do l=2, Nl

                    anu = tab_Vij_cont(l,kr,nact)
                    anu1 = tab_Vij_cont(l-1,kr,nact)

                    wl = (tab_lambda_nm(i0+l) - tab_lambda_nm(i0+l-1))

                    Jbar_up = Jbar_up + 0.5 * (anu*Ieff(l)/tab_lambda_nm(i0+l)+anu1*Ieff(l-1)/tab_lambda_nm(i0+l-1)) * wl

                    ehnukt = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l))/tab_lambda_nm(i0+l)
                    ehnukt1 = exp(-hc_k/T(icell)/tab_lambda_nm(i0+l-1))/tab_lambda_nm(i0+l-1)

                    Jbar_down = jbar_down + 0.5 * ( anu*Ieff(l)*ehnukt + anu1*Ieff(l-1)*ehnukt1 ) * wl

                enddo

                xcc_down = sum(chi_up(Nb:Nr,i,nact,id)*Uji_down(Nb:Nr,j,nact,id)*psi(Nb:Nr,1,id))

                atom%continua(kr)%Rij(id) = atom%continua(kr)%Rij(id) + dOmega*fourpi_h * Jbar_up
                !init at tab_Uji_down(kr,nact,icell) <=> Aji
                atom%continua(kr)%Rji(id) = atom%continua(kr)%Rji(id) + dOmega*fourpi_h * Jbar_down * ni_on_nj_star


                !cross-coupling terms
                xcc_up = 0.0_dp!
                do krr=1, atom%Nline
                    ip = atom%lines(krr)%i
                    jp = atom%lines(krr)%j
                    Nbp = atom%lines(krr)%Nb
                    Nrp = atom%lines(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo
                do krr=1, atom%Ncont
                    ip = atom%continua(krr)%i
                    jp = atom%continua(krr)%j
                    Nbp = atom%continua(krr)%Nb
                    Nrp = atom%continua(krr)%Nr
                    if (jp==i) then !i upper level of another transitions
                        xcc_up = xcc_up + sum(chi_down(Nbp:Nrp,j,nact,id)*psi(Nbp:Nrp,1,id)*Uji_down(Nbp:Nrp,i,nact,id))
                    endif
                enddo

                atom%continua(kr)%Rji(id) = atom%continua(kr)%Rji(id) - xcc_down * dOmega
                atom%continua(kr)%Rij(id) = atom%continua(kr)%Rij(id) + xcc_up * dOmega

            enddo
            atom => NULL()

        end do atom_loop

        return
    end subroutine accumulate_radrates_mali

    subroutine update_populations(id, icell, iterate_ne, delta)
    ! --------------------------------------------------------------------!
    ! Performs a solution of SEE for each atom and simultneously solves
    ! the charge conservation equations for hydrogen and helium.          !
    ! To do: implements Ng's acceleration iteration here
    ! --------------------------------------------------------------------!
        integer, intent(in) :: id, icell
        logical, intent(in) :: iterate_ne
        type(AtomType), pointer :: atom
        integer :: nact, nact_start
        real(kind=dp) :: dM
        real(kind=dp), intent(out) :: delta
        real(kind=dp) :: dne

        delta = 0.0_dp

        if (iterate_ne) then
            !-> currenlty H is always active if He active
            write(*,*) "iterate ne not yet!"
            stop
        else
            nact_start = 1
            atom_loop : do nact=nact_start, Nactiveatoms

                atom => activeatoms(nact)%p
                call SEE_atom(id, icell, atom, dM)
                ! call calc_delta_Tex_atom(icell, atom, dT, Tex, Tion,.false.)

                !max for all atoms
                delta = max(delta, dM)

                atom => null()

            enddo atom_loop

        endif
        return
    end subroutine update_populations


    subroutine radrate_matrix_atom(id, icell, atom, gr, dgrdne)

    !Rji continua are divided by ne !

        integer, intent(in) :: icell, id
        type (AtomType), intent(in) :: atom
        real(kind=dp), intent(out), dimension(atom%Nlevel, atom%Nlevel) :: gr, dgrdne
        integer :: kr, i, j

        gr(:,:) = 0.0_dp
        dgrdne(:,:) = 0.0_dp !derivative of recombinaison rate matrix to dne

        do kr=1, atom%Nline
            i = atom%lines(kr)%i; j = atom%lines(kr)%j
            gr(j,i) = gr(j,i) - atom%lines(kr)%Rij(id)
            gr(i,j) = gr(i,j) - atom%lines(kr)%Rji(id)
        enddo

        do kr=1, atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            gr(j,i) = gr(j,i) - atom%continua(kr)%Rij(id)
            gr(i,j) = gr(i,j) - atom%continua(kr)%Rji(id) / ne(icell)
            dgrdne(i,j) = dgrdne(i,j) - atom%continua(kr)%Rji(id) / ne(icell)
        enddo

        do i=1, atom%Nlevel
            gr(i,i) = -sum(gr(:,i))
            dgrdne(i,i) = -sum(dgrdne(:,i))
        enddo


        return
    end subroutine radrate_matrix_atom

    subroutine colrate_matrix_atom(id, icell, atom, gc, dgcdne)

    !see Hubeny & Mihalas 2014 eq. 14.8b

        integer, intent(in) :: id, icell
        type (AtomType), intent(in) :: atom!, pointer
        real(kind=dp), intent(out), dimension(atom%Nlevel,atom%Nlevel) :: gc, dgcdne
        integer :: kr, i, j

        gc(:,:) = 0.0_dp
        dgcdne(:,:) = 0.0_dp

        do kr=1, atom%Nline
            i = atom%lines(kr)%i; j = atom%lines(kr)%j
            gc(j,i) = gc(j,i) - atom%lines(kr)%Cij(id)
            gc(i,j) = gc(i,j) - atom%lines(kr)%Cji(id)
            dgcdne(i,j) = dgcdne(i,j) - 2.0 * atom%lines(kr)%Cji(id) / ne(icell)
            dgcdne(j,i) = dgcdne(j,i) - atom%lines(kr)%Cij(id) / ne(icell)
        enddo

        do kr=1, atom%Ncont
            i = atom%continua(kr)%i; j = atom%continua(kr)%j
            gc(j,i) = gc(j,i) - atom%continua(kr)%Cij(id)
            gc(i,j) = gc(i,j) - atom%continua(kr)%Cji(id)
            dgcdne(j,i) = dgcdne(j,i) - atom%continua(kr)%Cij(id) / ne(icell)
            dgcdne(i,j) = dgcdne(i,j) - 2.0 * atom%continua(kr)%Cji(id) / ne(icell)
        enddo

        do i=1, atom%Nlevel
            gc(i,i) = -sum(gc(:,i))
            dgcdne(i,i) = -sum(dgcdne(:,i))
        enddo

        return
    end subroutine colrate_matrix_atom

!-> adapt if H is not active
  subroutine particle_conservation_H_and_He (icell, Neq, x, f, df)
    integer, intent(in) :: icell, Neq
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(inout) :: f(neq), df(neq, neq)
    integer :: imaxpop, nlev, l, lp
    real(kind=dp) :: ntotal, ntotal_He

    ntotal = nHtot(icell)

    nlev = Neq - 1 !remove electrons

    imaxpop = locate(x(1:hydrogen%Nlevel), maxval(x(1:hydrogen%Nlevel)))

    !mass conservation: ntotal = sum_i n(i)
    !eq imaxpop is replaced by:
    !f(imaxpop) = ntotal - sum_i n(i) = 0
    !or
    !f(imaxpop) = 1.0 - sum_i n(i) / ntotal = 0

    f(imaxpop) = 1.0_dp
    f(imaxpop) = f(imaxpop) -sum(x(1:hydrogen%Nlevel)) / ntotal !~0, better to set to 0?

    df(imaxpop,:) = 0.0_dp
    df(imaxpop,1:hydrogen%Nlevel) = -1.0_dp / ntotal

    !-> H+He ?
    ! 		df(imaxpop,1:neq-1) = -1.0_dp / (ntotal+ntotal_atom(icell, helium))

    if (helium_is_active) then
        ntotal_He = helium%Abund * nHtot(icell)

       !relative to the sub array !!
       imaxpop = locate(x(hydrogen%Nlevel+1:nlev), maxval(x(hydrogen%Nlevel+1:nlev)))
       imaxpop = hydrogen%Nlevel+imaxpop

       f(imaxpop) = 1.0_dp
       f(imaxpop) = f(imaxpop) -sum(x(hydrogen%Nlevel+1:nlev)) / ntotal_He !~0 better to set 0 ?

       df(imaxpop,:) = 0.0_dp
       df(imaxpop,hydrogen%Nlevel+1:nlev) = -1.0_dp / ntotal_He

    endif

    return
  end subroutine particle_conservation_H_and_He

  subroutine multivariate_newton_raphson (neq, df, f, x)

    !matmul(Jacobian, delta) = - f_equation_vector

    integer, intent(in) :: neq
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
    integer :: ieq, jvar


    do ieq=1, neq
       f(ieq) = f(ieq) * -1.0_dp
       do jvar=1, neq
          df(ieq,jvar) = df(ieq,jvar) * x(jvar)
       enddo
    enddo

    call GaussSlv(df,f,neq)

    return
  end subroutine multivariate_newton_raphson

!-> adapt if H is not active
  !F_cc = 1.0 - 1/ne * (np + nHeII + 2*nHeIII) = 0
!   subroutine non_lte_charge_conservation_H_and_He (icell, neq, x, f, df)
!   !TO DO: move to elecdensity module
!     use atom_type, only : find_continuum
!     !also computes the derivative of the charge conservation equation with
!     !ne and total population of stages with stage > 0 (like nH(Nlevel)=nHII)
!     integer, intent(in) :: icell, neq
!     real(kind=dp), intent(in) :: x(neq)
!     real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
!     integer :: n, j, n_start, j0, j1
!     type (element) :: elem
!     real(kind=dp) :: akj, np, nheII, nheIII
!     real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk
!     !Separate in two parts: 1 non-LTE atoms and 2 electrons from LTE elements

!     !first identify proton levels (hydrogen%Nlevel)
!     !heII levels and heIII levels

!     !derivative of CC (charge conservation) to ne
!     !x = (n1,n2,n3 .... ne)
!     !Ionisation stage of nHII = np is 1, 1 for nHeII and 2 for nHeIII
!     np = hydrogen%n(hydrogen%Nlevel,icell)
!     if (helium_is_active) then
!        !I'm assuming I'm starting from HeI
!        j0 = find_continuum(helium, 1) !ground state of HeII
!        j1 = find_continuum(helium, j0) !heIII
!        nHeII = sum(helium%stage(j0:j1-1)*helium%n(j0:j1-1,icell))
!        nHeIII = helium%stage(j1) * helium%n(j1,icell)

!        n_start = 3
!        df(neq,neq) = 1/x(neq)**2 * (np + nheII + nheIII) !the factor stage is included!!
!        !the sum should return 0 * nHeI + 1 * nHeII + 2 * nHeIII
!        !df(neq,neq) = 1/x(neq)**2 * (np + sum(helium%stage(:) * helium%n(:,icell)))
!        f(neq) = 1.0 - (1.0/x(neq)) * (np + nheII + nheIII) !factor 2 in nheIII included
!        !derivative to protons
!        df(neq,hydrogen%Nlevel) = -1.0 / x(neq)
!        !derivative to nheIII
!        df(neq, neq-1) = -2.0 / x(neq)
!        !derivative to nHeII!
!        !or derivative to nHeII_i individual levels
!        do j=1, j1-j0
!           df(neq, neq-1-j) = -1/x(neq) !or for each level of nHeII I have the same derivative?
!        enddo
!        ! 			stop
!     else
!        n_start = 2
!        df(neq,neq) = 1/x(neq)**2 * np
!        f(neq) = 1.0 - (1.0/x(neq)) * np
!        !derivative to protons
!        df(neq,hydrogen%Nlevel) = -1.0 / x(neq)
!     endif

!     !check:
!     !if (hydrogen%n(hydrogen%Nlevel,icell) /= sum(hydrogen%n(:,icell)*hydrogen%stage(:)))
!     !to do


!     !now contribution from LTE atoms.
!     do n=n_start, 26 !for 26 elements

!         !if H is passive should be computed here
!         !but if helium is active, He should be avoided.
!         !if H is active, n_start avoid or non He if it is active/passive.

!        elem = Elements(n)

!        call ionisation_frac_lte(elem, icell, x(neq), fjk, dfjk)

!        do j=2, elem%Nstage
!           !pure LTE term j = 1 corresponds to the first ionisation stage always
!           !unlike in solve_ne where this loop is also for non-LTE models
!           akj = (j-1) * elem%abund * nhtot(icell) !(j-1) = 0 for neutrals and 1 for singly ionised
!           f(neq) = f(neq) - akj * fjk(j) / x(neq)
!           df(neq,neq) = df(neq,neq) + akj * fjk(j) / x(neq)**2 - akj / x(neq) * dfjk(j)
!        end do

!     enddo


!     return
!   end subroutine non_lte_charge_conservation_H_and_He


! !-> adapt if H is not active
!   subroutine see_ionisation_nonlte_H_and_He(id, icell, dM, dne)

!     integer, intent(in) :: id, icell
!     real(kind=dp), intent(out) :: dM, dne
!     integer :: nact, Neq, nact_start, n_iter, i, kr
!     !smoothing parameter default if damp_char. If native pops, damp the iterations by a facteur damp_scaling larger
!     real(kind=dp), parameter :: damp_scaling = 10.0, damp_char = 5.0, precision = 1d-5
!     real(kind=dp) :: d_damp, nedag, delta_f, dfpop, dfne, dM_he, dfpop_he, dfpop_h
!     real(kind=dp), allocatable :: gamma_r(:,:,:), ndag(:), RR(:,:), CC(:,:), deriv(:,:), dgamrdne(:,:), dgamcdne(:,:), gamma_tot(:,:)
!     real(kind=dp), allocatable :: fvar(:), dfvar(:,:), xvar(:) !Jacobian
!     logical :: lconverged, neg_pops, rest_damping, lsecond_try, verbose
!     integer, parameter :: maxIter = 1000
!     integer :: max_iter, Nlevel_atoms, ihel=0, ih=0, n_active

!     verbose = .false. !debug mode
!     lconverged = .false.
!     rest_damping = .false.
!     d_damp = damp_char
!     lsecond_try = .false. !set to .true. to avoid starting again if the first iterative scheme failed


!     !number of equations = Nlevel_atom1 + Nlevel_atom2 + ... + 1 for ne
!     Nlevel_atoms = hydrogen%Nlevel
!     n_active = 1
!     ih = 1
!     if (helium_is_active) then
!        ihel = Nlevel_atoms + 1 !position of SEE of helium among
!                                !all levels
!        Nlevel_atoms = Nlevel_atoms + helium%Nlevel !sum not max
!        !SEE(H) -> ih : hydrogen%Nlevel
!        !SEE(He) -> ihel : ihel + helium%Nlevel
!        n_active = n_active + 1
!     endif
!     !add the contributions of charge conservation.
!     Neq = Nlevel_atoms + 1
!     !currently allocate and deallocate for each cell
!     !futur fixed size
!     !unknowns = populations of atomic levels + electron density
!     allocate(xvar(neq)); xvar(:) = 0.0_dp
!     !equations = Nlevel SEE per atom  + charge conservation (CC)
!     allocate(fvar(neq)); fvar(:) = 0.0_dp
!     !Jacobian of the systems = derivative of i equations to xvar of sum_j partial(fvar(i))/partial(xvar(j))
!     allocate(dfvar(neq,neq)); dfvar(:,:) = 0.0_dp
!     !old ne values to rest for the other cell. The new value is stored in ne_new for the next iterations
!     nedag = ne(icell)

!     allocate(ndag(Nlevel_atoms)) !hydrogen%Nlevel if only H
!     ndag(1:hydrogen%Nlevel) = hydrogen%n(:,icell)
!     if (helium_is_active) then
!        ndag(ihel:Nlevel_atoms) = helium%n(:,icell)
!     endif

!     !temporary rate deriv matrix
!     allocate(deriv(Nlevel_Atoms, Nlevel_atoms))


!     !local radiative rates and collisional rates matrix to build total rate matrix
!     !for current estimate of ne.
!     allocate(RR(Nlevel_atoms,Nlevel_atoms), CC(Nlevel_atoms,Nlevel_atoms), dgamcdne(Nlevel_atoms,Nlevel_atoms))
!     !ne dependent radiative rate matrix for each atom + derivative to ne of rate matrix for all equations (of all atoms)
!     allocate(gamma_r(n_active, Nlevel_atoms,Nlevel_atoms),dgamrdne(Nlevel_atoms,Nlevel_atoms))
!     dgamrdne = 0.0_dp; dgamcdne = 0.0_dp
!     call radrate_matrix_atom(id, icell, hydrogen, gamma_r(1,1:hydrogen%Nlevel,1:hydrogen%Nlevel), &
!          deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel))
!     dgamrdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel) = deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel)
!     if (helium_is_active) then
!        deriv = 0.0_dp
!        call radrate_matrix_atom(id, icell, helium, gamma_r(n_active,1:helium%Nlevel,1:helium%Nlevel), &
!             deriv(1:helium%Nlevel,1:helium%Nlevel))
!        dgamrdne(ihel:Nlevel_atoms,ihel:Nlevel_atoms) = deriv(1:helium%Nlevel,1:helium%Nlevel)
!     endif

!     !start Newton-Raphson iterations
!     if (verbose) write(*,*) "T = ", T(icell), " nHtot = ", nHtot(icell)

!     n_iter = 0
!     max_iter = maxIter

!     iterative_loop : do while (.not.lconverged)
!        delta_f = 0.0_dp
!        dfpop = 0.0_dp !only populations
!        dfne = 0.0_dp
!        dfpop_he = 0.0_dp
!        dfpop_h = 0.0_dp

!        n_iter = n_iter + 1

!        !always the last one
!        xvar(Neq) = ne(icell)
!        fvar(:) = 0.0
!        dfvar(:,:) = 0.0

!        !set unknowns to old values
!        !First hydrogen
!        xvar(1:hydrogen%Nlevel) = hydrogen%n(:,icell)
!        RR(:,:) = gamma_r(1,:,:) !init to radiative rates
!        !we multiply ne by zero if lforce_lte.
!        do kr=1,hydrogen%Ncont
!           RR(hydrogen%continua(kr)%i,hydrogen%continua(kr)%j) = RR(hydrogen%continua(kr)%i,hydrogen%continua(kr)%j) * ne(icell)
!        enddo
 
!        call collision_rates_hydrogen_loc(id,icell)
!        call collrate_matrix_atom(id, icell, hydrogen, CC(1:hydrogen%Nlevel,1:hydrogen%Nlevel), &
!             deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel), dgamcdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel))

!        ! 			gamma_tot(1:hydrogen%Nlevel,1:hydrogen%nlevel) = RR(1:hydrogen%Nlevel,1:hydrogen%Nlevel) + CC(1:hydrogen%Nlevel,1:hydrogen%Nlevel)
!        hydrogen%Gamma(:,:,id) = RR(1:hydrogen%Nlevel,1:hydrogen%Nlevel) + CC(1:hydrogen%Nlevel,1:hydrogen%Nlevel)

!     !    !now compute the diagonal elements
!     !    do i=1, hydrogen%Nlevel !diagonal hence full rate matrix
!     !       ! 				gamma_tot(i,i) = sum(-gamma_tot(1:hydrogen%Nlevel,i)) !positive
!     !       hydrogen%gamma(i,i,id) = -sum(hydrogen%gamma(:,i,id))
!     !    enddo
!        !stacking equations for helium
!        !now helium
!        if (helium_is_active) then
!           xvar(ihel:Nlevel_atoms) = helium%n(:,icell)
!           RR(:,:) = gamma_r(n_active,:,:) !init to radiative rates
!           !we multiply ne by zero if lforce_lte.
!           do kr=1,helium%Ncont
!              !Rji in rate matrix is at rate(i,j)
!              RR(helium%continua(kr)%i,helium%continua(kr)%j) = RR(helium%continua(kr)%i,helium%continua(kr)%j) * ne(icell)
!           enddo
!           !cumulate the derivative for helium equations
!           call collision_rates_atom_loc(id, icell, helium)
!           call collrate_matrix_atom(id, icell, helium, CC(1:helium%Nlevel,1:helium%Nlevel), deriv(1:helium%Nlevel,1:helium%Nlevel), &
!                dgamcdne(ihel:Nlevel_atoms,ihel:Nlevel_atoms))

!           ! 				gamma_tot(ihel:Nlevel_atoms,ihel:Nlevel_atoms) = RR(1:helium%Nlevel,1:helium%Nlevel) + CC(1:helium%Nlevel,1:helium%Nlevel)
!           helium%gamma(:,:,id) = RR(1:helium%Nlevel,1:helium%Nlevel) + CC(1:helium%Nlevel,1:helium%Nlevel)
!         !   do i=1, helium%Nlevel !diagonal hence full rate matrix
!         !      ! 					gamma_tot(ihel+i - 1,ihel+i - 1) = sum(-gamma_tot(ihel:Nlevel_atoms,ihel+i - 1)) !positive
!         !      helium%gamma(i,i,id) = sum(-helium%gamma(:,i,id)) !positive
!         !   enddo
!        endif
!        !gamma_tot is 0 for SEE of H if N>Nlevel_H
!        !so that the product with all pops (including nHe pops) is 0.

!        !rate equation for this atom stored in f and df !
!        !an equation per level dn_i/dt = 0 = sum_lp n_l * Gamma_lp_l
!        call rate_equations_H_and_He(id, icell, Neq, ihel, dgamrdne(:,:), dgamcdne(:,:), xvar, fvar, dfvar)

!        !charge conservation!
!        call non_lte_charge_conservation_H_and_He (icell, neq, xvar, fvar, dfvar)

!        !replace one equation of SEE by particle number conservation
!        !particule conservation!

!        call particle_conservation_H_and_He (icell, Neq, xvar, fvar, dfvar)

!        !newton raphson!
!        call multivariate_newton_raphson (neq, dfvar, fvar, xvar)
!        ! 			call multivariate_newton_raphson (neq-1, dfvar(1:neq-1,1:neq-1), fvar(1:neq-1), xvar(1:neq-1))

!        !update atomic populations and ne
!        neg_pops = .false.
!        do i=1, hydrogen%Nlevel
!           hydrogen%n(i,icell) = hydrogen%n(i,icell) * ( 1.0 + fvar(i)/(1.0 + d_damp * abs(fvar(i))) )
!           if (hydrogen%n(i,icell) < 0.0) neg_pops = .true.
!        enddo
!        if (helium_is_active) then
!           do i=1, helium%Nlevel
!              helium%n(i,icell) = helium%n(i,icell) * ( 1.0 + fvar(ihel+i-1)/(1.0 + d_damp * abs(fvar(ihel+i-1))) )
!              if (helium%n(i,icell) < 0.0) neg_pops = .true.
!           enddo
!        endif
!        !keep ne constant for tests
!        ne(icell) = ne(icell) * ( 1.0 + fvar(neq)/(1.0 + d_damp * abs(fvar(neq))) )


!        if ( (ne(icell) < 1d-16 * hydrogen%n(hydrogen%Nlevel,icell)).or.(neg_pops) ) rest_damping = .true.

!        if (rest_damping .and. d_damp < (damp_char + 1.0)) then !restart with more iterations and larger damping (more stable, slower convergence)
!           ! 				do i=1, hydrogen%Nlevel
!           ! 					xvar(i) = ndag(i)
!           ! 				enddo
!           ! 				if (helium_is_active) then
!           ! 					do i=1, helium%Nlevel
!           ! 						xvar(ihel+i) = ndag(ihel+i)
!           ! 					enddo
!           ! 				endif
!           xvar(1:neq-1) = ndag(:)
!           xvar(neq) = nedag
!           d_damp = damp_char * damp_scaling
!           rest_damping = .false.
!           n_iter = 0
!           max_iter = damp_scaling * maxIter !can be large !

!        elseif (rest_damping) then
!           neg_pops = .false.
!           lconverged = .false.
!        endif

!        !should be fractional
!        delta_f = max(delta_f, maxval(abs(fvar)))
!        dfpop = max(dfpop, maxval(abs(fvar(1:Neq-1))))
!        dfne = max(dfne, abs(fvar(neq)))
!        dfpop_h = max(dfpop_h, maxval(abs(fvar(1:hydrogen%Nlevel)))) !same as dfpop if only H
!        if (helium_is_active) dfpop_he = max(dfpop_he,  maxval(abs(fvar(ihel:Neq-1))))

!        if (verbose) then
!           write(*,'("niter #"(1I4))') n_iter
!           write(*,'("non-LTE ionisation delta="(1ES17.8E3)" dfH="(1ES17.8E3)" dfne="(1ES17.8E3)" dfHe="(1ES17.8E3) )') &
!                delta_f, dfpop_H, dfne, dfpop_he
!        endif

!        if (n_iter > max_iter) then
!           if (lsecond_try) then
!              lconverged = .true.
!              write(*,*) "Not enough iterations to converge", max_iter
!              write(*,*) "err = ", delta_f
!              write(*,*) " cell", icell, " id", id
!              write(*,*) ""
!           else
!              lsecond_try = .true.
!              lconverged = .false.
!              n_iter = 0
!              max_iter = damp_scaling * maxIter
!              d_damp = damp_scaling * damp_char
!              xvar(1:neq-1) = ndag(:) !should include helium
!              xvar(neq) = nedag
!           endif
!        endif

!        if ((delta_f < precision).and..not.rest_damping) then
!           lconverged = .true.
!           exit iterative_loop
!        endif

!     enddo iterative_loop !on convergence

!     dne = abs(1.0 - nedag/ne(icell))
!     dM = maxval(abs(1.0 - ndag(1:hydrogen%Nlevel)/hydrogen%n(:,icell)))
!     if (helium_is_active) then
!        dM_he = maxval(abs(1.0 - ndag(ihel:Nlevel_atoms)/helium%n(:,icell)))
!        !!write(*,*) maxval(ndag(ihel:Nlevel_atoms)), maxval(helium%n(:,icell))
!        if (verbose) write(*,'("(DELTA) non-LTE ionisation dM_H="(1ES17.8E3)" dne="(1ES17.8E3)" dM_He="(1ES17.8E3) )') &
!             dM, dne, dM_he
!        dM = max(dM, dM_he)
!     else
!        if (verbose) write(*,'("(DELTA) non-LTE ionisation dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne
!     endif

!     n_new(1,1:hydrogen%Nlevel,icell) = hydrogen%n(:,icell)
!     ne_new(icell) = ne(icell) !will be set to ne once the new populations on all cells (with old quantities) have bee computed
!     hydrogen%n(:,icell) = ndag(1:hydrogen%Nlevel) !the first value before iterations
!     !needed for global convergence

!     ne(icell) = nedag !reset because we have to loop over all cells
!     !and we don't want to change the old values by the new one
!     !until all  cells are treated
!     if (helium_is_active) then
!        n_new(helium%activeindex, 1:helium%Nlevel,icell) = helium%n(:,icell)
!        helium%n(:,icell) = ndag(ihel:Nlevel_atoms)
!     endif

!     deallocate(xvar,fvar,dfvar)
!     deallocate(gamma_r, ndag, RR, CC, deriv, dgamrdne, dgamcdne)
!     !deallocate(gamma_tot)

!     if (verbose) stop
!     return
!   end subroutine see_ionisation_nonlte_H_and_He

!  subroutine rate_equations_H_and_He(id,icell, neq, ihel, dgrdne, dgcdne, x, f, df)
!     integer, intent(in) :: id, icell, neq, ihel
!     real(kind=dp), intent(in) :: x(neq)
!     real(kind=dp), intent(in), dimension(neq-1,neq-1) :: dgrdne, dgcdne
!     real(kind=dp), intent(out) ::  f(neq), df(neq,neq)
!     integer :: i, j, l, lp, Nlevels
!     real(kind=dp) :: part1, part2
!     !derivative of SEE to ne
!     !dCij/dne = Cij/ne
!     !dCji/dne = 2*Cji/ne
!     !dRji_cont/dne = Rji_cont/ne

!     !Have to do better, to much avoidable loops!!

!     Nlevels = neq - 1 !should be hydrogen%Nlevel + helium%Nlevel

!     !f and df 0 where helium or H levels stops

!     !equations and derivative for H
!     f(1:hydrogen%Nlevel) = matmul(hydrogen%Gamma(:,:,id), x(1:hydrogen%Nlevel))
!     !derivative to levels n
!     df(1:hydrogen%Nlevel,1:hydrogen%Nlevel) = hydrogen%gamma(:,:,id)
!     !derivative to ne
!     df(1:hydrogen%Nlevel,neq) = matmul( dgcdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel), x(1:hydrogen%Nlevel) )

!     ! do i=1, hydrogen%Nlevel
!     ! 	write(*,*) i, " eq deriv:", (df(i, j), j=1, neq)
!     ! enddo
!     if (helium_is_active) then
!        !He
!        f(ihel:Nlevels) = matmul(helium%gamma(:,:,id), x(ihel:Nlevels))
!        df(ihel:Nlevels,ihel:Nlevels) = helium%Gamma(:,:,id)
!        df(ihel:Nlevels,neq) = matmul (dgcdne(ihel:Nlevels,ihel:Nlevels), x(ihel:Nlevels))
!     endif

!     ! do i=1, helium%Nlevel
!     ! 	write(*,*) i, " eq deriv:", (df(i-1+ihel, j), j=1, neq)
!     ! enddo
!     ! stop
!     if (lforce_lte) return

!     !radiative part
!     !H
!     df(1:hydrogen%Nlevel,neq) = df(1:hydrogen%Nlevel,neq) + matmul( dgrdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel), &
!          x(1:hydrogen%Nlevel) )
!     if (helium_is_active) then
!        !He
!        df(ihel:Nlevels,neq) = df(ihel:Nlevels,neq) + matmul(dgrdne(ihel:Nlevels,ihel:Nlevels), x(ihel:Nlevels))
!     endif

!     return

!     return
!   end subroutine rate_equations_H_and_He

 
!    subroutine opacity_atom_bb_subiter(id, icell, iray, chi, Snu)
!    !obtain bound-bound opacities
!       integer, intent(in) :: id, icell, iray
!       real(kind=dp), intent(inout), dimension(n_lambda) :: chi, Snu
!       integer :: nat, Nred, Nblue, kr, i, j, Nlam
!       type(AtomType), pointer :: atom

!       atom_loop : do nat = 1, N_Atoms
!          atom => Atoms(nat)%p

!          tr_loop : do kr = 1,atom%Nline

!             if (.not.atom%lines(kr)%lcontrib) cycle


!             Nred = atom%lines(kr)%Nr; Nblue = atom%lines(kr)%Nb
!             Nlam = atom%lines(kr)%Nlambda
!             i = atom%lines(kr)%i; j = atom%lines(kr)%j

!             if ((atom%n(i,icell) - atom%n(j,icell)*atom%lines(kr)%gij) <= 0.0_dp) then
!                cycle tr_loop
!             endif
            
!             chi(Nblue:Nred) = chi(Nblue:Nred) + &
!                hc_fourPI * atom%lines(kr)%Bij * phi_loc(1:Nlam,kr,nact,iray,id) * (n_new(i,nact,icell) - atom%lines(kr)%gij*n_new(j,nact,icell))

!             Snu(Nblue:Nred) = Snu(Nblue:Nred) + &
!                hc_fourPI * atom%lines(kr)%Aji * phi_loc(1:Nlam,kr,nact,iray,id) * n_new(j,nact,icell)

!          end do tr_loop

!          atom => null()

!       end do atom_loop


!       return
!    end subroutine opacity_atom_bb_subiter

end module see