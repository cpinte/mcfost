MODULE statequil_atoms

  use atmos_type, only : ne, T, nTotal_atom, ActiveAtoms, Atoms, Natom, Nactiveatoms, nHmin, ds, elements, hydrogen, helium, &
       nHtot, helium_is_active
  use atom_type
  use spectrum_type, only					: lambda, Nlambda, Nlambda_cont, lambda_cont, Jnu_cont, Itot, &
       dk, dk_min, dk_max, &
       Jnu, eta_c, sca_c, chi_c, chi_c_nlte, eta_c_nlte, eta0_bb, chi0_bb
  use constant
  use opacity, only 						: eta_atoms, Uji_down, chi_down, chi_up, R_xcc
  use math, only		: locate, any_nan_infinity_matrix, any_nan_infinity_vector, is_nan_infinity, solve_lin
  use parametres, only 					: ldissolve, lelectron_scattering, n_cells, lforce_lte, lsobolev_regime
  use collision, only						: collision_rates_atom_loc!, CollisionRate_old
  use impact, only						: Collision_Hydrogen
  use getlambda, only						: hv, Nlambda_max_trans
  use occupation_probability, only 		: D_i, wocc_n
  use opacity, only						: prec_pops, frac_limit_pops, frac_ne_limit
  use lte, only 							: getPartitionFunctionk, sahaeq!, phi_T
  use solvene, only 						: Max_ionisation_stage
  use background_opacity, only			: Thomson

  use mcfost_env, only					: dp
  use constantes, only					: tiny_dp, sigma
  use messages, only 						: error, warning
  use utils, only 						: gaussslv, interp
  use input, only 						: nb_proc

  IMPLICIT NONE
  logical :: ldamp_jacobi = .false.
  !Nlambda,  Nproc not stored for all ray since no sub it and ray-ray building of rate matrix
  real(kind=dp), allocatable :: psi(:,:,:), omega_sor_atom(:)
  real(kind=dp), parameter :: Tmax = 1d10 !value of temperature of transition if ni = nj*gij
  real(kind=dp), allocatable :: n_new(:,:,:), ne_new(:)
  character(len=50), parameter :: invpop_file = "inversion_populations.txt"
  character(len=50), parameter :: profiles_file = "line_profile.txt"
  character(len=50), parameter :: convergence_file = "thresh.txt", sub_convergence_file = "subthresh.txt"
  integer, parameter :: unit_invfile = 20, unit_profiles = 25, unit_conver = 30, unit_subconver = 35

CONTAINS


  !n_new must not be allocated to fill atom%n
  !see see_atom (otherwise new pop = n_new and atom%n=nstar
  subroutine radiation_free_pops_atom(id, icell, atom, verbose)

    !still with Collisions
    !Rij = 0 for both line and cont
    !Rji = Aji for lines
    !Rji = 4pi/h int_nu (dnu exp(-hnu/kt) * alpha_nu * towhnu3_c2) * nstar(i)/nstar(j) for continua

    integer, intent(in) :: id, icell
    type (AtomType), intent(inout) :: atom
    logical, intent(in) :: verbose
    integer :: kr, kc, l, Nl, icell_d, Nblue, Nred, i, j
    real(kind=dp) :: ehnukt, twohnu3_c2, dM, wl, sigma, JJ

    !if (allocated(n_new)) call error("n_new must not be allocated with ZERO_RAD method")

    call collision_matrix_atom(id, icell, atom)
    call initGamma_atom(id, atom)

    tr_loop : do kr=1, atom%Ntr
       kc = atom%at(kr)%ik

       SELECT CASE (atom%at(kr)%trtype)

       CASE ("ATOMIC_LINE")
          i = atom%lines(kc)%i; j = atom%lines(kc)%j
          ! 				wj = 1.0
          ! 				wi = 1.0

          atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - atom%lines(kc)%Aji

       CASE ("ATOMIC_CONTINUUM")
          i = atom%continua(kc)%i; j = atom%continua(kc)%j
          Nblue = atom%continua(kc)%Nb; Nred = atom%continua(kc)%Nr
          Nl = Nred-Nblue + 1

          JJ = 0.0

          icell_d = 1
          if (ldissolve) then
             if (atom%ID=="H") icell_d = icell
          endif

          do l=1, Nl
             if (l==1) then
                wl = 0.5*(lambda(Nblue+l)-lambda(Nblue)) / lambda(Nblue)
             elseif (l==Nl) then
                wl = 0.5*(lambda(Nred)-lambda(Nred-1))  / lambda(Nred)
             else
                wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-2)) / lambda(Nblue+l-1)
             endif

             sigma = atom%continua(kc)%alpha_nu(l,icell_d)

             ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-1))
             twohnu3_c2 = twohc/lambda(Nblue+l-1)**3

             JJ = JJ + wl * twohnu3_c2  * ehnukt * sigma

          enddo
          atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - fourpi_h * atom%nstar(i,icell)/atom%nstar(j,icell) * JJ

       CASE DEFAULT

          CALL Error("Unkown transition type", atom%at(kr)%trtype)

       END SELECT

    end do tr_loop

    call see_atom(id, icell, atom, dM)

    if (verbose) then
       write(*,*) icell, " T=", T(icell)
       write(*,*) "nstar:", (atom%nstar(l,icell), l=1, atom%Nlevel)
       write(*,*) " n:", (atom%n(l,icell), l=1, atom%Nlevel)
       write(*,*) "ratio:", (atom%n(l,icell)/atom%nstar(l,icell), l=1, atom%Nlevel)
       write(*,*) "diff:", dM
    endif

    return
  end subroutine radiation_free_pops_atom


  SUBROUTINE fill_collision_matrix(id, icell)
    integer, intent(in) :: icell, id
    type (AtomType), pointer :: atom
    integer :: nact

    do nact=1, NactiveAtoms
       atom => ActiveAtoms(nact)%ptr_atom

       CALL collision_matrix_atom(id, icell, atom)

       atom => NULL()
    enddo

    RETURN
  END SUBROUTINE fill_collision_matrix



  !Occupa prob ?
  SUBROUTINE collision_matrix_atom(id, icell, atom, deriv)
    integer, intent(in) :: icell, id
    type (AtomType), intent(inout) :: atom  !, pointer
    real(kind=dp), intent(out), dimension(atom%Nlevel,atom%Nlevel), optional :: deriv


    atom%C(:,:,id) = 0d0
    if (atom%ID=="H") then
       !7/03/2021, factorize out ne
       atom%C(:,:,id) = ne(icell) * Collision_Hydrogen(icell,deriv)
       ! 			if (cswitch > 1.0) atom%C(:,:,id) = atom%C(:,:,id) * cswitch
    else
       ! 			write(*,*) "Collision RH not set for large atoms"
       ! 			stop
       ! 			if (id==1) write(*,*) " Using RH routine for collision atom ", atom%ID
       !write(*,*) "First test it to retrive lte results and compare with RH!"
       !stop
       ! 			atom%C(:,:,id) = CollisionRate_old(icell, atom)
       atom%C(:,:,id) = collision_rates_atom_loc(icell, atom, deriv)
    end if
    if (atom%cswitch > 1.0) atom%C(:,:,id) = atom%C(:,:,id) * atom%cswitch


    if (any_nan_infinity_matrix(atom%C(:,:,id))>0) then
       write(*,*) atom%C(:,:,id)
       write(*,*) "Bug at compute_collision_matrix", id, icell, atom%ID, atom%n(:,icell)
       stop
    endif


    RETURN
  END SUBROUTINE collision_matrix_atom


  SUBROUTINE initGamma_atom(id, atom)
    !see Hubeny & Mihalas 2014 eq. 14.8b

    integer, intent(in) :: id
    type (AtomType), intent(inout) :: atom!, pointer
    integer :: l,lp
    ! 		real(kind=dp) :: Ctest(atom%Nlevel,atom%nlevel)

    !because atom%C(i,j) is the collision rate from state i to j
    !but rate matrix A(i,j) = -C(j,i)
    atom%Gamma(:,:,id) = 0.0_dp
    !atom%Gamma(:,:,id) = -1.0_dp * transpose(atom%C(:,:,id))
    do l=1,atom%Nlevel-1
       do lp=l+1,atom%Nlevel
          atom%Gamma(l,lp,id) = -atom%C(lp,l,id)
          atom%Gamma(lp,l,id) = -atom%C(l,lp,id)
       enddo
    enddo
    ! 		Ctest = -transpose(atom%C(:,:,id))
    ! 		do l=1,atom%nlevel
    ! 			write(*,*) (atom%Gamma(l,lp,id),lp=1,atom%nlevel)
    ! 			write(*,*) (atom%C(l,lp,id),lp=1,atom%nlevel)
    ! 			write(*,*) (Ctest(l,lp),lp=1,atom%Nlevel)
    ! 		enddo
    ! 		stop

    call check_elements(id, id, atom%Nlevel,atom%Gamma(:,:,id), "in initGamma_atom")

    RETURN
  END SUBROUTINE initGamma_atom

  SUBROUTINE initGamma(id, init_bb_only)
    ! ------------------------------------------------------ !
    !for each active atom, allocate and init Gamma matrix
    !deallocated in freeAtom()
    ! Collision matrix has to be allocated
    ! if already allocated, simply set Gamma to its new icell
    ! value.
    !
    ! n(l)*Cl->lp = n(lp)*Clp->l
    ! e.g.,
    ! Cij = Cji * (nj/ni)^star, with Cji = C(ij) =
    ! colision rate from j to i.
    !
    ! (ij = j->i = (i-1)*N + j)
    ! (ji = i->j = (j-1)*N + i)
    ! ------------------------------------------------------ !
    integer :: nact, Nlevel, lp, l, ij, ji, nati, natf
    integer, intent(in) :: id
    logical, intent(in), optional :: init_bb_only
    type (AtomType), pointer :: atom
    logical :: init_bb_o

    if (present(init_bb_only)) then
       init_bb_o = init_bb_only
    else
       init_bb_o = .false.
    endif

    nati = 1; natf = Nactiveatoms
    do nact = nati, natf
       atom => ActiveAtoms(nact)%ptr_atom

       call initGamma_atom(id, atom)
       if (init_bb_o) then
          call init_bb_rates_atom(id, atom)
       else
          call init_rates_atom(id,atom)
       endif

       NULLIFY(atom)
    enddo
    RETURN
  END SUBROUTINE initGamma


  SUBROUTINE init_rates_atom(id, atom)
    integer, intent(in) :: id
    type (AtomType), intent(inout) :: atom
    integer :: k

    do k=1, atom%Nline

       atom%lines(k)%Rij(id) = 0.0
       atom%lines(k)%Rji(id) = atom%lines(k)%Aji

    enddo

    do k=1, atom%Ncont

       atom%continua(k)%Rij(id) = 0.0
       atom%continua(k)%Rji(id) = 0.0

    enddo


    RETURN
  END SUBROUTINE init_rates_atom

  SUBROUTINE init_bb_rates_atom(id, atom)
    integer, intent(in) :: id
    type (AtomType), intent(inout) :: atom
    integer :: k

    do k=1, atom%Nline

       atom%lines(k)%Rij(id) = 0.0
       atom%lines(k)%Rji(id) = atom%lines(k)%Aji

    enddo

    RETURN
  END SUBROUTINE init_bb_rates_atom

  !building inclued n_new, new wave integ
  subroutine calc_rates_draft(id, icell, iray, n_rayons)
    integer, intent(in) :: id, icell, iray, n_rayons
    integer :: iray_p
    real(kind=dp) :: chicont, etacont
    real(kind=dp), dimension(Nlambda) :: chi_tot, eta_tot, Ieff
    integer :: nact, Nred, Nblue, kc, kr, i, j, l, a0, Nl, icell_d, imu, iphi
    real(kind=dp) :: etau, wphi, JJ, JJb, j1, j2, I1, I2, ehnukt, twohnu3_c2
    real(kind=dp) :: wl, a1, a2
    type(AtomType), pointer :: aatom, atom

    real(kind=dp) :: wi, wj, nn, dk0

    !should also contains lines
    chi_tot(:) = chi0_bb(:,icell)
    !add continuum in
    eta_tot(:) = eta0_bb(:,icell)


    if (lelectron_scattering) then
    	write(*,*) "calc_rates_draft, elec scatt not available!"
    	stop
!        eta_tot(:) = eta_tot(:) + thomson(icell) * Jnu(:,id) !eta_es(:,icell)
    endif


    atom_loop : do nact = 1, Natom
       atom => Atoms(nact)%ptr_atom

       tr_loop : do kr = 1, atom%Ntr_line


          kc = atom%at(kr)%ik !relative index of a transition among continua or lines

          Nred = atom%lines(kc)%Nred;Nblue = atom%lines(kc)%Nblue
          i = atom%lines(kc)%i;j = atom%lines(kc)%j

          wj = 1.0; wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn

                wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1),hydrogen%n(1,icell)) !1 for H
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1),hydrogen%n(1,icell))
             endif
          endif

          !if ((atom%n(i,icell)*wj/wi - atom%n(j,icell)*atom%lines(kc)%gij) > 0.0_dp) then
          if (n_new(atom%activeindex,i,icell)*wj/wi - n_new(atom%activeindex,j,icell)*atom%lines(kc)%gij > 0.0_dp) then

             chi_tot(Nblue+dk_min:Nred+dk_max) = chi_tot(Nblue+dk_min:Nred+dk_max) + &
                  hc_fourPI * atom%lines(kc)%Bij * atom%lines(kc)%phi_loc(:,iray,id) * &
                  ( n_new(atom%activeindex,i,icell)*wj/wi - n_new(atom%activeindex,j,icell)*atom%lines(kc)%gij )
             !(atom%n(i,icell)*wj/wi - atom%lines(kc)%gij*atom%n(j,icell))

             eta_tot(Nblue+dk_min:Nred+dk_max)= eta_tot(Nblue+dk_min:Nred+dk_max) + &
                  hc_fourPI * atom%lines(kc)%Aji * atom%lines(kc)%phi_loc(:,iray,id) * n_new(atom%activeindex,j,icell)
             !atom%n(j,icell)


          else !neg or null

             eta_tot(Nblue+dk_min:Nred+dk_max)= eta_tot(Nblue+dk_min:Nred+dk_max) + &
                  atom%lines(kc)%Aji * hc_fourPI * atom%lines(kc)%phi_loc(:,iray,id) * n_new(atom%activeindex,j,icell)
             !atom%n(j,icell)
             ! 					eta_tot(Nblue+dk_min:Nred+dk_max)= eta_tot(Nblue+dk_min:Nred+dk_max) + 0.0_dp
             chi_tot(Nblue+dk_min:Nred+dk_max) = chi_tot(Nblue+dk_min:Nred+dk_max) - &
                  hc_fourPI * atom%lines(kc)%Bij * atom%lines(kc)%phi_loc(:,iray,id) * &
                  ( n_new(atom%activeindex,i,icell)*wj/wi - n_new(atom%activeindex,j,icell)*atom%lines(kc)%gij )
             !(atom%n(i,icell)*wj/wi - atom%lines(kc)%gij*atom%n(j,icell))
             !
             !
          endif


       end do tr_loop

       atom => NULL()

    end do atom_loop
    Ieff = ( Itot(:,iray,id) * exp(-ds(iray,id)*chi_tot(:)) + (eta_tot(:)/chi_tot(:)) * (1.0_dp -  exp(-ds(iray,id)*chi_tot(:))) )
    if (any_nan_infinity_vector(Ieff) /= 0.0) then
       write(*,*) " Ieff", Ieff, " Ieff"
       write(*,*) " chi0_bb", chi0_bb(:,icell), 'chi0_bb'
       write(*,*) " eta0_bb", eta0_bb(:,icell), 'eta0_bb'
       stop
    endif
    !now start integrates bound rates
    aatom_loop : do nact = 1, Nactiveatoms
       aatom => ActiveAtoms(nact)%ptr_atom

       atr_loop : do kr = 1, aatom%Ntr_line

          kc = aatom%at(kr)%ik

          Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
          i = aatom%lines(kc)%i
          j = aatom%lines(kc)%j

          Nl = Nred-dk_min+dk_max-Nblue+1
          a0 = Nblue+dk_min-1

          JJ = 0.0
          wphi = 0.0

          do l=2, Nl

             wphi = wphi + 0.5 * (aatom%lines(kc)%phi_loc(l,iray,id)+aatom%lines(kc)%phi_loc(l-1,iray,id)) * 1e3*hv

             J1 = Ieff(a0+l) * aatom%lines(kc)%phi_loc(l,iray,id)
             J2 = Ieff(a0+l-1) * aatom%lines(kc)%phi_loc(l-1,iray,id)

             JJ = JJ + 0.5 * (J1 + J2) * hv * 1d3
             !write(*,*) iray, l, J1, J2, wphi

          enddo

          JJ = JJ / wphi / n_rayons


          !init at Aji
          aatom%lines(kc)%Rji(id) = aatom%lines(kc)%Rji(id) + JJ * aatom%lines(kc)%Bji
          aatom%lines(kc)%Rij(id) = aatom%lines(kc)%Rij(id) + JJ * aatom%lines(kc)%Bij

          ! 				if ((any_nan_infinity_vector(aatom%lines(kc)%Rji) /= 0.0)) then
          ! 					write(*,*) "icell=",icell," id=",id," nb_proc=", nb_proc, aatom%lines(kc)%Rji(id)
          ! 					write(*,*) "Rji=",aatom%lines(kc)%Rji
          ! 					stop
          ! 				endif
          ! 				if ((any_nan_infinity_vector(aatom%lines(kc)%Rij) /= 0.0)) then
          ! 					write(*,*) "icell=",icell," id=",id," nb_proc=", nb_proc, aatom%lines(kc)%Rij(id)
          ! 					write(*,*) " Rij=", aatom%lines(kc)%Rij
          ! 					stop
          ! 				endif

       end do atr_loop

       atrc_loop : do kr = aatom%Ntr_line+1, aatom%Ntr

          kc = aatom%at(kr)%ik

          i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
          ! 				chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
          ! 				neff = aatom%stage(j) * sqrt(aatom%Rydberg / (aatom%E(j) - aatom%E(i)) )


          Nblue = aatom%continua(kc)%Nb; Nred = aatom%continua(kc)%Nr
          Nl = Nred-Nblue + 1

          JJ = 0.0
          JJb = 0.0

          icell_d = 1
          if (ldissolve) then
             if (aatom%ID=="H") icell_d = icell
          endif

          do l=2, Nl
             wl = ( lambda(Nblue+l-1) - lambda(Nblue+l-2) )

             a1 = aatom%continua(kc)%alpha_nu(l,icell_d)
             a2 = aatom%continua(kc)%alpha_nu(l-1,icell_d)

             J1 = Ieff(Nblue+l-1) * a1 / lambda(Nblue+l-1)
             J2 = Ieff(Nblue+l-2) * a2 / lambda(Nblue+l-2)

             JJ = JJ + 0.5 * (J1 + J2) * wl

             ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-1))
             twohnu3_c2 = twohc/lambda(Nblue+l-1)**3
             J1 = ( Ieff(Nblue+l-1) + twohnu3_c2 ) * ehnukt * a1 / lambda(Nblue+l-1)

             ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-2))
             twohnu3_c2 = twohc/lambda(Nblue+l-2)**3
             J2 = ( Ieff(Nblue+l-2) + twohnu3_c2 ) * ehnukt * a2 / lambda(Nblue+l-2)

             JJb = JJb + 0.5 * (J1 + J2) * wl

          enddo

          aatom%continua(kc)%Rij(id) = aatom%continua(kc)%Rij(id) + fourpi_h * ( JJ / n_rayons )
          aatom%continua(kc)%Rji(id) = aatom%continua(kc)%Rji(id) + fourpi_h * ( JJb / n_rayons) * &
               ( aatom%nstar(i,icell)/aatom%nstar(j,icell) )


       end do atrc_loop

       aatom => NULL()

    end do aatom_loop


    return
  end subroutine calc_rates_draft

  !occupa prob
  !Define a negative Tex for inversion of Temperature ? because the ratio of ni/nj exists
  !even with inversion of pops.
  !If negative Tex, T = Tmax, just like if ni = nj*gij
  !Opacity is then 0 and only emissivity contributes.
!!!!!!! Can be zero if the conditions not met and lead to infinity !!!
!!!!!
  SUBROUTINE calc_delta_Tex_atom(icell, atom, dT, Tex, Tion, write_neg_Tex)
    integer, intent(in) :: icell
    type(AtomType), intent(inout) :: atom
    logical, intent(in) :: write_neg_tex
    real(kind=dp), intent(out) :: dT, Tex, Tion
    integer :: nact, kr, kc, i, j
    real(kind=dp) :: deltaE_k, ratio, Tdag, gij, dT_loc, Tion_loc, Tex_loc
    real(kind=dp) :: wi, wj, ni_on_nj_star


    dT = 0.0_dp !for all transitions of one atom
    Tex = 0.0
    Tion = 0.0

    tr_loop : do kr=1, atom%Ntr
       kc = atom%at(kr)%ik

       SELECT CASE (atom%at(kr)%trtype)

       CASE ("ATOMIC_LINE")

          i = atom%lines(kc)%i; j = atom%lines(kc)%j
          wj = 1.0; wi = 1.0
          ! 				if (ldissolve) then
          ! 					if (atom%ID=="H") then
          ! 												!nn
          ! 						wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1))
          ! 						wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1))
          !
          ! 					endif
          ! 				endif

          !Same goes for the continuum:
          !I check populations inversions, and do not compute temperature in that case
          !or if populations of one level is below a threshold I skip two.
          !The latter does not prevent negative Tex and Tion
          !The following check is similar to the one for computing opacities

          !condition on values of individual populations ? (wrt a threshold)
          if ((atom%n(i,icell) - atom%lines(kc)%gij * atom%n(j,icell)) <= 0.0 ) cycle tr_loop

          ! 				if ( (atom%n(i,icell)/ntotal_atom(icell,atom) < frac_limit_pops) .or. (atom%n(j,icell)/ntotal_atom(icell,atom) < frac_limit_pops) ) cycle tr_loop


          Tdag = atom%lines(kc)%Tex(icell)
          deltaE_k = (atom%E(j)-atom%E(i)) / KBOLTZMANN

          ratio = dlog(wi * atom%n(j,icell) * atom%lines(kc)%gij / (atom%n(i,icell)*wj))

          !write(*,*) "line"
          !write(*,*) "nstar:", atom%nstar(i,icell), atom%nstar(j,icell)
          !write(*,*) "n:", atom%n(i,icell), atom%n(j,icell)
          !write(*,*) "T:", T(icell), Tdag, -deltaE_k/ratio

          if (ratio < 0.0_dp) then

             atom%lines(kc)%Tex(icell) = -deltaE_k / ratio
             Tex_loc = atom%lines(kc)%Tex(icell)

          else

             ! 					atom%lines(kc)%Tex(icell) = Tmax
             ! 					Tex_loc = atom%lines(kc)%Tex(icell)
             if (ratio == 0.0) then
                atom%lines(kc)%Tex(icell) = Tmax
                Tex_loc = atom%lines(kc)%Tex(icell)
             else !ratio > 0, T < 0
                atom%lines(kc)%Tex(icell) = -deltaE_k / ratio
                Tex_loc = atom%lines(kc)%Tex(icell)
             endif

             ! 					if (write_neg_tex) then
             ! 					write(unit_invfile,*) "-------------------------------------------------------------------------------"
             ! 					write(unit_invfile,"('icell = '(1I9), ' atom '(1A2), ' T='(1F12.5), ' lam='(1F12.5)' nm' )") icell, atom%ID, T(icell), atom%lines(kc)%lambda0
             ! 					write(unit_invfile, "(' -> line '(1I2)' -> '(1I2), ' d(ni-njxgij)/ni='(1ES20.7E3) )") i, j, (atom%n(i,icell)*wi/wj - atom%n(j,icell) * gij)/atom%n(j,icell)
             ! 					write(unit_invfile,"(' Tex='(1ES20.7E3) )") Tex_loc
             ! 					write(unit_invfile, "('log(ratio) = '(1ES20.7E3), ' ratio = '(1ES20.7E3) )") ratio, exp(ratio)
             ! 					write(unit_invfile,"( 'w(i)='(1ES14.7), ' w(j)='(1ES20.7E3) )") wi, wj
             ! 					write(unit_invfile,"( 'n(i)='(1ES14.7), ' n(j)='(1ES20.7E3), ' n(j)xgij='(1ES20.7E3) )") atom%n(i,icell), atom%n(j,icell), atom%n(j,icell) * gij
             ! 					write(unit_invfile,"( 'n*(i)='(1ES20.7E3), ' n*(j)='(1ES20.7E3) )") atom%nstar(i,icell), atom%nstar(j,icell)
             ! 					write(unit_invfile,"( 'ne='(1ES20.7E3), ' ntotal='(1ES20.7E3) )") ne(icell), ntotal_atom(icell,atom)
             ! 					write(unit_invfile,*) "-------------------------------------------------------------------------------"
             ! 					endif
          endif
          !write(*,*) "loc, line:", T(icell), Tex_loc

          dT_loc = abs(Tdag-Tex_loc)/(tiny_dp + Tex_loc)
          !        			dT_loc = abs(Tdag-Tex_loc)/(tiny_dp + Tdag)
          dT = max(dT, dT_loc)
          Tex = max(Tex, Tex_loc)


       CASE ("ATOMIC_CONTINUUM")

          i = atom%continua(kc)%i; j = atom%continua(kc)%j
          wj = 1.0; wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(j)),hydrogen%n(1,icell))

             endif
          endif

          !condition on values of individual populations ? (wrt a threshold)


          Tdag = atom%continua(kc)%Tex(icell)

          deltaE_k = (atom%E(j)-atom%E(i)) / KBOLTZMANN

          !at threshold
          !i.e., deltaE is hnu0

          ! 				ni_on_nj_star = ne(icell) * phi_T(icell, atom%g(i)/atom%g(j), atom%E(j)-atom%E(i))
          ni_on_nj_star = atom%nstar(i,icell)/atom%nstar(j,icell)


          gij = ni_on_nj_star * exp(-hc_k/atom%continua(kc)%lambda0/T(icell))
          ratio = log( wj*atom%n(i,icell)  / ( wi * atom%n(j,icell) * gij ) )

          if ((atom%n(i,icell) - gij * atom%n(j,icell)) <= 0.0 ) cycle tr_loop
          ! 				if ( (atom%n(i,icell)/ntotal_atom(icell,atom) < frac_limit_pops) .or. (atom%n(j,icell)/ntotal_atom(icell,atom) < frac_limit_pops) ) cycle tr_loop


          !write(*,*) "cont"
          !write(*,*) "nstar:", atom%nstar(i,icell), atom%nstar(j,icell)
          !write(*,*) "n:", atom%n(i,icell), atom%n(j,icell)
          !write(*,*) "T:", T(icell), Tdag, deltaE_k/ratio

          if (ratio > 0.0_dp) then

             !ionisation temperature
             atom%continua(kc)%Tex(icell) = deltaE_k / ratio
             Tion_loc = atom%continua(kc)%Tex(icell)

          else

             ! 					atom%continua(kc)%Tex(icell) = Tmax
             ! 					Tion_loc = atom%continua(kc)%Tex(icell)
             if (ratio == 0.0) then
                atom%continua(kc)%Tex(icell) = Tmax
                Tion_loc = atom%continua(kc)%Tex(icell)
             else
                atom%continua(kc)%Tex(icell) = deltaE_k / ratio
                Tion_loc = atom%continua(kc)%Tex(icell)
             endif

             ! 					if (write_neg_tex) then
             ! 					write(unit_invfile,*) "-------------------------------------------------------------------------------"
             ! 					write(unit_invfile,"('icell = '(1I9), ' atom '(1A2), ' T='(1F12.5), ' lam='(1F12.5)' nm' ) ") icell, atom%ID, T(icell), atom%continua(kc)%lambda0
             ! 					write(unit_invfile, "(' -> cont '(1I2)' -> '(1I2), ' d(ni-njxgij)/ni='(1ES20.7E3) )") i, j, (atom%n(i,icell)*wi/wj - atom%n(j,icell) * gij)/atom%n(j,icell)
             ! 					write(unit_invfile,"(' Tion='(1ES20.7E3) )") Tion_loc
             ! 					write(unit_invfile, "('log(ratio) = '(1ES20.7E3), ' ratio = '(1ES20.7E3) )") ratio, exp(ratio)
             ! 					write(unit_invfile,"( 'w(i)='(1ES14.7), ' w(j)='(1ES20.7E3) )") wi, wj
             ! 					write(unit_invfile,"( 'n(i)='(1ES14.7), ' n(j)='(1ES20.7E3), ' n(j)xgij='(1ES20.7E3) )") atom%n(i,icell), atom%n(j,icell), atom%n(j,icell) * gij
             ! 					write(unit_invfile,"( 'n*(i)='(1ES20.7E3), ' n*(j)='(1ES20.7E3) )") atom%nstar(i,icell), atom%nstar(j,icell)
             ! 					write(unit_invfile,"( 'ne='(1ES20.7E3), ' ntotal='(1ES20.7E3) )") ne(icell), ntotal_atom(icell,atom)
             ! 					write(unit_invfile,*) "-------------------------------------------------------------------------------"
             ! 					endif
          endif
          !write(*,*) "loc, cont:", T(icell), Tion_loc

          dT_loc = abs(Tdag-Tion_loc)/(Tion_loc + tiny_dp)
          !        			dT_loc = abs(Tdag-Tion_loc)/(Tdag + tiny_dp)
          dT = max(dT, dT_loc)
          Tion = max(Tion_loc, Tion)


       CASE DEFAULT

          CALL Error("Unkown transition type", atom%at(kr)%trtype)

       END SELECT

    end do tr_loop

    RETURN
  END SUBROUTINE calc_delta_Tex_atom


  SUBROUTINE calc_Tex(icell) !for all atoms
    integer, intent(in) :: icell
    type(AtomType), pointer :: atom
    integer :: nact
    real(kind=dp) :: dT, Tex, Tion

    do nact=1, NactiveAtoms
       atom => ActiveAtoms(nact)%ptr_atom

       CALL calc_delta_Tex_atom(icell, atom, dT, Tex, Tion, .false.)

       atom => NULL()
    enddo

    RETURN
  END SUBROUTINE calc_Tex



  SUBROUTINE calc_rate_matrix(id, icell, switch_lte)
    integer, intent(in) :: id, icell
    logical, intent(in) :: switch_lte
    integer :: nact

    do nact=1, NactiveAtoms
       call rate_matrix_atom(id, icell, ActiveAtoms(nact)%ptr_atom, switch_lte)
    enddo

    RETURN
  END SUBROUTINE calc_rate_matrix

  subroutine print_gamma(id, atom)
    integer, intent(in) :: id
    type (AtomType), intent(in) :: atom
    integer :: l, lp

    do l=1, atom%Nlevel
       write(*, '(1I3":", *(ES14.5E3))') l, (atom%Gamma(l,lp,id), lp=1, atom%Nlevel)
    enddo

    return
  end subroutine print_gamma

  subroutine eliminate_delta(id, atom)
    !see Hubeny & Mihalas 2014 eq. 14.8a
    !calc the term in delta(l,l') (diagonal elements of rate matrix)
    !delta(l,l') = Gamma(l,l) = sum_l' R_ll' + C_ll' = sum_l' -Gamma(l',l)
    integer, intent(in) :: id
    type (AtomType), intent(inout) :: atom
    integer :: l

    do l = 1, atom%Nlevel
       atom%Gamma(l,l,id) = 0.0_dp
       !Gamma(j,i) = -(Rij + Cij); Gamma(i,j) = -(Rji + Cji)
       !Gamma(i,i) = Rij + Cij = -sum(Gamma(j,i))
       atom%Gamma(l,l,id) = sum(-atom%Gamma(:,l,id)) !positive
    end do


    return
  end subroutine eliminate_delta


  subroutine check_elements(id, icell, N,A,printed_out_message)
    !check that digonal of A is >= 0 and the off diagonal are < 0
    real :: sign = -1.0 !save attribute but unchanged
    character(len=*), intent(in) :: printed_out_message
    integer :: i,j, l, lp
    integer, intent(in) :: N, id, icell
    real(kind=dp), intent(in) :: A(N,N)
    real(kind=dp) :: diag! = 0.0_dp !save attribute if initialized and changed !

    diag = 0.0_dp

    do i=1,N
       diag = diag + A(i,i)
    enddo


    do i=1,N
       do j=1,N
          if (i==j) then
             if (A(i,i) < 0.0) then
                write(*,*) printed_out_message
                write(*,*) "id=",id," icell=",icell
                do l=1, N
                   write(*,*) l, A(l,l)
                enddo
                call error("The diagonal of the rate matrix should be positive!")
             endif
          else
             if (A(i,j) > 0.0) then
                write(*,*) printed_out_message
                write(*,*) "id=",id," icell=",icell
                do l=1, N
                   write(*,*) l
                   do lp=1,N
                      if (lp /= l) write(*,'(2ES17.8E3, 1X)', advance="no") A(l,lp)
                   enddo
                enddo
                call error("off-diagonal elements of the rate matrix should be negative")
             endif
          endif
       enddo
    enddo

    return
  end subroutine check_elements


  !occupa prob
  SUBROUTINE rate_matrix_atom(id, icell, atom, switch_to_lte)
    !see Hubeny & Mihalas 2014 eq. 14.8a to 14.8c for the rate matrix elements
    !cannot directly remove diagonal here because rate matrix is already filled with
    ! collision rates. Needs to add collision rates here otherwise.
    integer, intent(in) :: icell, id
    type (AtomType), intent(inout) :: atom
    logical, optional, intent(in) :: switch_to_lte
    integer :: kr, kc, i, j, l, Nblue, Nred, Nl
    real(kind=dp) :: wj, wi, neff

    if (present(switch_to_lte)) then
       !because Gamma = C
       if (switch_to_lte) return
    end if

    tr_loop : do kr=1, atom%Ntr
       kc = atom%at(kr)%ik

       SELECT CASE (atom%at(kr)%trtype)

       CASE ("ATOMIC_LINE")
          i = atom%lines(kc)%i; j = atom%lines(kc)%j
          wj = 1.0
          wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1),hydrogen%n(1,icell))
                wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1),hydrogen%n(1,icell))
                ! 					else
                !
                ! 						!neff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(find_continuum(atom,j)) - atom%E(i)) )

             endif
          endif

          atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - atom%lines(kc)%Rij(id) * wj/wi
          atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - atom%lines(kc)%Rji(id)



       CASE ("ATOMIC_CONTINUUM")
          i = atom%continua(kc)%i; j = atom%continua(kc)%j
          wj = 1.0
          wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1),hydrogen%n(1,icell))
                wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1),hydrogen%n(1,icell))
                ! 					else
                !
                ! 						!neff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(find_continuum(atom,j)) - atom%E(i)) )

             endif
          endif

          atom%Gamma(j,i,id) = atom%Gamma(j,i,id) - atom%continua(kc)%Rij(id) * wj/wi
          atom%Gamma(i,j,id) = atom%Gamma(i,j,id) - atom%continua(kc)%Rji(id)

       CASE DEFAULT

          CALL Error("Unkown transition type", atom%at(kr)%trtype)

       END SELECT

    end do tr_loop


    RETURN
  END SUBROUTINE rate_matrix_atom

  !to do occupa prob
  SUBROUTINE radrate_matrix_atom(id, icell, atom, gam_rad, dgam_raddne)

    !Rji continua are divided by ne !

    integer, intent(in) :: icell, id
    type (AtomType), intent(in) :: atom
    real(kind=dp), intent(out), dimension(atom%Nlevel, atom%Nlevel) :: gam_rad, dgam_raddne
    integer :: kr, kc, i, j

    gam_rad(:,:) = 0.0_dp
    dgam_raddne(:,:) = 0.0_dp !derivative of recombinaison rate matrix to dne

    if (lforce_lte) return

    tr_loop : do kr=1, atom%Ntr
       kc = atom%at(kr)%ik

       SELECT CASE (atom%at(kr)%trtype)

       CASE ("ATOMIC_LINE")
          i = atom%lines(kc)%i; j = atom%lines(kc)%j

          gam_rad(j,i) = gam_rad(j,i) - atom%lines(kc)%Rij(id)
          gam_rad(i,j) = gam_rad(i,j) - atom%lines(kc)%Rji(id)



       CASE ("ATOMIC_CONTINUUM")
          i = atom%continua(kc)%i; j = atom%continua(kc)%j


          gam_rad(j,i) = gam_rad(j,i) - atom%continua(kc)%Rij(id)
          gam_rad(i,j) = gam_rad(i,j) - atom%continua(kc)%Rji(id) / ne(icell)
          dgam_raddne(i,j) = dgam_raddne(i,j) - atom%continua(kc)%Rji(id) / ne(icell)

       CASE DEFAULT

          CALL Error("Unkown transition type", atom%at(kr)%trtype)

       END SELECT

    end do tr_loop

    do i=1, atom%Nlevel
       dgam_raddne(i,i) = -sum(dgam_raddne(:,i))
    enddo


    RETURN
  END SUBROUTINE radrate_matrix_atom

  subroutine partial_radrate_matrix_partial_ne_atom(id,icell,atom, dgamma_dne)
    integer, intent(in) :: id, icell
    type (atomtype), intent(in) :: atom
    real(kind=dp), dimension(atom%Nlevel,atom%Nlevel), intent(out) :: dgamma_dne

    integer :: kr, kc, i, j

    dgamma_dne(:,:) = 0.0_dp
    if (lforce_lte) return

    tr_loop : do kr=atom%Ntr_line+1, atom%Ntr
       kc = atom%at(kr)%ik

       SELECT CASE (atom%at(kr)%trtype)

       CASE ("ATOMIC_LINE")
          call error("No bound-bound transition allowed for partial derivative of Gamma_rad to ne!")
       CASE ("ATOMIC_CONTINUUM")
          i = atom%continua(kc)%i; j = atom%continua(kc)%j

          dgamma_dne(i,j) = dgamma_dne(i,j) - atom%continua(kc)%Rji(id) / ne(icell) !current estimate of ne

       CASE DEFAULT

          CALL Error("Unkown transition type", atom%at(kr)%trtype)

       END SELECT

    end do tr_loop

    do i=1,atom%Nlevel
       dgamma_dne(i,i) = -sum(dgamma_dne(:,i))
    enddo

    return
  end subroutine partial_radrate_matrix_partial_ne_atom

  subroutine collrate_matrix_atom (id, icell, atom, gamma_col, deriv_rate, dgamma_dne)

    !see Hubeny & Mihalas 2014 eq. 14.8b

    integer, intent(in) :: id, icell
    type (AtomType), intent(in) :: atom!, pointer
    real(kind=dp), intent(in), dimension(atom%Nlevel,atom%Nlevel) :: deriv_rate
    real(kind=dp), intent(out), dimension(atom%Nlevel,atom%Nlevel) :: gamma_col, dgamma_dne
    integer :: l,lp

    !see initgamma_atom
    gamma_col(:,:) = 0.0_dp
    dgamma_dne(:,:) = 0.0_dp
    !I need the temporary deriv_Rate juste like atom%C because
    !otherwise dgamma_dne is set to zero dans not properly computed.
    !alternative dgamma_dne = transpose(-dgamma_dne) but invoke internal loop.

    do l=1,atom%Nlevel-1
       do lp=l+1,atom%Nlevel
          gamma_col(l,lp) = -atom%C(lp,l,id)
          gamma_col(lp,l) = -atom%C(l,lp,id)
          dgamma_dne(l,lp) = -deriv_rate(lp,l)!dgamma_dne(lp,l)!-2*atom%C(lp,l,id)/ne(icell)
          dgamma_dne(lp,l) = -deriv_rate(l,lp)!dgamma_dne(l,lp)!-atom%C(l,lp,id)/ne(icell)
       enddo
    enddo

    do l=1,atom%Nlevel
       dgamma_dne(l,l) = -sum(dgamma_dne(:,l))
    enddo

    return
  end subroutine collrate_matrix_atom


  subroutine ionisation_frac_lte(elem, k, ne, fjk, dfjk)
    real(kind=dp), intent(in) :: ne
    integer, intent(in) :: k
    type (Element), intent(in) :: Elem
    real(kind=dp), dimension(:), intent(inout) :: fjk, dfjk
    real(kind=dp) :: Uk, Ukp1, sum1, sum2
    integer :: j,  Nstage
    !return Nj / Ntot
    !for neutral j * Nj/Ntot = 0

    !fjk(1) is N(1) / ntot !N(1) = sum_j=1 sum_i=1^N(j) n(i)
    fjk(1)=1.
    dfjk(1)=0.
    sum1 = 1.
    sum2 = 0.
    Uk = getPartitionFunctionk(elem,1,k)
    do j=2,Elem%Nstage
       Ukp1 = getPartitionFunctionk(elem,j,k)
       fjk(j) = Sahaeq(k,fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)
       !fjk(j) = ne * cste, der = cste  fjk(j) / ne

       dfjk(j) = -(j-1)*fjk(j)/ne

       sum1 = sum1 + fjk(j)
       sum2 = sum2 + dfjk(j)

       Uk = Ukp1
    end do

    fjk(:)=fjk(:)/sum1
    dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1


    return
  end subroutine ionisation_frac_lte

  !should replace an equation for H SEE and He SEE
  !-> do we need the same equation for Ntotal = H + He
  !or one equation per SEE is to be replace by notal_H for H
  !and notal_He for He ?
  !-> why 1 - sum(pops)/ntot works and not to set 0, while it should be equivalent ?
  subroutine particle_conservation_H_and_He (icell, Neq, x, f, df)
    integer, intent(in) :: icell, Neq
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(inout) :: f(neq), df(neq, neq)
    integer :: imaxpop, nlev, l, lp
    real(kind=dp) :: ntotal

    ntotal = ntotal_atom(icell,hydrogen)

    nlev = Neq - 1 !remove electrons

    imaxpop = locate(x(1:hydrogen%Nlevel), maxval(x(1:hydrogen%Nlevel)))

    !mass conservation: ntotal = sum_i n(i)
    !eq imaxpop is replaced by:
    !f(imaxpop) = ntotal - sum_i n(i) = 0
    !or
    !f(imaxpop) = 1.0 - sum_i n(i) / ntotal = 0

    f(imaxpop) = 1.0_dp
    f(imaxpop) = f(imaxpop) -sum(x(1:hydrogen%Nlevel)) / ntotal !~0, better to set to 0?

    ! 		df(imaxpop,neq) = 0.0_dp
    df(imaxpop,:) = 0.0_dp
    df(imaxpop,1:hydrogen%Nlevel) = -1.0_dp / ntotal

    !-> H+He ?
    ! 		df(imaxpop,1:neq-1) = -1.0_dp / (ntotal+ntotal_atom(icell, helium))

    if (helium_is_active) then

       !relative to the sub array !!
       imaxpop = locate(x(hydrogen%Nlevel+1:nlev), maxval(x(hydrogen%Nlevel+1:nlev)))
       imaxpop = hydrogen%Nlevel+imaxpop
       ! 			write(*,*) "imaxpop = ", imaxpop
       f(imaxpop) = 1.0_dp
       f(imaxpop) = f(imaxpop) -sum(x(hydrogen%Nlevel+1:nlev)) / ntotal_atom(icell,helium) !~0 better to set 0 ?

       df(imaxpop,:) = 0.0_dp
       df(imaxpop,hydrogen%Nlevel+1:nlev) = -1.0_dp / ntotal_atom(icell,helium)

       ! 			write(*,*) ntotal_atom(icell,helium), sum(x(hydrogen%Nlevel+1:nlev)), f(imaxpop)
       ! 			write(*,*) f(hydrogen%Nlevel+1:neq-1)
       ! 			stop

    endif
    ! 		write(*,*) "H, He=", -1/ntotal_atom(icell, hydrogen), -1.0/ntotal_atom(icell, helium)
    ! 		write(*,*) f(1:neq-1)
    ! 		do l=1, neq-1
    ! 			write(*,*) l, (df(l,lp), lp=1,neq-1)
    ! 		enddo
    ! 		stop

    return
  end subroutine particle_conservation_H_and_He

  !F_cc = 1.0 - 1/ne * (np + nHeII + 2*nHeIII) = 0
  subroutine non_lte_charge_conservation_H_and_He (icell, neq, x, f, df)
    use atom_type, only : find_continuum
    !also computes the derivative of the charge conservation equation wuth
    !ne and total population of stages with stage > 0 (like nH(Nlevel)=nHII)
    integer, intent(in) :: icell, neq
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
    integer :: n, j, n_start, j0, j1
    type (element), pointer :: elem
    real(kind=dp) :: akj, np, nheII, nheIII
    real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk
    !Separate in two parts: 1 non-LTE atoms and 2 electrons from LTE elements

    !first identify proton levels (hydrogen%Nlevel)
    !heII levels and heIII levels

    !derivative of CC (charge conservation) to ne
    !x = (n1,n2,n3 .... ne)
    !Ionisation stage of nHII = np is 1, 1 for nHeII and 2 for nHeIII
    np = hydrogen%n(hydrogen%Nlevel,icell)
    if (helium_is_active) then
       !I'm assuming I'm starting from HeI
       j0 = find_continuum(helium, 1) !ground state of HeII
       j1 = find_continuum(helium, j0) !heIII
       nHeII = sum(helium%stage(j0:j1-1)*helium%n(j0:j1-1,icell))
       nHeIII = helium%stage(j1) * helium%n(j1,icell)
       ! 			write(*,*) " nHeIII=", nHeIII, 2 * helium%n(helium%Nlevel,icell)
       ! 			write(*,*) nHeII + nheIII/2. + sum(helium%n(1:j0-1,icell)), ntotal_atom(icell, helium)
       ! 			write(*,*) sum(helium%stage(:) * helium%n(:,icell)), nHeII + nHeIII
       !->  check that match the atomic file
       ! 			write(*,*) "HeII ground level index (fortran)", j0, " he III", j1
       ! 			write(*,*) size(helium%stage(1:j0-1)), size(helium%stage(j0:j1-1)), 1
       n_start = 3
       df(neq,neq) = 1/x(neq)**2 * (np + nheII + nheIII) !the factor stage is included!!
       !the sum should return 0 * nHeI + 1 * nHeII + 2 * nHeIII
       !df(neq,neq) = 1/x(neq)**2 * (np + sum(helium%stage(:) * helium%n(:,icell)))
       f(neq) = 1.0 - (1.0/x(neq)) * (np + nheII + nheIII) !factor 2 in nheIII included
       !derivative to protons
       df(neq,hydrogen%Nlevel) = -1.0 / x(neq)
       !derivative to nheIII
       df(neq, neq-1) = -2.0 / x(neq)
       !derivative to nHeII!
       !or derivative to nHeII_i individual levels
       do j=1, j1-j0
          df(neq, neq-1-j) = -1/x(neq) !or for each level of nHeII I have the same derivative?
       enddo
       ! 			stop
    else
       n_start = 2
       df(neq,neq) = 1/x(neq)**2 * np
       f(neq) = 1.0 - (1.0/x(neq)) * np
       !derivative to protons
       df(neq,hydrogen%Nlevel) = -1.0 / x(neq)
    endif

    !check:
    !if (hydrogen%n(hydrogen%Nlevel,icell) /= sum(hydrogen%n(:,icell)*hydrogen%stage(:)))
    !to do


    !now contribution from LTE atoms.
    do n=n_start, 26 !for 26 elements
       !should avoid H and He if present !
       elem => Elements(n)%ptr_elem

       call ionisation_frac_lte(elem, icell, x(neq), fjk, dfjk)

       do j=2, elem%Nstage
          !pure LTE term j = 1 corresponds to the first ionisation stage always
          !unlike in solve_ne where this loop is also for non-LTE models
          akj = (j-1) * elem%abund * nhtot(icell) !(j-1) = 0 for neutrals and 1 for singly ionised
          f(neq) = f(neq) - akj * fjk(j) / x(neq)
          df(neq,neq) = df(neq,neq) + akj * fjk(j) / x(neq)**2 - akj / x(neq) * dfjk(j)
       end do

       elem => null()

    enddo

    ! 	write(*,*) "df(neq,neq)=", df(neq, neq)
    ! 	write(*,*) "df(:,neq)=", df(:,neq)
    ! 	write(*,*) "df(neq,:)=",df(neq,:)
    ! 	write(*,*) "total=",df


    return
  end subroutine non_lte_charge_conservation_H_and_He

  subroutine rate_equations_H_and_He(id,icell, neq, ihel, dgrdne, dgcdne, x, f, df)
    integer, intent(in) :: id, icell, neq, ihel
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(in), dimension(neq-1,neq-1) :: dgrdne, dgcdne
    real(kind=dp), intent(out) ::  f(neq), df(neq,neq)
    integer :: i, j, l, lp, Nlevels
    real(kind=dp) :: part1, part2
    !derivative of SEE to ne
    !dCij/dne = Cij/ne
    !dCji/dne = 2*Cji/ne
    !dRji_cont/dne = Rji_cont/ne

    !Have to do better, to much avoidable loops!!

    Nlevels = neq - 1 !should be hydrogen%Nlevel + helium%Nlevel

    !f and df 0 where helium or H levels stops

    !equations and derivative for H
    f(1:hydrogen%Nlevel) = matmul(hydrogen%Gamma(:,:,id), x(1:hydrogen%Nlevel))
    !derivative to levels n
    df(1:hydrogen%Nlevel,1:hydrogen%Nlevel) = hydrogen%gamma(:,:,id)
    !derivative to ne
    df(1:hydrogen%Nlevel,neq) = matmul( dgcdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel), x(1:hydrogen%Nlevel) )

    ! do i=1, hydrogen%Nlevel
    ! 	write(*,*) i, " eq deriv:", (df(i, j), j=1, neq)
    ! enddo
    if (helium_is_active) then
       !He
       f(ihel:Nlevels) = matmul(helium%gamma(:,:,id), x(ihel:Nlevels))
       df(ihel:Nlevels,ihel:Nlevels) = helium%Gamma(:,:,id)
       df(ihel:Nlevels,neq) = matmul (dgcdne(ihel:Nlevels,ihel:Nlevels), x(ihel:Nlevels))
    endif

    ! do i=1, helium%Nlevel
    ! 	write(*,*) i, " eq deriv:", (df(i-1+ihel, j), j=1, neq)
    ! enddo
    ! stop
    if (lforce_lte) return

    !radiative part
    !H
    df(1:hydrogen%Nlevel,neq) = df(1:hydrogen%Nlevel,neq) + matmul( dgrdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel), &
         x(1:hydrogen%Nlevel) )
    if (helium_is_active) then
       !He
       df(ihel:Nlevels,neq) = df(ihel:Nlevels,neq) + matmul(dgrdne(ihel:Nlevels,ihel:Nlevels), x(ihel:Nlevels))
    endif

    return

    return
  end subroutine rate_equations_H_and_He

  !H active if He active at the moment !
  subroutine see_ionisation_nonlte_H_and_He(id, icell, dM, dne)
    !TO DO :positivity tests + zero value below prec_pops
    !Stack all equations and derivative in the same matrix

    integer, intent(in) :: id, icell
    real(kind=dp), intent(out) :: dM, dne
    integer :: nact, Neq, nact_start, n_iter, i, kr
    !smoothing parameter default if damp_char. If native pops, damp the iterations by a facteur damp_scaling larger
    real(kind=dp), parameter :: damp_scaling = 10.0, damp_char = 5.0, precision = 1d-5
    real(kind=dp) :: d_damp, nedag, delta_f, dfpop, dfne, dM_he, dfpop_he, dfpop_h
    real(kind=dp), allocatable :: gamma_r(:,:,:), ndag(:), RR(:,:), CC(:,:), deriv(:,:), dgamrdne(:,:), dgamcdne(:,:), gamma_tot(:,:)
    real(kind=dp), allocatable :: fvar(:), dfvar(:,:), xvar(:) !Jacobian
    logical :: lconverged, neg_pops, rest_damping, lsecond_try, verbose
    integer, parameter :: maxIter = 1000
    integer :: max_iter, Nlevel_atoms, ihel=0, ih=0, n_active

    verbose = .false. !debug mode
    lconverged = .false.
    rest_damping = .false.
    d_damp = damp_char
    lsecond_try = .false. !set to .true. to avoid starting again if the first iterative scheme failed


    !Hydrogen always active here.
    !number of equations = Nlevel_atom1 + Nlevel_atom2 + ... + 1 for ne
    Nlevel_atoms = hydrogen%Nlevel
    n_active = 1
    ih = 1
    if (helium_is_active) then
       ihel = Nlevel_atoms + 1    !position of SEE of helium among
       !all levels
       Nlevel_atoms = Nlevel_atoms + helium%Nlevel !sum not max
       !SEE(H) -> ih : hydrogen%Nlevel
       !SEE(He) -> ihel : ihel + helium%Nlevel
       n_active = n_active + 1
    endif
    Neq = Nlevel_atoms + 1
    !currently allocate and deallocate for each cell
    !futur fixed size
    !unknowns = populations of atomic levels + electron density
    allocate(xvar(neq)); xvar(:) = 0.0_dp
    !equations = Nlevel SEE per atom  + charge conservation (CC)
    allocate(fvar(neq)); fvar(:) = 0.0_dp
    !Jacobian of the systems = derivative of i equations to xvar of sum_j partial(fvar(i))/partial(xvar(j))
    allocate(dfvar(neq,neq)); dfvar(:,:) = 0.0_dp
    !old ne values to rest for the other cell. The new value is stored in ne_new for the next iterations
    nedag = ne(icell)

    allocate(ndag(Nlevel_atoms)) !hydrogen%Nlevel if only H
    ndag(1:hydrogen%Nlevel) = hydrogen%n(:,icell)
    if (helium_is_active) then
       ndag(ihel:Nlevel_atoms) = helium%n(:,icell)
    endif

    !temporary rate deriv matrix
    allocate(deriv(Nlevel_Atoms, Nlevel_atoms))


    !local radiative rates and collisional rates matrix to build total rate matrix
    !for current estimate of ne.
    allocate(RR(Nlevel_atoms,Nlevel_atoms), CC(Nlevel_atoms,Nlevel_atoms), dgamcdne(Nlevel_atoms,Nlevel_atoms))
    !ne dependent radiative rate matrix for each atom + derivative to ne of rate matrix for all equations (of all atoms)
    allocate(gamma_r(n_active, Nlevel_atoms,Nlevel_atoms),dgamrdne(Nlevel_atoms,Nlevel_atoms))
    dgamrdne = 0.0_dp; dgamcdne = 0.0_dp
    !!allocate(gamma_tot(Nlevel_atoms, Nlevel_atoms)); gamma_tot = 0.0_dp
    call radrate_matrix_atom(id, icell, hydrogen, gamma_r(1,1:hydrogen%Nlevel,1:hydrogen%Nlevel), &
         deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel))
    dgamrdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel) = deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel)
    if (helium_is_active) then
       deriv = 0.0_dp
       call radrate_matrix_atom(id, icell, helium, gamma_r(n_active,1:helium%Nlevel,1:helium%Nlevel), &
            deriv(1:helium%Nlevel,1:helium%Nlevel))
       dgamrdne(ihel:Nlevel_atoms,ihel:Nlevel_atoms) = deriv(1:helium%Nlevel,1:helium%Nlevel)
    endif
    !gamma_r(i,j) = -Rji, diagonal not included!
    !derivative of gamma_r wrt ne is included (with the diagonal!!). It is only non-zero
    !for recombinason rate and is independent of ne as it is Rji/ne (Rji = Rji/ne_old * ne_new at each iteration so derivative to ne_new is always Rji/ne_old)
    ! 		if (lforce_lte) then
    ! 			gamma_r(:,:,:) = 0.0_dp
    ! 			dgamrdne(:,:) = 0.0_dp
    ! 		endif

    !start Newton-Raphson iterations
    if (verbose) write(*,*) "T = ", T(icell), " nHtot = ", nHtot(icell)

    n_iter = 0
    max_iter = maxIter

    iterative_loop : do while (.not.lconverged)
       delta_f = 0.0_dp
       dfpop = 0.0_dp !only populations
       dfne = 0.0_dp
       dfpop_he = 0.0_dp
       dfpop_h = 0.0_dp

       n_iter = n_iter + 1

       !always the last one
       xvar(Neq) = ne(icell)
       fvar(:) = 0.0
       dfvar(:,:) = 0.0

       !set unknowns to old values
       !First hydrogen
       xvar(1:hydrogen%Nlevel) = hydrogen%n(:,icell)
       RR(:,:) = gamma_r(1,:,:) !init to radiative rates
       !we multiply ne by zero if lforce_lte.
       do kr=1,hydrogen%Ncont
          RR(hydrogen%continua(kr)%i,hydrogen%continua(kr)%j) = RR(hydrogen%continua(kr)%i,hydrogen%continua(kr)%j) * ne(icell)
       enddo
       call collision_matrix_atom(id, icell, hydrogen, deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel))
       ! 			!fill atom%C and its derivative
       !derivative independent of the new value of ne but Collision need to be updated.
       !build collisional rate matrix and its equivalent derivative
       call collrate_matrix_atom(id, icell, hydrogen, CC(1:hydrogen%Nlevel,1:hydrogen%Nlevel), &
            deriv(1:hydrogen%Nlevel,1:hydrogen%Nlevel), dgamcdne(1:hydrogen%Nlevel,1:hydrogen%Nlevel))

       ! 			gamma_tot(1:hydrogen%Nlevel,1:hydrogen%nlevel) = RR(1:hydrogen%Nlevel,1:hydrogen%Nlevel) + CC(1:hydrogen%Nlevel,1:hydrogen%Nlevel)
       hydrogen%Gamma(:,:,id) = RR(1:hydrogen%Nlevel,1:hydrogen%Nlevel) + CC(1:hydrogen%Nlevel,1:hydrogen%Nlevel)

       !now compute the diagonal elements
       do i=1, hydrogen%Nlevel !diagonal hence full rate matrix
          ! 				gamma_tot(i,i) = sum(-gamma_tot(1:hydrogen%Nlevel,i)) !positive
          hydrogen%gamma(i,i,id) = -sum(hydrogen%gamma(:,i,id))
       enddo
       !stacking equations for helium
       !now helium
       if (helium_is_active) then
          xvar(ihel:Nlevel_atoms) = helium%n(:,icell)
          RR(:,:) = gamma_r(n_active,:,:) !init to radiative rates
          !we multiply ne by zero if lforce_lte.
          do kr=1,helium%Ncont
             !Rji in rate matrix is at rate(i,j)
             RR(helium%continua(kr)%i,helium%continua(kr)%j) = RR(helium%continua(kr)%i,helium%continua(kr)%j) * ne(icell)
          enddo
          !cumulate the derivative for helium equations
          call collision_matrix_atom(id, icell, helium, deriv(1:helium%Nlevel,1:helium%Nlevel))
          call collrate_matrix_atom(id, icell, helium, CC(1:helium%Nlevel,1:helium%Nlevel), deriv(1:helium%Nlevel,1:helium%Nlevel), &
               dgamcdne(ihel:Nlevel_atoms,ihel:Nlevel_atoms))

          ! 				gamma_tot(ihel:Nlevel_atoms,ihel:Nlevel_atoms) = RR(1:helium%Nlevel,1:helium%Nlevel) + CC(1:helium%Nlevel,1:helium%Nlevel)
          helium%gamma(:,:,id) = RR(1:helium%Nlevel,1:helium%Nlevel) + CC(1:helium%Nlevel,1:helium%Nlevel)
          do i=1, helium%Nlevel !diagonal hence full rate matrix
             ! 					gamma_tot(ihel+i - 1,ihel+i - 1) = sum(-gamma_tot(ihel:Nlevel_atoms,ihel+i - 1)) !positive
             helium%gamma(i,i,id) = sum(-helium%gamma(:,i,id)) !positive
          enddo
       endif
       !gamma_tot is 0 for SEE of H if N>Nlevel_H
       !so that the product with all pops (including nHe pops) is 0.

       !rate equation for this atom stored in f and df !
       !an equation per level dn_i/dt = 0 = sum_lp n_l * Gamma_lp_l
       call rate_equations_H_and_He(id, icell, Neq, ihel, dgamrdne(:,:), dgamcdne(:,:), xvar, fvar, dfvar)

       !charge conservation!
       call non_lte_charge_conservation_H_and_He (icell, neq, xvar, fvar, dfvar)

       !replace one equation of SEE by particle number conservation
       !particule conservation!

       call particle_conservation_H_and_He (icell, Neq, xvar, fvar, dfvar)

       !newton raphson!
       call multivariate_newton_raphson (neq, dfvar, fvar, xvar)
       ! 			call multivariate_newton_raphson (neq-1, dfvar(1:neq-1,1:neq-1), fvar(1:neq-1), xvar(1:neq-1))

       !update atomic populations and ne
       neg_pops = .false.
       do i=1, hydrogen%Nlevel
          hydrogen%n(i,icell) = hydrogen%n(i,icell) * ( 1.0 + fvar(i)/(1.0 + d_damp * abs(fvar(i))) )
          if (hydrogen%n(i,icell) < 0.0) neg_pops = .true.
       enddo
       if (helium_is_active) then
          do i=1, helium%Nlevel
             helium%n(i,icell) = helium%n(i,icell) * ( 1.0 + fvar(ihel+i-1)/(1.0 + d_damp * abs(fvar(ihel+i-1))) )
             if (helium%n(i,icell) < 0.0) neg_pops = .true.
          enddo
       endif
       !keep ne constant for tests
       ne(icell) = ne(icell) * ( 1.0 + fvar(neq)/(1.0 + d_damp * abs(fvar(neq))) )


       if ( (ne(icell) < 1d-16 * hydrogen%n(hydrogen%Nlevel,icell)).or.(neg_pops) ) rest_damping = .true.

       if (rest_damping .and. d_damp < (damp_char + 1.0)) then !restart with more iterations and larger damping (more stable, slower convergence)
          ! 				do i=1, hydrogen%Nlevel
          ! 					xvar(i) = ndag(i)
          ! 				enddo
          ! 				if (helium_is_active) then
          ! 					do i=1, helium%Nlevel
          ! 						xvar(ihel+i) = ndag(ihel+i)
          ! 					enddo
          ! 				endif
          xvar(1:neq-1) = ndag(:)
          xvar(neq) = nedag
          d_damp = damp_char * damp_scaling
          rest_damping = .false.
          n_iter = 0
          max_iter = damp_scaling * maxIter !can be large !

       elseif (rest_damping) then
          neg_pops = .false.
          lconverged = .false.
       endif

       !should be fractional
       delta_f = max(delta_f, maxval(abs(fvar)))
       dfpop = max(dfpop, maxval(abs(fvar(1:Neq-1))))
       dfne = max(dfne, abs(fvar(neq)))
       dfpop_h = max(dfpop_h, maxval(abs(fvar(1:hydrogen%Nlevel)))) !same as dfpop if only H
       if (helium_is_active) dfpop_he = max(dfpop_he,  maxval(abs(fvar(ihel:Neq-1))))

       if (verbose) then
          write(*,'("niter #"(1I4))') n_iter
          write(*,'("non-LTE ionisation delta="(1ES17.8E3)" dfH="(1ES17.8E3)" dfne="(1ES17.8E3)" dfHe="(1ES17.8E3) )') &
               delta_f, dfpop_H, dfne, dfpop_he
       endif

       if (n_iter > max_iter) then
          if (lsecond_try) then
             lconverged = .true.
             write(*,*) "Not enough iterations to converge", max_iter
             write(*,*) "err = ", delta_f
             write(*,*) " cell", icell, " id", id
             write(*,*) ""
          else
             lsecond_try = .true.
             lconverged = .false.
             n_iter = 0
             max_iter = damp_scaling * maxIter
             d_damp = damp_scaling * damp_char
             xvar(1:neq-1) = ndag(:) !should include helium
             xvar(neq) = nedag
          endif
       endif

       if ((delta_f < precision).and..not.rest_damping) then
          lconverged = .true.
          exit iterative_loop
       endif

    enddo iterative_loop !on convergence

    dne = abs(1.0 - nedag/ne(icell))
    dM = maxval(abs(1.0 - ndag(1:hydrogen%Nlevel)/hydrogen%n(:,icell)))
    if (helium_is_active) then
       dM_he = maxval(abs(1.0 - ndag(ihel:Nlevel_atoms)/helium%n(:,icell)))
       !!write(*,*) maxval(ndag(ihel:Nlevel_atoms)), maxval(helium%n(:,icell))
       if (verbose) write(*,'("(DELTA) non-LTE ionisation dM_H="(1ES17.8E3)" dne="(1ES17.8E3)" dM_He="(1ES17.8E3) )') &
            dM, dne, dM_he
       dM = max(dM, dM_he)
    else
       if (verbose) write(*,'("(DELTA) non-LTE ionisation dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne
    endif

    n_new(1,1:hydrogen%Nlevel,icell) = hydrogen%n(:,icell)
    ne_new(icell) = ne(icell) !will be set to ne once the new populations on all cells (with old quantities) have bee computed
    hydrogen%n(:,icell) = ndag(1:hydrogen%Nlevel) !the first value before iterations
    !needed for global convergence

    ne(icell) = nedag !reset because we have to loop over all cells
    !and we don't want to change the old values by the new one
    !until all  cells are treated
    if (helium_is_active) then
       n_new(helium%activeindex, 1:helium%Nlevel,icell) = helium%n(:,icell)
       helium%n(:,icell) = ndag(ihel:Nlevel_atoms)
    endif

    deallocate(xvar,fvar,dfvar)
    deallocate(gamma_r, ndag, RR, CC, deriv, dgamrdne, dgamcdne)
    !deallocate(gamma_tot)

    if (verbose) stop
    return
  end subroutine see_ionisation_nonlte_H_and_He

  subroutine see_ionisation_nonlte_helium(id, icell, dM, dne)
    !positivity tests + zero value below prec_pops

    integer, intent(in) :: id, icell
    real(kind=dp), intent(out) :: dM, dne
    integer :: nact, Neq, nact_start, n_iter, i, kr
    !smoothing parameter default if damp_char. If native pops, damp the iterations by a facteur damp_scaling larger
    real(kind=dp), parameter :: damp_scaling = 10.0, damp_char = 5.0, precision = 1d-5
    real(kind=dp) :: d_damp, nedag, delta_f, dfpop, dfne
    real(kind=dp), allocatable :: gamma_r(:,:), ndag(:), RR(:,:), CC(:,:), dgamrdne(:,:), dgamcdne(:,:), deriv(:,:)
    real(kind=dp), allocatable :: fvar(:), dfvar(:,:), xvar(:) !Jacobian
    logical :: lconverged, neg_pops, rest_damping, lsecond_try, verbose
    integer, parameter :: maxIter = 1000
    integer :: max_iter

    verbose = .false. !debug mode
    lconverged = .false.
    rest_damping = .false.
    d_damp = damp_char
    lsecond_try = .false. !set to .true. to avoid starting again if the first iterative scheme failed


    Neq = helium%Nlevel + 1
    allocate(xvar(neq)); xvar(:) = 0.0_dp
    !equations = Nlevel SEE per atom  + charge conservation (CC)
    allocate(fvar(neq)); fvar(:) = 0.0_dp
    !Jacobian of the systems = derivative of i equations to xvar of sum_j partial(fvar(i))/partial(xvar(j))
    allocate(dfvar(neq,neq)); dfvar(:,:) = 0.0_dp
    !old ne values to rest for the other cell. The new value is stored in ne_new for the next iterations
    nedag = ne(icell)
    !old populations of atoms for this cell
    !should include all levels of all atoms
    allocate(ndag(helium%Nlevel))
    ndag(:) = helium%n(:,icell)

    !Radiative rate matrix and collisional rate matrix
    allocate(RR(helium%Nlevel,helium%Nlevel), CC(helium%Nlevel,helium%Nlevel), dgamcdne(helium%Nlevel,helium%Nlevel), &
         deriv(helium%Nlevel,helium%Nlevel))
    !normalised to ne radiative rate matrix
    allocate(gamma_r(helium%Nlevel,helium%Nlevel),dgamrdne(helium%Nlevel,helium%Nlevel))
    !Gamma_r(i_cont,j_cont) = R_jcont_icont / ne
    call radrate_matrix_atom(id, icell, helium, gamma_r, dgamrdne)

    !start Newton-Raphson iterations
    if (verbose) write(*,*) "T = ", T(icell), " nHtot = ", nHtot(icell)

    n_iter = 0
    max_iter = maxIter

    iterative_loop : do while (.not.lconverged)
       delta_f = 0.0_dp
       dfpop = 0.0_dp !only populations
       dfne = 0.0_dp

       n_iter = n_iter + 1

       !set unknowns to old valu
       !only H at the moment
       xvar(1:helium%Nlevel) = helium%n(:,icell)

       !always the last one
       xvar(Neq) = ne(icell)
       fvar(:) = 0.0
       dfvar(:,:) = 0.0

       !multiply downward continuum rate by current ne
       RR = gamma_r !init to radiative rates
       !we multiply ne by zero if lforce_lte.
       do kr=1,helium%Ncont
          RR(helium%continua(kr)%i,helium%continua(kr)%j) = RR(helium%continua(kr)%i,helium%continua(kr)%j) * ne(icell)
       enddo
       !the derivative is independent of the value of ne
       ! 			if (.not.lforce_lte) call partial_radrate_matrix_partial_ne_atom(id,icell,hydrogen, dgamrdne)

       !get collisional matrix
       !
       !
       !but n(i)/nHII is held constant while it changes with ne.
       !
       !
       !Derivative is okay (indep of ne)
       !if (n_iter==1) update deriv otherwise no
       call collision_matrix_atom(id, icell, helium, deriv)
       ! 			!fill it into CC
       !derivative independent of the new value of ne but Collision need to be updated.
       !! deriv and dgmacdne cannot be the same at the moment ! dgamcdne is set to zero in the following.
       call collrate_matrix_atom(id, icell, helium, CC, deriv, dgamcdne)

       !overwrite the value accumulated. But for the first iteration it is the same (because of ne).
       helium%Gamma(:,:,id) = RR(:,:) + CC(:,:)
       !now compute the diagonal elements
       do i=1, helium%Nlevel !diagonal hence full rate matrix
          helium%Gamma(i,i,id) = sum(-helium%Gamma(:,i,id)) !positive
       enddo
       !equivalent
       ! 			call eliminate_delta(id, hydrogen)


       !rate equation for this atom stored in f and df !
       !an equation per level dn_i/dt = 0 = sum_lp n_l * Gamma_lp_l
       call rate_equations_atom(id, icell, Neq, helium, RR(:,:), dgamrdne(:,:), CC(:,:), dgamcdne(:,:), xvar, fvar, dfvar)

       !charge conservation!
       call non_lte_charge_conservation_He (icell, neq, xvar, fvar, dfvar)

       !replace one equation of SEE by particle number conservation
       !particule conservation!
       call particle_conservation (icell, Neq, helium, xvar, fvar, dfvar)

       !newton raphson!
       call multivariate_newton_raphson (neq, dfvar, fvar, xvar)
       ! 			call multivariate_newton_raphson (neq-1, dfvar(1:neq-1,1:neq-1), fvar(1:neq-1), xvar(1:neq-1))


       !update atomic populations and ne
       neg_pops = .false.
       do i=1, helium%Nlevel
          helium%n(i,icell) = helium%n(i,icell) * ( 1.0 + fvar(i)/(1.0 + d_damp * abs(fvar(i))) )
          if (helium%n(i,icell) < 0.0) neg_pops = .true.
       enddo
       !keep ne constant for tests
       ne(icell) = ne(icell) * ( 1.0 + fvar(neq)/(1.0 + d_damp * abs(fvar(neq))) )


       if ( (ne(icell) < 1d-16 * helium%n(helium%Nlevel,icell)).or.(neg_pops) ) rest_damping = .true.

       if (rest_damping .and. d_damp < (damp_char + 1.0)) then !restart with more iterations and larger damping (more stable, slower convergence)
          xvar(1:helium%Nlevel) = ndag(:)
          xvar(neq) = nedag
          d_damp = damp_char * damp_scaling
          rest_damping = .false.
          n_iter = 0
          max_iter = damp_scaling * maxIter !can be large !

       elseif (rest_damping) then
          neg_pops = .false.
          lconverged = .false.
       endif

       !should be fractional
       delta_f = max(delta_f, maxval(abs(fvar)))
       dfpop = max(dfpop, maxval(abs(fvar(1:Neq-1))))
       dfne = max(dfne, abs(fvar(neq)))

       if (verbose) then
          write(*,'("niter #"(1I4))') n_iter
          write(*,'("non-LTE ionisation delta="(1ES17.8E3)" dfn="(1ES17.8E3)" dfne="(1ES17.8E3) )') delta_f, dfpop, dfne
       endif


       if (n_iter > max_iter) then
          if (lsecond_try) then
             lconverged = .true.
             write(*,*) "Not enough iterations to converge", max_iter
             write(*,*) "err = ", delta_f
             write(*,*) " cell", icell, " id", id
             write(*,*) ""
          else
             lsecond_try = .true.
             lconverged = .false.
             n_iter = 0
             max_iter = damp_scaling * maxIter
             d_damp = damp_scaling * damp_char
             xvar(1:helium%Nlevel) = ndag(:)
             xvar(neq) = nedag
          endif
       endif

       if ((delta_f < precision).and..not.rest_damping) then
          lconverged = .true.
          exit iterative_loop
       endif

    enddo iterative_loop !on convergence

    dne = abs(1.0 - nedag/ne(icell))
    dM = maxval(abs(1.0 - ndag(:)/helium%n(:,icell)))
    if (verbose) write(*,'("(DELTA) non-LTE ionisation dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne

    n_new(helium%activeindex,1:helium%Nlevel,icell) = helium%n(:,icell)
    ne_new(icell) = ne(icell) !will be set to ne once the new populations on all cells (with old quantities) have bee computed
    helium%n(:,icell) = ndag(:) !the first value before iterations
    !needed for global convergence
    ne(icell) = nedag !reset because we have to loop over all cells
    !and we don't want to change the old values by the new one
    !until all  cells are treated

    deallocate(xvar,fvar,dfvar)
    deallocate(gamma_r, ndag, RR, CC, deriv, dgamrdne, dgamcdne)

    !almost the same routine as for hydrogen only, only the ionisation routine changes

    return
  end subroutine see_ionisation_nonlte_helium

  subroutine non_lte_charge_conservation_He (icell, neq, x, f, df)
    use atom_type, only : find_continuum
    integer, intent(in) :: icell, neq
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
    integer :: n, j, n_start, j0, j1
    type (element), pointer :: elem
    real(kind=dp) :: akj, nHeII, nHeIII
    real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk

    j0 = find_continuum(helium, 1) !ground state of HeII
    j1 = find_continuum(helium, j0) !heIII
    nHeII = sum(helium%stage(j0:j1-1)*helium%n(j0:j1-1,icell))
    nHeIII = helium%stage(j1) * helium%n(j1,icell)
    df(neq,neq) = 1/x(neq)**2 * (nheII + nheIII) !the factor stage is included!!
    f(neq) = 1.0 - (1.0/x(neq)) * (nheII + nheIII) !factor 2 in nheIII included
    !derivative to nheIII
    df(neq, neq-1) = -2.0 / x(neq)
    !derivative to nHeII!
    !or derivative to nHeII_i individual levels
    do j=1, j1-j0
       df(neq, neq-1-j) = -1/x(neq) !or for each level of nHeII I have the same derivative?
    enddo


    !now contribution from LTE atoms.
    do n=1, 26 !for 26 elements avoinding helium here
       !need to avoid other non-LTE elements
       elem => Elements(n)%ptr_elem
       if (elem%ID=="He") cycle

       call ionisation_frac_lte(elem, icell, x(neq), fjk, dfjk)

       do j=2, elem%Nstage
          !pure LTE term j = 1 corresponds to the first ionisation stage always
          !unlike in solve_ne where this loop is also for non-LTE models
          akj = (j-1) * elem%abund * nhtot(icell) !(j-1) = 0 for neutrals and 1 for singly ionised
          f(neq) = f(neq) - akj * fjk(j) / x(neq)
          df(neq,neq) = df(neq,neq) + akj * fjk(j) / x(neq)**2 - akj / x(neq) * dfjk(j)
       end do

       elem => null()

    enddo


    return
  end subroutine non_lte_charge_conservation_He

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

    call solve_lin(df, f, neq, .true.)

    return
  end subroutine multivariate_newton_raphson

  !store atom%dCdne like atom%C
  !or atom%lines%Cij, %Cji and for cont ??
  subroutine see_ionisation_nonlte_hydrogen(id, icell, dM, dne)
    !positivity tests + zero value below prec_pops

    integer, intent(in) :: id, icell
    real(kind=dp), intent(out) :: dM, dne
    integer :: nact, Neq, nact_start, n_iter, i, kr
    !smoothing parameter default if damp_char. If native pops, damp the iterations by a facteur damp_scaling larger
    real(kind=dp), parameter :: damp_scaling = 10.0, damp_char = 5.0, precision = 1d-5
    real(kind=dp) :: d_damp, nedag, delta_f, dfpop, dfne
    real(kind=dp), allocatable :: gamma_r(:,:), ndag(:), RR(:,:), CC(:,:), dgamrdne(:,:), dgamcdne(:,:), deriv(:,:)
    real(kind=dp), allocatable :: fvar(:), dfvar(:,:), xvar(:) !Jacobian
    logical :: lconverged, neg_pops, rest_damping, lsecond_try, verbose
    integer, parameter :: maxIter = 1000
    integer :: max_iter

    verbose = .false. !debug mode
    lconverged = .false.
    rest_damping = .false.
    d_damp = damp_char
    lsecond_try = .false. !set to .true. to avoid starting again if the first iterative scheme failed

    !list of indexes for non-LTE atoms included ??

    !number of equations = Nlevel_atom1 + Nlevel_atom2 + ... + 1 for ne
    Neq = hydrogen%Nlevel + 1
    !currently allocate and deallocate for each cell
    !futur fixed size
    !unknowns = populations of atomic levels + electron density
    allocate(xvar(neq)); xvar(:) = 0.0_dp
    !equations = Nlevel SEE per atom  + charge conservation (CC)
    allocate(fvar(neq)); fvar(:) = 0.0_dp
    !Jacobian of the systems = derivative of i equations to xvar of sum_j partial(fvar(i))/partial(xvar(j))
    allocate(dfvar(neq,neq)); dfvar(:,:) = 0.0_dp
    !old ne values to rest for the other cell. The new value is stored in ne_new for the next iterations
    nedag = ne(icell)
    !old populations of atoms for this cell
    !should include all levels of all atoms
    allocate(ndag(hydrogen%Nlevel))
    ndag(:) = hydrogen%n(:,icell)

    !Radiative rate matrix and collisional rate matrix
    allocate(RR(hydrogen%Nlevel,hydrogen%Nlevel), CC(hydrogen%Nlevel,hydrogen%Nlevel), dgamcdne(hydrogen%Nlevel,hydrogen%Nlevel), &
         deriv(hydrogen%Nlevel,hydrogen%Nlevel))
    !normalised to ne radiative rate matrix
    allocate(gamma_r(hydrogen%Nlevel,hydrogen%Nlevel),dgamrdne(hydrogen%Nlevel,hydrogen%Nlevel))
    !Gamma_r(i_cont,j_cont) = R_jcont_icont / ne
    call radrate_matrix_atom(id, icell, hydrogen, gamma_r, dgamrdne)
    !gamma_r(i,j) = -Rji, diagonal not included!
    !derivative of gamma_r wrt ne is included (with the diagonal!!). It is only non-zero
    !for recombinason rate and is independent of ne as it is Rji/ne (Rji = Rji/ne_old * ne_new at each iteration so derivative to ne_new is always Rji/ne_old)
    ! 		if (lforce_lte) then
    ! 			gamma_r = 0.0_dp
    ! 			dgamrdne = 0.0_dp
    ! 		endif

    !start Newton-Raphson iterations
    if (verbose) write(*,*) "T = ", T(icell), " nHtot = ", nHtot(icell)

    n_iter = 0
    max_iter = maxIter

    iterative_loop : do while (.not.lconverged)
       delta_f = 0.0_dp
       dfpop = 0.0_dp !only populations
       dfne = 0.0_dp

       n_iter = n_iter + 1

       !set unknowns to old valu
       !only H at the moment
       xvar(1:hydrogen%Nlevel) = hydrogen%n(:,icell)

       !always the last one
       xvar(Neq) = ne(icell)
       fvar(:) = 0.0
       dfvar(:,:) = 0.0

       !multiply downward continuum rate by current ne
       RR = gamma_r !init to radiative rates
       !we multiply ne by zero if lforce_lte.
       do kr=1,hydrogen%Ncont
          RR(hydrogen%continua(kr)%i,hydrogen%continua(kr)%j) = RR(hydrogen%continua(kr)%i,hydrogen%continua(kr)%j) * ne(icell)
       enddo
       !the derivative is independent of the value of ne
       ! 			if (.not.lforce_lte) call partial_radrate_matrix_partial_ne_atom(id,icell,hydrogen, dgamrdne)

       !get collisional matrix
       !
       !
       !but n(i)/nHII is held constant while it changes with ne.
       !
       !
       !Derivative is okay (indep of ne)
       !if (n_iter==1) update deriv otherwise no
       call collision_matrix_atom(id, icell, hydrogen, deriv)
       ! 			!fill it into CC
       !derivative independent of the new value of ne but Collision need to be updated.
       !! deriv and dgmacdne cannot be the same at the moment ! dgamcdne is set to zero in the following.
       call collrate_matrix_atom(id, icell, hydrogen, CC, deriv, dgamcdne)

       !overwrite the value accumulated. But for the first iteration it is the same (because of ne).
       hydrogen%Gamma(:,:,id) = RR(:,:) + CC(:,:)
       !now compute the diagonal elements
       do i=1, hydrogen%Nlevel !diagonal hence full rate matrix
          hydrogen%Gamma(i,i,id) = sum(-hydrogen%Gamma(:,i,id)) !positive
       enddo
       !equivalent
       ! 			call eliminate_delta(id, hydrogen)


       !rate equation for this atom stored in f and df !
       !an equation per level dn_i/dt = 0 = sum_lp n_l * Gamma_lp_l
       call rate_equations_atom(id, icell, Neq, hydrogen, RR(:,:), dgamrdne(:,:), CC(:,:), dgamcdne(:,:), xvar, fvar, dfvar)

       !charge conservation!
       call non_lte_charge_conservation_H (icell, neq, xvar, fvar, dfvar)

       !replace one equation of SEE by particle number conservation
       !particule conservation!
       call particle_conservation (icell, Neq, hydrogen, xvar, fvar, dfvar)

       !newton raphson!
       call multivariate_newton_raphson (neq, dfvar, fvar, xvar)
       ! 			call multivariate_newton_raphson (neq-1, dfvar(1:neq-1,1:neq-1), fvar(1:neq-1), xvar(1:neq-1))


       !update atomic populations and ne
       neg_pops = .false.
       do i=1, hydrogen%Nlevel
          hydrogen%n(i,icell) = hydrogen%n(i,icell) * ( 1.0 + fvar(i)/(1.0 + d_damp * abs(fvar(i))) )
          if (hydrogen%n(i,icell) < 0.0) neg_pops = .true.
       enddo
       !keep ne constant for tests
       ne(icell) = ne(icell) * ( 1.0 + fvar(neq)/(1.0 + d_damp * abs(fvar(neq))) )


       if ( (ne(icell) < 1d-16 * hydrogen%n(hydrogen%Nlevel,icell)).or.(neg_pops) ) rest_damping = .true.

       if (rest_damping .and. d_damp < (damp_char + 1.0)) then !restart with more iterations and larger damping (more stable, slower convergence)
          do i=1, hydrogen%Nlevel
             xvar(i) = ndag(i)
          enddo
          xvar(neq) = nedag
          d_damp = damp_char * damp_scaling
          rest_damping = .false.
          n_iter = 0
          max_iter = damp_scaling * maxIter !can be large !

       elseif (rest_damping) then
          neg_pops = .false.
          lconverged = .false.
       endif

       !should be fractional
       delta_f = max(delta_f, maxval(abs(fvar)))
       dfpop = max(dfpop, maxval(abs(fvar(1:Neq-1))))
       dfne = max(dfne, abs(fvar(neq)))

       if (verbose) then
          write(*,'("niter #"(1I4))') n_iter
          write(*,'("non-LTE ionisation delta="(1ES17.8E3)" dfn="(1ES17.8E3)" dfne="(1ES17.8E3) )') delta_f, dfpop, dfne
       endif


       if (n_iter > max_iter) then
          if (lsecond_try) then
             lconverged = .true.
             write(*,*) "Not enough iterations to converge", max_iter
             write(*,*) "err = ", delta_f
             write(*,*) " cell", icell, " id", id
             write(*,*) ""
          else
             lsecond_try = .true.
             lconverged = .false.
             n_iter = 0
             max_iter = damp_scaling * maxIter
             d_damp = damp_scaling * damp_char
             do i=1, hydrogen%Nlevel
                xvar(i) = ndag(i)
             enddo
             xvar(neq) = nedag
          endif
       endif

       if ((delta_f < precision).and..not.rest_damping) then
          lconverged = .true.
          exit iterative_loop
       endif

    enddo iterative_loop !on convergence

    dne = abs(1.0 - nedag/ne(icell))
    dM = maxval(abs(1.0 - ndag(:)/hydrogen%n(:,icell)))
    if (verbose) write(*,'("(DELTA) non-LTE ionisation dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne

    n_new(1,1:hydrogen%Nlevel,icell) = hydrogen%n(:,icell)
    ne_new(icell) = ne(icell) !will be set to ne once the new populations on all cells (with old quantities) have bee computed
    hydrogen%n(:,icell) = ndag(:) !the first value before iterations
    !needed for global convergence
    ne(icell) = nedag !reset because we have to loop over all cells
    !and we don't want to change the old values by the new one
    !until all  cells are treated

    deallocate(xvar,fvar,dfvar)
    deallocate(gamma_r, ndag, RR, CC, deriv, dgamrdne, dgamcdne)


    return
  end subroutine see_ionisation_nonlte_hydrogen

  subroutine non_lte_charge_conservation_H (icell, neq, x, f, df)
    !also computes the derivative of the charge conservation equation wuth
    !ne and total population of stages with stage > 0 (like nH(Nlevel)=nHII)
    integer, intent(in) :: icell, neq
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(inout) :: df(neq,neq), f(neq)
    integer :: n, j, n_start
    type (element), pointer :: elem
    real(kind=dp) :: akj, np_j
    real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk
    !Separate in two parts: 1 non-LTE atoms and 2 electrons from LTE elements

    !derivative of CC (charge conservation) to ne
    !x = (n1,n2,n3 .... ne)
    np_j = 1.0_dp * hydrogen%n(hydrogen%Nlevel,icell) !the 1.0_dp represents the ionisation stage of nHII (2 for nHeIII)
    df(neq,neq) = 1/x(neq)**2 * np_j!(sum(hydrogen%n(:,icell)*hydrogen%stage(:)))
    !in principle if stage=0 we only add 0.
    !for H it is only hydrogen%(nlevel,icell)
    !more atoms can be added
    !check:
    !if (hydrogen%n(hydrogen%Nlevel,icell) /= sum(hydrogen%n(:,icell)*hydrogen%stage(:)))
    !derivative to ions currently only nHII so the equation before neq is nHII
    df(neq,neq-1) = -1.0 / x(neq) !pos_HII = neq - 1 if only hydrogen
    !for helium II it is -1/x(neq) df(neq, pos_heII)
    !for helium III it is -2/x(neq) df(neq, pos_heIII)

    !non-LTE part of the charge conservation
    f(neq) = 1.0 - (1.0 / x(neq)) * np_j!(sum(hydrogen%n(:,icell)*hydrogen%stage(:)))
    !for H only should be 1 - hydrogen%n(hydrogen%Nlevel,icell)/x(eq) that is 1 - x(eq-1)/x(eq)

    !now contribution from LTE atoms.
    do n=2, 26 !for 26 elements avoinding hydrogen here
       !need to avoid other non-LTE elements
       elem => Elements(n)%ptr_elem

       call ionisation_frac_lte(elem, icell, x(neq), fjk, dfjk)

       do j=2, elem%Nstage
          !pure LTE term j = 1 corresponds to the first ionisation stage always
          !unlike in solve_ne where this loop is also for non-LTE models
          akj = (j-1) * elem%abund * nhtot(icell) !(j-1) = 0 for neutrals and 1 for singly ionised
          f(neq) = f(neq) - akj * fjk(j) / x(neq)
          df(neq,neq) = df(neq,neq) + akj * fjk(j) / x(neq)**2 - akj / x(neq) * dfjk(j)
       end do

       elem => null()

    enddo

    ! 	write(*,*) "df(neq,neq)=", df(neq, neq)
    ! 	write(*,*) "df(:,neq)=", df(:,neq)
    ! 	write(*,*) "df(neq,:)=",df(neq,:)
    ! 	write(*,*) "total=",df
    ! 	stop

    return
  end subroutine non_lte_charge_conservation_H

  subroutine rate_equations_atom(id,icell, neq, atom, gamma_rad, dgrdne, gamma_col,  dgcdne, x, f, df)
    integer, intent(in) :: id, icell, neq
    type (atomtype), intent(in) :: atom !contains total rate matrix without diagonal in atom%Gamma
    real(kind=dp), intent(in) :: x(neq)
    real(kind=dp), intent(in), dimension(atom%Nlevel,atom%Nlevel) :: gamma_rad, dgrdne, gamma_col, dgcdne
    real(kind=dp), intent(out) ::  f(neq), df(neq,neq)
    integer :: i, j, l, lp
    real(kind=dp) :: part1, part2
    !individual equation for each level

    !also computes the derivative of each equation by each level

    !-> here Gamma diag is filled (with eliminate_delta or explicit summation gamma(:,l,id) for each level l)
    f(1:atom%Nlevel) = matmul( atom%gamma(:,:,id), x(1:atom%Nlevel) )
    !derivative of l row for each lp column
    !derivative of each SEE (a l row of Gamma) for all n(i)
    !it is simply Gamma(i,l) since Gamma does not depend on n(i)
    do l=1, atom%Nlevel
       do lp=1,atom%nlevel
          df(l,lp) = atom%gamma(l,lp,id)!hydrogen%gamma(l,lp,id)
          !for row l computes the lp derivatives of SEE(l)=f(k) with respect to n(lp)
       enddo
    enddo

    !derivative of SEE to ne
    !dCij/dne = Cij/ne
    !dCji/dne = 2*Cji/ne
    !dRji_cont/dne = Rji_cont/ne

    df(1:neq-1,neq) = matmul( dgcdne, x(1:atom%Nlevel) )!hydrogen%
    if (.not.lforce_lte) then
       df(1:neq-1,neq) = df(1:neq-1,neq) + matmul( dgrdne, x(1:atom%Nlevel) )
    endif

    return

    write(*,*) "should not be here"
    !-> really needed to avoid some temporary storage.
    !To Do
    !-> here are more optimised code that we need to avoid many loops over levels

    !-> here diag(Gamma) must be zero
    do i=1, atom%Nlevel
       part1 = -sum(atom%Gamma(:,i,id))
       f(i) = x(i) * part1!SEE for the first level
       do j=1, atom%Nlevel
          f(i) = f(i) + x(j) * atom%Gamma(i,j,id)
          if (i==j) then
             df(i,j) = part1
          else
             df(i,j) = atom%Gamma(i,j,id) !index ?
          endif
       enddo

    enddo

    do i=1, atom%Nlevel

       !add derivative to ne for recombinaison rates and collision rates

       if (i==atom%Nlevel) then
          part1 = -sum(gamma_rad(:,i)) / x(neq)
          part2 = -sum(gamma_col(:,i)) / x(neq)
          df(i,neq) = x(i) * (part1 + 2.0 * part2)
          do j=1, atom%Nlevel

             df(i,neq) = df(i,neq) + x(j)/x(neq) * gamma_col(i,j)

          enddo
       else

          part1 = -sum(gamma_col(:,i)) / x(neq)
          df(i,neq) = x(i) * part1

          do j=1, atom%Nlevel

             df(i,neq) = df(i,neq) + x(j) * gamma_col(i,j) / x(neq)

          enddo

          !recombination
          !x(atom%Nlevel) is np in that case
          df(i,neq) = df(i,neq) + (x(i) / x(neq)) * gamma_col(i,atom%Nlevel) + &
               (x(atom%Nlevel)/x(neq))*(gamma_rad(atom%Nlevel,i)+2.0*gamma_col(atom%Nlevel,i))
       endif
    enddo


    return
  end subroutine rate_equations_atom

  subroutine particle_conservation (icell, Neq, atom, x, f, df)
    integer, intent(in) :: icell, Neq
    real(kind=dp), intent(in) :: x(neq)
    type(AtomType), intent(in) :: atom
    real(kind=dp), intent(inout) :: f(neq), df(neq, neq)
    integer :: imaxpop, nlev
    real(kind=dp) :: ntotal

    ntotal = ntotal_atom(icell,atom)!hydrogen

    nlev = Neq - 1 !remove electrons

    imaxpop = locate(x(1:nlev), maxval(x(1:nlev)))

    !mass conservation: ntotal = sum_i n(i)
    !eq imaxpop is replaced by:
    !f(imaxpop) = ntotal - sum_i n(i) = 0
    !or
    !f(imaxpop) = 1.0 - sum_i n(i) / ntotal = 0

    f(imaxpop) = 1.0_dp
    f(imaxpop) = f(imaxpop) -sum(x(1:nlev)) / ntotal

    df(imaxpop,neq) = 0.0_dp
    df(imaxpop,1:nlev) = -1.0_dp / ntotal

    return
  end subroutine particle_conservation

  subroutine update_populations_and_electrons(id, icell, delta, verbose, nit, iterate_ne)
    ! --------------------------------------------------------------------!
    ! Performs a solution of SEE for each atom and simultneously solves
    ! the charge conservation equations for hydrogen
    ! H only at the moment
    !
    ! To do: implements Ng's acceleration iteration here
    ! --------------------------------------------------------------------!


    integer, intent(in) :: id, icell, nit
    logical, intent(in) :: verbose
    logical, intent(in) :: iterate_ne
    type(AtomType), pointer :: atom
    integer :: nact, nact_start
    logical :: nonlte_ionisation
    real(kind=dp) :: dM, dT, dTex, dpop, Tion, Tex
    real(kind=dp), intent(out) :: delta
    real(kind=dp) :: dne

    dpop = 0.0_dp
    dtex = 0.0_dp
    ! if (icell /= n_cells) return
    !electron density not iterated ?
    !at the moment direcly call the old routine!
    !works also for rest iterations (iterations where ne not iterated) ?
    if (.not.iterate_ne) then
       call update_populations(id, icell, dM, verbose, nit)
       return
    endif

    !here iterate_ne is .true. otherwise we exit after entering in update_populations
    ! 		if (helium_is_active .or. hydrogen%active) nonlte_ionisation = .true.
    !global variable helium_is_active is true if helium is associated and active.
    !it avoids to double test the associate and the active statuses.

    !-> currenlty H is always active if H active

    if (verbose) write(*,'(" niter #"(1I4)" id #"(1I4))') nit, id


    nact = 1 !future test on all elements that can be included in the equation
    ! 		if (hydrogen%active .and. activeatoms(nact)%ptr_atom%id=='H') then
    if (hydrogen%active) then
       if (verbose) write(*,*) " Solving SEE and charge conservation for H+He+ne"
       !actually includes all non-LTE atom part of the equation
       !loop for all non-LTE atoms included in the equations
       !currently only H and He

       if (helium_is_active) then
          call see_ionisation_nonlte_H_and_He(id, icell, dM, dne)
          if (verbose) then
             write(*,'(" --> non-LTE ionisation (H+He+ne) dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne
          endif
          dpop = max(dpop, dM)
          call calc_delta_Tex_atom(icell, hydrogen, dT, Tex, Tion,.false.)
          if (verbose) then
             write(*,*) " Hydrogen:"
             write(*,'(" --> dT="(1ES17.8E3)" Te="(1ES17.8E3)" K")'), dT, T(icell)
             write(*,'(" >> Tex="(1ES17.8E3)" K Tion="(1ES17.8E3)" K")') Tex, Tion
          endif
          dTex = max(dTex, dT)
          call calc_delta_Tex_atom(icell, helium, dT, Tex, Tion,.false.)
          dTex = max(dTex, dT)
          if (verbose) then
             write(*,*) " Helium:"
             write(*,'(" --> dT="(1ES17.8E3)" Te="(1ES17.8E3)" K")'), dT, T(icell)
             write(*,'(" >> Tex="(1ES17.8E3)" K Tion="(1ES17.8E3)" K")') Tex, Tion
          endif
       else
          call see_ionisation_nonlte_hydrogen(id, icell, dM, dne)
          if (verbose) then
             write(*,'(" --> non-LTE ionisation (H+ne) dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne
          endif
          dpop = max(dpop, dM)
          call calc_delta_Tex_atom(icell, hydrogen, dT, Tex, Tion,.false.)
          !For one atoms = for all transitions
          if (verbose) then
             write(*,'(" --> dT="(1ES17.8E3)" Te="(1ES17.8E3)" K")'), dT, T(icell)
             write(*,'(" >> Tex="(1ES17.8E3)" K Tion="(1ES17.8E3)" K")') Tex, Tion
          endif
          dTex = max(dTex, dT)
       endif

    else
       if (helium_is_active) then
          call see_ionisation_nonlte_Helium(id, icell, dM, dne)
          if (verbose) then
             write(*,'(" --> non-LTE ionisation (He+ne) dM="(1ES17.8E3)" dne="(1ES17.8E3) )') dM, dne
          endif
          dpop = max(dpop, dM)
          call calc_delta_Tex_atom(icell, helium, dT, Tex, Tion,.false.)
          if (verbose) then
             write(*,'(" --> dT="(1ES17.8E3)" Te="(1ES17.8E3)" K")'), dT, T(icell)
             write(*,'(" >> Tex="(1ES17.8E3)" K Tion="(1ES17.8E3)" K")') Tex, Tion
          endif
          dTex = max(dTex, dT)
       endif

       !at the moment separate hydrogen and the others
       !hydrogen always the first element. So if hydrogen is active, it is the first one!
       !-> beware, if hydrogen is not active , nact should be 1
       !-> check for helium too.
       if (hydrogen%active) then
          nact_start = 2
       else
          nact_start = 1
       endif
       aatom_loop : do nact=nact_start, Nactiveatoms
          !skip helium if active as it is included in the nonlte ionisation
          if (activeatoms(nact)%ptr_atom%ID=="He") then
             cycle aatom_loop
          endif

          if (verbose) write(*,*) " Solving SEE for atom ", atom%ID
          atom => activeatoms(nact)%ptr_atom
          call SEE_atom(id, icell, atom, dM)
          call calc_delta_Tex_atom(icell, atom, dT, Tex, Tion,.false.)

          !For one atoms = for all transitions
          if (verbose) then
             write(*,'(" --> dM="(1ES17.8E3))') dM
             write(*,'(" --> dT="(1ES17.8E3)" Te="(1ES17.8E3)" K")'), dT, T(icell)
             write(*,'(" >> Tex="(1ES17.8E3)" K Tion="(1ES17.8E3)" K")') Tex, Tion
          endif

          !max for all atoms
          dpop = max(dpop, dM)
          ! and all transitions...
          dTex = max(dTex, dT)

          atom => NULL()

       enddo aatom_loop
    endif

    !flag the one with a "*"
    !delta = dTex
    delta = dM

    !Compare all atoms
    if (verbose) then
       write(*,'(" <:> cell #"(1I4)" dTmax="(1ES17.8E3)" dnmax="(1ES17.8E3))') icell, dT, dpop
    endif

    return
  end subroutine update_populations_and_electrons

  subroutine see_atom(id, icell,atom, dM)
    ! --------------------------------------------------------------------!
    ! For atom atom solves for the Statistical Equilibrium Equations (SEE)
    !
    ! write matrix here
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
    real(kind=dp) :: n0 = 0.0_dp, ntotal, Gamma_dag(atom%Nlevel,atom%Nlevel) !debug
    real(kind=dp) :: dn_n

    if (n0 /= 0.0_dp) call error("n0 should be zero at the moment !")
    ntotal = ( ntotal_atom(icell,atom) - n0 )

    ndag(:) = atom%n(:,icell) / ntotal
    imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
    !imaxpop = atom%Nlevel

    call eliminate_delta(id, atom)
    !check that diagonal is > 0 and off-diagonals < 0
    ! 		call check_elements(id, icell, atom%Nlevel, atom%Gamma(:,:,id), "in see_atom: "//atom%ID)

    atom%n(:,icell) = 0d0
    atom%n(imaxpop,icell) = 1.0_dp
    ! 		!atom%Gamma(imaxpop,:,id) = 1.0_dp


    Gamma_dag = atom%Gamma(:,:,id) !with diagonal
    atom%Gamma(imaxpop,:,id) = 1.0_dp

	if (ldamp_jacobi) then
    !-> solve for residual. By default omega_sor_atom is 1.0
    	delta = atom%n(:,icell) - matmul(atom%Gamma(:,:,id), ndag)
    	call GaussSlv(atom%Gamma(:,:,id), delta(:), atom%Nlevel)
   		atom%n(:,icell) = ndag(:) + omega_sor_atom(atom%activeindex) * delta(:)
   	else !Jacobi
    	call GaussSlv(atom%Gamma(:,:,id), atom%n(:,icell), atom%Nlevel)
   	endif
!     call GaussSlv(atom%Gamma(:,:,id), atom%n(:,icell), atom%Nlevel)
!     !call solve_lin(atom%Gamma(:,:,id), atom%n(:,icell), atom%Nlevel, .true.)

    if ((maxval(atom%n(:,icell)) < 0.0)) then
       !raise warning or error if all populations are negative. Otherwise, handle
       !the negative populations to not stop the calculations ??
       ! 		if ((maxval(atom%n(:,icell)) <= 0.0).or.(minval(atom%n(:,icell))<0.0)) then
       write(*,*) atom%ID, id, icell
       write(*,'("nstar: "*(ES14.5E3))') (ndag(l),l=1,atom%Nlevel) !*ntotal
       write(*,'("n: "*(ES14.5E3))') (atom%n(l,icell),l=1,atom%Nlevel) !*ntotal
       write(*,'("b: "*(ES14.5E3))') (atom%n(l,icell)/ndag(l),l=1,atom%Nlevel) !*ntotal
       do l=1, atom%Nlevel
          write(*, '(1I3, *(ES14.5E3))') l, (atom%Gamma(l,lp,id), lp=1, atom%Nlevel)
       enddo
       ! 			call warning("Negative pops after SEE!")
       call warning("All pops are negative after SEE!")
       !or error ?
    endif

    if ((any_nan_infinity_vector(atom%n(:,icell))>0)) then
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
       write(*,*) "Radiative rates"
       do l=1, atom%Nline
          write(*,*) " line ", atom%lines(l)%i, atom%lines(l)%j
          write(*,*) "-> Rij=",atom%lines(l)%Rij(id)," Rji=",atom%lines(l)%Rji(id)
       enddo
       do l=1, atom%Ncont
          write(*,*) " cont ", atom%continua(l)%i, atom%continua(l)%j
          write(*,*) "-> Rij=",atom%continua(l)%Rij(id)," Rji=",atom%continua(l)%Rji(id)
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

       !Small populations are kept but not used in the convergence test
       !and in Tex. Populations below prec_pops (including negative) are set
       !to zero.
       ! 			if (atom%n(l,icell) < frac_ne_limit * ne(icell)) then
       if (atom%n(l,icell) < frac_limit_pops * ntotal) then
          Nsmall_pops = Nsmall_pops + 1
          level_index(lp) = l
          lp = lp + 1
          if (atom%n(l,icell) <= prec_pops * ntotal) then
             atom%n(l,icell) = 0.0_dp
             Nneg_or_null_pops = Nneg_or_null_pops + 1
          endif
       else
          dn_n = (1.0_dp - ndag(l))/atom%n(l,icell)
          dM = max(dM, abs(1.0_dp - ndag(l)/atom%n(l,icell)))
          ! 				dM = max(dM, abs(atom%n(l,icell)-ndag(l))/ndag(l))
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


  SUBROUTINE update_populations(id, icell, delta, verbose, nit)
    ! --------------------------------------------------------------------!
    ! Performs a solution of SEE for each atom.
    !
    ! to do: implements Ng's acceleration iteration here
    ! --------------------------------------------------------------------!

    integer, intent(in) :: id, icell, nit
    logical, intent(in) :: verbose
    type(AtomType), pointer :: atom
    integer :: nat, nati, natf
    logical :: accelerate = .false.
    real(kind=dp) :: dM, dT, dTex, dpop, Tion, Tex
    real(kind=dp), intent(out) :: delta


    dpop = 0.0_dp
    dTex = 0.0_dp

    if (verbose) write(*,'(" niter #"(1I4)" id #"(1I4))') nit, id

    nati = 1; natf = Nactiveatoms
    do nat=nati,natf !loop over each active atoms

       atom => ActiveAtoms(nat)%ptr_atom

       call SEE_atom(id, icell, atom, dM)

       !compute dTex and new values of Tex
       call calc_delta_Tex_atom(icell, atom, dT, Tex, Tion,.true.)


       !For one atoms = for all transitions
       if (verbose) then
          write(*,'(" --> dM="(1ES17.8E3))') dM
          write(*,'(" --> dT="(1ES17.8E3)" Te="(1ES17.8E3)" K")'), dT, T(icell)
          write(*,'(" >> Tex="(1ES17.8E3)" K Tion="(1ES17.8E3)" K")') Tex, Tion
       endif

       !max for all atoms
       dpop = max(dpop, dM)
       ! and all transitions...
       dTex = max(dTex, dT)

       atom => NULL()

    enddo

    !flag the one with a "*"
    !delta = dTex
    delta = dM

    !Compare all atoms
    if (verbose) then
       write(*,'(" <:> cell #"(1I4)" dTmax="(1ES17.8E3)" dnmax="(1ES17.8E3))') icell, dT, dpop
    endif

    RETURN
  END SUBROUTINE update_populations


  subroutine calc_rates_mali(id, icell, iray, waq)
    integer, intent(in) :: id, icell, iray
    real(kind=dp), intent(in) :: waq
    real(kind=dp) :: chicont, etacont
    !takes the max shift of lines
    real(kind=dp), dimension(Nlambda_max_trans) :: Ieff
    integer :: nact, Nred, Nblue, kc, kcc, kr, krr, i, j, ip, jp, l, a0, Nl, icell_d
    integer :: N1, N2
    real(kind=dp) :: etau, wphi, JJ, JJb, j1, j2, I1, I2, ehnukt, twohnu3_c2
    real(kind=dp) :: wl, a1, a2
    type(AtomType), pointer :: aatom, atom
    real(kind=dp) :: wi, wj, nn, dk0, xcc_ij, xcc_ji, ni_on_nj_star

    aatom_loop : do nact = 1, Nactiveatoms
       aatom => ActiveAtoms(nact)%ptr_atom


       atr_loop : do kr = 1, aatom%Ntr_line

          kc = aatom%at(kr)%ik
          !all lines contribute, so don't evaluate lcontrib_to_opac

          Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
          i = aatom%lines(kc)%i
          j = aatom%lines(kc)%j

          !dk should be zero now
          Nl = Nred-dk_min+dk_max-Nblue+1
          a0 = Nblue+dk_min-1

          JJ = 0.0
          wphi = 0.0
          xcc_ji = 0.0

          Ieff(1:Nl) = Itot(Nblue+dk_min:Nred+dk_max,iray,id) - Psi(Nblue+dk_min:Nred+dk_max,iray,id) * &
               eta_atoms(Nblue+dk_min:Nred+dk_max,nact,id)
          if (lsobolev_regime) then
            Ieff(1:Nl) = 0.0
          endif

          ! 				enddo
          do l=1,Nl
             if (l==1) then
                wl = 0.5*(1d3*hv)
                !wl = 0.5*(lambda(a0+l+1)-lambda(a0+l)) * clight / aatom%lines(kc)%lambda0
             elseif (l==Nl) then
                wl = 0.5*(1d3*hv)
                !wl = 0.5*(lambda(a0+l)-lambda(a0+l-1)) * clight / aatom%lines(kc)%lambda0
             else
                wl = 1d3*hv
                !wl = 0.5*(lambda(a0+l+1)-lambda(a0+l-1)) * clight / aatom%lines(kc)%lambda0
             endif

             JJ = JJ + wl * Ieff(l) * aatom%lines(kc)%phi_loc(l,iray,id)
             wphi = wphi + wl * aatom%lines(kc)%phi_loc(l,iray,id)
             xcc_ji = xcc_ji + chi_up(a0+l,i,nact,id)*psi(a0+l,iray,id)*Uji_down(a0+l,j,nact,id)

          enddo

		  !Trick here, it is however dangerous!
		  !wphi should never be 0 or < 0 ! or there is a problem with the local profile.
          if (wphi > 0.0) then	  
          	JJ = waq * JJ / wphi
          else
          	JJ = 0.0
          endif

          !init at Aji
          aatom%lines(kc)%Rji(id) = aatom%lines(kc)%Rji(id) + JJ * aatom%lines(kc)%Bji
          aatom%lines(kc)%Rij(id) = aatom%lines(kc)%Rij(id) + JJ * aatom%lines(kc)%Bij

          !cross-coupling terms
          xcc_ij = 0.0_dp!
          sub_atr_loop : do krr=1, aatom%Ntr

             kcc = aatom%at(krr)%ik

             select case (aatom%at(krr)%trtype)

             case ("ATOMIC_LINE")
                ip = aatom%lines(kcc)%i
                jp = aatom%lines(kcc)%j
                N1 = aatom%lines(kcc)%Nblue+dk_min
                N2 = aatom%lines(kcc)%Nred+dk_max
                if (jp==i) then !i upper level of another transitions
                   xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
                endif

             case ("ATOMIC_CONTINUUM")
                ip = aatom%continua(kcc)%i
                jp = aatom%continua(kcc)%j
                N1 = aatom%continua(kcc)%Nb
                N2 = aatom%continua(kcc)%Nr
                if (jp==i) then !i upper level of another transitions
                   xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
                endif

             case default
                call error("Transition type unknown cross-coupling (line)", aatom%at(krr)%trtype)
             end select


          enddo sub_atr_loop

          aatom%lines(kc)%Rji(id) = aatom%lines(kc)%Rji(id) - xcc_ji * waq
          aatom%lines(kc)%Rij(id) = aatom%lines(kc)%Rij(id) + xcc_ij * waq

       end do atr_loop

       atrc_loop : do kr = aatom%Ntr_line+1, aatom%Ntr

          kc = aatom%at(kr)%ik
          !all transitions contribue so I don't test lcontrib_to_opac

          i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
          ! 				chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
          ! 				neff = aatom%stage(j) * sqrt(aatom%Rydberg / (aatom%E(j) - aatom%E(i)) )

          ! 				ni_on_nj_star = ne(icell) * phi_T(icell, aatom%g(i)/aatom%g(j), aatom%E(j)-aatom%E(i))
          ni_on_nj_star = aatom%nstar(i,icell)/aatom%nstar(j,icell)

          Nblue = aatom%continua(kc)%Nb; Nred = aatom%continua(kc)%Nr
          Nl = Nred-Nblue + 1

          JJ = 0.0
          JJb = 0.0
          xcc_ji = 0.0

          Ieff(1:Nl) = Itot(Nblue:Nred,iray,id) - Psi(Nblue:Nred,iray,id) * eta_atoms(Nblue:Nred,nact,id)


          icell_d = 1
          if (ldissolve) then
             if (aatom%ID=="H") icell_d = icell
          endif

          do l=1, Nl
             if (l==1) then
                wl = 0.5*(lambda(Nblue+l)-lambda(Nblue)) / lambda(Nblue)
             elseif (l==Nl) then
                wl = 0.5*(lambda(Nred)-lambda(Nred-1))  / lambda(Nred)
             else
                wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-2)) / lambda(Nblue+l-1)
             endif

             a1 = aatom%continua(kc)%alpha_nu(l,icell_d)

             JJ = JJ + wl  * a1 * Ieff(l)

             ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-1))
             twohnu3_c2 = twohc/lambda(Nblue+l-1)**3
             JJb = JJb + wl * ( Ieff(l) + twohnu3_c2 ) * ehnukt * a1

             xcc_ji = xcc_ji + chi_up(Nblue+l-1,i,nact,id)*Uji_down(Nblue+l-1,j,nact,id)*psi(Nblue+l-1,iray,id)

          enddo

          aatom%continua(kc)%Rij(id) = aatom%continua(kc)%Rij(id) + waq*fourpi_h * JJ
          aatom%continua(kc)%Rji(id) = aatom%continua(kc)%Rji(id) + waq*fourpi_h * JJb * ni_on_nj_star

          xcc_ij = 0.0_dp
          sub_atrc_loop : do krr=1, aatom%Ntr

             kcc = aatom%at(krr)%ik

             select case (aatom%at(krr)%trtype)
             case ("ATOMIC_LINE")
                ip = aatom%lines(kcc)%i
                jp = aatom%lines(kcc)%j
                N1 = aatom%lines(kcc)%Nblue+dk_min
                N2 = aatom%lines(kcc)%Nred+dk_max
                if (jp==i) then !i upper level of another transitions
                   xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
                endif
             case ("ATOMIC_CONTINUUM")
                ip = aatom%continua(kcc)%i
                jp = aatom%continua(kcc)%j
                N1 = aatom%continua(kcc)%Nb
                N2 = aatom%continua(kcc)%Nr
                if (jp==i) then !i upper level of another transitions
                   xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
                endif

             case default
                call error("Transition type unknown cross-coupling (continuum)", aatom%at(krr)%trtype)
             end select


          enddo sub_atrc_loop

          aatom%continua(kc)%Rji(id) = aatom%continua(kc)%Rji(id) - xcc_ji * waq
          aatom%continua(kc)%Rij(id) = aatom%continua(kc)%Rij(id) + xcc_ij * waq

       end do atrc_loop

       aatom => NULL()

    end do aatom_loop

    return
  end subroutine calc_rates_mali

  ! 	subroutine calc_rates_mali_old(id, icell, iray, waq)
  ! 		integer, intent(in) :: id, icell, iray
  ! 		real(kind=dp), intent(in) :: waq
  ! 		real(kind=dp) :: chicont, etacont
  ! 		real(kind=dp), dimension(Nlambda_max_trans) :: Ieff
  ! 		!real(kind=dp), dimension(Nlambda) :: Ieff!chi_tot, eta_tot,
  ! 		integer :: nact, Nred, Nblue, kc, kcc, kr, krr, i, j, ip, jp, l, a0, Nl, icell_d
  ! 		integer :: N1, N2
  ! 		real(kind=dp) :: etau, wphi, JJ, JJb, j1, j2, I1, I2, ehnukt, twohnu3_c2
  ! 		real(kind=dp) :: wl, a1, a2
  ! 		type(AtomType), pointer :: aatom, atom
  ! 		real(kind=dp) :: wi, wj, nn, dk0, xcc_ij, xcc_ji
  !
  ! 		aatom_loop : do nact = 1, Nactiveatoms
  ! 			aatom => ActiveAtoms(nact)%ptr_atom
  !
  ! 			!move inside transition Ieff(1:Nlam) = Itot(Nb:Nred) ...
  ! 			Ieff = Itot(:,iray,id) - Psi(:,iray,id) * eta_atoms(:,nact,id)
  ! 			if (any_nan_infinity_vector(Ieff)>0) then
  ! 				write(*,*) id, icell, iray, waq
  ! 				write(*,*) "Ieff=", Ieff
  ! 				call error("Infinity in Ieff")
  ! 			endif
  !
  ! 			atr_loop : do kr = 1, aatom%Ntr_line
  !
  ! 				kc = aatom%at(kr)%ik
  !
  ! 				Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
  ! 				i = aatom%lines(kc)%i
  ! 				j = aatom%lines(kc)%j
  !
  ! 				Nl = Nred-dk_min+dk_max-Nblue+1
  ! 				a0 = Nblue+dk_min-1
  !
  ! 				JJ = 0.0
  ! 				wphi = 0.0
  ! 				xcc_ji = 0.0
  !
  ! ! 				enddo
  ! 				do l=1,Nl
  ! 					if (l==1) then
  ! 						wl = 0.5*(1d3*hv)
  ! 						!wl = 0.5*(lambda(a0+l+1)-lambda(a0+l)) * clight / aatom%lines(kc)%lambda0
  ! 					elseif (l==Nl) then
  ! 						wl = 0.5*(1d3*hv)
  ! 						!wl = 0.5*(lambda(a0+l)-lambda(a0+l-1)) * clight / aatom%lines(kc)%lambda0
  ! 					else
  ! 						wl = 1d3*hv
  ! 						!wl = 0.5*(lambda(a0+l+1)-lambda(a0+l-1)) * clight / aatom%lines(kc)%lambda0
  ! 					endif
  !
  ! 					JJ = JJ + wl * Ieff(a0+l) * aatom%lines(kc)%phi_loc(l,iray,id)
  ! 					wphi = wphi + wl * aatom%lines(kc)%phi_loc(l,iray,id)
  ! 					xcc_ji = xcc_ji + chi_up(a0+l,i,nact,id)*psi(a0+l,iray,id)*Uji_down(a0+l,j,nact,id)
  !
  ! 				enddo
  !
  ! 				JJ = waq * JJ / wphi
  !
  ! 				!init at Aji
  ! 				aatom%lines(kc)%Rji(id) = aatom%lines(kc)%Rji(id) + JJ * aatom%lines(kc)%Bji
  ! 				aatom%lines(kc)%Rij(id) = aatom%lines(kc)%Rij(id) + JJ * aatom%lines(kc)%Bij
  !
  ! 				!cross-coupling terms
  ! 				xcc_ij = 0.0_dp!
  ! 					sub_atr_loop : do krr=1, aatom%Ntr
  !
  ! 						kcc = aatom%at(krr)%ik
  !
  ! 						select case (aatom%at(krr)%trtype)
  !
  ! 						case ("ATOMIC_LINE")
  ! 							ip = aatom%lines(kcc)%i
  ! 							jp = aatom%lines(kcc)%j
  ! 							N1 = aatom%lines(kcc)%Nblue+dk_min
  ! 							N2 = aatom%lines(kcc)%Nred+dk_max
  ! 							if (jp==i) then !i upper level of another transitions
  ! 								xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
  ! 							endif
  !
  ! 						case ("ATOMIC_CONTINUUM")
  ! 							ip = aatom%continua(kcc)%i
  ! 							jp = aatom%continua(kcc)%j
  ! 							N1 = aatom%continua(kcc)%Nb
  ! 							N2 = aatom%continua(kcc)%Nr
  ! 							if (jp==i) then !i upper level of another transitions
  ! 								xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
  ! 							endif
  !
  ! 						case default
  ! 							call error("Transition type unknown cross-coupling (line)", aatom%at(krr)%trtype)
  ! 						end select
  !
  !
  ! 					enddo sub_atr_loop
  !
  ! 				aatom%lines(kc)%Rji(id) = aatom%lines(kc)%Rji(id) - xcc_ji * waq
  ! 				aatom%lines(kc)%Rij(id) = aatom%lines(kc)%Rij(id) + xcc_ij * waq
  !
  ! 			end do atr_loop
  !
  ! 			atrc_loop : do kr = aatom%Ntr_line+1, aatom%Ntr
  !
  ! 	 			kc = aatom%at(kr)%ik
  !
  ! 				i = aatom%continua(kc)%i; j = aatom%continua(kc)%j
  ! ! 				chi_ion = Elements(aatom%periodic_table)%ptr_elem%ionpot(aatom%stage(j))
  ! ! 				neff = aatom%stage(j) * sqrt(aatom%Rydberg / (aatom%E(j) - aatom%E(i)) )
  !
  !
  ! 				Nblue = aatom%continua(kc)%Nb; Nred = aatom%continua(kc)%Nr
  ! 				Nl = Nred-Nblue + 1
  !
  ! 				JJ = 0.0
  ! 				JJb = 0.0
  ! 				xcc_ji = 0.0
  !
  ! 				icell_d = 1
  ! 				if (ldissolve) then
  ! 					if (aatom%ID=="H") icell_d = icell
  ! 				endif
  !
  ! 				do l=1, Nl
  ! 					if (l==1) then
  ! 						wl = 0.5*(lambda(Nblue+l)-lambda(Nblue)) / lambda(Nblue)
  ! 					elseif (l==Nl) then
  ! 						wl = 0.5*(lambda(Nred)-lambda(Nred-1))  / lambda(Nred)
  ! 					else
  ! 						wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-2)) / lambda(Nblue+l-1)
  ! 					endif
  !
  ! 					a1 = aatom%continua(kc)%alpha_nu(l,icell_d)
  !
  ! 					JJ = JJ + wl  * a1 * Ieff(Nblue+l-1)
  !
  ! 					ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-1))
  ! 					twohnu3_c2 = twohc/lambda(Nblue+l-1)**3
  ! 					JJb = JJb + wl * ( Ieff(Nblue+l-1) + twohnu3_c2 ) * ehnukt * a1
  !
  ! 					xcc_ji = xcc_ji + chi_up(Nblue+l-1,i,nact,id)*Uji_down(Nblue+l-1,j,nact,id)*psi(Nblue+l-1,iray,id)
  !
  ! 				enddo
  !
  ! 				aatom%continua(kc)%Rij(id) = aatom%continua(kc)%Rij(id) + waq*fourpi_h * JJ
  ! 				aatom%continua(kc)%Rji(id) = aatom%continua(kc)%Rji(id) + waq*fourpi_h * JJb * ( aatom%nstar(i,icell)/aatom%nstar(j,icell) )
  !
  ! 				xcc_ij = 0.0_dp
  ! 					sub_atrc_loop : do krr=1, aatom%Ntr
  !
  ! 						kcc = aatom%at(krr)%ik
  !
  ! 						select case (aatom%at(krr)%trtype)
  ! 						case ("ATOMIC_LINE")
  ! 							ip = aatom%lines(kcc)%i
  ! 							jp = aatom%lines(kcc)%j
  ! 							N1 = aatom%lines(kcc)%Nblue+dk_min
  ! 							N2 = aatom%lines(kcc)%Nred+dk_max
  ! 							if (jp==i) then !i upper level of another transitions
  ! 								xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
  ! 							endif
  ! 						case ("ATOMIC_CONTINUUM")
  ! 							ip = aatom%continua(kcc)%i
  ! 							jp = aatom%continua(kcc)%j
  ! 							N1 = aatom%continua(kcc)%Nb
  ! 							N2 = aatom%continua(kcc)%Nr
  ! 							if (jp==i) then !i upper level of another transitions
  ! 								xcc_ij = xcc_ij + sum(chi_down(N1:N2,j,nact,id)*psi(N1:N2,iray,id)*Uji_down(N1:N2,i,nact,id))
  ! 							endif
  !
  ! 						case default
  ! 							call error("Transition type unknown cross-coupling (continuum)", aatom%at(krr)%trtype)
  ! 						end select
  !
  !
  ! 					enddo sub_atrc_loop
  !
  ! 				aatom%continua(kc)%Rji(id) = aatom%continua(kc)%Rji(id) - xcc_ji * waq
  ! 				aatom%continua(kc)%Rij(id) = aatom%continua(kc)%Rij(id) + xcc_ij * waq
  !
  ! 			end do atrc_loop
  !
  ! 			aatom => NULL()
  !
  ! 		end do aatom_loop
  !
  ! 	return
  ! 	end subroutine calc_rates_mali_old


  subroutine store_radiative_rates_mali(id, icell, init, waq, Nmaxtr, Rij, Rji)
    !continua radiative rates computed also with total I
    integer, intent(in) :: icell, Nmaxtr, id
    logical, intent(in) :: init
    real(kind=dp), intent(in) :: waq
    real(kind=dp), dimension(NactiveAtoms, Nmaxtr), intent(inout) :: Rij, Rji
    integer :: kc, kr, nact, l, imu,iphi,iray_p
    integer :: i, j, Nl, Nblue, Nred, a0, icell_d
    real(kind=dp) :: a1, JJb, a2, di1, di2, ehnukt, twohnu3_c2, wl
    real(kind=dp) :: chi_ion, neff, Ieff1, Ieff2, JJ, wphi, J1, J2
    real(kind=dp) :: ni_on_nj_star
    type (AtomType), pointer :: atom

    if (init) then
       Rij(:,:) = 0.0_dp
       Rji(:,:) = 0.0_dp
    endif
    
    !if lforce_lte we only store the Aji terms in the line radiative rates
    !and 0 elsewhere.
    if (lforce_lte) then
    	do nact=1, NactiveAtoms

       		atom => Activeatoms(nact)%ptr_atom

       		do kc=1, atom%Ntr_line

        		kr = atom%at(kc)%ik

        		if (init) Rji(nact,kc) = atom%lines(kr)%Aji
        	
        	end do
        	
        	atom => NULL()
        	
        enddo
             
    	return
    endif

    do nact=1, NactiveAtoms

       atom => Activeatoms(nact)%ptr_atom

       do kc=1, atom%Ntr

          kr = atom%at(kc)%ik

          select case (atom%at(kc)%trtype)

          case ("ATOMIC_LINE")

             Nred = atom%lines(kr)%Nred
             Nblue = atom%lines(kr)%Nblue
             i = atom%lines(kr)%i
             j = atom%lines(kr)%j

             Nl = Nred-dk_min+dk_max-Nblue+1
             a0 = Nblue+dk_min-1
             if (init) Rji(nact,kc) = atom%lines(kr)%Aji

             JJ = 0.0
             wphi = 0.0

             do l=1,Nl
                if (l==1) then
                   wl = 0.5*(1d3*hv)
                   !wl = 0.5*(lambda(a0+l+1)-lambda(a0+l)) * clight / atom%lines(kr)%lambda0
                elseif (l==Nl) then
                   !wl = 0.5*(lambda(a0+l)-lambda(a0+l-1)) * clight / atom%lines(kr)%lambda0
                   wl = 0.5*(1d3*hv)
                else
                   wl = 1d3*hv
                   !wl = 0.5*(lambda(a0+l+1)-lambda(a0+l-1)) * clight / atom%lines(kr)%lambda0
                endif

                JJ = JJ + wl * Itot(a0+l,1,id) * atom%lines(kr)%phi_loc(l,1,id)
                wphi = wphi + wl * atom%lines(kr)%phi_loc(l,1,id)

             enddo
             !init at Aji
             Rji(nact,kc) = Rji(nact,kc) + atom%lines(kr)%Bji * JJ * waq / wphi
             Rij(nact,kc) = Rij(nact,kc) + atom%lines(kr)%Bij * JJ * waq / wphi



          case ("ATOMIC_CONTINUUM")


             i = atom%continua(kr)%i; j = atom%continua(kr)%j
             chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
             neff = atom%stage(j) * sqrt(atom%Rydberg / (atom%E(j) - atom%E(i)) )


             Nblue = atom%continua(kr)%Nb; Nred = atom%continua(kr)%Nr
             Nl = Nred-Nblue + 1

             icell_d = 1
             if (ldissolve) then
                if (atom%ID=="H") icell_d = icell
             endif

             ! 					ni_on_nj_star = ne(icell) * phi_T(icell, atom%g(i)/atom%g(j), atom%E(j)-atom%E(i))
             ni_on_nj_star = atom%nstar(i,icell)/atom%nstar(j,icell)

             JJ = 0.0
             JJb = 0.0
             do l=1, Nl

                if (l==1) then
                   wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-1)) / lambda(Nblue+l-1)
                elseif (l==Nl) then
                   wl = 0.5*(lambda(Nblue+l-1)-lambda(Nblue+l-2))  / lambda(Nblue+l-1)
                else
                   wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-2)) / lambda(Nblue+l-1)
                endif

                a1 = atom%continua(kr)%alpha_nu(l,icell_d)

                JJ = JJ + wl * Itot(Nblue+l-1,1,id) * a1
                JJb = JJb + wl * a1 * exp(-hc_k/T(icell)/lambda(Nblue+l-1)) * &
                     (twohc/lambda(Nblue+l-1)**3 + Itot(Nblue+l-1,1,id))
             enddo

             Rij(nact,kc) = Rij(nact,kc) + waq*fourpi_h * JJ
             Rji(nact,kc) = Rji(nact,kc) + waq*fourpi_h * JJb * ni_on_nj_star

          case default
             call error("transition type unknown", atom%at(l)%trtype)

          end select

       enddo

       atom => NULL()
    enddo

    return
  end subroutine store_radiative_rates_mali

  !used only to write radiative rates computed with the converged radiation field
  !not used in the iterative scheme
  ! 	subroutine store_radiative_rates(id, icell, n_rayons, Nmaxtr, Rij, Rji, Jr, angle_quad)
  ! 	!continua radiative rates computed also with total I
  ! 		integer, intent(in) :: icell, n_rayons, Nmaxtr, id
  ! 		logical, intent(in) :: angle_quad
  ! 		real(kind=dp), dimension(NactiveAtoms, Nmaxtr), intent(inout) :: Rij, Rji
  ! 		real(kind=dp), dimension(Nlambda) :: Jr
  ! 		integer :: kc, kr, nact, l, iray,imu,iphi,iray_p
  ! 		integer :: i, j, Nl, Nblue, Nred, a0, icell_d
  ! 		real(kind=dp) :: a1, JJb, a2, di1, di2, ehnukt, twohnu3_c2, wl
  ! 		real(kind=dp) :: chi_ion, neff, Ieff1, Ieff2, JJ, wphi, J1, J2
  ! 		type (AtomType), pointer :: atom
  !
  ! 		Rij(:,:) = 0.0_dp
  ! 		Rji(:,:) = 0.0_dp
  !
  ! 		do nact=1, NactiveAtoms
  !
  ! 			atom => Activeatoms(nact)%ptr_atom
  !
  ! 			do kc=1, atom%Ntr
  !
  ! 				kr = atom%at(kc)%ik
  !
  ! 				select case (atom%at(kc)%trtype)
  !
  ! 				case ("ATOMIC_LINE")
  !
  ! 					Nred = atom%lines(kr)%Nred
  ! 					Nblue = atom%lines(kr)%Nblue
  ! 					i = atom%lines(kr)%i
  ! 					j = atom%lines(kr)%j
  !
  ! 					Nl = Nred-dk_min+dk_max-Nblue+1
  ! 					a0 = Nblue+dk_min-1
  ! 					Rji(nact,kc) = atom%lines(kr)%Aji
  !
  !
  ! 					if (angle_quad) then
  ! 						write(*,*) " Angle_quad not implemented yet"
  ! 						!check iray if lmali
  ! 						stop " mali"
  !
  ! 						else !pure mc
  ! 							do iray=1, n_rayons
  ! 								JJ = 0.0
  ! 								wphi = 0.0
  !
  ! 								do l=1,Nl
  !
  ! 									if (l==1) then
  ! 										wl = 0.5*(lambda(a0+l+1)-lambda(a0+l)) * clight / atom%lines(kr)%lambda0
  ! 									elseif (l==Nl) then
  ! 										wl = 0.5*(lambda(a0+l)-lambda(a0+l-1)) * clight / atom%lines(kr)%lambda0
  ! 									else
  ! 										wl = 0.5*(lambda(a0+l+1)-lambda(a0+l-1)) * clight / atom%lines(kr)%lambda0
  ! 									endif
  !
  ! 									JJ = JJ + wl * Itot(a0+l,iray,id) * atom%lines(kr)%phi_loc(l,iray,id)
  ! 									wphi = wphi + wl * atom%lines(kr)%phi_loc(l,iray,id)
  !
  ! 								enddo
  ! 								!init at Aji
  ! 								Rji(nact,kc) = Rji(nact,kc) + atom%lines(kr)%Bji * JJ / n_rayons / wphi
  ! 								Rij(nact,kc) = Rij(nact,kc) + atom%lines(kr)%Bij * JJ / n_rayons / wphi
  !
  ! 							enddo
  !
  ! 						endif !over angle_quad
  !
  !
  ! 				case ("ATOMIC_CONTINUUM")
  !
  !
  ! 					i = atom%continua(kr)%i; j = atom%continua(kr)%j
  ! 					chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
  ! 					neff = atom%stage(j) * sqrt(atom%Rydberg / (atom%E(j) - atom%E(i)) )
  !
  !
  ! 					Nblue = atom%continua(kr)%Nb; Nred = atom%continua(kr)%Nr
  ! 					Nl = Nred-Nblue + 1
  !
  ! 					icell_d = 1
  ! 					if (ldissolve) then
  ! 						if (atom%ID=="H") icell_d = icell
  ! 					endif
  !
  ! 					JJ = 0.0
  ! 					JJb = 0.0
  !
  ! 					do l=1, Nl
  !
  ! 						if (l==1) then
  ! 							wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-1)) / lambda(Nblue+l-1)
  ! 						elseif (l==Nl) then
  ! 							wl = 0.5*(lambda(Nblue+l-1)-lambda(Nblue+l-2))  / lambda(Nblue+l-1)
  ! 						else
  ! 							wl = 0.5*(lambda(Nblue+l)-lambda(Nblue+l-2)) / lambda(Nblue+l-1)
  ! 						endif
  !
  ! 						a1 = atom%continua(kr)%alpha_nu(l,icell_d)
  !
  ! 						JJ = JJ + wl * Jr (Nblue+l-1) * a1
  ! 						JJb = JJb + wl * a1 * exp(-hc_k/T(icell)/lambda(Nblue+l-1)) * &
  ! 							(twohc/lambda(Nblue+l-1)**3 + Jr(Nblue+l-1))
  ! 					enddo
  !
  !
  ! 					Rij(nact,kc) = fourpi_h * JJ
  ! 					Rji(nact,kc) = fourpi_h * JJb * atom%nstar(i,icell)/atom%nstar(j,icell)
  !
  ! 				case default
  ! 					call error("transition type unknown", atom%at(l)%trtype)
  !
  ! 				end select
  !
  ! 			enddo
  !
  ! 			atom => NULL()
  ! 		enddo
  !
  ! 	return
  ! 	end subroutine store_radiative_rates

  subroutine store_rate_matrices(id, icell, Nmaxlevel, Aij)
    !atom%Gamma should be filled with the value before exciting subiterations
	!built rate matrix with allocated radiative rate
	!These are the MALI rates at the moment (not Rij_all and Rji_all as it should rather be.)
    integer, intent(in) :: icell, Nmaxlevel, id
    real(kind=dp), dimension(NactiveAtoms, Nmaxlevel, Nmaxlevel), intent(inout) :: Aij
    integer :: nact
    type (AtomType), pointer :: atom

    Aij(:,:,:) = 0.0_dp

    do nact=1, NactiveAtoms

       atom => Activeatoms(nact)%ptr_atom
       !!call collision_matrix_atom(id, icell, atom)
       call initGamma_atom(id, atom) !init with collisions
       !if lforce_lte it will be equivalent to the collision matrix!
       !lforce_lte is set to .false. here so that, even at LTE, Gamma includes Aji terms
       !in the radiative rates and the collision matrix is the pure Gamma matrix!
       !In that case Rij_all contains only Aji
       call rate_matrix_atom(id, icell, atom, .false.)
       call eliminate_delta(id, atom)
       Aij(nact,1:atom%Nlevel,1:atom%Nlevel) = atom%Gamma(:,:,id)
       atom => NULL()

    enddo

    return
  end subroutine store_rate_matrices


  subroutine calc_bb_rates_hoger(id, icell, iray, n_rayons)

    integer, intent(in) :: id, icell, iray, n_rayons
    real(kind=dp) :: chicont, etacont
    real(kind=dp), dimension(Nlambda) :: chi_tot, eta_tot
    integer :: nact, Nred, Nblue, kc, kr, i, j, l, a0, Nl
    real(kind=dp) :: etau, wphi, JJ, JJb, j1, j2, Ieff1, Ieff2
    real(kind=dp) :: dissolve, chi_ion, neff, wl
    type(AtomType), pointer :: aatom, atom

    real(kind=dp) :: wi, wj, nn, dk0

    !only with lines, continuum added during integration
    chi_tot(:) = chi0_bb(:,icell)!0.0
    eta_tot(:) = eta0_bb(:,icell)!0.0
    if (lelectron_scattering) then
    	write(*,*) "elec scatt not included in calc_bb_rates_hoger!"
    	stop
!        eta_tot(:) = Jnu(:,id) * thomson(icell)
    endif

    atom_loop : do nact = 1, Natom!Nactiveatoms
       atom => Atoms(nact)%ptr_atom!ActiveAtoms(nact)%ptr_atom

       tr_loop : do kr = 1, atom%Ntr_line


          kc = atom%at(kr)%ik !relative index of a transition among continua or lines

          Nred = atom%lines(kc)%Nred;Nblue = atom%lines(kc)%Nblue
          i = atom%lines(kc)%i;j = atom%lines(kc)%j

          wj = 1.0; wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn
                wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1),hydrogen%n(1,icell)) !1 for H
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1),hydrogen%n(1,icell))
             endif
          endif

          if ((atom%n(i,icell)*wj/wi - atom%n(j,icell)*atom%lines(kc)%gij) > 0.0_dp) then


             chi_tot(Nblue+dk_min:Nred+dk_max) = chi_tot(Nblue+dk_min:Nred+dk_max) + &
                  hc_fourPI * atom%lines(kc)%Bij * atom%lines(kc)%phi_loc(:,iray,id) * &
                  (atom%n(i,icell)*wj/wi - atom%lines(kc)%gij*atom%n(j,icell))

             eta_tot(Nblue+dk_min:Nred+dk_max)= eta_tot(Nblue+dk_min:Nred+dk_max) + &
                  atom%lines(kc)%twohnu3_c2 * atom%lines(kc)%gij * hc_fourPI * atom%lines(kc)%Bij * &
                  atom%lines(kc)%phi_loc(:,iray,id) * atom%n(j,icell)


          else !neg or null

             ! 					eta_tot(Nblue+dk_min:Nred+dk_max)= eta_tot(Nblue+dk_min:Nred+dk_max) + &
             ! 						atom%lines(kc)%twohnu3_c2 * atom%lines(kc)%gij * hc_fourPI * atom%lines(kc)%Bij * atom%lines(kc)%phi_loc(:,iray,id) * atom%n(j,icell)
             eta_tot(Nblue+dk_min:Nred+dk_max)= eta_tot(Nblue+dk_min:Nred+dk_max) + 0.0_dp
             chi_tot(Nblue+dk_min:Nred+dk_max) = chi_tot(Nblue+dk_min:Nred+dk_max) - &
                  hc_fourPI * atom%lines(kc)%Bij * atom%lines(kc)%phi_loc(:,iray,id) * &
                  (atom%n(i,icell)*wj/wi - atom%lines(kc)%gij*atom%n(j,icell))


          endif


       end do tr_loop

       atom => NULL()

    end do atom_loop


    !now start integrates bound rates
    aatom_loop : do nact = 1, Nactiveatoms
       aatom => ActiveAtoms(nact)%ptr_atom

       atr_loop : do kr = 1, aatom%Ntr_line

          !relative index of a transition among continua or lines
          kc = aatom%at(kr)%ik

          Nred = aatom%lines(kc)%Nred;Nblue = aatom%lines(kc)%Nblue
          i = aatom%lines(kc)%i
          j = aatom%lines(kc)%j

          Nl = Nred-dk_min+dk_max-Nblue+1
          a0 = Nblue+dk_min-1

          JJ = 0.0
          wphi = 0.0

          !!continuum around this line
          !!suppose that in the region where the line lies (which can be far from lambda0 if large velocity)
          !!the continuum is constant. Might not be true
          ! 				if (iray==1) then
          ! 					!which is fastest / accurate ? eta_es already added if any
          ! 					!call continuum_line (icell, aatom%lines(kc)%lambda0,  aatom%lines(kc)%chi0,  aatom%lines(kc)%eta0)
          ! 					aatom%lines(kc)%chi0 = interp(chi_c(:,icell)+chi_c_nlte(:,icell), lambda_cont, aatom%lines(kc)%lambda0)
          ! 					aatom%lines(kc)%eta0 = interp(eta_c(:,icell)+eta_c_nlte(:,icell), lambda_cont, aatom%lines(kc)%lambda0)
          ! 				endif

          do l=2, Nl
             !
             wphi = wphi + 0.5 * (aatom%lines(kc)%phi_loc(l,iray,id)+aatom%lines(kc)%phi_loc(l-1,iray,id)) * 1e3*hv

             ! 					etau = exp(-(chi_tot(a0+l)+aatom%lines(kc)%chi0)*ds(iray,id))
             ! 					Ieff1 = Itot(a0+l,iray,id) * etau + ( (eta_tot(a0+l)+aatom%lines(kc)%eta0)/(chi_tot(a0+l)+aatom%lines(kc)%chi0 + tiny_dp) ) * (1.0 - etau)
             !
             ! 					etau = exp(-(chi_tot(a0+l-1)+aatom%lines(kc)%chi0)*ds(iray,id))
             ! 					Ieff2 = Itot(a0+l-1,iray,id) * etau + ( (eta_tot(a0+l-1)+aatom%lines(kc)%eta0)/(chi_tot(a0+l-1)+aatom%lines(kc)%chi0+tiny_dp) ) * (1.0 - etau)

             etau = exp(-chi_tot(a0+l)*ds(iray,id))
             Ieff1 = Itot(a0+l,iray,id) * etau +  (eta_tot(a0+l) / (chi_tot(a0+l) + tiny_dp) ) * (1.0 - etau)

             etau = exp(-chi_tot(a0+l-1)*ds(iray,id))
             Ieff2 = Itot(a0+l-1,iray,id) * etau + (eta_tot(a0+l-1) /(chi_tot(a0+l-1)+tiny_dp)) * (1.0 - etau)

             J1 = Ieff1 * aatom%lines(kc)%phi_loc(l,iray,id)
             J2 = Ieff2 * aatom%lines(kc)%phi_loc(l-1,iray,id)
             !
             JJ = JJ + 0.5 * (J1 + J2) * hv * 1d3


          enddo

          JJ = JJ / n_rayons / wphi
          !
          ! 				Renormalise profile ?
          ! 				if ((wphi < 0.95) .or. (wphi > 1.05)) then
          ! 					write(*,*) "icell = ", icell, " id = ", id
          ! 					write(*,*) " --> Beware, profile not well normalized for line ", i, j, " area = ", wphi
          ! 				endif

          !init at Aji
          aatom%lines(kc)%Rji(id) = aatom%lines(kc)%Rji(id) + JJ * aatom%lines(kc)%Bji
          aatom%lines(kc)%Rij(id) = aatom%lines(kc)%Rij(id) + JJ * aatom%lines(kc)%Bij

          if (any_nan_infinity_vector(aatom%lines(kc)%Rij) /= 0.0) then
             write(*,*) "icell=",icell," id=",id," nb_proc=", nb_proc, aatom%lines(kc)%Rij(id)
             write(*,*) "Rij=", aatom%lines(kc)%Rij
             write(*,*) "I=",Itot(:,iray,id)
             write(*,*) "chi0_bb=",chi0_bb(:,icell)
             write(*,*) "eta0_bb=",eta0_bb(:,icell)
             write(*,*) "phi(iray)=", aatom%lines(kc)%phi_loc(:,iray,id)
             stop
          endif
          if (any_nan_infinity_vector(aatom%lines(kc)%Rji) /= 0.0) then
             write(*,*) "icell=",icell," id=",id," nb_proc=", nb_proc, aatom%lines(kc)%Rji(id)
             write(*,*) "Rji=",aatom%lines(kc)%Rji
             write(*,*) "I=",Itot(:,iray,id)
             write(*,*) "chi0_bb=",chi0_bb(:,icell)
             write(*,*) "eta0_bb=",eta0_bb(:,icell)
             write(*,*) "phi(iray)=", aatom%lines(kc)%phi_loc(:,iray,id)
             stop
          endif

       end do atr_loop

       aatom => NULL()

    end do aatom_loop

    return
  end subroutine calc_bb_rates_hoger

  SUBROUTINE calc_bf_rates(id, icell, Jr)
    integer, intent(in) :: id, icell
    integer :: nact
    real(kind=dp), dimension(Nlambda_cont), intent(in) :: Jr

    do nact=1, NactiveAtoms

       ! 			call calc_bf_rates_atom(id, icell, ActiveAtoms(nact)%ptr_atom, Jr)
       call calc_bf_rates_atom_2(id, icell, ActiveAtoms(nact)%ptr_atom, Jr)

    enddo

    RETURN
  END SUBROUTINE calc_bf_rates

  !only depends on Jr(size Jnu_cont), therefore no ray integ
  SUBROUTINE calc_bf_rates_atom(id, icell, atom, Jr)

    integer, intent(in) :: icell, id
    type (AtomType), intent(inout) :: atom
    real(kind=dp), dimension(Nlambda_cont), intent(in) :: Jr
    integer :: kr, kc, i, j, l, Nblue, Nred, Nl, iray
    real(kind=dp) :: Ieff1, Ieff2, wphi, JJ, JJb, j1, j2, twohnu3_c2
    real(kind=dp) :: di1, di2, chi_ion, neff, wl, ehnukt, test_integ

    tr_loop : do kr=atom%Ntr_line+1, atom%Ntr
       kc = atom%at(kr)%ik


       i = atom%continua(kc)%i; j = atom%continua(kc)%j
       chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
       neff = atom%stage(j) * sqrt(atom%Rydberg / (atom%E(j) - atom%E(i)) )


       Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
       Nl = atom%continua(kc)%Nlambda

       !not useful
       atom%continua(kc)%Rij(id) = 0.0
       atom%continua(kc)%Rji(id) = 0.0

       JJ = 0.0
       JJb = 0.0

       !test_integ = 0
       do l=2, Nl
          wl = lambda_cont(Nblue+l-1) - lambda_cont(Nblue+l-2)
          di1 = D_i(icell, neff, real(atom%stage(i)), 1.0, lambda_cont(Nblue+l-1), atom%continua(kc)%lambda0, chi_ion)
          di2 = D_i(icell, neff, real(atom%stage(i)), 1.0, lambda_cont(Nblue+l-2), atom%continua(kc)%lambda0, chi_ion)

          ! 					J1 = Jnu_cont(Nblue+l-1,icell) * atom%continua(kc)%alpha(l)*wl/lambda_cont(Nblue+l-1)
          ! 					J2 = Jnu_cont(Nblue+l-2,icell) * atom%continua(kc)%alpha(l-1)*wl/lambda_cont(Nblue+l-2)
          J1 = Jr(Nblue+l-1) * di1 * atom%continua(kc)%alpha(l)/lambda_cont(Nblue+l-1)
          J2 = Jr(Nblue+l-2) * di2 * atom%continua(kc)%alpha(l-1)/lambda_cont(Nblue+l-2)

          JJ = JJ + 0.5 * (J1 + J2) * wl

          ehnukt = exp(-hc_k/T(icell)/lambda_cont(Nblue+l-1))
          twohnu3_c2 = twohc/lambda_cont(Nblue+l-1)**3
          ! 					J1 = ( Jnu_cont(Nblue+l-1,icell) + atom%continua(kc)%twohnu3_c2(l)) * ehnukt * atom%continua(kc)%alpha(l) * wl/lambda_cont(Nblue+l-1)
          J1 = ( Jr(Nblue+l-1) + twohnu3_c2 ) * ehnukt * di1 * atom%continua(kc)%alpha(l) / lambda_cont(Nblue+l-1)
          ehnukt = exp(-hc_k/T(icell)/lambda_cont(Nblue+l-2))
          twohnu3_c2 = twohc/lambda_cont(Nblue+l-2)**3
          ! 					J2 = ( Jnu_cont(Nblue+l-2,icell) + atom%continua(kc)%twohnu3_c2(l-1)) * ehnukt * atom%continua(kc)%alpha(l-1) * wl/lambda_cont(Nblue+l-2)
          J2 = ( Jr(Nblue+l-2) + twohnu3_c2 ) * ehnukt * di2 * atom%continua(kc)%alpha(l-1) / lambda_cont(Nblue+l-2)

          JJb = JJb + 0.5 * (J1 + J2) * wl
          !test_integ = test_integ + 0.5 * (	exp(-hc_k/T(icell)/lambda_cont(Nblue+l-1)) + exp(-hc_k/T(icell)/lambda_cont(Nblue+l-2))) * clight * (1.0/lambda_cont(Nblue+l-2) - 1.0/lambda_cont(Nblue+l-1)) * 1e9
       enddo
       !write(*,*) test_integ, KBOLTZMANN*T(icell)/Hplanck * (exp(-hc_k/T(icell)/lambda_cont(Nl+Nblue-1)) - exp(-hc_k/T(icell)/lambda_cont(Nblue)))
       !stop
       atom%continua(kc)%Rij(id) = fourpi_h * JJ
       atom%continua(kc)%Rji(id) = fourpi_h * JJb * atom%nstar(i,icell)/atom%nstar(j,icell)

    end do tr_loop


    RETURN
  END SUBROUTINE calc_bf_rates_atom

  !Jr size Itot(:,1,1)
  SUBROUTINE calc_bf_rates_atom_2(id, icell, atom, Jr)

    integer, intent(in) :: icell, id
    type (AtomType), intent(inout) :: atom
    real(kind=dp), dimension(Nlambda), intent(in) :: Jr
    integer :: kr, kc, i, j, l, Nblue, Nred, Nl, iray
    real(kind=dp) :: Ieff1, Ieff2, wphi, JJ, JJb, j1, j2, twohnu3_c2, a1, a2
    real(kind=dp) :: di1, di2, chi_ion, neff, wl, ehnukt, test_integ

    tr_loop : do kr=atom%Ntr_line+1, atom%Ntr
       kc = atom%at(kr)%ik


       i = atom%continua(kc)%i; j = atom%continua(kc)%j
       chi_ion = Elements(atom%periodic_table)%ptr_elem%ionpot(atom%stage(j))
       neff = atom%stage(j) * sqrt(atom%Rydberg / (atom%E(j) - atom%E(i)) )


       Nblue = atom%continua(kc)%Nb; Nred = atom%continua(kc)%Nr
       Nl = Nred-Nblue + 1


       JJ = 0.0
       JJb = 0.0

       do l=2, Nl
          wl = lambda(Nblue+l-1) - lambda(Nblue+l-2)
          di1 = D_i(icell, neff, real(atom%stage(i)), 1.0, lambda(Nblue+l-1), atom%continua(kc)%lambda0, chi_ion)
          di2 = D_i(icell, neff, real(atom%stage(i)), 1.0, lambda(Nblue+l-2), atom%continua(kc)%lambda0, chi_ion)
          a1 = interp(atom%continua(kc)%alpha,lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), lambda(Nblue+l-1))
          a2 = interp(atom%continua(kc)%alpha,lambda_cont(atom%continua(kc)%Nblue:atom%continua(kc)%Nred), lambda(Nblue+l-2))

          J1 = Jr(Nblue+l-1) * di1 * a1 / lambda(Nblue+l-1)
          J2 = Jr(Nblue+l-2) * di2 * a2 / lambda(Nblue+l-2)

          JJ = JJ + 0.5 * (J1 + J2) * wl

          ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-1))
          twohnu3_c2 = twohc/lambda(Nblue+l-1)**3
          J1 = ( Jr(Nblue+l-1) + twohnu3_c2 ) * ehnukt * di1 * a1 / lambda(Nblue+l-1)

          ehnukt = exp(-hc_k/T(icell)/lambda(Nblue+l-2))
          twohnu3_c2 = twohc/lambda(Nblue+l-2)**3
          J2 = ( Jr(Nblue+l-2) + twohnu3_c2 ) * ehnukt * di2 * a2 / lambda(Nblue+l-2)

          JJb = JJb + 0.5 * (J1 + J2) * wl
          ! 					if (i==1 .and. j==6) then
          ! 						write(*,*) "lamcont1=", lambda_cont(atom%continua(kc)%Nblue), lambda_cont(atom%continua(kc)%Nred)
          ! 						write(*,*) atom%continua(kc)%alpha(1), atom%continua(kc)%alpha(atom%continua(kc)%Nlambda)
          ! 						write(*,*) l, lambda(Nblue+l-1), lambda(Nblue+l-2)
          ! 						write(*,*) "alpha=", a1, " alpha2=",a2, " diss1&2=",di1, di2
          ! 						write(*,*) "J1&2=",Jr(Nblue+l-1), Jr(Nblue+l-2)
          ! 					endif
       enddo

       atom%continua(kc)%Rij(id) = fourpi_h * JJ
       atom%continua(kc)%Rji(id) = fourpi_h * JJb * atom%nstar(i,icell)/atom%nstar(j,icell)

    end do tr_loop


    RETURN
  END SUBROUTINE calc_bf_rates_atom_2

  !occupa prob
  SUBROUTINE rate_matrix_atom_old(id, icell, atom, switch_to_lte)
    !Gamma init to collision rates
    integer, intent(in) :: icell, id
    type (AtomType), intent(inout) :: atom
    logical, optional, intent(in) :: switch_to_lte
    integer :: kr, kc, i, j, l, Nblue, Nred, Nl
    real(kind=dp) :: wj, wi, neff

    if (present(switch_to_lte)) then
       if (switch_to_lte) return
    end if

    tr_loop : do kr=1, atom%Ntr
       kc = atom%at(kr)%ik

       SELECT CASE (atom%at(kr)%trtype)

       CASE ("ATOMIC_LINE")
          i = atom%lines(kc)%i; j = atom%lines(kc)%j
          wj = 1.0
          wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn
                !neff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(find_continuum(atom,j)) - atom%E(i)) )
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1),hydrogen%n(1,icell))
                wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1),hydrogen%n(1,icell))

             endif
          endif

          atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + atom%lines(kc)%Rji(id)
          atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + atom%lines(kc)%Rij(id) * wj/wi
          !!write(*,*) "l", icell, id, i, j, " Rij=", atom%lines(kc)%Rij(id), " Rji=",atom%lines(kc)%Rji(id)



       CASE ("ATOMIC_CONTINUUM")
          i = atom%continua(kc)%i; j = atom%continua(kc)%j
          wj = 1.0
          wi = 1.0
          if (ldissolve) then
             if (atom%ID=="H") then
                !nn
                !neff = (atom%stage(i)+1) * sqrt(atom%Rydberg / (atom%E(find_continuum(atom,j)) - atom%E(i)) )
                wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1),hydrogen%n(1,icell))
             endif
          endif


          atom%Gamma(i,j,id) = atom%Gamma(i,j,id) + atom%continua(kc)%Rij(id) * wj/wi
          atom%Gamma(j,i,id) = atom%Gamma(j,i,id) + atom%continua(kc)%Rji(id)
          !!write(*,*) "c", icell, id, i, j, " Rij=", atom%continua(kc)%Rij(id), " Rji=",atom%continua(kc)%Rji(id)

       CASE DEFAULT

          CALL Error("Unkown transition type", atom%at(kr)%trtype)

       END SELECT

    end do tr_loop


    RETURN
  END SUBROUTINE rate_matrix_atom_old


  SUBROUTINE remove_delta_old(id, atom)
    !Remove therm in delta(l,l')
    integer, intent(in) :: id
    type (AtomType), intent(inout) :: atom
    integer :: l

    do l = 1, atom%Nlevel
       atom%Gamma(l,l,id) = 0.0_dp
       atom%Gamma(l,l,id) = -sum(atom%Gamma(l,:,id))
    end do

    RETURN
  END SUBROUTINE remove_delta_old

  subroutine see_atom_old(id, icell,atom, dM)
    ! --------------------------------------------------------------------!
    ! For atom atom solves for the Statistical Equilibrium Equations (SEE)
    !
    ! We solve for :
    !  Sum_l' Gamma_l'l n_l' = 0 (m^-3 s^-1)
    !
    ! which is equivalent in matrix notation to:
    !
    ! GG(l,l') dot n_l' = 0 with l' is on the columns and l on the rows
    !
    ! In particular for a 2-level atom the two equations are:
    !
    ! n1*G_11 + n2*G_21 = 0
    ! n1*G12 + n2*G22 = 0,
    ! with one of these equations has to be replaced by
    ! n1 + n2 = N
    !
    ! For numerical stability, the row with the largest populations is
    ! removed (i.e., it is replaced by the mass conservation law).
    ! --------------------------------------------------------------------!

    integer, intent(in) :: icell, id
    type(AtomType), intent(inout) :: atom
    integer :: lp, imaxpop, l
    real(kind=dp), intent(out) :: dM
    real(kind=dp), dimension(atom%Nlevel) :: ndag
    real(kind=dp) :: n0
    real(kind=dp), dimension(atom%Nlevel, atom%Nlevel) :: Aij

    n0 = 0.0_dp
    !-> not yet
    !!if (atom%ID=="H") n0 = nHmin(icell)
    !!if (n0 >= ntotal_atom(icell,atom)) call error("Error n0")

    CALL remove_delta_old(id, atom)


    ndag(:) = atom%n(:,icell) / ( ntotal_atom(icell,atom) - n0 )
    atom%n(:,icell) = 0d0

    imaxpop = locate(atom%n(:,icell), maxval(atom%n(:,icell)))
    !imaxpop = atom%Nlevel

    atom%n(imaxpop,icell) = 1d0

    !Sum_l'_imaxpop * n_l' = N
    !(G11 G21)  (n1)  (0)
    !(       ) .(  ) =( )
    !(1    1 )  (n2)  (N)


    atom%Gamma(:,imaxpop,id) = 1d0 !all columns of the last row for instance
    Aij = transpose(atom%Gamma(:,:,id))

    if ((any_nan_infinity_matrix(atom%Gamma(:,:,id))>0)) then
       write(*,*) "BUG Gamma", id, icell
       write(*,*) atom%Gamma(:,:,id)
       write(*,*) "n = ", atom%n(:,icell)
       write(*,*) "ndag=", ndag
       stop
    end if

    CALL GaussSlv(Aij, atom%n(:,icell),atom%Nlevel)

    if ((any_nan_infinity_vector(atom%n(:,icell))>0)) then
       write(*,*) "BUG pops", id, icell, atom%n(:,icell)
       write(*,*) "ndag", ndag
       write(*,*) atom%Gamma(:,:,id)
       stop
    end if

    dM = 0.0_dp
    do l=1,atom%Nlevel

       if (atom%n(l,icell) <= prec_pops) then !relative to ntotal
          atom%n(l,icell) = 0.0_dp
       else !denormalise
          dM = max(dM, abs(atom%n(l,icell) - ndag(l))/atom%n(l,icell))
          atom%n(l,icell) = atom%n(l,icell) * ( ntotal_atom(icell,atom) - n0 )
       endif


    enddo

    return
  end subroutine see_atom_old


END MODULE statequil_atoms
