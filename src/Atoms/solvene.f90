! ------------------------------------------------------------------- !
! ------------------------------------------------------------------- !
! Module that solves for electron density for given
! temperature grid, Hydrogen total populations
! elements and their abundance and their LTE or NLTE populations.
!
!
! If ne_intial_slution is "NEMODEL" then first guess is computed using
! the value stored in ne as first guess.
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

  use atmos_type, only : Nelem, Elements, Atoms, Natom, nHtot, ne, T, icompute_atomRT, nTotal_atom
  use atom_type, only : Element, AtomType
  use math, only : interp1D
  use constant
  use lte
  use messages, only : Error, Warning
  use input, only : nb_proc
  use mcfost_env, only : dp
  use constantes, only : tiny_dp, huge_dp
  use parametres, only : n_cells
  !$ use omp_lib

  IMPLICIT NONE

  integer, parameter :: N_negative_ions = 0 !currently only H- at index 0
  real(kind=dp), parameter :: MAX_ELECTRON_ERROR=1d-6
  integer, parameter :: N_MAX_ELECTRON_ITERATIONS=50!100!50
  integer, parameter :: N_MAX_ELEMENT=26 !100 is the max, corresponding the atomic number
  real(kind=dp), parameter :: ne_min_limit = 1d-100!if ne < ne_min_limt, ne = 0
  real(kind=dp), parameter :: T_lim_non_lte = 5000.0_dp

  integer :: Max_ionisation_stage
  real(kind=dp) :: min_f_HII, max_f_HII

  !Introducing a lock for ne iterations with icompute_atomRT == 2
  !the transfer (and opacity) is solved for icompute_atomRT > 0
  !but if icompute_atomRT==2 electron density has the init value.
  !then electron iterations skipped for icompute_atomRT /= 1


CONTAINS

  function get_max_nstage()
    integer :: get_max_nstage
    integer :: ell

    get_max_nstage = 0
    do ell=1, Nelem

       get_max_nstage = max(get_max_nstage, elements(ell)%ptr_elem%Nstage)

    enddo

    write(*,*) " Maximum number of ionisation stages among all elements: ", get_max_nstage

    return
  end function get_max_nstage

  ! ----------------------------------------------------------------------!
  ! do not forget to use the good value of parition function and potential
  ! parition function read in logarithm and used is 10**(U)
  ! potential in cm-1, converted in the routines in J.
  ! ----------------------------------------------------------------------!

  SUBROUTINE ne_Hionisation0 (k, U0, U1, ne)
    ! ----------------------------------------------------------------------!
    ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
    ! of H.
    ! Njl = Nj1l * ne * phi_jl
    ! ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
    ! eq. 17.77 and 17.78
    ! ne(H) = (sqrt(N*phiH + 1)-1)/phiH
    ! ----------------------------------------------------------------------!

    integer, intent(in) :: k
    real(kind=dp), intent(in) :: U0, U1
    real(kind=dp) :: phiH
    real(kind=dp), intent(out) :: ne
    phiH = phi_jl(k, U0, U1, Elements(1)%ptr_elem%ionpot(1))

    !if no free e-, Nt = NH + NH+ with NH+=ne
    if (phiH>0) then
       ne = (sqrt(nHtot(k)*phiH*4. + 1)-1)/(2.*phiH) !without H minus
    else
       ne = 0d0
    endif

    RETURN
  END SUBROUTINE ne_Hionisation0

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
    real(kind=dp), intent(in) :: U0, U1
    real(kind=dp) :: phiH
    real(kind=dp), intent(out) :: ne

    phiH = phi_jl(k, U0, U1, Elements(1)%ptr_elem%ionpot(1))

    !if no free e-, Nt = NH + NH+ with NH+=ne
    !ne = (sqrt(nHtot(k)*phiH*4. + 1)-1)/(2.*phiH) !without H minus
    !if free e-, Nt=Nhot + NH+ + ne
    if (phiH>0) then
       ne = (sqrt(nHtot(k)*phiH + 1)-1)/(phiH)
    else
       ne = 0d0
    endif

    RETURN
  END SUBROUTINE ne_Hionisation

  SUBROUTINE ne_Metal(k, U0, U1, chi, A, ne)
    ! ----------------------------------------------------------------------!
    ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
    ! of a single metal.
    ! ----------------------------------------------------------------------!

    integer, intent(in) :: k
    real(kind=dp), intent(in) :: U0, U1, chi, A
    real(kind=dp) :: phiM, alphaM
    real(kind=dp), intent(out) :: ne

    phiM = phi_jl(k, U0, U1, chi)
    alphaM = A ! relative to H, for instance 1.-6 etc

    if (phiM > 0) then
       ne = (sqrt(alphaM*nHtot(k)*phiM +0.25*(1+alphaM)**2) - 0.5*(1+alphaM) ) / phiM
    else
       ne = 0d0
    endif

    RETURN
  END SUBROUTINE ne_Metal


  !  FUNCTION getPartitionFunctionk(elem, stage, k) result (Uk)
  !  ! ----------------------------------------------------------------------!
  !   ! Interpolate the partition function of Element elem in ionisation stage
  !   ! stage at cell point k.
  !  ! ----------------------------------------------------------------------!
  !
  !   type(Element) :: elem
  !   integer, intent(in) :: stage, k
  !   real(kind=dp) :: Uk, part_func(Npf)
  !
  !   part_func = elem%pf(stage,:)
  !   Uk = Interp1D(Tpf,part_func,T(k))
  !   !!Uk = (10.d0)**(Uk) !! 29/12/2019 changed log10 in log
  !   Uk = exp(Uk)
  !  RETURN
  ! END FUNCTION getPartitionFunctionk

  !should be min max for all cells.
  subroutine show_electron_given_per_elem(id, k, max_fjk)
    integer, intent(in) :: k, id
    real(kind=dp), intent(in), dimension(-N_negative_ions:N_MAX_ELEMENT) :: max_fjk
    integer :: ell
    integer, dimension(-N_negative_ions:9) :: list_elem!only six elements

    list_elem(1) = 1 !H
    list_elem(2) = 2 !He
    list_elem(3) = 3 !Li
    list_elem(4) = 6 !C
    list_elem(5) = 11 !Na
    list_elem(6) = 12 !Mg
    list_elem(7) = 13 !Al
    list_elem(8) = 19 !K
    list_elem(9) = 20 !Ca

    !locally
    if (k > 0 .and. k < n_cells + 1) then
       write(*,*) " -- Electron given for that cell --"
       write(*,'("  -- cell #"(1I5)," id#"(1I2)" T="(1F14.4)" K, nHtot="(1ES17.8E3)" m-3")') k, id, T(k), nHtot(k)
    else!for all
       write(*,*) " -- Max electron given for all cells --"
    endif

    do ell=-N_negative_ions, 9

       if (ell==0) then

          write(*,'("  >> H- fjk="(1ES17.8E3))') max_fjk(ell)
          !(1F18.1)" (x10^10)
       else

          write(*,'("  >> "(1A)" fjk="(1ES17.8E3))') elements(list_elem(ell))%ptr_elem%id, max_fjk(list_elem(ell))

       endif
    enddo


    !not working properly in term of display
    ! 	write(*,'(" Element          :  ", (A,1x),*(A,10x))') "H-", (elements(list_elem(ell))%ptr_elem%id, ell=1,9)!N_max_element)
    ! 	write(*,'(" fjk (x10^10)     :  ", (1F14.1,1x),*(1F14.1,1x))') 1d10*max_fjk(0), (1d10*max_fjk(list_elem(ell)), ell=1,9)!N_max_element)

    return
  end subroutine show_electron_given_per_elem

  subroutine calc_ionisation_frac(elem, k, ne, fjk, dfjk, n0)
    !return j * Nj / Ntot and N0/Ntot

    ! ------------------------------------------------------------------------!
    ! fractional population f_j(ne,T)=N_j/N for element Elem
    ! and its partial derivative with ne. If Elem is an element with
    ! detailed model and if NLTE populations for it exist, there are used
    ! instead of LTE.
    !
    ! ::new:: n0 is population of the lowest ionisation stage (like HI)
    ! Should in principle be a list for the calculation of ionisation from negative ions.
    ! currently only H-
    !
    ! to test: If detailed_model is not applied only to active atoms
    ! it means that some atoms will contribute their LTE values, therefore
    ! they have to be updated. Is it okey or not ?
    ! OR it better to finder a consistant solution with LTE atoms (elements or not)
    ! ------------------------------------------------------------------------!

    real(kind=dp), intent(in) :: ne
    integer, intent(in) :: k
    type (Element), intent(in) :: Elem
    real(kind=dp), dimension(:), intent(inout) :: fjk, dfjk
    real(kind=dp), intent(out) :: n0
    real(kind=dp) :: Uk, Ukp1, sum1, sum2
    logical :: detailed_model
    integer :: j, i, Nstage, min_j

    Nstage = elem%Nstage
    fjk(:) = 0.0_dp
    dfjk(:) = 0.0_dp

    detailed_model = .false.

    if (associated(elem%model)) then
       !otherwise, the first call (wihtout knowing ne at all)
       !would not work since atom%n = 0.0 (or atom%nstar)
       !.true. before the maxval
       !detailed_model = (maxval(elem%model%n)>0.0).and.(elem%model%active)
       !-> two much time to test n > 0 for large model
       !or just nltepops ? for passive atoms ?
       detailed_model = (elem%model%NLTEpops .and. elem%model%active)
       !can add condition on active or not, but check bckgr opac and eval of lte pops after
    endif


    if (detailed_model) then
       if (.not.elem%model%active) call warning("Beware passive atoms used!")
       ! 		fjk = 0d0
       ! 		dfjk = 0d0
       !derivative is 0 = constant for subsequent iterations, for non-LTE loop.

       n0 = 0.0_dp
       min_j = minval(elem%model%stage)
       ! 		n0 = sum(elem%model%n(:,k),mask=elem%model%stage==min_j)

       do i = 1, elem%model%Nlevel


          fjk(elem%model%stage(i)+1) = fjk(elem%model%stage(i)+1)+( elem%model%stage(i)*elem%model%n(i,k) )
          if (elem%model%stage(i) == min_j) n0 = n0 + elem%model%n(i,k)


       end do

       fjk(1:Nstage) = fjk(1:Nstage)/nTotal_atom(k, elem%model)
       ! 		n0 = n0 / ntotal_atom(k,elem%model)
       n0 = 0.0_dp

       if (elem%model%id == "H") then
          min_f_HII = min(min_f_HII, fjk(2))
          max_f_HII = max(max_f_HII, fjk(2))
       endif

    else !no model use LTE
       !fjk(1) is N(1) / ntot !N(1) = sum_j=1 sum_i=1^N(j) n(i)
       fjk(1)=1.
       dfjk(1)=0.
       sum1 = 1.
       sum2 = 0.
       Uk = getPartitionFunctionk(elem,1,k)
       do j=2,Elem%Nstage
          Ukp1 = getPartitionFunctionk(elem,j,k)
          fjk(j) = Sahaeq(k,fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)

          ! 			if (ne>0) then
          ! 				dfjk(j) = -(j-1)*fjk(j)/ne
          ! 			else
          ! 				dfjk(j) = 0d0
          ! 			endif
          !-> fjk is zero if ne=0 in Sahaeq
          dfjk(j) = -(j-1)*fjk(j)/( ne + ne_min_limit)

          sum1 = sum1 + fjk(j)
          sum2 = sum2 + dfjk(j)

          Uk = Ukp1
       end do

       fjk(:)=fjk(:)/sum1
       dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1

       n0 = fjk(1)

       !0 for neutrals => j==1 (elements which do not have a detailed model are ordered by neutral to ionised
       !and fjk is Nj / Ntot
       do j=1,elem%Nstage
          fjk(j) = (j-1) * fjk(j)
          dfjk(j) = (j-1) * dfjk(j)
       enddo

       if (elem%id == "H") then
          min_f_HII = min(min_f_HII, fjk(2))
          max_f_HII = max(max_f_HII, fjk(2))
       endif

    endif

    return
  end subroutine calc_ionisation_frac

  !try without para ! for epsilon
  !get cells that have not changed with subsequent iteration in non-LTE loop
  subroutine solve_electron_density(ne_initial_solution, verbose, epsilon)
    ! ----------------------------------------------------------------------!
    ! Solve for electron density for a set of elements
    ! stored in Elements. Elements up to N_MAX_ELEMENT are
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
    !
    ! Also returns epsilon, the change in electron density from the intial guess
    ! to the converged values. THIS IS NOT the convergence threshold used inside
    ! the electron loop.
    ! Epsilon is used to check for the convergence of the electron density,
    ! between subsequent iterations of the non-LTE loop.
    !
    ! epsilon is shared, thats' okay ?
    ! ----------------------------------------------------------------------!
    character(len=20), intent(in) :: ne_initial_solution
    logical, intent(in) :: verbose
    real(kind=dp), intent(inout) :: epsilon !difference wrt the initial solution, not criterion of convergence (dne)
    real(kind=dp) :: delta, ne_old, akj, sum, Uk, Ukp1, dne, ne0
    real(kind=dp):: ne_oldM, UkM, PhiHmin, n0
    real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk
    real(kind=dp), dimension(-N_negative_ions:N_MAX_ELEMENT) :: max_fjk, min_fjk !Negative ions from -N_neg to 0 (H-), then from 1 to Nelem positive ions
    integer :: n, k, niter, j, ZM, id
    type (Element), pointer :: elem
    integer :: unconverged_cells(nb_proc), ik_max

    !!if (verbose) write(*,*) "Initial solution for electron loop:", ne_initial_solution

    unconverged_cells(:) = 0
    id = 1
    ik_max = 0

    epsilon = 0.0

    min_f_HII = 1d50
    max_f_HII = 0.0_dp
    max_fjk(:) = 0.0_dp

    !$omp parallel &
    !$omp default(none) &
    !$omp private(k,n,j,fjk,dfjk,ne_old,niter,delta,sum,PhiHmin,Uk,Ukp1,ne_oldM) &
    !$omp private(dne, akj, id, ne0, elem, n0) &
    !$omp shared(n_cells, Elements, ne_initial_solution,Hydrogen, ZM, unconverged_cells, Nelem) &
    !$omp shared(ne, T, icompute_atomRT, nHtot, epsilon, max_f_HII, min_f_HII, ik_max, max_fjk)
    !$omp do
    do k=1,n_cells
       !$ id = omp_get_thread_num() + 1

       ! 		if (icompute_atomRT(k) <= 0) cycle !transparent or dark
       ! 		if (icompute_atomRT(k)==2) then
       ! 			write(*,*) "id",id," skipping cell ", k, " for electron density! fixed values, ne = ", ne(k)
       ! 		endif
       if (icompute_atomRT(k) /= 1) cycle  !<=0 transparent or dark
       !==2 ne locked to model value
       !==1 compute with RT

       !! Initial solution for this cell!!

       !compute that value in any case of the starting solution
       !Metal
       Zm = 11
       elem => Elements(ZM)%ptr_elem
       Uk = getPartitionFunctionk(Elements(ZM)%ptr_elem, 1, k)
       Ukp1 = getPartitionFunctionk(Elements(ZM)%ptr_elem, 2, k)
       CALL ne_Metal(k, Uk, Ukp1, elements(ZM)%ptr_elem%ionpot(1),elements(ZM)%ptr_elem%Abund, ne_oldM)

       if (ne_initial_solution.eq."N_PROTON") then
          ne_old = Hydrogen%n(Hydrogen%Nlevel,k)!np(k)
       else if (ne_initial_solution.eq."NE_MODEL") then
          ne_old = ne(k)
       else !"H_IONISATION" or unkown

          !Initial solution ionisation of H
          elem => Elements(1)%ptr_elem
          Uk = getPartitionFunctionk(elem, 1, k)
          Ukp1 = getPartitionFunctionk(elem, 2, k)
          elem => NULL()

          if (Ukp1 /= 1.0_dp) then
             CALL Warning("Partition function of H+ should be 1")
             write(*,*) Uk, Ukp1
             stop
          end if
          CALL ne_Hionisation (k, Uk, Ukp1, ne_old)


          !if Abund << 1. and chiM << chiH then
          ! ne (H+M) = ne(H) + ne(M)
          ne_old = ne_old + ne_oldM
          !to avoid having very low values, between say tiny_dp and 1e-100
          !just set to 0 if < 0: crude ?
          ! 			if (ne_oldM < ne_min_limit) ne_oldM = 0d0
          if (ne_old < ne_min_limit) ne_old = 0.0_dp

       end if

       !Loop starts
       ne(k) = ne_old
       ne0 = ne_old
       niter=0
       do while (niter < N_MAX_ELECTRON_ITERATIONS)
          delta = ne_old/nHtot(k)
          sum = 0.0

          do n=1,Nelem
             if (n > N_MAX_ELEMENT) exit
             !do n=1, N_MAX_ELEMENT
             elem => Elements(n)%ptr_elem

             !times stage !!
             call calc_ionisation_frac(elem, k, ne_old, fjk, dfjk, n0)


             if (n.eq.1)  then ! H minus for H
                !2 = partition function of HI, should be replace by getPartitionFunctionk(elements(1)%ptr_elem, 1, icell)
                PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
                ! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)

                !in non-LTE Z * fjk is stored in fjk meaning 0 for fjk(1) !neutrals
                !How to deal with non-LTE populations with negative ions ?
                !For instance, fjk(1) = 0 for Hydrogen during the non-LTE (stage * sum(n) / Ntot)
                !but not at LTE.

                !H- must be multiply their by 1 or by sum(hydrogen%n(1:hydrogen%Nlevel-1))/nHtot
                !! replace fjk(1) by n0, should be the same at LTE
                !negative contribution to the electron density
                delta = delta + ne_old*PhiHmin*n0!
                sum = sum-(n0 + ne_old*dfjk(1))*PhiHmin ! (fjk(1) + ne_old*dfjk(1))

                max_fjk(0) = max(max_fjk(0),ne_old*PhiHmin*n0)

             end if

             !neutrals do not contribute
             !j-1 = 0 for neutrals
             !j-1 = 1 for singly ionised ions.

             !avoiding neutrals
             !fjk are Nj / Not with Nj the total population in stage j evaluated at LTE
             ! 				do j=2, elem%Nstage
             ! 					akj = elem%Abund*(j-1) !because j starts at 0 for neutrals, 1 for singly ionised etc
             ! 					!positive contribution to the electron density
             ! 					delta = delta -akj*fjk(j)
             ! 					sum = sum + akj*dfjk(j)
             ! 					max_fjk(n) = max(max_fjk(n), akj*fjk(j))
             ! 				end do
             akj = elem%Abund
             !new the term j-1 already included in fjk and dfjk
             do j=1, elem%Nstage !start at 1 but fjk and dfjk are 0 if stage 1 = neutral
                !this allows a better match with detailed models that can have ionised level for j=1 (like Ca II)
                !positive contribution to the electron density
                delta = delta -akj*fjk(j)
                sum = sum + akj*dfjk(j)
                max_fjk(n) = max(max_fjk(n), akj*fjk(j))
             end do

             elem => NULL()
          end do !loop over elem


          ne(k) = ne_old - nHtot(k)*delta / (1.0-nHtot(k)*sum)

          ! 			if (ne(k) == 0.0) then
          ! 				write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, sum
          ! 				call Warning("electron density is 0 setting to metal value")
          ! 				ne(k) = ne_oldM
          ! 			endif

          ! 			dne = abs((ne(k)-ne_old)/ne_old)
          dne = abs ( 1.0_dp - ne_old / ne(k) )
          ne_old = ne(k)

          ! 			test epsilon en para et pas para!!!
          !can I do that in parallel ? epsilon is shared
          ! 			compare to the initial solution (the previous solution in non-LTE)
          if (abs(1.0_dp - ne0 / ne(k)) > epsilon) then
             epsilon = abs(1.0_dp - ne0 / ne(k))
             ! 			if (abs((ne(k) - ne0) / ne0) > epsilon) then
             ! 				epsilon = abs((ne(k) - ne0) / ne0)
             ik_max = k
          endif

          if (is_nan_infinity(ne(k))>0) then
             write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, " sum=", sum
             call error("electron density is nan or inf!")
          else if (ne(k) < 0.0) then
             write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, sum
             call error("Negative electron density!")
          endif

          niter = niter + 1
          if (dne <= MAX_ELECTRON_ERROR) then
             if (ne(k) < ne_min_limit) then
                ne(k) = 0.0_dp !or ne(k) = ne_min_limit
                write(*,*) " (Solve ne) ne < ne_min_limit, setting cell transparent!", k
                icompute_atomRT(k) = 0
             endif
             exit
          else if (niter >= N_MAX_ELECTRON_ITERATIONS) then
             if (dne > 0.0) then !shows the warning only if dne is actually greater than 1
                CALL Warning("Electron density has not converged for this cell")
                write(*,*) "icell=",k,"maxIter=",N_MAX_ELECTRON_ITERATIONS,"dne=",dne, &
                     "max(err)=", MAX_ELECTRON_ERROR, "ne=",ne(k), "T=",T(k)," nH=",nHtot(k)
                !set articially ne to some value ?
                if (dne >= 1.0) then
                   ne(k) = ne_oldM !already tested if < ne_min_limit
                   write(*,*) " -> Setting ne to ionisation of ", elements(ZM)%ptr_elem%ID, ne_oldM
                endif
             end if
             unconverged_cells(id) = unconverged_cells(id) + 1 !but count each cell for with dne > MAX_ELECTRON_ERROR
          end if !convergence test

       end do !while loop

    end do !loop over spatial points
    !$omp end do
    !$omp end parallel

    if (verbose) then
       write(*,*) " ------------------------------------------------ "
       write(*,'("ne(min)="(1ES17.8E3)" m^-3 ;ne(max)="(1ES17.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
       write(*,'("   >>>  epsilon="(1ES17.8E3)" at cell "(1I5))') epsilon, ik_max
       write(*,*) " T = ", T(ik_max)," nH = ", nHtot(ik_max)
       write(*,*) " "
       write(*,'("Ionisation fraction of HII "(1ES17.8E3, 1ES17.8E3))') max_f_HII, min_f_HII
       ! 			call show_electron_given_per_elem(0, 0, max_fjk)
       write(*,*) " ------------------------------------------------ "
    endif

    !that's because I use the sum as a variable so the function doesn't exist.
    do k=2, nb_proc
       unconverged_cells(1) = unconverged_cells(1) + unconverged_cells(k)
    end do
    if (unconverged_cells(1) > 0) then
       write(*,*) "Found ", unconverged_cells(1)," unconverged cells (%):", 100*real(unconverged_Cells(1))/real(n_cells)
    endif

    return
  end subroutine solve_electron_density

  SUBROUTINE getfjk (Elem, ne, k, fjk, dfjk)
    ! ----------------------------------------------------------------------!
    ! fractional population f_j(ne,T)=N_j/N for element Elem
    ! and its partial derivative with ne. If Elem is an element with
    ! detailed model and if NLTE populations for it exist, their are used
    ! instead of LTE.
    ! ----------------------------------------------------------------------!

    real(kind=dp), intent(in) :: ne
    integer, intent(in) :: k
    type (Element), intent(in) :: Elem
    real(kind=dp), dimension(:), intent(inout) :: fjk, dfjk
    real(kind=dp) :: Uk, Ukp1, sum1, sum2
    logical :: has_nlte_pops
    integer :: j, i

    has_nlte_pops = .false. !reset for all atom

    if (associated(elem%model)) then
       !it also check that it is active ?
       !because a model can be associated but it doesn't mean it has nlte pops !
       if (elem%model%NLTEpops) then
          has_nlte_pops = .true.
       endif
    endif

    fjk(:) = 0d0
    dfjk(:) = 0d0

    !may be active without NLTEpops or passive with read NLTE pops
    if (has_nlte_pops) then
       fjk = 0d0
       dfjk = 0d0
       !For Nlevel, Nlevel stages
       !fj = Nj/Ntot
       !first, Nj = Sum_i nij
       do i = 1, elem%model%Nlevel
          fjk(elem%model%stage(i)+1) = fjk(elem%model%stage(i)+1)+(elem%model%stage(i))*elem%model%n(i,k)
       end do                                        !is there a + 1 here so that fjk(1)=1*atom%n(1)/Ntot
       !Divide by Ntotal and retrieve fj = Nj/N for all j
       fjk(:) = fjk(:)/(nHtot(k)*elem%model%Abund)
    else !not active or active but first iteration of the NLTEloop so that
       !NLTEpos has been set to .false., whateveeer, use LTE
       fjk(1)=1.
       dfjk(1)=0.
       sum1 = 1.
       sum2 = 0.
       Uk = getPartitionFunctionk(elem,1,k)
       do j=2,Elem%Nstage !-->vectorize ?
          Ukp1 = getPartitionFunctionk(elem,j,k)
          ! fjk(j) = Nj / Sum_j Nj
          ! Nj = Nj-1/(phi*ne) Saha Equation
          ! fj = Nj/N = Nj/N0 / N/N0
          ! -> first computes Nj/N0 using Nj-1/N0 relative to N0
          ! -> sum up = 1 + N1/N0 + Nj/N0
          ! -> divide Nj/N0 by this sum and retrive Nj/N for all j
          !  --> Nj/N0 / (1+Nj/N0) = Nj/(N0*(1+Nj/N0)) = Nj / (N0 + Nj) = fj
          !-> Check that again
          fjk(j) = Sahaeq(k,fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)
          !write(*,*) "j=",j," fjk=",fjk(j)
          !write(*,*) fjk(j-1), fjk(j), fjk(j)/(ne+tiny_dp)

          !why without this I have nan ???
          if (ne>0) then
             dfjk(j) = -(j-1)*fjk(j)/ne  !-> fjk(j) should be 0 if ne = 0
             !see Saheq.
             !dfjk(j) = -(j-1) * fjk(j) * min(huge_dp,1/ne)
          else
             dfjk(j) = 0d0
          endif

          sum1 = sum1 + fjk(j)
          sum2 = sum2 + dfjk(j)

          !write(*,*) "j=",j," dfjk=",dfjk(j), "fjk=", fjk(j), sum1, sum2


          Uk = Ukp1
       end do
       fjk(:)=fjk(:)/sum1
       dfjk(:)=(dfjk(:)-fjk(:)*sum2)/sum1

       if (any_nan_infinity_vector(dfjk)>0) then
          write(*,*) "j=",j," dfjk=",dfjk(j), "fjk=", fjk(j), ne, T(k), sum1, sum2
          stop
       endif
    end if

    RETURN
  END SUBROUTINE getfjk

  SUBROUTINE Solve_Electron_Density_old(ne_initial_solution)
    ! ----------------------------------------------------------------------!
    ! Solve for electron density for a set of elements
    ! stored in Elements. Elements up to N_MAX_ELEMENT are
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
    real(kind=dp) :: error, ne_old, akj, sum, Uk, dne, Ukp1
    real(kind=dp):: ne_oldM, UkM, PhiHmin
    !   real(kind=dp), dimension(n_cells) :: np
    real(kind=dp), dimension(:), allocatable :: fjk, dfjk
    integer :: Nmaxstage=0, n, k, niter, j, ZM, id, ninit, nfini
    integer :: unconverged_cells(nb_proc)

    if (.not. present(ne_initial_solution)) then
       initial="H_IONISATION"!use HIONISAtion
    else
       initial=ne_initial_solution
       !     if ((initial .neq. "N_PROTON") .neq. (initial /= "NE_MODEL")) then
       !      write(*,*) 'NE initial solution unkown, set to H_IONISATION'
       !      initial = "H_IONISATION"
       !     end if
    end if
    write(*,*) "Initial solution for electron loop:", initial

    unconverged_cells(:) = 0
    id = 1
    do n=1,Nelem
       Nmaxstage=max(Nmaxstage,Elements(n)%ptr_elem%Nstage)
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
    write(*,*) "Test Nelem",Nelem

    !$omp parallel &
    !$omp default(none) &
    !$omp private(k,n,j,fjk,dfjk,ne_old,niter,error,sum,PhiHmin,Uk,Ukp1,ne_oldM) &
    !$omp private(dne, akj, id, ninit, nfini) &
    !$omp shared(n_cells, Elements, initial,Hydrogen, ZM, unconverged_cells, Nelem) &
    !$omp shared(ne, T, icompute_atomRT, nHtot)
    !$omp do
    do k=1,n_cells
       !$ id = omp_get_thread_num() + 1
       if (icompute_atomRT(k) <= 0) cycle !transparent or dark


       !write(*,*) "The thread,", omp_get_thread_num() + 1," is doing the cell ", k
       if (initial.eq."N_PROTON") then
          ne_old = Hydrogen%n(Hydrogen%Nlevel,k)!np(k)
       else if (initial.eq."NE_MODEL") then
          ne_old = ne(k)
       else !"H_IONISATION" or unkown

          !Initial solution ionisation of H
          Uk = getPartitionFunctionk(Elements(1)%ptr_elem, 1, k)
          Ukp1 = getPartitionFunctionk(Elements(1)%ptr_elem, 2, k)

          if (Ukp1 /= 1d0) then
             CALL Warning("Partition function of H+ should be 1")
             write(*,*) Uk, Ukp1
             stop
          end if
          CALL ne_Hionisation (k, Uk, Ukp1, ne_old)


          if (T(k) >= 20d3) then
             ZM = 2 !He
             !else if (T(k) <= 1000d0) then
          else
             ZM = 11 ! Na !trade off between higher A and lower ionpot
             !else
             !  ZM = 26 ! Fe
          end if
          !add Metal
          Uk = getPartitionFunctionk(Elements(ZM)%ptr_elem, 1, k)
          Ukp1 = getPartitionFunctionk(Elements(ZM)%ptr_elem, 2, k)
          CALL ne_Metal(k, Uk, Ukp1, elements(ZM)%ptr_elem%ionpot(1),elements(ZM)%ptr_elem%Abund, ne_oldM)
          !write(*,*) "neMetal=",ne_oldM
          !if Abund << 1. and chiM << chiH then
          ! ne (H+M) = ne(H) + ne(M)
          ne_old = ne_old + ne_oldM

!!!!Check here !!!
          !to avoid having very low values, between say tiny_dp and 1e-100
          !just set to 0 if < 0: crude ?
          if (ne_oldM < ne_min_limit) ne_oldM = 0d0
          if (ne_old < ne_min_limit) ne_old = 0d0
!!!!!!!!!!!!!!!!!
       end if

       ne(k) = ne_old
       niter=0
       do while (niter < N_MAX_ELECTRON_ITERATIONS)
          error = ne_old/nHtot(k)
          sum = 0.

          !ninit = (1. * (id-1)) / nb_proc * Nelem + 1
          !nfini = (1. * id) / nb_proc * Nelem
          do n=1,Nelem
             if (n > N_MAX_ELEMENT) exit

             CALL getfjk(Elements(n)%ptr_elem,ne_old,k,fjk,dfjk)!

             if (n.eq.1)  then ! H minus for H
                !2 = partition function of HI, should be replace by getPartitionFunctionk(elements(1)%ptr_elem, 1, icell)
                PhiHmin = phi_jl(k, 1d0, 2d0, E_ION_HMIN)
                ! = 1/4 * (h^2/(2PI m_e kT))^3/2 exp(Ediss/kT)
                error = error + ne_old*fjk(1)*PhiHmin
                sum = sum-(fjk(1)+ne_old*dfjk(1))*PhiHmin
                !write(*,*) "phiHmin=",PhiHmin,error, sum
             end if
             !neutrals do not contribute right
             do j=2,elements(n)%ptr_elem%Nstage
                akj = elements(n)%ptr_elem%Abund*(j-1) !because j starts at 0 for neutrals, 1 for singly ionised etc
                error = error -akj*fjk(j)
                sum = sum + akj*dfjk(j)
                !write(*,*) n-1, j-1, akj, error, sum, dfjk(j)
             end do
          end do !loop over elem

          ne(k) = ne_old - nHtot(k)*error /&
               (1.-nHtot(k)*sum)
          dne = dabs((ne(k)-ne_old)/(ne_old+tiny_dp))
          ne_old = ne(k)

          if (is_nan_infinity(ne(k))>0) then
             write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), &
                  "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, sum
             stop
          endif
          niter = niter + 1
          if (dne <= MAX_ELECTRON_ERROR) then
             !write(*,*) niter, "icell=",k," T=",T(k)," ne=",ne(k), " dne=", dne
             if (ne(k) < ne_min_limit) ne(k) = 0d0 !tiny_dp or < tiny_dp
             exit
          else if (niter >= N_MAX_ELECTRON_ITERATIONS) then
             if (dne >= 1) then !shows the warning only if dne is actually greater than 1
                CALL Warning("Electron density has not converged for this cell")
                write(*,*) "icell=",k,"maxIter=",N_MAX_ELECTRON_ITERATIONS,"dne=",dne,"max(err)=", MAX_ELECTRON_ERROR, &
                     "ne=",ne(k), "T=",T(k)," nH=",nHtot(k)
                !set articially ne to some value ?
                ne(k) = ne_oldM !already tested if < ne_min_limit
                write(*,*) "Setting ne to ionisation of ", elements(ZM)%ptr_elem%ID, ne_oldM
                !       else
                !           write(*,*) "icell=",k,"maxIter=",N_MAX_ELECTRON_ITERATIONS,"dne=",dne,"max(err)=", MAX_ELECTRON_ERROR, &
                !           "ne=",ne(k), "T=",T(k)," nH=",nHtot(k)
             end if
             unconverged_cells(id) = unconverged_cells(id) + 1 !but count each cell for with dne > MAX_ELECTRON_ERROR
          end if !convergence test

       end do !while loop
    end do !loop over spatial points
    !$omp end do
    !$omp end parallel

    deallocate(fjk, dfjk)

    write(*,*) "Maximum/minimum Electron density in the model (m^-3):"
    write(*,*) MAXVAL(ne),MINVAL(ne,mask=icompute_atomRT>0)

    !that's because I use the sum as a variable so the function doesn't exist.
    do k=2, nb_proc
       unconverged_cells(1) = unconverged_cells(1) + unconverged_cells(k)
    end do
    if (unconverged_cells(1) > 0) write(*,*) "Found ", unconverged_cells(1)," unconverged cells (%):", &
         100*real(unconverged_Cells(1))/real(n_cells)

    RETURN
  END SUBROUTINE Solve_Electron_Density_old


END MODULE solvene
