! ------------------------------------------------------------------- !
! ------------------------------------------------------------------- !
! Module that solves for electron density for given
! temperature grid, Hydrogen total populations
! elements and their abundance and their LTE or NLTE populations.
!
!
! If ne_intial is 1 then first guess is computed using
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

module elecdensity

   use atom_type, only : atoms, AtomType, hydrogen
   use elements_type
   use constantes
   use messages, only : Error, Warning
   use input, only : nb_proc
   use mcfost_env, only : dp
   use parametres, only : n_cells
   use utils, only : progress_bar, is_nan_infinity, interp_dp
   use grid, only : T, ne, nHtot, icompute_atomRT
   use fits_utils
   !$ use omp_lib

   implicit none

   integer, parameter :: N_negative_ions = 0 !currently only H- at index 0
   real(kind=dp), parameter :: MAX_ELECTRON_ERROR=1d-6
   integer, parameter :: N_MAX_ELECTRON_ITERATIONS=50
   integer, parameter :: N_MAX_ELEMENT=26!stops at Iron.
   real(kind=dp) :: min_f_HII, max_f_HII
   character(len=10) :: ne_filename = "ne.fits.gz"
   real(kind=dp), parameter :: ne_min_limit = 1d-10 ! ne_min_limit * nHtot = electron per cubic meters.
   !if ne is below ne_min_limit, we set it to 1 electron per cc (ne_one) !
   real(kind=dp), parameter :: ne_one = 1.0_dp !one electron per cubic meters.

   !Introducing a lock for ne iterations with icompute_atomRT == 2
   !the transfer (and opacity) is solved for icompute_atomRT > 0
   !but if icompute_atomRT==2 electron density has the init value.
   !then electron iterations skipped for icompute_atomRT /= 1

   contains


   function ne_Hionisation0 (temp, N, U0, U1)
   ! ----------------------------------------------------------------------!
   ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
   ! of H.
   ! Njl = Nj1l * ne * phi_jl
   ! ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
   ! eq. 17.77 and 17.78
   ! ne(H) = (sqrt(N*phiH + 1)-1)/phiH
   ! ----------------------------------------------------------------------!

   real(kind=dp), intent(in) :: temp, N
   real(kind=dp), intent(in) :: U0, U1
   real(kind=dp) :: phiH
   real(kind=dp) :: ne_Hionisation0
   phiH = phi_jl(temp, U0, U1, elems(1)%ionpot(1))

   !if no free e-, Nt = NH + NH+ with NH+=ne
   if (phiH>0) then
      ne_Hionisation0 = (sqrt(N*phiH*4. + 1)-1)/(2.*phiH) !without H minus
   else
      ne_Hionisation0 = 0.0
   endif

   return
   end function ne_Hionisation0

   function ne_Hionisation (temp, N, U0, U1)
    ! ----------------------------------------------------------------------!
    ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
    ! of H.
    ! Njl = Nj1l * ne * phi_jl
    ! ne(H) = NH-1 / (NH * phi_-1l) chi = ionpot0
    ! eq. 17.77 and 17.78
    ! ne(H) = (sqrt(N*phiH + 1)-1)/phiH
    ! ----------------------------------------------------------------------!

    real(kind=dp), intent(in) :: temp, N
    real(kind=dp), intent(in) :: U0, U1
    real(kind=dp) :: phiH
    real(kind=dp) :: ne_Hionisation

    phiH = phi_jl(temp, U0, U1, Elems(1)%ionpot(1))

    !if no free e-, Nt = NH + NH+ with NH+=ne
    !ne = (sqrt(nHtot(k)*phiH*4. + 1)-1)/(2.*phiH) !without H minus
    !if free e-, Nt=Nhot + NH+ + ne
    if (phiH>0) then
      ne_Hionisation = (sqrt(N*phiH + 1)-1)/(phiH)
    else
      ne_Hionisation = 0.0
    endif

    return
   end function ne_Hionisation

   function ne_Metal(temp, N, U0, U1, chi, A)
    ! ----------------------------------------------------------------------!
    ! Application of eq. 4.35 of Hubeny & Mihalas to ionisation
    ! of a single metal.
    ! ----------------------------------------------------------------------!

    real(kind=dp), intent(in) :: temp, N
    real(kind=dp), intent(in) :: U0, U1, chi, A
    real(kind=dp) :: phiM, alphaM
    real(kind=dp) :: ne_Metal

    phiM = phi_jl(temp, U0, U1, chi)
    alphaM = A ! relative to H, for instance 1.-6 etc

    if (phiM > 0) then
      ne_Metal = (sqrt(alphaM*N*phiM +0.25*(1+alphaM)**2) - 0.5*(1+alphaM) ) / phiM
    else
      ne_Metal = 0.0
    endif

    RETURN
   end function ne_Metal

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

          write(*,'("  >> "(1A)" fjk="(1ES17.8E3))') elems(list_elem(ell))%id, max_fjk(list_elem(ell))

       endif
    enddo


    !not working properly in term of display
    ! 	write(*,'(" Element          :  ", (A,1x),*(A,10x))') "H-", (elems(list_elem(ell))%id, ell=1,9)!N_max_element)
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
    type(AtomType), pointer :: atom
    real(kind=dp) :: ntotal

    Nstage = elem%Nstage
    fjk(:) = 0.0_dp
    dfjk(:) = 0.0_dp

    detailed_model = .false.

    if (elem%nm>0) then
       !otherwise, the first call (wihtout knowing ne at all)
       !would not work since atom%n = 0.0 (or atom%nstar)
       !.true. before the maxval
       !detailed_model = (maxval(elem%model%n)>0.0).and.(elem%model%active)
       !-> two much time to test n > 0 for large model
       !or just nltepops ? for passive atoms ?
       detailed_model = (Atoms(elem%nm)%p%NLTEpops .and. Atoms(elem%nm)%p%active)
      !can add condition on active or not, but check bckgr opac and eval of lte pops after
    endif


    !model and active (non-LTE)
    if (detailed_model) then
       atom => Atoms(elem%nm)%p
       ntotal = atom%Abund * nHtot(k)
       ! 		fjk = 0d0
       ! 		dfjk = 0d0
       !derivative is 0 = constant for subsequent iterations, for non-LTE loop.

       n0 = 0.0_dp
       min_j = minval(atom%stage)
       !n0 = sum(elem%model%n(:,k),mask=elem%model%stage==min_j)

       do i = 1,atom%Nlevel


          fjk(atom%stage(i)+1) = fjk(atom%stage(i)+1)+( atom%stage(i)*atom%n(i,k) )
          if (atom%stage(i) == min_j) n0 = n0 + atom%n(i,k)


       end do

       fjk(1:Nstage) = fjk(1:Nstage)/ntotal
       ! 		n0 = n0 / ntotal
       n0 = 0.0_dp

       if (atom%id == "H") then
          min_f_HII = min(min_f_HII, fjk(2))
          max_f_HII = max(max_f_HII, fjk(2))
       endif

    else !no model use LTE
       !fjk(1) is N(1) / ntot !N(1) = sum_j=1 sum_i=1^N(j) n(i)
       fjk(1)=1.
       dfjk(1)=0.
       sum1 = 1.
       sum2 = 0.
       Uk = get_pf(elem,1,T(k))
       do j=2,Elem%Nstage
          Ukp1 = get_pf(elem,j,T(k))
          fjk(j) = Sahaeq(T(k),fjk(j-1),Ukp1,Uk,elem%ionpot(j-1),ne)

          ! 			if (ne>0) then
          ! 				dfjk(j) = -(j-1)*fjk(j)/ne
          ! 			else
          ! 				dfjk(j) = 0d0
          ! 			endif
          !-> fjk is zero if ne=0 in Sahaeq
          if (ne > 0.0) then
            dfjk(j) = -(j-1)*fjk(j)/ne
          else
            dfjk(j) = 0.0
          endif

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

  subroutine solve_ne_loc(k,ne_init)
   !solve for electronic density locally (one point).
   integer, intent(in) :: k
   real(kind=dp), intent(in) :: ne_init
   real(kind=dp) :: delta, ne_old, akj, sum, Uk, Ukp1, dne
   real(kind=dp):: PhiHmin, n0
   real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk
   !!real(kind=dp), dimension(-N_negative_ions:N_MAX_ELEMENT) :: max_fjk, min_fjk !Negative ions from -N_neg to 0 (H-), then from 1 to Nelem positive ions
   integer :: n, niter, j


   !ne_old is updated inside
   !difference between ne(k) and ne_init gives the different to the initial solution
   !not the criterion of convergence
   ne_old = ne_init
   
   !Loop starts
   ne(k) = ne_old
   niter=0
   do while (niter < N_MAX_ELECTRON_ITERATIONS)
      delta = ne_old/nHtot(k)
      sum = 0.0

      do n=1,Nelem
         if (n > N_MAX_ELEMENT) exit
            !do n=1, N_MAX_ELEMENT

         !times stage !!
         call calc_ionisation_frac(elems(n), k, ne_old, fjk, dfjk, n0)


         if (n.eq.1)  then ! H minus for H
               !2 = partition function of HI, should be replace by getPartitionFunctionk(elements(1)%ptr_elem, 1, icell)
            PhiHmin = phi_jl(T(k), 1d0, 2d0, E_ION_HMIN)
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

            !!max_fjk(0) = max(max_fjk(0),ne_old*PhiHmin*n0)

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
         akj = elems(n)%Abund
            !new the term j-1 already included in fjk and dfjk
         do j=1, elems(n)%Nstage !start at 1 but fjk and dfjk are 0 if stage 1 = neutral
               !this allows a better match with detailed models that can have ionised level for j=1 (like Ca II)
               !positive contribution to the electron density
            delta = delta -akj*fjk(j)
            sum = sum + akj*dfjk(j)
            !!max_fjk(n) = max(max_fjk(n), akj*fjk(j))
         end do

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



      if (is_nan_infinity(ne(k))>0) then
         write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, " sum=", sum
         call error("electron density is nan or inf!")
      else if (ne(k) <= 0.0) then
         write(*,*) niter, "icell=",k," T=",T(k)," nH=",nHtot(k), " ne0=", ne_init
         write(*,*) "dne = ",dne, " ne=",ne(k), " nedag = ", ne_old, "sum=", sum
         call error("Negative electron density!")
      endif

      niter = niter + 1
      if (dne <= MAX_ELECTRON_ERROR) then
         if (ne(k) < nHtot(k)*ne_min_limit) then
               write(*,*) " (Solve ne) ne < ne_min_limit, setting cell transparent!", k
               write(*,*) "T=", T(k), ' nHtot=', nHtot(k), " ne=", ne(k)
               write(*,*) " -> setting ne to ", ne_one, " m^-3"
               ne(k) = ne_one
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
                  ne(k) = ne_init
                  write(*,*) " -> Setting ne to initial variable "
               endif
         end if
      end if !convergence test

   end do !while loop

   return
  end subroutine solve_ne_loc

  
  subroutine solve_ne(initial, verbose, epsilon)
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
    integer, intent(in) :: initial
    logical, intent(in) :: verbose
    !difference wrt the initial solution, not criterion of convergence (dne)
    real(kind=dp), intent(inout) :: epsilon
    real(kind=dp):: ne_oldM, UkM, PhiHmin, ne0, Uk, Ukp1, eps_id(nb_proc)
    real(kind=dp), dimension(max_ionisation_stage) :: fjk, dfjk
    real(kind=dp), dimension(-N_negative_ions:N_MAX_ELEMENT) :: max_fjk, min_fjk !Negative ions from -N_neg to 0 (H-), then from 1 to Nelem positive ions
    integer :: k, ZM, id, ibar, n_cells_done, n_cells_skipped
    integer :: unconverged_cells(nb_proc), ik_max, ik_max_id(nb_proc)

    ibar = 0
    n_cells_done = 0
    n_cells_skipped = 0!size(pack(icompute_atomRT,mask=icompute_atomRT /= 1))
    
    unconverged_cells(:) = 0
    eps_id(:) = 0.0
    ik_max_id = 0.0
    id = 1

    !difference with initial solution
    epsilon = 0.0
    ik_max = 0

    min_f_HII = 1d50
    max_f_HII = 0.0_dp
    max_fjk(:) = 0.0_dp


    call progress_bar(0)
    !$omp parallel &
    !$omp default(none) &
    !$omp private(k, ne_oldM, id, ne0, Uk, Ukp1) &
    !$omp shared(n_cells, initial, Hydrogen, ZM, unconverged_cells, Elems)&
    !$omp shared(ik_max_id, eps_id, ne, T, icompute_atomRT, nHtot) &
    !$omp shared(max_f_HII, min_f_HII, ik_max, max_fjk, ibar, n_cells_done, n_cells_skipped)
    !$omp do schedule(dynamic)
    do k=1,n_cells
       !$ id = omp_get_thread_num() + 1

       if (icompute_atomRT(k) /= 1) cycle
       !<=0 transparent or dark
       !==2 ne locked to model value
       !==1 compute with RT

       !compute that value in any case of the starting solution
       !Metal
       Zm = 11
       Uk = get_pf(Elems(ZM), 1, T(k))
       Ukp1 = get_pf(Elems(ZM), 2, T(k))
       ne_oldM = ne_Metal(T(k),nHtot(k), Uk, Ukp1, elems(ZM)%ionpot(1),elems(ZM)%Abund)

       if (initial == 2) then
          ne0 = Hydrogen%n(Hydrogen%Nlevel,k)!np(k)
       else if (initial == 1) then
          ne0 = ne(k)
       else !"H_IONISATION" or unkown

          !Initial solution ionisation of H
          Uk = get_pf(elems(1), 1, T(k))
          Ukp1 = get_pf(elems(1), 2, T(k))

          if (Ukp1 /= 1.0_dp) then
             CALL Warning("Partition function of H+ should be 1")
             write(*,*) Uk, Ukp1
             stop
          end if
          ne0 = ne_Hionisation (T(k),nHtot(k), Uk, Ukp1)


          !if Abund << 1. and chiM << chiH then
          ! ne (H+M) = ne(H) + ne(M)
          ne0 = ne0 + ne_oldM
          if (ne0 < nHtot(k)*ne_min_limit) ne0 = ne_one

       end if

       if (t(k) > 1d6) then
       !fully ionised
         ne(k) = 1.2 * nHtot(k)
       else
       !Loop starts
         call solve_ne_loc(k, ne0)
       endif
      
 
      if (abs(1.0_dp - ne0 / ne(k)) > eps_id(id)) then
         eps_id(id) = abs(1.0_dp - ne0 / ne(k))
         ik_max_id(id) = k
      endif

       
       ! Progress bar
       !$omp atomic
       n_cells_done = n_cells_done + 1
       if (real(n_cells_done) > 0.02*ibar*(n_cells-n_cells_skipped)) then
       	call progress_bar(ibar)
       	!$omp atomic
       	ibar = ibar+1
       endif

    end do !loop over spatial points
    !$omp end do
    !$omp end parallel
    call progress_bar(50)
    epsilon = maxval(eps_id)
    ik_max = ik_max_id(locate(eps_id, epsilon))

    if (verbose) then
       write(*,*) " ---------------------------------------------------- "
       write(*,'("ne(min)="(1ES16.8E3)" m^-3 ;ne(max)="(1ES16.8E3)" m^-3")') minval(ne,mask=icompute_atomRT>0), maxval(ne)
       write(*,'("   >>>  Diff to previous solution="(1ES13.5E3)" at cell "(1I7))') epsilon, ik_max
       write(*,*) " T = ", T(ik_max)," nH = ", nHtot(ik_max)
       write(*,*) " "
       write(*,'("Ionisation fraction of HII "(1ES13.5E3, 1ES13.5E3))') max_f_HII, min_f_HII
      !  write(*,*) nHtot(locate(ne/(1d-50+nHtot),maxval(ne/(1d-50+nHtot))))
       write(*,'("nH/ne "(1ES13.5E3, 1ES13.5E3))') maxval(nHtot/ne,mask=ne>0), minval(nHtot/ne,mask=ne>0)
       ! 			call show_electron_given_per_elem(0, 0, max_fjk)
       write(*,*) " ---------------------------------------------------- "
    endif

    !that's because I use the sum as a variable so the function doesn't exist.
    do k=2, nb_proc
       unconverged_cells(1) = unconverged_cells(1) + unconverged_cells(k)
    end do
    if (unconverged_cells(1) > 0) then
       write(*,*) "Found ", unconverged_cells(1)," unconverged cells (%):", 100*real(unconverged_Cells(1))/real(n_cells)
    endif

    return
  end subroutine solve_ne


  



  subroutine write_Electron
   ! ------------------------------------ !
   ! write Electron density if computed
   ! by the code.
   ! The electron density, ne, is in m^-3
   ! ------------------------------------ !
   integer :: unit, EOF = 0, blocksize, naxes(4), naxis,group, bitpix, fpixel
   logical :: extend, simple
   integer :: nelements, nfirst,sys_status
   character(len=512) :: cmd


   !check if data file already exist, can be the case if initial solutions is OLD_POPULATIONS
   cmd = "ls "//trim(ne_filename)
   call appel_syst(cmd, sys_status)
   if (sys_status == 0) then !means the file exist

      cmd = "mv "//trim(ne_filename)//" "//"ne_old.fits.gz"
      call appel_syst(cmd, sys_status)
      if (sys_status /= 0) then
         call error("Error in copying old ne!")
      endif

   endif

   !get unique unit number
   call ftgiou(unit,EOF)
   blocksize=1
   simple = .true. !Standard fits
   group = 1 !??
   fpixel = 1
   extend = .true.
   bitpix = -64

   if (lVoronoi) then
      naxis = 1
      naxes(1) = n_cells
      nelements = naxes(1)
   else
      if (l3D) then
         naxis = 3
         naxes(1) = n_rad
         naxes(2) = 2*nz
         naxes(3) = n_az
         nelements = naxes(1) * naxes(2) * naxes(3)
      else
         naxis = 2
         naxes(1) = n_rad
         naxes(2) = nz
         nelements = naxes(1) * naxes(2)
      end if
   end if


   call ftinit(unit,trim(ne_filename),blocksize,EOF)
   !  Initialize parameters about the FITS image
   call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,EOF)
   ! Additional optional keywords
   call ftpkys(unit, "UNIT", "m^-3", ' ', EOF)
   !write data
   call ftpprd(unit,group,fpixel,nelements,ne,EOF)

   call ftclos(unit, EOF)
   call ftfiou(unit, EOF)

   if (EOF > 0) call print_error(EOF)

   return
 end subroutine write_Electron

 subroutine read_electron(lelectron_read)
   ! ------------------------------------ !
   ! read Electron density from file
   ! ------------------------------------ !
   logical, intent(out) :: lelectron_read
   integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
   logical :: extend, simple, anynull
   integer :: nelements, naxis2(4), sys_status, naxis_found, hdutype
   character(len=512) :: cmd, some_comments
   real :: nullval

   !check if data file already exist, otherwise lead and return .false.
   cmd = "ls "//trim(ne_filename)
   call appel_syst(cmd, sys_status)
   if (sys_status == 0) then !means the file exist
      write(*,*) " Reading old ne.fits.gz file"
      lelectron_read = .true.
   else
      write(*,*) " found no ne.fits.gz"
      lelectron_read = .false.
      return
   endif

   status = 0
   call ftgiou(unit,status)

   call ftopen(unit, trim(ne_filename), 0, blocksize, status)
   if (status > 0) then
      write(*,*) "Cannot open electron file! ", ne_filename
      call print_error(status)
   endif


   simple = .true. !Standard fits
   group = 1
   fpixel = 1
   extend = .true.
   bitpix = -64
   nullval = -999

   call ftmahd(unit,1,hdutype,status)
   if (status > 0) then
      write(*,*) "Cannot read ne! "
      call print_error(status)
      stop
   endif

   if (lVoronoi) then
      naxis = 1

      ! number of axes
      call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
      if (status > 0) then
         write(*,*) "error reading number of axis (naxis)"
         call print_error(status)
         stop
      endif

      call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading ncells from file (naxis1)"
         call print_error(status)
         stop
      endif
      nelements = naxis_found

      if (nelements /= n_cells) then
         write(*,*) " found ",nelements," but need ",n_cells,": read_electron does not do interpolation yet!"
         call error(" Model read does not match simulation box")
      endif

      call ftgpvd(unit,group,1,nelements,nullval,ne,anynull,status)
      if (status > 0) then
         write(*,*) "error reading ne density"
         call print_error(status)
         stop
      endif

   else

      if (l3D) then
         naxis = 3
      else
         naxis = 2
      end if

      !Number of axis
      call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
      if (status > 0) then
         write(*,*) "error reading number of axis (naxis)"
         call print_error(status)
         stop
      endif

      !R
      call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading nrad from file (naxis1)"
         call print_error(status)
         stop
      endif
      nelements = naxis_found

      !Z
      call ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading nz from file (naxis2)"
         call print_error(status)
         stop
      endif
      nelements = nelements * naxis_found


      !phi axis ?
      if (l3D) then
         call ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
         if (status > 0) then
            write(*,*) "error reading naz from file (naxis3)"
            call print_error(status)
            stop
         endif
         nelements = nelements * naxis_found
      endif

      if (nelements /= n_cells) then
         write(*,*) " read_electron does not do interpolation yet!"
         call Error (" Model read does not match simulation box")
      endif

      call FTG2Dd(unit,group,-999,shape(ne),n_cells,1,ne,anynull,status)
      !call ftgpvd(unit,group,1,n_cells,-999,ne,anynull,status)

      if (status > 0) then
         write(*,*) "error reading ne density"
         call print_error(status)
         stop
      endif

   endif !lvoronoi

   call ftclos(unit, status) !close
   if (status > 0) then
      write(*,*) "error cannot close file in ", trim(ne_filename)
      call print_error(status)
      stop
   endif

   call ftfiou(unit, status) !free
   if (status > 0) then
      write(*,*) "error cannot free file unit!"
      call print_error(status)
      ! 		stop
   endif

   !call FTG2Dd(unit,1,-999,shape(ne),naxes(1),naxes(2),ne,anynull,EOF)

   write(*,'("  -- min(ne)="(1ES20.7E3)" m^-3; max(ne)="(1ES20.7E3)" m^-3")') , minval(ne,mask=(ne>0)), maxval(ne)

   return
 end subroutine read_electron

END MODULE elecdensity
