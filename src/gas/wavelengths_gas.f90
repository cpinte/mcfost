module wavelengths_gas


   use atom_type
   use grid, only : v_char, B_char, vturb, T
   use constantes
   use utils, only : locate, span, spanl, spanl_dp, span_dp
   use parametres
   use utils, only : span, spanl, spanl_dp, span_dp
   use sort, only : index_bubble_sort
   use messages
 
   implicit none
 
   !Number of points for each transition
   !continuum wavelength double for level dissolution !
   integer, parameter :: Nlambda_cont = 101!61!41! 141 !continuum, linear
   integer, parameter :: Nlambda_cont_log = 31!61!31!ontinuum log scaled
   integer, parameter :: Nlambda_line_w = 24, Nlambda_line_c_log = 51 !12 et 31
   integer, parameter :: Nlambda_line_c = 51!line linear1
   real, parameter    :: hvel_nlte = 6.0!for line in km/s, 1-3 for static models
   real, parameter	 :: delta_lambda_cont = 5.0 !nm
   real               :: hv = hvel_nlte !can change due to image grid
   integer            :: Nlambda_max_line, Nlambda_max_cont, Nlambda_max_trans !number max of lambda for all lines

   !real(kind=dp) :: group
 
   contains

   subroutine make_sub_wavelength_grid_cont_log_nu(cont, lambdamin, lambdamax)
   ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   ! ----------------------------------------------------------------- !
      type (AtomicContinuum), intent(inout) :: cont
      real(kind=dp), intent(in) :: lambdamin, lambdamax
      real(kind=dp) :: resol
      integer :: la, N1, N2
      real(kind=dp) :: l0, l1
      real(kind=dp) :: nu1, nu0, nu2, dnu

      !write(*,*) "Atom for which the continuum belongs to:", cont%atom%ID

      l1 = lambdamax
      l0 = lambdamin

      if (lambdamax > cont%lambda0) then
         N2 = Nlambda_cont_log + 1
         N1 = Nlambda_cont_log
      else
         N2 = 0
         N1 = Nlambda_cont_log
      endif
      cont%Nlambda = N1 + N2
      ! if (lambdamax > cont%lambda0) then
      !    N2 = cont%Nlambda + 1
      !    N1 = cont%Nlambda
      !    cont%Nlambda = N1 + N2
      ! else
      !    N2 = 0
      !    N1 = cont%Nlambda
      ! endif
      allocate(cont%lambda(cont%Nlambda))


      nu0 = (M_TO_NM * C_LIGHT / cont%lambda0)/1e15
      nu1 = (M_TO_NM * C_LIGHT / lambdamin)/1e15
      nu2 = (M_TO_NM * C_LIGHT / lambdamax)/1e15


      cont%lambda(N1:1:-1) = 1e-15 * c_light / spanl_dp(nu0,nu1,N1,1) * m_to_nm

      if (N2 > 0) then
         dnu = abs(cont%lambda(N1)-cont%lambda(N1-1))
         cont%lambda(N2+N1:N1+1:-1) = 1e-15 * c_light / spanl_dp(nu1+dnu,nu2,N2,1) * m_to_nm
      endif


      do la=1,cont%Nlambda
         if (cont%lambda(la) < 0) then
            write(*,*) "Error, lambda negative"
            write(*,*) "cont log"

            stop
         endif
      end do

      return
   end subroutine make_sub_wavelength_grid_cont_log_nu
 
 
   subroutine make_sub_wavelength_grid_cont_lin(cont, lambdamin, lambdamax)
   ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   ! ----------------------------------------------------------------- !
     type (AtomicContinuum), intent(inout) :: cont
     real(kind=dp), intent(in) :: lambdamin, lambdamax
     real(kind=dp) :: resol
     integer :: la, N1, N2, Nmid
     real(kind=dp) :: l0, l1, dl
 
     l1 = lambdamax
     l0 = lambdamin
 
     if (lambdamax > cont%lambda0) then
        N2 = Nlambda_cont + 1
        N1 = Nlambda_cont
        Nmid = N1 + N2
     else
        N2 = 0
        N1 = Nlambda_cont
        Nmid = N1
     endif
     cont%Nlambda = N1 + N2
     allocate(cont%lambda(cont%Nlambda))
 
     cont%lambda(1:N1) = span_dp(l0, cont%lambda0,N1,-1)!from lam0 to 1
     !still result in lamin to lamax
 
     ! N1+1 points are the same
     if (N2>0) then
        dl = abs(cont%lambda(N1) - cont%lambda(N1-1)) !2-1
        cont%lambda(N1+1:N2+N1) = span_dp(cont%lambda0+dl,l1,N2,1)
     endif
 
     do la=1,cont%Nlambda
        if (cont%lambda(la) < 0) then
           write(*,*) "Error, lambda negative"
           write(*,*) "cont lin"
           stop
        endif
     end do

     return
   end subroutine make_sub_wavelength_grid_cont_lin
 

   subroutine make_sub_wavelength_grid_cont_linlog(cont, lambdamin, lambdamax)
   ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   ! linear from 0 to lambda0 then log from lambda0 to lambdamax
   ! ----------------------------------------------------------------- !
     type (AtomicContinuum), intent(inout) :: cont
     real(kind=dp), intent(in) :: lambdamin, lambdamax
     real(kind=dp) :: resol
     integer :: la, N1, N2
     real(kind=dp) :: l0, l1, dl
 
 
     l1 = lambdamax
     l0 = lambdamin
 
     if (lambdamax > cont%lambda0) then
        N2 = Nlambda_cont_log + 1
        N1 =  Nlambda_cont
     else
        N2 = 0
        N1 =  Nlambda_cont
     endif
     cont%Nlambda = N1 + N2
     allocate(cont%lambda(cont%Nlambda))
 
 
     cont%lambda(1:N1) = span_dp(l0, cont%lambda0,N1,-1)!from lam0 to 1
 
 
     if (N2 > 0) then
        dl = abs(cont%lambda(N1)-cont%lambda(N1-1))
        cont%lambda(1+N1:N2+N1) = spanl_dp(cont%lambda0+dl,l1,N2,1)
     endif
 
 
     do la=1,cont%Nlambda
        if (cont%lambda(la) < 0) then
           write(*,*) "Error, lambda negative"
           write(*,*) "cont log"
 
           stop
        endif
     end do
 
     return
   end subroutine make_sub_wavelength_grid_cont_linlog
 
   function line_u_grid(k, line, N)
      !return for line line and cell k, the paramer:
      ! u = (lambda - lambda0) / lambda0 * c_light / vth
      integer, parameter :: Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
      integer, intent(in) :: k
      integer, intent(out) :: N
      real(kind=dp), dimension(Nlambda) :: line_u_grid
      type (atomicline), intent(in) :: line
      real(kind=dp) ::  vcore, vwing, v0, v1, vth, vB
      integer :: la, Nmid

      vth = vbroad(T(k), line%atom%weight, vturb(k))
      vB = 0.0_dp
      ! if (line%polarizable) then
      !    !replace b_char by magnetic field modulus here
      !    vB = B_char * LARMOR * (line%lambda0*NM_TO_M) * abs(line%g_lande_eff)
      ! else
      !    vB = 0.0_dp
      ! endif
      vwing = line%qwing * (vth + vB)

      vcore = 2.5 * vth!wing_to_core * vwing
  
      v0 = - vwing
      v1 = + vwing
      line_u_grid = 0.0_dp

  
      line_u_grid(Nlambda_line_w:1:-1) = -spanl_dp(vcore, vwing, Nlambda_line_w, -1)

      line_u_grid(Nlambda_line_w:Nlambda_line_c_log+Nlambda_line_w-1) = span_dp(-vcore, 0.0_dp, Nlambda_line_c_log, 1)
    
      Nmid = Nlambda/2 + 1
    
      line_u_grid(1:Nmid) = line_u_grid(1:Nmid)
      line_u_grid(Nmid+1:Nlambda) = -line_u_grid(Nmid-1:1:-1)

      line_u_grid(:) = line_u_grid(:) / vth

      N = nlambda     
 
      return
   end function line_u_grid
 
   
 
!    subroutine  define_local_gauss_profile_grid (atom)
!       type (AtomType), intent(inout) :: atom
!       real(kind=dp), dimension(2*(Nlambda_line_c_log+Nlambda_line_w-1)-1) :: vel
!       real(kind=dp) ::  vcore, vwing, v0, v1, vbroad
!       integer :: Nlambda, la, Nmid, kr
!  !     real, parameter :: fw = 7.0, fc = 4.0
!       real, parameter :: fw = 3.0, fc = 1.0
 
!       vbroad = maxval(atom%vbroad)
!       vwing = fw * vbroad
!       vcore = fc * vbroad
 
!       v0 = - vwing
!       v1 = + vwing
!       vel = 0.0_dp
 
 
!       vel(Nlambda_line_w:1:-1) = -spanl_dp(vcore, vwing, Nlambda_line_w, -1)
!       vel(Nlambda_line_w:Nlambda_line_c_log+Nlambda_line_w-1) = span_dp(-vcore, 0.0_dp, Nlambda_line_c_log, 1)
 
!       Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
 
!       Nmid = Nlambda/2 + 1
 
!       allocate(atom%ug(Nlambda)); atom%ug = 0.0
 
!       atom%ug(1:Nmid) = vel(1:Nmid)
!       atom%ug(Nmid+1:Nlambda) = -vel(Nmid-1:1:-1)
 
!       return
!    end subroutine define_local_gauss_profile_grid
  
 
   !TO DO: case without lines
!    subroutine make_wavelength_group()
 
!       integer :: Ntrans, Ncont, Nlines, Nlambda_cont
!       real(kind=dp) :: delta_v
!       integer :: n, kr, Nlam
!       type (Atom_Type), pointer :: atom
 
!       Ntrans = 0
!       delta_v = 0.0
 
!       !maximum and minimum wavelength for only lines, including max velocity field
!       !Count Number of transitions and number of lines
!       Nlines = 0
!       Nlambda_cont = 0
!       Ncont = 0
!       do n=1, N_atoms
!          atom => atoms(n)%p
!          do kr=1,atom%Ncont
!             Nlambda_cont = Nlambda_cont + atom%continua(kr)%Nlambda
!          enddo
!          Ntrans = Ntrans + atom%Ntr
      
!          Nlines = Nlines + atom%Nline
!          Ncont = Ncont + atom%Ncont
!       enddo
!       atom => null()
 
!     !However, it is very likely that the reference wavelength is not added at the moment
!     !if it falls in the group region
!      if (wl_ref > 0) then
!         Nlambda_cont = Nlambda_cont + 1
!         write(*,*) " Adding reference wavelength at ", wl_ref," nm!"
!      endif
 
!      allocate(cont_waves(Nlambda_cont), stat=alloc_status)
!      if (alloc_status>0) then
!         write(*,*) "Allocation error cont_waves"
!         stop
!      endif
 
!      if (wl_ref > 0) then
!         lac = 1
!         cont_waves(1) = wl_ref
!      else
!         lac = 0
!      endif
 
!      do n=1, Natom
!         atom => atoms(n)%ptr_atom
!         do kr=1,atom%Ncont
!            do la=1, atom%continua(kr)%Nlambda
!               lac = lac + 1
!               cont_waves(lac) = atom%continua(kr)%lambda(la)
!            enddo
!            !not used anymore
!            deallocate(atom%continua(kr)%lambda)
!            !lambda_file (and alpha_file) kept if allocated for explicit continua
!         enddo
!      enddo
 
!      !sort continuum frequencies
!      Nmore_cont_freq = 0.0
!      allocate(sorted_indexes(Nlambda_cont),stat=alloc_status)
!      if (alloc_status > 0) call error ("Allocation error sorted_indexes (cont)")
!      sorted_indexes = index_bubble_sort(cont_waves)
!      cont_waves(:) = cont_waves(sorted_indexes)
!      deallocate(sorted_indexes)
 
!      !remove duplicates
!      allocate(tmp_grid(Nlambda_cont), stat=alloc_status)
!      if (alloc_status > 0) call error ("Allocation error tmp_grid (cont)")
!      tmp_grid(2:Nlambda_cont) = 0.0
!      tmp_grid(1) = cont_waves(1)
!      Nremoved = 0
!      do la = 2, Nlambda_cont
!         if (cont_waves(la) > cont_waves(la-1)) then
!            tmp_grid(la) = cont_waves(la)
!         else
!            Nremoved = Nremoved + 1
!         endif
!      enddo
 
!      !write(*,*) "Total continuum frequencies, before merging : ", Nlambda_cont - Nremoved!lac
!      if (Nremoved > 0) then
!         write(*,*) " ->", Nremoved, " duplicate frequencies"
!         deallocate(cont_waves)
!         allocate(cont_waves(Nlambda_cont-Nremoved), stat=alloc_status)
!         cont_waves(:) = Pack(tmp_grid, tmp_grid > 0)
!      endif
!      deallocate(tmp_grid)
!      max_cont = maxval(cont_waves)!used to test if we need to add points below lines
!      Nlambda_cont = Nlambda_cont - Nremoved
 
!      !-> lines + cont
!      lthere_is_lines = (Nlam > 0)
!      if (lthere_is_lines) then
!         allocate(all_lamin(Nlam), all_lamax(Nlam), all_lam0(Nlam), stat=alloc_status)
!         if (alloc_status > 0) then
!            write(*,*) "Allocation error all_lam*"
!            stop
!         endif
 
!         !Store the maximum and minimum extent of each line including max velocity field
!         !lines can be red/blue-shifted at a max/min wavelength of red/blue *(1+/- delta_v/clight)
!         Nlam = 0
!         do n=1, Natom
!            atom => atoms(n)%ptr_atom
 
!            do kr=1,atom%Nline
!               Nlam = Nlam + 1
!               all_lamin(Nlam) = atom%lines(kr)%lambdamin * (1.0 - delta_v / clight)
!               all_lamax(Nlam) = atom%lines(kr)%lambdamax * ( 1.0 + 2*delta_v / clight)
!               all_lam0(Nlam) = atom%lines(kr)%lambda0
!            enddo
 
!         enddo
 
!         allocate(sorted_indexes(Nlam),stat=alloc_status)
!         if (alloc_status > 0) then
!            write(*,*) "Allocation error sorted_indexes(Nlam)"
!            stop
!         endif
!        !sort lines so that all_lamin(1) is always the first line
!        sorted_indexes = index_bubble_sort(all_lamin)
!        all_lamin(:) = all_lamin(sorted_indexes)
!        !->not necessarily ordered by min to max, but follows the order of lamin
!        !so that lmax(1) is associated to lamin(1)  which is important.
!        !If lines overlap, the lmax(1) could be associated to lmin(2) for instance.
!        all_lamax(:) = all_lamax(sorted_indexes)
!        deallocate(sorted_indexes)
 
!        group_blue(:) = -1.0
!        group_red(:) = -1.0
!        group_mid(:) = -1.0
 
!        Ngroup = 1
!        group_blue(Ngroup) = all_lamin(1)
!        group_red(Ngroup) = all_lamax(1)
!        !Find group of lines, and store for each group the lambda_blue and lambda_red of each group
!        !if a line overlaps with the previous line, add it to the same group and check the next line.
!        !Stop counting lines in a group if the next line does not overlap with the previous line. In
!        !the latter case, create a new group and start again.
!      ! Note: the first and last lines of a group may not overlap.
!        Nline_per_group(:) = 0
!        Nline_per_group(1) = 1
!        do Nlam = 2, size(all_lamin)
 
!            !Is the line overlapping the previous line ?
 
!            !Yes, add it to the same group
!            if (((all_lamin(Nlam) >= group_blue(Ngroup)).and.&
!             (all_lamin(Nlam) <= group_red(Ngroup))).or.&
!             ((all_lamax(Nlam) >= group_blue(Ngroup)).and.&
!             (all_lamax(Nlam) <= group_red(Ngroup)))) then
 
!               group_blue(Ngroup) = min(all_lamin(Nlam), group_blue(Ngroup))
!               group_red(Ngroup) = max(all_lamax(Nlam), group_red(Ngroup))
 
!               Nline_per_group(Ngroup) = Nline_per_group(Ngroup) + 1
 
!               !no, create a new group, starting with this line at first element
!            else
!               Ngroup = Ngroup + 1
!               if (Ngroup > MAX_GROUP_OF_LINES) then
!                  write(*,*) " Error, Ngroup > MAX_GROUP_OF_LINES", Ngroup
!                  stop
!               endif
!               group_blue(Ngroup) = all_lamin(Nlam)
!               group_red(Ngroup) = all_lamax(Nlam)
!               Nline_per_group(Ngroup) = 1
!            endif
 
!        enddo
 
!        write(*,*) " Found ", Ngroup, " groups of lines"
!  !       do Nlam=1, Ngroup
!  !       	write(*,*) " group ", Nlam, " lamin = ", group_blue(Nlam), ' lamax = ', group_red(Nlam), &
!  !       		" lamid = ", 0.5 * (group_blue(Nlam)+group_red(Nlam))
!  !       enddo
!  !       write(*,*) " -> ", sum(Nline_per_group), " lines"
!  !       write(*,*) " Found ", sum(Nline_per_group) - Ngroup, " overlapping regions for", sum(Nline_per_group), " lines"
!  !       write(*,*) " --> ", Ngroup, " groups of lines"
!        allocate(Nlambda_per_group(Ngroup), stat=alloc_status)
!        if (alloc_status > 0) then
!            write(*,*) "Allocation error Nlambda_per_group"
!            stop
!        endif
 
!        !sample each group with constant resol in km/s from group_blue(n) to group_red(n)
!        allocate(line_waves(1000000),stat=alloc_status)
!        if (alloc_status > 0) call error("Allocation error line_waves!")
!        line_waves = -1.0
 
!        !add all edges of each line ? + ref wavelength. Break the regularity of the grid ?
!        lac = 1
!  !       if (wl_ref > 0.0) then
!  !       	lac = 2
!  !       	line_waves(1) = wl_ref
!  !       else
!  !       	lac = 1
!  !       endif
!  !       do n=1, Natom
!  !       	atom => atoms(n)%ptr_atom
!  !
!  !       	do kr=1,atom%Nline
!  !         	line_waves(lac) = atom%lines(kr)%lambda0
!  !         	lac = lac + 1
!  !         	line_waves(lac) = atom%lines(kr)%lambdamin
!  !         	lac = lac + 1
!  !         	line_waves(lac) = atom%lines(kr)%lambdamax
!  !         	lac = lac + 1
!  !       	enddo
!  !
!  !       enddo
 
!        shift = lac!1
!        do n=1, Ngroup
 
!           line_waves(shift) = group_blue(n)
!  !       	write(*,*) "n=", n, " shift = ", shift
!           lac = 2
!           Nlambda_per_group(n) = 1
!           inf : do
!              la = lac + shift - 1
!              line_waves(la) = line_waves(la-1) * (1.0 + 1d3 * hv / clight)
!              Nlambda_per_group(n) = Nlambda_per_group(n) + 1
!              if (line_waves(la) >= group_red(n)) exit inf
!              if (la >= size(line_waves)) call error('la larger than size line_waves!')
!              lac = lac + 1
!  !       		write(*,*) lac, la, line_waves(la)
!           enddo inf
!           !change group_red!!
 
!  !       	write(*,*) 'gr change from to ', group_red(n), line_waves(la), la
!           group_red(n) = line_waves(la)!make the continuum below lines the same
 
!           shift = shift + Nlambda_per_group(n)
!  !       	write(*,*) ' Nlambda = ', Nlambda_per_group(n)
 
!        enddo
!  !       write(*,*) Nlambda_per_group
!        Nspec_line = sum(Nlambda_per_group)
!        line_waves = pack(line_waves,mask=line_waves > 0)
!        !write(*,*) " Nspec_line = ", Nspec_line, " check=", size(line_waves)
 
!        deallocate(all_lamin, all_lamax, all_lam0)
 
 
!        !-> Add continuum points beyond last continuum
!        !In case they are lines beyond the last continuum I add at least3 points per line for the continuum in this region
!        !->cont end		this is heavy for nothing but should work !
!        !finalise continuum here by reckoning how much freq we need
!        !make the new points go farther than line_waves for interpolation.
!        Nmore_cont_freq = 0
!        do n=1, Ngroup
!            l0 = group_blue(n)
!            l1 = group_red(n)
!            if (l0 > max_cont) then
!                  Nmore_cont_freq = Nmore_cont_freq + 1
!  !              	write(*,*) "adding more continuum (1)"
!  !              	write(*,*) l0, max_cont
!            endif
!            if (l1 > max_cont) then
!                  Nmore_cont_freq = Nmore_cont_freq + 1
!  !              	write(*,*) "adding more continuum (2)"
!  !              	write(*,*) l1, max_cont
!            endif
!            if (0.5*(l0+l1) > max_cont) then
!                  Nmore_cont_freq = Nmore_cont_freq + 1
!  !              	write(*,*) "adding more continuum (3)"
!  !              	write(*,*) 0.5*(l0+l1), max_cont
!            endif
 
!        enddo
 
!        check_new_freq = Nmore_cont_freq
!        if (Nmore_cont_freq > 0) then
!            write(*,*) "Adding new wavelength points for lines beyond continuum max!"
!            write(*,*) "  -> Adding ", Nmore_cont_freq," points"
!            allocate(tmp_grid(Nlambda_cont))
!            tmp_grid = cont_waves
!            deallocate(cont_waves)
!            allocate(cont_waves(Nlambda_cont + Nmore_cont_freq))
!            cont_waves(1:Nlambda_cont) = tmp_grid(:)
!            deallocate(tmp_grid)
!            allocate(tmp_grid(Nmore_cont_freq))
!            tmp_grid(:) = 0.0_dp
 
 
!            Nmore_cont_freq = 0
!            do n=1, Ngroup
!  !              	write(*,*) "n=",n
!                  l0 = group_blue(n)
!                  l1 = group_red(n)
!                  if (l0 > max_cont) then
!                     Nmore_cont_freq = Nmore_cont_freq + 1
!                     tmp_grid(Nmore_cont_freq) = l0
!  !              		write(*,*) "(1)", Nmore_cont_freq , "l0=",l0
!                  endif
!                  if (0.5*(l0+l1) > max_cont) then
!                     Nmore_cont_freq = Nmore_cont_freq + 1
!                     tmp_grid(Nmore_cont_freq) = 0.5 * (l0+l1)
!  !              		write(*,*) "(2)", Nmore_cont_freq , "lmid=", 0.5*(l0+l1)
!                  endif
!                  if (l1 > max_cont) then
!                     Nmore_cont_freq = Nmore_cont_freq + 1
!                     tmp_grid(Nmore_cont_freq) = l1
!  !              		write(*,*) "(3)", Nmore_cont_freq , "l1=",l1
!                  endif
 
!            enddo
!            if (Nmore_cont_freq /= check_new_freq) then
!                  call Warning("There are probably some frequency missing!")
!                  write(*,*) "Nmore_freq: ",check_new_freq," Nfreq_added: ", Nmore_cont_freq
!            endif
 
!            allocate(sorted_indexes(Nmore_cont_freq))
!            sorted_indexes(:) = index_bubble_sort(tmp_grid)
!            tmp_grid(:) = tmp_grid(sorted_indexes)
!            cont_waves(Nlambda_cont+1:Nlambda_cont + Nmore_cont_freq) = tmp_grid(:)
!            Nlambda_cont = Nlambda_cont + Nmore_cont_freq
!            deallocate(tmp_grid, sorted_indexes)
!            allocate(cont_grid(Nlambda_cont),stat=alloc_status)
!            if (alloc_status > 0) call error("Allocation error cont_grid")
!            cont_grid(:) = cont_waves(:)
!        else
!             allocate(cont_grid(Nlambda_cont),stat=alloc_status)
!             if (alloc_status > 0) call error("Allocation error cont_grid")
!             cont_grid(:) = cont_waves(:)
!        endif
!  !->cont end
 
!        !write(*,*) "bounds:"
!        !write(*,*) "  -> max lam:", maxval(line_waves), " max_cont:", maxval(cont_grid)
!            ! " max reddest line:", maxval(group_red,mask=group_red>-1),
 
!        !add lines + continua frequencies
!        Nspec_cont = size(cont_waves)
!        if (Nspec_cont /= Nlambda_cont) then
!              write(*,*) " Something went wrong with Nlambda cont"
!              stop
!        endif
 
!        !initiate with lines
!        allocate(tmp_grid(Nspec_line+Nspec_cont), stat=alloc_status)
!        if (alloc_status > 0) call error ("Allocation error tmp_grid")
!        tmp_grid(:) = -99
!        do la=1,Nspec_line
!              tmp_grid(la) = line_waves(la)
!        enddo
 
!        Nwaves = Nspec_line
!        !add continuum wavlengths (including reference wavelength), only outside line groups
!        !First values below or beyond first and last groups
!        la = 0
!        do lac=Nspec_line+1, Nspec_cont+Nspec_line
!              if ((cont_waves(lac-Nspec_line) < group_blue(1)) .or. (cont_waves(lac-Nspec_line) > group_red(Ngroup))) then
!                    tmp_grid(lac) = cont_waves(lac-Nspec_line)
!                    Nwaves = Nwaves + 1
!                    if (cont_waves(lac-Nspec_line) < group_blue(1)) la = lac
!              endif
!        enddo
 
!        !now values between groups
!        do lac=la+1, Nspec_cont+Nspec_line
!              group_loop : do n=2, Ngroup
!                    if ((cont_waves(lac-Nspec_line) > group_red(n-1)).and.(cont_waves(lac-Nspec_line) < group_blue(n))) then
!                          Nwaves = Nwaves + 1
!                          tmp_grid(lac) = cont_waves(lac-Nspec_line)
!                  !else
!                  ! be smart and cycle to accelerate
!                    endif
!              enddo group_loop
!        enddo
!        !write(*,*) " Check Nwaves:", Nwaves,  size(pack(tmp_grid, tmp_grid > 0))
!        deallocate(Nlambda_per_group, line_waves)
!        deallocate(cont_waves)
 
!        !continuum frequencies are sorted and so are the line frequencies
!        !but they are added at the end, so sorted is needed, but I can improve the previous
!        !loop to fill the tmp_frid in the ascending order of wavelengths
!        allocate(outgrid(Nwaves),stat=alloc_status)
!        tmp_grid = tmp_grid(index_bubble_sort(tmp_grid))
!        outgrid(:) = -99.0 !check
!        outgrid = pack(tmp_grid, mask=tmp_grid > 0)
 
!        do lac=2, Nwaves
!              if (outgrid(lac) <= outgrid(lac-1)) then
!                    write(*,*) lac, "lambda = ", outgrid(lac), outgrid(lac-1), minval(outgrid)
!                    call error("Sorted problem")
!              endif
!        enddo
 
!        deallocate(tmp_grid)
 
!      else !pure cont
!        Nspec_cont = size(cont_waves)
!        if (Nspec_cont /= Nlambda_cont) then
!              write(*,*) " Something went wrong with Nlambda cont"
!              stop
!        endif
!        allocate(cont_grid(Nlambda_cont),stat=alloc_status)
!        Nwaves = Nlambda_cont
!        if (alloc_status > 0) call error("Allocation error cont_grid")
!        cont_grid(:) = cont_waves(:)
!        write(*,*) "  -> pure cont max_cont:", maxval(cont_grid)
!        Nspec_line = 0
!        allocate(outgrid(Nwaves))
!        outgrid(:) = cont_grid(:)
!        deallocate(cont_waves)
!      endif !there is lines
 
 
!      write(*,*) Nwaves, " unique wavelengths" !they are no eliminated lines
!      write(*,*) Nspec_line, " line wavelengths"
!      write(*,*) Nspec_cont, " continuum wavelengths"
!      if (lthere_is_lines) then
!        write(*,*) "Mean number of lines per group:", real(sum(Nline_per_group))/real(Ngroup)
!        write(*,*) "Mean number of wavelengths per group:", real(Nspec_line)/real(Ngroup)
!        write(*,*) "Mean number of wavelengths per line:", real(Nspec_line)/real(Ntrans-Ncont)
!        write(*,*) "Resolution of line's groups (km/s):", hv
!      endif
!      write(*,*) "Mean number of wavelengths per continuum:", real(Nwaves - Nspec_line) / real(Ncont)
 
!      lam_unit = "nm"
!      l0 = minval(outgrid); l1 = maxval(outgrid)
!      if (l1 > 1e6) then
!        l1 = 10000000./l1
!        lam_unit = "cm^-1"
!      else if (l1 > 1e6) then
!        l1 = l1 * 1e-9 * 1e3
!        lam_unit = "mm"
!      else if (l1 > 1e7) then
!        l1 = l1 * 1e-9 * 1e2
!        lam_unit = "cm"
!      endif
!      write(*,'("Wavelength grid: "(1F12.4)" nm to",(1F12.4)" "(1A15))') l0, l1, lam_unit
 
!      !allocate indexes on the grid
!      Nlambda_max_line = 0
!      Nlambda_max_cont = 0
!      Nlambda_max_trans = 0
!      do n=1,Natom
!         atom => Atoms(n)%ptr_atom
!         do kr=1,atom%Ncont !only on the cont grid !
 
!            atom%continua(kr)%Nb = locate(outgrid, atom%continua(kr)%lambdamin)
!            atom%continua(kr)%Nr = locate(outgrid, atom%continua(kr)%lambdamax)
 
!            atom%continua(kr)%Nblue = locate(cont_grid, atom%continua(kr)%lambdamin)
!            atom%continua(kr)%Nred = locate(cont_grid, atom%continua(kr)%lambdamax)
!            atom%continua(kr)%Nmid = locate(cont_grid, 0.5*(atom%continua(kr)%lambdamin+atom%continua(kr)%lambdamax))
!            atom%continua(kr)%N0 = locate(cont_grid, atom%continua(kr)%lambda0)
!            atom%continua(kr)%Nlambda = atom%continua(kr)%Nred - atom%continua(kr)%Nblue + 1
 
!            !in any problem of grid resolution etc or locate approximation.
!            !We take Nred-1 to be sure than the lambda_cont(Nred) <= lambda0.
!            !Only if not dissolution.
!            !We just need to avoind having cont_lambda(Nred)>lambda0, since the cross section is in (lambda/lambda0)**3
!            if (cont_grid(atom%continua(kr)%N0) /= atom%continua(kr)%lambda0) then
!            !!write(*,*) " Beware, cont%lambda0 might not be on the grid", kr, atom%continua(kr)%i, atom%continua(kr)%j
!               if (.not.ldissolve) then
!                  if (cont_grid(atom%continua(kr)%Nred) > atom%continua(kr)%lambda0) then
!                     call Warning("continuum Nred larger than lambda0 !")
!                     write(*,*) atom%continua(kr)%lambda0, cont_grid(atom%continua(kr)%Nred)
!                     if (cont_grid(atom%continua(kr)%Nred-1) <= atom%continua(kr)%lambda0) then
!                        write(*,*) " ++++ adjusting Nred", " lambda0=",atom%continua(kr)%lambda0
!                        !To do while until <= lambda0
!                        atom%continua(kr)%Nred = atom%continua(kr)%Nred-1
!                        atom%continua(kr)%N0 = atom%continua(kr)%Nred
!                        atom%continua(kr)%Nlambda = atom%continua(kr)%Nred - atom%continua(kr)%Nblue + 1
!                        atom%continua(kr)%Nr = locate(outgrid,cont_grid(atom%continua(kr)%Nred))
!                        write(*,*) "new val at Nred:", outgrid(atom%continua(kr)%Nr), cont_grid(atom%continua(kr)%Nred)
!                     endif
!                  endif
!               endif
!            endif
!            !sur la grille totale
!            Nlambda_max_cont = max(Nlambda_max_cont,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1)
 
!         enddo
 
!        !The indexes on the grid do not include the maximum shift. So lambdamax/min += max shift < / > max/min lambda grid
!         do kr=1,atom%Nline
!            atom%lines(kr)%Nblue = locate(outgrid, atom%lines(kr)%lambdamin)
!            atom%lines(kr)%Nred = locate(outgrid, atom%lines(kr)%lambdamax)
!            atom%lines(kr)%Nmid = locate(outgrid, atom%lines(kr)%lambda0)
!            atom%lines(kr)%Nlambda = atom%lines(kr)%Nred - atom%lines(kr)%Nblue + 1
!            !write(*,*) "line", kr, " lam0=",atom%lines(kr)%lambda0, atom%lines(kr)%lambdamin, atom%lines(kr)%lambdamax
!            !write(*,*) " -> bounds on the grid:", outgrid(atom%lines(kr)%Nblue), outgrid(atom%lines(kr)%Nred)
!            !write(*,*) " -> max extent:", outgrid(atom%lines(kr)%Nblue)*(1.0 - delta_v/clight), outgrid(atom%lines(kr)%Nred)*(1.0 + delta_v/clight)
!            if (atom%lines(kr)%Nblue-dshift < 1) then
!               write(*,*) "Error for line ", kr, " of atom ", atom%ID
!               write(*,*) " Nblue below 1"
!               write(*,*) "Nb=",atom%lines(kr)%Nblue, " shift=",-dshift
!               write(*,*) " -> sum :", atom%lines(kr)%Nred-nint(1e-3 * delta_v/hv)
!               write(*,*) " lambdamin (CMF) = ", atom%lines(kr)%lambdamin
!               write(*,*) " max(lambda) grid  = ", minval(outgrid)
!               write(*,*) outgrid(atom%lines(kr)%Nblue)*(1.0 - delta_v/clight)
!               stop
!            else if(atom%lines(kr)%Nred+dshift > Nwaves) then
!               write(*,*) "Error for line ", kr, " of atom ", atom%ID
!               write(*,*) " Nred larger than Nwaves", Nwaves, size(outgrid)
!               write(*,*) "Nr=",atom%lines(kr)%Nred, " shift=",dshift
!               write(*,*) " -> sum :", atom%lines(kr)%Nred+nint(1e-3 * delta_v/hv)
!               write(*,*) " lambdamax (CMF) = ", atom%lines(kr)%lambdamax
!               write(*,*) " max(lambda) grid  = ", maxval(outgrid)
!               write(*,*) outgrid(atom%lines(kr)%Nred)*(1.0 + delta_v/clight)
!               write(*,*) "reducing shift"
!               dshift = Nwaves - atom%lines(kr)%Nred
!               write(*,*) dshift
!  !              stop
!            endif
!            !does not take the shift into account
!            Nlambda_max_line = max(Nlambda_max_line, atom%lines(kr)%Nlambda)
!         enddo
 
!         atom => NULL()
!      enddo
!      write(*,*) "Number of max freq points for all lines at this resolution :", Nlambda_max_line
!      write(*,*) "Number of max freq points for all cont at this resolution :", Nlambda_max_cont
!      !takes the shift into account, maximum size needed for a frequency array for line opacity.
!  !     Nlambda_max_trans = max(Nlambda_max_line+2*int( sign(1.0_dp, delta_v) * ( 1e-3 * abs(delta_v) / hv + 0.5 ) ),Nlambda_max_cont)
!      Nlambda_max_trans = max(Nlambda_max_line+2*dshift,Nlambda_max_cont)
!      write(*,*) "Number of max freq points for all trans at this resolution :", Nlambda_max_trans
!  ! stop
!      return
!    end subroutine make_wavelength_group
 
 
 end module wavelengths_gas
 
