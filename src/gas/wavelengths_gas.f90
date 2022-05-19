module wavelengths_gas

!TO DO:
   ! merge all wavelengths grids

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
   integer, parameter :: Nlambda_cont = 101
   integer, parameter :: Nlambda_cont_log = 31
   integer, parameter :: Nlambda_line_w = 24
   integer, parameter :: Nlambda_line_c_log = 14
   real               :: hv
   !number max of lambda for all lines
   integer            :: Nlambda_max_line, Nlambda_max_cont, Nlambda_max_trans

   !group of overlapping transitions
   integer :: N_groups, n_lambda_cont
   real(kind=dp), dimension(:), allocatable :: group_blue, group_red
   integer, dimension(:), allocatable :: Nline_per_group, Nlambda_per_group
   !non-lte freq grid
   real(kind=dp), allocatable :: tab_lambda_cont(:), tab_lambda_nm(:)

   !generator of grid for specific transitions
   procedure(make_sub_wavelength_grid_cont_log_nu), pointer :: subgrid_cont => make_sub_wavelength_grid_cont_log_nu
   procedure(line_lambda_grid), pointer :: subgrid_line => line_lambda_grid!_dv
 
   contains

   function make_sub_wavelength_grid_cont_log_nu(cont)
   ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   ! ----------------------------------------------------------------- !
      real(kind=dp), dimension(Nlambda_cont_log) :: make_sub_wavelength_grid_cont_log_nu
      type (AtomicContinuum), intent(in) :: cont
      real(kind=dp) :: resol
      integer :: la, N1, N2
      real(kind=dp) :: l0, l1
      real(kind=dp) :: nu1, nu0, nu2, dnu

      l1 = cont%lambdamax
      l0 = cont%lambdamin

      if (cont%lambdamax > cont%lambda0) then
         N2 = Nlambda_cont_log + 1
         N1 = Nlambda_cont_log
      else
         N2 = 0
         N1 = Nlambda_cont_log
      endif


      nu0 = (M_TO_NM * C_LIGHT / cont%lambda0)/1e15
      nu1 = (M_TO_NM * C_LIGHT / l0)/1e15
      nu2 = (M_TO_NM * C_LIGHT / l1)/1e15


      make_sub_wavelength_grid_cont_log_nu(N1:1:-1) = 1e-15 * c_light / spanl_dp(nu0,nu1,N1,1) * m_to_nm

      if (N2 > 0) then
         dnu = abs(make_sub_wavelength_grid_cont_log_nu(N1)-make_sub_wavelength_grid_cont_log_nu(N1-1))
         make_sub_wavelength_grid_cont_log_nu(N2+N1:N1+1:-1) = 1e-15 * c_light / spanl_dp(nu1+dnu,nu2,N2,1) * m_to_nm
      endif

      return
   end function make_sub_wavelength_grid_cont_log_nu

   function line_u_grid_loc(k, line, N)
      !return for line line and cell k, the paramer:
      ! u = (lambda - lambda0) / lambda0 * c_light / vth
      integer, parameter :: Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
      integer, intent(in) :: k
      integer, intent(out) :: N
      real(kind=dp), dimension(Nlambda) :: line_u_grid_loc
      type (atomicline), intent(in) :: line
      real(kind=dp) ::  vcore, vwing, vth, vB
      integer :: la, Nmid

      vth = vbroad(T(k), line%atom%weight, vturb(k))
      vB = 0.0_dp
      if (line%polarizable) then
         !replace b_char by magnetic field modulus here
         vB = B_char * LARMOR * (line%lambda0*NM_TO_M) * abs(line%g_lande_eff)
      else
         vB = 0.0_dp
      endif
      vwing = line%qwing * (vth + vB)

      vcore = 2.5 * vth!wing_to_core * vwing

      line_u_grid_loc = 0.0_dp

  
      line_u_grid_loc(Nlambda_line_w:1:-1) = -spanl_dp(vcore, vwing, Nlambda_line_w, -1)

      line_u_grid_loc(Nlambda_line_w:Nlambda_line_c_log+Nlambda_line_w-1) = span_dp(-vcore, 0.0_dp, Nlambda_line_c_log, 1)
    
      Nmid = Nlambda/2 + 1
    
      line_u_grid_loc(1:Nmid) = line_u_grid_loc(1:Nmid)
      line_u_grid_loc(Nmid+1:Nlambda) = -line_u_grid_loc(Nmid-1:1:-1)

      line_u_grid_loc(:) = line_u_grid_loc(:) / vth

      N = nlambda     
 
      return
   end function line_u_grid_loc

   subroutine compute_line_bound(line,limage)
      ! ------------------------------------------------------------ !
      ! Compute the line bounds : lamndamin and lambdamax
      ! the total extent of line in the atom's frame (no velocity).
      ! ------------------------------------------------------------ !
      type (AtomicLine), intent(inout) :: line
      logical, intent(in) :: limage
      real(kind=dp) :: vB, vmax, vth
      
         
      if (line%polarizable) then
         vB = B_char * LARMOR * (line%lambda0*NM_TO_M) * abs(line%g_lande_eff)
      else
         vB = 0.0_dp
      endif
      vth = vbroad(maxval(T), line%atom%weight, maxval(vturb))
      vmax = line%qwing * (vth + vB)
      if (limage) vmax = line%atom%vmax_rt*1d3
       
       
      line%lambdamin = line%lambda0*(1.0-vmax/C_LIGHT)
      line%lambdamax = line%lambda0*(1.0+vmax/C_LIGHT)
      line%vmax = vmax

      !avoid to print too much info. all lines are kept for lmodel_1d spectra.
      !to remove after debug ?
      if (limage.and..not.lmodel_1d) then
         write(*,'("line "(1I2)"->"(1I2)"; Vmax (Ray-Trace)="(1F12.3)" km/s")') line%j, line%i, real(line%vmax)*1d-3
         write(*,'(" lamin="(1F12.3)" lam0="(1F12.3)" lamax="(1F12.3))') line%lambdamin, line%lambda0, line%lambdamax
      ! else
      !    write(*,'("line "(1I2)"->"(1I2)"; Vmax (non-LTE)="(1F12.3)" km/s")') line%j, line%i, real(line%vmax)*1d-3
      !    write(*,'(" lamin="(1F12.3)" lam0="(1F12.3)" lamax="(1F12.3))') line%lambdamin, line%lambda0, line%lambdamax
      endif

      return
   end subroutine compute_line_bound

   function line_lambda_grid(line,Nlambda)
      !integer, parameter :: Nlambda = (2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1)
      integer, intent(in) :: nlambda
      real(kind=dp), dimension(Nlambda) :: line_lambda_grid
      type (atomicline), intent(inout) :: line
      real(kind=dp) ::  vcore, vwing, v0, v1, vth, vB
      integer :: la, Nmid

      !not the exact width at a given cell, but the maximum extent!
      vth = vbroad(maxval(T), line%atom%weight, maxval(vturb))
      vB = 0.0_dp
      ! if (line%polarizable) then
      !    !replace b_char by magnetic field modulus here
      !    vB = B_char * LARMOR * (line%lambda0*NM_TO_M) * abs(line%g_lande_eff)
      ! else
      !    vB = 0.0_dp
      ! endif
      vwing = line%qwing * (vth + vB)
      vwing = c_light * (line%lambdamax - line%lambda0) / line%lambda0 !line%vmax

      vcore = 2.5 * vth!wing_to_core * vwing
  
      v0 = - vwing
      v1 = + vwing
      line_lambda_grid = 0.0_dp

  
      line_lambda_grid(Nlambda_line_w:1:-1) = -spanl_dp(vcore, vwing, Nlambda_line_w, -1)

      line_lambda_grid(Nlambda_line_w:Nlambda_line_c_log+Nlambda_line_w-1) = span_dp(-vcore, 0.0_dp, Nlambda_line_c_log, 1)
    
      Nmid = Nlambda/2 + 1
    
      line_lambda_grid(1:Nmid) = line_lambda_grid(1:Nmid)
      line_lambda_grid(Nmid+1:Nlambda) = -line_lambda_grid(Nmid-1:1:-1)

      line_lambda_grid = (1.0 + line_lambda_grid/c_light) * line%lambda0

      ! line%lambdamin = (1.0 + minval(line_lambda_grid)/c_light) * line%lambda0
      ! line%lambdamax = (1.0 + maxval(line_lambda_grid)/c_light) * line%lambda0
 
      return
   end function line_lambda_grid

   function line_lambda_grid_dv(line,Nlambda)
      !integer, parameter :: Nlambda = nint( (2*line%vmax)/hv + 1 )
      integer, intent(in) :: nlambda
      real(kind=dp), dimension(Nlambda) :: line_lambda_grid_dv
      type (atomicline), intent(inout) :: line
      real(kind=dp) ::  vwing, vth
      integer :: la, Nmid

      !not the exact width at a given cell, but the maximum extent!
 
      vwing = c_light * (line%lambdamax - line%lambda0) / line%lambda0
      vwing = line%vmax

      line_lambda_grid_dv = 0.0_dp
      line_lambda_grid_dv(1) = line%lambdamin
      do la=2, nlambda

         line_lambda_grid_dv(la) = line_lambda_grid_dv(la-1) * (1.0 + 1d3 * hv/c_light)

      enddo
 
      return
   end function line_lambda_grid_dv

   subroutine deallocate_wavelengths_gasrt(lambda)
      real(kind=dp), intent(inout), allocatable, dimension(:) :: lambda
      if (allocated(lambda)) deallocate(lambda)
      if (allocated(tab_lambda_nm)) deallocate(tab_lambda_cont, tab_lambda_nm)
      !deallocate atom grid..
      !deallocate group..
      if (allocated(group_blue)) deallocate(group_blue, group_red)
      if (allocated(Nline_per_group)) deallocate(Nline_per_group, Nlambda_per_group)
      return
   end subroutine deallocate_wavelengths_gasrt

   subroutine make_wavelengths_nlte(lambda,Vmax_overlap)
      !Total frequency grid, concatenation of all
      !individual grids.
      ! used only for the non-LTE loop.
      real(kind=dp), optional :: vmax_overlap
      real(kind=dp), intent(out), allocatable :: lambda(:)
      integer :: Ntrans, Ncont, Nlines, Nlambda_cont, Nlambda
      integer :: n, kr, Nlam, Nmore_cont_freq, Nremoved, Nwaves, check_new_freq
      type (AtomType), pointer :: atom
      integer :: alloc_status, lac, la, nb, krr, max_Nlambda_indiv
      real(kind=dp), dimension(:), allocatable :: tmp_grid, tmp_grid2, all_l0, all_l1
      integer, parameter :: Ngroup_max = 1000
      real(kind=dp), dimension(Ngroup_max) :: group_blue_tmp, group_red_tmp, Nline_per_group_tmp
      integer, dimension(:), allocatable :: sorted_indexes
      real(kind=dp) :: max_cont, l0, l1, dvmin
      real :: hv_loc !like hv
      real, parameter :: prec_vel = 0.0 !remove points is abs(hv_loc-hv)>prec_vel
      logical :: check_for_overlap

      check_for_overlap = .false.
      if (present(vmax_overlap)) then
         check_for_overlap = (abs(vmax_overlap) > 0.0) 
         if (check_for_overlap) write(*,*) "*** Searching for maximum overlap between lines!"
      endif
 
      Ntrans = 0
      !Compute number of transitions in total, the number of wavelengths for the continua
      ! and computes the optimal resolution for all lines including all atoms.
      Nlines = 0
      Nlambda_cont = 0
      Ncont = 0
      hv = 1d30
      do n=1, N_atoms
         atom => atoms(n)%p
         !km/s
         hv = min(hv, 1d-3 * 0.46*vbroad(minval(T,mask=T>0),atom%weight,minval(vturb,mask=T>0)))
         do kr=1,atom%Ncont
            if (atom%continua(kr)%hydrogenic) then
               atom%continua(kr)%Nlambda = Nlambda_cont_log
               if (atom%continua(kr)%lambdamax > atom%continua(kr)%lambda0) then
                  atom%continua(kr)%Nlambda = atom%continua(kr)%Nlambda + Nlambda_cont_log + 1
               endif               
            endif
            Nlambda_cont = Nlambda_cont + atom%continua(kr)%Nlambda
         enddo
         do kr=1,atom%Nline
            if (associated(subgrid_line,line_lambda_grid)) then
               atom%lines(kr)%Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
            ! elseif (associated(subgrid_line,line_lambda_grid_dv)) then
            !    atom%lines(kr)%Nlambda = nint(2 * line%vmax / hv + 1)
            ! else
            !    call error("wrong association for subline_grid!")
            endif
         enddo 
         Ntrans = Ntrans + atom%Ntr
      
         Nlines = Nlines + atom%Nline
         Ncont = Ncont + atom%Ncont
      enddo
      atom => null()
      if (associated(subgrid_line,line_lambda_grid_dv)) then
         if (art_hv > 0.0) then
            write(*,'("R="(1F7.3)" km/s; optimal:"(1F7.3)" km/s")') art_hv, hv
            hv = art_hv
         else
            write(*,'("R="(1F7.3)" km/s")') hv        
         endif
      else
         hv = 0.0
      endif
      ! ********************** Pure continuum RT ********************** !
      !Even if continuum RT is removed, that part is important
      allocate(tab_lambda_cont(Nlambda_cont), stat=alloc_status)
      if (alloc_status>0) then
         call error("Allocation error cont_waves")
      endif

      lac = 1
      do n=1, N_atoms
         atom => atoms(n)%p
         do kr=1,atom%Ncont
            if (atom%continua(kr)%hydrogenic) then
               tab_lambda_cont(lac:(lac-1)+atom%continua(kr)%Nlambda) = subgrid_cont(atom%continua(kr))
            else 
               tab_lambda_cont(lac:lac-1+atom%continua(kr)%Nlambda) = atom%continua(kr)%lambda_file(:)
            endif
            lac = lac + atom%continua(kr)%Nlambda
        enddo
        atom => null()
      enddo
      ! write(*,*) "lac end=", lac-1
 
      !sort continuum frequencies
      Nmore_cont_freq = 0.0
      allocate(sorted_indexes(Nlambda_cont),stat=alloc_status)
      if (alloc_status > 0) call error ("Allocation error sorted_indexes (cont)")
      sorted_indexes = index_bubble_sort(tab_lambda_cont)
      tab_lambda_cont(:) = tab_lambda_cont(sorted_indexes)
      deallocate(sorted_indexes)
 
      !remove duplicates
      allocate(tmp_grid(Nlambda_cont), stat=alloc_status)
      if (alloc_status > 0) call error ("Allocation error tmp_grid (cont)")
      tmp_grid(2:Nlambda_cont) = 0.0
      tmp_grid(1) = tab_lambda_cont(1)
      Nremoved = 0
      do la = 2, Nlambda_cont
         if (tab_lambda_cont(la) > tab_lambda_cont(la-1)) then
            tmp_grid(la) = tab_lambda_cont(la)
         else
            Nremoved = Nremoved + 1
         endif
      enddo
 
      if (Nremoved > 0) then
         write(*,*) " ->", Nremoved, " duplicate frequencies"
         deallocate(tab_lambda_cont)
         allocate(tab_lambda_cont(Nlambda_cont-Nremoved), stat=alloc_status)
         tab_lambda_cont(:) = Pack(tmp_grid, tmp_grid > 0)
      endif
      deallocate(tmp_grid)
      max_cont = maxval(tab_lambda_cont)
      Nlambda_cont = Nlambda_cont - Nremoved
      ! ********************** ***************** ********************** !

      ! ********************** cont + line    RT ********************** !
      if (Nlines > 0) then !(.true.) then!lines are present in some atoms
         allocate(all_l0(Nlines), all_l1(Nlines), stat=alloc_status)

         Nlam = 0
         do n=1, N_atoms
            atom => atoms(n)%p
 
            do kr=1,atom%Nline
               Nlam = Nlam + 1
               all_l0(Nlam) = atom%lines(kr)%lambdamin
               all_l1(Nlam) = atom%lines(kr)%lambdamax
            enddo
 
         enddo
         atom => null()
 
         allocate(sorted_indexes(Nlam),stat=alloc_status)
         if (alloc_status > 0) then
            call error("Allocation error sorted_indexes(Nlam)")
         endif
         sorted_indexes = index_bubble_sort(all_l0)
         !I sort the lambdamin, and lambdamax array is ordered so that
         !the index in all_l0 and in all_l1 correspond to the same pair of 
         !lambda min/max
         all_l0(:) = all_l0(sorted_indexes)
         all_l1(:) = all_l1(sorted_indexes)
         deallocate(sorted_indexes)

         group_blue_tmp(:) = -1.0
         group_red_tmp(:) = -1.0
 
         N_groups = 1
         group_blue_tmp(N_groups) = all_l0(1)
         group_red_tmp(N_groups) = all_l1(1)
         !Find group of lines, and store for each group the lambda_blue and lambda_red of each group
         !if a line overlaps with the previous line, add it to the same group and check the next line.
         !Stop counting lines in a group if the next line does not overlap with the previous line. In
         !the latter case, create a new group and start again.
         ! Note: the first and last lines of a group may not overlap.
         Nline_per_group_tmp(:) = 0
         Nline_per_group_tmp(N_groups) = 1
         do Nlam = 2, Nlines
 
            !Is the line overlapping the previous line ?
 
            !Yes, add it to the same group
            if (((all_l0(Nlam) >= group_blue_tmp(N_groups)).and.&
               (all_l0(Nlam) <= group_red_tmp(N_groups))).or.&
               ((all_l1(Nlam) >= group_blue_tmp(N_groups)).and.&
               (all_l1(Nlam) <= group_red_tmp(N_groups)))) then
 
               group_blue_tmp(N_groups) = min(all_l0(Nlam), group_blue_tmp(N_groups))
               group_red_tmp(N_groups) = max(all_l1(Nlam), group_red_tmp(N_groups))
 
               Nline_per_group_tmp(N_groups) = Nline_per_group_tmp(N_groups) + 1
 
            !no, create a new group, starting with this line at first element
            else
               N_groups = N_groups + 1
               if (N_groups > Ngroup_max) then
                  call error(" Error, Ngroup > Ngroup_max")
               endif
               group_blue_tmp(N_groups) = all_l0(Nlam)
               group_red_tmp(N_groups) = all_l1(Nlam)
               Nline_per_group_tmp(N_groups) = 1
           endif
 
         enddo !end loop to create group
         allocate(Nlambda_per_group(N_groups), group_red(N_groups), &
         group_blue(N_groups), Nline_per_group(N_groups), stat=alloc_status)
         if (alloc_status > 0) then
           call error("Allocation error groups")
         endif
         group_blue(:) = group_blue_tmp(1:N_groups)
         group_red(:) = group_red_tmp(1:N_groups)
         Nline_per_group(:) = Nline_per_group_tmp(1:N_groups)

         !some statistics
         write(*,*) " Found ", N_groups, " groups of lines for ", Nlines, " lines in total."
         write(*,*) "   -> ", sum(Nline_per_group) - N_groups, " overlapping regions"! for", sum(Nline_per_group), " lines"
         ! do Nlam=1, N_groups
         !    write(*,*) " group ", Nlam, " lamin = ", group_blue(Nlam), ' lamax = ', group_red(Nlam)
         ! enddo

         allocate(tmp_grid(1000000),stat=alloc_status)
         if (alloc_status > 0) call error("Allocation error tmp_grid (line)!")
         tmp_grid = -1.0
 
         !groups are not really used here
         lac = 1
         max_Nlambda_indiv = 0
         do n=1, N_atoms
            atom => atoms(n)%p
 
            do kr=1,atom%Nline

               if (hv>0.0) then
                  atom%lines(kr)%Nlambda = nint(2 * 1d-3 * atom%lines(kr)%vmax / hv + 1)
                  ! write(*,*) "Nlambda line hv>0:", nint(2 * 1d-3 * atom%lines(kr)%vmax / hv + 1)
               ! else
                  ! write(*,*) "Nlambda line:", atom%lines(kr)%Nlambda
               endif
               max_Nlambda_indiv = max(max_Nlambda_indiv,atom%lines(kr)%Nlambda)
               tmp_grid(lac:atom%lines(kr)%Nlambda+lac-1) = subgrid_line(atom%lines(kr),atom%lines(kr)%Nlambda)
               lac = lac + atom%lines(kr)%Nlambda

            enddo

         enddo
         atom => null()
         tmp_grid = pack(tmp_grid,tmp_grid > 0)
         Nlambda = size(tmp_grid)
         allocate(sorted_indexes(Nlambda),stat=alloc_status)
         if (alloc_status > 0) call error ("Allocation error sorted_indexes (line)")
         sorted_indexes = index_bubble_sort(tmp_Grid)
         tmp_grid(:) = tmp_grid(sorted_indexes)
         deallocate(sorted_indexes)
         !count the number of wavelengths per groups ? 
         !remove duplicates or points below resolution.
         allocate(tmp_grid2(Nlambda), stat=alloc_status)
         if (alloc_status > 0) call error ("Allocation error tmp_grid2")
         tmp_grid2(2:Nlambda) = 0.0
         tmp_grid2(1) = tmp_grid(1)
         Nremoved = 0

         !here we will remove points that are also below hv so as to preserve the grid resolutin.
         !however we keep points with hv_loc above hv as it means we are comparing lambda of different
         !line groups.
         !TO DO: check that when hv_loc >> hv it is indeed  because we change group.
         if (hv>0.0) then
            do la = 2, Nlambda
               hv_loc = ceiling ( real (1d-3 * c_light * (tmp_grid(la)-tmp_grid(la-1)) / tmp_grid(la) ) )
               if (abs(hv_loc - hv) >= prec_vel) then !hv or more.
                  tmp_grid2(la) = tmp_grid(la)
               else!below grid resolution, remove.
                  Nremoved = Nremoved + 1
               endif
            enddo
         else
            do la = 2, Nlambda
               if (tmp_grid(la) > tmp_grid(la-1)) then
                  tmp_grid2(la) = tmp_grid(la)
               else
                  Nremoved = Nremoved + 1
               endif
            enddo
         endif
         if (Nremoved > 0) then
            write(*,*) " ->", Nremoved, " duplicate frequencies in lines"
            deallocate(tmp_grid)
            ! allocate(tmp_grid(Nlambda-Nremoved), stat=alloc_status)
            tmp_grid(:) = Pack(tmp_grid2, tmp_grid2 > 0)
         endif
         deallocate(tmp_grid2)
         Nlambda = size(tmp_grid)!Nlambda - Nremoved
         ! do n=1,N_groups
         !    write(*,*) n, size(pack(tmp_grid,&
         !    (tmp_grid >= group_blue(n)).and.tmp_grid <= group_red(n)))
         !    Nlambda_per_group(n) = size(pack(tmp_grid,&
         !       (tmp_grid >= group_blue(n)).and.tmp_grid <= group_red(n)))
         ! enddo
         ! write(*,*) " ** BUG HERE ?? ** "
         ! write(*,*) "There are in total ", Nlambda, " wavelength points for the lines grid.", sum(Nlambda_per_group)
         ! write(*,*) "** ?? ** "

         ! ----------- small cont grid start ------------- !
         !-> This is only useful for the pure-continuum RT !!!!
         !and should disappear at some point.
         !-> Add continuum points beyond last "bound-free" continuum
         !In case they are lines beyond the last continuum I add at least3 points per line for the continuum in this region
         !->cont end		this is heavy for nothing but should work !
         !finalise continuum here by reckoning how much freq we need
         !make the new points go farther than line_waves for interpolation.
         Nmore_cont_freq = 0
         do n=1, N_groups
            l0 = group_blue(n)
            l1 = group_red(n)
            if (l0 > max_cont) then
               Nmore_cont_freq = Nmore_cont_freq + 1
            endif
            if (l1 > max_cont) then
               Nmore_cont_freq = Nmore_cont_freq + 1
            endif
            if (0.5*(l0+l1) > max_cont) then
               Nmore_cont_freq = Nmore_cont_freq + 1
            endif
 
         enddo
 
         check_new_freq = Nmore_cont_freq
         if (Nmore_cont_freq > 0) then
            write(*,*) "Adding new wavelength points for lines beyond continuum max!"
            write(*,*) "  -> Adding ", Nmore_cont_freq," points"
            write(*,*) "max cont, max line", max_cont, maxval(tmp_grid)
            allocate(tmp_grid2(Nlambda_cont))
            tmp_grid2 = tab_lambda_cont
            deallocate(tab_lambda_cont)
            allocate(tab_lambda_cont(Nlambda_cont + Nmore_cont_freq))
            tab_lambda_cont(1:Nlambda_cont) = tmp_grid2(:)
            deallocate(tmp_grid2)
            allocate(tmp_grid2(Nmore_cont_freq))
            tmp_grid2(:) = 0.0_dp
 
            Nmore_cont_freq = 0
            do n=1, N_groups
               l0 = group_blue(n)
               l1 = group_red(n)
               if (l0 > max_cont) then
                  Nmore_cont_freq = Nmore_cont_freq + 1
                  tmp_grid2(Nmore_cont_freq) = l0
               endif
               if (0.5*(l0+l1) > max_cont) then
                  Nmore_cont_freq = Nmore_cont_freq + 1
                  tmp_grid2(Nmore_cont_freq) = 0.5 * (l0+l1)
               endif
               if (l1 > max_cont) then
                  Nmore_cont_freq = Nmore_cont_freq + 1
                  tmp_grid2(Nmore_cont_freq) = l1
               endif
            enddo
            if (Nmore_cont_freq /= check_new_freq) then
               call Warning("There are probably some frequency missing!")
               write(*,*) "Nmore_freq: ",check_new_freq," Nfreq_added: ", Nmore_cont_freq
            endif
 
            allocate(sorted_indexes(Nmore_cont_freq))
            sorted_indexes(:) = index_bubble_sort(tmp_grid2)
            tmp_grid2(:) = tmp_grid2(sorted_indexes)
            tab_lambda_cont(Nlambda_cont+1:Nlambda_cont + Nmore_cont_freq) = tmp_grid2(:)
            Nlambda_cont = Nlambda_cont + Nmore_cont_freq
            deallocate(tmp_grid2, sorted_indexes)
         endif
 
         if (size(tab_lambda_cont) /= Nlambda_cont) then
             write(*,*) " Something went wrong with Nlambda cont"
             stop
         endif

         ! ----------- small cont grid end ------------- !

         !initiate with lines
         allocate(tmp_grid2(nlambda+Nlambda_cont), stat=alloc_status)
         if (alloc_status > 0) call error ("Allocation error tmp2_grid (final)")
         tmp_grid2(:) = -99
         tmp_grid2(1:Nlambda) = tmp_grid(:)
 
         Nwaves = Nlambda
         !add continuum wavlengths only outside line groups
         !First values below or beyond first and last groups
         la = 0
         do lac=Nlambda+1, Nlambda+Nlambda_cont
            if ((tab_lambda_cont(lac-Nlambda) < group_blue(1)) .or. (tab_lambda_cont(lac-Nlambda) > group_red(N_groups))) then
                   tmp_grid2(lac) = tab_lambda_cont(lac-Nlambda)
                   Nwaves = Nwaves + 1
                   if (tab_lambda_cont(lac-Nlambda) < group_blue(1)) la = lac
             endif
         enddo
 
         !now values between groups
         do lac=la+1, Nlambda_cont+Nlambda
            group_loop : do n=2, N_groups
               if ((tab_lambda_cont(lac-Nlambda) > group_red(n-1)).and.(tab_lambda_cont(lac-Nlambda) < group_blue(n))) then
                  Nwaves = Nwaves + 1
                  tmp_grid2(lac) = tab_lambda_cont(lac-nlambda)
                 !else
                 ! be smart and cycle to accelerate
               endif
            enddo group_loop
         enddo
         ! write(*,*) " Check Nwaves:", Nwaves,  size(pack(tmp_grid2, tmp_grid2 > 0))
         deallocate(tmp_grid)
         !
         ! If I remove continuum RT !
         !
         ! deallocate(tab_lambda_cont)
 
         !continuum frequencies are sorted and so are the line frequencies
         !but they are added at the end, so sorted is needed, but I can improve the previous
         !loop to fill the tmp_frid in the ascending order of wavelengths
         allocate(lambda(Nwaves),stat=alloc_status)
         tmp_grid2 = tmp_grid2(index_bubble_sort(tmp_grid2))
         lambda(:) = -99.0 !check
         lambda = pack(tmp_grid2, mask=tmp_grid2 > 0)
 
         do lac=2, Nwaves
            if (lambda(lac) <= lambda(lac-1)) then
               write(*,*) lac, "lambda = ", lambda(lac), lambda(lac-1), minval(lambda)
               call error("Sorted problem")
            endif
         enddo
         deallocate(tmp_grid2)
 
      else !pure cont
         Nwaves = Nlambda_cont
         allocate(lambda(Nwaves),stat=alloc_status)
         if (alloc_status>0) call error("allocation error lambda, in pure cont!")
         lambda = tab_lambda_cont
         deallocate(tab_lambda_cont)
         Nlambda = 0
      endif !there is lines

      write(*,*) Nwaves, " unique wavelengths" !they are no eliminated lines
      write(*,*) Nlambda, " line wavelengths"
      write(*,*) Nlambda_cont, " continuum wavelengths"
      write(*,*) max_Nlambda_indiv, " max individual wavelengths for a line" 
      ! if (Nlines>1) then
      !    write(*,*) "Mean number of lines per group:", real(sum(Nline_per_group))/real(Ngroup)
      !    write(*,*) "Mean number of wavelengths per group:", real(Nspec_line)/real(Ngroup)
      !    write(*,*) "Mean number of wavelengths per line:", real(Nspec_line)/real(Ntrans-Ncont)
      !    write(*,*) "Resolution of line's groups (km/s):", hv
      ! endif 
      ! ************************** ************ ************************* !

      !Now indexes of each transition on the lambda grid
      Nlambda_max_line = 0
      Nlambda_max_cont = 0
      Nlambda_max_trans = 0
      do n=1,N_atoms
         atom => Atoms(n)%p
         do kr=1,atom%Ncont

            atom%continua(kr)%Nbc = locate(tab_lambda_cont, atom%continua(kr)%lambdamin)
            atom%continua(kr)%Nrc = locate(tab_lambda_cont, atom%continua(kr)%lambdamax)
            atom%continua(kr)%Nlambdac = atom%continua(kr)%Nrc - atom%continua(kr)%Nbc + 1

            atom%continua(kr)%Nb = locate(lambda, atom%continua(kr)%lambdamin)
            atom%continua(kr)%Nr = locate(lambda, atom%continua(kr)%lambdamax)
            atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
            l0 = lambda(locate(lambda, atom%continua(kr)%lambda0))
 
           !in any problem of grid resolution etc or locate approximation.
           !We take Nred-1 to be sure than the tab_lambda_cont(Nred) <= lambda0.
           !Only if not dissolution.
           !We just need to avoind having cont_lambda(Nred)>lambda0, since the cross section is in (lambda/lambda0)**3
            if (l0 /= atom%continua(kr)%lambda0) then
               if (.not.ldissolve) then
                  if (lambda(atom%continua(kr)%Nr) > atom%continua(kr)%lambda0) then
                     call Warning("continuum Nred larger than lambda0 !")
                     write(*,*) " ++++ adjusting Nred", " lambda0=",atom%continua(kr)%lambda0, ' Nr=', atom%continua(kr)%Nr
                     do while (lambda(atom%continua(kr)%Nr) > atom%continua(kr)%lambda0)
                        atom%continua(kr)%Nr = atom%continua(kr)%Nr-1
                        atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
                        if (atom%continua(kr)%Nlambda == 1) call error('accomadting bound-free egdes')
                     enddo
                     write(*,*) "new val at Nred:", lambda(atom%continua(kr)%Nr), ' Nr=', atom%continua(kr)%Nr
                  !    write(*,*) atom%continua(kr)%lambda0
                  !    if (lambda(atom%continua(kr)%Nr-1) <= atom%continua(kr)%lambda0) then
                  !       write(*,*) " ++++ adjusting Nred", " lambda0=",atom%continua(kr)%lambda0
                  !       !To do while until <= lambda0
                  !       atom%continua(kr)%Nr = atom%continua(kr)%Nr-1
                  !       atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
                  !       write(*,*) "new val at Nred:", lambda(atom%continua(kr)%Nr)
                  !   endif
                     atom%continua(kr)%Nbc = locate(tab_lambda_cont, atom%continua(kr)%lambdamin)
                     atom%continua(kr)%Nrc = locate(tab_lambda_cont, lambda(atom%continua(kr)%Nr))
                     atom%continua(kr)%Nlambdac = atom%continua(kr)%Nrc - atom%continua(kr)%Nbc + 1
                  endif
               endif
            endif
            Nlambda_max_cont = max(Nlambda_max_cont,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1)
 
         enddo
 
         do kr=1,atom%Nline
            atom%lines(kr)%Nb = locate(lambda, atom%lines(kr)%lambdamin)
            atom%lines(kr)%Nr = locate(lambda, atom%lines(kr)%lambdamax)
            atom%lines(kr)%Nlambda = atom%lines(kr)%Nr - atom%lines(kr)%Nb + 1
            ! write(*,*) "line", kr, " lam0=",atom%lines(kr)%lambda0, atom%lines(kr)%lambdamin, atom%lines(kr)%lambdamax
            ! write(*,*) " -> bounds on the grid:", lambda(atom%lines(kr)%Nb), lambda(atom%lines(kr)%Nr)
            Nlambda_max_line = max(Nlambda_max_line, atom%lines(kr)%Nlambda)
            !second loop over all other lines of all other atoms to check for overlap !
            !init there are no overlaps
            ! *** natural overlaps are already taken into account into Nblue and Nred ***
            atom%lines(kr)%Nover_inf = atom%lines(kr)%Nb
            atom%lines(kr)%Nover_sup = atom%lines(kr)%Nr
            !Does not change Nlambda (?)
            if (check_for_overlap) then

               inner_atom_loop : do nb = 1, n_atoms
                  inner_line_loop : do krr=1, atoms(nb)%p%Nline

                     !Nover_inf/sup = Nb/r
                     if ( (nb==n).and.(kr==krr) ) cycle inner_line_loop

                     l0 = atom%lines(kr)%lambda0
                     l1 = atoms(nb)%p%lines(krr)%lambda0

                     write(*,*) l0, l1, abs(vmax_overlap)/c_light
                     write(*,*) abs(l1-l0)/l0, abs(vmax_overlap)/c_light

                     if ( abs(l1-l0)/l0 <= abs(vmax_overlap)/c_light ) then
                        dvmin = c_light * abs(l1-l0)/l0
                        ! atom%lines(kr)%dvmin = max(dvmin,atom%lines(kr)%dvmin)
                        atom%lines(kr)%Nover_sup = max(atom%lines(kr)%Nover_sup, locate(lambda, atom%lines(kr)%lambda0*(1.0 +  dvmin/c_light)))
                        atom%lines(kr)%Nover_inf = min(atom%lines(kr)%Nover_inf, locate(lambda, atom%lines(kr)%lambda0*(1.0 -  dvmin/c_light)))
                        write(*,*) "overlap of line", kr, atom%lines(kr)%lambda0, " of atom ", atom%ID, &
                           " with line", krr, atoms(nb)%p%lines(krr)%lambda0, " of atom ", atoms(nb)%p%ID, " N0 = ", &
                           locate(lambda, atoms(nb)%p%lines(krr)%lambda0), " dvmin=", dvmin *1d-3
                     endif


                  enddo inner_line_loop
               enddo inner_atom_loop 

               if (atom%lines(kr)%Nover_sup /= atom%lines(kr)%Nr .or. atom%lines(kr)%Nover_inf /= atom%lines(kr)%Nb) then
                  write(*,*) "Consistency of overlap must be checked!"
                  write(*,*) "line", kr, " is overlapping",  atom%lines(kr)%Nover_inf,  atom%lines(kr)%Nb, &
                                 atom%lines(kr)%Nover_sup,  atom%lines(kr)%Nr
                  !re calc Nlambda ??
               endif
            endif !check for overlap
         enddo
 
         atom => NULL()
      enddo
      write(*,*) "Number of max freq points for all lines resolution :", Nlambda_max_line
      write(*,*) "Number of max freq points for all cont resolution :", Nlambda_max_cont
      Nlambda_max_trans = max(Nlambda_max_line,Nlambda_max_cont)
      write(*,*) "Number of max freq points for all trans :", Nlambda_max_trans

      n_lambda_cont = size(tab_lambda_cont)
      ! !output grid in micron
      ! lambda = lambda * m_to_km
      !now lambda will be stored in micron for compatibility
      !but a lots of atomic variables use nm.
      write(*,'("Wavelength grid: "(1F12.4)" nm to",(1F12.4)" nm")') minval(lambda),maxval(lambda)

      return
   end subroutine make_wavelengths_nlte

!    subroutine make_wavelength_group()
      ! Divide the frequency interval in groups where lines overlap.
      ! if no overlap there is Nline group + Ngroup_cont.
      ! For the continua, groups are added between lines
      ! and similarly the cont groups include the overlap of continua.
      
      
      ! Store also the information on the transitions of each atom
      ! belonging to each group.
      ! Start with the same structure as
 
!      return
!    end subroutine make_wavelength_group

  subroutine make_wavelengths_raytracing(lambda)
   !Make a wavelength grid for images and flux that is different from
   !the non-LTE wavelength grid (the egdes of the lines take the velocity shift).
   !The grid is linearly sampled (in velocity) from -vmax_rt to +vmax_rt with n_speed_rt
   ! 1) if limage:
   !    keep only ray-traced lines of each atom and, the transitions (continua, other lines) with which they overlap 
   !    in the opacity (still the grid is defined only around those lines).
   !
   ! 2) if .not. limage:
   !    a) if lsed, use the lambda table in the parameter file to compute a sed.
   !    No image will be written. If line%Nlambda < 3 the line contirbution to the opacity is neglected for the flux.
   !    b) it .not. lsed, all transitions are kept taking into account the maximum velocity shift.
   !       This represent a high resolution spectrum on a grid similar to the non-LTE grid (except that it takes velocity shifts).
      real(kind=dp), intent(out), allocatable :: lambda(:)
      integer :: Ntrans, Ncont, Nlines, Nlambda_cont, Nlambda
      integer :: n, kr, Nlam, Nremoved, Nwaves
      type (AtomType), pointer :: atom
      integer :: alloc_status, lac, la, nb, nr, krr, Noverlap
      real(kind=dp), dimension(:), allocatable :: tmp_grid, tmp_grid2, all_l0, all_l1
      integer, parameter :: Ngroup_max = 1000
      real(kind=dp), dimension(Ngroup_max) :: group_blue_tmp, group_red_tmp, Nline_per_group_tmp
      integer, dimension(:), allocatable :: sorted_indexes
      real :: hv_loc !like hv
      real, parameter :: prec_vel = 0.0 !remove points is abs(hv_loc-hv)>prec_vel
   
      Ntrans = 0
      Nlines = 0 !ray-tracing
      Nlambda_cont = 0
      Ncont = 0
      do n=1, N_atoms
         atom => atoms(n)%p
         !km/s
         hv = 2.0 * atom%vmax_rt / (real(atom%n_speed_rt) - 1.0)
         do kr=1,atom%Ncont
            atom%continua(kr)%lcontrib = .false.
         enddo
         do kr=1,atom%Nline
            atom%lines(kr)%lcontrib = .false.
         enddo
         do kr=1, atom%nTrans_rayTracing
            krr = atom%ij_to_trans(atom%i_Trans_rayTracing(kr),atom%j_Trans_rayTracing(kr))
            atom%lines(krr)%lcontrib = .true.
            atom%lines(krr)%Nlambda = atom%n_speed_rt
            !recomputes line bound
            call compute_line_bound(atom%lines(krr),.true.)
            Nlines = Nlines + 1
         enddo

         Ntrans = Ntrans + atom%Ntr
         ! Nlines = Nlines + atom%Nline
         Ncont = Ncont + atom%Ncont
      enddo
      atom => null()


      if( Nlines==0 ) then
         call error("It seems that you want to compute an image but there is no lines!")
      endif

      !tab_lambda_cont = lambda here

   
      allocate(all_l0(Nlines), all_l1(Nlines), stat=alloc_status)

      Nlam = 0
      do n=1, N_atoms
         atom => atoms(n)%p
 
         do kr=1,atom%Nline
            if (.not.atom%lines(kr)%lcontrib) cycle
            Nlam = Nlam + 1
            all_l0(Nlam) = atom%lines(kr)%lambdamin
            all_l1(Nlam) = atom%lines(kr)%lambdamax
         enddo
 
      enddo
      atom => null()
 
      allocate(sorted_indexes(Nlam),stat=alloc_status)
      if (alloc_status > 0) then
         call error("Allocation error sorted_indexes(Nlam)")
      endif
      sorted_indexes = index_bubble_sort(all_l0)
      !I sort the lambdamin, and lambdamax array is ordered so that
      !the index in all_l0 and in all_l1 correspond to the same pair of 
      !lambda min/max
      all_l0(:) = all_l0(sorted_indexes)
      all_l1(:) = all_l1(sorted_indexes)
      deallocate(sorted_indexes)

      group_blue_tmp(:) = -1.0
      group_red_tmp(:) = -1.0
 
      N_groups = 1
      group_blue_tmp(N_groups) = all_l0(1)
      group_red_tmp(N_groups) = all_l1(1)

      Nline_per_group_tmp(:) = 0
      Nline_per_group_tmp(N_groups) = 1
      do Nlam = 2, Nlines
 
         !Is the line overlapping the previous line ?
 
         !Yes, add it to the same group
         if (((all_l0(Nlam) >= group_blue_tmp(N_groups)).and.&
            (all_l0(Nlam) <= group_red_tmp(N_groups))).or.&
            ((all_l1(Nlam) >= group_blue_tmp(N_groups)).and.&
            (all_l1(Nlam) <= group_red_tmp(N_groups)))) then
 
            group_blue_tmp(N_groups) = min(all_l0(Nlam), group_blue_tmp(N_groups))
            group_red_tmp(N_groups) = max(all_l1(Nlam), group_red_tmp(N_groups))
 
            Nline_per_group_tmp(N_groups) = Nline_per_group_tmp(N_groups) + 1
 
         !no, create a new group, starting with this line at first element
         else
            N_groups = N_groups + 1
            if (N_groups > Ngroup_max) then
               call error(" Error, Ngroup > Ngroup_max")
            endif
            group_blue_tmp(N_groups) = all_l0(Nlam)
            group_red_tmp(N_groups) = all_l1(Nlam)
            Nline_per_group_tmp(N_groups) = 1
         endif
 
      enddo !end loop to create group
      allocate(Nlambda_per_group(N_groups), group_red(N_groups), &
      group_blue(N_groups), Nline_per_group(N_groups), stat=alloc_status)
      if (alloc_status > 0) then
         call error("Allocation error groups")
      endif
      group_blue(:) = group_blue_tmp(1:N_groups)
      group_red(:) = group_red_tmp(1:N_groups)
      Nline_per_group(:) = Nline_per_group_tmp(1:N_groups)


      allocate(tmp_grid(1000000),stat=alloc_status)
      if (alloc_status > 0) call error("Allocation error tmp_grid (line)!")
      tmp_grid = -1.0
 
      !groups are not really used here
      lac = 1
      do n=1, N_atoms
         atom => atoms(n)%p
 
         do kr=1,atom%Nline
            if (.not.atom%lines(kr)%lcontrib) cycle

            tmp_grid(lac:atom%lines(kr)%Nlambda+lac-1) = line_lambda_grid_dv(atom%lines(kr),atom%lines(kr)%Nlambda)
            lac = lac + atom%lines(kr)%Nlambda

         enddo

      enddo
      atom => null()
      tmp_grid = pack(tmp_grid,tmp_grid > 0)
      Nlambda = size(tmp_grid)
      allocate(sorted_indexes(Nlambda),stat=alloc_status)
      if (alloc_status > 0) call error ("Allocation error sorted_indexes (line)")
      sorted_indexes = index_bubble_sort(tmp_Grid)
      tmp_grid(:) = tmp_grid(sorted_indexes)
      deallocate(sorted_indexes)
      !count the number of wavelengths per groups ? 
      !remove duplicates or points below resolution.
      allocate(tmp_grid2(Nlambda), stat=alloc_status)
      if (alloc_status > 0) call error ("Allocation error tmp_grid2")
      tmp_grid2(2:Nlambda) = 0.0
      tmp_grid2(1) = tmp_grid(1)
      Nremoved = 0

      !here we will remove points that are also below hv so as to preserve the grid resolutin.
      !however we keep points with hv_loc above hv as it means we are comparing lambda of different
      !line groups.
      !TO DO: check that when hv_loc >> hv it is indeed  because we change group.
      if (hv>0.0) then
         do la = 2, Nlambda
            hv_loc = ceiling ( real (1d-3 * c_light * (tmp_grid(la)-tmp_grid(la-1)) / tmp_grid(la) ) )
            if (abs(hv_loc - hv) >= prec_vel) then !hv or more.
               tmp_grid2(la) = tmp_grid(la)
            else!below grid resolution, remove.
               Nremoved = Nremoved + 1
            endif
         enddo
      else
         do la = 2, Nlambda
            if (tmp_grid(la) > tmp_grid(la-1)) then
               tmp_grid2(la) = tmp_grid(la)
            else
               Nremoved = Nremoved + 1
            endif
         enddo
      endif
      if (Nremoved > 0) then
         write(*,*) " ->", Nremoved, " duplicate frequencies in lines"
         deallocate(tmp_grid)
         ! allocate(tmp_grid(Nlambda-Nremoved), stat=alloc_status)
         tmp_grid(:) = Pack(tmp_grid2, tmp_grid2 > 0)
      endif
      deallocate(tmp_grid2)
      Nlambda = size(tmp_grid)!Nlambda - Nremoved
      Nwaves = Nlambda

      allocate(lambda(Nwaves),stat=alloc_status)
      tmp_grid = tmp_grid(index_bubble_sort(tmp_grid))
      lambda(:) = -99.0 !check
      lambda = pack(tmp_grid, mask=tmp_grid > 0)
 
      do lac=2, Nwaves
         if (lambda(lac) <= lambda(lac-1)) then
            write(*,*) lac, "lambda = ", lambda(lac), lambda(lac-1), minval(lambda)
            call error("Sorted problem")
         endif
      enddo
      deallocate(tmp_grid)

      allocate(tab_lambda_cont(size(lambda)))
      tab_lambda_cont = lambda
      n_lambda_cont = size(tab_lambda_cont)
      write(*,*) Nlambda, " line wavelengths"


      !allocate indexes on the grid
      Nlambda_max_line = 0
      Nlambda_max_cont = 0
      Nlambda_max_trans = 0
      do n=1,n_atoms
         atom => Atoms(n)%p
         cont_loop : do kr=1,atom%Ncont


            atom%continua(kr)%Nb = 2*size(lambda)
            atom%continua(kr)%Nr = 0
            Noverlap = 0       
            loop_group_a : do lac=1,n_groups
               if ( (atom%continua(kr)%lambdamax >= group_blue(lac)).and.(atom%continua(kr)%lambdamin <= group_red(lac)) ) then
                  atom%continua(kr)%lcontrib = .true.
                  Noverlap = Noverlap + 1
    
                  nb = locate(lambda,group_blue(lac))
                  nr = locate(lambda,group_red(lac))
                  atom%continua(kr)%Nb = min(atom%continua(kr)%Nb,nb)
                  atom%continua(kr)%Nr = max(atom%continua(kr)%Nr,nr)
        
                  atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
               endif
 
            enddo loop_group_a
            atom%continua(kr)%Nbc = atom%continua(kr)%Nb
            atom%continua(kr)%Nrc = atom%continua(kr)%Nr
            atom%continua(kr)%Nlambdac = atom%continua(kr)%Nlambda
            if (Noverlap > 0) then
               if (Noverlap==1) then
                  write(*,*) " *** cont transition ", kr, " of atom ", atom%ID, " overlaps with ",Noverlap," line" 
               else
                  write(*,*) " *** cont transition ", kr, " of atom ", atom%ID, " overlaps with ",Noverlap," lines" 
               endif
               write(*,*) "    -> bounds ", lambda(atom%continua(kr)%Nb), lambda(atom%continua(kr)%Nr)
            endif

            Nlambda_max_cont = max(Nlambda_max_cont,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1)
         enddo cont_loop


         line_loop : do kr=1,atom%Nline
            atom%lines(kr)%Nb = locate(lambda, atom%lines(kr)%lambdamin)
            atom%lines(kr)%Nr = locate(lambda, atom%lines(kr)%lambdamax)
            atom%lines(kr)%Nlambda = atom%lines(kr)%Nr - atom%lines(kr)%Nb + 1


            !line is not a raytraced line ? lcontrib is .false. (since begining)
            if (.not.atom%lines(kr)%lcontrib) then
               !just check that this line is not overlapping, that is nblue and nred falls onto the grid.
               !The grid has be generated to include only lines for the raytracing, so we test if the other lines partially fall
               !in these regions.


               !testing if a line overlap at least one of the raytraced line
               !the line bounds are unchanged (look at the continuum)
               loop_group_b : do lac=1,n_groups
               !the shift is added to check if the line will overlap a group
                  if ( ((atom%lines(kr)%lambdamax < group_red(lac)).and.&
                        (atom%lines(kr)%lambdamax > group_blue(lac))).or.&
                        ((atom%lines(kr)%lambdamin > group_blue(lac)).and.&
                           (atom%lines(kr)%lambdamin < group_red(lac))) ) then
                     atom%lines(kr)%lcontrib = .true.
                     write(*,*) " *** line transition ", kr, " of atom ", &
                        atom%ID, " overlaps with another line" 
                     write(*,*) "    -> in group ", lac, group_blue(lac), group_red(lac)
                     exit loop_group_b
                  endif
               enddo loop_group_b
               !find nothing probably out otherwise lcontrib is true (even if it was false first)
            endif

            Nlambda_max_line = max(Nlambda_max_line, atom%lines(kr)%Nlambda)
         enddo line_loop

         atom => NULL()
      enddo
      write(*,*) "Number of max freq points for all lines at this resolution :", Nlambda_max_line
      write(*,*) "Number of max freq points for all cont at this resolution :", Nlambda_max_cont
      Nlambda_max_trans = max(Nlambda_max_line,Nlambda_max_cont)
      write(*,*) "Number of max freq points for all trans at this resolution :", Nlambda_max_trans

      return
   end subroutine make_wavelengths_raytracing

   subroutine make_wavelengths_flux(lambda,lfrom_file)
      real(kind=dp), intent(inout), allocatable :: lambda(:) !can be the one defined in the parameter file.
      logical, intent(in) :: lfrom_file
      integer :: Ntrans, Ncont, Nlines, Nlambda_cont, Nlambda
      integer :: n, kr, Nlam, Nmore_cont_freq, Nremoved, Nwaves, check_new_freq
      type (AtomType), pointer :: atom
      integer :: alloc_status, lac, la, nb, krr
      real(kind=dp), dimension(:), allocatable :: tmp_grid, tmp_grid2, all_l0, all_l1
      integer, parameter :: Ngroup_max = 1000
      real(kind=dp), dimension(Ngroup_max) :: group_blue_tmp, group_red_tmp, Nline_per_group_tmp
      integer, dimension(:), allocatable :: sorted_indexes
      real(kind=dp) :: max_cont, l0, l1, lr
      real :: hv_loc !like hv
      real, parameter :: prec_vel = 0.0 !remove points is abs(hv_loc-hv)>prec_vel

      if (.not.lfrom_file) then
 
         Ntrans = 0
         Nlines = 0 !ray-tracing
         Nlambda_cont = 0
         Ncont = 0
         do n=1, N_atoms
            atom => atoms(n)%p
         !km/s
            hv = 2.0 * atom%vmax_rt / (real(atom%n_speed_rt) - 1.0)
            do kr=1, atom%nline
               atom%lines(kr)%Nlambda = atom%n_speed_rt
               call compute_line_bound(atom%lines(kr),.true.)
            enddo

            do kr=1,atom%Ncont
               if (atom%continua(kr)%hydrogenic) then
                  atom%continua(kr)%Nlambda = Nlambda_cont_log
                  if (atom%continua(kr)%lambdamax > atom%continua(kr)%lambda0) then
                     atom%continua(kr)%Nlambda = atom%continua(kr)%Nlambda + Nlambda_cont_log + 1
                  endif  
               else
                  !needs a reset
                  atom%continua(kr)%Nlambda = size(atom%continua(kr)%lambda_file)             
               endif
               Nlambda_cont = Nlambda_cont + atom%continua(kr)%Nlambda
            enddo

            Ntrans = Ntrans + atom%Ntr
            Nlines = Nlines + atom%Nline
            Ncont = Ncont + atom%Ncont
         enddo
         atom => null()

         ! ********************** Pure continuum RT ********************** !
         !Even if continuum RT is removed, that part is important
         allocate(tab_lambda_cont(Nlambda_cont), stat=alloc_status)
         if (alloc_status>0) then
            call error("Allocation error cont_waves")
         endif

         lac = 1
         do n=1, N_atoms
            atom => atoms(n)%p
            do kr=1,atom%Ncont
               if (atom%continua(kr)%hydrogenic) then
                  tab_lambda_cont(lac:(lac-1)+atom%continua(kr)%Nlambda) = subgrid_cont(atom%continua(kr))
               else 
                  tab_lambda_cont(lac:lac-1+atom%continua(kr)%Nlambda) = atom%continua(kr)%lambda_file(:)
               endif
               lac = lac + atom%continua(kr)%Nlambda
            enddo
            atom => null()
         enddo
         ! write(*,*) "lac end=", lac-1
 
         !sort continuum frequencies
         Nmore_cont_freq = 0.0
         allocate(sorted_indexes(Nlambda_cont),stat=alloc_status)
         if (alloc_status > 0) call error ("Allocation error sorted_indexes (cont)")
         sorted_indexes = index_bubble_sort(tab_lambda_cont)
         tab_lambda_cont(:) = tab_lambda_cont(sorted_indexes)
         deallocate(sorted_indexes)
         !remove duplicates
         allocate(tmp_grid(Nlambda_cont), stat=alloc_status)
         if (alloc_status > 0) call error ("Allocation error tmp_grid (cont)")
         tmp_grid(2:Nlambda_cont) = 0.0
         tmp_grid(1) = tab_lambda_cont(1)
         Nremoved = 0
         do la = 2, Nlambda_cont
            if (tab_lambda_cont(la) > tab_lambda_cont(la-1)) then
               tmp_grid(la) = tab_lambda_cont(la)
          else
               Nremoved = Nremoved + 1
            endif
         enddo
 
         if (Nremoved > 0) then
            write(*,*) " ->", Nremoved, " duplicate frequencies"
            deallocate(tab_lambda_cont)
            allocate(tab_lambda_cont(Nlambda_cont-Nremoved), stat=alloc_status)
            tab_lambda_cont(:) = Pack(tmp_grid, tmp_grid > 0)
         endif
         deallocate(tmp_grid)
         max_cont = maxval(tab_lambda_cont)
         Nlambda_cont = Nlambda_cont - Nremoved
         ! ********************** ***************** ********************** !

         ! ********************** cont + line    RT ********************** !
         if (Nlines > 0) then !(.true.) then!lines are present in some atoms
            allocate(all_l0(Nlines), all_l1(Nlines), stat=alloc_status)

            Nlam = 0
            do n=1, N_atoms
               atom => atoms(n)%p
 
               do kr=1,atom%Nline
                  Nlam = Nlam + 1
                  all_l0(Nlam) = atom%lines(kr)%lambdamin
                  all_l1(Nlam) = atom%lines(kr)%lambdamax
               enddo
 
            enddo
            atom => null()
 
            allocate(sorted_indexes(Nlam),stat=alloc_status)
            if (alloc_status > 0) then
               call error("Allocation error sorted_indexes(Nlam)")
            endif
            sorted_indexes = index_bubble_sort(all_l0)
            !I sort the lambdamin, and lambdamax array is ordered so that
            !the index in all_l0 and in all_l1 correspond to the same pair of 
            !lambda min/max
            all_l0(:) = all_l0(sorted_indexes)
            all_l1(:) = all_l1(sorted_indexes)
            deallocate(sorted_indexes)

            group_blue_tmp(:) = -1.0
            group_red_tmp(:) = -1.0
 
            N_groups = 1
            group_blue_tmp(N_groups) = all_l0(1)
            group_red_tmp(N_groups) = all_l1(1)
            !Find group of lines, and store for each group the lambda_blue and lambda_red of each group
            !if a line overlaps with the previous line, add it to the same group and check the next line.
            !Stop counting lines in a group if the next line does not overlap with the previous line. In
            !the latter case, create a new group and start again.
            ! Note: the first and last lines of a group may not overlap.
            Nline_per_group_tmp(:) = 0
            Nline_per_group_tmp(N_groups) = 1
            do Nlam = 2, Nlines
 
            !Is the line overlapping the previous line ?
 
            !Yes, add it to the same group
               if (((all_l0(Nlam) >= group_blue_tmp(N_groups)).and.&
                  (all_l0(Nlam) <= group_red_tmp(N_groups))).or.&
                  ((all_l1(Nlam) >= group_blue_tmp(N_groups)).and.&
                  (all_l1(Nlam) <= group_red_tmp(N_groups)))) then
 
                  group_blue_tmp(N_groups) = min(all_l0(Nlam), group_blue_tmp(N_groups))
                  group_red_tmp(N_groups) = max(all_l1(Nlam), group_red_tmp(N_groups))
 
                  Nline_per_group_tmp(N_groups) = Nline_per_group_tmp(N_groups) + 1
 
            !no, create a new group, starting with this line at first element
               else
                  N_groups = N_groups + 1
                  if (N_groups > Ngroup_max) then
                     call error(" Error, Ngroup > Ngroup_max")
                  endif
                  group_blue_tmp(N_groups) = all_l0(Nlam)
                  group_red_tmp(N_groups) = all_l1(Nlam)
                  Nline_per_group_tmp(N_groups) = 1
               endif
 
            enddo !end loop to create group
            allocate(Nlambda_per_group(N_groups), group_red(N_groups), &
            group_blue(N_groups), Nline_per_group(N_groups), stat=alloc_status)
            if (alloc_status > 0) then
               call error("Allocation error groups")
            endif
            group_blue(:) = group_blue_tmp(1:N_groups)
            group_red(:) = group_red_tmp(1:N_groups)
            Nline_per_group(:) = Nline_per_group_tmp(1:N_groups)

            !some statistics
            write(*,*) " Found ", N_groups, " groups of lines for ", Nlines, " lines in total."
            write(*,*) "   -> ", sum(Nline_per_group) - N_groups, " overlapping regions"! for", sum(Nline_per_group), " lines"
            ! do Nlam=1, N_groups
            !    write(*,*) " group ", Nlam, " lamin = ", group_blue(Nlam), ' lamax = ', group_red(Nlam)
            ! enddo

            allocate(tmp_grid(1000000),stat=alloc_status)
            if (alloc_status > 0) call error("Allocation error tmp_grid (line)!")
            tmp_grid = -1.0
 
            !groups are not really used here
            lac = 1
            do n=1, N_atoms
               atom => atoms(n)%p
 
               do kr=1,atom%Nline

                  tmp_grid(lac:atom%lines(kr)%Nlambda+lac-1) = line_lambda_grid_dv(atom%lines(kr),atom%lines(kr)%Nlambda)
                  lac = lac + atom%lines(kr)%Nlambda

               enddo

            enddo
            atom => null()
            tmp_grid = pack(tmp_grid,tmp_grid > 0)
            Nlambda = size(tmp_grid)
            allocate(sorted_indexes(Nlambda),stat=alloc_status)
            if (alloc_status > 0) call error ("Allocation error sorted_indexes (line)")
            sorted_indexes = index_bubble_sort(tmp_Grid)
            tmp_grid(:) = tmp_grid(sorted_indexes)
            deallocate(sorted_indexes)
            !count the number of wavelengths per groups ? 
            !remove duplicates or points below resolution.
            allocate(tmp_grid2(Nlambda), stat=alloc_status)
            if (alloc_status > 0) call error ("Allocation error tmp_grid2")
            tmp_grid2(2:Nlambda) = 0.0
            tmp_grid2(1) = tmp_grid(1)
            Nremoved = 0

            !here we will remove points that are also below hv so as to preserve the grid resolutin.
            !however we keep points with hv_loc above hv as it means we are comparing lambda of different
            !line groups.
            !TO DO: check that when hv_loc >> hv it is indeed  because we change group.
            if (hv>0.0) then
               do la = 2, Nlambda
                  hv_loc = ceiling ( real (1d-3 * c_light * (tmp_grid(la)-tmp_grid(la-1)) / tmp_grid(la) ) )
                  if (abs(hv_loc - hv) >= prec_vel) then !hv or more.
                     tmp_grid2(la) = tmp_grid(la)
                  else!below grid resolution, remove.
                     Nremoved = Nremoved + 1
                  endif
               enddo
            else
               do la = 2, Nlambda
                  if (tmp_grid(la) > tmp_grid(la-1)) then
                     tmp_grid2(la) = tmp_grid(la)
                  else
                     Nremoved = Nremoved + 1
                  endif
               enddo
            endif
            if (Nremoved > 0) then
               write(*,*) " ->", Nremoved, " duplicate frequencies in lines"
               deallocate(tmp_grid)
               ! allocate(tmp_grid(Nlambda-Nremoved), stat=alloc_status)
               tmp_grid(:) = Pack(tmp_grid2, tmp_grid2 > 0)
            endif
            deallocate(tmp_grid2)
            Nlambda = size(tmp_grid)!Nlambda - Nremoved
            ! do n=1,N_groups
            !    write(*,*) n, size(pack(tmp_grid,&
            !    (tmp_grid >= group_blue(n)).and.tmp_grid <= group_red(n)))
            !    Nlambda_per_group(n) = size(pack(tmp_grid,&
            !       (tmp_grid >= group_blue(n)).and.tmp_grid <= group_red(n)))
            ! enddo

            ! ----------- small cont grid start ------------- !
            !-> This is only useful for the pure-continuum RT !!!!
            !and should disappear at some point.
            !-> Add continuum points beyond last "bound-free" continuum
            !In case they are lines beyond the last continuum I add at least3 points per line for the continuum in this region
            !->cont end		this is heavy for nothing but should work !
            !finalise continuum here by reckoning how much freq we need
            !make the new points go farther than line_waves for interpolation.
            Nmore_cont_freq = 0
            do n=1, N_groups
               l0 = group_blue(n)
               l1 = group_red(n)
               if (l0 > max_cont) then
                  Nmore_cont_freq = Nmore_cont_freq + 1
               endif
               if (l1 > max_cont) then
                  Nmore_cont_freq = Nmore_cont_freq + 1
               endif
               if (0.5*(l0+l1) > max_cont) then
                  Nmore_cont_freq = Nmore_cont_freq + 1
               endif
 
            enddo
            check_new_freq = Nmore_cont_freq
            if (Nmore_cont_freq > 0) then
               write(*,*) "Adding new wavelength points for lines beyond continuum max!"
               write(*,*) "  -> Adding ", Nmore_cont_freq," points"
               write(*,*) "max cont, max line", max_cont, maxval(tmp_grid)
               allocate(tmp_grid2(Nlambda_cont))
               tmp_grid2 = tab_lambda_cont
               deallocate(tab_lambda_cont)
               allocate(tab_lambda_cont(Nlambda_cont + Nmore_cont_freq))
               tab_lambda_cont(1:Nlambda_cont) = tmp_grid2(:)
               deallocate(tmp_grid2)
               allocate(tmp_grid2(Nmore_cont_freq))
               tmp_grid2(:) = 0.0_dp
 
               Nmore_cont_freq = 0
               do n=1, N_groups
                  l0 = group_blue(n)
                  l1 = group_red(n)
                  if (l0 > max_cont) then
                     Nmore_cont_freq = Nmore_cont_freq + 1
                     tmp_grid2(Nmore_cont_freq) = l0
                  endif
                  if (0.5*(l0+l1) > max_cont) then
                     Nmore_cont_freq = Nmore_cont_freq + 1
                     tmp_grid2(Nmore_cont_freq) = 0.5 * (l0+l1)
                  endif
                  if (l1 > max_cont) then
                     Nmore_cont_freq = Nmore_cont_freq + 1
                     tmp_grid2(Nmore_cont_freq) = l1
                  endif
               enddo
               if (Nmore_cont_freq /= check_new_freq) then
                  call Warning("There are probably some frequency missing!")
                  write(*,*) "Nmore_freq: ",check_new_freq," Nfreq_added: ", Nmore_cont_freq
               endif
 
               allocate(sorted_indexes(Nmore_cont_freq))
               sorted_indexes(:) = index_bubble_sort(tmp_grid2)
               tmp_grid2(:) = tmp_grid2(sorted_indexes)
               tab_lambda_cont(Nlambda_cont+1:Nlambda_cont + Nmore_cont_freq) = tmp_grid2(:)
               Nlambda_cont = Nlambda_cont + Nmore_cont_freq
               deallocate(tmp_grid2, sorted_indexes)
            endif
 
            if (size(tab_lambda_cont) /= Nlambda_cont) then
               write(*,*) " Something went wrong with Nlambda cont"
               stop
            endif

            ! ----------- small cont grid end ------------- !

            !initiate with lines
            allocate(tmp_grid2(nlambda+Nlambda_cont), stat=alloc_status)
            if (alloc_status > 0) call error ("Allocation error tmp2_grid (final)")
            tmp_grid2(:) = -99
            tmp_grid2(1:Nlambda) = tmp_grid(:)
 
            Nwaves = Nlambda
            !add continuum wavlengths only outside line groups
            !First values below or beyond first and last groups
            la = 0
            do lac=Nlambda+1, Nlambda+Nlambda_cont
               if ((tab_lambda_cont(lac-Nlambda) < group_blue(1)) .or. (tab_lambda_cont(lac-Nlambda) > group_red(N_groups))) then
                   tmp_grid2(lac) = tab_lambda_cont(lac-Nlambda)
                   Nwaves = Nwaves + 1
                   if (tab_lambda_cont(lac-Nlambda) < group_blue(1)) la = lac
               endif
            enddo
 
            !now values between groups
            do lac=la+1, Nlambda_cont+Nlambda
               group_loop : do n=2, N_groups
                  if ((tab_lambda_cont(lac-Nlambda) > group_red(n-1)).and.(tab_lambda_cont(lac-Nlambda) < group_blue(n))) then
                     Nwaves = Nwaves + 1
                     tmp_grid2(lac) = tab_lambda_cont(lac-nlambda)
                 !else
                 ! be smart and cycle to accelerate
                  endif
               enddo group_loop
            enddo
            ! write(*,*) " Check Nwaves:", Nwaves,  size(pack(tmp_grid2, tmp_grid2 > 0))
            deallocate(tmp_grid)
 
            !continuum frequencies are sorted and so are the line frequencies
            !but they are added at the end, so sorted is needed, but I can improve the previous
            !loop to fill the tmp_frid in the ascending order of wavelengths
            allocate(lambda(Nwaves),stat=alloc_status)
            tmp_grid2 = tmp_grid2(index_bubble_sort(tmp_grid2))
            lambda(:) = -99.0 !check
            lambda = pack(tmp_grid2, mask=tmp_grid2 > 0)
 
            do lac=2, Nwaves
               if (lambda(lac) <= lambda(lac-1)) then
                  write(*,*) lac, "lambda = ", lambda(lac), lambda(lac-1), minval(lambda)
                  call error("Sorted problem")
               endif
            enddo
            deallocate(tmp_grid2)
 
         else !pure cont
            Nwaves = Nlambda_cont
            allocate(lambda(Nwaves),stat=alloc_status)
            if (alloc_status>0) call error("allocation error lambda, in pure cont!")
            lambda = tab_lambda_cont
            deallocate(tab_lambda_cont)
            Nlambda = 0
         endif !there is lines
      else
         !lambda is in micron and is read from the parameter file
         Nwaves = size(lambda)
         Nlambda = Nwaves
         Nlambda_cont = Nwaves
         allocate(tab_lambda_nm(size(lambda)), tab_lambda_cont(size(lambda)))
         tab_lambda_nm = lambda * km_to_m !micron to nm
         tab_lambda_cont = tab_lambda_nm
         !we just need to check the bounds of the transitions on this grid
         !in particular bound-free egdes and lines.
         do n=1,n_atoms
            do kr=1,atoms(n)%p%Nline
               call compute_line_bound(atoms(n)%p%lines(kr),.true.)
            enddo
         enddo
      endif

      write(*,*) Nwaves, " unique wavelengths" !they are no eliminated lines
      if (.not.lfrom_file) then
         write(*,*) Nlambda, " line wavelengths"
         write(*,*) Nlambda_cont, " continuum wavelengths"
      endif

      !Now indexes of each transition on the lambda grid
      Nlambda_max_line = 0
      Nlambda_max_cont = 0
      Nlambda_max_trans = 0
      do n=1,N_atoms
         atom => Atoms(n)%p
         do kr=1,atom%Ncont

            atom%continua(kr)%Nbc = locate(tab_lambda_cont, atom%continua(kr)%lambdamin)
            atom%continua(kr)%Nrc = locate(tab_lambda_cont, atom%continua(kr)%lambdamax)
            atom%continua(kr)%Nlambdac = atom%continua(kr)%Nrc - atom%continua(kr)%Nbc + 1

            if (lfrom_file) then
               atom%continua(kr)%Nb = locate(tab_lambda_nm, atom%continua(kr)%lambdamin)
               atom%continua(kr)%Nr = locate(tab_lambda_nm, atom%continua(kr)%lambdamax)
               atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
               l0 = tab_lambda_nm(locate(tab_lambda_nm, atom%continua(kr)%lambda0))
               lr = tab_lambda_nm(atom%continua(kr)%Nr)
            else
               atom%continua(kr)%Nb = locate(lambda, atom%continua(kr)%lambdamin)
               atom%continua(kr)%Nr = locate(lambda, atom%continua(kr)%lambdamax)
               atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
               l0 = lambda(locate(lambda, atom%continua(kr)%lambda0))
               lr = lambda(atom%continua(kr)%Nr)
            endif
            if (l0 /= atom%continua(kr)%lambda0) then
               if (.not.ldissolve) then
                  if (lr > atom%continua(kr)%lambda0) then
                     call Warning("continuum Nred larger than lambda0 !")
                     write(*,*) " ++++ adjusting Nred", " lambda0=",atom%continua(kr)%lambda0, ' Nr=', atom%continua(kr)%Nr
                     if (lfrom_file) then
                        do while (tab_lambda_nm(atom%continua(kr)%Nr) > atom%continua(kr)%lambda0)
                           write(*,*) tab_lambda_nm(atom%continua(kr)%Nr), atom%continua(kr)%lambda0
                           atom%continua(kr)%Nr = atom%continua(kr)%Nr-1
                           atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
                           if (atom%continua(kr)%Nlambda == 1) call error('accomadting bound-free egdes')
                        enddo
                        atom%continua(kr)%Nbc = atom%continua(kr)%Nb
                        atom%continua(kr)%Nrc = atom%continua(kr)%Nr
                        atom%continua(kr)%Nlambdac = atom%continua(kr)%Nlambda
                        write(*,*) "new val at Nred:", tab_lambda_nm(atom%continua(kr)%Nr), ' Nr=', atom%continua(kr)%Nr
                     else
                        do while (lambda(atom%continua(kr)%Nr) > atom%continua(kr)%lambda0)
                           atom%continua(kr)%Nr = atom%continua(kr)%Nr-1
                           atom%continua(kr)%Nlambda = atom%continua(kr)%Nr - atom%continua(kr)%Nb + 1
                           if (atom%continua(kr)%Nlambda == 1) call error('accomadting bound-free egdes')
                        enddo
                        atom%continua(kr)%Nbc = locate(tab_lambda_cont, atom%continua(kr)%lambdamin)
                        atom%continua(kr)%Nrc = locate(tab_lambda_cont, lambda(atom%continua(kr)%Nr))
                        atom%continua(kr)%Nlambdac = atom%continua(kr)%Nrc - atom%continua(kr)%Nbc + 1
                        write(*,*) "new val at Nred:", lambda(atom%continua(kr)%Nr), ' Nr=', atom%continua(kr)%Nr
                     endif
                  endif
               endif
            endif
            Nlambda_max_cont = max(Nlambda_max_cont,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1)
 
         enddo
 
         do kr=1,atom%Nline
            atom%lines(kr)%Nb = locate(lambda, atom%lines(kr)%lambdamin)
            atom%lines(kr)%Nr = locate(lambda, atom%lines(kr)%lambdamax)
            atom%lines(kr)%Nlambda = atom%lines(kr)%Nr - atom%lines(kr)%Nb + 1
            ! write(*,*) "line", kr, " lam0=",atom%lines(kr)%lambda0, atom%lines(kr)%lambdamin, atom%lines(kr)%lambdamax
            ! write(*,*) " -> bounds on the grid:", lambda(atom%lines(kr)%Nb), lambda(atom%lines(kr)%Nr)
            Nlambda_max_line = max(Nlambda_max_line, atom%lines(kr)%Nlambda)
            !second loop over all other lines of all other atoms to check for overlap !
            !init there are no overlaps
            ! *** natural overlaps are already taken into account into Nblue and Nred ***
            !Does not change Nlambda (?)
            atom%lines(kr)%Nover_inf = atom%lines(kr)%Nb
            atom%lines(kr)%Nover_sup = atom%lines(kr)%Nr
            if (lfrom_file) then
               !check edges
               if (atom%lines(kr)%Nlambda < 3) atom%lines(kr)%lcontrib = .false.
            endif
         enddo
 
         atom => NULL()
      enddo
      write(*,*) "Number of max freq points for all lines resolution :", Nlambda_max_line
      write(*,*) "Number of max freq points for all cont resolution :", Nlambda_max_cont
      Nlambda_max_trans = max(Nlambda_max_line,Nlambda_max_cont)
      write(*,*) "Number of max freq points for all trans :", Nlambda_max_trans

      n_lambda_cont = size(tab_lambda_cont)
      ! !output grid in micron
      ! lambda = lambda * m_to_km
      !now lambda will be stored in micron for compatibility
      !but a lots of atomic variables use nm.
      if (lfrom_file) then
         write(*,'("Wavelength grid: "(1F12.4)" nm to ",(1F12.4)" nm")') minval(tab_lambda_nm),maxval(tab_lambda_nm)
         write(*,'(" -> "(1F12.4)" mum to ",(1F12.4)" mum")') minval(lambda),maxval(lambda)
      else
         write(*,'("Wavelength grid: "(1F12.4)" nm to",(1F12.4)" nm")') minval(lambda),maxval(lambda)
      endif

      return
   end subroutine make_wavelengths_flux

 end module wavelengths_gas
 
 
 
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

    ! subroutine make_sub_wavelength_grid_cont_lin(cont, lambdamin, lambdamax)
   ! ! ----------------------------------------------------------------- !
   ! ! Make an individual wavelength grid for the AtomicContinuum cont.
   ! !  -> cross-section is extrapolated beyond edge to be use with
   ! ! level's dissolution, if lambdamax > lambda0
   ! !
   ! ! Allocate cont%lambda.
   ! ! ----------------------------------------------------------------- !
   !   type (AtomicContinuum), intent(inout) :: cont
   !   real(kind=dp), intent(in) :: lambdamin, lambdamax
   !   real(kind=dp) :: resol
   !   integer :: la, N1, N2, Nmid
   !   real(kind=dp) :: l0, l1, dl
 
   !   l1 = lambdamax
   !   l0 = lambdamin
 
   !   if (lambdamax > cont%lambda0) then
   !      N2 = Nlambda_cont + 1
   !      N1 = Nlambda_cont
   !      Nmid = N1 + N2
   !   else
   !      N2 = 0
   !      N1 = Nlambda_cont
   !      Nmid = N1
   !   endif
   !   cont%Nlambda = N1 + N2
   !   allocate(cont%lambda(cont%Nlambda))
 
   !   cont%lambda(1:N1) = span_dp(l0, cont%lambda0,N1,-1)!from lam0 to 1
   !   !still result in lamin to lamax
 
   !   ! N1+1 points are the same
   !   if (N2>0) then
   !      dl = abs(cont%lambda(N1) - cont%lambda(N1-1)) !2-1
   !      cont%lambda(N1+1:N2+N1) = span_dp(cont%lambda0+dl,l1,N2,1)
   !   endif
 
   !   do la=1,cont%Nlambda
   !      if (cont%lambda(la) < 0) then
   !         write(*,*) "Error, lambda negative"
   !         write(*,*) "cont lin"
   !         stop
   !      endif
   !   end do

   !   return
   ! end subroutine make_sub_wavelength_grid_cont_lin
 

   ! subroutine make_sub_wavelength_grid_cont_linlog(cont, lambdamin, lambdamax)
   ! ! ----------------------------------------------------------------- !
   ! ! Make an individual wavelength grid for the AtomicContinuum cont.
   ! !  -> cross-section is extrapolated beyond edge to be use with
   ! ! level's dissolution, if lambdamax > lambda0
   ! !
   ! ! Allocate cont%lambda.
   ! ! linear from 0 to lambda0 then log from lambda0 to lambdamax
   ! ! ----------------------------------------------------------------- !
   !   type (AtomicContinuum), intent(inout) :: cont
   !   real(kind=dp), intent(in) :: lambdamin, lambdamax
   !   real(kind=dp) :: resol
   !   integer :: la, N1, N2
   !   real(kind=dp) :: l0, l1, dl
 
 
   !   l1 = lambdamax
   !   l0 = lambdamin
 
   !   if (lambdamax > cont%lambda0) then
   !      N2 = Nlambda_cont_log + 1
   !      N1 =  Nlambda_cont
   !   else
   !      N2 = 0
   !      N1 =  Nlambda_cont
   !   endif
   !   cont%Nlambda = N1 + N2
   !   allocate(cont%lambda(cont%Nlambda))
 
 
   !   cont%lambda(1:N1) = span_dp(l0, cont%lambda0,N1,-1)!from lam0 to 1
 
 
   !   if (N2 > 0) then
   !      dl = abs(cont%lambda(N1)-cont%lambda(N1-1))
   !      cont%lambda(1+N1:N2+N1) = spanl_dp(cont%lambda0+dl,l1,N2,1)
   !   endif
 
 
   !   do la=1,cont%Nlambda
   !      if (cont%lambda(la) < 0) then
   !         write(*,*) "Error, lambda negative"
   !         write(*,*) "cont log"
 
   !         stop
   !      endif
   !   end do
 
   !   return
   ! end subroutine make_sub_wavelength_grid_cont_linlog
 
