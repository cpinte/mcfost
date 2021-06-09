module getlambda

  use atom_type, only : AtomicContinuum, AtomicLine, AtomType
  use atmos_type, only : atomPointerArray, Atoms, Natom, ActiveAtoms, Nactiveatoms, v_char, B_char
  use constant
  use getline, only : getnextline, MAX_LENGTH

  use parametres
  use utils, only : span, spanl, bubble_sort, spanl_dp, span_dp
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
  integer 			 :: Nlambda_max_line, Nlambda_max_cont, Nlambda_max_trans !number max of lambda for all lines
  		

  contains

  subroutine Read_wavelengths_table(lambda_table, Nlam_I)
  ! -------------------------------------------------------- !
  ! Read a wavelenth table of the form:
  ! Nregion
  ! Nlam_I(1)
  ! Nlam_I(2)
  ! ...
  ! Nlam_I(Nregion)
  ! lambda(1) (first region)
  ! ....
  ! lambda(sum(Nlam_I)) (last region, last point)
  ! -------------------------------------------------------- !

   real(kind=dp), dimension(:), allocatable, intent(inout) :: lambda_table
   integer, dimension(:), allocatable, intent(out) :: Nlam_I
   character(len=MAX_LENGTH) :: inputline, FormatLine
   integer :: Nread = 0, Nl, k, Nr, Nt

   if (.not.ltab_wavelength_image) return

   if (allocated(lambda_table)) deallocate(lambda_table) !should not happen
   write(FormatLine,'("(1"A,I3")")') "A", 256

   open(unit=1, file=TRIM(tab_wavelength_image), status='old')
   CALL getnextline(1, "#", FormatLine, inputline, Nread)
   read(inputline,*) Nl
   allocate(Nlam_I(Nl)) !wavelength per regions

   do k=1,Nl
    CALL getnextline(1, "#", FormatLine, inputline, Nread)
    read(inputline,*) Nr
    Nlam_I(k) = Nr
   end do
   Nt = sum(Nlam_I)

   !read wavelengths
   allocate(lambda_table(Nt))
   do k=1,Nt
      CALL getnextline(1, "#", FormatLine, inputline, Nread)
      read(inputline,*) lambda_table(k)
   end do

   close(1)

  return
  end subroutine Read_wavelengths_table

! - > bug format write, same for write_wavelengths_table_NLTE_lines
  subroutine write_lines_grid()
  ! -------------------------------------------------------- !
  ! write individual lines grid in the format read by
  ! write_wavelength_table
  ! for all atoms
  !
  ! This version writes the individual grid of each line
  ! after they have been created before, creating the wavelength
  ! grid of the transfer
  ! -------------------------------------------------------- !
   !!use atmos_type, only : Natom, Atoms
   type (AtomType), pointer :: atom
   type (AtomicLine) :: line
   real(kind=dp), dimension(:,:), allocatable :: lambda_table
   character(len=12) :: line_waves = "line_waves.s"
   character(len=MAX_LENGTH) :: inputline, FormatLine
   integer :: Nl, k, la, na, maxL
   integer, dimension(1000) :: Nl2
   !cound lines
   Nl = 0 !total lines of all atoms
   maxL = 0

  if (.not.allocated(Atoms(1)%ptr_atom%lines(1)%lambda)) then
   write(*,*) " Error, lambda grid of lines may not been allocated yet or anymore"
   stop
  endif

   do na=1, Natom!NActiveatoms
    atom => atoms(na)%ptr_atom!Activeatoms(na)%ptr_atom
    do k=1, atom%Nline
     Nl = Nl + 1
     if (Nl>1000) then
      write(*,*) " too many lines in write_Wavelengths_table"
      stop
     endif
     Nl2(Nl) = atom%lines(k)%Nlambda
     maxL = max(Nl2(Nl), maxL)
    enddo
    atom=>NULL()
   enddo
   !write(*,*) maxL, Nl, Nl2(1:Nl)

   allocate(lambda_table(Nl, maxL))
   write(FormatLine,'("(1"A,I3")")') "A", 256

   Nl = 0
   do na=1, Natom!NActiveatoms
    atom => atoms(na)%ptr_atom!Activeatoms(na)%ptr_atom
    do k=1, atom%Nline
     Nl = Nl + 1
    lambda_table(Nl,1:Nl2(Nl)) = atom%lines(k)%lambda

    enddo
    atom => NULL()
   enddo

   open(unit=1, file=TRIM(line_waves), status='unknown')
   write(1,*) Nl
   !write(*,*) Nl

   do k=1,Nl
    write(1,*) Nl2(Nl)
    !write(*,*) Nl2(Nl)
   end do

   Nl = 0
   do na=1, Natom!NActiveatoms
    atom => atoms(na)%ptr_atom!Activeatoms(na)%ptr_atom
    do k=1, atom%Nline
     Nl = Nl + 1
      do la=1, Nl2(Nl)
       write(1,'(1F6.3)') lambda_table(Nl, la)
       !write(*,*) Nl, NL2(Nl), lambda_table(Nl, la)
      enddo
    enddo
    atom=>NULL()
   enddo

   close(1)
   deallocate(lambda_table)

  return
  end subroutine write_lines_grid

  subroutine write_wavelengths_table_NLTE_lines(waves)
  ! -------------------------------------------------------- !
  ! write individual lines grid in the format read by
  ! write_wavelength_table
  ! for all atoms
  !
  ! This version if for the line grid taken on the final
  ! waves grid
  ! -------------------------------------------------------- !
   !!use atmos_type, only : Natom, Atoms, NactiveAtoms
   type (AtomType), pointer :: atom
   real(kind=dp), dimension(:), intent(in) :: waves
   type (AtomicLine) :: line
   real(kind=dp), dimension(:,:), allocatable :: lambda_table
   character(len=12) :: line_waves = "line_waves.s"
   character(len=MAX_LENGTH) :: inputline, FormatLine
   integer :: Nl, k, la, na, maxL
   integer, dimension(1000) :: Nl2

   if (NactiveAtoms==0) then
    write(*,*) " NLTE lines wavelength table not written, N active atoms = ", NactiveAtoms
    return
   endif

   !cound lines
   Nl = 0 !total lines of all atoms
   maxL = 0
   do na=1, NActiveatoms
    atom => Activeatoms(na)%ptr_atom
    do k=1, atom%Nline
     Nl = Nl + 1
     if (Nl>1000) then
      write(*,*) " too many lines in write_Wavelengths_table"
      stop
     endif
     Nl2(Nl) = atom%lines(k)%Nlambda
     maxL = max(Nl2(Nl), maxL)
    enddo
    atom=>NULL()
   enddo

   allocate(lambda_table(Nl, maxL))
   write(FormatLine,'("(1"A,I3")")') "A", 256

   Nl = 0
   do na=1, NActiveatoms
    atom => Activeatoms(na)%ptr_atom
    do k=1, atom%Nline
     Nl = Nl + 1
    lambda_table(k,1:Nl2(Nl)) = waves(atom%lines(k)%Nblue:atom%lines(k)%Nred)

    enddo
    atom => NULL()
   enddo

   open(unit=1, file=TRIM(line_waves), status='unknown')
   write(1,*) Nl

   do k=1,Nl
    write(1,*) Nl2(Nl)
   end do

  Nl = 0
   do na=1, NActiveatoms
    atom => Activeatoms(na)%ptr_atom
    do k=1, atom%Nline
     Nl = Nl + 1
      do la=1, Nl2(Nl)
       write(1,'(1F6.3)') lambda_table(Nl, la)
      enddo
    enddo
    atom=>NULL()
   enddo

   close(1)

   deallocate(lambda_table)

  return
  end subroutine write_wavelengths_table_NLTE_lines
  

  subroutine make_sub_wavelength_grid_cont_old(cont, lambdamin)
  ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   ! The resolution is constant in nm.
   ! lambda must be lower that lambda0 and lambda(Nlambda)=lambda0.
   ! Allocate cont%lambda.
   ! cont%alpha (cross-section of photoionisation) is not used.
  ! ----------------------------------------------------------------- !
   type (AtomicContinuum), intent(inout) :: cont
   real(kind=dp), intent(in) :: lambdamin
   real(kind=dp) :: resol
   integer :: la
   real(kind=dp) :: l0, l1

   !write(*,*) "Atom for which the continuum belongs to:", cont%atom%ID

   l1 = cont%lambda0 !cannot be larger than lambda0 ! frequency for photoionisation
   l0 = lambdamin
   cont%Nlambda = Nlambda_cont
   allocate(cont%lambda(cont%Nlambda))
   resol = (l1 - l0) / real(cont%Nlambda - 1, kind=dp)
!    write(*,*) "Continuum:", cont%lambda0, cont%j,"->",cont%i, &
!               " Resolution (nm):", resol, " lambdamin =", lambdamin

   cont%lambda(1) = l0
   do la=2,cont%Nlambda
    if (cont%lambda(la-1) < 0) then
     write(*,*) "Error, lambda negative"
     stop
    endif
    cont%lambda(la) = cont%lambda(la-1) + resol
   end do

   !does not allocate cross-section, here
  return
  end subroutine make_sub_wavelength_grid_cont_old
  

 !-> preferred to be used without level dissolution
  subroutine make_sub_wavelength_grid_cont(cont, lambdamin, lambdamax)
  ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   !linearly spaced to lambdamin, lambdamax from lambda0
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

!    cont%lambda(N1:1:-1) = (1e-15 * CLIGHT) / span(nu0,nu1,N1) * M_TO_NM
!                                                                     !Because otherwise N1 and
!                                                                     ! N1+1 points are the same
!    if (N2>0) cont%lambda(N2+N1:N1+1:-1) = (1e-15 * CLIGHT) / span(nu2,nu0*(1-0.05),N2) * M_TO_NM
!    !write(*,*) cont%lambda(N1+1:N2+N1)

!   cont%lambda(1:N1) = span_dp(l0, cont%lambda0,N1,1)!from 1 to lam0
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
   
! stop
   !does not allocate cross-section, here
  return
  end subroutine make_sub_wavelength_grid_cont

  !with level dissolution
  !linear from 0 to lambda0 then log from lambda0 to lambdamax
  subroutine make_sub_wavelength_grid_cont_linlog(cont, lambdamin, lambdamax)
  ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   ! cont%alpha (cross-section of photoionisation) is not used.
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

  !not working properly with level dissolution
  subroutine make_sub_wavelength_grid_cont_log_nu(cont, lambdamin, lambdamax)
  ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   !  -> cross-section is extrapolated beyond edge to be use with
   ! level's dissolution, if lambdamax > lambda0
   !
   ! Allocate cont%lambda.
   ! cont%alpha (cross-section of photoionisation) is not used.
   ! It is logarithmic from lambda0 to lambdamin/lambdamax
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
   allocate(cont%lambda(cont%Nlambda))

   
   nu0 = (M_TO_NM * CLIGHT / cont%lambda0)/1e15
   nu1 = (M_TO_NM * CLIGHT / lambdamin)/1e15
   nu2 = (M_TO_NM * CLIGHT / lambdamax)/1e15


  
  cont%lambda(N1:1:-1) = 1e-15 * clight / spanl_dp(nu0,nu1,N1,1) * m_to_nm

  if (N2 > 0) then
   dnu = abs(cont%lambda(N1)-cont%lambda(N1-1))
   cont%lambda(N2+N1:N1+1:-1) = 1e-15 * clight / spanl_dp(nu1+dnu,nu2,N2,1) * m_to_nm
  endif


   do la=1,cont%Nlambda
    if (cont%lambda(la) < 0) then
     write(*,*) "Error, lambda negative"
     write(*,*) "cont log"

     stop
    endif
   end do
! stop
   !does not allocate cross-section, here
  return
  end subroutine make_sub_wavelength_grid_cont_log_nu

	subroutine compute_line_bound(line, maxV)
	! ------------------------------------------------------------ !
	! Compute the line bound, from lamndamin and lambdamax
	! the total extent of the unperturbed (no velocity ) line
	!
	! to do: Bmag
	! 	Actually it is more complicated for magnetic fields
	!    And it's not possible to store the profile ??
	!    but we can first include the magnetic field to 
	!    to the max extent of the line
	! ------------------------------------------------------------ !
		type (AtomicLine), intent(inout) :: line
		real(kind=dp), intent(in) :: maxV
		real :: v_char, Nlam
		real(kind=dp) :: vB
		integer :: la, v0
		real, parameter :: Ld = 1.0 ! L% of max extent of the line without velocity field
   							   !In general, 2 * vD * (1 + aD) is enough
   							   !2.5 is for security
	!Not used if maxvel is already the maximum extent in velocity !
  
		if (line%polarizable) then
			vB = B_char * LARMOR * (line%lambda0*NM_TO_M) * dabs(line%g_lande_eff)
		else
			vB = 0.0_dp
		endif
		v0 = int(maxV)
		v_char =  Ld * real(v0 + int(vB)) !m/s
   
   
		Nlam = 2 * ( 1e-3 * v_char / hv ) + 1
		line%Nlambda = nint(Nlam)
		!but not use at this point. Will be replaced when the full grid
   							 !is computed
   
		!!Actual boundaries of the line, in absence of velocity shift and magnetic field
		!! This will be used to fine the index of the line on the global grid
		line%lambdamin = line%lambda0*(1-v_char/CLIGHT)
		line%lambdamax = line%lambda0*(1+v_char/CLIGHT)
		write(*,*) line%lambdamin, line%lambdamax, line%Nlambda


	return
	end subroutine compute_line_bound
	
	subroutine  define_local_profile_grid (line)
		type (AtomicLine), intent(inout) :: line
   		real(kind=dp), dimension(2*(Nlambda_line_c_log+Nlambda_line_w-1)-1) :: vel
		real(kind=dp) ::  vcore, vwing, v0, v1, vbroad
		integer :: Nlambda, la, Nmid
		real, parameter :: wing_to_core = 0.6 !0.3


		vbroad = maxval(line%atom%vbroad) !takes time, but done only once per line!
		vwing = line%qwing * vbroad
		
		!I compute the line bound up to qwing * vwing. But I expand the local profile
		!farther
		!This avoid allocating to much points when velocity are present
		!while allowing a large range for the interpolation ?
		line%lambdamin = line%lambda0*(1-vwing/CLIGHT)
		line%lambdamax = line%lambda0*(1+vwing/CLIGHT)
		line%Nlambda = nint(2 * ( 1e-3 * vwing / hv )  ) - 1

		!leaving here to not allocate u in case the profile is not interpolated
		!is impossible since we need to check the association profile, local_profile_interp
		!which will raise a circular compiler erros.
		!To Do: global parameter determining the type of profiles.

! 		if (line%voigt) vwing = line%qwing * vbroad * 10
		!otherwise for Gaussians not needed
		vcore = 2.5 * vbroad!wing_to_core * vwing
		!vcore = min(vcore, wing_to_core * vwing)

		v0 = - vwing
		v1 = + vwing
		vel = 0.0_dp

		!from -vwing to 0
! 		dvw = (vwing-vcore)/real(Nlambda_line_w-1,kind=dp)
! 		dvc = vcore/real(Nlambda_line_c_log-1,kind=dp)

		vel(Nlambda_line_w:1:-1) = -spanl_dp(vcore, vwing, Nlambda_line_w, -1)
! 		vel(Nlambda_line_w:1:-1) = -span_dp(vcore, vwing, Nlambda_line_w, -1)
! 		write(*,*) v0*1e-3, -vcore*1e-3
! 		write(*,*) vel(1:Nlambda_line_w)*1e-3
! 		stop
		vel(Nlambda_line_w:Nlambda_line_c_log+Nlambda_line_w-1) = span_dp(-vcore, 0.0_dp, Nlambda_line_c_log, 1)
! 		write(*,*) -vcore*1e-3, 0.0
! 		write(*,*) vel*1e-3
! 		stop		
		
		Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
   
		Nmid = Nlambda/2 + 1
		
		allocate(line%u(Nlambda)); line%u = 0.0

		line%u(1:Nmid) = vel(1:Nmid)
		line%u(Nmid+1:Nlambda) = -vel(Nmid-1:1:-1)
		
		!!full linear with Nlambda from -vwing to vwing
		!!line%u = span_dp(-vwing,vwing,Nlambda, 1)
		
! 		line%lambdamin = line%lambda0*(1+line%u(1)/CLIGHT)
! 		line%lambdamax = line%lambda0*(1+line%u(Nlambda)/CLIGHT)
		
! 		write(*,*) "line bounds:", "Nlam_local=",Nlambda
! 		write(*,*) line%lambdamin, line%lambdamax
! 		write(*,*) line%u(1)*1e-3, line%u(Nlambda)*1e-3
		
! 		line%Nlambda = nint(2 * ( 1e-3 * vwing / hv )  ) - 1
! 		write(*,*) " Nlam_grid=", line%Nlambda," maxv=", vwing*1e-3
! 
! 		do la=1,Nlambda
! 			write(*,*) la, line%u(la), vel(la)
! 		enddo
! 		stop
	
	return
	end subroutine define_local_profile_grid
	
	subroutine  define_local_gauss_profile_grid (atom)
		type (AtomType), intent(inout) :: atom
   		real(kind=dp), dimension(2*(Nlambda_line_c_log+Nlambda_line_w-1)-1) :: vel
		real(kind=dp) ::  vcore, vwing, v0, v1, vbroad
		integer :: Nlambda, la, Nmid, kr
		real, parameter :: wing_to_core = 0.6 !0.3


		vbroad = maxval(atom%vbroad)
		vwing = 7.0 * vbroad
		vcore = 4.0 * vbroad

		v0 = - vwing
		v1 = + vwing
		vel = 0.0_dp


		vel(Nlambda_line_w:1:-1) = -spanl_dp(vcore, vwing, Nlambda_line_w, -1)
		vel(Nlambda_line_w:Nlambda_line_c_log+Nlambda_line_w-1) = span_dp(-vcore, 0.0_dp, Nlambda_line_c_log, 1)	
		
		Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
		
		Nmid = Nlambda/2 + 1
		
		allocate(atom%ug(Nlambda)); atom%ug = 0.0

		atom%ug(1:Nmid) = vel(1:Nmid)
		atom%ug(Nmid+1:Nlambda) = -vel(Nmid-1:1:-1)
		
		!changing Nlambda for each gauss line
! 		do kr=1, atom%Nline
! 			if (.not.atom%lines(kr)%voigt) then
! 				write(*,*) " Nlambda for gauss line changed from ", atom%lines(kr)%Nlambda, " to ", Nlambda
! 				atom%lines(kr)%Nlambda = Nlambda
! 			endif
! 		enddo
		
	
	return
	end subroutine define_local_gauss_profile_grid

	subroutine make_sub_wavelength_grid_line_lin(line, vD, aD)
	! ------------------------------------------------------------ !
	! Make an individual wavelength grid for the AtomicLine line.
	! The wavelength grid is symmetric wrt lambda0.
	! v_char need to be sufficiently large enough to encompass shift
	! of the order of the maximum velocity.
	! ------------------------------------------------------------ !
		type (AtomicLine), intent(inout) :: line
		real(kind=dp), intent(in) :: vD, aD !maximum thermal width of the atom in m/s
		real(kind=dp) :: v_char, dlam
		real(kind=dp) :: vB
		real(kind=dp) :: lam0, lam1
		integer :: la, Nlambda, Nmid
		real, parameter :: L = 10.01 ! L% of max extent contains the line
		real(kind=dp), dimension(Nlambda_line_c) :: xlam !dimension(2*(Nlambda_line_c+Nlambda_line_w-1)-1)
  

		vB = 0d0
		if (line%polarizable) vB =  &
			B_char * LARMOR * (line%lambda0*NM_TO_M) * dabs(line%g_lande_eff)


		v_char =  min(L * vD * (1.0 + aD**(-1.0)) + L * vB,3e5)


		xlam = 0d0
		!lam0 = line%lambda0*(1-v_char/CLIGHT)
		!lam1 = line%lambda0*(1+v_char/CLIGHT)
		line%lambdamin = line%lambda0*(1-v_char/CLIGHT)
		line%lambdamax = line%lambda0*(1+v_char/CLIGHT)

		Nlambda = Nlambda_line_c
		line%Nlambda = Nlambda
		if (mod(line%Nlambda,2)==0) line%Nlambda = line%Nlambda + 1

		allocate(line%lambda(line%Nlambda))
		line%lambda(1) = line%lambdamin
		dlam = (line%lambdamax-line%lambdamin) / (real(line%Nlambda-1))
		do la=2,line%Nlambda
		if (line%lambda(la-1) < 0) then
			write(*,*) "Error lambda nagative"
			write(*,*) "line lin"
			stop
		endif
		line%lambda(la) = line%lambda(la-1) + dlam
		enddo


	return
	end subroutine make_sub_wavelength_grid_line_lin

  subroutine make_sub_wavelength_grid_line(line, vD, aD)
  ! ------------------------------------------------------------ !
   ! Make an individual wavelength grid for the AtomicLine line.
   ! The wavelength grid is symmetric wrt lambda0.
   ! It is by default, logarithmic in the wing and linear in the
   ! core.
   !
   ! I recommend to not use it in case of velocity shifts larger
   ! than the doppler width of the line !
   !
  ! ------------------------------------------------------------ !
   type (AtomicLine), intent(inout) :: line
   real(kind=dp), intent(in) :: vD !maximum thermal width of the atom in m/s
   real(kind=dp), intent(in) :: aD !maximum damping wrt vbroad
   real(kind=dp) :: v_char, dvc, dvw
   real(kind=dp) :: vcore, vB, v0, v1!km/s
   integer :: la, Nlambda, Nmid
   real, parameter :: wing_to_core = 0.01, L = 70 !0.01, 70
   		!if it is the max, then L is close to 1, if it is the min, L >> 1, if it is the mean etc..
   !!integer, parameter :: Nc = 51, Nw = 7 !ntotal = 2*(Nc + Nw - 1) - 1
   real(kind=dp), dimension(2*(Nlambda_line_c_log+Nlambda_line_w-1)-1) :: vel
   !!real(kind=dp), dimension(2*(Nc+Nw-1)-1) :: vel !Size should be 2*(Nc+Nw-1)-1
   													 !if error try, 2*(Nc+Nw)


   vB = 0d0
   if (line%polarizable) vB =  &
   				B_char * LARMOR * (line%lambda0*NM_TO_M) * dabs(line%g_lande_eff)

   !!v_char = L * (v_char + 2d0*vD + vB) !=maximum extension of a line
   !!v_char = v_char + L * (vB  + vD)
   !!vcore = v_char * wing_to_core
   !transition between wing and core in velocity
   !!vcore = L * v_char * wing_to_core ! == fraction of line extent

   vcore = 20 * (vD + aD**(-1.0)) 
   v_char = max(v_char + L * vD * (1.0 + aD**(-1.0)) + vB, 2000e3)
   !!vcore = v_char * wing_to_core

   !for compatibiliy
   line%lambdamin = line%lambda0 * (1 - v_char / clight)
   line%lambdamax = line%lambda0 * (1 + v_char / clight)
   
   v0 = -v_char !* L
   v1 = +v_char !* L
   vel = 0d0

   !from -v_char to 0
   dvw = (v_char-vcore)/real(Nlambda_line_w-1,kind=dp) !(L * v_char-vcore)/real(Nw-1,kind=dp), old
   dvc = vcore/real(Nlambda_line_c_log-1,kind=dp)


!! Log wing
   !should be log from vcore to v0 not the opposite
   !v0 is < 0 but spanl takes the abs
!   vel(1:Nw) = -real(spanl(real(v0), real(vcore), Nw),kind=dp)
   vel(Nlambda_line_w:1:-1) = -real(spanl(real(vcore), real(v0), Nlambda_line_w),kind=dp)
!! end scale of wing points
!   vel(Nw:Nw-1+Nc) = -real(span(real(vcore), real(0.), Nc+1),kind=dp)
!   write(*,*) Nw, vel(Nw), vcore
   !vel(Nw) = -vcore!should be okey at the numerical precision
   !la goes from 1 to Nw + Nc -1 total number of points.
   !if Nc=101 and Nw =11, 111 velocity points,because the point vel(11) of the wing grid
   !is shared with the point vel(1) of the core grid.
   do la=Nlambda_line_w+1, Nlambda_line_c_log+Nlambda_line_w-1 !Half line core
    vel(la) = vel(la-1) + dvc
   end do

   !! Just a check here, maybe forcing the mid point to be zero is brutal
   !! but by construction it should be zero !
   !if (dabs(vel(Nw+Nc-1)) <= 1d-7) vel(Nw+Nc-1) = 0d0
   if (dabs(vel(Nlambda_line_w+Nlambda_line_c_log-1)) /= 0d0) vel(Nlambda_line_w+Nlambda_line_c_log-1) = 0d0
   if (vel(Nlambda_line_w+Nlambda_line_c_log-1) /= 0) write(*,*) 'Vel(Nw+Nc-1) should be 0.0'

  !number of points from -vchar to 0 is Nw+Nc-1, -1 because otherwise we count
  ! 2 times vcore which is shared by the wing (upper boundary) and core (lower boundary) grid.
  ! Total number of points is 2*(Nw+Nc-1) but here we count 2 times lambda0., therefore
  ! we remove 1 point.
   Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c_log - 1) - 1
   line%Nlambda = Nlambda
   Nmid = Nlambda/2 + 1 !As Nlambda is odd '1--2--3', Nmid = N/2 + 1 = 2, because 3/2 = 1
   						!because the division of too integers is the real part.
   allocate(line%lambda(line%Nlambda))

   line%lambda(1:Nmid) = line%lambda0*(1d0 + vel(1:Nmid)/CLIGHT)
   line%lambda(Nmid+1:Nlambda) = line%lambda0*(1d0 -vel(Nmid-1:1:-1)/CLIGHT)

   if (line%lambda(Nmid) /= line%lambda0) write(*,*) 'Lambda(Nlambda/2+1) should be lambda0'

   do la=1,line%Nlambda
     if (line%lambda(la) < 0) then
      write(*,*) "maxv=",maxval(vel), "vcore=",vcore, "v0",v0, "Nlam=",Nlambda_line_w
      write(*,*) "Error lambda negative"
      write(*,*) "line log", vD, aD, la, line%lambda(la)
      stop
     endif
   enddo

  return
  end subroutine make_sub_wavelength_grid_line

!   subroutine make_wavelength_grid(Natom, Atoms, inoutgrid, Ntrans, wl_ref)
!   use math, only : locate
!   ! --------------------------------------------------------------------------- !
!    ! construct and sort a wavelength grid for atomic line radiative transfer.
!    ! Computes also the edge of a line profile: Nblue and Nred.
!    ! The grid is built by merging the individual wavelength grids of each
!    ! transitions. Duplicates are removed and therefore changes the value of
!    ! line%Nlambda and continuum%Nlambda.
!    ! line%lambda and continuum%lambda are useless now. Except that
!    ! continuu%lambda is still used in .not.continuum%Hydrogenic !!
!   ! --------------------------------------------------------------------------- !
!    type (atomPointerArray), intent(inout), dimension(Natom) :: Atoms
!    integer, intent(in) :: Natom
!    real(kind=dp), intent(in) :: wl_ref
!    integer, intent(out) :: Ntrans !Total number of transitions (cont + line)
!    ! output grid. May contain values that are added to the final list before
!    ! deallocating the array. It is reallocated when the final list is known.
!    real(kind=dp), allocatable, dimension(:), intent(inout) :: inoutgrid
!    ! temporary storage for transitions, to count and add them.
!    integer, parameter :: MAX_TRANSITIONS = 50000
!    type (AtomicLine), allocatable, dimension(:) :: alllines
!    type (AtomicContinuum), allocatable, dimension(:) :: allcont
!    integer :: kr, kc, n, Nspect, Nwaves, Nlinetot, Nctot
!    integer :: la, nn, nnn !counters: wavelength, number of wavelengths, line wavelength
!    integer :: Nred, Nblue, Nlambda_original!read from model
!    real(kind=dp), allocatable, dimension(:) :: tempgrid
!    integer, allocatable, dimension(:) :: sorted_indexes
!    real(kind=dp) :: l0, l1 !ref wavelength of each transitions
!    character(len=10) :: lam_unit
! 
!    !write(*,*) ' Defining the nLTE wavelength grid, using ', Nlambda_cont,' points for each continuum, and ', &
!    ! 2*(Nlambda_line_w+Nlambda_line_c-1)-1, " points for each line."
! 
!    nn = 0
!    allocate(alllines(MAX_TRANSITIONS), allcont(MAX_TRANSITIONS))!stored on heap
!    ! if allocated inoutgrid then add to the number of
!    ! wavelength points to the points of the inoutgrid.
!    Nspect = 0
!    if (allocated(inoutgrid)) Nspect = size(inoutgrid)
! 
!    Nlinetot = 0
!    Nctot = 0
!    do n=1,Natom
!    !unlike RH, even passive atoms have dedicated wavelength grids
!    ! Count number of total wavelengths
!    do kr=1,Atoms(n)%ptr_atom%Nline
!     Nspect = Nspect + Atoms(n)%ptr_atom%lines(kr)%Nlambda
!     Nlinetot = Nlinetot + 1
!     if (Nlinetot > MAX_TRANSITIONS) then
!      write(*,*) "too many transitions"
!      stop
!     end if
!     alllines(Nlinetot) = Atoms(n)%ptr_atom%lines(kr)
!    end do
!    do kc=1,Atoms(n)%ptr_atom%Ncont
!     Nspect = Nspect + Atoms(n)%ptr_atom%continua(kc)%Nlambda
!     Nctot = Nctot + 1
!     if (Nctot > MAX_TRANSITIONS) then
!      write(*,*) "too many transitions"
!      stop
!     end if
!     allcont(Nctot) = Atoms(n)%ptr_atom%continua(kc)
!    end do
!   end do ! end loop over atoms
! 
! 
!   ! add ref wavelength if any and allocate temp array
!   if (wl_ref > 0.) then
!    Nspect = Nspect + 1
!    nn = 1
!    allocate(tempgrid(Nspect))
!    tempgrid(1)=wl_ref
!   else
!    allocate(tempgrid(Nspect))
!   end if
! 
!   Ntrans = Nlinetot + Nctot
!   write(*,*) "Adding ", Nspect," wavelengths for a total of", Ntrans, &
!              " transitions."
!   if (allocated(inoutgrid)) write(*,*) "  ->", size(inoutgrid)," input wavelengths"
! 
! !   allocate(Nred_array(Ntrans), Nblue_array(Ntrans), Nmid_array(Ntrans))
! 
!   ! add wavelength from mcfost inoutgrid if any
!   ! and convert it to nm, then deallocate
!   if (allocated(inoutgrid)) then
!    do la=1, size(inoutgrid)
!     nn = nn + 1
!     tempgrid(nn) = inoutgrid(la) !nm or convert to nm here
!    end do
!    deallocate(inoutgrid)
!   end if
! 
!   ! start adding continua and lines wavelength grid
!   nnn = 0!it is just an indicative counter
!   		 ! because nn is dependent on wl_ref and on the inoutgrid
!   		 ! it is filled.
!   do kc=1,Nctot
!    do la=1,allcont(kc)%Nlambda
!     nn = nn + 1
!     nnn = nnn + 1
!     tempgrid(nn) = allcont(kc)%lambda(la)
!    end do
!   end do
!   write(*,*) "  ->", nnn," continuum wavelengths"
!   nnn = 0
!   do kr=1,Nlinetot
!    do la=1,alllines(kr)%Nlambda
!     nn = nn + 1
!     nnn = nnn + 1
!     tempgrid(nn) = alllines(kr)%lambda(la)
!    end do
!   end do
!   write(*,*) "  ->", nnn," line wavelengths"
!   if (wl_ref > 0)   write(*,*) "  ->", 1," reference wavelength at", wl_ref, 'nm'
!   ! sort wavelength
!   !!CALL sort(tempgrid, Nspect)
!   !this should work ?
!   allocate(sorted_indexes(Nspect))
!   sorted_indexes = bubble_sort(tempgrid)
!   tempgrid(:) = tempgrid(sorted_indexes)
! 
!   !check for dupplicates
!   !tempgrid(1) already set
!   Nwaves = 2
! !   write(*,*) "l0", tempgrid(1)
!   do la=2,Nspect
! !   write(*,*) "dlam>0", tempgrid(la)-tempgrid(la-1)
!    if (tempgrid(la) > tempgrid(la-1)) then
!     tempgrid(Nwaves) = tempgrid(la)
! !     write(*,*) tempgrid(la), tempgrid(la-1), Nwaves+1, la, Nspect, tempgrid(Nwaves)
!     Nwaves = Nwaves + 1
!    end if
!   end do
! 
!   Nwaves = Nwaves-1 ! I don't understand yet but it works
!   					! Maybe, if only 1 wavelength, Nwaves = 2 - 1 = 1 (and not 2 as it
!   					! is implied by the above algorithm...)
!   					! or if 2 wavelengths, Nwaves = 2 = 3 - 1 (not 3 again etc)
!   					! and so on, actually I have 1 wavelength more than I should have.
!   write(*,*) Nwaves, " unique wavelengths: ", Nspect-Nwaves," eliminated lines"
! 
!   ! allocate and store the final grid
!   allocate(inoutgrid(Nwaves))
!   do la=1,Nwaves
!    inoutgrid(la) = tempgrid(la)
!   end do
! 
!   lam_unit = "nm"
!   l0 = minval(inoutgrid); l1 = maxval(inoutgrid)
!   if (l1 > 1500.) then
!    l1 = l1 *1e-4
!    lam_unit = "microns"
! !   else if (l1 > 1e6) then
! !    l1 = 10000000./l1
! !    lam_unit = "cm^-1"
!   else if (l1 > 1e6) then
!    l1 = l1 * 1e-9 * 1e3
!    lam_unit = "mm"
!   else if (l1 > 1e7) then
!    l1 = l1 * 1e-9 * 1e2
!    lam_unit = "cm"
!   endif
!   write(*,*) "Wavelength grid:", nint(l0)," nm",nint(l1),lam_unit  
! !   write(*,*) "Wavelength grid:", real(minval(inoutgrid)),' nm to',real(maxval(inoutgrid)), ' nm'
! 
!   !!should not dot that but error somewhere if many atoms
!   !!CALL sort(inoutgrid, Nwaves)
! 
!   !free some space
!   deallocate(tempgrid, alllines, allcont, sorted_indexes)
! 
! !   Now replace the line%Nlambda and continuum%Nlambda by the new values.
! !   we do that even for PASSIVE atoms
! !   deallocate line%lambda because it does not correspond to the new
! !   Nblue and Nlambda anymore.
! ! However, cont%lambda is need for interpolation of cont%alpha on the new grid if
! ! cont is not hydrogenic
! !   nn = 1
!   do n=1,Natom
!    !first continuum transitions
! !   write(*,*) " ------------------------------------------------------------------ "
!    do kc=1,Atoms(n)%ptr_atom%Ncont
!     Nlambda_original = Atoms(n)%ptr_atom%continua(kc)%Nlambda
!     l0 = Atoms(n)%ptr_atom%continua(kc)%lambda(1)
!     l1 = Atoms(n)%ptr_atom%continua(kc)%lambda(Nlambda_original) !on the subgrid
! !    Nred = locate(inoutgrid,l1)
! !     Nblue = locate(inoutgrid,l0)
! !      write(*,*) locate(inoutgrid,l0), locate(inoutgrid,l1), l0, l1!, Nblue, Nred
!     Atoms(n)%ptr_atom%continua(kc)%Nblue = locate(inoutgrid,l0)
!     Atoms(n)%ptr_atom%continua(kc)%Nred = locate(inoutgrid,l1)
!     Atoms(n)%ptr_atom%continua(kc)%Nlambda = Atoms(n)%ptr_atom%continua(kc)%Nred - &
!                                     Atoms(n)%ptr_atom%continua(kc)%Nblue + 1
! !      write(*,*) Atoms(n)%ID, " continuum:",kr, " Nlam_ori:", Nlambda_original, &
! !      " l0:", l0, " l1:", l1, " Nred:",  Atoms(n)%continua(kc)%Nred, &
! !        " Nblue:", Atoms(n)%continua(kc)%Nblue, " Nlambda:", Atoms(n)%continua(kc)%Nlambda, &
! !        " Nblue+Nlambda-1:", Atoms(n)%continua(kc)%Nblue + Atoms(n)%continua(kc)%Nlambda - 1
! !! For continuum transitions, lambda0 is at Nred, check definition of the wavelength grid
! !! which means that cont%Nmid = locate(inoutgrid, lam(Nred)+lam(Nb)/(Nlambda))
! !! and l1, lam(Nlambda) = lambda0
!     Atoms(n)%ptr_atom%continua(kc)%Nmid = locate(inoutgrid,0.5*(l0+l1))
!     Atoms(n)%ptr_atom%continua(kc)%N0 = locate(inoutgrid, Atoms(n)%ptr_atom%continua(kc)%lambda0)
! 
!     !if (Atoms(n)%ptr_atom%continua(kc)%Hydrogenic) & 
!     deallocate(atoms(n)%ptr_atom%continua(kc)%lambda)
!     !table of photoion kept on %lambda_file and alpha_file
! !!deprecated
! !     CALL fillPhotoionisationCrossSection(Atoms(n)%ptr_atom, kc, &
! !     	Atoms(n)%ptr_atom%continua(kc)%lambda,Nwaves, inoutgrid)
! !
! !     !TEST
! !     Atoms(n)%ptr_atom%continua(kc)%Nlambda = size( Atoms(n)%ptr_atom%continua(kc)%alpha)
!       !!-> also deallocates cont%lambda which is the original lambda read from file
!     !allocate(Atoms(n)%continua(kc)%lambda(Atoms(n)%continua(kc)%Nlambda))
!     !Atoms(n)%continua(kc)%lambda(Atoms(n)%continua(kc)%Nblue:Atoms(n)%continua(kc)%Nred) &
!     ! = inoutgrid(Atoms(n)%continua(kc)%Nblue:Atoms(n)%continua(kc)%Nred)
! !!!     Nred_array(nn) = Atoms(n)%continua(kc)%Nred
! !!!     Nmid_array(nn) = Atoms(n)%continua(kc)%Nmid
! !!!     Nblue_array(nn) = Atoms(n)%continua(kc)%Nblue
! !!!     nn= nn + 1
!    end do
!    !then bound-bound transitions
!    do kr=1,Atoms(n)%ptr_atom%Nline
!     Nlambda_original = Atoms(n)%ptr_atom%lines(kr)%Nlambda
!     l0 = Atoms(n)%ptr_atom%lines(kr)%lambda(1)
!     l1 = Atoms(n)%ptr_atom%lines(kr)%lambda(Nlambda_original)
! !    Nred = locate(inoutgrid,l1)
! !    Nblue = locate(inoutgrid,l0)
! !     write(*,*) locate(inoutgrid,l0), locate(inoutgrid,l1), l0, l1!, Nblue, Nred
!     Atoms(n)%ptr_atom%lines(kr)%Nblue = locate(inoutgrid,l0)!Nblue
!     Atoms(n)%ptr_atom%lines(kr)%Nred = locate(inoutgrid,l1)
!     Atoms(n)%ptr_atom%lines(kr)%Nlambda = Atoms(n)%ptr_atom%lines(kr)%Nred - &
!                                  Atoms(n)%ptr_atom%lines(kr)%Nblue + 1
! !     write(*,*) Atoms(n)%ID, " line:",kr, " Nlam_ori:", Nlambda_original, &
! !     " l0:", l0, " l1:", l1, " Nred:",  Atoms(n)%lines(kr)%Nred, &
! !       " Nblue:", Atoms(n)%lines(kr)%Nblue, " Nlambda:", Atoms(n)%lines(kr)%Nlambda, &
! !       " Nblue+Nlambda-1:", Atoms(n)%lines(kr)%Nblue + Atoms(n)%lines(kr)%Nlambda - 1
!     Atoms(n)%ptr_atom%lines(kr)%Nmid = locate(inoutgrid,Atoms(n)%ptr_atom%lines(kr)%lambda0)
!     deallocate(Atoms(n)%ptr_atom%lines(kr)%lambda) !does not correpond to the new grid, indexes might be wrong
!     !allocate(Atoms(n)%ptr_atom%lines(kr)%lambda(Atoms(n)%ptr_atom%lines(kr)%Nlambda))
!     !Atoms(n)%lines(kr)%lambda(Atoms(n)%lines(kr)%Nblue:Atoms(n)%lines(kr)%Nred) &
!     ! = inoutgrid(Atoms(n)%lines(kr)%Nblue:Atoms(n)%lines(kr)%Nred)
! !!!     Nred_array(nn) = Atoms(n)%lines(kr)%Nred
! !!!     Nmid_array(nn) = Atoms(n)%lines(kr)%Nmid
! !!!     Nblue_array(nn) = Atoms(n)%lines(kr)%Nblue
! !!!     nn = nn + 1
!    end do
! !     write(*,*) " ------------------------------------------------------------------ "
!   end do !over atoms
! 
!   return
!   end subroutine make_wavelength_grid

  subroutine adjust_wavelength_grid(old_grid, lambda, Lam_region)
   ! ------------------------------------------ !
    ! Reallocate wavelengths and indexes arrays
    ! to compute images on a user defined grid
   ! ------------------------------------------ !
   use math, only : locate
   use atmos_type, only : realloc_Transitions
   real(kind=dp), dimension(:), intent(in) :: old_grid
   integer, dimension(:), intent(in) :: lam_region
   !!type (atomPointerArray), dimension(:), intent(inout) :: Atoms
   real(kind=dp), dimension(:), intent(inout) :: lambda
   real(kind=dp), dimension(size(lambda)) :: lambda_us
   integer :: Nwaves, n, kr, kc, Nlambda_original, Nblue, Nred, Natom, ll, lll, alloc_status
   real(kind=dp) :: l0, l1 !ref wavelength of each transitions
   real(kind=dp) :: x0, x1,maxV !bound of the new grid
   integer, dimension(:), allocatable :: sorted_indexes
   logical, dimension(:), allocatable :: trans_contribute
   type (AtomicContinuum), dimension(:), allocatable :: conta
   type (AtomicLine), dimension(:), allocatable :: lines
   logical :: in_chan

   Nwaves = size(lambda)
   !check lambda is sorted ?
   !--> moved after the test over transitions now
   allocate(sorted_indexes(Nwaves), stat=alloc_status)
   if (alloc_status > 0) then
    call error("Allocation error, sorted_indexes")
   endif
   lambda_us(:) = lambda(:)
   sorted_indexes = bubble_sort(lambda)
   lambda(:) = lambda(sorted_indexes)
   x0 = minval(lambda); x1 = maxval(lambda)
   hv = 1e3 !m/s
   do ll=2, Nwaves
   	hv = min(hv, real((lambda(ll) - lambda(ll-1)) / lambda(ll-1) * clight) * 1e-3)
   enddo
   write(*,*) " -> New grid resolution for image in km/s", hv !should be integer or half integer and constant..
   !Realloc space for atoms
   !we need to test if a transition is on the new grid or not. Because the final grid
   !is not the sum of the individual grid, some transitions can be neglected because
   !they are out of range
   do n=1,Natom

    allocate(trans_contribute(atoms(n)%ptr_atom%Ntr)); trans_contribute(:)=.true.!by default

    do kc=1,Atoms(n)%ptr_atom%Ncont
     !on the old_Grid (i.e., the grid for NLTE which is built using all transitions)
     Nlambda_original = Atoms(n)%ptr_atom%continua(kc)%Nlambda !on the old_grid
     Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue
     Nred = Atoms(n)%ptr_atom%continua(kc)%Nred
     l0 = Atoms(n)%ptr_atom%continua(kc)%lambdamin
     !new with level's dissolution, equivalent to old case if lambdamax=lambda0
     l1 = Atoms(n)%ptr_atom%continua(kc)%lambdamax!Atoms(n)%ptr_atom%continua(kc)%lambda0

     in_chan = .false. !equivalent of trans_contribute so be smarter please
     ll = 0
     region_loop : do kr=1, size(lam_region)
     !relative index of regions
     ll = 1 + ll; lll = sum(lam_region(1:kr))
     x0 = minval(lambda_us(ll:lll)); x1 = maxval(lambda_us(ll:lll))
     if (l1 <= x0.or. l0 >= x1) then
      in_chan = .false.
     else
      in_chan = .true.
      exit region_loop !because if in one region no need to test the others
     end if
     ll = ll + sum(lam_region(1:kr))
     end do region_loop
     if (in_chan) then
      Atoms(n)%ptr_atom%continua(kc)%Nred = locate(lambda,l1) ! closest value return by locate is Nwaves if l1>lambda(Nwaves)
      Atoms(n)%ptr_atom%continua(kc)%Nblue = locate(lambda,l0) !
      Nred = Atoms(n)%ptr_atom%continua(kc)%Nred; Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue
     else
      Atoms(n)%ptr_atom%continua(kc)%Nred = -99
      Atoms(n)%ptr_atom%continua(kc)%Nblue = -99
     end if


     Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue; Nred = Atoms(n)%ptr_atom%continua(kc)%Nred
     if (Nred==-99.and.Nblue==-99) then
      Atoms(n)%ptr_atom%continua(kc)%Nlambda = -99

      Atoms(n)%ptr_atom%continua(kc)%Nmid = -99

      !!Atoms(n)%ptr_atom%continua(kc)%lcontrib_to_opac=.false.
      trans_contribute(atoms(n)%ptr_atom%Nline+kc) = .false.
      Atoms(n)%ptr_atom%at(atoms(n)%ptr_atom%Nline+kc)%lcontrib_to_opac=.false.
      write(*,*) " :: b-f transition", Atoms(n)%ptr_atom%continua(kc)%j,"->",Atoms(n)%ptr_atom%continua(kc)%i,&
       " for atom ",Atoms(n)%ptr_atom%ID, l0,"-",l1," not counted." !, " removed."
     else
      Atoms(n)%ptr_atom%continua(kc)%Nlambda = Atoms(n)%ptr_atom%continua(kc)%Nred - &
                                 Atoms(n)%ptr_atom%continua(kc)%Nblue + 1

      Atoms(n)%ptr_atom%continua(kc)%Nmid = locate(lambda,lambda(Atoms(n)%ptr_atom%continua(kc)%Nlambda)/2+1)
      Atoms(n)%ptr_atom%continua(kc)%N0 = locate(lambda, Atoms(n)%ptr_atom%continua(kc)%lambda0)

     end if
    end do

    !then bound-bound transitions
    do kr=1,Atoms(n)%ptr_atom%Nline
    
     Nlambda_original = Atoms(n)%ptr_atom%lines(kr)%Nlambda !on the old_grid
     Nblue = Atoms(n)%ptr_atom%lines(kr)%Nblue
     Nred = Atoms(n)%ptr_atom%lines(kr)%Nred
     l0 = Atoms(n)%ptr_atom%lines(kr)%lambdamin!old_grid(Nblue)
     l1 = Atoms(n)%ptr_atom%lines(kr)%lambdamax!old_grid(Nred)
     
     !recompute bounds with new shift ? 
     !call compute_line_bound(Atoms(n)%ptr_atom%lines(kr), maxV)
     maxV = hvel_nlte * dabs(l0-Atoms(n)%ptr_atom%lines(kr)%lambda0)/Atoms(n)%ptr_atom%lines(kr)%lambda0
     l0 = l0 + maxV/hvel_nlte - maxV/hv
     l1 = l1 - maxV/hvel_nlte + maxV/hv


     in_chan = .false.
     ll = 0
     !!write(*,*) size(lam_region), l0, l1, Atoms(n)%ptr_atom%lines(kr)%lambda0
     region_loop_l : do kc=1, size(lam_region)
     !relative index of regions
     ll = 1 + ll; lll = sum(lam_region(1:kc))
     !!write(*,*) kc, ll, lll, lam_region(kc), lll-ll+1
     x0 = minval(lambda_us(ll:lll)); x1 = maxval(lambda_us(ll:lll))
     !!write(*,*) x0, x1
     if (l1 <= x0.or. l0 >= x1) then
      in_chan = .false.
     else
      in_chan = .true.
      exit region_loop_l !because if in one region no need to test the others
     end if
     ll = sum(lam_region(1:kc))
     end do region_loop_l
     !!write(*,*) in_chan

     if (in_chan) then
      Atoms(n)%ptr_atom%lines(kr)%Nblue = locate(lambda,l0) ! closest value return by locate is Nwaves if l1>lambda(Nwaves)
      Atoms(n)%ptr_atom%lines(kr)%Nred = locate(lambda,l1) !
     else
      Atoms(n)%ptr_atom%lines(kr)%Nblue = -99
      Atoms(n)%ptr_atom%lines(kr)%Nred = -99
     end if

     Nblue = Atoms(n)%ptr_atom%lines(kr)%Nblue; Nred = Atoms(n)%ptr_atom%lines(kr)%Nred
     if (Nred==-99.and.Nblue==-99) then
      Atoms(n)%ptr_atom%lines(kr)%Nlambda = -99

      Atoms(n)%ptr_atom%lines(kr)%Nmid = -99

      !!Atoms(n)%ptr_atom%lines(kr)%lcontrib_to_opac=.false.
      trans_contribute(kr) = .false.
      Atoms(n)%ptr_atom%at(kr)%lcontrib_to_opac=.false.
      write(*,*) " :: b-b transition", Atoms(n)%ptr_atom%lines(kr)%j,"->",Atoms(n)%ptr_atom%lines(kr)%i,&
       " for atom " ,Atoms(n)%ptr_atom%ID, l0,"-",l1," not counted."!, " removed."
     else
      Atoms(n)%ptr_atom%lines(kr)%Nlambda = Atoms(n)%ptr_atom%lines(kr)%Nred - &
                                 Atoms(n)%ptr_atom%lines(kr)%Nblue + 1

      Atoms(n)%ptr_atom%lines(kr)%Nmid = locate(lambda,Atoms(n)%ptr_atom%lines(kr)%lambda0)
     end if
    if (allocated(Atoms(n)%ptr_atom%lines(kr)%lambda)) &
    	deallocate(Atoms(n)%ptr_atom%lines(kr)%lambda)
    end do

     CALL realloc_transitions(Atoms(n)%ptr_atom, count(trans_contribute), trans_contribute)
     deallocate(trans_contribute)
   end do !over atoms

  return
  end subroutine adjust_wavelength_grid
  
  !separate cont grid and total grid and so are Nblue, Nred, Nlambda and N0 ect, define wrt their grid, not the same !
  !Nblue (cont_grid) cannot be used on total_grid for instance
	subroutine make_wavelength_grid_new(wl_ref, dvmax, outgrid, Ntrans, cont_grid)
		use math, only : locate
 
		! Create a wavelength grid around group of lines (at least 1 line in a group)
		! with a constant velocity spacing for each group min(lambda(group)) - dvmax to max(lambda(group))+dvmax.
		! Continuum points are added to the grid, outside lines group.
		!This allows to not break the line sampling inside each group.
   
		!A group of lines are lines that may overlap because they are close.
		!Take into account dvmax to compute the overlap.

		!type (atomPointerArray), intent(inout), dimension(Natom) :: Atoms
		!integer, intent(in) :: Natom
		real(kind=dp), intent(in) :: wl_ref, dvmax
		integer, intent(out) :: Ntrans
		real(kind=dp), allocatable, dimension(:), intent(out) :: outgrid, cont_grid
		real(kind=dp), dimension(:), allocatable :: cont_waves, line_waves
		integer, parameter :: MAX_GROUP_OF_LINES = 1000
		integer :: Nwaves
		integer :: Ngroup, Nlam, Nlambda_cont, Ncont, Nline_per_group(MAX_GROUP_OF_LINES)
		real(kind=dp), dimension(MAX_GROUP_OF_LINES) :: group_blue, group_red
		real(kind=dp), dimension(:), allocatable :: all_lamin, all_lamax, tmp_grid, delta_lam
		integer, dimension(:), allocatable :: sorted_indexes, Nlambda_per_group
		integer :: Nspec_cont, Nspec_line, Nremoved, Nwaves_cont
		integer :: n, kr, la, lac, shift, alloc_status, Nmore_cont_freq, check_new_freq
		real(kind=dp) :: lambda_max, lambda_min, l0, l1, max_cont, delta_v
		real(kind=dp) :: vmin, vmax
		real(kind=dp), dimension(:), allocatable :: more_cont_waves, corr_hv
		type (AtomType), pointer :: atom
		character(len=15) :: lam_unit
		logical :: add_cont, lthere_is_lines = .false.

		if (allocated(outgrid)) then
			write(*,*) " Cannot use non-empty grid for this wavelength grid !"
			deallocate(outgrid)
		endif
		Ntrans = 0
		!maximum extent of lines with or without dvmax (i.e., in rest frame)
		delta_v = dvmax +  0 * hv * 1e3 !m/s

		!maximum and minimum wavelength for only lines, including max velocity field
		!Count Number of transitions and number of lines
		!lambda_max = 0
		!lambda_min = 1d100
		Nlam = 0
		Nlambda_cont = 0
		Ncont = 0
		do n=1, Natom
			atom => atoms(n)%ptr_atom
			do kr=1,atom%Ncont
				Nlambda_cont = Nlambda_cont + atom%continua(kr)%Nlambda
! 				lambda_max = min(lambda_max, atom%continua(kr)%lambdamax)
! 				lambda_min = min(lambda_min, atom%continua(kr)%lambdamin)

			enddo
			Ntrans = Ntrans + atom%Ntr!atom%Ncont + atom%Nline
    
			Nlam = Nlam + atom%Nline
			Ncont = Ncont + atom%Ncont

			do kr=1,atom%Nline
				if (allocated(atom%lines(kr)%lambda)) deallocate(atom%lines(kr)%lambda)
! 				lambda_max = max(lambda_max, atom%lines(kr)%lambdamax)
! 				lambda_min = min(lambda_min, atom%lines(kr)%lambdamin)
				
			enddo

		enddo
		
		!This is not used, but can be useful for informative purpose.
! 		lambda_min = lambda_min * (1.0 - dvmax/clight) !max(1.0, lambda_min * (1.0 - dvmax/clight))
! 		lambda_max = lambda_max * (1 + dvmax/clight) * 1.1
! 		write(*,*) "min/max line grid", lambda_min,lambda_max, Nlam

		if (wl_ref > 0) Nlambda_cont = Nlambda_cont + 1
		allocate(cont_waves(Nlambda_cont), stat=alloc_status)
		if (alloc_status>0) then
			write(*,*) "Allocation error cont_waves"
			stop
		endif
		if (wl_ref > 0) then
			lac = 1
			cont_waves(1) = wl_ref
		else
			lac = 0
		endif
		do n=1, Natom
			atom => atoms(n)%ptr_atom
			do kr=1,atom%Ncont
				do la=1, atom%continua(kr)%Nlambda
					lac = lac + 1
					cont_waves(lac) = atom%continua(kr)%lambda(la)
				enddo
				!not used anymore
				deallocate(atom%continua(kr)%lambda)
				!lambda_file (and alpha_file) kept if allocated for explicit continua
			enddo
		enddo
		!sort continuum frequencies
		Nmore_cont_freq = 0.0
		allocate(sorted_indexes(Nlambda_cont),stat=alloc_status)
		if (alloc_status > 0) call error ("Allocation error sorted_indexes (cont)")
		sorted_indexes = bubble_sort(cont_waves)
		cont_waves(:) = cont_waves(sorted_indexes)
		deallocate(sorted_indexes)
		
		!remove duplicates
		allocate(tmp_grid(Nlambda_cont), stat=alloc_status)
		if (alloc_status > 0) call error ("Allocation error tmp_grid (cont)")
		tmp_grid(2:Nlambda_cont) = 0.0
		tmp_grid(1) = cont_waves(1)
		Nremoved = 0
		do la = 2, Nlambda_cont
			if (cont_waves(la) > cont_waves(la-1)) then
				tmp_grid(la) = cont_waves(la)
			else
		 		Nremoved = Nremoved + 1
		 		!!write(*,*) " Removing from the cont grid:", cont_waves(la)
			endif
		enddo

		write(*,*) "Total continuum frequencies, before merging : ", Nlambda_cont - Nremoved!lac
		if (Nremoved > 0) then 
			write(*,*) " ->", Nremoved, " duplicate frequencies"
			deallocate(cont_waves)
			allocate(cont_waves(Nlambda_cont-Nremoved), stat=alloc_status)
			cont_waves(:) = Pack(tmp_grid, tmp_grid > 0)
		endif
		deallocate(tmp_grid)
		max_cont = maxval(cont_waves)
		Nlambda_cont = Nlambda_cont - Nremoved
		
		!Keep, sorted, unique, sole continuum wavelengths
! 		allocate(cont_grid(Nlambda_cont),stat=alloc_status)
! 		if (alloc_status > 0) call error("Allocation error cont_grid")
! 		cont_grid(:) = cont_waves(:)

!-> lines + cont
		lthere_is_lines = (Nlam > 0)
		if (lthere_is_lines) then
			allocate(all_lamin(Nlam), all_lamax(Nlam),stat=alloc_status)
			if (alloc_status > 0) then
				write(*,*) "Allocation error all_lam*"
				stop
			endif
		
		!Store the maximum and minimum extent of each line including max velocity field
		!add the reference wavelength as a line with no width ?? (lambdamin = lambdamax = ref)
			Nlam = 0
			do n=1, Natom
				atom => atoms(n)%ptr_atom

				do kr=1,atom%Nline
    				Nlam = Nlam + 1
					all_lamin(Nlam) = atom%lines(kr)%lambdamin * (1.0 - delta_v/clight)
					all_lamax(Nlam) = atom%lines(kr)%lambdamax * ( 1.0 + delta_v/clight)
				enddo

			enddo
		
			allocate(sorted_indexes(Nlam),stat=alloc_status)
			if (alloc_status > 0) then
				write(*,*) "Allocation error sorted_indexes(Nlam)"
				stop
			endif	
		!sort lines so that all_lamin(1) is always the first line
			sorted_indexes = bubble_sort(all_lamin)
			all_lamin(:) = all_lamin(sorted_indexes)
		!->not necessarily ordered by min to max, but follows the order of lamin
		!so that lmax(1) is associated to lamin(1)  which is important.
		!If lines overlap, the lmax(1) could be associated to lmin(2) for instance.
			all_lamax(:) = all_lamax(sorted_indexes)
			deallocate(sorted_indexes)
		
			group_blue(:) = -1.0
			group_red(:) = -1.0
		
			Ngroup = 1
			group_blue(Ngroup) = all_lamin(1)
			group_red(Ngroup) = all_lamax(1)
		!Find group of lines, and store for each group the lambda_blue and lambda_red of each group
		!if a line overlaps with the previous line, add it to the same group and check the next line.
		!Stop counting lines in a group if the next line does not overlap with the previous line. In
		!the latter case, create a new group and start again.
		! Note: the first and last lines of a group may not overlap. 
			Nline_per_group(:) = 0
			Nline_per_group(1) = 1
			do Nlam = 2, size(all_lamin)
		
				!Is the line overlapping the previous line ? 
				
				!Yes, add it to the same group
					if (((all_lamin(Nlam) >= group_blue(Ngroup)).and.&
					(all_lamin(Nlam) <= group_red(Ngroup))).or.&
					((all_lamax(Nlam) >= group_blue(Ngroup)).and.&
					(all_lamax(Nlam) <= group_red(Ngroup)))) then
				
						group_blue(Ngroup) = min(all_lamin(Nlam), group_blue(Ngroup))
						group_red(Ngroup) = max(all_lamax(Nlam), group_red(Ngroup))
					
						Nline_per_group(Ngroup) = Nline_per_group(Ngroup) + 1
					
				!no, create a new group, starting with this line at first element								 
					else
						Ngroup = Ngroup + 1
						if (Ngroup > MAX_GROUP_OF_LINES) then
							write(*,*) " Error, Ngroup > MAX_GROUP_OF_LINES", Ngroup
							stop
					 	endif
						group_blue(Ngroup) = all_lamin(Nlam)
						group_red(Ngroup) = all_lamax(Nlam)
						Nline_per_group(Ngroup) = 1
					endif		
		
			enddo
! 			write(*,*) Ngroup, Nline_per_group, sum(Nline_per_group)
		
		!write(*,*) " Found ", Ngroup, " groups of lines"
		!write(*,*) " -> ", sum(Nline_per_group), " lines"
			write(*,*) " Found ", sum(Nline_per_group) - Ngroup, " overlapping regions for", sum(Nline_per_group), " lines"
			write(*,*) " --> ", Ngroup, " groups of lines"
			allocate(Nlambda_per_group(Ngroup),corr_hv(Ngroup), stat=alloc_status)
			if (alloc_status > 0) then
				write(*,*) "Allocation error Nlambda_per_group"
				stop
			endif		
		!find number of points to cover lambda red lambda blue
		!!write(*,*) group_blue(1:Ngroup), group_red(1:Ngroup)
		!!stop
			do la=1,Ngroup
				Nlambda_per_group(la) = nint(1 + 1e-3 * clight / hv * 2.0*(group_red(la)-group_blue(la))/(group_red(la)+group_blue(la)))
				l0 = 0.5*(group_red(la)+group_blue(la))
				vmin = 1d-3 * clight*(group_blue(la)-l0)/l0  !km/s
				vmax = 1d-3 * clight*(group_red(la)-l0)/l0
				!equivalent
				Nlambda_per_group(la) = 1 + nint( (vmax-vmin)/hv )!max(1 + nint( (vmax-vmin)/hv ), 2)
! 			write(*,*) Nlambda_per_group(la)
				corr_hv(la) = (vmax-vmin)/(Nlambda_per_group(la)-1) - hv
			enddo

		!Now gather and creates grid for lines
			Nspec_line = sum(Nlambda_per_group)
			allocate(line_waves(Nspec_line), stat=alloc_status)
			if (alloc_status > 0) then
				write(*,*) "Allocation error line_waves"
				stop
			endif
			line_waves = 0.0

			shift = 1
			la = 0
			do n=1, Ngroup	
				line_waves(shift) = group_blue(n)
			!!write(*,*) "start=", n, shift, line_waves(shift)
				do lac=2, Nlambda_per_group(n)
					la = lac + (shift - 1)
 					line_waves(la) = line_waves(la-1) * (1.0 + 1d3 * hv / clight)
! 				write(*,*) n, shift, la, line_waves(la), "gr=", group_red(n)
				enddo
				shift = shift + Nlambda_per_group(n)
			enddo

			deallocate(all_lamin, all_lamax, corr_hv)
		
		!-> Add continuum points beyond last continuum
!In case they are lines beyond the last continuum I add at least3 points per line for the continuum in this region			
!->cont end		this is heavy for nothing but should work !
		!finalise continuum here by reckoning how much freq we need
		
		!make the new points go farther than line_waves for interpolation.
			Nmore_cont_freq = 0
			do n=1, Ngroup

				l0 = group_blue(n)
				l1 = l0 * (1.0 + 1d3 * hv / clight)**(Nlambda_per_group(n))!group_red(n)
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
				allocate(tmp_grid(Nlambda_cont))
				tmp_grid = cont_waves
				deallocate(cont_waves)
				allocate(cont_waves(Nlambda_cont + Nmore_cont_freq))
				cont_waves(1:Nlambda_cont) = tmp_grid(:)
				deallocate(tmp_grid)
				allocate(tmp_grid(Nmore_cont_freq))
				tmp_grid(:) = 0.0_dp
	

				Nmore_cont_freq = 0
				do n=1, Ngroup
! 				write(*,*) "n=",n
					l0 = group_blue(n)
					l1 = l0 * (1.0 + 1d3 * hv / clight)**(Nlambda_per_group(n))!group_red(n)
					if (l0 > max_cont) then
						Nmore_cont_freq = Nmore_cont_freq + 1
						tmp_grid(Nmore_cont_freq) = l0
! 					write(*,*) Nmore_cont_freq , "l0=",l0
					endif
					if (0.5*(l0+l1) > max_cont) then
						Nmore_cont_freq = Nmore_cont_freq + 1
						tmp_grid(Nmore_cont_freq) = 0.5 * (l0+l1)
! 					write(*,*) Nmore_cont_freq , "lmid=", 0.5*(l0+l1)
					endif
					if (l1 > max_cont) then
						Nmore_cont_freq = Nmore_cont_freq + 1
						tmp_grid(Nmore_cont_freq) = l1
! 					write(*,*) Nmore_cont_freq , "l1=",l1
					endif
	
				enddo

				if (Nmore_cont_freq	/= check_new_freq) then
					call Warning("There are probably some frequency missing!")
					write(*,*) "Nmore_freq: ",check_new_freq," Nfreq_added: ", Nmore_cont_freq
				endif				
				allocate(sorted_indexes(Nmore_cont_freq))
				sorted_indexes(:) = bubble_sort(tmp_grid)
				tmp_grid(:) = tmp_grid(sorted_indexes)
				cont_waves(Nlambda_cont+1:Nlambda_cont + Nmore_cont_freq) = tmp_grid(:)
				Nlambda_cont = Nlambda_cont + Nmore_cont_freq
				deallocate(tmp_grid, sorted_indexes)
				allocate(cont_grid(Nlambda_cont),stat=alloc_status)
				if (alloc_status > 0) call error("Allocation error cont_grid")
				cont_grid(:) = cont_waves(:)
			else
				allocate(cont_grid(Nlambda_cont),stat=alloc_status)
				if (alloc_status > 0) call error("Allocation error cont_grid")
				cont_grid(:) = cont_waves(:)
			endif
!->cont end	

			write(*,*) "bounds:"
			write(*,*) "  -> max lam:", maxval(line_waves), " max_cont:", maxval(cont_grid)
		! " max reddest line:", maxval(group_red,mask=group_red>-1),
		else !pure cont
			allocate(cont_grid(Nlambda_cont),stat=alloc_status)
			if (alloc_status > 0) call error("Allocation error cont_grid")
			cont_grid(:) = cont_waves(:)
			write(*,*) "  -> pure cont max_cont:", maxval(cont_grid)
			Nspec_line = 0
		endif !there is lines
		!add lines + continua frequencies
		Nspec_cont = size(cont_waves)
		if (Nspec_cont /= Nlambda_cont) then
		 write(*,*) " Something went wrong with Nlambda cont"
		 stop
		endif
		
		!initiate with lines
		allocate(tmp_grid(Nspec_line+Nspec_cont), stat=alloc_status)
! 		allocate(tmp_grid(Nspec_line), stat=alloc_status)
		if (alloc_status > 0) call error ("Allocation error tmp_grid")
		tmp_grid(:) = -99
		do la=1,Nspec_line
			tmp_grid(la) = line_waves(la)
		enddo

		Nwaves = Nspec_line
		!add continuum wavlengths (including reference wavelength), only outside line groups		
! 		!First values below or beyond first and last groups
		la = 0
		do lac=Nspec_line+1, Nspec_cont+Nspec_line
			if ((cont_waves(lac-Nspec_line) < group_blue(1)) .or. (cont_waves(lac-Nspec_line) > group_red(Ngroup))) then
				tmp_grid(lac) = cont_waves(lac-Nspec_line)
				Nwaves = Nwaves + 1
				if (cont_waves(lac-Nspec_line) < group_blue(1)) la = lac
			endif
		enddo

		!now values between groups
		do lac=la+1, Nspec_cont+Nspec_line
			group_loop : do n=2, Ngroup
				if ((cont_waves(lac-Nspec_line) > group_red(n-1)).and.(cont_waves(lac-Nspec_line) < group_blue(n))) then
					Nwaves = Nwaves + 1
					tmp_grid(lac) = cont_waves(lac-Nspec_line)
				!else
					! be smart and cycle to accelerate
				endif
			enddo group_loop
		enddo
		!write(*,*) maxval(cont_grid), minval(cont_grid)
		!write(*,*) maxval(tmp_grid), minval(tmp_grid,mask=tmp_grid /= -99)
		
		!!Nwaves = size(pack(tmp_grid, tmp_grid > 0))
		if (lthere_is_lines) deallocate(Nlambda_per_group, line_waves)
		deallocate(cont_waves)

		
		!continuum frequencies are sorted and so are the line frequencies
		!but they are added at the end, so sorted is needed, but I can improve the previous
		!loop to fill the tmp_frid in the ascending order of wavelengths
		allocate(outgrid(Nwaves),stat=alloc_status)
		tmp_grid = tmp_grid(bubble_sort(tmp_grid))
		outgrid(:) = -99.0 !check
		outgrid = pack(tmp_grid, mask=tmp_grid > 0)
		!write(*,*) maxval(outgrid), minval(outgrid)
		!write(*,*) "Nwaves=",Nwaves," la=",1, " lambda=",outgrid(1)
		do lac=2, Nwaves
			!write(*,*) "la=",lac, " lambda=",outgrid(lac)
			if (outgrid(lac) <= outgrid(lac-1)) then
				write(*,*) lac, "lambda = ", outgrid(lac), outgrid(lac-1), minval(outgrid)
				call error("Sorted problem")
			endif
		enddo
! 		write(*,*) outgrid(1), minval(group_blue,mask=group_blue >= 0)
! 		write(*,*) outgrid(Nwaves), maxval(group_red,mask=group_red >= 0), maxval(cont_grid)
! stop
		deallocate(tmp_grid)

		write(*,*) Nwaves, " unique wavelengths" !they are no eliminated lines
		write(*,*) Nspec_line, " line wavelengths"
		write(*,*) Nspec_cont, " continuum wavelengths"
		if (lthere_is_lines) then
			write(*,*) "Mean number of lines per group:", real(sum(Nline_per_group))/real(Ngroup)
			write(*,*) "Mean number of wavelengths per group:", real(Nspec_line)/real(Ngroup)
			write(*,*) "Mean number of wavelengths per line:", real(Nspec_line)/real(Ntrans-Ncont)
			write(*,*) "Resolution of line's groups (km/s):", hv
		endif
		write(*,*) "Mean number of wavelengths per continuum:", real(Nwaves - Nspec_line) / real(Ncont)

		lam_unit = "nm"
		l0 = minval(outgrid); l1 = maxval(outgrid)
		if (l1 > 1500.) then
			l1 = l1 *1e-4
			lam_unit = "mum"
		else if (l1 > 1e6) then
			l1 = 10000000./l1
			lam_unit = "cm^-1"
		else if (l1 > 1e6) then
			l1 = l1 * 1e-9 * 1e3
			lam_unit = "mm"
		else if (l1 > 1e7) then
			l1 = l1 * 1e-9 * 1e2
			lam_unit = "cm"
		endif
		write(*,*) "Wavelength grid:", l0," nm", l1,lam_unit  
			
		!allocate indexes on the grid
		Nlambda_max_line = 0
		Nlambda_max_cont = 0
		Nlambda_max_trans = 0
		do n=1,Natom
			atom => Atoms(n)%ptr_atom
			do kr=1,atom%Ncont !only on the cont grid !

				atom%continua(kr)%Nb = locate(outgrid, atom%continua(kr)%lambdamin)
				atom%continua(kr)%Nr = locate(outgrid, atom%continua(kr)%lambdamax)

				atom%continua(kr)%Nblue = locate(cont_grid, atom%continua(kr)%lambdamin)
				atom%continua(kr)%Nred = locate(cont_grid, atom%continua(kr)%lambdamax)
				atom%continua(kr)%Nmid = locate(cont_grid, 0.5*(atom%continua(kr)%lambdamin+atom%continua(kr)%lambdamax))
				atom%continua(kr)%N0 = locate(cont_grid, atom%continua(kr)%lambda0)
				atom%continua(kr)%Nlambda = atom%continua(kr)%Nred - atom%continua(kr)%Nblue + 1
				
				!in any problem of grid resolution etc or locate approximation.
				!We take Nred-1 to be sure than the lambda_cont(Nred) <= lambda0.
				!Only if not dissolution.
				!We just need to avoind having cont_lambda(Nred)>lambda0, since the cross section is in (lambda/lambda0)**3
 				if (cont_grid(atom%continua(kr)%N0) /= atom%continua(kr)%lambda0) then
					!!write(*,*) " Beware, cont%lambda0 might not be on the grid", kr, atom%continua(kr)%i, atom%continua(kr)%j
					if (.not.ldissolve) then
						if (cont_grid(atom%continua(kr)%Nred) > atom%continua(kr)%lambda0) then
							call Warning("continuum Nred larger than lambda0 !")
							write(*,*) atom%continua(kr)%lambda0, cont_grid(atom%continua(kr)%Nred)
							if (cont_grid(atom%continua(kr)%Nred-1) <= atom%continua(kr)%lambda0) then
								write(*,*) " ++++ adjusting Nred", " lambda0=",atom%continua(kr)%lambda0
								!To do while until <= lambda0
								atom%continua(kr)%Nred = atom%continua(kr)%Nred-1
								atom%continua(kr)%N0 = atom%continua(kr)%Nred
								atom%continua(kr)%Nlambda = atom%continua(kr)%Nred - atom%continua(kr)%Nblue + 1
								atom%continua(kr)%Nr = locate(outgrid,cont_grid(atom%continua(kr)%Nred))
								write(*,*) "new val at Nred:", outgrid(atom%continua(kr)%Nr), cont_grid(atom%continua(kr)%Nred)
							endif
						endif
					endif
				endif
				!sur la grille totale
				Nlambda_max_cont = max(Nlambda_max_cont,atom%continua(kr)%Nr-atom%continua(kr)%Nb+1) 
				
			enddo
			
			do kr=1,atom%Nline
				atom%lines(kr)%Nblue = locate(outgrid, atom%lines(kr)%lambdamin)
				atom%lines(kr)%Nred = locate(outgrid, atom%lines(kr)%lambdamax)
				atom%lines(kr)%Nmid = locate(outgrid, atom%lines(kr)%lambda0)
				atom%lines(kr)%Nlambda = atom%lines(kr)%Nred - atom%lines(kr)%Nblue + 1
				if (atom%lines(kr)%Nblue-nint(1e-3 * dvmax/hv) < 1) then
					write(*,*) "Error for line ", kr, " of atom ", atom%ID
					write(*,*) " Nblue below 1"
					write(*,*) "Nb=",atom%lines(kr)%Nblue, " dvmin=",-nint(1e-3 * dvmax/hv)
					write(*,*) " -> sum :", atom%lines(kr)%Nred-nint(1e-3 * dvmax/hv)
					stop
				else if(atom%lines(kr)%Nred+nint(1e-3 * dvmax/hv) > Nwaves) then
					write(*,*) "Error for line ", kr, " of atom ", atom%ID
					write(*,*) " Nred larger than Nwaves", Nwaves, size(outgrid)
					write(*,*) "Nr=",atom%lines(kr)%Nred, " dvmax=",nint(1e-3 * dvmax/hv)
					write(*,*) " -> sum :", atom%lines(kr)%Nred+nint(1e-3 * dvmax/hv)
					stop
				endif
! 				write(*,*) "line", kr, " lam0=",atom%lines(kr)%lambda0, atom%lines(kr)%lambdamin, atom%lines(kr)%lambdamax
! 				write(*,*) " -> bounds on the grid:", outgrid(atom%lines(kr)%Nblue), outgrid(atom%lines(kr)%Nred)
! 				write(*,*) " -> max extent:", outgrid(atom%lines(kr)%Nblue)*(1.0 - delta_v/clight), outgrid(atom%lines(kr)%Nred)*(1.0 + delta_v/clight)
				
				!does not take the shift into account
				Nlambda_max_line = max(Nlambda_max_line, atom%lines(kr)%Nlambda)
			enddo

			atom => NULL()
		enddo
		write(*,*) "Number of max freq points for all lines at this resolution :", Nlambda_max_line 
		write(*,*) "Number of max freq points for all cont at this resolution :", Nlambda_max_cont
		!takes the shift into account
		Nlambda_max_trans = max(Nlambda_max_line+2*int( sign(1.0_dp, delta_v) * ( 1e-3 * abs(delta_v) / hv + 0.5 ) ),Nlambda_max_cont)
		write(*,*) "Number of max freq points for all trans at this resolution :", Nlambda_max_trans

	return
	end subroutine make_wavelength_grid_new

  end module getlambda
