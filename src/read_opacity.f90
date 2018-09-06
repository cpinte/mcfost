module read_opacity

  use parametres
  use grains
  use wavelengths
  use mcfost_env
  use utils
  use read_DustEM


  implicit none

  character(len=512), dimension(:), allocatable :: sh_file

  contains

subroutine read_opacity_file(pop)

  integer, intent(in) :: pop
  character(len=512) :: filename

  integer :: l
  type(dust_pop_type), pointer :: dp

  dp => dust_pop(pop)

  ! We look for the directory and check if the file exists in  get_opacity_file_dim
  ! We do not need to redo it here
  filename = trim(dp%indices(1))
  write(*,*) "Reading opacity file "//trim(filename) ;

  ! load optical data from file
  if (dp%is_Misselt_opacity_file) then
     call misselt_load(filename, pop)
  else if (dp%is_DustEM_opacity_file) then
     call read_DustEM_cross_sections(pop)
  else
     dp%component_rho1g(1) = 2.5
     dp%rho1g_avg = 2.5
     call draine_load(filename, pop)
     !call draine_load(filename, op_file_n_lambda, op_file_na, 10, 1, &
    !      tmp_lambda, tmp_r_grain,  tmp_Q_ext, tmp_Q_abs, tmp_Q_sca, tmp_g, 4)
  endif

  ! abs car le fichier de lambda pour les PAHs est a l'envers
  op_file_delta_lambda(1,pop) = abs(op_file_lambda(2,pop) - op_file_lambda(1,pop))
  if (n_lambda > 1) then
     op_file_delta_lambda(n_lambda,pop) = abs(op_file_lambda(n_lambda,pop) - op_file_lambda(n_lambda-1,pop))
  endif
  do l=2,op_file_n_lambda(pop)-1
     op_file_delta_lambda(l,pop) = 0.5* abs(op_file_lambda(l+1,pop) - op_file_lambda(l-1,pop))
  enddo

  return

end subroutine read_opacity_file

!*************************************************************

subroutine get_opacity_file_dim(pop)

  integer, intent(in) :: pop
  character(len=512) :: filename, dir
  type(dust_pop_type), pointer :: dp

  integer :: ios

  dp => dust_pop(pop)


  filename = trim(dp%indices(1))

  if (.not.dp%is_DustEM_opacity_file) then
     dir = in_dir(filename, dust_dir,  status=ios)
     if (ios /=0) call error("dust file cannot be found:"//trim(filename))
     filename = trim(dir)//trim(filename) ;

     ! Updating indice filename with directory
     dp%indices(1) = filename ;

     if (dp%is_Misselt_opacity_file) then
        call get_misselt_dim(filename, pop)
     else
        dp%component_rho1g(1) = 2.5
        dp%rho1g_avg = 2.5
        call get_draine_dim(filename, pop)
     endif

  else ! DustEM file
     call get_DustEM_dim(pop)
  endif

  return

end subroutine get_opacity_file_dim

!*************************************************************

subroutine alloc_mem_opacity_file()

  integer :: alloc_status, n_lambda, na, nT, pop

  na = maxval(op_file_na(:))
  n_lambda = maxval(op_file_n_lambda(:))

  allocate(op_file_Qext(n_lambda,na,n_pop), op_file_Qsca(n_lambda,na,n_pop), &
       op_file_g(n_lambda,na,n_pop), op_file_log_r_grain(na,n_pop), &
       op_file_lambda(n_lambda,n_pop), op_file_delta_lambda(n_lambda,n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error opacity file')

  if (lnRE) then
     do pop = 1, n_pop
        if (dust_pop(pop)%methode_chauffage == 3) then
           if (dust_pop(pop)%is_Misselt_opacity_file) call get_Misselt_specific_heat_dim(pop)
           if (dust_pop(pop)%is_DustEM_opacity_file) call get_DustEM_specific_heat_dim(pop)
        endif
     enddo

     nT = maxval(file_sh_nT(:))
     allocate(file_sh_T(nT,n_pop), file_sh(nT,n_pop))
     if (alloc_status > 0) call error('Allocation error specific heat capacity file')
  endif

  return

end subroutine alloc_mem_opacity_file

!*************************************************************

subroutine free_mem_opacity_file()
  ! Pas utilise pour le moment

  deallocate(op_file_Qext, op_file_Qsca, op_file_g, op_file_log_r_grain, &
       op_file_lambda, op_file_delta_lambda)

  if (allocated(file_sh_T)) deallocate(file_sh_T, file_sh)

  return

end subroutine free_mem_opacity_file

!*************************************************************

subroutine draine_load(file,pop)

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7, nh2 = 1, ns=1

  real :: fbuffer, a
  integer :: i,j,l

  open(unit=1,file=file)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) op_file_na(pop)
  read(1,*) op_file_n_lambda(pop)

  do l=1, nh2
     read(1,*)
  enddo

  ! lecture
  do i=1, op_file_na(pop)
     read(1,*) a ; op_file_log_r_grain(i,pop) = log(a)
     read(1,*)
     do j=1, op_file_n_lambda(pop)
        read(1,*) op_file_lambda(j,pop), op_file_Qext(j,i,pop), fbuffer, op_file_Qsca(j,i,pop), op_file_g(j,i,pop)
     enddo
     if (i < op_file_na(pop)) then
        ! skip ns ligne
        do l=1,ns
           read(1,*)
        enddo
     endif
  enddo
  close(unit=1)

  return

end subroutine draine_load

!*************************************************************

subroutine get_draine_dim(file,pop)

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7
  integer :: l

  open(unit=1,file=file)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) op_file_na(pop)
  read(1,*) op_file_n_lambda(pop)

  close(unit=1)

  return

end subroutine get_draine_dim

!*************************************************************

subroutine misselt_load(file, pop)
  ! Allocate and fills the arrays op_file_XXXX
  ! C. Pinte
  ! 27/06/2014

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7, nh2 = 2, ns=1

  character(len=5) :: sbuffer

  real :: fbuffer, fbuffer2, a, material_density
  integer :: i,j,l, ibuffer

  open(unit=1,file=file)
  read(1,*)  sh_file(pop)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) sbuffer, op_file_na(pop)
  read(1,*) sbuffer, op_file_n_lambda(pop)
  read(1,*) ibuffer, material_density
  dust_pop(pop)%component_rho1g(1) = material_density
  dust_pop(pop)%rho1g_avg = material_density

  do l=1, nh2
     read(1,*)
  enddo

  ! lecture
  do i=1,op_file_na(pop)
     read(1,*) a ; op_file_log_r_grain(i,pop) = log(a)
     read(1,*)
     do j=1,op_file_n_lambda(pop)
        read(1,*) fbuffer, op_file_lambda(j,pop), fbuffer2, op_file_Qsca(j,i,pop), op_file_Qext(j,i,pop), op_file_g(j,i,pop)
     enddo
     if (i < op_file_na(pop)) then
        ! skip ns ligne
        do l=1,ns
           read(1,*)
        enddo
     endif
  enddo
  close(unit=1)

  return

end subroutine misselt_load

!*************************************************************

subroutine get_misselt_dim(file, pop)

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7

  character(len=5) :: sbuffer
  real :: material_density
  integer :: l, ibuffer

  open(unit=1,file=file)
  read(1,*)  sh_file(pop)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) sbuffer, op_file_na(pop)
  read(1,*) sbuffer, op_file_n_lambda(pop)
  read(1,*) ibuffer, material_density
  dust_pop(pop)%component_rho1g(1) = material_density
  dust_pop(pop)%rho1g_avg = material_density

  close(unit=1)

  return

end subroutine get_misselt_dim

!*************************************************************

subroutine read_Misselt_specific_heat(pop)

  integer, intent(in) :: pop

  integer :: ios, status, n_comment, k
  real :: fbuffer

  open(unit=1, file=sh_file(pop), status='old', iostat=ios)
  if (ios/=0) call error("cannot open SH file "//trim(sh_file(pop)))
  write(*,*) "Reading specific heat capacity: "//trim(sh_file(pop))

  ! On elimine les lignes avec des commentaires
  status = 1
  n_comment = 0
  do while (status /= 0)
     n_comment = n_comment + 1
     read(1,*,iostat=status) fbuffer
  enddo
  n_comment = n_comment - 1

  ! Lecture proprement dite
  rewind(1)
  ! On passe les commentaires
  do k=1, n_comment
     read(1,*)
  enddo

  ! Lecture specific heat capacity
  do k=1,file_sh_nT(pop)
     read(1,*) file_sh_T(k,pop), fbuffer, file_sh(k,pop)
  enddo

  return

end subroutine read_Misselt_specific_heat

!******************************************************

subroutine get_Misselt_specific_heat_dim(pop)

  integer, intent(in) :: pop

  integer :: ios, status, n_comment, n_Temperatures
  real :: fbuffer

  open(unit=1, file=sh_file(pop), status='old', iostat=ios)
  if (ios/=0) call error("cannot open SH file "//trim(sh_file(pop)))
  write(*,*) "Reading specific heat capacity: "//trim(sh_file(pop))

  ! On elimine les lignes avec des commentaires
  status = 1
  n_comment = 0
  do while (status /= 0)
     n_comment = n_comment + 1
     read(1,*,iostat=status) fbuffer
  enddo
  n_comment = n_comment - 1

  ! On compte les lignes avec des donnees
  status=0
  n_Temperatures=1 ! On a deja lu une ligne en cherchant les commentaires
  do while(status==0)
     n_Temperatures=n_Temperatures+1
     read(1,*,iostat=status)
  enddo
  n_Temperatures = n_Temperatures - 1

  file_sh_nT(pop) = n_Temperatures

  return

end subroutine get_Misselt_specific_heat_dim

!*************************************************************


end module read_opacity
