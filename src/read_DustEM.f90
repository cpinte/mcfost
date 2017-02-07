module read_DustEM

  use grains


  implicit none

  character(len=512) :: DustEM_dir

contains

  subroutine get_DustEM_dim(pop)

    integer, intent(in) :: pop

    character(len=512) :: filename, the_char
    real :: material_density

    call get_environment_variable('DUSTEM_DIR',DustEM_dir)
    if (DustEM_dir == "") then
       write(*,*) "ERROR: environnement variable DUSTEM_DIR is not defined."
       write(*,*) "Exiting."
       stop
    endif

    ! 1) open lambda file
    filename = trim(DustEM_dir)//"/oprop/LAMBDA.DAT"
    open(unit=1, file=filename,status='old')
    the_char = read_com(1)
    read(the_char, fmt=*) op_file_n_lambda(pop)
    close(unit=1)

    ! 2) open opacity file : DustEM/oprop/Q_<NAME>.DAT
    filename = trim(DustEM_dir)//"/oprop/Q_"//trim(dust_pop(pop)%indices(1))//".DAT"
    open(unit=1, file=filename,status='old')
    the_char = read_com(1)
    read(the_char, fmt=*) op_file_na(pop)
    close(unit=1)

    ! Material densities are parameters provided in data directory in DustEM
    ! for simplicity, they are hard-coded here
    select case(trim(dust_pop(pop)%indices(1)))
    case('BG_DBP90')
       material_density = 3.0
    case('Gra')
       material_density = 2.24
    case('PAH0')
       material_density = 2.24
    case('PAH0_DBP90')
       material_density = 2.25
    case('PAH0_DL01')
       material_density = 2.24
    case('PAH0_DL07')
       material_density = 2.24
    case('PAH1')
       material_density = 2.24
    case('PAH1_DBP90')
       material_density = 2.25
    case('PAH1_DL01')
       material_density = 2.24
    case('PAH1_DL07')
       material_density = 2.24
    case('VSG_DBP90')
       material_density = 2.30
    case('aSil')
       material_density = 3.5
    case('amCBE')
       material_density = 1.81
    case('amCBEx')
       material_density = 1.81
    case default
       write(*,*) "Material density is unknown for the DustEM material: "//trim(dust_pop(pop)%indices(1))
       write(*,*) "Exiting"
       stop
    end select
    dust_pop(pop)%component_rho1g(1) = material_density
    dust_pop(pop)%rho1g_avg = material_density

    return

  end subroutine get_DustEM_dim

  !*********************************************************************************

  subroutine read_DustEM_cross_sections(pop)

    integer, intent(in) :: pop

    character(len=512) :: filename, the_char
    integer :: i,j, na

    ! opacity file : DustEM/oprop/C_<NAME>.DAT
    ! g file : DustEM/oprop/G_<NAME>.DAT
    ! heating capacity file : DustEM/hcap/C_<NAME>.DAT

    !1) reading lambda grid for all Q files (LAMBDA.DAT)
    filename = trim(DustEM_dir)//"/oprop/LAMBDA.DAT"
    open(unit=1, file=filename,status='old')
    the_char = read_com(1)
    read(the_char, fmt=*) op_file_n_lambda(pop)
    do j=1,op_file_n_lambda(pop)
       read(1,*) op_file_lambda(j,pop)
    enddo
    close(unit=1)

    !2) reading Qabs, Qsca (Q_*.DAT) and G_*.DAT
    filename = trim(DustEM_dir)//"/oprop/Q_"//trim(dust_pop(pop)%indices(1))//".DAT"
    open(unit=1, file=filename,status='old')
    the_char = read_com(1)
    read(the_char, fmt=*) op_file_na(pop)

    ! get grain sizes
    read(1,*) (op_file_log_r_grain(i,pop), i=1, op_file_na(pop))
    ! convert in log
    do i=1, op_file_na(pop)
       op_file_log_r_grain(i,pop) = log(op_file_log_r_grain(i,pop))
    enddo

    ! get Qabs
    the_char = read_com(1)
    backspace(1)
    do j=1,op_file_n_lambda(pop)
       read(1,*) (op_file_Qext(j,i,pop), i=1, op_file_na(pop)) ! this is actually Qabs
    enddo

    ! get Qsca
    the_char = read_com(1)
    backspace(1)
    do j=1,op_file_n_lambda(pop)
       read(1,*) (op_file_Qsca(j,i,pop), i=1, op_file_na(pop))
    enddo

    ! Qext = Qabs + Qsca
    op_file_Qext(:,:,pop) = op_file_Qext(:,:,pop) + op_file_Qsca(:,:,pop)
    close(unit=1)

    !3) read g
    filename = trim(DustEM_dir)//"/oprop/G_"//trim(dust_pop(pop)%indices(1))//".DAT"
    open(unit=1, file=filename,status='old')
    the_char = read_com(1)
    read(the_char, fmt=*) na
    if (na /= op_file_na(pop)) then
       write(*,*) "DustEM: G file must have same dimension as opacity file"
       write(*,*) "Exiting"
       stop
    endif

    ! get grain sizes
    read(1,*) (op_file_log_r_grain(i,pop), i=1, op_file_na(pop))
    ! convert in log
    do i=1, op_file_na(pop)
       op_file_log_r_grain(i,pop) = log(op_file_log_r_grain(i,pop))
    enddo

    ! get g
    the_char = read_com(1)
    backspace(1)
    do j=1,op_file_n_lambda(pop)
       read(1,*) (op_file_g(j,i,pop), i=1, op_file_na(pop))
    enddo
    close(1)

    return

  end subroutine read_DustEM_cross_sections

  !*********************************************************************************

   subroutine get_DustEM_specific_heat_dim(pop)

    use utils, only : is_file

    integer, intent(in) :: pop

    character(len=512) :: filename, the_char
    integer :: na

    ! heating capacity file : DustEM/hcap/C_<NAME>.DAT
    filename = trim(DustEM_dir)//"/hcap/C_"//trim(dust_pop(pop)%indices(1))//".DAT"
    if (is_file(filename)) then
       open(unit=1, file=filename,status='old')
       the_char = read_com(1)
       read(the_char, fmt=*) na
       read(1,*) ! skip 1 line
       read(1,*) file_sh_nT(pop)
       close(unit=1)
    else
       write(*,*) "ERROR: "//trim(filename)//" does not exist"
       write(*,*) "Are you sure you wish to perform stochastic heating ?"
       write(*,*) "Exiting"
       stop
    endif

    return

  end subroutine get_DustEM_specific_heat_dim

  !*********************************************************************************

  subroutine read_DustEM_specific_heat(pop)

    integer, intent(in) :: pop

    character(len=512) :: filename, the_char
    real :: C
    integer :: k, na

    filename = trim(DustEM_dir)//"/hcap/C_"//trim(dust_pop(pop)%indices(1))//".DAT"

    write(*,*) "Reading specific heat file: "//trim(filename)

    open(unit=1, file=filename,status='old')
    the_char = read_com(1)
    read(the_char, fmt=*) na
    if (na /= op_file_na(pop)) then
       write(*,*) "DustEM: heat capacity file must have same dimension as opacity file"
       write(*,*) "Exiting"
       stop
    endif
    read(1,*) ! skip 1 line : I already have the size
    read(1,*) file_sh_nT(pop) ! already done but ok

    do k=1, file_sh_nT(pop)
       read(1,*) file_sh_T(k,pop), C
       ! convert log T to Temperature
       file_sh_T(k,pop) = 10**file_sh_T(k,pop)
       ! convert from log and to specific heat capacity (ie volumic to massic)
       file_sh(k,pop) = 10**C / dust_pop(pop)%rho1g_avg ! erg/K/g
    enddo ! k

    close(unit=1)

    return

  end subroutine read_DustEM_specific_heat

  !*********************************************************************************

  function read_com(file_input) result(line)
    ! reads comment lines beginning with #
    ! returns 1st line that is not a comment as a string
    ! Adapted from DustEM

    integer, intent(in)    :: file_input
    logical                :: comment_found
    character (len=512)    :: line

    comment_found = .true.
    do while(comment_found)
       read(file_input,*) line
       if (line(1:1) /= '#') then
          comment_found = .false.
          return
       endif
    enddo

    return

  end function read_com

end module read_DustEM
