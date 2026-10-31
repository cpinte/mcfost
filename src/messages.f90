module messages

  implicit none

  public :: warning, error

  private

contains

  subroutine warning(msg,msg2,msg3)

    character(len=*), intent(in)           :: msg
    character(len=*), intent(in), optional :: msg2, msg3

    ! explicit format: list-directed output is wrapped at 80 columns by ifort
    write(*,'(1x,a)') "WARNING: "//trim(msg)
    if (present(msg2)) write(*,'(1x,a)') trim(msg2)
    if (present(msg3)) write(*,'(1x,a)') trim(msg3)

    return

  end subroutine warning

  !---------------------------------

  subroutine error(msg,msg2,msg3,ierr)

    character(len=*), intent(in)           :: msg
    character(len=*), intent(in), optional :: msg2, msg3
    integer, intent(out), optional         :: ierr

    ! explicit format: list-directed output is wrapped at 80 columns by ifort
    write(*,'(1x,a)') "ERROR: "//trim(msg)
    if (present(msg2)) write(*,'(1x,a)') trim(msg2)
    if (present(msg3)) write(*,'(1x,a)') trim(msg3)

    if (present(ierr)) then
       ierr = 1
       return
    else
       write(*,*) "Exiting."
       call exit(1)
    endif

  end subroutine error

end module messages
