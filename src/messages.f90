module messages

  implicit none

  public :: warning, error

  private

contains

  subroutine warning(msg,msg2,msg3)

    character(len=*), intent(in)           :: msg
    character(len=*), intent(in), optional :: msg2, msg3

    write(*,*) "WARNING: "//trim(msg)
    if (present(msg2)) write(*,*) trim(msg2)
    if (present(msg3)) write(*,*) trim(msg3)

    return

  end subroutine warning

  !---------------------------------

  subroutine error(msg,msg2,msg3)

    character(len=*), intent(in)           :: msg
    character(len=*), intent(in), optional :: msg2, msg3

    write(*,*) "ERROR: "//trim(msg)
    if (present(msg2)) write(*,*) trim(msg2)
    if (present(msg3)) write(*,*) trim(msg3)
    write(*,*) "Exiting."
    call exit(1)

  end subroutine error

end module messages
