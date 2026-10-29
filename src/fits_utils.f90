module fits_utils

  use messages, only : error

  implicit none

  public :: cfitsWrite, print_error, fits_overwrite

  private

contains

subroutine cfitsWrite(filename,tab,dim)

  implicit none

  character(len=*), intent(in) :: filename
  real, dimension(*), intent(in) :: tab
  integer, dimension(:), intent(in) :: dim ! dim == shape(tab)

  integer :: status,unit,blocksize,bitpix,naxis
  integer :: group,fpixel,nelements
  logical :: simple, extend

  !  Get an unused Logical Unit Number to use to open the FITS file.
  status=0
  call ftgiou (unit,status)

  !  Create the new empty FITS file.
  blocksize=1
  call ftinit(unit,fits_overwrite(filename),blocksize,status)

  simple=.true.
  ! le signe - signifie que l'on ecrit des reels dans le fits
  bitpix=-32
  extend=.true.


  !  Write the required header keywords.
  naxis=size(dim)
  call ftphpr(unit,simple,bitpix,naxis,dim,0,1,extend,status)

  !  Write the array to the FITS file.
  group=1
  fpixel=1
  nelements=product(dim)

  call ftppre(unit,group,fpixel,nelements,tab,status)

  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages and stop
  if (status > 0) call print_error(status)

  return

end subroutine cfitsWrite

!***************************************************


subroutine print_error(status)
  ! PRINT_ERROR prints out the FITSIO error messages to the user.

  integer status
  character ( len = 30 ) errtext
  character ( len = 80 ) errmessage

  !  Check if status is OK (no error); if so, simply return.
  if (status <= 0) then
     return
  end if

  !  Get the text string which describes the error
  call ftgerr(status,errtext)
  print *,'FITSIO Error Status =',status,': ',errtext

  !  Read and print out all the error messages on the FITSIO stack
  call ftgmsg(errmessage)
  do while (errmessage .ne. ' ')
     print *,errmessage
     call ftgmsg(errmessage)
  end do

  ! A failed read or write would leave missing, incomplete or stale outputs
  call error("FITS file input/output failed (see the FITSIO messages above)")

  return
end subroutine print_error

!***************************************************

function fits_overwrite(filename) result(name)
  ! FITSIO name to create filename, replacing an existing file ("!" prefix).
  ! Outputs from a previous run are protected by the data directory backup
  ! (data_*_old) and -no_backup: within a run, files are rewritten (e.g. iterations)

  character(len=*), intent(in) :: filename
  character(len=:), allocatable :: name

  if (index(filename, "!") == 1) then
     name = trim(filename)
  else
     name = "!"//trim(filename)
  endif

end function fits_overwrite

end module fits_utils
