MODULE getline

  !Routine to read formatted input files

  IMPLICIT NONE
  integer, parameter :: MAX_LENGTH=512 !maximum character read on a line
  integer, parameter :: MAX_KEYWORD_SIZE = 15

  CONTAINS

  SUBROUTINE getnextLine(unit, commentChar, FMT, line, Nread)

  !Read next line which is not a comment line nor an empty line
  !return that line and the len of the line Nread

  character(len=MAX_LENGTH), intent(out) :: line
  character(len=5) :: trline
  character(len=1), intent(in) :: commentChar
  integer, intent(out) :: Nread
  integer, intent(in) :: unit
  integer :: EOF, n, nMax, count
  character(len=*), intent(in) :: FMT

  Nread = 0
  n = 0
  EOF = 0
  nMax = 30 !if error stop at n = 15
  count = 0
  !write(*,*) "formatline used=", FMT
  do while (n.eq.0 .and. count.lt.nMax)
   read(unit, FMT, IOSTAT=EOF) line !'(512A)'
   !!write(*,*) n, "count=",count, line
   if (EOF.ne.0) then
    !write(*,*) "EOF error, skipping"
    !stop
   else if (len(trim(line)).eq.0) then
    !write(*,*) "Empty line, skipping"
    !stop
   else if (line(1:1).eq.commentChar) then
    !write(*,*) "Commented line, skipping"
    !stop
   else
    n = 1
    Nread = len(line)!len(trim(line))
   end if
   count = count + 1
  end do

  if (n.eq.0) then
    write(*,*) "Error when reading line",line, count
    write(*,*) "exiting..."
    stop
  end if

  RETURN
  END SUBROUTINE getnextLine

  SUBROUTINE KompressIndex (str, i, j)
   character(len=*), intent(in) :: str
   integer, intent(out) :: i, j
   integer :: N, k
   N = len(str)

   do k=1,N
    if (str(k:k).ne." ") then
     i = k
     exit
    end if
   end do
   do k=N,1,-1
    if (str(k:k).ne." ") then
     j = k
     exit
    end if
   end do

  RETURN
  END SUBROUTINE KompressIndex


  END MODULE getline
