MODULE uplow

  ! Function to change from upper/lower case string to
  ! lower/upper case string.

  IMPLICIT NONE
  character(len=26), parameter, private :: low  = "abcdefghijklmnopqrstuvwxyz"
  character(len=26), parameter, private :: high = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
contains

  function to_upper(s) result(t)
    ! returns upper case of s
    implicit none
    character(len=*), intent(in) :: s
    character(len=len(s))        :: t

    character(len=1), save       :: convtable(0:255)
    logical, save                :: first = .true.
    integer                      :: i

    if(first) then
       do i=0,255
          convtable(i) = char(i)
       enddo
       do i=1,len(low)
          convtable(iachar(low(i:i))) = char(iachar(high(i:i)))
       enddo
       first = .false.
    endif

    t = s

    do i=1,len_trim(s)
       t(i:i) = convtable(iachar(s(i:i)))
    enddo

  end function to_upper

  function to_lower(s) result(t)
    ! returns lower case of s
    implicit none
    character(len=*), intent(in) :: s
    character(len=len(s))        :: t

    character(len=1), save :: convtable(0:255)
    logical, save          :: first = .true.
    integer                :: i

    if(first) then
       do i=0,255
          convtable(i) = char(i)
       enddo
       do i = 1,len(low)
          convtable(iachar(high(i:i))) = char(iachar(low(i:i)))
       enddo
       first = .false.
    endif

    t = s

    do i=1,len_trim(s)
       t(i:i) = convtable(iachar(s(i:i)))
    enddo

  end function to_lower



  FUNCTION strcompress( input_string, n ) RESULT ( output_string )

    ! -- Arguments
    CHARACTER( * ), INTENT( IN )  :: input_string
    INTEGER,        INTENT( OUT ) :: n

    ! -- Function result
    CHARACTER( LEN( input_string ) ) :: output_string

    ! -- Local parameters
    INTEGER,        PARAMETER :: IACHAR_SPACE = 32, &
         IACHAR_TAB   = 9

    ! -- Local variables
    INTEGER :: i
    INTEGER :: iachar_character

    ! -- Initialise output string
    output_string = ' '

    ! -- Initialise output string "useful" length counter
    n = 0

    ! -- Loop over string elements
    DO i = 1, LEN( input_string )

       ! -- Convert the current character to its position
       ! -- in the ASCII collating sequence
       iachar_character = IACHAR( input_string( i:i ) )

       ! -- If the character is NOT a space ' ' or a tab '->|'
       ! -- copy it to the output string.
       IF ( iachar_character /= IACHAR_SPACE .AND. &
            iachar_character /= IACHAR_TAB         ) THEN
          n = n + 1
          output_string( n:n ) = input_string( i:i )
       END IF

    END DO

  END FUNCTION strcompress


END MODULE uplow
