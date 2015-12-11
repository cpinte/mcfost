!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2015 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://users.monash.edu.au/~dprice/phantom                               !
!--------------------------------------------------------------------------!
!+
!  MODULE: dump_utils
!
!  DESCRIPTION: Utility routines used when reading and writing the
!   sphNG/Phantom dump file format
!
!  REFERENCES: None
!
!  OWNER: Daniel Price
!
!  $Id: 46275e053845ccce42b676028be688c5d01b851c $
!
!  RUNTIME PARAMETERS: None
!
!  DEPENDENCIES: None
!+
!--------------------------------------------------------------------------
module dump_utils
 implicit none
 public :: open_dumpfile_w, open_dumpfile_r, get_error_text
 public :: tag,check_tag
 public :: skipblock,read_blockheader,skip_headerblock
 public :: get_dumpname
 public :: add_to_header,add_to_rheader,add_to_iheader
 public :: num_in_header,reset_header,extract
 public :: read_array_from_file
 integer, parameter, public :: lentag = 16    ! tag length
 integer, parameter, public :: lenid  = 100
 integer, parameter, public :: maxphead = 256
 !
 ! Format version. Change this ONLY if the dump format is incompatible
 ! with previous formats
 !
 integer, parameter, public :: iversion = 1

 ! magic numbers
 integer(kind=4), parameter, public :: int1=060769,int2=060878
 integer(kind=4), parameter, public :: int1o=690706,int2o=780806

 ! data types
 integer, parameter, public :: ndatatypes = 8
 integer, parameter, public :: i_int   = 1, &
                               i_int1  = 2, &
                               i_int2  = 3, &
                               i_int4  = 4, &
                               i_int8  = 5, &
                               i_real  = 6, &
                               i_real4 = 7, &
                               i_real8 = 8

 ! error codes
 integer, parameter, public :: ierr_fileopen = 1,&
                               ierr_endian   = 2, &
                               ierr_version  = 3, &
                               ierr_realsize = 4, &
                               ierr_intsize  = 5, &
                               ierr_notags   = 6, &
                               ierr_unknown  = 7

 type dump_h
  integer :: nums(ndatatypes)
  character(len=lentag), allocatable :: inttags(:)
  character(len=lentag), allocatable :: int1tags(:)
  character(len=lentag), allocatable :: int2tags(:)
  character(len=lentag), allocatable :: int4tags(:)
  character(len=lentag), allocatable :: int8tags(:)
  character(len=lentag), allocatable :: realtags(:)
  character(len=lentag), allocatable :: real4tags(:)
  character(len=lentag), allocatable :: real8tags(:)
  integer,         allocatable :: intvals(:)
  integer(kind=1), allocatable :: int1vals(:)
  integer(kind=2), allocatable :: int2vals(:)
  integer(kind=4), allocatable :: int4vals(:)
  integer(kind=8), allocatable :: int8vals(:)
  real,            allocatable :: realvals(:)
  real(kind=4),    allocatable :: real4vals(:)
  real(kind=8),    allocatable :: real8vals(:)
 end type dump_h

 public :: dump_h

 public :: write_header, read_header
 public :: allocate_header, free_header

 ! generic interface to extract quantities from header
 interface extract
  module procedure extract_int4, extract_int8, &
                   extract_real4, extract_real8, &
                   extract_int4arr, extract_int8arr, &
                   extract_real4arr, extract_real8arr
  module procedure extracthdr_int4, extracthdr_int8, &
   extracthdr_real4, extracthdr_real8, &
   extracthdr_int4arr, extracthdr_int8arr, &
   extracthdr_real4arr, extracthdr_real8arr
 end interface extract

 ! generic interface for writing values to header
 interface add_to_header
  module procedure add_to_header_int4, add_to_header_int8, &
    add_to_header_int4arr,  add_to_header_int8arr, &
    add_to_header_real4,    add_to_header_real8, &
    add_to_header_real4arr, add_to_header_real8arr
 end interface add_to_header

 ! add to the default real section of the header
 interface add_to_rheader
  module procedure add_to_rheader, add_to_rheader_arr
 end interface add_to_rheader

 ! add to the default int section of the header
 interface add_to_iheader
  module procedure add_to_iheader, add_to_iheader_arr
 end interface add_to_iheader

 ! generic interface for reset of header
 interface reset_header
  module procedure reset_header_real
 end interface reset_header

 private

contains
!--------------------------------------------------------------------
!+
!  filenaming convention for split MPI dumps
!+
!--------------------------------------------------------------------
function get_dumpname(filename,id)
 character(len=*), intent(in) :: filename
 character(len=len_trim(filename)+8) :: get_dumpname
 integer,          intent(in) :: id

 write(get_dumpname,"(a,a5,i3.3)") trim(filename),'_part',id+1

 return
end function get_dumpname

!--------------------------------------------------------------------
!+
!  small utility to skip an entire block in a file
!+
!-------------------------------------------------------------------
subroutine skipblock(iunit,nums1,nums2,nums3,nums4,ierr)
 integer, intent(in)  :: iunit
 integer, intent(in)  :: nums1(ndatatypes),nums2(ndatatypes),nums3(ndatatypes),nums4(ndatatypes)
 integer, intent(out) :: ierr
 integer :: i

 ierr = 0
 do i=1,sum(nums1)+sum(nums2)+sum(nums3)+sum(nums4)
    read(iunit,iostat=ierr)
 enddo

 return
end subroutine skipblock

!--------------------------------------------------------------------
!+
!  small utility to skip the whole single variable header
!+
!-------------------------------------------------------------------
subroutine skip_headerblock(iunit,ierr)
 integer, intent(in)  :: iunit
 integer, intent(out) :: ierr
 integer :: number

 read(iunit,iostat=ierr) number
 if (ierr /= 0) return
 if (number > 0) then
    read(iunit,iostat=ierr) ! skip tags
    if (ierr /= 0) return
    read(iunit,iostat=ierr) ! skip variables
    if (ierr /= 0) return
 endif

end subroutine skip_headerblock

!---------------------------------------------------------------------
!+
! Construct tag, padded or truncated to 16-characters based on
!  input string
!+
!---------------------------------------------------------------------
elemental function tag(label)
 character(len=lentag) :: tag
 character(len=*), intent(in) :: label

 tag = adjustl(label)

end function tag

!-------------------------------------------
!+
! Check tag against an expected value
!+
!-------------------------------------------
subroutine check_tag(tag,expectedtag)
 character(len=*), intent(in) :: tag, expectedtag

 if (trim(tag)/=trim(expectedtag)) then
    print "(a)",' ERROR reading file: expecting '//trim(expectedtag)//' but got '//trim(tag)
 endif

end subroutine check_tag

!------------------------------------------------
!+
!  Read single value from the header structure
!  For default int and default real we first
!  check the default header, but then check
!  the other headers of the same type
!+
!------------------------------------------------
subroutine extracthdr_int4(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 integer(kind=4),  intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer(kind=4),  intent(in), optional :: default
 integer(kind=4) :: def
 integer :: ival,idef

 if (present(default)) then
    def = default
    idef = default
 else
    def = 0
    idef = 0
 endif
 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.,default=idef)
    if (ierr==0) then
       val = ival
    else
       call extract_int4(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr,default=def)
    endif
 else
    call extract_int4(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr,default=def)
 endif

end subroutine extracthdr_int4

subroutine extracthdr_int8(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 integer(kind=8),  intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer(kind=8),  intent(in), optional :: default
 integer(kind=8) :: def
 integer :: ival,idef

 if (present(default)) then
    def = default
    idef = int(default)
 else
    def = 0_8
    idef = 0
 endif
 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.,default=idef)
    if (ierr==0) then
       val = ival
    else
       call extract_int8(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr,default=def)
    endif
 else
    call extract_int8(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr,default=def)
 endif

end subroutine extracthdr_int8

subroutine extracthdr_real4(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 real(kind=4),     intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real(kind=4),     intent(in), optional :: default
 real(kind=4) :: def
 real :: rval,rdef

 if (present(default)) then
    def = default
    rdef = default
 else
    def = 0.d0
    rdef = 0.
 endif
 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,default=rdef,q=.true.)
    if (ierr==0) then
       val = real(rval,kind=4)
    else
       call extract_real4(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr,default=def)
    endif
 else
    call extract_real4(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr,default=def)
 endif
end subroutine extracthdr_real4

subroutine extracthdr_real8(tag,val,hdr,ierr,default)
 character(len=*), intent(in)  :: tag
 real(kind=8),     intent(out) :: val
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real(kind=8),     intent(in), optional :: default
 real(kind=8) :: def
 real :: rval,rdef
 
 if (present(default)) then
    def = default
    rdef = default
 else
    def = 0.d0
    rdef = 0.
 endif
 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,q=.true.,default=rdef)
    if (ierr==0) then
       val = rval
    else
       call extract_real8(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr,default=def)
    endif
 else
    call extract_real8(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr,default=def)
 endif

end subroutine extracthdr_real8

!------------------------------------------------
!+
!  Read array from the header structure
!+
!------------------------------------------------
subroutine extracthdr_int4arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 integer(kind=4),  intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer :: ival(size(val))

 ierr = 0
 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.)
    if (ierr==0) then
       val = ival
    else
       call extract_int4arr(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr)
    endif
 else
    call extract_int4arr(tag,val,hdr%int4vals,hdr%int4tags,hdr%nums(i_int4),ierr)
 endif

end subroutine extracthdr_int4arr

subroutine extracthdr_int8arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 integer(kind=8),  intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 integer :: ival(size(val))

 if (kind(ival)==kind(val)) then
    call extract(tag,ival,hdr%intvals,hdr%inttags,hdr%nums(i_int),ierr,q=.true.)
    if (ierr==0) then
       val = ival
    else
       call extract_int8arr(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr)
    endif
 else
    call extract_int8arr(tag,val,hdr%int8vals,hdr%int8tags,hdr%nums(i_int8),ierr)
 endif

end subroutine extracthdr_int8arr

subroutine extracthdr_real4arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 real(kind=4),     intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real :: rval(size(val))

 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,q=.true.)
    if (ierr==0) then
       val = real(rval,kind=4)
    else
       call extract_real4arr(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr)
    endif
 else
    call extract_real4arr(tag,val,hdr%real4vals,hdr%real4tags,hdr%nums(i_real4),ierr)
 endif
end subroutine extracthdr_real4arr

subroutine extracthdr_real8arr(tag,val,hdr,ierr)
 character(len=*), intent(in)  :: tag
 real(kind=8),     intent(out) :: val(:)
 type(dump_h),     intent(in)  :: hdr
 integer,          intent(out) :: ierr
 real :: rval(size(val))

 if (kind(rval)==kind(val)) then
    call extract(tag,rval,hdr%realvals,hdr%realtags,hdr%nums(i_real),ierr,q=.true.)
    if (ierr==0) then
       val = rval
    else
       call extract_real8arr(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr)
    endif
 else
    call extract_real8arr(tag,val,hdr%real8vals,hdr%real8tags,hdr%nums(i_real8),ierr)
 endif

end subroutine extracthdr_real8arr

!-------------------------------------------
! Extraction of int*4 variables from header
!-------------------------------------------
subroutine extract_int4(tag,ival,intarr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=4),       intent(out) :: ival
 integer,               intent(in)  :: ntags
 integer(kind=4),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 integer(kind=4),       intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ! default if not found
 if (present(default)) then
    ival = default
 else
    ival = 0
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i) then
          ival = intarr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int4

!-------------------------------------------
! Extraction of int*8 variables from header
!-------------------------------------------
subroutine extract_int8(tag,ival,intarr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=8),       intent(out) :: ival
 integer,               intent(in)  :: ntags
 integer(kind=8),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 integer(kind=8),       intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ! default if not found
 if (present(default)) then
    ival = default
 else
    ival = 0_8
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i) then
          ival = intarr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int8

!---------------------------------------------------
! Extraction of single real*8 variables from header
!---------------------------------------------------
subroutine extract_real8(tag,rval,r8arr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 real(kind=8),          intent(out) :: rval
 real(kind=8),          intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 real(kind=8),          intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 if (present(default)) then
    rval = default
 else
    rval = 0.d0 ! default if not found
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r8arr) >= i) then
          rval = r8arr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real8

!---------------------------------------------------
! Extraction of single real*4 variables from header
!---------------------------------------------------
subroutine extract_real4(tag,rval,r4arr,tags,ntags,ierr,default,q)
 character(len=*),      intent(in)  :: tag
 real(kind=4),          intent(out) :: rval
 real(kind=4),          intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 real(kind=4),          intent(in), optional :: default
 logical,               intent(in), optional :: q
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 if (present(default)) then
    rval = default
 else
    rval = 0. ! default if not found
 endif
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r4arr) >= i) then
          rval = r4arr(i)
          matched = .true.
       endif
       exit over_tags  ! only match first occurrence
    endif
 enddo over_tags
 if (matched) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real4

!------------------------------------------
! Extraction of int*4 arrays frmo header
!------------------------------------------
subroutine extract_int4arr(tag,ival,intarr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=4),       intent(out) :: ival(:)
 integer,               intent(in)  :: ntags
 integer(kind=4),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 ival(:) = 0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i .and. size(ival) > nmatched) then
          nmatched = nmatched + 1
          ival(nmatched) = intarr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(ival)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int4arr

!------------------------------------------
! Extraction of int*8 arrays from header
!------------------------------------------
subroutine extract_int8arr(tag,ival,intarr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 integer(kind=8),       intent(out) :: ival(:)
 integer,               intent(in)  :: ntags
 integer(kind=8),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 ival(:) = 0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(intarr) >= i .and. size(ival) > nmatched) then
          nmatched = nmatched + 1
          ival(nmatched) = intarr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(ival)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_int8arr

!------------------------------------------
! Extraction of real*8 arrays from header
!------------------------------------------
subroutine extract_real8arr(tag,rval,r8arr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 real(kind=8),          intent(out) :: rval(:)
 real(kind=8),          intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 rval = 0.d0 ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r8arr) >= i .and. size(rval) > nmatched) then
          nmatched = nmatched + 1
          rval(nmatched) = r8arr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(rval)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real8arr

!------------------------------------------
! extraction of real*4 arrays from header
!------------------------------------------
subroutine extract_real4arr(tag,rval,r4arr,tags,ntags,ierr,q)
 character(len=*),      intent(in)  :: tag
 real(kind=4),          intent(out) :: rval(:)
 real(kind=4),          intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical,               intent(in), optional :: q
 integer :: i,nmatched

 ierr = 1
 nmatched = 0
 rval = 0. ! default if not found
 over_tags: do i=1,min(ntags,size(tags))
    if (trim(tags(i))==trim(adjustl(tag))) then
       if (size(r4arr) >= i .and. size(rval) > nmatched) then
          nmatched = nmatched + 1
          rval(nmatched) = r4arr(i)
       endif
    endif
 enddo over_tags
 if (nmatched==size(rval)) ierr = 0
 if (ierr /= 0 .and. .not.present(q)) then
    print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'
 endif

end subroutine extract_real4arr

!------------------------------------------
! reset header to all blank entries
!------------------------------------------
subroutine reset_header_real(rheader,rtags)
 real,                  intent(out)   :: rheader(:)
 character(len=lentag), intent(inout) :: rtags(:)
 integer :: i

 do i=1,size(rheader)
    rheader(i) = 0.
    rtags(i)   = ''
 enddo

end subroutine reset_header_real

!------------------------------------------
! add item to int header
!------------------------------------------
subroutine add_to_header_int4(ival,tag,hdr,ierr)
 integer(kind=4),  intent(in)    :: ival
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_int4) + 1
 if (i > size(hdr%int4vals)) then
    ierr = 1
 else
    hdr%nums(i_int4) = i
    hdr%int4vals(i)  = ival
    hdr%int4tags(i)  = tag
 endif

end subroutine add_to_header_int4

!------------------------------------------
! add item to int header
!------------------------------------------
subroutine add_to_header_int8(ival,tag,hdr,ierr)
 integer(kind=8),  intent(in)    :: ival
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_int8) + 1
 if (i > size(hdr%int8vals)) then
    ierr = 1
 else
    hdr%nums(i_int8) = i
    hdr%int8vals(i)  = ival
    hdr%int8tags(i)  = tag
 endif

end subroutine add_to_header_int8

!------------------------------------------
! add array to integer header
!------------------------------------------
subroutine add_to_header_int4arr(ival,tag,hdr,ierr)
 integer(kind=4),  intent(in)    :: ival(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_int4) + 1
 do j=1,size(ival)
    if (i < size(hdr%int4vals)) then
       hdr%nums(i_int4) = i
       hdr%int4vals(i) = ival(j)
       hdr%int4tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_int4arr

!------------------------------------------
! add array to integer*8 header
!------------------------------------------
subroutine add_to_header_int8arr(ival,tag,hdr,ierr)
 integer(kind=8),  intent(in)    :: ival(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_int8) + 1
 do j=1,size(ival)
    if (i < size(hdr%int8vals)) then
       hdr%nums(i_int8) = i
       hdr%int8vals(i) = ival(j)
       hdr%int8tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_int8arr

!------------------------------------------
! add item to real*4 header
!------------------------------------------
subroutine add_to_header_real4(rval,tag,hdr,ierr)
 real(kind=4),     intent(in)    :: rval
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_real4) + 1
 if (i > size(hdr%real4vals)) then
    ierr = 1
 else
    hdr%nums(i_real4) = i
    hdr%real4vals(i)  = rval
    hdr%real4tags(i)  = tag
 endif

end subroutine add_to_header_real4

!------------------------------------------
! add item to real*8 header
!------------------------------------------
subroutine add_to_header_real8(rval,tag,hdr,ierr)
 real(kind=8),     intent(in)    :: rval
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_real8) + 1
 if (i > size(hdr%real8vals)) then
    ierr = 1
 else
    hdr%nums(i_real8) = i
    hdr%real8vals(i)  = rval
    hdr%real8tags(i)  = tag
 endif

end subroutine add_to_header_real8

!------------------------------------------
! add array to real*4 header
!------------------------------------------
subroutine add_to_header_real4arr(rval,tag,hdr,ierr)
 real(kind=4),     intent(in)    :: rval(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_real4) + 1
 do j=1,size(rval)
    if (i < size(hdr%real4vals)) then
       hdr%nums(i_real4) = i
       hdr%real4vals(i) = rval(j)
       hdr%real4tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_real4arr

!------------------------------------------
! add array to real*8 header
!------------------------------------------
subroutine add_to_header_real8arr(rval,tag,hdr,ierr)
 real(kind=8),     intent(in)    :: rval(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_real8) + 1
 do j=1,size(rval)
    if (i < size(hdr%real8vals)) then
       hdr%nums(i_real8) = i
       hdr%real8vals(i) = rval(j)
       hdr%real8tags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_header_real8arr

!------------------------------------------
! add item to default real header
!------------------------------------------
subroutine add_to_rheader(rval,tag,hdr,ierr)
 real,             intent(in)    :: rval
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_real) + 1
 if (i > size(hdr%realvals)) then
    ierr = 1
 else
    hdr%nums(i_real) = i
    hdr%realvals(i)  = rval
    hdr%realtags(i)  = tag
 endif

end subroutine add_to_rheader

!------------------------------------------
! add array to default real header
!------------------------------------------
subroutine add_to_rheader_arr(rval,tag,hdr,ierr)
 real,             intent(in)    :: rval(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_real) + 1
 do j=1,size(rval)
    if (i < size(hdr%realvals)) then
       hdr%nums(i_real) = i
       hdr%realvals(i) = rval(j)
       hdr%realtags(i) = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_rheader_arr

!------------------------------------------
! add item to default int header
!------------------------------------------
subroutine add_to_iheader(ival,tag,hdr,ierr)
 integer,          intent(in)    :: ival
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i

 i = hdr%nums(i_int) + 1
 if (i > size(hdr%intvals)) then
    ierr = 1
 else
    hdr%nums(i_int) = i
    hdr%intvals(i)  = ival
    hdr%inttags(i)  = tag
 endif

end subroutine add_to_iheader

!------------------------------------------
! add array to default int header
!------------------------------------------
subroutine add_to_iheader_arr(ival,tag,hdr,ierr)
 integer,          intent(in)    :: ival(:)
 character(len=*), intent(in)    :: tag
 type(dump_h),     intent(inout) :: hdr
 integer,          intent(inout) :: ierr
 integer :: i,j

 i = hdr%nums(i_int) + 1
 do j=1,size(ival)
    if (i < size(hdr%intvals)) then
       hdr%nums(i_int) = i
       hdr%intvals(i)  = ival(j)
       hdr%inttags(i)  = tag
    else
       ierr = 1
    endif
    i = i + 1
 enddo

end subroutine add_to_iheader_arr

!------------------------------------------
! add item to real header
!------------------------------------------
integer function num_in_header(tags)
 character(len=lentag), intent(in) :: tags(:)
 integer :: i

 ! cycle through header backwards until non-blank value is found
 i = size(tags)
 do while(len_trim(tags(i))==0 .and. i > 0)
    i = i - 1
 enddo
 num_in_header = i

end function num_in_header

!----------------------------------------
! open a dump file and write the file id
! and other generic header information
!----------------------------------------
subroutine open_dumpfile_w(iunit,filename,fileid,ierr,singleprec)
 integer,              intent(in)  :: iunit
 character(len=*),     intent(in)  :: filename
 character(len=lenid), intent(in)  :: fileid
 integer,              intent(out) :: ierr
 logical,              intent(in), optional :: singleprec
 integer :: i1
 logical :: r4
 real    :: r1

 open(unit=iunit,file=filename,status='replace',form='unformatted',iostat=ierr)
 if (ierr /= 0) return

 i1 = int1
 r1 = real(int2)

 r4 = .false.
 if (present(singleprec)) r4 = singleprec

 if (r4) then
    write(iunit,iostat=ierr) int1,real(r1,kind=4),int2,iversion,int1o
 else
    write(iunit,iostat=ierr) int1,r1,int2,iversion,int1o
 endif
 if (ierr /= 0) return

 write(iunit,iostat=ierr) fileid

end subroutine open_dumpfile_w

!-----------------------------------------
! open a dump file and read the file id
! and generic header information
!-----------------------------------------
subroutine open_dumpfile_r(iunit,filename,fileid,ierr,singleprec,requiretags)
 integer,              intent(in)  :: iunit
 character(len=*),     intent(in)  :: filename
 character(len=lenid), intent(out) :: fileid
 integer,              intent(out) :: ierr
 logical,              intent(in), optional :: singleprec,requiretags
 integer(kind=4) :: int1i,int2i,int3i
 integer         :: iversion_file,ierr1
 logical         :: r4,must_have_tags
 real(kind=4)    :: r1s
 real            :: r1i

 r4 = .false.
 must_have_tags = .false.
 if (present(singleprec)) r4 = singleprec
 if (present(requiretags)) must_have_tags = requiretags
!
!--open dump file
!
 open(unit=iunit,file=filename,status='old',form='unformatted',iostat=ierr)
 if (ierr /= 0) then
    ierr = ierr_fileopen
    return
 endif
!
!--read output file
!
 if (r4) then
    read (iunit, iostat=ierr1) int1i,r1s,int2i,iversion_file,int3i
 else
    read (iunit, iostat=ierr1) int1i,r1i,int2i,iversion_file,int3i
 endif
 if (int1i /= int1 .and. int1i /= int1o) then
    ierr = ierr_endian
    return
 endif

! handle version numbers
! version 0 had iversion=690706

 if (iversion_file==int1o) iversion_file = 0
 if (iversion_file > iversion) then
    !write (*,"(a,i2,a,i2)") 'error 4 in readdump: format version is ',iversion, &
    !   ' but this version of Phantom can only read ',maxversion
    ierr = ierr_version
 endif

 read (iunit, iostat=ierr1) fileid

 if (int2i /= int2 .and. int2i /= int2o) then
    ierr = ierr_realsize
    return
 endif
 if (int3i /= int1o) then
    ierr = ierr_intsize
    return
 endif

 ! generic read error, only return this if other errors have been ruled out first
 if (ierr1 /= 0) then
    ierr = ierr_unknown
    return
 endif

 ! return error if not tagged format, if this is required
 if (must_have_tags) then
    if (fileid(2:2) /= 'T' .and. fileid(2:2) /= 't') then
       ierr = ierr_notags
       return
    endif
 endif

 return
end subroutine open_dumpfile_r

!-------------------------------------------------------
!+
!  error handling routine so that errors can be handled
!  outside of this library
!+
!-------------------------------------------------------
character(len=60) function get_error_text(ierr)
 integer, intent(in) :: ierr

 select case(ierr)
 case(ierr_fileopen)
   get_error_text = 'error opening file'
 case(ierr_endian)
   get_error_text = 'wrong endian?'
 case(ierr_version)
   get_error_text = 'file format version newer than current code can read'
 case(ierr_realsize)
   get_error_text = 'default real size wrong'
 case(ierr_intsize)
   get_error_text = 'default int size wrong'
 case(ierr_notags)
   get_error_text = 'routine requires tagged format but not detected'
 case default
   get_error_text = 'unknown error'
 end select

end function get_error_text

!-------------------------------------------------------
!+
!  read the file header into the dump_header structure
!+
!-------------------------------------------------------
subroutine read_header(iunit,hdr,tagged,ierr,singleprec)
 integer,      intent(in) :: iunit
 type(dump_h), intent(out) :: hdr
 logical,      intent(in)  :: tagged
 integer,      intent(out) :: ierr
 logical,      intent(in), optional :: singleprec
 logical :: convert_prec
 integer :: i,n
 real(kind=4), allocatable :: dumr4(:)

 convert_prec = .false.
 if (present(singleprec)) convert_prec = singleprec

 do i=1,ndatatypes
    read (iunit, iostat=ierr) n
    if (n < 0) n = 0
    hdr%nums(i) = n
    select case(i)
    case(i_int)
       allocate(hdr%inttags(n),hdr%intvals(n),stat=ierr)
       hdr%inttags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%inttags(1:n)
          read(iunit, iostat=ierr) hdr%intvals(1:n)
       endif
    case(i_int1)
       allocate(hdr%int1tags(n),hdr%int1vals(n),stat=ierr)
       hdr%int1tags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%int1tags(1:n)
          read(iunit, iostat=ierr) hdr%int1vals(1:n)
       endif
    case(i_int2)
       allocate(hdr%int2tags(n),hdr%int2vals(n),stat=ierr)
       hdr%int2tags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%int2tags(1:n)
          read(iunit, iostat=ierr) hdr%int2vals(1:n)
       endif
    case(i_int4)
       allocate(hdr%int4tags(n),hdr%int4vals(n),stat=ierr)
       hdr%int4tags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%int4tags(1:n)
          read(iunit, iostat=ierr) hdr%int4vals(1:n)
       endif
    case(i_int8)
       allocate(hdr%int8tags(n),hdr%int8vals(n),stat=ierr)
       hdr%int8tags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%int8tags(1:n)
          read(iunit, iostat=ierr) hdr%int8vals(1:n)
       endif
    case(i_real)
       allocate(hdr%realtags(n),hdr%realvals(n),stat=ierr)
       hdr%realtags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%realtags(1:n)
          if (convert_prec .and. kind(0.) /= 4) then
             allocate(dumr4(n),stat=ierr)
             read(iunit, iostat=ierr) dumr4(1:n)
             hdr%realvals(1:n) = real(dumr4(1:n))
             deallocate(dumr4)
          else
             read(iunit, iostat=ierr) hdr%realvals(1:n)
          endif
       endif
    case(i_real4)
       allocate(hdr%real4tags(n),hdr%real4vals(n),stat=ierr)
       hdr%real4tags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%real4tags(1:n)
          read(iunit, iostat=ierr) hdr%real4vals(1:n)
       endif
    case(i_real8)
       allocate(hdr%real8tags(n),hdr%real8vals(n),stat=ierr)
       hdr%real8tags(:) = ''
       if (n > 0) then
          if (tagged) read(iunit, iostat=ierr) hdr%real8tags(1:n)
          read(iunit, iostat=ierr) hdr%real8vals(1:n)
       endif
    end select
 enddo

end subroutine read_header

!-------------------------------------------------------
!+
!  allocate the dump header structure for writing
!  IN: (all are optional, default size is 256)
!     nint,nint2 : size of buffer for each data type
!               e.g. (/256,0,0,256,256,256,256,256/)
!  OUT:
!     hdr : header structure allocated for each
!           data type
!+
!-------------------------------------------------------
function allocate_header(nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8,err) result(hdr)
 integer, intent(in),  optional :: nint,nint1,nint2,nint4,nint8,nreal,nreal4,nreal8
 integer, intent(out), optional :: err
 type(dump_h) :: hdr
 integer      :: size(ndatatypes)
 integer      :: ierrs(ndatatypes)
 integer      :: ierr

 ! make sure header is deallocated first
 call free_header(hdr,ierr)

 size(:) = maxphead
 if (present(nint))   size(i_int)  = nint
 if (present(nint1))  size(i_int1) = nint1
 if (present(nint2))  size(i_int2) = nint2
 if (present(nint4))  size(i_int4) = nint4
 if (present(nreal))  size(i_real) = nreal
 if (present(nreal4)) size(i_real4) = nreal4
 if (present(nreal8)) size(i_real8) = nreal8

 if (present(err)) err = 0
 ierrs(:) = 0
 hdr%nums(:) = 0
 if (size(i_int) > 0)  then
    allocate(hdr%inttags(size(i_int)),hdr%intvals(size(i_int)),stat=ierrs(1))
    if (ierrs(1)==0) hdr%inttags(:) = ''
 endif
 if (size(i_int1) > 0) then
    allocate(hdr%int1tags(size(i_int1)),hdr%int1vals(size(i_int1)),stat=ierrs(2))
    if (ierrs(2)==0) hdr%int1tags(:) = ''
 endif
 if (size(i_int2) > 0) then
    allocate(hdr%int2tags(size(i_int2)),hdr%int2vals(size(i_int2)),stat=ierrs(3))
    if (ierrs(3)==0) hdr%int2tags(:) = ''
 endif
 if (size(i_int4) > 0) then
    allocate(hdr%int4tags(size(i_int4)),hdr%int4vals(size(i_int4)),stat=ierrs(4))
    if (ierrs(4)==0) hdr%int4tags(:) = ''
 endif
 if (size(i_int8) > 0) then
    allocate(hdr%int8tags(size(i_int8)),hdr%int8vals(size(i_int8)),stat=ierrs(5))
    if (ierrs(5)==0) hdr%int8tags(:) = ''
 endif
 if (size(i_real) > 0)  then
    allocate(hdr%realtags(size(i_real)),hdr%realvals(size(i_real)),stat=ierrs(6))
    if (ierrs(6)==0) hdr%realtags(:) = ''
 endif
 if (size(i_real4) > 0)  then
    allocate(hdr%real4tags(size(i_real4)),hdr%real4vals(size(i_real4)),stat=ierrs(7))
    if (ierrs(7)==0) hdr%real4tags(:) = ''
 endif
 if (size(i_real8) > 0)  then
    allocate(hdr%real8tags(size(i_real8)),hdr%real8vals(size(i_real8)),stat=ierrs(8))
    if (ierrs(8)==0) hdr%real8tags(:) = ''
 endif

 if (present(err) .and. any(ierrs /= 0)) err = 1

end function allocate_header

!-------------------------------------------------------
!+
!  allocate the dump header structure for writing
!+
!-------------------------------------------------------
subroutine free_header(hdr,ierr)
 type(dump_h), intent(inout) :: hdr
 integer,      intent(out)   :: ierr

 if (allocated(hdr%inttags))   deallocate(hdr%inttags)
 if (allocated(hdr%int1tags))  deallocate(hdr%int1tags)
 if (allocated(hdr%int2tags))  deallocate(hdr%int2tags)
 if (allocated(hdr%int4tags))  deallocate(hdr%int4tags)
 if (allocated(hdr%int8tags))  deallocate(hdr%int8tags)
 if (allocated(hdr%realtags))  deallocate(hdr%realtags)
 if (allocated(hdr%real4tags)) deallocate(hdr%real4tags)
 if (allocated(hdr%real8tags)) deallocate(hdr%real8tags)

 if (allocated(hdr%intvals))   deallocate(hdr%intvals)
 if (allocated(hdr%int1vals))  deallocate(hdr%int1vals)
 if (allocated(hdr%int2vals))  deallocate(hdr%int2vals)
 if (allocated(hdr%int4vals))  deallocate(hdr%int4vals)
 if (allocated(hdr%realvals))  deallocate(hdr%realvals)
 if (allocated(hdr%real4vals)) deallocate(hdr%real4vals)
 if (allocated(hdr%real8vals)) deallocate(hdr%real8vals)

end subroutine free_header

!-------------------------------------------------------
!+
!  write the header to file
!+
!-------------------------------------------------------
subroutine write_header(iunit,hdr,ierr,singleprec)
 integer,      intent(in)  :: iunit
 type(dump_h), intent(in)  :: hdr
 integer,      intent(out) :: ierr
 logical,      intent(in), optional :: singleprec
 integer :: idum,number,i,j,ierrs(17)
 integer(kind=1), parameter :: idum1 = 0_1
 integer(kind=2), parameter :: idum2 = 0_2
 integer(kind=4), parameter :: idum4 = 0_4
 integer(kind=8), parameter :: idum8 = 0_8
 real, parameter :: dum = 0.
 real(kind=4), parameter :: dum4 = 0._4
 real(kind=8), parameter :: dum8 = 0._8
 logical :: sing_prec

 ! optional argument to write real header in single precision
 sing_prec = .false.
 if (present(singleprec)) sing_prec = singleprec

 ierrs = 0
 do i=1,ndatatypes
    number = hdr%nums(i)
    write(iunit,iostat=ierrs(1)) number
    if (number > 0) then
       select case(i)
       case(i_int)
          if (allocated(hdr%inttags) .and. allocated(hdr%intvals)) then
             write(iunit,iostat=ierrs(2)) hdr%inttags(1:number)
             write(iunit,iostat=ierrs(3)) hdr%intvals(1:number)
          else
             write(iunit,iostat=ierrs(2)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(3)) (idum,j=1,number)
          endif
       case(i_int1)
          if (allocated(hdr%int1tags) .and. allocated(hdr%int1vals)) then
             write(iunit,iostat=ierrs(4)) hdr%int1tags(1:number)
             write(iunit,iostat=ierrs(5)) hdr%int1vals(1:number)
          else
             write(iunit,iostat=ierrs(4)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(5)) (idum1,j=1,number)
          endif
       case(i_int2)
          if (allocated(hdr%int2tags) .and. allocated(hdr%int2vals)) then
             write(iunit,iostat=ierrs(6)) hdr%int2tags(1:number)
             write(iunit,iostat=ierrs(7)) hdr%int2vals(1:number)
          else
             write(iunit,iostat=ierrs(6)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(7)) (idum2,j=1,number)
          endif
       case(i_int4)
          if (allocated(hdr%int4tags) .and. allocated(hdr%int4vals)) then
             write(iunit,iostat=ierrs(8)) hdr%int4tags(1:number)
             write(iunit,iostat=ierrs(9)) hdr%int4vals(1:number)
          else
             write(iunit,iostat=ierrs(8)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(9)) (idum4,j=1,number)
          endif
       case(i_int8)
          if (allocated(hdr%int8tags) .and. allocated(hdr%int8vals)) then
             write(iunit,iostat=ierrs(10)) hdr%int8tags(1:number)
             write(iunit,iostat=ierrs(11)) hdr%int8vals(1:number)
          else
             write(iunit,iostat=ierrs(10)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(11)) (idum8,j=1,number)
          endif
       case(i_real)
          if (allocated(hdr%realtags) .and. allocated(hdr%realvals)) then
             write(iunit,iostat=ierrs(12)) hdr%realtags(1:number)
             if (sing_prec) then
                write(iunit,iostat=ierrs(13)) real(hdr%realvals(1:number),kind=4)
             else
                write(iunit,iostat=ierrs(13)) hdr%realvals(1:number)
             endif
          else
             write(iunit,iostat=ierrs(12)) (tag('unknown'),j=1,number)
             if (sing_prec) then
                write(iunit,iostat=ierrs(13)) (real(dum,kind=4),j=1,number)
             else
                write(iunit,iostat=ierrs(13)) (dum,j=1,number)
             endif
          endif
       case(i_real4)
          if (allocated(hdr%real4tags) .and. allocated(hdr%real4vals)) then
             write(iunit,iostat=ierrs(14)) hdr%real4tags(1:number)
             write(iunit,iostat=ierrs(15)) hdr%real4vals(1:number)
          else
             write(iunit,iostat=ierrs(14)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(15)) (dum4,j=1,number)
          endif
       case(i_real8)
          if (allocated(hdr%real8tags) .and. allocated(hdr%real8vals)) then
             write(iunit,iostat=ierrs(16)) hdr%real8tags(1:number)
             write(iunit,iostat=ierrs(17)) hdr%real8vals(1:number)
          else
             write(iunit,iostat=ierrs(16)) (tag('unknown'),j=1,number)
             write(iunit,iostat=ierrs(17)) (dum8,j=1,number)
          endif
       end select
    endif
 enddo
 if (any(ierrs /= 0)) ierr = 1

end subroutine write_header

!-----------------------------------------------------
!+
!  The following routine can be used to read a single
!  array matching a particular tag from the main blocks
!  in the file
!+
!-----------------------------------------------------
subroutine read_array_from_file(iunit,filename,tag,array,ierr)
 integer,               intent(in) :: iunit
 character(len=*),      intent(in) :: filename
 character(len=*),      intent(in) :: tag
 real,    intent(out) :: array(:)
 integer, intent(out) :: ierr
 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,nblockarrays,number
 integer :: intarr(maxphead)
 character(len=lentag) :: tagarr(maxphead)
 character(len=lenid)  :: fileid

 array = 0.

 ! open file for read
 call open_dumpfile_r(iunit,filename,fileid,ierr,requiretags=.true.)
 if (ierr /= 0) return

 ! read nblocks from int header
 read(iunit,iostat=ierr) number
 nblocks = 1
 if (number >= 5) then
    if (number > maxphead) number = maxphead
    read(iunit,iostat=ierr) tagarr(1:number)
    read(iunit,iostat=ierr) intarr(1:number)
    call extract('nblocks',nblocks,intarr,tagarr,number,ierr)
    if (ierr /= 0) nblocks = 1
 elseif (number > 0) then
    nblocks = 1
    read(iunit,iostat=ierr)
    read(iunit,iostat=ierr)
 endif

 ! no need to read rest of header
 do i=1,ndatatypes-1
    call skip_headerblock(iunit,ierr)
    if (ierr /= 0) print*,' error skipping header block'
    if (ierr /= 0) return
 enddo

 read (iunit, iostat=ierr) number
 if (ierr /= 0) return
 narraylengths = number/nblocks
! print*,' got nblocks = ',nblocks,' narraylengths = ',narraylengths

 ! skip each block that is too small
 nblockarrays = narraylengths*nblocks
 do iblock = 1,nblocks
    call read_blockheader(iunit,narraylengths,number8,nums,ierr)
    do j=1,narraylengths
       !print*,'block ',j
       do i=1,ndatatypes
          !print*,' data type ',i,' arrays = ',nums(i,j)
          do k=1,nums(i,j)
             if (i==i_real) then
                read(iunit, iostat=ierr) tagarr(1)
                if (trim(tagarr(1))==trim(tag)) then
                   read(iunit, iostat=ierr) array(1:min(int(number8(j)),size(array)))
                   print*,'->',tagarr(1)
                else
                   print*,'  ',tagarr(1)
                   read(iunit, iostat=ierr)
                endif
             else
                read(iunit, iostat=ierr) tagarr(1) ! tag
                !print*,tagarr(1)
                read(iunit, iostat=ierr) ! array
             endif
          enddo
       enddo
    enddo
 enddo

 close(iunit)

end subroutine read_array_from_file

!--------------------------------------------------------------
!+
!  internal routine to read the header arrays for each block
!+
!--------------------------------------------------------------
subroutine read_blockheader(iunit,narraylengths,number8,nums,ierr)
 integer,         intent(in)  :: iunit,narraylengths
 integer(kind=8), intent(out) :: number8(narraylengths)
 integer,         intent(out) :: nums(ndatatypes,narraylengths)
 integer,         intent(out) :: ierr
 integer :: i,j

 do i=1,narraylengths
    read(iunit,iostat=ierr) number8(i), (nums(j,i), j=1,ndatatypes)
 enddo

end subroutine read_blockheader

end module dump_utils
