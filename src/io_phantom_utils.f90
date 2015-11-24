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
!  $Id: 3ae57c9c8e683de7f0b4902adf0cd92d8a7fb3bd $
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
 public :: add_to_header,num_in_header,reset_header,extract
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

 ! generic interface to extract quantities from header
 interface extract
  module procedure extract_int4, extract_int8, &
                   extract_real4, extract_real8, &
                   extract_int4arr, extract_int8arr, &
                   extract_real4arr, extract_real8arr
 end interface extract

 ! generic interface for writing values to header
 interface add_to_header
  module procedure add_to_header_real, add_to_header_realarr
 end interface add_to_header

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

!-------------------------------------------
! Extraction of int*4 variables from header
!-------------------------------------------
subroutine extract_int4(tag,ival,intarr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 integer(kind=4),       intent(out) :: ival
 integer,               intent(in)  :: ntags
 integer(kind=4),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ival = 0 ! default if not found
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_int4

!-------------------------------------------
! Extraction of int*8 variables from header
!-------------------------------------------
subroutine extract_int8(tag,ival,intarr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 integer(kind=8),       intent(out) :: ival
 integer,               intent(in)  :: ntags
 integer(kind=8),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 ival = 0 ! default if not found
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_int8

!---------------------------------------------------
! Extraction of single real*8 variables from header
!---------------------------------------------------
subroutine extract_real8(tag,rval,r8arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real(kind=8),          intent(out) :: rval
 real(kind=8),          intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 rval = 0.d0 ! default if not found
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real8

!---------------------------------------------------
! Extraction of single real*4 variables from header
!---------------------------------------------------
subroutine extract_real4(tag,rval,r4arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real(kind=4),          intent(out) :: rval
 real(kind=4),          intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
 logical :: matched
 integer :: i

 ierr = 1
 matched = .false.
 rval = 0. ! default if not found
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real4

!------------------------------------------
! Extraction of int*4 arrays frmo header
!------------------------------------------
subroutine extract_int4arr(tag,ival,intarr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 integer(kind=4),       intent(out) :: ival(:)
 integer,               intent(in)  :: ntags
 integer(kind=4),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_int4arr

!------------------------------------------
! Extraction of int*8 arrays from header
!------------------------------------------
subroutine extract_int8arr(tag,ival,intarr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 integer(kind=8),       intent(out) :: ival(:)
 integer,               intent(in)  :: ntags
 integer(kind=8),       intent(in)  :: intarr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(out) :: ierr
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_int8arr

!------------------------------------------
! Extraction of real*8 arrays from header
!------------------------------------------
subroutine extract_real8arr(tag,rval,r8arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real(kind=8),          intent(out) :: rval(:)
 real(kind=8),          intent(in)  :: r8arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

end subroutine extract_real8arr

!------------------------------------------
! extraction of real*4 arrays from header
!------------------------------------------
subroutine extract_real4arr(tag,rval,r4arr,tags,ntags,ierr)
 character(len=*),      intent(in)  :: tag
 real(kind=4),          intent(out) :: rval(:)
 real(kind=4),          intent(in)  :: r4arr(:)
 character(len=lentag), intent(in)  :: tags(:)
 integer,               intent(in)  :: ntags
 integer,               intent(out) :: ierr
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
 if (ierr /= 0) print "(a)",' ERROR: could not find '//trim(adjustl(tag))//' in header'

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
! add item to real header
!------------------------------------------
subroutine add_to_header_real(rval,tag,rheader,rtags,ierr)
 real,                  intent(in)    :: rval
 character(len=*),      intent(in)    :: tag
 real,                  intent(inout) :: rheader(:)
 character(len=lentag), intent(inout) :: rtags(:)
 integer,               intent(inout) :: ierr
 integer :: i

 !ierr = 0
 ! cycle through header until non-blank value is found
 i = 1
 do while(len_trim(rtags(i)) > 0 .and. i < size(rtags))
    i = i + 1
 enddo
 if (i==size(rheader) .and. (len_trim(rtags(i)) > 0)) then
    ierr = 1
 else
    rheader(i) = rval
    rtags(i)   = tag
 endif

end subroutine add_to_header_real

!------------------------------------------
! add array to real header
!------------------------------------------
subroutine add_to_header_realarr(rval,tag,rheader,rtags,ierr)
 real,                  intent(in)    :: rval(:)
 character(len=*),      intent(in)    :: tag
 real,                  intent(inout) :: rheader(:)
 character(len=lentag), intent(inout) :: rtags(:)
 integer,               intent(inout) :: ierr
 integer :: i,j

 ! cycle through header until non-blank value is found
 i = 1
 do while(len_trim(rtags(i)) > 0 .and. i < size(rtags))
    i = i + 1
 enddo
 if (i==size(rheader) .and. (len_trim(rtags(i)) > 0)) then
    ierr = 1
 else
    do j=1,size(rval)
       if (i <= size(rheader)) then
          rheader(i) = rval(j)
          rtags(i)   = tag
       else
          ierr = 2
       endif
       i = i + 1
    enddo
 endif

end subroutine add_to_header_realarr

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
                   read(iunit, iostat=ierr) array(1:min(number8(j),size(array))
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
