module io_phantom

  use parametres
  use dump_utils

  implicit none

  contains

subroutine read_phantom_file(iunit,filename,x,y,z,rhogas,rhodust,ncells,ierr)
 integer,               intent(in) :: iunit
 character(len=*),      intent(in) :: filename
 real(db), intent(out), dimension(:), allocatable :: x,y,z,rhogas,rhodust
 integer, intent(out) :: ncells,ierr
 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,nblockarrays,number
 integer :: intarr(maxphead)
 character(len=lentag) :: tagarr(maxphead)
 character(len=lenid)  :: fileid
 integer :: np,ntypes
 integer, parameter :: maxtypes = 10
 integer :: npartoftype(maxtypes)
 real(db) :: massoftype(maxtypes),realarr(maxphead),hfact
 integer, allocatable, dimension(:) :: itype
 real(4), allocatable, dimension(:) :: tmp,dustfrac
 real(db), allocatable, dimension(:,:) :: xyzh

 logical :: got_h, got_dustfrac

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

 call extract('nparttot',np,intarr,tagarr,number,ierr)
 call extract('ntypes',ntypes,intarr,tagarr,number,ierr)
 call extract('npartoftype',npartoftype(1:ntypes),intarr,tagarr,number,ierr)

 print*,' npart = ',np,' ntypes = ',ntypes
 print*,' npartoftype = ',npartoftype(1:ntypes)

 allocate(xyzh(4,np),itype(np),dustfrac(np),tmp(np))
 itype = 1

 ! no need to read rest of header
 do i=i_int+1,i_real-1
    call skip_headerblock(iunit,ierr)
    if (ierr /= 0) print*,' error skipping header block'
    if (ierr /= 0) return
 enddo

 ! read real header
 read(iunit,iostat=ierr) number
 number = min(number,maxphead)
 if (number > 0) then
    read(iunit,iostat=ierr) tagarr(1:number)
    read(iunit,iostat=ierr) realarr(1:number)
 endif
 !print*,tagarr(1:number)

 ! skip rest of header
 do i=i_real+1,ndatatypes
    call skip_headerblock(iunit,ierr)
 enddo

 ! extract info from real header
 call extract('massoftype',massoftype(1:ntypes),realarr,tagarr,number,ierr)
 call extract('hfact',hfact,realarr,tagarr,number,ierr)
 print*,' hfact = ',hfact
 print*,' massoftype = ',massoftype(1:ntypes)

 read (iunit, iostat=ierr) number
 if (ierr /= 0) return
 narraylengths = number/nblocks

 got_h = .false.
 got_dustfrac = .false.
 ! skip each block that is too small
 nblockarrays = narraylengths*nblocks
 do iblock = 1,nblocks
    call read_blockheader(iunit,narraylengths,number8,nums,ierr)
    do j=1,narraylengths
       !print*,'block ',j
       do i=1,ndatatypes
          !print*,' data type ',i,' arrays = ',nums(i,j)
          do k=1,nums(i,j)
             if (j==1 .and. number8(j)==np) then
                read(iunit, iostat=ierr) tagarr(1)
                print*,' ',trim(tagarr(1))
                if (i==i_real .or. i==i_real8) then
                   select case(trim(tagarr(1)))
                   case('x')
                      print*,'->',tagarr(1)
                      read(iunit,iostat=ierr) xyzh(1,1:np)
                   case('y')
                      print*,'->',tagarr(1)
                      read(iunit,iostat=ierr) xyzh(2,1:np)
                   case('z')
                      print*,'->',tagarr(1)
                      read(iunit,iostat=ierr) xyzh(3,1:np)
                   case('h')
                      print*,'->',tagarr(1)
                      read(iunit,iostat=ierr) xyzh(4,1:np)
                      got_h = .true.
                   case default
                      read(iunit,iostat=ierr)
                   end select
                elseif (i==i_real4) then
                   select case(trim(tagarr(1)))
                   case('h')
                      print*,'->',tagarr(1)
                      read(iunit,iostat=ierr) tmp(1:np)
                      xyzh(4,1:np) = real(tmp(1:np),kind=db)
                      got_h = .true.
                   case('dustfrac')
                      print*,'->',tagarr(1)
                      read(iunit,iostat=ierr) dustfrac(1:np)
                      got_dustfrac = .true.
                   case default
                      read(iunit,iostat=ierr)
                   end select
                elseif (i==i_int .or. i==i_int4) then
                   select case(trim(tagarr(1)))
                   case('itype')
                      read(iunit,iostat=ierr) itype(1:np)
                   case default
                      read(iunit,iostat=ierr)
                   end select
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

 if (.not.got_dustfrac) then
    dustfrac = 0.
 endif

 if (got_h) then
    call phantom2mcfost(np,xyzh,itype,dustfrac,ntypes,massoftype(1:ntypes),hfact,x,y,z,rhogas,rhodust,ncells)
 else
    ncells = 0
    print*,' ERROR reading h from file'
 endif

 deallocate(xyzh,itype,dustfrac,tmp)

end subroutine read_phantom_file

!*************************************************************************

pure subroutine phantom_2_mcfost(np,xyzh,iphase,dustfrac,ntypes,massoftype,hfact,x,y,z,rhogas,rhodust,ncells)

  integer, intent(in) :: np, ntypes
  real(db), dimension(4,np), intent(in) :: xyzh
  integer, dimension(np), intent(in) :: iphase
  real(sl), dimension(np), intent(in) :: dustfrac
  real(db), dimension(ntypes), intent(in) :: massoftype
  real(db), intent(in) :: hfact

  real(db), dimension(:), allocatable, intent(out) :: x,y,z,rhogas,rhodust
  integer, intent(out) :: ncells

  integer :: i,j,itypei
  real(db) :: xi, yi, zi, hi, rhoi

 ! convert to dust and gas density
 j = 0
 do i=1,np
    if (xyzh(4,i) > 0. .and. abs(iphase(i))==1)  j = j + 1
 enddo
 ncells = j
 allocate(x(ncells),y(ncells),z(ncells),rhogas(ncells),rhodust(ncells))

 j = 0
 do i=1,np
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    itypei = abs(iphase(i))
    if (hi > 0. .and. itypei==1) then
       j = j + 1
       x(j) = xi
       y(j) = yi
       z(j) = zi
       rhoi = massoftype(itypei)*(hfact/hi)**3
       rhogas(j) = (1 - dustfrac(i))*rhoi
       rhodust(j) = dustfrac(i)*rhoi
    endif
 enddo
 ncells = j

 return

end subroutine phantom_2_mcfost

subroutine read_phantom_input_file(filename,iunit,grainsize,graindens,ierr)
  use infile_utils
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: filename
  real, intent(out) :: grainsize,graindens
  integer,  intent(out) :: ierr
  type(inopts), allocatable :: dbin(:)

  call open_db_from_file(dbin,filename,iunit,ierr)
  call read_inopt(graindens,'graindens',dbin,ierr)
  call read_inopt(grainsize,'grainsize',dbin,ierr)
  call close_db(dbin)

end subroutine read_phantom_input_file

end module io_phantom
