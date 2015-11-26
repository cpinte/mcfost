module io_phantom

  use parametres
  use dump_utils

  implicit none

  contains

subroutine read_phantom_file(iunit,filename,x,y,z,rhogas,rhodust,ndusttypes,ncells,ierr)
 integer,               intent(in) :: iunit
 character(len=*),      intent(in) :: filename
 real(db), intent(out), dimension(:),   allocatable :: x,y,z,rhogas
 real(db), intent(out), dimension(:,:), allocatable :: rhodust
 integer, intent(out) :: ndusttypes,ncells,ierr
 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,nblockarrays,number
 integer :: intarr(maxphead)
 character(len=lentag) :: tagarr(maxphead)
 character(len=lenid)  :: fileid
 integer :: np,ntypes,nptmass,ipos,ngrains
 integer, parameter :: maxtypes = 10
 integer :: npartoftype(maxtypes)
 real(db) :: massoftype(maxtypes),realarr(maxphead),hfact,umass,utime,udist
 integer, allocatable, dimension(:) :: itype
 real(4), allocatable, dimension(:) :: tmp
 real(db), allocatable, dimension(:) :: grainsize
 real(db) :: graindens
 real(db), allocatable, dimension(:,:) :: xyzh,dustfrac,xyzmh_ptmass

 logical :: got_h, got_dustfrac

 ! open file for read
 call open_dumpfile_r(iunit,filename,fileid,ierr,requiretags=.true.)
 print "(a)",trim(fileid)
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

 nptmass = 0
 call extract('nparttot',np,intarr,tagarr,number,ierr)
 call extract('ntypes',ntypes,intarr,tagarr,number,ierr)
 call extract('npartoftype',npartoftype(1:ntypes),intarr,tagarr,number,ierr)
 call extract('ndusttypes',ndusttypes,intarr,tagarr,number,ierr)
 call extract('nptmass',nptmass,intarr,tagarr,number,ierr)
 if (ndusttypes==0) ndusttypes = 1

 !print*,' npart = ',np,' ntypes = ',ntypes
 !print*,' npartoftype = ',npartoftype(1:ntypes)

 allocate(xyzh(4,np),itype(np),dustfrac(ndusttypes,np),grainsize(ndusttypes),tmp(np))
 allocate(xyzmh_ptmass(5,nptmass))
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

 ! extract info from real header
 call extract('massoftype',massoftype(1:ntypes),realarr,tagarr,number,ierr)
 call extract('hfact',hfact,realarr,tagarr,number,ierr)
 call extract('grainsize',grainsize(1:ndusttypes),realarr,tagarr,number,ierr)
 call extract('graindens',graindens,realarr,tagarr,number,ierr)
 print*,' hfact = ',hfact
 print*,' massoftype = ',massoftype(1:ntypes)

 ! skip real*4 header
 call skip_headerblock(iunit,ierr)

 ! read units from real*8 header
 read(iunit,iostat=ierr) number
 number = min(number,maxphead)
 if (number > 0) then
    read(iunit,iostat=ierr) tagarr(1:number)
    read(iunit,iostat=ierr) realarr(1:number)
 endif
 call extract('umass',umass,realarr,tagarr,number,ierr)
 call extract('utime',utime,realarr,tagarr,number,ierr)
 call extract('udist',udist,realarr,tagarr,number,ierr)

 read (iunit, iostat=ierr) number
 if (ierr /= 0) return
 narraylengths = number/nblocks

 got_h = .false.
 got_dustfrac = .false.
 ! skip each block that is too small
 nblockarrays = narraylengths*nblocks

 ngrains = 0
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
                   case('dustfrac')
                      ngrains = ngrains + 1
                      read(iunit,iostat=ierr) dustfrac(ngrains,1:np)
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
                      ngrains = ngrains + 1
                      read(iunit,iostat=ierr) tmp(1:np)
                      dustfrac(ngrains,1:np) = real(tmp(1:np),kind=db)
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
             elseif (j==1 .and. number8(j)==nptmass) then
                read(iunit,iostat=ierr) tagarr(1)
                if (i==i_real .or. i==i_real8) then
                   ipos = 0
                   select case(trim(tagarr(1)))
                   case('x')
                      ipos = 1
                   case('y')
                      ipos = 2
                   case('z')
                      ipos = 3
                   case('m')
                      ipos = 4
                   case('h')
                      ipos = 5
                   end select
                   if (ipos > 0) then
                      read(iunit,iostat=ierr) xyzmh_ptmass(ipos,1:nptmass)
                   endif
                else
                   read(iunit,iostat=ierr)
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
    call phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,xyzh,itype,grainsize,dustfrac,&
         massoftype(1:ntypes),xyzmh_ptmass,hfact,umass,utime,udist,graindens,x,y,z,rhogas,rhodust,ncells)
 else
    ncells = 0
    print*,' ERROR reading h from file'
 endif

 deallocate(xyzh,itype,dustfrac,tmp,xyzmh_ptmass)

end subroutine read_phantom_file

!*************************************************************************

subroutine phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,xyzh,iphase,grainsize,dustfrac,&
     massoftype,xyzmh_ptmass,hfact,umass,utime,udist,graindens,x,y,z,rhogas,rhodust,ncells)

  use constantes, only : au_to_cm, Msun_to_g
  use prop_star

  integer, intent(in) :: np, nptmass, ntypes, ndusttypes
  real(db), dimension(4,np), intent(in) :: xyzh
  integer, dimension(np), intent(in) :: iphase
  real(db), dimension(ndusttypes,np), intent(in) :: dustfrac
  real(db), dimension(ndusttypes),    intent(in) :: grainsize
  real(db), intent(in) :: graindens
  real(db), dimension(ntypes), intent(in) :: massoftype
  real(db), intent(in) :: hfact,umass,utime,udist
  real(db), dimension(:,:), intent(in) :: xyzmh_ptmass

  real(db), dimension(:),   allocatable, intent(out) :: x,y,z,rhogas
  real(db), dimension(:,:), allocatable, intent(out) :: rhodust
  integer, intent(out) :: ncells

  integer :: i,j,k,itypei
  real(db) :: xi, yi, zi, hi, rhoi, udens, ulength, usolarmass, dustfraci

  udens = umass/udist**3
  ulength = udist/au_to_cm
  usolarmass = umass/Msun_to_g

 ! convert to dust and gas density
 j = 0
 do i=1,np
    if (xyzh(4,i) > 0. .and. abs(iphase(i))==1)  j = j + 1
 enddo
 ncells = j
 allocate(x(ncells),y(ncells),z(ncells),rhogas(ncells),rhodust(ndusttypes,ncells))

 j = 0
 do i=1,np
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    itypei = abs(iphase(i))
    if (hi > 0. .and. itypei==1) then
       j = j + 1
       x(j) = xi * ulength
       y(j) = yi * ulength
       z(j) = zi * ulength
       rhoi = massoftype(itypei)*(hfact/hi)**3  * udens ! g/cm**3
       dustfraci = sum(dustfrac(:,i))
       rhogas(j) = (1 - dustfraci)*rhoi
       do k=1,ndusttypes
          rhodust(k,j) = dustfrac(k,i)*rhoi
       enddo
    endif
 enddo
 ncells = j


 allocate(etoile(nptmass))
 do i=1,nptmass
    etoile(i)%x = xyzmh_ptmass(1,i) * ulength
    etoile(i)%y = xyzmh_ptmass(2,i) * ulength
    etoile(i)%z = xyzmh_ptmass(3,i) * ulength

    etoile(i)%M = xyzmh_ptmass(4,i) * usolarmass

    ! T, fUV, slope_UV, lb_body, spectre, ri, zj, phik
 enddo

 return

end subroutine phantom_2_mcfost

subroutine read_phantom_input_file(filename,iunit,graintype,graindens,ierr)
  use infile_utils
  integer, intent(in) :: iunit
  character(len=*), intent(in) :: filename
  real, intent(out) :: graintype,graindens
  integer,  intent(out) :: ierr
  type(inopts), allocatable :: dbin(:)

  call open_db_from_file(dbin,filename,iunit,ierr)
  call read_inopt(graindens,'graindens',dbin,ierr)
  call read_inopt(graintype,'graintype',dbin,ierr)
  call close_db(dbin)

end subroutine read_phantom_input_file

end module io_phantom
