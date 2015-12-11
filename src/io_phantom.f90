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
 character(len=lentag) :: tag
 character(len=lenid)  :: fileid
 integer :: np,ntypes,nptmass,ipos,ngrains
 integer, parameter :: maxtypes = 10
 integer :: npartoftype(maxtypes)
 real(db) :: massoftype(maxtypes),hfact,umass,utime,udist
 integer(kind=1), allocatable, dimension(:) :: itype
 real(4),  allocatable, dimension(:) :: tmp
 real(db), allocatable, dimension(:) :: grainsize
 real(db) :: graindens
 real(db), allocatable, dimension(:,:) :: xyzh,dustfrac,xyzmh_ptmass
 type(dump_h) :: hdr
 logical :: got_h, got_dustfrac, tagged, matched

 ! open file for read
 call open_dumpfile_r(iunit,filename,fileid,ierr,requiretags=.true.)
 if (ierr /= 0) then
    write(*,"(/,a,/)") ' *** ERROR - '//trim(get_error_text(ierr))//' ***'
    return
 endif
 print "(1x,a)",trim(fileid)

 call read_header(iunit,hdr,tagged,ierr)
 if (.not.tagged) then
    write(*,"(/,a,/)") ' *** ERROR - Phantom dump too old to be read by MCFOST ***'
    ierr = 1000
    return
 endif
 !print*,hdr%inttags(:)
 !print*,hdr%realtags(:)
 
 ! get nblocks
 call extract('nblocks',nblocks,hdr,ierr,default=1)
 call extract('nparttot',np,hdr,ierr)
 call extract('ntypes',ntypes,hdr,ierr)
 call extract('npartoftype',npartoftype(1:ntypes),hdr,ierr)
 call extract('ndusttypes',ndusttypes,hdr,ierr,default=1)
 call extract('nptmass',nptmass,hdr,ierr,default=0)

 print*,' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
 print*,' npartoftype = ',npartoftype(1:ntypes)
 if (npartoftype(2) > 0) then
    write(*,"(/,a,/)") ' *** WARNING: Phantom dump contains two-fluid dust particles, will be discarded ***'
 endif

 allocate(xyzh(4,np),itype(np),dustfrac(ndusttypes,np),grainsize(ndusttypes),tmp(np))
 allocate(xyzmh_ptmass(5,nptmass))
 itype = 1

 ! extract info from real header
 call extract('massoftype',massoftype(1:ntypes),hdr,ierr)
 call extract('hfact',hfact,hdr,ierr)
 call extract('grainsize',grainsize(1:ndusttypes),hdr,ierr)
 call extract('graindens',graindens,hdr,ierr)
 !print*,' hfact = ',hfact
 !print*,' massoftype = ',massoftype(1:ntypes)

 call extract('umass',umass,hdr,ierr)
 call extract('utime',utime,hdr,ierr)
 call extract('udist',udist,hdr,ierr)

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
                read(iunit, iostat=ierr) tag
                write(*,"(1x,a)",advance='no') trim(tag)
                matched = .true.
                if (i==i_real .or. i==i_real8) then
                   select case(trim(tag))
                   case('x')
                      read(iunit,iostat=ierr) xyzh(1,1:np)
                   case('y')
                      read(iunit,iostat=ierr) xyzh(2,1:np)
                   case('z')
                      read(iunit,iostat=ierr) xyzh(3,1:np)
                   case('h')
                      read(iunit,iostat=ierr) xyzh(4,1:np)
                      got_h = .true.
                   case('dustfrac')
                      ngrains = ngrains + 1
                      read(iunit,iostat=ierr) dustfrac(ngrains,1:np)
                   case default
                      matched = .false.
                      read(iunit,iostat=ierr)
                   end select
                elseif (i==i_real4) then
                   select case(trim(tag))
                   case('h')
                      read(iunit,iostat=ierr) tmp(1:np)
                      xyzh(4,1:np) = real(tmp(1:np),kind=db)
                      got_h = .true.
                   case('dustfrac')
                      ngrains = ngrains + 1
                      read(iunit,iostat=ierr) tmp(1:np)
                      dustfrac(ngrains,1:np) = real(tmp(1:np),kind=db)
                      got_dustfrac = .true.
                   case default
                      matched = .false.
                      read(iunit,iostat=ierr)
                   end select
                elseif (i==i_int1) then
                   select case(trim(tag))
                   case('iphase')
                      matched = .true.
                      read(iunit,iostat=ierr) itype(1:np)
                   case default
                      read(iunit,iostat=ierr)
                   end select
                else
                   read(iunit,iostat=ierr)
                endif
                if (matched) then
                   write(*,"(a)") '->'//trim(tag)
                else
                   write(*,"(a)")
                endif
             elseif (j==1 .and. number8(j)==nptmass) then
                read(iunit,iostat=ierr) tag
                matched = .true.
                if (i==i_real .or. i==i_real8) then
                   ipos = 0
                   select case(trim(tag))
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
                   case default
                      matched = .false.
                   end select
                   if (ipos > 0) then
                      read(iunit,iostat=ierr) xyzmh_ptmass(ipos,1:nptmass)
                   endif
                else
                   matched = .false.
                   read(iunit,iostat=ierr)
                endif
                if (matched) then
                   write(*,"(a)") '->',trim(tag)
                else
                   write(*,"(a)")
                endif
             else
                read(iunit, iostat=ierr) tag ! tag
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
    write(*,"(a,i8,a)") ' Using ',ncells,' particles from Phantom file'
 else
    ncells = 0
    write(*,*) ' ERROR reading h from file'
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
  integer(kind=1), dimension(np), intent(in) :: iphase
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

 if (allocated(etoile)) deallocate(etoile)
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
