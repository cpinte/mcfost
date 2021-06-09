module read_phantom

  use parametres
  use dump_utils
  use messages
  use constantes
  use mess_up_SPH

  implicit none

  contains

subroutine read_phantom_bin_files(iunit,n_files,filenames,x,y,z,h,vx,vy,vz,particle_id,massgas,massdust,&
      rhogas,rhodust,Tgas,extra_heating,ndusttypes,SPH_grainsizes,mask,n_SPH,ierr)

 integer,               intent(in) :: iunit, n_files
 character(len=*),dimension(n_files), intent(in) :: filenames
 real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,vx,vy,vz,rhogas,Tgas,massgas,SPH_grainsizes
 integer,  intent(out), dimension(:),   allocatable :: particle_id
 real(dp), intent(out), dimension(:,:), allocatable :: rhodust,massdust
 logical, dimension(:), allocatable, intent(out) :: mask
 real, intent(out), dimension(:), allocatable :: extra_heating
 integer, intent(out) :: ndusttypes,n_SPH,ierr

 integer, parameter :: maxarraylengths = 12
 integer, parameter :: nsinkproperties = 17
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,nblockarrays,number,idust
 character(len=lentag) :: tag
 character(len=lenid)  :: fileid
 integer :: np,ntypes,nptmass,ipos,ngrains,dustfluidtype,ndudt
 integer, parameter :: maxtypes = 6
 integer, parameter :: maxfiles = 3
 integer, parameter :: maxinblock = 128
 integer, allocatable, dimension(:) :: npartoftype
 real(dp), allocatable, dimension(:,:) :: massoftype !(maxfiles,maxtypes)
 real(dp) :: hfact,umass,utime,udist
 integer(kind=1), allocatable, dimension(:) :: itype, ifiles
 real(4),  allocatable, dimension(:) :: tmp
 real(dp), allocatable, dimension(:) :: grainsize, graindens
 real(dp), allocatable, dimension(:) :: dudt, tmp_dp
 real(dp), allocatable, dimension(:,:) :: xyzh,xyzmh_ptmass,vxyz_ptmass,dustfrac,vxyzu
 type(dump_h) :: hdr
 logical :: got_h,got_dustfrac,got_itype,tagged,matched

 integer :: ifile, np0, ntypes0, np_tot, ntypes_tot, ntypes_max, ndustsmall, ndustlarge

 ! We first read the number of particules in each phantom file
 np_tot = 0
 ntypes_tot = 0
 ntypes_max = 0
 do ifile=1, n_files
    write(*,*) "---- Reading header file #", ifile
    !call read_phantom_bin_header()
    ! open file for read
    call open_dumpfile_r(iunit,filenames(ifile),fileid,ierr,requiretags=.true.)
    if (ierr /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - '//trim(get_error_text(ierr))//' ***'
       return
    endif

    if (fileid(2:2)=='T') then
       tagged = .true.
    else
       tagged = .false.
    endif
    call read_header(iunit,hdr,tagged,ierr)
    if (.not.tagged) then
       write(*,"(/,a,/)") ' *** ERROR - Phantom dump too old to be read by MCFOST ***'
       ierr = 1000
       return
    endif

    call extract('nparttot',np,hdr,ierr)
    call extract('ndusttypes',ndusttypes,hdr,ierr,default=0)
    if (ierr /= 0) then
       ! ndusttypes is for pre-largegrain multigrain headers
       call extract('ndustsmall',ndustsmall,hdr,ierr,default=0)
       call extract('ndustlarge',ndustlarge,hdr,ierr,default=0)
       ! ndusttype must be the same for all files : todo : add a test
       ndusttypes = ndustsmall + ndustlarge

       ! If ndusttypes, ndustlarge and ndustsmall are all missing, manually count grains
       if (ndusttypes==0 .and. ierr/=0) then
             ! For older files where ndusttypes is not output to the header
          idust = 0
          do i = 1,maxinblock
             if (hdr%realtags(i)=='grainsize') idust = idust + 1
          enddo
          write(*,"(a)")    ' Warning! Could not find ndusttypes in header'
          write(*,"(a,I4)") '          ...counting grainsize arrays...ndusttypes =',idust
          ndusttypes = idust
       endif

    endif
    call extract('ntypes',ntypes,hdr,ierr)

    np_tot = np_tot + np
    ntypes_tot = ntypes_tot + ntypes
    ntypes_max = max(ntypes_max,ntypes)

    call free_header(hdr, ierr)

    close(iunit)
    write(*,*) "---- Done"
 enddo ! ifile

 allocate(xyzh(4,np_tot),itype(np_tot),vxyzu(4,np_tot))
 allocate(dustfrac(ndusttypes,np_tot),grainsize(ndusttypes),graindens(ndusttypes))
 allocate(dudt(np_tot),ifiles(np_tot),massoftype(n_files,ntypes_max),npartoftype(ntypes_tot))


 np0 = 0
 ntypes0 = 0
 do ifile=1, n_files
    write(*,*) "---- Reading data file #", ifile

    ! open file for read
    call open_dumpfile_r(iunit,filenames(ifile),fileid,ierr,requiretags=.true.)
    if (ierr /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - '//trim(get_error_text(ierr))//' ***'
       return
    endif
    print "(1x,a)",trim(fileid)

    if (fileid(2:2)=='T') then
       tagged = .true.
    else
       tagged = .false.
    endif
    call read_header(iunit,hdr,tagged,ierr)
    if (.not.tagged) then
       write(*,"(/,a,/)") ' *** ERROR - Phantom dump too old to be read by MCFOST ***'
       ierr = 1000
       return
    endif

    ! get nblocks
    call extract('nblocks',nblocks,hdr,ierr,default=1)
    call extract('nparttot',np,hdr,ierr)
    call extract('ntypes',ntypes,hdr,ierr)
    call extract('npartoftype',npartoftype(ntypes0+1:ntypes0+ntypes),hdr,ierr)
    call extract('ndusttypes',ndusttypes,hdr,ierr,default=0)
    if (ierr /= 0) then
       ! ndusttypes is for pre-largegrain multigrain headers
       call extract('ndustsmall',ndustsmall,hdr,ierr,default=0)
       if (ierr /= 0) call warning("could not read ndustsmall")
       call extract('ndustlarge',ndustlarge,hdr,ierr,default=0)

       ! ndusttype must be the same for all files : todo : add a test
       ndusttypes = ndustsmall + ndustlarge

       ! If ndusttypes, ndustlarge and ndustsmall are all missing, manually count grains
       if (ndusttypes==0 .and. ierr/=0) then
          ! For older files where ndusttypes is not output to the header
          idust = 0
          do i = 1,maxinblock
             if (hdr%realtags(i)=='grainsize') idust = idust + 1
          enddo
          write(*,"(a)")    ' Warning! Could not find ndusttypes in header'
          write(*,"(a,I4)") '          ...counting grainsize arrays...ndusttypes =',idust
          ndusttypes = idust
      endif
    endif
    call extract('ntypes',ntypes,hdr,ierr)
    call extract('nptmass',nptmass,hdr,ierr,default=0)
    !call extract('isink',isink,hdr,ierr,default=0)

    write(*,*) ' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
    !write(*,*) ' npartoftype = ',npartoftype(ntypes0+1:ntypes0+ntypes)
    write(*,*) ' nptmass = ', nptmass

    allocate(tmp(np), tmp_dp(np))

    if (npartoftype(2) > 0) then
       dustfluidtype = 2
       write(*,"(/,a,/)") ' *** WARNING: Phantom dump contains two-fluid dust particles, may be discarded ***'
    else
       dustfluidtype = 1
    endif

    ! extract info from real header
    call extract('massoftype',massoftype(ifile,1:ntypes),hdr,ierr)
    call extract('hfact',hfact,hdr,ierr)
    if (ndusttypes > 0) then
       call extract('grainsize',grainsize(1:ndusttypes),hdr,ierr) ! code units here
       call extract('graindens',graindens(1:ndusttypes),hdr,ierr)
    endif

    call extract('umass',umass,hdr,ierr)
    call extract('utime',utime,hdr,ierr)
    call extract('udist',udist,hdr,ierr)

    read (iunit, iostat=ierr) number
    if (ierr /= 0) return
    narraylengths = number/nblocks

    got_h = .false.
    got_dustfrac = .false.
    got_itype = .false.
    ! skip each block that is too small
    nblockarrays = narraylengths*nblocks

    ndudt = 0
    ngrains = 0
    do iblock = 1,nblocks
       call read_block_header(narraylengths,number8,nums,iunit,ierr)
       if (ierr /= 0) call error('Reading block header')
       do j=1,narraylengths
          do i=1,ndatatypes
             !write(*,*) ' data type ',i,' arrays = ',nums(i,j)
             do k=1,nums(i,j)
                !write(*,*) "k=", k, np0, np

                if (j==1 .and. number8(j)==np) then
                   read(iunit, iostat=ierr) tag
                   if (ierr /= 0) call error('Reading tag')
                   !write(*,"(1x,a)",advance='no') trim(tag)

                   matched = .true.
                   if (i==i_real .or. i==i_real8) then
                      select case(trim(tag))
                      case('x')
                         read(iunit,iostat=ierr) tmp_dp ; xyzh(1,np0+1:np0+np) = tmp_dp
                      case('y')
                         read(iunit,iostat=ierr) tmp_dp ; xyzh(2,np0+1:np0+np) = tmp_dp
                      case('z')
                         read(iunit,iostat=ierr) tmp_dp ; xyzh(3,np0+1:np0+np) = tmp_dp
                      case('h')
                         read(iunit,iostat=ierr) tmp_dp ; xyzh(4,np0+1:np0+np) = tmp_dp
                         got_h = .true.
                      case('vx')
                         read(iunit,iostat=ierr) tmp_dp ; vxyzu(1,np0+1:np0+np) = tmp_dp
                      case('vy')
                         read(iunit,iostat=ierr) tmp_dp ; vxyzu(2,np0+1:np0+np) = tmp_dp
                      case('vz')
                         read(iunit,iostat=ierr) tmp_dp ; vxyzu(3,np0+1:np0+np) = tmp_dp
                      case('u')
                         read(iunit,iostat=ierr) tmp_dp ; vxyzu(4,np0+1:np0+np) = tmp_dp
                      case('dustfrac')
                         ngrains = ngrains + 1
                         if (ngrains > ndusttypes) then
                            write(*,*) "ERROR ngrains > ndusttypes:", ngrains, ndusttypes
                            ierr = 1 ;
                            return
                         endif
                         read(iunit,iostat=ierr) dustfrac(ngrains,1:np)
                         got_dustfrac = .true.
                      case default
                         matched = .false.
                         read(iunit,iostat=ierr)
                      end select
                   elseif (i==i_real4) then
                      select case(trim(tag))
                      case('h')
                         read(iunit,iostat=ierr) tmp(1:np)
                         xyzh(4,np0+1:np0+np) = real(tmp(1:np),kind=dp)
                         got_h = .true.
                      case('dustfrac')
                         ngrains = ngrains + 1
                         read(iunit,iostat=ierr) tmp(1:np)
                         dustfrac(ngrains,np0+1:np0+np) = real(tmp(1:np),kind=dp)
                         got_dustfrac = .true.
                      case('luminosity')
                         read(iunit,iostat=ierr) tmp(1:np)
                         dudt(np0+1:np0+np) = real(tmp(1:np),kind=dp)
                         ndudt = np
                      case default
                         matched = .false.
                         read(iunit,iostat=ierr)
                      end select
                   elseif (i==i_int1) then
                      select case(trim(tag))
                      case('itype')
                         matched = .true.
                         got_itype = .true.
                         read(iunit,iostat=ierr) itype(np0+1:np0+np)
                         itype(np0+1:np0+np) = itype(np0+1:np0+np) + ntypes0 ! shifting types for succesive files
                         write(*,*) "MAX", maxval(itype(np0+1:np0+np)), ntypes0
                      case default
                         read(iunit,iostat=ierr)
                      end select
                   else
                      read(iunit,iostat=ierr)
                   endif
                   if (ierr /= 0) call error("Error reading tag: "//trim(tag))
                !elseif (j==1 .and. number8(j)==nptmass) then
                elseif (j==2) then ! HACK : what is j exactly anyway ? and why would we need to test for j==1
                   nptmass = number8(j) ! HACK
                   if (.not.allocated(xyzmh_ptmass)) then
                      allocate(xyzmh_ptmass(nsinkproperties,nptmass), vxyz_ptmass(3,nptmass)) !HACK
                      xyzmh_ptmass(:,:) = 0.
                      vxyz_ptmass(:,:) = 0.
                   endif

                   read(iunit,iostat=ierr) tag
                   matched = .true.
                   if (i==i_real .or. i==i_real8) then
                      select case(trim(tag))
                      case('x')
                         read(iunit,iostat=ierr) xyzmh_ptmass(1,1:nptmass)
                      case('y')
                         read(iunit,iostat=ierr) xyzmh_ptmass(2,1:nptmass)
                      case('z')
                         read(iunit,iostat=ierr) xyzmh_ptmass(3,1:nptmass)
                      case('m')
                         read(iunit,iostat=ierr) xyzmh_ptmass(4,1:nptmass)
                      case('h')
                         read(iunit,iostat=ierr) xyzmh_ptmass(5,1:nptmass)
                      case('mdotav')
                         read(iunit,iostat=ierr) xyzmh_ptmass(16,1:nptmass)
                      case('vx')
                         read(iunit,iostat=ierr) vxyz_ptmass(1,1:nptmass)
                      case('vy')
                         read(iunit,iostat=ierr) vxyz_ptmass(2,1:nptmass)
                      case('vz')
                         read(iunit,iostat=ierr) vxyz_ptmass(3,1:nptmass)
                      case default
                         matched = .false.
                      end select
                   else
                      matched = .false.
                      read(iunit,iostat=ierr)
                   endif
                   !     if (matched) then
                   !        write(*,"(a)") '->',trim(tag)
                   !     else
                   !        write(*,"(a)")
                   !     endif
                else
                   read(iunit, iostat=ierr) tag ! tag
                   !print*,tagarr(1)
                   read(iunit, iostat=ierr) ! array
                endif
             enddo
          enddo
       enddo
    enddo ! block

    deallocate(tmp, tmp_dp)
    call free_header(hdr, ierr)
    write(*,*) "---- Done"
    close(iunit)

    ifiles(np0+1:np0+np) = ifile ! index of the file

    if (.not. got_itype) then
       itype(np0+1:np0+np) =  1
    endif

    np0 = np0 + np
    ntypes0 = ntypes0 + ntypes
 enddo ! ifile

 if ((ndusttypes > 0) .and. .not. got_dustfrac) then
    dustfrac = 0.
 endif

 if (ndudt == 0) then
    dudt = 0.
 endif

 write(*,*) "Found", nptmass, "point masses in the phantom file"

 if (got_h) then
    call modify_dump(np, nptmass, xyzh, vxyzu, xyzmh_ptmass, udist, mask)

    call phantom_2_mcfost(np_tot,nptmass,ntypes_max,ndusttypes,n_files,dustfluidtype,xyzh,&
         vxyzu,itype,grainsize,dustfrac,massoftype,xyzmh_ptmass,vxyz_ptmass,&
         hfact,umass,utime,udist,graindens,ndudt,dudt,ifiles, &
         n_SPH,x,y,z,h,vx,vy,vz,particle_id, &
         SPH_grainsizes,massgas,massdust,rhogas,rhodust,Tgas,extra_heating)
    write(*,"(a,i8,a)") ' Using ',n_SPH,' particles from Phantom file'
 else
    n_SPH = 0
    write(*,*) "ERROR reading h from file"
 endif

 write(*,*) "Phantom dump file processed ok"
 deallocate(xyzh,itype,vxyzu)
 if (allocated(xyzmh_ptmass)) deallocate(xyzmh_ptmass,vxyz_ptmass)

end subroutine read_phantom_bin_files

!*************************************************************************

subroutine read_phantom_hdf_files(iunit,n_files, filenames, x,y,z,h,vx,vy,vz,  &
                                  particle_id,massgas,massdust,rhogas,rhodust, &
                                  Tgas,extra_heating,ndusttypes,SPH_grainsizes,     &
                                  mask,n_SPH,ierr)

 use utils_hdf5, only:open_hdf5file,    &
                      close_hdf5file,   &
                      open_hdf5group,   &
                      close_hdf5group,  &
                      read_from_hdf5,   &
                      HID_T

 integer, intent(in) :: iunit, n_files
 character(len=*),dimension(n_files), intent(in) :: filenames
 real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,        &
                                                       vx,vy,vz,       &
                                                       rhogas,massgas, &
                                                       Tgas,SPH_grainsizes
 integer,  intent(out), dimension(:),   allocatable :: particle_id
 real(dp), intent(out), dimension(:,:), allocatable :: rhodust,massdust
 logical, dimension(:), allocatable, intent(out) :: mask
 real, intent(out), dimension(:), allocatable :: extra_heating
 integer, intent(out) :: ndusttypes,n_SPH,ierr

 character(len=200) :: filename

 logical :: got_dustfrac,got_itype

 integer :: ifile, np0, ntypes0, np_tot, ntypes_tot, ntypes_max
 integer :: np,ntypes,nptmass,dustfluidtype,ndudt
 integer :: error,ndustsmall,ndustlarge

 integer, parameter :: maxtypes = 100
 integer, parameter :: nsinkproperties = 17

 integer(kind=1), allocatable, dimension(:) :: itype, ifiles
 real(4),  allocatable, dimension(:)   :: tmp, tmp_header
 real(dp), allocatable, dimension(:)   :: dudt, tmp_dp
 real(dp), allocatable, dimension(:)   :: grainsize, graindens
 real(dp), allocatable, dimension(:,:) :: xyzh,xyzmh_ptmass,vxyz_ptmass,dustfrac,vxyzu
 integer, allocatable, dimension(:) :: npartoftype
 real(dp), allocatable, dimension(:,:) :: massoftype !(maxfiles,maxtypes)
 real(dp) :: hfact,umass,utime,udist

 integer(HID_T) :: hdf5_file_id
 integer(HID_T) :: hdf5_group_id

 logical :: got

  error = 0
  np_tot = 0
  ntypes_tot = 0
  ntypes_max = 0

  ! Read file headers to get particle numbers
  file_headers: do ifile=1, n_files

    write(*,*) "---- Reading header file #", ifile

    filename = trim(filenames(ifile))

    ! open file
    call open_hdf5file(filename,hdf5_file_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot open Phantom HDF file ***'
    endif

    ! open header group
    call open_hdf5group(hdf5_file_id,'header',hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot open Phantom HDF header group ***'
    endif

    ! read header values
    call read_from_hdf5(np,'nparttot',hdf5_group_id,got,error)
    call read_from_hdf5(ntypes,'ntypes',hdf5_group_id,got,error)
    call read_from_hdf5(ndusttypes,'ndusttypes',hdf5_group_id,got,error)
    if (.not. got) then
       ! ndusttypes is for pre-largegrain multigrain headers
       call read_from_hdf5(ndustsmall,'ndustsmall',hdf5_group_id,got,error)
       call read_from_hdf5(ndustlarge,'ndustlarge',hdf5_group_id,got,error)
       ! ndusttype must be the same for all files : todo : add a test
       ndusttypes = ndustsmall + ndustlarge
    endif

    ! close the header group
    call close_hdf5group(hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot close Phantom HDF header group ***'
    endif

    ! close file
    call close_hdf5file(hdf5_file_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot close Phantom HDF file ***'
    endif

    np_tot = np_tot + np
    ntypes_tot = ntypes_tot + ntypes
    ntypes_max = max(ntypes_max,ntypes)

  enddo file_headers

  ! Allocate arrays for SPH data
  allocate(xyzh(4,np_tot),                        &
           itype(np_tot),                         &
           vxyzu(4,np_tot),                       &
           dustfrac(ndusttypes,np_tot),           &
           grainsize(ndusttypes),                 &
           graindens(ndusttypes),                 &
           dudt(np_tot),                          &
           ifiles(np_tot),                        &
           massoftype(n_files,ntypes_max),        &
           npartoftype(ntypes_tot))

  ! Read file data
  np0 = 0
  ntypes0 = 0
  file_data: do ifile=1, n_files

    ! open file
    filename = trim(filenames(ifile))
    call open_hdf5file(filename,hdf5_file_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot open Phantom HDF file ***'
    endif

    ! open header group
    call open_hdf5group(hdf5_file_id,'header',hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot open Phantom HDF header group ***'
    endif

    ! read from header
    call read_from_hdf5(np,'nparttot',hdf5_group_id,got,error)
    call read_from_hdf5(ntypes,'ntypes',hdf5_group_id,got,error)
    call read_from_hdf5(npartoftype,'npartoftype',hdf5_group_id,got,error)
    call read_from_hdf5(nptmass,'nptmass',hdf5_group_id,got,error)
    call read_from_hdf5(ndusttypes,'ndusttypes',hdf5_group_id,got,error)
    if (.not. got) then
       ! ndusttypes is for pre-largegrain multigrain headers
       call read_from_hdf5(ndustsmall,'ndustsmall',hdf5_group_id,got,error)
       call read_from_hdf5(ndustlarge,'ndustlarge',hdf5_group_id,got,error)
       ! ndusttype must be the same for all files : todo : add a test
       ndusttypes = ndustsmall + ndustlarge
    endif
    call read_from_hdf5(massoftype(ifile,1:ntypes),'massoftype',hdf5_group_id,got,error)
    call read_from_hdf5(hfact,'hfact',hdf5_group_id,got,error)
    if (ndusttypes > 0) then
       allocate(tmp_header(np))
       call read_from_hdf5(tmp_header,'grainsize',hdf5_group_id,got,error)
       grainsize(1:ndusttypes) = tmp_header(1:ndusttypes)
       call read_from_hdf5(tmp_header,'graindens',hdf5_group_id,got,error)
       graindens(1:ndusttypes) = tmp_header(1:ndusttypes)
       deallocate(tmp_header)
    endif
    call read_from_hdf5(umass,'umass',hdf5_group_id,got,error)
    call read_from_hdf5(utime,'utime',hdf5_group_id,got,error)
    call read_from_hdf5(udist,'udist',hdf5_group_id,got,error)

    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot read Phantom HDF header group ***'
    endif

    ! close the header group
    call close_hdf5group(hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot close Phantom HDF header group ***'
    endif

    write(*,*) ' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
    !write(*,*) ' npartoftype = ',npartoftype(ntypes0+1:ntypes0+ntypes)
    write(*,*) ' nptmass = ', nptmass

    allocate(tmp(np), tmp_dp(np))

    if (npartoftype(2) > 0) then
       dustfluidtype = 2
       write(*,"(/,a,/)") ' *** WARNING: Phantom dump contains two-fluid dust particles, may be discarded ***'
    else
       dustfluidtype = 1
    endif

    ! open particles group
    call open_hdf5group(hdf5_file_id,'particles',hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot open Phantom HDF particles group ***'
    endif

    got_dustfrac = .false.
    got_itype = .false.
    ndudt = 0

    ! TODO: read particle arrays
    call read_from_hdf5(itype(np0+1:np0+np),'itype',hdf5_group_id,got,error)
    if (got) got_itype = .true.
    call read_from_hdf5(xyzh(1:3,np0+1:np0+np),'xyz',hdf5_group_id,got,error)
    call read_from_hdf5(xyzh(4,np0+1:np0+np),'h',hdf5_group_id,got,error)
    call read_from_hdf5(vxyzu(1:3,np0+1:np0+np),'vxyz',hdf5_group_id,got,error)
    call read_from_hdf5(dustfrac(1:ndusttypes,np0+1:np0+np),'dustfrac',hdf5_group_id,got,error)
    if (got) got_dustfrac = .true.

    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot read Phantom HDF particles group ***'
    endif

    ! close the particles group
    call close_hdf5group(hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot close Phantom HDF particles group ***'
    endif

    ! open sinks group
    call open_hdf5group(hdf5_file_id,'sinks',hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot open Phantom HDF sinks group ***'
    endif

    ! read sinks
    allocate(xyzmh_ptmass(nsinkproperties,nptmass), vxyz_ptmass(3,nptmass))
    call read_from_hdf5(xyzmh_ptmass(1:3,:),'xyz',hdf5_group_id,got,error)
    call read_from_hdf5(xyzmh_ptmass(4,:),'m',hdf5_group_id,got,error)
    call read_from_hdf5(xyzmh_ptmass(5,:),'h',hdf5_group_id,got,error)
    call read_from_hdf5(vxyz_ptmass(:,:),'vxyz',hdf5_group_id,got,error)

    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot read Phantom HDF sinks group ***'
    endif

    call read_from_hdf5(xyzmh_ptmass(16,:),'mdotav',hdf5_group_id,got,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** WARNMING : mdotav not present'
       error =0
    endif

    ! close the sinks group
    call close_hdf5group(hdf5_group_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot close Phantom HDF sinks group ***'
    endif

    ! close file
    call close_hdf5file(hdf5_file_id,error)
    if (error /= 0) then
       write(*,"(/,a,/)") ' *** ERROR - cannot close Phantom HDF file ***'
    endif

    ifiles(np0+1:np0+np) = ifile ! index of the file

    if (.not. got_itype) then
       itype(np0+1:np0+np) =  1
    endif

    np0 = np0 + np
    ntypes0 = ntypes0 + ntypes

  enddo file_data

  if ((ndusttypes > 0) .and. .not. got_dustfrac) then
     dustfrac = 0.
  endif

  if (ndudt == 0) then
     dudt = 0.
  endif

  write(*,*) "Found", nptmass, "point masses in the phantom file"

  call modify_dump(np, nptmass, xyzh, vxyzu, xyzmh_ptmass, udist, mask)

  call phantom_2_mcfost(np_tot,nptmass,ntypes_max,ndusttypes,n_files,      &
                        dustfluidtype,xyzh,vxyzu,itype,grainsize,dustfrac, &
                        massoftype,xyzmh_ptmass,vxyz_ptmass,hfact,umass,       &
                        utime,udist,graindens,ndudt,dudt,ifiles,   &
                        n_SPH,x,y,z,h,vx,vy,vz,particle_id,SPH_grainsizes,     &
                        massgas,massdust,rhogas,rhodust,Tgas,extra_heating)

  write(*,"(a,i8,a)") ' Using ',n_SPH,' particles from Phantom file'

  write(*,*) "Phantom dump file processed ok"
  ierr = 0
  deallocate(xyzh,itype,vxyzu)
  if (allocated(xyzmh_ptmass)) deallocate(xyzmh_ptmass,vxyz_ptmass)

end subroutine read_phantom_hdf_files

!*************************************************************************

subroutine modify_dump(np, nptmass, xyzh, vxyzu, xyzmh_ptmass, udist, mask)

  integer, intent(in) :: np, nptmass
  real(kind=dp), dimension(:,:), intent(inout) :: xyzh, vxyzu
  real(kind=dp), dimension(:,:), intent(inout) :: xyzmh_ptmass
  real(kind=dp), intent(in) :: udist

  real(kind=dp), dimension(3) :: centre

  logical, dimension(:), allocatable, intent(out) :: mask

  integer :: i

  ! Modifying SPH dump
  if (ldelete_Hill_sphere) then
     allocate(mask(np))
     call mask_Hill_sphere(np, nptmass, xyzh, xyzmh_ptmass,udist, mask)
  endif
  if (lrandomize_azimuth)     call randomize_azimuth(np, xyzh, vxyzu, mask)
  if (lrandomize_gap)         call randomize_gap(np, nptmass, xyzh, vxyzu, xyzmh_ptmass,udist, gap_factor, .true.)
  if (lrandomize_outside_gap) call randomize_gap(np, nptmass, xyzh, vxyzu, xyzmh_ptmass,udist, gap_factor, .false.)

  if (lcentre_on_sink) then
     write(*,*) "Recentering on sink #", isink_centre, nptmass, np
     centre = xyzmh_ptmass(1:3,isink_centre)
     do i=1,nptmass
        xyzmh_ptmass(1:3,i) = xyzmh_ptmass(1:3,i) - centre(:)
     enddo
     do i=1, np
        xyzh(1:3,i) = xyzh(1:3,i) - centre(:)
     enddo
  endif

  return

end subroutine modify_dump

!*************************************************************************


subroutine phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,n_files,dustfluidtype,xyzh, &
     vxyzu,iphase,grainsize,dustfrac,massoftype,xyzmh_ptmass,vxyz_ptmass,hfact,umass, &
     utime, udist,graindens,ndudt,dudt,ifiles, &
     n_SPH,x,y,z,h,vx,vy,vz,particle_id, &
     SPH_grainsizes,massgas,massdust,rhogas,rhodust,Tgas,extra_heating,T_to_u)

  ! Convert phantom quantities & units to mcfost quantities & units
  ! x,y,z are in au
  ! rhodust & rhogas are in g/cm3
  ! extra_heating is in W

  use constantes, only : au_to_cm,Msun_to_g,erg_to_J,m_to_cm,Lsun,cm_to_mum,deg_to_rad,Ggrav,cm_to_m,g_to_kg
  use parametres, only : ldudt_implicit,ufac_implicit, lplanet_az, planet_az, lfix_star, RT_az_min, RT_az_max, RT_n_az
  use parametres, only : lscale_length_units,scale_length_units_factor,lscale_mass_units,scale_mass_units_factor

  integer, intent(in) :: np,nptmass,ntypes,ndusttypes,dustfluidtype, n_files
  real(dp), dimension(4,np), intent(in) :: xyzh,vxyzu
  integer(kind=1), dimension(np), intent(in) :: iphase, ifiles
  real(dp), dimension(ndusttypes,np), intent(in) :: dustfrac
  real(dp), dimension(ndusttypes),    intent(in) :: grainsize ! code units
  real(dp), dimension(ndusttypes),    intent(in) :: graindens
  real(dp), dimension(n_files,ntypes), intent(in) :: massoftype
  real(dp), intent(in) :: hfact,umass,utime,udist
  real(dp), dimension(:,:), intent(in) :: xyzmh_ptmass, vxyz_ptmass
  integer, intent(in) :: ndudt
  real(dp), dimension(:), intent(in) :: dudt

  ! MC
  integer, intent(out) :: n_SPH
  real(dp), dimension(:),   allocatable, intent(out) :: x,y,z,h,vx,vy,vz,rhogas,massgas,Tgas ! massgas [Msun]
  integer, dimension(:),    allocatable, intent(out) :: particle_id
  real(dp), dimension(:,:), allocatable, intent(out) :: rhodust,massdust
  real(dp), dimension(:), allocatable, intent(out) :: SPH_grainsizes ! mum
  real, dimension(:), allocatable, intent(out) :: extra_heating
  real(dp), intent(in), optional :: T_to_u


  integer  :: i,j,k,itypei,alloc_status,i_etoile, n_etoiles_old, ifile
  real(dp) :: xi,yi,zi,hi,vxi,vyi,vzi,rhogasi,rhodusti,gasfraci,dustfraci,totlum,qtermi
  real(dp) :: udist_scaled, umass_scaled, utime_scaled,udens,uerg_per_s,uWatt,ulength_au,usolarmass,uvelocity
  real(dp) :: vphi, vr, phi, cos_phi, sin_phi, r_cyl, r_cyl2, r_sph, G_phantom

  logical :: use_dust_particles = .false. ! 2-fluid: choose to use dust


  ! We check the units by recomputing G
  G_phantom = udist**3 / (utime**2 * umass)
  if (abs(G_phantom - Ggrav*1e3) > 1e-2 * G_phantom) call error("Phantom units are not consistent")

  udist_scaled = udist
  umass_scaled = umass
  if (lscale_length_units) then
     write(*,*) 'Lengths are rescaled by ', scale_length_units_factor
     udist_scaled = udist * scale_length_units_factor
  else
     scale_length_units_factor = 1.0
  endif
  if (lscale_mass_units) then
     write(*,*) 'Mass are rescaled by ', scale_mass_units_factor
     umass_scaled = umass * scale_mass_units_factor
  else
     scale_mass_units_factor = 1.0
  endif
  utime_scaled = sqrt(udist_scaled**3/(G_phantom*umass_scaled))

  udens = umass_scaled/udist_scaled**3
  uerg_per_s = umass_scaled*udist_scaled**2/utime_scaled**3
  uWatt = uerg_per_s * erg_to_J
  ulength_au = udist_scaled / (au_to_cm)
  uvelocity =  udist_scaled / (m_to_cm) / utime_scaled ! m/s
  usolarmass = umass_scaled / Msun_to_g

 if (dustfluidtype == 1) then
    ! 1-fluid: always use gas particles for Voronoi mesh
    use_dust_particles = .false.
 endif

 write(*,*) ''
 if (use_dust_particles) then
    write(*,*) '*** Using dust particles for Voronoi mesh ***'
 else
    write(*,*) '*** Using gas particles for Voronoi mesh ***', n_files
 endif

 ! convert to dust and gas density
 j = 0
 do i=1,np
    if (xyzh(4,i) > 0.) then
       if (.not. use_dust_particles .and. abs(iphase(i))==1)  j = j + 1
       if (      use_dust_particles .and. abs(iphase(i))==2)  j = j + 1
    endif
 enddo
 n_SPH = j

 ! TODO : use mcfost quantities directly rather that these intermediate variables
 ! Voronoi()%x  densite_gaz & densite_pous
 alloc_status = 0
 allocate(rhodust(ndusttypes,n_SPH),massdust(ndusttypes,n_SPH),SPH_grainsizes(ndusttypes),particle_id(n_SPH),&
      x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),massgas(n_SPH),rhogas(n_SPH),Tgas(n_SPH),stat=alloc_status)
 if (alloc_status /=0) then
    write(*,*) "Allocation error in phanton_2_mcfost"
    write(*,*) "Exiting"
 endif

 ! Dust grain sizes in microns
 SPH_grainsizes(:) = grainsize(:) * udist_scaled * cm_to_mum
 ! graindens * udens is in g/cm3

 if (lemission_mol .or. lemission_atom) then
    allocate(vx(n_SPH),vy(n_SPH),vz(n_SPH),stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error velocities in phanton_2_mcfost"
       write(*,*) "Exiting"
    endif
 endif

 j = 0
 do i=1,np
    ifile = ifiles(i)

    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)

    vxi = vxyzu(1,i)
    vyi = vxyzu(2,i)
    vzi = vxyzu(3,i)

    itypei = abs(iphase(i))
    if (hi > 0.) then
       if (use_dust_particles .and. dustfluidtype==2 .and. ndusttypes==1 .and. itypei==2) then
          j = j + 1
          particle_id(j) = i
          x(j) = xi * ulength_au
          y(j) = yi * ulength_au
          z(j) = zi * ulength_au
          h(j) = hi * ulength_au
          if (lemission_mol) then
             vx(j) = vxi * uvelocity
             vy(j) = vyi * uvelocity
             vz(j) = vzi * uvelocity
          endif
          rhodusti = massoftype(ifile,itypei) * (hfact/hi)**3  * udens ! g/cm**3
          gasfraci = dustfrac(1,i)
          rhodust(1,j) = rhodusti
          massdust(1,j) = massoftype(ifile,itypei) * usolarmass ! Msun
          rhogas(j) = gasfraci*rhodusti
          massgas(j) = gasfraci*massdust(1,j)
       elseif (.not. use_dust_particles .and. itypei==1) then
          j = j + 1
          particle_id(j) = i
          x(j) = xi * ulength_au
          y(j) = yi * ulength_au
          z(j) = zi * ulength_au
          h(j) = hi * ulength_au
          if (lemission_mol) then
             vx(j) = vxi * uvelocity
             vy(j) = vyi * uvelocity
             vz(j) = vzi * uvelocity
          endif
          rhogasi = massoftype(ifile,itypei) *(hfact/hi)**3  * udens ! g/cm**3
          dustfraci = sum(dustfrac(:,i))
          if (dustfluidtype==1) then
             rhogas(j) = rhogasi
          else
             rhogas(j) = (1 - dustfraci)*rhogasi
          endif
          massgas(j) =  massoftype(ifile,itypei) * usolarmass ! Msun
          do k=1,ndusttypes
             rhodust(k,j) = dustfrac(k,i)*rhogasi
             massdust(k,j) = dustfrac(k,i)*massgas(j)
          enddo
          if (lemission_atom) then
             !
             ! solve for the gas temperature from the thermal energy
             ! this should only be done if the temperature is NOT read from phantom
             ! also: mu here should ideally be consistent with the abundances used for
             ! the atomic transfer and the ionisation state, at the moment
             ! we just assume everything is fully ionised to compute Tgas
             !
             rhogasi = rhogasi*g_to_kg/cm_to_m**3
             Tgas(j) = get_temp_from_u(vxyzu(4,i)*uvelocity**2,rhogasi,mu=0.6_dp)
             if (mod(j,100000).eq.0) print*,i,Tgas(j),vxyzu(4,i),rhogasi/g_to_kg*cm_to_m**3,'g/cm**3'
          endif
       endif
    endif
 enddo
 n_SPH = j

 if (sum(massgas) == 0.) call error('Using old phantom dumpfile without gas-to-dust ratio on dust particles.', &
      msg2='Set use_dust_particles = .false. (read_phantom.f90) and compile and run again.')

 write(*,*) ''
 write(*,*) 'SPH gas mass is  ', real(sum(massgas)), 'Msun'
 write(*,*) 'SPH dust mass is ', real(sum(massdust)),'Msun'
 write(*,*) ''

 if (lemission_atom) write(*,*) ' max temperature is ',maxval(Tgas),' min temperature is ',minval(Tgas)

 lextra_heating = .false. ; ldudt_implicit = .false.
 if (.not.lno_internal_energy) then
    if (ndudt == np) then
       write(*,*) "Computing energy input"
       lextra_heating = .true.
       allocate(extra_heating(n_SPH), stat=alloc_status)
       if (alloc_status /=0) call error("Allocation error in phanton_2_mcfost")

       totlum = 0.
       do i=1,n_SPH
          qtermi = dudt(particle_id(i)) * uWatt
          totlum = totlum + qtermi
          extra_heating(i) = qtermi
       enddo

       write(*,*) "Total energy input = ",real(totlum),' W'
       write(*,*) "Total energy input = ",real(totlum/Lsun),' Lsun'
    endif

    if (present( T_to_u )) then
       ufac_implicit = T_to_u * uWatt
       ldudt_implicit = .true.
       write(*,*) "Implicit heating scheme"
    else
       ufac_implicit = 0.
       ldudt_implicit = .false.
    endif
 endif

 write(*,*) "# Sink particles:"
 n_etoiles_old = n_etoiles
 n_etoiles = 0
 do i=1,nptmass
    n_etoiles = n_etoiles + 1
    if (real(xyzmh_ptmass(4,i)) * scale_mass_units_factor > 0.013) then
       write(*,*) "Sink #", i, "xyz=", real(xyzmh_ptmass(1:3,i) * scale_length_units_factor), "au, M=", &
            real(xyzmh_ptmass(4,i) * scale_mass_units_factor), "Msun, Mdot=", &
            real(xyzmh_ptmass(16,i) * usolarmass / utime_scaled * year_to_s ), "Msun/yr"
    else
       write(*,*) "Sink #", i, "xyz=", real(xyzmh_ptmass(1:3,i) * scale_length_units_factor), "au, M=", &
            real(xyzmh_ptmass(4,i) * GxMsun/GxMjup * scale_mass_units_factor), "Mjup, Mdot=", &
            real(xyzmh_ptmass(16,i) * usolarmass / utime_scaled * year_to_s ), "Msun/yr"
    endif
    if (i>1) write(*,*)  "       distance=", real(norm2(xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,1))*scale_length_units_factor), "au"
 enddo


 if (lupdate_velocities) then
    if (lno_vz) vz(:) = 0.

    if (lno_vr) then
       do i=1, n_SPH
          phi = atan2(y(i),x(i)) ; cos_phi = cos(phi) ; sin_phi = sin(phi)
          vphi = - vx(i) * sin_phi + vy(i) * cos_phi
          ! vr = vx(i) * cos_phi + vy(i) * sin_phi
          ! vx = vr cos phi - vphi sin phi
          ! vy = vr sin phi + vphi cos phi
          vx(i) = - vphi * sin_phi !  vr = 0
          vy(i) = vphi * cos_phi
       enddo
    endif

    if (lVphi_kep) then
       do i=1, n_SPH
          r_cyl2 = x(i)**2 + y(i)**2
          r_cyl = sqrt(r_cyl2)
          cos_phi = x(i)/r_cyl ; sin_phi = y(i)/r_cyl
          r_sph = sqrt(r_cyl2 + z(i)**2)

          !U_r = (/ x(i)/r_cyl, y(i)/r_cyl /)
          !U_phi = (/ -y(i)/r_cyl, x(i)/r_cyl /)
          !v = (/vx(i), vy(i)/)

          vr = vx(i) * cos_phi   + vy(i) * sin_phi
          !vphi = - vx(i) * sin_phi + vy(i) * cos_phi

          ! Keplerian vphi
          vphi = sqrt(Ggrav * xyzmh_ptmass(4,1) * scale_mass_units_factor * Msun_to_kg  &
               * (r_cyl * AU_to_m)**2 /  (r_sph * AU_to_m)**3 )

          vx(i) = vr * cos_phi - vphi * sin_phi
          vy(i) = vr * sin_phi + vphi * cos_phi
       enddo
    endif
 endif


 if (lplanet_az) then
    if ((nptmass /= 2).and.(which_planet==0)) then
       call error("option -planet_az: you need to specify the sink particle with -planet option")
    endif
    if (nptmass == 2) which_planet=2
    if (which_planet > nptmass) call error("specified sink particle does not exist")

    RT_n_az = 1
    RT_az_min = planet_az + atan2(xyzmh_ptmass(2,which_planet) - xyzmh_ptmass(2,1), &
                                  xyzmh_ptmass(1,which_planet) - xyzmh_ptmass(1,1)) &
                            / deg_to_rad
    RT_az_max = RT_az_min
    write(*,*) "Moving sink particle #", which_planet, "to azimuth =", planet_az
    write(*,*) "WARNING: updating the azimuth to:", RT_az_min
 endif

 if (lfix_star) then
    write(*,*) ""
    write(*,*) "Stellar parameters will not be updated, only the star positions and masses"
    if (n_etoiles > n_etoiles_old) then
       write(*,*) "WARNING: sink with id >", n_etoiles_old, "will be ignored in the RT"
       n_etoiles = n_etoiles_old
    endif
    do i_etoile = 1, n_etoiles
       etoile(i_etoile)%x = xyzmh_ptmass(1,i_etoile) * ulength_au
       etoile(i_etoile)%y = xyzmh_ptmass(2,i_etoile) * ulength_au
       etoile(i_etoile)%z = xyzmh_ptmass(3,i_etoile) * ulength_au

       etoile(i_etoile)%vx = vxyz_ptmass(1,i_etoile) * uvelocity
       etoile(i_etoile)%vy = vxyz_ptmass(2,i_etoile) * uvelocity
       etoile(i_etoile)%vz = vxyz_ptmass(3,i_etoile) * uvelocity

       etoile(i_etoile)%M = xyzmh_ptmass(4,i_etoile) * usolarmass

       etoile(i_etoile)%Mdot = 0.
    enddo
 else
    write(*,*) ""
    write(*,*) "Updating the stellar properties:"
    write(*,*) "There are now", n_etoiles, "stars in the model"
    if (allocated(etoile)) deallocate(etoile)
    allocate(etoile(n_etoiles))

    etoile(:)%x = xyzmh_ptmass(1,:) * ulength_au
    etoile(:)%y = xyzmh_ptmass(2,:) * ulength_au
    etoile(:)%z = xyzmh_ptmass(3,:) * ulength_au

    etoile(:)%vx = vxyz_ptmass(1,:) * uvelocity
    etoile(:)%vy = vxyz_ptmass(2,:) * uvelocity
    etoile(:)%vz = vxyz_ptmass(3,:) * uvelocity

    etoile(:)%M = xyzmh_ptmass(4,:) * usolarmass

    etoile(:)%Mdot = xyzmh_ptmass(16,:) * usolarmass / utime_scaled * year_to_s ! Accretion rate is in Msun/year

    etoile(:)%find_spectrum = .true.
 endif

 return

end subroutine phantom_2_mcfost

!*************************************************************************
! routine to to compute temperature from
! internal energy assuming a mix of gas and radiation
! pressure, where Trad = Tgas. That is, we solve the
! quartic equation
!
!  a*T^4 + 3/2*rho*kb*T/mu = rho*u
!
! to determine the temperature from the supplied density
! and internal energy (rho, u).
! INPUT:
!    rho - density [kg/m3]
!    u - internal energy [J/kg]
! OUTPUT:
!    temp - temperature [K]
!*************************************************************************
real(dp) function get_temp_from_u(rho,u,mu) result(temp)
 use constantes, only:kb_on_mh,radconst
 real(dp), intent(in) :: rho,u,mu
 real(dp) :: ft,dft,dt
 real(dp), parameter :: tol = 1.e-8
 integer :: its

 temp = (u*rho/radconst)**0.25
 return
 ! Take minimum of gas and radiation temperatures as initial guess
 temp = min(u*mu/(1.5*kb_on_mh),(u*rho/radconst)**0.25)
  print*,'rho,u  = ',rho,u,' T =',u*mu/(1.5*kb_on_mh),(u*rho/radconst)**0.25,radconst,kb_on_mh

 dt = huge(0.)
 its = 0
 do while (abs(dt) > tol*temp .and. its < 500)
    its = its + 1
    ft = u*rho - 1.5*kb_on_mh*temp*rho/mu - radconst*temp**4
    dft = - 1.5*kb_on_mh*rho/mu - 4.*radconst*temp**3
    dt = ft/dft ! Newton-Raphson
    if (temp - dt > 1.2*temp) then
       temp = 1.2*temp
    elseif (temp - dt < 0.8*temp) then
       temp = 0.8*temp
    else
       temp = temp - dt
    endif
 enddo
 print*,'converged to T=',temp

end function get_temp_from_u

end module read_phantom
