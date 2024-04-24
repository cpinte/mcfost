module read_phantom

  use parametres
  use io_phantom_utils
  use messages
  use constantes
  use mess_up_SPH

  implicit none

contains

  subroutine read_phantom_bin_files(iunit,n_files, filenames, x,y,z,h,vx,vy,vz,T_gas,particle_id,massgas,massdust,&
       rhogas,rhodust,extra_heating,ndusttypes,SPH_grainsizes,mask,n_SPH,ldust_moments,dust_moments,mass_per_h,ierr)

    integer,               intent(in) :: iunit, n_files
    character(len=*),dimension(n_files), intent(in) :: filenames
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,vx,vy,vz
    real(dp), intent(out), dimension(:),   allocatable :: rhogas,massgas,SPH_grainsizes,T_gas
    integer,  intent(out), dimension(:),   allocatable :: particle_id
    real(dp), intent(out), dimension(:,:), allocatable :: rhodust,massdust, dust_moments
    integer, dimension(:), allocatable, intent(out) :: mask
    real, intent(out), dimension(:), allocatable :: extra_heating
    integer, intent(out) :: ndusttypes,n_SPH,ierr
    logical, intent(out) :: ldust_moments
    real(dp), intent(out) :: mass_per_H
    integer, parameter :: maxarraylengths = 12
    integer, parameter :: nsinkproperties = 17
    integer(kind=8) :: number8(maxarraylengths)
    integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
    integer :: nblocks,narraylengths,nblockarrays,number,idust,ieos
    character(len=lentag) :: tag
    character(len=lenid)  :: fileid
    integer :: np,ntypes,nptmass,ngrains,dustfluidtype,ndudt,nptmass0,nptmass_j,nptmass_found
    integer, parameter :: maxtypes = 6
    integer, parameter :: maxfiles = 3
    integer, parameter :: maxinblock = 128
    integer, parameter :: n_nucleation = 6
    integer, allocatable, dimension(:) :: npartoftype
    real(dp), allocatable, dimension(:,:) :: massoftype !(maxfiles,maxtypes)
    real(dp) :: hfact,umass,utime,ulength,gmw,x2
    integer(kind=1), allocatable, dimension(:) :: itype,ifiles
    real(4),  allocatable, dimension(:) :: tmp
    real(dp), allocatable, dimension(:) :: grainsize,graindens
    real(dp), allocatable, dimension(:) :: dudt,tmp_dp,gastemperature
    real(dp), allocatable, dimension(:,:) :: xyzh,xyzmh_ptmass,vxyz_ptmass,dustfrac,vxyzu,nucleation
    type(dump_h) :: hdr
    logical :: got_h,got_dustfrac,got_itype,tagged,matched
    logical :: got_temperature,got_u,lpotential,do_nucleation
    integer :: ifile, np0, ntypes0, np_tot, ntypes_tot, ntypes_max, ndustsmall, ndustlarge

    ! We first read the number of particules in each phantom file
    np_tot = 0
    ntypes_tot = 0
    ntypes_max = 0
    do ifile=1, n_files
       write(*,*) "---- Reading header file #", ifile

       ! open file for read
       call open_dumpfile_r(iunit,filenames(ifile),fileid,ierr,requiretags=.true.)
       if (ierr /= 0) call error("cannot open "//trim(filenames(ifile)))

       if (fileid(2:2)=='T') then
          tagged = .true.
       else
          tagged = .false.
       endif
       call read_header(iunit,hdr,tagged,ierr)
       if (.not.tagged) call error("Phantom dump too old to be read by MCFOST")

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

       call extract('mass_per_H',mass_per_h,hdr,ierr)
       if (ierr /= 0) then
          ! No dust nucleation
          ldust_moments = .false.
       else
          ldust_moments = .true.
       endif

       call free_header(hdr, ierr)

       close(iunit)
       write(*,*) "---- Done"
    enddo ! ifile

    allocate(xyzh(4,np_tot),itype(np_tot),vxyzu(4,np_tot),gastemperature(np_tot))
    allocate(dustfrac(ndusttypes,np_tot),grainsize(ndusttypes),graindens(ndusttypes))
    allocate(dudt(np_tot),ifiles(np_tot),massoftype(n_files,ntypes_max),npartoftype(ntypes_tot))
    if (ldust_moments) allocate(nucleation(n_nucleation,np_tot))

    np0 = 0
    ntypes0 = 0
    do ifile=1, n_files
       write(*,*) "---- Reading data file #", ifile

       ! open file for read
       call open_dumpfile_r(iunit,filenames(ifile),fileid,ierr,requiretags=.true.)
       if (ierr /= 0) call error("cannot open "//trim(filenames(ifile)))

       if (fileid(2:2)=='T') then
          tagged = .true.
       else
          tagged = .false.
       endif
       call read_header(iunit,hdr,tagged,ierr)
       if (.not.tagged) call error("Phantom dump too old to be read by MCFOST")

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

       write(*,*) ' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
       write(*,*) ' nptmass = ', nptmass

       write(*,*) "Found", nptmass, "point masses"

       ! testing for binary potential
       call extract('x2',x2,hdr,ierr)
       if (ierr == 0) then
          nptmass = nptmass + 2
          write(*,*) "Found 2 gravitational potentials"
          lpotential = .true.
       else
          lpotential = .false.
          ierr = 0
       endif
       write(*,*) "Converting them to", nptmass, "stars/planets"

       if (nptmass > 0) then
          allocate(xyzmh_ptmass(nsinkproperties,nptmass), vxyz_ptmass(3,nptmass)) !HACK
          xyzmh_ptmass(:,:) = 0.
          vxyz_ptmass(:,:) = 0.
       endif

       allocate(tmp(np), tmp_dp(np))

       if (npartoftype(2) > 0) then
          dustfluidtype = 2
          call warning("Phantom dump contains two-fluid dust particles, may be discarded")
       else
          dustfluidtype = 1
       endif

       ! extract info from real header
       call extract('massoftype',massoftype(ifile,1:ntypes),hdr,ierr)
       call extract('hfact',hfact,hdr,ierr)
       call extract('gmw',gmw,hdr,ierr)
       call extract('ieos',ieos,hdr,ierr)
       if (ierr /= 0) then
          write(*,*) "Using mcfost value instead : ", mu
          gmw=mu
       endif
       if (ndusttypes > 0) then
          call extract('grainsize',grainsize(1:ndusttypes),hdr,ierr) ! code units here
          call extract('graindens',graindens(1:ndusttypes),hdr,ierr)
       endif

       call extract('umass',umass,hdr,ierr)
       call extract('utime',utime,hdr,ierr)
       call extract('udist',ulength,hdr,ierr)
       call extract('time',simu_time,hdr,ierr)

       read (iunit, iostat=ierr) number
       if (ierr /= 0) return
       narraylengths = number/nblocks

       got_h = .false.
       got_dustfrac = .false.
       got_itype = .false.
       got_temperature = .false.
       got_u = .false.
       ! skip each block that is too small
       nblockarrays = narraylengths*nblocks

       ndudt = 0
       ngrains = 0
       nptmass_found = 0
       do iblock = 1,nblocks
          call read_block_header(narraylengths,number8,nums,iunit,ierr)
          if (ierr /= 0) call error('Reading block header')
          do j=1,narraylengths
             do i=1,ndatatypes
                do k=1,nums(i,j)
                   if (j==1 .and. number8(j)==np) then
                      read(iunit, iostat=ierr) tag
                      if (ierr /= 0) call error('Reading tag')
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
                            got_u = .true.
                         case('temperature')
                            read(iunit,iostat=ierr) tmp_dp ; gastemperature(np0+1:np0+np) = tmp_dp
                            got_temperature = .true.
                         case('dustfrac')
                            ngrains = ngrains + 1
                            if (ngrains > ndusttypes) then
                               write(*,*) "ERROR ngrains > ndusttypes:", ngrains, ndusttypes
                               ierr = 1 ;
                               return
                            endif
                            read(iunit,iostat=ierr) dustfrac(ngrains,:)
                            got_dustfrac = .true.
                         ! Nucleation
                         case('K0')
                            read(iunit,iostat=ierr) tmp_dp ; nucleation(1,np0+1:np0+np) = tmp_dp
                         case('K1')
                            read(iunit,iostat=ierr) tmp_dp ; nucleation(2,np0+1:np0+np) = tmp_dp
                         case('K2')
                            read(iunit,iostat=ierr) tmp_dp ; nucleation(3,np0+1:np0+np) = tmp_dp
                         case('K3')
                            read(iunit,iostat=ierr) tmp_dp ; nucleation(4,np0+1:np0+np) = tmp_dp
                         case('mu')
                            read(iunit,iostat=ierr) tmp_dp ; nucleation(5,np0+1:np0+np) = tmp_dp
                         case('kappa')
                            read(iunit,iostat=ierr) tmp_dp ; nucleation(6,np0+1:np0+np) = tmp_dp
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
                         case default
                            read(iunit,iostat=ierr)
                         end select
                      else
                         read(iunit,iostat=ierr)
                      endif
                      if (ierr /= 0) call error("Error reading tag: "//trim(tag))
                   elseif (j==2 .and. number8(j) > 0) then ! sink particles
                      read(iunit,iostat=ierr) tag
                      matched = .true.
                      if (i==i_real .or. i==i_real8) then
                         select case(trim(tag))
                         case('x')
                            nptmass_j = number8(j) ! We rely on x being the 1st data array for sink particles
                            nptmass0 = nptmass_found+1
                            nptmass_found = nptmass_found + nptmass_j
                            read(iunit,iostat=ierr) xyzmh_ptmass(1,nptmass0:nptmass_found)
                         case('y')
                            read(iunit,iostat=ierr) xyzmh_ptmass(2,nptmass0:nptmass_found)
                         case('z')
                            read(iunit,iostat=ierr) xyzmh_ptmass(3,nptmass0:nptmass_found)
                         case('m')
                            read(iunit,iostat=ierr) xyzmh_ptmass(4,nptmass0:nptmass_found)
                         case('h')
                            read(iunit,iostat=ierr) xyzmh_ptmass(5,nptmass0:nptmass_found)
                         case('hsoft')
                            read(iunit,iostat=ierr) xyzmh_ptmass(6,nptmass0:nptmass_found)
                         case('lum')
                            read(iunit,iostat=ierr) xyzmh_ptmass(12,nptmass0:nptmass_found)
                         case('Teff')
                            read(iunit,iostat=ierr) xyzmh_ptmass(13,nptmass0:nptmass_found)
                         case('Reff')
                            read(iunit,iostat=ierr) xyzmh_ptmass(14,nptmass0:nptmass_found)
                         case('mdotav')
                            read(iunit,iostat=ierr) xyzmh_ptmass(16,nptmass0:nptmass_found)
                         case('vx')
                            read(iunit,iostat=ierr) vxyz_ptmass(1,nptmass0:nptmass_found)
                         case('vy')
                            read(iunit,iostat=ierr) vxyz_ptmass(2,nptmass0:nptmass_found)
                         case('vz')
                            read(iunit,iostat=ierr) vxyz_ptmass(3,nptmass0:nptmass_found)
                         case default
                            matched = .false.
                            read(iunit,iostat=ierr)
                         end select
                      else
                         matched = .false.
                         read(iunit,iostat=ierr)
                      endif
                   else
                      read(iunit, iostat=ierr) tag ! tag
                      read(iunit, iostat=ierr) ! array
                   endif
                enddo
             enddo
          enddo
       enddo ! block

       if (lpotential) then
          call extract('x2',xyzmh_ptmass(1,nptmass_found+1),hdr,ierr)
          call extract('y2',xyzmh_ptmass(2,nptmass_found+1),hdr,ierr)
          call extract('z2',xyzmh_ptmass(3,nptmass_found+1),hdr,ierr)
          call extract('m2',xyzmh_ptmass(4,nptmass_found+1),hdr,ierr)
          call extract('h2',xyzmh_ptmass(5,nptmass_found+1),hdr,ierr)
          call extract('vx2',vxyz_ptmass(1,nptmass_found+1),hdr,ierr)
          call extract('vy2',vxyz_ptmass(2,nptmass_found+1),hdr,ierr)
          call extract('vz2',vxyz_ptmass(3,nptmass_found+1),hdr,ierr)

          call extract('x1',xyzmh_ptmass(1,nptmass_found+2),hdr,ierr)
          call extract('y1',xyzmh_ptmass(2,nptmass_found+2),hdr,ierr)
          call extract('z1',xyzmh_ptmass(3,nptmass_found+2),hdr,ierr)
          call extract('m1',xyzmh_ptmass(4,nptmass_found+2),hdr,ierr)
          call extract('h1',xyzmh_ptmass(5,nptmass_found+2),hdr,ierr)
          call extract('vx1',vxyz_ptmass(1,nptmass_found+2),hdr,ierr)
          call extract('vy1',vxyz_ptmass(2,nptmass_found+2),hdr,ierr)
          call extract('vz1',vxyz_ptmass(3,nptmass_found+2),hdr,ierr)
       endif

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

    if (got_temperature) then
       write(*,*) 'Using temperature from phantom temperature writeout'
    elseif(got_u) then
       write(*,*) 'Calculating temperature from u'
       ! Phantom uses erg/g for u, and mcfost uses J/K for kb.
       gastemperature = vxyzu(4,:) * 2./3.*gmw*amu/kb*erg_to_J * (ulength/utime)**2
       write(*,*) 'Maximum temperature = ', maxval(gastemperature)
    else
       write(*,*) 'Gas temperature not found, setting to T=cmb to avoid dust sublimation'
       gastemperature = 2.74
    endif

    if (got_h) then
       call modify_dump(np, nptmass, xyzh, vxyzu, xyzmh_ptmass, ulength, mask)

       call phantom_2_mcfost(np_tot,nptmass,ntypes_max,ndusttypes,ldust_moments,n_files,dustfluidtype,xyzh,&
            vxyzu,gastemperature,itype,grainsize,dustfrac,nucleation,massoftype,xyzmh_ptmass,vxyz_ptmass,&
            hfact,umass,utime,ulength,graindens,ndudt,dudt,ifiles, &
            n_SPH,x,y,z,h,vx,vy,vz,T_gas,particle_id, &
            SPH_grainsizes,massgas,massdust,rhogas,rhodust,dust_moments,extra_heating,ieos)
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

  subroutine read_phantom_hdf_files(iunit,n_files, filenames, x,y,z,h,vx,vy,vz,T_gas,&
       particle_id,massgas,massdust,rhogas,rhodust, &
       extra_heating,ndusttypes,SPH_grainsizes,     &
       mask,n_SPH,ldust_moments,dust_moments,mass_per_H,ierr)

    use utils_hdf5, only:open_hdf5file,    &
         close_hdf5file,   &
         open_hdf5group,   &
         close_hdf5group,  &
         read_from_hdf5,   &
         HID_T

    integer, intent(in) :: iunit, n_files
    character(len=*),dimension(n_files), intent(in) :: filenames
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,        &
         vx,vy,vz,T_gas, &
         rhogas,massgas, &
         SPH_grainsizes
    integer,  intent(out), dimension(:),   allocatable :: particle_id
    real(dp), intent(out), dimension(:,:), allocatable :: rhodust,massdust,dust_moments
    integer, dimension(:), allocatable, intent(out) :: mask
    real, intent(out), dimension(:), allocatable :: extra_heating
    integer, intent(out) :: ndusttypes,n_SPH,ierr
    logical, intent(out) :: ldust_moments
    real(dp), intent(out) :: mass_per_H

    character(len=200) :: filename

    logical :: got_dustfrac,got_itype

    integer :: ifile, np0, ntypes0, np_tot, ntypes_tot, ntypes_max
    integer :: np,ntypes,nptmass,dustfluidtype,ndudt
    integer :: ndustsmall,ndustlarge,ieos

    integer, parameter :: maxtypes = 100
    integer, parameter :: nsinkproperties = 17

    integer(kind=1), allocatable, dimension(:) :: itype, ifiles
    real(4),  allocatable, dimension(:)   :: tmp, tmp_header
    real(dp), allocatable, dimension(:)   :: dudt, tmp_dp,gastemperature
    real(dp), allocatable, dimension(:)   :: grainsize, graindens
    real(dp), allocatable, dimension(:,:) :: xyzh,xyzmh_ptmass,vxyz_ptmass,dustfrac,vxyzu,nucleation
    integer, allocatable, dimension(:) :: npartoftype
    real(dp), allocatable, dimension(:,:) :: massoftype !(maxfiles,maxtypes)
    real(dp) :: hfact,umass,utime,ulength

    integer(HID_T) :: hdf5_file_id
    integer(HID_T) :: hdf5_group_id

    logical :: got

    ldust_moments=.false.

    ierr = 0
    np_tot = 0
    ntypes_tot = 0
    ntypes_max = 0

    ! Read file headers to get particle numbers
    file_headers: do ifile=1, n_files

       write(*,*) "---- Reading header file #", ifile

       filename = trim(filenames(ifile))

       ! open file
       call open_hdf5file(filename,hdf5_file_id,ierr)
       if (ierr /= 0) call error("cannot open Phantom HDF file: "//trim(filename))

       ! open header group
       call open_hdf5group(hdf5_file_id,'header',hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot open Phantom HDF header group in "//trim(filename))

       ! read header values
       call read_from_hdf5(np,'nparttot',hdf5_group_id,got,ierr)
       call read_from_hdf5(ntypes,'ntypes',hdf5_group_id,got,ierr)
       call read_from_hdf5(ndusttypes,'ndusttypes',hdf5_group_id,got,ierr)
       if (.not. got) then
          ! ndusttypes is for pre-largegrain multigrain headers
          call read_from_hdf5(ndustsmall,'ndustsmall',hdf5_group_id,got,ierr)
          call read_from_hdf5(ndustlarge,'ndustlarge',hdf5_group_id,got,ierr)
          ! ndusttype must be the same for all files : todo : add a test
          ndusttypes = ndustsmall + ndustlarge
       endif

       ! close the header group
       call close_hdf5group(hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot close Phantom HDF header group in "//trim(filename))

       ! close file
       call close_hdf5file(hdf5_file_id,ierr)
       if (ierr /= 0) call error("cannot close Phantom HDF file: "//trim(filename))

       np_tot = np_tot + np
       ntypes_tot = ntypes_tot + ntypes
       ntypes_max = max(ntypes_max,ntypes)

    enddo file_headers

    ! Allocate arrays for SPH data
    allocate(xyzh(4,np_tot),                        &
         itype(np_tot),                         &
         vxyzu(4,np_tot),                       &
         gastemperature(np_tot),                &
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
       call open_hdf5file(filename,hdf5_file_id,ierr)
       if (ierr /= 0) call error("cannot open Phantom HDF file: "//trim(filename))

       ! open header group
       call open_hdf5group(hdf5_file_id,'header',hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot open Phantom HDF header group in "//trim(filename))

       ! read from header
       call read_from_hdf5(np,'nparttot',hdf5_group_id,got,ierr)
       call read_from_hdf5(ntypes,'ntypes',hdf5_group_id,got,ierr)
       call read_from_hdf5(npartoftype,'npartoftype',hdf5_group_id,got,ierr)
       call read_from_hdf5(nptmass,'nptmass',hdf5_group_id,got,ierr)
       call read_from_hdf5(ndusttypes,'ndusttypes',hdf5_group_id,got,ierr)
       if (.not. got) then
          ! ndusttypes is for pre-largegrain multigrain headers
          call read_from_hdf5(ndustsmall,'ndustsmall',hdf5_group_id,got,ierr)
          call read_from_hdf5(ndustlarge,'ndustlarge',hdf5_group_id,got,ierr)
          ! ndusttype must be the same for all files : todo : add a test
          ndusttypes = ndustsmall + ndustlarge
       endif
       call read_from_hdf5(massoftype(ifile,1:ntypes),'massoftype',hdf5_group_id,got,ierr)
       call read_from_hdf5(hfact,'hfact',hdf5_group_id,got,ierr)
       if (ndusttypes > 0) then
          allocate(tmp_header(np))
          call read_from_hdf5(tmp_header,'grainsize',hdf5_group_id,got,ierr)
          grainsize(1:ndusttypes) = tmp_header(1:ndusttypes)
          call read_from_hdf5(tmp_header,'graindens',hdf5_group_id,got,ierr)
          graindens(1:ndusttypes) = tmp_header(1:ndusttypes)
          deallocate(tmp_header)
       endif
       call read_from_hdf5(umass,'umass',hdf5_group_id,got,ierr)
       call read_from_hdf5(utime,'utime',hdf5_group_id,got,ierr)
       call read_from_hdf5(ulength,'udist',hdf5_group_id,got,ierr)

       if (ierr /= 0) call error("cannot read Phantom HDF header group in "//trim(filename))

       ! close the header group
       call close_hdf5group(hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot close Phantom HDF header group in "//trim(filename))

       write(*,*) ' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
       !write(*,*) ' npartoftype = ',npartoftype(ntypes0+1:ntypes0+ntypes)
       write(*,*) ' nptmass = ', nptmass

       allocate(tmp(np), tmp_dp(np))

       if (npartoftype(2) > 0) then
          dustfluidtype = 2
          call warning("Phantom dump contains two-fluid dust particles, may be discarded")
       else
          dustfluidtype = 1
       endif

       ! open particles group
       call open_hdf5group(hdf5_file_id,'particles',hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot open Phantom HDF particles group in "//trim(filename))

       got_dustfrac = .false.
       got_itype = .false.
       ndudt = 0

       ! TODO: read particle arrays
       call read_from_hdf5(itype(np0+1:np0+np),'itype',hdf5_group_id,got,ierr)
       if (got) got_itype = .true.
       call read_from_hdf5(xyzh(1:3,np0+1:np0+np),'xyz',hdf5_group_id,got,ierr)
       call read_from_hdf5(xyzh(4,np0+1:np0+np),'h',hdf5_group_id,got,ierr)
       call read_from_hdf5(vxyzu(1:3,np0+1:np0+np),'vxyz',hdf5_group_id,got,ierr)
       call read_from_hdf5(dustfrac(1:ndusttypes,np0+1:np0+np),'dustfrac',hdf5_group_id,got,ierr)
       if (got) got_dustfrac = .true.

       if (ierr /= 0) call error("cannot read Phantom HDF particles group in "//trim(filename))

       ! close the particles group
       call close_hdf5group(hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot close Phantom HDF particles group in "//trim(filename))

       ! open sinks group
       call open_hdf5group(hdf5_file_id,'sinks',hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot open Phantom HDF sinks group in "//trim(filename))

       ! read sinks
       allocate(xyzmh_ptmass(nsinkproperties,nptmass), vxyz_ptmass(3,nptmass))
       call read_from_hdf5(xyzmh_ptmass(1:3,:),'xyz',hdf5_group_id,got,ierr)
       call read_from_hdf5(xyzmh_ptmass(4,:),'m',hdf5_group_id,got,ierr)
       call read_from_hdf5(xyzmh_ptmass(5,:),'h',hdf5_group_id,got,ierr)
       call read_from_hdf5(xyzmh_ptmass(6,:),'hsoft',hdf5_group_id,got,ierr)
       call read_from_hdf5(vxyz_ptmass(:,:),'vxyz',hdf5_group_id,got,ierr)

       if (ierr /= 0)  call error("cannot read Phantom HDF sinks group in "//trim(filename))

       call read_from_hdf5(xyzmh_ptmass(16,:),'mdotav',hdf5_group_id,got,ierr)
       if (ierr /= 0) then
          call warning("mdotav not present")
          ierr =0
       endif

       ! close the sinks group
       call close_hdf5group(hdf5_group_id,ierr)
       if (ierr /= 0) call error("cannot close Phantom HDF sinks group in "//trim(filename))

       ! close file
       call close_hdf5file(hdf5_file_id,ierr)
       if (ierr /= 0) call error("cannot close Phantom HDF file "//trim(filename))

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

    gastemperature = 2.74
    write(*,*) 'Gas temperature not read from hdf files, setting to T=cmb to avoid dust sublimation'

    write(*,*) "Found", nptmass, "point masses in the phantom file"

    call modify_dump(np, nptmass, xyzh, vxyzu, xyzmh_ptmass, ulength, mask)

    call phantom_2_mcfost(np_tot,nptmass,ntypes_max,ndusttypes,ldust_moments,n_files,      &
         dustfluidtype,xyzh,vxyzu,gastemperature,itype,grainsize,dustfrac,nucleation, &
         massoftype,xyzmh_ptmass,vxyz_ptmass,hfact,umass,       &
         utime,ulength,graindens,ndudt,dudt,ifiles,   &
         n_SPH,x,y,z,h,vx,vy,vz,T_gas,particle_id,SPH_grainsizes,     &
         massgas,massdust,rhogas,rhodust,dust_moments,extra_heating,ieos)

    write(*,"(a,i8,a)") ' Using ',n_SPH,' particles from Phantom file'

    write(*,*) "Phantom dump file processed ok"
    ierr = 0
    deallocate(xyzh,itype,vxyzu)
    if (allocated(xyzmh_ptmass)) deallocate(xyzmh_ptmass,vxyz_ptmass)

  end subroutine read_phantom_hdf_files

  !*************************************************************************


  subroutine modify_dump(np, nptmass, xyzh, vxyzu, xyzmh_ptmass, ulength, mask)

    integer, intent(in) :: np, nptmass
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: xyzh, vxyzu
    real(kind=dp), allocatable, dimension(:,:), intent(inout) :: xyzmh_ptmass
    real(kind=dp), intent(in) :: ulength
    real(kind=dp), dimension(3) :: centre
    integer, dimension(:), allocatable, intent(out) :: mask

    integer :: i

    ! Modifying SPH dump
    if (ldelete_Hill_sphere .or. ldelete_inside_rsph .or. ldelete_outside_rsph .or. ldelete_above_theta) then
       allocate(mask(np))
       mask(:) = 0
    endif

    if (ldelete_Hill_sphere)  call mask_Hill_sphere(np, nptmass, xyzh, xyzmh_ptmass,ulength, mask)
    if (ldelete_inside_rsph)  call mask_inside_rsph(np, xyzh, ulength, rsph_min, mask)
    if (ldelete_outside_rsph) call mask_outside_rsph(np, xyzh, ulength, rsph_max, mask)
    if (ldelete_above_theta)  call mask_above_theta(np, xyzh, ulength, theta_max, mask)

    if (lrandomize_azimuth)     call randomize_azimuth(np, xyzh, vxyzu, mask)
    if (lrandomize_gap)         call randomize_gap(np, nptmass, xyzh, vxyzu, xyzmh_ptmass,ulength, gap_factor, .true.)
    if (lrandomize_outside_gap) call randomize_gap(np, nptmass, xyzh, vxyzu, xyzmh_ptmass,ulength, gap_factor, .false.)

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

  subroutine phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,ldust_moments,n_files,dustfluidtype,xyzh, &
       vxyzu,gastemperature,iphase,grainsize,dustfrac,nucleation,massoftype,xyzmh_ptmass,vxyz_ptmass,hfact,umass, &
       utime, ulength,graindens,ndudt,dudt,ifiles, &
       n_SPH,x,y,z,h,vx,vy,vz,T_gas,particle_id, &
       SPH_grainsizes, massgas,massdust, rhogas,rhodust,dust_moments,extra_heating,ieos)

    ! Convert phantom quantities & units to mcfost quantities & units
    !
    ! Outputs :
    ! x,y,z are in au
    ! rhodust & rhogas are in g/cm3
    ! extra_heating is in W

    use constantes, only : au_to_cm,Msun_to_g,erg_to_J,m_to_cm, Lsun, cm_to_mum, deg_to_rad, Ggrav
    use parametres, only : ldudt_implicit,ufac_implicit, lplanet_az, planet_az, lfix_star, RT_az_min, RT_az_max, RT_n_az
    use parametres, only : lscale_length_units,scale_length_units_factor,lscale_mass_units,scale_mass_units_factor

    ! Input arguments
    integer, intent(in) :: np,nptmass,ntypes,ndusttypes, n_files,dustfluidtype
    logical, intent(in) :: ldust_moments
    real(dp), dimension(4,np), intent(in) :: xyzh,vxyzu
    integer(kind=1), dimension(np), intent(in) :: iphase, ifiles
    real(dp), dimension(ndusttypes,np), intent(in) :: dustfrac
    real(dp), dimension(ndusttypes),    intent(in) :: grainsize ! code units
    real(dp), dimension(ndusttypes),    intent(in) :: graindens
    real(dp), dimension(n_files,ntypes), intent(in) :: massoftype
    real(dp), intent(in) :: hfact,umass,utime,ulength
    real(dp), dimension(:,:), intent(in) :: xyzmh_ptmass, vxyz_ptmass
    real(dp), dimension(:,:), intent(inout), allocatable :: nucleation
    real(dp), dimension(:), intent(in) :: gastemperature
    integer, intent(in) :: ndudt
    real(dp), dimension(:), intent(in) :: dudt

    ! Output arguments
    integer, intent(out) :: n_SPH, ieos
    real(dp), dimension(:),   allocatable, intent(out) :: x,y,z,h,vx,vy,vz,T_gas,rhogas,massgas ! massgas [Msun]
    integer, dimension(:),    allocatable, intent(out) :: particle_id
    real(dp), dimension(:,:), allocatable, intent(out) :: rhodust,massdust,dust_moments
    real(dp), dimension(:), allocatable, intent(out) :: SPH_grainsizes ! mum
    real, dimension(:), allocatable, intent(out) :: extra_heating


    type(star_type), dimension(:), allocatable :: etoile_old
    integer  :: i,j,k,itypei,alloc_status,i_etoile, n_etoiles_old, ifile
    real(dp) :: xi,yi,zi,hi,vxi,vyi,vzi,T_gasi,rhogasi,rhodusti,gasfraci,dustfraci,totlum,qtermi
    real(dp) :: ulength_scaled, umass_scaled, utime_scaled,udens,uerg_per_s,uWatt,ulength_au,ulength_m,usolarmass,uvelocity
    real(dp) :: vphi, vr, phi, cos_phi, sin_phi, r_cyl, r_cyl2, r_sph, G_phantom

    logical :: use_dust_particles = .false. ! 2-fluid: choose to use dust
    logical :: lupdate_photosphere

    ! We check the units by recomputing G
    G_phantom = ulength**3 / (utime**2 * umass)
    if (abs(G_phantom - Ggrav*1e3) > 1e-2 * G_phantom) call error("Phantom units are not consistent")

    ulength_scaled = ulength
    umass_scaled = umass
    if (lscale_length_units) then
       write(*,*) 'Lengths are rescaled by ', real(scale_length_units_factor)
       ulength_scaled = ulength * scale_length_units_factor
    endif
    if (lscale_mass_units) then
       write(*,*) 'Mass are rescaled by ', real(scale_mass_units_factor)
       umass_scaled = umass * scale_mass_units_factor
    endif
    utime_scaled = sqrt(ulength_scaled**3/(G_phantom*umass_scaled))

    udens = umass_scaled/ulength_scaled**3
    uerg_per_s = umass_scaled*ulength_scaled**2/utime_scaled**3
    uWatt = uerg_per_s * erg_to_J
    ulength_m = ulength_scaled / (m_to_cm)
    ulength_au = ulength_scaled / (au_to_cm)
    uvelocity =  ulength_scaled / (m_to_cm) / utime_scaled ! m/s
    usolarmass = umass_scaled / Msun_to_g

    simu_time = simu_time * utime_scaled

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
         x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),massgas(n_SPH),rhogas(n_SPH),T_gas(n_SPH),stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error in phanton_2_mcfost"
       write(*,*) "Exiting"
    endif

    ! Dust grain sizes in microns
    SPH_grainsizes(:) = grainsize(:) * ulength_scaled * cm_to_mum
    ! graindens * udens is in g/cm3

    if (lemission_mol .or. lemission_atom) then
       allocate(vx(n_SPH),vy(n_SPH),vz(n_SPH),stat=alloc_status)
       if (alloc_status /=0) then
          write(*,*) "Allocation error velocities in phanton_2_mcfost"
          write(*,*) "Exiting"
       endif
    endif

    if (ldust_moments) then
       allocate(dust_moments(4,n_SPH), stat=alloc_status)
       if (alloc_status /=0) then
          write(*,*) "Allocation error dust_moments in phanton_2_mcfost"
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
       T_gasi = gastemperature(i)

       itypei = abs(iphase(i))
       if (hi > 0.) then
          if (use_dust_particles .and. dustfluidtype==2 .and. ndusttypes==1 .and. itypei==2) then
             j = j + 1
             particle_id(j) = i
             x(j) = xi * ulength_au
             y(j) = yi * ulength_au
             z(j) = zi * ulength_au
             h(j) = hi * ulength_au
             if (lemission_mol .or. lemission_atom) then
                vx(j) = vxi * uvelocity
                vy(j) = vyi * uvelocity
                vz(j) = vzi * uvelocity
             endif
             T_gas(j) = T_gasi
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
             if (lemission_mol .or. lemission_atom) then
                vx(j) = vxi * uvelocity
                vy(j) = vyi * uvelocity
                vz(j) = vzi * uvelocity
             endif

             if (ldust_moments) dust_moments(:,j) = nucleation(1:4,i) ! indexing is different from phantom as I read starting at k0

             T_gas(j) = T_gasi
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

    lextra_heating = .false. ; ldudt_implicit = .false.
    if (.not.lno_internal_energy) then
       if (ndudt == np) then
          write(*,*) "Computing energy input"
          lextra_heating = .true.
          allocate(extra_heating(n_SPH), stat=alloc_status)
          if (alloc_status /=0) call error("Allocation error in phantom_2_mcfost")

          totlum = 0.
          do i=1,n_SPH
             qtermi = dudt(particle_id(i)) * uWatt
             totlum = totlum + qtermi
             extra_heating(i) = qtermi
          enddo

          write(*,*) "Total energy input = ",real(totlum),' W'
          write(*,*) "Total energy input = ",real(totlum/Lsun),' Lsun'
       endif

       ufac_implicit = 0.
       ldudt_implicit = .false.
    endif

    n_etoiles_old = n_etoiles
    if (lignore_sink) then
       n_etoiles = 0
    else
       write(*,*) "# Stars/planets:"
       n_etoiles = 0
       do i=1,nptmass
          n_etoiles = n_etoiles + 1
          if (real(xyzmh_ptmass(4,i)) * scale_mass_units_factor > 0.013) then
             write(*,*) "Star   #", i, "xyz=", real(xyzmh_ptmass(1:3,i) * ulength_au), "au, M=", &
                  real(xyzmh_ptmass(4,i) * usolarmass), "Msun, Mdot=", &
                  real(xyzmh_ptmass(16,i) * usolarmass / utime_scaled * year_to_s ), "Msun/yr"
          else
             write(*,*) "Planet #", i, "xyz=", real(xyzmh_ptmass(1:3,i) * ulength_au), "au, M=", &
                  real(xyzmh_ptmass(4,i) * GxMsun/GxMjup * usolarmass), "Mjup, Mdot=", &
                  real(xyzmh_ptmass(16,i) * usolarmass / utime_scaled * year_to_s ), "Msun/yr"
          endif
          if (i>1) write(*,*)  "       distance=", real(norm2(xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,1))*ulength_au), "au"
       enddo
    endif

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
             vphi = sqrt(Ggrav *  (xyzmh_ptmass(4,1) + xyzmh_ptmass(4,2)) * scale_mass_units_factor * Msun_to_kg  &
                  * (r_cyl * AU_to_m)**2 /  (r_sph * AU_to_m)**3 )

             vx(i) = vr * cos_phi - vphi * sin_phi
             vy(i) = vr * sin_phi + vphi * cos_phi
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
             T_gas(j) = get_temp_from_u(rhogasi,vxyzu(4,i)*uvelocity**2,0.6_dp,ieos)
             if (mod(j,100000).eq.0) print*,i,T_gas(j),'vel=',vxi*uvelocity,vyi*uvelocity,vzi*uvelocity,'m/s'
          endif
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
       write(*,*) "Stellar parameters will not be updated, only the star positions, velocities and masses"
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
       enddo
    else

       write(*,*) ""
       write(*,*) "Updating the stellar properties:"
       write(*,*) "There are now", n_etoiles, "stars in the model"

       ! Saving if the accretion rate was forced
       allocate(etoile_old(n_etoiles_old))
       if (allocated(etoile)) then
          etoile_old(:) = etoile(:)
          deallocate(etoile)
       endif
       allocate(etoile(n_etoiles))
       do i=1, min(n_etoiles, n_etoiles_old)
          etoile(i)%force_Mdot = etoile_old(i)%force_Mdot
          etoile(i)%Mdot = etoile_old(i)%Mdot
       enddo
       ! If we have new stars
       do i=n_etoiles_old+1,n_etoiles
          etoile(i)%force_Mdot = .false.
          etoile(i)%Mdot = 0.
       enddo
       deallocate(etoile_old)

       if (n_etoiles > 0) then
          etoile(:)%x = xyzmh_ptmass(1,:) * ulength_au
          etoile(:)%y = xyzmh_ptmass(2,:) * ulength_au
          etoile(:)%z = xyzmh_ptmass(3,:) * ulength_au

          etoile(:)%vx = vxyz_ptmass(1,:) * uvelocity
          etoile(:)%vy = vxyz_ptmass(2,:) * uvelocity
          etoile(:)%vz = vxyz_ptmass(3,:) * uvelocity

          etoile(:)%M = xyzmh_ptmass(4,:) * usolarmass

          do i=1, n_etoiles
             if (.not.etoile(i)%force_Mdot) then
                etoile(i)%Mdot = xyzmh_ptmass(16,i) * usolarmass / utime_scaled * year_to_s ! Accretion rate is in Msun/year
             endif
          enddo
          etoile(:)%find_spectrum = .true.

          lupdate_photosphere = .false.
          do i=1,n_etoiles
             if (xyzmh_ptmass(13,i) > tiny_real) then
                ! sink particle Teff is defined
                write(*,*) "Using Teff and Lum from phantom:"
                lupdate_photosphere = .true.
                etoile(i)%T = max(xyzmh_ptmass(13,i),100.)
                ! deriving radius from luminosity (Reff is not correct in dump)
                etoile(i)%r = sqrt(xyzmh_ptmass(12,i) * uWatt / Lsun) * (Tsun/etoile(i)%T)**2
             else if (xyzmh_ptmass(14,i) > tiny_real) then
                ! sink particle Reff is defined
                write(*,*) "Using Reff and Lum from phantom:"
                lupdate_photosphere = .true.
                etoile(i)%r = xyzmh_ptmass(14,i) * ulength_m / Rsun
                etoile(i)%T =  Tsun * (xyzmh_ptmass(12,i) * uWatt / Lsun)**0.25 / sqrt(etoile(i)%r)
             else if (xyzmh_ptmass(12,i) > tiny_real) then
                ! We use the luminosity and a radius estimate
                write(*,*) "Using Lum from phantom + effective radius from sink:"
                etoile(i)%T = 100.
                etoile(i)%r = 0.001
                if (xyzmh_ptmass(12,i) > tiny_real) then
                   etoile(i)%r = max(xyzmh_ptmass(5,i),xyzmh_ptmass(6,i)) * ulength_m / Rsun
                   if (etoile(i)%r > 0.) etoile(i)%T =  Tsun * (xyzmh_ptmass(12,i) * uWatt / Lsun)**0.25 / sqrt(etoile(i)%r)
                endif
             else
                if (lupdate_photosphere) then ! if we updated the 1st sink
                   etoile(i)%r = max(max(xyzmh_ptmass(5,i),xyzmh_ptmass(6,i)) * ulength_m / Rsun, 0.001)
                   etoile(i)%T = 100
                endif
             endif
          enddo

          if (lupdate_photosphere) then
             lfix_star = .true.
             do i=1,n_etoiles
                write(*,*) "Star #",i,"  Teff=", etoile(i)%T, "K, r=", real(etoile(i)%r), "Rsun"
             enddo
             ! Convert radius to au
             etoile(:)%r = etoile(:)%r * Rsun_to_AU

             etoile(:)%lb_body=.true.
          endif
       endif
    endif

    return

  end subroutine phantom_2_mcfost

  !*************************************************************************

  real(dp) function get_temp_from_u(rho,u,mu,ieos) result(temp)
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

    use constantes, only:kb_on_mh,radconst
    real(dp), intent(in) :: rho,u,mu
    integer, intent(in) :: ieos
    real(dp) :: ft,dft,dt
    real(dp), parameter :: tol = 1.e-8
    integer :: its

    if (ieos==2) then
       temp = u*mu/(1.5*kb_on_mh)
       ! temp = (u*rho/radconst)**0.25
    elseif (ieos==12) then
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
    endif

  end function get_temp_from_u

end module read_phantom
