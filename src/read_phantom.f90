module read_phantom

  use parametres
  use dump_utils
  use messages

  implicit none

  contains

    subroutine read_phantom_files(iunit,n_files, filenames, x,y,z,h,vx,vy,vz,particle_id,massgas,massdust,&
         rhogas,rhodust,extra_heating,ndusttypes,SPH_grainsizes,n_SPH,ierr)

 integer,               intent(in) :: iunit, n_files
 character(len=*),dimension(n_files), intent(in) :: filenames
 real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h, vx,vy,vz, rhogas,massgas,SPH_grainsizes
 integer,  intent(out), dimension(:),   allocatable :: particle_id
 real(dp), intent(out), dimension(:,:), allocatable :: rhodust,massdust

 real, intent(out), dimension(:), allocatable :: extra_heating
 integer, intent(out) :: ndusttypes,n_SPH,ierr

 integer, parameter :: maxarraylengths = 12
 integer(kind=8) :: number8(maxarraylengths)
 integer :: i,j,k,iblock,nums(ndatatypes,maxarraylengths)
 integer :: nblocks,narraylengths,nblockarrays,number
 character(len=lentag) :: tag
 character(len=lenid)  :: fileid
 integer :: np,ntypes,nptmass,ipos,ngrains,dustfluidtype,ndudt
 integer, parameter :: maxtypes = 6
 integer, parameter :: maxfiles = 3
 integer, allocatable, dimension(:) :: npartoftype
 real(dp), allocatable, dimension(:,:) :: massoftype !(maxfiles,maxtypes)
 real(dp) :: hfact,umass,utime,udist
 integer(kind=1), allocatable, dimension(:) :: itype, ifiles
 real(4),  allocatable, dimension(:) :: tmp
 real(dp), allocatable, dimension(:) :: grainsize, graindens
 real(dp), allocatable, dimension(:) :: dudt, tmp_dp
 real(dp), allocatable, dimension(:,:) :: xyzh,xyzmh_ptmass,dustfrac,vxyzu
 type(dump_h) :: hdr
 logical :: got_h,got_dustfrac,got_itype,tagged,matched

 integer :: ifile, np0, ntypes0, np_tot, ntypes_tot, np_max, ntypes_max, ndustsmall, ndustlarge

 ! We first read the number of particules in each phantom file
 np_tot = 0
 np_max = 0
 ntypes_tot = 0
 ntypes_max = 0
 do ifile=1, n_files
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
    endif
    call extract('ntypes',ntypes,hdr,ierr)

    np_tot = np_tot + np
    np_max = max(np_max,np)
    ntypes_tot = ntypes_tot + ntypes
    ntypes_max = max(ntypes_max,ntypes)
    close(iunit)
 enddo ! ifile

 allocate(xyzh(4,np_tot),itype(np_tot),tmp(np_max),vxyzu(4,np_tot),tmp_dp(np_max))
 allocate(dustfrac(ndusttypes,np_tot),grainsize(ndusttypes),graindens(ndusttypes))
 allocate(dudt(np_tot),ifiles(np_tot),massoftype(n_files,ntypes_max),npartoftype(ntypes_tot))


 np0 = 0
 ntypes0 = 0
 do ifile=1, n_files
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
       call extract('ndustlarge',ndustlarge,hdr,ierr,default=0)
       ! ndusttype must be the same for all files : todo : add a test
       ndusttypes = ndustsmall + ndustlarge
    endif
    call extract('nptmass',nptmass,hdr,ierr,default=0)
    !call extract('isink',isink,hdr,ierr,default=0)

    write(*,*) ' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
    write(*,*) ' npartoftype = ',npartoftype(ntypes0+1:ntypes0+ntypes)
    write(*,*) ' nptmass = ', nptmass

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
       do j=1,narraylengths
          !write(*,*) 'block ',iblock, j, number8(j)
          do i=1,ndatatypes
             !print*,' data type ',i,' arrays = ',nums(i,j)
             do k=1,nums(i,j)
                if (j==1 .and. number8(j)==np) then
                   read(iunit, iostat=ierr) tag
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
                !elseif (j==1 .and. number8(j)==nptmass) then
                elseif (j==2) then ! HACK : what is j exactly anyway ? and why would we need to test for j==1
                   nptmass = number8(j) ! HACK
                   if (.not.allocated(xyzmh_ptmass)) allocate(xyzmh_ptmass(5,nptmass)) !HACK

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
    call phantom_2_mcfost(np_tot,nptmass,ntypes_max,ndusttypes,n_files,dustfluidtype,xyzh,&
         vxyzu,itype,grainsize,dustfrac,massoftype,xyzmh_ptmass,&
         hfact,umass,utime,udist,graindens,ndudt,dudt,ifiles, &
         n_SPH,x,y,z,h,vx,vy,vz,particle_id, &
         SPH_grainsizes,massgas,massdust,rhogas,rhodust,extra_heating)
    write(*,"(a,i8,a)") ' Using ',n_SPH,' particles from Phantom file'
 else
    n_SPH = 0
    write(*,*) "ERROR reading h from file"
 endif

 write(*,*) "Phantom dump file processed ok"
 deallocate(xyzh,itype,tmp,vxyzu,tmp_dp)
 if (allocated(xyzmh_ptmass)) deallocate(xyzmh_ptmass)

end subroutine read_phantom_files

!*************************************************************************

subroutine phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,n_files,dustfluidtype,xyzh, &
     vxyzu,iphase,grainsize,dustfrac,massoftype,xyzmh_ptmass,hfact,umass, &
     utime, udist,graindens,ndudt,dudt,ifiles, &
     n_SPH,x,y,z,h,vx,vy,vz,particle_id, &
     SPH_grainsizes, massgas,massdust, rhogas,rhodust,extra_heating,T_to_u)

  ! Convert phantom quantities & units to mcfost quantities & units
  ! x,y,z are in au
  ! rhodust & rhogas are in g/cm3
  ! extra_heating is in W

  use constantes, only : au_to_cm,Msun_to_g,erg_to_J,m_to_cm, Lsun, cm_to_mum, deg_to_rad
  use parametres, only : ldudt_implicit,ufac_implicit, lplanet_az, planet_az, lfix_star, RT_az_min, RT_az_max, RT_n_az
  use parametres, only : lscale_units,scale_units_factor

  integer, intent(in) :: np,nptmass,ntypes,ndusttypes,dustfluidtype, n_files
  real(dp), dimension(4,np), intent(in) :: xyzh,vxyzu
  integer(kind=1), dimension(np), intent(in) :: iphase, ifiles
  real(dp), dimension(ndusttypes,np), intent(in) :: dustfrac
  real(dp), dimension(ndusttypes),    intent(in) :: grainsize ! code units
  real(dp), dimension(ndusttypes),    intent(in) :: graindens
  real(dp), dimension(n_files,ntypes), intent(in) :: massoftype
  real(dp), intent(in) :: hfact,umass,utime,udist
  real(dp), dimension(:,:), intent(in) :: xyzmh_ptmass
  integer, intent(in) :: ndudt
  real(dp), dimension(:), intent(in) :: dudt

  ! MC
  integer, intent(out) :: n_SPH
  real(dp), dimension(:),   allocatable, intent(out) :: x,y,z,h,vx,vy,vz,rhogas,massgas ! massgas [Msun]
  integer, dimension(:),    allocatable, intent(out) :: particle_id
  real(dp), dimension(:,:), allocatable, intent(out) :: rhodust,massdust
  real(dp), dimension(:), allocatable, intent(out) :: SPH_grainsizes ! mum
  real, dimension(:), allocatable, intent(out) :: extra_heating
  real(dp), intent(in), optional :: T_to_u


  integer  :: i,j,k,itypei,alloc_status,i_etoiles, ifile
  real(dp) :: xi,yi,zi,hi,vxi,vyi,vzi,rhogasi,rhodusti,gasfraci,dustfraci,totlum,qtermi
  real(dp) :: udist_scaled, umass_scaled, udens,uerg_per_s,uWatt,ulength_au,usolarmass,uvelocity

  logical :: use_dust_particles = .false. ! 2-fluid: choose to use dust


 if (lscale_units) then
    write(*,*) 'Lengths are rescaled by ', scale_units_factor
    udist_scaled = udist * scale_units_factor
    umass_scaled = umass * scale_units_factor**3
 else
    udist_scaled = udist
    umass_scaled = umass
 endif

  udens = umass_scaled/udist_scaled**3
  uerg_per_s = umass_scaled*udist_scaled**2/utime**3
  uWatt = uerg_per_s * erg_to_J
  ulength_au = udist_scaled/ (au_to_cm)
  uvelocity =  udist_scaled / (m_to_cm) / utime ! m/s
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
      x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),massgas(n_SPH),rhogas(n_SPH),stat=alloc_status)
 if (alloc_status /=0) then
    write(*,*) "Allocation error in phanton_2_mcfost"
    write(*,*) "Exiting"
 endif

 ! Dust grain sizes in microns
 SPH_grainsizes(:) = grainsize(:) * udist_scaled * cm_to_mum
 ! graindens * udens is in g/cm3

 if (lemission_mol) then
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
 n_etoiles = 0
 do i=1,nptmass
    write(*,*) "Sink #", i, "xyz=", real(xyzmh_ptmass(1:3,i)), "au, M=", real(xyzmh_ptmass(4,i)), "Msun"
    if (i>1) write(*,*)  "       distance=", real(norm2(xyzmh_ptmass(1:3,i) - xyzmh_ptmass(1:3,1))), "au"
    if (xyzmh_ptmass(4,i) > 0.0124098) then ! 13 Jupiter masses
       n_etoiles = n_etoiles + 1
    endif
 enddo

 if (lplanet_az) then
    if (nptmass /= 2) call error("option -planet_az only works with 2 sink particles")
    RT_n_az = 1
    RT_az_min = planet_az + atan2(xyzmh_ptmass(2,2) - xyzmh_ptmass(2,1),xyzmh_ptmass(1,2) - xyzmh_ptmass(1,1)) / deg_to_rad
    RT_az_max = RT_az_min
    write(*,*) "WARNING: updating the azimuth to:", RT_az_min
 endif

 if (lfix_star) then
    write(*,*) ""
    write(*,*) "Stellar parameters will not be updated"
 else
    write(*,*) ""
    write(*,*) "Updating the stellar properties:"
    write(*,*) "There are", n_etoiles, "stars in the model"
    if (allocated(etoile)) deallocate(etoile)
    allocate(etoile(n_etoiles))

    i_etoiles = 0
    do i=1,nptmass
       if (xyzmh_ptmass(4,i) > 0.0124098) then ! 13 Jupiter masses
          i_etoiles = i_etoiles + 1
          etoile(i_etoiles)%x = xyzmh_ptmass(1,i) * ulength_au
          etoile(i_etoiles)%y = xyzmh_ptmass(2,i) * ulength_au
          etoile(i_etoiles)%z = xyzmh_ptmass(3,i) * ulength_au
          etoile(i_etoiles)%M = xyzmh_ptmass(4,i) * usolarmass
       endif
    enddo
 endif

 return

end subroutine phantom_2_mcfost

!*************************************************************************

end module read_phantom
