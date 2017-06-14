module read_phantom

  use parametres
  use dump_utils

  implicit none

  contains

subroutine read_phantom_file(iunit,filename,x,y,z,particle_id,massgas,massdust,&
      rhogas,rhodust,extra_heating,ndusttypes,grainsize,n_SPH,ierr)

 integer,               intent(in) :: iunit
 character(len=*),      intent(in) :: filename
 real(dp), intent(out), dimension(:),   allocatable :: x,y,z,rhogas,massgas,grainsize
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
 integer, parameter :: maxtypes = 10
 integer :: npartoftype(maxtypes)
 real(dp) :: massoftype(maxtypes),hfact,umass,utime,udist
 integer(kind=1), allocatable, dimension(:) :: itype
 real(4),  allocatable, dimension(:) :: tmp
 real(dp) :: graindens
 real(dp), allocatable, dimension(:) :: dudt
 real(dp), allocatable, dimension(:,:) :: xyzh,xyzmh_ptmass,dustfrac,vxyzu
 type(dump_h) :: hdr
 logical :: got_h,got_dustfrac,got_itype,tagged,matched

 ! open file for read
 call open_dumpfile_r(iunit,filename,fileid,ierr,requiretags=.true.)
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
 !print*,hdr%inttags(:)
 !print*,hdr%realtags(:)

 ! get nblocks
 call extract('nblocks',nblocks,hdr,ierr,default=1)
 call extract('nparttot',np,hdr,ierr)
 call extract('ntypes',ntypes,hdr,ierr)
 call extract('npartoftype',npartoftype(1:ntypes),hdr,ierr)
 call extract('ndusttypes',ndusttypes,hdr,ierr,default=0)
 call extract('nptmass',nptmass,hdr,ierr,default=0)
 !call extract('isink',isink,hdr,ierr,default=0)

 write(*,*) ' npart = ',np,' ntypes = ',ntypes, ' ndusttypes = ',ndusttypes
 write(*,*) ' npartoftype = ',npartoftype(1:ntypes)
 write(*,*) ' nptmass = ', nptmass

 if (npartoftype(2) > 0) then
    dustfluidtype = 2
    write(*,"(/,a,/)") ' *** WARNING: Phantom dump contains two-fluid dust particles, may be discarded ***'
 else
    dustfluidtype = 1
 endif

 allocate(xyzh(4,np),itype(np),tmp(np),vxyzu(4,np))
 allocate(dustfrac(ndusttypes,np),grainsize(ndusttypes))
 allocate(dudt(np))

 !allocate(xyzmh_ptmass(5,nptmass)) ! HACK : Bug :  nptmass not defined yet, the keyword does not exist in the dump
 ! itype = 1

 ! extract info from real header
 call extract('massoftype',massoftype(1:ntypes),hdr,ierr)
 call extract('hfact',hfact,hdr,ierr)
 if (ndusttypes > 0) then
    call extract('grainsize',grainsize(1:ndusttypes),hdr,ierr)
    call extract('graindens',graindens,hdr,ierr)
 endif
 !write(*,*) ' hfact = ',hfact
 !write(*,*) ' massoftype = ',massoftype(1:ntypes)

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
                      read(iunit,iostat=ierr) xyzh(1,1:np)
                   case('y')
                      read(iunit,iostat=ierr) xyzh(2,1:np)
                   case('z')
                      read(iunit,iostat=ierr) xyzh(3,1:np)
                   case('h')
                      read(iunit,iostat=ierr) xyzh(4,1:np)
                      got_h = .true.
                   case('vx')
                      read(iunit,iostat=ierr) vxyzu(1,1:np)
                   case('vy')
                      read(iunit,iostat=ierr) vxyzu(2,1:np)
                   case('vz')
                      read(iunit,iostat=ierr) vxyzu(3,1:np)
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
                      xyzh(4,1:np) = real(tmp(1:np),kind=dp)
                      got_h = .true.
                   case('dustfrac')
                      ngrains = ngrains + 1
                      read(iunit,iostat=ierr) tmp(1:np)
                      dustfrac(ngrains,1:np) = real(tmp(1:np),kind=dp)
                      got_dustfrac = .true.
                   case('luminosity')
                      read(iunit,iostat=ierr) tmp(1:np)
                      dudt(1:np) = real(tmp(1:np),kind=dp)
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
                      read(iunit,iostat=ierr) itype(1:np)
                   case default
                      read(iunit,iostat=ierr)
                   end select
                else
                   read(iunit,iostat=ierr)
                endif
              !  if (matched) then
              !     write(*,"(a)") '->'//trim(tag)
              !  else
              !     write(*,"(a)")
              !  endif
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
 enddo

 close(iunit)

 if (.not. got_itype) then
    itype = 1
 endif

 if ((ndusttypes > 0) .and. .not. got_dustfrac) then
    dustfrac = 0.
 endif

 if (ndudt == 0) then
    dudt = 0.
 endif

 write(*,*) "Found", nptmass, "point masses in the phantom file"

 if (got_h) then
    call phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,dustfluidtype,xyzh,&
         vxyzu,itype,grainsize,dustfrac,massoftype(1:ntypes),xyzmh_ptmass,&
         hfact,umass,utime,udist,graindens,ndudt,dudt,n_SPH,x,y,z,particle_id,&
         massgas,massdust,rhogas,rhodust,extra_heating)
    write(*,"(a,i8,a)") ' Using ',n_SPH,' particles from Phantom file'
 else
    n_SPH = 0
    write(*,*) "ERROR reading h from file"
 endif

 write(*,*) "Phantom dump file processed ok"
 deallocate(xyzh,itype,tmp,vxyzu)
 if (allocated(xyzmh_ptmass)) deallocate(xyzmh_ptmass)

end subroutine read_phantom_file

!*************************************************************************

subroutine phantom_2_mcfost(np,nptmass,ntypes,ndusttypes,dustfluidtype,xyzh,&
     vxyzu,iphase,grainsize,dustfrac,massoftype,xyzmh_ptmass,hfact,umass,utime,&
     udist,graindens,ndudt,dudt,n_SPH,x,y,z,particle_id,massgas,massdust,&
     rhogas,rhodust,extra_heating,T_to_u)

  ! Convert phantom quantities & units to mcfost quantities & units
  ! x,y,z are in au
  ! rhodust & rhogas are in g/cm3
  ! extra_heating is in W

  use constantes, only : au_to_cm,Msun_to_g,erg_to_J
  use prop_star
  use parametres, only : ldudt_implicit,ufac_implicit

  integer, intent(in) :: np,nptmass,ntypes,ndusttypes,dustfluidtype
  real(dp), dimension(4,np), intent(in) :: xyzh,vxyzu
  integer(kind=1), dimension(np), intent(in) :: iphase
  real(dp), dimension(ndusttypes,np), intent(in) :: dustfrac
  real(dp), dimension(ndusttypes),    intent(in) :: grainsize
  real(dp), intent(in) :: graindens
  real(dp), dimension(ntypes), intent(in) :: massoftype
  real(dp), intent(in) :: hfact,umass,utime,udist
  real(dp), dimension(:,:), intent(in) :: xyzmh_ptmass
  integer, intent(in) :: ndudt
  real(dp), dimension(:), intent(in) :: dudt

  real(dp), dimension(:),   allocatable, intent(out) :: x,y,z,rhogas,massgas
  integer, dimension(:),    allocatable, intent(out) :: particle_id
  real(dp), dimension(:,:), allocatable, intent(out) :: rhodust,massdust
  real, dimension(:), allocatable, intent(out) :: extra_heating
  integer, intent(out) :: n_SPH
  real(dp), intent(in), optional :: T_to_u
  integer  :: i,j,k,itypei,alloc_status,i_etoiles
  real(dp) :: xi,yi,zi,hi,rhogasi,rhodusti,udens,uerg_per_s,uWatt,ulength_au,usolarmass
  real(dp) :: gasfraci,dustfraci,totlum,qtermi

  logical :: use_dust_particles = .false. ! 2-fluid: choose to use dust particles for Voronoi mesh

  real, parameter :: Lsun = 3.839e26 ! W

  udens = umass/udist**3
  uerg_per_s = umass*udist**2/utime**3
  uWatt = uerg_per_s * erg_to_J
  ulength_au = udist/ (au_to_cm )
  usolarmass = umass/Msun_to_g

 if (dustfluidtype == 1) then
    ! 1-fluid: always use gas particles for Voronoi mesh
    use_dust_particles = .false.
 endif

 write(*,*) ''
 if (use_dust_particles) then
    write(*,*) '*** Using dust particles for Voronoi mesh ***'
 else
    write(*,*) '*** Using gas particles for Voronoi mesh ***'
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
 allocate(rhodust(ndusttypes,n_SPH),massdust(ndusttypes,n_SPH),particle_id(n_SPH),&
      x(n_SPH),y(n_SPH),z(n_SPH),massgas(n_SPH),rhogas(n_SPH),stat=alloc_status)
 if (alloc_status /=0) then
    write(*,*) "Allocation error in phanton_2_mcfost"
    write(*,*) "Exiting"
 endif

 j = 0
 do i=1,np
    xi = xyzh(1,i)
    yi = xyzh(2,i)
    zi = xyzh(3,i)
    hi = xyzh(4,i)
    itypei = abs(iphase(i))
    if (hi > 0.) then
       if (use_dust_particles .and. &
           dustfluidtype==2 .and. ndusttypes==1 .and. itypei==2) then
          j = j + 1
          particle_id(j) = i
          x(j) = xi * ulength_au
          y(j) = yi * ulength_au
          z(j) = zi * ulength_au
          rhodusti = massoftype(itypei)*(hfact/hi)**3  * udens ! g/cm**3
          gasfraci = dustfrac(1,i)
          rhodust(1,j) = rhodusti
          massdust(1,j) = massoftype(itypei) * usolarmass ! Msun
          rhogas(j) = gasfraci*rhodusti
          massgas(j) = gasfraci*massdust(1,j)
       elseif (.not. use_dust_particles .and. itypei==1) then
          j = j + 1
          particle_id(j) = i
          x(j) = xi * ulength_au
          y(j) = yi * ulength_au
          z(j) = zi * ulength_au
          rhogasi = massoftype(itypei)*(hfact/hi)**3  * udens ! g/cm**3
          dustfraci = sum(dustfrac(:,i))
          if (dustfluidtype==2) then
             rhogas(j) = rhogasi
          else
             rhogas(j) = (1 - dustfraci)*rhogasi
          endif
          massgas(j) =  massoftype(itypei) * usolarmass ! Msun
          do k=1,ndusttypes
             rhodust(k,j) = dustfrac(k,i)*rhogasi
             massdust(k,j) = dustfrac(k,i)*massgas(j)
          enddo
       endif
    endif
 enddo
 n_SPH = j

 if (sum(massgas) == 0.) then
    write(*,*) ''
    write(*,*) 'Using old phantom dumpfile without gas-to-dust ratio on dust particles.'
    write(*,*) 'Set use_dust_particles = .false. (read_phantom.f90) and compile and run again.'
    write(*,*) 'Exiting...'
    stop
 endif

 write(*,*) ''
 write(*,*) 'Gas mass is      ', sum(massgas)
 write(*,*) 'SPH dust mass is ', sum(massdust)
 write(*,*) ''

 if (.not.lno_internal_energy) then
    if (ndudt == np) then
       write(*,*) "Computing energy input"
       lextra_heating = .true.
       allocate(extra_heating(n_SPH), stat=alloc_status)
       if (alloc_status /=0) then
          write(*,*) "Allocation error in phanton_2_mcfost"
          write(*,*) "Exiting"
       endif

       totlum = 0.
       do i=1,n_SPH
          qtermi = dudt(particle_id(i)) * uWatt
          totlum = totlum + qtermi
          extra_heating(i) = qtermi
       enddo

       write(*,*) "Total energy input = ",real(totlum),' W'
       write(*,*) "Total energy input = ",real(totlum/Lsun),' Lsun'
    endif
 endif
 if (present( T_to_u )) then
   ufac_implicit = T_to_u * uWatt
   ldudt_implicit = .true.
 else
   ufac_implicit = 0.
   ldudt_implicit = .false.
 endif

 write(*,*) "Updating the stellar properties:"
 n_etoiles = 0
 do i=1,nptmass
    if (xyzmh_ptmass(4,i) > 0.0124098) then ! 13 Jupiter masses
       n_etoiles = n_etoiles + 1
    endif
 enddo
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

 return

end subroutine phantom_2_mcfost

!*************************************************************************

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

end module read_phantom
