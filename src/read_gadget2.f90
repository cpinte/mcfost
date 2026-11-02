!******************************************************************************
! Module: read_gadget2
! Purpose: Read classic Gadget-2 binary dumps and GIZMO/Gadget-style HDF5
!          snapshots for the SPH to Voronoi pipeline.
!******************************************************************************
module read_gadget2

  use parameters
  use constants
  use utils
  use messages
  use mcfost_env, only : sp

  implicit none

contains

  !***************************************************************************
  ! Dispatcher: choose HDF5 or binary reader from the filename extension
  !***************************************************************************
  subroutine read_gadget2_file(iunit,filename,x,y,z,h,vx,vy,vz,particle_id,massgas,rhogas, &
       rhodust,ndusttypes,n_SPH,ierr)

    integer,               intent(in) :: iunit
    character(len=*),      intent(in) :: filename
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,vx,vy,vz,rhogas,massgas
    integer,  intent(out), dimension(:),   allocatable :: particle_id
    real(dp), intent(out), dimension(:,:), allocatable :: rhodust
    integer, intent(out) :: ndusttypes,n_SPH,ierr

    integer :: ilen
    logical :: is_hdf5

    ierr = 0
    is_hdf5 = .false.

    ! Detect Gadget/GIZMO HDF5 dumps from the filename suffix
    ilen = index(filename,'.',back=.true.)
    if (ilen > 0) then
       if (filename(ilen:len_trim(filename)) == '.h5' .or. &
            filename(ilen:len_trim(filename)) == '.hdf5') then
          is_hdf5 = .true.
       endif
    endif

    if (is_hdf5) then
       call read_gadget_hdf5_file(filename,x,y,z,h,vx,vy,vz,particle_id,massgas,rhogas, &
            rhodust,ndusttypes,n_SPH,ierr)
    else
       call read_gadget_binary_file(iunit,filename,x,y,z,h,vx,vy,vz,particle_id,massgas, &
            rhogas,rhodust,ndusttypes,n_SPH,ierr)
    endif

    return

  end subroutine read_gadget2_file

  !***************************************************************************
  ! Read a GIZMO / Gadget HDF5 snapshot (PartType0 gas + optional PartType5)
  !***************************************************************************
  subroutine read_gadget_hdf5_file(filename,x,y,z,h,vx,vy,vz,particle_id,massgas,rhogas, &
       rhodust,ndusttypes,n_SPH,ierr)

    use HDF5_utils, only: open_hdf5file, close_hdf5file,   &
         open_hdf5group, close_hdf5group,                   &
         hdf_read_attribute, read_from_hdf5, HID_T

    character(len=*),      intent(in) :: filename
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,vx,vy,vz,rhogas,massgas
    integer,  intent(out), dimension(:),   allocatable :: particle_id
    real(dp), intent(out), dimension(:,:), allocatable :: rhodust
    integer, intent(out) :: ndusttypes,n_SPH,ierr

    integer(HID_T) :: file_id, header_id, group_id
    integer :: npart(6), n_sinks, alloc_status, i
    real(dp) :: mass_table(6), unit_length_cgs, unit_mass_cgs, unit_vel_cgs
    real(dp) :: ulength_au, usolarmass, uvelocity, x0, y0, z0
    real(dp), allocatable :: coords(:,:), vels(:,:), sink_coords(:,:), sink_vels(:,:)
    real(dp), allocatable :: sink_mass(:)
    real(sp), allocatable :: hsml_sp(:), mass_sp(:)
    logical :: got, recentred

    ierr = 0
    ndusttypes = 0
    recentred = .false.

    write(*,*) ''
    write(*,*) ' Opening Gadget/GIZMO HDF5 file: ', trim(filename)

    call open_hdf5file(filename,file_id,ierr)
    if (ierr /= 0) call error("cannot open Gadget/GIZMO HDF5 file: "//trim(filename))

    call open_hdf5group(file_id,'Header',header_id,ierr)
    if (ierr /= 0) call error("cannot open Header group in "//trim(filename))

    ! Header quantities are stored as attributes (Gadget/GIZMO convention)
    call hdf_read_attribute(header_id,"","NumPart_ThisFile",npart)
    call hdf_read_attribute(header_id,"","MassTable",mass_table)
    call hdf_read_attribute(header_id,"","UnitLength_In_CGS",unit_length_cgs)
    call hdf_read_attribute(header_id,"","UnitMass_In_CGS",unit_mass_cgs)
    call hdf_read_attribute(header_id,"","UnitVelocity_In_CGS",unit_vel_cgs)

    call close_hdf5group(header_id,ierr)
    if (ierr /= 0) call error("cannot close Header group in "//trim(filename))

    n_SPH = npart(1)
    n_sinks = npart(6)
    if (n_SPH <= 0) call error("No gas particles (PartType0) in "//trim(filename))

    ulength_au = unit_length_cgs / AU_to_cm
    usolarmass = unit_mass_cgs / Msun_to_g
    ! Match phantom: store velocities in m/s for Voronoi / line transfer
    uvelocity = unit_vel_cgs * cm_to_m

    write(*,*) '  UnitLength_In_CGS = ', real(unit_length_cgs), ' cm (', real(ulength_au), ' AU)'
    write(*,*) '  UnitMass_In_CGS   = ', real(unit_mass_cgs), ' g (', real(usolarmass), ' Msun)'
    write(*,*) '  UnitVelocity_In_CGS = ', real(unit_vel_cgs), ' cm/s (', real(uvelocity), ' m/s)'
    write(*,*) '  Gas particles     = ', n_SPH
    write(*,*) '  Sink particles    = ', n_sinks

    ! GIZMO stores Coordinates/Velocities as C (npart,3); the Fortran HDF5
    ! interface reports that as (3,npart). Match phantom: read into (3,n).
    allocate(x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),vx(n_SPH),vy(n_SPH),vz(n_SPH), &
         massgas(n_SPH),rhogas(n_SPH),particle_id(n_SPH),coords(3,n_SPH),vels(3,n_SPH), &
         hsml_sp(n_SPH),stat=alloc_status)
    if (alloc_status /= 0) call error("Allocation error reading Gadget/GIZMO gas arrays")

    call open_hdf5group(file_id,'PartType0',group_id,ierr)
    if (ierr /= 0) call error("cannot open PartType0 in "//trim(filename))

    ! Coordinates are double precision in this dump; h and masses are single
    call read_from_hdf5(coords,'Coordinates',group_id,got,ierr)
    if (.not.got .or. ierr /= 0) call error("cannot read PartType0/Coordinates")

    call read_from_hdf5(vels,'Velocities',group_id,got,ierr)
    if (.not.got .or. ierr /= 0) call error("cannot read PartType0/Velocities")

    call read_from_hdf5(hsml_sp,'SmoothingLength',group_id,got,ierr)
    if (.not.got .or. ierr /= 0) call error("cannot read PartType0/SmoothingLength")

    if (mass_table(1) == 0.0_dp) then
       allocate(mass_sp(n_SPH),stat=alloc_status)
       if (alloc_status /= 0) call error("Allocation error reading Gadget/GIZMO Masses")
       call read_from_hdf5(mass_sp,'Masses',group_id,got,ierr)
       if (.not.got .or. ierr /= 0) call error("cannot read PartType0/Masses")
       massgas(:) = real(mass_sp(:),kind=dp) * usolarmass
       deallocate(mass_sp)
    else
       massgas(:) = mass_table(1) * usolarmass
    endif

    call close_hdf5group(group_id,ierr)
    if (ierr /= 0) call error("cannot close PartType0 in "//trim(filename))

    x(:) = coords(1,:)
    y(:) = coords(2,:)
    z(:) = coords(3,:)
    vx(:) = vels(1,:)
    vy(:) = vels(2,:)
    vz(:) = vels(3,:)
    h(:) = real(hsml_sp(:),kind=dp)
    deallocate(coords,vels,hsml_sp)

    ! Optional sinks / black holes (PartType5); recenter on the first one
    if (n_sinks > 0) then
       allocate(sink_coords(3,n_sinks),sink_vels(3,n_sinks),sink_mass(n_sinks),stat=alloc_status)
       if (alloc_status /= 0) call error("Allocation error reading Gadget/GIZMO sinks")

       call open_hdf5group(file_id,'PartType5',group_id,ierr)
       if (ierr /= 0) call error("cannot open PartType5 in "//trim(filename))

       call read_from_hdf5(sink_coords,'Coordinates',group_id,got,ierr)
       if (.not.got .or. ierr /= 0) call error("cannot read PartType5/Coordinates")

       call read_from_hdf5(sink_vels,'Velocities',group_id,got,ierr)
       if (.not.got .or. ierr /= 0) call error("cannot read PartType5/Velocities")

       if (mass_table(6) == 0.0_dp) then
          allocate(mass_sp(n_sinks),stat=alloc_status)
          if (alloc_status /= 0) call error("Allocation error reading PartType5/Masses")
          call read_from_hdf5(mass_sp,'Masses',group_id,got,ierr)
          if (.not.got .or. ierr /= 0) call error("cannot read PartType5/Masses")
          sink_mass(:) = real(mass_sp(:),kind=dp)
          deallocate(mass_sp)
       else
          sink_mass(:) = mass_table(6)
       endif

       call close_hdf5group(group_id,ierr)
       if (ierr /= 0) call error("cannot close PartType5 in "//trim(filename))

       ! Recentre so the first sink sits at the origin (mcfost disc convention)
       x0 = sink_coords(1,1)
       y0 = sink_coords(2,1)
       z0 = sink_coords(3,1)
       x(:) = x(:) - x0
       y(:) = y(:) - y0
       z(:) = z(:) - z0
       sink_coords(1,:) = sink_coords(1,:) - x0
       sink_coords(2,:) = sink_coords(2,:) - y0
       sink_coords(3,:) = sink_coords(3,:) - z0
       recentred = .true.

       write(*,*) '  Recentring on first sink at code coords: ', real(x0), real(y0), real(z0)

       ! Update stellar positions / masses, preserving para-file T, radius, spectrum
       if (n_sinks /= n_stars) then
          write(*,*) '  Updating number of stars from', n_stars, 'to', n_sinks
          call resize_stars(n_sinks, lkeep_all=.true.)
       endif

       do i = 1, n_stars
          star(i)%x = sink_coords(1,i) * ulength_au
          star(i)%y = sink_coords(2,i) * ulength_au
          star(i)%z = sink_coords(3,i) * ulength_au
          star(i)%vx = sink_vels(1,i) * uvelocity
          star(i)%vy = sink_vels(2,i) * uvelocity
          star(i)%vz = sink_vels(3,i) * uvelocity
          star(i)%M = sink_mass(i) * usolarmass
          write(*,*) '  Sink', i, 'mass =', real(star(i)%M), 'Msun'
       enddo

       deallocate(sink_coords,sink_vels,sink_mass)
    endif

    call close_hdf5file(file_id,ierr)
    if (ierr /= 0) call error("cannot close Gadget/GIZMO HDF5 file: "//trim(filename))

    ! Convert lengths to AU and velocities to m/s; masses already in Msun
    x(:) = x(:) * ulength_au
    y(:) = y(:) * ulength_au
    z(:) = z(:) * ulength_au
    h(:) = h(:) * ulength_au
    vx(:) = vx(:) * uvelocity
    vy(:) = vy(:) * uvelocity
    vz(:) = vz(:) * uvelocity
    rhogas(:) = 0.0_dp

    do i = 1, n_SPH
       particle_id(i) = i
    enddo

    write(*,*) '  Total gas mass    = ', real(sum(massgas)), 'Msun'
    if (recentred) then
       write(*,*) '  Coordinates recentred on first sink (AU)'
    endif
    write(*,*) ' Closing Gadget/GIZMO HDF5 file'
    write(*,*) ''

    return

  end subroutine read_gadget_hdf5_file

  !***************************************************************************
  ! Read a classic Gadget-2 unformatted binary dump
  !***************************************************************************
  subroutine read_gadget_binary_file(iunit,filename,x,y,z,h,vx,vy,vz,particle_id,massgas, &
       rhogas,rhodust,ndusttypes,n_SPH,ierr)

    integer,               intent(in) :: iunit
    character(len=*),      intent(in) :: filename
    real(dp), intent(out), dimension(:),   allocatable :: x,y,z,h,vx,vy,vz,rhogas,massgas
    integer,  intent(out), dimension(:),   allocatable :: particle_id
    real(dp), intent(out), dimension(:,:), allocatable :: rhodust
    integer, intent(out) :: ndusttypes,n_SPH,ierr

    integer :: i, j, testlen
    logical :: exists
    character*4 :: blockid
    ! npart(0) = SPH gas; npart(5) = N-body / sink particles
    integer :: npart(0:5), ntot, len, flag_sfr, flag_feedback, unused(33), nall(0:5)
    integer :: n_mass, n_sinks, alloc_status, i_star, i_mass, i_part, k
    double precision :: massarr(0:5), time, redshift
    real, allocatable :: pos(:,:), vel(:,:), mass_block(:), u(:), rho(:), hsml(:)
    integer, allocatable :: id(:)
    real, allocatable :: nh(:), ne(:)

    ! Legacy hard-coded units for older disc binary dumps
    real, parameter :: ulength =  3.29120e+12
    real, parameter :: ulength_au = ulength / AU_to_cm
    real, parameter :: uspeed = 1.20524e+07 ! cm/s
    real, parameter :: uvelocity = uspeed * cm_to_m ! m/s
    real, parameter :: udens =  2.00932e-04
    real, parameter :: umass = udens * ulength**3
    real, parameter :: usolarmass = umass / Msun_to_g

    ierr = 0
    ndusttypes = 0

    inquire(file = filename(1:len_trim(filename)), exist = exists)
    if (.not.exists) call error("File "//trim(filename)//" does not exist")

    write(*,*) ''
    write(*,*) '      Opening Gadget-2 file: ', filename(1:len_trim(filename))
    open(unit = iunit, file = filename(1:len_trim(filename)), form = 'unformatted')

    read(iunit) blockid, len
    testlen = 264
    if (len /= testlen) then
       write(*,*) "ERROR: cannot read header of Gadget dump file"
       write(*,*) 'Incorrect length: ', len
       write(*,*) '       Should be: ', testlen
       write(*,*) '        On block: ', blockid
       call error("Invalid Gadget-2 binary header")
    end if
    read(iunit) npart, massarr, time, redshift, flag_sfr, flag_feedback, nall, unused

    write(*,*) ''
    write(*,*) '        Total gas particles: ', npart(0)
    write(*,*) '          Gas particle mass: ', massarr(0)
    if (massarr(0) == 0.) then
       write(*,*) "            Total gas mass : individual masses (MASS block)"
    else
       write(*,*) "            Total gas mass : ", npart(0) * massarr(0) * usolarmass
    endif
    write(*,*) '           particle #5 mass: ', massarr(:)
    write(*,*) '     Total N-body particles: ', npart(5)
    write(*,*) '             Time in orbits: ', time / (2. * pi)
    write(*,*) ''

    n_SPH = npart(0)
    ntot = sum(npart)

    ! Particles with zero MassTable entry store individual masses in the MASS block, in
    ! the order of the types. Stars: the non-gas particles with individual masses (the
    ! gas was counted as stars when it had individual masses)
    n_mass = 0
    do j = 0, 5
       if (massarr(j) == 0.) n_mass = n_mass + npart(j)
    end do
    n_sinks = 0
    do j = 1, 5
       if (massarr(j) == 0.) n_sinks = n_sinks + npart(j)
    end do

    allocate(pos(3,ntot), vel(3,ntot), id(ntot), mass_block(n_mass), &
         u(npart(0)), ne(npart(0)), nh(npart(0)))
    allocate(hsml(npart(0)), rho(npart(0)))

    loop : do
       read(iunit) blockid, len
       if (blockid == 'POS ') then
          testlen = (4 * 3 * ntot) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) pos
          write(*,*) '                 Block read: ', blockid, '        (Positions)'
       else if (blockid == 'VEL ') then
          testlen = (4 * 3 * ntot) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) vel
          write(*,*) '                 Block read: ', blockid, '        (Velocities)'
       else if (blockid == 'ID  ') then
          testlen = (4 * ntot) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) id
          write(*,*) '                 Block read: ', blockid, '        (IDs)'
       else if (blockid == 'MASS') then
          testlen = (4 * n_mass) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) mass_block
          write(*,*) '                 Block read: ', blockid, '        (Masses)'
       else if (blockid == 'U   ') then
          testlen = (4 * npart(0)) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) u
          write(*,*) '                 Block read: ', blockid, '        (Internal energies)'
       else if (blockid == 'RHO ') then
          testlen = (4 * npart(0)) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) rho
          write(*,*) '                 Block read: ', blockid, '        (Densities)'
       else if (blockid == 'NE  ') then
          testlen = (4 * npart(0)) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) ne
          write(*,*) '                 Block read: ', blockid, '        (Electron abundences: not used)'
       else if (blockid == 'NH  ') then
          testlen = (4 * npart(0)) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) nh
          write(*,*) '                 Block read: ', blockid, '        (Hydrogen fractions: not used)'
       else if (blockid == 'HSML') then
          testlen = (4 * npart(0)) + 8
          if (len /= testlen) then
             ierr = 3
             exit loop
          end if
          read(iunit) hsml
          write(*,*) '                 Block read: ', blockid, '        (Smoothing lengths)'
          exit
       else
          write(*,*) "Block skipped ", blockid
       endif
    end do loop
    write(*,*) ''

    write(*,*) '      Closing Gadget-2 file: ', filename(1:len_trim(filename))
    write(*,*) ''
    close(iunit)

    if (ierr == 3) then
       write(*,*) ''
       write(*,*) 'Error:'
       write(*,*) 'Incorrect length: ', len
       write(*,*) '       Should be: ', testlen
       write(*,*) '        On block: ', blockid
       call error("Invalid block length in Gadget-2 binary dump")
    endif

    alloc_status = 0
    allocate(x(n_SPH),y(n_SPH),z(n_SPH),h(n_SPH),vx(n_SPH),vy(n_SPH),vz(n_SPH), &
         massgas(n_SPH),rhogas(n_SPH),particle_id(n_SPH),stat=alloc_status)
    if (alloc_status /= 0) call error("Allocation error in gadget_2_mcfost")

    x(:) = pos(1,1:n_SPH) * ulength_au
    y(:) = pos(2,1:n_SPH) * ulength_au
    z(:) = pos(3,1:n_SPH) * ulength_au
    h(:) = hsml(1:n_SPH) * ulength_au
    vx(:) = vel(1,1:n_SPH) * uvelocity
    vy(:) = vel(2,1:n_SPH) * uvelocity
    vz(:) = vel(3,1:n_SPH) * uvelocity
    if (massarr(0) == 0.) then ! gas particles with individual masses: first in the MASS block
       massgas(:) = mass_block(1:n_SPH) * usolarmass
    else
       massgas(:) = massarr(0) * usolarmass
    endif
    rhogas(:) = 0.0_dp

    do i = 1, n_SPH
       particle_id(i) = i
    enddo

    ! As in the HDF5 reader: without star in the file, the stars of the parameter file are
    ! used; otherwise only the positions, velocities and masses are updated (para-file T,
    ! radius and spectrum preserved, as with -fix_star)
    write(*,*) "Found", n_sinks, "stars in the Gadget-2 file"
    if (n_sinks > 0) then
       write(*,*) "Updating the stellar positions, velocities and masses"
       if (n_sinks /= n_stars) then
          write(*,*) '  Updating number of stars from', n_stars, 'to', n_sinks
          call resize_stars(n_sinks, lkeep_all=.true.)
       endif

       ! particles in the order of the types (positions), masses in the MASS block for the
       ! types with individual masses only (after the gas if it has individual masses)
       i_star = 0 ; i_part = npart(0) ; i_mass = 0
       if (massarr(0) == 0.) i_mass = npart(0)
       do j=1, 5
          do k=1, npart(j)
             i_part = i_part + 1
             if (massarr(j) == 0.) then
                i_star = i_star + 1 ; i_mass = i_mass + 1
                star(i_star)%x = pos(1,i_part) * ulength_au
                star(i_star)%y = pos(2,i_part) * ulength_au
                star(i_star)%z = pos(3,i_part) * ulength_au
                star(i_star)%vx = vel(1,i_part) * uvelocity
                star(i_star)%vy = vel(2,i_part) * uvelocity
                star(i_star)%vz = vel(3,i_part) * uvelocity
                star(i_star)%M =  mass_block(i_mass) * usolarmass
                write(*,*) i_star, star(i_star)%M
             endif
          enddo
       enddo
    endif

    return

  end subroutine read_gadget_binary_file

end module read_gadget2
