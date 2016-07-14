module read_gadget2

  use parametres
  use constantes
  use utils

  implicit none

contains

  subroutine read_gadget2_file(iunit,filename,x,y,z,massgas,rhogas,rhodust,ndusttypes,n_SPH,ierr)

    use prop_star

    integer,               intent(in) :: iunit
    character(len=*),      intent(in) :: filename
    real(db), intent(out), dimension(:),   allocatable :: x,y,z,rhogas,massgas
    real(db), intent(out), dimension(:,:), allocatable :: rhodust
    integer, intent(out) :: ndusttypes,n_SPH

    ! dump number, error flag
    integer :: i, j, ierr, testlen
    ! logical operator for file check
    logical :: exists
    ! id for reading data blocks
    character*4 :: blockid
    ! number of particles per type: npart(0) gives SPH particles, npart(5) gives N-body particles.
    ! ntot is total particle number
    ! others are required for data read but not used
    integer :: npart(0:5), ntot, len, flag_sfr, flag_feedback, unused(33), nall(0:5), n_stars, alloc_status
    ! massarr gives particle mass per type (as with npart) for particle
    ! types of uniform mass (only type 0 SPH particles here)
    double precision :: massarr(0:5), time, redshift
    ! position, velocity, mass, internal energy, density and smoothing length arrays
    real, allocatable :: pos(:,:), vel(:,:), mass_stars(:), u(:), rho(:), hsml(:)
    ! particle id array
    integer, allocatable :: id(:)
    ! other SPH data (unused in current simulations)
    real, allocatable :: nh(:), ne(:)


    ! For positions in cm, multiply code units by:    3.29120e+12
    ! For velocities in cm/s, multiply code units by: 1.20524e+07
    ! For densities in g/cm^3, multiply code units by:        2.00932e-04
    real, parameter :: Msun = 1.9892e+33 ! en g
    real, parameter :: AU = 1.49598e+13 ! en cm
    real, parameter :: ulength =  3.29120e+12
    real, parameter :: ulength_au = ulength / AU_to_cm ! factor ulength to convert in cm
    real, parameter :: uspeed = 1.20524e+07
    real, parameter :: udens =  2.00932e-04
    real, parameter :: umass = udens * ulength**3
    real, parameter :: usolarmass = umass / Msun_to_g  ! for mass in g, multiply the code by umass

    ierr = 0
    ndusttypes = 1 ! only gas for gadget

    ! check dump file exists and open
    inquire(file = filename(1:len_trim(filename)), exist = exists)
    if(exists) then
       write(*,*) ''
       write(*,*) '      Opening Gadget-2 file: ', filename(1:len_trim(filename))
       open(unit = 1, file = filename(1:len_trim(filename)), form = 'unformatted')
       ! read header data
       read(1) blockid, len
       testlen = 264
       if(len /= testlen)then
          write(*,*) "ERROR: cannot read header of Gadget dump file"
          write(*,*) 'Incorrect length: ', len
          write(*,*) '       Should be: ', testlen
          write(*,*) '        On block: ', blockid
          write(*,*) "Exiting"
          stop
       end if
       read(1) npart, massarr, time, redshift, flag_sfr, flag_feedback, nall, unused
       ! print selected header data
       write(*,*) ''
       write(*,*) '        Total gas particles: ', npart(0)
       write(*,*) '          Gas particle mass: ', massarr(0)
       write(*,*) "            Total gas mass : ", npart(0) * massarr(0) * usolarmass
       write(*,*) '           particle #5 mass: ', massarr(:)
       write(*,*) '     Total N-body particles: ', npart(5)
       write(*,*) '             Time in orbits: ', time / (2. * pi)
       write(*,*) ''
    else
       write(*,*) "File "//trim(filename)//" does not exist"
       write(*,*) "Exiting"
       stop
    end if

    n_SPH = npart(0)
    n_stars = npart(5) ! number of stellar objects

    ! set array sizes for reading
    ntot = sum(npart)
    !write(*,*) "Ntot =", ntot

    n_stars = 0
    do j = 0, 5
       if(massarr(j) == 0.) n_stars = n_stars + npart(j)
    end do

    ! allocate arrays for particle data (and deallocate if necessary)
    allocate(pos(3,ntot), vel(3,ntot), id(ntot), mass_stars(n_stars), u(npart(0)), ne(npart(0)), nh(npart(0)))
    allocate(hsml(npart(0)), rho(npart(0)))

    ! cycle through blocks, check length, and read data
    loop : do
       read(1) blockid, len
       if(blockid == 'POS ')then
          testlen = (4 * 3 * ntot) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) pos
          write(*,*) '                 Block read: ', blockid, '        (Positions)'
       else if(blockid == 'VEL ')then
          testlen = (4 * 3 * ntot) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) vel
          write(*,*) '                 Block read: ', blockid, '        (Velocities)'
       else if(blockid == 'ID  ')then
          testlen = (4 * ntot) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) id
          write(*,*) '                 Block read: ', blockid, '        (IDs)'
       else if(blockid == 'MASS')then
          testlen = (4 * n_stars) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) mass_stars
          write(*,*) '                 Block read: ', blockid, '        (Masses)'
       else if(blockid == 'U   ')then
          testlen = (4 * npart(0)) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) u
          write(*,*) '                 Block read: ', blockid, '        (Internal energies)'
       else if(blockid == 'RHO ')then
          testlen = (4 * npart(0)) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) rho
          write(*,*) '                 Block read: ', blockid, '        (Densities)'
       else if(blockid == 'NE  ')then
          testlen = (4 * npart(0)) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) ne
          write(*,*) '                 Block read: ', blockid, '        (Electron abundences: not used)'
       else if(blockid == 'NH  ')then
          testlen = (4 * npart(0)) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) nh
          write(*,*) '                 Block read: ', blockid, '        (Hydrogen fractions: not used)'
       else if(blockid == 'HSML')then
          testlen = (4 * npart(0)) + 8
          if(len /= testlen)then
             ierr = 3
             exit loop
          end if
          read(1) hsml
          write(*,*) '                 Block read: ', blockid, '        (Smoothing lengths)'
          exit
       else
          write(*,*) "Block skipped ", blockid
       endif
    end do loop
    write(*,*) ''

    ! Close file
    write(*,*) '      Closing Gadget-2 file: ', filename(1:len_trim(filename))
    write(*,*) ''
    close(1)

    ! Check errors
    if (ierr == 3) then
       write(*,*) ''
       write(*,*) 'Error:'
       write(*,*) 'Incorrect length: ', len
       write(*,*) '       Should be: ', testlen
       write(*,*) '        On block: ', blockid
       write(*,*) "Exiting"
       stop
    endif

    alloc_status = 0
    allocate(x(n_SPH),y(n_SPH),z(n_SPH),massgas(n_SPH),rhogas(n_SPH),rhodust(ndusttypes,n_SPH), stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error in phanton_2_mcfost"
       write(*,*) "Exiting"
    endif

    x(:) = pos(1,1:n_SPH) * ulength_au
    y(:) = pos(2,1:n_SPH) * ulength_au
    z(:) = pos(3,1:n_SPH) * ulength_au
    massgas(:) = massarr(0) * usolarmass ! en Msun

    rhogas(:) = 0.
    rhodust(:,:) = 0.

    write(*,*) "Found", n_stars, "stars in the Gadget-2 file"
    if (n_stars > 0) then
       write(*,*) "Updating the stellar properties"
       if (allocated(etoile)) deallocate(etoile)
       allocate(etoile(n_stars))

       do i=1,n_stars
          etoile(i)%x = pos(1,n_SPH + i) * ulength_au
          etoile(i)%y = pos(2,n_SPH + i) * ulength_au
          etoile(i)%z = pos(3,n_SPH + i) * ulength_au

          etoile(i)%M =  mass_stars(i) * usolarmass
          write(*,*) i, etoile(i)%M
       enddo
    endif

    return

  end subroutine read_gadget2_file

end module read_gadget2
