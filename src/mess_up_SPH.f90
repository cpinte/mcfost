module mess_up_SPH

  implicit none

  public :: mask_Hill_sphere, randomize_azimuth, randomize_gap

  private

  integer, parameter :: i_star = 1
  integer, parameter :: dp = selected_real_kind(p=13,r=200)
  real(dp), parameter :: pi = 4._dp*atan(1._dp)

contains

  subroutine mask_Hill_sphere(np, nptmass, xyzh, xyzmh_ptmass,udist,mask)
    ! Create and return a mask with the particles that are not within the Hill spheres
    ! of all the planets in a dump
    ! mask is set to True is particle is inside the Hill sphere

    integer, intent(in) :: np, nptmass
    real(kind=dp), dimension(4,np), intent(in) :: xyzh
    real(kind=dp), dimension(:,:), intent(in) :: xyzmh_ptmass
    real(kind=dp), intent(in) :: udist

    logical, dimension(np), intent(out) :: mask

    integer :: i_planet, i, n_delete
    real(kind=dp) :: d2, r_Hill2, r_hill, dx, dy, dz

    ! We assume that the 1st sink particle is the actual star
    ! and the following sink particles are planets
    mask(:) = .false.
    do i_planet=2, nptmass
       n_delete = 0

       r_Hill2 = Hill_radius2(nptmass, xyzmh_ptmass, i_planet)
       r_Hill = sqrt(r_Hill2)
       write(*,*) "Sink particle #", i_planet, "Hill radius =", r_Hill * udist, "au"

       particle_loop : do i=1, np
          ! We first exlude particles that are not with a cube around the sink particle
          dx = xyzh(1,i) - xyzmh_ptmass(1,i_planet)
          if (dx > r_Hill) cycle particle_loop
          dy = xyzh(2,i) - xyzmh_ptmass(2,i_planet)
          if (dy > r_Hill) cycle particle_loop
          dz = xyzh(3,i) - xyzmh_ptmass(3,i_planet)
          if (dz > r_Hill) cycle particle_loop

          ! We then test on the sphere itself
          d2 = dx**2 + dy**2 + dz**2
          if (d2 < r_Hill2) then ! particle is in Hill sphere
             mask(i) = .true.
             n_delete = n_delete + 1
          endif
       enddo particle_loop

       write(*,*) "Deleting", n_delete, "particles in Hill sphere of sink particle #", i_planet
    enddo

    return

  end subroutine mask_Hill_sphere

  !*********************************************************
  ! Compute the square of distance between star and sink particle i_planet

  function d2_from_star(nptmass, xyzmh_ptmass, i_planet) result(d2)
    integer, intent(in) :: nptmass, i_planet
    real(kind=dp), dimension(:,:), intent(in) :: xyzmh_ptmass
    real(kind=dp) :: d2

    d2 = (xyzmh_ptmass(1,i_planet) - xyzmh_ptmass(1,i_star))**2 + &
         (xyzmh_ptmass(2,i_planet) - xyzmh_ptmass(2,i_star))**2 + &
         (xyzmh_ptmass(3,i_planet) - xyzmh_ptmass(3,i_star))**2

  end function d2_from_star

  !*********************************************************

  function Hill_radius2(nptmass, xyzmh_ptmass, i_planet)
    ! Compute the square of the Hill radius for sink particle i_planet
    ! Assuming the star is sink particle #1
    ! Units : [length code units**2]

    integer, intent(in) :: nptmass, i_planet
    real(kind=dp), dimension(:,:), intent(in) :: xyzmh_ptmass
    real(kind=dp) :: d2, Hill_radius2

    d2 = d2_from_star(nptmass,xyzmh_ptmass,i_planet)
    Hill_radius2 = d2 * (xyzmh_ptmass(4,i_planet) / (3*xyzmh_ptmass(4,i_star)))**(2./3)

    return

  end function Hill_radius2

  !*********************************************************

  subroutine randomize_azimuth(np, xyzh, vxyzu, mask)
    ! Randomly rotate all the particles around the z axis
    ! Only rotates particles that are masked (ie where mask == .true.)

    integer, intent(in) :: np
    real(kind=dp), dimension(4,np), intent(inout) :: xyzh, vxyzu
    logical, dimension(np), intent(in), optional :: mask
    real(kind=dp) :: cos_phi, sin_phi, phi, x_tmp, y_tmp
    integer :: i

    particle_loop : do i=1, np
       if (present(mask)) then
          if (.not.mask(i)) cycle particle_loop
       endif
       call random_number(phi)
       phi = 2.*pi*phi
       cos_phi = cos(phi) ; sin_phi = sin(phi)

       !-- position
       x_tmp = xyzh(1,i) * cos_phi + xyzh(2,i) * sin_phi
       y_tmp = - xyzh(1,i) * sin_phi +  xyzh(2,i) * cos_phi
       xyzh(1,i) = x_tmp ;  xyzh(2,i) = y_tmp

       !-- velocities
       x_tmp = vxyzu(1,i) * cos_phi + vxyzu(2,i) * sin_phi
       y_tmp = -vxyzu(1,i) * sin_phi + vxyzu(2,i) * cos_phi
       vxyzu(1,i) = x_tmp ; vxyzu(2,i) = y_tmp
    enddo particle_loop

    return

  end subroutine randomize_azimuth

  !*********************************************************

  subroutine randomize_gap(np, nptmass, xyzh, vxyzu, xyzmh_ptmass,udist)

    integer, intent(in) :: np, nptmass
    real(kind=dp), dimension(4,np), intent(inout) :: xyzh, vxyzu
    real(kind=dp), dimension(:,:), intent(in) :: xyzmh_ptmass
    real(kind=dp), intent(in) :: udist

    logical, dimension(np) :: mask ! true for the particles in the gaps

    integer :: i_planet, i, n_delete
    real(kind=dp) :: r_Hill2, r_hill, r_planet, r_minus2, r_plus2, r2

    ! We assume that the 1st sink particle is the actual star
    ! and the following sink particles are planets
    mask(:) = .false.
    do i_planet=2, nptmass
       n_delete = 0

       r_Hill2 = Hill_radius2(nptmass, xyzmh_ptmass, i_planet)
       r_Hill = sqrt(r_Hill2)
       write(*,*) "Sink particle #", i_planet, "Hill radius =", r_Hill * udist, "au"

       r_planet = sqrt((xyzmh_ptmass(1,i_planet) - xyzmh_ptmass(1,i_star))**2 + &
                       (xyzmh_ptmass(2,i_planet) - xyzmh_ptmass(2,i_star))**2)
       r_plus2 = (r_planet + 3.*r_Hill)**2
       r_minus2 = (r_planet - 3.*r_Hill)**2

       particle_loop : do i=1, np
          r2 = (xyzh(1,i)-xyzmh_ptmass(1,i_star))**2 + (xyzh(2,i)-xyzmh_ptmass(2,i_star))**2

          if (r2 < r_plus2) then
             if (r2 > r_minus2) then
                mask(i) = .true.
                n_delete = n_delete + 1
             endif
          end if
       enddo particle_loop

       write(*,*) "Randomizing", n_delete, "particles in gap of sink particle #", i_planet
    enddo

    call randomize_azimuth(np, xyzh, vxyzu, mask)

    return

  end subroutine randomize_gap

  !*********************************************************

end module mess_up_SPH
