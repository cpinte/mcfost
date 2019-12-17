module mess_up_SPH

  implicit none

  public :: delete_Hill_sphere, randomize_azimuth

  private

  integer, parameter :: dp = selected_real_kind(p=13,r=200)

  contains

  subroutine delete_Hill_sphere(np, nptmass, xyzh, xyzmh_ptmass,udist,mask)

    integer, intent(in) :: np, nptmass
    real(kind=dp), dimension(4,np), intent(inout) :: xyzh
    real(kind=dp), dimension(5,nptmass), intent(inout) :: xyzmh_ptmass
    real(kind=dp), intent(in) :: udist

    logical, dimension(np), intent(out) :: mask

    integer :: istar, i, n_delete
    real(kind=dp) :: d2, r_Hill2, r_hill, dx, dy, dz

    ! We assume that the 1st sink particle is the actual star
    ! and the following sink particles are planets
    mask(:) = .true.
    do istar=2, nptmass
       n_delete = 0

       d2 = (xyzmh_ptmass(1,istar) - xyzmh_ptmass(1,1))**2 + (xyzmh_ptmass(3,istar) - xyzmh_ptmass(3,1))**2 + (xyzmh_ptmass(3,istar) - xyzmh_ptmass(3,1))**2
       r_Hill2 = d2 * (xyzmh_ptmass(4,istar) / (3*xyzmh_ptmass(4,istar)))**(2./3)
       r_Hill = sqrt(r_Hill2)

       write(*,*) "Sink particle #", istar, "Hill radius =", r_Hill * udist, "au"

       particle_loop : do i=1, np
          dx = xyzh(1,i) - xyzmh_ptmass(1,istar)
          if (dx > r_Hill) cycle particle_loop
          dy = xyzh(2,i) - xyzmh_ptmass(2,istar)
          if (dy > r_Hill) cycle particle_loop
          dz = xyzh(3,i) - xyzmh_ptmass(3,istar)
          if (dz > r_Hill) cycle particle_loop

          d2 = dx**2 + dy**2 + dz**2
          if (d2 < r_Hill2) then ! particle is in Hill sphere
             mask(i) = .false.
             n_delete = n_delete + 1
          endif
       enddo particle_loop

       write(*,*) "Deleting", n_delete, "particles in Hill sphere of sink particle #", istar
    enddo

    return

  end subroutine delete_Hill_sphere

  !*********************************************************

  subroutine randomize_azimuth(np, xyzh, vxyzu)

    integer, intent(in) :: np
    real(kind=dp), dimension(np,4), intent(inout) :: xyzh, vxyzu

    real(kind=dp) :: cos_phi, sin_phi, phi, x_tmp, y_tmp
    integer :: i

    particle_loop : do i=1, np
       call random_number(phi)
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




end module mess_up_SPH
