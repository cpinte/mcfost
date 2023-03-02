module MRW
  use mcfost_env, only : dp
  use constantes, only : pi
  use utils, only : interp

  implicit none

  integer, parameter :: n = 10000
  real(dp), dimension(n), save :: zeta, y_MRW

contains

 subroutine initialize_cumulative_zeta()
  ! Compute equation (7) of Min et al (2009):
  ! P(t) = 2 * Sum_{n=1}^{infinity} (-1)^(n+1) * y^(n^2)
  !
  ! The function requires high n values close to y=1,
  ! so we pre-compute it here and we will then interpolate

  implicit none

  integer :: i,j
  real(kind=8) :: term

  zeta = 0._dp

  do i=1,n
     y_MRW(i) = real(i-1, dp)/real(n-1, dp)

     if(i==n) then
        zeta(i) = 0.5_dp
     else
        j = 0
        do
           j = j + 1
           term = y_MRW(i)**(j**2)
           if(term == 0._dp) exit
           if(mod(j, 2)==0) then
              zeta(i) = zeta(i) - term
           else
              zeta(i) = zeta(i) + term
           endif
        enddo
     endif
  enddo

  zeta = zeta * 2._dp

  return

 end subroutine initialize_cumulative_zeta

 !----------------------------------------------

 real(dp) function sample_zeta(id)

   use naleat, only : stream
#include "sprng_f.h"

   integer, intent(in) :: id

   real(dp) :: zeta_random

   zeta_random = sprng(stream(id))
   sample_zeta = interp(y_MRW, zeta, zeta_random)

 end function sample_zeta

 !----------------------------------------------

 subroutine make_MRW_step(id,icell,lambda, x,y,z,R0, E)

   integer, intent(in) :: id, icell
   integer, intent(inout) :: lambda
   real(kind=dp), intent(inout) :: x,y,z
   real(kind=dp), intent(in) :: R0, E

   real(dp) :: u,v,w, ct

   ! Place photon randomly on sphere of radius R0 around current position
   call random_isotropic_direction(id, u,v,w)
   x = x + u*R0
   y = y + v* R0
   z = z + w*R0

   ! Select new direction : isotropic or emitting sphere ?

   ! sample P0
   y = sample_zeta(id)

   ! Solve Eq. 8 of Min et al 2009
   ct = -log(y) / diff_coeff * (R0/pi)**2.

   ! Deposit energy using Planck mean opacity

   ! Update temperature and dissuion coeeficient ?

   ! Select new wavelength

   return

 end subroutine make_MRW_step


end module MRW
