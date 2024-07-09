module MRW
  use mcfost_env, only : dp
  use constantes, only : pi
  use utils, only : interp

  implicit none

  integer, parameter :: n = 10000
  real(dp), dimension(n), save :: zeta, y_MRW

  real, parameter :: gamma_MRW = 2.0
  real(dp), parameter :: cst_ct = 3 / pi**2

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

 subroutine make_MRW_step(id,icell, x,y,z,E, d, rec_Planck_opacity)

   use naleat, only : random_isotropic_direction

   integer, intent(in) :: id, icell
   real(kind=dp), intent(inout) :: x,y,z
   real(kind=dp), intent(in) :: E, d, rec_Planck_opacity

   real(dp) :: u,v,w, ct

   ! Place photon randomly on sphere of radius R0 around current position
   call random_isotropic_direction(id, u,v,w)
   x = x + u*d
   y = y + v*d
   z = z + w*d

   ! Select new direction : isotropic or emitting sphere ?

   ! sample P0
   y = sample_zeta(id)

   ! Solve Eq. 8 of Min et al 2009
   ! D = 1/(3*rec_Planck_opacity)
   ! cst_ct = 3/pi**3
   ct = -log(y) * cst_ct * rec_Planck_opacity * d**2.


   ! Steps that needs to be done
   ! Deposit energy using Planck mean opacity

   ! Update temperature and diffusion coefficient ?


   ! Only at end of MRW, to switch to MC
   ! Select new wavelength

   ! Select new photon direction


   return

 end subroutine make_MRW_step


end module MRW
