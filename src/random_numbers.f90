module naleat
!**********************************************************
! Parametres des nombres aleatoires
!**********************************************************
! SPRNG renvoie un nbre aleatoire dans [0,1[
! faire attention a conversion dp -> sp sur machine 64 bits
! qui peut donner 1.0 -> BUG
!**********************************************************

  use mcfost_env, only : dp
  use constantes, only : pi

  implicit none

  save

#include "sprng_f.h"

  integer :: seed = 269753

  integer :: gtype = 2
  ! 0 -> Additive Lagged Fibonacci Generator
  ! 1 -> 48 bit Linear Congruential Generator with Prime Addend
  ! 2 -> 64 bit Linear Congruential Generator with Prime Addend
  ! 3 -> Combined multiple recursive generator
  ! 4 -> Multiplicative Lagged Fibonacci Generator

  SPRNG_POINTER, dimension(:), allocatable :: stream

contains

  subroutine random_isotropic_direction(id, u,v,w)
    ! Generate a randomly distributed vector (u,w,w)

    integer, intent(in) :: id
    real(kind=dp), intent(out) :: u,v,w

    real(kind=dp) :: uv, phi
    real :: rand

    rand = sprng(stream(id))
    w = 2.0_dp * rand - 1.0_dp
    uv = sqrt (1.0 - w*w)
    rand = sprng(stream(id))
    phi = pi * (2.0_dp * rand - 1.0_dp)
    u = uv * cos(phi)
    v = uv * sin(phi)

    return

  end subroutine random_isotropic_direction

end module naleat
