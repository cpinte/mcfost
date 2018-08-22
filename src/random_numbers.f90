module naleat
!**********************************************************
! Parametres des nombres aleatoires
!**********************************************************
! SPRNG renvoie un nbre aleatoire dans [0,1[
! faire attention a conversion dp -> sp sur machine 64 bits
! qui peut donner 1.0 -> BUG
!**********************************************************

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

end module naleat
