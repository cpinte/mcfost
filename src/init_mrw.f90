module zeta_MRW
    use mcfost_env, only : dp
    implicit none

    integer, parameter :: n = 10000
    real(dp), dimension(n) :: zeta, y_MRW 

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
           end if
        end do
     end if

!        write(*,*) i, j, y_MRW(i), 2*zeta(i)
     end do

  zeta = zeta * 2._dp

 end subroutine initialize_cumulative_zeta

end module zeta_MRW
