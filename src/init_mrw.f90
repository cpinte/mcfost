program initialize_cumulative_zeta
  ! Compute equation (7) of Min et al (2009):
  ! P(t) = 2 * Sum_{n=1}^{infinity} (-1)^(n+1) * y^(n^2)
  !
  ! The function requires high n values close to y=1,
  ! so we pre-compute it here and we will then interpolate

  implicit none

  integer, parameter :: sp = selected_real_kind(p=6,r=37)
  integer, parameter :: dp = selected_real_kind(p=13,r=200)

  integer, parameter :: n = 10000

  real(kind=8), dimension(n) :: x, y
  integer :: i,j
  real(kind=8) :: term

  y = 0._dp

  do i=1,n
     x(i) = real(i-1, dp)/real(n-1, dp)

     if(i==n) then
        y(i) = 0.5_dp
     else
        j = 0
        do
           j = j + 1
           term = x(i)**(j**2)
           if(term == 0._dp) exit
           if(mod(j, 2)==0) then
              y(i) = y(i) - term
           else
              y(i) = y(i) + term
           end if
        end do
     end if

        write(*,*) i, j, x(i), 2*y(i)
     end do

  y = y * 2._dp

end program initialize_cumulative_zeta
