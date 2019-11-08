MODULE accelerate
! Adapted from RH H. Uitenbroek
  !
  ! Implements Ng structure similar to RH for convergence
  ! acceleration
  
  use mcfost_env, only : dp
  use utils, only :  GaussSlv

  IMPLICIT NONE
  
  
  logical, parameter :: improve_sol = .true.


  TYPE Ng
   integer :: N, N2, Ndelay, Norder, Nperiod, i
   real(8), allocatable, dimension(:) :: b
   real(8), allocatable, dimension(:,:) :: previous, A
  END TYPE Ng

  CONTAINS
  
  

  SUBROUTINE initNg(N, Ndelay, Norder, Nperiod, Ngs)!solution, 
   integer, intent(in) :: N, Ndelay, Norder, Nperiod
   !real(kind=dp), dimension(N), intent(in) :: solution
   type (Ng), intent(inout) :: Ngs

   Ngs%N = N
   Ngs%Norder = Norder !number of previous iterations kept
   Ngs%Nperiod = Nperiod !number of time we wait before a new acceleration
   !Because if Norder is 1, Ndelay should be at least 3
   Ngs%Ndelay = MAX(Ndelay, Norder+2) !number of time we wait until the first acceleration

   !We test outside that Ngorder > 0, and if it is not, we never enter here
   if (Norder > 0) then
    allocate(Ngs%A(Norder,Norder))
    allocate(Ngs%b(Norder))
    allocate(Ngs%previous(Norder+2,N)) !need to store Norder + 2 iterations
    Ngs%i= 0
   end if


  RETURN
  END SUBROUTINE initNg
  
  FUNCTION Acceleration(Ngs, solution)
   logical :: acceleration
   type(Ng), intent(inout) :: Ngs
   real(kind=dp), intent(inout), dimension(Ngs%N) :: solution!
   integer :: i, j, k
   real(kind=dp) :: weight, delta, di, dj!, Aij(Ngs%Norder, Ngs%Norder)
   
    acceleration = .false.
    !! at the moment
    !We test outside that Ngorder > 0
    !and outside that we start after the Ndelay
    !Again, we test outside the Nperiod
    !!
   
    !init a at 0
    Ngs%i = Ngs%i + 1
    Ngs%previous(Ngs%i,:) = solution(:)
    if (Ngs%i < Ngs%Norder+2) return
    
    acceleration = .true.
    
    !write(*,*) Ngs%Norder, Ngs%i," accumulated solutions, extrapolating..."

    Ngs%A(:,:) = 0d0
    Ngs%b(:) = 0d0
    
    do k=1, Ngs%N
    	if (solution(k)==0d0) cycle !anyway the residual should be 0 if the cell is empty
    	weight = 1d0 / (dabs(solution(k))) !never be zero
    	!delta(i==1) = latest stored solution
    	delta = Ngs%previous(Ngs%i,k) - Ngs%previous(Ngs%i-1,k) !r_ip+m
    	do i=1, Ngs%Norder
    	    !at Norder; we are at i=1, the first stored solution
        	di = Ngs%previous(Ngs%i-i,k) - Ngs%previous(Ngs%i-i-1,k)
        	Ngs%b(i) = Ngs%b(i) + weight * delta * (di - delta)
        	do j=1,Ngs%Norder
         		dj = Ngs%previous(Ngs%i,k) - Ngs%previous(Ngs%i-j-1,k)
         		Ngs%A(i,j) = Ngs%A(i,j) + weight * (di-delta)*(dj-delta)
        	enddo
        
      	enddo
      
    	enddo !end loop over space

    !maybe a transpose here
!     Aij = transpose(Ngs%A)
    CALL  GaussSlv(Ngs%A, Ngs%b, Ngs%Norder)
    
    do i=1, Ngs%Norder
     solution(:) = solution(:) - Ngs%b(i) * (Ngs%previous(Ngs%i,:) - Ngs%previous(Ngs%i-i,:))
    enddo
    
!     write(*,*) maxval(solution), minval(solution)
! stop
    !reset counter for next block of iterations
    Ngs%i = 0
   
  RETURN
  END FUNCTION Acceleration
  
  SUBROUTINE freeNg(Ngs)
    type (Ng), intent(inout) :: Ngs
    if (Ngs%Norder > 0) then
     deallocate(Ngs%A)
     deallocate(Ngs%b)
     deallocate(Ngs%previous)
    end if
  RETURN
  END SUBROUTINE freeNg


END MODULE accelerate

