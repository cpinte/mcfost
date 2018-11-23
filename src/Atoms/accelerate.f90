MODULE accelerate
! Adapted from RH H. Uitenbroek
  !
  ! Implements Ng structure similar to RH for convergence
  ! acceleration
  !
  
  !use some_modules, only : linear_eq_solver ??

  IMPLICIT NONE
  
  
 logical, parameter :: improve_sol = .true.


  TYPE Ng
   integer :: N, Ndelay, Norder, Nperiod, count
   real(8), allocatable, dimension(:) :: b, theStorage
   real(8), allocatable, dimension(:,:) :: previous, A
  END TYPE Ng

  CONTAINS
  
  

  SUBROUTINE initNg(N, Ndelay, Norder, Nperiod, solution, Ngs)
   integer, intent(in) :: N, Ndelay, Norder, Nperiod
   integer :: k
   real(8), allocatable, dimension(:), intent(in) :: solution
   type (Ng), intent(inout) :: Ngs

   Ngs%N = N
   Ngs%Norder = Norder
   Ngs%Nperiod = Nperiod
   Ngs%Ndelay = MAX(Ndelay, Norder+2)

   if (Norder.lt.0) then
    allocate(Ngs%A(Norder,Norder))
    allocate(Ngs%b(Norder))
   end if

   allocate(Ngs%previous(Norder+2,N))
   Ngs%theStorage = Ngs%previous(1,:)

   do k=1,N
    Ngs%previous(1,k)=solution(k)
   end do
   Ngs%count = 1

  RETURN
  END SUBROUTINE initNg

  SUBROUTINE NgAcceleration(Ngs, solution, accel)
   !Ngs already initialised
   integer :: i, j, k, ip, ipp, i0
   logical, intent(inout) :: accel
   real(8), allocatable, dimension(:) :: weight
   real(8), allocatable, dimension(:,:) ::  delta
   type (Ng), intent(inout) :: Ngs
   real(8), dimension(Ngs%N), intent(inout) :: solution
   !logical :: improve_sol = .true.

   i = MOD(Ngs%count, (Ngs%Norder+2))
   do k=1,Ngs%N
    Ngs%previous(i,k) = solution(k)
   end do
   Ngs%count = Ngs%count + 1 !increament for next iteration

   if ((Ngs%Norder.gt.0).and.(Ngs%count.ge.Ngs%Ndelay) &
       .and.(MOD(Ngs%count - Ngs%Ndelay, Ngs%Nperiod).gt.0) ) then

      allocate(Delta(Ngs%Norder+1,Ngs%N))
      allocate(weight(Ngs%N))

      do i=1,Ngs%Norder
       ip = MOD(Ngs%count -1 -i, Ngs%Norder+2)
       ipp = MOD(Ngs%count - 2 -i, Ngs%Norder+2)
       do k=1,Ngs%N
         Delta(i,k) = Ngs%previous(ip,k) - Ngs%previous(ipp,k)
       end do
      end do

      do k=1,Ngs%N
       weight(k) = 1/dabs(solution(k))
      end do
      do i=1,Ngs%Norder
       Ngs%b(i) = 0.
       do j=1,Ngs%Norder
        Ngs%A(i,j) = 0.
       end do
      end do


     do j=1,Ngs%Norder
       do k=1,Ngs%N
        Ngs%b(j) = Ngs%b(j) + weight(k) * Delta(1,k)* &
                      (Delta(1,k) - Delta(j+1,k))
       end do
       do i=1,Ngs%Norder
        do k=1,Ngs%N
         Ngs%A(i,j) = Ngs%A(i,j) + weight(k)*&
                     (Delta(j+1,k) - Delta(1,k)) * &
                     (Delta(i+1,k) - Delta(1,k))
        end do
       end do
     end do
   
     !CALL some_linear_solver(Ngs%Norder,Ngs%A,Ngs%b, improve_sol)

     i0 = MOD(Ngs%count -1, Ngs%Norder+2)
     do i=1,Ngs%Norder
      ip = MOD(Ngs%count-2-i, Ngs%Norder+2)
      do k=1,Ngs%N
       solution(k) = solution(k) + Ngs%b(i) * &
             (Ngs%previous(ip,k) - Ngs%previous(i0,k))
      end do
     end do

     do k=1,Ngs%N
      Ngs%previous(i0,k) = solution(k)
     end do

     accel = .true.

     deallocate(Delta)
     deallocate(weight)
   else
     accel = .false.
   end if

  RETURN
  END SUBROUTINE NgAcceleration

  SUBROUTINE freeNg(Ngs)
    type (Ng), intent(inout) :: Ngs
    if (Ngs%Norder.gt.0) then
    deallocate(Ngs%A)
    deallocate(Ngs%b)
    end if
  RETURN
  END SUBROUTINE freeNg

  SUBROUTINE MaxChange(Ngs, text, verbose, dM)
    type(Ng), intent(in) :: Ngs
    logical, intent(in) :: verbose
    character(len=*), intent(in) :: text
    real, intent(out) :: dM
    integer :: k
    real(8), allocatable, dimension(:) :: old, new
    dM = 0.

    allocate(old(MOD(Ngs%count-2,Ngs%Norder+2)))
    allocate(new(MOD(Ngs%count-1,Ngs%Norder+2)))

    do k=1,Ngs%N
     if (new(k).gt.0.) then
      dM = MAX(dM, dabs((new(k)-old(k))/new(k)))
     end if
    end do

    if (verbose) then
     write(*,*) text," delta = ", dM
    end if

  RETURN
  END SUBROUTINE MaxChange


END MODULE accelerate

