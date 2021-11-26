module sort

  use mcfost_env, only : sp, dp
  use utils, only : indgen
  use messages, only : error

  implicit none

contains

function index_bubble_sort(data_in)
  ! Implementation of Bubble sort
  ! Warning : this is N^2, only for small arrays
  ! Same behaviour as yorick to allow ordering of mutiple arrays
  ! Return the order of data, the sorted array would be data_in(order)
  !
  ! C. Pinte
  ! 02/05/11

  real(kind=dp), dimension(:), intent(in) :: data_in
  real(kind=dp), dimension(:), allocatable :: data
  integer, dimension(:), allocatable :: index_bubble_sort ! indices

  integer :: i, pass, n, tmp_i
  logical :: sorted
  real(kind=dp) ::  temp

  n = size(data_in)
  allocate(index_bubble_sort(n),data(n))
  data = data_in

  index_bubble_sort = indgen(n)

  pass = 1
  sorted = .false.

  do while(.not.sorted)
     sorted = .true.

     do i = 1,n-pass
        if(data(i) > data(i+1)) then
           temp = data(i)
           data(i) = data(i+1)
           data(i+1) = temp
           sorted = .false.

           ! same operation on indices
           tmp_i = index_bubble_sort(i)
           index_bubble_sort(i) = index_bubble_sort(i+1)
           index_bubble_sort(i+1) = tmp_i
        endif
     end do ! i
     pass = pass +1
  end do ! while

  return

end function index_bubble_sort

!************************************************************

subroutine quicksort(a)
  ! Non-recursive version of the quicksort algorithm using a stack.
  ! By always pushing the larger "half" to the stack, and moving directly to calculate
  ! the smaller "half", it ensures that the stack needs no more than log_2(N) entries

  real(dp), dimension(:), intent(inout):: a
  real(dp) :: temp,pivot
  integer :: i,j,left_s,right_s,left,right,stack_ptr
  integer, dimension(2,64) :: stack

  integer, parameter :: Qsort_limit = 50

  left=1
  right=size(a)
  stack_ptr=1

  do
     if (right-left < Qsort_limit) then ! use insertion sort on small arrays
        do i=left+1,right
           temp=a(i)
           do j=i-1,left,-1
              if (a(j) <= temp) exit
              a(j+1)=a(j)
           enddo
           a(j+1)=temp
        enddo
        ! pop from stack
        if (stack_ptr == 1) return
        stack_ptr=stack_ptr-1
        left=stack(1,stack_ptr)
        right=stack(2,stack_ptr)
     else
        ! find median of three pivots and place sentinels at first and last elements
        temp=a((left+right)/2)
        a((left+right)/2)=a(left+1)
        if (temp > a(right)) then
           a(left+1)=a(right)
           a(right)=temp
        else
           a(left+1)=temp
        endif
        if (a(left) > a(right)) then
           temp=a(left)
           a(left)=a(right)
           a(right)=temp
        endif
        if (a(left) > a(left+1)) then
           temp=a(left)
           a(left)=a(left+1)
           a(left+1)=temp
        endif
        pivot=a(left+1)
        left_s=left+2
        right_s=right-1
        do
           do while(a(left_s) < pivot)
              left_s=left_s+1
           enddo
           do while(a(right_s) > pivot)
              right_s=right_s-1
           enddo
           if (left_s >= right_s) exit
           temp=a(left_s)
           a(left_s)=a(right_s)
           a(right_s)=temp
           left_s=left_s+1
           right_s=right_s-1
        enddo
        if (left_s == right_s) left_s=left_s+1
        if (left_s < (left+right)/2) then
           stack(1,stack_ptr)=left_s
           stack(2,stack_ptr)=right
           stack_ptr=stack_ptr+1
           right=left_s-1
        else
           stack(1,stack_ptr)=left
           stack(2,stack_ptr)=left_s-1
           stack_ptr=stack_ptr+1
           left=left_s
        endif
     endif
  enddo

  return

end subroutine quicksort

!************************************************************

function index_quicksort(a) result(idx)
  ! Use the same non-recursive quicksort as above, but return an array with the order of the indices
  ! Does not alter the orginal array

  real(kind=dp), dimension(:), intent(in)  :: a
  integer, dimension(:), allocatable :: idx

  real(dp) :: temp,pivot
  integer :: i,j,left_s,right_s,left,right,stack_ptr,i_tmp,k
  integer, dimension(2,64) :: stack

  integer, parameter :: Qsort_limit = 50

  left=1
  right=size(a)
  stack_ptr=1

  allocate(idx(right))
  idx = indgen(right)

  do
     if (right-left < Qsort_limit) then ! use insertion sort on small arrays
        do i=left+1,right
           i_tmp = idx(i)
           temp=a(i_tmp)
           do j=i-1,left,-1
              if (a(idx(j)) <= temp) exit
              idx(j+1)=idx(j)
           enddo
           idx(j+1)=i_tmp
        enddo
        ! pop from stack
        if (stack_ptr == 1) return
        stack_ptr=stack_ptr-1
        left=stack(1,stack_ptr)
        right=stack(2,stack_ptr)
     else
        ! find median of three pivots and place sentinels at first and last elements
        k=(left+right)/2
        i_tmp=idx(k)
        idx(k)=idx(left+1)
        if (a(i_tmp) > a(idx(right))) then
           idx(left+1)=idx(right)
           idx(right)=i_tmp
        else
           idx(left+1)=i_tmp
        endif
        if (a(idx(left)) > a(idx(right))) then
           i_tmp=idx(left)
           idx(left)=idx(right)
           idx(right)=i_tmp
        endif
        if (a(idx(left)) > a(idx(left+1))) then
           i_tmp=idx(left)
           idx(left)=idx(left+1)
           idx(left+1)=i_tmp
        endif
        pivot=a(idx(left+1))
        left_s=left+2
        right_s=right-1
        do
           do while(a(idx(left_s)) < pivot)
              left_s=left_s+1
           enddo
           do while(a(idx(right_s)) > pivot)
              right_s=right_s-1
           enddo
           if (left_s >= right_s) exit
           i_tmp=idx(left_s)
           idx(left_s)=idx(right_s)
           idx(right_s)=i_tmp
           left_s=left_s+1
           right_s=right_s-1
        enddo
        if (left_s == right_s) left_s=left_s+1
        if (left_s < (left+right)/2) then
           stack(1,stack_ptr)=left_s
           stack(2,stack_ptr)=right
           stack_ptr=stack_ptr+1
           right=left_s-1
        else
           stack(1,stack_ptr)=left
           stack(2,stack_ptr)=left_s-1
           stack_ptr=stack_ptr+1
           left=left_s
        endif
     endif
  enddo

  return

end function index_quicksort

!************************************************************

subroutine test_quicksort()

  integer, parameter :: n=100000

  real(kind=dp), dimension(n) :: A, B
  integer, dimension(n) :: order
  integer :: seed_size
  integer, dimension(:), allocatable :: seed
  integer :: i
  real :: r

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed(:) = 42
  call random_seed(put=seed)

  do i=1,n
     call random_number(r)
     A(i) = r
  enddo
  B = A

  write(*,*) "testing quicksort ..."
  call quicksort(A)
  do i=2,n
     if (a(i) < a(i-1)) call error("quicksort")
  enddo
  write(*,*) "ok"

  write(*,*) "testing index_quicksort ..."
  order = index_quicksort(B)
  do i=2,n
     if (b(order(i)) < b(order(i-1))) call error("index_quicksort")
  enddo
  write(*,*) "ok"
  stop

  return

end subroutine test_quicksort

!************************************************************

function find_kth_smallest(k,array)
  ! Returns the kth smallest value in the array.
  ! The input array will be rearranged to have this value in location array(k),
  ! with all smaller elements moved to arr(1:k-1) (in arbitrary order) and
  ! all larger elements in arr(k+1:) (also in arbitrary order).

  integer, intent(in) :: k
  real(sp), dimension(:), intent(inout) :: array
  real(sp) :: find_kth_smallest

  integer :: i,r,j,l
  real(sp) :: a

  l=1
  r=size(array)
  do
     if (r-l <= 1) then ! Active partition contains 1 or 2 elements.
        if (r-l == 1) then
           if (array(l)>array(r)) call swap(array(l),array(r))  ! Active partition contains 2 elements.
        endif
        find_kth_smallest=array(k)
        return
     else
        ! Choose median of left, center, and right elements as partitioning element a.
        ! Also rearrange so that array(l) <= array(l+1) <= array(r).
        i=(l+r)/2
        call swap(array(i),array(l+1))
        if (array(l)>array(r))   call swap(array(l),array(r))
        if (array(l+1)>array(r)) call swap(array(l+1),array(r))
        if (array(l)>array(l+1)) call swap(array(l),array(l+1))
        i=l+1 ! Initialize pointers for partitioning.
        j=r
        a=array(l+1) ! Partitioning element.
        do
           do ! Scan up to find element > a.
              i=i+1
              if (array(i) >= a) exit
           enddo
           do ! Scan down to find element < a.
              j=j-1
              if (array(j) <= a) exit
           enddo
           if (j < i) exit ! Pointers crossed. Exit with partitioning complete.
           call swap(array(i),array(j)) ! Exchange elements.
        enddo
        array(l+1)=array(j) ! Insert partitioning element.
        array(j)=a
        if (j >= k) r=j-1 ! Keep active the partition that contains the kth element.
        if (j <= k) l=i
     endif
  enddo

  return

end function find_kth_smallest

!************************************************************

function find_kth_smallest_inplace(k,array)
  ! Returns the kth smallest value in the array, without altering the input array.
  ! In Fortran 90's assumed memory-rich environment, we just call select in scratch space.
  ! C.P. : this is an issue for very large arrays: using allocatable array instead

  integer, intent(in) :: k
  real(sp), dimension(:), intent(in) :: array
  real(sp) :: find_kth_smallest_inplace

  real(sp), dimension(:), allocatable :: tmp_array

  allocate(tmp_array(size(array)))

  tmp_array=array
  find_kth_smallest_inplace=find_kth_smallest(k,tmp_array)
  deallocate(tmp_array)

  return

end function find_kth_smallest_inplace

!************************************************************

subroutine swap(a,b)
  real(sp), intent(inout) :: a,b
  real(sp) :: tmp

  tmp=a ; a=b ; b=tmp
  return

end subroutine swap

!***************************************************

subroutine Knuth_shuffle(a)

  integer, intent(inout) :: a(:)

  integer :: seed_size
  integer, dimension(:), allocatable :: seed
  integer :: i, randpos, temp
  real :: r

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  seed(:) = 42
  call random_seed(put=seed)

  do i = size(a), 2, -1
     call random_number(r)
     randpos = int(r * i) + 1
     temp = a(randpos)
     a(randpos) = a(i)
     a(i) = temp
  enddo

  return

end subroutine Knuth_Shuffle

end module sort
