module Voronoi_grid

  implicit none

  type Voronoi_cell
     real :: x, y, z, V
     integer :: id, first_neighbour, last_neighbour
  end type Voronoi_cell

  type(Voronoi_cell), dimension(:), allocatable :: Voronoi
  integer, dimension(:), allocatable :: neighbours_list


  contains

  subroutine read_Voronoi(n)

    integer, intent(in) :: n

    integer :: i, ios, n_neighbours, n_neighbours_tot, ifirst


    allocate(Voronoi(n))

    n_neighbours_tot = 0
    open(unit=1, file="Voronoi.txt", status='old', iostat=ios)
    do i=1, n
       read(1,*) Voronoi(i)%id, Voronoi(i)%x, Voronoi(i)%y, Voronoi(i)%z, Voronoi(i)%V, n_neighbours
       if (i>1) then
          Voronoi(i)%first_neighbour = Voronoi(i-1)%last_neighbour + 1
          Voronoi(i)%last_neighbour  = Voronoi(i-1)%last_neighbour + n_neighbours
       else
          Voronoi(i)%first_neighbour = 1
          Voronoi(i)%last_neighbour =  n_neighbours
       endif
       n_neighbours_tot = n_neighbours_tot + n_neighbours
    enddo

    write(*,*)  n_neighbours_tot, "voisins en tout"
    write(*,*)  sum(Voronoi%V)
    rewind(1)

    allocate(neighbours_list(n_neighbours_tot))

    do i=1, n
       read(1,*) Voronoi(i)%id, Voronoi(i)%x, Voronoi(i)%y, Voronoi(i)%z, Voronoi(i)%V, n_neighbours, neighbours_list(Voronoi(i)%first_neighbour:Voronoi(i)%last_neighbour)
    enddo


    close(unit=1)

    return

  end subroutine read_Voronoi

  !----------------------------------------

  subroutine cross_Voronoi_cell(id, x,y,z, u,v,w, next_cell, s)

    integer, intent(in) :: id
    real, intent(in) :: x,y,z, u,v,w

    real, intent(out) ::  s
    integer, intent(out) :: next_cell

    real :: s_tmp, den
    integer :: i

    real :: n, p ! vectors

    s = 1e30 !huge_real
    next_cell = 0
    nb_loop : do i=Voronoi(id)%first_neighbour, Voronoi(id)%first_neighbour

       ! unnormalized vector to plane
       n = Voronoi(id)%x - Voronoi(i)%x

       ! point on the plane
       p = 0.5 * (Voronoi(id)%x + Voronoi(i)%x)

       ! test direction
       den = u * p
       if (den <= 0) cycle nb_loop ! dot product

       s_tmp = n * (p -x) ! dot product
       if (s_tmp < s) then
          s = s_tmp
          next_cell = i
       endif

    enddo nb_loop

    return

  end subroutine cross_Voronoi_cell


end module Voronoi_grid

!----------------------------------------

program Voronoi

  use Voronoi_grid

  implicit none

  call read_Voronoi(999997)

end program Voronoi
