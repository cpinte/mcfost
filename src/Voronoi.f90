module Voronoi_grid

  use parametres
  use utils, only : bubble_sort

  implicit none

  integer, parameter :: max_wall_neighbours = 100000

  type Voronoi_cell
     real :: x, y, z, V
     integer :: id, first_neighbour, last_neighbour
  end type Voronoi_cell

  type Voronoi_wall
     character(len=10) :: type
     real :: x1, x2, x3, x4, x5, x6, x7
     integer :: n_neighbours
     integer, dimension(max_wall_neighbours) :: neighbour_list ! Warning hard coded
  end type Voronoi_wall


  type(Voronoi_cell), dimension(:), allocatable :: Voronoi
  type(Voronoi_wall), dimension(:), allocatable :: wall
  integer, dimension(:), allocatable :: neighbours_list

  integer :: n_walls

  contains

  subroutine read_Voronoi(n)

    integer, intent(in) :: n

    integer :: i, j, k, ios, n_neighbours, n_neighbours_tot, ifirst


    n_walls = 6
    write(*,*) "Reading ", n_walls, "walls"


    call init_Voronoi_walls()

    allocate(Voronoi(n))

    write(*,*) "Reading Voronoi :", n, "cells"


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

    rewind(1)
    write(*,*) "neighbours list size =", n_neighbours_tot
    write(*,*)  "Voronoi volume =", sum(Voronoi%V)
    write(*,*) "Trying to allocate", 4*n_neighbours_tot/ 1024.**2, "MB for neighbours list"
    allocate(neighbours_list(n_neighbours_tot))



    do i=1, n
       read(1,*) Voronoi(i)%id, Voronoi(i)%x, Voronoi(i)%y, Voronoi(i)%z, Voronoi(i)%V, n_neighbours, neighbours_list(Voronoi(i)%first_neighbour:Voronoi(i)%last_neighbour)

       ! todo : find the cells touching the walls
       do k=1, n_neighbours
          j = neighbours_list(Voronoi(i)%first_neighbour + k-1)
          if (j < 0) then ! wall
             j = -j ! wall index
             wall(j)%n_neighbours = wall(j)%n_neighbours+1
             if (wall(j)%n_neighbours > max_wall_neighbours) then
                write(*,*) "ERROR : Voronoi wall", j, "max number of neighbours reached"
             endif
             wall(j)%neighbour_list(wall(j)%n_neighbours) = i
          endif ! wall
       enddo ! k
    enddo ! i

    close(unit=1)

    do k=1, n_walls
       write(*,*) "wall", k, wall(k)%n_neighbours, "voisins"
       !if (k==1) then
       !   do i=1,  wall(k)%n_neighbours
       !      write(*,*) i, wall(k)%neighbour_list(i)
       !   enddo
       !endif
    enddo


    return

  end subroutine read_Voronoi

  !----------------------------------------

  subroutine cross_Voronoi_cell(icell, previous_cell, x,y,z, u,v,w, next_cell, s)

    integer, intent(in) :: icell, previous_cell
    real, intent(in) :: x,y,z, u,v,w

    real, intent(out) ::  s
    integer, intent(out) :: next_cell

    real :: s_tmp, den
    integer :: i, in

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real, dimension(3) :: n, p, r, k

    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    s = 1e30 !huge_real
    next_cell = 0
    nb_loop : do i=Voronoi(icell)%first_neighbour, Voronoi(icell)%first_neighbour
       in = neighbours_list(i) ! id du voisin

       if (in==previous_cell) cycle nb_loop

       if (in > 0) then ! cellule
          ! unnormalized vector to plane
          n(1) = Voronoi(in)%x - Voronoi(icell)%x
          n(2) = Voronoi(in)%y - Voronoi(icell)%y
          n(3) = Voronoi(in)%z - Voronoi(icell)%z

          ! test direction
          den = dot_product(n, k)
          if (den <= 0) cycle nb_loop

          ! point on the plane
          p(1) = 0.5 * (Voronoi(in)%x + Voronoi(icell)%x)
          p(2) = 0.5 * (Voronoi(in)%y + Voronoi(icell)%y)
          p(3) = 0.5 * (Voronoi(in)%z + Voronoi(icell)%z)

          s_tmp = dot_product(n, p-r) / den
       else ! i < 0 ; le voisin est un wall
          s_tmp = distance_to_wall(x,y,z, u,v,w, -in) ;
       endif

       if (s_tmp < s) then
          s = s_tmp
          next_cell = in
       endif
    enddo nb_loop ! i

    return

  end subroutine cross_Voronoi_cell

  !----------------------------------------

  subroutine init_Voronoi_walls()

    integer :: iwall

    real, parameter :: xmin = 0 , xmax = 1, ymin = 0, ymax = 1, zmin = 0, zmax = 1

    allocate(wall(n_walls))

    ! initialise plane walls
    do iwall=1, n_walls
       wall(iwall)%type = "plane"
       wall(iwall)%n_neighbours = 0
    enddo

    ! test pour localiser les murs par defaut
    !if (j==6) then
    !   write(*,*) Voronoi(i)%x, Voronoi(i)%y, Voronoi(i)%z
    !endif

    ! j=1 ---> x = xmin ; n = (-1,0,0) : normal towards outside
    ! j=2 ---> x = xmax
    ! j=3 ---> y = ymin
    ! j=4 ---> y = ymax
    ! j=5 ---> z = zmin
    ! j=6 ---> z = zmax

    wall(1)%x1 = -1 ; wall(1)%x2 = 0  ; wall(1)%x3 = 0  ; wall(1)%x4 = -xmin
    wall(2)%x1 =  1 ; wall(2)%x2 = 0  ; wall(2)%x3 = 0  ; wall(2)%x4 =  xmax
    wall(3)%x1 =  0 ; wall(3)%x2 = -1 ; wall(3)%x3 = 0  ; wall(3)%x4 = -ymin
    wall(4)%x1 =  0 ; wall(4)%x2 = 1  ; wall(4)%x3 = 0  ; wall(4)%x4 =  ymax
    wall(5)%x1 =  0 ; wall(5)%x2 = 0  ; wall(5)%x3 = -1 ; wall(5)%x4 = -zmin
    wall(6)%x1 =  0 ; wall(6)%x2 = 0  ; wall(6)%x3 = 1  ; wall(6)%x4 =  zmax


    return


  end subroutine init_Voronoi_walls

  !----------------------------------------


  real function distance_to_wall(x,y,z, u,v,w, iwall)
    ! Mur plan pour le moment : meme algorithme que cross Voronoi cell
    ! renvoie une valeur :
    ! - negative si on rentre dans le modele (car la normale du mur est vers l'exterieur)
    ! - positive is on sort du modele

    real, intent(in) :: x,y,z, u,v,w
    integer, intent(in) :: iwall

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real, dimension(3) :: n, p, r, k

    real :: den

    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    n(1) = wall(iwall)%x1 ; n(2) = wall(iwall)%x2 ;  n(3) = wall(iwall)%x3 ;
    p = wall(iwall)%x4 * n  ! todo : verifier le signe

    den = dot_product(n, k)

    if (abs(den) > 0) then
       distance_to_wall = dot_product(n, p-r) / den
    endif

    return

  end function distance_to_wall

  !----------------------------------------

  subroutine move_to_Voronoi_grid(x,y,z, u,v,w, s,icell)

    real, intent(in) :: x,y,z,u,v,w
    real, intent(out) :: s
    integer, intent(out) :: icell

    logical, dimension(n_walls) :: intersect
    real, dimension(n_walls) :: s_walls
    integer, dimension(n_walls) :: order

    real :: l
    integer :: iwall

    ! Find out which plane we are crossing first
    ! and test boundaries of the plane
    s_walls(:) = huge(1.0)
    intersect(:) = .false.
    do iwall=1, n_walls
       l = -distance_to_wall(x,y,z, u,v,w, iwall) ! signe - car on rentre dans le volume
       if (l >= 0) then
          intersect(iwall) = .true.
          s_walls(iwall) = l
       endif
    enddo

    order = bubble_sort(real(s_walls,kind=db))

    ! Move to the plane

    ! Find out the closest cell

    ! Move to the cell (if wall is approximate)

    return

  end subroutine move_to_Voronoi_grid

  !----------------------------------------

!  subroutine length_Voronoi(id,lambda,Stokes,id,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
!    !Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
!    ! C. Pinte
!
!    integer, intent(in) :: id,lambda
!    integer, intent(inout) :: id
!    real(kind=db), dimension(4), intent(in) :: Stokes
!    logical, intent(in) :: flag_star, flag_direct_star
!    real(kind=db), intent(inout) :: u,v,w
!    real, intent(in) :: extrin
!    real(kind=db), intent(inout) :: xio,yio,zio
!    real, intent(out) :: ltot
!    logical, intent(out) :: flag_sortie
!
!    return
!
!  end subroutine length_Voronoi


end module Voronoi_grid

!----------------------------------------

program Voronoi

  use Voronoi_grid

  implicit none

  call read_Voronoi(999997)

end program Voronoi
