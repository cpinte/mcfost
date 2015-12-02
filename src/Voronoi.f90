module Voronoi_grid

  use constantes
  use parametres
  use utils, only : bubble_sort, appel_syst
  use naleat, only : seed, stream, gtype

  implicit none

  integer, parameter :: max_wall_neighbours = 100000
  real(kind=db), parameter :: prec = 1.e-6_db

  type Voronoi_cell
     real :: x, y, z, V
     integer :: id, first_neighbour, last_neighbour
  end type Voronoi_cell

  type Voronoi_wall
     character(len=10) :: type
     real :: x1, x2, x3, x4, x5, x6, x7
     integer :: n_neighbours
     integer, dimension(max_wall_neighbours) :: neighbour_list ! Warning hard coded

     ! Plane wall :
     ! x1, x2, x3 : normal
     ! x4 : displacement along the normal

  end type Voronoi_wall


  type(Voronoi_cell), dimension(:), allocatable :: Voronoi
  type(Voronoi_wall), dimension(:), allocatable :: wall
  integer, dimension(:), allocatable :: neighbours_list

  integer :: n_walls

  contains

  subroutine read_Voronoi(n)

    integer, intent(in) :: n

    integer :: i, j, k, ios, n_neighbours, n_neighbours_tot, ifirst

    !integer, dimension(:), allocatable :: id_list

    ! For testing purposes
    real :: x, y, z, u, v, w, s, norme
    integer :: icell

    integer :: id, lambda
    real :: tau, lvol

    real(kind=db), dimension(4) :: Stokes
    logical :: flag_star, flag_direct_star, flag_sortie

    ! Fin testing

    id = 1

    n_walls = 6
    write(*,*) "Finding ", n_walls, "walls"


    call init_Voronoi_walls()
    allocate(Voronoi(n))

    write(*,*) "Finding:", n, " Voronoi cells"

    n_neighbours_tot = 0
    open(unit=1, file="Voronoi.txt", status='old', iostat=ios)
    do i=1, n
       read(1,*) Voronoi(i)%id, Voronoi(i)%x, Voronoi(i)%y, Voronoi(i)%z, Voronoi(i)%V, n_neighbours
       Voronoi(i)%id = i ! id a un PB car Voronoi fait sauter des points
       !write(*,*) "Voronoi id = ", Voronoi(i)%id

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
       ! id a un PB car Voronoi fait sauter des points quand bcp de ponts
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
       write(*,*) "wall", k, wall(k)%n_neighbours, "neighbours"
       !if (k==1) then
       !   do i=1,  wall(k)%n_neighbours
       !      write(*,*) i, wall(k)%neighbour_list(i)
       !   enddo
       !endif
    enddo


    write(*,*) "Testing radiative transfer routines on Voronoi grid"
    x = -200.0 ; y = -200.0 ; z = -200.0 ;
    !u = 1.2 ; v = 1.0 ; w = 1.1 ;
    !u = 1.45 ; v = 1.2 ; w = 1.0 ;

    u = 1.0 ; v = 1.0 ; w = 1.0 ;
    norme = sqrt(u*u + v*v + w*w)
    u = u / norme ; v = v / norme ;  w = w / norme ;

    call move_to_Voronoi_grid(x,y,z, u,v,w, s, icell)

    x = x + s*u ; y = y + s*v ; z = z+s*w ! a mettre dans move ??

    write(*,*) "Packet is entering volume in cell", icell

    id = 1 ; lambda = 1 !TODO
    Stokes(1) = 1.0 ; Stokes(2:4) = 0.0

    flag_star = .true.
    flag_direct_star = .true.

    tau = 500 ;

    write(*,*) "Testing length_Voronoi with tau=", tau

    if (icell > 0) then
       call length_Voronoi(id,lambda,Stokes, icell,x,y,z, u,v,w, flag_star,flag_direct_star, tau, lvol,flag_sortie)

       write(*,*) "lvol=", lvol
       write(*,*) "Did I exit ?", flag_sortie
       write(*,*) "Last cell", icell
       write(*,*) "Final position", x, y, z
    else
       write(*,*) "Packet did not reach model volume"
    endif

    write(*,*) "TEST DONE"

    ! OK jusqu'ici
    !call test_emission()

    return

  end subroutine read_Voronoi

!----------------------------------------

  subroutine Voronoi_tesselation(n_cells, x,y,z)

    integer, intent(in) :: n_cells
    real(kind=db), dimension(n_cells), intent(in) :: x, y, z

    character(len=512) :: cmd
    integer :: i, syst_status

    open(unit=1, file="particles.txt", status="replace")
    do i=1, n_cells
       write(unit=1,fmt="(i5,f15.6,f15.6,f15.6)") i, real(x(i)), real(y(i)), real(z(i))
    enddo
    close(unit=1)

    ! Run voro++ command line for now
    ! while I fix the c++/fortran interface and make all the tests
    write(*,*) "Performing Voronoi tesselation on ", n_cells, "SPH particles"
    cmd = "~/codes/voro++-0.4.6/src/voro++  -v -o -g -c '%i %q %v %s %n' -150 150 -150 150 -150 150 particles.txt ; mv particles.txt.vol Voronoi.txt"
    call appel_syst(cmd,syst_status)
    write(*,*) "Voronoi Tesselation done"

    write(*,*) "TMP : filtering out 10 cells for safety, will do it better later"

    call read_Voronoi(n_cells-10)

    write(*,*) "Tesselation finished"
    return

  end subroutine Voronoi_tesselation

!----------------------------------------

  subroutine test_emission()

#include "sprng_f.h"

    integer, parameter :: n_sample = 10000000

    real(kind=db) :: rand, rand2, rand3
    real :: x, y, z
    integer :: icell, np_proc, i, id, k


    nb_proc = 1 ; id = 1
    allocate(stream(nb_proc))

    stream = 0.0
    do i=1, nb_proc
       !write(*,*) gtype, i-1,nb_proc,seed,SPRNG_DEFAULT
       !init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
       stream(i) = init_sprng(gtype, i-1,nb_proc,seed,SPRNG_DEFAULT)
    enddo


    ! Testing pos
    icell = 1
    write(*,*) "testing emission position", Voronoi(icell)%x, Voronoi(icell)%y, Voronoi(icell)%z

    do k=1, n_sample
       rand  = sprng(stream(id))
       rand2 = sprng(stream(id))
       rand3 = sprng(stream(id))

       call pos_em_cellule_Voronoi(icell,rand,rand2,rand3, x,y,z)
       write(*,*) x,y,z
    enddo

  end subroutine test_emission

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

    nb_loop : do i=Voronoi(icell)%first_neighbour, Voronoi(icell)%last_neighbour
       in = neighbours_list(i) ! id du voisin
       if (in==previous_cell) cycle nb_loop

       if (in > 0) then ! cellule
          ! unnormalized vector to plane
          n(1) = Voronoi(in)%x - Voronoi(icell)%x
          n(2) = Voronoi(in)%y - Voronoi(icell)%y
          n(3) = Voronoi(in)%z - Voronoi(icell)%z

          ! test direction
          den = dot_product(n, k)
          if (den <= 0) cycle nb_loop ! car s_tmp sera < 0

          ! point on the plane
          p(1) = 0.5 * (Voronoi(in)%x + Voronoi(icell)%x)
          p(2) = 0.5 * (Voronoi(in)%y + Voronoi(icell)%y)
          p(3) = 0.5 * (Voronoi(in)%z + Voronoi(icell)%z)

          s_tmp = dot_product(n, p-r) / den
       else ! i < 0 ; le voisin est un wall
          s_tmp = distance_to_wall(x,y,z, u,v,w, -in) ;

          ! si c'est le wall d'entree : peut-etre a faire sauter en sauvegardabt le wall d'entree
          if (s_tmp < 0.) s_tmp = huge(1.0)
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

    real, parameter :: xmin = -150, xmax = 150
    real, parameter :: ymin = xmin, ymax = xmax
    real, parameter :: zmin = xmin, zmax = xmax

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

    wall(1)%x1 = -1 ; wall(1)%x2 = 0  ; wall(1)%x3 = 0  ; wall(1)%x4 = abs(xmin)
    wall(2)%x1 =  1 ; wall(2)%x2 = 0  ; wall(2)%x3 = 0  ; wall(2)%x4 = abs(xmax)
    wall(3)%x1 =  0 ; wall(3)%x2 = -1 ; wall(3)%x3 = 0  ; wall(3)%x4 = abs(ymin)
    wall(4)%x1 =  0 ; wall(4)%x2 = 1  ; wall(4)%x3 = 0  ; wall(4)%x4 = abs(ymax)
    wall(5)%x1 =  0 ; wall(5)%x2 = 0  ; wall(5)%x3 = -1 ; wall(5)%x4 = abs(zmin)
    wall(6)%x1 =  0 ; wall(6)%x2 = 0  ; wall(6)%x3 = 1  ; wall(6)%x4 = abs(zmax)

    return

  end subroutine init_Voronoi_walls

  !----------------------------------------


  real function distance_to_wall(x,y,z, u,v,w, iwall)
    ! Mur plan pour le moment : meme algorithme que cross Voronoi cell

    real, intent(in) :: x,y,z, u,v,w
    integer, intent(in) :: iwall

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real, dimension(3) :: n, p, r, k

    real :: den

    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    n(1) = wall(iwall)%x1 ; n(2) = wall(iwall)%x2 ;  n(3) = wall(iwall)%x3 ;

    p = wall(iwall)%x4 * n

    den = dot_product(n, k) ! le signe depend du sens de propagation par rapport a la normale

    if (abs(den) > 0) then
       distance_to_wall = dot_product(n, p-r) / den
    else
       distance_to_wall = huge(1.0)
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

    real :: l, x_test, y_test, z_test
    integer :: i, iwall


    ! Find out which plane we are crossing first
    ! and test boundaries of the plane
    s_walls(:) = huge(1.0)
    intersect(:) = .false.
    do iwall=1, n_walls
       l = distance_to_wall(x,y,z, u,v,w, iwall) ! signe - car on rentre dans le volume

       if (l >= 0) then
          intersect(iwall) = .true.
          s_walls(iwall) = l * (1.0_db + prec)
       else
          s_walls(iwall) = huge(1.0)
       endif
    enddo

    order = bubble_sort(real(s_walls,kind=db))

    ! Move to the closest plane & check the packet is in the model
    check_wall : do i = 1, n_walls
       iwall = order(i)
       l = s_walls(iwall)

       x_test = x + l*u
       y_test = y + l*v
       z_test = z + l*w

       if (is_in_model(x_test,y_test,z_test)) then
          s = l ! distance to the closest wall
          exit check_wall
       endif

       if (i==n_walls) then
          ! The packet does not reach the model
          icell = 0
          s = 0.0
          return
       endif
    enddo check_wall

    ! Find out the closest cell
    icell = find_Voronoi_cell(iwall, x_test, y_test, z_test)

    ! Move to the cell (if wall is approximate)

    return

  end subroutine move_to_Voronoi_grid

  !----------------------------------------

  subroutine length_Voronoi(id,lambda,Stokes,cell_io,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
    ! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
    ! C. Pinte

    integer, intent(in) :: id,lambda
    integer, intent(inout) :: cell_io
    real(kind=db), dimension(4), intent(in) :: Stokes
    logical, intent(in) :: flag_star, flag_direct_star
    real, intent(inout) :: u,v,w
    real, intent(in) :: extrin
    real, intent(inout) :: xio,yio,zio
    real, intent(out) :: ltot
    logical, intent(out) :: flag_sortie


    logical :: lstop

    real(kind=db) :: extr, tau, opacite !, correct_moins, correct_plus
    integer :: previous_cell, cell, next_cell

    real :: x, y, z, l

    !correct_moins = 1.0_db - prec_grille
    !correct_plus = 1.0_db + prec_grille

    extr = extrin

    previous_cell = 0
    cell = cell_io

    lstop = .false.

    x = xio ; y = yio ; z=zio
    ltot = 0.0

    ! Boucle infinie sur les cellules
    do
       !write(*,*) "I am in cell ", cell, "position", real(x), real(y), real(z)
       call cross_Voronoi_cell(cell, previous_cell, x,y,z, u,v,w, next_cell, l)
       opacite=1.0 !kappa(lambda,cell,1,1) ! TODO !!!

       ! Calcul longeur de vol et profondeur optique dans la cellule
       tau=l*opacite ! opacite constante dans la cellule

       ! Comparaison integrale avec tau
       ! et ajustement longueur de vol evntuellement
       if(tau > extr) then ! On a fini d'integrer
          lstop = .true.
          l = l * (extr/tau) ! on rescale l pour que tau=extr
          ltot=ltot+l
       else ! Il reste extr - tau a integrer dans la cellule suivante
          extr=extr-tau
          ltot=ltot+l
       endif

       ! Update position
       x=x+l*u ; y=y+l*v ; z=z+l*w

       ! Test si on on sort de la routine ou pas
       if (lstop) then ! On a fini d'integrer
          flag_sortie = .false.
          cell_io = cell
          xio=x ; yio=y ;zio=z
          return
       else ! On passe a la cellule suivante
          previous_cell = cell
          cell = next_cell

          if (cell < 0) then ! on sort du volume
             flag_sortie = .true.
             return
          endif
       endif
    enddo ! Boucle sur cellules

    return

  end subroutine length_Voronoi

!----------------------------------------

logical function is_in_model(x,y,z)

  real, intent(in) :: x,y,z

  is_in_model = .false.

  if ((x > wall(1)%x4 *  wall(1)%x1).and.(x < wall(2)%x4 * wall(2)%x1)) then
     if ((y > wall(3)%x4 * wall(3)%x2).and.(y < wall(4)%x4 * wall(4)%x2)) then
        if ((z > wall(5)%x4 * wall(5)%x3).and.(z < wall(6)%x4 * wall(6)%x3)) then
           is_in_model = .true.
        endif
     endif
  endif

  return

end function is_in_model


!----------------------------------------

integer function find_Voronoi_cell(iwall, x,y,z)
  ! Methode debile : boucle sur toutes les cellules pour test

  integer, intent(in) :: iwall
  real, intent(in) :: x, y, z

  real :: dist2, dist2_min, i
  integer :: icell, icell_min

  dist2_min = huge(1.0)
  do i=1, wall(iwall)%n_neighbours
     icell = wall(iwall)%neighbour_list(i)
     dist2 = (Voronoi(icell)%x - x)**2 + (Voronoi(icell)%y - y)**2 + (Voronoi(icell)%z - z)**2

     if (dist2 < dist2_min) then
        icell_min = icell
        dist2_min = dist2
     endif
  enddo

  find_Voronoi_cell = icell_min

  return

end function find_Voronoi_cell

!----------------------------------------

subroutine pos_em_cellule_Voronoi(icell,aleat1,aleat2,aleat3,x,y,z)
! Choisit la position d'emission uniformement dans la cellule
! C. Pinte
! 20/05/14

  implicit none

  integer, intent(in) :: icell
  real(kind=db), intent(in) :: aleat1, aleat2, aleat3
  real, intent(out) :: x,y,z

  real(kind=db) :: u, v, w, srw2, argmt
  integer :: previous_cell, next_cell
  real :: l

  ! Direction aleatoire
  w = 2.0_db * aleat1 - 1.0_db
  srw2 = sqrt(1.0_db-w*w)
  argmt = pi*(2.0_db*aleat2-1.0_db)
  u = srw2 * cos(argmt)
  v = srw2 * sin(argmt)

  ! Distance jusqu'au bord de la cellule
  previous_cell = 0
  x = Voronoi(icell)%x ; y = Voronoi(icell)%y ; z = Voronoi(icell)%z
  call cross_Voronoi_cell(icell, previous_cell, x,y,z, real(u),real(v),real(w), next_cell, l)

  ! Repartition uniforme selon cette direction
  l = l * aleat3**(1./3)
  x=x+l*u ; y=y+l*v ; z=z+l*w

  return

end subroutine pos_em_cellule_Voronoi

end module Voronoi_grid

!----------------------------------------

!program Voronoi
!
!  use Voronoi_grid
!
!  implicit none
!
!  call read_Voronoi(1000000)
!
!end program Voronoi
