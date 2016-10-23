module Voronoi_grid

  use constantes
  use parametres
  use utils, only : bubble_sort, appel_syst
  use naleat, only : seed, stream, gtype
  use opacity, only : volume
  use disk, only : density_file
  use prop_star
  use kdtree2_module

  implicit none

  integer, parameter :: max_wall_neighbours = 100000
  real(kind=db), parameter :: prec = 1.e-6_db

  type Voronoi_cell
     real(kind=db), dimension(3) :: xyz
     integer :: id, first_neighbour, last_neighbour
     logical :: exist
  end type Voronoi_cell

  type Voronoi_wall
     character(len=10) :: type
     real :: x1, x2, x3, x4, x5, x6, x7
     integer :: n_neighbours
     integer, dimension(max_wall_neighbours) :: neighbour_list ! Warning hard coded

     ! Fortran trick to make an array of pointers
     type(kdtree2), pointer :: tree

     ! Plane wall :
     ! x1, x2, x3 : normal
     ! x4 : displacement along the normal
  end type Voronoi_wall

  real(kind=db), dimension(:,:,:), allocatable :: wall_cells ! 3, n_cells_wall, n_walls

  type(Voronoi_cell), dimension(:), allocatable :: Voronoi
  type(Voronoi_wall), dimension(:), allocatable :: wall
  integer, dimension(:), allocatable :: neighbours_list

  integer :: n_walls

  interface
     subroutine voro(n_points, max_neighbours, limits,x,y,z, &
          n_in, volume, first_neighbours,last_neighbours, n_neighbours_tot,neighbours_list, ierr) bind(C, name='voro_C')
       use, intrinsic :: iso_c_binding

       integer(c_int), intent(in), value :: n_points, max_neighbours
       real(c_double), dimension(6), intent(in) :: limits
       real(c_double), dimension(n_points), intent(in) :: x,y,z

       integer(c_int), intent(out) :: n_in, n_neighbours_tot,  ierr
       real(c_double), dimension(n_points), intent(out) :: volume
       integer(c_int), dimension(n_points), intent(out) :: first_neighbours,last_neighbours
       integer(c_int), dimension(n_points * 20), intent(out) :: neighbours_list


     end subroutine voro
  end interface

  contains

    subroutine define_Voronoi_grid()

      return

    end subroutine define_Voronoi_grid

  !************************************************************************

  subroutine Voronoi_tesselation_cmd_line(n_points, x,y,z,  limits)

    integer, intent(in) :: n_points
    real(kind=db), dimension(n_points), intent(in) :: x, y, z
    real(kind=db), dimension(6), intent(in) :: limits

    character(len=512) :: cmd
    integer :: i, syst_status, time1, time2, itime, iVoronoi, alloc_status, ios, icell
    integer :: n_neighbours, n_neighbours_tot, first_neighbour, last_neighbour, k, j, nVoronoi
    real :: time
    real(kind=db) :: x_tmp, y_tmp, z_tmp, vol

    character(len=128) :: slimits

    logical, parameter :: lrun = .true.


    open(unit=1, file="particles.txt", status="replace")
    iVoronoi = 0
    do i=1, n_points
       if ((x(i) > limits(1)).and.(x(i) < limits(2))) then
          if ((y(i) > limits(3)).and.(y(i) < limits(4))) then
             if ((z(i) > limits(5)).and.(z(i) < limits(6))) then
                iVoronoi = iVoronoi + 1
                write(unit=1,fmt="(i7,f15.6,f15.6,f15.6)") iVoronoi, real(x(i)), real(y(i)), real(z(i))
             endif
          endif
       endif
    enddo
    close(unit=1)
    n_cells = iVoronoi

    allocate(Voronoi(n_cells), volume(n_cells), stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error Voronoi structure"
       write(*,*) "Exiting"
       stop
    endif
    Voronoi(:)%exist = .false. ! cells are not defined yet
    volume(:) = 0.0

    write(*,*) n_cells, "particles are in the volume"
    write(slimits,fmt="(f15.6,f15.6,f15.6,f15.6,f15.6,f15.6)") limits(1), limits(2), &
         limits(3), limits(4), limits(5), limits(6)

    call system_clock(time1)
    if (lrun) then
       ! Run voro++ command line for now
       ! while I fix the c++/fortran interface and make all the tests
       write(*,*) "Performing Voronoi tesselation on ", n_cells, "SPH particles"
       cmd = "~/codes/voro++-0.4.6/src/voro++  -v -o -c '%i %q %v %s %n' "//&
            trim(slimits)//&
            " particles.txt ; mv particles.txt.vol Voronoi.txt"
       write(*,*) trim (cmd)
       call appel_syst(cmd,syst_status)
    else
       write(*,*) "Using previous Voronoi tesselation"
    endif

    cmd = "rm -rf nVoronoi.txt ; wc -l Voronoi.txt > nVoronoi.txt"
    call appel_syst(cmd,syst_status)

    open(unit=1,file="nVoronoi.txt",status="old")
    read(1,*) nVoronoi
    close(unit=1)
    write(*,*) "Mesh was tesselated with ", nVoronoi, "cells"
    if (nVoronoi /= n_cells) then
       write(*,*) "*****************************************"
       write(*,*) "WARNING : some particles are not the mesh"
       write(*,*) "*****************************************"
    endif

    call system_clock(time2)
    time=(time2 - time1)/real(time_tick)
    if (time > 60) then
       itime = int(time)
       write (*,'(" Tesselation Time = ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
    else
       write (*,'(" Tesselation Time = ", F5.2, "s")')  time
    endif
    write(*,*) "Voronoi Tesselation done"

    !************************************
    ! We have to read the file now
    !************************************
    n_walls = 6
    write(*,*) "Finding ", n_walls, "walls"
    call init_Voronoi_walls(n_walls, limits)
    write(*,*) "Reading", nVoronoi, " Voronoi cells"

    n_neighbours_tot = 0
    open(unit=1, file="Voronoi.txt", status='old', iostat=ios)

    ! The Voronoi mesh is not defined to start with
    Voronoi(:)%exist = .false.
    Voronoi(:)%first_neighbour = 0
    Voronoi(:)%last_neighbour = 0

    do i=1, nVoronoi
       ! We read the file to figure out the neighbours list
       read(1,*) icell , x_tmp, y_tmp, z_tmp, vol, n_neighbours

       ! We use temporary variables first_neighbour and last_neighbour
       ! as the id of the Voronoi cells might not be consecutive
       if (i == 1) then
          first_neighbour = 1
          last_neighbour = n_neighbours
       else
          first_neighbour = last_neighbour + 1
          last_neighbour = last_neighbour + n_neighbours
       endif
       Voronoi(icell)%first_neighbour = first_neighbour
       Voronoi(icell)%last_neighbour  = last_neighbour

       if (Voronoi(icell)%first_neighbour < 0) then
          write(*,*) "ERROR in Voronoi cell map"
          write(*,*) "icell=", icell,Voronoi(icell)%first_neighbour, Voronoi(icell)%last_neighbour
          write(*,*) "icell-1=", icell-1,Voronoi(icell-1)%first_neighbour, Voronoi(icell-1)%last_neighbour
          write(*,*) "Exiting"
          stop
       endif
       n_neighbours_tot = n_neighbours_tot + n_neighbours
    enddo

    write(*,*) "Neighbours list size =", n_neighbours_tot
    write(*,*) "Trying to allocate", 4*n_neighbours_tot/ 1024.**2, "MB for neighbours list"
    allocate(neighbours_list(n_neighbours_tot))

    ! We read the file again with all the information
    rewind(1)
    do i=1, nVoronoi
       read(1,*) icell, x_tmp, y_tmp, z_tmp, vol, n_neighbours, &
            neighbours_list(Voronoi(icell)%first_neighbour:Voronoi(icell)%last_neighbour)

       Voronoi(icell)%exist = .true.
       Voronoi(icell)%xyz(1) = x_tmp
       Voronoi(icell)%xyz(2) = y_tmp
       Voronoi(icell)%xyz(3) = z_tmp
       Volume(icell) = vol
       ! todo : find the cells touching the walls
       do k=1, n_neighbours
          j = neighbours_list(Voronoi(icell)%first_neighbour + k-1)
          if (j < 0) then ! wall
             j = -j ! wall index
             wall(j)%n_neighbours = wall(j)%n_neighbours+1
             if (wall(j)%n_neighbours > max_wall_neighbours) then
                write(*,*) "ERROR : Voronoi wall", j, "max number of neighbours reached"
                write(*,*) wall(j)%n_neighbours, max_wall_neighbours
                write(*,*) "Exiting"
                stop
             endif
             wall(j)%neighbour_list(wall(j)%n_neighbours) = icell
          endif ! wall
       enddo ! k
    enddo ! i

    close(unit=1)

    write(*,*) "Voronoi volume =", sum(volume(1:nVoronoi))

    write(*,*) "Tesselation finished"

    return

  end subroutine Voronoi_tesselation_cmd_line

!----------------------------------------

  subroutine Voronoi_tesselation(n_points, x,y,z,  limits, check_previous_tesselation)

    integer, intent(in) :: n_points
    real(kind=db), dimension(n_points), intent(in) :: x, y, z
    real(kind=db), dimension(6), intent(in) :: limits
    logical, intent(in) :: check_previous_tesselation

    integer, parameter :: max_neighbours = 20  ! maximum number of neighbours per cell (to build neighbours list)

    real(kind=db), dimension(:), allocatable :: x_tmp, y_tmp, z_tmp
    integer, dimension(:), allocatable :: SPH_id
    real :: time
    integer :: n_in, n_neighbours_tot, ierr, alloc_status, k, j, time1, time2, itime, i, icell, istar, n_sublimate, n_missing_cells
    integer, dimension(:), allocatable :: first_neighbours,last_neighbours

    logical :: is_outside_stars, lcompute

    real(kind=db), dimension(n_etoiles) :: deuxr2_star
    real(kind=db) :: dx, dy, dz, dist2

    n_walls = 6
    write(*,*) "Finding ", n_walls, "walls"
    call init_Voronoi_walls(n_walls, limits)

    ! For initial poistion of rays in ray-tracing mode
    Rmax = sqrt( (limits(2)-limits(1))**2 + (limits(4)-limits(3))**2 + (limits(6)-limits(5))**2 )

    allocate(x_tmp(n_points+n_etoiles), y_tmp(n_points+n_etoiles), z_tmp(n_points+n_etoiles), &
         SPH_id(n_points+n_etoiles),stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error Voronoi temp arrays"
       write(*,*) "Exiting"
       stop
    endif

    do istar=1, n_etoiles
       deuxr2_star(istar) = (2*etoile(istar)%r)**2
    enddo

    ! Filtering particles outside the limits and inside stars
    icell = 0
    n_sublimate = 0
    do i=1, n_points
       ! We test if the point is in the model volume
       if ((x(i) > limits(1)).and.(x(i) < limits(2))) then
          if ((y(i) > limits(3)).and.(y(i) < limits(4))) then
             if ((z(i) > limits(5)).and.(z(i) < limits(6))) then

                ! We also test if the edge of the cell can be inside the star
                ! We test for the edge, so the center needs to be at twice the radius
                is_outside_stars = .true.
                loop_stars : do istar=1, n_etoiles
                   dx = x(i) - etoile(istar)%x
                   dy = y(i) - etoile(istar)%y
                   dz = z(i) - etoile(istar)%z

                   if (min(dx,dy,dz) < 2*etoile(istar)%r) then
                      dist2 = dx**2 + dy**2 + dz**2
                      if (dist2 < deuxr2_star(istar)) then
                         is_outside_stars = .false.
                         n_sublimate = n_sublimate + 1
                         exit loop_stars
                      endif
                   endif
                enddo loop_stars

                if (is_outside_stars) then
                   icell = icell + 1
                   SPH_id(icell) = i
                   x_tmp(icell) = x(i) ; y_tmp(icell) = y(i) ; z_tmp(icell) = z(i) ;
                endif
             endif
          endif
       endif
    enddo

    if (n_sublimate > 0) then
       write(*,*) n_sublimate, "particles have been sublimated"
       write(*,*) "Not implemented yet : MCFOST will probably crash !!!!"
    endif

    ! Filtering stars outside the limits
    etoile(:)%out_model = .true.
    etoile(:)%icell = 0
    do i=1, n_etoiles
       if ((etoile(i)%x > limits(1)).and.(etoile(i)%x < limits(2))) then
          if ((etoile(i)%y > limits(3)).and.(etoile(i)%y < limits(4))) then
             if ((etoile(i)%z > limits(5)).and.(etoile(i)%z < limits(6))) then
                icell = icell + 1
                x_tmp(icell) = etoile(i)%x ; y_tmp(icell) = etoile(i)%y ; z_tmp(icell) = etoile(i)%z ;
                etoile(i)%out_model = .false.
                etoile(i)%icell = icell
                SPH_id(icell) = 0
             endif
          endif
       endif
    enddo
    n_cells = icell

    allocate(Voronoi(n_cells), volume(n_cells), first_neighbours(n_cells),last_neighbours(n_cells), stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Allocation error Voronoi structure"
       write(*,*) "Exiting"
       stop
    endif
    volume(:) = 0.0 ; first_neighbours(:) = 0 ; last_neighbours(:) = 0 ;
    Voronoi(:)%exist = .true. ! we filter before, so all the cells should exist now
    Voronoi(:)%first_neighbour = 0
    Voronoi(:)%last_neighbour = 0

    write(*,*) "Trying to allocate", 4*n_cells * max_neighbours/ 1024.**2, "MB for neighbours list"
    allocate(neighbours_list(n_cells * max_neighbours), stat=alloc_status)
    if (alloc_status /=0) then
       write(*,*) "Error when allocating neighbours list"
       write(*,*) "Exiting"
    endif
    neighbours_list = 0

    do icell=1, n_cells
       Voronoi(icell)%xyz(1) = x_tmp(icell)
       Voronoi(icell)%xyz(2) = y_tmp(icell)
       Voronoi(icell)%xyz(3) = z_tmp(icell)
       Voronoi(icell)%id     = SPH_id(icell)
    enddo

    call system_clock(time1)
    write(*,*) "Performing Voronoi tesselation on ", n_cells, "SPH particles"
    if (check_previous_tesselation) then
       call read_saved_Voronoi_tesselation(n_cells,max_neighbours, limits, &
            lcompute, n_in,first_neighbours,last_neighbours,n_neighbours_tot,neighbours_list)
    else
       lcompute = .true.
    endif
    if (lcompute) then
       call voro(n_cells,max_neighbours,limits,x_tmp,y_tmp,z_tmp,  &
            n_in,volume,first_neighbours,last_neighbours,n_neighbours_tot,neighbours_list,ierr)
       if (ierr /= 0) then
          write(*,*) "Voro++ excited with an error", ierr
          write(*,*) "Exiting"
          stop
       endif
       if (check_previous_tesselation) then
          call save_Voronoi_tesselation(limits, n_in, n_neighbours_tot,first_neighbours,last_neighbours,neighbours_list)
       endif
    else
       write(*,*) "Reading previous Voronoi tesselation"
    endif

    ! Conversion to Fortran indices
    Voronoi(:)%first_neighbour = first_neighbours(:) + 1
    Voronoi(:)%last_neighbour = last_neighbours(:) + 1
    deallocate(first_neighbours,last_neighbours)

    ! Setting-up the walls
    n_missing_cells = 0
    do icell=1, n_cells
       if (volume(icell) < tiny_real) then
          n_missing_cells = n_missing_cells + 1
       endif

       ! todo : find the cells touching the walls
       do k=Voronoi(icell)%first_neighbour,Voronoi(icell)%last_neighbour
          j = neighbours_list(k)
          if (j < 0) then ! wall
             j = -j ! wall index
             wall(j)%n_neighbours = wall(j)%n_neighbours+1
             if (wall(j)%n_neighbours > max_wall_neighbours) then
                write(*,*) "ERROR : Voronoi wall", j, "max number of neighbours reached"
                write(*,*) wall(j)%n_neighbours, max_wall_neighbours
                write(*,*) "Exiting"
                stop
             endif
             wall(j)%neighbour_list(wall(j)%n_neighbours) = icell
          endif ! wall
       enddo ! k
    enddo ! i

    write(*,*) "Building the kd-trees for the model walls"
    call build_wall_kdtrees()

    if (n_missing_cells > 0) then
       write(*,*) "*******************************************"
       write(*,*) "WARNING:", n_missing_cells, "are missing"
       write(*,*) "*******************************************"
    endif
    !write(*,*) "icell #", icell, "is missing", x_tmp(icell), y_tmp(icell), z_tmp(icell)

    call system_clock(time2)
    time=(time2 - time1)/real(time_tick)
    if (time > 60) then
       itime = int(time)
       write (*,'(" Tesselation Time = ", I3, "h", I3, "m", I3, "s")')  itime/3600, mod(itime/60,60), mod(itime,60)
    else
       write (*,'(" Tesselation Time = ", F5.2, "s")')  time
    endif
    write(*,*) "Found", n_in, "cells"

    if (n_in /= n_cells) then
       write(*,*) "*******************************************"
       write(*,*) "ERROR : some particles are not in the mesh"
       write(*,*) "*******************************************"
       stop
    endif

    write(*,*) "Neighbours list size =", n_neighbours_tot
    write(*,*) "Average number of neighbours =", real(n_neighbours_tot)/n_cells
    write(*,*) "Voronoi volume =", sum(volume)

    deallocate(x_tmp, y_tmp, z_tmp, SPH_id)
    !i = 1
    !write(*,*) i, real(Voronoi(i)%xyz(1)), real(Voronoi(i)%xyz(2)), real(Voronoi(i)%xyz(3)), real(volume(i)), Voronoi(i)%last_neighbour-Voronoi(i)%first_neighbour+1, neighbours_list(Voronoi(i)%first_neighbour:Voronoi(i)%last_neighbour)

    return

  end subroutine Voronoi_tesselation

  !**********************************************************

  subroutine save_Voronoi_tesselation(limits, n_in, n_neighbours_tot, first_neighbours,last_neighbours,neighbours_list)

    real(kind=db), intent(in), dimension(6) :: limits
    integer, intent(in) :: n_in, n_neighbours_tot
    integer, dimension(:), intent(in) :: first_neighbours,last_neighbours, neighbours_list
    character(len=512) :: filename
    character(len=40) :: voronoi_sha1

    call get_voronoi_sha1(density_file, voronoi_sha1)

    filename = "_voronoi.tmp"
    open(1,file=filename,status='replace',form='unformatted')
    ! todo : add id for the SPH file : filename + sha1 ??  + limits !!
    write(1) voronoi_sha1, limits, n_in, n_neighbours_tot, volume, first_neighbours,last_neighbours, neighbours_list
    close(1)

    return

  end subroutine save_Voronoi_tesselation


  !**********************************************************

  subroutine read_saved_Voronoi_tesselation(n_cells,max_neighbours, limits, &
       lcompute, n_in,first_neighbours,last_neighbours,n_neighbours_tot,neighbours_list)

    integer, intent(in) :: n_cells,max_neighbours
    real(kind=db), intent(in), dimension(6) :: limits

    logical, intent(out) :: lcompute
    integer, intent(out) :: n_in, n_neighbours_tot
    integer, dimension(n_cells), intent(out) :: first_neighbours,last_neighbours
    integer, dimension(n_cells*max_neighbours), intent(out) :: neighbours_list

    character(len=512) :: filename
    integer :: ios

    character(len=40) :: voronoi_sha1, voronoi_sha1_saved
    real(kind=db), dimension(6) :: limits_saved

    lcompute = .true.

    call get_voronoi_sha1(density_file, voronoi_sha1)

    filename = "_voronoi.tmp"
    ! check if there is a Voronoi file
    ios = 0
    open(1,file=filename,status='old',form='unformatted',iostat=ios)
    if (ios /= 0)  then
       close(unit=1)
       return
    endif

    ! read the saved Voronoi mesh
    read(1,iostat=ios) voronoi_sha1_saved, limits_saved, n_in, n_neighbours_tot, volume, &
         first_neighbours,last_neighbours, neighbours_list
    close(unit=1)
    if (ios /= 0) then ! if some dimension changed
       return
    endif

    ! We are using the same file with the same limits
    if (voronoi_sha1 == voronoi_sha1_saved) then
       if (maxval(abs(limits-limits_saved)) < 1e-3) then
          lcompute = .false.
       endif
    endif

    return

  end subroutine read_saved_Voronoi_tesselation

  !**********************************************************

  subroutine get_voronoi_sha1(filename, voronoi_sha1)

    use system

    character(len=512), intent(in) :: filename
    character(len=40), intent(out) :: voronoi_sha1
    character(len=512) :: cmd
    integer :: ios, syst_status

    if (operating_system=="Linux ") then
       cmd = "sha1sum  "//trim(filename)//" > voronoi.sha1"
    else if (operating_system=="Darwin") then
       cmd = "openssl sha1 "//trim(filename)//" | awk '{print $2}' > voronoi.sha1"
    else
       write(*,*) "Unknown operatinf system"
       write(*,*) "Can't compute sha1 of "//trim(filename)
       return
    endif
    call appel_syst(cmd, syst_status)
    open(unit=1, file="voronoi.sha1", status='old',iostat=ios)
    read(1,*,iostat=ios) voronoi_sha1
    close(unit=1,status="delete",iostat=ios)

    return

  end subroutine get_voronoi_sha1

  !**********************************************************

  subroutine test_walls()

    integer :: iwall

    do iwall=1, n_walls
       write(*,*) "wall#", iwall, wall(iwall)%n_neighbours
    enddo

    stop
    return

  end subroutine test_walls

  !----------------------------------------


  subroutine test_emission()

#include "sprng_f.h"

    integer, parameter :: n_sample = 10000000

    real :: rand, rand2, rand3
    real(kind=db) :: x, y, z
    integer :: icell, i, id, k


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
    write(*,*) "testing emission position", Voronoi(icell)%xyz(1), Voronoi(icell)%xyz(2), Voronoi(icell)%xyz(3)

    do k=1, n_sample
       rand  = sprng(stream(id))
       rand2 = sprng(stream(id))
       rand3 = sprng(stream(id))

       call pos_em_cellule_Voronoi(icell,rand,rand2,rand3, x,y,z)
       write(*,*) x,y,z
    enddo

  end subroutine test_emission

  !----------------------------------------

  subroutine cross_Voronoi_cell(x,y,z, u,v,w, icell, previous_cell, x1,y1,z1, next_cell, s)

    integer, intent(in) :: icell, previous_cell
    real(kind=db), intent(in) :: x,y,z, u,v,w

    real(kind=db), intent(out) ::  s
    integer, intent(out) :: next_cell

    real(kind=db) :: s_tmp, den
    integer :: i, id_n

    real(kind=db), intent(out) :: x1, y1, z1

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real, dimension(3) :: n, p, r, k

    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    s = 1e30 !huge_real
    next_cell = 0

    nb_loop : do i=Voronoi(icell)%first_neighbour, Voronoi(icell)%last_neighbour
       id_n = neighbours_list(i) ! id du voisin

       if (id_n==previous_cell) cycle nb_loop

       if (id_n > 0) then ! cellule
          ! unnormalized vector to plane
          n(:) = Voronoi(id_n)%xyz(:) - Voronoi(icell)%xyz(:)

          ! test direction
          den = dot_product(n, k)
          if (den <= 0.) cycle nb_loop ! car s_tmp sera < 0

          ! point on the plane
          p(:) = 0.5 * (Voronoi(id_n)%xyz(:) + Voronoi(icell)%xyz(:))

          s_tmp = dot_product(n, p-r) / den

          if (s_tmp < 0.) s_tmp = huge(1.0)
       else ! i < 0 ; le voisin est un wall
          s_tmp = distance_to_wall(x,y,z, u,v,w, -id_n) ;

          ! si c'est le wall d'entree : peut-etre a faire sauter en sauvegardant le wall d'entree
          if (s_tmp < 0.) s_tmp = huge(1.0)
       endif

       if (s_tmp < s) then
          s = s_tmp
          next_cell = id_n
       endif

    enddo nb_loop ! i

    x1 = x + u*s
    y1 = y + v*s
    z1 = z + w*s

    if (next_cell == 0) then ! there is a rounding-off error somewehere
       ! We correct the cell index and do not move the packet
       x1 = x ; y1 = y ; z1 = z ; s = 0.0
       if (is_in_volume(x,y,z)) then
          call indice_cellule_voronoi(x,y,z, next_cell)
          if (icell == next_cell) then ! that means we are out of the grid already
             next_cell = -1 ! the exact index does not matter
          endif
       else
          next_cell = -1 ! the exact index does not matter
       endif
    endif

    return

  end subroutine cross_Voronoi_cell

  !----------------------------------------

  subroutine init_Voronoi_walls(n_walls, limits)

    integer, intent(in) :: n_walls
    real(kind=db), dimension(n_walls), intent(in) :: limits
    integer :: iwall

    if (n_walls /= 6) then
       write(*,*) "ERROR: n_walls must be 6 in Voronoi grid"
       write(*,*) "Exiting"
       stop
    endif

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

    wall(1)%x1 = -1 ; wall(1)%x2 = 0  ; wall(1)%x3 = 0  ; wall(1)%x4 = limits(1)  !xmin
    wall(2)%x1 =  1 ; wall(2)%x2 = 0  ; wall(2)%x3 = 0  ; wall(2)%x4 = limits(2)  !xmax
    wall(3)%x1 =  0 ; wall(3)%x2 = -1 ; wall(3)%x3 = 0  ; wall(3)%x4 = limits(3)  !ymin
    wall(4)%x1 =  0 ; wall(4)%x2 = 1  ; wall(4)%x3 = 0  ; wall(4)%x4 = limits(4)  !ymax
    wall(5)%x1 =  0 ; wall(5)%x2 = 0  ; wall(5)%x3 = -1 ; wall(5)%x4 = limits(5)  !zmin
    wall(6)%x1 =  0 ; wall(6)%x2 = 0  ; wall(6)%x3 = 1  ; wall(6)%x4 = limits(6)  !zmax

    return

  end subroutine init_Voronoi_walls

  !----------------------------------------


  real(kind=db) function distance_to_wall(x,y,z, u,v,w, iwall)
    ! Mur plan pour le moment : meme algorithme que cross Voronoi cell

    real(kind=db), intent(in) :: x,y,z, u,v,w
    integer, intent(in) :: iwall

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real(kind=db), dimension(3) :: n, p, r, k

    real :: den

    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    n(1) = wall(iwall)%x1 ; n(2) = wall(iwall)%x2 ; n(3) = wall(iwall)%x3

    p = wall(iwall)%x4 * abs(n)

    den = dot_product(n, k) ! le signe depend du sens de propagation par rapport a la normale

    if (abs(den) > tiny_real) then
       distance_to_wall = dot_product(n, p-r) / den
    else
       distance_to_wall = huge(1.0)
    endif

    return

  end function distance_to_wall

  !----------------------------------------

  subroutine move_to_grid_Voronoi(id, x,y,z, u,v,w, icell, lintersect)

    integer, intent(in) :: id
    real(kind=db), intent(inout) :: x,y,z
    real(kind=db), intent(in) :: u,v,w

    integer, intent(out) :: icell
    logical, intent(out) :: lintersect

    logical, dimension(n_walls) :: intersect
    real(kind=db), dimension(n_walls) :: s_walls
    integer, dimension(n_walls) :: order

    real(kind=db) :: s, l, x_test, y_test, z_test
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

    ! Move to the closest plane & check that the packet is in the model
    check_wall : do i = 1, n_walls
       iwall = order(i)
       l = s_walls(iwall)

       x_test = x + l*u
       y_test = y + l*v
       z_test = z + l*w

       if (is_in_volume(x_test,y_test,z_test)) then
          s = l ! distance to the closest wall
          exit check_wall
       endif

       if (i==n_walls) then
          ! The packet does not reach the model
          icell = 0
          lintersect = .false.
          return
       endif
    enddo check_wall

    lintersect = .true.
    x = x_test ; y = y_test ; z = z_test

    ! Find out the closest cell
    icell = find_Voronoi_cell(id, iwall, x, y, z)

    return

  end subroutine move_to_grid_Voronoi

  !----------------------------------------

  subroutine length_Voronoi(id,lambda,Stokes,cell_io,xio,yio,zio,u,v,w,flag_star,flag_direct_star,extrin,ltot,flag_sortie)
    ! Ne met a jour xio, ... que si le photon ne sort pas de la nebuleuse (flag_sortie=1)
    ! C. Pinte

    integer, intent(in) :: id,lambda
    integer, intent(inout) :: cell_io
    real(kind=db), dimension(4), intent(in) :: Stokes
    logical, intent(in) :: flag_star, flag_direct_star
    real(kind=db), intent(inout) :: u,v,w
    real(kind=db), intent(in) :: extrin
    real(kind=db), intent(inout) :: xio,yio,zio
    real(kind=db), intent(out) :: ltot
    logical, intent(out) :: flag_sortie


    logical :: lstop

    real(kind=db) :: extr, tau, opacite !, correct_moins, correct_plus
    integer :: previous_cell, cell, next_cell

    real(kind=db) :: x, y, z, l, x1,y1,z1

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
       call cross_Voronoi_cell(x,y,z, u,v,w, cell, previous_cell, x1,y1,z1, next_cell, l)
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
       x=x1 ; y=y1 ; z=z1

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

pure logical function test_exit_grid_Voronoi(icell, x,y,z)

  integer, intent(in) :: icell
  real(kind=db), intent(in) :: x,y,z

  if (icell < 0) then
     test_exit_grid_Voronoi = .true.
  else
     test_exit_grid_Voronoi = .false.
  endif

  return

end function test_exit_grid_Voronoi

!----------------------------------------

pure logical function is_in_volume(x,y,z)

  real(kind=db), intent(in) :: x,y,z

  is_in_volume  = .false.
  if ((x > wall(1)%x4).and.(x < wall(2)%x4)) then
     if ((y > wall(3)%x4).and.(y < wall(4)%x4)) then
        if ((z > wall(5)%x4).and.(z < wall(6)%x4)) then
           is_in_volume  = .true.
        endif
     endif
  endif

  return

end function is_in_volume

!----------------------------------------

integer function find_Voronoi_cell_brute_force(id, iwall, x,y,z)
  ! Methode debile : boucle sur toutes les cellules pour test

  integer, intent(in) :: id, iwall
  real(kind=db), intent(in) :: x, y, z

  real :: dist2, dist2_min
  integer :: icell, icell_min, i

  dist2_min = huge(1.0)
  do i=1, wall(iwall)%n_neighbours
     icell = wall(iwall)%neighbour_list(i)
     dist2 = (Voronoi(icell)%xyz(1) - x)**2 + (Voronoi(icell)%xyz(2) - y)**2 + (Voronoi(icell)%xyz(3) - z)**2

     if (dist2 < dist2_min) then
        icell_min = icell
        dist2_min = dist2
     endif
  enddo

  find_Voronoi_cell_brute_force = icell_min

  return

end function find_Voronoi_cell_brute_force

!----------------------------------------

subroutine pos_em_cellule_Voronoi(icell,aleat1,aleat2,aleat3, x,y,z)
! Choisit la position d'emission uniformement dans la cellule
! C. Pinte
! 20/05/14

  implicit none

  integer, intent(in) :: icell
  real, intent(in) :: aleat1, aleat2, aleat3
  real(kind=db), intent(out) :: x,y,z

  real(kind=db) :: u, v, w, srw2, argmt, x1,y1,z1
  integer :: previous_cell, next_cell
  real(kind=db) :: l

  ! Direction aleatoire
  w = 2.0_db * aleat1 - 1.0_db
  srw2 = sqrt(1.0_db-w*w)
  argmt = pi*(2.0_db*aleat2-1.0_db)
  u = srw2 * cos(argmt)
  v = srw2 * sin(argmt)

  ! Distance jusqu'au bord de la cellule
  previous_cell = 0
  x = Voronoi(icell)%xyz(1) ; y = Voronoi(icell)%xyz(2) ; z = Voronoi(icell)%xyz(3)
  call cross_Voronoi_cell(x,y,z, u,v,w, icell, previous_cell, x1,y1,z1, next_cell, l)

  ! Repartition uniforme selon cette direction
  l = l * aleat3**(1./3)
  !x=x+l ; y=y+l*v ; z=z+l*w ! emission au centre cellule si cette ligne est commentee

  return

end subroutine pos_em_cellule_Voronoi


!----------------------------------------

subroutine indice_cellule_Voronoi(xin,yin,zin, icell)

    implicit none

    real(kind=db), intent(in) :: xin,yin,zin
    integer, intent(out) :: icell

    integer :: i
    real :: dist2_min, dist2

    dist2_min = huge(1.0)
    do i=1, n_cells
       if (Voronoi(i)%exist) then ! testing only on the actual Voronoi cells
          dist2 = (Voronoi(i)%xyz(1) - xin)**2 + (Voronoi(i)%xyz(2) - yin)**2 + (Voronoi(i)%xyz(3) - zin)**2

          if (dist2 < dist2_min) then
             icell = i
             dist2_min = dist2
          endif
       endif
    enddo

    return

  end subroutine indice_cellule_Voronoi

!----------------------------------------

subroutine deallocate_Voronoi()

  if (allocated(Voronoi)) deallocate(Voronoi)
  if (allocated(volume)) deallocate(volume)
  if (allocated(neighbours_list)) deallocate(neighbours_list)
  if (allocated(wall)) deallocate(wall)

  return

end subroutine deallocate_Voronoi

!----------------------------------------

subroutine build_wall_kdtrees()

  integer :: iwall, i, icell, alloc_status

  call allocate_kdtree2_search(nb_proc)

  allocate(wall_cells(3,maxval(wall(:)%n_neighbours),n_walls), stat=alloc_status)
  if (alloc_status /=0) then
     write(*,*) "Allocation error wall kdtrees"
     write(*,*) "Exiting"
     stop
  endif
  wall_cells = 0.

  do iwall=1, n_walls
     ! Filling up the data array for the kdtree
     do i=1, wall(iwall)%n_neighbours
        icell = wall(iwall)%neighbour_list(i)
        wall_cells(:,i,iwall) = Voronoi(icell)%xyz(:)
     enddo
  enddo

  do iwall=1, n_walls
     ! Create the tree, ask for internally rearranged data for speed,
     ! and for output NOT sorted by increasing distance from the query vector
     wall(iwall)%tree => wall_kdtree2_create(wall_cells,iwall,n_points=wall(iwall)%n_neighbours,&
          rearrange=.false.,sort=.false.)
  end do

  return

end subroutine build_wall_kdtrees

!----------------------------------------

integer function find_Voronoi_cell(id, iwall, x,y,z)

  integer, intent(in) :: id, iwall
  real(kind=db), intent(in) :: x, y, z

  real(kind=db), dimension(3), target :: qv
  integer, parameter :: NN = 1 ! we only need the closest neighbour

  type(kdtree2_result), dimension(NN) :: results

  qv(1) = x ; qv(2) = y ; qv(3) = z

  ! Find the first vectors in the tree nearest to 'qv' in euclidean norm
  call kdtree2_n_nearest(id, wall(iwall)%tree,qv,NN,results)

  ! Convert the neighbour index to a Voronoi cell index
  find_Voronoi_cell = wall(iwall)%neighbour_list(results(1)%idx)

  return

end function find_Voronoi_cell

!----------------------------------------


!  subroutine kdtree2_example
!
!    !type(kdtree2), pointer :: tree
!    integer :: N,d
!    real(kind=db), allocatable :: mydata(:,:)
!
!    type(kdtree2_result), allocatable :: results(:)
!    ! user sets d, the dimensionality of the Euclidean space
!    ! and N, the number of points in the set.
!    allocate(mydata(d,N))
!    ! note order, d is first, N second.
!    ! read in vectors into mydata(j,i) for j=1..d, i=1..N
!
!    allocate(wall_tree(6))
!
!    wall_tree(1)%tree => kdtree2_create(mydata,rearrange=.true.,sort=.true.)
!    ! Create the tree, ask for internally rearranged data for speed,
!    ! and for output sorted by increasing distance from the
!    ! query vector
!    allocate(results(20))
!    !call kdtree2_n_nearest_around_point(tree,idxin=100,nn=20,correltime=50,results)
!    !call kdtree2_n_nearest_around_point(tree,100,20,50,results)
!
!
!    ! Now the 20 nearest neighbors to mydata(*,100) are in results(:) except
!    ! that points within 50 time units of idxin=50 are not considered as valid neighbors.
!    !
!    write (*,*) "The first 10 near neighbor distances are: ", results(1:10)%dis
!    write (*,*) "The first 10 near neighbor indexes   are: ", results(1:10)%idx
!
!    return
!
!  end subroutine kdtree2_example

end module Voronoi_grid
