module Voronoi_grid

  use constantes
  use mcfost_env
  use parametres
  use utils, only : appel_syst
  use sort, only : Knuth_shuffle, index_quicksort
  use naleat, only : seed, stream, gtype
  use cylindrical_grid, only : volume
  use kdtree2_module
  use messages
  use os

  implicit none
  save

  integer, parameter :: max_wall_neighbours = 100000
  ! n_saved_neighbours = 30 is fine when there is no randomization of particles in mess_up_SPH, otherwise 40.
  ! Todo:  0 does not change speed ... we probably do not even need it
  integer, parameter :: n_saved_neighbours = 0
  real(kind=dp), parameter :: prec = 1.e-6_dp

  type Voronoi_cell
     real(kind=dp), dimension(3) :: xyz, vxyz
     real(kind=dp) :: h ! SPH smoothing lengths or typical cell size in grid based codes
     integer :: id, original_id, first_neighbour, last_neighbour
     logical(kind=lp) :: exist, is_star, was_cut
     logical :: is_star_neighbour!.true. if .not.is_star and a star's neighbour! default is .false..
     logical :: masked
  end type Voronoi_cell

  type Voronoi_wall
     character(len=10) :: type
     real :: x1, x2, x3, x4, x5, x6, x7
     ! Plane wall :
     ! x1, x2, x3 : normal
     ! x4 : displacement along the normal

     integer :: n_neighbours
     integer, dimension(max_wall_neighbours) :: neighbour_list ! Warning max_wall_neighbours is hard coded

     ! Fortran trick to make an array of pointers
     type(kdtree2), pointer :: tree
  end type Voronoi_wall

  type Platonic_Solid
     integer :: n_faces

     ! distance face as a function of h
     real(kind=dp) :: cutting_distance_o_h
     ! distance to vertex as a function of h
     real(kind=dp) :: radius_o_h

     ! Normalized vectors to faces
     real(kind=dp), dimension(:,:), allocatable :: vectors
  end type Platonic_Solid

  real(kind=dp), dimension(:,:,:), allocatable :: wall_cells ! 3, n_cells_wall, n_walls

  type(Voronoi_cell), dimension(:), allocatable :: Voronoi
  real, dimension(:,:), allocatable :: Voronoi_xyz
  real, dimension(:,:,:), allocatable :: Voronoi_neighbour_xyz
  type(Voronoi_wall), dimension(:), allocatable :: wall
  integer, dimension(:), allocatable :: neighbours_list

  type(Platonic_Solid) :: PS

  integer :: n_walls

  interface
     subroutine voro(n_points, max_neighbours, limits,x,y,z,h, threshold, n_vectors, cutting_vectors, cutting_distance_o_h, &
          icell_start,icell_end, cpu_id, n_cpu, n_points_per_cpu, &
          n_in, volume, first_neighbours,last_neighbours,n_neighbours,neighbours_list, was_cell_cut, &
          n_stars, stars_cell, stars_radius, stellar_neighb, ierr) bind(C, name='voro_C')
       use, intrinsic :: iso_c_binding

       integer(c_int), intent(in), value :: n_points, max_neighbours, icell_start, icell_end
       integer(c_int), intent(in), value :: cpu_id, n_cpu, n_points_per_cpu, n_vectors, n_stars
       real(c_double), dimension(6), intent(in) :: limits
       real(c_double), dimension(n_points), intent(in) :: x,y,z,h
       real(c_double), intent(in), value :: threshold, cutting_distance_o_h ! defines at which value we decide to cut the cell, and at how many h the cell will be cut
       ! normal vectors to the cutting plane (need to be normalized)
       real(c_double), dimension(3,n_vectors), intent(in) :: cutting_vectors  ! access order is wrong, but only way to pass a 2d array to C with an unfixed dim
       integer(c_int), intent(in) :: stars_cell(n_stars)
       real(c_double), intent(in) :: stars_radius(n_stars)

       integer(c_int), intent(out) :: n_in,  ierr
       integer(c_int), dimension(n_cpu), intent(out) ::  n_neighbours
       real(c_double), dimension(n_points_per_cpu), intent(out) :: volume
       integer(c_int), dimension(n_points_per_cpu), intent(out) :: first_neighbours,last_neighbours
       integer(c_int), dimension(n_points_per_cpu * max_neighbours * n_cpu), intent(out) :: neighbours_list
       logical(c_bool), dimension(n_points_per_cpu), intent(out) :: was_cell_cut
       logical(c_bool), dimension(n_points_per_cpu), intent(out) :: stellar_neighb

     end subroutine voro
  end interface

  contains

    subroutine define_Voronoi_grid()
      ! This is an empty routine as a target for define_grid

      return

    end subroutine define_Voronoi_grid

    !************************************************************************

    subroutine init_Platonic_Solid(n_faces, radius_o_h)

      integer, intent(in) :: n_faces
      real(kind=dp), intent(in) ::  radius_o_h

      real(kind=dp), parameter :: Phi = (1+sqrt(5.0))/2. ! Golden ratio
      real(kind=dp), parameter :: psi = (1-sqrt(5.0))/2. ! Golden ratio


      real(kind=dp) :: radius, dist_to_face,  face_o_edge, f, fPhi, fpsi, volume, volume_relative_to_circumsphere

      if ((n_faces /= 12).and.(n_faces /= 20)) call error("Only dodecahedon and icosahedron are implemented so far")

      PS%n_faces = n_faces
      PS%radius_o_h = radius_o_h ! We save the input value just in case
      if (.not.allocated(PS%vectors)) allocate(PS%vectors(3,n_faces))

      if (n_faces == 12) then ! regular dodecahedron
         ! a is the edge length
         radius = sqrt(3.0)/2. * Phi ! x a  radius of circumscribed sphere
         dist_to_face = Phi**3 / (2*sqrt(Phi**2+1)) ! x a

         volume = (15 + 7*sqrt(5.0))/4. ! x a**3 ! filling factor ~ 0.665

         f = 1.0/sqrt(1.0+Phi*Phi)  ! 1/Norm of vector with components (0,1,Phi) (in any order)
         fPhi = f * Phi

         ! Normalized vectors perpendicular to faces
         PS%vectors(:,1)  = (/0.0_dp,fPhi,f/)
         PS%vectors(:,2)  = (/0.0_dp,-fPhi,f/)
         PS%vectors(:,3)  = (/0.0_dp,fPhi,-f/)
         PS%vectors(:,4)  = (/0.0_dp,-fPhi,-f/)
         PS%vectors(:,5)  = (/f,0.0_dp,fPhi/)
         PS%vectors(:,6)  = (/-f,0.0_dp,fPhi/)
         PS%vectors(:,7)  = (/f,0.0_dp,-fPhi/)
         PS%vectors(:,8)  = (/-f,0.0_dp,-fPhi/)
         PS%vectors(:,9)  = (/fPhi,f,0.0_dp/)
         PS%vectors(:,10) = (/-fPhi,f,0.0_dp/)
         PS%vectors(:,11) = (/fPhi,-f,0.0_dp/)
         PS%vectors(:,12) = (/-fPhi,-f,0.0_dp/)

      else if (n_faces == 20) then ! regular icosahedron
         ! a is the edge length
         radius = 0.5 * sqrt(Phi * sqrt(5.0))  ! x a
         dist_to_face = Phi**2 / (2*sqrt(3.0)) ! x a

         volume = 5. * (3+sqrt(5.0))/12. ! a**3   filling factor is 0.605

         f = 1.0/sqrt(3.0)  ! 1/Norm of vector with components (1,1,1) or (0,Phi,psi) (in any order)
         fPhi = f * Phi
         fpsi = f * psi

         ! Normalized vectors perpendicular to faces
         PS%vectors(:,1)  = (/f,f,f/)
         PS%vectors(:,2)  = (/-f,f,f/)
         PS%vectors(:,3)  = (/f,-f,f/)
         PS%vectors(:,4)  = (/-f,-f,f/)
         PS%vectors(:,5)  = (/f,f,-f/)
         PS%vectors(:,6)  = (/-f,f,-f/)
         PS%vectors(:,7)  = (/f,-f,-f/)
         PS%vectors(:,8)  = (/-f,-f,-f/)
         PS%vectors(:,9)  = (/0.0_dp,fpsi,fPhi/)
         PS%vectors(:,10) = (/0.0_dp,fpsi,-fPhi/)
         PS%vectors(:,11) = (/0.0_dp,-fpsi,fPhi/)
         PS%vectors(:,12) = (/0.0_dp,-fpsi,-fPhi/)
         PS%vectors(:,13) = (/fPhi,0.0_dp,fpsi/)
         PS%vectors(:,14) = (/fPhi,0.0_dp,-fpsi/)
         PS%vectors(:,15) = (/-fPhi,0.0_dp,fpsi/)
         PS%vectors(:,16) = (/-fPhi,0.0_dp,-fpsi/)
         PS%vectors(:,17) = (/fpsi,fPhi,0.0_dp/)
         PS%vectors(:,18) = (/fpsi,-fPhi,0.0_dp/)
         PS%vectors(:,19) = (/-fpsi,fPhi,0.0_dp/)
         PS%vectors(:,20) = (/-fpsi,-fPhi,0.0_dp/)
      endif

      face_o_edge = dist_to_face/radius
      volume_relative_to_circumsphere = volume / (4.*pi/3 * radius**3)

      ! distance of a face / smoothing length (if edges are at radius_o_h x h from the center)
      !PS%cutting_distance_o_h = radius_o_h * face_o_edge ! PS is exactly included in the sphere
      PS%cutting_distance_o_h = radius_o_h * face_o_edge / (volume_relative_to_circumsphere)**(1./3) ! PS has the same volume as sphere of radius h*radius_o_h

      return

    end subroutine init_Platonic_Solid

    !************************************************************************

  subroutine Voronoi_tesselation(n_points, particle_id, x,y,z,h, vx,vy,vz, limits, check_previous_tesselation)

    use iso_fortran_env
    use, intrinsic :: iso_c_binding, only : c_bool
    !$ use omp_lib

    integer, intent(in) :: n_points
    real(kind=dp), dimension(n_points), intent(in) :: x, y, z, h
    real(kind=dp), dimension(:), allocatable, intent(in) :: vx,vy,vz
    integer, dimension(n_points), intent(in) :: particle_id
    real(kind=dp), dimension(6), intent(in) :: limits
    logical, intent(in) :: check_previous_tesselation


    integer, parameter :: max_neighbours = 30  ! maximum number of neighbours per cell (to build temporary neighbours list)

    real(kind=dp), dimension(:), allocatable :: x_tmp, y_tmp, z_tmp, h_tmp
    integer, dimension(:), allocatable :: SPH_id, SPH_original_id
    real :: time, mem
    integer :: n_in, n_neighbours_tot, ierr, alloc_status, k, j, n, time1, time2, itime, i
    integer :: icell, istar, n_sublimate, n_missing_cells, n_cells_per_cpu, nb_proc_voro
    real(kind=dp), dimension(:), allocatable :: V_tmp
    integer, dimension(:), allocatable :: first_neighbours,last_neighbours
    integer, dimension(:), allocatable :: neighbours_list_loc
    integer, dimension(:), allocatable :: n_neighbours ! nb_proc
    logical(c_bool), dimension(:), allocatable :: was_cell_cut, was_cell_cut_tmp
    logical(c_bool), dimension(:), allocatable :: star_neighb_tmp, star_neighb

    logical :: is_outside_stars, lcompute

    real(kind=dp) :: dx, dy, dz, dist2

    integer :: icell_start, icell_end, id, row, l, n_cells_before_stars, n_average_neighb_stars

    real(kind=dp), parameter :: threshold = 3 ! defines at how many h cells will be cut
    character(len=2) :: unit

    logical, parameter :: lrandom = .true.
!     logical :: lrandom
    integer, dimension(:), allocatable :: order,SPH_id2,SPH_original_id2
    real(kind=dp), dimension(:), allocatable :: x_tmp2,y_tmp2,z_tmp2,h_tmp2

!     lrandom = lrandomize_voronoi
!     if (lrandom) then
!     	write(*,*) " Randomized particles, lrandom = True"
!     else
!     	write(*,*) " lrandom = .false."
!     endif

    if (nb_proc > 16) write(*,*) "Using 16 cores for Voronoi tesselation" ! Overheads dominate above 16 cores

    nb_proc_voro = min(16,nb_proc)
    allocate(n_neighbours(nb_proc_voro))

    ! Defining Platonic solid that will be used to cut the wierly shaped Voronoi cells
    call init_Platonic_Solid(12, threshold)

    n_walls = 6
    write(*,*) "Finding ", n_walls, "walls"
    call init_Voronoi_walls(n_walls, limits)

    ! For initial position of rays in ray-tracing mode
    Rmax = sqrt( (limits(2)-limits(1))**2 + (limits(4)-limits(3))**2 + (limits(6)-limits(5))**2 )

    alloc_status = 0
    allocate(x_tmp(n_points+n_etoiles), y_tmp(n_points+n_etoiles), z_tmp(n_points+n_etoiles), h_tmp(n_points+n_etoiles), &
         SPH_id(n_points+n_etoiles), SPH_original_id(n_points+n_etoiles), stat=alloc_status)
    if (alloc_status /=0) call error("Allocation error Voronoi temp arrays")

    ! Filtering particles outside the limits
    icell = 0
    n_sublimate = 0
    do i=1, n_points
       ! We test if the point is in the model volume
       !-> Voronoi cells are now cut at the surface of the star. We only need
       ! to test if a particle is below Rstar.
       if ((x(i) > limits(1)).and.(x(i) < limits(2))) then
          if ((y(i) > limits(3)).and.(y(i) < limits(4))) then
             if ((z(i) > limits(5)).and.(z(i) < limits(6))) then

                ! We also test if the edge of the cell can be inside the star
                is_outside_stars = .true.

                loop_stars : do istar=1, n_etoiles
                   dx = x(i) - etoile(istar)%x
                   dy = y(i) - etoile(istar)%y
                   dz = z(i) - etoile(istar)%z

                   if (min(dx,dy,dz) < etoile(istar)%r) then
                      dist2 = dx**2 + dy**2 + dz**2
                      if (dist2 < etoile(istar)%r**2) then
                         is_outside_stars = .false.
                         n_sublimate = n_sublimate + 1
                         exit loop_stars
                      endif
                   endif
                enddo loop_stars

                if (is_outside_stars) then
                   icell = icell + 1
                   SPH_id(icell) = i ; SPH_original_id(icell) = particle_id(i)
                   x_tmp(icell) = x(i) ; y_tmp(icell) = y(i) ; z_tmp(icell) = z(i) ;  h_tmp(icell) = h(i)
                endif
             endif
          endif
       endif
    enddo
    n_cells_before_stars = icell
    n_cells = icell
    n_cells_per_cpu = (1.0*n_cells) / nb_proc_voro + 1

    ! Randomizing particles
    if (lrandom) then
       allocate(x_tmp2(n_cells), y_tmp2(n_cells), z_tmp2(n_cells), h_tmp2(n_cells), &
            SPH_id2(n_cells), SPH_original_id2(n_cells), order(n_cells), stat=alloc_status)
       if (alloc_status /=0) call error("Allocation error Voronoi temp2 arrays")

       order=[(i,i=1,n_cells)]
       call Knuth_shuffle(order)

       do icell=1, n_cells
          i = order(icell)
          x_tmp2(icell) = x_tmp(i)
          y_tmp2(icell) = y_tmp(i)
          z_tmp2(icell) = z_tmp(i)
          h_tmp2(icell) = h_tmp(i)
          SPH_id2(icell) = SPH_id(i)
          SPH_original_id2(icell) = SPH_original_id(i)
       enddo

       x_tmp(1:n_cells) = x_tmp2(1:n_cells) ; y_tmp(1:n_cells) = y_tmp2 ; z_tmp(1:n_cells) = z_tmp2 ; h_tmp(1:n_cells) = h_tmp2
       SPH_id(1:n_cells) = SPH_id2(1:n_cells) ; SPH_original_id(1:n_cells) = SPH_original_id2
       deallocate(order,x_tmp2,y_tmp2,z_tmp2,h_tmp2,SPH_id2,SPH_original_id2)
    endif

    if (n_sublimate > 0) then
       write(*,*) n_sublimate, "particles have been sublimated"
       write(*,*) "Not implemented yet : MCFOST will probably crash !!!!"
    endif

    ! Filtering stars outside the limits
    etoile(:)%out_model = .true.
    etoile(:)%icell = 0
    icell = n_cells_before_stars
    do i=1, n_etoiles
       ! We test is the star is in the model
       if ((etoile(i)%x > limits(1)).and.(etoile(i)%x < limits(2))) then
          if ((etoile(i)%y > limits(3)).and.(etoile(i)%y < limits(4))) then
             if ((etoile(i)%z > limits(5)).and.(etoile(i)%z < limits(6))) then
                icell = icell + 1
                SPH_id(icell) = 0 ; SPH_original_id(icell) = 0
                x_tmp(icell) = etoile(i)%x ; y_tmp(icell) = etoile(i)%y ; z_tmp(icell) = etoile(i)%z ; h_tmp(icell) = huge_real ;
                etoile(i)%out_model = .false.
                etoile(i)%icell = icell
             endif
          endif
       endif
    enddo
    n_cells = icell

    alloc_status = 0
    allocate(Voronoi(n_cells), Voronoi_xyz(3,n_cells), volume(n_cells), first_neighbours(n_cells),last_neighbours(n_cells), &
         was_cell_cut(n_cells), star_neighb(n_cells), stat=alloc_status)
    if (alloc_status /=0) call error("Allocation error Voronoi structure")
    volume(:) = 0.0 ; first_neighbours(:) = 0 ; last_neighbours(:) = 0 ; was_cell_cut(:) = .false.
    Voronoi(:)%exist = .true. ! we filter before, so all the cells should exist now
    Voronoi(:)%first_neighbour = 0
    Voronoi(:)%last_neighbour = 0
    Voronoi(:)%is_star = .false.
    star_neighb(:) = .false.

    do icell=1, n_cells
       Voronoi(icell)%xyz(1) = x_tmp(icell)
       Voronoi(icell)%xyz(2) = y_tmp(icell)
       Voronoi(icell)%xyz(3) = z_tmp(icell)
       Voronoi(icell)%h      = h_tmp(icell)
       Voronoi_xyz(:,icell)  = Voronoi(icell)%xyz(:)
       Voronoi(icell)%id     = SPH_id(icell) ! this is the id of the particles passed to mcfost
       Voronoi(icell)%original_id = SPH_original_id(icell) ! this is the particle id in phantom
    enddo

    !*************************
    ! Velocities
    !*************************
    if ((lemission_mol).or.(lemission_atom)) then
       do icell=1,n_cells_before_stars
          i = SPH_id(icell)
          Voronoi(icell)%vxyz(1) = vx(i)
          Voronoi(icell)%vxyz(2) = vy(i)
          Voronoi(icell)%vxyz(3) = vz(i)
       enddo
    endif

    do i=1, n_etoiles
       if (etoile(i)%icell > 0) Voronoi(etoile(i)%icell)%is_star = .true.
    enddo

    call system_clock(time1)
    write(*,*) "Performing Voronoi tesselation on ", n_cells, "SPH particles"
    mem =  (1.0*n_cells) * (5*2 + 1) * 4.
    if (mem > 1e9) then
       mem = mem / 1024**3
       unit = "GB"
    else
       mem = mem / 1024**2
       unit = "MB"
    endif
    write(*,*) "mcfost will require ~", mem, unit//" of temporary memory for the tesselation" ! 5 double + 2 int arrays

    if (operating_system == "Darwin" .and. (n_cells > 2e6)) then
       call warning("Voronoi tesselation will likely crash with that many particle on a Mac. Switch to linux")
    endif

    if (check_previous_tesselation) then
       call read_saved_Voronoi_tesselation(n_cells,max_neighbours, limits, &
            lcompute, n_in,first_neighbours,last_neighbours,n_neighbours_tot,neighbours_list,was_cell_cut,star_neighb)
    else
       lcompute = .true.
    endif


    if (lcompute) then
       mem = 4. * n_cells * max_neighbours
       if (mem > 1e9) then
          mem = mem / 1024**3
          unit = "GB"
       else
          mem = mem / 1024**2
          unit = "MB"
       endif
       write(*,*) "Trying to allocate", mem, unit//" for temporary neighbours list"
       alloc_status = 0
       allocate(neighbours_list_loc(n_cells_per_cpu * max_neighbours * nb_proc_voro), stat=alloc_status)
       if (alloc_status /=0) call error("Allocation error neighbours_list_loc")
       neighbours_list_loc = 0

       n_in = 0 ! We initialize value at 0 as we have a reduction + clause
       !$omp parallel default(none) num_threads(nb_proc_voro) &
       !$omp shared(n_cells,limits,x_tmp,y_tmp,z_tmp,h_tmp,nb_proc_voro,n_cells_per_cpu,etoile) &
       !$omp shared(first_neighbours,last_neighbours,neighbours_list_loc,n_neighbours,PS) &
       !$omp private(id,n,icell_start,icell_end,ierr) &
       !$omp private(V_tmp,was_cell_cut_tmp,alloc_status,star_neighb_tmp) &
       !$omp shared(volume,was_cell_cut,star_neighb, n_etoiles) &
       !$omp reduction(+:n_in)
       id = 1
       !$ id = omp_get_thread_num() + 1
       icell_start = (1.0 * (id-1)) / nb_proc_voro * n_cells + 1
       icell_end = (1.0 * (id)) / nb_proc_voro * n_cells

       ! Allocating results array
       alloc_status = 0
       allocate(V_tmp(n_cells), was_cell_cut_tmp(n_cells), star_neighb_tmp(n_cells), stat=alloc_status)
       if (alloc_status /=0) call error("Allocation error before voro++ call")

       call voro(n_cells,max_neighbours,limits,x_tmp,y_tmp,z_tmp,h_tmp, threshold, PS%n_faces, &
            PS%vectors, PS%cutting_distance_o_h, icell_start-1,icell_end-1, id-1,nb_proc_voro,n_cells_per_cpu, &
            n_in,V_tmp,first_neighbours,last_neighbours,n_neighbours,neighbours_list_loc,was_cell_cut_tmp,&
            n_etoiles, etoile(:)%icell-1,etoile(:)%r,star_neighb_tmp, ierr) ! icell & id shifted by 1 for C
       if (ierr /= 0) then
          write(*,*) "Voro++ excited with an error", ierr, "thread #", id
          write(*,*) "Exiting"
          call exit(1)
       endif

       ! Merging outputs from various cores
       n=icell_end-icell_start+1
       volume(icell_start:icell_end) = V_tmp(1:n)
       was_cell_cut(icell_start:icell_end) = was_cell_cut_tmp(1:n)
       star_neighb(icell_start:icell_end) = star_neighb_tmp(1:n)

       deallocate(V_tmp, was_cell_cut_tmp,star_neighb_tmp)
       !$omp end parallel

       !-----------------------------------------------------------
       ! Merging the neighbour arrays from the different threads
       !-----------------------------------------------------------
       ! We need to shift the indices of the neighbours
       do id=2, nb_proc_voro
          icell_start = (1.0 * (id-1)) / nb_proc_voro * n_cells + 1
          icell_end = (1.0 * (id)) / nb_proc_voro * n_cells

          ! Pointers to the first and last neighbours of the cell
          first_neighbours(icell_start:icell_end) = first_neighbours(icell_start:icell_end) + last_neighbours(icell_start-1) + 1
          last_neighbours(icell_start:icell_end)  = last_neighbours(icell_start:icell_end)  + last_neighbours(icell_start-1) + 1
       enddo

       ! Compacting neighbours list
       n_neighbours_tot = sum(n_neighbours)
       mem = 4. * n_neighbours_tot
       if (mem > 1e9) then
          mem = mem / 1024**3
          unit = "GB"
       else
          mem = mem / 1024**2
          unit = "MB"
       endif
       write(*,*) "Trying to allocate", mem, unit//" for neighbours list"
       allocate(neighbours_list(n_neighbours_tot), stat=alloc_status)
       if (alloc_status /=0) call error("Allocation error neighbours list")

       row = n_cells_per_cpu * max_neighbours ;
       k = 0 ;
       do id=1, nb_proc_voro
          do i=1, n_neighbours(id)
             k = k+1
             neighbours_list(k) = neighbours_list_loc(row * (id-1) + i)
          enddo
       enddo
       write(*,*) "Freeing temporary neighbours list"
       deallocate(neighbours_list_loc)

       if (check_previous_tesselation) then
          call save_Voronoi_tesselation(limits, n_in, n_neighbours_tot, first_neighbours, last_neighbours, &
               neighbours_list,was_cell_cut,star_neighb)
       endif
    else !lcompute
       write(*,*) "Reading previous Voronoi tesselation"
    endif
    write(*,*) "Tesselation done"

    ! Conversion to Fortran indices
    Voronoi(:)%first_neighbour = first_neighbours(:) + 1
    Voronoi(:)%last_neighbour = last_neighbours(:) + 1
    deallocate(first_neighbours,last_neighbours)

    Voronoi(:)%was_cut = was_cell_cut(:)
    deallocate(was_cell_cut)
    Voronoi(:)%is_star_neighbour = star_neighb(:)
    deallocate(star_neighb)

    n_average_neighb_stars = 0
    do icell=1,n_cells
       	if (Voronoi(icell)%is_star_neighbour) then
       		n_average_neighb_stars = n_average_neighb_stars + 1
       	endif
    enddo
    if (n_average_neighb_stars>0) write(*,'("mean number of stellar neighbours="(1I5))') &
    	int(real(size(pack(Voronoi(:)%is_star_neighbour,mask=Voronoi(:)%is_star_neighbour)))/real(n_etoiles))


    ! Saving position of the first neighbours to save time on memory access
    ! Warning this seems to cost a fair bit of time
    allocate(Voronoi_neighbour_xyz(3,n_saved_neighbours,n_cells)) ! tried 4 to align vectors, does not help
    do icell=1, n_cells
       l=0
       do k=Voronoi(icell)%first_neighbour,Voronoi(icell)%last_neighbour
          l = l+1
          if (l <= n_saved_neighbours) then
             j = neighbours_list(k)
             if (j> n_cells) then
                write(*,*) "icell =", icell, "k=", k, "j=", j, "n_cells=", n_cells
                call error("Voronoi neighbour list index is invalid")
             endif
             if (j>0) Voronoi_neighbour_xyz(1:3,l,icell) = Voronoi_xyz(:,j)  ! this is slow
          endif
       enddo
    enddo

    ! Setting-up the walls
    n_missing_cells = 0
    do icell=1, n_cells
       ! We check first that there is no issue in the tesselation
       if (volume(icell) < tiny_real) then
          n_missing_cells = n_missing_cells + 1
          write(*,*) "WARNING: cell #", icell, "is missing"
          write(*,*) "original id =", SPH_original_id(icell)
          write(*,*) "xyz=", x_tmp(icell), y_tmp(icell), z_tmp(icell)
          write(*,*) "volume =", volume(icell)
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
                call exit(1)
             endif
             wall(j)%neighbour_list(wall(j)%n_neighbours) = icell
          endif ! wall
       enddo ! k
    enddo ! icell

    if (n_missing_cells > 0) then
       write(*,*) "*******************************************"
       if (n_missing_cells == 1) then
           write(*,*) "WARNING: 1 cell is missing"
       else
          write(*,*) "WARNING:", n_missing_cells, "cells are missing"
       endif
       write(*,*) "*******************************************"
    endif

    write(*,*) "Building the kd-trees for the model walls"
    call build_wall_kdtrees()

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
       write(*,*) "n_cells =", n_cells
       write(*,*) "n in tesselation =", n_in
       call error("some particles are not in the mesh")
    endif

    write(*,*) "Neighbours list size =", n_neighbours_tot
    write(*,*) "Average number of neighbours =", real(n_neighbours_tot)/n_cells
    write(*,*) "Voronoi volume =", sum(volume)

    !i = 1
    !write(*,*) i, real(Voronoi(i)%xyz(1)), real(Voronoi(i)%xyz(2)), real(Voronoi(i)%xyz(3)), real(volume(i)), Voronoi(i)%last_neighbour-Voronoi(i)%first_neighbour+1, neighbours_list(Voronoi(i)%first_neighbour:Voronoi(i)%last_neighbour)

    return

  end subroutine Voronoi_tesselation

  !**********************************************************

  subroutine save_Voronoi_tesselation(limits, n_in, n_neighbours_tot, first_neighbours,&
   last_neighbours,neighbours_list,was_cell_cut,star_neighb)

    use, intrinsic :: iso_c_binding, only : c_bool

    real(kind=dp), intent(in), dimension(6) :: limits
    integer, intent(in) :: n_in, n_neighbours_tot
    integer, dimension(:), intent(in) :: first_neighbours,last_neighbours, neighbours_list
    logical(c_bool), dimension(:), intent(in) :: was_cell_cut, star_neighb
    character(len=512) :: filename
    character(len=40) :: voronoi_sha1

    call get_voronoi_sha1(density_files(1), voronoi_sha1)

    filename = trim(tmp_dir)//"_voronoi.tmp"
    open(1,file=filename,status='replace',form='unformatted')
    ! todo : add id for the SPH file : filename + sha1 ??  + limits !!
    write(1) voronoi_sha1, limits, n_in, n_neighbours_tot
    write(1) volume, first_neighbours,last_neighbours, neighbours_list, was_cell_cut, star_neighb
    close(1)

    return

  end subroutine save_Voronoi_tesselation


  !**********************************************************

  subroutine read_saved_Voronoi_tesselation(n_cells,max_neighbours, limits, &
       lcompute, n_in,first_neighbours,last_neighbours,n_neighbours_tot,neighbours_list,was_cell_cut,star_neighb)

    use, intrinsic :: iso_c_binding, only : c_bool

    integer, intent(in) :: n_cells,max_neighbours
    real(kind=dp), intent(in), dimension(6) :: limits

    logical, intent(out) :: lcompute
    integer, intent(out) :: n_in, n_neighbours_tot
    integer, dimension(n_cells), intent(out) :: first_neighbours,last_neighbours
    integer, dimension(:), allocatable, intent(out) :: neighbours_list
    logical(c_bool), dimension(n_cells), intent(out) :: was_cell_cut, star_neighb

    character(len=512) :: filename
    integer :: ios, alloc_status
    real :: mem

    character(len=40) :: voronoi_sha1, voronoi_sha1_saved
    real(kind=dp), dimension(6) :: limits_saved
    character(len=2) :: unit

    lcompute = .true.
    if (lrandomize_azimuth) return

    call get_voronoi_sha1(density_files(1), voronoi_sha1)

    filename = trim(tmp_dir)//"_voronoi.tmp"
    ! check if there is a Voronoi file
    ios = 0
    open(1,file=filename,status='old',form='unformatted',iostat=ios)
    if (ios /= 0)  then
       close(unit=1)
       return
    endif

    ! read the saved Voronoi mesh
    read(1,iostat=ios) voronoi_sha1_saved, limits_saved, n_in, n_neighbours_tot
    mem = 4. * n_neighbours_tot
    if (mem > 1e9) then
       mem = mem / 1024**3
       unit = "GB"
    else
       mem = mem / 1024**2
       unit = "MB"
    endif

    alloc_status = 0
    write(*,*) "Trying to allocate", mem, unit//" for neighbours list"
    allocate(neighbours_list(n_neighbours_tot), stat=alloc_status)
    if (alloc_status /=0) call error("Allocating neighbours list in read_saved_tesselation")
    read(1,iostat=ios) volume, first_neighbours,last_neighbours, neighbours_list, was_cell_cut, star_neighb
    close(unit=1)

    if (ios == 0) then ! some dimension
       ! We are using the same file with the same limits
       if (voronoi_sha1 == voronoi_sha1_saved) then
          if (maxval(abs(limits-limits_saved)) < 1e-3) then
             lcompute = .false.
          else
             write(*,*) "Limits have changed"
          endif
       else
          write(*,*) "sha1 has changed"
       endif
    else
       write(*,*) "Problem reading file"
    endif

    if (lcompute) then
       write(*,*) "Tesselation needs to be recomputed, deallocating neighbours list"
       deallocate(neighbours_list)
    endif

    return

  end subroutine read_saved_Voronoi_tesselation

  !**********************************************************

  subroutine get_voronoi_sha1(filename, voronoi_sha1)

    use os

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

    call exit(0)
    return

  end subroutine test_walls

  !----------------------------------------


  subroutine test_emission()

#include "sprng_f.h"

    integer, parameter :: n_sample = 10000000

    real :: rand, rand2, rand3
    real(kind=dp) :: x, y, z
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

  subroutine cross_Voronoi_cell(x,y,z, u,v,w, icell, previous_cell, x1,y1,z1, next_cell, s, s_contrib, s_void_before)

    integer, intent(in) :: icell, previous_cell
    real(kind=dp), intent(in) :: x,y,z, u,v,w

    real(kind=dp), intent(out) :: x1, y1, z1, s, s_contrib, s_void_before
    integer, intent(out) :: next_cell

    real(kind=dp) :: s_tmp, den, num, s_entry, s_exit
    integer :: i, id_n, l, ifirst, ilast

    integer :: i_star
    real(kind=dp) :: d_to_star
    logical :: is_a_star_neighbour

    real(kind=dp) :: b, c, delta, rac, s1, s2, h
    real(kind=dp), dimension(3) :: delta_r

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real, dimension(3) :: n, p, r, k, r_cell, r_neighbour
    logical :: was_cut

    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    s = 1e30 !huge_real
    next_cell = 0

    ! We do all the access to Voronoi(icell) now
    r_cell(:) = Voronoi(icell)%xyz(:)
    ifirst = Voronoi(icell)%first_neighbour
    ilast = Voronoi(icell)%last_neighbour
    was_cut = Voronoi(icell)%was_cut
    h = Voronoi(icell)%h
    !Store the id of the neighbouring stars instead ?
    is_a_star_neighbour = Voronoi(icell)%is_star_neighbour

    l=0
    nb_loop : do i=ifirst,ilast
       l = l+1
       id_n = neighbours_list(i) ! id du voisin

       if (id_n==previous_cell) cycle nb_loop

       if (id_n > 0) then ! cellule
          if (l <= n_saved_neighbours) then ! we used an ordered array to limit cache misses
             r_neighbour(:) = Voronoi_neighbour_xyz(:,l,icell)
          else
             r_neighbour(:) = Voronoi_xyz(:,id_n)
          endif

          ! unnormalized vector to plane
          n(:) = r_neighbour(:) - r_cell(:)

          ! test direction
          den = dot_product(n, k)
          if (den <= 0.) cycle nb_loop ! car s_tmp sera < 0

          ! point on the plane
          p(:) = 0.5 * (r_neighbour(:) + r_cell(:))

          s_tmp = dot_product(n, p-r) / den

          if (s_tmp < 0.) s_tmp = huge(1.0)
       else ! id_n < 0 ; le voisin est un wall
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

    if (was_cut) then
       ! We check where the packet intersect the sphere of radius Voronoi(icell)%h * PS%cutting_distance_o_h
       ! centered on the center of the cell
       delta_r = r - r_cell(:)

       	b = dot_product(delta_r,k)
       	c = dot_product(delta_r,delta_r) - (h * PS%cutting_distance_o_h)**2
       	delta = b*b - c

       	if (delta < 0.) then ! the packet never encounters the sphere
         	 s_void_before = s
          	s_contrib = 0.0_dp
       	else ! the packet encounters the sphere
          	rac = sqrt(delta)
          	s1 = -b - rac
          	s2 = -b + rac

          	if (s1 < 0) then ! we already entered the sphere
             	if (s2 < 0) then ! we already exited the sphere
                	s_void_before = s
                	s_contrib = 0.0_dp
             	else ! We are still in the sphere and will exit it
                	s_void_before = 0.0_dp
                	s_contrib = min(s2,s)
             	endif
          	else ! We will enter in the sphere (both s1 and s2 are > 0)
             	if (s1 < s) then ! We will enter the sphere in this cell
                	s_void_before = s1
                	s_contrib = min(s2,s) - s1
             	else ! We will not enter the sphere in this sphere
                	s_void_before = s
                	s_contrib = 0.0_dp
             	endif
          	endif
       	endif ! delta < 0
    else ! the cell was not cut
       s_void_before = 0.0_dp
       s_contrib = s
    endif


    if (is_a_star_neighbour) then
       d_to_star = distance_to_star(x,y,z,u,v,w,i_star)
       ! It is a neighbour so if the star is intersected (i_star>0) it might be the next cell.
       ! Otherwise, i_star == 0 (does not intersect the star in this direction (u,v,w)).
       if (i_star > 0) then
          if (d_to_star < s) then ! indeed a star, we use d_to_stars and set next_cell
             s_contrib = d_to_star
             next_cell = etoile(i_star)%icell
          endif
       endif
    endif

    return

  end subroutine cross_Voronoi_cell

  !----------------------------------------

  subroutine cross_Voronoi_cell_vect(x,y,z, u,v,w, icell, previous_cell, x1,y1,z1, next_cell, s, s_contrib, s_void_before)

    integer, intent(in) :: icell, previous_cell
    real(kind=dp), intent(in) :: x,y,z, u,v,w

    real(kind=dp), intent(out) :: x1, y1, z1, s, s_contrib, s_void_before
    integer, intent(out) :: next_cell

    integer :: i, id_n, l


    integer, dimension(100) :: tab_id, tab_id_wall
    real, dimension(3,100) :: n, p_r
    real, dimension(100) :: den, num, s_tmp

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real, dimension(3) :: r, k, r_cell, r_neighbour

    real, dimension(3) :: n1
    real(kind=dp), dimension(3) :: delta_r
    real :: s_tmp_wall, den1
    real(kind=dp) :: b, c, delta, rac, s1, s2, h

    integer :: n_neighbours, n_wall, ifirst, ilast, id_min, id_max

    logical :: lhit_wall, was_cut


    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    ! We do all the access to Voronoi(icell) now
    r_cell = Voronoi(icell)%xyz(:)
    ifirst = Voronoi(icell)%first_neighbour
    ilast = Voronoi(icell)%last_neighbour
    was_cut = Voronoi(icell)%was_cut
    h = Voronoi(icell)%h

    s = 1e30 !huge_real
    next_cell = 0

    ! Counting number of neighbouring cells and walls
    n_neighbours = 0
    n_walls=0
    lhit_wall = .false.

    l=0
    do i=ifirst,ilast
       l = l+1
       id_n = neighbours_list(i) ! id du voisin : values are all over the place : many cache misses

       if (id_n==previous_cell) cycle

       if (id_n > 0) then ! cellule
          if (l <= n_saved_neighbours) then ! we used an ordered array to limit cache misses
             r_neighbour(:) = Voronoi_neighbour_xyz(:,l,icell)
          else
             r_neighbour(:) = Voronoi_xyz(:,id_n)
          endif

          n1(:) = r_neighbour(:) - r_cell(:) ! 1 normal vector
          den1 = dot_product(n1, k)
          if (den1 > 0) then ! we keep that cell and start filling some arrays
             n_neighbours = n_neighbours + 1
             tab_id(n_neighbours) = id_n

             n(:,n_neighbours) = n1(:)
             den(n_neighbours) = den1

             p_r(:,n_neighbours) = 0.5 * (r_neighbour(:) + r_cell(:)) - r(:) ! we can save an operation here
          endif
       else
          lhit_wall = .true.
          n_walls = n_walls + 1
          tab_id_wall(n_walls) = id_n
       endif
    enddo

    ! Allocate arrays
    !allocate(n(3,n_neighbours), p_r(3,n_neighbours), den(n_neighbours), num(n_neighbours), s_tmp(n_neighbours))

    ! Array of vectors
  !  do i=1,n_neighbours
  !     id_n = tab_id(i)
  !     !n(i,:) = Voronoi(id_n)%xyz(:) - Voronoi(icell)%xyz(:)
  !     p_r(i,:) = 0.5 * (Voronoi(id_n)%xyz(:) + r_cell(:)) - r(:) ! we can save an operation here
  !  enddo

    num(1:n_neighbours) = n(1,1:n_neighbours) * p_r(1,1:n_neighbours) + &
         n(2,1:n_neighbours) * p_r(2,1:n_neighbours) + n(3,1:n_neighbours) * p_r(3,1:n_neighbours)
    s_tmp(1:n_neighbours) = num(1:n_neighbours)/den(1:n_neighbours)

    s = huge(1.0)
    do i=1,n_neighbours
       if (s_tmp(i) > 0) then
          if (s_tmp(i) < s) then
             s = s_tmp(i)
             next_cell = tab_id(i)
          endif
       endif
    enddo

    ! We test for walls now
    if (lhit_wall) then
       do i=1, n_walls
          id_n = tab_id_wall(i)

          s_tmp_wall = distance_to_wall(x,y,z, u,v,w, -id_n) ;

          ! si c'est le wall d'entree : peut-etre a faire sauter en sauvegardant le wall d'entree
          if (s_tmp_wall < 0.) s_tmp_wall = huge(1.0)

          if (s_tmp_wall < s) then
             s = s_tmp_wall
             next_cell = id_n
          endif
       enddo
    endif

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

    if (was_cut) then
       ! We check where the packet intersect the sphere of radius Voronoi(icell)%h * PS%cutting_distance_o_h
       ! centered on the center of the cell
       delta_r = r(:) - r_cell(:)
       b = dot_product(delta_r,k)
       c = dot_product(delta_r,delta_r) - (h * PS%cutting_distance_o_h)**2
       delta = b*b - c

       if (delta < 0.) then ! the packet never encounters the sphere
          s_void_before = s
          s_contrib = 0.0_dp
       else ! the packet encounters the sphere
          rac = sqrt(delta)
          s1 = -b - rac
          s2 = -b + rac

          if (s1 < 0) then ! we already entered the sphere
             if (s2 < 0) then ! we already exited the sphere
                s_void_before = s
                s_contrib = 0.0_dp
             else ! We are still in the sphere and will exit it
                s_void_before = 0.0_dp
                s_contrib = min(s2,s)
             endif
          else ! We will enter in the sphere (both s1 and s2 are > 0)
             if (s1 < s) then ! We will enter the sphere in this cell
                s_void_before = s1
                s_contrib = min(s2,s) - s1
             else ! We will not enter the sphere in this sphere
                s_void_before = s
                s_contrib = 0.0_dp
             endif
          endif
       endif ! delta < 0
    else ! the cell was not cut
       s_void_before = 0.0_dp
       s_contrib = s
    endif

    return

  end subroutine cross_Voronoi_cell_vect

  !----------------------------------------

  subroutine init_Voronoi_walls(n_walls, limits)

    integer, intent(in) :: n_walls
    real(kind=dp), dimension(n_walls), intent(in) :: limits
    integer :: iwall, alloc_status

    if (n_walls /= 6) call error("n_walls must be 6 in Voronoi grid")

    allocate(wall(n_walls),stat=alloc_status)
    if (alloc_status /= 0) call error("Allocation error Voronoi wall")

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


  real(kind=dp) function distance_to_wall(x,y,z, u,v,w, iwall)
    ! Mur plan pour le moment : meme algorithme que cross Voronoi cell

    real(kind=dp), intent(in) :: x,y,z, u,v,w
    integer, intent(in) :: iwall

    ! n = normale a la face, p = point sur la face, r = position du photon, k = direction de vol
    real(kind=dp), dimension(3) :: n, p, r, k

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

  function distance_to_star(x,y,z,u,v,w,i_star)
  ! This routine implies that a star is in a unique cell
  ! return the index of a star (i_star > 0) if a star will be intersected
  ! by a ray going in the direction (u,v,w).
  ! Computes also the distance between the point (x,y,z) and the star in that direction.

    real(kind=dp) :: distance_to_star
    real(kind=dp), intent(in) :: x,y,z, u,v,w
    integer, intent(out) :: i_star

    real(kind=dp), dimension(3) :: r, k, delta_r
    real(kind=dp) :: b,c, delta, rac, s1, s2
    integer :: i


    r(1) = x ; r(2) = y ; r(3) = z
    k(1) = u ; k(2) = v ; k(3) = w

    distance_to_star = huge(1.0_dp)

    i_star = 0
    star_loop : do i = 1, n_etoiles
     	delta_r(:)  = r(:) - (/etoile(i)%x, etoile(i)%y, etoile(i)%z/)
     	b = dot_product(delta_r,k)
     	c = dot_product(delta_r,delta_r) - (etoile(i)%r)**2
     	delta = b*b - c

     	if (delta >= 0.) then ! the packet will encounter (or has encoutered) the star
        	rac = sqrt(delta)
        	s1 = -b - rac

        	if (s1 < 0) then ! we already entered the star
           ! We can probably skip that test, s1 must be positive as we must be outside the star
           	s2 = -b + rac
           	if (s2 > 0) then ! for s2 < 0: we already exited the star
              ! We are still in the sphere and will exit it
              ! This means that we had a round-off error somewhere
              	distance_to_star = 0.0_dp
              	i_star = i
           	endif

        else ! We will enter in the star
           	if (s1 < distance_to_star) then
              	distance_to_star = s1
              	i_star = i
           		endif
        	endif ! s1 < 0

     	endif ! delta < 0

    enddo star_loop


    return
  end function distance_to_star

  !----------------------------------------

  subroutine move_to_grid_Voronoi(id, x,y,z, u,v,w, icell, lintersect)

    integer, intent(in) :: id
    real(kind=dp), intent(inout) :: x,y,z
    real(kind=dp), intent(in) :: u,v,w

    integer, intent(out) :: icell
    logical, intent(out) :: lintersect

    logical, dimension(n_walls) :: intersect
    real(kind=dp), dimension(n_walls) :: s_walls
    integer, dimension(n_walls) :: order

    real(kind=dp) :: s, l, x_test, y_test, z_test
    integer :: i, iwall

    ! Find out which plane we are crossing first
    ! and test boundaries of the plane
    s_walls(:) = huge(1.0)
    intersect(:) = .false.
    do iwall=1, n_walls
       l = distance_to_wall(x,y,z, u,v,w, iwall) ! signe - car on rentre dans le volume

       if (l >= 0) then
          intersect(iwall) = .true.
          s_walls(iwall) = l * (1.0_dp + prec)
       else
          s_walls(iwall) = huge(1.0)
       endif
    enddo

    order = index_quicksort(real(s_walls,kind=dp))

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

pure logical function test_exit_grid_Voronoi(icell, x,y,z)

  integer, intent(in) :: icell
  real(kind=dp), intent(in) :: x,y,z

  if (icell < 0) then
     test_exit_grid_Voronoi = .true.
  else
     test_exit_grid_Voronoi = .false. ! .not.Voronoi(next_cell)%is_star
  endif

  return

end function test_exit_grid_Voronoi

!----------------------------------------

pure logical function is_in_volume(x,y,z)

  real(kind=dp), intent(in) :: x,y,z

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
  real(kind=dp), intent(in) :: x, y, z

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
  real(kind=dp), intent(out) :: x,y,z

  real(kind=dp) :: u, v, w, srw2, argmt, x1,y1,z1
  integer :: previous_cell, next_cell
  real(kind=dp) :: l, l_contrib, l_void_before

  ! Direction aleatoire
  w = 2.0_dp * aleat1 - 1.0_dp
  srw2 = sqrt(1.0_dp-w*w)
  argmt = pi*(2.0_dp*aleat2-1.0_dp)
  u = srw2 * cos(argmt)
  v = srw2 * sin(argmt)

  ! Distance jusqu'au bord de la cellule
  previous_cell = 0
  x = Voronoi(icell)%xyz(1) ; y = Voronoi(icell)%xyz(2) ; z = Voronoi(icell)%xyz(3)
  call cross_Voronoi_cell(x,y,z, u,v,w, icell, previous_cell, x1,y1,z1, next_cell, l, l_contrib, l_void_before)

  ! Repartition uniforme selon cette direction
  l = l_contrib * aleat3**(1./3)
  !x=x+l ; y=y+l*v ; z=z+l*w ! emission au centre cellule si cette ligne est commentee

  return

end subroutine pos_em_cellule_Voronoi


!----------------------------------------

subroutine indice_cellule_Voronoi(xin,yin,zin, icell)

    implicit none

    real(kind=dp), intent(in) :: xin,yin,zin
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
  if (allocated(Voronoi_xyz)) deallocate(Voronoi_xyz)
  if (allocated(Voronoi_neighbour_xyz)) deallocate(Voronoi_neighbour_xyz)
  if (allocated(volume)) deallocate(volume)
  if (allocated(neighbours_list)) deallocate(neighbours_list)
  if (allocated(wall)) deallocate(wall)

  if (allocated(wall_cells)) deallocate(wall_cells)
  call deallocate_kdtree2_search()

  return

end subroutine deallocate_Voronoi

!----------------------------------------

subroutine build_wall_kdtrees()

  integer :: iwall, i, icell, alloc_status

  call allocate_kdtree2_search(nb_proc)

  allocate(wall_cells(3,maxval(wall(:)%n_neighbours),n_walls), stat=alloc_status)
  if (alloc_status /=0) call error("Allocation error wall kdtrees")
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
  real(kind=dp), intent(in) :: x, y, z

  real(kind=dp), dimension(3), target :: qv
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
!    real(kind=dp), allocatable :: mydata(:,:)
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
