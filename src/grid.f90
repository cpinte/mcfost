module grid

  use parametres
  use constantes
  use grains
  use utils
  use sort, only : index_quicksort
  use mcfost_env, only : dp

  use cylindrical_grid
  use spherical_grid
  use Voronoi_grid

  implicit none

  procedure(cross_cylindrical_cell), pointer :: cross_cell => null()
  procedure(pos_em_cellule_cyl), pointer :: pos_em_cellule => null()
  procedure(move_to_grid_cyl), pointer :: move_to_grid => null()
  procedure(indice_cellule_cyl), pointer :: indice_cellule => null()
  procedure(test_exit_grid_cyl), pointer :: test_exit_grid => null()
  procedure(define_cylindrical_grid), pointer :: define_grid => null()

  real(kind=dp) :: v_char, B_char
  logical :: lcalc_ne, lmagnetized
  !real(kind=dp), dimension(:), allocatable :: vfield_v1, vfield_v2, vfield_v3
  real(kind=dp), dimension(:), allocatable :: ne, nHtot, T, nHmin, vturb
  integer, dimension(:), allocatable :: icompute_atomRT
  real, dimension(:,:), allocatable :: vfield3d ! n_cells x 3
  real, dimension(:), allocatable :: vfield ! n_cells

  contains

  subroutine alloc_atomrt_grid()
   integer(kind=8) :: mem_alloc_local = 0


   allocate(nHtot(n_cells), ne(n_cells), T(n_cells))
   allocate(nHmin(n_cells), vturb(n_cells))

   !here vfield_x,_y and _z are either R,z, phi or r,theta, phi.
   !works onlt if lvelocity_file = .false.
   if (allocated(vfield_x)) deallocate(vfield_x, vfield_y, vfield_z)
   if (.not.lvoronoi) then
   !    allocate(vfield_v1(n_cells)); vfield_v1 = 0.0_dp
   !    allocate(vfield_v2(n_cells)); vfield_v2 = 0.0_dp
   !    allocate(vfield_v3(n_cells)); vfield_v3 = 0.0_dp
       allocate(vfield_x(n_cells)); vfield_x = 0.0_dp
       allocate(vfield_y(n_cells)); vfield_y = 0.0_dp
       allocate(vfield_z(n_cells)); vfield_z = 0.0_dp
   end if

   lcalc_ne = .false.
   T = 0.0
   nHtot = 0.0
   ne = 0.0
   nHmin = 0.0
   vturb = 0.0 !m/s
   

   allocate(icompute_atomRT(n_cells))
   icompute_atomRT = 0 !everything transparent at init.


   ! mem_alloc_tot = mem_alloc_tot + mem_alloc_local
   ! write(*,'("Total memory allocated in alloc_atomic_atmos:"(1ES17.8E3)" MB")') mem_alloc_local / 1024./1024.


   return
 end subroutine alloc_atomrt_grid

 subroutine dealloc_atomrt_grid

   deallocate(nHtot, ne, T, vturb)
   if (allocated(nHmin)) deallocate(nHmin)

   if (.not.lvoronoi) then
      deallocate(vfield_x,vfield_y,vfield_z)
   end if


   deallocate(icompute_atomRT)

   return
 end subroutine dealloc_atomrt_grid

 subroutine order_zones()
   ! Order the various zones according to their Rin
   ! C. Pinte
   ! 04/05/11

   integer, dimension(n_zones) :: order
   type(disk_zone_type), dimension(n_zones) :: disk_zone_tmp
   integer, dimension(n_pop) :: Izone_tmp

   integer :: i, ipop

   ! Save arrays to order
   disk_zone_tmp(:) = disk_zone(:)
   Izone_tmp(:) =  dust_pop(:)%zone

   ! order following Rin
   order = index_quicksort(disk_zone(:)%Rmin)

   ! Reordering zones
   do i=1, n_zones
      disk_zone(i) =  disk_zone_tmp(order(i))
      do ipop=1,n_pop  ! reordering zone index in pops
         if (Izone_tmp(ipop) == order(i))  dust_pop(ipop)%zone = i
      enddo
   enddo

   return

 end subroutine order_zones

!******************************************************************************

subroutine define_physical_zones()
  ! Recheche les connections de zone 2 a 2
  ! on doit pouvoir faire mieux que ca
  ! C. Pinte
  ! 03/05/11

  integer :: i, j, index, i_region, iter, ir, k
  logical, dimension(n_zones) :: zone_scanned
  real(kind=dp) :: r1, r2, minR, maxR
  character(len=10) :: n1, n2

  logical :: test_j_in_i, test_i_in_j

  ! Detecting connected zones
  zone_scanned(:) = .false.
  index = 0
  do i=1, n_zones
     if (.not.zone_scanned(i)) then
        index = index + 1
        disk_zone(i)%region = index
        zone_scanned(i) = .true.

        ! Current minimum & maximum radii of region
        minR = disk_zone(i)%Rmin
        maxR = disk_zone(i)%Rmax

        ! Besoin d'iterer au cas ou les connections entre zones sont multiples
        ! probablement pas autant mais ca ne coute rien en calcul
        do iter=1, n_zones-1
           do j=i+1, n_zones

              r1 = disk_zone(j)%Rmin
              r2 = disk_zone(j)%Rmax

              ! Test if the 2 zones are imbrigated
              test_j_in_i = ((r1 > minR).and.(r1 < maxR)) .or. ((r2 > minR).and.(r2 <= maxR))
              test_i_in_j = ((minR > r1).and.(minR < r2)) .or. ((maxR > r1).and.(maxR <= r2))

              if ( test_j_in_i .or. test_i_in_j ) then
                 if (.not.zone_scanned(j)) then
                    i_region = index
                 else
                    i_region = disk_zone(j)%region
                 endif ! zone_scanned

                 disk_zone(j)%region = i_region
                 zone_scanned(j) = .true.

                 ! Updating minimum & maximum radii of region
                 minR = min(minR,r1)
                 maxR = max(maxR,r2)
              endif ! test rayon

           enddo ! j
        enddo ! iter
     endif !.not.zone_scanned(i)
  enddo !i

  n_regions = maxval(disk_zone(:)%region)

  allocate(regions(n_regions))
  do ir=1,n_regions
     k = 0
     do i=1, n_zones
        if (disk_zone(i)%region == ir)  k=k+1
     enddo ! i
     regions(ir)%n_zones = k
     allocate(regions(ir)%zones(regions(ir)%n_zones))

     k = 0
     do i=1, n_zones
        if (disk_zone(i)%region == ir)  then
           k=k+1
           regions(ir)%zones(k) = i
        endif
     enddo ! i
  enddo ! ir


  do ir = 1, n_regions
     regions(ir)%Rmin = 1e30
     regions(ir)%Rmax = 0
     do i=1, n_zones
        if (disk_zone(i)%region == ir) then
           regions(ir)%Rmin = min(regions(ir)%Rmin,disk_zone(i)%Rmin)
           regions(ir)%Rmax = max(regions(ir)%Rmax,disk_zone(i)%Rmax)
        endif
     enddo !i
  enddo !ir

  ! Adjusting grid for prolate or oblate envelope
  if (z_scaling_env > 1.) then
     regions(:)%Rmax = regions(:)%Rmax * z_scaling_env
  endif
  if (z_scaling_env < 1.) then
     regions(:)%Rmin = regions(:)%Rmin * z_scaling_env
  endif

  Rmin = minval(regions(:)%Rmin)
  Rmax = maxval(regions(:)%Rmax)

  if (Rmin < 0.0) call error("R_min < 0.0")
  do i=1, n_etoiles
     if ( (abs(etoile(i)%x) < tiny_real).and.(abs(etoile(i)%x) < tiny_real).and.(abs(etoile(i)%x) < tiny_real) ) then
        if (etoile(i)%r > Rmin) then
           write(*,*) "Rstar =", etoile(i)%r, "Rmin=", Rmin
           call error("inner disk radius is smaller than stellar radius")
        endif
     endif
  enddo

  write(*,fmt='(" Number of regions detected:",i2)') n_regions

  do i=1, n_zones
     R1 = real(regions(disk_zone(i)%region)%Rmin)
     R2 = real(regions(disk_zone(i)%region)%Rmax)
     ! Format
     if ((R1 <= 1e-2).or.(R1>=1e6)) then
        n1 = "es8.2"
     else
        n1 = "f"//achar(int(abs(log10(R1))+1)+iachar('3'))//".2"
     endif
     if ((R2 <= 1e-2).or.(R2>=1e6)) then
        n2 = "es8.2"
     else
        n2 = "f"//achar(int(abs(log10(R2))+1)+iachar('3'))//".2"
     endif
     write(*,fmt='(" zone",i2," --> region=",i2," : R=",'//trim(n1)//'," to ",'//trim(n2)//'," AU")') &
          i, disk_zone(i)%region, R1, R2
  enddo

  return

end subroutine define_physical_zones

!******************************************************************************

subroutine setup_grid()

  logical, save :: lfirst = .true.
  integer :: mem_size

  if (.not.lVoronoi) then
     nrz = n_rad * nz
     if (l3D) then
        n_cells = 2*nrz*n_az
     else
        n_cells = nrz
     endif
  endif

  if (n_cells < 1e6) then
     write(*,*) "Using", n_cells, "cells"
  else
     write(*,*) "Using", real(n_cells)/1e6, "million cells"
  endif

  if (lvariable_dust) then
     p_n_cells = n_cells
  else
     p_n_cells = 1
  endif

  if (lVoronoi) then
     write(*,*) "Using a Voronoi mesh"
     lcylindrical = .false.
     lspherical = .false.
     cross_cell => cross_Voronoi_cell
     !cross_cell => cross_Voronoi_cell_vect
     pos_em_cellule => pos_em_cellule_Voronoi
     move_to_grid => move_to_grid_Voronoi
     indice_cellule => indice_cellule_Voronoi
     test_exit_grid => test_exit_grid_Voronoi
     define_grid => define_Voronoi_grid
  else
     if (lvariable_dust) then
        p_n_rad=n_rad ; p_nz = nz
     else
        p_n_rad=1 ;  p_nz=1
     endif

     if (l3D) then
        j_start = -nz
        if (lvariable_dust) then
           p_n_az = n_az
        else
           p_n_az = 1
        endif
     else
        j_start = 1
        p_n_az = 1
     endif

     if ((p_nz /= 1).and.l3D) then
        pj_start = -nz
     else
        pj_start = 1
     endif

     if (grid_type == 1) then
        lcylindrical = .true.
        lspherical = .false.
        cross_cell => cross_cylindrical_cell
        pos_em_cellule => pos_em_cellule_cyl
        move_to_grid => move_to_grid_cyl
        indice_cellule => indice_cellule_cyl
        test_exit_grid => test_exit_grid_cyl
        define_grid => define_cylindrical_grid
     else if (grid_type == 2) then
        lcylindrical = .false.
        lspherical = .true.
        cross_cell => cross_spherical_cell
        pos_em_cellule => pos_em_cellule_sph
        move_to_grid => move_to_grid_sph
        indice_cellule => indice_cellule_sph
        test_exit_grid => test_exit_grid_sph
        define_grid => define_cylindrical_grid ! same routine at the moment
     else
        call error("Unknown grid type")
     endif

     if (lfirst) then
        call build_cylindrical_cell_mapping() ! same for spherical
        lfirst = .false.
     endif
  endif

  return

end subroutine setup_grid

!to do
!subroutine to find the neighbours of a cell.

subroutine nnk_smooth(arr)
	!Nearest Neighbours Kernel smoother for spherical and Voronoi grids.
	!TO DO:
	! 	check Voronoi implementation
	!   weight the value of the neighbours depending of their distance to icell
  	real(kind=dp), intent(inout) :: arr(n_cells)
  	real(kind=dp) :: A(n_cells)
  	integer :: i,j,k, icell
  	integer :: n_az_end
  	real(kind=dp) :: w0, w1, w2
  	integer :: ifirst, ilast, l, id_n
  	
  	!in principle the weights can be function of the dimension
  	w0 = 1.0
  	w1 = w0 / 6.0
  	w2 = w1 / 6.0
  	
  	!does not work if grid is cylindrical yet
  	if (lcylindrical) return
  	
  	!building
  	if (lvoronoi) then
  		!for each Voronoi cell, loop over the nearest neighbours 
  		do icell=1,n_cells
			ifirst = Voronoi(icell)%first_neighbour
			ilast = Voronoi(icell)%last_neighbour
			l=0
			A(icell) = w0 * arr(icell)
			nb_loop : do i=ifirst,ilast
				l = l+1
				id_n = neighbours_list(i)
				if (id_n==icell) cycle nb_loop
				
				!is that correct with Voronoi ? id_n is the neigbour ?
				!There should also be a way of using a weight depending on the distance of the neighbour
       			if (id_n > 0) then ! cellule
					A(icell) = A(icell) + w0 * arr(id_n)
				endif !no wall
				
				
			enddo nb_loop	
			A(icell) = A(icell) / (w0 + l * w0)	
  		enddo  	
  		arr = A
  	else !spherical only at the moment

  		if (l3d) then
  			!To Do. Computes the value of the edges
    		A = arr !copy array to temporarily handle the edges (which are not smoothed)!
    		if (n_az==1) then
    			!2.5d because z < 0 but phi = 0 everywhere
    			n_az_end = 1
    		else
    			!full 3d
    			n_az_end = n_az-1
    		endif
			do i=2, n_rad-1
				do j=j_start+1,nz-1
					if (j==0) cycle
					do k=1, n_az_end
						icell = cell_map(i,j,k)
						A(icell) = ( w0 * arr(icell) &
						+ w1*arr(cell_map(i,j,k+1)) &
						+ w1*arr(cell_map(i,j,k-1)) &
						+ w1*arr(cell_map(i,j+1,k)) &
						+ w2*arr(cell_map(i,j+1,k+1)) &
						+ w2*arr(cell_map(i,j+1,k-1)) &
						+ w1*arr(cell_map(i,j-1,k)) &
						+ w2*arr(cell_map(i,j-1,k+1)) &
						+ w2*arr(cell_map(i,j-1,k-1)) &
						+ w1*arr(cell_map(i+1,j,k)) &
						+ w1*arr(cell_map(i-1,j,k)) &
						+ w2*arr(cell_map(i+1,j,k+1)) &
						+ w2*arr(cell_map(i+1,j,k-1)) &
						+ w2*arr(cell_map(i+1,j+1,k)) &
						+ w2*arr(cell_map(i+1,j+1,k+1)) &
						+ w2*arr(cell_map(i+1,j+1,k-1)) &
						+ w2*arr(cell_map(i+1,j-1,k)) &
						+ w2*arr(cell_map(i+1,j-1,k+1)) &
						+ w2*arr(cell_map(i+1,j-1,k-1)) &
						+ w2*arr(cell_map(i-1,j,k+1)) &
						+ w2*arr(cell_map(i-1,j,k-1)) &
						+ w2*arr(cell_map(i-1,j+1,k)) &
						+ w2*arr(cell_map(i-1,j+1,k+1)) &
						+ w2*arr(cell_map(i-1,j+1,k-1)) &
						+ w2*arr(cell_map(i-1,j-1,k)) &
						+ w2*arr(cell_map(i-1,j-1,k+1)) &
						+ w2*arr(cell_map(i-1,j-1,k-1)) ) / (w0 + 6*w1 + 20*w2)
					enddo
				enddo
			enddo
		else
			k = 1
			do i=2, n_rad-1
				do j=j_start+1,nz-1
					if (j==0) cycle
					icell = cell_map(i,j,k)
					A(icell) = ( w0*arr(icell) &
						+ w1*arr(cell_map(i,j+1,k)) &
						+ w1*arr(cell_map(i,j-1,k)) &
						+ w1*arr(cell_map(i+1,j,k)) &
						+ w1*arr(cell_map(i-1,j,k)) &
						+ w2*arr(cell_map(i+1,j+1,k)) &
						+ w2*arr(cell_map(i+1,j-1,k)) &
						+ w2*arr(cell_map(i-1,j+1,k)) &
						+ w2*arr(cell_map(i-1,j-1,k)) ) / (w0 + 4 * w1 + 3 * w2)
				enddo
				!for all i, the edges j points
				A(cell_map(i,j_start,k)) = ( w0*arr(cell_map(i,j_start,k)) + &
					w1*arr(cell_map(i,j_start+1,k))+w1*arr(cell_map(i+1,j_start,k))+&
					w2*arr(cell_map(i+1,j_start+1,k))+w1*arr(cell_map(i-1,j_start,k))+&
					w2*arr(cell_map(i-1,j_start+1,k)) ) / (w0 + 3*w1+2*w2)

				A(cell_map(i,nz,k)) = ( w0*arr(cell_map(i,nz,k)) + &
					w1*arr(cell_map(i,nz-1,k))+w1*arr(cell_map(i+1,nz,k))+&
					w2*arr(cell_map(i+1,nz-1,k))+w1*arr(cell_map(i-1,nz,k))+&
					w2*arr(cell_map(i-1,nz-1,k)) ) / (w0 + 3*w1+2*w2)
			enddo
			!for all j, the edges in i
			do j=j_start+1,nz-1
				if (j==0) cycle
				icell = cell_map(1,j,k)
				A(icell) = ( w0*arr(icell) &
					+ w1*arr(cell_map(1,j+1,k)) &
					+ w1*arr(cell_map(1,j-1,k)) &
					+ w1*arr(cell_map(1+1,j,k)) &
					+ w2*arr(cell_map(1+1,j+1,k)) &
					+ w2*arr(cell_map(1+1,j-1,k)) ) / (w0 + 3*w1 + 2*w2)

				icell = cell_map(n_rad,j,k)					
				A(icell) = ( w0*arr(icell) &
						+ w1*arr(cell_map(n_rad,j+1,k)) &
						+ w1*arr(cell_map(n_rad,j-1,k)) &
						+ w1*arr(cell_map(n_rad-1,j,k)) &
						+ w2*arr(cell_map(n_rad-1,j+1,k)) &
						+ w2*arr(cell_map(n_rad-1,j-1,k)) ) / (w0 + 3*w1 + 2*w2)
			enddo
			i = 1
			!for edges in i and j
			A(cell_map(i,j_start,k)) = ( w0*arr(cell_map(i,j_start,k)) + &
					w1*arr(cell_map(i,j_start+1,k))+w1*arr(cell_map(i+1,j_start,k))+&
					w2*arr(cell_map(i+1,j_start+1,k)) ) / (w0 + 2*w1 + w2)
			A(cell_map(i,nz,k)) = ( w0*arr(cell_map(i,nz,k)) + &
					w1*arr(cell_map(i,nz-1,k))+w1*arr(cell_map(i+1,nz,k))+&
					w2*arr(cell_map(i+1,nz-1,k)) ) / (w0 + 2*w1 + w2)
					
			i = n_rad
			A(cell_map(i,j_start,k)) = ( w0*arr(cell_map(i,j_start,k)) + &
					w1*arr(cell_map(i,j_start+1,k))+w1*arr(cell_map(i-1,j_start,k))+&
					w2*arr(cell_map(i-1,j_start+1,k)) ) / (w0 + 2*w1 + w2)
			A(cell_map(i,nz,k)) = ( w0*arr(cell_map(i,nz,k)) + &
					w1*arr(cell_map(i,nz-1,k))+w1*arr(cell_map(i-1,nz,k))+&
					w2*arr(cell_map(i-1,nz-1,k)) ) / (w0 + 2*w1 + w2)
		endif !2d
		arr = A !lspherical
	endif

  return
  end subroutine nnk_smooth

!******************************************************************************

end module grid
