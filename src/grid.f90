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
  real(kind=dp), dimension(:), allocatable :: ne, nHtot, T, nHmin, vturb ! n_cells
  integer, dimension(:), allocatable :: icompute_atomRT !n_cells
  real, dimension(:,:), allocatable :: vfield3d ! n_cells x 3
  real, dimension(:), allocatable :: vfield ! n_cells

  contains

  subroutine alloc_atomrt_grid()
    integer(kind=8) :: mem_alloc_local = 0

    !merge vturb and v_turb (molecular emission)
    !TO DO: move vturb in molecular emission in grid.f90
    !  use that vturb for atomRT
    allocate(nHtot(n_cells), ne(n_cells), T(n_cells))
    allocate(nHmin(n_cells), vturb(n_cells))

    !here vfield3d are either R,z, phi or r,theta, phi.
    !works onlt if lvelocity_file = .false.
    if (allocated(vfield3d)) deallocate(vfield3d)
    if (.not.lvoronoi) then
       allocate(vfield3d(n_cells,3)); vfield3d = 0.0_dp
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
   if (.not.lvoronoi) deallocate(vfield3d)
   deallocate(icompute_atomRT)

   return

 end subroutine dealloc_atomrt_grid

  subroutine check_for_zero_electronic_density()
   integer :: N_fixed_ne, icell

   lcalc_ne = .false.
   icell_loop : do icell=1,n_cells
      !check that in filled cells there is electron density otherwise we need to compute it
      !from scratch.
      if (icompute_atomRT(icell) > 0) then

         if (ne(icell) <= 0.0_dp) then
            write(*,*) "  ** No electron density found in the model! ** "
            lcalc_ne = .true.
            exit icell_loop
         endif

      endif
   enddo icell_loop
   N_fixed_ne = size(pack(icompute_atomRT,mask=(icompute_atomRT==2)))
   if (N_fixed_ne > 0) then
      write(*,'("Found "(1I5)" cells with fixed electron density values! ("(1I3)" %)")') &
           N_fixed_ne, nint(real(N_fixed_ne) / real(n_cells) * 100)
   endif

   return

 end subroutine check_for_zero_electronic_density

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

!******************************************************************************

end module grid
