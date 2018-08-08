module grid

  use parametres
  use constantes
  use opacity
  use grains
  use em_th
  use prop_star
  use utils
  use cylindrical_grid
  use spherical_grid
  use Voronoi_grid
  use ray_tracing, only : lscatt_ray_tracing

  implicit none

  procedure(cross_cylindrical_cell), pointer :: cross_cell => null()
  procedure(pos_em_cellule_cyl), pointer :: pos_em_cellule => null()
  procedure(move_to_grid_cyl), pointer :: move_to_grid => null()
  procedure(indice_cellule_cyl), pointer :: indice_cellule => null()
  procedure(test_exit_grid_cyl), pointer :: test_exit_grid => null()
  procedure(define_cylindrical_grid), pointer :: define_grid => null()

  contains

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
   order = bubble_sort(disk_zone(:)%Rmin)

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

  use dust_ray_tracing, only : select_scattering_method

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

  ! parametrage methode de diffusion en fonction de la taille de la grille
  ! 1 : per dust grain
  ! 2 : per cell
  call select_scattering_method(p_n_cells)

  lMueller_pos_multi = .false.
  if (lmono) then
     p_n_lambda_pos = 1
  else
     if (scattering_method==1) then
        p_n_lambda_pos = 1
     else
        p_n_lambda_pos = n_lambda
        lMueller_pos_multi = .true.
     endif
  endif

  return

end subroutine setup_grid

!******************************************************************************

subroutine init_lambda()
  ! Initialisation table de longueurs d'onde

  integer :: i, alloc_status

  allocate(tab_lambda(n_lambda), tab_lambda_inf(n_lambda), tab_lambda_sup(n_lambda), tab_delta_lambda(n_lambda), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_lambda')
  tab_lambda=0.0 ;  tab_lambda_inf=0.0 ; tab_lambda_sup=0.0 ; tab_delta_lambda=0.0

  allocate(tab_amu1(n_lambda, n_pop), tab_amu2(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_amu')
  tab_amu1=0.0; tab_amu2=0.0;

  allocate(tab_amu1_coating(n_lambda, n_pop), tab_amu2_coating(n_lambda, n_pop), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_amu_coating')
  tab_amu1_coating=0.0; tab_amu2_coating=0.0;

  if (lmono0) then
     ! Lecture longueur d'onde
     read(band,*) tab_lambda(1)
     tab_delta_lambda(1) = 1.0
     tab_lambda_inf(1) = tab_lambda(1)
     tab_lambda_sup(1) = tab_lambda(1)

  else
     ! Initialisation longueurs d'onde
     !delta_lambda = (lambda_max/lambda_min)**(1.0/real(n_lambda))
     delta_lambda =  exp( (1.0_dp/real(n_lambda,kind=dp)) * log(lambda_max/lambda_min) )

     tab_lambda_inf(1) = lambda_min
     tab_lambda(1)=lambda_min*sqrt(delta_lambda)
     tab_lambda_sup(1) = lambda_min*delta_lambda
     do i=2, n_lambda
        tab_lambda(i)= tab_lambda(i-1)*delta_lambda
        tab_lambda_sup(i)= tab_lambda_sup(i-1)*delta_lambda
        tab_lambda_inf(i)= tab_lambda_sup(i-1)
     enddo

     do i=1, n_lambda
        tab_delta_lambda(i) = tab_lambda_sup(i) - tab_lambda_inf(i)
     enddo

  endif

end subroutine init_lambda

!**********************************************************************

subroutine init_lambda2()
  ! Initialisation table en lambda sed

  implicit none

  integer :: i

  n_lambda=n_lambda2
  do i=1, n_lambda2
     tab_lambda(i)= tab_lambda2(i)
     tab_lambda_inf(i)= tab_lambda2_inf(i)
     tab_lambda_sup(i)= tab_lambda2_sup(i)
     tab_delta_lambda(i)= tab_delta_lambda2(i)
  enddo

  return

end subroutine init_lambda2

!**********************************************************************


end module grid
