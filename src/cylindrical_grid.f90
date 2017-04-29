module cylindrical_grid

  use mem

  implicit none

contains

  subroutine build_cylindrical_cell_mapping() ! work also in sph

    integer :: i,j,k,icell, ntot, ntot2, alloc_status
    integer :: istart,iend,jstart,jend,kstart,kend, istart2,iend2,jstart2,jend2,kstart2,kend2

    icell_ref = 1

    istart = 1
    iend = n_rad

    jstart = j_start
    jend = nz

    kstart=1
    kend = n_az

    if (j_start < 0) then
       ntot = (iend - istart + 1) * (jend - jstart) * (kend - kstart + 1)
    else
       ntot = (iend - istart + 1) * (jend - jstart +1) * (kend - kstart + 1)
    endif

    if (ntot /= n_cells) then
       write(*,*) "ERROR in 'build_cylindrical_cell_mapping'"
       write(*,*) "The number of cells is not matching :"
       write(*,*) "ntot=", ntot, "should be", n_cells
       write(*,*) "Exiting."
       stop
    endif

    istart2 = 0
    iend2 = n_rad + 1

    jstart2 = min(1,j_start)-1
    jend2 = nz+1

    kstart2=1
    kend2 = n_az

    if (jstart2 < 0) then
       ntot2 = (iend2 - istart2 + 1) * (jend2 - jstart2) * (kend2 - kstart2 + 1)
    else
       ntot2 = (iend2 - istart2 + 1) * (jend2 - jstart2 +1) * (kend2 - kstart2 + 1)
    endif
    allocate(cell_map(istart2:iend2,jstart2:jend2,kstart2:kend2))
    allocate(cell_map_i(ntot2), cell_map_j(ntot2), cell_map_k(ntot2))

    ! Actual cells
    icell = 0
    do k=kstart, kend
       bz : do j=j_start, jend
          if (j==0) cycle bz
          do i=istart, iend

             icell = icell+1
             if (icell > ntot) then
                write(*,*) "ERROR : there is an issue in the cell mapping"
                write(*,*) "Exiting"
                stop
             endif

             cell_map_i(icell) = i
             cell_map_j(icell) = j
             cell_map_k(icell) = k

             cell_map(i,j,k) = icell
          enddo
       enddo bz
    enddo

    if (icell /= ntot) then
       write(*,*) "Something went wrong in the call mapping"
       write(*,*) "I am missing some real cells"
       write(*,*) icell, ntot
       write(*,*)
       stop
    endif


    ! Virtual cell indices for when the packets are just around the grid

    ! Can the packet exit from this cell : 0 -> no, 1 -> radially, 2 -> vertically
    ! Warning lexit_cells == 2 works only in cylindrical
    allocate(lexit_cell(1:ntot2), stat=alloc_status)
    if (alloc_status > 0) then
       write(*,*) 'Allocation error lexit_cell'
       stop
    endif
    lexit_cell(:) = 0

    ! Cases j=0 and j=nz+1
    do k=kstart, kend
       do j = jstart2, jend2, jend2 - jstart2
          do i=istart2, iend2

             icell = icell+1
             if (icell > ntot2) then
                write(*,*) "ERROR : there is an issue in the cell mapping"
                write(*,*) "Exiting"
                stop
             endif

             if (abs(j)==jend2) lexit_cell(icell) = 2
             if (i==iend2) lexit_cell(icell) = 1

             cell_map_i(icell) = i
             cell_map_j(icell) = j
             cell_map_k(icell) = k

             cell_map(i,j,k) = icell
          enddo
       enddo
    enddo

    ! Cases i=0 and i=n_rad+1 (except j=0 and j=nz+1 done above)
    do k=kstart, kend
       bz2 : do j = jstart, jend
          if (j==0) cycle bz2
          do i=istart2,iend2, iend2-istart2

             icell = icell+1
             if (icell > ntot2) then
                write(*,*) "ERROR : there is an issue in the cell mapping"
                write(*,*) "Extra cells:", icell, ntot2
                write(*,*) i,j,k
                write(*,*) "Exiting"
                stop
             endif

             if (i==iend2) lexit_cell(icell) = 1

             cell_map_i(icell) = i
             cell_map_j(icell) = j
             cell_map_k(icell) = k

             cell_map(i,j,k) = icell
          enddo
       enddo bz2
    enddo

    if (icell /= ntot2) then
       write(*,*) "Something went wrong in the cell mapping"
       write(*,*) "I am missing some virtual cells"
       write(*,*) icell, ntot2
       write(*,*)
       stop
    endif

    !if (cell_map(1,1,1) /= 1) then
    !   write(*,*) "WARNING : mapping of cell (1,1,1) is not 1"
    !   write(*,*) "(1,1,1) --->", cell_map(1,1,1)
    !   write(*,*) "MCFOST might crash"
    !   !write(*,*) "Exiting"
    !   !stop
    !endif

    return

  end subroutine build_cylindrical_cell_mapping

!******************************************************************************

subroutine define_cylindrical_grid()
  ! Definit la grille du code
  ! Calcule les tableaux zmax, volume, r_lim, r_lim_2, z_lim
  ! et la variable Rmax2
  ! Version 4 gere les subdivisions pour les zones multiples
  ! C. Pinte
  ! 03/05/11, version 3 :  27/04/05

  real, parameter :: pi = 3.1415926535
  real(kind=dp) :: rcyl, puiss, rsph, w, uv, p, rcyl_min, rcyl_max, frac
  real :: phi
  integer :: i,j,k, izone, ii, ii_min, ii_max, icell

  !tab en cylindrique ou spherique suivant grille
  real(kind=dp), dimension(n_rad) :: V
  real(kind=dp), dimension(n_rad+1) :: tab_r, tab_r2, tab_r3
  real(kind=dp) ::   r_i, r_f, dr, fac, r0, H, hzone
  real(kind=dp) :: delta_r, ln_delta_r, delta_r_in, ln_delta_r_in
  integer :: ir, iz, n_cells_tmp, n_rad_region, n_rad_in_region, n_empty, istart, alloc_status

  type(disk_zone_type) :: dz

  real(kind=dp), dimension(:,:), allocatable :: r_grid_tmp, z_grid_tmp
  real(kind=dp), dimension(:), allocatable :: phi_grid_tmp

  if (l3D) then
     allocate(r_grid_tmp(n_rad,-nz:nz), z_grid_tmp(n_rad,-nz:nz), phi_grid_tmp(n_az), stat=alloc_status)
  else
     allocate(r_grid_tmp(n_rad,nz), z_grid_tmp(n_rad,nz), phi_grid_tmp(n_az), stat=alloc_status)
  endif


  ! **************************************************
  ! Tableaux relatifs a la grille
  ! **************************************************
  allocate(r_lim(0:n_rad), r_lim_2(0:n_rad), r_lim_3(0:n_rad), &
  delta_z(n_rad), dr2_grid(n_rad), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_lim'
     stop
  endif
  r_lim = 0.0 ; r_lim_2=0.0; r_lim_3=0.0 ; delta_z=0.0 ; dr2_grid=0.0

  allocate(r_grid(n_cells), z_grid(n_cells), phi_grid(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_lim'
     stop
  endif
  r_grid=0.0; z_grid=0.0 ; phi_grid = 0.0

  allocate(tab_region(n_rad), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_region'
     stop
  endif
  tab_region = 0

  allocate(z_lim(n_rad,nz+2), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error z_lim'
     stop
  endif
  z_lim = 0.0

  allocate(w_lim(0:nz),  theta_lim(0:nz),tan_theta_lim(0:nz),tan_phi_lim(n_az), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tan_phi_lim'
     stop
  endif
  w_lim = 0.0
  theta_lim=0.0
  tan_theta_lim = 0.0
  tan_phi_lim = 0.0

  allocate(zmax(n_rad),volume(n_cells), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error zmax, volume'
     stop
  endif
  zmax = 0.0 ; volume=0.0
  ! end allocation


  Rmax2 = Rmax*Rmax

  n_rad_in = max(n_rad_in,1) ! in case n_rad_in is set to 0 by user

  ! Definition du nombre de chaques cellules
  n_empty = 3
  n_rad_region = (n_rad - (n_regions -1) * n_empty) / n_regions
  n_rad_in_region = n_rad_in

  n_cells_tmp = 0

  istart = 1
  tab_r(:) = 0.0_dp
  do ir=1, n_regions
     regions(ir)%iRmin = istart ; regions(ir)%iRmax = istart+n_rad_region-1 ;

     if (ir == n_regions) then
        n_rad_region = n_rad - n_cells_tmp ! On prend toutes les celles restantes
     endif

     ! Pour eviter d'avoir 2 cellules a la meme position si les regions se touchent
     R0 =  regions(ir)%Rmin
     if (ir > 1) then
        if (regions(ir)%Rmin == regions(ir-1)%Rmax) then
           R0 =  regions(ir)%Rmin * 1.00001_dp
        endif
     endif

     ! Grille log avec subdivision cellule interne
     !delta_r = (rout/rmin)**(1.0/(real(n_rad-n_rad_in+1)))
     ln_delta_r = (1.0_dp/real(n_rad_region-n_rad_in_region+1,kind=dp))*log(regions(ir)%Rmax/R0)
     delta_r = exp(ln_delta_r)

     ln_delta_r_in = (1.0_dp/real(n_rad_in_region,kind=dp))*log(delta_r)
     delta_r_in = exp(ln_delta_r_in)

     ! Selection de la zone correpondante : pente la plus forte
     puiss = 0.0_dp
     do iz=1, n_zones
        if (disk_zone(iz)%region == ir) then
           p=1+dz%surf-dz%exp_beta
           if (p > puiss) then
              puiss = p
           endif
        endif
     enddo

     ! Calcul recursif hors boucle //
     ! Calcul les rayons separant les cellules de (1 a n_rad + 1)
     tab_r(istart) = R0
     tab_r2(istart) = tab_r(istart) * tab_r(istart)
     tab_r3(istart) = tab_r2(istart) * tab_r(istart)

     if (puiss == 0.0) then
        do i=istart+1, istart + n_rad_in_region
           tab_r(i) = exp(log(R0) - (log(R0)-log(R0*delta_r))*(2.0**(i-istart)-1.0)/(2.0**n_rad_in_region-1.0))
           tab_r2(i) = tab_r(i) * tab_r(i)
           tab_r3(i) = tab_r2(i) * tab_r(i)
        enddo
     else
        r_i = exp(puiss*log(R0))
        r_f = exp(puiss*log(R0*delta_r))
        dr=r_f-r_i
        fac = 1.0/(2.0**(n_rad_in_region+1)-1.0)
        do i=istart+1, istart + n_rad_in_region
           tab_r(i) = (R0**puiss - (R0**puiss-(R0*delta_r)**puiss) &
                *(2.0**(i-istart+1)-1.0)/(2.0**(n_rad_in_region+1)-1.0))**(1.0/puiss)
           !     tab_rcyl(i) = exp( 1.0/puiss * log(r_i + dr * (2.0**(i)-1.0) * fac) )
           !if (tab_rcyl(i) - tab_rcyl(i-1) < 1.0d-15*tab_rcyl(i-1)) then
           if (tab_r(i) - tab_r(i-1) < prec_grille*tab_r(i-1)) then
              write(*,*) "Error : spatial grid resolution too high"
              write(*,*) "Differences between two cells are below double precision"
              stop
           endif
           tab_r2(i) = tab_r(i) * tab_r(i)
           tab_r3(i) = tab_r2(i) * tab_r(i)
        enddo
     endif

     ! Grille log apres subdivision "1ere" cellule
     do i=istart + n_rad_in_region+1, istart+n_rad_region
        tab_r(i) = tab_r(i-1) * delta_r
        tab_r2(i) = tab_r(i) * tab_r(i)
        tab_r3(i) = tab_r2(i) * tab_r(i)
     enddo

     n_cells_tmp = istart+n_rad_region

     ! Cellules vides
     if (ir < n_regions) then
        if ( (regions(ir+1)%Rmin > regions(ir)%Rmax) ) then
           ln_delta_r = (1.0_dp/real(n_empty+1,kind=dp))*log(regions(ir+1)%Rmin/regions(ir)%Rmax)
           delta_r = exp(ln_delta_r)
           do i=istart+n_rad_region+1, istart+n_rad_region+n_empty
              tab_r(i) = tab_r(i-1) * delta_r
              tab_r2(i) = tab_r(i) * tab_r(i)
              tab_r3(i) = tab_r2(i) * tab_r(i)
           enddo
           n_cells_tmp = n_cells_tmp + n_empty
        endif
     endif

     istart = n_cells_tmp+1
  enddo ! ir

  r_lim(0)= rmin
  r_lim_2(0)= rmin**2
  r_lim_3(0) = rmin**3 ! only for spherical
  do i=1, n_rad
     r_lim(i)=tab_r(i+1)
     r_lim_2(i)= tab_r2(i+1)
     r_lim_3(i)= tab_r3(i+1)
     if (r_lim(i) < r_lim(i-1)) then
        write(*,*) "ERROR in gridding: this is likely to be a bug"
        write(*,*) "i", i, r_lim(i), r_lim(i-1)
        write(*,*) "Exiting"
        stop
     endif
  enddo !i

  if (lcylindrical) then
     ! Calcul volume des cellules (pour calculer leur masse)
     ! On prend ici le rayon au milieu de la cellule
     ! facteur 2 car symétrie
     ! tab_r est en cylindrique ici

     do i=1, n_rad
        rcyl = 0.5*(r_lim(i) +r_lim(i-1))
        r_grid_tmp(i,:) = rcyl!sqrt(r_lim(i) +r_lim(i-1)))

        ! Estimation du zmax proprement
        ! Recherche de l'echelle de hauteur max des zones pertinentes au rayon donne
        H = 0.
        do izone=1,n_zones
           dz=disk_zone(izone)
           if ((dz%rmin < rcyl).and.(rcyl < dz%Rmax)) then
              hzone = dz%sclht * (rcyl/dz%rref)**dz%exp_beta
              if (hzone > H) H = hzone
           endif ! test rcyl
        enddo ! izone
        zmax(i) = cutoff * H
     enddo ! i

     do i=1, n_rad
        ! Interpolation pour les cellules ou H n'est pas defini (ie entre les zones)
        if (zmax(i) < tiny_real)  then
           search_min: do ii = i-1, 1, -1
              if (zmax(ii) > tiny_real) then
                 ii_min = ii
                 exit search_min
              endif
           enddo search_min !ii

           search_max: do ii = i+1, n_rad
              if (zmax(ii) > tiny_real) then
                 ii_max = ii
                 exit search_max
              endif
           enddo search_max !ii

           ! Interpolation lineaire en log(r)
           rcyl = r_grid_tmp(i,1) ; rcyl_min =  r_grid_tmp(ii_min,1)  ; rcyl_max =  r_grid_tmp(ii_max,1)
           frac = (log(rcyl) - log(rcyl_min)) / (log(rcyl_max) - log(rcyl_min))
           zmax(i) = exp(log(zmax(ii_max)) * frac + log(zmax(ii_min)) * (1.0 - frac))
        endif ! zmax(i) < tiny_real
     enddo !i

     do i=1, n_rad
        if ((tab_r2(i+1)-tab_r2(i)) > 1.0e-6*tab_r2(i)) then
           V(i)=2.0_dp*pi*(tab_r2(i+1)-tab_r2(i)) * zmax(i)/real(nz)
           dr2_grid(i) = tab_r2(i+1)-tab_r2(i)
        else
           rcyl = r_grid_tmp(i,1)
           V(i)=4.0_dp*pi*rcyl*(tab_r(i+1)-tab_r(i)) * zmax(i)/real(nz)
           dr2_grid(i) = 2.0_dp * rcyl*(tab_r(i+1)-tab_r(i))
        endif

        delta_z(i)=zmax(i)/real(nz)
        ! Pas d'integration = moitie + petite dimension cellule
        z_lim(i,nz+1)=zmax(i)

        do j=1,nz
           z_lim(i,j) = (real(j,kind=dp)-1.0_dp)*delta_z(i)
           z_grid_tmp(i,j) = (real(j,kind=dp)-0.5_dp)*delta_z(i)
        enddo
     enddo

     z_lim(:,nz+2)=1.0e30
     zmaxmax = maxval(zmax)

  else !lspherical
     izone=1
     dz=disk_zone(izone)


     ! tab_r est en spherique ici
     w_lim(0) = 0.0_dp
     theta_lim(0) = 0.0_dp
     tan_theta_lim(0) = 1.0e-10_dp

     w_lim(nz) = 1.0_dp
     theta_lim(nz) = pi/2.
     tan_theta_lim(nz) = 1.e30_dp

     do j=1, nz-1
        ! repartition uniforme en cos
        w= real(j,kind=dp)/real(nz,kind=dp)
        w_lim(j) = w
        tan_theta_lim(j) = w / sqrt(1.0_dp - w*w)
        theta_lim(j) = atan(tan_theta_lim(j))
     enddo

     do i=1, n_rad
        !rsph = 0.5*(r_lim(i) +r_lim(i-1))
        rsph = sqrt(r_lim(i) * r_lim(i-1))

        do j=1,nz
           w = (real(j,kind=dp)-0.5_dp)/real(nz,kind=dp)
           uv = sqrt(1.0_dp - w*w)
           r_grid_tmp(i,j)=rsph * uv
           z_grid_tmp(i,j)=rsph * w
        enddo

        if (rsph > dz%Rmax) then
           izone = izone +1
           dz=disk_zone(izone)
        endif

        if ((tab_r3(i+1)-tab_r3(i)) > 1.0e-6*tab_r3(i)) then
           V(i)=4.0/3.0*pi*(tab_r3(i+1)-tab_r3(i)) /real(nz)
        else
           V(i)=4.0*pi*rsph**2*(tab_r(i+1)-tab_r(i)) /real(nz)
        endif
     enddo

  endif ! cylindrique ou spherique

  ! Version 3D
  if (l3D) then
     do k=1, n_az
        phi_grid_tmp(k) = 2.0*pi*real(k)/real(n_az)
        phi = phi_grid_tmp(k)
        if (abs(modulo(phi-0.5*pi,pi)) < 1.0e-6) then
           tan_phi_lim(k) = 1.0d300
        else
           tan_phi_lim(k) = tan(phi)
        endif
     enddo !k

     V(:) = V(:) * 0.5 / real(n_az)

     do j=1,nz
        r_grid_tmp(:,-j) = r_grid_tmp(:,j)
        z_grid_tmp(:,-j) = -z_grid_tmp(:,j)
     enddo
  endif

  ! Determine the zone for each cell
  do ir = 1, n_regions
     do i=1, n_rad
        if ((r_grid_tmp(i,1) >  regions(ir)%Rmin).and.(r_grid_tmp(i,1) <  regions(ir)%Rmax)) then
           tab_region(i) = ir
        endif
     enddo
  enddo

  ! Volume and cell arrays with 1D index
  do icell=1, n_cells
     i = cell_map_i(icell)
     j = cell_map_j(icell)
     k = cell_map_k(icell)
     volume(icell) = V(i)

     r_grid(icell) = r_grid_tmp(i,j)
     z_grid(icell) = z_grid_tmp(i,j)
     phi_grid(icell) = phi_grid_tmp(k)
  enddo


  ! Pour Sebastien Charnoz
  if (lSeb_Charnoz) then
     write(*,*) "# n_rad nz"
     write(*,*) n_rad, nz
     write(*,*) "# ir	iz	Rmin		deltaR			Zmin		deltaZ"
     j = 1
     do i=1, n_rad
        do j=1, nz
           write(*,'(I3,3X,I3,3X,ES16.9,3X,ES16.9,3X,ES16.9,3X,ES16.9)') &
                i, j, r_lim(i-1), r_lim(i) - r_lim(i-1), z_lim(i,j),  z_lim(i,j+1) -  z_lim(i,j)
        enddo
     enddo
     stop
  endif ! lSeb_Charnoz

  deallocate(r_grid_tmp,z_grid_tmp,phi_grid_tmp)

  return

end subroutine define_cylindrical_grid

!******************************************************************************

  pure logical function test_exit_grid_cyl(icell, x, y, z)

    integer, intent(in) :: icell
    real(kind=dp), intent(in) :: x,y,z

    if (icell <= n_cells) then
       test_exit_grid_cyl = .false.
       return
    endif

    if (lexit_cell(icell)==0) then
       test_exit_grid_cyl = .false.
    else if (lexit_cell(icell)==1) then ! radial
       test_exit_grid_cyl = .true.
    else ! 2 --> vertical
       if (abs(z) > zmaxmax) then
          test_exit_grid_cyl = .true.
       else
          test_exit_grid_cyl = .false.
       endif
    endif

    return

  end function test_exit_grid_cyl

  !******************************************************************************


  subroutine test_convert()

    integer :: i, j, k, icell
    integer :: i2,j2,k2


    write(*,*)
    write(*,*) "TEST CONVERT"

    do k=1, n_az
       do j=1, nz+1
          do i=0, n_rad

             icell = cell_map(i,j,k)
             write(*,*) "convert", i,j,k, "-->", icell

             call cell2cylindrical(icell, i2,j2,k2)
             if (i>0) then
                if ((i/=i2).or.(j/=j2).or.(k2/=k)) then
                   write(*,*) "PB test convert"
                   write(*,*) i,j,k, "-->", icell
                   write(*,*) icell, "-->", i2,j2,k2
                   stop
                endif
             else
                if ((i/=i2)) then ! seul i est defini ds la cas 0
                   write(*,*) "PB test convert"
                   write(*,*) i,j,k, "-->", icell
                   write(*,*) icell, "-->", i2,j2,k2
                   stop
                endif
             endif
          enddo
       enddo
    enddo

    write(*,*) "DONE"
    stop
    return


  end subroutine test_convert

  !******************************************************************************
  !pure subroutine cylindrical2cell(i,j,k, icell)
  !
  !  integer, intent(in) :: i,j,k
  !  integer, intent(out) :: icell
  !
  !  icell = cell_map(i,j,k)
  !
  !  return
  !
  !end subroutine cylindrical2cell

  !******************************************************************************

  pure subroutine cell2cylindrical(icell, i,j,k) ! work also in sph

    integer, intent(in) :: icell
    integer, intent(out) :: i,j,k

    i = cell_map_i(icell)
    j = cell_map_j(icell)
    k = cell_map_k(icell)

    return

  end subroutine cell2cylindrical

  !******************************************************************************

  subroutine cylindrical2cell_old(i,j,k, icell)
    ! icell is between 1 and n_rad * (n_z+1) * n_az

    integer, intent(in) :: i,j,k
    integer, intent(out) :: icell

    if ((i==0).and.(j==0)) then
       icell = 0
    else if (j>nz+1) then
       icell = -i
    else
       icell = i + n_rad * ( j-1 + nz * (k-1))
    endif

    return

  end subroutine cylindrical2cell_old

  !******************************************************************************

  subroutine cell2cylindrical_old(icell, i,j,k)

    integer, intent(in) :: icell
    integer, intent(out) :: i,j,k

    integer :: ij ! indice combine i et j, ie : i + (j-1) * n_rad

    if (icell==0) then
       i=0
       j=0
       k=1
    else if (icell < 0) then
       i = -icell
       j = nz+2
       k = 1
    else
       k = (icell-1)/nrz + 1 ; if (k > n_az) k=n_az

       ij = icell - (k-1)*nrz
       j = (ij-1)/n_rad + 1 ; if (j > nz+1) j=nz+1

       i = ij - (j-1)*n_rad

       !write(*,*) "TEST ij", ij,  ij/n_rad
       !write(*,*) "i,j", i, j
    endif

    return

  end subroutine cell2cylindrical_old

  !******************************************************************************

  subroutine indice_cellule_cyl(xin,yin,zin, icell)

    implicit none

    real(kind=dp), intent(in) :: xin,yin,zin
    integer, intent(out) :: icell

    real(kind=dp) :: r2, phi
    integer :: ri, ri_min, ri_max, ri_out, zj_out, phik_out

    r2 = xin*xin+yin*yin

    if (r2 < r_lim_2(0)) then
       ri_out=0
       if (abs(zin) > zmax(1)) then
          zj_out = nz+1
          if (l3D) then
             if (zin < 0.0)  zj_out = -zj_out
          endif
       else
          zj_out=1
       endif
       phik_out=1
    else if (r2 > Rmax2) then
       ri_out=n_rad+1
       zj_out=1
       phik_out=1
    else
       ri_min=0
       ri_max=n_rad
       ri=(ri_min+ri_max)/2

       do while((ri_max-ri_min) > 1)
          if(r2 > r_lim_2(ri)) then
             ri_min=ri
          else
             ri_max=ri
          endif
          ri=(ri_min+ri_max)/2
       enddo
       ri_out=ri+1

       zj_out = floor(min(real(abs(zin)/zmax(ri_out) * nz),max_int))+1

       if (l3D) then
          if (zj_out > nz) zj_out = nz+1
          if (zin < 0.0)  zj_out = -zj_out
          if (zin /= 0.0) then
             phi=modulo(atan2(yin,xin),2*real(pi,kind=dp))
             phik_out=floor(phi/(2*pi)*real(N_az))+1
             if (phik_out==n_az+1) phik_out=n_az
          else
             phik_out=1
          endif
       else ! 2D
          if (zj_out > nz) zj_out = nz + 1
          phik_out=1
       endif
    endif

    icell = cell_map(ri_out,zj_out,phik_out)

    return

  end subroutine indice_cellule_cyl

  !******************************************************************************

  subroutine indice_cellule_3D_phi(xin,yin,zin,phik_out)
    ! ok : not necessary anymore, included directly dans cross_cylindrical_cell

    implicit none

    real(kind=dp), intent(in) :: xin,yin,zin
    integer, intent(out) :: phik_out

    real(kind=dp) :: phi

    if (zin /= 0.0) then
       phi=modulo(atan2(yin,xin),2*real(pi,kind=dp))
       phik_out=floor(phi/(2*pi)*real(N_az))+1
       if (phik_out==n_az+1) phik_out=n_az
    else
       phik_out=1
    endif

    return

  end subroutine indice_cellule_3D_phi

  !******************************************************************************

  subroutine cross_cylindrical_cell(x0,y0,z0, u,v,w,  cell, previous_cell, x1,y1,z1, next_cell, l)

    integer, intent(in) :: cell, previous_cell
    real(kind=dp), intent(in) :: x0,y0,z0
    real(kind=dp), intent(in) :: u,v,w ! Todo : check that

    real(kind=dp), intent(out) :: x1, y1, z1
    integer, intent(out) :: next_cell
    real(kind=dp), intent(out) :: l

    ! Variables to be sorted out
    integer :: ri0,zj0,k0, k0m1
    integer ::  delta_rad, delta_zj, delta_phi, ri1, zj1, k1

    real(kind=dp) :: inv_a, a, b, c, s, rac, t, t_phi, delta, inv_w, r_2, den, tan_angle_lim
    real(kind=dp) :: phi, delta_vol, zlim, dotprod
    real(kind=dp) :: correct_moins, correct_plus


    ! TODO: Can be calculated outside
    correct_moins = 1.0_dp - prec_grille
    correct_plus = 1.0_dp + prec_grille

    a=u*u+v*v
    if (a > tiny_real) then
       inv_a=1.0_dp/a
    else
       inv_a=huge_real
    endif

    if (abs(w) > tiny_real) then
       inv_w=1.0_dp/w
    else
       inv_w=sign(huge_dp,w) ! huge_real avant
    endif
    ! End : TODO : Can be calculated outside

    ! 3D cell indices
    call cell2cylindrical(cell, ri0,zj0,k0)

    ! Detection interface
    r_2=x0*x0+y0*y0
    b=(x0*u+y0*v)*inv_a

    if (ri0==0) then
       ! Si on est avant le bord interne,  on passe forcement par rmin
       ! et on cherche forcement la racine positive (unique)
       c=(r_2-r_lim_2(0))*inv_a
       delta=b*b-c
       rac=sqrt(delta)
       s = (-b+rac) * correct_plus
       t=huge_real
       t_phi= huge_real
       delta_rad=1
    else
       ! 1) position interface radiale
       ! on avance ou recule en r ? -> produit scalaire
       dotprod=u*x0+v*y0  ! ~ b
       if (dotprod < 0.0_dp) then
          ! on recule : on cherche rayon inférieur
          c=(r_2-r_lim_2(ri0-1)*correct_moins)*inv_a
          delta=b*b-c
          if (delta < 0.0_dp) then ! on ne rencontre pas le rayon inférieur
             ! on cherche le rayon supérieur
             c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
             delta=max(b*b-c,0.0_dp) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
             delta_rad=1
          else
             delta_rad=-1
          endif
       else
          ! on avance : on cherche le rayon supérieur
          c=(r_2-r_lim_2(ri0)*correct_plus)*inv_a
          delta=max(b*b-c,0.0_dp) ! on force 0.0 si pb de precision qui donnerait delta=-epsilon
          delta_rad=1
       endif !dotprod
       rac=sqrt(delta)
       s=(-b-rac) * correct_plus
       if (s < 0.0_dp) then
          s=(-b+rac) * correct_plus
       else if (s==0.0_dp) then
          s=prec_grille
       endif


       ! 2) position interface verticale
       ! on monte ou on descend par rapport au plan équatorial ?
       dotprod=w*z0
       if (dotprod == 0.0_dp) then
          t=1.0e10
       else
          if (dotprod > 0.0_dp) then
             ! on s'eloigne du midplane (ou on monte en 2D)
             if (abs(zj0)==nz+1) then
                delta_zj=0
                zlim=sign(1.0e10_dp,z0)
             else
                zlim= sign(z_lim(ri0,abs(zj0)+1)*correct_plus, z0)  ! BUUG HERE TODO
                delta_zj=1
                if (l3D.and.(z0 < 0.0))  delta_zj=-1
             endif
          else
             !  on se rapproche du midplane (ou on descend en 2D)
             if (l3D) then
                if (z0 > 0.0) then
                   zlim=z_lim(ri0,abs(zj0))*correct_moins
                   delta_zj=-1
                   if (zj0==1) delta_zj=-2 ! pas d'indice 0
                else
                   zlim=-z_lim(ri0,abs(zj0))*correct_moins
                   delta_zj=1
                   if (zj0==-1) delta_zj=2 ! pas d'indice 0
                endif
             else ! 2D
                if (zj0==1) then
                   ! on traverse le plan eq donc on va remonter
                   ! et z va changer de signe
                   delta_zj=1
                   if (z0 > 0.0_dp) then
                      zlim=-z_lim(ri0,2)*correct_moins
                   else
                      zlim=z_lim(ri0,2)*correct_moins
                   endif
                else !(zj0==1)
                   ! on ne traverse pas z=0.
                   if (z0 > 0.0_dp) then
                      zlim=z_lim(ri0,zj0)*correct_moins
                   else
                      zlim=-z_lim(ri0,zj0)*correct_moins
                   endif
                   delta_zj=-1
                endif !(zj0==1)
             endif ! 3D
          endif ! monte ou descend
          t=(zlim-z0)*inv_w
          ! correct pb precision
          if (t < 0.0_dp) t=prec_grille
       endif !dotprod=0.0


       ! 3) position interface azimuthale
       if (l3D) then
          dotprod =  x0*v - y0*u
          if (abs(dotprod) < 1.0e-10) then
             ! on ne franchit pas d'interface azimuthale
             t_phi = 1.0e30
          else
             ! Quelle cellule on va franchir
             if (dotprod > 0.0) then
                tan_angle_lim = tan_phi_lim(k0)
                delta_phi=1
             else
                k0m1=k0-1
                if (k0m1==0) k0m1=N_az
                tan_angle_lim = tan_phi_lim(k0m1)
                delta_phi=-1
             endif
             ! Longueur av interserction
             if (tan_angle_lim > 1.0d299) then
                t_phi = -x0/u
             else
                den= v-u*tan_angle_lim
                if (abs(den) > 1.0e-6) then
                   t_phi = -(y0-x0*tan_angle_lim)/den
                else
                   t_phi = 1.0e30
                endif
             endif
             if (t_phi < 0.0) t_phi = 1.0e30
          endif !dotprod = 0.0
       else ! l3D
          t_phi = huge_real
       endif
    endif ! ri0==0


    ! 4) interface en r ou z ?
    if ((s < t).and.(s < t_phi)) then ! r
       l=s
       delta_vol=s
       ! Position au bord de la cellule suivante
       x1=x0+delta_vol*u
       y1=y0+delta_vol*v
       z1=z0+delta_vol*w
       ri1=ri0+delta_rad
       if ((ri1<1).or.(ri1>n_rad)) then
          zj1=zj0
       else
          zj1= floor(min(real(abs(z1)/zmax(ri1)*nz),real(max_int))) + 1
          if (zj1>nz) zj1=nz+1
          if (l3D.and.(z1 < 0.0)) zj1=-zj1
       endif

       k1=k0
       if (l3D) then
          ! We need to find the azimuth when we enter the disc
          ! It can be different from the initial azimuth if the star is not centered
          ! so we need to compute it here
          if (ri0==0) then
             phi=modulo(atan2(y1,x1),2*real(pi,kind=dp))
             k1=floor(phi*un_sur_deux_pi*real(N_az))+1
             if (k1==n_az+1) k1=n_az
          endif
       endif

    else if (t < t_phi) then ! z
       l=t
       delta_vol=t
       ! Position au bord de la cellule suivante
       x1=x0+delta_vol*u
       y1=y0+delta_vol*v
       z1=z0+delta_vol*w
       ri1=ri0
       zj1=zj0+delta_zj
       k1=k0
    else ! phi --> only happens in 3D
       l=t_phi
       delta_vol=correct_plus*t_phi
       ! Position au bord de la cellule suivante
       x1=x0+delta_vol*u
       y1=y0+delta_vol*v
       z1=z0+delta_vol*w
       ri1=ri0
       zj1= floor(abs(z1)/zmax(ri1)*nz) + 1
       if (zj1>nz) zj1=nz+1
       if (z1 < 0.0) zj1=-zj1
       k1=k0+delta_phi
       if (k1 == 0) k1=N_az
       if (k1 == N_az+1) k1=1
    endif

    ! Correction if z1==0, otherwise dotprod (in z) will be 0 at the next iteration
    if (z1 == 0.0_dp) then
       if (l3D) then
          z1 = sign(prec_grille,w)
       else
          z1 = prec_grille
       endif
    endif

    !call cylindrical2cell(ri1,zj1,1, next_cell)
    next_cell = cell_map(ri1,zj1,k1)

    return

  end subroutine cross_cylindrical_cell

  !***********************************************************

  subroutine verif_cell_position_cyl(icell, x, y, z)

    real(kind=dp), intent(inout) :: x,y,z
    integer, intent(inout) :: icell

    integer :: ri, zj, ri0, zj0, tmp_k
    real(kind=dp) :: factor, correct_moins, correct_plus

    correct_moins = 1.0_dp - prec_grille
    correct_plus = 1.0_dp + prec_grille

    ! todo : tmp :
    call cell2cylindrical(icell, ri0,zj0, tmp_k) ! converting current cell index

    ! locate current cell index
    call indice_cellule_cyl(x,y,z, icell)
    ri = cell_map_i(icell)

    ! Patch pour eviter BUG sur position radiale
    ! a cause de limite de precision
    if (ri==0) then
       factor = rmin/ sqrt(x*x+y*y) * correct_plus
       x = x * factor
       y = y * factor
       z = z * factor

       ! On verifie que c'est OK maintenant
       call indice_cellule_cyl(x,y,z, icell)
       ri = cell_map_i(icell)
       if (ri==0) then
          write(*,*) "BUG in verif_cell_position_cyl"
          write(*,*) "Exiting"
          stop
       endif
    endif

    if (l_dark_zone(icell)) then ! Petit test de securite
       zj = cell_map_j(icell)
       ! On resort le paquet
       if (zj < zj0) then
          zj = zj0
          z = z_lim(ri0,zj0)*correct_plus
       endif
       if (ri < ri0) then
          ri = ri0
          x = x * correct_plus
          y = y * correct_plus
       else if (ri > ri0) then
          ri = ri0
          x = x * correct_moins
          y = y * correct_moins
       endif
    endif

  end subroutine verif_cell_position_cyl

!**********************************************************************

  subroutine move_to_grid_cyl(id, x,y,z,u,v,w, icell,lintersect)
    ! Calcule la position au bord de la grille dans
    ! la direction donnee pour grille cylindrique
    ! C. Pinte
    ! 19/09/07

    implicit none

    integer, intent(in) :: id
    real(kind=dp), intent(inout) :: x,y,z
    real(kind=dp), intent(in) :: u,v,w
    integer, intent(out) :: icell
    logical, intent(out) :: lintersect

    real(kind=dp) :: x0, y0, z0, z1, a, inv_a, r_2, b, c, delta, rac, s1, s2, dotprod, t1, t2
    real(kind=dp) :: zlim, zlim2, delta_vol, inv_w, correct_moins, correct_plus

    correct_moins = 1.0_dp - 1.0e-10_dp
    correct_plus = 1.0_dp + 1.0e-10_dp

    x0=x ; y0=y ; z0=z

    a=u*u+v*v
    if (a > tiny_real) then
       inv_a=1.0_dp/a
    else
       inv_a=huge_real
    endif

    if (abs(w) > tiny_real) then
       inv_w=1.0_dp/w
    else
       inv_w=sign(huge_dp,w) ! huge_real avant
    endif

    ! Longueur de vol pour atteindre le rayon cylindrique rout
    r_2=x0*x0+y0*y0
    b=(x0*u+y0*v)*inv_a

    c=(r_2-r_lim_2(n_rad)*correct_moins)*inv_a
    delta=b*b-c
    if (delta < 0.0_dp) then
       ! On ne rencontre pas le cylindre
       s1 = huge_real
       s2 = huge_real
    else
       ! On rencontre le cylindre
       rac=sqrt(delta)

       ! Les deux racines doivent etre positives sinon BUG !!
       ! s1 < s2
       s1=-b-rac
       s2=-b+rac

       ! TMP : BUG : ca plante si rayon vertical !!!
       !  if (s1 < 0.0) then
       !     write(*,*) "Bug dans ray tracing !!!", s1, s2
       !  endif
       ! END TMP
    endif

    ! longueur de vol pour atteindre zmax
    ! le rayon monte ou descend ?
    dotprod=w*z0

    if (abs(dotprod) < tiny_real) then ! rayon horizontal
       t1=huge_real
       t2=huge_real
    else
       if (z0 > 0.0_dp) then
          zlim=zmaxmax*correct_moins
          zlim2=-zmaxmax*correct_moins
       else
          zlim=-zmaxmax*correct_moins
          zlim2=zmaxmax*correct_moins
       endif
       t1=(zlim-z0)*inv_w
       t2=(zlim2-z0)*inv_w
    endif !dotprod=0.0

    ! On ne rencontre ni le cylindre ni les plans
    if (t1 > 1e20) then
       if (s1 > 1e20) then
          lintersect = .false.
          return
       endif
    endif

    ! On rentre ?? et si oui, par le dessus ou par rout ??
    if (t1 > s1) then ! On rentre d'abord dans le cylindre
       if (t1 > s2) then ! on ressort du cylindre avant de croiser la tranche
          ! On se place au bord du cylindre
          delta_vol = s1
          z1 = z0 + delta_vol * w
          if (abs(z1) > zmaxmax) then ! On est constamment en dehors des 2 plans
             lintersect=.false.
             return
          else ! on est constamment entre les 2 plans
             lintersect = .true.
          endif
       else  ! on va croiser la surface dans le cylindre
          lintersect = .true.
          delta_vol = t1
       endif

    else  ! on croise d'abord le plan
       if (t2 < s1) then
          ! on ressort de la tranche avant de croiser le cylindre
          lintersect=.false.
          return
       else
          lintersect = .true.
          delta_vol = s1
       endif
    endif

    ! Position au bord de la grille
    x=x0+delta_vol*u!*correct_plus
    y=y0+delta_vol*v!*correct_plus
    z=z0+delta_vol*w!*correct_plus

    ! Determination de l'indice de la premiere cellule traversee
    ! pour initialiser la propagation
    call indice_cellule_cyl(x,y,z, icell)

    return

  end subroutine move_to_grid_cyl

  !**********************************************************************

  subroutine pos_em_cellule_cyl(icell,aleat1,aleat2,aleat3, x,y,z)
    ! Choisit la position d'emission uniformement
    ! dans la cellule (ri,zj)
    ! Geometrie cylindrique
    ! C. Pinte
    ! 04/02/05

    implicit none

    integer, intent(in) :: icell
    real, intent(in) :: aleat1, aleat2, aleat3
    real(kind=dp), intent(out) :: x,y,z

    real(kind=dp) :: r,phi
    integer :: ri, zj, phik

    ri = cell_map_i(icell)
    zj = cell_map_j(icell)
    phik = cell_map_k(icell)

    ! Position aleatoire dans cellule
    ! Position radiale
    !  r=r_lim(ri-1)+aleat1*(r_lim(ri)-r_lim(ri-1))
    !  r=sqrt(r_lim(ri-1)**2+aleat1*(r_lim(ri)**2-r_lim(ri-1)**2))

    r=sqrt(r_lim_2(ri-1)+aleat1*(r_lim_2(ri)-r_lim_2(ri-1)))
    ! Position verticale
    if (l3D) then ! signe de z = signe de zj
       if (zj > 0) then
          z=z_lim(ri,zj)+aleat2*(z_lim(ri,zj+1)-z_lim(ri,zj))
       else
          z= -(z_lim(ri,-zj)+aleat2*(z_lim(ri,-zj+1)-z_lim(ri,-zj)))
       endif
    else ! 2D : choix aléatoire du signe
       if (aleat2 > 0.5_dp) then
          z=z_lim(ri,zj)+(2.0_dp*(aleat2-0.5_dp))*(z_lim(ri,abs(zj)+1)-z_lim(ri,zj))
       else
          z=-(z_lim(ri,zj)+(2.0_dp*aleat2)*(z_lim(ri,zj+1)-z_lim(ri,zj)))
       endif
    endif

    ! Position azimuthale
    !phi=(2.0*aleat3-1.0)*pi
    phi = 2.0_dp*pi * (real(phik,kind=dp)-1.0_dp+aleat3)/real(n_az,kind=dp)

    ! x et y
    x=r*cos(phi)
    y=r*sin(phi)

    return

  end subroutine pos_em_cellule_cyl

end module cylindrical_grid
