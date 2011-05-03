module grid

  use parametres
  use constantes
  use disk
  use opacity
  use grains
  use em_th
  use prop_star
  use mem

  implicit none

  contains


!******************************************************************************

subroutine define_physical_zones()
  ! Recheche les connections de zone 2 a 2
  ! on doit pouvoir faire mieux que ca

  implicit none

  integer :: i, j, index, i_region, iter, ir
  type(disk_zone_type), dimension(n_zones) :: disk_zone_tmp
  type(disk_zone_type) :: dz

  logical, dimension(n_zones) :: zone_scanned
  
  real(kind=db) :: r1, r2, minR, maxR

  logical :: test_j_in_i, test_i_in_j

  allocate(region(n_zones)) 

  ! Detecting connected zones 
  zone_scanned(:) = .false.
  index = 0 
  do i=1, n_zones
     if (.not.zone_scanned(i)) then
        index = index + 1 
        region(i) = index 
        zone_scanned(i) = .true.

        ! Current minimum & maximum radii of region
        minR = disk_zone(i)%rin
        maxR = disk_zone(i)%rout

        ! Besoin d'iterer au cas ou les connections entre zones sont multiples
        ! probablement pas autant mais ca ne coute rien en calcul
        do iter=1, n_zones-1 
           do j=i+1, n_zones
           
              r1 = disk_zone(j)%rin
              r2 = disk_zone(j)%rout
              
              ! Test if the 2 zones are imbrigated
              test_j_in_i = ((r1 > minR).and.(r1 < maxR)) .or. ((r2 > minR).and.(r2 < maxR))
              test_i_in_j = ((minR > r1).and.(minR < r2)) .or. ((minR > r1).and.(maxR < r2))

              if ( test_j_in_i .or. test_i_in_j ) then
                 if (.not.zone_scanned(j)) then
                    i_region = index
                 else
                    i_region = region(j) 
                 endif ! zone_scanned
                 
                 region(j) = i_region
                 zone_scanned(j) = .true.
                 
                 ! Updating minimum & maximum radii of region
                 minR = min(minR,r1)
                 maxR = max(maxR,r2)
              endif ! test rayon
              
           enddo ! j
        enddo ! iter
     endif !.not.zone_scanned(i)
  enddo !i

  

  n_regions = maxval(region(:))

  allocate(Rmin_region(n_regions),Rmax_region(n_regions))

  do ir = 1, n_regions
     Rmin_region(ir) = 1e30 
     Rmax_region(ir) = 0
     do i=1, n_zones
        if (region(i) == ir) then
           Rmin_region(ir) = min(Rmin_region(ir),disk_zone(i)%rin)
           Rmax_region(ir) = max(Rmax_region(ir),disk_zone(i)%rout)
        endif
     enddo !i
  enddo !ir

  !write(*,*) "Number of regions detected =", n_regions
  !do i=1, n_zones
  !   write(*,*) "zone=", i, "region=", region(i), real(Rmin_region(region(i))), real(Rmax_region(region(i)))
  !enddo
    
  return

end subroutine define_physical_zones

!******************************************************************************

subroutine define_grid4()
  implicit none

  real, parameter :: pi = 3.1415926535
  real(kind=db) :: rcyl, puiss, rsph, w, uv, p, rcyl_min, rcyl_max, frac
  real :: phi
  integer :: i,j,k, izone, i_subdivide, iz_subdivide, ii, ii_min, ii_max

  !tab en cylindrique ou spherique suivant grille
  real(kind=db), dimension(n_rad+1) :: tab_r, tab_r2, tab_r3 
  real(kind=db) ::   r_i, r_f, dr, fac, r0, rout_cell, H, hzone
  integer :: ir, iz, n_cells, n_rad_region, n_rad_in_region, n_empty, istart

  real(kind=db), dimension(n_rad-n_rad_in+2) :: tab_r_tmp
  real(kind=db), dimension(n_rad_in+1) :: tab_r_tmp2

  type(disk_zone_type) :: dz

  logical, parameter :: lprint = .false. ! TEMPORARY : the time to validate and test the new routine

 ! calcul des parametres de la table log
!  delta_r = (rout/(rmin))**(1.0/(real(resol)-1.0))

 ! rmin2=rmin*rmin
  rout2=rout*rout
 ! rmin_1 = 1.0/rmin
 ! rmin2_1 = 1.0/(rmin2)

  if (grid_type == 1) then
     lcylindrical = .true.
     lspherical = .false.
  else
     lcylindrical = .false.
     lspherical = .true.
  endif
       
  if ((.not.lcylindrical).and.(is_there_disk)) then
     !w = 0.5_db / real(nz,kind=db)
!     w = 0.15_db  ! Empirique pour echantilloner suffisamment le bord interne verticalement
!    grid_rmin = rmin / sqrt(1.0_db -w*w)     ! ---> Ca cree un bug dans indice_cellule_sph    
     grid_rmin=rmin
  else
     grid_rmin=rmin
  endif

  if (llinear_grid) then
     
     do i=1, n_rad+1
        tab_r(i) = grid_rmin + (rout - grid_rmin) * real(i-1)/real(n_rad)
        tab_r2(i) = tab_r(i) * tab_r(i)
        tab_r3(i) = tab_r2(i) * tab_r(i) 
     enddo

  else 

     ! Definition du nombre de chaques cellules
     n_empty = 3 
     n_rad_region = (n_rad - (n_regions -1) * n_empty) / n_regions
     n_rad_in_region = n_rad_in

     n_cells = 0 

     istart = 1
     tab_r(:) = 0.0_db
     do ir=1, n_regions
        if (lprint) write(*,*) "**********************"
        if (lprint) write(*,*) "New region", ir 
        if (lprint) write(*,*) "istart", istart, n_rad_in_region, n_rad_in
        if (lprint) write(*,*) "R=", Rmin_region(ir), Rmax_region(ir)


        if (ir == n_regions) then
           n_rad_region = n_rad - n_cells ! On prend toutes les celles restantes
        endif

        ! Pour eviter d'avoir 2 cellules a la meme position si les regions se touchent
        R0 =  Rmin_region(ir)
        if (ir > 1) then
           if (Rmin_region(ir) == Rmax_region(ir-1)) then
              R0 =  Rmin_region(ir) * 1.00001_db
           endif
        endif

        ! Grille log avec subdivision cellule interne
        !delta_r = (rout/rmin)**(1.0/(real(n_rad-n_rad_in+1)))
        ln_delta_r = (1.0_db/real(n_rad_region-n_rad_in_region+1,kind=db))*log(Rmax_region(ir)/R0) 
        delta_r = exp(ln_delta_r)

        ln_delta_r_in = (1.0_db/real(n_rad_in_region,kind=db))*log(delta_r)
        delta_r_in = exp(ln_delta_r_in)

        if (lprint) write(*,*) "Delta_r", delta_r, delta_r_in

        ! Selection de la zone correpondante : pente la plus forte
        puiss = -1.e30
        do iz=1, n_zones
           if (region(iz) == ir) then
              p=1+dz%surf-dz%exp_beta
              if (p > puiss) then
                 puiss = p
              endif
           endif
        enddo

        if (lprint) write(*,*) "n_rad_in, puiss=", puiss
        
        ! Calcul recursif hors boucle //
        ! Calcul les rayons separant les cellules de (1 a n_rad + 1)

        tab_r(istart) = R0
        tab_r2(istart) = tab_r(istart) * tab_r(istart)
        tab_r3(istart) = tab_r2(istart) * tab_r(istart)

         if (lprint) write(*,*) istart, ir, tab_r(istart)

        if (puiss == 0.0) then
           do i=istart+1, istart + n_rad_in_region
              tab_r(i) = exp(log(R0) - (log(R0)-log(R0*delta_r))*(2.0**(i-istart)-1.0)/(2.0**n_rad_in_region-1.0))
              tab_r2(i) = tab_r(i) * tab_r(i)
              tab_r3(i) = tab_r2(i) * tab_r(i)

              if (lprint) write(*,*) i, ir, tab_r(i)
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

              if (lprint) write(*,*) i, ir, tab_r(i)
           enddo
        endif
     
        if (lprint) write(*,*) "n_rad"


        do i=istart + n_rad_in_region+1, istart+n_rad_region
           tab_r(i) = tab_r(i-1) * delta_r
           tab_r2(i) = tab_r(i) * tab_r(i)
           tab_r3(i) = tab_r2(i) * tab_r(i)

           if (lprint) write(*,*) i, ir, tab_r(i)
        enddo
        
        
        n_cells = istart+n_rad_region

        ! Cellules vides
        if (ir < n_regions) then
           if ( (Rmin_region(ir+1) > Rmax_region(ir)) ) then
              if (lprint) write(*,*) "empty cells"
              ln_delta_r = (1.0_db/real(n_empty+1,kind=db))*log(Rmin_region(ir+1)/Rmax_region(ir)) 
              delta_r = exp(ln_delta_r)
              do i=istart+n_rad_region+1, istart+n_rad_region+n_empty
                 tab_r(i) = tab_r(i-1) * delta_r
                 tab_r2(i) = tab_r(i) * tab_r(i)
                 tab_r3(i) = tab_r2(i) * tab_r(i)

                 if (lprint) write(*,*) i, ir, tab_r(i)
              enddo
              n_cells = n_cells + n_empty
           endif
        endif

        istart = n_cells+1
     enddo ! ir
     
  endif ! llinear_grid


  
  r_lim(0)= grid_rmin
  r_lim_2(0)= grid_rmin**2
  r_lim_3(0) = grid_rmin**3
  do i=1, n_rad
     r_lim(i)=tab_r(i+1)
     r_lim_2(i)= r_lim(i)**2
     r_lim_3(i)= r_lim(i)**3
  enddo !i


  if (lcylindrical) then
     ! Calcul volume des cellules (pour calculer leur masse)
     ! On prend ici le rayon au milieu de la cellule
     ! facteur 2 car symétrie
     ! tab_r est en cylindrique ici

     do i=1, n_rad
        rcyl = 0.5*(r_lim(i) +r_lim(i-1))        
        r_grid(i,:) = rcyl!sqrt(r_lim(i) +r_lim(i-1)))

        ! Estimation du zmax proprement
        ! Recherche de l'echelle de hauteur max des zones pertinentes au rayon donne
        H = 0.
        do izone=1,n_zones
           dz=disk_zone(izone)
           if ((dz%rmin < rcyl).and.(rcyl < dz%rout)) then
              hzone = dz%sclht * (rcyl/dz%rref)**dz%exp_beta
              if (hzone > H) H = hzone
           endif ! test rcyl
        enddo ! izone
        zmax(i) = cutoff * H
     enddo ! i

     do i=1, n_rad
        ! Interpolation pour les cellules ou H n'est pas defini
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
           rcyl = r_grid(i,1) ; rcyl_min =  r_grid(ii_min,1)  ; rcyl_max =  r_grid(ii_max,1)         
           frac = (log(rcyl) - log(rcyl_min)) / (log(rcyl_max) - log(rcyl_min))
           zmax(i) = exp(log(zmax(ii_max)) * frac + log(zmax(ii_min)) * (1.0 - frac))
        endif ! zmax(i) < tiny_real
     enddo !i

     ! Version basique et initiale : prend la premiere zone pour determiner echelle de hauteur
     !do i=1, n_rad
     !   rcyl = r_grid(i,1)
     !   dz=disk_zone(1)
     !   zmax(i) = cutoff * dz%sclht * (rcyl/dz%rref)**dz%exp_beta
     !enddo !i


     do i=1, n_rad
        if ((tab_r2(i+1)-tab_r2(i)) > 1.0e-6*tab_r2(i)) then
           volume(i)=2.0_db*pi*(tab_r2(i+1)-tab_r2(i)) * zmax(i)/real(nz)
           dr2_grid(i) = tab_r2(i+1)-tab_r2(i)
        else
           volume(i)=4.0_db*pi*rcyl*(tab_r(i+1)-tab_r(i)) * zmax(i)/real(nz)
           dr2_grid(i) = 2.0_db * rcyl*(tab_r(i+1)-tab_r(i)) 
        endif
        ! Conversion en cm**3 
        ! 1 AU = 1.49597870691e13 cm
        !     volume(i)=volume(i)*3.347929e39

        delta_z(i)=zmax(i)/real(nz)
        ! Pas d'integration = moitie + petite dimension cellule 
        if (linteg_dic) delta0(i) = 0.5*min(tab_r(i+1)-tab_r(i),delta_z(i))
        z_lim(i,nz+1)=zmax(i)

        do j=1,nz
           z_lim(i,j) = (real(j,kind=db)-1.0_db)*delta_z(i) 
           z_grid(i,j) = (real(j,kind=db)-0.5_db)*delta_z(i) 
        enddo
     enddo

     z_lim(:,nz+2)=1.0e30
     zmaxmax = maxval(zmax)

     if (linteg_dic) delta0(0)=0.5*grid_rmin

  else !lspherical
     izone=1
     dz=disk_zone(izone)


     ! tab_r est en spherique ici
     w_lim(0) = 0.0_db
     theta_lim(0) = 0.0_db
     tan_theta_lim(0) = 1.0e-10_db
     
     w_lim(nz) = 1.0_db
     theta_lim(nz) = pi/2.
     tan_theta_lim(nz) = 1.e30_db
     
     do j=1, nz-1
        ! repartition uniforme en cos
        w= real(j,kind=db)/real(nz,kind=db)
        w_lim(j) = w
        tan_theta_lim(j) = w / sqrt(1.0_db - w*w)
        theta_lim(j) = atan(tan_theta_lim(j))
 !       write(*,*) "tan_theta_lim", j, w, tan_theta_lim(j)
     enddo
  !   stop

     
     do i=1, n_rad
        !rsph = 0.5*(r_lim(i) +r_lim(i-1))
        rsph = sqrt(r_lim(i) * r_lim(i-1))
        
        do j=1,nz
           w = (real(j,kind=db)-0.5_db)/real(nz,kind=db)
           uv = sqrt(1.0_db - w*w)
           r_grid(i,j)=rsph * uv 
           z_grid(i,j)=rsph * w
        enddo

        if (rsph > dz%rout) then
           izone = izone +1
           dz=disk_zone(izone)
        endif
        
        if ((tab_r3(i+1)-tab_r3(i)) > 1.0e-6*tab_r3(i)) then
           volume(i)=4.0/3.0*pi*(tab_r3(i+1)-tab_r3(i)) /real(nz)
        else
           volume(i)=4.0*pi*rsph**2*(tab_r(i+1)-tab_r(i)) /real(nz)
        endif
     enddo

  endif ! cylindrique ou spherique

  ! Version 3D
  if (l3D) then
     do k=1, n_az
        phi_grid(k) = 2.0*pi*real(k)/real(n_az)
        phi = phi_grid(k)
        if (abs(modulo(phi-0.5*pi,pi)) < 1.0e-6) then
           tan_phi_lim(k) = 1.0d300
        else
           tan_phi_lim(k) = tan(phi)
        endif
     enddo !pk

     volume(:) = volume(:) * 0.5 / real(n_az) 

     do j=1,nz
        z_grid(:,-j) = z_grid(:,j)
     enddo
  endif

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

  return

end subroutine define_grid4

!******************************************************************************

subroutine define_grid3()
! Definit la grille du code
! Calcule les tableaux zmax, volume, r_lim, r_lim_2, z_lim et delta0
! C. Pinte 
! 27/04/05 

  implicit none

  real, parameter :: pi = 3.1415926535
  real(kind=db) :: rcyl, puiss, rsph, w, uv, p, rcyl_min, rcyl_max, frac
  real :: phi
  integer :: i,j,k, izone, i_subdivide, iz_subdivide, ii, ii_min, ii_max

  !tab en cylindrique ou spherique suivant grille
  real(kind=db), dimension(n_rad+1) :: tab_r, tab_r2, tab_r3 
  real(kind=db) ::   r_i, r_f, dr, fac, r0, rout_cell, H, hzone

  real(kind=db), dimension(n_rad-n_rad_in+2) :: tab_r_tmp
  real(kind=db), dimension(n_rad_in+1) :: tab_r_tmp2

  type(disk_zone_type) :: dz

 ! calcul des parametres de la table log
!  delta_r = (rout/(rmin))**(1.0/(real(resol)-1.0))

 ! rmin2=rmin*rmin
  rout2=rout*rout
 ! rmin_1 = 1.0/rmin
 ! rmin2_1 = 1.0/(rmin2)

  if (grid_type == 1) then
     lcylindrical = .true.
     lspherical = .false.
  else
     lcylindrical = .false.
     lspherical = .true.
  endif
       
  if ((.not.lcylindrical).and.(is_there_disk)) then
     !w = 0.5_db / real(nz,kind=db)
!     w = 0.15_db  ! Empirique pour echantilloner suffisamment le bord interne verticalement
!    grid_rmin = rmin / sqrt(1.0_db -w*w)     ! ---> Ca cree un bug dans indice_cellule_sph    
     grid_rmin=rmin
  else
     grid_rmin=rmin
  endif

  if (llinear_grid) then
     
     do i=1, n_rad+1
        tab_r(i) = grid_rmin + (rout - grid_rmin) * real(i-1)/real(n_rad)
        tab_r2(i) = tab_r(i) * tab_r(i)
        tab_r3(i) = tab_r2(i) * tab_r(i) 
     enddo

  else if (lr_subdivide) then
     ln_delta_r = (1.0_db/real(n_rad-n_rad_in+1,kind=db))*log(rout/grid_rmin) 
     delta_r = exp(ln_delta_r)

     ! Repartition des cellules "entieres" en log
     tab_r_tmp(1) = grid_rmin
     do i=2, n_rad - n_rad_in + 2
        tab_r_tmp(i) = tab_r_tmp(i-1) * delta_r
     enddo

     ! Recherche de la cellule a diviser
     search : do i=1, n_rad
        !rcyl = 0.5 * (tab_r_tmp(i) + tab_r_tmp(i+1))
        !if (rcyl > r_subdivide) then
        if (tab_r_tmp(i+1) > r_subdivide) then
           i_subdivide = i
           exit
        endif
     enddo search

     ! Selection de la zone correpondante : pente la plus forte
     iz_subdivide = 0
     puiss = -1.e30
     do i=1, n_zones
        dz=disk_zone(i)
        if ((dz%rin <= r_subdivide).and.(dz%rout >= r_subdivide)) then
           p=1+dz%surf-dz%exp_beta
           if (p > puiss) then
              puiss = p
              iz_subdivide = i
           endif
        endif
     enddo
     if (iz_subdivide == 0) then
        write(*,*) "Error in cell subdivision : no zone at this radius"
        write(*,*) "Exiting."
        stop
     endif
     dz = disk_zone(iz_subdivide)

     write(*,*) "Subdividing cell", i_subdivide," in zone", iz_subdivide
     write(*,*) "Inner cell radius=", tab_r_tmp(i_subdivide)
     write(*,*) "Outer cell radius=", tab_r_tmp(i_subdivide+1)
     

     ! Subdivision
!     ln_delta_r_in = (1.0_db/real(n_rad_in-1,kind=db))*log(delta_r)
!     delta_r_in = exp(ln_delta_r_in)

     r0 = tab_r_tmp(i_subdivide)
     tab_r_tmp2(1) = r0
     tab_r_tmp2(2) = r_subdivide
     rout_cell = tab_r_tmp(i_subdivide+1)
     if (puiss == 0.0) then
        do i=3, n_rad_in+1
           tab_r_tmp2(i) = exp(log(r_subdivide) - (log(r_subdivide)-log(rout_cell))*(2.0**(i-1)-1.0)/(2.0**(n_rad_in-1)-1.0))
        enddo
     else
        r_i = exp(puiss*log(r_subdivide))
        r_f = exp(puiss*log(rout_cell))
        dr=r_f-r_i
        fac = 1.0/(2.0**(n_rad_in)-1.0)
        do i=3, n_rad_in+1
           tab_r_tmp2(i) = (r_subdivide**puiss - (r_subdivide**puiss-(rout_cell)**puiss) & 
                *(2.0**(i)-1.0)/(2.0**(n_rad_in)-1.0))**(1.0/puiss)
           if (tab_r_tmp2(i) - tab_r_tmp2(i-1) < prec_grille*tab_r_tmp2(i-1)) then
              write(*,*) "Error : spatial grid resolution too high"
              write(*,*) "Differences between two cells are below double precision"
              stop
           endif
        enddo
     endif

     ! On remet le tout dans l'ordre
     do i=1,i_subdivide
        tab_r(i) = tab_r_tmp(i)
     enddo

     do i=i_subdivide+1,i_subdivide+n_rad_in
        tab_r(i) = tab_r_tmp2(i - i_subdivide +1)
     enddo

     do i=i_subdivide+n_rad_in+1,n_rad+1
        tab_r(i) = tab_r_tmp(i - n_rad_in + 1)
     enddo

     do i=1,n_rad+1
        tab_r2(i) = tab_r(i) * tab_r(i)
        tab_r3(i) = tab_r2(i) * tab_r(i) 
     enddo

     
  else ! Grille log avec subdivision cellule interne
     !delta_r = (rout/rmin)**(1.0/(real(n_rad-n_rad_in+1)))
     ln_delta_r = (1.0_db/real(n_rad-n_rad_in+1,kind=db))*log(rout/grid_rmin) 
     delta_r = exp(ln_delta_r)
  
     !delta_r_in = delta_r**(1.0/real(n_rad_in))
     ln_delta_r_in = (1.0_db/real(n_rad_in,kind=db))*log(delta_r)
     delta_r_in = exp(ln_delta_r_in)
     
     !   delta_0 = (delta_r - 1.0) * rmin pour length4. Il faut prendre le moitie pour length3
     !delta_0 = (delta_r - 1.0) * rmin !* 0.5
     !delta_test = 2*delta_0
     
     ! Selection de la 1er zone
     dz=disk_zone(1) ! Les zones doivent etre dans l'ordre

     puiss=1+dz%surf-dz%exp_beta
     
     ! Calcul recursif hors boucle //
     ! Calcul les rayons separant les cellules de (1 a n_rad + 1)
     tab_r(1) = grid_rmin
     tab_r2(1) = grid_rmin*grid_rmin
     tab_r3(1) = grid_rmin*grid_rmin*grid_rmin
     if (puiss == 0.0) then
        do i=2, n_rad_in+1
           tab_r(i) = exp(log(grid_rmin) - (log(grid_rmin)-log(grid_rmin*delta_r))*(2.0**(i-1)-1.0)/(2.0**n_rad_in-1.0))
           tab_r2(i) = tab_r(i) * tab_r(i)
           tab_r3(i) = tab_r2(i) * tab_r(i)
        enddo
     else
        r_i = exp(puiss*log(grid_rmin))
        r_f = exp(puiss*log(grid_rmin*delta_r))
        dr=r_f-r_i
        fac = 1.0/(2.0**(n_rad_in+1)-1.0)
        do i=2, n_rad_in+1
           tab_r(i) = (grid_rmin**puiss - (grid_rmin**puiss-(grid_rmin*delta_r)**puiss) & 
                *(2.0**(i)-1.0)/(2.0**(n_rad_in+1)-1.0))**(1.0/puiss)
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
     do i=n_rad_in+2, n_rad+1
        tab_r(i) = tab_r(i-1) * delta_r
        tab_r2(i) = tab_r(i) * tab_r(i)
        tab_r3(i) = tab_r2(i) * tab_r(i)
     enddo
     
  endif ! llinear_grid

  
  r_lim(0)= grid_rmin
  r_lim_2(0)= grid_rmin**2
  r_lim_3(0) = grid_rmin**3
  do i=1, n_rad
     r_lim(i)=tab_r(i+1)
     r_lim_2(i)= r_lim(i)**2
     r_lim_3(i)= r_lim(i)**3
  enddo !i


  if (lcylindrical) then
     ! Calcul volume des cellules (pour calculer leur masse)
     ! On prend ici le rayon au milieu de la cellule
     ! facteur 2 car symétrie
     ! tab_r est en cylindrique ici

     do i=1, n_rad
        rcyl = 0.5*(r_lim(i) +r_lim(i-1))        
        r_grid(i,:) = rcyl!sqrt(r_lim(i) +r_lim(i-1)))

        ! Estimation du zmax proprement
        ! Recherche de l'echelle de hauteur max des zones pertinentes au rayon donne
        H = 0.
        do izone=1,n_zones
           dz=disk_zone(izone)
           if ((dz%rmin < rcyl).and.(rcyl < dz%rout)) then
              hzone = dz%sclht * (rcyl/dz%rref)**dz%exp_beta
              if (hzone > H) H = hzone
           endif ! test rcyl
        enddo ! izone
        zmax(i) = cutoff * H
     enddo ! i

     do i=1, n_rad
        ! Interpolation pour les cellules ou H n'est pas defini
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
           rcyl = r_grid(i,1) ; rcyl_min =  r_grid(ii_min,1)  ; rcyl_max =  r_grid(ii_max,1)         
           frac = (log(rcyl) - log(rcyl_min)) / (log(rcyl_max) - log(rcyl_min))
           zmax(i) = exp(log(zmax(ii_max)) * frac + log(zmax(ii_min)) * (1.0 - frac))
        endif ! zmax(i) < tiny_real
     enddo !i

     ! Version basique et initiale : prend la premiere zone pour determiner echelle de hauteur
     !do i=1, n_rad
     !   rcyl = r_grid(i,1)
     !   dz=disk_zone(1)
     !   zmax(i) = cutoff * dz%sclht * (rcyl/dz%rref)**dz%exp_beta
     !enddo !i


     do i=1, n_rad
        if ((tab_r2(i+1)-tab_r2(i)) > 1.0e-6*tab_r2(i)) then
           volume(i)=2.0_db*pi*(tab_r2(i+1)-tab_r2(i)) * zmax(i)/real(nz)
           dr2_grid(i) = tab_r2(i+1)-tab_r2(i)
        else
           volume(i)=4.0_db*pi*rcyl*(tab_r(i+1)-tab_r(i)) * zmax(i)/real(nz)
           dr2_grid(i) = 2.0_db * rcyl*(tab_r(i+1)-tab_r(i)) 
        endif
        ! Conversion en cm**3 
        ! 1 AU = 1.49597870691e13 cm
        !     volume(i)=volume(i)*3.347929e39

        delta_z(i)=zmax(i)/real(nz)
        ! Pas d'integration = moitie + petite dimension cellule 
        if (linteg_dic) delta0(i) = 0.5*min(tab_r(i+1)-tab_r(i),delta_z(i))
        z_lim(i,nz+1)=zmax(i)

        do j=1,nz
           z_lim(i,j) = (real(j,kind=db)-1.0_db)*delta_z(i) 
           z_grid(i,j) = (real(j,kind=db)-0.5_db)*delta_z(i) 
        enddo
     enddo

     z_lim(:,nz+2)=1.0e30
     zmaxmax = maxval(zmax)

     if (linteg_dic) delta0(0)=0.5*grid_rmin

  else !lspherical
     izone=1
     dz=disk_zone(izone)


     ! tab_r est en spherique ici
     w_lim(0) = 0.0_db
     theta_lim(0) = 0.0_db
     tan_theta_lim(0) = 1.0e-10_db
     
     w_lim(nz) = 1.0_db
     theta_lim(nz) = pi/2.
     tan_theta_lim(nz) = 1.e30_db
     
     do j=1, nz-1
        ! repartition uniforme en cos
        w= real(j,kind=db)/real(nz,kind=db)
        w_lim(j) = w
        tan_theta_lim(j) = w / sqrt(1.0_db - w*w)
        theta_lim(j) = atan(tan_theta_lim(j))
 !       write(*,*) "tan_theta_lim", j, w, tan_theta_lim(j)
     enddo
  !   stop

     
     do i=1, n_rad
        !rsph = 0.5*(r_lim(i) +r_lim(i-1))
        rsph = sqrt(r_lim(i) * r_lim(i-1))
        
        do j=1,nz
           w = (real(j,kind=db)-0.5_db)/real(nz,kind=db)
           uv = sqrt(1.0_db - w*w)
           r_grid(i,j)=rsph * uv 
           z_grid(i,j)=rsph * w
        enddo

        if (rsph > dz%rout) then
           izone = izone +1
           dz=disk_zone(izone)
        endif
        
        if ((tab_r3(i+1)-tab_r3(i)) > 1.0e-6*tab_r3(i)) then
           volume(i)=4.0/3.0*pi*(tab_r3(i+1)-tab_r3(i)) /real(nz)
        else
           volume(i)=4.0*pi*rsph**2*(tab_r(i+1)-tab_r(i)) /real(nz)
        endif
     enddo

  endif ! cylindrique ou spherique

  ! Version 3D
  if (l3D) then
     do k=1, n_az
        phi_grid(k) = 2.0*pi*real(k)/real(n_az)
        phi = phi_grid(k)
        if (abs(modulo(phi-0.5*pi,pi)) < 1.0e-6) then
           tan_phi_lim(k) = 1.0d300
        else
           tan_phi_lim(k) = tan(phi)
        endif
     enddo !pk

     volume(:) = volume(:) * 0.5 / real(n_az) 

     do j=1,nz
        z_grid(:,-j) = z_grid(:,j)
     enddo
  endif

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

  return

end subroutine define_grid3

!******************************************************************************

subroutine indice_cellule(xin,yin,zin,ri_out,zj_out)

  implicit none
  
  real(kind=db), intent(in) :: xin,yin,zin
  integer, intent(out) :: ri_out, zj_out

  real(kind=db) :: r2
  integer :: ri, ri_min, ri_max

  r2 = xin*xin+yin*yin

    
  if (r2 < r_lim_2(0)) then
     ri_out=0
     zj_out=1
     return
  else if (r2 > rout2) then
     ri_out=n_rad
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
  endif
  
  zj_out = floor(min(real(abs(zin)/zmax(ri_out) * nz),max_int))+1
  if (zj_out > nz) then 
     zj_out = nz + 1
  endif

  return

end subroutine indice_cellule

!******************************************************************************

subroutine indice_cellule_sph(xin,yin,zin,ri_out,thetaj_out)

  implicit none
  
  real(kind=db), intent(in) :: xin,yin,zin
  integer, intent(out) :: ri_out, thetaj_out

  real(kind=db) :: r2, r02, tan_theta
  integer :: ri, ri_min, ri_max, thetaj, thetaj_min, thetaj_max

  r02 = xin*xin+yin*yin
  r2 = r02+zin*zin

  ! Peut etre un bug en 0, 0 due a la correction sur grid_rmin dans define_grid3
  if (r2 < r_lim_2(0)) then
     ri_out=0
  else if (r2 > rout2) then
     ri_out=n_rad
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
  endif
  
  ! thetaj_out
  if (r02 > tiny_db) then
     tan_theta = abs(zin) / sqrt(r02)
  else
     tan_theta = 1.0e30
  endif

  thetaj_min = 0
  thetaj_max = nz
  thetaj=(thetaj_min+thetaj_max)/2

  do while((thetaj_max-thetaj_min) > 1)
     if(tan_theta > tan_theta_lim(thetaj)) then
        thetaj_min=thetaj
     else
        thetaj_max=thetaj
     endif
     thetaj=(thetaj_min+thetaj_max)/2
  enddo
  thetaj_out=thetaj+1
  
  return

end subroutine indice_cellule_sph

!******************************************************************************

subroutine indice_cellule_sph_theta(xin,yin,zin,thetaj_out)

  implicit none
  
  real(kind=db), intent(in) :: xin,yin,zin
  integer, intent(out) :: thetaj_out

  real(kind=db) :: r02, tan_theta
  integer :: thetaj, thetaj_min, thetaj_max

  r02 = xin*xin+yin*yin
  
  ! thetaj_out
  if (r02 > tiny_db) then
     tan_theta = abs(zin) / sqrt(r02)
  else
     tan_theta = 1.0e30
  endif

  thetaj_min = 0
  thetaj_max = nz
  thetaj=(thetaj_min+thetaj_max)/2

  do while((thetaj_max-thetaj_min) > 1)
     if(tan_theta > tan_theta_lim(thetaj)) then
        thetaj_min=thetaj
     else
        thetaj_max=thetaj
     endif
     thetaj=(thetaj_min+thetaj_max)/2
  enddo
  thetaj_out=thetaj+1
  
  return

end subroutine indice_cellule_sph_theta

!******************************************************************************

subroutine indice_cellule_3D(xin,yin,zin,ri_out,zj_out,phik_out)

  implicit none
  
  real(kind=db), intent(in) :: xin,yin,zin
  integer, intent(out) :: ri_out, zj_out, phik_out

  real(kind=db) :: r2, phi
  integer :: ri, ri_min, ri_max

  r2 = xin*xin+yin*yin

    
  if (r2 < r_lim_2(0)) then
     ri_out=0
  else if (r2 > rout2) then
     ri_out=n_rad
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
  endif
  
  zj_out = floor(min(real(abs(zin)/zmax(ri_out) * nz),max_int))+1
  if (zj_out > nz) zj_out = nz
  if (zin < 0.0)  zj_out = -zj_out


  if (zin /= 0.0) then
     phi=modulo(atan2(yin,xin),2*real(pi,kind=db))
     phik_out=floor(phi/(2*pi)*real(N_az))+1
     if (phik_out==n_az+1) phik_out=n_az
  else
     phik_out=1
  endif

  return

end subroutine indice_cellule_3D

!******************************************************************************

subroutine init_lambda()
  ! Initialisation table de longueurs d'onde

  implicit none

  real :: facteur_largeur
  integer :: i

  if (lmono0) then
     ! Lecture longueur d'onde
     read(band,*) tab_lambda(1)
     tab_delta_lambda(1) = 1.0

  else
     ! Initialisation longueurs d'onde
     !delta_lambda = (lambda_max/lambda_min)**(1.0/real(n_lambda))
     delta_lambda =  exp((1.0_db/real(n_lambda,kind=db)) * log(lambda_max/lambda_min))
     
     tab_lambda(1)=lambda_min*sqrt(delta_lambda)
     do i=2, n_lambda
        tab_lambda(i)= tab_lambda(i-1)*delta_lambda
     enddo
     
     facteur_largeur=sqrt(delta_lambda)-1.0/sqrt(delta_lambda)
     do i=1, n_lambda
        tab_delta_lambda(i) = tab_lambda(i) * facteur_largeur
     enddo

  endif

end subroutine init_lambda

!**********************************************************************

subroutine init_lambda2()
  ! Initialisation table en lambda sed

  implicit none

  integer :: i
  
  ! reorganisation memoire
  call realloc_step2
  
  n_lambda=n_lambda2
  do i=1, n_lambda2
     tab_lambda(i)= tab_lambda2(i)
     tab_delta_lambda(i)= tab_delta_lambda2(i)
  enddo
  

  return
  
end subroutine init_lambda2

!**********************************************************************

subroutine select_cellule(lambda,aleat,ri,zj, phik)
! Sélection de la cellule qui va émettre le photon 
! C. Pinte 
! 04/02/05
! Modif 3D 10/06/05

  implicit none

  integer, intent(in) :: lambda
  real, intent(in) :: aleat
  integer, intent(out) :: ri,zj, phik
  integer :: k, k2, kmin, kmax, nz2
  
  
  ! Dichotomie
  kmin=0
  if (l3D) then
     nz2=2*nz
     kmax=n_rad*nz2*n_az
  else
     kmax=n_rad*nz
  endif
  k=(kmin+kmax)/2

  do while ((kmax-kmin) > 1)
     if (prob_E_cell(lambda,k) < aleat) then
        kmin = k
     else
        kmax = k
     endif       
     k = (kmin + kmax)/2
   enddo   ! while
   k=kmax

   if (l3D) then
      ! Reconstruction indices 3D (et 2D)
      phik = mod(k,n_az)
      if (phik==0) phik=n_az
      k2=(k-phik)/n_az + 1
      zj=mod(k2,nz2)
      if (zj==0) zj=nz2
      ri=1+(k2-zj)/nz2
      ! recentrage par rapport à 0
      zj = zj -nz  
      if (zj <= 0) zj = zj -1
   else
      ! Reconstruction indices 2D
      phik=1
      zj=mod(k,nz)
      if (zj==0) zj=nz
      ri=1+(k-zj)/nz
   endif

   return

end subroutine select_cellule

!**********************************************************************

subroutine pos_em_cellule(ri,zj,phik,aleat1,aleat2,aleat3,x,y,z)
! Choisit la poistion d'emission uniformement dans la cellule
! C. Pinte
! 8/06/07

  implicit none

  integer, intent(in) :: ri, zj, phik
  real, intent(in) :: aleat1, aleat2, aleat3
  real(kind=db), intent(out) :: x,y,z

  if (lspherical) then
     call pos_em_cellule_sph(ri,zj,phik,aleat1,aleat2,aleat3,x,y,z)
  else
     call pos_em_cellule_cyl(ri,zj,phik,aleat1,aleat2,aleat3,x,y,z)
  endif

  return

end subroutine pos_em_cellule

!**********************************************************************

subroutine pos_em_cellule_cyl(ri,zj,phik,aleat1,aleat2,aleat3,x,y,z)
! Choisit la position d'emission uniformement
! dans la cellule (ri,zj)
! Geometrie cylindrique
! C. Pinte 
! 04/02/05

  implicit none

  integer, intent(in) :: ri, zj, phik
  real, intent(in) :: aleat1, aleat2, aleat3
  real(kind=db), intent(out) :: x,y,z

  real(kind=db) :: r,phi

  ! Position aleatoire dans cellule
  ! Position radiale
!  r=r_lim(ri-1)+aleat1*(r_lim(ri)-r_lim(ri-1))
!  r=sqrt(r_lim(ri-1)**2+aleat1*(r_lim(ri)**2-r_lim(ri-1)**2))
  
  ! La premiere cellule ne peut plus etre dans la zone sombre maintenant
 ! if (ri==1) then
 !    r=sqrt(rmin2+aleat1*(r_in_opacite2(zj,phik)-rmin2))
 ! else
     r=sqrt(r_lim_2(ri-1)+aleat1*(r_lim_2(ri)-r_lim_2(ri-1)))
 ! endif
  ! Position verticale
  if (l3D) then ! signe de z = signe de zj
     if (zj > 0) then
        z=z_lim(ri,zj)+aleat2*(z_lim(ri,zj+1)-z_lim(ri,zj))
     else
        z= -(z_lim(ri,-zj)+aleat2*(z_lim(ri,-zj+1)-z_lim(ri,-zj)))
     endif
  else ! 2D : choix aléatoire du signe
     if (aleat2 > 0.5_db) then
        z=z_lim(ri,zj)+(2.0_db*(aleat2-0.5_db))*(z_lim(ri,abs(zj)+1)-z_lim(ri,zj))
     else
        z=-(z_lim(ri,zj)+(2.0_db*aleat2)*(z_lim(ri,zj+1)-z_lim(ri,zj)))
     endif
  endif

  
  ! Position azimuthale
  !phi=(2.0*aleat3-1.0)*pi
  phi = 2.0_db*pi * (real(phik,kind=db)-1.0_db+aleat3)/real(n_az,kind=db)

  ! x et y
  x=r*cos(phi)
  y=r*sin(phi)

  return

end subroutine pos_em_cellule_cyl

!***********************************************************

subroutine pos_em_cellule_sph(ri,thetaj,phik,aleat1,aleat2,aleat3,x,y,z)
! Choisit la position d'emission uniformement
! dans la cellule (ri,thetaj)
! Geometrie spherique
! C. Pinte 
! 08/06/07

  implicit none

  integer, intent(in) :: ri, thetaj, phik
  real, intent(in) :: aleat1, aleat2, aleat3
  real(kind=db), intent(out) :: x,y,z

  real(kind=db) :: r, theta, phi, r_cos_theta, w

  integer :: ri_tmp, thetaj_tmp


  ! Position radiale
  r=(r_lim_3(ri-1)+aleat1*(r_lim_3(ri)-r_lim_3(ri-1)))**un_tiers
  
  ! Position theta
  if (aleat2 > 0.5_db) then
     theta=theta_lim(thetaj-1)+(2.0_db*(aleat2-0.5_db))*(theta_lim(thetaj)-theta_lim(thetaj-1))
  else
     theta=-(theta_lim(thetaj-1)+(2.0_db*aleat2)*(theta_lim(thetaj)-theta_lim(thetaj-1)))
  endif

  theta=theta_lim(thetaj-1)+aleat2*(theta_lim(thetaj)-theta_lim(thetaj-1))

  
  ! BUG ??? : ca doit etre uniforme en w, non ??


  ! Position azimuthale
  phi = 2.0_db*pi * (real(phik)-1.0_db+aleat3)/real(n_az)

  ! x et y
  z=r*sin(theta)
  r_cos_theta = r*cos(theta)
  x=r_cos_theta*cos(phi)
  y=r_cos_theta*sin(phi)



!!$  ! Position theta
!!$  if (aleat2 > 0.5) then
!!$     w=w_lim(thetaj-1)+(2.0_db*(aleat2-0.5_db))*(w_lim(thetaj)-w_lim(thetaj-1))
!!$  else
!!$     w=-(w_lim(thetaj-1)+(2.0_db*aleat2)*(w_lim(thetaj)-w_lim(thetaj-1)))
!!$  endif
!!$
!!$  ! Position azimuthale
!!$  phi = 2.0_db*pi * (real(phik)-1.0_db+aleat3)/real(n_az)
!!$
!!$  ! x et y
!!$  z=r*w
!!$  r_cos_theta = r*sqrt(1.0_db-w*w)
!!$  x=r_cos_theta*cos(phi)
!!$  y=r_cos_theta*sin(phi)

!  z = z_grid(ri,thetaj)
!  r = r_grid(ri,thetaj)
!  x = r*cos(phi)
!  y = r*sin(phi)

 ! call indice_cellule_sph(x,y,z,ri_tmp,thetaj_tmp)
 ! if (ri /= ri_tmp) then
 !    write(*,*) "Bug ri", ri, ri_tmp, sqrt(x*x+y*y)
 !    read(*,*)
 ! else if (thetaj /= thetaj_tmp) then
 !    write(*,*) "Bug zj", w, thetaj, thetaj_tmp
 !    call indice_cellule_sph(x,y,-z,ri_tmp,thetaj_tmp)
 !    write(*,*) -z, thetaj_tmp
 !    read(*,*)
 ! endif


  return

end subroutine pos_em_cellule_sph

!***********************************************************

subroutine verif_cellule_1(lambda)
! Verifie opacite cellule 1
  implicit none

  integer, intent(in) :: lambda
  integer :: n
  real :: tau

  n = 1 + floor(log(tau-1.)/log(2.))

end subroutine verif_cellule_1

!***********************************************************

subroutine angle_disque()

  implicit none

  integer :: i
  logical :: l_1etoile
  real :: r, zr, rmax, zrmax, zzmax
  
  ! test si le systeme est axisymetrique
  if (n_etoiles > 1) then
     l_1etoile=.false.
  else
     if ((abs(etoile(1)%x) > 1.0e-6).or.(abs(etoile(1)%y) > 1.0e-6).or.(abs(etoile(1)%z) > 1.0e-6))  then
        l_1etoile=.false.
     else
        l_1etoile=.true.
     endif
  endif

  if ((l_sym_axiale).and.(l_1etoile)) then
     ! On cherche le zmax / r le plus grand
     zzmax = zmax(1)
     zrmax = zzmax / rmin
     rmax = rmin
     do i = 1,n_rad
        ! On prend le rayon au bord interne de la cellule
        r= r_lim(i-1)
        zr = zmax(i) / r
        if (zr > zrmax) then
           zrmax = zr
           zzmax = zmax(i)
           rmax = r
        endif
     enddo !i
     
     ! On calcule la hauteur correspondante a rmin
     cos_max2 = rmax**2/(rmax**2+zzmax**2)
  
     ! On place le photon juste apres le bord interne (ie dans une cellule non vide)
     r_bord2 = (rmin*1.0001)**2
  else
     cos_max2 = 0.0
  endif

  return

end subroutine angle_disque

end module grid
