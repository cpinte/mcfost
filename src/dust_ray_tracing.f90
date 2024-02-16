module dust_ray_tracing
  ! RT 1 : calcul l'intensite specifique diffusee dans les directions souhaitees. : mode prefere pour SED et pour 3D
  ! RT 2 : sauvegarde l'intensite specifique et sa dependence angulaire : mode preferee pour image 2D

  use parametres
  use constantes
  use stars, only : E_stars
  use scattering
  use Temperature
  use cylindrical_grid
  use dust_prop
  !$ use omp_lib

  implicit none
  save

  real, dimension(:), allocatable :: tab_RT_incl, tab_RT_az
  real(kind=dp), dimension(:), allocatable :: tab_uv_rt, tab_w_rt
  real(kind=dp), dimension(:,:), allocatable :: tab_u_rt, tab_v_rt

  real(kind=dp), dimension(:,:), allocatable :: n_phot_envoyes

  ! Sauvegarde champ de radiation pour rt2
  integer ::  n_phi_I,  n_theta_I ! 15 et 9 ok avec 30 et 30 en mode SED
  ! Pour rt 2 : nbre d'angle de visee en azimuth. TODO : calculer automatiquement en fct de la fct de phase + interpolation
  integer :: nang_ray_tracing, nang_ray_tracing_star

  ! rt1 : dimension of specific intensity array
  integer :: n_az_rt, n_theta_rt

  ! intensite specifique
  real, dimension(:), allocatable :: J_th ! n_cells
  real, dimension(:,:,:,:,:,:), allocatable ::  xI_scatt ! 4, RT_n_incl * RT_n_az, n_cells, n_az_rt, n_theta_rt, ncpus
  real(kind=dp), dimension(:,:,:), allocatable ::  I_scatt ! 4, n_az_rt, 2

  ! methode RT 1 : saving scattered specific intensity (SED + image 3D)
  ! todo faire sauter le 2 pour gagner une dimension et rester sous la limite de 7

  integer, dimension(:,:,:), allocatable :: itheta_rt1 ! RT_n_incl,RT_n_az,nb_proc
  real(kind=dp), dimension(:,:,:), allocatable ::  sin_omega_rt1, cos_omega_rt1, sin_scatt_rt1 ! RT_n_incl,RT_n_az,nb_proc
  real(kind=dp), dimension(:,:,:,:), allocatable ::  eps_dust1 !N_type_flux, n_cells, n_az_rt,n_theta_rt

  ! methode RT 2 : saving specific intensity (image 2D)
  real, dimension(:,:,:,:,:), allocatable :: I_spec ! 4, n_theta_I, n_phi_I, n_cells, ncpus
  real, dimension(:,:), allocatable :: I_spec_star ! n_cells, ncpus

  ! Fonction source: Ok en simple
  real, dimension(:,:,:,:), allocatable ::  I_sca2 ! n_type_flux, nang_ray_tracing, 2, n_cells
  real, dimension(:,:,:,:), allocatable ::  eps_dust2 ! n_type_flux, nang_ray_tracing, 2, n_rad, nz
  real, dimension(:,:,:,:), allocatable ::  eps_dust2_star ! n_type_flux, nang_ray_tracing, 2, n_rad, nz

  real, dimension(:,:,:), allocatable ::  cos_thet_ray_tracing, omega_ray_tracing ! nang_ray_tracing, 2 (+z et -z), nb_proc
  real, dimension(:,:,:), allocatable ::  cos_thet_ray_tracing_star, omega_ray_tracing_star ! nang_ray_tracing, 2 (+z et -z), nb_proc

  real, dimension(:,:,:,:,:,:,:), allocatable :: Stokes_ray_tracing ! n_lambda, nx, ny, RT_n_incl, RT_n_az, n_type_flux, ncpus
  real, dimension(:,:,:,:), allocatable :: star_position ! n_stars, RT_n_incl, RT_n_az, 2
  real, dimension(:,:,:), allocatable :: star_vr

  real, dimension(:,:,:,:,:,:), allocatable :: tau_surface_map ! nx, ny, RT_n_incl, RT_n_az, 3, ncpus
  real, dimension(:,:,:,:,:), allocatable :: tau_map ! nx, ny, RT_n_incl, RT_n_az, ncpus
  real, dimension(:,:,:), allocatable :: stars_map ! nx, ny, 4

  contains


integer function RT2d_to_RT1d(ibin, iaz)
  ! Combines the 2d RT direction indices into 1 index
  ! RT2d_to_RT1d is between 1 and RT_n_incl * RT_n_az

  integer, intent(in) :: ibin,iaz

  RT2d_to_RT1d = ibin + RT_n_incl * (iaz-1)

  return

end function RT2d_to_RT1d

!******************************************************************************

subroutine alloc_ray_tracing()
  ! Alloue les tableaux pour le ray-tracing
  ! C. Pinte
  ! 13/10/08

  integer :: alloc_status
  real :: mem_size

  ! Ok, this is the trick which makes mcfost rt2 in 2D very fast
  ! in 2D, we save the scattered radiation field for various azimuth (as the cell is a ring) + 2 elevation directions (up & down)
  ! in 3D, each cell is already divided in azimuth, so there is no need to divide the radiation field in azymuth (it is too memory expensive anyway)

  ! Variables for rt1
  if (l3D) then
     n_az_rt = 1
     n_theta_rt = 1
  else
     n_az_rt = 45
     n_theta_rt = 2
  endif

  ! Variables for rt2
  ! n_phi_I and n_theta_I are the angular dimensions of the stored radiation field for rt2
  ! 15 15 90 90 OK pour benchmark
  n_phi_I = 15 ; ! 17/11/10 semble mieux si egal a nang_ray_tracing
  n_theta_I = 15 ;

  ! Number of angular directions were the the scattered radiation field will be calculated
  ! We use a higher direction for the stellar radiation field, as it the radiation field is directive and not 2 or 3D
  nang_ray_tracing = 15 ;
  nang_ray_tracing_star = 1000 ;

  if (lisotropic) then
     n_phi_I = 1 ;
     n_theta_I = 1 ;
     nang_ray_tracing = 1 ;
  endif

  ! datacube images
  if (lsed.and.(RT_sed_method == 1)) then
      allocate(Stokes_ray_tracing(n_lambda,1,1,RT_n_incl,RT_n_az,N_type_flux,nb_proc), stars_map(1,1,1), stat=alloc_status)
  else
     allocate(Stokes_ray_tracing(n_lambda,npix_x,npix_y,RT_n_incl,RT_n_az,N_type_flux,nb_proc), stat=alloc_status)
     if (lsepar_pola.and.llimb_darkening)then
        allocate(stars_map(npix_x,npix_y,3), stat=alloc_status)
     else
        allocate(stars_map(npix_x,npix_y,1), stat=alloc_status)
     endif
  endif
  if (alloc_status > 0) call error('Allocation error Stokes_ray_tracing')
  Stokes_ray_tracing = 0.0 ; stars_map = 0.0

  if (ltau_surface) then
     allocate(tau_surface_map(npix_x,npix_y,RT_n_incl,RT_n_az,3,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tau_surface_mao')
     tau_surface_map = 0.0
  endif

  if (ltau_map) then
     allocate(tau_map(npix_x,npix_y,RT_n_incl,RT_n_az,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error tau_map')
     tau_map = 0.0
  endif

  allocate(J_th(n_cells), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error J_th')
  J_th = 0.0_dp

  if (lscatt_ray_tracing1) then
     mem_size = (1.0*N_type_flux + 2) * RT_n_incl * RT_n_az * n_cells * n_az_rt * n_theta_rt * nb_proc * 4 / 1024.**3
     if (mem_size > 0.5) write(*,*) "Trying to allocate", mem_size, "GB for ray-tracing"

     allocate(xI_scatt(n_az_rt,n_theta_rt,N_type_flux,RT_n_incl*RT_n_az,n_cells,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error xI_scatt')
     xI_scatt = 0.0_dp

     allocate(I_scatt(n_az_rt,n_theta_rt,N_type_flux), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error I_scatt')
     I_scatt = 0.0_dp

     allocate(eps_dust1(n_az_rt,n_theta_rt,N_type_flux,n_cells), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error eps_dust1')
     eps_dust1 =0._dp

     allocate(itheta_rt1(RT_n_incl,RT_n_az,nb_proc), sin_scatt_rt1(RT_n_incl,RT_n_az,nb_proc), &
          sin_omega_rt1(RT_n_incl,RT_n_az,nb_proc), cos_omega_rt1(RT_n_incl,RT_n_az,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error itheta_rt1')
     itheta_rt1 = 1
     sin_scatt_rt1 = 0.0
     sin_omega_rt1 = 0.5
     cos_omega_rt1 = 0.5

  else !rt 2
     allocate(I_spec_star(n_cells,nb_proc), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error I_spec_star')
     I_spec_star = 0.0_dp

     allocate(eps_dust2(N_type_flux,nang_ray_tracing,0:1,n_cells), &
          I_sca2(N_type_flux,nang_ray_tracing,0:1,n_cells), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error eps_dust2')
     eps_dust2 =0.0_dp ; I_sca2 = 0.0_dp ;


     allocate(eps_dust2_star(n_Stokes,nang_ray_tracing_star,0:1,n_cells), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error eps_dust2_star')
     eps_dust2_star =0.0_dp


     allocate(cos_thet_ray_tracing(nang_ray_tracing,0:1,nb_proc), omega_ray_tracing(nang_ray_tracing,0:1,nb_proc),&
          stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error cos_thet_ray_tracing')
     cos_thet_ray_tracing = 0._dp
     omega_ray_tracing = 0._dp

     allocate(cos_thet_ray_tracing_star(nang_ray_tracing_star,0:1,nb_proc), &
          omega_ray_tracing_star(nang_ray_tracing_star,0:1,nb_proc),&
          stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error cos_thet_ray_tracing_star')
     cos_thet_ray_tracing_star = 0._dp
     omega_ray_tracing_star = 0._dp

     allocate(I_spec(N_type_flux,n_theta_I,n_phi_I,n_cells,nb_proc),stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error xI')
     I_spec = 0.0_dp
  endif !lscatt_ray_tracing1

  return

end subroutine alloc_ray_tracing

!***********************************************************

subroutine dealloc_ray_tracing()
  ! desalloue les tableaux pour le ray-tracing
  ! C. Pinte
  ! 13/10/08

  if (allocated(Stokes_ray_tracing)) then
     deallocate(Stokes_ray_tracing,stars_map,J_th)

     if (lscatt_ray_tracing1) then
        deallocate(xI_scatt,I_scatt,eps_dust1,sin_scatt_rt1,sin_omega_rt1,cos_omega_rt1)
     else
        deallocate(I_spec_star,eps_dust2, eps_dust2_star, I_sca2,cos_thet_ray_tracing,omega_ray_tracing,I_spec)
     endif
     deallocate(tab_RT_incl,tab_RT_az,tab_uv_rt,tab_u_rt,tab_v_rt,tab_w_rt)
  endif

  return

end subroutine dealloc_ray_tracing

!***********************************************************

subroutine init_directions_ray_tracing()
  ! initialise les directions (inclinaisons) pour le ray-tracing
  ! C. Pinte
  ! 09/09/08

  real(kind=dp) :: cos_min, cos_max
  integer :: ibin, iaz, alloc_status

  alloc_status = 0

  if (.not.allocated(tab_RT_incl)) then
     allocate(tab_RT_incl(RT_n_incl),tab_RT_az(RT_n_az), tab_uv_rt(RT_n_incl), &
          tab_u_rt(RT_n_incl,RT_n_az),tab_v_rt(RT_n_incl,RT_n_az),tab_w_rt(RT_n_incl))
  endif

  if (RT_n_incl==1) then
     tab_RT_incl(1) = RT_imin
  else
     cos_min = cos(RT_imin * deg_to_rad)
     cos_max = cos(RT_imax * deg_to_rad)

     if (lRT_i_centered) then
        do ibin=1, RT_n_incl
           tab_RT_incl(ibin) = acos(  cos_min + (real(ibin) -0.5)/real(RT_n_incl) * (cos_max - cos_min) ) / deg_to_rad
        enddo
     else
        do ibin=1, RT_n_incl
           tab_RT_incl(ibin) = acos(  cos_min + (real(ibin) - 1)/(real(RT_n_incl) -1) * (cos_max - cos_min) ) / deg_to_rad
        enddo
     endif
  endif

  if (RT_n_az==1) then
     tab_RT_az(1) = RT_az_min
  else
     do iaz=1, RT_n_az
        tab_RT_az(iaz) = (RT_az_min +  (real(iaz) - 1)/(real(RT_n_az) -1) * (RT_az_max - RT_az_min) )
     enddo
  endif


  do ibin=1, RT_n_incl
     ! 0 est remplace par un epsilon car il faut donner un axe de reference
     ! pour les differentes directions de ray-tracing utilisees dans le RT2
     tab_uv_rt(ibin) = sin( sign(max(abs(tab_RT_incl(ibin)),1e-20), tab_RT_incl(ibin)) *deg_to_rad ) ! uv_rt mais v_rt = 0 ici
     tab_w_rt(ibin) = sqrt(1.0_dp - tab_uv_rt(ibin)*tab_uv_rt(ibin))
     if (abs(tab_RT_incl(ibin)) > 90.) tab_w_rt(ibin) = - tab_w_rt(ibin)

     ! phi = 0 correspond a -v
     do iaz =1, RT_n_az
        tab_u_rt(ibin,iaz) =   tab_uv_rt(ibin) * sin(tab_RT_az(iaz)*deg_to_rad)
        tab_v_rt(ibin,iaz) =   - tab_uv_rt(ibin) * cos(tab_RT_az(iaz)*deg_to_rad)
     enddo
  enddo

  if (.not. allocated(star_position)) then
     allocate(star_position(n_etoiles,RT_n_incl,RT_n_az,2), star_vr(n_etoiles,RT_n_incl,RT_n_az), stat=alloc_status)
     if (alloc_status > 0) call error('Allocation error star_position')
  endif

  return

end subroutine init_directions_ray_tracing

!***********************************************************

subroutine angles_scatt_rt2(id,ibin,x,y,z,u,v,w,lstar)
  ! Calcul les cosinus des angles de diffusion vers toutes les directions
  ! utilisees pour le ray tracing
  ! Calcul aussi les matrices de rotation pour le vecteur de Stokes
  ! Remplit les tableaux cos_thet_ray_tracing et omega_ray_tracing
  ! RT2
  ! C. Pinte
  ! 19/01/08

  ! NB: cos_thet_ray_tracing n'est plus utilise  pour la composante diffusee
  ! plusieurs fois depuis le super-echantillonage

  integer, intent(in) :: id, ibin ! RT2 : la direction d'azimuth n'a pas d'importance ici
  real(kind=dp), intent(in) :: x,y,z, u, v, w  ! z est inutilise
  logical, intent(in) :: lstar
  real(kind=dp) :: uv0, w0, phi, phi_pos,u_ray_tracing,  v_ray_tracing, w_ray_tracing, prod1, prod2
  real(kind=dp) :: w2, v1pi, v1pj, v1pk, xnyp, costhet, theta, omega
  integer :: iscatt, direction, N

  ! Direction observateur dans repere refence
  uv0 = tab_uv_rt(ibin) ; w0 = tab_w_rt(ibin) ! epsilon needed here

  ! Calcul angle phi a la position
  phi_pos = modulo(atan2(x,y) + deux_pi,deux_pi)

  ! Nouvelle orientation: les directions de ray-tracing sont definies par rapport a l'axe de ref y
  ! il faut se replacer dans ce cas en fct de la position
  ! ie il faut soustraire le phi de la position du paquet

  if (lstar) then
     N = nang_ray_tracing_star
  else
     N = nang_ray_tracing
  endif

  ! creation des tableaux de u,v,w
  do iscatt=1, N
     ! if (l_sym_ima) then
     !    phi = pi * real(iscatt-1) / real(nang_ray_tracing-1) + pi
     ! else
     phi = deux_pi * real(iscatt) / real(N)
     ! endif

     phi = phi - phi_pos

     u_ray_tracing = uv0 * sin(phi)
     v_ray_tracing = - uv0 * cos(phi)
     w_ray_tracing = w0

     prod1 = u_ray_tracing * u  + v_ray_tracing * v

     do direction=0,1
        if (direction==1) then
           w2 = w
        else
           w2 = -w
        endif

        prod2 = w_ray_tracing * w2
        if (lstar) then
           cos_thet_ray_tracing_star(iscatt,direction,id) =  prod1 + prod2
        else
           cos_thet_ray_tracing(iscatt,direction,id) =  prod1 + prod2
        endif

        if (lsepar_pola) then
           call rotation(u,v,w2,-u_ray_tracing,-v_ray_tracing,-w_ray_tracing,v1pi,v1pj,v1pk)

           xnyp = sqrt(v1pk*v1pk + v1pj*v1pj)
           if (xnyp < 1e-10) then
              xnyp = 0.0
              costhet = 1.0_dp ! omega=0
           else
              costhet = v1pj / xnyp
           endif

           ! calcul de l'angle entre la normale et l'axe z (theta)
           theta = acos(costhet)
           if (theta >= pi) theta = 0.0

           !     le plan de diffusion est a +ou- 90deg de la normale
           ! Ne doit pas etre utilise en ray-tracing !!! teste OK sans
           !theta = theta  + pi_sur_deux

           !----dans les matrices de rotation l'angle est omega = 2 * theta-----
           omega = 2.0_dp * theta ! A verifier: le moins est pour corriger un bug de signe trouve par Marshall (supprime)
           !     prochain if car l'arccos va de 0 a pi seulement
           !     le +/- pour faire la difference dans le sens de rotation
           if (v1pk < 0.0) omega = -1.0_dp * omega

           if (lstar) then
              omega_ray_tracing_star(iscatt,direction,id) = omega
           else
              omega_ray_tracing(iscatt,direction,id) = omega
           endif
        endif
     enddo ! direction
  enddo !iscatt

  return

end subroutine angles_scatt_rt2

!***********************************************************

subroutine angles_scatt_rt1(id,u,v,w)
  ! Calcul le (ou les) angle(s) de diffusion pour le RT1
  ! dans une direction de diffusion a 1 position
  ! lance 1 seule pour unr direction de paquet
  !
  ! renvoie
  ! - itheta : indice de l'angle de diffusion
  ! - cosw et sinw, cos et sin de l'angle de la matrice de rotation
  ! qui sont utilises par new_stokes_rt1
  ! C. Pinte
  ! 13/09/09, version intiale  19/01/08

  integer, intent(in) :: id
  real(kind=dp), intent(in) :: u, v, w

  real(kind=dp) :: v1pi, v1pj, v1pk, xnyp, costhet, theta, omega, cosw, sinw
  real :: cos_scatt
  integer :: k, ibin, iaz

  ! Direction observateur dans repere refence
  do ibin=1,RT_n_incl
     do iaz=1, RT_n_az
        cos_scatt = tab_u_rt(ibin,iaz) * u +  tab_v_rt(ibin,iaz) * v + tab_w_rt(ibin) * w

        k = nint(acos(cos_scatt) * real(nang_scatt)/pi) ! TODO : + 1 pour test
        if (k > nang_scatt) k = nang_scatt
        if (k < 1) k = 1
        itheta_rt1(ibin,iaz,id) = k
        sin_scatt_rt1(ibin,iaz,id) = sqrt(1.0 - cos_scatt**2)

        if (lsepar_pola) then
           ! Matrice de rotation pour se recaller sur le Nord celeste
           call rotation(u,v,w,-tab_u_rt(ibin,iaz),-tab_v_rt(ibin,iaz),-tab_w_rt(ibin),v1pi,v1pj,v1pk)
           xnyp = sqrt(v1pk*v1pk + v1pj*v1pj)
           if (xnyp < 1e-10) then
              xnyp = 0.0_dp
              costhet = 1.0_dp ! omega = 0
           else
              costhet = -1.0_dp*v1pj / xnyp
           endif

           ! calcul de l'angle entre la normale et l'axe z (theta)
           theta = acos(costhet)
           if (theta >= pi) theta = 0.0_dp

           !     le plan de diffusion est a +ou- 90deg de la normale
           theta = theta + pi_sur_deux

           !----dans les matrices de rotation l'angle est omega = 2 * theta-----
           omega = 2.0_dp * theta
           !     prochain if car l'arccos va de 0 a pi seulement
           !     le +/- pour faire la difference dans le sens de rotation
           if (v1pk < 0.0) omega = -1.0_dp * omega

           cosw = cos(omega)
           sinw = sin(omega)
           if (abs(cosw) < 1e-06) cosw = 0.0_dp
           if (abs(sinw) < 1e-06) sinw = 0.0_dp

           cos_omega_rt1(ibin,iaz,id) = cosw
           sin_omega_rt1(ibin,iaz,id) = sinw
        endif ! lsepar_pola
     enddo ! iaz
  enddo ! ibin

  return

end subroutine angles_scatt_rt1

!***********************************************************

subroutine calc_xI_scatt(id,lambda,p_lambda,icell, phik,psup,l,stokes,flag_star)
  ! RT1
  ! Calcul les matrices de Mueller pour une direction de diffusion
  ! lance dans chaque cellule traversee
  ! utilise les resultats de angles_scatt_rt1
  ! C. Pinte
  ! 13/09/09, version intiale  19/01/08
  ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
  ! 20/04/2023

  implicit none
  ! pour chaque cellule stocke le champ diffusee par un paquet dans les directions
  ! de ray-tracing ibin et dir en fonction de atan2(x,y) : ie un paquet
  ! n'allume qu'un seul endroit de l'anneau au lieu de partout avant,
  ! mais la taille des tableaux est la meme
  ! pour le ray-tracing, il ne faut plus interpoler entre les directions de ray-tracing
  ! mais entre les positions dans le disque

  real(kind=dp), intent(in) :: Stokes, l
  logical, intent(in) :: flag_star
  integer, intent(in) :: id, lambda, p_lambda, icell, phik, psup

  real(kind=dp) :: flux
  integer :: ibin, iaz, iRT, it, p_icell

  if (lvariable_dust) then
     p_icell = icell
  else
     p_icell = icell_ref
  endif

  do ibin = 1, RT_n_incl
     do iaz=1, RT_n_az
        it = itheta_rt1(ibin,iaz,id)
        if (lmueller) then
           flux = l * stokes * tab_mueller_pos(1,1,it,p_icell,p_lambda)
        else
           flux = l * stokes * tab_s11_pos(it,p_icell,p_lambda) !* sin_scatt_rt1(ibin,id)
        endif
        ! TODO : est-ce qu'il ne faut pas moyenner par le sin de l'angle de scatt dans la cellule azimuthale ???

        iRT = RT2d_to_RT1d(ibin, iaz)
        xI_scatt(phik,psup,1,iRT,icell,id) =  xI_scatt(phik,psup,1,iRT,icell,id) + flux

        if (lsepar_contrib) then
           if (flag_star) then
              xI_scatt(phik,psup,n_Stokes+2,iRT,icell,id) =  xI_scatt(phik,psup,n_Stokes+2,iRT,icell,id) + flux
           else
              xI_scatt(phik,psup,n_Stokes+4,iRT,icell,id) =  xI_scatt(phik,psup,n_Stokes+4,iRT,icell,id) + flux
           endif
        endif
     enddo ! iaz
  enddo ! ibin

  return

end subroutine calc_xI_scatt

!***********************************************************

subroutine calc_xI_scatt_pola(id,lambda,p_lambda,icell,phik,psup,l,stokes,flag_star)
  ! RT1
  ! Calcul les matrices de Mueller pour une direction de diffusion
  ! lance dans chaque cellule traversee
  ! utilise les resultats de angles_scatt_rt1
  ! C. Pinte
  ! 13/09/09, version intiale  19/01/08
  ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
  ! 20/04/2023

  ! TODO : a mettre a jour pour RT 3D !

  implicit none

  real(kind=dp), dimension(4), intent(in) :: Stokes
  real(kind=dp), intent(in) :: l
  logical, intent(in) :: flag_star
  integer, intent(in) :: id, lambda, p_lambda, icell, phik, psup

  real(kind=dp), dimension(4) :: C, D, S
  real(kind=dp), dimension(4,4) ::  M, ROP, RPO
  real(kind=dp) :: cosw, sinw, flux
  real :: s11, s12, s33, s34
  integer :: ibin, iaz, iRT, it, p_icell


  ROP = 0.0_dp
  RPO = 0.0_dp
  RPO(1,1) = 1.0_dp
  ROP(1,1) = 1.0_dp
  RPO(4,4) = 1.0_dp
  ROP(4,4) = 1.0_dp

  M = 0.0_dp

  if (lvariable_dust) then
     p_icell = icell
  else
     p_icell = icell_ref
  endif

  do ibin = 1, RT_n_incl
     do iaz = 1, RT_n_az
        ! Matrice de Mueller
        it = itheta_rt1(ibin,iaz,id)

        if (lmueller) then
           s11 = tab_mueller_pos(1,1,it,p_icell,p_lambda)
           M(:,:) = tab_mueller_pos(:,:,it,p_icell,p_lambda)*(-s11)
           M(:,4) = -M(:,4)
           M(4,:) = -M(4,:)
           M(2,2) = -M(2,2)
           M(1,1) = s11
           !On passe à l'opposé la quatrième collonne et ligne à cause de la convention
           !Le - M(2,2) est pour être cohérent avec les opérations effectuées dans le cas mie

        else

           s11 = tab_s11_pos(it,p_icell,p_lambda)
           s12 = -s11 * tab_s12_o_s11_pos(it,p_icell,p_lambda)
           s33 = -s11 * tab_s33_o_s11_pos(it,p_icell,p_lambda)
           s34 = -s11 * tab_s34_o_s11_pos(it,p_icell,p_lambda)
           M(1,1) = s11 ; M(2,2) = s11 ; M(1,2) = s12 ; M(2,1) = s12
           M(3,3) = s33 ; M(4,4) = s33 ; M(3,4) = -s34 ; M(4,3) = s34

        endif

        ! Matrices de rotation
        cosw = cos_omega_rt1(ibin,iaz,id)
        sinw = sin_omega_rt1(ibin,iaz,id)

        RPO(2,2) = -cosw
        ROP(2,2) = cosw
        RPO(2,3) = -sinw
        ROP(2,3) = -sinw
        RPO(3,2) = -sinw
        ROP(3,2) = sinw
        RPO(3,3) = cosw
        ROP(3,3) = cosw

        !  STOKE FINAL = RPO * M * ROP * STOKE INITIAL

        if (lmueller) then
           ! 1ere rotation
           C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
           C(1)=stokes(1)
           C(4)=stokes(4)

           ! multiplication matrice Mueller
           D=matmul(M,C)

           ! 2nde rotation
           S(2:3)=matmul(RPO(2:3,2:3),D(2:3))
           S(1)=D(1)
           S(4)=D(4)

        else
           ! 1ere rotation
           C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
           C(1)=stokes(1)
           C(4)=stokes(4)

           ! multiplication matrice Mueller par bloc
           D(1:2)=matmul(M(1:2,1:2),C(1:2))
           D(3:4)=matmul(M(3:4,3:4),C(3:4))

           ! 2nde rotation
           S(2:3)=matmul(RPO(2:3,2:3),D(2:3))
           S(1)=D(1)
           S(4)=D(4)
        endif

        iRT = RT2d_to_RT1d(ibin, iaz)
        xI_scatt(phik,psup,1:4,iRT,icell,id) =  xI_scatt(phik,psup,1:4,iRT,icell,id) + l * S(:)

        if (lsepar_contrib) then
           flux = l * S(1)
           if (flag_star) then
              !n_Stokes+2 = 6
              xI_scatt(phik,psup,6,iRT,icell,id) =  xI_scatt(phik,psup,6,iRT,icell,id) + flux
           else
              !n_Stokes+4 = 8
              xI_scatt(phik,psup,8,iRT,icell,id) =  xI_scatt(phik,psup,8,iRT,icell,id) + flux
           endif
        endif
     enddo ! iaz
  enddo ! ibin

  return

end subroutine calc_xI_scatt_pola

!***********************************************************

subroutine init_dust_source_fct1(lambda,ibin,iaz)
  ! RT1

  implicit none

  integer, intent(in) :: lambda, ibin, iaz

  integer :: itype, iRT
  integer :: icell, p_icell
  real(kind=dp) :: facteur, energie_photon, n_photons_envoyes, kappa_sca, kappa_ext

  iRT = RT2d_to_RT1d(ibin, iaz)

  eps_dust1(:,:,:,:) = 0.0_dp

  ! Contribution emission thermique directe
  call calc_Jth(lambda)
  !unpolarized_Stokes = 0. ;   unpolarized_Stokes(1) = 1.

  ! TODO : la taille de eps_dust1 est le facteur limitant pour le temps de calcul

  if (lmono0) then ! image
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (real(nbre_photons_loop)*real(nbre_photons_image) *  AU_to_cm * pi)
  else ! SED
     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (n_photons_envoyes *  AU_to_cm * pi) !lambda.F_lambda
  endif

  !write(*,*) "TATA", lambda, kappa(1,lambda), tab_albedo_pos(1,lambda)

  ! Intensite specifique diffusion
  !$omp parallel &
  !$omp default(none) &
  !$omp private(icell,p_icell,facteur,kappa_sca,kappa_ext,itype,I_scatt) &
  !$omp shared(n_cells,energie_photon,volume,n_az_rt,n_theta_rt,kappa,kappa_factor,N_type_flux,lambda,iRT) &
  !$omp shared(J_th,xI_scatt,eps_dust1,lsepar_pola,lsepar_contrib,n_Stokes,nb_proc,tab_albedo_pos,lvariable_dust,icell_ref)

  p_icell = icell_ref

  !$omp do schedule(static,n_cells/nb_proc)
  do icell=1, n_cells
     if (lvariable_dust) p_icell = icell
     facteur = energie_photon / volume(icell) * n_az_rt * n_theta_rt ! n_az_rt * n_theta_rt car subdivision virtuelle des cellules
     kappa_ext = kappa(p_icell,lambda) * kappa_factor(icell)
     kappa_sca = kappa_ext * tab_albedo_pos(p_icell,lambda)

     !write(*,*) icell,  kappa_sca, kappa(icell,lambda), tab_albedo_pos(icell,lambda)

     ! TODO : les lignes suivantes sont tres cheres en OpenMP
     if (kappa_ext > tiny_dp) then
        do itype=1,N_type_flux
           I_scatt(:,:,itype) = sum(xI_scatt(:,:,itype,iRT,icell,:),dim=3) * facteur * kappa_sca
        enddo ! itype

        eps_dust1(:,:,1,icell) =  (I_scatt(:,:,1) +  J_th(icell)) / kappa_ext

        ! TODO : there is might a problem in polarization in 2D
        if (lsepar_pola) then
           eps_dust1(:,:,2:4,icell) =  I_scatt(:,:,2:4) / kappa_ext
        endif

        if (lsepar_contrib) then
           eps_dust1(:,:,n_Stokes+2,icell) = I_scatt(:,:,n_Stokes+2) / kappa_ext
           eps_dust1(:,:,n_Stokes+3,icell) = J_th(icell) / kappa_ext
           eps_dust1(:,:,n_Stokes+4,icell) = I_scatt(:,:,n_Stokes+4) / kappa_ext
        endif ! lsepar_contrib
     else
        eps_dust1(:,:,:,icell) = 0.0_dp
     endif ! kappa > 0
  enddo ! icell
  !$omp enddo
  !$omp end parallel

  return

end subroutine init_dust_source_fct1

!***********************************************************

subroutine init_dust_source_fct2(lambda,p_lambda,ibin)
  ! RT2
  ! calcule la fct source pour integrer par ray-tracing
  ! Cette version ne marche qu'en 2D
  ! C. Pinte
  ! 25/09/08


  integer, intent(in) :: lambda, p_lambda,ibin
  integer :: iscatt, dir, icell, p_icell
  real(kind=dp) :: kappa_ext

  real :: Q, U

  if (lmono0) then
     write(*,*) "i=", tab_RT_incl(ibin)
     write(*,*) "Scattered specific intensity ..."
  endif

  I_sca2 = 0.0_dp
  eps_dust2_star = 0.0_dp
  eps_dust2 = 0.0_dp

  ! Ajout du champ de radiation stellaire diffuse 1 seule fois
  call calc_Isca_rt2_star(lambda, p_lambda, ibin) ! this is the slow line  : 1.6s sur 2.4s  ---> 0.54s apres parallelization

  ! Contribution lumiere diffusee (y compris multiple et thermique diffusee)
  call calc_Isca_rt2(lambda, p_lambda, ibin)  ! ~ 0.7s sur 2.4s

  ! Contribution emission thermique directe
  call calc_Jth(lambda)   !---> 0.05s sur 2.4s

  ! Fonction source, indices : pola, iscatt, dir, i, j
  !$omp parallel &
  !$omp default(none) &
  !$omp shared(lambda,n_cells,nang_ray_tracing,eps_dust2,I_sca2,J_th,kappa,kappa_factor,eps_dust2_star,lsepar_pola) &
  !$omp shared(lsepar_contrib,n_Stokes,nang_ray_tracing_star,nb_proc,tab_albedo_pos,icell_ref,lvariable_dust) &
  !$omp private(icell,p_icell,dir,iscatt,Q,U,kappa_ext)
  p_icell = icell_ref
  !$omp do schedule(static,n_cells/nb_proc)
  do icell=1, n_cells
     if (lvariable_dust) p_icell = icell

     kappa_ext = kappa(p_icell,lambda) * kappa_factor(icell)
     if (kappa_ext > tiny_dp) then
        ! Boucle sur les directions de ray-tracing
        do dir=0,1
           do iscatt = 1, nang_ray_tracing
              eps_dust2(1,iscatt,dir,icell) =  ( I_sca2(1,iscatt,dir,icell)  +  J_th(icell) ) / kappa_ext

              if (lsepar_pola) then
                 eps_dust2(2:4,iscatt,dir,icell) =  I_sca2(2:4,iscatt,dir,icell)  / kappa_ext

                 ! Transforming to P*I and 2*theta for new interpolation scheme in dust_source_fct()
                 Q = eps_dust2(2,iscatt,dir,icell) ; U = eps_dust2(3,iscatt,dir,icell)
                 eps_dust2(2,iscatt,dir,icell) = sqrt(Q**2 + U**2)
                 eps_dust2(3,iscatt,dir,icell) = atan2(U,Q)
              endif

              if (lsepar_contrib) then
                 eps_dust2(n_Stokes+2,iscatt,dir,icell) =    I_sca2(n_Stokes+2,iscatt,dir,icell) / kappa_ext
                 eps_dust2(n_Stokes+3,iscatt,dir,icell) =    J_th(icell) / kappa_ext
                 eps_dust2(n_Stokes+4,iscatt,dir,icell) =    I_sca2(n_Stokes+4,iscatt,dir,icell) / kappa_ext
              endif ! lsepar_contrib
           enddo ! iscatt

           do iscatt = 1, nang_ray_tracing_star
              eps_dust2_star(:,iscatt,dir,icell) = eps_dust2_star(:,iscatt,dir,icell) / kappa_ext

              if (lsepar_pola) then
                 ! Transforming to P*I and 2*theta for new interpolation scheme in dust_source_fct()
                 Q = eps_dust2_star(2,iscatt,dir,icell) ; U = eps_dust2_star(3,iscatt,dir,icell)
                 eps_dust2_star(2,iscatt,dir,icell) = sqrt(Q**2 + U**2)
                 eps_dust2_star(3,iscatt,dir,icell) = atan2(U,Q)
              endif
           enddo ! iscatt
        enddo ! dir
     else
        eps_dust2(:,:,:,icell) = 0.0_dp
        eps_dust2_star(:,:,:,icell) = 0.0_dp
     endif
  enddo !icell
  !$omp end do
  !$omp end parallel

  if (lmono0) write(*,*)  "Done"

  return

end subroutine init_dust_source_fct2

!***********************************************************

subroutine calc_Jth(lambda)
  ! calcul emissivite thermique
  ! C. Pinte
  ! 25/09/08

  integer, intent(in) :: lambda
  integer :: l, T, icell, p_icell
  real(kind=dp) ::  Temp, cst_wl, wl, coeff_exp, cst_E

  ! longueur d'onde en metre
  wl = tab_lambda(lambda)*1.e-6

  J_th(:) = 0.0


  ! Intensite specifique emission thermique
  if ((l_em_disk_image).or.(lsed)) then
     if (lRE_LTE) then
        cst_E=2.0*hp*c_light**2
        !$omp parallel &
        !$omp default(none) &
        !$omp shared(n_cells,Tdust,wl,lambda,kappa_abs_LTE,kappa_factor,cst_E,J_th) &
        !$omp shared(lvariable_dust,icell_ref) &
        !$omp private(icell,p_icell,Temp,cst_wl,coeff_exp)
        p_icell = icell_ref
        !$omp do
        do icell=1, n_cells
           if (lvariable_dust) p_icell = icell
           Temp=Tdust(icell) ! que LTE pour le moment
           if (Temp*wl > 3.e-4) then
              cst_wl=cst_th/(Temp*wl)
              coeff_exp=exp(cst_wl)
              J_th(icell) = cst_E/((wl**5)*(coeff_exp-1.0)) * wl * kappa_abs_LTE(p_icell,lambda) * kappa_factor(icell) ! Teste OK en mode SED avec echantillonnage lineaire du plan image
           else
              J_th(icell) = 0.0_dp
           endif
        enddo ! icell
        !$omp end do
        !$omp end parallel
     endif !lRE_LTE

     if (lRE_nLTE) then
        cst_E=2.0*hp*c_light**2
        do icell=1,n_cells
           do l=grain_RE_nLTE_start,grain_RE_nLTE_end
              Temp=Tdust_1grain(l,icell)
              if (Temp*wl > 3.e-4) then
                 cst_wl=cst_th/(Temp*wl)
                 coeff_exp=exp(cst_wl)
                 J_th(icell) = J_th(icell) + cst_E/((wl**5)*(coeff_exp-1.0)) * wl * &
                      C_abs_norm(l,lambda) * densite_pouss(l,icell)
              endif
           enddo ! l
        enddo ! icell
     endif !lRE_nLTE

     if (lnRE) then
        do icell=1,n_cells
           do l=grain_nRE_start,grain_nRE_end
              if (l_RE(l,icell)) then ! le grain a une temperature
                 Temp=Tdust_1grain_nRE(l,icell) ! WARNING : TODO : this does not work in 3D
                 if (Temp*wl > 3.e-4) then
                    cst_wl=cst_th/(Temp*wl)
                    coeff_exp=exp(cst_wl)
                    J_th(icell) = J_th(icell) + cst_E/((wl**5)*(coeff_exp-1.0)) * wl * &
                         C_abs_norm(l,lambda)*densite_pouss(l,icell)
                 endif !cst_wl
              else ! ! la grain a une proba de T
                 do T=1,n_T
                    temp=tab_Temp(T)
                    if (Temp*wl > 3.e-4) then
                       cst_wl=cst_th/(Temp*wl)
                       coeff_exp=exp(cst_wl)
                       J_th(icell) = J_th(icell) + cst_E/((wl**5)*(coeff_exp-1.0)) * wl * &
                            C_abs_norm(l,lambda)*densite_pouss(l,icell) * Proba_Tdust(T,l,icell)
                    endif !cst_wl
                 enddo ! T
              endif ! l_RE
           enddo ! l
        enddo ! icell
     endif !lnRE

  endif ! lsed or lemission_disk

  return

end subroutine calc_Jth

!***********************************************************

subroutine calc_Isca_rt2(lambda,p_lambda,ibin)
  ! Version acceleree
  ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
  ! 20/04/2023

  integer, intent(in) :: lambda, p_lambda, ibin

  real(kind=dp), dimension(4) :: Stokes, S, C, D
  real(kind=dp) :: x, y, z, u, v, w, phi, w02, facteur, energie_photon, kappa_sca
  real(kind=dp) :: u_ray_tracing, v_ray_tracing, w_ray_tracing, uv0, w0

  real(kind=dp), dimension(:,:,:,:), allocatable :: Inu ! sum of I_spec over cpus
  real(kind=dp), dimension(4,4) ::  M, ROP, RPO

  integer :: theta_I, phi_I, dir, iscatt, id
  integer :: k, alloc_status, i1, i2, icell, p_icell
  real :: cos_scatt, sum_sin, f1, f2, sin_scatt, phi_scatt
  real(kind=dp) :: omega, sinw, cosw, n_photons_envoyes, v1pi, v1pj, v1pk, xnyp, costhet, theta

  integer, parameter :: N_super = 5 ! number of elements (N_super^2) to average the phase function
  ! 5 cree un leger surcout dans le cas avec strat (qq 10 sec par inclinaison et par lambda)
  ! 15 cree un important surcout

  real :: s11, s12, s33, s34

  ! Many dimensions but small numbers (<1MB for default values)
  real, dimension(N_super,N_super,n_theta_I,n_phi_I,nang_ray_tracing,0:1) :: tab_sin_scatt_norm
  integer, dimension(N_super,N_super,n_theta_I,n_phi_I,nang_ray_tracing,0:1) :: tab_k

  real, dimension(n_theta_I,n_phi_I,nang_ray_tracing,0:1) :: s11_save
  real(kind=dp), dimension(n_theta_I,n_phi_I,nang_ray_tracing,0:1) :: tab_sinw, tab_cosw

  ! Direction observateur dans repere refence
  uv0 = tab_uv_rt(ibin) ; w0 = tab_w_rt(ibin) ! epsilon needed here

  ! Allocation dynamique pour passer en stack
  ! Todo : we do not need to compute this for each ibin !!!
  allocate(Inu(N_type_flux,n_theta_I,n_phi_I,n_cells), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error I_nu')
  Inu = 0.0_dp

  id = 1

  if (lmono0) then ! image
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (real(nbre_photons_loop)*real(nbre_photons_image) *  AU_to_cm * pi)
  else ! SED
     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (n_photons_envoyes *  AU_to_cm * pi)
  endif

  ! Position : virtuelle !!
  x = 0.0_dp
  y = 1.0_dp ! doit juste etre non nul
  z = 1.0_dp ! sert a rien

  ! Champ de radiation
  do icell=1, n_cells
     do phi_I=1,n_phi_I
        do theta_I=1,n_theta_I
           Inu(:,theta_I,phi_I,icell) = sum(I_spec(:,theta_I,phi_I,icell,:),dim=2)
        enddo ! phi_I
     enddo ! theta_I
  enddo !icell


  ! Precalcul des directions ou on va calculer Inu * s11
  ! Permet d'etre aussi rapide que la version initiale calc_Isca2

  ! TODO : tab_k, tab_sin_scatt_norm, tab_cosw and tab_cosw can be saved
  ! They only to be compured once per ibin

  do dir=0,1
     do iscatt = 1, nang_ray_tracing
        phi_scatt = deux_pi * real(iscatt) / real(nang_ray_tracing)  ! todo : precalculate this section

        ! Direction of the various rays that will sampled a "cylindrical" cell
        u_ray_tracing = uv0 * sin(phi_scatt)
        v_ray_tracing = - uv0 * cos(phi_scatt)
        w_ray_tracing = w0

        ! Direction of the stored specific intensity
        do phi_I=1,n_phi_I
           do theta_I=1,n_theta_I
              ! Moyennage de la fct de phase sur le  bin
              sum_sin = 0.
              do i2=1, N_super
                 do i1 = 1, N_super
                    ! direction de vol moyenne du bin
                    f1 = real(i1) / (N_super + 1)
                    f2 = real(i2) / (N_super + 1)

                    w = (2.0_dp * ((real(theta_I,kind=dp) - f1) / real(n_theta_I,kind=dp) ) - 1.0_dp) * (2*dir-1)
                    phi = deux_pi * (real(phi_I,kind=dp) - f2) / real(n_phi_I,kind=dp)

                    w02 = sqrt(1.0_dp-w*w)
                    u = w02 * sin(phi)
                    v = - w02 * cos(phi)
                    ! BUG : u_ray_tracing

                    ! Angle de diffusion --> n'utilise plus cos_thet_ray_tracing depuis super-echantillonage
                    cos_scatt = u_ray_tracing * u + v_ray_tracing * v + w_ray_tracing * w

                    k = nint(acos(cos_scatt) * real(nang_scatt)/pi)
                    if (k > nang_scatt) k = nang_scatt
                    if (k < 0) k = 0
                    tab_k(i1,i2,theta_I,phi_I,iscatt,dir) = k

                    sin_scatt = sqrt(1.0_dp - cos_scatt*cos_scatt)
                    tab_sin_scatt_norm(i1,i2,theta_I,phi_I,iscatt,dir) = sin_scatt

                    sum_sin = sum_sin + sin_scatt
                 enddo !i1
              enddo !i2

              ! Normalization du facteur sin ici
              ! tab_sin_scatt depends on ibin !!!
              tab_sin_scatt_norm(:,:,theta_I,phi_I,iscatt,dir) = tab_sin_scatt_norm(:,:,theta_I,phi_I,iscatt,dir) / sum_sin

              if (lsepar_pola) then ! On calcule les s12, s33, s34 et la matrice de rotation

                 ! On prend le milieu du bin uniquement pour la pola, car les fct sont plus smooth
                 w = (2.0_dp * ((real(theta_I,kind=dp) - 0.5) / real(n_theta_I,kind=dp) ) - 1.0_dp) * (2*dir-1)
                 phi = deux_pi * (real(phi_I,kind=dp) - 0.5) / real(n_phi_I,kind=dp)

                 w02 = sqrt(1.0_dp-w*w)
                 u = w02 * sin(phi)
                 v = - w02 * cos(phi)

                 ! Angle de diffusion
                 cos_scatt = u_ray_tracing * u + v_ray_tracing * v + w_ray_tracing * w

                 k = nint(acos(cos_scatt) * real(nang_scatt)/pi)
                 if (k > nang_scatt) k = nang_scatt
                 if (k < 0) k = 0

                 ! Calcul de l'angle omega
                 call rotation(u,v,w,-u_ray_tracing,-v_ray_tracing,-w_ray_tracing,v1pi,v1pj,v1pk)
                 xnyp = sqrt(v1pk*v1pk + v1pj*v1pj)
                 if (xnyp < 1e-10) then
                    xnyp = 0.0
                    costhet = 1.0_dp ! omega = 0
                 else
                    costhet = v1pj / xnyp
                 endif

                 ! calcul de l'angle entre la normale et l'axe z (theta)
                 theta = acos(costhet)
                 if (theta >= pi) theta = 0.0

                 !     le plan de diffusion est a +ou- 90deg de la normale
                 ! Ne doit pas etre utilise en ray-tracing !!! teste OK sans
                 !theta = theta  + pi_sur_deux

                 !----dans les matrices de rotation l'angle est omega = 2 * theta-----
                 omega = 2.0_dp * theta ! A verifier: le moins est pour corriger un bug de signe trouve par Marshall (supprime)
                 !     prochain if car l'arccos va de 0 a pi seulement
                 !     le +/- pour faire la difference dans le sens de rotation
                 if (v1pk < 0.0) omega = -1.0_dp * omega

                 ! Matrices de rotations
                 cosw = cos(omega)
                 sinw = sin(omega)
                 if (abs(cosw) < 1e-06) cosw = 0.0_dp
                 if (abs(sinw) < 1e-06) sinw = 0.0_dp

                 tab_cosw(theta_I,phi_I,iscatt,dir) = cosw
                 tab_sinw(theta_I,phi_I,iscatt,dir) = sinw
              endif ! lsepar_pola
           enddo !theta_I
        enddo !phi_I

     enddo ! iscatt
  enddo ! dir


!$omp parallel &
  !$omp default(none) &
  !$omp shared(lvariable_dust,Inu,I_sca2,n_cells,tab_s11_pos,uv0,w0,n_Stokes, tab_mueller_pos, lmueller) &
  !$omp shared(tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos,icell_ref,energie_photon,volume) &
  !$omp shared(lsepar_pola,tab_k,tab_sin_scatt_norm,lambda,p_lambda,n_phi_I,n_theta_I,nang_ray_tracing,lsepar_contrib) &
  !$omp shared(s11_save,tab_cosw,tab_sinw,nb_proc,kappa,kappa_factor,tab_albedo_pos) &
  !$omp private(iscatt,id,u_ray_tracing,v_ray_tracing,w_ray_tracing,theta_I,phi_I,i1,i2,u,v,w,cos_scatt,sin_scatt) &
  !$omp private(sum_sin,icell,p_icell,stokes,s11,k,alloc_status,dir,phi_scatt) &
  !$omp private(s12,s33,s34,M,ROP,RPO,v1pi,v1pj,v1pk,xnyp,costhet,theta,omega,cosw,sinw,C,D,S,facteur,kappa_sca)
  id = 1 ! pour code sequentiel

  ! Matrice de Mueller
  M = 0.0_dp

  ! matrices de rotation
  ROP = 0.0_dp
  RPO = 0.0_dp

  RPO(1,1) = 1.0_dp
  ROP(1,1) = 1.0_dp
  RPO(4,4) = 1.0_dp
  ROP(4,4) = 1.0_dp

  if (.not.lvariable_dust) then   ! we precalculate the s11 as they are all the same
     icell = icell_ref

     do dir=0,1
        do iscatt = 1, nang_ray_tracing
           ! Loop over the specific intensity bins
           do phi_I=1,n_phi_I
              do theta_I=1,n_theta_I
                 ! Average of the phase function over the bin
                 s11 = 0. ;
                 do i2=1, N_super
                    do i1 = 1, N_super
                       k = tab_k(i1,i2,theta_I,phi_I,iscatt,dir)
                       if (lmueller) then
                          s11 = s11 + tab_mueller_pos(1,1,k,icell,p_lambda) * tab_sin_scatt_norm(i1,i2,theta_I,phi_I,iscatt,dir)
                       else
                          s11 = s11 + tab_s11_pos(k,icell,p_lambda) * tab_sin_scatt_norm(i1,i2,theta_I,phi_I,iscatt,dir)
                       endif !lmueller
                    enddo ! i1
                 enddo !i2
                 s11_save(theta_I,phi_I,iscatt,dir) = s11

              enddo ! theta_I
           enddo ! phi_I

        enddo ! iscatt
     enddo ! dir
  endif ! .not.lvariable_dust

  p_icell = icell_ref
  !$omp do schedule(static, n_cells/nb_proc)
  do icell = 1, n_cells
     !$ id = omp_get_thread_num() + 1
     if (lvariable_dust) p_icell = icell

     ! Boucle sur les directions de ray-tracing
     do dir=0,1
        do iscatt = 1, nang_ray_tracing

           do phi_I=1,n_phi_I
              do theta_I=1,n_theta_I

                 if (lvariable_dust) then
                    ! Average of the phase over the bin
                    s11 = 0.
                    do i2=1, N_super
                       do i1 = 1, N_super
                          k = tab_k(i1,i2,theta_I,phi_I,iscatt,dir)
                          if (lmueller) then
                             s11 = s11 + tab_mueller_pos(1,1,k,p_icell,p_lambda) * &
                                  tab_sin_scatt_norm(i1,i2,theta_I,phi_I,iscatt,dir)
                          else
                             s11 = s11 + tab_s11_pos(k,p_icell,p_lambda) * tab_sin_scatt_norm(i1,i2,theta_I,phi_I,iscatt,dir)
                          endif !lmueller
                       enddo ! i1
                    enddo !i2
                 else ! does not depend on icell, we use the stored value
                    s11 = s11_save(theta_I,phi_I,iscatt,dir)
                 endif


                 if (lsepar_pola) then ! On calcule les s12, s33, s34 et la matrice de rotation
                    ! The polarisability phase function is not as steep as the phase function,
                    ! so we only use the center of the bin for the polarization to make it faster
                    i1 = N_super/2+1 ; i2 = i1
                    k = tab_k(i1,i2,theta_I,phi_I,iscatt,dir)

                    cosw = tab_cosw(theta_I,phi_I,iscatt,dir)
                    sinw = tab_sinw(theta_I,phi_I,iscatt,dir)

                    RPO(2,2) = cosw
                    ROP(2,2) = cosw
                    RPO(2,3) = sinw
                    ROP(2,3) = -1.0_dp * sinw
                    RPO(3,2) = -1.0_dp * sinw
                    ROP(3,2) = sinw
                    RPO(3,3) = cosw
                    ROP(3,3) = cosw

                   if (lmueller) then

                       ! Champ de radiation
                       stokes(:) = Inu(1:4,theta_I,phi_I,icell)

                       M(:,:) = tab_mueller_pos(:,:,k,p_icell,p_lambda)*(-s11)
                       M(:,4) = -M(:,4)
                       M(4,:) = -M(4,:)
                       M(2,2) = -M(2,2)
                       M(1,1) = s11
                       !On passe à l'opposé la quatrième collonne et ligne à cause de la convention
                       !Le - M(2,2) est pour être cohérent avec les opérations effectuées dans le cas mie

                    else

                       s12 = - s11 * tab_s12_o_s11_pos(k,p_icell,p_lambda)
                       s33 = - s11 * tab_s33_o_s11_pos(k,p_icell,p_lambda)
                       s34 = - s11 * tab_s34_o_s11_pos(k,p_icell,p_lambda)


                       ! Champ de radiation
                       stokes(:) = Inu(1:4,theta_I,phi_I,icell)

                       M(1,1) = s11 ; M(2,2) = s11 ; M(1,2) = s12 ; M(2,1) = s12
                       M(3,3) = s33 ; M(4,4) = s33 ; M(3,4) = -s34 ; M(4,3) = s34

                    endif

                    !  STOKE FINAL = RPO * M * ROP * STOKE INITIAL

                    if (lmueller) then
                       ! 1ere rotation
                       C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
                       C(1)=stokes(1)
                       C(4)=stokes(4)

                       ! multiplication matrice Mueller
                       D=matmul(M,C)

                       ! 2nde rotation
                       S(2:3)=matmul(RPO(2:3,2:3),D(2:3))
                       S(1)=D(1)
                       S(4)=D(4)
                       S(3)=-S(3)

                    else

                       ! 1ere rotation
                       C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
                       C(1)=stokes(1)
                       C(4)=stokes(4)

                       ! multiplication matrice Mueller par bloc
                       D(1:2)=matmul(M(1:2,1:2),C(1:2))
                       D(3:4)=matmul(M(3:4,3:4),C(3:4))

                       ! 2nde rotation
                       S(2:3)=matmul(RPO(2:3,2:3),D(2:3))
                       S(1)=D(1)
                       S(4)=D(4)
                       S(3)=-S(3)
                    endif

                    I_sca2(1:4,iscatt,dir,icell) = I_sca2(1:4,iscatt,dir,icell) + S(:)

                 else ! lsepar_pola
                    ! Champ de radiation
                    stokes(1) = Inu(1,theta_I,phi_I,icell)
                    I_sca2(1,iscatt,dir,icell) = I_sca2(1,iscatt,dir,icell) + s11 * stokes(1)
                 endif  ! lsepar_pola

                 if (lsepar_contrib) then
                    I_sca2(n_Stokes+2,iscatt,dir,icell) = I_sca2(n_Stokes+2,iscatt,dir,icell) + &
                         s11 * Inu(n_Stokes+2,theta_I,phi_I,icell)

                    I_sca2(n_Stokes+4,iscatt,dir,icell) = I_sca2(n_Stokes+4,iscatt,dir,icell) + &
                         s11 * Inu(n_Stokes+4,theta_I,phi_I,icell)
                 endif ! lsepar_contrib

              enddo ! theta_I
           enddo ! phi_I

        enddo ! iscatt
     enddo ! dir

     ! Normalisation
     facteur = energie_photon / volume(icell)
     kappa_sca = kappa(p_icell,lambda) * kappa_factor(icell) * tab_albedo_pos(p_icell,lambda)
     I_sca2(:,:,:,icell) =  I_sca2(:,:,:,icell) *  facteur * kappa_sca
  enddo ! icell
  !$omp enddo
  !$omp end parallel

  deallocate(Inu)
  return

end subroutine calc_Isca_rt2

!***********************************************************

subroutine calc_Isca_rt2_star(lambda,p_lambda,ibin)
  ! Version acceleree !TMP
  ! Routine liee a la precedente
  ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
  ! 20/04/2023

  integer, intent(in) :: lambda, p_lambda, ibin

  integer :: p_ri, p_zj, dir, iscatt, id
  real(kind=dp), dimension(4) :: Stokes, S, C, D
  real(kind=dp) :: x, y, z, u, v, w, facteur, kappa_sca

  real(kind=dp), dimension(4,4) ::  M, ROP, RPO

  integer :: k, icell, p_icell
  real :: s11, s12, s33, s34, cos_scatt
  real(kind=dp) :: omega, sinw, cosw, norme, energie_photon, n_photons_envoyes

  real(kind=dp), parameter :: prec = 0._dp

  if (lmono0) then ! image
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (real(nbre_photons_loop)*real(nbre_photons_image) *  AU_to_cm * pi)
  else ! SED
     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (n_photons_envoyes *  AU_to_cm * pi)
  endif

  if (n_etoiles > 1) then
     write(*,*) "WARNING : RT2 needs to be updated for multiple stars"
     write(*,*) " if extra stars/planets emit significantly"
  endif

  ! Matrice de Mueller
  M = 0.0_dp

  ! matrices de rotation
  ROP = 0.0_dp
  RPO = 0.0_dp

  RPO(1,1) = 1.0_dp
  ROP(1,1) = 1.0_dp
  RPO(4,4) = 1.0_dp
  ROP(4,4) = 1.0_dp

  stokes(:) = 0.0_dp
  p_ri = 1
  p_zj = 1

  id = 1 ! sequentiel

  ! Boucle sur les cellules
  !$omp parallel &
  !$omp default(none) &
  !$omp shared(n_cells,lambda,p_lambda,ibin,I_spec_star,nang_ray_tracing_star,cos_thet_ray_tracing_star) &
  !$omp shared(lvariable_dust,r_grid,z_grid,icell_ref,omega_ray_tracing_star,lsepar_pola) &
  !$omp shared(tab_s11_pos,tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos,eps_dust2_star) &
  !$omp shared(energie_photon,volume,tab_albedo_pos,kappa,kappa_factor, tab_mueller_pos, lmueller) &
  !$omp private(id,icell,p_icell,stokes,x,y,z,norme,u,v,w,dir,iscatt,cos_scatt,k,omega,cosw,sinw,RPO,ROP,S) &
  !$omp private(s11,s12,s33,s34,C,D,M,facteur,kappa_sca)
  p_icell = icell_ref
  !$omp do
  do icell=1, n_cells
     !$ id = omp_get_thread_num() + 1
     if (lvariable_dust) p_icell = icell

     ! Champ de radiation
     stokes(1) = sum(I_spec_star(icell,:))

     if (stokes(1) < 1.e-30_dp) cycle

     ! Direction de vol : lumiere stellaire directe : direction == position
     ! WARNING : only works for single star at x,y,z = 0 !!!
     x = 0.0_dp
     y = r_grid(icell) ! doit juste etre non nul
     z = z_grid(icell) ! sert a rien

     norme = sqrt(y**2 + z**2)
     u = 0.0_dp
     v = y / norme
     w = z / norme

     ! cos_scatt_ray_tracing correspondants
     call angles_scatt_rt2(id,ibin,x,y,z,u,v,w,.true.)

     ! Boucle sur les directions de ray-tracing
     do dir=0,1
        do iscatt = 1, nang_ray_tracing_star

           ! Angle de diffusion
           cos_scatt = cos_thet_ray_tracing_star(iscatt,dir,id)

           ! Dichotomie
!!$      kmin=0 ; kmax=nang_scatt ; k=(kmin+kmax)/2
!!$      do while ((kmax-kmin) > 1)
!!$         if (tab_cos_scatt(k) > cos_scatt) then ! cos decroissant
!!$            kmin = k
!!$         else
!!$            kmax = k
!!$         endif
!!$         k = (kmin + kmax)/2
!!$      enddo   ! while
!!$      k=kmax
!!$
!!$      ! Interpolation lineaire (en cos) des elements de la matrice de Mueller
!!$      frac =  (tab_cos_scatt(k-1) - cos_scatt) / (tab_cos_scatt(k-1) - tab_cos_scatt(k))
!!$      frac_m1 = 1.0_dp - frac

           ! 2 methode ---> semble OK
           ! Plus propre que l'interpolation avec la dicotomie ci-dessus
           ! en particulier quand nang_ray_tracing_star est faible
           ! la dichotomie + interpolation donne n'inporte quoi
           ! ---> c'est pas une bonne idee d'interpoler en cos car c'est rend la fct de phase encore plus pique pres de 0

           k = nint(acos(cos_scatt) * real(nang_scatt)/pi)
           if (k > nang_scatt) k = nang_scatt
           if (k < 1) k = 1 ! pour resoudre bug --> ne change pour le bench
           !frac = 1. ; frac_m1 = 0.

           ! Interpolation en lineaire ne marche pas non plus
           ! sans doute parce que le fct de phase de phase est tres pique
           ! et que c'est de toute facon pas une bonne idee d'interpoler
           ! --> le pb ne vient pas de la dichotomie mais de l'interpolation

           !   angle = acos(cos_scatt) * real(nang_scatt)/pi
           !   k = ceiling(angle)
           !   if (k > nang_scatt) then
           !      k = nang_scatt
           !      frac = 1. ; frac_m1 = 0.
           !   else
           !      frac =  (k - angle)
           !      frac_m1 = 1.0 - frac
           !   endif

           if (lsepar_pola) then
              ! Rotation pour se recaller sur le Nord celeste
              omega = omega_ray_tracing_star(iscatt,dir,id)

              cosw = cos(omega)
              sinw = sin(omega)
              if (abs(cosw) < 1e-06) cosw = 0.0_dp
              if (abs(sinw) < 1e-06) sinw = 0.0_dp

              RPO(2,2) = COSW
              ROP(2,2) = COSW
              RPO(2,3) = SINW
              ROP(2,3) = -1.0_dp * SINW
              RPO(3,2) = SINW
              ROP(3,2) = SINW
              RPO(3,3) = -COSW
              ROP(3,3) = COSW

              ! Matrice de Mueller

              if (lmueller) then
                 s11 = tab_mueller_pos(1,1,k,p_icell,p_lambda)
                 M(:,:) = tab_mueller_pos(:,:,k,p_icell,p_lambda)*(-s11)
                 M(:,4) = -M(:,4)
                 M(4,:) = -M(4,:)
                 M(2,2) = -M(2,2)
                 M(1,1) = s11
                 !On passe à l'opposé la quatrième collonne et ligne à cause de la convention
                 !Le - M(2,2) est pour être cohérent avec les opérations effectuées dans le cas mie

              else
                 s11 = tab_s11_pos(k,p_icell,p_lambda)
                 s12 = - s11 * tab_s12_o_s11_pos(k,p_icell,p_lambda)
                 s33 = - s11 * tab_s33_o_s11_pos(k,p_icell,p_lambda)
                 s34 = - s11 * tab_s34_o_s11_pos(k,p_icell,p_lambda)

                 M(1,1) = s11
                 M(2,2) = s11
                 M(1,2) = s12
                 M(2,1) = s12

                 M(3,3) = s33
                 M(4,4) = s33
                 M(3,4) = -s34
                 M(4,3) = s34
                 M(1,1) = s11
              endif !lmueller


              !  STOKE FINAL = RPO * M * ROP * STOKE INITIAL

              if (lmueller) then
                 ! 1ere rotation
                 C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
                 C(1)=stokes(1)
                 C(4)=stokes(4)

                 ! multiplication matrice Mueller par bloc
                 D=matmul(M,C)

                 ! 2nde rotation
                 S(2:3)=matmul(RPO(2:3,2:3),D(2:3))
                 S(1)=D(1)
                 S(4)=D(4)

              else

                 ! 1ere rotation
                 C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
                 C(1)=stokes(1)
                 C(4)=stokes(4)

                 ! multiplication matrice Mueller par bloc
                 D(1:2)=matmul(M(1:2,1:2),C(1:2))
                 D(3:4)=matmul(M(3:4,3:4),C(3:4))

                 ! 2nde rotation
                 S(2:3)=matmul(RPO(2:3,2:3),D(2:3))
                 S(1)=D(1)
                 S(4)=D(4)

              endif

              eps_dust2_star(:,iscatt,dir,icell) =  eps_dust2_star(:,iscatt,dir,icell) + S(:)
           else ! .not.lsepar_pola
              if (lmueller) then
                 s11 = tab_mueller_pos(1,1,k,p_icell,p_lambda)
              else
                 s11 = tab_s11_pos(k,p_icell,p_lambda)
              endif
              eps_dust2_star(1,iscatt,dir,icell) =  eps_dust2_star(1,iscatt,dir,icell) + s11 * stokes(1)
           endif ! lsepar_pola

        enddo ! iscatt
     enddo ! dir

     ! Normalisation
     facteur = energie_photon / volume(icell)
     kappa_sca = kappa(p_icell,lambda) * kappa_factor(icell) * tab_albedo_pos(p_icell,lambda)
     eps_dust2_star(:,:,:,icell) =  eps_dust2_star(:,:,:,icell) *  facteur * kappa_sca

  enddo ! icell
  !$omp end do
  !$omp end parallel

  return

end subroutine calc_Isca_rt2_star

!***********************************************************

function dust_source_fct(icell, x,y,z)
  ! La direction du rayon est maintenant fixee le long de l'axe x
  ! l'angle de diffusion ne depend que de la position x, y, z

  integer, intent(in) :: icell
  real(kind=dp), intent(in) :: x, y, z

  real(kind=dp), dimension(N_type_flux) :: dust_source_fct, SF1, SF2, SF3, SF4

  real(kind=dp) :: phi_pos, frac, un_m_frac, xiscatt, frac_r, frac_z, r
  integer :: iscatt1, iscatt2, dir, psup, ri1, zj1, ri2, zj2
  integer :: k, icell_tmp, ri, zj, phik

  SF1 = 0 ; SF2 = 0 ; SF3 = 0 ; SF4 = 0

  ! Interpolations angulaires
  if (lscatt_ray_tracing1) then
     if (l3D) then
        psup = 1
        k = 1
     else
        if (z > 0.0_dp) then
           psup=1
        else
           psup=2
        endif
        ! Cette methode est OK mais on voit les cellules
        phi_pos = atan2(x,y)
        k = floor(modulo(phi_pos, deux_pi) / deux_pi * n_az_rt) + 1
        if (k > n_az_rt) k = n_az_rt
     endif

     ! OK : lisse les images par rapport a la methode en dessous
!---     xi_az =  modulo(phi_pos, deux_pi) / deux_pi * n_az_rt + 0.5
!---     phi_k =  floor(xi_az)
!---     frac = xi_az - phi_k
!---     un_m_frac = 1.0_dp - frac
!---
!---     phi_k_p1 = phi_k + 1
!---     if (phi_k <=0) phi_k = 1
!---     if (phi_k_p1 > n_az_rt) phi_k_p1=n_az_rt
!---
!---     dust_source_fct(:) = eps_dust1(:,ri,zj,phi_k_p1,psup) * frac + eps_dust1(:,ri,zj,phi_k,psup) * un_m_frac
     dust_source_fct(:) = eps_dust1(k,psup,:,icell)  ! ??? ce n'est pas lineaire

  else ! Methode 2 --> only 2D, so phik should be 1
     ri = cell_map_i(icell)
     zj = cell_map_j(icell)
     phik = cell_map_k(icell)

     ! Pour interpolations spatiales
     r = sqrt(x*x + y*y)

     if (r > r_grid(icell)) then
        ri1 = ri
        ri2 = ri + 1
     else
        ri1 = ri - 1
        ri2 = ri
     endif

     if (abs(z) > z_grid(icell)) then
        zj1 = zj
        zj2 = zj + 1
     else
        zj1 = zj - 1
        zj2 = zj
     endif

     if (ri2 > n_rad) then
        ri2 = n_rad
        frac_r = 0._dp
     else  if (ri1 < 1) then
        ri1 = 1
        frac_r = 0._dp
     else
        frac_r = (log(r_grid(cell_map(ri2,zj,phik))) - log(r)) / &
             (log(r_grid(cell_map(ri2,zj,phik))) - log(r_grid(cell_map(ri1,zj,phik))))
     endif

     if (zj2 > nz) then
        zj2 = nz
        frac_z = 1.0_dp
     else  if (zj1 < 1) then
        zj1 = 1
        frac_z = 1.0_dp
     else
        frac_z = (z_grid(cell_map(ri,zj2,phik)) - abs(z)) / &
             (z_grid(cell_map(ri,zj2,phik)) - z_grid(cell_map(ri,zj1,phik)))
     endif

     ! Ok si je decommente les 2
     ri1 = ri    ! OK si je decommente seulement r ---> le pb est l'interpolation en r
     !zj1 = zj    ! Pb si je decommente seulement z


     ! Patch : dans les cas ou frac_z est > 1 ou < 0  !!!!!
     ! ca veut surement dire que l'on est pas dans la bonne cellule
     frac_z = max(min(1.0,frac_z),0.0)
     frac_r = max(min(1.0,frac_r),0.0)


     ! Calcul angle phi a la position
     phi_pos = modulo(atan2(x,y) + deux_pi,deux_pi)

     if (z > 0.0_dp) then
        dir=1
     else
        dir=0
     endif

     !----------------------------------------------------
     ! Emissivite du champ diffuse plusieurs fois
     !----------------------------------------------------
     xiscatt = max(phi_pos/deux_pi  * nang_ray_tracing, 0.0_dp)

     iscatt1 = floor(xiscatt)
     frac = xiscatt - iscatt1
     un_m_frac = 1.0_dp - frac
     iscatt2 = iscatt1 + 1

     ! Les limites periodiques
     !    if (iscatt2 > nang_ray_tracing) iscatt2 = 1
     !    if (iscatt1 < 1) iscatt1 = nang_ray_tracing
     iscatt1 = modulo(iscatt1,nang_ray_tracing)
     if (iscatt1==0) iscatt1 = nang_ray_tracing

     iscatt2 = modulo(iscatt2,nang_ray_tracing)
     if (iscatt2==0) iscatt2 = nang_ray_tracing

     ! Fct source des cellules
     icell_tmp = cell_map(ri1,zj1,phik)
     SF1(1) = eps_dust2(1,iscatt2,dir,icell_tmp) * frac + eps_dust2(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF1(2:3) = interpolate_Stokes_QU(eps_dust2(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif
     if (lsepar_contrib) then
        SF1(n_Stokes+1:N_type_flux) = eps_dust2(n_Stokes+1:N_type_flux,iscatt2,dir,icell_tmp) * frac + &
             eps_dust2(n_Stokes+1:N_type_flux,iscatt1,dir,icell_tmp) * un_m_frac
     endif

     icell_tmp = cell_map(ri2,zj1,phik)
     SF2(1) = eps_dust2(1,iscatt2,dir,icell_tmp) * frac + eps_dust2(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF2(2:3) = interpolate_Stokes_QU(eps_dust2(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif
     if (lsepar_contrib) then
        SF2(n_Stokes+1:N_type_flux) = eps_dust2(n_Stokes+1:N_type_flux,iscatt2,dir,icell_tmp) * frac + &
             eps_dust2(n_Stokes+1:N_type_flux,iscatt1,dir,icell_tmp) * un_m_frac
     endif

     icell_tmp = cell_map(ri1,zj2,phik)
     SF3(1) = eps_dust2(1,iscatt2,dir,icell_tmp) * frac + eps_dust2(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF3(2:3) = interpolate_Stokes_QU(eps_dust2(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif
     if (lsepar_contrib) then
        SF3(n_Stokes+1:N_type_flux) = eps_dust2(n_Stokes+1:N_type_flux,iscatt2,dir,icell_tmp) * frac + &
             eps_dust2(n_Stokes+1:N_type_flux,iscatt1,dir,icell_tmp) * un_m_frac
     endif

     icell_tmp = cell_map(ri2,zj2,phik)
     SF4(1) = eps_dust2(1,iscatt2,dir,icell_tmp) * frac + eps_dust2(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF4(2:3) = interpolate_Stokes_QU(eps_dust2(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif
     if (lsepar_contrib) then
        SF4(n_Stokes+1:N_type_flux) = eps_dust2(n_Stokes+1:N_type_flux,iscatt2,dir,icell_tmp) * frac + &
             eps_dust2(n_Stokes+1:N_type_flux,iscatt1,dir,icell_tmp) * un_m_frac
     endif


     !----------------------------------------------------
     ! Emissivite du champ stellaire diffuse 1 fois
     !----------------------------------------------------
     xiscatt =  max(phi_pos/deux_pi  * nang_ray_tracing_star, 0.0_dp)

     iscatt1 = floor(xiscatt)
     frac = xiscatt - iscatt1
     un_m_frac = 1.0_dp - frac
     iscatt2 = iscatt1 + 1

     ! Les limites periodiques
     iscatt1 = modulo(iscatt1,nang_ray_tracing_star)
     if (iscatt1==0) iscatt1 = nang_ray_tracing_star

     iscatt2 = modulo(iscatt2,nang_ray_tracing_star)
     if (iscatt2==0) iscatt2 = nang_ray_tracing_star


     ! Fct source des cellules
!     SF1(:) = SF1(:) + exp(log(eps_dust2_star(:,iscatt2,dir,ri1,zj1)) * frac + log(eps_dust2_star(:,iscatt1,dir,ri1,zj1)) * un_m_frac)
     ! TODO : passer en log ameliore un peu la forme a longue distance mais cree des bug si 0. !!

     icell_tmp = cell_map(ri1,zj1,phik)
     SF1(1) = SF1(1) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
          + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac

     if (lsepar_pola) then
        SF1(2:3) = SF1(2:3) + interpolate_Stokes_QU(eps_dust2_star(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2_star(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif

     icell_tmp = cell_map(ri2,zj1,phik)
     SF2(1) = SF2(1) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
          + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF2(2:3) = SF2(2:3) + interpolate_Stokes_QU(eps_dust2_star(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2_star(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif

     icell_tmp = cell_map(ri1,zj2,phik)
     SF3(1) = SF3(1) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
          + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF3(2:3) = SF3(2:3) + interpolate_Stokes_QU(eps_dust2_star(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2_star(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif

     icell_tmp = cell_map(ri2,zj2,phik)
     SF4(1) = SF4(1) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
          + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
     if (lsepar_pola) then
        SF4(2:3) = SF4(2:3) + interpolate_Stokes_QU(eps_dust2_star(2:3,iscatt1,dir,icell_tmp), &
             eps_dust2_star(2:3,iscatt2,dir,icell_tmp),un_m_frac)
     endif

     if (lsepar_contrib) then
        icell_tmp = cell_map(ri1,zj1,phik)
        SF1(n_Stokes+2) = SF1(n_Stokes+2) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
             + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
        icell_tmp = cell_map(ri2,zj1,phik)
        SF2(n_Stokes+2) = SF2(n_Stokes+2) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
             + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
        icell_tmp = cell_map(ri1,zj2,phik)
        SF3(n_Stokes+2) = SF3(n_Stokes+2) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
             + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
        icell_tmp = cell_map(ri2,zj2,phik)
        SF4(n_Stokes+2) = SF4(n_Stokes+2) + eps_dust2_star(1,iscatt2,dir,icell_tmp) * frac &
             + eps_dust2_star(1,iscatt1,dir,icell_tmp) * un_m_frac
     endif

     frac_r = 1.0
     !frac_z = 1.0
     ! Fct source interpolee
     dust_source_fct(:) = &
          frac_r * frac_z * SF1(:) + &
          (1.0_dp - frac_r) * frac_z * SF2(:) + &
          frac_r * (1.0_dp - frac_z) * SF3(:) + &
          (1.0_dp - frac_r) * (1.0_dp - frac_z) * SF4(:)

  endif !lscatt_ray_tracing

  return

end function dust_source_fct

!***********************************************************

function interpolate_Stokes_QU(PI_deuxtheta1,PI_deuxtheta2,frac1)
  ! This fonction interpolates Stokes Q and U
  ! by a linear interpolation in polarisation intensity and angle
  ! as interpolating Q and U directly creates artefacts
  ! (reduced pola when interpolating 2 Stokes vectors that are not aligned)
  ! C. Pinte 1/11/2017

  real, dimension(2), intent(in) :: PI_deuxtheta1,PI_deuxtheta2
  real(dp), intent(in) :: frac1
  real, dimension(2) :: interpolate_Stokes_QU

  real :: PxI1, PxI2, PxI, deux_theta1, deux_theta2, deux_theta

  !PxI1 = sqrt(QU1(1)**2 + QU1(2)**2) ; deux_theta1 = atan2(QU1(2),QU1(1))
  !PxI2 = sqrt(QU2(1)**2 + QU2(2)**2) ; deux_theta2 = atan2(QU2(2),QU2(1))

  ! Using PxI and deux_theta
  PxI1 = PI_deuxtheta1(1) ; deux_theta1 = PI_deuxtheta1(2)
  PxI2 = PI_deuxtheta2(1) ; deux_theta2 = PI_deuxtheta2(2)

  ! Linear interpolation of polarized intensity
  PxI = PxI2 * (1.0_dp-frac1) + PxI1 * frac1

  ! Making sure there is no wrapping before the linear interpolation in angle
  if (abs(deux_theta2 - deux_theta1) >= pi) then
     if (deux_theta2 > deux_theta1) then
        deux_theta1 = deux_theta1 + deux_pi
     else
        deux_theta2 = deux_theta2 + deux_pi
     endif
  endif
  deux_theta = deux_theta2 * (1.0_dp-frac1) + deux_theta1 * frac1

  interpolate_Stokes_QU = PxI * [cos(deux_theta),-sin(deux_theta)] ! - sign to match IAU convention

  return

end function interpolate_Stokes_QU

!***********************************************************

end module dust_ray_tracing
