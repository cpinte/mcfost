! TODO : benchmark 2009 
! large g
! optimization memoire (ie, dimension 1 de eps_dust1)

module dust_ray_tracing

  use parametres
  use constantes
  use em_th
  use disk
  use opacity
  use prop_star
  use resultats
  use ray_tracing
  !$ use omp_lib

  use scattering

  implicit none

  contains


subroutine alloc_ray_tracing()
  ! Alloue les tableaux pour le ray-tracing
  ! C. Pinte
  ! 13/10/08
  
  integer :: alloc_status

  ! 15 15 90 90 OK pour benchmark

  n_phi_I = 15 ; ! 17/11/10 semble mieux si egal a nang_ray_tracing
  n_theta_I = 15 ;
  nang_ray_tracing = 15 ;
  nang_ray_tracing_star = 1000 ; ! Bug dans les cartes de pola si ce n'est pas le meme nbre

  if (lisotropic) then
     n_phi_I = 1 ;
     n_theta_I = 1 ;
     nang_ray_tracing = 1 ;
  endif

  allocate(kappa_sca(n_lambda,n_rad,nz+1,1), tab_s11_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt),stat=alloc_status)
  if (alloc_status > 0) then 
     write(*,*) 'Allocation error kappa_sca, tab_s11_ray_tracing'
     stop
  endif
  kappa_sca = 0. 
  tab_s11_ray_tracing = 0.

  if (lsepar_pola) then
     allocate(tab_s12_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt),  &
          tab_s33_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt),&
          tab_s34_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt),&
          tab_s12_o_s11_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt),  &
          tab_s33_o_s11_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt), &
          tab_s34_o_s11_ray_tracing(n_lambda,p_n_rad,p_nz,0:nang_scatt), &
          stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error tab_s12_ray_tracing'
        stop
     endif
     tab_s12_ray_tracing = 0.
     tab_s33_ray_tracing = 0.
     tab_s34_ray_tracing = 0.
     tab_s12_o_s11_ray_tracing = 0.
     tab_s33_o_s11_ray_tracing = 0.
     tab_s34_o_s11_ray_tracing = 0.
  endif

  if (lsed.and.(RT_sed_method == 1)) then
      allocate(Stokes_ray_tracing(n_lambda,1,1,RT_n_ibin,N_type_flux,nb_proc), stat=alloc_status)
  else
     allocate(Stokes_ray_tracing(n_lambda,igridx,igridy,RT_n_ibin,N_type_flux,nb_proc), stat=alloc_status)
  endif
  if (alloc_status > 0) then  
     write(*,*) 'Allocation error Stokes_ray_tracing'
     stop
  endif
  Stokes_ray_tracing = 0.0_db 

  allocate(J_th(n_rad,nz), stat=alloc_status)
  if (alloc_status > 0) then 
     write(*,*) 'Allocation error J_th'
     stop
  endif
  J_th = 0.0_db
  
  if (lscatt_ray_tracing1) then
     allocate(xI_scatt(N_type_flux,RT_n_ibin,n_rad,nz,n_az_rt,0:1,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error xI_scatt'
        stop
     endif
     xI_scatt = 0.0_db

     allocate(xsin_scatt(RT_n_ibin,n_rad,nz,n_az_rt,0:1,nb_proc), &
          xN_scatt(RT_n_ibin,n_rad,nz,n_az_rt,0:1,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error xsin_scatt'
        stop
     endif
     xsin_scatt = 0.0_db
     xN_scatt = 0.0_db

     !I_scatt = sum(xI_scatt(:,ibin,i,j,:,:,:),dim=4) * facteur * kappa_sca(lambda,i,j,1)
     allocate(I_scatt(N_type_flux,n_az_rt,0:1), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error I_scatt'
        stop
     endif
     I_scatt = 0.0_db
 

     allocate(eps_dust1(N_type_flux,n_rad,nz,n_az_rt,0:1), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error eps_dust1'
        stop
     endif
     eps_dust1 =0._db 

     allocate(itheta_rt1(RT_n_ibin,nb_proc), sin_scatt_rt1(RT_n_ibin,nb_proc), sin_omega_rt1(RT_n_ibin,nb_proc), &
          cos_omega_rt1(RT_n_ibin,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error itheta_rt1'
        stop
     endif
     itheta_rt1 = 1
     sin_scatt_rt1 = 0.0
     sin_omega_rt1 = 0.5
     cos_omega_rt1 = 0.5

  else !rt 2
     allocate(xI_star(n_rad,nz,nb_proc), xw_star(n_rad,nz,nb_proc), xl_star(n_rad,nz,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error xI_star'
        stop
     endif
     xI_star = 0.0_db ;  xl_star = 0.0_db 
     
     allocate(eps_dust2(N_type_flux,nang_ray_tracing,0:1,n_rad,nz), &
          I_sca2(N_type_flux,nang_ray_tracing,0:1,n_rad,nz,nb_proc), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error eps_dust2'
        stop
     endif
     eps_dust2 =0.0_db ; I_sca2 = 0.0_db ;
     
     allocate(eps_dust2_star(4,nang_ray_tracing_star,0:1,n_rad,nz), stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error eps_dust2_star'
        stop
     endif
     eps_dust2_star =0.0_db 
          
         
     allocate(cos_thet_ray_tracing(nang_ray_tracing,0:1,nb_proc), omega_ray_tracing(nang_ray_tracing,0:1,nb_proc),&
          stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error cos_thet_ray_tracing'
        stop
     endif
     cos_thet_ray_tracing = 0._db
     omega_ray_tracing = 0._db

     allocate(cos_thet_ray_tracing_star(nang_ray_tracing_star,0:1,nb_proc), &
          omega_ray_tracing_star(nang_ray_tracing_star,0:1,nb_proc),&
          stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error cos_thet_ray_tracing_star'
        stop
     endif
     cos_thet_ray_tracing_star = 0._db
     omega_ray_tracing_star = 0._db
     
     allocate(xI(N_type_flux,n_theta_I,n_phi_I,n_rad,nz,nb_proc),stat=alloc_status)
     if (alloc_status > 0) then 
        write(*,*) 'Allocation error xI'
        stop
     endif
     xI = 0.0_db
  endif

  return

end subroutine alloc_ray_tracing

!***********************************************************

subroutine dealloc_ray_tracing()
  ! desalloue les tableaux pour le ray-tracing
  ! C. Pinte
  ! 13/10/08
  
     
  deallocate(eps_dust2, eps_dust2_star, I_sca2)
  deallocate(cos_thet_ray_tracing,omega_ray_tracing)
  deallocate(xI)

  return

end subroutine dealloc_ray_tracing

!***********************************************************

subroutine init_directions_ray_tracing()
  ! initialise les directions (inclinaisons) pour le ray-tracing
  ! C. Pinte 
  ! 09/09/08
 
  real(kind=db) :: cos_min, cos_max
  integer :: ibin
  
  allocate(tab_RT_incl(RT_n_ibin),tab_u_rt(RT_n_ibin),tab_v_rt(RT_n_ibin),tab_w_rt(RT_n_ibin))

  if (RT_n_ibin==1) then
     tab_RT_incl(1) = RT_imin
  else
     cos_min = cos(RT_imin / 180.0_db * pi)
     cos_max = cos(RT_imax / 180.0_db * pi)

     if (lRT_i_centered) then
        do ibin=1, RT_n_ibin
           tab_RT_incl(ibin) = acos(  cos_min + (real(ibin) -0.5)/real(RT_n_ibin) * (cos_max - cos_min) ) /pi * 180.0_db
        enddo
     else
        do ibin=1, RT_n_ibin
           tab_RT_incl(ibin) = acos(  cos_min + (real(ibin) - 1)/(real(RT_n_ibin) -1) * (cos_max - cos_min) ) /pi * 180.0_db
        enddo
     endif
  endif

  do ibin=1, RT_n_ibin
     ! 0 est remplace par un epsilon car il faut donner un axe de reference
     ! pour les differentes directions de ray-tracing utilisees dans le RT2
     tab_u_rt(ibin) = sin(max(tab_RT_incl(ibin),1e-20)/180._db*pi) ! uv_rt mais v_rt = 0 ici
     tab_w_rt(ibin) = sqrt(1.0_db - tab_u_rt(ibin)*tab_u_rt(ibin))
  enddo
  tab_v_rt(:) = 0.0

  return

end subroutine init_directions_ray_tracing

!***********************************************************

subroutine angles_scatt_ray_tracing(id,ibin,x,y,z,u,v,w,lstar)
  ! Calcul les cosinus des angles de diffusion vers toutes les directions
  ! utilisees pour le ray tracing
  ! Calcul aussi les matrices de rotation pour le vecteur de Stokes 
  ! Remplit les tableaux cos_thet_ray_tracing et omega_ray_tracing
  ! RT2
  ! C. Pinte
  ! 19/01/08

  ! TODO: cos_thet_ray_tracing n'est plus utilise depuis le super-echantillonage

  integer, intent(in) :: id, ibin
  real(kind=db), intent(in) :: x,y,z, u, v, w  ! z est inutilise
  logical, intent(in) :: lstar
  real(kind=db) :: uv0, w0, phi, phi_pos,u_ray_tracing,  v_ray_tracing, w_ray_tracing, prod1, prod2
  real(kind=db) :: w2, v1pi, v1pj, v1pk, xnyp, costhet, theta, omega
  integer :: iscatt, direction, N

  ! Direction observateur dans repere refence
  uv0 = tab_u_rt(ibin) ; w0 = tab_w_rt(ibin) ! epsilon needed here

  ! Calcul angle phi a la position
  phi_pos = modulo(atan2(y,x) + deux_pi,deux_pi)

  ! les directions de ray-tracing sont definis pr rapport a la ref u0
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
     
     phi = phi + phi_pos
        
     u_ray_tracing = uv0 * cos(phi)
     v_ray_tracing = uv0 * sin(phi)
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
              costhet = 1.0
           else
              costhet = -1.0*v1pj / xnyp
           endif
              
           ! calcul de l'angle entre la normale et l'axe z (theta)
           theta = acos(costhet)
           if (theta >= pi) theta = 0.0
     
           !     le plan de diffusion est a +ou- 90deg de la normale
           ! Ne doit pas etre utilise en ray-tracing !!! teste OK sans
           !theta = theta  + pi_sur_deux  
              
           !----dans les matrices de rotation l'angle est omega = 2 * theta-----
           omega = 2.0_db * theta ! A verifier: le moins est pour corriger un bug de signe trouve par Marshall (supprime)
           !     prochain if car l'arccos va de 0 a pi seulement
           !     le +/- pour faire la difference dans le sens de rotation
           if (v1pk < 0.0) omega = -1.0_db * omega
           
           if (lstar) then
              omega_ray_tracing_star(iscatt,direction,id) = omega
           else
              omega_ray_tracing(iscatt,direction,id) = omega
           endif
        endif
     enddo ! direction
  enddo !iscatt

  return

end subroutine angles_scatt_ray_tracing

!***********************************************************

subroutine angle_scatt_rt(id,u,v,w)
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
  real(kind=db), intent(in) :: u, v, w

  real(kind=db) :: v1pi, v1pj, v1pk, xnyp, costhet, theta, omega, cosw, sinw
  real :: cos_scatt
  integer :: k, ibin, dir
  
  ! Direction observateur dans repere refence
  do ibin=1,RT_n_ibin        
     cos_scatt = tab_u_rt(ibin) * u +  tab_w_rt(ibin) * w ! v_ray_tracing = 0

     k = nint(acos(cos_scatt) * real(nang_scatt)/pi) ! TODO : + 1 pour test
     if (k > nang_scatt) k = nang_scatt
     if (k < 1) k = 1 
     
     itheta_rt1(ibin,id) = k
     sin_scatt_rt1(ibin,id) = sqrt(1.0 - cos_scatt**2) 

     if (lsepar_pola) then
        ! Matrice de rotation pour se recaller sur le Nord celeste
        call rotation(u,v,w,-tab_u_rt(ibin),-tab_v_rt(ibin),-tab_w_rt(ibin),v1pi,v1pj,v1pk)
        xnyp = sqrt(v1pk*v1pk + v1pj*v1pj)
        if (xnyp < 1e-10) then
           xnyp = 0.0_db
           costhet = 1.0_db
        else
           costhet = -1.0_db*v1pj / xnyp
        endif
        
        ! calcul de l'angle entre la normale et l'axe z (theta)
        theta = acos(costhet)
        if (theta >= pi) theta = 0.0_db
     
        !     le plan de diffusion est a +ou- 90deg de la normale
        theta = theta + pi_sur_deux
  
        !----dans les matrices de rotation l'angle est omega = 2 * theta-----
        omega = 2.0_db * theta
        !     prochain if car l'arccos va de 0 a pi seulement
        !     le +/- pour faire la difference dans le sens de rotation
        if (v1pk < 0.0) omega = -1.0_db * omega
           
        cosw = cos(omega)
        sinw = sin(omega)
        if (abs(cosw) < 1e-06) cosw = 0.0_db
        if (abs(sinw) < 1e-06) sinw = 0.0_db
        
        cos_omega_rt1(ibin,id) = cosw
        sin_omega_rt1(ibin,id) = sinw
     endif ! lsepar_pola
  enddo ! ibin

  return

end subroutine angle_scatt_rt

!***********************************************************

subroutine calc_xI_scatt(id,lambda,ri,zj,phi_k,psup,l,stokes,flag_star)
  ! RT1
  ! Calcul les matrices de Mueller pour une direction de diffusion
  ! lance dans chaque cellule traversee
  ! utilise les resultats de angle_scatt_rt1
  ! C. Pinte
  ! 13/09/09, version intiale  19/01/08

  implicit none
  ! pour chaque cellule stocke le champ diffusee par un paquet dans les directions
  ! de ray-tracing ibin et dir en fonction de atan2(x,y) : ie un paquet 
  ! n'allume qu'un seul endroit de l'anneau au lieu de partout avant, 
  ! mais la taille des tableaux est la meme
  ! pour le ray-tracing, il ne faut plus interpoler entre les directions de ray-tracing
  ! mais entre les positions dans le disque

  real(kind=db), intent(in) :: Stokes, l
  logical, intent(in) :: flag_star
  integer, intent(in) :: id, lambda, ri, zj, phi_k, psup

  real(kind=db) :: flux
  integer :: ibin, it, p_ri, p_zj
  
  if (lstrat) then
     p_ri = ri
     p_zj = zj
  else
     p_ri = 1
     p_zj = 1
  endif

  do ibin = 1, RT_n_ibin
     it = itheta_rt1(ibin,id)
     flux = l * stokes * tab_s11_ray_tracing(lambda,p_ri,p_zj,it) !* sin_scatt_rt1(ibin,id) 
     ! TODO : est-ce qu'il ne faut pas moyenner par le sin de l'angle de scatt dans la cellule azimuthale ???
     
     xI_scatt(1,ibin,ri,zj,phi_k,psup,id) =  xI_scatt(1,ibin,ri,zj,phi_k,psup,id) + flux
     xsin_scatt(ibin,ri,zj,phi_k,psup,id) =  xsin_scatt(ibin,ri,zj,phi_k,psup,id) + 1.0_db !sin_scatt_rt1(ibin,id) 
     xN_scatt(ibin,ri,zj,phi_k,psup,id) =  xN_scatt(ibin,ri,zj,phi_k,psup,id) + 1.0_db

     if (lsepar_contrib) then
        if (flag_star) then
           xI_scatt(n_Stokes+2,ibin,ri,zj,phi_k,psup,id) =  xI_scatt(n_Stokes+2,ibin,ri,zj,phi_k,psup,id) + flux
        else
           xI_scatt(n_Stokes+4,ibin,ri,zj,phi_k,psup,id) =  xI_scatt(n_Stokes+4,ibin,ri,zj,phi_k,psup,id) + flux
        endif
     endif
  enddo ! ibin

  return

end subroutine calc_xI_scatt

!***********************************************************

subroutine calc_xI_scatt_pola(id,lambda,ri,zj,phi_k,psup,l,stokes)
  ! RT1
  ! Calcul les matrices de Mueller pour une direction de diffusion
  ! lance dans chaque cellule traversee
  ! utilise les resultats de angle_scatt_rt1
  ! C. Pinte
  ! 13/09/09, version intiale  19/01/08

  implicit none

  real(kind=db), dimension(4), intent(in) :: Stokes
  real(kind=db), intent(in) :: l
  integer, intent(in) :: id, lambda, ri, zj, phi_k, psup

  real(kind=db), dimension(4) :: C, D
  real(kind=db), dimension(4,4) ::  M, ROP, RPO
  real(kind=db) :: cosw, sinw
  real :: s11, s12, s33, s34
  integer :: ibin, it

  ROP = 0.0_db
  RPO = 0.0_db
  RPO(1,1) = 1.0_db
  ROP(1,1) = 1.0_db
  RPO(4,4) = 1.0_db
  ROP(4,4) = 1.0_db

  M = 0.0_db 

  do ibin = 1, RT_n_ibin
     ! Matrice de Mueller
     it = itheta_rt1(ibin,id)

     s11 = tab_s11_ray_tracing(lambda,ri,zj,it)
     s12 = tab_s12_ray_tracing(lambda,ri,zj,it)
     s33 = tab_s33_ray_tracing(lambda,ri,zj,it)
     s34 = tab_s34_ray_tracing(lambda,ri,zj,it)

     M(1,1) = s11 ; M(2,2) = s11 ; M(1,2) = s12 ; M(2,1) = s12
     M(3,3) = s33 ; M(4,4) = s33 ; M(3,4) = -s34 ; M(4,3) = s34

     ! Matrices de rotation 
     cosw = cos_omega_rt1(ibin,id)
     sinw = sin_omega_rt1(ibin,id)

     RPO(2,2) = cosw
     ROP(2,2) = cosw
     RPO(2,3) = sinw
     ROP(2,3) = -1.0_db * sinw
     RPO(3,2) = -1.0_db * sinw
     ROP(3,2) = sinw
     RPO(3,3) = cosw
     ROP(3,3) = cosw

     !  STOKE FINAL = RPO * M * ROP * STOKE INITIAL
     ! 1ere rotation
     C(2:3) = matmul(ROP(2:3,2:3),stokes(2:3))
     C(1)=stokes(1)
     C(4)=stokes(4)

     ! multiplication matrice Mueller par bloc
     D(1:2)=matmul(M(1:2,1:2),C(1:2))
     D(3:4)=matmul(M(3:4,3:4),C(3:4))
     
     ! 2nde rotation
     xI_scatt(2:3,ibin,ri,zj,phi_k,psup,id) =  xI_scatt(2:3,ibin,ri,zj,phi_k,psup,id) + l * matmul(RPO(2:3,2:3),D(2:3))
     xI_scatt(1,ibin,ri,zj,phi_k,psup,id) =  xI_scatt(1,ibin,ri,zj,phi_k,psup,id) + l * D(1)
     xI_scatt(4,ibin,ri,zj,phi_k,psup,id) =  xI_scatt(4,ibin,ri,zj,phi_k,psup,id) + l * D(4)
  enddo ! ibin

  ! TODO : 

  return

end subroutine calc_xI_scatt_pola


!***********************************************************

subroutine init_dust_source_fct1(lambda,ibin)
  ! RT1
  
  implicit none

  integer, intent(in) :: lambda, ibin

  integer :: i,j, alloc_status, l, k, itype
  real(kind=db) :: facteur, energie_photon, Temp, cst_wl, wl, cst, coeff_exp, cst_E, n_photons_envoyes
  real(kind=db), dimension(n_az_rt,0:1) :: norme

  if (lmono0) write(*,*) "i=", tab_RT_incl(ibin)

  eps_dust1(:,:,:,:,:) = 0.0_db 

  ! Contribution emission thermique directe
  call calc_Ith(lambda)
  !unpolarized_Stokes = 0. ;   unpolarized_Stokes(1) = 1. 

  ! TODO : la taille de eps_dust1 est le facteur limitant pour le temps de calcul
  
  if (lmono0) then ! image
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (real(nbre_photons_loop)*real(nbre_photons_image) *  AU_to_cm * pi)
  else ! SED
     !n_photons_envoyes = sum(n_phot_envoyes(:,lambda))
     !energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
     !(n_photons_envoyes *  AU_to_cm * (8*pi))

     ! OK pour test
     !n_photons_envoyes = sum(n_phot_envoyes(:,lambda))
     !energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda)) / n_photons_envoyes &
     !    * tab_lambda(lambda) * 1.0e-6 * 100!lambda.F_lambda

     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = hp * c_light**2 / 2. * (E_stars(lambda) + E_disk(lambda) ) / n_photons_envoyes &
         * tab_lambda(lambda) *1e-2/14.!lambda.F_lambda
  endif

  ! Intensite specifique diffusion
  do i=1,n_rad
     facteur = energie_photon / volume(i) * n_az_rt * 2 ! n_az_rt * 2 car subdivision virtuelle des cellules
     do j=1,nz
        if (kappa(lambda,i,j,1) > tiny_db) then
           ! TODO : pb de pola 
           do itype=1,N_type_flux
              norme = sum(xN_scatt(ibin,i,j,:,:,:),dim=3) / max(sum(xsin_scatt(ibin,i,j,:,:,:),dim=3),tiny_db) 
              I_scatt(itype,:,:) = sum(xI_scatt(itype,ibin,i,j,:,:,:),dim=3) * norme  * facteur * kappa_sca(lambda,i,j,1)
           enddo ! itype

           eps_dust1(1,i,j,:,:) =  (  I_scatt(1,:,:) +  J_th(i,j) ) / kappa(lambda,i,j,1) 
           if (lsepar_contrib) then
              eps_dust1(n_Stokes+2,i,j,:,:) =    I_scatt(n_Stokes+2,:,:) / kappa(lambda,i,j,1) 
              eps_dust1(n_Stokes+3,i,j,:,:) =    J_th(i,j) / kappa(lambda,i,j,1) 
              eps_dust1(n_Stokes+4,i,j,:,:) =    I_scatt(n_Stokes+4,:,:) / kappa(lambda,i,j,1) 
           endif ! lsepar_contrib
        else
           eps_dust1(:,i,j,:,:) = 0.0_db
        endif ! kappa > 0
     enddo ! j
  enddo !i
  
  
!  write(*,*) maxval(eps_dust1(1,:,:,:,:)), maxval(eps_dust1(2,:,:,:,:)), maxval(eps_dust1(3,:,:,:,:)), maxval(eps_dust1(4,:,:,:,:)), maxval(eps_dust1(5,:,:,:,:))

  ! Fouction source
!---  do i=1,n_rad
!---     do j=1,nz
!---         if (l_sym_ima) then
!---           do k = 1,n_az_rt/2
!---              eps_dust1(:,i,j,k,:) = 0.5 * (eps_dust1(:,i,j,k,:) + eps_dust1(i,j,n_az_rt-k+1,:,:))
!---              ! symetrique mais pas besoin
!---           enddo
!---        endif ! lsym_ima
!---     enddo ! j
!---  enddo !i

  return

end subroutine init_dust_source_fct1

!***********************************************************

subroutine init_dust_source_fct2(lambda,ibin)
  ! RT2
  ! calcule la fct source pour integrer par ray-tracing
  ! C. Pinte
  ! 25/09/08

  integer, intent(in) :: lambda, ibin
  integer :: i,j, iscatt, dir

  if (lmono0) write(*,*) "i=", tab_RT_incl(ibin)
  if (lmono0) write(*,'(a33, $)') " Scattered specific intensity ..."

  I_sca2 = 0.0_db 
  eps_dust2_star = 0.0_db
  eps_dust2 = 0.0_db

  ! Ajout du champ de radiation stellaire diffuse 1 seule fois
  call calc_Isca2_star(lambda, ibin) 

  ! Contribution lumiere diffusee (y compris multiple et thermique diffusee)
  call calc_Isca2_new(lambda, ibin) 

  ! Contribution emission thermique directe
  call calc_Ith(lambda)

  ! Fouction source, indices : pola, iscatt, dir, i, j
  do j=1,nz
     do i=1,n_rad
        if (kappa(lambda,i,j,1) > tiny_db) then
           ! Boucle sur les directions de ray-tracing
           do dir=0,1 
              do iscatt = 1, nang_ray_tracing 
                 eps_dust2(1,iscatt,dir,i,j) =  ( sum(I_sca2(1,iscatt,dir,i,j,:))  +  J_th(i,j) ) / kappa(lambda,i,j,1)
                                  
                 if (lsepar_pola) then
                    eps_dust2(2:4,iscatt,dir,i,j) =  sum(I_sca2(2:4,iscatt,dir,i,j,:),dim=2)  / kappa(lambda,i,j,1)
                 endif
              
                 if (lsepar_contrib) then
                    eps_dust2(n_Stokes+2,iscatt,dir,i,j) =    sum(I_sca2(n_Stokes+2,iscatt,dir,i,j,:)) / kappa(lambda,i,j,1) 
                    eps_dust2(n_Stokes+3,iscatt,dir,i,j) =    J_th(i,j) / kappa(lambda,i,j,1) 
                    eps_dust2(n_Stokes+4,iscatt,dir,i,j) =    sum(I_sca2(n_Stokes+4,iscatt,dir,i,j,:)) / kappa(lambda,i,j,1) 
                 endif ! lsepar_contrib
              enddo ! iscatt

              do iscatt = 1, nang_ray_tracing_star 
                 eps_dust2_star(:,iscatt,dir,i,j) = eps_dust2_star(:,iscatt,dir,i,j) / kappa(lambda,i,j,1)
              enddo ! iscatt
           enddo ! dir
        else
           eps_dust2(:,:,:,i,j) = 0.0_db
           eps_dust2_star(:,:,:,i,j) = 0.0_db
        endif
     enddo !i
  enddo !j

  if (lmono0) write(*,*)  "Done"
  
  return

end subroutine init_dust_source_fct2

!***********************************************************

subroutine calc_Ith(lambda)
  ! calcul emissivite thermique 
  ! C. Pinte
  ! 25/09/08

  integer, intent(in) :: lambda
  integer :: i, j, l, T
  real(kind=db) ::  Temp, cst_wl, wl, coeff_exp, cst_E


  ! longueur d'onde en metre
  wl = tab_lambda(lambda)*1.e-6

  J_th(:,:) = 0.0
 
  ! Intensite specifique emission thermique
  if ((l_em_disk_image).or.(lsed)) then
     if (lRE_LTE) then
        cst_E=2.0*hp*c_light**2 
        do i=1,n_rad
           do j=1,nz
              Temp=Temperature(i,j,1) ! que LTE pour le moment
              cst_wl=cst_th/(Temp*wl)
              if (cst_wl < 500.0) then
                 coeff_exp=exp(cst_wl)
                 J_th(i,j) = cst_E/((wl**5)*(coeff_exp-1.0)) * wl * kappa_abs_eg(lambda,i,j,1) ! Teste OK en mode SED avec echantillonnage lineaire du plan image
              else
                 J_th(i,j) = 0.0_db
              endif
           enddo
        enddo
     endif !lRE_LTE

     if (lRE_nLTE) then
        cst_E=2.0*hp*c_light**2 
        do i=1,n_rad
           do j=1,nz
              do l=grain_RE_nLTE_start,grain_RE_nLTE_end
                 Temp=Temperature_1grain(i,j,l) 
                 cst_wl=cst_th/(Temp*wl)
                 if (cst_wl < 500.0) then
                    coeff_exp=exp(cst_wl)
                    J_th(i,j) = J_th(i,j) + cst_E/((wl**5)*(coeff_exp-1.0)) * wl * &
                         q_abs(lambda,l)*densite_pouss(i,j,1,l) 
                 endif
              enddo
           enddo
        enddo
     endif !lRE_nLTE

     if (lnRE) then
        do i=1,n_rad
           do j=1,nz
              do l=grain_nRE_start,grain_nRE_end
                 if (l_RE(i,j,l)) then ! le grain a une temperature
                    Temp=Temperature_1grain_nRE(i,j,l)
                    cst_wl=cst_th/(Temp*wl)
                    if (cst_wl < 500.) then
                       coeff_exp=exp(cst_wl)
                       J_th(i,j) = J_th(i,j) + cst_E/((wl**5)*(coeff_exp-1.0)) * wl * &
                            q_abs(lambda,l)*densite_pouss(i,j,1,l) 
                    endif !cst_wl
                 else ! ! la grain a une proba de T
                    do T=1,n_T
                       temp=tab_Temp(T)
                       cst_wl=cst_th/(Temp*wl)
                       if (cst_wl < 500.) then
                          coeff_exp=exp(cst_wl)
                          J_th(i,j) = J_th(i,j) + cst_E/((wl**5)*(coeff_exp-1.0)) * wl * &
                               q_abs(lambda,l)*densite_pouss(i,j,1,l) * Proba_Temperature(T,i,j,l)
                       endif !cst_wl
                    enddo ! T
                 endif ! l_RE
              enddo !
           enddo !j
        enddo ! i
     endif !lnRE

  endif ! lsed or lemission_disk
  
  return

end subroutine calc_Ith

!***********************************************************

subroutine calc_Isca2_new(lambda,ibin)
  ! Version acceleree

  integer, intent(in) :: lambda, ibin

  integer :: ri, zj, p_ri, p_zj, theta_I, phi_I, dir, iscatt, id, dir_I
  real(kind=db), dimension(4) :: Stokes, S, C, D
  real(kind=db) :: x, y, z, u, v, w, phi, w02, facteur, energie_photon
  real(kind=db) :: u_ray_tracing, v_ray_tracing, w_ray_tracing, uv0, w0

  real(kind=db), dimension(:,:,:,:,:), allocatable :: Inu
  real(kind=db), dimension(n_rad,nz) :: Inu_max
  real(kind=db), dimension(4,4) ::  M, ROP, RPO

  integer :: k, kmin, kmax, alloc_status, i1, i2, correct_w
  real :: frac, frac_m1, cos_scatt, sum_sin, f1, f2, sin_scatt, phi_scatt
  real(kind=db) :: omega, sinw, cosw, n_photons_envoyes, v1pi, v1pj, v1pk, xnyp, costhet, theta


  integer, parameter :: N_super = 5 
  ! 5 cree un leger surcout dans le cas avec strat (qq 10 sec par inclinaison et par lambda)
  ! 15 cree un important surcout 

  real, dimension(N_super,N_super,n_theta_I,n_phi_I) :: tab_u, tab_v, tab_w

  real, dimension(:,:), allocatable :: s11, sum_s11, s12, s33, s34

  integer :: syst_status

  real :: rnd

  ! Direction observateur dans repere refence
  uv0 = tab_u_rt(ibin) ; w0 = tab_w_rt(ibin) ! epsilon needed here

  ! Allocation dynamique pour passer en stack
  allocate(Inu(N_type_flux,n_theta_I,n_phi_I,n_rad,nz), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error I_nu'
     stop
  endif
  Inu = 0.0_db

  id = 1
  
  if (lmono0) then ! image
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (real(nbre_photons_loop)*real(nbre_photons_image) *  AU_to_cm * pi)
  else ! SED
     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (n_photons_envoyes *  AU_to_cm * pi)
  endif

  ! Position : virtuel !!
  x = 1.0_db ! doit juste etre non nul
  y = 0.0_db 
  z = 1.0_db ! sert a rien

  ! Champ de radiation
  do zj=1,nz
     do ri=1, n_rad
        do phi_I=1,n_phi_I
           do theta_I=1,n_theta_I
              Inu(:,theta_I,phi_I,ri,zj) = sum(xI(:,theta_I,phi_I,ri,zj,:),dim=2)
           enddo ! phi_I
        enddo ! theta_I

        Inu_max(ri,zj) = maxval(Inu(1,:,:,ri,zj))
        
     enddo !zj
  enddo !ri


  ! Precalcul des directions ou on va calculer Inu * s11
  ! Permet d'etre aussi rapide que calc_Isca2
  do theta_I=1,n_theta_I
     do phi_I=1,n_phi_I
        ! Moyennage de la fct de phase sur le  bin
        do i2=1, N_super
           do i1 = 1, N_super
              ! direction de vol moyenne du bin
              f1 = real(i1) / (N_super + 1)
              f2 = real(i2) / (N_super + 1)
              
              w = 2.0_db * ((real(theta_I,kind=db) - f1) / real(n_theta_I,kind=db) ) - 1.0_db
              tab_w(i1,i2,theta_I,phi_I) =  w 
              ! if (l_sym_ima) then
              !    if (dir == 0) then
              !       phi = pi * (real(n_phi_I - phi_I + 1,kind=db) - f2) / real(n_phi_I,kind=db)  
              !    else
              !       phi = pi * (real(phi_I,kind=db) - f2) / real(n_phi_I,kind=db)  
              !    endif
              ! else
              phi = deux_pi * (real(phi_I,kind=db) - f2) / real(n_phi_I,kind=db)  
              !endif
              
              w02 = sqrt(1.0_db-w*w)
              tab_u(i1,i2,theta_I,phi_I) = w02 * cos(phi)
              tab_v(i1,i2,theta_I,phi_I) = w02 * sin(phi)
           enddo !i1
        enddo !i2
     enddo !phi_I
  enddo !theta_I



  !$omp parallel &
  !$omp default(none) &
  !$omp shared(lstrat,Inu,I_sca2,n_rad,nz,tab_s11_ray_tracing,uv0,w0,n_Stokes) &
  !$omp shared(tab_s12_o_s11_ray_tracing,tab_s33_o_s11_ray_tracing,tab_s34_o_s11_ray_tracing) &
  !$omp shared(lsepar_pola,tab_u,tab_v,tab_w,lambda,n_phi_I,n_theta_I,nang_ray_tracing,lsepar_contrib) &
  !$omp private(iscatt,id,u_ray_tracing,v_ray_tracing,w_ray_tracing,theta_I,phi_I,sum_s11,i1,i2,u,v,w,cos_scatt,sin_scatt) &
  !$omp private(sum_sin,p_ri,p_zj,ri,zj,stokes,s11,k,alloc_status,dir,correct_w,phi_scatt) &
  !$omp private(s12,s33,s34,M,ROP,RPO,v1pi,v1pj,v1pk,xnyp,costhet,theta,omega,cosw,sinw,C,D,S) 
  id = 1 ! pour code sequentiel

  if (lstrat) then
     allocate(sum_s11(n_rad,nz),s11(n_rad,nz), stat=alloc_status)
  else
     allocate(sum_s11(1,1),s11(1,1), stat=alloc_status)
  endif
  if (alloc_status > 0) then
     write(*,*) 'Allocation error sum_s11'
     stop
  endif
  sum_s11 = 0.0
  s11 = 0.0
  
  if (lsepar_pola) then
     if (lstrat) then
        allocate(s12(n_rad,nz),s33(n_rad,nz),s34(n_rad,nz), stat=alloc_status)
     else
        allocate(s12(1,1),s33(1,1),s34(1,1), stat=alloc_status)
     endif
     if (alloc_status > 0) then
        write(*,*) 'Allocation error sum_s11'
        stop
     endif
     s12 = 0.0
     s33 = 0.0
     s34 = 0.0
  endif


  ! Matrice de Mueller
  M = 0.0_db
  
  ! matrices de rotation
  ROP = 0.0_db
  RPO = 0.0_db

  RPO(1,1) = 1.0_db
  ROP(1,1) = 1.0_db
  RPO(4,4) = 1.0_db
  ROP(4,4) = 1.0_db

  ! Boucle sur les directions de ray-tracing
  do dir=0,1 
     correct_w = (2 * dir - 1)

     !$omp do 
     do iscatt = 1, nang_ray_tracing 
        !$ id = omp_get_thread_num() + 1
        phi_scatt = deux_pi * real(iscatt) / real(nang_ray_tracing)
        u_ray_tracing = uv0 * cos(phi_scatt)
        v_ray_tracing = uv0 * sin(phi_scatt)
        w_ray_tracing = w0

        do theta_I=1,n_theta_I
           do phi_I=1,n_phi_I
              ! Moyennage de la fct de phase sur le  bin
              sum_s11(:,:) = 0.0
              sum_sin = 0.0
              do i2=1, N_super
                 do i1 = 1, N_super
                    u = tab_u(i1,i2,theta_I,phi_I)
                    v = tab_v(i1,i2,theta_I,phi_I)
                    w = tab_w(i1,i2,theta_I,phi_I) * correct_w
              
                    ! Angle de diffusion --> n'utilise plus cos_thet_ray_tracing depuis super-echantillonage
                    cos_scatt = u_ray_tracing * u + v_ray_tracing * v + w_ray_tracing * w
                                        
                    k = nint(acos(cos_scatt) * real(nang_scatt)/pi)
                    if (k > nang_scatt) k = nang_scatt
                    if (k < 0) k = 0

                    sin_scatt = sqrt(1.0_db - cos_scatt*cos_scatt)
                    
                    sum_sin = sum_sin + sin_scatt
                    
                    if (lstrat) then
                       do ri=1,n_rad
                          do zj=1,nz
                             sum_s11(ri,zj) = sum_s11(ri,zj) + tab_s11_ray_tracing(lambda,ri,zj,k) * sin_scatt
                          enddo !zj
                       enddo !ri
                    else ! pas de strat
                       sum_s11(1,1) = sum_s11(1,1) + tab_s11_ray_tracing(lambda,1,1,k) * sin_scatt
                    endif

                 enddo ! i1
              enddo !i2

              if (sum_sin > 0.0) then
                 s11(:,:) = sum_s11(:,:) / sum_sin
              else
                 s11(:,:) = 0.0
              endif

              if (lsepar_pola) then ! On calcule les s12, s33, s34 et la matrice de rotation
                 
                 ! On prend le milieu du bin uniquement pour la pola, car les fct sont plus smooth
                 i1 = N_super/2 + 1 
                 i2 = i1 
                 u = tab_u(i1,i2,theta_I,phi_I)
                 v = tab_v(i1,i2,theta_I,phi_I)
                 w = tab_w(i1,i2,theta_I,phi_I) * correct_w
                 
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
                    costhet = 1.0
                 else
                    costhet = -1.0*v1pj / xnyp
                 endif
              
                 ! calcul de l'angle entre la normale et l'axe z (theta)
                 theta = acos(costhet)
                 if (theta >= pi) theta = 0.0
     
                 !     le plan de diffusion est a +ou- 90deg de la normale
                 ! Ne doit pas etre utilise en ray-tracing !!! teste OK sans
                 !theta = theta  + pi_sur_deux  
              
                 !----dans les matrices de rotation l'angle est omega = 2 * theta-----
                 omega = 2.0_db * theta ! A verifier: le moins est pour corriger un bug de signe trouve par Marshall (supprime)
                 !     prochain if car l'arccos va de 0 a pi seulement
                 !     le +/- pour faire la difference dans le sens de rotation
                 if (v1pk < 0.0) omega = -1.0_db * omega

                 ! Matrices de rotations 
                 cosw = cos(omega)
                 sinw = sin(omega)
                 if (abs(cosw) < 1e-06) cosw = 0.0_db
                 if (abs(sinw) < 1e-06) sinw = 0.0_db

                 RPO(2,2) = cosw
                 ROP(2,2) = cosw
                 RPO(2,3) = sinw
                 ROP(2,3) = -1.0_db * sinw
                 RPO(3,2) = -1.0_db * sinw
                 ROP(3,2) = sinw
                 RPO(3,3) = cosw
                 ROP(3,3) = cosw

                 if (lstrat) then
                    do ri=1,n_rad
                       do zj=1,nz
                          s12(ri,zj) = s11(ri,zj) * tab_s12_o_s11_ray_tracing(lambda,ri,zj,k)
                          s33(ri,zj) = s11(ri,zj) * tab_s33_o_s11_ray_tracing(lambda,ri,zj,k)
                          s34(ri,zj) = s11(ri,zj) * tab_s34_o_s11_ray_tracing(lambda,ri,zj,k)
                       enddo ! zj
                    enddo ! ri
                 else ! pas de strat
                    s12(1,1) = s11(1,1) * tab_s12_o_s11_ray_tracing(lambda,1,1,k)
                    s33(1,1) = s11(1,1) * tab_s33_o_s11_ray_tracing(lambda,1,1,k)
                    s34(1,1) = s11(1,1) * tab_s34_o_s11_ray_tracing(lambda,1,1,k)
                 endif ! lstrat
              endif ! lsepar_pola
              
              !write(*,*) "s12"                
                       
              ! Boucle sur les cellules pour calculer l'intensite diffusee
              p_ri = 1 ; p_zj =1
              do ri=1, n_rad
                 do zj=1,nz
                    if (lstrat) then
                       p_ri = ri
                       p_zj = zj
                    endif

                    if (lsepar_pola) then ! on calcule les s12, s33 et s34
                       ! Champ de radiation
                       stokes(:) = Inu(1:4,theta_I,phi_I,ri,zj)

                        M(1,1) = s11(p_ri,p_zj)
                        M(2,2) = s11(p_ri,p_zj)
                        M(1,2) = s12(p_ri,p_zj)
                        M(2,1) = s12(p_ri,p_zj)
                       
                        M(3,3) = s33(p_ri,p_zj)
                        M(4,4) = s33(p_ri,p_zj)
                        M(3,4) = -s34(p_ri,p_zj)
                        M(4,3) = s34(p_ri,p_zj)
                        
                        !  STOKE FINAL = RPO * M * ROP * STOKE INITIAL

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

                        I_sca2(:,iscatt,dir,ri,zj,id) = I_sca2(:,iscatt,dir,ri,zj,id) + S(:)

                     else ! lsepar_pola

                       ! Champ de radiation
                       stokes(1) = Inu(1,theta_I,phi_I,ri,zj)
                       I_sca2(1,iscatt,dir,ri,zj,id) = I_sca2(1,iscatt,dir,ri,zj,id) + s11(p_ri,p_zj) * stokes(1)
                    endif  ! lsepar_pola

                    if (lsepar_contrib) then
                       I_sca2(n_Stokes+2,iscatt,dir,ri,zj,id) = I_sca2(n_Stokes+2,iscatt,dir,ri,zj,id) + &
                               s11(p_ri,p_zj) * Inu(n_Stokes+2,theta_I,phi_I,ri,zj)
                       
                       I_sca2(n_Stokes+4,iscatt,dir,ri,zj,id) = I_sca2(n_Stokes+4,iscatt,dir,ri,zj,id) + &
                               s11(p_ri,p_zj) * Inu(n_Stokes+4,theta_I,phi_I,ri,zj)
                    endif ! lsepar_contrib
                    

                 enddo ! zj
              enddo !ri

              !write(*,*) "Isca"

           enddo ! phi_I
        enddo ! theta_I

     enddo ! iscatt
     !$omp enddo
  enddo ! dir
  deallocate(s11,sum_s11)
  if (lsepar_pola) deallocate(s12,s33,s34)
  !$omp end parallel
  


  ! Normalisation
  ! Boucle sur les cellules
  do ri=1, n_rad
     do zj=1,nz
        facteur = energie_photon / volume(ri) 
        I_sca2(:,:,:,ri,zj,:) =  I_sca2(:,:,:,ri,zj,:) *  facteur * kappa_sca(lambda,ri,zj,1)
     enddo
  enddo

  deallocate(Inu)
  return

end subroutine calc_Isca2_new

!***********************************************************

subroutine calc_Isca2_star(lambda,ibin)
  ! Version acceleree !TMP
  ! Routine liee a la precedente

  integer, intent(in) :: lambda, ibin

  integer :: ri, zj, p_ri, p_zj, dir, iscatt, id, dir_I
  real(kind=db), dimension(4) :: Stokes, S, C, D
  real(kind=db) :: x, y, z, u, v, w, phi, w02, facteur

  real(kind=db), dimension(n_theta_I,n_phi_I) :: tab_u, tab_v, tab_w
  real(kind=db), dimension(n_rad,nz) :: Inu
  real(kind=db), dimension(4,4) ::  M, ROP, RPO

  integer :: k, kmin, kmax, alloc_status
  real :: frac, frac_m1, s11, s12, s33, s34, cos_scatt
  real(kind=db) :: omega, sinw, cosw, norme, energie_photon, n_photons_envoyes

  real(kind=db), parameter :: prec = 0._db
  real :: angle

  id = 1 ! TMP

  
  if (lmono0) then ! image
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (real(nbre_photons_loop)*real(nbre_photons_image) *  AU_to_cm * pi)
  else ! SED
     n_photons_envoyes = sum(n_phot_envoyes(lambda,:))
     energie_photon = (E_stars(lambda) + E_disk(lambda)) * tab_lambda(lambda) * 1.0e-6  / &
          (n_photons_envoyes *  AU_to_cm * pi)
  endif

  ! Matrice de Mueller
  M = 0.0_db
  
  ! matrices de rotation
  ROP = 0.0_db
  RPO = 0.0_db

  RPO(1,1) = 1.0_db
  ROP(1,1) = 1.0_db
  RPO(4,4) = 1.0_db
  ROP(4,4) = 1.0_db

  ! Champ de radiation
  do ri=1, n_rad
     do zj=1,nz
        Inu(ri,zj) = sum(xI_star(ri,zj,:))
     enddo !zj
  enddo !ri

  stokes(:) = 0.0_db
  p_ri = 1
  p_zj = 1

  ! Boucle sur les cellules
  do ri=1, n_rad
     do zj=1,nz
        ! Champ de radiation
        stokes(1) = Inu(ri,zj)

        if (stokes(1) < 1.e-30_db) cycle

        ! Direction de vol 
        x = r_grid(ri,zj) ! doit juste etre non nul
        y = 0.0_db
        z = z_grid(ri,zj) ! sert a rien

        norme = sqrt(x**2 + z**2)
        u = x / norme
        v = 0.0_db
        w = z / norme

        !! Ca change quasi rien : rassurant !!!!!
!        if (stokes(1) > 0.0_db) then
!           w = sum(xw_star(ri,zj,:)) / stokes(1) 
!           u = sqrt(1.0_db - w*w)
!           v = 0.0_db
!        else
!           cycle
!        endif

        ! cos_scatt_ray_tracing correspondants
        call angles_scatt_ray_tracing(id,ibin,x,y,z,u,v,w,.true.)

        if (lstrat) then
           p_ri = ri
           p_zj = zj
        endif

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
        !!$      frac_m1 = 1.0_db - frac

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
                 if (abs(cosw) < 1e-06) cosw = 0.0_db
                 if (abs(sinw) < 1e-06) sinw = 0.0_db
                 
                 RPO(2,2) = COSW
                 ROP(2,2) = COSW
                 RPO(2,3) = SINW
                 ROP(2,3) = -1.0_db * SINW
                 RPO(3,2) = -1.0_db * SINW
                 ROP(3,2) = SINW
                 RPO(3,3) = COSW
                 ROP(3,3) = COSW
                 
                 ! Matrice de Mueller
                 s11 = tab_s11_ray_tracing(lambda,p_ri,p_zj,k)
                 s12 = tab_s12_ray_tracing(lambda,p_ri,p_zj,k)
                 s33 = tab_s33_ray_tracing(lambda,p_ri,p_zj,k)
                 s34 = tab_s34_ray_tracing(lambda,p_ri,p_zj,k)
                 
                 M(1,1) = s11
                 M(2,2) = s11
                 M(1,2) = s12
                 M(2,1) = s12
                 
                 M(3,3) = s33
                 M(4,4) = s33
                 M(3,4) = -s34
                 M(4,3) = s34
                 M(1,1) = s11
                       
                       
                 !  STOKE FINAL = RPO * M * ROP * STOKE INITIAL
                 
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
                 
                 eps_dust2_star(:,iscatt,dir,ri,zj) =  eps_dust2_star(:,iscatt,dir,ri,zj) + S(:)
              else ! .not.lsepar_pola
                 s11 = tab_s11_ray_tracing(lambda,p_ri,p_zj,k)
                 eps_dust2_star(1,iscatt,dir,ri,zj) =  eps_dust2_star(1,iscatt,dir,ri,zj) + s11 * stokes(1)
              endif ! lsepar_pola
        
           enddo ! iscatt
        enddo ! dir

     enddo ! zj
  enddo ! ri

  ! Normalisation
  ! Boucle sur les cellules
  do ri=1, n_rad
     do zj=1,nz
        facteur = energie_photon / volume(ri) 
        eps_dust2_star(:,:,:,ri,zj) =  eps_dust2_star(:,:,:,ri,zj) *  facteur * kappa_sca(lambda,ri,zj,1) 
     enddo
  enddo

  return

end subroutine calc_Isca2_star

!***********************************************************

function dust_source_fct(lambda,ri,zj,x,y,z)
  ! La direction du rayon est maintenant fixee le long de l'axe x
  ! l'angle de diffusion ne depend que de la position x, y, z

  integer, intent(in) :: lambda,ri, zj
  real(kind=db), intent(in) :: x, y, z

  real(kind=db), dimension(N_type_flux) :: dust_source_fct, SF1, SF2, SF3, SF4

  real(kind=db) :: phi_pos, frac, un_m_frac, xiscatt, xi_az, frac_r, frac_z, r
  integer :: iscatt1, iscatt2, dir, i_az, i_az_p1, psup, ri1, zj1, ri2, zj2
  integer :: phi_k, phi_k_p1, n_pola

  SF1 = 0 ; SF2 = 0 ; SF3 = 0 ; SF4 = 0

!  write(*,*) "test"


  ! Interpolations angulaires
  if (lscatt_ray_tracing1) then
     if (z > 0.0_db) then
        psup=1
     else
        psup=0
     endif

     phi_pos = atan2(y,x) 

     ! OK : lisse les images par rapport a la methode en dessous
!---     xi_az =  modulo(phi_pos, deux_pi) / deux_pi * n_az_rt + 0.5
!---     phi_k =  floor(xi_az)
!---     frac = xi_az - phi_k
!---     un_m_frac = 1.0_db - frac
!---
!---     phi_k_p1 = phi_k + 1
!---     if (phi_k <=0) phi_k = 1 
!---     if (phi_k_p1 > n_az_rt) phi_k_p1=n_az_rt 
!---
!---     dust_source_fct(:) = eps_dust1(:,ri,zj,phi_k_p1,psup) * frac + eps_dust1(:,ri,zj,phi_k,psup) * un_m_frac

     ! Cette methode est OK mais on voit les cellules
     phi_k =  floor(modulo(phi_pos, deux_pi) / deux_pi * n_az_rt) + 1
     if (phi_k > n_az_rt) phi_k = n_az_rt
     dust_source_fct(:) = eps_dust1(:,ri,zj,phi_k,psup)  ! ??? ce n'est pas lineaire

  else ! Methode 2 : la seule actuelle
     ! Pour interpolations spatiales
     r = sqrt(x*x + y*y)

     if (r > r_grid(ri,zj)) then
        ri1 = ri
        ri2 = ri + 1
     else
        ri1 = ri - 1
        ri2 = ri
     endif
     
     if (abs(z) > z_grid(ri,zj)) then
        zj1 = zj
        zj2 = zj + 1
     else
        zj1 = zj - 1
        zj2 = zj
     endif

     if (ri2 > n_rad) then
        ri2 = n_rad
        frac_r = 0._db
     else  if (ri1 < 1) then
        ri1 = 1
        frac_r = 0._db
     else
        frac_r = (log(r_grid(ri2,zj)) - log(r)) / (log(r_grid(ri2,zj)) - log(r_grid(ri1,zj))) 
     endif
  
     if (zj2 > nz) then
        zj2 = nz
        frac_z = 1.0_db
     else  if (zj1 < 1) then
        zj1 = 1
        frac_z = 1.0_db
     else
        frac_z = (z_grid(ri,zj2) - abs(z)) / (z_grid(ri,zj2) - z_grid(ri,zj1)) 
     endif


     ! Calcul angle phi a la position
     phi_pos = modulo(atan2(y,x) + deux_pi,deux_pi)

     if (z > 0.0_db) then
        dir=1
     else
        dir=0
     endif

     !----------------------------------------------------
     ! Emissivite du champ diffuse plusieurs fois
     !----------------------------------------------------
     xiscatt =   max(phi_pos/deux_pi  * nang_ray_tracing, 0.0_db) 
     
     iscatt1 = floor(xiscatt)
     frac = xiscatt - iscatt1
     un_m_frac = 1.0_db - frac
     iscatt2 = iscatt1 + 1

     ! Les limites periodiques
     !    if (iscatt2 > nang_ray_tracing) iscatt2 = 1
     !    if (iscatt1 < 1) iscatt1 = nang_ray_tracing
     iscatt1 = modulo(iscatt1,nang_ray_tracing)
     if (iscatt1==0) iscatt1 = nang_ray_tracing

     iscatt2 = modulo(iscatt2,nang_ray_tracing)
     if (iscatt2==0) iscatt2 = nang_ray_tracing

     ! Fct source des cellules
     SF1(:) = eps_dust2(:,iscatt2,dir,ri1,zj1) * frac + eps_dust2(:,iscatt1,dir,ri1,zj1) * un_m_frac 
     SF2(:) = eps_dust2(:,iscatt2,dir,ri2,zj1) * frac + eps_dust2(:,iscatt1,dir,ri2,zj1) * un_m_frac
     SF3(:) = eps_dust2(:,iscatt2,dir,ri1,zj2) * frac + eps_dust2(:,iscatt1,dir,ri1,zj2) * un_m_frac
     SF4(:) = eps_dust2(:,iscatt2,dir,ri2,zj2) * frac + eps_dust2(:,iscatt1,dir,ri2,zj2) * un_m_frac

     !----------------------------------------------------
     ! Emissivite du champ stellaire diffuse 1 fois
     !----------------------------------------------------
     xiscatt =   max(phi_pos/deux_pi  * nang_ray_tracing_star, 0.0_db) 
     
     iscatt1 = floor(xiscatt)
     frac = xiscatt - iscatt1
     un_m_frac = 1.0_db - frac
     iscatt2 = iscatt1 + 1

     ! Les limites periodiques
     iscatt1 = modulo(iscatt1,nang_ray_tracing_star)
     if (iscatt1==0) iscatt1 = nang_ray_tracing_star

     iscatt2 = modulo(iscatt2,nang_ray_tracing_star)
     if (iscatt2==0) iscatt2 = nang_ray_tracing_star

      
     ! Fct source des cellules
!     SF1(:) = SF1(:) + exp(log(eps_dust2_star(:,iscatt2,dir,ri1,zj1)) * frac + log(eps_dust2_star(:,iscatt1,dir,ri1,zj1)) * un_m_frac)  
     ! TODO : passer en log ameliore un peu la forme a longue distance mais cree des bug si 0. !!

     if (lsepar_pola) then
        n_pola = 4
     else
        n_pola = 1
     endif
     
     ! TODO : petit bug ici --> on peut ajuster la taille de eps_dust2_star au lieu de la fixer a 4
     SF1(1:n_pola) = SF1(1:n_pola) + eps_dust2_star(1:n_pola,iscatt2,dir,ri1,zj1) * frac &
          + eps_dust2_star(1:n_pola,iscatt1,dir,ri1,zj1) * un_m_frac 
     SF2(1:n_pola) = SF2(1:n_pola) + eps_dust2_star(1:n_pola,iscatt2,dir,ri2,zj1) * frac &
          + eps_dust2_star(1:n_pola,iscatt1,dir,ri2,zj1) * un_m_frac
     SF3(1:n_pola) = SF3(1:n_pola) + eps_dust2_star(1:n_pola,iscatt2,dir,ri1,zj2) * frac &
          + eps_dust2_star(1:n_pola,iscatt1,dir,ri1,zj2) * un_m_frac
     SF4(1:n_pola) = SF4(1:n_pola) + eps_dust2_star(1:n_pola,iscatt2,dir,ri2,zj2) * frac &
          + eps_dust2_star(1:n_pola,iscatt1,dir,ri2,zj2) * un_m_frac

     if (lsepar_contrib) then
         SF1(n_pola+2) = SF1(n_pola+2) + eps_dust2_star(1,iscatt2,dir,ri1,zj1) * frac &
              + eps_dust2_star(1,iscatt1,dir,ri1,zj1) * un_m_frac 
         SF2(n_pola+2) = SF2(n_pola+2) + eps_dust2_star(1,iscatt2,dir,ri2,zj1) * frac &
              + eps_dust2_star(1,iscatt1,dir,ri2,zj1) * un_m_frac 
         SF3(n_pola+2) = SF3(n_pola+2) + eps_dust2_star(1,iscatt2,dir,ri1,zj2) * frac &
              + eps_dust2_star(1,iscatt1,dir,ri1,zj2) * un_m_frac 
         SF4(n_pola+2) = SF4(n_pola+2) + eps_dust2_star(1,iscatt2,dir,ri2,zj2) * frac &
              + eps_dust2_star(1,iscatt1,dir,ri2,zj2) * un_m_frac 
     endif

     ! Fct source interpolee
     dust_source_fct(:) = &
          frac_r * frac_z * SF1(:) + &
          (1.0_db - frac_r) * frac_z * SF2(:) + &
          frac_r * (1.0_db - frac_z) * SF3(:) + &
          (1.0_db - frac_r) * (1.0_db - frac_z) * SF4(:) 

  endif !lscatt_ray_tracing

  return
  
end function dust_source_fct

!***********************************************************

end module dust_ray_tracing



