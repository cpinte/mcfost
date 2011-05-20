module dust

  use parametres
  use grains
  use disk, only : n_zones
  use em_th
  use constantes
  use opacity
  use ray_tracing
  use dust_ray_tracing

  use scattering
  use input
  use output

  implicit none

  contains

subroutine taille_grains()
! Distribution en taille des grains
! Sortie de prop_grains pour cas multi-lambda
! C. Pinte 14/01/05  

  implicit none
  integer :: k, i, j
  real :: a, alfa, qext=0.0, qsca=0.0, M_tot, nbre_tot_grains
  real(kind=db) :: exp_grains
  real :: masse_pop

  type(dust_pop_type), pointer :: dp

!* --- masse moyenne d'un grain ( XMG ) ---
!*     Volume = 4/3*pi * a**3 * 1E-12  ou a = rayon en micron
!*     Rho1g = about 3 g/cm3
!*
!*      XMG = 4.1887E-12 * A**3 * RHO1G
!*
!* Avec une distribution en taille des grains, la formule est
!* modifiee: il faut calculer la masse moyenne des grains, ponderee
!* suivant la distribution. Avec une distribution du type n(a) = a**(-p),
!* cela devient:

  

  ! Boucle sur les populations de grains
  do i=1, n_pop
     dp => dust_pop(i)

     if (dp%aexp < 0) then
        write(*,*) "****************************************"
        write(*,*) "Warning: slope grains size negative !!!!"
        write(*,*) "****************************************"
     endif

     if (abs(dp%amin - dp%amax) < 1.0e-5 * dp%amax) then
        a=dp%amin
        dp%xmg = 4.1887E-12 * a**3 * dp%rho1g
     else
        if (abs(dp%aexp - 4.) > 1.0e-5)  then
           if (abs(dp%aexp - 1.) > 1.0e-5) then
              dp%xmg = 4.1887E-12 * dp%rho1g *(1-dp%aexp)/(4-dp%aexp) *(dp%amax**(4-dp%aexp)-dp%amin**(4-dp%aexp)) / &
                   (dp%amax**(1-dp%aexp)-dp%amin**(1-dp%aexp))
           else
              dp%xmg = 4.1887E-12 * dp%rho1g /(4-dp%aexp) *(dp%amax**(4-dp%aexp)-dp%amin**(4-dp%aexp)) / &
                   (log(dp%amax)-log(dp%amin))
           endif
        else
           dp%xmg = 4.1887E-12 * dp%rho1g *(1-dp%aexp)*(log(dp%amax)-log(dp%amin)) / &
                (dp%amax**(1-dp%aexp)-dp%amin**(1-dp%aexp))
        endif
     endif
 

     ! Proprietes des grains
     !exp_grains = (amax/amin)**(1./real(n_grains_tot))

     if ((dp%n_grains==1).and.(abs(dp%amax-dp%amin) > 1.0e-3 * dp%amin)) then
        write(*,*) "You have specified 1 grain size but amin != amax. Are you sure ?"
        write(*,*) "If yes, press return"
        read(*,*)
     endif

     exp_grains =  exp((1.0_db/real(dp%n_grains,kind=db)) * log(dp%amax/dp%amin))
     a = dp%amin*sqrt(exp_grains)
     tab_a(dp%ind_debut) = a
     nbre_grains(dp%ind_debut) = a**(-dp%aexp) * a
     grain(dp%ind_debut)%methode_chauffage = dp%methode_chauffage
     grain(dp%ind_debut)%zone = dp%zone
     grain(dp%ind_debut)%pop = i
        
     ! Taille des grains (recursif)
     masse_pop = nbre_grains(dp%ind_debut)
     do  k=dp%ind_debut+1, dp%ind_fin
        a= tab_a(k-1) * exp_grains
        tab_a(k) = a
        ! Multiplication par a car da = a.dln(a)
        nbre_grains(k) = a**(-dp%aexp) * a
        grain(k)%methode_chauffage = dp%methode_chauffage
        grain(k)%zone = dp%zone
        grain(k)%pop = i
        masse_pop = masse_pop + nbre_grains(k)
     enddo !k

     masse_pop = masse_pop * dp%xmg

     ! Normalisation du nombre de grains pour atteindre la bonne masse
     nbre_grains(dp%ind_debut:dp%ind_fin) = nbre_grains(dp%ind_debut:dp%ind_fin) * dp%masse/masse_pop

     ! Normalisation de tous les grains au sein d'une pop
     nbre_tot_grains = 0.0
     do k=dp%ind_debut,dp%ind_fin
        nbre_tot_grains =  nbre_tot_grains + nbre_grains(k)
     enddo
     
     ! Fraction de grains de taille k au sein d'une pop
     do  k=dp%ind_debut,dp%ind_fin
        nbre_grains(k) = nbre_grains(k)/nbre_tot_grains
     enddo !k

  enddo !i

  return

end subroutine taille_grains

!********************************************************************

subroutine init_indices_optiques()
  ! calcule les tableaux tab_amu1 et tab_amu2
  ! aux longueurs d'onde de tab_lambda
  ! C. Pinte
  ! 18/10/07 (separation depuis init_lambda)

  implicit none

  complex :: m1, m2, mavg, eps1, eps2, eavg1, eavg2, eavg, a, b, c, delta
  real :: frac, f1, f2, fbuffer
  integer :: n, i, k, syst_status, alloc_status, pop, status, n_ind, buffer, ios, n_comment

  character(len=512) :: filename

  real, dimension(:), allocatable :: tab_l, tab_a1, tab_a2

  logical :: l_ordre_decroissant

  do pop=1, n_pop
     if (.not.dust_pop(pop)%is_PAH) then 

        ! Lecture fichier indices
        filename = trim(dust_dir)//trim(dust_pop(pop)%indices)

        open(unit=1,file=filename, status='old', iostat=ios)
        if (ios /=0) then
           write(*,*) "ERROR: dust file cannot be opened:",trim(filename)
           write(*,*) "Exiting"
           stop
        endif

        ! On elimine les lignes avec des commentaires
        status = 1
        n_comment = 0 
        do while (status /= 0)
           n_comment = n_comment + 1
           read(1,*,iostat=status) fbuffer
        enddo
        n_comment = n_comment - 1

        ! On compte les lignes avec des donnees
        status=0
        n_ind=1 ! On a deja lu une ligne en cherchant les commentaires
        do while(status==0)
           n_ind=n_ind+1
           read(1,*,iostat=status)
        enddo
        n_ind = n_ind - 1

        ! On enleve les 2 premieres lignes
        n_ind = n_ind - 2 
        
        ! Allocation de tab
        allocate(tab_l(n_ind), tab_a1(n_ind), tab_a2(n_ind), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error zsup'
           stop
        endif
        tab_l=0.0 ; tab_a1=0.0 ; tab_a2=0.0 
        
        
        ! Lecture proprement dite
        rewind(1)
        ! On passe les commentaires
        do i=1, n_comment
           read(1,*)
        enddo

        ! lecture densite
        read(1,*) dust_pop(pop)%rho1g
        if (dust_pop(pop)%rho1g > 10.) then
           write(*,*) "ERROR: optical indices file has the wrong format"
           write(*,*) "Exiting"
           stop
        endif
        ! ligne vide
        read(1,*)

        ! Lecture indices
        do i=1,n_ind
           read(1,*) tab_l(n_ind-i+1), tab_a1(n_ind-i+1), tab_a2(n_ind-i+1)
        enddo
           
        close(unit=1)

        if (tab_l(2) > tab_l(1)) then
           l_ordre_decroissant = .false. 
        else
           l_ordre_decroissant = .true.
        endif

        if (dust_pop(pop)%porosity > tiny_real) then 
           f2 = dust_pop(pop)%porosity
           f1 = 1. - f2

           ! On corrige la densite 
           dust_pop(pop)%rho1g  = dust_pop(pop)%rho1g * f1

           ! On corrige les indices optiques
           do i=1, n_ind
              m1 = cmplx(tab_a1(i),tab_a2(i))
              m2 = cmplx(1.0,0.0) ! Vide
              eps1 = m1*m1
              eps2 = m2*m2
              
              ! Effective medium theory (Bruggeman rule)
              b = 2*f1*eps1 -f1*eps2 + 2*f2*eps2 -f2*eps1
              c = eps1 * eps2
              a = -2.   ! --> 1/2a  = -0.25
              delta = b*b - 4*a*c 
              eavg1 = (-b - sqrt(delta)) * (-0.25)
              eavg2 = (-b + sqrt(delta)) * (-0.25) 

              ! Selection (et verif) de l'unique racine positive (a priori c'est eavg1)
              if (aimag(eavg1) > 0) then
                 if (aimag(eavg2) > 0) then
                    write(*,*) "Error in Bruggeman EMT rule !!!"
                    write(*,*) "All imaginary parts are > 0"
                    write(*,*) "Exiting."
                    stop
                 endif
                 eavg = eavg1
              else
                 if (aimag(eavg2) < 0) then
                    write(*,*) "Error in Bruggeman EMT rule !!!"
                    write(*,*) "All imaginary parts are < 0"
                    write(*,*) "Exiting."
                    stop
                 endif
                 eavg = eavg2
              endif

              mavg = sqrt(eavg)
         !     write(*,*) "verif Bruggeman", (abs(f1 * (eps1 - eavg)/(eps1 + 2*eavg) + &
         !     f2 * (eps2 - eavg)/(eps2 + 2*eavg))) ;

              ! Moyennage a la Mathis & Whieffen 89
              ! Valeur max pour indices Voschnikov 2005, A&A, 429, 371 : p380
          !    eavg = f1 * eps1 + f2 * eps2
          !    mavg = sqrt(eavg)
              
              ! Valeur min pour indices Voschnikov 2005, A&A, 429, 371 : p380
              ! eavg = 1.0 / (f1/eps1 + f2/eps2)
              ! mavg = sqrt(eavg)


              ! Evite la chute d'opacite a grand lambda : proche M&W 
              ! Opcaite un peu plus faible dans le mm
              ! mavg = f1 * m1 + f2 * m2 

              tab_a1(i) = real(mavg)
              tab_a2(i) = aimag(mavg)
           
 !             write(*,*) tab_lambda(i), tab_amu1(i,pop), tab_amu2(i,pop)

           enddo ! n_ind

        endif ! porosity

        ! Calcul des indices par interpolation lineaire
        do i=1, n_lambda
           if (l_ordre_decroissant) then
              k=2 ! deplace ici car les lambda ne sont plus dans l'ordre pour l'emission moleculaire
              do while((tab_lambda(i) < tab_l(k)).and.(k <= n_ind-1))
                 k=k+1
              enddo
           
              frac=(log(tab_l(k))-log(tab_lambda(i)))/(log(tab_l(k))-log(tab_l(k-1)))
              tab_amu1(i,pop)=exp(log(tab_a1(k-1))*frac+log(tab_a1(k))*(1.0-frac))
              tab_amu2(i,pop)=exp(log(tab_a2(k-1))*frac+log(tab_a2(k))*(1.0-frac))
              
           else ! ordre croissant
              k=2 ! deplace ici car les lambda ne sont plus dans l'ordre pour l'emission moleculaire
              do while((tab_lambda(i) > tab_l(k)).and.(k <= n_ind-1))
                 k=k+1
              enddo
           
              frac=(log(tab_l(k))-log(tab_lambda(i)))/(log(tab_l(k))-log(tab_l(k-1)))
              tab_amu1(i,pop)=exp(log(tab_a1(k-1))*frac+log(tab_a1(k))*(1.0-frac))
              tab_amu2(i,pop)=exp(log(tab_a2(k-1))*frac+log(tab_a2(k))*(1.0-frac))
           endif ! croissant ou decroissant
           
        enddo ! lambda
        
        deallocate(tab_l,tab_a1,tab_a2)
     else ! notPAH       
        dust_pop(pop)%rho1g = 2.5
     endif ! notPAH       
     
  enddo ! pop

  return

end subroutine init_indices_optiques

!******************************************************************************

subroutine prop_grains(lambda, p_lambda)
  ! Calcule les propriétés des grains
  ! Masse, fraction en nombre, sections efficaces, matrices de muellers
  ! sortie des routines opacite
  ! C. Pinte 7/01/05

  implicit none

  integer, intent(in) :: lambda, p_lambda
  real, parameter :: pi = 3.1415926535
  real :: a, alfa, qext, qsca, fact, gsca, amu1, amu2
  integer :: k, alloc_status, i, pop, l

  type(dust_pop_type) :: dp

  ! PAH
  real, dimension(PAH_n_lambda,PAH_n_rad) :: tmp_Q_ext, tmp_Q_abs, tmp_Q_sca, tmp_g
  character(len=512) :: opt_file
  real, dimension(PAH_n_lambda) :: tmp_PAH_lambda
  real, dimension(PAH_n_rad) :: tmp_PAH_rad


  qext=0.0
  qsca=0.0

  ! Longueur d'onde
  wavel=tab_lambda(lambda)

  if (lnRE) then
     if (.not.lread_PAH) then ! variable logique pour ne lire qu'une seule fois le fichier de prop de PAH
        allocate(PAH_Q_ext(PAH_n_lambda,PAH_n_rad,n_pop), &
             PAH_Q_sca(PAH_n_lambda,PAH_n_rad,n_pop), &
             PAH_g(PAH_n_lambda,PAH_n_rad,n_pop), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error PAH_Q_ext'
           stop
        endif
        
        do i=1, n_pop
           dp = dust_pop(i)
           if (dp%is_PAH) then 
              opt_file = trim(dust_dir)//"/"//trim(dp%indices)
              ! load optical data from file
              call draine_load(opt_file, PAH_n_lambda, PAH_n_rad, 10, 1, &
                   tmp_PAH_lambda, tmp_PAH_rad,  tmp_Q_ext, tmp_Q_abs, tmp_Q_sca, tmp_g, 4)
              PAH_Q_ext(:,:,i) = tmp_Q_ext
              PAH_Q_sca(:,:,i) = tmp_Q_sca
              PAH_g(:,:,i) = tmp_g
              PAH_lambda = tmp_PAH_lambda
              log_PAH_rad = log(tmp_PAH_rad)
           endif
        enddo !i
        
        ! abs car le fichier de lambda pour les PAHs est a l'envers
        PAH_delta_lambda(1) = abs(PAH_lambda(2) - PAH_lambda(1))
        if (n_lambda > 1) then
           PAH_delta_lambda(n_lambda) = abs(PAH_lambda(n_lambda) - PAH_lambda(n_lambda-1))
        endif
        do l=2,PAH_n_lambda-1
           PAH_delta_lambda(l) = 0.5* abs(PAH_lambda(l+1) - PAH_lambda(l-1))
        enddo

        lread_PAH = .true.
     endif
  endif

  ! Prop optiques
  !$omp parallel &
  !$omp default(none) &
  !$omp private(k,a,alfa,qext,qsca,fact,gsca,amu1,amu2,pop) &
  !$omp shared(tab_a,q_ext,q_sca,q_abs,wavel,aexp,tab_albedo,lambda,p_lambda,tab_g,grain) &
  !$omp shared(laggregate,tab_amu1,tab_amu2,n_grains_tot,is_pop_PAH,is_grain_PAH) &
  !$omp shared(dust_pop) 
  !$omp do schedule(dynamic,10)
  ! on fait la boucle a l'envers pour optimiser la parallelisation
  ! et savoir des le debut si l'alloc. mem. ds bhmie passe ou pas
  do  k=1,n_grains_tot
     pop = grain(k)%pop
     a = tab_a(k)
     if (.not.dust_pop(pop)%is_PAH) then ! theorie de mie ou gmm
        amu1=tab_amu1(lambda,pop)
        amu2=tab_amu2(lambda,pop)
        alfa = 2.0 * pi * a / wavel
        if (laggregate) then
           call mueller_gmm(p_lambda,k,alfa,qext,qsca,gsca)
        else
           call mueller2(p_lambda,k,alfa,amu1,amu2,qext,qsca,gsca)
        endif
     else ! grain de PAH
        is_grain_PAH(k) = .true.
        call mueller_PAH(lambda,p_lambda,k,qext,qsca,gsca)
     endif
     tab_albedo(lambda,k)=qsca/qext
     tab_g(lambda,k) = gsca
     ! tau est sans dimension : [kappa * lvol = density * a² * lvol] 
     ! a² microns² -> 1e-8 cm²             \
     ! density en cm-3                      > reste facteur 149595.0 
     ! longueur de vol en AU = 1.5e13 cm   / 
     fact =  pi * a * a * 149595.0
     !q_geo(k) = pi * a * a * 1.e-12 ! en m^2
     q_ext(lambda,k) = qext * fact
     q_sca(lambda,k) = qsca * fact
     q_abs(lambda,k) = q_ext(lambda,k) - q_sca(lambda,k)     
  enddo !k
  !$omp enddo
  !$omp end parallel

  return

end subroutine prop_grains

!******************************************************************************

subroutine opacite2(lambda)
! Calcule la table d'opacite et probsizecumul
! Inclus stratification empirique
! Utilise les resultats des routine densite et prop_grains
! Doit etre utilise juste apres prop_grain : lambda ne doit pas changer entre 2 
! (dans le cas multi-longueurs d'onde) 
! update : 12/09/06

  implicit none

  integer, intent(in) :: lambda

  real, parameter :: G = 6.672e-8

  integer :: i,j, pk, k, ii, jj, kk, l, k_min, thetaj
  real(kind=db) ::  density, proba, q_sca_tot, q_ext_tot,norme
  logical :: lcompute_obs
  
  ! Pour spline
  real, dimension(1:n_prob) :: xa,ya,y2
  real :: delta_prob1, delta_probn, yp1, ypn

  real :: yp_p1, yp_m1, ypp, somme

  real :: z, rcyl

  logical :: ldens0

  real(kind=db) :: kappa_abs_RE, kappa_abs_tot, angle

  ! Attention : dans le cas no_strat, il ne faut pas que la cellule (1,1,1) soit vide.
  ! on la met à nbre_grains et on effacera apres
  ! c'est pour les prop de diffusion en relatif donc la veleur exacte n'a pas d'importante
  ldens0 = .false.
  if (.not.lstrat) then
     if (maxval(densite_pouss(1,1,1,:)) < tiny_real) then
        ldens0 = .true.
        densite_pouss(1,1,1,:) = nbre_grains(:)
     endif
  endif
        
  if (lmono0) then
     lcompute_obs = .true.
  else
     if (lsed_complete) then
        lcompute_obs = .true.
     else
        lcompute_obs = .not.letape_th
     endif
  endif

!***************
! Calcul proba cumulee en dessous d'une taille de grains 
! et opacite pour chaque position
!
!* probsizecumul(i) represente la probabilite cumulee en-dessous d'une
!* certaine taille de grain. Ce tableau est utilise pour le tirage
!* aleatoire de la taille du grain diffuseur, puisqu'elle doit prendre
!* en compte le nombre de grains en meme temps que leur probabilite
!* individuelle de diffuser (donnee par qsca*pi*a**2). 


  if (scattering_method == 2) then
     tab_s11_pos(lambda,:,:,:,:)=0.0
     if (lsepar_pola) then
        tab_s12_pos(lambda,:,:,:,:)=0.0
        tab_s33_pos(lambda,:,:,:,:)=0.0
        tab_s34_pos(lambda,:,:,:,:)=0.0
     endif
  endif

  ! Calcul opacite et probabilite de diffusion
  do i=1, n_rad
     bz : do j=j_start,nz
        if (j==0) cycle bz
        do pk=1, n_az
           kappa(lambda,i,j,pk) = 0.0
           kappa_abs_tot = 0.0
           kappa_abs_RE = 0.0
           do  k=1,n_grains_tot
              density=densite_pouss(i,j,pk,k)
              kappa(lambda,i,j,pk) = kappa(lambda,i,j,pk) + q_ext(lambda,k) * density
              kappa_abs_tot = kappa_abs_tot + q_abs(lambda,k) * density
           enddo !k

           if (lRE_LTE) then 
              kappa_abs_eg(lambda,i,j,pk) = 0.0
              do k=grain_RE_LTE_start,grain_RE_LTE_end
                 density=densite_pouss(i,j,pk,k)
                 kappa_abs_eg(lambda,i,j,pk) =  kappa_abs_eg(lambda,i,j,pk) + q_abs(lambda,k) * density
                 kappa_abs_RE = kappa_abs_RE + q_abs(lambda,k) * density
              enddo
            endif

           if (lRE_nLTE) then 
              do k=grain_RE_nLTE_start,grain_RE_nLTE_end
                 density=densite_pouss(i,j,pk,k)
                 kappa_abs_RE = kappa_abs_RE + q_abs(lambda,k) * density
              enddo
           endif

           if (lcompute_obs.and.lscatt_ray_tracing) then
              kappa_sca(lambda,i,j,pk) = 0.0
              do k=1,n_grains_tot
                 density=densite_pouss(i,j,pk,k)
                 kappa_sca(lambda,i,j,pk) =  kappa_sca(lambda,i,j,pk) + q_sca(lambda,k) * density
              enddo
           endif
       
           if (kappa_abs_tot > tiny_db) proba_abs_RE(lambda,i,j,pk) = kappa_abs_RE/kappa_abs_tot
        enddo !pk
     enddo bz !j
  enddo !i

  proba_abs_RE(lambda,:,nz+1,:) = proba_abs_RE(lambda,:,nz,:)
  
  ! proba absorption sur une taille donnée
  if (lRE_nLTE) then
     if (.not.lmono) then
        do i=1, n_rad
           do j=1, nz
              prob_kappa_abs_1grain(lambda,i,j,0)=0.0  
              do  k=1, n_grains_RE_nLTE
                 density=densite_pouss(i,j,1,k)
                 prob_kappa_abs_1grain(lambda,i,j,k)=prob_kappa_abs_1grain(lambda,i,j,k-1) + &
                      q_abs(lambda,k) * density
              enddo !k
              prob_kappa_abs_1grain(lambda,i,j,:) =  prob_kappa_abs_1grain(lambda,i,j,:)/&
                   prob_kappa_abs_1grain(lambda,i,j,n_grains_RE_nLTE) 
           enddo !j
        enddo !i
     endif !.not.lmono
  endif !lnLTE


  if (lkappa_abs_grain) then
     if (.not.lmono) then
        write(*,*) "Error : lkappa_abs_grain only monochrmatic"
        stop
     endif
     
     allocate(kappa_abs_tmp(n_rad,nz,n_grains_tot))

     kappa_abs_tmp =0._db
     do i=1, n_rad
        do j=1,nz
           pk=1
           do  k=1,n_grains_tot
              kappa_abs_tmp(i,j,k) =  q_abs(lambda,k) * densite_pouss(i,j,pk,k)
           enddo !k
        enddo
     enddo

     write(*,*) "Writing kappa_abs_grain.fits.gz"
     call cfitsWrite("kappa_abs_grain.fits.gz",kappa_abs_tmp,shape(kappa_abs_tmp))
     deallocate(kappa_abs_tmp)
     write(*,*) "Done"
     stop
  endif !lkappa_abs_grain

  if (lscatt_ray_tracing) then
     do thetaj=0,nang_scatt
        angle = real(thetaj)/real(nang_scatt)*pi
        tab_cos_scatt(thetaj) = cos(angle)
     enddo
  endif


  !$omp parallel &
  !$omp default(none) &
  !$omp shared(tab_s11_pos,tab_s12_pos,tab_s33_pos,tab_s34_pos,lcompute_obs) &
  !$omp shared(tab_s11,tab_s12,tab_s33,tab_s34,lambda,n_grains_tot) &
  !$omp shared(tab_albedo_pos,prob_s11_pos,valeur_prob,amax_reel,somme) &
  !$omp private(i,j,k,pk,density,k_min,proba,q_sca_tot,q_ext_tot,norme,angle)&
  !$omp shared(zmax,kappa,kappa_abs_eg,probsizecumul,ech_prob,p_n_rad,p_nz,p_n_az,j_start,pj_start) &
  !$omp shared(q_ext,q_sca,densite_pouss,scattering_method,tab_g_pos,aniso_method,tab_g,lisotropic) &
  !$omp shared(lscatt_ray_tracing,tab_s11_ray_tracing,tab_s12_ray_tracing,tab_s33_ray_tracing,tab_s34_ray_tracing) &
  !$omp shared(tab_s12_o_s11_ray_tracing,tab_s33_o_s11_ray_tracing,tab_s34_o_s11_ray_tracing,lsepar_pola)
  !$omp do schedule(dynamic,1) 
  do i=1, p_n_rad
     bz2 : do j=pj_start,p_nz
        if (j==0) cycle bz2
        do pk=1, p_n_az
           q_sca_tot=0.0
           q_ext_tot=0.0
     
           if (scattering_method == 2) then
              if (aniso_method==1) then
                 tab_s11_pos(lambda,i,j,pk,:) = 0.
                 if (lsepar_pola) then
                    tab_s12_pos(lambda,i,j,pk,:) = 0.
                    tab_s33_pos(lambda,i,j,pk,:) = 0.
                    tab_s34_pos(lambda,i,j,pk,:) = 0.
                 endif
              endif
           else
              probsizecumul(lambda,i,j,pk,0)=0.0
           endif
           
           somme=0.0
           do  k=1,n_grains_tot
              density=densite_pouss(i,j,pk,k)
              q_sca_tot = q_sca_tot + q_sca(lambda,k)*density
              q_ext_tot = q_ext_tot + q_ext(lambda,k)*density
              tab_g_pos(lambda,i,j,pk) = tab_g_pos(lambda,i,j,pk) + q_sca(lambda,k)*density * tab_g(lambda,k)
              !somme = somme + q_sca(lambda,k)*density
         
              if (scattering_method == 2) then
                 if (aniso_method==1) then
                    ! Moyennage matrice de mueller (long) (dernier indice : angle)
                    ! TODO : indice 1 suppose que mueller (via prop_grain) a ete calcule juste avant
                    tab_s11_pos(lambda,i,j,pk,:) = tab_s11(1,k,:) * density +  tab_s11_pos(lambda,i,j,pk,:)
                    if (lsepar_pola) then
                       tab_s12_pos(lambda,i,j,pk,:) = tab_s12(1,k,:) * density +  tab_s12_pos(lambda,i,j,pk,:)
                       tab_s33_pos(lambda,i,j,pk,:) = tab_s33(1,k,:) * density +  tab_s33_pos(lambda,i,j,pk,:)
                       tab_s34_pos(lambda,i,j,pk,:) = tab_s34(1,k,:) * density +  tab_s34_pos(lambda,i,j,pk,:)
                    endif
                 endif !aniso_method 
              else
                 ! Au choix suivant que l'on considère un albedo par cellule ou par grain
                 ! albedo par cellule :
                 probsizecumul(lambda,i,j,pk,k) = probsizecumul(lambda,i,j,pk,k-1) + q_sca(lambda,k)*density
                 ! albedo par grains :
                 !              probsizecumul(lambda,i,j,k) = probsizecumul(lambda,i,j,k-1) + q_ext(lambda,k)*density
              endif !scattering_method
           enddo !k
           
           
           if (lcompute_obs.and.lscatt_ray_tracing) then
              if (scattering_method == 1) then
                 write(*,*) "ERROR: ray-tracing is incompatible with scattering method 1"
                 stop
              endif

              ! Normalisation : total energie diffusee sur [0,pi] en theta et [0,2pi] en phi = 1
              norme = 0.0_db
              do thetaj=0,nang_scatt
                 angle = real(thetaj)/real(nang_scatt)*pi
                 norme=norme + tab_s11_pos(lambda,i,j,1,thetaj) * sin(angle)
              enddo
              
              if (abs(norme) > 0.0_db) then
                 norme = 1.0_db / (norme * deux_pi)

                 tab_s11_ray_tracing(lambda,i,j,:) =  tab_s11_pos(lambda,i,j,1,:) * norme
                 if (lsepar_pola) then 
                    ! Signe moins pour corriger probleme de signe pola decouvert par Gaspard 
                    ! Le transfer est fait a l'envers (direction de propagation inversee), il faut donc changer 
                    ! le signe de la matrice de Mueller
                    ! (--> supprime le signe dans dust_ray_tracing pour corriger le bug trouve par Marshall)
                    tab_s12_ray_tracing(lambda,i,j,:) =  - tab_s12_pos(lambda,i,j,1,:) * norme
                    tab_s33_ray_tracing(lambda,i,j,:) =  - tab_s33_pos(lambda,i,j,1,:) * norme
                    tab_s34_ray_tracing(lambda,i,j,:) =  - tab_s34_pos(lambda,i,j,1,:) * norme
                 endif
              else
                 tab_s11_ray_tracing(lambda,i,j,:) =  0.0_db
                 if (lsepar_pola) then
                    tab_s12_ray_tracing(lambda,i,j,:) =  0.0_db
                    tab_s33_ray_tracing(lambda,i,j,:) =  0.0_db
                    tab_s34_ray_tracing(lambda,i,j,:) =  0.0_db
                 endif
              endif ! norme

              if (lsepar_pola) then
                 tab_s12_o_s11_ray_tracing(lambda,i,j,:) = tab_s12_ray_tracing(lambda,i,j,:) / tab_s11_ray_tracing(lambda,i,j,:)
                 tab_s33_o_s11_ray_tracing(lambda,i,j,:) = tab_s33_ray_tracing(lambda,i,j,:) / tab_s11_ray_tracing(lambda,i,j,:)
                 tab_s34_o_s11_ray_tracing(lambda,i,j,:) = tab_s34_ray_tracing(lambda,i,j,:) / tab_s11_ray_tracing(lambda,i,j,:)
              endif
                 
              if (lisotropic) tab_s11_ray_tracing(lambda,i,j,:) = q_sca_tot/q_ext_tot * 1.e-4!sum(tab_s11_ray_tracing(lambda,i,j,:)) / (nang_scatt + 1) 

           endif !lscatt_ray_tracing

           if (sum(densite_pouss(i,j,pk,:)) > tiny_real) then

              if (q_ext_tot > tiny_real) tab_albedo_pos(lambda,i,j,pk) = q_sca_tot/q_ext_tot  ! NEW TEST STRAT LAURE
              if (q_sca_tot > tiny_real) tab_g_pos(lambda,i,j,pk) = tab_g_pos(lambda,i,j,pk)/q_sca_tot ! NEW TEST STRAT LAURE
              
              if (scattering_method==2) then
                 if (aniso_method==1) then
                    ! Propriétés optiques des cellules
                    prob_s11_pos(lambda,i,j,pk,0)=0.0
                    do l=1,180
                       prob_s11_pos(lambda,i,j,pk,l)=prob_s11_pos(lambda,i,j,pk,l-1)+ &
                            tab_s11_pos(lambda,i,j,pk,l)*sin(real(l)/180.*pi)  ! *pi/(180.)
                    enddo
                    if (prob_s11_pos(lambda,i,j,pk,180) > tiny_real) then ! NEW TEST STRAT LAURE
                       do l=1,180
                          prob_s11_pos(lambda,i,j,pk,l)=prob_s11_pos(lambda,i,j,pk,l)/prob_s11_pos(lambda,i,j,pk,180)
                       enddo
                    endif
                    ! Normalisation (idem que dans mueller2)
                    do l=0,180
                       if (tab_s11_pos(lambda,i,j,pk,l) > tiny_real) then ! NEW TEST STRAT LAURE
                          norme=1.0/tab_s11_pos(lambda,i,j,pk,l)
                          tab_s11_pos(lambda,i,j,pk,l)=tab_s11_pos(lambda,i,j,pk,l)*norme
                          if (lsepar_pola) then
                             tab_s12_pos(lambda,i,j,pk,l)=tab_s12_pos(lambda,i,j,pk,l)*norme
                             tab_s33_pos(lambda,i,j,pk,l)=tab_s33_pos(lambda,i,j,pk,l)*norme
                             tab_s34_pos(lambda,i,j,pk,l)=tab_s34_pos(lambda,i,j,pk,l)*norme
                          endif
                       endif
                    enddo
                 else !aniso_method
                    tab_s11_pos(lambda,i,j,pk,:)=1.0
                    if (lsepar_pola) then
                       tab_s12_pos(lambda,i,j,pk,:)=0.0
                       tab_s33_pos(lambda,i,j,pk,:)=0.0
                       tab_s34_pos(lambda,i,j,pk,:)=0.0
                    endif
                 endif !aniso_method
                 
              else !scattering_method
                 if  (probsizecumul(lambda,i,j,pk,n_grains_tot) > tiny_real) then
                    do k=1, n_grains_tot
                       probsizecumul(lambda,i,j,pk,k)= probsizecumul(lambda,i,j,pk,k)/ &
                            probsizecumul(lambda,i,j,pk,n_grains_tot)
                    enddo !k
                    ech_prob(lambda,i,j,pk,0) = 1
                    ech_prob(lambda,i,j,pk,n_prob+1)=n_grains_tot
                    k_min = 1
                    ! Cas particulier proba inf
                    ech_proba_inf : do k=1, n_grains_tot
                       if (probsizecumul(lambda,i,j,pk,k) >= proba_resol) then
                          k_min = k
                          ech_prob(lambda,i,j,pk,1) = k
                          valeur_prob(lambda,i,j,pk,1)=probsizecumul(lambda,i,j,pk,k)
                          exit ech_proba_inf
                       endif
                    enddo ech_proba_inf
                    ! cas general
                    do l=2,n_prob-1
                       proba =  real(l-1)/real(n_prob-1)
                       echantillonage : do k=k_min, n_grains_tot
                          if (probsizecumul(lambda,i,j,pk,k) >= proba) then
                             k_min = k
                             ech_prob(lambda,i,j,pk,l) = k
                             valeur_prob(lambda,i,j,pk,l)=probsizecumul(lambda,i,j,pk,k)
                             exit echantillonage
                          end if
                       end do echantillonage
                    enddo !l
                    ! Cas particulier proba sup
                    ech_proba_sup : do k=k_min, n_grains_tot
                       if (probsizecumul(lambda,i,j,pk,k) >= (1.0-proba_resol)) then
                          k_min = k
                          ech_prob(lambda,i,j,pk,n_prob) = k
                          valeur_prob(lambda,i,j,pk,n_prob)=probsizecumul(lambda,i,j,pk,k)
                          exit  ech_proba_sup
                       endif
                    enddo ech_proba_sup
                    ! Cas particulier proba=1.0
                    ech_proba1 : do k=k_min, n_grains_tot
                       if ((1.0 - probsizecumul(lambda,i,j,pk,k)) <  1.e-6) then
                          amax_reel(lambda,i,j,pk) = k
                          ech_prob(lambda,i,j,pk,n_prob+1) = k
                          valeur_prob(lambda,i,j,pk,n_prob+1) = 1.0
                          exit  ech_proba1
                       endif
                    enddo  ech_proba1 !k
                 else
                    ! a la surface, on peut avoir une proba de 0.0 partout
                    ! dans ce cas, on decide qu'il n'y a que les plus petits grains
                    ! rq : en pratique, la densite est trop faible pour qu'il y ait
                    ! une diffusion a cet endroit.
                    do  k=1,n_grains_tot
                       probsizecumul(lambda,i,j,pk,k) = 1.0
                    enddo !k
                    do l=0,n_prob+1
                       ech_prob(lambda,i,j,pk,l) = 1
                    enddo
                 endif
              endif !scattering_method
           else !densite_pouss = 0.0
              tab_albedo_pos(lambda,i,j,pk)=0.0
              if (scattering_method ==2) then
                 ! Propriétés optiques des cellules
                 prob_s11_pos(lambda,i,j,pk,:)=1.0
                 prob_s11_pos(lambda,i,j,pk,0)=0.0
                 ! Normalisation (idem que dans mueller2)
                 tab_s11_pos(lambda,i,j,pk,:)=1.0
                 if (lsepar_pola) then
                    tab_s12_pos(lambda,i,j,pk,:)=0.0
                    tab_s33_pos(lambda,i,j,pk,:)=0.0
                    tab_s34_pos(lambda,i,j,pk,:)=0.0
                 endif
              else !scattering_method
                 probsizecumul(lambda,i,j,pk,:)=1.0
                 probsizecumul(lambda,i,j,pk,0)=0.0
                 ech_prob(lambda,i,j,pk,:) = 1
                 valeur_prob(lambda,i,j,pk,:) = 1.0
              endif ! scattering_method
              
           endif !densite_pouss = 0.0
        enddo !pk
     enddo bz2 !j
  enddo !i
  !$omp enddo
  !$omp end parallel

  ! Supression scattering
  if (lno_scattering) then
     kappa = kappa_abs_eg
     tab_albedo_pos = 0.0_db
  endif

  ! scattering = abs
  if (lqsca_equal_qabs) then
     kappa = 2.0_db * kappa_abs_eg
     tab_albedo_pos = 0.5_db
  endif

  ! On remet la densite à zéro si besoin
  if (ldens0) densite_pouss(1,1,1,:) = 0.0_db
        
!  write(*,*) "lambda", tab_lambda(lambda), tab_g_pos(lambda,1,1,1)
!  write(*,*) tab_s11_ray_tracing(lambda,1,1,:)

!!$! Spline
!!$  do i=1,n_rad
!!$     do j=1,nz
!!$        do l=1,n_prob
!!$           xa(l)=valeur_prob(lambda,i,j,l) ! Proba
!!$           ya(l)=real(ech_prob(lambda,i,j,l)) ! Taille du grain
!!$        enddo !l
!!$        k=ya(1)
!!$        delta_prob1 = probsizecumul(lambda,i,j,k+1) - probsizecumul(lambda,i,j,k)
!!$        k=ya(n_prob)
!!$        delta_probn = probsizecumul(lambda,i,j,k) - probsizecumul(lambda,i,j,k-1)
!!$ !       write(*,*) probsizecumul(lambda,i,j,n_grains_tot-1), probsizecumul(lambda,i,j,n_grains_tot)
!!$ !       read(*,*)
!!$
!!$        if (delta_prob1 > 0.0) then
!!$           yp1=1.0/delta_prob1 ! 1.0 = 2-1
!!$!           yp1=(ya(2)-ya(1))/(xa(2)-xa(1))           
!!$        else
!!$           yp1=1.e30
!!$        endif
!!$        if (delta_probn > 0.0) then
!!$           ypn=1.0/delta_probn ! 1.0 = n_grains_tot -(n_grains_tot -1)
!!$!           ypn=(ya(n_prob)-ya(n_prob-1))/(xa(n_prob)-xa(n_prob-1))
!!$        else
!!$           ypn=1.e30
!!$        endif
!!$
!!$
!!$!        read(*,*)
!!$!        write(*,*) i,j,yp1,ypn
!!$!        read(*,*)
!!$
!!$!        write(*,*) 'SPLINE'
!!$        call spline(xa,ya,yp1,ypn,y2)

!!!!!!!!!!!!!!!
!!$        ! Calcul derivees a la main !CA BUGGG
!!$        do l=1,n_prob
!!$        ! derivee 1 y'=dk/dproba avec dk=1 k indice taille du grain
!!$           k=ya(l)
!!$!           yp_p1 = 1.0/(probsizecumul(lambda,i,j,k+1)-probsizecumul(lambda,i,j,k))
!!$!           yp_m1 = 1.0/(probsizecumul(lambda,i,j,k)-probsizecumul(lambda,i,j,k-1))
!!$
!!$        ! derivee 2
!!$           ypp = 1.0/(probsizecumul(lambda,i,j,k+2)-probsizecumul(lambda,i,j,k-2))
!!$!           write(*,*) l, xa(l), ya(l), y2(l), ypp
!!$           y2(l)=ypp
!!$        enddo !l
!        read(*,*)

!!!!!!!!!!!!!!!

!!$!        write(*,*) 'SPLINE'
!!$
!!$        do l=1,n_prob
!!$           xspline(lambda,i,j,l)=y2(l)
!!$        enddo !l
!!$     enddo !j
!!$  enddo !i

!  write(*,*) "test"
!  write(*,*) tab_s11_pos(1,1,143), tab_s12_pos(1,1,143), tab_s33_pos(1,1,143), tab_s34_pos(1,1,143)

! write(*,*) tab_lambda(lambda), kappa(lambda, 1, 40,1) !, kappa_abs_eg(lambda, 1, 1,1), densite_pouss(1,1,1,1)

!  write(*,*) sum(densite_pouss),   q_abs(lambda,1), tab_a(1)
!  stop

  if ((ldust_prop).and.(lambda == n_lambda)) then 
     ! Only do it after the last pass through the wavelength table
     ! in order to populate the tab_s11_pos and tab_s12_pos tables first!
     allocate(kappa_lambda(n_lambda))
     allocate(albedo_lambda(n_lambda))
     allocate(g_lambda(n_lambda))
     allocate(pol_lambda_theta(n_lambda,nang_scatt))

     kappa_lambda=real((kappa(:,1,1,1)/AU_to_cm)/(masse(1,1,1)/(volume(1)*AU_to_cm**3)))
     albedo_lambda=tab_albedo_pos(:,1,1,1)
     g_lambda=tab_g_pos(:,1,1,1)
         
     call cfitsWrite("data_dust/lambda.fits.gz",tab_lambda,shape(tab_lambda))
     call cfitsWrite("data_dust/kappa.fits.gz",kappa_lambda,shape(kappa_lambda))
     call cfitsWrite("data_dust/albedo.fits.gz",albedo_lambda,shape(albedo_lambda))
     call cfitsWrite("data_dust/g.fits.gz",g_lambda,shape(g_lambda))
     
     if (lsepar_pola) then
        do i=1,n_lambda
           pol_lambda_theta(i,:)=-tab_s12_pos(i,1,1,1,:)/tab_s11_pos(i,1,1,1,:)
        enddo
        call cfitsWrite("data_dust/polar.fits.gz",pol_lambda_theta,shape(pol_lambda_theta))
     endif

     deallocate(kappa_lambda,albedo_lambda,g_lambda,pol_lambda_theta)
     write(*,*) "Writing dust prop."
     write(*,*) "Exiting"
     stop
  endif !ldust_prop

  return

end subroutine opacite2

!******************************************************************************

end module dust
