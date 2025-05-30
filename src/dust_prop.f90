module dust_prop

  use mcfost_env
  use parametres
  use grains
  use constantes
  use wavelengths
  use density
  use cylindrical_grid
  use utils
  use scattering
  use coated_sphere
  use read_opacity

  implicit none

  real(dp), dimension(:), allocatable :: kappa_factor
  real(dp), dimension(:,:), allocatable :: kappa !n_cells, n_lambda
  real(dp), dimension(:,:), allocatable :: kappa_abs_LTE ! n_cells, n_lambda
  real(dp), dimension(:,:), allocatable :: kappa_abs_nLTE, kappa_abs_RE ! n_cells, n_lambda
  real(dp), dimension(:,:), allocatable :: proba_abs_RE, proba_abs_RE_LTE, proba_abs_RE_LTE_p_nLTE
  real(dp), dimension(:,:,:), allocatable :: kabs_nLTE_CDF, kabs_nRE_CDF ! 0:n_grains, n_cells, n_lambda

  real(dp), dimension(:,:,:), allocatable :: ksca_CDF ! 0:n_grains, n_cells, n_lambda
  !* ksca_CDF(i) represente la probabilite cumulee en-dessous d'une
  !* certaine taille de grain. Ce tableau est utilise pour le tirage
  !* aleatoire de la taille du grain diffuseur, puisqu'elle doit prendre
  !* en compte le nombre de grains en meme temps que leur probabilite
  !* individuelle de diffuser (donnee par qsca*pi*a**2).

  contains

subroutine build_grain_size_distribution()
! Distribution en taille des grains
! Sortie de prop_grains pour cas multi-lambda
! C. Pinte 14/01/05

  implicit none

  integer :: k, pop
  real :: a, nbre_tot_grains
  real(kind=dp) :: exp_grains, sqrt_exp_grains
  real :: masse_pop
  real :: correct_fact_r

  type(dust_pop_type), pointer :: d_p

  integer :: ios, status, n_comment, n_grains, alloc_status
  real :: fbuffer

  ! **************************************************
  ! Tableaux relatifs aux grains
  ! **************************************************
  allocate(nbre_grains(n_grains_tot), r_grain(n_grains_tot),  r_grain_min(n_grains_tot), r_grain_max(n_grains_tot), &
       S_grain(n_grains_tot), M_grain(n_grains_tot), r_core(n_grains_tot), grain(n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error r_grain')
  nbre_grains = 0.0   ; r_core=0.0
  r_grain=0.0 ; r_grain_min=0.0 ; r_grain_max=0.0 ; S_grain=0.0 ; M_grain=0.0

  ! Boucle sur les populations de grains
  do pop=1, n_pop
     d_p => dust_pop(pop)

      if (lread_grain_size_distrib) then
        if (n_pop > 1) call error("you cannot provide a grain size distribution with more than 1 population")

        open(unit=1, file=grain_size_file, status='old', iostat=ios)
        if (ios/=0) call error("cannot open grain size file"//trim(grain_size_file))
        write(*,*) "Reading "//trim(grain_size_file)

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
        n_grains=1 ! On a deja lu une ligne en cherchant les commentaires
        do while(status==0)
           n_grains=n_grains+1
           read(1,*,iostat=status)
        enddo
        n_grains = n_grains - 1

        if (n_grains /= n_grains_tot) call error("the number of grains must be the same as in the parameter file.")

        ! Lecture proprement dite
        rewind(1)
        ! On passe les commentaires
        do k=1, n_comment
           read(1,*)
        enddo

        ! Lecture indices
        masse_pop = 0.0
        do k=1,n_grains
           read(1,*) a, nbre_grains(k)

           r_grain(k) = a ! micron
           S_grain(k) = pi * a**2 ! micron^2
           M_grain(k) = quatre_tiers_pi * (a*mum_to_cm)**3 * d_p%rho1g_avg ! masse en g

           ! Multiplication par a car da = a.dln(a)
           nbre_grains(k) =  nbre_grains(k) * a
           grain(k)%methode_chauffage = d_p%methode_chauffage
           grain(k)%zone = d_p%zone
           grain(k)%pop = pop
           masse_pop = masse_pop + nbre_grains(k)

           if (d_p%is_PAH) grain(k)%is_PAH = .true.
        enddo ! k

        d_p%avg_grain_mass = sum(M_grain(:) * nbre_grains(:)) / sum(nbre_grains(:))

     else ! lread_grain_size_distribution

        if (d_p%aexp < 0) then
           write(*,*) "****************************************"
           write(*,*) "Warning: slope grains size negative !!!!"
           write(*,*) "****************************************"
        endif

        if (abs(d_p%amin - d_p%amax) < 1.0e-5 * d_p%amax) then
           a=d_p%amin
           d_p%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * a**3 * d_p%rho1g_avg
        else
           if (abs(d_p%aexp - 4.) > 1.0e-5) then
              if (abs(d_p%aexp - 1.) > 1.0e-5) then
                 d_p%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * d_p%rho1g_avg * &
                      (1-d_p%aexp)/(4-d_p%aexp)*(d_p%amax**(4-d_p%aexp)-d_p%amin**(4-d_p%aexp)) / &
                      (d_p%amax**(1-d_p%aexp)-d_p%amin**(1-d_p%aexp))
              else
                 d_p%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * d_p%rho1g_avg /(4-d_p%aexp) * &
                      (d_p%amax**(4-d_p%aexp)-d_p%amin**(4-d_p%aexp)) / &
                      (log(d_p%amax)-log(d_p%amin))
              endif
           else
              d_p%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * d_p%rho1g_avg *&
                   (1-d_p%aexp)*(log(d_p%amax)-log(d_p%amin)) / &
                   (d_p%amax**(1-d_p%aexp)-d_p%amin**(1-d_p%aexp))
           endif
        endif

        ! Proprietes des grains
        !exp_grains = (amax/amin)**(1./real(n_grains_tot))
        if ((d_p%n_grains==1).and.(abs(d_p%amax-d_p%amin) > 1.0e-3 * d_p%amin)) then
           write(*,*) "You have specified 1 grain size but amin != amax. Are you sure ?"
           write(*,*) "If yes, press return"
           read(*,*)
        endif

        ! Taille des grains (recursif)
        masse_pop = nbre_grains(d_p%ind_debut)
        exp_grains =  exp((1.0_dp/real(d_p%n_grains,kind=dp)) * log(d_p%amax/d_p%amin))
        sqrt_exp_grains = sqrt(exp_grains)
        do  k=d_p%ind_debut, d_p%ind_fin
           if (k==d_p%ind_debut) then
              a = d_p%amin*sqrt_exp_grains
           else
              a= r_grain(k-1) * exp_grains
           endif

           r_grain(k) = a ! micron
           r_grain_min(k) = a/sqrt_exp_grains ! taille min
           r_grain_max(k) = a*sqrt_exp_grains ! taille max

           S_grain(k) = pi * a**2 ! micron^2
           M_grain(k) = quatre_tiers_pi * (a*mum_to_cm)**3 * d_p%rho1g_avg ! masse en g

           ! Multiplication par a car da = a.dln(a)
           nbre_grains(k) = a**(-d_p%aexp) * a
           grain(k)%methode_chauffage = d_p%methode_chauffage
           grain(k)%zone = d_p%zone
           grain(k)%pop = pop
           masse_pop = masse_pop + nbre_grains(k)

           if (d_p%is_PAH) grain(k)%is_PAH = .true.
        enddo !k

     endif ! lread_grain_size_distribution

     masse_pop = masse_pop * d_p%avg_grain_mass


     ! Normalisation du nombre de grains pour atteindre la bonne masse
     nbre_grains(d_p%ind_debut:d_p%ind_fin) = nbre_grains(d_p%ind_debut:d_p%ind_fin) * d_p%masse/masse_pop

     ! Normalisation de tous les grains au sein d'une pop
     nbre_tot_grains = 0.0
     do k=d_p%ind_debut,d_p%ind_fin
        nbre_tot_grains =  nbre_tot_grains + nbre_grains(k)
     enddo

     ! Fraction de grains de taille k au sein d'une pop
     do  k=d_p%ind_debut,d_p%ind_fin
        nbre_grains(k) = nbre_grains(k)/nbre_tot_grains
     enddo !k

     ! Total radius is kept constant with coating
     if (d_p%lcoating) then
        correct_fact_r =  (1-d_p%component_volume_fraction(2))**(1./3)  ! r_core = r_tot * correct_fact_r
        do k=d_p%ind_debut,d_p%ind_fin
           r_core(k) = r_grain(k)  * correct_fact_r
        enddo ! k
     endif
  enddo !i

  return

end subroutine build_grain_size_distribution

!********************************************************************

subroutine init_indices_optiques()
  ! calcule les tableaux tab_amu1 et tab_amu2
  ! aux longueurs d'onde de tab_lambda
  ! C. Pinte
  ! 18/10/07 (separation depuis init_lambda)

  implicit none

  complex :: mavg
  complex, dimension(:), allocatable :: m

  real :: frac, fbuffer
  real, dimension(:), allocatable :: f
  integer :: i, k, ii, alloc_status, pop, status, n_ind, ios, n_comment, n_components

  character(len=512) :: filename, dir

  real, dimension(:), allocatable :: tab_l, tab_n, tab_k
  real, dimension(:,:), allocatable :: tab_tmp_amu1, tab_tmp_amu2

  logical :: l_ordre_decroissant, lread_opacity_file

  logical, save :: lfirst = .true.

  if (lfirst) then
     lread_opacity_file = .false.
  endif

  do pop=1, n_pop
     if (.not.dust_pop(pop)%is_opacity_file) then

        n_components = dust_pop(pop)%n_components
        if (dust_pop(pop)%porosity > tiny_real) then
           n_components = n_components + 1
           if (dust_pop(pop)%mixing_rule /= 1) then
              dust_pop(pop)%mixing_rule = 1 ! forcing EMT when there is only 1 component
              if (n_components > 2) then ! mcfost does not know what to do when there is more than 1 component
                 write(*,*) "dust population #", pop
                 call error("cannot use porosity with coating")
              endif
           endif
        endif
        allocate(tab_tmp_amu1(n_lambda,n_components), tab_tmp_amu2(n_lambda,n_components), stat=alloc_status)
        if (alloc_status > 0) call error('Allocation error tab_tmp_amu1')
        tab_tmp_amu1 = 0. ; tab_tmp_amu2 = 0.

        ! Lecture fichier indices
        do k=1, dust_pop(pop)%n_components
           filename = trim(dust_pop(pop)%indices(k))

           dir = in_dir(filename, dust_dir,  status=ios)
           if (ios /=0) call error("dust file cannot be found: "//trim(filename))

           filename = trim(dir)//trim(filename) ;
           write(*,*) "Reading "//trim(filename) ;

           open(unit=1,file=filename, status='old', iostat=ios)
           if (ios /=0) call error("dust file cannot be opened: "//trim(filename))

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
           allocate(tab_l(n_ind), tab_n(n_ind), tab_k(n_ind), stat=alloc_status)
           if (alloc_status > 0) call error('Allocation error zsup')
           tab_l=0.0 ; tab_n=0.0 ; tab_k=0.0

           ! Lecture proprement dite
           rewind(1)
           ! On passe les commentaires
           do i=1, n_comment
              read(1,*)
           enddo

           ! lecture densite
           read(1,*) dust_pop(pop)%component_rho1g(k), dust_pop(pop)%component_T_sub(k)
           if (dust_pop(pop)%component_rho1g(k) > 10.) call error("optical indices file has the wrong format")
           ! ligne vide
           read(1,*)

           ! Lecture indices
           do i=1,n_ind
              read(1,*) tab_l(n_ind-i+1), tab_n(n_ind-i+1), tab_k(n_ind-i+1)
           enddo

           if (tab_l(2) > tab_l(1)) then
              l_ordre_decroissant = .false.
           else
              l_ordre_decroissant = .true.
           endif
           close(unit=1)

           ! Interpolation des indices optiques aux longeurs d'onde de MCFOST
           do i=1, n_lambda
              if (l_ordre_decroissant) then
                 ! Pas d'extrapolation aux courtes longueurs d'onde
                 if (tab_lambda(i) < tab_l(n_ind)) then
                    tab_tmp_amu1(i,k)=tab_n(n_ind)
                    tab_tmp_amu2(i,k)=tab_k(n_ind)
                 else
                    ii=2 ! deplace ici car les lambda ne sont plus dans l'ordre pour l'emission moleculaire
                    do while((tab_lambda(i) < tab_l(ii)).and.(ii <= n_ind-1))
                       ii=ii+1
                    enddo

                    frac=(log(tab_l(ii))-log(tab_lambda(i)))/(log(tab_l(ii))-log(tab_l(ii-1)))
                    tab_tmp_amu1(i,k)=exp(log(tab_n(ii-1))*frac+log(tab_n(ii))*(1.0-frac))
                    tab_tmp_amu2(i,k)=exp(log(tab_k(ii-1))*frac+log(tab_k(ii))*(1.0-frac))
                 endif
              else ! ordre croissant
                 ! Pas d'extrapolation aux courtes longueurs d'onde
                 if (tab_lambda(i) < tab_l(1)) then
                    tab_tmp_amu1(i,k)=tab_n(1)
                    tab_tmp_amu2(i,k)=tab_k(1)
                 else
                    ii=2 ! deplace ici car les lambda ne sont plus dans l'ordre pour l'emission moleculaire
                    do while((tab_lambda(i) > tab_l(ii)).and.(ii <= n_ind-1))
                       ii=ii+1
                    enddo

                    frac=(log(tab_l(ii))-log(tab_lambda(i)))/(log(tab_l(ii))-log(tab_l(ii-1)))
                    tab_tmp_amu1(i,k)=exp(log(tab_n(ii-1))*frac+log(tab_n(ii))*(1.0-frac))
                    tab_tmp_amu2(i,k)=exp(log(tab_k(ii-1))*frac+log(tab_k(ii))*(1.0-frac))
                 endif
              endif ! croissant ou decroissant

           enddo ! lambda

           deallocate(tab_l,tab_n,tab_k)

        enddo ! k , boucle components

        ! Add vacuum as last component
        if (dust_pop(pop)%porosity > tiny_real) then
           tab_tmp_amu1(:,n_components) = 1.0_dp
           tab_tmp_amu2(:,n_components) = 0.0_dp
        endif

        if (n_components == 1) then
           tab_amu1(:,pop) = tab_tmp_amu1(:,1)
           tab_amu2(:,pop) = tab_tmp_amu2(:,1)
           dust_pop(pop)%T_sub = dust_pop(pop)%component_T_sub(1)
        else
           if (dust_pop(pop)%mixing_rule == 1) then ! Regle de melange
              allocate(m(n_components), f(n_components))
              do i=1, n_lambda
                 m = cmplx(tab_tmp_amu1(i,:),tab_tmp_amu2(i,:))
                 f(1:dust_pop(pop)%n_components) = dust_pop(pop)%component_volume_fraction(1:dust_pop(pop)%n_components) &
                      * (1.0_dp - dust_pop(pop)%porosity)
                 if (dust_pop(pop)%porosity > tiny_real)  f(n_components) = dust_pop(pop)%porosity

                 mavg = Bruggeman_EMT(i, m, f)

                 tab_amu1(i,pop) = real(mavg)
                 tab_amu2(i,pop) = aimag(mavg)
              enddo
              deallocate(m,f)
           else ! coating : 2 composants max pour coating
              !write(*,*) "Applying coating for pop.", pop
              if (n_components /= 2) then
                 write(*,*) "there is", n_components, "components in pop #", pop
                 call error("coating can only be computed with 2 components")
              endif
              tab_amu1(:,pop) = tab_tmp_amu1(:,1)
              tab_amu2(:,pop) = tab_tmp_amu2(:,1)
              tab_amu1_coating(:,pop) = tab_tmp_amu1(:,2)
              tab_amu2_coating(:,pop) = tab_tmp_amu2(:,2)
           endif

           ! On suppose que le grain est detruit des que le premier composant se vaporise
           write(*,*) "WARNING : T_sub is set to the minimum T_sub of the individual components"
           dust_pop(pop)%T_sub = minval(dust_pop(pop)%component_T_sub(1:dust_pop(pop)%n_components))
        endif ! n_components

        deallocate(tab_tmp_amu1,tab_tmp_amu2)

        ! Compute average material density
        dust_pop(pop)%rho1g_avg = 0.0
        do k=1, dust_pop(pop)%n_components
           dust_pop(pop)%rho1g_avg = dust_pop(pop)%rho1g_avg  + &
                dust_pop(pop)%component_rho1g(k) * dust_pop(pop)%component_volume_fraction(k)
        enddo
        dust_pop(pop)%rho1g_avg = dust_pop(pop)%rho1g_avg * (1.0-dust_pop(pop)%porosity)

        !write (*,*) "Material average density",pop,dust_pop(pop)%rho1g_avg

     else ! fichier d'opacite
        if (lfirst) then
           ! We allocate the arrays to save dimension the 1st time we enter this section
           if (.not.lread_opacity_file) then
              allocate(op_file_na(n_pop), op_file_n_lambda(n_pop), sh_file(n_pop), file_sh_nT(n_pop))
              op_file_na = 0 ; op_file_n_lambda = 0 ; sh_file = ""; file_sh_nT = 0
           endif
           lread_opacity_file = .true.

           call get_opacity_file_dim(pop)
        endif

        if (dust_pop(pop)%n_components > 1) call error("cannot mix dust with opacity file with other component")
     endif ! fichier d'opacite
  enddo ! pop

  if (lfirst.and.lread_opacity_file) call alloc_mem_opacity_file()
  lfirst = .false.

  return

end subroutine init_indices_optiques

!******************************************************************************

function Bruggeman_EMT(lambda,m,f) result(m_eff)
  ! Effective medium theory following Bruggeman rule
  !  - for 2 components : find 2nd degree polynomial root
  !  - more than 2 components : iterative solution
  !
  ! C. Pinte
  ! 28/09/11
  ! 27/10/11 : more than 2 components
  ! 18/09/13 : changing iteration for a much more stable scheme

  implicit none

  integer, intent(in) :: lambda
  complex, dimension(:), intent(in) :: m ! optical indices of the components
  real,    dimension(:), intent(in) :: f ! volume fractions of the components
  complex ::  m_eff

  integer, parameter :: n_iter_max = 100
  real, parameter :: accuracy = 1.e-6

  complex, dimension(size(m)) :: eps

  integer :: n, iter, k
  complex :: eps_eff, eps_eff1, eps_eff2, a, b, c, delta, S

  ! Number of components
  n = size(m) ;

  ! Sanity check
  if (n /= size(f)) then
     write(*,*) "Bruggeman rule error"
     write(*,*) "Incorrect number of components"
     write(*,*) "Exiting"
  endif

  ! Constantes dielectriques
  eps = m*m

  if (n==2) then ! equation du degre 2
     b = 2*f(1)*eps(1) - f(1)*eps(2) + 2*f(2)*eps(2) -f(2)*eps(1)
     c = eps(1) * eps(2)
     a = -2.   ! --> 1/2a  = -0.25
     delta = b*b - 4*a*c
     eps_eff1 = (-b - sqrt(delta)) * (-0.25)
     eps_eff2 = (-b + sqrt(delta)) * (-0.25)

     ! Selection (et verif) de l'unique racine positive (a priori c'est eps_eff1)
     if (aimag(eps_eff1) > 0) then
        if (aimag(eps_eff2) > 0) then
           call error("Bruggeman EMT rule: all imaginary parts are > 0")
        endif
        eps_eff = eps_eff1
     else
        if (aimag(eps_eff2) < 0) then
           call error("Bruggeman EMT rule: all imaginary parts are < 0")
        endif
        eps_eff = eps_eff2
     endif

  else ! n > 2, pas de solution analytique, methode iterative
     eps_eff = (sum(f(:) * eps(:)**(1./3)))**3 ! Initial guess : Landau, Lifshitz, Looyenga rule
     ! Marche aussi avec eps_eff = 1.0 mais necessite quelques iterations de plus

     eps_eff1 = eps_eff ; k=0 ! in case the iteration does not work
     iteration : do iter=1, n_iter_max
        k=k+1
        S = sum( (eps(:)-eps_eff)/(eps(:)+2*eps_eff) * f(:) )
        eps_eff =  eps_eff * (2*S+1)/(1-S)

        if (abs(S) < accuracy) exit iteration
     enddo iteration

     if (k==n_iter_max) then
        write(*,*) "****************************************************"
        write(*,*) "WARNING: lambda=",lambda, "wl=",tab_lambda(lambda)
        write(*,*) "Bruggeman not converging"
        write(*,*) "Using Landau, Lifshitz, Looyenga instead"
        write(*,*) "****************************************************"
        eps_eff = eps_eff1
     endif

  endif ! n > 2

  ! Indices optiques
  m_eff = sqrt(eps_eff)

  ! Verification finale ---> OK
  !S = sum( (eps(:)-eps_eff)/(eps(:)+2*eps_eff) * f(:) )
  !write(*,*) "Verif Bruggeman", lambda, iter, m_eff, abs(S)

  return

end function Bruggeman_EMT

!******************************************************************************

subroutine prop_grains(lambda)
  ! Calcule les propriétés des grains
  ! Masse, fraction en nombre, sections efficaces, matrices de muellers
  ! sortie des routines opacite
  ! C. Pinte 7/01/05
  ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
  ! 20/04/2023

  !use benchmarks, only : read_Pascucci_cross_sections

  implicit none

  integer, intent(in) :: lambda
  real :: a, wavel, x, qext, qsca, gsca, amu1, amu2, amu1_coat, amu2_coat
  integer :: k, pop

  qext=0.0
  qsca=0.0

  ! Longueur d'onde
  wavel=tab_lambda(lambda)

  ! Prop optiques
  ! Une premiere boucle pour les grains definis par un fichier d'indice
  !$omp parallel &
  !$omp default(none) &
  !$omp private(k,a,x,qext,qsca,gsca,amu1,amu2,pop) &
  !$omp shared(r_grain,C_ext,C_sca,C_abs,C_abs_norm,wavel,aexp,tab_albedo,lambda,tab_g,grain) &
  !$omp shared(laggregate, lFresnel, lFresnel_per_size, tab_amu1,tab_amu2,n_grains_tot,S_grain) &
  !$omp shared(tab_amu1_coating,tab_amu2_coating,amu1_coat,amu2_coat) &
  !$omp shared(dust_pop)
  !$omp do schedule(dynamic,1)
  ! on fait la boucle a l'envers pour optimiser la parallelisation
  ! et savoir des le debut si l'alloc. mem. ds bhmie passe ou pas
  do  k=n_grains_tot,1,-1
     pop = grain(k)%pop
     if (.not.dust_pop(pop)%is_opacity_file) then
        a = r_grain(k)
        amu1=tab_amu1(lambda,pop)
        amu2=tab_amu2(lambda,pop)
        x = 2.0 * pi * a / wavel
        if (laggregate) then
           call Mueller_GMM(lambda,k,qext,qsca,gsca)
        else if (lFresnel) then
           ! on désactive la parallélisation pour lire le fichier
           !$omp critical (read)
           if (lFresnel_per_size) then
	      call Fresnel_input_size(lambda,k,qext,qsca,gsca)
           else
	      call Fresnel_input(lambda,k,qext,qsca,gsca)
	   endif
	   !$omp end critical (read)
        else
           if ((dust_pop(pop)%type=="Mie").or.(dust_pop(pop)%type=="mie").or.(dust_pop(pop)%type=="MIE")) then
              if (dust_pop(pop)%lcoating) then
                 amu1_coat=tab_amu1_coating(lambda,pop)
                 amu2_coat=tab_amu2_coating(lambda,pop)
                 call Mueller_coated_sphere(lambda,k,wavel,amu1,amu2,amu1_coat,amu2_coat, qext,qsca,gsca)
              else
                 call Mueller_Mie(lambda,k,x,amu1,amu2, qext,qsca,gsca)
              endif
           else if ((dust_pop(pop)%type=="DHS").or.(dust_pop(pop)%type=="dhs")) then
              if (x < 1e4) then
                 call Mueller_DHS(lambda,k,wavel,amu1,amu2, qext,qsca,gsca)
              else
                 call Mueller_Mie(lambda,k,x,amu1,amu2, qext,qsca,gsca)
              endif
           else
              write(*,*) "Unknow dust type : ", dust_pop(pop)%type
              call error("")
           endif
        endif ! laggregate
        tab_albedo(k,lambda)=qsca/qext
        tab_g(k,lambda) = gsca

        ! section efficace
        C_ext(k,lambda) = qext * S_grain(k) ! micron^2
        C_sca(k,lambda) = qsca * S_grain(k)
        C_abs(k,lambda) = C_ext(k,lambda) - C_sca(k,lambda)

        ! Normalisation des opacites pour etre en AU^-1 pour fichier thermal_emission.f90
        ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
        ! a² microns² -> 1e-8 cm²             \
        ! density en cm-3                      > reste facteur AU_to_cm * mum_to_cm**2 = 149595.0
        ! longueur de vol en AU = 1.5e13 cm   /
        C_abs_norm(k,lambda) = C_abs(k,lambda) * AU_to_cm * mum_to_cm**2
     endif ! is_opacity_file
  enddo !k
  !$omp enddo
  !$omp end parallel

  ! Opacity files
  ! On refait exactement la meme boucle en sequentielle (pour eviter pb d'allocation en parallel)
  ! pour les grains definis par un fichier d'opacite
  ! Boucle a l'endroit cette fois
  do  k=1,n_grains_tot
     pop = grain(k)%pop
     if (dust_pop(pop)%is_opacity_file) then
        a = r_grain(k)
        x = 2.0 * pi * a / wavel
        call Mueller_opacity_file(lambda,k, qext,qsca,gsca)

        tab_albedo(k,lambda)=qsca/qext
        tab_g(k,lambda) = gsca

        C_ext(k,lambda) = qext * S_grain(k)
        C_sca(k,lambda) = qsca * S_grain(k)
        C_abs(k,lambda) = C_ext(k,lambda) - C_sca(k,lambda)
        C_abs_norm(k,lambda) = C_abs(k,lambda) * AU_to_cm * mum_to_cm**2
     endif ! is_opacity_file
  enddo !k

  if (loverwrite_s12) call overwrite_s12(Pmax)

  ! Verif normalization
  !-- do  k=1,n_grains_tot
  !--    norme = 0.0
  !--    dtheta = pi/real(nang_scatt)
  !--    do j=2,nang_scatt ! probabilite de diffusion jusqu'a l'angle j, on saute j=0 car sin(theta) = 0
  !--       theta = real(j)*dtheta
  !--       norme =  norme +   tab_s11(j,k,lambda)*sin(theta)*dtheta
  !--    enddo
  !--    qsca= C_sca(lambda,k)/S_grain(k)
  !--    write(*,*) "Verif phase function (a<<lambda) : ", tab_lambda(lambda), r_grain(k), norme/qsca, qsca
  !-- enddo

  return

end subroutine prop_grains

!******************************************************************************

subroutine save_dust_prop(letape_th)
 ! Ajout du cas ou les matrices de Mueller sont donnees en entrees
 ! 20/04/2023

  logical, intent(in) :: letape_th
  character(len=512) :: filename

  if (letape_th) then
     filename=trim(tmp_dir)//"/_dust_prop_th.tmp" ;
  else
     filename=trim(tmp_dir)//"/_dust_prop_SED.tmp" ;
  endif

  open(1,file=filename,status='replace',form='unformatted')
  write(1) para_version, scattering_method, dust_pop, grain, &
       n_lambda, lambda_min, lambda_max, tab_wavelength, &
       C_ext, C_sca, C_abs, C_abs_norm, tab_g, tab_albedo, &
       prob_s11, tab_s11, tab_s12, tab_s22, tab_s33, tab_s34, tab_s44
  close(unit=1)

  return

end subroutine save_dust_prop

!******************************************************************************

subroutine read_saved_dust_prop(letape_th, lcompute)
! Ajout du cas ou les matrices de Mueller sont donnees en entrees
! 20/04/2023

  logical, intent(in) :: letape_th
  logical, intent(out) :: lcompute

  type(dust_pop_type), dimension(n_pop) :: dust_pop_save
  type(dust_grain_type), dimension(n_grains_tot) :: grain_save
  character(len=512) :: filename, tab_wavelength_save
  integer :: n_lambda_save, scattering_method_save
  real :: lambda_min_save, lambda_max_save, para_version_save

  integer :: i, pop, ios
  logical :: ok

  lcompute = .true.
  if (lread_Misselt.or.lread_DustEM) return

  if (letape_th) then
     filename=trim(tmp_dir)//"/_dust_prop_th.tmp" ;
  else
     filename=trim(tmp_dir)//"_dust_prop_SED.tmp" ;
  endif

  ! check if there is a dust population file
  ios = 0
  open(1,file=filename,status='old',form='unformatted',iostat=ios)
  if (ios /= 0)  then
     close(unit=1)
     return
  endif

  ok = .true.
  ! read the saved dust properties
  read(1,iostat=ios) para_version_save, scattering_method_save, dust_pop_save, grain_save, &
       n_lambda_save, lambda_min_save, lambda_max_save, tab_wavelength_save, &
       C_ext, C_sca, C_abs, C_abs_norm, tab_g, tab_albedo, &
       prob_s11, tab_s11, tab_s12, tab_s22, tab_s33, tab_s34, tab_s44
  close(unit=1)
  if (ios /= 0) return ! if some dimension changed

  ! The parameter version has to be the same for safety
  ! in case something important was modified between 2 versions
  if (para_version /= para_version_save) return

  ! If the scattering method changed, the normalization is off
  if ( scattering_method /=  scattering_method_save) return

  ! check if the dust population has changed
  ok = (n_lambda == n_lambda_save) .and. (real_equality(lambda_min,lambda_min_save)) .and. &
       (real_equality(lambda_max,lambda_max_save)) .and. (tab_wavelength == tab_wavelength_save)
  if (.not.ok) return

  do pop=1, n_pop
     ok = ok .and. (dust_pop(pop)%type  == dust_pop_save(pop)%type)
     ok = ok .and. (dust_pop(pop)%n_grains  == dust_pop_save(pop)%n_grains)
     ok = ok .and. (dust_pop(pop)%methode_chauffage  == dust_pop_save(pop)%methode_chauffage)
     ok = ok .and. (dust_pop(pop)%mixing_rule  == dust_pop_save(pop)%mixing_rule)
     ok = ok .and. (dust_pop(pop)%n_components  == dust_pop_save(pop)%n_components)
     do i=1, dust_pop(pop)%n_components
        ok = ok .and. (dust_pop(pop)%indices(i)  == dust_pop_save(pop)%indices(i))
     enddo ! i
     ok = ok .and. (real_equality(dust_pop(pop)%amin,dust_pop_save(pop)%amin))
     ok = ok .and. (real_equality(dust_pop(pop)%amax,dust_pop_save(pop)%amax))
     ok = ok .and. (real_equality(dust_pop(pop)%aexp,dust_pop_save(pop)%aexp))
     ok = ok .and. (real_equality(dust_pop(pop)%frac_mass,dust_pop_save(pop)%frac_mass))
     ok = ok .and. (real_equality(dust_pop(pop)%porosity,dust_pop_save(pop)%porosity))
     ok = ok .and. (real_equality(dust_pop(pop)%dhs_maxf,dust_pop_save(pop)%dhs_maxf))
     if (.not.ok) return
  enddo ! pop

  ! for some reason workd better with a temporary structure
  grain = grain_save

  lcompute = .false.
  return

end subroutine read_saved_dust_prop

!******************************************************************************

subroutine opacite(lambda, p_lambda, no_scatt)
! Calcule la table d'opacite et ksca_CDF
! Inclus stratification empirique
! Utilise les resultats des routine densite et prop_grains
! Doit etre utilise juste apres prop_grain : lambda ne doit pas changer entre 2
! (dans le cas multi-longueurs d'onde)
! update : 12/09/06
! Ajout du cas ou les matrices de Mueller sont donnees en entrees
! 20/04/2023

  implicit none

  integer, intent(in) :: lambda, p_lambda
  logical, intent(in), optional :: no_scatt

  integer :: icell, k
  real(kind=dp) ::  density, fact, k_abs_RE, k_abs_LTE, k_abs_tot, k_sca_tot, rho0
  logical :: lcompute_obs,  ldens0, compute_scatt


  if (present(no_scatt)) then
     compute_scatt = .not.no_scatt
  else
     compute_scatt = .true.
  endif

  ! Attention : dans le cas no_strat, il ne faut pas que la cellule (1,1,1) soit vide.
  ! on la met à nbre_grains et on effacera apres
  ! c'est pour les prop de diffusion en relatif donc la veleur exacte n'a pas d'importante
  ldens0 = .false.
  if (.not.lvariable_dust) then
     if (icell_not_empty > icell1) then ! icell1==1
        ldens0 = .true.
        densite_pouss(:,icell1) = densite_pouss(:,icell_not_empty)
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

  ! Calcul opacite et probabilite de diffusion
  do icell=1, p_n_cells
     kappa(icell,lambda) = 0.0
     k_sca_tot = 0.0
     !k_abs_tot = 0.0 ! we do not need the normalisation here --> next loop now
     !k_abs_RE = 0.0

     do  k=1,n_grains_tot ! Expensive when n_cells is large
        density=densite_pouss(k,icell)
        kappa(icell,lambda) = kappa(icell,lambda) + C_ext(k,lambda) * density

        k_sca_tot = k_sca_tot + C_sca(k,lambda) * density
        !k_abs_tot = k_abs_tot + C_abs(k,lambda) * density
     enddo !k

     if (kappa(icell,lambda) > tiny_real) tab_albedo_pos(icell,lambda) = k_sca_tot/kappa(icell,lambda)

     if (aniso_method==2) then
        tab_g_pos(icell,lambda) = 0.0
        do  k=1,n_grains_tot ! Expensive when n_cells is large
           density=densite_pouss(k,icell)
           tab_g_pos(icell,lambda) = tab_g_pos(icell,lambda) + C_sca(k,lambda) * density * tab_g(k,lambda)
        enddo ! k
        if (k_sca_tot > tiny_real) tab_g_pos(icell,lambda) = tab_g_pos(icell,lambda)/k_sca_tot
     endif

     if (lRE_LTE) then
        kappa_abs_LTE(icell,lambda) = 0.0
        do k=grain_RE_LTE_start,grain_RE_LTE_end   ! Expensive when n_cells is large
           kappa_abs_LTE(icell,lambda) =  kappa_abs_LTE(icell,lambda) + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo
        !k_abs_RE = k_abs_RE + kappa_abs_LTE(icell,lambda)
     endif

     if (lRE_nLTE) then
        kappa_abs_nLTE(icell,lambda) = 0.0
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           kappa_abs_nLTE(icell,lambda) =  kappa_abs_nLTE(icell,lambda) + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo
        !k_abs_RE = k_abs_RE + kappa_abs_nLTE(icell,lambda)
     endif

     ! This has been moved to next loop :
     ! nRE opacities are updated live and per cell (as grains are flagged in equilibrium), so we can not use a cell pointer here
   !  if (letape_th) then
   !     if (lnRE.and.(k_abs_tot > tiny_dp)) then
   !        kappa_abs_RE(icell,lambda) = k_abs_RE
   !        proba_abs_RE(icell,lambda) = k_abs_RE/k_abs_tot
   !     endif
   !
   !     if (.not. (lonly_LTE.or.lonly_nLTE)) then
   !        if (k_abs_RE > tiny_dp) then
   !           Proba_abs_RE_LTE(icell,lambda) = kappa_abs_LTE(icell,lambda) / (k_abs_RE)
   !        else ! the cell is probably empty
   !           Proba_abs_RE_LTE(icell,lambda) = 1.0
   !        endif
   !     endif
   !     if (lRE_nLTE) Proba_abs_RE_LTE_p_nLTE(icell,lambda) = 1.0 ! so far, might be updated if nRE --> qRE grains
   !  endif ! letape_th

  enddo ! p_icell

  ! nRE opacities and probabilities are updated live and per cell
  ! (as grains are flagged in quasi-equilibrium), so we can not use a cell pointer here
  if ((letape_th) .and. (.not. lonly_LTE)) then
     do icell = 1, n_cells
        ! Calculating normalising opacities
        k_abs_tot = 0.0
        k_abs_RE = 0.0
        k_abs_LTE = 0.0
        do k=1, n_grains_tot
           k_abs_tot = k_abs_tot + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo
        do k=grain_RE_LTE_start,grain_RE_LTE_end
           k_abs_LTE =  k_abs_LTE + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo
        k_abs_RE = k_abs_LTE
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           k_abs_RE =  k_abs_RE + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo

        ! Computing probabilities
        if (lnRE.and.(k_abs_tot > tiny_dp)) then
           kappa_abs_RE(icell,lambda) = k_abs_RE
           proba_abs_RE(icell,lambda) = k_abs_RE / k_abs_tot
        endif

        if (.not. (lonly_LTE.or.lonly_nLTE)) then
           if (k_abs_RE > tiny_dp) then
              Proba_abs_RE_LTE(icell,lambda) = k_abs_LTE / k_abs_RE
           else ! the cell is probably empty
              Proba_abs_RE_LTE(icell,lambda) = 1.0
           endif
        endif
        if (lRE_nLTE) Proba_abs_RE_LTE_p_nLTE(icell,lambda) = 1.0 ! so far, might be updated if nRE --> qRE grains
     enddo
  endif ! letape_th and not onlyLTE


  ! proba absorption sur une taille donnée
  if (lRE_nLTE .and. (.not.low_mem_th_emission_nLTE)) then
     if (letape_th) then
        do icell=1, n_cells
           kabs_nLTE_CDF(grain_RE_nLTE_start-1,icell,lambda)=0.0
           do  k=grain_RE_nLTE_start, grain_RE_nLTE_end
              density=densite_pouss(k,icell)
              kabs_nLTE_CDF(k,icell,lambda) = kabs_nLTE_CDF(k-1,icell,lambda) + &
                   C_abs(k,lambda) * density
           enddo !k
           if (kabs_nLTE_CDF(grain_RE_nLTE_end,icell,lambda) > tiny_real) then
              kabs_nLTE_CDF(:,icell,lambda) =  kabs_nLTE_CDF(:,icell,lambda)/&
                   kabs_nLTE_CDF(grain_RE_nLTE_end,icell,lambda)
           endif
        enddo !icell
     endif !.not.lmono
  endif !lnLTE

  if (icell_not_empty <= 0) call error("could not find a non empty cell")
  rho0 = masse(icell_not_empty)/volume(icell_not_empty) ! normalising by density in a non-empty cell
  if (rho0 < tiny_dp) call error("cannot normalise by density in first non-empty cell")

  ! We apply a corrective factor per cell --> to get kappa, we need to do kappa(icell,lambda) * kappa_factor(icell)
  if (lvariable_dust) then
     kappa_factor(:) = 1.0_dp
  else
     kappa_factor(1:n_cells) = masse(1:n_cells)/volume(1:n_cells) / rho0 ! ie rho / rho(icell_not_empty)
  endif

  ! Normalisation des opacites kappa_abs pour etre en AU^-1
  ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
  ! a² microns² -> 1e-8 cm²             \
  ! density en cm-3                      > reste facteur AU_to_cm * mum_to_cm**2 = 149595.0
  ! longueur de vol en AU = 1.5e13 cm   /
  ! fact =  pi * a * a * 149595.0
  ! les k_abs_XXX n'ont pas besoin d'etre normalise car tout est relatif
  fact = AU_to_cm * mum_to_cm**2

  kappa(:,lambda) = kappa(:,lambda) * fact ! this is kappa in cell # icell_not_empty or in all cells if lvariable_dust
  if (lRE_LTE) kappa_abs_LTE(:,lambda) = kappa_abs_LTE(:,lambda) * fact
  if (lRE_nLTE) kappa_abs_nLTE(:,lambda) = kappa_abs_nLTE(:,lambda) * fact
  if (letape_th.and.lnRE) kappa_abs_RE(:,lambda) =  kappa_abs_RE(:,lambda) * fact

  if (compute_scatt) then
     if (scattering_method==2) then
        call calc_local_scattering_matrices(lambda, p_lambda)
     else ! scattering_method ==1
        if (.not.low_mem_scattering) then ! we precompute the CDF (otherwise, we recalculate it live)
           do icell=1, p_n_cells
              ksca_CDF(0,icell,p_lambda)=0.0

              do  k=1,n_grains_tot
                 ksca_CDF(k,icell,p_lambda) = ksca_CDF(k-1,icell,p_lambda) + C_sca(k,lambda) * densite_pouss(k,icell)
              enddo !k

              if  (ksca_CDF(n_grains_tot,icell,p_lambda) > tiny_real) then
                 ksca_CDF(:,icell,p_lambda)= ksca_CDF(:,icell,p_lambda)/ ksca_CDF(n_grains_tot,icell,p_lambda)
              else
                 ! a la surface, on peut avoir une proba de 0.0 partout
                 ! dans ce cas, on decide qu'il n'y a que les plus petits grains
                 ! rq : en pratique, la densite est trop faible pour qu'il y ait
                 ! une diffusion a cet endroit.
                 ksca_CDF(:,icell,p_lambda) = 1.0
              endif
           enddo ! icell
        endif ! .not.low_mem_scattering
     endif
  endif

  ! Supression scattering
  if (lno_scattering) then
     kappa = kappa_abs_LTE
     tab_albedo_pos = 0.0_dp
  endif

  ! scattering = abs
  if (lqsca_equal_qabs) then
     kappa = 2.0_dp * kappa_abs_LTE
     tab_albedo_pos = 0.5_dp
  endif

  ! On remet la densite à zéro si besoin
  if (ldens0) then
     densite_pouss(:,icell1) = 0.0_sp
  endif

  if ((ldust_prop).and.(lambda == n_lambda)) then
     call write_dust_prop()

     if (lstop_after_init) then
        write(*,*) "Exiting"
        call exit(0)
     else
        ! Re-Normalisation S11
        ! la normalisation n'a pas eu lieu dans le cas ldust_prop pour sauver S11 dans le fichier fits
        ! on l'a fait donc maintenant pour comtinuer les calculs
        if ((scattering_method==2).and.(aniso_method==1)) tab_s11_pos = 1.0
     endif

  endif !ldust_prop

  return

end subroutine opacite

!******************************************************************************

subroutine calc_local_scattering_matrices(lambda, p_lambda)

  integer, intent(in) :: lambda, p_lambda

  real(kind=dp), parameter :: dtheta = pi/real(nang_scatt)
  real(kind=dp) :: density, theta, norme, fact, k_sca_tot
  real :: mu, g, g2

  integer :: icell, k, l
  logical :: ldens0

  fact = AU_to_cm * mum_to_cm**2
  !write(*,*) "Computing local scattering properties", lambda, p_lambda

  ! see opacite()
  ! Attention : dans le cas no_strat, il ne faut pas que la cellule (1,1,1) soit vide.
  ! on la met à nbre_grains et on effacera apres
  ! c'est pour les prop de diffusion en relatif donc la veleur exacte n'a pas d'importante
  ldens0 = .false.
  if (.not.lvariable_dust) then
     if (icell_not_empty > icell1) then
        ldens0 = .true.
        densite_pouss(:,icell1) = densite_pouss(:,icell_not_empty)
     endif
  endif

  !$omp parallel &
  !$omp default(none) &
  !$omp shared(tab_s11_pos,tab_s12_o_s11_pos,tab_s22_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos,tab_s44_o_s11_pos) &
  !$omp shared(tab_s11,tab_s12,tab_s22,tab_s33,tab_s34,tab_s44,lambda,p_lambda,n_grains_tot,tab_albedo_pos,prob_s11_pos) &
  !$omp shared(zmax,kappa,kappa_abs_LTE,ksca_CDF,p_n_cells,fact) &
  !$omp shared(C_ext,C_sca,densite_pouss,S_grain,scattering_method,tab_g_pos,aniso_method,tab_g,lisotropic,low_mem_scattering) &
  !$omp shared(lscatt_ray_tracing,letape_th,lsepar_pola,ldust_prop,lphase_function_file,s11_file,loverwrite_s12,Pmax) &
  !$omp private(icell,k,density,norme,theta,k_sca_tot,mu,g,g2)
  !$omp do schedule(dynamic,1)
  do icell=1, p_n_cells
     if (aniso_method==1) then
        tab_s11_pos(:,icell,p_lambda) = 0.
        if (lsepar_pola) then
           tab_s12_o_s11_pos(:,icell,p_lambda) = 0.
           tab_s22_o_s11_pos(:,icell,p_lambda) = 0.
           tab_s33_o_s11_pos(:,icell,p_lambda) = 0.
           tab_s34_o_s11_pos(:,icell,p_lambda) = 0.
           tab_s44_o_s11_pos(:,icell,p_lambda) = 0.
        endif
     endif

     do  k=1,n_grains_tot
        density=densite_pouss(k,icell)
        if (aniso_method==1) then
           ! Moyennage matrice de mueller (long en cpu ) (le dernier indice est l'angle)
           ! tab_s11 est normalisee a Qsca --> facteur S_grain * density pour que
           ! tab_s11_pos soit normalisee a k_sca_tot
           tab_s11_pos(:,icell,p_lambda) = tab_s11_pos(:,icell,p_lambda) + tab_s11(:,k,lambda) * S_grain(k) * density
           if (lsepar_pola) then
              tab_s12_o_s11_pos(:,icell,p_lambda) = tab_s12_o_s11_pos(:,icell,p_lambda) + tab_s12(:,k,lambda) * S_grain(k) * density
              tab_s22_o_s11_pos(:,icell,p_lambda) = tab_s22_o_s11_pos(:,icell,p_lambda) + tab_s22(:,k,lambda) * S_grain(k) * density
              tab_s33_o_s11_pos(:,icell,p_lambda) = tab_s33_o_s11_pos(:,icell,p_lambda) + tab_s33(:,k,lambda) * S_grain(k) * density
              tab_s34_o_s11_pos(:,icell,p_lambda) = tab_s34_o_s11_pos(:,icell,p_lambda) + tab_s34(:,k,lambda) * S_grain(k) * density
              tab_s44_o_s11_pos(:,icell,p_lambda) = tab_s44_o_s11_pos(:,icell,p_lambda) + tab_s44(:,k,lambda) * S_grain(k) * density
           endif
        endif !aniso_method
     enddo !k

     k_sca_tot = kappa(icell,lambda) * tab_albedo_pos(icell,lambda) / fact ! We renormalize to remove the factor from opacity

     ! Over-riding phase function
     if (lphase_function_file)  then
        tab_s11_pos(:,icell,p_lambda) = s11_file(:)

        ! Normalisation of phase-function to k_sca_tot, ie like the internal routines
        norme = 0.0
        do l=1,nang_scatt-1 ! on saute j=0 & nang_scatt car sin(theta) = 0
           theta = real(l)*dtheta
           norme = norme + tab_s11_pos(l,icell,p_lambda)*sin(theta)*dtheta
        enddo
        tab_s11_pos(:,icell,p_lambda) = tab_s11_pos(:,icell,p_lambda) * k_sca_tot / norme
     endif

     if (k_sca_tot > tiny_real) then

        if (aniso_method==1) then ! full phase function
           ! Propriétés optiques des cellules
           prob_s11_pos(0,icell,p_lambda)=0.0

           do l=2,nang_scatt ! probabilite de diffusion jusqu'a l'angle j, on saute j=0 car sin(theta) = 0
              theta = real(l)*dtheta
              prob_s11_pos(l,icell,p_lambda)=prob_s11_pos(l-1,icell,p_lambda)+ &
                   tab_s11_pos(l,icell,p_lambda)*sin(theta)*dtheta
           enddo

           ! tab_s11_pos est calculee telle que la normalisation soit: k_sca_tot
           prob_s11_pos(1:nang_scatt,icell,p_lambda) = prob_s11_pos(1:nang_scatt,icell,p_lambda) + &
                k_sca_tot - prob_s11_pos(nang_scatt,icell,p_lambda)

           ! Normalisation de la proba cumulee a 1
           prob_s11_pos(:,icell,p_lambda)=prob_s11_pos(:,icell,p_lambda)/k_sca_tot

           ! Normalisation des matrices de Mueller (idem que dans mueller_Mie)
           do l=0,nang_scatt
              if (tab_s11_pos(l,icell,p_lambda) > tiny_real) then
                 norme=1.0/tab_s11_pos(l,icell,p_lambda)
                 if (lsepar_pola) then
                    tab_s12_o_s11_pos(l,icell,p_lambda)=tab_s12_o_s11_pos(l,icell,p_lambda)*norme
                    tab_s22_o_s11_pos(l,icell,p_lambda)=tab_s22_o_s11_pos(l,icell,p_lambda)*norme
                    tab_s33_o_s11_pos(l,icell,p_lambda)=tab_s33_o_s11_pos(l,icell,p_lambda)*norme
                    tab_s34_o_s11_pos(l,icell,p_lambda)=tab_s34_o_s11_pos(l,icell,p_lambda)*norme
                    tab_s44_o_s11_pos(l,icell,p_lambda)=tab_s44_o_s11_pos(l,icell,p_lambda)*norme
                 endif
              endif
           enddo

           ! Normalisation : on veut que l'energie total diffusee sur [0,pi] en theta et [0,2pi] en phi = 1
           ! (for ray-tracing)
           tab_s11_pos(:,icell,p_lambda) = tab_s11_pos(:,icell,p_lambda) * dtheta / (k_sca_tot * deux_pi)

           if (lsepar_pola .and. loverwrite_s12) then
              do l=0,nang_scatt
                 theta = real(l)*dtheta
                 tab_s12_o_s11_pos(l,icell,p_lambda) = - Pmax * (1-(cos(theta))**2)
              enddo
           endif

           ! -- ! aniso_method = 2 --> HG
           ! -- gsca = tab_g_pos(icell,lambda)
           ! -- tab_s11_ray_tracing(:,icell,lambda) =  tab_s11_ray_tracing(:,icell,lambda) / (k_sca_tot * deux_pi)
           ! --
           ! -- if (lisotropic) tab_s11_ray_tracing(:,icell,lambda) = 1.0 / (4.* nang_scatt)
        else !aniso_method == 2 : HG

           ! Normalisation : on veut que l'energie total diffusee sur [0,pi] en theta et [0,2pi] en phi = 1
           ! (for ray-tracing)
           do l=0,nang_scatt
              g = tab_g_pos(icell,p_lambda) ; g2 = g**2
              mu = cos((real(l))/real(nang_scatt)*pi) ! cos(theta)
              tab_s11_pos(l,icell,p_lambda) = 1.0/quatre_pi * (1-g2) * (1+g2-2*g*mu)**(-1.5) * dtheta
           enddo

           if (lsepar_pola) then
              tab_s12_o_s11_pos(:,icell,p_lambda)=0.0
              tab_s22_o_s11_pos(:,icell,p_lambda)=0.0
              tab_s33_o_s11_pos(:,icell,p_lambda)=0.0
              tab_s34_o_s11_pos(:,icell,p_lambda)=0.0
              tab_s44_o_s11_pos(:,icell,p_lambda)=0.0
           endif
        endif !aniso_method

        ! todo :
        ! utiliser tab_s11_pos pour raie tracing et tab_s12_o_s11
        ! 1) normaliser les tab_s12_o_s11
        ! 2) normaliser les tab_s11 avec  norme = dtheta / (k_sca_tot * deux_pi)

        ! Multi ou mono-lambda : mono-lambda si 180 * p_n_cells * p_n_lambda >> 1
        ! 1) p_n_lambda = n_lambda au step 1
        ! 2) p_n_lambda = 1 au step 2 si p_n_cells > 1, sinon on peut utiliser n_lambda
        ! 3) on ne recalcule pas si fait au step 1

     else ! k_sca_tot = 0.0
        tab_albedo_pos(icell,lambda)=0.0
        ! Propri\E9t\E9s optiques des cellules
        prob_s11_pos(:,icell,p_lambda)=1.0
        prob_s11_pos(0,icell,p_lambda)=0.0
        ! Normalisation (idem que dans mueller)
        tab_s11_pos(:,icell,p_lambda)=1.0
        if (lsepar_pola) then
           tab_s12_o_s11_pos(:,icell,p_lambda)=0.0
           tab_s22_o_s11_pos(:,icell,p_lambda)=0.0
           tab_s33_o_s11_pos(:,icell,p_lambda)=0.0
           tab_s34_o_s11_pos(:,icell,p_lambda)=0.0
           tab_s44_o_s11_pos(:,icell,p_lambda)=0.0
        endif

     endif ! k_sca_tot > or = 0
  enddo !icell
  !$omp enddo
  !$omp end parallel

  ! see opacite()
  ! On remet la densite à zéro si besoin
  if (ldens0) then
     densite_pouss(:,icell1) = 0.0_sp
  endif

  return

end subroutine calc_local_scattering_matrices

!******************************************************************************

integer function select_grainsize_high_mem(lambda,aleat, icell) result(k)
!-----------------------------------------------------
!  Nouvelle version : janvier 04
!  la dichotomie se fait en comparant les indices
!  et non plus en utilisant la taille du pas
!  --> plus precis, important quand on reduit n_grains_tot
!
!  23 Fevrier 04 : ajout d'une prediction avant la
!  premiere iteration de la dichotomie
!-----------------------------------------------------

  implicit none

  integer, intent(in) :: lambda, icell
  real, intent(in) :: aleat
  real :: prob
  integer :: kmin, kmax

  prob = aleat ! ksca_CDF(n_grains_tot,icell,lambda) est normalise a 1.0

  ! dichotomie
  kmin = 0
  kmax = n_grains_tot
  k=(kmin + kmax)/2

  do while (ksca_CDF(k,icell,lambda) /= prob)
     if (ksca_CDF(k,icell,lambda) < prob) then
        kmin = k
     else
        kmax = k
     endif

     k = (kmin + kmax)/2
     if ((kmax-kmin) <= 1) then
        exit ! Sortie du while
     endif
  enddo   ! while
  k=kmax

  return

end function select_grainsize_high_mem

!***************************************************

integer function select_scattering_grain(lambda,icell, aleat) result(k)
  ! This routine will select randomly the scattering grain from the CDF of ksca
  ! Because we cannot store all the CDF for all cells (n_grains x ncells x n_lambda),
  ! the CDF is recomputed on the fly here.
  ! The normalization is saved via kappa_sca, so we do not need to compute the whole CDF.

  implicit none

  integer, intent(in) :: lambda, icell
  real, intent(in) :: aleat
  real :: prob, CDF, norm

  if (low_mem_scattering) then
     ! We scale the random number so that it is between 0 and kappa_sca (= last value of CDF)
     norm =  kappa(icell,lambda) * tab_albedo_pos(icell,lambda) / (AU_to_cm * mum_to_cm**2)

     if (aleat < 0.5) then ! We start from first grain
        prob = aleat * norm
        CDF = 0.0
        do k=1, n_grains_tot
           CDF = CDF + C_sca(k,lambda) * densite_pouss(k,icell)
           if (CDF > prob) exit
        enddo
     else ! We start from the end of the grain size distribution
        prob = (1.0-aleat) * norm
        CDF = 0.0
        do k=n_grains_tot, 1, -1
           CDF = CDF + C_sca(k,lambda) * densite_pouss(k,icell)
           if (CDF > prob) exit
        enddo
     endif
  else
     k = select_grainsize_high_mem(lambda,aleat, icell) ! is this ever used ?
  endif

  return

end function select_scattering_grain

!*******************************************************************

subroutine write_dust_prop()
! Ajout du cas ou les matrices de Mueller sont donnees en entrees
! 20/04/2023

  use fits_utils, only : cfitsWrite
  use density, only : masse

  integer :: icell, l, p_icell

  real, dimension(:), allocatable :: kappa_lambda,albedo_lambda,g_lambda
  real, dimension(:,:), allocatable :: S11_lambda_theta, pol_lambda_theta, kappa_grain

  write(*,*) "Writing dust properties"
  ! Rewrite step2 on top of step1 (we still get step 1 if step 2 does not finish)

  ! Only do it after the last pass through the wavelength table
  ! in order to populate the tab_s11_pos and tab_s12_pos tables first!
  allocate(kappa_lambda(n_lambda))
  allocate(albedo_lambda(n_lambda))
  allocate(g_lambda(n_lambda))
  allocate(S11_lambda_theta(n_lambda,0:nang_scatt),pol_lambda_theta(n_lambda,0:nang_scatt))
  allocate(kappa_grain(n_lambda,n_grains_tot))

  if (lvariable_dust) then
     write(*,*) "Warning: dust is not uniform, picking cell #",  icell_not_empty
     icell = icell_not_empty
     p_icell = icell_not_empty
  else
     icell  = icell_not_empty
     p_icell = 1
  endif

  kappa_lambda=real((kappa(icell,:)*kappa_factor(icell)/AU_to_cm)/(masse(icell)/(volume(icell)*AU_to_cm**3))) ! cm^2/g
  albedo_lambda=tab_albedo_pos(icell,:)

  call cfitsWrite("!data_dust/lambda.fits.gz",real(tab_lambda),shape(tab_lambda))
  call cfitsWrite("!data_dust/kappa.fits.gz",kappa_lambda,shape(kappa_lambda))
  call cfitsWrite("!data_dust/albedo.fits.gz",albedo_lambda,shape(albedo_lambda))

  if (aniso_method==2) then
     g_lambda=tab_g_pos(icell,:)
     call cfitsWrite("!data_dust/g.fits.gz",g_lambda,shape(g_lambda))
  endif

  do l=1, n_lambda
     kappa_grain(l,:) = C_abs(:,l) * mum_to_cm**2 / M_grain(:) ! cm^2/g
  enddo
  call cfitsWrite("!data_dust/kappa_grain.fits.gz",kappa_grain,shape(kappa_grain)) ! lambda, n_grains

  do l=1, n_lambda
     S11_lambda_theta(l,:)= tab_s11_pos(:,p_icell,l)
  enddo
  call cfitsWrite("!data_dust/phase_function.fits.gz",S11_lambda_theta,shape(S11_lambda_theta))

  if (lsepar_pola) then
     do l=1, n_lambda
        pol_lambda_theta(l,:) = -tab_s12_o_s11_pos(:,p_icell,l) ! Deja normalise par S11
     enddo
     call cfitsWrite("!data_dust/polarizability.fits.gz",pol_lambda_theta,shape(pol_lambda_theta))
  endif

  deallocate(kappa_lambda,albedo_lambda,g_lambda,S11_lambda_theta,pol_lambda_theta,kappa_grain)

  return

end subroutine write_dust_prop

end module dust_prop
