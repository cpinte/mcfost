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
  use coated_sphere
  use input
  use output

  implicit none

  contains

subroutine build_grain_size_distribution()
! Distribution en taille des grains
! Sortie de prop_grains pour cas multi-lambda
! C. Pinte 14/01/05

  implicit none

  integer :: k, pop
  real :: a, nbre_tot_grains
  real(kind=db) :: exp_grains, sqrt_exp_grains
  real :: masse_pop
  real :: correct_fact_r

  type(dust_pop_type), pointer :: dp

  integer :: ios, status, n_comment, n_grains, alloc_status
  real :: fbuffer

  ! **************************************************
  ! Tableaux relatifs aux grains
  ! **************************************************
  allocate(nbre_grains(n_grains_tot), r_grain(n_grains_tot),  r_grain_min(n_grains_tot), r_grain_max(n_grains_tot), &
       S_grain(n_grains_tot), M_grain(n_grains_tot), r_core(n_grains_tot), grain(n_grains_tot), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error r_grain'
     stop
  endif
  nbre_grains = 0.0   ; r_core=0.0
  r_grain=0.0 ; r_grain_min=0.0 ; r_grain_max=0.0 ; S_grain=0.0 ; M_grain=0.0

  ! Boucle sur les populations de grains
  do pop=1, n_pop
     dp => dust_pop(pop)

      if (lread_grain_size_distrib) then
        if (n_pop > 1) then
           write(*,*) "ERROR: you cannot provide a grain size distribution with more than 1 population"
           write(*,*) "Exiting"
           stop
        endif

        open(unit=1, file=grain_size_file, status='old', iostat=ios)
        if (ios/=0) then
           write(*,*) "ERROR: cannot open grain size file"//trim(grain_size_file)
           write(*,*) "Exiting"
           stop
        endif
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

        if (n_grains /= n_grains_tot) then
           write(*,*) "ERROR: the number of grains must be the same as in the parameter file"
           write(*,*) "I will correct that later", n_grains, n_grains_tot
           write(*,*) "Exiting"
           stop
        endif


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
           M_grain(k) = quatre_tiers_pi * (a*mum_to_cm)**3 * dp%rho1g_avg ! masse en g

           ! Multiplication par a car da = a.dln(a)
           nbre_grains(k) =  nbre_grains(k) * a
           grain(k)%methode_chauffage = dp%methode_chauffage
           grain(k)%zone = dp%zone
           grain(k)%pop = pop
           masse_pop = masse_pop + nbre_grains(k)

           if (dp%is_PAH) grain(k)%is_PAH = .true.
        enddo ! k

        dp%avg_grain_mass = sum(M_grain(:) * nbre_grains(:)) / sum(nbre_grains(:))

     else ! lread_grain_size_distribution

        if (dp%aexp < 0) then
           write(*,*) "****************************************"
           write(*,*) "Warning: slope grains size negative !!!!"
           write(*,*) "****************************************"
        endif

        if (abs(dp%amin - dp%amax) < 1.0e-5 * dp%amax) then
           a=dp%amin
           dp%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * a**3 * dp%rho1g_avg
        else
           if (abs(dp%aexp - 4.) > 1.0e-5) then
              if (abs(dp%aexp - 1.) > 1.0e-5) then
                 dp%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * dp%rho1g_avg * &
                      (1-dp%aexp)/(4-dp%aexp)*(dp%amax**(4-dp%aexp)-dp%amin**(4-dp%aexp)) / &
                      (dp%amax**(1-dp%aexp)-dp%amin**(1-dp%aexp))
              else
                 dp%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * dp%rho1g_avg /(4-dp%aexp) * &
                      (dp%amax**(4-dp%aexp)-dp%amin**(4-dp%aexp)) / &
                      (log(dp%amax)-log(dp%amin))
              endif
           else
              dp%avg_grain_mass = quatre_tiers_pi * mum_to_cm**3 * dp%rho1g_avg *&
                   (1-dp%aexp)*(log(dp%amax)-log(dp%amin)) / &
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

        ! Taille des grains (recursif)
        masse_pop = nbre_grains(dp%ind_debut)
        exp_grains =  exp((1.0_db/real(dp%n_grains,kind=db)) * log(dp%amax/dp%amin))
        sqrt_exp_grains = sqrt(exp_grains)
        do  k=dp%ind_debut, dp%ind_fin
           if (k==dp%ind_debut) then
              a = dp%amin*sqrt_exp_grains
           else
              a= r_grain(k-1) * exp_grains
           endif

           r_grain(k) = a ! micron
           r_grain_min(k) = a/sqrt_exp_grains ! taille min
           r_grain_max(k) = a*sqrt_exp_grains ! taille max

           S_grain(k) = pi * a**2 ! micron^2
           M_grain(k) = quatre_tiers_pi * (a*mum_to_cm)**3 * dp%rho1g_avg ! masse en g

           ! Multiplication par a car da = a.dln(a)
           nbre_grains(k) = a**(-dp%aexp) * a
           grain(k)%methode_chauffage = dp%methode_chauffage
           grain(k)%zone = dp%zone
           grain(k)%pop = pop
           masse_pop = masse_pop + nbre_grains(k)

           if (dp%is_PAH) grain(k)%is_PAH = .true.
        enddo !k

     endif ! lread_grain_size_distribution

     masse_pop = masse_pop * dp%avg_grain_mass


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

     ! Total radius is kept constant with coating
     if (dp%lcoating) then
        correct_fact_r =  (1-dp%component_volume_fraction(2))**(1./3)  ! r_core = r_tot * correct_fact_r
        do k=dp%ind_debut,dp%ind_fin
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
        if (dust_pop(pop)%porosity > tiny_real) n_components = n_components + 1
        allocate(tab_tmp_amu1(n_lambda,n_components), tab_tmp_amu2(n_lambda,n_components), stat=alloc_status)
        if (alloc_status > 0) then
           write(*,*) 'Allocation error tab_tmp_amu1'
           stop
        endif
        tab_tmp_amu1 = 0. ; tab_tmp_amu2 = 0.

        ! Lecture fichier indices
        do k=1, dust_pop(pop)%n_components
           filename = trim(dust_pop(pop)%indices(k))

           dir = in_dir(filename, dust_dir,  status=ios)
           if (ios /=0) then
              write(*,*) "ERROR: dust file cannot be found:",trim(filename)
              write(*,*) "Exiting"
              stop
           else
              filename = trim(dir)//trim(filename) ;
              write(*,*) "Reading "//trim(filename) ;
           endif

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
           allocate(tab_l(n_ind), tab_n(n_ind), tab_k(n_ind), stat=alloc_status)
           if (alloc_status > 0) then
              write(*,*) 'Allocation error zsup'
              stop
           endif
           tab_l=0.0 ; tab_n=0.0 ; tab_k=0.0

           ! Lecture proprement dite
           rewind(1)
           ! On passe les commentaires
           do i=1, n_comment
              read(1,*)
           enddo

           ! lecture densite
           read(1,*) dust_pop(pop)%component_rho1g(k), dust_pop(pop)%component_T_sub(k)
           if (dust_pop(pop)%component_rho1g(k) > 10.) then
              write(*,*) "ERROR: optical indices file has the wrong format"
              write(*,*) "Exiting"
              stop
           endif
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
                 ii=2 ! deplace ici car les lambda ne sont plus dans l'ordre pour l'emission moleculaire
                 do while((tab_lambda(i) < tab_l(ii)).and.(ii <= n_ind-1))
                    ii=ii+1
                 enddo

                 frac=(log(tab_l(ii))-log(tab_lambda(i)))/(log(tab_l(ii))-log(tab_l(ii-1)))
                 tab_tmp_amu1(i,k)=exp(log(tab_n(ii-1))*frac+log(tab_n(ii))*(1.0-frac))
                 tab_tmp_amu2(i,k)=exp(log(tab_k(ii-1))*frac+log(tab_k(ii))*(1.0-frac))

              else ! ordre croissant
                 ii=2 ! deplace ici car les lambda ne sont plus dans l'ordre pour l'emission moleculaire
                 do while((tab_lambda(i) > tab_l(ii)).and.(ii <= n_ind-1))
                    ii=ii+1
                 enddo

                 frac=(log(tab_l(ii))-log(tab_lambda(i)))/(log(tab_l(ii))-log(tab_l(ii-1)))
                 tab_tmp_amu1(i,k)=exp(log(tab_n(ii-1))*frac+log(tab_n(ii))*(1.0-frac))
                 tab_tmp_amu2(i,k)=exp(log(tab_k(ii-1))*frac+log(tab_k(ii))*(1.0-frac))
              endif ! croissant ou decroissant

           enddo ! lambda

           deallocate(tab_l,tab_n,tab_k)

        enddo ! k , boucle components

        ! Add vacuum as last component
        if (dust_pop(pop)%porosity > tiny_real) then
           tab_tmp_amu1(:,n_components) = 1.0_db
           tab_tmp_amu2(:,n_components) = 0.0_db
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
                      * (1.0_db - dust_pop(pop)%porosity)
                 if (dust_pop(pop)%porosity > tiny_real)  f(n_components) = dust_pop(pop)%porosity


                 mavg = Bruggeman_EMT(i, m, f)

                 tab_amu1(i,pop) = real(mavg)
                 tab_amu2(i,pop) = aimag(mavg)
              enddo
              deallocate(m,f)
           else ! coating : 2 composants max pour coating
              !write(*,*) "Applying coating for pop.", pop
              if (n_components /= 2) then
                 write(*,*) "ERROR: coating can only be computed with 2 components"
                 write(*,*) "there is", n_components, "components in pop #", pop
                 write(*,*) "Exiting"
                 stop
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

        if (dust_pop(pop)%n_components > 1) then
           write(*,*) "ERROR : cannot mix dust with opacity file with other component"
           write(*,*) "Exiting"
           stop
        endif
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
           write(*,*) "Error in Bruggeman EMT rule !!!"
           write(*,*) "All imaginary parts are > 0"
           write(*,*) "Exiting."
           stop
        endif
        eps_eff = eps_eff1
     else
        if (aimag(eps_eff2) < 0) then
           write(*,*) "Error in Bruggeman EMT rule !!!"
           write(*,*) "All imaginary parts are < 0"
           write(*,*) "Exiting."
           stop
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

  implicit none

  integer, intent(in) :: lambda
  real :: a, wavel, x, qext, qsca, gsca, amu1, amu2, amu1_coat, amu2_coat
  integer :: k, pop

  !-- real :: norme, dtheta, theta
  !-- integer :: j

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
  !$omp shared(laggregate,tab_amu1,tab_amu2,n_grains_tot,S_grain) &
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
              write(*,*) "Exiting"
              stop
           endif
           !           write(*,*) wavel, qext,qsca,gsca
        endif ! laggregate
        tab_albedo(k,lambda)=qsca/qext
        tab_g(k,lambda) = gsca

        ! section efficace
        C_ext(k,lambda) = qext * S_grain(k)
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
  !-- stop

  return

end subroutine prop_grains

!******************************************************************************

subroutine save_dust_prop(letape_th)

  logical, intent(in) :: letape_th
  character(len=512) :: filename

  if (letape_th) then
     filename="_dust_prop_th.tmp" ;
  else
     filename="_dust_prop_SED.tmp" ;
  endif

  open(1,file=filename,status='replace',form='unformatted')
  write(1) para_version, scattering_method, dust_pop, grain, &
       n_lambda, lambda_min, lambda_max, tab_wavelength, &
       C_ext, C_sca, C_abs, C_abs_norm, tab_g, tab_albedo, prob_s11, tab_s11, tab_s12, tab_s33, tab_s34
  close(unit=1)

  return

end subroutine save_dust_prop

!******************************************************************************

subroutine read_saved_dust_prop(letape_th, lcompute)

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
  if (lread_Misselt) return

  if (letape_th) then
     filename="_dust_prop_th.tmp" ;
  else
     filename="_dust_prop_SED.tmp" ;
  endif

  ! check if there is a dust population file
  ios = 0
  open(1,file=filename,status='old',form='unformatted',iostat=ios)
  if (ios /= 0)  then
     close(unit=1)
     return
  endif

  ! read the saved dust properties
  read(1,iostat=ios) para_version_save, scattering_method_save, dust_pop_save, grain_save, &
       n_lambda_save, lambda_min_save, lambda_max_save, tab_wavelength_save, &
       C_ext, C_sca, C_abs, C_abs_norm, tab_g, tab_albedo, prob_s11, tab_s11, tab_s12, tab_s33, tab_s34
  close(unit=1)
  if (ios /= 0) then ! if some dimension changed
     return
  endif

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

subroutine opacite(lambda, p_lambda)
! Calcule la table d'opacite et ksca_CDF
! Inclus stratification empirique
! Utilise les resultats des routine densite et prop_grains
! Doit etre utilise juste apres prop_grain : lambda ne doit pas changer entre 2
! (dans le cas multi-longueurs d'onde)
! update : 12/09/06

  implicit none

  integer, intent(in) :: lambda, p_lambda

  integer :: icell, k, thetaj
  real(kind=db) ::  density, fact, k_abs_RE, k_abs_tot, angle
  logical :: lcompute_obs,  ldens0

  ! Attention : dans le cas no_strat, il ne faut pas que la cellule (1,1,1) soit vide.
  ! on la met à nbre_grains et on effacera apres
  ! c'est pour les prop de diffusion en relatif donc la veleur exacte n'a pas d'importante
  ldens0 = .false.
  if (.not.lvariable_dust) then
     icell = icell_ref
     if (maxval(densite_pouss(:,icell)) < tiny_real) then
        ldens0 = .true.
        densite_pouss(:,icell) = nbre_grains(:)
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
  do icell=1, n_cells
     kappa(icell,lambda) = 0.0
     kappa_sca(icell,lambda) = 0.0
     k_abs_tot = 0.0
     k_abs_RE = 0.0

     do  k=1,n_grains_tot
        density=densite_pouss(k,icell)
        kappa(icell,lambda) = kappa(icell,lambda) + C_ext(k,lambda) * density
        kappa_sca(icell,lambda) = kappa_sca(icell,lambda) + C_sca(k,lambda) * density
        tab_g_pos(icell,lambda) = tab_g_pos(icell,lambda) + C_sca(k,lambda) * density * tab_g(k,lambda)

        k_abs_tot = k_abs_tot + C_abs(k,lambda) * density
     enddo !k

     if (kappa(icell,lambda) > tiny_real) tab_albedo_pos(icell,lambda) = kappa_sca(icell,lambda)/kappa(icell,lambda)
     if (kappa_sca(icell,lambda) > tiny_real) tab_g_pos(icell,lambda) = tab_g_pos(icell,lambda)/kappa_sca(icell,lambda)

     if (lRE_LTE) then
        kappa_abs_LTE(icell,lambda) = 0.0
        do k=grain_RE_LTE_start,grain_RE_LTE_end
           kappa_abs_LTE(icell,lambda) =  kappa_abs_LTE(icell,lambda) + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo
        k_abs_RE = k_abs_RE + kappa_abs_LTE(icell,lambda)
     endif

     if (lRE_nLTE) then
        kappa_abs_nLTE(icell,lambda) = 0.0
        do k=grain_RE_nLTE_start,grain_RE_nLTE_end
           kappa_abs_nLTE(icell,lambda) =  kappa_abs_nLTE(icell,lambda) + C_abs(k,lambda) * densite_pouss(k,icell)
        enddo
        k_abs_RE = k_abs_RE + kappa_abs_nLTE(icell,lambda)
     endif

     if (letape_th) then
        if (lnRE.and.(k_abs_tot > tiny_db)) then
           kappa_abs_RE(icell,lambda) = k_abs_RE
           proba_abs_RE(icell,lambda) = k_abs_RE/k_abs_tot
        endif

        if (k_abs_RE > tiny_db) then
           Proba_abs_RE_LTE(icell,lambda) = kappa_abs_LTE(icell,lambda) / (k_abs_RE)
        else ! the cell is probably empty
           Proba_abs_RE_LTE(icell,lambda) = 1.0
        endif
        if (lRE_nLTE) Proba_abs_RE_LTE_p_nLTE(icell,lambda) = 1.0 ! so far, might be updated if nRE --> qRE grains
     endif ! letape_th
  enddo !icell

  ! proba absorption sur une taille donnée
  if (lRE_nLTE.and. (.not.low_mem_th_emission_nLTE)) then
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

  ! Normalisation des opacites pour etre en AU^-1
  ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
  ! a² microns² -> 1e-8 cm²             \
  ! density en cm-3                      > reste facteur AU_to_cm * mum_to_cm**2 = 149595.0
  ! longueur de vol en AU = 1.5e13 cm   /
  ! fact =  pi * a * a * 149595.0
  ! les k_abs_XXX n'ont pas besoin d'etre normalise car tout est relatif
  fact = AU_to_cm * mum_to_cm**2
  kappa(:,lambda) = kappa(:,lambda) * fact
  kappa_sca(:,lambda) = kappa_sca(:,lambda) * fact

  if (lRE_LTE) kappa_abs_LTE(:,lambda) = kappa_abs_LTE(:,lambda) * fact
  if (lRE_nLTE) kappa_abs_nLTE(:,lambda) = kappa_abs_nLTE(:,lambda) * fact
  if (letape_th.and.lnRE) kappa_abs_RE(:,lambda) =  kappa_abs_RE(:,lambda) * fact

  if (lscatt_ray_tracing) then
     do thetaj=0,nang_scatt
        angle = real(thetaj)/real(nang_scatt)*pi
        tab_cos_scatt(thetaj) = cos(angle)
     enddo
  endif


  if (scattering_method==2) then
     call calc_local_scattering_matrices(lambda, p_lambda)
  else ! scattering_method ==1
     if (.not.low_mem_scattering) then ! we precompute the CDF (otherwise, we recalculate it live)
        do icell=1, p_n_cells
           ksca_CDF(0,icell,lambda)=0.0

           do  k=1,n_grains_tot
              ksca_CDF(k,icell,lambda) = ksca_CDF(k-1,icell,lambda) + C_sca(k,lambda) * densite_pouss(k,icell)
           enddo !k

           if  (ksca_CDF(n_grains_tot,icell,lambda) > tiny_real) then
              ksca_CDF(:,icell,lambda)= ksca_CDF(:,icell,lambda)/ ksca_CDF(n_grains_tot,icell,lambda)
           else
              ! a la surface, on peut avoir une proba de 0.0 partout
              ! dans ce cas, on decide qu'il n'y a que les plus petits grains
              ! rq : en pratique, la densite est trop faible pour qu'il y ait
              ! une diffusion a cet endroit.
              ksca_CDF(:,icell,lambda) = 1.0
           endif
        enddo ! icell
     endif ! .not.low_mem_scattering
  endif

  ! Supression scattering
  if (lno_scattering) then
     kappa = kappa_abs_LTE
     tab_albedo_pos = 0.0_db
  endif

  ! scattering = abs
  if (lqsca_equal_qabs) then
     kappa = 2.0_db * kappa_abs_LTE
     tab_albedo_pos = 0.5_db
  endif

  ! On remet la densite à zéro si besoin
  if (ldens0) then
     icell = icell_ref
     densite_pouss(:,icell) = 0.0_db
     kappa(icell,lambda) = 0.0_db
     if (lRE_LTE) then
        kappa_abs_LTE(icell,lambda) = 0.0_db
     endif
     if (lcompute_obs.and.lscatt_ray_tracing.or.lProDiMo2mcfost) then
        kappa_sca(icell,lambda) = 0.0_db
     endif
     if (lRE_nLTE) then
        kabs_nLTE_CDF(:,icell,lambda) = 0.0
     endif
  endif

  if ((ldust_prop).and.(lambda == n_lambda)) then
     call write_dust_prop()

     if (lstop_after_init) then
        write(*,*) "Exiting"
        stop
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

  real(kind=db), parameter :: dtheta = pi/real(nang_scatt)
  real(kind=db) :: density, theta, norme, fact, k_sca_tot

  integer :: icell, k, l

  fact = AU_to_cm * mum_to_cm**2
  !write(*,*) "Computing local scattering properties", lambda, p_lambda

  !$omp parallel &
  !$omp default(none) &
  !$omp shared(tab_s11_pos,tab_s12_o_s11_pos,tab_s33_o_s11_pos,tab_s34_o_s11_pos) &
  !$omp shared(tab_s11,tab_s12,tab_s33,tab_s34,lambda,p_lambda,n_grains_tot,tab_albedo_pos,prob_s11_pos) &
  !$omp shared(zmax,kappa,kappa_sca,kappa_abs_LTE,ksca_CDF,p_n_cells,fact) &
  !$omp shared(C_ext,C_sca,densite_pouss,S_grain,scattering_method,tab_g_pos,aniso_method,tab_g,lisotropic,low_mem_scattering) &
  !$omp shared(lscatt_ray_tracing,letape_th,lsepar_pola,ldust_prop) &
  !$omp private(icell,k,density,norme,theta,k_sca_tot)
  !$omp do schedule(dynamic,1)
  do icell=1, p_n_cells
     if (aniso_method==1) then
        tab_s11_pos(:,icell,p_lambda) = 0.
        if (lsepar_pola) then
           tab_s12_o_s11_pos(:,icell,p_lambda) = 0.
           tab_s33_o_s11_pos(:,icell,p_lambda) = 0.
           tab_s34_o_s11_pos(:,icell,p_lambda) = 0.
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
              tab_s12_o_s11_pos(:,icell,p_lambda) = tab_s12_o_s11_pos(:,icell,p_lambda) + &
                   tab_s12(:,k,lambda) * S_grain(k) * density
              tab_s33_o_s11_pos(:,icell,p_lambda) = tab_s33_o_s11_pos(:,icell,p_lambda) + &
                   tab_s33(:,k,lambda) * S_grain(k) * density
              tab_s34_o_s11_pos(:,icell,p_lambda) = tab_s34_o_s11_pos(:,icell,p_lambda) + &
                   tab_s34(:,k,lambda) * S_grain(k) * density
           endif
        endif !aniso_method
     enddo !k

     k_sca_tot = kappa_sca(icell,lambda) / fact ! We renormalize to remove the factor from opacity

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
                    tab_s33_o_s11_pos(l,icell,p_lambda)=tab_s33_o_s11_pos(l,icell,p_lambda)*norme
                    tab_s34_o_s11_pos(l,icell,p_lambda)=tab_s34_o_s11_pos(l,icell,p_lambda)*norme
                 endif
              endif
           enddo

           ! Normalisation : on veut que l'energie total diffusee sur [0,pi] en theta et [0,2pi] en phi = 1
           ! (for ray-tracing)
           tab_s11_pos(:,icell,p_lambda) = tab_s11_pos(:,icell,p_lambda) * dtheta / (k_sca_tot * deux_pi)


           ! -- ! aniso_method = 2 --> HG
           ! -- gsca = tab_g_pos(icell,lambda)
           ! -- tab_s11_ray_tracing(:,icell,lambda) =  tab_s11_ray_tracing(:,icell,lambda) / (k_sca_tot * deux_pi)
           ! --
           ! -- if (lisotropic) tab_s11_ray_tracing(:,icell,lambda) = 1.0 / (4.* nang_scatt)
        else !aniso_method == 2 : HG
           tab_s11_pos(:,icell,p_lambda) = 1.0
           if (lsepar_pola) then
              tab_s12_o_s11_pos(:,icell,p_lambda)=0.0
              tab_s33_o_s11_pos(:,icell,p_lambda)=0.0
              tab_s34_o_s11_pos(:,icell,p_lambda)=0.0
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
        ! Propriétés optiques des cellules
        prob_s11_pos(:,icell,p_lambda)=1.0
        prob_s11_pos(0,icell,p_lambda)=0.0
        ! Normalisation (idem que dans mueller)
        tab_s11_pos(:,icell,p_lambda)=1.0
        if (lsepar_pola) then
           tab_s12_o_s11_pos(:,icell,p_lambda)=0.0
           tab_s33_o_s11_pos(:,icell,p_lambda)=0.0
           tab_s34_o_s11_pos(:,icell,p_lambda)=0.0
        endif

     endif ! k_sca_tot > or = 0
  enddo !icell
  !$omp enddo
  !$omp end parallel

  return

end subroutine calc_local_scattering_matrices

!******************************************************************************

end module dust
