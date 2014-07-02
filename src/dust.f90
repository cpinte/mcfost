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

subroutine taille_grains()
! Distribution en taille des grains
! Sortie de prop_grains pour cas multi-lambda
! C. Pinte 14/01/05

  implicit none
  integer :: k, i, j
  real :: a, alfa, qext=0.0, qsca=0.0, M_tot, nbre_tot_grains
  real(kind=db) :: exp_grains
  real :: masse_pop
  real :: correct_fact_r, correct_fact_S, correct_fact_M

  type(dust_pop_type), pointer :: dp

  integer :: ios, status, n_comment, n_grains
  real :: fbuffer


  ! Boucle sur les populations de grains
  do i=1, n_pop
     dp => dust_pop(i)

     if (lread_grain_size_distrib) then
        if (n_pop > 1) then
           write(*,*) "ERROR : you cannot provide a grain size distribution with more than 1 population"
           write(*,*) "Exiting"
           stop
        endif

        open(unit=1, file=grain_size_file, status='old', iostat=ios)
        if (ios/=0) then
           write(*,*) "ERROR : cannot open "//trim(grain_size_file)
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
           write(*,*) "ERROR : the number of grains must be the same as in the parameter file"
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
           grain(k)%pop = i
           masse_pop = masse_pop + nbre_grains(k)
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
           do  k=dp%ind_debut, dp%ind_fin
              if (k==dp%ind_debut) then
                 a = dp%amin*sqrt(exp_grains)
              else
                 a= r_grain(k-1) * exp_grains
              endif

              r_grain(k) = a ! micron
              S_grain(k) = pi * a**2 ! micron^2
              M_grain(k) = quatre_tiers_pi * (a*mum_to_cm)**3 * dp%rho1g_avg ! masse en g

              ! Multiplication par a car da = a.dln(a)
              nbre_grains(k) = a**(-dp%aexp) * a
              grain(k)%methode_chauffage = dp%methode_chauffage
              grain(k)%zone = dp%zone
              grain(k)%pop = i
              masse_pop = masse_pop + nbre_grains(k)
           enddo !k

        endif ! lread_grain_size_distribution

        masse_pop = masse_pop * dp%avg_grain_mass


        ! Normalisation du nombre de grains pour atteindre la bonne masse
        nbre_grains(dp%ind_debut:dp%ind_fin) = nbre_grains(dp%ind_debut:dp%ind_fin) * dp%masse/masse_pop
        ! we need to correct and divide by rho1g when reading the material density afterwards

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

end subroutine taille_grains

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
  integer :: n, i, k, ii, syst_status, alloc_status, pop, status, n_ind, buffer, ios, n_comment, n_components

  character(len=512) :: filename, dir

  real, dimension(:), allocatable :: tab_l, tab_n, tab_k
  real, dimension(:,:), allocatable :: tab_tmp_amu1, tab_tmp_amu2

  logical :: l_ordre_decroissant


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
                 write(*,*) "ERROR : coating can only be computed with 2 components"
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

        ! Force a temporary density, we will correct for it in read_opacity_file
        dust_pop(pop)%component_rho1g(1) = 1.0
        dust_pop(pop)%rho1g_avg = 1.0

        if (dust_pop(pop)%n_components > 1) then
           write(*,*) "ERROR : cannot mix PAH with other component"
           write(*,*) "Exiting"
           stop
        endif
     endif ! fin test fichier opacite

  enddo ! pop

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

  integer :: n, i, iter, k
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

subroutine prop_grains(lambda, p_lambda)
  ! Calcule les propriétés des grains
  ! Masse, fraction en nombre, sections efficaces, matrices de muellers
  ! sortie des routines opacite
  ! C. Pinte 7/01/05

  implicit none

  integer, intent(in) :: lambda, p_lambda
  real, parameter :: pi = 3.1415926535
  real :: a, wavel, alfa, qext, qsca, fact, gsca, amu1, amu2, amu1_coat, amu2_coat
  integer :: k, i, pop, l

  qext=0.0
  qsca=0.0

  ! Longueur d'onde
  wavel=tab_lambda(lambda)

  ! Prop optiques
  ! Une premiere boucle pour les grains definis par un fichier d'indice
  !$omp parallel &
  !$omp default(none) &
  !$omp private(k,a,alfa,qext,qsca,fact,gsca,amu1,amu2,pop) &
  !$omp shared(r_grain,q_ext,q_sca,q_abs,wavel,aexp,tab_albedo,lambda,p_lambda,tab_g,grain) &
  !$omp shared(laggregate,tab_amu1,tab_amu2,n_grains_tot,is_pop_PAH,is_grain_PAH) &
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
        alfa = 2.0 * pi * a / wavel
        if (laggregate) then
           call mueller_gmm(p_lambda,k,alfa,qext,qsca,gsca)
        else
           if ((dust_pop(pop)%type=="Mie").or.(dust_pop(pop)%type=="mie").or.(dust_pop(pop)%type=="MIE")) then
              call mueller2(p_lambda,k,alfa,amu1,amu2, qext,qsca,gsca)
              if (dust_pop(pop)%lcoating) then
                 amu1_coat=tab_amu1_coating(lambda,pop)
                 amu2_coat=tab_amu2_coating(lambda,pop)
                 call mueller_coated_sphere(p_lambda,k,wavel,amu1,amu2,amu1_coat,amu2_coat, qext,qsca,gsca)
              endif
           else if ((dust_pop(pop)%type=="DHS").or.(dust_pop(pop)%type=="dhs")) then
              call mueller_DHS(p_lambda,k,wavel,amu1,amu2, qext,qsca,gsca)
           else
              write(*,*) "Unknow dust type : ", dust_pop(pop)%type
              write(*,*) "Exiting"
              stop
           endif
           !           write(*,*) wavel, qext,qsca,gsca
        endif ! laggregate
        tab_albedo(lambda,k)=qsca/qext
        tab_g(lambda,k) = gsca
        ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
        ! a² microns² -> 1e-8 cm²             \
        ! density en cm-3                      > reste facteur 149595.0
        ! longueur de vol en AU = 1.5e13 cm   /
        fact =  pi * a * a * 149595.0
        !q_geo(k) = pi * a * a * 1.e-12 ! en m^2
        q_ext(lambda,k) = qext * fact ! todo : renommer C_ext
        q_sca(lambda,k) = qsca * fact
        q_abs(lambda,k) = q_ext(lambda,k) - q_sca(lambda,k) ! section efficace
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
        if (dust_pop(pop)%is_PAH) is_grain_PAH(k) = .true.
        call mueller_opacity_file(lambda,p_lambda,k,qext, qsca,gsca)

        tab_albedo(lambda,k)=qsca/qext
        tab_g(lambda,k) = gsca
        ! tau est sans dimension : [kappa * lvol = density * a² * lvol]
        ! a² microns² -> 1e-8 cm²             \
        ! density en cm-3                      > reste facteur 149595.0
        ! longueur de vol en AU = 1.5e13 cm   /
        fact =  pi * a * a * 149595.0
        !q_geo(k) = pi * a * a * 1.e-12 ! en m^2
        q_ext(lambda,k) = qext * fact ! todo : renommer C_ext
        q_sca(lambda,k) = qsca * fact
        q_abs(lambda,k) = q_ext(lambda,k) - q_sca(lambda,k)
     endif ! is_opacity_file
  enddo !k

  return

end subroutine prop_grains

!******************************************************************************

subroutine save_dust_prop(letape_th)

  logical, intent(in) :: letape_th

  type(dust_pop_type), dimension(n_pop) :: dust_pop_save
  character(len=512) :: filename, tab_wavelength_save
  integer :: n_lambda_save
  real :: lambda_min_save, lambda_max_save

  dust_pop_save = dust_pop
  n_lambda_save = n_lambda
  lambda_min_save = lambda_min
  lambda_max_save = lambda_max
  tab_wavelength_save = tab_wavelength

  if (letape_th) then
     filename=".dust_prop1.tmp" ;
  else
     filename=".dust_prop2.tmp" ;
  endif

  open(1,file=filename,status='replace',form='unformatted')
  write(1) dust_pop_save, q_ext, q_sca, q_abs, tab_g, tab_albedo, prob_s11, tab_s11, tab_s12, tab_s33, tab_s34, &
       n_lambda_save, lambda_min_save, lambda_max_save, tab_wavelength_save
  close(unit=1)

  return

end subroutine save_dust_prop

!******************************************************************************

subroutine read_saved_dust_prop(letape_th, lcompute)

  logical, intent(in) :: letape_th
  logical, intent(out) :: lcompute

  type(dust_pop_type), dimension(n_pop) :: dust_pop_save
  character(len=512) :: filename, tab_wavelength_save
  integer :: n_lambda_save
  real :: lambda_min_save, lambda_max_save

  integer :: i, pop, ios
  logical :: ok

  lcompute = .true.

  if (letape_th) then
     filename=".dust_prop1.tmp" ;
  else
     filename=".dust_prop2.tmp" ;
  endif

  ! check if there is a dust population file
  ios = 0
  open(1,file=filename,status='old',form='unformatted',iostat=ios)
  if (ios /= 0)  then
     close(unit=1)
     return
  endif

  ! read the saved dust properties
  read(1,iostat=ios)  dust_pop_save, q_ext, q_sca, q_abs, tab_g, tab_albedo, prob_s11, tab_s11, tab_s12, tab_s33, tab_s34, &
       n_lambda_save, lambda_min_save, lambda_max_save, tab_wavelength_save
  close(unit=1)
  if (ios /= 0)  return

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

  lcompute = .false.
  return

end subroutine read_saved_dust_prop

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
!  real, dimension(1:n_prob) :: xa,ya,y2
!  real :: delta_prob1, delta_probn, yp1, ypn
!  real :: yp_p1, yp_m1, ypp



  real :: z, rcyl, somme, gsca

  logical :: ldens0

  real(kind=db) :: kappa_abs_RE, kappa_abs_tot, angle

  real, dimension(:), allocatable :: kappa_lambda,albedo_lambda,g_lambda
  real, dimension(:,:), allocatable :: S11_lambda_theta, pol_lambda_theta

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

           if (lcompute_obs.and.lscatt_ray_tracing.or.lProDiMo2mcfost) then
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
  !$omp private(i,j,k,pk,density,k_min,proba,q_sca_tot,q_ext_tot,norme,angle,gsca)&
  !$omp shared(zmax,kappa,kappa_abs_eg,probsizecumul,ech_prob,p_n_rad,p_nz,p_n_az,j_start,pj_start) &
  !$omp shared(q_ext,q_sca,densite_pouss,scattering_method,tab_g_pos,aniso_method,tab_g,lisotropic) &
  !$omp shared(lscatt_ray_tracing,tab_s11_ray_tracing,tab_s12_ray_tracing,tab_s33_ray_tracing,tab_s34_ray_tracing) &
  !$omp shared(tab_s12_o_s11_ray_tracing,tab_s33_o_s11_ray_tracing,tab_s34_o_s11_ray_tracing,lsepar_pola,ldust_prop)
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


           if (q_ext_tot > tiny_real) tab_albedo_pos(lambda,i,j,pk) = q_sca_tot/q_ext_tot
           if (q_sca_tot > tiny_real) tab_g_pos(lambda,i,j,pk) = tab_g_pos(lambda,i,j,pk)/q_sca_tot



           if (lcompute_obs.and.lscatt_ray_tracing) then
              if (scattering_method == 1) then
                 write(*,*) "ERROR: ray-tracing is incompatible with scattering method 1"
                 stop
              endif

              if (aniso_method == 1) then
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
                    tab_s12_o_s11_ray_tracing(lambda,i,j,:) = tab_s12_ray_tracing(lambda,i,j,:) / &
                         max(tab_s11_ray_tracing(lambda,i,j,:),tiny_real)
                    tab_s33_o_s11_ray_tracing(lambda,i,j,:) = tab_s33_ray_tracing(lambda,i,j,:) / &
                         max(tab_s11_ray_tracing(lambda,i,j,:),tiny_real)
                    tab_s34_o_s11_ray_tracing(lambda,i,j,:) = tab_s34_ray_tracing(lambda,i,j,:) / &
                         max(tab_s11_ray_tracing(lambda,i,j,:),tiny_real)
                 endif
              else ! aniso_method = 2 --> HG
                 gsca = tab_g_pos(lambda,i,j,pk)

                 norme = 0.0
                 do thetaj=0,nang_scatt
                    angle = real(thetaj)/real(nang_scatt)*pi
                    tab_s11_ray_tracing(lambda,i,j,thetaj) =((1-gsca**2)/(2.0))* &
                         (1+gsca**2-2*gsca*cos((real(j))/real(nang_scatt)*pi))**(-1.5)
                    norme=norme + tab_s11_ray_tracing(lambda,i,j,thetaj) * sin(angle)
                 enddo
                 tab_s11_ray_tracing(lambda,i,j,:) =  tab_s11_ray_tracing(lambda,i,j,:) / (norme * deux_pi)
              endif

              if (lisotropic) tab_s11_ray_tracing(lambda,i,j,:) = 1.0 / (4.* nang_scatt)

            !  ! Verification normalization
            !  norme = 0.0
            !  do thetaj=0,nang_scatt
            !     angle = real(thetaj)/real(nang_scatt)*pi
            !     norme=norme + tab_s11_ray_tracing(lambda,i,j,thetaj) * sin(angle)
            !  enddo
            !  write(*,*) "test norme ", norme * deux_pi, gsca ! doit faire 1
           endif !lscatt_ray_tracing

           if (sum(densite_pouss(i,j,pk,:)) > tiny_real) then

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
                          if (.not.ldust_prop) then
                             tab_s11_pos(lambda,i,j,pk,l)= 1.0 !tab_s11_pos(lambda,i,j,pk,l)*norme
                          endif
                          if (lsepar_pola) then
                             tab_s12_pos(lambda,i,j,pk,l)=tab_s12_pos(lambda,i,j,pk,l)*norme
                             tab_s33_pos(lambda,i,j,pk,l)=tab_s33_pos(lambda,i,j,pk,l)*norme
                             tab_s34_pos(lambda,i,j,pk,l)=tab_s34_pos(lambda,i,j,pk,l)*norme
                          endif
                       endif
                    enddo
                 else !aniso_method
                    tab_s11_pos(lambda,i,j,pk,:) = 1.0
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

!  write(*,*) sum(densite_pouss),   q_abs(lambda,1), r_grain(1)
!  stop

  if ((ldust_prop).and.(lambda == n_lambda)) then
     write(*,*) "Writing dust propreties"

     ! Only do it after the last pass through the wavelength table
     ! in order to populate the tab_s11_pos and tab_s12_pos tables first!
     allocate(kappa_lambda(n_lambda))
     allocate(albedo_lambda(n_lambda))
     allocate(g_lambda(n_lambda))
     allocate(S11_lambda_theta(n_lambda,0:nang_scatt),pol_lambda_theta(n_lambda,0:nang_scatt))

     kappa_lambda=real((kappa(:,1,1,1)/AU_to_cm)/(masse(1,1,1)/(volume(1)*AU_to_cm**3))) ! cm^2/g
     albedo_lambda=tab_albedo_pos(:,1,1,1)
     g_lambda=tab_g_pos(:,1,1,1)

     call cfitsWrite("data_dust/lambda.fits.gz",tab_lambda,shape(tab_lambda))
     call cfitsWrite("data_dust/kappa.fits.gz",kappa_lambda,shape(kappa_lambda))
     call cfitsWrite("data_dust/albedo.fits.gz",albedo_lambda,shape(albedo_lambda))
     call cfitsWrite("data_dust/g.fits.gz",g_lambda,shape(g_lambda))

     S11_lambda_theta(:,:)= tab_s11_pos(:,1,1,1,:)
     call cfitsWrite("data_dust/phase_function.fits.gz",S11_lambda_theta,shape(S11_lambda_theta))

     if (lsepar_pola) then
        pol_lambda_theta(:,:)=-tab_s12_pos(:,1,1,1,:)/tab_s11_pos(:,1,1,1,:)
        call cfitsWrite("data_dust/polarizability.fits.gz",pol_lambda_theta,shape(pol_lambda_theta))
     endif

     deallocate(kappa_lambda,albedo_lambda,g_lambda,pol_lambda_theta)

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

end subroutine opacite2

!******************************************************************************

end module dust
