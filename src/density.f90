module density

  use parametres
  use constantes
  use molecular_emission
  use opacity
  use em_th
  use prop_star
  use grains
  use disk
  use nr, only : rkqs
  use nrtype, only : sp
  use ode_data
  use wall
  use grid
  use utils
  use output

  implicit none

  contains

subroutine define_density()

  call define_gas_density()
  call define_dust_density()

  return

end subroutine define_density

!*************************************************************

subroutine define_gas_density()
  ! Calcule la table de densite: densite_gaz et masse_gaz
  ! C. Pinte : re-ecrit le 27/04/2013

  implicit none

  integer :: i,j, k, izone, alloc_status
  real(kind=db), dimension(n_zones) :: cst_gaz
  real(kind=db) :: z, density, fact_exp, rsph, mass, puffed, facteur, z0, phi, surface, H, C

  type(disk_zone_type) :: dz

  ! Tableau temporaire pour densite gaz dans 1 zone (pour renormaliser zone par zone)
  ! Pas besoin dans la poussiere car a chaque pop, il y a des tailles de grains independantes
  real(kind=db), dimension(:,:,:), allocatable :: densite_gaz_tmp


  if (l3D) then
     allocate(densite_gaz_tmp(n_rad,-nz:nz,n_az), stat=alloc_status)
  else
     allocate(densite_gaz_tmp(n_rad,0:nz,1), stat=alloc_status)
  endif
  densite_gaz_tmp = 0.0
  densite_gaz = 0.0 ;

  do izone=1, n_zones
     dz = disk_zone(izone)
     if (dz%geometry <= 2) then ! Disque
        ! La constante n'est OK que pour les lois de puissance
        ! mais sera renormalisee pour les tappered-edge, inner-edge, etc
        if (abs(dz%surf+2.0) > 1.0e-5) then
           C =(((2.0+dz%surf))/(((2*pi)**1.5)*dz%sclht*(dz%rin**2)*&
                ((dz%rmax/dz%rin)**(2+dz%surf)-1)))*((dz%rref/dz%rin)**dz%surf)
        else
           C =1.0/((2*pi)**1.5*dz%sclht*(dz%rref**2) * log(dz%rmax/dz%rin))
        endif
     else ! envelope
        if (abs(dz%surf+3.0) > 1.0e-5) then
           C =((3.0+dz%surf))/ (4*pi*(dz%rmax**(3.+dz%surf) - dz%rin**(3.+dz%surf)))
        else
           C = 1.0/( 4*pi * log(dz%rmax/dz%rin) )
        endif
     end if

     ! Facteur multiplicatif pour passer en masse de gaz
     ! puis en nH2/AU**3 puis en nH2/m**3
     cst_gaz(izone) = C * dz%diskmass * dz%gas_to_dust / masse_mol_gaz * Msun_to_g / AU3_to_m3
  enddo

  do izone=1, n_zones
     dz=disk_zone(izone)
     densite_gaz_tmp = 0.0

     if (dz%geometry <= 2) then ! Disque
        do i=1, n_rad
           bz : do j=min(0,j_start),nz
              !if (j==0) cycle bz

              ! On calcule la densite au milieu de la cellule
              if (j==0) then
                 rcyl = r_grid(i,1)
                 z = 0.0_db
              else
                 rcyl = r_grid(i,j)
                 z = z_grid(i,j)
              endif

              H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta

              ! Puffed-up rim analytique
              if (lpuffed_rim) then
                 puffed = 1.0 + (puffed_rim_h - 1.0) / (exp((rcyl - puffed_rim_r)/puffed_rim_delta_r) + 1.0)
              else
                 puffed = 1.0
              endif

              if (dz%geometry == 1) then ! power-law
                 fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta)
              else ! tappered-edge : dz%surf correspond a -gamma
                 fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta) * exp( -(rcyl/dz%rc)**(2+dz%moins_gamma_exp) )
              endif
              coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))

              do k=1, n_az
                 phi = phi_grid(k)
                 ! Warp analytique
                 if (lwarp) then
                    z0 = z_warp * (rcyl/dz%rref)**3 * cos(phi)
                 else
                    z0 = 0.0
                 endif

                 !  Densite du gaz (lsetup_gas avant)
                 if (rcyl > dz%rmax) then
                    density = 0.0
                 else if (rcyl < dz%rmin) then
                    density = 0.0
                 else if (rcyl < dz%rin) then
                    density = cst_gaz(izone)*fact_exp * exp(-(((z-z0)/(dz%sclht*puffed))**2)/(coeff_exp))*&
                         exp(-((rcyl-dz%rin)**2)/(2.*dz%edge**2))
                 else
                    density = cst_gaz(izone)*fact_exp * exp(-(((z-z0)/(dz%sclht*puffed))**2)/(coeff_exp))
                 endif
                 densite_gaz_tmp(i,j,k) = density

              enddo !k
           enddo bz !j
        enddo ! i

     else if (dz%geometry == 3) then ! enveloppe : 2D uniquement pour le moment
        do i=1, n_rad
           do j=0,nz
              ! On calcule la densite au milieu de la cellule
              if (j==0) then
                 rcyl = r_grid(i,1)
                 z = 0.0_db
              else
                 rcyl = r_grid(i,j)
                 z = z_grid(i,j)
              endif
              rsph = sqrt(rcyl**2+z**2)


              ! Setup densite du gaz
              if (rsph > dz%rmax) then
                 density = 0.0
              else if (rsph < dz%rmin) then
                 density = 0.0
              else if (rsph < dz%rin) then
                 density = cst_gaz(izone) * rsph**(dz%surf)  * exp(-((rsph-dz%rin)**2)/(2.*dz%edge**2))
              else
                 density = cst_gaz(izone) * rsph**(dz%surf)
              endif
              densite_gaz_tmp(i,j,1) = density
           enddo !j
        enddo ! i

     endif ! dz%geometry

     !----------------------------------------
     ! Normalisation masse de gaz par zone
     !----------------------------------------

     ! Calcul de la masse de gaz de la zone
     mass = 0.
     do i=1,n_rad
        bz_gas_mass : do j=j_start,nz
           if (j==0) cycle bz_gas_mass
           do k=1,n_az
              mass = mass + densite_gaz_tmp(i,j,k) *  masse_mol_gaz * volume(i)
           enddo  !k
        enddo bz_gas_mass
     enddo !i
     mass =  mass * AU3_to_m3 * g_to_Msun

        ! Normalisation
     if (mass > 0.0) then ! pour le cas ou gas_to_dust = 0.
        facteur = dz%diskmass * dz%gas_to_dust / mass
        !     write(*,*) "VERIF gas mass: zone ",  izone, dz%diskmass * dz%gas_to_dust, mass, facteur

        ! Somme sur les zones pour densite finale
        do i=1,n_rad
           bz_gas_mass2 : do j=min(0,j_start),nz
              do k=1, n_az
                 densite_gaz(i,j,k) = densite_gaz(i,j,k) + densite_gaz_tmp(i,j,k) * facteur
              enddo !k
           enddo bz_gas_mass2
        enddo ! i
     endif
  enddo ! n_zones

  ! Ajout cavite vide
  if (lcavity) then
     do i=1, n_rad
        do j = 1, nz
           surface = cavity%sclht * (r_grid(i,j) / cavity%rref)**cavity%exp_beta
           if (z_grid(i,j) > surface) then
              densite_gaz(i,j,1) = 0.0_db
           endif
        enddo
     enddo
  endif

  ! Tableau de masse de gaz
  do i=1,n_rad
     facteur = masse_mol_gaz * volume(i) * AU3_to_m3
     bz_gas_mass3 : do j=j_start,nz
        if (j==0) cycle bz_gas_mass3
        do k=1, n_az
           masse_gaz(i,j,k) =  densite_gaz(i,j,k) * facteur
        enddo !k
     enddo bz_gas_mass3
  enddo ! i
  write(*,*) 'Total  gas mass in model:', real(sum(masse_gaz) * g_to_Msun),' Msun'

  if (lcorrect_density) then
     write(*,*) "Correcting density ..."
     do i=1, n_rad
        if ((r_grid(i,1) >= correct_density_Rin).and.(r_grid(i,1) <= correct_density_Rout)) then
           densite_gaz(i,:,:) = densite_gaz(i,:,:) *  correct_density_factor
           masse_gaz(i,:,:) = masse_gaz(i,:,:) *  correct_density_factor
        endif
     enddo
  endif

  return

end subroutine define_gas_density

!*************************************************************

subroutine define_dust_density()
! Calcule la table de densite
! Inclus stratification analytique
! Calcule les tableaux densite_pouss et masse
! et indique ri_not_empty, zj_not_empty et phik_not_empty
! C. Pinte : re-ecrit le 27/04/2013

  implicit none

  integer :: i,j, k, ii, jj, kk, l, k_min, izone, pop
  real(kind=db), dimension(n_pop) :: cst, cst_pous
  real(kind=db) :: lrin,  lrout, ledge, rcyl, rsph, mass
  real(kind=db) :: z, z_demi, fact_exp, coeff_exp, density, OmegaTau, h_H2, coeff_strat, proba
  real(kind=db) :: puffed, facteur, z0, phi, surface, norme, S

  real(kind=db), dimension(n_grains_tot) :: correct_strat, N_tot, N_tot2

  real(kind=db) :: rho, rho0, ztilde, dtilde, h, hd

  ! Pour puffed-up inner rim
  real(kind=db) :: correct_H

  real(kind=db) :: s_opt

  type(disk_zone_type) :: dz
  type(dust_pop_type), pointer :: dp

  ! Pour simus Seb
  real, parameter :: Sc = 1.5 ! nbre de Schmidt

  ! Pour simus Dubrulle
  real, parameter :: gamma = 2.0 ! exposant de turbulence, 2 pour turbulence compressible

  logical :: lwarning

  densite_pouss = 0.0; masse = 0.0

  ! Coefficient de diffusion constant
  Dtilde = alpha / Sc

  ! Cste densite disque : va eventuellement etre renormalise
  ! mais donne la constante pour la pop de poussiere dans une zone donnee
  do i=1, n_pop
     izone=dust_pop(i)%zone
     dz=disk_zone(izone)

     if ((1.0_db+1e-10_db) * dz%rin >=  dz%rmax) then
        write(*,*) "ERROR: Rout must be larger than than Rin in zone", izone
        write(*,*) "Exiting"
        stop
     endif

     if (dz%geometry == 4) lwall = .true.

     if (dz%geometry <= 2) then ! Disque
        if (abs(dz%surf+2.0) > 1.0e-5) then
           cst(i)=(((2.0+dz%surf)*dust_pop(i)%masse)/(((2*pi)**1.5)*dz%sclht*(dz%rin**2)*&
                ((dz%rmax/dz%rin)**(2+dz%surf)-1)))*((dz%rref/dz%rin)**dz%surf)
        else
           cst(i)=dust_pop(i)%masse/((2*pi)**1.5*dz%sclht*(dz%rref**2) * log(dz%rmax/dz%rin))
        endif
     else ! envelope
        if (abs(dz%surf+3.0) > 1.0e-5) then
           cst(i)=((3.0+dz%surf)*dust_pop(i)%masse)/ (4*pi*(dz%rmax**(3.+dz%surf) - dz%rin**(3.+dz%surf)))
        else
           cst(i)= dust_pop(i)%masse/( 4*pi * log(dz%rmax/dz%rin) )
        endif
     end if
  enddo !i pop

  ! facteur multiplicatif pour passer en g/cm**3:
  cst=cst * Msun_to_g / AU3_to_cm3

  ! facteur multiplicatif pour passer en part/cm**3: avg_grain_mass
  do i=1, n_pop
     cst_pous(i) = cst(i)/dust_pop(i)%avg_grain_mass
  enddo

  do  l=1,n_grains_tot
     ! Correction stratification
     if (lstrat.and.(settling_type == 1)) then
        ! loi de puissance
        if (r_grain(l) > a_strat) then
           correct_strat(l) = (r_grain(l)/a_strat)**exp_strat   ! (h_gas/h_dust)^2
        else
           correct_strat(l) = 1.0
        endif
        ! loi exponentielle (Garaud , Barriere 2004)
        !        correct_strat(k) = exp(r_grain(k)*fact_strat)/exp(amin*fact_strat)
     else
        correct_strat(l) = 1.0
     endif
  enddo !k

  ! Pour normalisation apres migration radiale
  N_tot(:) = 0.0 ; N_tot2(:) = 0.0 ;

  ! Boucle pop a l'exterieur pour pouvoir superposer differente pop
  do pop=1, n_pop
     izone=dust_pop(pop)%zone
     dz=disk_zone(izone)

     if (dz%geometry <= 2) then ! Disque

        do i=1, n_rad
           rho0 = densite_gaz(i,0,1) ! midplane density (j=0)

           !write(*,*) "     ", rcyl, rho0*masse_mol_gaz*cm_to_m**2, dust_pop(pop)%rho1g_avg
           !write(*,*) "s_opt", rcyl, s_opt/1000.

           bz : do j=j_start,nz
              if (j==0) cycle bz

              ! On calcule la densite au milieu de la cellule
              rcyl = r_grid(i,j)
              z = z_grid(i,j)

              H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta

              ! Puffed-up rim analytique
              if (lpuffed_rim) then
                 puffed = 1.0 + (puffed_rim_h - 1.0) / (exp((rcyl - puffed_rim_r)/puffed_rim_delta_r) + 1.0)
              else
                 puffed = 1.0
              endif

              if (dz%geometry == 1) then ! power-law
                 fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta)
              else ! tappered-edge : dz%surf correspond a -gamma
                 fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta) * exp( -(rcyl/dz%rc)**(2+dz%moins_gamma_exp) )
              endif
              coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))

              do k=1, n_az
                 phi = phi_grid(k)
                 ! Warp analytique
                 if (lwarp) then
                    z0 = z_warp * (rcyl/dz%rref)**3 * cos(phi)
                 else
                    z0 = 0.0
                 endif

                 ! Densite de la poussiere
                 do  l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    ! Settling a la Dubrulle
                    if (lstrat.and.(settling_type == 2)) then

                       !h_H=(1d0/(1d0+gamma))**(0.25)*sqrt(alpha/(Omega*tau_f)) ! echelle de hauteur du rapport gaz/poussiere / H_gaz
                       !hd_H=h_H*(1d0+h_H**2)**(-0.5)                           ! echelle de hauteur de la poussiere / H_gaz
                       OmegaTau = omega_tau(rho0,H,l)

                       h_H2= sqrt(1./(1.+gamma)) * alpha/OmegaTau
                       correct_strat(l) = (1 + h_H2) / h_H2 ! (h_gas/h_dust)^2
                    endif


                    if (rcyl > dz%rmax) then
                       density = 0.0
                    else if (rcyl < dz%rmin) then
                       density = 0.0
                    else if (rcyl < dz%rin) then
                       ! Strat radiale dans edge (pour GG Tau)
                       ! density = nbre_grains(l) * sqrt(correct_strat(l)) *  &
                       !      cst_pous(pop)*fact_exp * exp(-((z(j)/dz%sclht)**2*(correct_strat(l)))/(coeff_exp))*&
                       !      exp(-((rcyl-dz%rin)**2*(correct_strat(l)))/(2.*dz%edge**2))

                       ! Pas de strat radial dans edge
                       density = nbre_grains(l) * sqrt(correct_strat(l)) *  &
                            cst_pous(pop)*fact_exp * exp(-(((z-z0)/(dz%sclht*puffed))**2*(correct_strat(l)))/(coeff_exp))*&
                            exp(-((rcyl-dz%rin)**2)/(2.*dz%edge**2))
                    else
                       density = nbre_grains(l) * sqrt(correct_strat(l)) * cst_pous(pop)*fact_exp * &
                            exp(-(((z-z0)/(dz%sclht*puffed))**2*(correct_strat(l)))/(coeff_exp))
                    endif
                    densite_pouss(i,j,k,l) = density
                 enddo !l
              enddo !k
           enddo bz !j


           if (lstrat.and.(settling_type == 2)) then
              if (lspherical) then
                 write(*,*) "ERROR: settling following Dubrulle's prescription is only"
                 write(*,*) "implemented on a spherical grid"
                 write(*,*) "Exiting"
                 stop
              endif

              if ((rcyl > dz%rmin).and.(rcyl < dz%rmax)) then
                 ! Renormalisation pour  les cas ou il y a peu de resolution en z
                 do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    ! normalization en z
                    norme = 0.0
                    do j=j_start,nz
                       if (j==0) cycle
                       norme = norme + densite_pouss(i,j,1,l)
                    enddo !j

                    ! Si tous les grains sont sedimentes, on les met dans le plan median
                    if (norme < tiny_db) then
                       densite_pouss(i,1,1,l)  = 1.0_db
                       norme = 1.0_db

                       write(*,*) "WARNING: Vertical settling unresolved for"
                       write(*,*) "grain larger than", r_grain(l), "at R > ", real(rcyl)
                    endif

                    do j=j_start,nz
                       if (j==0) cycle
                       if (norme > tiny_db) densite_pouss(i,j,1,l) = densite_pouss(i,j,1,l) / norme * rho0 * nbre_grains(l)
                    enddo !j
                 enddo ! l
              endif ! test r
           endif ! settling==2

        enddo ! i

        if (lstrat.and.(settling_type == 3)) then
           ! Si strat a la Seb. Fromang : on ecrase le tableau de densite
           ! je ne code que la dependence en z dans un premier temps puis normalise et ajoute la dependence en R et taille de grain

           if (lspherical) then
              write(*,*) "ERROR: settling following Fromang's prescription is only"
              write(*,*) "implemented on a spharical grid"
              write(*,*) "Exiting"
              stop
           endif

           do i=1, n_rad
              lwarning = .true.
              rho0 = densite_gaz(i,0,1) ! pour dependance en R : pb en coord sperique

              rcyl = r_grid(i,1)
              H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta

              if ((rcyl > dz%rmin).and.(rcyl < dz%rmax)) then
                 ! Renormalisation pour  les cas ou il y a peu de resolution en z
                 do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    !calculate omega_tau in the disk midplane
                    OmegaTau = omega_tau(rho0,H,l)

                    do j=j_start,nz ! dependence en z uniquement ici !!!!
                       if (j==0) cycle

                       !calculate h & z/h
                       z = z_grid(i,j)
                       Ztilde=z/H

                       ! Fit Gaussien du profile de densite
                       !densite_pouss(i,j,1,l)=  exp(-(1+OmegaTau/dtilde) * (Ztilde**2/2.))

                       ! Coefficient de diffusion constant
                       densite_pouss(i,j,1,l)=  exp( -OmegaTau/dtilde * (exp(Ztilde**2/2.)-1) - Ztilde**2/2 )  ! formule 19
                    enddo!j

                    ! normalization en z
                    norme = 0.0
                    do j=j_start,nz
                       if (j==0) cycle
                       norme = norme + densite_pouss(i,j,1,l)
                    enddo !j

                    ! Si tous les grains sont sedimentes, on les met dans le plan median
                    if (norme < tiny_db) then
                       densite_pouss(i,1,1,l)  = 1.0_db
                       norme = 1.0_db

                       if (lwarning) then
                          write(*,*)
                          write(*,*) "WARNING : Vertical settling unresolved for"
                          write(*,*) "grain larger than", r_grain(l), "at R > ", real(rcyl)
                          lwarning = .false. ! on ne fait un warning qu'1 fois par rayon
                       endif
                    endif

                    do j=j_start,nz
                       if (j==0) cycle
                       if (norme > tiny_db) densite_pouss(i,j,1,l) = densite_pouss(i,j,1,l) / norme * rho0 * nbre_grains(l)
                    enddo !j

                 enddo ! l
              endif ! test r
           enddo ! i
        endif ! Settling Fromang

        lwarning = .true.
        if (lmigration) then
           do i=1, n_rad
              rho0 = densite_gaz(i,0,1) ! pour dependance en R : pb en coord sperique
              !s_opt = rho_g * cs / (rho * Omega)    ! cs = H * Omega ! on doit trouver 1mm vers 50AU
              !omega_tau= dust_pop(ipop)%rho1g_avg*(r_grain(l)*mum_to_cm) / (rho * masse_mol_gaz/m_to_cm**3 * H*AU_to_cm)

              rcyl = r_grid(i,1)
              H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta
              s_opt = (rho0*masse_mol_gaz*cm_to_m**3  /dust_pop(pop)%rho1g_avg) *  H * AU_to_m * m_to_mum

              if ((s_opt < dust_pop(pop)%amin).and.(lwarning)) then
                 write(*,*)
                 write(*,*) "WARNING: a_migration = ", s_opt
                 write(*,*) "is smaller than amin for dust pop #", pop
                 write(*,*) "MCFOST will exit with an error is there are no smaller grains"
                 if (s_opt < tiny_db) write(*,*) "is your gas-to-dust ratio = 0 ?"
                 lwarning = .false.
              endif

              !write(*,*) "s_opt", i, rcyl, s_opt

              do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                 S = sum(densite_pouss(i,:,:,l))
                 N_tot(l) = N_tot(l) + S
                 if (r_grain(l) > s_opt) then ! grains plus gros que taille optimale de migration
                    !write(*,*) "migration", i, rcyl, l, r_grain(l)
                    densite_pouss(i,:,:,l) = 0.0
                 else
                    N_tot2(l) = N_tot2(l) + S
                 endif
              enddo ! l
           enddo !i

           ! Renormalisation : on garde le meme nombre de grains par taille que avant la migration
           do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
              if (N_tot2(l) > tiny_db) then
                 densite_pouss(:,:,:,l) = densite_pouss(:,:,:,l) * N_tot(l)/N_tot2(l)
              endif
           enddo ! l
        endif ! migration

     else if (dz%geometry == 3) then ! enveloppe : 2D uniquement pour le moment
        do i=1, n_rad
           do j=1,nz
              ! On calcule la densite au milieu de la cellule
              rcyl = r_grid(i,j)
              z = z_grid(i,j)
              rsph = sqrt(rcyl**2+z**2)

              do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                 if (rsph > dz%rmax) then
                    density = 0.0
                 else if (rsph < dz%rmin) then
                    density = 0.0
                 else if (rsph < dz%rin) then
                    density = nbre_grains(l) * cst_pous(pop) * rsph**(dz%surf)  * exp(-((rsph-dz%rin)**2)/(2.*dz%edge**2))
                 else
                    density = nbre_grains(l) * cst_pous(pop) * rsph**(dz%surf)
                 endif
                 densite_pouss(i,j,1,l) = density
              enddo !l
           enddo !j
        enddo ! i

     endif ! dz%geometry

  enddo ! pop

  ! Ajout cavite vide
  if (lcavity) then
     do i=1, n_rad
        do j = 1, nz
           surface = cavity%sclht * (r_grid(i,j) / cavity%rref)**cavity%exp_beta
           if (z_grid(i,j) > surface) then
              densite_pouss(i,j,1,:) = 0.0_db
           endif
        enddo
     enddo
  endif

  if (lgap_ELT) then
     do i=1, n_rad
        densite_pouss(i,:,:,:) = densite_pouss(i,:,:,:) * &
             (1.0 - exp(-0.5 * ((r_grid(i,1) - r_gap_ELT) / sigma_gap_ELT)**2 ))
     enddo
  endif


  search_not_empty : do l=1,n_grains_tot
     do j=1,nz
        do i=1,n_rad
           if (densite_pouss(i,j,1,l) > 0.0_db) then
              ri_not_empty = i
              zj_not_empty = j
              phik_not_empty = 1
              exit search_not_empty
           endif
        enddo
     enddo
  enddo search_not_empty



  ! Normalisation poussiere: re-calcul masse totale par population a partir de la densite (utile quand edge /= 0)
  do pop=1, n_pop
     izone=dust_pop(pop)%zone
     dz=disk_zone(izone)

     if (dz%geometry /= 4) then ! pas de wall ici
        dp => dust_pop(pop)
        mass = 0.0

        do i=1,n_rad
           bz2 : do j=j_start,nz
              if (j==0) cycle bz2

              do k=1,n_az
                 do l=dp%ind_debut,dp%ind_fin
                    mass=mass + densite_pouss(i,j,k,l) * M_grain(l) * volume(i)
                 enddo !l
              enddo !k
           enddo bz2
        enddo !i
        mass =  mass * AU3_to_cm3 * g_to_Msun

        if (mass < tiny_db) then
           write(*,*)
           write(*,*) "ERROR : something went wrong, there is no dust in the disk"
           write(*,*) "Exiting." ; stop
        endif

        facteur = dp%masse / mass

        do i=1,n_rad
           bz3 : do j=j_start,nz
              if (j==0) cycle bz3
              do k=1, n_az
                 do l=dp%ind_debut,dp%ind_fin
                    densite_pouss(i,j,k,l) = densite_pouss(i,j,k,l) * facteur
                    masse(i,j,k) = masse(i,j,k) + densite_pouss(i,j,k,l) * M_grain(l) * volume(i)
                 enddo !l
              enddo !k
           enddo bz3
        enddo ! i

     endif ! test wall
  enddo ! pop

  masse(:,:,:) = masse(:,:,:) * AU3_to_cm3
  write(*,*) 'Total dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'

  if (ldust_gas_ratio) then
     allocate(dust_gas_ratio(n_rad,nz))
     dust_gas_ratio = masse(:,:,1) / masse_gaz(:,:,1)

     write(*,*) "Writing dust_gas_ratio"
     call cfitsWrite("dust_gas_ratio.fits.gz",dust_gas_ratio,shape(dust_gas_ratio))
     write(*,*) "Done"
     stop
  endif

  if (lcorrect_density) then
     write(*,*) "Correcting density ..."
     do i=1, n_rad
        if ((r_grid(i,1) >= correct_density_Rin).and.(r_grid(i,1) <= correct_density_Rout)) then
           densite_pouss(i,:,:,:) = densite_pouss(i,:,:,:) * correct_density_factor
           masse(i,:,:) = masse(i,:,:) *  correct_density_factor
        endif
     enddo

      write(*,*) 'Total corrected dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'
  endif

  ! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  ! Valable que dans le cas cylindrique mais pas de pb dans le cas spherique
  if (lcylindrical) densite_pouss(:,nz+1,:,:) = densite_pouss(:,nz,:,:)

  return

end subroutine define_dust_density

!*************************************************************

subroutine define_density_wall3D()
  ! superpose un mur avec une forme en cosinus
  ! sur le disque
  ! C. Pinte  25/01/2011

  integer :: pop, i, j, k, l, izone, alloc_status
  type(disk_zone_type) :: dz
  type(dust_pop_type), pointer :: dp

  real(kind=db) :: rcyl, z, phi, density, facteur, hh, mass

  real(kind=db), dimension(:,:,:,:), allocatable :: density_wall
  real(kind=db), dimension(:,:,:), allocatable :: masse_wall


  write(*,*) "*********************************************************"
  write(*,*) "Adding 3D wall structure ...."

  allocate(density_wall(n_rad,-nz-1:nz+1,n_az,n_grains_tot), stat=alloc_status)
  allocate(masse_wall(n_rad,-nz:nz,n_az), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error wall'
     stop
  endif
  density_wall = 0.
  masse_wall = 0.

  do pop=1, n_pop
     izone=dust_pop(pop)%zone
     dz=disk_zone(izone)


     if (dz%geometry == 3) then ! wall
        if (abs(dz%exp_beta) > 1e-6) then
           write(*,*) "ERROR: the wall must have Beta = 0"
           write(*,*) "Exiting"
           stop
        endif

        write(*,*) "Wall between", real(dz%rin), "and", real(dz%rmax), "AU"
        h_wall = real(dz%sclht)
        write(*,*) "h_wall =", real(h_wall)


        do i=1, n_rad
           bz : do j=j_start,nz
              if (j==0) cycle bz

              ! On calcule la densite au milieu de la cellule
              rcyl = r_grid(i,j)
              z = z_grid(i,j)

              do k=1, n_az
                 phi = phi_grid(k)

                 do  l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    if (rcyl > dz%rmax) then
                       density = 0.0
                    else if (rcyl < dz%rin) then
                       density = 0.0
                    else
                       density = nbre_grains(l) ! densite constante dans le mur
                    endif

                    ! Variation de hauteur du mur en cos(phi/2)
                    hh = h_wall * (1.+cos(phi+pi))/2.

                    if ((z > 0.).and.(z < hh)) then
                       density_wall(i,j,k,l) = density
                       !if (density > 0.) write(*,*) i,j, k, l, density_wall(i,j,k,l)
                    else
                       density_wall(i,j,k,l) = 0.0
                    endif
                 enddo ! l

              enddo !k
           enddo bz !j
        enddo !i

     endif ! wall
  enddo ! pop

  ! Normalisation de la masse du mur
  masse_wall(:,:,:) = 0.0

  do pop=1, n_pop
     izone=dust_pop(pop)%zone
     dz=disk_zone(izone)

     if (dz%geometry == 3) then ! wall
        dp => dust_pop(pop)
        mass = 0.0

        do i=1,n_rad
           bz2 : do j=j_start,nz
              if (j==0) cycle bz2

              do k=1,n_az
                 do l=dp%ind_debut,dp%ind_fin
                    mass=mass + density_wall(i,j,k,l) * M_grain(l) * (volume(i) * AU3_to_cm3)
                 enddo !l
              enddo !k
           enddo bz2
        enddo !i
        mass =  mass*g_to_Msun

        facteur = dp%masse / mass

        do i=1,n_rad
           bz3 : do j=j_start,nz
              if (j==0) cycle bz3
              do k=1, n_az
                 do l=dp%ind_debut,dp%ind_fin
                    density_wall(i,j,k,l) = density_wall(i,j,k,l) * facteur
                    masse_wall(i,j,k) = masse_wall(i,j,k) + density_wall(i,j,k,l) * M_grain(l) * volume(i)
                 enddo !l
              enddo !k
           enddo bz3
        enddo ! i

     endif ! zone = wall
  enddo ! pop

  masse_wall(:,:,:) = masse_wall(:,:,:) * AU3_to_cm3
  write(*,*) 'Wall dust mass:', real(sum(masse_wall)*g_to_Msun),' Msun'

  ! superposition du mur sur le disque
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) + density_wall(:,:,:,:)
  masse(:,:,:) = masse(:,:,:) + masse_wall(:,:,:)

  write(*,*) 'Total dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'
  write(*,*) "Done"
  write(*,*) "*********************************************************"


  return

end subroutine define_density_wall3D

!*************************************************************

subroutine densite_data2()
! Lecture du 2eme jeux de simulation de Laure Barriere
! La densite de surface n'est plus une loi de puissance
! Premier fits de GG Tau
! Mars 2005

  implicit none

  real, parameter :: G = 6.672e-8

  integer :: i,j, k, ii, jj, kk, l, k_min, nbre_a, alloc_status
  real ::  cst, cst_pous, cst_gaz, lrin,  lrout, ledge, rcyl, M_star, a, test
  real :: z_demi, fact_exp, coeff_exp, density, coeff_strat, proba, r0

  real, dimension(nz) :: z

  real, dimension(0:n_grains_tot) :: tab_cst, tab_surf, tab_beta, tab_h0
  real(kind=db) :: exp_grains, h, dsigma, ntot_grains, somme

  real, dimension(5) ::  a_sph,a1, a2, a0, b0, b1, b2, b3, b4
  real, dimension(n_grains_tot) ::  a1_int, a2_int, a0_int, b0_int, b1_int, b2_int, b3_int, b4_int

  real(kind=db) :: mass

  nbre_a=5

  ! Lecture données
  open(unit=1,file="/gagax1/ur2/cpinte/Laure/resultats_smooth/fits.dat",status="old")
!  open(unit=1,file="/gagax1/ur2/cpinte/Laure/resultats_smallsizes/fits.dat",status="old")
  do i=1,nbre_a
     read(1,*) a_sph(i)
     read(1,*) a2(i), a1(i), a0(i)
     read(1,*) b4(i),b3(i),b2(i), b1(i), b0(i)
  enddo
  close(unit=1)

  ! Interpolations
  l=0

  if (lstrat) then
     do k=1,n_grains_tot
        ! Petits grains
        if (r_grain(k) < a_sph(1)) then
           a0_int(k) = a0(1)
           a1_int(k) = a1(1)
           a2_int(k) = a2(1)
           b4_int(k) = b4(1)
           b3_int(k) = b3(1)
           b2_int(k) = b2(1)
           b1_int(k) = b1(1)
           b0_int(k) = b0(1)
           ! Gros grains
        else if (r_grain(k) > a_sph(nbre_a)) then
           a0_int(k) = a0(5)
           a1_int(k) = a1(5)
           a2_int(k) = a2(5)
           b4_int(k) = b4(5)
           b3_int(k) = b3(5)
           b2_int(k) = b2(5)
           b1_int(k) = b1(5)
           b0_int(k) = b0(5)
        else
           ! Autres grains : interpolation
           if (r_grain(k) > a_sph(l+1)) l = l+1
           a0_int(k) = a0(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(a0(l+1)-a0(l))
           a1_int(k) = a1(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(a1(l+1)-a1(l))
           a2_int(k) = a2(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(a2(l+1)-a2(l))
           b4_int(k) = b4(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(b4(l+1)-b4(l))
           b3_int(k) = b3(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(b3(l+1)-b3(l))
           b2_int(k) = b2(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(b2(l+1)-b2(l))
           b1_int(k) = b1(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(b1(l+1)-b1(l))
           b0_int(k) = b0(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(b0(l+1)-b0(l))
        endif
     enddo
  else !lstrat
     a0_int(:) = a0(1)
     a1_int(:) = a1(1)
     a2_int(:) = a2(1)
     b4_int(:) = b4(1)
     b3_int(:) = b3(1)
     b2_int(:) = b2(1)
     b1_int(:) = b1(1)
     b0_int(:) = b0(1)
  endif

  ! Constantes
  ! On divise ici par (AU/cm)**3 pour eviter overflow
  ntot_grains=diskmass*Msun_to_g/AU3_to_cm3 /dust_pop(1)%avg_grain_mass

  do i=1, n_rad
!     rcyl = i*1.0/real(resol)
     ! On calcule la densite au milieu de la cellule
     rcyl = sqrt(r_lim(i) * r_lim(i-1))

     ! Calcul opacite et probabilite de diffusion
     !$omp parallel &
     !$omp default(none) &
     !$omp shared(i,rcyl, tab_beta, tab_surf, tab_h0, tab_cst,nz,n_grains_tot) &
     !$omp shared(amax_reel,densite_pouss,z_lim,a1_int, a2_int, a0_int)&
     !$omp shared(b0_int, b1_int, b2_int, b3_int, b4_int,ntot_grains) &
     !$omp private(j,k,z,density,k_min,proba,h,dSigma, somme)&
     !$omp shared(zmax,kappa,probsizecumul,ech_prob,nbre_grains,cst_pous,q_ext,q_sca,delta_z)
     !$omp do schedule(dynamic,10)
     do j=1,nz
        z(j) = (real(j)-0.5)*delta_z(i)
        do  k=1,n_grains_tot
           H = 100.0*(a2_int(k)*(0.01*rcyl)**2 + a1_int(k) * (0.01*rcyl) + a0_int(k))
           dSigma =exp(b4_int(k)*(0.01*rcyl)**4+b3_int(k)*(0.01*rcyl)**3+ &
                b2_int(k)*(0.01*rcyl)**2+b1_int(k)*(0.01*rcyl)+b0_int(k))
           densite_pouss(i,j,1,k) = dSigma/(sqrt(2*pi)*H) * exp(-z(j)**2/(2*H**2))
!           write(*,*) densite_pouss(i,j,k)
        enddo !k
     enddo !j
     !$omp end do
     !$omp end parallel
  enddo !i


  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do k=1,n_grains_tot
     somme=0.0
     do i=1,n_rad
        do j=1,nz
           if (densite_pouss(i,j,1,k) < 0.0) then
              write(*,*) i,j,k
              read(*,*)
           endif
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo

  ! Normalisation : Calcul masse totale
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  ! Verif
!  masse = 0.0
!  do i=1,n_rad
!     do j=1,nz
!        do k=1,n_grains_tot
!           masse=masse + densite_pouss(i,j,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
!        enddo
!     enddo
!  enddo
!  masse =  masse/Msun_to_g
!  write(*,*) "masse totale =", masse


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_data2

!***********************************************************

subroutine densite_data_SPH_binaire()
! Lecture du 2eme jeux de simulation de Laure Barriere
! La densite de surface n'est plus une loi de puissance
! 2eme fits de GG Tau
! cree tableaux densite_pouss et masse
! 04/05/05

  implicit none

  real, parameter :: G = 6.672e-8

  integer :: i,j, k, ii, jj, kk, l, k_min, nbre_a, alloc_status
  real ::  cst, cst_pous, cst_gaz, lrin,  lrout, ledge, rcyl, M_star, a, test
  real :: z_demi, fact_exp, coeff_exp, density, coeff_strat, proba, r0

  real, dimension(nz) :: z

  real, dimension(0:n_grains_tot) :: tab_cst, tab_surf, tab_beta, tab_h0
  real(kind=db) :: exp_grains, h, dsigma, ntot_grains, somme

  real, dimension(5) ::  a_sph,a1, a2, a0, b0, b1, b2
  real, dimension(n_grains_tot) ::  a1_int, a2_int, a0_int, b0_int, b1_int, b2_int

  real(kind=db) :: mass, fact

  nbre_a=5

  write(*,*) "Strat SPH binaire"

  ! Lecture données
  open(unit=1,file=trim(home)//"/Laure/fits_GG_Tau/fits.dat",status="old")
  do i=1,nbre_a
     read(1,*) a_sph(i)
     read(1,*) a2(i), a1(i), a0(i)
     read(1,*) b2(i), b1(i), b0(i)
  enddo
  close(unit=1)

  ! Interpolations
  l=0

  if (.not.lno_strat_SPH_bin) then
     do k=1,n_grains_tot
        ! Petits grains
        if (r_grain(k) < a_sph(1)) then
           a0_int(k) = a0(1)
           a1_int(k) = a1(1)
           a2_int(k) = a2(1)
           b2_int(k) = b2(1)
           b1_int(k) = b1(1)
           b0_int(k) = b0(1)
           ! Gros grains
        else if (r_grain(k) > a_sph(nbre_a)) then
           a0_int(k) = a0(5)
           a1_int(k) = a1(5)
           a2_int(k) = a2(5)
           b2_int(k) = b2(5)
           b1_int(k) = b1(5)
           b0_int(k) = b0(5)
        else
           ! Autres grains : interpolation
           if (r_grain(k) > a_sph(l+1)) l = l+1
           !fact = (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))
           fact = (log(r_grain(k))-log(a_sph(l)))/(log(a_sph(l+1))-log(a_sph(l)))
           a0_int(k) = a0(l) + fact * (a0(l+1)-a0(l))
           a1_int(k) = a1(l) + fact * (a1(l+1)-a1(l))
           a2_int(k) = a2(l) + fact * (a2(l+1)-a2(l))
           b2_int(k) = b2(l) + fact * (b2(l+1)-b2(l))
           b1_int(k) = b1(l) + fact * (b1(l+1)-b1(l))
           b0_int(k) = b0(l) + fact * (b0(l+1)-b0(l))
        endif
     enddo
  else !lstrat
     a0_int(:) = a0(1)
     a1_int(:) = a1(1)
     a2_int(:) = a2(1)
     b2_int(:) = b2(1)
     b1_int(:) = b1(1)
     b0_int(:) = b0(1)
  endif

  ! Constantes
  ! On divise ici par (AU/cm)**3 pour eviter overflow
  ntot_grains=diskmass*Msun_to_g/AU3_to_cm3 /dust_pop(1)%avg_grain_mass

  do i=1, n_rad
!     rcyl = i*1.0/real(resol)
     ! On calcule la densite au milieu de la cellule
     rcyl = 0.5*(r_lim(i) + r_lim(i-1))

     ! Calcul opacite et probabilite de diffusion
     !$omp parallel &
     !$omp default(none) &
     !$omp shared(i,rcyl, tab_beta, tab_surf, tab_h0, tab_cst) &
     !$omp shared(amax_reel,densite_pouss,z_lim,a1_int, a2_int, a0_int)&
     !$omp shared(b0_int, b1_int, b2_int,ntot_grains,nz,n_grains_tot) &
     !$omp private(j,k,z,density,k_min,proba,h,dSigma, somme)&
     !$omp shared(zmax,kappa,probsizecumul,ech_prob,nbre_grains,cst_pous,q_ext,q_sca,delta_z)
     !$omp do schedule(dynamic,10)
     do j=1,nz
        z(j) = (real(j)-0.5)*delta_z(i)
        do  k=1,n_grains_tot
           H = 100.0*(a2_int(k)*(0.01*rcyl)**2 + a1_int(k) * (0.01*rcyl) + a0_int(k))
           dSigma = b2_int(k) * exp(-(log(0.01*rcyl)-b1_int(k)-log(0.01))**2/(2*b0_int(k)**2))
           densite_pouss(i,j,1,k) = dSigma/(sqrt(2*pi)*H) * exp(-z(j)**2/(2*H**2))
           if (densite_pouss(i,j,1,k) < 0.0) densite_pouss(i,j,1,k) = 0.0
!           write(*,*) densite_pouss(i,j,k)
           enddo !k
     enddo !j
     !$omp end do
     !$omp end parallel
  enddo !i


  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do k=1,n_grains_tot
     somme=0.0
     do i=1,n_rad
        do j=1,nz
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = densite_pouss(:,:,1,k)*nbre_grains(k)
  enddo

  ! Normalisation : Calcul masse totale
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  ! Verif
!  masse = 0.0
!  do i=1,n_rad
!     do j=1,nz
!        do k=1,n_grains_tot
!           masse=masse + densite_pouss(i,j,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
!        enddo
!     enddo
!  enddo
!  masse =  masse/Msun_to_g
!  write(*,*) "masse totale =", masse


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_data_SPH_binaire

!*****************************************************

subroutine densite_data_SPH_TTauri()
! Lecture du 2eme jeux de simulation de Laure Barriere
! La densite de surface n'est plus une loi de puissance
! 2eme fits de T Tauri
! 22 juin 2005

  implicit none

  real, parameter :: G = 6.672e-8

  integer, parameter :: nbre_a = 6
  integer :: i,j, k, ii, jj, kk, l, k_min, alloc_status
  real ::  cst, cst_pous, cst_gaz, lrin,  lrout, ledge, rcyl, M_star, a, test
  real :: z_demi, fact_exp, coeff_exp, density, coeff_strat, proba, r0

  real :: z

  real, dimension(0:n_grains_tot) :: tab_cst, tab_surf, tab_beta, tab_h0
  real(kind=db) :: exp_grains, ntot_grains, somme

  real, dimension(nbre_a) ::  a_sph, a0, a1, a2, a3, a4, b0, b1, b2, b3, b4
  real, dimension(n_grains_tot) ::  a1_int, a2_int, a0_int, a3_int, a4_int, b0_int, b1_int, b2_int, b3_int, b4_int
  real(kind=db), dimension(n_grains_tot) :: h, dsigma
  real(kind=db) :: mass, fact

  ! Lecture données
    write(*,*) "Reading SPH data ..."
!  open(unit=1,file=trim(home)//"/Laure/Donnees_TTauri_lr/fits.dat",status="old")
  open(unit=1,file=trim(home)//"/Laure/TTauri_hr/fits.dat",status="old")

!  open(unit=1,file="~/Laure/TTauri_hr/fits.dat",status="old")
  do i=1,nbre_a
     read(1,*) a_sph(i)
     a_sph(i) = a_sph(i) *1.0e6
     read(1,*) a4(i), a3(i), a2(i), a1(i), a0(i)
     read(1,*) b4(i),b3(i),b2(i), b1(i), b0(i)
  enddo
  close(unit=1)

  ! Interpolations
  l=0

  if (.not.lno_strat_SPH) then
     do k=1,n_grains_tot
        ! Petits grains
        if (r_grain(k) < a_sph(1)) then
           a0_int(k) = a0(1)
           a1_int(k) = a1(1)
           a2_int(k) = a2(1)
           a3_int(k) = a3(1)
           a4_int(k) = a4(1)
           b4_int(k) = b4(1)
           b3_int(k) = b3(1)
           b2_int(k) = b2(1)
           b1_int(k) = b1(1)
           b0_int(k) = b0(1)
           ! Gros grains
        else if (r_grain(k) > a_sph(nbre_a)) then
           a0_int(k) = a0(nbre_a)
           a1_int(k) = a1(nbre_a)
           a2_int(k) = a2(nbre_a)
           a3_int(k) = a3(nbre_a)
           a4_int(k) = a4(nbre_a)
           b4_int(k) = b4(nbre_a)
           b3_int(k) = b3(nbre_a)
           b2_int(k) = b2(nbre_a)
           b1_int(k) = b1(nbre_a)
           b0_int(k) = b0(nbre_a)
        else
           ! Autres grains : interpolation
           if (r_grain(k) > a_sph(l+1)) l = l+1
           fact = (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))
           !fact = (log(r_grain(k))-log(a_sph(l)))/(log(a_sph(l+1))-log(a_sph(l)))
           a0_int(k) = a0(l) + fact * (a0(l+1)-a0(l))
           a1_int(k) = a1(l) + fact * (a1(l+1)-a1(l))
           a2_int(k) = a2(l) + fact * (a2(l+1)-a2(l))
           a3_int(k) = a3(l) + fact * (a3(l+1)-a3(l))
           a4_int(k) = a4(l) + fact * (a4(l+1)-a4(l))
           b4_int(k) = b4(l) + fact * (b4(l+1)-b4(l))
           b3_int(k) = b3(l) + fact * (b3(l+1)-b3(l))
           b2_int(k) = b2(l) + fact * (b2(l+1)-b2(l))
           b1_int(k) = b1(l) + fact * (b1(l+1)-b1(l))
           b0_int(k) = b0(l) + fact * (b0(l+1)-b0(l))
        endif
     enddo
  else !lstrat
     a0_int(:) = a0(1)
     a1_int(:) = a1(1)
     a2_int(:) = a2(1)
     a3_int(:) = a3(1)
     a4_int(:) = a4(1)
     b4_int(:) = b4(1)
     b3_int(:) = b3(1)
     b2_int(:) = b2(1)
     b1_int(:) = b1(1)
     b0_int(:) = b0(1)
  endif


!  write(*,*) a0(1), a1(1), a2(1), a3(1), a4(1)


  ! Constantes
  ! On divise ici par (AU/cm)**3 pour eviter overflow
  ntot_grains=diskmass*Msun_to_g/AU3_to_cm3 /dust_pop(1)%avg_grain_mass

  do i=1, n_rad
!     rcyl = i*1.0/real(resol)
     ! On calcule la densite au milieu de la cellule
     rcyl = sqrt(r_lim(i) * r_lim(i-1))

     do k=1, n_grains_tot
!        H(k) = 100.0*exp(a4_int(k)*(0.01*rcyl)**4+a3_int(k)*(0.01*rcyl)**3+ &
!             a2_int(k)*(0.01*rcyl)**2 + a1_int(k) * (0.01*rcyl) + a0_int(k))
! Modif donnees 16/01/05
        H(k) = 100.0*(a4_int(k)*(0.01*rcyl)**4+a3_int(k)*(0.01*rcyl)**3+ &
             a2_int(k)*(0.01*rcyl)**2 + a1_int(k) * (0.01*rcyl) + a0_int(k))
        dSigma(k) =exp(b4_int(k)*(0.01*rcyl)**4+b3_int(k)*(0.01*rcyl)**3+ &
             b2_int(k)*(0.01*rcyl)**2+b1_int(k)*(0.01*rcyl)+b0_int(k))
!        H(k) = 10.0*(rcyl/100.0)**1.125
!        dSigma(k) = rcyl**(-1.0)
     enddo


     ! Calcul opacite et probabilite de diffusion
     !$omp parallel &
     !$omp default(none) &
     !$omp shared(i,rcyl, tab_beta, tab_surf, tab_h0, tab_cst,nz) &
     !$omp shared(amax_reel,densite_pouss,z_lim,a1_int, a2_int, a0_int, a3_int, a4_int)&
     !$omp shared(b0_int, b1_int, b2_int, b3_int, b4_int,ntot_grains,h,dSigma, n_grains_tot) &
     !$omp private(j,k,z,density,k_min,proba,somme)&
     !$omp shared(zmax,kappa,probsizecumul,ech_prob,nbre_grains,cst_pous,q_ext,q_sca,delta_z)
     !$omp do schedule(dynamic,10)
     do j=1,nz
        z = (real(j)-0.5)*delta_z(i)
        do k=1, n_grains_tot
           if (H(k) < 0.0) then
              densite_pouss(i,j,1,k) = 0.0
           else
              densite_pouss(i,j,1,k) = dSigma(k)/(sqrt(2*pi)*H(k)) * exp(-z**2/(2*H(k)**2))
           endif
        enddo !k
     enddo !j
     !$omp end do
     !$omp end parallel
  enddo !i


  ! Normalisation : on a 1 grain de chaque taille dans le disque
  somme=0.0
  do k=1,n_grains_tot
     do i=1,n_rad
        do j=1,nz
           if (densite_pouss(i,j,1,k) < 0.0) then
              write(*,*) "Pb densite"
              write(*,*) i,j,k
              read(*,*)
           endif
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo


  ! Normalisation : Calcul masse totale
  masse = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/sum(masse) * Msun_to_g
  masse =  masse* diskmass/sum(masse) * Msun_to_g


 ! Verif
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass
  write(*,*) "SPH data OK"
  write(*,*) "masse totale =", mass/ Msun_to_g


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_data_SPH_TTauri

!***********************************************************

subroutine densite_data_SPH_TTauri_1()
! Lecture du premier jeu de données de Laure Barriere
! Mai-Juin 2004

  implicit none

  real, parameter :: G = 6.672e-8

  integer :: i,j, k, ii, jj, kk, l, k_min, nbre_a, alloc_status
  real ::  cst, cst_pous, cst_gaz, lrin,  lrout, ledge, rcyl, M_star, a, test
  real :: z_demi, fact_exp, coeff_exp, density, coeff_strat, proba, r0

  real, dimension(nz) :: z

  real, dimension(0:n_grains_tot) :: tab_cst, tab_surf, tab_beta, tab_h0
  real :: exp_grains

  real, dimension(:,:), allocatable :: sph
  real, dimension(:), allocatable ::  a_sph, beta_sph, surf_sph, h0_sph, sigma_sph

  real ::  q_sca_tot, q_ext_tot,norme

  real(kind=db) :: somme, mass

  type(disk_zone_type) :: dz

  dz=disk_zone(1)

100 format(3X,I2,3X,F7.3)
101 format(e7.1,4(2X,F8.5))

  if (lno_strat_SPH) then
     open(unit=1,file=trim(home)//"/Laure/no_strat.dat",status='old')
  else
     open(unit=1,file=trim(home)//"/Laure/strat.dat",status='old')
  endif
  read(1,100) nbre_a, r0
  write(*,*) nbre_a, r0
  ! allocation des tab
  if(r0/=disk_zone(1)%rref) stop
  allocate(sph(nbre_a,5), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  sph = 0.0

  allocate(a_sph(nbre_a), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  a_sph = 0.0

  allocate(beta_sph(nbre_a), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  beta_sph = 0.0

  allocate(surf_sph(nbre_a), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  surf_sph = 0.0

  allocate(h0_sph(nbre_a), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  h0_sph = 0.0

  allocate(sigma_sph(nbre_a), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error'
     stop
  endif
  sigma_sph = 0.0

! On lit surf, beta, h0, r0 pour chaque taille de grains
  ! nbre de taille de grains et r0 (a=0.0 -> gaz)
  do l=1, nbre_a
     read(1,101) sph(l,1), sph(l,2), sph(l,3), sph(l,4), sph(l,5)
     write(*,*) sph(l,1), sph(l,2), sph(l,3), sph(l,4), sph(l,5)
     a_sph(l) = 1.e6*sph(l,1)
     beta_sph(l) = -sph(l,2)/2.0
     surf_sph(l) = sph(l,4) - sph(l,2)/2.0
     h0_sph(l) = 100./sqrt(2*exp(sph(l,3)))
     sigma_sph(l) = exp(sph(l,5) - sph(l,3)/2.0)
  enddo

  write(*,*) ' '
  write(*,*) ' a               beta          surf         h0      h0(rin) '
  do l=1,nbre_a
     write(*,*) a_sph(l), beta_sph(l), surf_sph(l), h0_sph(l), h0_sph(l)*(dz%rin/r0)**beta_sph(l)
  enddo

! Correction echelle de hauteur : au rayon interne, les grains les plus petits doivent etre les plus haut
  do l=1+1,nbre_a
     if(h0_sph(l)*(dz%rin/r0)**beta_sph(l) > h0_sph(1)*(dz%rin/r0)**beta_sph(1)) then
        h0_sph(l) = h0_sph(1)*(dz%rin/r0)**(beta_sph(1)-beta_sph(l))
     endif
  enddo

  write(*,*) ' '
  write(*,*) ' a               beta          surf         h0      h0(rin) '
  do l=1,nbre_a
     write(*,*) a_sph(l), beta_sph(l), surf_sph(l), h0_sph(l), h0_sph(l)*(dz%rin/r0)**beta_sph(l)
  enddo



! On interpole et on stocke dans des tables
  l = 0
  do k=1, n_grains_tot
     ! Petits grains
     if (r_grain(k) < a_sph(1)) then
        tab_beta(k) = beta_sph(1)
        tab_surf(k) = surf_sph(1)
        tab_h0(k) = h0_sph(1)
     ! Gros grains
     else if (r_grain(k) > a_sph(nbre_a)) then
        tab_beta(k) = beta_sph(nbre_a)
        tab_surf(k) = surf_sph(nbre_a)
        tab_h0(k) = h0_sph(nbre_a)
     else
        ! Autres grains : interpolation
        if (r_grain(k) > a_sph(l+1)) l = l+1
!        write(*,*) l, k
        tab_beta(k) = beta_sph(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(beta_sph(l+1)-beta_sph(l))
        tab_surf(k) = surf_sph(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(surf_sph(l+1)-surf_sph(l))
        tab_h0(k) = h0_sph(l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))*(h0_sph(l+1)-h0_sph(l))
     endif
!        write(*,*) k, r_grain(k), tab_beta(k), tab_surf(k),  tab_h0(k)
  enddo

  dz%exp_beta=beta_sph(1)
  dz%surf=surf_sph(1)
  dz%sclht=h0_sph(1)

  ! ATTENTION : cst et cst_pouss dependent de k
  do k=1, n_grains_tot
     ! Cste densite disque (depend de a par surf, ... et nbre_grains(k))
     if (abs(tab_surf(k)+2.0) > 1.0e-5) then
        cst=(((2+tab_surf(k))*diskmass)/(((2*pi)**1.5)*tab_h0(k)*(dz%rin**2) &
             *((dz%rmax/dz%rin)**(2+tab_surf(k))-1)))*((dz%rref/dz%rin)**tab_surf(k))
     else
        cst=diskmass/((2*pi)**1.5*tab_h0(k)*(dz%rref**2)) * log(dz%rin/dz%rmax)
     endif

     ! facteur multiplicatif pour passer en g/cm**3:
     cst=cst*Msun_to_g / AU3_to_cm3
     ! facteur multiplicatif pour passer en part/cm**3: AVG_GRAIN_MASS
     cst_pous = cst/dust_pop(1)%avg_grain_mass
     ! Facteur multiplicatif pour prendre en compte que les grains de taille a
     tab_cst(k) = cst_pous*nbre_grains(k)
  enddo

  somme = 0.0
  do i=1, n_rad
     ! On calcule la densite au milieu de la cellule
     rcyl = 0.5*(r_lim(i) + r_lim(i-1))
     fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta)
     coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))

     ! Calcul opacite et probabilite de diffusion
     ! $ omp parallel &
     ! $ omp default(none) &
     ! $ omp shared(i,rcyl, tab_beta, tab_surf, tab_h0, tab_cst) &
     ! $ omp shared(z_lim,densite_pouss,rref,nz) &
     ! $ omp private(j,k,z,density,k_min,proba,q_sca_tot,q_ext_tot,norme)&
     ! $ omp shared(zmax,kappa,probsizecumul,ech_prob,nbre_grains,cst_pous,q_ext,q_sca,delta_z,n_grains_tot)
     ! $ omp do schedule(dynamic,10)
     do j=1,nz
        z(j) = (real(j)-0.5)*delta_z(i)
        do  k=1,n_grains_tot
           densite_pouss(i,j,1,k) = tab_cst(k)*(rcyl/dz%rref)**(tab_surf(k)-tab_beta(k)) &
           *exp(-((z(j)/(tab_h0(k)))**2)/(2*(rcyl/dz%rref)**(2*tab_beta(k))))
            somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo !k
     enddo !j
     ! $ omp enddo
     ! $ omp end parallel
  enddo !i
  close(unit=1)

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo


  ! Normalisation : Calcul masse totale
  masse = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/sum(masse) * Msun_to_g
  masse =  masse* diskmass/sum(masse) * Msun_to_g


 ! Verif
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass
  write(*,*) "SPH data OK"
  write(*,*) "masse totale =", mass/ Msun_to_g


  !***************
  ! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_data_SPH_TTauri_1

!**********************************************************************

subroutine densite_data_SPH_TTauri_2()
! Lecture du 2eme jeux de simulation de Laure Barriere
! La densite de surface n'est plus une loi de puissance
! 3eme fits de T Tauri : calcul sur la grille de mcfost
! 13 Octobre 2006

  implicit none

  real, parameter :: G = 6.672e-8

  integer, parameter :: nbre_a = 6
  integer :: i,j, k, ii, jj, kk, l, k_min, alloc_status
  real ::  cst, cst_pous, cst_gaz, lrin,  lrout, ledge, rcyl, M_star, a, test
  real :: z_demi, fact_exp, coeff_exp, density, coeff_strat, proba, r0

  real :: z

  real, dimension(0:n_grains_tot) :: tab_cst, tab_surf, tab_beta, tab_h0
  real(kind=db) :: exp_grains, ntot_grains, somme

  real, dimension(nbre_a) ::  a_sph
  real, dimension(n_rad,nz,nbre_a) :: rho
  real :: HH, SS, rr, zz

  real, dimension(n_grains_tot) ::  a1_int, a2_int, a0_int, a3_int, a4_int, b0_int, b1_int, b2_int, b3_int, b4_int
  real(kind=db), dimension(n_grains_tot) :: h, dsigma
  real(kind=db) :: mass, fact

  character(len=512) :: s, dir

  ! Lecture données
  write(*,*) "Reading SPH data ..."
   do l=1, nbre_a
     write(s,fmt='(I1)') 7-l
     a_sph(l) = 10**(l-1);
     dir=trim(home)//"TTauri_hr_61107_nomin_nomax"
     !dir=trim(home)//"TTauri_hr_61015"
     write(*,*) trim(dir)//"/fitC_1d-"//trim(s)//"m.out"
     open(unit=1,file=trim(dir)//"/fitC_1d-"//trim(s)//"m.out",status="old")
     read(1,*)
     do i=1, n_rad
        do j=1, nz
           read(1,*) ii, jj, rr, zz, rho(i,j,l), HH, SS
        enddo
     enddo
     close(unit=1)
  enddo


  l=0
  ! Interpolation
  if (.not.lno_strat_SPH) then
     do k=1,n_grains_tot
        ! Petits grains
        if (r_grain(k) < a_sph(1)) then
           do i=1, n_rad
              do j=1,nz
                 densite_pouss(i,j,1,k) = rho(i,j,1)
              enddo
           enddo
           ! Gros grains
        else if (r_grain(k) > a_sph(nbre_a)) then
           do i=1, n_rad
              do j=1,nz
                 densite_pouss(i,j,1,k) = rho(i,j,1)
              enddo
           enddo
        else
           ! Autres grains : interpolation
           if (r_grain(k) > a_sph(l+1)) l = l+1
           !write(*,*) r_grain(k), l, a_sph(l), a_sph(l+1)
           !fact = (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))
           fact = (log(r_grain(k))-log(a_sph(l)))/(log(a_sph(l+1))-log(a_sph(l)))
           do i=1, n_rad
              do j=1,nz
                 densite_pouss(i,j,1,k) = exp(log(rho(i,j,l))  + fact* (log(rho(i,j,l+1)) - log(rho(i,j,l))) )
              enddo
           enddo
        endif
     enddo
  else !lstrat
     do k=1,n_grains_tot
        do i=1, n_rad
           do j=1,nz
              densite_pouss(i,j,1,k) = rho(i,j,1)
           enddo
        enddo
     enddo
  endif

  ! Normalisation : on a 1 grain de chaque taille dans le disque
  somme=0.0
  do k=1,n_grains_tot
     do i=1,n_rad
        do j=1,nz
           if (densite_pouss(i,j,1,k) < 0.0) then
              write(*,*) "Pb densite", densite_pouss(i,j,1,k)
              write(*,*) i,j,k
              read(*,*)
           endif
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo


  ! Normalisation : Calcul masse totale
  masse = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,k) * (dust_pop(grain(k)%pop)%rho1g_avg * 4.*pi/3. * &
                (r_grain(k)*1.e-4)**3) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/sum(masse) * Msun_to_g
  masse =  masse* diskmass/sum(masse) * Msun_to_g


 ! Verif
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass
  write(*,*) "SPH data OK"
  write(*,*) "masse totale =", mass/ Msun_to_g


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_data_SPH_TTauri_2

!**********************************************************************

subroutine densite_data_LAURE_SED()

! redefinit la grille pour suivre celle de Laure
! Calcule les tableaux zmax, volume, r_lim, r_lim_2, z_lim et delta0

! Calcule la table de densite
! Inclus stratification empirique
! Calcule les tableaux densite_pouss et masse
! + densite_gaz et masse_gaz
! et indique ri_not_empty, zj_not_empty et phik_not_empty

  integer, parameter :: nr = 340
  integer, parameter :: nvert = 50
  real, dimension(nr,nvert) :: rho, T

  real :: TT, tautau, z_r ! buffer
  real ::  M

  integer :: i, j, k, l, pop, sys_status
  real(kind=db) :: facteur, mass

  character(len=512) :: cmd

  type(dust_pop_type), pointer :: dp

  ! Tests dimensions
  if ((n_rad/=nr).or.(nz/=nvert).or.(n_rad_in/=1)) then
     write(*,*) "ERROR in spatial grid"
     write(*,*) "n_rad, nz and n_rad_in must be 340,",nvert," and 1"
     write(*,*) "Exiting"
     stop
  endif

  write(*,*) "*********************************"
  write(*,*) "*    reading Laure's files      *"
  write(*,*) "*  Updating grid and density    *"
  write(*,*) "*********************************"

  allocate(Temperature_Laure_SED(n_rad,nz,1))


  ! Lecture du fichier
  open(unit=1,file=trim(Laure_SED_filename),status="old")
  read(1,*)
  read(1,*)
  do i=1, n_rad
     do j=1, nz
        read(1,*) r_grid(i,j), z_r, T(i,j), rho(i,j), tautau
        z_grid(i,j) = r_grid(i,j) * z_r
        Temperature_Laure_SED(i,j,1) = T(i,j)
     enddo

     zmax(i) = z_grid(i,nz) + 0.5 * (z_grid(i,nz) - z_grid(i,nz-1))
     delta_z(i)=zmax(i)/real(nz)
  enddo !i
  close(unit=1)
  zmaxmax = maxval(zmax)

  call cfitsWrite("data_th/r_grid.fits", real(r_grid), shape(r_grid))
  call cfitsWrite("data_th/z_grid.fits", real(z_grid), shape(z_grid))
  call cfitsWrite("data_th/T_Laure.fits", real(T), shape(T))

  ! limite superieure des cellules
  do i=1, n_rad-1
     r_lim(i) = sqrt(r_grid(i,1) * r_grid(i+1,1))
  enddo
  r_lim(0) = r_grid(1,1) / sqrt(r_grid(2,1)/r_grid(1,1))
  r_lim(n_rad) = r_grid(n_rad,1) * sqrt(r_grid(n_rad,1)/r_grid(n_rad-1,1))
  rmax = r_lim(n_rad)
  rmin = r_lim(0)


  r_lim_2(:) = r_lim(:) * r_lim(:)
  r_lim_3(:) = r_lim_2(:) * r_lim(:)

  do i=1, n_rad
     delta_z(i)=zmax(i)/real(nz)
     z_lim(i,nz+1)=zmax(i)
     do j=1,nz
        z_lim(i,j) = (real(j,kind=db)-1.0_db)*delta_z(i)
     enddo
  enddo

  ! Volume
  do i=1,n_rad
     volume(i)=2.0_db*pi*(r_lim_2(i)-r_lim_2(i-1)) * zmax(i)/real(nz)
  enddo

  ! Masse totale
  M = 0
  do i=1, n_rad
     do j=1, nz
        M = M + Rho(i,j) * (volume(i) * AU_to_cm**3) / Msun_to_g
     enddo
  enddo
  !write(*,*) "M", M


  ! Table de densite
  do i=1, n_rad
     do j=1, nz
        do l=1,n_grains_tot
           densite_pouss(i,j,1,l) = Rho(i,j) * nbre_grains(l)
        enddo ! k
        !if (i==1) write(*,*) j, sum(densite_pouss(i,j,1,:))
     enddo !j
  enddo !i

  ! Normalisation : Calcul masse totale
  do pop=1, n_pop
     dp => dust_pop(pop)
     mass = 0.0
     do i=1,n_rad
        bz2 : do j=j_start,nz
           if (j==0) cycle bz2

           do k=1,n_az
              do l=dp%ind_debut,dp%ind_fin
                 mass=mass + densite_pouss(i,j,k,l) * M_grain(l) * (volume(i) * AU3_to_cm3)
              enddo !l
           enddo !k
        enddo bz2
     enddo !i
     mass =  mass*g_to_Msun

     facteur = M / mass * 0.01 ! masse de poussiere

     masse = 0.0
     do i=1,n_rad
        bz3 : do j=j_start,nz
           if (j==0) cycle bz3
           do k=1, n_az
              do l=dp%ind_debut,dp%ind_fin
                 densite_pouss(i,j,k,l) = densite_pouss(i,j,k,l) * facteur
                 masse(i,j,k) = masse(i,j,k) + densite_pouss(i,j,k,l) * M_grain(l) * volume(i)
              enddo !l
           enddo !k
        enddo bz3
     enddo ! i

  enddo ! pop

  masse(:,:,:) = masse(:,:,:) * AU3_to_cm3

  write(*,*) 'Total dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'

  cmd = "rm -rf T_Laure.fits.gz"
  call appel_syst(cmd, sys_status)
  call cfitsWrite("T_Laure.fits.gz",T,shape(T))

  cmd = "rm -rf Dens_Laure.fits.gz"
  call appel_syst(cmd, sys_status)
  call cfitsWrite("Dens_Laure.fits.gz",Rho,shape(Rho))


  ri_not_empty = 1
  zj_not_empty = 1
  phik_not_empty = 1

  if (lcylindrical) densite_pouss(:,nz+1,:,:) = 0. !densite_pouss(:,nz,:,:)

  return

end subroutine densite_data_LAURE_SED

!**********************************************************************

subroutine force_T_LAURE_SED()

  write(*,*) "Forcing temperature"
  temperature(:,:,:) = temperature_Laure_SED(:,:,:)

  return

end subroutine force_T_LAURE_SED

!**********************************************************************

subroutine densite_eqdiff()

  implicit none

  real, parameter :: G = 6.672e-8
  real, parameter :: gas_dust = 100

  integer :: i,j, k, l, k_min, jj
  real :: cst, cst_pous, cst_gaz, M_star, a
  real :: fact_exp, density,  c_sound, proba, coeff_grav, omega, D0, eps, pas_z, somme1, somme2, correct,test

  real, dimension(1) :: y

  real, dimension(nz) :: z

  real, dimension(n_grains_tot, nz) :: correct_strat
  real, dimension(nz) :: rho

  type(disk_zone_type) :: dz

  if (n_zones > 1) then
     write(*,*) "Error : n_zones must be set to 1 when using densite_eqdiff"
     stop
  endif

  dz=disk_zone(1)

  ! Cste densite disque
  if (abs(dz%surf+2.0) > 1.0e-5) then
     cst=(((2.0+dz%surf)*dz%diskmass)/(((2*pi)**1.5)*dz%sclht*(dz%rin**2)*&
          ((dz%rmax/dz%rin)**(2+dz%surf)-1)))*((dz%rref/dz%rin)**dz%surf)
  else
     cst=diskmass/((2*pi)**1.5*dz%sclht*(dz%rref**2)) * log(dz%rin/dz%rmax)
  endif
  ! facteur multiplicatif pour passer en g/cm**3:
  cst=cst * Msun_to_g / AU3_to_cm3

  ! facteur multiplicatif pour passer en part/cm**3: AVG_GRAIN_MASS
  cst_pous = cst/dust_pop(1)%avg_grain_mass

  ! Facteur multiplicatif pour passer en masse de gaz
  cst_gaz = cst*gas_dust


!******************************************************************************
! Unites
!******************************************************************************
  ! Masse etoile en gramme
  M_star = sum(etoile(:)%M)*Msun_to_g
  coeff_grav = sum(etoile(:)%M) * Msun_to_g/AU3_to_cm3


!***************
! Calcul proba cumulee en dessous d'une taille de grains
! et opacite pour chaque position
!
!* probsizecumul(i) represente la probabilite cumulee en-dessous d'une
!* certaine taille de grain. Ce tableau est utilise pour le tirage
!* aleatoire de la taille du grain diffuseur, puisqu'elle doit prendre
!* en compte le nombre de grains en meme temps que leur probabilite
!* individuelle de diffuser (donnee par qsca*pi*a**2).

  do i=1, n_rad
!     rcyl = i*1.0/real(resol)
     ! On calcule la densite au milieu de la cellule
     rcyl = sqrt((r_lim(i) * r_lim(i-1)))
     fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta)
     coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))
     c_sound = v_sound * 100. * (rcyl/rref_v)**(dz%exp_beta - 3./2.) ! 100 -> cm
     omega =  sqrt(coeff_grav/(rcyl)**3)
     D0 = alpha*c_sound*dz%sclht*(rcyl/dz%rref)**dz%exp_beta*1.5e13 !1.5e13 = schlt en cm

     do j=1,nz
        z(j) = delta_z(i) * (real(j)-0.5)
        if (rcyl < dz%rin) then
           rho(j) = cst_pous*fact_exp * exp(-((z(j)/dz%sclht)**2)/(coeff_exp)) * exp(-((rcyl-dz%rin)**2)/(2*dz%edge**2))
        else
           rho(j) = cst_pous*fact_exp *  exp(-((z(j)/dz%sclht)**2)/(coeff_exp))
        endif
     enddo

     ! Densite de colonne (integ + analytique pour verifier si nz suffisant)
     ! Verifie Ok par rapport a Dullemond et al 2004
     if (i==140) then
        test = 0.0
        do j=1,nz
           test = test + rho(j)*delta_z(i)*1.5e13
        enddo
        write(*,*) 'Colonne densite (g/cm^2) a r=', rcyl
        write(*,*) 'Integration / Analytique'
        write(*,*) 2.*test*cst_gaz/cst_pous, fact_exp*cst_gaz*sqrt(coeff_exp*pi)*dz%sclht*1.5e13
     endif

    ! Resolution eq_diff pour stratification
     if (lstrat) then
        eps = 1.e-2
        pas_z = delta_z(i)
        do k=1, n_grains_tot
           ! Coeff dependant de a
           coeff1 = (omega**2*dust_pop(1)%rho1g_avg*r_grain(k)*mum_to_cm)/(c_sound*D0*cst_gaz*fact_exp)
           coeff2 = (omega*dust_pop(1)%rho1g_avg*r_grain(k)*mum_to_cm)/(c_sound**cst_gaz*fact_exp)
           ! Integration main loop
           y(1) = 1.0
           correct_strat(k,1) = 1.0
           vertical : do j=1,nz-1
              call odeint(y,z(j),z(j+1),eps,pas_z,0.,derivs,rkqs)
              if (y(1) < 1.e-15) then
                 ! On met tout ce qu'il y a au dessus a 0 et on ne reintegre pas
                 ! influe le tps de calcul et bug si precision trop grande
                 do jj = j,nz-1
                    correct_strat(k,jj+1) = 0.0
                 enddo !jj
                 exit vertical
!                 y(1) = 0.0
              endif
              correct_strat(k,j+1) = y(1)
           enddo vertical !j
           ! Conservation de la masse
           somme1=0.0
           somme2=0.0
           do j=1,nz
              somme1 = somme1 + rho(j)
              somme2 = somme2 + rho(j)* correct_strat(k,j)
           enddo
           correct = somme1/somme2
           do j=1,nz
              correct_strat(k,j) = correct_strat(k,j)*correct
           enddo
!           if ((rcyl > 10.0).and.(r_grain(k)>0.1)) then
!              do j=1,nz
!                 write(*,*) z(j)/rcyl,  correct_strat(k,j)
!              enddo
!              stop
!           endif
        enddo !k
     else!lstrat
        correct_strat = 1.0
     endif!lstrat
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     ! Calcul opacite et probabilite de diffusion
     !$omp parallel &
     !$omp default(none) &
     !$omp shared(i,rcyl,nz,n_grains_tot) &
     !$omp shared(amax_reel,densite_pouss) &
     !$omp private(j,k,z,density,k_min,proba)&
     !$omp shared(zmax,kappa,probsizecumul,ech_prob,nbre_grains,cst_pous,q_ext,q_sca,rho,correct_strat, delta_z)
     !$omp do schedule(dynamic,10)
     do j=1,nz
        do  k=1,n_grains_tot
           densite_pouss(i,j,1,k) = nbre_grains(k) * correct_strat(k,j) * rho(j)
        enddo !k
     enddo !j
     !$omp enddo
     !$omp end parallel
  enddo !i


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_eqdiff

!***********************************************************

subroutine derivs(x,y,dydx)

  implicit none

  REAL(SP), INTENT(IN) :: x
  REAL(SP), DIMENSION(:), INTENT(IN) :: y
  REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
  real :: rho_gaz, schlt

  type(disk_zone_type) :: dz

  dz= disk_zone(1)

  ! x en AU
  if (rcyl < dz%rin) then
     rho_gaz = exp(-((x/dz%sclht)**2)/(coeff_exp)) * exp(-((rcyl-dz%rin)**2)/(2*dz%edge**2))
  else
     rho_gaz = exp(-((x/dz%sclht)**2)/(coeff_exp))
  endif

  ! (1.5e13)^2 = 2.25e26 -> 1 AU pour x et 1 AU pour dx
  dydx = -y * coeff1 * x/rho_gaz * (1.0 + coeff2/rho_gaz) *  2.25e26

end subroutine derivs

!***********************************************************

subroutine densite_data_hd32297()
! Lecture table densite JC
! 3/05/05

  implicit none


!  integer, parameter :: nr=1150
  integer :: i, j, k, kmin, longueur, status, nr, alloc_status
  real, dimension(:), allocatable :: rayon, dens
  real, dimension(n_rad) :: density
  real, dimension(nz) ::  z
  real :: somme, mass, rcyl, coeff_exp, a, b, gamma

  character(len=62) :: buf

  type(disk_zone_type) :: dz

  dz=disk_zone(1)

  longueur = len_trim(para)

  open(unit=1,file=trim(para(1:longueur-4))//"txt",status='old')

   ! On compte d'abord les lignes
  status=0
  nr=0
  do while(status==0)
     nr=nr+1
     read(1,*,IOSTAT=status)
  enddo
  nr = nr-1

  ! on va lire 8 lignes a part
  nr = nr - 8

  ! Ollocation des tableaux
  allocate(rayon(nr), dens(nr), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  rayon=0.0
  dens=0.0

  ! Lecture
  rewind(1)

  read(1,*)
  read(1,*) buf,buf, gamma
  read(1,*)  buf, buf, dz%exp_beta
  read(1,*)  buf, buf, dz%sclht
  write(*,*) gamma, dz%exp_beta, dz%sclht
  do i=1,4
     read(1,*)
  enddo

  ! On peut redefinir la grille
  call define_grid3

  do i=1,nr
     read(1,*) rayon(i), dens(i)!, a, b
     if (dens(i) < tiny_real) dens(i) = 1.0e-20
  enddo
  close(unit=1)

  kmin=1

  do i=1,n_rad
     rcyl=0.5*(r_lim(i-1)+r_lim(i))
     !     coeff_exp = (2*(rcyl/rref)**(2*exp_beta))
     coeff_exp = ((rcyl/dz%rref)**(2*dz%exp_beta)) ! Gaussienne definie sans le 1/2
     ! Interpolation densite midplane
     inter : do k=kmin, nr
        if (rayon(k) > rcyl) then
           kmin=k-1
           exit inter
        endif
     enddo inter
     density(i) = dens(kmin) + (dens(kmin+1)-dens(kmin)) * (rcyl-rayon(kmin))/(rayon(kmin+1)-rayon(kmin))
   !  write(*,*) i,  rayon(kmin), rcyl,  rayon(kmin+1), dens(kmin),  density(i), dens(kmin+1)


     do j=1,nz
     ! On calcule la densite au milieu de la cellule
        z(j) = (real(j)-0.5)*delta_z(i)
        do  k=1,n_grains_tot
           densite_pouss(i,j,1,k) = density(i) * nbre_grains(k) * exp(-((z(j)/dz%sclht)**gamma)/(coeff_exp))
        enddo !k
     enddo !j
  enddo


!***************
! Remplissage a zero pour z > zmax
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do k=1,n_grains_tot
     somme=0.0
     do i=1,n_rad
        do j=1,nz
           if (densite_pouss(i,j,1,k) < 0.0) then
              write(*,*) i,j,k, densite_pouss(i,j,1,k)
              read(*,*)
           endif
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo

  ! Normalisation : Calcul masse totale
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  return

end subroutine densite_data_hd32297

!***********************************************************

subroutine densite_data_gap()
! Lecture table desnite Peggy
! 5/01/06

  implicit none

  integer :: i, j, k, kmin, longueur, status, nr, alloc_status
  real, dimension(:), allocatable :: rayon, dens
  real, dimension(n_rad) :: density
  real, dimension(nz) ::  z
  real :: somme, mass, rcyl, coeff_exp, a, b, gamma

  character(len=62) :: buf

  integer :: ind_debut, ind_fin
  real :: r_debut, r_fin

  type(disk_zone_type) :: dz

  dz = disk_zone(1)

  open(unit=1,file="gasdens4000.txt",status='old')

  ! On compte d'abord les lignes
  status=0
  nr=0
  do while(status==0)
     nr=nr+1
     read(1,*,IOSTAT=status)
  enddo
  nr = nr-1

  ! Allocation des tableaux
  allocate(rayon(nr), dens(nr), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error crea_noyau'
     stop
  endif
  rayon=0.0
  dens=0.0

  ! Lecture
  rewind(1)

  do i=1,nr
     read(1,*) rayon(i), dens(i)
     if (dens(i) < tiny_real) dens(i) = 1.0e-20
  enddo
  close(unit=1)

  kmin=1

  ind_fin = 100000

  ! Interpolation densite
  do i=1,n_rad
     rcyl=0.5*(r_lim(i-1)+r_lim(i))
     if (rcyl < rayon(1)) then
        ind_debut=i
     elseif (rcyl > rayon(nr)) then
        if (ind_fin > i) ind_fin=i
     else
        ! Interpolation densite midplane
        inter : do k=kmin, nr
           if (rayon(k) > rcyl) then
              kmin=k-1
              exit inter
           endif
        enddo inter
        density(i) = dens(kmin) + (dens(kmin+1)-dens(kmin)) * (rcyl-rayon(kmin))/(rayon(kmin+1)-rayon(kmin))
     end if
  enddo

  ! Loi de puissance pour region interne
  r_debut=0.5*(r_lim(ind_debut)+r_lim(ind_debut+1))
  do i=1,ind_debut
     rcyl=0.5*(r_lim(i-1)+r_lim(i))
     density(i) =  density(ind_debut+1) * (rcyl/r_debut)**dz%surf
  enddo

  ! Loi de puissance pour region externe
  r_fin=0.5*(r_lim(ind_fin-2)+r_lim(ind_fin-1))
  do i=ind_fin,n_rad
     rcyl=0.5*(r_lim(i-1)+r_lim(i))
     density(i) =  density(ind_fin-1) * (rcyl/r_fin)**dz%surf
  enddo

  write(*,*) density(ind_fin-1)

  ! Profil vertical
  do i=1,n_rad
     rcyl=0.5*(r_lim(i-1)+r_lim(i))
     coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))
     do j=1,nz
     ! On calcule la densite au milieu de la cellule
        z(j) = (real(j)-0.5)*delta_z(i)
        do  k=1,n_grains_tot
           densite_pouss(i,j,1,k) = density(i) * nbre_grains(k) * exp(-((z(j)/dz%sclht)**2)/(coeff_exp))
        enddo !k
     enddo !j
  enddo


!***************
! Remplissage a zero pour z > zmax
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do k=1,n_grains_tot
     somme=0.0
     do i=1,n_rad
        do j=1,nz
           if (densite_pouss(i,j,1,k) < 0.0) then
              write(*,*) i,j,k, densite_pouss(i,j,1,k)
              read(*,*)
           endif
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo

  ! Normalisation : Calcul masse totale
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  return

end subroutine densite_data_gap

!***********************************************************

subroutine init_opacity_wall()
  ! Calcule la densite du mur pour avoir la bonne opacite
  ! ATTENTION : ca n'est qu'un mur en opacite pas en emission !!!
  ! C. Pinte
  ! 1/05/07 et oui, c'est pas ferie en Angleterre !!! :-(

  implicit none

  ! Opacite du mur qui a pour largeur celle de la premiere cellule
  kappa_wall = tau_wall / (r_lim(1) - r_lim(0))

  ! On ne vire pas les photons au-dessus de l'echelle de hauteur du disque
  ! Au cas ou ils taperrait dans le mur
  cos_max2 = 0.0

  write (*,*) 'Masse mur ='!,8*delta_r*r_wall*h_wall*rho_wall*avg_grain_mass/0.594098e-6

  return

end subroutine init_opacity_wall

!********************************************************************

subroutine densite_gap_laure()
  ! Remplace par densite_gap_laure2

  implicit none

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j, hdunum, hdutype
  integer :: nullval
  integer, dimension(4) :: naxes
  logical :: anynull

  integer, parameter :: nbre_a = 4

  real, dimension(4) :: a_sph
  integer :: k, l, i
  real(kind=db) :: somme, mass

  real, dimension(n_rad,nz,nbre_a) :: sph_dens


  sph_dens = 1.

  ! Lecture donnees
  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  write(*,*) "Reading density file : "//trim(density_file)

  readwrite=0
  call ftopen(unit,density_file,readwrite,blocksize,status)
  if (status /= 0) then ! le fichier temperature n'existe pas
     write(*,*) "ERROR : density file needed"
     stop
  endif

  group=1
  firstpix=1
  nullval=-999

  !  determine the size of density file
  call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
  if (nfound /= 3) then
     write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
     write(*,*) 'of '//trim(density_file)//' file. Exiting.'
     stop
  endif

  if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
     write(*,*) "Error : "//trim(density_file)//" does not have the"
     write(*,*) "right dimensions. Exiting."
     write(*,*) naxes(1), n_rad
     write(*,*) naxes(2), nz
     stop
  endif

  npixels=naxes(1)*naxes(2)*naxes(3)

  nbuffer=npixels
  ! read_image
  call ftgpve(unit,group,firstpix,nbuffer,nullval,sph_dens,anynull,status)

  call ftclos(unit, status)
  call ftfiou(unit, status)

  ! Interpolation en tailles
  a_sph(1) = 10. ! densite des grains <= 10 microns == celle du gas
  a_sph(2) = 100.
  a_sph(3) = 1000.
  a_sph(4) = 10000.


  if (lstrat) then
     write(*,*) "Differential gap"
     l=1
     do k=1,n_grains_tot
        if (r_grain(k) < a_sph(1)) then  ! Petits grains
           densite_pouss(:,1:nz,1,k) = sph_dens(:,:,1)
        else if (r_grain(k) > a_sph(nbre_a)) then ! Gros grains
           densite_pouss(:,1:nz,1,k) = sph_dens(:,:,nbre_a)
        else  ! Autres grains : interpolation
           if (r_grain(k) > a_sph(l+1)) l = l+1
           densite_pouss(:,1:nz,1,k) = sph_dens(:,:,l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l)) * &
                ( sph_dens(:,:,l+1) -  sph_dens(:,:,l) )
        endif
     enddo
  else ! Tous les grains suivent le gas
     write(*,*) "Constant Gap"
     do k=1,n_grains_tot
        densite_pouss(:,1:nz,1,k) = sph_dens(:,:,1)
     enddo
  endif  !lstrat

  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do k=1,n_grains_tot
     somme=0.0
     do i=1,n_rad
        do j=1,nz
           if (densite_pouss(i,j,1,k) < 0.0) densite_pouss(i,j,1,k) = 1.0e-20
           somme=somme+densite_pouss(i,j,1,k)*volume(i)
        enddo
     enddo
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)
  enddo

  ! Normalisation : on a 1 grain en tout dans le disque
  do k=1,n_grains_tot
     densite_pouss(:,:,1,k) = (densite_pouss(:,:,1,k)/somme)*nbre_grains(k)
  enddo

  ! Normalisation : Calcul masse totale
  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  write(*,*) "Done"

  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,k) * M_grain(k) * volume(i)
        enddo !k
     enddo !j
  enddo ! i

  masse(:,:,:) = masse(:,:,:) * AU3_to_cm3

  write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'

  return

end subroutine densite_gap_laure

!**********************************************************

subroutine densite_gap_laure2()
  ! Nouvelle routine pour lire les grilles de densite
  ! calculees par Yorick directement a partir des donnees SPH
  ! Les donnees sont directement lissees sur la grille de MCFOST
  ! C. Pinte
  ! 12/0/09

  implicit none

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j, hdunum, hdutype
  integer :: nullval, stat
  integer, dimension(4) :: naxes
  logical :: anynull
  character(len=80) :: comment

  integer :: k, l, i, n_a
  real(kind=db) :: somme, mass
  real :: tmp
  character :: s

  real, dimension(:,:,:,:), allocatable :: sph_dens ! (n_rad,nz,n_az,n_a)
  real, dimension(:), allocatable :: a_sph ! n_a


  sph_dens = 1.

  ! Lecture donnees
  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  write(*,*) "Reading density file : "//trim(density_file)

  readwrite=0
  call ftopen(unit,density_file,readwrite,blocksize,status)
  if (status /= 0) then ! le fichier temperature n'existe pas
     write(*,*) "ERROR : density file needed"
     stop
  endif

  group=1
  firstpix=1
  nullval=-999

  !  determine the size of temperature file
  call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
  if (nfound /= 4) then
     write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
     write(*,*) 'of '//trim(density_file)//' file. Exiting.'
     stop
  endif

  if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= n_az) ) then
     write(*,*) "Error : "//trim(density_file)//" does not have the"
     write(*,*) "right dimensions. Exiting."
     write(*,*) "# fits_file vs mcfost_grid"
     write(*,*) naxes(1), n_rad
     write(*,*) naxes(2), nz
     write(*,*) naxes(3), n_az
     !write(*,*) naxes(4), n_a
     stop
  endif
  n_a = naxes(4)
  write(*,*) n_a, "grain sizes found"
  npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
  nbuffer=npixels

  allocate(sph_dens(n_rad,nz,n_az,n_a), a_sph(n_a))

  ! read_image
  call ftgpve(unit,group,firstpix,nbuffer,nullval,sph_dens,anynull,status)

  ! Au cas ou
  sph_dens = sph_dens + 1e-20

  ! Lecture des tailles de grains (en microns)
  if (n_a > 9) then
     write(*,*) "ERROR : max 9 grain sizes at the moment"
     write(*,*) "code must be updated if you need more"
     write(*,*) "Exiting."
     stop
  endif
  do i=1,n_a
     write(s,'(i1)') i ; call ftgkye(unit,'grain_size_'//s,tmp,comment,stat)
     a_sph(i) = tmp ! cannot read directly into an array element
     write(*,*) i, a_sph(i), "microns"
  enddo

  ! On verifie que les grains sont tries
  do i=1, n_a-1
     if (a_sph(i) > a_sph(i+1)) then
        write(*,*) "ERROR : grains must be ordered from small to large"
        write(*,*) "Exiting"
        stop
     endif
  enddo

  call ftclos(unit, status)
  call ftfiou(unit, status)


  if (lstrat) then
     write(*,*) "Differential gap"
     l=1
     do k=1,n_grains_tot
        if (r_grain(k) < a_sph(1)) then  ! Petits grains
           densite_pouss(:,1:nz,:,k) = sph_dens(:,:,:,1)
           !write(*,*) k, r_grain(k), "p"
        else if (r_grain(k) > a_sph(n_a)) then ! Gros grains
           densite_pouss(:,1:nz,:,k) = sph_dens(:,:,:,n_a)
           !write(*,*) k, r_grain(k), "g"
        else  ! Autres grains : interpolation
           !write(*,*) k, r_grain(k), "i"
           if (r_grain(k) > a_sph(l+1)) l = l+1
           densite_pouss(:,1:nz,:,k) = sph_dens(:,:,:,l) + (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l)) * &
                ( sph_dens(:,:,:,l+1) -  sph_dens(:,:,:,l) )
        endif
     enddo
  else ! Tous les grains suivent le gas
     write(*,*) "Constant Gap"
     do k=1,n_grains_tot
        densite_pouss(:,1:nz,:,k) = sph_dens(:,:,:,1)
     enddo
  endif  !lstrat

  ! Symetrie verticale en z
  do j=1,nz
     densite_pouss(:,-j,:,:) = densite_pouss(:,j,:,:)
  enddo

  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do l=1,n_grains_tot
     somme=0.0
     do i=1,n_rad
        do j=-nz,nz
           if (j==0) cycle
           do k=1,n_az
              if (densite_pouss(i,j,k,l) <= 0.0) densite_pouss(i,j,k,l) = 1.0e-20
              somme=somme+densite_pouss(i,j,k,l)*volume(i)
           enddo !k
        enddo !j
     enddo !i
     densite_pouss(:,:,:,l) = (densite_pouss(:,:,:,l)/somme)
  enddo !l

  ! Normalisation : on a 1 grain en tout dans le disque
  do l=1,n_grains_tot
     densite_pouss(:,:,:,l) = (densite_pouss(:,:,:,l)/somme)*nbre_grains(l)
  enddo

  search_not_empty : do l=1,n_grains_tot
     do k=1,n_az
        do j=1,nz
           do i=1,n_rad
              if (densite_pouss(i,j,k,l) > 0.0_db) then
                 ri_not_empty = i
                 zj_not_empty = j
                 phik_not_empty = k
                 exit search_not_empty
              endif
           enddo !i
        enddo !j
     enddo !k
  enddo search_not_empty


  ! Normalisation : Calcul masse totale
  mass = 0.0
  do i=1,n_rad
     do j=-nz,nz
        if (j==0) cycle
        do k=1,n_az
           do l=1,n_grains_tot
              mass=mass + densite_pouss(i,j,k,l) * M_grain(l) * (volume(i) * AU3_to_cm3)
           enddo !l
        enddo !k
     enddo !j
  enddo !i
  mass =  mass/Msun_to_g
  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  write(*,*) "Done"

  do i=1,n_rad
     !write(*,*) i, r_grid(i,1), sum(densite_pouss(i,:,:,3))
     do j=-nz,nz
        if (j==0) cycle
        do k=1,n_az
           do l=1,n_grains_tot
              masse(i,j,k) = masse(i,j,k) + densite_pouss(i,j,k,l) * M_grain(l) * volume(i)
           enddo !l
        enddo !k
     enddo !j
  enddo ! i

  masse(:,:,:) = masse(:,:,:) * AU3_to_cm3

  write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'
  deallocate(sph_dens,a_sph)

  return

end subroutine densite_gap_laure2

!**********************************************************

subroutine densite_debris
! Lecture du 1er jeu de simulation de Holly Manness
! G.D., 29 Fevrier 2008

  implicit none

  integer :: i, j, k, l

  real, dimension(n_rad,nz,n_grains_tot) :: rho

  real(kind=db) :: mass

  ! Lecture données
  write(*,*) "Reading debris disk dynamical structure file... "
  open(unit=1,file=debris_file,status='old')
  do l=1,37
     read(1,*)
  enddo
  do i=1, n_rad
     do j=1, nz
        do l=1, n_grains_tot
           read(1,*) rho(i,j,l)
        enddo
     enddo
  enddo
  close(unit=1)
  write(*,*) " ... done!"

  ! Populating grid
  do k=1,n_grains_tot
     do i=1, n_rad
        do j=1,nz
           densite_pouss(i,j,1,k) = rho(i,j,k) !/ (volume(i) * AU3_to_cm3)
        enddo
     enddo
!     write(6,*) r_grain(k),tab_densite(1,1,1,k)
  enddo

  ! Normalisation : Calcul masse totale
  masse = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass=sum(masse)

  write(*,*) "Debris disk data OK"
  write(*,*) "masse totale =", mass/ Msun_to_g


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return
end subroutine densite_debris

!**********************************************************************

subroutine densite_fits
! Lecture du 2e jeu de simulation de Holly Manness, using a FITS input file
! Input file is assumed to be N(r,z,a) in units of cm^-3
! G.D., 28 May 2009

  implicit none

  integer :: i, j, k, l, unit, status, readwrite, blocksize
  integer :: npixels, group, firstpix, nullval, nbuffer, alloc_status
  logical ::  anynull

  real, allocatable, dimension(:,:,:,:) :: rho_tmp

  type(dust_pop_type), pointer :: dp

  real(kind=db) :: mass

  ! Lecture données
  write(*,'(a49, $)') "Reading input FITS file for disk structure ... "
  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  readwrite=0
  call ftopen(unit,struct_fits_file,readwrite,blocksize,status)
  if (status /= 0)then
     write(*,*) 'ERROR opening FITS file '//trim(struct_fits_file)
     stop
  endif

  dp => dust_pop(1)

  allocate(rho_tmp(n_rad,nz,dp%n_grains,n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error rho_tmp'
     stop
  endif

  !  initialize variables.
  npixels=n_rad*nz*n_grains_tot
  group=1
  firstpix=1
  nullval=-999
  nbuffer=npixels

  ! read input file
  call ftgpve(unit,group,firstpix,nbuffer,nullval,rho_tmp,anynull,status)

  do i=1, n_rad
     do j=1, nz
        do k=1, dp%n_grains
           do l=0, n_pop-1
              densite_pouss(i,j,1,k+l*dp%n_grains)=rho_tmp(i,j,k,l+1)
           enddo
        enddo
     enddo
  enddo

  ! closing input file and freeing unit number and temporary memory
  call ftclos(unit, status)
  call ftfiou(unit, status)
  deallocate(rho_tmp)

  write(*,*) " Done"

  ! Normalisation : Calcul masse totale
  masse = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,k) * M_grain(l) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass=sum(masse)

  write(*,*) "Total dust mass =", mass/ Msun_to_g," Msun"


!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(i,nz+1,1,k) = densite_pouss(i,nz,1,k)
     enddo !k
  enddo !i

  return

end subroutine densite_fits

!**********************************************************************

real(kind=db) function omega_tau(rho,H,l)
  ! Pour les calculs de Sebastien Fromang
  ! rho doit etre en g.cm-3

  real(kind=db), intent(in) :: rho, H
  integer, intent(in) :: l

  integer :: ipop

  ipop = grain(l)%pop
  omega_tau= dust_pop(ipop)%rho1g_avg*(r_grain(l)*mum_to_cm) / (rho * masse_mol_gaz/m_to_cm**3 * H*AU_to_cm)

  return

end function omega_tau

!**********************************************************************

subroutine densite_Seb_Charnoz()

  integer :: Nr_Seb, Nz_Seb, Na_Seb
  real(kind=db), dimension(n_grains_tot) :: density_Seb, taille_grains_Seb, N_grains
  real(kind=db) :: Rmin, Dr, Zmin, Dz
  real(kind=db) :: Rmin_mcfost, Dr_mcfost, Zmin_mcfost, Dz_mcfost

  real(kind=db) :: Somme
  integer :: ii, jj, i, j, k, l

  write(*,*) "***********************************************"
  write(*,*) "Reading results from Sebastien Charnoz ..."
  open(unit=1,file="/Volumes/home/cpinte/Seb_C/twhydra_simuturb_mcfost.dat",status="old")
  read(1,*)
  read(1,*) Nr_Seb, Nz_Seb, Na_Seb

  if ((Nr_Seb /= n_rad).or.(Nz_Seb /= nz)) then
     write(*,*) "ERROR: Spatial grid does not match!"
     write(*,*) "Exiting"
     stop
  endif

  if (Na_Seb /= n_grains_tot) then
     write(*,*) "ERROR: Grain size grid does not match!"
     write(*,*) "Exiting"
     stop
  endif

  read(1,*)
  read(1,*) taille_grains_Seb
  r_grain(:) = taille_grains_Seb(:) * m_to_mum ! conversion en micron

  read(1,*)
  Somme = 0.
  do i=1,n_rad
     do j=1, nz
        Rmin_mcfost = r_lim(i-1)
        Dr_mcfost = r_lim(i) - r_lim(i-1)
        Zmin_mcfost = z_lim(i,j)
        Dz_mcfost = z_lim(i,j+1) -  z_lim(i,j)

        read(1,*)  ii, jj, Rmin, Dr, Zmin, Dz, density_Seb(:)

        if (is_diff(Rmin,Rmin_mcfost).and.(j==1)) write(*,*) "Pb R cell", i, Rmin, Rmin_mcfost
        if (is_diff(Dr,Dr_mcfost).and.(j==1)) write(*,*) "Pb Dr cell", i, Dr, Dr_mcfost
        if (is_diff(Zmin,Zmin_mcfost)) write(*,*) "Pb Z cell", i,j, Zmin, Zmin_mcfost
        if (is_diff(Dz,Dz_mcfost))  write(*,*) "Pb Dz cell", i,j


        densite_pouss(i,j,1,:) = density_Seb(:) / (volume(i)*AU3_to_cm3) ! comversion en densite volumique
        Somme = Somme +  1.6 * 4.*pi/3. *  (mum_to_cm)**3 * sum( density_Seb(:) * r_grain(:)**3 )
     enddo ! j
  enddo !i
  write(*,*) "Dust mass from Seb's file :", real(Somme * g_to_Msun), "Msun"

  search_not_empty : do l=1,n_grains_tot
     do j=1,nz
        do i=1,n_rad
           if (densite_pouss(i,j,1,l) > 0.0_db) then
              ri_not_empty = i
              zj_not_empty = j
              phik_not_empty = 1
              exit search_not_empty
           endif
        enddo
     enddo
  enddo search_not_empty


  ! Re-population du tableau de grains
  ! les methodes de chauffages etc, ne changent pas

  do k=1, n_grains_tot
     N_grains(k) = sum(densite_pouss(:,:,1,k))
  enddo
  nbre_grains(:) = N_grains(:)/sum(N_grains(:))


  do i=1,n_rad
     bz : do j=j_start,nz
        if (j==0) cycle bz
        do k=1, n_az
           do l=1,n_grains_tot
              masse(i,j,k) = masse(i,j,k) + densite_pouss(i,j,k,l) * M_grain(l) * volume(i)
           enddo !l
        enddo !k
     enddo bz
  enddo ! i

  masse(:,:,:) = masse(:,:,:) * AU3_to_cm3 * 1600./3500 ! TMP

  write(*,*) 'Total dust mass in model  :', real(sum(masse)*g_to_Msun),'Msun'

  write(*,*) "Done"
  write(*,*) "***********************************************"

end subroutine densite_Seb_Charnoz

!**********************************************************

subroutine densite_Seb_Charnoz2()

  implicit none

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j, hdunum, hdutype
  integer :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull

  integer :: k, l, i
  real(kind=db) :: somme, mass, somme2

  real, dimension(n_rad,nz) :: dens

  dens = 0.

  ! Lecture donnees
  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  write(*,*) "Density structure from Seb. Charnoz" ;
  write(*,*) "Reading density file : "//trim(density_file)

  readwrite=0
  call ftopen(unit,density_file,readwrite,blocksize,status)
  if (status /= 0) then ! le fichier de densite n'existe pas
     write(*,*) "ERROR : density file needed"
     stop
  endif

  group=1
  firstpix=1
  nullval=-999

  !  determine the size of density file
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  if (nfound /= 2) then
     write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
     write(*,*) 'of '//trim(density_file)//' file. Exiting.'
     stop
  endif

  if ((naxes(1) /= n_rad).or.(naxes(2) /= nz) ) then
     write(*,*) "Error : "//trim(density_file)//" does not have the"
     write(*,*) "right dimensions. Exiting."
     write(*,*) "# fits_file vs mcfost_grid"
     write(*,*) naxes(1), n_rad
     write(*,*) naxes(2), nz
     stop
  endif

  npixels=naxes(1)*naxes(2)

  nbuffer=npixels
  ! read_image
  call ftgpve(unit,group,firstpix,nbuffer,nullval,dens,anynull,status)

  call ftclos(unit, status)
  call ftfiou(unit, status)

  ! Au cas ou
  dens = dens + 1e-30

  somme = 0
  somme2 = 0
  do j=1, nz
     !write(*,*) j, z_lim(1,j+1)-z_lim(1,j)
     somme2 = somme2 + dens(1,j) * (z_lim(1,j+1)-z_lim(1,j)) * AU_to_m * 2. ! SI : kg/m^2
  enddo
  write(*,*) "Surface density at 1st radius (g/cm-2) =", somme2 * 1000. / 100**2 ! g/cm-2
  ! END VERIF


  ! Tous les grains suivent le gas
  ! dens est la densite volumique de poussiere en SI. kg/m^3
  do k=1,n_grains_tot
     densite_pouss(:,1:nz,1,k) = dens(:,:)
  enddo

  ! kg/m^3  ---> part/cm^3
  densite_pouss = densite_pouss / ( (cm_to_m)**3  * dust_pop(1)%avg_grain_mass * 1e3)

  ! BUG correction : il manque un facteur
  densite_pouss = densite_pouss * 1e-6


  do l=1,n_grains_tot
     densite_pouss(:,:,:,l) = densite_pouss(:,:,:,l)*nbre_grains(l)
  enddo

!
!  ! Normalisation : on a 1 grain de chaque taille dans le disque
!  do l=1,n_grains_tot
!     somme=0.0
!     do i=1,n_rad
!        do j=1,nz
!           if (densite_pouss(i,j,1,l) <= 0.0) densite_pouss(i,j,1,l) = 1.0e-30
!           somme=somme+densite_pouss(i,j,1,l)*volume(i)
!        enddo !j
!     enddo !i
!     densite_pouss(:,:,:,l) = (densite_pouss(:,:,:,l)/somme)
!  enddo !l
!
!  ! Normalisation : on a 1 grain en tout dans le disque
!  do l=1,n_grains_tot
!     densite_pouss(:,:,:,l) = (densite_pouss(:,:,:,l)/somme)*nbre_grains(l)
!  enddo

  search_not_empty : do l=1,n_grains_tot
     do k=1,n_az
        do j=1,nz
           do i=1,n_rad
              if (densite_pouss(i,j,k,l) > 0.0_db) then
                 ri_not_empty = i
                 zj_not_empty = j
                 phik_not_empty = k
                 exit search_not_empty
              endif
           enddo !i
        enddo !j
     enddo !k
  enddo search_not_empty

!  ! Normalisation : Calcul masse totale
!  mass = 0.0
!  do i=1,n_rad
!     do j=1,nz
!        do l=1,n_grains_tot
!           mass=mass + densite_pouss(i,j,1,l) * M_grain(l) * (volume(i) * AU3_to_cm3)
!        enddo !l
!     enddo !j
!  enddo !i
!  mass =  mass/Msun_to_g
!  densite_pouss(:,:,:,:) = densite_pouss(:,:,:,:) * diskmass/mass

  write(*,*) "Done"

  do i=1,n_rad
     do j=1,nz
        do l=1,n_grains_tot
           masse(i,j,1) = masse(i,j,1) + densite_pouss(i,j,1,l) * M_grain(l) * volume(i)
        enddo !l
     enddo !j
  enddo ! i

  masse(:,:,:) = masse(:,:,:) * AU3_to_cm3

  write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'
  write(*,*) "Density from Seb. Charnoz set up OK"

  return

end subroutine densite_Seb_Charnoz2

!**********************************************************************

subroutine remove_specie()

  implicit none

  integer :: i, j, pk, k
  real :: mass

  write(*,*) "Removing specie", specie_removed, "where T >", T_rm

  do i=1,n_rad
     do j=1,nz
        do pk=1,n_az
           do k=1,n_grains_tot
              if (grain(k)%pop==specie_removed) then
                 if (Temperature(i,j,pk) > T_rm) then
                    densite_pouss(i,j,pk,k) = 0.0
                 endif
              endif
           enddo
        enddo
     enddo
  enddo

  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        do k=1,n_grains_tot
           mass=mass + densite_pouss(i,j,1,k) * M_grain(k) * (volume(i) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g

  write(*,*) 'New total dust mass in model :', mass,' Msun'

  return

end subroutine remove_specie

!************************************************************

end module density
