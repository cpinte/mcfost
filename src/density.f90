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
  use grid
  use utils
  use output

  implicit none

  real :: coeff_exp, coeff1, coeff2, rcyl

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

  integer :: i,j, k, izone, alloc_status, icell
  real(kind=dp), dimension(n_zones) :: cst_gaz
  real(kind=dp) :: z, density, fact_exp, rsph, mass, puffed, facteur, z0, phi, surface, H, C, somme

  type(disk_zone_type) :: dz

  ! Tableau temporaire pour densite gaz dans 1 zone (pour renormaliser zone par zone)
  ! Pas besoin dans la poussiere car a chaque pop, il y a des tailles de grains independantes
  real(kind=dp), dimension(:), allocatable :: densite_gaz_tmp, densite_gaz_midplane_tmp


  allocate(densite_gaz_tmp(n_cells), densite_gaz_midplane_tmp(n_rad), stat=alloc_status)
  densite_gaz_tmp = 0.0 ; densite_gaz_midplane_tmp = 0.0
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

              do k=1, n_az
                 ! On calcule la densite au milieu de la cellule
                 if (j==0) then
                    icell = cell_map(i,1,k)
                    rcyl = r_grid(icell)
                    z = 0.0_dp ! utilisation indice 0 pour densite plan median
                 else
                    icell = cell_map(i,j,k)
                    rcyl = r_grid(icell)
                    z = z_grid(icell)
                 endif
                 phi = phi_grid(icell)

                 H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta

                 ! Puffed-up rim analytique
                 if (lpuffed_rim) then
                    puffed = 1.0 + (puffed_rim_h - 1.0) / (exp((rcyl - puffed_rim_r)/puffed_rim_delta_r) + 1.0)
                 else
                    puffed = 1.0
                 endif

                 if (dz%geometry == 1) then ! power-law
                    fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta)
                 else  if (dz%geometry == 2) then ! tappered-edge : dz%surf correspond a -gamma
                    if (dz%rc < tiny_dp) then
                       write(*,*) "ERROR : tappered-edge structure with Rc = 0."
                       write(*,*) "Exiting"
                       stop
                    endif
                    fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta) * exp( -(rcyl/dz%rc)**(2+dz%moins_gamma_exp) )
                 endif
                 coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))


                 ! Warp analytique
                 if (lwarp) then
                    z0 = z_warp * (rcyl/dz%rref)**3 * cos(phi)
                 else if (ltilt) then
                    if (izone==1) then
                       z0 = rcyl * cos(phi) * tan(tilt_angle * deg_to_rad)
                    else
                       z0 = 0.0
                    endif
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
                 if (j>0) then
                    densite_gaz_tmp(icell) = density
                 else
                    densite_gaz_midplane_tmp(i) = density
                 endif

              enddo !k
           enddo bz !j

           if ((lSigma_file).and.(izone==1)) then
              ! Normalisation pour densite de surface dans fichier
              ! todo : only works for k = 1
              somme = 0.0
              bz2 : do j=min(1,j_start),nz
                 somme = somme + densite_gaz_tmp(cell_map(i,j,1)) *  (z_lim(i,j+1) - z_lim(i,j))
              enddo bz2
              if (somme > tiny_dp) then
                 do j=min(1,j_start),nz
                    densite_gaz_tmp(cell_map(i,j,1)) = densite_gaz_tmp(cell_map(i,j,1)) * Surface_density(i)/somme
                 enddo ! j
                 densite_gaz_midplane_tmp(i) = densite_gaz_midplane_tmp(i) * Surface_density(i)/somme
              endif
           endif
        enddo ! i

     else if (dz%geometry == 3) then ! envelope
        do i=1, n_rad
           bz_env : do j=j_start,nz
              if (j==0) cycle bz_env
              do k=1, n_az
                 if (j==0) then
                    icell = cell_map(i,1,k)
                    rcyl = r_grid(icell)
                    z = 0.0_dp
                 else
                    icell = cell_map(i,j,k)
                    rcyl = r_grid(icell)
                    z = z_grid(icell)
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
              if (j/=0) then
                 densite_gaz_tmp(icell) = density
              else
                 densite_gaz_midplane_tmp(i) = density
              endif
              enddo !k
           enddo bz_env !j
        enddo !i

     else if (dz%geometry == 4) then ! debris
        k=1
        do i=1, n_rad
           bz_debris : do j=j_start,nz
              if (j==0) cycle bz_debris
              do k=1, n_az
                 icell = cell_map(i,j,k)
                 rcyl = r_grid(icell)
                 z = z_grid(icell)
                 phi = phi_grid(icell)

                 ! Warp analytique
                 if (lwarp) then
                    z0 = z_warp * (rcyl/dz%rref)**3 * cos(phi)
                 else if (ltilt) then
                    if (izone==1) then
                       z0 = rcyl * cos(phi) * tan(tilt_angle * deg_to_rad)
                    else
                       z0 = 0.0
                    endif
                 else
                    z0 = 0.0
                 endif

                 H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta
                 if (rcyl > dz%rmax) then
                    density = 0.0
                 else if (rcyl < dz%rmin) then
                    density = 0.0
                 else
                    density = cst_gaz(izone) * &
                         ( (rcyl/dz%Rc)**(-2*dz%surf) + (rcyl/dz%Rc)**(-2*dz%moins_gamma_exp) )**(-0.5) * &
                         exp( - (abs(z-z0)/h)**dz%vert_exponent)
                 endif
                 densite_gaz_tmp(icell) = density
              enddo ! k
           enddo bz_debris !j
        enddo !i

     endif ! dz%geometry

     if (lgap_Gaussian) then
        do icell=1, n_cells
           densite_gaz_tmp(icell) = densite_gaz_tmp(icell) * (1.0 - f_gap_Gaussian * &
                exp(-0.5 * ((r_grid(icell) - r_gap_Gaussian) / sigma_gap_Gaussian)**2 ))
        enddo
     endif

     !----------------------------------------
     ! Normalisation masse de gaz par zone
     !----------------------------------------

     ! Calcul de la masse de gaz de la zone
     mass = 0.
     do icell = 1, n_cells
        mass = mass + densite_gaz_tmp(icell) *  masse_mol_gaz * volume(icell)
     enddo
     mass =  mass * AU3_to_m3 * g_to_Msun

     ! Normalisation
     if (mass > 0.0) then ! pour le cas ou gas_to_dust = 0.
        facteur = dz%diskmass * dz%gas_to_dust / mass
        !     write(*,*) "VERIF gas mass: zone ",  izone, dz%diskmass * dz%gas_to_dust, mass, facteur

        ! Somme sur les zones pour densite finale
        do i=1,n_rad
           bz_gas_mass2 : do j=min(1,j_start),nz
              if (j==0) cycle
              do k=1, n_az
                 icell = cell_map(i,j,k)
                 densite_gaz(icell) = densite_gaz(icell) + densite_gaz_tmp(icell) * facteur
              enddo !k
              densite_gaz_midplane(i) = densite_gaz_midplane_tmp(i) + densite_gaz_midplane_tmp(i) * facteur
           enddo bz_gas_mass2
        enddo ! i
     endif
  enddo ! n_zones

  ! Ajout cavite vide
  if (lcavity) then
     do icell=1, n_cells
        surface = cavity%sclht * (r_grid(icell) / cavity%rref)**cavity%exp_beta
        if (abs(z_grid(icell)) > surface) then
           densite_gaz(icell) = 0.0_dp
        endif
     enddo
  endif

  ! Tableau de masse de gaz
  do icell=1,n_cells
     masse_gaz(icell) =  densite_gaz(icell) * masse_mol_gaz * volume(icell) * AU3_to_m3
  enddo
  write(*,*) 'Total  gas mass in model:', real(sum(masse_gaz) * g_to_Msun),' Msun'

  if (lcorrect_density) then
     write(*,*) "Correcting density ..."
     do icell=1, n_cells
        if ((r_grid(icell) >= correct_density_Rin).and.(r_grid(icell) <= correct_density_Rout)) then
           densite_gaz(icell) = densite_gaz(icell) *  correct_density_factor
           masse_gaz(icell) = masse_gaz(icell) *  correct_density_factor
        endif
     enddo
  end if

  return

end subroutine define_gas_density

!*************************************************************

subroutine define_dust_density()
! Calcule la table de densite
! Inclus stratification analytique
! Calcule les tableaux densite_pouss et masse
! et indique icell_not_empty
! C. Pinte : re-ecrit le 27/04/2013

  implicit none

  integer :: i,j, k, icell, l, izone, pop
  real(kind=dp), dimension(n_pop) :: cst, cst_pous
  real(kind=dp) :: rcyl, rsph, mass
  real(kind=dp) :: z, fact_exp, coeff_exp, density, OmegaTau, h_H2
  real(kind=dp) :: puffed, facteur, z0, phi, surface, norme

  real(kind=dp), dimension(n_grains_tot) :: correct_strat, N_tot, N_tot2

  real(kind=dp) :: rho0, ztilde, Dtilde, h, s_opt, somme

  type(disk_zone_type) :: dz
  type(dust_pop_type), pointer :: d_p

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

     if ((1.0_dp+1e-10_dp) * dz%rin >=  dz%rmax) then
        write(*,*) "ERROR: Rout must be larger than than Rin in zone", izone
        write(*,*) "Exiting"
        stop
     endif

     if (dz%geometry == 5) lwall = .true.

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
     if (lvariable_dust.and.(settling_type == 1)) then
        ! loi de puissance
        a_strat = max(a_strat,minval(r_grain))
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
           rho0 = densite_gaz_midplane(i) ! midplane density (j=0)

           !write(*,*) "     ", rcyl, rho0*masse_mol_gaz*cm_to_m**2, dust_pop(pop)%rho1g_avg
           !write(*,*) "s_opt", rcyl, s_opt/1000.

           bz : do j=j_start,nz
              if (j==0) cycle bz

              icell = cell_map(i,j,1) ! to compute factors independent of azimuth
              ! On calcule la densite au milieu de la cellule
              rcyl = r_grid(icell)
              z = z_grid(icell)

              H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta

              ! Puffed-up rim analytique
              if (lpuffed_rim) then
                 puffed = 1.0 + (puffed_rim_h - 1.0) / (exp((rcyl - puffed_rim_r)/puffed_rim_delta_r) + 1.0)
              else
                 puffed = 1.0
              endif

              if (dz%geometry == 1) then ! power-law
                 fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta)
              else if (dz%geometry == 2) then ! tappered-edge : dz%surf correspond a -gamma
                 fact_exp = (rcyl/dz%rref)**(dz%surf-dz%exp_beta) * exp( -(rcyl/dz%rc)**(2+dz%moins_gamma_exp) )
              endif
              coeff_exp = (2*(rcyl/dz%rref)**(2*dz%exp_beta))

              do k=1, n_az
                 icell = cell_map(i,j,k)
                 phi = phi_grid(icell)

                 ! Warp analytique
                 if (lwarp) then
                    z0 = z_warp * (rcyl/dz%rref)**3 * cos(phi)
                 else if (ltilt) then
                    if (izone==1) then
                       z0 = rcyl * cos(phi) * tan(tilt_angle * deg_to_rad)
                    else
                       z0 = 0.0
                    endif
                 else
                    z0 = 0.0
                 endif

                 ! Densite de la poussiere
                 do  l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    ! Settling a la Dubrulle
                    if (lvariable_dust.and.(settling_type == 2)) then
                       !h_H=(1d0/(1d0+gamma))**(0.25)*sqrt(alpha/(Omega*tau_f)) ! echelle de hauteur du rapport gaz/poussiere / H_gaz
                       !hd_H=h_H*(1d0+h_H**2)**(-0.5)                           ! echelle de hauteur de la poussiere / H_gaz
                       if (rho0 > tiny_dp) then
                          OmegaTau = omega_tau(rho0,H,l)
                          h_H2= sqrt(1./(1.+gamma)) * alpha/OmegaTau
                          correct_strat(l) = (1 + h_H2) / h_H2 ! (h_gas/h_dust)^2
                       else
                          correct_strat(l) = 1.0
                       endif
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
                    densite_pouss(l,icell) = density
                 enddo !l
              enddo !k
           enddo bz !j

           if ((lSigma_file).and.(izone==1)) then
              ! Normalisation pour densite de surface dans fichier
              ! todo : only works for k = 1
              do k=1,n_az
                 do  l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    somme = 0.0
                    do j=j_start,nz
                       if (j==0) cycle
                       icell = cell_map(i,j,k)
                       somme = somme + densite_pouss(l,icell)  *  (z_lim(i,j+1) - z_lim(i,j))
                    enddo ! j
                    if (somme > tiny_dp) then
                       do j=j_start,nz
                          if (j==0) cycle
                          icell = cell_map(i,j,k)
                          densite_pouss(l,icell) = densite_pouss(l,icell)  * Surface_density(i)/somme * nbre_grains(l)
                       enddo ! j
                    endif
                 enddo ! l
              enddo ! k
           endif


           if (lvariable_dust.and.(settling_type == 2)) then
              if (lspherical) then
                 write(*,*) "ERROR: settling following Dubrulle's prescription is only"
                 write(*,*) "implemented on a cylindrical grid so far"
                 write(*,*) "Exiting"
                 stop
              endif

              if ((rcyl > dz%rmin).and.(rcyl < dz%rmax)) then
                 ! Renormalisation pour  les cas ou il y a peu de resolution en z
                 do k=1, n_az
                    do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                       ! normalization en z
                       norme = 0.0
                       do j=j_start,nz
                          if (j==0) cycle
                          icell = cell_map(i,j,k)
                          norme = norme + densite_pouss(l,icell)
                       enddo !j

                       ! Si tous les grains sont sedimentes, on les met dans le plan median
                       if (norme < 1.0e-200_dp) then
                          icell = cell_map(i,1,k)
                          densite_pouss(l,icell)  = 1.0_dp
                          norme = 1.0_dp

                          write(*,*) "WARNING: Vertical settling unresolved for"
                          write(*,*) "grain larger than", r_grain(l), "at R > ", real(rcyl)
                       endif

                       do j=j_start,nz
                          if (j==0) cycle
                          icell = cell_map(i,j,k)
                          if (norme > tiny_dp) densite_pouss(l,icell) = densite_pouss(l,icell) / norme * rho0 * nbre_grains(l)
                       enddo !j
                    enddo ! l
                 enddo ! k
              endif ! test r
           endif ! settling==2

        enddo ! i

        if (lvariable_dust.and.(settling_type == 3)) then
           ! Si strat a la Seb. Fromang : on ecrase le tableau de densite
           ! je ne code que la dependence en z dans un premier temps puis normalise et ajoute la dependence en R et taille de grain

           if (lspherical) then
              write(*,*) "ERROR: settling following Fromang's prescription is only"
              write(*,*) "implemented on a cylindrical grid so far"
              write(*,*) "Exiting"
              stop
           endif

           do i=1, n_rad
              lwarning = .true.
              rho0 = densite_gaz_midplane(i) ! pour dependance en R : pb en coord sperique
              icell = cell_map(i,1,1)
              rcyl = r_grid(icell)
              H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta

              if ((rcyl > dz%rmin).and.(rcyl < dz%rmax)) then
                 do k=1, n_az

                    ! Renormalisation pour  les cas ou il y a peu de resolution en z
                    do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                       !calculate omega_tau in the disk midplane
                       OmegaTau = omega_tau(rho0,H,l)

                       do j=j_start,nz ! dependence en z uniquement ici !!!!
                          if (j==0) cycle
                          icell = cell_map(i,j,k)

                          !calculate h & z/h
                          z = z_grid(icell)
                          Ztilde=z/H

                          ! Fit Gaussien du profile de densite
                          !densite_pouss(l,icell)=  exp(-(1+OmegaTau/Dtilde) * (Ztilde**2/2.))

                          ! Coefficient de diffusion constant
                          densite_pouss(l,icell)=  exp( -OmegaTau/Dtilde * (exp(Ztilde**2/2.)-1) - Ztilde**2/2 )  ! formule 19
                       enddo!j

                       ! normalization en z
                       norme = 0.0
                       do j=j_start,nz
                          if (j==0) cycle
                          icell = cell_map(i,j,k)
                          norme = norme + densite_pouss(l,icell)
                       enddo !j

                       ! Si tous les grains sont sedimentes, on les met dans le plan median
                       if (norme < 1e-200_dp) then
                          icell = cell_map(i,1,k)
                          densite_pouss(l,icell)  = 1.0_dp
                          norme = 1.0_dp

                          if (lwarning) then
                             write(*,*)
                             write(*,*) "WARNING : Vertical settling unresolved for"
                             write(*,*) "grain larger than", r_grain(l), "at R > ", real(rcyl)
                             lwarning = .false. ! on ne fait un warning qu'1 fois par rayon
                          endif
                       endif

                       do j=j_start,nz
                          if (j==0) cycle
                          icell = cell_map(i,j,k)
                          if (norme > tiny_dp) densite_pouss(l,icell) = densite_pouss(l,icell) / norme * rho0 * nbre_grains(l)
                       enddo !j

                    enddo ! l
                 enddo ! k
              endif ! test r
           enddo ! i
        endif ! Settling Fromang

        lwarning = .true.
        if (lmigration) then
           ! distribution en taille de grains avant la migration
           do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
              do icell=1,n_cells
                 N_tot(l) = N_tot(l) + densite_pouss(l,icell) * volume(icell)
              enddo ! icell
           enddo !l

           do i=1, n_rad
              do k=1, n_az
                 rho0 = densite_gaz(cell_map(i,1,k)) ! pour dependance en R : pb en coord sperique
                 !s_opt = rho_g * cs / (rho * Omega)    ! cs = H * Omega ! on doit trouver 1mm vers 50AU
                 !omega_tau= dust_pop(ipop)%rho1g_avg*(r_grain(l)*mum_to_cm) / (rho * masse_mol_gaz/m_to_cm**3 * H*AU_to_cm)
                 icell = cell_map(i,1,k)
                 rcyl = r_grid(icell)
                 H = dz%sclht * (rcyl/dz%rref)**dz%exp_beta
                 s_opt = (rho0*masse_mol_gaz*cm_to_m**3  /dust_pop(pop)%rho1g_avg) *  H * AU_to_m * m_to_mum

                 write(*,*) "r=", rcyl, "a_migration =", s_opt

                 if ((s_opt < dust_pop(pop)%amin).and.(lwarning)) then
                    write(*,*)
                    write(*,*) "WARNING: a_migration = ", s_opt
                    write(*,*) "is smaller than amin for dust pop #", pop
                    write(*,*) "MCFOST will exit with an error as there are no smaller grains"
                    if (s_opt < tiny_dp) write(*,*) "is your gas-to-dust ratio = 0 ?"
                    lwarning = .false.
                 endif

                 do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    if (r_grain(l) > s_opt) then ! grains plus gros que taille optimale de migration
                       do j=j_start,nz
                          if (j==0) cycle
                          icell = cell_map(i,j,k)
                          densite_pouss(l,icell) = 0.0
                       enddo !j
                    endif
                 enddo ! l
              enddo ! k
           enddo !i

           ! distribution en taille de grains apres la migration
           do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
              do icell=1,n_cells
                 N_tot2(l) = N_tot2(l) + densite_pouss(l,icell) * volume(icell)
              enddo ! i
           enddo ! l

           ! Renormalisation : on garde le meme nombre de grains par taille que avant la migration
           do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
              if (N_tot2(l) > tiny_dp) then
                 densite_pouss(l,:) = densite_pouss(l,:) * N_tot(l)/N_tot2(l)
              endif
           enddo ! l

        endif ! migration

     else if (dz%geometry == 3) then ! envelope
        do i=1, n_rad
           bz_env : do j=j_start,nz
              if (j==0) cycle bz_env
              do k=1, n_az
                 icell = cell_map(i,j,k)
                 ! On calcule la densite au milieu de la cellule
                 rcyl = r_grid(icell)
                 z = z_grid(icell)
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
                    densite_pouss(l,icell) = density
                 enddo !l

              enddo ! k
           enddo bz_env !j
        enddo ! i

     else if (dz%geometry == 4) then ! disque de debris
        do i=1, n_rad
           bz_debris : do j=j_start,nz
              if (j==0) cycle bz_debris
              do k=1, n_az
                 icell = cell_map(i,j,k)
                 ! On calcule la densite au milieu de la cellule
                 rcyl = r_grid(icell)
                 z = z_grid(icell)
                 phi = phi_grid(icell)

                 h = dz%sclht * (rcyl/dz%Rref)**dz%exp_beta

                 !R(r) = (  (r/rc)^-2alpha_in + (r/rc)^-2alpha_out )^-1/2
                 !Z(r,z) =  exp( - (abs(z)/h(r))^gamma  )


                 ! Warp analytique
                 if (lwarp) then
                    z0 = z_warp * (rcyl/dz%rref)**3 * cos(phi)
                 else if (ltilt) then
                    if (izone==1) then
                       z0 = rcyl * cos(phi) * tan(tilt_angle * deg_to_rad)
                    else
                       z0 = 0.0
                    endif
                 else
                    z0 = 0.0
                 endif

                 do l=dust_pop(pop)%ind_debut,dust_pop(pop)%ind_fin
                    if (rcyl > dz%rmax) then
                       density = 0.0
                    else if (rcyl < dz%rmin) then
                       density = 0.0
                    else
                       density = nbre_grains(l) * cst_pous(pop) * &
                            ( (rcyl/dz%Rc)**(-2*dz%surf) + (rcyl/dz%Rc)**(-2*dz%moins_gamma_exp) )**(-0.5) * &
                            exp( - (abs(z -z0)/h)**dz%vert_exponent)
                    endif
                    densite_pouss(l,icell) = density
                 enddo ! l


              enddo !k
           enddo bz_debris !j
        enddo !i

     endif ! dz%geometry

  enddo ! pop

  ! Ajout cavite vide
  if (lcavity) then
     do icell=1,n_cells
        surface = cavity%sclht * (r_grid(icell) / cavity%rref)**cavity%exp_beta
        if (abs(z_grid(icell)) > surface) then
           densite_pouss(:,icell) = 0.0_dp
        endif
     enddo
  endif

  if (lgap_Gaussian) then
     do icell=1, n_cells
        densite_pouss(:,icell) = densite_pouss(:,icell) * (1.0 - f_gap_Gaussian * &
             exp(-0.5 * ((r_grid(icell) - r_gap_Gaussian) / sigma_gap_Gaussian)**2 ))
     enddo
  endif

  search_not_empty : do l=1,n_grains_tot
     do icell=1,n_cells
        if (densite_pouss(l,icell) > 0.0_dp) then
           icell_not_empty = icell
           exit search_not_empty
        endif
     enddo
  enddo search_not_empty


  ! Normalisation poussiere: re-calcul masse totale par population a partir de la densite (utile quand edge /= 0)
  do pop=1, n_pop
     izone=dust_pop(pop)%zone
     dz=disk_zone(izone)

     if (dz%geometry /= 5) then ! pas de wall ici
        d_p => dust_pop(pop)
        mass = 0.0

        do icell=1,n_cells
           do l=d_p%ind_debut,d_p%ind_fin
              mass=mass + densite_pouss(l,icell) * M_grain(l) * volume(icell)
           enddo !l
        enddo !icell
        mass =  mass * AU3_to_cm3 * g_to_Msun

        if (mass < tiny_dp) then
           write(*,*)
           write(*,*) "ERROR : something went wrong, there is no dust in the disk"
           write(*,*) "Exiting." ; stop
        endif

        facteur = d_p%masse / mass

        do icell=1,n_cells
           do l=d_p%ind_debut,d_p%ind_fin
              densite_pouss(l,icell) = densite_pouss(l,icell) * facteur
              masse(icell) = masse(icell) + densite_pouss(l,icell) * M_grain(l) * volume(icell)
           enddo !l
        enddo ! icell

     endif ! test wall
  enddo ! pop

  masse(:) = masse(:) * AU3_to_cm3
  write(*,*) 'Total dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'

  if (lcorrect_density) then
     write(*,*) "Correcting density ..."
     do i=1, n_rad
        icell = cell_map(i,1,1)
        if ((r_grid(icell) >= correct_density_Rin).and.(r_grid(icell) <= correct_density_Rout)) then
           do j=j_start,nz
              if (j==0) cycle
              do k=1, n_az
                 icell = cell_map(i,j,k)
                 densite_pouss(:,icell) = densite_pouss(:,icell) * correct_density_factor
                 masse(icell) = masse(icell) *  correct_density_factor
              enddo !k
           enddo ! j
        endif
     enddo

      write(*,*) 'Total corrected dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'
  endif

  ! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  ! Valable que dans le cas cylindrique mais pas de pb dans le cas spherique
  !if (lcylindrical) densite_pouss(:,nz+1,:,:) = densite_pouss(:,nz,:,:)

  return

end subroutine define_dust_density

!*************************************************************

subroutine define_density_wall3D()
  ! superpose un mur avec une forme en cosinus
  ! sur le disque
  ! C. Pinte  25/01/2011

  integer :: pop, l, izone, alloc_status, icell
  type(disk_zone_type) :: dz
  type(dust_pop_type), pointer :: d_p

  real(kind=dp) :: rcyl, z, phi, density, facteur, hh, mass, h_wall

  real(kind=dp), dimension(:,:), allocatable :: density_wall
  real(kind=dp), dimension(:), allocatable :: masse_wall

  write(*,*) "*********************************************************"
  write(*,*) "Adding 3D wall structure ...."

  allocate(density_wall(n_cells,n_grains_tot), stat=alloc_status)
  allocate(masse_wall(n_cells), stat=alloc_status)
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
        h_wall = dz%sclht
        write(*,*) "h_wall =", real(h_wall)


        do icell=1, n_cells
           ! On calcule la densite au milieu de la cellule
           rcyl = r_grid(icell)
           z = z_grid(icell)
           phi = phi_grid(icell)

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
                 density_wall(icell,l) = density
                 !if (density > 0.) write(*,*) i,j, k, l, density_wall(i,j,k,l)
              else
                 density_wall(icell,l) = 0.0
              endif
           enddo ! l

        enddo !icell

     endif ! wall
  enddo ! pop

  ! Normalisation de la masse du mur
  masse_wall(:) = 0.0

  do pop=1, n_pop
     izone=dust_pop(pop)%zone
     dz=disk_zone(izone)

     if (dz%geometry == 3) then ! wall
        d_p => dust_pop(pop)
        mass = 0.0

        do icell=1,n_cells
           do l=d_p%ind_debut,d_p%ind_fin
              mass=mass + density_wall(icell,l) * M_grain(l) * (volume(icell) * AU3_to_cm3)
           enddo !l
        enddo !icell
        mass =  mass*g_to_Msun

        facteur = d_p%masse / mass

        do icell=1,n_cells
           do l=d_p%ind_debut,d_p%ind_fin
              density_wall(icell,l) = density_wall(icell,l) * facteur
              masse_wall(icell) = masse_wall(icell) + density_wall(icell,l) * M_grain(l) * volume(icell)
           enddo !l
        enddo ! icell

     endif ! zone = wall
  enddo ! pop

  masse_wall(:) = masse_wall(:) * AU3_to_cm3
  write(*,*) 'Wall dust mass:', real(sum(masse_wall)*g_to_Msun),' Msun'

  ! superposition du mur sur le disque
  densite_pouss(:,:) = densite_pouss(:,:) + density_wall(:,:)
  masse(:) = masse(:) + masse_wall(:)

  write(*,*) 'Total dust mass in model:', real(sum(masse)*g_to_Msun),' Msun'
  write(*,*) "Done"
  write(*,*) "*********************************************************"


  return

end subroutine define_density_wall3D

!*************************************************************

subroutine densite_eqdiff()
  ! only 2D : old routine, has not been used in a long time

  implicit none

  real, parameter :: G = 6.672e-8
  real, parameter :: gas_dust = 100

  integer :: i,j, k, jj, icell
  real :: cst, cst_pous, cst_gaz, M_star
  real :: fact_exp, c_sound, coeff_grav, omega, D0, eps, pas_z, somme1, somme2, correct,test

  real, dimension(1) :: y

  real, dimension(nz) :: z

  real, dimension(n_grains_tot, nz) :: correct_strat
  real, dimension(nz) :: rho

  type(disk_zone_type) :: dz

  ! Vitesses (en m/s) (Dullemond et Dubrulle)
  real, parameter :: v_sound = 380.0 ! m/s
  ! Rayon de definition des vitesses
  real, parameter :: rref_v = 50. ! au


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


  do i=1, n_rad
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
     if (lvariable_dust) then
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
        enddo !k
     else!lvariable_dust
        correct_strat = 1.0
     endif!lvariable_dust

     ! Calcul opacite et probabilite de diffusion
     do j=1,nz
        icell = cell_map(i,j,1) ! only 2D
        do  k=1,n_grains_tot
           densite_pouss(k,icell) = nbre_grains(k) * correct_strat(k,j) * rho(j)
        enddo !k
     enddo !j
  enddo !i

!***************
! Remplissage a zero pour z > zmax que l'en envoie sur l'indice j=0
  do i=1,n_rad
     do  k=1,n_grains_tot
        densite_pouss(k,cell_map(i,nz+1,1)) = densite_pouss(k,cell_map(i,nz,1))
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
  real :: rho_gaz

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

!********************************************************************

subroutine densite_file()
  ! Nouvelle routine pour lire les grilles de densite
  ! calculees par Yorick directement a partir des donnees SPH (ou autre)
  ! Les donnees sont directement lissees sur la grille de MCFOST
  ! C. Pinte
  ! 12/0/09
  ! Mise a jour 9/11/13

  use grains
  use utils

  implicit none

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,npixels,j, hdutype, bitpix
  integer :: nullval
  integer, dimension(4) :: naxes
  logical :: anynull, l3D_file
  character(len=80) :: comment

  integer :: k, l, i, n_a, read_n_a, jj, icell, phik
  real(kind=dp) :: somme, mass, facteur
  real :: a, tmp

  real, dimension(:,:,:,:), allocatable :: sph_dens ! (n_rad,nz,n_az,n_a)
  real, dimension(:), allocatable :: a_sph, n_a_sph, log_a_sph, log_n_a_sph ! n_a

  real(kind=dp), dimension(:,:,:,:), allocatable :: sph_dens_dp
  real(kind=dp), dimension(:), allocatable :: a_sph_dp
  real(kind=dp) :: f

  type(disk_zone_type) :: dz

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
  call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
  if (nfound /= 4) then
     write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
     write(*,*) 'of '//trim(density_file)//' file. Exiting.'
     write(*,*) "I found", nfound, "axis instead of 4"
     stop
  endif

  if ((naxes(1) /= n_rad).or.((naxes(2) /= nz).and.(naxes(2) /= 2*nz+1)).or.(naxes(3) /= n_az) ) then
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

  if (naxes(2) == 2*nz+1) then
     l3D_file = .true.
  else
     write(*,*) "The density file only has > 0 z, making it symmetric"
     l3D_file = .false.
  endif


  if (l3D_file) then
     allocate(sph_dens(n_rad,-nz:nz,n_az,n_a), a_sph(n_a), n_a_sph(n_a))
  else
     allocate(sph_dens(n_rad,nz,n_az,n_a), a_sph(n_a), n_a_sph(n_a))
  endif
  sph_dens = 0.0 ; a_sph = 0.0 ; n_a_sph = 0.0

  bitpix = 0
  call ftgkyj(unit,"bitpix",bitpix,comment,status)

  ! read_image
  if (bitpix==-32) then
     call ftgpve(unit,group,firstpix,npixels,nullval,sph_dens,anynull,status)
  else if (bitpix==-64) then
     if (l3D_file) then
        allocate(sph_dens_dp(n_rad,-nz:nz,n_az,n_a))
     else
        allocate(sph_dens_dp(n_rad,nz,n_az,n_a))
     endif
     sph_dens_dp = 0.0_dp
     call ftgpvd(unit,group,firstpix,npixels,nullval,sph_dens_dp,anynull,status)
     sph_dens = real(sph_dens_dp,kind=sp)
     deallocate(sph_dens_dp)
  else
     write(*,*) "ERROR: cannot read bitpix in fits file"
     stop
  endif

  write(*,*) "Density range:", minval(sph_dens), maxval(sph_dens)

  ! Au cas ou : on elimine les valeurs a 0
  sph_dens = sph_dens/maxval(sph_dens) ! normalization avant d'ajouter une constante
  sph_dens = max(sph_dens,1e10*tiny_real)

  if (n_a > 1) then

     read_n_a = 0
     call ftgkyj(unit,"read_n_a",read_n_a,comment,status)

     write(*,*) "read_n_a", read_n_a

     !---------------------------------------------------------
     ! HDU 2 : grain sizes
     !---------------------------------------------------------
     !  move to next hdu
     call ftmrhd(unit,1,hdutype,status)
     nfound=1
     ! Check dimensions
     naxes(:) = 0
     call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
     if (nfound /= 1) then
        write(*,*) 'READ_IMAGE did not find 1 dimension in HDU 2'
        write(*,*) 'HDU 2 has', nfound, 'dimensions.'
        write(*,*) 'Exiting.'
        write(*,*)
        stop
     endif
     if ((naxes(1) /= n_a)) then
        write(*,*) "Error : HDU 2 does not have the right dimension"
        write(*,*) "It has ", naxes(1), "instead of ", n_a
        write(*,*) "Exiting."
        stop
     endif
     npixels=naxes(1)

     bitpix=0
     call ftgkyj(unit,"bitpix",bitpix,comment,status)

     ! read_image
     if (bitpix==-32) then
        call ftgpve(unit,group,firstpix,npixels,nullval,a_sph,anynull,status)
     else if (bitpix==-64) then
        allocate(a_sph_dp(n_a)) ; a_sph_dp = 0.0_dp
        call ftgpvd(unit,group,firstpix,npixels,nullval,a_sph_dp,anynull,status)
        a_sph = real(a_sph_dp,kind=sp)
        deallocate(a_sph_dp)
     else
        write(*,*) "ERROR: cannot read bitpix in fits file"
        stop
     endif

     ! read_image

     ! On verifie que les grains sont tries
     do i=1, n_a-1
        if (a_sph(i) >= a_sph(i+1)) then
           write(*,*) "ERROR : grains must be ordered from small to large"
           write(*,*) "I found the follwing grain sizes in the fits :"
           do j=1, n_a
              write(*,*) a_sph(j)
           enddo
           write(*,*) "Exiting"
           stop
        endif
     enddo


     ! On lit au besoin la distribution en taille (dn(a) / da)
     if (read_n_a==1) then
        write(*,*) "Reading grain size distribution from fits file"

        !---------------------------------------------------------
        ! HDU 3 : nombre de grains
        !---------------------------------------------------------
        !  move to next hdu
        call ftmrhd(unit,1,hdutype,status)

        nfound = 0 ; naxes = 0 ;
        ! Check dimensions
        call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
        if (nfound /= 1) then
           write(*,*) 'READ_IMAGE did not find 1 dimension in HDU 2'
           write(*,*) 'HDU 3 has', nfound, 'dimensions.'
           write(*,*) 'Exiting.'
           stop
        endif
        if ((naxes(1) /= n_a)) then
           write(*,*) "Error : HDU 3 does not have the right dimension"
           write(*,*) "It has ", naxes(1), "instead of ", n_a
           write(*,*) "Exiting."
           stop
        endif
        npixels=naxes(1)

        bitpix = 0
        call ftgkyj(unit,"bitpix",bitpix,comment,status)

        ! read_image
        if (bitpix==-32) then
           call ftgpve(unit,group,firstpix,npixels,nullval,n_a_sph,anynull,status)
        else if (bitpix==-64) then
           allocate(a_sph_dp(n_a)) ; a_sph_dp = 0.0_dp
           call ftgpvd(unit,group,firstpix,npixels,nullval,a_sph_dp,anynull,status)
           n_a_sph = real(a_sph_dp,kind=sp)
           deallocate(a_sph_dp)
        else
           write(*,*) "ERROR: cannot read bitpix in fits file"
           stop
        endif

        tmp = sum(n_a_sph)
        write(*,*) "The following grain sizes were found in the fits file:"
        do i=1,n_a
           write(*,*) i, a_sph(i), "microns", n_a_sph(i) / tmp
        enddo
        write(*,*) "They will be used to set the integrated grain size distribution"


        if (n_pop > 1) then
           write(*,*) "ERROR : density fits interface only works for 1 dust pop"
           stop
        endif

        allocate(log_a_sph(n_a), log_n_a_sph(n_a))
        log_a_sph = log(a_sph) ; log_n_a_sph = log(n_a_sph)

        ! Multiplication par a car da = a.dln(a)
        do k=1, n_grains_tot
           write(*,*) k
           a = r_grain(k)
           ! todo : peut etre optimise sans interp
           nbre_grains(k) = exp( interp(log_n_a_sph, log_a_sph, log(a)) )  * a
        enddo !k

        ! Normalisation de tous les grains au sein d'une pop
        nbre_grains = nbre_grains / sum(nbre_grains)
     else
        write(*,*) "The following grain sizes were found in the fits file:"
        do i=1,n_a
           write(*,*) i, a_sph(i), "microns"
        enddo
        write(*,*) "Using values from parameter file to set the integrated grain size distribution"
     endif

  else ! n_a > 1
     write(*,*) "Using grain size distribution from parameter file"
     a_sph(1) = 1.0
  endif

  call ftclos(unit, status)
  call ftfiou(unit, status)

  ! Densite du gaz : gaz = plus petites particules
  dz = disk_zone(1)
  !densite_gaz(:,1:nz,:) = sph_dens(:,:,:,1) ! marche pas, bizarre ???
  do k=1, n_az
     do j=j_start,nz
        if (l3D_file) then
           jj = j
        else
           jj = abs(j)
        endif
        if (j==0) then
           !densite_gaz(cell_map(i,j,k)) =0.0
        else
           do i=1, n_rad
              densite_gaz(cell_map(i,j,k)) = sph_dens(i,jj,k,1) ! gaz = plus petites particules
           enddo
        endif
     enddo
  enddo

  ! Calcul de la masse de gaz de la zone
  mass = 0.
  do icell=1,n_cells
     mass = mass + densite_gaz(icell) *  masse_mol_gaz * volume(icell)
  enddo !icell
  mass =  mass * AU3_to_m3 * g_to_Msun

  ! Normalisation
  if (mass > 0.0) then ! pour le cas ou gas_to_dust = 0.
     facteur = dz%diskmass * dz%gas_to_dust / mass
     ! Somme sur les zones pour densite finale
     do icell=1,n_cells
        densite_gaz(icell) = densite_gaz(icell) * facteur
     enddo ! icell
  endif

  ! Tableau de masse de gaz
  do icell=1,n_cells
     do k=1, n_az
        masse_gaz(icell) =  densite_gaz(icell) * masse_mol_gaz * volume(icell) * AU3_to_m3
     enddo !k
  enddo ! icell

  if (lvariable_dust) then
     write(*,*) "Differential spatial distribution"
     l=1
     do k=1,n_grains_tot
        if (r_grain(k) < a_sph(1)) then  ! Petits grains
           do phik=1, n_az
              do j=j_start,nz
                 if (j==0) cycle
                 if (l3D_file) then
                    jj = j
                 else
                    jj = abs(j)
                 endif
                 do i=1, n_rad
                    densite_pouss(k,cell_map(i,j,phik)) = sph_dens(i,jj,phik,1)
                 enddo ! phik
              enddo ! j
           enddo ! i
        else if (r_grain(k) > a_sph(n_a)) then ! Gros grains

           do phik=1, n_az
              do j=j_start,nz
                 if (l3D_file) then
                    jj = j
                 else
                    jj = abs(j)
                 endif
                 if (j==0) cycle
                 do i=1, n_rad
                    densite_pouss(k,cell_map(i,j,phik)) = sph_dens(i,jj,phik,n_a)
                 enddo ! phik
              enddo ! j
           enddo ! i
        else  ! Autres grains : interpolation
           if (r_grain(k) > a_sph(l+1)) l = l+1
           f = (r_grain(k)-a_sph(l))/(a_sph(l+1)-a_sph(l))

           do phik=1, n_az
              do j=j_start,nz
                 if (j==0) cycle
                 if (l3D_file) then
                    jj = j
                 else
                    jj = abs(j)
                 endif
                 do i=1, n_rad
                    densite_pouss(k,cell_map(i,j,phik)) = sph_dens(i,jj,phik,l) + f * &
                         ( sph_dens(i,jj,phik,l+1) -  sph_dens(i,jj,phik,l) )
                 enddo ! phik
              enddo ! j
           enddo ! i
        endif
     enddo
  else ! Tous les grains suivent le gas
     write(*,*) "Constant spatial distribution"
     do k=1,n_grains_tot
        do phik=1, n_az
           do j=j_start,nz
              if (j==0) cycle
              if (l3D_file) then
                 jj = j
              else
                 jj = abs(j)
              endif
              do i=1, n_rad
                 densite_pouss(k,cell_map(i,j,phik)) = sph_dens(i,jj,phik,1)
              enddo ! phik
           enddo ! j
        enddo ! i
     enddo ! k
  endif  !lvariable_dust

  ! Normalisation : on a 1 grain de chaque taille dans le disque
  do l=1,n_grains_tot
     somme=0.0

     do icell=1,n_cells
        if (densite_pouss(l,icell) <= 0.0) densite_pouss(l,icell) = tiny_dp
        somme=somme+densite_pouss(l,icell)*volume(icell)
     enddo !icell
     densite_pouss(l,:) = (densite_pouss(l,:)/somme)
  enddo !l

  ! Normalisation : on a 1 grain en tout dans le disque
  do l=1,n_grains_tot
     densite_pouss(l,:) = (densite_pouss(l,:)/somme)*nbre_grains(l)
  enddo

  search_not_empty : do l=1,n_grains_tot
     do icell=1, n_cells
        if (densite_pouss(l,icell) > 0.0_dp) then
           icell_not_empty = icell
           exit search_not_empty
        endif
     enddo !icell
  enddo search_not_empty


  ! Normalisation : Calcul masse totale
  mass = 0.0
  do icell=1,n_cells
     do l=1,n_grains_tot
        mass=mass + densite_pouss(l,icell) * M_grain(l) * (volume(icell) * AU3_to_cm3)
     enddo !l
  enddo !icell
  mass =  mass/Msun_to_g
  densite_pouss(:,:) = densite_pouss(:,:) * diskmass/mass

  do icell=1,n_cells
     do l=1,n_grains_tot
        masse(icell) = masse(icell) + densite_pouss(l,icell) * M_grain(l) * volume(icell)
     enddo !l
  enddo ! icell
  masse(:) = masse(:) * AU3_to_cm3

  write(*,*) 'Total  gas mass in model:', real(sum(masse_gaz) * g_to_Msun),' Msun'
  write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'
  deallocate(sph_dens,a_sph)

  !write(*,*) "MODIFYING 3D DENSITY !!!"
  !k = 36
  !densite_pouss(:,:,k,:) = densite_pouss(:,:,k,:)* 1e20
  !densite_gaz(:,:,k) = densite_gaz(:,:,k)* 1e20

  return

end subroutine densite_file

!**********************************************************

subroutine read_Sigma_file()
  ! Nouvelle routine pour lire une densite de surface
  ! C. Pinte
  ! 22/01/15

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,npixels, bitpix
  integer :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull
  character(len=80) :: comment
  real, dimension(:), allocatable :: sigma_sp

  ! Lecture donnees
  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)
  write(*,*) "Reading surface density file : "//trim(sigma_file)

  if (n_zones > 1) then
     write(*,*) "This surface density will be applied to zone 1"
     ! todo : a verifier apres reordering
  endif

  readwrite=0
  call ftopen(unit,sigma_file,readwrite,blocksize,status)
  if (status /= 0) then ! le fichier temperature n'existe pas
     write(*,*) "ERROR : surface density file needed"
     stop
  endif

  group=1
  firstpix=1
  nullval=-999

  ! determine the size of density file
  call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
  if (nfound > 2) then
     write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
     write(*,*) 'of '//trim(density_file)//' file. Exiting.'
     write(*,*) "nfound = ", nfound, "instead of 1 or 2"
     stop
  endif

  if ((naxes(1) /= n_rad)) then
     write(*,*) "Error : "//trim(sigma_file)//" does not have the"
     write(*,*) "right dimensions. Exiting."
     write(*,*) "Axis #1 (radius) :  fits_file vs mcfost_grid"
     write(*,*) naxes(1), n_rad
     stop
  endif

  if (nfound==2) then
     if ((naxes(2) /= n_az)) then
        write(*,*) "Error : "//trim(sigma_file)//" does not have the"
        write(*,*) "right dimensions. Exiting."
        write(*,*) "Axis #2 (azimuth):  fits_file vs mcfost_grid"
        write(*,*) naxes(1), n_az
        stop
     endif
     npixels = naxes(1) * naxes(2)
  else
     npixels = naxes(1)
  endif

  bitpix = 0
  call ftgkyj(unit,"bitpix",bitpix,comment,status)

  ! read_image
  if (bitpix==-32) then
     allocate(sigma_sp(n_rad))
     sigma_sp = 0.0_dp
     call ftgpve(unit,group,firstpix,npixels,nullval,sigma_sp,anynull,status)
     surface_density = real(sigma_sp,kind=dp)
     deallocate(sigma_sp)
  else if (bitpix==-64) then
     call ftgpvd(unit,group,firstpix,npixels,nullval,surface_density,anynull,status)
  else
     write(*,*) "ERROR: cannot read bitpix in fits file"
     stop
  endif

  ! Au cas ou
  surface_density = max(surface_density,tiny_dp)

  call ftclos(unit, status)
  call ftfiou(unit, status)

  return

end subroutine read_Sigma_file

!**********************************************************************

real(kind=dp) function omega_tau(rho,H,l)
  ! Pour les calculs de Sebastien Fromang
  ! rho doit etre en g.cm-3

  real(kind=dp), intent(in) :: rho, H
  integer, intent(in) :: l

  integer :: ipop

  ipop = grain(l)%pop
  !write(*,*) ipop, dust_pop(ipop)%rho1g_avg, rho
  if (rho > tiny_dp) then
     omega_tau = dust_pop(ipop)%rho1g_avg*(r_grain(l)*mum_to_cm) / (rho * masse_mol_gaz/m_to_cm**3 * H*AU_to_cm)
  else
     omega_tau = huge_dp
  endif

  return

end function omega_tau

!**********************************************************************

subroutine densite_Seb_Charnoz()

  integer :: Nr_Seb, Nz_Seb, Na_Seb
  real(kind=dp), dimension(n_grains_tot) :: density_Seb, taille_grains_Seb, N_grains
  real(kind=dp) :: Rmin, Dr, Zmin, Dz
  real(kind=dp) :: Rmin_mcfost, Dr_mcfost, Zmin_mcfost, Dz_mcfost

  real(kind=dp) :: Somme
  integer :: ii, jj, i, j, k, l, icell

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

        icell = cell_map(i,j,1) ! only 2D
        densite_pouss(:,icell) = density_Seb(:) / (volume(icell)*AU3_to_cm3) ! comversion en densite volumique
        Somme = Somme +  1.6 * 4.*pi/3. *  (mum_to_cm)**3 * sum( density_Seb(:) * r_grain(:)**3 )
     enddo ! j
  enddo !i
  write(*,*) "Dust mass from Seb's file :", real(Somme * g_to_Msun), "Msun"

  search_not_empty : do l=1,n_grains_tot
     do icell=1, n_cells
        if (densite_pouss(l,icell) > 0.0_dp) then
           icell_not_empty = icell
           exit search_not_empty
        endif
     enddo
  enddo search_not_empty

  ! Re-population du tableau de grains
  ! les methodes de chauffages etc, ne changent pas

  do k=1, n_grains_tot
     N_grains(k) = sum(densite_pouss(k,:))
  enddo
  nbre_grains(:) = N_grains(:)/sum(N_grains(:))


  do icell=1,n_cells
     do l=1,n_grains_tot
        masse(icell) = masse(icell) + densite_pouss(l,icell) * M_grain(l) * volume(icell)
     enddo !l
  enddo ! icell

  masse(:) = masse(:) * AU3_to_cm3 * 1600./3500 ! TMP

  write(*,*) 'Total dust mass in model  :', real(sum(masse)*g_to_Msun),'Msun'

  write(*,*) "Done"
  write(*,*) "***********************************************"

end subroutine densite_Seb_Charnoz

!**********************************************************

subroutine densite_Seb_Charnoz2()

  implicit none

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j
  integer :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull

  integer :: k, l, i, icell
  real(kind=dp) :: somme, somme2

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
  call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
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
  dens = dens + tiny_real

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
     do i=1, n_rad
        do j=1,nz
           densite_pouss(k,cell_map(i,j,1)) = dens(i,j) ! only 2D
        enddo
     enddo
  enddo

  ! kg/m^3  ---> part/cm^3
  densite_pouss = densite_pouss / ( (cm_to_m)**3  * dust_pop(1)%avg_grain_mass * 1e3)

  ! BUG correction : il manque un facteur
  densite_pouss = densite_pouss * 1e-6


  do l=1,n_grains_tot
     densite_pouss(l,:) = densite_pouss(l,:)*nbre_grains(l)
  enddo

  search_not_empty : do l=1,n_grains_tot
     do icell=1, n_cells
        if (densite_pouss(l,icell) > 0.0_dp) then
           icell_not_empty = icell
           exit search_not_empty
        endif
     enddo !icell
  enddo search_not_empty

  write(*,*) "Done"

  do icell=1,n_cells
     do l=1,n_grains_tot
        masse(icell) = masse(icell) + densite_pouss(l,icell) * M_grain(l) * volume(icell)
     enddo !l
  enddo ! icell

  masse(:) = masse(:) * AU3_to_cm3

  write(*,*) 'Total dust mass in model :', real(sum(masse)*g_to_Msun),' Msun'
  write(*,*) "Density from Seb. Charnoz set up OK"

  return

end subroutine densite_Seb_Charnoz2

!**********************************************************************

subroutine remove_specie()

  implicit none

  integer :: k, icell
  real :: mass

  write(*,*) "Removing specie", specie_removed, "where T >", T_rm

  do icell=1,n_cells
     do k=1,n_grains_tot
        if (grain(k)%pop==specie_removed) then
           if (Temperature(icell) > T_rm) densite_pouss(k,icell) = 0.0
        endif
     enddo
  enddo

  mass = 0.0
  do icell=1,n_cells
     do k=1,n_grains_tot
        mass=mass + densite_pouss(k,icell) * M_grain(k) * (volume(icell) * AU3_to_cm3)
     enddo
  enddo
  mass =  mass/Msun_to_g

  write(*,*) 'New total dust mass in model :', mass,' Msun'

  return

end subroutine remove_specie

!************************************************************

end module density
