module disk_physics

  use grains
  use opacity
  use em_th
  use prop_star
  use molecular_emission, only : densite_gaz
  use constantes

contains


subroutine compute_othin_sublimation_radius()
  ! Dans le cas optiquement mince, ne dépend que de la température de l'étoile

  implicit none

  real(kind=db) :: E_dust, E_etoile
  real :: cst, wl, delta_wl
  integer :: lambda, i, icell
  real(kind=db) :: sublimation_radius, coeff_exp, cst_wl

  ! TODO : pb de normalization spectre_etoiles si Teff n'est pas celle du spectre en mode non-bb


  !  if (lvariable_dust) then
  !     write(*,*) "Sublimation radius calculation not implemented"
  !     write(*,*) "in case of statification"
  !     stop
  !  endif
  ! --> kappa moyenne verticalement

  !  Emission poussière
  E_dust = 0.0
  cst=cst_th/dust_pop(1)%T_sub
  cst=cst_th/1500.

  icell = cell_map(1,1,1)

  do lambda=1, n_lambda
     ! longueur d'onde en metre
     wl = tab_lambda(lambda)*1.e-6
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     cst_wl=cst/wl
     if (cst_wl < 500.0) then
        coeff_exp=exp(cst_wl)
        E_dust = E_dust + 4.0 * kappa_abs_LTE(icell,lambda)/((wl**5)*(coeff_exp-1.0)) *delta_wl
     endif
  enddo
  E_dust = E_dust * 2.0*pi*hp*c_light**2

  ! emission étoile : BB seulement pour le moment
  E_etoile = 0.0
  cst=cst_th/etoile(1)%T
  do lambda=1, n_lambda
     ! longueur d'onde en metre
     wl = tab_lambda(lambda)*1.e-6
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     cst_wl=cst/wl
     if (cst_wl < 500.0) then
        coeff_exp=exp(cst_wl)
!        E_etoile = E_etoile + sum(kappa_abs_LTE(lambda,1,:,1)) /((wl**5)*(coeff_exp-1.0)) *delta_wl
        E_etoile = E_etoile + kappa_abs_LTE(icell,lambda) * spectre_etoiles(lambda) / ( 4*pi * AU_to_m**2)


       ! write(*,*)  2.0*pi*hp*c_light**2  * 4*pi*etoile(1)%r**2 * AU_to_m**2 / ((wl**5)*(coeff_exp-1.0)) * delta_wl  /  spectre_etoiles(lambda) !----> OK, c'est la bonne valeur de spectre etoile pour 1BB quand n_lambda est grand (binnage negligeable)
     endif
  enddo

  if (E_dust < tiny_real) then
     write(*,*) "Sublimation radius : something is wrong"
     write(*,*) "Opacity is not defined yet"
     write(*,*) "Maybe the parameter file is old ?"
     write(*,*) "Exiting"
     stop
  endif
  sublimation_radius = sqrt(E_etoile/E_dust)



  ! -------------
!---
!---  log_frac_E_abs=log(J_abs*n_phot_L_tot + E0(ri,zj,phik))
!---
!---
!---  Cst0 = 2.0*hp*c_light**2 * 1e-6 * pi*Rsun_to_AU**2  / (pc_to_AU**2)
!---
!---  ! pour le spectre de MCFOST
!---  tab_spectre(i,l) = Cst0/ ( ((exp(cst_w)) -1.)) * (wl**5))  ;
!---  terme = (surface / Cst0) * exp(interp(log_spectre, log_wl_spectre, log_lambda(lambda)))
!4*pi*(etoile(i)%r**2)
!---  terme = terme / N * (surface / Cst0)
!---
!---  spectre = spectre + terme * delta_wl

!---  spectre_etoiles(:) =  spectre_etoiles(:) * cst_spectre_etoiles
! 2.0*pi*hp*c_light**2 * (AU_to_m)**2
!---   surface=4*pi*(etoile(i)%r**2)
!---     cst_spectre_etoiles = 2.0*pi*hp*c_light**2 * (AU_to_m)**2 ! cst BB + r_etoile est en AU
!---  ---> spectre_etoile = B * surface * cst_spectre_etoiles * delta_wl
!---
  ! ---------------


  !  sublimation_radius =  (etoile(1)%T/T_max)**2 * etoile(1)%r
  write(*,*) "Optically thin sublimation radius =", real(sublimation_radius), "AU"

  sublimation_radius = sublimation_radius * 1.6

  do i=1,n_zones
     !write(*,*) "zone", i,sublimation_radius, disk_zone(i)%rmin, sublimation_radius > disk_zone(i)%rmin
     if (sublimation_radius > disk_zone(i)%rmin) then
        !write(*,*) "SUBLIMATING DUST IN ZONE", i
        disk_zone(i)%rmin = sublimation_radius
        disk_zone(i)%rin = disk_zone(i)%rmin !+ 5* disk_zone(1)%edge
        disk_zone(i)%edge = 0.0
     endif
  enddo !i
  rmin = minval(disk_zone(:)%rmin)
  write(*,*) "New minimum radius = ", rmin

  do i = 1, n_regions
     regions(i)%Rmin = max(regions(i)%Rmin,sublimation_radius)
  enddo !i

  return

end subroutine compute_othin_sublimation_radius

!***********************************************************

subroutine sublimate_dust()
  ! Supprime les grains dont la temperature depasse la temperature de sublimation
  ! C. Pinte
  ! 07/08/12

  integer :: i, j, pk, icell, k, ipop
  real :: mass

  write(*,*) "Sublimating dust"

  do i=1,n_rad
     do j=1,nz
        do pk=1,n_az
           icell = cell_map(i,j,pk)

           ! Cas LTE
           do k=1,n_grains_tot
              ipop = grain(k)%pop

              if (.not.dust_pop(ipop)%is_PAH) then
                 if (Temperature(icell) > dust_pop(ipop)%T_sub) then
                    densite_pouss(k,icell) = 0.0
                 endif
              endif
           enddo


        enddo !pk
     enddo !j
  enddo !i

  mass = 0.0
  do i=1,n_rad
     do j=1,nz
        icell = cell_map(i,j,1)
        do k=1,n_grains_tot
           mass=mass + densite_pouss(k,icell) * M_grain(k) * (volume(icell) * AU3_to_cm3)
        enddo
     enddo
  enddo
  mass =  mass/Msun_to_g

  write(*,*) 'New total dust mass in model :', mass,' Msun'

end subroutine sublimate_dust


!**********************************************************

subroutine equilibre_hydrostatique()
  ! Calcul l'equilibre hydrostatique pour chaque rayon
  ! Equation 2.4.3 de la these (page 38, 52 du pdf)
  ! Valable pour disque de gaz parfait, non-autogravitant, geometriquement mince
  !
  ! C. Pinte
  ! 25/09/07

  implicit none

  real, dimension(nz) :: rho, ln_rho
  real :: dz, dz_m1, dTdz, fac, fac1, fac2, M_etoiles, M_mol, somme, cst
  integer :: i,j, k, icell, icell_m1

  real, parameter :: gas_dust = 100

  M_etoiles = sum(etoile(:)%M) * Msun_to_kg
  M_mol = masse_mol_gaz * g_to_kg

  cst = Ggrav * M_etoiles * M_mol / (kb * AU_to_m**2)

  do k=1, n_az
     do i=1, n_rad
        ln_rho(1) = 0.
        rho(1) = 1.
        dz = delta_z(i)
        dz_m1 = 1.0/dz
        somme = rho(1)
        do j = 2, nz
           icell = cell_map(i,j,k)
           icell_m1 = cell_map(i,j-1,k)
           dTdz = (Temperature(icell)-Temperature(icell_m1)) * dz_m1
           fac1 = cst * z_grid(icell)/ (r_grid(icell)**3)
           fac2 = -1.0 * (dTdz + fac1) / Temperature(icell)
           ln_rho(j) = ln_rho(j-1) + fac2 * dz
           rho(j) = exp(ln_rho(j))
           somme = somme + rho(j)
        enddo !j

        ! Renormalisation
        do j = 1, nz
           icell = cell_map(i,j,k)
           fac = gas_dust * masse_rayon(i,k) / (volume(icell) * somme) ! TODO : densite est en particule, non ???
           densite_gaz(icell) =  rho(j) * fac
        enddo

     enddo !i
  enddo !k

  return

end subroutine equilibre_hydrostatique

!**********************************************************

end module disk_physics
