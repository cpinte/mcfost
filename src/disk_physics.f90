module disk_physics

  use grains
  use mcfost_env
  use opacity
  use em_th
  use molecular_emission, only : densite_gaz
  use constantes
  use stars, only : spectre_etoiles
  use messages
  use wavelengths

  implicit none

  character(len=32) :: sublimationFile = "sublimation_radius.txt"

contains


subroutine compute_othin_sublimation_radius()
  ! Dans le cas optiquement mince, ne dépend que de la température (et spectre) de l'étoile

  implicit none

  real(kind=dp) :: E_dust, E_etoile, coeff_exp, cst_wl, sublimation_radius
  real :: cst, wl, delta_wl
  integer :: lambda, icell, i

  E_dust = 0.0
  cst=cst_th/dust_pop(1)%T_sub
  cst=cst_th/1500.

  icell = icell_ref

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

  ! Emission étoiles
  do i=1, n_etoiles
     E_etoile = 0.0
     do lambda=1, n_lambda
        E_etoile = E_etoile + kappa_abs_LTE(icell,lambda) * spectre_etoiles(lambda) / ( 4*pi * AU_to_m**2)
     enddo

     if (E_dust < tiny_real) then
        call error("Sublimation radius : opacity is not defined yet", &
             msg2="Maybe the parameter file is old ?")
     endif
     etoile(i)%othin_sublimation_radius = sqrt(E_etoile/E_dust)
  enddo

  if (.not.lVoronoi) then
     sublimation_radius = real(etoile(1)%othin_sublimation_radius)
     write(*,*) "Optically thin sublimation radius =", real(sublimation_radius), "AU"
     sublimation_radius = sublimation_radius * 1.6

     open(unit=1,file=trim(data_dir)//"/"//trim(sublimationFile),status="replace")
     write(1,*) etoile(:)%othin_sublimation_radius
     close(1)

     call set_sublimation_radius(sublimation_radius)
  endif

  return

end subroutine compute_othin_sublimation_radius

!***********************************************************

subroutine set_sublimation_radius(sublimation_radius)

  real(kind=dp), intent(in) :: sublimation_radius

  integer :: i

  do i=1,n_zones
     !write(*,*) "zone", i,sublimation_radius, disk_zone(i)%rmin, sublimation_radius > disk_zone(i)%rmin
     if (sublimation_radius < disk_zone(i)%rmin) then
        !write(*,*) "SUBLIMATING DUST IN ZONE", i
        disk_zone(i)%rmin = sublimation_radius
        disk_zone(i)%rin = disk_zone(i)%rmin + 5* disk_zone(1)%edge
        disk_zone(i)%edge = 0.0
     endif
  enddo !i
  rmin = minval(disk_zone(:)%rmin)
  write(*,*) "New minimum radius = ", rmin, "AU"

  do i = 1, n_regions
     regions(i)%Rmin = max(regions(i)%Rmin,sublimation_radius)
  enddo !i

end subroutine set_sublimation_radius

!***********************************************************

subroutine read_sublimation_radius()

  real(kind=dp) :: sublimation_radius

  write(*,*) "Reading sublimation file : ./data_th/"//trim(sublimationfile)

  open(unit=1,file=trim(root_dir)//"/"//trim(seed_dir)//"/data_th/"//trim(sublimationfile), status="old")
  read(1,*) sublimation_radius
  close(unit=1)

  call set_sublimation_radius(sublimation_radius)

  return

end subroutine read_sublimation_radius

!**********************************************************************

subroutine sublimate_dust()
  ! Supprime les grains dont la temperature depasse la temperature de sublimation
  ! C. Pinte
  ! 07/08/12

  integer :: i, j, pk, icell, k, ipop, l
  real :: mass

  write(*,*) "Sublimating dust"

  do i=1,n_rad
     do j=j_start,nz
        if (j==0) cycle
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
     do j=j_start,nz
        if (j==0) cycle
        do k=1,n_az
           icell = cell_map(i,j,k)
           do l=1,n_grains_tot
              mass=mass + densite_pouss(l,icell) * M_grain(l) * (volume(icell) * AU3_to_cm3)
           enddo
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
