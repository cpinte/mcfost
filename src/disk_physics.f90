module disk_physics

  use grains
  use opacity
  use em_th
  use prop_star
  use constantes

contains
  

subroutine compute_othin_sublimation_radius()
  ! Dans le cas optiquement mince, ne dépend que de la température de l'étoile
  
  implicit none

  real(kind=db) :: E_dust, E_etoile
  real :: cst, cst_wl, coeff_exp, wl, delta_wl
  integer :: lambda, k, i
  real(kind=db) :: sublimation_radius, Rsub, Hsub

  if (n_etoiles > 1) then
     write(*,*) "Sublimation radius calculation not implemented"
     write(*,*) "in case of several stars"
     stop
  endif

  if (n_zones > 1) then
     write(*,*) "Sublimation radius calculation not implemented"
     write(*,*) "in case of several zones"
     stop
  endif
     

  !  if (lstrat) then
  !     write(*,*) "Sublimation radius calculation not implemented"
  !     write(*,*) "in case of statification"
  !     stop
  !  endif
  ! --> kappa moyenne verticalement

  !  Emission poussière
  E_dust = 0.0
  cst=cst_th/T_max
  do lambda=1, n_lambda
     ! longueur d'onde en metre
     wl = tab_lambda(lambda)*1.e-6
     delta_wl=tab_delta_lambda(lambda)*1.e-6
     cst_wl=cst/wl
     if (cst_wl < 500.0) then
        coeff_exp=exp(cst_wl)
        E_dust = E_dust + 4.0*sum(kappa_abs_eg(lambda,1,:,1))/((wl**5)*(coeff_exp-1.0))*delta_wl 
     endif
  enddo

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
        E_etoile = E_etoile + sum(kappa_abs_eg(lambda,1,:,1))/((wl**5)*(coeff_exp-1.0))*delta_wl 
     endif
  enddo

  ! sublimation_radius = etoile(1)%r * sqrt(E_etoile/E_dust) * sqrt(1 + disk_zone(1)%sclht / disk_zone(1)%rref)
 
  !Rsub = disk_zone(1)%rmin
  !do i=1,10 ! pour iteration, converge tres vite
  !   Hsub = disk_zone(1)%sclht * (Rsub/disk_zone(1)%rref)**disk_zone(1)%exp_beta    
  !   Rsub = etoile(1)%r * sqrt(E_etoile/E_dust) * sqrt(1 +  min(Hsub/Rsub,1.0_db))
  !enddo

  Rsub = etoile(1)%r * sqrt(E_etoile/E_dust) * sqrt(2.)

  sublimation_radius = Rsub
  
  !  sublimation_radius =  (etoile(1)%T/T_max)**2 * etoile(1)%r
  write(*,*) "Optically thin sublimation radius =", real(sublimation_radius), "AU"

  do i=1,n_zones
     if (sublimation_radius > disk_zone(i)%rmin) then
        disk_zone(i)%rmin = sublimation_radius
        disk_zone(i)%rin = disk_zone(i)%rmin !+ 5* disk_zone(1)%edge
        disk_zone(i)%edge = 0.0
     endif
  enddo !i
  rmin = minval(disk_zone(:)%rmin)

  do i = 1, n_regions
     Rmin_region(i) = max(Rmin_region(i),sublimation_radius)
  enddo !i
  
  return

end subroutine compute_othin_sublimation_radius

!***********************************************************

subroutine sublimate_dust()
  ! Supprime les grains dont la temperature depasse la temperature de sublimation
  ! C. Pinte
  ! 07/08/12

  integer :: i, j, pk, k, ipop
  real :: mass

  write(*,*) "Sublimating dust"
    
  do i=1,n_rad
     do j=1,nz
        do pk=1,n_az
           
           ! Cas LTE
           do k=1,n_grains_tot
              ipop = grain(k)%pop
  
              if (.not.dust_pop(ipop)%is_PAH) then
                 if (Temperature(i,j,pk) > dust_pop(ipop)%T_sub) then
                    densite_pouss(i,j,pk,k) = 0.0
                 endif
              endif
           enddo
           
           
        enddo !pk
     enddo !j
  enddo !i
  
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
    
end subroutine sublimate_dust

end module disk_physics
