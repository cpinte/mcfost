MODULE PROFILES

 use atmos_type, only                : atmos, B_project, VBROAD_atom
 use constant
 use atom_type
 use spectrum_type, only			 : NLTEspec
 use voigtfunctions, only 			 : Voigt
 use broad, only 					 : Damping
 use math

 ! MCFOST's original
 use mcfost_env, only : dp
 use molecular_emission, only		 : v_proj
 use parametres
 use input
 use constantes, only				 : tiny_dp, huge_dp

 IMPLICIT NONE

 PROCEDURE(Iprofile), pointer :: Profile => null()

 CONTAINS

 !TO DO: add the Dimension N of the output profile, could be usefull
 !if we do wavelength by wavelength opacities.
 SUBROUTINE Iprofile (line, icell,x,y,z,x1,y1,z1,u,v,w,l, P, phi, psi)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
  integer, intent(in) 							            :: icell
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  type (AtomicLine), intent(in)								:: line
  double precision, dimension(line%Nlambda)					:: vvoigt, F, vv
  integer, parameter										:: NvspaceMax = 101
  double precision, dimension(NvspaceMax)					:: omegav
  integer													:: Nvspace, nv, Nred, Nblue, i, j
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv, vbroad
  double precision, intent(out), dimension(:)               :: P
  double precision, intent(out), dimension(:,:), optional   :: phi, psi

  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   v0 = v_proj(icell,x,y,z,u,v,w)
   omegav(1) = v0
  end if
  
  vbroad = VBROAD_atom(icell,line%atom)
   
  if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0)
     Nvspace = max(2,nint(20*dv/vbroad))
     Nvspace = min(Nvspace,NvspaceMax)
     omegav(Nvspace) = v1
    do nv=2,Nvspace-1
      delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
      xphi=x+delta_vol_phi*u
      yphi=y+delta_vol_phi*v
      zphi=z+delta_vol_phi*w
      omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
    end do 
  end if


  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue

  P = 0d0
  !allocate(vv(line%Nlambda), F(line%Nlambda), vvoigt(line%Nlambda))
  vv = 0d0
  vv(:) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * vbroad)


  if (line%voigt) then
  !Now we have a pointer to atom in line. atom(n)%lines(kr)%atom => atom(n) 
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
          vvoigt(:) = vv(:) - omegav(nv) / vbroad

          P(:) = P(:) + &
             Voigt(line%Nlambda, line%adamp,vvoigt(:), F, VoigtMethod) / Nvspace

      end do
  else !Gaussian !only for checking
      do nv=1, Nvspace
      
         vvoigt(:) = vv(:) - omegav(nv) / vbroad
         P(:) = P(:) + dexp(-(vvoigt(:))**2) / Nvspace 
         
      end do
 end if !line%voigt
 P(:) = P(:) / (SQRTPI * vbroad)
 if (minval(P) < 0) CALL error("P should not be negative")

 !deallocate(vv, vvoigt, F)
 RETURN
 END SUBROUTINE IProfile
 
 SUBROUTINE ZProfile (line, icell,x,y,z,x1,y1,z1,u,v,w,l, P, phi, psi)
  integer, intent(in) 							            :: icell
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  type (AtomicLine), intent(in)								:: line
  double precision, dimension(line%Nlambda)                 :: vvoigt, vv, F, LV, vvoigt_b
  integer, parameter										:: NvspaceMax = 101, NbspaceMax=101
  double precision, dimension(NvspaceMax)					:: omegav
  double precision, dimension(NbspaceMax)					:: omegaB, gamma, chi
  integer													:: Nvspace, nv, Nred, Nblue, nc, &
  															   Nbspace, nb, Nzc, i, j,qz
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv, b0, b1,g1,c1,dB,vbroad
  double precision, intent(out), dimension(:)               :: P
  double precision, intent(out), dimension(:,:) 		    :: phi, psi !eta_QUV; rho_QUV
  double precision, dimension(3,line%Nlambda) 				:: phi_zc, psi_zc!Sigma_b, PI, sigma_r
  logical 													:: B_flag = .true.
  !or allocate deallocate only on Nlambda. Lower arrays dimension but took time to allocate

  omegaB = 0d0
  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   v0 = v_proj(icell,x,y,z,u,v,w)
   omegav(1) = v0
  end if
  
  vbroad = VBROAD_atom(icell,line%atom)

  b0 = B_project(icell,x,y,z,u,v,w,g1,c1)
  omegaB(1) = b0; Nbspace = 1
  gamma(1) = g1; chi(1)=c1

  if (maxval(abs(atmos%Bxyz(icell,:))) == 0d0) B_flag = .false.
   
  if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0) 
     Nvspace = max(2,nint(20*dv/vbroad))
     Nvspace = min(Nvspace,NvspaceMax)
     omegav(Nvspace) = v1
    do nv=2,Nvspace-1
      delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
      xphi=x+delta_vol_phi*u
      yphi=y+delta_vol_phi*v
      zphi=z+delta_vol_phi*w
      omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
    end do 
  end if


  if (.not.lvoronoi .and. B_flag) then
      b1 = B_project(icell,x1,y1,z1,u,v,w,g1,c1)
      !Nbspace = NbspaceMax
      dB = (b1-b0)
      dB = dabs(dB * line%g_lande_eff) * LARMOR * line%lambda0 * NM_TO_M
      Nbspace = max(2,nint(20*dB/vbroad))
      Nbspace = min(Nbspace,NbspaceMax)
      !write(*,*) Nbspace, b1*1e4, b0*1e4, dB
      omegaB(Nbspace) = b1
      gamma(Nbspace) = g1; chi(Nbspace)=c1
      do nv=2,Nbspace-1
       delta_vol_phi = (real(nv,kind=dp))/(real(Nbspace,kind=dp)) * l
       xphi=x+delta_vol_phi*u
       yphi=y+delta_vol_phi*v
       zphi=z+delta_vol_phi*w
       omegaB(nv) = B_project(icell,xphi,yphi,zphi,u,v,w,g1,c1)
       gamma(nv) = g1; chi(nv)=c1
      end do      
  end if

  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue
  P = 0d0
  
  !allocate(vv(line%Nlambda), vvoigt(line%Nlambda))

  vv(:) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * vbroad)

  Nzc = line%zm%Ncomponent
  if (.not.line%voigt) then !unpolarised line assumed even if line%polarizable
      do nv=1, Nvspace
      
         vvoigt(:) = vv(:) - omegav(nv) / vbroad
         P(:) = P(:) + dexp(-(vvoigt(:))**2) / Nvspace
      !derivative of Gaussian:
!          F(Nblue:Nred) = F(Nblue:Nred) - &
!            2d0 * dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace * &
!            NLTEspec%lambda(Nblue:Nred) * CLIGHT / (line%atom%vbroad(icell) * line%lambda0)

      end do
      P(:) = P(:) / (SQRTPI * vbroad)
!       F(Nblue:Nred) = F(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))
      !deallocate(vv, vvoigt)
      RETURN
  end if
     
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)

  !allocate(LV(line%Nlambda), F(line%Nlambda), psi_zc(3,line%Nlambda),phi_zc(3,line%Nlambda))
  LV = 0d0
  F = 0d0
  psi_zc(:,:) = 0d0; phi_zc(:,:) = 0d0
  psi(:,:) = 0d0; phi(:,:) = 0d0
  !Should work also for unpolarised voigt line because Ncz=1,S=0,q=0,shift=0
       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
        vvoigt(:) = vv(:) - omegav(nv) / vbroad
        
         do nb=1,Nbspace !Nbspace=1 if Voronoi, or magnetic field is 0d0.But present

          do nc=1,Nzc
             ! the splitting is 0 if unpolarized 'cause zm%shift(nc=Nzc=1)=0d0
             !there is a + omegaB because, -deltaLam^JJp_MMp=splitting = lamB * (gp*Mp - g*M)
             vvoigt_b(:) = vvoigt(:) + line%zm%shift(nc) * omegaB(nb) * &
                               LARMOR * line%lambda0 * NM_TO_M / vbroad
             !The normalisation by Nzc is done when compute the strength*profiles.
             !In case of unpolarised line, Nzc = 1 and there is no need to normalised
             LV(:) = Voigt(line%Nlambda, line%adamp,vvoigt_b,F, VoigtMethod)
             !qz = 1 if line%zm%q = - 1, 2 if line%zm%q=0, 3 if line%zm%q = 1
             !with negative index, qz = -line%zm%q
             SELECT CASE (line%zm%q(nc))
              CASE (-1)
               qz = 1
              CASE (+0)
               qz = 2
              CASE (+1)
               qz = 3
              CASE DEFAULT
               CALL ERROR("line%zm%q should be (/-1,0,1/)!")
             END SELECT
             !!if (abs(line%zm%q(nc)) > 1) CALL Error("BUG") !in the SELECT CASE
             psi_zc(qz,:) = psi_zc(qz,:) + &
              		line%zm%strength(nc) * F(:) / (SQRTPI * vbroad * Nbspace*Nvspace)
             phi_zc(qz,:) = phi_zc(qz,:) + &
              		line%zm%strength(nc) * LV(:) / (SQRTPI * vbroad * Nbspace*Nvspace)
             LV(:) = 0d0; F(:) = 0d0
             !write(*,*) Nzc, qz, line%zm%q(nc), vvoigt_b(line%Nlambda/2+1), line%zm%shift(nc), line%zm%strength(nc), omegab(nb)
          end do !components 
          !the output, for the other we store chi_pol/chi_I, rho_pol/chi_I etc
          !write(*,*) dsin(gamma(nb)), dcos(gamma(nb)), dsin(2*chi(nb)), dcos(2*chi(nb))
          P(:) = P(:) + 5d-1 *(phi_zc(2,:) * dsin(gamma(nb))*dsin(gamma(nb)) + \
            5d-1 *(1d0+dcos(gamma(nb))*dcos(gamma(nb))) * (phi_zc(1,:)+phi_zc(3,:))) ! profile in chiI, etaI
          !rhoQ/chiI
!           psi(1,:) = psi(1,:) + &
!           		0.5*(real(psi_zc(2,:))-0.5*real(psi_zc(1,:)+psi_zc(3,:)))*cos(2*chi(nb))*sin(gamma(nb))**2	
!           !rhoU/chiI
!           psi(2,:) = psi(2,:) + &
!           		0.5*(real(psi_zc(2,:))-0.5*real(psi_zc(1,:)+psi_zc(3,:)))*sin(2*chi(nb))*sin(gamma(nb))**2
!           !rhoV/chiI
!           psi(3,:) = psi(3,:) + 0.5*(real(psi_zc(3,:))-real(psi_zc(1,:)))*cos(gamma(nb))
!           !chiQ/chiI
!           phi(1,:) = phi(1,:) + &
!           		0.5*(real(phi_zc(2,:))-0.5*real(phi_zc(1,:)+phi_zc(3,:)))*cos(2*chi(nb))*sin(gamma(nb))**2
!           !chiU/chiI
!           phi(2,:) = phi(2,:) + &
!           		0.5*(real(phi_zc(2,:))-0.5*real(phi_zc(1,:)+phi_zc(3,:)))*sin(2*chi(nb))*sin(gamma(nb))**2
!           !chiV/chiI
!           phi(3,:) = phi(3,:) + 0.5*(real(phi_zc(3,:))-real(phi_zc(1,:)))*cos(gamma(nb))
          !rhoQ/chiI
          psi(1,:) = psi(1,:) + &
          		0.5*(psi_zc(2,:)-0.5*(psi_zc(1,:)+psi_zc(3,:)))*cos(2*chi(nb))*sin(gamma(nb))**2	
          !rhoU/chiI
          psi(2,:) = psi(2,:) + &
          		0.5*(psi_zc(2,:)-0.5*(psi_zc(1,:)+psi_zc(3,:)))*sin(2*chi(nb))*sin(gamma(nb))**2
          !rhoV/chiI
          psi(3,:) = psi(3,:) + 0.5*(psi_zc(3,:)-psi_zc(1,:))*cos(gamma(nb))
          !chiQ/chiI
          phi(1,:) = phi(1,:) + &
          		0.5*(phi_zc(2,:)-0.5*(phi_zc(1,:)+phi_zc(3,:)))*cos(2*chi(nb))*sin(gamma(nb))**2
          !chiU/chiI
          phi(2,:) = phi(2,:) + &
          		0.5*(phi_zc(2,:)-0.5*(phi_zc(1,:)+phi_zc(3,:)))*sin(2*chi(nb))*sin(gamma(nb))**2
          !chiV/chiI
          phi(3,:) = phi(3,:) + 0.5*(phi_zc(3,:)-phi_zc(1,:))*cos(gamma(nb))
        end do !magnetic field     
        
       end do !velocity
       !check that if line is not polarised the Zeeman components are 0
       !write(*,*) "tpt", allocated(psi_zc),allocated(phi_zc), allocated(LV), allocated(F), &
       ! allocated(vv), allocated(vvoigt)
  !deallocate(psi_zc, phi_zc, LV, F, vv, vvoigt) 
 RETURN
 END SUBROUTINE ZProfile


END MODULE PROFILES
