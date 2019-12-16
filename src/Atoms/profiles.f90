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
 !real(kind=dp), dimension(:,:,:) :: line_profiles

 CONTAINS
 
!if projection done before, we do not need x,y,z,l ect

 SUBROUTINE Iprofile (line,icell,x,y,z,x1,y1,z1,u,v,w,l,id, Nvspace, Omegav)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
  integer, intent(in) 							            :: icell, id, Nvspace
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  real(kind=dp), intent(in) 								:: Omegav(:)
  type (AtomicLine), intent(inout)								:: line
  real(kind=dp), dimension(line%Nlambda)					:: vvoigt
  integer													::  Nred, Nblue, i, j, nv
  real(kind=dp)												:: norm


  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue
  
  line%phi_loc(:,id) = 0d0
  !line_profiles(:,kr,id) = 0d0
  
  norm =  SQRTPI * line%atom%vbroad(icell) * Nvspace

  if (line%voigt) then
  !Now we have a pointer to atom in line. atom(n)%lines(kr)%atom => atom(n)
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
            !                 2) Voronoi grid is used                 
                        
          vvoigt(:) = ( line%u(:) - omegav(nv) ) / line%atom%vbroad(icell)


!           write(*,*) nv, "0loc=", locate(vvoigt, 0d0)
!           write(*,*) maxval(vvoigt), vbroad
          line%phi_loc(:,id) = line%phi_loc(:,id) + &
          						Voigt(line%Nlambda, line%a(icell),vvoigt(:)) !/ Nvspace

      end do
  else !Gaussian !only for checking
      do nv=1, Nvspace

         vvoigt(:) = ( line%u(:) - omegav(nv) ) / line%atom%vbroad(icell)
         line%phi_loc(:,id) = line%phi_loc(:,id) + dexp(-(vvoigt(:))**2) !/ Nvspace 

      end do
 end if !line%voigt
 line%phi_loc(:,id) = line%phi_loc(:,id) / norm !/ (SQRTPI * line%atom%vbroad(icell))

!  if (any_nan_infinity_vector(line%phi_loc(:,id))>0 .or. minval(line%phi_loc(:,id)) < 0) then
!   write(*,*) line%Nlambda, icell, id
!   if (line%voigt) write(*,*) "Damping = ", line%a(icell)
!   write(*,*) "vv=",vvoigt !unitless
!   write(*,*) "vv2(km/s)=", vvoigt(:) * line%atom%vbroad(icell)*1d-3
!   write(*,*) " Error with Profile"
!   write(*,*) line%phi_loc(:,id)
!   stop
!  end if

 RETURN
 END SUBROUTINE IProfile
 
 !-> TO be included as a Voigt procedure
 SUBROUTINE Iprofile_thomson (line,icell,x,y,z,x1,y1,z1,u,v,w,l,id, Nvspace, Omegav)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
  integer, intent(in) 							            :: icell, id, Nvspace
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  real(kind=dp), intent(in) 								:: Omegav(:)
  type (AtomicLine), intent(inout)								:: line
  real(kind=dp), dimension(line%Nlambda)					:: vvoigt
  integer													::  Nred, Nblue, i, j, nv

  nv = 1 !=Nvspace

  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue
  
  line%phi_loc(:,id) = 0d0
  !line_profiles(:,kr,id) = 0d0
  
!   r = line%atom%vbroad(icell)*line%a(icell)/line%aeff(icell)
!   eta = 1.36603*r - 0.47719*r*r + 0.11116*r*r*r

  !normed
!   line%phi_loc(:,id) = eta*(line%aeff(icell) / ((line%u(:)-omegav(nv))**2.+line%aeff(icell)**2.) / pi) &
!                        + (1-eta)*dexp(-(vvoigt(:)/line%aeff(icell))**2)/SQRTPI/line%aeff(icell)
                       
  line%phi_loc(:,id) = line%r(icell)*line%aeff(icell) / ((line%u(:)-omegav(nv))**2.+line%aeff(icell)**2.) &
                       + line%r1(icell)*dexp(-(line%u(:)/line%aeff(icell))**2)


 RETURN
 END SUBROUTINE IProfile_thomson

 !interpolation or shifting, building
 SUBROUTINE IProfile_cmf_to_obs(line,icell,x,y,z,x1,y1,z1,u,v,w,l, id, Nvspace, Omegav)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
  integer, intent(in) 							            :: icell, id, Nvspace
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  real(kind=dp), intent(in) 								:: Omegav(:)
  type (AtomicLine), intent(inout)								:: line
  real(kind=dp), dimension(line%Nlambda)					:: vvoigt
  integer													::  Nred, Nblue, i, j, nv
  real(kind=dp), dimension(line%Nlambda)                    :: u1, u1p
 

  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue

  !temporary here
  line%phi_loc(:,id) = 0d0
 
  u1(:) = line%u(:)/line%atom%vbroad(icell)

 do nv=1, Nvspace 
 
         u1p(:) = u1(:) - omegav(nv)/line%atom%vbroad(icell)
             
         line%phi_loc(:,id) = line%phi_loc(:,id) + &
                 linear_1D_sorted(line%Nlambda,u1,line%phi(:,icell),line%Nlambda,u1p) / Nvspace /sqrtpi / line%atom%vbroad(icell)
         	

 enddo
 


 RETURN
 END SUBROUTINE IProfile_cmf_to_obs
 
 SUBROUTINE ZProfile (line, icell,x,y,z,x1,y1,z1,u,v,w,l,id, Nvspace, Omegav)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
!->>>>>>
integer :: iray = 1 !futur deprecation
!<<<<<<-
  integer, intent(in) 							            :: icell, id, Nvspace
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  real(kind=dp), intent(in) 								:: Omegav(:)
  type (AtomicLine), intent(inout)								:: line
  real(kind=dp), dimension(line%Nlambda)                 :: vvoigt, F, LV, vvoigt_b
  integer, parameter										:: NbspaceMax=101
  real(kind=dp), dimension(NbspaceMax)					:: omegaB, gamma, chi
  integer													:: nv, Nred, Nblue, nc, &
  															   Nbspace, nb, Nzc, i, j,qz
  real(kind=dp) 											:: b0, b1,g1,c1,dB,norm

  real(kind=dp), dimension(3,line%Nlambda) 				:: phi_zc, psi_zc!Sigma_b, PI, sigma_r
  logical 													:: B_flag = .true.
  !or allocate deallocate only on Nlambda. Lower arrays dimension but took time to allocate

  omegaB = 0d0
  
  CALL error("Do the  projection of B in metal_bb or NLTE_bound_bound line v_proj")

  b0 = B_project(icell,x,y,z,u,v,w,g1,c1)
  omegaB(1) = b0; Nbspace = 1
  gamma(1) = g1; chi(1)=c1
  Nbspace = 1



!   if (.not.lvoronoi .and. B_flag) then
!       b1 = B_project(icell,x1,y1,z1,u,v,w,g1,c1)
!       !Nbspace = NbspaceMax
!       dB = (b1-b0)
!       dB = dabs(dB * line%g_lande_eff) * LARMOR * line%lambda0 * NM_TO_M
!       Nbspace = max(2,nint(20*dB/vbroad))
!       Nbspace = min(Nbspace,NbspaceMax)
!       !write(*,*) Nbspace, b1*1e4, b0*1e4, dB
!       omegaB(Nbspace) = b1
!       gamma(Nbspace) = g1; chi(Nbspace)=c1
!       do nv=2,Nbspace-1
!        delta_vol_phi = (real(nv,kind=dp))/(real(Nbspace,kind=dp)) * l
!        xphi=x+delta_vol_phi*u
!        yphi=y+delta_vol_phi*v
!        zphi=z+delta_vol_phi*w
!        omegaB(nv) = B_project(icell,xphi,yphi,zphi,u,v,w,g1,c1)
!        gamma(nv) = g1; chi(nv)=c1
!       end do      
!   end if

  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue
  line%phi_loc(:,id) = 0d0
  
  !allocate(vv(line%Nlambda), vvoigt(line%Nlambda))


  Nzc = line%zm%Ncomponent
  if (.not.line%voigt) then !unpolarised line assumed even if line%polarizable
     norm  = line%atom%vbroad(icell) * Nvspace * SQRTPI

      do nv=1, Nvspace
      
         vvoigt(:) = (line%u - omegav(nv)) / line%atom%vbroad(icell)
         line%phi_loc(:,id) = line%phi_loc(:,id) + dexp(-(vvoigt(:))**2) !/ Nvspace

      !derivative of Gaussian:
!          F(Nblue:Nred) = F(Nblue:Nred) - &
!            2d0 * dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace * &
!            NLTEspec%lambda(Nblue:Nred) * CLIGHT / (line%atom%vbroad(icell) * line%lambda0)

      end do
      line%phi_loc(:,id) = line%phi_loc(:,id) / norm !/ (SQRTPI * vbroad)
!       F(Nblue:Nred) = F(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))
      !deallocate(vv, vvoigt)
      RETURN
  end if
  norm  = line%atom%vbroad(icell) * Nvspace * Nbspace * SQRTPI

  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)

  !allocate(LV(line%Nlambda), F(line%Nlambda), psi_zc(3,line%Nlambda),phi_zc(3,line%Nlambda))
  LV = 0d0
  F = 0d0
  psi_zc(:,:) = 0d0; phi_zc(:,:) = 0d0
  line%psi(:,:,iray) = 0d0; line%phiZ(:,iray,id) = 0d0
  !Should work also for unpolarised voigt line because Ncz=1,S=0,q=0,shift=0
       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
      
        vvoigt(:) = (line%u - omegav(nv)) / line%atom%vbroad(icell)

         do nb=1,Nbspace !Nbspace=1 if Voronoi, or magnetic field is 0d0.But present

          do nc=1,Nzc
             ! the splitting is 0 if unpolarized 'cause zm%shift(nc=Nzc=1)=0d0
             !there is a + omegaB because, -deltaLam^JJp_MMp=splitting = lamB * (gp*Mp - g*M)
             vvoigt_b(:) = vvoigt(:) + line%zm%shift(nc) * omegaB(nb) * &
                               LARMOR * line%lambda0 * NM_TO_M / line%atom%vbroad(icell)
             !The normalisation by Nzc is done when compute the strength*profiles.
             !In case of unpolarised line, Nzc = 1 and there is no need to normalised
             LV(:) = Voigt(line%Nlambda, line%adamp,vvoigt_b, F)
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
              		line%zm%strength(nc) * F(:) / norm !/ (SQRTPI * vbroad * Nbspace*Nvspace)
             phi_zc(qz,:) = phi_zc(qz,:) + &
              		line%zm%strength(nc) * LV(:) / norm !/ (SQRTPI * vbroad * Nbspace*Nvspace)
             LV(:) = 0d0; F(:) = 0d0
             !write(*,*) Nzc, qz, line%zm%q(nc), vvoigt_b(line%Nlambda/2+1), line%zm%shift(nc), line%zm%strength(nc), omegab(nb)
          end do !components
          !the output, for the other we store chi_pol/chi_I, rho_pol/chi_I etc
          !write(*,*) dsin(gamma(nb)), dcos(gamma(nb)), dsin(2*chi(nb)), dcos(2*chi(nb))
          line%phi_loc(:,id) = line%phi_loc(:,id) + 5d-1 *(phi_zc(2,:) * dsin(gamma(nb))*dsin(gamma(nb)) + \
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
          
!! Should be wore like psi(Nlambda, ncompo, id) because we do not care about keeping it for every direction

          line%psi(1,:,iray) = line%psi(1,:,iray) + &
          		0.5*(psi_zc(2,:)-0.5*(psi_zc(1,:)+psi_zc(3,:)))*cos(2*chi(nb))*sin(gamma(nb))**2	

          !rhoU/chiI
          line%psi(2,:,iray) = line%psi(2,:,iray) + &
          		0.5*(psi_zc(2,:)-0.5*(psi_zc(1,:)+psi_zc(3,:)))*sin(2*chi(nb))*sin(gamma(nb))**2
          !rhoV/chiI
          line%psi(3,:,iray) = line%psi(3,:,iray) + 0.5*(psi_zc(3,:)-psi_zc(1,:))*cos(gamma(nb))
          !chiQ/chiI
          line%phiZ(1,:,iray) = line%phiZ(1,:,iray) + &
          		0.5*(phi_zc(2,:)-0.5*(phi_zc(1,:)+phi_zc(3,:)))*cos(2*chi(nb))*sin(gamma(nb))**2
          !chiU/chiI
          line%phiZ(2,:,iray) = line%phiZ(2,:,iray) + &
          		0.5*(phi_zc(2,:)-0.5*(phi_zc(1,:)+phi_zc(3,:)))*sin(2*chi(nb))*sin(gamma(nb))**2
          !chiV/chiI
          line%phiZ(3,:,iray) = line%phiZ(3,:,iray) + 0.5*(phi_zc(3,:)-phi_zc(1,:))*cos(gamma(nb))
        end do !magnetic field     

       end do !velocity
       !check that if line is not polarised the Zeeman components are 0
       !write(*,*) "tpt", allocated(psi_zc),allocated(phi_zc), allocated(LV), allocated(F), &
       ! allocated(vv), allocated(vvoigt)
  !deallocate(psi_zc, phi_zc, LV, F, vv, vvoigt)
 RETURN
 END SUBROUTINE ZProfile


END MODULE PROFILES
