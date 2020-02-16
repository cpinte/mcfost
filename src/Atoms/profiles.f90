MODULE PROFILES

	use atmos_type, only				: atmos, B_project, VBROAD_atom
	use constant
	use atom_type
	use spectrum_type, only				: NLTEspec
	use voigtfunctions, only			: Voigt
	use broad, only						: Damping
	use math
	use getlambda, only 				: hv

 ! MCFOST's original
	use mcfost_env, only				: dp
	use molecular_emission, only		: v_proj
	use parametres
	use input
	use constantes, only				: tiny_dp, huge_dp

	IMPLICIT NONE

	PROCEDURE(Iprofile), pointer :: Profile => null()

	CONTAINS
 
	SUBROUTINE compute_shift_index(id, icell, iray, x, y, z, x1, y1, z1, u, v, w, l)!, Nvspace, Omegav)
		integer, intent(in) :: id, icell, iray
		real(kind=dp), intent(in) :: x, y, z, x1, y1, z1, u, v, w, l
		!integer, optional, intent(in) :: Nvspace
		!real(kind=dp), intent(inout), optional :: Omegav(Nvspace)
		integer :: kr, nact
		real(kind=dp) :: v0
		
! 		if (present(Nvsapce)) then
! 		
! 			....
! 		
! 		endif

		v0 = v_proj(icell,(x+x1)*0.5,(y+y1)*0.5,(z+z1)*0.5,u,v,w)


		do nact=1,atmos%NactiveAtoms
		
			do kr=1,atmos%ActiveAtoms(nact)%ptr_atom%Nline
			
				atmos%ActiveAtoms(nact)%ptr_atom%lines(kr)%dk(iray,id) = nint(1e-3 * v0 / hv)
			
			enddo
		
		enddo
		
	
	
	RETURN
	END SUBROUTINE compute_shift_index
 
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
  
  line%phi(:,id) = 0d0
  !line_profiles(:,kr,id) = 0d0
  
  norm =  SQRTPI * line%atom%vbroad(icell) * Nvspace
  
  !project magnetic field if any
  

!   if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
!      v1 = v_proj(icell,x1,y1,z1,u,v,w)
!      dv = dabs(v1-v0)
!      Nvspace = max(2,nint(20*dv/line%atom%vbroad(icell)))
!      Nvspace = min(Nvspace,NvspaceMax)
!      omegav(Nvspace) = v1
!     do nv=2,Nvspace-1
!       delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
!       xphi=x+delta_vol_phi*u
!       yphi=y+delta_vol_phi*v
!       zphi=z+delta_vol_phi*w
!       omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
!     end do 
!   end if
! write(*,*) "Nvspace=", Nvspace

  if (line%voigt) then
  !Now we have a pointer to atom in line. atom(n)%lines(kr)%atom => atom(n)
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields
            !                 2) Voronoi grid is used                 
                        
          vvoigt(:) = ( line%u(:) - omegav(nv) ) / line%atom%vbroad(icell)


!           write(*,*) nv, "0loc=", locate(vvoigt, 0d0)
!           write(*,*) maxval(vvoigt), vbroad
          line%phi(:,id) = line%phi(:,id) + &
          						Voigt(line%Nlambda, line%a(icell),vvoigt(:)) !/ Nvspace

      end do
  else !Gaussian !only for checking
      do nv=1, Nvspace

         vvoigt(:) = ( line%u(:) - omegav(nv) ) / line%atom%vbroad(icell)
         line%phi(:,id) = line%phi(:,id) + dexp(-(vvoigt(:))**2) !/ Nvspace 

      end do
 end if !line%voigt
 line%phi(:,id) = line%phi(:,id) / norm !/ (SQRTPI * line%atom%vbroad(icell))


 RETURN
 END SUBROUTINE IProfile
 
 !-> TO be included as a Voigt procedure
 !Missing velocity shift here
 !building
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
  
  line%phi(:,id) = 0d0
  !line_profiles(:,kr,id) = 0d0
  
!   r = line%atom%vbroad(icell)*line%a(icell)/line%aeff(icell)
!   eta = 1.36603*r - 0.47719*r*r + 0.11116*r*r*r

  !normed
!   line%phi(:,id) = eta*(line%aeff(icell) / ((line%u(:)-omegav(nv))**2.+line%aeff(icell)**2.) / pi) &
!                        + (1-eta)*dexp(-(vvoigt(:)/line%aeff(icell))**2)/SQRTPI/line%aeff(icell)
                       
  line%phi(:,id) = line%r(icell)*line%aeff(icell) / ((line%u(:)-omegav(nv))**2.+line%aeff(icell)**2.) &
                       + line%r1(icell)*dexp(-(line%u(:)/line%aeff(icell))**2)


 RETURN
 END SUBROUTINE IProfile_thomson


 SUBROUTINE IProfile_cmf_to_obs(line,icell,x,y,z,x1,y1,z1,u,v,w,l, id, Nvspace, Omegav)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
  integer, intent(in) 							            :: icell, id, Nvspace
  real(kind=dp), intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  real(kind=dp), intent(in) 								:: Omegav(:)
  type (AtomicLine), intent(inout)								:: line
  integer													::  Nred, Nblue, i, j, nv
  real(kind=dp), dimension(line%Nlambda)                    :: u1, u1p, phi0
 

  i = line%i; j = line%j
  Nred = line%Nred; Nblue = line%Nblue

  !temporary here
  phi0 = 0d0
 write(*,*) " Interp not ready for profile"
 stop
 u1(:) = line%u(:)/line%atom%vbroad(icell)

 do nv=1, Nvspace 
 
         u1p(:) = u1(:) - omegav(nv)/line%atom%vbroad(icell)
         
         if (omegav(nv) == 0.0) then
         
          phi0 = phi0 + line%phi(:,icell)
         else

          phi0 = phi0 + &
               linear_1D_sorted(line%Nlambda,u1,line%phi(:,icell),line%Nlambda,u1p)
         endif
         	
 enddo
 
 phi0 = phi0 / Nvspace /sqrtpi / line%atom%vbroad(icell)


 RETURN
 END SUBROUTINE IProfile_cmf_to_obs
 
 !building
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
  line%phi(:,id) = 0d0
  
  !allocate(vv(line%Nlambda), vvoigt(line%Nlambda))


  Nzc = line%zm%Ncomponent
  if (.not.line%voigt) then !unpolarised line assumed even if line%polarizable
     norm  = line%atom%vbroad(icell) * Nvspace * SQRTPI

      do nv=1, Nvspace
      
         vvoigt(:) = (line%u - omegav(nv)) / line%atom%vbroad(icell)
         line%phi(:,id) = line%phi(:,id) + dexp(-(vvoigt(:))**2) !/ Nvspace

      !derivative of Gaussian:
!          F(Nblue:Nred) = F(Nblue:Nred) - &
!            2d0 * dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace * &
!            NLTEspec%lambda(Nblue:Nred) * CLIGHT / (line%atom%vbroad(icell) * line%lambda0)

      end do
      line%phi(:,id) = line%phi(:,id) / norm !/ (SQRTPI * vbroad)
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
       do nv=1, Nvspace !one iteration if 1) No velocity fields
      
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
          line%phi(:,id) = line%phi(:,id) + 5d-1 *(phi_zc(2,:) * dsin(gamma(nb))*dsin(gamma(nb)) + \
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
 
 SUBROUTINE write_profiles_ascii(unit, atom, delta_k)
  type(AtomType), intent(in) :: atom
  integer, intent(in) :: unit
  integer, intent(in), optional :: delta_k
  integer :: dk, kr, l, la, icell, Np
  type(AtomType), pointer :: HH
  
  write(*,*)  " Writing profiles for atom ", atom%ID
  HH => atmos%Atoms(1)%ptr_atom
  
  if (present(delta_k)) then
   dk = delta_k
   if (dk <= 0 .or. dk > atmos%Nspace) then
    dk = 1
    write(*,*) "delta_k cannot be out bound!"
   endif
  else
   dk = 1
  endif
  
  Np = n_cells
  Np = int((n_cells-1)/dk + 1)
  if (Np /= n_cells) then
   write(*,*) " Effective number of depth points written:", Np, n_cells
  endif
  
  open(unit, file=trim(atom%ID)//"_profiles.txt", status="unknown")
  write(unit,*) Np, atom%Nline
  do icell=1, n_cells, dk
   if (atmos%icompute_atomRT(icell) > 0) then
     write(unit,*) icell, atmos%T(icell), atmos%ne(icell), HH%n(1,icell), HH%n(HH%Nlevel,icell)
     write(unit,*) atom%vbroad(icell)*1e-3
     do kr=1, atom%Nline
      write(unit, *) kr, atom%lines(kr)%Nlambda, atom%lines(kr)%a(icell)
      write(unit, *) atom%lines(kr)%lambdamin, atom%lines(kr)%lambda0, atom%lines(kr)%lambdamax
      do la=1, atom%lines(kr)%Nlambda
       l = atom%lines(kr)%Nblue - 1 + la
       write(unit,*) NLTEspec%lambda(l), atom%lines(kr)%phi(la,icell)
      enddo
     enddo
   endif
  enddo

  close(unit)
 
 RETURN
 END SUBROUTINE write_profiles_ascii


END MODULE PROFILES
