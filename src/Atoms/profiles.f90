MODULE PROFILES

	use atmos_type, only				: B_project, VBROAD_atom
	use constant
	use atom_type
	use spectrum_type, only				: lambda, dk_min, dk_max
	use voigtfunctions, only			: Voigt, dirac_line, gate_line
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

	PROCEDURE(local_profile_v), pointer :: profile => null()
	integer, parameter :: NvspaceMax = 3, NbspaceMax = 50

	CONTAINS
	
	function gradv (icell, x,y,z,x1,y1,z1,u,v,w,l,dk)
		real(kind=dp) :: gradv, v0, v1
		real(kind=dp), intent(in) :: x,y,z,x1,y1,z1,u,v,w,l
		integer, intent(in) :: icell
		integer, intent(out) :: dk
		
		dk = 0
		v0 = v_proj(icell,x,y,z,u,v,w)
		v1 = v_proj(icell,x1,y1,z1,u,v,w)
		gradv = (v1-v0)/l

		dk = int(1e-3 * max(v1,v0) / hv + 0.5)
	
	return
	end function gradv

	!-> subroutine to call different line profile for different lines
	!more general and easier than having a line=> Voigt() or a line=>profile() ?

! 	subroutine line_profile(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l, prof)
! 		! phi = Voigt / sqrt(pi) / vbroad(icell)
! 		integer, intent(in) 							            :: icell, N
! 		logical, intent(in)											:: lsubstract_avg
! 		real(kind=dp), dimension(N), intent(in)						:: lambda
! 		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
!                                 				               			x1,y1,z1, &      ! velocity field and magnetic field
!                                 				               			l !physical length of the cell
! 		type (AtomicLine), intent(in)								:: line
! 		real(kind=dp), intent(out) 									:: prof
! 		
! 		
! 		if (line%voigt) then
! 			prof = local_profile_v(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
! 		elseif (line%pvoigt) then
! 			prof = local_profile_thomson(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
! 		elseif (line%voigt_interp .or. line%gauss_interp) then
! 			prof = local_profile_interp(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
! 		else!only gauss
! 			prof = local_profile_v(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
! 		endif
! 			
! 	return
! 	end subroutine line_profile 

 	!Might be better because for all lines of an atom or even for all lines of all atoms should be equivalent
	!if projection done before, we do not need x,y,z,l ect
 
 
	function local_profile_v(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
		! phi = Voigt / sqrt(pi) / vbroad(icell)
		integer, intent(in) 							            :: icell, N
		logical, intent(in)											:: lsubstract_avg
		real(kind=dp), dimension(N), intent(in)						:: lambda
		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
                                				               			x1,y1,z1, &      ! velocity field and magnetic field
                                				               			l !physical length of the cell
		integer 													:: Nvspace
		real(kind=dp), dimension(NvspaceMax) 						:: Omegav
		real(kind=dp) 												:: norm
		real(kind=dp) 												:: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
																		dv, omegav_mean
		type (AtomicLine), intent(in)								:: line
		integer														:: Nred, Nblue, i, j, nv
		real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_v
 

		Nvspace = NvspaceMax
		i = line%i; j = line%j
		Nred = line%Nred; Nblue = line%Nblue

		local_profile_v = 0d0
		u1(:) = ( (lambda - line%lambda0)/line%lambda0 ) * ( clight/line%atom%vbroad(icell) )
		
		Omegav = 0d0
		v0 = v_proj(icell,x,y,z,u,v,w)
		omegav(1) = v0
		v1 = v_proj(icell,x1,y1,z1,u,v,w)

		dv = abs(v1-v0)
		Nvspace = min(max(2,nint(dv/line%atom%vbroad(icell)*20.)),NvspaceMax)

		do nv=2, Nvspace-1
			delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
			xphi=x+delta_vol_phi*u
			yphi=y+delta_vol_phi*v
			zphi=z+delta_vol_phi*w
			omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
		enddo
		omegav(Nvspace) = v1

		omegav_mean = 0.0_dp
		!!if (lsubstract_avg) omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)

		norm = Nvspace * line%atom%vbroad(icell) * sqrtpi
		
        if (line%voigt) then

			do nv=1, Nvspace
 
				u1p(:) = u1(:) - (omegav(nv) - omegav_mean)/line%atom%vbroad(icell)
         
				local_profile_v(:) = local_profile_v(:) + Voigt(N, line%a(icell), u1p)
			enddo
			
		else
			do nv=1, Nvspace
			
				u1p(:) = u1(:) - (omegav(nv) - omegav_mean)/line%atom%vbroad(icell)

				local_profile_v(:) = local_profile_v(:) + exp(-u1p**2)
				
			enddo
		endif
			

		local_profile_v(:) = local_profile_v(:) / norm

	return
	end function local_profile_v
 
 
	function local_profile_thomson(line,icell,lsubstract_avg, N, lambda, x,y,z,x1,y1,z1,u,v,w,l)
		! phi = Voigt / sqrt(pi) / vbroad(icell)
		integer, intent(in) 							            :: icell, N
		logical, intent(in) 										:: lsubstract_avg
		real(kind=dp), dimension(N), intent(in)						:: lambda
		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
                                				               			x1,y1,z1, &      ! velocity field and magnetic field
                                				               			l !physical length of the cell
		integer 													:: Nvspace
		real(kind=dp), dimension(NvspaceMax) 						:: Omegav
		real(kind=dp) 												:: v0, v1, dv, delta_vol_phi, xphi, yphi, zphi, aeff, eta, ratio, aL, vbroad
		type (AtomicLine), intent(in)								:: line
		integer														::  Nred, Nblue, i, j, nv
		real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_thomson
 

		Nvspace = NvspaceMax
		i = line%i; j = line%j
		Nred = line%Nred; Nblue = line%Nblue
		vbroad = line%atom%vbroad(icell)


		local_profile_thomson = 0d0
		u1(:) = (lambda - line%lambda0)/line%lambda0 * clight
		
		Omegav = 0d0
		v0 = v_proj(icell,x,y,z,u,v,w)
		omegav(1) = v0
		v1 = v_proj(icell,x1,y1,z1,u,v,w)

		dv = abs(v1-v0)
		Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)
		
		do nv=2, Nvspace-1
      		delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
			xphi=x+delta_vol_phi*u
			yphi=y+delta_vol_phi*v
			zphi=z+delta_vol_phi*w
			omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
		enddo
		omegav(Nvspace) = v1

		
        if (line%voigt) then

          	!to do optimize:
          	!Can store them on the grid instead ! but it is fast to evaluate ?
		
			aL = line%a(icell) * vbroad !(m/s), adamp in doppler units
		
			aeff = (vbroad**5. + 2.69269*vbroad**4. * aL + 2.42843*vbroad**3. * aL**2. + &
					4.47163*vbroad**2. *aL**3. + 0.07842*vbroad*aL**4. + aL**5.)**(0.2)
          		
          	
			ratio = aL/aeff
			eta = 1.36603*ratio - 0.47719*(ratio*ratio) + 0.11116*(ratio*ratio*ratio)
			
			do nv=1, Nvspace
 
				u1p(:) = ( u1(:) - omegav(nv) ) 
         
				local_profile_thomson(:) = local_profile_thomson(:) + &
					eta * ( aeff/pi * (u1p(:)**2 + aeff**2)**(-1.0) ) + &
					(1.0_dp - eta) * exp(-(u1p(:)/aeff)**2) / aeff / sqrtpi
					
			enddo

			local_profile_thomson(:) = local_profile_thomson(:) / Nvspace
			
		else !pure Gauss, no approximation
			do nv=1, Nvspace
			
				u1p(:) = ( u1(:) - omegav(nv) ) / vbroad

				local_profile_thomson(:) = local_profile_thomson(:) + exp(-u1p**2)
				
			enddo
			local_profile_thomson(:) = local_profile_thomson(:) / Nvspace /sqrtpi / vbroad

		endif
 

	return
	end function local_profile_thomson

	!gaussian are also interpolated
	function local_profile_interp(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
		! phi = Voigt / sqrt(pi) / vbroad(icell)
		integer, intent(in) 							            :: icell, N
		logical, intent(in)											:: lsubstract_avg
		real(kind=dp), intent(in), dimension(N)						:: lambda
		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
                                				               			x1,y1,z1, &      ! velocity field and magnetic field
                                				               			l !physical length of the cell
		integer 													:: Nvspace
		real(kind=dp), dimension(NvspaceMax) 						:: Omegav
		real(kind=dp) :: v0, v1, delta_vol_phi, xphi, yphi, zphi, t, vbroad, dv
		type (AtomicLine), intent(in)								:: line
		integer														::  Nred, Nblue, i, j, nv, la, j0, i0
		real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_interp
! 		real(kind=dp), dimension(N)									:: locprof
 

		Nvspace = NvspaceMax
		i = line%i; j = line%j
		Nred = line%Nred; Nblue = line%Nblue
		vbroad = line%atom%vbroad(icell)

		local_profile_interp = 0d0
		u1(:) = (lambda - line%lambda0)/line%lambda0 * clight/vbroad
		
		Omegav = 0d0
		v0 = v_proj(icell,x,y,z,u,v,w)
		omegav(1) = v0 / vbroad
		v1 = v_proj(icell,x1,y1,z1,u,v,w)
		dv = abs(v1-v0)
		Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)
! 		write(*,*) "Nv = ", max(2,nint(dv/line%atom%vbroad(icell)*20.))
! 		write(*,*) v0/1e3, v1/1e3, dv/1e3, NvspaceMax
		do nv=2, Nvspace-1
      		delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
			xphi=x+delta_vol_phi*u
			yphi=y+delta_vol_phi*v
			zphi=z+delta_vol_phi*w
			omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w) / vbroad
		enddo
		omegav(Nvspace) = v1 / vbroad


!  					if (i==2 .and. j==3 .and. icell==n_cells .and. abs(v0) > 45d3 ) then
!  						write(*,*) v0
! 						open(unit=32, file="profile.txt",status="unknown")
! 						do la=1, size(line%u)
! 						
! 							write(32,*) line%u(la), line%phi(la,icell)
! 							
! 						enddo
! 						close(32)
! 
! 					endif

! 		if (line%voigt) then !faster to interpolate ?
!  						locprof = 0
			do nv=1, Nvspace
 
				u1p(:) = u1(:) - omegav(nv)


				local_profile_interp(:) = local_profile_interp(:) + linear_1D_sorted(size(line%u),line%u/vbroad,line%phi(:,icell),N,u1p)
! locprof = locprof + voigt(N, line%a(icell), u1p)/sqrtpi/vbroad/nvspace
			enddo

! 		else !gaussian, faster to evaluate in place
! 			do nv=1, Nvspace
!  
! 				u1p(:) = u1(:) - omegav(nv)
! 
!               	local_profile_interp(:) = local_profile_interp(:) + exp(-u1p*u1p) / sqrtpi / vbroad
!               	
! 			enddo		
! 		endif 
 
		local_profile_interp(:) = local_profile_interp(:) / Nvspace


!  					if (i==2 .and. j==3 .and. icell==n_cells .and. abs(v0) > 45d3) then
! 
! ! 						open(unit=32, file="profile_t.txt",status="unknown")
! ! 						do la=1, N
! ! 						
! ! 							write(32,*) vbroad*u1(la), locprof(la)
! ! 							
! ! 						enddo
! ! 						close(32)
! 
!  						write(*,*) line%a(icell), vbroad
! 						open(unit=32, file="profile_i.txt",status="unknown")
! 						do la=1, N
! 						
! 							write(32,*) vbroad*u1(la), local_profile_interp(la), u1p(la)
! 							
! 						enddo
! 						close(32)
! 						stop
! 					endif

	return
	end function local_profile_interp
	
	function local_profile_dirac(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
		! phi = Voigt / sqrt(pi) / vbroad(icell)
		integer, intent(in) 							            :: icell, N
		logical, intent(in)											:: lsubstract_avg
		real(kind=dp), intent(in), dimension(N)						:: lambda
		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
                                				               			x1,y1,z1, &      ! velocity field and magnetic field
                                				               			l !physical length of the cell
		integer 													:: Nvspace
		real(kind=dp), dimension(NvspaceMax) 						:: Omegav
		real(kind=dp) :: v0, v1, delta_vol_phi, xphi, yphi, zphi, t, vbroad, dv
		type (AtomicLine), intent(in)								:: line
		integer														::  Nred, Nblue, i, j, nv, la, j0, i0
		real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_dirac
 

		Nvspace = NvspaceMax
		i = line%i; j = line%j
		Nred = line%Nred; Nblue = line%Nblue
		vbroad = line%atom%vbroad(icell)

		local_profile_dirac = 0d0
		u1(:) = (lambda - line%lambda0)/line%lambda0 * clight/vbroad
		
		Omegav = 0d0
		v0 = v_proj(icell,x,y,z,u,v,w)
		omegav(1) = v0 / vbroad
		v1 = v_proj(icell,x1,y1,z1,u,v,w)
		dv = abs(v1-v0)
		Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)

		do nv=2, Nvspace-1
      		delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
			xphi=x+delta_vol_phi*u
			yphi=y+delta_vol_phi*v
			zphi=z+delta_vol_phi*w
			omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w) / vbroad
		enddo
		omegav(Nvspace) = v1 / vbroad


			do nv=1, Nvspace
 
				u1p(:) = u1(:) - omegav(nv)


				local_profile_dirac(:) = local_profile_dirac(:) + dirac_line(N, u1p)
			enddo

 
		local_profile_dirac(:) = local_profile_dirac(:) / Nvspace / vbroad / sqrtpi


	return
	end function local_profile_dirac
	
	function local_profile_gate(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
		! phi = Voigt / sqrt(pi) / vbroad(icell)
		integer, intent(in) 							            :: icell, N
		logical, intent(in)											:: lsubstract_avg
		real(kind=dp), intent(in), dimension(N)						:: lambda
		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
                                				               			x1,y1,z1, &      ! velocity field and magnetic field
                                				               			l !physical length of the cell
		integer 													:: Nvspace
		real(kind=dp), dimension(NvspaceMax) 						:: Omegav
		real(kind=dp) :: v0, v1, delta_vol_phi, xphi, yphi, zphi, t, vbroad, dv, max_u
		type (AtomicLine), intent(in)								:: line
		integer														::  Nred, Nblue, i, j, nv, la, j0, i0
		real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_gate
 

		Nvspace = NvspaceMax
		i = line%i; j = line%j
		Nred = line%Nred; Nblue = line%Nblue
		vbroad = line%atom%vbroad(icell)

		local_profile_gate = 0d0
		u1(:) = (lambda - line%lambda0)/line%lambda0 * clight/vbroad
		max_u = (lambda(Nred) - line%lambda0) / line%lambda0 * clight / vbroad
		
		Omegav = 0d0
		v0 = v_proj(icell,x,y,z,u,v,w)
		omegav(1) = v0 / vbroad
		v1 = v_proj(icell,x1,y1,z1,u,v,w)
		dv = abs(v1-v0)
		Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)

		do nv=2, Nvspace-1
      		delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
			xphi=x+delta_vol_phi*u
			yphi=y+delta_vol_phi*v
			zphi=z+delta_vol_phi*w
			omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w) / vbroad
		enddo
		omegav(Nvspace) = v1 / vbroad


			do nv=1, Nvspace
 
				u1p(:) = u1(:) - omegav(nv)


				local_profile_gate(:) = local_profile_gate(:) + gate_line(N, u1p, max_u)
			enddo

 
		local_profile_gate(:) = local_profile_gate(:) / Nvspace / vbroad / sqrtpi


	return
	end function local_profile_gate

	function local_profile_dk(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l)
	!comoving (local) profile is shifted on the observed grid depending on the velocity.
	!The profile is defined on a N size grid which encompass the maximum possible displacement
	!due to the velocity, but the local profile is shifted only from i1:i2 (Nblue and Nred on this
	!lambda(N) grid)
		! phi = Voigt / sqrt(pi) / vbroad(icell)
		integer, intent(in) 							            :: icell, N
		logical, intent(in)											:: lsubstract_avg
		real(kind=dp), dimension(N), intent(in)						:: lambda
		real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
                                				               			x1,y1,z1, &      ! velocity field and magnetic field
                                				               			l !physical length of the cell
		integer 													:: Nvspace
		real(kind=dp) 												:: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
																		dv
		type (AtomicLine), intent(in)								:: line
		integer														:: i1, i2, i, j, nv, dk_mean, dk(NvspaceMax)
		real(kind=dp), dimension(N)			                    	:: local_profile_dk
 
		Nvspace = NvspaceMax
		i = line%i; j = line%j
		i1 = locate(lambda, line%lambdamin)
		i2 = locate(lambda, line%lambdamax)

		local_profile_dk = 0d0
		
		dk = 0
		v0 = v_proj(icell,x,y,z,u,v,w)
		dk(1) = int(1e-3 * v0/hv + 0.5)
		v1 = v_proj(icell,x1,y1,z1,u,v,w)

		dv = abs(v1-v0)
		Nvspace = min(max(2,nint(dv/line%atom%vbroad(icell)*20.)),NvspaceMax)
		do nv=2, Nvspace-1
      		delta_vol_phi = (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l
			xphi=x+delta_vol_phi*u
			yphi=y+delta_vol_phi*v
			zphi=z+delta_vol_phi*w
			dk(nv) = int(1e-3*v_proj(icell,xphi,yphi,zphi,u,v,w)/hv + 0.5)
		enddo
		dk(Nvspace) = int(1e-3 * v1/hv + 0.5)
		dk_mean = 0
		!!if (lsubstract_avg) dk_mean = sum(dk(1:Nvspace))/Nvspace

		do nv=1, Nvspace
			local_profile_dk(i1+dk(nv):i2+dk(nv)) = local_profile_dk(i1+dk(nv):i2+dk(nv)) + line%phi(:,icell)
		enddo

		local_profile_dk(:) = local_profile_dk(:) / Nvspace

	return
	end function local_profile_dk
	
	
! 	function local_zprofile_vd !voigt + dispersion profiles
! 	
! 	return
! 	end function local_zprofile_vd
! 	
! 	function local_zprofile_thomson
! 	return
! 	end function local_zprofile_thomson
! 	
! 	function local_zprofile_interp
! 	return
! 	end function local_zprofile_interp
 
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
         line%phi(:,id) = line%phi(:,id) + exp(-(vvoigt(:))**2) !/ Nvspace

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
 
 SUBROUTINE write_profile(unit, icell, line, kc, wphi)
 	integer, intent(in) :: unit, icell, kc
 	type (AtomicLine), intent(in) :: line
 	real(kind=dp), intent(in) :: wphi
 	real(kind=dp) :: damp
 	!!CALL Damping(icell, line%atom, kc, damp)
 	!damp not computed, to much time ?
 	write(unit, *) " icell = ", icell, " atom = ", line%atom%ID, " vbroad = ", line%atom%vbroad(icell)
 	write(unit, *) " l0 = ", line%lambda0, " lmin = ", line%lambdamin, " lmax = ", line%lambdamax
 	write(unit, *) " resol (nm) = ", lambda(line%Nblue+1)-lambda(line%Nblue), " resol(km/s) = ",1d-3 * clight*(lambda(line%Nblue+1)-lambda(line%Nblue))/lambda(line%Nblue)
 	if (allocated(line%a)) then
 		write(unit, *) " Area = ", wphi," damping = ", line%a(icell)
 	else
  		write(unit, *) " Area = ", wphi
	endif
 	!!write(unit,*) " Vd (km/s) = ", 1e-3*line%atom%vbroad(icell),  " a = ", damp
 
 RETURN
 END SUBROUTINE 
 
!  SUBROUTINE write_profiles_ascii(unit, atom, delta_k)
!	use atmos, only : ....
!   type(AtomType), intent(in) :: atom
!   integer, intent(in) :: unit
!   integer, intent(in), optional :: delta_k
!   integer :: dk, kr, l, la, icell, Np
!   type(AtomType), pointer :: HH
!   
!   write(*,*)  " Writing profiles for atom ", atom%ID
!   HH => atmos%Atoms(1)%ptr_atom
!   
!   if (present(delta_k)) then
!    dk = delta_k
!    if (dk <= 0 .or. dk > atmos%Nspace) then
!     dk = 1
!     write(*,*) "delta_k cannot be out bound!"
!    endif
!   else
!    dk = 1
!   endif
!   
!   Np = n_cells
!   Np = int((n_cells-1)/dk + 1)
!   if (Np /= n_cells) then
!    write(*,*) " Effective number of depth points written:", Np, n_cells
!   endif
!   
!   open(unit, file=trim(atom%ID)//"_profiles.txt", status="unknown")
!   write(unit,*) Np, atom%Nline
!   do icell=1, n_cells, dk
!    if (atmos%icompute_atomRT(icell) > 0) then
!      write(unit,*) icell, atmos%T(icell), atmos%ne(icell), HH%n(1,icell), HH%n(HH%Nlevel,icell)
!      write(unit,*) atom%vbroad(icell)*1e-3
!      do kr=1, atom%Nline
!       write(unit, *) kr, atom%lines(kr)%Nlambda, atom%lines(kr)%a(icell)
!       write(unit, *) atom%lines(kr)%lambdamin, atom%lines(kr)%lambda0, atom%lines(kr)%lambdamax
!       do la=1, atom%lines(kr)%Nlambda
!        l = atom%lines(kr)%Nblue - 1 + la
!        write(unit,*) NLTEspec%lambda(l), atom%lines(kr)%phi(la,icell)
!       enddo
!      enddo
!    endif
!   enddo
! 
!   close(unit)
!  
!  RETURN
!  END SUBROUTINE write_profiles_ascii



END MODULE PROFILES
