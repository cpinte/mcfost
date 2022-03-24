module PROFILES

   use atmos_type, only             : B_project_from_vc, B_project_angles, VBROAD_atom
   use constant
   use atom_type
   use spectrum_type, only          : lambda, dk_min, dk_max
   use voigtfunctions, only         : Voigt, dirac_line, gate_line
   use math
   use getlambda, only              : hv
   use mcfost_env, only             : dp
   use molecular_emission, only     : v_proj
   use parametres
   use input
   use constantes, only             : tiny_dp, huge_dp

   implicit none

   procedure(local_profile_interp), pointer :: profile => null()
   integer, parameter :: NvspaceMax = 151
   integer, parameter :: NbspaceMax = 17

   contains



   function local_profile_v(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib,vmean)
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp), intent(inout)           :: vmean
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                                dv, omegav_mean
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: u1, u1p, local_profile_v


      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nred; Nblue = line%Nblue

      local_profile_v = 0d0
      u1(:) = ( (lambda - line%lambda0)/line%lambda0 ) * ( clight/line%atom%vbroad(icell) )

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
         omegav_mean = v0
      else

         Omegav = 0.0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/line%atom%vbroad(icell)*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      if (lsubstract_avg) then
         vmean = omegav_mean
      !    omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
      ! else
      !    omegav(1:Nvspace) = omegav(1:Nvspace) - vmean
      endif

      norm = Nvspace * line%atom%vbroad(icell) * sqrtpi

      if (line%voigt) then

         do nv=1, Nvspace

            u1p(:) = u1(:) - omegav(nv)/line%atom%vbroad(icell)

            local_profile_v(:) = local_profile_v(:) + Voigt(N, line%a(icell), u1p)

         enddo

      else
         do nv=1, Nvspace

            u1p(:) = u1(:) - omegav(nv)/line%atom%vbroad(icell)

            local_profile_v(:) = local_profile_v(:) + exp(-u1p**2)

         enddo
      endif


      local_profile_v(:) = local_profile_v(:) / norm

      return
   end function local_profile_v

  !gaussian are NOT interpolated (because it is not much faster but would be costly in memory to store all gaussian lines)
   function local_profile_interp(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib, vmean)
   ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                          x1,y1,z1, &      ! velocity field and magnetic field
                                          l_void_before,l_contrib !physical length of the cell
      integer 											:: Nvspace
      real(kind=dp), intent(inout)           :: vmean
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: norm
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                             dv, omegav_mean, vbroad
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: u1, u1p, local_profile_interp
      real(kind=dp), dimension(size(line%u)) :: u0

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nred; Nblue = line%Nblue
      vbroad = line%atom%vbroad(icell)
      !u0 = line%u / vbroad !avoid creating temporary array in interpolation
      !move below

      local_profile_interp = 0.0
      u1(:) = (lambda - line%lambda0)/line%lambda0 * clight/vbroad


      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0 / vbroad
         Nvspace = 1
         omegav_mean = v0 / vbroad
      else

         Omegav = 0.0_dp
         omegav(1) = v0 / vbroad
         v1 = v_proj(icell,x1,y1,z1,u,v,w)
         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)
         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w) / vbroad
         enddo
         omegav(Nvspace) = v1 / vbroad
         omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)
      endif
      if (lsubstract_avg) then
         vmean = omegav_mean * vbroad
      !    omegav(1:Nvspace) = omegav(1:Nvspace) - omegav_mean
      ! else
      !    omegav(1:Nvspace) = omegav(1:Nvspace) - vmean / vbroad
      endif


      if (line%voigt) then
         u0 = line%u / vbroad
         do nv=1, Nvspace

         u1p(:) = u1(:) - omegav(nv)

         local_profile_interp(:) = local_profile_interp(:) + linear_1D_sorted(size(line%u),u0,line%phi(:,icell),N,u1p)
         enddo

      else !Common gauss profile interpolated!
         u0 = line%atom%ug(:) / vbroad
         do nv=1, Nvspace
            u1p(:) = u1(:) - omegav(nv)
            local_profile_interp(:) = local_profile_interp(:) + &
               linear_1D_sorted(size(line%atom%ug),u0,line%atom%gauss_prof(:,icell),N,u1p)
         enddo
      endif

      local_profile_interp(:) = local_profile_interp(:) / Nvspace

      return
   end function local_profile_interp

   function local_profile_thomson(line,icell,lsubstract_avg, N, lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer                                :: Nvspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
      real(kind=dp)                          :: v0, v1, dv, delta_vol_phi, xphi, yphi, zphi, aeff, eta, ratio, aL, vbroad
      type (AtomicLine), intent(in)          :: line
      integer                                ::  Nred, Nblue, i, j, nv
      real(kind=dp), dimension(N)            :: u1, u1p, local_profile_thomson


      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nred; Nblue = line%Nblue
      vbroad = line%atom%vbroad(icell)


      local_profile_thomson = 0.0
      u1(:) = (lambda - line%lambda0)/line%lambda0 * clight

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
      else
         Omegav = 0d0
         omegav(1) = v0
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before+(real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
      endif


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



   function local_profile_dk(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
      !comoving (local) profile is shifted on the observed grid depending on the velocity.
      !The profile is defined on a N size grid which encompass the maximum possible displacement
      !due to the velocity, but the local profile is shifted only from i1:i2 (Nblue and Nred on this
      !lambda(N) grid)
      ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer                                :: Nvspace
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, dv
      type (AtomicLine), intent(in)          :: line
      integer                                :: i1, i2, i, j, nv, dk_mean, dk(NvspaceMax)
      real(kind=dp), dimension(N)            :: local_profile_dk

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      i1 = locate(lambda, line%lambdamin)
      i2 = locate(lambda, line%lambdamax)

      local_profile_dk = 0d0

      dk = 0
      v0 = v_proj(icell,x,y,z,u,v,w)

      if (lvoronoi) then
         dk(1) = int(1e-3 * v0/hv + 0.5)
         Nvspace = 1
      else
         dk(1) = int(1e-3 * v0/hv + 0.5)
         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/line%atom%vbroad(icell)*20.)),NvspaceMax)
         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            dk(nv) = int(1e-3*v_proj(icell,xphi,yphi,zphi,u,v,w)/hv + 0.5)
         enddo
         dk(Nvspace) = int(1e-3 * v1/hv + 0.5)
      endif
      dk_mean = 0
      !!if (lsubstract_avg) dk_mean = sum(dk(1:Nvspace))/Nvspace

      do nv=1, Nvspace
         local_profile_dk(i1+dk(nv):i2+dk(nv)) = local_profile_dk(i1+dk(nv):i2+dk(nv)) + line%phi(:,icell)
      enddo

      local_profile_dk(:) = local_profile_dk(:) / Nvspace

      return
   end function local_profile_dk


  !Currently, the magnetic field is assumed constant inside the cell and there is not projection
  !To Do split the magnetic field like the velocity field
   subroutine local_profile_zv(line,icell,lsubstract_avg,N,lambda, phi0, phiz, psiz, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
   ! phi = Voigt / sqrt(pi) / vbroad(icell)
      integer, intent(in)                    :: icell, N
      logical, intent(in)                    :: lsubstract_avg
      real(kind=dp), dimension(N), intent(in):: lambda
      real(kind=dp), intent(in)              :: x,y,z,u,v,w,& !positions and angles used to project
                                             x1,y1,z1, &      ! velocity field and magnetic field
                                             l_void_before,l_contrib !physical length of the cell
      integer                                :: Nvspace, Nzc, Nbspace
      real(kind=dp), dimension(NvspaceMax)   :: Omegav
!     real(kind=dp), dimension(NbspaceMax)   :: Omegab
      real(kind=dp)                          :: norm, vbroad, admp, cog, s2c, c2c, B, sigsq
      real(kind=dp)                          :: v0, v1, delta_vol_phi, xphi, yphi, zphi, &
                                             dv, omegav_mean, dummy
      type (AtomicLine), intent(in)          :: line
      integer                                :: Nred, Nblue, i, j, nv, nc
      real(kind=dp), dimension(N)            :: u1, u1p, ub
      real(kind=dp), intent(out)             :: phi0(N), phiZ(N,3), psiZ(N,3)
      real(kind=dp)                          :: H(N), F(N), psi(N,-1:1), phi(N,-1:1)
      logical                                :: lnot_magnetized
    

      B = B_project_angles(icell,x,y,z,u,v,w,cog,sigsq,c2c,s2c,lnot_magnetized)

      !Output arrays correspond to I, Q, U, V with phi0 for I and psiZ(:,i) for i=Q,U,V

      if (lnot_magnetized.or..not.(line%polarizable)) then
         !The test could be done elsewhere though
         phi0 = local_profile_v(line,icell,lsubstract_avg,N,lambda,x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib,dummy)
         return
      endif
    
      !tmp
      if (line%voigt) then
         admp = line%a(icell)
      else
         admp = 0.0_dp
         !to handle dispersion profile with gaussian (Voigt(a=0))
      endif

      Nvspace = NvspaceMax
      i = line%i; j = line%j
      Nred = line%Nred; Nblue = line%Nblue
      vbroad = line%atom%vbroad(icell)
      Nzc = line%zm%Ncomponent

      u1(:) = ( (lambda - line%lambda0)/line%lambda0 ) * ( clight/vbroad )

      v0 = v_proj(icell,x,y,z,u,v,w)
      if (lvoronoi) then
         omegav(1) = v0
         Nvspace = 1
      else
         Omegav = 0d0
         omegav(1) = v0

         v1 = v_proj(icell,x1,y1,z1,u,v,w)

         dv = abs(v1-v0)
         Nvspace = min(max(2,nint(dv/line%atom%vbroad(icell)*20.)),NvspaceMax)

         do nv=2, Nvspace-1
            delta_vol_phi = l_void_before + (real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
            xphi=x+delta_vol_phi*u
            yphi=y+delta_vol_phi*v
            zphi=z+delta_vol_phi*w
            omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w)
         enddo
         omegav(Nvspace) = v1
      endif
      omegav_mean = 0.0_dp
      !!if (lsubstract_avg) omegav_mean = sum(omegav(1:Nvspace))/real(Nvspace,kind=dp)

      norm = Nvspace * vbroad * sqrtpi
      psi = 0.0_dp; phi = 0.0_dp

      !psi and phi are the Zeeman components with -1, 0 and +1 depending on the deltaM

      !if (line%voigt) then

      do nv=1, Nvspace

         u1p(:) = u1(:) - (omegav(nv) - omegav_mean)/vbroad

         do nc=1,Nzc

            ub = u1p(:) - line%zm%shift(nc) * B * LARMOR * line%lambda0 * NM_TO_M / vbroad

            !not line%Nlambda. line%Nlambda==N only if not shift in index! (0 fields)
            H = Voigt(N, admp, ub, F)

            psi(:,line%zm%q(nc)) = psi(:,line%zm%q(nc)) + line%zm%strength(nc) * F(:) / norm

            phi(:,line%zm%q(nc)) = phi(:,line%zm%q(nc)) + line%zm%strength(nc) * H(:) / norm
             
         end do !components

      enddo
       
      !-> Can au Gaussian line be polarized ??
      !-> need a version for a = 0 that returns gaussian + equivalent dispersion profile with a = 0

!     else
!        do nv=1, Nvspace
! 
!           u1p(:) = u1(:) - (omegav(nv) - omegav_mean)/vbroad
! 
!           do nc=1,Nzc
!              ! the splitting is 0 if unpolarized 'cause zm%shift(nc=Nzc=1)=0d0
!              !there is a + omegaB because, -deltaLam^JJp_MMp=splitting = lamB * (gp*Mp - g*M)
!              ub = u1p(:) + line%zm%shift(nc) * B * LARMOR * line%lambda0 * NM_TO_M / vbroad
! 
!              !H = exp(-ub**2)
!              !F = -2 * ub(:) * H !?
!              ! H = Voigt(line%Nlambda, 0.0_dp, ub, F)
! 
!              psi(:,line%zm%q(nc)) = psi(:,line%zm%q(nc)) + line%zm%strength(nc) * F(:) / norm
! 
!              phi(:,line%zm%q(nc)) = phi(:,line%zm%q(nc)) + line%zm%strength(nc) * H(:) / norm
!           end do !components
! 
!        enddo
!     endif

      phi0(:) = 0.5 *(phi(:,0) * sigsq + 0.5 *(1.0+cog*cog) * (phi(:,-1)+phi(:,1)))

      !chiQ/chiI
      phiZ(:,1) = 0.5*(phi(:,0)-0.5*(phi(:,-1)+phi(:,1)))*c2c*sigsq
      !chiU/chiI
      phiZ(:,2) = 0.5*(phi(:,0)-0.5*(phi(:,-1)+phi(:,1)))*s2c*sigsq
      !chiV/chiI
      phiZ(:,3) = 0.5*(phi(:,-1)-phi(:,1))*cog

      !rhoQ/chiI
      psiz(:,1) = 0.5*(psi(:,0)-0.5*(psi(:,-1)+psi(:,1)))*c2c*sigsq
      !rhoU/chiI
      psiz(:,2) = 0.5*(psi(:,0)-0.5*(psi(:,-1)+psi(:,1)))*s2c*sigsq
      !rhoV/chiI
      psiz(:,3) = 0.5*(psi(:,-1)-psi(:,1))*cog

      return
   end subroutine local_profile_zv

end module profiles
!   SUBROUTINE write_profile(unit, icell, line, kc, wphi)
!     integer, intent(in) :: unit, icell, kc
!     type (AtomicLine), intent(in) :: line
!     real(kind=dp), intent(in) :: wphi
!     real(kind=dp) :: damp
!     write(unit, *) " icell = ", icell, " atom = ", line%atom%ID, " vbroad = ", line%atom%vbroad(icell)
!     write(unit, *) " l0 = ", line%lambda0, " lmin = ", line%lambdamin, " lmax = ", line%lambdamax
!     write(unit, *) " resol (nm) = ", lambda(line%Nblue+1)-lambda(line%Nblue), &
!          " resol(km/s) = ",1d-3 * clight*(lambda(line%Nblue+1)-lambda(line%Nblue))/lambda(line%Nblue)
!     if (allocated(line%a)) then
!        write(unit, *) " Area = ", wphi," damping = ", line%a(icell)
!     else
!        write(unit, *) " Area = ", wphi
!     endif
!     !!write(unit,*) " Vd (km/s) = ", 1e-3*line%atom%vbroad(icell),  " a = ", damp

!     RETURN
!   END SUBROUTINE write_profile

!   function local_profile_dirac(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
!    ! phi = Voigt / sqrt(pi) / vbroad(icell)
!    integer, intent(in) 							            :: icell, N
!    logical, intent(in)											:: lsubstract_avg
!    real(kind=dp), intent(in), dimension(N)						:: lambda
!    real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
!         x1,y1,z1, &      ! velocity field and magnetic field
!         l_void_before,l_contrib !physical length of the cell
!    integer 													:: Nvspace
!    real(kind=dp), dimension(NvspaceMax) 						:: Omegav
!    real(kind=dp) :: v0, v1, delta_vol_phi, xphi, yphi, zphi, t, vbroad, dv
!    type (AtomicLine), intent(in)								:: line
!    integer														::  Nred, Nblue, i, j, nv, la, j0, i0
!    real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_dirac


!    Nvspace = NvspaceMax
!    i = line%i; j = line%j
!    Nred = line%Nred; Nblue = line%Nblue
!    vbroad = line%atom%vbroad(icell)

!    local_profile_dirac = 0d0
!    u1(:) = (lambda - line%lambda0)/line%lambda0 * clight/vbroad

!    v0 = v_proj(icell,x,y,z,u,v,w)
!    if (lvoronoi) then
!       omegav(1) = v0 / vbroad
!       Nvspace = 1
!    else
!       Omegav = 0.0
!       omegav(1) = v0 / vbroad
!       v1 = v_proj(icell,x1,y1,z1,u,v,w)
!       dv = abs(v1-v0)
!       Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)

!       do nv=2, Nvspace-1
!          delta_vol_phi = l_void_before+(real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
!          xphi=x+delta_vol_phi*u
!          yphi=y+delta_vol_phi*v
!          zphi=z+delta_vol_phi*w
!          omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w) / vbroad
!       enddo
!       omegav(Nvspace) = v1 / vbroad
!    endif


!    do nv=1, Nvspace

!       u1p(:) = u1(:) - omegav(nv)


!       local_profile_dirac(:) = local_profile_dirac(:) + dirac_line(N, u1p)
!    enddo


!    local_profile_dirac(:) = local_profile_dirac(:) / Nvspace / vbroad / sqrtpi


!    return
!  end function local_profile_dirac

!  function local_profile_gate(line,icell,lsubstract_avg,N,lambda, x,y,z,x1,y1,z1,u,v,w,l_void_before,l_contrib)
!    ! phi = Voigt / sqrt(pi) / vbroad(icell)
!    integer, intent(in) 							            :: icell, N
!    logical, intent(in)											:: lsubstract_avg
!    real(kind=dp), intent(in), dimension(N)						:: lambda
!    real(kind=dp), intent(in) 					            	:: x,y,z,u,v,w,& !positions and angles used to project
!         x1,y1,z1, &      ! velocity field and magnetic field
!         l_void_before,l_contrib !physical length of the cell
!    integer 													:: Nvspace
!    real(kind=dp), dimension(NvspaceMax) 						:: Omegav
!    real(kind=dp) :: v0, v1, delta_vol_phi, xphi, yphi, zphi, t, vbroad, dv, max_u
!    type (AtomicLine), intent(in)								:: line
!    integer														::  Nred, Nblue, i, j, nv, la, j0, i0
!    real(kind=dp), dimension(N)			                    	:: u1, u1p, local_profile_gate


!    Nvspace = NvspaceMax
!    i = line%i; j = line%j
!    Nred = line%Nred; Nblue = line%Nblue
!    vbroad = line%atom%vbroad(icell)

!    local_profile_gate = 0d0
!    u1(:) = (lambda - line%lambda0)/line%lambda0 * clight/vbroad
!    max_u = (lambda(Nred) - line%lambda0) / line%lambda0 * clight / vbroad

!    v0 = v_proj(icell,x,y,z,u,v,w)
!    if (lvoronoi) then
!       omegav(1) = v0 / vbroad
!       Nvspace = 1
!    else
!       Omegav = 0d0
!       omegav(1) = v0 / vbroad
!       v1 = v_proj(icell,x1,y1,z1,u,v,w)
!       dv = abs(v1-v0)
!       Nvspace = min(max(2,nint(dv/vbroad*20.)),NvspaceMax)

!       do nv=2, Nvspace-1
!          delta_vol_phi = l_void_before+(real(nv,kind=dp))/(real(Nvspace,kind=dp)) * l_contrib
!          xphi=x+delta_vol_phi*u
!          yphi=y+delta_vol_phi*v
!          zphi=z+delta_vol_phi*w
!          omegav(nv) = v_proj(icell,xphi,yphi,zphi,u,v,w) / vbroad
!       enddo
!       omegav(Nvspace) = v1 / vbroad
!    endif


!    do nv=1, Nvspace

!       u1p(:) = u1(:) - omegav(nv)


!       local_profile_gate(:) = local_profile_gate(:) + gate_line(N, u1p, max_u)
!    enddo


!    local_profile_gate(:) = local_profile_gate(:) / Nvspace / vbroad / sqrtpi


!    return
!  end function local_profile_gate

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