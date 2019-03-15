MODULE PROFILES

 use atmos_type, only                : atmos, B_project
 use constant
 use atom_type
 use spectrum_type, only			 : NLTEspec
 use voigtfunctions, only 			 : Voigt
 use broad, only 					 : Damping

 ! MCFOST's original
 use mcfost_env, only : dp
 use molecular_emission, only		 : v_proj
 use parametres
 use input
 use constantes, only				 : tiny_dp, huge_dp

 IMPLICIT NONE

 PROCEDURE(Iprofile), pointer :: Profile => null()
 PROCEDURE(Iprofile_lambda), pointer :: Profile_lambda => null()

 CONTAINS


 SUBROUTINE Iprofile (line, icell,x,y,z,x1,y1,z1,u,v,w,l, P, F)
 ! phi = Voigt / sqrt(pi) / vbroad(icell)
  integer, intent(in) 							            :: icell
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  type (AtomicLine), intent(in)								:: line
  double precision, dimension(NLTEspec%Nwaves)              :: vvoigt, phip, vv
  integer, parameter										:: NvspaceMax = 101
  double precision, dimension(NvspaceMax)					:: omegav
  integer													:: Nvspace, nv, Nred, Nblue, i, j
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: P
  double precision, intent(out), dimension(NLTEspec%Nwaves), optional :: F


  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   v0 = v_proj(icell,x,y,z,u,v,w)
   omegav(1) = v0
  end if

   
  if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0) 
     Nvspace = max(2,nint(20*dv/line%atom%vbroad(icell)))
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
  vv = 0d0
  vv(Nblue:Nred) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * line%atom%vbroad(icell))


  if (line%voigt) then
  !Now we have a pointer to atom in line. atom(n)%lines(kr)%atom => atom(n) 
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
          vvoigt(Nblue:Nred) = vv(Nblue:Nred) - &
                                       omegav(nv) / line%atom%vbroad(icell)

          P(Nblue:Nred) = P(Nblue:Nred) + &
            Voigt(line%Nlambda, line%adamp,vvoigt(Nblue:Nred), &
                  phip, VoigtMethod) / Nvspace

      end do
  else !Gaussian !only for checking
      do nv=1, Nvspace
      
         vvoigt(Nblue:Nred) = vv(Nblue:Nred) - omegav(nv) / line%atom%vbroad(icell)
         P(Nblue:Nred) = P(Nblue:Nred) + dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace 
         
      end do
 end if !line%voigt
 P(Nblue:Nred) = P(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))

 RETURN
 END SUBROUTINE IProfile
 
 SUBROUTINE ZProfile (line, icell,x,y,z,x1,y1,z1,u,v,w,l, P, F)
  integer, intent(in) 							            :: icell
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  type (AtomicLine), intent(in)								:: line
  double precision, dimension(NLTEspec%Nwaves)              :: vvoigt, phip, vv
  integer, parameter										:: NvspaceMax = 101, NbspaceMax=15
  double precision, dimension(NvspaceMax)					:: omegav
  double precision, dimension(NbspaceMax)					:: omegaB, gamma, chi
  integer													:: Nvspace, nv, Nred, Nblue, nc, &
  															   Nbspace, nb, Nzc, i, j
  double precision 											:: delta_vol_phi, xphi, yphi, zphi,&
  															   v0, v1, dv, b0, b1,g1,c1,dB
  double precision, intent(out), dimension(NLTEspec%Nwaves) :: P
  double precision, intent(out), dimension(NLTEspec%Nwaves), optional :: F


  NLTEspec%S_QUV = 0d0
  omegaB = 0d0
  ! v_proj in m/s at point icell
  omegav = 0d0
  Nvspace = 1
  if (.not.lstatic) then
   v0 = v_proj(icell,x,y,z,u,v,w)
   omegav(1) = v0
  end if


  b0 = B_project(icell,x,y,z,u,v,w,g1,c1)
  omegaB(1) = b0
  gamma(1) = g1; chi(1)=c1

   
  if (.not.lstatic .and. .not.lVoronoi) then ! velocity is varying across the cell
     v1 = v_proj(icell,x1,y1,z1,u,v,w)
     dv = dabs(v1-v0) 
     Nvspace = max(2,nint(20*dv/line%atom%vbroad(icell)))
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


  if (.not.lvoronoi) then
      b1 = B_project(icell,x1,y1,z1,u,v,w,g1,c1)
      Nbspace = NbspaceMax
!       dB = dabs(b1-b0) * LARMOR * (line%lambda0 * NM_TO_M) **2
!       Nbspace = max(2,nint(20*dB/atom%vbroad(icell)))
!       Nbspace = min(Nbspace,NbspaceMax)
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
  F = 0d0

  vv(Nblue:Nred) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * line%atom%vbroad(icell))
           
  Nzc = line%zm%Ncomponent
  if (.not.line%voigt) then !unpolarised line assumed
      do nv=1, Nvspace
      
         vvoigt(Nblue:Nred) = vv(Nblue:Nred) - omegav(nv) / line%atom%vbroad(icell)
         P(Nblue:Nred) = P(Nblue:Nred) + dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace
      !derivative of Gaussian:
!          F(Nblue:Nred) = F(Nblue:Nred) - &
!            2d0 * dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace * &
!            NLTEspec%lambda(Nblue:Nred) * CLIGHT / (line%atom%vbroad(icell) * line%lambda0)

      end do
      P(Nblue:Nred) = P(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))
!       F(Nblue:Nred) = F(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))
      RETURN
  end if
     
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)

  !Should work also for unpolarised voigt line because Ncz=1,S=0,q=0,shift=0
       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
        vvoigt(Nblue:Nred) = vv(Nblue:Nred) - &
                                       omegav(nv) / line%atom%vbroad(icell)
         do nb=1,Nbspace
         
          do nc=1,Nzc
             ! the splitting is 0 if unpolarized 'cause zm%shift(nc=Nzc=1)=0d0
             vvoigt(Nblue:Nred) = vvoigt(Nblue:Nred) - omegaB(nb) / line%atom%vbroad(icell) * &
                                  LARMOR * (line%lambda0 * NM_TO_M)**2 * line%zm%shift(nc)
             !The normalisation by Nzc is done when compute the strength*profiles.
             !In case of unpolarised line, Nzc = 1 and there is no need to normalised
             P(Nblue:Nred) = P(Nblue:Nred) + &
          					Voigt(line%Nlambda, line%adamp,vvoigt(Nblue:Nred), &
                  			phip, VoigtMethod) / Nvspace / Nbspace
             F(Nblue:Nred) = F(Nblue:Nred) + phip(Nblue:Nred)/Nvspace/Nbspace
             phip = 0d0
             !do something here
             CALL Error("Full profile not implemented")
          end do !components 
          
        end do !magnetic field     
        
       end do !velocity
       !check that if line is not polarised the Zeeman components are 0

  P(Nblue:Nred) = P(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))
  !F is 0 if unpolarised Voigt profile
  F(Nblue:Nred) = F(Nblue:Nred) / (SQRTPI * line%atom%vbroad(icell))

 RETURN
 END SUBROUTINE ZProfile

 
 SUBROUTINE WEAKFIELD_FLUX(ipix, jpix, ibin, iaz)
  use math, only : cent_deriv
  !Should be fine for non-overlapping lines
  integer, intent(in) :: ipix, jpix, ibin, iaz
  integer :: k, nat, Nred, Nblue
  type (AtomicLine) :: line
  type (AtomType), pointer :: atom
  double precision :: dlamB, Ipol(3,NLTEspec%Nwaves), dlamB21, dlamB22
  
  do nat=1,atmos%Natom
   atom => atmos%Atoms(nat)%ptr_atom
   do k=1,atom%Nline
     line = atom%lines(k)
     Nred = line%Nred; Nblue=line%Nblue
     if (.not.line%polarizable) CYCLE !Should be Bl at the "surface"
     dlamB = -line%g_lande_eff * maxval(atmos%Bxyz) * LARMOR * (line%lambda0**2) * NM_TO_M!can be 0 anyway
     Ipol = 0d0
     CALL cent_deriv(line%Nlambda,NLTEspec%lambda(Nblue:Nred),&
              NLTEspec%Flux(Nblue:Nred,ipix, jpix,ibin, iaz)*dlamB, Ipol(3,Nblue:Nred))
     NLTEspec%F_QUV(3,:,ipix,jpix,ibin,iaz) = NLTEspec%F_QUV(3,:,ipix,jpix,ibin,iaz) + &
                   Ipol(3,:)
   end do
   NULLIFY(atom)
  end do
 RETURN
 END SUBROUTINE

 !--> I should include the zeeman lines in phi, because it plays a role in opacity
 !and tau, so in the map calculations it might matter. But perhaps negligible
 SUBROUTINE Iprofile_lambda (line, icell,x,y,z,x1,y1,z1,u,v,w,l, P, F)
  integer, intent(in) 							            :: icell
  double precision, intent(in) 					            :: x,y,z,u,v,w,& !positions and angles used to project
                                				               x1,y1,z1, &      ! velocity field and magnetic field
                                				               l !physical length of the cell
  character(len=20)							                :: VoigtMethod = "HUMLICEK"
  integer 													:: i, j
  type (AtomicLine), intent(in)								:: line
  double precision, intent(out), dimension(1) :: P
  double precision, intent(out), dimension(1), optional :: F

  P(1) = 0d0

 RETURN
 END SUBROUTINE IProfile_lambda


END MODULE PROFILES
