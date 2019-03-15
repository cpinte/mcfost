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
  															   v0, v1, dv, dlamB, b0, b1,g1,c1,dB
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
  Nzc = 0
  if (line%polarizable) Nzc = line%zm%Ncomponent
  if (.not.line%voigt) then !futur use weak field only when line%Gauss
      do nv=1, Nvspace
         vvoigt(Nblue:Nred) = vv(Nblue:Nred) - omegav(nv) / line%atom%vbroad(icell)
         P(Nblue:Nred) = P(Nblue:Nred) + dexp(-(vvoigt(Nblue:Nred))**2) / Nvspace

      end do
      F = 0d0
      CALL Warning("Warning this line is not polarised because only Voigt profile for Zeeman calculation!")
      RETURN
  end if

  vv(Nblue:Nred) = (NLTEspec%lambda(Nblue:Nred)-line%lambda0) * &
           CLIGHT / (line%lambda0 * line%atom%vbroad(icell))
     
  !Computed before or change damping to use only line
  !CALL Damping(icell, line%atom, kr, line%adamp)

       ! init for this line of this atom accounting for Velocity fields
       do nv=1, Nvspace !one iteration if 1) No velocity fields or lstatic
                        !                 2) Voronoi grid is used                 
                        
        vvoigt(Nblue:Nred) = vv(Nblue:Nred) - &
                                       omegav(nv) / line%atom%vbroad(icell)
         do nb=1,Nbspace
           if (Nzc == 0) & !line not polarizable or weak field
                   P(Nblue:Nred) = P(Nblue:Nred) + &
          					Voigt(line%Nlambda, line%adamp,vvoigt(Nblue:Nred), &
                  			phip, VoigtMethod) / Nvspace / Nbspace
          do nc=1,Nzc
             vvoigt(Nblue:Nred) = vvoigt(Nblue:Nred) - omegaB(nb) / line%atom%vbroad(icell) * &
                                  LARMOR * (line%lambda0 * NM_TO_M)**2 * line%zm%shift(nc)
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
     
     if (line%polarizable) then
      if (line%ZeemanPattern==0) then
       !gamma(1) = !Bl = B*cos(gamma)
       !chi(1) = 0d0
       !derivative after, because we want dI/dlambda, I = exp(-tau)*(1.-exp(-dtau))*S
       dlamB = -line%zm%shift(1) * 1*LARMOR * (line%lambda0)**2 * NM_TO_M !result in nm
        NLTEspec%S_QUV(3,line%Nblue:line%Nred) = &
         NLTEspec%S_QUV(3,line%Nblue:line%Nred) + dlamB
      else
        CALL ERROR("Zeeman Opac not implemented yet")
      end if
     end if

 RETURN
 END SUBROUTINE ZProfile

 !--> I should include the zeeman lines in phi, because it plays a role in opacity
 !and tau
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
