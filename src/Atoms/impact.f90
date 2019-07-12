MODULE IMPACT
! ------------------------------------------------------------------------------------- ! 
 ! Hydrogen Collisional excitation / ionisation.
 !Collision coefficients C(j,i) = ne * CE * gi/gj * sqrt(T) in s^-1 with 
 ! CE in s^-1 K^-1/2 m^3.
 !
 ! - Effective collision strength 
 ! 1s - 2s; 2s - 2p
 !    T. T Scholz, H.R.J Walters, P.G. Burke, and M.P. Scott 1990,
 !    MNRAS 242, 692-697
 !
 ! - Electron impact excitation rates
 !   C. Giovanardi, A. Natta, and F. Palla 1987,
 !   Astron. Astrophys. Suppl. 70, 269-280
 !   C. Giovanardi, A, and F. Palla 1989
 !   Astron. Astrophys. Suppl. 77, 157-160 
 !   
 !   and also:
 !    E.S Chang, E.H. Avrett, and R. Loeser 1991, 
 !    A&A 247, 580-583
 !
 ! Other Elements
 !
 !
! ------------------------------------------------------------------------------------- !

 use atmos_type, only : atmos
 use mcfost_env, only : dp
 use constant
 
 IMPLICIT NONE
 
 !abscissa and weights for Gauss Laguerre integration of function int_0înf f(x)*exp(-x) dx
 real(kind=dp), dimension(8) :: x_laguerre, w_laguerre
 !integ = sum(y(x_laguerre)*w_laguerre) with y = f(x)*exp(-x)
 data x_laguerre / 0.170279632305d0,  0.903701776799d0,  2.251086629866d0,  4.266700170288d0, &
                   7.045905402393d0, 10.758516010181d0, 15.740678641278d0, 22.863131736889d0 /
                   
 data w_laguerre / 3.69188589342d-01, 4.18786780814d-01, 1.75794986637d-01, 3.33434922612d-02, & 
                   2.79453623523d-03, 9.07650877336d-05, 8.48574671627d-07, 1.04800117487d-09 /

 !Giovanardi coefficients, Cl are for low temperature , see table 1.
 !it is Cl_1s(1,:) = 1s_2s, C0, C1, C2, C3
 real(kind=dp), Cl_1s(5, 4), Cl_2s(3, 4), Cl_2p(3,4)
 data Cl_1s(1,:) / 2.302d-1, 5.248d-6, -1.144d-10, 8.24d-16 /
 
 
 
 CONTAINS
 
 CONTAINS
 
 FUNCTION Collision_Hydrogen(icell) result(Cji)
  integer :: icell, i, j
  double precision :: Cji(Hydrogen%Nlevel, Hydrogen%Nlevel)
  double precision :: nr_ji, CI(Hydrogen%Nlevel,Hydrogen%Nlevel), CE(Hydrogen%Nlevel,Hydrogen%Nlevel)
  
   Cji(:,:) = 0d0; CI = 0d0; CE(:,:) = 0d0 
   CALL Johnson_CI(icell, CI(:,Hydrogen%Nlevel)) !bound-free i->Nlevel
   CALL Johnson_CE(icell, CE) !among all levels

   do j=1,Hydrogen%Nlevel
    do i=1,Hydrogen%Nlevel
     nr_ji = Hydrogen%nstar(i,icell)/Hydrogen%nstar(j,icell)
     Cji(j,i) = CE(j,i) +  CI(i,j)  * nr_ji
     Cji(i,j) = CE(j,i)/nr_ji + CI(i,j)
    end do
   end do
   Cji(:,:) = Cji(:,:) * atmos%ne(icell)

 RETURN
 END FUNCTION Collision_Hydrogen

 SUBROUTINE Johnson_CI(icell, Cje)
 ! --------------------------------------------------- !
  ! Ionisation rate coefficient for
  ! Hydrogen atom, from
  ! L.C Johnson 
  ! ApJ 74:227-236, 1972 May 15; eq. 39
  !
  ! ne factorised
  !
  ! return C(i,j) with j = Nlevel (bound-free)
 ! --------------------------------------------------- !
   integer, intent(in) :: icell
   double precision, intent(out), dimension(:) :: Cje
   integer :: i, j, Nl
   double precision :: C0, pia0sq, rn, bn, n, An, En, yn, zn, S, Bnp
   type (AtomType) :: atom
   
   atom = atmos%Atoms(1)%ptr_atom

   C0 = dsqrt(8.*KBOLTZMANN*atmos%T(icell) / pi / M_ELECTRON)
   pia0sq = 2d0 * pi * RBOHR**2
   
   Nl = atom%Nlevel
   !Hydrogen level are ordered by n increasing, except for the continuum level
   !n = 1., but stops before Nlevel
   
   do i=1, Nl-1 !collision from neutral states to the ground state of H+
    n = real(i,kind=dp)
    if (i==1) then !n=i
     rn = 0.45
     bn = -0.603
    else 
     rn = 1.94*n**(-1.57)
     bn = 1d0/n * (4. - 18.63/n + 36.24/n/n - 28.09/n/n/n)
    end if

    ! in Joules
    En = E_RYDBERG / n / n !energie of level with different quantum number in 13.6eV: En = 13.6/n**2
    yn = En / KBOLTZMANN / atmos%T(icell)
    An = 32. / 3. / dsqrt(3d0) / pi * n  * (g0(n)/3. + g1(n)/4. + g2(n)/5.)
    Bnp = 2./3. * n*n * (5. + bn)
    zn = rn + yn

    S = C0 * pia0sq * (n*yn)**2 * (An*(E1(yn)/yn - E1(zn)/zn) + &
   			(Bnp - An*log(2*n*n))*(ksi_johnson(yn)-ksi_johnson(zn)))
   	!!write(*,*) i-1, Nl-1, "S=", S, S*dexp(yn)/dsqrt(atmos%T(icell)) 
    !check that otherwise multiply by dexp(yn)
    Cje(i) = S !RH -> dexp(yn) / dsqrt(atmos%T(icell)) !per ne
    		   !we compute it at icell so in fact we could avoid / sqrt(atmos%icell)
    		   !as col = Cje * sqrt(T) * exp(-de). Normally dE == En/kT=y so
    		   !it is not useful to multiply except to take into account the slightly diff
    		   !between En and (atom%E(Nl)-atom%E(i))
    !!write(*,*) En, (atom%E(Nl)-atom%E(i)) !Should be similar, because E(j)=13.6/j**2
   end do
 RETURN
 END SUBROUTINE Johnson_CI

 
 SUBROUTINE Johnson_CE(icell, Cje)
 ! ----------------------------------------------------- !
  ! Excitation rate coefficient for
  ! Hydrogen atom, from
  ! ApJ 74:227-236, 1972 May 15; eq. 36
  !
  ! CE = S = C(i,j) all transitions from j, i and i, j
  ! 
  ! ( -> transform C(i,j) to specifically C(j,i)
  ! S * dexp(y) * atom%g(i)/atom%g(j)
  ! = C(i,j) * exp(hnu/kt)*gi/gj = C(i,j) * ni/nj = C(j,i)
  ! at LTE: (gi/gj * nj/ni)  = exp(-hnu/kT) )
 ! ----------------------------------------------------- !
   integer, intent(in) :: icell
   double precision, intent(out), dimension(:,:) :: Cje
   integer :: i, j, Nl
   double precision :: C0, pia0sq, rn, bn, n, Ennp, y, z, S, Bnnp, En
   double precision :: np, x, fnnp, rnnp, Annp, Gaunt_bf
   type (AtomType) :: atom
   
   atom = atmos%Atoms(1)%ptr_atom

   C0 = dsqrt(8.*KBOLTZMANN*atmos%T(icell) / pi / M_ELECTRON)
   pia0sq = 2d0 * pi * RBOHR**2
   
   Nl = atom%Nlevel
   !Hydrogen level are ordered by n increasing, except for the continuum level
   !n = 1., but stops before Nlevel
   
   do i=1, Nl-1 !collision between neutral states, n to n'
    n = real(i,kind=dp)
    if (i==1) then !n=i
     rn = 0.45
     bn = -0.603
    else 
     rn = 1.94*n**(-1.57)
     bn = 1d0/n * (4. - 18.63/n + 36.24/n/n - 28.09/n/n/n)
    end if
    
    do j=i+1, Nl-1
     np = dble(j)!n'
     x = 1d0 - (n/np)**2 ! = Enn'/Rdybg
     !Gauntfactor * 32./3./dsqrt(3.)/pi * n/np**3 /x**3
     Gaunt_bf = g0(n) + g1(n)/x + g2(n)/x/x
     fnnp = Gaunt_bf * 32./3./dsqrt(3d0)/pi * n / np / np /np / x / x / x
     rnnp = rn * x
     Annp = 2d0 * n*n*fnnp/x
    ! in Joules
     En = E_RYDBERG / n / n !energie of level with different quantum number in 13.6eV = ionisation E of n
     y = x * En / KBOLTZMANN / atmos%T(icell) !x = ratio of E/En
     Bnnp = 4d0 * (n**4)/(np**3) / x / x * (1. + 4./3. /x + bn/x/x)
     z = rnnp + y
   
     S = C0 * pia0sq * n*n*y*y/x * (Annp*((1./y + 0.5)*E1(y)-(1./z + 0.5)*E1(z))+&
     	(Bnnp-Annp*dlog(2*n*n/x))*(E2(y)/y - E2(z)/z))
     
     Cje(j,i) = S * dexp(y) * atom%g(i)/atom%g(j)
     !!Cje(i,j) = S
     !write(*,*) atom%E(j) - atom%E(i), En*x !because x = deltaE/En
    end do !over np
   end do  !over n

 RETURN
 END SUBROUTINE Johnson_CE
 
 FUNCTION ksi_johnson(t) result(y)
  double precision :: t, y
  !E0
  y = dexp(-t)/t - 2d0*E1(t) + E2(t)
   
 RETURN
 END FUNCTION ksi_johnson
 
 FUNCTION g0 (n) result(g)
  double precision :: g, n
  
  SELECT CASE (int(n))
  CASE (1)
   g = 1.1330
  CASE (2)
   g = 1.0785
  CASE DEFAULT
   g = 0.9935 + 0.2328/n - 0.1296 / n / n
  END SELECT
 
 RETURN
 END FUNCTION g0
 
 FUNCTION g1 (n) result(g)
  double precision :: g, n
  
  SELECT CASE (int(n))
  CASE (1)
   g = -0.4059
  CASE (2)
   g = -0.2319
  CASE DEFAULT
   g = -1d0/n * (0.6282 - 0.5598/n + 0.5299 / n / n)
  END SELECT  
 
 RETURN
 END FUNCTION g1
 
 FUNCTION g2 (n) result(g)
  double precision :: g, n
 
  SELECT CASE (int(n))
  CASE (1)
   g = 0.07014
  CASE (2)
   g = 0.02947
  CASE DEFAULT
   g = 1d0/n/n * (0.3887 - 1.181 / n + 1.470 /n / n)
  END SELECT  
 
 RETURN
 END FUNCTION g2
 
 SUBROUTINE Scholz_et_al(icell, C)
  !define the output, and check results
  integer, intent(in) :: icell
  real(kind=dp), intent(inout) :: C !to define
  
  real(kind=dp) :: gamma_1s_2s, gamma_1s_2p, Omega_2s, Omega_2p
  real(kind=dp) :: b(8), c(6), x
  real(kind=dp), parameter: expdE = dexp(-7.5d-1), g0 = 2d0 !g value for ground-state of Hydrogen is 2
  
  if (atmos%T(icell) < 5d3 .or. atmos%T(icell) > 5d5) RETURN
  
  data b / 4.5168d-02,  2.8056d+01,  7.2945d+00, 2.4805d-01,  1.0044d-01, -1.1143d-02,  &
          -1.3432d-03,  3.7570d-04  / 
  
  data c / 3.6177d-01,  1.3891d+00,  5.0866d-01, -3.8011d-01,  1.0158d-01, -1.0072d-02 /
  
  x = (KBOLTZMANN * atmos%T(icell)) / E_RYDBERG
  dE_kT = E_RYDBERG
  !natural log
  gamma_1s_2s = b(1) * log(b(2)*x) * dexp(-b(3)*x) + b(4) + b(5)*x + b(6)*x*x * &
  				b(7)*x*x*x + b(8)*x*x*x*x
  				 !c0   !c1*x     !c2*x**3   !c3*x**3      !c4*x**4        !c5*x**5
  gamma_1s_2p = c(1) + c(2)*x + c(3)*x*x + c(4)*x*x*x + c(5)*x*x*x*x + c(6)*x*x*x*x*x
 
  Omega_2s = 8.63d-6 * gamma_1s_2s * expdE / gi / dsqrt(atmos%T(icell))
  Omega_2p = 8.63d-6 * gamma_1s_2p * expdE / gi / dsqrt(atmos%T(icell))

   
 RETURN
 END SUBROUTINE Scholz_et_al
 
 
 
 SUBROUTINE Giovanardi_et_al()
  

 RETURN
 END SUBROUTINE Giovanardi_et_al

END MODULE IMPACT