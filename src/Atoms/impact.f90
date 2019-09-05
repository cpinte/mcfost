MODULE IMPACT
! ----------------------------------------------------------------------
!Note: C in Unit of number density is m^-3 and atom%Ckij and atom%C in s^-1.
!
!Convention: C_ij = C[i][j] represents the
!transition j --> i = Cul
! C_ij = C[ith ligne][jth column]
! fortran is column row
!     ij = (i-1)*atom%Nlevel +j : j-->i
!     ji = (j-1)*atom%Nlevel +i : i-->j
! ----------------------------------------------------------------------

! ------------------------------------------------------------------------------------- ! 
 ! Hydrogen Collisional excitation / ionisation.
 ! Collision coefficients C(j,i) = ne * CE * gi/gj * sqrt(T) in s^-1 with 
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
 ! TODO: Integration of the cross-section itself. 
! ------------------------------------------------------------------------------------- !

 use atmos_type, only : atmos, Hydrogen
 use atom_type!, only : AtomicLine, AtomicContinuum, AtomType
 use mcfost_env, only : dp
 use constant
 !!use special_functions
 use math
 use utils, only : interp_sp
 
 IMPLICIT NONE
 
 !abscissa and weights for Gauss Laguerre integration of function int_0^inf f(x)*exp(-x) dx
 real(kind=dp), dimension(8) :: x_laguerre, w_laguerre
 !integ = sum(y(x_laguerre)*w_laguerre) with y = f(x)*exp(-x)
 data x_laguerre / 0.170279632305d0,  0.903701776799d0,  2.251086629866d0,  4.266700170288d0, &
                   7.045905402393d0, 10.758516010181d0, 15.740678641278d0, 22.863131736889d0 /
                   
 data w_laguerre / 3.69188589342d-01, 4.18786780814d-01, 1.75794986637d-01, 3.33434922612d-02, & 
                   2.79453623523d-03, 9.07650877336d-05, 8.48574671627d-07, 1.04800117487d-09 /

 !Impact Parameter approximation of Seaton 1962 vol 72
 !I interpolate using its calculations, could be better to compute everything direcly
 !but beta involves calculations of transcendental equation using modified Bessel function K1 and K0..
 real, dimension(24) :: beta, zeta, phi
 data beta  / 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, &
              0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 / 
 !these are two functions of beta
 data zeta  / 1.0, 1.015, 1.030, 1.035, 1.035, 1.025, 1.010, 0.963, 0.9, 0.829, 0.754, &
 			  0.68, 0.608, 0.540, 0.476, 0.418, 0.365, 0.318, 0.276, 0.239, 0.206, 0.177, 0.152, 0.130 /
              
 data phi  / 3.1, 2.392, 1.972, 1.674, 1.444, 1.259, 0.974, 0.766, 0.608, 0.486, 0.390, &
 			 0.314, 0.253, 0.2049, 0.1661, 0.1347, 0.1095, 0.0889, 0.0723, 0.0589, 0.0480, &
 			 0.0390, 0.0319 /
              
 
 
 CONTAINS
 
 FUNCTION Seaton_IP(Sji)
 !
 ! Seaton Impact Parameter approximation
 ! eq. 17 + eq. 19 and eq. 
 !
  real(kind=dp) :: Seaton_IP, Sji
  
  Seaton_IP = 4d0 * pia0squarex2 * E_RYDBERG**2  * Sji
  
  RETURN
 END FUNCTION Seaton_IP
 
 !Check addition of Collision ! probably a missing factor<

 
 SUBROUTINE Collision_with_electrons(icell, atom, C)
 !Landolt-Boernstein, Volume 2, Astronomy and Astrophysics,
 !       subvolume b, Stars and Star clusters, Springer-Verlag, 98-100).
 ! For optically forbidden transition, set correct = 0.1, not really optimal
  integer, intent(in) :: icell
  type(AtomType), intent(in) :: atom
  real(kind=dp), intent(out) :: C(atom%Nlevel, atom%Nlevel)
  integer :: i, j, kr, ic
  type(AtomicLine) :: line
  type(AtomicContinuum) :: cont
  real :: gbar_i, correct = 1.!, xions, xneutrals
  real(kind=dp) :: C2_ions, C2_neutrals, dE, C_bf, deltam, R0, beta0, ri, rj, neff, Sji, Wi
  logical :: forbidden
  
  
   deltam = 1. + M_ELECTRON / (atom%weight * AMU)
   !xions = 0.
   !xneutrals = 0.68
   C2_ions = 3.96d-6 !s^-1 m^3 K^3/2
   C2_neutrals = 2.15d-6
   C_bf = 1.55d11 ! s^-1 K^1/2 m^1
   
   C(:,:) = 0d0
   
   do i=1,atom%Nlevel !between pair of levels resulting in a b-b transition
    Wi = atmos%T(icell) * KBOLTZMANN
    do j=i+1,atom%Nlevel
     
      dE = (atom%E(j) - atom%E(i)) / Wi !Joules I hope
     
      !eq. 34
      if (atom%stage(j) == atom%stage(i)) then !collisions between bound-bound transitions
       
       forbidden = .false.
       correct = 1.
       lines_loop : do kr=1, atom%Nline
         line = atom%lines(kr)
          if (line%i == i .and. line%j == j) then
           forbidden = .true.
           correct = 0.1
           exit lines_loop
          endif
       enddo lines_loop

       if (atom%stage(i) == 0) then !neutral
       
         if (forbidden) then !Modified Van Regemorter with correct.
 
            C(i,j) = correct * C2_neutrals * atmos%T(icell)**(-1.5) * dexp(-dE) * line%fosc * &
                         dE**(-1.68) !power of 1+xneutrals
         else !optically allowed: could also use Van Regemorter with correct = 1.
            ic = locate(atom%stage*1d0, 1d0*atom%stage(i)+1d0)
            if (ic == i) then
             write(*,*) "Error, couldn't find continuum for level", atom%label(i)
            end if
            neff = (atom%stage(i)+1) * dsqrt(E_RYDBERG/deltam / (atom%E(ic)-atom%E(i)))
            write(*,*) "(n,l,Z)", neff, atom%Lorbit(i), atom%stage(i)+1
            ri = atomic_orbital_radius(neff, atom%Lorbit(i),atom%stage(i)+1)
            ic = locate(atom%stage*1d0, 1d0*atom%stage(j)+1d0)
            if (ic == j) then
             write(*,*) "Error, couldn't find continuum for level", atom%label(j)
            end if
            neff = (atom%stage(j)+1) * dsqrt(E_RYDBERG/deltam / (atom%E(ic)-atom%E(j)))
            write(*,*) "(n,l,Z)", neff, atom%Lorbit(j), atom%stage(j)+1
            ri = atomic_orbital_radius(neff, atom%Lorbit(j),atom%stage(j)+1)
            R0 = min(ri, rj)
            beta0 = R0 * dsqrt(Wi*deltam/E_RYDBERG) * dE
            write(*,*) "R0 (a0) = ", R0," beta0 = ", real(beta0), beta0
            Sji = line%fosc*interp_sp(phi,beta, real(beta0)) / deltam / deltam / dE / Wi**2
            C(i,j) = Seaton_IP(Sji)
         
         end if !forbidden neutral stage 
       else !ions !Van Regemorter
          C(i,j) = correct * C2_ions * atmos%T(icell)**(-1.5) * dexp(-dE) * line%fosc * dE
       endif

      endif !eq. 34
    
    enddo
   enddo
   
   !now bound-free eq. 26
   !or eq. 9.60 Hubeny Mihalas
   do kr=1,atom%Ncont !only over continua
    cont = atom%continua(kr)
    i = cont%i; j = cont%j
    dE = (atom%E(j) - atom%E(i)) / (KBOLTZMANN * atmos%T(icell))
    
    select case (atom%stage(i))
     case (0)
      gbar_i = 0.1 ! stage neutral
     case (1) 
      gbar_i = 0.2 ! Z = 2
     case default
      gbar_i = 0.3 ! Z > 2
    end select
    
    C(i,j) = C_bf * atmos%T(icell)**(-0.5) * cont%alpha0 * gbar_i * (dE**-1) * &
              dexp(-dE)
    
   enddo
 
 RETURN
 END SUBROUTINE Collision_with_electrons
 
 FUNCTION Psi_Drawin(u)
  real(kind=dp) :: Psi_Drawin, u

  Psi_Drawin = dexp(-u) / (1d0 + 2d0/u)
  
 RETURN
 END FUNCTION Psi_Drawin

 FUNCTION Psi_Drawin_ionisation(u)
  real(kind=dp) :: Psi_Drawin_ionisation, u

  Psi_Drawin_ionisation = (1d0 + 2d0/u) * dexp(-u)
  Psi_drawin_ionisation = Psi_drawin_ionisation / &
     1d0 + 2d0 * M_ELECTRON / (u * (M_ELECTRON + Hydrogen%weight*AMU))
  
 RETURN
 END FUNCTION Psi_Drawin_ionisation
 
 SUBROUTINE Collision_with_HI(icell, atom, C)
  !Hubeny Mihalas, eq. 9.62 Drawin approximation
  integer :: icell
  type (AtomType) :: atom
  real(kind=dp), intent(out) :: C(atom%Nlevel, atom%Nlevel)
  type (AtomicContinuum) :: cont
  type (AtomicLine) :: line
  real(kind=dp) :: psi0, dE, mu, u0, C0, qij, f, psi, deltam
  integer :: kr, i, j
  
  !reduce mass of the system
  !mu = atomic_weights(atom%periodic_table)*atomic_weights(1) /&
  !    (atomic_weights(atom%periodic_table) + atomic_weights(1))
  mu = atom%weight * Hydrogen%weight / (atom%weight + Hydrogen%weight)
  deltam = 1. + M_ELECTRON/ (atom%weight * AMU) !correction to E_RYDBERG = Rinf
  C0 = 8* pia0squarex2 * (2d0 * KBOLTZMANN * atmos%T(icell) / PI / mu)**(0.5)
        !16 * pi * RBHOR**2
  do kr=1, atom%Ntr
   
   select case (atom%at(kr)%trtype)
   
    case ("ATOMIC_LINE") !collisional excitation
     line = atom%lines(atom%at(kr)%ik)
     i = line%i; j = line%j
     dE = atom%E(j) - atom%E(i)
     u0 = dE / KBOLTZMANN / atmos%T(icell)
     f = line%fosc
     psi = Psi_drawin(u0)
     
    case ("ATOMIC_CONTNUUM") !collisional ionisation
     cont = atom%continua(atom%at(kr)%ik)
     i = cont%i; j = cont%j
     dE = atom%E(j) - atom%E(i)
     f = real(i)/2. ! HERE IT IS NOT CORRECT
     u0 = dE / KBOLTZMANN / atmos%T(icell)
     psi = Psi_Drawin_ionisation(u0)
    
    case default
     write(*,*) atom%at(kr)%trtype, " unkown!"
     stop
   end select
   
   !Drawin approx
   qij = (E_RYDBERG/deltaM/dE)**2 * f * atom%weight/Hydrogen%weight * &
     M_ELECTRON / (M_ELECTRON + AMU * Hydrogen%weight) * psi
     
   C(i,j) = qij
  end do
  
 RETURN
 END SUBROUTINE Collision_with_HI
 
 SUBROUTINE Charge_transfer_with_H()
 
 RETURN
 END SUBROUTINE Charge_transfer_with_H
 
 FUNCTION Collision_Hydrogen(icell) result(Cji)
  integer :: icell, i, j
  real(kind=dp) :: Cji(Hydrogen%Nlevel, Hydrogen%Nlevel)
  real(kind=dp) :: nr_ji, CI(Hydrogen%Nlevel,Hydrogen%Nlevel), CE(Hydrogen%Nlevel,Hydrogen%Nlevel)
  
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
   real(kind=dp), intent(out), dimension(:) :: Cje
   integer :: i, j, Nl
   real(kind=dp) :: C0, rn, bn, n, An, En, yn, zn, S, Bnp, deltaM
   type (AtomType) :: atom
   
   atom = atmos%Atoms(1)%ptr_atom
   deltam = 1. + M_ELECTRON/ (atom%weight * AMU)

   C0 = dsqrt(8.*KBOLTZMANN*atmos%T(icell) / pi / M_ELECTRON)
   
   Nl = atom%Nlevel
   !Hydrogen level are ordered by n increasing, except for the continuum level
   !n = 1., but stops before Nlevel
   
   do i=1, Nl-2 !collision from neutral states to the ground state of H+
    n = real(i,kind=dp)
    if (i==1) then !n=i
     rn = 0.45
     bn = -0.603
    else 
     rn = 1.94*n**(-1.57)
     bn = 1d0/n * (4. - 18.63/n + 36.24/n/n - 28.09/n/n/n)
    end if

    ! in Joules
    En = E_RYDBERG /deltam / n / n !energie of level with different quantum number in 13.6eV: En = 13.6/n**2
    yn = En / KBOLTZMANN / atmos%T(icell)
    An = 32. / 3. / dsqrt(3d0) / pi * n  * (g0(n)/3. + g1(n)/4. + g2(n)/5.)
    Bnp = 2./3. * n*n * (5. + bn)
    zn = rn + yn

    S = C0 * pia0squarex2 * (n*yn)**2 * (An*(E1(yn)/yn - E1(zn)/zn) + &
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
   real(kind=dp), intent(out), dimension(:,:) :: Cje
   integer :: i, j, Nl
   real(kind=dp) :: C0, rn, bn, n, Ennp, y, z, S, Bnnp, En
   real(kind=dp) :: np, x, fnnp, rnnp, Annp, Gaunt_bf, deltam
   type (AtomType) :: atom
   
   atom = atmos%Atoms(1)%ptr_atom
   deltam = 1. + M_ELECTRON/ (atom%weight * AMU)

   C0 = dsqrt(8.*KBOLTZMANN*atmos%T(icell) / pi / M_ELECTRON)
   
   Nl = atom%Nlevel
   !Hydrogen level are ordered by n increasing, except for the continuum level
   !n = 1., but stops before Nlevel
   
   do i=1, Nl-2 !collision between neutral states, n to n'
    n = real(i,kind=dp)
    if (i==1) then !n=i
     rn = 0.45
     bn = -0.603
    else 
     rn = 1.94*n**(-1.57)
     bn = 1d0/n * (4. - 18.63/n + 36.24/n/n - 28.09/n/n/n)
    end if
    
    do j=i+1, Nl-2
     np = dble(j)!n'
     x = 1d0 - (n/np)**2 ! = Enn'/Rdybg
     !Gauntfactor * 32./3./dsqrt(3.)/pi * n/np**3 /x**3
     Gaunt_bf = g0(n) + g1(n)/x + g2(n)/x/x
     fnnp = Gaunt_bf * 32./3./dsqrt(3d0)/pi * n / np / np /np / x / x / x
     rnnp = rn * x
     Annp = 2d0 * n*n*fnnp/x
    ! in Joules
     En = E_RYDBERG /deltam / n / n !energie of level with different quantum number in 13.6eV = ionisation E of n
     y = x * En / KBOLTZMANN / atmos%T(icell) !x = ratio of E/En
     Bnnp = 4d0 * (n**4)/(np**3) / x / x * (1. + 4./3. /x + bn/x/x)
     z = rnnp + y
   
     S = C0 * pia0squarex2 * n*n*y*y/x * (Annp*((1./y + 0.5)*E1(y)-(1./z + 0.5)*E1(z))+&
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
  real(kind=dp) :: g, n
  
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
  real(kind=dp) :: g, n
  
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
  real(kind=dp) :: g, n
 
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
  real(kind=dp) :: b(8), cc(6), x, dE_kT
  real(kind=dp), parameter :: expdE = dexp(-7.5d-1), gi = 2d0 !g value for ground-state of Hydrogen is 2
  
  if (atmos%T(icell) < 5d3 .or. atmos%T(icell) > 5d5) RETURN
  
  data b / 4.5168d-02,  2.8056d+01,  7.2945d+00, 2.4805d-01,  1.0044d-01, -1.1143d-02,  &
          -1.3432d-03,  3.7570d-04  / 
  
  data cc / 3.6177d-01,  1.3891d+00,  5.0866d-01, -3.8011d-01,  1.0158d-01, -1.0072d-02 /
  
  x = (KBOLTZMANN * atmos%T(icell)) / E_RYDBERG / (1d0 + M_ELECTRON/hydrogen%weight*AMU)
  dE_kT = E_RYDBERG
  !natural log
  gamma_1s_2s = b(1) * log(b(2)*x) * dexp(-b(3)*x) + b(4) + b(5)*x + b(6)*x*x * &
  				b(7)*x*x*x + b(8)*x*x*x*x
  				 !c0   !c1*x     !c2*x**3   !c3*x**3      !c4*x**4        !c5*x**5
  gamma_1s_2p = cc(1) + cc(2)*x + cc(3)*x*x + cc(4)*x*x*x + cc(5)*x*x*x*x + cc(6)*x*x*x*x*x
 
  Omega_2s = 8.63d-6 * gamma_1s_2s * expdE / gi / dsqrt(atmos%T(icell))
  Omega_2p = 8.63d-6 * gamma_1s_2p * expdE / gi / dsqrt(atmos%T(icell))

   
 RETURN
 END SUBROUTINE Scholz_et_al
 
  !Giovanardi coefficients, Cl are for low temperature , see table 1.
 !it is Cl_1s(1,:) = 1s_2s, C0, C1, C2, C3
!  real(kind=dp), Cl_1s(5, 4), Cl_2s(3, 4), Cl_2p(3,4)
!  data Cl_1s(1,:) / 2.302d-1, 5.248d-6, -1.144d-10, 8.24d-16 /
 
 SUBROUTINE Giovanardi_et_al()
  

 RETURN
 END SUBROUTINE Giovanardi_et_al

END MODULE IMPACT