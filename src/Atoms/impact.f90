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
 !   and also:
 !    E.S Chang, E.H. Avrett, and R. Loeser 1991,
 !    A&A 247, 580-583
 !
 ! Other Elements
 !
 ! TODO: Integration of the cross-section itself.
! ------------------------------------------------------------------------------------- !

 use atmos_type, only : ne, T, Atoms, Natom, Hydrogen
 use atom_type!, only : AtomicLine, AtomicContinuum, AtomType
 use mcfost_env, only : dp
 use constant
 !!use special_functions
 use math
 use utils, only : interp_sp
 use messages, only : error, warning
 use parametres, only : ldissolve
 use occupation_probability, only : wocc_n

 IMPLICIT NONE

 !abscissa and weights for Gauss Laguerre integration of function int_0^inf f(x)*exp(-x) dx
 real(kind=dp), dimension(8) :: x_laguerre, w_laguerre
 !integ = sum(y(x_laguerre)*w_laguerre) with y = f(x)*exp(-x)
 data x_laguerre / 0.170279632305d0,  0.903701776799d0,  2.251086629866d0,  4.266700170288d0, &
                   7.045905402393d0, 10.758516010181d0, 15.740678641278d0, 22.863131736889d0 /

 data w_laguerre / 3.69188589342d-01, 4.18786780814d-01, 1.75794986637d-01, 3.33434922612d-02, &
                   2.79453623523d-03, 9.07650877336d-05, 8.48574671627d-07, 1.04800117487d-09 /

 !Impact Parameter approximation of Seaton 1962 vol 72
 real, dimension(24) :: beta, zeta, phi
 data beta  / 0.0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.4, 0.5, 0.6, 0.7, &
              0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0 /
 !these are two functions of beta
 data zeta  / 1.0, 1.015, 1.030, 1.035, 1.035, 1.025, 1.010, 0.963, 0.9, 0.829, 0.754, &
 			  0.68, 0.608, 0.540, 0.476, 0.418, 0.365, 0.318, 0.276, 0.239, 0.206, 0.177, 0.152, 0.130 /
              
 data phi  / 0.0, 3.1, 2.392, 1.972, 1.674, 1.444, 1.259, 0.974, 0.766, 0.608, 0.486, 0.390, &
 			 0.314, 0.253, 0.2049, 0.1661, 0.1347, 0.1095, 0.0889, 0.0723, 0.0589, 0.0480, &
 			 0.0390, 0.0319 /

 CONTAINS

 !Check addition of Collision ! probably a missing factor<
 !and the rates from j->i missing only i->j and i->ionisation (not recombination) present yet

 FUNCTION collision_atom(icell, atom) result(Cji)
  integer :: i, j, icell
  type (AtomType) :: atom
  real(kind=dp), dimension(atom%Nlevel, atom%Nlevel) :: Cji

   Cji(:,:) = 0d0

 RETURN
 END FUNCTION collision_atom

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
  real(kind=dp) :: C2_ions, C2_neutrals, dE, C_bf, deltam, R0, beta0, ri, rj, neff, Qij, Wi, Qij_strong
  logical :: forbidden


   deltam = 1. + M_ELECTRON / (atom%weight * AMU)
   !xions = 0.
   !xneutrals = 0.68
   C2_ions = 3.96d-6 !s^-1 m^3 K^3/2
   C2_neutrals = 2.15d-6
   C_bf = 1.55d11 ! s^-1 K^1/2 m^1

   C(:,:) = 0d0

   do i=1,atom%Nlevel !between pair of levels resulting in a b-b transition
    Wi = T(icell) * KBOLTZMANN
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

            C(i,j) = correct * C2_neutrals * T(icell)**(-1.5) * exp(-dE) * line%fosc * &
                         dE**(-1.68) !power of 1+xneutrals
          else !optically allowed: could also use Van Regemorter with correct = 1.
          ! Seaton Impact Paramater approximations !
          ! min of (Qij_weak, Qij_strong)

            ic = find_continuum(atom, i)
		    if (ic==0) CALL Error ("Impact")
            neff = (atom%stage(i)+1) * sqrt(E_RYDBERG/deltam / (atom%E(ic)-atom%E(i)))
            write(*,*) "(n,l,Z)", neff, atom%Lorbit(i), atom%stage(i)+1
            ri = atomic_orbital_radius(neff, atom%Lorbit(i),atom%stage(i)+1)
!             ic = locate(atom%stage*1d0, 1d0*atom%stage(j)+1d0)
!             if (ic == j) then
!              write(*,*) "Error, couldn't find continuum for level", atom%label(j)
!             end if
		    ic = find_continuum(atom, j)
		    if (ic==0) CALL Error ("Impact")
            neff = (atom%stage(j)+1) * sqrt(E_RYDBERG/deltam / (atom%E(ic)-atom%E(j)))
            write(*,*) "(n,l,Z)", neff, atom%Lorbit(j), atom%stage(j)+1
            ri = atomic_orbital_radius(neff, atom%Lorbit(j),atom%stage(j)+1)
            R0 = min(ri, rj)
            !! or neff = min(neff_i, neff_u) and R0 = 5*neff*2 + 1 / 4.
            beta0 = R0 * sqrt(Wi*deltam/E_RYDBERG) * dE
            write(*,*) "R0 (a0) = ", R0," beta0 = ", real(beta0), beta0
            !in the weak coupling
            Qij = 4d0 * pia0squarex2 * E_RYDBERG**2 * &
            	line%fosc*interp_sp(phi,beta, real(beta0)) / deltam / deltam / dE / Wi**2

            !Estimate of Qij_strong. The real calculation involved to solve for beta1 un beta1**2 * A = f(beta1)
            Qij_Strong = 2d0 * pia0squarex2 * (E_RYDBERG/deltam)**2 * line%fosc / dE/ Wi
            C(i,j) = min(Qij, Qij_strong)

          end if !forbidden neutral stage
       else !ions !Van Regemorter
          C(i,j) = correct * C2_ions * T(icell)**(-1.5) * exp(-dE) * line%fosc * dE
       endif

      endif !eq. 34

    enddo
   enddo

   !now bound-free eq. 26
   !or eq. 9.60 Hubeny Mihalas
   do kr=1,atom%Ncont !only over continua
    cont = atom%continua(kr)
    i = cont%i; j = cont%j
    dE = (atom%E(j) - atom%E(i)) / (KBOLTZMANN * T(icell))

    select case (atom%stage(i))
     case (0)
      gbar_i = 0.1 ! stage neutral
     case (1)
      gbar_i = 0.2 ! Z = 2
     case default
      gbar_i = 0.3 ! Z > 2
    end select

    C(i,j) = C_bf * T(icell)**(-0.5) * cont%alpha0 * gbar_i * exp(-dE) / dE

   enddo

 RETURN
 END SUBROUTINE Collision_with_electrons

 FUNCTION Psi_Drawin(u)
  real(kind=dp) :: Psi_Drawin, u

  Psi_Drawin = exp(-u) / (1d0 + 2d0/u)

 RETURN
 END FUNCTION Psi_Drawin

 FUNCTION Psi_Drawin_ionisation(u)
  real(kind=dp) :: Psi_Drawin_ionisation, u

  Psi_Drawin_ionisation = (1d0 + 2d0/u) * exp(-u)
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
  C0 = 8* pia0squarex2 * (2d0 * KBOLTZMANN * T(icell) / PI / mu)**(0.5)
        !16 * pi * RBHOR**2
  do kr=1, atom%Ntr

   select case (atom%at(kr)%trtype)

    case ("ATOMIC_LINE") !collisional excitation
     line = atom%lines(atom%at(kr)%ik)
     i = line%i; j = line%j
     dE = atom%E(j) - atom%E(i)
     u0 = dE / KBOLTZMANN / T(icell)
     f = line%fosc
     psi = Psi_drawin(u0)

    case ("ATOMIC_CONTNUUM") !collisional ionisation
     cont = atom%continua(atom%at(kr)%ik)
     i = cont%i; j = cont%j
     dE = atom%E(j) - atom%E(i)
     f = real(i)/2. ! HERE IT IS NOT CORRECT
     u0 = dE / KBOLTZMANN / T(icell)
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

 FUNCTION Collision_Hydrogen(icell) result(Cij)
  !Matrix containing the collision rates of hydrogen C(i,j) = collrate from i->j
  !This is not the collision rate matrix!
  integer :: icell, i, j
  real(kind=dp) :: Cij(Hydrogen%Nlevel, Hydrogen%Nlevel)
  real(kind=dp) :: nr_ji, CI(hydrogen%Nlevel), CE(Hydrogen%Nlevel,Hydrogen%Nlevel)!, CI(Hydrogen%Nlevel,Hydrogen%Nlevel)
  real(kind=dp) :: wi, wj

   Cij(:,:) = 0d0; CI = 0d0; CE(:,:) = 0d0
	wj = 1.0
	wi = 1.0
! 	if (ldissolve) then
! 		if (atom%ID=="H") then
! 												!nn
! 			wi = wocc_n(icell, real(i,kind=dp), real(atom%stage(i)), real(atom%stage(i)+1))
! 			wj = wocc_n(icell, real(j,kind=dp), real(atom%stage(j)), real(atom%stage(j)+1))
! 			
! 		endif
! 	endif
	
   CALL Johnson_CI(icell, CI) !bound-free i->Nlevel
   CALL Johnson_CE(icell, CE) !among all levels
   
   do i=1,Hydrogen%Nlevel
    Cij(i,hydrogen%Nlevel) = CI(i)
    Cij(hydrogen%Nlevel,i) = CI(i) * hydrogen%nstar(i,icell) / hydrogen%nstar(hydrogen%Nlevel,icell)
    do j=i+1,Hydrogen%Nlevel-1
    	!write(*,*), i, j, CE(i,j),CE(i,j) * hydrogen%nstar(i,icell)/hydrogen%nstar(j,icell)
 		Cij(i,j) = Cij(i,j) + CE(i,j)
 		Cij(j,i) = Cij(j,i) + CE(i,j) * hydrogen%nstar(i,icell)/hydrogen%nstar(j,icell)
    end do
   end do

   Cij(:,:) = Cij(:,:) * ne(icell)

!    if (T(icell)==T(55)) then
!    	do i=1, hydrogen%nlevel
!    		write(*,*) "****",i, "****"
!    		!write(*,*) (j, Cij(i,j), j=1,hydrogen%nlevel)
!    		do j=i+1, hydrogen%nlevel-1
!    			write(*,*) i,j, CE(i,j)
!    		enddo
!    	enddo
!    	stop
!    endif

 RETURN
 END FUNCTION Collision_Hydrogen


 SUBROUTINE Johnson_CI(icell, Cik)
 ! --------------------------------------------------- !
  ! Ionisation rate coefficient for Hydrogen atom
  ! from lower level i to the continuum level k
  ! C(i,k) = C_{i->k}
  !
  ! from L.C Johnson
  ! ApJ 74:227-236, 1972 May 15; eq. 39
  !
  ! ne factorised
  !
  ! return C(i,k) with k = Nlevel (bound-free)
 ! --------------------------------------------------- !
   integer, intent(in) :: icell
   real(kind=dp), intent(out), dimension(:) :: Cik
   integer :: i, j, Nl
   !real(kind=dp) :: x0 = 1 - (n/n0)**2 with n0 -> infinity
   real(kind=dp) :: C0, rn, bn, n, An, En, yn, zn, S, Bnp, deltaM

   deltam = 1. + M_ELECTRON/ (hydrogen%weight * AMU)

   C0 = sqrt(8.*KBOLTZMANN*T(icell) / pi / M_ELECTRON)
   
   Cik(:) = 0.0_dp

   Nl = hydrogen%Nlevel
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
    En = E_RYDBERG /deltam / n / n !energie of level with different quantum number in 13.6eV: En = 13.6/n**2
    yn = En / KBOLTZMANN / T(icell)
!     An = 32. / 3. / sqrt(3d0) / pi * n  * (g0(n)/3. + g1(n)/4. + g2(n)/5.)
    An = 32. / 3. / sqrt(3d0) / pi * n  * (g0(i)/3. + g1(i)/4. + g2(i)/5.)
    Bnp = 2./3. * n*n * (5. + bn)
    zn = rn + yn

    S = C0 * pia0squarex2 * (n*yn)**2 * (An*(E1(yn)/yn - E1(zn)/zn) + &
   			(Bnp - An*log(2*n*n))*(ksi_johnson(yn)-ksi_johnson(zn)))
   	!write(*,*) icell, i, "S=", S
    !check that otherwise multiply by exp(yn)
    Cik(i) = S !RH -> exp(yn) / sqrt(T(icell)) !per ne
    		   !we compute it at icell so in fact we could avoid / sqrt(icell)
    		   !as col = Cje * sqrt(T) * exp(-de). Normally dE == En/kT=y so
    		   !it is not useful to multiply except to take into account the slightly diff
    		   !between En and (hydrogen%E(Nl)-hydrogen%E(i))
    !!write(*,*) En, (hydrogen%E(Nl)-hydrogen%E(i)) !Should be similar, because E(j)=13.6/j**2
   end do

 RETURN
 END SUBROUTINE Johnson_CI


 SUBROUTINE Johnson_CE(icell, Cij)
 ! ----------------------------------------------------- !
  ! Excitation rate coefficient for
  ! Hydrogen atom, from
  ! ApJ 74:227-236, 1972 May 15; eq. 36
 ! ----------------------------------------------------- !
   integer, intent(in) :: icell
   real(kind=dp), intent(out), dimension(:,:) :: Cij
   integer :: i, j, Nl
   real(kind=dp) :: C0, rn, bn, n, Ennp, y, z, S, Bnnp, En
   real(kind=dp) :: np, x, fnnp, rnnp, Annp, Gaunt_bf, deltam

   deltam = 1. + M_ELECTRON/ (hydrogen%weight * AMU)

   C0 = sqrt(8.*KBOLTZMANN*T(icell) / pi / M_ELECTRON)

   Nl = hydrogen%Nlevel
   !Hydrogen level are ordered by n increasing, except for the continuum level
   !n = 1., but stops before Nlevel

   do i=1, Nl-1 !collision between neutral states, n to n'
    n = real(i,kind=dp)
    if (i==1) then !n=i
     rn = 0.45
     bn = -0.603
    else
     rn = 1.94/(i**1.57)
     bn = (1.0/n) * (4.0 - 18.63/n + 36.24/n/n - 28.09/n/n/n)
    end if

    do j=i+1, Nl-1
     np = real(j,kind=dp)!n'
     x = 1d0 - (n/np)**2 ! = Enn'/Rdybg
!      Gaunt_bf = g0(n) + g1(n)/x + g2(n)/x/x
     Gaunt_bf = g0(i) + g1(i)/x + g2(i)/x/x

     fnnp = Gaunt_bf * 32.0/ (3.0 * sqrt(3d0) * pi) * n / np / np /np / x / x / x
     rnnp = rn * x
     Annp = 2.0 * n*n*fnnp/x
    ! in Joules
     En = E_RYDBERG /deltam / n / n !energie of level with different quantum number in 13.6eV = ionisation E of n
     y = x * En / KBOLTZMANN / T(icell) !x = ratio of E/En
     Bnnp = 4d0 * (n**4)/(np**3) / x / x * (1. + 4./(3.0 * x) + bn/x/x)
     z = rnnp + y

     S = C0 * pia0squarex2 * n*n*y*y/x * (Annp*((1./y + 0.5)*E1(y)-(1./z + 0.5)*E1(z))+&
     	(Bnnp-Annp*dlog(2*n*n/x))*(E2(y)/y - E2(z)/z))

     Cij(i,j) = S
!      if (T(icell)==T(55)) then
!      	write(*,*) i, j
!      	write(*,*) n, np, S
!      	stop
!      endif
    end do !over np
   end do  !over n

 RETURN
 END SUBROUTINE Johnson_CE

 FUNCTION ksi_johnson(t) result(y)
  double precision :: t, y
  !E0
  y = exp(-t)/t - 2d0*E1(t) + E2(t)

 RETURN
 END FUNCTION ksi_johnson

 FUNCTION g0 (n) result(g)
  real(kind=dp) :: g
  integer :: n

  SELECT CASE (n)
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
  real(kind=dp) :: g
  integer :: n

  SELECT CASE (n)
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
  real(kind=dp) :: g
  integer :: n

  SELECT CASE (n)
  CASE (1)
   g = 0.07014
  CASE (2)
   g = 0.02947
  CASE DEFAULT
   g = 1.0/(n*n) * (0.3887 - 1.181 / n + 1.470 / n / n)
  END SELECT

 RETURN
 END FUNCTION g2

 SUBROUTINE Scholz_et_al(icell, C)
  !define the output, and check results
  integer, intent(in) :: icell
  real(kind=dp), intent(inout) :: C !to define

  real(kind=dp) :: gamma_1s_2s, gamma_1s_2p, Omega_2s, Omega_2p
  real(kind=dp) :: b(8), cc(6), x, dE_kT
  real(kind=dp), parameter :: expdE = exp(-7.5d-1), gi = 2d0 !g value for ground-state of Hydrogen is 2

  if (T(icell) < 5d3 .or. T(icell) > 5d5) RETURN

  data b / 4.5168d-02,  2.8056d+01,  7.2945d+00, 2.4805d-01,  1.0044d-01, -1.1143d-02,  &
          -1.3432d-03,  3.7570d-04  /

  data cc / 3.6177d-01,  1.3891d+00,  5.0866d-01, -3.8011d-01,  1.0158d-01, -1.0072d-02 /

  x = (KBOLTZMANN * T(icell)) / E_RYDBERG / (1d0 + M_ELECTRON/hydrogen%weight*AMU)
  dE_kT = E_RYDBERG
  !natural log
  gamma_1s_2s = b(1) * log(b(2)*x) * exp(-b(3)*x) + b(4) + b(5)*x + b(6)*x*x * &
  				b(7)*x*x*x + b(8)*x*x*x*x
  				 !c0   !c1*x     !c2*x**3   !c3*x**3      !c4*x**4        !c5*x**5
  gamma_1s_2p = cc(1) + cc(2)*x + cc(3)*x*x + cc(4)*x*x*x + cc(5)*x*x*x*x + cc(6)*x*x*x*x*x

  Omega_2s = 8.63d-6 * gamma_1s_2s * expdE / gi / sqrt(T(icell))
  Omega_2p = 8.63d-6 * gamma_1s_2p * expdE / gi / sqrt(T(icell))


 RETURN
 END SUBROUTINE Scholz_et_al


END MODULE IMPACT
