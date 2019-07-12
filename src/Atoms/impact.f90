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