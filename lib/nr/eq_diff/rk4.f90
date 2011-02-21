SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
! Given values for the N variables y and their derivatives dydx known at x, use
! the fourthorder Runge-Kutta method to advance the solution over an interval 
! h and return the incremented variables as yout, which need not be a distinct 
! array from y. y, dydx and yout are all of length N. The user supplies the 
! subroutine derivs(x,y,dydx), which returns derivatives dydx at x.
  USE nrtype; USE nrutil, ONLY : assert_eq
  IMPLICIT NONE
  REAL(SP), DIMENSION(:), INTENT(IN) :: y,dydx
  REAL(SP), INTENT(IN) :: x,h
  REAL(SP), DIMENSION(:), INTENT(OUT) :: yout
  INTERFACE
     SUBROUTINE derivs(x,y,dydx)
       USE nrtype
       IMPLICIT NONE
       REAL(SP), INTENT(IN) :: x
       REAL(SP), DIMENSION(:), INTENT(IN) :: y
       REAL(SP), DIMENSION(:), INTENT(OUT) :: dydx
     END SUBROUTINE derivs
  END INTERFACE
  INTEGER(I4B) :: ndum
  REAL(SP) :: h6,hh,xh
  REAL(SP), DIMENSION(size(y)) :: dym,dyt,yt
  ! assert sert a verifier que les sizes sont egales
  ndum=assert_eq(size(y),size(dydx),size(yout),'rk4')
  hh=h*0.5_sp
  h6=h/6.0_sp
  xh=x+hh
  yt=y+hh*dydx
  call derivs(xh,yt,dyt)
  yt=y+hh*dyt
  call derivs(xh,yt,dym)
  yt=y+h*dym
  dym=dyt+dym
  call derivs(x+h,yt,dyt)
  yout=y+h6*(dydx+dyt+2.0_sp*dym)
END SUBROUTINE rk4
