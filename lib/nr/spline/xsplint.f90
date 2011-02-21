PROGRAM xsplint
  !	driver for routine splint, which calls spline
  USE nrtype; USE nrutil
  USE nr
  IMPLICIT NONE
  INTEGER(I4B), PARAMETER :: NP=10
  INTEGER(I4B) :: i,nfunc
  REAL(SP) :: f,x,y,yp1,ypn
  REAL(SP), DIMENSION(NP) :: xa,ya,y2
  do nfunc=1,2
     if (nfunc == 1) then
        write(*,*) 'sinus'

        xa(1:NP)=arth(1,1,NP)*PI/NP
        ya(:)=sin(xa(:))
        yp1=cos(xa(1))
        ypn=cos(xa(NP))
     else if (nfunc == 2) then
        write(*,*) 'exp'
        xa(1:NP)=arth(1.0_sp,1.0_sp,NP)/NP
        ya(:)=exp(xa(:))
        yp1=exp(xa(1))
        ypn=exp(xa(NP))
     else
        stop
     end if
     !	call SPLINE to get second derivatives
     call spline(xa,ya,yp1,ypn,y2)
     !	call SPLINT for interpolations
     write(*,*) '        x         f(x)       interpol'
     do i=1,10
        if (nfunc == 1) then
           x=(-0.05_sp+i/10.0_sp)*PI
           f=sin(x)
        else if (nfunc == 2) then
           x=-0.05_sp+i/10.0_sp
           f=exp(x)
        end if
        y=splint(xa,ya,y2,x)
        write(*,'(1x,3f12.6)') x,f,y
     end do
     read(*,*)
  end do
END PROGRAM xsplint
