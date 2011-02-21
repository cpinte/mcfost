!-----------------------------------------------------------------------------
  subroutine WRITE_MCFOST
!-----------------------------------------------------------------------------
  use mcfost2ProDiMo
  use PARAMETERS,ONLY: Rin,Rout
  use GRID,ONLY: NXX,NZZ,Tdust,nmol
  use ProDiMo2mcfost,ONLY: varP,CIIpop,OIpop,deltaV_CII,deltaV_OI
  use Nature,ONLY: AU,pi
  implicit none
  integer :: i,ix,iz,iCII,iOI,SPindex
  real,dimension(n_rad,nz) :: nCII,nOI,delV_CII,delV_OI
  real :: CIIpopMF(2,n_rad,nz),OIpopMF(3,n_rad,nz)
  real :: x,z,Ntot,dV,rr(0:n_rad),dz,ProDiMoFIT

  iCII = SPindex("C+")
  iOI  = SPindex("O")
  varP(:,:) = nmol(iCII,:,:)            ! particle density of C+ [1/cm3]
  do ix=1,n_rad
    do iz=1,nz
      x = mcfost%r_grid(ix,iz)*AU
      z = mcfost%z_grid(ix,iz)*AU
      nCII(ix,iz) = ProDiMoFIT(x,z,.true.,.false.)
    enddo
  enddo
  rr(0) = Rin
  rr(n_rad) = Rout
  do ix=1,n_rad-1
    rr(ix) = sqrt(mcfost%r_grid(ix,1)*mcfost%r_grid(ix+1,1))*AU
  enddo
  Ntot = 0.0
  do ix=1,n_rad
    dz = (mcfost%z_grid(ix,2)-mcfost%z_grid(ix,1))*AU
    do iz=1,nz
      dV = pi*(rr(ix)**2-rr(ix-1)**2)*dz
      Ntot = Ntot + nCII(ix,iz)*dV
    enddo
  enddo
  Ntot = 2.0*Ntot   ! add lower side of the disk
  write(*,*) "N_CII = ",Ntot

  varP(:,:) = nmol(iOI,:,:)             ! particle density of O [1/cm3]
  do ix=1,n_rad
    do iz=1,nz
      x = mcfost%r_grid(ix,iz)*AU
      z = mcfost%z_grid(ix,iz)*AU
      nOI(ix,iz) = ProDiMoFIT(x,z,.true.,.false.)
    enddo
  enddo
  varP(:,:) = deltaV_CII(:,:)           ! line broadening parameter [km/s]
  do ix=1,n_rad
    do iz=1,nz
      x = mcfost%r_grid(ix,iz)*AU
      z = mcfost%z_grid(ix,iz)*AU
      delV_CII(ix,iz) = ProDiMoFIT(x,z,.false.,.false.)
    enddo
  enddo
  varP(:,:) = deltaV_OI(:,:)            ! line broadening parameter [km/s]
  do ix=1,n_rad
    do iz=1,nz
      x = mcfost%r_grid(ix,iz)*AU
      z = mcfost%z_grid(ix,iz)*AU
      delV_OI(ix,iz) = ProDiMoFIT(x,z,.false.,.false.)
    enddo
  enddo
  do i=1,2
    varP(:,:) = CIIpop(i,:,:)           ! relative pop. number level i
    do ix=1,n_rad
      do iz=1,nz
        x = mcfost%r_grid(ix,iz)*AU
        z = mcfost%z_grid(ix,iz)*AU
        CIIpopMF(i,ix,iz) = ProDiMoFIT(x,z,.false.,.false.)
      enddo
    enddo
  enddo
  do i=1,3
    varP(:,:) = OIpop(i,:,:)           ! relative pop. number level i
    do ix=1,n_rad
      do iz=1,nz
        x = mcfost%r_grid(ix,iz)*AU
        z = mcfost%z_grid(ix,iz)*AU
        OIpopMF(i,ix,iz) = ProDiMoFIT(x,z,.false.,.false.)
      enddo
    enddo
  enddo

  open(unit=12,file='forMCFOST.dat',status='replace')
  write(12,*) n_rad,nz
  write(12,'(2a4,99(a16))') "ix","iz","r/AU","z/AU", &
          & "nCII[cm-3]","dvCII[km/s]","CIIpop_1","CIIpop_2", &
          & "nOI[cm-3]","dvOI[km/s]","OIpop_1","OIpop_2","OIpop_3"
  do ix=1,n_rad
    do iz=1,nz
      x = mcfost%r_grid(ix,iz)
      z = mcfost%z_grid(ix,iz)
      write(12,'(2(i4),99(1pE16.7))') ix,iz,x,z, &
                & nCII(ix,iz),delV_CII(ix,iz),CIIpopMF(1:2,ix,iz), &
                & nOI(ix,iz), delV_OI(ix,iz), OIpopMF(1:3,ix,iz)
    enddo
  enddo
  close(12)

  end subroutine WRITE_MCFOST
