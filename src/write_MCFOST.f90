module ProDiMo2mcfost
  use GRID,ONLY: NXX,NZZ
  real :: varP(NXX,NZZ)              ! working variable for 2D interpolation
  integer,parameter :: MC_NSP=5      ! number of species
  integer,parameter :: MC_LEVMAX=10  ! maximum number of levels
  TYPE LPOP
    character(len=10) :: name(MC_NSP)       ! species name
    integer :: isp(MC_NSP)                  ! species indices
    integer :: lmax(MC_NSP)                 ! maximum level recorded
    real :: amass(MC_NSP)                   ! mass of atom/molecule
    real :: deltaV(MC_NSP,NXX,NZZ)          ! turbulent+thermal line widths
    real :: pops(MC_NSP,MC_LEVMAX,NXX,NZZ)  ! level populations
  END TYPE LPOP
  TYPE(LPOP) :: lpops 
end module

!-----------------------------------------------------------------------
  SUBROUTINE INIT_MCFOST 
!-----------------------------------------------------------------------
    use ProDiMo2mcfost,ONLY: MC_NSP,lpops
    use Nature,ONLY: amu
    integer :: i,SPindex
    lpops%name(1) = "C+"
    lpops%name(2) = "O"
    lpops%name(3) = "CO"
    lpops%name(4) = "H2O"
    lpops%name(5) = "H2O"
    do i=1,MC_NSP
      lpops%isp(i) = SPindex(lpops%name(i))
    enddo
    lpops%name(4) = "oH2O"
    lpops%name(5) = "pH2O"  
    lpops%amass(1) = 12.0*amu
    lpops%amass(2) = 16.0*amu
    lpops%amass(3) = 28.0*amu
    lpops%amass(4) = 18.0*amu
    lpops%amass(5) = 18.0*amu
    lpops%lmax(1) = 2                ! C+
    lpops%lmax(2) = 3                ! O
    lpops%lmax(3) = 10               ! CO
    lpops%lmax(4) = 10               ! oH2O
    lpops%lmax(5) = 9                ! pH2O
  end SUBROUTINE INIT_MCFOST
  
!-----------------------------------------------------------------------
  REAL FUNCTION ProDiMoFIT(x,z,elog,extra) 
!-----------------------------------------------------------------------
    use GRID,ONLY: Nx=>NXX,Nz=>NZZ,xx,zz
    use ProDiMo2mcfost,ONLY: var=>varP
    use NATURE,ONLY: AU
    implicit none
    real,intent(in) :: x,z
    logical,intent(in) :: elog,extra
    real :: zpos,f11,f12,f21,f22,x1,x2,z1,z2
    real,save :: xold=-1.0,zold=-1.0,xfac,zfac
    integer,save :: ix=1,iz=1
    
    if ((x.ne.xold).or.(z.ne.zold)) then
      xold = x
      zold = z
      !---------------------------------------------------------------
      ! ***  search for index ix to achieve xx(ix) <= x < x(ix+1)  ***
      !---------------------------------------------------------------
      do
        if ((xx(ix).le.x).or.(ix.le.1)) exit
        ix=ix-1
      enddo
      do
        if ((xx(ix+1).gt.x).or.(ix.ge.Nx-1)) exit
        ix=ix+1
      enddo
      x1 = LOG(xx(ix))
      x2 = LOG(xx(ix+1))
      xfac = (LOG(x)-x1)/(x2-x1)           ! x=x1 => xfac=0
      !write(*,'(i4,4(0pF10.5))') ix,EXP(x1)/AU,x/AU,EXP(x2)/AU,xfac

      !------------------------------------------------------------------------
      ! ***  search for index iz to achieve zpos(x,iz) <= z < zpos(x,iz+1)  ***
      !------------------------------------------------------------------------
      do
        zpos = zz(ix,iz) + xfac*(zz(ix+1,iz)-zz(ix,iz))
        if ((zpos.le.z).or.(iz.le.1)) exit
        iz=iz-1
      enddo
      do
        zpos = zz(ix,iz+1) + xfac*(zz(ix+1,iz+1)-zz(ix,iz+1))
        if ((zpos.gt.z).or.(iz.ge.Nz-1)) exit
        iz=iz+1
      enddo
      z1 = zz(ix,iz)   + xfac*(zz(ix+1,iz)  -zz(ix,iz)  )
      z2 = zz(ix,iz+1) + xfac*(zz(ix+1,iz+1)-zz(ix,iz+1))
      zfac = (z-z1)/(z2-z1)                ! z=z1 => zfac=0
      !write(*,'(i4,4(0pF10.5))') iz,z1/AU,z/AU,z2/AU,zfac
      !stop

      !---------------------------------
      ! ***  disable extrapolation?  ***
      !---------------------------------
      if(.not.extra) then
        xfac = MAX(0.0,MIN(1.0,xfac))
        zfac = MAX(0.0,MIN(1.0,zfac))
      endif  

    endif  
      
    f11 = var(ix,iz)
    f12 = var(ix,iz+1)
    f21 = var(ix+1,iz)
    f22 = var(ix+1,iz+1)
    if ((.not.elog).or.(f11.le.0.0).or.(f12.le.0.0) &
                 & .or.(f21.le.0.0).or.(f22.le.0.0)) then
      ! --- linear interpolation --- 
      ProDiMoFIT = f11*(1.0-xfac)*(1.0-zfac) + f12*(1.0-xfac)*zfac &
      &          + f21*   xfac   *(1.0-zfac) + f22*   xfac   *zfac
    else 
      ! --- log interpolation --- 
      f11 = LOG(f11)
      f12 = LOG(f12)
      f21 = LOG(f21)
      f22 = LOG(f22)
      ProDiMoFIT = f11*(1.0-xfac)*(1.0-zfac) + f12*(1.0-xfac)*zfac &
      &          + f21*   xfac   *(1.0-zfac) + f22*   xfac   *zfac
      ProDiMoFIT = EXP(ProDiMoFIT)
    endif

    !write(*,*) xfac,zfac
    !write(*,'(1pE10.3,10x,1pE10.3)') var(ix,iz+1),var(ix+1,iz+1)
    !write(*,'(10x,1pE10.3)') ProDiMoFIT
    !write(*,'(1pE10.3,10x,1pE10.3)') var(ix,iz),var(ix+1,iz)

  end FUNCTION ProDiMoFIT


!-----------------------------------------------------------------------------
  subroutine WRITE_MCFOST
!-----------------------------------------------------------------------------
  use mcfost2ProDiMo
  use ProDiMo2mcfost,ONLY: MC_NSP,MC_LEVMAX,varP,lpops
  use PARAMETERS,ONLY: Rin,Rout
  use GRID,ONLY: NXX,NZZ,Tgas,Tdust,nmol
  use Nature,ONLY: AU,pi,Msun
  implicit none
  integer :: i,l,ix,iz
  real,dimension(MC_NSP,n_rad,nz) :: nsp,delV
  real,dimension(MC_NSP,MC_LEVMAX,n_rad,nz) :: pops
  real :: x,z,dV,rr(0:n_rad),dz,fortho,fpara,ProDiMoFIT
  real :: Ntot(MC_NSP)

  ! for the fits interface
  integer :: status,unit,blocksize,bitpix,naxis
  integer, dimension(5) :: naxes
  logical :: simple, extend
  character(len=512) :: filename

  do i=1,MC_NSP

    varP(:,:) = nmol(lpops%isp(i),:,:)       ! particle density [1/cm^3]
    if (lpops%name(i).eq."oH2O") then
      do ix=1,NXX
        do iz=1,NZZ
          call H2O_ORTHO_PARA(Tgas(ix,iz),fortho,fpara) 
          varP(ix,iz) = fortho*varP(ix,iz)
        enddo
      enddo  
    else if (lpops%name(i).eq."pH2O") then
      do ix=1,NXX
        do iz=1,NZZ
          call H2O_ORTHO_PARA(Tgas(ix,iz),fortho,fpara) 
          varP(ix,iz) = fpara*varP(ix,iz)
        enddo
      enddo  
    endif   
    do ix=1,n_rad
      do iz=1,nz
        x = mcfost%r_grid(ix,iz)*AU
        z = mcfost%z_grid(ix,iz)*AU
        nsp(i,ix,iz) = ProDiMoFIT(x,z,.true.,.false.) 
      enddo
    enddo  
    rr(0) = Rin
    rr(n_rad) = Rout
    do ix=1,n_rad-1
      rr(ix) = sqrt(mcfost%r_grid(ix,1)*mcfost%r_grid(ix+1,1))*AU  
    enddo
    Ntot(i) = 0.0
    do ix=1,n_rad
      dz = (mcfost%z_grid(ix,2)-mcfost%z_grid(ix,1))*AU 
      do iz=1,nz
        dV = pi*(rr(ix)**2-rr(ix-1)**2)*dz 
        Ntot(i) = Ntot(i) + nsp(i,ix,iz)*dV 
      enddo
    enddo  
    Ntot(i) = 2.0*Ntot(i)   ! add lower side of the disk
    varP(:,:) = lpops%deltaV(i,:,:)      ! line broadening parameter [km/s]
    do ix=1,n_rad
      do iz=1,nz
        x = mcfost%r_grid(ix,iz)*AU
        z = mcfost%z_grid(ix,iz)*AU
        delV(i,ix,iz) = ProDiMoFIT(x,z,.false.,.false.) 
      enddo
    enddo  

    do l=1,lpops%lmax(i)
      varP(:,:) = lpops%pops(i,l,:,:)    ! relative population numbers
      do ix=1,n_rad
        do iz=1,nz
          x = mcfost%r_grid(ix,iz)*AU
          z = mcfost%z_grid(ix,iz)*AU
          pops(i,l,ix,iz) = ProDiMoFIT(x,z,.false.,.false.) 
        enddo  
      enddo
    enddo  

  enddo  ! end of loop over output species


  ! New fits interface
  ! to be checked as I can't compile it and it is easy to make small 
  ! mistakes with cfitsio

   filename = "forMCFOST.fits.gz"
   ! Initialize parameters about the FITS image
   simple=.true.
   extend=.true.
   group=1
   fpixel=1

   !  Get an unused Logical Unit Number to use to open the FITS file.
   status=0
   call ftgiou(unit,status)

   !  Create the new empty FITS file.
   blocksize=1
   call ftinit(unit,trim(filename),blocksize,status)

  !------------------------------------------------------------------------------
  ! HDU 1 : grid, to make a test inside mcfost
  ! Format is the same as for the mcfost2ProDiMo interface
  !------------------------------------------------------------------------------
  bitpix=-64
  naxis=3
  naxes(1)=n_rad
  naxes(2)=nz
  naxes(3)=2
  nelements=naxes(1)*naxes(2)*naxes(3)

  !  Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  !  Write the array to the FITS file.
  do ix=1, n_rad
     grid(ix,:,1) = mcfost%r_grid(ix,:)
     do iz=1,nz
        grid(ix,iz,2) = mcfost%z_grid(ix,iz)
     enddo
  enddo
  call ftpprd(unit,group,fpixel,nelements,grid,status)

  !------------------------------------------------------------------------------
  ! HDU 2 :  particle density [1/cm^3]
  !------------------------------------------------------------------------------
  bitpix = -32
  naxis=3
  naxes(1)=MC_NSP
  naxes(2)=n_rad
  naxes(3)=nz
  nelements=naxes(1)*naxes(2)*naxes(3)

  ! create new hdu
  call ftcrhd(unit, status)

  ! Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,nsp,status)

  !------------------------------------------------------------------------------
  ! HDU 3 :  line broadening parameter [km/s]
  !------------------------------------------------------------------------------
  bitpix = -32
  naxis=3
  naxes(1)=MC_NSP
  naxes(2)=n_rad
  naxes(3)=nz
  nelements=naxes(1)*naxes(2)*naxes(3)

  ! create new hdu
  call ftcrhd(unit, status)

  ! Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,delV,status)

  !------------------------------------------------------------------------------
  ! HDU 4 : relative population numbers
  !------------------------------------------------------------------------------
  bitpix = -32
  naxis=3
  naxes(1)=MC_NSP
  naxes(2)=lpops%lmax(i)
  naxes(3)=n_rad
  naxes(4)=nz
  nelements=naxes(1)*naxes(2)*naxes(3)*naxes(4)

  ! create new hdu
  call ftcrhd(unit, status)

  ! Write the required header keywords.
  call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

  ! Write the array to the FITS file.
  call ftppre(unit,group,fpixel,nelements,pops,status)
  
  !  Close the file and free the unit number.
  call ftclos(unit, status)
  call ftfiou(unit, status)

  !  Check for any error, and if so print out error messages
  if (status > 0) then
     call print_error(status)
     stop
  end if
  
!$--  open(unit=12,file='forMCFOST.dat',status='replace')
!$--  write(12,*) n_rad,nz 
!$--  write(12,'(2a4,99(a16))') "ix","iz","r/AU","z/AU", &
!$--          & "nCII[cm^-3]","dvCII[km/s]","CIIpop_1","CIIpop_2", &
!$--          & "nOI[cm^-3]","dvOI[km/s]","OIpop_1","OIpop_2","OIpop_3", &
!$--          & "nCO[cm^-3]","dvCO[km/s]","COpop_1","COpop_2", "COpop_3", &
!$--          &              "COpop_4","COpop_5","COpop_6","COpop_7", &
!$--          &              "COpop_8","COpop_9","COpop_10", &
!$--          & "noH2O[cm^-3]","dvoH2O[km/s]","oH2Opop_1","oH2Opop_2","oH2Opop_3",&
!$--          &              "oH2Opop_4","oH2Opop_5","oH2Opop_6","oH2Opop_7", &
!$--          &              "oH2Opop_8","oH2Opop_9","oH2Opop_10", &
!$--          & "npH2O[cm^-3]","dvpH2O[km/s]","pH2Opop_1","pH2Opop_2","pH2Opop_3",&
!$--          &              "pH2Opop_4","pH2Opop_5","pH2Opop_6","pH2Opop_7", &
!$--          &              "pH2Opop_8","pH2Opop_9"
!$--  do ix=1,n_rad
!$--    do iz=1,nz
!$--      x = mcfost%r_grid(ix,iz)
!$--      z = mcfost%z_grid(ix,iz)
!$--      write(12,'(2(i4),99(1pE16.7))') ix,iz,x,z, &
!$--                & (   nsp(i,ix,iz), &
!$--                &    delV(i,ix,iz), &
!$--                    (pops(i,l,ix,iz),l=1,lpops%lmax(i)),i=1,MC_NSP)
!$--    enddo
!$--  enddo
!$--  write(12,*)
!$--  do i=1,MC_NSP
!$--    write(*,'(a10,"Ntotal=",1pE12.4,"   Mtot/Msun=",1pE12.4)') &
!$--            & lpops%name(i),Ntot(i),Ntot(i)*lpops%amass(i)/Msun
!$--    write(12,'(a10,"Ntotal=",1pE12.4,"   Mtot/Msun=",1pE12.4)') &
!$--            & lpops%name(i),Ntot(i),Ntot(i)*lpops%amass(i)/Msun
!$--  enddo  
!$--  close(12)

  end subroutine WRITE_MCFOST
