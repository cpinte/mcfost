! ----------------------------------------------------------------------
! Adapted from RH H. Uitenbroek
! ----------------------------------------------------------------------
! --- Reads collisional data from atomic data file and computes
! collisional rates. ---
!
!
!Input line can be either
!
!TEMP  Nitem  T[0]     ...   T[Nitem-1]
!
!to read in the temperature grid for collisional coefficients, or
!
!KEYWORD  i1  i2   coeff[0]    ...   coeff[Nitem-1]
!
!to read the coefficients for transitions between i1 and i2.
!Multiple entries of the temperature grid are allowed, but
!at least one entry with the correct number of grid points
!has to precede an entry with coefficients.
!
!Allowed keywords are:
!
!Keyword    Transition type
!----------------------------------------------------
!
!TEMP  -->  Temperature grid
!
!OMEGA -->  Collisional de-excitation of ions by electrons
!CE    -->  Collisional de-excitation of neutrals by electrons
!CI    -->  Collisional ionization by electrons
!CP    -->  Collisional de-excitation by protons
!
!CH0   -->  Charge exchange of ion with neutral hydrogen
!CH+   -->  Charge exchange of neutral with protons
!
! See also Hubeny Mihalas stellar atmospheres,
! chapter on collisional processes
!
! AR85-CEA  -->  Collisional auto-ionisation following Arnaud &
!Rothenflug (1985, ApJS 60)
!AR85-CDI  ---> Collisional ionisation following Arnaud & Rothenflug
!(1985, ApJS 60)
!AR85-CHP  -->  Charge exchange with ionised hydrogen following
!Arnaud & Rothenflug (1985, ApJS 60)
!AR85-CHH  -->  Charge exchange with neutral hydrogen following
!Arnaud & Rothenflug (1985, ApJS 60)
!BURGESS   -->  Collisional ionisation from excited states following
!Burgess & Chidichimo (1983, MNRAS 203, 1269)
!BADNELL   -->  Dielectronic recombination following the Badnell
!recipe (Badnell 2006)
!SHULL82   -->  Coefficients for collisional ionization, radiative
!recombination, and dielectronic recombination following
!Shull & van Steenberg (1982, ApJS, 48, 95)
!END   -->  End of input data
!----------------------------------------------------
!
!Note: C in Unit of number density is m^-3 and atom%Ckij and atom%C in s^-1.
!
!Convention: C_ij = C[i][j] represents the
!transition j --> i = Cul
! C_ij = C[ith ligne][jth column]
! fortran is column row
!     ij = (i-1)*atom%Nlevel +j : j-->i
!     ji = (j-1)*atom%Nlevel +i : i-->j
! ----------------------------------------------------------------------


MODULE collision
 use math, only : E1, E2, bezier3_interp, SQ, CUBE, dPOW, interp1D
 use constant
 use atom_type, only : AtomType, ATOM_ID_WIDTH, PTrowcol, CollisionData
 use readatom,only : path_to_atoms
 use atmos_type, only : atomZnumber, atmos, Hydrogen, Helium
 use getline
 use mcfost_env

 IMPLICIT NONE

 real(8), parameter :: BRK=4.0
 integer, parameter :: MSHELL=5
 character, parameter :: COMMENT_CHAR="#"


 CONTAINS

 FUNCTION fone(x) result(y)
 ! --------------------------------------------
 ! f_1 function of
 !  Arnaud & Rothenflug, 1985, A&ASS, 60, 425
 ! --------------------------------------------

  double precision :: x, y

  if (x.le.50.0) then
   y = exp(x)*E1(x)
  else
   y = 1./x
  end if

 RETURN
 END FUNCTION fone

 FUNCTION ftwo(x) result(y)
 ! --------------------------------------------
 ! f_2 function of
 !  Arnaud & Rothenflug, 1985, A&ASS, 60, 425
 !
 !  Improved description when x < 4 from:
 !    Hummer, 1983, jqsrt, 30 281
 ! --------------------------------------------

  integer :: i
  double precision :: x, y, p(15), q(15), px, xfact
  double precision :: qx, gamma, f0x, count, fact, term

  p = (/1.0000e+00, 2.1658e+02, 2.0336e+04, 1.0911e+06,&
        3.7114e+07, 8.3963e+08, 1.2889e+10, 1.3449e+11,&
        9.4002e+11, 4.2571e+12, 1.1743e+13, 1.7549e+13,&
        1.0806e+13, 4.9776e+11, 0.0000/)

  q = (/1.0000e+00, 2.1958e+02, 2.0984e+04, 1.1517e+06,&
        4.0349e+07, 9.4900e+08, 1.5345e+10, 1.7182e+11,&
        1.3249e+12, 6.9071e+12, 2.3531e+13, 4.9432e+13,&
        5.7760e+13, 3.0225e+13, 3.3641e+12/)

  if (x.gt.BRK) then
   px=p(1)
   xfact=1.
   do i=2,15
    xfact = xfact/x
    px = px+p(i)*xfact
   end do
   qx=q(1)
   xfact=1.
   do i=2,15
    xfact = xfact/x
    qx = qx+q(i)*xfact
   end do
   y = px/(qx*SQ(x))
  else
   gamma = 0.5772156649
   f0x = SQ(PI) / 12.
   term = 1.
   count = 0.0
   fact = 1.0
   xfact = 1.0

   do while (dabs(term/f0x).gt.1d-8)
    count = count+1
    fact = fact*count
    xfact = xfact*(-x)
    term = xfact/(SQ(count)*fact)
    f0x = f0x+term
    if (count.gt.100.0) then
     write(*,*) "Error in ftwo function in collision.f90"
     write(*,*) "code not stopped..."
    end if
   end do
   y = exp(x) * ((log(x) + gamma) * &
        (log(x) + gamma)*0.5 + f0x)
  end if
 RETURN
 END FUNCTION ftwo

 ! -----------------------------------------------
 ! REMEMBER -> stage is like in C, it goes from 0
 ! for neutrals, to Nstage-1 for the most ionised
 ! stage (instead of 1 to Nstage)
 ! Furhter, in readatom.f90, i is i+1 and j is j+1
 ! to be in aggreement with the index convention of
 ! fortran 90. First level is i(j) = 1, not 0 !
 ! and last level is i(j) = Nlevel, not Nlevel-1.
 ! -----------------------------------------------

 FUNCTION ar85cea(i,j,k, atom) result(cup)
 ! -----------------------------------------------
 ! Computes collisional autoionisation rates
 ! using the formalism from Arnaud and Rothenflug
 ! 1985, A&ASS, 60, 425 (a.k.a ar85)
 !
 ! -----------------------------------------------
  integer :: i, j, k, iz, ichrge, isoseq
  type (AtomType) :: atom
  double precision :: cup
  character(len=ATOM_ID_WIDTH) :: cseq
  double precision :: zz, bkt, b, zeff, iea, y
  double precision :: f1y,a,g,cea

  ! -- Initialisation --
  cea = 0.
  y = 0.
  f1y = 0.
  cup = 0.

  ! -- Search for element --
  iz = atomZnumber(atom%ID)!Hydrogen 1
  zz = iz
  if (iz.lt.1 .and. iz.gt.92) then
   write(*,*) "Limits of the routine in terms of ", &
    "elements exceeded, returing..."
   RETURN
  end if

  ! -- get iso-electronic sequence --
  ichrge = atom%stage(i) !stages from 0 for neutrals to N-1
  isoseq = iz-ichrge
  if (isoseq.lt.29) cseq = atmos%Elements(isoseq)%ID

  ! -- Temperature in eV --
  bkt = KBOLTZMANN * atmos%T(k) / EV

  ! All ID in atmos%Elements are of the form Xx
  ! just as in atom%ID
  ! -- Lithium sequence --
  if (cseq.eq."Li") then
   iea = 13.6*(dpow(zz-0.835d0,2.d0)-&
      0.25*dpow(zz-1.62d0,2.d0))
   b = 1./(1.+2.0d-4 * dpow(zz, 3.d0))
   zeff = zz-0.43
   y = iea/bkt
   f1y = fone(y)
   g = 2.22*f1y + 0.67*(1.-y*f1y)+0.49*y*f1y + &
     1.2*y*(1.-y*f1y)
   cup = (1.6d-7 * 1.2*b) / (dpow(zeff,2.d0)*&
      dsqrt(bkt)) * exp(-y)*g
   if (atom%ID.eq."C ") then
    ! -- C IV, app A ar85 --
    cup = cup*0.6
   else if (atom%ID.eq."N ") then
    ! -- N V, app A ar85 --
    cup = cup*0.8
   else if (atom%ID.eq."O ") then
    ! -- O VI, app A ar85 --
    cup = cup*1.25
   end if
  else if (cseq.eq."Na") then
   ! -- Sodium sequence --
   if (iz.le.16) then
    iea = 26.*(zz-10.)
    a = 2.9d-17 * dpow(zz-11.,-7d-1)
    y = iea/bkt
    f1y = fone(y)
    cup = 6.69d7 * a * iea / dsqrt(bkt) * &
         exp(-y) * (1.-y*f1y)
   else if (iz.ge.18 .and. iz.le.28) then
    iea = 11.*(zz-10.)*sqrt(zz-10.)
    a = 1.4e-14 * dpow(zz-10., -3.73d0)
    y = iea/bkt
    f1y = fone(y)
    cup = 6.69d7 * a * iea/dsqrt(bkt)  * exp(-y) * &
       (1.-0.5*(y-SQ(y)+SQ(y)*y*f1y))
   else
    cup = 0.
   end if
  end if
  ! -- Magnesium-sulfur sequences --
  if ((cseq.eq."Mg") .or. (cseq.eq."Al") .or. &
      (cseq.eq."Si") .or. (cseq.eq."P ") .or. &
      (cseq.eq."S ")) then

   if (cseq.eq."Mg") iea=10.3*dpow(zz-10.,1.52d0)
   if (cseq.eq."Al") iea=18.0*dpow(zz-11.,1.33d0)
   if (cseq.eq."Si") iea=18.4*dpow(zz-12.,1.36d0)
   if (cseq.eq."P ") iea=23.7 *dpow(zz-13., 1.29d0)
   if (cseq.eq."S ") iea=40.1*dpow(zz-14.,1.1d0)

   a = 4.d-13 / (SQ(zz)*iea)
   y = iea/bkt
   f1y = fone(y)
   cup = 6.69d7 *a*iea / dsqrt(bkt) * exp(-y) * &
     (1.-0.5*(y-SQ(y)+SQ(y)*y*f1y))

  end if
  ! -- Special cases --
  if ((atom%ID.eq."Ca").and.(ichrge.eq.0)) then
   iea = 25.
   a = 9.8d-17
   b = 1.12
   cup = 6.69d7 * a * iea / dsqrt(bkt) * exp(-y) * &
    (1. + b*f1y)
  else if ((atom%ID.eq."Ca").and.(ichrge.eq.1)) then
   iea = 25.0
   a = 6.0d-17
   b = 1.12
   cup = 6.69d7 * a * iea / dsqrt(bkt) * exp(-y) * &
    (1. + b*f1y)
  else if ((atom%ID.eq."Fe").and.(ichrge.eq.3)) then
   a = 1.8d-17
   iea = 60.0
   b = 1.0
   cup = 6.69d7 * a * iea / dsqrt(bkt) * exp(-y) * &
    (1. + b*f1y)
  else if ((atom%ID.eq."Fe").and.(ichrge.eq.4)) then
   a = 5.d-17
   iea = 73.
   b = 1.
   cup = 6.69d7 * a * iea / dsqrt(bkt) * exp(-y) * &
    (1.+b*f1y)
  end if

  cup = cup * CUBE(CM_TO_M) ! in m3

 RETURN
 END FUNCTION

 FUNCTION summers(i, j, nne, atom) result(y)
 ! -----------------------------------------------------
 ! --- Density sensitive dielectronic recombination ----
 ! the term adi may be multiplied by a density-sensitive
 ! factor if needed, this is crucial for Li and B-like
 ! ions collinding with impacting electrons.
 !
 ! This simple formulation was derived from a study of
 ! the dependence of the dielectronic "bump" in the figures
 ! of Summers 1974 and fitting according to Ne/Z^7
 ! Author, Uitenbroek, Judge, Leenaarts

  integer :: i, j, isoseq, row, col, iz
  double precision :: y, zz, rho0, rhoq, x, beta, nne
  type (AtomType) :: atom
  character(len=ATOM_ID_WIDTH) :: cseq

  iz = atomZnumber(atom%ID)
  if (iz.lt.1 .or. iz.gt.92) then
   write(*,*) "In summers, out of bound"
  end if

  ! -- charge of recombining ion --
  zz = atom%stage(j)
  isoseq = iz - atom%stage(i)
  if (isoseq.lt.29) cseq = atmos%Elements(isoseq)%ID

  !write(*,*) i, j, zz, isoseq, cseq

  ! -- find ROW and COLumn in periodic table --
  ! Note that for the three first lines we have the
  ! particular cases:
  ! - if Z=2 (He) returns, 1, 18 because
  !   He is on the first row at the 18th column
  !   if we look at the maximum number of columns
  ! - if Z>=5 and Z<=10 or Z >= 13 and Z<= 18
  !   we add +10 to the column number in the subroutine
  !  since they are 10 empty columns for theses rows in
  !  in the periodic table.
  !  if you want the relative positions:
  !   col(He) = col(He)-16 (He->2)
  !   col(Z3) = col(Z3)-10 (B->3)


  CALL PTrowcol(isoseq, row, col)
  if (isoseq.eq.2) col = col-16 ! Helium, relative col=2
  if ( ((isoseq.ge.5).and.(isoseq.le.10)).or.&
       ((isoseq.ge.13).and.(isoseq.le.18))) col=col-10
  !write(*,*) "row=",row, "col=",col

  rhoq = nne*CUBE(CM_TO_M) / dpow(zz, 7.d0)
  x = (0.5*zz + (col-1.)) * row/3.
  beta = -0.2 / log(x+2.71828)
  rho0 = 30.0 + 50.0*x;
  y = dpow(1. + rhoq/rho0,beta)

  !write(*,*) "rhoq=",rhoq," x=",x," beta=", beta," y=", y

 RETURN
 END FUNCTION summers
 
 SUBROUTINE openCollisionFile(atom)
  type (AtomType), intent(in) :: atom
  integer :: checkfseek, fseek
  
   !re-open atomic model, but set the cursor at the position
   !of collisional data.
  write(*,*) "Openning collision file for atom ", atom%ID
  open(unit=atom%colunit,file=trim(mcfost_utils)//path_to_atoms//trim(atom%inputFile),status="old")
  write(*,*) "Warning in CollisionRate "&
    ": fseek is different in gfortran and ifort"
  !gfortran ->
  !!CALL FSEEK(colunit, atom%offset_coll,0,checkfseek)
  !ifort ->
  checkfseek = fseek(atom%colunit,atom%offset_coll, 0)
  if (checkfseek.gt. 0 ) then 
    write(*,*) 'fseek error'
    stop 
  end if
  return
 END SUBROUTINE
 
 SUBROUTINE closeCollisionFile(atom)
  type (AtomType), intent(in) :: atom
  write(*,*) "Closing collision file for atom ", atom%ID, atom%colunit
  close(atom%colunit)
  return
 END SUBROUTINE closeCollisionFile
!  
!  SUBROUTINE initCollision(atom)
!   type (AtomType), intent(inout) :: atom
!   integer :: ij, k, Nread, countline=0, colunit, checkfseek, fseek
!   integer :: NTMP, n, m, ii, Nitem, i1, i2, i, j, ji, Nitem2, Nkey
!   character(len=8) :: END_OF_FILE="END     ", key
!   double precision :: C, np
!   real(8) :: deltaE, C0, Cdown, Cup, gij, xj, fac, fxj
!   integer :: Ncoef, Nrow
!   character(len=MAX_LENGTH) :: inputline, FormatLine
!   double precision, dimension(:), allocatable :: coeff
!   real(8), dimension(:,:), allocatable :: badi, cdi
!   real(8) :: acolsh, tcolsh, aradsh, xradsh, adish, bdish
!   real(8) :: t0sh, t1sh, summrs, tg, cdn, ccup
!   real(8) :: ar85t1, ar85t2, ar85a, ar85b, ar85c, ar85d, t4
!   real(8) :: de,zz,betab,cbar,dekt,dekti,wlog,wb,sumscl=0.0
! 
!   colunit = atom%colunit
!   write(*,*) "Collision()->colunit = ", colunit, atom%colunit
!   write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
!   
!   C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *&
!         sqrt(8.0/(PI*KBOLTZMANN))
! 
!   key="        "
!   Nkey = 0
!   ! -- Go throught the remaining lines in the file --
!   ! read collisional data depending the cases of recipes.
!   
!   !First count how many data to read
! 
!   do while(key.ne.END_OF_FILE)
!    countline = countline + 1
!    CALL getnextline(colunit, COMMENT_CHAR, &
!          FormatLine, inputline, Nread)
! 
!    key = adjustl(inputline(1:len(key)))
!    if (key.eq."TEMP") then
!     read(inputline(len(key)+1:len(key)+3), *) atom%col_data%NTMP
!     allocate(atom%col_data%TGRID(atom%col_data%NTMP))
! 
!     read(inputline(len(key)+3:Nread),*) (TGRID(k), k=1,NTMP)
!     !write(*,*) (TGRID(k), k=1,NTMP)
!     Nitem = atom%col_data%NTMP
!     Nitem2 = atom%col_data%TGRID
!     NTMP = atom%col_data%NTMP
! 
!    else if ((key.eq."OMEGA") .or. (key.eq."CE") .or.&
!         (key.eq."CI") .or. (key.eq."CP") .or. &
!         (key.eq."CH0") .or. (key.eq."CH+") .or. &
!         (key.eq."CH") .or. (key.eq."CR")) then
!     Nkey = Nkey + 1
!     if (key == "OMGEA") then
!      atom%col_data%ir(Nkey) = 1
!     else if (key == "CE") then
!      atom%col_data%ir(Nkey) = 2
!     else if (key == "CI") then
!      atom%col_data%ir(Nkey) = 3
!     else if (key == "CP") then
!      atom%col_data%ir(Nkey) = 4
!     else if (key == "CH0") then
!      atom%col_data%ir(Nkey) = 5
!     else if (key == "CH+") then
!      atom%col_data%ir(Nkey) = 6
!     else if (key == "CH") then
!      atom%col_data%ir(Nkey) = 7
!     else if (key == "CR") then
!      atom%col_data%ir(Nkey) = 8
!     else
!      write(*,*) Key, "collision recipe unknonw"
!      stop
!     end if 
!     ! -- read level indices and collision coefficients --
!     Nitem = NTMP
!     allocate(atom%col_data%coeffs_1(Nitem))
!     read(inputline(len(key)+1:),*) i1, i2, &
!          (coeff(k),k=1,Nitem)
! 
!     Nitem2 = size(atom%col_data%coeffs_1)
! 
!     i = MIN(i1,i2) + 1
!     j = MAX(i1,i2) + 1 !from 1 to Nlevel (in C, 0 to Nlevel-1)
!     ij = (i-1)*atom%Nlevel + j!j->i
!     ji = (j-1)*atom%Nlevel + i !i->j
!     atom%col_data%ij(atom%col_data%ir(Nkey)) = ij
!     atom%col_data%ji(atom%col_data%ir(Nkey)) = ij
!     ! -- Note:
!     ! Cf = C(Nlevel,Nlevel).flatten, of size Nlevel*Nlevel
!     ! C(i,j) = Cf(ij) = Cf(i*Nlevel + j)
!     ! C(1,1) = Cf(1) -> i*Nlevel +j = 1 -> i=i-1
! 
!    else if ((key.eq."AR85-CHP").or.(key.eq."AR85-CHH")) then
!     Nitem = 6
!     if (allocated(coeff)) deallocate(coeff)
!     allocate(coeff(Nitem))
!     read(inputline(len(key)+1:),*) i1, i2, &
!         (coeff(k),k=1,Nitem)
!     !write(*,*) inputline(len(key)+1:)
! 
!     Nitem2 = size(coeff)
! 
!     i = MIN(i1,i2) + 1
!     j = MAX(i1,i2) + 1
!     ij = (i-1)*atom%Nlevel + j
!     ji = (j-1)*atom%Nlevel + i
! 
!    else if ((key.eq.'AR85-CEA').or.(key.eq."BURGESS")) then
!      Nitem = 1
!      if (allocated(coeff)) deallocate(coeff)
!      allocate(coeff(Nitem))
! 
!      read(inputline(len(key)+1:),*) i1, i2, coeff(1)
!      i = MIN(i1,i2) + 1
!      j = MAX(i1,i2) + 1
!      ij = (i-1)*atom%Nlevel + j
!      ji = (j-1)*atom%Nlevel + i
! 
!      Nitem2 = size(coeff)
!      !write(*,*) "CEA or BURGESS", i1, i2, i, ij
! 
!    else if (key.eq."SHULL82") then
!      Nitem = 8
!      if (allocated(coeff)) deallocate(coeff)
!      allocate(coeff(Nitem))
! 
!      read(inputline(len(key)+1:),*) i1, i2,&
!          (coeff(k),k=1,Nitem)
!      i = MIN(i1,i2) + 1
!      j = MAX(i1,i2) + 1
!      ij = (i-1)*atom%Nlevel+j
!      ji = (j-1)*atom%Nlevel +i
! 
!      Nitem2 = size(coeff)
! 
! 
!    else if (key.eq."BADNELL") then
!     ! -- BADNELL formula for dielectronic recombination --
! 
!     !write(*,*) inputline, Nread
!     read(inputline(len(key)+1:),*) i1, i2, Ncoef
!     Nrow = 2
!     Nitem = Nrow*Ncoef
!     allocate(badi(Nrow,Ncoef))
!     ! they are two lines of Nitem elements to read
!     !write(*,*) i1, i2, Ncoef
!     do m=1,Nrow
!       CALL getnextline(colunit, COMMENT_CHAR, &
!              FormatLine, inputline, Nread)
!       !write(*,*) 'row=',m, inputline, Nread
!       read(inputline,*) (badi(m,k),k=1,Ncoef)
!     end do
!     i = MIN(i1,i2) + 1
!     j = MAX(i1,i2) + 1
!     ij = atom%Nlevel*(i-1) + j
!     ji = atom%Nlevel*(j-1) + i
! 
!     Nitem2 = size(badi(1,:)) * size(badi(:,1))
! 
!    else if (key.eq."SUMMERS") then
!     ! -- Switch for density dependent DR coefficient
!     !
!     ! Give the default multiplication factor of summers
!     ! density dependence of dielectronic recombination:
!     ! sumscl = 0. -> no density dependence
!     ! sumscl = 1. -> full density dependence
! 
!     Nitem = 1
!     read(inputline(len(key)+1:), *) sumscl
! 
!     Nitem2 = 1
! 
!    else if (key.eq."AR85-CDI") then
!     read(inputline(len(key)+1:),*) i1, i2, Nrow
! 
!     if (Nrow.gt.MSHELL) then
!      write(*,*) "Nrow is greater than MSHELL, exiting..."
!      stop
!     end if
!     Nitem = Nrow*MSHELL
!     allocate(cdi(Nrow, MSHELL))
!     do m=1,Nrow
!       CALL getnextline(colunit, COMMENT_CHAR, &
!        FormatLine, inputline, Nread)
!       read(inputline,*) (cdi(m,k),k=1,MSHELL)
!     end do
! 
!     Nitem2 = size(cdi(1,:)) * size(cdi(:,1))
! 
!     i = MIN(i1,i2) + 1
!     j = MAX(i1,i2) + 1
!     ij = (i-1)*atom%Nlevel +j
!     ji = (j-1)*atom%Nlevel +i
! 
!    else if (key.eq.END_OF_FILE) then
!     exit
!    else
!     write(*,*) "Keyword '", key,"' unknown!"
!     write(*,*) "exiting..."
!     stop
!    end if ! end over possible cases key
! 
!  RETURN
!  END SUBROUTINE initCollision

 SUBROUTINE CollisionRate(icell, atom)
 ! --------------------------------------------------
 ! Computes collisional rates and fills the C matrix.
 ! The routine is made such that new recipe can be
 ! implemented easily (I think).
 ! C_ij = C(i,j) = C j->i
 ! Remeber, atom%C(Nlevel*Nlevel) now
 ! we define indexes ij, and ji to run in the matrix
 !
 !
 ! (1) atomic file is opened with openCollisionFile at the
 ! atom%offset_coll position in the file. Before leaving
 ! the subroutine, the cursor is reset to atom%offset_coll position
 ! for next cell point. 
 ! The file is closed at the end of the simulation. 
 ! This trick avoid closing/reopening the file at each cell point.
 ! --------------------------------------------------
  type (AtomType), intent(inout) :: atom
  integer, intent(in) :: icell
  integer :: ij, k, Nread, countline=0, colunit, checkfseek, fseek
  integer :: NTMP, n, m, ii, Nitem, i1, i2, i, j, ji, Nitem2
  character(len=8) :: END_OF_FILE="END     ", key
  real(8), dimension(:), allocatable :: TGRID, coeff
  double precision :: C, np
  real(8) :: deltaE, C0, Cdown, Cup, gij, xj, fac, fxj
  integer :: Ncoef, Nrow
  character(len=MAX_LENGTH) :: inputline, FormatLine
  real(8), dimension(:,:), allocatable :: badi, cdi
  real(8) :: acolsh, tcolsh, aradsh, xradsh, adish, bdish
  real(8) :: t0sh, t1sh, summrs, tg, cdn, ccup
  real(8) :: ar85t1, ar85t2, ar85a, ar85b, ar85c, ar85d, t4
  real(8) :: de,zz,betab,cbar,dekt,dekti,wlog,wb,sumscl=0.0

  ! -- Nitem2 is here to check
  ! that Nitem (read) is equals to the acutla number
  ! of item read Nitem2
  ! Actually the error will be in read(inputline,*)
  ! because it reads in an array of size Nitem, so
  ! Nitem are expected to be read and an error will be
  ! raise otherwise.

  ! file is already opened, and will be close at the end of the RT run
  colunit = atom%colunit
  !write(*,*) "Collision()->colunit = ", colunit, atom%colunit
 
  write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
  
  C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) *&
        sqrt(8.0/(PI*KBOLTZMANN))

  ! -- Initialise the collisional matrix at each cell--
  atom%C(:) = 0d0
  !temporary
  atom%Ckij(icell,:) = 0d0 !C j->i = C(j,i)



  key="        "
  ! -- Go throught the remaining lines in the file --
  ! read collisional data depending the cases of recipes.

  do while(key.ne.END_OF_FILE)
   countline = countline + 1
   CALL getnextline(colunit, COMMENT_CHAR, &
         FormatLine, inputline, Nread)

   ! do not go further, because after there is
   ! the number of grid points, which is in general
   ! one or two digits.
   ! if you intend to increase the number of grid points
   ! for the interpolation of the coefficients
   ! you have to modify the file to flush right
   ! everything after the TEMP keyword.

   key = adjustl(inputline(1:len(key)))
   !write(*,*) trim(inputline)
   !write(*,*) "Line = ", countline, " key=",key,"."
   if (key.eq."TEMP") then
    read(inputline(len(key)+1:len(key)+3), *) NTMP
    !write(*,*) "NTMP = ", NTMP
    if (allocated(TGRID)) deallocate(TGRID)
    allocate(TGRID(NTMP))
    ! Nread is len(NTMP) in string format, offset
    ! inputline to read TGRID only, after NTMP.

    read(inputline(len(key)+3:Nread),*) (TGRID(k), k=1,NTMP)
    !write(*,*) (TGRID(k), k=1,NTMP)
    Nitem = NTMP
    Nitem2 = size(TGRID)

   else if ((key.eq."OMEGA") .or. (key.eq."CE") .or.&
        (key.eq."CI") .or. (key.eq."CP") .or. &
        (key.eq."CH0") .or. (key.eq."CH+") .or. &
        (key.eq."CH") .or. (key.eq."CR")) then

    ! -- read level indices and collision coefficients --
    Nitem = NTMP
    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(Nitem))
    ! offset inputline to avoid re-reading key.
    ! but read(inputline,*) key, i1, i2, coeff is OK
    !write(*,*) key, inputline(len(key)+1:)
    read(inputline(len(key)+1:),*) i1, i2, &
         (coeff(k),k=1,Nitem)

    Nitem2 = size(coeff)

    i = MIN(i1,i2) + 1
    j = MAX(i1,i2) + 1 !from 1 to Nlevel (in C, 0 to Nlevel-1)
    ij = (i-1)*atom%Nlevel + j!j->i
    ji = (j-1)*atom%Nlevel + i !i->j
    ! -- Note:
    ! Cf = C(Nlevel,Nlevel).flatten, of size Nlevel*Nlevel
    ! C(i,j) = Cf(ij) = Cf(i*Nlevel + j)
    ! C(1,1) = Cf(1) -> i*Nlevel +j = 1 -> i=i-1

   else if ((key.eq."AR85-CHP").or.(key.eq."AR85-CHH")) then
    Nitem = 6
    if (allocated(coeff)) deallocate(coeff)
    allocate(coeff(Nitem))
    read(inputline(len(key)+1:),*) i1, i2, &
        (coeff(k),k=1,Nitem)
    !write(*,*) inputline(len(key)+1:)

    Nitem2 = size(coeff)

    i = MIN(i1,i2) + 1
    j = MAX(i1,i2) + 1
    ij = (i-1)*atom%Nlevel + j
    ji = (j-1)*atom%Nlevel + i

   else if ((key.eq.'AR85-CEA').or.(key.eq."BURGESS")) then
     Nitem = 1
     if (allocated(coeff)) deallocate(coeff)
     allocate(coeff(Nitem))

     read(inputline(len(key)+1:),*) i1, i2, coeff(1)
     i = MIN(i1,i2) + 1
     j = MAX(i1,i2) + 1
     ij = (i-1)*atom%Nlevel + j
     ji = (j-1)*atom%Nlevel + i

     Nitem2 = size(coeff)
     !write(*,*) "CEA or BURGESS", i1, i2, i, ij

   else if (key.eq."SHULL82") then
     Nitem = 8
     if (allocated(coeff)) deallocate(coeff)
     allocate(coeff(Nitem))

     read(inputline(len(key)+1:),*) i1, i2,&
         (coeff(k),k=1,Nitem)
     i = MIN(i1,i2) + 1
     j = MAX(i1,i2) + 1
     ij = (i-1)*atom%Nlevel+j
     ji = (j-1)*atom%Nlevel +i

     Nitem2 = size(coeff)


   else if (key.eq."BADNELL") then
    ! -- BADNELL formula for dielectronic recombination --

    !write(*,*) inputline, Nread
    read(inputline(len(key)+1:),*) i1, i2, Ncoef
    Nrow = 2
    Nitem = Nrow*Ncoef
    allocate(badi(Nrow,Ncoef))
    ! they are two lines of Nitem elements to read
    !write(*,*) i1, i2, Ncoef
    do m=1,Nrow
      CALL getnextline(colunit, COMMENT_CHAR, &
             FormatLine, inputline, Nread)
      !write(*,*) 'row=',m, inputline, Nread
      read(inputline,*) (badi(m,k),k=1,Ncoef)
    end do
    i = MIN(i1,i2) + 1
    j = MAX(i1,i2) + 1
    ij = atom%Nlevel*(i-1) + j
    ji = atom%Nlevel*(j-1) + i

    Nitem2 = size(badi(1,:)) * size(badi(:,1))

   else if (key.eq."SUMMERS") then
    ! -- Switch for density dependent DR coefficient
    !
    ! Give the default multiplication factor of summers
    ! density dependence of dielectronic recombination:
    ! sumscl = 0. -> no density dependence
    ! sumscl = 1. -> full density dependence

    Nitem = 1
    read(inputline(len(key)+1:), *) sumscl

    Nitem2 = 1

   else if (key.eq."AR85-CDI") then
    read(inputline(len(key)+1:),*) i1, i2, Nrow

    if (Nrow.gt.MSHELL) then
     write(*,*) "Nrow is greater than MSHELL, exiting..."
     stop
    end if
    Nitem = Nrow*MSHELL
    allocate(cdi(Nrow, MSHELL))
    do m=1,Nrow
      CALL getnextline(colunit, COMMENT_CHAR, &
       FormatLine, inputline, Nread)
      read(inputline,*) (cdi(m,k),k=1,MSHELL)
    end do

    Nitem2 = size(cdi(1,:)) * size(cdi(:,1))

    i = MIN(i1,i2) + 1
    j = MAX(i1,i2) + 1
    ij = (i-1)*atom%Nlevel +j
    ji = (j-1)*atom%Nlevel +i

   else if (key.eq.END_OF_FILE) then
    exit
   else
    write(*,*) "Keyword '", key,"' unknown!"
    write(*,*) "exiting..."
    stop
   end if ! end over possible cases key

   if (Nitem2 .ne. Nitem) then
    ! -- Actually, the error should be in the reading.
    write(*,*) "Error, read ", Nitem2," items->expected ", &
     Nitem
    write(*,*) "This should never happen, here."
    write(*,*) "Exiting..."
    stop
   end if

   ! -- END of reading, filling now ...

   if ((key.eq."OMEGA") .or. (key.eq."CE") .or. &
       (key.eq."CI") .or. (key.eq."CP") .or. &
       (key.eq."CH0") .or. (key.eq."CH+") .or. &
       (key.eq."CH") .or. (key.eq."CR") ) then


                      !here Nitem is NTMP !!
    C = interp1D(TGRID, coeff, atmos%T(icell))
    !-> because we compute cell by cell now
    !CALL bezier3_interp(Nitem, TGRID, coeff, &
    !     atmos%Nspace, atmos%T, C)
    !CALL bezier2_interp(Nitem, TGRID, coeff, &
    !     atmos%Nspace, atmos%T, C)

    !write(*,*) "Coefficients interpolated onto the T grid"
    !write(*,*) (C(k),k=1,atmos%Nspace)
   end if
   if (key.eq."OMEGA") then
    ! -- Collisional excitation of ions

    !!do k=1, atmos%Nspace !! cell by cell
     !! C is a constant, replace by C(k) if for all grid at once
     !! and atom%C(ij) by atom%C(ij,k)
     Cdown = C0 * atmos%ne(icell) * C / &
      (atom%g(j)*dsqrt(atmos%T(icell)))
     atom%C(ij) = atom%C(ij) + Cdown
     ! remember the relation between Cij and Cji
     ! which takes place at LTE.
     atom%C(ji) = atom%C(ji) + Cdown * &
      atom%nstar(j,icell)/atom%nstar(i,icell)
    !!end do !! cell by cell
   else if (key.eq."CE") then
    ! -- Collisional excitation of neutrals
    gij = atom%g(i)/atom%g(j)
    !!do k=1,atmos%Nspace
     Cdown = C * atmos%ne(icell) * gij * dsqrt(atmos%T(icell))
     !write(*,*) key, "k=",k, "Cdown = ", Cdown, C(k)
     !write(*,*) "ne=",atmos%ne(k), gij, "sqrt(T)=",dsqrt(atmos%T(k))
     atom%C(ij) = atom%C(ij) + Cdown
     atom%C(ji) = atom%C(ji) + Cdown * &
      atom%nstar(j,icell)/atom%nstar(i,icell)
    !!end do
   else if (key.eq."CI") then
    ! -- Collisional ionisation
    deltaE = atom%E(j) - atom%E(i)
    !!do k=1,atmos%Nspace
     Cup = C * atmos%ne(icell) * &
       exp(-deltaE/(KBOLTZMANN*atmos%T(icell))) *dsqrt(atmos%T(icell))
     !write(*,*) key, "k=",k, "Cup = ", Cup, C(k)
     !write(*,*) "dE=",deltaE," exp()=",exp(-deltaE/(KBOLTZMANN*atmos%T(k)))
     atom%C(ji) = atom%C(ji) + Cup
     atom%C(ij) = atom%C(ij) + Cup * &
      atom%nstar(i,icell)/atom%nstar(j,icell)
    !!end do
   else if (key.eq."CR") then
    ! -- Collisional de-excitation by electrons
    !!do k=1,atmos%Nspace
     Cdown = atmos%ne(icell) * C
     atom%C(ij) = atom%C(ij) + Cdown
    !!end do
   else if (key.eq."CP") then
    ! -- collisions with protons
    ! protons are the last level of Hydrogen atoms
    np = Hydrogen%n(Hydrogen%Nlevel,icell)
    !!do k=1,atmos%Nspace
     Cdown = np * C
     atom%C(ij) = atom%C(ij) + Cdown
     atom%C(ji) = atom%C(ji) + Cdown * &
     atom%nstar(j,icell) / atom%nstar(i,icell)
    !!end do
   else if (key.eq."CH") then
    ! -- Collisions with neutral hydrogen
    !!do k=1,atmos%Nspace
     Cup = Hydrogen%n(1,icell) * C
     atom%C(ji) = atom%C(ji) + Cup
     atom%C(ij) = atom%C(ij) + Cup * &
      atom%nstar(i,icell) / atom%nstar(j,icell)
   !! end do
   else if (key.eq."CH0") then
    ! -- Charge exchange with neutral hydrogen
    !!do k=1,atmos%Nspace
     atom%C(ij) = atom%C(ij) + Hydrogen%n(1,icell)*C
    !!end do
   else if (key.eq."CH+") then
    ! -- charge exchange with protons
    np = Hydrogen%n(Hydrogen%Nlevel,icell)
    !!do k=1,atmos%Nspace
     atom%C(ji) = atom%C(ji) + np*C
    !!end do
   else if (key.eq."SHULL82") then
    acolsh = coeff(1)
    tcolsh = coeff(2)
    aradsh = coeff(3)
    xradsh = coeff(4)
    adish = coeff(5)
    bdish = coeff(6)
    t0sh = coeff(7)
    t1sh = coeff(8)
    !!do k=1,atmos%Nspace
     summrs = sumscl*summers(i,j,atmos%ne(icell),atom)
     summrs = summrs + (1.-sumscl)
     tg = atmos%T(icell)
     cdn = aradsh * dpow(tg/1.d4, -xradsh) + &
       summrs*adish / tg / dsqrt(tg) * exp(-t0sh/tg) * &
       (1. + (bdish*exp(-t1sh/tg)))
     cup = acolsh * dsqrt(tg) * exp(-tcolsh/tg) / &
      (1. + 0.1*tg/tcolsh)
     ! -- convert from cm3 /s to m3/s
     cdn = cdn * atmos%ne(icell) * CUBE(CM_TO_M)
     cup = cup*atmos%ne(icell)*CUBE(CM_TO_M)
     ! -- 3-body recombination (high density limit)
     cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)
     !write(*,*) "k=",k, " cdn = ", cdn
     atom%C(ij) = atom%C(ij) + cdn
     atom%C(ji) = atom%C(ji) + cup
    !!end do
   else if (key.eq."BADNELL") then
    ! -- Fit for dielectronic recombination form Badnell
    ! First line coefficients are energies in K (from Chianti)
    ! Second line coefficients are coefficients (from Chianti)
    ! See Badnell 2006 for more details

    !!do k =1,atmos%Nspace
     summrs = sumscl*summers(i,j,atmos%ne(icell),atom) + &
       (1.-sumscl)
     tg = atmos%T(icell)
     cdn = 0.
     do ii=1,Ncoef
      cdn = cdn + badi(2,ii) * exp(-badi(1,ii)/tg)
     end do
     cdn = cdn * dpow(tg, -1.5d0)
     !write(*,*) "k=",k, " cdn = ", cdn, " summrs = ",summrs, "cdn=", cdn
     ! -- convert from cm3/s to m3/s
     cdn = cdn *atmos%ne(icell) * summrs * CUBE(CM_TO_M)
     cup = cdn * atom%nstar(j,icell)/atom%nstar(i,icell)

     cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)

     atom%C(ij) = atom%C(ij) + cdn
     atom%C(ji) = atom%C(ji) + cup
     !write(*,*) "BADNELL", cdn, cup
    !!end do
    deallocate(badi)
   else if (key.eq."AR85-CDI") then
    ! -- Direct collisional ionisation
    !!do k=1,atmos%Nspace
     cup = 0.
     tg = atmos%T(icell)
     do m=1,Nrow
      xj = cdi(m,1) * EV/ (KBOLTZMANN*tg)
      fac = exp(-xj) * dsqrt(xj)
      fxj = cdi(m,2) + cdi(m,3) * (1.+xj) + &
       (cdi(m,4) - xj*(cdi(m,2)+cdi(m,3)*(2.+xj)))*&
        fone(xj) + cdi(m,5)*xj*ftwo(xj)
      fxj = fxj * fac
      !write(*,*) fxj, ftwo(xj), fone(xj)
      fac = 6.69d-7 / dpow(cdi(m,1),1.5d0)
      cup = cup + fac*fxj*CUBE(CM_TO_M)
     end do
     if (cup.lt.0.) cup = 0.
     cup = cup * atmos%ne(icell)
     cdn = cup * atom%nstar(i,icell)/atom%nstar(j,icell)
     atom%C(ij) = atom%C(ij) + cdn
     atom%C(ji) = atom%C(ji) + cup
     !write(*,*) "CDU", cdn, cup
     !write(*,*) "CDI: line=",countline,ij, k, " C[ij,k]=",atom%C(ij,k), &
     ! " C[ji,k]=",atom%C(ji,k)
    !!end do
    deallocate(cdi)
   else if (key.eq."AR85-CEA") then
    ! -- Autoionisation
    !!do k=1,atmos%Nspace
     fac = ar85cea(i,j,icell,atom)
     !write(*,*) "fac=", fac
     cup = coeff(1)*fac*atmos%ne(icell)
     atom%C(ji) = atom%C(ji) + cup
     !write(*,*) "AR85-CEA, cup=", cup
     !write(*,*) "CEA: line=",countline,ij, k, " C[ij,k]=",atom%C(ij,k), &
     !  " C[ji,k]=",atom%C(ji,k)
    !!end do
   else if (key.eq."AR85-CHP") then
    ! -- Charge transfer with ionised Hydrogen
    ar85t1 = coeff(1)
    ar85t2 = coeff(2)
    ar85a = coeff(3)
    ar85b = coeff(4)
    ar85c = coeff(5)
    ar85d = coeff(6)
    !!do k =1,atmos%Nspace
     if ((atmos%T(icell).ge.ar85t1).and.&
          (atmos%T(icell).le.ar85t2)) then
      t4 = atmos%T(icell)/1d4
      cup = ar85a * 1d-9 * dpow(t4,ar85b) * &
        exp(-ar85c*t4) * &
        exp(-ar85d*EV/KBOLTZMANN/atmos%T(icell)) * &
        Hydrogen%n(Hydrogen%Nlevel,icell) * &
        CUBE(CM_TO_M)
      atom%C(ji) = atom%C(ji) + cup
     end if
    !!end do
   else if (key.eq."AR85-CHH") then
    ! Charge transfer with neutral Hydrogen
    ar85t1 = coeff(1)
    ar85t2 = coeff(2)
    ar85a = coeff(3)
    ar85b = coeff(4)
    ar85c = coeff(5)
    ar85d = coeff(6)
    !!do k=1,atmos%Nspace
     if ((atmos%T(icell).ge.ar85t1).and.&
         (atmos%T(icell).le.ar85t2)) then
      t4 = atmos%T(icell)/1d4
      cdn = ar85a * 1d-9 * dpow(t4, ar85b) * &
       (1. + ar85c*exp(ar85d*t4)) * &
       Hydrogen%n(1,icell) * CUBE(CM_TO_M)
      atom%C(ij) = atom%C(ij) + cdn
     end if
    !!end do
   else if (key.eq."BURGESS") then
    write(*,*) "BURGESS NOT CHECK"
    ! -- Electron impact ionisation from Burgess Chidichimo
    ! 1983, MNRAS, 203, 1269-1280
    de = (atom%E(j)-atom%E(i))/EV
    zz = atom%stage(i) !0 for neutrals
    betab = 0.25*(dsqrt((100.*zz+91)/(4.*zz+3.))-5.)
    cbar = 2.3
    !!do k=1,atmos%Nspace
     dekt = de*EV / (KBOLTZMANN*atmos%T(k))
     dekt = MIN(500., dekt)
     dekti = 1./dekt
     wlog = log(1.+dekti)
     wb = dpow(wlog, betab/(1.+dekti))
     cup = 2.1715d-8 * cbar * dpow(13.6/de, 1.5d0) * &
         dsqrt(dekt) * E1(dekt) * wb * atmos%ne(icell) * &
         CUBE(CM_TO_M)
     ! -- Fudge factor
     cup = cup * coeff(1)
     cdn = cup * atom%nstar(i,icell) / atom%nstar(j,icell)
     write(*,*) "BRUGESS, cdn=", cdn, " cup=", cup
     atom%C(ji) = atom%C(ji) + cup
     atom%C(ij) = atom%C(ij) + cdn
    !!end do
   end if
  end do !loop over (remaining) lines (of the atomic model)

 
 !!close(colunit) !! closed latter
 !! reset the cursor for next cell point !!
 checkfseek = fseek(atom%colunit,atom%offset_coll, 0)
 if (checkfseek.gt. 0 ) then 
   write(*,*) 'fseek error'
   stop 
 end if

 atom%Ckij(icell,:) = atom%C(:)
 deallocate(TGRID)
 !!deallocate(C) !! not an array anymore
 deallocate(coeff)

 RETURN
 END SUBROUTINE CollisionRate

END MODULE collision
