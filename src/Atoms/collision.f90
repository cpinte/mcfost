! ---------------------------------------------------------------------- !
! Convention: C_ij = C[i][j] represents the
! transition j --> i = Cul
! C_ij = C[ith ligne][jth column]
! fortran is column row
!     ij = (i-1)*atom%Nlevel +j : j-->i
!     ji = (j-1)*atom%Nlevel +i : i-->j
! ---------------------------------------------------------------------- !

module collision
  use math, only : E1, E2
  use constant
  use atom_type, only : AtomType, ATOM_ID_WIDTH, PTrowcol, atomZnumber, Nmax_line_per_collision
  use atmos_type, only : Hydrogen, Helium, ne, Atoms, T, Elements
  use messages
  use mcfost_env, only : dp
  use utils, only : interp_dp

  implicit none

  real, parameter :: BRK=4.0
  integer, parameter :: MSHELL=5
  character, parameter :: COMMENT_CHAR="#"
  real(kind=dp), parameter :: C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*(RBOHR)**2) * sqrt(8.0/(PI*KBOLTZMANN))
  integer, parameter :: Max_size_elements = 100
  integer, parameter :: Nmax_lines = 10000 !number max of lines for all collision data

contains

  subroutine read_collisions(unit, atom)

    integer, intent(in) :: unit
    type (AtomType), intent(inout) :: atom
    character(len=Nmax_line_per_collision), dimension(Nmax_lines) :: lines_in_file !check len char matches the one in atom%
    integer :: N !real number of lines
    integer :: Nread, status
    character(len=3) :: END_OF_FILE="END", key
    character(len=Nmax_line_per_collision*10) :: inputline

    n = 0
    !it is still important to read an END in the file
    key = "   "
    do while(key /= END_OF_FILE)
      read(unit,'(1A<Nmax_line_per_collision*10>)',iostat=status) inputline
      key = adjustl(inputline)
      Nread = len(trim(inputline))
      if (inputline(1:1)=='#' .or. Nread==0) cycle
      if (Nread > Nmax_line_per_collision) call error("Collision: Nread > Nmax_line_per_collision")
      n = n + 1
      if (n > Nmax_lines) call error("Collision: not enough number of lines for all collision data")
      lines_in_file(n) = inputline(1:Nread)
    enddo

    allocate(atom%collision_lines(N))
    atom%collision_lines(:) = lines_in_file(1:N)
   !  do nread=1,N
   !    write(*,*) atom%collision_lines(nread)
   !  enddo
   !  if (atom%ID=="He")stop

    return
  end subroutine read_collisions

  function fone(x) result(y)
    ! --------------------------------------------
    ! f_1 function of
    !  Arnaud & Rothenflug, 1985, A&ASS, 60, 425
    ! --------------------------------------------

    real(kind=dp) :: x, y

    if (x <= 50.0) then
       y = exp(x)*E1(x)
    else
       y = 1./x
    end if

    return

  end function fone

  function ftwo(x) result(y)
    ! --------------------------------------------
    ! f_2 function of
    !  Arnaud & Rothenflug, 1985, A&ASS, 60, 425
    !
    !  Improved description when x < 4 from:
    !    Hummer, 1983, jqsrt, 30 281
    ! --------------------------------------------

    integer :: i
    real(kind=dp) :: x, y, p(15), q(15), px, xfact
    real(kind=dp) :: qx, gamma, f0x, count, fact, term

    p = (/1.0000d+00, 2.1658d+02, 2.0336d+04, 1.0911d+06,&
         3.7114d+07, 8.3963d+08, 1.2889d+10, 1.3449d+11,&
         9.4002d+11, 4.2571d+12, 1.1743d+13, 1.7549d+13,&
         1.0806d+13, 4.9776d+11, 0d0/)

    q = (/1.0000d+00, 2.1958d+02, 2.0984d+04, 1.1517d+06,&
         4.0349d+07, 9.4900d+08, 1.5345d+10, 1.7182d+11,&
         1.3249d+12, 6.9071d+12, 2.3531d+13, 4.9432d+13,&
         5.7760d+13, 3.0225d+13, 3.3641d+12/)

    if (x > BRK) then
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
       y = px/(qx*x*x)
    else
       gamma = 0.5772156649
       f0x = pi*pi / 12.
       term = 1.
       count = 0.0
       fact = 1.0
       xfact = 1.0

       do while (abs(term/f0x) > 1d-8)
          count = count+1
          fact = fact*count
          xfact = xfact*(-x)
          term = xfact/(count*count*fact)
          f0x = f0x+term
          if (count > 100.0) then
             write(*,*) "Error in ftwo function in collision.f90"
             write(*,*) "code not stopped..."
          end if
       end do
       y = exp(x) * ((log(x) + gamma) * (log(x) + gamma)*0.5 + f0x)
    end if
    return
  end function ftwo

  ! -----------------------------------------------
  ! REMEMBER -> stage is like in C, it goes from 0
  ! for neutrals, to Nstage-1 for the most ionised
  ! stage (instead of 1 to Nstage)
  ! Furhter, in readatom.f90, i is i+1 and j is j+1
  ! to be in aggreement with the index convention of
  ! fortran 90. First level is i(j) = 1, not 0 !
  ! and last level is i(j) = Nlevel, not Nlevel-1.
  ! -----------------------------------------------

  function ar85cea(i,j,k, atom) result(cup)
    ! -----------------------------------------------
    ! Computes collisional autoionisation rates
    ! using the formalism from Arnaud and Rothenflug
    ! 1985, A&ASS, 60, 425 (a.k.a ar85)
    !
    ! -----------------------------------------------
    integer :: i, j, k, iz, ichrge, isoseq
    type (AtomType) :: atom
    real(kind=dp) :: cup
    character(len=ATOM_ID_WIDTH) :: cseq
    real(kind=dp) :: zz, bkt, b, zeff, iea, y
    real(kind=dp) :: f1y,a,g,cea

    ! -- Initialisation --
    cea = 0.
    y = 0.
    f1y = 0.
    cup = 0.

    ! -- Search for element --
    iz = atomZnumber(atom)!Hydrogen 1 !old atomID
    zz = iz
    if (iz < 1 .or. iz > .92) then
       write(*,*) k, "(ar85cea): Limits of the routine in terms of elements exceeded, returing..."
       return
    end if

    ! -- get iso-electronic sequence --
    ichrge = atom%stage(i) !stages from 0 for neutrals to N-1
    isoseq = iz-ichrge
    if (isoseq < 29) cseq = Elements(isoseq)%ptr_elem%ID

    ! -- Temperature in eV --
    bkt = KBOLTZMANN * T(k) / EV

    ! All ID in Elements are of the form Xx
    ! just as in atom%ID
    ! -- Lithium sequence --
    if (cseq.eq."Li") then
       iea = 13.6*((zz-0.835d0)**(2.d0)-0.25*(zz-1.62d0)**(2.d0))
       b = 1./(1.+2.0d-4 * (zz)**(3.d0))
       zeff = zz-0.43
       y = iea/bkt
       f1y = fone(y)
       g = 2.22*f1y + 0.67*(1.-y*f1y)+0.49*y*f1y + 1.2*y*(1.-y*f1y)
       cup = (1.6d-7 * 1.2*b) / ((zeff)**(2.d0)*sqrt(bkt)) * exp(-y)*g
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
       if (iz <= 16) then
          iea = 26.*(zz-10.)
          a = 2.9d-17 * (zz-11.)**(-7d-1)
          y = iea/bkt
          f1y = fone(y)
          cup = 6.69d7 * a * iea / sqrt(bkt) * exp(-y) * (1.-y*f1y)
       else if (iz >= 18 .and. iz <= 28) then
          iea = 11.*(zz-10.)*sqrt(zz-10.)
          a = 1.4e-14 * (zz-10.)**(-3.73d0)
          y = iea/bkt
          f1y = fone(y)
          cup = 6.69d7 * a * iea/sqrt(bkt)  * exp(-y) * (1.-0.5*(y-y*y+y*y*y*f1y))
       else
          cup = 0.
       end if
    end if
    ! -- Magnesium-sulfur sequences --
    if ((cseq.eq."Mg") .or. (cseq.eq."Al") .or. &
         (cseq.eq."Si") .or. (cseq.eq."P ") .or. &
         (cseq.eq."S ")) then

       if (cseq.eq."Mg") iea=10.3*(zz-10.)**(1.52d0)
       if (cseq.eq."Al") iea=18.0*(zz-11.)**(1.33d0)
       if (cseq.eq."Si") iea=18.4*(zz-12.)**(1.36d0)
       if (cseq.eq."P ") iea=23.7*(zz-13.)**(1.29d0)
       if (cseq.eq."S ") iea=40.1*(zz-14.)**(1.1d0)

       a = 4.d-13 / (zz*zz*iea)
       y = iea/bkt
       f1y = fone(y)
       cup = 6.69d7 *a*iea / sqrt(bkt) * exp(-y) * (1.-0.5*(y-y*y+y*y*y*f1y))

    end if
    ! -- Special cases --
    if ((atom%ID.eq."Ca").and.(ichrge==0)) then
       iea = 25.
       a = 9.8d-17
       b = 1.12
       cup = 6.69d7 * a * iea / sqrt(bkt) * exp(-y) * (1. + b*f1y)
    else if ((atom%ID.eq."Ca").and.(ichrge==1)) then
       iea = 25.0
       a = 6.0d-17
       b = 1.12
       cup = 6.69d7 * a * iea / sqrt(bkt) * exp(-y) * (1. + b*f1y)
    else if ((atom%ID.eq."Fe").and.(ichrge==3)) then
       a = 1.8d-17
       iea = 60.0
       b = 1.0
       cup = 6.69d7 * a * iea / sqrt(bkt) * exp(-y) * (1. + b*f1y)
    else if ((atom%ID.eq."Fe").and.(ichrge==4)) then
       a = 5.d-17
       iea = 73.
       b = 1.
       cup = 6.69d7 * a * iea / sqrt(bkt) * exp(-y) * (1.+b*f1y)
    end if

    cup = cup * (CM_TO_M)**3 ! in m3

    return
  end function ar85cea

  function summers(i, j, nne, atom) result(y)
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
    real(kind=dp) :: y, zz, rho0, rhoq, x, beta, nne
    type (AtomType) :: atom
    character(len=ATOM_ID_WIDTH) :: cseq

    iz = atomZnumber(atom)
    if (iz < 1 .or. iz > 92) then
       write(*,*) "(summers): Limits of the routine in terms of elements exceeded, returing..."
    end if

    ! -- charge of recombining ion --
    zz = atom%stage(j)
    isoseq = iz - atom%stage(i)
    if (isoseq < 29) cseq = Elements(isoseq)%ptr_elem%ID


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


    call PTrowcol(isoseq, row, col)
    if (isoseq==2) col = col-16 ! Helium, relative col=2
    if ( ((isoseq>=5).and.(isoseq<=10)).or.&
         ((isoseq>=13).and.(isoseq<=18))) col=col-10
    !write(*,*) "row=",row, "col=",col

    rhoq = nne*(CM_TO_M)**3 / (zz)**(7.d0)
    x = (0.5*zz + (col-1.)) * row/3.
    beta = -0.2 / log(x+2.71828)
    rho0 = 30.0 + 50.0*x;
    y = (1. + rhoq/rho0)**(beta)

    return
  end function summers

  function collision_rates_atom_loc(icell, atom, deriv) result(C)
    ! --------------------------------------------------------------------------- !
    ! Computes collisional rates and fills the C matrix.
    ! The routine is made such that new recipe can be
    ! implemented easily (I think).
    !
    ! C(i,j) = Cji = C i -> j
    !
    ! C in m^-3 s^-1 collisional rate units.
    ! it is ne(icell) * Cij
    !
    ! Collision data are stored in an array as characters.
    ! The function then reads this array at each cell to compute
    ! the collisional rates. This is not very efficient as interpolation
    ! are done on the fly. Still, it avoids storing the full collision matrix
    ! in memory and the time to compute the collisional rates is negligible
    ! compared to the time to integrate the RTE.
    ! --------------------------------------------------------------------------- !
    type (AtomType), intent(inout) :: atom
    integer, intent(in) :: icell
    real(kind=dp), dimension(atom%Nlevel,atom%Nlevel) :: C
    real(kind=dp), intent(out), optional :: deriv(atom%Nlevel, atom%Nlevel)
    integer :: ij, k, Nread, countline=0, colunit
    integer :: NTMP, n, m, ii, Nitem, i1, i2, i, j, ji, Nitem2
    character(len=3) :: END_OF_FILE="END"
    integer, parameter :: max_len_key = 8!set by the maximum length of a keyword!!
    character(len=max_len_key) :: key
    ! 		real(kind=dp), dimension(:), allocatable :: TGRID, coeff
    ! 		real(kind=dp), dimension(:,:), allocatable :: badi, cdi
    real(kind=dp), dimension(Max_size_elements) :: TGRID, coeff
    real(kind=dp), dimension(Max_size_elements,Max_size_elements) :: badi, cdi
    real(kind=dp) :: np, CC
    real(kind=dp) :: deltaE, Cdown, Cup, gij, xj, fac, fxj
    integer :: Ncoef, Nrow, Nlines, k1
    character(len=Nmax_line_per_collision) :: inputline, FormatLine
    real(kind=dp) :: acolsh, tcolsh, aradsh, xradsh, adish, bdish
    real(kind=dp) :: t0sh, t1sh, summrs, tg, cdn, ccup
    real(kind=dp) :: ar85t1, ar85t2, ar85a, ar85b, ar85c, ar85d, t4
    real(kind=dp) :: de,zz,betab,cbar,dekt,dekti,wlog,wb,sumscl
    real(kind=dp) :: ni_on_nj

    sumscl = 0.0

    ! -- Nitem2 is here to check
    ! that Nitem (read) is equals to the acutla number
    ! of item read Nitem2
    ! Actually the error will be in read(inputline,*)
    ! because it reads in an array of size Nitem, so
    ! Nitem are expected to be read and an error will be
    ! raise otherwise.

    write(FormatLine,'("(1"A,I5")")') "A", Nmax_line_per_collision


    ! -- Initialise the collisional matrix at each cell -- !
    C(:,:) = 0d0
    !derivative to ne
    if (present(deriv)) deriv(:,:) = 0.0_dp

    Nlines = size(atom%collision_lines)

    ! -- Go throught the remaining lines in the file -- !
    ! read collisional data depending the cases of recipes.

    loop_lines_in_file : do k1=1, Nlines
       countline = countline + 1
       inputline = atom%collision_lines(k1)
       Nread = len(inputline)!already trimed
       ! do not go further, because after there is
       ! the number of grid points, which is in general
       ! one or two digits.
       ! if you intend to increase the number of grid points
       ! for the interpolation of the coefficients
       ! you have to modify the file to flush right
       ! everything after the TEMP keyword.

       key = adjustl(inputline(1:max_len_key))
       !write(*,*) trim(inputline)
       !write(*,*) "Line = ", countline, " key=",key,"."
       if (key.eq."TEMP") then
          read(inputline(max_len_key+1:max_len_key+3), *) NTMP
          !write(*,*) "NTMP = ", NTMP
          ! 				if (allocated(TGRID)) deallocate(TGRID)
          ! 				allocate(TGRID(NTMP))
          ! Nread is len(NTMP) in string format, offset
          ! inputline to read TGRID only, after NTMP.
          Tgrid(1:NTMP) = 0.0_dp
          read(inputline(max_len_key+3:Nread),*) (TGRID(k), k=1,NTMP)
          !write(*,*) (TGRID(k), k=1,NTMP)
          Nitem = NTMP
          if (Nitem > max_size_elements) call error("Number temperature points for collision data larger than the available size!")
          ! 				Nitem2 = size(TGRID) !does not work if not allocatable

       else if ((key.eq."OMEGA") .or. (key.eq."CE") .or.&
            (key.eq."CI") .or. (key.eq."CP") .or. &
            (key.eq."CH0") .or. (key.eq."CH+") .or. &
            (key.eq."CH") .or. (key.eq."CR")) then

          ! -- read level indices and collision coefficients --
          Nitem = NTMP
          if (Nitem > max_size_elements) call error("Number of coeffs for collision data larger than the available size!")
          ! 				if (allocated(coeff)) deallocate(coeff)
          ! 				allocate(coeff(Nitem))
          ! offset inputline to avoid re-reading key.
          ! but read(inputline,*) key, i1, i2, coeff is OK
          !write(*,*) key, inputline(len(key)+1:)
          coeff(1:Nitem) = 0.0_dp
          read(inputline(max_len_key+1:),*) i1, i2,(coeff(k),k=1,Nitem)
          ! 				Nitem2 = size(coeff)

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
          ! 				if (allocated(coeff)) deallocate(coeff)
          ! 				allocate(coeff(Nitem))
          coeff(1:Nitem) = 0.0_dp
          read(inputline(max_len_key+1:),*) i1, i2, (coeff(k),k=1,Nitem)
          !write(*,*) inputline(max_len_key+1:)
          ! 				Nitem2 = size(coeff)

          i = MIN(i1,i2) + 1
          j = MAX(i1,i2) + 1
          ij = (i-1)*atom%Nlevel + j
          ji = (j-1)*atom%Nlevel + i

       else if ((key.eq.'AR85-CEA').or.(key.eq."BURGESS")) then
          Nitem = 1
          ! 				if (allocated(coeff)) deallocate(coeff)
          ! 				allocate(coeff(Nitem))
          coeff(1:Nitem) = 0.0_dp
          read(inputline(max_len_key+1:),*) i1, i2, coeff(1)
          i = MIN(i1,i2) + 1
          j = MAX(i1,i2) + 1
          ij = (i-1)*atom%Nlevel + j
          ji = (j-1)*atom%Nlevel + i

          ! 				Nitem2 = size(coeff)
          !write(*,*) "CEA or BURGESS", i1, i2, i, ij

       else if (key.eq."SHULL82") then
          Nitem = 8
          ! 				if (allocated(coeff)) deallocate(coeff)
          ! 				allocate(coeff(Nitem))
          coeff(1:Nitem) = 0.0_dp
          read(inputline(max_len_key+1:),*) i1, i2,(coeff(k),k=1,Nitem)
          i = MIN(i1,i2) + 1
          j = MAX(i1,i2) + 1
          ij = (i-1)*atom%Nlevel + j
          ji = (j-1)*atom%Nlevel + i

          ! 				Nitem2 = size(coeff)


       else if (key.eq."BADNELL") then
          ! -- BADNELL formula for dielectronic recombination --

          !write(*,*) inputline, Nread
          read(inputline(max_len_key+1:),*) i1, i2, Ncoef
          Nrow = 2
          Nitem = Nrow*Ncoef
          if (Nitem > max_size_elements) call error("Number of coeffs for BADNELL collision data larger than the available size!")
          badi(1:Nrow, 1:Ncoef) = 0.0_dp
          ! 			allocate(badi(Nrow,Ncoef))
          ! they are two lines of Nitem elements to read
          !write(*,*) i1, i2, Ncoef
          do m=1,Nrow
             inputline = atom%collision_lines(k1+m)
             !write(*,*) 'row=',m, inputline, Nread
             read(inputline,*) (badi(m,k),k=1,Ncoef)
          end do
          i = MIN(i1,i2) + 1
          j = MAX(i1,i2) + 1
          ij = atom%Nlevel*(i-1) + j
          ji = atom%Nlevel*(j-1) + i

          ! 				Nitem2 = size(badi(1,:)) * size(badi(:,1))

       else if (key.eq."SUMMERS") then
          ! -- Switch for density dependent DR coefficient
          !
          ! Give the default multiplication factor of summers
          ! density dependence of dielectronic recombination:
          ! sumscl = 0. -> no density dependence
          ! sumscl = 1. -> full density dependence

          Nitem = 1
          read(inputline(max_len_key+1:), *) sumscl
          ! 				Nitem2 = 1

       else if (key.eq."AR85-CDI") then
          read(inputline(max_len_key+1:),*) i1, i2, Nrow

          if (Nrow > MSHELL) then
             call error("(Collision AR85-CDI) Nrow is greater than MSHELL, exiting...")
          end if
          Nitem = Nrow*MSHELL
          if (Nitem > max_size_elements) call error("Number of coeffs for AR85-CDI collision data larger than the available size!")
          cdi(1:Nrow, 1:MSHELL) = 0.0_dp
          ! 				allocate(cdi(Nrow, MSHELL))
          do m=1,Nrow
             inputline = atom%collision_lines(k1+m)
             read(inputline,*) (cdi(m,k),k=1,MSHELL)
          end do

          ! 				Nitem2 = size(cdi(1,:)) * size(cdi(:,1))

          i = MIN(i1,i2) + 1
          j = MAX(i1,i2) + 1
          ij = (i-1)*atom%Nlevel +j
          ji = (j-1)*atom%Nlevel +i

       else if (key.eq.END_OF_FILE) then
          exit loop_lines_in_file
       else
          write(*,*) "Error in collision! : Keyword ", key," unknown!"
          stop
       end if ! end over possible cases key

       !-> this error should be handled in reading
       ! 			if (Nitem2 .ne. Nitem) then
       ! 			! -- Actually, the error should be in the reading.
       ! 				write(*,*) "Error, read ", Nitem2," items->expected ", Nitem
       ! 				write(*,*) "This should never happen, here."
       ! 				write(*,*) "Exiting..."
       ! 				stop
       ! 			end if


       ! -- END of reading data for key, now computing collision rates -- !
       ! Still in the loop over the lines in the atomic file.

       if ((key.eq."OMEGA") .or. (key.eq."CE") .or. &
            (key.eq."CI") .or. (key.eq."CP") .or. &
            (key.eq."CH0") .or. (key.eq."CH+") .or. &
            (key.eq."CH").or.(key.eq."CR")) then


          !here Nitem is NTMP !
          if (Nitem /= NTMP) call error("Error, Nitem should be NTMP!")

          !             	CC = interp_dp(coeff, TGRID, T(icell))
          !if coeff and TGRID are not allocatables,
          !the size of coeff and Tgrid must be given (otherwise interpolation errors appear) !
          CC = interp_dp(coeff(1:Nitem), TGRID(1:Nitem), T(icell))

       end if

       !Compute also the derivative if present.
       !note that in principle ni_bb / nj_cont is prop to ne and
       !ne * phi(T) should be use because LTE populations are updated afterwards
       !not during the electron + SEE iterations. See impact.f90 and see_ionisation_nonlte in statequil_atom.f90 for more details.

       if (key.eq."OMEGA") then
          ! -- Collisional excitation of ions
          !Cdown means from j->i
          Cdown = C0 * ne(icell) * CC / (atom%g(j)*sqrt(T(icell)))
          C(j,i) = C(j,i) + Cdown

          ! remember the relation between Cij and Cji
          ! which takes place at LTE.
          C(i,j) = C(i,j) + Cdown * atom%nstar(j,icell)/atom%nstar(i,icell)
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 * Cdown/ne(icell)
             deriv(i,j) = deriv(i,j) + Cdown * atom%nstar(j,icell)/atom%nstar(i,icell) / ne(icell)
          endif
       else if (key.eq."CE") then
          ! -- Collisional excitation of neutrals
          gij = atom%g(i)/atom%g(j)
          !write(*,*) CC * gij * sqrt(T(icell))
          !exp(de/kT)
          Cdown = CC * ne(icell) * gij * sqrt(T(icell)) !CE(RH)=CE*gj/gi*ni/nj / sqrt(T)=CC
          !write(*,*) key, "k=",k, "Cdown = ", Cdown, C(k)
          !write(*,*) "ne=",ne(k), gij, "sqrt(T)=",sqrt(T(k))
          C(j,i) = C(j,i) + Cdown
          C(i,j) = C(i,j) + Cdown * atom%nstar(j,icell)/atom%nstar(i,icell)
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 * Cdown/ne(icell)
             deriv(i,j) = deriv(i,j) + Cdown * atom%nstar(j,icell)/atom%nstar(i,icell) / ne(icell)
          endif

       else if (key.eq."CI") then
          ! -- Collisional ionisation
          deltaE = atom%E(j) - atom%E(i)
          Cup = CC * ne(icell) * exp(-deltaE/(KBOLTZMANN*T(icell))) *sqrt(T(icell))
          !write(*,*) key, "k=",k, "Cup = ", Cup, C(k)
          !write(*,*) "dE=",deltaE," exp()=",exp(-deltaE/(KBOLTZMANN*T(k)))
          C(i,j) = C(i,j) + Cup
          C(j,i) = C(j,i) + Cup * atom%nstar(i,icell)/atom%nstar(j,icell)

          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * Cup * atom%nstar(i,icell)/atom%nstar(j,icell)
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif

       else if (key.eq."CR") then
          ! -- Collisional de-excitation by electrons
          Cdown = ne(icell) * CC
          C(j,i) = C(j,i) + Cdown !here

          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * cdown
          endif

       else if (key.eq."CP") then
          ! -- collisions with protons
          ! protons are the last level of Hydrogen atoms
          np = Hydrogen%n(Hydrogen%Nlevel,icell)
          Cdown = np * CC
          C(j,i) = C(j,i) + Cdown
          C(i,j) = C(i,j) + Cdown * atom%nstar(j,icell) / atom%nstar(i,icell)
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 * Cdown/ne(icell)
             deriv(i,j) = deriv(i,j) + Cdown * atom%nstar(j,icell)/atom%nstar(i,icell) / ne(icell)
          endif
       else if (key.eq."CH") then
          ! -- Collisions with neutral hydrogen
          Cup = Hydrogen%n(1,icell) * CC
          C(i,j) = C(i,j) + Cup
          C(j,i) = C(j,i) + Cup * atom%nstar(i,icell) / atom%nstar(j,icell)
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * Cup * atom%nstar(i,icell)/atom%nstar(j,icell)
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif
       else if (key.eq."CH0") then
          ! -- Charge exchange with neutral hydrogen
          C(j,i) = C(j,i) + Hydrogen%n(1,icell)*CC
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * Hydrogen%n(1,icell)*CC
          endif
       else if (key.eq."CH+") then
          ! -- charge exchange with protons
          np = Hydrogen%n(Hydrogen%Nlevel,icell)
          C(i,j) = C(i,j) + np*CC
          if (present(deriv)) then
             deriv(i,j) = deriv(i,j) + np*CC/ne(icell)
          endif
       else if (key.eq."SHULL82") then
          acolsh = coeff(1)
          tcolsh = coeff(2)
          aradsh = coeff(3)
          xradsh = coeff(4)
          adish = coeff(5)
          bdish = coeff(6)
          t0sh = coeff(7)
          t1sh = coeff(8)
          summrs = sumscl*summers(i,j,ne(icell),atom)
          summrs = summrs + (1.-sumscl)
          tg = T(icell)
          cdn = aradsh * (tg/1.d4)**(-xradsh) + summrs*adish / tg / sqrt(tg) * exp(-t0sh/tg) * (1. + (bdish*exp(-t1sh/tg)))
          cup = acolsh * sqrt(tg) * exp(-tcolsh/tg) / (1. + 0.1*tg/tcolsh)
          ! -- convert from cm3 /s to m3/s
          cdn = cdn * ne(icell) * (CM_TO_M)**3
          cup = cup*ne(icell) * (CM_TO_M)**3
          ! -- 3-body recombination (high density limit)
          cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)
          !write(*,*) "k=",k, " cdn = ", cdn
          C(j,i) = C(j,i) + cdn
          C(i,j) = C(i,j) + cup
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0*cdn / ne(icell)
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif
       else if (key.eq."BADNELL") then
          ! -- Fit for dielectronic recombination form Badnell
          ! First line coefficients are energies in K (from Chianti)
          ! Second line coefficients are coefficients (from Chianti)
          ! See Badnell 2006 for more details

          summrs = sumscl*summers(i,j,ne(icell),atom) + (1.-sumscl)
          tg = T(icell)
          cdn = 0.
          do ii=1,Ncoef
             cdn = cdn + badi(2,ii) * exp(-badi(1,ii)/tg)
          end do
          cdn = cdn * (tg)**(-1.5d0)
          !write(*,*) "k=",k, " cdn = ", cdn, " summrs = ",summrs, "cdn=", cdn
          ! -- convert from cm3/s to m3/s
          cdn = cdn *ne(icell) * summrs * (CM_TO_M)**3
          cup = cdn * atom%nstar(j,icell)/atom%nstar(i,icell)

          cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)

          C(j,i) = C(j,i) + cdn
          C(i,j) = C(i,j) + cup
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * cdn
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif

          ! 				deallocate(badi)
       else if (key.eq."AR85-CDI") then
          ! -- Direct collisional ionisation
          cup = 0.
          tg = T(icell)
          do m=1,Nrow
             xj = cdi(m,1) * EV/ (KBOLTZMANN*tg)
             fac = exp(-xj) * sqrt(xj)
             fxj = cdi(m,2) + cdi(m,3) * (1.+xj) + (cdi(m,4) - xj*(cdi(m,2)+cdi(m,3)*(2.+xj)))*&
                  fone(xj) + cdi(m,5)*xj*ftwo(xj)
             fxj = fxj * fac
             !write(*,*) fxj, ftwo(xj), fone(xj)
             fac = 6.69d-7 / (cdi(m,1))**(1.5d0)
             cup = cup + fac*fxj*(CM_TO_M)**3
          end do
          if (cup < 0.0_dp) cup = 0.0_dp
          cup = cup * ne(icell)
          cdn = cup * atom%nstar(i,icell)/atom%nstar(j,icell)
          C(j,i) = C(j,i) + cdn
          C(i,j) = C(i,j) + cup
          !write(*,*) "CDU", cdn, cup
          !write(*,*) "CDI: line=",countline,ij, k, " C[ij,k]=",C(ij,k), &
          ! " C[ji,k]=",C(ji,k)
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * cdn
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif
          ! 				deallocate(cdi)
       else if (key.eq."AR85-CEA") then
          ! -- Autoionisation
          fac = ar85cea(i,j,icell,atom)
          !write(*,*) "fac=", fac
          cup = coeff(1)*fac*ne(icell)
          C(i,j) = C(i,j) + cup
          !write(*,*) "AR85-CEA, cup=", cup
          !write(*,*) "CEA: line=",countline,ij, k, " C[ij,k]=",C(ij,k), &
          !  " C[ji,k]=",C(ji,k)
          if (present(deriv)) then
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif
       else if (key.eq."AR85-CHP") then
          ! -- Charge transfer with ionised Hydrogen
          ar85t1 = coeff(1)
          ar85t2 = coeff(2)
          ar85a = coeff(3)
          ar85b = coeff(4)
          ar85c = coeff(5)
          ar85d = coeff(6)
          if ((T(icell) >= ar85t1).and.(T(icell) <= ar85t2)) then
             t4 = T(icell)/1d4
             cup = ar85a * 1d-9 * (t4)**(ar85b) * exp(-ar85c*t4) * exp(-ar85d*EV/KBOLTZMANN/T(icell)) * &
                  Hydrogen%n(Hydrogen%Nlevel,icell) * (CM_TO_M)**3
             C(i,j) = C(i,j) + cup
             if (present(deriv)) then
                deriv(i,j) = deriv(i,j) + Cup / ne(icell)
             endif
          end if
       else if (key.eq."AR85-CHH") then
          ! Charge transfer with neutral Hydrogen
          ar85t1 = coeff(1)
          ar85t2 = coeff(2)
          ar85a = coeff(3)
          ar85b = coeff(4)
          ar85c = coeff(5)
          ar85d = coeff(6)
          if ((T(icell) >= ar85t1).and.(T(icell) <= ar85t2)) then
             t4 = T(icell)/1d4
             cdn = ar85a * 1d-9 * (t4)**(ar85b) * (1. + ar85c*exp(ar85d*t4)) * Hydrogen%n(1,icell) * (CM_TO_M)**3
             C(j,i) = C(j,i) + cdn
             if (present(deriv)) then
                deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * cdn
             endif
          end if
       else if (key.eq."BURGESS") then
          write(*,*) "BURGESS NOT CHECK"
          ! -- Electron impact ionisation from Burgess Chidichimo
          ! 1983, MNRAS, 203, 1269-1280
          de = (atom%E(j)-atom%E(i))/EV
          zz = atom%stage(i) !0 for neutrals
          betab = 0.25*(sqrt((100.*zz+91)/(4.*zz+3.))-5.)
          cbar = 2.3
          dekt = de*EV / (KBOLTZMANN*T(k))
          dekt = MIN(500., dekt)
          dekti = 1./dekt
          wlog = log(1.+dekti)
          wb = (wlog)**(betab/(1.+dekti))
          cup = 2.1715d-8 * cbar * (13.6/de)**(1.5d0) * sqrt(dekt) * E1(dekt) * wb * ne(icell) * (CM_TO_M)**3
          ! -- Fudge factor
          cup = cup * coeff(1)
          cdn = cup * atom%nstar(i,icell) / atom%nstar(j,icell)
          write(*,*) "BRUGESS, cdn=", cdn, " cup=", cup
          C(i,j) = C(i,j) + cup
          C(j,i) = C(j,i) + cdn
          if (present(deriv)) then
             deriv(j,i) = deriv(j,i) + 2.0 / ne(icell) * cdn
             deriv(i,j) = deriv(i,j) + Cup / ne(icell)
          endif
       end if

    end do loop_lines_in_file

    ! 		deallocate(TGRID, coeff)

    return
  end function collision_rates_atom_loc

  ! FUNCTION CollisionRate_old(icell, atom) result(C)
  !  ! -------------------------------------------------------------
  !  ! Computes collisional rates and fills the C matrix.
  !  ! The routine is made such that new recipe can be
  !  ! implemented easily (I think).
  !  ! C_ij = C(i,j) = C j->i
  !  !
  !  !
  !  !
  !  ! This version reads lines in a array atom%collision_lines,
  !  ! avoid reading in the atomic file. Essentially the same
  !  ! as CollisionRate_File, which is RH original.
  !  ! -------------------------------------------------------------
  !  use math, only : sq, cube, dpow
  !   type (AtomType), intent(inout) :: atom
  !   integer, intent(in) :: icell
  !   double precision, dimension(atom%Nlevel,atom%Nlevel) :: C
  !   integer :: ij, k, Nread, countline=0, colunit
  !   integer :: NTMP, n, m, ii, Nitem, i1, i2, i, j, ji, Nitem2
  !   character(len=8) :: END_OF_FILE="END     ", key
  !   real(8), dimension(:), allocatable :: TGRID, coeff
  !   double precision :: np, CC
  !   real(8) :: deltaE, C0, Cdown, Cup, gij, xj, fac, fxj
  !   integer :: Ncoef, Nrow, Nlines, k1
  !   character(len=MAX_LENGTH) :: inputline, FormatLine
  !   real(8), dimension(:,:), allocatable :: badi, cdi
  !   real(8) :: acolsh, tcolsh, aradsh, xradsh, adish, bdish
  !   real(8) :: t0sh, t1sh, summrs, tg, cdn, ccup
  !   real(8) :: ar85t1, ar85t2, ar85a, ar85b, ar85c, ar85d, t4
  !   real(8) :: de,zz,betab,cbar,dekt,dekti,wlog,wb,sumscl=0.0
  !
  !   ! -- Nitem2 is here to check
  !   ! that Nitem (read) is equals to the acutla number
  !   ! of item read Nitem2
  !   ! Actually the error will be in read(inputline,*)
  !   ! because it reads in an array of size Nitem, so
  !   ! Nitem are expected to be read and an error will be
  !   ! raise otherwise.
  !
  !   write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
  !
  !   C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) * sqrt(8.0/(PI*KBOLTZMANN))
  !
  !   ! -- Initialise the collisional matrix at each cell--
  !   C(:,:) = 0d0
  !   Nlines = size(atom%collision_lines)
  !
  !   key="        "
  !   ! -- Go throught the remaining lines in the file --
  !   ! read collisional data depending the cases of recipes.
  !
  !   !do while(key.ne.END_OF_FILE)
  !   do k1=1, Nlines
  !    countline = countline + 1
  !    !CALL getnextline(colunit, COMMENT_CHAR, FormatLine, inputline, Nread)
  !    inputline = atom%collision_lines(k1)
  !    Nread = len(inputline)
  !    ! do not go further, because after there is
  !    ! the number of grid points, which is in general
  !    ! one or two digits.
  !    ! if you intend to increase the number of grid points
  !    ! for the interpolation of the coefficients
  !    ! you have to modify the file to flush right
  !    ! everything after the TEMP keyword.
  !
  !    key = adjustl(inputline(1:len(key)))
  !    !write(*,*) trim(inputline)
  !    !write(*,*) "Line = ", countline, " key=",key,"."
  !    if (key.eq."TEMP") then
  !     read(inputline(len(key)+1:len(key)+3), *) NTMP
  !     !write(*,*) "NTMP = ", NTMP
  !     if (allocated(TGRID)) deallocate(TGRID)
  !     allocate(TGRID(NTMP))
  !     ! Nread is len(NTMP) in string format, offset
  !     ! inputline to read TGRID only, after NTMP.
  !
  !     read(inputline(len(key)+3:Nread),*) (TGRID(k), k=1,NTMP)
  !     !write(*,*) (TGRID(k), k=1,NTMP)
  !     Nitem = NTMP
  !     Nitem2 = size(TGRID)
  !
  !    else if ((key.eq."OMEGA") .or. (key.eq."CE") .or.&
  !         (key.eq."CI") .or. (key.eq."CP") .or. &
  !         (key.eq."CH0") .or. (key.eq."CH+") .or. &
  !         (key.eq."CH") .or. (key.eq."CR")) then
  !
  !     ! -- read level indices and collision coefficients --
  !     Nitem = NTMP
  !     if (allocated(coeff)) deallocate(coeff)
  !     allocate(coeff(Nitem))
  !     ! offset inputline to avoid re-reading key.
  !     ! but read(inputline,*) key, i1, i2, coeff is OK
  !     !write(*,*) key, inputline(len(key)+1:)
  !     read(inputline(len(key)+1:),*) i1, i2,(coeff(k),k=1,Nitem)
  !
  !     Nitem2 = size(coeff)
  !
  !     i = MIN(i1,i2) + 1
  !     j = MAX(i1,i2) + 1 !from 1 to Nlevel (in C, 0 to Nlevel-1)
  !     ij = (i-1)*atom%Nlevel + j!j->i
  !     ji = (j-1)*atom%Nlevel + i !i->j
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
  !       inputline = atom%collision_lines(k1+m)
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
  !       inputline = atom%collision_lines(k1+m)
  !
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
  !    if (Nitem2 .ne. Nitem) then
  !     ! -- Actually, the error should be in the reading.
  !     write(*,*) "Error, read ", Nitem2," items->expected ", &
  !      Nitem
  !     write(*,*) "This should never happen, here."
  !     write(*,*) "Exiting..."
  !     stop
  !    end if
  !
  !    ! -- END of reading, filling now ...
  !
  !    if ((key.eq."OMEGA") .or. (key.eq."CE") .or. &
  !        (key.eq."CI") .or. (key.eq."CP") .or. &
  !        (key.eq."CH0") .or. (key.eq."CH+") .or. &
  !        (key.eq."CH") .or. (key.eq."CR") ) then
  !
  !
  !                       !here Nitem is NTMP !!
  !     !CC = interp1D(TGRID, coeff, T(icell))
  !     CC = interp_dp(coeff, TGRID, T(icell))
  !
  !    end if
  !    if (key.eq."OMEGA") then
  !     ! -- Collisional excitation of ions
  !      Cdown = C0 * ne(icell) * CC / &
  !       (atom%g(j)*sqrt(T(icell)))
  !      C(j,i) = C(j,i) + Cdown
  !      ! remember the relation between Cij and Cji
  !      ! which takes place at LTE.
  !      C(i,j) = C(i,j) + Cdown * &
  !       atom%nstar(j,icell)/atom%nstar(i,icell)
  !     !!end do !! cell by cell
  !    else if (key.eq."CE") then
  !     ! -- Collisional excitation of neutrals
  !     gij = atom%g(i)/atom%g(j)
  !     !!do k=1,Nspace
  !     !write(*,*) CC * gij * sqrt(T(icell))
  !     																	  !exp(de/kT)
  !      Cdown = CC * ne(icell) * gij * sqrt(T(icell)) !CE(RH)=CE*gj/gi*ni/nj / sqrt(T)=CC
  !      !write(*,*) key, "k=",k, "Cdown = ", Cdown, C(k)
  !      !write(*,*) "ne=",ne(k), gij, "sqrt(T)=",sqrt(T(k))
  !      C(j,i) = C(j,i) + Cdown
  !      C(i,j) = C(i,j) + Cdown * &
  !       atom%nstar(j,icell)/atom%nstar(i,icell)
  !     !!end do
  !    else if (key.eq."CI") then
  !     ! -- Collisional ionisation
  !     deltaE = atom%E(j) - atom%E(i)
  !     !!do k=1,Nspace
  !      Cup = CC * ne(icell) * &
  !        exp(-deltaE/(KBOLTZMANN*T(icell))) *sqrt(T(icell))
  !      !write(*,*) key, "k=",k, "Cup = ", Cup, C(k)
  !      !write(*,*) "dE=",deltaE," exp()=",exp(-deltaE/(KBOLTZMANN*T(k)))
  !      C(i,j) = C(i,j) + Cup
  !      C(j,i) = C(j,i) + Cup * &
  !       atom%nstar(i,icell)/atom%nstar(j,icell)
  !     !!end do
  !    else if (key.eq."CR") then
  !     ! -- Collisional de-excitation by electrons
  !     !!do k=1,Nspace
  !      Cdown = ne(icell) * CC
  !      C(j,i) = C(j,i) + Cdown
  !     !!end do
  !    else if (key.eq."CP") then
  !     ! -- collisions with protons
  !     ! protons are the last level of Hydrogen atoms
  !     np = Hydrogen%n(Hydrogen%Nlevel,icell)
  !     !!do k=1,Nspace
  !      Cdown = np * CC
  !      C(j,i) = C(j,i) + Cdown
  !      C(i,j) = C(i,j) + Cdown * &
  !      atom%nstar(j,icell) / atom%nstar(i,icell)
  !     !!end do
  !    else if (key.eq."CH") then
  !     ! -- Collisions with neutral hydrogen
  !     !!do k=1,Nspace
  !      Cup = Hydrogen%n(1,icell) * CC
  !      C(i,j) = C(i,j) + Cup
  !      C(j,i) = C(j,i) + Cup * &
  !       atom%nstar(i,icell) / atom%nstar(j,icell)
  !    !! end do
  !    else if (key.eq."CH0") then
  !     ! -- Charge exchange with neutral hydrogen
  !     !!do k=1,Nspace
  !      C(j,i) = C(j,i) + Hydrogen%n(1,icell)*CC
  !     !!end do
  !    else if (key.eq."CH+") then
  !     ! -- charge exchange with protons
  !     np = Hydrogen%n(Hydrogen%Nlevel,icell)
  !     !!do k=1,Nspace
  !      C(i,j) = C(i,j) + np*CC
  !     !!end do
  !    else if (key.eq."SHULL82") then
  !     acolsh = coeff(1)
  !     tcolsh = coeff(2)
  !     aradsh = coeff(3)
  !     xradsh = coeff(4)
  !     adish = coeff(5)
  !     bdish = coeff(6)
  !     t0sh = coeff(7)
  !     t1sh = coeff(8)
  !     !!do k=1,Nspace
  !      summrs = sumscl*summers(i,j,ne(icell),atom)
  !      summrs = summrs + (1.-sumscl)
  !      tg = T(icell)
  !      cdn = aradsh * dpow(tg/1.d4, -xradsh) + &
  !        summrs*adish / tg / sqrt(tg) * exp(-t0sh/tg) * &
  !        (1. + (bdish*exp(-t1sh/tg)))
  !      cup = acolsh * sqrt(tg) * exp(-tcolsh/tg) / &
  !       (1. + 0.1*tg/tcolsh)
  !      ! -- convert from cm3 /s to m3/s
  !      cdn = cdn * ne(icell) * CUBE(CM_TO_M)
  !      cup = cup*ne(icell)*CUBE(CM_TO_M)
  !      ! -- 3-body recombination (high density limit)
  !      cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)
  !      !write(*,*) "k=",k, " cdn = ", cdn
  !      C(j,i) = C(j,i) + cdn
  !      C(i,j) = C(i,j) + cup
  !     !!end do
  !    else if (key.eq."BADNELL") then
  !     ! -- Fit for dielectronic recombination form Badnell
  !     ! First line coefficients are energies in K (from Chianti)
  !     ! Second line coefficients are coefficients (from Chianti)
  !     ! See Badnell 2006 for more details
  !
  !     !!do k =1,Nspace
  !      summrs = sumscl*summers(i,j,ne(icell),atom) + &
  !        (1.-sumscl)
  !      tg = T(icell)
  !      cdn = 0.
  !      do ii=1,Ncoef
  !       cdn = cdn + badi(2,ii) * exp(-badi(1,ii)/tg)
  !      end do
  !      cdn = cdn * dpow(tg, -1.5d0)
  !      !write(*,*) "k=",k, " cdn = ", cdn, " summrs = ",summrs, "cdn=", cdn
  !      ! -- convert from cm3/s to m3/s
  !      cdn = cdn *ne(icell) * summrs * CUBE(CM_TO_M)
  !      cup = cdn * atom%nstar(j,icell)/atom%nstar(i,icell)
  !
  !      cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)
  !
  !      C(j,i) = C(j,i) + cdn
  !      C(i,j) = C(i,j) + cup
  !      !write(*,*) "BADNELL", cdn, cup
  !     !!end do
  !     deallocate(badi)
  !    else if (key.eq."AR85-CDI") then
  !     ! -- Direct collisional ionisation
  !     !!do k=1,Nspace
  !      cup = 0.
  !      tg = T(icell)
  !      do m=1,Nrow
  !       xj = cdi(m,1) * EV/ (KBOLTZMANN*tg)
  !       fac = exp(-xj) * sqrt(xj)
  !       fxj = cdi(m,2) + cdi(m,3) * (1.+xj) + &
  !        (cdi(m,4) - xj*(cdi(m,2)+cdi(m,3)*(2.+xj)))*&
  !         fone(xj) + cdi(m,5)*xj*ftwo(xj)
  !       fxj = fxj * fac
  !       !write(*,*) fxj, ftwo(xj), fone(xj)
  !       fac = 6.69d-7 / dpow(cdi(m,1),1.5d0)
  !       cup = cup + fac*fxj*CUBE(CM_TO_M)
  !      end do
  !      if (cup.lt.0.) cup = 0.
  !      cup = cup * ne(icell)
  !      cdn = cup * atom%nstar(i,icell)/atom%nstar(j,icell)
  !      C(j,i) = C(j,i) + cdn
  !      C(i,j) = C(i,j) + cup
  !      !write(*,*) "CDU", cdn, cup
  !      !write(*,*) "CDI: line=",countline,ij, k, " C[ij,k]=",C(ij,k), &
  !      ! " C[ji,k]=",C(ji,k)
  !     !!end do
  !     deallocate(cdi)
  !    else if (key.eq."AR85-CEA") then
  !     ! -- Autoionisation
  !     !!do k=1,Nspace
  !      fac = ar85cea(i,j,icell,atom)
  !      !write(*,*) "fac=", fac
  !      cup = coeff(1)*fac*ne(icell)
  !      C(i,j) = C(i,j) + cup
  !      !write(*,*) "AR85-CEA, cup=", cup
  !      !write(*,*) "CEA: line=",countline,ij, k, " C[ij,k]=",C(ij,k), &
  !      !  " C[ji,k]=",C(ji,k)
  !     !!end do
  !    else if (key.eq."AR85-CHP") then
  !     ! -- Charge transfer with ionised Hydrogen
  !     ar85t1 = coeff(1)
  !     ar85t2 = coeff(2)
  !     ar85a = coeff(3)
  !     ar85b = coeff(4)
  !     ar85c = coeff(5)
  !     ar85d = coeff(6)
  !     !!do k =1,Nspace
  !      if ((T(icell).ge.ar85t1).and.&
  !           (T(icell).le.ar85t2)) then
  !       t4 = T(icell)/1d4
  !       cup = ar85a * 1d-9 * dpow(t4,ar85b) * &
  !         exp(-ar85c*t4) * &
  !         exp(-ar85d*EV/KBOLTZMANN/T(icell)) * &
  !         Hydrogen%n(Hydrogen%Nlevel,icell) * &
  !         CUBE(CM_TO_M)
  !       C(i,j) = C(i,j) + cup
  !      end if
  !     !!end do
  !    else if (key.eq."AR85-CHH") then
  !     ! Charge transfer with neutral Hydrogen
  !     ar85t1 = coeff(1)
  !     ar85t2 = coeff(2)
  !     ar85a = coeff(3)
  !     ar85b = coeff(4)
  !     ar85c = coeff(5)
  !     ar85d = coeff(6)
  !     !!do k=1,Nspace
  !      if ((T(icell).ge.ar85t1).and.&
  !          (T(icell).le.ar85t2)) then
  !       t4 = T(icell)/1d4
  !       cdn = ar85a * 1d-9 * dpow(t4, ar85b) * &
  !        (1. + ar85c*exp(ar85d*t4)) * &
  !        Hydrogen%n(1,icell) * CUBE(CM_TO_M)
  !       C(j,i) = C(j,i) + cdn
  !      end if
  !     !!end do
  !    else if (key.eq."BURGESS") then
  !     write(*,*) "BURGESS NOT CHECK"
  !     ! -- Electron impact ionisation from Burgess Chidichimo
  !     ! 1983, MNRAS, 203, 1269-1280
  !     de = (atom%E(j)-atom%E(i))/EV
  !     zz = atom%stage(i) !0 for neutrals
  !     betab = 0.25*(sqrt((100.*zz+91)/(4.*zz+3.))-5.)
  !     cbar = 2.3
  !     !!do k=1,Nspace
  !      dekt = de*EV / (KBOLTZMANN*T(k))
  !      dekt = MIN(500., dekt)
  !      dekti = 1./dekt
  !      wlog = log(1.+dekti)
  !      wb = dpow(wlog, betab/(1.+dekti))
  !      cup = 2.1715d-8 * cbar * dpow(13.6/de, 1.5d0) * &
  !          sqrt(dekt) * E1(dekt) * wb * ne(icell) * &
  !          CUBE(CM_TO_M)
  !      ! -- Fudge factor
  !      cup = cup * coeff(1)
  !      cdn = cup * atom%nstar(i,icell) / atom%nstar(j,icell)
  !      write(*,*) "BRUGESS, cdn=", cdn, " cup=", cup
  !      C(i,j) = C(i,j) + cup
  !      C(j,i) = C(j,i) + cdn
  !     !!end do
  !    end if
  !   end do !loop over (remaining) lines (of the atomic model)
  !
  !  deallocate(TGRID)
  !  deallocate(coeff)
  !
  !  RETURN
  !  END FUNCTION CollisionRate_old
end module collision

!-> Futur deprecation, will be computed on the fly, but his routine will be kept
!-> Or data should be kept in memory but seems diffcult to handle all cases
!  SUBROUTINE CollisionRate_File(icell, atom)
!  ! --------------------------------------------------
!  ! Computes collisional rates and fills the C matrix.
!  ! The routine is made such that new recipe can be
!  ! implemented easily (I think).
!  ! C_ij = C(i,j) = C j->i
!  ! Remeber, atom%C(Nlevel*Nlevel) now
!  ! we define indexes ij, and ji to run in the matrix
!  !
!  !
!  ! (1) atomic file is opened with openCollisionFile at the
!  ! atom%offset_coll position in the file. Before leaving
!  ! the subroutine, the cursor is reset to atom%offset_coll position
!  ! for next cell point.
!  ! The file is closed at the end of the simulation.
!  ! This trick avoid closing/reopening the file at each cell point.
!  ! --------------------------------------------------
!   type (AtomType), intent(inout) :: atom
!   integer, intent(in) :: icell
!   integer :: ij, k, Nread, countline=0, colunit, checkfseek, fseek
!   integer :: NTMP, n, m, ii, Nitem, i1, i2, i, j, ji, Nitem2
!   character(len=8) :: END_OF_FILE="END     ", key
!   real(8), dimension(:), allocatable :: TGRID, coeff
!   double precision :: C(atom%Nlevel*atom%Nlevel), np, CC
!   real(8) :: deltaE, C0, Cdown, Cup, gij, xj, fac, fxj
!   integer :: Ncoef, Nrow
!   character(len=MAX_LENGTH) :: inputline, FormatLine
!   real(8), dimension(:,:), allocatable :: badi, cdi
!   real(8) :: acolsh, tcolsh, aradsh, xradsh, adish, bdish
!   real(8) :: t0sh, t1sh, summrs, tg, cdn, ccup
!   real(8) :: ar85t1, ar85t2, ar85a, ar85b, ar85c, ar85d, t4
!   real(8) :: de,zz,betab,cbar,dekt,dekti,wlog,wb,sumscl=0.0
!
!   ! -- Nitem2 is here to check
!   ! that Nitem (read) is equals to the acutla number
!   ! of item read Nitem2
!   ! Actually the error will be in read(inputline,*)
!   ! because it reads in an array of size Nitem, so
!   ! Nitem are expected to be read and an error will be
!   ! raise otherwise.
!
!   ! file is already opened, and will be close at the end of the RT run
!   colunit = atom%colunit
!   !write(*,*) "Collision()->colunit = ", colunit, atom%colunit
!
!   write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH
!
!   C0 = ((E_RYDBERG/sqrt(M_ELECTRON)) * PI*SQ(RBOHR)) * sqrt(8.0/(PI*KBOLTZMANN))
!
!   ! -- Initialise the collisional matrix at each cell--
!   C(:) = 0d0
!   !temporary
!   atom%Ckij(icell,:) = 0d0 !C j->i = C(j,i)
!
!
!
!   key="        "
!   ! -- Go throught the remaining lines in the file --
!   ! read collisional data depending the cases of recipes.
!
!   do while(key.ne.END_OF_FILE)
!    countline = countline + 1
!    CALL getnextline(colunit, COMMENT_CHAR, FormatLine, inputline, Nread)
!
!    ! do not go further, because after there is
!    ! the number of grid points, which is in general
!    ! one or two digits.
!    ! if you intend to increase the number of grid points
!    ! for the interpolation of the coefficients
!    ! you have to modify the file to flush right
!    ! everything after the TEMP keyword.
!
!    key = adjustl(inputline(1:len(key)))
!    !write(*,*) trim(inputline)
!    !write(*,*) "Line = ", countline, " key=",key,"."
!    if (key.eq."TEMP") then
!     read(inputline(len(key)+1:len(key)+3), *) NTMP
!     !write(*,*) "NTMP = ", NTMP
!     if (allocated(TGRID)) deallocate(TGRID)
!     allocate(TGRID(NTMP))
!     ! Nread is len(NTMP) in string format, offset
!     ! inputline to read TGRID only, after NTMP.
!
!     read(inputline(len(key)+3:Nread),*) (TGRID(k), k=1,NTMP)
!     !write(*,*) (TGRID(k), k=1,NTMP)
!     Nitem = NTMP
!     Nitem2 = size(TGRID)
!
!    else if ((key.eq."OMEGA") .or. (key.eq."CE") .or.&
!         (key.eq."CI") .or. (key.eq."CP") .or. &
!         (key.eq."CH0") .or. (key.eq."CH+") .or. &
!         (key.eq."CH") .or. (key.eq."CR")) then
!
!     ! -- read level indices and collision coefficients --
!     Nitem = NTMP
!     if (allocated(coeff)) deallocate(coeff)
!     allocate(coeff(Nitem))
!     ! offset inputline to avoid re-reading key.
!     ! but read(inputline,*) key, i1, i2, coeff is OK
!     !write(*,*) key, inputline(len(key)+1:)
!     read(inputline(len(key)+1:),*) i1, i2, &
!          (coeff(k),k=1,Nitem)
!
!     Nitem2 = size(coeff)
!
!     i = MIN(i1,i2) + 1
!     j = MAX(i1,i2) + 1 !from 1 to Nlevel (in C, 0 to Nlevel-1)
!     ij = (i-1)*atom%Nlevel + j!j->i
!     ji = (j-1)*atom%Nlevel + i !i->j
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
!    if (Nitem2 .ne. Nitem) then
!     ! -- Actually, the error should be in the reading.
!     write(*,*) "Error, read ", Nitem2," items->expected ", &
!      Nitem
!     write(*,*) "This should never happen, here."
!     write(*,*) "Exiting..."
!     stop
!    end if
!
!    ! -- END of reading, filling now ...
!
!    if ((key.eq."OMEGA") .or. (key.eq."CE") .or. &
!        (key.eq."CI") .or. (key.eq."CP") .or. &
!        (key.eq."CH0") .or. (key.eq."CH+") .or. &
!        (key.eq."CH") .or. (key.eq."CR") ) then
!
!
!                       !here Nitem is NTMP !!
!     CC = interp_dp(coeff, TGRID, T(icell))
!
!    end if
!    if (key.eq."OMEGA") then
!     ! -- Collisional excitation of ions
!
!     !!do k=1, Nspace !! cell by cell
!      !! C is a constant, replace by C(k) if for all grid at once
!      !! and C(ij) by C(ij,k)
!      Cdown = C0 * ne(icell) * CC / &
!       (atom%g(j)*sqrt(T(icell)))
!      C(ij) = C(ij) + Cdown
!      ! remember the relation between Cij and Cji
!      ! which takes place at LTE.
!      C(ji) = C(ji) + Cdown * &
!       atom%nstar(j,icell)/atom%nstar(i,icell)
!     !!end do !! cell by cell
!    else if (key.eq."CE") then
!     ! -- Collisional excitation of neutrals
!     gij = atom%g(i)/atom%g(j)
!     !!do k=1,Nspace
!      Cdown = CC * ne(icell) * gij * sqrt(T(icell))
!      !write(*,*) key, "k=",k, "Cdown = ", Cdown, C(k)
!      !write(*,*) "ne=",ne(k), gij, "sqrt(T)=",sqrt(T(k))
!      C(ij) = C(ij) + Cdown
!      C(ji) = C(ji) + Cdown * &
!       atom%nstar(j,icell)/atom%nstar(i,icell)
!     !!end do
!    else if (key.eq."CI") then
!     ! -- Collisional ionisation
!     deltaE = atom%E(j) - atom%E(i)
!     !!do k=1,Nspace
!      Cup = CC * ne(icell) * &
!        exp(-deltaE/(KBOLTZMANN*T(icell))) *sqrt(T(icell))
!      !write(*,*) key, "k=",k, "Cup = ", Cup, C(k)
!      !write(*,*) "dE=",deltaE," exp()=",exp(-deltaE/(KBOLTZMANN*T(k)))
!      C(ji) = C(ji) + Cup
!      C(ij) = C(ij) + Cup * &
!       atom%nstar(i,icell)/atom%nstar(j,icell)
!     !!end do
!    else if (key.eq."CR") then
!     ! -- Collisional de-excitation by electrons
!     !!do k=1,Nspace
!      Cdown = ne(icell) * CC
!      C(ij) = C(ij) + Cdown
!     !!end do
!    else if (key.eq."CP") then
!     ! -- collisions with protons
!     ! protons are the last level of Hydrogen atoms
!     np = Hydrogen%n(Hydrogen%Nlevel,icell)
!     !!do k=1,Nspace
!      Cdown = np * CC
!      C(ij) = C(ij) + Cdown
!      C(ji) = C(ji) + Cdown * &
!      atom%nstar(j,icell) / atom%nstar(i,icell)
!     !!end do
!    else if (key.eq."CH") then
!     ! -- Collisions with neutral hydrogen
!     !!do k=1,Nspace
!      Cup = Hydrogen%n(1,icell) * CC
!      C(ji) = C(ji) + Cup
!      C(ij) = C(ij) + Cup * &
!       atom%nstar(i,icell) / atom%nstar(j,icell)
!    !! end do
!    else if (key.eq."CH0") then
!     ! -- Charge exchange with neutral hydrogen
!     !!do k=1,Nspace
!      C(ij) = C(ij) + Hydrogen%n(1,icell)*CC
!     !!end do
!    else if (key.eq."CH+") then
!     ! -- charge exchange with protons
!     np = Hydrogen%n(Hydrogen%Nlevel,icell)
!     !!do k=1,Nspace
!      C(ji) = C(ji) + np*CC
!     !!end do
!    else if (key.eq."SHULL82") then
!     acolsh = coeff(1)
!     tcolsh = coeff(2)
!     aradsh = coeff(3)
!     xradsh = coeff(4)
!     adish = coeff(5)
!     bdish = coeff(6)
!     t0sh = coeff(7)
!     t1sh = coeff(8)
!     !!do k=1,Nspace
!      summrs = sumscl*summers(i,j,ne(icell),atom)
!      summrs = summrs + (1.-sumscl)
!      tg = T(icell)
!      cdn = aradsh * dpow(tg/1.d4, -xradsh) + &
!        summrs*adish / tg / sqrt(tg) * exp(-t0sh/tg) * &
!        (1. + (bdish*exp(-t1sh/tg)))
!      cup = acolsh * sqrt(tg) * exp(-tcolsh/tg) / &
!       (1. + 0.1*tg/tcolsh)
!      ! -- convert from cm3 /s to m3/s
!      cdn = cdn * ne(icell) * CUBE(CM_TO_M)
!      cup = cup*ne(icell)*CUBE(CM_TO_M)
!      ! -- 3-body recombination (high density limit)
!      cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)
!      !write(*,*) "k=",k, " cdn = ", cdn
!      C(ij) = C(ij) + cdn
!      C(ji) = C(ji) + cup
!     !!end do
!    else if (key.eq."BADNELL") then
!     ! -- Fit for dielectronic recombination form Badnell
!     ! First line coefficients are energies in K (from Chianti)
!     ! Second line coefficients are coefficients (from Chianti)
!     ! See Badnell 2006 for more details
!
!     !!do k =1,Nspace
!      summrs = sumscl*summers(i,j,ne(icell),atom) + &
!        (1.-sumscl)
!      tg = T(icell)
!      cdn = 0.
!      do ii=1,Ncoef
!       cdn = cdn + badi(2,ii) * exp(-badi(1,ii)/tg)
!      end do
!      cdn = cdn * dpow(tg, -1.5d0)
!      !write(*,*) "k=",k, " cdn = ", cdn, " summrs = ",summrs, "cdn=", cdn
!      ! -- convert from cm3/s to m3/s
!      cdn = cdn *ne(icell) * summrs * CUBE(CM_TO_M)
!      cup = cdn * atom%nstar(j,icell)/atom%nstar(i,icell)
!
!      cdn = cdn + cup*atom%nstar(i,icell)/atom%nstar(j,icell)
!
!      C(ij) = C(ij) + cdn
!      C(ji) = C(ji) + cup
!      !write(*,*) "BADNELL", cdn, cup
!     !!end do
!     deallocate(badi)
!    else if (key.eq."AR85-CDI") then
!     ! -- Direct collisional ionisation
!     !!do k=1,Nspace
!      cup = 0.
!      tg = T(icell)
!      do m=1,Nrow
!       xj = cdi(m,1) * EV/ (KBOLTZMANN*tg)
!       fac = exp(-xj) * sqrt(xj)
!       fxj = cdi(m,2) + cdi(m,3) * (1.+xj) + &
!        (cdi(m,4) - xj*(cdi(m,2)+cdi(m,3)*(2.+xj)))*&
!         fone(xj) + cdi(m,5)*xj*ftwo(xj)
!       fxj = fxj * fac
!       !write(*,*) fxj, ftwo(xj), fone(xj)
!       fac = 6.69d-7 / dpow(cdi(m,1),1.5d0)
!       cup = cup + fac*fxj*CUBE(CM_TO_M)
!      end do
!      if (cup.lt.0.) cup = 0.
!      cup = cup * ne(icell)
!      cdn = cup * atom%nstar(i,icell)/atom%nstar(j,icell)
!      C(ij) = C(ij) + cdn
!      C(ji) = C(ji) + cup
!      !write(*,*) "CDU", cdn, cup
!      !write(*,*) "CDI: line=",countline,ij, k, " C[ij,k]=",C(ij,k), &
!      ! " C[ji,k]=",C(ji,k)
!     !!end do
!     deallocate(cdi)
!    else if (key.eq."AR85-CEA") then
!     ! -- Autoionisation
!     !!do k=1,Nspace
!      fac = ar85cea(i,j,icell,atom)
!      !write(*,*) "fac=", fac
!      cup = coeff(1)*fac*ne(icell)
!      C(ji) = C(ji) + cup
!      !write(*,*) "AR85-CEA, cup=", cup
!      !write(*,*) "CEA: line=",countline,ij, k, " C[ij,k]=",C(ij,k), &
!      !  " C[ji,k]=",C(ji,k)
!     !!end do
!    else if (key.eq."AR85-CHP") then
!     ! -- Charge transfer with ionised Hydrogen
!     ar85t1 = coeff(1)
!     ar85t2 = coeff(2)
!     ar85a = coeff(3)
!     ar85b = coeff(4)
!     ar85c = coeff(5)
!     ar85d = coeff(6)
!     !!do k =1,Nspace
!      if ((T(icell).ge.ar85t1).and.&
!           (T(icell).le.ar85t2)) then
!       t4 = T(icell)/1d4
!       cup = ar85a * 1d-9 * dpow(t4,ar85b) * &
!         exp(-ar85c*t4) * &
!         exp(-ar85d*EV/KBOLTZMANN/T(icell)) * &
!         Hydrogen%n(Hydrogen%Nlevel,icell) * &
!         CUBE(CM_TO_M)
!       C(ji) = C(ji) + cup
!      end if
!     !!end do
!    else if (key.eq."AR85-CHH") then
!     ! Charge transfer with neutral Hydrogen
!     ar85t1 = coeff(1)
!     ar85t2 = coeff(2)
!     ar85a = coeff(3)
!     ar85b = coeff(4)
!     ar85c = coeff(5)
!     ar85d = coeff(6)
!     !!do k=1,Nspace
!      if ((T(icell).ge.ar85t1).and.&
!          (T(icell).le.ar85t2)) then
!       t4 = T(icell)/1d4
!       cdn = ar85a * 1d-9 * dpow(t4, ar85b) * &
!        (1. + ar85c*exp(ar85d*t4)) * &
!        Hydrogen%n(1,icell) * CUBE(CM_TO_M)
!       C(ij) = C(ij) + cdn
!      end if
!     !!end do
!    else if (key.eq."BURGESS") then
!     write(*,*) "BURGESS NOT CHECK"
!     ! -- Electron impact ionisation from Burgess Chidichimo
!     ! 1983, MNRAS, 203, 1269-1280
!     de = (atom%E(j)-atom%E(i))/EV
!     zz = atom%stage(i) !0 for neutrals
!     betab = 0.25*(sqrt((100.*zz+91)/(4.*zz+3.))-5.)
!     cbar = 2.3
!     !!do k=1,Nspace
!      dekt = de*EV / (KBOLTZMANN*T(k))
!      dekt = MIN(500., dekt)
!      dekti = 1./dekt
!      wlog = log(1.+dekti)
!      wb = dpow(wlog, betab/(1.+dekti))
!      cup = 2.1715d-8 * cbar * dpow(13.6/de, 1.5d0) * &
!          sqrt(dekt) * E1(dekt) * wb * ne(icell) * &
!          CUBE(CM_TO_M)
!      ! -- Fudge factor
!      cup = cup * coeff(1)
!      cdn = cup * atom%nstar(i,icell) / atom%nstar(j,icell)
!      write(*,*) "BRUGESS, cdn=", cdn, " cup=", cup
!      C(ji) = C(ji) + cup
!      C(ij) = C(ij) + cdn
!     !!end do
!    end if
!   end do !loop over (remaining) lines (of the atomic model)
!
!
!  !!close(colunit) !! closed latter
!  !! reset the cursor for next cell point !!
!  checkfseek = fseek(atom%colunit,atom%offset_coll, 0)
!  if (checkfseek.gt. 0 ) then
!    write(*,*) 'fseek error'
!    stop
!  end if
!
!  atom%Ckij(icell,:) = C(:)
!  deallocate(TGRID)
!  !!deallocate(C) !! not an array anymore
!  deallocate(coeff)
!
!  RETURN
!  END SUBROUTINE CollisionRate_File

!  FUNCTION Collision_Hydrogen(icell) result(Cji)
!   integer :: icell, i, j
!   double precision :: Cji(Hydrogen%Nlevel, Hydrogen%Nlevel)
!   double precision :: nr_ji, CI(Hydrogen%Nlevel,Hydrogen%Nlevel), CE(Hydrogen%Nlevel,Hydrogen%Nlevel)
!
!    Cji(:,:) = 0d0; CI = 0d0; CE(:,:) = 0d0
!    CALL Johnson_CI(icell, CI(:,Hydrogen%Nlevel)) !bound-free i->Nlevel
!    CALL Johnson_CE(icell, CE) !among all levels
!
!    do j=1,Hydrogen%Nlevel
!     do i=1,Hydrogen%Nlevel
!      nr_ji = Hydrogen%nstar(i,icell)/Hydrogen%nstar(j,icell)
!      Cji(j,i) = CE(j,i) +  CI(i,j)  * nr_ji
!      Cji(i,j) = CE(j,i)/nr_ji + CI(i,j)
!     end do
!    end do
!    Cji(:,:) = Cji(:,:) * ne(icell)
!
!  RETURN
!  END FUNCTION Collision_Hydrogen
!
!  SUBROUTINE Johnson_CI(icell, Cje)
!  ! --------------------------------------------------- !
!   ! Ionisation rate coefficient for
!   ! Hydrogen atom, from
!   ! L.C Johnson
!   ! ApJ 74:227-236, 1972 May 15; eq. 39
!   !
!   ! ne factorised
!   !
!   ! return C(i,j) with j = Nlevel (bound-free)
!  ! --------------------------------------------------- !
!    integer, intent(in) :: icell
!    double precision, intent(out), dimension(:) :: Cje
!    integer :: i, j, Nl
!    double precision :: C0, pia0sq, rn, bn, n, An, En, yn, zn, S, Bnp
!    type (AtomType) :: atom
!
!    atom = Atoms(1)%ptr_atom
!
!    C0 = sqrt(8.*KBOLTZMANN*T(icell) / pi / M_ELECTRON)
!    pia0sq = 2d0 * pi * RBOHR**2
!
!    Nl = atom%Nlevel
!    !Hydrogen level are ordered by n increasing, except for the continuum level
!    !n = 1., but stops before Nlevel
!
!    do i=1, Nl-1 !collision from neutral states to the ground state of H+
!     n = real(i,kind=dp)
!     if (i==1) then !n=i
!      rn = 0.45
!      bn = -0.603
!     else
!      rn = 1.94*n**(-1.57)
!      bn = 1d0/n * (4. - 18.63/n + 36.24/n/n - 28.09/n/n/n)
!     end if
!
!     ! in Joules
!     En = E_RYDBERG / n / n !energie of level with different quantum number in 13.6eV: En = 13.6/n**2
!     yn = En / KBOLTZMANN / T(icell)
!     An = 32. / 3. / sqrt(3d0) / pi * n  * (g0(n)/3. + g1(n)/4. + g2(n)/5.)
!     Bnp = 2./3. * n*n * (5. + bn)
!     zn = rn + yn
!
!     S = C0 * pia0sq * (n*yn)**2 * (An*(E1(yn)/yn - E1(zn)/zn) + &
!    			(Bnp - An*log(2*n*n))*(ksi_johnson(yn)-ksi_johnson(zn)))
!    	!!write(*,*) i-1, Nl-1, "S=", S, S*exp(yn)/sqrt(T(icell))
!     !check that otherwise multiply by exp(yn)
!     Cje(i) = S !RH -> exp(yn) / sqrt(T(icell)) !per ne
!     		   !we compute it at icell so in fact we could avoid / sqrt(icell)
!     		   !as col = Cje * sqrt(T) * exp(-de). Normally dE == En/kT=y so
!     		   !it is not useful to multiply except to take into account the slightly diff
!     		   !between En and (atom%E(Nl)-atom%E(i))
!     !!write(*,*) En, (atom%E(Nl)-atom%E(i)) !Should be similar, because E(j)=13.6/j**2
!    end do
!  RETURN
!  END SUBROUTINE Johnson_CI
!
!
!  SUBROUTINE Johnson_CE(icell, Cje)
!  ! ----------------------------------------------------- !
!   ! Excitation rate coefficient for
!   ! Hydrogen atom, from
!   ! ApJ 74:227-236, 1972 May 15; eq. 36
!   !
!   ! CE = S = C(i,j) all transitions from j, i and i, j
!   !
!   ! ( -> transform C(i,j) to specifically C(j,i)
!   ! S * exp(y) * atom%g(i)/atom%g(j)
!   ! = C(i,j) * exp(hnu/kt)*gi/gj = C(i,j) * ni/nj = C(j,i)
!   ! at LTE: (gi/gj * nj/ni)  = exp(-hnu/kT) )
!  ! ----------------------------------------------------- !
!    integer, intent(in) :: icell
!    double precision, intent(out), dimension(:,:) :: Cje
!    integer :: i, j, Nl
!    double precision :: C0, pia0sq, rn, bn, n, Ennp, y, z, S, Bnnp, En
!    double precision :: np, x, fnnp, rnnp, Annp, Gaunt_bf
!    type (AtomType) :: atom
!
!    atom = Atoms(1)%ptr_atom
!
!    C0 = sqrt(8.*KBOLTZMANN*T(icell) / pi / M_ELECTRON)
!    pia0sq = 2d0 * pi * RBOHR**2
!
!    Nl = atom%Nlevel
!    !Hydrogen level are ordered by n increasing, except for the continuum level
!    !n = 1., but stops before Nlevel
!
!    do i=1, Nl-1 !collision between neutral states, n to n'
!     n = real(i,kind=dp)
!     if (i==1) then !n=i
!      rn = 0.45
!      bn = -0.603
!     else
!      rn = 1.94*n**(-1.57)
!      bn = 1d0/n * (4. - 18.63/n + 36.24/n/n - 28.09/n/n/n)
!     end if
!
!     do j=i+1, Nl-1
!      np = dble(j)!n'
!      x = 1d0 - (n/np)**2 ! = Enn'/Rdybg
!      !Gauntfactor * 32./3./sqrt(3.)/pi * n/np**3 /x**3
!      Gaunt_bf = g0(n) + g1(n)/x + g2(n)/x/x
!      fnnp = Gaunt_bf * 32./3./sqrt(3d0)/pi * n / np / np /np / x / x / x
!      rnnp = rn * x
!      Annp = 2d0 * n*n*fnnp/x
!     ! in Joules
!      En = E_RYDBERG / n / n !energie of level with different quantum number in 13.6eV = ionisation E of n
!      y = x * En / KBOLTZMANN / T(icell) !x = ratio of E/En
!      Bnnp = 4d0 * (n**4)/(np**3) / x / x * (1. + 4./3. /x + bn/x/x)
!      z = rnnp + y
!
!      S = C0 * pia0sq * n*n*y*y/x * (Annp*((1./y + 0.5)*E1(y)-(1./z + 0.5)*E1(z))+&
!      	(Bnnp-Annp*dlog(2*n*n/x))*(E2(y)/y - E2(z)/z))
!
!      Cje(j,i) = S * exp(y) * atom%g(i)/atom%g(j)
!      !!Cje(i,j) = S
!      !write(*,*) atom%E(j) - atom%E(i), En*x !because x = deltaE/En
!     end do !over np
!    end do  !over n
!
!  RETURN
!  END SUBROUTINE Johnson_CE
!
!  FUNCTION ksi_johnson(t) result(y)
!   double precision :: t, y
!   !E0
!   y = exp(-t)/t - 2d0*E1(t) + E2(t)
!
!  RETURN
!  END FUNCTION ksi_johnson
!
!  FUNCTION g0 (n) result(g)
!   double precision :: g, n
!
!   SELECT CASE (int(n))
!   CASE (1)
!    g = 1.1330
!   CASE (2)
!    g = 1.0785
!   CASE DEFAULT
!    g = 0.9935 + 0.2328/n - 0.1296 / n / n
!   END SELECT
!
!  RETURN
!  END FUNCTION g0
!
!  FUNCTION g1 (n) result(g)
!   double precision :: g, n
!
!   SELECT CASE (int(n))
!   CASE (1)
!    g = -0.4059
!   CASE (2)
!    g = -0.2319
!   CASE DEFAULT
!    g = -1d0/n * (0.6282 - 0.5598/n + 0.5299 / n / n)
!   END SELECT
!
!  RETURN
!  END FUNCTION g1
!
!  FUNCTION g2 (n) result(g)
!   double precision :: g, n
!
!   SELECT CASE (int(n))
!   CASE (1)
!    g = 0.07014
!   CASE (2)
!    g = 0.02947
!   CASE DEFAULT
!    g = 1d0/n/n * (0.3887 - 1.181 / n + 1.470 /n / n)
!   END SELECT
!
!  RETURN
!  END FUNCTION g2
