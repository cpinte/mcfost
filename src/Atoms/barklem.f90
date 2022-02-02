!See: - Anstee & O'Mara 1995, MNRAS 276, 859-866 (s-p, p-s)
!  - Barklem & O'Mara 1997, MNRAS 290, 102 (p-d, d-p)
!  - Barklem, O'Mara & Ross 1998, MNRAS 296, 1057-1060 (d-f, f-d)
!  - Barklem, O'Mara 1998, MNRAS 300, 863-871

! Adapted from RH H. Uitenbroek


MODULE barklem

  use atom_type, only : AtomType, determinate, ATOM_LABEL_WIDTH
  use getline
  use math, only : interp2D, locate, gammln
  use constant

  !MCFOST's originals
  use mcfost_env, only : dp, mcfost_utils! convert from the relative location of atomic data
  ! to mcfost's environnement folders.

  IMPLICIT NONE

  character(len=*), parameter :: BARKLEM_SP_DATA="/Atoms/Barklem_spdata.dat"
  integer, parameter :: BARKLEM_SP_NS=21
  integer, parameter :: BARKLEM_SP_NP=18
  real, parameter :: BARKLEM_SP_NEFF1=1.0
  real, parameter :: BARKLEM_SP_NEFF2=1.3

  character(len=*), parameter :: BARKLEM_PD_DATA="/Atoms/Barklem_pddata.dat"
  integer, parameter :: BARKLEM_PD_NP=18
  integer, parameter :: BARKLEM_PD_ND=18
  real, parameter :: BARKLEM_PD_NEFF1=1.3
  real, parameter :: BARKLEM_PD_NEFF2=2.3

  character(len=*), parameter :: BARKLEM_DF_DATA="/Atoms/Barklem_dfdata.dat"
  integer, parameter :: BARKLEM_DF_ND=18
  integer, parameter :: BARKLEM_DF_NF=18
  real, parameter :: BARKLEM_DF_NEFF1=2.3
  real, parameter :: BARKLEM_DF_NEFF2=3.3

  real, parameter :: BARKLEM_DELTA_NEFF=0.1


  TYPE BarklemType
     integer :: N1, N2
     real(kind=dp), allocatable, dimension(:) :: neff1, neff2
     real(kind=dp), allocatable, dimension(:,:) :: cross, alpha
     character(len=2) :: barklem_type
  END TYPE BarklemType

CONTAINS

  !cross sections data in 10e4 m/s
  !exponent alpha of cross section unitless
  SUBROUTINE initBarklem(btype, bs, res)
    character(len=2), intent(in) :: btype
    TYPE (BarklemType), intent(inout) :: bs
    logical, intent(out) :: res
    real(kind=dp) :: neff1_0, neff2_0
    integer :: j, i, n, Nread
    character(len=MAX_LENGTH) :: inputline, FomatLine !len=120
    character(len=50) :: barklem_data
    character(len=1) :: COMMENT_CHAR = "c"
    FomatLine = "(1A120)"

    SELECT CASE (btype)
    CASE ("SP")
       !write(*,*) "SP Barklem's type"
       bs%N1 = BARKLEM_SP_NS
       bs%N2 = BARKLEM_SP_NP
       neff1_0 = BARKLEM_SP_NEFF1
       neff2_0 = BARKLEM_SP_NEFF2
       barklem_data = BARKLEM_SP_DATA
    CASE ("PD")
       !write(*,*) "PD Barklem's type"
       bs%N1 = BARKLEM_PD_NP
       bs%N2 = BARKLEM_PD_ND
       neff1_0 = BARKLEM_PD_NEFF1
       neff2_0 = BARKLEM_PD_NEFF2
       barklem_data = BARKLEM_PD_DATA
    CASE ("DF")
       !write(*,*) "DF Barklem's type"
       bs%N1 = BARKLEM_DF_ND
       bs%N2 = BARKLEM_DF_NF
       neff1_0 = BARKLEM_DF_NEFF1
       neff2_0 = BARKLEM_DF_NEFF2
       barklem_data = BARKLEM_DF_DATA
    CASE DEFAULT
       write(*,'((1A2) " Unrocognised Barklem transition")') btype
       res = .false.
       RETURN
    END SELECT

    allocate(bs%neff1(bs%N1))
    allocate(bs%neff2(bs%N2))

    do n=1,bs%N1
       bs%neff1(n) = neff1_0 + n*BARKLEM_DELTA_NEFF
    end do

    do n=1,bs%N2
       bs%neff2(n) = neff2_0 + n*BARKLEM_DELTA_NEFF
    end do

    allocate(bs%cross(bs%N1,bs%N2))
    allocate(bs%alpha(bs%N1,bs%N2))

    ! read data
    open(unit=2,file=trim(mcfost_utils)//trim(barklem_data),status="old")

    do i=1,bs%N1+bs%N1
       if (i.le.bs%N1) then
          CALL getnextline(2, COMMENT_CHAR,FomatLine,&
               inputline, Nread)
          read(inputline,*) (bs%cross(i,j), j=1,bs%N2)
          !write(*,*) "Barklem cross-sections: " , &
          !(bs%cross(i,j), j=1,bs%N2)
       else
          CALL getnextline(2, COMMENT_CHAR,FomatLine,&
               inputline, Nread)
          read(inputline,*) (bs%alpha(i-bs%N1,j), j=1,bs%N2)
          !write(*,*) "Barklem alpha: " , &
          !(bs%alpha(i,j), j=1,bs%N2)
       end if
    end do

    close(2)
    res = .true.
    RETURN
  END SUBROUTINE initBarklem



  SUBROUTINE getBarklem(atom, kr, res)
    integer, intent(in) :: kr
    TYPE (AtomType), intent(inout) :: atom
    TYPE (BarklemType) :: bs
    logical, intent(out) :: res
    logical :: determined
    real :: Si, Sj
    real :: Jj, Ji
    real(kind=dp) :: neff1, neff2, nefftmp
    real(kind=dp) :: findex1, findex2, reducemass, meanvelocity
    real(kind=dp) :: crossmean, E_Rydberg2, deltaEi, deltaEj
    integer :: Li, Lj, i, j, Z, ic, k, index, index1
    character(len=ATOM_LABEL_WIDTH) :: label_i, label_j
    real(kind=dp) :: testcc

    i = Atom%lines(kr)%i
    j = Atom%lines(kr)%j

    label_i = atom%label(i)
    label_j = atom%label(j)

    !Barklem only for neutrals
    ! stage is 0 for neutrals, 1 for singly ionised
    if (atom%stage(i).lt.0) then
       res = .false.
       RETURN
    end if

    CALL determinate(label_i, atom%g(i), Si, Li,Ji, determined)
    CALL determinate(label_j, atom%g(j), Sj, Lj,Jj, determined)
    if (.not. determined) res = .false.

    !write(*,'("Si="(1F2.2)", Li="(1I2)", Ji="(1F2.2))') Si, Li, Ji
    !write(*,'("Sj="(1F2.2)", Lj="(1I2)", Jj="(1F2.2))') Sj, Lj, Jj
    write(*,*) "Si=",Si," Li=",Li," Ji=",Ji
    write(*,*) "Sj=",Sj," Lj=",Lj," Jj=",Jj

    if (determined) then
       if ((Li.eq.0 .and. Lj.eq.1) .or. &
            (Lj.eq.0 .and. Li.eq.1)) then!S and P orbits
          CALL initBarklem("SP", bs, res)
       else if ((Li.eq.1 .and. lj.eq.2) .or. &
            (Lj.eq.1 .and. li.eq.2)) then !P and D orbits
          CALL initBarklem("PD", bs, res)
       else if ((Li.eq.2 .and. lj.eq.3) .or. &
            (Lj.eq.2 .and. li.eq.3)) then !D and F orbits
          CALL initBarklem("DF", bs, res)
       end if
    end if

    if (.not. res .and. .not. determined) then
       res=.false.
       RETURN
    end if

    !index of the appropriate continuum level
    Z = atom%stage(j)+1
    ! here, we do not add + 1 wrt to the C indexing
    ! because it is not an index but a physical Z ?

    ic = j + 1
    do while (atom%stage(ic) .lt. Z)
       ic = ic + 1
    end do
    !write(*,*) i, j, Z, ic, size(atom%E)

    deltaEi = atom%E(ic) - atom%E(i)
    deltaEj = atom%E(ic) - atom%E(j)
    E_Rydberg2 = E_RYDBERG / &
         (1. + M_ELECTRON/(atom%weight*AMU) )
    !write(*,*) deltaEI, deltaEj, E_Rydberg2

    ! this explains why we do not add +1 to the Z wrt to the C indexing
    neff1 = Z * sqrt(E_Rydberg2/deltaEi)
    neff2 = Z * sqrt(E_Rydberg2/deltaEj)
    ! REMEMBER -> dE=(-) (Zeff/neff)**2 * E_RYDBERG(=13.6eV)
    ! so neff = sqrt(E_RYDBERG/dE) * Zeff
    !https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy-term

    if (Li.gt.Lj) then
       nefftmp = neff1
       neff1 = neff2
       neff2 = nefftmp
    end if


    if ((neff1 .lt. bs%neff1(1)) .or. &
         neff1 .gt. (bs%neff1(bs%N1))) then
       !write(*,*) neff1, bs%neff1(1), bs%neff1(bs%N1)
       !write(*,*) neff2, bs%neff2(1), bs%neff2(bs%N2)
       write(*,*) "neff outside domain, use Unsold"
       res = .false.
       RETURN
    end if

    ! the Following findexes are used in cubicconvolution
    ! algorithm where interpolation is done over pixels of
    ! the input table. They are define so that the pixel
    ! is close to the real value neff at a given error delta_neff.
    ! Wth spline interpolation, interpolation is done
    ! between neff1 and neff2.

    index = locate(bs%neff1,neff1)!minloc(abs(bs%neff1 - neff1), 1)
    !!CALL Hunt(bs%neff1, neff1, index)
    findex1 = index + &
         (neff1 - bs%neff1(index))/BARKLEM_DELTA_NEFF

    index1 = index

    if ((neff2 .lt.bs%neff2(1)) .or. &
         (neff2 .gt. bs%neff2(bs%N2))) then
       res = .false.
       RETURN
    end if

    index = locate(bs%neff2,neff2)
    findex2 = index + &
         (neff2 - bs%neff2(index))/BARKLEM_DELTA_NEFF

    atom%lines(kr)%cvdWaals(1) = &
         interp2D(bs%neff1,bs%neff2,bs%cross,neff1,neff2)

    atom%lines(kr)%cvdWaals(2) = &
         interp2D(bs%neff1,bs%neff2,bs%alpha,neff1,neff2)

    write(*,*) "neff1_t = ", bs%neff1(index1), " neff2_t = ", &
         bs%neff2(index)
    write(*,*) "index1 = ", index1, " index2 = ", index
    write(*,*) "neff1 = ", neff1, " neff2 = ", neff2
    write(*,*) "findex1 = ", findex1 , " findex2 = ", findex2
    write(*,*) "c = ", bs%cross(index1,index), " a = ", bs%alpha(index1,index)
    write(*,*) "c_i = ", atom%lines(kr)%cvdWaals(1), &
         " a_i = ", atom%lines(kr)%cvdWaals(2)
    !-->RH
    !neff1_t = 1.600000, neff2_t = 2.100000
    !index1 = 6, index2=8
    !neff1 = 1.627083, neff2 = 2.116624
    !findex1 = 6.270830, findex2 = 8.166237
    !c = 396.000000, a = 0.283000
    !c_i = 407.162966, a_i = 0.271340


    reducemass = AMU/(1.0/atomic_weights(1) + &
         1./atom%weight)
    meanvelocity = sqrt(8.*KBOLTZMANN / (PI*reducemass))

    crossmean = (RBOHR)**2 * &
         (meanvelocity/1d4)**(-atom%lines(kr)%cvdWaals(2))

    write(*,*) "reducemass, meanvelo, crossmean"
    write(*,*) reducemass, meanvelocity, crossmean

    atom%lines(kr)%cvdWaals(1) = atom%lines(kr)%cvdWaals(1)*&
         2.*(4./PI)**(atom%lines(kr)%cvdWaals(2)) * &
         exp(gammln(4.-atom%lines(kr)%cvdWaals(2)/2.)) *&
         meanvelocity * crossmean

    ! Use UNSOLD for the contribution of Helium atoms

    atom%lines(kr)%cvdWaals(3) = 1.0;
    atom%lines(kr)%cvdWaals(4) = 0.0;

    write(*,*) "VdWaals coeffs"
    write(*,*) (atom%lines(kr)%cvdWaals(j), j=1,4)

    res = .true.
    RETURN
  END SUBROUTINE getBarklem

END MODULE barklem
