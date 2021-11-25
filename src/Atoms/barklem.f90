!See: - Anstee & O'Mara 1995, MNRAS 276, 859-866 (s-p, p-s)
!  - Barklem & O'Mara 1997, MNRAS 290, 102 (p-d, d-p)
!  - Barklem, O'Mara & Ross 1998, MNRAS 296, 1057-1060 (d-f, f-d)
!  - Barklem, O'Mara 1998, MNRAS 300, 863-871

! The routine to read the data file is adapated from RH H. Uitenbroek


module barklem

   use atom_type, only  : AtomType, parse_label, ATOM_LABEL_WIDTH
   use getline
   use math, only       : interp2D, locate, gammln
   use constant
   use mcfost_env, only : dp, mcfost_utils! convert from the relative location of atomic data
                                          ! to mcfost's environnement folders.

   implicit none

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

   type Barklem_cross_data
      integer :: N1, N2
      real(kind=dp), allocatable, dimension(:) :: neff1, neff2
      real(kind=dp), allocatable, dimension(:,:) :: cross, alpha
      character(len=2) :: barklem_type
   end type Barklem_cross_data

   contains

   subroutine dealloc_barklem(barklem)
      type(barklem_cross_data), intent(inout) :: barklem
      deallocate(barklem%neff1, barklem%neff2, barklem%cross,barklem%alpha)
      return 
   end subroutine dealloc_barklem


   subroutine init_Barklem_cross_data(btype, bs, res)
      !cross sections data for relative velocity in 1e4 m/s
      !exponent alpha of cross section unitless.
      !if data cannot be read, res is false.
      character(len=2), intent(in) :: btype
      type (barklem_cross_data), intent(inout) :: bs
      logical, intent(out) :: res
      real(kind=dp) :: neff1_0, neff2_0
      integer :: j, i, n, Nread
      character(len=MAX_LENGTH) :: inputline, FomatLine !len=120
      character(len=50) :: barklem_data
      character(len=1) :: COMMENT_CHAR = "c"
      FomatLine = "(1A120)"

      select case (btype)
         case ("SP")
         !write(*,*) "SP Barklem's type"
            bs%N1 = BARKLEM_SP_NS
            bs%N2 = BARKLEM_SP_NP
            neff1_0 = BARKLEM_SP_NEFF1
            neff2_0 = BARKLEM_SP_NEFF2
            barklem_data = BARKLEM_SP_DATA
         case ("PD")
         !write(*,*) "PD Barklem's type"
            bs%N1 = BARKLEM_PD_NP
            bs%N2 = BARKLEM_PD_ND
            neff1_0 = BARKLEM_PD_NEFF1
            neff2_0 = BARKLEM_PD_NEFF2
         barklem_data = BARKLEM_PD_DATA
         case ("DF")
       !write(*,*) "DF Barklem's type"
            bs%N1 = BARKLEM_DF_ND
            bs%N2 = BARKLEM_DF_NF
            neff1_0 = BARKLEM_DF_NEFF1
            neff2_0 = BARKLEM_DF_NEFF2
            barklem_data = BARKLEM_DF_DATA
         case default
            write(*,'((1A2) " Unrocognised Barklem transition")') btype
            res = .false.
            return
      end select

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
            call getnextline(2, COMMENT_CHAR,FomatLine,&
               inputline, Nread)
            read(inputline,*) (bs%cross(i,j), j=1,bs%N2)
         else
            call getnextline(2, COMMENT_CHAR,FomatLine,&
               inputline, Nread)
            read(inputline,*) (bs%alpha(i-bs%N1,j), j=1,bs%N2)
         end if
      end do

      close(2)
      res = .true.
      return
   end subroutine init_Barklem_cross_data


   subroutine get_Barklem_cross_data(atom, kr, res)
      !if atom%lines(kr)%cvdWaals(1) > 1.0 = xxx.yyy
      ! use xxx as the cross section and yyy as the exponent.
      ! Otherwise attempt to obtaine these values from the Barklem data.
      ! If not possible, Unsold is used (like for Helium).
      integer, intent(in) :: kr
      type (AtomType), intent(inout) :: atom
      type (barklem_cross_data) :: bs
      logical, intent(out) :: res
      logical :: determined, lread_from_table
      real :: Si, Sj
      real :: Jj, Ji
      real(kind=dp) :: vref = 1d4 !m/s
      real(kind=dp) :: neff1, neff2, nefftmp
      real(kind=dp) :: findex1, findex2, mu, vmean, sqa0
      real(kind=dp) :: cross, E_Rydberg2, deltaEi, deltaEj, expo
      integer :: Li, Lj, i, j, Z, ic, k, index, index1
      character(len=ATOM_LABEL_WIDTH) :: label_i, label_j
      real(kind=dp) :: testcc

      i = Atom%lines(kr)%i
      j = Atom%lines(kr)%j

      label_i = atom%label(i)
      label_j = atom%label(j)

      ! Init Unsold for He, does not use what's in the file
      atom%lines(kr)%cvdWaals(3) = 1.0
      atom%lines(kr)%cvdWaals(4) = 0.0

      mu = AMU/(1.0/atomic_weights(1) + 1./atom%weight)
      vmean = sqrt(8.0*KBOLTZMANN/pi/mu) !missing a sqrt(T)

      sqa0 = RBOHR**2

      lread_from_table = .false.
      if (atom%lines(kr)%cvdWaals(1) > 1.0) then
         cross = real( int(atom%lines(kr)%cvdWaals(1)) )
         expo = real( atom%lines(kr)%cvdWaals(1)-int(atom%lines(kr)%cvdWaals(1)) )
         atom%lines(kr)%cvdWaals(2) = expo
      else !read from table
         !after because atom%lines(kr)%cvdWaals(1) might be > 1
         atom%lines(kr)%cvdWaals(2) = 0.0
         atom%lines(kr)%cvdWaals(1) = 1.0

         !Barklem only for neutrals
         ! stage is 0 for neutrals, 1 for singly ionised
         if (atom%stage(i).lt.0) then
            res = .false. !use Unsold
            return
         end if

         !parse_label
         call parse_label(label_i, atom%g(i), Si, Li,Ji, determined)
         call parse_label(label_j, atom%g(j), Sj, Lj,Jj, determined)
         if (.not. determined) then
            write(*,*) "WARNING: Barklem, cannot parse label, using Unsold !"
            res = .false.
            return
         endif 
         !write(*,'("Si="(1F2.2)", Li="(1I2)", Ji="(1F2.2))') Si, Li, Ji
         !write(*,'("Sj="(1F2.2)", Lj="(1I2)", Jj="(1F2.2))') Sj, Lj, Jj

         if ((Li.eq.0 .and. Lj.eq.1) .or. &
               (Lj.eq.0 .and. Li.eq.1)) then!S and P orbits
            call init_Barklem_cross_data("SP", bs, res)
         else if ((Li.eq.1 .and. lj.eq.2) .or. &
            (Lj.eq.1 .and. li.eq.2)) then !P and D orbits
            call init_Barklem_cross_data("PD", bs, res)
         else if ((Li.eq.2 .and. lj.eq.3) .or. &
            (Lj.eq.2 .and. li.eq.3)) then !D and F orbits
            call init_Barklem_cross_data("DF", bs, res)
         end if

         if (.not. res ) then
            write(*,*) "WARNING: Barklem, cannot get data, using Unsold !"          
            return
         endif
         lread_from_table = .true.
         !index of the appropriate continuum level
         Z = atom%stage(j)+1

         ic = j + 1
         do while (atom%stage(ic) .lt. Z)
            ic = ic + 1
         end do

         deltaEi = atom%E(ic) - atom%E(i)
         deltaEj = atom%E(ic) - atom%E(j)
         E_Rydberg2 = E_RYDBERG / (1. + M_ELECTRON/(atom%weight*AMU) )

         neff1 = Z * sqrt(E_Rydberg2/deltaEi)
         neff2 = Z * sqrt(E_Rydberg2/deltaEj)
         ! REMEMBER -> dE=(-) (Zeff/neff)**2 * E_RYDBERG(=13.6eV)
         ! so neff = sqrt(E_RYDBERG/dE) * Zeff
         !https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas-notation-data-and-formulas/atomic-spectroscopy-term

         !swamp double
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
            return       
         end if

         index = locate(bs%neff1,neff1)!minloc(abs(bs%neff1 - neff1), 1)
         findex1 = index + &
            (neff1 - bs%neff1(index))/BARKLEM_DELTA_NEFF

         index1 = index

         if ((neff2 .lt.bs%neff2(1)) .or. &
            (neff2 .gt. bs%neff2(bs%N2))) then
            res = .false.
            return
         end if

         index = locate(bs%neff2,neff2)
         findex2 = index + &
            (neff2 - bs%neff2(index))/BARKLEM_DELTA_NEFF

         !get cross section in atomic unit
         atom%lines(kr)%cvdWaals(1) = &
            interp2D(bs%neff1,bs%neff2,bs%cross,neff1,neff2)

         !get exponent
         atom%lines(kr)%cvdWaals(2) = &
            interp2D(bs%neff1,bs%neff2,bs%alpha,neff1,neff2)

      endif !read from table

      !comptue the actual cross-section with T factorised out.

      cross = 2.0*sqa0 * (4.0/pi)**(0.5*atom%lines(kr)%cvdWaals(2)) *  &
         vmean * (vmean/vref)**(-atom%lines(kr)%cvdWaals(2)) 

      atom%lines(kr)%cvdWaals(1) = atom%lines(kr)%cvdWaals(1) * cross * & 
         exp(gammln(2.0-0.5*atom%lines(kr)%cvdWaals(2)))

      res = .true.
      if (lread_from_table) call dealloc_barklem(bs)
      return
   end subroutine get_Barklem_cross_data

end module barklem