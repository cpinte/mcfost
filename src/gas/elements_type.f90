module elements_type
    !Small module to keep track of data of elements either bound to molecules or in atoms
    !sets the abundance of the model, used for atom rt and chemical equilibrium if 
    !molecules and atoms overlap.

    use mcfost_env, only : dp, mcfost_utils
    use constantes, only : hp, c_light, cm_to_m
    use utils !do not forget to add the math module un int

    implicit none 

    type Element

        character(len=2) :: ID
        logical :: abundance_set
        integer :: Nstage, Nmolecule ! number of Molecules having an element
        integer, allocatable, dimension(:)  :: mol_index !track molecules in which Element is present
        real(kind=dp) :: weight, abund, massf
        real(kind=dp), allocatable, dimension(:)  :: ionpot
        real(kind=dp), allocatable, dimension(:,:)  :: pf
        integer :: nm = 0 !-> index in Atoms, of a given elements if any (>0)!

    end type Element

    type elemPointerArray
        type(Element), pointer :: p => NULL()
    end type elemPointerArray  


    ! type (elemPointerArray), dimension(:), allocatable :: Elements
    type (Element), dimension(:), allocatable :: elems

    !file where the abundance of all elements are
    integer, parameter :: NELEM_WEIGHTS = 99
    character(len=50), parameter :: ABUNDANCE_FILE="/Atoms/abundance.input"
    !partition function file, by default Kurucz
    character(len=50), parameter :: KURUCZ_PF_FILE="/Atoms/pf_Kurucz.fits.gz"

    real :: metallicity = 0.0
    real(kind=dp) :: totalAbund, avgWeight, wght_per_H
    real(kind=dp), dimension(:), allocatable :: Tpf
    integer :: Nelem, Npf
    real(kind=dp), parameter :: phi_min_limit = 1d-100!tiny_dp!1d-100 !1d-50, tiny_dp
    integer :: Max_ionisation_stage


  ! Atomic properties constants

    real, dimension(NELEM_WEIGHTS) :: atomic_weights
    !starts at 1 for H, ends at NELEM_WEIGHTS
    DATA atomic_weights /1.008,4.003,6.939,9.013,10.81,12.01,  &
         14.01,16.0,19.0,20.18,22.990,24.310,26.98,28.09,30.98,  &
         32.07,35.45,39.95,39.1,40.08,44.96,47.9,50.94,52.0,   &
         54.94,55.85,58.94,58.71,63.55,65.37,69.72,72.6,74.92, &
         78.96,79.91,83.8,85.48,87.63,88.91,91.22,92.91,95.95, &
         99.0,101.1,102.9,106.4,107.9,112.4,114.8,118.7,121.8, &
         127.6,126.9,131.3,132.9,137.4,138.9,140.1,140.9,      &
         144.3,147.0,150.4,152.0,157.3,158.9,162.5,164.9,      &
         167.3,168.9,173.0,175.0,178.5,181.0,183.9,186.3,      &
         190.2,192.2,195.1,197.0,200.6,204.4,207.2,209.0,      &
         210.0,211.0,222.0,223.0,226.1,227.1,232.0,231.0,      &
         238.0,237.0,244.0,243.0,247.0,247.0,251.0, 254.0/
  
    character(len=2), dimension(NELEM_WEIGHTS) :: elemental_ID
  
    DATA elemental_ID /'H ','He','Li','Be','B ','C ','N ','O ',    &
         'F ','Ne','Na','Mg','Al','Si','P ','S ','Cl', &
         'Ar','K ','Ca','Sc','Ti','V','Cr','Mn','Fe',  &
         'Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br', &
         'Kr','Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru', &
         'Rh','Pd','Ag','Cd','In','Sn','Sb','Te','I ', &
         'Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm', &
         'Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu', &
         'Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
         'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac', &
         'Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es'/
  


    contains

    subroutine dealloc_elements()
        integer :: n

        deallocate(Tpf)
        deallocate(elems)
        !individual deallocation ?
        ! do n=1,size(elements)

        ! enddo
        return
    end subroutine dealloc_elements

    !to do, add a more flexible partition function.
    subroutine read_abundance()
        !This routine reads the abundance of elements listed in the
        ! abundance.input file, and fill the elems structure.
        !
        !Elements in Abundance file are read in order of Z, so that H is 1, He is 2 ect
        !
        ! Also read partition function

        integer :: EOF, n, blocksize, unit, i, j, syst_status
        integer :: NAXIST(1), naxis_found, hdutype,  Nread
        character(len=256) :: some_comments
        logical :: anynull
        character(len=4) :: charID
        character(len=10) :: inputline, FormatLine
        real :: A
        integer :: Nstage, Z, shape_pf(2)
        real(kind=dp), allocatable, dimension(:,:) :: data_krz
        integer, parameter :: Nstage_max = 100
        integer :: NAXISPF(2)
    
        eof = 0
        !Start reading partition function by first reading the grid
        !get unique unit number
        call ftgiou(unit,EOF)
        ! open fits file in readmode'
        call ftopen(unit, TRIM(mcfost_utils)//TRIM(KURUCZ_PF_FILE), 0, blocksize, EOF)
        !some check about the first axis!
        !get first axis
        call ftgknj(unit, 'NAXIS', 1, 1, NAXIST, naxis_found, EOF)
        !read dimension of partition function table
        call ftgkyj(unit, "NPOINTS", Npf, some_comments, EOF)
        if (NAXIST(1) /= Npf) then
           write(*,*) "Error in reading pf data !"
           write(*,*) "NAXIST=",NAXIST," NPOINTS=", Npf
           write(*,*) "exiting... !"
           stop
        end if
    
        allocate(Tpf(Npf))
        !now read temperature grid
        call ftgpvd(unit,1,1,Npf,-999,Tpf,anynull,EOF) !d because double !!
        shape_pf(1) = Npf+1

        allocate(data_krz(Npf+1,Nstage_max))
        data_krz = 0.0

        !do not close yet
        !will be closed after all data are read!
    
    
        !read abundances
        write(FormatLine,'("(1"A,I3")")') "A", 10
    
        open(unit=1, file=TRIM(mcfost_utils)//TRIM(ABUNDANCE_FILE),status="old")
        call read_line(1, FormatLine, inputline, Nread)
        read(inputline,*) Nelem
        write(*,*) " Reading abundances of ", Nelem, " elements..."
        if (Nelem < 3) write(*,*) " WARNING:, not that using Nelem < 3 will cause problem", &
             " in solvene.f90 or other functions assuming abundances of metal of to 26 are known!!"


        allocate(elems(Nelem))
        totalAbund = 0.0
        avgWeight = 0.0
        do n=1, Nelem !H is first, then He ect
    
           call read_line(1, FormatLine, inputline, Nread)
           read(inputline,*) charID, A
    
           elems(n)%weight=atomic_weights(n)
           elems(n)%ID = charID(1:2) !supposed to be lowercases!!!!!
           elems(n)%abund = 10.**(A-12.0)
    
           if (A <= -99) elems(n)%abund = 0.0
    
           elems(n)%abundance_set = .true.

           totalAbund = totalAbund+elems(n)%abund
           avgWeight = avgWeight+elems(n)%abund*elems(n)%weight
    
           !read pf data, the fits is not closed yet. 
           call FTMAHD(unit,n+1,hdutype,EOF)
           if (EOF.ne.0) then
                write(*,*) "error reading partition function"
                write(*,*) "error to change HDU at ", n+1
                write(*,*) "exiting..."
                stop
            end if
            call ftgkyj(unit, "NSTAGE", Nstage, some_comments, EOF) !replace j with d if read double instead of int
            call ftgkyj(unit, "Z", Z, some_comments, EOF)
            !j'arrive pas a lire ce string, but it is not necessary
            !but should asks christophe to know how to do
            !call ftgkyj(unit, "ELEM", AtomID, commentfits, EOF)
     
            !write(*,*) "Nstage", Nstage, "Z=", Z, "for elem:", elemental_ID(code)
            !write(*,*) "Npf points = ", Npf
            !now read partition function
            call ftgknj(unit, 'NAXIS', 1, 2, NAXISPF, naxis_found, EOF)
            if (NAXISPF(1).ne.Npf.and.NAXISPF(2).ne.Nstage) then
                write(*,*) "error reading partition function"
                write(*,*) "NAXISPF(len=2)=",NAXISPF," NPOINTS=", &
                    Npf, " NSTAGE=", Nstage
                write(*,*) "exiting"
                stop
            end if
            !data are actually transposed from my fits to fortran
            elems(n)%Nstage = Nstage
            allocate(elems(n)%ionpot(elems(n)%Nstage))
            !-> all T for a given ionisation stage are contiguous in memory!
            allocate(elems(n)%pf(Npf,elems(n)%Nstage))
            !now read the value of pf for that atom
            ! remember: pf(:,1) =  ion potentials
            !           pf(:,2:) = partition functions for all ion pots.

            shape_pf(2) = elems(n)%Nstage
            call FTG2Dd(unit,1,-999,shape_pf,Npf+1,elems(n)%Nstage,data_krz(:,1:elems(n)%Nstage),anynull,EOF)
            !do i=1,Nstage
            ! write(*,*) "potential in cm-1 for elem ", elemental_ID(code),":", &
            !            data_krz(1,i)
            !fill ionpot array and pf array for that elem
            ! store value in Joule
            elems(n)%ionpot(:) = data_krz(1,:) * hp*c_light / (CM_TO_M)
            !write(*,*) "chi(i) = ", data_krz(1,i)
            !do not forget to write partition function has (Nstage,Npf)
            !instead of (Npf, Nstage) as read by FTG2D
            !store value in logarithm of partition function, for interpolation
            do i=1,Nstage
                do j=2,Npf+1
               !!if (code.eq.1) write(*,*) 'pf=', data_krz(j,i)
     !!!pf(i,j-1) = LOG10(data_krz(j,i)) !29/12/2019, using neperien log
                    elems(n)%pf(j-1,i) = LOG(data_krz(j,i))
               !!if (code.eq.1) write(*,*) 'pf10powLog=', 10**(pf(i,j-1))
                end do
            end do
        end do !over elem read
        close(unit=1)
        call ftclos(unit, EOF)
        ! free the unit number
        call ftfiou(unit, EOF)
        deallocate(data_krz)
    
        wght_per_H = avgWeight
        avgWeight = avgWeight/totalAbund
        write(*,*) "Total Abundance in the atmosphere = ", totalAbund
        write(*,*) "Total average weight = ", avgWeight
        write(*,*) "Weight per Hydrogen = ", wght_per_H !i.e., per AMU = total mass
        write(*,*) ""
    
        !store also the mass fraction of each element
        do n=1, Nelem
           elems(n)%massf = elems(n)%weight*elems(n)%Abund/wght_per_H
        enddo

        Max_ionisation_stage =  maxval(elems(:)%Nstage)
    
        return
    end subroutine read_abundance

    subroutine atom_pos(Z, row, col)
        ! ---------------------------------------------------
        ! Position in the periodic table of an atom according
        ! to this Z number.
        ! min row is 1, max row is 87 Francium
        ! No special ceses yet for Lanthanides and Actinides
        ! min col is 1, max col is 18
        ! beware He, is on the 18 column, because they are 16
        ! empty columns between H and He
        ! ---------------------------------------------------
     
           integer i
           integer, intent(in) :: Z
           integer, intent(out) :: row, col
           integer :: istart(7)
     
           row = 0
           col = 0
     
           istart = (/1, 3, 11, 19, 37, 55, 87/) !first Z of each rows
     
           ! -- find row
           do i=1,6
              if ((Z .ge. istart(i)) .and. (Z .lt. istart(i + 1))) then
                 ! we are on the row istart(i)
                 row=i
                 exit
              end if
           end do
           ! find column on the row
           col = Z - istart(i) + 1
     
           ! for relative position just comment out the following lines
           ! or take col=col-10 (or-18 for He) for these cases.
     
           ! special case of Z=2 for Helium because there are 16
           ! empty columns between Hydrogen and Helium
           if (Z.eq.2) then
              col = col + 16
              ! ten empty lines between Be and B
              ! ten empty lines between Mg and Al
           else if (((istart(i).eq.3).or.(istart(i)).eq.11).and.&
              (Z.gt.istart(i)+1)) then
              !row = istart(i)
              col = col+10
           end if
           return
    end subroutine atom_pos

    function get_pf(elem, j, temp) result (Uk)
        ! ----------------------------------------------------------------------!
        ! Interpolate the partition function of Element elem in ionisation stage
        ! j at temperature temp
        ! 24/02/2023: added linear extrapolation for points outside range.
        ! ----------------------------------------------------------------------!
    
        type(Element), intent(in) :: elem
        integer, intent(in) :: j
        real(kind=dp), intent(in) :: temp
        real(kind=dp) :: Uk
        real(kind=dp) :: Uka(1), tp(1)

        ! Uk = exp( interp(elem%pf(:,j),Tpf,temp) )
        ! elem%pf(:,j) is ln(U)
    
        !->faster
        !out of bound the function return 0 not the inner (outer) bound.
          tp(1) = temp
          Uka(:) = linear_1D_sorted(Npf, Tpf, elem%pf(:,j), 1, tp)
        !   Uk = exp(Uka(1))
        !   if( temp < Tpf(1) ) Uk = exp(elem%pf(1,j))
        !   if (temp > Tpf(Npf)) Uk = exp(elem%pf(Npf,j))

          !linear extrapolation
          ! TO DO: in linear_1D_sorted
          if (temp < Tpf(1)) then
            ! write(*,*) "xi=",temp, "x1=",Tpf(1), "x2=",Tpf(2), "y1=",elem%pf(1,j), "y2=",elem%pf(2,j)
            Uk = exp( elem%pf(2,j) + (temp - Tpf(2))/(Tpf(1)-Tpf(2)) * (elem%pf(1,j)-elem%pf(2,j)) )
          elseif (temp > Tpf(Npf)) then
            ! write(*,*) "xi=",temp, "x1=",Tpf(Npf-1), "x2=",Tpf(Npf), "y1=",elem%pf(Npf-1,j), "y2=",elem%pf(Npf,j)
            Uk = exp( elem%pf(Npf-1,j) + (temp - Tpf(Npf-1))/(Tpf(Npf)-Tpf(Npf-1)) * (elem%pf(Npf,j)-elem%pf(Npf-1,j)) )
          else ! in range
            Uk = exp(Uka(1))
          endif
    
        return
    end function get_pf


    function phi_jl(temp, Ujl, Uj1l, ionpot)
        ! -------------------------------------------------------------- !
        ! Hubeny & Mihalas (2014)
        ! "Theory of Stellar Atmospheres", p.  94, eq. 4.35
        !
        !
        ! ionisation state j of element l ! (j1=j+1)
        ! ionpot = ionpot of ion j of element l (ionpot0l(H) = 13.6 eV)
        ! ionpot in J
        !
        ! Njl = Nj1l * ne * phi_jl
        ! --> SahaEq: Nj1l = Njl / (phi_jl * ne)
        ! --> NH- = NH * ne * phi_jl
        ! NH = NH+ * ne * phi_jl
        ! -------------------------------------------------------------- !
        real(kind=dp), intent(in) :: temp
        real(kind=dp) :: phi_jl
        real(kind=dp) ::  C1, expo
        real(kind=dp), intent(in) :: Ujl, Uj1l, ionpot
        C1 = 0.5*(HP**2 / (2.*PI*mel*kb))**1.5
        !C1 = (HPLANCK**2 / (2.*PI*M_ELECTRON*KBOLTZMANN))**1.5
    
        !!ionpot = ionpot * 100.*HPLANCK*CLIGHT !cm1->J
        !! -> avoiding dividing by big numbers causing overflow.
        !!maximum argument is 600, exp(600) = 3.77e260
        expo = exp(min(600d0,ionpot/(kb*temp)))
        if (ionpot/(kb*temp) >= 600d0) expo = huge_dp
    
        !if exp(300) it means phi is "infinity", exp(300) == 2e130 so that's enough
        phi_jl = C1 * (Ujl / Uj1l) * expo / (temp**1.5 + tiny_dp)
        if (phi_jl < phi_min_limit) phi_jl = 0d0 !tiny_dp ! or phi = 0d0 should be more correct ?
        ! but test to not divide by 0
    
        if (is_nan_infinity(phi_jl)>0) then
           write(*,*) "error, phi_jl", phi_jl, temp
           stop
        endif
    
        return
    end function phi_jl

    !Inverse Saha
    function SahaEq(temp, NI, UI1, UI, chi, ne) result(NI1)
        ! -------------------------------------------------------------- !
        ! Hubeny & Mihalas (2014)
        ! "Theory of Stellar Atmospheres" eq. 9.2
        ! (here, the inverse, i.e., NI+1/NI)
        !
        ! NI1 = NI * 2.*(h**2/(2pimkT))^-1.5 * &
        !                   UI1_UI * np.exp(-chiI/kT) / ne
        ! NI1 = NI/(ne*phijl)
        ! where NI, NI1 are the total number of particles of an element
        ! in the stage I, I+1
        ! -------------------------------------------------------------- !
        real(kind=dp), intent(in) :: temp
        real(kind=dp), intent(in) :: NI, UI1, UI, chi, ne
        real(kind=dp) :: phi, NI1
        phi = phi_jl(temp, UI, UI1, chi) ! chi in J
    
        !phi should be between phi_min_limit and exp(600)
        !and ne between ne_limit and nemax, so never nan nor infinity
        !further in General, a larg phi implies a ne close to 0, so the product remains small
    
        if (ne > 0.0) then
           NI1 = NI/(phi*ne)
        else !all in ground state, phi->infinity if T->0 and phi->0 if T->infinty
           !AND IF	ne -> 0, T likely goes to 0, all in neutral states
           !meaning if ne->0 phi goes to infinity
           NI1 = 0.0
        endif
    
        RETURN
    end function SahaEq

    !to do...

    !subroutine H_pf()
    !end subroutine H_pf

    !subroutine read_pf(type="....")
    !end subroutine

end module