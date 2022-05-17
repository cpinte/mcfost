module io_atom

   use atom_type
   use elements_type
   !use zeeman, only              : Lande_eff, ZeemanMultiplet
   use constantes
   use uplow
   use abo, only             : get_Barklem_cross_data
   use collision_atom, only      : read_collisions
   use messages
   use mcfost_env, only          : mcfost_utils
   use parametres, only          : art_hv
   use grid, only                : T, vturb, B_char
   use fits_utils
   use wavelengths_gas, only : compute_line_bound

   implicit none


   character(len=*), parameter :: path_to_atoms = "/Atoms/"
   real(kind=dp), parameter :: C1 = deux_pi * (electron_charge/EPSILON_0) * (electron_charge/mel / C_light)

   contains

   subroutine read_model_atom(atomunit, atom)
   !
   ! read independent atomic model
   ! Iinitialize the atom values
   !
      integer, intent(in) :: atomunit
      type (AtomType), intent(inout), target :: atom
      real(kind=dp) :: dummy
      character(len=2) :: IDread
      character(len=512) :: inputline, FormatLine
      logical :: find_abund, res, setup_common_gauss_prof
      real, allocatable, dimension(:) :: levelNumber
      logical, dimension(:), allocatable :: determined, parse_labs
      character(len=20) :: shapeChar, vdWChar, nuDepChar
      real(kind=dp) :: f, lambdaji, lambdamin
      integer :: Nread, i,j, EOF, nll, nc, kr, la

      
      EOF = 0
      res = .false.
      setup_common_gauss_prof = .false.



      open(unit=atomunit,file=trim(mcfost_utils)//trim(path_to_atoms)//trim(atom%filename),status="old")
      !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
      write(FormatLine,'("(1"A,I3")")') "A", 512

      !Read ID and fill atom abundance and weight
      CALL read_line(atomunit, FormatLine, inputline, Nread)
      read(inputline,*) IDread

      IDread(2:2) = to_lower(IDread(2:2))
      atom%ID = IDread
      write(*,*) "Reading atomic model of atom ", atom%ID
      find_abund = .false.
      do nll=1,Nelem
         if (Elems(nll)%ID == atom%ID) then
            atom%massf = Elems(nll)%massf
            write(*,*) "Abundance of atom ",atom%ID,": A =",Elems(nll)%Abund
            if (atom%ID=="H" .or. atom%ID=="He") then
               write(*,*) " -> mass fraction (%) = ", 100.*real(atom%massf)
            else
               write(*,*) " -> mass fraction (m/m(Fe) %) = ", 100.*real(atom%massf/Elems(26)%massf)
            endif
            if (Elems(nll)%abundance_set) then
               atom%periodic_table = nll
               atom%Abund = Elems(nll)%Abund
               atom%weight = Elems(nll)%weight
               atom%Rydberg = rydberg_atom(atom)
               find_abund=.true.
            end if
            exit
         end if
      end do

      if (.not.find_abund) then
         write(*,*) "Error no abundance found for atom ", atom%filename, " ", atom%ID
         stop
      end if

      !read Nlevel, Nline, Ncont
      call read_line(atomunit, FormatLine, inputline, Nread)
      read(inputline,*) atom%Nlevel, atom%Nline, atom%Ncont
      write(*,*) "Nlevel=",atom%Nlevel," Nline=",atom%Nline,&
                  " Ncont=", atom%Ncont

      !atomic level spectroscopic term
      allocate(atom%label(atom%Nlevel))
      !atomic level Energy w/ ground level
      allocate(atom%E(atom%Nlevel))
      !atomic level statistical weight
      allocate(atom%g(atom%Nlevel))
      !atomic level ionisation stage (i.e., 0 for neutral, 1 for singly ionised ...)
      allocate(atom%stage(atom%Nlevel))
      !atomic level id w/ other levels, RELATIVE TO that model.
      !index id starts at 0 for first level of the model.
      allocate(levelNumber(atom%Nlevel))
      !atomic level orbital quantum number
      allocate(atom%Lorbit(atom%Nlevel))
      !atomic level Spin (multiplicity)
      allocate(atom%qS(atom%Nlevel))
      !atomic level total angular momentum quantum number
      allocate(atom%qJ(atom%Nlevel))
      !is the level well identified ?
      allocate(determined(atom%Nlevel))
      !same as above
      allocate(parse_labs(atom%Nlevel))
      !atomic level's population at LTE (n^*)
      allocate(atom%nstar(atom%Nlevel,n_cells))
      !total number of transitions for this model
      atom%Ntr = atom%Nline + atom%Ncont
      !-> to remove tab_trans ??
      ! allocate(atom%tab_trans(atom%Ntr))
      allocate(atom%ij_to_trans(atom%Nlevel,atom%Nlevel))
      atom%ij_to_trans(:,:) = -99
      allocate(atom%i_trans(atom%Ntr),atom%j_trans(atom%Ntr))
      atom%Ntr_line = atom%Nline
      atom%nstar(:,:) = 0d0
      !default
      atom%qJ(:) = -99.
      atom%qS(:) = -99.0
      atom%Lorbit(:) = -99
      determined(:) = .false.
      parse_labs(:) = .false.

      !Start reading energy levels
      do i=1,atom%Nlevel
         call read_line(atomunit, FormatLine, inputline, Nread)
         read(inputline,*) atom%E(i), atom%g(i), atom%label(i), &
                           atom%stage(i), levelNumber(i)
         atom%E(i) = atom%E(i) * HP*C_LIGHT / (CM_TO_M)
      end do 

      ! Check if there is at least one continuum transition
      if (atom%stage(atom%Nlevel) /= atom%stage(atom%Nlevel-1)+1) then
         write(*,*) atom%stage
         write(*,*) atom%stage(atom%Nlevel), atom%stage(atom%Nlevel-1)+1
         write(*,*) "Atomic model does not have an overlying continuum"
         write(*,*) "exiting..."
         stop
      end if

      !Starting from now, i is the index of the lower level
      !it cannot be used as a loop index :-)

      !read bounb-bound (line) transitions
      write(*,*) "MAGNETIC POLARIZATION REMOVED AT THE MOMENT"
      allocate(atom%lines(atom%Nline))
      do kr=1,atom%Nline

         ! atom%tab_trans(kr) = kr
         atom%lines(kr)%lcontrib = .true. !init, all transitions contributes to opacity
                                          !in images, lines can be removed by setting lcontrib to .false. (expect if overlap)

         atom%lines(kr)%atom => atom
         atom%lines(kr)%polarizable = .false.
         atom%lines(kr)%g_lande_eff = -99.0
         atom%lines(kr)%glande_i = -99.0; atom%lines(kr)%glande_j = -99.0

         call read_line(atomunit, FormatLine, inputline, Nread)                       
         Nread = len(trim(inputline))
         !!write(*,*) Nread, inputline(1:Nread)
         read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%qwing, vdWChar,&
               atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
               atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
               atom%lines(kr)%Grad, atom%lines(kr)%cStark
 
         !indexes in the atomic model are C indexes!
         i = i + 1
         j = j + 1

         atom%lines(kr)%i = min(i,j)
         atom%lines(kr)%j = max(i,j)
         atom%ij_to_trans(atom%lines(kr)%i,atom%lines(kr)%j) = kr
         atom%i_trans(kr) = atom%lines(kr)%i
         atom%j_trans(kr) = atom%lines(kr)%j

         if (atom%lines(kr)%qwing < 3.0) then

            call Warning("qwing for line lower than 3! setting to 3")
            atom%lines(kr)%qwing = 3.0

         endif

         !because levels correspond to different transitions, 
         !we need to test if the level has already been indentified
         if (.not.parse_labs(atom%lines(kr)%i)) then
            call parse_label(atom%label(atom%lines(kr)%i),&
               atom%g(atom%lines(kr)%i),&
               atom%qS(atom%lines(kr)%i),&
               atom%Lorbit(atom%lines(kr)%i),&
               atom%qJ(atom%lines(kr)%i), &
               determined(atom%lines(kr)%i))
            parse_labs(atom%lines(kr)%i) = .true.
         end if
         if (.not.parse_labs(atom%lines(kr)%j)) then
            call parse_label(atom%label(atom%lines(kr)%j),&
               atom%g(atom%lines(kr)%j),&
               atom%qS(atom%lines(kr)%j),&
               atom%Lorbit(atom%lines(kr)%j),&
               atom%qJ(atom%lines(kr)%j), &
               determined(atom%lines(kr)%j))
            parse_labs(atom%lines(kr)%i) = .true. !even if determined is false.
         end if

         ! if (lmagnetized) then
         !    if ((abs(atom%qJ(atom%lines(kr)%i) - atom%qJ(atom%lines(kr)%j)) <= 1.) .and. &
         !       (atom%lines(kr)%g_Lande_eff <= -99) .and. &
         !       (determined(atom%lines(kr)%j)) .and. (determined(atom%lines(kr)%i))) then !
         !       ! do not compute geff if term is not
         !       ! determined or if the geff is read from file
         !       ! ie if g_lande_eff > -99
         !       ! fill lande_g_factor if determined and not given
         !       ! for b-b transitions
         !       call Lande_eff(atom, kr)
         !    end if
         !    atom%lines(kr)%polarizable = (atom%lines(kr)%g_lande_eff > -99).and.(atom%lines(kr)%ZeemanPattern /= 0)
         !    if (atom%lines(kr)%polarizable) then
         !       call ZeemanMultiplet(atom%lines(kr))
         !    endif
         ! endif !lmagnetized
         ! oscillator strength saved
         atom%lines(kr)%fosc = f
         lambdaji = (HP * C_LIGHT) / (atom%E(j) - atom%E(i))
         atom%lines(kr)%Aji = C1 / (lambdaji**2) * (atom%g(i) / atom%g(j)) * f
         atom%lines(kr)%Bji = (lambdaji**3) / (2.0 * HP * C_LIGHT) *atom%lines(kr)%Aji
         atom%lines(kr)%Bij = (atom%g(j) / atom%g(i)) * atom%lines(kr)%Bji
         atom%lines(kr)%lambda0 = lambdaji / NM_TO_M

         !gi * Bij = gj * Bji
         !gi/gj = Bji/Bij so that niBij - njBji = Bij (ni - gij nj)
         atom%lines(kr)%gij = atom%lines(kr)%Bji / atom%lines(kr)%Bij !gi/gj
         atom%lines(kr)%twohnu3_c2 = atom%lines(kr)%Aji / atom%lines(kr)%Bji

         !profile function for lines
         select case(shapechar)
            case ("GAUSS")
               atom%lines(kr)%Voigt = .false.
               if (.not.setup_common_gauss_prof) then
                  setup_common_gauss_prof = .true.
                  atom%lgauss_prof = .true. !set only once per atom!
               endif
            case ("VOIGT")
               atom%lines(kr)%Voigt = .true.
            ! case ("THOMSON")
            !    atom%lines(kr)%Voigt = .false.
            !    atom%lines(kr)%pVoigt = .true.
            case default
               write(*,*) "Line profile shape", shapechar, " unknown, using voigt."
               atom%lines(kr)%Voigt = .true.
         end select


         !Van der Waaks collision method
         !line%a allocated in opacity_atom.f90
         atom%lines(kr)%cvdWaals(4) = 0.
         atom%lines(kr)%cvdWaals(2) = 0.
         atom%lines(kr)%cvdWaals(1) = 0.
         atom%lines(kr)%cvdWaals(3) = 0.
         select case (vdwchar)
            case ("BARKLEM")
               atom%lines(kr)%vdWaals = "BARKLEM"
               call get_Barklem_cross_data(atom, kr, res)
               if (.not. res) then
                  write(*,*) &
                       "Line <atom%lines(kr)%j>->atom%lines(kr)%i>", &
                       " cannot be treated with Barklem type", &
                       " broadening."
                  write(*,*) "using UNSOLD"
                  atom%lines(kr)%vdWaals = "UNSOLD"
                  atom%lines(kr)%cvdWaals(4) = 0.
                  atom%lines(kr)%cvdWaals(2) = 0.
               end if
            case default
               atom%lines(kr)%vdWaals = "UNSOLD"
         end select

         !define the extension of a line for non-LTE loop and images.
         !here, limage is .false. vmax is computed from the atomic file
         !and from the "natural" width of the line.
         call compute_line_bound(atom%lines(kr),.false.)
         !if limage is .true., vmax is atom%vmax_rt.
         !with limage, line%Nlambda is also fixed by atom%n_speed_rt.


      end do !end loop over bound-bound transitions

      ! ----------------------------------------- !
      ! starts reading bound-free transitions
      ! cross-sections allocated once the final
      ! wavelength grid is known.
      ! ----------------------------------------- !
      allocate(atom%continua(atom%Ncont))
      do kr=1,atom%Ncont

         ! atom%tab_trans(kr+atom%Ntr_line) = kr+atom%Ntr_line
         atom%continua(kr)%lcontrib = .true.

         atom%continua(kr)%atom => atom

         call read_line(atomunit,FormatLine, inputline, Nread)
         read(inputline, *) j, i, atom%continua(kr)%alpha0,&
            atom%continua(kr)%Nlambda, nuDepChar, atom%continua(kr)%lambdamin
         j = j + 1
         i = i + 1
         atom%continua(kr)%j = max(i,j)
         atom%continua(kr)%i = min(i,j)
         atom%ij_to_trans(atom%continua(kr)%i,atom%continua(kr)%j) = kr+atom%Ntr_line
         atom%i_trans(kr+atom%Ntr_line) = atom%continua(kr)%i
         atom%j_trans(kr+atom%Ntr_line) = atom%continua(kr)%j

         lambdaji = (HP*C_LIGHT)/ (atom%E(j)-atom%E(i))
         atom%continua(kr)%lambda0 = lambdaji/NM_TO_M !nm
         atom%continua(kr)%lambdamax = atom%continua(kr)%lambda0

         !no need to allocate cont%lambda and cont%alpha.
         !only cont%alpha_file and cont%lambda_file are needed if present in the file.
         select case (nudepchar)
            case ("EXPLICIT")
               allocate(atom%continua(kr)%alpha_file(atom%continua(kr)%Nlambda))
               allocate(atom%continua(kr)%lambda_file(atom%continua(kr)%Nlambda))
               atom%continua(kr)%hydrogenic=.false.
               !written from increasing frequencies so decreasing wavelengths
               do la=atom%continua(kr)%Nlambda,1,-1
                  call read_line(atomunit,FormatLine, inputline, Nread)
                  read(inputline,*) atom%continua(kr)%lambda_file(la), &
                       atom%continua(kr)%alpha_file(la)
               end do

               do la=2,atom%continua(kr)%Nlambda
                  if (atom%continua(kr)%lambda_file(la) < &
                       atom%continua(kr)%lambda_file(la-1)) then
                     write(*,*) "continuum wavelength not monotonous"
                     write(*,*) "exiting..."
                     stop
                  end if
               end do
               !not extrapolated if explicit ?
               !should consider neglecting occupation probability for that case
               atom%continua(kr)%lambdamax = maxval(atom%continua(kr)%lambda_file)
            case ("HYDROGENIC")
               atom%continua(kr)%hydrogenic=.true.
               !!tmp
               !        atom%continua(kr)%lambdamin = 5.0_dp
               !        atom%continua(kr)%lambdamin = 0.05 * atom%continua(kr)%lambda0
               !        atom%continua(kr)%lambdamin = max(10.0_dp, 1d-2 * atom%continua(kr)%lambda0)
               !!tmp
               write(*,'(" Continuum "(1I3)" -> "(1I3)" at "(1F12.5)" nm")') &
                  atom%continua(kr)%i, atom%continua(kr)%j, atom%continua(kr)%lambda0
               write(*,'(" -> lower edge cut at "(1F12.5)" nm !")'), atom%continua(kr)%lambdamin
     
               if (atom%continua(kr)%lambdamin>=atom%continua(kr)%lambda0) then
                  write(*,*) "Minimum wavelength for continuum is larger than continuum edge."
                  write(*,*) kr, atom%continua(kr)%lambda0, atom%continua(kr)%lambdamin
                  write(*,*) "exiting..."
                  stop
               end if

               if (ldissolve) then
                  if (atom%ID=="H") then
                     call search_cont_lambdamax (atom%continua(kr), atom%Rydberg, atom%stage(i)+1,atom%E(j),atom%E(i))
                  endif
               endif
               !the Nlambda from file is not used except if not hydrogenic
               atom%continua(kr)%Nlambda = 0
            case default

               call error("nudepchar!")
            
         end select
      enddo !bound-free

      atom%set_ltepops = .true. !by default compute lte populations
      !non-LTE pops in electronic density ? write non-LTE pops to file ? 
      atom%NLTEpops = .false. !set to true during non-LTE loop if initial /= 0

      ! allocate some space
      if (atom%initial==2) then
         if (.not.atom%active) then
            write(*,*) atom%ID, " is passive! cannot use ZERO_RADIATION solution, set to LTE."
            atom%initial=0
         endif
      endif

      if (atom%active) then

         atom%cswitch = 1.0_dp

         ! reading collision rates of RH
         if (atom%ID /= "H") then
            write(*,*) "  -> Reading collision data from RH for atom ", atom%ID
            call read_collisions(atomunit, atom)
            call warning("  ** the data for each transition must be consistent with transitions in the file")
            write(*,*) "a value of -99 is given to missing transitions and are skipped **"
         endif

         allocate(atom%n(atom%Nlevel,n_cells))
         atom%n = 0.0_dp
         if (atom%initial == 0) then!(initial==0 .or. initial==3)
            atom%n = atom%nstar !still need to be computed
            atom%set_ltepops = .true.
            !if (initial==0) then
            write(*,*) " -> Setting initial solution to LTE "
            !else
            !  cswitch
            !endif


         else if (atom%initial == 1) then

            write(*,*) " -> Reading (non-LTE AND LTE) populations from file..."
            call read_pops_atom(atom)
            atom%NLTEpops = .true.
            atom%set_ltepops = .false. !read and USE also LTE populations from file!!

         else if (atom%initial == 2) then

            call error("initial solution == 2 (opt thin) not implemented yet!")
            atom%n = atom%nstar
            atom%set_ltepops = .true.
            !nlte pops is false and are set after first electron density and lte pops.

         else if (atom%initial == 3) then
            atom%n = atom%nstar !still need to be computed
            !like initial==0
            atom%set_ltepops = .true.
            if (.not. lforce_lte) then !otherwise cswitch is not LTE
               write(*,*) " -> Setting initial solution to LTE with CSWITCH "
               atom%cswitch = cswitch_val
               if (.not. lcswitch_enabled) lcswitch_enabled = .true.!we need at least one
            endif

         else if (atom%initial == 4) then

            call error("initial solution == 4 (Sobolev) not implemented yet!")
            atom%n = atom%nstar
            atom%set_ltepops = .true.
            !nlte pops is false and are set after first electron density and lte pops.

         end if
      else !not active = PASSIVE
         !other initials do not matter, always LTE ! (except if read)
         if (atom%initial == 1) then

            allocate(atom%n(atom%Nlevel,n_cells)) !not allocated if passive, n->nstar
            write(*,*) " -> Reading (non-LTE AND LTE) populations from file for passive atom..."
            call read_pops_atom(atom)
            !by default. Need a case with only LTE ?
            atom%NLTEpops = .true.
            atom%set_ltepops = .false. !read and USE also LTE populations from file!!

          !atom%NLTEpops = .false. still false at this point as we need pops to do electron densities
         else !pure passive without nlte pops from previous run
          !        atom%NLTEpops=.false.  !-> default values
          !        atom%set_ltepops = .true.
            atom%n => atom%nstar !initialised to zero
          !atom%n is an alias for nstar in this case
         end if
      end if !end is active

      deallocate(levelNumber)
      deallocate(determined)
      deallocate(parse_labs)
      !close atomic file
      close(unit=atomunit) !later it will be open again for

      return
   end subroutine read_Model_Atom

   subroutine read_atomic_models()
      integer :: EOF=0,Nread, nmet, mmet, nblancks, nact, npass
      integer :: kr, k, imax, ic
      real(kind=dp) :: min_resol, max_resol
      integer :: unit = 1
      character(len=15) :: FormatLine
      character(len=512) :: popsfile, filename, inputline
      character(len=2) :: IDread
      type (AtomType), pointer :: atom


      write(FormatLine,'("(1"A,I3")")') "A", 512
      if (N_atoms > 1) then
         write(*,*) "Reading ", N_atoms, " atoms"
      else 
         write(*,*) "Reading ", N_atoms, " atom"
      endif

      !go through file and read the different atomic models
      do nmet = 1, N_atoms

         open(unit=unit+nmet,file=trim(mcfost_utils)//TRIM(path_to_atoms)//trim(atoms(nmet)%p%filename),status="old")
         call read_line(unit+nmet, FormatLine, inputline, Nread)
         read(inputline,*) IDread
         if (nmet==1 .and. IDread/="H ") then
            write(*,*) "first atom is ", IDread
            write(*,*) "First atomic model read has to be Hydrogen"
            write(*,*) "Exting..."
            stop
         end if
         close(unit+nmet)
         ! check duplicates
         IDread(2:2) = to_lower(IDread(2:2))
         do mmet = 1,nmet-1 !compare the actual (IDread) with previous
         !write(*,*) mmet, nmet
            if (Atoms(mmet)%p%ID == IDread) then
               write(*,*) "Already read a model for this atom ", IDread
               write(*,*) "exiting..."
               stop
            end if
         end do

         !sets atom%ID too.
         call read_model_atom(unit+nmet, Atoms(nmet)%p)
         ! if (Atoms(nmet)%p%lgauss_prof) then
         !    write(*,*) "**-> allocating gauss profile for atom", Atoms(nmet)%p%id
         !    call alloc_common_gauss_profile(Atoms(nmet)%p) !deallocated in opacity
         ! endif

      end do !over atoms

      !Temporary
      !check that if helium is active and electron are iterated H must be active at the moment !!
      check_helium : do nmet=1, N_atoms
         if (atoms(nmet)%p%id=="He") then

            if ((n_iterate_ne > 0).and.(atoms(nmet)%p%active)) then

               if (atoms(nmet)%p%stage(1) > 0) then
                  write(*,*) " !!!!!!!!! "
                  call warning("Helium ground state is not in neutral stage ! Must be He I")
                  write(*,*) atoms(nmet)%p%stage(1), atoms(nmet)%p%label(1), &
                     atoms(nmet)%p%E(1), atoms(nmet)%p%g(1)
                  write(*,*) " !!!!!!!!! "
                  stop
               endif

               !H always present and the first.
               !force to be active if helium active ?
               !at the moment print a warning!
               if (.not.atoms(1)%p%active) then
                  write(*,*) " !!!!!!!!!!!!!!!!!!!!!!! "
                  call WARNING(" Hydrogen is passive, while helium is active and n_iterate_ne > 0!!")
                  write(*,*) " !!!!!!!!!!!!!!!!!!!!!!! "
                  exit check_helium
               endif

            endif

         endif
      enddo check_helium


      ! Alias to the most importent one
      !always exits !!
      Hydrogen=>Atoms(1)%p
      if (.not.associated(Hydrogen, Atoms(1)%p)) CALL Error(" Hydrogen alias not associated to atomic model!")

      ! Aliases to active atoms
      nact = 0; npass = 0
      if (Nactiveatoms > 0) then
         allocate(ActiveAtoms(Nactiveatoms))
      end if
      if (n_atoms - Nactiveatoms == Npassiveatoms) then
         allocate(PassiveAtoms(Npassiveatoms))
      else
         write(*,*) "Error, number of passive atoms is not N Nactive"
         stop
      end if

      ! keep a duplicate in Elements
      write(*,*) "order#     ID   periodic-table#    ACTIVE    #lines   #continua"
      do nmet=1,n_atoms
         write(*,*) nmet, Atoms(nmet)%p%ID, &
            Atoms(nmet)%p%periodic_table, Atoms(nmet)%p%active, &
            Atoms(nmet)%p%Nline, Atoms(nmet)%p%Ncont

         Elems(Atoms(nmet)%p%periodic_table)%nm = nmet

         !Check Nstage
         if (Elems(Atoms(nmet)%p%periodic_table)%Nstage < maxval(Atoms(nmet)%p%stage) + 1) then
            write(*,*) Atoms(nmet)%p%id, maxval(Atoms(nmet)%p%stage) + 1
            write(*,*) "Ns pf = ", Elems(Atoms(nmet)%p%periodic_table)%Nstage
            call error("Model has more ionisation stages than the one in the partition function!")
         endif

         if (Atoms(nmet)%p%periodic_table==2)  then
            Helium => Atoms(nmet)%p
            helium_is_active = helium%active
            write(*,*) "Helium pointers associated to atom in the table", nmet
            if (helium_is_active) then
               write(*,*) " And it is active !", Atoms(nmet)%p%active
            endif
            if (.not.associated(Helium,Atoms(nmet)%p)) &
               call Warning(" Helium alias is not associated to an atomic model!")
         end if

         if (allocated(ActiveAtoms)) then
            if (Atoms(nmet)%p%active) then
               nact = nact + 1 !got the next index of active atoms
               ActiveAtoms(nact)%p => Atoms(nmet)%p
               Atoms(nmet)%p%activeindex = nact
            end if
         end if

         if (allocated(PassiveAtoms)) then
            if (.not.Atoms(nmet)%p%active) then
               npass = npass+1
               PassiveAtoms(npass)%p => Atoms(nmet)%p
            end if
         end if
      end do

      return
   end subroutine read_Atomic_Models


   subroutine search_cont_lambdamax (cont, Rinf, Z, Ej, Ei)
      !Search the lambdamax = first lambda for which Gaunt < 0
      !Because I do not have an extrapolation routine, there is no need to extrapolate beyond
      !lambdamax if gaunt is negative.
      !
      ! If extrapolation, we need to find where alpha * wocc goes to zero
      !  approximately
      !

      real(kind=dp), intent(in) :: Ej, Ei, Rinf
      integer, intent(in) :: Z
      type (AtomicContinuum), intent(inout) :: cont
      real, parameter :: l1 = 1e3, dlam = 1.0 !nm !0.01
      real(kind=dp) :: n_eff,  u ! = n_Eff**2 * eps = hnu/Z/Z/E_RYDBERG - 1
      real(kind=dp) :: Gaunt_bf, x

      x = l1 * cont%lambda0
      cont%lambdamax = cont%lambda0
      n_eff = Z * sqrt(Rinf / (Ej-Ei))
      Gaunt_bf = 1.0_dp
      do while(cont%lambdamax < x)

         u = n_eff**2 * HP*C_LIGHT / (NM_TO_M * cont%lambdamax) / Z*Z / E_RYDBERG - 1

         Gaunt_bf = 1d0 + 0.1728 * (n_eff**(-2./3.)) * (u+1d0)**(-2./3.) * (u-1d0) &
            - 0.0496*(n_eff**(-4./3.)) * (u+1d0)**(-4./3.) * (u*u + 4./3. * u + 1d0)

         if (Gaunt_bf < 0) exit

         cont%lambdamax = cont%lambdamax + dlam

      enddo

      cont%lambdamax = min(cont%lambdamax,1d5)
      write(*,*) cont%atom%ID, " cont, n=",real(n_eff), cont%j, cont%i, real(cont%lambda0),"nm"
      write(*,*) " -> pseudo-continuum: ", real(cont%lambdamax), " nm: frac=", real(cont%lambdamax/cont%lambda0)

      return
   end subroutine search_cont_lambdamax


   subroutine write_pops_atom(atom,iter,step)
   ! ----------------------------------------------------------------- !
   ! write Atom populations
   ! First, NLTE populations, then LTE populations.
   ! ----------------------------------------------------------------- !
   type (AtomType), intent(in) :: atom
   integer,  optional :: iter, step
   integer :: unit, blocksize, naxes(5), naxis,group, bitpix, fpixel
   logical :: extend, simple, lte_only
   integer :: nelements, hdutype, status, k, sys_status
   character(len=512) :: cmd, popsF
   character(len=10000) :: step_c, iter_c

   !lte_only = .not.atom%active
   lte_only = (.not.atom%NLTEpops) !can be true if atom ACTIVE or INITIAL_SOLUTION==OLD_POPULATIONS
   !because maybe NLTE pops are read but atom is considered as PASSIVE (i.e., no nlte iterations)

   status = 0
   !get unique unit number
   call ftgiou(unit,status)
   blocksize=1
   simple = .true. !Standard fits
   group = 1 !??
   fpixel = 1
   extend = .true.
   bitpix = -64

   popsF = trim(atom%ID)//".fits.gz"


   if (present(iter)) then
      !suppose that the file does not exist
      write(iter_c, '(i0)') iter
      if (.not.present(step)) step=0
      write(step_c, '(i0)') step
      !folder iterations must exist
      popsF = "iterations/"//trim(atom%ID)//'_iter'//trim(iter_c)//'_step'//trim(step_c)//'.fits.gz'

   else
      !check if data file already exist, can be the case if initial solutions is OLD_POPULATIONS
      cmd = "ls "//popsF!trim(atom%ID)//".fits.gz"
      call appel_syst(cmd, sys_status)
      if (sys_status == 0) then !means the file exist

         cmd = "mv "//trim(atom%ID)//".fits.gz"//" "//trim(atom%ID)//"_oldpop.fits.gz"
         call appel_syst(cmd, sys_status)
         if (sys_status /= 0) then
            call error("Error in copying old pops!")
         endif

      endif
   endif

   call ftinit(unit, trim(root_dir)//"/"//TRIM(popsF), blocksize, status)
   if (status > 0) then
      write(*,*) "Cannot create fits file ", popsF
      call print_error(status)
   endif
   naxes(1) = atom%Nlevel

   if (lVoronoi) then
      naxis = 2
      naxes(2) = n_cells
      nelements = naxes(1)*naxes(2)
   else
      if (l3D) then
         naxis = 4
         naxes(2) = n_rad
         naxes(3) = 2*nz
         naxes(4) = n_az
         nelements = naxes(1) * naxes(2) * naxes(3) * naxes(4)
      else
         naxis = 3!4
         naxes(2) = n_rad
         naxes(3) = nz
         ! 			naxes(4) = 1 !write naz which is one
         nelements = naxes(1) * naxes(2) * naxes(3)
      end if
   end if


   if (lte_only) then

      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
      if (status > 0) then
         write(*,*) "Error writing LTE pops to file"
         call print_error(status)
      endif

      ! Additional optional keywords
      call ftpkys(unit, "UNIT", "m^-3", "(LTE)", status)
      if (status > 0) then
         write(*,*) "Error writing LTE pops to file (2)"
         call print_error(status)
      endif
      !write data
      call ftpprd(unit,group,fpixel,nelements,atom%nstar,status)
      if (status > 0) then
         write(*,*) "Error writing LTE pops to file (3)"
         call print_error(status)
      endif
   else !NLTE + LTE
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
      if (status > 0) then
         write(*,*) "Error writing nLTE pops to file"
         call print_error(status)
      endif
      ! Additional optional keywords
      call ftpkys(unit, "UNIT", "m^-3", "(NLTE) ", status)
      if (status > 0) then
         write(*,*) "Error writing nLTE pops to file (2)"
         call print_error(status)
      endif

      call ftpprd(unit,group,fpixel,nelements,atom%n,status)
      if (status > 0) then
         write(*,*) "Error writing nLTE pops to file (3)"
         call print_error(status)
      endif

      call ftcrhd(unit, status)
      if (status > 0) then
         write(*,*) "Error writing nLTE pops to file (4)"
         call print_error(status)
      endif

      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)
      if (status > 0) then
         write(*,*) "Error writing LTE pops to file (5)"
         call print_error(status)
      endif
      ! Additional optional keywords
      call ftpkys(unit, "UNIT", "m^-3", "(LTE)",status)
      if (status > 0) then
         write(*,*) "Error writing LTE pops to file (6)"
         call print_error(status)
      endif
      !write data
      call ftpprd(unit,group,fpixel,nelements,atom%nstar,status)
      if (status > 0) then
         write(*,*) "Error writing LTE pops to file (7)"
         call print_error(status)
      endif

   end if


   call ftclos(unit,status) !close
   call ftfiou(unit,status) !free

   if (status > 0) call print_error(status)

   return
   end subroutine write_pops_atom

   subroutine read_pops_atom(atom)
   ! ---------------------------------------------------------- !
   ! read Atom populations.
   ! First, NLTE populations, then LTE populations.
   !
   ! To DO:
   ! If the read populations are on a different grid than the
   ! actual grid, they are interpolated.
   ! However: there is not extrapolation at the moment.
   ! ---------------------------------------------------------- !
   type (AtomType), intent(inout) :: atom
   integer :: unit, status, blocksize, naxis,group, bitpix, fpixel
   logical :: extend, simple, anynull, show_warning = .true.
   integer :: nelements, naxis2(4), Nl, naxis_found, hdutype, l, icell
   character(len=256) :: some_comments, popsF

   status = 0
   call ftgiou(unit,status)

   popsF = TRIM(atom%ID)//".fits.gz"

   call ftopen(unit, trim(popsF), 0, blocksize, status)
   if (status > 0) then
      write(*,*) "Cannot open file with populations! ", atom%id
      call print_error(status)
   endif


   simple = .true. !Standard fits
   group = 1
   fpixel = 1
   extend = .true.
   bitpix = -64

   !NLTE
   call ftmahd(unit,1,hdutype,status)
   if (status > 0) then
      write(*,*) "Cannot read nlte populations ! ", atom%id
      call print_error(status)
      stop
   endif

   if (lVoronoi) then
      naxis = 2
      call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
      if (status > 0) then
         write(*,*) "error reading number of axis (naxis)"
         call print_error(status)
         stop
      endif

      !Nlevel
      call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading Nlevel from file (naxis1)"
         call print_error(status)
         stop
      endif

      if (naxis_found /= atom%Nlevel) then
         if (naxis_found /= naxis2(1)) then
            write(*,*) "Nlevel read does not match atom !", atom%Nlevel
            stop
         endif
      endif
      nelements = naxis_found

      !n_cells
      call ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading nrad from file (naxis2)"
         call print_error(status)
         stop
      endif

      nelements = nelements * naxis_found

   else

      if (l3D) then
         naxis = 4
      else
         !no naz
         naxis = 3!4
      end if

      !Number of axis
      call ftgknj(unit, 'NAXIS', 1, naxis, naxis2, naxis_found, status)
      !naxis2 = (Nlevel,nrad, naz, nphi) = (Naxis1,Naxis2,Naxis3,Naxis4)
      !naxis4 is 0 if not l3D
      if (status > 0) then
         write(*,*) "error reading number of axis (naxis)"
         call print_error(status)
         stop
      endif

      !Nlevel
      call ftgkyj(unit, "NAXIS1", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading Nlevel from file (naxis1)"
         call print_error(status)
         stop
      endif

      if (naxis_found /= atom%Nlevel) then
         if (naxis_found /= naxis2(1)) then
            write(*,*) "Nlevel read does not match atom !", atom%Nlevel
            stop
         endif
      endif
      nelements = naxis_found

      !R
      call ftgkyj(unit, "NAXIS2", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading nrad from file (naxis2)"
         call print_error(status)
         stop
      endif
      ! 		if (naxis_found /= nrad) then
      ! 			write(*,*) "nrad read does not match grid !", naxis_found
      ! 			stop
      ! 		endif
      nelements = nelements * naxis_found

      !z
      call ftgkyj(unit, "NAXIS3", naxis_found, some_comments, status)
      if (status > 0) then
         write(*,*) "error reading nz from file (naxis3)"
         call print_error(status)
         stop
      endif
      ! 		if (naxis_found /= nz) then
      ! 			write(*,*) "nz read does not match grid !", naxis_found
      ! 			stop
      ! 		endif
      nelements = nelements * naxis_found

      !phi axis ?
      if (l3D) then
         call ftgkyj(unit, "NAXIS4", naxis_found, some_comments, status)
         if (status > 0) then
            write(*,*) "error reading naz from file (naxis4)"
            call print_error(status)
            stop
         endif
         ! 			if (naxis_found /= n_az) then
         ! 				write(*,*) "n_az read does not match grid !", naxis_found
         ! 				stop
         ! 			endif
         nelements = nelements * naxis_found
      endif

   endif !lvoronoi

   if (nelements /= atom%Nlevel * n_cells) then
      write(*,*) " read_pops_atom does not do interpolation yet!"
      call Error (" Model read does not match simulation box")
   endif

   !READ NLTE
   call FTG2Dd(unit,1,-999,shape(atom%n),atom%Nlevel,n_cells,atom%n,anynull,status)
   if (status > 0) then
      write(*,*) "error reading non-LTE populations ", atom%id
      call print_error(status)
      stop
   endif

   !now LTE
   call FTMAHD(unit,2,hdutype,status)
   if (status > 0) then
      write(*,*) "error opening LTE hdu ", atom%id
      call print_error(status)
      stop
   endif

   atom%nstar(:,:) = 0.0_dp
   call FTG2Dd(unit,1,-999,shape(atom%nstar),atom%Nlevel,n_cells,atom%nstar,anynull,status)
   if (status > 0) then
      write(*,*) "error reading LTE populations ", atom%id
      call print_error(status)
      stop
   endif

   call ftclos(unit, status) !close
   if (status > 0) then
      write(*,*) "error cannot close file in ", trim(popsF)
      call print_error(status)
      stop
   endif

   call ftfiou(unit, status) !free
   if (status > 0) then
      write(*,*) "error cannot free file unit!"
      call print_error(status)
      ! 		stop
   endif

   ! write(*,*) " min/max pops for each level:"
   do l=1,atom%Nlevel
      write(*,"('Level #'(1I3))") l
      write(*,'("  -- min(n)="(1ES20.7E3)" m^-3; max(n)="(1ES20.7E3)" m^-3")') , &
           minval(atom%n(l,:),mask=(atom%n(l,:)>0)), maxval(atom%n(l,:))
      write(*,'("  -- min(nstar)="(1ES20.7E3)" m^-3; max(nstar)="(1ES20.7E3)" m^-3")')  &
           minval(atom%nstar(l,:),mask=(atom%nstar(l,:)>0)), maxval(atom%nstar(l,:))
   enddo

   return

   end subroutine read_pops_atom

end module io_atom