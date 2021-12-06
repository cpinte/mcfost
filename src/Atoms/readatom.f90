!
! Read ascii model atom in RH format (H. Uitenbroak).
! However, some values in the model atom are not used compared to RH.
!
! To Do: mcfost atomic model atom format
!        include in ref.para files.
!
module readatom

   use atom_type, only           : AtomicLine, AtomicContinuum, AtomType, Element, parse_label, &
                                 rydberg_atom, n_eff, find_continuum, atomZnumber, ATOM_ID_WIDTH
   use atmos_type, only          : Nelem, Hydrogen, Helium, Elements, T, ne, vturb, lmagnetized, &
                                 icompute_atomRT, Natom, NpassiveAtoms, NactiveAtoms, Atoms,     &
                                 PassiveAtoms, ActiveAtoms, helium_is_active
   use zeeman, only              : Lande_eff, ZeemanMultiplet
   use getlambda
   use constant
   use uplow
   use getline
   use barklem, only             : get_Barklem_cross_data
   use io_atomic_pops, only      : read_pops_atom
   use collision, only           : read_collisions
   use solvene, only             : Max_ionisation_stage, get_max_nstage
   !$ use omp_lib
   use messages
   use mcfost_env, only          : mcfost_utils
   use parametres, only          : art_hv

   implicit none

   character, parameter :: COMMENT_CHAR="#"
   character(len=*), parameter :: ATOMS_INPUT = "./atoms.input"!"/Atoms/atoms.input"
   character(len=*), parameter :: path_to_atoms = "/Atoms/"
   !real, parameter :: MAX_ABUND_ERROR=0.001

   integer, parameter :: cwitch_niter = 12!nint(ceil(log(cswitch_val)/log(cswitch_down_scaling_factor)))
   real(kind=dp), parameter :: cswitch_val = 1d12
   real(kind=dp), parameter :: cswitch_down_scaling_factor = 10.0!ceil(exp(log(cswitch_val)/cswitch_niter))
   logical :: cswitch_enabled = .false.

   !Until an extrapolation routine is found, the lambdamax is found where Gaunt(lambda) < 0
   real, parameter :: fact_pseudo_cont = 100.0 ! the hydrogenic bound-free will be extrapolated
   !beyond the edge limit (lambda0) up to lambda0 * fact_pseudo_cont
   !if it is 1.0, no extrapolation i.e., lambdamax=lambda0

   contains

   subroutine read_Model_Atom(atomunit, atom, atom_file)
   !
   ! read independent atomic model
   ! Iinitialize the atom values
   !
      integer, intent(in) :: atomunit
      integer :: kr, k, la, alloc_status
      type (AtomType), intent(inout), target :: atom
      character(len=*), intent(in) :: atom_file
      character(len=MAX_LENGTH) :: inputline, FormatLine
      integer :: Nread, i,j, EOF, nll, nc, Nfixed !deprecation future
      real, allocatable, dimension(:) :: levelNumber
      logical :: Debeye, match, res, setup_common_gauss_prof
      logical, dimension(:), allocatable :: determined, parse_labs
      real(kind=dp), dimension(:), allocatable :: old_nHtot
      character(len=20) :: shapeChar, symmChar, optionChar, vdWChar, nuDepChar
      character(len=2) :: IDread
      real(kind=dp) :: C1, vDoppler, f, lambdaji
      real(kind=dp) :: lambdamin
      real :: geff, gamma_j, gamma_i
      EOF = 0
      res = .false.
      setup_common_gauss_prof = .false.


      C1 = 2.*PI * (Q_ELECTRON/EPSILON_0) * (Q_ELECTRON/M_ELECTRON / CLIGHT)

      open(unit=atomunit,file=atom_file,status="old")
      !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
      write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH

      !Read ID and fill atom abundance and weight
      CALL getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
      read(inputline,*) IDread

      IDread(2:2) = to_lower(IDread(2:2))
      atom%ID = IDread
      write(*,*) "Reading atomic model of atom ", atom%ID
      match = .false.
      do nll=1,Nelem
         if (Elements(nll)%ptr_elem%ID.eq.atom%ID) then
            atom%massf = Elements(nll)%ptr_elem%massf
            write(*,*) "Abundance of atom ",atom%ID,": A =",Elements(nll)%ptr_elem%Abund
            if (atom%ID=="H" .or. atom%ID=="He") then
               write(*,*) " -> mass fraction (%) = ", 100.*real(atom%massf)
            else
               write(*,*) " -> mass fraction (m/m(Fe) %) = ", 100.*real(atom%massf/Elements(26)%ptr_elem%massf)
            endif
            if (Elements(nll)%ptr_elem%abundance_set) then
               atom%periodic_table=nll
               atom%Abund=Elements(nll)%ptr_elem%Abund
               atom%weight=Elements(nll)%ptr_elem%weight
               atom%Rydberg = rydberg_atom(atom)
               match=.true.
            end if
            exit
         end if
      end do

      if (.not.match) then
         write(*,*) "Error no abundance found for atom ", atom_file, " ", atom%ID
         stop
      end if

      !read Nlevel, Nline, Ncont and Nfixed transitions
      !fixed transitions are read for compatibility with rh.
      call getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
      read(inputline,*) atom%Nlevel, atom%Nline, atom%Ncont, Nfixed
      write(*,*) "Nlevel=",atom%Nlevel," Nline=",atom%Nline,&
                  " Ncont=", atom%Ncont!, " Nfixed=", atom%Nfixed

      atom%cswitch = 1.0_dp

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
      allocate(atom%at(atom%Ntr))
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
         call getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)
         read(inputline,*) atom%E(i), atom%g(i), atom%label(i), &
                           atom%stage(i), levelNumber(i)
         atom%E(i) = atom%E(i) * HPLANCK*CLIGHT / (CM_TO_M)
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

      allocate(atom%vbroad(n_cells))
      atom%vbroad = sqrt(vtherm/atom%weight * T + vturb**2) !vturb in m/s
      !atom%ntotal = atom%Abund * nHtot
      !note: for Hydrogen, ntotal is actually (nHtot - nHmin)

      VDoppler = sqrt(vtherm/atom%weight * maxval(T) + maxval(vturb)**2)

      !read bounb-bound (line) transitions
      allocate(atom%lines(atom%Nline))
      do kr=1,atom%Nline
         atom%lines(kr)%atom => atom
         atom%at(kr)%trtype = atom%lines(kr)%trtype; atom%at(kr)%ik=kr
         atom%at(kr)%lcontrib_to_opac=.true.
         atom%lines(kr)%symmetric = .false. !not used
         atom%lines(kr)%polarizable = .false.
         atom%lines(kr)%isotope_frac = 1.
         atom%lines(kr)%g_lande_eff = -99.0
         atom%lines(kr)%glande_i = -99.0; atom%lines(kr)%glande_j = -99.0
         !atom%lines(kr)%trtype="ATOMIC_LINE"
         call getnextline(atomunit, COMMENT_CHAR, FormatLine, inputline, Nread)                       
         Nread = len(trim(inputline))
         atom%lines(kr)%ZeemanPattern = 1 !should be read in file
         if (.not.lmagnetized) atom%lines(kr)%ZeemanPattern = 0
         ! -1 = effective triplet, +1
         if (Nread <= 115) then
            read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
               symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
               atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
               atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
               atom%lines(kr)%Grad, atom%lines(kr)%cStark
         else
            write(*,*) " ->Read aditional g_lande_eff for that line", kr
            read(inputline(1:Nread),*) j, i, f, shapeChar, atom%lines(kr)%Nlambda, &
               symmChar, atom%lines(kr)%qcore,atom%lines(kr)%qwing, vdWChar,&
               atom%lines(kr)%cvdWaals(1), atom%lines(kr)%cvdWaals(2), &
               atom%lines(kr)%cvdWaals(3), atom%lines(kr)%cvdWaals(4), &
               atom%lines(kr)%Grad, atom%lines(kr)%cStark, &
               atom%lines(kr)%g_Lande_eff, atom%lines(kr)%glande_j, atom%lines(kr)%glande_i
            !if glande_eff given, we need to read a value for gi and gj even if it is 0.
            !if glane is <= -99, gi, and gj and geff are computed eventually.
            !landé upper / lower levels in case the coupling scheme is not accurate
            if (atom%lines(kr)%g_lande_eff > -99) atom%lines(kr)%ZeemanPattern = -1 !effective T assumed
            if (atom%lines(kr)%g_lande_eff <= -99) atom%lines(kr)%ZeemanPattern = 1 !Full using gi and gj
            if (atom%lines(kr)%glande_j <= -99 .or. atom%lines(kr)%glande_i<= -99) then
               call Warning("Unable to use read lande factors, try to compute them..")
               atom%lines(kr)%g_lande_eff = -99 !force calculation
            end if
            if ((atom%lines(kr)%glande_j == 0.0 .and. atom%lines(kr)%glande_i == 0.0).or.&
               (atom%lines(kr)%g_lande_eff == 0.0)) then
               atom%lines(kr)%ZeemanPattern = 0
               write(*,*) " ++-> line", i, j, " unpolarized!"
            endif
         end if
         i = i + 1
         j = j + 1 !because RH uses C indexes

         atom%lines(kr)%i = min(i,j)
         atom%lines(kr)%j = max(i,j)

         !Temporary attributes which lines are stored for images
         !to be moved to parameter files
         if (atom%ID=="H") then

            !Ly alpha
            if (atom%g(i)==2 .and. atom%g(j)==8) atom%lines(kr)%write_flux_map =.true.
            !H alpha
            if (atom%g(i)==8 .and. atom%g(j)==18) atom%lines(kr)%write_flux_map=.true.
            !H beta
            if (atom%g(i)==8 .and. atom%g(j)==32) atom%lines(kr)%write_flux_map =.true.
            !H gamma
            if (atom%g(i)==8 .and. atom%g(j)==50) atom%lines(kr)%write_flux_map =.true.
            !Pa beta
            if (atom%g(i)==18 .and. atom%g(j)==50) atom%lines(kr)%write_flux_map =.true.
            !Br gamma
            if (atom%g(i)==32 .and. atom%g(j)==98) atom%lines(kr)%write_flux_map =.true.

         elseif (atom%ID=="He") then
            if ((kr == 5).or.(kr==15)) then
               atom%lines(kr)%write_flux_map =.true.
            endif
         else

            if (kr <= 4) then
               atom%lines(kr)%write_flux_map =.true.
            endif

         endif

         if (atom%lines(kr)%qwing < 2.0) then

            call Warning("qwing for line lower than 2! setting to 2")
            atom%lines(kr)%qwing = 2.0

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

         if (lmagnetized) then
            if ((abs(atom%qJ(atom%lines(kr)%i) - atom%qJ(atom%lines(kr)%j)) <= 1.) .and. &
               (atom%lines(kr)%g_Lande_eff <= -99) .and. &
               (determined(atom%lines(kr)%j)) .and. (determined(atom%lines(kr)%i))) then !
               ! do not compute geff if term is not
               ! determined or if the geff is read from file
               ! ie if g_lande_eff > -99
               ! fill lande_g_factor if determined and not given
               ! for b-b transitions
               call Lande_eff(atom, kr)
            end if
            atom%lines(kr)%polarizable = (atom%lines(kr)%g_lande_eff > -99).and.(atom%lines(kr)%ZeemanPattern /= 0)
            if (atom%lines(kr)%polarizable) then
               call ZeemanMultiplet(atom%lines(kr))
            endif
         endif !lmagnetized
         ! oscillator strength saved
         atom%lines(kr)%fosc = f
         lambdaji = (HPLANCK * CLIGHT) / (atom%E(j) - atom%E(i))
         atom%lines(kr)%Aji = C1 / (lambdaji**2) * (atom%g(i) / atom%g(j)) * f
         atom%lines(kr)%Bji = (lambdaji**3) / (2.0 * HPLANCK * CLIGHT) *atom%lines(kr)%Aji
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
            case ("PFR")
               write(*,*) "PFR not yet"
               stop
               atom%lines(kr)%Voigt = .true.
            case default
               write(*,*) "Line profile shape", shapechar, " unknown, using voigt."
               atom%lines(kr)%Voigt = .true.
         end select


         !Van der Waaks collision method
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

      end do !end loop over bound-bound transitions

      ! ----------------------------------------- !
      ! starts reading bound-free transitions
      ! cross-sections allocated once the final
      ! wavelength grid is known.
      ! ----------------------------------------- !
      allocate(atom%continua(atom%Ncont))
      do kr=1,atom%Ncont
         atom%continua(kr)%isotope_frac=1.
         atom%continua(kr)%atom => atom
         atom%at(kr)%lcontrib_to_opac=.true.
         atom%at(atom%Nline+kr)%trtype = atom%continua(kr)%trtype; atom%at(kr+atom%Nline)%ik=kr

         call getnextline(atomunit, COMMENT_CHAR, &
            FormatLine, inputline, Nread)
         read(inputline, *) j, i, atom%continua(kr)%alpha0,&
            atom%continua(kr)%Nlambda, nuDepChar, atom%continua(kr)%lambdamin
         j = j + 1
         i = i + 1
         atom%continua(kr)%j = max(i,j)
         atom%continua(kr)%i = min(i,j)
         lambdaji = (HPLANCK*CLIGHT)/ (atom%E(j)-atom%E(i))
         atom%continua(kr)%lambda0 = lambdaji/NM_TO_M !nm
         atom%continua(kr)%lambdamax = atom%continua(kr)%lambda0


         select case (nudepchar)
            case ("EXPLICIT")
               allocate(atom%continua(kr)%alpha_file(atom%continua(kr)%Nlambda))
               allocate(atom%continua(kr)%lambda_file(atom%continua(kr)%Nlambda))
               allocate(atom%continua(kr)%lambda(atom%continua(kr)%Nlambda))
               atom%continua(kr)%hydrogenic=.false.
               do la=atom%continua(kr)%Nlambda,1,-1
                  call getnextline(atomunit, COMMENT_CHAR, &
                       FormatLine, inputline, Nread)
                  read(inputline,*) atom%continua(kr)%lambda_file(la), &
                       atom%continua(kr)%alpha_file(la)
               end do
               atom%continua(kr)%lambda(:) = atom%continua(kr)%lambda_file(:)
               do la=2,atom%continua(kr)%Nlambda
                  if (atom%continua(kr)%lambda_file(la).lt.&
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
               !used to extrapolate the b-f cross-section for level dissolution
               !if I go to far in lambda0 actually g_bf is < 0 and atm I set g_bf = 0 if g_bf < 0.
               !to go further in lambda (i.e. for lambda with g_bf(lambda)->0), I need to properly extrapolated g_bf.
               !!atom%continua(kr)%lambdamax = atom%continua(kr)%lambda0 * fact_pseudo_cont!7
               !Meanwhile, I search the lambda for which g_bf is < 0 (then 0)
               if (ldissolve) then !only if we actually want to extrapolate
                  if (atom%ID=="H") then ! .or. atom%ID=="He") then
                     call search_cont_lambdamax (atom%continua(kr), atom%Rydberg, atom%stage(i)+1,atom%E(j),atom%E(i))
                     call make_sub_wavelength_grid_cont_linlog(atom%continua(kr), &
                   atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
                     !!CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
                  else
                     !there is dissolve but not for this atom
                     call make_sub_wavelength_grid_cont(atom%continua(kr), &
                          atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
                     ! 			call make_sub_wavelength_grid_cont_log_nu(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
                  endif
               else !no dissolve
                  ! call make_sub_wavelength_grid_cont(atom%continua(kr), &
                     !   atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
                  	call make_sub_wavelength_grid_cont_log_nu(atom%continua(kr), &
                        atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
               endif
               ! %lambda allocated inside the routines.
               !        CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
               !make it log if not occupation probability ?
               !        CALL make_sub_wavelength_grid_cont_log(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
               !Can be done elsewhere, with atomic lines ? But unlike some lines grid, this does not depend on the local condition ATM
            case default

               call error("nudepchar!")
            
         end select
      enddo !bound-free


      if (Nfixed.gt.0) call error("Fixed transitions no handled !")


      atom%set_ltepops = .true. !by default compute lte populations
      atom%NLTEpops = .false.

      ! allocate some space
      if (atom%initial_solution=="ZERO_RADIATION") then
         if (.not.atom%active) then
            write(*,*) atom%ID, " is passive! cannot use ZERO_RADIATION solution, set to LTE."
            atom%initial_solution="LTE_POPULATIONS"
         endif
      endif

      if (atom%active) then


         ! reading collision rates of RH
         if (atom%ID /= "H") then
            write(*,*) "  -> Reading collision data from RH for atom ", atom%ID
            ! 		atom%colunit = atom%periodic_table*2 + 1
            call read_collisions(atomunit, atom)
         endif

         !!allocate(atom%C(atom%Nlevel*atom%Nlevel,n_cells))
         !!now Collision matrix is constructed cell by cell, therefore allocated elsewhere
         allocate(atom%n(atom%Nlevel,n_cells))
         atom%n = 0d0
         if (atom%initial_solution == "LTE_POPULATIONS") then
            atom%n = atom%nstar !still need to be computed
            atom%set_ltepops = .true.
            write(*,*) " -> Setting initial solution to LTE "

         else if (atom%initial_solution == "CSWITCH") then
            atom%n = atom%nstar !still need to be computed
            atom%set_ltepops = .true.
            if (.not. lforce_lte) then !otherwise cswitch is not LTE
               write(*,*) " -> Setting initial solution to LTE with CSWITCH "
               atom%cswitch = cswitch_val
               if (.not. cswitch_enabled) cswitch_enabled = .true.!we need at least one
            endif
         else if (atom%initial_solution == "ZERO_RADIATION") then

            atom%n = atom%nstar
            atom%set_ltepops = .true.
            !nlte pops is false

         else if (atom%initial_solution == "OLD_POPULATIONS") then

            write(*,*) " -> Reading (non-LTE AND LTE) populations from file..."
            call read_pops_atom(atom)
            atom%NLTEpops = .true.
            atom%set_ltepops = .false. !read and USE also LTE populations from file!!

         end if
      else !not active = PASSIVE
         if (atom%initial_solution == "OLD_POPULATIONS") then

            allocate(atom%n(atom%Nlevel,n_cells)) !not allocated if passive, n->nstar
            write(*,*) " -> Reading (non-LTE AND LTE) populations from file for passive atom..."
            call read_pops_atom(atom)
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

   subroutine read_Atomic_Models(unit)
      !Read all atomic files present in atoms.input file
      !successive call of readModelAtom()
      integer :: EOF=0,Nread, nmet, mmet, nblancks, nact, npass
      integer :: kr, k, imax, ic
      real, parameter :: epsilon = 5e-3
      real :: eps
      real(kind=dp) :: epsilon_l_max !if epsilon > 1/pi/adamp, the value of xwing_lorentz is negative
      real(kind=dp) :: max_adamp, adamp, maxvel, vel, min_resol, max_resol
      integer, intent(in) :: unit
      character(len=MAX_LENGTH) :: inputline
      character(len=15) :: FormatLine
      character(len=MAX_LENGTH) :: popsfile, filename
      character(len=MAX_KEYWORD_SIZE) :: actionKey, popsKey
      character(len=2) :: IDread
      type (AtomType), pointer :: atom

      !create formatline
      !FormatLine = "(1A<MAX_LENGTH>)" !not working with ifort
      !write(FormatLine,'("("I3,A")")') MAX_LENGTH,"A"
      write(FormatLine,'("(1"A,I3")")') "A", MAX_LENGTH

    !   if (fact_pseudo_cont <= 0.0) then
    !    write(*,*) "fact_pseudo_cont=", fact_pseudo_cont
    !    call error("Wrong value for fact_pseudo_cont!")
    !   endif
    !   if (fact_pseudo_cont > 1.0) then
    !    Write(*,*) " Hydrogenic continua extrapolated up to lambda0 x ", fact_pseudo_cont
    !   endif

      Nactiveatoms = 0
      Npassiveatoms = 0
      open(unit=unit,file=TRIM(ATOMS_INPUT), status="old")!mcfost_utils)//TRIM(ATOMS_INPUT)

      !get number of atomic models to read
      call getnextline(unit, COMMENT_CHAR,FormatLine, &
         inputline, Nread)
      read(inputline,*) Natom
      if (Natom > 1) then
         write(*,*) "Reading ", Natom, " atom"
      else 
         write(*,*) "Reading ", Natom, " atoms"
      endif
      !array containing informations about all atoms
      allocate(Atoms(Natom))

      !go through file and read the different atomic models
      do nmet = 1, Natom
         call getnextline(unit, COMMENT_CHAR, FormatLine, inputline, Nread)

         allocate(Atoms(nmet)%ptr_atom)

         !the population file is not read anymore.
         ! read(inputline,'(1A28, 1A7, 1A22, 1A20)') filename, actionKey, popsKey, popsFile
         read(inputline,*) filename, actionKey, popsKey!, popsFile

         Atoms(nmet)%ptr_atom%initial_solution=adjustl(popsKey)
         Atoms(nmet)%ptr_atom%inputFile=trim(filename)


         ! would be pssoible in the future to read all J value also
         if (Atoms(nmet)%ptr_atom%initial_solution/="OLD_POPULATIONS" &
            .and.Atoms(nmet)%ptr_atom%initial_solution/="LTE_POPULATIONS"&
            .and.Atoms(nmet)%ptr_atom%initial_solution/="ZERO_RADIATION"&
            .and.Atoms(nmet)%ptr_atom%initial_solution/="CSWITCH")then
            write(*,*) "Initial solution ", Atoms(nmet)%ptr_atom%initial_solution,&
               " unkown!"
            write(*,*) "Exiting..."
            stop
         end if
         !!Atoms(nmet)%ptr_atom%dataFile = trim(popsFile) ! now the file atomID.fits.gz contains
         ! all informations including populations.


         !Active atoms are treated in NLTE
         if (adjustl(actionKey)=="ACTIVE") then
            Atoms(nmet)%ptr_atom%active=.true.
            !write(*,*) "atom is active"
            Nactiveatoms = Nactiveatoms+1
         else
            Atoms(nmet)%ptr_atom%active=.false.
            NpassiveAtoms = NpassiveAtoms+1
            !write(*,*) "atom is passive"
         end if
         !   ! just opoen the model to check that Hydrogen is the first model
         !   !
         open(unit=unit+nmet,file=trim(mcfost_utils)//TRIM(path_to_atoms)//trim(filename),status="old")
         call getnextline(unit+nmet, COMMENT_CHAR, FormatLine, inputline, Nread)
         read(inputline,*) IDread
         if (nmet==1 .and. IDread/="H ") then
            write(*,*) "First atomic model read has to be Hydrogen"
            write(*,*) "Exting..."
            stop
         end if
         close(unit+nmet)
         ! check duplicates
         IDread(2:2) = to_lower(IDread(2:2))
         do mmet = 1,nmet-1 !compare the actual (IDread) with previous
         !write(*,*) mmet, nmet
            if (Atoms(mmet)%ptr_atom%ID.eq.IDread) then
               write(*,*) "Already read a model for this atom ", IDread
               write(*,*) "exiting..."
               stop
            end if
         end do

         ! read and fill atomic structures saved in atoms
         ! create an alias for Hydrogen
         call read_Model_Atom(unit+nmet, Atoms(nmet)%ptr_atom, &
            trim(mcfost_utils)//TRIM(path_to_atoms)//trim(filename))

      end do
      close(unit)

      !Temporary
      !check that if helium is active and electron are iterated H must be active at the moment !!
      check_helium : do nmet=1, Natom
         if (atoms(nmet)%ptr_atom%id=="He") then

            if ((n_iterate_ne > 0).and.(atoms(nmet)%ptr_atom%active)) then

               if (atoms(nmet)%ptr_atom%stage(1) > 0) then
                  write(*,*) " !!!!!!!!! "
                  call warning("Helium ground state is not in neutral stage ! Must be He I")
                  write(*,*) atoms(nmet)%ptr_atom%stage(1), atoms(nmet)%ptr_atom%label(1), &
                     atoms(nmet)%ptr_atom%E(1), atoms(nmet)%ptr_atom%g(1)
                  write(*,*) " !!!!!!!!! "
                  stop
               endif

               !H always present and the first.
               !force to be active if helium active ?
               !at the moment print a warning!
               if (.not.atoms(1)%ptr_atom%active) then
                  write(*,*) " !!!!!!!!!!!!!!!!!!!!!!! "
                  call WARNING(" Hydrogen is passive, while helium is active and n_iterate_ne > 0!!")
                  write(*,*) " !!!!!!!!!!!!!!!!!!!!!!! "
                  exit check_helium
               endif

            endif

         endif
      enddo check_helium

      if (lfix_backgrnd_opac) then
         do nmet=1, Natom
            if (atoms(nmet)%ptr_atom%initial_solution.eq."OLD_POPULATIONS") then
               call warning(" Using previous NLTE populations with fixed bacgkrnd!!!!")
               write(*,*) " ----->n check consistency of the result"
               exit
            endif
         enddo
      endif


      ! Alias to the most importent one
      !always exits !!
      Hydrogen=>Atoms(1)%ptr_atom
      if (.not.associated(Hydrogen, Atoms(1)%ptr_atom)) CALL Error(" Hydrogen alias not associated to atomic model!")

      ! Aliases to active atoms
      nact = 0; npass = 0
      if (Nactiveatoms > 0) then
         allocate(ActiveAtoms(Nactiveatoms))
      end if
      if (Natom - Nactiveatoms == Npassiveatoms) then
         allocate(PassiveAtoms(Npassiveatoms))
      else
         write(*,*) "Error, number of passive atoms is not N Nactive"
         stop
      end if

      ! keep a duplicate in Elements
      write(*,*) "order#     ID   periodic-table#    ACTIVE    #lines   #continua"
      do nmet=1,Natom
         write(*,*) nmet, Atoms(nmet)%ptr_atom%ID, &
            Atoms(nmet)%ptr_atom%periodic_table, Atoms(nmet)%ptr_atom%active, &
            Atoms(nmet)%ptr_atom%Nline, Atoms(nmet)%ptr_atom%Ncont
         ! create alias in Elements for elements that have
         ! a model atom. It means all elements here.
         !atom => Atoms(nmet)
         Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%model => Atoms(nmet)%ptr_atom

         if (.not.associated(Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%model, &
            Atoms(nmet)%ptr_atom)) then
            write(*,*) Atoms(nmet)%ptr_atom%id
            call error(" Elemental model not associated to atomic model!")
         endif

         !Check Nstage
         if (Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%Nstage < maxval(Atoms(nmet)%ptr_atom%stage) + 1) then
            write(*,*) Atoms(nmet)%ptr_atom%id, maxval(Atoms(nmet)%ptr_atom%stage) + 1
            write(*,*) "Ns pf = ", Elements(Atoms(nmet)%ptr_atom%periodic_table)%ptr_elem%Nstage
            call error("Model has more ionisation stages than the one in the partition function!")
         endif

         if (Atoms(nmet)%ptr_atom%periodic_table==2)  then
            NULLIFY(Helium) !Because, it is associated to an Elem by default
            Helium => Atoms(nmet)%ptr_atom
            helium_is_active = helium%active
            write(*,*) "Helium pointers associated to atom in the table", nmet
            if (helium_is_active) then
               write(*,*) " And it is active !", Atoms(nmet)%ptr_atom%active
            endif
            if (.not.associated(Helium,Atoms(nmet)%ptr_atom)) &
               call Warning(" Helium alias is not associated to an atomic model!")
         end if

         if (allocated(ActiveAtoms)) then
            if (Atoms(nmet)%ptr_atom%active) then
               nact = nact + 1 !got the next index of active atoms
               ActiveAtoms(nact)%ptr_atom => Atoms(nmet)%ptr_atom
               Atoms(nmet)%ptr_atom%activeindex = nact
            end if
         end if

         if (allocated(PassiveAtoms)) then
            if (.not.Atoms(nmet)%ptr_atom%active) then
               npass = npass+1
               PassiveAtoms(npass)%ptr_atom => Atoms(nmet)%ptr_atom
            end if
         end if
      end do

      min_Resol = 1d30
      max_resol = 0.0
      do nmet=1,Natom
         do k=1,n_cells
            if (icompute_atomRT(k)>0) then
               min_resol = min(Atoms(nmet)%ptr_atom%vbroad(k), min_resol)
               max_resol = max(Atoms(nmet)%ptr_atom%vbroad(k), max_resol)
            endif
         enddo
      enddo
      hv = 0.46 * real(min_resol) * 1e-3

      if (art_hv > 0.0) then
         hv = art_hv
      endif
      write(*,'("R="(1F7.3)" km/s; min(Vth)="(1F7.3)" km/s; max(Vth)="(1F7.3)" km/s")') hv, min_resol * 1d-3, max_resol * 1d-3


      write(*,*) " Generating sub wavelength grid and lines boundary for all atoms..."
      do nmet=1, Natom
         atom => Atoms(nmet)%ptr_atom
         !!write(*,*) "ID:", Atoms(nmet)%ptr_atom%ID
         do kr=1, Atoms(nmet)%ptr_atom%Nline

            maxvel = Atoms(nmet)%ptr_atom%lines(kr)%qwing * maxval(atom%vbroad)

 
            call define_local_profile_grid (Atoms(nmet)%ptr_atom%lines(kr))


            !-> depends if we interpolate profile on finer grid !
            !!-> linear
            !      CALL make_sub_wavelength_grid_line_lin(Atoms(nmet)%ptr_atom%lines(kr),&
            !                                         maxval(Atoms(nmet)%ptr_atom%vbroad), max_adamp)
            !!-> logarithmic
            !      CALL make_sub_wavelength_grid_line(Atoms(nmet)%ptr_atom%lines(kr),&
            !                                         maxval(Atoms(nmet)%ptr_atom%vbroad), max_adamp)


         enddo !over lines
         if (atom%lgauss_prof) then
            write(*,*) " -> gauss profile for that atom ", atom%id
            call define_local_gauss_profile_grid(atom)
         endif
         atom => null()
      enddo !over atoms
      write(*,*) "..done"
      Max_ionisation_stage = get_max_nstage()

      return
   end subroutine read_Atomic_Models

   subroutine search_cont_lambdamax (cont, Rinf, Z, Ej, Ei)
      !Search the lambdamax = first lambda for which Gaunt < 0
      !Because I do not have an extrapolation routine, there is no need to extrapolate beyond
      !lambdamax if gaunt is negative.

      !See Gaunt_bf in Hydrogen.f90
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

         u = n_eff**2 * HPLANCK*CLIGHT / (NM_TO_M * cont%lambdamax) / Z*Z / E_RYDBERG - 1

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


   !But why, a cswitch per atom ? It is going at the same speed for all atoms right ?
   function maxval_cswitch_atoms ()
      !for all atoms, check the maximum value of the cswitch
      integer :: n
      real(kind=dp) ::  maxval_cswitch_atoms

      maxval_cswitch_atoms = 1.0_dp
      do n=1, Natom

         maxval_cswitch_atoms = max(maxval_cswitch_atoms, atoms(n)%ptr_atom%cswitch)

      enddo

      return
   end function maxval_cswitch_atoms

   !could be done for only one common parameter by the way
   subroutine adjust_cswitch_atoms ()
      !for all active atoms, decreases the cswitch value from cswitch_down_scaling_factor
      integer :: n
      logical :: print_message = .false.
      real(kind=dp) :: new_cs

      print_message = .false.

      do n=1, NactiveAtoms

         if (activeatoms(n)%ptr_atom%cswitch > 1.0) then
            activeatoms(n)%ptr_atom%cswitch = max(1.0_dp, min(activeatoms(n)%ptr_atom%cswitch, &
               activeatoms(n)%ptr_atom%cswitch/cswitch_down_scaling_factor))
            if (.not. print_message) then
               print_message = .true.
               new_cs = activeatoms(n)%ptr_atom%cswitch
            endif
         endif
      enddo

      if (print_message) write(*,'(" cswitch for next iteration: "(1ES17.8E3))') new_cs


      return
   end subroutine adjust_cswitch_atoms

end module readatom
!    !This is not used anymore., can be use to check atomZnumber in atom_type.f90
! function atomZnumber_old(atomID) result(Z)
!    !--------------------------------------
!    !return the atomic number of an atom
!    !with ID = atomID.
!    !Hydrogen is 1
!    !--------------------------------------
!     character(len=ATOM_ID_WIDTH) :: atomID
!     integer :: Z, i

!     Z = 1
!     do i=1,Nelem
!      if (Elements(i)%ptr_elem%ID == atomID) then
!        Z = i
!        exit
!      end if
!     enddo
!     write(*,*) "atomZnumber_old found:", Z
!     Z = 1
!     do while (Elements(Z)%ptr_elem%ID /= atomID)
!      Z=Z+1
!     end do
!     write(*,*) "atomZnumber_old found(bis):", Z


!    return
! end function atomZnumber_old
 !     if (trim(nuDepChar).eq."EXPLICIT") then
   !        ! Nlambda set in atomic file
   !        allocate(atom%continua(kr)%alpha_file(atom%continua(kr)%Nlambda))
   !        allocate(atom%continua(kr)%lambda_file(atom%continua(kr)%Nlambda))
   !        allocate(atom%continua(kr)%lambda(atom%continua(kr)%Nlambda))
   !        atom%continua(kr)%hydrogenic=.false.
   !        ! rearanging them in increasing order
   !        ! because in the atomic file they are
   !        ! given in decreasing order !
   !        do la=atom%continua(kr)%Nlambda,1,-1
   !           CALL getnextline(atomunit, COMMENT_CHAR, &
   !                FormatLine, inputline, Nread)
   !           read(inputline,*) atom%continua(kr)%lambda_file(la), &
   !                atom%continua(kr)%alpha_file(la)
   !           ! though they are printed in decreasing order
   !           !write(*,*) "l = ",atom%continua(kr)%lambda(la), &
   !           !    " a = ", atom%continua(kr)%alpha(la)
   !        end do
   !        atom%continua(kr)%lambda(:) = atom%continua(kr)%lambda_file(:)
   !        do la=2,atom%continua(kr)%Nlambda
   !           if (atom%continua(kr)%lambda_file(la).lt.&
   !                atom%continua(kr)%lambda_file(la-1)) then
   !              write(*,*) "continuum wavelength not monotonous"
   !              write(*,*) "exiting..."
   !              stop
   !           end if
   !        end do
   !        !not extrapolated if explicit ?
   !        !should consider neglecting occupation probability for that case
   !        atom%continua(kr)%lambdamax = maxval(atom%continua(kr)%lambda_file)
   !     else if (trim(nuDepChar).eq."HYDROGENIC") then
   !        atom%continua(kr)%hydrogenic=.true.
   !        !!tmp
   !        !        atom%continua(kr)%lambdamin = 5.0_dp
   !        !        atom%continua(kr)%lambdamin = 0.05 * atom%continua(kr)%lambda0
   !        !        atom%continua(kr)%lambdamin = max(10.0_dp, 1d-2 * atom%continua(kr)%lambda0)
   !        !!tmp
   !        write(*,'(" Continuum "(1I3)" -> "(1I3)" at "(1F12.5)" nm")') &
   !             atom%continua(kr)%i, atom%continua(kr)%j, atom%continua(kr)%lambda0
   !        write(*,'(" -> lower edge cut at "(1F12.5)" nm !")'), atom%continua(kr)%lambdamin

   !        if (atom%continua(kr)%lambdamin>=atom%continua(kr)%lambda0) then
   !           write(*,*) "Minimum wavelength for continuum is larger than continuum edge."
   !           write(*,*) kr, atom%continua(kr)%lambda0, atom%continua(kr)%lambdamin
   !           write(*,*) "exiting..."
   !           stop
   !        end if
   !        !used to extrapolate the b-f cross-section for level dissolution
   !        !if I go to far in lambda0 actually g_bf is < 0 and atm I set g_bf = 0 if g_bf < 0.
   !        !to go further in lambda (i.e. for lambda with g_bf(lambda)->0), I need to properly extrapolated g_bf.
   !        !!atom%continua(kr)%lambdamax = atom%continua(kr)%lambda0 * fact_pseudo_cont!7
   !        !Meanwhile, I search the lambda for which g_bf is < 0 (then 0)
   !        if (ldissolve) then !only if we actually want to extrapolate
   !           if (atom%ID=="H") then ! .or. atom%ID=="He") then
   !      	CALL search_cont_lambdamax (atom%continua(kr), atom%Rydberg, atom%stage(i)+1,atom%E(j),atom%E(i))
   !      	CALL make_sub_wavelength_grid_cont_linlog(atom%continua(kr), &
   !            atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !              !!CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !           else
   !              !there is dissolve but not for this atom
   !              CALL make_sub_wavelength_grid_cont(atom%continua(kr), &
   !                   atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !              ! 			call make_sub_wavelength_grid_cont_log_nu(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !           endif
   !        else !no dissolve
   !           CALL make_sub_wavelength_grid_cont(atom%continua(kr), &
   !                atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !           ! 			call make_sub_wavelength_grid_cont_log_nu(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !        endif
   !        ! %lambda allocated inside the routines.
   !        !        CALL make_sub_wavelength_grid_cont(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !        !make it log if not occupation probability ?
   !        !        CALL make_sub_wavelength_grid_cont_log(atom%continua(kr), atom%continua(kr)%lambdamin,atom%continua(kr)%lambdamax)
   !        !Can be done elsewhere, with atomic lines ? But unlike some lines grid, this does not depend on the local condition ATM
   !     else
   !        write(*,*) "Invalid continuum type : ", trim(nuDepChar)
   !        write(*,*) "exiting..."
   !        stop
   !     end if
   !  end do !end loop over bound-free transitions

    ! now fixed transitions
    ! fixed transitions are usefull to describe NLTE problem
    ! in the solar chromosphere in 2D,3D without not so much
    ! computational cost.
    ! They are not implemented because only relevent for solar
    ! application
         !      ic = find_continuum(atom,atom%lines(kr)%i)
          !
          !      max_adamp = 1d100
          !      maxvel = 0
          !      adamp = 0.0
          !      imax = 1
          !      !!write(*,*) " line @ ", Atoms(nmet)%ptr_atom%lines(kr)%lambda0,'nm'
          !      do k=1, n_cells
          !        if (icompute_atomRT(k)>0) then
          !
          !         vel = Atoms(nmet)%ptr_atom%vbroad(k)
          !
          !         if (Atoms(nmet)%ptr_atom%lines(kr)%voigt) then
          !          !no damping if Gaussian
          !          CALL Damping(k, Atoms(nmet)%ptr_atom, kr, adamp)
          !          !-> LTE populations not known, we cannot use damping here
          !          !Because some of the damping are zero if prop to populations
          !          !Typically there are only Stark and Natural depending on electrons and lines.
          !          !If only Adamp its to short
          !          !adamp = atom%lines(kr)%Aji * (NM_TO_M*atom%lines(kr)%lambda0) / (4.*pi) / atom%vbroad(k)
          !          !max_adamp = max(max_adamp, adamp)
          !          max_adamp = min(max_adamp, adamp) ! min actually
          !
          !          epsilon_l_max = 1.0/pi/adamp
          !
          !          if (n_eff(atom%Rydberg, atom%E(ic), atom%E(atom%lines(kr)%i), atom%stage(ic)) <= 3.0) then
          !           eps = 1e-4
          !          else
          !           eps = epsilon
          !          endif
          ! !-> with vD
          ! !         maxvel = max(maxvel, &
          ! !          vel * sqrt(abs(-log(epsilon))), vel * sqrt(min(adamp, 0.9*1.0/pi/epsilon)/pi/min(0.9*epsilon_l_max, epsilon) - min(adamp,0.9*1.0/pi/epsilon)**2) )
          ! !-> only from damping
          ! ! 		 maxvel = max(maxvel, vel * sqrt(min(adamp, 0.9*1.0/pi/eps)/pi/min(0.9*epsilon_l_max, eps) - min(adamp,0.9*1.0/pi/eps)**2))
          ! 		else
          ! 		 !-> exact
          ! ! 		 maxvel = max(maxvel, vel * sqrt(abs(-log(eps))))
          !         endif
          !
          !   		maxvel = max(maxvel, vel)
          !
          !
          !        endif
          !
          !      enddo
          !! maxvel = maxvel_atom_transfer * 1e3 * ( 1e-3 * vel/1.0)

          ! 	 maxvel = maxval(atom%vbroad)
          !write(*,*) "maxvel for line ", kr, maxvel * 1e-3
          !!-> group of lines, linear in v
          !for Gauss and voigt and sets line bounds!line%u deallocated if not needed