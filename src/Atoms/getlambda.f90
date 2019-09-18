MODULE getlambda

  use atom_type, only : AtomicContinuum, AtomicLine, AtomType
  use atmos_type, only : atmos, atomPointerArray
  use constant
  use getline, only : getnextline, MAX_LENGTH
  
  use parametres
  use utils, only : span, spanl, bubble_sort
  use messages
  
  IMPLICIT NONE
  
  !Number of points for each transition
  integer, parameter :: Nlambda_cont = 20, &
  						Nlambda_line_w = 9, Nlambda_line_c = 101
  						!Nwing from -vchar to vcore; Ncore from vcore to 0

  CONTAINS
  
  SUBROUTINE Read_wavelengths_table(lambda_table, Nlam_I)
  ! -------------------------------------------------------- !
  ! Read a wavelenth table of the form:
  ! Nregion
  ! Nlam_I(1)
  ! Nlam_I(2)
  ! ...
  ! Nlam_I(Nregion)
  ! lambda(1) (first region)
  ! ....
  ! lambda(sum(Nlam_I)) (last region, last point)
  ! -------------------------------------------------------- !

   real(kind=dp), dimension(:), allocatable, intent(inout) :: lambda_table
   integer, dimension(:), allocatable, intent(out) :: Nlam_I
   character(len=MAX_LENGTH) :: inputline, FormatLine
   integer :: Nread = 0, Nl, k, Nr, Nt
   
   if (.not.ltab_wavelength_image) RETURN
   
   if (allocated(lambda_table)) deallocate(lambda_table) !should not happen
   write(FormatLine,'("(1"A,I3")")') "A", 256

   open(unit=1, file=TRIM(tab_wavelength_image), status='old')
   CALL getnextline(1, "#", FormatLine, inputline, Nread)
   read(inputline,*) Nl
   allocate(Nlam_I(Nl)) !wavelength per regions
   
   do k=1,Nl
    CALL getnextline(1, "#", FormatLine, inputline, Nread)
    read(inputline,*) Nr
    Nlam_I(k) = Nr
   end do
   Nt = sum(Nlam_I)
 
   !read wavelengths
   allocate(lambda_table(Nt))
   do k=1,Nt
      CALL getnextline(1, "#", FormatLine, inputline, Nread)
      read(inputline,*) lambda_table(k)
   end do
   
   close(1)
   
  RETURN
  END SUBROUTINE Read_wavelengths_table
  
  SUBROUTINE make_sub_wavelength_grid_cont(cont, lambdamin)
  ! ----------------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicContinuum cont.
   ! The resolution is constant in nm.
   ! lambda must be lower that lambda0 and lambda(Nlambda)=lambda0.
   ! Allocate cont%lambda.
   ! cont%alpha (cross-section of photoionisation) is not used.
  ! ----------------------------------------------------------------- !
   type (AtomicContinuum), intent(inout) :: cont
   real(kind=dp), intent(in) :: lambdamin
   real(kind=dp) :: resol
   integer :: la
   real(kind=dp) :: l0, l1
   
   !write(*,*) "Atom for which the continuum belongs to:", cont%atom%ID

   l1 = cont%lambda0 !cannot be larger than lambda0 ! minimum frequency for photoionisation
   l0 = lambdamin
   cont%Nlambda = Nlambda_cont
   allocate(cont%lambda(cont%Nlambda))
   resol = (l1 - l0) / real(cont%Nlambda - 1, kind=dp)
!    write(*,*) "Continuum:", cont%lambda0, cont%j,"->",cont%i, &
!               " Resolution (nm):", resol, " lambdamin =", lambdamin  
   
   cont%lambda(1) = l0
   do la=2,cont%Nlambda
    cont%lambda(la) = cont%lambda(la-1) + resol
   end do
   !does not allocate cross-section, here
  RETURN
  END SUBROUTINE make_sub_wavelength_grid_cont
  
  !Probleme passe pas part lambda0
  SUBROUTINE make_sub_wavelength_grid_line_lin(line, vD)
  ! ------------------------------------------------------------ !
   ! Make an individual wavelength grid for the AtomicLine line.
   ! The wavelength grid is symmetric wrt lambda0.
   ! It is by default, logarithmic in the wing and linear in the
   ! core.
  ! ------------------------------------------------------------ !
   type (AtomicLine), intent(inout) :: line
   real(kind=dp), intent(in) :: vD !maximum thermal width of the atom in m/s
   real(kind=dp) :: v_char, dlam
   real(kind=dp) :: vB
   real(kind=dp) :: lam0, lam1
   integer :: la, Nlambda, Nmid
   real(kind=dp) :: adamp_char = 0d0
   real, parameter :: L = 0.55 ! L% of max extent contains the line
   real(kind=dp), dimension(Nlambda_line_c) :: xlam !dimension(2*(Nlambda_line_c+Nlambda_line_w-1)-1)

   adamp_char = 1d3 * line%Grad/ (4d0 * PI) * (NM_TO_M*line%lambda0) / vD
   vB = 0d0
   if (line%polarizable) vB =  &
   				2d0*atmos%B_char * LARMOR * (line%lambda0*NM_TO_M) * dabs(line%g_lande_eff)
   								 
   v_char = L * (atmos%v_char + 2d0*vD*(1. + adamp_char) + vB) !=maximum extension of a line

   xlam = 0d0
   lam0 = line%lambda0*(1-v_char/CLIGHT)
   lam1 = line%lambda0*(1+v_char/CLIGHT)

   !if (dabs(vel(Nlambda_line_w+Nlambda_line_c-1)) /= 0d0) vel(Nlambda_line_w+Nlambda_line_c-1) = 0d0
   !if (vel(Nlambda_line_w+Nlambda_line_c-1) /= 0) write(*,*) 'Vel(Nw+Nc-1) should be 0.0'
   
   Nlambda = Nlambda_line_c!2 * (Nlambda_line_w + Nlambda_line_c - 1) - 1 
   line%Nlambda = Nlambda
   if (mod(line%Nlambda,2)==0) line%Nlambda = line%Nlambda + 1
   !Nmid = Nlambda/2 + 1 !As Nlambda is odd '1--2--3', Nmid = N/2 + 1 = 2, because 3/2 = 1
   !						!because the division of too integers is the real part. 
   allocate(line%lambda(line%Nlambda))
   line%lambda(1) = lam0
   dlam = (lam1-lam0) / (real(line%Nlambda-1))
   do la=2,line%Nlambda
     line%lambda(la) = line%lambda(la-1) + dlam
   enddo

   !line%lambda(1:Nmid) = line%lambda0*(1d0 + vel(1:Nmid)/CLIGHT)
   !line%lambda(Nmid+1:Nlambda) = line%lambda0*(1d0 -vel(Nmid-1:1:-1)/CLIGHT)

   !if (line%lambda(Nmid) /= line%lambda0) write(*,*) 'Lambda(Nlambda/2+1) should be lambda0'

  RETURN
  END SUBROUTINE make_sub_wavelength_grid_line_lin
  
  SUBROUTINE make_sub_wavelength_grid_line(line, vD)
  ! ------------------------------------------------------------ !
   ! Make an individual wavelength grid for the AtomicLine line.
   ! The wavelength grid is symmetric wrt lambda0.
   ! It is by default, logarithmic in the wing and linear in the
   ! core.
  ! ------------------------------------------------------------ !
   type (AtomicLine), intent(inout) :: line
   real(kind=dp), intent(in) :: vD !maximum thermal width of the atom in m/s
   real(kind=dp) :: v_char, dvc, dvw
   real(kind=dp) :: vcore, vB, v0, v1!km/s
   integer :: la, Nlambda, Nmid
   real(kind=dp) :: adamp_char = 0d0
   real, parameter :: wing_to_core = 0.5, L = 2. !only if velocity field
   		!if it is the max, then L is close to 1, if it is the min, L >> 1, if it is the mean etc..
   !!integer, parameter :: Nc = 51, Nw = 7 !ntotal = 2*(Nc + Nw - 1) - 1
   real(kind=dp), dimension(2*(Nlambda_line_c+Nlambda_line_w-1)-1) :: vel
   !!real(kind=dp), dimension(2*(Nc+Nw-1)-1) :: vel !Size should be 2*(Nc+Nw-1)-1
   													 !if error try, 2*(Nc+Nw)

   adamp_char = 1d3 * line%Grad/ (4d0 * PI) * (NM_TO_M*line%lambda0) / vD
   !write(*,*) adamp_char*vD/1000.
   vB = 0d0
   if (line%polarizable) vB =  &
   				2d0*atmos%B_char * LARMOR * (line%lambda0*NM_TO_M) * dabs(line%g_lande_eff)
   								 
   v_char = L * (atmos%v_char + 2d0*vD*(1. + adamp_char) + vB) !=maximum extension of a line
   vcore = (atmos%v_char + 2d0*vD*(1. + adamp_char) + vB) * wing_to_core
   !transition between wing and core in velocity
   !!vcore = L * v_char * wing_to_core ! == fraction of line extent
   vcore = v_char * wing_to_core !with *L if velocity field

   v0 = -v_char !* L
   v1 = +v_char !* L
   vel = 0d0

   !from -v_char to 0
   dvw = (v_char-vcore)/real(Nlambda_line_w-1,kind=dp) !(L * v_char-vcore)/real(Nw-1,kind=dp), old
   dvc = vcore/real(Nlambda_line_c-1,kind=dp)
   
!! Linear wing
   !vel(1:Nw) = -real(span(real(v1), real(vcore), Nw),kind=dp)
   !vel(1) = v0 !half wing
   !!write(*,*) 1, vel(1), v0, L*v_char/1d3
   !wing loop
!    do la=2, Nw
!     vel(la) = vel(la-1) + dvw
!    !!write(*,*) la, vel(la)
!    end do
!! Log wing
   !should be log from vcore to v0 not the opposite
   !v0 is < 0 but spanl takes the abs
!   vel(1:Nw) = -real(spanl(real(v0), real(vcore), Nw),kind=dp)
   vel(Nlambda_line_w:1:-1) = -real(spanl(real(vcore), real(v0), Nlambda_line_w),kind=dp)
!    do la=1,Nw
!     write(*,*) la, vel(la)/1000., vcore, -v_char, adamp_char*vD / 1000.
!    end do
!! end scale of wing points
!   vel(Nw:Nw-1+Nc) = -real(span(real(vcore), real(0.), Nc+1),kind=dp) 
!   write(*,*) Nw, vel(Nw), vcore
   !vel(Nw) = -vcore!should be okey at the numerical precision 
   !la goes from 1 to Nw + Nc -1 total number of points.
   !if Nc=101 and Nw =11, 111 velocity points,because the point vel(11) of the wing grid
   !is shared with the point vel(1) of the core grid.
   do la=Nlambda_line_w+1, Nlambda_line_c+Nlambda_line_w-1 !Half line core
    vel(la) = vel(la-1) + dvc
    !write(*,*) la-Nw+1,la, vel(la) !relative index of core grid starts at 2 because vel(Nw)
    							   ! is the first point.
   end do   
   
   !! Just a check here, maybe forcing the mid point to be zero is brutal
   !! but by construction it should be zero !
   !if (dabs(vel(Nw+Nc-1)) <= 1d-7) vel(Nw+Nc-1) = 0d0
   if (dabs(vel(Nlambda_line_w+Nlambda_line_c-1)) /= 0d0) vel(Nlambda_line_w+Nlambda_line_c-1) = 0d0
   if (vel(Nlambda_line_w+Nlambda_line_c-1) /= 0) write(*,*) 'Vel(Nw+Nc-1) should be 0.0'

  !number of points from -vchar to 0 is Nw+Nc-1, -1 because otherwise we count
  ! 2 times vcore which is shared by the wing (upper boundary) and core (lower boundary) grid.
  ! Total number of points is 2*(Nw+Nc-1) but here we count 2 times lambda0., therefore
  ! we remove 1 point.
   Nlambda = 2 * (Nlambda_line_w + Nlambda_line_c - 1) - 1 
   line%Nlambda = Nlambda
   Nmid = Nlambda/2 + 1 !As Nlambda is odd '1--2--3', Nmid = N/2 + 1 = 2, because 3/2 = 1
   						!because the division of too integers is the real part. 
   allocate(line%lambda(line%Nlambda))

   line%lambda(1:Nmid) = line%lambda0*(1d0 + vel(1:Nmid)/CLIGHT)
   line%lambda(Nmid+1:Nlambda) = line%lambda0*(1d0 -vel(Nmid-1:1:-1)/CLIGHT)

   if (line%lambda(Nmid) /= line%lambda0) write(*,*) 'Lambda(Nlambda/2+1) should be lambda0'
   
!    if (line%j==3 .and. line%i==2) then
!    do la=1,Nlambda
!     write(*,*) la, line%lambda(la), line%lambda0
!    end do
!    !stop
!    end if

  RETURN
  END SUBROUTINE make_sub_wavelength_grid_line  

  SUBROUTINE make_wavelength_grid(Natom, Atoms, inoutgrid, Ntrans, wl_ref)
  use math, only : locate
  ! --------------------------------------------------------------------------- !
   ! construct and sort a wavelength grid for atomic line radiative transfer.
   ! Computes also the edge of a line profile: Nblue and Nred.
   ! The grid is built by merging the individual wavelength grids of each 
   ! transitions. Duplicates are removed and therefore changes the value of
   ! line%Nlambda and continuum%Nlambda.
   ! line%lambda and continuum%lambda are useless now. Except that 
   ! continuu%lambda is still used in .not.continuum%Hydrogenic !!
  ! --------------------------------------------------------------------------- !
   type (atomPointerArray), intent(inout), dimension(Natom) :: Atoms
   integer, intent(in) :: Natom
   real(kind=dp), intent(in) :: wl_ref
   integer, intent(out) :: Ntrans !Total number of transitions (cont + line)
   ! output grid. May contain values that are added to the final list before
   ! deallocating the array. It is reallocated when the final list is known.
   real(kind=dp), allocatable, dimension(:), intent(inout) :: inoutgrid
   ! temporary storage for transitions, to count and add them.
   integer, parameter :: MAX_TRANSITIONS = 50000
   type (AtomicLine), allocatable, dimension(:) :: alllines
   type (AtomicContinuum), allocatable, dimension(:) :: allcont
   integer :: kr, kc, n, Nspect, Nwaves, Nlinetot, Nctot
   integer :: la, nn, nnn !counters: wavelength, number of wavelengths, line wavelength
   integer :: Nred, Nblue, Nlambda_original!read from model
   real(kind=dp), allocatable, dimension(:) :: tempgrid, sorted_indexes
   real(kind=dp) :: l0, l1 !ref wavelength of each transitions
   
   write(*,*) ' Defining the nLTE wavelength grid, using ', Nlambda_cont,' points for each continuum, and ', &
    2*(Nlambda_line_w+Nlambda_line_c-1)-1, " points for each line."

   nn = 0
   allocate(alllines(MAX_TRANSITIONS), allcont(MAX_TRANSITIONS))!stored on heap
   ! if allocated inoutgrid then add to the number of
   ! wavelength points to the points of the inoutgrid.
   Nspect = 0
   if (allocated(inoutgrid)) Nspect = size(inoutgrid)

   Nlinetot = 0
   Nctot = 0
   do n=1,Natom
   !unlike RH, even passive atoms have dedicated wavelength grids
   ! Count number of total wavelengths
   do kr=1,Atoms(n)%ptr_atom%Nline
    Nspect = Nspect + Atoms(n)%ptr_atom%lines(kr)%Nlambda
    Nlinetot = Nlinetot + 1
    if (Nlinetot > MAX_TRANSITIONS) then
     write(*,*) "too many transitions"
     stop
    end if
    alllines(Nlinetot) = Atoms(n)%ptr_atom%lines(kr)
   end do
   do kc=1,Atoms(n)%ptr_atom%Ncont
    Nspect = Nspect + Atoms(n)%ptr_atom%continua(kc)%Nlambda
    Nctot = Nctot + 1
    if (Nctot > MAX_TRANSITIONS) then
     write(*,*) "too many transitions"
     stop
    end if
    allcont(Nctot) = Atoms(n)%ptr_atom%continua(kc)
   end do
  end do ! end loop over atoms


  ! add ref wavelength if any and allocate temp array
  if (wl_ref > 0.) then
   Nspect = Nspect + 1
   nn = 1
   allocate(tempgrid(Nspect))
   tempgrid(1)=wl_ref
  else
   allocate(tempgrid(Nspect))
  end if

  Ntrans = Nlinetot + Nctot
  write(*,*) "Adding ", Nspect," wavelengths for a total of", Ntrans, &
             " transitions."
  if (allocated(inoutgrid)) write(*,*) "  ->", size(inoutgrid)," input wavelengths"

!   allocate(Nred_array(Ntrans), Nblue_array(Ntrans), Nmid_array(Ntrans))

  ! add wavelength from mcfost inoutgrid if any
  ! and convert it to nm, then deallocate
  if (allocated(inoutgrid)) then
   do la=1, size(inoutgrid)
    nn = nn + 1
    tempgrid(nn) = inoutgrid(la) !nm or convert to nm here
   end do
   deallocate(inoutgrid)
  end if

  ! start adding continua and lines wavelength grid
  nnn = 0!it is just an indicative counter
  		 ! because nn is dependent on wl_ref and on the inoutgrid
  		 ! it is filled. 
  do kc=1,Nctot
   do la=1,allcont(kc)%Nlambda
    nn = nn + 1
    nnn = nnn + 1
    tempgrid(nn) = allcont(kc)%lambda(la)
   end do
  end do
  write(*,*) "  ->", nnn," continuum wavelengths"
  nnn = 0
  do kr=1,Nlinetot
   do la=1,alllines(kr)%Nlambda
    nn = nn + 1
    nnn = nnn + 1
    tempgrid(nn) = alllines(kr)%lambda(la)
   end do
  end do
  write(*,*) "  ->", nnn," line wavelengths"
  if (wl_ref > 0)   write(*,*) "  ->", 1," reference wavelength at", wl_ref, 'nm'
  ! sort wavelength
  !!CALL sort(tempgrid, Nspect)
  !this should work ?
  allocate(sorted_indexes(Nspect))
  sorted_indexes = bubble_sort(tempgrid)
  tempgrid(:) = tempgrid(sorted_indexes)

  !check for dupplicates
  !tempgrid(1) already set
  Nwaves = 2
  !write(*,*) "l0", tempgrid(1)
  do la=2,Nspect
  !write(*,*) "dlam>0", tempgrid(la)-tempgrid(la-1)
   if (tempgrid(la) > tempgrid(la-1)) then
    tempgrid(Nwaves) = tempgrid(la)
    !write(*,*) tempgrid(la), tempgrid(la-1), Nwaves+1, la, Nspect, tempgrid(Nwaves)
    Nwaves = Nwaves + 1
   end if
  end do
  Nwaves = Nwaves-1 ! I don't understand yet but it works
  					! Maybe, if only 1 wavelength, Nwaves = 2 - 1 = 1 (and not 2 as it 
  					! is implied by the above algorithm...)
  					! or if 2 wavelengths, Nwaves = 2 = 3 - 1 (not 3 again etc)
  					! and so on, actually I have 1 wavelength more than I should have.
  write(*,*) Nwaves, " unique wavelengths: ", Nspect-Nwaves,&
   " eliminated lines"

  ! allocate and store the final grid
  allocate(inoutgrid(Nwaves))
  do la=1,Nwaves
   inoutgrid(la) = tempgrid(la)
   ! write(*,*) "lam(la) =", inoutgrid(la)
  end do

  !!should not dot that but error somewhere if many atoms
  !!CALL sort(inoutgrid, Nwaves)

  !free some space
  deallocate(tempgrid, alllines, allcont, sorted_indexes)

!   Now replace the line%Nlambda and continuum%Nlambda by the new values.
!   we do that even for PASSIVE atoms
!   deallocate line%lambda because it does not correspond to the new
!   Nblue and Nlambda anymore.
! However, cont%lambda is need for interpolation of cont%alpha on the new grid if
! cont is not hydrogenic
!   nn = 1
  do n=1,Natom
   !first continuum transitions
!   write(*,*) " ------------------------------------------------------------------ "
   do kc=1,Atoms(n)%ptr_atom%Ncont
    Nlambda_original = Atoms(n)%ptr_atom%continua(kc)%Nlambda
    l0 = Atoms(n)%ptr_atom%continua(kc)%lambda(1)
    l1 = Atoms(n)%ptr_atom%continua(kc)%lambda(Nlambda_original) !on the subgrid
!    Nred = locate(inoutgrid,l1)
!     Nblue = locate(inoutgrid,l0)
!      write(*,*) locate(inoutgrid,l0), locate(inoutgrid,l1), l0, l1!, Nblue, Nred
    Atoms(n)%ptr_atom%continua(kc)%Nblue = locate(inoutgrid,l0)
    Atoms(n)%ptr_atom%continua(kc)%Nred = locate(inoutgrid,l1)
    Atoms(n)%ptr_atom%continua(kc)%Nlambda = Atoms(n)%ptr_atom%continua(kc)%Nred - &
                                    Atoms(n)%ptr_atom%continua(kc)%Nblue + 1
!      write(*,*) Atoms(n)%ID, " continuum:",kr, " Nlam_ori:", Nlambda_original, &
!      " l0:", l0, " l1:", l1, " Nred:",  Atoms(n)%continua(kc)%Nred, & 
!        " Nblue:", Atoms(n)%continua(kc)%Nblue, " Nlambda:", Atoms(n)%continua(kc)%Nlambda, & 
!        " Nblue+Nlambda-1:", Atoms(n)%continua(kc)%Nblue + Atoms(n)%continua(kc)%Nlambda - 1
!! For continuum transitions, lambda0 is at Nred, check definition of the wavelength grid
!! which means that cont%Nmid = locate(inoutgrid, lam(Nred)+lam(Nb)/(Nlambda))
!! and l1, lam(Nlambda) = lambda0
    Atoms(n)%ptr_atom%continua(kc)%Nmid = locate(inoutgrid,0.5*(l0+l1))
    if (Atoms(n)%ptr_atom%continua(kc)%Hydrogenic) deallocate(atoms(n)%ptr_atom%continua(kc)%lambda) !kept only if we interpolate
!!deprecated
!     CALL fillPhotoionisationCrossSection(Atoms(n)%ptr_atom, kc, &
!     	Atoms(n)%ptr_atom%continua(kc)%lambda,Nwaves, inoutgrid)
!     	
!     !TEST 
!     Atoms(n)%ptr_atom%continua(kc)%Nlambda = size( Atoms(n)%ptr_atom%continua(kc)%alpha)
      !!-> also deallocates cont%lambda which is the original lambda read from file
    !allocate(Atoms(n)%continua(kc)%lambda(Atoms(n)%continua(kc)%Nlambda))
    !Atoms(n)%continua(kc)%lambda(Atoms(n)%continua(kc)%Nblue:Atoms(n)%continua(kc)%Nred) &
    ! = inoutgrid(Atoms(n)%continua(kc)%Nblue:Atoms(n)%continua(kc)%Nred)
!!!     Nred_array(nn) = Atoms(n)%continua(kc)%Nred
!!!     Nmid_array(nn) = Atoms(n)%continua(kc)%Nmid
!!!     Nblue_array(nn) = Atoms(n)%continua(kc)%Nblue
!!!     nn= nn + 1
   end do
   !then bound-bound transitions
   do kr=1,Atoms(n)%ptr_atom%Nline
    Nlambda_original = Atoms(n)%ptr_atom%lines(kr)%Nlambda
    l0 = Atoms(n)%ptr_atom%lines(kr)%lambda(1)
    l1 = Atoms(n)%ptr_atom%lines(kr)%lambda(Nlambda_original)
!    Nred = locate(inoutgrid,l1)
!    Nblue = locate(inoutgrid,l0)
!     write(*,*) locate(inoutgrid,l0), locate(inoutgrid,l1), l0, l1!, Nblue, Nred
    Atoms(n)%ptr_atom%lines(kr)%Nblue = locate(inoutgrid,l0)!Nblue
    Atoms(n)%ptr_atom%lines(kr)%Nred = locate(inoutgrid,l1)
    Atoms(n)%ptr_atom%lines(kr)%Nlambda = Atoms(n)%ptr_atom%lines(kr)%Nred - &
                                 Atoms(n)%ptr_atom%lines(kr)%Nblue + 1
!     write(*,*) Atoms(n)%ID, " line:",kr, " Nlam_ori:", Nlambda_original, &
!     " l0:", l0, " l1:", l1, " Nred:",  Atoms(n)%lines(kr)%Nred, & 
!       " Nblue:", Atoms(n)%lines(kr)%Nblue, " Nlambda:", Atoms(n)%lines(kr)%Nlambda, & 
!       " Nblue+Nlambda-1:", Atoms(n)%lines(kr)%Nblue + Atoms(n)%lines(kr)%Nlambda - 1
    Atoms(n)%ptr_atom%lines(kr)%Nmid = locate(inoutgrid,Atoms(n)%ptr_atom%lines(kr)%lambda0)
    deallocate(Atoms(n)%ptr_atom%lines(kr)%lambda) !does not correpond to the new grid, indexes might be wrong
    !allocate(Atoms(n)%ptr_atom%lines(kr)%lambda(Atoms(n)%ptr_atom%lines(kr)%Nlambda))
    !Atoms(n)%lines(kr)%lambda(Atoms(n)%lines(kr)%Nblue:Atoms(n)%lines(kr)%Nred) &
    ! = inoutgrid(Atoms(n)%lines(kr)%Nblue:Atoms(n)%lines(kr)%Nred)
!!!     Nred_array(nn) = Atoms(n)%lines(kr)%Nred
!!!     Nmid_array(nn) = Atoms(n)%lines(kr)%Nmid
!!!     Nblue_array(nn) = Atoms(n)%lines(kr)%Nblue
!!!     nn = nn + 1
   end do
!     write(*,*) " ------------------------------------------------------------------ "
  end do !over atoms
  
  RETURN
  END SUBROUTINE make_wavelength_grid

  SUBROUTINE adjust_wavelength_grid(old_grid, lambda, Lam_region, Atoms)
   ! ------------------------------------------ !
    ! Reallocate wavelengths and indexes arrays
    ! to compute images on a user defined grid
   ! ------------------------------------------ !
   use math, only : locate
   use atmos_type, only : realloc_Transitions
   real(kind=dp), dimension(:), intent(in) :: old_grid
   integer, dimension(:), intent(in) :: lam_region
   type (atomPointerArray), dimension(:), intent(inout) :: Atoms
   real(kind=dp), dimension(:), intent(inout) :: lambda
   real(kind=dp), dimension(size(lambda)) :: lambda_us
   integer :: Nwaves, n, kr, kc, Nlambda_original, Nblue, Nred, Natom, ll, lll
   real(kind=dp) :: l0, l1 !ref wavelength of each transitions
   real(kind=dp) :: x0, x1 !bound of the new grid
   integer, dimension(:), allocatable :: sorted_indexes
   logical, dimension(:), allocatable :: trans_contribute
   type (AtomicContinuum), dimension(:), allocatable :: conta
   type (AtomicLine), dimension(:), allocatable :: lines
   logical :: in_chan

   Natom = size(Atoms)
   Nwaves = size(lambda)
   !check lambda is sorted ?
   !--> moved after the test over transitions now
   allocate(sorted_indexes(Nwaves))
   lambda_us(:) = lambda(:)
   sorted_indexes = bubble_sort(lambda)
   lambda(:) = lambda(sorted_indexes)
   x0 = minval(lambda); x1 = maxval(lambda)

   !Realloc space for atoms
   !we need to test if a transition is on the new grid or not. Because the final grid
   !is not the sum of the individual grid, some transitions can be neglected because
   !they are out of range
   do n=1,Natom
       
    allocate(trans_contribute(atoms(n)%ptr_atom%Ntr)); trans_contribute(:)=.true.!by default
    
    do kc=1,Atoms(n)%ptr_atom%Ncont
     !on the old_Grid (i.e., the grid for NLTE which is built using all transitions)
     Nlambda_original = Atoms(n)%ptr_atom%continua(kc)%Nlambda !on the old_grid
     Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue
     Nred = Atoms(n)%ptr_atom%continua(kc)%Nred
     l0 = Atoms(n)%ptr_atom%continua(kc)%lambdamin!old_grid(Nblue)
     l1 = Atoms(n)%ptr_atom%continua(kc)%lambda0!old_grid(Nred)

     in_chan = .false. !equivalent of trans_contribute so be smarter please
     ll = 0
     region_loop : do kr=1, size(lam_region)
     !relative index of regions
     ll = 1 + ll; lll = sum(lam_region(1:kr))
     x0 = minval(lambda_us(ll:lll)); x1 = maxval(lambda_us(ll:lll))
     if (l1 <= x0.or. l0 >= x1) then
      in_chan = .false.
     else
      in_chan = .true.
      exit region_loop !because if in one region no need to test the others
     end if
     ll = ll + sum(lam_region(1:kr))
     end do region_loop
     if (in_chan) then
      Atoms(n)%ptr_atom%continua(kc)%Nred = locate(lambda,l1) ! closest value return by locate is Nwaves if l1>lambda(Nwaves)
      Atoms(n)%ptr_atom%continua(kc)%Nblue = locate(lambda,l0) !
      Nred = Atoms(n)%ptr_atom%continua(kc)%Nred; Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue
     else
      Atoms(n)%ptr_atom%continua(kc)%Nred = -99
      Atoms(n)%ptr_atom%continua(kc)%Nblue = -99
     end if

     
     Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue; Nred = Atoms(n)%ptr_atom%continua(kc)%Nred
     if (Nred==-99.and.Nblue==-99) then
      Atoms(n)%ptr_atom%continua(kc)%Nlambda = -99

      Atoms(n)%ptr_atom%continua(kc)%Nmid = -99
      
      !!Atoms(n)%ptr_atom%continua(kc)%lcontrib_to_opac=.false.
      trans_contribute(atoms(n)%ptr_atom%Nline+kc) = .false.
      Atoms(n)%ptr_atom%at(atoms(n)%ptr_atom%Nline+kc)%lcontrib_to_opac=.false.
      write(*,*) " :: b-f transition", Atoms(n)%ptr_atom%continua(kc)%j,"->",Atoms(n)%ptr_atom%continua(kc)%i,&
       " for atom ",Atoms(n)%ptr_atom%ID, l0,"-",l1," not counted." !, " removed."
     else
      Atoms(n)%ptr_atom%continua(kc)%Nlambda = Atoms(n)%ptr_atom%continua(kc)%Nred - &
                                 Atoms(n)%ptr_atom%continua(kc)%Nblue + 1

      Atoms(n)%ptr_atom%continua(kc)%Nmid = locate(lambda,lambda(Atoms(n)%ptr_atom%continua(kc)%Nlambda)/2+1)

     end if                          
    end do

    !then bound-bound transitions
    do kr=1,Atoms(n)%ptr_atom%Nline
     Nlambda_original = Atoms(n)%ptr_atom%lines(kr)%Nlambda !on the old_grid
     Nblue = Atoms(n)%ptr_atom%lines(kr)%Nblue
     Nred = Atoms(n)%ptr_atom%lines(kr)%Nred
     l0 = old_grid(Nblue)
     l1 = old_grid(Nred)

     in_chan = .false.
     ll = 0
     !!write(*,*) size(lam_region), l0, l1, Atoms(n)%ptr_atom%lines(kr)%lambda0
     region_loop_l : do kc=1, size(lam_region)
     !relative index of regions
     ll = 1 + ll; lll = sum(lam_region(1:kc))
     !!write(*,*) kc, ll, lll, lam_region(kc), lll-ll+1
     x0 = minval(lambda_us(ll:lll)); x1 = maxval(lambda_us(ll:lll))
     !!write(*,*) x0, x1
     if (l1 <= x0.or. l0 >= x1) then
      in_chan = .false.
     else
      in_chan = .true.
      exit region_loop_l !because if in one region no need to test the others
     end if
     ll = sum(lam_region(1:kc))
     end do region_loop_l
     !!write(*,*) in_chan

     if (in_chan) then
      Atoms(n)%ptr_atom%lines(kr)%Nblue = locate(lambda,l0) ! closest value return by locate is Nwaves if l1>lambda(Nwaves)
      Atoms(n)%ptr_atom%lines(kr)%Nred = locate(lambda,l1) !
     else
      Atoms(n)%ptr_atom%lines(kr)%Nblue = -99
      Atoms(n)%ptr_atom%lines(kr)%Nred = -99
     end if
     
     Nblue = Atoms(n)%ptr_atom%lines(kr)%Nblue; Nred = Atoms(n)%ptr_atom%lines(kr)%Nred
     if (Nred==-99.and.Nblue==-99) then
      Atoms(n)%ptr_atom%lines(kr)%Nlambda = -99

      Atoms(n)%ptr_atom%lines(kr)%Nmid = -99

      !!Atoms(n)%ptr_atom%lines(kr)%lcontrib_to_opac=.false.
      trans_contribute(kr) = .false.
      Atoms(n)%ptr_atom%at(kr)%lcontrib_to_opac=.false.
      write(*,*) " :: b-b transition", Atoms(n)%ptr_atom%lines(kr)%j,"->",Atoms(n)%ptr_atom%lines(kr)%i,&
       " for atom " ,Atoms(n)%ptr_atom%ID, l0,"-",l1," not counted."!, " removed."
     else
      Atoms(n)%ptr_atom%lines(kr)%Nlambda = Atoms(n)%ptr_atom%lines(kr)%Nred - &
                                 Atoms(n)%ptr_atom%lines(kr)%Nblue + 1

      Atoms(n)%ptr_atom%lines(kr)%Nmid = locate(lambda,Atoms(n)%ptr_atom%lines(kr)%lambda0)
     end if                               
    if (allocated(Atoms(n)%ptr_atom%lines(kr)%lambda)) &
    	deallocate(Atoms(n)%ptr_atom%lines(kr)%lambda)
    end do
    
     CALL realloc_transitions(Atoms(n)%ptr_atom, count(trans_contribute), trans_contribute)
     deallocate(trans_contribute)
   end do !over atoms

  RETURN
  END SUBROUTINE adjust_wavelength_grid

  END MODULE getlambda