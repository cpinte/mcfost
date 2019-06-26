MODULE getlambda

  use atom_type, only : AtomicContinuum, AtomicLine, AtomType
  use atmos_type, only : atmos, atomPointerArray
  use constant
  use getline, only : getnextline, MAX_LENGTH
  
  use parametres
  use utils, only : span, spanl, bubble_sort
  use messages
  
  IMPLICIT NONE
  
!   double precision, dimension(:), allocatable :: Nred_array, Nblue_array, Nmid_array

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

   double precision, dimension(:), allocatable, intent(inout) :: lambda_table
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
   double precision, intent(in) :: lambdamin
   double precision :: resol
   integer, parameter :: Nlambda = 71!101
   integer :: la
   double precision :: l0, l1
   
   !write(*,*) "Atom for which the continuum belongs to:", cont%atom%ID

   l1 = cont%lambda0 !cannot be larger than lambda0 ! minimum frequency for photoionisation
   l0 = lambdamin
   cont%Nlambda = Nlambda
   allocate(cont%lambda(cont%Nlambda))
   resol = (l1 - l0) / real(Nlambda - 1, kind=dp)
!    write(*,*) "Continuum:", cont%lambda0, cont%j,"->",cont%i, &
!               " Resolution (nm):", resol, " lambdamin =", lambdamin  
   
   cont%lambda(1) = l0
   do la=2,cont%Nlambda
    cont%lambda(la) = cont%lambda(la-1) + resol
   end do
   !does not allocate cross-section, here
  RETURN
  END SUBROUTINE make_sub_wavelength_grid_cont
  
  SUBROUTINE make_sub_wavelength_grid_line(line, vD)
  ! ------------------------------------------------------------ !
   ! Make an individual wavelength grid for the AtomicLine line.
   ! The wavelength grid is symmetric wrt lambda0.
   ! It is by default, logarithmic in the wing and linear in the
   ! core.
  ! ------------------------------------------------------------ !
   type (AtomicLine), intent(inout) :: line
   double precision, intent(in) :: vD !maximum thermal width of the atom in m/s
   double precision :: v_char, dvc, dvw
   double precision :: vcore, vB, v0, v1!km/s
   integer :: la, Nlambda, Nmid
   double precision :: adamp_char = 0d0
   double precision, parameter :: wing_to_core = 0.5, L = 5d0 !only if velocity field
   		!if it is the max, then L is close to 1, if it is the min, L >> 1, if it is the mean etc..
   integer, parameter :: Nc = 45, Nw = 7 !ntotal = 2*(Nc + Nw - 1) - 1
   double precision, dimension(2*(Nc+Nw-1)-1) :: vel !Size should be 2*(Nc+Nw-1)-1
   													 !if error try, 2*(Nc+Nw)

   adamp_char = 1d3 * L * line%Grad/ (4d0 * PI) * (NM_TO_M*line%lambda0) / vD
   !write(*,*) adamp_char*vD/1000.
   vB = 0d0
   if (line%polarizable) vB =  &
   				2d0*atmos%B_char * LARMOR * (line%lambda0*NM_TO_M) * dabs(line%g_lande_eff)
   								 
   v_char = (atmos%v_char*L + 2d0*vD*(1. + adamp_char) + vB) !=maximum extension of a line
   vcore = (atmos%v_char + 2d0*vD*(1. + adamp_char) + vB) * wing_to_core
   !transition between wing and core in velocity
   !!vcore = L * v_char * wing_to_core ! == fraction of line extent
   vcore = v_char * wing_to_core !with *L if velocity field

   v0 = -v_char !* L
   v1 = +v_char !* L
   vel = 0d0

   !from -v_char to 0
   dvw = (v_char-vcore)/real(Nw-1,kind=dp) !(L * v_char-vcore)/real(Nw-1,kind=dp), old
   dvc = vcore/real(Nc-1,kind=dp)
   
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
   vel(Nw:1:-1) = -real(spanl(real(vcore), real(v0), Nw),kind=dp)
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
   do la=Nw+1, Nc+Nw-1 !Half line core
    vel(la) = vel(la-1) + dvc
    !write(*,*) la-Nw+1,la, vel(la) !relative index of core grid starts at 2 because vel(Nw)
    							   ! is the first point.
   end do
   !! Just a check here, maybe forcing the mid point to be zero is brutal
   !! but by construction it should be zero !
   !if (dabs(vel(Nw+Nc-1)) <= 1d-7) vel(Nw+Nc-1) = 0d0
   if (dabs(vel(Nw+Nc-1)) /= 0d0) vel(Nw+Nc-1) = 0d0
   if (vel(Nw+Nc-1) /= 0) write(*,*) 'Vel(Nw+Nc-1) should be 0.0'

  !number of points from -vchar to 0 is Nw+Nc-1, -1 because otherwise we count
  ! 2 times vcore which is shared by the wing (upper boundary) and core (lower boundary) grid.
  ! Total number of points is 2*(Nw+Nc-1) but here we count 2 times lambda0., therefore
  ! we remove 1 point.
   Nlambda = 2 * (Nw + Nc - 1) - 1 
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
   double precision, intent(in) :: wl_ref
   integer, intent(out) :: Ntrans !Total number of transitions (cont + line)
   ! output grid. May contain values that are added to the final list before
   ! deallocating the array. It is reallocated when the final list is known.
   double precision, allocatable, dimension(:), intent(inout) :: inoutgrid
   ! temporary storage for transitions, to count and add them.
   integer, parameter :: MAX_TRANSITIONS = 50000
   type (AtomicLine), allocatable, dimension(:) :: alllines
   type (AtomicContinuum), allocatable, dimension(:) :: allcont
   integer :: kr, kc, n, Nspect, Nwaves, Nlinetot, Nctot
   integer :: la, nn, nnn !counters: wavelength, number of wavelengths, line wavelength
   integer :: Nred, Nblue, Nlambda_original!read from model
   double precision, allocatable, dimension(:) :: tempgrid, sorted_indexes
   double precision :: l0, l1 !ref wavelength of each transitions

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
    CALL fillPhotoionisationCrossSection(Atoms(n)%ptr_atom, kc, &
    	Atoms(n)%ptr_atom%continua(kc)%lambda,Nwaves, inoutgrid)
    	
    !TEST 
    Atoms(n)%ptr_atom%continua(kc)%Nlambda = size( Atoms(n)%ptr_atom%continua(kc)%alpha)

    	
    	
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
!--> Futur deprecation of this. Two memory consuming  
  SUBROUTINE fillPhotoionisationCrossSection(atom, kc, old_waves, Nwaves, waves)
   !computes photoionisation Xsections on the waves grid.
   !Eventually interpolate the values of the old grid old_waves onto waves
   use atmos_type, only : Hydrogen
   use math, only : bezier3_interp, Gaunt_bf
    type(AtomType), intent(inout) :: atom
    integer, intent(in) :: kc, Nwaves
    double precision, dimension(Nwaves) :: waves
    integer :: i, j, Nblue, Nred
    double precision, dimension(:), intent(in) :: old_waves
    double precision, dimension(Nwaves) :: uu, g_bf
    double precision, dimension(:), allocatable :: old_alpha
    double precision :: n_eff, Z, gbf_0(1), uu0(1), lambdaEdge, sigma0
    
    !pour H only?
    sigma0 = (32.)/(PI*3.*dsqrt(3d0)) * EPSILON_0 * &
          (HPLANCK**(3d0)) / (CLIGHT * &
          (M_ELECTRON*Q_ELECTRON)**(2d0))
    
    i = atom%continua(kc)%i; j = atom%continua(kc)%j
    Z = real(atom%stage(i)+1); lambdaEdge = atom%continua(kc)%lambda0
    Nblue = atom%continua(kc)%Nblue; Nred = atom%continua(kc)%Nred
    if (Nblue==-99 .and. Nred==-99) then
     !if (allocated(atom%continua(kc)%lambda)) deallocate(atom%continua(kc)%lambda)
     if (allocated(atom%continua(kc)%alpha)) deallocate(atom%continua(kc)%alpha)
     RETURN !no need to compute
    end if

    if (atom%continua(kc)%hydrogenic) then
     
     if (atom%ID == "H ") then
      n_eff = dsqrt(Hydrogen%g(i)/2.)  !only for Hydrogen !
     else
     !obtained_n = getPrincipal(metal%label(continuum%i), n_eff)
     !if (.not.obtained_n) &
        n_eff = Z*dsqrt(E_RYDBERG / (atom%E(j) - atom%E(i)))
     end if
     if (allocated(atom%continua(kc)%alpha)) deallocate(atom%continua(kc)%alpha)
     allocate(atom%continua(kc)%alpha(Nwaves)) !it is allocated in readatom only for non
     									   ! hydrogenic continua
     atom%continua(kc)%alpha(:) = 0d0
     g_bf = 0d0
     uu(Nblue:Nred) = n_eff*n_eff*HPLANCK*CLIGHT/ & 
        (NM_TO_M*waves(Nblue:Nred)) / (Z*Z) / E_RYDBERG - 1.
     uu0 = n_eff*n_eff * HPLANCK*CLIGHT / (NM_TO_M * lambdaEdge) / Z / Z / E_RYDBERG - 1.
       
     gbf_0 = Gaunt_bf(1, uu0, n_eff)
    
     g_bf(Nblue:Nred) = &
      Gaunt_bf(Nred-Nblue+1, uu(Nblue:Nred), n_eff)
      
     atom%continua(kc)%alpha(Nblue:Nred) = &
        atom%continua(kc)%alpha0 * g_bf(Nblue:Nred) * (waves(Nblue:Nred)/lambdaEdge)**3  / gbf_0(1)!m^2

    else !cont%alpha is allocated and filled with read values
     !!CALL Warning(" Beware, memory error if tow many cross-sections")
     write(*,*) "Interpolating photoionisation cross-section for atom ", &
      j,'->',i,atom%continua(kc)%atom%ID, " (", atom%ID,')'
     
     !that is because in this case, %alpha and %lambda have their size taken from the atomic
     !file with is not %Nlambda anymore at this point.
     allocate(old_alpha(size(old_waves)))
     old_alpha(:) = atom%continua(kc)%alpha(:)
!      if (size(waves)==1) then
!      write(*,*) size(atom%continua(kc)%alpha(:)), size(old_waves)
!      write(*,*) size(old_alpha)
!      write(*,*) waves(Nblue:Nred), atom%continua(kc)%Nlambda
!      end if
     deallocate(atom%continua(kc)%alpha)
     allocate(atom%continua(kc)%alpha(Nwaves))

     atom%continua(kc)%alpha(:) = 0d0
     CALL bezier3_interp(size(old_alpha),old_waves, old_alpha, & !Now the interpolation grid
             atom%continua(kc)%Nlambda,  waves(Nblue:Nred), &
             atom%continua(kc)%alpha(Nblue:Nred)) !end
     deallocate(old_alpha)
    end if
    if (allocated(atom%continua(kc)%lambda)) deallocate(atom%continua(kc)%lambda) !not used
    !!if (atom%ID=="He") read(*,*)
  RETURN
  END SUBROUTINE fillPhotoionisationCrossSection

  SUBROUTINE adjust_wavelength_grid(old_grid, lambda, Lam_region, Atoms)
   ! ------------------------------------------ !
    ! Reallocate wavelengths and indexes arrays
    ! to compute images on a user defined grid
   ! ------------------------------------------ !
   use math, only : locate
   use atmos_type, only : realloc_continuum_Transitions, realloc_line_Transitions
   double precision, dimension(:), intent(in) :: old_grid
   integer, dimension(:), intent(in) :: lam_region
   type (atomPointerArray), dimension(:), intent(inout) :: Atoms
   double precision, dimension(:), intent(inout) :: lambda
   double precision, dimension(size(lambda)) :: lambda_us
   integer :: Nwaves, n, kr, kc, Nlambda_original, Nblue, Nred, Natom, ll, lll
   double precision :: l0, l1 !ref wavelength of each transitions
   double precision :: x0, x1 !bound of the new grid
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
    allocate(trans_contribute(Atoms(n)%ptr_atom%Ncont)); trans_contribute(:)=.true.
    do kc=1,Atoms(n)%ptr_atom%Ncont
     !on the old_Grid (i.e., the grid for NLTE which is built using all transitions)
     Nlambda_original = Atoms(n)%ptr_atom%continua(kc)%Nlambda !on the old_grid
     Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue
     Nred = Atoms(n)%ptr_atom%continua(kc)%Nred
     l0 = Atoms(n)%ptr_atom%continua(kc)%lambdamin!old_grid(Nblue)
     l1 = Atoms(n)%ptr_atom%continua(kc)%lambda0!old_grid(Nred)
     !!if l0 out of the new grid, remove it, otherwise get the new index
     !if (l1 < minval(lambda)) then
!      if (l1 <= x0.or. l0 >= x1) then
!       !if l1 < x0 then automatically l0 < x0;
! 									!because l0<l1.
!      								!and conversely, if l0 > x1 -> l1 > x1
!       Atoms(n)%ptr_atom%continua(kc)%Nred = -99
!       Atoms(n)%ptr_atom%continua(kc)%Nblue = -99
!      else
!       !Otherwise, if l1 in )x0, x1] and l0 < x0 then l0==x0
!       !and if l0 in [x0, x1( and l1 > x1 then l1 == x1
!       Atoms(n)%ptr_atom%continua(kc)%Nred = locate(lambda,l1) ! closest value return by locate is Nwaves if l1>lambda(Nwaves)
!       Atoms(n)%ptr_atom%continua(kc)%Nblue = locate(lambda,l0) !
!      end if
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
!      if (l0 > maxval(lambda)) then
!       Atoms(n)%ptr_atom%continua(kc)%Nblue = -99
!      else
!       Atoms(n)%ptr_atom%continua(kc)%Nblue = locate(lambda,l0)
!      end if
!      
!      if ((Atoms(n)%ptr_atom%continua(kc)%Nblue==-99).and.(Atoms(n)%ptr_atom%continua(kc)%Nred/=-99)) then
!       Atoms(n)%ptr_atom%continua(kc)%Nred=-99
!      end if
!      if ((Atoms(n)%ptr_atom%continua(kc)%Nblue/=-99).and.(Atoms(n)%ptr_atom%continua(kc)%Nred==-99)) then
!       Atoms(n)%ptr_atom%continua(kc)%Nblue=-99
!      end if
     
     Nblue = Atoms(n)%ptr_atom%continua(kc)%Nblue; Nred = Atoms(n)%ptr_atom%continua(kc)%Nred
     if (Nred==-99.and.Nblue==-99) then
      Atoms(n)%ptr_atom%continua(kc)%Nlambda = -99

      Atoms(n)%ptr_atom%continua(kc)%Nmid = -99
      
      trans_contribute(kc) = .false.
      write(*,*) " :: b-f transition", Atoms(n)%ptr_atom%continua(kc)%j,"->",Atoms(n)%ptr_atom%continua(kc)%i,&
       " for atom ",Atoms(n)%ptr_atom%ID, l0,"-",l1, " removed."
     else
      Atoms(n)%ptr_atom%continua(kc)%Nlambda = Atoms(n)%ptr_atom%continua(kc)%Nred - &
                                 Atoms(n)%ptr_atom%continua(kc)%Nblue + 1

      Atoms(n)%ptr_atom%continua(kc)%Nmid = locate(lambda,lambda(Atoms(n)%ptr_atom%continua(kc)%Nlambda)/2+1)
     !cont%alpha presently defined on old_lambda, and then on new lambda lambda.
     !cont%lambda deallocated
      CALL fillPhotoionisationCrossSection(Atoms(n)%ptr_atom, kc, old_grid, Nwaves, lambda)
      !!-> also deallocates cont%lambda for hydrogenic atoms
      !! for all transitions even if we removed them?
     end if                          
    end do
    !If the transition is removed no need for testing Nred and Nlblue in the opacities
    CALL realloc_continuum_transitions(Atoms(n)%ptr_atom, count(trans_contribute), trans_contribute)
    !write(*,*) maxval(Atoms(n)%ptr_atom%continua(1)%alpha)
    deallocate(trans_contribute)
    !then bound-bound transitions
    allocate(trans_contribute(Atoms(n)%ptr_atom%Nline)); trans_contribute(:)=.true.
    do kr=1,Atoms(n)%ptr_atom%Nline
     Nlambda_original = Atoms(n)%ptr_atom%lines(kr)%Nlambda !on the old_grid
     Nblue = Atoms(n)%ptr_atom%lines(kr)%Nblue
     Nred = Atoms(n)%ptr_atom%lines(kr)%Nred
     l0 = old_grid(Nblue)
     l1 = old_grid(Nred)
     !if l0 out of the new grid, remove it, otherwise get the new index
     !if (l0 > maxval(lambda)) then
     !relative index of regions
!      if (l1 <= x0 .or. l0 >= x1) then
!       Atoms(n)%ptr_atom%lines(kr)%Nblue = -99
!       Atoms(n)%ptr_atom%lines(kr)%Nred = -99
!      else
!       Atoms(n)%ptr_atom%lines(kr)%Nblue = locate(lambda,l0)
!       Atoms(n)%ptr_atom%lines(kr)%Nred = locate(lambda,l1)
!       exit region_loop_l !because if in one region no need to test the others
!      end if
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
!      if (l1 < minval(lambda)) then
!       Atoms(n)%ptr_atom%lines(kr)%Nred = -99
!      else
!       Atoms(n)%ptr_atom%lines(kr)%Nred = locate(lambda,l1)
!      end if
!      
!      if ((Atoms(n)%ptr_atom%lines(kr)%Nblue==-99).and.(Atoms(n)%ptr_atom%lines(kr)%Nred/=-99)) then
!       Atoms(n)%ptr_atom%lines(kr)%Nred=-99
!      end if
!      if ((Atoms(n)%ptr_atom%lines(kr)%Nblue/=-99).and.(Atoms(n)%ptr_atom%lines(kr)%Nred==-99)) then
!       Atoms(n)%ptr_atom%lines(kr)%Nblue=-99
!      end if
     
     Nblue = Atoms(n)%ptr_atom%lines(kr)%Nblue; Nred = Atoms(n)%ptr_atom%lines(kr)%Nred
     if (Nred==-99.and.Nblue==-99) then
      Atoms(n)%ptr_atom%lines(kr)%Nlambda = -99

      Atoms(n)%ptr_atom%lines(kr)%Nmid = -99

      trans_contribute(kr) = .false.
      write(*,*) " :: b-b transition", Atoms(n)%ptr_atom%lines(kr)%j,"->",Atoms(n)%ptr_atom%lines(kr)%i,&
       " for atom " ,Atoms(n)%ptr_atom%ID, l0,"-",l1, " removed."
     else
      Atoms(n)%ptr_atom%lines(kr)%Nlambda = Atoms(n)%ptr_atom%lines(kr)%Nred - &
                                 Atoms(n)%ptr_atom%lines(kr)%Nblue + 1

      Atoms(n)%ptr_atom%lines(kr)%Nmid = locate(lambda,Atoms(n)%ptr_atom%lines(kr)%lambda0)
     end if                               
    if (allocated(Atoms(n)%ptr_atom%lines(kr)%lambda)) &
    	deallocate(Atoms(n)%ptr_atom%lines(kr)%lambda)
    end do
!      Atoms(n)%ptr_atom%lines = PACK(Atoms(n)%ptr_atom%lines, trans_contribute)
!      Atoms(n)%ptr_atom%Nline = count(trans_contribute)
     CALL realloc_line_transitions(Atoms(n)%ptr_atom, count(trans_contribute), trans_contribute)
     deallocate(trans_contribute)
   end do !over atoms
   if (allocated(trans_contribute)) deallocate(trans_contribute) !in case

  RETURN
  END SUBROUTINE adjust_wavelength_grid

  END MODULE getlambda
