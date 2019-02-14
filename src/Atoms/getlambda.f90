MODULE getlambda

  use atom_type, only : AtomicContinuum, AtomicLine, AtomType
  use atmos_type, only : atmos
  use constant
  
  use parametres
  use utils, only : span, spanl
  
  IMPLICIT NONE
  
  double precision, dimension(:), allocatable :: Nred_array, Nblue_array, Nmid_array

  CONTAINS
  
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
   integer, parameter :: Nlambda = 41
   integer :: la
   double precision :: l0, l1

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
   !do not allocate cross-section
   !allocate(cont%alpha(cont%Nlambda)) -> not used in this case
   ! necessary only if the explicit wavelength dependence is given
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
   double precision :: vcore, v0, v1!km/s
   integer :: la, Nlambda, Nmid
   double precision, parameter :: wing_to_core = 0.3, L = 10d0!50d0
   integer, parameter :: Nc = 71, Nw = 11 !ntotal = 2*(Nc + Nw - 1) - 1
   double precision, dimension(2*(Nc+Nw-1)-1) :: vel !Size should be 2*(Nc+Nw-1)-1
   													 !if error try, 2*(Nc+Nw)
   
   v_char = (atmos%v_char + vD) !=maximum extension of a line
   !atmos%v_char is minimum of Vfield and vD is minimum of atom%vbroad presently
   v0 = -v_char * L
   v1 = +v_char * L
   vel = 0d0
   !transition between wing and core in velocity
   !vcore = L * v_char * wing_to_core ! == fraction of line extent
   vcore = v_char * 1.5d0

   !from -v_char to 0
   dvw = (L * v_char-vcore)/real(Nw-1,kind=dp)
   dvc = vcore/real(Nc-1,kind=dp)
!    write(*,*) "line:", line%lambda0,line%j,"->",line%i, &
!               " Resolution wing,core (km/s):", dvw/1d3,dvc/1d3
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
   vel(1:Nw) = -real(spanl(real(v0), real(vcore), Nw),kind=dp)
!! end scale of wing points
!   vel(Nw:Nw-1+Nc) = -real(span(real(vcore), real(0.), Nc+1),kind=dp) 
!   write(*,*) Nw, vel(Nw), vcore
   !vel(Nw) = -vcore!should be okey at the numerical precision 
   !la goes from 1 to Nw + Nc -1 total number of points.
   !if Nc=101 and Nw =11, 111 velocity points,because the point vel(11) of the wing grid
   !is shared with the point vel(1) of the core grid.
   do la=Nw+1, Nc+Nw-1 !Half line core
    vel(la) = vel(la-1) + dvc
!    write(*,*) la-Nw+1,la, vel(la) !relative index of core grid starts at 2 because vel(Nw)
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
   
!    do la=1,Nlambda
!     write(*,*) la, line%lambda(la), line%lambda0
!    end do

  RETURN
  END SUBROUTINE make_sub_wavelength_grid_line
  
  FUNCTION IntegrationWeightLine(line, la) result (wlam)
  ! ------------------------------------------------------- !
   ! result in m/s
   ! A 1/vbroad term is missing.
  ! ------------------------------------------------------- !
   double precision :: wlam, Dopplerwidth
   type (AtomicLine) :: line
   integer :: la
   
   Dopplerwidth = CLIGHT/line%lambda0
   if (line%symmetric) Dopplerwidth = Dopplerwidth*2
   
   if (la == 1) then
    wlam = 0.5*(line%lambda(la+1)-line%lambda(la))
   else if (la == line%Nlambda) then
    wlam = 0.5*(line%lambda(la)-line%lambda(la-1))
   else
    wlam = 0.5*(line%lambda(la+1)-line%lambda(la-1))
   end if
   
   wlam = wlam * DopplerWidth

  RETURN
  END FUNCTION IntegrationWeightLine

  FUNCTION IntegrationWeightCont(cont, la) result(wlam)
  ! ------------------------------------------------------- !
   ! result in nm
  ! ------------------------------------------------------- !
   type (AtomicContinuum) :: cont
   integer :: la
   double precision :: wlam

   if (la == 1) then
    wlam = 0.5*(cont%lambda(la+1)-cont%lambda(la))
   else if (la == cont%Nlambda-1) then
    wlam = 0.5*(cont%lambda(la)-cont%lambda(la-1))
   else
    wlam = 0.5*(cont%lambda(la+1)-cont%lambda(la-1))
   end if

  RETURN
  END FUNCTION

  SUBROUTINE make_wavelength_grid(Natom, Atoms, inoutgrid, Ntrans, wl_ref)
  use math, only : locate, sort
  ! --------------------------------------------------------------------------- !
   ! construct and sort a wavelength grid for atomic line radiative transfer.
   ! Computes also the edge of a line profile: Nblue and Nred.
   ! The grid is built by merging the individual wavelength grids of each 
   ! transitions. Duplicates are removed and therefore changes the value of
   ! line%Nlambda and continuum%Nlambda.
   ! line%lambda and continuum%lambda are useless now. Except that 
   ! continuu%lambda is still used in .not.continuum%Hydrogenic !!
  ! --------------------------------------------------------------------------- !
   type (AtomType), intent(inout), dimension(Natom) :: Atoms
   integer, intent(in) :: Natom
   double precision, intent(in) :: wl_ref
   integer, intent(out) :: Ntrans !Total number of transitions (cont + line)
   ! output grid. May contain values that are added to the final list before
   ! deallocating the array. It is reallocated when the final list is known.
   double precision, allocatable, dimension(:), intent(inout) :: inoutgrid
   ! temporary storage for transitions, to count and add them.
   integer, parameter :: MAX_TRANSITIONS = 10000
   type (AtomicLine), allocatable, dimension(:) :: alllines
   type (AtomicContinuum), allocatable, dimension(:) :: allcont
   integer :: kr, kc, n, Nspect, Nwaves, Nlinetot, Nctot
   integer :: la, nn, nnn !counters: wavelength, number of wavelengths, line wavelength
   integer :: Nred, Nblue, Nlambda_original!read from model
   double precision, allocatable, dimension(:) :: tempgrid
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
   do kr=1,Atoms(n)%Nline
    Nspect = Nspect + Atoms(n)%lines(kr)%Nlambda
    Nlinetot = Nlinetot + 1
    if (Nlinetot > MAX_TRANSITIONS) then
     write(*,*) "too many transitions"
     stop
    end if
    alllines(Nlinetot) = Atoms(n)%lines(kr)
   end do
   do kc=1,Atoms(n)%Ncont
    Nspect = Nspect + Atoms(n)%continua(kc)%Nlambda
    Nctot = Nctot + 1
    if (Nctot > MAX_TRANSITIONS) then
     write(*,*) "too many transitions"
     stop
    end if
    allcont(Nctot) = Atoms(n)%continua(kc)
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

  allocate(Nred_array(Ntrans), Nblue_array(Ntrans), Nmid_array(Ntrans))

  ! add wavelength from mcfost inoutgrid if any
  ! and convert it to nm, then deallocate
  if (allocated(inoutgrid)) then
   do la=1, size(inoutgrid)
    nn = nn + 1
    tempgrid(nn) = inoutgrid(la) * 1d3 !micron to nm
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
  CALL sort(tempgrid, Nspect)

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

  !should not dot that but error somewhere if many atoms
  !CALL sort(inoutgrid, Nwaves)

  !free some space
  deallocate(tempgrid, alllines, allcont)

  ! Now replace the line%Nlambda and continuum%Nlambda by the new values.
  ! we do that even for PASSIVE atoms
  ! deallocate line%lambda and continuum%lambda because they do not correspond to the new
  ! Nblue and Nlambda anymore.
  nn = 1
  do n=1,Natom
   !first continuum transitions
!   write(*,*) " ------------------------------------------------------------------ "
   do kc=1,Atoms(n)%Ncont
    Nlambda_original = Atoms(n)%continua(kc)%Nlambda
    l0 = Atoms(n)%continua(kc)%lambda(1)
    l1 = Atoms(n)%continua(kc)%lambda(Nlambda_original)
!    Nred = locate(inoutgrid,l1)
!     Nblue = locate(inoutgrid,l0)
!     write(*,*) locate(inoutgrid,l0), locate(inoutgrid,l1), l0, l1!, Nblue, Nred
    Atoms(n)%continua(kc)%Nblue = locate(inoutgrid,l0)
    Atoms(n)%continua(kc)%Nred = locate(inoutgrid,l1)
    Atoms(n)%continua(kc)%Nlambda = Atoms(n)%continua(kc)%Nred - &
                                    Atoms(n)%continua(kc)%Nblue + 1
!     write(*,*) Atoms(n)%ID, " continuum:",kr, " Nlam_ori:", Nlambda_original, &
!     " l0:", l0, " l1:", l1, " Nred:",  Atoms(n)%continua(kc)%Nred, & 
!       " Nblue:", Atoms(n)%continua(kc)%Nblue, " Nlambda:", Atoms(n)%continua(kc)%Nlambda, & 
!       " Nblue+Nlambda-1:", Atoms(n)%continua(kc)%Nblue + Atoms(n)%continua(kc)%Nlambda - 1
!! For continuum transitions, lambda0 is at Nred, check definition of the wavelength grid
!! which means that cont%Nmid = locate(inoutgrid, lam(Nred)+lam(Nb)/(Nlambda))
!! and l1, lam(Nlambda) = lambda0
    Atoms(n)%continua(kc)%Nmid = locate(inoutgrid,0.5*(l0+l1))!locate(inoutgrid,Atoms(n)%continua(kc)%lambda0)
    deallocate(Atoms(n)%continua(kc)%lambda)
    !allocate(Atoms(n)%continua(kc)%lambda(Atoms(n)%continua(kc)%Nlambda))
    !Atoms(n)%continua(kc)%lambda(Atoms(n)%continua(kc)%Nblue:Atoms(n)%continua(kc)%Nred) &
    ! = inoutgrid(Atoms(n)%continua(kc)%Nblue:Atoms(n)%continua(kc)%Nred)
    Nred_array(nn) = Atoms(n)%continua(kc)%Nred
    Nmid_array(nn) = Atoms(n)%continua(kc)%Nmid
    Nblue_array(nn) = Atoms(n)%continua(kc)%Nblue
    nn= nn + 1
   end do
   !then bound-bound transitions
   do kr=1,Atoms(n)%Nline
    Nlambda_original = Atoms(n)%lines(kr)%Nlambda
    l0 = Atoms(n)%lines(kr)%lambda(1)
    l1 = Atoms(n)%lines(kr)%lambda(Nlambda_original)
!    Nred = locate(inoutgrid,l1)
!    Nblue = locate(inoutgrid,l0)
!     write(*,*) locate(inoutgrid,l0), locate(inoutgrid,l1), l0, l1!, Nblue, Nred
    Atoms(n)%lines(kr)%Nblue = locate(inoutgrid,l0)!Nblue
    Atoms(n)%lines(kr)%Nred = locate(inoutgrid,l1)
    Atoms(n)%lines(kr)%Nlambda = Atoms(n)%lines(kr)%Nred - &
                                 Atoms(n)%lines(kr)%Nblue + 1
!     write(*,*) Atoms(n)%ID, " line:",kr, " Nlam_ori:", Nlambda_original, &
!     " l0:", l0, " l1:", l1, " Nred:",  Atoms(n)%lines(kr)%Nred, & 
!       " Nblue:", Atoms(n)%lines(kr)%Nblue, " Nlambda:", Atoms(n)%lines(kr)%Nlambda, & 
!       " Nblue+Nlambda-1:", Atoms(n)%lines(kr)%Nblue + Atoms(n)%lines(kr)%Nlambda - 1
    Atoms(n)%lines(kr)%Nmid = locate(inoutgrid,Atoms(n)%lines(kr)%lambda0)
    deallocate(Atoms(n)%lines(kr)%lambda) !does not correpond to the new grid, indexes might be wrong
    !allocate(Atoms(n)%lines(kr)%lambda(Atoms(n)%lines(kr)%Nlambda))
    !Atoms(n)%lines(kr)%lambda(Atoms(n)%lines(kr)%Nblue:Atoms(n)%lines(kr)%Nred) &
    ! = inoutgrid(Atoms(n)%lines(kr)%Nblue:Atoms(n)%lines(kr)%Nred)
    Nred_array(nn) = Atoms(n)%lines(kr)%Nred
    Nmid_array(nn) = Atoms(n)%lines(kr)%Nmid
    Nblue_array(nn) = Atoms(n)%lines(kr)%Nblue
    nn = nn + 1
   end do
!     write(*,*) " ------------------------------------------------------------------ "
  end do !over atoms

  RETURN
  END SUBROUTINE make_wavelength_grid

  END MODULE getlambda
