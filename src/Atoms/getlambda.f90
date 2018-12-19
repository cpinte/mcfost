MODULE getlambda

  use atom_type, only : AtomicContinuum, AtomicLine, AtomType
  use atmos_type, only : atmos
  use constant
  
  IMPLICIT NONE

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
   integer, parameter :: Nlambda = 25
   integer :: la
   double precision :: l0, l1

   l1 = cont%lambda0 !cannot be larger than lambda0 ! minimum frequency for photoionisation
   l0 = lambdamin
   cont%Nlambda = Nlambda
   allocate(cont%lambda(cont%Nlambda))
   resol = (l1 - l0) / (Nlambda - 1)
   write(*,*) "Continuum:", cont%lambda0, cont%j,"->",cont%i, &
              " Resolution (nm):", resol, " lambdamin =", lambdamin  
   
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
  ! ------------------------------------------------------------ !
   type (AtomicLine), intent(inout) :: line
   double precision, intent(in) :: vD !maximum thermal width of the atom in m/s
   double precision :: v_char, dvc, dvw
   double precision :: vcore, v0, v1!km/s
   integer :: la, Nlambda, Nmid
   double precision, parameter :: core_to_wing=5d-1, L=7d0
   integer, parameter :: Nc = 101, Nw = 11 !ntotal = 2*(Nc + Nw - 1)
   double precision, dimension(5*(Nc+Nw)) :: vel
   
   v_char = L * (atmos%v_char + vD)
   v0 = -v_char
   v1 = +v_char
   vel = 0d0
   vcore = 2d0 * core_to_wing * v_char /L !converge to 2/L*core_to_wing*vD in case of static atm.
   
   !from -v_char to 0
   dvw = (v_char-vcore)/(Nw-1)
   dvc = vcore/(Nc-1)
   write(*,*) "line:", line%lambda0,line%j,"->",line%i, &
              " Resolution wing,core (km/s):", dvw/1d3,dvc/1d3
   vel(1) = v0 !half wing
!   write(*,*) 1, vel(1), v0, v_char/1d3
   !wing loop
   do la=2, Nw
    vel(la) = vel(la-1) + dvw
!    write(*,*) la, vel(la)
   end do
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
   if (dabs(vel(Nw+Nc-1)) <= 1d-7) vel(Nw+Nc-1) = 0d0
   if (vel(Nw+Nc-1) /= 0) write(*,*) 'Vel(Nw+Nc-1) should be 0.0'
!stop
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
   
!    write(*,*) 1, line%lambda(1)
!    do la=2,Nlambda
!     write(*,*) la, line%lambda(la), line%lambda0
!    end do
!    stop
  RETURN
  END SUBROUTINE make_sub_wavelength_grid_line
  
  !building
  SUBROUTINE make_sub_wavelength_grid_asymm(line, vD)
  ! ------------------------------------------------------- !
   ! Make an individual wavelength grid for the AtomicLine
   ! line with constant resolution in km/s. When dlambda 
   ! is lower than vsub, the resolution is increased.
   ! Allocate line%lambda.
  ! ------------------------------------------------------- !
   type (AtomicLine), intent(inout) :: line
   double precision, intent(in) :: vD !maximum thermal width of the atom in m
   double precision :: v_char
   double precision :: subresol, vsub, resol!km/s
   integer :: la, Nlambda, Nmid
   integer, parameter :: NlambdaMax = 30000
   double precision, parameter :: L = 2d0, core_to_wing = 0.3d0, R = 4.6
   double precision :: lam_grid(NlambdaMax), dlambda, l0, l1 !nm
   !the line domain falls in [l0*(1+resol), l1*(1+resol)] with, 
   !l0 = (lambda0 - lambda0*v_char/c - vWing) and l1 = (lambda0 + lambda0*v_char/c + vWing).
   v_char = L * (atmos%v_char + vD + R*1d3)
   resol = R !km/s
   dlambda = v_char /CLIGHT * line%lambda0
   vsub = core_to_wing * (R + vD * 1d-3)!km/s
   subresol = 1d-2 * resol
   l1 = (line%lambda0 + dlambda) * (1 + 1d3 * resol/CLIGHT)
   l0 = (line%lambda0 - dlambda) * (1 - 1d3 * resol/CLIGHT)
   !symmetric boundaries
!    write(*,*) v_char/1d3, dlambda, l1, l0, vsub, subresol
!    stop
   
   Nlambda = 2
   lam_grid(1) = line%lambda0
   infinie : do ! from lambda0 to l1
    lam_grid(Nlambda) = lam_grid(Nlambda - 1) * (1d0 + 1d3 * resol / CLIGHT)
    ! increase resolution in the line core
    if (CLIGHT * dabs(lam_grid(Nlambda) - line%lambda0)/line%lambda0 <= 1d3 * vsub) then
     sub_infinie : do
     lam_grid(Nlambda) = lam_grid(Nlambda-1) * (1d0 + 1d3 * subresol/CLIGHT)

     if ((Nlambda > NlambdaMax).or.(lam_grid(Nlambda)>l1)) exit infinie
     if (CLIGHT * dabs(lam_grid(Nlambda) - line%lambda0)/line%lambda0 > 1d3 * vsub) exit sub_infinie
     write(*,*) Nlambda, lam_grid(Nlambda)
     Nlambda = Nlambda + 1
     end do sub_infinie
    end if
     write(*,*) Nlambda, lam_grid(Nlambda)

    if ((Nlambda > NlambdaMax).or.(lam_grid(Nlambda)>l1)) exit infinie
    Nlambda = Nlambda + 1
   end do infinie

!    stop
   
   Nlambda = 2 * Nlambda
   line%Nlambda = Nlambda
   Nmid = Nlambda/2
   allocate(line%lambda(line%Nlambda))   
   do la=Nmid, Nlambda
    line%lambda(la) = lam_grid(la)
    line%lambda(Nlambda - la + 1) = lam_grid(la) - line%lambda0
   end do
   do la=1,Nlambda
    write(*,*) la, line%lambda(la)
   end do
   stop
   
  RETURN
  END SUBROUTINE make_sub_wavelength_grid_asymm

  SUBROUTINE getLambdaCont(continuum, lambdamin)
   type (AtomicContinuum), intent(inout) :: continuum
   double precision, intent(in) :: lambdamin
   integer :: la, Nlambda
   double precision :: dlambda
   ! from lambdamin to lambdaEdge or lambda ionisation
   ! or from maximum frequency to frequency threshold.

   ! Nlambda set in the atomic file
   allocate(continuum%lambda(continuum%Nlambda))
   !allocate(continuum%alpha(continuum%Nlambda))   !not used in this case
   
   Nlambda =  continuum%Nlambda
   dlambda = (continuum%lambda0 - lambdamin) / (Nlambda-1)
   continuum%lambda(1) = lambdamin
   do la=2,Nlambda
    continuum%lambda(la) = continuum%lambda(la-1) + dlambda
   end do

  RETURN
  END SUBROUTINE getLambdaCont
  
  SUBROUTINE getLambdaLine(line)
  use math, only : SQ, CUBE, locate, sort
  ! Adapted from RH H. Uitenbroek
  ! Nlambda, qcore and qwing read from file.
  ! --- Construct a wavelength grid that is approximately equidistant
  !in the core (q <= qcore) and equidistant in log(q) in the wings
  !(qcore < q <= qwing).
  !
  !Look for a function of the form: q[n] = a*(n + (exp(b*n)-1));
  !n=0, N-1 satisfying the following conditions:
  !
  ! -- q[0]   = 0          (this is true for all a and b)
  ! -- q[(N-1)/2] = qcore
  ! -- q[N-1] = qwing
   type (AtomicLine), intent(inout) :: line
   integer :: la, n, Nlambda, NB=0, Nmid
   real(8) :: beta, a, b, y, q_to_lambda, qB_char, qB_shift,dlambda
   real(8), dimension(:), allocatable :: q
   real(8) :: g_lande_eff, lambda0

   if ((line%qcore.le.0.).or.(line%qwing.le.0.)) then
    write(*,*) "Either qcore or qwing <= 0, for line ", &
        line%j,"->",line%i
    write(*,*) "exiting..."
    stop
   end if
   
   !write(*,*) line%j, "->", line%i, line%qcore, line%qwing

   !store in unit of lambda0 instead of adimensional values
   !deltaLambda / lambda0 = deltav/c
   if (atmos%v_char <= 0.) then
   write(*,*) "Error with getlambda, set v_char /= 0"
   stop
   end if
   q_to_lambda = line%lambda0 * (atmos%v_char/CLIGHT)

   g_lande_eff = line%g_lande_eff

   !--- Create full (instead of half) line profiles in case of:

  !1) a set keyword ASYMM in this_atom.atom for this line
  !2) a moving grid wich is the case for accretion disk
  !3) a magnetic field is present and the g_Lande is
  !greater than -99 for this line
  !4) the line has more than one components

  ! Because atmos will be often moving, line
  ! will be always asymmetric
  if ((atmos%Magnetized .and. line%polarizable) .or. &
     atmos%moving) then
     line%symmetric=.false.
  end if

  ! Make sure the half-line profile always has an odd number of
  ! grid points

  if (line%symmetric) then
   if (MOD(line%Nlambda,2) >= 1.) then
     Nlambda=line%Nlambda !odd
   else
     Nlambda = line%Nlambda+1 !result of MOD is 0, even
   end if
  else
   if (MOD(line%Nlambda,2) >= 1.) then
     Nlambda = line%Nlambda/2 !odd
   else
     Nlambda = (line%Nlambda+1)/2 !result of MOD is 0, even
   end if
  end if
  !write(*,*) "Nlambda, line%Nlambda = ",Nlambda, line%Nlambda

  if (line%qwing <= 2.*line%qcore) then
   !linear scale out to qwing
   write(*,*) "Ratio of qwing/qcore < 1 for line ", &
          line%j,"->",line%i
   beta = 1
  else
   beta = line%qwing/(2.*line%qcore)
  end if
  y = beta+sqrt(SQ(beta) + (beta-1)*Nlambda + 2.-3.*beta)
  b = 2.*log(y) / (Nlambda-1)
  a = line%qwing / (Nlambda-2. + SQ(y))

  allocate(q(Nlambda))
  do la=1,Nlambda
    q(la) = a*((la-1)+(exp(b*(la-1))-1))
    !write(*,*) q(la)*q_to_lambda, " nm"
  end do
  if (line%polarizable) then
      write(*,*) "line is polarizable and B is present"
      !characteristic Zeeman splitting
      qB_char = g_lande_eff * (Q_ELECTRON/(4*PI*M_ELECTRON))*&
           (line%lambda0*NM_TO_M) * atmos%B_char / &
            atmos%v_char

    if (qB_char/2. >= line%qcore) then
      write(*,*) "Characteristic Zeeman Splitting >= qcore ",&
      "for line ", line%j,"->",line%i
    end if
     NB = locate(q, qB_char/2.)
     qB_shift = 2.*q(NB)
     !!CALL realloc(q, Nlambda+2*NB)
     deallocate(q)
     allocate(q(Nlambda+2*NB))

     do la=NB+1,2*NB+1
       q(la) = qB_shift - a*(2*NB-la-1 + (exp(b*(2*NB-la-1)) &
           - 1.))
     end do

     do la=2*NB+1,Nlambda+2*NB
       q(la) = qB_shift + a*(la-1 -2*NB+(exp(b*(-2*NB+la-1)) &
           - 1.))
     end do

     Nlambda = Nlambda+2*NB

   end if !over polarized line

   !beware because the grid is most of the time moving
   !accretion disk, stellar hydro models etc..
   !we rarely inter in the case line%symmetric
   !even if the line is defined symmetric in the model atom.
   if (line%symmetric) then
    line%Nlambda = Nlambda
    allocate(line%lambda(Nlambda))
    do la=1,Nlambda
     line%lambda(la) = line%lambda0+q_to_lambda*q(la)
    end do
   else !not symmetric
    !write(*,*) "Setting up grid for asymmetric line"
    ! minus - 1 to force it to be odd
    line%Nlambda = (2*Nlambda-1)! * line%Ncomponent
    !!write(*,*) Nlambda, 2*Nlambda-1, line%Nlambda
    allocate(line%lambda(line%Nlambda))
    !do n=1,line%Ncomponent
    n = 1
     !in c, n=0,Ncompo - 1
     Nmid = (n-1)*(2*Nlambda-1) + Nlambda -1
     lambda0 = line%lambda0! + line%c_shift(n)
     line%lambda(Nmid) = lambda0
      do la=1,Nlambda !la=Nlambda->
        ! Nlambda-1 +Nlambda = 2*Nlambda-1, last point
        dlambda = q_to_lambda*q(la)
        line%lambda(Nmid-la) = lambda0-dlambda
        line%lambda(Nmid+la) = lambda0+dlambda
      end do
      !do la=1,line%Nlambda
      ! write(*,*) line%lambda(la)
      !end do
    !end do

!     if (line%Ncomponent.gt.1) then
!      CALL sort(line%lambda, line%Nlambda)
!      !!CALL hpsort(line%Nlambda, line%lambda) ! NumRec
!     end if

   end if

  deallocate(q)
  RETURN
  END SUBROUTINE getLambdaLine

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
   !!-> why it is not working if I want to remove the dynamic allocation? ...
   allocate(alllines(MAX_TRANSITIONS), allcont(MAX_TRANSITIONS))
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
  write(*,*) "Adding ", Nspect," lines and continua, for a total of", Ntrans, &
             " transitions."

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
  write(*,*) Nwaves, " unique wavelengths, ", Nspect-Nwaves,&
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
  !deallocate line%lambda and continuum%lambda because they do not correspond to the new
  ! Nblue and Nlambda anymore.
  do n=1,Natom
   Nred = 0
   !first continuum transitions
!   write(*,*) " ------------------------------------------------------------------ "
   do kc=1,Atoms(n)%Ncont
    Nlambda_original = Atoms(n)%continua(kc)%Nlambda
    l0 = Atoms(n)%continua(kc)%lambda(1)
    l1 = Atoms(n)%continua(kc)%lambda(Nlambda_original)
    Nred = locate(inoutgrid,l1)
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
    Atoms(n)%continua(kc)%Nmid = locate(inoutgrid,Atoms(n)%continua(kc)%lambda0)
   end do
   !then bound-bound transitions
   do kr=1,Atoms(n)%Nline
    Nlambda_original = Atoms(n)%lines(kr)%Nlambda
    l0 = Atoms(n)%lines(kr)%lambda(1)
    l1 = Atoms(n)%lines(kr)%lambda(Nlambda_original)
    Nred = locate(inoutgrid,l1)
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
   end do
!     write(*,*) " ------------------------------------------------------------------ "
  end do !over atoms

  RETURN
  END SUBROUTINE make_wavelength_grid

  END MODULE getlambda
