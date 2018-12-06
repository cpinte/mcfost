MODULE getlambda

  use atom_type, only : AtomicContinuum, AtomicLine, AtomType
  use atmos_type, only : atmos
  use constant
  use math, only : SQ, CUBE, locate, sort

  IMPLICIT NONE

  integer, parameter :: MAX_ACTIVE_TRANSITIONS = 1000 !maximum number of ACTIVE transitions


  CONTAINS


  SUBROUTINE getLambdaCont(continuum, lambdamin)
   type (AtomicContinuum), intent(inout) :: continuum
   double precision, intent(in) :: lambdamin
   integer :: la, Nlambda
   double precision :: dlambda
   ! from lambdamin to lambdaEdge or lambda ionisation
   ! or from maximum frequency to frequency threshold.
   Nlambda =  continuum%Nlambda
   dlambda = (continuum%lambda0 - lambdamin) / (Nlambda-1)
   continuum%lambda(1) = lambdamin
   do la=2,Nlambda
    continuum%lambda(la) = continuum%lambda(la-1) + dlambda
   end do

  RETURN
  END SUBROUTINE getLambdaCont
  
  SUBROUTINE make_sub_wavelength_grid(line, vWing)
   type (AtomicLine), intent(inout) :: line
   double precision, intent(in) :: vWing
   double precision :: v_char
   double precision, parameter :: resol = 1.78888d0, subresol = 0.2, vsub = 4.6 !km/s
   integer :: la, Nlambda
   integer, parameter :: NlambdaMax = 10000
   double precision :: lam_grid(NlambdaMax), dlambda, l0, l1
   
   v_char = max(atmos%v_char, vWing)
   
   Nlambda = 2
   dlambda = v_char/CLIGHT * line%lambda0 !v_char / CLIGHT * line%lambda0
   l1 = (line%lambda0 + dlambda) * (1 + 1d3 * resol/CLIGHT)
   l0 = (line%lambda0 - dlambda) * (1 - 1d3 * resol/CLIGHT)
   lam_grid(1) = l0
    infinie : do !from l0 to l1
    lam_grid(Nlambda) = lam_grid(Nlambda-1) * (1d0 + 1d3 * resol / CLIGHT)
    ! increase resolution in the line core
    if (CLIGHT * dabs(lam_grid(Nlambda) - line%lambda0)/line%lambda0 <= 1d3 * vsub) then
      sub_infinie : do
       lam_grid(Nlambda) = lam_grid(Nlambda-1) * (1d0 + 1d3 * subresol/CLIGHT)
       if (CLIGHT * dabs(lam_grid(Nlambda) - line%lambda0)/line%lambda0 > 1d3 * vsub) & 
          exit sub_infinie
       Nlambda = Nlambda + 1
       if (Nlambda > NlambdaMax) exit !infinie exit all
      end do sub_infinie
    end if
     if (lam_grid(Nlambda) >= l1) exit infinie
     Nlambda = Nlambda + 1
     if (Nlambda > NlambdaMax) exit infinie
   end do infinie

   line%Nlambda = Nlambda
   allocate(line%lambda(line%Nlambda))
   line%lambda(1:line%Nlambda) = lam_grid(1:line%Nlambda) 
!     do la=1,line%Nlambda
!     write(*,*) la, line%lambda(la)
!     end do
  RETURN
  END SUBROUTINE make_sub_wavelength_grid

  SUBROUTINE getLambdaLine(line)
! Adapted from RH H. Uitenbroek

 ! --- Construct a wavelength grid that is approximately equidistant
 !in the core (q <= qcore) and equidistant in log(q) in the wings
 !(qcore < q <= qwing).

 !Look for a function of the form: q[n] = a*(n + (exp(b*n)-1));
 !n=0, N-1 satisfying the following conditions:

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
   if (MOD(line%Nlambda,2).ge.1.) then
     Nlambda=line%Nlambda !odd
   else
     Nlambda = line%Nlambda+1 !result of MOD is 0, even
   end if
  else
   if (MOD(line%Nlambda,2).ge.1.) then
     Nlambda = line%Nlambda/2 !odd
   else
     Nlambda = (line%Nlambda+1)/2 !result of MOD is 0, even
   end if
  end if
  !write(*,*) "Nlambda, line%Nlambda = ",Nlambda, line%Nlambda

  if (line%qwing.le.2.*line%qcore) then
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

    if (qB_char/2. .ge.line%qcore) then
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
   double precision :: wlam, Dopplerwidth
   type (AtomicLine) :: line
   integer :: la
   !Return wavelength interval.
   !In case of atomic bound-bound
   ! return interval
   !in Doppler units
   ! Thermal velocity  not factorized here !!
   ! result in m/s

   Dopplerwidth = CLIGHT/line%lambda0
   if (line%symmetric) Dopplerwidth = Dopplerwidth*2
   if (la == 1) then
    wlam = 0.5*(line%lambda(la+1)-line%lambda(la))
   else if (la == line%Nlambda) then
    wlam = 0.5*(line%lambda(la)-line%lambda(la-1))
   else
    wlam = 0.5*(line%lambda(la+1)-line%lambda(la-1))
   end if

  RETURN
  END FUNCTION IntegrationWeightLine

  FUNCTION IntegrationWeightCont(cont, la) result(wlam)
  ! result in nm
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

  SUBROUTINE make_wavelength_grid(Natom, Atoms, inoutgrid, wl_ref)
  ! construct un sort a wavelength grid for NLTE line transfer
  ! after the grid has been constructed, saves Nblue for continuum
  ! and bound-bound transitions: The bluest absolute index of the extent of the
  ! transitions.
  ! Note that because duplicate wavelengths are removed, a new %Nlambda for
  ! each line and continuum are recomputed.
   type (AtomType), intent(inout), dimension(Natom) :: Atoms
   integer, intent(in) :: Natom
   double precision, intent(in) :: wl_ref
   type (AtomicLine), allocatable, dimension(:) :: alllines
   type (AtomicContinuum), allocatable, dimension(:) :: allcont
   double precision, allocatable, dimension(:), intent(inout) :: inoutgrid
   integer :: kr, kc, n, Nspect, Nwaves, Nlinetot, Nctot, idx
   integer :: la, nn, compar, Nred, Nblue, Nlambda_original!read from model
   double precision, allocatable, dimension(:) :: tempgrid
   double precision :: l0, l1 !ref wavelength of each transitions

   nn = 0
   allocate(alllines(MAX_ACTIVE_TRANSITIONS))
   allocate(allcont(MAX_ACTIVE_TRANSITIONS))

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
    if (Nlinetot.gt.MAX_ACTIVE_TRANSITIONS) then
     write(*,*) "too many Active transitions"
     stop
    end if
    alllines(Nlinetot) = Atoms(n)%lines(kr)
   end do
   do kc=1,Atoms(n)%Ncont
    Nspect = Nspect + Atoms(n)%continua(kc)%Nlambda
    Nctot = Nctot + 1
    if (Nctot.gt.MAX_ACTIVE_TRANSITIONS) then
     write(*,*) "too many Active transitions"
     stop
    end if
    allcont(Nctot) = Atoms(n)%continua(kc)
   end do
  end do ! end loop over atoms

  ! add ref wavelength if any and allocate temp array
  if (wl_ref.gt.0.) then
   Nspect = Nspect + 1
   nn = 1
   allocate(tempgrid(Nspect))
   tempgrid(1)=wl_ref
  else
   allocate(tempgrid(Nspect))
  end if

  write(*,*) "Adding ", Nspect," lines and continua, for a total of", Nlinetot+Nctot,&
             " transitions."

  ! add wavelength from mcfost inoutgrid if any
  ! convert it to nm
  if (allocated(inoutgrid)) then
   do la=1, size(inoutgrid)
    nn = nn + 1
    tempgrid(nn) = inoutgrid(la)
   end do
   deallocate(inoutgrid)
  end if

  ! start adding continua and lines wavelength grid
  do kc=1,Nctot
   do la=1,allcont(kc)%Nlambda
    nn = nn + 1
    tempgrid(nn) = allcont(kc)%lambda(la)
   end do
  end do
  do kr=1,Nlinetot
   do la=1,alllines(kr)%Nlambda
    nn = nn + 1
    tempgrid(nn) = alllines(kr)%lambda(la)
   end do
  end do

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
  Nwaves = Nwaves-1
  write(*,*) Nwaves, " unique wavelengths, ", Nspect-Nwaves,&
   " eliminated lines"

  allocate(inoutgrid(Nwaves))
  do la=1,Nwaves
   inoutgrid(la) = tempgrid(la)
   ! write(*,*) "lam(la) =", inoutgrid(la)
  end do

  !should not dot that but error somewhere if many atoms
  !CALL sort(inoutgrid, Nwaves)

  deallocate(tempgrid)
  deallocate(alllines)
  deallocate(allcont)

  ! Now replace the line%Nlambda and continuum%Nlambda by the new values.
  ! we do that even for PASSIVE atoms
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
    end do
!     write(*,*) " ------------------------------------------------------------------ "
  end do !over atoms

!deallocate line%lambda ?continuum%lambda because it will not correspond to the new
! Nblue and Nlambda

  RETURN
  END SUBROUTINE make_wavelength_grid

  END MODULE getlambda
