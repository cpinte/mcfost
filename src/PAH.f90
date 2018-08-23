module PAH

  use parametres
  use grains
  use constantes
  use dust_prop
  use utils, only : interp

  implicit none

contains

function specific_heat(T, taille_grain)
  ! C. Pinte
  ! 30/01/07

  integer, intent(in) :: taille_grain ! indice taille de grain
  real(kind=dp), dimension(:), intent(in) :: T
  real(kind=dp), dimension(size(T)) :: specific_heat
  integer :: pop

  pop =  grain(taille_grain)%pop

  if (dust_pop(pop)%is_Misselt_opacity_file) then
     specific_heat = Misselt_specific_heat(T,taille_grain)
  else if (dust_pop(pop)%is_DustEM_opacity_file) then
     specific_heat = DustEM_specific_heat(T,taille_grain)
  else if (grain(taille_grain)%is_PAH) then
     specific_heat = PAH_specific_heat(T, taille_grain)
  else
     specific_heat = astrosil_specific_heat(T, taille_grain)
  endif

  return

end function specific_heat

!**********************************************************************

function astrosil_specific_heat(T, taille_grain)
  ! C. Pinte
  ! 26/01/07

  integer, intent(in) :: taille_grain ! indice taille de grain
  real(kind=dp), dimension(:), intent(in) :: T
  real(kind=dp), dimension(size(T)) :: astrosil_specific_heat

  real :: Na

  Na = get_astrosil_Na(taille_grain)

  ! specific heat [J/K]
  astrosil_specific_heat = (Na-2.)* kb * (2.*sh_helper(T/500., 2) + sh_helper(T/1500., 3))

end function astrosil_specific_heat

!**********************************************************************

function PAH_specific_heat(T,taille_grain)
  ! return the specific heat capacity for PAHs, in [erg/K]
  ! input is Temperature [K] and grain radius [cm]
  ! optional output are mode energies hbarw [erg] and modes/energy g
  ! C. Pinte
  ! 30/01/07

  integer, intent(in) :: taille_grain
  real(kind=dp), dimension(:), intent(in) :: T
  real(kind=dp), dimension(size(T)) :: PAH_specific_heat

  real(kind=dp), dimension(:), allocatable :: hbarw, g

  integer :: NC, NH, n, Nm, i, nT
  real(kind=dp) :: ThetaD, beta, a

  integer :: N_CCop, N_CCip, N_H, N_tot


  real(kind=dp), dimension(:), allocatable :: deltaj, g_CCop, g_CCip,  hbarw_CCop, hbarw_CCip, x
  integer, dimension(:), allocatable :: s

  real, dimension(3) :: wavenum, hbarw_CH, g_CH

  logical, save :: l_first_time = .true.

  a = r_grain(taille_grain)

  nT=size(T)

  NC = nint((a*1.e3)**3*468.)   ! number of Carbon atoms (Draine & Li 2001 Eq 8+)
  NH = get_NH(NC)               ! number of Hydrogen atoms

  if (l_first_time) then
     write(*,*) "PAH, NC=", NC
     l_first_time = .false.
  endif

  if (NC > limite_stack) then
     call error("PAH grain size is too large, reaching stack limit", &
          msg2="Cannot compute the mode spectra")
  endif

  !compute mode spectra
  n = 2                        ! PAH C-C

  ! out-of-plane C-C
  Nm = NC-2   ; N_CCop = Nm    ! number of modes
  ThetaD = 863.                ! Debye temperature [K]

  allocate(deltaj(Nm),hbarw_CCop(Nm))
  deltaj=0.5 ; deltaj(2:3) =1.0 ! PAH C-C (DL01 eqs. 5,6)
  beta = get_beta(NC, Nm)
  ! do i=1, nT
  !    if (T(i) > ThetaD) then
  !       write(*,*) 'temp greater than Debye!'
  !    endif
  ! enddo

  hbarw_CCop = mode_spectrum(ThetaD, Nm, beta, deltaj) ! mode spectrum

  allocate(g_CCop(Nm)) ; g_CCop = 1. ! number of oscillators per energy
  deallocate(deltaj)

  ! in-plane C-C
  Nm = 2*NC-2  ; N_CCip = Nm    ! number of modes
  ThetaD = 2504.                ! Debye temperature [K]

  allocate(deltaj(Nm),hbarw_CCip(Nm))
  deltaj = 0.5  ; deltaj(2:3) = 1. ! PAH C-C (DL01 eqs. 5,6)
  beta = get_beta(NC, Nm)
  ! do i=1,NT
  !    if (T(i) > ThetaD) then
  !       write(*,*) 'temp greater than Debye!'
  !    endif
  ! enddo
  hbarw_CCip = mode_spectrum(ThetaD, Nm, beta, deltaj) ! mode spectrum
  allocate(g_CCip(Nm)) ;  g_CCip = 1.   ! number of oscillators per energy
  deallocate(deltaj)

  ! C-H modes (out-of-plane bending, in-plane bending, stretching)
  N_H=3
  ! mode wavenumbers
  wavenum(1) = 688.
  wavenum(2) = 1161.
  wavenum(3) = 3030.
  wavenum = wavenum * 1e2  ! [cm^-1] -> [m^-1]
  hbarw_CH = hp*c_light*wavenum  ! mode spectrum [erg] -> [J]
  g_CH = NH                ! number of oscillators per energy


  ! combine mode spectra
  N_tot = N_CCop +  N_CCip + N_H

  if (N_tot > limite_stack) then
     call error("PAH grain size is too large, reaching stack limit", &
          msg2="Cannot compute the mode spectra")
  endif

  allocate(hbarw(N_tot),g(N_tot),s(N_tot), x(N_tot))
  do i=1,N_CCop
     hbarw(i)=hbarw_CCop(i)
     g(i) = g_CCop(i)
  enddo
  do i=1,N_CCip
     hbarw(N_CCop+i)=hbarw_CCip(i)
     g(N_CCop+i) = g_CCip(i)
  enddo
  do i=1,N_H
     hbarw(N_CCop+N_CCip+i)=hbarw_CH(i)
     g(N_CCop+N_CCip+i) = g_CH(i)
  enddo

  !hbarw = [hbarw_CCop, hbarw_CCip, hbarw_CH]
  !g = [g_CCop, g_CCip, g_CH]

  s = bubble_sort(hbarw) ! renvoie les indices tries (=sort en Yorick), slow
  hbarw = hbarw(s)
  g = g(s)

  ! Reproduction Fig 1 Draine 2001 avec test_PAH_specific_heat
  ! open(unit=1,file="hbarw.txt")
  ! do i=1, n_tot
  !    write(1,*) hbarw(i)/(hp*c_light*1e2) , g(i)
  ! enddo
  ! close(unit=1)

  ! compute heat capacity
  do i = 1, nT
     x = hbarw/(kb*T(i))
     PAH_specific_heat(i) = kb*sum(g*exp(-x)*(x/(1. -exp(-x)))**2)
  end do

  return

end function PAH_specific_heat

!******************************************************

subroutine test_PAH_specific_heat()


  implicit none

  integer :: nc
  real :: a
  real(kind=dp), dimension(1) :: T, C

  nc=24
  a=1d-3*(nc/468.)**(1/3.)
  r_grain(1) = a
  T(1) = 500.

  C = PAH_specific_heat(T,1)

  call exit(0)

end subroutine test_PAH_specific_heat

!******************************************************

function Misselt_specific_heat(T,taille_grain)
  ! return the specific heat capacity in [erg/K]

  integer, intent(in) :: taille_grain ! indice taille de grain
  real(kind=dp), dimension(:), intent(in) :: T
  real(kind=dp), dimension(size(T)) :: Misselt_specific_heat

  integer :: k, pop

  pop =  grain(taille_grain)%pop

  ! todo : interpoler a l'avance dans read_file_specific_heat
  do k=1,size(T)
     Misselt_specific_heat(k) = interp(file_sh(:,pop), file_sh_T(:,pop), T(k)) * M_grain(taille_grain)  / 1.e7  * 0.5  ! todo : c'est des Joules l'unite non ?
     ! todo 2 : facteur 0.5 ??? donne un meilleur accord mais je ne capte pas pourquoi
  enddo

  return

end function Misselt_specific_heat

!**********************************************************************

function DustEM_specific_heat(T,taille_grain)
  ! return the specific heat capacity in [erg/K]

  integer, intent(in) :: taille_grain ! indice taille de grain
  real(kind=dp), dimension(:), intent(in) :: T
  real(kind=dp), dimension(size(T)) :: DustEM_specific_heat

  integer :: k, pop

  pop =  grain(taille_grain)%pop

  ! todo : interpoler a l'avance dans read_file_specific_heat
  do k=1,size(T)
     DustEM_specific_heat(k) = interp(file_sh(:,pop), file_sh_T(:,pop), T(k)) * M_grain(taille_grain)
  enddo

  return

end function DustEM_specific_heat

!**********************************************************************

real function get_astrosil_Na(taille_grain)
  ! C. Pinte
  ! 26/01/07

  implicit none

  integer, intent(in) :: taille_grain
  real :: a

  a=r_grain(taille_grain)

  get_astrosil_Na=4.*pi/3. * a**3 * 3.7e10

  !Na = 4.d*!dpi/3.*a^3 * 3.7d22

end function get_astrosil_Na

!******************************************************

function mode_spectrum(ThetaD,Nm,beta,deltaj)
  ! Draine & Lee eq 4 with n=2
  ! C. Pinte
  ! 30/01/07

  implicit none

  real(kind=dp), intent(in) :: ThetaD, beta
  integer, intent(in) :: Nm
  real(kind=dp), dimension(Nm) :: mode_spectrum
  real(kind=dp), dimension(Nm), intent(in) :: deltaj

  integer :: j
  real :: fact

  fact = (1.-beta)/Nm

  do j=1,Nm
     mode_spectrum(j) = sqrt(fact*(j-deltaj(j))+beta)
  enddo
  mode_spectrum= mode_spectrum * kb * ThetaD

  return

end function mode_spectrum

!******************************************************

function get_beta(NC, Nm)
  ! helper for specific heat calc (DL01 eq. 7)
  ! C. Pinte
  ! 30/01/07
  implicit none

  integer, intent(in) :: NC, Nm
  real :: get_beta

  if (NC <= 54) then
     get_beta = 0 ! PAH C-C (DL01 eq. 7)
  else if (NC <= 102) then
     get_beta = (NC-52.)/52./(2.*Nm-1)
  else
     get_beta = ((NC-54.)/52.*(102./NC)**(2./3.)-1.)/(2.*Nm-1)
  end if

  return

end function get_beta

!**********************************************************************

function sh_helper(x, n)
  ! helper for the specific heats
  ! returns f'n(x)  (DL01, eq. 10)  nn = 100
  ! C. Pinte
  ! 26/01/07

  implicit none

  real(kind=dp), dimension(:), intent(in) :: x
  integer, intent(in) :: n
  real, dimension(size(x)) :: sh_helper

  integer, parameter :: nn=100
  integer :: i, nx, j
  real :: dy
  real, dimension(nn) :: y
  real(kind=dp) :: eyx, yx

  do j=1, nn
     y(j)=(real(j)-0.5)/real(nn);
  enddo
  dy = 1.0/real(nn);

  nx=size(x)

  sh_helper=0.0
  do i=1, nx
     do j=1, nn
        yx = y(j)/x(i)
        if (yx < 350.) then ! limite double precision representable
           eyx = exp(yx)
           sh_helper(i)= sh_helper(i) + (y(j)**(n+1)*eyx/(eyx-1.)**2)
        endif
     enddo
  enddo

  sh_helper=sh_helper*dy*real(n)/x**2;

end function sh_helper

!************************************************************

end module PAH
