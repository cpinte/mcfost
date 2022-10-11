module utils

  use mcfost_env
  use naleat, only : stream
  use constantes
  use sha
  use messages

  implicit none

  real, parameter :: VACUUM_TO_AIR_LIMIT=200.0000
  real, parameter :: AIR_TO_VACUUM_LIMIT=199.9352

  public :: interp
  interface interp
     module procedure  interp_sp
     module procedure  interp_dp
  end interface

contains

function span(xmin,xmax,n)

  implicit none

  real, intent(in) :: xmin, xmax
  integer, intent(in) :: n
  real, dimension(n) :: span

  integer :: i
  real :: delta_x

  delta_x = (xmax-xmin)/real(n-1)
  span(1) = xmin
  do i=2,n
     span(i) = span(i-1) + delta_x
  enddo

  return

end function span

function span_dp(xmin,xmax,n, dk)

  implicit none
  integer, intent(in) ::dk
  real(kind=dp), intent(in) :: xmin, xmax
  integer, intent(in) :: n
  real(kind=dp), dimension(n) :: span_dp
  real(kind=dp) :: x1, x2
  integer :: istart, iend, i0

  integer :: i
  real(kind=dp) :: delta_x

  delta_x = (xmax-xmin)/real(n-1,kind=dp)
  
  if (dk < 0) then
  	x1 = xmax
  	x2 = xmin
  	istart = N-1
  	iend = 1
  	i0 = N
  else
  	x1 = xmin
  	x2 = xmax
  	istart = 2
  	iend = n
  	i0 = 1
  endif
  
  span_dp(i0) = x1
  do i=istart,iend,dk
     span_dp(i) = span_dp(i-dk) + dk * delta_x
  enddo

  return

end function span_dp

!************************************************************

function spanl(xmin,xmax,n)

  implicit none

  real, intent(in) :: xmin, xmax
  integer, intent(in) :: n
  real, dimension(n) :: spanl

  spanl = exp(span(log(abs(xmin)), log(abs(xmax)), n))

end function spanl

function spanl_dp(xmin,xmax,n,dk)

  implicit none
  integer, intent(in) :: dk
  real(kind=dp), intent(in) :: xmin, xmax
  integer, intent(in) :: n
  real(kind=dp), dimension(n) :: spanl_dp

  spanl_dp = exp(span_dp(log(abs(xmin)), log(abs(xmax)), n, dk))

end function spanl_dp

!************************************************************

real(kind=sp) function interp_sp(y, x, xp)
! interpole lineaire entre les points
! fait une extrapolation lineaire entre les 2 premiers ou 2 derniers points
! en cas de point cp en dehors de x
! Modif : pas d'extrapolation : renvoie valeur a xmin ou xmax
! C. Pinte
! 26/01/07

  implicit none

  real(kind=sp), dimension(:), intent(in) :: x, y
  real(kind=sp) :: xp, frac

  integer :: n, j

  n=size(x) !; ny=size(y)

  !if (n /= ny) then
  !   write(*,*) "Error in interp : y and x must have same dim"
  !endif

  ! PAS D'EXTRAPOLATION
  if (xp < minval(x)) then
     if (x(n) > x(1)) then ! xcroissant
        interp_sp = y(1)
     else
        interp_sp = y(n)
     endif
     return
  endif
  if (xp > maxval(x)) then
     if (x(n) > x(1)) then ! xcroissant
        interp_sp = y(n)
     else
        interp_sp = y(1)
     endif
     return
  endif

  ! Suivant la croissance ou decroissance de x
  if (x(n) > x(1)) then ! xcroissant
     search : do j=2,n-1
        if (x(j) > xp) exit search
     enddo search
     frac = (xp-x(j-1))/(x(j)-x(j-1))
     interp_sp = y(j-1) * (1.-frac)   + y(j) * frac
  else ! x decroissant
     search2 : do j=2,n-1
        if (x(j) < xp) exit search2
     enddo search2
     frac = (xp-x(j))/(x(j-1)-x(j))
     interp_sp = y(j) * (1.-frac)   + y(j-1) * frac
  endif

  return

end function interp_sp

!----------------------------

real(kind=dp) function interp_dp(y, x, xp)
! interpole lineaire entre les points
! fait une extrapolation lineaire entre les 2 premiers ou 2 derniers points
! en cas de point cp en dehors de x
! Modif : pas d'extrapolation : renvoie valeur a xmin ou xmax
! C. Pinte
! 26/01/07

  implicit none

  real(kind=dp), dimension(:), intent(in) :: x, y
  real(kind=dp) :: xp, frac

  integer :: n, j

  n=size(x) !; ny=size(y)

  !if (n /= ny) then
  !   write(*,*) "Error in interp : y and x must have same dim"
  !endif


  ! PAS D'EXTRAPOLATION
  if (xp < minval(x)) then
     if (x(n) > x(1)) then ! xcroissant
        interp_dp = y(1)
     else
        interp_dp = y(n)
     endif
     return
  endif
  if (xp > maxval(x)) then
     if (x(n) > x(1)) then ! xcroissant
        interp_dp = y(n)
     else
        interp_dp = y(1)
     endif
     return
  endif

  ! Suivant la croissance ou decroissance de x
  if (x(n) > x(1)) then ! xcroissant
     search : do j=2,n-1
        if (x(j) > xp) exit search
     enddo search
     frac = (xp-x(j-1))/(x(j)-x(j-1))
     interp_dp = y(j-1) * (1.-frac)   + y(j) * frac
  else ! x decroissant
     search2 : do j=2,n-1
        if (x(j) < xp) exit search2
     enddo search2
     frac = (xp-x(j))/(x(j-1)-x(j))
     interp_dp = y(j) * (1.-frac)   + y(j-1) * frac
  endif

  return

end function interp_dp

!**********************************************************************

subroutine GaussSlv(a, b, n)
  ! Resolution d'un systeme d'equation par methode de Gauss
  ! Non optimise : meme pas de recherche du pivot max !!!
  ! C. Pinte
  ! 22/09/07

  implicit none

  integer, intent(in) :: n
  real(kind=dp), dimension(n,n), intent(inout) :: a
  real(kind=dp), dimension(n), intent(inout) :: b

  real(kind=dp) :: factor
  integer :: i, j, k, l

  if (n == 1) then
     b(n) = b(n) / a(n,n)
     return
  end if

  ! Triangularisation de la matrice
  do i=1,n-1 ! Boucle sur colonnes
     do k=i+1,n ! Boucle sur lignes
        factor = a(k,i)/a(i,i)
        ! Operation sur la ligne
        do j=i+1,n
           a(k,j) = a(k,j) - a(i,j) * factor
        end do
        b(k) = b(k) - b(i) * factor
     end do
  end do
  b(n) = b(n) / a(n,n)

  ! Resolution de la matrice triangulaire
  do i=1,n-1
     k = n-i
     l = k+1
     do j = l, n
        b(k) = b(k) - b(j)*a(k,j)
     end do
     b(k) = b(k) / a(k,k)
  end do

  return

end subroutine GaussSlv

!***********************************************************

subroutine rotation(xinit,yinit,zinit,u1,v1,w1,xfin,yfin,zfin)
  ! Effectue une rotation du vecteur (xinit,yinit,zinit)
  ! Le resultat (xfin,yfin,zfin) est dans
  ! le systeme de coordonnees o� le vecteur (u1,v1,w1)=(1,0,0).
  ! ie applique la meme rotation que celle qui transforme
  ! (1,0,0) en (u1,v1,w1) sur (xinit,yinit,zinit)
  ! C. Pinte : 1/03/06
  ! Nouvelle version d'une routine de F. M�nard

  implicit none

  real(kind=dp), intent(in) ::  xinit, yinit, zinit, u1, v1, w1
  real(kind=dp), intent(out) :: xfin, yfin, zfin

  real(kind=dp) :: cost, sint, sing, prod, theta

  if (w1 > 0.999999999_dp) then
     cost = 1.0_dp
     sint = 0.0_dp
     sing = 0.0_dp
  else
     if (abs(u1) < tiny_real) then
        cost=0.0_dp
        sint=1.0_dp
        sing=sqrt(1.0 - w1*w1)
     else
   !     x=v1/u1
   !     !c'est pas un atan mais un atan2 d'on le sign(u1)
   !     cost=sign(1.0/sqrt(1.0+x*x),u1)        !cos(atan(x)) = 1/(1+sqrt(x*x))
   !     sint=x*cost                            !sin(atan(x)) = x/(1+sqrt(x*x))
        ! Equivalent  �  (~ meme tps cpu):
         theta=atan2(v1,u1)
         cost=cos(theta)
         sint=sin(theta)
         sing = sqrt(1.0_dp - w1*w1)
      endif
  endif

  prod=cost*xinit + sint*yinit

  xfin = sing * prod + w1*zinit
  yfin = cost*yinit - sint*xinit
  zfin = sing*zinit - w1 * prod

  return

end subroutine rotation

!***********************************************************

function calc_mu0(mu,a)
  ! Resout l'equation de degre 3:   mu0**3 +(a-1)*mu0 - a*mu == 0
  ! Forme speciale de l'equation:  x**3 +a1*x**2 +a2*x +a3,
  ! Utilisation des formules de Tartaglia-Cardan
  ! C. Pinte
  ! 01/09/08

  real(kind=dp), intent(in) :: mu, a
  real(kind=dp) :: calc_mu0

  real(kind=dp) :: a1, a2, a3, q, r, q3, qh, denom, theta, mu01, mu02, mu03, factor, factor3


  a1 = 0.0_dp
  a2 = a - 1.0_dp
  a3 = -a * mu

  q = (a1**2 - 3.0_dp*a2)/9.0_dp
  r = (2.0_dp*a1**3 - 9.0_dp*a1*a2 + 27.0_dp*a3)/54.0_dp
  q3 = q**3

  if ((q3 - r**2) >= 0.0_dp) then
     qh = sqrt(q)
     denom = sqrt(q3)
     theta = acos(r/denom)
     mu01 = -2.0_dp*qh*cos(theta/3.0_dp) - a1/3.0_dp
     mu02 = -2._dp*qh*cos((theta + deux_pi)/3.0_dp) - a1/3.0_dp
     mu03 = -2._dp*qh*cos((theta + quatre_pi)/3.0_dp) - a1/3.0_dp
     if ((mu01 > 0.0_dp).and.(mu02 > 0.0_dp)) then
        call error("calc_mu0: 2 positives roots, mu01, mu02")
     else if ((mu01 > 0.0_dp).and.(mu03 > 0.0_dp)) then
        call error("calc_mu0: 2 positives roots, mu01, mu03")
     elseif ((mu02 > 0.0_dp).and.(mu03 > 0.0_dp)) then
        call error("calc_mu0: 2 positives roots, mu02, mu03")
     endif

     if (mu01 > 0.0_dp) calc_mu0 = mu01
     if (mu02 > 0.0_dp) calc_mu0 = mu02
     if (mu03 > 0.0_dp) calc_mu0 = mu03

  else
     factor = sqrt(r**2 - q3) + abs(r)
     factor3 = factor**un_tiers
     calc_mu0 = -1.0_dp*sign(1.0_dp,r)*(factor3 + q/factor3) - a1/3.0_dp
  endif

  return

end function calc_mu0

!***********************************************************

function Bnu(nu,T)
! Loi de Planck
! Bnu en SI : W.m-2.Hz-1.sr-1
! nu en Hz
! C. Pinte
! 17/10/7

  implicit none

  real(kind=dp) ,intent(in) :: nu
  real, intent(in) :: T
  real(kind=dp) :: Bnu

  real(kind=dp) :: hnu_kT

  hnu_kT = (hp * nu) / (kb * T)

  if (hnu_kT > 100._dp) then
     Bnu = 0.0_dp
  else
     Bnu = 2.0_dp*hp/c_light**2 * nu**3 / (exp(hnu_kT)-1.0_dp)
  endif

  return

end function Bnu

!***********************************************************

!-> note. Elemental function are too slow!
Function Bpnu (N,lambda,T)
! -----------------------------------------------------
! Return an array of planck functions at all wavelengths
! for the cell temperature.
! Bnu in W/m2/Hz/sr
! lambda in nm
! -----------------------------------------------------
integer, intent(in) :: N
real(kind=dp), intent(in) :: T, lambda(N)
real(kind=dp) :: hnu_kT, twohnu3_c2, Bpnu(N)
integer la

do la=1, N
   twohnu3_c2 = 2.*HC / (NM_TO_M * lambda(la))**3
   hnu_kT = hc_k / lambda(la) / T

   if (hnu_kT > 100.0) then
      Bpnu(la) = 0.0
   else
      Bpnu(la) = twohnu3_c2 / (exp(hnu_kT)-1.0)
   end if 
enddo

return
end function bpnu

!***********************************************************

function Blambda(wl,T)
! Loi de Planck
! Blambda en SI : W.m-2.m-1.sr-1
! wl en m
! C. Pinte
! 23/05/09

  implicit none

  real ,intent(in) :: wl
  real, intent(in) :: T
  real :: Blambda

  real(kind=dp) :: hnu_kT

  hnu_kT = (hp * c_light)/ (kb * T * wl)

  if (hnu_kT > 100.) then
     Blambda = 0.0
  else
     Blambda = 2.0*hp*c_light**2 / wl**5 / (exp(hnu_kT)-1.0)
  endif

  return

end function Blambda

!******************************************************

function Blambda_dp(wl,T)
! Loi de Planck
! Blambda en SI : W.m-2.s-1.sr-1
! wl en m
! C. Pinte
! 23/05/09

  implicit none

  real(kind=dp) ,intent(in) :: wl
  real, intent(in) :: T
  real(kind=dp) :: Blambda_dp

  real(kind=dp) :: hnu_kT

  hnu_kT = (hp * c_light)/ (kb * T * wl)

  if (hnu_kT > 700.) then
     Blambda_dp = 0.0_dp
  else
     Blambda_dp = 2.0_dp*hp*c_light**2 / wl**5 / (exp(hnu_kT)-1.0_dp)
  endif

  return

end function Blambda_dp

!******************************************************

function get_NH(NC)
  ! returns the number of Hydrogen atoms given the number of Carbon
  ! atoms (DL01 eq. 8)
  ! C. Pinte
  ! 30/01/07
  implicit none

  integer, intent(in) :: NC
  integer :: get_NH

  if (NC <= 25) then
     get_NH = floor(0.5*NC+0.5)
  else if (NC <= 100) then
     get_NH = floor(2.5*sqrt(real(NC))+0.5)
  else
     get_NH = floor(0.25*NC+0.5)
  endif

  return

end function get_NH

!***********************************************************

subroutine mcfost_setup()

  character(len=512) :: cmd, s
  integer ::  syst_status, mcfost_git

  call get_utils()
  ! Write date of the last time an update was search for
  cmd = "rm -rf "//trim(mcfost_utils)//"/.last_update"//" ; date +%s > "&
       //trim(mcfost_utils)//"/.last_update"
  call appel_syst(cmd, syst_status)

  ! We do not need to download the parameter file if we use the sources
  call get_environment_variable('MCFOST_GIT',s)
  if (s/="") then
     read(s,*) mcfost_git
     if (mcfost_git /= 1) call mcfost_get_ref_para()
  endif

  write(*,*) "MCFOST set-up was sucessful"

  write(*,*) "You can find the full documentation at:"
  write(*,*) trim(doc_webpage)

  return

end subroutine mcfost_setup

!***********************************************************

function mcfost_update(lforce_update, lmanual, n_days)

  use os

  logical, intent(in) :: lforce_update, lmanual
  integer, intent(in), optional :: n_days
  logical :: lupdate, mcfost_update

  character(len=512) :: cmd, url, url_sha1, last_version, current_binary, s
  character(len=40) :: mcfost_sha1, mcfost_update_sha1
  integer ::  syst_status, ios


  ! Last version
  write(*,*) "Checking last version ..."
  syst_status = 0
  !cmd = "wget "//trim(webpage)//"/version.txt -q -T 5 -t 3"
  cmd = "curl "//trim(webpage)//"version.txt -O -s"
  call appel_syst(cmd, syst_status)
  if (syst_status/=0) then
     write(*,*) "ERROR: Cannot connect to MCFOST server."
     call exit_update(lmanual, n_days, lupdate)
     return
  endif

  ios=0
  open(unit=1, file="version.txt", status='old',iostat=ios)
  read(1,*,iostat=ios) last_version
  close(unit=1,status="delete",iostat=ios)
  if ( (ios/=0) .or. (.not.is_digit(last_version(1:1)))) then
     write(*,*) "ERROR: Cannot get MCFOST last version number."
     write(*,*) "Cannot read version file."
     call exit_update(lmanual, n_days, lupdate)
     return
  endif

  ! Do we need to update ?
  if (last_version == mcfost_release) then
     write(*,*) "MCFOST is up-to-date"
     lupdate = .false.
  else ! Updating
     write(*,*) "New MCFOST version available: ", trim(last_version)
     lupdate = .true.
  endif
  if (lforce_update) then
     write(*,*) "Forcing update"
     lupdate = .true.
  endif
  write(*,*) " "

  if (lupdate) then
     ! get the correct url corresponding to the system
     if (operating_system=="Linux ") then
        url = trim(webpage)//"linux/mcfost_bin.tgz"
        url_sha1 = trim(webpage)//"linux/mcfost.sha1"
     else if (operating_system=="Darwin") then
        url = trim(webpage)//"macos/mcfost_bin.tgz"
        url_sha1 = trim(webpage)//"macos/mcfost.sha1"
     else
        write(*,*) "Unknown operating system : error 2"
        write(*,*) "Cannot download new binary"
        call exit_update(lmanual, n_days, lupdate)
        return
     endif

     write(*,*) "Your system is ", operating_system

     ! Download
     write(*,'(a32, $)') "Downloading the new version ..."
     cmd = "curl "//trim(url)//" -o mcfost_bin.tgz -s"    ; call appel_syst(cmd, syst_status)
     cmd = "curl "//trim(url_sha1)//" -o mcfost.sha1 -s" ; call appel_syst(cmd, syst_status)
     if (syst_status==0) then
        write(*,*) "Done"
     else
        cmd = "rm -rf mcfost_bin.tgz"
        call appel_syst(cmd, syst_status)
        write(*,*) "ERROR during download. MCFOST has not been updated."
        call exit_update(lmanual, n_days, lupdate)
        return
     endif

     ! check sha
     write(*,'(a20, $)') "Checking binary ..."
     if (operating_system=="Linux ") then
        cmd = "sha1sum  mcfost_bin.tgz > mcfost_update.sha1"
     else if (operating_system=="Darwin") then
        cmd = "openssl sha1 mcfost_bin.tgz | awk '{print $2}' > mcfost_update.sha1"
     endif
     call appel_syst(cmd, syst_status)

     open(unit=1, file="mcfost.sha1", status='old',iostat=ios)
     read(1,*,iostat=ios) mcfost_sha1
     close(unit=1,status="delete",iostat=ios)

     open(unit=1, file="mcfost_update.sha1", status='old',iostat=ios)
     read(1,*,iostat=ios) mcfost_update_sha1
     close(unit=1,status="delete",iostat=ios)

     if ( (ios/=0) .or. (mcfost_sha1/=mcfost_update_sha1)) then
        !cmd = "rm -rf mcfost_update" ; call appel_syst(cmd, syst_status)
        write(*,*) " "
        write(*,*) "ERROR: binary sha1 is incorrect. MCFOST has not been updated."
        write(*,*) "The downloaded file has sha1: "
        write(*,*) mcfost_update_sha1
        write(*,*) "It should be:"
        write(*,*) mcfost_sha1
        call exit_update(lmanual, n_days, lupdate)
        return
     else
        write(*,*) "Done"
        write(*,'(a25, $)') "Decompressing binary ..."
        cmd = "tar xzf mcfost_bin.tgz ; rm -rf mcfost_bin.tgz"
        call appel_syst(cmd, syst_status)
        write(*,*) "Done"
     endif

     ! check where is the current binary
     call get_command_argument(0,current_binary)
     if (current_binary(1:1)/=".") then

        write(*,'(a28, $)') "Locating current binary ..."
        cmd = "rm -rf which_mcfost_binary.txt && which "//trim(current_binary)// &
             " | awk '{print "//' "\"" $NF "\""'//"}' > which_mcfost_binary.txt"
        call appel_syst(cmd, syst_status)

        ios=0
        open(unit=1, file="which_mcfost_binary.txt", status='old',iostat=ios)
        read(1,*,iostat=ios) current_binary
        close(unit=1,iostat=ios)

        if ( (ios/=0) .or. (.not.is_digit(last_version(1:1)))) then
           write(*,*) ""
           write(*,*) "ERROR: Cannot locate current MCFOST binary,"
           write(*,*) "the new binary will downloaded in the current directory."
        else
           write(*,*) "Done"
        endif
     endif

     ! make binary executable
     write(*,'(a20, $)') "Updating binary ..."
     cmd = "chmod a+x mcfost_update ; mv mcfost_update "//trim(current_binary)
     call appel_syst(cmd, syst_status)
     if (syst_status /= 0) then
        write(*,*) "ERROR : the update failed for some unknown reason"
        write(*,*) "You may want to have a look at the file named mcfost_update"
        call exit_update(lmanual, n_days, lupdate)
        return
     endif
     write(*,*) "Done"
     write(*,*) "MCFOST has been updated"
  endif ! lupdate

  ! Write date of the last time an update was search for
  cmd = "rm -rf "//trim(mcfost_utils)//"/.last_update"//" ; date +%s > "//trim(mcfost_utils)//"/.last_update"
  call appel_syst(cmd, syst_status)

  mcfost_update = lupdate

  return

end function mcfost_update

!***********************************************************

subroutine exit_update(lmanual, n_days, lupdate)
  ! Exit update properly :
  !  - exit mcfost if it is a manual update
  !  - simply return if it is an auto-update

  logical, intent(in) :: lmanual
  integer, intent(in) :: n_days
  logical, intent(out) :: lupdate

  character(len=512) :: cmd, s
  integer :: syst_status

  if (lmanual) then
     write(*,*) "Exiting."
     call exit(1)
  else ! if it is an auto-update, we don't exit
     lupdate = .false.
     ! We try again tomorrow : Write date of the last time an update was search for - (mcfost_auto_update -1) days
     write(*,*) "WARNING: Skiping auto-update. MCFOST will try again tomorrow."
     write(s,*) (n_days-1) * 3600 * 24
     cmd = "rm -rf "//trim(mcfost_utils)//"/.last_update"//" ; expr `date +%s` - "&
          //trim(s)//" > "//trim(mcfost_utils)//"/.last_update"
     call appel_syst(cmd, syst_status)
  endif

  return

end subroutine exit_update

!***********************************************************

subroutine mcfost_history()

  character(len=512) :: cmd
  integer ::  syst_status

  ! Last version
  write(*,*) "Getting MCFOST history ..."
  cmd = "curl "//trim(webpage)//"history.txt"
  call appel_syst(cmd, syst_status)
  if (syst_status/=0) call error("Cannot get MCFOST history")
  write(*,*) " "

  return

end subroutine mcfost_history

!***********************************************************

subroutine mcfost_get_ref_para()

  character(len=512) :: cmd
  character(len=12) :: ref_file
  character(len=18) :: ref_file_multi
  character(len=15) :: ref_file_3D
  integer ::  syst_status

  if (mcfost_release(5:5) == ".") then
     ref_file = "ref"//mcfost_release(1:4)//".para"
     ref_file_multi = "ref"//mcfost_release(1:4)//"_multi.para"
     ref_file_3D = "ref"//mcfost_release(1:4)//"_3D.para"
  else if (mcfost_release(4:4) == ".") then
     ref_file = "ref"//mcfost_release(1:3)//".para"
     ref_file_multi = "ref"//mcfost_release(1:3)//"_multi.para"
     ref_file_3D = "ref"//mcfost_release(1:3)//"_3D.para"
  else
     write(*,*) "Cannot parse "//trim(mcfost_release)//" to find parameter files"
     return
  endif


  write(*,*) "Getting MCFOST reference files: "//ref_file//" & "//ref_file_multi//" & "//ref_file_3D
  cmd = "curl "//trim(webpage)//ref_file//" -O -s"
  call appel_syst(cmd, syst_status)
  cmd = "curl "//trim(webpage)//ref_file_multi//" -O -s"
  call appel_syst(cmd, syst_status)
  cmd = "curl "//trim(webpage)//ref_file_3D//" -O -s"
  call appel_syst(cmd, syst_status)
  if (syst_status/=0) call error("Cannot get MCFOST reference file")
  write(*,*) "Done"

  return

end subroutine mcfost_get_ref_para

!***********************************************************

subroutine mcfost_get_yorick()

  character(len=512) :: cmd
  character(len=128) :: yo_dir, yo_file1, yo_file2
  integer ::  syst_status

  yo_dir = "yorick/"
  yo_file1 = "mcfost_struct.i"
  yo_file2 = "mcfost_utils.i"

  write(*,*) "Getting MCFOST yorick files: ", trim(yo_file1), " ", trim(yo_file2)
  cmd = "curl "//trim(webpage)//trim(yo_dir)//trim(yo_file1)//" -O -s ; "//&
       "curl "//trim(webpage)//trim(yo_dir)//trim(yo_file2)//" -O -s"

  call appel_syst(cmd, syst_status)
  if (syst_status/=0) call error("Cannot get MCFOST yorick package")
  write(*,*) "Done"

  return

end subroutine mcfost_get_yorick

!***********************************************************

subroutine mcfost_v()

  character(len=512) :: cmd, last_version
  integer ::  syst_status, line_number, ios
  character(len=128) :: sline_number

  syst_status = 0

  write(*,*) "Binary compiled the ",__DATE__," at ",__TIME__
#if defined (__INTEL_COMPILER)
  write(*,fmt='(" with INTEL compiler version ",i4)')  __INTEL_COMPILER
#endif
#if defined (__GFORTRAN__)
  write(*,fmt='(" with GFORTRAN compiler version ",i2,".",i1,"."i1)') __GNUC__ , __GNUC_MINOR__,  __GNUC_PATCHLEVEL__
#endif
#if defined (__G95__)
  write(*,fmt='(" with G95 compiler version ",i1,".",i2)') __G95__, __G95_MINOR__
#endif
  write(*,*) " "

  ! Last version
  write(*,*) "Checking last version ..."
  cmd = "curl "//trim(webpage)//"version.txt -O -s"
  call appel_syst(cmd, syst_status)
  if (syst_status/=0) &
       call error("Cannot get MCFOST last version number (Error 1)","'"//trim(cmd)//"' did not run as expected.")

  open(unit=1, file="version.txt", status='old',iostat=ios)
  read(1,*,iostat=ios) last_version
  close(unit=1,status="delete",iostat=ios)

  if ((ios/=0) .or. (.not.is_digit(last_version(1:1)))) &
       call error("Cannot get MCFOST last version number (Error 2)","Cannot read new version file")

  ! Do we have the last version ?
  if (last_version == mcfost_release) then
     write(*,*) "MCFOST is up-to-date"
  else ! Print the history of new version
     write(*,*) "A new version of MCFOST is available: ", trim(last_version)

     ! Last version
     write(*,*)
     write(*,*) "Getting MCFOST history ..."
     cmd = "curl "//trim(webpage)//"history.txt -O -s"
     call appel_syst(cmd, syst_status)
     if (syst_status/=0) call error("Cannot get MCFOST history")

     ! Getting line number of current version in history
     cmd = "grep -n "//trim(mcfost_release)//" history.txt | awk -F : '{print $1}' > line_number.txt"
     call appel_syst(cmd, syst_status)
     open(unit=1, file="line_number.txt", status='old',iostat=ios)
     read(1,*,iostat=ios) line_number
     close(unit=1,status="delete",iostat=ios)

     ! Printing history
     if (ios==0) then
        if (line_number > 1) then
           write(sline_number,*)  line_number-1
           write(*,*) "Changes since your version:"
           cmd = "head -"//trim(adjustl(sline_number))//" history.txt ; rm -rf history.txt"
           call appel_syst(cmd, syst_status)
        else
           write(*,*) "No changes found since your version (?????)"
        endif
        write(*,*) " "
        write(*,*) "Full history is available with mcfost -h"
        write(*,*) "The new version can be downloaded with mcfost -u"
     else
        cmd = "rm -rf history.txt"
        call appel_syst(cmd, syst_status)
        write(*,*) "ERROR: I did not find any modifications in the history"
        write(*,*) "Exiting"
     endif ! ios
  endif ! last_version = mcfost_release

  return

end subroutine mcfost_v

!***********************************************************

subroutine update_utils(lforce_update)

  logical, intent(in) :: lforce_update
  real :: utils_current_version, last_version

  character(len=512) :: cmd, s_last_version
  integer ::  syst_status, ios
  logical :: lupdate

  ! Last version
  utils_current_version =  get_mcfost_utils_version()
  write(*,*) "Version ", utils_current_version
  write(*,*) "Checking last version of MCFOST UTILS ..."

  cmd = "curl "//trim(utils_webpage)//"Version -O -s"
  call appel_syst(cmd, syst_status)
  if (syst_status/=0) call error("Cannot get MCFOST UTILS last version number (Error 1)")
  open(unit=1, file="Version", status='old',iostat=ios)
  read(1,*,iostat=ios) s_last_version
  close(unit=1,status="delete",iostat=ios)

  if ( (ios/=0) .or. (.not.is_digit(s_last_version(1:1)))) then
     call error("Cannot get MCFOST UTILS last version number (Error 2)")
  endif
  read(s_last_version,*) last_version


  ! Do we need to update ?
  if (last_version == required_utils_version) then
     !Ok we have the correct version of mcfost
     if (last_version == utils_current_version) then
        write(*,*) "MCFOST UTILS is up-to-date"
        lupdate = .false.
     else ! we update
        lupdate = .true.
     endif
  else ! We need to update MCFOST first
     write(*,*) "New MCFOST UTILS version available: ", trim(s_last_version)
     write(*,*) "Please update mcfost first with mcfost -u"
     lupdate = .false.
  endif

  if (lforce_update) then
     write(*,*) "Forcing update"
     lupdate = .true.
  endif

  if (lupdate) call get_utils()

  return

end subroutine update_utils

!***********************************************************

subroutine get_utils()

  character(len=512) :: cmd
  integer :: syst_status

  write(*,*) "Downloading MCFOST UTILS (this may take a while) ..."
  write(*,*) "from "//trim(utils_webpage)//"mcfost_utils.tgz"

  cmd = "mkdir -p "//trim(mcfost_utils)//&
  " ; curl "//trim(utils_webpage)//"mcfost_utils.tgz | tar xzf - -C"//trim(mcfost_utils)
  call appel_syst(cmd, syst_status)
  if (syst_status == 0) then
     write(*,*) "Done"
  else
     write(*,*) "Error during download."
  endif

  return

end subroutine get_utils

!***********************************************************

real function get_mcfost_utils_version()

  integer :: ios

  ! Check utils directory
  get_mcfost_utils_version = 0.0
  open(unit=1, file=trim(mcfost_utils)//"/Version", status='old',iostat=ios)
  if (ios /= 0) then
     write(*,*) "I could not find the file: "//trim(mcfost_utils)//"/Version"
     get_mcfost_utils_version = 0.0
  endif
  read(1,*,iostat=ios) get_mcfost_utils_version
  close(unit=1,iostat=ios)
  if (ios /= 0) then
     write(*,*) "I could not read the version number in file: "//trim(mcfost_utils)//"/Version"
     get_mcfost_utils_version = 0.0
  endif

  if (abs(get_mcfost_utils_version) < 1e-6) call error("I could not find current MCFOST utils version")

  return

end function get_mcfost_utils_version

!***********************************************************

function indgen(n)

  integer, intent(in) :: n
  integer, dimension(n) :: indgen

  integer :: i

  do i=1,n
     indgen(i) = i
  enddo

  return

end function indgen

!************************************************************

function is_digit(ch) result(res)
  ! Returns .true. if ch is a digit (0,1,...,9) and .false. otherwise
  ! From stringmod.f90 by George Benthien
  ! http://www.gbenthien.net/strings/index.html

  character, intent(in) :: ch
  logical :: res

  select case(ch)
  case('0':'9')
     res=.true.
  case default
     res=.false.
  end select

  return

end function is_digit

!************************************************************

function is_diff(a,b)

  logical is_diff
  real(kind=dp), intent(in) :: a,b

  if (abs(a-b) > 0.5e-5 * abs(a+b)) then
     is_diff = .true.
  else
     is_diff = .false.
  endif

end function is_diff

!************************************************************

subroutine lgnd(lmax,x,p)
  ! Subroutine to generate Legendre polynomials P_L(X)
  ! for L = 0,1,...,LMAX with given X.

  ! This computer program is part of the book, "An Introduction to
  ! Computational Physics," written by Tao Pang and published and
  ! copyrighted by Cambridge University Press in 1997.
  ! Program 5.4 : originally named LGND

  ! 29/08/11 : C. Pinte : translated in fortran 90

  integer, intent(in) :: lmax
  real(kind=dp), intent(in) :: x
  real(kind=dp), dimension(0:lmax), intent(out) :: P

  integer :: l

  P(0) = 1._dp
  P(1) = x
  do l = 1, lmax-1
     P(l+1) = ((2.0*l+1)*x*P(l)-l*P(l-1))/(l+1)
  enddo

  return

end subroutine lgnd

!************************************************************

logical function is_dir(filename)

  character(len=*), intent(in) :: filename

#if defined (__INTEL_COMPILER)
  inquire(directory=trim(filename),exist=is_dir)
#else
  inquire(file=trim(filename),exist=is_dir)
#endif

  return

end function is_dir

!************************************************************

logical function is_file(filename)

  character(len=*), intent(in) :: filename

  inquire(file=trim(filename),exist=is_file)

  return

end function is_file

!************************************************************

function in_dir(filename,dir_list, status)

  character(len=*), intent(in) :: filename
  character(len=*), dimension(:), intent(in) :: dir_list
  character(len=128) :: in_dir

  integer, intent(out), optional  :: status

  logical :: is_file
  integer :: i

  in_dir = "None"
  if (present(status)) status = 1
  do i=1,size(dir_list)
     inquire(file=trim(dir_list(i))//"/"//trim(filename),exist=is_file)
     if (is_file) then
        in_dir = dir_list(i)
        if (present(status)) status = 0
        return
     endif
  enddo !i

  return

end function in_dir

!************************************************************

subroutine appel_syst(cmd, status)

#if defined (__INTEL_COMPILER)
  use IFPORT, only : system
#endif

  implicit none

  character(len=*), intent(in) :: cmd
  integer, intent(out), optional :: status

  status = system(cmd) ! linux
  !call system(cmd, status) ! aix

  return

end subroutine appel_syst

!************************************************************

subroutine Gauss_Legendre_quadrature(x1,x2,n, x,w)
  ! Computes the abscissas and weights of the
  ! Gauss-Legendre n-point quadrature formula.
  !
  ! Inputs :
  ! - x1 and x2 are the lower and upper limits of integration,
  ! - n is the number of points
  !
  ! Outputs :
  !  - arrays x and w of length n, containing abscissas and weights

  ! - 20/06/05 by RGA.
  ! - modified 31/01/13 by C. Pinte

  implicit none

  integer, intent(in) :: n  ! Order of the quadrature

  real(kind=dp), intent(in) :: x1, x2 ! Lower and upper integration limit
  real(kind=dp), dimension(n), intent(out) :: x, w ! abscissas and weights

  real(kind=dp), parameter :: eps = 3.0d-14
  real(kind=dp) :: dn, dj, p1, p2, p3, pp, xl, xm, z, z1
  integer :: i, j, m
  logical :: conv

  ! The roots are symmetric in the interval, so we only have to find half of them.
  m = (n+1)/2
  xm = 0.5_dp * (x2+x1)
  xl = 0.5_dp * (x2-x1)

  dn = real(n,kind=dp)

  do i=1,m ! Loop over the desired roots
     z = cos(pi*(real(i,kind=dp)-0.25_dp)/(dn+0.5_dp))
     ! Starting with the above approximation to the ith root, we
     ! enter the main loop of refinement by Newton's method
     conv = .false.
     do while (.not.conv)
        p1 = 1.0_dp ; p2 = 0.0_dp

        ! Loop up the recurrence relation to get the Legendre polynomial evaluated at z
        do j = 1, n
           dj=real(j,kind=dp)
           p3 = p2
           p2 = p1
           p1 = ((2.0_dp*dj-1.0_dp)*z*p2 - (dj-1.0_dp)*p3)/dj
        enddo ! j

        ! p1 is now the desired Legendre polynomial. We next compute pp, its
        ! derivative, by a standard relation involving also p2, the
        ! polynomial of one lower order.
        pp = dn * (z * p1 - p2) / (z * z - 1.0_dp)

        ! Newton's method.
        z1 = z
        z = z1 - p1 / pp
        conv = abs(z-z1) <= eps
     enddo ! conv

     ! Scale the root to the desired interval, and put in its symmetric counterpart.
     x(i) = xm - xl * z
     x(n+1-i) = xm + xl * z

     ! Compute the weight and its symmetric counterpart.
     w(i) = 2.0_dp * xl / ((1.0_dp - z * z) * pp * pp)
     w(n+1-i) = w(i)
  enddo ! i

  return

end subroutine Gauss_Legendre_quadrature

!************************************************************

logical function real_equality(x,y)

  real, intent(in) :: x, y

  real_equality = .false.
  if (abs(x) < tiny_real) then
     if (abs(y) < tiny_real) then
        real_equality = .true.
     endif
  else
     real_equality = abs(x-y) < 1e-5 * abs(x)
  endif

  return

end function real_equality

!************************************************************

function rotation_3d(axis,angle,v)
  ! Rotates a vector around an axis vector in 3D.
  !
  !  Parameters:
  !    - axis : the axis vector for the rotation. Must be normalized !
  !    - angle :the angle, in degrees, of the rotation.
  !    - v : the vector to be rotated.
  !  Output: the rotated vector.
  !
  ! C. Pinte 17/07/14


  real(kind=dp), dimension(3), intent(in) :: axis, v
  real(kind=dp), intent(in) :: angle
  real(kind=dp), dimension(3) :: rotation_3d

  real(kind=dp), dimension(3) :: vn, vn2, vp
  real(kind=dp) :: norm

  ! Compute the parallel component of the vector.
  vp(:) = dot_product(v, axis) * axis(:)

  ! Compute the normal component of the vector.
  vn(:) = v(:) - vp(:)

  norm = sqrt(sum(vn*vn))
  if (norm < tiny_dp) then
     rotation_3d(:) = vp(:)
     return
  endif
  vn(:) = vn(:)/norm

  ! Compute a second vector, lying in the plane, perpendicular
  ! to vn, and forming a right-handed system.
  vn2(:) = cross_product(axis, vn)

  ! Rotate the normal component by the angle.
  vn(:) = norm * (cos(angle * deg_to_rad) * vn(:) + sin(angle * deg_to_rad) * vn2(:))

  ! The rotated vector is the parallel component plus the rotated component.
  rotation_3d(:) = vp(:) + vn(:)

  return

end function rotation_3d

!************************************************************

function cross_product(v1, v2)

  real(kind=dp), dimension(3), intent(in) :: v1, v2
  real(kind=dp), dimension(3) :: cross_product

  cross_product(1) = v1(2) * v2(3) - v1(3) * v2(2)
  cross_product(2) = v1(3) * v2(1) - v1(1) * v2(3)
  cross_product(3) = v1(1) * v2(2) - v1(2) * v2(1)

  return

end function cross_product

!************************************************************

subroutine progress_bar(j)
  ! progress bar with steps of 2%
  ! j must be between 0 and 50

  integer, intent(in) :: j
  integer :: k
  character(len=58) :: bar

  bar = " ???% |                                                  |"
  write(unit=bar(2:4),fmt="(i3)") 2*j
  do k=1,j
     bar(7+k:7+k)="="
  enddo

  ! print the progress bar.
  write(unit=6,fmt="(a1,a58)",advance="no") char(13), bar
  if (j/=50) then
     flush(unit=6)
  else
     write(unit=6,fmt=*)
  endif

  return

end subroutine progress_bar

!************************************************************

subroutine cdapres(cospsi, phi, u0, v0, w0, u1, v1, w1)
!***************************************************
!*--------COSINUS DIRECTEURS APRES LA DIFFUSION-----
!*
!*        U0,V0 ET W0 ; COSINUS DIRECTEURS AVANT
!*        U1,V1 ET W1 ; COSINUS DIRECTEURS APRES
!*        U ; SELON X
!*        V ; SELON Y    X VERS L'OBSERVATEUR
!*        W ; SELON Z
!*
!*        COSPSI = COSINUS DE L'ANGLE DE DIFFUSION
!*        PHI = AZIMUTH
!*
!*        FRANCOIS MENARD, 21 NOVEMBRE 1988, UDEM
!*
!***************************************************
! 06/12/05 : - passage double car bug sous Icare
!            - reduction nbre d'op�rations
!            (C. Pinte)
!***************************************************

  implicit none

  real(kind=dp), intent(in) :: cospsi, phi, u0, v0, w0
  real(kind=dp), intent(out) :: u1, v1, w1

  real(kind=dp) :: cpsi, spsi, cphi, sphi, a, b,c, Aw0, Cm1

  cpsi = cospsi
  spsi = sqrt(1.0_dp - cpsi*cpsi)
  sphi = sin(phi)
  cphi = cos(phi)

  a = spsi * cphi
  b = spsi * sphi

  ! Calcul de u1,v1 et w1
  if ( abs(w0) <= 0.999999 ) then
     c = sqrt(1.0_dp - w0*w0)
     cm1 = 1.0_dp/c
     aw0 = a*w0
     u1 = ( aw0*u0 - b*v0 ) * cm1 + cpsi*u0
     v1 = ( aw0*v0 + b*u0 ) * cm1 + cpsi*v0
     w1 =  cpsi*w0 - a*c
  else
     u1 = a
     v1 = b
     w1 = cpsi!*w0
  endif

  return

end subroutine cdapres

!***********************************************************

subroutine read_comments(iunit)
  ! Skip comment lines starting with a #

  integer, intent(in) :: iunit
  logical :: comment_found
  character(len=128) :: line

  comment_found = .true.
  do while(comment_found)
     read (iunit, fmt=*) line
     line = adjustl(line)
     if (line(1:1) /= '#') comment_found = .false.
  enddo
  backspace(iunit)

  return

end subroutine read_comments

subroutine read_line(unit,FMT,line,Nread,commentchar)
   !
   !Get next line which is not a comment line nor an empty line
   !return that line and the len of the line Nread

   character(len=*), intent(out) :: line
   character(len=*), optional :: commentchar
   character(len=1) :: com
   integer, intent(out) :: Nread
   integer, intent(in) :: unit
   integer :: EOF
   character(len=*), intent(in) :: FMT

   Nread = 0
   EOF = 0

   com = "#"
   if (present(commentchar)) then
      com = commentchar
   endif

   do while (EOF == 0)
      read(unit, FMT, IOSTAT=EOF) line !'(512A)'
      Nread = len(trim(line))
      !comment or empty ? -> go to next line
      if ((line(1:1).eq.com).or.(Nread==0)) cycle 
      ! line read exit ! to process it
      exit
   enddo ! if EOF > 0 reaches end of file, leaving

   return
end subroutine read_line

function is_nan_infinity(y)
   real(kind=dp), intent(in) :: y
   integer :: is_nan_infinity

   is_nan_infinity = 0
   if (y /= y) then
      write(*,*) "(Nan):", y
      is_nan_infinity = 1
      return
   else if (y > 0 .and. (y==y*10)) then
      write(*,*) "(infinity):", y
      is_nan_infinity = 2
      return
   end if

   return
end function is_nan_infinity

function is_nan_infinity_matrix(y)
    real(kind=dp), intent(in) :: y(:,:)
    integer :: is_nan_infinity_matrix, i, j

    is_nan_infinity_matrix = 0
    do i=1,size(y(:,1))
       do j=1, size(y(1,:))
          if (y(i,j) /= y(i,j)) then
             write(*,*) "(Nan):", y(i,j), " i=", i, " j=",j
             is_nan_infinity_matrix = 1
             return
          else if (y(i,j) > 0 .and. (y(i,j)==y(i,j)*10)) then
             write(*,*) "(infinity):", y(i,j), y(i,j)*(1+0.1), " i=", i, " j=",j
             is_nan_infinity_matrix = 2
             return
          end if
       end do
    end do
    return
end function is_nan_infinity_matrix

function is_nan_infinity_vector(y)
    real(kind=dp), intent(in) :: y(:)
    integer :: is_nan_infinity_vector, i

    is_nan_infinity_vector = 0
    do i=1,size(y)
       if (y(i) /= y(i)) then
          write(*,*) "(Nan):", y(i), " i=", i, y(i)/=y(i)
          is_nan_infinity_vector = 1
          return
       else if (y(i)>0 .and. (y(i)==y(i)*10)) then
          write(*,*) "(infinity):", y(i), y(i)*(1+0.1), " i=", i, (y(i)==y(i)*10)
          is_nan_infinity_vector = 2
          return
       end if
    end do
    return
end function is_nan_infinity_vector

function vacuum2air(Nlambda, lambda_vac) result(lambda_air)
   !wavelength in nm
   integer, intent(in) :: Nlambda
   real(kind=dp), dimension(Nlambda), intent(in) :: lambda_vac
   real(kind=dp), dimension(Nlambda) :: lambda_air
   real(kind=dp), dimension(Nlambda) :: sqwave, reduction

   where (lambda_vac >= VACUUM_TO_AIR_LIMIT)
      sqwave = 1_dp/(lambda_vac**2)
      reduction = 1.0 + 2.735182d-4 + &
           (1.314182 + 2.76249d4 * sqwave) * sqwave
      lambda_air = lambda_vac / reduction
   else where(lambda_vac < VACUUM_TO_AIR_LIMIT)
      lambda_air = lambda_vac
   end where


   return
 end function vacuum2air

 function air2vacuum(Nlambda, lambda_air) result(lambda_vac)
   !wavelength in nm
   integer, intent(in) :: Nlambda
   real(kind=dp), dimension(Nlambda), intent(in) :: lambda_air
   real(kind=dp), dimension(Nlambda) :: lambda_vac
   real(kind=dp), dimension(Nlambda) :: sqwave, increase

   where (lambda_air >= AIR_TO_VACUUM_LIMIT)
      sqwave = (1.0d7 / lambda_air)**2
      increase = 1.0000834213d+00 + &
           2.406030d6/(1.30d10 - sqwave) + &
           1.5997d4/(3.89d9 - sqwave)
      lambda_vac = lambda_air * increase
   else where(lambda_air < AIR_TO_VACUUM_LIMIT)
      lambda_vac = lambda_air
   end where

   return
 end function air2vacuum


function locate(xx,x,mask)
   !wrapper function to locate the position of x in array xx.
   !the closest position is returned.
   real(kind=dp), dimension(:), intent(in) :: xx
   real(kind=dp), intent(in) :: x
   logical, intent(in), dimension(:), optional :: mask
   integer :: locate

   if (present(mask)) then
      locate = minloc((xx-x)**2,1,mask=mask)
   else
      ! 1D array
      locate = minloc((xx-x)**2,1) !(xx(:)-x)*(xx(:)-x)
   end if

   return
end function locate

 function bilinear(N,xi,i0,M,yi,j0,f,xo,yo)
   !bilinear interpolation of the function f(N,M)
   !defined on points xi(N), yi(M) at real xo, real yo.
   !too slow ? f***ck
   real(kind=dp) :: bilinear
   integer, intent(in) :: N, M
   real(kind=dp), intent(in) :: xi(N),yi(M),f(N,M)
   real(kind=dp), intent(in) :: xo,yo 
   integer, intent(in) :: i0, j0
   integer :: i, j
   real(kind=dp) :: norm, f11, f21, f12, f22

   !find closest point in i0 and j0
   ! i0 = max(locate(xi,xo),2)
   ! j0 = max(locate(yi,yo),2)

   ! write(*,*) i0, j0

   norm = ((xi(i0) - xi(i0-1)) * (yi(j0) - yi(j0-1)))
   f11 = f(i0-1,j0-1)
   f21 = f(i0,j0-1)
   f12 = f(i0-1,j0)
   f22 = f(i0,j0)

   bilinear = ((f11 * (xi(i0) - xo) * (yi(j0) - yo) + &
      f21 * (xo - xi(i0-1)) * (yi(j0) - yo) + &
      f12 * (xi(i0) - xo) * (yo - yi(j0-1)) + &
      f22 * (xo - xi(i0-1)) * (yo - yi(j0-1)))) / norm

   return
 end function bilinear

 function linear_1D_sorted(n,x,y, np,xp)
   ! assumes that both x and xp are in increasing order
   ! We only loop once over the initial array, and we only perform 1 test per element
   ! TO DO:
   !  - the bounds are not well handled

   integer, intent(in)                      :: n, np
   real(kind=dp), dimension(n),  intent(in) :: x,y
   real(kind=dp), dimension(np), intent(in) :: xp
   real(kind=dp), dimension(np)             :: linear_1D_sorted

   real(kind=dp) :: t
   integer :: i, j, i0, j0

   linear_1D_sorted(:) = 0.0_dp

   j0=np+1
   do j=1, np
      if (xp(j) > x(1)) then
         j0 = j
         exit
      endif
   enddo
   !write(*,*) "jstar=",j0, xp(j0), x(1)

   ! We perform the 2nd pass where we do the actual interpolation
   ! For points larger than x(n), value will stay at 0
   i0 = 2
   do j=j0, np
      loop_i : do i=i0, n
         if (x(i) > xp(j)) then
            t = (xp(j) - x(i-1)) / (x(i) - x(i-1))
            !write(*,*) j,i, t
            linear_1D_sorted(j) = (1.0_dp - t) * y(i-1)  + t * y(i)
            i0 = i
            exit loop_i
         endif
      enddo loop_i
   enddo


   return
 end function linear_1D_sorted

 function E1 (x)
   !First exponential integral
   real(kind=dp), dimension(6) :: a6
   real(kind=dp), dimension(4) :: a4, b4
   real(kind=dp), intent(in) :: x
   real(kind=dp) :: E1

   a6(:) = 1.0_dp * (/-0.57721566,  0.99999193, -0.24991055,0.05519968, -0.00976004,  0.00107857 /)

   a4(:) = 1.0_dp * (/ 8.5733287401, 18.0590169730, 8.6347608925,  0.2677737343 /)
   b4(:) = 1.0_dp * (/9.5733223454, 25.6329561486, 21.0996530827,  3.9584969228/)

   !Error here, should return an error!
   if (x<=0.0) then
      E1 = 0.0_dp
      return
   else if (x > 0 .and. x <= 1.0_dp) then
      E1 = -log(x) + a6(1) + x*(a6(2) + x*(a6(3) + x*(a6(4) + x*(a6(5) + x*a6(6)))))
   else if (x > 1.0_dp .and. x <= 80.0_dp) then
      E1  = a4(4)/x +  a4(3) + x*(a4(2) + x*(a4(1) + x))
      E1 = E1 / ( b4(4) + x*(b4(3) + x*(b4(2) + x*(b4(1) + x))) )
      E1 = E1 * exp(-x);
   else
      E1 = 0.0_dp
   end if

   return
 end function E1

 function E2(x)
   ! second exponential integral, using recurrence relation
   ! En+1 = 1/n * ( exp(-x) -xEn(x))
   real(kind=dp), intent(in) :: x
   real(kind=dp) :: E2

   if (x <= 0.0) then
      E2 = 0.0_dp
      return
   endif

   if (x > 0.0 .and. x <= 80.0) then
      E2 = 1.0_dp * (exp(-x) - x * E1(x))
   else
      E2 = 0.0_dp
   end if

   return
 end function E2

end module utils
