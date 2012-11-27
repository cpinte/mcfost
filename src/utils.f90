module utils

  use parametres
  use naleat, only : stream, gauss_random_saved, lgauss_random_saved
  use constantes
  use sha

  implicit none

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

!************************************************************

function spanl(xmin,xmax,n)

  implicit none

  real, intent(in) :: xmin, xmax
  integer, intent(in) :: n
  real, dimension(n) :: spanl

  spanl = exp(span(log(abs(xmin)), log(abs(xmax)), n))

end function spanl

!************************************************************

function gauss_random(id)
  ! Retourne un nombre aleatoire distribue suivant une gaussienne
  ! avec moyenne nulle et variance 1
  ! Adapte de numerical recipes pour fonctionner en parallel avec sprng
  ! C. Pinte
  ! 19/10/07

  implicit none

#include "sprng_f.h"

  integer, intent(in) :: id
  real(kind=db) :: gauss_random
  real(kind=db) :: rsq,rand1,rand2

  if (lgauss_random_saved(id)) then
     gauss_random=gauss_random_saved(id)
     lgauss_random_saved(id)=.false.
  else
     loop : do
        rand1 = sprng(stream(id))
        rand2 = sprng(stream(id))
        rand1 = 2.0_db * rand1 - 1.0_db
        rand2 = 2.0_db * rand2 - 1.0_db
        rsq=rand1**2+rand2**2
        if (rsq > 0.0_db .and. rsq < 1.0_db) exit loop
     enddo loop
     rsq=sqrt(-2.0_db*log(rsq)/rsq)
     gauss_random = rand1 * rsq
     gauss_random_saved(id)= rand2 * rsq
     lgauss_random_saved(id)=.true.
  endif
	
  return

end function gauss_random

!**********************************************************

subroutine polint(xa,ya,n,x,y,dy)
! Interpolation polynomiale
! xa, ya : tab des abscisses et ordonnees
! n : degre du polynome + 1 = taille xa, ya 
! x : abscisse du point considere
! y, dy : valeur et erreur
  implicit none

  integer, intent(in) :: n
  integer, parameter :: nmax = 10
  real, intent(in) :: x
  real, dimension(n) , intent(in) :: xa, ya
  real, intent(out) :: y, dy

  integer :: i,m,ns
  real :: den,dif,dift,ho,hp,w
  real, dimension(nmax) :: c,d

  ns=1
  dif=abs(x-xa(1))
  do i=1,n
     dift=abs(x-xa(i))
     if (dift.lt.dif) then
        ns=i
        dif=dift
     endif
     c(i)=ya(i)
     d(i)=ya(i)
  enddo
  y=ya(ns)
  ns=ns-1
  do m=1,n-1
     do i=1,n-m
        ho=xa(i)-x
        hp=xa(i+m)-x
        w=c(i+1)-d(i)
        den=ho-hp
        if(den.eq.0.) then
           write(*,*) 'failure in polint'
           stop
        endif
        den=w/den
        d(i)=hp*den
        c(i)=ho*den
     enddo !i
     if (2*ns.lt.n-m)then
        dy=c(ns+1)
     else
        dy=d(ns)
        ns=ns-1
     endif
     y=y+dy
  enddo !m
  return
end subroutine polint

!***************************************************

real(kind=sl) function interp_sp(y, x, xp)
! interpole lineaire entre les points
! fait une extrapolation lineaire entre les 2 premiers ou 2 derniers points
! en cas de point cp en dehors de x
! Modif : pas d'extrapolation : renvoie valeur mini en dehors de l'intervalle
! C. Pinte
! 26/01/07
    
  implicit none
  
  real(kind=sl), dimension(:), intent(in) :: x, y
  real(kind=sl) :: xp

  integer :: n, np, j, ny
    
  real :: frac

  n=size(x) !; ny=size(y)
    
  !if (n /= ny) then
  !   write(*,*) "Error in interp : y and x must have same dim"
  !endif


  ! PAS D'EXTRAPOLATION
  if ((xp < minval(x)) .or. (xp > maxval(x))) then
     interp_sp = minval(y) 
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

real(kind=db) function interp_dp(y, x, xp)
! interpole lineaire entre les points
! fait une extrapolation lineaire entre les 2 premiers ou 2 derniers points
! en cas de point cp en dehors de x
! Modif : pas d'extrapolation : renvoie valeur mini en dehors de l'intervalle
! C. Pinte
! 26/01/07
    
  implicit none
  
  real(kind=db), dimension(:), intent(in) :: x, y
  real(kind=db) :: xp

  integer :: n, np, j, ny
    
  real :: frac

  n=size(x) !; ny=size(y)
    
  !if (n /= ny) then
  !   write(*,*) "Error in interp : y and x must have same dim"
  !endif


  ! PAS D'EXTRAPOLATION
  if ((xp < minval(x)) .or. (xp > maxval(x))) then
     interp_dp = minval(y) 
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
  ! C'est bourrin de chez bourrin : meme pas de recherche du pivot max !!!
  ! C. Pinte
  ! 22/09/07

  implicit none
  
  integer, intent(in) :: n 
  real(kind=db), dimension(n,n), intent(inout) :: a
  real(kind=db), dimension(n), intent(inout) :: b    

  real(kind=db) :: factor
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
  ! le systeme de coordonnees où le vecteur (u1,v1,w1)=(1,0,0).
  ! ie applique la meme rotation que celle qui transforme
  ! (1,0,0) en (u1,v1,w1) sur (xinit,yinit,zinit)
  ! C. Pinte : 1/03/06
  ! Nouvelle version d'une routine de F. Ménard

  implicit none

  real(kind=db), intent(in) ::  xinit, yinit, zinit, u1, v1, w1
  real(kind=db), intent(out) :: xfin, yfin, zfin

  real(kind=db) :: cost, sint, sing, prod, x, theta

  if (w1 > 0.999999999_db) then 
     cost = 1.0_db
     sint = 0.0_db
     sing = 0.0_db
  else
     if (abs(u1) < tiny_real) then
        cost=0.0_db
        sint=1.0_db
        sing=sqrt(1.0 - w1*w1)
     else
   !     x=v1/u1
   !     !c'est pas un atan mais un atan2 d'on le sign(u1)
   !     cost=sign(1.0/sqrt(1.0+x*x),u1)        !cos(atan(x)) = 1/(1+sqrt(x*x))
   !     sint=x*cost                            !sin(atan(x)) = x/(1+sqrt(x*x))
        ! Equivalent  à  (~ meme tps cpu):
         theta=atan2(v1,u1)
         cost=cos(theta)
         sint=sin(theta)
         sing = sqrt(1.0_db - w1*w1)
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

  real(kind=db), intent(in) :: mu, a
  real(kind=db) :: calc_mu0

  real(kind=db) :: a1, a2, a3, q, r, q3, qh, denom, theta, mu01, mu02, mu03, factor, factor3


  a1 = 0.0_db
  a2 = a - 1.0_db
  a3 = -a * mu

  q = (a1**2 - 3.0_db*a2)/9.0_db
  r = (2.0_db*a1**3 - 9.0_db*a1*a2 + 27.0_db*a3)/54.0_db
  q3 = q**3
  
  if ((q3 - r**2) >= 0.0_db) then
     qh = sqrt(q)
     denom = sqrt(q3)
     theta = acos(r/denom)
     mu01 = -2.0_db*qh*cos(theta/3.0_db) - a1/3.0_db
     mu02 = -2.d0*qh*cos((theta + deux_pi)/3.0_db) - a1/3.0_db
     mu03 = -2.d0*qh*cos((theta + quatre_pi)/3.0_db) - a1/3.0_db
     if ((mu01 > 0.0_db).and.(mu02 > 0.0_db)) then
        write(*,*) "ERREUR: calc_mu0: 2 racines: mu01,mu02",mu01,mu02
        stop
     else if ((mu01 > 0.0_db).and.(mu03 > 0.0_db)) then
        write(*,*) "ERREUR: calc_mu0: 2 racines: mu01,mu03",mu01,mu03
        stop 
     elseif ((mu02 > 0.0_db).and.(mu03 > 0.0_db)) then
        write(*,*) "ERREUR: calc_mu0: 2 racines: mu02,mu03",mu02,mu03
        stop
     endif

     if (mu01 > 0.0_db) calc_mu0 = mu01
     if (mu02 > 0.0_db) calc_mu0 = mu02
     if (mu03 > 0.0_db) calc_mu0 = mu03
 
  else
     factor = sqrt(r**2 - q3) + abs(r)
     factor3 = factor**un_tiers
     calc_mu0 = -1.0_db*sign(1.0_db,r)*(factor3 + q/factor3) - a1/3.0_db   
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

  real(kind=db) ,intent(in) :: nu
  real, intent(in) :: T
  real(kind=db) :: Bnu

  real(kind=db) :: hnu_kT

  hnu_kT = (hp * nu) / (kb * T)
 
  if (hnu_kT > 100._db) then
     Bnu = 0.0_db
  else
     Bnu = 2.0_db*hp/c_light**2 * nu**3 / (exp(hnu_kT)-1.0_db)
  endif
     
  return

end function Bnu

!***********************************************************

function Blambda(wl,T)
! Loi de Planck 
! Blambda en SI : W.m-2.s-1.sr-1
! wl en m
! C. Pinte
! 23/05/09

  implicit none

  real ,intent(in) :: wl
  real, intent(in) :: T
  real :: Blambda

  real(kind=db) :: hnu_kT

  hnu_kT = (hp * c_light)/ (kb * T * wl)
 
  if (hnu_kT > 100.) then
     Blambda = 0.0
  else
     Blambda = 2.0*hp*c_light**2 / wl**5 / (exp(hnu_kT)-1.0)
  endif
     
  return

end function Blambda

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

subroutine mcfost_update(lforce_update)

  logical, intent(in) :: lforce_update
  logical :: lupdate

  character(len=512) :: cmd, url, last_version, machtype, ostype, system
  integer ::  syst_status, ios

  write(*,*) "Version ", mcfost_release
  
  ! Last version
  write(*,*) "Checking last version ..."
  !cmd = "wget "//trim(webpage)//"/version.txt -q -T 5 -t 3"
  cmd = "curl "//trim(webpage)//"version.txt -O -s"
  call appel_syst(cmd, syst_status)           
  if (syst_status/=0) then 
     write(*,*) "ERROR: Cannot get MCFOST last version number."
     write(*,*) "Current version file not found."
     write(*,*) "Exiting"
     stop
  endif

  open(unit=1, file="version.txt", status='old',iostat=ios)
  read(1,*,iostat=ios) last_version
  close(unit=1,status="delete",iostat=ios)
  if ( (ios/=0) .or. (.not.is_digit(last_version(1:1)))) then 
     write(*,*) "ERROR: Cannot get MCFOST last version number."
     write(*,*) "Cannot read version file."
     write(*,*) "Exiting"
     stop
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
     ! get system info
     call get_environment_variable('OSTYPE',ostype)
     call get_environment_variable('MACHTYPE',machtype)
     system = trim(ostype)//" "//trim(machtype)

     
     ! get the correct url corresponding to the system
     if (ostype(1:5)=="linux") then
        if (machtype(1:4)=="i386") then
           url = trim(webpage)//"linux_32bits/mcfost"
        else if (machtype(1:6)=="x86_64") then
           url = trim(webpage)//"linux_64bits/mcfost"
        else
           write(*,*) "Linux system but unknown architecture"
           write(*,*) "Cannot download new binary"
           write(*,*) "Exiting."
           stop
        endif
         system = trim(ostype)//" "//trim(machtype)
     else if (ostype(1:6)=="darwin") then
        system = "MacOS X x86_64" 
        url = trim(webpage)//"macos_intel_64bits/mcfost"
     else
        write(*,*) "Unknown operating system"
        write(*,*) "Cannot download new binary"
        write(*,*) " "
        write(*,*) "Detected configuration:"
        write(*,*) "OSTYPE: ", trim(ostype)
        write(*,*) "MACHTYPE: ", trim(machtype)
        write(*,*) " "
        write(*,*) "You can try to update the environement variables OSTYPE & MACHTYPE"
        write(*,*) "Known values:"
        write(*,*) " - OSTYPE: linux or darwin"
        write(*,*) " - MACHTYPE: i386 or x86_64"
        
        write(*,*) "Exiting."
        stop
     endif

     write(*,*) "Your system ", trim(system)
     
     
     ! Download
     write(*,'(a32, $)') "Downloading the new version ..."
     !cmd = "wget -q "//trim(url)//" -O mcfost_update"
     cmd = "curl "//trim(url)//" -o mcfost_update -s"
     call appel_syst(cmd, syst_status)  
     if (syst_status==0) then 
        write(*,*) "Done"
        cmd = "chmod a+x mcfost_update ; mv mcfost mcfost_"//trim(mcfost_release)//" ; mv mcfost_update mcfost"
        call appel_syst(cmd, syst_status)    
        write(*,*) " "
        write(*,*) "MCFOST has been updated"
        write(*,*) "The previous version has been saved as mcfost_"//trim(mcfost_release)
     else
        cmd = "rm -rf mcfost_update"
        call appel_syst(cmd, syst_status) 
        write(*,*) "Download error"
        write(*,*) "Exiting"
     endif
        
  endif ! lupdate

  return

end subroutine mcfost_update

!***********************************************************

subroutine mcfost_history()

  character(len=512) :: cmd
  integer ::  syst_status
  
  ! Last version
  write(*,*) "Getting MCFOST history ..."
  cmd = "curl "//trim(webpage)//"history.txt"
  call appel_syst(cmd, syst_status)           
  if (syst_status/=0) then 
     write(*,*) "Cannot get MCFOST history"
     write(*,*) "Exiting"
     stop
  endif
  write(*,*) " "

  return

end subroutine mcfost_history

!***********************************************************

subroutine mcfost_get_ref_para()

  character(len=512) :: cmd
  character(len=12) :: ref_file
  integer ::  syst_status
  
  ref_file = "ref"//mcfost_release(1:4)//".para"
  
  write(*,*) "Getting MCFOST reference file: "//ref_file
  cmd = "curl "//trim(webpage)//ref_file//" -O -s"
  call appel_syst(cmd, syst_status)           
  if (syst_status/=0) then 
     write(*,*) "Cannot get MCFOST reference file"
     write(*,*) "Exiting"
     stop
  endif
  write(*,*) "Done"

  return

end subroutine mcfost_get_ref_para

!***********************************************************

subroutine mcfost_v()

  character(len=512) :: cmd, last_version
  integer ::  syst_status, line_number, ios
  character(len=128) :: sline_number

  write(*,*) "This is MCFOST version: ", mcfost_release
  write(*,fmt="(A24, F5.2)") "Parameter file version ", mcfost_version
  write(*,*) "SHA-1 id = ", sha_id
  !write(*,*) 1.0/0.0 , __LINE__, __FILE__
  write(*,*) "Binary compiled the ",__DATE__," at ",__TIME__
#if defined (__INTEL_COMPILER)
  write(*,fmt='(" with INTEL compiler version ",i4)')  __INTEL_COMPILER
#endif
#if defined (__GFORTRAN__)
  write(*,fmt='(" with GFORTAN compiler version ",i1,".",i1,"."i1)') __GNUC__ , __GNUC_MINOR__,  __GNUC_PATCHLEVEL__
#endif
#if defined (__G95__)
  write(*,fmt='(" with G95 compiler version ",i1,".",i2)') __G95__, __G95_MINOR__
#endif
  write(*,*) " "
  

  ! Last version
  write(*,*) "Checking last version ..."
  cmd = "curl "//trim(webpage)//"version.txt -O -s"

  call appel_syst(cmd, syst_status)           
  if (syst_status/=0) then 
     write(*,*) "ERROR: Cannot get MCFOST last version number (Error 1)"
     write(*,*) "Exiting"
     stop
  endif
  open(unit=1, file="version.txt", status='old',iostat=ios)
  read(1,*,iostat=ios) last_version
  close(unit=1,status="delete",iostat=ios)

  if ( (ios/=0) .or. (.not.is_digit(last_version(1:1)))) then 
     write(*,*) "ERROR: Cannot get MCFOST last version number (Error 2)"
     write(*,*) "Exiting"
     stop
  endif

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
     if (syst_status/=0) then 
        write(*,*) "Cannot get MCFOST history"
        write(*,*) "Exiting"
        stop
     endif

     ! Getting line number of current version in history
     cmd = "grep -n "//trim(mcfost_release)//" history.txt | awk -F : '{print $1}' > line_number.txt"
     call appel_syst(cmd, syst_status)
     open(unit=1, file="line_number.txt", status='old',iostat=ios)
     read(1,*,iostat=ios) line_number
     close(unit=1,status="delete",iostat=ios) 
     
     ! Printing history
     if (ios==0) then
        write(sline_number,*)  line_number-1
        write(*,*) "Changes since your version:"
        cmd = "head -"//trim(adjustl(sline_number))//" history.txt ; rm -rf history.txt"
        call appel_syst(cmd, syst_status)

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
  if (syst_status/=0) then 
     write(*,*) "ERROR: Cannot get MCFOST UTILS last version number (Error 1)"
     write(*,*) "Exiting"
     stop
  endif
  open(unit=1, file="Version", status='old',iostat=ios)
  read(1,*,iostat=ios) s_last_version
  close(unit=1,status="delete",iostat=ios)

  if ( (ios/=0) .or. (.not.is_digit(s_last_version(1:1)))) then 
     write(*,*) "ERROR: Cannot get MCFOST UTILS last version number (Error 2)"
     write(*,*) "Exiting"
     stop
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
     write(*,*) "New MCFOST version available: ", trim(s_last_version)
     write(*,*) "Please update mcfost first with mcfost -u"
     lupdate = .false.
  endif

  if (lforce_update) then
     write(*,*) "Forcing update"
     lupdate = .true.
  endif

  if (lupdate) then
     write(*,*) "Updating MCFOST UTILS ..."
     call get_utils()
  endif

  return

end subroutine update_utils

!***********************************************************

subroutine get_utils()

  character(len=512) :: cmd
  integer :: syst_status

  cmd = "curl "//trim(utils_webpage)//"mcfost_utils.tgz -s | tar xzvf - -C"//trim(mcfost_utils)
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
  read(1,*,iostat=ios) get_mcfost_utils_version
  close(unit=1,iostat=ios)
  if (ios /= 0) get_mcfost_utils_version = 0.0

  if (abs(get_mcfost_utils_version) < 1e-6) then
     write(*,*) "ERROR : could not find current MCFOST utils version"
     write(*,*) "Exiting"
     stop
  endif

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

function bubble_sort(data_in)
  ! Implementation of Bubble sort
  ! Warning : this is N^2, only for small arrays 
  ! Same behaviour as yorick to allow ordering of mutiple arrays
  ! Return the order of data, the sorted array would be data_in(order)
  !
  ! C. Pinte
  ! 02/05/11

  real(kind=db), dimension(:), intent(in) :: data_in   
  real(kind=db), dimension(:), allocatable :: data 
  integer, dimension(:), allocatable :: bubble_sort ! indices

  integer :: i, pass, n, tmp_i
  logical :: sorted
  real ::  temp

  n = size(data_in)
  allocate(bubble_sort(n),data(n)) 
  data = data_in

  bubble_sort = indgen(n) ;

  pass = 1
  sorted = .false.
  
  do while(.not.sorted) 
     sorted = .true.
     
     do i = 1,n-pass
        if(data(i) > data(i+1)) then
           temp = data(i)
           data(i) = data(i+1)
           data(i+1) = temp
           sorted = .false.

           ! same operation on indices
           tmp_i = bubble_sort(i) 
           bubble_sort(i) = bubble_sort(i+1)
           bubble_sort(i+1) = tmp_i
           
        endif
     end do ! i
     pass = pass +1

  end do ! while

  return

end function bubble_sort

!************************************************************

function is_diff(a,b)

  logical is_diff
  real(kind=db), intent(in) :: a,b
  
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

  ! 29/08/11 : C. Pinte : translated in fortran 90 TODO : regarder NR

  integer, intent(in) :: lmax
  real(kind=db), intent(in) :: x
  real(kind=db), dimension(0:lmax), intent(out) :: P

  integer :: l

  P(0) = 1._db
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

subroutine appel_syst(cmd, status)

#if defined (__INTEL_COMPILER)
  use IFPORT, only : system
#endif

  implicit none

  character(len=*), intent(in) :: cmd
  integer, intent(out), optional :: status

  status = system(cmd) ! limux
  !call system(cmd, status) ! aix

  return

end subroutine appel_syst

end module utils
