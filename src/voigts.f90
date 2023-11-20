! --------------------------------------------------------------------------- !
! Implements Voigt function procedures
! --------------------------------------------------------------------------- !
module voigts

    use constantes, only : pi, sqrtpi
    use mcfost_env, only : dp
    
    implicit none

    real(kind=dp), dimension(7) :: Aj, Bj, Cj
    data Aj /122.607931777104326,214.382388694706425,181.928533092181549,&
        93.155580458138441, 30.180142196210589,5.912626209773153,0.564189583562615/
    data Bj /122.607931773875350,352.730625110963558,457.334478783897737,&
        348.703917719495792,170.354001821091472,53.992906912940207,10.479857114260399/
    data Cj /0.5641641,0.8718681,1.474395,-19.57826,802.4513,-4850.316,8031.468/

    !At the moment, the procedure does not work because the different functions
    !do not take the same parameters as input.
    procedure(VoigtHumlicek), pointer :: Voigt => VoigtHumlicek

    contains

    function VoigtHumlicek(N, a, v, F)
    !--------------------------------------------------------------
    ! Generates Voigt and anomalous dispersion profiles
    ! See Humlicek (1982) JQSRT 27, 437
    ! W4 = Voigt + i*Faraday-Voigt (dispersion profile)
    !--------------------------------------------------------------
        integer, intent(in) :: N
        real(kind=dp), intent(in) :: a, v(N)
        real(kind=dp) :: VoigtHumlicek(N)
        real(kind=dp), intent(out), optional :: F(N)
        complex(kind=dp) :: w4, z, u, v4
        real(kind=dp) :: s
        integer :: i

        do i = 1, n

            z = cmplx(a,-v(i))
            s = abs(v(i)) + a
            u = z*z

            if (s >= 15.d0) then
          ! Region I
                w4 = z * 0.5641896d0 / (0.5d0+u)
            elseif (s >= 5.5) then
          ! Region II
                w4 = z*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
            elseif (a >= 0.195d0*abs(v(i))-0.176d0) then
          ! Region III
                w4 = (16.4955d0+z*(20.20933d0+z*(11.96482d0+z*(3.778987d0+z*0.5642236d0)))) / &
                    (16.4955d0+z*(38.82363d0+z*(39.27121d0+z*(21.69274d0+z*(6.699398d0+z)))))
            else
          ! Region IV
                w4 = z*(36183.31d0-u*(3321.9905d0-u*(1540.787d0-u*(219.0313d0-u*(35.76683d0-&
                    u*(1.320522d0-u*0.56419d0))))))
                v4 = (32066.6d0-u*(24322.84d0-u*(9022.228d0-u*(2186.181d0-u*(364.2191d0-&
                    u*(61.57037d0-u*(1.841439d0-u)))))))
                w4 = exp(u) - w4/v4 !cdexp, can be optimized by splitting the exp on real and imag parts
            endif
            VoigtHumlicek(i) = real(w4,kind=dp)!w4%re
            if (present(F)) then 
                F(i) = aimag(w4)!w4%im
            endif

        enddo
    
        return
    end function VoigtHumlicek

 
    function VoigtThomson(N, a, v, vbroad)
        integer, intent(in) :: N
        real(kind=dp), intent(in) :: a, vbroad, v(N)
        real(kind=dp) :: VoigtThomson(N)
        integer :: i
        real(kind=dp) :: f, rr, eta, ap

        ap = a * vbroad

        f = (vbroad**5. + 2.69269*vbroad**4. * ap + 2.42843*vbroad**3. * ap**2. + &
         4.47163*vbroad**2. *ap**3. + 0.07842*vbroad*ap**4. + ap**5.)**(0.2)

        rr = ap/f
        eta = 1.36603*rr - 0.47719*rr*rr + 0.11116*rr*rr*rr

        ! VoigtThomson = eta * ( f/pi / ((vbroad*v)**2 + f**2) ) + &
        !     (1.0_dp - eta) * exp(-(v*vbroad/f)**2) / f / sqrtpi

        !because Thomson is in units of 1/vbroad/sqrtpi already
        VoigtThomson = eta * ( vbroad * f / pi / sqrtpi / ((vbroad*v)**2 + f**2) ) + &
            (1.0_dp - eta) * exp(-(v*vbroad/f)**2) * vbroad / f
    
        return
    end function VoigtThomson

   function max_voigt_profile(vth,a)
   !Evaluate the maximum of the voigt function, using the Thomson approximation.
      real(kind=dp) :: max_voigt_profile
      real(kind=dp), intent(in) :: vth, a
      real(kind=dp) :: f, ap, eta

      ap = vth * a
      f = (vth**5 + 2.69269*vth**4 * ap + 2.42843*vth**3 * ap**2 + &
         4.47163*vth**2*ap**3 + 0.07842*vth*ap**4 + ap**5)**(1./5.)
	
      eta = 1.36603*(ap/f) - 0.47719*(ap/f)**2+0.11116*(ap/f)**3
      max_voigt_profile = eta / f / pi + (1.0-eta) / sqrt(pi) / f

      return 
   end function max_voigt_profile

   function dmax_voigt(vth,a,eps)
   !Evaluate the distance in units of vth, at which the voigt profile
   ! reaches "eps * d" value : the value of the voigt profile in units of 
   ! its peak. The Thomson approximation is used.
    real(kind=dp) :: dmax_voigt
    real, intent(in) :: eps
    real(kind=dp), intent(in) :: vth, a
    real(kind=dp) :: peak_frac
    real(kind=dp) :: f, ap, eta

    ap = vth * a
    peak_frac = max_voigt_profile(vth,a) * eps

	f = (vth**5 + 2.69269*vth**4 * ap + 2.42843*vth**3 * ap**2 + &
	    4.47163*vth**2*ap**3 + 0.07842*vth*ap**4 + ap**5)**(1./5.)
	
	eta = 1.36603*(ap/f) - 0.47719*(ap/f)**2+0.11116*(ap/f)**3
	
	dmax_voigt = eta * sqrt(f/pi/peak_frac - f**2) + (1.0 - eta) * f * sqrt(-log(sqrt(pi)*f*peak_frac))
    ! write(*,*) "eta=", eta
    ! write(*,*) "xp(L) = ", sqrt(f/pi/peak_frac - f**2)
    ! write(*,*) "xp(G) = ", f * sqrt(-log(sqrt(pi)*f*peak_frac))

    return
   end function dmax_voigt

   function dmax_lorentz(vth, a, eps)
   !assuming at large x, the Voigt profile collapse to a Lorentzian
   !even for small a
    real(kind=dp) :: dmax_lorentz
    real(kind=dp), intent(in) :: vth, a
    real, intent(in) :: eps
    real(kind=dp) :: peak_frac

    !note the peak of a Voigt and not a Lorentzian
    peak_frac = max_voigt_profile(vth,a) * eps

    dmax_lorentz = vth * sqrt(a/sqrtpi/vth/peak_frac - a*a)

    return
   end function dmax_lorentz

    !building
    function VoigtAller(N,a,v)
        integer, intent(in) :: N
        real(kind=dp), intent(in) :: a, v(N)
        real(kind=dp) :: VoigtAller(N)
        integer :: i,j
        real(kind=dp) :: sum1, base
        complex(kind=dp) :: ksi, sum2
                
        if (a <= 0.1) then
            do i=1,N
                if (abs(v(i)) >= 2.5) then
                    base = (1.0+a*a*(1-2*v(i)**2))*exp(-v(i)**2)
                    sum1 = 0
                    do j=1,7
                        sum1 = sum1 + cj(j)*v(i)**(-2*(j-1))
                    enddo
                    VoigtAller(i) = a/v(i)/v(i) * sum1 + base
                else
                    ksi = cmplx(a, -v(i))
                    sum1 = 0
                    sum2 = ksi**7
                    do j=1,7
                        sum1 = sum1 + Aj(j)*ksi**(j-1)
                        sum2 = sum2 + Bj(j)*ksi**(j-1) 
                    enddo
                    VoigtAller(i) = real(sum1/sum2, kind=dp)
                endif	
            enddo
        else
            do i=1,N
                ksi = cmplx(a, -v(i))
                sum1 = 0
                sum2 = ksi**7
                do j=1,7
                    sum1 = sum1 + Aj(j)*ksi**(j-1)
                    sum2 = sum2 + Bj(j)*ksi**(j-1) 
                enddo
                VoigtAller(i) = real(sum1/sum2, kind=dp)
            enddo
        endif
    
        return
    end function VoigtAller

end module voigts