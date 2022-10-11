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
            s = dabs(v(i)) + a
            u = z*z

            if (s >= 15.d0) then
          ! Region I
                w4 = z * 0.5641896d0 / (0.5d0+u)
            elseif (s >= 5.5) then
          ! Region II
                w4 = z*(1.410474d0+u*0.5641896d0)/(0.75d0+u*(3.d0+u))
            elseif (a >= 0.195d0*dabs(v(i))-0.176d0) then
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
        real(kind=dp) :: aeff, ratio, eta, al

        al = a * vbroad

        aeff = (vbroad**5. + 2.69269*vbroad**4. * aL + 2.42843*vbroad**3. * aL**2. + &
         4.47163*vbroad**2. *aL**3. + 0.07842*vbroad*aL**4. + aL**5.)**(0.2)



        ratio = aL/aeff
        eta = 1.36603*ratio - 0.47719*(ratio*ratio) + 0.11116*(ratio*ratio*ratio)


        VoigtThomson = eta * ( aeff/pi * (v**2 + aeff**2)**(-1.0) ) + &
         (1.0_dp - eta) * exp(-(v/aeff)**2) / aeff / sqrtpi


        return
    end function VoigtThomson

    
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