module mcfost2ProDiMo

  use constantes, only : pi
  
  implicit none

  ! Spatial grid of mcfost
  integer, parameter :: n_rad = 100
  integer, parameter :: nz = 40
  ! Wavelength grid of mcfost
  integer, parameter :: n_lambda = 39

  ! Data file
  character(len=512), parameter :: filename = "data_th/forProDiMo.fits.gz"

  type mcfost_model
     
     ! parameters
     real :: diskmass, Teff, Rstar, Mstar, fUV, fUV_slope, rin, rout, amin, amax, aexp, rho_grain, a_settle

     ! Saptial grid (Unit: AU)
     real, dimension(n_rad,nz) :: r_grid, z_grid
  
     ! Temperature structure (Unit: K)
     real, dimension(n_rad,nz) :: temperature

     ! Wavelength bins (Unit: microns)
     real, dimension(n_lambda) :: wavelengths

     ! input Stellar spectrum (lambda.F_lambda) (Units : W.m^-2)
     ! sum(F_lambda x dlambda = sigma Teff^4 in case of 1 star without UV excess)
     real, dimension(n_lambda) :: J_star_input

     ! input ISM spectrum (lambda.F_lambda) (Units : W.m^-2)
     real, dimension(n_lambda) :: J_ISM_input

     ! Radiation field
     ! lambda.J_lambda (Units : W.m^-2)
     real, dimension(n_rad,nz,n_lambda) :: J, J_ISM
     
     ! Statictic of the radiation field (number of packets)
     real, dimension(n_rad,nz,n_lambda) :: nJ, nJ_ISM
       
     ! Gas density assuming a gas/dust mass ratio of 100 (Unit: part.m-3)
     real, dimension(n_rad,nz) :: density

     ! 0th, 1st, 2nd and 3rd moment of the grain size distribution
     ! 0th moment: unit = part.m-3
     ! 1th, 2nd and 3rd moment: units = micron, micron^2, micron^3
     real, dimension(n_rad,nz,0:3) :: grain_size

     ! Dust opacities (Units : AU^-1)
     real, dimension(n_rad,nz,n_lambda) :: kappa_ext, kappa_abs
  end type mcfost_model

contains

  function read_mcfost() result(mcfost)

    integer :: status, readwrite, unit, blocksize, nfound, group, firstpix, nbuffer, npixels, hdunum, hdutype
    real :: nullval
    logical :: anynull

    integer, dimension(4) :: naxes

    character(len=30) errtext
    character (len=80) errmessage, comment

    real, dimension(n_rad,nz,2) :: grid
    real, dimension(n_rad,nz,2,n_lambda) :: opacity
    
    type(mcfost_model) :: mcfost


    !  Get an unused Logical Unit Number to use to open the FITS file.
    status=0
    call ftgiou (unit,status)

    ! Open file
    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status /= 0) then ! le fichier temperature n'existe pas
       write(*,*) "ERROR : ProDiMo.fits.gz file needed"
       stop
    endif

    group=1
    firstpix=1
    nullval=-999

    !------------------------------------------------------------------------------
    ! Spatial grid
    !------------------------------------------------------------------------------
    ! hdu 1

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 1 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= 2)) then
       write(*,*) "Error : HDU 1 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! Keywords
    ! FTGKY[EDJKLS](unit,keyword, > keyval,comment,status)
    call ftgkye(unit,'disk_dust_mass',mcfost%diskmass,comment,status)

    call ftgkye(unit,'Teff',mcfost%Teff,comment,status)
    call ftgkye(unit,'Rstar',mcfost%Rstar,comment,status)
    call ftgkye(unit,'Mstar',mcfost%Mstar,comment,status)
    call ftgkye(unit,'fUV',mcfost%fUV,comment,status)
    call ftgkye(unit,'slope_UV',mcfost%slope_UV,comment,status)
 
    call ftgkye(unit,'Rin',mcfost%rin,comment,status)
    call ftgkye(unit,'Rout',mcfost%rout,comment,status)

    call ftgkye(unit,'amin',mcfost%amin,comment,status)
    call ftgkye(unit,'amax',mcfost%amax,comment,status)
    call ftgkye(unit,'aexp',mcfost%aexp,comment,status)
    call ftgkye(unit,'a_settle',mcfost%a_settle,comment,status)
    call ftgkye(unit,'rho_grain',mcfost%rho_grain,comment,status)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,grid,anynull,status)
    mcfost%r_grid(:,:) = grid(:,:,1)
    mcfost%z_grid(:,:) = grid(:,:,2)
    
    !------------------------------------------------------------------------------
    ! HDU 2: Temperature 
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)    

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    if (nfound /= 2) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 2 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
       write(*,*) "Error : HDU 2 does not have the"
       write(*,*) "right dimensions. Exiting."
       write(*,*) naxes(1), n_lambda
       stop
    endif
    npixels=naxes(1)*naxes(2)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%temperature,anynull,status)


    !------------------------------------------------------------------------------
    ! HDU 3: Wavelengths
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)
    

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 3 file. Exiting.'
       stop
    endif
    if (naxes(1) /= n_lambda) then
       write(*,*) "Error : HDU 3 does not have the"
       write(*,*) "right dimensions. Exiting."
       write(*,*) naxes(1), n_lambda
       stop
    endif
    npixels=naxes(1)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%wavelengths,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 4 : Stellar input radiation field
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)
    

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 3 file. Exiting.'
       stop
    endif
    if (naxes(1) /= n_lambda) then
       write(*,*) "Error : HDU 3 does not have the"
       write(*,*) "right dimensions. Exiting."
       write(*,*) naxes(1), n_lambda
       stop
    endif
    npixels=naxes(1)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%J_star_input,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 5 : ISM input radiation field
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)
    

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,1,naxes,nfound,status)
    if (nfound /= 1) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 3 file. Exiting.'
       stop
    endif
    if (naxes(1) /= n_lambda) then
       write(*,*) "Error : HDU 3 does not have the"
       write(*,*) "right dimensions. Exiting."
       write(*,*) naxes(1), n_lambda
       stop
    endif
    npixels=naxes(1)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%J_ISM_input,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 6 : Radiation fied
    !------------------------------------------------------------------------------  
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 4 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= n_lambda)) then
       write(*,*) "Error : HDU 4 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%J,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 7 : Statistic of the radiation fied
    !------------------------------------------------------------------------------  
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 4 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= n_lambda)) then
       write(*,*) "Error : HDU 4 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%nJ,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 8 : ISM radiation fied
    !------------------------------------------------------------------------------  
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 4 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= n_lambda)) then
       write(*,*) "Error : HDU 4 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%J_ISM,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 9 : Statistic of the ISM radiation fied
    !------------------------------------------------------------------------------  
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 4 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= n_lambda)) then
       write(*,*) "Error : HDU 4 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%nJ_ISM,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 10 : Gaz density for a gas to dust mass ratio of 100
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
    if (nfound /= 2) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 5 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
       write(*,*) "Error : HDU 5 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%density,anynull,status)

    !------------------------------------------------------------------------------
    ! HDU 11 : Opacities
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
    if (nfound /= 4) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 6 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= 2).or.(naxes(4) /= n_lambda)) then
       write(*,*) "Error : HDU 6 does not have the"
       write(*,*) "right dimensions. Exiting."
  !     write(*,*) naxes(1:4)
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,opacity,anynull,status)
    mcfost%kappa_ext(:,:,:) = opacity(:,:,1,:)
    mcfost%kappa_abs(:,:,:) = opacity(:,:,2,:)

    !------------------------------------------------------------------------------
    ! HDU 12 : Moments of the grain size distribution
    !------------------------------------------------------------------------------
    !  move to next hdu
    call ftmrhd(unit,1,hdutype,status)

    ! Check dimensions
    call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
    if (nfound /= 3) then
       write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
       write(*,*) 'of HDU 7 file. Exiting.'
       stop
    endif
    if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= 4)) then
       write(*,*) "Error : HDU 7 does not have the"
       write(*,*) "right dimensions. Exiting."
       stop
    endif
    npixels=naxes(1)*naxes(2)*naxes(3)

    ! read_image
    call ftgpve(unit,group,firstpix,npixels,nullval,mcfost%grain_size,anynull,status)

    !  Close the file and free the unit number.
    call ftclos(unit, status)
    call ftfiou(unit, status)

    !  Check for any error, and if so print out error messages
    !  Get the text string which describes the error
    if (status > 0) then
       call ftgerr(status,errtext)
       print *,'FITSIO Error Status =',status,': ',errtext

       !  Read and print out all the error messages on the FITSIO stack
       call ftgmsg(errmessage)
       do while (errmessage .ne. ' ')
          print *,errmessage
          call ftgmsg(errmessage)
       end do
    endif

    return

  end function read_mcfost

end module mcfost2ProDiMo

!****************************************

program test

  use mcfost2ProDiMo

  implicit none
  
  type(mcfost_model) :: mcfost
  
  real :: fact_lambda, delta_lambda, wl, delta_wl, B, Sum_J, Sum_B, cst_th, T, Sum
  integer :: i, j, l

  real, parameter :: hp = 6.626e-34
  real, parameter :: c = 3.0e8
  real, parameter :: kb = 1.38e-23
  integer, parameter :: db = selected_real_kind(p=13,r=200)
 
  integer :: iz
  real :: delz, int_dust1, int_dust2, int_dz, rho_gas, rho_dust1, n_dust, dust_to_gas_c, rho_gr, m_dust, rho_dust2
  real, dimension(nz) :: dz

  dust_to_gas_c = 0.01


  cst_th = hp*c/kb 

  mcfost =  read_mcfost()

  write(*,*) "Teff", mcfost%Teff
  write(*,*) "Rstar", mcfost%Rstar

  write(*,*) "disk mass", mcfost%diskmass
  write(*,*) "Rin", mcfost%rin
  write(*,*) "Rout", mcfost%rout

  write(*,*) "amin", mcfost%amin
  write(*,*) "amax", mcfost%amax
  write(*,*) "aexp", mcfost%aexp
  write(*,*) "rho grain", mcfost%rho_grain

  write(*,*) "r_grid", minval(mcfost%r_grid), maxval(mcfost%r_grid)
  write(*,*) "z_grid", minval(mcfost%z_grid), maxval(mcfost%z_grid)
  write(*,*) "wavelengths", minval(mcfost%wavelengths), maxval(mcfost%wavelengths)
  write(*,*) "J_star", minval(mcfost%J_star_input), maxval(mcfost%J_star_input)
  write(*,*) "J", minval(mcfost%J), maxval(mcfost%J)
  write(*,*) "nJ", minval(mcfost%nJ), maxval(mcfost%nJ)
  write(*,*) "J_ISM", minval(mcfost%J_ISM), maxval(mcfost%J_ISM)
  write(*,*) "nJ_ISM", minval(mcfost%nJ_ISM), maxval(mcfost%nJ_ISM)
  write(*,*) "density", minval(mcfost%density), maxval(mcfost%density)
  write(*,*) "grain_size", minval(mcfost%grain_size), maxval(mcfost%grain_size)
  write(*,*) "opacity", minval(mcfost%kappa_ext), maxval(mcfost%kappa_ext)


  ! Verification radiation field vs temperature
!  Sum = 0.0
!  fact_lambda = (mcfost%wavelengths(n_lambda) /  mcfost%wavelengths(1))**(1./(n_lambda-1))
!  delta_lambda = fact_lambda - 1./fact_lambda
!  do i=1, n_rad
!     Sum_J = 0. ;  Sum_B = 0. 
!     
!     do j=1,nz
!        T = mcfost%temperature(i,j)
!
!        do l=1,n_lambda
!           wl = mcfost%wavelengths(l) * 1.e-6
!           delta_wl = wl * delta_lambda
!           B = 2.*hp*c**2 *  1./ ( (exp(cst_th/(T*wl)) -1.) * wl**5 ) 
!        
!           ! mcfost%J is lambda.J_lambda
!           Sum_J = Sum_J + mcfost%kappa_abs(i,j,l) * mcfost%J(i,j,l)/wl * delta_wl 
!           Sum_B = Sum_B + mcfost%kappa_abs(i,j,l) * B * delta_wl 
!        enddo
!     
!     enddo !j
!     !   write(*,*) " "
!     write(*,*) "Verification radiation field vs T", i, Sum_J, Sum_J/Sum_B ;
!     Sum = Sum +  Sum_J/Sum_B
!  enddo
!  write(*,*) Sum / N_rad

     !-----------------------------------------
    ! ***  sum up dust in vertical column  ***
    !-----------------------------------------
    dz(:) = 0.0
    do iz=2,nz
      delz = mcfost%z_grid(1,iz)-mcfost%z_grid(1,iz-1)
      dz(iz) = dz(iz)+0.5*delz
      dz(iz-1) = dz(iz-1)+0.5*delz
    enddo
    

    int_dust1 = 0.0
    int_dust2 = 0.0
    int_dz = 0.0
    do iz=nz,1,-1
      rho_gas   = mcfost%density(1,iz)*0.01/dust_to_gas_c
      rho_dust1 = rho_gas*dust_to_gas_c
      n_dust    = mcfost%grain_size(1,iz,0)
      m_dust    = 4.0/3.0 * pi * mcfost%rho_grain * mcfost%grain_size(1,iz,3)
    
      rho_dust2 = n_dust*m_dust * 1e-18
      
      int_dz    = int_dz+dz(iz)
      int_dust1 = int_dust1 + rho_dust1*dz(iz)
      int_dust2 = int_dust2 + rho_dust2*dz(iz)
   !   write(*,*) rho_dust1,rho_dust2 
    enddo

    write(*,*) "TEST", sum(mcfost%grain_size(1,:,0))
    write(*,*) "TEST", sum(mcfost%density(1,:))


    write(*,*) "TEST_new", sum(mcfost%grain_size(1,:,3) * mcfost%grain_size(1,:,0))

    

    write(*,*) "TEST_density", sum(mcfost%grain_size(1,:,3) * mcfost%grain_size(1,:,0)) * 4./3. * pi * mcfost%rho_grain * 100., sum(mcfost%density(1,:)) *1e18
 
 !   write(*,*) int_dz,mcfost%z_grid(1,nz)-mcfost%z_grid(1,1)
    write(*,*) int_dust1,int_dust2

    write(*,*) dz


end program test
