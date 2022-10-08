module input

  use parametres
  use mcfost_env
  use constantes
  use molecular_emission
  use grains
  use read_DustEM
  use utils, only : in_dir
  use messages
  use wavelengths
  use temperature

  implicit none

  character(len=512) :: Tfile = "./data_th/Temperature.fits.gz"
  character(len=512) :: Tfile_nLTE = "./data_th/Temperature_nLTE.fits.gz"
  character(len=512) :: Tfile_Diff_approx = "./data_th/Temperature_DA.fits.gz"
  character(len=512) :: Tfile_nRE = "./data_th/Temperature_nRE.fits.gz"

  contains


subroutine read_molecules_names(imol)

  integer, intent(in) :: imol

  character(len=80) :: junk
  character(len=512) :: filename, dir
  integer :: ios

  filename = trim(mol(imol)%filename)
  dir = in_dir(filename, mol_dir,  status=ios)
  if (ios /=0) then
     call error("molecule file cannot be found: "//trim(filename))
  else
     filename = trim(dir)//trim(filename) ;
     write(*,*) "Reading "//trim(filename) ;
  endif

  open(unit=1, file=filename, status="old")

  read(1,*) junk
  read(1,'(a)') mol(imol)%name
  close(unit=1)

  return

end subroutine read_molecules_names

!*************************************************************


subroutine readmolecule(imol)
  ! Lit les parametres de la molecule etudiee
  ! remplit Aul, Bul, Blu + les f,   Level_energy, transfreq,
  ! iTransUpper, iTransLower + les donnees de collisions

  integer, intent(in) :: imol
  integer, parameter :: nCollTemp_max = 50

  character(len=515) :: filename, dir
  character(len=80) :: junk
  integer :: i, j, iLow, iUp, iPart, ios
  real :: a, freq, eu, molecularWeight
  real, dimension(nCollTemp_max) :: collrates_tmp, colltemps_tmp

  write(*,*) ""
  write(*,*) "-------------------------------------------------------"

  filename = trim(mol(imol)%filename)
  dir = in_dir(filename, mol_dir,  status=ios)
  if (ios /=0) then
     call error("molecule file cannot be found: "//trim(filename))
  else
     filename = trim(dir)//trim(filename) ;
     write(*,*) "Reading "//trim(filename) ;
  endif

  open(unit=1, file=filename, status="old")

  read(1,*) junk
  read(1,'(a)') mol(imol)%name

  read(1,*) junk
  read(1,*) molecularWeight
  masse_mol = masseH * molecularWeight

  read(1,*) junk
  read(1,*) nLevels

  allocate(Level_energy(nLevels),poids_stat_g(nLevels),j_qnb(nLevels))

  read(1,*) junk
  do i = 1, nLevels
     read(1,*) j, Level_energy(i), poids_stat_g(i) , j_qnb(i)
 !    read(1,*) j, Level_energy(i), poids_stat_g(i) !, buffer
     Level_energy(i) = Level_energy(i) / 8065.541  ! per cm to ev
  enddo

  read(1,*) junk
  read(1,*) nTrans_tot

  allocate(Aul(1:nTrans_tot),fAul(nTrans_tot))
  allocate(Bul(1:nTrans_tot),fBul(nTrans_tot))
  allocate(Blu(1:nTrans_tot),fBlu(nTrans_tot))
  allocate(transfreq(1:nTrans_tot))
  allocate(itransUpper(1:nTrans_tot))
  allocate(itransLower(1:nTrans_tot))
  read(1,*) junk
  do i = 1, nTrans_tot
     read(1,*) j, iUp, iLow, a, freq, eu

     Aul(i) = a ! en s-1
     transfreq(i) = freq*1.d9 ! passage en Hz
     itransUpper(i) = iUp
     itransLower(i) = iLow

     ! Transformation Aul -> Bul
     Bul(i) = a * (c_light**2)/(2.d0*hp*(transfreq(i))**3)
     ! Transformation Bul -> Blu
     Blu(i) = Bul(i) * poids_stat_g(iUp)/poids_stat_g(iLow)
  enddo

  fAul(:) = Aul(:) * hp * transfreq(:)/(4*pi)
  fBul(:) = Bul(:) * hp * transfreq(:)/(4*pi)
  fBlu(:) = Blu(:) * hp * transfreq(:)/(4*pi)

  read(1,*) junk
  read(1,*) nCollPart

  allocate(nCollTrans(1:nCollPart))
  allocate(nCollTemps(1:nCollPart))
  allocate(collTemps(1:nCollPart, 1:nCollTemp_max))
  allocate(collBetween(1:nCollPart))

  do iPart = 1, nCollPart

     read(1,*) junk
     read(1,*) collBetween(iPart)

     read(1,*) junk
     read(1,*) nCollTrans(iPart)

     read(1,*) junk
     read(1,*) nCollTemps(iPart)


     read(1,*) junk
     read(1,*) collTemps_tmp(1:nCollTemps(ipart))
     collTemps(iPart,1:nCollTemps(ipart)) = colltemps_tmp(1:nCollTemps(ipart))

     if (iPart == 1) then
        allocate(collRates(1:nCollPart, 1:nCollTrans(iPart) + 50, 1:nCollTemps(ipart) + 50)) ! TODO : passage par des pointeurs, c'est crade
        allocate(iCollUpper(1:nCollPart, 1:nCollTrans(iPart) + 50))
        allocate(iCollLower(1:nCollPart, 1:nCollTrans(iPart) + 50))
        collRates = 0.d0
     endif

     read(1,*) junk
     do j = 1, nCollTrans(iPart)
        read(1,*) i, iCollUpper(ipart,j),  iCollLower(ipart,j), collRates_tmp(1:nCollTemps(iPart))
        collRates(iPart, j, 1:nCollTemps(iPart)) = collRates_tmp(1:nCollTemps(iPart))
     enddo
  enddo

  close(unit=1)

  write(*,*) ""
  write(*,*) "-------------------------------------------------------"
  write(*,*) trim(mol(imol)%name)//" molecular file read successfully"

  return

end subroutine readmolecule

!*************************************************************

subroutine lect_Temperature()

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, hdutype
  real :: nullval
  integer, dimension(5) :: naxes
  logical :: anynull

  if (lemission_atom) return   !B. Tessore., Temporary

  if (lRE_LTE) then
     status=0
     !  Get an unused Logical Unit Number to use to open the FITS file.
     call ftgiou(unit,status)

     Tfile = trim(root_dir)//"/"//trim(Tfile)
     write(*,*) "Reading temperature file : "//trim(Tfile)

     readwrite=0
     call ftopen(unit,Tfile,readwrite,blocksize,status)
     if (status /= 0) call error("temperature file needed")

     group=1
     firstpix=1
     nullval=-999

     !  determine the size of temperature file
     call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
     if (lVoronoi) then
        if (nfound /= 1) then
           write(*,*) "# fits file,   required"
           write(*,*) nfound, 1
           call error('failed to read the NAXISn keywords of '//trim(Tfile))
        endif
        if ((naxes(1) /= n_cells)) then
           write(*,*) "# fits file,   required"
           write(*,*) naxes(1), n_cells
           call error("Temperature.fits.gz does not have the right dimensions.")
        endif
        npixels=naxes(1)
     else
        if (l3D) then
           if (nfound /= 3) then
              write(*,*) "# fits file,   required"
              write(*,*) nfound, 3
              call error('failed to read the NAXISn keywords of '//trim(Tfile))
           endif
           if ((naxes(1) /= n_rad).or.(naxes(2) /= 2*nz).or.(naxes(3) /= n_az)) then
              write(*,*) "# fits file,   required"
              write(*,*) naxes(1), n_rad
              write(*,*) naxes(2), 2*nz
              write(*,*) naxes(3), n_az
              call error("Temperature.fits.gz does not have the right dimensions.")
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)
        else
           if (nfound /= 2) then
              write(*,*) "# fits file,   required"
              write(*,*) nfound, 2
              call error('failed to read the NAXISn keywords of '//trim(Tfile))
           endif

           if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
              call error("Temperature.fits.gz does not have the right dimensions.")
           endif
           npixels=naxes(1)*naxes(2)
        endif ! l3D
     endif ! lVoronoi

     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Tdust,anynull,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
  endif

  if (lRE_nLTE) then
     status=0
     !  Get an unused Logical Unit Number to use to open the FITS file.
     call ftgiou(unit,status)

     Tfile_nLTE = trim(root_dir)//"/"//trim(Tfile_nLTE)
     write(*,*) "Reading temperature file : "//trim(Tfile_nLTE)

     readwrite=0
     call ftopen(unit,Tfile_nLTE,readwrite,blocksize,status)
     if (status /= 0) call error("temperature file needed")

     group=1
     firstpix=1
     nullval=-999

     ! determine the size of temperature file
     call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
     if (lVoronoi) then
        if (nfound /= 2) then
           write(*,*) "Found", nfound, "axes", naxes
           call error('failed to read the NAXISn keywords of '//trim(Tfile_nLTE))
        endif
        npixels=naxes(1)*naxes(2)
     else
        if (l3D) then
           if (nfound /= 4) then
              write(*,*) "Found", nfound, "axes", naxes
              call error('failed to read the NAXISn keywords of '//trim(Tfile_nLTE)//' file.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           if (nfound /= 3) then
              write(*,*) "Found", nfound, "axes", naxes
              call error('failed to read the NAXISn keywords of '//trim(Tfile_nLTE)//' file.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)
        endif
     endif

     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Tdust_1grain,anynull,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
  endif



  if (lnRE) then
     status=0
     !  Get an unused Logical Unit Number to use to open the FITS file.
     call ftgiou(unit,status)

     Tfile_nRE=trim(root_dir)//"/"//trim(Tfile_nRE)
     write(*,*) "Reading temperature file : "//trim(Tfile_nRE)
     !call init_tab_Temp()

     readwrite=0
     call ftopen(unit,Tfile_nRE,readwrite,blocksize,status)
     if (status /= 0) call error("ERROR : temperature file needed")

     group=1
     firstpix=1
     nullval=-999

     ! HDU 1 : Teq
     call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
     if (lVoronoi) then
        if (nfound /= 2) then
           write(*,*) "Found", nfound, "axes", naxes
           call error('failed to read the NAXISn keywords of '//trim(Tfile_nLTE)//' file.')
        endif
        npixels=naxes(1)*naxes(2)
     else
        if (l3D) then
           if (nfound /= 4) then
              call error('failed to read the NAXISn keywords of Temperature_nRE.fits.gz file HDU 1.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           if (nfound /= 3) then
              call error('failed to read the NAXISn keywords of Temperature_nRE.fits.gz file HDU 1.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)
        endif
     endif
     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Tdust_1grain_nRE,anynull,status)

     ! HDU 2 : is_eq
     call ftmahd(unit,2,hdutype,status)
     call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
     if (lVoronoi) then
        if (nfound /= 2) then
           write(*,*) "Found", nfound, "axes", naxes
           call error('failed to read the NAXISn keywords of '//trim(Tfile_nLTE)//' file.')
        endif
        npixels=naxes(1)*naxes(2)
     else
        if (l3D) then
           if (nfound /= 4) then
              call error('failed to read the NAXISn keywords of Temperature_nRE.fits.gz file HDU 2.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        else
           if (nfound /= 3) then
              call error('failed to read the NAXISn keywords of Temperature_nRE.fits.gz file HDU 2.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)
        endif
     endif
     nbuffer=npixels
     ! read_image
     call ftgpvj(unit,group,firstpix,nbuffer,nullval,l_RE,anynull,status)

     ! HDU 4 : proba temperature
     call ftmahd(unit,4,hdutype,status)

     call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)
     if (lVoronoi) then
        if (nfound /= 3) then
           write(*,*) "Found", nfound, "axes", naxes
           call error('failed to read the NAXISn keywords of '//trim(Tfile_nLTE)//' file.')
        endif
        npixels=naxes(1)*naxes(2)*naxes(3)
     else
        if (l3D) then
           if (nfound /= 5) then
              call error('failed to read the NAXISn keywords of Temperature_nRE.fits.gz file HDU 4.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)*naxes(5)
        else
           if (nfound /= 4) then
              call error('failed to read the NAXISn keywords of Temperature_nRE.fits.gz file HDU 4.')
           endif
           npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
        endif
     endif
     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Proba_Tdust,anynull,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
  endif

  return

end subroutine lect_Temperature

!**********************************************************************

subroutine read_abundance(imol)

  ! C.pinte   1/10/07

  integer, intent(in) :: imol

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels
  real :: nullval
  integer, dimension(4) :: naxes
  logical :: anynull
  character(len=512) abundance_file

  abundance_file = mol(imol)%abundance_file

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  write(*,*) "Reading abundance file : "//trim(abundance_file)

  readwrite=0
  call ftopen(unit,abundance_file,readwrite,blocksize,status)
  if (status /= 0)call error("abundance file needed")

  group=1
  firstpix=1
  nullval=-999

  !  determine the size of temperature file
  call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)

  if (lVoronoi) then
     if (nfound /= 1) then
        write(*,*) "# fits file,   required"
        write(*,*) nfound, 1
        call error('failed to read the NAXISn keywords of '//trim(abundance_file)//' file.')
     endif
     if ((naxes(1) /= n_cells)) then
        write(*,*) "# fits file,   required"
        write(*,*) naxes(1), n_cells
        call error(trim(abundance_file)//" does not have the right dimensions.")
     endif
     npixels=naxes(1)
  else
     if (l3D) then
        if (nfound /= 3) then
           write(*,*) "# fits file,   required"
           write(*,*) nfound, 3
           call error('failed to read the NAXISn keywords of '//trim(abundance_file)//' file.')
        endif
        if ((naxes(1) /= n_rad).or.(naxes(2) /= 2*nz).or.(naxes(3) /= n_az)) then
           write(*,*) "# fits file,   required"
           write(*,*) naxes(1), n_rad
           write(*,*) naxes(2), 2*nz
           write(*,*) naxes(3), n_az
           call error(trim(abundance_file)//" does not have the right dimensions.")
        endif
        npixels=naxes(1)*naxes(2)*naxes(3)
     else
        if (nfound /= 2) then
           write(*,*) "# fits file,   required"
           write(*,*) nfound, 2
           call error('failed to read the NAXISn keywords of '//trim(abundance_file)//' file.')
        endif
        if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
           call error(trim(abundance_file)//" does not have the right dimensions.")
        endif
        npixels=naxes(1)*naxes(2)
     endif
  endif
  nbuffer=npixels

  ! read_image
  call ftgpve(unit,group,firstpix,nbuffer,nullval,tab_abundance,anynull,status)

  call ftclos(unit, status)
  call ftfiou(unit, status)

  return

end subroutine read_abundance

!**********************************************************************

subroutine lect_lambda()
  ! Remplit la variable lambda_filename
  ! et le tableau de lambda pour sed2

  integer :: alloc_status, lambda, status, n_comment, i, ios
  real :: fbuffer

  real :: wl_factor = 1.0005

  character(len=512) :: dir

  lambda_filename = trim(tab_wavelength)

  dir = in_dir(lambda_filename, lambda_dir,  status=ios)
  if (ios /=0) then
     call error("lambda file cannot be found: "//trim(lambda_filename))
  else
     lambda_filename = trim(dir)//trim(lambda_filename) ;
     write(*,*) "Reading "//trim(lambda_filename) ;
  endif

  open(unit=1,file=lambda_filename,status="old",iostat=status)
  if (status /= 0) then
     call warning(trim(lambda_filename)//"' does not exit.", &
          msg2="Looking into working directory ...")

     lambda_filename = trim(tab_wavelength)
     open(unit=1,file=lambda_filename,status="old",iostat=status)
     if (status /= 0) then
        call error(trim(lambda_filename)//" does not exit.")
     else
        write(*,*) "Using local file: "//trim(lambda_filename)
     endif
  endif

  ! On elimine les lignes avec des commentaires
  status = 1
  n_comment = 0
  do while (status /= 0)
     n_comment = n_comment + 1
     read(1,*,iostat=status) fbuffer
  enddo
  n_comment = n_comment - 1

  ! On compte d'abord le nombre de lignes
  status=0
  n_lambda2=1 ! On a deja lu une ligne en cherchant les commentaires
  do while(status==0)
     n_lambda2=n_lambda2+1
     read(1,*,iostat=status)
  enddo
  n_lambda2 = n_lambda2-1

  rewind(1)

  ! On passe les commentaires
  do i=1, n_comment
     read(1,*)
  enddo

  ! Allocations des tab
  allocate(tab_lambda2(n_lambda2), tab_lambda2_inf(n_lambda2),tab_lambda2_sup(n_lambda2), &
       tab_delta_lambda2(n_lambda2), stat=alloc_status)
  if (alloc_status > 0) call error('Allocation error tab_lambda2')

  ! Lecture proprement dite
  do lambda=1, n_lambda2
     read(1,*) tab_lambda2(lambda)!, tab_delta_lambda2(lambda)
     tab_lambda2_inf(lambda) = tab_lambda2(lambda) / wl_factor
     tab_lambda2_sup(lambda) = tab_lambda2(lambda) * wl_factor
     tab_delta_lambda2(lambda) = tab_lambda2_sup(lambda) - tab_lambda2_inf(lambda)
  enddo
  close(unit=1)

  if (lProDiMo) then
     if (abs(tab_lambda2(1) - 0.0912) > 1e-6) then
        call warning("matching step2 wavelength bins to ProDiMo set-up")
        tab_lambda2_inf(1) = 0.0912
        tab_lambda2_sup(1) = tab_lambda2(1) * tab_lambda2(1)/tab_lambda2_inf(1)
        tab_delta_lambda2(1) = tab_lambda2_sup(1) - tab_lambda2_inf(1)
        do lambda=2, n_lambda2
           tab_lambda2_inf(lambda) = tab_lambda2_sup(lambda-1)
           tab_lambda2_sup(lambda) = tab_lambda2(lambda)  * tab_lambda2(lambda)/tab_lambda2_inf(lambda)
           tab_delta_lambda2(lambda) = tab_lambda2_sup(lambda) - tab_lambda2_inf(lambda)
        enddo

        do lambda=1, n_lambda2
           if (tab_delta_lambda2(lambda) <= 0.0) then
              write(*,*) "ERROR in ProDiMo wavelength grid."
              write(*,*) "Exiting"
              write(*,*) lambda, "wl=", tab_lambda2(lambda), "delta=", tab_delta_lambda2(lambda)
              write(*,*) "            min, max =", tab_lambda2_inf(lambda), tab_lambda2_sup(lambda)
              call exit(1)
           endif
        enddo
     endif ! test tab_lambda2 == 0.0912
  endif !  lProDiMo

  do lambda=1, n_lambda2
     tab_delta_lambda2(lambda) = 0.01 * tab_lambda2(lambda)
  enddo

  return

end subroutine lect_lambda

!***************************************************

subroutine read_limb_darkening_file(lambda)

  integer, intent(in) :: lambda

  integer :: lb_n_wl, lb_n_mu
  real, dimension(:,:), allocatable :: Imu, Q_o_Imu ! lb_n_mu, lb_n_wl
  real, dimension(:), allocatable :: lb_wavelength ! lb_n_wl
  real, dimension(:), allocatable :: tmp ! lb_n_mu

  character(len=4) :: buffer4
  character(len=20) :: buffer8
  real :: fbuffer

  real :: frac
  integer :: i, j, alloc_status


  write(*,*) "Reading limb darkening & stellar polarisation file: "//trim(limb_darkening_file)

  open(unit=1,file=limb_darkening_file)
  ! line 1 : skipping header line
  read(1,*) buffer4

  if (buffer4 == "Teff") then
     lb_n_wl = 20 ; lb_n_mu = 159
  else if (buffer4 == "Rmax") then
     lb_n_wl = 10 ; lb_n_mu = 99
  else
     call error("Unknown limb darkening file format")
  endif

  allocate(limb_darkening(lb_n_mu), pola_limb_darkening(lb_n_mu), mu_limb_darkening(lb_n_mu), stat=alloc_status)
  if (alloc_status /= 0) call error("Allocation error in limb_darkening")

  allocate(Imu(lb_n_mu, lb_n_wl), Q_o_Imu(lb_n_mu, lb_n_wl), lb_wavelength(lb_n_wl), tmp(lb_n_mu), stat=alloc_status)
  if (alloc_status /= 0) call error("Allocation error in limb_darkening (local variables)")

  if (buffer4 == "Teff") then
     ! line 2 : reading 20 wavelengths & converting to mum
     read(1,*) lb_wavelength ! the wavelengths are from long to short
     lb_wavelength = lb_wavelength * A_to_mum

     if ( (tab_lambda(lambda) > lb_wavelength(1)).or.(tab_lambda(lambda) < lb_wavelength(lb_n_wl)) ) then
        write(*,*) "ERROR : wavelength is outside limb darkening wavelength range"
        write(*,*) "valod values are between", lb_wavelength(lb_n_wl), "and", lb_wavelength(1), "mum"
        write(*,*) "Exiting"
        call exit(1)
     endif

     ! line 3 : sKipping Flux line
     read(1,*)

     ! line 4 : reading 159 mu values
     read(1,*) mu_limb_darkening

     ! line 5 : skipping the theta values
     read(1,*)

     do i=1, lb_n_wl
        read(1,*) tmp ; Imu(:,i) = tmp(:)
        read(1,*) ! skipping Q
        read(1,*) tmp ; Q_o_Imu(:,i) = tmp(:)
     enddo

     ! Interpolation
     do i= lb_n_wl-1, 1, -1 ! the wavelengths are from long to short
        if (lb_wavelength(i) > tab_lambda(lambda)) then
           frac = (tab_lambda(lambda) - lb_wavelength(i+1)) / (lb_wavelength(i) - lb_wavelength(i+1))
           exit
        endif
     enddo

     limb_darkening(:) = frac * Imu(:,i) + (1.0 - frac) * Imu(:,i+1)
     pola_limb_darkening(:) = -(frac * Q_o_Imu(:,i) + (1.0 - frac) * Q_o_Imu(:,i+1)) ! signe - in file

  else if (buffer4 == "Rmax") then
     do i=1, lb_n_wl
        read(1,*) buffer8, buffer8, fbuffer
        lb_wavelength(i) = fbuffer

        if (i==1) read(1,*) ! header line (only 1st block)
        do j=1, lb_n_mu
           read(1,*) mu_limb_darkening(j), fbuffer, Imu(j,i), Q_o_Imu(j,i)
        enddo ! j
     enddo ! i

     ! Interpolation
     do i= 2, lb_n_wl ! the wavelengths are froms short to long
        if (lb_wavelength(i) > tab_lambda(lambda)) then
           frac = (lb_wavelength(i) - tab_lambda(lambda)) / (lb_wavelength(i) - lb_wavelength(i-1))
           exit
        endif
     enddo

     limb_darkening(:) = frac * Imu(:,i) + (1.0 - frac) * Imu(:,i-1)
     pola_limb_darkening(:) = frac * Q_o_Imu(:,i) + (1.0 - frac) * Q_o_Imu(:,i-1)

  else
     write(*,*) "mcfost dhould never be here"
     write(*,*) "Exiting"
  endif

  ! Normalizing limb_darkening
  limb_darkening(:) =  limb_darkening(:) / limb_darkening(lb_n_mu)

  write(*,*) "Maximum limb darkening = ", limb_darkening(1)
  write(*,*) "Maximum pola limb darkening = ", pola_limb_darkening(1)

  deallocate(Imu, Q_o_Imu, lb_wavelength, tmp)

  return

end subroutine read_limb_darkening_file

!*************************************************************

subroutine read_phase_function(phase_function_file)

  character(len=512), intent(in) :: phase_function_file

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, hdutype, i
  real :: nullval
  integer, dimension(2) :: naxes
  logical :: anynull

  real, dimension(0:nang_scatt,2) :: PF

  allocate(s11_file(0:nang_scatt))

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  write(*,*) "Reading phase function file : "//trim(phase_function_file)

  readwrite=0
  call ftopen(unit,phase_function_file,readwrite,blocksize,status)
  if (status /= 0) call error("phase function file needed")

  group=1
  firstpix=1
  nullval=-999

  !  determine the size of temperature file
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  if ((naxes(1) /= nang_scatt+1)) then
     write(*,*) "AXIS #1"
     write(*,*) "# fits file,   required"
     write(*,*) naxes(1), nang_scatt+1
     call error("phase function file does not have the right dimensions.")
  endif
  if ((naxes(2) /= 2)) then
     write(*,*) "AXIS #2"
     write(*,*) "# fits file,   required"
     write(*,*) naxes(1), nang_scatt+1
     call error("phase function file does not have the right dimensions.")
  endif
  npixels=naxes(1)*naxes(2)

  ! read_image
  nbuffer=npixels
  call ftgpve(unit,group,firstpix,nbuffer,nullval,PF,anynull,status)

  call ftclos(unit, status)
  call ftfiou(unit, status)

  ! Checking angles, will replace by interpolation later
  do i=0,nang_scatt
     if (abs(PF(i,1) - i) > 1e-5) then
        write(*,*) PF(i,1)
        write(*,*) "ERROR: phase function file"
        write(*,*) "Angle #", i, "should be equal to", i
        write(*,*) "Exiting"
        call exit(1)
     endif
  enddo
  s11_file(:) = PF(:,2)

  return

end subroutine read_phase_function


end module input
