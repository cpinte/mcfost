module input

  use parametres
  use constantes
  use molecular_emission
  use em_th
  use grains
  use disk
  use utils, only : in_dir

  implicit none

  contains

subroutine read_opacity_file(pop)

  integer, intent(in) :: pop
  character(len=512) :: filename, dir

  integer :: ios, l
  type(dust_pop_type), pointer :: dp

  dp => dust_pop(pop)

  ! We look for the directory and check if the file exists in  get_opacity_file_dim
  ! We do not need to redo it here
  filename = trim(dp%indices(1))
  write(*,*) "Reading opacity file "//trim(filename) ;

  ! load optical data from file
  if (dp%is_Misselt_opacity_file) then
     call misselt_load(filename, pop)
  else
     dp%component_rho1g(1) = 2.5
     dp%rho1g_avg = 2.5
     call draine_load(filename, pop)
     !call draine_load(filename, op_file_n_lambda, op_file_na, 10, 1, &
    !      tmp_lambda, tmp_r_grain,  tmp_Q_ext, tmp_Q_abs, tmp_Q_sca, tmp_g, 4)
  endif

  ! abs car le fichier de lambda pour les PAHs est a l'envers
  op_file_delta_lambda(1,pop) = abs(op_file_lambda(2,pop) - op_file_lambda(1,pop))
  if (n_lambda > 1) then
     op_file_delta_lambda(n_lambda,pop) = abs(op_file_lambda(n_lambda,pop) - op_file_lambda(n_lambda-1,pop))
  endif
  do l=2,op_file_n_lambda(pop)-1
     op_file_delta_lambda(l,pop) = 0.5* abs(op_file_lambda(l+1,pop) - op_file_lambda(l-1,pop))
  enddo

  return

end subroutine read_opacity_file

!*************************************************************

subroutine get_opacity_file_dim(pop)

  integer, intent(in) :: pop
  character(len=512) :: filename, dir
  type(dust_pop_type), pointer :: dp

  integer :: ios

  dp => dust_pop(pop)

  filename = trim(dp%indices(1))
  dir = in_dir(filename, dust_dir,  status=ios)
  if (ios /=0) then
     write(*,*) "ERROR: dust file cannot be found:",trim(filename)
     write(*,*) "Exiting"
     stop
  endif
  filename = trim(dir)//trim(filename) ;

  ! Updating indice filename with directory
  dp%indices(1) = filename ;

  if (dp%is_Misselt_opacity_file) then
     call get_misselt_dim(filename, pop)
  else
     dp%component_rho1g(1) = 2.5
     dp%rho1g_avg = 2.5
     call get_draine_dim(filename, pop)
  endif

  return

end subroutine get_opacity_file_dim

!*************************************************************

subroutine alloc_mem_opacity_file()

  integer :: alloc_status, n_lambda, na, nT, pop

  na = maxval(op_file_na(:))
  n_lambda = maxval(op_file_n_lambda(:))

  allocate(op_file_Qext(n_lambda,na,n_pop), op_file_Qsca(n_lambda,na,n_pop), &
       op_file_g(n_lambda,na,n_pop), op_file_log_r_grain(na,n_pop), &
       op_file_lambda(n_lambda,n_pop), op_file_delta_lambda(n_lambda,n_pop), stat=alloc_status)
  if (alloc_status > 0) then
     write(*,*) 'Allocation error opacity file'
     stop
  endif

  if (lnRE) then
     do pop = 1, n_pop
        if (lread_Misselt.and.dust_pop(pop)%is_opacity_file) call get_file_specific_heat_dim(pop)
     enddo

     nT = maxval(file_sh_nT(:))
     allocate(file_sh_T(nT,n_pop), file_sh(nT,n_pop))
     if (alloc_status > 0) then
        write(*,*) 'Allocation error specific heat capacity file'
        stop
     endif
  endif

  return

end subroutine alloc_mem_opacity_file

!*************************************************************

subroutine free_mem_opacity_file()
  ! Pas utilise pour le moment

  deallocate(op_file_Qext, op_file_Qsca, op_file_g, op_file_log_r_grain, &
       op_file_lambda, op_file_delta_lambda)

  if (allocated(file_sh_T)) deallocate(file_sh_T, file_sh)

  return

end subroutine free_mem_opacity_file

!*************************************************************

subroutine draine_load(file,pop)

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7, nh2 = 1, ns=1

  real :: fbuffer, a
  integer :: i,j,l

  open(unit=1,file=file)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) op_file_na(pop)
  read(1,*) op_file_n_lambda(pop)

  do l=1, nh2
     read(1,*)
  enddo

  ! lecture
  do i=1, op_file_na(pop)
     read(1,*) a ; op_file_log_r_grain(i,pop) = log(a)
     read(1,*)
     do j=1, op_file_n_lambda(pop)
        read(1,*) op_file_lambda(j,pop), op_file_Qext(j,i,pop), fbuffer, op_file_Qsca(j,i,pop), op_file_g(j,i,pop)
     enddo
     if (i < op_file_na(pop)) then
        ! skip ns ligne
        do l=1,ns
           read(1,*)
        enddo
     endif
  enddo
  close(unit=1)

  return

end subroutine draine_load

!*************************************************************

subroutine get_draine_dim(file,pop)

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7
  integer :: l

  open(unit=1,file=file)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) op_file_na(pop)
  read(1,*) op_file_n_lambda(pop)

  close(unit=1)

  return

end subroutine get_draine_dim

!*************************************************************

subroutine misselt_load(file, pop)
  ! Allocate and fills the arrays op_file_XXXX
  ! C. Pinte
  ! 27/06/2014

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7, nh2 = 2, ns=1

  character(len=5) :: sbuffer

  real :: fbuffer, fbuffer2, a, material_density
  integer :: i,j,l, ibuffer

  open(unit=1,file=file)
  read(1,*)  sh_file(pop)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) sbuffer, op_file_na(pop)
  read(1,*) sbuffer, op_file_n_lambda(pop)
  read(1,*) ibuffer, material_density
  dust_pop(pop)%component_rho1g(1) = material_density
  dust_pop(pop)%rho1g_avg = material_density

  do l=1, nh2
     read(1,*)
  enddo

  ! lecture
  do i=1,op_file_na(pop)
     read(1,*) a ; op_file_log_r_grain(i,pop) = log(a)
     read(1,*)
     do j=1,op_file_n_lambda(pop)
        read(1,*) fbuffer, op_file_lambda(j,pop), fbuffer2, op_file_Qsca(j,i,pop), op_file_Qext(j,i,pop), op_file_g(j,i,pop)
     enddo
     if (i < op_file_na(pop)) then
        ! skip ns ligne
        do l=1,ns
           read(1,*)
        enddo
     endif
  enddo
  close(unit=1)

  return

end subroutine misselt_load

!*************************************************************

subroutine get_misselt_dim(file, pop)

  character(len=512), intent(in) :: file
  integer, intent(in) :: pop

  ! nh : n header ; ns : n_skip
  integer, parameter :: nh1 = 7

  character(len=5) :: sbuffer
  real :: material_density
  integer :: l, ibuffer

  open(unit=1,file=file)
  read(1,*)  sh_file(pop)
  ! header
  do l=1, nh1
     read(1,*)
  enddo
  read(1,*) sbuffer, op_file_na(pop)
  read(1,*) sbuffer, op_file_n_lambda(pop)
  read(1,*) ibuffer, material_density
  dust_pop(pop)%component_rho1g(1) = material_density
  dust_pop(pop)%rho1g_avg = material_density

  close(unit=1)

  return

end subroutine get_misselt_dim

!*************************************************************

subroutine read_file_specific_heat(pop)

  integer, intent(in) :: pop

  integer :: ios, status, n_comment, k
  real :: fbuffer

  open(unit=1, file=sh_file(pop), status='old', iostat=ios)
  if (ios/=0) then
     write(*,*) "ERROR : cannot open SH file "//trim(sh_file(pop))
     write(*,*) "Exiting"
     stop
  endif
  write(*,*) "Reading specific heat capacity: "//trim(sh_file(pop))

  ! On elimine les lignes avec des commentaires
  status = 1
  n_comment = 0
  do while (status /= 0)
     n_comment = n_comment + 1
     read(1,*,iostat=status) fbuffer
  enddo
  n_comment = n_comment - 1

  ! Lecture proprement dite
  rewind(1)
  ! On passe les commentaires
  do k=1, n_comment
     read(1,*)
  enddo

  ! Lecture specific heat capacity
  do k=1,file_sh_nT(pop)
     read(1,*) file_sh_T(k,pop), fbuffer, file_sh(k,pop)
  enddo

  return

end subroutine read_file_specific_heat

!******************************************************

subroutine get_file_specific_heat_dim(pop)

  integer, intent(in) :: pop

  integer :: ios, status, n_comment, n_Temperatures, k
  real :: fbuffer

  open(unit=1, file=sh_file(pop), status='old', iostat=ios)
  if (ios/=0) then
     write(*,*) "ERROR : cannot open SH file "//trim(sh_file(pop))
     write(*,*) "Exiting"
     stop
  endif
  write(*,*) "Reading specific heat capacity: "//trim(sh_file(pop))

  ! On elimine les lignes avec des commentaires
  status = 1
  n_comment = 0
  do while (status /= 0)
     n_comment = n_comment + 1
     read(1,*,iostat=status) fbuffer
  enddo
  n_comment = n_comment - 1

  ! On compte les lignes avec des donnees
  status=0
  n_Temperatures=1 ! On a deja lu une ligne en cherchant les commentaires
  do while(status==0)
     n_Temperatures=n_Temperatures+1
     read(1,*,iostat=status)
  enddo
  n_Temperatures = n_Temperatures - 1

  file_sh_nT(pop) = n_Temperatures

  return

end subroutine get_file_specific_heat_dim

!*************************************************************

subroutine read_molecules_names(imol)

  integer, intent(in) :: imol

  character(len=80) :: junk
  character(len=512) :: filename, dir
  integer :: ios

  filename = trim(mol(imol)%filename)
  dir = in_dir(filename, mol_dir,  status=ios)
  if (ios /=0) then
     write(*,*) "ERROR: molecule file cannot be found:",trim(filename)
     write(*,*) "Exiting"
     stop
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
  real :: a, freq, eu
  real, dimension(nCollTemp_max) :: collrates_tmp, colltemps_tmp

  character(len=10) :: buffer

  filename = trim(mol(imol)%filename)
  dir = in_dir(filename, mol_dir,  status=ios)
  if (ios /=0) then
     write(*,*) "ERROR: molecule file cannot be found:",trim(filename)
     write(*,*) "Exiting"
     stop
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
  !   read(1,*) j, Level_energy(i), poids_stat_g(i) , j_qnb(i)
     read(1,*) j, Level_energy(i), poids_stat_g(i) !, buffer
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

  write(*,*) trim(mol(imol)%name)//" molecular file read successfully"

  return

end subroutine readmolecule

!*************************************************************

subroutine lect_Temperature()

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j, hdunum, hdutype
  real :: nullval
  integer, dimension(4) :: naxes
  logical :: anynull

  real, dimension(n_rad,nz,n_az) :: tmp_temp

  if (lRE_LTE) then
     status=0
     !  Get an unused Logical Unit Number to use to open the FITS file.
     call ftgiou(unit,status)

     write(*,*) "Reading temperature file : "//trim(Tfile)

     readwrite=0
     call ftopen(unit,Tfile,readwrite,blocksize,status)
     if (status /= 0) then ! le fichier temperature n'existe pas
        write(*,*) "ERROR : temperature file needed"
        stop
     endif

     group=1
     firstpix=1
     nullval=-999

     !  determine the size of temperature file
     if (l3D) then
        call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
        if (nfound /= 3) then
           write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
           write(*,*) 'of Temperature.fits.gz file. Exiting.'
           stop
        endif
        if ((naxes(1) /= n_rad).or.(naxes(2) /= nz).or.(naxes(3) /= n_az)) then
           write(*,*) "Error : Temperature.fits.gz does not have the"
           write(*,*) "right dimensions. Exiting."
           write(*,*) naxes(1), n_rad
           write(*,*) naxes(2), nz
           write(*,*) naxes(3), n_az

           stop
        endif
        npixels=naxes(1)*naxes(2)*naxes(3)
     else
        call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
        if (nfound /= 2) then
           write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
           write(*,*) 'of Temperature.fits.gz file. Exiting.'
           stop
        endif


        if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
           write(*,*) "Error : Temperature.fits.gz does not have the"
           write(*,*) "right dimensions. Exiting."
           stop
        endif
        npixels=naxes(1)*naxes(2)
     endif



     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,tmp_temp,anynull,status)

     if (l3D) then
        do j=1, nz
           temperature(:,-j,:) = tmp_temp(:,j,:)
           temperature(:,j,:) = tmp_temp(:,j,:)
        enddo
        Temperature(:,0,:) = T_min
     else
        Temperature(:,:,1) = tmp_temp(:,:,1)
     endif

       call ftclos(unit, status)
       call ftfiou(unit, status)
  endif

  if (lRE_nLTE) then
     status=0
     !  Get an unused Logical Unit Number to use to open the FITS file.
     call ftgiou(unit,status)

     write(*,*) "Reading temperature file : "//trim(Tfile)

     readwrite=0
     call ftopen(unit,Tfile,readwrite,blocksize,status)
     if (status /= 0) then ! le fichier temperature n'existe pas
        write(*,*) "ERROR : temperature file needed"
        stop
     endif

     group=1
     firstpix=1
     nullval=-999

     !  determine the size of temperature file
     call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
     if (nfound /= 3) then
        write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
        write(*,*) 'of Temperature.fits.gz file. Exiting.'
        stop
     endif
     npixels=naxes(1)*naxes(2)*naxes(3)
     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Temperature_1grain,anynull,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
  endif



  if (lnRE) then
     status=0
     !  Get an unused Logical Unit Number to use to open the FITS file.
     call ftgiou(unit,status)


     write(*,*) "Reading temperature file : "//trim(Tfile_nRE)
     !call init_tab_Temp()

     readwrite=0
     call ftopen(unit,Tfile_nRE,readwrite,blocksize,status)
     if (status /= 0) then ! le fichier temperature n'existe pas
        write(*,*) "ERROR : temperature file needed"
        stop
     endif

     group=1
     firstpix=1
     nullval=-999


!     call ftcrhd(unit, hdunum,status)
!     write(*,*) hdunum
!     read(*,*)


     !  determine the size of temperature file
     call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
     if (nfound /= 3) then
        write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
        write(*,*) 'of Temperature_nRE.fits.gz file HDU 1. Exiting.'
        stop
     endif
     npixels=naxes(1)*naxes(2)*naxes(3)
     nbuffer=npixels
     ! read_image
     call ftgpvj(unit,group,firstpix,nbuffer,nullval,l_RE,anynull,status)

     ! 2eme HDU
     call ftmahd(unit,2,hdutype,status)

     call ftgknj(unit,'NAXIS',1,4,naxes,nfound,status)
     if (nfound /= 4) then
        write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
        write(*,*) 'of Temperature_nRE.fits.gz file HDU 2. Exiting.'
        stop
     endif
     npixels=naxes(1)*naxes(2)*naxes(3)*naxes(4)
     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Proba_Temperature,anynull,status)

     ! 3eme HDU
     call ftmahd(unit,3,hdutype,status)

     call ftgknj(unit,'NAXIS',1,3,naxes,nfound,status)
     if (nfound /= 3) then
        write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
        write(*,*) 'of Temperature_nRE.fits.gz file HDU 3. Exiting.'
        stop
     endif
     npixels=naxes(1)*naxes(2)*naxes(3)
     nbuffer=npixels
     ! read_image
     call ftgpve(unit,group,firstpix,nbuffer,nullval,Temperature_1grain_nRE,anynull,status)

     call ftclos(unit, status)
     call ftfiou(unit, status)
  endif

  return

end subroutine lect_Temperature

!**********************************************************************

subroutine read_abundance(imol)

  ! C.pinte   1/10/07

  integer, intent(in) :: imol

  integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels,j, hdunum, hdutype
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
  if (status /= 0) then ! le fichier temperature n'existe pas
     write(*,*) "ERROR : abundance file needed"
     stop
  endif

  group=1
  firstpix=1
  nullval=-999

  !  determine the size of temperature file
  call ftgknj(unit,'NAXIS',1,2,naxes,nfound,status)
  if (nfound /= 2) then
     write(*,*) 'READ_IMAGE failed to read the NAXISn keywords'
     write(*,*) 'of '//trim(abundance_file)//' file. Exiting.'
     stop
  endif

  if ((naxes(1) /= n_rad).or.(naxes(2) /= nz)) then
     write(*,*) "Error : "//trim(abundance_file)//" does not have the"
     write(*,*) "right dimensions. Exiting."
     stop
  endif

  npixels=naxes(1)*naxes(2)
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

  real, parameter :: wl_factor = 1.025

  character(len=512) :: dir

  lambda_filename = trim(tab_wavelength)

  dir = in_dir(lambda_filename, lambda_dir,  status=ios)
  if (ios /=0) then
     write(*,*) "ERROR: lambda file cannot be found:",trim(lambda_filename)
     write(*,*) "Exiting"
     stop
  else
     lambda_filename = trim(dir)//trim(lambda_filename) ;
     write(*,*) "Reading "//trim(lambda_filename) ;
  endif

  open(unit=1,file=lambda_filename,status="old",iostat=status)
  if (status /= 0) then
     write(*,*) "WARNING: '"//trim(lambda_filename)//"' does not exit."
     write(*,*) "Looking into working directory ..."

     lambda_filename = trim(tab_wavelength)
     open(unit=1,file=lambda_filename,status="old",iostat=status)
     if (status /= 0) then
        write(*,*) "ERROR: '"//trim(lambda_filename)//"' does not exit."
        write(*,*) "Exiting"
        stop
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
  if (alloc_status > 0) then
     write(*,*) 'Allocation error tab_lambda2'
     stop
  endif

  ! Lecture proprement dite
  do lambda=1, n_lambda2
     read(1,*) tab_lambda2(lambda)!, tab_delta_lambda2(lambda)
     tab_lambda2_inf(lambda) = tab_lambda2(lambda) / wl_factor
     tab_lambda2_sup(lambda) = tab_lambda2(lambda) * wl_factor
     tab_delta_lambda2(lambda) = tab_lambda2_sup(lambda) - tab_lambda2_inf(lambda)
  enddo
  close(unit=1)


  if (lProDiMo) then
     write(*,*) "WARNING:  matching step2 wavelength bins to ProDiMo set-up"
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
           stop
        endif
     enddo
  endif ! lProDiMo

  return

end subroutine lect_lambda

!**********************************************************************

subroutine read_struct_fits_file()
  ! Need to read the FITS file once here to read off the header keywords necessary to define the grid.
  ! Minimal set of keywords: AMIN, AMAX [in micron], RIN, ROUT, ZMAX, RREF [in AU] and BETA, plus the NAXIS suite.
  ! Grid structure has to follow H(r) = H0 x (r/r0)^beta, and ZMAX is the total height at radius RREF.
  ! Format of input FITS file: n_rad x nz x n_grains [x n_species]

  integer :: unit, status, blocksize, readwrite, nfound_dim
  integer, dimension(4) :: naxes
  character(len=80) :: comment

  grid_type=1

  status=0
  !  Get an unused Logical Unit Number to use to open the FITS file.
  call ftgiou(unit,status)

  readwrite=0
  call ftopen(unit,struct_fits_file,readwrite,blocksize,status)
  if (status /= 0)then
     write(*,*) 'ERROR opening fits file '//trim(struct_fits_file)
     stop
  endif

  !  determine the size of the input file (Not more than 4 dimensions!)
  call ftgknj(unit,'NAXIS',1,4,naxes,nfound_dim,status)
  !  check that it found at least 3 NAXISn keywords
  if (nfound_dim < 3) then
     write(*,*) 'ERROR! File '//trim(struct_fits_file)//' does has not enough dimensions (3 are required)'
     stop
  endif

  if (nfound_dim == 4) then
     struct_file_nspecies=naxes(4)
  else
     struct_file_nspecies=1
  endif
  write(*,*) "Input structure file contains ",struct_file_nspecies," dust population(s)"

  ! determine spatial extent of grid (rin, rout, zmax, beta, rref)
  status=0
  call ftgkye(unit,'RIN',struct_file_rin,comment,status)
  if (status /= 0) then
     write(*,*) 'ERROR! Input file does not contain RIN keyword'
     stop
  else
     status=0
     call ftgkye(unit,'ROUT',struct_file_rout,comment,status)
     if (status /= 0) then
        write(*,*) 'ERROR! Input file does not contain ROUT keyword'
        stop
     else
        status=0
        call ftgkye(unit,'ZMAX',struct_file_zmax,comment,status)
        if (status /= 0) then
           write(*,*) 'ERROR! Input file does not contain ZMAX keyword'
           stop
        else
           status=0
           call ftgkye(unit,'BETA',struct_file_beta,comment,status)
           if (status /= 0) then
              write(*,*) 'ERROR! Input file does not contain BETA keyword'
              stop
           else
              status=0
              call ftgkye(unit,'RREF',struct_file_rref,comment,status)
              if (status /= 0) then
                 write(*,*) 'ERROR! Input file does not contain RREF keyword'
                 stop
              endif
           endif
        endif
     endif
  endif

  ! determine extent of grain size distribution
  status=0
  call ftgkye(unit,'AMIN',struct_file_amin,comment,status)
  if (status /= 0) then
     write(*,*) 'ERROR! Input file does not contain AMIN keyword'
     stop
  else
     status=0
     call ftgkye(unit,'AMAX',struct_file_amax,comment,status)
     if (status /= 0) then
        write(*,*) 'ERROR! Input file does not contain AMAX keyword'
        stop
     endif
  endif

  ! closing input file and freeing unit number
  call ftclos(unit, status)
  call ftfiou(unit, status)

  n_rad=naxes(1) ; nz=naxes(2) ; n_az=1 ; n_rad_in=1
  struct_file_n_grains = naxes(3)

  return

end subroutine read_struct_fits_file

!***************************************************

subroutine init_tab_Temp()

  real(kind=db) :: delta_T
  integer :: t

  tab_Temp=0.0
  ! Echantillonage temperature
  !delta_T=(T_max)**(1.0/(n_T-1))
  delta_T=exp((1.0_db/(real(n_T,kind=db)))*log(T_max/T_min))
  tab_Temp(1)=T_min*sqrt(delta_T)
   do t=2,n_T
     tab_Temp(t)=delta_T*tab_Temp(t-1)
  enddo

end subroutine init_tab_Temp

end module input
