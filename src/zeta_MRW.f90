module zeta_MRW

    implicit none

    real, dimension(:, :), allocatable :: zeta_y_data
    real, dimension(:), allocatable :: zeta, y_MRW 
    character(len=512) :: zeta_file = "./zeta_MRW.fits"

 contains

    subroutine read_zeta()

        integer :: status, readwrite, unit, blocksize,nfound,group,firstpix,nbuffer,npixels, hdutype
        real :: nullval
        integer, dimension(2) :: naxes
        logical :: anynull

        zeta_file = trim(".")//"/"//trim(zeta_file)
        write(*,*) "Reading zeta file : "//trim(zeta_file)

        status=0
        !  Get an unused Logical Unit Number to use to open the FITS
        !  file.
        call ftgiou(unit,status)

        readwrite=0
        call ftopen(unit,zeta_file,readwrite,blocksize,status)

        group=1
        firstpix=1
        nullval=-999

        call ftgknj(unit,'NAXIS',1,10,naxes,nfound,status)

        npixels=naxes(1)*naxes(2)

        nbuffer=npixels
        ! read_image
        allocate(zeta_y_data(naxes(1),naxes(2)))
        allocate(zeta(naxes(2)))
        allocate(y_MRW(naxes(2)))
        call ftgpve(unit,group,firstpix,nbuffer,nullval,zeta_y_data,anynull,status)

        zeta(:) = zeta_y_data(2,:)
        y_MRW(:) = zeta_y_data(1,:)

        ! write(*,*) status
        ! write(*,*) zeta

    end subroutine read_zeta

end module zeta_MRW

