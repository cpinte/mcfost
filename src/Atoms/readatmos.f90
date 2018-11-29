!
!
! This module reads external model "atmosphere" file from MHD simulations.
!
!
! atmos structure is allocated here for the RT Calculations
!
! Note: The word "atmosphere" is a general designation for accretion disk models,
! or whatever. It only requires the total Hydrogen densities (nHtot), eventually
! electron densities , Temperature, velocity fields, magnetic fields etc. Further
! the conditions of atomic line transfer should apply.
!
!
!
!
!
MODULE READATMOS

 use atmos_type, only : atmos
 use project
 use solvene, only : SolveElectronDensity
 use constant

 !MCFOST's original modules
 use parametres
 use grid
 use density
 use constantes, only : AU_to_Rsun, Rsun


 IMPLICIT NONE

 CONTAINS

 SUBROUTINE readPLUTO()
 
 RETURN
 END SUBROUTINE readPLUTO

 SUBROUTINE readatmos_1D(model1D)
  ! Reads 1D spherically symmetric models in fits.gz format
  ! then the model is mapped to mcfost grid
  character(len=*), intent(in) :: model1D
  integer :: unit, EOF, blocksize, hdutype
  integer :: naxis, Ndep, naxis1, i, j, k, icell
  logical :: anynull
  character(len=256) :: some_comments
  double precision :: logg, Rstar, rcyl, z, r, rmax, rmin
  double precision, allocatable, dimension(:) :: rmod, nHmod, Tmod, Vzmod

  EOF = 0

  CALL ftgiou(unit,EOF)
  !open to the first HDU
  CALL ftopen(unit, TRIM(mcfost_utils)//TRIM(model1D), 0, blocksize, EOF)
  CALL ftgknj(unit, 'NAXIS', 1, 1, naxis, naxis1, EOF)
  CALL ftgkyj(unit, "NAXIS1", Ndep, some_comments, EOF)
  CALL ftgkyd(unit, "RSTAR", Rstar, some_comments, EOF)
  CALL ftgkyd(unit, "g", logg, some_comments, EOF)


  logg = log10(logg) !g in m/s2

  allocate(rmod(Ndep))
  allocate(nHmod(Ndep))
  allocate(Tmod(Ndep))
  allocate(Vzmod(Ndep))

  CALL ftgpvd(unit,1,1,Ndep,-999,rmod,anynull,EOF)
  rmod = rmod * Rstar !Rsun
  CALL FTMAHD(unit,2,hdutype,EOF)
  CALL ftgpvd(unit,1,1,Ndep,-999,nHmod,anynull,EOF) !m^-3
  CALL FTMAHD(unit,3,hdutype,EOF)
  CALL ftgpvd(unit,1,1,Ndep,-999,Tmod,anynull,EOF) !K
  CALL FTMAHD(unit,4,hdutype,EOF)
  CALL ftgpvd(unit,1,1,Ndep,-999,Vzmod,anynull,EOF) !km/s
  Vzmod = Vzmod * 1000d0 !m/s

!   write(*,*) Rstar, logg, Ndep
!   write(*,*) "R:", maxval(rmod), minval(rmod)
!   write(*,*) "T:", maxval(Tmod), minval(Tmod)
!   write(*,*) "nH:", maxval(nHmod), minval(nHmod)
!   write(*,*) "Vz:", maxval(Vzmod), minval(Vzmod)

  CALL ftclos(unit, EOF)
  CALL ftfiou(unit, EOF)

  !map the 1D atmos to the mcfost's grid
  rmax = maxval(rmod) * Rstar * Rsun !m
  rmin = minval(rmod) * Rstar * Rsun
  do i=1, n_rad
   do j=j_start,nz
    do k=1, n_az
     if (j==0) then !pourquoi ici
      icell = cell_map(i,1,k)
      rcyl = r_grid(icell) * AU_TO_RSUN * Rsun
      z = 0d0
     else
      icell = cell_map(i,j,k)
      rcyl = r_grid(icell) * AU_TO_RSUN * Rsun
      z = z_grid(icell) * AU_TO_RSUN * Rsun
     endif
     r = sqrt(rcyl**2 + z**2) !m
     if ((r <= rmax).and.(r>=rmin).and.(j/=0)) then

     end if
    end do !k
   end do !j
  end do !i
  !CALL init_atomic_atmos(n_cells, Tmap, NHtotmap, nemap)
  !atmos%Vmap = 0d0

 RETURN
 END SUBROUTINE readatmos_1D

END MODULE READATMOS
