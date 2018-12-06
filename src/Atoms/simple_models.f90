MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos
 use input
 use parametres
 use grid
 use density
 use constantes
 
 
 IMPLICIT NONE

 CONTAINS
 
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid
   ! Values are from a Solar model atmosphere by Fontenla et al.
   ! Namely the FAL C model.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne

   nHtot =  1d12
   T = 3000d0
   ne = 1d-2 * nHtot

!    !idk = 10
!    nHtot =  2.27414200581936d16
!    T = 45420d0
!    ne = 2.523785d16

!    !idk = 75
!    T=7590d0
!    ne = 4.446967d20
!    nHtot = 1.259814d23

!    !idk = 81   
!     T=9400d0
!     ne = 3.831726d21
!     nHtot = 1.326625d23

!    !idk = 0 
!     T = 100000d0
!     ne = 1.251891d16
!     nHtot = 1.045714d16
    
   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving = .true.
   atmos%vturb = 0d0
   !atmos%vturb = 9.506225d3 !m/s !idk=10
   !atmos%vturb = 1.696164d3 !idk = 75
   !atmos%vturb = 1.806787d3 !idk=81
   !atmos%vturb = 10.680960d3 !idk=0
   atmos%velocity_law = 0 !keplerian
   if (atmos%moving.and.atmos%velocity_law /= 0) then
     if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
     atmos%Vmap = 0d0
   end if

  RETURN
  END SUBROUTINE uniform_law_model
  
  SUBROUTINE prop_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! proportional to dust density and temperature.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne

   nHtot = 1d15 * densite_gaz/MAXVAL(densite_gaz) + 1d12
   !nHtot =  1d0 *  1d3 * densite_gaz * masse_mol_gaz / m3_to_cm3 / masseH + 1d12
   T = Tdust * 1d0 + 3000d0!10d0, depends on the stellar flux
   ne = 1d-2 * maxval(nHtot) * nHtot/maxval(nHtot)!1d-2 * nHtot
   
!   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving=.true.
   atmos%velocity_law = 0 !keplerian
   if (atmos%moving.and.atmos%velocity_law /= 0) then
     if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
     atmos%Vmap = 0d0
   end if
   
  RETURN
  END SUBROUTINE prop_law_model
  
  SUBROUTINE spherically_symmetric_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple spherically symmetric model of envelop-
   ! -e in expansion or in contraction.
   ! the gradient of velocity: dlog(Vr)/dlog(r) = alpha, drives
   ! the shape of the emerging flux.
   ! 
   ! alpha > 1: 
   ! See: Bertout & Magnan 1997
  ! ----------------------------------------------------------- !
  integer :: izone, i, j, k, icell, n_zones=1 !Only one component firstly
  double precision, dimension(n_cells) :: Vr, Tr,vturbr, nHtotr, ner
  type(disk_zone_type) :: dz ! to define the properties of the model
  double precision :: r, ne0, T0, vt, rcyl, z, v0, vph, rho0, Mdot, r0
  double precision :: alpha, yr_to_sec = 3.154d7, nH0 !s/yr
  
  nH0 = 1d12 !m^-3 !0d0
  if (nH0 > 0d0) then
   rho0 = nH0 * masseH / 1d3
  else
    rho0 = 7.107d-6!kg/m3
  end if

  ne0 = 1d-2 * 1d3/masseH * rho0
  T0 = 2000d0
  vt = 2d0
  v0 = 1000d0 !km/s
  r0 = 5d0 !AU
  alpha = 2d0

  Mdot = 1d-6 * Msun_to_kg / yr_to_sec ! Msun/yr -> kg/s
  vph = Mdot / (4*PI * (r0*AU_to_m)**2 * rho0) !m/s ! at the top of the stellar photosphere
  write(*,*) " ------ Initial parameters --------"
  write(*,*) Mdot*yr_to_sec / Msun_to_kg, v0/1d3, rho0 * 1d3/masseH, T0, ne0, alpha
  
  Tr = 0d0
  nHtotr = 0d0
  Vr = 0d0
  vturbr = 0d0
  ner = 0d0
  
  write(*,*) " ------ parameters --------"
  do izone=1, n_zones
   !dz = disk_zone(izon) ! for now the parameters are hardcoded
    do i=1, n_rad
     do j=j_start,nz
      do k=1, n_az
       if (j==0) then !midplane
        icell = cell_map(i,1,k)
        rcyl = r_grid(icell) !AU
        z = 0.0_dp
       else
        icell = cell_map(i,j,k)
        rcyl = r_grid(icell)
        z = z_grid(icell)/z_scaling_env
       end if
       r = sqrt(rcyl**2 + z**2)
       !!if ((r <= 300.).and.(r>=r0)) then
        Tr(icell) = T0 * sqrt(r0/r)
        Vr(icell) = 1d3 * (vph + (v0-vph) * (1.-r0/r)**alpha) !m/s
        !Vr(icell) = 1d3 * v0 * (r/r0)**alpha !m/s
        !nHtotr(icell) = 1d3/masseH * Mdot / (4*pi * (r*AU_to_m)**2 * Vr(icell)) !m^-3
        nHtotr(icell) = 1d3/masseH *rho0 * sqrt(r0/r)
        ner(icell) = ne0 * sqrt(r0/r)!nHtotr(icell) !m^-3
        vturbr(icell) = vt
        write(*,*) r/r0, Tr(icell), nHtotr(icell), ner(icell), Vr(icell)/1d3
       !!end if
      end do !k
     end do !j
    end do !i
  end do !over n_zones
  write(*,*) " ------ ------- --------"


  CALL init_atomic_atmos(n_cells, Tr, ner, nHtotr)
  atmos%moving=.true.
  atmos%vturb = vt*1d3
  if (atmos%moving .and. .not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
  atmos%Vmap = 0d0
  atmos%Vmap = Vr !and set keyword lkeplerian=.false. and linfall=.true.
  				  ! otherwise vinfall/vkep = cte = expected.
  atmos%velocity_law = -1
  RETURN
  END SUBROUTINE spherically_symmetric_model
 

 END MODULE simple_models