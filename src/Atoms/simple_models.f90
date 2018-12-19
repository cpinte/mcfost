MODULE simple_models

 use atmos_type, only : atmos, init_atomic_atmos
 
 ! MCFOST's modules
 use input
 use parametres
 use grid
 use density
 use constantes
 use fits_utils, only : print_error
 
 IMPLICIT NONE

 CONTAINS
 
  SUBROUTINE uniform_law_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple model with density and electron density
   ! constant on the grid
   ! Values width idk are from a Solar model atmosphere by 
   ! Fontenla et al. namely the FAL C model.
  ! ----------------------------------------------------------- !
   double precision, dimension(n_cells) :: nHtot, T, ne
   double precision :: Vkepmax
   !idk = 10
   nHtot =  2.27414200581936d16
   T = 45420d0
   ne = 2.523785d16

   !idk = 75
!    T=7590d0
!    ne = 4.446967d20
!    nHtot = 1.259814d23

    !idk = 81   
!     T=9400d0
!     ne = 3.831726d21
!     nHtot = 1.326625d23

    !idk = 0 
!     T = 100000d0
!     ne = 1.251891d16
!     nHtot = 1.045714d16
    
   !!more or less the same role as init_molecular_disk
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving = .false.
   !atmos%vturb = 0d0
   atmos%vturb = 9.506225d3 !m/s !idk=10
   !atmos%vturb = 1.696164d3 !idk = 75
   !atmos%vturb = 1.806787d3 !idk=81
   !atmos%vturb = 10.680960d3 !idk=0
   atmos%velocity_law = 0 !keplerian = 0
   Vkepmax = 0d0
   atmos%v_char = atmos%v_char + Vkepmax + maxval(atmos%vturb)
   write(*,*) "Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
   if (atmos%moving.and.atmos%velocity_law == 0) then
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
   double precision :: Vkepmax
   
   nHtot = 1d14! * densite_gaz/MAXVAL(densite_gaz) + 1d12
   T = 3000d0!+Tdust
   ne = 1d-2 * maxval(nHtot)! * nHtot/maxval(nHtot)
   
   CALL init_atomic_atmos(n_cells, T, ne, nHtot)
   atmos%moving=.true.
   atmos%velocity_law = 0 !keplerian
   if (atmos%moving) then
    Vkepmax = dsqrt(Ggrav * sum(etoile%M) * Msun_to_kg * Rmin**2 / &
                ((Rmin**2 + minval(z_grid)**2)**1.5 * AU_to_m) )

    atmos%v_char = Vkepmax + maxval(atmos%vturb)
    write(*,*) " >-< Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
    if (atmos%velocity_law == 0) then
      if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
      atmos%Vmap = 0d0
    end if
   end if

  RETURN
  END SUBROUTINE prop_law_model
  
  SUBROUTINE spherically_symmetric_shells_model()
  ! ----------------------------------------------------------- !
   ! Implements a simple spherically symmetric model of shells
   ! in expansion or in contraction.
  ! ----------------------------------------------------------- !
  integer :: izone, i, j, k, icell, n_zones=1 !Only one component firstly
  double precision, dimension(n_cells) :: Vr, Tr,vturbr, nHtotr, ner
  type(disk_zone_type) :: dz ! to define the properties of the model
  double precision :: r, ne0, T0, vt, rcyl, z, v0, rho0, Mdot, r0, beta = 5d-1
  double precision :: alpha, yr_to_sec = 3.154d7, nH0, r1, Rlimit, Vlimit !s/yr
  

  Mdot = 1d-4 * Msun_to_kg / yr_to_sec ! Msun/yr -> kg/s
  nH0 = 1d15 !m^-3 !the inner radius should be opaque.
  alpha = 0.5d0
  rho0 = nH0 * masseH * 1d-3
  ne0 = 1d-2 * nH0
  T0 = 7d3
  vt = 0d0
  v0 = 100d0 !km/s
  r0 = 1.5*Rmin!AU
  !Vlimit = 1d9 !km/s, stop the calculation at this limit
  !Rlimit = r0 * (Vlimit/v0) ** (1./alpha) !depending on the velocity law
  Rlimit = Rmax
  
  Tr = 0d0
  nHtotr = 0d0
  Vr = 0d0
  vturbr = 0d0
  ner = 0d0
  
  main : do izone=1, n_zones
   !dz = disk_zone(izone) ! for now the parameters are hardcoded
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
       !Rin = Rmin if edge = 0
       r = sqrt(rcyl**2 + z**2)
       if (r>= Rmin .and. r<=Rlimit) then
 
        Tr(icell) = T0 * dsqrt(Rmin/r)
        nHtotr(icell) = nH0 * dsqrt(Rmin/r)
        !! if ner = 0d0, it is computed in solvene.f90
        !!ner(icell) = 1d-2 * nHtotr(icell)!ne0 * dsqrt(Rmin/r)
        Vr(icell) = 1d3 * v0 * (r/r0)**alpha !m/s

        
        !Tr(icell) = T0 * sqrt(Rmin/r)
        !Vr(icell) = 1d3 * v0 * (1d0-Rmin/r)**beta !m/s
        !nHtotr(icell) = 1d-3/masseH * Mdot / (4*pi * (r*AU_to_m)**2 * Vr(icell)) !m^-3
        
        
        vturbr(icell) = vt*1d3 !m/s   
        !define an opaque core
        if (r <= r0) then
         Tr(icell) = 1d4
         nHtotr(icell) = 1d23
         Vr(icell) = 0d0
        end if   
!         write(*,*) icell, nHtotr(icell), ner(icell), Tr(icell), Vr(icell)
       end if
      end do !k
     end do !j
    end do !i
  end do main !over n_zones
  
    
  CALL init_atomic_atmos(n_cells, Tr, ner, nHtotr)
  atmos%moving=.true.
  atmos%vturb = vturbr
  if (atmos%moving) then
   if (.not.allocated(atmos%Vmap)) allocate(atmos%Vmap(n_cells))
   atmos%Vmap = Vr !and set keyword lkeplerian=.false. and linfall=.true.
  				  ! otherwise vinfall/vkep = cte = expected.
   atmos%velocity_law = -1
   atmos%v_char = maxval(abs(atmos%Vmap)) + maxval(atmos%vturb)
   write(*,*) "Macroscopic motions typical velocity (km/s):", atmos%v_char/1d3
  end if
  RETURN
  END SUBROUTINE spherically_symmetric_shells_model

 END MODULE simple_models