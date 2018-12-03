! This module contains functions to project the velocity
! and the magnetic fields according the parameters of the cell
! in which a (local) line profile is computed.
! Use MCFOST routines.


MODULE project

 use atmos_type, only : atmos


 IMPLICIT NONE

 CONTAINS

 FUNCTION vproject(icell, x, y, z, u, v, w, V_x, V_y, V_z, lvoronoi, lkeplerian) result(v_proj)

  logical, intent(in) :: lvoronoi, lkeplerian
  integer, intent(in) :: icell
  double precision , intent(in) :: x,y,z,u,v,w, V_x, V_y, V_z
  double precision :: vitesse, vx, vy, vz, norme, r
  logical :: ldensity_file
  double precision :: v_proj
  ! if not from file, V_x=Vfield and V_y=V_z=0
  
  if (V_y + V_z > 0d0) ldensity_file = .true.

  if (lVoronoi) then
!      vx = Voronoi(icell)%vxyz(1)
!      vy = Voronoi(icell)%vxyz(2)
!      vz = Voronoi(icell)%vxyz(3)
! 
!      v_proj = vx * u + vy * v + vz * w
  else
     if (ldensity_file) then
        vx = V_x ; vy = V_y; vz = V_z
        v_proj = vx * u + vy * v + vz * w
     else ! Using Analytical velocity field
        vitesse = V_x

        if (lkeplerian) then
           r = sqrt(x*x+y*y)
           if (r > 0d0) then
              norme = 1.d0/r
              vx = -y * norme * vitesse
              vy = x * norme * vitesse
              vz = 0.
              v_proj = vx * u + vy * v
           else
              v_proj = 0d0
           endif
        endif
     endif ! ldensity_file
  endif
 RETURN
 END FUNCTION

 FUNCTION Bproject(icell) result(y)
  integer :: icell
  double precision :: y

  y = 0.

 END FUNCTION

END MODULE project
