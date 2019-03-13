MODULE zeeman
 ! see LL04
 use atom_type, only : AtomType, ZeemanType, determinate, getorbital, AtomicLine
 use atmos_type, only : atmos
 use messages

 IMPLICIT NONE
 CONTAINS

 FUNCTION Lande(S, L, J) result(g)
  real(8) :: S, J
  integer :: L
  real(8) :: g

  if (J .eq. 0.0) then
    g = 0.0
  else
    g = 1.5 + (S*(S + 1.0) - L*(L + 1)) / (2.0*J*(J + 1.0))
  end if

  RETURN
  END FUNCTION Lande

 SUBROUTINE Lande_eff(atom, kr)
   type (AtomType), intent(inout) :: atom
   integer :: kr
   integer :: i, j
   real(8) g, gi, gj

   i = atom%lines(kr)%i
   j = atom%lines(kr)%j

   !write(*,*) i, j

   gi = Lande(atom%qS(i), atom%Lorbit(i), atom%qJ(i));
   gj = Lande(atom%qS(j), atom%Lorbit(j), atom%qJ(j));
   !write(*,*) "landes", gi, gj

   g = 0.5*(gj + gi) + 0.25*(gj - gi) * &
        (atom%qJ(j)*(atom%qJ(j) + 1.0) - &
        atom%qJ(i)*(atom%qJ(i) + 1.0))

   atom%lines(kr)%g_lande_eff = g
   atom%lines(kr)%glande_i = gi
   atom%lines(kr)%glande_j = gj

 RETURN
 END SUBROUTINE Lande_eff

 SUBROUTINE ZeemanStrength()


 RETURN
 END SUBROUTINE ZeemanStrength

 SUBROUTINE ZeemanMultiplet(line) !Called only once per line
  type(AtomicLine), intent(inout) :: line
  type(ZeemanType), pointer :: zm
  
  !init the Zeeman component for this line at this position

  if (atmos%B_SOLUTION == "WEAK_FIELD") then
   if (line%g_Lande_eff <= -99) then
    CALL Warning("Landé factor not parsed for this line")
    RETURN
   end if
   allocate(line%zm%q(3), line%zm%strength(3), line%zm%shift(3))
   zm => line%zm 
   
   zm%Ncomponent = 3
   zm%q = (/-1, 0, 1/)
   zm%strength = 1d0
   zm%shift = zm%q * line%g_Lande_eff
  else
   CALL Error("Full Zeeman components not implemented yet!")   
  end if

 RETURN
 END SUBROUTINE ZeemanMultiplet

END MODULE Zeeman
