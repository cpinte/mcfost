MODULE zeeman
 ! see LL04
 use atom_type, only : AtomType, ZeemanType, determinate, getorbital

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

 RETURN
 END SUBROUTINE Lande_eff

 SUBROUTINE ZeemanStrength()

 RETURN
 END SUBROUTINE ZeemanStrength

 SUBROUTINE ZeemanMultiplet(zm)
  type(ZeemanType), intent(out) :: zm
  
  zm%Ncomponent = 0
  zm%shift = 0d0
  zm%strength = 0d0
  zm%q = 0

 RETURN
 END SUBROUTINE ZeemanMultiplet

END MODULE Zeeman
