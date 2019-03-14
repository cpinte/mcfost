MODULE zeeman
 ! see LL04
 use atom_type, only : AtomType, ZeemanType, determinate, getorbital, AtomicLine
 use atmos_type, only : atmos
 use messages
 use parametres

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
   if (atmos%magnetized) & !otherwise we don't care seeing that
   	write(*,*) " -> line gi = ", gi, " line gj = ", gj," line geff = ", g

 RETURN
 END SUBROUTINE Lande_eff

 SUBROUTINE ZeemanStrength()


 RETURN
 END SUBROUTINE ZeemanStrength

 SUBROUTINE ZeemanMultiplet(line) !Called only once per line
  type(AtomicLine), intent(inout) :: line  
  
  
!    if (line%g_Lande_eff <= -99) then
!     write(*,*) "Should not happen anymore"
!     CALL Warning("  -> Landé factor not parsed for this line")
!     !write(*,*) line%g_lande_eff, line%atom%qJ(line%i), line%atom%qJ(line%j)
!     RETURN
!    end if

  if (line%ZeemanPattern == 0) then !Weak Field
   allocate(line%zm%q(1), line%zm%strength(1), line%zm%shift(1)) 
   
   line%zm%Ncomponent = 0
   line%zm%q = 1d0
   line%zm%strength = 1d0 !not used in this case.
   line%zm%shift = line%zm%q*line%g_lande_eff !unique shift = deltaLamB/B
   
   !The Stokes parameters are only prop to g_lande_eff * dI/dlambda

  else if (line%ZeemanPattern == -1) then
   line%zm%Ncomponent = 3
   line%zm%q = (/-1, 0, 1/)
   line%zm%strength = 1d0 !Here all components have the same strength
   line%zm%shift = line%zm%q * line%g_Lande_eff !same shift
   CALL Error("Effective Zeeman Triplet (EZT) not implemented yet!")   
   write(*,*) "  Line ", line%j,"->",line%i," has", line%zm%Ncomponent,&
   			  " Zeeman components, geff=", line%g_lande_eff
   			  
  else if (line%ZeemanPattern == 1) then
    !Strength relative of all components
   CALL Error("Full Zeeman components not implemented yet!")
  else
    CALL Error("Zeeman components Recipe unknown!")
  end if

 RETURN
 END SUBROUTINE ZeemanMultiplet

END MODULE Zeeman
