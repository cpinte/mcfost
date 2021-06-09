MODULE zeeman
  ! see LL04
  use atom_type, only : AtomType, ZeemanType, determinate, getorbital, AtomicLine
  use atmos_type, only : lmagnetized
  use messages
  use parametres

  use mcfost_env, only : dp

  IMPLICIT NONE
CONTAINS

  FUNCTION Lande(S, L, J) result(g)
    real(kind=dp) :: S
    integer :: L!, J
    real(kind=dp) :: g, J

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
    real(kind=dp) g, gi, gj

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
    if (lmagnetized) & !otherwise we don't care seeing that
         write(*,*) " -> line gi = ", gi, " line gj = ", gj!," line geff = ", g

    RETURN
  END SUBROUTINE Lande_eff

  FUNCTION ZeemanStrength(Ji, Mi, Jj, Mj)
    use math, only : w3js
    real(kind=dp) :: ZeemanStrength, dM
    !   integer, intent(in) :: Ji, Jj, Mi, Mj
    real(kind=dp), intent(in) :: Ji, Jj, Mi, Mj
    integer :: q

    !q = -1 = sB, 0 spi, +1 sr
    !dM = 1 -> sr; dM = +1 sb
    dM = Mj - Mi
    q = -int(dM)
    if (abs(dM) > 1) then
       write(*,*) dM, " is not satisfying selection rules!"
       ZeemanStrength = 0d0
       RETURN !should not happen
    end if

    ZeemanStrength = 3.0 * w3js(int(2*Jj),int(2*Ji),2,&
         -int(2*Mj),int(2*Mi),-2*q)**2
    RETURN
  END FUNCTION ZeemanStrength

  SUBROUTINE ZeemanMultiplet(line) !Called only once per line
    type(AtomicLine), intent(inout) :: line
    integer :: nc, i1, i2
    real(kind=dp) :: Mi, Mj
    !, norm(3) !sum of -1, 0 and +1 components
    !not need, j-symbols normalised

    if (line%ZeemanPattern == -1 .and. line%polarizable) then
       line%zm%Ncomponent = 3
       allocate(line%zm%q(3), line%zm%strength(3), line%zm%shift(3))
       line%zm%q = (/-1, 0, 1/)
       line%zm%strength = 1d0 !Here all components have the same strength
       line%zm%shift = line%zm%q * line%g_Lande_eff !same shift
       write(*,*) "  Line ", line%j,"->",line%i," has", line%zm%Ncomponent,&
            " Zeeman components, geff=", line%g_lande_eff
       !    write(*,*) line%zm%q
       !    write(*,*) line%zm%shift
       !    write(*,*) line%zm%strength

    else if (line%ZeemanPattern == 1 .and. line%polarizable) then
       !Strength relative of all components
       !First count number of components
       line%zm%Ncomponent = 0
       !    do i1=1,2*line%atom%qJ(line%j)+1
       !     Mj = line%atom%qJ(line%j) + 1 - i1
       !     do i2=1,2*line%atom%qJ(line%i)+1
       !      Mi = line%atom%qJ(line%i) + 1 - i2
       do Mj=-line%atom%qJ(line%j),line%atom%qJ(line%j)
          do Mi=-line%atom%qJ(line%i),line%atom%qJ(line%i)
             if (abs(Mi-Mj) <= 1) line%zm%Ncomponent = line%zm%Ncomponent + 1
          end do
       end do

       allocate(line%zm%q(line%zm%Ncomponent), &
            line%zm%strength(line%zm%Ncomponent), &
            line%zm%shift(line%zm%Ncomponent))

       write(*,*) "  Line ", line%j,"->",line%i," has", line%zm%Ncomponent,&
            " Zeeman components, geff=", line%g_lande_eff
       write(*,*) "J' = ", line%atom%qJ(line%j), " J = ", line%atom%qJ(line%i)
       nc = 0
       do Mi=-line%atom%qJ(line%i),line%atom%qJ(line%i)
          do Mj=-line%atom%qJ(line%j),line%atom%qJ(line%j)
             if (abs(Mj-Mi) <= 1) then
                nc = nc + 1
                line%zm%q(nc) = -int(Mj - Mi)
                line%zm%shift(nc) = line%glande_i * Mi - line%glande_j * Mj
                line%zm%strength(nc) = ZeemanStrength(line%atom%qJ(line%i),Mi,line%atom%qJ(line%j), Mj)
                write(*,*) line%zm%q(nc), line%zm%shift(nc), "Strength = ",line%zm%strength(nc)
             end if
          end do
       end do

    else if (.not.line%polarizable) then !unpolarized line
       allocate(line%zm%q(1), line%zm%strength(1), line%zm%shift(1))
       line%zm%Ncomponent = 1
       line%zm%q = 0d0
       line%zm%strength = 0d0 !not used in this case.
       line%zm%shift = 0d0

    else
       CALL Error("Zeeman components Recipe unknown!")
    end if

    RETURN
  END SUBROUTINE ZeemanMultiplet

END MODULE Zeeman
