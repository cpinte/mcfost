module zeeman
  ! see LL04
    use atom_type, only : AtomType, ZeemanType, parse_label, getorbital, AtomicLine
    use atmos_type, only : lmagnetized
    use messages
    use parametres
    use math, only : w3js

    use mcfost_env, only : dp

    implicit none


    contains


    function Lande(S, L, J)
	!Return the Landé factor assuming LS coupling
        real, intent(in) :: S
        integer, intent(in) :: L!, J
        real, intent(in) :: J
        real :: Lande

        if (J == 0.0) then
            Lande = 0.0
        else
            Lande = 1.5 + (S*(S + 1.0) - L*(L + 1)) / (2.0*J*(J + 1.0))
        end if

	    return
    end function Lande

    subroutine Lande_eff(atom, kr)
	    type (AtomType), intent(inout) :: atom
        integer :: kr
        integer :: i, j
        real g, gi, gj

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
        if (lmagnetized) then
            write(*,'(" -> line g(lower)="(1F12.5)" g(upper)="(1F12.3)" geff="(1F12.5))') &
         	    gi, gj,  g
        endif

        return
    end subroutine Lande_eff
    
    function ZeemanStrength(Ji, Mi, Jj, Mj)
        real(kind=dp) :: ZeemanStrength
        real :: dM
        real, intent(in) :: Ji, Jj, Mi, Mj
        integer :: q

    !q = -1 = sB, 0 spi, +1 sr
    !dM = 1 -> sr; dM = +1 sb
        dM = Mj - Mi
        q = -int(dM)
        if (abs(dM) > 1) then
            call warning(" dM is not satisfying selection rules!")
            write(*,*) q, Ji, Jj, Mi, Mj, dM
            ZeemanStrength = 0.0
            return !should not happen
        end if

        ZeemanStrength = 3.0 * w3js(int(2*Jj),int(2*Ji),2,&
            -int(2*Mj), int(2*Mi),-2*q)**2
         
         
     !Using table 3.1 of LL04
     !... 
         
	    return
    end function ZeemanStrength

    subroutine ZeemanMultiplet(line)
        type(AtomicLine), intent(inout) :: line
        integer :: nc, i1, i2
        real :: Mi, Mj

	!Line is always polarizable
	!The unpolarized profile function is used instead of a ZeemanMulitplet
	!with only 1 component! (Sum of all components but unshifted actually)
	    if (.not.line%polarizable) call error("Line must be polarizable to have ZM!")

	!Effective Zeeman triplet using effective Landé
        if (line%ZeemanPattern == -1) then
    
            line%zm%Ncomponent = 3
            allocate(line%zm%q(3), line%zm%strength(3), line%zm%shift(3))
            line%zm%q = (/-1, 0, 1/)
            line%zm%strength = 1.0 !Here all components have the same strength
            line%zm%shift = line%zm%q * line%g_Lande_eff
            !or line%zm%shift = (/-line%g_i, 0.0, line%g_j/)
            write(*,'(" -> line "(1I2)" -> "(1I2)" has"(1I4)" components (EZT)!")') &
                line%i, line%j, line%zm%Ncomponent
            !    write(*,*) line%zm%q
            !    write(*,*) line%zm%shift
            !    write(*,*) line%zm%strength

            !Full Zeeman pattern
        else if (line%ZeemanPattern == 1) then
            !First count number of components
            line%zm%Ncomponent = 0
            do i1=1,int(2*line%atom%qJ(line%i)+1)
                Mi = line%atom%qJ(line%i) + 1 - i1
                do i2=1,int(2*line%atom%qJ(line%j)+1)
                    Mj = line%atom%qJ(line%j) + 1 - i2
                    if (abs(Mj - Mi) <= 1.0) line%zm%Ncomponent = line%zm%Ncomponent + 1
                enddo
            enddo

            allocate(line%zm%q(line%zm%Ncomponent), &
                    line%zm%strength(line%zm%Ncomponent), &
                    line%zm%shift(line%zm%Ncomponent))

            write(*,'(" -> line "(1I2)" -> "(1I2)" has"(1I4)" components (FZM)!")') &
                line%i, line%j, line%zm%Ncomponent
            write(*,'("   "(1A2)"="(1F12.2),(1A2)"="(1F12.2))') "J'", line%atom%qJ(line%j), "J", line%atom%qJ(line%i)
            nc = 0

            do i1=1,int(2*line%atom%qJ(line%i)+1)
            Mi = line%atom%qJ(line%i) + 1 - i1
                do i2=1,int(2*line%atom%qJ(line%j)+1)
                    Mj = line%atom%qJ(line%j) + 1 - i2
                    if (abs(Mj - Mi) <= 1.0) then
                        nc = nc + 1
                        line%zm%q(nc) = -int(Mj - Mi)
                        line%zm%shift(nc) = line%glande_j * Mj - line%glande_i * Mi
                        line%zm%strength(nc) = ZeemanStrength(line%atom%qJ(line%i),Mi,&
                                            line%atom%qJ(line%j), Mj)
                        write(*,'(" +++> q="(1I2)" shift="(1F12.6)" S="(1F12.3))') &
                                line%zm%q(nc), line%zm%shift(nc), line%zm%strength(nc)
                    end if
                end do
            end do

!     else if (.not.line%polarizable) then
!        allocate(line%zm%q(1), line%zm%strength(1), line%zm%shift(1))
!        line%zm%Ncomponent = 1
!        line%zm%q = 0.0
!        line%zm%strength = 0.0 !not used in this case.
!        line%zm%shift = 0.0

        else
            call Error("Zeeman pattern type unknown!")
            write(*,*) line%ZeemanPattern, line%i, line%j
        end if

        return
    end subroutine ZeemanMultiplet

!to do
    subroutine write_zeeman_multiplet(atom)
        !write ZM structure for all polarizable lines of an atom
        !to ascii file atom.zm.
        !This is use for debug!
        type(atomtype), intent(in) :: atom
	
	
        return
    end subroutine write_zeeman_multiplet

end module Zeeman