program test_stokes_scattering
  ! Unit test: Stokes vector after successive scatterings (scattering, update_Stokes)
  ! against an independent reference in Jones calculus.
  !
  ! The reference carries the coherency matrix J of the beam in an explicit basis
  ! (e_a = projected z axis, e_b = k x e_a) and applies the complex scattering
  ! amplitudes (S2 parallel, S1 perpendicular to the scattering plane) in the basis
  ! (e_par, e_perp), e_perp = k_in x k_out normalised, e_par = e_perp x k.
  ! update_Stokes receives the Mueller matrix built from the same amplitudes with the
  ! MCFOST definitions (mueller_Mie). Quantities independent of the reference frames
  ! are compared after each scattering: I, the linearly polarised flux sqrt(Q^2+U^2)
  ! and |V| (the sign of V depends on the convention).
  ! Built by "make unit_tests" in src/ and run by test_suite/test_unit.py.
  ! Exits with status 1 on failure.

  use mcfost_env, only : dp
  use scattering, only : update_Stokes

  implicit none

  integer, parameter :: n_seq = 50, n_scatt = 5
  real(kind=dp), parameter :: tol = 1e-10_dp

  integer :: n_fail, n_test, iseq, i
  real(kind=dp), dimension(3,0:n_scatt) :: k
  complex(kind=dp), dimension(n_scatt) :: S1, S2
  complex(kind=dp), dimension(2,2) :: J
  real(kind=dp), dimension(4) :: S_ref, S_mc
  real(kind=dp), dimension(4,4) :: M
  real(kind=dp) :: err_I, err_P, err_V
  integer, allocatable :: seed(:)

  n_fail = 0 ; n_test = 0
  err_I = 0.0_dp ; err_P = 0.0_dp ; err_V = 0.0_dp

  call random_seed(size=i)
  allocate(seed(i)) ; seed = 12345
  call random_seed(put=seed)

  do iseq=1, n_seq
     ! the scattering plane is undefined for forward or backward scattering: such
     ! (measure zero) draws are redone
     k(:,0) = random_direction()
     do i=1, n_scatt
        k(:,i) = random_direction()
        do while (abs(dot_product(k(:,i), k(:,i-1))) > 0.999999_dp)
           k(:,i) = random_direction()
        enddo
     enddo
     do i=1, n_scatt
        S1(i) = random_complex() ; S2(i) = random_complex()
     enddo

     ! unpolarised beam, I = 1
     J = (0.0_dp,0.0_dp) ; J(1,1) = (0.5_dp,0.0_dp) ; J(2,2) = (0.5_dp,0.0_dp)
     S_mc = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

     do i=1, n_scatt
        call scatter_reference(J, k(:,i-1), k(:,i), S1(i), S2(i))
        S_ref = stokes(J)

        call mcfost_mueller(S1(i), S2(i), M)
        call update_Stokes(S_mc, k(1,i-1),k(2,i-1),k(3,i-1), k(1,i),k(2,i),k(3,i), M)

        err_I = max(err_I, abs(S_mc(1) - S_ref(1)))
        err_P = max(err_P, abs(hypot(S_mc(2),S_mc(3)) - hypot(S_ref(2),S_ref(3))))
        err_V = max(err_V, abs(abs(S_mc(4)) - abs(S_ref(4))))
     enddo
  enddo

  call check("intensity", err_I)
  call check("linear polarisation", err_P)
  call check("circular polarisation |V|", err_V)

  write(*,'(a,i0,a,i0)') "test_stokes_scattering: failures = ", n_fail, " / ", n_test
  if (n_fail > 0) stop 1

contains

  subroutine check(name, err)
    character(len=*), intent(in) :: name
    real(kind=dp), intent(in) :: err

    n_test = n_test + 1
    if (err > tol) then
       n_fail = n_fail + 1
       write(*,'(a,a,a,es10.3,a,i0,a,i0,a)') "FAIL ", name, ": max error ", err, &
            " over ", n_seq, " sequences of ", n_scatt, " scatterings"
    endif
  end subroutine check

  function random_direction() result(v)
    real(kind=dp), dimension(3) :: v
    real(kind=dp) :: r(2), cost, phi, sint
    call random_number(r)
    cost = 2.0_dp*r(1) - 1.0_dp
    phi = 2.0_dp * acos(-1.0_dp) * r(2)
    sint = sqrt(1.0_dp - cost**2)
    v = [sint*cos(phi), sint*sin(phi), cost]
  end function random_direction

  function random_complex() result(c)
    complex(kind=dp) :: c
    real(kind=dp) :: r(2)
    call random_number(r)
    c = cmplx(2.0_dp*r(1) - 1.0_dp, 2.0_dp*r(2) - 1.0_dp, kind=dp)
  end function random_complex

  function cross(a, b) result(c)
    real(kind=dp), dimension(3), intent(in) :: a, b
    real(kind=dp), dimension(3) :: c
    c = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), a(1)*b(2) - a(2)*b(1)]
  end function cross

  subroutine beam_basis(kv, ea, eb)
    ! e_a : projected z axis, e_b = k x e_a
    real(kind=dp), dimension(3), intent(in) :: kv
    real(kind=dp), dimension(3), intent(out) :: ea, eb
    if (abs(kv(3)) > 0.999999_dp) then ! along z: any perpendicular axis
       ea = [1.0_dp, 0.0_dp, 0.0_dp] - kv(1) * kv
    else
       ea = [0.0_dp, 0.0_dp, 1.0_dp] - kv(3) * kv
    endif
    ea = ea / norm2(ea)
    eb = cross(kv, ea)
  end subroutine beam_basis

  subroutine scatter_reference(J, k_in, k_out, S1, S2)
    complex(kind=dp), dimension(2,2), intent(inout) :: J
    real(kind=dp), dimension(3), intent(in) :: k_in, k_out
    complex(kind=dp), intent(in) :: S1, S2

    real(kind=dp), dimension(3) :: a_in, b_in, a_out, b_out, perp, par_in, par_out
    real(kind=dp), dimension(2,2) :: P_in, P_out
    complex(kind=dp), dimension(2,2) :: A, T

    call beam_basis(k_in, a_in, b_in)
    call beam_basis(k_out, a_out, b_out)
    perp = cross(k_in, k_out) ; perp = perp / norm2(perp)
    par_in = cross(perp, k_in)
    par_out = cross(perp, k_out)

    ! (a,b) of the incoming beam -> (par, perp)
    P_in = reshape([dot_product(par_in,a_in), dot_product(perp,a_in), &
         dot_product(par_in,b_in), dot_product(perp,b_in)], [2,2])
    ! (par, perp) -> (a,b) of the outgoing beam
    P_out = reshape([dot_product(a_out,par_out), dot_product(b_out,par_out), &
         dot_product(a_out,perp), dot_product(b_out,perp)], [2,2])
    A = (0.0_dp,0.0_dp) ; A(1,1) = S2 ; A(2,2) = S1

    T = matmul(P_out, matmul(A, P_in))
    J = matmul(T, matmul(J, conjg(transpose(T))))
    J = J / real(J(1,1) + J(2,2), kind=dp) ! normalised to I = 1, as update_Stokes
  end subroutine scatter_reference

  function stokes(J) result(S)
    complex(kind=dp), dimension(2,2), intent(in) :: J
    real(kind=dp), dimension(4) :: S
    S(1) = real(J(1,1) + J(2,2), kind=dp)
    S(2) = real(J(1,1) - J(2,2), kind=dp)
    S(3) = 2.0_dp * real(J(1,2), kind=dp)
    S(4) = -2.0_dp * aimag(J(1,2))
  end function stokes

  subroutine mcfost_mueller(S1, S2, M)
    ! Mueller matrix from the amplitudes with the definitions of mueller_Mie, normalised
    ! by s11 and arranged as in get_Mueller_matrix_per_cell
    complex(kind=dp), intent(in) :: S1, S2
    real(kind=dp), dimension(4,4), intent(out) :: M
    real(kind=dp) :: s11, s12, s33, s34
    complex(kind=dp) :: s

    s11 = 0.5_dp * (abs(S2)**2 + abs(S1)**2)
    s12 = 0.5_dp * (abs(S2)**2 - abs(S1)**2)
    s = S2 * conjg(S1)
    s33 = real(s, kind=dp) ; s34 = aimag(s)

    M = 0.0_dp
    M(1,1) = 1.0_dp ; M(2,2) = 1.0_dp
    M(1,2) = s12/s11 ; M(2,1) = M(1,2)
    M(3,3) = s33/s11 ; M(4,4) = s33/s11
    M(3,4) = -s34/s11 ; M(4,3) = -M(3,4)
  end subroutine mcfost_mueller

end program test_stokes_scattering
