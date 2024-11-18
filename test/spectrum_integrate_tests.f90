module spectrum_integrate_tests
    use iso_fortran_env
    use fortran_test_helper
    use spectrum
    implicit none

contains
! ------------------------------------------------------------------------------
function test_integrate() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: tol = 1.0d-6

    ! Local Variables
    integer(int32) :: i
    real(real64) :: t(n), x(n), v(n), ans(n)

    ! Initialization
    rst = .true.
    t = (/ (i * dt, i = 0, n - 1) /)
    v = cos(5.0d0 * t)
    ans = 0.2d0 * sin(5.0d0 * t)

    ! Test
    x = integrate(dt, v)
    if (.not.assert(ans, x, tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_integrate -1"
    end if
end function

! ------------------------------------------------------------------------------
end module