module spectrum_diff_tests
    use iso_fortran_env
    use fortran_test_helper
    use spectrum
    implicit none

contains
! ------------------------------------------------------------------------------
function test_finite_difference() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    integer(int32) :: i
    real(real64) :: ti, t(n), x(n), dxdt(n), ans(n)

    ! Initialization
    rst = .true.
    ti = 0.0d0
    do i = 1, n
        t(i) = ti
        x(i) = sin(2.0d0 * ti)
        ans(i) = 2.0d0 * cos(2.0d0 * ti)
        ti = ti + dt
    end do

    ! Test
    dxdt = finite_difference(dt, x)
    if (.not.assert(ans(2:n-1), dxdt(2:n-1), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_finite_difference -1"
    end if

    dxdt = finite_difference(t, x)
    if (.not.assert(ans(2:n-1), dxdt(2:n-1), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_finite_difference -2"
    end if
end function

! ------------------------------------------------------------------------------
function test_stencil_diff() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    integer(int32) :: i
    real(real64) :: ti, t(n), x(n), dxdt(n), ans(n)

    ! Initialization
    rst = .true.
    ti = 0.0d0
    do i = 1, n
        t(i) = ti
        x(i) = sin(2.0d0 * ti)
        ans(i) = 2.0d0 * cos(2.0d0 * ti)
        ti = ti + dt
    end do

    ! Test
    dxdt = stencil_diff_5(dt, x)
    if (.not.assert(ans(3:n-2), dxdt(3:n-2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_stencil_diff -1"
    end if
end function

! ------------------------------------------------------------------------------
function test_stencil_diff_2() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: tol = 1.0d-4

    ! Local Variables
    integer(int32) :: i
    real(real64) :: ti, t(n), x(n), dxdt(n), ans(n)

    ! Initialization
    rst = .true.
    ti = 0.0d0
    do i = 1, n
        t(i) = ti
        x(i) = sin(2.0d0 * ti)
        ans(i) = -4.0d0 * sin(2.0d0 * ti)
        ti = ti + dt
    end do

    ! Test
    dxdt = stencil_second_diff_5(dt, x)
    if (.not.assert(ans(3:n-2), dxdt(3:n-2), tol)) then
        rst = .false.
        print "(A)", "TEST FAILED: test_stencil_diff_2 -1"
    end if
end function

! ------------------------------------------------------------------------------
end module