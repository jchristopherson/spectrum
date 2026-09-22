! MIT License
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
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
function test_tvr_derivative_sparse() result(rst)
    logical :: rst

    integer(int32), parameter :: n = 1000
    real(real64), parameter :: dt = 1.0d-3
    real(real64), parameter :: alpha = 1.0d0
    real(real64), parameter :: tol = 1.0d-5
    real(real64) :: x(n), dxdt(n)
    integer(int32) :: i, niter

    do i = 1, n
        x(i) = 2.0d0 * real(i - 1, real64) * dt + 3.0d0
    end do

    dxdt = tvr_derivative(dt, x, alpha, use_sparse = .true., niter = niter)
    rst = niter > 0 .and. all(abs(dxdt - 2.0d0) < tol)
    if (.not.rst) print "(A)", "TEST FAILED: test_tvr_derivative_sparse -1"
end function

! ------------------------------------------------------------------------------
end module