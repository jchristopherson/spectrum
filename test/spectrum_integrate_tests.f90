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
function test_integrate_boundaries() result(rst)
    logical :: rst

    real(real64), parameter :: dt = 0.5d0
    real(real64), parameter :: tol = 1.0d-12
    real(real64), parameter :: x1(1) = [2.0d0]
    real(real64), parameter :: x2(2) = [2.0d0, 2.0d0]
    real(real64), parameter :: x3(3) = [2.0d0, 2.0d0, 2.0d0]
    real(real64), parameter :: expected2(2) = [1.0d0, 2.0d0]
    real(real64), parameter :: expected3(3) = [1.0d0, 2.0d0, 3.0d0]
    real(real64), allocatable :: y1(:), y2(:), y3(:), y6(:)

    y1 = integrate(dt, x1, 1.0d0)
    y2 = integrate(dt, x2, 1.0d0)
    y3 = integrate(dt, x3, 1.0d0)
    y6 = integrate(dt, [x3, x3])
    rst = assert(y1, [1.0d0], tol) .and. assert(y2, expected2, tol) .and. &
        assert(y3, expected3, tol) .and. maxval(abs(y6 - &
        [0.0d0, 1.0d0, 2.0d0, 3.0d0, 4.0d0, 5.0d0])) < tol
    if (.not.rst) print "(A)", "TEST FAILED: test_integrate_boundaries -1"
end function

! ------------------------------------------------------------------------------
end module