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
module spectrum_convolution_tests
    use iso_fortran_env
    use spectrum
    use fortran_test_helper
    implicit none
contains

function test_convolution() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64) :: u1(3), v1(2), u2(7), v2(4), u3(3), v3(7), &
        ans1(4), ans2(7), ans3(9)
    real(real64), allocatable, dimension(:) :: w1, w2, w3, w4

    ! Initialization
    rst = .true.
    u1 = [1.0d0, 0.0d0, 1.0d0]
    v1 = [2.0d0, 7.0d0]
    u2 = [-1.0d0, 2.0d0, 3.0d0, -2.0d0, 0.0d0, 1.0d0, 2.0d0]
    v2 = [2.0d0, 4.0d0, -1.0d0, 1.0d0]
    u3 = [1.0d0, 1.0d0, 1.0d0]
    v3 = [1.0d0, 1.0d0, 0.0d0, 0.0d0, 0.0d0, 1.0d0, 1.0d0]

    ! Define the solutions
    ans1 = [2.0d0, 7.0d0, 2.0d0, 7.0d0]
    ans2 = [15.0d0, 5.0d0, -9.0d0, 7.0d0, 6.0d0, 7.0d0, -1.0d0]
    ans3 = [1.0d0, 2.0d0, 2.0d0, 1.0d0, 0.0d0, 1.0d0, 2.0d0, 2.0d0, 1.0d0]

    ! Test 1
    w1 = convolve(u1, v1)
    if (.not.assert(ans1, w1)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_convolution 1-1"
    end if

    ! Test 2
    w2 = convolve(u2, v2, method = SPCTRM_CENTRAL_CONVOLUTION)
    if (.not.assert(ans2, w2)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_convolution 1-2"
    end if

    ! Test 3
    w3 = convolve(u3, v3)
    if (.not.assert(ans3, w3)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_convolution 1-3"
    end if

    w4 = deconvolve(ans1, v1)
    if (.not.assert([1.0d0, 0.0d0, 1.0d0, 0.0d0, 0.0d0], w4, 1.0d-10)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_convolution 1-4"
    end if
end function

end module