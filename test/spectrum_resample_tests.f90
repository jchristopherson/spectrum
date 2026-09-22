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
module spectrum_resample_tests
    use iso_fortran_env
    use spectrum
    use fortran_test_helper
    implicit none
contains

function test_resample() result(rst)
    logical :: rst

    integer(int32), parameter :: npts = 8
    real(real64), parameter :: fs = 8.0d0
    real(real64), parameter :: tol = 1.0d-10
    real(real64) :: x(npts)
    real(real64), allocatable :: y(:), z(:), identity(:)

    x = 1.0d0
    y = upsample(2_int32, fs, x)
    z = downsample(2_int32, fs, x)
    identity = upsample(1_int32, fs, x)

    rst = allocated(y) .and. allocated(z) .and. allocated(identity)
    if (rst) rst = size(y) == 2 * npts .and. size(z) == npts / 2
    if (rst) rst = maxval(abs(y - 1.0d0)) < tol .and. &
        maxval(abs(z - 1.0d0)) < tol .and. assert(identity, x, tol)
    if (.not.rst) print '(A)', "TEST FAILED: test_resample -1"
end function

end module
