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
