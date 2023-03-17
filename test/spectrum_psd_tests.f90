module spectrum_psd_tests
    use iso_fortran_env
    use spectrum
    use fortran_test_helper
    implicit none

contains

function test_psd() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64), parameter :: tol = 1.0d-6
    integer(int32), parameter :: npts = 10000
    integer(int32), parameter :: winsize = 1024
    real(real64), parameter :: sample_rate_hz = 1.024d3
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: freq1 = 2.5d1
    real(real64), parameter :: amp1 = 2.0d0
    integer(int32) :: i
    real(real64) :: t, dt, ans, df, p
    real(real64) :: x(npts)
    real(real64), allocatable :: pwr(:)
    type(hann_window) :: win

    ! Initialization
    rst = .true.
    win%size = winsize

    ! Construct a sinusoidal signal
    t = 0.0d0
    dt = 1.0d0 / sample_rate_hz
    x(1) = 0.0d0
    do i = 2, npts
        t = t + dt
        x(i) = amp1 * sin(2.0d0 * pi * freq1 * t)
    end do

    ! Compute the power spectrum
    pwr = psd(win, x, sample_rate_hz)

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, winsize)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_psd 1-1"
    end if
end function




function test_csd() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64), parameter :: tol = 1.0d-6
    integer(int32), parameter :: npts = 10000
    integer(int32), parameter :: winsize = 1024
    real(real64), parameter :: sample_rate_hz = 1.024d3
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: freq1 = 2.5d1
    real(real64), parameter :: amp1 = 2.0d0
    integer(int32) :: i
    real(real64) :: t, dt, ans, df, p
    real(real64) :: x(npts)
    real(real64), allocatable :: pwr(:)
    type(hann_window) :: win

    ! Initialization
    rst = .true.
    win%size = winsize

    ! Construct a sinusoidal signal
    t = 0.0d0
    dt = 1.0d0 / sample_rate_hz
    x(1) = 0.0d0
    do i = 2, npts
        t = t + dt
        x(i) = amp1 * sin(2.0d0 * pi * freq1 * t)
    end do

    ! Compute the CSD of x with itself.  This is equivalent to the PSD of x.
    pwr = csd(win, x, x, sample_rate_hz)

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, winsize)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_csd 1-1"
    end if
end function

end module