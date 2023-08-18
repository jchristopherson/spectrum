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
    integer(int32), parameter :: nfft = 4096
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

    ! Zero Padding
    deallocate(pwr)
    pwr = psd(win, x, fs = sample_rate_hz, nfft = nfft)

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, nfft)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_psd 1-2"
    end if
end function



function test_periodogram() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64), parameter :: tol = 1.0d-6
    integer(int32), parameter :: npts = 1024
    integer(int32), parameter :: nfft = 4096
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
    win%size = npts

    ! Construct a sinusoidal signal
    t = 0.0d0
    dt = 1.0d0 / sample_rate_hz
    x(1) = 0.0d0
    do i = 2, npts
        t = t + dt
        x(i) = amp1 * sin(2.0d0 * pi * freq1 * t)
    end do

    ! Compute the power spectrum
    ! pwr = psd(win, x, sample_rate_hz)
    pwr = periodogram(win, x, sample_rate_hz)

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, npts)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_periodogram 1-1"
    end if

    ! Test 2 - pad with zeros
    deallocate(pwr)
    pwr = periodogram(win, x, nfft = nfft)

    ! Compute the expected solution
    ans = amp1**2

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_periodogram 1-2"
    end if
end function



function test_csd() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64), parameter :: tol = 1.0d-6
    integer(int32), parameter :: npts = 10000
    integer(int32), parameter :: winsize = 1024
    integer(int32), parameter :: nfft = 4096
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
    pwr = abs(csd(win, x, x, sample_rate_hz))

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, winsize)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_csd 1-1"
    end if

    ! Zero padding
    deallocate(pwr)
    pwr = abs(csd(win, x, x, fs = sample_rate_hz, nfft = nfft))

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, nfft)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_csd 1-2"
    end if
end function



function test_cross_periodogram() result(rst)
    ! Arguments
    logical :: rst

    ! Local Variables
    real(real64), parameter :: tol = 1.0d-6
    integer(int32), parameter :: npts = 1024
    integer(int32), parameter :: nfft = 4096
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
    win%size = npts

    ! Construct a sinusoidal signal
    t = 0.0d0
    dt = 1.0d0 / sample_rate_hz
    x(1) = 0.0d0
    do i = 2, npts
        t = t + dt
        x(i) = amp1 * sin(2.0d0 * pi * freq1 * t)
    end do

    ! Compute the CSD of x with itself.  This is equivalent to the PSD of x.
    pwr = abs(cross_periodogram(win, x, x, sample_rate_hz))

    ! Compute the expected solution
    df = frequency_bin_width(sample_rate_hz, npts)
    ans = amp1**2 / df

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_cross_periodogram 1-1"
    end if

    ! Test zero padding
    deallocate(pwr)
    pwr = abs(cross_periodogram(win, x, x, nfft = nfft))

    ! Compute the expected solution
    ans = amp1**2

    ! Test
    p = maxval(pwr)
    if (.not.assert(p, ans, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_cross_periodogram 1-2"
    end if
end function



function test_spectrogram() result(rst)
    ! Arguments
    logical :: rst

    ! Parameters
    integer(int32), parameter :: window_size = 512
    real(real64), parameter :: fs = 2048.0d0
    real(real64), parameter :: f0 = 1.0d2
    real(real64), parameter :: f1 = 1.0d3
    real(real64), parameter :: duration = 50.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: tol = 1.0d-1

    ! Local Variables
    integer(int32) :: i, npts
    real(real64) :: k, check, dt
    complex(real64), allocatable, dimension(:,:) :: r
    real(real64), allocatable, dimension(:) :: t, x
    real(real64), allocatable, dimension(:,:) :: mag
    type(flat_top_window) :: win

    ! Initialization
    rst = .true.

    ! Create the exponential chirp signal
    npts = floor(duration * fs) + 1
    allocate(t(npts))
    dt = duration / (npts - 1.0d0)
    t = (/ (dt * i, i = 0, npts - 1) /)
    k = (f1 / f0)**(1.0 / duration)
    x = sin(2.0d0 * pi * f0 * (k**t - 1.0d0) / log(k))

    ! Define the window
    win%size = window_size

    ! Compute the spectrogram of x
    r = spectrogram(win, x)

    ! Compute the magnitude
    mag = abs(r)

    ! Ensure the magnitude of the largest component of each transform is one
    do i = 1, size(mag, 2)
        check = maxval(mag(:,i))
        if (.not.assert(check, 1.0d0, tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_spectrogram 1-2"
        end if
    end do
end function

end module