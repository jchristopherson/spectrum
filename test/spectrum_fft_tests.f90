module spectrum_fft_tests
    use iso_fortran_env
    use spectrum
    use fortran_test_helper
    implicit none

    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: tol = 1.0d-10

contains
! ------------------------------------------------------------------------------
pure function reference_dft(x) result(rst)
    !! A direct, brute-force evaluation of the positive half of the discrete
    !! Fourier transform used as an independent reference for rfft.
    real(real64), intent(in), dimension(:) :: x
    complex(real64), allocatable, dimension(:) :: rst

    integer(int32) :: i, j, n, m

    n = size(x)
    m = compute_transform_length(n)
    allocate(rst(m), source = (0.0d0, 0.0d0))
    do i = 1, m
        do j = 1, n
            rst(i) = rst(i) + x(j) * exp(cmplx(0.0d0, &
                -2.0d0 * pi * (i - 1.0d0) * (j - 1.0d0) / n, real64))
        end do
    end do
end function

! ------------------------------------------------------------------------------
function assert_transform(ans, x, id) result(rst)
    !! Compares the real and imaginary components of two complex-valued arrays.
    complex(real64), intent(in), dimension(:) :: ans, x
    character(len = *), intent(in) :: id
    logical :: rst

    rst = .true.
    if (size(ans) /= size(x)) then
        rst = .false.
        print '(A)', "TEST FAILED: " // id // " (size mismatch)"
        return
    end if
    if (.not.assert(real(ans, real64), real(x, real64), tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: " // id // " (real component)"
    end if
    if (.not.assert(aimag(ans), aimag(x), tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: " // id // " (imaginary component)"
    end if
end function

! ------------------------------------------------------------------------------
function test_rfft_even_length() result(rst)
    !! Tests rfft for an even-length signal.  The transform of an N-point
    !! signal has N / 2 + 1 points, and both the DC and Nyquist terms are
    !! purely real.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 8
    integer(int32), parameter :: m = n / 2 + 1

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), a, b
    complex(real64), allocatable, dimension(:) :: y
    complex(real64) :: ans(m)

    ! Initialization
    rst = .true.
    call random_number(x)

    ! Test 1 - Compare against a brute-force DFT
    y = rfft(x)
    ans = reference_dft(x)
    if (.not.assert_transform(ans, y, "test_rfft_even_length 1-1")) rst = .false.

    ! Test 2 - The DC term is the sum of the signal, and the Nyquist term is
    ! the alternating sum.  Both must be purely real.
    if (.not.assert(sum(x), real(y(1), real64), tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_even_length 1-2a"
    end if
    a = 0.0d0
    do i = 1, n
        a = a + x(i) * (-1.0d0)**(i - 1)
    end do
    if (.not.assert(a, real(y(m), real64), tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_even_length 1-2b"
    end if
    if (.not.assert(0.0d0, aimag(y(1)), tol) .or. &
        .not.assert(0.0d0, aimag(y(m)), tol)) &
    then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_even_length 1-2c"
    end if

    ! Test 3 - A pure sinusoid landing exactly on bin 2 must produce a single
    ! non-zero, non-DC term of magnitude a * n / 2.
    a = 3.0d0
    do i = 1, n
        x(i) = a * cos(2.0d0 * pi * 2.0d0 * (i - 1.0d0) / n)
    end do
    y = rfft(x)
    ans = (0.0d0, 0.0d0)
    ans(3) = cmplx(0.5d0 * a * n, 0.0d0, real64)
    if (.not.assert_transform(ans, y, "test_rfft_even_length 1-3")) rst = .false.

    ! Test 4 - A sine on bin 2 produces a purely imaginary term of -a * n / 2.
    b = 2.0d0
    do i = 1, n
        x(i) = b * sin(2.0d0 * pi * 2.0d0 * (i - 1.0d0) / n)
    end do
    y = rfft(x)
    ans = (0.0d0, 0.0d0)
    ans(3) = cmplx(0.0d0, -0.5d0 * b * n, real64)
    if (.not.assert_transform(ans, y, "test_rfft_even_length 1-4")) rst = .false.
end function

! ------------------------------------------------------------------------------
function test_rfft_odd_length() result(rst)
    !! Tests rfft for an odd-length signal.  The transform of an N-point signal
    !! has (N + 1) / 2 points, and there is no Nyquist term.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 7
    integer(int32), parameter :: m = (n + 1) / 2

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), a, b
    complex(real64), allocatable, dimension(:) :: y
    complex(real64) :: ans(m)

    ! Initialization
    rst = .true.
    call random_number(x)

    ! Test 1 - Compare against a brute-force DFT
    y = rfft(x)
    ans = reference_dft(x)
    if (.not.assert_transform(ans, y, "test_rfft_odd_length 1-1")) rst = .false.

    ! Test 2 - The DC term is the sum of the signal, and must be purely real
    if (.not.assert(sum(x), real(y(1), real64), tol) .or. &
        .not.assert(0.0d0, aimag(y(1)), tol)) &
    then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_odd_length 1-2"
    end if

    ! Test 3 - A sinusoid landing on the highest available bin retains both
    ! real and imaginary components as there is no Nyquist bin.
    a = 1.5d0
    b = -0.5d0
    do i = 1, n
        x(i) = a * cos(2.0d0 * pi * (m - 1.0d0) * (i - 1.0d0) / n) - &
            b * sin(2.0d0 * pi * (m - 1.0d0) * (i - 1.0d0) / n)
    end do
    y = rfft(x)
    ans = (0.0d0, 0.0d0)
    ans(m) = cmplx(0.5d0 * a * n, 0.5d0 * b * n, real64)
    if (.not.assert_transform(ans, y, "test_rfft_odd_length 1-3")) rst = .false.
end function

! ------------------------------------------------------------------------------
function test_rfft_padded_length() result(rst)
    !! Tests the optional length argument of rfft for both zero-padding and
    !! truncation of the input signal.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: n = 5
    integer(int32), parameter :: npad = 16
    integer(int32), parameter :: ntrunc = 4

    ! Local Variables
    real(real64) :: x(n), xpad(npad)
    complex(real64), allocatable, dimension(:) :: y, ans

    ! Initialization
    rst = .true.
    call random_number(x)
    xpad = 0.0d0
    xpad(1:n) = x

    ! Test 1 - Zero-padding to a longer length must match the transform of the
    ! explicitly padded signal
    y = rfft(x, npad)
    if (size(y) /= npad / 2 + 1) then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_padded_length 1-1a"
        return
    end if
    ans = reference_dft(xpad)
    if (.not.assert_transform(ans, y, "test_rfft_padded_length 1-1b")) &
        rst = .false.

    ! Test 2 - Zero-padding preserves the DC term
    if (.not.assert(sum(x), real(y(1), real64), tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_padded_length 1-2"
    end if

    ! Test 3 - Requesting the actual signal length must match the unpadded
    ! transform
    y = rfft(x, n)
    ans = rfft(x)
    if (.not.assert_transform(ans, y, "test_rfft_padded_length 1-3")) &
        rst = .false.

    ! Test 4 - A length shorter than the signal truncates the signal
    y = rfft(x, ntrunc)
    if (size(y) /= ntrunc / 2 + 1) then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_padded_length 1-4a"
        return
    end if
    ans = reference_dft(x(1:ntrunc))
    if (.not.assert_transform(ans, y, "test_rfft_padded_length 1-4b")) &
        rst = .false.

    ! Test 5 - A zero-padded transform can be inverted back to the padded
    ! signal (scaled by the transform length as irfft is unnormalized)
    if (.not.assert(npad * xpad, irfft(rfft(x, npad)), tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_rfft_padded_length 1-5"
    end if
end function

! ------------------------------------------------------------------------------
function reference_stft(win, x, nk) result(rst)
    !! Builds the expected STFT by explicitly windowing each segment and
    !! applying rfft, then applying the amplitude scaling documented for stft.
    class(window), intent(in) :: win
    real(real64), intent(in), dimension(:) :: x
    integer(int32), intent(in) :: nk
    complex(real64), allocatable, dimension(:,:) :: rst

    integer(int32) :: i, k, m, nx, nxfrm, i1
    real(real64) :: del, sumw, fac
    real(real64), allocatable, dimension(:) :: seg, w
    complex(real64), allocatable, dimension(:) :: xfrm

    nx = size(x)
    m = win%size
    nxfrm = compute_transform_length(m)
    allocate(rst(nxfrm, nk), seg(m), w(m))
    do k = 1, m
        w(k) = win%evaluate(k - 1)
    end do
    sumw = sum(w)
    fac = 2.0d0 / sumw    ! (m / sumw) * (2 / m)
    if (nk > 1) then
        del = (nx - m) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if

    do i = 1, nk
        i1 = int((i - 1) * del + 0.5d0, int32)
        seg = w * x(i1+1:i1+m)
        xfrm = rfft(seg)
        rst(:,i) = fac * xfrm
        rst(1,i) = 0.5d0 * fac * cmplx(real(xfrm(1), real64), 0.0d0, real64)
        if (mod(m, 2) == 0) then
            rst(nxfrm,i) = 0.5d0 * fac * &
                cmplx(real(xfrm(nxfrm), real64), 0.0d0, real64)
        end if
    end do
end function

! ------------------------------------------------------------------------------
function check_stft(win, x, nk, id) result(rst)
    !! Verifies the transform, the size of the output, and the segment offsets.
    class(window), intent(in) :: win
    real(real64), intent(in), dimension(:) :: x
    integer(int32), intent(in) :: nk
    character(len = *), intent(in) :: id
    logical :: rst

    integer(int32) :: i, m, nxfrm
    type(stft_result) :: y
    complex(real64), allocatable, dimension(:,:) :: ans

    rst = .true.
    m = win%size
    nxfrm = compute_transform_length(m)
    y = stft(win, x)

    if (size(y%stft, 1) /= nxfrm .or. size(y%stft, 2) /= nk .or. &
        size(y%offsets) /= nk) &
    then
        rst = .false.
        print '(A)', "TEST FAILED: " // id // " (output size)"
        return
    end if

    ! The segments must span the signal from the first to the last point
    if (y%offsets(1) /= 1 .or. y%offsets(nk) + m - 1 /= size(x)) then
        rst = .false.
        print '(A)', "TEST FAILED: " // id // " (offsets)"
    end if

    ans = reference_stft(win, x, nk)
    do i = 1, nk
        if (.not.assert_transform(ans(:,i), y%stft(:,i), id)) rst = .false.
    end do
end function

! ------------------------------------------------------------------------------
function test_stft_even_window() result(rst)
    !! Tests stft using an even-sized window.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: m = 8
    integer(int32), parameter :: nx = 32
    integer(int32), parameter :: nxfrm = m / 2 + 1
    integer(int32), parameter :: nk = (nx - 1) / nxfrm

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(nx), a, ans(nxfrm)
    type(rectangular_window) :: rect
    type(hamming_window) :: hamm
    type(stft_result) :: y

    ! Initialization
    rst = .true.
    rect%size = m
    hamm%size = m
    call random_number(x)

    ! Test 1 - A rectangular window against an explicitly computed reference
    if (.not.check_stft(rect, x, nk, "test_stft_even_window 1-1")) rst = .false.

    ! Test 2 - A non-trivial window must be normalized such that the result is
    ! independent of the window's amplitude
    if (.not.check_stft(hamm, x, nk, "test_stft_even_window 1-2")) rst = .false.

    ! Test 3 - A bin-aligned sinusoid must report its actual amplitude
    a = 2.5d0
    do i = 1, nx
        x(i) = a * cos(2.0d0 * pi * 2.0d0 * (i - 1.0d0) / m)
    end do
    y = stft(rect, x)
    ans = 0.0d0
    ans(3) = a
    do i = 1, nk
        if (.not.assert(ans, abs(y%stft(:,i)), tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_stft_even_window 1-3"
        end if
    end do

    ! Test 4 - The DC and Nyquist terms must be purely real
    x = 1.0d0
    y = stft(rect, x)
    ans = 0.0d0
    ans(1) = 1.0d0
    do i = 1, nk
        if (.not.assert(ans, abs(y%stft(:,i)), tol) .or. &
            .not.assert(0.0d0, aimag(y%stft(1,i)), tol) .or. &
            .not.assert(0.0d0, aimag(y%stft(nxfrm,i)), tol)) &
        then
            rst = .false.
            print '(A)', "TEST FAILED: test_stft_even_window 1-4"
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_stft_odd_window() result(rst)
    !! Tests stft using an odd-sized window.  There is no Nyquist term in this
    !! case, so the highest-frequency bin is treated as any other bin.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: m = 7
    integer(int32), parameter :: nx = 30
    integer(int32), parameter :: nxfrm = (m + 1) / 2
    integer(int32), parameter :: nk = 7

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(nx), a, ans(nxfrm)
    type(rectangular_window) :: rect
    type(hann_window) :: hann
    type(stft_result) :: y

    ! Initialization
    rst = .true.
    rect%size = m
    hann%size = m
    call random_number(x)

    ! Test 1 - A rectangular window against an explicitly computed reference
    if (.not.check_stft(rect, x, nk, "test_stft_odd_window 1-1")) rst = .false.

    ! Test 2 - A tapered window
    if (.not.check_stft(hann, x, nk, "test_stft_odd_window 1-2")) rst = .false.

    ! Test 3 - A sinusoid on the highest available bin must report its actual
    ! amplitude
    a = 1.75d0
    do i = 1, nx
        x(i) = a * cos(2.0d0 * pi * (nxfrm - 1.0d0) * (i - 1.0d0) / m)
    end do
    y = stft(rect, x)
    ans = 0.0d0
    ans(nxfrm) = a
    do i = 1, nk
        if (.not.assert(ans, abs(y%stft(:,i)), tol)) then
            rst = .false.
            print '(A)', "TEST FAILED: test_stft_odd_window 1-3"
        end if
    end do

    ! Test 4 - A constant-valued signal yields a purely real, unity DC term
    x = 1.0d0
    y = stft(rect, x)
    ans = 0.0d0
    ans(1) = 1.0d0
    do i = 1, nk
        if (.not.assert(ans, abs(y%stft(:,i)), tol) .or. &
            .not.assert(0.0d0, aimag(y%stft(1,i)), tol)) &
        then
            rst = .false.
            print '(A)', "TEST FAILED: test_stft_odd_window 1-4"
        end if
    end do
end function

! ------------------------------------------------------------------------------
function test_stft_odd_input_size_guard() result(rst)
    !! Tests that stft returns an empty result for invalid window sizes.
    logical :: rst

    ! Local Variables
    real(real64) :: x(16)
    type(rectangular_window) :: win
    type(stft_result) :: y

    ! Initialization
    rst = .true.
    x = 1.0d0

    ! Test 1 - A zero-sized window
    win%size = 0
    y = stft(win, x)
    if (allocated(y%stft) .or. allocated(y%offsets)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_stft_odd_input_size_guard 1-1"
    end if

    ! Test 2 - A window longer than the signal
    win%size = size(x) + 1
    y = stft(win, x)
    if (allocated(y%stft) .or. allocated(y%offsets)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_stft_odd_input_size_guard 1-2"
    end if
end function

! ------------------------------------------------------------------------------
function test_irfft_odd_input_size() result(rst)
    !! Tests irfft for an input array of odd length.  An odd-length input of
    !! length M is interpreted as the positive half of the transform of an
    !! even-length signal of length N = 2 * (M - 1).
    logical :: rst

    ! Parameters
    integer(int32), parameter :: m = 5     ! odd-length input
    integer(int32), parameter :: n = 2 * (m - 1)

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), ans(n), a, b
    real(real64), allocatable, dimension(:) :: y
    complex(real64) :: xfrm(m)

    ! Initialization
    rst = .true.

    ! Test 1 - Verify the output length
    xfrm = (0.0d0, 0.0d0)
    y = irfft(xfrm)
    if (size(y) /= n) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_odd_input_size 1-1"
        return
    end if

    ! Test 2 - A DC-only transform must produce a constant-valued signal
    a = 2.5d0
    xfrm = (0.0d0, 0.0d0)
    xfrm(1) = cmplx(a, 0.0d0, real64)
    y = irfft(xfrm)
    ans = a
    if (.not.assert(ans, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_odd_input_size 1-2"
    end if

    ! Test 3 - A single non-DC, non-Nyquist term.  As FFTPACK's backward
    ! transform is unnormalized, the expected result is
    ! y(k) = 2 * (Re(X) * cos(2 * pi * j * k / n) - Im(X) * sin(2 * pi * j * k / n))
    ! for the j-th bin.
    a = 1.5d0
    b = -0.75d0
    xfrm = (0.0d0, 0.0d0)
    xfrm(2) = cmplx(a, b, real64)
    y = irfft(xfrm)
    do i = 1, n
        ans(i) = 2.0d0 * (a * cos(2.0d0 * pi * (i - 1.0d0) / n) - &
            b * sin(2.0d0 * pi * (i - 1.0d0) / n))
    end do
    if (.not.assert(ans, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_odd_input_size 1-3"
    end if

    ! Test 4 - The Nyquist term is only present for even-length signals, and
    ! must produce an alternating +/- sequence.
    a = 3.0d0
    xfrm = (0.0d0, 0.0d0)
    xfrm(m) = cmplx(a, 0.0d0, real64)
    y = irfft(xfrm)
    do i = 1, n
        ans(i) = a * (-1.0d0)**(i - 1)
    end do
    if (.not.assert(ans, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_odd_input_size 1-4"
    end if

    ! Test 5 - Round-trip through rfft.  The inverse transform is unnormalized
    ! such that irfft(rfft(x)) = n * x.
    call random_number(x)
    y = irfft(rfft(x))
    if (.not.assert(n * x, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_odd_input_size 1-5"
    end if
end function

! ------------------------------------------------------------------------------
function test_irfft_even_input_size() result(rst)
    !! Tests irfft for an input array of even length.  An even-length input of
    !! length M is interpreted as the positive half of the transform of an
    !! odd-length signal of length N = 2 * M - 1.
    logical :: rst

    ! Parameters
    integer(int32), parameter :: m = 4     ! even-length input
    integer(int32), parameter :: n = 2 * m - 1

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n), ans(n), a, b
    real(real64), allocatable, dimension(:) :: y
    complex(real64) :: xfrm(m)

    ! Initialization
    rst = .true.

    ! Test 1 - Verify the output length
    xfrm = (0.0d0, 0.0d0)
    y = irfft(xfrm)
    if (size(y) /= n) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_even_input_size 1-1"
        return
    end if

    ! Test 2 - A DC-only transform must produce a constant-valued signal
    a = -1.25d0
    xfrm = (0.0d0, 0.0d0)
    xfrm(1) = cmplx(a, 0.0d0, real64)
    y = irfft(xfrm)
    ans = a
    if (.not.assert(ans, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_even_input_size 1-2"
    end if

    ! Test 3 - A single non-DC term
    a = 0.8d0
    b = 2.2d0
    xfrm = (0.0d0, 0.0d0)
    xfrm(3) = cmplx(a, b, real64)
    y = irfft(xfrm)
    do i = 1, n
        ans(i) = 2.0d0 * (a * cos(4.0d0 * pi * (i - 1.0d0) / n) - &
            b * sin(4.0d0 * pi * (i - 1.0d0) / n))
    end do
    if (.not.assert(ans, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_even_input_size 1-3"
    end if

    ! Test 4 - The highest-frequency term of an odd-length signal retains both
    ! real and imaginary components as there is no Nyquist bin.
    a = 1.0d0
    b = -2.0d0
    xfrm = (0.0d0, 0.0d0)
    xfrm(m) = cmplx(a, b, real64)
    y = irfft(xfrm)
    do i = 1, n
        ans(i) = 2.0d0 * (a * cos(2.0d0 * pi * (m - 1.0d0) * (i - 1.0d0) / n) - &
            b * sin(2.0d0 * pi * (m - 1.0d0) * (i - 1.0d0) / n))
    end do
    if (.not.assert(ans, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_even_input_size 1-4"
    end if

    ! Test 5 - Round-trip through rfft
    call random_number(x)
    y = irfft(rfft(x))
    if (.not.assert(n * x, y, tol)) then
        rst = .false.
        print '(A)', "TEST FAILED: test_irfft_even_input_size 1-5"
    end if
end function

! ------------------------------------------------------------------------------
end module
