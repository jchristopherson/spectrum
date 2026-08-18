module spectrum_filter_design
    use iso_fortran_env
    use spectrum_filter, only : LOW_PASS_FILTER, HIGH_PASS_FILTER, &
        BAND_PASS_FILTER, BAND_STOP_FILTER
    use spectrum_windows, only : window, hamming_window
    implicit none
    private
    public :: butterworth_filter_order
    public :: design_fir_filter
    public :: design_iir_filter
    public :: filter_frequency_response

contains
! ------------------------------------------------------------------------------
pure function butterworth_filter_order(fp, fs, fstop, apass, astop, ftype) result(n)
    !! Computes the minimum order of a Butterworth low-pass or high-pass filter.
    !!
    !! Frequencies are specified in Hz and attenuation values are specified in dB.
    !! The frequency prewarping matches [[design_iir_filter]]. An order of zero
    !! indicates an invalid specification.
    real(real64), intent(in) :: fp
        !! The passband edge frequency, in Hz.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    real(real64), intent(in) :: fstop
        !! The stopband edge frequency, in Hz.
    real(real64), intent(in) :: apass
        !! The maximum passband attenuation, in dB.
    real(real64), intent(in) :: astop
        !! The minimum stopband attenuation, in dB.
    integer(int32), intent(in), optional :: ftype
        !! LOW_PASS_FILTER (default) or HIGH_PASS_FILTER.
    integer(int32) :: n
        !! The minimum filter order.

    integer(int32) :: filter_type
    real(real64) :: pi, wp, ws, numerator, denominator, order_value

    n = 0_int32
    filter_type = LOW_PASS_FILTER
    if (present(ftype)) filter_type = ftype
    if (fs <= 0.0d0 .or. fp <= 0.0d0 .or. fstop <= 0.0d0 .or. &
        fp >= 0.5d0 * fs .or. fstop >= 0.5d0 * fs .or. &
        apass <= 0.0d0 .or. astop <= apass .or. &
        (filter_type /= LOW_PASS_FILTER .and. filter_type /= HIGH_PASS_FILTER)) return
    if ((filter_type == LOW_PASS_FILTER .and. fstop <= fp) .or. &
        (filter_type == HIGH_PASS_FILTER .and. fstop >= fp)) return

    pi = acos(-1.0d0)
    wp = 2.0d0 * fs * tan(pi * fp / fs)
    ws = 2.0d0 * fs * tan(pi * fstop / fs)
    if (filter_type == HIGH_PASS_FILTER) then
        ws = wp
        wp = 2.0d0 * fs * tan(pi * fstop / fs)
    end if

    numerator = log10(10.0d0 ** (astop / 10.0d0) - 1.0d0)
    denominator = log10(10.0d0 ** (apass / 10.0d0) - 1.0d0)
    order_value = (numerator - denominator) / (2.0d0 * log10(ws / wp))
    if (order_value > 0.0d0) n = int(ceiling(order_value), int32)
end function

! ------------------------------------------------------------------------------
pure subroutine design_fir_filter(n, fc, fs, b, a, fc2, ftype, win)
    !! Designs a windowed, linear-phase FIR filter.
    !!
    !! The returned coefficients use the same convention as filter: b contains
    !! the numerator coefficients in increasing delay order and a is [1.0].
    integer(int32), intent(in) :: n
        !! The requested number of taps. Even values are increased by one.
    real(real64), intent(in) :: fc
        !! The first cutoff frequency, in Hz.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    real(real64), allocatable, intent(out) :: b(:)
        !! The numerator coefficients of the designed filter.
    real(real64), allocatable, intent(out) :: a(:)
        !! The denominator coefficients. This is always [1.0].
    real(real64), intent(in), optional :: fc2
        !! The second cutoff frequency for band-pass and band-stop filters.
    integer(int32), intent(in), optional :: ftype
        !! The filter type. The default is LOW_PASS_FILTER.
    class(window), intent(inout), optional, target :: win
        !! The window function. The default is a Hamming window.

    integer(int32) :: i, ntaps, filter_type, midpoint
    real(real64) :: cutoff2, offset, sum_b
    type(hamming_window), target :: default_window
    class(window), pointer :: selected_window

    ntaps = n
    if (mod(ntaps, 2) == 0) ntaps = ntaps + 1
    filter_type = LOW_PASS_FILTER
    if (present(ftype)) filter_type = ftype
    cutoff2 = 0.0d0
    if (present(fc2)) cutoff2 = fc2

    if (ntaps < 1 .or. fs <= 0.0d0 .or. fc <= 0.0d0 .or. &
        fc >= 0.5d0 * fs) return
    if (filter_type < LOW_PASS_FILTER .or. filter_type > BAND_STOP_FILTER) return
    if (filter_type == BAND_PASS_FILTER .or. filter_type == BAND_STOP_FILTER) then
        if (.not. present(fc2) .or. cutoff2 <= fc .or. cutoff2 >= 0.5d0 * fs) return
    end if

    allocate(b(ntaps), a(1))
    a = 1.0d0
    midpoint = (ntaps - 1) / 2
    if (present(win)) then
        selected_window => win
    else
        selected_window => default_window
    end if
    selected_window%size = ntaps

    do i = 1, ntaps
        offset = real(i - 1 - midpoint, real64)
        b(i) = windowed_kernel(fc, fs, offset, ntaps, selected_window)
    end do

    select case (filter_type)
    case (HIGH_PASS_FILTER)
        b = -b
        b(midpoint + 1) = b(midpoint + 1) + 1.0d0
    case (BAND_PASS_FILTER, BAND_STOP_FILTER)
        do i = 1, ntaps
            offset = real(i - 1 - midpoint, real64)
            b(i) = b(i) - windowed_kernel(cutoff2, fs, offset, ntaps, selected_window)
        end do
        if (filter_type == BAND_STOP_FILTER) then
            b(midpoint + 1) = b(midpoint + 1) + 1.0d0
        else
            b = -b
        end if
    end select

    sum_b = sum(b)
    if (filter_type == LOW_PASS_FILTER .or. filter_type == BAND_STOP_FILTER) then
        if (abs(sum_b) > epsilon(1.0d0)) b = b / sum_b
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine design_iir_filter(n, fc, fs, b, a, ftype)
    !! Designs an nth-order Butterworth IIR low-pass or high-pass filter.
    !!
    !! The returned coefficients are normalized such that a(1) is one and use
    !! the increasing-delay convention required by filter.
    integer(int32), intent(in) :: n
        !! The Butterworth filter order.
    real(real64), intent(in) :: fc
        !! The cutoff frequency, in Hz.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    real(real64), allocatable, intent(out) :: b(:)
        !! The numerator coefficients of the designed filter.
    real(real64), allocatable, intent(out) :: a(:)
        !! The denominator coefficients of the designed filter.
    integer(int32), intent(in), optional :: ftype
        !! LOW_PASS_FILTER (default) or HIGH_PASS_FILTER.

    integer(int32) :: i, j, filter_type
    real(real64) :: analog_cutoff, gain, pi, phase, numerator, denominator_sum
    complex(real64) :: pole, root
    complex(real64), allocatable :: denominator(:)

    filter_type = LOW_PASS_FILTER
    if (present(ftype)) filter_type = ftype
    if (n < 1 .or. fs <= 0.0d0 .or. fc <= 0.0d0 .or. &
        fc >= 0.5d0 * fs) return
    if (filter_type /= LOW_PASS_FILTER .and. &
        filter_type /= HIGH_PASS_FILTER) return

    pi = acos(-1.0d0)
    analog_cutoff = 2.0d0 * fs * tan(pi * fc / fs)
    allocate(b(n + 1), a(n + 1), denominator(n + 1))
    b = 0.0d0
    denominator = (0.0d0, 0.0d0)
    b(1) = 1.0d0
    denominator(1) = (1.0d0, 0.0d0)

    do i = 1, n
        phase = pi * real(2 * i + n - 1, real64) / real(2 * n, real64)
        pole = analog_cutoff * exp(cmplx(0.0d0, phase, real64))
        root = (2.0d0 * fs + pole) / (2.0d0 * fs - pole)
        do j = i + 1, 2, -1
            denominator(j) = denominator(j) - root * denominator(j - 1)
            if (filter_type == LOW_PASS_FILTER) then
                b(j) = b(j) + b(j - 1)
            else
                b(j) = b(j) - b(j - 1)
            end if
        end do
    end do
    a = real(denominator, real64)

    if (filter_type == LOW_PASS_FILTER) then
        gain = sum(a) / sum(b)
    else
        numerator = 0.0d0
        denominator_sum = 0.0d0
        do i = 1, n + 1
            numerator = numerator + (-1.0d0) ** (i - 1) * a(i)
            denominator_sum = denominator_sum + (-1.0d0) ** (i - 1) * b(i)
        end do
        gain = numerator / denominator_sum
    end if
    b = gain * b
end subroutine

! ------------------------------------------------------------------------------
pure function filter_frequency_response(b, a, f, fs) result(rst)
    !! Computes the complex frequency response of a rational filter.
    !!
    !! The coefficient ordering is the same as [[filter]]: \(b(i)\) and \(a(i)\) 
    !! where \(z = \exp{\frac{2 \pi i f}{f_s}}\). The returned values are 
    !! therefore \(H(f) = \frac{B(z)}{A(z)}\) at each requested frequency.
    real(real64), intent(in) :: b(:)
        !! The numerator coefficients of the rational transfer function.
    real(real64), intent(in) :: a(:)
        !! The denominator coefficients of the rational transfer function.
    real(real64), intent(in) :: f(:)
        !! The frequencies at which to evaluate the response, in Hz.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    complex(real64), allocatable :: rst(:)
        !! The complex frequency response at each frequency in f.

    integer(int32) :: i, j, n
    real(real64) :: pi, phase
    complex(real64) :: z, numerator, denominator

    n = size(f)
    if (size(a) < 1 .or. fs <= 0.0d0) return
    allocate(rst(n), source = (0.0d0, 0.0d0))
    pi = acos(-1.0d0)

    do i = 1, n
        phase = 2.0d0 * pi * f(i) / fs
        z = exp(cmplx(0.0d0, -phase, real64))
        numerator = (0.0d0, 0.0d0)
        denominator = (0.0d0, 0.0d0)
        do j = size(b), 1, -1
            numerator = numerator * z + b(j)
        end do
        do j = size(a), 1, -1
            denominator = denominator * z + a(j)
        end do
        if (abs(denominator) > tiny(1.0d0)) then
            rst(i) = numerator / denominator
        end if
    end do
end function

! ------------------------------------------------------------------------------
pure function lowpass_kernel(fc, fs, offset) result(value)
    real(real64), intent(in) :: fc, fs, offset
    real(real64) :: value, pi, argument

    pi = acos(-1.0d0)
    argument = 2.0d0 * fc / fs * offset
    if (abs(argument) < epsilon(1.0d0)) then
        value = 2.0d0 * fc / fs
    else
        value = 2.0d0 * fc / fs * sin(pi * argument) / (pi * argument)
    end if
end function

! ------------------------------------------------------------------------------
pure function windowed_kernel(fc, fs, offset, ntaps, win) result(value)
    real(real64), intent(in) :: fc, fs, offset
    integer(int32), intent(in) :: ntaps
    class(window), intent(in) :: win
    real(real64) :: value, pi, index, window

    pi = acos(-1.0d0)
    value = lowpass_kernel(fc, fs, offset)
    if (ntaps > 1) then
        index = offset + 0.5d0 * real(ntaps - 1, real64)
        window = win%evaluate(int(index, int32))
        value = window * value
    end if
end function

! ------------------------------------------------------------------------------
end module