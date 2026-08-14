module spectrum_filter_design
    use iso_fortran_env
    use spectrum_filter, only : LOW_PASS_FILTER, HIGH_PASS_FILTER, &
        BAND_PASS_FILTER, BAND_STOP_FILTER
    implicit none
    private
    public :: design_fir_filter
    public :: filter_frequency_response

contains
! ------------------------------------------------------------------------------
pure subroutine design_fir_filter(n, fc, fs, b, a, fc2, ftype)
    !! Designs a Hamming-windowed, linear-phase FIR filter.
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

    integer(int32) :: i, ntaps, filter_type, midpoint
    real(real64) :: cutoff2, offset, sum_b

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

    do i = 1, ntaps
        offset = real(i - 1 - midpoint, real64)
        b(i) = windowed_kernel(fc, fs, offset, ntaps)
    end do

    select case (filter_type)
    case (HIGH_PASS_FILTER)
        b = -b
        b(midpoint + 1) = b(midpoint + 1) + 1.0d0
    case (BAND_PASS_FILTER, BAND_STOP_FILTER)
        do i = 1, ntaps
            offset = real(i - 1 - midpoint, real64)
            b(i) = b(i) - windowed_kernel(cutoff2, fs, offset, ntaps)
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
pure function filter_frequency_response(b, a, f, fs) result(rst)
    !! Computes the complex frequency response of a rational filter.
    !!
    !! The coefficient ordering is the same as filter: b(i) and a(i) multiply
    !! z**(-(i - 1)), where z = exp(2*pi*i*f/fs). The returned values are
    !! therefore H(f) = B(z) / A(z) at each requested frequency.
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
pure function windowed_kernel(fc, fs, offset, ntaps) result(value)
    real(real64), intent(in) :: fc, fs, offset
    integer(int32), intent(in) :: ntaps
    real(real64) :: value, pi, index, window

    pi = acos(-1.0d0)
    value = lowpass_kernel(fc, fs, offset)
    if (ntaps > 1) then
        index = offset + 0.5d0 * real(ntaps - 1, real64)
        window = 0.54d0 - 0.46d0 * cos(2.0d0 * pi * index / real(ntaps - 1, real64))
        value = window * value
    end if
end function

! ------------------------------------------------------------------------------
end module