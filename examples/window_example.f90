program example
    use iso_fortran_env
    use ieee_arithmetic
    use spectrum
    use fftpack
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: winsize = 64
    integer(int32), parameter :: fftsize = 1024

    ! Local Variables
    integer(int32) :: i, j
    real(real64), dimension(winsize) :: samples, hw, hmw, bw, ww
    real(real64) :: freq(fftsize)
    type(hann_window) :: hannwin
    type(hamming_window) :: hamwin
    type(blackman_harris_window) :: bwin
    type(welch_window) :: wwin
    complex(real64), dimension(fftsize) :: hwfft, hmwfft, bwfft, wwfft

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2, pd3, pd4
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd
    class(terminal), pointer :: term

    ! Build the windows
    hannwin%size = winsize
    hamwin%size = winsize
    bwin%size = winsize
    wwin%size = winsize
    do i = 1, winsize
        j = i - 1
        samples(i) = real(j, real64)
        hw(i) = hannwin%evaluate(j)
        hmw(i) = hamwin%evaluate(j)
        bw(i) = bwin%evaluate(j)
        ww(i) = wwin%evaluate(j)
    end do

    ! Compute FFT's of the window functions
    hwfft = fftshift( fft(cmplx(hw, 0.0d0, real64), fftsize) ) / &
        (0.5d0 * winsize)
    hmwfft = fftshift( fft(cmplx(hmw, 0.0d0, real64), fftsize) ) / &
        (0.5d0 * winsize)
    bwfft = fftshift( fft(cmplx(bw, 0.0d0, real64), fftsize) ) / &
        (0.5d0 * winsize)
    wwfft = fftshift( fft(cmplx(ww, 0.0d0, real64), fftsize) ) / &
        (0.5d0 * winsize)

    ! Plot the windows
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()
    term => plt%get_terminal()

    call xAxis%set_title("Sample")
    call yAxis%set_title("Amplitude")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, real(winsize, real64))
    call lgnd%set_is_visible(.true.)
    call lgnd%set_draw_border(.false.)
    call lgnd%set_draw_inside_axes(.false.)
    call term%set_window_width(800)

    call pd1%define_data(samples, hw)
    call pd1%set_name("Hann (Hanning)")
    call pd1%set_line_width(2.0)
    call plt%push(pd1)

    call pd2%define_data(samples, hmw)
    call pd2%set_name("Hamming")
    call pd2%set_line_width(2.0)
    call pd2%set_line_style(LINE_DASHED)
    call plt%push(pd2)

    call pd3%define_data(samples, bw)
    call pd3%set_name("Blackman-Harris")
    call pd3%set_line_width(2.0)
    call pd3%set_line_style(LINE_DASH_DOTTED)
    call plt%push(pd3)

    call pd4%define_data(samples, ww)
    call pd4%set_name("Welch")
    call pd4%set_line_width(2.0)
    call pd4%set_line_style(LINE_DOTTED)
    call plt%push(pd4)
    
    call plt%draw()
    call plt%clear_all()

    ! Plot the window FFT's
    freq = linspace(-0.5d0, 0.5d0, fftsize)
    call xAxis%set_title("Normalized Frequency")
    call xAxis%set_limits(-0.5d0, 0.5d0)
    call yAxis%set_title("Magnitude [dB]")
    call yAxis%set_autoscale(.false.)
    call yAxis%set_limits(-1.2d2, 0.0d0)

    call pd1%define_data(freq, to_db(hwfft))
    call pd1%set_line_width(1.0)
    call plt%push(pd1)

    call pd2%define_data(freq, to_db(hmwfft))
    call pd2%set_line_width(1.0)
    call plt%push(pd2)

    call pd3%define_data(freq, to_db(bwfft))
    call pd3%set_line_width(1.0)
    call plt%push(pd3)

    call pd4%define_data(freq, to_db(wwfft))
    call pd4%set_line_width(1.0)
    call plt%push(pd4)

    call plt%draw()

contains
    pure function to_db(x) result(rst)
        complex(real64), intent(in) :: x(:)
        real(real64) :: rst(size(x))
        real(real64) :: absx(size(x)), mx
        absx = abs(x)
        mx = maxval(absx)
        rst = 2.0d1 * log10(absx / mx)
        where (.not.ieee_is_finite(rst)) rst = ieee_value(mx, ieee_quiet_nan)
    end function
end program