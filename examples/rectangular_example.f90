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
    real(real64) :: freq(fftsize), samples(winsize), w(winsize)
    complex(real64) :: wfft(fftsize)
    type(rectangular_window) :: win

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Build the window
    win%size = winsize
    do i = 1, winsize
        j = i - 1
        samples(i) = real(j, real64)
        w(i) = win%evaluate(j)
    end do

    ! Compute the FFT
    wfft = fftshift( fft(cmplx(w, 0.0d0, real64), fftsize) ) / (0.5d0 * winsize)

    ! Plot the window
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("Sample")
    call yAxis%set_title("Amplitude")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, real(winsize, real64))

    call pd%define_data(samples, w)
    call pd%set_line_width(2.0)
    call plt%push(pd)

    call plt%draw()
    call plt%clear_all()

    ! Plot the FFT
    freq = linspace(-0.5d0, 0.5d0, fftsize)
    call xAxis%set_title("Normalized Frequency")
    call xAxis%set_limits(-0.5d0, 0.5d0)
    call yAxis%set_title("Magnitude [dB]")
    call yAxis%set_autoscale(.false.)
    call yAxis%set_limits(-1.2d2, 0.0d0)

    call pd%define_data(freq, to_db(wfft))
    call plt%push(pd)

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