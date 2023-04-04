program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 20000
    integer(int32), parameter :: winsize = 1024
    real(real64), parameter :: sample_rate = 2.048d3
    real(real64), parameter :: freq1 = 5.0d1
    real(real64), parameter :: freq2 = 2.5d2
    real(real64), parameter :: phase1 = 0.0d0
    real(real64), parameter :: phase2 = 4.5d1
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: amp1 = 1.5d0
    real(real64), parameter :: amp2 = 7.5d-1

    ! Local Variables
    integer(int32) :: i, nxfrm
    real(real64) :: df, dt, t(npts), x(npts), y(npts)
    real(real64), allocatable, dimension(:) :: xfrm, freq
    type(hann_window) :: win

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Build the signal - assume units of V
    dt = 1.0d0 / sample_rate
    t = (/ (dt * i, i = 0, npts - 1) /)
    x = amp1 * sin(2.0d0 * pi * freq1 * t - pi * phase1 / 1.8d2) + &
        0.25d0 * amp1 * sin(2.0d0 * pi * freq2 * t)
    y = amp2 * sin(2.0d0 * pi * freq2 * t - pi * phase2 / 1.8d2) + &
        0.5d0 * amp1 * sin(2.0d0 * pi * freq1 * t)

    ! Define the window
    win%size = winsize

    ! Compute the
    xfrm = csd(win, x, y)

    ! BUild a corresponding array of frequency values
    df = frequency_bin_width(sample_rate, winsize)
    allocate(freq(size(xfrm)))
    freq = (/ (df * i, i = 0, size(xfrm) - 1) /)

    ! Plot the spectrum
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("X * Y [V^{2}]")

    call xAxis%set_is_log_scaled(.true.)
    call xAxis%set_use_default_tic_label_format(.false.)
    call xAxis%set_tic_label_format("%0.0e")

    call pd%define_data(freq, xfrm)
    call pd%set_line_width(2.0)
    call plt%push(pd)
    call plt%draw()
end program