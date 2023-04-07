! Model: SDOF mechanical vibrating system:
! x" + 2 z wn x' + wn**2 x = 2 z wn y' + wn**2 y
!
! Input: Delta Function (Y(s) = 1 / s)
! Output:
! x(t) = 1 + exp(-wn * z * t) * [(z / sqrt(1 - z**2)) * sin(a) - cos(a)] 
! where,
! a = t * wn * sqrt(1 - z**2)
program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 10000
    integer(int32), parameter :: winsize = 512
    integer(int32), parameter :: nxfrm = winsize / 2 + 1
    real(real64), parameter :: sample_rate = 1024.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: zeta = 0.05d0
    real(real64), parameter :: nat_freq_hz = 7.5d1
    real(real64), parameter :: nat_freq = 2.0d0 * pi * nat_freq_hz

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dt, zarg, arg, t(npts), x(npts), y(npts)
    real(real64) :: df, freq(nxfrm), mag(nxfrm), phase(nxfrm)
    complex(real64) :: tf(nxfrm)
    type(rectangular_window) :: win

    ! Plot Variables
    type(multiplot) :: mplt
    type(plot_2d) :: plt, plt1, plt2
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Define the time vector, input signal, and the resulting output signal
    dt = 1.0d0 / sample_rate
    t = (/ (dt * i, i = 0, npts - 1) /)
    y(1) = 1.0d0
    y(2:) = 0.0d0
    zarg = sqrt(1.0d0 - zeta**2)
    arg = nat_freq * zarg
    x = 1.0d0 + exp(-nat_freq * zeta * t) * &
        ((zeta / zarg) * sin(arg * t) - cos(arg * t))

    ! Compute the transfer function
    win%size = winsize
    tf = siso_transfer_function(win, y, x)

    ! Compute the frequency, magnitude, and phase
    df = frequency_bin_width(sample_rate, winsize)
    freq = (/ (df * i, i = 0, nxfrm - 1) /)
    mag = 2.0d1 * log10(abs(tf))    ! Convert to dB
    phase = atan2(aimag(tf), real(tf))

    ! Unwrap the phase
    call unwrap(phase, pi / 2.0d0)

    ! Set up the plot objects
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    ! Plot the time signal
    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")

    call pd%define_data(t, x)
    call pd%set_line_width(2.0)
    call plt%Push(pd)
    call plt%draw()

    ! Plot the transfer function
    call mplt%initialize(2, 1)
    call plt1%initialize()
    call plt2%initialize()

    ! Plot the magnitude component
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("|X / Y| [dB]")

    call pd%define_data(freq(:nxfrm-1), mag(:nxfrm-1))
    call plt1%push(pd)

    ! Plot the phase component
    xAxis => plt2%get_x_axis()
    yAxis => plt2%get_y_axis()
    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("{/Symbol f} [deg]")

    call pd%define_data(freq, phase * 1.8d2 / pi)
    call plt2%push(pd)

    call mplt%set(1, 1, plt1)
    call mplt%set(2, 1, plt2)
    call mplt%draw()
end program