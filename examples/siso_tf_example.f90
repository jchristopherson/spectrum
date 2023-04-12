program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    use csv_module
    implicit none

    ! Parameters
    integer(int32), parameter :: winsize = 1024
    integer(int32), parameter :: nxfrm = winsize / 2 + 1
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    complex(real64), parameter :: j = (0.0d0, 1.0d0)
    real(real64), parameter :: zeta = 1.0d-2
    real(real64), parameter :: fn = 2.5d2
    real(real64), parameter :: wn = 2.0d0 * pi * fn

    ! Local Variables
    logical :: ok
    integer(int32) :: i
    real(real64) :: dt, fs, df, freq(nxfrm)
    real(real64), allocatable, dimension(:) :: t, x, y, mag, phase, maga, pa
    complex(real64), allocatable, dimension(:) :: tf, s, tfa
    type(hamming_window) :: win
    type(csv_file) :: file

    ! Plot Variables
    type(multiplot) :: mplt
    type(plot_2d) :: plt, plt1, plt2
    type(plot_data_2d) :: pd, pda
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Import time history data
    call file%read("files/chirp.csv", status_ok = ok)
    call file%get(1, t, ok)
    call file%get(2, y, ok)
    call file%get(3, x, ok)

    ! Determine the sample rate
    dt = t(2) - t(1)
    fs = 1.0d0 / dt
    print 100, "Sample Rate: ", fs, " Hz"

    ! Compute the transfer function
    win%size = winsize
    tf = siso_transfer_function(win, y, x)

    ! Compute the frequency, magnitude, and phase
    df = frequency_bin_width(fs, winsize)
    freq = (/ (df * i, i = 0, nxfrm - 1) /)
    mag = 2.0d1 * log10(abs(tf))    ! Convert to dB
    phase = atan2(aimag(tf), real(tf))
    call unwrap(phase, pi / 2.0d0)

    ! Compare to the analytical solution
    s = j * (2.0d0 * pi * freq)
    tfa = (2.0d0 * zeta * wn * s + wn**2) / &
        (s**2 + 2.0d0 * zeta * wn * s + wn**2)
    maga = 2.0d1 * log10(abs(tfa))
    pa = atan2(aimag(tfa), real(tfa))
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
    call plt%push(pd)
    call plt%draw()

    ! Plot the transfer function
    call mplt%initialize(2, 1)
    call plt1%initialize()
    call plt2%initialize()

    ! Plot the magnitude component
    xAxis => plt1%get_x_axis()
    yAxis => plt1%get_y_axis()
    lgnd => plt1%get_legend()
    call xAxis%set_title("f [Hz]")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 5.0d2)
    call yAxis%set_title("|X / Y| [dB]")
    call lgnd%set_is_visible(.true.)

    call pd%define_data(freq, mag)
    call pd%set_name("Computed")
    call plt1%push(pd)

    call pda%define_data(freq, maga)
    call pda%set_name("Analytical")
    call pda%set_line_style(LINE_DASHED)
    call pda%set_line_width(2.5)
    call plt1%push(pda)

    ! Plot the phase component
    xAxis => plt2%get_x_axis()
    yAxis => plt2%get_y_axis()
    call xAxis%set_title("f [Hz]")
    call xAxis%set_autoscale(.false.)
    call xAxis%set_limits(0.0d0, 5.0d2)
    call yAxis%set_title("{/Symbol f} [deg]")

    call pd%define_data(freq, phase * 1.8d2 / pi)
    call plt2%push(pd)

    call pda%define_data(freq, pa * 1.8d2 / pi)
    call plt2%push(pda)

    call mplt%set(1, 1, plt1)
    call mplt%set(2, 1, plt2)
    call mplt%draw()

    ! Formatting
100 format(A, F6.1, A)
end program