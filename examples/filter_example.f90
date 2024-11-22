program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: winsize = 10
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: f1 = 5.0d0
    real(real64), parameter :: f2 = 1.5d1
    real(real64), parameter :: cutoff = 2.0d1
    real(real64), parameter :: alpha = 3.0d0
    
    ! Local Variables
    integer(int32) :: i, k
    real(real64) :: fs, t(npts), x(npts), y(npts), b(winsize), a(1), ys(npts)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2, d3
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Build the signal
    t = linspace(0.0d0, 1.0d0, npts)
    fs = 1.0d0 / (t(2) - t(1))
    call random_number(x)
    x = 0.25d0 * (x - 0.5d0) + sin(2.0d0 * pi * f1 * t) + &
        0.5 * sin(2.0d0 * pi * f2 * t)

    ! Define the filter coefficients.  An averaging-type filter of window size 
    ! winsize is defined.  This is an FIR type filter.
    b = 1.0d0 / winsize ! all values in the array are the same
    a = 1.0d0

    ! Apply the filter
    y = filter(b, a, x)

    ! Apply a sinc filter
    ys = sinc_filter(cutoff, fs, x)

    ! Plot the results
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call lgnd%set_is_visible(.true.)

    call d1%define_data(t, x)
    call d1%set_name("Original")
    call plt%push(d1)

    call d2%define_data(t, y)
    call d2%set_name("Smoothed")
    call plt%push(d2)

    call d3%define_data(t, ys)
    call d3%set_name("Sinc Filtered")
    call plt%push(d3)

    call plt%draw()
end program