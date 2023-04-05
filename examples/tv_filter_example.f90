program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 10000
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: f1 = 5.0d0
    real(real64), parameter :: f2 = 1.5d1
    real(real64), parameter :: alpha = 3.0d0
    
    ! Local Variables
    integer(int32) :: i, k
    real(real64) :: t(npts), x(npts), y(npts)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: d1, d2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Build the signal
    t = linspace(0.0d0, 1.0d0, npts)
    call random_number(x)
    x = 0.25d0 * (x - 0.5d0) + sin(2.0d0 * pi * f1 * t) + &
        0.5 * sin(2.0d0 * pi * f2 * t)

    ! Apply the filter
    y = tv_filter(x, 0.5d0)

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

    call plt%draw()
end program