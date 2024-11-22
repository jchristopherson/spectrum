program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Variables
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: dt = 1.0d-3
    integer(int32) :: i
    real(real64) :: t(npts), x(npts), xm(npts)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Generate the signal
    t = (/ (i * dt, i = 0, npts - 1) /)
    call random_number(x)

    ! Remove the mean
    xm = remove_mean(x)

! ------------------------------------------------------------------------------
    ! Plot
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(t, x)
    call pd1%set_name("Original")
    call plt%push(pd1)

    call pd2%define_data(t, xm)
    call pd2%set_name("Mean Removed")
    call plt%push(pd2)

    call plt%draw()
end program