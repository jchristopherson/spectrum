module plot_helper
    use iso_fortran_env
    use fplot_core
    implicit none

contains
    subroutine plot_time_responses(t, x, xf)
        real(real64), intent(in), dimension(:) :: t
        real(real64), intent(in), dimension(size(t)) :: x
        real(real64), intent(in), dimension(size(t)) :: xf

        type(plot_2d) :: plt
        type(plot_data_2d) :: pd1, pd2
        class(legend), pointer :: lgnd
        class(plot_axis), pointer :: xAxis, yAxis

        call plt%initialize()
        lgnd => plt%get_legend()
        xAxis => plt%get_x_axis()
        yAxis => plt%get_y_axis()
        call xAxis%set_title("t")
        call yAxis%set_title("x(t)")
        call lgnd%set_is_visible(.true.)

        call pd1%define_data(t, x)
        call pd1%set_name("Original")
        call plt%push(pd1)

        call pd2%define_data(t, xf)
        call pd2%set_name("Filtered")
        call plt%push(pd2)

        call plt%draw()
    end subroutine

    subroutine plot_spectral_response(f, x, xf)
        real(real64), intent(in), dimension(:) :: f
        real(real64), intent(in), dimension(size(f)) :: x
        real(real64), intent(in), dimension(size(f)) :: xf

        type(plot_2d) :: plt
        type(plot_data_2d) :: pd1, pd2
        class(legend), pointer :: lgnd
        class(plot_axis), pointer :: xAxis, yAxis

        call plt%initialize()
        lgnd => plt%get_legend()
        xAxis => plt%get_x_axis()
        yAxis => plt%get_y_axis()
        call xAxis%set_title("f")
        call yAxis%set_title("X(f)^2")
        call yAxis%set_is_log_scaled(.true.)
        call lgnd%set_is_visible(.true.)

        call pd1%define_data(f, x)
        call pd1%set_line_width(2.0)
        call pd1%set_name("Original")
        call plt%push(pd1)

        call pd2%define_data(f, xf)
        call pd2%set_line_width(2.0)
        call pd2%set_name("Filtered")
        call plt%push(pd2)

        call plt%draw()
    end subroutine
end module

program example
    use iso_fortran_env
    use spectrum
    use plot_helper
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 1000
    integer(int32), parameter :: ntaps = 10
    real(real64), parameter :: cutoff_hz = 5.0d1
    real(real64), parameter :: sample_hz = 1.024d3
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: f1 = 1.0d1
    real(real64), parameter :: f2 = 1.25d2

    ! Local Variables
    integer(int32) :: i
    real(real64) :: dt, t0, df, t(npts), x(npts), xf(npts)
    real(real64), allocatable, dimension(:) :: alp, blp
    real(real64), allocatable, dimension(:) :: f, p1, p1f
    type(hamming_window) :: win

    ! Design the filter
    call design_fir_filter(ntaps, cutoff_hz, sample_hz, blp, alp, ftype = LOW_PASS_FILTER)

    ! Display the coefficients
    print "(A)", "FILTER COEFFICIENTS:"
    print "(A)", "B = "
    do i = 1, size(blp)
        print "(A, A, I0, A, F8.3)", achar(9), "b(", i, ") = ", blp(i)
    end do
    print "(A)", "A = "
    do i = 1, size(alp)
        print "(A, A, I0, A, F8.3)", achar(9), "a(", i, ") = ", alp(i)
    end do

    ! Define the signal
    call random_number(x)
    dt = 1.0d0 / sample_hz
    t0 = -dt
    do i = 1, npts
        t(i) = t0 + dt
        t0 = t(i)
        x(i) = 0.5d0 * (x(i) - 0.5d0) + sin(2.0d0 * pi * f1 * t(i)) + &
            0.5d0 * sin(2.0d0 * pi * f2 * t(i))
    end do

    ! Apply the filter
    xf = filter(blp, alp, x)

    ! Plot the time histories
    call plot_time_responses(t, x, xf)

    ! Plot the spectral responses
    win%size = 512
    p1 = psd(win, x, fs = sample_hz)
    p1f = psd(win, xf, fs = sample_hz)
    allocate(f(size(p1)))
    df = frequency_bin_width(sample_hz, win%size)
    f = (/ (i * df, i = 0, size(p1) - 1) /)
    call plot_spectral_response(f, p1, p1f)
end program