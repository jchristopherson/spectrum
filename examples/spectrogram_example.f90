program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: window_size = 512
    real(real64), parameter :: fs = 2048.0d0
    real(real64), parameter :: f0 = 1.0d2
    real(real64), parameter :: f1 = 1.0d3
    real(real64), parameter :: duration = 50.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i, npts
    integer(int32), allocatable, dimension(:) :: offsets
    real(real64) :: k, df
    complex(real64), allocatable, dimension(:,:) :: rst
    real(real64), allocatable, dimension(:) :: t, x, f, s
    real(real64), allocatable, dimension(:,:) :: mag
    real(real64), allocatable, dimension(:,:,:) :: xy
    type(hamming_window) :: win

    ! Plot Variables
    type(surface_plot) :: plt
    type(surface_plot_data) :: pd
    class(plot_axis), pointer :: xAxis, yAxis
    type(rainbow_colormap) :: map

    ! Create the exponential chirp signal
    npts = floor(duration * fs) + 1
    t = linspace(0.0d0, duration, npts)
    k = (f1 / f0)**(1.0 / duration)
    x = sin(2.0d0 * pi * f0 * (k**t - 1.0d0) / log(k))

    ! Determine sampling frequency parameters
    df = frequency_bin_width(fs, window_size)

    ! Define the window
    win%size = window_size

    ! Compute the spectrogram of x
    rst = spectrogram(win, x, offsets)

    ! Compute the magnitude, along with each frequency and time point
    mag = abs(rst)

    allocate(f(size(mag, 1)))
    f = (/ (df * i, i = 0, size(f) - 1) /)

    allocate(s(size(mag, 2)))
    do i = 1, size(s)
        if (i == 1) then
            s(i) = 0.5d0 * window_size / fs
        else
            s(i) = i * (offsets(i) - offsets(i-1)) / fs
        end if
    end do
    xy = meshgrid(s, f)

    ! Plot the results
    call plt%initialize()
    call plt%set_colormap(map)
    call plt%set_use_map_view(.true.)
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("Time [s]")
    call yAxis%set_title("Frequency [Hz]")
    call yAxis%set_autoscale(.false.)
    call yAxis%set_limits(0.0d0, f(size(mag, 1)))
    
    call pd%define_data(xy(:,:,1), xy(:,:,2), mag)
    call plt%push(pd)
    call plt%draw()
end program