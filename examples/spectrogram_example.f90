program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: npts = 16384
    integer(int32), parameter :: window_size = 128
    integer(int32), parameter :: noverlap = 64
    real(real64), parameter :: f0 = 1.0d1
    real(real64), parameter :: f1 = 1.0d2
    real(real64), parameter :: duration = 5.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i
    real(real64) :: t(npts), x(npts), k
    complex(real64), allocatable, dimension(:,:) :: rst
    real(real64), allocatable, dimension(:,:) :: mag
    type(hamming_window) :: win

    ! Create the exponential chirp signal
    t = linspace(0.0d0, duration, npts)
    k = (f1 / f0)**(1.0 / duration)
    x = sin(2.0d0 * pi * f0 * (k**t - 1.0d0) / log(k))

    ! Define the window
    win%size = window_size

    ! Compute the spectrogram of x
    rst = spectrogram(win, x, noverlap)

    ! Compute the magnitude, along with each frequency and time point
    mag = abs(rst)
end program