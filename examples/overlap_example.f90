program example
    use iso_fortran_env
    use spectrum
    use fftpack
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
    integer(int32) :: i, nxfrm, noverlap
    real(real64) :: dt, t(npts), x(npts), buffer(winsize), w(winsize), sumw, df
    real(real64), allocatable, dimension(:) :: freq, pwr
    complex(real64) :: xfrm(winsize)
    type(hamming_window) :: win

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd
    class(plot_axis), pointer :: xAxis, yAxis

    ! Build the signal
    dt = 1.0d0 / sample_rate
    t = (/ (dt * i, i = 0, npts - 1) /)
    x = amp1 * sin(2.0d0 * pi * freq1 * t - pi * phase1 / 1.8d2) + &
        amp2 * sin(2.0d0 * pi * freq2 * t - pi * phase2 / 1.8d2)

    ! Define the window
    win%size = winsize
    w = (/ (win%evaluate(i), i = 0, winsize - 1) /)
    sumw = sum(w) / 2.0d0

    ! Determine how many overlapped segments
    noverlap = compute_overlap_segment_count(npts, winsize)

    ! Compute the transforms of each overlapped segment
    nxfrm = compute_transform_length(winsize)
    allocate(pwr(nxfrm), source = 0.0d0)
    do i = 1, noverlap
        ! Extract the relevant portion of the original signal
        call overlap(x, i, winsize, buffer)

        ! Apply the window
        buffer = buffer * w

        ! Compute the transform and add the result to the previous segment.
        ! Be sure to scale each transform to account for the windowing.
        xfrm = fft(cmplx(buffer, 0.0d0, real64)) / sumw

        ! Accumulate the power spectrum - only keeping a 1-sided spectrum
        pwr = pwr + abs(xfrm(1:nxfrm))**2
    end do

    ! Average the accumulated transforms
    pwr = pwr / noverlap

    ! Compute the frequency, and get the magnitude and phase information
    df = frequency_bin_width(sample_rate, winsize)
    allocate(freq(nxfrm))
    freq = (/ (df * i, i = 0, nxfrm - 1) /)

    ! Set up the plots
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("f [Hz]")
    call yAxis%set_title("|X|^{2}")

    call xAxis%set_is_log_scaled(.true.)
    
    call xAxis%set_use_default_tic_label_format(.false.)
    call xAxis%set_tic_label_format("%0.0e")

    call pd%define_data(freq, pwr)
    call pd%set_line_width(2.0)
    call plt%push(pd)
    call plt%draw()
end program