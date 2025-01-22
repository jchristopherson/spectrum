program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: winsize = 1024
    integer(int32), parameter :: n = 10000
    real(real64), parameter :: fs = 1.0d3
    real(real64), parameter :: fc1 = 1.0d2
    real(real64), parameter :: fc2 = 3.0d2

    ! Local Variables
    integer(int32) :: i
    real(real64) :: df, x(n), xlp(n), xhp(n), xbp(n), xbs(n)
    real(real64), allocatable, dimension(:) :: freq, p, plp, php, pbp, pbs
    type(hamming_window) :: win

    ! Plot Variables
    type(multiplot) :: freqPlot
    type(plot_2d) :: plt, pltlp, plthp, pltbp, pltbs
    class(terminal), pointer :: term
    type(plot_data_2d) :: pd, pdlp, pdhp, pdbp, pdbs

    ! Build the signal and remove the mean
    call random_number(x)
    x = remove_mean(x)

    ! Filter the signals
    xlp = sinc_filter(fc1, fs, x, ftype = LOW_PASS_FILTER)
    xhp = sinc_filter(fc1, fs, x, ftype = HIGH_PASS_FILTER)
    xbp = sinc_filter(fc1, fs, x, ftype = BAND_PASS_FILTER, fc2 = fc2)
    xbs = sinc_filter(fc1, fs, x, ftype = BAND_STOP_FILTER, fc2 = fc2)

    ! Compute PSD's
    win%size = winsize
    p = psd(win, x, fs = fs)
    plp = psd(win, xlp, fs = fs)
    php = psd(win, xhp, fs = fs)
    pbp = psd(win, xbp, fs = fs)
    pbs = psd(win, xbs, fs = fs)

    ! Build a corresponding array of frequency values
    df = frequency_bin_width(fs, winsize)
    allocate(freq(size(p)))
    freq = (/ (df * i, i = 0, size(p) - 1) /)

    ! Create the plot
    call freqPlot%initialize(5, 1)
    call plt%initialize()
    call pltlp%initialize()
    call plthp%initialize()
    call pltbp%initialize()
    call pltbs%initialize()
    term => freqPlot%get_terminal()
    call term%set_window_height(1200)
    
    call plt%set_title("Original Signal")
    call pltlp%set_title("Low Pass Filtered Signal")
    call plthp%set_title("High Pass Filtered Signal")
    call pltbp%set_title("Band Pass Filtered Signal")
    call pltbs%set_title("Band Stop Filtered Signal")

    call pd%define_data(freq, p)
    call plt%push(pd)
    call freqPlot%set(1, 1, plt)

    call pdlp%define_data(freq, plp)
    call pltlp%push(pdlp)
    call freqPlot%set(2, 1, pltlp)

    call pdhp%define_data(freq, php)
    call plthp%push(pdhp)
    call freqPlot%set(3, 1, plthp)

    call pdbp%define_data(freq, pbp)
    call pltbp%push(pdbp)
    call freqPlot%set(4, 1, pltbp)

    call pdbs%define_data(freq, pbs)
    call pltbs%push(pdbs)
    call freqPlot%set(5, 1, pltbs)

    call freqPlot%draw()
end program