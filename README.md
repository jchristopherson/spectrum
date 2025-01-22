# spectrum
Spectrum is a library containing signal analysis routines with a focus towards spectral routines.

## Status
[![CMake](https://github.com/jchristopherson/spectrum/actions/workflows/cmake.yml/badge.svg)](https://github.com/jchristopherson/spectrum/actions/workflows/cmake.yml)
[![Actions Status](https://github.com/jchristopherson/spectrum/workflows/fpm/badge.svg)](https://github.com/jchristopherson/spectrum/actions)

## Available Operations
- Spectral Estimators (power spectral density, cross spectral density, etc.)
- Short Time Fourier Transform (STFT)
- Windowing
- Convolution
- Offset Removal (Mean Removal)
- Filtering
    - Low-Pass, High-Pass, Band-Pass, & Band-Stop
    - Total Variation
    - Gaussian
    - Moving Average
- Upsampling & Downsampling
- Transfer Function Estimation
- Differentiation & Integration

## Documentation
The documentation can be found [here](https://jchristopherson.github.io/spectrum/).

## Building Spectrum
[CMake](https://cmake.org/)This library can be built using CMake.  For instructions see [Running CMake](https://cmake.org/runningcmake/).

[FPM](https://github.com/fortran-lang/fpm) can also be used to build this library using the provided fpm.toml.
```txt
fpm build
```
The SPECTRUM library can be used within your FPM project by adding the following to your fpm.toml file.
```toml
[dependencies]
spectrum = { git = "https://github.com/jchristopherson/spectrum" }
```

## External Libraries
The FPLOT library depends upon the following libraries.
- [FERROR](https://github.com/jchristopherson/ferror)
- [FFTPACK](https://github.com/fortran-lang/fftpack)
- [BLAS](http://www.netlib.org/blas/)
- [LAPACK](http://www.netlib.org/lapack/)
- [LINALG](https://github.com/jchristopherson/linalg)

Spectrum contains a large selection of signal processing routines.  The following examples illustrate some of the spectral capabilities in addition to some of the filtering options available.

## Filtering & PSD Example
This example illustrates some available filtering options on a randomly generated signal by comparing the power spectrums of the unfiltered vs. filtered signals.  This example utilizes the [fplot](https://github.com/jchristopherson/fplot) library for plotting.
```fortran
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
    call term%set_window_height(900)
    
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
```
This example produces the following plot.

![](images/filter_comparison.png?raw=true)

## Transfer Function Example
This example illustrates the calculation of a single-input/single-output (SISO) transfer function from a signal stored in a CSV file.  This example utilizes the [fplot](https://github.com/jchristopherson/fplot) library for plotting and the [fortran-csv-module](https://github.com/jacobwilliams/fortran-csv-module) library for reading the CSV file.
```fortran
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
```
This example produces the following plots.

![](images/siso_tf_example_time_history.png?raw=true)
![](images/siso_tf_example_bode.png?raw=true)

## STFT Example
This example illustrates the STFT capabilities.
```fortran
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
    type(hann_window) :: win

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
    rst = stft(win, x, offsets)

    ! Compute the magnitude, along with each frequency and time point
    mag = abs(rst)

    allocate(f(size(mag, 1)))
    f = (/ (df * i, i = 0, size(f) - 1) /)

    allocate(s(size(mag, 2)))
    do i = 1, size(s)
        if (i == 1) then
            s(i) = offsets(i) / fs
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
```
This example produces the following plot.

![](images/spectrogram_example_1.png?raw=true)

## References
1. Welch, P.D. (1967). The Use of Fast Fourier Transform for the Estimation of Power Spectra: A Method Based on Time Averaging Over Short, Modified Periodograms. IEEE Transactions on Audio and Electroacoustics, AU-15 (2): 70-73.
2. Stoica, Petre, and Randolph Moses. Spectral Analysis of Signals. Upper Saddle River, NJ: Prentice Hall, 2005.
3. William H. Press, Saul A. Teukolsky, William T. Vetterling, and Brian P. Flannery. 2007. Numerical Recipes 3rd Edition: The Art of Scientific Computing (3rd. ed.). Cambridge University Press, USA.
4. Millioz, Fabien & Martin, Nadine. (2011). Circularity of the STFT and Spectral Kurtosis for Time-Frequency Segmentation in Gaussian Environment. Signal Processing, IEEE Transactions on. 59. 515 - 524. 10.1109/TSP.2010.2081986. 
5. Wikipedia Contributors. (2022, February 9). Welch’s method. Wikipedia; Wikimedia Foundation. https://en.wikipedia.org/wiki/Welch%27s_method
6. Spectral density. (2023, April 9). Wikipedia. https://en.wikipedia.org/wiki/Spectral_density#Cross-spectral_density
7. Wikipedia Contributors. (2019, April 12). Short-time Fourier transform. Wikipedia; Wikimedia Foundation. https://en.wikipedia.org/wiki/Short-time_Fourier_transform
8. Window function. (2020, December 12). Wikipedia. https://en.wikipedia.org/wiki/Window_function
9. Wikipedia Contributors. (2019, March 22). Gaussian filter. Wikipedia; Wikimedia Foundation. https://en.wikipedia.org/wiki/Gaussian_filter
10. Selesnick Andilker Bayram, I. (2010). Total Variation Filtering. https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/TV_filtering.pdf
11. Unavane, T., & Panse. (2015). New Method for Online Frequency Response Function Estimation Using Circular Queue. INTERNATIONAL JOURNAL for RESEARCH in EMERGING SCIENCE and TECHNOLOGY, 2. https://ijrest.net/downloads/volume-2/issue-6/pid-ijrest-26201530.pdf

‌