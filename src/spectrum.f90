! TO DO:
! - windowed FFT (https://en.wikipedia.org/wiki/Short-time_Fourier_transform)
! - wavelet transforms???
! - a routine for creating a frequency vector

!> @brief
module spectrum
    use iso_fortran_env
    use ferror
    implicit none
    private
    public :: window
    public :: window_function
    public :: rectangular_window
    public :: hann_window
    public :: hamming_window
    public :: welch_window
    public :: blackman_harris_window
    public :: flat_top_window
    public :: compute_transform_length
    public :: frequency_bin_width
    public :: is_power_of_two
    public :: next_power_of_two
    public :: psd
    public :: csd
    public :: spectrogram
    public :: convolve
    public :: gaussian_filter
    public :: tv_filter
    public :: compute_overlap_segment_count
    public :: overlap
    public :: SPCTRM_MEMORY_ERROR
    public :: SPCTRM_INVALID_INPUT_ERROR
    public :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR
    public :: SPCTRM_FULL_CONVOLUTION
    public :: SPCTRM_CENTRAL_CONVOLUTION

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: SPCTRM_MEMORY_ERROR = 10000
    integer(int32), parameter :: SPCTRM_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR = 10002

    !> A flag for requesting a full convolution.
    integer(int32), parameter :: SPCTRM_FULL_CONVOLUTION = 50000
    !> A flag for requesting the central portion of the convolution that is the
    !! same length as the input signal.
    integer(int32), parameter :: SPCTRM_CENTRAL_CONVOLUTION = 50001

! ******************************************************************************
! SPECTRUM_WINDOWS.F90
! ------------------------------------------------------------------------------
    !> @brief Defines a window.
    type, abstract :: window
        ! Window size.
        integer(int32), public :: size = 0
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate(class(window) this, integer(int32) bin)
        !! @endcode
        !!
        !! @param[in] this The @ref window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure(window_function), public, deferred, pass :: evaluate
    end type

    interface
        pure function window_function(this, bin) result(rst)
            use iso_fortran_env, only : real64, int32
            import window
            class(window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function
    end interface

    !> @brief Defines a rectangular window
    type, extends(window) :: rectangular_window
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate( &
        !!  class(rectangular_window) this, &
        !!  integer(int32) bin &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref rectangular_window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure, public :: evaluate => rw_eval
    end type

    !> @brief Defines a Hann window.
    !!
    !! @par
    !! A Hann window is defined as follows.
    !! 
    !! @par
    !! \f$ w(j) = \frac{1}{2} \left( 1 - \cos \left( \frac{2 \pi j}{n} \right) \right) \f$
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    type, extends(window) :: hann_window
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate( &
        !!  class(hann_window) this, &
        !!  integer(int32) bin &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref hann_window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure, public :: evaluate => hann_eval
    end type

    !> @brief Defines a Hamming window.
    !!
    !! @par
    !! A Hamming window is defined as follows.
    !! 
    !! @par
    !! \f$ w(j) = 0.54 - 0.46 \cos \left( \frac{2 \pi j}{n} \right) \f$
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    type, extends(window) :: hamming_window
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate( &
        !!  class(hamming_window) this, &
        !!  integer(int32) bin &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref hamming_window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure, public :: evaluate => hamming_eval
    end type

    !> @brief Defines a Welch window.
    !!
    !! @par
    !! A Welch window is defined as follows.
    !! 
    !! @par
    !! \f$ w(j) = 1 - \left( \frac{j - \frac{n}{2} }{ \frac{n}{2} } \right)^2 \f$
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    type, extends(window) :: welch_window
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate( &
        !!  class(welch_window) this, &
        !!  integer(int32) bin &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref welch_window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure, public :: evaluate => welch_eval
    end type

    !> @brief Defines a Blackman-Harris window.
    !!
    !! @par
    !! A Blackman-Harris window is defined as follows.
    !! 
    !! @par
    !! \f$ w(j) = 0.3635819 - 0.4891775 \cos \left( \frac{2 \pi j}{n} \right)
    !! + 0.1365995 \cos \left( \frac{4 \pi j}{n} \right) - 0.0106411
    !! \cos \left( \frac{6 \pi j}{n} \right) \f$
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    type, extends(window) :: blackman_harris_window
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate( &
        !!  class(blackman_harris_window) this, &
        !!  integer(int32) bin &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref blackman_harris_window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure, public :: evaluate => bhw_eval
    end type

    !> @brief Defines a flat top window.
    !!
    !! @par
    !! A flat top window is defined as follows.
    !!
    !! @par
    !! \f$ w(j) = 0.21557895 - 
    !! 0.41663158 \cos \left( \frac{2 \pi j}{N} \right) +
    !! 0.277263158 \cos \left( \frac{4 \pi j}{N} \right) -
    !! 0.083578947 \cos \left( \frac{6 \pi j}{N} \right)  +
    !! 0.006947368 \cos \left( \frac{8 \pi j}{N} \right) \f$
    !!
    !! @par See Also
    !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    type, extends(window) :: flat_top_window
    contains
        !> @brief Evaluates the window function.
        !!
        !! @par Syntax
        !! @code{.f90}
        !! real(real64) pure function evaluate( &
        !!  class(flat_top_window) this, &
        !!  integer(int32) bin &
        !! )
        !! @endcode
        !!
        !! @param[in] this The @ref flat_top_window object.
        !! @param[in] bin The index or bin number (0 <= @p bin <= n) where n
        !!  is the size of the window.
        !!
        !! @return The window function value.
        procedure, public :: evaluate => ftw_eval
    end type

    interface
        pure module function rw_eval(this, bin) result(rst)
            class(rectangular_window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function

        pure module function hann_eval(this, bin) result(rst)
            class(hann_window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function

        pure module function hamming_eval(this, bin) result(rst)
            class(hamming_window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function

        pure module function welch_eval(this, bin) result(rst)
            class(welch_window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function

        pure module function bhw_eval(this, bin) result(rst)
            class(blackman_harris_window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function

        pure module function ftw_eval(this, bin) result(rst)
            class(flat_top_window), intent(in) :: this
            integer(int32), intent(in) :: bin
            real(real64) :: rst
        end function
    end interface

! ******************************************************************************
! SPECTRUM_ROUTINES.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the length of the positive half of a discrete Fourier
    !! transform for a specific signal length.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! integer(int32) pure function compute_transform_length(integer(int32) n)
    !! @endcode
    !!
    !! @param[in] n The signal length (length of the signal put forth to the
    !!  Fourier transform).
    !! @return The length of the positive half of the discrete Fourier
    !!  transform.
    !!
    !! @par Remarks
    !! The length of the discrete Fourier transform is determined as follows.
    !! @par
    !! \f$ m = \frac{n_{signal}}{2} + 1 \f$ if \f$ n_{signal} \f$ is even; else,
    !! \f$ m = \frac{n_{signal} + 1}{2} \f$ if \f$ n_{signal} \f$ is odd.
    interface compute_transform_length
        module procedure :: compute_xfrm_length_1
    end interface

    !> @brief Computes the bin width for a discrete frequency spectrum based
    !! upon the data sample rate.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! real(real64) pure elemental function frequency_bin_width( &
    !!  real(real64) fs, &
    !!  integer(int32) n &
    !! )
    !! @endcode
    !!
    !! @param[in] fs The rate at which the signal was sampled.  The units of
    !!  this value will be the units of the output.
    !! @param[in] n The signal length (length of the signal put forth to the
    !!  Fourier transform).
    !!
    !! @par Remarks
    !! The frequency bin width (\f$ \Delta f \f$) is computed as follows.
    !! \f$ \Delta f = \frac{f_{s}}{2 \left(m - 1 \right)} \f$
    !! where,
    !! \f$ m = \frac{n_{window}}{2} + 1 \f$ if \f$ n_{window} \f$ is even; else,
    !! \f$ m = \frac{n_{window} + 1}{2} \f$ if \f$ n_{window} \f$ is odd.
    !!
    !! @par
    !! \f$ m \f$ can be computed using @ref compute_transform_length
    interface frequency_bin_width
        module procedure :: frequency_bin_width_1
    end interface

    !> @brief Tests to see if a value is an integer power of two.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! logical pure pure elemental function is_power_of_two(integer(int32) n)
    !! @endcode
    !!
    !! @param[in] n The value to test.
    !! @return Returns true if @p n is an integer power of two; else, returns 
    !!  false.
    interface is_power_of_two
        module procedure :: is_power_of_two_1
    end interface

    !> @brief Provides the next higher integer power of two.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! integer(int32) pure elemental function next_power_of_two(integer(int32) n)
    !! @endcode
    !!
    !! @param[in] n The value to test.
    !! @return The next power of two higher than @p n.  If @p n is already
    !! a power of two, it's value is simply returned.  For instance, if @p
    !! is set to 128, then a value of 7 is returned.  However, if a value
    !! of 129 is supplied, then a value of 8 is returned.
    interface next_power_of_two
        module procedure :: next_power_of_two_1
    end interface

    interface
        pure module function compute_xfrm_length_1(n) result(rst)
            integer(int32), intent(in) :: n
            integer(int32) :: rst
        end function

        pure module elemental function frequency_bin_width_1(fs, n) result(rst)
            real(real64), intent(in) :: fs
            integer(int32), intent(in) :: n
            real(real64) :: rst
        end function

        pure module elemental function is_power_of_two_1(n) result(rst)
            integer(int32), intent(in) :: n
            logical :: rst
        end function

        pure module elemental function next_power_of_two_1(n) result(rst)
            integer(int32), intent(in) :: n
            integer(int32) :: rst
        end function
    end interface

! ******************************************************************************
! SPECTRUM_PSD.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the power spectral density (PSD) of a signal via Welch's
    !! averaged method.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function psd( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  optional real(real64) fs, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The @ref window to apply.  The size of the window must be
    !!  non-zero and positive-valued.
    !! @param[in] x The signal to transform.
    !! @param[in] fs An optional input, that if supplied, allows for 
    !!  normalization of the computed PSD by the frequency resolution.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p win is not sized 
    !!      appropriately.
    !!
    !! @return An array containing the discrete PSD estimate.  The estimate is
    !!  returned at discrete frequency intervals that can be determined as 
    !!  follows.
    !!
    !! @par
    !! \f$ \Delta f = \frac{f_{s}}{2 \left(m - 1 \right)} \f$
    !! where,
    !! \f$ m = \frac{n_{window}}{2} + 1 \f$ if \f$ n_{window} \f$ is even; else,
    !! \f$ m = \frac{n_{window} + 1}{2} \f$ if \f$ n_{window} \f$ is odd.
    !!
    !! @par
    !! \f$ m \f$ can be computed using @ref compute_transform_length, and
    !! \f$ \Delta f \f$ can be computed using @ref frequency_bin_width.
    !!
    !! @par References
    !! - Welch, P.D. "The Use of Fast Fourier Transform for the Estimation of 
    !!  Power Spectra: A Method Based on Time Averaging Over Short, Modified 
    !!  Periodograms." IEEE Transactions on Audio and Electroacoustics, 
    !!  AU-15 (2): 70-73, 1967.
    !! - [Wikipedia - Welch's Method](https://en.wikipedia.org/wiki/Welch%27s_method)
    interface psd
        module procedure :: psd_welch
    end interface

    interface
        module function psd_welch(win, x, fs, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:)
            real(real64), intent(in), optional :: fs
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_CSD.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the cross spectral density (CSD) of a signal via Welch's
    !! averaged method.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function csd( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  real(real64) y(:), &
    !!  optional real(real64) fs, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The @ref window to apply.  The size of the window must be
    !!  non-zero and positive-valued.
    !! @param[in] x The first N-element signal to transform.
    !! @param[in] y The second N-element signal to transform.
    !! @param[in] fs An optional input, that if supplied, allows for 
    !!  normalization of the computed CSD by the frequency resolution.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p win is not sized 
    !!      appropriately.
    !!
    !! @return An array containing the discrete CSD estimate.  The estimate is
    !!  returned at discrete frequency intervals that can be determined as 
    !!  follows.
    !!
    !! @par
    !! \f$ \Delta f = \frac{f_{s}}{2 \left(m - 1 \right)} \f$
    !! where,
    !! \f$ m = \frac{n_{window}}{2} + 1 \f$ if \f$ n_{window} \f$ is even; else,
    !! \f$ m = \frac{n_{window} + 1}{2} \f$ if \f$ n_{window} \f$ is odd.
    !!
    !! @par
    !! \f$ m \f$ can be computed using @ref compute_transform_length, and
    !! \f$ \Delta f \f$ can be computed using @ref frequency_bin_width.
    !!
    !! @par References
    !! - Welch, P.D. "The Use of Fast Fourier Transform for the Estimation of 
    !!  Power Spectra: A Method Based on Time Averaging Over Short, Modified 
    !!  Periodograms." IEEE Transactions on Audio and Electroacoustics, 
    !!  AU-15 (2): 70-73, 1967.
    !! - [Wikipedia - Welch's Method](https://en.wikipedia.org/wiki/Welch%27s_method)
    interface csd
        module procedure :: csd_welch
    end interface

    interface
        module function csd_welch(win, x, y, fs, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:), y(:)
            real(real64), intent(in), optional :: fs
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_FFT.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the spectrogram of a signal.  Only the positive half of
    !! the transform is returned.  The symmetry of the transform may be 
    !! exploited if the both the negative and positive frequency components
    !! are desired.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable complex(real64) rst(:,:) function spectrogram( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  optional allocatable integer(int32) offsets(:), &
    !!  optional logical par, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The window to apply.
    !! @param[in] x The signal to analyze.  The signal must be longer than the
    !!  size of the window @p win.
    !! @param[out] offsets An optional allocatable array that, if supplied, will
    !!  be filled with the starting indices of each window segment.
    !! @param[in] par An optional input that will utilize OpenMP to parallelize
    !!  the transform operations if set to true.  If set to false, each 
    !!  transform will be computed in a serial manner.  The default is false 
    !!  unless the number of transforms exceeds 50.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if the signal in @p x is too short
    !!      relative to the window size in @p win.
    !!
    !! @return An M-by-N matrix containing the M-element complex-valued 
    !!  transforms for each of the N time points studied.  M is the size of the
    !!  positive half of the transform, and N is the total number of transformed
    !!  segments.
    !!
    !! @par Example
    !! The following example illustrates how to compute the spectrogram of an
    !! exponential chirp signal.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use spectrum
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: window_size = 512
    !!     real(real64), parameter :: fs = 2048.0d0
    !!     real(real64), parameter :: f0 = 1.0d2
    !!     real(real64), parameter :: f1 = 1.0d3
    !!     real(real64), parameter :: duration = 50.0d0
    !!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i, npts
    !!     integer(int32), allocatable, dimension(:) :: offsets
    !!     real(real64) :: k, df
    !!     complex(real64), allocatable, dimension(:,:) :: rst
    !!     real(real64), allocatable, dimension(:) :: t, x, f, s
    !!     real(real64), allocatable, dimension(:,:) :: mag
    !!     real(real64), allocatable, dimension(:,:,:) :: xy
    !!     type(hann_window) :: win
    !!
    !!     ! Plot Variables
    !!     type(surface_plot) :: plt
    !!     type(surface_plot_data) :: pd
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!     type(rainbow_colormap) :: map
    !!
    !!     ! Create the exponential chirp signal
    !!     npts = floor(duration * fs) + 1
    !!     t = linspace(0.0d0, duration, npts)
    !!     k = (f1 / f0)**(1.0 / duration)
    !!     x = sin(2.0d0 * pi * f0 * (k**t - 1.0d0) / log(k))
    !!
    !!     ! Determine sampling frequency parameters
    !!     df = frequency_bin_width(fs, window_size)
    !!
    !!     ! Define the window
    !!     win%size = window_size
    !!
    !!     ! Compute the spectrogram of x
    !!     rst = spectrogram(win, x, offsets)
    !!
    !!     ! Compute the magnitude, along with each frequency and time point
    !!     mag = abs(rst)
    !!
    !!     allocate(f(size(mag, 1)))
    !!     f = (/ (df * i, i = 0, size(f) - 1) /)
    !!
    !!     allocate(s(size(mag, 2)))
    !!     do i = 1, size(s)
    !!         if (i == 1) then
    !!             s(i) = offsets(i) / fs
    !!         else
    !!             s(i) = i * (offsets(i) - offsets(i-1)) / fs
    !!         end if
    !!     end do
    !!     xy = meshgrid(s, f)
    !!
    !!     ! Plot the results
    !!     call plt%initialize()
    !!     call plt%set_colormap(map)
    !!     call plt%set_use_map_view(.true.)
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!
    !!     call xAxis%set_title("Time [s]")
    !!     call yAxis%set_title("Frequency [Hz]")
    !!     call yAxis%set_autoscale(.false.)
    !!     call yAxis%set_limits(0.0d0, f(size(mag, 1)))
    !!
    !!     call pd%define_data(xy(:,:,1), xy(:,:,2), mag)
    !!     call plt%push(pd)
    !!     call plt%draw()
    !! end program
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html spectrogram_example_1.png
    interface spectrogram
        module procedure :: stft
    end interface

    interface
        module function stft(win, x, offsets, par, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:)
            integer(int32), intent(out), optional, allocatable :: offsets(:)
            logical, intent(in), optional :: par
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable :: rst(:,:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_CONVOLVE.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the convolution of a signal and kernel.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) convolve( &
    !!  real(real64) x(:), &
    !!  real(real64) y(:), &
    !!  optional integer(int32) method, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] x The N-element signal.
    !! @param[in] y The M-element kernel.
    !! @param[in] method An optional input that dictates the expected
    !!  convolution result.  The following options are available.
    !!  - SPCTRM_FULL_CONVOLUTION: The full convolution results are provided, 
    !!      including the portions polluted courtesy of the zero-padding and
    !!      the corresponding wrap-around effects.  The length of this output
    !!      is N + M - 1.
    !!  - SPCTRM_CENTRAL_CONVOLUTION: The N-element result containing the 
    !!      convolved signal not poluted by the zero-padding and corresponding
    !!      wrap-around effects.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!
    !! @return The convolved result.
    !!
    !! @par Example
    !! The following example illustrates how to use convolution to apply a 
    !! Gaussian filter to noisy data.  This example utilizes a window that
    !! spans 3 sigma.
    !! @code{.f90}
    !! program example
    !!     use iso_fortran_env
    !!     use spectrum
    !!     use fplot_core
    !!     implicit none
    !!
    !!     ! Parameters
    !!     integer(int32), parameter :: nkernel = 21
    !!     integer(int32), parameter :: npts = 10000
    !!     real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    !!     real(real64), parameter :: f1 = 5.0d0
    !!     real(real64), parameter :: f2 = 1.5d1
    !!     real(real64), parameter :: alpha = 3.0d0
    !!
    !!     ! Local Variables
    !!     integer(int32) :: i, k
    !!     real(real64) :: t(npts), x(npts), y(npts), g(nkernel), sumg
    !!
    !!     ! Plot Variables
    !!     type(plot_2d) :: plt
    !!     type(plot_data_2d) :: d1, d2
    !!     class(plot_axis), pointer :: xAxis, yAxis
    !!     class(legend), pointer :: lgnd
    !!
    !!     ! Build the signal
    !!     t = linspace(0.0d0, 1.0d0, npts)
    !!     call random_number(x)
    !!     x = 0.25d0 * (x - 0.5d0) + sin(2.0d0 * pi * f1 * t) + &
    !!         0.5 * sin(2.0d0 * pi * f2 * t)
    !!
    !!     ! Define the Gaussian kernel
    !!     ! Let sigma = (NKERNEL - 1) / (2 alpha) such that:
    !!     ! exp(-k**2 / (2 sigma**2)) = exp(-(1/2) (alpha k / (NKERNEL - 1) / 2)**2)
    !!     k = -(nkernel - 1) / 2  ! NKERNEL should be odd, or made to odd, for a symmetric window
    !!     sumg = 0.0d0
    !!     do i = 1, nkernel 
    !!         g(i) = exp(-0.5d0 * (2.0d0 * alpha * k / (nkernel - 1.0d0))**2)
    !!         k = k + 1
    !!         sumg = sumg + g(i)
    !!     end do
    !!
    !!     ! Normalize the kernel to have a sum of one
    !!     g = g / sumg
    !!
    !!     ! Compute the convolution and keep only the non-poluted data
    !!     y = convolve(x, g, SPCTRM_CENTRAL_CONVOLUTION)
    !!
    !!     ! Create the plot
    !!     call plt%initialize()
    !!     xAxis => plt%get_x_axis()
    !!     yAxis => plt%get_y_axis()
    !!     lgnd => plt%get_legend()
    !!
    !!     call xAxis%set_title("t")
    !!     call yAxis%set_title("x(t)")
    !!     call lgnd%set_is_visible(.true.)
    !!
    !!     call d1%define_data(t, x)
    !!     call d1%set_name("Original")
    !!     call plt%push(d1)
    !!
    !!     call d2%define_data(t, y)
    !!     call d2%set_name("Smoothed")
    !!     call plt%push(d2)
    !!
    !!     call plt%draw()
    !!
    !!     ! Note: The derivative of the signal may also be computed by noting:
    !!     ! d/dt( g * x ) = d/dt(g) * x
    !!     ! So, only the time derivative of the kernel is needed, and then the
    !!     ! convolution proceeds as per normal.
    !! end program
    !! @endcode
    !! The above program produces the following plot using the 
    !! [FPLOT](https://github.com/jchristopherson/fplot) library.
    !! @image html gaussian_filter_example_1.png
    interface convolve
        module procedure :: convolve_1
    end interface

    interface
        module function convolve_1(x, y, method, err) result(rst)
            real(real64), intent(in) :: x(:), y(:)
            integer(int32), intent(in), optional :: method
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_FILTER.F90
! ------------------------------------------------------------------------------
    !> @brief Applies a Gaussian filter to a signal.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function gaussian_filter( &
    !!  real(real64) x(:), &
    !!  real(real64) alpha, &
    !!  integer(int32) k, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in] alpha A parameter that specifies the number of standard
    !!  deviations \f$ \sigma \f$ desired in the kernel.  This parameter is
    !!  related to the standard deviation by \f$ \sigma = \frac{k - 1}{2 \alpha}
    !!  \f$.
    !! @param[in] k The kernel size.  This value must be a positive, non-zero
    !!  integer value less than N.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p k is not within the proper
    !!      bounds.
    !!
    !! @return An N-element array containing the filtered signal.
    interface gaussian_filter
        module procedure :: gaussian_filter_1
    end interface

    !> @brief Applies a total-variation filter to a signal.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function tv_filter( &
    !!  real(real64) x(:),
    !!  real(real64) lambda,
    !!  optional integer(int32) niter,
    !!  optional class(errors), intent(inout) err
    !! )
    !! @endcode
    !!
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in] lambda The regularization parameter.  The actual value to use
    !!  is problem dependent, but the noisier the data, the larger this value
    !!  should be.  A good starting point is typically 0.3 - 0.5; however, the
    !!  actual value is problem dependent.
    !! @param[in] niter An optional parameter controlling the number of 
    !!  iterations performed.  The default limit is 10 iterations.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p niter is less than one.
    !!
    !! @return An N-element array containing the filtered signal.
    interface tv_filter
        module procedure :: filter_tv_1
    end interface

    interface
        module function gaussian_filter_1(x, alpha, k, err) result(rst)
            real(real64), intent(in) :: x(:), alpha
            integer(int32), intent(in) :: k
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function

        module function filter_tv_1(x, lambda, niter, err) result(rst)
            real(real64), intent(in) :: x(:), lambda
            integer(int32), intent(in), optional :: niter
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_OVERLAP.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the number of overlapped signals using a nominal 50% 
    !! overlap.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! integer(int32) pure function compute_overlap_segment_count( &
    !!  integer(int32) n, &
    !!  integer(int32) winsize &
    !! )
    !! @endcode
    !!
    !! @param[in] n The total length of the signal being overlapped.
    !! @param[in] winsize The window size.
    !! @return The number of segments.
    interface compute_overlap_segment_count
        module procedure :: overlap_segment_count_1
    end interface

    !> @brief Extracts a segment from a signal allowing for a nominally 50% 
    !! overlap.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine overlap( &
    !!  real(real64) x(:), &
    !!  integer(int32) seg, &
    !!  integer(int32) winsize, &
    !!  real(real64) buffer(:), &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] x An N-element array containing the entire signal.
    !! @param[in] seg The one-based index of the segment to extract.
    !! @param[in] winsize The size of the window (segment).  If this value is
    !!  less than N, the end of the segment will be padded with zeros.
    !! @param[out] buffer A @p winsize array where the segment will be written.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p winsize is less than one.
    interface overlap
        module procedure :: fill_overlap_buffer_1
    end interface

    interface
        pure module function overlap_segment_count_1(n, winsize) result(rst)
            integer(int32), intent(in) :: n, winsize
            integer(int32) :: rst
        end function
        module subroutine fill_overlap_buffer_1(x, seg, winsize, buffer, err)
            real(real64), intent(in) :: x(:)
            integer(int32), intent(in) :: seg, winsize
            real(real64), intent(out) :: buffer(:)
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface

! ------------------------------------------------------------------------------
end module