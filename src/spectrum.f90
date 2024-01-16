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
    public :: filter
    public :: moving_average_filter
    public :: compute_overlap_segment_count
    public :: overlap
    public :: unwrap
    public :: siso_transfer_function
    public :: periodogram
    public :: cross_periodogram
    public :: SPCTRM_MEMORY_ERROR
    public :: SPCTRM_INVALID_INPUT_ERROR
    public :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR
    public :: SPCTRM_ARRAY_SIZE_ERROR
    public :: SPCTRM_FULL_CONVOLUTION
    public :: SPCTRM_CENTRAL_CONVOLUTION
    public :: SPCTRM_H1_ESTIMATOR
    public :: SPCTRM_H2_ESTIMATOR

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: SPCTRM_MEMORY_ERROR = 10000
    integer(int32), parameter :: SPCTRM_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR = 10002
    integer(int32), parameter :: SPCTRM_ARRAY_SIZE_ERROR = 10003

    !> A flag for requesting a full convolution.
    integer(int32), parameter :: SPCTRM_FULL_CONVOLUTION = 50000
    !> A flag for requesting the central portion of the convolution that is the
    !! same length as the input signal.
    integer(int32), parameter :: SPCTRM_CENTRAL_CONVOLUTION = 50001
    !> A flag for requesting an H1 transfer function estimator.  An H1 estimator
    !! is best used when noise is uncorrelated with the input, and results in
    !! a transfer function estimate of
    !! @par
    !! \f$ H_{1} = \frac{P_{yx}}{P_{xx}} \f$.
    integer(int32), parameter :: SPCTRM_H1_ESTIMATOR = 50002
    !> A flag for requesting an H2 transfer function estimator.  An H2 estimator
    !! is best used when noise is uncorrelated with the output, and results in
    !! a transfer function estimate of
    !! @par
    !! \f$ H_{2} = \frac{P_{yy}}{P_{xy}} \f$.
    integer(int32), parameter :: SPCTRM_H2_ESTIMATOR = 50003

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

    !> @brief Defines a rectangular window.
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

    !> @brief Shifts phase angle arrays to deal with jumps greater than or equal
    !! to @p tol by adding multiples of +/- 2 \f$ \pi \f$ until
    !! the jump is less than @p tol.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! subroutine unwrap( &
    !!  real(real64) x(:), &
    !!  optional real(real64) tol &
    !! )
    !! @endcode
    !!
    !! @param[in] x On input, the phase angle array.  On ouptut, the modified
    !!  phase angle array.
    !! @param[in] tol An optional input that controls the tolerated jump size.
    !!  The default value is \f$ \pi \f$.
    interface unwrap
        module procedure :: unwrap_1
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

        module subroutine unpack_real_transform(x, cx, fac)
            real(real64), intent(in) :: x(:)
            complex(real64), intent(out) :: cx(:)
            real(real64), intent(in), optional :: fac
        end subroutine

        module subroutine unwrap_1(x, tol)
            real(real64), intent(inout) :: x(:)
            real(real64), intent(in), optional :: tol
        end subroutine
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
    !!  optional integer(int32) nfft, &
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
        module function psd_welch(win, x, fs, nfft, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:)
            real(real64), intent(in), optional :: fs
            integer(int32), intent(in), optional :: nfft
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_CSD.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the cross spectral density (CSD) of a signal via Welch's
    !! averaged method (sometimes referred to as cross power spectral density
    !! or CPSD).
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function csd( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  real(real64) y(:), &
    !!  optional real(real64) fs, &
    !!  integer(int32) nfft, &
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
    !!  - SPCTRM_ARRAY_SIZE_MISMATCH_ERROR: Occurs if @p x and @p y are not the
    !!      same size.
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
    !! - [Wikipedia - Cross Power Spectral Density](https://en.wikipedia.org/wiki/Spectral_density#Cross-spectral_density)
    interface csd
        module procedure :: csd_welch
    end interface

    interface
        module function csd_welch(win, x, y, fs, nfft, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:), y(:)
            real(real64), intent(in), optional :: fs
            integer(int32), intent(in), optional :: nfft
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable :: rst(:)
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
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The window to apply.
    !! @param[in] x The signal to analyze.  The signal must be longer than the
    !!  size of the window @p win.
    !! @param[out] offsets An optional allocatable array that, if supplied, will
    !!  be filled with the starting indices of each window segment.
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
    !! @par References
    !! - [Wikipedia - Short Time Fourier Transform](https://en.wikipedia.org/wiki/Short-time_Fourier_transform)
    interface spectrogram
        module procedure :: stft
    end interface

    interface
        module function stft(win, x, offsets, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:)
            integer(int32), intent(out), optional, allocatable :: offsets(:)
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
    !!
    !! @par Remarks
    !! The algorithm used by this routine is based upon the algorithm presented 
    !! by [Selesnick and Bayram]
    !! (https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/TV_filtering.pdf).
    interface tv_filter
        module procedure :: filter_tv_1
    end interface

    !> @brief Applies the specified filter to a signal.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function filter( &
    !!  real(real64) b(:), &
    !!  real(real64) a(:), &
    !!  real(real64) x(:), &
    !!  optional allocatable real(real64) delays(:), &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] b The numerator coefficients of the rational transfer 
    !!  function.
    !! @param[in] a The denominator coefficients of the ration transfer 
    !!  function.  In the case of an FIR filter, this parameter should be set
    !!  to a one-element array with a value of one.  Regardless, the value of
    !!  a(1) must be non-zero.
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in,out] delays An optional array of length 
    !!  MAX(size(@p a), size(@p b)) - 1 that, on input, provides the initial 
    !!  conditions for filter delays, and on ouput, the final conditions for 
    !!  filter delays.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if a(1) is zero.
    !!  - SPCTRM_ARRAY_SIZE_ERROR: Occurs if @p is not sized correctly, or if
    !!      @p delays is not sized correctly.
    !!
    !! @return An N-element array containing the filtered signal.
    !!
    !! @par Remarks
    !! The description of the filter in the Z-transform domain is a rational
    !! transfer function of the form:
    !! @par
    !! \f$ Y(z) = \frac{b(1) + b(2) z^{-1} + ... + b(n_b + 1)z^{-n_b}}
    !! {1 + a(2) z^{-1} + ... + a(n_a + 1) z^{-n_a}} X(z) \f$,
    !! which handles both IIR and FIR filters. The above form assumes a
    !! normalization of a(1) = 1; however, the routine will appropriately 
    !! handle the situation where a(1) is not set to one.
    interface filter
        module procedure :: filter_1
    end interface

    !> @brief Applies a moving average filter to a signal.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function moving_average_filter( &
    !!  integer(int32), intent(in) :: navg, &
    !!  real(real64), intent(in) :: x(:), &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] navg The size of the averaging window.  This parameter must
    !!  be positive and non-zero.
    !! @param[in] x An N-element array containing the signal to filter.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p navg is less than one.
    interface moving_average_filter
        module procedure :: avg_filter_1
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

        module function filter_1(b, a, x, delays, err) result(rst)
            real(real64), intent(in) :: b(:), a(:), x(:)
            real(real64), intent(inout), optional, target :: delays(:)
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function

        module function avg_filter_1(navg, x, err) result(rst)
            integer(int32), intent(in) :: navg
            real(real64), intent(in) :: x(:)
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

! ******************************************************************************
! SPECTRUM_TF.F90
! ------------------------------------------------------------------------------
    !> @brief Estimates the transfer function for a single-input/single-output
    !! (SISO) system.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function siso_transfer_function( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  real(real64) y(:), &
    !!  optional integer(int32) xfrmtype, &
    !!  class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The window object.
    !! @param[in] x An N-element array containing the input signal.
    !! @param[in] y An N-element array containing the output signal.
    !! @param[in] etype An optional input that, if supplied, denotes the
    !!  estimator to use.  If no value is specified, an H1 estimator is used.
    !!  The following options are supported.
    !!  - @ref SPCTRM_H1_ESTIMATOR: Uses an H1 estimate.
    !!  - @ref SPCTRM_H2_ESTIMATOR: Uses an H2 estimate.
    !!  If an unrecognized value is provided, the routine defaults to an H1 
    !!  estimator.
    !! @param[in,out] err
    !!
    !! @return Returns the complex-valued transfer function estimate.
    interface siso_transfer_function
        module procedure :: siso_xfrm
    end interface

    interface
        module function siso_xfrm(win, x, y, etype, nfft, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:), y(:)
            integer(int32), intent(in), optional :: etype
            integer(int32), intent(in), optional :: nfft
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable :: rst(:)
        end function
    end interface

! ******************************************************************************
! SPECTRUM_PERIODOGRAM.F90
! ------------------------------------------------------------------------------
    !> @brief Computes the periodogram of a signal.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function periodogram( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  optional real(real64) fs, &
    !!  optional integer(int32) nfft, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The @ref window to apply.  The size of the window must be
    !!  non-zero and positive-valued.
    !! @param[in] x The signal to transform.  The length of this array must be
    !!  the same size as the window @p win.
    !! @param[in] fs An optional input, that if supplied, allows for 
    !!  normalization of the computed spectrum by the frequency resolution.
    !! @param[in] nfft An optional input, that if supplied, allows for zero
    !!  padding the FFT to the desired length.  If not supplied, this parameter
    !!  defaults to the window size.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p win is not sized 
    !!      appropriately.
    !!  - SPCTRM_ARRAY_SIZE_MISMATCH_ERROR: Occurs if @p x and @p y are not the 
    !!      same size as the window @p win.
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
    interface periodogram
        module procedure :: periodogram_1
    end interface

    !> @brief Computes the cross-spectral periodogram of two signals.
    !!
    !! @par Syntax
    !! @code{.f90}
    !! allocatable real(real64)(:) function csd( &
    !!  class(window) win, &
    !!  real(real64) x(:), &
    !!  real(real64) y(:), &
    !!  optional real(real64) fs, &
    !!  optional integer(int32) nfft, &
    !!  optional class(errors) err &
    !! )
    !! @endcode
    !!
    !! @param[in] win The @ref window to apply.  The size of the window must be
    !!  non-zero and positive-valued.
    !! @param[in] x The first N-element signal to transform.
    !! @param[in] y The second N-element signal to transform.
    !! @param[in] fs An optional input, that if supplied, allows for 
    !!  normalization of the computed spectrum by the frequency resolution.
    !! @param[in] nfft An optional input, that if supplied, allows for zero
    !!  padding the FFT to the desired length.  If not supplied, this parameter
    !!  defaults to the window size.
    !! @param[in,out] err An optional errors-based object that if provided can
    !!  be used to retrieve information relating to any errors encountered 
    !!  during execution.  If not provided, a default implementation of the 
    !!  errors class is used internally to provide error handling.  Possible 
    !!  errors and warning messages that may be encountered are as follows.
    !!  - SPCTRM_MEMORY_ERROR: Occurs if there is insufficient memory 
    !!      available.
    !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if @p win is not sized 
    !!      appropriately.
    !!  - SPCTRM_ARRAY_SIZE_MISMATCH_ERROR: Occurs if @p x and @p y are not the
    !!      same size, or if @p x and @p y are not the same size as the window
    !!      @p win.
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
    interface cross_periodogram
        module procedure :: cross_periodogram_1
    end interface

    interface
        module function periodogram_1(win, x, fs, nfft, err) result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:)
            real(real64), intent(in), optional :: fs
            integer(int32), intent(in), optional :: nfft
            class(errors), intent(inout), optional, target :: err
            real(real64), allocatable :: rst(:)
        end function

        module function cross_periodogram_1(win, x, y, fs, nfft, err) &
            result(rst)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:), y(:)
            real(real64), intent(in), optional :: fs
            integer(int32), intent(in), optional :: nfft
            class(errors), intent(inout), optional, target :: err
            complex(real64), allocatable :: rst(:)
        end function
        
        module subroutine periodogram_driver(win, x, xfrm, fs, nfft, work, &
            initxfrm, cwork, err)
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:)
            real(real64), intent(out) :: xfrm(:)
            real(real64), intent(in), optional :: fs
            integer(int32), intent(in), optional :: nfft
            real(real64), intent(out), optional, target :: work(:)
            logical, intent(in), optional :: initxfrm
            complex(real64), intent(out), optional, target :: cwork(:)
            class(errors), intent(inout), optional, target :: err
        end subroutine

        module subroutine cross_periodogram_driver(win, x, y, xfrm, fs, nfft, &
            work, initxfrm, cwork, err)
            ! Arguments
            class(window), intent(in) :: win
            real(real64), intent(in) :: x(:), y(:)
            complex(real64), intent(out) :: xfrm(:)
            real(real64), intent(in), optional :: fs
            integer(int32), intent(in), optional :: nfft
            real(real64), intent(out), optional, target :: work(:)
            logical, intent(in), optional :: initxfrm
            complex(real64), intent(out), optional, target :: cwork(:)
            class(errors), intent(inout), optional, target :: err
        end subroutine
    end interface


! ------------------------------------------------------------------------------
end module