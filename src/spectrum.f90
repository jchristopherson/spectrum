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
    public :: compute_transform_length
    public :: frequency_bin_width
    public :: is_power_of_two
    public :: next_power_of_two
    public :: psd
    public :: csd
    public :: SPCTRM_MEMORY_ERROR
    public :: SPCTRM_INVALID_INPUT_ERROR
    public :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR

! ******************************************************************************
! CONSTANTS
! ------------------------------------------------------------------------------
    integer(int32), parameter :: SPCTRM_MEMORY_ERROR = 10000
    integer(int32), parameter :: SPCTRM_INVALID_INPUT_ERROR = 10001
    integer(int32), parameter :: SPCTRM_ARRAY_SIZE_MISMATCH_ERROR = 10002

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

! ------------------------------------------------------------------------------
end module