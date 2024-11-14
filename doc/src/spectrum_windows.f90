module spectrum_windows
    use iso_fortran_env
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

    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    
    type, abstract :: window
        !! Defines the structure of a window.
        integer(int32), public :: size = 0
            !! The window size.
    contains
        procedure(window_function), public, deferred, pass :: evaluate
    end type

    interface
        pure function window_function(this, bin) result(rst)
            !! Evaluates the window function.
            use iso_fortran_env, only : real64, int32
            import window
            class(window), intent(in) :: this
                !! The window object.
            integer(int32), intent(in) :: bin
                !! The index or bin number [0, n], where n is the window size.
            real(real64) :: rst
                !! The function value.
        end function
    end interface

    !> @brief Defines a rectangular window.
    type, extends(window) :: rectangular_window
        !! Defines a rectangular window.
    contains
        procedure, public :: evaluate => rw_eval
    end type

    type, extends(window) :: hann_window
        !! Defines a Hann window.
        !!
        !! $$ w(j) = \frac{1}{2} \left( 1 - \cos \left( \frac{2 \pi j}{n} \right) \right) $$.
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    contains
        procedure, public :: evaluate => hann_eval
    end type

    type, extends(window) :: hamming_window
        !! Defines a Hamming window.
        !!
        !! $$ w(j) = 0.54 - 0.46 \cos \left( \frac{2 \pi j}{n} \right) $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    contains
        procedure, public :: evaluate => hamming_eval
    end type

    type, extends(window) :: welch_window
        !! Defines a Welch window.
        !!
        !! $$ w(j) = 1 - \left( \frac{j - \frac{n}{2} }{ \frac{n}{2} } \right)^2 $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    contains
        procedure, public :: evaluate => welch_eval
    end type

    type, extends(window) :: blackman_harris_window
        !! Defines a Blackman-Harris window.
        !!
        !! $$ w(j) = 0.3635819 - 0.4891775 \cos \left( \frac{2 \pi j}{n} \right)
        !! + 0.1365995 \cos \left( \frac{4 \pi j}{n} \right) - 0.0106411
        !! \cos \left( \frac{6 \pi j}{n} \right) $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    contains
        procedure, public :: evaluate => bhw_eval
    end type

    type, extends(window) :: flat_top_window
        !! Defines a flat-top window.
        !!
        !! $$ w(j) = 0.21557895 - 
        !! 0.41663158 \cos \left( \frac{2 \pi j}{N} \right) +
        !! 0.277263158 \cos \left( \frac{4 \pi j}{N} \right) -
        !! 0.083578947 \cos \left( \frac{6 \pi j}{N} \right)  +
        !! 0.006947368 \cos \left( \frac{8 \pi j}{N} \right) $$
        !!
        !! See Also
        !!
        !! - [Wikipedia](https://en.wikipedia.org/wiki/Window_function)
    contains
        procedure, public :: evaluate => ftw_eval
    end type

contains
! ------------------------------------------------------------------------------
pure function rw_eval(this, bin) result(rst)
    !! Evaluates the window function.
    class(rectangular_window), intent(in) :: this
        !! The rectangular_window object.
    integer(int32), intent(in) :: bin
        !! The index or bin number [0, n], where n is the window size.
    real(real64) :: rst
        !! The function value.
    rst = 1.0d0
end function

! ------------------------------------------------------------------------------
pure function hann_eval(this, bin) result(rst)
    !! Evaluates the window function.
    class(hann_window), intent(in) :: this
        !! The hann_window object.
    integer(int32), intent(in) :: bin
        !! The index or bin number [0, n], where n is the window size.
    real(real64) :: rst
        !! The function value.
    rst = 0.5d0 * (1.0d0 - cos(2.0d0 * pi * bin / this%size))
end function

! ------------------------------------------------------------------------------
pure function hamming_eval(this, bin) result(rst)
    !! Evaluates the window function.
    class(hamming_window), intent(in) :: this
        !! The hamming_window object.
    integer(int32), intent(in) :: bin
        !! The index or bin number [0, n], where n is the window size.
    real(real64) :: rst
        !! The function value.
    rst = 0.54d0 - 0.46d0 * cos(2.0d0 * pi * bin / this%size)
end function

! ------------------------------------------------------------------------------
pure function welch_eval(this, bin) result(rst)
    !! Evaluates the window function.
    class(welch_window), intent(in) :: this
        !! The welch_window object.
    integer(int32), intent(in) :: bin
        !! The index or bin number [0, n], where n is the window size.
    real(real64) :: rst
        !! The function value.

    real(real64) :: halfsize
    halfsize = 0.5d0 * this%size
    rst = 1.0d0 - ((bin - halfsize) / halfsize)**2
end function

! ------------------------------------------------------------------------------
pure function bhw_eval(this, bin) result(rst)
    !! Evaluates the window function.
    class(blackman_harris_window), intent(in) :: this
        !! The blackman_harris_window object.
    integer(int32), intent(in) :: bin
        !! The index or bin number [0, n], where n is the window size.
    real(real64) :: rst
        !! The function value.
    rst = 0.35875d0 - &
        0.48829d0 * cos(2.0d0 * pi * bin / this%size) + &
        0.14128d0 * cos(4.0d0 * pi * bin / this%size) - &
        0.01168d0 * cos(6.0d0 * pi * bin / this%size)
end function

! ------------------------------------------------------------------------------
pure function ftw_eval(this, bin) result(rst)
    !! Evaluates the window function.
    class(flat_top_window), intent(in) :: this
        !! The flat_top_window object.
    integer(int32), intent(in) :: bin
        !! The index or bin number [0, n], where n is the window size.
    real(real64) :: rst
        !! The function value.
    
    real(real64), parameter :: a0 = 0.21557895d0
    real(real64), parameter :: a1 = 0.41663158d0
    real(real64), parameter :: a2 = 0.277263158d0
    real(real64), parameter :: a3 = 8.3578947d-2
    real(real64), parameter :: a4 = 6.947368d-3

    real(real64) :: arg
    arg = pi * bin / this%size

    rst = a0 - a1 * cos(2.0d0 * arg) + a2 * cos(4.0d0 * arg) - &
        a3 * cos(6.0d0 * arg) + a4 * cos(8.0d0 * arg)
end function

! ------------------------------------------------------------------------------
end module