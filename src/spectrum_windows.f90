submodule (spectrum) spectrum_windows
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
contains
! ------------------------------------------------------------------------------
pure module function rw_eval(this, bin) result(rst)
    class(rectangular_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    rst = 1.0d0
end function

! ------------------------------------------------------------------------------
pure module function hann_eval(this, bin) result(rst)
    class(hann_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    rst = 0.5d0 * (1.0d0 - cos(2.0d0 * pi * bin / this%size))
end function

! ------------------------------------------------------------------------------
pure module function hamming_eval(this, bin) result(rst)
    class(hamming_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    rst = 0.54d0 - 0.46d0 * cos(2.0d0 * pi * bin / this%size)
end function

! ------------------------------------------------------------------------------
pure module function welch_eval(this, bin) result(rst)
    class(welch_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst, halfsize
    halfsize = 0.5d0 * this%size
    rst = 1.0d0 - ((bin - halfsize) / halfsize)**2
end function

! ------------------------------------------------------------------------------
pure module function bhw_eval(this, bin) result(rst)
    class(blackman_harris_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    rst = 0.35875d0 - &
        0.48829d0 * cos(2.0d0 * pi * bin / this%size) + &
        0.14128d0 * cos(4.0d0 * pi * bin / this%size) - &
        0.01168d0 * cos(6.0d0 * pi * bin / this%size)
end function

! ------------------------------------------------------------------------------
pure module function ftw_eval(this, bin) result(rst)
    class(flat_top_window), intent(in) :: this
    integer(int32), intent(in) :: bin
    real(real64) :: rst
    
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
end submodule