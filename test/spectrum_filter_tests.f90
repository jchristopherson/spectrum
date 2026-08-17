module spectrum_filter_tests
    use iso_fortran_env
    use spectrum
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_sinc_filter() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: n = 1000
        real(real64), parameter :: fs = 1.0d3
        real(real64), parameter :: fc = 1.5d2
        real(real64), parameter :: threshold = 1.0d-8

        ! Local Variables
        integer(int32) :: i, nxfrm
        real(real64) :: df, x(n), xf(n)
        real(real64), allocatable, dimension(:) :: freq, pwr
        type(rectangular_window) :: win

        ! Initialization
        rst = .true.
        call random_number(x)
        win%size = n

        ! Filter the signal
        xf = sinc_filter(fc, fs, x)

        ! Compute the PSD of the filtered signal
        pwr = periodogram(win, xf)
        nxfrm = size(pwr)

        ! Construct the frequency vector
        df = frequency_bin_width(fs, n)
        allocate(freq(nxfrm))
        freq = (/ (i * df, i = 0, nxfrm - 1) /)

        ! Cycle through the frequency array checking the power term
        do i = 2, nxfrm ! no need to look at the DC term, start at 2
            if (freq(i) > fc) then
                if (pwr(i) > threshold) then
                    rst = .false.
                    print "(A)", "TEST FAILED: test_sinc_filter -1"
                    exit
                end if
            end if
        end do
    end function

! ------------------------------------------------------------------------------
function test_filter_frequency_response() result(rst)
    logical :: rst

    real(real64), parameter :: b(2) = [1.0d0, -1.0d0]
    real(real64), parameter :: a(1) = [1.0d0]
    real(real64), parameter :: f(2) = [0.0d0, 0.5d0]
    complex(real64), allocatable :: response(:)

    response = filter_frequency_response(b, a, f, 1.0d0)
    rst = abs(response(1)) < 1.0d-12 .and. &
        abs(abs(response(2)) - 2.0d0) < 1.0d-12
    if (.not.rst) print "(A)", "TEST FAILED: test_filter_frequency_response -1"
end function

! ------------------------------------------------------------------------------
function test_design_iir_filter() result(rst)
    logical :: rst

    real(real64), parameter :: fs = 1.0d3
    real(real64), parameter :: fc = 1.0d2
    real(real64), parameter :: f(3) = [0.0d0, fc, 0.45d0 * fs]
    real(real64), allocatable :: b(:), a(:)
    complex(real64), allocatable :: response(:)

    call design_iir_filter(4_int32, fc, fs, b, a)
    response = filter_frequency_response(b, a, f, fs)
    rst = size(a) == 5 .and. abs(a(1) - 1.0d0) < 1.0d-12 .and. &
        abs(abs(response(1)) - 1.0d0) < 1.0d-12 .and. &
        abs(abs(response(2)) - 1.0d0 / sqrt(2.0d0)) < 1.0d-10 .and. &
        abs(response(3)) < 1.0d-3

    call design_iir_filter(4_int32, fc, fs, b, a, HIGH_PASS_FILTER)
    response = filter_frequency_response(b, a, f, fs)
    rst = rst .and. abs(response(1)) < 1.0d-12 .and. &
        abs(abs(response(2)) - 1.0d0 / sqrt(2.0d0)) < 1.0d-10 .and. &
        abs(abs(response(3)) - 1.0d0) < 1.0d-3
    if (.not.rst) print "(A)", "TEST FAILED: test_design_iir_filter -1"
end function

! ------------------------------------------------------------------------------
function test_butterworth_filter_order() result(rst)
    logical :: rst

    integer(int32) :: lowpass_order, highpass_order, invalid_order

    lowpass_order = butterworth_filter_order(100.0d0, 1000.0d0, 200.0d0, &
        1.0d0, 40.0d0)
    highpass_order = butterworth_filter_order(200.0d0, 1000.0d0, 100.0d0, &
        1.0d0, 40.0d0, HIGH_PASS_FILTER)
    invalid_order = butterworth_filter_order(200.0d0, 1000.0d0, 100.0d0, &
        1.0d0, 40.0d0)
    rst = lowpass_order == 7 .and. highpass_order == 7 .and. invalid_order == 0
    if (.not.rst) print "(A)", "TEST FAILED: test_butterworth_filter_order -1"
end function

! ------------------------------------------------------------------------------
function test_design_fir_filter() result(rst)
    logical :: rst

    real(real64), allocatable :: b(:), a(:)
    type(rectangular_window) :: win

    call design_fir_filter(5_int32, 100.0d0, 1000.0d0, b, a, win = win)
    rst = size(b) == 5 .and. size(a) == 1 .and. win%size == 5 .and. &
        abs(sum(b) - 1.0d0) < 1.0d-12 .and. abs(b(1) - b(5)) < 1.0d-12
    if (.not.rst) print "(A)", "TEST FAILED: test_design_fir_filter -1"
end function

! ------------------------------------------------------------------------------
end module