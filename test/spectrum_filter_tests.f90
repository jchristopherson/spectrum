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
end module