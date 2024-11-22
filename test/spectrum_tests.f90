program test
    use spectrum_psd_tests
    use spectrum_convolution_tests
    use spectrum_diff_tests
    use spectrum_integrate_tests
    use spectrum_filter_tests
    implicit none

    ! Local Variables
    integer(int32) :: flag
    logical :: local

    ! Initialization
    flag = 0

    ! Tests
    local = test_psd()
    if (.not.local) flag = 1

    local = test_periodogram()
    if (.not.local) flag = 2

    local = test_csd()
    if (.not.local) flag = 3

    local = test_convolution()
    if (.not.local) flag = 4

    local = test_spectrogram()
    if (.not.local) flag = 5

    local = test_finite_difference()
    if (.not.local) flag = 6

    local = test_stencil_diff()
    if (.not.local) flag = 7

    local = test_stencil_diff_2()
    if (.not.local) flag = 8

    local = test_integrate()
    if (.not.local) flag = 9

    local = test_sinc_filter()
    if (.not.local) flag = 10

    ! Output
    stop flag
end program