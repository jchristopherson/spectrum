program test
    use spectrum_psd_tests
    use spectrum_convolution_tests
    use spectrum_diff_tests
    use spectrum_integrate_tests
    use spectrum_filter_tests
    use spectrum_tf_tests
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

    local = test_spectral_endpoint_scaling()
    if (.not.local) flag = 4

    local = test_convolution()
    if (.not.local) flag = 5

    local = test_spectrogram()
    if (.not.local) flag = 6

    local = test_stft_scaling()
    if (.not.local) flag = 7

    local = test_finite_difference()
    if (.not.local) flag = 8

    local = test_stencil_diff()
    if (.not.local) flag = 9

    local = test_stencil_diff_2()
    if (.not.local) flag = 10

    local = test_tvr_derivative_sparse()
    if (.not.local) flag = 11

    local = test_integrate()
    if (.not.local) flag = 12

    local = test_sinc_filter()
    if (.not.local) flag = 13

    local = test_filter_frequency_response()
    if (.not.local) flag = 13

    local = test_design_iir_filter()
    if (.not.local) flag = 14

    local = test_butterworth_filter_order()
    if (.not.local) flag = 15

    local = test_design_fir_filter()
    if (.not.local) flag = 16

    local = test_siso_transfer_function()
    if (.not.local) flag = 17

    local = test_mimo_transfer_function()
    if (.not.local) flag = 18

    ! Output
    stop flag
end program