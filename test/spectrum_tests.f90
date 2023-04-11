program test
    use spectrum_psd_tests
    use spectrum_convolution_tests
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

    ! Output
    stop flag
end program