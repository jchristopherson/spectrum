program test
    use spectrum_psd_tests
    implicit none

    ! Local Variables
    integer(int32) :: flag
    logical :: local

    ! Initialization
    flag = 0

    ! Tests
    local = test_psd()
    if (.not.local) flag = 1

    local = test_csd()
    if (.not.local) flag = 2

    ! Output
    stop flag
end program