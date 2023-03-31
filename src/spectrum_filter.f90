submodule (spectrum) spectrum_filter
contains
! ******************************************************************************
! GAUSSIAN FILTER
! ------------------------------------------------------------------------------
module function gaussian_filter_1(x, alpha, k, err) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:), alpha
    integer(int32), intent(in) :: k
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: kappa, nk, flag
    real(real64) :: sumg
    real(real64), allocatable :: g(:)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (mod(k, 2) == 0) then
        nk = k + 1
    else
        nk = k
    end if
    kappa = -(nk - 1) / 2
    sumg = 0.0d0

    ! Input Checking
    if (nk > size(x) .or. k < 1) go to 20
    if (alpha <= 0.0d0) go to 30

    ! Memory Allocation
    allocate(g(nk), stat = flag)
    if (flag /= 0) go to 10

    ! Define the kernel
    do i = 1, nk
        g(i) = exp(-0.5d0 * (2.0d0 * alpha * kappa / (nk - 1.0d0))**2)
        kappa = kappa + 1
        sumg = sumg + g(i)
    end do

    ! Normalize the kernel to have a sum of one
    g = g / sumg

    ! Compute the convolution and keep only the non-poluted data
    rst = convolve(x, g, SPCTRM_CENTRAL_CONVOLUTION, err = errmgr)

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("gaussian_filter_1", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Kernel Size Error
20  continue
    call errmgr%report_error("gaussian_filter_1", "The kernel size must " // &
        "be a positive valued integer less than the size of the signal " // &
        "being filtered.", SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Invalid Kernel Parameter
30  continue
    call errmgr%report_error("gaussian_filter_1", "The kernal parameter " // &
        "alpha must be positive-valued.", SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule