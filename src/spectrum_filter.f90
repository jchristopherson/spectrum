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

! ******************************************************************************
! TOTAL VARIATION FILTERING
! ------------------------------------------------------------------------------
! REF:
! https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/TV_filtering.pdf
module function filter_tv_1(x, lambda, niter, err) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:), lambda
    integer(int32), intent(in), optional :: niter
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Parameters
    real(real64), parameter :: alpha = 4.0d0

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    integer(int32) :: i, k, nit, n, flag
    real(real64) :: t
    real(real64), allocatable :: z(:), work(:), dx(:)
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(niter)) then
        nit = niter
    else
        nit = 10
    end if
    n = size(x)

    ! Input Checking
    if (nit < 1) go to 20

    ! Memory Allocations
    allocate(rst(n), stat = flag)
    if (flag == 0) allocate(z(n - 1), stat = flag, source = 0.0d0)
    if (flag == 0) allocate(dx(n - 1), stat = flag)
    if (flag == 0) allocate(work(n), stat = flag)
    if (flag /= 0) go to 10

    ! Process
    t = 0.5d0 * lambda
    do i = 1, nit
        call difference(z, work(2:n-1))
        work(1) = z(1)
        work(n) = -z(n - 1)
        rst = x + work

        call difference(rst, dx)
        do k = 1, n - 1
            z(k) = z(k) + dx(k) / alpha
            z(k) = max(min(z(k), t), -t)
        end do
    end do

    ! End
    return

    ! Memory Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("filter_tv_1", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Invalid Input Error
20  continue
    call errmgr%report_error("filter_tv_1", "The number of input " // &
        "iterations must be at least 1.", SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
    100 format(A, I0, A)
    101 format(A)
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------

! ******************************************************************************
! HELPER ROUTINES
! ------------------------------------------------------------------------------
! Computes the difference between elements in an array.
!
! - x: An N-element array.
! - dx: The N-1 element results array.
subroutine difference(x, dx)
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64), intent(out) :: dx(:)

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    do i = 1, n - 1
        dx(i) = x(i+1) - x(i)
    end do
end subroutine

! ------------------------------------------------------------------------------
end submodule