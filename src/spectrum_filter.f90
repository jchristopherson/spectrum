submodule (spectrum) spectrum_filter
    implicit none
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
    integer(int32) :: i, kappa, nk, flag
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
        ! call difference(z, work(2:n-1))
        work(2:n-1) = difference(z)
        work(1) = z(1)
        work(n) = -z(n - 1)
        rst = x + work

        ! call difference(rst, dx)
        dx = difference(rst)
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
module function filter_1(b, a, x, delays, err) result(rst)
    ! Arguments
    real(real64), intent(in) :: b(:), a(:), x(:)
    real(real64), intent(inout), optional, target :: delays(:)
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Parameters
    real(real64), parameter :: tol = 2.0d0 * epsilon(2.0d0)

    ! Local Variables
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    integer(int32) :: i, m, na, nb, n, nx, flag
    real(real64), allocatable :: aa(:), bb(:)
    real(real64), allocatable, target :: zdef(:)
    real(real64), pointer :: zptr(:)
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nx = size(x)
    na = size(a)
    nb = size(b)
    if (na > nb) then
        n = na
    else
        n = nb
    end if

    ! Input Checking
    if (na < 1) go to 20
    if (abs(a(1)) < tol) go to 30

    ! Memory Allocations
    if (present(delays)) then
        if (size(delays) /= n - 1) go to 40
        zptr(1:n-1) => delays
    else
        allocate(zdef(n - 1), stat = flag, source = 0.0d0)
        if (flag /= 0) go to 10
        zptr(1:n-1) => zdef
    end if
    allocate(aa(n), stat = flag, source = 0.0d0)
    if (flag == 0) allocate(bb(n), stat = flag, source = 0.0d0)
    if (flag == 0) allocate(rst(nx), stat = flag, source = 0.0d0)
    if (flag /= 0) go to 10

    ! Copy over A & B and scale such that A(1) = 1
    if (abs(a(1) - 1.0d0) > tol) then
        bb(1:nb) = b / a(1)
        aa(1:na) = a / a(1)
    else
        bb(1:nb) = b
        aa(1:na) = a
    end if

    ! Process
    if (na > 1) then ! IIR
        do m = 1, nx
            rst(m) = bb(1) * x(m) + zptr(1)
            do i = 2, n - 1
                zptr(i-1) = bb(i) * x(m) + zptr(i) - aa(i) * rst(m)
            end do
            zptr(n-1) = bb(n) * x(m) - aa(n) * rst(m)
            ! Omit z(n), which is always zero
        end do
    else ! FIR
        do m = 1, nx
            rst(m) = bb(1) * x(m) + zptr(1)
            do i = 2, n - 1
                zptr(i-1) = bb(i) * x(m) + zptr(i)
            end do
            ! Omit z(n), which is always zero
        end do
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("filter_1", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Filter Coefficient Array A Size Error Handler
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The filter coefficient array 'A' must have at " // &
        "least 1 element, but was found to have ", na, " elements."
    call errmgr%report_error("filter_1", trim(errmsg), SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Filter Coefficient Array Value Error Handler
30  continue
    call errmgr%report_error("filter_1", &
        "The 'A(1)' coefficient must be non-zero.", &
        SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Input delays array is not sized correctly
40  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The filter delay array was expected to be of size ", &
        n - 1, "; however, was found to be of size ", size(delays), "."
    call errmgr%report_error("filter_1", trim(errmsg), SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0)
end function

! ------------------------------------------------------------------------------
module function avg_filter_1(navg, x, err) result(rst)
    ! Arguments
    integer(int32), intent(in) :: navg
    real(real64), intent(in) :: x(:)
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Local Variables
    real(real64) :: a(1), b(navg)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if

    ! Input Check
    if (navg < 1) go to 10

    ! Process
    a = 1.0d0
    b = 1.0d0 / navg
    rst = filter(b, a, x, err = errmgr)

    ! End
    return

    ! NAVG invalid value (< 1)
10  continue
    call errmgr%report_error("avg_filter_1", &
        "The averaging window must have at least 1 element.", &
        SPCTRM_INVALID_INPUT_ERROR)
    return
end function

! ------------------------------------------------------------------------------
end submodule