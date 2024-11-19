module spectrum_filter
    use iso_fortran_env
    use ferror
    use spectrum_errors
    use spectrum_convolve
    use spectrum_routines
    use fftpack
    implicit none
    private
    public :: gaussian_filter
    public :: tv_filter
    public :: filter
    public :: moving_average_filter
    public :: sinc_filter

contains
! ******************************************************************************
! GAUSSIAN FILTER
! ------------------------------------------------------------------------------
function gaussian_filter(x, alpha, k, err) result(rst)
    !! Applies a Gaussian filter to a signal.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the signal to filter.
    real(real64), intent(in) :: alpha
        !! A parameter that specifies the number of standard deviations 
        !! \( \sigma \) desired in the kernel.  This parameter is related to the 
        !! standard deviation by \( \sigma = \frac{k - 1}{2 \alpha} \).
    integer(int32), intent(in) :: k
        !! The kernel size.  This value must be a positive, non-zero
        !! integer value less than N.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
        !!      available.
        !! 
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if k is not within the proper
        !!      bounds.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

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
    call errmgr%report_error("gaussian_filter", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Kernel Size Error
20  continue
    call errmgr%report_error("gaussian_filter", "The kernel size must " // &
        "be a positive valued integer less than the size of the signal " // &
        "being filtered.", SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Invalid Kernel Parameter
30  continue
    call errmgr%report_error("gaussian_filter", "The kernal parameter " // &
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
function tv_filter(x, lambda, niter, err) result(rst)
    !! Applies a total-variation filter to a signal.
    !!
    !! The algorithm used by this routine is based upon the algorithm presented 
    !! by [Selesnick and Bayram]
    !! (https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/TV_filtering.pdf).
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the signal to filter.
    real(real64), intent(in) :: lambda
        !! The regularization parameter.  The actual value to use
        !! is problem dependent, but the noisier the data, the larger this value
        !! should be.  A good starting point is typically 0.3 - 0.5; however, the
        !! actual value is problem dependent.
    integer(int32), intent(in), optional :: niter
        !! An optional parameter controlling the number of iterations performed.
        !! The default limit is 10 iterations.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
        !!      available.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if niter is less than one.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

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
    call errmgr%report_error("tv_filter", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Invalid Input Error
20  continue
    call errmgr%report_error("tv_filter", "The number of input " // &
        "iterations must be at least 1.", SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
    100 format(A, I0, A)
    101 format(A)
end function

! ------------------------------------------------------------------------------
function filter(b, a, x, delays, err) result(rst)
    !! Applies the specified filter to a signal.
    !!
    !! The description of the filter in the Z-transform domain is a rational
    !! transfer function of the form:
    !! 
    !! \( Y(z) = \frac{b(1) + b(2) z^{-1} + ... + b(n_b + 1)z^{-n_b}}
    !! {1 + a(2) z^{-1} + ... + a(n_a + 1) z^{-n_a}} X(z) \),
    !! which handles both IIR and FIR filters. The above form assumes a
    !! normalization of a(1) = 1; however, the routine will appropriately 
    !! handle the situation where a(1) is not set to one.
    real(real64), intent(in) :: b(:)
        !! The numerator coefficients of the rational transfer function.
    real(real64), intent(in) :: a(:)
        !! The denominator coefficients of the ration transfer function.  In 
        !! the case of an FIR filter, this parameter should be set to a 
        !! one-element array with a value of one.  Regardless, the value of
        !! a(1) must be non-zero.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the signal to filter.
    real(real64), intent(inout), optional, target :: delays(:)
        !! An optional array of length 
        !! MAX(size(a), size(b)) - 1 that, on input, provides the initial 
        !! conditions for filter delays, and on ouput, the final conditions for 
        !! filter delays.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
        !!      available.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if a(1) is zero.
        !!
        !!  - SPCTRM_ARRAY_SIZE_ERROR: Occurs if a is not sized correctly, or if
        !!      delays is not sized correctly.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

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
    call errmgr%report_error("filter", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Filter Coefficient Array A Size Error Handler
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "The filter coefficient array 'A' must have at " // &
        "least 1 element, but was found to have ", na, " elements."
    call errmgr%report_error("filter", trim(errmsg), SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Filter Coefficient Array Value Error Handler
30  continue
    call errmgr%report_error("filter", &
        "The 'A(1)' coefficient must be non-zero.", &
        SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Input delays array is not sized correctly
40  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The filter delay array was expected to be of size ", &
        n - 1, "; however, was found to be of size ", size(delays), "."
    call errmgr%report_error("filter", trim(errmsg), SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0)
end function

! ------------------------------------------------------------------------------
function moving_average_filter(navg, x, err) result(rst)
    !! Applies a moving average filter to a signal.
    integer(int32), intent(in) :: navg
        !! The size of the averaging window.  This parameter must be positive
        !! and non-zero.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the signal to filter.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_OUT_OF_MEMORY_ERROR: Occurs if there is insufficient memory 
        !!      available.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if navg is less than one.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

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
    call errmgr%report_error("moving_average_filter", &
        "The averaging window must have at least 1 element.", &
        SPCTRM_INVALID_INPUT_ERROR)
    return
end function

! ******************************************************************************
! V1.1.3 ADDITIONS
! ------------------------------------------------------------------------------
function sinc_filter(fc, fs, x, err) result(rst)
    !! Applies a sinc-in-time filter (rectangular frequency response).
    real(real64), intent(in) :: fc
        !! The filter cutoff frequency, in Hz.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    real(real64), intent(in), dimension(:) :: x
        !! The signal to filter.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !! - SPCTRM_INVALID_INPUT_ERROR: Occurs if the cutoff frequency is 
        !!      greater than or equal to the sampling frequency, or if either 
        !!      the cutoff or sampling frequency is zero or negative-valued.
    real(real64), allocatable, dimension(:) :: rst
        !! The filtered signal.

    ! Local Variables
    integer(int32) :: start, n, nw, flag
    real(real64) :: df
    real(real64), allocatable, dimension(:) :: wsave
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    nw = 2 * n + 15

    ! Input Checking
    if (fs <= 0.0d0) then
        call errmgr%report_error("sinc_filter", &
            "The sampling frequency must be non-zero and positive-valued.", &
            SPCTRM_INVALID_INPUT_ERROR)
        return
    end if
    if (fc <= 0.0d0) then
        call errmgr%report_error("sinc_filter", &
            "The cutoff frequency must be non-zero and positive-valued.", &
            SPCTRM_INVALID_INPUT_ERROR)
        return
    end if
    if (fc >= fs) then
        call errmgr%report_error("sinc_filter", &
            "The cutoff frequency must be less than the sampling frequency.", &
            SPCTRM_INVALID_INPUT_ERROR)
        return
    end if

    ! Memory Allocations
    allocate(rst(n), stat = flag, source = x)
    if (flag == 0) allocate(wsave(nw), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("sinc_filter", &
            "Memory allocation error.", SPCTRM_MEMORY_ERROR)
        return
    end if

    ! Initialize and compute the Fourier transform
    call dffti(n, wsave)
    call dfftf(n, rst, wsave)

    ! Zero out anything above the cutoff frequency
    df = frequency_bin_width(fs, n)
    start = floor(fc / df) * 2
    rst(start:n) = 0.0d0

    ! Compute the inverse transform to retrieve the filtered signal
    call dfftb(n, rst, wsave)
    rst = rst / n
end function

! ------------------------------------------------------------------------------
end module