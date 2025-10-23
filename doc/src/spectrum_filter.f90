module spectrum_filter
    use iso_fortran_env
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
    public :: LOW_PASS_FILTER
    public :: HIGH_PASS_FILTER
    public :: BAND_PASS_FILTER
    public :: BAND_STOP_FILTER

    integer(int32), parameter :: LOW_PASS_FILTER = 1
        !! Denotes a low-pass filter.
    integer(int32), parameter :: HIGH_PASS_FILTER = 2
        !! Denotes a high-pass filter.
    integer(int32), parameter :: BAND_PASS_FILTER = 3
        !! Denotes a band-pass filter.
    integer(int32), parameter :: BAND_STOP_FILTER = 4
        !! Denotes a band-stop filter.


contains
! ******************************************************************************
! GAUSSIAN FILTER
! ------------------------------------------------------------------------------
pure function gaussian_filter(x, alpha, k) result(rst)
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
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

    ! Local Variables
    integer(int32) :: i, kappa, nk
    real(real64) :: sumg
    real(real64), allocatable :: g(:)
    
    ! Initialization
    if (mod(k, 2) == 0) then
        nk = k + 1
    else
        nk = k
    end if
    kappa = -(nk - 1) / 2
    sumg = 0.0d0

    ! Input Checking
    if (nk > size(x) .or. k < 1) return
    if (alpha <= 0.0d0) return

    ! Memory Allocation
    allocate(g(nk))

    ! Define the kernel
    do i = 1, nk
        g(i) = exp(-0.5d0 * (2.0d0 * alpha * kappa / (nk - 1.0d0))**2)
        kappa = kappa + 1
        sumg = sumg + g(i)
    end do

    ! Normalize the kernel to have a sum of one
    g = g / sumg

    ! Compute the convolution and keep only the non-poluted data
    rst = convolve(x, g, SPCTRM_CENTRAL_CONVOLUTION)
end function

! ******************************************************************************
! TOTAL VARIATION FILTERING
! ------------------------------------------------------------------------------
! REF:
! https://eeweb.engineering.nyu.edu/iselesni/lecture_notes/TV_filtering.pdf
pure function tv_filter(x, lambda, niter) result(rst)
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
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

    ! Parameters
    real(real64), parameter :: alpha = 4.0d0

    ! Local Variables
    integer(int32) :: i, k, nit, n
    real(real64) :: t
    real(real64), allocatable :: z(:), work(:), dx(:)
    
    ! Initialization
    if (present(niter)) then
        nit = niter
    else
        nit = 10
    end if
    n = size(x)

    ! Input Checking
    if (nit < 1) return

    ! Memory Allocations
    allocate(rst(n))
    allocate(z(n - 1), source = 0.0d0)
    allocate(dx(n - 1))
    allocate(work(n))

    ! Process
    t = 0.5d0 * lambda
    do i = 1, nit
        work(2:n-1) = difference(z)
        work(1) = z(1)
        work(n) = -z(n - 1)
        rst = x + work

        dx = difference(rst)
        do k = 1, n - 1
            z(k) = z(k) + dx(k) / alpha
            z(k) = max(min(z(k), t), -t)
        end do
    end do
end function

! ------------------------------------------------------------------------------
pure function filter(b, a, x, delays) result(rst)
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
    real(real64), intent(in), optional, target :: delays(:)
        !! An optional array of length MAX(size(a), size(b)) - 1 that provides 
        !! the initial conditions for filter delays, and on ouput, the final 
        !! conditions for filter delays.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

    ! Parameters
    real(real64), parameter :: tol = 2.0d0 * epsilon(2.0d0)

    ! Local Variables
    integer(int32) :: i, m, na, nb, n, nx
    real(real64), allocatable, dimension(:) :: aa, bb, z
    
    ! Initialization
    nx = size(x)
    na = size(a)
    nb = size(b)
    if (na > nb) then
        n = na
    else
        n = nb
    end if

    ! Input Checking
    if (na < 1) return
    if (abs(a(1)) < tol) return

    ! Memory Allocations
    if (present(delays)) then
        ! if (size(delays) /= n - 1) return
        ! zptr(1:n-1) => delays
        allocate(z(n - 1), source = delays)
    else
        ! allocate(zdef(n - 1), source = 0.0d0)
        ! zptr(1:n-1) => zdef
        allocate(z(n - 1), source = 0.0d0)
    end if
    allocate(aa(n), bb(n), rst(nx), source = 0.0d0)

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
            rst(m) = bb(1) * x(m) + z(1)
            do i = 2, n - 1
                z(i-1) = bb(i) * x(m) + z(i) - aa(i) * rst(m)
            end do
            z(n-1) = bb(n) * x(m) - aa(n) * rst(m)
            ! Omit z(n), which is always zero
        end do
    else ! FIR
        do m = 1, nx
            rst(m) = bb(1) * x(m) + z(1)
            do i = 2, n - 1
                z(i-1) = bb(i) * x(m) + z(i)
            end do
            ! Omit z(n), which is always zero
        end do
    end if
end function

! ------------------------------------------------------------------------------
pure function moving_average_filter(navg, x) result(rst)
    !! Applies a moving average filter to a signal.
    integer(int32), intent(in) :: navg
        !! The size of the averaging window.  This parameter must be positive
        !! and non-zero.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the signal to filter.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the filtered signal.

    ! Local Variables
    real(real64) :: a(1), b(navg)

    ! Input Check
    if (navg < 1) return

    ! Process
    a = 1.0d0
    b = 1.0d0 / navg
    rst = filter(b, a, x)
end function

! ******************************************************************************
! V1.1.3 ADDITIONS
! ------------------------------------------------------------------------------
pure function sinc_filter(fc, fs, x, fc2, ftype) result(rst)
    !! Applies a sinc-in-time filter (rectangular frequency response).
    real(real64), intent(in) :: fc
        !! The filter cutoff frequency, in Hz.
    real(real64), intent(in) :: fs
        !! The sampling frequency, in Hz.
    real(real64), intent(in), dimension(:) :: x
        !! The signal to filter.
    real(real64), intent(in), optional :: fc2
        !! The second cutoff frequency for band-pass and band-stop filters.
    integer(int32), intent(in), optional :: ftype
        !! The filter type.  This parameter must be one of the following values:
        !!
        !!  - LOW_PASS_FILTER: Denotes a low-pass filter.
        !!
        !!  - HIGH_PASS_FILTER: Denotes a high-pass filter.
        !!
        !!  - BAND_PASS_FILTER: Denotes a band-pass filter.
        !!
        !!  - BAND_STOP_FILTER: Denotes a band-stop filter.
        !!
        !! The default value is LOW_PASS_FILTER.
    real(real64), allocatable, dimension(:) :: rst
        !! The filtered signal.

    ! Local Variables
    integer(int32) :: start, finish, n, nw, filterType
    real(real64) :: df
    real(real64), allocatable, dimension(:) :: wsave
    
    ! Initialization
    if (present(ftype)) then
        filterType = ftype
    else
        filterType = LOW_PASS_FILTER
    end if
    n = size(x)
    nw = 2 * n + 15

    ! Input Checking
    if (fs <= 0.0d0) then
        return
    end if
    if (fc <= 0.0d0) then
        return
    end if
    if (fc >= 0.5d0 * fs) then
        return
    end if
    if (filterType < LOW_PASS_FILTER .or. filterType > BAND_STOP_FILTER) then
        return
    end if
    if (filterType == BAND_PASS_FILTER .or. filterType == BAND_STOP_FILTER) then
        if (.not. present(fc2)) then
            return
        end if
        if (fc2 <= 0.0d0) then
            return
        end if
        if (fc2 >= 0.5d0 * fs) then
            return
        end if
    end if

    ! Memory Allocations
    allocate(rst(n), source = x)
    allocate(wsave(nw))

    ! Initialize and compute the Fourier transform
    call dffti(n, wsave)
    call dfftf(n, rst, wsave)

    ! Zero out anything above the cutoff frequency
    df = frequency_bin_width(fs, n)
    start = floor(fc / df) * 2
    if (present(fc2)) then
        finish = floor(fc2 / df) * 2
    end if

    ! Apply the appropriate filter
    select case (filterType)
    case (LOW_PASS_FILTER)
        ! Zero out anything above the cutoff frequency
        rst(start:n) = 0.0d0
    case (HIGH_PASS_FILTER)
        ! Zero out anything below the cutoff frequency
        rst(1:start) = 0.0d0
    case (BAND_PASS_FILTER)
        ! Zero out anything below the first cutoff frequency
        rst(1:start) = 0.0d0
        ! Zero out anything above the second cutoff frequency
        rst(finish:n) = 0.0d0
    case (BAND_STOP_FILTER)
        ! Zero out anything between the two cutoff frequencies
        rst(start:finish) = 0.0d0
    end select

    ! Compute the inverse transform to retrieve the filtered signal
    call dfftb(n, rst, wsave)
    rst = rst / n
end function

! ------------------------------------------------------------------------------
end module