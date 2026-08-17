module spectrum_convolve
    use iso_fortran_env
    use fftpack
    implicit none
    private
    public :: convolve
    public :: deconvolve
    public :: SPCTRM_FULL_CONVOLUTION
    public :: SPCTRM_CENTRAL_CONVOLUTION

    integer(int32), parameter :: SPCTRM_FULL_CONVOLUTION = 50000
        !! A flag for requesting a full convolution.
    integer(int32), parameter :: SPCTRM_CENTRAL_CONVOLUTION = 50001
        !! A flag for requesting the central portion of the convolution that is 
        !! the same length as the input signal.

contains
! ------------------------------------------------------------------------------
pure function convolve(x, y, method) result(rst)
    !! Computes the convolution of a signal and kernel.
    real(real64), intent(in) :: x(:)
        !! The N-element signal.
    real(real64), intent(in) :: y(:)
        !! The M-element kernel.
    integer(int32), intent(in), optional :: method
        !! An optional input that dictates the expected
        !! convolution result.  The following options are available.
        !!
        !! - SPCTRM_FULL_CONVOLUTION: The full convolution results are provided, 
        !!     including the portions polluted courtesy of the zero-padding and
        !!     the corresponding wrap-around effects.  The length of this output
        !!     is N + M - 1.
        !!
        !!   SPCTRM_CENTRAL_CONVOLUTION: The N-element result containing the 
        !!     convolved signal not poluted by the zero-padding and 
        !!     corresponding wrap-around effects.
    real(real64), allocatable :: rst(:)
        !! The convolved result.

    ! Local Variables
    integer(int32) :: mth, n1, n2
    real(real64), allocatable :: temp(:)
    complex(real64), allocatable :: xmod(:), ymod(:)
    
    ! Initialization
    if (present(method)) then
        mth = method
    else
        mth = SPCTRM_FULL_CONVOLUTION
    end if

    ! Input Checking
    if (mth /= SPCTRM_FULL_CONVOLUTION .and. &
        mth /= SPCTRM_CENTRAL_CONVOLUTION)  &
    then
        ! Reset
        mth = SPCTRM_FULL_CONVOLUTION
    end if

    ! Prepare the signals
    if (size(x) < 1 .or. size(y) < 1) return
    call prepare_conv(x, y, xmod, ymod, n1, n2)

    ! Compute the convolution
    if (mth == SPCTRM_CENTRAL_CONVOLUTION) then
        call convolve_driver(xmod, ymod, 1, temp)
        rst = temp(n1:n2)
    else
        call convolve_driver(xmod, ymod, 1, rst)
    end if
end function

! ------------------------------------------------------------------------------
pure function deconvolve(x, y, method) result(rst)
    !! Computes the frequency-domain deconvolution of a signal and kernel.
    !!
    !! The result uses the same output selection as [[convolve]]. Full
    !! deconvolution returns N + M - 1 samples; central deconvolution returns
    !! the central portion with the length of the longer input.
    real(real64), intent(in) :: x(:)
        !! The signal to deconvolve.
    real(real64), intent(in) :: y(:)
        !! The kernel to remove from the signal.
    integer(int32), intent(in), optional :: method
        !! SPCTRM_FULL_CONVOLUTION (default) or SPCTRM_CENTRAL_CONVOLUTION.
    real(real64), allocatable :: rst(:)
        !! The deconvolved result.

    integer(int32) :: mth, n1, n2
    real(real64), allocatable :: temp(:)
    complex(real64), allocatable :: xmod(:), ymod(:)

    if (present(method)) then
        mth = method
    else
        mth = SPCTRM_FULL_CONVOLUTION
    end if
    if (mth /= SPCTRM_FULL_CONVOLUTION .and. &
        mth /= SPCTRM_CENTRAL_CONVOLUTION) mth = SPCTRM_FULL_CONVOLUTION
    if (size(x) < 1 .or. size(y) < 1) return

    call prepare_conv(x, y, xmod, ymod, n1, n2)
    if (mth == SPCTRM_CENTRAL_CONVOLUTION) then
        call convolve_driver(xmod, ymod, -1, temp)
        rst = temp(n1:n2)
    else
        call convolve_driver(xmod, ymod, -1, rst)
    end if
end function

! ------------------------------------------------------------------------------
! Prepares two signals for convolution.
!
! - x: The original signal of length N.
! - y: The original response signal of length M where M <= N.
! - xmod: The padded X signal of length M + N - 1.
! - ymod: The padded Y signal of length M + N - 1.
! - n2: The ending index of the unspoiled data
pure subroutine prepare_conv(x, y, xmod, ymod, n1, n2)
    ! Arguments
    real(real64), intent(in) :: x(:), y(:)
    complex(real64), intent(out), allocatable :: xmod(:), ymod(:)
    integer(int32), intent(out) :: n1, n2

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: nx, ny, n, nbig, dn, m
    logical :: evenY
    ! Initialization
    nx = size(x)
    ny = size(y)
    n = nx + ny - 1

    ! Memory Allocation
    allocate(xmod(n), source = zero)
    allocate(ymod(n), source = zero)

    ! Copy over X & Y
    if (nx >= ny) then
        xmod(1:nx) = cmplx(x, 0.0d0, real64)
        ymod(1:ny) = cmplx(y, 0.0d0, real64)
        nbig = nx
    else
        xmod(1:ny) = cmplx(y, 0.0d0, real64)
        ymod(1:nx) = cmplx(x, 0.0d0, real64)
        nbig = ny
    end if

    ! Locate the unspoiled data points
    dn = n - nbig
    if (mod(dn, 2) == 0) then
        n1 = dn / 2 + 1
    else
        n1 = (dn + 1) / 2 + 1
    end if
    n2 = n1 + nbig - 1
end subroutine

! ------------------------------------------------------------------------------
! Computes the convolution or deconvolution of two signals.
!
! - x: The zero-padded original signal.
! - y: The padded and wrapped response signal.
! - method: +1 for convolution, -1 for deconvolution.
! - rst: The convolved signal.
pure subroutine convolve_driver(x, y, method, rst)
    ! Arguments
    complex(real64), intent(inout) :: x(:), y(:)
    integer(int32), intent(in) :: method
    real(real64), intent(out), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: i, n, lwork, nend
    real(real64) :: fac
    real(real64), allocatable, dimension(:) :: work
    complex(real64), allocatable, dimension(:) :: cxy
    ! Generate the workspace array for the transforms & allocate memory
    n = size(x)
    lwork = 4 * n + 15
    allocate(work(lwork))
    allocate(rst(n))
    allocate(cxy(n))

    ! Compute the FFT's of both signals
    call zffti(n, work)
    call zfftf(n, x, work)
    call zfftf(n, y, work)

    ! Compute the scaling factor
    fac = 1.0d0 / real(n, real64)

    ! Perform the convolution operation
    if (method == 1) then
        ! Multiply the FFT's to convolve
        cxy = x * y
    else
        ! Divide the FFT's to deconvolve
        cxy = (0.0d0, 0.0d0)
        do i = 1, n
            if (abs(y(i)) > tiny(1.0d0)) cxy(i) = x(i) / y(i)
        end do
    end if

    ! Compute the inverse transform
    call zfftb(n, cxy, work)
    rst = real(fac * cxy)
end subroutine

! ------------------------------------------------------------------------------
end module