submodule (spectrum) spectrum_convolve
    use fftpack
contains
! ------------------------------------------------------------------------------
module function convolve_1(x, y, method, err) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:), y(:)
    integer(int32), intent(in), optional :: method
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: mth, n1, n2
    real(real64), allocatable :: temp(:)
    complex(real64), allocatable :: xmod(:), ymod(:)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
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
    call prepare_conv(x, y, xmod, ymod, n1, n2, errmgr)
    if (errmgr%has_error_occurred()) return

    ! Compute the convolution
    if (mth == SPCTRM_CENTRAL_CONVOLUTION) then
        call convolve_driver(xmod, ymod, 1, temp, errmgr)
        if (errmgr%has_error_occurred()) return
        rst = temp(n1:n2)
    else
        call convolve_driver(xmod, ymod, 1, rst, errmgr)
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
subroutine prepare_conv(x, y, xmod, ymod, n1, n2, err)
    ! Arguments
    real(real64), intent(in) :: x(:), y(:)
    complex(real64), intent(out), allocatable :: xmod(:), ymod(:)
    integer(int32), intent(out) :: n1, n2
    class(errors), intent(inout) :: err

    ! Parameters
    complex(real64), parameter :: zero = (0.0d0, 0.0d0)

    ! Local Variables
    integer(int32) :: nx, ny, n, nbig, dn, m, flag
    logical :: evenY
    character(len = :), allocatable :: errmsg

    ! Initialization
    nx = size(x)
    ny = size(y)
    n = nx + ny - 1

    ! Memory Allocation
    allocate(xmod(n), stat = flag, source = zero)
    if (flag == 0) allocate(ymod(n), stat = flag, source = zero)
    if (flag /= 0) then
        ! ERROR:
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Memory allocation error flag ", flag, "."
        call err%report_error("prepare_conv", trim(errmsg), &
            SPCTRM_MEMORY_ERROR)
    end if

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

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
! Computes the convolution or deconvolution of two signals.
!
! - x: The zero-padded original signal.
! - y: The padded and wrapped response signal.
! - method: +1 for convolution, -1 for deconvolution.
! - rst: The convolved signal.
! - err: An error handling object.
subroutine convolve_driver(x, y, method, rst, err)
    ! Arguments
    complex(real64), intent(inout) :: x(:), y(:)
    integer(int32), intent(in) :: method
    real(real64), intent(out), allocatable :: rst(:)
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: i, n, lwork, flag, nend
    real(real64) :: fac
    real(real64), allocatable, dimension(:) :: work
    complex(real64), allocatable, dimension(:) :: cxy
    character(len = :), allocatable :: errmsg

    ! Generate the workspace array for the transforms & allocate memory
    n = size(x)
    lwork = 4 * n + 15
    allocate(work(lwork), stat = flag)
    if (flag == 0) allocate(rst(n), stat = flag)
    if (flag /= 0) then
        ! ERROR:
        allocate(character(len = 256) :: errmsg)
        write(errmsg, 100) "Memory allocation error flag ", flag, "."
        call err%report_error("convolve_driver", trim(errmsg), &
            SPCTRM_MEMORY_ERROR)
    end if

    ! Compute the FFT's of both signals
    call zffti(n, work)
    call zfftf(n, x, work)
    call zfftf(n, y, work)

    ! Compute the scaling factor
    fac = 2.0d0 / real(n, real64)

    ! Perform the convolution operation
    if (method == 1) then
        ! Multiply the FFT's to convolve
        cxy = fac**2 * x * y
    else
        ! Divide the FFT's to deconvolve
        cxy = x / y
    end if

    ! Compute the inverse transform
    call zfftb(n, cxy, work)
    rst = real(cxy)

    ! Formatting
100 format(A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end submodule