module spectrum_resample
    use iso_fortran_env
    use ferror
    use spectrum_filter
    use spectrum_errors
    implicit none
    private
    public :: upsample
    public :: downsample

contains
! ------------------------------------------------------------------------------
function upsample(n, fs, x, err) result(rst)
    !! Upsamples an evenly sampled signal by the specified factor.
    integer(int32), intent(in) :: n
        !! The upsample factor.  This value must be non-zero and positive 
        !! valued.
    real(real64), intent(in) :: fs
        !! The original signal sample rate, in Hz.
    real(real64), intent(in), dimension(:) :: x
        !! The signal to upsample.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !! - SPCTRM_INVALID_INPUT_ERROR: Occurs if n is negative or zero-valued.
    real(real64), allocatable, dimension(:) :: rst
        !! The upsampled signal.

    ! Local Variables
    integer(int32) :: i, npts, nnew, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    nnew = n * npts

    ! Input Checking
    if (n < 1) then
        call errmgr%report_error("upsample", &
            "The upsampling factor must be non-zero and positive-valued.", &
            SPCTRM_INVALID_INPUT_ERROR)
        return
    end if

    ! Memory Allocations
    allocate(rst(nnew), stat = flag, source = 0.0d0)
    if (flag /= 0) then
        call errmgr%report_error("upsample", "Memory allocation error.", &
            SPCTRM_MEMORY_ERROR)
        return
    end if

    ! Quick Return
    if (n == 1) then
        rst = x
        return
    end if

    ! Populate the "upsampled" signal
    do i = 1, npts
        rst(n * (i - 1) + 1) = x(i)
    end do

    ! Filter the upsampled frequency at the sample rate of the old signal
    rst = sinc_filter(0.5d0 * fs, fs * n, rst, err = errmgr) * n
end function

! ------------------------------------------------------------------------------
function downsample(n, fs, x, err) result(rst)
    !! Downsamples an evenly sampled signal by the specified factor.
    integer(int32), intent(in) :: n
        !! The downsample factor.  The value must be non-zer and positive
        !! valued.
    real(real64), intent(in) :: fs
        !! The original signal sample rate, in Hz.
    real(real64), intent(in), dimension(:) :: x
        !! The signal to downsample.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !! - SPCTRM_INVALID_INPUT_ERROR: Occurs if n is negative or zero-valued.
    real(real64), allocatable, dimension(:) :: rst
        !! The downsampled signal.

    ! Local Variables
    integer(int32) :: i, npts, nnew, flag
    real(real64), allocatable, dimension(:) :: xf
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    npts = size(x)
    nnew = npts / n

    ! Input Checking
    if (n < 1) then
        call errmgr%report_error("downsample", &
            "The downsampling factor must be non-zero and positive-valued.", &
            SPCTRM_INVALID_INPUT_ERROR)
        return
    end if

    ! Memory Allocations
    allocate(rst(nnew), stat = flag, source = 0.0d0)
    if (flag /= 0) then
        call errmgr%report_error("downsample", "Memory allocation error.", &
            SPCTRM_MEMORY_ERROR)
        return
    end if

    ! Quick Return
    if (n == 1) then
        rst = x
        return
    end if

    ! Filter the signal at the downsampled frequency - prevents aliasing issues
    xf = sinc_filter(0.5d0 * fs / n, fs, x, err = errmgr)
    if (errmgr%has_error_occurred()) return

    ! Resample
    do i = 1, nnew
        rst(i) = xf(n * (i - 1) + 1)
    end do
end function

! ------------------------------------------------------------------------------
end module