module spectrum_resample
    use iso_fortran_env
    use spectrum_filter
    implicit none
    private
    public :: upsample
    public :: downsample

contains
! ------------------------------------------------------------------------------
pure function upsample(n, fs, x) result(rst)
    !! Upsamples an evenly sampled signal by the specified factor.
    integer(int32), intent(in) :: n
        !! The upsample factor.  This value must be non-zero and positive 
        !! valued.
    real(real64), intent(in) :: fs
        !! The original signal sample rate, in Hz.
    real(real64), intent(in), dimension(:) :: x
        !! The signal to upsample.
    real(real64), allocatable, dimension(:) :: rst
        !! The upsampled signal.

    ! Local Variables
    integer(int32) :: i, npts, nnew
    
    ! Initialization
    npts = size(x)

    ! Input Checking
    if (n < 1 .or. fs <= 0.0d0 .or. npts < 1) return
    nnew = n * npts

    ! Memory Allocations
    allocate(rst(nnew), source = 0.0d0)

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
    rst = sinc_filter(0.5d0 * fs, fs * n, rst) * n
end function

! ------------------------------------------------------------------------------
pure function downsample(n, fs, x) result(rst)
    !! Downsamples an evenly sampled signal by the specified factor.
    integer(int32), intent(in) :: n
        !! The downsample factor.  The value must be non-zer and positive
        !! valued.
    real(real64), intent(in) :: fs
        !! The original signal sample rate, in Hz.
    real(real64), intent(in), dimension(:) :: x
        !! The signal to downsample.
    real(real64), allocatable, dimension(:) :: rst
        !! The downsampled signal.

    ! Local Variables
    integer(int32) :: i, npts, nnew
    real(real64), allocatable, dimension(:) :: xf
    
    ! Initialization
    npts = size(x)

    ! Input Checking
    if (n < 1 .or. fs <= 0.0d0 .or. npts < 1) return
    nnew = npts / n
    if (nnew < 1) return

    ! Memory Allocations
    allocate(rst(nnew), source = 0.0d0)

    ! Quick Return
    if (n == 1) then
        rst = x
        return
    end if

    ! Filter the signal at the downsampled frequency - prevents aliasing issues
    xf = sinc_filter(0.5d0 * fs / n, fs, x)

    ! Resample
    do i = 1, nnew
        rst(i) = xf(n * (i - 1) + 1)
    end do
end function

! ------------------------------------------------------------------------------
end module