module spectrum_fft
    use iso_fortran_env
    use fftpack
    use spectrum_windows
    use spectrum_routines
    use ferror
    implicit none
    private
    public :: stft

contains
! ------------------------------------------------------------------------------
function stft(win, x, offsets, err) result(rst)
    !! Computes the short time Fourier transform of a signal.
    !!
    !! See Also
    !!
    !! - [Wikipedia - Short Time Fourier Transform](https://en.wikipedia.org/wiki/Short-time_Fourier_transform)
    class(window), intent(in) :: win
        !! The window to apply.
    real(real64), intent(in) :: x(:)
        !! The signal to analyze.  The signal must be longer than the size of 
        !! the window.
    integer(int32), intent(out), optional, allocatable :: offsets(:)
        !! An optional allocatable array that, if supplied, will be filled with 
        !! the starting indices of each window segment.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if the signal in x is too short
        !!      relative to the window size in win.
    complex(real64), allocatable :: rst(:,:)
        !! An M-by-N matrix containing the M-element complex-valued 
        !! transforms for each of the N time points studied.  M is the size of 
        !! the positive half of the transform, and N is the total number of 
        !! transformed segments.

    ! Local Variables
    integer(int32) :: i, j, k, m, nx, nxfrm, nk, lwork, flag, i1, nend
    real(real64) :: del, sumw, w, fac, scale
    real(real64), allocatable, dimension(:) :: work, buffer
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nx = size(x)
    m = win%size          ! # of rows in the output matrix (transform length)
    nxfrm = compute_transform_length(m)
    lwork = 2 * m + 15
    nk = compute_overlap_segment_count(nx, m)
    if (nk > 1) then
        del = (nx - m) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if
    if (mod(m, 2) == 0) then
        nend = nxfrm - 1
        scale = 2.0d0 / m
    else
        nend = nxfrm
        scale = 2.0d0 / (m - 1.0d0)
    end if

    ! Input Checking
    if (m > nx) go to 20

    ! Memory Allocation
    allocate(rst(nxfrm, nk), stat = flag)
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (flag == 0) allocate(buffer(m), stat = flag)
    if (present(offsets)) then
        if (flag == 0) allocate(offsets(nk), stat = flag, source = 0)
    end if
    if (flag /= 0) go to 10

    ! Initialize the trig terms for the transforms.  As all transforms are the 
    ! same length, this array can be shared.
    call dffti(m, work)

    ! Compute each transform
    do i = 1, nk
        ! Compute the offset
        i1 = int((i - 1) * del + 0.5d0, int32)
        if (present(offsets)) offsets(i) = i1 + 1

        ! Apply the window
        sumw = 0.0d0
        j = 0
        do k = 1, m
            w = win%evaluate(j)
            j = j + 1
            sumw = sumw + w
            buffer(k) = w * x(k + i1)
        end do
        fac = m / sumw

        ! Compute the transform
        call dfftf(m, buffer, work)

        ! Scale the transform
        rst(1,i) = fac * scale * cmplx(buffer(1), 0.0d0, real64)
        do k = 2, nend
            rst(k,i) = fac * scale * cmplx(buffer(2*k-1), buffer(2*k), real64)
        end do
        if (nend /= nxfrm) then
            rst(nxfrm,i) = fac * scale * cmplx(buffer(m), 0.0d0, real64)
        end if
    end do

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("stft", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Window size exceeds signal length
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "Window size (", m, ") exceeds the signal length (", &
        nx, ")."
    call errmgr%report_error("stft", trim(errmsg), SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end function

! ------------------------------------------------------------------------------
end module