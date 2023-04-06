submodule (spectrum) spectrum_csd
    use fftpack
contains
! ------------------------------------------------------------------------------
module function csd_welch(win, x, y, fs, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:), y(:)
    real(real64), intent(in), optional :: fs
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: nx, ny, nwin, nxfrm, noverlaps, lwork, flag
    real(real64) :: fres, fac
    real(real64), allocatable :: work(:)
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    nwin = win%size
    nx = size(x)
    ny = size(y)
    nxfrm = compute_transform_length(nwin)
    lwork = 6 * nwin + 15
    fres = 1.0d0
    if (present(fs)) then
        fres = frequency_bin_width(fs, nwin)
    end if

    ! Input Check
    if (nwin < 2) go to 20
    if (nx /= ny) go to 30

    ! Memory Allocation
    allocate(rst(nxfrm), stat = flag)
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (flag /= 0) go to 10

    ! Overlap and compute the transforms
    call overlap_segments(win, x, y, rst, noverlaps, work)

    ! Normalize
    fac = 1.0d0 / (fres * noverlaps)
    rst = fac * rst

    ! Account for the symmetry of the transform
    rst(2:nxfrm - 1) = 2.0d0 * rst(2:nxfrm - 1)

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("csd_welch", trim(errmsg), SPCTRM_MEMORY_ERROR)
    return

    ! Window Size Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) &
        "The window must have at least 2 points, but was found to have ", &
        nx, "."
    call errmgr%report_error("csd_welch", trim(errmsg), &
        SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Array Size Mismatch Error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The two arrays must be the same length.  " // &
        "Array X was found to be of length ", nx, ", and array Y was " // &
        "found to be of length ", ny, "."
    call errmgr%report_error("csd_welch", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_MISMATCH_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end function

! ------------------------------------------------------------------------------
! Compiles a series of overlapped transforms.
!
! - win: The window object with a window size of N.
! - x: The P-element array containing the first signal to transform.  If P is
!       less than N, P is padded with zeros to be N-elements long; else, it is 
!       overlapped.
! - y: The P-element array containing the second signal to transform.  If P is
!       less than N, P is padded with zeros to be N-elements long; else, it is 
!       overlapped.
! - xfrms: An M-element array where the transformed data will be written.
! - noverlaps: The number of overlapped segments.
! - work: A 6 * N + 15 element workspace array.
subroutine overlap_segments(win, x, y, xfrms, noverlaps, work)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:), y(:)
    real(real64), intent(out) :: xfrms(:)
    integer(int32), intent(out) :: noverlaps
    real(real64), intent(out), target :: work(:)

    ! Local Variables
    integer(int32) :: k, nx, n, lw, lwork, nk
    real(real64), pointer :: zx(:), zy(:), w(:)

    ! Initialization
    xfrms = 0.0d0
    nx = size(x)
    n = win%size
    nk = compute_overlap_segment_count(nx, n)
    noverlaps = 0
    lw = 4 * n + 15
    lwork = lw + 2 * n

    ! Assign Pointers
    zx(1:n) => work(1:n)
    zy(1:n) => work(n + 1:2 * n)
    w(1:lw) => work(2 * n + 1:lwork)

    ! Process
    do k = 1, nk
        call overlap(x, k, n, zx)
        call overlap(y, k, n, zy)
        call buffer_segment(win, zx, zy, xfrms, noverlaps, w)
    end do
end subroutine

! ------------------------------------------------------------------------------
! Adds another segment to a pre-existing buffer by computing the windowed FFT
! of the two signals of equal length, and adding the product of their transforms
! to the buffer.
!
! Arguments:
! - win: The window object
! - xnew: An N-element array containing the first signal to analyze.
! - ynew: An N-element array containing the second signal to analyze.
! - buffer: An M-element array containing the already buffered signals, and on
!       output, the updated buffer.
! - counter: A tracking variable determining how many signals have been 
!       buffered.
! - work: A 4 * N + 15 element workspace array.
subroutine buffer_segment(win, xnew, ynew, buffer, counter, work)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: xnew(:), ynew(:)
    real(real64), intent(inout) :: buffer(:)
    integer(int32), intent(inout) :: counter
    real(real64), intent(inout), target :: work(:)

    ! Local Variables
    integer(int32) :: i, j, m, n, flag, lwork, nend
    real(real64) :: sumw, w, fac, scale
    real(real64), pointer :: wxfrm(:), xxfrm(:), yxfrm(:)
    complex(real64) :: cx, cy

    ! Initialization
    n = win%size
    m = compute_transform_length(n)
    lwork = 2 * n + 15

    ! Get pointers to the workspace arrays
    wxfrm(1:lwork) => work(1:lwork)
    xxfrm(1:n) => work(lwork + 1:n + lwork)
    yxfrm(1:n) => work(lwork + n + 1:lwork + 2 * n)

    ! Initialize the FFT workspace - it can be shared between both transforms
    if (counter < 1) then
        call dffti(n, wxfrm)
    end if

    ! Window both signals
    sumw = 0.0d0
    do i = 1, n
        j = i - 1
        w = win%evaluate(j)
        sumw = sumw + w
        xxfrm(i) = w * xnew(i)
        yxfrm(i) = w * ynew(i)
    end do
    fac = n / sumw

    ! Compute the FFT's of the windowed signals
    call dfftf(n, xxfrm, wxfrm)
    call dfftf(n, yxfrm, wxfrm)

    ! Buffer the result
    if (mod(n, 2) == 0) then
        scale = 2.0d0 / n   ! scale the transform by it's length
        nend = m - 1
    else
        scale = 2.0d0 / (n - 1.0d0) ! scale the transform by it's length
        nend = m
    end if

    buffer(1) = fac * abs((scale * xxfrm(1)) * (scale * yxfrm(1))) + buffer(1)
    do i = 2, nend
        cx = scale * cmplx(xxfrm(2 * i - 2), xxfrm(2 * i - 1), real64)
        cy = scale * cmplx(yxfrm(2 * i - 2), yxfrm(2 * i - 1), real64)
        buffer(i) = fac * abs(cx * cy) + buffer(i)
    end do
    if (nend /= m) then
        buffer(nend) = fac * abs((scale * xxfrm(n)) * (scale * yxfrm(n))) + &
            buffer(nend)
    end if

    ! Increment the counter
    counter = counter + 1
end subroutine

! ------------------------------------------------------------------------------
end submodule