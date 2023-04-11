submodule (spectrum) spectrum_csd
    use fftpack
contains
module function csd_welch(win, x, y, fs, nfft, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:), y(:)
    real(real64), intent(in), optional :: fs
    integer(int32), intent(in), optional :: nfft
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable :: rst(:)

    ! Local Variables
    logical :: init
    integer(int32) :: nx, ny, nw, nk, nf, lwork, flag
    real(real64) :: fres, fac
    real(real64), allocatable, dimension(:) :: work, xw, yw
    complex(real64), allocatable, dimension(:) :: cwork, buffer
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
    ny = size(y)
    nw = win%size
    if (present(nfft)) then
        nf = nfft
    else
        nf = nw
    end if
    nxfrm = compute_transform_length(nf)
    nk = compute_overlap_segment_count(nx, nw)
    lwork = 4 * nf + 15

    ! Input Checking
    if (nx < 2 .or. ny < 2) go to 20

    ! Memory Allocation
    allocate(rst(nxfrm), stat = flag, source = (0.0d0, 0.0d0))
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (flag == 0) allocate(xw(nw), stat = flag)
    if (flag == 0) allocate(yw(nw), stat = flag)
    if (flag == 0) allocate(buffer(nxfrm), stat = flag)
    if (flag == 0) allocate(cwork(2 * nxfrm), stat = flag)
    if (flag /= 0) go to 10
    
    ! Cycle over each segment
    init = .true.
    do i = 1, nk
        call overlap(x, i, nw, xw)
        call overlap(y, i, nw, yw)
        call cross_periodogram_driver(win, xw, yw, buffer, fs, nf, work, &
            init, cwork, errmgr)
        if (errmgr%has_error_occurred()) return
        rst = rst + buffer
        init = .false.
    end do

    ! Average the result
    rst = rst / real(nk, real64)

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

end submodule