submodule (spectrum) spectrum_psd
    use fftpack
contains
module function psd_welch(win, x, fs, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:)
    real(real64), intent(in), optional :: fs
    class(errors), intent(inout), optional, target :: err
    real(real64), allocatable :: rst(:)

    ! Local Variables
    logical :: init
    integer(int32) :: i, nx, nxfrm, nw, nk, lwork, flag
    real(real64) :: fres, fac
    real(real64), allocatable, dimension(:) :: work, xw, buffer
    complex(real64), allocatable :: cwork(:)
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
    nw = win%size
    nxfrm = compute_transform_length(nw)
    nk = compute_overlap_segment_count(nx, nw)
    lwork = 3 * nw + 15

    ! Input Checking
    if (size(x) < 2) go to 20

    ! Memory Allocation
    allocate(rst(nxfrm), stat = flag, source = 0.0d0)
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (flag == 0) allocate(xw(nw), stat = flag)
    if (flag == 0) allocate(buffer(nxfrm), stat = flag)
    if (flag == 0) allocate(cwork(nxfrm), stat = flag)
    if (flag /= 0) go to 10
    
    ! Cycle over each segment
    init = .true.
    do i = 1, nk
        call overlap(x, i, nw, xw)
        call periodogram_driver(win, xw, buffer, fs, work, init, cwork, errmgr)
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
    call errmgr%report_error("psd_welch", trim(errmsg), SPCTRM_MEMORY_ERROR)
    return

    ! Window Size Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) &
        "The window must have at least 2 points, but was found to have ", &
        nx, "."
    call errmgr%report_error("psd_welch", trim(errmsg), &
        SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

end submodule