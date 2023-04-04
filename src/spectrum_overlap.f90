submodule (spectrum) spectrum_overlap
contains

pure module function overlap_segment_count_1(n, winsize) result(rst)
    ! Arguments
    integer(int32), intent(in) :: n, winsize
    integer(int32) :: rst

    ! Local Variables
    integer(int32) :: nxfrm

    ! Process
    nxfrm = compute_transform_length(winsize)
    rst = (n - 1) / nxfrm
end function

module subroutine fill_overlap_buffer_1(x, seg, winsize, buffer, err)
    ! Arguments
    real(real64), intent(in) :: x(:)
    integer(int32), intent(in) :: seg, winsize
    real(real64), intent(out) :: buffer(:)
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    integer(int32) :: k, nx, nxfrm, nk, noff, i1, i2
    real(real64) :: del
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
    nxfrm = compute_transform_length(winsize)
    nk = (nx - 1) / m
    if (nk > 1) then
        del = (nx - winsize) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if
    
    ! Input Checking
    if (winsize < 1) go to 10

    ! Slice up the array
    if (winsize > nx) then
        buffer(1:nx) = x
        buffer(nx+1:winsize) = 0.0d0
    else
        k = seg - 1
        noff = int(k * del + 0.5d0, int32)
        i1 = noff + 1
        i2 = i1 + winsize - 1
        buffer = x(i1:i2)
    end if

    ! End
    return

    ! Window Size Error
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) &
        "The window size must be greater than 1, but was found to be ", &
        winsize, &
        "."
    call errmgr%report_error("fill_overlap_buffer_1", trim(errmsg), &
        SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end subroutine

end submodule