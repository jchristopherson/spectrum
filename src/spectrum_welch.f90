submodule (spectrum) spectrum_welch
    use fftpack
contains
! ------------------------------------------------------------------------------
module function welch_periodogram_1(win, x, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:)
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: nwin, nx, nxfrm, flag, lw
    real(real64), allocatable, dimension(:) :: work

    ! Initialization
    nwin = win%size
    nx = size(x)
    nxfrm = compute_transform_length(nwin)
    noverlaps = compute_overlap_segment_count(nx, nwin)
    lw = 2 * nwin + 15

    ! Input Checking

    ! Memory Allocation
    allocate(rst(nxfrm), stat = flag, source = (0.0d0, 0.0d0))
    if (flag == 0) allocate(work(lw), stat = flag)
    if (flag /= 0) go to 10

    ! Initialize the Fourier transform workspace array
    call dffti(nwin, work)



    ! End
    return

    ! Memory Error Handling
10  continue
    return
end function

! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
end submodule