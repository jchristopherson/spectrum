submodule (spectrum) spectrum_welch
    use fftpack
contains
! ------------------------------------------------------------------------------
module subroutine periodogram_driver(win, x, xfrm, fs, work, initxfrm, cwork, err)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:)
    real(real64), intent(out) :: xfrm(:)
    real(real64), intent(in), optional :: fs
    real(real64), intent(out), optional, target :: work(:)
    logical, intent(in), optional :: initxfrm
    complex(real64), intent(out), optional, target :: cwork(:)
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    logical :: init
    integer(int32) :: i, j, nx, nxfrm, lw, lwork, flag
    real(real64) :: wval, wsum, scale, fac, df
    real(real64), allocatable, target, dimension(:) :: wdef
    real(real64), pointer, dimension(:) :: w, xw
    complex(real64), allocatable, target, dimension(:) :: cwdef
    complex(real64), pointer, dimension(:) :: cw
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
    nxfrm = compute_transform_length(nx)
    lw = 2 * nx + 15
    lwork = lw + nx
    if (present(work) .and. present(initxfrm)) then
        init = initxfrm
    else
        init = .true.
    end if
    if (mod(nx, 2) == 0) then
        scale = 2.0d0 / nx
    else
        scale = 2.0d0 / (nx - 1.0d0)
    end if

    ! Input Checking
    if (win%size /= nx) go to 20
    if (size(xfrm) /= nxfrm) go to 50

    ! Workspace
    if (present(work)) then
        if (size(work) < lwork) go to 30
        w(1:lw) => work(1:lw)
        xw(1:nx) => work(lw+1:lwork)
    else
        allocate(wdef(lwork), stat = flag)
        if (flag /= 0) go to 10
        w(1:lw) => wdef(1:lw)
        xw(1:nx) => wdef(lw+1:lwork)
    end if
    if (init) call dffti(nx, w)

    if (present(cwork)) then
        if (size(cwork) < nxfrm) go to 40
        cw(1:nxfrm) => cwork(1:nxfrm)
    else
        allocate(cwdef(nxfrm), stat = flag)
        if (flag /= 0) go to 10
        cw(1:nxfrm) => cwdef(1:nxfrm)
    end if

    ! Apply the window
    wsum = 0.0d0
    do i = 1, nx
        j = i - 1
        wval = win%evaluate(j)
        wsum = wsum + wval
        xw(i) = wval * x(i)
    end do
    fac = nx / wsum

    ! Compute the transform
    call dfftf(nx, xw, w)
    call unpack_real_transform(xw, cw, scale)

    ! Compute the power from the windowed transform
    xfrm = 2.0d0 * fac * real(cw * conjg(cw))

    ! Normalize by frequency
    if (present(fs)) then
        df = frequency_bin_width(fs, nx)
        xfrm = xfrm / df
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("periodogram_driver", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Window Size Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The window size (", win%size, &
        ") must be the same size as the signal array (", nx, ")."
    call errmgr%report_error("periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_MISMATCH_ERROR)
    return

    ! Workspace Too Small Error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The workspace array was expected to have " // &
        "a length of ", lwork, " elements, but was found to have ", &
        size(work), " elements."
    call errmgr%report_error("periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Complex-Valued Workspace Too Small Error
40  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The complex-valued workspace array was expected " // &
        "to have a length of ", nxfrm, " elements, but was found to have ", &
        size(cwork), " elements."
    call errmgr%report_error("periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Output Array Size Error
50  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The output array was expected to have a length of ", &
        nxfrm, " elements, but was found to have ", size(xfrm), " elements."
    call errmgr%report_error("periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
end submodule