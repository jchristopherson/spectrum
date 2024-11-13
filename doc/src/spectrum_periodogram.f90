module spectrum_periodogram
    use iso_fortran_env
    use spectrum_windows
    use spectrum_routines
    use fftpack
    use ferror
    implicit none
    private
    public :: psd
    public :: csd
    public :: periodogram
    public :: cross_periodogram

contains
! ------------------------------------------------------------------------------
function psd(win, x, fs, nfft, err) result(rst)
    !! Computes the power spectral density (PSD) of a signal via Welch's
    !! method.
    !!
    !! References
    !!
    !! - Welch, P.D. "The Use of Fast Fourier Transform for the Estimation of 
    !!  Power Spectra: A Method Based on Time Averaging Over Short, Modified 
    !!  Periodograms." IEEE Transactions on Audio and Electroacoustics, 
    !!  AU-15 (2): 70-73, 1967.
    !!
    !! - [Wikipedia - Welch's Method](https://en.wikipedia.org/wiki/Welch%27s_method)
    class(window), intent(in) :: win
        !! The window to apply.  The size of the window must be non-zero and 
        !! positive-valued.
    real(real64), intent(in) :: x(:)
        !! The signal to transform.
    real(real64), intent(in), optional :: fs
        !! An optional input, that if supplied, allows for normalization of 
        !! the computed spectrum by the frequency resolution.
    integer(int32), intent(in), optional :: nfft
        !! An optional input that can be used to force the length of each
        !! individual DFT operation by padding any remaining space with zeros.
        !! If not supplied, the window size is used to determine the size of
        !! the DFT.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if win is not sized 
        !!      appropriately.
    real(real64), allocatable :: rst(:)
        !! An array containing the discrete PSD estimate.  The estimate is
        !! returned at discrete frequency intervals that can be determined by
        !! a call to frequency_bin_width.

    ! Local Variables
    logical :: init
    integer(int32) :: i, nx, nxfrm, nw, nk, nf, lwork, flag
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
    if (present(nfft)) then
        nf = nfft
    else
        nf = nw
    end if
    nxfrm = compute_transform_length(nf)
    nk = compute_overlap_segment_count(nx, nw)
    lwork = 3 * nf + 15

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
        call periodogram_driver(win, xw, buffer, fs, nf, work, init, cwork, &
            errmgr)
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

! ------------------------------------------------------------------------------
function csd(win, x, y, fs, nfft, err) result(rst)
    !! Computes the cross spectral density (CSD) of a signal via Welch's
    !! method (sometimes referred to as cross power spectral density
    !! or CPSD).
    !!
    !! References
    !!
    !! - Welch, P.D. "The Use of Fast Fourier Transform for the Estimation of 
    !!  Power Spectra: A Method Based on Time Averaging Over Short, Modified 
    !!  Periodograms." IEEE Transactions on Audio and Electroacoustics, 
    !!  AU-15 (2): 70-73, 1967.
    !!
    !! - [Wikipedia - Welch's Method](https://en.wikipedia.org/wiki/Welch%27s_method)
    !!
    !! - [Wikipedia - Cross Power Spectral Density](https://en.wikipedia.org/wiki/Spectral_density#Cross-spectral_density)
    class(window), intent(in) :: win
        !! The window to apply.  The size of the window must be non-zero and 
        !! positive-valued.
    real(real64), intent(in) :: x(:)
        !! The first N-element signal to transform.
    real(real64), intent(in) :: y(:)
        !! The second N-element signal to transform.
    real(real64), intent(in), optional :: fs
        !! An optional input, that if supplied, allows for normalization of 
        !! the computed spectrum by the frequency resolution.
    integer(int32), intent(in), optional :: nfft
        !! An optional input that can be used to force the length of each
        !! individual DFT operation by padding any remaining space with zeros.
        !! If not supplied, the window size is used to determine the size of
        !! the DFT.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if win is not sized 
        !!      appropriately.
        !!
        !! - SPCTRM_ARRAY_SIZE_MISMATCH_ERROR: Occurs if x and y are not the 
        !!      same size.
    complex(real64), allocatable :: rst(:)
        !! An array containing the discrete CSD estimate.  The estimate is
        !! returned at discrete frequency intervals that can be determined by
        !! a call to frequency_bin_width.

    ! Local Variables
    logical :: init
    integer(int32) :: i, nx, ny, nw, nk, nf, lwork, flag, nxfrm
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

! ------------------------------------------------------------------------------
module function periodogram(win, x, fs, nfft, err) result(rst)
    !! Computes the periodogram of a signal.
    class(window), intent(in) :: win
        !! The window to apply.  The size of the window must be non-zero and 
        !! positive-valued.
    real(real64), intent(in) :: x(:)
        !! The signal to transform.
    real(real64), intent(in), optional :: fs
        !! An optional input, that if supplied, allows for normalization of 
        !! the computed spectrum by the frequency resolution.
    integer(int32), intent(in), optional :: nfft
        !! An optional input that can be used to force the length of each
        !! individual DFT operation by padding any remaining space with zeros.
        !! If not supplied, the window size is used to determine the size of
        !! the DFT.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if win is not sized 
        !!      appropriately.
    real(real64), allocatable :: rst(:)
        !! An array containing the discrete PSD estimate.  The estimate is
        !! returned at discrete frequency intervals that can be determined by
        !! a call to frequency_bin_width.

    ! Local Variables
    integer(int32) :: n, nxfrm, nf, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    character(len = :), allocatable :: errmsg
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    if (present(nfft)) then
        nf = nfft
    else
        nf = n
    end if
    nxfrm = compute_transform_length(nf)

    ! Memory Allocation
    allocate(rst(nxfrm), stat = flag)
    if (flag /= 0) go to 10

    ! Process
    call periodogram_driver(win, x, rst, fs = fs, nfft = nfft, err = errmgr)
    
    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("periodogram", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------
module function cross_periodogram(win, x, y, fs, nfft, err) result(rst)
    !! Computes the cross-spectral periodogram of two signals.
    class(window), intent(in) :: win
        !! The window to apply.  The size of the window must be non-zero and 
        !! positive-valued.
    real(real64), intent(in) :: x(:)
        !! The first N-element signal to transform.
    real(real64), intent(in) :: y(:)
        !! The second N-element signal to transform.
    real(real64), intent(in), optional :: fs
        !! An optional input, that if supplied, allows for normalization of 
        !! the computed spectrum by the frequency resolution.
    integer(int32), intent(in), optional :: nfft
        !! An optional input that can be used to force the length of each
        !! individual DFT operation by padding any remaining space with zeros.
        !! If not supplied, the window size is used to determine the size of
        !! the DFT.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if win is not sized 
        !!      appropriately.
        !!
        !!  - SPCTRM_ARRAY_SIZE_MISMATCH_ERROR: Occurs if x and y are not the
        !!      same size, or if x and y are not the same size as the window
        !!      win.
    complex(real64), allocatable :: rst(:)
        !! An array containing the discrete cross periodogram estimate.  The 
        !! estimate is returned at discrete frequency intervals that can be 
        !! determined by a call to frequency_bin_width.

    ! Local Variables
    integer(int32) :: nx, nxfrm, nf, flag
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
    if (present(nfft)) then
        nf = nfft
    else
        nf = nx
    end if
    nxfrm = compute_transform_length(nf)

    ! Memory Allocation
    allocate(rst(nxfrm), stat = flag)
    if (flag /= 0) go to 10

    ! Process
    call cross_periodogram_driver(win, x, y, rst, fs = fs, nfft = nfft, &
        err = errmgr)

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("cross_periodogram", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------
subroutine periodogram_driver(win, x, xfrm, fs, nfft, work, initxfrm, &
    cwork, err)
    ! Arguments
    class(window), intent(in) :: win        ! size = n
    real(real64), intent(in) :: x(:)        ! size = n
    real(real64), intent(out) :: xfrm(:)    ! size = m = (nfft + 1) / 2 or nfft / 2 + 1
    real(real64), intent(in), optional :: fs
    integer(int32), intent(in), optional :: nfft ! defaults to n, must be at least n
    real(real64), intent(out), optional, target :: work(:)  ! size = 3 * nfft + 15
    logical, intent(in), optional :: initxfrm
    complex(real64), intent(out), optional, target :: cwork(:)  ! size = m
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    logical :: init
    integer(int32) :: i, j, nx, nxfrm, nf, lw, lwork, flag
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
    if (present(nfft)) then
        nf = nfft
    else
        nf = nx
    end if
    nxfrm = compute_transform_length(nf)
    lw = 2 * nf + 15
    lwork = lw + nf
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
    if (nf < nx) go to 60

    ! Workspace
    if (present(work)) then
        if (size(work) < lwork) go to 30
        w(1:lw) => work(1:lw)
        xw(1:nf) => work(lw+1:lwork)
    else
        allocate(wdef(lwork), stat = flag)
        if (flag /= 0) go to 10
        w(1:lw) => wdef(1:lw)
        xw(1:nf) => wdef(lw+1:lwork)
    end if
    if (init) call dffti(nf, w)

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

    ! Pad with zeros
    if (nf > nx) xw(nx+1:nf) = 0.0d0

    ! Compute the transform
    call dfftf(nf, xw, w)
    call unpack_real_transform(xw, cw, fac * scale)

    ! Compute the power from the windowed transform
    xfrm = real(cw * conjg(cw))

    ! Normalize by frequency
    if (present(fs)) then
        df = frequency_bin_width(fs, nf)
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

    ! Too small of NFFT
60  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The length of the FFT (", nf, &
        ") must be at least the size of the window (", nx, ")."
    call errmgr%report_error("periodogram_driver", trim(errmsg), &
        SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
subroutine cross_periodogram_driver(win, x, y, xfrm, fs, nfft, work, &
    initxfrm, cwork, err)
    ! Arguments
    class(window), intent(in) :: win        ! size = n
    real(real64), intent(in) :: x(:), y(:)  ! size = n
    complex(real64), intent(out) :: xfrm(:)    ! size = m = (nfft + 1) / 2 or nfft / 2 + 1
    real(real64), intent(in), optional :: fs
    integer(int32), intent(in), optional :: nfft ! defaults to n, must be at least n
    real(real64), intent(out), optional, target :: work(:)  ! size = 4 * nfft + 15
    logical, intent(in), optional :: initxfrm
    complex(real64), intent(out), optional, target :: cwork(:)  ! size = 2 * m
    class(errors), intent(inout), optional, target :: err

    ! Local Variables
    logical :: init
    integer(int32) :: i, j, nx, ny, nxfrm, nf, lw, lwork, flag
    real(real64) :: wval, wsum, scale, fac, df
    real(real64), allocatable, target, dimension(:) :: wdef
    real(real64), pointer, dimension(:) :: w, xw, yw
    complex(real64), allocatable, target, dimension(:) :: cwdef
    complex(real64), pointer, dimension(:) :: cx, cy
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
    if (present(nfft)) then
        nf = nfft
    else
        nf = nx
    end if
    nxfrm = compute_transform_length(nf)
    lw = 2 * nf + 15
    lwork = lw + 2 * nf
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
    if (nx /= ny) go to 60
    if (win%size /= nx) go to 20
    if (size(xfrm) /= nxfrm) go to 50
    if (nf < nx) go to 70

    ! Workspace
    if (present(work)) then
        if (size(work) < lwork) go to 30
        w(1:lw) => work(1:lw)
        xw(1:nf) => work(lw+1:lw+nf)
        yw(1:nf) => work(lw+nf+1:lwork)
    else
        allocate(wdef(lwork), stat = flag)
        if (flag /= 0) go to 10
        w(1:lw) => wdef(1:lw)
        xw(1:nf) => wdef(lw+1:lw+nf)
        yw(1:nf) => wdef(lw+nf+1:lwork)
    end if
    if (init) call dffti(nf, w)

    if (present(cwork)) then
        if (size(cwork) /= 2 * nxfrm) go to 40
        cx(1:nxfrm) => cwork(1:nxfrm)
        cy(1:nxfrm) => cwork(nxfrm+1:2*nxfrm)
    else
        allocate(cwdef(2 * nxfrm), stat = flag)
        if (flag /= 0) go to 10
        cx(1:nxfrm) => cwdef(1:nxfrm)
        cy(1:nxfrm) => cwdef(nxfrm+1:2*nxfrm)
    end if

    ! Apply the window
    wsum = 0.0d0
    do i = 1, nx
        j = i - 1
        wval = win%evaluate(j)
        wsum = wsum + wval
        xw(i) = wval * x(i)
        yw(i) = wval * y(i)
    end do
    fac = nx / wsum

    ! Pad with zeros
    if (nf > nx) then
        xw(nx+1:nf) = 0.0d0
        yw(ny+1:nf) = 0.0d0
    end if

    ! Compute the transforms
    call dfftf(nf, xw, w)
    call dfftf(nf, yw, w)
    call unpack_real_transform(xw, cx, fac * scale)
    call unpack_real_transform(yw, cy, fac * scale)

    ! Compute the cross spectrum
    xfrm = cx * conjg(cy)

    ! Normalize by frequency
    if (present(fs)) then
        df = frequency_bin_width(fs, nf)
        xfrm = xfrm / df
    end if

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Window Size Error
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The window size (", win%size, &
        ") must be the same size as the signal array (", nx, ")."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_MISMATCH_ERROR)
    return

    ! Workspace Too Small Error
30  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The workspace array was expected to have " // &
        "a length of ", lwork, " elements, but was found to have ", &
        size(work), " elements."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Complex-Valued Workspace Too Small Error
40  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The complex-valued workspace array was expected " // &
        "to have a length of ", 2 * nxfrm, &
        " elements, but was found to have ", size(cwork), " elements."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! Output Array Size Error
50  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The output array was expected to have a length of ", &
        nxfrm, " elements, but was found to have ", size(xfrm), " elements."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_ERROR)
    return

    ! X & Y are not the same size
60  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The first input array (size = ", nx, &
        "), and the second input array (size = ", ny, &
        ") must be the same size."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_ARRAY_SIZE_MISMATCH_ERROR)

! Too small of NFFT
70  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "The length of the FFT (", nf, &
        ") must be at least the size of the window (", nx, ")."
    call errmgr%report_error("cross_periodogram_driver", trim(errmsg), &
        SPCTRM_INVALID_INPUT_ERROR)
    return
    
    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end subroutine

! ------------------------------------------------------------------------------
end module