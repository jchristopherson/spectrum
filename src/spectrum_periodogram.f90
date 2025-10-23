module spectrum_periodogram
    use iso_fortran_env
    use spectrum_windows
    use spectrum_routines
    use fftpack
    implicit none
    private
    public :: psd
    public :: csd
    public :: periodogram
    public :: cross_periodogram

contains
! ------------------------------------------------------------------------------
pure function psd(win, x, fs, nfft) result(rst)
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
    real(real64), allocatable :: rst(:)
        !! An array containing the discrete PSD estimate.  The estimate is
        !! returned at discrete frequency intervals that can be determined by
        !! a call to frequency_bin_width.

    ! Local Variables
    logical :: init
    integer(int32) :: i, nx, nxfrm, nw, nk, nf, lwork
    real(real64) :: fres, fac
    real(real64), allocatable, dimension(:) :: work, xw, buffer
    complex(real64), allocatable :: cwork(:)
    
    ! Initialization
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
    if (size(x) < 2) return

    ! Memory Allocation
    allocate(rst(nxfrm), source = 0.0d0)
    allocate(work(lwork))
    allocate(xw(nw))
    allocate(buffer(nxfrm))
    allocate(cwork(nxfrm))
    
    
    ! Cycle over each segment
    init = .true.
    do i = 1, nk
        call overlap(x, i, nw, xw)
        call periodogram_driver(win, xw, buffer, fs, nf, work, init, cwork)
        rst = rst + buffer
        init = .false.
    end do
    
    ! Average the result
    rst = rst / real(nk, real64)
end function

! ------------------------------------------------------------------------------
pure function csd(win, x, y, fs, nfft) result(rst)
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
    complex(real64), allocatable :: rst(:)
        !! An array containing the discrete CSD estimate.  The estimate is
        !! returned at discrete frequency intervals that can be determined by
        !! a call to frequency_bin_width.

    ! Local Variables
    logical :: init
    integer(int32) :: i, nx, ny, nw, nk, nf, lwork, nxfrm
    real(real64) :: fres, fac
    real(real64), allocatable, dimension(:) :: work, xw, yw
    complex(real64), allocatable, dimension(:) :: cwork, buffer
    
    ! Initialization
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
    if (nx < 2 .or. ny < 2) return

    ! Memory Allocation
    allocate(rst(nxfrm), source = (0.0d0, 0.0d0))
    allocate(work(lwork))
    allocate(xw(nw))
    allocate(yw(nw))
    allocate(buffer(nxfrm))
    allocate(cwork(2 * nxfrm))
    
    ! Cycle over each segment
    init = .true.
    do i = 1, nk
        call overlap(x, i, nw, xw)
        call overlap(y, i, nw, yw)
        call cross_periodogram_driver(win, xw, yw, buffer, fs, nf, work, &
            init, cwork)
        rst = rst + buffer
        init = .false.
    end do

    ! Average the result
    rst = rst / real(nk, real64)
end function

! ------------------------------------------------------------------------------
pure function periodogram(win, x, fs, nfft) result(rst)
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
    real(real64), allocatable :: rst(:)
        !! An array containing the discrete PSD estimate.  The estimate is
        !! returned at discrete frequency intervals that can be determined by
        !! a call to frequency_bin_width.

    ! Local Variables
    integer(int32) :: n, nxfrm, nf
    
    ! Initialization
    n = size(x)
    if (present(nfft)) then
        nf = nfft
    else
        nf = n
    end if
    nxfrm = compute_transform_length(nf)

    ! Memory Allocation
    allocate(rst(nxfrm))

    ! Process
    call periodogram_driver(win, x, rst, fs = fs, nfft = nfft)
end function

! ------------------------------------------------------------------------------
pure function cross_periodogram(win, x, y, fs, nfft) result(rst)
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
    complex(real64), allocatable :: rst(:)
        !! An array containing the discrete cross periodogram estimate.  The 
        !! estimate is returned at discrete frequency intervals that can be 
        !! determined by a call to frequency_bin_width.

    ! Local Variables
    integer(int32) :: nx, nxfrm, nf
    
    ! Initialization
    nx = size(x)
    if (present(nfft)) then
        nf = nfft
    else
        nf = nx
    end if
    nxfrm = compute_transform_length(nf)

    ! Memory Allocation
    allocate(rst(nxfrm))

    ! Process
    call cross_periodogram_driver(win, x, y, rst, fs = fs, nfft = nfft)
end function

! ------------------------------------------------------------------------------
pure subroutine periodogram_driver(win, x, xfrm, fs, nfft, work, initxfrm, &
    cwork)
    ! Arguments
    class(window), intent(in) :: win        ! size = n
    real(real64), intent(in) :: x(:)        ! size = n
    real(real64), intent(out) :: xfrm(:)    ! size = m = (nfft + 1) / 2 or nfft / 2 + 1
    real(real64), intent(in), optional :: fs
    integer(int32), intent(in), optional :: nfft ! defaults to n, must be at least n
    real(real64), intent(out), optional, target :: work(:)  ! size = 3 * nfft + 15
    logical, intent(in), optional :: initxfrm
    complex(real64), intent(out), optional, target :: cwork(:)  ! size = m

    ! Local Variables
    logical :: init
    integer(int32) :: i, j, nx, nxfrm, nf, lw, lwork
    real(real64) :: wval, wsum, scale, fac, df
    real(real64), allocatable, target, dimension(:) :: wdef
    real(real64), pointer, dimension(:) :: w, xw
    complex(real64), allocatable, target, dimension(:) :: cwdef
    complex(real64), pointer, dimension(:) :: cw
    
    ! Initialization
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
    if (win%size /= nx) return
    if (size(xfrm) /= nxfrm) return
    if (nf < nx) return

    ! Workspace
    if (present(work)) then
        if (size(work) < lwork) return
        w(1:lw) => work(1:lw)
        xw(1:nf) => work(lw+1:lwork)
    else
        allocate(wdef(lwork))
        w(1:lw) => wdef(1:lw)
        xw(1:nf) => wdef(lw+1:lwork)
    end if
    if (init) call dffti(nf, w)

    if (present(cwork)) then
        if (size(cwork) < nxfrm) return
        cw(1:nxfrm) => cwork(1:nxfrm)
    else
        allocate(cwdef(nxfrm))
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
    cw = unpack_real_transform(xw, fac * scale)

    ! Compute the power from the windowed transform
    xfrm = real(cw * conjg(cw))

    ! Normalize by frequency
    if (present(fs)) then
        df = frequency_bin_width(fs, nf)
        xfrm = xfrm / df
    end if
end subroutine

! ------------------------------------------------------------------------------
pure subroutine cross_periodogram_driver(win, x, y, xfrm, fs, nfft, work, &
    initxfrm, cwork)
    ! Arguments
    class(window), intent(in) :: win        ! size = n
    real(real64), intent(in) :: x(:), y(:)  ! size = n
    complex(real64), intent(out) :: xfrm(:)    ! size = m = (nfft + 1) / 2 or nfft / 2 + 1
    real(real64), intent(in), optional :: fs
    integer(int32), intent(in), optional :: nfft ! defaults to n, must be at least n
    real(real64), intent(out), optional, target :: work(:)  ! size = 4 * nfft + 15
    logical, intent(in), optional :: initxfrm
    complex(real64), intent(out), optional, target :: cwork(:)  ! size = 2 * m

    ! Local Variables
    logical :: init
    integer(int32) :: i, j, nx, ny, nxfrm, nf, lw, lwork
    real(real64) :: wval, wsum, scale, fac, df
    real(real64), allocatable, target, dimension(:) :: wdef
    real(real64), pointer, dimension(:) :: w, xw, yw
    complex(real64), allocatable, target, dimension(:) :: cwdef
    complex(real64), pointer, dimension(:) :: cx, cy

    ! Initialization
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
    if (nx /= ny) return
    if (win%size /= nx) return
    if (size(xfrm) /= nxfrm) return
    if (nf < nx) return

    ! Workspace
    if (present(work)) then
        if (size(work) < lwork) return
        w(1:lw) => work(1:lw)
        xw(1:nf) => work(lw+1:lw+nf)
        yw(1:nf) => work(lw+nf+1:lwork)
    else
        allocate(wdef(lwork))
        w(1:lw) => wdef(1:lw)
        xw(1:nf) => wdef(lw+1:lw+nf)
        yw(1:nf) => wdef(lw+nf+1:lwork)
    end if
    if (init) call dffti(nf, w)

    if (present(cwork)) then
        if (size(cwork) /= 2 * nxfrm) return
        cx(1:nxfrm) => cwork(1:nxfrm)
        cy(1:nxfrm) => cwork(nxfrm+1:2*nxfrm)
    else
        allocate(cwdef(2 * nxfrm))
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
    cx = unpack_real_transform(xw, fac * scale)
    cy = unpack_real_transform(yw, fac * scale)

    ! Compute the cross spectrum
    xfrm = cx * conjg(cy)

    ! Normalize by frequency
    if (present(fs)) then
        df = frequency_bin_width(fs, nf)
        xfrm = xfrm / df
    end if
end subroutine

! ------------------------------------------------------------------------------
end module