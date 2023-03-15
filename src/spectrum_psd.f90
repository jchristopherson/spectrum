submodule (spectrum) spectrum_psd
    use fftpack
contains
! ------------------------------------------------------------------------------
module function psd_welch(win, x, fs, err) result(rst)
    ! Arguments
end function

! ------------------------------------------------------------------------------
! Compiles a series of overlapped transforms (sum of the square of the transform 
! magnitude - sum of the power for each overlapped segment).
!
! - win: The window object with a window size of N.
! - x: The P-element array containing the signal to transform.  If P is less 
!       than N, P is padded with zeros to be N-elements long; else, it is 
!       overlapped.
! - xfrms: An M-element array where the transformed data will be written.
! - noverlaps: The number of overlapped segments.
! - work: A 4 * N + 15 element workspace array.
subroutine overlap_segments(win, x, xfrms, noverlaps, work)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:)
    real(real64), intent(out) :: xfrms(:)
    integer(int32), intent(out) :: noverlaps
    real(real64), intent(out), target :: work(:)

    ! Local Variables
    integer(int32) :: i, k, nk, nx, m, n, noff, lw, lwork
    real(real64) :: del
    real(real64), pointer :: z(:), w(:)

    ! Initialization
    xfrms = 0.0d0
    nx = size(x)
    n = win%size
    m = compute_transform_length(n)
    nk = (nx - 1) / m
    if (nk > 1) then
        del = (nx - n) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if
    noverlaps = 0
    lw = 3 * n + 15
    lwork = lw + n

    ! Assign Pointers
    z(1:n) => work(1:n)
    w(1:lw) => work(n+1:lwork)

    ! Process
    if (nx < n) then
        ! Pad with zeros
        z(1:nx) = x
        z(nx+1:n) = 0.0d0

        ! Add to the buffer
        call add_segment(win, z, xfrms, noverlaps, w)
    else
        ! Overlap
        do k = 0, nk - 1
            noff = int(k * del + 0.5d0, int32)
            do i = 1, n
                z(i) = x(noff + i)
            end do
            call add_segment(win, z, xfrms, noverlaps, w)
        end do
    end if
end subroutine

! ------------------------------------------------------------------------------
! Adds another segment to a pre-existing buffer by computing it's windowed FFT,
! and then squaring it and adding the squared trasform to the buffer. 
!
! Arguments:
! - win: The window object
! - xnew: An N-element array containing the new signal to buffer
! - buffer: An M-element array containing the already buffered signals, and on
!       output, the updated buffer.
! - counter: A tracking variable determining how many signals have been 
!       buffered.
! - work: A 3 * N + 15 element workspace array.
subroutine buffer_segment(win, xnew, buffer, counter, work)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: xnew(:)
    real(real64), intent(inout) :: buffer(:)
    integer(int32), intent(inout) :: counter
    real(real64), intent(inout), target :: work(:)

    ! Local Variables
    integer(int32) :: i, j, m, n, flag, lwork, nend
    real(real64) :: sumw, w, fac, scale
    real(real64), pointer :: wxfrm(:), xfrm(:)
    complex(real64) :: cx

    ! Initialization
    n = win%size
    m = compute_transform_length(n)
    lwork = 2 * n + 15

    ! Get pointers to the workspace arrays
    wxfrm(1:lwork) => work(1:lwork)
    xfrm(1:n) => work(lwork+1:n+lwork)

    ! Initialize the FFT workspace if needed
    if (counter < 1) then
        call dffti(n, wxfrm)
    end if

    ! Apply the window
    sumw = 0.0d0
    do i = 1, n
        j = i - 1
        w = win%evaluate(j)
        sumw = sumw + w
        xfrm(i) = w * xnew(i)
    end do
    fac = n / sumw

    ! Compute the FFT of the windowed signal
    call dfftf(n, xfrm, wxfrm)

    ! Buffer the result
    if (mod(n, 2) == 0) then
        scale = 2.0d0 / n   ! scale the transform by it's length
        nend = m - 1
    else
        scale = 2.0d0 / (n - 1.0d0) ! scale the transform by it's length
        nend = m
    end if
    buffer(1) = fac * abs(scale * xfrm(1))**2
    do i = 2, nend
        cx = scale * cmplx(xfrm(2 * i - 2), xfrm(2 * i - 1), real64)
        buffer(i) = fac * abs(cx)**2
    end do
    if (nend /= m) then
        buffer(nend) = fac * abs(scale * xfrm(n))**2
    end if

    ! Increment the counter
    counter = counter + 1
end subroutine

! ------------------------------------------------------------------------------
end submodule