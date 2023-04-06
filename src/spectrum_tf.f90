submodule (spectrum) spectrum_tf
    use fftpack
contains
! ------------------------------------------------------------------------------
! REF: 
! - https://github.com/giuliovv/tfest/blob/main/tfest/tfest.py
! - https://dsp.stackexchange.com/questions/71811/understanding-the-h1-and-h2-estimators
module function siso_xfrm_h1(win, x, y, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:), y(:)
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable :: rst

    ! Local Variables
    integer(int32) :: m, n, nx, lt, flag
    real(real64), allocatable, dimension(:) :: trig, work
    complex(real64), allocatable, dimension(:) :: cwork
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
    n = win%size
    m = compute_transform_length(n)
    lw = 2 * n + 15

    ! Input Checking
    if (size(y) /= nx) go to 20
    if (n < 1) go to 30

    ! Memory Allocations
    allocate(rst(m), stat = flag, source = (0.0d0, 0.0d0))
    if (flag == 0) allocate(cwork(m), stat = flag)
    if (flag == 0) allocate(trig(lt), stat = flag)
    if (flag == 0) allocate(work(4 * n), stat = flag)
    if (flag /= 0) go to 10

    ! Initialize the transforms
    call dffti(n, trig)

    ! Accumulate and overlap to compute the transforms

    ! End
    return

    ! Memory Error Handling
10  continue
    return

    ! X & Y are not the same size
20  continue
    return

    ! Invalid window size
30  continue
    return

    ! Formatting
100 format(A, I0, A)
end function

! ------------------------------------------------------------------------------

! ******************************************************************************
! HELPER ROUTINES
! ------------------------------------------------------------------------------
! Segments, windows, and overlaps the supplied signals by 50% to accumulate an
! averaged transform.
!
! Inputs:
! - win: The window object.
! - x: An N-element array containing the input signal.
! - y: An N-element array containing the output signal.
! - trig: A 2*N+15 element array containing the trig terms for the FFT's.
! 
! Outputs:
! - xfrms: An M-element array containing the averaged transform.
!       M = N / 2 + 1, when N is even; else M = (N + 1) / 2.
! - work: A 4*N-element workspace array.
! - cwork: An M-element complex-valued workspace array.
subroutine overlap_segments(win, x, y, trig, xfrms, work, cwork)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:), y(:), trig(:)
    real(real64), intent(out), target :: work(:)
    complex(real64), intent(out) :: xfrms(:), cwork(:)

    ! Local Variables
    integer(int32) :: k, n, nk, lt
    real(real64), pointer, dimension(:) :: zx, zy, w

    ! Initialization
    xfrms = (0.0d0, 0.0d0)
    n = size(x)
    nk = compute_overlap_segment_count(n, win%size)

    ! Assign Pointers
    zx(1:n) => work(1:n)
    zy(1:n) => work(n+1:2*n)
    w(1:2*n) => work(2*n+1:4*n)

    ! Process
    do k = 1, nk
        call overlap(x, k, n, zx)
        call overlap(y, k, n, zy)
        call buffer_h1_segment(win, zx, zy, trig, xfrms, w, cwork)
    end do

    ! Average the transforms
    if (nk > 1) xfrms = xfrms / nk
end subroutine

! ------------------------------------------------------------------------------
! Adds a new set of transforms to the complex-valued buffer.
!
! Inputs:
! - win: The window object.
! - xnew: An N-element array containing the input signal.
! - ynew: An N-element array containing the output signal.
! - trig: A 2*N+15 element array containing the trig terms used for the FFT's.
! 
! Outputs:
! - buffer: An M-element array containing the complex-valued transforms.
!       M = N / 2 + 1, when N is even; else M = (N + 1) / 2.
! - work: A 2*N-element workspace array.
! - cwork: An M-element complex-valued workspace array.
subroutine buffer_h1_segment(win, xnew, ynew, trig, buffer, work, cwork)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: xnew(:), ynew(:), trig(:)
    complex(real64), intent(inout) :: buffer(:)
    real(real64), intent(out), target :: work(:)
    complex(real64), intent(out) :: cwork(:)

    ! Local Variables
    integer(int32) :: i, j, n, m
    real(real64) :: w, sumw, fac
    real(real64), pointer, dimension(:) :: xwin, ywin

    ! Initialization
    n = win%size
    m = compute_transform_length(n)

    ! Get pointers to the workspace arrays
    xwin(1:n) => work(1:n)
    ywin(1:n) => work(n + 1:2 * n)

    ! Apply the window to both signals
    sumw = 0.0d0
    do i = 1, n
        j = i - 1
        w = win%evaluate(j)
        sumw = sumw + w
        xwin(i) = w * xnew(i)
        ywin(i) = w * ynew(i)
    end do
    fac = n / sumw

    ! Compute the transforms
    call h1_transform(xwin, ywin, trig, cwork)
    buffer = buffer + cwork
end subroutine

! ------------------------------------------------------------------------------
! Computes the cross-spectrum and power spectrum transforms necessary for an H1
! transfer function estimate.
!
! Inputs:
! - x: An N-element array describing the system input.
! - y: An N-element array describing the system output.
! - trig: A 2*N+15 element workspace array containing the trig terms necessary
!       for the real-valued FFT algorithm.
!
! Outputs:
! - x: An N-element array containing the FFT of X.
! - y: An N-element array containing the FFT of Y.
! - xfrm: An M-element array where the transform will be written.  M = N / 2 + 1
!       if N is even; else, M = (N + 1) / 2.
subroutine h1_transform(x, y, trig, xfrm)
    ! Arguments
    real(real64), intent(inout) :: x(:), y(:)
    real(real64), intent(in) :: trig(:)
    complex(real64), intent(out) :: xfrm(:)

    ! Local Variables
    integer(int32) :: i, m, n, nend
    complex(real64) :: cx, cy

    ! Initialization
    n = size(x)
    if (mod(n, 2) == 0) then
        m = n / 2 + 1
        nend = m - 1
    else
        m = (n + 1) / 2
        nend = m
    end if

    ! Compute the transforms
    call dfftf(n, x, trig)
    call dfftf(n, y, trig)

    ! Compute the complex-valued transfer function
    xfrm(1) = y(1) / x(1)
    do i = 2, nend
        cx = cmplx(x(2*i-2), x(2*i-1), real64)
        cy = cmplx(y(2*i-2), y(2*i-1), real64)
        xfrm(i) = cy * conjg(cx) / (cx * conjg(cx))
    end do
    if (nend /= m) then
        xfrm(nend) = y(n) / x(n)
    end if
end subroutine

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end submodule