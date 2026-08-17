module spectrum_fft
    use iso_fortran_env
    use fftpack, only : dffti, dfftf, dfftb
    use spectrum_windows
    use spectrum_routines
    implicit none
    private
    public :: stft
    public :: rfft
    public :: irfft
    public :: stft_result

    type stft_result
        !! A container for STFT output.
        complex(real64), allocatable, dimension(:,:) :: stft
            !! An M-by-N matrix containing the M-element complex-valued 
            !! transforms for each of the N time points studied.  M is the 
            !! size of the positive half of the transform, and N is the total 
            !! number of transformed segments.
        integer(int32), allocatable, dimension(:) :: offsets
            !! The starting indices of each window segment.
    end type

contains
! ------------------------------------------------------------------------------
pure function rfft(x, n) result(rst)
    !! Computes the Fourier transform of a real-valued data set.  Only the 
    !! positive half of the transform is returned.
    real(real64), intent(in), dimension(:) :: x
        !! The real-valued array to transform.
    integer(int32), intent(in), optional :: n
        !! An optional input that can be used to specify the length of the 
        !! transform.  If less than the length of x, x is truncated; however, 
        !! if greater than the length of x, x is padded with zeros.
    complex(real64), allocatable, dimension(:) :: rst
        !! The complex-valued result of the transform.

    ! Local Variables
    integer(int32) :: i, m, nx, nxfrm, lensav
    real(real64), allocatable, dimension(:) :: wsave, rrst

    ! Initialization
    nx = size(x)
    if (present(n)) then
        nxfrm = n
        if (nxfrm <= nx) then
            rrst = x(1:nxfrm)
        else
            rrst = [x, (0.0d0, i = 1, nxfrm - nx)]
        end if
    else
        nxfrm = nx
        rrst = x
    end if
    lensav = 2 * nxfrm + 15
    allocate(wsave(lensav))
    call dffti(nxfrm, wsave)

    ! Process
    call dfftf(nxfrm, rrst, wsave)
    m = compute_transform_length(nxfrm)
    allocate(rst(m))
    rst = unpack_real_transform(rrst)
end function

! ------------------------------------------------------------------------------
pure function irfft(x, n) result(rst)
    !! Computes the inverse Fourier transform for a real-valued data set.
    complex(real64), intent(in), dimension(:) :: x
        !! The positive half of the transform.
    integer(int32), intent(in), optional :: n
        !! An optional input that can be used to specify the length of the 
        !! transform.  If less than the length of x, x is truncated; however, 
        !! if greater than the length of x, x is padded with zeros.
    real(real64), allocatable, dimension(:) :: rst
        !! The real-valued result of the inverse transform.

    ! Local Variables
    integer(int32) :: i, m, nx, nxfrm, lensav, mn
    real(real64), allocatable, dimension(:) :: wsave

    ! Initialization
    nx = size(x)
    if (mod(nx, 2) == 0) then
        m = 2 * nx - 1
    else
        m = 2 * (nx - 1)
    end if
    allocate(rst(m), source = 0.0d0)
    if (present(n)) then
        nxfrm = n
    else
        nxfrm = nx
    end if
    mn = min((m + 1) / 2, nx)
    rst(1) = real(x(1), real64)
    do i = 2, mn
        rst(2*i-2) = real(x(i), real64)
        rst(2*i-1) = aimag(x(i))
    end do
    lensav = 2 * nxfrm + 15
    allocate(wsave(lensav))

    ! Process
    call dfftb(nxfrm, rst, wsave)
end function

! ------------------------------------------------------------------------------
pure function stft(win, x) result(rst)
    !! Computes the short time Fourier transform of a signal.
    !!
    !! See Also
    !!
    !! - [Wikipedia - Short Time Fourier Transform](https://en.wikipedia.org/wiki/Short-time_Fourier_transform)
    class(window), intent(in) :: win
        !! The window to apply.
    real(real64), intent(in) :: x(:)
        !! The signal to analyze.  The signal must be longer than the size of 
        !! the window.
    type(stft_result) :: rst
        !! A container filled with the transform results.

    ! Local Variables
    integer(int32) :: i, j, k, m, nx, nxfrm, nk, lwork, i1, nend
    real(real64) :: del, sumw, w, fac, scale
    real(real64), allocatable, dimension(:) :: work, buffer
    
    ! Initialization
    nx = size(x)
    m = win%size          ! # of rows in the output matrix (transform length)
    if (m < 1 .or. nx < m) return
    nxfrm = compute_transform_length(m)
    lwork = 2 * m + 15
    nk = compute_overlap_segment_count(nx, m)
    if (nk > 1) then
        del = (nx - m) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if
    if (mod(m, 2) == 0) then
        nend = nxfrm - 1
        scale = 2.0d0 / m
    else
        nend = nxfrm
        scale = 2.0d0 / (m - 1.0d0)
    end if

    allocate(work(lwork), buffer(m))
    allocate(rst%stft(nxfrm, nk))
    allocate(rst%offsets(nk), source = 0)

    ! Initialize the the transform.  As all transforms are the 
    ! same length, this array can be shared.
    call dffti(m, work)

    ! Compute each transform
    do i = 1, nk
        ! Compute the offset
        i1 = int((i - 1) * del + 0.5d0, int32)
        ! if (present(offsets)) offsets(i) = i1 + 1
        rst%offsets(i) = i1 + 1

        ! Apply the window
        sumw = 0.0d0
        j = 0
        do k = 1, m
            w = win%evaluate(j)
            j = j + 1
            sumw = sumw + w
            buffer(k) = w * x(k + i1)
        end do
        if (abs(sumw) <= tiny(1.0d0)) return
        fac = m / sumw

        ! Compute the transform
        call dfftf(m, buffer, work)

        ! Scale the transform
        if (m == 1) then
            rst%stft(1,i) = fac * cmplx(buffer(1), 0.0d0, real64)
        else
            rst%stft(1,i) = 0.5d0 * fac * scale * &
                cmplx(buffer(1), 0.0d0, real64)
        end if
        do k = 2, nend
            rst%stft(k,i) = fac * scale * cmplx(buffer(2*k-1), buffer(2*k), real64)
        end do
        if (nend /= nxfrm) then
            rst%stft(nxfrm,i) = 0.5d0 * fac * scale * &
                cmplx(buffer(m), 0.0d0, real64)
        end if
    end do
end function

! ------------------------------------------------------------------------------
end module