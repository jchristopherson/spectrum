module spectrum_routines
    use iso_fortran_env
    implicit none
    private
    public :: compute_transform_length
    public :: frequency_bin_width
    public :: is_power_of_two
    public :: next_power_of_two
    public :: unpack_real_transform
    public :: cumulative_sum
    public :: unwrap
    public :: difference
    public :: compute_overlap_segment_count
    public :: overlap

contains
! ------------------------------------------------------------------------------
pure function compute_transform_length(n) result(rst)
    !! Computes the length of the positive half of a discrete Fourier
    !! transform for a specific signal length.
    !!
    !! The length of the discrete Fourier transform is determined as follows.
    !! 
    !! If n is even,
    !! $$ m = \frac{n}{2} + 1 $$
    !!
    !! else;
    !! $$ m = \frac{n + 1}{2} $$.
    integer(int32), intent(in) :: n
        !! The signal length (length of the signal put forth to the Fourier 
        !! transform).
    integer(int32) :: rst
        !! The length of the positive half of the discrete Fourier transform.
    if (mod(n, 2) == 0) then
        rst = n / 2 + 1
    else
        rst = (n + 1) / 2
    end if
end function

! ------------------------------------------------------------------------------
pure elemental function frequency_bin_width(fs, n) result(rst)
    !! Computes the bin width for a discrete frequency spectrum based
    !! upon the data sample rate.
    !!
    !! The frequency bin width is computed as follows.
    !!
    !! $$ \Delta f = \frac{f_{s}}{2 \left(m - 1 \right)} $$
    !!
    !! where m is computed via compute_transform_length.
    real(real64), intent(in) :: fs
        !! The rate at which the signal was sampled.  The units of this value 
        !! will be the units of the output.
    integer(int32), intent(in) :: n
        !! The signal length (length of the signal put forth to the
        !! Fourier transform).
    real(real64) :: rst
        !! The frequency bin width.

    ! Local Variables
    integer(int32) :: m

    ! Process
    m = compute_transform_length(n)
    rst = 0.5d0 * fs / (m - 1.0d0)
end function

! ------------------------------------------------------------------------------
pure elemental function is_power_of_two(n) result(rst)
    !! Tests to see if a value is an integer power of two.
    integer(int32), intent(in) :: n
        !! The value to test.
    logical :: rst
        !! Returns true if n is a power of two; else, false.

    ! Process
    rst = (n /= 0) .and. (iand(n, n - 1) == 0)
    ! This is equivalent to the C form : n != 0 && ((n & (n - 1)) == 0)
end function

! ------------------------------------------------------------------------------
pure elemental function next_power_of_two(n) result(rst)
    !! Provides the next higher integer power of two.
    integer(int32), intent(in) :: n
        !! The value to test.
    integer(int32) :: rst
        !! The next power of two higher than n.  If n is already a power of two,
        !! its value is simply returned.  For instance, if n is set to 128, 
        !! then a value of 7 is returned ( \(2^7 = 128 \) ).

    ! Process
    rst = ceiling(log(real(n)) / log(2.0), int32)
end function

! ------------------------------------------------------------------------------
pure subroutine unpack_real_transform(x, cx, fac)
    !! Unpacks a real-valued transform into its complex-valued format.
    real(real64), intent(in) :: x(:)
        !! The complex-valued signal stored in a real-valued array.  This array
        !! is assumed to be of length N.
    complex(real64), intent(out) :: cx(:)
        !! The M-element array where the complex form of x is written.  M is
        !! determined by calling compute_transform_length(N).
    real(real64), intent(in), optional :: fac
        !! An optional scaling input.  The default is 1 such that no scaling
        !! is performed.

    ! Local Variables
    integer(int32) :: i, nx, nxfrm
    real(real64) :: f
    logical :: is_even

    ! Initialization
    nx = size(x)
    is_even = mod(nx, 2) == 0
    nxfrm = compute_transform_length(nx)
    if (present(fac)) then
        f = fac
    else
        f = 1.0d0
    end if

    ! Quick Return
    if (f == 0.0d0) then
        cx = (0.0d0, 0.0d0)
        return
    end if

    ! Process
    if (f == 1.0d0) then
        cx(1) = cmplx(x(1), 0.0d0, real64)  ! the DC term is always real
        if (is_even) then
            do i = 2, nxfrm - 1
                cx(i) = cmplx(x(2*i-2), x(2*i-1), real64)
            end do
            cx(nxfrm) = cmplx(x(nx), 0.0d0, real64) ! always real for even-lengths
        else
            do i = 2, nxfrm
                cx(i) = cmplx(x(2*i-2), x(2*i-1), real64)
            end do
        end if
    else
        cx(1) = f * cmplx(x(1), 0.0d0, real64)  ! the DC term is always real
        if (is_even) then
            do i = 2, nxfrm - 1
                cx(i) = f * cmplx(x(2*i-2), x(2*i-1), real64)
            end do
            cx(nxfrm) = f * cmplx(x(nx), 0.0d0, real64) ! always real for even-lengths
        else
            do i = 2, nxfrm
                cx(i) = f * cmplx(x(2*i-2), x(2*i-1), real64)
            end do
        end if
    end if
end subroutine

! ------------------------------------------------------------------------------
! Rounds a number to the required precision (p), but rounds 0.5 down to the
! next lowest integer.
pure elemental function round(x, p) result(rst)
    ! Arguments
    real(real64), intent(in) :: x
    integer(int32), intent(in) :: p
    real(real64) :: rst

    ! Local Variables
    real(real64) :: scale, val

    ! Process
    scale = 10.0d0**(-p)
    val = x * scale + 0.49d0
    rst = floor(val) / scale
end function

! ------------------------------------------------------------------------------
pure function difference(x) result(rst)
    !! Computes the difference between each element in an array.
    real(real64), intent(in) :: x(:)
        !! The N-element array on which to operate.
    real(real64), allocatable :: rst(:)
        !! The N-1 element array containing the difference between each element
        !! in x.

    ! Local Variables
    integer(int32) :: i, n

    ! Initialization
    n = size(x)
    if (n == 0) return

    ! Process
    if (n == 1) then
        allocate(rst(1), source = 0.0d0)
        return
    end if
    allocate(rst(n - 1))
    do i = 1, n - 1
        rst(i) = x(i+1) - x(i)
    end do
end function

! ------------------------------------------------------------------------------
pure function cumulative_sum(x) result(rst)
    !! Computes the cumulative sum of an array.
    real(real64), intent(in) :: x(:)
        !! The N-element array on which to operate.
    real(real64), allocatable :: rst(:)
        !! An N-element array containing the cumulative sum of each element
        !! in x (e.g. cumulative_sum(x) = [x(1), x(1) + x(2), ...]).

    ! Local Variables
    integer(int32) :: i, n
    
    ! Initialization
    n = size(x)
    allocate(rst(n))

    ! Process
    rst(1) = x(1)
    do i = 2, n
        rst(i) = rst(i-1) + x(i)
    end do
end function

! ------------------------------------------------------------------------------
subroutine unwrap(x, tol)
    !! Shifts phase angle arrays to deal with jumps greater than or equal
    !! to tol by adding multiples of \(  \pm 2 \pi \) until
    !! the jump is less than tol.
    real(real64), intent(inout) :: x(:)
        !! On input, the phase array.  On output, the unwrapped phase array.
    real(real64), intent(in), optional :: tol
        !! The tolerance value.

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i, n
    real(real64), allocatable, dimension(:) :: dp, dpcorr
    real(real64) :: cutoff

    ! Initialization
    n = size(x)
    if (present(tol)) then
        cutoff = tol
    else
        cutoff = pi
    end if

    ! Compute the incremental variations
    dp = difference(x)

    ! Compute how many time 2*pi we're off, and round to the nearest integer 
    ! with the tie-breaker rounding n+0.5 down to n.
    dpcorr = round(dp / (2.0d0 * pi), -1)

    ! Stop the jump from happening if dp < cutoff
    do i = 1, n - 1
        if (abs(dp(i)) < cutoff) dpcorr(i) = 0.0d0
    end do

    ! Apply the corrections
    x(2:n) = x(2:n) - (2.0d0 * pi) * cumulative_sum(dpcorr)
end subroutine

! ------------------------------------------------------------------------------
pure function compute_overlap_segment_count(n, winsize) result(rst)
    !! Computes the number of overlapped signals using a nominal 50% overlap.
    integer(int32), intent(in) :: n
        !! The total length of the signal being overlapped.
    integer(int32), intent(in) :: winsize
        !! The window size.
    integer(int32) :: rst
        !! The number of segments.

    ! Local Variables
    integer(int32) :: nxfrm

    ! Process
    nxfrm = compute_transform_length(winsize)
    rst = max((n - 1) / nxfrm, 1)
end function

! ------------------------------------------------------------------------------
pure subroutine overlap(x, seg, winsize, buffer)
    !! Extracts a segment from a signal allowing for a nominally 50% overlap.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the entire signal.
    integer(int32), intent(in) :: seg
        !! The one-based index of the segment to extract.
    integer(int32), intent(in) :: winsize
        !! The size of the window (segment).  If this value is less than N, the 
        !! end of the segment will be padded with zeros.
    real(real64), intent(out) :: buffer(:)
        !! A winsize array where the segment will be written.

    ! Local Variables
    integer(int32) :: k, nx, nxfrm, nk, noff, i1, i2
    real(real64) :: del
    
    ! Initialization
    nx = size(x)
    nxfrm = compute_transform_length(winsize)
    nk = (nx - 1) / nxfrm
    if (nk > 1) then
        del = (nx - winsize) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if

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
end subroutine

! ------------------------------------------------------------------------------
end module