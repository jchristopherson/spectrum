submodule (spectrum) spectrum_routines
contains
! ------------------------------------------------------------------------------
pure module function compute_xfrm_length_1(n) result(rst)
    integer(int32), intent(in) :: n
    integer(int32) :: rst
    if (mod(n, 2) == 0) then
        rst = n / 2 + 1
    else
        rst = (n + 1) / 2
    end if
end function

! ------------------------------------------------------------------------------
pure module elemental function frequency_bin_width_1(fs, n) result(rst)
    ! Arguments
    real(real64), intent(in) :: fs
    integer(int32), intent(in) :: n
    real(real64) :: rst

    ! Local Variables
    integer(int32) :: m

    ! Process
    m = compute_transform_length(n)
    rst = 0.5d0 * fs / (m - 1.0d0)
end function

! ------------------------------------------------------------------------------
pure module elemental function is_power_of_two_1(n) result(rst)
    ! Arguments
    integer(int32), intent(in) :: n
    logical :: rst

    ! Process
    rst = (n /= 0) .and. (iand(n, n - 1) == 0)
    ! This is equivalent to the C form : n != 0 && ((n & (n - 1)) == 0)
end function

! ------------------------------------------------------------------------------
pure module elemental function next_power_of_two_1(n) result(rst)
    ! Arguments
    integer(int32), intent(in) :: n
    integer(int32) :: rst

    ! Process
    rst = ceiling(log(real(n)) / log(2.0), int32)
end function

! ------------------------------------------------------------------------------
module subroutine unpack_real_transform(x, cx, fac)
    ! Arguments
    real(real64), intent(in) :: x(:)
    complex(real64), intent(out) :: cx(:)
    real(real64), intent(in), optional :: fac

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
pure function diff(x) result(rst)
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64), allocatable :: rst(:)

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
    ! Arguments
    real(real64), intent(in) :: x(:)
    real(real64), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    allocate(rst(n))
    rst(1) = x(1)
    do i = 2, n
        rst(i) = rst(i-1) + x(i)
    end do
end function

! ------------------------------------------------------------------------------
module subroutine unwrap_1(x, tol)
    ! Arguments
    real(real64), intent(inout) :: x(:)
    real(real64), intent(in), optional :: tol

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
    dp = diff(x)

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
end submodule