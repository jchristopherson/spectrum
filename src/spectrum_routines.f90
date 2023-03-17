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
end submodule