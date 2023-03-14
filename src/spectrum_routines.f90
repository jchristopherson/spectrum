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
end submodule