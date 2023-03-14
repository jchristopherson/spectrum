submodule (spectrum) spectrum_overlap
    use fftpack
contains
! ------------------------------------------------------------------------------
subroutine buffer_segment(win, xnew, xfrmsum, counter, work, cwork, err)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: xnew(:)
    real(real64), intent(inout) :: xfrmsum(:)
    integer(int32), intent(inout) :: counter
    real(real64), intent(inout), target :: work(:)
    complex(real64), intent(out) :: cwork(:)
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: i, j, n, flag
    real(real64) :: sumw, w, fac

    ! Initialization
    n = win%size
end subroutine

! ------------------------------------------------------------------------------
end submodule