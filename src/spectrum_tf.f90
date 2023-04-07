submodule (spectrum) spectrum_tf
    use fftpack
contains
! ------------------------------------------------------------------------------
! REF: 
! - https://github.com/giuliovv/tfest/blob/main/tfest/tfest.py
! - https://dsp.stackexchange.com/questions/71811/understanding-the-h1-and-h2-estimators
! - https://github.com/epezent/etfe/blob/main/include/ETFE.hpp
module function siso_xfrm(win, x, y, etype, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:), y(:)
    integer(int32), intent(in), optional :: etype
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable :: rst(:)

    ! Local Variables
    integer(int32) :: est
    complex(real64), allocatable, dimension(:) :: pcross
    real(real64), allocatable, dimension(:) :: pwr
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(etype)) then
        est = etype
    else
        est = SPCTRM_H1_ESTIMATOR
    end if
    if (est /= SPCTRM_H1_ESTIMATOR .and. est /= SPCTRM_H2_ESTIMATOR) then
        est = SPCTRM_H1_ESTIMATOR
    end if

    ! Process
    select case (est)
    case (SPCTRM_H1_ESTIMATOR)
        pcross = csd(win, x, y, err = errmgr)
        if (errmgr%has_error_occurred()) return
        pwr = psd(win, x, err = errmgr)
        if (errmgr%has_error_occurred()) return
        rst = pcross / pwr
    case (SPCTRM_H2_ESTIMATOR)
        pcross = csd(win, y, x, err = errmgr)
        if (errmgr%has_error_occurred()) return
        pwr = psd(win, y, err = errmgr)
        if (errmgr%has_error_occurred()) return
        rst = pwr / pcross
    end select
end function

! ------------------------------------------------------------------------------
end submodule