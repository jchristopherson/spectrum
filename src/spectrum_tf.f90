module spectrum_tf
    use iso_fortran_env
    use spectrum_periodogram
    use spectrum_windows
    use spectrum_errors
    use ferror
    implicit none
    private
    public :: siso_transfer_function
    public :: SPCTRM_H1_ESTIMATOR
    public :: SPCTRM_H2_ESTIMATOR

    integer(int32), parameter :: SPCTRM_H1_ESTIMATOR = 50002
        !! A flag for requesting an H1 transfer function estimator.  An H1 
        !! estimator is best used when noise is uncorrelated with the input, 
        !! and results in a transfer function estimate of the form
        !! 
        !! $$ H_{1} = \frac{P_{yx}}{P_{xx}} $$.
    integer(int32), parameter :: SPCTRM_H2_ESTIMATOR = 50003
        !! A flag for requesting an H2 transfer function estimator.  An H2 
        !! estimator is best used when noise is uncorrelated with the output, 
        !! and results in a transfer function estimate of the form
        !! 
        !! $$ H_{2} = \frac{P_{yy}}{P_{xy}} $$.

contains
! ------------------------------------------------------------------------------
! REF: 
! - https://github.com/giuliovv/tfest/blob/main/tfest/tfest.py
! - https://dsp.stackexchange.com/questions/71811/understanding-the-h1-and-h2-estimators
! - https://github.com/epezent/etfe/blob/main/include/ETFE.hpp
function siso_transfer_function(win, x, y, etype, nfft, err) result(rst)
    !! Estimates the transfer function for a single-input/single-output
    !! (SISO) system.
    class(window), intent(in) :: win
        !! The window object.
    real(real64), intent(in) :: x(:)
        !! An N-element array containing the input signal.
    real(real64), intent(in) :: y(:)
        !! An N-element array containing the output signal.
    integer(int32), intent(in), optional :: etype
        !! An optional input that, if supplied, denotes the
        !! estimator to use.  If no value is specified, an H1 estimator is used.
        !! The following options are supported.
        !!
        !!  - SPCTRM_H1_ESTIMATOR: Uses an H1 estimate.
        !!
        !!  - SPCTRM_H2_ESTIMATOR: Uses an H2 estimate.
        !!
        !! If an unrecognized value is provided, the routine defaults to an 
        !! H1 estimator.
    integer(int32), intent(in), optional :: nfft
        !! An optional input that can be used to force the length of each
        !! individual DFT operation by padding any remaining space with zeros.
        !! If not supplied, the window size is used to determine the size of
        !! the DFT.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_INVALID_INPUT_ERROR: Occurs if win is not sized 
        !!      appropriately.
    complex(real64), allocatable :: rst(:)
        !! Returns the complex-valued transfer function estimate.

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
        pcross = csd(win, y, x, nfft = nfft, err = errmgr)
        if (errmgr%has_error_occurred()) return
        pwr = psd(win, x, nfft = nfft, err = errmgr)
        if (errmgr%has_error_occurred()) return
        rst = pcross / pwr
    case (SPCTRM_H2_ESTIMATOR)
        pcross = csd(win, x, y, nfft = nfft, err = errmgr)
        if (errmgr%has_error_occurred()) return
        pwr = psd(win, y, err = errmgr)
        if (errmgr%has_error_occurred()) return
        rst = pwr / pcross
    end select
end function

! ------------------------------------------------------------------------------
end module