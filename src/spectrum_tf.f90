module spectrum_tf
    use iso_fortran_env
    use spectrum_periodogram
    use spectrum_routines
    use spectrum_windows
    implicit none
    private
    public :: siso_transfer_function
    public :: mimo_transfer_function
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
pure function siso_transfer_function(win, x, y, etype, nfft) result(rst)
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
    complex(real64), allocatable :: rst(:)
        !! Returns the complex-valued transfer function estimate.

    ! Local Variables
    integer(int32) :: est
    complex(real64), allocatable, dimension(:) :: pcross
    real(real64), allocatable, dimension(:) :: pwr
    
    ! Initialization
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
        pcross = csd(win, y, x, nfft = nfft)
        pwr = psd(win, x, nfft = nfft)
        allocate(rst(size(pcross)))
        rst = (0.0d0, 0.0d0)
        where (abs(pwr) > tiny(1.0d0))
            rst = pcross / pwr
        end where
    case (SPCTRM_H2_ESTIMATOR)
        pcross = csd(win, x, y, nfft = nfft)
        pwr = psd(win, y, nfft = nfft)
        allocate(rst(size(pcross)))
        rst = (0.0d0, 0.0d0)
        where (abs(pcross) > tiny(1.0d0))
            rst = pwr / pcross
        end where
    end select
end function

! ------------------------------------------------------------------------------
pure function mimo_transfer_function(win, x, y, etype, nfft) result(rst)
    !! Estimates the transfer function matrix for a multiple-input/
    !! multiple-output (MIMO) system.
    !!
    !! The input signal matrix is expected to have dimensions (n, nin) where
    !! the second index denotes the input channel.  The output signal matrix is
    !! expected to have dimensions (n, nout) where the second index denotes the
    !! output channel.  The returned transfer matrix has dimensions
    !! (nout, nin, m) where m is the transform length.
    class(window), intent(in) :: win
        !! The window object.
    real(real64), intent(in) :: x(:, :)
        !! An N-by-NIN array containing the input signals.
    real(real64), intent(in) :: y(:, :)
        !! An N-by-NOUT array containing the output signals.
    integer(int32), intent(in), optional :: etype
        !! An optional input that, if supplied, denotes the estimator to use.
        !! Supported values are the same as for SISO transfer estimates.
    integer(int32), intent(in), optional :: nfft
        !! An optional input that can be used to force the length of each
        !! individual DFT operation by padding any remaining space with zeros.
    complex(real64), allocatable :: rst(:, :, :)
        !! Returns the complex-valued transfer function estimate matrix.

    ! Local Variables
    integer(int32) :: est, nx, ny, nin, nout, i, j, nf
    complex(real64), allocatable, dimension(:) :: pcross
    real(real64), allocatable, dimension(:) :: pwr

    ! Initialization
    nx = size(x, 1)
    ny = size(y, 1)
    nin = size(x, 2)
    nout = size(y, 2)
    if (present(etype)) then
        est = etype
    else
        est = SPCTRM_H1_ESTIMATOR
    end if
    if (est /= SPCTRM_H1_ESTIMATOR .and. est /= SPCTRM_H2_ESTIMATOR) then
        est = SPCTRM_H1_ESTIMATOR
    end if
    if (nx /= ny .or. nin < 1 .or. nout < 1) return
    if (present(nfft)) then
        nf = nfft
    else
        nf = win%size
    end if
    allocate(rst(nout, nin, compute_transform_length(nf)))
    rst = (0.0d0, 0.0d0)

    ! Process
    select case (est)
    case (SPCTRM_H1_ESTIMATOR)
        do j = 1, nout
            do i = 1, nin
                pcross = csd(win, y(:, j), x(:, i), nfft = nfft)
                pwr = psd(win, x(:, i), nfft = nfft)
                where (abs(pwr) > tiny(1.0d0))
                    rst(j, i, :) = pcross / pwr
                elsewhere
                    rst(j, i, :) = (0.0d0, 0.0d0)
                end where
            end do
        end do
    case (SPCTRM_H2_ESTIMATOR)
        do j = 1, nout
            do i = 1, nin
                pcross = csd(win, x(:, i), y(:, j), nfft = nfft)
                pwr = psd(win, y(:, j), nfft = nfft)
                where (abs(pcross) > tiny(1.0d0))
                    rst(j, i, :) = pwr / pcross
                elsewhere
                    rst(j, i, :) = (0.0d0, 0.0d0)
                end where
            end do
        end do
    end select
end function

! ------------------------------------------------------------------------------
end module