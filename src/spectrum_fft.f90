submodule (spectrum) spectrum_fft
    use fftpack
contains
! REFERENCES:
! - https://en.wikipedia.org/wiki/Short-time_Fourier_transform
! - https://en.wikipedia.org/wiki/Overlap%E2%80%93add_method
!
! Overlap-Add References:
! - https://ccrma.stanford.edu/~jos/sasp/Example_Overlap_Add_Convolution.html
! - https://www.mathworks.com/matlabcentral/fileexchange/41338-block-convolution-using-overlap-add-method
! - https://github.com/vishwas1234567/Octave-scilab/blob/master/ovelapadd%20and%20ovelap%20save%20in%20matlab.m

! Other Options:
! - Frequency Domain Windowing:
!   - https://www.embedded.com/dsp-tricks-frequency-domain-windowing/


module function stft(win, x, offsets, err) result(rst)
    ! Arguments
    class(window), intent(in) :: win
    real(real64), intent(in) :: x(:)
    integer(int32), intent(out), optional, allocatable :: offsets(:)
    class(errors), intent(inout), optional, target :: err
    complex(real64), allocatable :: rst(:,:)

    ! Local Variables
    integer(int32) :: i, j, k, m, nx, nxfrm, nend, nk, lwork, flag, i1
    real(real64) :: fac, w, sumw, del, scale
    real(real64), allocatable :: work(:), buffer(:,:)
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
    m = win%size          ! # of rows in the output matrix (transform length)
    nxfrm = compute_transform_length(m)
    if (mod(m, 2) == 0) then
        nend = nxfrm - 1
        scale = 2.0d0 / m
    else
        nend = nxfrm
        scale = 2.0d0 / (m - 1.0d0)
    end if
    lwork = 2 * m + 15
    nk = (nx - 1) / nxfrm
    if (nk > 1) then
        del = (nx - m) / (nk - 1.0d0)
    else
        del = 0.0d0
    end if

    ! Input Checking
    if (m > nx) go to 20

    ! Memory Allocation
    allocate(rst(nxfrm, nk), stat = flag)
    if (flag == 0) allocate(buffer(m, nk), stat = flag, source = 0.0d0)
    if (flag == 0) allocate(work(lwork), stat = flag)
    if (present(offsets)) then
        if (flag == 0) allocate(offsets(nk), stat = flag, source = 0)
    end if
    if (flag /= 0) go to 10

    ! Initialize the trig terms for the transforms.  As all transforms are the 
    ! same length, this array can be shared.
    call dffti(m, work)

    ! Compute each transform
!$OMP PARALLEL PRIVATE(i1, j, k, w, sumw, fac)
!$OMP DO
    do i = 1, nk
        ! Locate the portion of the signal on which to operate
        i1 = int((i - 1) * del + 0.5d0, int32)
        if (present(offsets)) offsets(i) = i1

        ! Apply the window
        j = 0
        sumw = 0.0d0
        do k = 1, m
            w = win%evaluate(j)
            j = j + 1
            sumw = sumw + w
            buffer(k,i) = w * x(i1+k)
        end do
        fac = m / sumw

        ! Compute the transform
        call dfftf(m, buffer(:,i), work)

        ! Scale the transform
        rst(1,i) = fac * scale * cmplx(buffer(1,i), 0.0d0, real64)
        do k = 2, nend
            rst(k,i) = fac * scale * cmplx(buffer(2*k-1,i), &
                buffer(2*k,i), real64)
        end do
        if (nend /= nxfrm) rst(nxfrm,i) = fac * scale * &
            cmplx(buffer(m,i), 0.0d0, real64)
    end do
!$OMP END DO
!$OMP END PARALLEL

    ! End
    return

    ! Memory Error Handling
10  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 100) "Memory allocation error flag ", flag, "."
    call errmgr%report_error("stft", trim(errmsg), &
        SPCTRM_MEMORY_ERROR)
    return

    ! Window size exceeds signal length
20  continue
    allocate(character(len = 256) :: errmsg)
    write(errmsg, 101) "Window size (", m, ") exceeds the signal length (", &
        nx, ")."
    call errmgr%report_error("stft", trim(errmsg), SPCTRM_INVALID_INPUT_ERROR)
    return

    ! Formatting
100 format(A, I0, A)
101 format(A, I0, A, I0, A)
end function

end submodule