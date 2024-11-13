! REF:
! - https://ejde.math.txstate.edu/conf-proc/21/k3/knowles.pdf
! - https://arxiv.org/pdf/2009.01911.pdf

module spectrum_diff
    use iso_fortran_env
    use blas
    use linalg
    use ferror
    use spectrum_errors
    implicit none
    private
    public :: finite_difference
    public :: tvr_derivative

    interface finite_difference
        module procedure :: finite_difference_1
        module procedure :: finite_difference_2
    end interface

    interface finite_difference_driver
        module procedure :: finite_difference_driver_1
        module procedure :: finite_difference_driver_2
    end interface

    ! BLAS Routines:
    interface
        ! subroutine dgbmv(trans, m, n, kl, ku, alpha, a, lda, x, incx, beta, y, incy)
        !     use iso_fortran_env, only : int32, real64
        !     character, intent(in) :: trans
        !     integer(int32), intent(in) :: m, n, kl, ku, lda, incx, incy
        !     real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
        !     real(real64), intent(inout) :: y(*)
        ! end subroutine

        ! subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        !     use iso_fortran_env, only : int32, real64
        !     character, intent(in) :: transa, transb
        !     integer(int32), intent(in) :: m, n, k, lda, ldb, ldc
        !     real(real64), intent(in) :: alpha, beta, a(lda,*), b(ldb,*)
        !     real(real64), intent(inout) :: c(ldc,*)
        ! end subroutine

        ! subroutine dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        !     use iso_fortran_env, only : int32, real64
        !     character, intent(in) :: trans
        !     integer(int32), intent(in) :: m, n, lda, incx, incy
        !     real(real64), intent(in) :: alpha, beta, a(lda,*), x(*)
        !     real(real64), intent(inout) :: y(*)
        ! end subroutine

        subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
            use iso_fortran_env, only : int32, real64
            integer(int32), intent(in) :: n, nrhs, lda, ldb
            real(real64), intent(inout) :: a(lda,*), b(ldb,*)
            integer(int32), intent(out) :: ipiv(*), info
        end subroutine
    end interface
contains
! ------------------------------------------------------------------------------
pure subroutine finite_difference_driver_1(dt, x, dxdt)
    ! Arguments
    real(real64), intent(in) :: dt
    real(real64), intent(in), dimension(:) :: x
    real(real64), intent(out), dimension(:) :: dxdt

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    if (n == 0) return
    if (n == 1) then
        dxdt = 0.0d0
        return
    end if
    dxdt(1) = (x(2) - x(1)) / dt
    do i = 2, n - 1
        dxdt(i) = 0.5d0 * (x(i + 1) - x(i - 1)) / dt
    end do
    dxdt(n) = (x(n) - x(n - 1)) / dt
end subroutine

! ------------------------------------------------------------------------------
pure subroutine finite_difference_driver_2(t, x, dxdt)
    ! Arguments
    real(real64), intent(in), dimension(:) :: t, x
    real(real64), intent(out), dimension(:) :: dxdt

    ! Local Variables
    integer(int32) :: i, n

    ! Process
    n = size(x)
    if (n == 0) return
    if (n == 1) then
        dxdt = 0.0d0
        return
    end if
    dxdt(1) = (x(2) - x(1)) / (t(2) - t(1))
    do i = 2, n - 1
        dxdt(i) = 0.5d0 * (x(i + 1) - x(i - 1)) / (t(i + 1) - t(i - 1))
    end do
    dxdt(n) = (x(n) - x(n - 1)) / (t(n) - t(n - 1))
end subroutine

! ------------------------------------------------------------------------------
function finite_difference_1(dt, x, err) result(rst)
    !! Estimates the derivative of a data set by means of a naive 
    !! implementation of a finite difference scheme based upon central 
    !! differences.
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the derivative estimate.

    ! Local Variables
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(x)
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("finite_difference_1", &
            "Memory allocation error.", SPCTRM_MEMORY_ERROR)
        return
    end if
    
    ! Process
    call finite_difference_driver(dt, x, rst)
end function

! ------------------------------------------------------------------------------
function finite_difference_2(t, x, err) result(rst)
    !! Computes an estimate to the derivative of an evenly-sampled data
    !! set using total variation regularization.
    real(real64), intent(in), dimension(:) :: t
        !! An N-element array containing the time points at which x was sampled.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_ARRAY_SIZE_MISMATCH_ERROR: Occurs if t and x are not the
        !!      same size.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the derivative estimate.

    ! Local Variables
    integer(int32) :: n, flag
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    n = size(t)

    ! Input Checking
    if (size(x) /= n) then
        call errmgr%report_error("finite_difference_2", &
            "The input arrays must be the same size.", SPCTRM_ARRAY_SIZE_ERROR)
        return
    end if

    ! Memory Allocation
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("finite_difference_2", &
            "Memory allocation error.", SPCTRM_MEMORY_ERROR)
        return
    end if

    ! Process
    call finite_difference_driver(t, x, rst)
end function

! ******************************************************************************
! TOTAL VARIATION REGULARIZATION
! ------------------------------------------------------------------------------
! REF: https://oliver-k-ernst.medium.com/how-to-differentiate-noisy-signals-2baf71b8bb65
! https://github.com/smrfeld/Total-Variation-Regularization-Derivative-Python/blob/main/python/diff_tvr.py
! https://github.com/florisvb/PyNumDiff/blob/master/pynumdiff/total_variation_regularization/__chartrand_tvregdiff__.py

! Constructs the N-by-N+1 D matrix:
!          | -1     1         |
! D = 1/dx |  0     -1      1 |
!          |  0     0      -1 |
subroutine make_d_full(dx, d)
    ! Arguments
    real(real64), intent(in) :: dx
    real(real64), intent(out), dimension(:,:) :: d

    ! Local Variables
    integer(int32) :: j, n
    real(real64) :: idx

    ! Process
    n = size(d, 1)
    idx = 1.0d0 / dx
    d = 0.0d0
    do j = 1, n + 1
        if (j > 1) d(j-1,j) = idx
        if (j <= n) d(j,j) = -idx
    end do
end subroutine

! ------------------------------------------------------------------------------
! Constructs the N-by-N+1 A matrix.
subroutine make_a_full(dx, a)
    ! Arguments
    real(real64), intent(in) :: dx
    real(real64), intent(out), dimension(:,:) :: a

    ! Local Variables
    integer(int32) :: j, n
    real(real64) :: hdx
    
    ! Process
    n = size(a, 1)
    hdx = 0.5d0 * dx
    do j = 1, n + 1
        if (j == 1) then
            a(:,j) = hdx
        else
            if (j > 1) a(j-1,j) = hdx
            if (j <= n) a(j:,j) = dx
        end if
    end do
end subroutine

! ------------------------------------------------------------------------------
! Constructs the N-by-N E matrix.  The matrix is a diagonal matrix with only the
! diagonal stored.
subroutine make_e(d, u, e)
    ! Arguments
    real(real64), intent(in), dimension(:,:) :: d
    real(real64), intent(in), dimension(:) :: u
    real(real64), intent(out), dimension(:) :: e

    ! Local Variables
    integer(int32) :: j, n, n1
    real(real64) :: eps
    
    ! Process
    eps = sqrt(epsilon(eps))
    n = size(d, 1)
    n1 = n + 1
    call dgemv("N", n, n1, 1.0d0, d, n, u, 1, 0.0d0, e, 1)
    do j = 1, n
        e(j) = 1.0d0 / sqrt(e(j)**2 + eps)
    end do
end subroutine

! ------------------------------------------------------------------------------
subroutine tvr_diff_small(alpha, dt, x, maxiter, dxdt, tol, niter, err)
    real(real64), intent(in) :: alpha ! variational parameter
    real(real64), intent(in) :: dt  ! time step
    real(real64), intent(in), dimension(:) :: x ! data array to differentiate
    integer(int32), intent(in) :: maxiter  ! max # of iterations
    real(real64), intent(out), dimension(:) :: dxdt ! derivative dx/dt
    real(real64), intent(in) :: tol ! tolerance on change in gradient
    integer(int32), intent(out) :: niter ! # of iterations taken
    class(errors), intent(inout) :: err

    ! Local Variables
    integer(int32) :: i, n, n1, flag
    integer(int32), allocatable, dimension(:) :: ipiv
    real(real64) :: offset, nrm, nrmold
    real(real64), allocatable, dimension(:,:) :: d, a, dte, l, ata, h
    real(real64), allocatable, dimension(:) :: e, u, atb, atau, lu, g

    ! Initialization
    n = size(x)
    n1 = n + 1
    nrmold = huge(nrmold)
    
    ! Memory Allocations
    allocate( &
        d(n, n1), &
        a(n, n1), &
        e(n), &
        u(n1), &
        atb(n), &
        dte(n1, n), &
        l(n1, n1), &
        ata(n1, n1), &
        atau(n1), &
        lu(n1), &
        g(n1), &
        h(n1, n1), &
        ipiv(n1), &
        stat = flag)
    if (flag /= 0) go to 10

    ! Construct matrices
    call make_d_full(dt, d)
    call make_a_full(dt, a)
    call dgemm("T", "N", n1, n1, n, 1.0d0, a, n, a, n, 0.0d0, ata, n1) ! A**T * A

    ! Provide a first estimate of the derivative
    u(1) = 0.0d0
    u(2:) = finite_difference(dt, x)

    ! Precompute A**T * (X(1) - X)
    call dgemv("T", n, n1, 1.0d0, a, n, offset - x, 1, 0.0d0, atb, 1)

    ! Iteration Process
    do i = 1, maxiter
        ! Compute E and L
        call make_e(d, u, e)
        call diag_mtx_mult(.false., .true., dt, e, d, 0.0d0, dte) ! dt * D**T * E
        call dgemm("N", "N", n1, n1, n, 1.0d0, dte, n1, d, n, 0.0d0, l, n1) ! L = (dx * D**T * E) * D

        ! Compute the gradient
        call dgemv("N", n1, n1, 1.0d0, ata, n1, u, 1, 0.0d0, atau, 1)
        call dgemv("N", n1, n1, alpha, l, n1, u, 1, 0.0d0, lu, 1)
        g = atau + atb + lu

        ! Compute H
        h = ata + alpha * l

        ! Solve H * s = g, for s - stored in g
        call dgesv(n1, 1, h, n1, ipiv, g, n1, flag)
        if (flag /= 0) go to 20

        ! Check the solution
        nrm = norm2(g)
        if (abs(nrm - nrmold) < tol) exit
        nrmold = nrm

        ! Update the derivative estimate
        u = u - g
    end do
    niter = min(i, maxiter)

    ! Extract the computed derivative
    dxdt = u(1:n)

    ! End
    return

    ! Memory Error
10  continue
    call err%report_error("tvr_diff_small", "Memory allocation error.", &
        SPCTRM_MEMORY_ERROR)
    return

    ! Solution Error - Singular matrix
20  continue
    call err%report_error("tvr_diff_small", "A singular Hessian matrix " // &
        "was encountered.  Check to ensure the problem is properly defined.", &
        SPCTRM_SINGULAR_MATRIX_ERROR)
    return
end subroutine

! ------------------------------------------------------------------------------
function tvr_derivative(dt, x, alpha, maxiter, tol, niter, err) result(rst)
    !! Computes an estimate to the derivative of an evenly-sampled data
    !! set using total variation regularization.
    !!
    !! See Also
    !!
    !! - van Breugel, Floris & Brunton, Bingni & Kutz, J.. (2020). Numerical 
    !!   differentiation of noisy data: A unifying multi-objective optimization 
    !!   framework. 
    !!
    !! - Oliver K. Ernst, Ph. D. (2021, February 16). How to differentiate 
    !!   noisy signals. Medium. https://oliver-k-ernst.medium.com/how-to-differentiate-noisy-signals-2baf71b8bb65 
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is
        !! to be estimated.
    real(real64), intent(in) :: alpha
        !! The regularization parameter.
    integer(int32), intent(in), optional :: maxiter
        !! The maximum number of iterations to allow.  The default is 20 
        !! iterations.
    real(real64), intent(in), optional :: tol
        !! The convergence tolerance to use.  The tolerance is 
        !! applied to the difference in Euclidean norms of the derivative update
        !! vector.  Once the norm of the update vector is changing less than 
        !! this tolerance, the iteration process will terminate.  The default
        !! is 1e-3.
    integer(int32), intent(out), optional :: niter
        !! The number of iterations actually performed.
    class(errors), intent(inout), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
        !!
        !!  - SPCTRM_SINGULAR_MATRIX_ERROR: Occurs if the internal Hessian 
        !!      estimate becomes singular.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the estimate of the derivative.

    ! Local Variables
    integer(int32) :: mi, n, flag, ni
    real(real64) :: gtol
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(maxiter)) then
        mi = maxiter
    else
        mi = 20
    end if
    if (present(tol)) then
        gtol = tol
    else
        gtol = 1.0d-3
    end if
    n = size(x)
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("tvr_derivative", "Memory allocation error.", &
            SPCTRM_MEMORY_ERROR)
        return
    end if

    ! Process
    call tvr_diff_small(alpha, dt, x, mi, rst, gtol, ni, errmgr)
    if (present(niter)) niter = ni
end function

! ------------------------------------------------------------------------------
end module