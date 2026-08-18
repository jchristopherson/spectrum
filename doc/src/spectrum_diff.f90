! REF:
! - https://ejde.math.txstate.edu/conf-proc/21/k3/knowles.pdf
! - https://arxiv.org/pdf/2009.01911.pdf

module spectrum_diff
    use iso_fortran_env
    use blas
    use linalg
    use lapack, only : dgesv
    implicit none
    private
    public :: finite_difference
    public :: tvr_derivative
    public :: stencil_diff_5
    public :: stencil_second_diff_5
    public :: filter_diff

    interface finite_difference
        module procedure :: finite_difference_1
        module procedure :: finite_difference_2
    end interface

    interface finite_difference_driver
        module procedure :: finite_difference_driver_1
        module procedure :: finite_difference_driver_2
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
        dxdt(i) = (x(i + 1) - x(i - 1)) / (t(i + 1) - t(i - 1))
    end do
    dxdt(n) = (x(n) - x(n - 1)) / (t(n) - t(n - 1))
end subroutine

! ------------------------------------------------------------------------------
pure function finite_difference_1(dt, x) result(rst)
    !! Estimates the derivative of a data set by means of a naive 
    !! implementation of a finite difference scheme based upon central 
    !! differences.
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the derivative estimate.

    ! Local Variables
    integer(int32) :: n
    
    ! Initialization
    n = size(x)
    allocate(rst(n))
    
    ! Process
    call finite_difference_driver(dt, x, rst)
end function

! ------------------------------------------------------------------------------
pure function finite_difference_2(t, x) result(rst)
    !! Computes an estimate to the derivative of an evenly-sampled data
    !! set using total variation regularization.
    real(real64), intent(in), dimension(:) :: t
        !! An N-element array containing the time points at which x was sampled.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the derivative estimate.

    ! Local Variables
    integer(int32) :: n
    
    ! Initialization
    n = size(t)

    ! Input Checking
    if (size(x) /= n) return

    ! Memory Allocation
    allocate(rst(n))

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
pure subroutine make_d_full(dx, d)
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
pure subroutine make_a_full(dx, a)
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
function tvr_diff_small(alpha, dt, x, maxiter, tol, niter) result(dxdt)
    real(real64), intent(in) :: alpha ! variational parameter
    real(real64), intent(in) :: dt  ! time step
    real(real64), intent(in), dimension(:) :: x ! data array to differentiate
    integer(int32), intent(in) :: maxiter  ! max # of iterations
    real(real64), intent(in) :: tol ! tolerance on change in gradient
    integer(int32), intent(out) :: niter ! # of iterations taken
    real(real64), allocatable, dimension(:) :: dxdt ! derivative dx/dt

    ! Local Variables
    integer(int32) :: i, n, n1, flag
    integer(int32), allocatable, dimension(:) :: ipiv
    real(real64) :: offset, nrm, nrmold
    real(real64), allocatable, dimension(:,:) :: d, a, dte, l, ata, h
    real(real64), allocatable, dimension(:) :: e, u, atb, atau, lu, g

    ! Initialization
    n = size(x)
    n1 = n + 1
    offset = x(1)
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
        dxdt(n) &
    )

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
        if (flag /= 0) return

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
end function

! ------------------------------------------------------------------------------
function tvr_derivative(dt, x, alpha, maxiter, tol, use_sparse, niter) result(rst)
    !! Computes an estimate to the derivative of an evenly-sampled data
    !! set using total variation regularization.
    !!
    !! This implementation solves the augmented dense formulation using an
    !! integration matrix and a first-difference regularization operator. The
    !! dense formulation uses the physical time step in its difference matrix
    !! and allocates dense matrices whose storage grows quadratically with the
    !! number of samples.
    !!
    !! When use_sparse is true, this routine dispatches to
    !! tvr_derivative_sparse. That solver uses a different data-fitting term,
    !! a second-difference operator, physical curvature scaling, and normalized
    !! IRLS weights. Consequently, alpha values are solver-specific and should
    !! not be compared directly.
    !!
    !! See Also
    !!
    !! - van Breugel, Floris & Brunton, Bingni & Kutz, J.. (2020). Numerical 
    !!   differentiation of noisy data: A unifying multi-objective optimization 
    !!   framework. 
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
        !! applied to the change in the update measure. The dense solver uses
        !! an absolute Euclidean norm, while the sparse solver uses a relative
        !! Euclidean norm. The default is 1e-3.
    logical, intent(in), optional :: use_sparse
        !! True if the sparse solver should be used vs. the dense solver.  This 
        !! is highly recommended when N is larger than ~1000.  The default is 
        !! true such that the sparse solver is used. The sparse solver has
        !! linear storage growth and is preferred for large data sets.
    integer(int32), intent(out), optional :: niter
        !! The number of iterations actually performed.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the estimate of the derivative.

    ! Local Variables
    logical :: sparse
    integer(int32) :: mi, ni
    real(real64) :: gtol
    
    ! Initialization
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
    if (present(use_sparse)) then
        sparse = use_sparse
    else
        sparse = .true.
    end if

    ! Process
    if (sparse) then
        rst = tvr_derivative_sparse(dt, x, alpha, maxiter = mi, tol = gtol, &
            niter = ni)
    else
        rst =  tvr_diff_small(alpha, dt, x, mi, gtol, ni)
    end if
    if (present(niter)) niter = ni
end function

! ------------------------------------------------------------------------------
function tvr_derivative_sparse(dt, x, alpha, maxiter, tol, niter) result(rst)
    !! Computes a total-variation-regularized derivative using sparse matrices.
    !!
    !! This routine regularizes the finite-difference derivative with a sparse
    !! second-difference operator and uses the derivative itself as the data
    !! fit. The second difference is scaled by dt**2 to represent physical
    !! curvature, then normalized before the IRLS weights are formed to keep
    !! the sparse systems well conditioned. The resulting iteratively
    !! reweighted systems contain at most five non-zero diagonals, so memory
    !! use grows linearly with the number of samples.
    !!
    !! This is a scalable alternative to tvr_derivative, not an algebraically
    !! equivalent sparse implementation of it. The two routines use different
    !! data-fitting and regularization formulations. The sparse alpha is a
    !! normalized regularization parameter and should not be compared directly
    !! with the dense alpha.
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is
        !! to be estimated.
    real(real64), intent(in) :: alpha
        !! The normalized regularization parameter. Larger values produce a
        !! smoother derivative estimate. Its value is not directly comparable
        !! with the alpha parameter used by the dense formulation.
    integer(int32), intent(in), optional :: maxiter
        !! The maximum number of reweighting iterations. The default is 20.
    real(real64), intent(in), optional :: tol
        !! The relative convergence tolerance applied to the derivative update
        !! norm. The default is 1e-3.
    integer(int32), intent(out), optional :: niter
        !! The number of iterations actually performed.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the estimate of the derivative.

    integer(int32) :: i, j, k, n, nrows, max_iterations, iterations
    real(real64) :: convergence_tol, change, curvature_scale, system_scale
    real(real64) :: data_coefficient, regularization_coefficient, dt2
    real(real64), allocatable :: estimate(:), candidate(:), weights(:), &
        second_difference(:), rhs(:)
    type(csr_matrix) :: d, h

    n = size(x)
    allocate(rst(n))
    if (n == 0) return
    if (n < 3 .or. dt == 0.0d0 .or. alpha < 0.0d0) then
        rst = finite_difference(dt, x)
        if (present(niter)) niter = 0
        return
    end if

    if (present(maxiter)) then
        max_iterations = maxiter
    else
        max_iterations = 20
    end if
    if (present(tol)) then
        convergence_tol = tol
    else
        convergence_tol = 1.0d-3
    end if
    if (max_iterations < 1) then
        rst = finite_difference(dt, x)
        if (present(niter)) niter = 0
        return
    end if

    nrows = n - 2
    allocate(estimate(n), candidate(n), rhs(n), weights(nrows), &
        second_difference(nrows))
    estimate = finite_difference(dt, x)
    rhs = estimate

    d = create_empty_csr_matrix(nrows, n, 3 * nrows)
    d%row_indices(1) = 1
    do i = 1, nrows
        k = 3 * (i - 1) + 1
        d%column_indices(k:k+2) = [i, i + 1, i + 2]
        d%values(k:k+2) = [1.0d0, -2.0d0, 1.0d0]
        d%row_indices(i + 1) = k + 3
    end do

    iterations = 0
    dt2 = dt**2
    do i = 1, max_iterations
        second_difference = matmul(d, estimate) / dt2
        curvature_scale = max(1.0d0, maxval(abs(second_difference)))
        weights = 1.0d0 / sqrt((second_difference / curvature_scale)**2 + &
            sqrt(epsilon(1.0d0)))

        ! Curvature is normalized before IRLS weighting. This preserves the
        ! relative weighting while avoiding large coefficients for flat data.
        system_scale = max(1.0d0, alpha * maxval(weights))
        data_coefficient = 1.0d0 / system_scale
        regularization_coefficient = alpha / system_scale

        h = create_empty_csr_matrix(n, n, 5 * n - 6)
        h%row_indices(1) = 1
        k = 1
        do j = 1, n
            do while (k <= h%row_indices(1) - 1)
                k = k + 1
            end do
            call fill_sparse_system_row(j, weights, data_coefficient, &
                regularization_coefficient, h, k)
            h%row_indices(j + 1) = k
        end do

        candidate = sparse_direct_solve(h, data_coefficient * rhs)
        change = norm2(candidate - estimate) / max(norm2(estimate), 1.0d0)
        estimate = candidate
        iterations = i
        if (change < convergence_tol) exit
    end do

    rst = estimate
    if (present(niter)) niter = iterations
end function

! ------------------------------------------------------------------------------
pure subroutine fill_sparse_system_row(row, weights, data_coefficient, &
    regularization_coefficient, h, offset)
    integer(int32), intent(in) :: row
    real(real64), intent(in) :: weights(:), data_coefficient, &
        regularization_coefficient
    type(csr_matrix), intent(inout) :: h
    integer(int32), intent(inout) :: offset

    integer(int32) :: col, first_row, last_row, difference_row
    real(real64) :: value

    first_row = max(1, row - 2)
    last_row = min(size(weights), row)
    do col = first_row, last_row + 2
        value = 0.0d0
        do difference_row = first_row, last_row
            value = value + second_difference_coefficient(difference_row, col) * &
                second_difference_coefficient(difference_row, row) * &
                weights(difference_row)
        end do
        value = value * regularization_coefficient
        if (col == row) value = value + data_coefficient
        h%column_indices(offset) = col
        h%values(offset) = value
        offset = offset + 1
    end do
end subroutine

! ------------------------------------------------------------------------------
pure function second_difference_coefficient(row, col) result(value)
    integer(int32), intent(in) :: row, col
    real(real64) :: value

    select case (col - row)
    case (0)
        value = 1.0d0
    case (1)
        value = -2.0d0
    case (2)
        value = 1.0d0
    case default
        value = 0.0d0
    end select
end function

! ******************************************************************************
! V1.1.2 ADDITIONS
! ------------------------------------------------------------------------------
pure function stencil_diff_5(dt, x) result(rst)
    !! Utilizes a 5-point stencil to estimate the derivative of a data set.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Five-point_stencil" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the derivative estimate.

    ! Local Variables
    integer(int32) :: i, n
    
    ! Initialization
    n = size(x)
    allocate(rst(n))

    ! Process
    ! Step in and out of the problem via finite differences; else, use
    ! a 5-point stencil of the form:
    !
    ! f'(x) = (-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)) / (12h)
    rst(1) = (x(2) - x(1)) / dt
    rst(2) = (x(3) - x(2)) / dt

    do i = 3, n - 2
        rst(i) = (-x(i + 2) + 8.0d0 * (x(i + 1) - x(i - 1)) + x(i - 2)) / &
            (12.0d0 * dt)
    end do

    rst(n-1) = (x(n-1) - x(n-2)) / dt
    rst(n) = (x(n) - x(n-1)) / dt
end function

! ------------------------------------------------------------------------------
pure function stencil_second_diff_5(dt, x) result(rst)
    !! Utilizes a 5-point stencil to estimate the second derivative of a data
    !! set.
    !!
    !! See Also
    !!
    !! - <a href="https://en.wikipedia.org/wiki/Five-point_stencil" target="_blank">Wikipedia</a>
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    real(real64), allocatable, dimension(:) :: rst
        !! An N-element array containing the derivative estimate.

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: h2
    
    ! Initialization
    n = size(x)
    allocate(rst(n))

    ! Process
    ! Step in and out of the problem via finite differences; else, use
    ! a 5-point stencil of the form:
    !
    ! f"(x) = (-f(x+2h) + 16f(x+h) - 30f(x) + 16f(x-h) - f(x-2h)) / (12h**2)
    h2 = dt**2
    rst(1) = (x(3) - 2.0d0 * x(2) + x(1)) / h2
    rst(2) = (x(4) - 2.0d0 * x(3) + x(2)) / h2

    do i = 3, n - 2
        rst(i) = (-x(i + 2) - 3.0d1 * x(i) + 1.6d1 * (x(i + 1) + x(i - 1)) - &
            x(i - 2)) / (1.2d1 * h2)
    end do

    rst(n - 1) = (x(n - 1) - 2.0d0 * x(n - 2) + x(n - 3)) / h2
    rst(n) = (x(n) - 2.0d0 * x(n - 1) + x(n - 2)) / h2
end function

! ******************************************************************************
! V1.1.3 ADDITIONS
! ------------------------------------------------------------------------------
pure function filter_diff(dt, x, fc) result(rst)
    !! Estimates the derivative of a signal by utilization of a second-order
    !! system as a filter.
    real(real64), intent(in) :: dt
        !! The time step between data points.
    real(real64), intent(in), dimension(:) :: x
        !! An N-element array containing the data whose derivative is to be 
        !! estimated.
    real(real64), intent(in) :: fc
        !! The filter cutoff frequency, in Hz.
    real(real64), allocatable, dimension(:,:) :: rst
        !! An N-element array containing the filtered signal in the first column
        !! and the derivative estimate in the second.

    ! Parameters
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: zeta = 0.5d0 * sqrt(2.0d0)

    ! Local Variables
    integer(int32) :: i, n
    real(real64) :: fs, wn
    
    ! Initialization
    n = size(x)
    fs = 1.0d0 / dt
    wn = 2.0d0 * pi * fc

    ! Input Checking
    if (fc >= 0.5d0 * fs .or. fc <= 0.0d0) return

    ! Memory Allocations
    allocate(rst(n,2))

    ! Define the initial conditions
    rst(1,1) = x(1)    ! output initial value is equivalent to the original value
    rst(1,2) = (x(2) - x(1)) / dt  ! finite difference estimate of the first point

    ! Perform the integration using Euler's method
    do i = 2, n
        ! Predictor Stage (Explicit Method)
        rst(i,:) = rst(i-1,:) + dt * fcn(rst(i-1,:), x(i-1), wn, zeta)

        ! Corrector Stage (Implicit Method)
        rst(i,:) = rst(i-1,:) + dt * fcn(Rst(i,:), x(i), wn, zeta)
    end do
end function

! ----------
pure function fcn(x, y, wn, zeta) result(dxdt)
    !! The second-order equations of motion.
    real(real64), intent(in) :: x(2)
        !! The current state variables.
    real(real64), intent(in) :: y
        !! The current value of the original signal.
    real(real64), intent(in) :: wn
        !! The second-order system natural frequency, in rad/s.
    real(real64), intent(in) :: zeta
        !! The second-order system damping ratio.
    real(real64) :: dxdt(2)
        !! The output derivative values.

    ! Equation of Motion:
    ! x" + 2 * zeta * wn * x' + wn**2 * x = wn**2 * y
    dxdt(1) = x(2)
    dxdt(2) = wn**2 * (y - x(1)) - 2.0d0 * zeta * wn * x(2)
end function

! ------------------------------------------------------------------------------
end module