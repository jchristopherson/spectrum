module spectrum_integrate
    use iso_fortran_env
    use ferror
    use spectrum_errors
    implicit none
    private
    public :: integrate

contains
! ------------------------------------------------------------------------------
function integrate(dt, x, iv, err) result(rst)
    !! Integrates a data set.
    real(real64), intent(in) :: dt
        !! The sample interval.
    real(real64), intent(in), dimension(:) :: x
        !! The data set to integrate.
    real(real64), intent(in), optional :: iv
        !! The initial value.  The default value is zero.
    class(errors), intent(in), optional, target :: err
        !! An optional errors-based object that if provided can
        !! be used to retrieve information relating to any errors encountered 
        !! during execution.  If not provided, a default implementation of the 
        !! errors class is used internally to provide error handling.  Possible 
        !! errors and warning messages that may be encountered are as follows.
        !!
        !!  - SPCTRM_MEMORY_ERROR: Occurs if a memory allocation error occurs.
    real(real64), allocatable, dimension(:) :: rst
        !! The integrated data set.

    ! Local Variables
    integer(int32) :: i, n, flag
    real(real64) :: init_val
    class(errors), pointer :: errmgr
    type(errors), target :: deferr
    
    ! Initialization
    if (present(err)) then
        errmgr => err
    else
        errmgr => deferr
    end if
    if (present(iv)) then
        init_val = iv
    else
        init_val = 0.0d0
    end if
    n = size(x)
    allocate(rst(n), stat = flag)
    if (flag /= 0) then
        call errmgr%report_error("integrate", "Memory allocation error.", &
            SPCTRM_MEMORY_ERROR)
        return
    end if

    ! Process
    !
    ! Utilize the fixed step-size Adams-Bashforth methods.
    rst(1) = init_val
    if (n > 2) rst(2) = rst(1) + dt * x(1)     ! Euler method
    if (n > 3) rst(3) = rst(2) + 0.5d0 * dt * (3.0d0 * x(2) - x(1))
    if (n > 4) rst(4) = rst(3) + (dt / 1.2d1) * (2.3d1 * x(3) - 1.6d1 * x(2) + &
        5.0d0 * x(1))
    if (n > 5) rst(5) = rst(4) + (dt / 2.4d1) * (5.5d1 * x(4) - 5.9d1 * x(3) + &
        3.7d1 * x(2) - 9.0d0 * x(1))
    if (n > 6) then
        do i = 6, n
            rst(i) = rst(i - 1) + (dt / 7.2d2) * (1.901d3 * x(i - 1) - &
                2.774d3 * x(i - 2) + 2.616d3 * x(i - 3) - 1.274d3 * x(i - 4) + &
                2.51d2 * x(i - 5))
        end do
    end if
end function

! ------------------------------------------------------------------------------

! ------------------------------------------------------------------------------
end module