module spectrum_tf_tests
    use iso_fortran_env
    use ieee_arithmetic
    use spectrum
    use fortran_test_helper
    implicit none
contains
! ------------------------------------------------------------------------------
    function test_siso_transfer_function() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: n = 256
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(n), y(n)
        complex(real64), allocatable :: tf(:)
        type(hamming_window) :: win

        ! Initialization
        rst = .true.
        win%size = n

        do i = 1, n
            x(i) = sin(2.0d0 * pi * real(i, real64) / real(n, real64))
        end do
        y = x

        tf = siso_transfer_function(win, x, y)

        if (size(tf) <= 0) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_transfer_function -1"
            return
        end if
        if (any(.not.ieee_is_finite(real(tf))) .or. any(.not.ieee_is_finite(aimag(tf)))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_transfer_function -2"
            return
        end if
        if (maxval(abs(tf)) <= 0.5d0) then
            rst = .false.
            print "(A)", "TEST FAILED: test_siso_transfer_function -3"
            return
        end if
    end function

! ------------------------------------------------------------------------------
    function test_mimo_transfer_function() result(rst)
        ! Arguments
        logical :: rst

        ! Parameters
        integer(int32), parameter :: n = 256
        integer(int32), parameter :: nout = 2
        integer(int32), parameter :: nin = 2
        real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

        ! Local Variables
        integer(int32) :: i
        real(real64) :: x(n, nin), y(n, nout)
        complex(real64), allocatable :: tf(:,:,:)
        type(hamming_window) :: win

        ! Initialization
        rst = .true.
        win%size = n

        ! Use orthogonal tones for each input channel so the diagonal entries are
        ! expected to be the dominant transfer terms and the off-diagonal terms
        ! remain near zero.
        do i = 1, n
            x(i,1) = sin(2.0d0 * pi * real(i, real64) / real(n, real64))
            x(i,2) = sin(2.0d0 * pi * 2.0d0 * real(i, real64) / real(n, real64))
        end do

        y(:,1) = x(:,1)
        y(:,2) = x(:,2)

        ! Compute the transfer matrix.
        tf = mimo_transfer_function(win, x, y)

        ! Validate the output shape and finiteness.
        if (size(tf, 1) /= nout .or. size(tf, 2) /= nin .or. size(tf, 3) /= n / 2 + 1) then
            rst = .false.
            print "(A)", "TEST FAILED: test_mimo_transfer_function -1"
            return
        end if
        if (any(.not.ieee_is_finite(real(tf))) .or. any(.not.ieee_is_finite(aimag(tf)))) then
            rst = .false.
            print "(A)", "TEST FAILED: test_mimo_transfer_function -2"
            return
        end if

        ! The diagonal terms should be finite and nontrivial for a matched
        ! MIMO identity mapping that preserves each input channel.
        if (maxval(abs(tf(1,1,:))) < 0.5d0 .or. maxval(abs(tf(2,2,:))) < 0.5d0) then
            rst = .false.
            print "(A)", "TEST FAILED: test_mimo_transfer_function -3"
            return
        end if

    end function
! ------------------------------------------------------------------------------
end module
