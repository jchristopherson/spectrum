program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Parameters
    integer(int32), parameter :: n = 256
    integer(int32), parameter :: nin = 2
    integer(int32), parameter :: nout = 2
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i
    real(real64) :: x(n, nin), y(n, nout), freq(n / 2 + 1), df
    complex(real64), allocatable :: tf(:,:,:)
    type(hamming_window) :: win

    ! Plot Variables
    type(multiplot) :: plt
    type(plot_2d) :: p1, p2, p3, p4
    type(plot_data_2d) :: pd

    ! Build two input channels with distinct tones
    do i = 1, n
        x(i,1) = sin(2.0d0 * pi * real(i, real64) / real(n, real64))
        x(i,2) = sin(2.0d0 * pi * 2.0d0 * real(i, real64) / real(n, real64))
    end do

    ! Make the outputs identical to their corresponding inputs
    y(:,1) = x(:,1)
    y(:,2) = x(:,2)

    ! Estimate the MIMO transfer matrix
    win%size = n
    tf = mimo_transfer_function(win, x, y)

    ! Frequency axis using the default sampling interval
    df = frequency_bin_width(1.0d0 / (x(2,1) - x(1,1)), n)
    freq = (/ (df * i, i = 0, size(tf, 3) - 1) /)

    ! Plots
    call plt%initialize(2, 2)
    call p1%initialize()
    call p2%initialize()
    call p3%initialize()
    call p4%initialize()

    call pd%define_data(freq, abs(tf(1,1,:)))
    call pd%set_line_width(2.0)
    call p1%push(pd)
    call plt%set(1, 1, p1)

    call pd%define_data(freq, abs(tf(2,1,:)))
    call p2%push(pd)
    call plt%set(2, 1, p2)

    call pd%define_data(freq, abs(tf(1,2,:)))
    call p3%push(pd)
    call plt%set(1, 2, p3)

    call pd%define_data(freq, abs(tf(2,2,:)))
    call p4%push(pd)
    call plt%set(2, 2, p4)

    call plt%draw()
end program
