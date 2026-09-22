! MIT License
!
! Copyright (c) 2023-2026 Jason Christopherson
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!
program example
    use iso_fortran_env
    use spectrum
    use fplot_core
    implicit none

    ! Variables
    integer(int32), parameter :: npts = 1000
    real(real64), parameter :: dt = 1.0d-3
    integer(int32) :: i
    real(real64) :: t(npts), x(npts), xm(npts)

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Generate the signal
    t = (/ (i * dt, i = 0, npts - 1) /)
    call random_number(x)

    ! Remove the mean
    xm = remove_mean(x)

! ------------------------------------------------------------------------------
    ! Plot
    call plt%initialize()
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(t, x)
    call pd1%set_name("Original")
    call plt%push(pd1)

    call pd2%define_data(t, xm)
    call pd2%set_name("Mean Removed")
    call plt%push(pd2)

    call plt%draw()
end program