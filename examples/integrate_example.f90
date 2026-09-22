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
    integer(int32) :: i
    real(real64), parameter :: dt = 1.0d-3
    real(real64) :: t(npts), x(npts), v(npts), act(npts), noise(npts)

    ! Plot Variables
    type(plot_2d) :: plt
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd
    type(plot_data_2d) :: pd1, pd2, pd3

    ! Build the signal and add some noise
    t = (/ (i * dt, i = 0, npts - 1) /)
    v = cos(5.0d0 * t)
    act = 0.2d0 * sin(5.0d0 * t)    ! actual solution
    call random_number(noise)
    noise = 0.1d0 * (noise - 0.5d0)
    v = v + noise

    ! Compute the integral numerically
    x = integrate(dt, v)    ! initial value is 0

! ------------------------------------------------------------------------------
    ! Plot
    call plt%initialize()
    call plt%set_use_y2_axis(.true.)
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()

    call xAxis%set_title("t")
    call yAxis%set_title("dx/dt")
    call y2Axis%set_title("x(t)")
    call lgnd%set_is_visible(.true.)

    ! Plot the data
    call pd1%define_data(t, v)
    call pd1%set_name("dx/dt")
    call plt%push(pd1)

    call pd2%define_data(t, x)
    call pd2%set_draw_against_y2(.true.)
    call pd2%set_name("Numerical")
    call plt%push(pd2)

    call pd3%define_data(t, act)
    call pd3%set_draw_against_y2(.true.)
    call pd3%set_name("Analytical")
    call pd3%set_line_style(LINE_DASHED)
    call plt%push(pd3)

    call plt%draw()
end program