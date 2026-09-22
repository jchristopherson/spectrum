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
    integer(int32), parameter :: hn = npts / 2
    real(real64), parameter :: slope = 1.0d0
    real(real64), parameter :: intercept = 1.0d0
    real(real64), parameter :: maxt = 1.0d0
    integer(int32) :: niter
    real(real64) :: dt, t(npts), x(npts), noise(npts), xf(npts,2)

    ! Plot Variables
    type(plot_2d) :: plt
    class(plot_axis), pointer :: xAxis, yAxis, y2Axis
    class(legend), pointer :: lgnd
    type(plot_data_2d) :: pd1, pd2, pd3

    ! Build the signal and add some noise
    t = linspace(0.0d0, maxt, npts)
    dt = t(2) - t(1)
    call random_number(noise)
    noise = 0.05d0 * (noise - 0.5d0)
    x(:hn) = slope * t(:hn)
    x(hn+1:) = -slope * t(hn+1:) + intercept
    x = x + noise

    ! Compute the derivative
    xf = filter_diff(dt, x, 0.025d0 / dt)

! ------------------------------------------------------------------------------
    ! Plot
    call plt%initialize()
    call plt%set_use_y2_axis(.true.)
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()
    y2Axis => plt%get_y2_axis()
    lgnd => plt%get_legend()
    call xAxis%set_title("t")
    call yAxis%set_title("x(t)")
    call y2Axis%set_title("dx/dt")
    call lgnd%set_is_visible(.true.)

    call pd1%define_data(t, x)
    call pd1%set_name("Original")
    call plt%push(pd1)

    call pd2%define_data(t, xf(:,1))
    call pd2%set_name("Filtered")
    call pd2%set_line_width(2.0)
    call plt%push(pd2)

    call pd3%define_data(t, xf(:,2))
    call pd3%set_draw_against_y2(.true.)
    call pd3%set_name("Derivative")
    call plt%push(pd3)

    call plt%draw()
end program