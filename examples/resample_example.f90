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
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)
    real(real64), parameter :: tmax = 1.0d0
    real(real64), parameter :: freq = 5.0d0
    real(real64), parameter :: fs1 = 40.0d0
    integer(int32), parameter :: nu = 4
    integer(int32), parameter :: nd = 2
    integer(int32) :: i, n
    real(real64) :: dt
    real(real64), allocatable, dimension(:) :: t, tu, td, x, xu, xd

    ! Plot Variables
    type(plot_2d) :: plt
    type(plot_data_2d) :: pd1, pd2, pd3
    class(plot_axis), pointer :: xAxis, yAxis
    class(legend), pointer :: lgnd

    ! Build the signal
    dt = 1.0d0 / fs1
    n = floor(tmax / dt) + 1
    allocate(t(n), x(n), tu(nu * n), td(n * nu / nd), xu(nu * n), xd(n * nu / nd))
    t = (/ (dt * i, i = 0, n - 1) /)
    x = sin(2.0d0 * pi * freq * t)

    ! Upsample the signal by a factor of 4
    xu = upsample(nu, fs1, x)
    tu = (/ (i * dt / nu, i = 0, n * nu - 1) /)

    ! Downsample the signal
    xd = downsample(nd, fs1, xu)
    td = (/ (i * dt * nd / nu, i = 0, n * nu / nd - 1) /)

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
    call pd1%set_draw_markers(.true.)
    call pd1%set_marker_style(MARKER_ASTERISK)
    call plt%push(pd1)

    call pd2%define_data(tu, xu)
    call pd2%set_name("Upsampled")
    call pd2%set_draw_markers(.true.)
    call pd2%set_marker_style(MARKER_EMPTY_NABLA)
    call plt%push(pd2)

    call pd3%define_data(td, xd)
    call pd3%set_name("Downsampled")
    call pd3%set_draw_markers(.true.)
    call pd3%set_marker_style(MARKER_EMPTY_SQUARE)
    call plt%push(pd3)

    call plt%draw()
end program