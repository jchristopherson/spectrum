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

    ! Parameters
    integer(int32), parameter :: window_size = 512
    real(real64), parameter :: fs = 2048.0d0
    real(real64), parameter :: f0 = 1.0d2
    real(real64), parameter :: f1 = 1.0d3
    real(real64), parameter :: duration = 50.0d0
    real(real64), parameter :: pi = 2.0d0 * acos(0.0d0)

    ! Local Variables
    integer(int32) :: i, npts
    integer(int32), allocatable, dimension(:) :: offsets
    real(real64) :: k, df
    complex(real64), allocatable, dimension(:,:) :: rst
    real(real64), allocatable, dimension(:) :: t, x, f, s
    real(real64), allocatable, dimension(:,:) :: mag
    real(real64), allocatable, dimension(:,:,:) :: xy
    type(hann_window) :: win
    type(stft_result) :: z

    ! Plot Variables
    type(surface_plot) :: plt
    type(surface_plot_data) :: pd
    class(plot_axis), pointer :: xAxis, yAxis
    type(rainbow_colormap) :: map

    ! Create the exponential chirp signal
    npts = floor(duration * fs) + 1
    t = linspace(0.0d0, duration, npts)
    k = (f1 / f0)**(1.0 / duration)
    x = sin(2.0d0 * pi * f0 * (k**t - 1.0d0) / log(k))

    ! Determine sampling frequency parameters
    df = frequency_bin_width(fs, window_size)

    ! Define the window
    win%size = window_size

    ! Compute the spectrogram of x
    z = stft(win, x)
    rst = z%stft
    offsets = z%offsets

    ! Compute the magnitude, along with each frequency and time point
    mag = abs(rst)

    allocate(f(size(mag, 1)))
    f = (/ (df * i, i = 0, size(f) - 1) /)

    allocate(s(size(mag, 2)))
    do i = 1, size(s)
        if (i == 1) then
            s(i) = offsets(i) / fs
        else
            s(i) = i * (offsets(i) - offsets(i-1)) / fs
        end if
    end do
    xy = meshgrid(s, f)

    ! Plot the results
    call plt%initialize()
    call plt%set_colormap(map)
    call plt%set_use_map_view(.true.)
    xAxis => plt%get_x_axis()
    yAxis => plt%get_y_axis()

    call xAxis%set_title("Time [s]")
    call yAxis%set_title("Frequency [Hz]")
    call yAxis%set_autoscale(.false.)
    call yAxis%set_limits(0.0d0, f(size(mag, 1)))
    
    call pd%define_data(xy(:,:,1), xy(:,:,2), mag)
    call plt%push(pd)
    call plt%draw()
end program