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
module spectrum
    use spectrum_windows
    use spectrum_routines
    use spectrum_periodogram
    use spectrum_fft
    use spectrum_convolve
    use spectrum_filter
    use spectrum_filter_design
    use spectrum_tf
    use spectrum_diff
    use spectrum_integrate
    use spectrum_resample
    use fftpack, only : zffti, zfftf, zfftb, dffti, dfftf, dfftb, fft, ifft, &
        fftshift, ifftshift
    implicit none
    private

    ! SPECTRUM_WINDOWS.F90
    public :: window
    public :: window_function
    public :: rectangular_window
    public :: hann_window
    public :: hamming_window
    public :: welch_window
    public :: blackman_harris_window
    public :: flat_top_window

    ! SPECTRUM_ROUTINES.F90
    public :: compute_transform_length
    public :: frequency_bin_width
    public :: is_power_of_two
    public :: next_power_of_two
    public :: unpack_real_transform
    public :: cumulative_sum
    public :: unwrap
    public :: difference
    public :: compute_overlap_segment_count
    public :: overlap
    public :: remove_mean

    ! SPECTRUM_PERIODOGRAM.F90
    public :: psd
    public :: csd
    public :: periodogram
    public :: cross_periodogram

    ! SPECTRUM_FFT.F90
    public :: rfft
    public :: irfft
    public :: stft
    public :: stft_result

    ! SPECTRUM_CONVOLVE.F90
    public :: convolve
    public :: deconvolve
    public :: SPCTRM_FULL_CONVOLUTION
    public :: SPCTRM_CENTRAL_CONVOLUTION

    ! SPECTRUM_FITLER.F90
    public :: gaussian_filter
    public :: tv_filter
    public :: filter
    public :: moving_average_filter
    public :: sinc_filter
    public :: LOW_PASS_FILTER
    public :: HIGH_PASS_FILTER
    public :: BAND_PASS_FILTER
    public :: BAND_STOP_FILTER

    ! SPECTRUM_FILTER_DESIGN.F90
    public :: butterworth_filter_order
    public :: design_fir_filter
    public :: design_iir_filter
    public :: filter_frequency_response

    ! SPECTRUM_TF.F90
    public :: siso_transfer_function
    public :: mimo_transfer_function
    public :: SPCTRM_H1_ESTIMATOR
    public :: SPCTRM_H2_ESTIMATOR

    ! SPECTRUM_DIFF.F90
    public :: finite_difference
    public :: tvr_derivative
    public :: stencil_diff_5
    public :: stencil_second_diff_5
    public :: filter_diff

    ! SPECTRUM_INTEGRATE.F90
    public :: integrate

    ! SPECTRUM_RESAMPLE.F90
    public :: upsample
    public :: downsample
 
    ! FFTPACK.F90
    public :: zffti
    public :: zfftf
    public :: zfftb
    public :: dffti
    public :: dfftf
    public :: dfftb
    public :: fft
    public :: ifft
    public :: fftshift
    public :: ifftshift
end module