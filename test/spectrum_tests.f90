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
program test
    use spectrum_psd_tests
    use spectrum_convolution_tests
    use spectrum_diff_tests
    use spectrum_integrate_tests
    use spectrum_filter_tests
    use spectrum_resample_tests
    use spectrum_tf_tests
    use spectrum_fft_tests
    implicit none

    ! Local Variables
    integer(int32) :: flag
    logical :: local

    ! Initialization
    flag = 0

    ! Tests
    local = test_psd()
    if (.not.local) flag = 1

    local = test_periodogram()
    if (.not.local) flag = 2

    local = test_csd()
    if (.not.local) flag = 3

    local = test_spectral_endpoint_scaling()
    if (.not.local) flag = 4

    local = test_convolution()
    if (.not.local) flag = 5

    local = test_spectrogram()
    if (.not.local) flag = 6

    local = test_stft_scaling()
    if (.not.local) flag = 7

    local = test_finite_difference()
    if (.not.local) flag = 8

    local = test_stencil_diff()
    if (.not.local) flag = 9

    local = test_stencil_diff_2()
    if (.not.local) flag = 10

    local = test_tvr_derivative_sparse()
    if (.not.local) flag = 11

    local = test_integrate()
    if (.not.local) flag = 12

    local = test_integrate_boundaries()
    if (.not.local) flag = 13

    local = test_sinc_filter()
    if (.not.local) flag = 14

    local = test_resample()
    if (.not.local) flag = 14

    local = test_filter_frequency_response()
    if (.not.local) flag = 16

    local = test_design_iir_filter()
    if (.not.local) flag = 17

    local = test_butterworth_filter_order()
    if (.not.local) flag = 18

    local = test_design_fir_filter()
    if (.not.local) flag = 19

    local = test_filter_boundaries()
    if (.not.local) flag = 20

    local = test_siso_transfer_function()
    if (.not.local) flag = 21

    local = test_mimo_transfer_function()
    if (.not.local) flag = 21

    local = test_irfft_odd_input_size()
    if (.not.local) flag = 22

    local = test_irfft_even_input_size()
    if (.not.local) flag = 23

    local = test_rfft_even_length()
    if (.not.local) flag = 24

    local = test_rfft_odd_length()
    if (.not.local) flag = 25

    local = test_rfft_padded_length()
    if (.not.local) flag = 26

    local = test_stft_even_window()
    if (.not.local) flag = 27

    local = test_stft_odd_window()
    if (.not.local) flag = 28

    local = test_stft_odd_input_size_guard()
    if (.not.local) flag = 29

    ! Output
    stop flag
end program