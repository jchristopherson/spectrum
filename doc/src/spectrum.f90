module spectrum
    use spectrum_windows
    use spectrum_routines
    use spectrum_periodogram
    use spectrum_fft
    use spectrum_convolve
    use spectrum_filter
    use spectrum_tf
    use spectrum_diff
    use spectrum_integrate
    use spectrum_resample
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
    public :: stft

    ! SPECTRUM_CONVOLVE.F90
    public :: convolve
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

    ! SPECTRUM_TF.F90
    public :: siso_transfer_function
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
 
end module