include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get Dependencies
add_subdirectory(fplot)
add_subdirectory(fortran-csv-module)

# Spectrogram Example
add_executable(spectrogram_example spectrogram_example.f90)
link_library(spectrogram_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(spectrogram_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Gaussian Filter Example
add_executable(gaussian_filter_example gaussian_filter_example.f90)
link_library(gaussian_filter_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(gaussian_filter_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

add_executable(gaussian_filter_example_2 gaussian_filter_example_2.f90)
link_library(gaussian_filter_example_2 ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(gaussian_filter_example_2 ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Overlap Example
add_executable(overlap_example overlap_example.f90)
link_library(overlap_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(overlap_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(overlap_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# PSD Example
add_executable(psd_example psd_example.f90)
link_library(psd_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(psd_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# CSD Example
add_executable(csd_example csd_example.f90)
link_library(csd_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(csd_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Window Example
add_executable(window_example window_example.f90)
link_library(window_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(window_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(window_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Hann Window Example
add_executable(hann_example hann_example.f90)
link_library(hann_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(hann_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(hann_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Hamming Window Example
add_executable(hamming_example hamming_example.f90)
link_library(hamming_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(hamming_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(hamming_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Welch Window Example
add_executable(welch_example welch_example.f90)
link_library(welch_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(welch_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(welch_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Blackman-Harris Window Example
add_executable(blackman_example blackman_example.f90)
link_library(blackman_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(blackman_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(blackman_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Flat-Top Window Example
add_executable(flat_top_example flat_top_example.f90)
link_library(flat_top_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(flat_top_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(flat_top_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Rectangular Window Example
add_executable(rectangular_example rectangular_example.f90)
link_library(rectangular_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(rectangular_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(rectangular_example ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# Filter Example
add_executable(filter_example filter_example.f90)
link_library(filter_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(filter_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Total Variation Filter Example
add_executable(tv_filter_example tv_filter_example.f90)
link_library(tv_filter_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(tv_filter_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# SISO Transfer Function Example
add_executable(siso_tf_example siso_tf_example.f90)
link_library(siso_tf_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(siso_tf_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(siso_tf_example ${fortran-csv-module_LIBRARY} ${fortran-csv-module_INCLUDE_DIR})

# TVR Differentiation Example
add_executable(tvr_diff_example tvr_diff_example.f90)
link_library(tvr_diff_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(tvr_diff_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Integration Example
add_executable(integrate_example integrate_example.f90)
link_library(integrate_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(integrate_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Remove Mean Example
add_executable(remove_mean_example remove_mean_example.f90)
link_library(remove_mean_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(remove_mean_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Resample Example
add_executable(resample_example resample_example.f90)
link_library(resample_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(resample_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Filtered Derivative Example
add_executable(filtered_derivative_example filtered_derivative_example.f90)
link_library(filtered_derivative_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(filtered_derivative_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
