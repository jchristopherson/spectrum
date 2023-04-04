include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Get FPLOT
add_subdirectory(fplot)

# Spectrogram Example
add_executable(spectrogram_example spectrogram_example.f90)
link_library(spectrogram_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(spectrogram_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Gaussian Filter Example
add_executable(gaussian_filter_example gaussian_filter_example.f90)
link_library(gaussian_filter_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(gaussian_filter_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})

# Overlap Example
add_executable(overlap_example overlap_example.f90)
link_library(overlap_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(overlap_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})

# PSD Example
add_executable(psd_example psd_example.f90)
link_library(psd_example ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(psd_example ${fplot_LIBRARY} ${fplot_INCLUDE_DIR})