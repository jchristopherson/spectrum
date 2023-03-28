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