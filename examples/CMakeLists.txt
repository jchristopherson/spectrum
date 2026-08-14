# Get Dependencies
add_subdirectory(fortran-csv-module)
include(FetchContent)
FetchContent_Declare(
    fplot
    GIT_REPOSITORY "https://github.com/jchristopherson/fplot"
)
FetchContent_MakeAvailable(fplot)
set(fplot_LIBRARY fplot)

# Spectrogram Example
add_executable(spectrogram_example spectrogram_example.f90)
target_link_libraries(spectrogram_example spectrum ${fplot_LIBRARY})

# Gaussian Filter Example
add_executable(gaussian_filter_example gaussian_filter_example.f90)
target_link_libraries(gaussian_filter_example spectrum ${fplot_LIBRARY})

add_executable(gaussian_filter_example_2 gaussian_filter_example_2.f90)
target_link_libraries(gaussian_filter_example_2 spectrum ${fplot_LIBRARY})

# Overlap Example
add_executable(overlap_example overlap_example.f90)
target_link_libraries(overlap_example spectrum ${fplot_LIBRARY})

# PSD Example
add_executable(psd_example psd_example.f90)
target_link_libraries(psd_example spectrum ${fplot_LIBRARY})

# CSD Example
add_executable(csd_example csd_example.f90)
target_link_libraries(csd_example spectrum ${fplot_LIBRARY})

# Window Example
add_executable(window_example window_example.f90)
target_link_libraries(window_example spectrum ${fplot_LIBRARY})

# Hann Window Example
add_executable(hann_example hann_example.f90)
target_link_libraries(hann_example spectrum ${fplot_LIBRARY})

# Hamming Window Example
add_executable(hamming_example hamming_example.f90)
target_link_libraries(hamming_example spectrum ${fplot_LIBRARY})

# Welch Window Example
add_executable(welch_example welch_example.f90)
target_link_libraries(welch_example spectrum ${fplot_LIBRARY})

# Blackman-Harris Window Example
add_executable(blackman_example blackman_example.f90)
target_link_libraries(blackman_example spectrum ${fplot_LIBRARY})

# Flat-Top Window Example
add_executable(flat_top_example flat_top_example.f90)
target_link_libraries(flat_top_example spectrum ${fplot_LIBRARY})

# Rectangular Window Example
add_executable(rectangular_example rectangular_example.f90)
target_link_libraries(rectangular_example spectrum ${fplot_LIBRARY})

# Filter Example
add_executable(filter_example filter_example.f90)
target_link_libraries(filter_example spectrum ${fplot_LIBRARY})

# Total Variation Filter Example
add_executable(tv_filter_example tv_filter_example.f90)
target_link_libraries(tv_filter_example spectrum ${fplot_LIBRARY})

# SISO Transfer Function Example
add_executable(siso_tf_example siso_tf_example.f90)
target_link_libraries(siso_tf_example spectrum ${fplot_LIBRARY} ${fortran-csv-module_LIBRARY})
target_include_directories(siso_tf_example PUBLIC ${fortran-csv-module_INCLUDE_DIR})

# MIMO Transfer Function Example
add_executable(mimo_tf_example mimo_tf_example.f90)
target_link_libraries(mimo_tf_example spectrum ${fplot_LIBRARY})

# TVR Differentiation Example
add_executable(tvr_diff_example tvr_diff_example.f90)
target_link_libraries(tvr_diff_example spectrum ${fplot_LIBRARY})

# Integration Example
add_executable(integrate_example integrate_example.f90)
target_link_libraries(integrate_example spectrum ${fplot_LIBRARY})

# Remove Mean Example
add_executable(remove_mean_example remove_mean_example.f90)
target_link_libraries(remove_mean_example spectrum ${fplot_LIBRARY})

# Resample Example
add_executable(resample_example resample_example.f90)
target_link_libraries(resample_example spectrum ${fplot_LIBRARY})

# Filtered Derivative Example
add_executable(filtered_derivative_example filtered_derivative_example.f90)
target_link_libraries(filtered_derivative_example spectrum ${fplot_LIBRARY})

# Filter Comparison Example
add_executable(filter_comparison_example filter_comparison_example.f90)
target_link_libraries(filter_comparison_example spectrum ${fplot_LIBRARY})

# Filter Design Example
add_executable(filter_design_example filter_design_example.f90)
target_link_libraries(filter_design_example spectrum ${fplot_LIBRARY})

# IIR Filter Design Example
add_executable(iir_filter_design_example iir_filter_design_example.f90)
target_link_libraries(iir_filter_design_example spectrum ${fplot_LIBRARY})

#
if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET spectrogram_example
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:spectrogram_example> 
        $<TARGET_FILE_DIR:spectrogram_example>
        COMMAND_EXPAND_LISTS
    )
endif()