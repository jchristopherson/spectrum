# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(SPECTRUM_SOURCES
    ${dir}/spectrum.f90
    ${dir}/spectrum_windows.f90
    ${dir}/spectrum_routines.f90
    ${dir}/spectrum_overlap.f90
    ${dir}/spectrum_psd.f90
    ${dir}/spectrum_csd.f90
    ${dir}/spectrum_fft.f90
    ${dir}/spectrum_convolve.f90
    ${dir}/spectrum_filter.f90
    ${dir}/spectrum_tf.f90
    ${dir}/spectrum_periodogram.f90
    ${dir}/spectrum_diff.f90
)
set(SPECTRUM_SOURCES ${SPECTRUM_SOURCES} PARENT_SCOPE)