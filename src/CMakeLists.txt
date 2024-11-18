# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
set(SPECTRUM_SOURCES
    ${dir}/spectrum.f90
    ${dir}/spectrum_windows.f90
    ${dir}/spectrum_routines.f90
    ${dir}/spectrum_fft.f90
    ${dir}/spectrum_convolve.f90
    ${dir}/spectrum_filter.f90
    ${dir}/spectrum_tf.f90
    ${dir}/spectrum_periodogram.f90
    ${dir}/spectrum_diff.f90
    ${dir}/spectrum_errors.f90
    ${dir}/spectrum_integrate.f90
)
set(SPECTRUM_SOURCES ${SPECTRUM_SOURCES} PARENT_SCOPE)