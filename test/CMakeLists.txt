include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

add_subdirectory(fortran_test_helper)

set(spectrum_test_sources
    spectrum_tests.f90
    spectrum_psd_tests.f90
)

add_executable(spectrum_tests ${spectrum_test_sources})
link_library(spectrum_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
link_library(spectrum_tests ${fortran_test_helper_LIBRARY} ${fortran_test_helper_INCLUDE_DIR})
add_test(
    NAME spectrum_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:spectrum_tests>
)