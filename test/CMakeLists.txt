set(spectrum_test_sources
    spectrum_tests.f90
)

add_executable(spectrum_tests ${spectrum_test_sources})
link_library(spectrum_tests ${PROJECT_NAME} ${PROJECT_INCLUDE_DIR})
add_test(
    NAME spectrum_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:spectrum_tests>
)