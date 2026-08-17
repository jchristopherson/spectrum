set(spectrum_test_sources
    spectrum_tests.f90
    spectrum_psd_tests.f90
    spectrum_convolution_tests.f90
    spectrum_diff_tests.f90
    spectrum_integrate_tests.f90
    spectrum_filter_tests.f90
    spectrum_resample_tests.f90
    spectrum_tf_tests.f90
)

find_package(fortran_test_helper QUIET)
if (NOT fortran_test_helper_FOUND)
    include(FetchContent)
    FetchContent_Declare(
        fortran_test_helper
        GIT_TAG "origin/main"
        GIT_REPOSITORY https://github.com/jchristopherson/fortran_test_helper
    )
    FetchContent_MakeAvailable(fortran_test_helper)
    set(FTH_LIBRARIES fortran_test_helper)
else()
    set(FTH_LIBRARIES fortran_test_helper::fortran_test_helper)
endif()

add_executable(spectrum_tests ${spectrum_test_sources})
target_link_libraries(
    spectrum_tests
    spectrum
    ${FTH_LIBRARIES}
)
add_test(
    NAME spectrum_tests
    WORKING_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}
    COMMAND $<TARGET_FILE:spectrum_tests>
)

if (${BUILD_SHARED_LIBS} AND WIN32)
    add_custom_command(
        TARGET spectrum_tests
        POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy 
        $<TARGET_RUNTIME_DLLS:spectrum_tests> 
        $<TARGET_FILE_DIR:spectrum_tests>
        COMMAND_EXPAND_LISTS
    )
endif()
