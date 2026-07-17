# Locate the source directory
set(dir ${CMAKE_CURRENT_SOURCE_DIR})

# Define the source files
add_subdirectory(fftpack)
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
    ${dir}/spectrum_integrate.f90
    ${dir}/spectrum_resample.f90
    ${dir}/fftpack/fftpack.f90
    ${dir}/fftpack/dfftb.f90
    ${dir}/fftpack/dfftf.f90
    ${dir}/fftpack/dffti.f90
    ${dir}/fftpack/rfftb1.f90
    ${dir}/fftpack/rfftf1.f90
    ${dir}/fftpack/rffti1.f90
    ${dir}/fftpack/zfftb.f90
    ${dir}/fftpack/zffti.f90
    ${dir}/fftpack/zfftf.f90
    ${FFTPACK_SOURCES}
)

# Dependencies
include(FetchContent)
find_package(ferror QUIET)
find_package(linalg 1.8.10 QUIET)

if (NOT ferror_FOUND)
    message(STATUS "FERROR could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        ferror
        GIT_REPOSITORY "https://github.com/jchristopherson/ferror"
    )
    FetchContent_MakeAvailable(ferror)
    set(FERROR_LIBRARIES ferror)
    set(BUILD_FERROR TRUE)
else()
    set(FERROR_LIBRARIES ferror::ferror)
    set(BUILD_FERROR FALSE)
endif()

if (NOT linalg_FOUND)
    message(STATUS "LINALG could not be found.  A reference version will be employed.")
    FetchContent_Declare(
        linalg
        GIT_REPOSITORY "https://github.com/jchristopherson/linalg"
    )
    FetchContent_MakeAvailable(linalg)
    set(LINALG_LIBRARIES linalg)
    set(BUILD_LINALG TRUE)
else()
    set(LINALG_LIBRARIES linalg::linalg)
    set(BUILD_LINALG FALSE)
endif()

# Build
add_library(${PROJECT_NAME} ${SPECTRUM_SOURCES})
target_link_libraries(
    ${PROJECT_NAME}
    PRIVATE
    ${FERROR_LIBRARIES}
    ${LINALG_LIBRARIES}
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    WINDOWS_EXPORT_ALL_SYMBOLS ON
)

set(SPECTRUM_MOD_DIR ${CMAKE_CURRENT_BINARY_DIR}/mod_files/)
if (NOT EXISTS "${SPECTRUM_MOD_DIR}")
    file(MAKE_DIRECTORY "${SPECTRUM_MOD_DIR}")
endif()

set_target_properties(
    ${PROJECT_NAME} PROPERTIES
    Fortran_MODULE_DIRECTORY ${SPECTRUM_MOD_DIR}
)
target_include_directories(
    ${PROJECT_NAME} PUBLIC
    $<BUILD_INTERFACE:${SPECTRUM_MOD_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_MODULEDIR}>
)

# Install
install(
    TARGETS ${PROJECT_NAME}
    EXPORT ${PROJECT_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
)
install(DIRECTORY ${SPECTRUM_MOD_DIR} DESTINATION ${CMAKE_INSTALL_MODULEDIR})

if (${BUILD_SHARED_LIBS})
    if (${BUILD_FERROR})
        install(TARGETS ${FERROR_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
else()
    if (${BUILD_FERROR})
        install(TARGETS ${FERROR_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
    if (${BUILD_LINALG})
        install(TARGETS ${LINALG_LIBRARIES}
            RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
            LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        )
    endif()
endif()