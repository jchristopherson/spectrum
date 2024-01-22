# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    spectrum
    LANGUAGES Fortran
    VERSION 1.0.1
)

# Get helper macros and functions
include("${PROJECT_SOURCE_DIR}/cmake/helper.cmake")

# Confgiure everything
add_subdirectory(configure)

# Deal with the dependencies
find_package(BLAS)
find_package(LAPACK)
add_subdirectory(dependencies)

# Source
add_subdirectory(src)
add_fortran_library(
    ${PROJECT_NAME}
    ${PROJECT_INCLUDE_DIR}
    ${CMAKE_INSTALL_INCLUDEDIR}
    ${PROJECT_VERSION}
    ${PROJECT_VERSION_MAJOR}
    ${SPECTRUM_SOURCES}
)
target_link_libraries(
  ${PROJECT_NAME} 
  ${BLAS_LIBRARIES} 
  ${LAPACK_LIBRARIES}
)
link_library(${PROJECT_NAME} ${ferror_LIBRARY} ${ferror_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${fftpack_LIBRARY} ${fftpack_INCLUDE_DIR})
link_library(${PROJECT_NAME} ${linalg_LIBRARY} ${linalg_INCLUDE_DIR})

# Installation
add_subdirectory(install)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build tests: ${BUILD_TESTING}")
if (BUILD_TESTING)
    enable_testing()
    add_subdirectory(test)
endif()

# Examples
option(BUILD_SPECTRUM_EXAMPLES "Build SPECTRUM examples")
message(STATUS "Build SPECTRUM examples: ${BUILD_SPECTRUM_EXAMPLES}")
if (BUILD_SPECTRUM_EXAMPLES)
    add_subdirectory(examples)
endif()