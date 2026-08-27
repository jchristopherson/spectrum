# Master CMAKE Build Script
cmake_minimum_required(VERSION 3.24)
project(
    spectrum
    LANGUAGES Fortran
    VERSION 1.4.1
)
set(CMAKE_Fortran_STANDARD 2018)
set(CMAKE_Fortran_STANDARD_REQUIRED TRUE)

# Confgiure everything
add_subdirectory(configure)

# Source
add_subdirectory(src)

# Install
install(
  EXPORT ${PROJECT_NAME}-targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}"
)

# Testing
option(BUILD_TESTING "Build tests")
include(CTest)
message(STATUS "Build SPECTRUM tests: ${BUILD_TESTING}")
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