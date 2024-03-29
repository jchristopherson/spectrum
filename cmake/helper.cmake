# helper.cmake
#
# A collection of macros and functions making life with CMake and Fortran a
# bit simpler.

# Use to include and export headers
function(include_headers lib dir install_dir)
    target_include_directories(
        ${lib}
        INTERFACE
        $<BUILD_INTERFACE:${dir}>
        $<INSTALL_INTERFACE:${install_dir}>
    )
endfunction()

# Use instead of add_library.
function(add_fortran_library lib_name mod_dir include_install_dir version major)
    add_library(${lib_name} ${ARGN})
    set_target_properties(
        ${lib_name}
        PROPERTIES
            POSITION_INDEPENDENT_CODE TRUE
            OUTPUT_NAME ${lib_name}
            VERSION ${version}
            SOVERSION ${major}
            Fortran_MODULE_DIRECTORY ${include_install_dir}
    )
    target_include_directories(
        ${lib_name}
        PUBLIC
        $<BUILD_INTERFACE:${mod_dir}>
        $<INSTALL_INTERFACE:${include_install_dir}>
    )
endfunction()

# Installs the library
function(install_library lib_name lib_install_dir bin_install_dir mod_dir install_dir)
    install(
        TARGETS ${lib_name}
        EXPORT ${lib_name}Targets
        RUNTIME DESTINATION ${bin_install_dir}
        LIBRARY DESTINATION ${lib_install_dir}
        ARCHIVE DESTINATION ${lib_install_dir}
        INCLUDES DESTINATION ${install_dir}/include
    )
    install(
        DIRECTORY ${mod_dir}
        DESTINATION ${install_dir}
    )
endfunction()

# Install the documentation files
function(install_documentation doc_dir install_dir)
    install(
        DIRECTORY ${doc_dir}
        DESTINATION ${install_dir}
    )
endfunction()

# Links the supplied library
function(link_library targ lib include_dir)
    target_link_libraries(${targ} ${lib})
    target_include_directories(${targ} PUBLIC $<BUILD_INTERFACE:${include_dir}>)
endfunction()

# Links with OpenMP
function(link_openmp targ)
    find_package(OpenMP)
    if (OpenMP_Fortran_FOUND)
        target_link_libraries(${targ} OpenMP::OpenMP_Fortran)
        target_compile_definitions(${targ} PUBLIC ${OpenMP_Fortran_FLAGS})
    endif()
endfunction()

# Enable parallelism with Do-Concurrent
function(enable_parallel_do targ nproc)
    if (CMAKE_Fortran_COMPILER_ID MATCHES "^Intel")
        # Intel - just requires OpenMP
        message(STATUS "IFORT: DO CONCURRENT will be parallelized.")
        link_openmp(${targ})
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        message(STATUS "GFORTRAN: DO CONCURRENT will be parallelized.")
        target_compile_definitions(${targ} PUBLIC "-ftree-parallelize-loops=${nproc}")
    else()
        # Unrecognized Compiler
        message(STATUS "Unrecognized compiler.  DO CONCURRENT not guaranteed to be parallelized.")
    endif()
endfunction()

# ------------------------------------------------------------------------------
# Helpful Macros
macro(print_all_variables)
    message(STATUS "---------- CURRENTLY DEFINED VARIABLES -----------")
    get_cmake_property(varNames VARIABLES)
    foreach(varName ${varNames})
        message(STATUS ${varName} = ${${varName}})
    endforeach()
    message(STATUS "---------- END ----------")
endmacro()
