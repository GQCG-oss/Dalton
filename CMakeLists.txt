cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# suppress Windows warning about policy
cmake_policy(SET CMP0010 OLD)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

project(DALTON Fortran C)

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"        OFF)
option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"  OFF)
option(ENABLE_UTILITIES      "Enable utilities"              OFF)
option(ENABLE_LOCAL_DISKS    "Enable local disks"            OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking" OFF)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/binary-info
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/parallel-environment
    ${CMAKE_SOURCE_DIR}/cmake/testing
    )

set(LIBS)

include(Sources)

include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(ConfigExternal)
include(ConfigMath)
include(ConfigParallelEnvironment)
include(ConfigSafeGuards)
include(ConfigTesting)
include(ConfigVersion)
include(GenericMacros)
include(BinaryInfo)

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/mm
    ${PROJECT_SOURCE_DIR}/dft
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

# set definitions

add_definitions(-DBUILD_GEN1INT)
add_definitions(-DVAR_MFDS)
add_definitions(-D_FILE_OFFSET_BITS=64)
add_definitions(-DIMPLICIT_NONE)
if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DVAR_INT64)
endif()
# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

set(WORK_MEM_WORDS
    "64000000"
    CACHE STRING
    "Work memory in words"
    )
add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})

set(2EL_MEM_WORDS
    "1"
    CACHE STRING
    "Static memory for storing 2-el integrals"
    )
add_definitions(-DINSTALL_MMWORK=${2EL_MEM_WORDS})

# find a good default scratch base directory
foreach(dir "/global/work/$ENV{USER}" "/work" "/scratch" "/scr" "/temp" "/tmp")
   if(EXISTS "${dir}/")
       set(DEFAULT_SCRATCH_BASE_DIR ${dir})
       break()
   endif()
endforeach()

# set scratch directory
set(SCRATCH_DIR
    "${DEFAULT_SCRATCH_BASE_DIR}/DALTON_scratch_$ENV{USER}"
    CACHE STRING
    "Default scratch directory"
    )
message("-- SCRATCH_DIR set to: ${SCRATCH_DIR}")

set(BASIS_DIR
    "${CMAKE_SOURCE_DIR}/basis"
    CACHE STRING
    "Basis set directory"
    FORCE
    )

add_subdirectory(gen1int)

configure_script(
    ${CMAKE_SOURCE_DIR}/dalton.in
    ${CMAKE_BINARY_DIR}/dalton
    )

add_library(
    dalton
    ${C_SOURCES}
    ${FREE_FORTRAN_SOURCES}
    ${FIXED_FORTRAN_SOURCES}
    )

add_dependencies(dalton gen1int_interface)

add_executable(
    dalton.x
    ${CMAKE_SOURCE_DIR}/abacus/dalton.F
    )

target_link_libraries(
    dalton.x
    dalton
    gen1int_interface
    ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
    ${LIBS}
    )

# get size of static allocations
add_custom_target(
    info
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_static_size.py dalton.x
    )

include(ConfigPackaging)
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR})
