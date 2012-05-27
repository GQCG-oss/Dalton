cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(DALTON Fortran C)

set(THIS_CODE "trunk")
if(THIS_CODE STREQUAL "linsca")
    add_definitions(-DVAR_LINSCA)
endif()

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/binary-info
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/parallel-environment
    )

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE
        Debug
        CACHE STRING
        "Choose the type of build, options are: None Debug Release."
        FORCE
        )
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if(    NOT cmake_build_type_tolower STREQUAL "debug"
		AND NOT cmake_build_type_tolower STREQUAL "release"
		AND NOT cmake_build_type_tolower STREQUAL "relwithdebinfo")
	message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
	Allowed values are Debug, Release, RelWithDebInfo (case-insensitive).")
endif()

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"        OFF)
option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"  OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int module"         OFF)
option(ENABLE_XCFUN          "Enable XCFun library"          OFF)
option(ENABLE_UTILITIES      "Enable utilities"              OFF)
option(ENABLE_LOCAL_DISKS    "Enable local disks"            OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking" OFF)

add_definitions(-DVAR_MFDS -D_FILE_OFFSET_BITS=64 -DIMPLICIT_NONE)

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
    "${DEFAULT_SCRATCH_BASE_DIR}/DALTON_scratch"
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

if(THIS_CODE STREQUAL "linsca")
    add_definitions(-DINSTALL_BASDIR="${BASIS_DIR}")
endif()

include(ConfigArchitecture)
include(ConfigVersion)
include(BinaryInfo)
include(Sources)
include(ConfigCompilerFlags)
include(ConfigParallelEnvironment)

set(LIBS)

include(ConfigMath)

if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DVAR_INT64)
    if(THIS_CODE STREQUAL "linsca")
        add_definitions(-DVAR_64BITS)
    endif()
endif()

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

# set Gen1Int library and interface
if(ENABLE_GEN1INT)
    message("-- Enable Gen1Int module")
    message("-- Gen1Int library at "${LIB_GEN1INT_PATH})
    add_definitions(-DBUILD_GEN1INT)
    # Gen1int library source code
    include(${LIB_GEN1INT_PATH}/cmake/Sources.cmake)
    # Gen1int library header files
    include_directories(${LIB_GEN1INT_PATH}/include)
    # Gen1Int interface source code
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        ${Gen1Int_LIB_SRCS}
        gen1int/gen1int_matrix.F90
        gen1int/gen1int_shell.F90
        gen1int/gen1int_cube.F90
        gen1int/gen1int_api.F90
        gen1int/gen1int_host.F90)
endif()

# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )
# set ctest own timeout
set(DART_TESTING_TIMEOUT
    "1200"
    CACHE STRING
    "Set timeout in seconds for every single test"
    )
include(Tests)
include(CTest)
enable_testing()

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/mm
    ${PROJECT_SOURCE_DIR}/dft
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

include(GenericMacros)
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

add_executable(
    dalton.x
    ${CMAKE_SOURCE_DIR}/abacus/dalton.F
)
target_link_libraries(
    dalton.x
    dalton
    ${LIBS}
    )

# get size of static allocations
add_custom_target(
    info
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/binary-info/get_static_size.py dalton.x
    )

include(ConfigPackaging)
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR})
