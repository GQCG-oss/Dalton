cmake_minimum_required(VERSION 2.6)

project(DALTON Fortran C)

# do not rebuild if rules (compiler flags) change
# this has advantages and disadvantages
# set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

set(THIS_CODE "trunk")
if(THIS_CODE STREQUAL "linsca")
    add_definitions(-DVAR_LINSCA)
endif()

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
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
# this has advantages and disadvantages
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"        OFF)
option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_EXTERNAL_MATH  "Enable external math library"  ON)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"  OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int module"         OFF)
option(ENABLE_XCFUN          "Enable XCFun library"          OFF)
option(ENABLE_UTILITIES      "Enable utilities"              OFF)
option(ENABLE_LOCAL_DISKS    "Enable local disks"            OFF)

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

include(Archs)
include(Sources)
include(Compilers)

set(LIBS)

include(Math)
if(NOT BLAS_FOUND)
    message("-- Using own BLAS implementation (slow)")
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        ${OWN_BLAS_SOURCES}
        )
endif()
if(NOT LAPACK_FOUND)
    message("-- Using own LAPACK implementation (slow)")
    set(FIXED_FORTRAN_SOURCES
        ${FIXED_FORTRAN_SOURCES}
        ${OWN_LAPACK_SOURCES}
        )
endif()

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_MPI)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif()
endif()

if(ENABLE_OMP)
    add_definitions(-DVAR_OMP)
endif()

if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DVAR_INT64)
    if(THIS_CODE STREQUAL "linsca")
        add_definitions(-DVAR_64BITS)
    endif()
endif()

if(ENABLE_GEN1INT)
    add_definitions(-DBUILD_GEN1INT)
    message(FATAL_ERROR "fix Gen1Int compilation in CMakeLists.txt")
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

configure_file(${CMAKE_SOURCE_DIR}/dalton.in
	${CMAKE_BINARY_DIR}/dalton)
execute_process(COMMAND
	chmod 755 ${CMAKE_BINARY_DIR}/dalton OUTPUT_QUIET)

add_executable(
    dalton.x
    ${C_SOURCES}
    ${FREE_FORTRAN_SOURCES}
    ${FIXED_FORTRAN_SOURCES}
    )

target_link_libraries(
    dalton.x
    ${LIBS}
    ) 
