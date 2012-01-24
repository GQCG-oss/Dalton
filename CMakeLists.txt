cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

project(DALTON Fortran C)

# set version
set(VERSION_MAJOR 2011)
set(VERSION_MINOR 0)
set(VERSION_PATCH 0)

set(THIS_CODE "trunk")
if(THIS_CODE STREQUAL "linsca")
    add_definitions(-DVAR_LINSCA)
endif()

add_definitions(-DTHIS_IS_CMAKE_BUILD)

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
option(ENABLE_INTERNAL_MATH  "Enable internal math library implementation"  OFF)
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

include(Archs)
include(Sources)
include(FCompilers)
include(CCompilers)

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
    set(Gen1Int_SRCS ${LIB_GEN1INT_PATH}/src/error_stop.F90
                     ${LIB_GEN1INT_PATH}/src/xtimer.F90
                     ${LIB_GEN1INT_PATH}/src/norm_contr_cgto.F90
                     ${LIB_GEN1INT_PATH}/src/norm_contr_sgto.F90
                     ${LIB_GEN1INT_PATH}/src/reorder_ints.F90
                     ${LIB_GEN1INT_PATH}/src/trace_ints.F90
                     ${LIB_GEN1INT_PATH}/src/dump_info.F90
                     ${LIB_GEN1INT_PATH}/src/binom_coeff.F90
                     ${LIB_GEN1INT_PATH}/src/next_permutation.F90
                     ${LIB_GEN1INT_PATH}/src/geom_total.F90
                     ${LIB_GEN1INT_PATH}/src/sort_cents.F90
                     ${LIB_GEN1INT_PATH}/src/shell_scatter.F90
                     ${LIB_GEN1INT_PATH}/src/geom_part_zero.F90
                     ${LIB_GEN1INT_PATH}/src/geom_part_one.F90
                     ${LIB_GEN1INT_PATH}/src/hgto_to_cgto.F90
                     ${LIB_GEN1INT_PATH}/src/hgto_to_lcgto.F90
                     ${LIB_GEN1INT_PATH}/src/hgto_to_sgto.F90
                     ${LIB_GEN1INT_PATH}/src/london_mom_hgto.F90
                     ${LIB_GEN1INT_PATH}/src/const_contr_ints.F90
                     ${LIB_GEN1INT_PATH}/src/carmom_deriv.F90
                     ${LIB_GEN1INT_PATH}/src/carmom_init.F90
                     ${LIB_GEN1INT_PATH}/src/carmom_to_hbra.F90
                     ${LIB_GEN1INT_PATH}/src/carmom_hrr_ket.F90
                     ${LIB_GEN1INT_PATH}/src/prim_hgto_carmom.F90
                     ${LIB_GEN1INT_PATH}/src/contr_cgto_carmom.F90
                     ${LIB_GEN1INT_PATH}/src/contr_sgto_carmom.F90
                     ${LIB_GEN1INT_PATH}/src/contr_csgto_carmom.F90
                     ${LIB_GEN1INT_PATH}/src/aux_boys_vec.F90
                     ${LIB_GEN1INT_PATH}/src/nucpot_nderv.F90
                     ${LIB_GEN1INT_PATH}/src/nucpot_hket.F90
                     ${LIB_GEN1INT_PATH}/src/nucpot_hbra.F90
                     ${LIB_GEN1INT_PATH}/src/nucpot_to_carmom.F90
                     ${LIB_GEN1INT_PATH}/src/nucpot_to_hgto.F90
                     ${LIB_GEN1INT_PATH}/src/prim_hgto_nucpot.F90
                     ${LIB_GEN1INT_PATH}/src/contr_cgto_nucpot.F90
                     ${LIB_GEN1INT_PATH}/src/contr_sgto_nucpot.F90
                     ${LIB_GEN1INT_PATH}/src/gaupot_nderv.F90
                     ${LIB_GEN1INT_PATH}/src/prim_hgto_gaupot.F90
                     ${LIB_GEN1INT_PATH}/src/contr_cgto_gaupot.F90
                     ${LIB_GEN1INT_PATH}/src/contr_sgto_gaupot.F90
                     ${LIB_GEN1INT_PATH}/src/gen1int.F90)
    # Gen1Int interface source code
    set(FREE_FORTRAN_SOURCES
        ${FREE_FORTRAN_SOURCES}
        ${Gen1Int_SRCS}
        gen1int/gen1int_matrix.F90
        gen1int/gen1int_shell.F90
        gen1int/gen1int_interface.F90)
    #find_library(GEN1INT_LIBRARIES NAMES libgen1int.a
    #             PATHS ${LIB_GEN1INT_PATH}
    #             PATH_SUFFIXES build lib obj src)
    #if(GEN1INT_LIBRARIES)
    #    message("-- Find Gen1Int library at "${GEN1INT_LIBRARIES})
    #    message("-- Enable Gen1Int module")
    #    add_definitions(-DBUILD_GEN1INT)
    #    set(LIBS
    #        ${LIBS}
    #        ${GEN1INT_LIBRARIES})
    #    set(FREE_FORTRAN_SOURCES
    #        ${FREE_FORTRAN_SOURCES}
    #        gen1int/gen1int_matrix.F90
    #        gen1int/gen1int_shell.F90
    #        gen1int/gen1int_interface.F90)
    #else()
    #    message("-- Disable Gen1Int module")
    #endif()
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

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/dalton.setup.in
    ${CMAKE_BINARY_DIR}/dalton.setup
    @ONLY
    )
execute_process(COMMAND	cat ${CMAKE_BINARY_DIR}/dalton.setup ${CMAKE_SOURCE_DIR}/dalton.gnr OUTPUT_FILE ${CMAKE_BINARY_DIR}/dalton)
execute_process(COMMAND	chmod 755 ${CMAKE_BINARY_DIR}/dalton OUTPUT_QUIET)

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

# get size of static allocations in dalton.x
add_custom_target(
    info
    COMMAND ${CMAKE_SOURCE_DIR}/cmake/static_size.py dalton.x
    )

include(CompilationInfo)
include(Install)
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR})
