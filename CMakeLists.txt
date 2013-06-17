

# exit if cmake version is below 2.8
# we need at least 2.8 for external projects support
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# set project name and languanges that are involved
project(DALTON Fortran C CXX)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake/binary-info
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/mpi
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_BINARY_DIR}
    )

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"                           OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                              OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"                             OFF)
option(ENABLE_MPI            "Enable MPI parallelization"                       OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    OFF)
option(ENABLE_RELEASE        "Test release mode which undefines MOD_UNRELEASED" OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"                     OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking"                  OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int library"                           ON)
option(ALWAYS_RESET_EXTERNAL "Always remove builds stamps of external projects" ON)
option(ENABLE_CSR            "Enable MKL compressed sparse row"                 OFF)
option(ENABLE_SCALAPACK      "Enable SCALAPACK"                                 OFF)
option(ENABLE_SCALASCA       "Enable scalasca profiler mode"                    OFF)
option(ENABLE_TIMINGS        "Enable TIMINGS"                                   OFF)
option(ENABLE_XCFUN          "Enable XCFUN"                                     OFF)
option(ENABLE_INTEREST       "Enable interest integrals"                        OFF)
option(ENABLE_CRAY_WRAPPERS  "Enable cray wrappers for BLAS/LAPACK and MPI"     OFF)
option(ENABLE_LSEEK          "Enable lseek"                                     OFF)
option(ENABLE_DEBUGPBC       "Enable DEBUG_PBC"                                 OFF)
option(ENABLE_MPI2_DETECTION "Enable detection of MPI-2 standard"               OFF)

set(LIBS)

include(ConfigVersion)
include(SourcesDALTON)
include(SourcesLSDALTON)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(ConfigExternal)
include(ConfigOMP)
include(ConfigMath)
include(ConfigMPI)
include(ConfigSafeGuards)
include(ConfigTesting)
include(GenericMacros)
include(BinaryInfo)
include(mergestaticlibs)

# set code coverage
if(ENABLE_CODE_COVERAGE)
    set(LIBS ${LIBS} gcov)
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/DALTON/include
    ${PROJECT_SOURCE_DIR}/DALTON/dft
    ${PROJECT_SOURCE_DIR}/LSDALTON/dft
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

# set definitions
include(Definitions)

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

# configure INSTALL_WRKMEM and INSTALL_MMWORK
include(ConfigWORKMEM)

configure_script(
    ${CMAKE_SOURCE_DIR}/DALTON/dalton.in
    ${CMAKE_BINARY_DIR}/dalton
    )
configure_script(
    ${CMAKE_SOURCE_DIR}/LSDALTON/lsdalton.in
    ${CMAKE_BINARY_DIR}/lsdalton
    )

# if BLAS and/or LAPACK not found, add own sources to the list of
# sources to compile
if(USE_OWN_BLAS)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${DALTON_OWN_BLAS_SOURCES}
        )
    set(LSDALTON_FIXED_FORTRAN_SOURCES
        ${LSDALTON_FIXED_FORTRAN_SOURCES}
        ${LSDALTON_OWN_BLAS_SOURCES}
        )
endif()
if(USE_OWN_LAPACK)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${DALTON_OWN_LAPACK_SOURCES}
        )
    set(LSDALTON_FIXED_FORTRAN_SOURCES
        ${LSDALTON_FIXED_FORTRAN_SOURCES}
        ${LSDALTON_OWN_LAPACK_SOURCES}
        )
endif()

# add binary_info.F90
set(GENERATED_FILES
    ${CMAKE_BINARY_DIR}/binary_info.F90
    )

# always invalidate generated files
add_custom_target(generated_files
    ALL
    DEPENDS
    ${GENERATED_FILES}
)

include(LibsDALTON)
include(LibsLSDALTON)

# copy basis/ to build directory
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR})

# this controlls "make install" target
include(ConfigMakeInstall)

# this controlls "make release" target
include(ConfigMakeRelease)

# give information about system, compiler flags, and size of static allocations
set(STATIC_MEM_INFO_BINARIES dalton)
include(ConfigInfo)
