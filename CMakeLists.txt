

# exit if cmake version is below 2.8
# we need at least 2.8 for external projects support
cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

# set not CYGWIN_WIN32, to quench warning messages from CYGWIN cmake
# will do no harm if not CYGWIN
set(CMAKE_LEGACY_CYGWIN_WIN32 0)

# set project name and languanges that are involved
project(DALTON Fortran C CXX)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# To link with fortran not c++
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake/binary-info
    ${CMAKE_SOURCE_DIR}/cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/mpi
    ${CMAKE_SOURCE_DIR}/cmake/fortran-standard/
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_BINARY_DIR}
    )

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"                           OFF)
option(ENABLE_OPENACC        "Enable OpenACC"                                   OFF)
option(ENABLE_COLLAPSE       "Enable COLLAPSE clause for reorderings"           OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"                              OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"                             OFF)
option(ENABLE_MPI            "Enable MPI parallelization"                       OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    OFF)
option(ENABLE_RELEASE        "Test release mode which undefines MOD_UNRELEASED" OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"                     OFF)
option(ENABLE_STATIC_LINKING "Enable static libraries linking"                  OFF)
option(ENABLE_GEN1INT        "Enable Gen1Int library"                           ON)
option(ENABLE_PELIB          "Enable PE library"                                ON)
option(ALWAYS_RESET_EXTERNAL "Always remove builds stamps of external projects" OFF)
option(ENABLE_CSR            "Enable MKL compressed sparse row"                 OFF)
option(ENABLE_SCALAPACK      "Enable SCALAPACK"                                 OFF)
option(ENABLE_SCALASCA       "Enable scalasca profiler mode"                    OFF)
option(ENABLE_TIMINGS        "Enable TIMINGS"                                   OFF)
option(ENABLE_XCFUN          "Enable XCFUN"                                     OFF)
option(ENABLE_INTEREST       "Enable interest integrals"                        OFF)
option(ENABLE_CRAY_WRAPPERS  "Enable cray wrappers for BLAS/LAPACK and MPI"     OFF)
option(ENABLE_LSEEK          "Enable lseek"                                     OFF)
option(ENABLE_DEBUGPBC       "Enable DEBUG_PBC"                                 OFF)
option(ENABLE_RSP            "Enable Response code in LSDALTON"                 ON)
option(ENABLE_ICHOR          "Enable ichor integrals"                           OFF)
option(ENABLE_OPENRSP        "Enable OpenRSP"                                   OFF)

option(ENABLE_BUILTIN_BLAS   "Enable builtin BLAS implementation (slow)"        OFF)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_BUILTIN_LAPACK "Enable builtin LAPACK implementation (slow)"      OFF)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)

set(EXTERNAL_LIBS)

include(ConfigVersion)
include(SourcesDALTON)
include(SourcesLSDALTON)
include(ConfigArchitecture)
include(ConfigCompilerFlags)
include(CheckFortranStandardCompatibility)
include(ConfigExternal)
include(ConfigOMP)

# math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

include(ConfigExplicitLibs)
include(ConfigMPI)
include(ConfigSafeGuards)
include(ConfigTesting)
include(GenericMacros)
include(BinaryInfo)
include(mergestaticlibs)

# set code coverage
if(ENABLE_CODE_COVERAGE)
    set(EXTERNAL_LIBS ${EXTERNAL_LIBS} gcov)
endif()

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

include_directories(
    ${CMAKE_SOURCE_DIR}/DALTON/include
    ${CMAKE_SOURCE_DIR}/DALTON/dft
    ${CMAKE_SOURCE_DIR}/DALTON/abacus
    ${CMAKE_SOURCE_DIR}/LSDALTON/dft
    ${CMAKE_SOURCE_DIR}/LSDALTON/geomopt
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    ${PROJECT_BINARY_DIR}
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

# set definitions
include(Definitions)
if(ENABLE_OPENRSP)
    include(DefinitionsOpenRSP)
endif()

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
configure_script(
    ${CMAKE_SOURCE_DIR}/DALTON/dalton_config.in
    ${CMAKE_BINARY_DIR}/dalton_config.h
    )

# if BLAS and/or LAPACK not found, add own sources to the list of
# sources to compile
if(USE_BUILTIN_BLAS)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${DALTON_OWN_BLAS_SOURCES}
        )
    set(LSDALTON_FIXED_FORTRAN_SOURCES
        ${LSDALTON_FIXED_FORTRAN_SOURCES}
        ${LSDALTON_OWN_BLAS_SOURCES}
        )
endif()
if(USE_BUILTIN_LAPACK)
    set(DALTON_FIXED_FORTRAN_SOURCES
        ${DALTON_FIXED_FORTRAN_SOURCES}
        ${DALTON_OWN_LAPACK_SOURCES}
        )
    set(LSDALTON_FIXED_FORTRAN_SOURCES
        ${LSDALTON_FIXED_FORTRAN_SOURCES}
        ${LSDALTON_OWN_LAPACK_SOURCES}
        )
endif()

include(LibsDALTON)
include(LibsLSDALTON)

# copy basis/ to build directory
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/basis ${CMAKE_BINARY_DIR})

# this controlls "make install" target
include(ConfigMakeInstall)

# this controlls "make release" target
include(ConfigMakeRelease)

# give information about system, compiler flags, and size of static allocations
set(STATIC_MEM_INFO_BINARIES dalton lsdalton)
include(ConfigInfo)
