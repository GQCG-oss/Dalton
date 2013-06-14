

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
    ${CMAKE_SOURCE_DIR}/../cmake/binary-info
    ${CMAKE_SOURCE_DIR}/../cmake/compilers
    ${CMAKE_SOURCE_DIR}/../cmake/math
    ${CMAKE_SOURCE_DIR}/../cmake/mpi
    ${CMAKE_SOURCE_DIR}/../cmake
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
include(GenericMacros)

# set definitions
include(Definitions)

# forward CPP directly to the code
set(CPP)
if(NOT "${CPP}" STREQUAL "")
    add_definitions(${CPP})
endif()

set(PARENT_DEFINITIONS "-DBUILD_GEN1INT -DVAR_LSDALTON")
set(ExternalProjectCMakeArgs
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/include
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    )

# if BLAS and/or LAPACK not found, add own sources to the list of
# sources to compile
if(USE_OWN_BLAS)
    set(LSDALTON_FIXED_FORTRAN_SOURCES
        ${LSDALTON_FIXED_FORTRAN_SOURCES}
        ${LSDALTON_OWN_BLAS_SOURCES}
        )
endif()
if(USE_OWN_LAPACK)
    set(LSDALTON_FIXED_FORTRAN_SOURCES
        ${LSDALTON_FIXED_FORTRAN_SOURCES}
        ${LSDALTON_OWN_LAPACK_SOURCES}
        )
endif()

if(ENABLE_SCALASCA)
    set(SCALASCA_INSTRUMENT ${CMAKE_Fortran_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/scalasca.in
        ${CMAKE_SOURCE_DIR}/scalascaf90.sh
        )
    set(SCALASCA_INSTRUMENT ${CMAKE_C_COMPILER})
    configure_script(
        ${CMAKE_SOURCE_DIR}/scalasca.in
        ${CMAKE_SOURCE_DIR}/scalascaCC.sh
        )
    unset(SCALASCA_INSTRUMENT)
    SET(CMAKE_Fortran_COMPILER "../scalascaf90.sh")
    SET(CMAKE_C_COMPILER "../scalascaCC.sh")
    MESSAGE(STATUS "Fortran Compiler " ${CMAKE_Fortran_COMPILER})
    MESSAGE(STATUS "C Compiler       " ${CMAKE_C_COMPILER})
endif()

# set code coverage
if(ENABLE_CODE_COVERAGE)
  set(LIBS ${LIBS} gcov)
endif()

# set cdash buildname
set(BUILDNAME
    "BUILDNAME-not-set"
    CACHE STRING
    "Name of build on the dashboard"
    )
# set ctest own timeout
if(cmake_build_type_tolower STREQUAL "debug")
    set(DART_TESTING_TIMEOUT
        "1100"
        CACHE STRING
        "Set timeout in seconds for every single test"
        )
else()
    set(DART_TESTING_TIMEOUT
        "400"
        CACHE STRING
        "Set timeout in seconds for every single test"
        )
endif()

include(BinaryInfo)

include(TestsLSDALTON)
include(CTest)
include(mergestaticlibs)
enable_testing()

set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)

get_directory_property( DirDefs DIRECTORY ${CMAKE_SOURCE_DIR} COMPILE_DEFINITIONS )
foreach( d ${DirDefs} )
    message( STATUS "Found Define: " ${d} )
endforeach()
message("-- COMPILE_DEFINITIONS   : ${DirDefs}" )

include_directories(
    ${PROJECT_SOURCE_DIR}/mm
    ${PROJECT_SOURCE_DIR}/dft
    ${CMAKE_Fortran_MODULE_DIRECTORY}
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

configure_script(
    ${CMAKE_SOURCE_DIR}/lsdalton.in
    ${CMAKE_BINARY_DIR}/lsdalton
    )

execute_process(COMMAND echo LSDALTON_script="${CMAKE_BINARY_DIR}/lsdalton" OUTPUT_FILE ${CMAKE_BINARY_DIR}/test.scipt)

# tell CMake the source won't be available until build time
SET_SOURCE_FILES_PROPERTIES(${PROJECT_BINARY_DIR}/binary_info.f90 PROPERTIES GENERATED 1)

add_library(
    lsutillib_precision
    ${LSUTIL_PRECISION_SOURCES}
    )

add_library(
    matrixmlib
    ${LSUTIL_MATRIXM_SOURCES}
    )
target_link_libraries(matrixmlib lsutillib_precision)

add_library(
    lsutillib_common
    ${LSUTIL_COMMON_SOURCES}
    )
target_link_libraries(lsutillib_common matrixmlib)

add_library(
    matrixolib
    ${LSUTIL_MATRIXO_SOURCES}
    ${LSUTIL_MATRIXO_C_SOURCES}
    )
target_link_libraries(matrixolib lsutillib_common)

add_library(
    matrixulib
    ${LSUTIL_MATRIXU_SOURCES}
    )

target_link_libraries(matrixulib matrixolib)

add_library(
    pdpacklib
    ${LSDALTON_FIXED_FORTRAN_SOURCES}
    )

target_link_libraries(pdpacklib matrixulib)

#automatially generate the manual_reorderdings.F90
execute_process(COMMAND python ${CMAKE_SOURCE_DIR}/lsutil/autogen/generate_man_reord.py nocollapse ${DirDefs})

add_library(
    lsutiltypelib_common
    ${LSUTIL_TYPE_SOURCES}
    )

target_link_libraries(lsutiltypelib_common pdpacklib)

add_library(
    lsutillib
    ${LSUTILLIB_SOURCES}
    ${PROJECT_BINARY_DIR}/binary_info.F90
    )
target_link_libraries(lsutillib lsutiltypelib_common)
# Make sure the source is generated before the lsutillib lib is built
add_dependencies(lsutillib generate_revinfo)
add_dependencies(lsutillib generate_compinfo)

set(PARENT_DEFINITIONS "-DVAR_LSDALTON")

set(ExternalProjectCMakeArgs
    -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/external
    -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
    -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
    -DPARENT_INCLUDE_DIR=${PROJECT_SOURCE_DIR}/include
    -DPARENT_MODULE_DIR=${PROJECT_BINARY_DIR}/modules
    -DPARENT_DEFINITIONS=${PARENT_DEFINITIONS}
    )

add_library(
    xcfun_interface
    xcfun_host/xcfun_host.F90
    )
add_dependencies(xcfun_interface lsutillib_precision)

if(ENABLE_XCFUN)
    add_external(xcfun)
    add_dependencies(xcfun_interface xcfun)
    add_definitions(-DVAR_XCFUN)
    set(LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun_f90_bindings.a
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        ${LIBS}
        )
endif()

if(NOT ENABLE_RELEASE)
if(ENABLE_INTEREST)
        add_definitions(-DVAR_INTEREST)
	add_library(
		interestlib
    		${INTERESTLIB_SOURCES}
		)
                if(ENABLE_XCFUN)
		target_link_libraries(interestlib xcfun_interface)
                endif()
endif()
endif()

add_library(
    fmmlib
    ${FMM_SOURCES}
    ${FMM_C_SOURCES}
    )
        if(ENABLE_XCFUN)
	target_link_libraries(fmmlib xcfun_interface)
        endif()
add_dependencies(fmmlib lsutillib_precision)
add_dependencies(fmmlib lsutillib_common)
add_dependencies(fmmlib lsutiltypelib_common)

if(NOT ENABLE_RELEASE)
if(ENABLE_INTEREST)
target_link_libraries(fmmlib interestlib)
endif()
endif()

add_library(
    dftfunclib
    ${DFTFUNC_SOURCES}
    ${DFTFUNC_F_SOURCES}
    )
target_link_libraries(dftfunclib fmmlib)

add_library(
    lsintlib
    ${LSINT_SOURCES}
    )
target_link_libraries(lsintlib dftfunclib)
add_dependencies(lsintlib pdpacklib)
add_dependencies(lsintlib lsutillib)
add_dependencies(lsintlib xcfun_interface)
add_library(
    pbclib
    ${PBC_FORTRAN_SOURCES}
    )

target_link_libraries(pbclib lsintlib)

add_library(
    ddynamlib
    ${DDYNAM_SOURCES}
    )

target_link_libraries(ddynamlib lsintlib)

add_library(
    declib
    ${DEC_SOURCES}
    ${DEC_C_SOURCES}
    )

target_link_libraries(declib lsintlib)

add_library(
    solverutillib
    ${SOLVERUTIL_SOURCES}
    )

target_link_libraries(solverutillib lsintlib)

add_library(
    rspsolverlib
    ${RSPSOLVER_SOURCES}
    )

target_link_libraries(rspsolverlib solverutillib)

add_library(
    linearslib
    ${LINEARS_SOURCES}
    )

target_link_libraries(linearslib rspsolverlib)

add_library(
    rsp_propertieslib
    ${RSP_PROPERTIES_SOURCES}
    )

target_link_libraries(rsp_propertieslib lsintlib)


add_library(
    geooptlib
    ${GEOOPT_SOURCES}
    )

target_link_libraries(geooptlib lsintlib)

add_library(
    lsdaltonmain 
    ${LSDALTONMAIN_FORTRAN_SOURCES}
    )

target_link_libraries(lsdaltonmain pbclib)
target_link_libraries(lsdaltonmain geooptlib)
target_link_libraries(lsdaltonmain linearslib)
target_link_libraries(lsdaltonmain declib)
target_link_libraries(lsdaltonmain ddynamlib)
target_link_libraries(lsdaltonmain rsp_propertieslib)
target_link_libraries(lsdaltonmain rspsolverlib)
target_link_libraries(lsdaltonmain xcfun_interface)

add_executable(
    lsdalton.x
    ${CMAKE_SOURCE_DIR}/lsdaltonsrc/lsdalton_wrapper.f90
    ${LINK_FLAGS}
    )

add_executable(
    lslib_tester.x
    ${LSLIB_SOURCES}
    ${LINK_FLAGS}
    )

if(MPI_FOUND)
    # Simen's magic fix for Mac/GNU/OpenMPI
    if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
        if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            SET_TARGET_PROPERTIES(lsdalton.x     PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
            SET_TARGET_PROPERTIES(lslib_tester.x PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
        endif()
    endif()
endif()

  MERGE_STATIC_LIBS(
    rsp_prop
    rsp_propertieslib
    )

if(NOT ENABLE_RELEASE)
  if(ENABLE_INTEREST)
    MERGE_STATIC_LIBS(
      lsint
      lsintlib
      interestlib
      )
  else()
    MERGE_STATIC_LIBS(
      lsint
      lsintlib
      )
  endif()
else()
  MERGE_STATIC_LIBS(
    lsint
    lsintlib
    )
endif()


MERGE_STATIC_LIBS(
    lsdalton
    lsutillib_precision
    matrixmlib
    lsutillib_common
    matrixolib
    matrixulib
    pdpacklib
    lsutiltypelib_common
    lsutillib
    fmmlib
    dftfunclib
    lsint
    pbclib
    ddynamlib
    declib
    solverutillib
    rspsolverlib
    linearslib
    rsp_prop
    geooptlib
    xcfun_interface
    lsdaltonmain 
    )

target_link_libraries(
    lsdalton
    lsdaltonmain
    ${LIBS}
    )

target_link_libraries(
   lsdalton.x
   lsdalton
   ) 

target_link_libraries(
   lslib_tester.x
   lsdalton
   ) 

# copy basis/ to build directory
execute_process(COMMAND cp -r ${CMAKE_SOURCE_DIR}/../basis ${CMAKE_BINARY_DIR})

# this controlls "make install" target
# include(ConfigMakeInstall)

# this controlls "make release" target
# include(ConfigMakeRelease)

# give information about system, compiler flags, and size of static allocations
set(STATIC_MEM_INFO_BINARIES lsdalton)
include(ConfigInfo)
