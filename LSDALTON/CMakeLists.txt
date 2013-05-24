cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(LSDALTON Fortran C CXX)

# do not rebuild if rules (compiler flags) change
set(CMAKE_SKIP_RULE_DEPENDENCY TRUE)

# flag to add a fix OpenMPI bug for gfortran/macOSX
set(CMAKE_GFORTRAN_OPENMPI_MAGIC_FIX FALSE)

add_definitions(-DTHIS_IS_CMAKE_BUILD)

# these are paths that CMake will search for cmake module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/../cmake
    ${CMAKE_SOURCE_DIR}/../cmake/compilers
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/../cmake/binary-info
    )

message("-- System: ${CMAKE_SYSTEM_NAME}")

include(ConfigVersion)

include(ConfigSafeGuards)
if(cmake_build_type_tolower STREQUAL "debug")
	add_definitions(-DVAR_DEBUGINT)
	add_definitions(-DVAR_DEBUG)
endif()

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"        OFF)
option(ENABLE_CSR            "Enable MKL compressed sparse row" OFF)
option(ENABLE_SCALAPACK      "Enable SCALAPACK"              OFF)
option(ENABLE_TIMINGS        "Enable TIMINGS"                OFF)
option(ENABLE_PROFILING      "Enable profiling"              OFF)
option(ENABLE_BOUNDS_CHECK   "Enable bounds check"           OFF)
option(ENABLE_CODE_COVERAGE  "Enable code coverage"          OFF)
option(ENABLE_INTERNAL_MATH  "Enable internal math library implementation"  OFF)
option(ENABLE_XCFUN          "Enable XCFUN"                  OFF)
option(ENABLE_MPI            "Enable MPI parallelization"    OFF)
option(ENABLE_MPI32          "Enable 32 bit integer MPI lib workaround" OFF)
option(ENABLE_OMP            "Enable OpenMP parallelization" OFF)
option(ENABLE_INTEREST       "Enable interest integrals"     OFF)
option(ENABLE_CUDA           "Enable CUDA GPU acceleration"  OFF)
option(ENABLE_OPENRSP        "Enable OpenRSP module"         OFF)
option(ENABLE_UTILITIES      "Enable utilities"              OFF)
option(ENABLE_LOCAL_DISKS    "Enable local disks"            OFF)
option(ENABLE_CRAYWRAPPER    "Enable cray wrapper"           OFF)
option(ENABLE_LSEEK          "Enable lseek"                  OFF)
option(ENABLE_DEBUGPBC       "Enable DEBUG_PBC"              OFF)
option(ENABLE_RELEASE        "Enable Release"                OFF)
#option(ENABLE_EXTERNLIB      "Create lsdaltonextern.a extern lib " OFF)
option(ENABLE_GIT_SUBMODULES "Download external projects using Git submodules at compilation time" ON)


# find a good default scratch base directory
foreach(dir "/work" "/scratch" "/scr" "/temp" "/tmp")
   if(EXISTS "${dir}/")
       set(DEFAULT_SCRATCH_BASE_DIR ${dir})
       break()
   endif()
endforeach()

# set scratch directory
set(SCRATCH_DIR
    "${DEFAULT_SCRATCH_BASE_DIR}/LSDALTON_scratch"
    CACHE STRING
    "Default scratch directory"
    )
message("-- SCRATCH_DIR set to: ${SCRATCH_DIR}")

set(BASIS_DIR
    "${CMAKE_SOURCE_DIR}/../basis"
    CACHE STRING
    "Basis set directory"
    )
message("-- BASIS_DIR set to: ${BASIS_DIR}")

    add_definitions(-DINSTALL_BASDIR="${BASIS_DIR}")

include(ConfigArchitecture)
include(Sources)
include(ConfigCompilerFlags)
include(ConfigExternal)

set(LIBS)

if(ENABLE_OPENRSP)
	set(PARENT_DEFINITIONS "-DBUILD_GEN1INT -DVAR_LSDALTON")
	add_definitions(-DBUILD_GEN1INT)
	set(CMAKE_Fortran_MODULE_DIRECTORY
	   ${PROJECT_BINARY_DIR}/modules
	   )
        set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    	   ${PROJECT_BINARY_DIR}/lib
    	   )
endif()

set(PARENT_DEFINITIONS "-DBUILD_GEN1INT -DVAR_LSDALTON")
add_definitions(-DVAR_LSDALTON)

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

if(NOT ENABLE_CRAYWRAPPER)
    include(Math)
    if(NOT BLAS_FOUND)
       message("-- Using own BLAS implementation (slow)")
       set(FIXED_FORTRAN_SOURCES
          ${FIXED_FORTRAN_SOURCES}
          ${OWN_BLAS_SOURCES}
          )
    endif()
    if(NOT LAPACK_FOUND)
       set(FIXED_FORTRAN_SOURCES
         ${FIXED_FORTRAN_SOURCES}
         ${OWN_LAPACK_SOURCES}
       )
    endif()
endif()

if(ENABLE_RELEASE)

else()
   add_definitions(-DMOD_UNRELEASED)
endif()

if(ENABLE_MPI)
    find_package(MPI)
    if(MPI_FOUND)
        add_definitions(-DVAR_LSMPI)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
        if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
          if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
            set(CMAKE_GFORTRAN_OPENMPI_MAGIC_FIX TRUE)
          endif()
        endif()
    endif()
endif()

if(ENABLE_MPI32)
    if(ENABLE_64BIT_INTEGERS)
        add_definitions(-DVAR_LSMPI_32)
    else()
        message("--mpi32 should only be enabled if --int64 is used.")
    endif()
endif()

if(ENABLE_OMP)
    add_definitions(-DVAR_OMP)
endif()

if(ENABLE_LSEEK)
    add_definitions(-DHAVE_NO_LSEEK64)
endif()

if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DVAR_INT64)
    add_definitions(-DVAR_64BITS)
endif()

if(ENABLE_CSR)
    add_definitions(-DVAR_MKL)
endif()

if(ENABLE_SCALAPACK)
    add_definitions(-DVAR_SCALAPACK)
endif()

if(ENABLE_TIMINGS)
    add_definitions(-DVAR_TIME)
endif()

if(ENABLE_DEBUGPBC)
  add_definitions(-DDEBUGPBC)
endif()
if(PROFILER STREQUAL "SCALASCA")
  SET(CMAKE_Fortran_COMPILER "../scalascaf90.sh")
  SET(CMAKE_C_COMPILER "../scalascaCC.sh")
  MESSAGE(STATUS "Fortran Compiler " ${CMAKE_Fortran_COMPILER})
  MESSAGE(STATUS "C Compiler       " ${CMAKE_C_COMPILER})
endif()

# set code coverage
if(ENABLE_CODE_COVERAGE)
  set(LIBS ${LIBS} gcov)
endif()

# gen1int library
if(ENABLE_OPENRSP)
   set(LIBS ${LIBS} ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a)
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

include(Tests)
include(CTest)
include(merge)
enable_testing()

if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
# Hack to make gfortran 4.6.2 work, also done for other gfortran versions
  if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    set(CMAKE_Fortran_MODULE_DIRECTORY
      ${PROJECT_BINARY_DIR}/modules
      )
  else()
    set(CMAKE_Fortran_MODULE_DIRECTORY
      ${PROJECT_BINARY_DIR}/modules
      )
  endif()
else()
    set(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/modules)
#    if(${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
#    endif()
endif()

if(cmake_build_type_tolower STREQUAL "debug")
       message("-- Fortran Compiler Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_DEBUG}")
       message("-- C Compiler Flags      : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_DEBUG}")
else()
       message("-- Fortran Compiler Flags: ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_RELEASE}")
       message("-- C Compiler Flags      : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_RELEASE}")
endif()

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

configure_file(
    ${CMAKE_SOURCE_DIR}/lsdalton.in
    ${CMAKE_BINARY_DIR}/lsdalton.in
    @ONLY
    )
execute_process(COMMAND	cat ${CMAKE_BINARY_DIR}/lsdalton.in OUTPUT_FILE ${CMAKE_BINARY_DIR}/lsdalton)
execute_process(COMMAND	chmod 755 ${CMAKE_BINARY_DIR}/lsdalton OUTPUT_QUIET)

execute_process(COMMAND echo LSDALTON_script="${CMAKE_BINARY_DIR}/lsdalton" OUTPUT_FILE ${CMAKE_BINARY_DIR}/test.scipt)

if(ENABLE_OPENRSP)
	add_subdirectory(gen1int)
endif()

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
    ${FIXED_FORTRAN_SOURCES}
    )

target_link_libraries(pdpacklib matrixulib)

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
#    add_subdirectory(xcfun_host)
    add_dependencies(xcfun_interface xcfun)
    add_definitions(-DVAR_XCFUN)
    set(LIBS
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun_f90_bindings.a
        ${PROJECT_BINARY_DIR}/external/lib/libxcfun.a
        ${LIBS}
        )
endif()

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

if(ENABLE_INTEREST)
target_link_libraries(fmmlib interestlib)
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

if(ENABLE_OPENRSP)
	add_dependencies(gen1int_interface lsutillib_precision)
	add_dependencies(gen1int_interface matrixmlib)
	add_dependencies(gen1int_interface lsutillib_common)
	add_dependencies(gen1int_interface matrixolib)
	add_dependencies(gen1int_interface matrixulib)
	add_dependencies(gen1int_interface lsutiltypelib_common)
	add_dependencies(lsintlib gen1int)
	add_dependencies(gen1int_interface lsintlib)
	add_dependencies(linearslib gen1int_interface)
	add_dependencies(lsdaltonmain gen1int_interface)
	add_dependencies(rsp_propertieslib gen1int_interface)
endif()


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

if(CMAKE_GFORTRAN_OPENMPI_MAGIC_FIX)
  SET_TARGET_PROPERTIES(lsdalton.x PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
  SET_TARGET_PROPERTIES(lslib_tester.x PROPERTIES LINK_FLAGS "-Wl,-commons,use_dylibs")
endif()

if(ENABLE_OPENRSP)
  MERGE_STATIC_LIBS(
    rsp_prop
    rsp_propertieslib
    gen1int_interface
#   ${PROJECT_BINARY_DIR}/external/lib/libgen1int.a
    )
else()
  MERGE_STATIC_LIBS(
    rsp_prop
    rsp_propertieslib
    )
endif()

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

include(Install)
