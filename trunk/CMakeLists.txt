cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH
    ${CMAKE_MODULE_PATH}
    ${CMAKE_SOURCE_DIR}/cmake
    )

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif ()

project(Dalton Fortran C)

option(DALTON_ENABLE_PROFILING     "Enable profiling"          OFF)
option(DALTON_ENABLE_BOUNDS_CHECK  "Enable bounds check"       OFF)
option(DALTON_ENABLE_CODE_COVERAGE "Enable code coverage"      OFF)
option(DALTON_ENABLE_EXTERNAL_MATH "Use external math library" ON)
option(DALTON_ENABLE_MPI           "Use MPI parallelization"   OFF)
option(DALTON_ENABLE_XCFUN         "Use XCFun library"         OFF)

set(DALTON_WORK_MEM_WORDS
    "100000000"
    CACHE STRING
    "DALTON work memory in words"
    )
add_definitions(-D'INSTALL_WRKMEM=${DALTON_WORK_MEM_WORDS}' -D'INSTALL_MMWORK=1' -D_FILE_OFFSET_BITS=64 -DIMPLICIT_NONE)

# now only Linux, later add more machine info
add_definitions(-DSYS_LINUX)

include(ConfigCompilerFlags)
if (CMAKE_COMPILER_IS_GNUC)
    set (CMAKE_C_FLAGS "-Wall")
    set (CMAKE_C_FLAGS_DEBUG "-O0 -g3")
    set (CMAKE_C_FLAGS_RELEASE "-g -O2 -Wno-unused")
elseif (CMAKE_C_COMPILER_ID MATCHES Intel)
    set (CMAKE_C_FLAGS "-wd981 -wd279 -wd383 -vec-report0 -wd1572 -wd177")
    set (CMAKE_C_FLAGS_DEBUG "-g -O0")
    set (CMAKE_C_FLAGS_RELEASE "-g -O2")
    set (CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -shared-intel")
endif ()
config_compiler_flags(C)

if (CMAKE_COMPILER_IS_GNUFORTRAN)
#   set(CMAKE_Fortran_FLAGS         "-Wall")
    set(CMAKE_Fortran_FLAGS         "     ")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g3")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -funroll-all-loops")
    add_definitions(-DVAR_GFORTRAN)
elseif (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    set(CMAKE_Fortran_FLAGS         "-w -assume byterecl")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    add_definitions(-DVAR_IFC)
endif ()
config_compiler_flags(Fortran)

set(DALTON_LIBRARIES_TO_LINK)

if (DALTON_ENABLE_EXTERNAL_MATH)
    find_package(BLAS)
    find_package(LAPACK)
    set(DALTON_LIBRARIES_TO_LINK
        ${DALTON_LIBRARIES_TO_LINK}
        blas
        lapack
        )
endif ()

if (DALTON_ENABLE_MPI)
    find_package(MPI)
    add_definitions(-DVAR_MPI)
    if (MPI_FOUND)
        set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${MPI_COMPILE_FLAGS}")
        include_directories(${MPI_INCLUDE_PATH})
    endif ()
endif ()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_BINARY_DIR}/modules
    )

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY
    ${PROJECT_BINARY_DIR}/lib
    )

set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/modules
    )

include(cmake/Sources.cmake)

if (DALTON_ENABLE_XCFUN)
    add_subdirectory(xcfun)
    set(DALTON_INTERFACES_LIST
        ${DALTON_INTERFACES_LIST}
        xcfun
        )
    set(DALTON_LIBRARIES_TO_LINK
        ${DALTON_LIBRARIES_TO_LINK}
        xcfun
        )
endif ()

add_executable(dalton.x
    ${DALTON_SOURCES}
    )

target_link_libraries(
    dalton.x
    ${DALTON_LIBRARIES_TO_LINK}
    )
