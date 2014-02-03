if(DEVELOPMENT_CODE AND NOT ENABLE_RELEASE)
    add_definitions(-DMOD_UNRELEASED)
endif()

if(ENABLE_GEN1INT)
    add_definitions(-DBUILD_GEN1INT)
#   currently not set, don't know why
#   add_definitions(-DBUILD_GEN1INT_LSDALTON)
endif()

if(ENABLE_CHEMSHELL)
    add_definitions(-DVAR_CHEMSHELL)
endif()

add_definitions(-DVAR_MFDS)
add_definitions(-D_FILE_OFFSET_BITS=64)
add_definitions(-DIMPLICIT_NONE)

add_definitions(-DBINARY_INFO_AVAILABLE)

if(cmake_build_type_tolower STREQUAL "debug")
    add_definitions(-DVAR_LSDEBUGINT)
    add_definitions(-DVAR_LSDEBUG)
endif()

add_definitions(-DINSTALL_BASDIR="${PROJECT_BINARY_DIR}/basis")

if(HAVE_MKL_LAPACK)
    add_definitions(-DVAR_MKL)
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
    add_definitions(-DVAR_CSR)
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
