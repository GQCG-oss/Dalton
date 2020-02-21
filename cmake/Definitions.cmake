if(ENABLE_GEN1INT)
    add_definitions(-DBUILD_GEN1INT)
endif()

if(ENABLE_CHEMSHELL)
    add_definitions(-DVAR_CHEMSHELL)
endif()

if(ENABLE_QFITLIB)
    add_definitions(-DBUILD_QFITLIB)
endif()

add_definitions(-DVAR_MFDS)
add_definitions(-D_FILE_OFFSET_BITS=64)
add_definitions(-DIMPLICIT_NONE)

add_definitions(-DBINARY_INFO_AVAILABLE)

add_definitions(-DINSTALL_BASDIR="${PROJECT_BINARY_DIR}/basis")

if(HAVE_MKL_LAPACK OR HAVE_MKL_BLAS)
    add_definitions(-DVAR_MKL)
endif()

if(ENABLE_LSEEK)
    add_definitions(-DHAVE_NO_LSEEK64)
endif()

if(ENABLE_64BIT_INTEGERS)
    add_definitions(-DVAR_INT64)
endif()

