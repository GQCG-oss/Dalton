if(DEVELOPMENT_CODE AND NOT ENABLE_RELEASE)
    add_definitions(-DMOD_UNRELEASED)
endif()

if(ENABLE_GEN1INT)
    add_definitions(-DBUILD_GEN1INT)
endif()


if(ENABLE_DEC)
  add_definitions(-DVAR_DEC)
  set(ENABLE_TENSORS ON)
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

if(ENABLE_TITANBUILD)
   add_definitions(-DVAR_HAVE_MPI3)
   if(CMAKE_Fortran_COMPILER_ID MATCHES Cray) 
      if(ENABLE_TITANBUILD)
          add_definitions(-DVAR_WORKAROUND_CRAY_MEM_ISSUE_LARGE_ASSIGN)
      endif()
   endif()
endif()

add_definitions(-DBINARY_INFO_AVAILABLE)

if(cmake_build_type_tolower STREQUAL "debug")
  add_definitions(-DVAR_LSDEBUGINT)
  add_definitions(-DVAR_LSDEBUG)
  set(reorder_definitions " --debug_version ${reorder_definitions}")
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
    #WARNING THIS IS TEMPORARY 
    #I know that the combi --int64 and MKL will result
    #in linking to a 64 bit integer lapack. To my 
    #knowledge no other combi will produce this
    if(HAVE_MKL_LAPACK)
      add_definitions(-DVAR_LAPACK_INT64)
    endif()
endif()

if(ENABLE_GPU)
    add_definitions(-DVAR_CUDA)
    if(ENABLE_CUBLAS)
        add_definitions(-DVAR_CUBLAS)
    endif()
endif()

if(ENABLE_REAL_SP)
    add_definitions(-DVAR_REAL_SP)
endif()

if(ENABLE_DEBUGPBC)
    add_definitions(-DDEBUGPBC)
endif()

if(ENABLE_QCMATRIX)
    add_definitions(-DENABLE_QCMATRIX)
endif()
