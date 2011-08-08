   subroutine print_compilation_info()
      implicit none
#ifdef THIS_IS_CMAKE_BUILD
#include "priunit.h"
#include "compinfo.h"
      write(lupri, *) 'CMAKE_SYSTEM_NAME          ', CMAKE_SYSTEM_NAME
      write(lupri, *)
      write(lupri, *) 'ENABLE_INTERNAL_MATH       ', ENABLE_INTERNAL_MATH
      write(lupri, *)
      write(lupri, *) 'ENABLE_64BIT_INTEGERS      ', ENABLE_64BIT_INTEGERS
      write(lupri, *)
      write(lupri, *) 'ENABLE_MPI                 ', ENABLE_MPI
      write(lupri, *)
      write(lupri, *) 'ENABLE_XCFUN               ', ENABLE_XCFUN
      write(lupri, *)
      write(lupri, *) 'CMAKE_BUILD_TYPE           ', CMAKE_BUILD_TYPE
      write(lupri, *)
      write(lupri, *) 'CMAKE_Fortran_COMPILER     ', CMAKE_Fortran_COMPILER
      write(lupri, *)
      write(lupri, *) 'FORTRAN_COMPILER_VERSION   ', FORTRAN_COMPILER_VERSION
      write(lupri, *)
      write(lupri, *) 'CMAKE_Fortran_FLAGS        ', CMAKE_Fortran_FLAGS
      write(lupri, *)
      write(lupri, *) 'CMAKE_Fortran_FLAGS_DEBUG  ', CMAKE_Fortran_FLAGS_DEBUG
      write(lupri, *)
      write(lupri, *) 'CMAKE_Fortran_FLAGS_RELEASE', CMAKE_Fortran_FLAGS_RELEASE
      write(lupri, *)
      write(lupri, *) 'CMAKE_C_COMPILER           ', CMAKE_C_COMPILER
      write(lupri, *)
      write(lupri, *) 'C_COMPILER_VERSION         ', C_COMPILER_VERSION
      write(lupri, *)
      write(lupri, *) 'CMAKE_C_FLAGS              ', CMAKE_C_FLAGS
      write(lupri, *)
      write(lupri, *) 'CMAKE_C_FLAGS_DEBUG        ', CMAKE_C_FLAGS_DEBUG
      write(lupri, *)
      write(lupri, *) 'CMAKE_C_FLAGS_RELEASE      ', CMAKE_C_FLAGS_RELEASE
      write(lupri, *)
      write(lupri, *) 'BLAS_LIBRARIES             ', BLAS_LIBRARIES
      write(lupri, *)
      write(lupri, *) 'LAPACK_LIBRARIES           ', LAPACK_LIBRARIES
      write(lupri, *)
      write(lupri, *) 'SVN_REVISION               ', SVN_REVISION
#endif
   end subroutine
