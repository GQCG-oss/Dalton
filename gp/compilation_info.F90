   subroutine print_compilation_info()
      implicit none
#ifdef THIS_IS_CMAKE_BUILD
!#include "priunit.h"
!#include "compinfo.h"
!     write(lupri, *) 'CMAKE_SYSTEM_NAME          ',
!     write(lupri, *)
!     write(lupri, *) 'ENABLE_INTERNAL_MATH       ',
!     write(lupri, *)
!     write(lupri, *) 'ENABLE_64BIT_INTEGERS      ',
!     write(lupri, *)
!     write(lupri, *) 'ENABLE_MPI                 ',
!     write(lupri, *)
!     write(lupri, *) 'ENABLE_XCFUN               ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_BUILD_TYPE           ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_Fortran_COMPILER     ',
!     write(lupri, *)
!     write(lupri, *) 'FORTRAN_COMPILER_VERSION   ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_Fortran_FLAGS        ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_Fortran_FLAGS_DEBUG  ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_Fortran_FLAGS_RELEASE',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_C_COMPILER           ',
!     write(lupri, *)
!     write(lupri, *) 'C_COMPILER_VERSION         ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_C_FLAGS              ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_C_FLAGS_DEBUG        ',
!     write(lupri, *)
!     write(lupri, *) 'CMAKE_C_FLAGS_RELEASE      ',
!     write(lupri, *)
!     write(lupri, *) 'BLAS_LIBRARIES             ',
!     write(lupri, *)
!     write(lupri, *) 'LAPACK_LIBRARIES           ',
!     write(lupri, *)
!     write(lupri, *) 'SVN_REVISION               ',
#endif
   end subroutine
