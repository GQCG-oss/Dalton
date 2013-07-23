program mpi_test

   implicit none

#include "mpif.h"

   integer :: ierr

   call mpi_init(ierr)
   call mpi_finalize(ierr)

end program
