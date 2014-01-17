program raboof
!  this program won't compile if compilers do not match
   implicit none
contains

   subroutine foo()

      use mpi

      integer :: request
      integer :: status_container(MPI_STATUS_SIZE)
      integer :: ierr

      ! intel mpi 64 currently (2013-04-22) does not support mpi.mod
      ! the following call will not resolve and make this test fail
      ! and force mpif.h
      call mpi_wait(request, status_container, ierr)

   end subroutine
end program
