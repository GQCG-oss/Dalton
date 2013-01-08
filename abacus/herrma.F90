!dalton_copyright_start
!
!
!dalton_copyright_end

! this module contains all functionality wrt RMA operations (MPI-one-sided communication) in hermit.
!
! written by sknecht for Dalton - linkoeping jan 2013
!
!
module herrma

#ifdef VAR_MPI
#ifndef VAR_USE_MPIF
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif


  public set_rma_windows
  public free_rma_windows

#include "priunit.h"
  integer                        :: ierr
  integer(kind=MPI_ADDRESS_KIND) :: lower_bound
  integer(kind=MPI_ADDRESS_KIND) :: size_dp
  integer                        :: istat(mpi_status_size)

contains

  subroutine set_rma_windows(rbuf,          &
                             my_win)
!*******************************************************************************
!
!     Communicate a scalar/vector/matrix via a remote memory access (RMA).
!     Data are put from the target memory to the origin.
!
!     output: 
!            rbuf:      updated origin buffer with jcount elements.
!     input: 
!            my_win:    memory window on itarget (must be initialized!)
!                       accessed at displacement idispl with jcount_t 
!                       elements.
!*******************************************************************************
  real(8), intent(inout)                     :: rbuf(*)
  integer, intent(in)                        :: my_win
!-------------------------------------------------------------------------------
  integer                                    :: datatype_out = mpi_real8
  integer                                    :: datatype_in  = mpi_real8
!-------------------------------------------------------------------------------
!
  end subroutine set_rma_windows
!*******************************************************************************

  subroutine free_rma_windows(rbuf,          &
                              my_win)
!*******************************************************************************
!
!     Communicate a scalar/vector/matrix via a remote memory access (RMA).
!     Data are put from the target memory to the origin.
!
!     output: 
!            rbuf:      updated origin buffer with jcount elements.
!     input: 
!            my_win:    memory window on itarget (must be initialized!)
!                       accessed at displacement idispl with jcount_t 
!                       elements.
!*******************************************************************************
  real(8), intent(inout)                     :: rbuf(*)
  integer, intent(in)                        :: my_win
!-------------------------------------------------------------------------------
  integer                                    :: datatype_out = mpi_real8
  integer                                    :: datatype_in  = mpi_real8
!-------------------------------------------------------------------------------
!
  end subroutine free_rma_windows
!*******************************************************************************

#endif
!*******************************************************************************

end module herrma
