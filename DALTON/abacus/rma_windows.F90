!dalton_copyright_start
!
!
!dalton_copyright_end

! RMA windows object (used for MPI-one-sided communication) and driver routines to 
! open and close memory windows.
!
! written by sknecht - linkoeping jan 2014
!
!
module rma_windows

#ifdef VAR_MPI
  use one_sided_communication_wrappers, only: &
      mpixwincreate,                          &
      mpixwinfree
#endif

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif
#else
  implicit none
#endif

#ifdef VAR_MPI
  public set_rma_window
  public free_rma_window
#endif

  private

  save

! rma-windows definition
! ----------------------------------------------------------------------------
  type, public :: rma_win

    integer ::                   &
      dmat_win,                  &             ! density matrix window (2e-integral codes)
      fmat_win                                 ! fock matrix window (2e-integral codes)

    integer ::                   &
      nmat_max_wo_win,           &             ! max number of matrices "outside" the memory window (no RMA communication needed for those)
      dmat_fh                                  ! local dmat shared file handle

    logical ::                   &
      rma_win_init    = .false.                ! status of the rma_win type
    logical ::                   &
      lock_win        = .true.                 ! lock windows for communication calls

  end type rma_win

! rma-windows object
  type(rma_win), public :: rma_win_info
! ----------------------------------------------------------------------------

#ifdef VAR_MPI

contains

  subroutine set_rma_window(rbuf,                &
                            nelement,            &
                            myid,                &
                            win_communicator,    &
                            my_win,              &
                            set_win_lock_info,   &
                            key,                 &
                            value,               &
                            external_info        &
                           )
!*******************************************************************************
!
!     open memory window to be used in one-sided MPI communication
!
!     INPUT:
!            array RBUF (should be allocated by mpi_alloc_mem to
!            assure lock functionality)
!            number of elements nelement of type real*8
!
!
!     OUTPUT:
!            new memory window handle my_win shared by all
!            processes in the communication group win_communicator.
!            extension of memory window on each process may be
!            unsymmetric.
!
!
!*******************************************************************************
  real(8), intent(inout)                     :: rbuf(*)
  integer, intent(in)                        :: myid
  integer, intent(in)                        :: win_communicator
  logical, intent(in)                        :: set_win_lock_info
  integer, intent(out)                       :: my_win
  integer(kind=MPI_ADDRESS_KIND), intent(in) :: nelement
  character*(*), intent(in), optional        :: key
  character*(*), intent(in), optional        :: value
  integer, intent(in),       optional        :: external_info
!-------------------------------------------------------------------------------
  integer                                    :: info_object
  integer(kind=MPI_INTEGER_KIND)             :: info_object_mpi, myid_mpi, win_comm_mpi
  integer(kind=MPI_INTEGER_KIND)             :: ierr
!-------------------------------------------------------------------------------
!
!     open memory window on each process shared by win_communicator
      info_object_mpi = mpi_info_null

      if(set_win_lock_info)then
        call mpi_info_create(info_object_mpi,ierr)
        call mpi_info_set(info_object_mpi,key,value,ierr)
      else
        if(present(external_info)) then
           info_object_mpi = external_info
        end if
      end if

      info_object = info_object_mpi
      call  mpixwincreate(rbuf,                &
                          nelement,            &
                          myid,                &
                          win_communicator,    &
                          info_object,         &
                          my_win)

      if(set_win_lock_info)then
        call mpi_info_free(info_object_mpi,ierr)
      end if
!
  end subroutine set_rma_window
!*******************************************************************************

  subroutine free_rma_window(my_win, myid)
!*******************************************************************************
!
!     close memory window used in one-sided MPI communication
!
!     INPUT: 
!            memory window handle my_win shared by all processes in
!            a communication group.
!
!*******************************************************************************
  integer, intent(in)    :: myid
  integer, intent(inout) :: my_win
!-------------------------------------------------------------------------------
!
!     close memory window on each process within the communication group
!
      call mpixwinfree(my_win,myid)
!
!
  end subroutine free_rma_window 
#endif
!*******************************************************************************

end module rma_windows
