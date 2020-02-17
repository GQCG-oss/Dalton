!dalton_copyright_start
!
!
!dalton_copyright_end

! this module contains macros for MPI-one-sided communication calls.
!
! re-written by sknecht for Dalton - linkoeping jan 2014
! original collection in Dirac - sknecht October 2007
!
!
! ideas for window hints: window for get operations in 2e-integral codes: no_locks == .true. (MPI-3)
!
!
#ifdef VAR_MPI
module one_sided_communication_wrappers

#ifdef USE_MPI_MOD_F90
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif


  public mpixget
  public mpixaccum
  public mpixwincreate
  public mpixwinfree

  private
#include "priunit.h"
  integer(kind=MPI_ADDRESS_KIND) :: lower_bound
  integer(kind=MPI_ADDRESS_KIND) :: size_dp
  !---
  integer(kind=MPI_INTEGER_KIND) :: ierr
  integer(kind=MPI_INTEGER_KIND) :: my_win_mpi
  integer(kind=MPI_INTEGER_KIND) :: istat(mpi_status_size)
  integer(kind=MPI_INTEGER_KIND) :: jcount_mpi, win_lock_mode_mpi, jcount_t_mpi
  integer(kind=MPI_INTEGER_KIND) :: itarget_mpi
  integer(kind=MPI_INTEGER_KIND) :: datatype_out = mpi_real8
  integer(kind=MPI_INTEGER_KIND) :: datatype_in  = mpi_real8

contains

  subroutine mpixget(rbuf,          &
                     jcount,        &
                     itarget,       &
                     idispl,        &
                     jcount_t,      &
                     myid,          &
                     win_lock_mode, &
                     lock_active,   &
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
  integer, intent(in)                        :: jcount
  integer, intent(in)                        :: jcount_t
  integer, intent(in)                        :: myid
  integer, intent(in)                        :: itarget
  integer, intent(in)                        :: win_lock_mode
  integer, intent(inout)                     :: my_win
  logical, intent(in)                        :: lock_active
  integer(kind=MPI_ADDRESS_KIND), intent(in) :: idispl
!-------------------------------------------------------------------------------
!
      my_win_mpi        = my_win
      win_lock_mode_mpi = win_lock_mode
      itarget_mpi       = itarget
      jcount_mpi        = jcount
      jcount_t_mpi      = jcount_t

!     lock window (MPI_LOCK_SHARED mode)
      if(lock_active) call mpi_win_lock(win_lock_mode_mpi,itarget_mpi,mpi_mode_nocheck,my_win_mpi,ierr)

!     print *, 'my_win, itarget, myid, jcount,idispl,jcount_t,datatype_in',&
!               my_win, itarget, myid, jcount,idispl,jcount_t,datatype_in
        
!
!     transfer data     
      call mpi_get(rbuf,jcount_mpi,datatype_out,itarget_mpi,idispl,jcount_t_mpi,datatype_in,my_win_mpi,ierr)
!
!     unlock
      if(lock_active) call mpi_win_unlock(itarget_mpi,my_win_mpi,ierr)
!
  end subroutine mpixget
!*******************************************************************************

  subroutine mpixaccum(rbuf,          &
                       jcount,        &
                       itarget,       &
                       idispl,        &
                       jcount_t,      &
                       myid,          &
                       win_lock_mode, &
                       my_win,        &
                       lock_active)
!*******************************************************************************
!
!     accumulate a scalar/vector/matrix via a remote memory access (RMA) 
!     routine. Data are put from the origin memory to the target memory.
!
!     output: 
!            my_win:   updated target memory window with jcount_t elements.
!     input: 
!            rbuf  :   memory buffer on origin with JCOUNT elements.
!
!     allowed OPERATION: mpi_sum
!
!*******************************************************************************
  real(8), intent(in)                        :: rbuf(*)
  integer, intent(in)                        :: jcount
  integer, intent(in)                        :: jcount_t
  integer, intent(in)                        :: myid
  integer, intent(in)                        :: itarget
  integer, intent(in)                        :: win_lock_mode
  integer, intent(in)                        :: my_win
  integer(kind=MPI_ADDRESS_KIND), intent(in) :: idispl
  logical, intent(in), optional              :: lock_active
!-------------------------------------------------------------------------------
  logical                                    :: lock
!-------------------------------------------------------------------------------
!
      my_win_mpi        = my_win
      win_lock_mode_mpi = win_lock_mode
      itarget_mpi       = itarget
      jcount_mpi        = jcount
      jcount_t_mpi      = jcount_t

      lock = .true.
      if(present(lock_active)) lock = lock_active
!
!     lock window (MPI_LOCK_SHARED mode)
      if(lock) call mpi_win_lock(win_lock_mode_mpi,itarget_mpi,mpi_mode_nocheck,my_win_mpi,ierr)

!     accumulate data
      call mpi_accumulate(rbuf,        &
                          jcount_mpi,  &
                          datatype_in, & 
                          itarget_mpi, &
                          idispl,      &
                          jcount_t_mpi,&
                          datatype_out,&
                          mpi_sum,     &
                          my_win_mpi,  &
                          ierr)

!     unlock window
      if(lock) call mpi_win_unlock(itarget_mpi,my_win_mpi,ierr)

  end subroutine mpixaccum
!*******************************************************************************

  subroutine mpixwincreate(rbuf,nelement,myid,win_communicator,win_info,my_win)
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
  integer, intent(in)                        :: win_info
  integer, intent(out)                       :: my_win
  integer(kind=MPI_ADDRESS_KIND), intent(in) :: nelement
  integer(kind=MPI_ADDRESS_KIND)             :: buf_len
!-------------------------------------------------------------------------------
  integer(kind=MPI_INTEGER_KIND)             :: size_dp_local
  integer(kind=MPI_INTEGER_KIND)             :: win_info_mpi
  integer(kind=MPI_INTEGER_KIND)             :: win_comm_mpi
  integer                                    :: memory_model
  logical                                    :: flag
!-------------------------------------------------------------------------------
!
!     open memory window on each process shared by win_communicator
!     isize_dp is scaling unit in memory window
!     --> scale by REAL*8
!
      call mpi_type_get_extent(mpi_real8,lower_bound,size_dp,ierr)

      buf_len       = size_dp * nelement
      size_dp_local = size_dp
      win_info_mpi  = win_info
      win_comm_mpi  = win_communicator
!
!     write(*,*) 'opening window: my_win,myid', my_win,myid
      call mpi_win_create(rbuf,buf_len,size_dp_local,win_info_mpi,win_comm_mpi,my_win_mpi,ierr)
      my_win = my_win_mpi
!mpi-3call mpi_win_get_attr(my_win, mpi_win_model, memory_model, flag, ierr)
!
!     write(*,*) ' window opened: my_win,myid', my_win,myid
!mpi-3write(lupri,*) ' window memory model:', memory_model
!
  end subroutine mpixwincreate
!*******************************************************************************

  subroutine mpixwinfree(my_win, myid)
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
!     write(*,*) ' closing window', my_win,myid
      my_win_mpi = my_win
      call mpi_win_free(my_win_mpi,ierr)
!
!     write(*,*) ' window closed', my_win,myid
!
  end subroutine mpixwinfree
!*******************************************************************************

end module one_sided_communication_wrappers
#else
subroutine one_sided_communication_wrappers
! dummy routine for non-mpi compilation
end
#endif
