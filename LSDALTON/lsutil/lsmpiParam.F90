module lsmpi_param
  use precision
#ifdef VAR_MPI
  use lsmpi_module
#endif
!!!!!!!!!!!!!!!!!!!!!!!!!
!Constants for MPIBUFFER!
!!!!!!!!!!!!!!!!!!!!!!!!!

  integer,parameter     :: LSMPIBROADCAST       = 1
  integer,parameter     :: LSMPIREDUCTION       = 2
  integer,parameter     :: LSMPIREDUCTIONmaster = 3
  integer,parameter     :: LSMPISENDRECV        = 4

#ifdef VAR_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!General MPI vars, aka junkbox!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer(kind=ls_mpik) :: MPI_COMM_LSDALTON = 0
  logical               :: LSMPIASYNCP                !contains environment value of async progress
  logical               :: lsmpi_enabled_comm_procs 

  !split mpi messages in case of 32bit mpi library to subparts, which are
  !describable by a 32bit integer and dividable by 8
  integer     :: SPLIT_MPI_MSG
  !split mpi one sided communication into 100MB chunks
  integer     :: MAX_SIZE_ONE_SIDED 

  !mpistatus
  integer(kind=ls_mpik) :: lsmpi_status(MPI_STATUS_SIZE) 

  type mpigroup
     integer(kind=ls_mpik)         :: groupsize
     integer(kind=ls_mpik),pointer :: ranks(:)
  end type mpigroup

#endif

  contains
    subroutine lsmpi_param_dummy()
      implicit none

    end subroutine lsmpi_param_dummy

end module lsmpi_param


