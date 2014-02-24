!> @file
!> Contains module for MPI integral evaluation settings
module infpar_module
use precision
#ifdef VAR_MPI

type infpar_struct
  integer :: nodeid(129), ncode,  iprpar, mtottk, ntask,&
       &nfmat,  ndegdi, mytid,  timing,   slave,&
       &debug
  integer :: inputBLOCKSIZE
  character(20) :: nodnam(128), myname
  !> Master rank
  integer(kind=ls_mpik) :: master
  !> Rank inside total (world) group
  integer(kind=ls_mpik) :: mynum
  !> Total number of nodes inside total (world) group
  integer(kind=ls_mpik) :: nodtot
  !> Rank inside local group where one of the "slaves" is a local master
  integer(kind=ls_mpik) :: lg_mynum
  !> Total number of nodes local group where one of the "slaves" is a local master
  integer(kind=ls_mpik) :: lg_nodtot
  !> Communicator for local group where one of the "slaves" is a local master
  integer(kind=ls_mpik) :: lg_comm
  !> Communicator structure to distinguish spawned thread from original thread,
  !for the parent this has the value MPI_COMM_NULL
  integer(kind=ls_mpik) :: parent_comm
  integer(kind=ls_mpik) :: child_comm
  !> intracommunicator of child and parent thread, the intercommunicator is
  !deleted
  integer(kind=ls_mpik) :: pc_comm
  !> rank in the parent child communicator
  integer(kind=ls_mpik) :: pc_mynum
  !> number of nodes in the parent child communicator
  integer(kind=ls_mpik) :: pc_nodtot

  !> Are there more local jobs?
  logical :: lg_morejobs
end type infpar_struct

type(infpar_struct) :: infpar
#endif
contains

!Added to avoid "has no symbols" linking warning
subroutine infpar_module_void()
end subroutine infpar_module_void

end module infpar_module
