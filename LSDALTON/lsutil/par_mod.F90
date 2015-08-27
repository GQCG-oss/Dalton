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
  !> number of nodes used in scalapack
  integer(kind=ls_mpik) :: ScalapackGroupSize
  !> Use scalapack workaround
  logical :: ScalapackWorkAround
  !> number of nodes used in PDMM
  integer(kind=ls_mpik) :: pdmmGroupSize

  !> Are there more local jobs?
  logical :: lg_morejobs
end type infpar_struct

type(infpar_struct) :: infpar
!> specifies if mpi_init should be called in lsmpi_init
logical, save :: call_mpi_init = .true.
#endif
contains

!Added to avoid "has no symbols" linking warning
subroutine infpar_module_void()
end subroutine infpar_module_void

end module infpar_module
