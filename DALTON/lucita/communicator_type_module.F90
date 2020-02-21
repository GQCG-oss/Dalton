!dalton_copyright_start
!
!
!dalton_copyright_end
!
! stefan dec 2011: this module contains the type definition for the
!                  MPI communcator handle used in the MCSCF-LUCITA interface + LUCITA program.
!
module communicator_type_module

  implicit none

  public communicator_type
  public communicator_init_lucipar
  public communicator_free_lucipar
  public communicator_switch_lucipar


! type definition
  type communicator_type

    integer ::                  &
      communicator_intranode,   &             ! intra-node communicator handle
      communicator_internode,   &             ! intra-node communicator handle
      communicator_sharedijkl,  &             ! shared-memory communicator handle: (ij|kl)
      communicator_sharedcvec                 ! shared-memory communicator handle: cvec

    logical ::                  &
      communicator_type_init = .false.        ! status of the communicator type

    integer, allocatable ::     &
      total_process_list(:),    &             ! list of all processes with their process ID
      intra_node_group_list(:)                ! list of all processes within an intra-node group and their process ID

  end type communicator_type

! ttss block type object
  type(communicator_type), public, save :: communicator_info

contains 

  subroutine communicator_init_lucipar(A, number_of_processes)

!   ----------------------------------------------------------------------------
    type(communicator_type)    :: A
    integer, intent(in)        :: number_of_processes
!   ----------------------------------------------------------------------------

!   reset old type information
    call communicator_free_lucipar(A)

    A%communicator_type_init  = .true.
    A%communicator_intranode  = -1
    A%communicator_internode  = -1
    A%communicator_sharedijkl = -1
    A%communicator_sharedcvec = -1

    allocate(A%total_process_list(number_of_processes))
    allocate(A%intra_node_group_list(number_of_processes))
    A%total_process_list      = -1
    A%intra_node_group_list   = -1

  end subroutine communicator_init_lucipar

  subroutine communicator_free_lucipar(A)

!   ----------------------------------------------------------------------------
    type(communicator_type) :: A
    integer                 :: ierr
!   ----------------------------------------------------------------------------

    if(.not. A%communicator_type_init) return

    A%communicator_type_init = .false.
    deallocate(A%total_process_list)
    deallocate(A%intra_node_group_list)

#ifdef VAR_MPI
    call mpi_comm_free(A%communicator_intranode,  ierr)
    call mpi_comm_free(A%communicator_internode,  ierr)
    call mpi_comm_free(A%communicator_sharedijkl, ierr)
    call mpi_comm_free(A%communicator_sharedcvec, ierr)
#endif

  end subroutine communicator_free_lucipar

  subroutine communicator_switch_lucipar(A,intra,inter,share_ijkl,share_cvec)

!   ----------------------------------------------------------------------------
    type(communicator_type)  :: A
    integer, intent(in) :: intra,inter,share_ijkl,share_cvec
!   ----------------------------------------------------------------------------

    A%communicator_intranode  = intra
    A%communicator_internode  = inter
    A%communicator_sharedijkl = share_ijkl
    A%communicator_sharedcvec = share_cvec

  end subroutine communicator_switch_lucipar

end module communicator_type_module
