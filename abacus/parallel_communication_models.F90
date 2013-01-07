!dalton_copyright_start
!
!
!dalton_copyright_end
!
! stefan jan 2013: this module contains the type definition for the
!                  MPI communication handles to be potentially used in DALTON.
!                  currently in use: 
!                                    - lucita (after merge with hjaa-srdft)
!                                    - mcscf  (after merge with hjaa-srdft)
!                                    - hermit
!
module parallel_communication_models

#ifdef VAR_MPI
#ifndef VAR_USE_MPIF
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif
#else
  implicit none
#endif

  public communication_type
  public communication_init
  public communication_free


! communication type definition
! ----------------------------------------------------------------------------
  type communication_type

    integer ::                   &
      communication_intranode,   &             ! intra-node communication handle
      communication_internode,   &             ! inter-node communication handle
      communication_sharedijkl,  &             ! shared-memory communication handle: integrals, fock matrices, density matrices, etc

    logical ::                   &
      communication_type_init = .false.        ! status of the communication type

    integer, allocatable ::      &
      total_process_list(:),     &             ! list of all processes with their process ID
      shared_mem_group_list(:),  &             ! list of all processes within a shared-memory group and their process ID
      intra_node_group_list(:)                 ! list of all processes within an intra-node group and their process ID

  end type communication_type

! communication type object
  type(communication_type), public, save :: communication_info
! ----------------------------------------------------------------------------

  integer, private                       :: ierr
#ifdef VAR_MPI
  integer, private                       :: istat(MPI_STATUS_SIZE)
#endif

contains 

!*******************************************************************************
  subroutine communication_init(A, nr_of_process_glb, my_process_id_glb, global_communicator)

!   ----------------------------------------------------------------------------
!                - provides: 
!                           a. communication handles and process-ids in each group
!                           b. arrays with group listings
!                - requires: 
!                           a. global (starting) communication handle
!                           b. number of processors
!
!    for a detailed description of the communication model see references:
!    S. Knecht, H. J. Aa. Jensen, and T. Fleig
!       JCP, 128, 014108 (2008)
!       JCP, 132, 014108 (2010)
!   ----------------------------------------------------------------------------
    type(communication_type)   :: A
    integer, intent(in)        :: nr_of_process_glb
    integer, intent(in)        :: global_communicator
    integer, intent(in)        :: my_process_id_glb
!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------

!   reset old type information
    call communication_free(A)

    A%communication_type_init  = .true.
    A%communication_intranode  = -1
    A%communication_internode  = -1
    A%communication_sharedijkl = -1

    allocate(A%total_process_list(nr_of_process_glb))
    allocate(A%shared_mem_group_list(nr_of_process_glb))
    allocate(A%intra_node_group_list(nr_of_process_glb))
    A%total_process_list      = -1
    A%intra_node_group_list   = -1

#ifdef VAR_MPI
!   start communication model setup
!   -------------------------------

!   1. determine # communication process groups and store group-IDs in total_process_list
    call set_communication_groups(A%total_process_list,     &
                                  my_process_id_glb,        &
                                  nr_of_process_glb,        &
                                  global_communicator)

!   2. setup communicators and process-id for the various communication levels
!      a. intra-node
!      b. inter-node
!      c. shared-memory
    call set_communication_levels(A%intra_node_group_list,                   &
                                  A%shared_mem_group_list,                   &
                                  A%total_process_list,                      &
                                  nr_of_process_glb,                         &
                                  my_process_id_glb,                         &
                                  communicator_glb,                          &
                                  my_intra_node_id,                          &
                                  my_inter_node_id,                          &
                                  my_shmem_ijkl_id,                          &
                                  intra_node_size,                           &
                                  inter_node_size,                           &
                                  shmem_ijkl_size,                           &
                                  intra_node_comm,                           &
                                  inter_node_comm,                           &
                                  shmem_ijkl_comm,                           &
                                  intra_node_group_id)                      
#endif

  end subroutine communication_init
!*******************************************************************************

  subroutine communication_free(A)

!   ----------------------------------------------------------------------------
    type(communication_type) :: A
    integer                 :: ierr
!   ----------------------------------------------------------------------------

    if(.not. A%communication_type_init) return

    A%communication_type_init = .false.
    deallocate(A%total_process_list)
    deallocate(A%shared_mem_group_list)
    deallocate(A%intra_node_group_list)

#ifdef VAR_MPI
    call mpi_comm_free(A%communication_intranode,  ierr)
    call mpi_comm_free(A%communication_internode,  ierr)
    call mpi_comm_free(A%communication_sharedijkl, ierr)
#endif

  end subroutine communication_free
!*******************************************************************************

  subroutine set_communication_groups(process_list_glb,        &
                                      my_process_id_glb,       &
                                      nr_of_process_glb,       &
                                      communicator_glb)
!*******************************************************************************
     integer, intent(out)             :: process_list_glb(nr_of_process_glb)
     integer, intent(in )             :: nr_of_process_glb
     integer, intent(in )             :: my_process_id_glb
     integer, intent(in )             :: communicator_glb
!-------------------------------------------------------------------------------
     integer                          :: process_name_length  
     integer                          :: local_counter_file_groups  
     integer                          :: current_process_id  
     integer                          :: proc_id  
     integer                          :: finished_loop  
     character (len=255)              :: process_name
     character (len=255)              :: scr_process_name
     integer,             allocatable :: scr_arr_name_length(:)
     character (len=255), allocatable :: scr_arr_process_name(:)
!-------------------------------------------------------------------------------
                                               
      process_list_glb = -1

!     find system-dependent unique process name
      call mpixprocname(process_name,process_name_length)

      allocate(scr_arr_name_length(nr_of_process_glb))
      scr_arr_name_length = 0

!     1. gather all name length
      call mpi_allgather(process_name_length,1,my_MPI_INTEGER,         &
     &                   scr_arr_name_length,1,my_MPI_INTEGER,         &
     &                   communicator_glb,ierr)

!     2. collect all names in temporary storage
      allocate(scr_arr_process_name(nr_of_process_glb*255))
      scr_arr_process_name(my_process_id_glb+1) = process_name

      do proc_id = 1, nr_of_process_glb

         scr_process_name = process_name
         call mpi_bcast(scr_process_name,                              &
                        scr_arr_name_length(proc_id),                  &
                        MPI_CHARACTER,                                 &
                        proc_id-1,                                     &
                        communicator_glb,                              &
                        ierr)

         if(my_process_id_glb /= proc_id -1)then
           scr_arr_process_name(proc_id) = scr_process_name(1:scr_arr_name_length(proc_id))
         end if

      end do

 
!     3. find all processors on the same deck and reorder (if necessary) 
!     to get the processors as close as possible, starting with the master (id == 0)
!     NOTE: reordering is currently NOT done
      current_process_id        = 1
      local_counter_file_groups = 0

      do ! indefinite loop - we check at the end if we meet an exit condition

        finished_loop = 1

        local_counter_file_groups = local_counter_file_groups + 1
        process_name = scr_arr_process_name(current_process_id)         &
                       (1:scr_arr_name_length(current_process_id))
 
        do proc_id = 1, nr_of_process_glb
 
          scr_process_name(1:scr_arr_name_length(proc_id)) =            &
          scr_arr_process_name(proc_id)                                 &
          (1:scr_arr_name_length(proc_id))
 
          if(scr_process_name(1:scr_arr_name_length(proc_id)) ==        &
             process_name(1:scr_arr_name_length(current_process_id)))   &
             process_list_glb(proc_id) = current_process_id
 
        end do
!
!       check if we are done and set finished_loop: done = 1, else 0
!       search for the next lowest cpu building the 'local group master'
!       which is current_process_id
!
        do proc_id = 1, nr_of_process_glb
          if(finished_loop /= 0 .and. process_list_glb(proc_id) == -1)then 
            finished_loop = 0
            current_process_id = proc_id
          end if
        end do
        if(finished_loop == 1) exit

      end do
 
      deallocate(scr_arr_process_name)
      deallocate(scr_arr_name_length)

!     master reports the count of intra-groups formed by all processes
      if(my_process_id_glb == 0)then
        write(*,'(/a,i4,1x,a)')                                             &
        '  *** final output from the group generator:',                     &
        local_counter_file_groups,' intra-node group(s) (is) are built. ***'
      end if

  end subroutine set_communication_groups
!*******************************************************************************

  subroutine set_communication_levels(group_list,                              &
                                      process_list_glb,                        &
                                      process_list_shared_mem_glb,             &
                                      nr_of_process_glb,                       &
                                      my_process_id_glb,                       &
                                      communicator_glb,                        &
                                      my_intra_node_id,                        &
                                      my_inter_node_id,                        &
                                      my_shmem_ijkl_id,                        &
                                      my_shmem_cvec_id,                        &
                                      intra_node_master,                       &
                                      shmem_master_ijkl,                       &
                                      intra_node_size,                         &
                                      inter_node_size,                         &
                                      shmem_ijkl_size,                         &
                                      intra_node_comm,                         &
                                      inter_node_comm,                         &
                                      shmem_ijkl_comm,                         &
                                      intra_node_group_id)                      
!********************************************************************************
!     purpose: setup communicators and process-id for the various               
!              communication levels:
!                a. intra-node
!                b. inter-node
!                c. shared-memory
!*******************************************************************************
     integer, intent(in )   :: nr_of_process_glb
     integer, intent(in )   :: communicator_glb
     integer, intent(in )   :: my_process_id_glb
     integer, intent(in )   :: process_list_glb(nr_of_process_glb)
     integer, intent(in )   :: process_list_shared_mem_glb(nr_of_process_glb)
     integer, intent(in )   :: shared_memory_lvl_ijkl
     integer, intent(in )   :: shared_memory_lvl_cvec
     integer, intent(out)   :: group_list(nr_of_process_glb)
     integer, intent(out)   :: my_intra_node_id
     integer, intent(out)   :: my_inter_node_id
     integer, intent(out)   :: my_shmem_ijkl_id
     integer, intent(out)   :: my_shmem_cvec_id
     integer, intent(out)   :: shmem_master_ijkl
     integer, intent(out)   :: intra_node_size
     integer, intent(out)   :: inter_node_size
     integer, intent(out)   :: shmem_ijkl_size
     integer, intent(out)   :: intra_node_comm
     integer, intent(out)   :: inter_node_comm
     integer, intent(out)   :: shmem_ijkl_comm
     integer, intent(out)   :: intra_node_group_id
     integer, intent(out)   :: intra_node_master
!-------------------------------------------------------------------------------
     integer                :: key
     integer                :: color
     integer                :: tmp_group_counter
     integer                :: current_proc
!-------------------------------------------------------------------------------

!     a. intra-node communicator
 
      key   = my_process_id_glb
      color = process_list_glb(my_process_id_glb+1)

      call build_new_communicator_group(communicator_glb,        &
                                        intra_node_comm,         &
                                        intra_node_size,         &
                                        my_intra_node_id,        &
                                        color,                   &
                                        key)

!     setup required information about each intra-node group:
!       - group id
!       - group list
      intra_node_group_id = color
      tmp_group_counter   = 1
      do current_proc = 1, nr_of_process_glb
        if(process_list_glb(current_proc) == intra_node_group_id)then
          group_list(tmp_group_counter) = current_proc      - 1
          tmp_group_counter             = tmp_group_counter + 1
        end if
      end do

!     b. inter-node communicator

      key   = process_list_glb(my_process_id_glb+1)
      color = 2
      if(my_intra_node_id /= 0) color = 3
 
      call build_new_communicator_group(communicator_glb,        &
                                        inter_node_comm,         &
                                        inter_node_size,         &
                                        my_inter_node_id,        &
                                        color,                   &
                                        key)
!     c. shared memory communicator

      key   = my_process_id_glb
      color = process_list_shared_mem_glb(my_process_id_glb+1)
      if(shared_memory_lvl_ijkl == 1) color = 1
 
      call build_new_communicator_group(communicator_glb,        &
                                        shmem_ijkl_comm,         &
                                        shmem_ijkl_size,         &
                                        my_shmem_ijkl_id,        &
                                        color,                   &
                                        key)

!     set new sub-group master (might be master of all CPU's)
      intra_node_master = 0
      shmem_master_ijkl = 0

  end subroutine set_communication_levels

!*******************************************************************************

  subroutine build_new_communication_group(old_communicator,      &
                                           new_communicator,      &
                                           group_size,            &
                                           my_id_group,           &
                                           color,                 &
                                           key)
!*******************************************************************************
!     purpose: split old communication group into new sub-groups 
!              using key and color. 
!              the new communication along with the process-ids in the new group 
!              and the group size are returned.
!*******************************************************************************
     integer, intent(in )   :: old_communication
     integer, intent(in )   :: color
     integer, intent(in )   :: key
     integer, intent(out)   :: group_size
     integer, intent(out)   :: my_id_group
     integer, intent(out)   :: new_communication
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      call mpi_comm_split(old_communication,color,key,new_communication,ierr)
 
!     collect required information about each group:
!       - group size
!       - process id in the group
      call mpi_comm_size(new_communication, group_size,ierr)
      call mpi_comm_rank(new_communication,my_id_group,ierr)

  end subroutine build_new_communication_group

end module parallel_communication_models
