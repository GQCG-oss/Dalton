!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module communication_model

! stefan: - this module provides all necessary functionality
!           to setup a communication model in parallel mcscf/ci 
!           calculations.
!
!           written by sknecht, may 2007 for DIRAC MCSCF/KR-CI/LUCITA
!           adapted for DALTON by sknecht, november 2010.
#ifndef VAR_USE_MPIF
  use mpi
  implicit none
#else
  implicit none
#include "mpif.h"
#endif

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

  public setup_communication_model
  public close_communication_model

  private

  save

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr

contains 
 
  subroutine setup_communication_model(nr_file_groups,                 &
                                       io_mode,                        &
                                       shared_memory_mode,             &
                                       shared_memory_lvl_ijkl,         &
                                       shared_memory_lvl_cvec,         &
                                       group_list,                     &
                                       process_list_glb,               &
                                       process_list_shared_mem_glb,    &
                                       my_process_id_glb,              &
                                       nr_of_process_glb,              &
                                       communicator_glb,               &
                                       my_intra_node_id,               &
                                       my_inter_node_id,               &
                                       my_shmem_ijkl_id,               &
                                       my_shmem_cvec_id,               &
                                       intra_node_size,                &
                                       inter_node_size,                &
                                       shmem_ijkl_size,                &
                                       shmem_cvec_size,                &
                                       intra_node_comm,                &
                                       inter_node_comm,                &
                                       shmem_ijkl_comm,                &
                                       shmem_cvec_comm,                &
                                       intra_node_group_id,            &                
                                       intra_node_master,              &
                                       shmem_master_ijkl,              &
                                       shmem_master_cvec,              &
                                       print_unit)
!******************************************************************************
!
!    purpose:  initialize communication model for parallel CI/MCSCF runs:
!                - provides: 
!                           a. communication handles and process-ids in each group
!                           b. arrays with group listings
!                - requires: 
!                           a. global (starting) communication handle
!                           a. number of processors
!                           b. information concerning shared-memory mode(s)
!
!    for a detailed description of the communication model see references:
!    S. Knecht, H. J. Aa. Jensen, and T. Fleig
!       JCP, 128, 014108 (2008)
!       JCP, 132, 014108 (2008)
!
!*******************************************************************************
     integer, intent(in )   :: nr_of_process_glb
     integer, intent(in )   :: communicator_glb
     integer, intent(in )   :: my_process_id_glb
     integer, intent(in )   :: print_unit
     integer, intent(out)   :: group_list(nr_of_process_glb)
     integer, intent(out)   :: process_list_glb(nr_of_process_glb)
     integer, intent(out)   :: process_list_shared_mem_glb(nr_of_process_glb)
     integer, intent(out)   :: nr_file_groups
     integer, intent(out)   :: io_mode
     integer, intent(out)   :: shared_memory_lvl_ijkl
     integer, intent(out)   :: shared_memory_lvl_cvec
     logical, intent(out)   :: shared_memory_mode
     integer, intent(out)   :: shmem_master_ijkl
     integer, intent(out)   :: shmem_master_cvec
     integer, intent(out)   :: intra_node_master
     integer, intent(out)   :: intra_node_group_id
     integer, intent(out)   :: my_intra_node_id
     integer, intent(out)   :: my_inter_node_id
     integer, intent(out)   :: my_shmem_ijkl_id
     integer, intent(out)   :: my_shmem_cvec_id
     integer, intent(out)   :: intra_node_size
     integer, intent(out)   :: inter_node_size
     integer, intent(out)   :: shmem_ijkl_size
     integer, intent(out)   :: shmem_cvec_size
     integer, intent(out)   :: intra_node_comm
     integer, intent(out)   :: inter_node_comm
     integer, intent(out)   :: shmem_ijkl_comm
     integer, intent(out)   :: shmem_cvec_comm
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     initialize number of process groups sharing, e.g. a c-vector file
      nr_file_groups              =  0
!     parallel i/o mode in general (1 = MPI-I/O)
      io_mode                     =  1
!     set default values (to be consistent with Dirac/KR-CI where these values are dynamic.) 
      shared_memory_lvl_ijkl      = -1
      shared_memory_lvl_cvec      = -1
      process_list_shared_mem_glb =  1 
      shared_memory_mode          = .false.

!     1. determine # process groups and store group-IDs in process_list_glb
      call set_file_groups(process_list_glb,          &
                           nr_file_groups,            &
                           my_process_id_glb,         &
                           nr_of_process_glb,         &
                           communicator_glb,          &
                           print_unit)

!     2. setup communicators and process-id for the various communication levels
!        a. intra-node
!        b. inter-node
!        c. shared-memory
      call set_communication_levels(group_list,                  &
                                    process_list_glb,            &
                                    process_list_shared_mem_glb, &
                                    nr_of_process_glb,           &
                                    my_process_id_glb,           &
                                    communicator_glb,            &
                                    shared_memory_lvl_ijkl,      &
                                    shared_memory_lvl_cvec,      &
                                    my_intra_node_id,            &
                                    my_inter_node_id,            &
                                    my_shmem_ijkl_id,            &
                                    my_shmem_cvec_id,            &
                                    intra_node_master,           &
                                    shmem_master_ijkl,           &
                                    shmem_master_cvec,           &
                                    intra_node_size,             &
                                    inter_node_size,             &
                                    shmem_ijkl_size,             &
                                    shmem_cvec_size,             &
                                    intra_node_comm,             &
                                    inter_node_comm,             &
                                    shmem_ijkl_comm,             &
                                    shmem_cvec_comm,             &
                                    intra_node_group_id)                      

  end subroutine setup_communication_model
!*******************************************************************************

  subroutine set_file_groups(process_list_glb,        &
                             nr_file_groups,          &
                             my_process_id_glb,       &
                             nr_of_process_glb,       &
                             communicator_glb,        &
                             print_unit)
!*******************************************************************************
     integer, intent(out)   :: process_list_glb(nr_of_process_glb)
     integer, intent(out)   :: nr_file_groups  
     integer, intent(in )   :: nr_of_process_glb
     integer, intent(in )   :: my_process_id_glb
     integer, intent(in )   :: communicator_glb
     integer, intent(in )   :: print_unit
!-------------------------------------------------------------------------------
     integer                          :: process_name_length  
     integer                          :: local_counter_file_groups  
     integer                          :: current_process_id  
     integer                          :: proc_id  
     integer                          :: finished_loop  
     character (len=100)              :: process_name
     character (len=100)              :: scr_process_name
     integer,             allocatable :: scr_arr_name_length(:)
     character (len=100), allocatable :: scr_arr_process_name(:)
!-------------------------------------------------------------------------------
                                               
      process_list_glb = -1

!     set this flag for parallel global filesystems - e.g. on IBM clusters
#ifdef VAR_PFS
      nr_file_groups   =  1
      process_list_glb =  1
#else
                                               
!     find system-dependent unique process name
      call mpixprocname(process_name,process_name_length)

      allocate(scr_arr_name_length(nr_of_process_glb))
      scr_arr_name_length = 0

!     1. gather all name length
      call mpi_allgather(process_name_length,1,my_MPI_INTEGER,         &
     &                   scr_arr_name_length,1,my_MPI_INTEGER,         &
     &                   communicator_glb,ierr)

!     2. collect all names in temporary storage
      allocate(scr_arr_process_name(nr_of_process_glb*100))
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
 

      nr_file_groups = local_counter_file_groups
 
!     final consistency check
      if(nr_file_groups == 1) process_list_glb = 1

      deallocate(scr_arr_process_name)
      deallocate(scr_arr_name_length)

#endif /* VAR_PFS */

!     report the count of intra-groups formed by all processes
      write(print_unit,'(/a,i4,1x,a)')                                    &
      '  *** final output from the group generator:',                     &
      local_counter_file_groups,' intra-node group(s) (is) are built. ***'

  end subroutine set_file_groups
!*******************************************************************************

  subroutine set_communication_levels(group_list,                              &
                                      process_list_glb,                        &
                                      process_list_shared_mem_glb,             &
                                      nr_of_process_glb,                       &
                                      my_process_id_glb,                       &
                                      communicator_glb,                        &
                                      shared_memory_lvl_ijkl,                  &
                                      shared_memory_lvl_cvec,                  &
                                      my_intra_node_id,                        &
                                      my_inter_node_id,                        &
                                      my_shmem_ijkl_id,                        &
                                      my_shmem_cvec_id,                        &
                                      intra_node_master,                       &
                                      shmem_master_ijkl,                       &
                                      shmem_master_cvec,                       &
                                      intra_node_size,                         &
                                      inter_node_size,                         &
                                      shmem_ijkl_size,                         &
                                      shmem_cvec_size,                         &
                                      intra_node_comm,                         &
                                      inter_node_comm,                         &
                                      shmem_ijkl_comm,                         &
                                      shmem_cvec_comm,                         &
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
     integer, intent(out)   :: shmem_master_cvec
     integer, intent(out)   :: intra_node_size
     integer, intent(out)   :: inter_node_size
     integer, intent(out)   :: shmem_ijkl_size
     integer, intent(out)   :: shmem_cvec_size
     integer, intent(out)   :: intra_node_comm
     integer, intent(out)   :: inter_node_comm
     integer, intent(out)   :: shmem_ijkl_comm
     integer, intent(out)   :: shmem_cvec_comm
     integer, intent(out)   :: intra_node_group_id
     integer, intent(out)   :: intra_node_master
!-------------------------------------------------------------------------------
     integer                :: key
     integer                :: color
     integer                :: tmp_group_counter
     integer                :: current_proc
!-------------------------------------------------------------------------------

!     a. intra-node communicator (primarly used as I/O communicator)
 
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

!     c.1. ijkl communicator
 
      key   = my_process_id_glb
      color = process_list_shared_mem_glb(my_process_id_glb+1)
      if(shared_memory_lvl_ijkl == 1) color = 1
 
      call build_new_communicator_group(communicator_glb,        &
                                        shmem_ijkl_comm,         &
                                        shmem_ijkl_size,         &
                                        my_shmem_ijkl_id,        &
                                        color,                   &
                                        key)

!     c.2. cvec communicator
 
      key   = my_process_id_glb
      color = process_list_shared_mem_glb(my_process_id_glb+1)
      if(shared_memory_lvl_cvec == 1) color = 1
 
      call build_new_communicator_group(communicator_glb,        &
                                        shmem_cvec_comm,         &
                                        shmem_cvec_size,         &
                                        my_shmem_cvec_id,        &
                                        color,                   &
                                        key)

!     set new sub-group master (might be master of all CPU's)
      intra_node_master = 0
      shmem_master_ijkl = 0
      shmem_master_cvec = 0

  end subroutine set_communication_levels

!*******************************************************************************

  subroutine build_new_communicator_group(old_communicator,      &
                                          new_communicator,      &
                                          group_size,            &
                                          my_id_group,           &
                                          color,                 &
                                          key)
!*******************************************************************************
!     purpose: split old communicator group into new sub-groups 
!              using key and color. 
!              the new communicator along with the process-ids in the new group 
!              and the group size are returned.
!*******************************************************************************
     integer, intent(in )   :: old_communicator
     integer, intent(in )   :: color
     integer, intent(in )   :: key
     integer, intent(out)   :: group_size
     integer, intent(out)   :: my_id_group
     integer, intent(out)   :: new_communicator
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      call mpi_comm_split(old_communicator,color,key,new_communicator,ierr)
 
!     collect required information about each group:
!       - group size
!       - process id in the group
      call mpi_comm_size(new_communicator, group_size,ierr)
      call mpi_comm_rank(new_communicator,my_id_group,ierr)

  end subroutine build_new_communicator_group
!*******************************************************************************

  subroutine close_communication_model(intra_node_comm,          &
                                       inter_node_comm,          &
                                       shmem_ijkl_comm,          &
                                       shmem_cvec_comm)
!*******************************************************************************
!
!     purpose: reset sub-group communicators.
!
!*******************************************************************************
     integer, intent(inout) :: intra_node_comm
     integer, intent(inout) :: inter_node_comm
     integer, intent(inout) :: shmem_ijkl_comm
     integer, intent(inout) :: shmem_cvec_comm
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

      call mpi_comm_free(intra_node_comm,ierr)
      call mpi_comm_free(inter_node_comm,ierr)
      call mpi_comm_free(shmem_ijkl_comm,ierr)
      call mpi_comm_free(shmem_cvec_comm,ierr)
 
  end subroutine close_communication_model
  
#else 
module dummy_comm_model
#endif
end module
