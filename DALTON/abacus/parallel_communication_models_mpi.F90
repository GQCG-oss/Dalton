!dalton_copyright_start
!
!
!dalton_copyright_end
!
! stefan jan 2014: this module contains the type definition for the
!                  MPI communication handles to be potentially used in modules.
!                  currently in use in Dalton: 
!                                    - lucita (after merge with hjaa-srdft)
!                                    - mcscf  (after merge with hjaa-srdft)
!                                    - hermit
!
module parallel_communication_models_mpi

#ifdef VAR_MPI
#ifdef USE_MPI_MOD_F90
  use mpi
#endif
#endif

  implicit none

  public communication_type_mpi
  public communication_init_mpi
  public communication_free_mpi


! communication type definition
! ----------------------------------------------------------------------------
  type communication_type_mpi

    integer ::                   &
      my_intra_node_id,          &             ! intra-node ID
      my_inter_node_id,          &             ! inter-node ID
      my_shmem_node_id,          &             ! shared-memory-node ID
      intra_node_size,           &             ! size of intra-node group
      inter_node_size,           &             ! size of inter-node group
      shmem_node_size,           &             ! size of shared-memory group
      intra_node_group_id,       &             ! intra-node group ID
      communication_intranode,   &             ! intra-node communication handle
      communication_internode,   &             ! inter-node communication handle
      communication_shmemnode,   &             ! shared-memory communication handle: integrals, fock matrices, density matrices, etc
      communication_glb_world                  ! global communicator handle

    logical ::                   &
      communication_type_init = .false.        ! status of the communication type

    integer, allocatable ::      &
      total_process_list(:),     &             ! list of all processes with their process ID
      shared_mem_group_list(:),  &             ! list of all processes within a shared-memory group and their process ID
      intra_node_group_list(:)                 ! list of all processes within an intra-node group and their process ID

  end type communication_type_mpi

! communication type object
  type(communication_type_mpi), public, save :: communication_info_mpi
! ----------------------------------------------------------------------------

  private

#ifdef VAR_MPI
#ifndef USE_MPI_MOD_F90
#include "mpif.h"
#endif
  integer(kind=MPI_INTEGER_KIND)         :: ierr_mpi
  integer(kind=MPI_INTEGER_KIND)         :: istat(MPI_STATUS_SIZE)
#endif

#if defined (VAR_INT64)
#define my_MPI_INTEGER MPI_INTEGER8
#else
#define my_MPI_INTEGER MPI_INTEGER4
#endif

contains 

!*******************************************************************************
  subroutine communication_init_mpi(A,                     & 
                                    my_process_id_glb,     &
                                    nr_of_process_glb,     &
                                    communicator_glb)

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
    type(communication_type_mpi) :: A
    integer, intent(in)          :: nr_of_process_glb
    integer, intent(in)          :: communicator_glb
    integer, intent(in)          :: my_process_id_glb
!   ----------------------------------------------------------------------------
!   ----------------------------------------------------------------------------

!   reset old type information
    call communication_free_mpi(A,nr_of_process_glb)

    A%communication_type_init  = .true.
    A%communication_intranode  = -1
    A%communication_internode  = -1
    A%communication_shmemnode  = -1
    A%communication_glb_world  = -1
    A%my_intra_node_id         = -1
    A%my_inter_node_id         = -1
    A%my_shmem_node_id         = -1
    A%intra_node_size          = -1
    A%inter_node_size          = -1
    A%shmem_node_size          = -1
    A%intra_node_group_id      = -1

    allocate(A%total_process_list(nr_of_process_glb))
    allocate(A%shared_mem_group_list(nr_of_process_glb))
    allocate(A%intra_node_group_list(nr_of_process_glb))
    A%total_process_list      = -1
    A%shared_mem_group_list   = -1
    A%intra_node_group_list   = -1

    if(nr_of_process_glb > 1)then
#ifdef VAR_MPI
!     start communication model setup
!     -------------------------------

!     1. determine # communication process groups and store group-IDs in total_process_list
      call set_communication_groups(A%total_process_list,     &
                                    my_process_id_glb,        &
                                    nr_of_process_glb,        &
                                    communicator_glb)

!     2. setup communicators and process-id for the various communication levels
!        a. intra-node
!        b. inter-node
!        c. shared-memory
      call set_communication_levels(A%intra_node_group_list,                   &
                                    A%shared_mem_group_list,                   &
                                    A%total_process_list,                      &
                                    nr_of_process_glb,                         &
                                    my_process_id_glb,                         &
                                    communicator_glb,                          &
                                    A%my_intra_node_id,                        &
                                    A%my_inter_node_id,                        &
                                    A%my_shmem_node_id,                        &
                                    A%intra_node_size,                         &
                                    A%inter_node_size,                         &
                                    A%shmem_node_size,                         &
                                    A%communication_intranode,                 &
                                    A%communication_internode,                 &
                                    A%communication_shmemnode,                 &
                                    A%communication_glb_world,                 &
                                    A%intra_node_group_id)                      
#endif
    end if

  end subroutine communication_init_mpi
!*******************************************************************************

  subroutine communication_free_mpi(A,nr_of_process_glb)

!   ----------------------------------------------------------------------------
    type(communication_type_mpi) :: A
    integer, intent(in)          :: nr_of_process_glb
!   ----------------------------------------------------------------------------
#ifdef VAR_MPI
    integer(kind=MPI_INTEGER_KIND) :: comm_mpi
#endif
!   ----------------------------------------------------------------------------

    if(.not. A%communication_type_init) return

    A%communication_type_init = .false.
    deallocate(A%total_process_list)
    deallocate(A%shared_mem_group_list)
    deallocate(A%intra_node_group_list)

    if(nr_of_process_glb > 1)then
#ifdef VAR_MPI
      comm_mpi = A%communication_intranode
      call mpi_comm_free(comm_mpi,  ierr_mpi)
      comm_mpi = A%communication_internode
      call mpi_comm_free(comm_mpi,  ierr_mpi)
      comm_mpi = A%communication_shmemnode
      call mpi_comm_free(comm_mpi,  ierr_mpi)
#endif
    end if

  end subroutine communication_free_mpi
!*******************************************************************************
#ifdef VAR_MPI

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
     integer                          :: finished_loop  
     character (len=255)              :: process_name
     character (len=255)              :: scr_process_name
     character (len=255), allocatable :: scr_arr_process_name(:)
     integer, allocatable             :: scr_arr_name_length(:)
!-------------------------------------------------------------------------------
     integer(kind=MPI_INTEGER_KIND)   :: proc_id, len_mpi, comm_glb_mpi, ierr_mpi
!-------------------------------------------------------------------------------

      comm_glb_mpi = communicator_glb
      process_list_glb = -1

!     find system-dependent unique process name
      call mpixprocname(process_name,process_name_length)

      allocate(scr_arr_name_length(nr_of_process_glb))
      scr_arr_name_length = 0

!     1. gather all name length
      len_mpi = 1
      call mpi_allgather(process_name_length,len_mpi,my_MPI_INTEGER,   &
     &                   scr_arr_name_length,len_mpi,my_MPI_INTEGER,   &
     &                   comm_glb_mpi,ierr_mpi)

!     2. collect all names in temporary storage
      allocate(scr_arr_process_name(nr_of_process_glb*255))
      scr_arr_process_name(my_process_id_glb+1) = process_name

      do proc_id = 0, nr_of_process_glb-1

         scr_process_name = process_name
         len_mpi          = scr_arr_name_length(proc_id+1)
         call mpi_bcast(scr_process_name,len_mpi, MPI_CHARACTER,       &
                        proc_id, comm_glb_mpi, ierr_mpi)

         if(my_process_id_glb /= proc_id)then
           scr_arr_process_name(proc_id+1) = scr_process_name(1:scr_arr_name_length(proc_id+1))
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
        write(*,'(/a)')                                                     &
        ' Output from the communication group generator:'
        if(local_counter_file_groups == 1)then
          write(*,'(i4,a/)')                                                &
          local_counter_file_groups,' intra-node group has been built.'
        else
          write(*,'(i4,a/)')                                                &
          local_counter_file_groups,' intra-node groups have been built.'
        end if
      end if

  end subroutine set_communication_groups
!*******************************************************************************

  subroutine set_communication_levels(group_list,                              &
                                      shared_mem_list,                         &
                                      process_list_glb,                        &
                                      nr_of_process_glb,                       &
                                      my_process_id_glb,                       &
                                      communicator_glb,                        &
                                      my_intra_node_id,                        &
                                      my_inter_node_id,                        &
                                      my_shmem_node_id,                        &
                                      intra_node_size,                         &
                                      inter_node_size,                         &
                                      shmem_node_size,                         &
                                      intra_node_comm,                         &
                                      inter_node_comm,                         &
                                      shmem_ijkl_comm,                         &
                                      total_area_comm,                         &
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
     integer, intent(out)   :: group_list(nr_of_process_glb)
! May 2016 hjaaj: removing intent on shared_mem_list, because it throws intel compile
!    warnings as shared_mem_list is not set in this routine (currently?)
!    integer, intent(out)   :: shared_mem_list(nr_of_process_glb)
     integer                :: shared_mem_list(nr_of_process_glb)
     integer, intent(out)   :: my_intra_node_id
     integer, intent(out)   :: my_inter_node_id
     integer, intent(out)   :: my_shmem_node_id
     integer, intent(out)   :: intra_node_size
     integer, intent(out)   :: inter_node_size
     integer, intent(out)   :: shmem_node_size
     integer, intent(out)   :: intra_node_comm
     integer, intent(out)   :: inter_node_comm
     integer, intent(out)   :: shmem_ijkl_comm
     integer, intent(out)   :: total_area_comm
     integer, intent(out)   :: intra_node_group_id
!-------------------------------------------------------------------------------
     integer                :: key
     integer                :: color
     integer                :: tmp_group_counter
     integer                :: current_proc
     integer                :: numa_procs
     integer                :: numa_counter
     integer                :: numa_nodes
     integer                :: i
     integer, allocatable   :: tmp_array(:)
     character(len=6)       :: numa_procs_env
!-------------------------------------------------------------------------------

!     0. global communicator
      total_area_comm = communicator_glb

!     a. intra-node communicator
 
      key   = my_process_id_glb
      color = process_list_glb(my_process_id_glb+1)

      call build_new_communication_group(communicator_glb,        &
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

!     b. shared memory communicator
      allocate(tmp_array(nr_of_process_glb))
      tmp_array(1:nr_of_process_glb) = process_list_glb(1:nr_of_process_glb) 
      key                            = my_process_id_glb
      color                          = tmp_array(my_process_id_glb+1)

!
!     IMPORTANT ADVICE from Peter xxx at the NSC Linkoeping - use e.g. hwloc to get the NUMA node master rather than the
!     intra-node master: NUMA ==> all processes with close memory == highest performance
!     todo: automatic scheme w/o required user intervention  -- sknecht feb 2013


!     !> NUMA node mode...
      numa_procs          = 0
      numa_procs_env(1:6) = '      '
      call getenv('NUMA_PROCS',numa_procs_env)
      read(numa_procs_env, '(i6)') numa_procs

!     if( numa_procs > 0)then
!       write(*,*) 'NUMA process distribution active for shared memory'
!       write(*,*) '# processes per NUMA node (user defined) ==>      ',numa_procs
!     end if

      if( numa_procs > 0)then

!       OpenMPIs --bysocket --bind-to-socket policy: round-robin fashion between the X sockets (often socket == NUMA node) 
!       thus, we will follow this strategy here...
        numa_nodes = intra_node_size/numa_procs
        if(mod(intra_node_size,numa_procs) /= 0 ) write(*,*) ' ** Warning: asymmetric NUMA node allocation'

        numa_counter = 0
        i            = 0
        do
          i = i + 1
          if( i             > nr_of_process_glb )exit
          if( tmp_array(i) ==             color )then
            numa_counter    = numa_counter + 1
            if( (i-1) == my_process_id_glb )then
              color         = color        + 10000*numa_counter ! offset for next NUMA node color: color + 10000*numa_counter
              tmp_array(i)  = color; exit
            end if
!           reset 'round-robin' counter
            if(numa_counter == numa_nodes) numa_counter = 0
          end if
        end do
!       write(*,*) 'pid, color, numa_counter, numa_nodes',my_process_id_glb,color, numa_counter, numa_nodes
      end if
!     !> end

      call build_new_communication_group(communicator_glb,        &
                                         shmem_ijkl_comm,         &
                                         shmem_node_size,         &
                                         my_shmem_node_id,        &
                                         color,                   &
                                         key)

      deallocate(tmp_array)

!     c. inter-node resp. inter-NUMA node communicator

      key   = process_list_glb(my_process_id_glb+1)
      color = 2
      if(my_shmem_node_id /= 0) color = 3
 
      call build_new_communication_group(communicator_glb,        &
                                         inter_node_comm,         &
                                         inter_node_size,         &
                                         my_inter_node_id,        &
                                         color,                   &
                                         key)
!     write(*,*) 'pid, # my_inter_node_id,my_shmem_node_id,inter_node_size',my_process_id_glb,&
!                        my_inter_node_id,my_shmem_node_id,inter_node_size

  end subroutine set_communication_levels

!*******************************************************************************

  subroutine build_new_communication_group(old_communication,     &
                                           new_communication,     &
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
     integer(kind=MPI_INTEGER_KIND)   :: old_comm_mpi
     integer(kind=MPI_INTEGER_KIND)   :: color_mpi
     integer(kind=MPI_INTEGER_KIND)   :: key_mpi
     integer(kind=MPI_INTEGER_KIND)   :: group_size_mpi
     integer(kind=MPI_INTEGER_KIND)   :: my_id_group_mpi
     integer(kind=MPI_INTEGER_KIND)   :: new_comm_mpi, ierr_mpi
!-------------------------------------------------------------------------------

      old_comm_mpi = old_communication
      color_mpi    = color
      key_mpi      = key

      call mpi_comm_split(old_comm_mpi,color_mpi,key_mpi,new_comm_mpi,ierr_mpi)
 
!     collect required information about each group:
!       - group size
!       - process id in the group
      call mpi_comm_size(new_comm_mpi,group_size_mpi,ierr_mpi)
      call mpi_comm_rank(new_comm_mpi,my_id_group_mpi,ierr_mpi)

      group_size        = group_size_mpi
      my_id_group       = my_id_group_mpi
      new_communication = new_comm_mpi

  end subroutine build_new_communication_group
#endif

end module parallel_communication_models_mpi
