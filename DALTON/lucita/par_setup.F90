!dalton_copyright_start
!
!
!dalton_copyright_end

module parallel_setup

! stefan: - this module provides all necessary functionality
!           to setup the parallel mcscf/ci model.
!
!           written by sknecht, may 2011 for DALTON LUCITA.

  use communicator_type_module
  use parallel_task_distribution_type_module
  use file_type_module

#ifdef VAR_MPI
  use communication_model
  use file_io_model

#ifdef USE_MPI_MOD_F90
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

  public lucita_setup_parallel_model

#endif
  public lucita_close_parallel_model
  public lucita_reset_parallel_model

  private

  save

contains

#ifdef VAR_MPI

  subroutine lucita_setup_parallel_model(block_list,par_dist_block_list,&
                                         rcctos,nblock,mxciv,nroot)
!******************************************************************************
!
!    purpose:  
!
!*******************************************************************************
#include "parluci.h"
#include "priunit.h"
     integer,    intent(in)                    :: nblock
     integer,    intent(inout)                 :: block_list(nblock)
     integer,    intent(inout)                 :: par_dist_block_list(nblock)
     integer,    intent(inout)                 :: rcctos(nblock)
     integer,    intent(in)                    :: mxciv
     integer,    intent(in)                    :: nroot
!-------------------------------------------------------------------------------
     integer                                   :: file_offset_off
     integer                                   :: my_mpi_comm_world
     integer,                      allocatable :: grouplist_shared_mem(:)
     real(8),                      allocatable :: tmp_block_scaling_fac(:)
!-------------------------------------------------------------------------------
 
      allocate(grouplist_shared_mem(luci_nmproc))
      grouplist_shared_mem = -1

      
!     all important variables (input parameters to module) 
!     are stored on common block /LUPARGROUP/

!     step 1: setup the communication model
!     -------------------------------------
      my_mpi_comm_world = mpi_comm_world ! from MPI_INTEGER_KIND to default integer kind
      call setup_communication_model(nflgrps,                         &
                                     iiomod,                          &
                                     shared_m,                        &
                                     it_shl,                          &
                                     ic_shl,                          &
                                     grouplist_shared_mem,            &
                                     luci_myproc,                     &
                                     luci_nmproc,                     &
                                     my_mpi_comm_world,               &
                                     mynew_id,                        &
                                     icomm_id,                        &
                                     mynew_id_sm,                     &
                                     mynew_id_sm_c,                   &
                                     newcomm_proc,                    &
                                     icomm_size,                      &
                                     newcomm_proc_sm,                 &
                                     newcomm_proc_sm_c,               &
                                     mynew_comm,                      &
                                     icomm,                           &
                                     mynew_comm_sm,                   &
                                     mynew_comm_sm_c,                 &
                                     my_groupn,                       &
                                     n_master,                        &
                                     n_master_sm,                     &
                                     n_master_sm_c,                   &
                                     lupri)
     
  
!     free currently unused shared-memory group list
      deallocate(grouplist_shared_mem)

!     step 2: setup static load balancing
!     -----------------------------------
!     block_list          == list of all blocks containing their length
!     par_dist_block_list == list of all blocks containing their assigned CPU
!     rcctos              == list of bvec <--> sigma connections

      if(.not.ptask_distribution%parallel_task_distribution_set)then
        allocate(tmp_block_scaling_fac(nblock))
        call block_distr_drv(nblock,block_list,par_dist_block_list,            &
                             rcctos,tmp_block_scaling_fac,l_combi,             &
                             communicator_info%total_process_list)
        ptask_distribution%parallel_task_distribution_set = .true.
        deallocate(tmp_block_scaling_fac)
      end if

!     step 3: organize MPI file I/O offsets with the following ordering: 
!     ------------------------------------------------------------------
!     iref, ilu1, ilu2, idia, iluc, ilu[3-7]

      if(.not.file_info%file_type_init)then
        call file_init_lucipar(file_info, mxciv, nroot)
        file_info%file_type_init = .true.
      else 
        return
      end if

      file_offset_off     = 0

      call set_file_io_offset(file_info%active_nr_f_lucipar,                                                 &
                              file_offset_off,                                                               &
                              file_info%file_offsets,                                                        &
                              file_info%facofffl,                                                            &
                              my_vec1_ioff,                                                                  &
                              my_vec2_ioff,                                                                  &
                              my_act_blk1,                                                                   &
                              my_act_blk2,                                                                   &
                              my_act_blk_all,                                                                &
                              .false.,                                                                       &
                              luci_myproc,                                                                   &
                              num_blocks2,                                                                   &
                              newcomm_proc,                                                                  &
                              communicator_info%intra_node_group_list,                                       &
                              par_dist_block_list,                                                           &
                              block_list)

!     save output in common block variables
      my_lur_off = file_info%file_offsets( 1)
      my_lu1_off = file_info%file_offsets( 2)
      my_lu2_off = file_info%file_offsets( 3)
      my_dia_off = file_info%file_offsets( 4)
      my_luc_off = file_info%file_offsets( 5)

      call file_set_list_lucipar(file_info, my_act_blk2, num_blocks2)

!     length for allocation of file arrays
      iall_lur   = file_info%max_list_length
      iall_lu1   = file_info%max_list_length
      iall_lu2   = file_info%max_list_length
      iall_luc   = file_info%max_list_length_bvec

!     step 5: setup the file i/o model
!     --------------------------------
      call setup_file_io_model(mynew_comm,                           &
                               file_info%active_nr_f_lucipar,        &
                               file_info%fh_lu,                      &
                               0,                                    &
                               my_groupn,                            &
                               newcomm_proc,                         &
                               'parci',                              &
                               lupri)
  
!     transfer file handles to common block /LUCIAPFILE/ (in parluci.h)
      ILUR = file_info%fh_lu( 1)
      ILU1 = file_info%fh_lu( 2)
      ILU2 = file_info%fh_lu( 3)
      IDIA = file_info%fh_lu( 4)
      ILUC = file_info%fh_lu( 5)

!     safe information for files used exclusively in parallel LUCITA-CI runs
      if(file_info%active_nr_f_lucipar > 5)then
        ILU3       = file_info%fh_lu( 6)
        ILU4       = file_info%fh_lu( 7)
        ILU5       = file_info%fh_lu( 8)
        ILU6       = file_info%fh_lu( 9)
        ILU7       = file_info%fh_lu(10)

        my_lu3_off = file_info%file_offsets( 6)
        my_lu4_off = file_info%file_offsets( 7)
        my_lu5_off = file_info%file_offsets( 8)
        my_lu6_off = file_info%file_offsets( 9)
        my_lu7_off = file_info%file_offsets(10)

        iall_lu3   = file_info%max_list_length
        iall_lu4   = file_info%max_list_length
        iall_lu5   = file_info%max_list_length
        iall_lu6   = file_info%max_list_length
        iall_lu7   = file_info%max_list_length
      end if

  end subroutine lucita_setup_parallel_model
!*******************************************************************************
#endif /* VAR_MPI */

  subroutine lucita_close_parallel_model(reset_communication)
!******************************************************************************
!
!    purpose: close down the parallel model
!         --> close parallel files and nullify file pointers
!         --> reset parallel communication groups
!
!*******************************************************************************
     logical, intent(in)               :: reset_communication
!-------------------------------------------------------------------------------
     integer                           :: fh_offset
!-------------------------------------------------------------------------------

!     step 1: close all open parallel files
#ifdef VAR_MPI
      fh_offset = 0
      if(file_info%file_type_init)then
        call close_file_io_model(file_info%active_nr_f_lucipar,               &
                                 fh_offset,                                   &
                                 file_info%fh_lu)
      end if
#endif

!     step 2: reset file lists and initialized variables
      call file_free_lucipar(file_info)

!     step 3: free communication handle(s)
      if(reset_communication) call communicator_free_lucipar(communicator_info)

  end subroutine lucita_close_parallel_model
!*******************************************************************************

  subroutine lucita_reset_parallel_model(block_list,par_dist_block_list,               &
                                         rcctos,nblock,mxciv,nroot,                    &
                                         reset_communication_model,set_mc_file_model)
!*******************************************************************************
!
!    purpose: close down the parallel model
!         --> close parallel files and nullify file pointers
!         --> reset parallel communication groups
!
!*******************************************************************************
     integer, intent(in)    :: nblock
     integer, intent(inout) :: block_list(nblock)
     integer, intent(inout) :: par_dist_block_list(nblock)
     integer, intent(inout) :: rcctos(nblock)
     integer, intent(in)    :: mxciv
     integer, intent(in)    :: nroot
     logical, intent(in)    :: reset_communication_model
     logical, intent(in)    :: set_mc_file_model
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

!     step 1: close parallel file (and possibly communication model)
      call lucita_close_parallel_model(reset_communication_model)

!     step 2: possibly switch to MCSCF file model
      file_info%file_type_mc = set_mc_file_model

#ifdef VAR_MPI
!     step 3: re-open file and communication model in parallel runs
      call lucita_setup_parallel_model(block_list,par_dist_block_list,                 &
                                       rcctos,nblock,mxciv,nroot)
#endif

  end subroutine lucita_reset_parallel_model
!*******************************************************************************
end module
