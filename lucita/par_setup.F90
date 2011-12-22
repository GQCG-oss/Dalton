!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module parallel_setup

! stefan: - this module provides all necessary functionality
!           to setup the parallel mcscf/ci model.
!
!           written by sknecht, may 2011 for DALTON LUCITA.

  use communicator_type_module
  use parallel_task_distribution_type_module
  use communication_model
  use file_io_model
  use file_type_module
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

  public lucita_setup_parallel_model
  public lucita_close_parallel_model

  private

  save

contains

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
     integer,                      allocatable :: grouplist_shared_mem(:)
     real(8),                      allocatable :: tmp_block_scaling_fac(:)
!-------------------------------------------------------------------------------
 
      allocate(grouplist_shared_mem(luci_nmproc))
      grouplist_shared_mem = -1

      
!     all important variables (input parameters to module) 
!     are stored on common block /LUPARGROUP/

!     step 1: setup the communication model
!     -------------------------------------
      call setup_communication_model(nflgrps,                         &
                                     iiomod,                          &
                                     shared_m,                        &
                                     it_shl,                          &
                                     ic_shl,                          &
                                     grouplist_shared_mem,            &
                                     luci_myproc,                     &
                                     luci_nmproc,                     &
                                     mpi_comm_world,                  &
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
!     idia, iluc, ilu[2-7], ilu1

      if(.not.file_info%file_type_init)then
        call file_init_lucipar(file_info, mxciv, nroot)
        file_info%file_type_init = .true.
      end if

      file_offset_off     = 0

      call set_file_io_offset(nr_files,                                                                      &
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
      my_lu3_off = file_info%file_offsets( 6)
      my_lu4_off = file_info%file_offsets( 7)
      my_lu5_off = file_info%file_offsets( 8)
      my_lu6_off = file_info%file_offsets( 9)
      my_lu7_off = file_info%file_offsets(10)

      call file_set_list_lucipar(file_info, my_act_blk2, num_blocks2)

!     length for allocation of file arrays
      iall_lur   = file_info%max_list_length
      iall_lu1   = file_info%max_list_length
      iall_lu2   = file_info%max_list_length
      iall_lu3   = file_info%max_list_length
      iall_lu4   = file_info%max_list_length
      iall_lu5   = file_info%max_list_length
      iall_lu6   = file_info%max_list_length
      iall_lu7   = file_info%max_list_length
      iall_luc   = file_info%max_list_length_bvec

!     step 5: setup the file i/o model
!     --------------------------------
      call setup_file_io_model(mynew_comm,                           &
                               nr_files,                             &
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
      ILU3 = file_info%fh_lu( 6)
      ILU4 = file_info%fh_lu( 7)
      ILU5 = file_info%fh_lu( 8)
      ILU6 = file_info%fh_lu( 9)
      ILU7 = file_info%fh_lu(10)

  end subroutine lucita_setup_parallel_model
!*******************************************************************************

  subroutine lucita_close_parallel_model(nr_files,lucita_ci_run_id,         &
                                         fh_array)
!******************************************************************************
!
!    purpose: close down the parallel model
!         --> close parallel files and nullify file pointers
!         --> reset parallel communication groups
!
!*******************************************************************************
     integer,            intent(in)    :: nr_files
     character (len=72), intent(in)    :: lucita_ci_run_id
     integer,            intent(inout) :: fh_array(nr_files)
!-------------------------------------------------------------------------------
     integer                           :: files_to_close
     integer                           :: fh_offset
!-------------------------------------------------------------------------------

!     step 1: close all remaining open parallel files
      files_to_close = nr_files
      fh_offset      = 0
      call close_file_io_model(files_to_close,                              &
                               fh_offset,                                   &
                               fh_array)

  end subroutine lucita_close_parallel_model
!*******************************************************************************
#else
module dummy_parallel_setup
#endif
end module
