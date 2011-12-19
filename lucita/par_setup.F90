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
  use communication_model
  use file_io_model
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
                                         nblock,                        &
                                         rcctos,kiluclist,              &
                                         kilu1list,kilu2list,kilu3list, &
                                         kilu4list,kilu5list,kilu6list, &
                                         kilu7list,fh_array,mxciv,nroot)
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
     integer   , intent(inout)                 :: rcctos(nblock)
     integer(8), intent(inout)                 :: kiluclist, kilu1list, kilu2list, kilu3list
     integer(8), intent(inout)                 :: kilu4list, kilu5list, kilu6list, kilu7list
     integer,    intent(inout)                 :: fh_array(nr_files)
     integer,    intent(in)                    :: mxciv
     integer,    intent(in)                    :: nroot
!-------------------------------------------------------------------------------
     integer                                   :: file_offset_off
     integer,                      allocatable :: grouplist_shared_mem(:)
     real(8),                      allocatable :: tmp_block_scaling_fac(:)
     integer,                      allocatable :: file_offset_fac_i4(:)
     integer(KIND=MPI_OFFSET_KIND),allocatable :: file_offset_fac(:)
     integer(KIND=MPI_OFFSET_KIND),allocatable :: file_offset_array(:)
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
     
!     step 2: setup the file i/o model
!     --------------------------------
      call  setup_file_io_model(mynew_comm,                           &
                                nr_files,                             &
                                fh_array,                             &
                                my_groupn,                            &
                                newcomm_proc,                         &
                                'parci',                              &
                                lupri)
  
!     transfer file handles to common block /LUCIAPFILE/ (in parluci.h)
!     ILU1 has to be the last in the row because it is supposed to be
!     closed after we have copied the final solution vector back to
!     sequential format. Before doing that we want to free disk space by
!     deleting all other MPI-I/O files.
      IDIA = fh_array(1)
      ILUC = fh_array(2)
      ILU7 = fh_array(3)
      ILU6 = fh_array(4)
      ILU5 = fh_array(5)
      ILU4 = fh_array(6)
      ILU3 = fh_array(7)
      ILU2 = fh_array(8)
      ILU1 = fh_array(9)
  
!     free currently unused shared-memory group list
      deallocate(grouplist_shared_mem)

!     step 3: setup static load balancing
!     -----------------------------------
!     block_list          == list of all blocks containing their length
!     par_dist_block_list == list of all blocks containing their assigned CPU

      allocate(tmp_block_scaling_fac(nblock))
      call block_distr_drv(nblock,block_list,par_dist_block_list,           &
                           rcctos,tmp_block_scaling_fac,                    &
                           l_combi,communicator_info%total_process_list)
      deallocate(tmp_block_scaling_fac)

!     step 4: organize MPI file I/O offsets with the following ordering: 
!     ------------------------------------------------------------------
!     idia, iluc, ilu[2-7], ilu1

      allocate(file_offset_fac_i4(nr_files))
      allocate(file_offset_fac(nr_files))
      allocate(file_offset_array(nr_files))
      file_offset_array = 0
      file_offset_fac   = 0
      file_offset_off   = 0

      file_offset_fac_i4(1) = 1
      file_offset_fac_i4(2) = 0
      file_offset_fac_i4(3) = mxciv + nroot
      file_offset_fac_i4(4) = mxciv + nroot
      file_offset_fac_i4(5) = mxciv + nroot
      file_offset_fac_i4(6) = mxciv + nroot
      file_offset_fac_i4(7) = 1
      file_offset_fac_i4(8) = mxciv
      file_offset_fac_i4(9) = nroot

      file_offset_fac(1) = file_offset_fac_i4(1)
      file_offset_fac(2) = file_offset_fac_i4(2)
      file_offset_fac(3) = file_offset_fac_i4(3)
      file_offset_fac(4) = file_offset_fac_i4(4)
      file_offset_fac(5) = file_offset_fac_i4(5)
      file_offset_fac(6) = file_offset_fac_i4(6)
      file_offset_fac(7) = file_offset_fac_i4(7)
      file_offset_fac(8) = file_offset_fac_i4(8)
      file_offset_fac(9) = file_offset_fac_i4(9)
      
      call set_file_io_offset(nr_files,                                 &
                              file_offset_off,                          &
                              file_offset_array,                        &
                              file_offset_fac,                          &
                              my_vec1_ioff,                             &
                              my_vec2_ioff,                             &
                              my_act_blk1,                              &
                              my_act_blk2,                              &
                              my_act_blk_all,                           &
                              .false.,                                  &
                              luci_myproc,                              &
                              num_blocks2,                              &
                              newcomm_proc,                             &
                              communicator_info%intra_node_group_list,  &
                              par_dist_block_list,                      &
                              block_list)

!     save output in common block variables
      my_dia_off = file_offset_array(1)
      my_luc_off = file_offset_array(2)
      my_lu2_off = file_offset_array(3)
      my_lu3_off = file_offset_array(4)
      my_lu4_off = file_offset_array(5)
      my_lu5_off = file_offset_array(6)
      my_lu6_off = file_offset_array(7)
      my_lu7_off = file_offset_array(8)
      my_lu1_off = file_offset_array(9)

      deallocate(file_offset_array)
      deallocate(file_offset_fac)

!     length for allocation of file arrays
      iall_luc =         1             * num_blocks2
      iall_lu2 = file_offset_fac_i4(3) * my_act_blk2
      iall_lu3 = file_offset_fac_i4(4) * my_act_blk2
      iall_lu4 = file_offset_fac_i4(5) * my_act_blk2
      iall_lu5 = file_offset_fac_i4(6) * my_act_blk2
      iall_lu6 = file_offset_fac_i4(7) * my_act_blk2
      iall_lu7 = file_offset_fac_i4(8) * my_act_blk2
      iall_lu1 = file_offset_fac_i4(9) * my_act_blk2

      deallocate(file_offset_fac_i4)

!     step 5: allocate file arrays - return pointers to calling subroutine
!     --------------------------------------------------------------------
      call memman(kiluclist,iall_luc,'ADDS  ',1,'LUCLST')
      call memman(kilu1list,iall_lu1,'ADDS  ',1,'LU1LST')
      call memman(kilu2list,iall_lu2,'ADDS  ',1,'LU2LST')
      call memman(kilu3list,iall_lu3,'ADDS  ',1,'LU3LST')
      call memman(kilu4list,iall_lu4,'ADDS  ',1,'LU4LST')
      call memman(kilu5list,iall_lu5,'ADDS  ',1,'LU5LST')
      call memman(kilu6list,iall_lu6,'ADDS  ',1,'LU6LST')
      call memman(kilu7list,iall_lu7,'ADDS  ',1,'LU7LST')

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
