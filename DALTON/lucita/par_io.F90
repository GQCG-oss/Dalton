!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
module par_mcci_io

! stefan: - this module provides all necessary functionality
!           to perform disk i/o operations in parallel mcscf/ci 
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

  public mcci_cp_vcd_batch

  private

  save

  integer, private                       :: istat(MPI_STATUS_SIZE)
  integer, private                       :: ierr

contains 

  subroutine mcci_cp_vcd_batch(luin,             &
                               luout,            &
                               xmat,             &
                               nbatch,           &
                               lbatch,           &
                               lebatch,          &
                               i1batch,          &
                               ibatch,           &
                               block_info,       &
                               my_ioff_luin,     &
                               my_ioff_luout,    &
                               luinlist,         &
                               luoutlist,        &
                               cs_fullvector_off,&
                               cs_fullactblk_off,&
                               joff)
!******************************************************************************
!
!    purpose:  copy c/sigma-vector from file LUIN to LUOUT batchwise:
!                - update the file lists
!                - active blocks on the MPI-files are flagged 
!                  by a nonzero length
!
!              note: joff = (ivec resp. iroot) - 1
!
!******************************************************************************
     real(8), intent(inout) :: xmat(*)
     integer, intent(in)    :: luin
!    integer, intent(in)    :: luout ! stefan: open mpi 1.5 with f90-bindings complains here
     integer                :: luout
     integer, intent(in)    :: nbatch
     integer, intent(in)    :: lbatch(*)
     integer, intent(in)    :: lebatch(*)
     integer, intent(in)    :: i1batch(*)
     integer, intent(in)    :: ibatch(8,*)
     integer, intent(in)    :: block_info(*)
     integer, intent(in)    :: luinlist(*)
     integer, intent(inout) :: luoutlist(*)
     integer, intent(in)    :: joff
!    individual offsets for each process
     integer(kind=MPI_OFFSET_KIND), intent(in)  :: my_ioff_luin
     integer(kind=MPI_OFFSET_KIND), intent(in)  :: my_ioff_luout
     integer(kind=MPI_OFFSET_KIND), intent(in)  :: cs_fullvector_off
     integer                      , intent(in)  :: cs_fullactblk_off
!------------------------------------------------------------------------------
     integer(kind=MPI_OFFSET_KIND)  :: ioffset_in  
     integer(kind=MPI_OFFSET_KIND)  :: ioffset_out
     integer(kind=MPI_OFFSET_KIND)  :: ioffset_scratch
     integer(kind=MPI_OFFSET_KIND)  :: my_ioffset_scratch
     integer(kind=MPI_OFFSET_KIND)  :: iscratch_special
     integer                        :: isbatch
     integer                        :: ioffset_int_in
     integer                        :: ioffset_int_out
     integer                        :: num_blk
     integer                        :: num_blk_cnt
     integer                        :: num_blk_cnt_act
!******************************************************************************
!
!     initialize scratch offsets
      my_ioffset_scratch = 0
      ioffset_scratch    = 0
      ioffset_in         = 0
      ioffset_out        = 0
      iscratch_special   = 0
      num_blk            = 0
      num_blk_cnt        = 0
      ioffset_int_in     = 0
      ioffset_int_out    = 0
      num_blk_cnt_act    = 0

!     transfer data from file "luin" batchwise to file "luout"
!     --------------------------------------------------------

      do isbatch = 1, nbatch

        call dzero(xmat,lebatch(isbatch))

!       set offset for data read operation
        ioffset_in = my_ioff_luin + (cs_fullvector_off * joff) + my_ioffset_scratch
        ioffset_int_in = (cs_fullactblk_off * joff) + num_blk_cnt + 1
 
        num_blk_cnt_act  = 0
        ioffset_int_out  = num_blk + 1

!       read batch 
        call read_batch_pario(xmat,                              &
                              luin,                              &
                              luinlist,                          &
                              luoutlist,                         &
                              lbatch(isbatch),                   &
                              ibatch(1,i1batch(isbatch)),        &
                              ioffset_in,                        &
                              ioffset_int_in,                    &
                              ioffset_int_out,                   &
                              num_blk_cnt_act)
 
!       set offset for data write operation
        ioffset_out = my_ioff_luout + ioffset_scratch
     
        iscratch_special = 0
 
        call write_batch_pario(xmat,                             &
                               luout,                            &
                               luoutlist,                        &
                               lbatch(isbatch),                  &
                               ibatch(1,i1batch(isbatch)),       &
                               block_info,                       &
                               ioffset_out,                      &
                               ioffset_int_out,                  &
                               iscratch_special)                  
 
!       keep track of correct offset
!       a. output
        ioffset_scratch    = ioffset_scratch + iscratch_special
        num_blk            = num_blk + lbatch(isbatch)
!       b. input
        num_blk_cnt        = num_blk_cnt + num_blk_cnt_act
        my_ioffset_scratch = my_ioffset_scratch + lebatch(isbatch)
 
      end do

  end subroutine mcci_cp_vcd_batch

!******************************************************************************
  subroutine read_batch_pario(xmat,                &
                              luin,                &
                              luinlist,            &
                              luoutlist,           &
                              num_blocks_in_batch, &
                              batch_info,          &
                              ioff_luin,           &
                              luinlist_offset,     &
                              luoutlist_offset,    &
                              num_blk_cnt_act)
!******************************************************************************
!
!    purpose:  read batch of c/sigma-vector from file LUIN:
!                - update the file lists
!                - active blocks on the MPI-files are flagged 
!                  by a nonzero length
!
!******************************************************************************
     real(8), intent(inout) :: xmat(*)
     integer, intent(in)    :: luin
     integer, intent(in)    :: luinlist(*)
     integer, intent(inout) :: luoutlist(*)
     integer, intent(in)    :: num_blocks_in_batch
     integer, intent(in)    :: batch_info(8,*)
     integer, intent(inout) :: num_blk_cnt_act
!    offset
     integer(kind=MPI_OFFSET_KIND), intent(inout) :: ioff_luin
     integer, intent(in)                          :: luinlist_offset
     integer, intent(in)                          :: luoutlist_offset
!------------------------------------------------------------------------------
     integer(kind=MPI_OFFSET_KIND)  :: blk_len
     integer                        :: mem_off
     integer                        :: is_blk
!******************************************************************************

      do is_blk = 1, num_blocks_in_batch

        blk_len = batch_info(8,is_blk)

        if(blk_len > 0)then

!         count the active blocks
          num_blk_cnt_act = num_blk_cnt_act + 1

          if(luinlist(luinlist_offset+num_blk_cnt_act-1) > 0)then

!           write length into file array for output file
            luoutlist(luoutlist_offset+is_blk-1) = batch_info(8,is_blk)

!           set memory offset
            mem_off = batch_info(6,is_blk)

!           read vector
            call mpi_file_read_at(luin,                    &
                                  ioff_luin,               &
                                  xmat(mem_off),           &
                                  batch_info(8,is_blk),    &
                                  MPI_REAL8,               &
                                  istat,                   &
                                  ierr)
 
          end if ! block is non-zero on input file
        end if ! block has non-zero length in general

!       keep track of offset
        ioff_luin = ioff_luin + blk_len

      end do ! loop over blocks in batch

   end subroutine read_batch_pario

!******************************************************************************
  subroutine write_batch_pario(xmat,                &
                               luout,               &
                               luoutlist,           &
                               num_blocks_in_batch, &
                               batch_info,          &
                               block_info,          &
                               ioff_luout,          &
                               luoutlist_offset,    &
                               ioff_luout_special)
!******************************************************************************
!
!    purpose:  write batch of c/sigma-vector to file luout:
!                - update the file lists
!                - active blocks on the MPI-files are flagged 
!                  by a nonzero length
!
!******************************************************************************
     real(8), intent(in)    :: xmat(*)
!    integer, intent(in)    :: luout ! stefan: open mpi 1.5 with f90-bindings complains here
     integer                :: luout
     integer, intent(in)    :: luoutlist(*)
     integer, intent(in)    :: num_blocks_in_batch
     integer, intent(in)    :: batch_info(8,*)
     integer, intent(in)    :: block_info(*)
!    offset
     integer(kind=MPI_OFFSET_KIND), intent(inout) :: ioff_luout
     integer(kind=MPI_OFFSET_KIND), intent(inout) :: ioff_luout_special
     integer, intent(in)                          :: luoutlist_offset
!------------------------------------------------------------------------------
     integer(kind=MPI_OFFSET_KIND)  :: blk_len
     integer                        :: mem_off
     integer                        :: is_blk
!******************************************************************************

      do is_blk = 1, num_blocks_in_batch

        blk_len = luoutlist(luoutlist_offset+is_blk-1)

        if(luoutlist(luoutlist_offset+is_blk-1) > 0)then

!           set memory offset
            mem_off = batch_info(6,is_blk)

!           write vector
            call mpi_file_write_at(luout,                                   &
                                   ioff_luout,                              &
                                   xmat(mem_off),                           &
                                   luoutlist(luoutlist_offset+is_blk-1),    &
                                   MPI_REAL8,                               &
                                   istat,                                   &
                                   ierr)
 
        end if ! block has non-zero length in general

!       keep track of offset
        blk_len            = block_info(luoutlist_offset+is_blk-1)
        ioff_luout         = ioff_luout         + blk_len
        ioff_luout_special = ioff_luout_special + blk_len

      end do ! loop over blocks in batch

   end subroutine write_batch_pario
!******************************************************************************

#else
module dummy_mcci_pario
#endif
end module
