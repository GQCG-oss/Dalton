!dalton_copyright_start
!
!
!dalton_copyright_end

#ifdef VAR_MPI
!**********************************************************************
!
!============================================================*
!  Common LUCITA and LUCIAREL routines for parallel purposes *
!                                                            *
!   collected by Stefan Knecht, Jan 18 - ???, 2007           *
!============================================================*
!
      subroutine copy_vector_mpi_2_sequential_io(luin,luout,vec1,
     &                                           my_ioff_luin,
     &                                           luinlist,
     &                                           par_dist_block_list,
     &                                           block_list,
     &                                           communicator_group,
     &                                           nblock,nvec,
     &                                           i_run_cplx)
C
C     Written by  S. Knecht         - June 11 2007
C
C**********************************************************************
C
C     copy vector from MPI file LUIN to 'normal file' LUOUT blockwise
C
C     NOTE: nvec = nvec
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      real(8), intent(inout) :: vec1(*) 
      integer, intent(in)    :: nblock, nvec, i_run_cplx
      integer, intent(in)    :: par_dist_block_list(nblock)
      integer, intent(in)    :: block_list(*)
      integer, intent(in)    :: communicator_group
      integer, intent(in)    :: luin, luout
      integer, intent(in)    :: luinlist(*)
      INTEGER(KIND=MPI_OFFSET_KIND), intent(in) :: MY_IOFF_LUIN
!-------------------------------------------------------------------------------
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
      integer sender_id
      integer ivec
!-------------------------------------------------------------------------------
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_INT_IN    = 0
C
C     initialize co-worker tag 
      sender_id = 0
C
#ifdef LUCI_DEBUG
      print *, 'I_RUN_CPLX',I_RUN_CPLX
      print *, 'NBLOCK',NBLOCK
      print *, 'par_dist list'
      call iwrtma(par_dist_block_list,1,NBLOCK,1,NBLOCK)
      print *, 'lbl_dist list'
      call iwrtma(block_list,1,NBLOCK,1,NBLOCK)
#endif
   
C
C     ================
C      COMPLEX VECTOR
C     ================

C     loop over vectors
      DO IVEC = 1, NVEC
C
        IOFFSET_SCRATCH = 0
        NUM_BLK         = 0
C
C       ...............
C       REAL PART FIRST 
C       ...............
C
C       set new offset
C
C       position in file is at the end of vector nvec - 1
C
C       note: real part: --> MY_VEC2_IOFF and MY_ACT_BLK2
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( ivec - 1 )   * MY_VEC2_IOFF
        IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                   ( ivec - 1 )   * MY_ACT_BLK2
C
        DO II = 1, NBLOCK ! loop over vector blocks

          IOFFSET_SCRATCH = 0
          NUM_BLK         = 0

          IF(LUCI_MYPROC .ne. LUCI_MASTER)THEN ! co-workers first ...

            IF(par_dist_block_list(II) .ne. LUCI_MYPROC) GOTO 500

          ELSE ! LUCI_MASTER

!           receive block from co-worker sender_id
            sender_id = par_dist_block_list(II)

            if(abs(sender_id) .gt. 0)then

              if(sender_id .eq. -2)then
!               block has zero length
                CALL ITODS(0,1,-1,LUOUT)
                IMZERO  = 1
                IAMPACK = 0
                CALL ZERORC(0,LUOUT,IAMPACK)
                GOTO 500
              end if

            else
 
!             LUCI_MASTER block
!             -----------------
              LBL = block_list(II)
   
              CALL ITODS(LBL,1,-1,LUOUT)
              IMZERO = 1
              IF(LUINLIST(IOFFSET_INT_IN) .gt. 0)THEN
C
C               read block in
                CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,vec1,LBL,
     &                                MPI_REAL8,my_STATUS,IERR)
#ifdef LUCI_DEBUG
                print *,'master - read-in: my vector piece'
                call wrtmatmn(vec1,1,lbl,1,lbl,luwrt)
#endif
                IMZERO = 0
              END IF
C
              IF(IMZERO.EQ.0) THEN
                 CALL TODSC_LUCI(vec1,LBL,-1,LUOUT)
              ELSE
                 CALL ZERORC(LBL,LUOUT,0)
              END IF
C
              IOFFSET_SCRATCH = LBL
              NUM_BLK         = 1
C
              GOTO 500
            END IF ! abs(sender_id) > 0?

          END IF ! co-worker / master switch
 
!         co-worker block
!         ---------------
          LBL = block_list(II)
 
          if(luci_myproc .ne. luci_master)then
C
            ISENDTOM = 0
            ISENDTOM = LUINLIST(IOFFSET_INT_IN)
C
            CALL MPI_SEND(ISENDTOM,1,my_MPI_INTEGER,LUCI_MASTER,90,
     &                    communicator_group,IERR)
C
            IF(ISENDTOM .gt. 0)THEN
C
              CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,vec1,LBL,
     &                              MPI_REAL8,my_STATUS,IERR)
              CALL MPI_SEND(vec1,LBL,MPI_REAL8,LUCI_MASTER,91,
     &                      communicator_group,IERR)
#ifdef LUCI_DEBUG
                 print *,'co-worker - read-in: my vector piece'
                 call wrtmatmn(vec1,1,lbl,1,lbl,luwrt)
#endif
C
            END IF
            IOFFSET_SCRATCH = LBL
            NUM_BLK         = 1

          else ! LUCI_MASTER
 
            CALL ITODS(LBL,1,-1,LUOUT)

            ISENDTOM = 0
            CALL MPI_RECV(ISENDTOM,1,my_MPI_INTEGER,sender_id,90,
     &                    communicator_group,my_STATUS,IERR)
C
            IF(ISENDTOM .gt. 0)THEN
C
              CALL MPI_RECV(vec1,LBL,MPI_REAL8,sender_id,91,
     &                       communicator_group,my_STATUS,IERR)
#ifdef LUCI_DEBUG
              print *,'master- recieve: my vector piece'
              call wrtmatmn(vec1,1,lbl,1,lbl,luwrt)
#endif
              CALL TODSC_LUCI(vec1,LBL,-1,LUOUT)
            ELSE
              CALL ZERORC(LBL,LUOUT,0)
            END IF
            IOFFSET_SCRATCH = 0
            NUM_BLK         = 0

          end if ! co-worker / master switch

 500      CONTINUE
        
!         keep track of correct offset
          IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
          IOFFSET_INT_IN  = IOFFSET_INT_IN  + NUM_BLK

        end do ! loop over vector blocks
C
        IF(LUCI_MYPROC. eq. LUCI_MASTER)THEN
C         write an EOF-mark
          CALL ITODS(-1,1,-1,LUOUT)
        ENDIF
C
C       ..............
C          IMAG PART 
C       ..............
C
        IOFFSET_SCRATCH = 0
        NUM_BLK         = 0
C
        IF(I_RUN_CPLX .eq. 2)THEN
C
C         set new offset
C
C         position in file is at the end of vector ivec - 1
C
C         note: complex part: --> MY_VEC1_IOFF and MY_ACT_BLK1
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                     ( ivec - 1 )   * MY_VEC2_IOFF    +
     &                                      MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                     ( ivec - 1 )   * MY_ACT_BLK2 + MY_ACT_BLK1
C
C LOOP OVER BLOCKS
C
         DO 601 II = 1, NBLOCK 
C
C           slaves first ...
C
           IF( LUCI_MYPROC .ne. LUCI_MASTER )THEN
             IF( par_dist_block_list(II) .ne. LUCI_MYPROC )THEN 
C
               IOFFSET_SCRATCH = 0
               NUM_BLK         = 0
C
               GOTO 600
             ELSE IF( par_dist_block_list(II) .eq. LUCI_MYPROC )THEN
               GOTO 102
             ENDIF
C
C          LUCI_MASTER follows ...
C
           ELSE IF( LUCI_MYPROC .eq. LUCI_MASTER )THEN
             IF( par_dist_block_list(II) .ne. LUCI_MYPROC )THEN
               IF(par_dist_block_list(II) .eq. -2 )THEN
C
C                block has zero length
C
                 LBL = 0
                 CALL ITODS(LBL,1,-1,LUOUT)
                 IMZERO = 1
                 IAMPACK = 0
                 CALL ZERORC(LBL,LUOUT,IAMPACK)
C
                 IOFFSET_SCRATCH = 0
                 NUM_BLK = 0
C
                 GOTO 600
               ELSE IF( par_dist_block_list(II) .gt. 0 )THEN
C
C                recieve block from node sender_id
C
                 sender_id = par_dist_block_list(II) 
                 GOTO 102 
               END IF
             ELSE IF( par_dist_block_list(II) .eq. LUCI_MYPROC )THEN
C
C              LUCI_MASTER block
C
               LBL = block_list(II)
C  
               CALL ITODS(LBL,1,-1,LUOUT)
C
               IF(LUINLIST(IOFFSET_INT_IN) .gt. 0)THEN
C
C                read block in
                 CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,vec1,LBL,
     &                                 MPI_REAL8,my_STATUS,IERR)
                 IMZERO = 0
               ELSE
                 IMZERO = 1
               END IF
C
               IF(IMZERO .eq. 0)THEN
                 call todsc_luci(vec1,LBL,-1,LUOUT)
               ELSE
                 CALL ZERORC(LBL,LUOUT,0)
               END IF
C
               NUM_BLK = 1
               IOFFSET_SCRATCH = LBL
C
               GOTO 600
             END IF ! par_dist_block_list eq 0 or ne 0
           END IF
C
C       ... continue
C
 102        CONTINUE
C
C           slave block
C
            LBL = block_list(II)
C
            IF( LUCI_MYPROC .eq. LUCI_MASTER ) THEN
C
               CALL ITODS(LBL,1,-1,LUOUT)
C  
            END IF
C
            IF( LUCI_MYPROC .ne. LUCI_MASTER ) THEN
C
              ISENDTOM = 0
              ISENDTOM = LUINLIST( IOFFSET_INT_IN )
C
              CALL MPI_SEND(ISENDTOM,1,my_MPI_INTEGER,LUCI_MASTER,93,
     &                      communicator_group,IERR)
C
              IF( ISENDTOM .gt. 0 ) THEN
C
                CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN_LUIN,vec1,LBL,
     &                                MPI_REAL8,my_STATUS,IERR)
                CALL MPI_SEND(vec1,LBL,MPI_REAL8,LUCI_MASTER,94,
     &                        communicator_group,IERR)
C
              END IF
              IOFFSET_SCRATCH = LBL
              NUM_BLK = 1
            END IF
C
C
            IF( LUCI_MYPROC .eq. LUCI_MASTER ) THEN
C
              ISENDTOM = 0
C
              CALL MPI_RECV(ISENDTOM,1,my_MPI_INTEGER,sender_id,93,
     &                      communicator_group,my_STATUS,IERR)
C
              IF(ISENDTOM .gt. 0)THEN
C
                CALL MPI_RECV(vec1,LBL,MPI_REAL8,sender_id,94,
     &                        communicator_group,my_STATUS,IERR)
                CALL TODSC_LUCI(vec1,LBL,-1,LUOUT)
              ELSE
                CALL ZERORC(LBL,LUOUT,0)
              END IF
              IOFFSET_SCRATCH = 0
              NUM_BLK         = 0
            END IF
C
 600        CONTINUE
C       
C           keep track of correct offset
            IOFFSET_IN_LUIN = IOFFSET_IN_LUIN + IOFFSET_SCRATCH
            IOFFSET_INT_IN  = IOFFSET_INT_IN  + NUM_BLK
 601      CONTINUE
C
          IF( LUCI_MYPROC. eq. LUCI_MASTER )THEN
C           write an EOF-mark
            CALL ITODS(-1,1,-1,LUOUT)
          ENDIF
C
        END IF ! COMPLEX PART

      END DO ! loop over vectors
C
      END
!**********************************************************************

      subroutine copy_vector_sequential_2_mpi_io(LUIN,LUOUT,VEC1,
     &                                           MY_IOFF_LUOUT,
     &                                           LUOUTLIST,
     &                                           par_dist_block_list,
     &                                           block_list,
     &                                           communicator_group,
     &                                           nblocks,NVEC,
     &                                           I_RUN_CPLX)
C
C     Written by  S. Knecht         - June 19 2007
C
C**********************************************************************
C
C     copy vector from sequential file LUIN to MPI-file LUOUT.
C     master sends a block to a co-worker according to par_dist_block_list.
C
C     NOTE: NVEC = NVEC
C
C     active blocks on the MPI-file are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
!------------------------------------------------------------------------------
      real(8), intent(inout) :: vec1(*)
      integer, intent(in)    :: nblocks, nvec, i_run_cplx, luin, luout
      integer, intent(out)   :: luoutlist(*)
      integer, intent(in)    :: par_dist_block_list(nblocks) 
      integer, intent(in)    :: block_list(*)
      integer, intent(in)    :: communicator_group
      integer(kind=MPI_OFFSET_KIND), intent(in) :: MY_IOFF_LUOUT
!------------------------------------------------------------------------------
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER                :: NUM_BLK     = 0 
      integer                :: receiver_id = 0
      integer                :: ilen        = 0
      integer, parameter     :: no_zeroing  = 1
!------------------------------------------------------------------------------

C
C     initialize scratch offsets
      IOFFSET_SCRATCH   = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_LUOUT = 0
#ifdef LUCI_DEBUG
      if(luci_myproc .eq. luci_master)then
        call rewino(LUIN)
        write(6,*) 'vectors to copy'
        do i = 1, NVEC
          call wrtvcd(vec1,luin,0,-1)
        end do
        call rewino(LUIN)
      end if
#endif
 
C
C     loop over vectors
      DO IVEC = 1, NVEC
C
        IOFFSET_SCRATCH = 0
        NUM_BLK         = 0
C
C       ================
C        COMPLEX VECTOR
C       ================
C
C       .................
C        REAL PART FIRST
C       .................
C
C       set new offset
C
C       position in file is at the end of vector IVEC - 1
C
C       note: real part: --> MY_VEC2_IOFF and MY_ACT_BLK2
C
        IOFFSET_LUOUT      = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                       ( IVEC - 1 )  * MY_VEC2_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK  +
     &                       ( IVEC - 1 ) * MY_ACT_BLK2
C
C       loop over blocks
C
        DO IBLK = 1, NBLOCKS
C
          ILEN = block_list(IBLK)

          IF(ILEN .GT. 0)THEN
C
C           get CPU id that will receive the block
            receiver_id = par_dist_block_list(IBLK)
C
            IF(LUCI_MYPROC .ne. LUCI_MASTER)THEN ! co-worker part...
C
              if(LUCI_MYPROC .ne. receiver_id) GOTO 500

              NUM_BLK                                    = NUM_BLK + 1
              LUOUTLIST(IOFFSET_INT_LUOUT + NUM_BLK - 1) = 0
              IMZERO                                     = 0

              CALL MPI_RECV(IMZERO,1,my_MPI_INTEGER,LUCI_MASTER,90,
     &                      communicator_group,my_STATUS,IERR)

              IF(IMZERO .eq. 0)THEN
C               recieve block and transfer to file
                CALL MPI_RECV(VEC1,ILEN,MPI_REAL8,LUCI_MASTER,91,
     &                        communicator_group,my_STATUS,IERR)
csk             WRITE(6,*) ' S will write this to disk for block',
csk     &       IOFFSET_INT_LUOUT + NUM_BLK - 1, IBLK
csk             CALL WRTMATMN(VEC1,1,ILEN,1,ILEN,6)
                CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_LUOUT,VEC1,ILEN,
     &                                 MPI_REAL8,my_STATUS,IERR)

                LUOUTLIST(IOFFSET_INT_LUOUT + NUM_BLK - 1) = 1
              END IF
C
C             keep track of correct offset
              IOFFSET_SCRATCH = ILEN
              IOFFSET_LUOUT   = IOFFSET_LUOUT + IOFFSET_SCRATCH

            ELSE ! LUCI_MASTER part...
C
C             read block in core and send to CPU receiver_id if necessary
              IMZERO = 0
              CALL IFRMDS(LBL,1,-1,LUIN)
              IF(LBL .ne. ILEN)THEN
                WRITE(LUWRT,*) ' Error in COP_REST_VEC_REL:'//
     &                         ' block sizes do not match. '
                call quit('Error in COP_REST_VEC_REL detected!')
              END IF
              CALL FRMDSC2(VEC1,ILEN,-1,LUIN,IMZERO,IAMPACK,NO_ZEROING)

              IF(LUCI_MYPROC .ne. receiver_id)THEN ! co-workers block
                CALL MPI_SEND(IMZERO,1,my_MPI_INTEGER,receiver_id,90,
     &                        communicator_group,IERR)
                IF(IMZERO .eq. 0)THEN
                  CALL MPI_SEND(VEC1,ILEN,MPI_REAL8,receiver_id,91,
     &                          communicator_group,IERR)
                END IF
              ELSE ! master block

C               count active blocks
                NUM_BLK                                    = NUM_BLK + 1
                LUOUTLIST(IOFFSET_INT_LUOUT + NUM_BLK - 1) = 0
C
                IF(IMZERO .eq. 0)THEN
csk               WRITE(LUWRT,*) ' M will write this to disk for block',
csk  &            IOFFSET_INT_LUOUT + NUM_BLK - 1, IBLK
csk               CALL WRTMATMN(VEC1,1,ILEN,1,ILEN,LUWRT)
                  CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_LUOUT,VEC1,ILEN,
     &                                 MPI_REAL8,my_STATUS,IERR)
                  LUOUTLIST( IOFFSET_INT_LUOUT + NUM_BLK - 1 ) = 1
                END IF 

!               keep track of correct offset
                IOFFSET_SCRATCH = ILEN
                IOFFSET_LUOUT   = IOFFSET_LUOUT + IOFFSET_SCRATCH
              END IF ! send id /= master?
            END IF ! co-workers / master switch
          ELSE
C           skip zero block on disk
            IF(LUCI_MYPROC .EQ. LUCI_MASTER)THEN
              CALL IFRMDS(LBL,1,-1,LUIN)
              CALL FRMDSC2(VEC1,LBL,-1,LUIN,IMZERO,IAMPACK,NO_ZEROING)
            END IF
          END IF ! length of block > 0

 500     CONTINUE
C
        END DO ! loop over vector blocks
C
C       ..............
C        COMPLEX PART
C       ..............
C
        IOFFSET_SCRATCH = 0
        NUM_BLK = 0
C
        IF( I_RUN_CPLX .eq. 2 ) THEN
C
          IF( LUCI_MYPROC .eq. LUCI_MASTER )THEN
C           skip marker on file 
            CALL IFRMDS(LBL,1,-1,LUIN)
C
          END IF
C
          IOFFSET_SCRATCH = 0
          NUM_BLK = 0
C
C         new offset for writing on LUOUT - complex part
C
          IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                          ( IVEC - 1 )  * MY_VEC2_IOFF +
     &                                          MY_VEC1_IOFF
C
          IOFFSET_INT_LUOUT  = 1 + NUM_BLK  +
     &                         ( IVEC - 1 ) * MY_ACT_BLK2 + MY_ACT_BLK1
C
C         loop over blocks
C
          DO IBLK = 1, NBLOCKS
C
            ILEN = block_list( IBLK )
CSK            WRITE(LUWRT,*) ' ILEN of IBLK',ILEN,IBLK
            IF( ILEN .GT. 0 ) THEN
C
C             get number of CPU that will recieve the block
              receiver_id = par_dist_block_list( IBLK )
CSK              WRITE(LUWRT,*) ' receiver_id of IBLK',receiver_id,IBLK
C
              IF( LUCI_MYPROC .ne. LUCI_MASTER .AND.
     &            LUCI_MYPROC .NE. receiver_id ) 
     &           GOTO 600
!
              IF(LUCI_MYPROC .ne. LUCI_MASTER)THEN ! co-workers
!
                NUM_BLK = NUM_BLK + 1
                IMZERO = 0
!
                CALL MPI_RECV(IMZERO,1,my_MPI_INTEGER,LUCI_MASTER,93,
     &                        communicator_group,my_STATUS,IERR)
!
                IF( IMZERO .eq. 0 ) THEN
!
!                 recieve block and transfer to file
                  CALL MPI_RECV(VEC1,ILEN,MPI_REAL8,LUCI_MASTER,94,
     &                          communicator_group,my_STATUS,IERR)
                  LUOUTLIST( IOFFSET_INT_LUOUT + NUM_BLK - 1 ) = 1
                  CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_LUOUT,VEC1,ILEN,
     &                                   MPI_REAL8,my_STATUS,IERR)
!
                ELSE
!
!                 mark block as zero block
                  LUOUTLIST( IOFFSET_INT_LUOUT + NUM_BLK - 1 ) = 0
!
                END IF
!
!               keep track of correct offset
                IOFFSET_SCRATCH = ILEN
                IOFFSET_LUOUT = IOFFSET_LUOUT + IOFFSET_SCRATCH
!
              ELSE ! LUCI_MASTER
!
!               read block in core and send to CPU receiver_id if necessary
                IMZERO = 0
                CALL IFRMDS(LBL,1,-1,LUIN)
                IF( LBL .ne. ILEN ) THEN
                  WRITE(LUWRT,*) '  Error in COP_REST_VEC_REL:'//
     &                           ' block sizes do not match! '
                  call quit('Error in COP_REST_VEC_REL detected!')
                END IF
                CALL FRMDSC2(VEC1,ILEN,-1,LUIN,IMZERO,IAMPACK,
     &                       NO_ZEROING)
                IF( LUCI_MYPROC .ne. receiver_id )THEN
                  CALL MPI_SEND(IMZERO,1,my_MPI_INTEGER,receiver_id,93,
     &                          communicator_group,IERR)
                  IF(IMZERO .eq. 0)THEN
                    CALL MPI_SEND(VEC1,ILEN,MPI_REAL8,receiver_id,94,
     &                            communicator_group,IERR)
                  END IF
                ELSE
!                 count active blocks
                  NUM_BLK = NUM_BLK + 1
!
                  IF( IMZERO .eq. 0 ) THEN
!
                    CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_LUOUT,VEC1,
     &                                  ILEN,MPI_REAL8,my_STATUS,IERR)
                    LUOUTLIST( IOFFSET_INT_LUOUT + NUM_BLK - 1 ) = 1
                  ELSE
                    LUOUTLIST( IOFFSET_INT_LUOUT + NUM_BLK - 1 ) = 0
                  END IF 
!                 keep track of correct offset
                  IOFFSET_SCRATCH = ILEN
                  IOFFSET_LUOUT = IOFFSET_LUOUT + IOFFSET_SCRATCH
!
                END IF ! needs send?
!
              END IF ! co-worker or master
           ELSE
!            skip zero block on disk
             IF( LUCI_MYPROC .EQ. LUCI_MASTER ) THEN
               CALL IFRMDS(LBL,1,-1,LUIN)
               CALL FRMDSC2(VEC1,LBL,-1,LUIN,IMZERO,IAMPACK,
     &                      NO_ZEROING)
             END IF
           END IF ! ILEN > 0?
 600       CONTINUE
!
          END DO ! loop over blocks
!
        END IF ! COMPLEX PART NEEDED? 
!
        IF( LUCI_MYPROC .eq. LUCI_MASTER ) THEN
!         skip end-of-vector marker on file LUIN
          CALL IFRMDS(LBL,1,-1,LUIN)
        END IF
!
      END DO ! loop over vectors
!
      END
***********************************************************************

      SUBROUTINE COPVCD_PAR_BDRIV5_REL(LUIN,LUOUT,SEGMNT,IBLOCKD,
     &                                 ISCALFAC,ISCALFAC_GROUP,
     &                                 IBLOCKL,NBLOCK,JCOMM,
     &                                 IGROUPLIST,IPROCLIST,IRILP)
C
C     Read blocks of the vector from disc file LUIN to array segment
C     and broadcast to nodes that need the block, save on diskfile luout 
C     if necessary. 
C     communication via communicator jcomm.
C
C     OUTPUT: update of LUOUT with new TTSS-c-blocks, 
C             modified ISCALFAC_GROUP.
C
C     based on the corresponding spinfree-routine
C
C     extended for an imaginary part of the c-vector
C     if NZ == IRILP == 2           S. Knecht - June 27 2007
C
C     Written by  S. Knecht         - May 23 2007
C
C**********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*), IBLOCKD(NBLOCK), IBLOCKL(NBLOCK)
      DIMENSION ISCALFAC(*), ISCALFAC_GROUP(*)
      DIMENSION IGROUPLIST(LUCI_NMPROC),IPROCLIST(LUCI_NMPROC)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCRATCH_P
      INTEGER IGROUPBLK, ROOTPROC, IOFF_BLOCK
      integer is_new_rootproc_tmp, is_new_rootproc
      CHARACTER*12 WALLTID, SECTID
C
      bcast_time = 0.0D0
      read_write = 0.0D0
C
      IBI_DIST_NBLOCKS = 0
      IBI_RECV_NBLOCKS = 0
C
      IOFF_BLOCK     = 0
      IOFF_SCRATCH_P = 0
csk   WRITE(LUWRT,*) ' MY GROUPLIST resp PROCLIST'
csk   CALL IWRTMAMN(IGROUPLIST,1,LUCI_NMPROC,1,LUCI_NMPROC,LUWRT)
csk   CALL IWRTMAMN(IPROCLIST,1,LUCI_NMPROC,1,LUCI_NMPROC,LUWRT)
C   
C     =========
C     REAL PART
C     =========
C
      ROOTPROC   = 0
C
C     file offset for each cpu is 0 - real part
C
      IOFFSET_IN  = 0
      IOFFSET_OUT = 0
C
C     LOOP OVER BLOCKS
C
      DO IBLK = 1, NBLOCK 
C
C       get length of block
C
        LBL = IBLOCKL(IBLK)
C
C       check if block is active
C
        IF( ISCALFAC(IBLK) .ne. 0 ) THEN
C
C         check if block belongs to own group
C
          IGROUPBLK = 0
          ISCRNODE = IBLOCKD(IBLK)
C
          DO IPROC = 1, NEWCOMM_PROC
            IF( ISCRNODE .eq. IGROUPLIST(IPROC) ) IGROUPBLK = 1
          END DO
          IF( IGROUPBLK .eq. 1 ) THEN
C
csk         WRITE(6,*) 'ACTIVE BLOCK ',IBLK
csk         WRITE(6,*) 'START READING AT IOFFSET_IN = ',IOFFSET_IN
            xrw_time = MPI_WTIME()
            CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN,SEGMNT,LBL,
     &                            MPI_REAL8,my_STATUS,IERR)
csk         WRITE(6,*) 'MY read-in result '
csk         CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,6)
            read_write = read_write + MPI_WTIME() - xrw_time
            IBI_DIST_NBLOCKS = IBI_DIST_NBLOCKS + 1
C
          END IF
C
C         find the sending root processor
          ROOTPROC = IPROCLIST(ISCRNODE+1) - 1
C
C         decide whether it is necessary to take part in communication
!
!         stefan: i commented the split/free code out because of
!         deadlock situations occuring with buggy mpi implementations on horseshoe
!         at sdu in odense.
!
          NNKEY   = LUCI_MYPROC + 1
          IF(LUCI_MYPROC .eq. ROOTPROC) NNKEY = 0
          NNCOLOR = 8
          IF(ISCALFAC_GROUP(IBLK) .ne. 0) NNCOLOR = 7

!         CALL MPI_COMM_SPLIT(JCOMM,NNCOLOR,NNKEY,IBLOCKCOMM,IERR)
 
!         IF( NNCOLOR .eq. 8 ) GOTO 200

C         broadcast the nonvanishing block
          xcast_time = MPI_WTIME()
!         CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,0,IBLOCKCOMM,IERR)

!         new code without mpi_comm_split
          is_new_rootproc     = 0
          is_new_rootproc_tmp = 0
          IF(LUCI_MYPROC .eq. ROOTPROC) is_new_rootproc_tmp = icomm_id
          CALL MPI_ALLREDUCE(is_new_rootproc_tmp,is_new_rootproc,1,
     &                       my_MPI_INTEGER,mpi_max,JCOMM,IERR)
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,is_new_rootproc,
     &                   JCOMM,IERR)
!         end of new code without mpi_comm_split

          bcast_time = bcast_time + MPI_WTIME() - xcast_time

C         transfer the recieved block to disk
#ifdef LUCI_DEBUG
          WRITE(LUWRT,*)'THIS IS WHAT I GOT',ICOMM_ID
          CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,LUWRT)
          WRITE(LUWRT,*) 'IOFFSET_OUT = ',IOFFSET_OUT
#endif
          if(NNCOLOR .eq. 7)then
            xrw_time = MPI_WTIME()
            CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_OUT,SEGMNT,LBL,
     &                             MPI_REAL8,my_STATUS,IERR)
            read_write = read_write + MPI_WTIME() - xrw_time
            IBI_RECV_NBLOCKS = IBI_RECV_NBLOCKS + 1
          end if
C
!200      CONTINUE
!         CALL MPI_COMM_FREE(IBLOCKCOMM,IERR)
C
        ELSE 
C         set local iscalfac to 0
          ISCALFAC_GROUP(IBLK) = 0
        END IF
C
C       keep track of correct offset
        IOFFSET_OUT = IOFFSET_OUT + LBL
        IOFFSET_IN  = IOFFSET_IN  + LBL
C
      END DO
C     ^ loop over blocks - real part
C
C
      IF( IRILP .eq. 1 ) GOTO 1001
C
C     =========
C     IMAG PART
C     =========
C
      ROOTPROC = 0
C
C     pointer for file offset: skip real part of the c-vector
C
      IOFF_SCRATCH_P = L_COMBI
C
      IOFFSET_IN  = IOFF_SCRATCH_P
      IOFFSET_OUT = IOFF_SCRATCH_P
C     set correct block counter for imag part
      IOFF_BLOCK = NBLOCK
CSK      WRITE(LUWRT,*) 'IOFF_BLOCK is',IOFF_BLOCK
C
C     LOOP OVER BLOCKS
C
      DO IBLK = 1, NBLOCK 
C
C       get length of block
C
        LBL = IBLOCKL( IBLK )
C
C       check if block is active
C
        IF( ISCALFAC( IBLK + IOFF_BLOCK ) .ne. 0 ) THEN
C
CSK          WRITE(LUWRT,*) 'ACTIVE BLOCK imag part',IBLK, 
CSK     &                    IBLK + IOFF_BLOCK
C         check if block belongs to own group
C
          IGROUPBLK = 0
          ISCRNODE = IBLOCKD( IBLK )
C
          DO IPROC = 1, NEWCOMM_PROC
            IF( ISCRNODE .eq. IGROUPLIST( IPROC ) ) IGROUPBLK = 1
          END DO
          IF( IGROUPBLK .eq. 1 ) THEN
C
CSK            WRITE(LUWRT,*) 'START READING AT IOFFSET_IN imag part = ',
CSK     &                      IOFFSET_IN
            xrw_time = MPI_WTIME()
            CALL MPI_FILE_READ_AT(LUIN,IOFFSET_IN,SEGMNT,LBL,
     &                            MPI_REAL8,my_STATUS,IERR)
            read_write = read_write + MPI_WTIME() - xrw_time
            IBI_DIST_NBLOCKS = IBI_DIST_NBLOCKS + 1
C
          END IF
C
C         find the sending root processor
          ROOTPROC = IPROCLIST(ISCRNODE+1) - 1
C
C         decide whether it is necessary to take part in communication
!
!         stefan: i commented the split/free code out because of
!         deadlock situations occuring with buggy mpi implementations on horseshoe
!         at sdu in odense.
!
          NNKEY   = LUCI_MYPROC + 1
          IF(LUCI_MYPROC .eq. ROOTPROC) NNKEY = 0
          NNCOLOR = 8
          IF(ISCALFAC_GROUP(IOFF_BLOCK+IBLK) .ne. 0) NNCOLOR = 7

!         CALL MPI_COMM_SPLIT(JCOMM,NNCOLOR,NNKEY,IBLOCKCOMM,IERR)

!         IF( NNCOLOR .eq. 8 ) GOTO 300

C         broadcast the nonvanishing block
          xcast_time = MPI_WTIME()
!         CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,0,IBLOCKCOMM,IERR)

!         new code without mpi_comm_split
          is_new_rootproc     = 0
          is_new_rootproc_tmp = 0
          IF(LUCI_MYPROC .eq. ROOTPROC) is_new_rootproc_tmp = icomm_id
          CALL MPI_ALLREDUCE(is_new_rootproc_tmp,is_new_rootproc,1,
     &                       my_MPI_INTEGER,mpi_max,JCOMM,IERR)
          CALL MPI_BCAST(SEGMNT,LBL,MPI_REAL8,is_new_rootproc,
     &                   JCOMM,IERR)
!         end of new code without mpi_comm_split
          bcast_time = bcast_time + MPI_WTIME() - xcast_time

C         transfer the recieved block to disk if necessary
          if(nncolor .eq. 7)then
#ifdef LUCI_DEBUG
            WRITE(LUWRT,*)'THIS IS WHAT I GOT from ROOTPROC',
     &                     LUCI_MYPROC,ROOTPROC
            CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,LUWRT)
            WRITE(LUWRT,*) 'IOFFSET_OUT imag part = ',IOFFSET_OUT
            WRITE(LUWRT,*) 'block to write imag part',IBLK,
     &                      IBLK + IOFF_BLOCK
#endif
            xrw_time = MPI_WTIME()
            CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET_OUT,SEGMNT,LBL,
     &                             MPI_REAL8,my_STATUS,IERR)
            read_write = read_write + MPI_WTIME() - xrw_time
            IBI_RECV_NBLOCKS = IBI_RECV_NBLOCKS + 1
          end if
C
!300      CONTINUE
!         CALL MPI_COMM_FREE(IBLOCKCOMM,IERR)
C
        ELSE 
C         set local iscalfac to 0
          ISCALFAC_GROUP( IBLK + IOFF_BLOCK ) = 0
        END IF
C
C       keep track of correct offset
        IOFFSET_OUT = IOFFSET_OUT + LBL
        IOFFSET_IN  = IOFFSET_IN  + LBL
C
      END DO
C     ^ loop over blocks - imag part
C
 1001 CONTINUE
C
      IF( TIMING_par ) THEN
C
C       print statistics
C
        WRITE(LUWRT,'(/A)') 
     &  '               coefficients exchange statistics  '
        WRITE(LUWRT,'(A/)')   
     &  '              __________________________________ '
        WRITE(LUWRT,'(2X,A,1X,I4)')   
     &  ' number of b_i blocks distributed         : ',IBI_DIST_NBLOCKS
        WRITE(LUWRT,'(2X,A,1X,I4)')   
     &  ' number of b_i blocks recieved            : ',IBI_RECV_NBLOCKS
        WALLTID = SECTID(bcast_time)
        WRITE(LUWRT,'(2X,A,1X,A)')   
     &  ' time spent in b_i communication (bcast)  : ', WALLTID
        WALLTID = SECTID(read_write)
        WRITE(LUWRT,'(2X,A,1X,A/)')   
     &  ' time spent in saving coefficients on disk: ', WALLTID
C
      END IF
C
      END
!**********************************************************************

      SUBROUTINE H0M1TD_REL_PAR(LUOUT,LUDIA,LUIN,SHIFT,VEC1,VEC2,
     &                          LUOUTLIST,LUINLIST,
     &                          NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                          IBATCH,MY_IOFF_LUOUT,MY_IOFF_LUDIA,
     &                          MY_IOFF_LUIN,INV)
C
C     calculate inverted general preconditioner matrix times vector
C
C     original written by Jeppe Olsen - September 1993
C
C     adaption of sequential routine for parallel purposes  
C                      by S. Knecht   -  March 13 2008
C
C     MPI file I/O version
C
C     vecout=  (H0 + shift )-1 * vecin
C
C      LUOUT       LUDIA        LUIN
C
C**********************************************************************
#include "implicit.h"
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
*
* =====
* Input
* =====
*
* LUOUT : File to contain output vector
* LUDIA : File containing diagonal of H0
* LUIN  : File containing input vector
* SHIFT : constant ADDED to diagonal
*
* ======
* Output
* ======
*
* LUOUT : contains output vector, not rewinded
*
* =======
* Scratch
* =======
*
* VEC1,VEC2 : Must each be able to hold largest segment of vector
C
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUDIA
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_IN_LUDIA
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_IN_LUDIA  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_LUOUT = 0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
        IOFFSET_IN_LUIN = MY_IOFF_LUIN + IOFFSET_SCRATCH
        IOFFSET_INT_IN  = 1 + NUM_BLK 
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK  &                  IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
        IOFFSET_IN_LUDIA = MY_IOFF_LUDIA + IOFFSET_SCRATCH
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUDIA',
CSK  &                  IOFFSET_IN_LUDIA
C
C       read in batch ISBATCH from LUDIA to VEC1
C
        CALL RDVEC_BATCH_DRV5(LUDIA,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUDIA)
C
CSK     WRITE(LUWRT,*) 'initial VEC1 on LUDIA'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       calculate inverse diagonal on VEC1
C
        ILEN_BATCH =  0
        ILEN_BATCH =  LEBATCH(ISBATCH)
C
        IF( ILEN_BATCH .gt. 0 )THEN
           IF( CSCREEN) THEN
C             set proper truncation factor
              THR_TRUNC  = TRUNC_FAC * RNORM_FAC
!             THR_ETRUNC = 1.0D-7 * THRES_E
! THRES_E is not defined
Chj           14-jun-07:   disable THR_ETRUNC
              THR_ETRUNC = -1.0D0
!testprint    WRITE(LUWRT,*) 'TRUNCATION FACTORS:',THR_TRUNC,THR_ETRUNC
              CALL DIAVC2_TRUNC(VEC1,VEC2,VEC1,SHIFT,LEBATCH(ISBATCH),
     &                          THR_TRUNC,THR_ETRUNC)
           ELSE
             CALL DIAVC2(VEC1,VEC2,VEC1,SHIFT,LEBATCH(ISBATCH))
           END IF
        END IF
CSK          WRITE(LUWRT,*) 'final VEC1'
CSK          WRITE(LUWRT,*) '  (D-E)-1 *( HX - EX ) '
CSK          WRITE(LUWRT,*) 'final VEC2'
CSK          CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                  LUWRT)
C
C       write VEC1 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = 1 + NUM_BLK
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE COPVCD_PP_CC_B(LUIN,LUOUT,SCR,NBATCH,LBATCH,LEBATCH,
     &                         I1BATCH,IBATCH,MY_IOFF_LUIN,
     &                         MY_IOFF_LUOUT,
     &                         LUINLIST,LUOUTLIST,IBLOCKL,JOFF)
C
C     Written by  S. Knecht         - May 18 2007
C
C**********************************************************************
C
C     copy c-vector from file LUIN to LUOUT batchwise
C     update the file lists.
C
C     NOTE: JOFF = (IVEC resp. IROOT) - 1
C
C
C     copy vector from file LUIN to file LUOUT --> (hint LUIN, LUOUT)
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*),LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*), IBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) ISCRATCH_SP
      INTEGER(KIND=MPI_OFFSET_KIND) ILEN_COMB
      INTEGER NUM_BLK
C
      I_RUN_CPLX = 1
C     some constants
      IONE = 1
      ITWO = 2
C     initialize scratch offsets
      NUM_BLK = 0
      MY_NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN  = 0
      IOFFSET_OUT  = 0
      IOFFSET_INT_IN  = 0
      IOFFSET_INT_OUT  = 0
      ILEN_COMB = 0
      ILEN_COMB = L_COMBI
      MY_IOFFSET_SCRATCH = 0
      ISCRATCH_SP = 0
      NUM_BLK_SP = 0
C
C
C     ================
C       REAL VECTOR
C     ================
C
C     position in file is at the end of vector JOFF
C
C     note: real part: --> MY_VEC2_IOFF, MY_ACT_BLK2
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(SCR,LEBATCH(ISBATCH))
C
        IOFFSET_IN = MY_IOFF_LUIN + ( MY_VEC2_IOFF * JOFF ) +
     &               MY_IOFFSET_SCRATCH
        IOFFSET_INT_IN = ( MY_ACT_BLK2 * JOFF ) + MY_NUM_BLK + 1
C
C       write active block array for LUOUT 
C
        IOFFSET_INT_OUT = NUM_BLK + 1
C
        NUM_BLK_SP  = 0
        ISCRATCH_SP = 0
C
csk     WRITE(LUWRT,*) ' THIS IS MY IOFFSET_IN',IOFFSET_IN
csk     WRITE(LUWRT,*) ' THIS IS MY IOFFSET_INT_IN',IOFFSET_INT_IN
csk     WRITE(LUWRT,*) ' THIS IS MY IOFFSET_INT_OUT',IOFFSET_INT_OUT
C
        CALL RDVEC_BATCH_DRV3(LUIN,SCR,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN,IOFFSET_INT_IN,
     &                       IOFFSET_INT_OUT,LUINLIST,
     &                       LUOUTLIST,NUM_BLK_SP)
C
        IOFFSET_OUT = MY_IOFF_LUOUT + IOFFSET_SCRATCH
     
csk     WRITE(LUWRT,*) ' SCR content is'
csk     CALL WRTMATMN(SCR,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
csk     WRITE(LUWRT,*) ' THIS IS MY IOFFSET_OUT',IOFFSET_OUT
        
        IOFFSET_INT_OUT = NUM_BLK + 1
C
        CALL WTVEC_BATCH_DRV3(SCR,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        LUOUT,IOFFSET_OUT,IOFFSET_INT_OUT,
     &                        LUOUTLIST,1,ISCRATCH_SP,
     &                        IBLOCKL)
C
C
C       count the length of the last copy
C
C       LUOUT
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + ISCRATCH_SP
        NUM_BLK         = NUM_BLK + LBATCH(ISBATCH)
C       LUIN
        MY_NUM_BLK         = MY_NUM_BLK + NUM_BLK_SP
        MY_IOFFSET_SCRATCH = MY_IOFFSET_SCRATCH + LEBATCH(ISBATCH)
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE RDVEC_BATCH_DRV3(LUIN,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET,
     &                            IOFFSET_INT_IN,IOFFSET_INT_OUT,
     &                            IVCOFF_IN,IVCOFF_OUT,
     &                            ICOUNT_ACT)
C
C     Written by  S. Knecht         -      May 21 2007
C
C**********************************************************************
C
C     read in a batch from a MPI-file LUIN to SCR - C-vector routine
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*) 
      DIMENSION IVCOFF_IN(*), IVCOFF_OUT(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
        DO IBLK = 1, NBATCH_BLK
C
C         check for need
C
          ILEN2 = NBATCH_INF(8,IBLK) 
C
          IF( ILEN2 .gt. 0 ) THEN
C
C           count the active blocks for this batch
C
            ICOUNT_ACT = ICOUNT_ACT + 1
C
            IF( IVCOFF_IN( IOFFSET_INT_IN + ICOUNT_ACT - 1 ) .gt.0)THEN
C
C             write ILEN2 into file array for LUOUT
C
              IVCOFF_OUT( IOFFSET_INT_OUT + IBLK - 1 ) = ILEN2
C
C             memory offset
C
              JOFF = NBATCH_INF(6,IBLK)
C
CSK          WRITE(LUWRT,*) ' JOFF for read in IBLK at offset',JOFF, IBLK,
CSK     &                 IOFFSET
C
C           read vector
C
              CALL MPI_FILE_READ_AT(LUIN,IOFFSET,SCR(JOFF),
     &                              ILEN2,MPI_REAL8,my_STATUS,IERR)
C
CSK              WRITE(LUWRT,'(2X,A,1X,I12,1X,A,1X,I6,1X,I6)') 
CSK     & 'Read-in from ILU1 at',IOFFSET,'for block',IOFFSET_INT_OUT +
CSK     &  IBLK - 1, JOFF
CSK              CALL WRTMATMN(SCR(JOFF),1,ILEN2,1,ILEN2,LUWRT)
C
            END IF 
C           ^ IVCOFF_IN .gt. 0 ?
          END IF
C         ^ ILEN2 .gt. 0 ?
C
          IOFFSET = IOFFSET + ILEN2
C
C
        ENDDO
C       ^ loop over blocks in batch
C
      END
***********************************************************************

      SUBROUTINE WTVEC_BATCH_DRV3(SCR,NBATCH_BLK,NBATCH_INF,
     &                            LUOUT,IOFFSET,IOFFSET_INT,
     &                            IVCOFF_OUT,NO_CHECK,ISCRATCH_SP,
     &                            IBLOCKL)
C
C     Written by  S. Knecht         -      May 21 2007
C
C**********************************************************************
C
C
C     write batch from SCR to MPI-file LUOUT - C-vector routine
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_OUT(*)
      DIMENSION IBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET, ISCRATCH_SP
      INTEGER JOFF, LENGTH
C
      JOFF   = 0
      LENGTH = 0
CSK      WRITE(LUWRT,*) 'calculating length start block is', IOFFSET_INT
C
C     loop over all blocks in that batch, do a norm check
C
      DO IBLK = 1, NBATCH_BLK
C
        JOFF = NBATCH_INF(6,IBLK)
csk     WRITE(LUWRT,*) 'calculating length for block',
csk  &                  IOFFSET_INT + IBLK - 1
        LENGTH = IBLOCKL( IOFFSET_INT + IBLK - 1 )
C       check norm of vector
        IF( NO_CHECK .eq. 0 ) THEN
          XXX = 0.0D0
          XXX = DDOT(LENGTH,SCR(JOFF),1,SCR(JOFF),1)
C
          IF( XXX .eq. 0.0D0 ) THEN
            IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) = 0
            GOTO 100
          END IF
        ELSE
C
          IF( IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) .eq. 0 ) GOTO 100
C
        END IF
C
csk     WRITE(LUWRT,'(2X,A,1X,I6,1X,A,1X,I14,1X,I6)') 'THIS IS block',
csk  &       IOFFSET_INT + IBLK - 1,'to go on pos (JOFF)',IOFFSET,JOFF
csk       CALL WRTMATMN(SCR(JOFF),1,LENGTH,1,LENGTH,LUWRT)
C
          CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET,SCR(JOFF),
     &                           LENGTH,MPI_REAL8,my_STATUS,IERR)
C
C
 100    CONTINUE  ! skip zero blocks on file
        IOFFSET = IOFFSET + LENGTH
        ISCRATCH_SP = ISCRATCH_SP + LENGTH
C
      END DO
C     ^ loop over blocks in a batch
C
      END

***********************************************************************

      SUBROUTINE DISTBLKND_1(NDIM,ICWEIGHTF,NPARBLOCK,NBLKWT,NPARBLKWT,
     &                       NVAR,ICCTOS,IBLOCKL,IPROCLIST,ICWEIGHT,
     &                       IABSOLUTE_WEIGHT)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
C     INPUT
      INTEGER ICOUNTABLK, NTEST,MWEIGHT
      INTEGER*8 ITOTBLCKL,MMPINFA
      DIMENSION ICWEIGHTF(NDIM),ICCTOS(NDIM,NDIM), IBLOCKL(NDIM)
      DIMENSION NPARBLOCK(NDIM),NBLKWT(NDIM),MMPINFA(LUCI_NMPROC)
      DIMENSION MWEIGHT(LUCI_NMPROC),IPROCLIST(LUCI_NMPROC)
      DIMENSION ICWEIGHT(NDIM)
C     scratch
C     MXSIZE: actual maximum size for one of all blocks
C     MXNUMB: number of the current largest block
C     MMPINFA(*,1): total number of blocks for the n-th proc
C     MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2,ILEN
      INTEGER*8 LBPROC,MXSZTMP2,MXSIZE2,ITEMPL2,IABSOLUTE_WEIGHT
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
      LOGICAL SUCCESS
C     initializtion
      SUCCESS = .FALSE.
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 0
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 0
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO I = 1, LUCI_NMPROC
        MMPINFA(I) = 0
        MWEIGHT(I) = 0
      END DO
C
C     count total number of active blocks
C     set need for non-zero block to -1
C
      DO ICBL = 1, NDIM
C
        ILEN = IBLOCKL(ICBL)
        IF(ILEN.NE.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
           ITOTBLCKL = ITOTBLCKL + ILEN
           NPARBLOCK(ICBL) = -1
        END IF
        IF(ILEN.GE.MXSIZE) MXSIZE = ILEN
C
      END DO
C
      ITEMPW1 = MXSIZE
C
      IF(NTEST.GT.0) THEN
        WRITE(luwrt,*) '  total number of processes to distribute on:',
     &                    LUCI_NMPROC
        WRITE(luwrt,*) '  total number of active blocks:', ICOUNTABLK
        WRITE(luwrt,*) '  overall active block length:', ITOTBLCKL
      END IF
*
      
      IAM_NOT_INV =   1
      IMINNP = LUCI_NMPROC
      NTEMPP = IMINNP
      IF(ICOUNTABLK.LT.LUCI_NMPROC) THEN
        IMINNP = MIN(ICOUNTABLK,IMINNP)
        write(luwrt,'(/a)') '  INFO from distblknd_1:'
        write(luwrt,'( a)')
     &    '  INFO: # of active blocks < total # of processes.' 
        write(luwrt,'(a,i6)') '  INFO: consider to decrease the'//
     &    ' # of processes to:',IMINNP
        write(luwrt,'(a)') '  INFO: alternatively change the GAS'//
     &    ' specification (more GA spaces).'
        write(luwrt,'(a/)') '  INFO: this ensures an increase'//
     &    ' in the number of active TTSS blocks.'
        NTEMPP = IMINNP
        IF(LUCI_MYPROC.GE.ICOUNTABLK) THEN
          IAM_NOT_INV =   0
          IMINUS2     = - 2
          CALL ISETVC(NPARBLOCK,IMINUS2,NDIM)
          GOTO 1001
        END IF
      END IF
      IAM_NOT_INV =   1
C
C     maximum number of determinants on a given cpu
C
      LBPROC = (( NVAR -1 ) / (LUCI_NMPROC)) + 1
C?      LBPROC = (( IABSOLUTE_WEIGHT  - 1 ) / (NTEMPP)) + 1
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'LBPROC is',LBPROC
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'NTEMPP is',NTEMPP
C
Csk?      IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
Csk?        CALL IWRTMA(NPARBLOCK,1,NDIM,1,NDIM)
Csk?      END IF
C
C     find optimal (?) c-block distribution 
C
C
C     start with the 'master group' and with the master
      IGROUP = 1
      IPROC  = 1
      ICPROC = 0
      IASSIGNED_BLK = 0
C
 100  CONTINUE
C
      MXSZTMP = 0
      ICPROC = ICPROC + 1
      IF( ICPROC .gt. NTEMPP ) GOTO 500
      IF( IASSIGNED_BLK .ge. ICOUNTABLK ) GOTO 1001
C     find largest unassigned block:
C     MXNUMB  = block number
C     MXSZTMP = block 'weight'
C
      DO II = 1, NDIM
C
C        ITEMPL = NBLKWT(II)
        ITEMPL = IBLOCKL(II)
        ITEMPD = NPARBLOCK(II)
        ITEMPN = II
C?        IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'ITEMPL for II',ITEMPL,II
C
        IF(ITEMPL.GT.0)THEN
          IF(ITEMPD.EQ.-1) THEN
            IF(ITEMPL.LE.MXSIZE)THEN
              IF(ITEMPL.GT.MXSZTMP) THEN
C?                 IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'MXSZTMP,ITEMPL',
C?     & MXSZTMP,ITEMPL,MXSIZE
                 MXSZTMP = ITEMPL
                 MXNUMB = ITEMPN
              END IF
            END IF
          END IF
        END IF
      END DO
C
 150  CONTINUE
C
C?        IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'assigned in 150',MXNUMB,IPROC
        MMPINFA(IPROC) = MMPINFA(IPROC) + MXSZTMP
        MWEIGHT(IPROC) = MWEIGHT(IPROC) + 1
        NPARBLOCK(MXNUMB) = IPROC - 1
        MXSIZE = MXSZTMP
        IASSIGNED_BLK = IASSIGNED_BLK + 1
C
C     assign all blocks that are connected to MXNUMB to a given CPU
C     as long as MMPINFA <= LBPROC !
C
      SUCCESS = .FALSE.
C
 200  CONTINUE 
C
C     find largest block connected to MXNUMB
C
C     1. find all connections
      JJJ = 0
      DO IJ = 1, NDIM
        IF(ICCTOS(IJ,MXNUMB).NE.0) THEN
           IF( NPARBLOCK(IJ) .eq. -1 ) THEN
             JJJ = JJJ + 1
             ICWEIGHT(JJJ) = IJ
           END IF
        ENDIF
      END DO
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'JJJ is',JJJ
C
      IRUN = 0
 250  CONTINUE
      MXSIZE2 = MXSIZE
      MXSZTMP = 0
      IF( JJJ .eq. 0 ) GOTO 290
C
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'current IPROC is',IPROC
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) THEN
C?        WRITE(6,*)'current ICWEIGHT'
C?        CALL IWRTMA(ICWEIGHT,1,JJJ,1,JJJ)
C?      END IF
C
C?      IRUN = 0
 275  CONTINUE
      MXSZTMP = 0
C     2. select the largest block not yet assigned
      DO IASB = 1, JJJ
C
        JBLK = ICWEIGHT(IASB)
        ITEMPL2 = IBLOCKL(JBLK)
C?       ITEMPL2 = NBLKWT(JBLK)
        ITEMPD2 = NPARBLOCK(JBLK)
        ITEMPN2 = JBLK
C
        IF(ITEMPD2.EQ.-1) THEN
          IF(ITEMPL2.LE.MXSIZE2)THEN
            IF(ITEMPL2.GT.MXSZTMP) THEN
               MXSZTMP  = ITEMPL2
               MXNUMB2  = ITEMPN2
               MXSZTMP2 = ITEMPL2
            END IF
          END IF
        END IF
      END DO
C
C     assign
C
      IF( ( MMPINFA(IPROC) + MXSZTMP2 ) .le. LBPROC ) THEN
C?        IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'assigned in 250',MXNUMB2,IPROC
        MMPINFA(IPROC) = MMPINFA(IPROC) + MXSZTMP2
        MWEIGHT(IPROC) = MWEIGHT(IPROC) + 1
        NPARBLOCK(MXNUMB2) = IPROC - 1
        MXSIZE2 = MXSZTMP2
        IASSIGNED_BLK = IASSIGNED_BLK + 1
      ELSE
C
C       CPU in same group?
 280    CONTINUE
C
        SUCCESS = .FALSE.
C
        DO ITCPU = 1, NTEMPP
C
          IF(.NOT.SUCCESS) THEN
            JGROUP = IPROCLIST(ITCPU)
            IF( JGROUP .eq. IGROUP ) THEN
              IF( ITCPU .gt. IPROC ) THEN
                IPROC = ITCPU
                SUCCESS = .TRUE.
              END IF
            END IF
          END IF
C
        END DO
C
        IF( .NOT. SUCCESS ) THEN
          IGROUP = IGROUP + 1
C?          IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'IGROUP,NFLGRPS',IGROUP,NFLGRPS
          IF( IGROUP .gt. NFLGRPS ) GOTO 500
C         find the lowest CPU number in the new group
          IPROC = 1
          GOTO 280
        ELSE
C         fresh cpu
          GOTO 100
        END IF
C
      END IF
      IRUN = IRUN + 1
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'IRUN IS ',IRUN
      IF( IRUN .lt. JJJ ) THEN
C?      IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'I WILL GO TO 250'
        GOTO 250
      ELSE 
        GOTO 290
      END IF
C       CPU in same group?
 290    CONTINUE
C
        SUCCESS = .FALSE.
C
        DO ITCPU = 1, NTEMPP
C
          IF(.NOT.SUCCESS) THEN
            JGROUP = IPROCLIST(ITCPU)
            IF( JGROUP .eq. IGROUP ) THEN
              IF( ITCPU .gt. IPROC ) THEN
                IPROC = ITCPU
                SUCCESS = .TRUE.
              END IF
            END IF
          END IF
C
        END DO
C
        IF( .NOT. SUCCESS ) THEN
          IGROUP = IGROUP + 1
C?          IF(LUCI_MYPROC.eq.LUCI_MASTER) WRITE(6,*)'2. IGROUP,NFLGRPS',
C?     &                         IGROUP,NFLGRPS
          IF( IGROUP .gt. NFLGRPS ) GOTO 500
C         find the lowest CPU number in the new group
          IPROC = 1
          GOTO 290
        END IF
C       fresh cpu
        GOTO 100
C
 500  CONTINUE
C
C     check for 'lonely' blocks
C
      LONBLK = 0
C
      DO IBLK = 1, NDIM
        IF( NPARBLOCK(IBLK) .eq. -1 ) THEN 
          LONBLK = LONBLK + 1
          ICWEIGHT(LONBLK) = IBLK
        END IF
      END DO
C?      IF (LUCI_MYPROC .eq. LUCI_MASTER) WRITE(6,*) 'LONBLK is ',LONBLK
C
      IF( LONBLK .eq. 0 ) GOTO 1001
C
C     good old normal block distribution
C
      IRUN = 0
      MXSIZE = ITEMPW1
C
 600  CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.LONBLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, LONBLK
*
          IBLK   = ICWEIGHT(II) 
          ITEMPL = IBLOCKL(IBLK)
          ITEMPD = NPARBLOCK(IBLK)
          ITEMPN = IBLK
C
C?          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
                END IF
              END IF
            END IF
C?          END IF
*
3000    CONTINUE
*
        DO 4000 IPR = 1, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR) 
CSK         IF(LUCI_MYPROC.EQ.0) 
CSK     & WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MYNEW_ID = 0 !
          IAMTPROC = IPR
CSK         IF(LUCI_MYPROC.EQ.0) WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
CSK         IF(LUCI_MYPROC.EQ.0) WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
C      now we should have found a proc and a block --> put both together !
C
C      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
C?       IF(NTEST.GT.0) THEN
C?         IF(LUCI_MYPROC.EQ.0) 
C?     & WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
C?       ENDIF
C
C      raising MMPINFA(proc,1) by 1
C      adding on MMPINFA(proc,2) the weight of the new block
C
       IAMTPROC = IAMGPROC-1
       MXSIZE = MXSZTMP
       MWEIGHT(IAMGPROC) = MWEIGHT(IAMGPROC) + 1
       MMPINFA(IAMGPROC) = MMPINFA(IAMGPROC) + MXSZTMP
       NPARBLOCK(MXNUMB)     = IAMTPROC
C
      ELSE
        GOTO 1001
      END IF
C     /\ IRUN !!!
C
      GOTO 600
C
1001  CONTINUE
CSK   IF(NTEST.GT.0) THEN
      WRITE(luwrt,'(3X,A,I4)')
     &'SUMMATION OF EVEN DISTRIBUTION OF LUCI_MYPROC:',LUCI_MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(luwrt,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MWEIGHT(ISTI),' blocks with a
     & total length of',MMPINFA(ISTI)
      END DO
CSK   END IF

      END
***********************************************************************

      SUBROUTINE DISTBLKND_2(NDIM,ICWEIGHTF,NPARBLOCK,IBLOCKL,
     &                       JSYM_DISTBLK)
*
*. Find distribution for (nonvanishing blocks) among the nodes
*. Non vanishing block has a nonzero blocklength
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
C     -----
C     INPUT
C     -----
C     NDIM: total number of c-blocks
C     IBLOCKL(NDIM): total length of each block
C     ICWEIGHTF(NDIM): total 'weight factors' for each block
C     JSYM_DISTBLK: symmetry irrep of right-hand side vector
C
C     ------
C     OUTPUT
C     ------
*     NPARBLOCK(NDIM): list of blocks with the corresponding assigned node  
*
      INTEGER   ICOUNTABLK, NTEST
      DIMENSION ICWEIGHTF(NDIM)
      DIMENSION NPARBLOCK(NDIM),IBLOCKL(NDIM)
C     -------
C     SCRATCH
C     -------
C     MXSIZE: actual maximum size for one of all blocks
C     MXNUMB: number of the current largest block
C     MMPINFA(*,1): total number of blocks for the n-th proc
C     MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
      INTEGER(KIND=MPI_OFFSET_KIND) ITLTND,MXSIZE,MXSZTMP
      INTEGER(KIND=MPI_OFFSET_KIND) ITEMPL, ITEMPW1, ITEMPW2
      INTEGER(KIND=MPI_OFFSET_KIND) ILEN, ITOTBLCKL, IILENGTH_MAX 
      INTEGER(KIND=MPI_OFFSET_KIND) IIWEIGHTBLK, IILENGTH 
*
      NTEST        = 0

      IIWEIGHTBLK  = 0
      IILENGTH     = 0
      ITLTND       = 0
      ICOUNTABLK   = 0
      MXSIZE       = 0
      MXSZTMP      = 0
      MXNUMB       = 0
      IRUN         = 0
      ITEMPW1      = 0
      ITEMPW2      = 0
      ITEMPN       = 0
      ITEMPL       = 0
      IAMGPROC     = 0
      IAMTPROC     = 0
      ITOTBLCKL    = 0
      ILEN         = 0
      IILENGTH_MAX = 0

*. end of initialization
      DO ICBL = 1, NDIM
C
        IIWEIGHTBLK = ICWEIGHTF(ICBL)
        IILENGTH    = IBLOCKL(ICBL) 
        ILEN        = IILENGTH * IIWEIGHTBLK
C
        IF( ILEN .lt. 0 ) THEN 
          WRITE(luwrt,*)'Attention, minus block detected',ICBL
          WRITE(luwrt,*)'ILEN = ',ILEN
          WRITE(luwrt,*)'IBLOCKL(ICBL) = ',IBLOCKL(ICBL)
          WRITE(luwrt,*)'ICWEIGHTF(ICBL) = ',ICWEIGHTF(ICBL)
        END IF

        IF(ILEN.ne.0) THEN
          ICOUNTABLK = ICOUNTABLK + 1
          ITOTBLCKL = ITOTBLCKL + ILEN
          NPARBLOCK(ICBL) = -1
        END IF

        MXSIZE       = max(MXSIZE,ILEN)
        IILENGTH_MAX = max(IILENGTH_MAX,IILENGTH)
      END DO
      ITEMPW1 = MXSIZE
!
!     print header for statistics
!
      WRITE(LUWRT,'(/7X,A)')
     &  '==================================================='
      WRITE(LUWRT,'(7X, A,I2)')
     &  ' parallel distribution setup for symmetry irrep ',JSYM_DISTBLK
      WRITE(LUWRT,'(7X,A/)')
     &  '==================================================='
      WRITE(LUWRT,'(A,I19)')
     &'  total number of processes to distribute on :',LUCI_NMPROC
      WRITE(LUWRT,'(A,I19)')
     &'  total number of blocks                     :',NDIM
      WRITE(LUWRT,'(A,I19)')
     &'  total number of active blocks              :',ICOUNTABLK
      WRITE(LUWRT,'(A,I19)')
     &'  size of largest TTSS block                 :',IILENGTH_MAX
      WRITE(LUWRT,'(A,I19)')
     &'  overall weighted active block length       :',ITOTBLCKL
      WRITE(LUWRT,'(A,I19)')
     &'  Maximum weighted block size                :',MXSIZE
*
      IAM_NOT_INV = 1
      IMINNP = LUCI_NMPROC
      NTEMPP = IMINNP
      IF( ICOUNTABLK .lt. LUCI_NMPROC) THEN
        IMINNP = MIN( ICOUNTABLK, IMINNP )
        write(luwrt,'(/a)') '  INFO from distblknd_2:'
        write(luwrt,'( a)')
     &    '  INFO: # of active blocks < total # of processes.' 
        write(luwrt,'(a,i6)') '  INFO: consider to decrease the'//
     &    ' # of processes to:',IMINNP
        write(luwrt,'(a)') '  INFO: alternatively change the GAS'//
     &    ' specification (more GA spaces).'
        write(luwrt,'(a/)') '  INFO: this ensures an increase'//
     &    ' in the number of active TTSS blocks.'
        NTEMPP = IMINNP
        IF( LUCI_MYPROC .ge. ICOUNTABLK ) THEN
          IAM_NOT_INV =   0
          IMINUS2     = - 2
          CALL ISETVC(NPARBLOCK,IMINUS2,NDIM)
          GOTO 101
        END IF
      END IF
      IAM_NOT_INV = 1
*
*. starting the treausure quest for the ?optimal? c-block distribution
*
100   CONTINUE
*
      IRUN = IRUN + 1
      IF(IRUN.LE.ICOUNTABLK) THEN
*
        MXSZTMP = 0
*
        DO 3000 II = 1, NDIM
*
          
          IIWEIGHTBLK = ICWEIGHTF(II)
          IILENGTH = IBLOCKL(II)
C
          ITEMPL = IILENGTH * IIWEIGHTBLK
C
          ITEMPD = NPARBLOCK(II)
          ITEMPN = II

CSK          IF(LUCI_MYPROC.EQ.0) 
CSK     & WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
CSK          IF(LUCI_MYPROC.EQ.0) WRITE(6,*) 'MXSIZE:',MXSIZE
CSK          IF(LUCI_MYPROC.EQ.0) WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.EQ.-1) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
CSK                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
CSK       IF(LUCI_MYPROC.EQ.0) WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 1, NTEMPP
*
C?          ITEMPW2 = MMPINFA(IPR) 
           ITEMPW2 = 0
           DO IBLK = 1, NDIM
            IF( NPARBLOCK(IBLK) .eq. IPR-1 ) THEN
              IIWEIGHTBLK = ICWEIGHTF(IBLK)
              IILENGTH = IBLOCKL(IBLK)
              ITEMPW2 = ITEMPW2 + ( IILENGTH * IIWEIGHTBLK )
            END IF
           END DO
          
CSK         IF(LUCI_MYPROC.EQ.0) 
CSK     & WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MYNEW_ID = 0 !
          IAMTPROC = IPR
CSK         IF(LUCI_MYPROC.EQ.0) WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
CSK         IF(LUCI_MYPROC.EQ.0) WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.     calculation of block MXNUMB by proc IAMGPROC -1
       IF(NTEST.GT.0) THEN
         WRITE(luwrt,*)'calculation of block MXNUMB by proc',
     &                  MXNUMB,IAMGPROC-1
       ENDIF
*
*.     raising MMPINFA(proc,1) by 1
*.     adding on MMPINFA(proc,2) the weight of the new block
       IAMTPROC = IAMGPROC-1
       MXSIZE = MXSZTMP
C?       MWEIGHT(IAMGPROC) = MWEIGHT(IAMGPROC) + 1
C?       MMPINFA(IAMGPROC) = MMPINFA(IAMGPROC) + MXSZTMP
       NPARBLOCK(MXNUMB)     = IAMTPROC
*
      ELSE
        GOTO 101
      END IF
*     /\ IRUN !!!
*
      GOTO 100
*
101   CONTINUE
*
      DO II = 1, NDIM
        IF(NPARBLOCK(II) .eq. -1)THEN
          WRITE(luwrt,*) 'AAAHHHH, block',II,' is not distributed'
          WRITE(luwrt,*) 'I do not know how to proceed, I will stop!'
          call quit('NE DISTBLKND_2: Block not distributed')
        END IF
      END DO

      WRITE(luwrt,'(/18X,A)')'================================'
      WRITE(luwrt,'( 18X,A)')' summation of even distribution '
      WRITE(luwrt,'(18X,A/)')'================================'

      DO ISTI = 1, NTEMPP
        ITLTND = 0
        ITLBND = 0
        DO IBLK = 1, NDIM
          IF( NPARBLOCK(IBLK) .eq. ISTI -1 ) THEN
C
             IIWEIGHTBLK = ICWEIGHTF(IBLK)
             IILENGTH = IBLOCKL(IBLK)
C
             ITLTND = ITLTND + ( IILENGTH * IIWEIGHTBLK )
             ITLBND = ITLBND + 1
          END IF
        END DO
        WRITE(luwrt,'(2X,A,I4,A,I5,1X,A,I17)')
     &  'CPU',ISTI-1,' computes',ITLBND,'blocks: total weight ==>',
     &  ITLTND
        call flshfo(luwrt)
      END DO
      
*
      END
***********************************************************************

      SUBROUTINE FIND_CCTOS(ISCALFAC_LOCAL,IBLOCKD,ICCTOS,NBLOCK)
*
*. Find all c-blocks connecting to a given sigma-block for each cpu
*  using the connection matrix ICCTOS.
*. Each cpu stores the information in ISCALFAC_LOCAL.
*. Connection is marked by 1.
*
      IMPLICIT REAL*8           ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
*     ==========
*       INPUT
*     ==========
      DIMENSION ICCTOS(NBLOCK,NBLOCK),IBLOCKD(NBLOCK)
*
*     ===========
*       OUTPUT
*     ===========
*
      DIMENSION ISCALFAC_LOCAL(NBLOCK)
*
      IONE = 1
*
      DO IBLK = 1, NBLOCK
*
        IF(IBLOCKD(IBLK).EQ.LUCI_MYPROC) THEN
*
          DO JBLK = 1,NBLOCK
*
            IF(ICCTOS(JBLK,IBLK).NE.0) ISCALFAC_LOCAL(JBLK) = IONE
*
          END DO
*
        END IF
*
      END DO
*
      END
***********************************************************************

      SUBROUTINE ISET_ARRAY_ACT_BLK(IAM_BLK_ACT,IBLOCKL,NDIM)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IBLOCKL(*), IAM_BLK_ACT(*)
C
        ICOUNT = 0
C
        DO IBLK = 1, NDIM
          IF( IBLOCKL(IBLK) .gt. 0 ) THEN
            ICOUNT = ICOUNT + 1
            IAM_BLK_ACT( IBLK ) = ICOUNT
          END IF
        END DO
C
      END 
***********************************************************************

      SUBROUTINE INFO_PRINT_BATCH_REL(LBATCH, LEBATCH,
     &                                I1BATCH,IBATCH,JBATCH,LUPRINT)
      IMPLICIT REAL*8(A-H,O-Z)
C     Input
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
C
C
C     
      DO ISBATCH = 1, JBATCH
        WRITE(LUPRINT,*)
        WRITE(LUPRINT,*) ' Info on batch ', ISBATCH
        WRITE(LUPRINT,*) ' *********************** '
        WRITE(LUPRINT,*)
        WRITE(LUPRINT,*) '      Number of blocks included ',
     &                          LBATCH(ISBATCH)
        WRITE(LUPRINT,*) '      Length of batch           ',
     &                          LEBATCH(ISBATCH)
        WRITE(LUPRINT,*) '      TTSS and offsets and lenghts of
     &                          each block '
        DO IBLOCK = I1BATCH(ISBATCH),I1BATCH(ISBATCH)+ LBATCH(ISBATCH)-1
          WRITE(LUPRINT,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
        END DO
      END DO
      END
***********************************************************************

      SUBROUTINE COPVCD_PP_CC_B_C(LUIN,LUOUT,SCR,NBATCH,LBATCH,LEBATCH,
     &                            I1BATCH,IBATCH,MY_IOFF_LUIN,
     &                            MY_IOFF_LUOUT,
     &                            LUINLIST,LUOUTLIST,IBLOCKL,JOFF)
C
C     Written by  S. Knecht         - May 18 2007
C
C**********************************************************************
C
C     copy c-vector from file LUIN to LUOUT batchwise
C     update the file lists.
C
C     NOTE: JOFF = (IVEC resp. IROOT) - 1
C
C
C     copy vector from file LUIN to file LUOUT --> (hint LUIN, LUOUT)
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     ready for complex double groups - S. Knecht, Odense, 22 June 2007
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*),LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*), IBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN, IOFFSET_OUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER(KIND=MPI_OFFSET_KIND) ISCRATCH_SP
      INTEGER(KIND=MPI_OFFSET_KIND) ILEN_COMB
      INTEGER NUM_BLK
C     some constants
      IONE = 1
      ITWO = 2
C     initialize scratch offsets
      NUM_BLK = 0
      MY_NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN  = 0
      IOFFSET_OUT  = 0
      IOFFSET_INT_IN  = 0
      IOFFSET_INT_OUT  = 0
      ILEN_COMB = 0
      ILEN_COMB = L_COMBI
      MY_IOFFSET_SCRATCH = 0
      ISCRATCH_SP = 0
      NUM_BLK_SP = 0
C
C     ================
C      COMPLEX VECTOR
C     ================
C
C     .................
C      REAL PART FIRST
C     .................
C
C     set new offset
C
C     position in file is at the end of vector JOFF
C
C     note: real part: --> MY_VEC2_IOFF, MY_ACT_BLK2
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(SCR,LEBATCH(ISBATCH))
C
        IOFFSET_IN = MY_IOFF_LUIN + ( MY_VEC2_IOFF * JOFF ) +
     &               MY_IOFFSET_SCRATCH
        IOFFSET_INT_IN = ( MY_ACT_BLK2 * JOFF ) + MY_NUM_BLK + 1
C
C       write active block array for LUOUT 
C
        IOFFSET_INT_OUT = NUM_BLK + 1
C
        NUM_BLK_SP  = 0
        ISCRATCH_SP = 0
C
CSK        WRITE(LUWRT,*) ' THIS IS MY IOFFSET_IN',IOFFSET_IN
CSK        WRITE(LUWRT,*) ' THIS IS MY IOFFSET_INT_IN',IOFFSET_INT_IN
CSK        WRITE(LUWRT,*) ' THIS IS MY IOFFSET_INT_OUT',IOFFSET_INT_OUT
C
        CALL RDVEC_BATCH_DRV3(LUIN,SCR,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN,IOFFSET_INT_IN,
     &                       IOFFSET_INT_OUT,LUINLIST,
     &                       LUOUTLIST,NUM_BLK_SP)
C
        IOFFSET_OUT = MY_IOFF_LUOUT + IOFFSET_SCRATCH
     
CSK        WRITE(LUWRT,*) ' THIS IS MY IOFFSET_OUT',IOFFSET_OUT
CSK        CALL WRTMATMN(SCR,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
CSK        WRITE(LUWRT,*) ' THIS IS MY IOFFSET_OUT',IOFFSET_OUT
        
        IOFFSET_INT_OUT = NUM_BLK + 1
C
        CALL WTVEC_BATCH_DRV3(SCR,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        LUOUT,IOFFSET_OUT,IOFFSET_INT_OUT,
     &                        LUOUTLIST,1,ISCRATCH_SP,
     &                        IBLOCKL)
C
C
C       count the length of the last copy
C
C       LUOUT
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + ISCRATCH_SP
        NUM_BLK         = NUM_BLK + LBATCH(ISBATCH)
C       LUIN
        MY_NUM_BLK         = MY_NUM_BLK + NUM_BLK_SP
        MY_IOFFSET_SCRATCH = MY_IOFFSET_SCRATCH + LEBATCH(ISBATCH)
C
      END DO
C
      IOFFSET_SCRATCH    = 0
      MY_IOFFSET_SCRATCH = 0
      NUM_BLK    = 0
      MY_NUM_BLK = 0
      NUM_BLK_SP = 0
      ISCRATCH_SP = 0
C
C     ..............
C      COMPLEX PART
C     ..............
C
C     set new offset
C
C     position in file is at the end of vector JOFF
C
C     note: complex part: --> MY_VEC1_IOFF and MY_ACT_BLK1
C
      DO JSBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(SCR,LEBATCH(JSBATCH))
C
        IOFFSET_IN = MY_IOFF_LUIN + ( MY_VEC2_IOFF * JOFF ) +
     &               MY_IOFFSET_SCRATCH +  MY_VEC1_IOFF 
C
        IOFFSET_INT_IN = ( MY_ACT_BLK2 * JOFF ) + MY_NUM_BLK + 
     &                   1 + MY_ACT_BLK1
C
C       write active block array for LUOUT - complex part 
C       --> MY_ACT_BLK_ALL needed!
C
        IOFFSET_INT_OUT = MY_ACT_BLK_ALL + NUM_BLK + 1
C
        NUM_BLK_SP  = 0
        ISCRATCH_SP = 0
C
CSK      WRITE(LUWRT,*) ' THIS IS MY IOFFSET_IN imag',IOFFSET_IN
CSK      WRITE(LUWRT,*) ' THIS IS MY IOFFSET_INT_IN imag',IOFFSET_INT_IN
CSK      WRITE(LUWRT,*) ' THIS IS MY IOFFSET_INT_OUT imag',IOFFSET_INT_OUT
C
        CALL RDVEC_BATCH_DRV3(LUIN,SCR,LBATCH(JSBATCH),
     &                        IBATCH(1,I1BATCH(JSBATCH)),
     &                        IOFFSET_IN,IOFFSET_INT_IN,
     &                        IOFFSET_INT_OUT,LUINLIST,
     &                        LUOUTLIST,NUM_BLK_SP)
C
        IOFFSET_OUT = MY_IOFF_LUOUT + ILEN_COMB + IOFFSET_SCRATCH
C     
CSK        WRITE(LUWRT,*) ' THIS IS MY IOFFSET_OUT imag',IOFFSET_OUT
CSK        CALL WRTMATMN(SCR,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
      
        IOFFSET_INT_OUT = MY_ACT_BLK_ALL + NUM_BLK + 1
C
        CALL WTVEC_BATCH_DRV3_C(SCR,LBATCH(JSBATCH),
     &                          IBATCH(1,I1BATCH(JSBATCH)),
     &                          LUOUT,IOFFSET_OUT,IOFFSET_INT_OUT,
     &                          LUOUTLIST,1,ISCRATCH_SP,
     &                          IBLOCKL)
C
C
C       count the length of the last copy
C
C       LUOUT
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + ISCRATCH_SP
        NUM_BLK         = NUM_BLK + LBATCH(JSBATCH)
C       LUIN
        MY_NUM_BLK         = MY_NUM_BLK + NUM_BLK_SP
        MY_IOFFSET_SCRATCH = MY_IOFFSET_SCRATCH + LEBATCH(JSBATCH)
C
      END DO
C
C
C
      END
***********************************************************************

      SUBROUTINE COPVCD_PP_B_RL(VEC1,LUINLIST,LUOUTLIST,
     &                          NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                          MY_IOFF_LUIN,MY_IOFF_LUOUT,
     &                          IROOT,LUIN,LUOUT)
C
C     Written by  S. Knecht         - June 11 2007
C
C**********************************************************************
C
C     copy vector from file LUIN to LUOUT batchwise
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_LUOUT = 0
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IROOT - 1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )   * MY_VEC1_IOFF
        IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                   ( IROOT - 1 )   * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                LUWRT)
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    +
     &                        ( IROOT - 1 ) * MY_ACT_BLK1
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE COPVCD_PP_B_CPX(VEC1,LUINLIST,LUOUTLIST,
     &                           NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                           MY_IOFF_LUIN,MY_IOFF_LUOUT,
     &                           IROOT,LUIN,LUOUT)
C
C     Written by  S. Knecht         - June 26 2007
C
C**********************************************************************
C
C     copy vector from file LUIN to LUOUT batchwise
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_LUOUT = 0
C
C
C     ==============
C     COMPLEX VECTOR
C     ==============
C
C     ---------
C     REAL PART
C     ---------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IROOT - 1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )   * MY_VEC2_IOFF
        IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                   ( IROOT - 1 )   * MY_ACT_BLK2
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                LUWRT)
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC2_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    +
     &                        ( IROOT - 1 ) * MY_ACT_BLK2
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C
C     zero scratch offsets
      NUM_BLK          = 0
      IOFFSET_SCRATCH  = 0
      NUM_ACTIVE_BATCH = 0
C
C     ---------
C     IMAG PART
C     ---------
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IROOT - 1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )  * MY_VEC2_IOFF    +
     &                                    MY_VEC1_IOFF
C
        IOFFSET_INT_IN   = 1 + NUM_BLK + MY_ACT_BLK1 +
     &                   ( IROOT - 1 ) * MY_ACT_BLK2
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                   IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                LUWRT)
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC2_IOFF    +
     &                                        MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    + MY_ACT_BLK1 +
     &                        ( IROOT - 1 ) * MY_ACT_BLK2
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE WTVEC_BATCH_DRV3_C(SCR,NBATCH_BLK,NBATCH_INF,
     &                              LUOUT,IOFFSET,IOFFSET_INT,
     &                              IVCOFF_OUT,NO_CHECK,ISCRATCH_SP,
     &                              IBLOCKL)
C
C     Written by  S. Knecht         -      June 22 2007
C
C**********************************************************************
C
C
C     write batch from SCR to MPI-file LUOUT - C-vector routine
C
C     COMPLEX vector part
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_OUT(*)
      DIMENSION IBLOCKL(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET, ISCRATCH_SP
      INTEGER JOFF, LENGTH
C
      JOFF   = 0
      LENGTH = 0
C
C     loop over all blocks in that batch, do a norm check
C
      DO IBLK = 1, NBATCH_BLK
C
        JOFF = NBATCH_INF(6,IBLK)
        LENGTH = IBLOCKL( IOFFSET_INT + IBLK - 1 - MY_ACT_BLK_ALL )
C       check norm of vector
        IF( NO_CHECK .eq. 0 ) THEN
          XXX = 0.0D0
          XXX = DDOT(LENGTH,SCR(JOFF),1,SCR(JOFF),1)
C
          IF( XXX .eq. 0.0D0 ) THEN
            IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) = 0
            GOTO 100
          END IF
        ELSE
C
          IF( IVCOFF_OUT( IOFFSET_INT + IBLK - 1 ) .eq. 0 ) GOTO 100
C
        END IF
C
CSK          WRITE(LUWRT,'(2X,A,1X,I6,1X,A,1X,I14,1X,I6)') 'THIS IS block',
CSK     &         IOFFSET_INT + IBLK - 1,'to go on pos (JOFF)',IOFFSET,JOFF
CSK          CALL WRTMATMN(SCR(JOFF),1,LENGTH,1,LENGTH,LUWRT)
C
          CALL MPI_FILE_WRITE_AT(LUOUT,IOFFSET,SCR(JOFF),
     &                           LENGTH,MPI_REAL8,my_STATUS,IERR)
C
C
 100    CONTINUE  ! skip zero blocks on file
        IOFFSET = IOFFSET + LENGTH
        ISCRATCH_SP = ISCRATCH_SP + LENGTH
C
      END DO
C     ^ loop over blocks in a batch
C
      END
***********************************************************************

      SUBROUTINE RDVEC_BATCH_DRV4(MY_LUIN,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET,IOFFSET_INT,IVCOFF_IN1,
     &                            NUM_A_B)
C
C     read in a batch from a MPI-file MY_LUIN to segment SCR
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
      NUM_A_B = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
            NUM_A_B = NUM_A_B + 1 
C
            ILEN  = IVCOFF_IN1( IOFFSET_INT + NUM_A_B -1 )
C
            IF( ILEN .gt. 0 ) THEN
C             memory offset
              JOFF = NBATCH_INF(6,IBLK)
C
CSK              WRITE(LUWRT,*) 'JOFF,ILEN2,ILEN again for block 
CSK     & IOFFSET_INT+ NUM_A_B -1',JOFF,ILEN2,ILEN,IOFFSET_INT+ NUM_A_B -1
              CALL MPI_FILE_READ_AT(MY_LUIN,IOFFSET,SCR(JOFF),
     &                              ILEN2,MPI_REAL8,my_STATUS,IERR)
CSK              WRITE(LUWRT,*) ' my block'
CSK              CALL WRTMATMN(SCR(JOFF),1,ILEN2,1,ILEN2,LUWRT)
            ENDIF 
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************

      SUBROUTINE WTVEC_BATCH_DRV4(MY_LUOUT,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET,IOFFSET_INT,IVCOFF_IN1,
     &                            NUM_A_B)
C
C     write batch from segment SCR to MPI-file MY_LUOUT
C
C     active blocks on the MPI-file are flagged by a nonzero entry in 
C     file list
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
      NUM_A_B = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
C           memory offset
            JOFF  = NBATCH_INF(6,IBLK)
            NUM_A_B = NUM_A_B + 1
C           check norm of vector
            XXX = 0.0D0
            XXX = DDOT(ILEN2,SCR(JOFF),1,SCR(JOFF),1)
C
            IF( XXX .eq. 0.0D0 ) THEN
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 0
              GOTO 100
            ELSE
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 1
            END IF
C
CSK              WRITE(LUWRT,*) 'JOFF,ILEN2,ILEN again for block 
CSK     & IOFFSET_INT+ NUM_A_B -1',JOFF,ILEN2,ILEN,IOFFSET_INT+ NUM_A_B -1
              CALL MPI_FILE_WRITE_AT(MY_LUOUT,IOFFSET,SCR(JOFF),
     &                               ILEN2,MPI_REAL8,my_STATUS,IERR)
C
 100        CONTINUE
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************

      SUBROUTINE WTVEC_BATCH_DRV4SP(MY_LUOUT,SCR,NBATCH_BLK,NBATCH_INF,
     &                              IOFFSET,IOFFSET_INT,IVCOFF_IN1,
     &                              NUM_A_B,IDEBUGPRNT)
C
C     write batch from segment SCR to MPI-file MY_LUOUT
C
C     active blocks on the MPI-file are flagged by a nonzero entry in 
C     file list
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*), IVCOFF_IN1(*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
      JOFF = 0
C
      NUM_A_B = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
C           memory offset
            JOFF  = NBATCH_INF(6,IBLK)
            NUM_A_B = NUM_A_B + 1
C           check norm of vector
            XXX = 0.0D0
            XXX = DDOT(ILEN2,SCR(JOFF),1,SCR(JOFF),1)
C
            IF( XXX .eq. 0.0D0 ) THEN
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 0
              GOTO 100
            ELSE
              IVCOFF_IN1( IOFFSET_INT + NUM_A_B - 1 ) = 1
            END IF
C
            IF( IDEBUGPRNT .ne. 0 ) THEN
              WRITE(LUWRT,'(A,1X,I6,1X,I6,1X,I6,1X,A,1X,I14)') 'JOFF, 
     & ILEN2 for block',JOFF,ILEN2,IOFFSET_INT+ NUM_A_B -1,'at off',
     & IOFFSET
              CALL WRTMATMN(SCR(JOFF),1,ILEN2,1,ILEN2,LUWRT)
            END IF
            CALL MPI_FILE_WRITE_AT(MY_LUOUT,IOFFSET,SCR(JOFF),
     &                             ILEN2,MPI_REAL8,my_STATUS,IERR)
C
 100        CONTINUE
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************

      SUBROUTINE RDVEC_BATCH_DRV5(MY_LUIN,SCR,NBATCH_BLK,NBATCH_INF,
     &                            IOFFSET)
C
C     read in a batch from a MPI-file MY_LUIN to segment SCR
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SCR(*), NBATCH_INF(8,*)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
      INTEGER JOFF
C
      JOFF = 0
C
        DO IBLK = 1, NBATCH_BLK
C
          ILEN2 = NBATCH_INF(8,IBLK) 
CSK          WRITE(LUWRT,*) 'ILEN2 for block IBLK',
CSK     &                    ILEN2,IBLK
          IF( ILEN2 .gt. 0 ) THEN
C
C           memory offset
C
            JOFF = NBATCH_INF(6,IBLK)
C
            CALL MPI_FILE_READ_AT(MY_LUIN,IOFFSET,SCR(JOFF),
     &                            ILEN2,MPI_REAL8,my_STATUS,IERR)
C
          ENDIF 
C
          IOFFSET = IOFFSET + ILEN2
C
        END DO
C       ^ loop over blocks in batch
C
      END
***********************************************************************

      SUBROUTINE TRAVC_B_RL_DRV(VEC1,VEC2,XMAT,LUIN1LIST,
     &                          LUOUTLIST,NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                          IBATCH,MY_IOFF_LUIN1,MY_IOFF_LUOUT,
     &                          NVEC,NVEC2,LUIN1,LUOUT,IALL_LUIN)
C
C     Written by  S. Knecht         - March 6 2008
C
C**********************************************************************
C
C     transforming vectors so that they become the actual approx. to the 
C     eigenvectors
C
C     NOTE: NVEC  = NVEC
C     NOTE: NVEC2 = NVEC2
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), XMAT(NVEC,NVEC2)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
C
C     transform: LUIN1 --> LUOUT
C
      DO IVECOUT = 1, NVEC2
         CALL TRAVC_B_RL(VEC1,VEC2,XMAT(1,IVECOUT),LUIN1LIST,
     &                   LUOUTLIST,NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                   IBATCH,MY_IOFF_LUIN1,MY_IOFF_LUOUT,NVEC,
     &                   IVECOUT,LUIN1,LUOUT)
      END DO
C
C     copy back: LUOUT --> LUIN1
C
      CALL IZERO(LUIN1LIST,IALL_LUIN)
C
      DO IROOT = 1, NVEC2
         CALL COPVCD_PP_B_RL(VEC1,LUOUTLIST,LUIN1LIST,
     &                       NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                       MY_IOFF_LUOUT,MY_IOFF_LUIN1,IROOT,LUOUT,
     &                       LUIN1)
      END DO
      END
***********************************************************************

      SUBROUTINE TRAVC_B_RL(VEC1,VEC2,FAC,LUIN1LIST,
     &                      LUOUTLIST,NBATCH,LBATCH,LEBATCH,I1BATCH,
     &                      IBATCH,MY_IOFF_LUIN1,MY_IOFF_LUOUT,
     &                      NVEC,ISVEC,LUIN1,LUOUT)
C
C     Written by  S. Knecht         - March 6 2008
C
C**********************************************************************
C
C     transforming vectors so that they become the actual approx. to the 
C     eigenvectors
C
C     NOTE: NVEC  = NVEC
C     NOTE: ISVEC = ISVEC
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), FAC(NVEC)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_LUOUT = 0
C
      FACTOR = 0.0D0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC
C
          FACTOR = FAC( IVEC )
C
C         set new offset
C
C         position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )  * MY_VEC1_IOFF
          IOFFSET_INT_IN1  = 1 + NUM_BLK   +
     &                     ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK     &                   IOFFSET_IN_LUIN1
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK            CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                           LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK            WRITE(LUWRT,*) 'initial VEC1 on LUIN1'
CSK            CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK     &                    LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
          END IF
C
  100   CONTINUE
C
C       write VEC2 to LUOUT wrt offset at ISVEC - 1
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + IOFFSET_SCRATCH + 
     &                      ( MY_VEC1_IOFF * ( ISVEC - 1 ) )
        IOFFSET_INT_LUOUT = 1 + NUM_BLK 
     &                        + ( MY_ACT_BLK1 * ( ISVEC - 1 ) )
csk        WRITE(LUWRT,*) 'This is my OFFSET for LUOUT',
csk     &                  IOFFSET_LUOUT
csk     WRITE(LUWRT,*) 'final VEC2 to write on LUOUT'
csk     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
csk  &                LUWRT)
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE ORTHG_VEC_BATCH_LUCI1(VEC1,VEC2,SUBSPH,SCRRED,
     &                                 SCRRED2,NROOT,JOFF,JOFF2,
     &                                 NBATCH,LBATCH,LEBATCH,
     &                                 I1BATCH,IBATCH,
     &                                 MY_IOFF_LUIN,MY_IOFF_LUIN2,
     &                                 MY_IOFF_LUIN3,MY_IOFF_LUOUT,
     &                                 LUINLIST,LUIN2LIST,
     &                                 LUIN3LIST,LUOUTLIST,
     &                                 LUIN,LUIN2,
     &                                 LUIN3,LUOUT)
C
C     calculate dot products between vector on file LUIN1 and LUIN2 
C     resp. LUIN1 and LUIN3. Store vector VEC2 on LUOUT 
C                            ==> LUIN1(VEC2) --> LUOUT(VEC2)
C
C     dot products SUM_batches ( VEC1 * VEC2 ) are reduced on every 
C     process in MPI_COMM_WORLD
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: JOFF  = JOFF
C           JOFF2 = JOFF2
C
C     Last revision:     S. Knecht       - March 2008
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION SUBSPH(*), SCRRED(*), SCRRED2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUIN2LIST(*), LUIN3LIST(*)
      DIMENSION LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK, ISTRED, ISTRED2
      CHARACTER*12 SECTID, WALLTID
C     initialize scratch vector offsets
      ISTRED  = 0
      ISTRED2 = 0
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH  = 0
      IOFFSET_IN_LUIN  = 0
      IOFFSET_IN_LUIN2 = 0
      IOFFSET_IN_LUIN3 = 0
      IOFFSET_LUOUT    = 0
      IOFFSET_INT_IN   = 0
      IOFFSET_INT_IN2  = 0
      IOFFSET_INT_IN3  = 0
      IOFFSET_INT_OUT  = 0
C     set scratch vectors to zero
      CALL DZERO(SCRRED,NROOT)
C
      IF( ( JOFF - 1 ) .gt. 0 ) CALL DZERO(SCRRED2, ( JOFF - 1) )
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector JOFF - 1
C       and JOFF2 - 1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( JOFF  - 1 )   * MY_VEC1_IOFF   +
     &                   ( JOFF2 - 1 )   * MY_VEC1_IOFF
        IOFFSET_INT_IN   = 1 + NUM_BLK                    +
     &                   ( JOFF  - 1 )   * MY_ACT_BLK1    +
     &                   ( JOFF2 - 1 )   * MY_ACT_BLK1
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK  &                  IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
        DO 100 IROOT = 1, NROOT
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IROOT - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN2  = 1 + NUM_BLK                     +
     &                     ( IROOT - 1 )   * MY_ACT_BLK1
C
C         read in batch ISBATCH from LUIN2 to VEC1
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) '1st VEC1 read from position',IROOT -1
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
          IF( IROOT .eq. 1 ) ISTRED = IROOT + ( JOFF - 1 ) * 2 * NROOT
C
C         SCRRED(IROOT) == SCRRED(IROOT) - SUM ( VEC1 * VEC2 )
          SCRRED(IROOT) =  SCRRED(IROOT)
     &                  -  DDOT( LEBATCH(ISBATCH), VEC1, 1, VEC2, 1 )
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 KVEC = 1, (JOFF - 1)
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( KVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN3  = 1 + NUM_BLK                     +
     &                     ( KVEC - 1 )    * MY_ACT_BLK1
C
          CALL RDVEC_BATCH_DRV4(LUIN3,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN3,IOFFSET_INT_IN3,
     &                          LUIN3LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) '2nd VEC1 read from position',KVEC -1
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
          IF( KVEC .eq. 1 ) ISTRED2 = NROOT + KVEC + (JOFF - 1 )*2*NROOT
C         SCRRED2(KVEC) == SCRRED2(KVEC) - SUM ( VEC1 * VEC2 )
          SCRRED2(KVEC) =  SCRRED2(KVEC) 
     &                  -  DDOT( LEBATCH(ISBATCH), VEC1, 1, VEC2, 1 )
C
  200   CONTINUE
C
C       write vector VEC2 from batch IBATCH to file LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )     * MY_VEC1_IOFF
        IOFFSET_INT_OUT    = 1 + NUM_BLK    +
     &                     ( JOFF - 1 )     * MY_ACT_BLK1
C
CSK     WRITE(LUWRT,*) 'final VEC2 to write on position',JOFF -1
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT,IOFFSET_INT_OUT,
     &                       LUOUTLIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches
C
C     individual dot products need to be communicated
C
      starttime = MPI_WTIME()
C
      CALL DZERO(SUBSPH(ISTRED),NROOT)
      CALL REDVEC(SCRRED,SUBSPH(ISTRED),NROOT,2,
     &            MPI_SUM,MPI_COMM_WORLD,-1)
C
      IF( (JOFF - 1) .gt. 0 ) THEN
        CALL DZERO(SUBSPH(ISTRED2),JOFF-1)
        CALL REDVEC(SCRRED2,SUBSPH(ISTRED2),JOFF-1,2,
     &              MPI_SUM,MPI_COMM_WORLD,-1)
      END IF
C
      endtime = MPI_WTIME()
      WALLTID = SECTID(endtime-starttime)
csk   WRITE(LUWRT,*) 'THIS IS WHAT I HAVE REDUCED TO'
csk   WRITE(LUWRT,*) ' ISTRED, ISTRED2, NROOT, JOFF',
csk  &                  ISTRED, ISTRED2, NROOT, JOFF-1
csk   CALL WRTMATMN(SUBSPH(ISTRED),1,NROOT,1,NROOT,LUWRT)
csk   CALL WRTMATMN(SUBSPH(ISTRED2),1,JOFF-1,1,JOFF-1,LUWRT)
      IF( TIMING_par )
     &WRITE(LUWRT,9450) WALLTID
 9450 FORMAT(' >>>  WALL TIME FOR REDUCING WORK                : ',A)
C
      END
***********************************************************************

      SUBROUTINE ORTHG_VEC_BATCH_LUCI2(VEC1,VEC2,FACIN1,FACIN2,NINVEC1,
     &                                 NINVEC2,NBATCH,LBATCH,LEBATCH,
     &                                 I1BATCH,IBATCH,MY_IOFF_LUIN,
     &                                 MY_IOFF_LUIN2,MY_IOFF_LUIN3,
     &                                 MY_IOFF_LUOUT,MY_IOFF_LUOUT2,
     &                                 SCAL_ARRAY,NSCALING,IADD,
     &                                 LUINLIST,LUIN2LIST,LUIN3LIST,
     &                                 LUOUTLIST,LUOUT2LIST,
     &                                 LUIN,LUIN2,LUIN3,LUOUT,LUOUT2)
C
C     add vectors on file LUIN1 and LUIN2 resp. LUIN2 and LUIN3. 
C     store resulting vector VEC2 on LUOUT with scaling factor FACOFSCAL.
C     LUOUT(VEC2) = LUIN2(VEC2) x FACOFSCAL
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - March 2008
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION SCAL_ARRAY(*), FACIN1(*), FACIN2(1)
      DIMENSION LUINLIST(*),LUIN2LIST(*), LUIN3LIST(*)
      DIMENSION LUOUTLIST(*), LUOUT2LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT, MY_IOFF_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT, IOFFSET_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
      REAL*8 REDFAC
C     initialize scratch factor 
      REDFAC = 0.0D0
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH  = 0
      IOFFSET_IN_LUIN  = 0
      IOFFSET_IN_LUIN2 = 0
      IOFFSET_IN_LUIN3 = 0
      IOFFSET_LUOUT    = 0
      IOFFSET_LUOUT2   = 0
      IOFFSET_INT_IN   = 0
      IOFFSET_INT_IN2  = 0
      IOFFSET_INT_IN3  = 0
      IOFFSET_INT_OUT  = 0
      IOFFSET_INT_OUT2 = 0
C
csk   WRITE(LUWRT,*) 'THIS IS WORK_SP(1+(JVEC-1)*2*NROOT) inside',myproc
csk   CALL WRTMATMN(FACIN1,1,NINVEC1,1,NINVEC1,LUWRT)
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       find offset 
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                   ( NINVEC2 - 1 )    * MY_VEC1_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK                     +
     &                   ( NINVEC2 - 1 )    * MY_ACT_BLK1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) ' VEC2 read from position',NINVEC2 - 1
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       calculate scalar products with vectors on LUIN
C   
        DO 100 IVEC = 1, NINVEC1
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN  + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK                     +
     &                     ( IVEC - 1 )    * MY_ACT_BLK1
C
          CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                          LUINLIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) ' VEC1 read from position',IVEC - 1
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          IF( IVEC .eq. 1 ) THEN
            CALL DAXPY(LEBATCH(ISBATCH),FACIN1(1),VEC1,1,VEC2,1)
          ELSE
            CALL DAXPY(LEBATCH(ISBATCH),FACIN1(IVEC),VEC1,1,VEC2,1)
          END IF
C
  100   CONTINUE
C
CSK     WRITE(LUWRT,*) ' VEC2 after 1st DAXPY part '
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       continue with vectors on LUIN3
C
        DO 200 IVEC = 1, ( NINVEC2 - 1)
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN3  = 1 + NUM_BLK                     +
     &                     ( IVEC - 1 )    * MY_ACT_BLK1
C
          CALL RDVEC_BATCH_DRV4(LUIN3,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN3,IOFFSET_INT_IN3,
     &                          LUIN3LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'VEC1 read from LUIN3 at IVEC - 1', IVEC - 1
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
C    &                  LUWRT)
C
C         VEC2 == VEC2 + VEC1 * FACIN2(IVEC)
C
          IF( IVEC .eq. 1 ) THEN
            CALL DAXPY(LEBATCH(ISBATCH),FACIN2(1),VEC1,1,VEC2,1)
          ELSE
            CALL DAXPY(LEBATCH(ISBATCH),FACIN2(IVEC),VEC1,1,VEC2,1)
          END IF
C
  200   CONTINUE
C
CSK     WRITE(LUWRT,*) 'VEC2 after 2nd DAXPY == after sec ort'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
C    &                LUWRT)
C
C       update REDFAC == REDFAC + SUM ( VEC2 * VEC2 )
        REDFAC = REDFAC + DDOT( LEBATCH(ISBATCH), VEC2, 1, VEC2, 1 )
CSK     WRITE(LUWRT,*) 'REDFAC = REDFAC + DDOT ...',REDFAC
C
C       write vector VEC2 from batch IBATCH to file LUIN2
C
        IOFFSET_IN_LUIN2 =  MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                   ( NINVEC2 - 1 )  * MY_VEC1_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK    +
     &                   ( NINVEC2 - 1 )  * MY_ACT_BLK1
C
CSK     WRITE(LUWRT,*) 'final VEC2 to write on position',NINVEC2 -1
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
        CALL WTVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                       LUIN2LIST,NUM_ACTIVE_BATCH)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches - part 1 -
C
C     communicate REDFAC to obtain correct scaling factor
      
      FACTOR = 0.0D0
      CALL REDVEC(REDFAC,FACTOR,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
csk   WRITE(LUWRT,*) 'REDUCED FACTOR IN PART 3.2',FACTOR
C
      SCALING = 0.0D0
      SCALING = 1.0D0/SQRT(FACTOR)
csk   WRITE(LUWRT,*) 'SCALING FACTOR IN PART 3.2',SCALING
C
C     scaling part 1
C
      CALL DSCAL(NSCALING, SCALING, SCAL_ARRAY, 1 )
C
csk   WRITE(LUWRT,*) 'THIS IS MY SCALED WORK PART'
csk   CALL WRTMATMN(SCAL_ARRAY,1,NSCALING,1,NSCALING,LUWRT)
C
C     initialize scratch offsets again
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
C
C     we loop again
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO( VEC2, LEBATCH(ISBATCH) )
C
C       read vector from file LUIN2
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                   ( NINVEC2 - 1 )    * MY_VEC1_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK                     +
     &                   ( NINVEC2 - 1 )    * MY_ACT_BLK1
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       scaling - part 2 -
C
        CALL DSCAL(LEBATCH(ISBATCH),SCALING,VEC2,1)
C
C       write scaled vector to files LUOUT and LUOUT2
C       writing to LUOUT1 only for ( NINVEC2 -1 + 1 .ne. IADD )
C
        IF( ( NINVEC2 - 1 + 1 ) .ne. IADD ) THEN
C
C         new offset for writing on LUOUT
C
          IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                       ( NINVEC2 - 1 )  * MY_VEC1_IOFF
C
          IOFFSET_INT_OUT    = 1 + NUM_BLK    +
     &                       ( NINVEC2 - 1 )  * MY_ACT_BLK1
C
          CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT,IOFFSET_INT_OUT,
     &                          LUOUTLIST,NUM_ACTIVE_BATCH)
C
        END IF
C
C       new offset for writing on LUOUT2
C
        IOFFSET_LUOUT2     =  MY_IOFF_LUOUT2 + IOFFSET_SCRATCH +
     &                     ( NINVEC2 - 1 )   * MY_VEC1_IOFF
C
        IOFFSET_INT_OUT2   = 1 + NUM_BLK     +
     &                     ( NINVEC2 - 1 )   * MY_ACT_BLK1
C
        CALL WTVEC_BATCH_DRV4(LUOUT2,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT2,IOFFSET_INT_OUT2,
     &                        LUOUT2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'This is what i have written'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches 2
      END
***********************************************************************

      SUBROUTINE ORTHG_VEC_BATCH_LUCI3(VEC1,VEC2,FACIN1,FACIN3,FACIN2,
     &                                 NROOT,JOFF,NBATCH,LBATCH,LEBATCH,
     &                                 I1BATCH,IBATCH,MY_IOFF_LUIN,
     &                                 MY_IOFF_LUIN2,MY_IOFF_LUIN3,
     &                                 MY_IOFF_LUOUT,LUINLIST,LUIN2LIST,
     &                                 LUIN3LIST,LUOUTLIST,JOFF2,
     &                                 LUIN,LUIN2,LUIN3,LUOUT)
C
C     calculate vector products between vector on file LUIN1 and LUIN2 
C     resp. LUIN1 and LUIN3. Store vector VEC2 on LUOUT 
C                            ==> LUIN1(VEC2) --> LUOUT(VEC2)
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: JOFF  = JOFF
C           JOFF2 = JOFF2
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION FACIN1(*), FACIN3(1)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUIN2LIST(*), LUIN3LIST(*)
      DIMENSION LUOUTLIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_IN_LUIN3  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_IN3   = 0
      IOFFSET_INT_OUT   = 0
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       position in file is at the end of vector JOFF - 1
C       and JOFF2 - 1
C
        IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                   ( JOFF  - 1 )   * MY_VEC1_IOFF   +
     &                   ( JOFF2 - 1 )   * MY_VEC1_IOFF
        IOFFSET_INT_IN   = 1 + NUM_BLK                    +
     &                   ( JOFF  - 1 )   * MY_ACT_BLK1    +
     &                   ( JOFF2 - 1 )   * MY_ACT_BLK1
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK  &                  IOFFSET_IN_LUIN
C
        CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                        LUINLIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
        DO 100 IROOT = 1, NROOT
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN2  = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                      ( IROOT - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN2   = 1 + NUM_BLK                     +
     &                      ( IROOT - 1 )   * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK  &                    IOFFSET_IN_LUIN2
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
C         scaling of VEC2
C
          IF( IROOT .eq. 1 ) THEN 
            CALL DSCAL(LEBATCH(ISBATCH),FACIN2,VEC2,1)
          END IF
C
C         VEC2 == VEC2 + VEC1 * FACIN1(IROOT)
C
          CALL DAXPY(LEBATCH(ISBATCH),FACIN1(IROOT),VEC1,1,VEC2,1)
CSK
CSK       WRITE(LUWRT,*) 'VEC1 and VEC2 after DAXPY == after first ort'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH))
CSK  &                  LUWRT)
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
  100   CONTINUE
C
C       continue with vectors on LUIN3
C
        DO 200 KVEC = 1, (JOFF - 1)
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3  = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                      ( KVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN3   = 1 + NUM_BLK                     +
     &                      ( KVEC - 1 )    * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN3',
CSK  &                    IOFFSET_IN_LUIN3
C
          CALL RDVEC_BATCH_DRV4(LUIN3,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN3,IOFFSET_INT_IN3,
     &                          LUIN3LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN3'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
C         VEC2 == VEC2 + VEC1 * FACIN3(KVEC)
C
          IF( KVEC .eq. 1 ) THEN
            CALL DAXPY(LEBATCH(ISBATCH),FACIN3(1),VEC1,1,VEC2,1)
          ELSE
            CALL DAXPY(LEBATCH(ISBATCH),FACIN3(KVEC),VEC1,1,VEC2,1)
          END IF
C
  200   CONTINUE
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT     =  MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                    ( JOFF - 1 )   * MY_VEC1_IOFF
C
        IOFFSET_INT_OUT   = 1 + NUM_BLK                      +
     &                    ( JOFF - 1 )   * MY_ACT_BLK1
C
        CALL WTVEC_BATCH_DRV4(LUOUT,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT,IOFFSET_INT_OUT,
     &                        LUOUTLIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'This is what I have written - sigma part'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
      END DO
C     ^ loop over batches
C
      END
***********************************************************************

      SUBROUTINE CALC_SUBSPACE_H_LUCI(VEC1,VEC2,SCRRED,IMUSTRED,IVEC,
     &                                NROOT,ISTRED,NBATCH,LBATCH,
     &                                LEBATCH,I1BATCH,IBATCH,ILOOP,
     &                                MY_IOFF_LUIN,MY_IOFF_LUIN2,
     &                                MY_IOFF_LUIN3,MY_IOFF_LUIN4,
     &                                LUINLIST,LUIN2LIST,LUIN3LIST,
     &                                LUIN4LIST,LUIN,LUIN2,LUIN3,LUIN4)
C
C     calculate dot products between vectors on file LUIN1 / LUIN2 
C     and LUIN3 resp. LUIN4
C
C     active blocks on the MPI-file are flagged by a nonzero length
C
C     Note: IVEC  = IVEC
C           ILOOP = ILOOP
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*),VEC2(*)
      DIMENSION SCRRED(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUIN2LIST(*), LUIN3LIST(*)
      DIMENSION LUIN4LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN, MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN4
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3, IOFFSET_IN_LUIN4
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH  = 0
      IOFFSET_IN_LUIN  = 0
      IOFFSET_IN_LUIN2 = 0
      IOFFSET_IN_LUIN3 = 0
      IOFFSET_IN_LUIN4 = 0
      IOFFSET_INT_IN   = 0
      IOFFSET_INT_IN2  = 0
      IOFFSET_INT_IN3  = 0
      IOFFSET_INT_IN4  = 0
C 
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        IF( ILOOP .le. 0 )THEN
C
          IOFFSET_IN_LUIN   = MY_IOFF_LUIN  + IOFFSET_SCRATCH +
     &                      ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN    = 1 + NUM_BLK                     +
     &                      ( IVEC - 1 )    * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN ',
CSK  &                    IOFFSET_IN_LUIN 
C
          CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                          LUINLIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC2 on LUIN '
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
        ELSE 
C
          IOFFSET_IN_LUIN2  = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                      ( ILOOP - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN2   = 1 + NUM_BLK                     +
     &                      ( ILOOP - 1 )   * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK  &                    IOFFSET_IN_LUIN2
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC2 on LUIN2'
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
        END IF
C
        DO 100 JVEC = 1, MIN(IVEC,NROOT)
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN3  = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                      ( JVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN3   = 1 + NUM_BLK                     +
     &                      ( JVEC - 1 )    * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN3',
CSK  &                    IOFFSET_IN_LUIN3
C
          CALL RDVEC_BATCH_DRV4(LUIN3,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN3,IOFFSET_INT_IN3,
     &                          LUIN3LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN3'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          IJ = IVEC * ( IVEC - 1 ) / 2 + JVEC
C         SCRRED(IJ) = SCRRED(IJ) + SUM ( VEC1 * VEC2 )
          SCRRED(IJ) = SCRRED(IJ) + 
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
          IF (ISBATCH .eq. 1 ) THEN
            IMUSTRED = IMUSTRED + 1
            IF( IVEC .eq. 1 .and. JVEC .eq. 1 ) ISTRED = IJ
          END IF
C
  100   CONTINUE
C
C       continue with vectors on LUIN4
C
        ILOOP2 = 0
C
        DO 200 JVEC = ( NROOT + 1), IVEC
C
          ILOOP2 = ILOOP2 + 1
C
C         set new offset and zero read-in vector
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          IOFFSET_IN_LUIN4  = MY_IOFF_LUIN4 + IOFFSET_SCRATCH +
     &                      ( ILOOP2 - 1 )  * MY_VEC1_IOFF
          IOFFSET_INT_IN4   = 1 + NUM_BLK                     +
     &                      ( ILOOP2 - 1 )  * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN4',
CSK  &                    IOFFSET_IN_LUIN4
C
          CALL RDVEC_BATCH_DRV4(LUIN4,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN4,IOFFSET_INT_IN4,
     &                          LUIN4LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN4'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          IJ = IVEC * ( IVEC - 1 ) / 2 + JVEC
C         SCRRED(IJ) = SCRRED(IJ) + SUM ( VEC1 * VEC2 )
          SCRRED(IJ) = SCRRED(IJ) +
     &                 DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
          IF (ISBATCH .eq. 1 ) IMUSTRED = IMUSTRED + 1
C
  200   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C     ^ loop over batches
C
      END
***********************************************************************

      SUBROUTINE PART_CIV_PAR1(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,IBLOCKD,NDIM)
C
C     Partition a CI vector into batches of blocks. The length of a
C     batch must be atmost MXLNG.
C     S-vector routine.
C
C     IF ICOMP. eq. 1: the complete CI vector is constructed in just one
C     batch.
C
C     OUTPUT
C     ======
C
C     NBATCH      : Number of batches
C     LBATCH(*)   : Number of blocks in a given batch
C     LEBATCH(*)  : Number of elements in a given batch ( packed ) !
C     I1BATCH(*)  : Number of first block in a given batch
C     IBATCH(8,*) : TTS blocks in Start of a given TTS block with respect to
C                   start
C     of batch --
C     IBATCH(1,*) : Alpha type
C     IBATCH(2,*) : Beta  type
C     IBATCH(3,*) : Sym of alpha
C     IBATCH(4,*) : Sym of beta
C     IBATCH(5,*) : Offset of block with respect to start of block in
C                   expanded form
C     IBATCH(6,*) : Offset of block with respect to start of block in
C                   packed form
C     IBATCH(7,*) : Length of block, expandend form
C     IBATCH(8,*) : Length of block, packed form
C    
C     original version : Jeppe Olsen     - August 1995
C     parallel adaption: S. Knecht       - March  2007 
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
C     Input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION IBLOCKD(NDIM)
C     Output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
*
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) '    PART_CIV_PAR1   '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(6,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(6,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(6,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
*
*. block zero
*
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
*. Loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
*. Loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
*. Should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
*. can this block be included
      IBSM = ISMOST(ISM)
      NSTA = NSSOA(ISM,IA)
      NSTB = NSSOB(IBSM,IB)
      LBLOCK= NSTA*NSTB
C     set unpacked length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. LUCI_MYPROC ) LBLOCK = 0
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
C     set packed length for unused block to zero
      IF( IBLOCKD(IBLOCK+1) .ne. LUCI_MYPROC ) LBLOCKP = 0
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
C       only s-blocks corresponding to distribution
        IF( IBLOCKD(IBLOCK) .ne. LUCI_MYPROC ) IBATCH(1,IBLOCK) = 0
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
        IBATCH(7,IBLOCK) = LBLOCK
        IBATCH(8,IBLOCK) = LBLOCKP
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single Block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV_PAR1'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise Buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV_PAR1. ' )
      ELSE
*. This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
*
      IF(NTEST.NE.0) THEN
        WRITE(luwrt,*) 'Output from PART_CIV_PAR1'
        WRITE(luwrt,*) '========================='
        WRITE(luwrt,*)
        WRITE(luwrt,*) ' Number of batches ', NBATCH
        DO JBATCH = 1, NBATCH
          WRITE(luwrt,*)
          WRITE(luwrt,*) ' Info on batch ', JBATCH
          WRITE(luwrt,*) ' *********************** '
          WRITE(luwrt,*)
          WRITE(luwrt,*) ' Number of blocks included ', LBATCH(JBATCH)
          WRITE(luwrt,*) ' TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(luwrt,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
      END IF
*
      END
***********************************************************************

      SUBROUTINE PART_CIV_PAR2(IDC,IBLTP,NSSOA,NSSOB,NOCTPA,NOCTPB,
     &                         NSMST,MXLNG,IOCOC,ISMOST,
     &                         NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                         ICOMP,ITTSS_ORD,SCLFAC,NDIM,IEXCLBLK)
C
C     Partition a CI vector into batches of blocks. The length of a
C     batch must be atmost MXLNG.
C     C-vector routine.
C
C     IF ICOMP. eq. 1: the complete CI vector is constructed in just one
C     batch.
C
C     OUTPUT
C     ======
C
C     NBATCH      : Number of batches
C     LBATCH(*)   : Number of blocks in a given batch
C     LEBATCH(*)  : Number of elements in a given batch ( packed ) !
C     I1BATCH(*)  : Number of first block in a given batch
C     IBATCH(8,*) : TTS blocks in Start of a given TTS block with respect to
C                   start
C     of batch --
C     IBATCH(1,*) : Alpha type
C     IBATCH(2,*) : Beta  type
C     IBATCH(3,*) : Sym of alpha
C     IBATCH(4,*) : Sym of beta
C     IBATCH(5,*) : Offset of block with respect to start of block in
C                   expanded form
C     IBATCH(6,*) : Offset of block with respect to start of block in
C                   packed form
C     IBATCH(7,*) : Length of block, expandend form
C     IBATCH(8,*) : Length of block, packed form
C    
C     original version : Jeppe Olsen     - August 1995
C     parallel adaption: S. Knecht       - March  2007 
C
C     Last revision:     S. Knecht       - March  2007
C
************************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
C     input
      INTEGER NSSOA(NSMST,*),NSSOB(NSMST,*)
      INTEGER IOCOC(NOCTPA,NOCTPB)
      INTEGER IBLTP(*)
      INTEGER ISMOST(*)
      DIMENSION SCLFAC(NDIM)
C     output
      INTEGER LBATCH(*)
      INTEGER LEBATCH(*)
      INTEGER I1BATCH(*)
      INTEGER IBATCH(8,*)
C     scratch
      INTEGER LBLOCK_SAVE, LBLOCKP_SAVE 
C
      LBLOCK_SAVE  = 0
      LBLOCKP_SAVE = 0 
C
      NTEST = 0000
      IF(NTEST.GE.100) THEN
        WRITE(6,*)
        WRITE(6,*) ' ==================='
        WRITE(6,*) '    PART_CIV_PAR2   '
        WRITE(6,*) ' ==================='
        WRITE(6,*) ' IDC = ', IDC
        WRITE(6,*)
        WRITE(6,*) ' IOCOC Array '
        CALL IWRTMA(IOCOC,NOCTPA,NOCTPB,NOCTPA,NOCTPB)
        if (NTEST.ge.500) then
          WRITE(6,*) ' NSSOA array ( input ) '
          CALL IWRTMA(NSSOA,NSMST,NOCTPA,NSMST,NOCTPA)
          WRITE(6,*) ' NSSOB array ( input ) '
          CALL IWRTMA(NSSOB,NSMST,NOCTPB,NSMST,NOCTPB)
          write(6,*) ' IBLTP array: '
          call iwrtma(IBLTP,1,NSMST,1,NSMST)
        end if
      END IF
C
C     block zero
C
      ILOOPBLK = 0
      IB = 1
      IA = 1
      ISM = 1
      IFRST = 1
      NBATCH = 0
      IBLOCK = 0
      IFINI = 0
C     loop over batches of blocks
 2000 CONTINUE
      NBATCH = NBATCH + 1
      LBATCH(NBATCH) = 0
      I1BATCH(NBATCH) = IBLOCK  + 1
      LENGTH = 0
      LENGTHP= 0
      NBLOCK = 0
      IFRST = 1
      IF( NBATCH .gt. 1 ) ILOOPBLK = ILOOPBLK - 1 
      
C     loop over blocks in batch
 1000 CONTINUE
      IF(IFRST.EQ.0) THEN
        call nxt_tts(ITTSS_ORD,IA,IB,ISM,IFINI,NOCTPA,NOCTPB,NSMST)
      END IF
      IFRST = 0
      IF (IFINI.EQ.1) GOTO 2002
C     should this block be included
      IF(IBLTP(ISM).EQ.0) GOTO 1000
      IF(IBLTP(ISM).EQ.2.AND.IA.LT.IB) GOTO 1000
      IF(IOCOC(IA,IB).EQ.0) GOTO 1000
      ILOOPBLK = ILOOPBLK + 1
C     can this block be included
      IBSM   = ISMOST(ISM)
      NSTA   = NSSOA(ISM,IA)
      NSTB   = NSSOB(IBSM,IB)
      LBLOCK = NSTA*NSTB
      LBLOCK_SAVE = LBLOCK
      IF(IBLTP(ISM).EQ.1.OR.(IBLTP(ISM).EQ.2.AND.IA.NE.IB)) THEN
        LBLOCKP = NSTA*NSTB
      ELSE IF (IBLTP(ISM) .EQ. 2.AND.IA.EQ.IB) THEN
        LBLOCKP = NSTA*(NSTA+1)/2
      END IF
      LBLOCKP_SAVE = LBLOCKP
C     check if we need this block, if not: 
C     LBLOCKP = 0
C     LBLOCK  = 0
CSK      IF(IEXCLBLK .eq. 1) THEN
CSK       WRITE(6,*) 'ILOOPBLK,SCLFAC(ILOOPBLK),IBLOCK',ILOOPBLK,
CSK     &             SCLFAC(ILOOPBLK),IBLOCK
CSK      END IF
      IF( IEXCLBLK .eq. 1 ) THEN
         IF( SCLFAC(ILOOPBLK) .eq. 0.0D0 ) THEN
           LBLOCKP = 0
           LBLOCK  = 0
           GOTO 1000
         END IF
      ELSE
        IF( SCLFAC(IBLOCK+1) .eq. 0.0D0 ) THEN
          LBLOCKP = 0
          LBLOCK  = 0
        END IF
      END IF
C?    write(6,*) ' IA IB ISM LBLOCK ', IA,IB,ISM,LBLOCK,myproc
      IF(LENGTH+LBLOCK.LE.MXLNG.OR.ICOMP.EQ.1) THEN
        NBLOCK = NBLOCK + 1
        IBLOCK = IBLOCK + 1
        LBATCH(NBATCH) = LBATCH(NBATCH)+1
        IBATCH(1,IBLOCK) = IA
        IBATCH(2,IBLOCK) = IB
        IBATCH(3,IBLOCK) = ISM
        IBATCH(4,IBLOCK) = IBSM
        IBATCH(5,IBLOCK) = LENGTH+1
        IBATCH(6,IBLOCK) = LENGTHP+1
C       keep length information, needed?
        IBATCH(7,IBLOCK) = LBLOCK_SAVE
        IBATCH(8,IBLOCK) = LBLOCKP_SAVE
C       all blocks are included, but only active blocks have a length
        LENGTH = LENGTH + LBLOCK
        LENGTHP= LENGTHP+ LBLOCKP
        LEBATCH(NBATCH) = LENGTHP
        GOTO 1000
      ELSE IF(ICOMP.EQ.0.AND.
     &  LENGTH+LBLOCK.GT. MXLNG .AND. NBLOCK.EQ.0) THEN
        WRITE(6,*) ' Not enough scratch space to include a single block'
        WRITE(6,*) ' Since I cannot procede I will stop '
        WRITE(6,*) ' Insufficient buffer detected in PART_CIV_PAR2'
        write(6,*) '  LENGTH,LBLOCK ',LENGTH,LBLOCK
        WRITE(6,*) ' Alter GAS space of raise buffer from ', MXLNG
        call quit( ' Insufficient buffer space in PART_CIV_PAR2. ' )
      ELSE
C       This batch is finished, goto next batch
        GOTO 2000
      END IF
 2002 CONTINUE
C
      IF(NTEST.NE.0) THEN
       IF(IEXCLBLK .eq. 1 ) THEN
        WRITE(6,*) 'Output from PART_CIV_PAR2'
        WRITE(6,*) '========================='
        WRITE(6,*)
        WRITE(6,*) ' Number of batches (LUCI_MYPROC)',
     &   NBATCH, LUCI_MYPROC
        DO JBATCH = 1, NBATCH
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*)
          WRITE(6,*) '      Length of batch           ', LEBATCH(JBATCH)
          WRITE(6,*) '      Number of blocks included ', LBATCH(JBATCH)
          WRITE(6,*) '      TTSS and offsets and lenghts of each block '
          DO IBLOCK = I1BATCH(JBATCH),I1BATCH(JBATCH)+ LBATCH(JBATCH)-1
            WRITE(6,'(10X,4I3,4I8)') (IBATCH(II,IBLOCK),II=1,8)
          END DO
        END DO
       END IF
      END IF
C
      END
***********************************************************************

      SUBROUTINE PRINT_BATCH_INFO(NCBATCH,LBATC,LEBATC,
     &           I1BATC,IBATC,ICCTOS,ICWEIGHT,NBLOCK,IBTOTW,ICWEIGHTF,
     &           IABSOLUTE_WEIGHT)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION LBATC(*),LEBATC(*),I1BATC(*),IBATC(8,*)
      DIMENSION ICCTOS(NBLOCK,NBLOCK),ICWEIGHT(NBLOCK)
      DIMENSION IBTOTW(NBLOCK),ICWEIGHTF(NBLOCK)
      INTEGER*8 IABSOLUTE_WEIGHT
*
      NTEST = 99
      IF(LUCI_MYPROC.EQ.LUCI_MASTER .and. NTEST .ge. 100 )THEN
        WRITE(6,*) ' OUTPUT from PRINT_BATCH_INFO'
        WRITE(6,*) ' ============================'
        WRITE(6,*) ' '
        WRITE(6,*) ' Number of batches ', NCBATCH
      END IF
      DO JBATCH = 1, NCBATCH
        IF(LUCI_MYPROC.EQ.LUCI_MASTER .and. NTEST .ge. 100 ) THEN
          WRITE(6,*)
          WRITE(6,*) ' Info on batch ', JBATCH
          WRITE(6,*) ' *********************** '
          WRITE(6,*) ' '
          WRITE(6,*) ' Number of blocks included ', LBATC(JBATCH)
        END IF
CSK        WRITE(6,*) '   TTSS and offsets and lenghts of each block '
        DO IBLOCK =I1BATC(JBATCH),I1BATC(JBATCH)+ LBATC(JBATCH)-1
CSK          WRITE(6,'(10X,4I3,4I8)') (IBATC(II,IBLOCK),II=1,8)
CSK          IF(IBATC(8,IBLOCK).NE.0.AND.LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10) 
CSK     &       WRITE(6,'(2X,A,1X,I5,1X,A,1X,I8)') 'Block number',IBLOCK,
CSK     & 'has length',(IBATC(II,IBLOCK),II=8,8)
          JJJ = 0
          DO IJ = 1, NBLOCK
            IF(ICCTOS(IJ,IBLOCK).NE.0) THEN
              JJJ = JJJ + 1
              ICWEIGHT(JJJ) = IJ
            ENDIF
          END DO
          ICWEIGHTF(IBLOCK) = JJJ
          IF(JJJ.GT.0.AND.LUCI_MYPROC.EQ.LUCI_MASTER
     &      .AND.NTEST.GE.100) THEN
            WRITE(6,'(2X,A,I4)')'Number of connections for this block:',
     &                           JJJ
            CALL IWRTMA(ICWEIGHT,1,JJJ,1,JJJ)
            
            WRITE(6,*)'  '
          ENDIF
          IBTOTW(IBLOCK) = (JJJ) * IBATC(8,IBLOCK)
        END DO
      END DO
      IF(LUCI_MYPROC.EQ.LUCI_MASTER.AND.NTEST.GE.10) THEN
CSK        WRITE(6,'(2X,A)')'TOTAL WEIGHT FOR EACH BLOCK'
CSK        CALL IWRTMA8(IBTOTW,1,NBLOCK,1,NBLOCK)
CSK        WRITE(6,'(2X,A)')'WEIGHT MATRIX'
CSK        CALL IWRTMA(ICWEIGHTF,1,NBLOCK,1,NBLOCK)
      END IF
      DO IBLK = 1, NBLOCK
        IABSOLUTE_WEIGHT = IABSOLUTE_WEIGHT + 
     &                     ( ICWEIGHTF(IBLK) * IBATC(8,IBLK) )
      END DO
*
      END
***********************************************************************

      SUBROUTINE PRSYMMN(A,MATDIM,LUWRT)
C PRINT LOWER HALF OF A SYMMETRIC MATRIX OF DIMENSION MATDIM.
C THE LOWER HALF OF THE MATRIX IS SUPPOSED TO BE IN VECTOR A.
#include "implicit.h"
      DIMENSION A(*)
      JSTART=1
      JSTOP=0
      DO 100 I=1,MATDIM
        JSTART=JSTART+I-1
        JSTOP=JSTOP +I
        WRITE(LUWRT,1010) I,(A(J),J=JSTART,JSTOP)
  100 CONTINUE
      RETURN
 1010 FORMAT(/I6,5E14.7,/,(6X,1P,5E14.7) )
      END
***********************************************************************

      SUBROUTINE REDVEC(VECIN,VECOUT,NDIM,NTYPE,NOP,NCOMM,NTARGET)
*
*. Reduce an array VECIN of dimension NDIM that is distributed over 
*. NCOMM to VECOUT on either a single process (NTARGET >= 0) or on 
*. all processes (NTARGET < 0). NTYPE determines the MPI_DATATYPE 
*. used in the reduce-call. NOP is the reduce-operation.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VECIN(*),VECOUT(*)
      INTEGER DATATYPE
*
      IF(NTYPE .EQ. 1) THEN
        DATATYPE = my_MPI_INTEGER
      ELSE IF (NTYPE .EQ. 2) THEN
C MPI_REAL8 is equivalent to MPI_DOUBLE_PRECISION
C the old LAM-MPI library does not support MPI_ALLREDUCE for MPI_REAL8
C so we use instead MPI_DOUBLE_PRECISION
CSK        DATATYPE = MPI_REAL8
        DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (NTYPE .EQ. 3)THEN
        DATATYPE = MPI_REAL
      ELSE IF (NTYPE .EQ. 4 .OR. NTYPE .EQ. 5) THEN
        DATATYPE  =  MPI_CHARACTER 
      ELSE
        WRITE(6,*)' REDVEC: TYPE "',NTYPE,'" does not exist!'
        CALL quit('ERROR in REDVEC: Nonexisting DATATYPE')
      END IF
*
      IF(NTARGET.GE.0)THEN
        CALL MPI_REDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,NTARGET,
     &                  NCOMM,IERR)
      ELSE IF(NTARGET.LT.0)THEN
        CALL MPI_ALLREDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,
     &                     NCOMM,IERR)
      ENDIF
*. return since everything is reduced now!
*
      END
***********************************************************************

      SUBROUTINE REDVEC_REL(VECIN,VECOUT,NDIM,NTYPE,NOP,NCOMM,NTARGET)
*
*. Reduce an array VECIN of dimension NDIM that is distributed over 
*. NCOMM to VECOUT on either a single process (NTARGET >= 0) or on 
*. all processes (NTARGET < 0). NTYPE determines the MPI_DATATYPE 
*. used in the reduce-call. NOP is the reduce-operation.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
      DIMENSION VECIN(*),VECOUT(*)
      INTEGER DATATYPE
*
      IF(NTYPE .EQ. 1) THEN
        DATATYPE = my_MPI_INTEGER
      ELSE IF (NTYPE .EQ. 2) THEN
C MPI_REAL8 is equivalent to MPI_DOUBLE_PRECISION
C the old LAM-MPI library does not support MPI_ALLREDUCE for MPI_REAL8
C so we use instead MPI_DOUBLE_PRECISION
CSK        DATATYPE = MPI_REAL8
        DATATYPE = MPI_DOUBLE_PRECISION
      ELSE IF (NTYPE .EQ. 3)THEN
        DATATYPE = MPI_REAL
      ELSE IF (NTYPE .EQ. 4 .OR. NTYPE .EQ. 5) THEN
        DATATYPE  =  MPI_CHARACTER 
      ELSE
        WRITE(6,*)' REDVEC: TYPE "',NTYPE,'" does not exist!'
        CALL quit('ERROR in REDVEC: Nonexisting DATATYPE')
      END IF
*
      IF(NTARGET.GE.0)THEN
        CALL MPI_REDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,NTARGET,
     &                  NCOMM,IERR)
      ELSE IF(NTARGET.LT.0)THEN
        CALL MPI_ALLREDUCE(VECIN,VECOUT,NDIM,DATATYPE,NOP,
     &                     NCOMM,IERR)
      ENDIF
*. return since everything is reduced now!
*
      END
***********************************************************************

      SUBROUTINE UPDATE_GEN_LIST(ISLIST1,ISLIST2,NDIM)
C
C     OUTPUT
C     ======
C
C     updated content of ISLIST1
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
#include "implicit.h"
      DIMENSION ISLIST1(*), ISLIST2(*)
      INTEGER NZERO
C
      NZERO = 0
C
      DO IBLK = 1, NDIM
C
        IF( ISLIST1( IBLK ) .ne. NZERO ) THEN
          IF( ISLIST2( IBLK ) .eq. NZERO ) THEN
            ISLIST1( IBLK ) = NZERO
          END IF
        END IF
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE IFACTOSFAC(ISLIST1,SLIST2,NDIM)
C
C     OUTPUT
C     ======
C
C     updated content of SLIST2
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
#include "implicit.h"
      DIMENSION ISLIST1(*), SLIST2(*)
      INTEGER NZERO
C
      NZERO = 0
      ONE  = 1.0D0
C
      DO IBLK = 1, NDIM
C
        IF( ISLIST1( IBLK ) .ne. NZERO ) THEN
           SLIST2( IBLK ) = ONE
        END IF
C
      END DO
C
      END
***********************************************************************

      SUBROUTINE WRSVCD_PAR(LU,LBLK,VEC1,IPLAC,VAL,NSCAT,
     &                      IBLOCKL,IBLOCKD,NBLOCK,JPACK,IXROOT,
     &                      LU1LIST)
*
* Write scattered vector to disc, every node writes only it's part.
* Information about the length of a blcok is stored in IBLOCKL.
* Information about which block belongs to which node is stored in
* IBLOCKD. 
* IXROOT is the current root to place a 1.0D0.
* LU1LIST is the MPI file list.
*
      IMPLICIT REAL*8 (A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
*     input
      DIMENSION IPLAC(*),VAL(*),IBLOCKL(NBLOCK),IBLOCKD(NBLOCK)
*     scratch
      DIMENSION VEC1(*)
      LOGICAL FOUND_ONE
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR1, IOFF_SCR2
      DIMENSION LU1LIST(*)
      INTEGER IBOTTO
*
      FOUND_ONE = .FALSE.
C
C     initialize offset
      IOFF_SCR1 = 0
      IOFF_SCR2 = 0
      INT_IOFF1  = 0
      INT_IOFF2  = 0
      IBOTTO = 0
C
      IOFF_SCR1 = MY_LU1_OFF + MY_VEC2_IOFF * ( IXROOT - 1 )
      INT_IOFF1 = 1 + MY_ACT_BLK2 * ( IXROOT    - 1 )
csk   WRITE(LUWRT,*) 'OFFSET at START',IOFF_SCR1
csk   WRITE(LUWRT,*) 'INT_IOFF1 at START ',INT_IOFF1
csk   WRITE(LUWRT,*) 'IBLOCKD '
csk   CALL IWRTMAMN(IBLOCKD,1,NBLOCK,1,NBLOCK,LUWRT)
csk   WRITE(LUWRT,*) 'IBLOCKL '
csk   CALL IWRTMAMN(IBLOCKL,1,NBLOCK,1,NBLOCK,LUWRT)
csk   WRITE(LUWRT,*) 'IBLOCKL '
csk   CALL IWRTMAMN(IBLOCKL,1,NBLOCK,1,NBLOCK,LUWRT)
csk   WRITE(LUWRT,*) ' VAL '
csk   CALL WRTMATMN(VAL,1,10,1,10,LUWRT)
csk   WRITE(LUWRT,*) ' IPLAC '
csk   CALL IWRTMAMN(IPLAC,1,10,1,10,LUWRT)

C
      IONE = 1
      NTEST = 0
      IBOTTO = 1
      LBL = 0
csk   WRITE(LUWRT,*) ' IBOTTO',IBOTTO
*     loop over blocks
      DO II = 1, NBLOCK
        IF( IBLOCKD(II) .ne. LUCI_MYPROC ) GOTO 300
*       length
        LBL = IBLOCKL(II)
csk     WRITE(LUWRT,*) ' LBL and IBOTTO for block II',LBL,IBOTTO,II

        CALL DZERO(VEC1,LBL)
C       *******************************************************
C       *******************************************************
C              Find the places where a 1.0 should be set 
C       *******************************************************
C       *******************************************************
        DO 200 IEFF = 1, NSCAT
csk       WRITE(LUWRT,*) ' IPLAC(IEFF), IEFF',IPLAC(IEFF), IEFF
          IF( IPLAC(IEFF).GE.IBOTTO.AND.IPLAC(IEFF).LE.IBOTTO+LBL-1)THEN
             VEC1(IPLAC(IEFF)-IBOTTO+1) = VAL(IEFF)
             FOUND_ONE = .TRUE.
          END IF
C
          NTEST = 00
C
          IF(NTEST.GE.10) THEN
            IF( IPLAC(IEFF).GE.IBOTTO.AND.IPLAC(IEFF).LE.IBOTTO+LBL-1)
     &      write(6,*) 'IBOT,IBOTTO+LBL-1',IBOTTO,IBOTTO+LBL-1
            IF( IPLAC(IEFF).GE.IBOTTO.AND.IPLAC(IEFF).LE.IBOTTO+LBL-1)
     &      write(6,*) ' Catch : IPLAC(IEFF) VAL(IEFF) ',
     &      IPLAC(IEFF),VAL(IEFF)
          END IF
*
  200     CONTINUE
C
C         new offset
C
          IOFF_SCR1 = IOFF_SCR1 + IOFF_SCR2
          INT_IOFF1 = INT_IOFF1 + INT_IOFF2
C
          IF( FOUND_ONE ) THEN
csk         WRITE(LUWRT,*) 'THIS IS WHAT I WILL WRITE',LUCI_MYPROC
csk         CALL WRTMATMN(VEC1,1,LBL,1,LBL,LUWRT)
csk         WRITE(LUWRT,*) 'OFFSET',IOFF_SCR1
csk         WRITE(LUWRT,*) 'INT_IOFF1 is ',INT_IOFF1
            LU1LIST( INT_IOFF1 ) = IONE
            CALL MPI_FILE_WRITE_AT(ILU1,IOFF_SCR1,VEC1,LBL,
     &                             MPI_REAL8,my_STATUS,IERR)
          END IF
C
          FOUND_ONE = .FALSE.
C         keep track of correct offset
          IOFF_SCR2 = LBL
          INT_IOFF2 = IONE
  300     CONTINUE
C         keep IBOT up-to-date for all nodes
          LBL = IBLOCKL( II )
          IBOTTO = IBOTTO + LBL
      END DO  
*
      END
***********************************************************************
*                                                                     *
* LUCIAREL, written by Timo Fleig and Jeppe Olsen                     *
*           parallelization by Stefan Knecht                          *
*                                                                     *
***********************************************************************
      SUBROUTINE INPROD_B_PAR_RL(LUIN1,LUIN2,VEC1,VEC2,SUBSPH,NBATCH,
     &                           LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                           MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                           LUIN1LIST,LUIN2LIST,JOFF)
C
C     Written by  S. Knecht         - May 29 2007
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: JOFF = IVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - May  2007
C
************************************************************************
      IMPLICIT REAL*8 ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK = 0
      IOFFSET_SCRATCH = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_INT_IN1  = 0
      IOFFSET_INT_IN2  = 0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C       position in file is at the end of vector JOFF - 1
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                     ( JOFF - 1 )  * MY_VEC1_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JOFF - 1 ) * MY_ACT_BLK1
C
csk     WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
csk  &                  IOFFSET_IN_LUIN1
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                       LUIN1LIST,NUM_ACTIVE_BATCH)
C
csk     WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
csk     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
        DO 100 IVEC = 1, JOFF
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( IVEC - 1 )  * MY_VEC1_IOFF
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( IVEC - 1 ) * MY_ACT_BLK1
C
csk     WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
csk  &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
csk     WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
csk     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = JOFF*(JOFF-1)/2 + IVEC
C          SUBSPH(IJ) == VEC1 * VEC2
C
           IJ = JOFF*(JOFF-1)/2 + IVEC
C
           SUBSPH(IJ) = SUBSPH(IJ) + 
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
#else
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE PAR_DUMMY2
      END 
***********************************************************************
#endif /* defined VAR_MPI */
