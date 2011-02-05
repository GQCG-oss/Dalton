!dirac_copyright_start
!      Copyright (c) 2010 by the authors of DIRAC.
!      All Rights Reserved.
!
!      This source code is part of the DIRAC program package.
!      It is provided under a written license and may be used,
!      copied, transmitted, or stored only in accordance to the
!      conditions of that written license.
!
!      In particular, no part of the source code or compiled modules may
!      be distributed outside the research group of the license holder.
!      This means also that persons (e.g. post-docs) leaving the research
!      group of the license holder may not take any part of Dirac,
!      including modified files, with him/her, unless that person has
!      obtained his/her own license.
!
!      For information on how to get a license, as well as the
!      author list and the complete list of contributors to the
!      DIRAC program, see: http://dirac.chem.vu.nl
!dirac_copyright_end

***********************************************************************
*
*========================================================*
*   Common LUCITA and LUCIAREL I/O routines              *
*                                                        *
*   Originals by Jeppe Olsen                             *
*   collected by Timo Fleig, Jul 12 - Sep, 2002          *
*========================================================*
*
***********************************************************************
      SUBROUTINE COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILE
C   LBLK .gt. 0 : LBL (length of next record) is equal to LBLK
C        .le. 0 : LBL (length of next record) is read from LUIN
*
* Structure of output file (LUOUT) is inherited from input file (LUIN):
* if input file is packed, so is output file
*
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        REWIND LUIN
        REWIND LUOUT
      END IF
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC_LUCI(SEGMNT,LBL,KBLK,LUOUT)
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            CALL ZERORC(LBL,LUOUT,IAMPACK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
      SUBROUTINE COPVCDP(LUIN,LUOUT,SEGMNT,IREW,LBLK)
C
C COPY VECTOR ON FILE LUIN TO FILE LUOUT
*
* Packed version 
C
C
C LBLK DEFINES STRUCTURE OF FILE
C Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        REWIND LUIN
        REWIND LUOUT
      END IF
 
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC_LUCI(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
      SUBROUTINE DMTVCD(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,LBLK)
C
C  IF( INV .NE. 0 ) THEN
C    V3(I) = (V1(I)+FAC)-1 * V2(I)
C    LU3      LU1            LU2
C  IF( INV .EQ. 0 ) THEN
C    V3(I) = (V1(I)+FAC) * V2(I)
C    LU3         LU1        LU2
C WHERE V1 AND V2 ARE VECTORS ON FILES LU1 AND LU2,
C AND LU3 IS WRITTEN ON FILE LU3
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
      DIMENSION  VEC1(*),VEC2(*)
C
      IF ( IREW .NE. 0 ) THEN
          REWIND LU1
          REWIND LU2
          REWIND LU3
      END IF
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(lupri,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVCD : '
     &                     , LBL1,LBL2
          Call Abend2( ' DIFFERENT BLOCKSIZES IN DMTVCD ' )
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC_LUCI(VEC1,LBL1,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC_LUCI(VEC2,LBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
             CALL DIAVC2(VEC2,VEC2,VEC1,FAC,LBL1)
            ELSE
             CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
             CALL VECSUM(VEC2,VEC1,VEC2,1.0D0,FAC,LBL1)
           END IF
C          CALL TODSC_LUCI(VEC2,LBL1,KBLK,LU3)
         END IF
         CALL TODSC_LUCI(VEC2,LBL1,KBLK,LU3)
      END IF
C
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
***********************************************************************
      SUBROUTINE DMTVDS(VEC1,VEC2,LU1,LU2,LU3,FAC,IREW,INV,
     &                  ISCAT,XSCAT,NSCAT,LBLK,XINOUT)
*
* Multiply/divide elements of two vectors residing on disk
* Elements corresponding to absolute adresses in ISCAT
* are set to the elements of XSCAT
*
* For elements not in ISCAT the operation is thus :
*
* For INV.NE. 0 :
*
*    V3(I) = (V1(I)+FAC)-1 * V2(I)
*    LU3      LU1            LU2
*
* For INV.EQ. 0 :
*
*    V3(I) = (V1(I)+FAC) * V2(I)
*    LU3         LU1        LU2
*
* LBLK defines structure of files
*
* Overlap between input and output vector is also calculated
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
      DIMENSION  VEC1(*),VEC2(*)
      DIMENSION  XSCAT(NSCAT),ISCAT(NSCAT)
*
      REAL*8 INPROD
*
      XINOUT = 0.0D0
*
      IF ( IREW .NE. 0 ) THEN
         REWIND LU1
         REWIND LU2
         REWIND LU3
      END IF
*
* Loop over blocks
*
        IOFF = 1
 1000 CONTINUE
        IF (LBLK .GT. 0 ) THEN
          LBL1 = LBLK
          LBL2 = LBLK
        ELSE
          READ(LU1) LBL1
          READ(LU2) LBL2
          WRITE(LU3) LBL1
        END IF
        IF(LBL1 .NE. LBL2 ) THEN
          WRITE(lupri,'(A,2I3)') ' DIFFERENT BLOCKSIZES IN DMTVSD : '
     &                     , LBL1,LBL2
          Call Abend2( ' DIFFERENT BLOCKSIZES IN DMTVSD ' )
        END IF
        IF(LBL1 .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL1
          ELSE
            KBLK = -1
          END IF
          LENGTH = LBL1
          CALL FRMDSC_LUCI(VEC1,LENGTH,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC_LUCI(VEC2,LENGTH,KBLK,LU2,IMZERO,IAMPACK)
          IF( LBL1 .GT. 0 )THEN
            IF(INV .NE. 0 ) THEN
              CALL DIAVC2(VEC1,VEC2,VEC1,FAC,LENGTH)
            ELSE
              CALL VVTOV(VEC1,VEC2,VEC1,LBL1)
              CALL VECSUM(VEC1,VEC1,VEC2,1.0D0,FAC,LENGTH)
            END IF
          END IF
*
          IF(NSCAT.NE.0) THEN
            IFIRST = IOFF
            ILAST = IOFF + LENGTH - 1
            DO 100 I = 1, NSCAT
              IF(IFIRST .LE. ISCAT(I) .AND. ISCAT(I) .LE. ILAST )
     &        VEC1(ISCAT(I)-IOFF+1) = XSCAT(I)
  100       CONTINUE
          END IF
*
          XINOUT = XINOUT + DDOT(LENGTH,VEC1,1,VEC2,1)
          CALL TODSC_LUCI(VEC1,LENGTH,KBLK,LU3)
          IOFF = IOFF + LENGTH
        END IF
*
      IF( LBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      RETURN
      END
***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS(LUIN,LBLK,BLK_A,SEGMNT)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      REWIND LUIN
*
      IBLK = 0
      NBLK_A = 0
*. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LUIN) LBL
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC_LUCI(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
      NBLK =  IBLK-1
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(lupri,'(A,I8,I8)')
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(lupri,*) ' Active blocks '
        CALL WRTMT_LU(BLK_A,1,NBLK,1,NBLK)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE FIND_ACTIVE_CLASSES
     &           (LUIN,LBLK,I_BLK_TO_CLS,ICLS_A,NCLS,SEGMNT)
*
*. Find the active (nonvanishing classes ) on LUIN)
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION I_BLK_TO_CLS(*)
*. Output
      DIMENSION ICLS_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*
      REWIND LUIN
      IZERO = 0
      CALL ISETVC(ICLS_A,IZERO,NCLS)
*
      IBLK = 0
      NCLS_A = 0
*. Loop over blocks
C?      write(6,*) ' ZAP_BLOCK_VEC :  LBLK = ', LBLK
 1000 CONTINUE
        IBLK = IBLK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LUIN) LBL
        END IF
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          NO_ZEROING = 1
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
C         CALL FRMDSC_LUCI(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN  
           IF(ICLS_A(I_BLK_TO_CLS(IBLK)).EQ.0) THEN
             NCLS_A = NCLS_A + 1
           END IF
           ICLS_A(I_BLK_TO_CLS(IBLK)) = 1
C          write(6,*) ' Active block and class ',IBLK,
C    &     I_BLK_TO_CLS(IBLK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*) ' FIND_A.... Number of active classes ',NCLS_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active classes '
        CALL IWRTMA(ICLS_A,1,NCLS,1,NCLS)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE FNDMND(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &                  SCR,ISCAT,SUBVAL,NTESTG)
*
* FIND NSUB LOWEST ELEMENTS OF VECTOR RESIDING ON FILE
* LU. ENSURE THAT NO DEGENERENCIES ARE SPLIT
*
*
* INPUT
*======
* LU : FILE WHERE VECTOR OF INTEREST IS RESIDING, REWOUND
* LBLK : DEFINES FILE STRUCTURE ON FILE LU
* NSUBMX: LARGEST ALLOWED NUMBER OF SORTED ELEMENTS
*
*OUTPUT
*======
* NSUB : ACTUAL NUMBER OF ELEMENTS OBTAINED. CAN BE SMALLER
*       THAN NSUBMX IF THE LAST ELEMENT BELONGS TO A DEGENERATE
*       SET
*ISCAT: SCATTERING ARRAY, ISCAT(I) GIVES FULL ADRESS OF SORTED
*       ELEMENT I
*       SUBVAL : VALUE OF SORTED ELEMENTS

      IMPLICIT REAL*8           ( A-H,O-Z)
#include "priunit.h"
#ifdef VAR_MPI
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#endif

      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
*
      REWIND LU
*
*.LOOP OVER BLOCKS
*
C     write(6,*) ' FNDMND NSUBMX = ', NSUBMX
*
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
*
      IBLK = 0
      IBASE = 1
      LSUB = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LU) LBL
        END IF
        IBLK = IBLK + 1
        IF(NTEST.GE.10) THEN
          WRITE(lupri,*) ' Info about block ',IBLK
          WRITE(lupri,*) ' Number of elements ', LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC_LUCI(SEGMNT,LBL ,-1,LU,IMZERO,IAMPACK)
           IF(NTEST.GE.100) THEN
             WRITE(lupri,*) ' Elements read in '
             CALL WRTMT_LU(SEGMNT,1,LBL,1,LBL)
           END IF
           IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
             MSUBMX = MIN(NSUBMX,LBL)
C     SUBROUTINE SORLOW(WRK,STVAL,ISTART,KZVAR,KEXSTV,JEXSTV,IPRT)
             IF(LBL.GE.1) THEN
               CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                     MSUBMX,MSUB,NTEST)
C              WRITE(lupri,*)
C    &       ' After SORLOW 1 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
              MSUB = 0
             END IF
             DO 10 I = 1, MSUB
   10        ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
C            Write(lupri,*)
C    &       ' After 10 , Scatter array for combined list '
C            call iwrtma(ISCR,1,LSUB+MSUB,1,LSUP+MSUP)
* SORT COMBINED LIST
             MSUBMX = MIN(NSUBMX,LSUB+MSUB)
             IF(MSUBMX.GT.0) THEN
             CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &         NTEST)
C              WRITE(lupri,*)
C    &       ' After SORLOW 2 , Scatter array for combined list '
C            call iwrtma(ISCR(1),1,LSUB+MSUB,1,LSUP+MSUP)
             ELSE
               LLSUB = 0
             END IF
             LSUB = LLSUB
             DO 20 I = 1, LSUB
               ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20        CONTINUE
*
             CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
             CALL COPVEC(SUBVAL,SCR,LSUB)
             IBASE = IBASE + LBL
             IF(NTEST .GE. 20 ) THEN
               WRITE(lupri,*)' Lowest elements and their original place'
               WRITE(lupri,*)' Number of elements obtained ', LSUB
               CALL WRTMT_LU(SUBVAL,1,LSUB,1,LSUB)
               CALL IWRTMA(ISCR,1,LSUB,1,LSUB)
             END IF
          END IF
*
        END IF
*
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      NSUB = LSUB
      CALL ICOPVE(ISCR,ISCAT,NSUB)
      IF(NTEST .GE. 10) THEN
        WRITE(lupri,*) ' Lowest elements and their original place '
        WRITE(lupri,*) ' Number of elements obtained ', NSUB
        CALL WRTMT_LU(SUBVAL,1,NSUB,1,NSUB)
        CALL IWRTMA(ISCAT,1,NSUB,1,NSUB)
      END IF
*
      END
***********************************************************************
      SUBROUTINE FRMDSC_LUCI(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED)
C
C     TRANSFER ARRAY FROM DISK FILE IFILE
C
*. Version allowing zero and packed blocks 
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY(*)
*
      integer ISCR(2)
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
      integer, intent(in) :: IFILE
 
      SAVE LBATCH
      DATA LBATCH/0/
C
      IPACK = 1
      ISCR(1) = 0
      ISCR(2) = 0
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
C       IF(MMBLOCK.GE.2) MMBLOCK = 2
C       CALL IFRMDS(ISCR,2,MMBLOCK,IFILE)
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO      = ISCR(1)
        I_AM_PACKED = ISCR(2)
        IF(I_AM_PACKED.NE.0) THEN
C         WRITE(6,*) ' File is packed, file number = ', IFILE
        END IF
        IF(IMZERO.EQ.1) THEN
          CALL dzero(ARRAY,NDIM)
          GOTO 1001
        END IF
      END IF
*
      IF(I_AM_PACKED.EQ.1) THEN
        CALL dzero(ARRAY,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements 
  999   CONTINUE
          NBATCH = NBATCH + 1
          LBATCHP = LBATCH
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT 
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE    
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              Call Abend2(' problem in FRMDSC ')
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
C
      END IF
*
 1001 CONTINUE
*
      END
***********************************************************************
      SUBROUTINE FRMDSC2(ARRAY,NDIM,MBLOCK,IFILE,IMZERO,I_AM_PACKED,
     &                   NO_ZEROING)
C
C     Transfer vector from disk file IFILE to array ARRAY
C
*. Version allowing zero and packed blocks 
*
* If NO_ZEROING = 1, the elements of zero blocks
*    are not set to zero, the routine just returns with 
*    IMZERO = 1
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION ARRAY(*)
*
      integer ISCR(2)
      PARAMETER(LPBLK=50 000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)

      SAVE LBATCH
      DATA LBATCH/0/
 
      IMZERO = 0
C
      IPACK = 1
      ISCR(1) = 0
      ISCR(2) = 0
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        MMBLOCK = MBLOCK
C       IF(MMBLOCK.GE.2) MMBLOCK = 2
C       CALL IFRMDS(ISCR,2,MMBLOCK,IFILE)
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO      = ISCR(1)
        I_AM_PACKED = ISCR(2)
        IF(IMZERO.EQ.1) THEN
          IF(NO_ZEROING.EQ.0) THEN
            CALL DZERO(ARRAY,NDIM)
          END IF
          GOTO 1001
        END IF
      END IF
*
      IF(I_AM_PACKED.EQ.1) THEN
        CALL DZERO(ARRAY,NDIM)
*. Loop over packed records of dimension LPBLK
      NBATCH = 0
C1000 CONTINUE
*. The next LPBLK elements 
  999   CONTINUE
          NBATCH = NBATCH + 1
          LBATCHP = LBATCH
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE) (IPAK(I),I=1, LBATCH)
            READ(IFILE) (XPAK(I),I=1, LBATCH)
          END IF
          READ(IFILE) ISTOP
          DO IELMNT = 1, LBATCH
            IF(IPAK(IELMNT).LE.0.OR.IPAK(IELMNT).GT.NDIM) THEN
              WRITE(6,*) ' FRMDSC : Problemo IELMNT = ',IELMNT 
              WRITE(6,*) ' IPAK(IELMNT) = ',IPAK(IELMNT )
              WRITE(6,*) ' LBATCH IFILE  = ',LBATCH,IFILE    
              IF(NBATCH.EQ.1) THEN
               WRITE(6,*) ' NBATCH = 1 '
              ELSE
               WRITE(6,*) ' NBATCH, LBATCHP', NBATCH,LBATCHP
              END IF
              WRITE(6,*) ' NDIM,IMZERO = ', NDIM,IMZERO
              Call Abend2(' problem in FRMDSC ')
            END IF
            ARRAY(IPAK(IELMNT)) = XPAK(IELMNT)
          END DO
        IF(ISTOP.EQ.0) GOTO 999
*. End of loop over records of truncated elements
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE) (ARRAY(IBASE+I),I=1,NBLOCK)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE) (ARRAY(IBASE+I),I=1,IREST)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
C
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
***********************************************************************
      SUBROUTINE FRMDSCN3(VEC,NREC,LBLK,LU,NO_ZEROING,I_AM_ZERO,
     &                    LBLOCK)
*
* Read  VEC as multiple record file, NREC records read
* If NO_ZEROING.EQ.0 then zero blocks are not
* set to zero;  a 1 is instead flagged in the relevant block
* of I_AM_ZERO
*
* IF LBLOCK(IBLK) is less than zero, then no space for this
* block is made
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION LBLOCK(*)
*. Output
      DIMENSION VEC(*),I_AM_ZERO(*)
*
      IOFF = 1
C?    write(6,*) ' FRMDSCN3: Number of records to be read', NREC
      DO IREC = 1, NREC
        READ (LU) LREC
        IF(LBLOCK(IREC).GE.0) THEN
          CALL FRMDSC2(VEC(IOFF),LREC,LBLK,LU,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          I_AM_ZERO(IREC) = IMZERO
          IOFF = IOFF + LREC
        ELSE
          I_AM_ZERO(IREC) = 1
          CALL SKPRCD2(LBLK,-1,LU)
        END IF

      END DO
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GATVCD(LU,LBLK,NGAT,IGAT,XGAT,SEGMNT,IPRT)
*
* Gather elements from a file LU
*
* XGAT(I) = Vector(IGAT(I))
*
* Jeppe Olsen, September 1993
*
      IMPLICIT REAL*8           (A-H,O-Z)
*. Input
      INTEGER IGAT(NGAT)
*. Output
      DIMENSION XGAT(NGAT)
*. Scratch
      DIMENSION SEGMNT(*)
*
      REWIND LU
*
      IBASE = 1
      IBLOCK = 0
*
*. Loop over blocks of file
*
 1000 CONTINUE
        IBLOCK = IBLOCK + 1
        CALL NEXREC(LU,LBLK,SEGMNT,IEND,LENGTH)
        IF(IPRT.GE.10)
     &  WRITE(6,*) LENGTH, ' elements in block ',IBLOCK
        IF(IEND.EQ.0) THEN
          IFIRST = IBASE
          ILAST = IBASE + LENGTH - 1
          DO 100 I = 1, NGAT
            IF(IFIRST .LE. IGAT(I) .AND. IGAT(I) .LE. ILAST )
     &      XGAT(I) = SEGMNT(IGAT(I)-IFIRST+1)
C?          IF(IFIRST .LE. IGAT(I) .AND. IGAT(I) .LE. ILAST )
C?   &      write(6,*) ' Catch I IGAT(I) XGAT(I) ',
C?   &                         I,IGAT(I),XGAT(I)
  100     CONTINUE
          IBASE = IBASE + LENGTH
      IF(LBLK.LT.0) GOTO 1000
        END IF
*
      NTEST = 0
      NTEST = MAX(IPRT,NTEST)
      IF(NTEST.GE.5) THEN
       WRITE(6,*) ' Gathered vector from GATVCD '
       CALL WRTMT_LU(XGAT,1,NGAT,1,NGAT)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GET_TTS_BATCHN(CTTS,NBLOCK,IBLOCK,NBLOCKC,IBLOCKC,
     &                  NOCTPA,NOCTPB,NSMST,NSASO,NSBSO,
     &                  IDC,LUC,IREW,ISCALE)
*
* Read in a batch of C blocks from file LUC.
*
* The complete file is defined by NBLOCKC,IBLOCKC,
* and the blocks of the actual batch is defined by NBLOCK,IBLOCK.
* Vector packing is defined by IDC
*
* Should be initialized with rewind on LUC i.e. IREW = 1
*
*. Feb 96 : Modified to accomodate packed files
*. March 97 : Shaved, assumes IBLOCK and IBLOCKC have identical ordering
*
*
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CTTS(*),NSASO(NSMST,*),NSBSO(NSMST,*)
      DIMENSION IBLOCKC(8,NBLOCKC),IBLOCK(8,NBLOCK)
*
      COMMON/HIDLUC/IBLK,IATP,IBTP,IASM,IBSM
*
      NTEST = 000
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Welcome to GET_TTS_BATCH '
        WRITE(6,*) ' ========================='
        WRITE(6,*) ' Number of blocks in batch ', NBLOCK
        WRITE(6,*)
        WRITE(6,*) ' IATP  IBTP  IASM  IBSM '
        WRITE(6,*) ' ====================== '
        DO JBLOCK = 1, NBLOCK
          WRITE(6,'(4I4)') (IBLOCK(II,JBLOCK),II = 1, 4 )
        END DO
      END IF
*
      IF(IREW.EQ.1) THEN
        REWIND (LUC)
        IBLK = 1
      END IF
*
      IAMPACK = 1
      IMZERO= 0
*. Ready to loop over blocks
*
*. Loop over blocks to be read
      DO JBLOCK = 1, NBLOCK
          IF(IBLOCK(1,JBLOCK).GT.0) THEN
            ISKIP = 0
          ELSE
            ISKIP = 1
          END IF
*
          IF(IBLK.GT.NBLOCKC) THEN
            REWIND LUC
            WRITE(6,*) ' Notice : LUC rewinded in GET_TTS'
            WRITE(6,*) ' Less than optimal programming ? '
            IBLK = 1
          END IF
*. Loop over blocks on file
  999   CONTINUE
*. is this block identical to next block on disk
          READ (LUC) LBL
          IF(ABS(IBLOCK(1,JBLOCK)).EQ.ABS(IBLOCKC(1,IBLK)) .AND.
     &       IBLOCK(2,JBLOCK).EQ.IBLOCKC(2,IBLK) .AND.
     &       IBLOCK(3,JBLOCK).EQ.IBLOCKC(3,IBLK) .AND.
     &       IBLOCK(4,JBLOCK).EQ.IBLOCKC(4,IBLK) ) THEN
             IF(ISKIP.EQ.0) THEN
*. Read me !
              IOFFO   = IBLOCK(6,JBLOCK)
              CALL FRMDSC_LUCI(CTTS(IOFFO),LBL,-1,LUC,IMZERO,IAMPACK)
              IBLK = IBLK + 1
            ELSE
*. Skip me, you are not interested in me at all
              CALL SKPRCD2(LBL,-1,LUC)
              IBLK = IBLK + 1
            END IF
          ELSE
*. Skip me
            CALL SKPRCD2(LBL,-1,LUC)
            IBLK = IBLK + 1
            GOTO 999
          END IF
C     END IF
*. End of loop over blocks on file
      END DO
*
*. Rescale from combination form to determinant
*
      IF(IDC.EQ.2.AND.ISCALE.NE.0) THEN
        CALL SCDTTS(CTTS,IBLOCK,NBLOCK,
     &  NSMST,NOCTPA,NOCTPB,NSASO,NSBSO,
     &  IDC,2,NTEST)
C     SCDTTS(BLOCKS,IBLOCK,NBLOCK,
C    &                  NSMST,NOCTPA,NOCTPB,
C    &                  NSASO,NSBSO,IDC,IWAY,IPRNT)
      END IF
*
      IF(NTEST.GE.1000) THEN
         WRITE(6,*) ' Batch of blocks read in '
         WRITE(6,*) ' ========================'
         WRITE(6,*)
         CALL WRTTTS(CTTS,IBLOCK,NBLOCK,NSMST,NOCTPA,NOCTPB,
     &               NSASO,NSBSO,IDC)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE GFLUSH(LU)
*
*. Outer routine for flushing UNITS !
*
*. RISC/6000
C      CALL FLUSH_(LU)
*. CRAY.
C      CALL FLUSH(LU)
*
      RETURN
      END
***********************************************************************
      SUBROUTINE IFRMDS(IARRAY,NDIM,MBLOCK,IFILE)
C
C     TRANSFER INTEGER ARRAY FROM DISK FILE IFILE
C
C If nblock .eq. 0 NBLOCK = NDIM
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IARRAY(*)
C
      NBLOCK = MBLOCK

        IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
          IF(IREST.GT.NBLOCK) THEN
            READ(IFILE) (IARRAY(IBASE+I),I=1,NBLOCK)
            IBASE=IBASE+NBLOCK
            IREST=IREST-NBLOCK
          ELSE
            READ(IFILE) (IARRAY(IBASE+I),I=1,IREST)
            IREST=0
          END IF
        IF( IREST .GT. 0 ) GOTO 100
      RETURN
      END
***********************************************************************
      REAL*8 FUNCTION INPRDD(VEC1,VEC2,LU1,LU2,IREW,LBLK)
C
C DISK VERSION OF INPROD
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 INPROD
      DIMENSION VEC1(*),VEC2(*)
      LOGICAL DIFVEC
      integer, intent(in) :: LU1, LU2, IREW, LBLK
      integer :: myLU1, myLU2, NBL1, NBL2
C
      X = 0.0D0
      myLU1 = LU1
      myLU2 = LU2

      DIFVEC =  .FALSE.
      IF(myLU1.ne.myLU2) DIFVEC = .TRUE.
C
      IF( IREW .NE. 0 ) THEN
        REWIND myLU1
        IF(DIFVEC) REWIND myLU2
      END IF
C
C LOOP OVER BLOCKS OF VECTORS
C
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE
          NBL1 = 0
          NBL2 = 0
          READ(myLU1) NBL1
          IF( DIFVEC) READ(myLU2) NBL2
        END IF
C
        IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
          CALL FRMDSC_LUCI(VEC1,NBL1,KBLK,myLU1,IMZERO,IAMPACK)
          IF( DIFVEC) THEN
            CALL FRMDSC_LUCI(VEC2,NBL1,KBLK,myLU2,IMZERO,IAMPACK)
            IF(NBL1 .GT. 0 )
     &      X = X + DDOT(NBL1,VEC1,1,VEC2,1)
C?          write(6,*) ' vec1 and vec2 in INPRDD '
C?         CALL WRTMT_LU(VEC1,1,NBL1,1,NBL1)
C?         CALL WRTMT_LU(VEC2,1,NBL1,1,NBL1)
          ELSE
          IF(NBL1 .GT. 0 )
     &    X = X + DDOT(NBL1,VEC1,1,VEC1,1) 
        END IF
      END IF
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      INPRDD = X
C
      END
***********************************************************************
      SUBROUTINE ITODS(IA,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY INTEGER IA(LENGTH NDIM) TO DISK FILE IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION IA(*)
      INTEGER START,STOP
*
      NBLOCK = MBLOCK
C
      IF(NBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (IA(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
C
C
      RETURN
      END
***********************************************************************
      SUBROUTINE MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
C
C ADD VECTORS ON FILE LUIN TIMES FACTOR AND STORE ON LUOUT
C
C LUOUT AND LUSCR ARE INITIALLY REWOUND
C
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(NVEC)
C
      IF( MOD(NVEC,2) .EQ. 0 ) THEN
        LLUOUT = LUSCR
        LLUSCR = LUOUT
      ELSE
        LLUOUT = LUOUT
        LLUSCR = LUSCR
      END IF
C
      IF(IREW .NE. 0 ) REWIND LUIN
C
      DO 100 IVEC = 1, NVEC
        REWIND LLUSCR
        REWIND LLUOUT
        IF( IVEC .EQ. 1 ) THEN
          CALL SCLVCD(LUIN,LLUOUT,FAC(IVEC),VEC1,0,LBLK)
        ELSE
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),1.0D0,LUIN,LLUSCR,LLUOUT,
     &                0,LBLK)
        END IF
C
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
  100 CONTINUE
C
      RETURN
      END
***********************************************************************
      SUBROUTINE MVCSMD2(LUIN,FAC,FACLUOUT,LUOUT,LUSCR,
     &                   VEC1,VEC2,NVEC,IREW,LBLK)
*
* LUOUT = Factor * LUOUT + sum(I=1,nvec)fac(I) LUIN(I)
*
* LUOUT and LUSCR are rewinded and input data 
* on these files are lost.
*
* MVCSMD2 is identical MVCSMD, except that
* the input vector on LUOUT can be included
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(1),VEC2(1)
      DIMENSION FAC(1)
*
      IF(IREW .NE. 0 ) REWIND LUIN
*
      LLUOUT = LUOUT
      LLUSCR = LUSCR
      DO IVEC = 1, NVEC
        REWIND LLUSCR
        REWIND LLUOUT
        IF( IVEC .EQ. 1 ) THEN
          IF(FACLUOUT.NE.0.0D0) THEN
            CALL VECSMD(VEC1,VEC2,FAC(1),FACLUOUT,LUIN,
     &                  LUOUT,LUSCR,0,LBLK)
          ELSE
            CALL SCLVCD(LUIN,LUSCR,FAC(1),VEC1,0,LBLK)
          END IF
        ELSE
          ONE = 1.0D0
          CALL VECSMD(VEC1,VEC2,FAC(IVEC),ONE,LUIN,LLUOUT,LLUSCR,
     &                0,LBLK)
        END IF
*
        LBUF = LLUOUT
        LLUOUT = LLUSCR
        LLUSCR = LBUF
      END DO
*
      IF(LLUOUT.NE.LUOUT) THEN
C            COPVCD(LUIN,LUOUT,SEGMNT,IREW,LBLK)
        CALL COPVCD(LLUOUT,LUOUT,VEC1,1,LBLK)
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE NEXREC(LU,LBLK,REC,IEND,LENGTH)
*
* OBTAIN NEXT RECORD ON FILE LU, IF
* AN END OF FILE IS ISSUED THE RECORD IS EMPTY
* AND IEND IS SET TO 1
*
      IMPLICIT REAL*8           (A-H,O-Z)
      DIMENSION REC(*)
*
      IF(LBLK .GT. 0 ) THEN
        IEND = 0
        LENGTH = LBLK
        CALL FRMDSC_LUCI(REC,LENGTH,-1,LU,IMZERO,IAMPACK)
      ELSE
        READ (LU) LENGTH
        IF(LENGTH.GE.0) THEN
          IEND = 0
          CALL FRMDSC_LUCI(REC,LENGTH,LBLK,LU,IMZERO,IAMPACK)
        ELSE
          IEND = 1
        END IF
      END IF
*
      RETURN
      END
***********************************************************************
      SUBROUTINE REWINE( LU ,LBLK )
* 
* LBLK used to select FASTIO or not, is now ignored.
      REWIND LU
*
      RETURN
      END
***********************************************************************
      SUBROUTINE SCLVCD(LUIN,LUOUT,SCALE,SEGMNT,IREW,LBLK)
C
C SCALE VECTOR ON FILE LUIN WITH FACTOR SCALE AND STORE ON LUOUT
C
C
C LBLK DEFINES STRUCTURE OF FILES
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) THEN
        REWIND LUIN
        REWIND LUOUT
      END IF
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LUIN) LBL
          WRITE(LUOUT) LBL
        END IF
C
        IF ( LBL .GE. 0 ) THEN
          IF(      LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C
          CALL FRMDSC_LUCI(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(LBL .GT. 0 )
     &    CALL DSCAL(LBL,SCALE,SEGMNT,1)
          CALL TODSC_LUCI(SEGMNT,LBL,KBLK,LUOUT)
        END IF
C
      IF( LBL .GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
***********************************************************************
      SUBROUTINE SKPRCD2(NDIM,MBLOCK,IFILE)
C
C     Skip record in file IFILE
C
*. Version allowing zero and packed blocks
*
      IMPLICIT REAL*8(A-H,O-Z)
*
      DIMENSION ISCR(2)
      PARAMETER(LPBLK=50000)

C
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Read if ARRAY is zero
        CALL IFRMDS(ISCR,2,2,IFILE)
        IMZERO=ISCR(1)
        I_AM_PACKED=ISCR(2)
CSK        WRITE(6,*) 'IMZERO is',IMZERO
        IF(IMZERO.EQ.1) THEN
          GOTO 1001
        END IF
      END IF
*
      IF(I_AM_PACKED.EQ.1) THEN
*. Loop over packed records of dimension LPBLK
*. The next LPBLK elements
  999   CONTINUE
*. Read next batch
          READ(IFILE) LBATCH
          IF(LBATCH.GT.0) THEN
            READ(IFILE)
            READ(IFILE)
          END IF
          READ(IFILE) ISTOP
        IF(ISTOP.EQ.0) GOTO 999
      ELSE IF ( I_AM_PACKED.EQ.0) THEN
        NBLOCK = MBLOCK
        IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
        IREST=NDIM
        IBASE=0
  100   CONTINUE
         IF(IREST.GT.NBLOCK) THEN
          READ(IFILE)
          IBASE=IBASE+NBLOCK
          IREST=IREST-NBLOCK
         ELSE
          READ(IFILE)
          IREST=0
         END IF
        IF( IREST .GT. 0 ) GOTO 100
C
      END IF
*
 1001 CONTINUE
*
      RETURN
      END
***********************************************************************
      SUBROUTINE SKPVCD(LU,NVEC,SEGMNT,IREW,LBLK)
C
C SKIP OVER NVEC VECTORS ON FILE LUIN
C
C LBLK DEFINES STRUCTURE OF FILE
C (see note on structure of files )
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) REWIND LU
      DO 1001 IVEC = 1, NVEC
C
C LOOP OVER BLOCKS OF GIVEN VECTOR
C
 1000   CONTINUE
C
          IF( LBLK .GT. 0 ) THEN
            LBL = LBLK
          ELSE
            READ(LU) LBL
          END IF
C
          IF( LBL .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = LBLK
            ELSE
              KBLK = -1
            END IF
            CALL FRMDSC_LUCI(SEGMNT,LBL,KBLK,LU,IMZERO,IAMPACK)
          END IF
        IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
 1001 CONTINUE
C
      RETURN
      END
***********************************************************************
      SUBROUTINE TODSC_LUCI(A,NDIM,MBLOCK,IFIL)
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO disk file IFIL IN
C RECORDS WITH LENGTH NBLOCK.
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
*
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = DDOT(NDIM,A,1,A,1)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. No packing 
        ISCR(2) = 0
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (A(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
*
 1001 CONTINUE
C
      RETURN
      END
***********************************************************************
      SUBROUTINE TODSC2(A,NDIM,MBLOCK,IFIL,IACTIVE)
C
C     TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO disk file IFIL IN
C     RECORDS WITH LENGTH NBLOCK.
C
C     'active' blocks are flagged by a nonzero entry for iactive
C
C     Last revision:     S. Knecht       - April  2007
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(*)
      INTEGER START,STOP
      INTEGER ISCR(2)
*
      IPACK = 1
      IF(IPACK.NE.0) THEN
        IF(IACTIVE .eq. 0 ) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. No packing 
         ISCR(2) = 0
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
      STOP=0
      NBACK=NDIM
C LOOP OVER RECORDS
  100 CONTINUE
       IF(NBACK.LE.NBLOCK) THEN
         NTRANS=NBACK
         NLABEL=-NTRANS
       ELSE
         NTRANS=NBLOCK
         NLABEL=NTRANS
       END IF
       START=STOP+1
       STOP=START+NBLOCK-1
       NBACK=NBACK-NTRANS
       WRITE(IFIL) (A(I),I=START,STOP),NLABEL
      IF(NBACK.NE.0) GOTO 100
*
 1001 CONTINUE
C
      END
***********************************************************************
      SUBROUTINE TODSCP(A,NDIM,MBLOCK,IFIL)
*
C TRANSFER ARRAY DOUBLE PRECISION  A(LENGTH NDIM) TO disk file IFIL IN
C RECORDS WITH LENGTH NBLOCK.
*
* Packed version : Store only nonzero elements
*. Small elements should be xeroed outside
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(1)
      INTEGER START,STOP
      REAL*8 INPROD
      INTEGER ISCR(2)
* 
      PARAMETER(LPBLK=50000)
      INTEGER IPAK(LPBLK)
      DIMENSION XPAK(LPBLK)
*
*
C?    write(6,*) ' entering TODSCP, file = ', IFIL
C?    CALL XFLUSH(6)
      IPACK = 1
      IF(IPACK.NE.0) THEN
*. Check norm of A before writing
        XNORM = DDOT(NDIM,A,1,A,1)
        IF(XNORM.EQ.0.0D0) THEN
          IMZERO = 1
        ELSE
          IMZERO = 0
        END IF
        MMBLOCK = MBLOCK
        IF(MMBLOCK.GT.2) MMBLOCK = 2
*
        ISCR(1) = IMZERO
*. Packing 
        ISCR(2) = 1
C       CALL ITODS(ISCR,2,MMBLOCK,IFIL)
        CALL ITODS(ISCR,2,2,IFIL)
        IF(IMZERO.EQ.1) GOTO 1001
      END IF
*
      NBLOCK = MBLOCK
      IF ( MBLOCK .LE. 0 ) NBLOCK = NDIM
*. Loop over packed records of dimension LPBLK
      IELMNT = 0
 1000 CONTINUE
*. The next LPBLK elements 
      LBATCH = 0
*. Obtain next batch of elemnts
  999 CONTINUE
       IF(NDIM.GE.1) THEN
       IELMNT = IELMNT+1
       IF(A(IELMNT).NE.0.0D0) THEN
         LBATCH=LBATCH+1
         IPAK(LBATCH) = IELMNT
         XPAK(LBATCH) = A(IELMNT)
       END IF
       END IF
       IF(LBATCH.EQ.LPBLK.OR.IELMNT.EQ.NDIM) goto 998
       GOTO 999
*. Send to DISK
 998   CONTINUE   
       WRITE(IFIL) LBATCH
       IF(LBATCH.GT.0) THEN 
         WRITE(IFIL) (IPAK(I),I=1, LBATCH)
         WRITE(IFIL) (XPAK(I),I=1, LBATCH)
       END IF
       IF(IELMNT.EQ.NDIM) THEN
         WRITE(IFIL) -1
       ELSE
         WRITE(IFIL) 0
         GOTO 1000
       END IF
*. End of loop over records of truncated elements
 1001 CONTINUE
*
      RETURN
      END
***********************************************************************
      SUBROUTINE TRAVCD(VEC1,VEC2,X,NVECIN,NVECOUT,LUIN,LUOUT,
     &                  ICOPY,LBLK,LUSCR1,LUSCR2)
*
* NVECIN vectors reside on LU1, Transform these vectors,
*  using LUSCR1 and  LUSCR2  as
* scratch files,  with matrix X to produce output file LUOUT.
*
* Since LUIN is accessed several times it is always
* REWOUND. LUOUT is written to from current start.
*
* I ICOPY .ne. 0 the transformed vectors are written back to LUIN
*
* Jeppe Olsen, April 1997
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Input
      DIMENSION X(NVECIN,NVECOUT)
*. Scratch
      DIMENSION VEC1(*),VEC2(*)
*
C             MVCSMD(LUIN,FAC,LUOUT,LUSCR,VEC1,VEC2,NVEC,IREW,LBLK)
      DO IVECOUT = 1, NVECOUT
         CALL MVCSMD(LUIN,X(1,IVECOUT), LUSCR1,LUSCR2,VEC1,VEC2,
     &               NVECIN,1,LBLK)
         REWIND LUSCR1
         CALL COPVCD(LUSCR1,LUOUT,VEC1,0,LBLK)
      END DO
*
      IF(ICOPY.EQ.1) THEN
        REWIND LUIN
        REWIND LUOUT
        DO IVECOUT = 1, NVECOUT
          CALL COPVCD(LUOUT,LUIN,VEC1,0,LBLK)
        END DO
      END IF
*
      RETURN
      END
***********************************************************************
      REAL*8 FUNCTION VCSMDN(VEC1,VEC2,FAC1,FAC2,LU1,LU2,IREW,LBLK)
*
* Norm of sum of two vectors residing on disk
*
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
*
      XNORM = 0.0D0
      IF(IREW .NE. 0 ) THEN
        REWIND  LU1
        REWIND  LU2
      END IF
*
* LOOP OVER BLOCKS OF VECTOR
*
 1000 CONTINUE
C
        IF( LBLK .GT. 0 ) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE
          READ(LU1) NBL1
          READ(LU2) NBL2
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
        WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &  NBL1,NBL2
        Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
      END IF
C
      IF(NBL1 .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = NBL1
          ELSE
            KBLK = -1
          END IF
        CALL FRMDSC_LUCI(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
        CALL FRMDSC_LUCI(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
        IF( NBL1 .GT. 0 ) THEN
          CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          XNORM = XNORM + DDOT(NBL1,VEC1,1,VEC1,1)
        END IF
      END IF
*
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
*
      VCSMDN = XNORM
      RETURN
      END
***********************************************************************
      SUBROUTINE VECSMD(VEC1,VEC2,FAC1,FAC2,LU1,LU2,LU3,IREW,LBLK)
C
C DISK VERSION OF VECSUM :
C
C      ADD BLOCKED VECTORS ON FILES LU1 AND LU2
C      AND STORE ON LU3
C
C LBLK DEFINES STRUCTURE OF FILE
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION VEC1(*),VEC2(*)
C
      IF(IREW .NE. 0 ) THEN
        REWIND  LU1
        REWIND  LU2
        REWIND  LU3
      END IF
C
C LOOP OVER BLOCKS OF VECTOR
C
 1000 CONTINUE
C
        IF (LBLK .GT. 0) THEN
          NBL1 = LBLK
          NBL2 = LBLK
        ELSE
          READ(LU1) NBL1
          READ(LU2) NBL2
          WRITE(LU3) NBL1
        END IF
        IF( NBL1 .NE. NBL2 ) THEN
          WRITE(6,'(A,2I5)') 'DIFFERENT BLOCKSIZES IN VECSMD ',
     &    NBL1,NBL2
          Call Abend2( ' INCOMPATIBLE BLOCKSIZES IN VECSMF ' )
        END IF
C
        IF(NBL1 .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = NBL1
            ELSE
              KBLK = -1
            END IF
          CALL FRMDSC_LUCI(VEC1,NBL1,KBLK,LU1,IMZERO,IAMPACK)
          CALL FRMDSC_LUCI(VEC2,NBL1,KBLK,LU2,IMZERO,IAMPACK)
          IF( NBL1 .GT. 0 )
     &    CALL VECSUM(VEC1,VEC1,VEC2,FAC1,FAC2,NBL1)
          IF(IAMPACK.EQ.0) THEN
            CALL TODSC_LUCI(VEC1,NBL1,KBLK,LU3)
          ELSE
            CALL TODSCP(VEC1,NBL1,KBLK,LU3)
          END IF
        END IF
C
      IF(NBL1.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
***********************************************************************
      SUBROUTINE WRTVCD(SEGMNT,LU,IREW,LBLK)
C
C PRINT VECTOR ON FILE LU
C
C LBLK DEFINES STRUCTURE OF FILES :
C
      IMPLICIT REAL*8(A-H,O-Z)
#include "priunit.h"
      DIMENSION SEGMNT(*)
C
      IF( IREW .NE. 0 ) REWIND LU
C LOOP OVER BLOCKS
C
      IBLK = 0
 1000 CONTINUE
        IF ( LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE
          READ(LU) LBL
        END IF
        IBLK = IBLK + 1
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
           CALL FRMDSC_LUCI(SEGMNT,LBL ,KBLK,LU,IMZERO,IAMPACK)
           IF(LBL .GT. 0 ) THEN
             WRITE(lupri,'(A,I3,A,I6)')
     &       ' Number of elements in segment ',IBLK,' IS ',LBL
             CALL WRTMT_LU(SEGMNT,1,LBL,1,LBL)
           END IF
        END IF
C
      IF( LBL.GE. 0 .AND. LBLK .LE. 0) GOTO 1000
C
      RETURN
      END
***********************************************************************
      SUBROUTINE ZERORC(MBLOCK,IFIL,IAMPACKED)
*
* A record was known to be identical zero
*
* Write corresponding info to file IFIL
*
* IAMPACKED added Oct. 98 / Jeppe Olsen
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER ISCR(2)
* Zero record 
      ISCR(1) = 1
*. Packed form
      ISCR(2) = IAMPACKED
*
      CALL ITODS(ISCR,2,2,IFIL)
*
      RETURN 
      END
