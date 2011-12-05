!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!

#ifdef VAR_MPI
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_RL_LUCI1(VEC1,VEC2,EIGAPR,RNRM,EIGSHF,
     &                             EIG,TEST,E_CONV,RTCNV,CONVER,ITER,
     &                             MAXIT,
     &                             IROOT,LUIN1LIST,LUIN2LIST,LUOUTLIST,
     &                             NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                             MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                             MY_IOFF_LUOUT,
     &                             SCRRED,LUIN1,LUIN2,LUOUT)
C
C     Written by  S. Knecht        - March 13 2008
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - March 2008
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUOUTLIST(*)
      DIMENSION RNRM(MAXIT,*), EIG(MAXIT,*)
      DIMENSION SCRRED(*)
      LOGICAL CONVER, RTCNV(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
      RNORM     = 0.0D0
      REDSCRVAR = 0.0D0
      CALL DZERO(SCRRED,IROOT)
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
        FACTOR = - EIGAPR
C
C       set new offset wrt IROOT
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )  * MY_VEC1_IOFF
        IOFFSET_INT_IN2  = 1 + NUM_BLK   +
     &                   ( IROOT - 1 )  * MY_ACT_BLK1
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK  &                 IOFFSET_IN_LUIN2
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       set offset for LUIN1 and zero read-in vector
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        IOFFSET_IN_LUIN1 = MY_IOFF_LUIN1 + IOFFSET_SCRATCH +
     &                   ( IROOT - 1 )    * MY_VEC1_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK   + 
     &                   ( IROOT - 1 )    * MY_ACT_BLK1
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN1',
CSK  &                  IOFFSET_IN_LUIN1
C
C       read in batch ISBATCH from LUIN1 to VEC2
C
        CALL RDVEC_BATCH_DRV4(LUIN1,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                        LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN1'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       VEC2 == VEC2 + VEC1 * FACTOR 
C
        CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
C       calculate partial RNORM
C
        REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC2,1)
C
C       write VEC2 to LUOUT
C
        IOFFSET_LUOUT     = MY_IOFF_LUOUT  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT = 1 + NUM_BLK
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUOUT',
CSK  &                  IOFFSET_LUOUT
CSK     WRITE(LUWRT,*) 'final VEC2 to write on LUOUT'
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
C     ^ loop over batches
C
C     communicate REDSCRVAR to get full RNORM
C
      CALL REDVEC_REL(REDSCRVAR,SCRRED,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
      CALL DCOPY(1,SCRRED,1,REDSCRVAR,1)
C
      RNORM = SQRT(REDSCRVAR)
C
      RNRM(ITER-1,IROOT) = RNORM
C
C     print norm and eigenvalue
C
      WRITE(LUWRT,'(A19,I8,1p,1E15.5,0p,3X,1F19.10)')
     &     ' Iter RNORM EIGAPR ', ITER-1,RNORM,EIGAPR+EIGSHF
      CALL FLSHFO(LUWRT)
C
C
C     screening of new trial vector
C
      RNORM_FAC = RNORM
C
      IF (TRUNC_FAC .GT. 0.1D0) THEN
          WRITE (LUWRT,*) 'TRUNC_FAC reset from ',TRUNC_FAC,' to',0.1D0
          TRUNC_FAC = 0.1D0
      END IF
C
C     check for convergence
C
      IF(RNORM .lt. TEST .OR. ( ITER .gt. 2 .and. 
     &  ABS(EIG(ITER-2,IROOT)-EIG(ITER-1,IROOT)).LT.E_CONV)) THEN
C
        RTCNV(IROOT) = .TRUE.
      ELSE
C
        RTCNV(IROOT) = .FALSE.
        CONVER       = .FALSE.
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_RL_LUCI2(VEC1,VEC2,SHIFT,IROOT,
     &                             LUINLIST,LUOUT1LIST,LUOUT2LIST,
     &                             LUOUT3LIST,
     &                             NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                             MY_IOFF_LUIN,MY_IOFF_LUOUT1,
     &                             MY_IOFF_LUOUT2,MY_IOFF_LUOUT3,
     &                             MY_IOFF_LUDIA,
     &                             LUIN,LUOUT1,LUOUT2,LUOUT3,LUDIA,INV)
C
C     Written by  S. Knecht         - March 13 2008
C
C**********************************************************************
C
C     part 1.2 of DAVIDSON-OLSEN algorithm in MPI-file I/O mode
C
C     NOTE: IROOT = IROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUT1LIST(*), LUOUT2LIST(*)
      DIMENSION LUOUT3LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT3
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUDIA
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN, IOFFSET_LUOUT1
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUDIA
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK            = 0
      IOFFSET_SCRATCH    = 0
      IOFFSET_IN_LUIN    = 0
      IOFFSET_LUOUT1     = 0
      IOFFSET_LUOUT2     = 0
      IOFFSET_LUOUT3     = 0
      IOFFSET_IN_LUDIA   = 0
      IOFFSET_INT_IN     = 0
      IOFFSET_INT_LUOUT1 = 0
      IOFFSET_INT_LUOUT2 = 0
      IOFFSET_INT_LUOUT3 = 0
C
      REDSCRVAR = 0.0D0
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
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
C       calculate inverse diagonal on VEC1
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C       new offset for file containing diagonal
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
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
        IF( LEBATCH(ISBATCH) .gt. 0 )THEN
          IF( CSCREEN) THEN
C           set proper truncation factor
            THR_TRUNC  = TRUNC_FAC * RNORM_FAC
Csk         WRITE(LUWRT,*) 'TRUNCATION FACTOR:',THR_TRUNC
Chj         14-jun-07:   disable THR_ETRUNC
Chj         THR_ETRUNC = 1.0D-6 * THRES_E
            THR_ETRUNC = -1.0D0
            CALL DIAVC2_TRUNC(VEC1,VEC2,VEC1,SHIFT,LEBATCH(ISBATCH),
     &                        THR_TRUNC,THR_ETRUNC)
          ELSE
            CALL DIAVC2(VEC1,VEC2,VEC1,SHIFT,LEBATCH(ISBATCH))
          END IF
        END IF
C
C       write VEC1 to LUOUT1 and VEC2 to LUOUT2
C
        IOFFSET_LUOUT1     = MY_IOFF_LUOUT1  + IOFFSET_SCRATCH
        IOFFSET_LUOUT2     = MY_IOFF_LUOUT2  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT1 = 1 + NUM_BLK
        IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
C       VEC1
        CALL WTVEC_BATCH_DRV4(LUOUT1,VEC1,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT1,IOFFSET_INT_LUOUT1,
     &                       LUOUT1LIST,NUM_ACTIVE_BATCH)
C       VEC2
        CALL WTVEC_BATCH_DRV4(LUOUT2,VEC2,LBATCH(ISBATCH),
     &                       IBATCH(1,I1BATCH(ISBATCH)),
     &                       IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                       LUOUT2LIST,NUM_ACTIVE_BATCH)
C
C
CSK     WRITE(LUWRT,*) 'THIS IS my partial REDSCRVAR',REDSCRVAR
CSK     WRITE(LUWRT,*) 'THIS IS LEBATCH(ISBATCH)',LEBATCH(ISBATCH)
C       calculate partial GAMMA
        IF( LEBATCH(ISBATCH) .gt. 0 ) THEN
          REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC1,1)
        END IF
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C     communicate REDSCRVAR to get full GAMMA
C
      CALL REDVEC_REL(REDSCRVAR,GAMMA,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
C
CSK   WRITE(LUWRT,*) 'THIS IS GAMMA',GAMMA
C
C     continue with VNORM ...
C
C     reset scratch offsets
      NUM_BLK            = 0
      IOFFSET_SCRATCH    = 0
      REDSCRVAR = 0.0D0
CSK   WRITE(LUWRT,*) 'THIS IS REDSCRVAR',REDSCRVAR
C
      DO ISBATCH = 1, NBATCH
C
        CALL DZERO(VEC1,LEBATCH(ISBATCH))
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       read VEC1 from LUOUT2 and VEC2 from LUOUT1
C
        IOFFSET_LUOUT1     = MY_IOFF_LUOUT1  + IOFFSET_SCRATCH
        IOFFSET_LUOUT2     = MY_IOFF_LUOUT2  + IOFFSET_SCRATCH
        IOFFSET_INT_LUOUT1 = 1 + NUM_BLK
        IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
C       VEC1
        CALL RDVEC_BATCH_DRV4(LUOUT2,VEC1,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                        LUOUT2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC1 on LUOUT2 in P1..._2 again'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C       VEC2
        CALL RDVEC_BATCH_DRV4(LUOUT1,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_LUOUT1,IOFFSET_INT_LUOUT1,
     &                        LUOUT1LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) ' VEC2 before DAXPY call'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
CSK     WRITE(LUWRT,*) ' VEC1 before DAXPY call'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       VEC2 == VEC2 + VEC1 * FACTOR
C
        CALL DAXPY(LEBATCH(ISBATCH),-GAMMA,VEC1,1,VEC2,1)
C
CSK     WRITE(LUWRT,*) ' VEC2 after DAXPY call'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
CSK     WRITE(LUWRT,*) ' VEC1 after DAXPY call'
CSK     CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       calculate partial VNORM_Q
C
        REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC2,1,VEC2,1)
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      VNORM_Q = 0.0D0
      VNORM   = 0.0D0
C
C     communicate REDSCRVAR to get full VNORM_Q
C
      CALL REDVEC_REL(REDSCRVAR,VNORM_Q,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
C
C     is X an eigen vector for (H0 - 1 ) - 1 ???
C
      VNORM = SQRT(VNORM_Q)
C
csk   WRITE(LUWRT,*) 'GAMMA ',GAMMA
csk   WRITE(LUWRT,*) 'VNORM ',VNORM
C
      IF( VNORM .GT. 1.0D-7 ) THEN
        IOLSAC = 1
      ELSE
        IOLSAC = 0
      END IF
      IF(IOLSAC .EQ. 1 ) THEN
C
CSK     WRITE(LUWRT,*) ' Olsen correction active'
        DELTA = 0.0D0 
C
C       continue with DELTA ...
C
C       reset scratch offsets
        NUM_BLK            = 0
        IOFFSET_SCRATCH    = 0
        REDSCRVAR = 0.0D0
C
        DO ISBATCH = 1, NBATCH
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
          CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C         read VEC1 from LUOUT2 and VEC2 from LUOUT3
C
          IOFFSET_LUOUT3     = MY_IOFF_LUOUT3  + IOFFSET_SCRATCH
          IOFFSET_LUOUT2     = MY_IOFF_LUOUT2  + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT3 = 1 + NUM_BLK
          IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
C         VEC1
          CALL RDVEC_BATCH_DRV4(LUOUT2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                          LUOUT2LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUOUT2'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C         VEC2
          CALL RDVEC_BATCH_DRV4(LUOUT3,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT3,IOFFSET_INT_LUOUT3,
     &                          LUOUT3LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC2 on LUOUT3'
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
C         calculate partial DELTA
C
          REDSCRVAR = REDSCRVAR + DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
C         keep track of correct offset
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
          NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
        END DO
C
C
C       communicate REDSCRVAR to get full DELTA
C
        CALL REDVEC_REL(REDSCRVAR,DELTA,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
C
C
CSK     WRITE(LUWRT,*) ' THIS IS DELTA'
C
        FACTOR = - DELTA / GAMMA
csk     WRITE(LUWRT,*) 'FACTOR, DELTA, GAMMA',FACTOR, DELTA, GAMMA
C
C       reset scratch offsets
        NUM_BLK            = 0
        IOFFSET_SCRATCH    = 0
C
        DO ISBATCH = 1, NBATCH
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
          CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C         read VEC1 from LUOUT1 and VEC2 from LUOUT3
C
          IOFFSET_LUOUT3     = MY_IOFF_LUOUT3  + IOFFSET_SCRATCH
          IOFFSET_LUOUT1     = MY_IOFF_LUOUT1  + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT3 = 1 + NUM_BLK
          IOFFSET_INT_LUOUT1 = 1 + NUM_BLK
C
C         VEC1
          CALL RDVEC_BATCH_DRV4(LUOUT1,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT1,IOFFSET_INT_LUOUT1,
     &                          LUOUT1LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUOUT1'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C         VEC2
          CALL RDVEC_BATCH_DRV4(LUOUT3,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT3,IOFFSET_INT_LUOUT3,
     &                          LUOUT3LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC2 on LUOUT3'
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
C         VEC2 == VEC2 + VEC1 * FACTOR
C
          CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
C         write VEC2 on LUOUT3
C      
          IOFFSET_LUOUT3     = MY_IOFF_LUOUT3  + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT3 = 1 + NUM_BLK
C
CSK       WRITE(LUWRT,*) 'final VEC2 to write on LUOUT3'
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          CALL WTVEC_BATCH_DRV4(LUOUT3,VEC2,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_LUOUT3,IOFFSET_INT_LUOUT3,
     &                         LUOUT3LIST,NUM_ACTIVE_BATCH)
C
C         keep track of correct offset
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
          NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
        END DO
C
      END IF
C     ^ IOLSAC ?
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE P1_B_PAR_RL_LUCI3(VEC1,VEC2,SUBSPH,
     &                             LUIN1LIST,LUIN2LIST,LUOUTLIST,
     &                             LUOUT2LIST,
     &                             NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                             MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                             MY_IOFF_LUOUT,MY_IOFF_LUOUT2,
     &                             SCRRED,NVEC,IADD,
     &                             LUIN1,LUIN2,LUOUT,LUOUT2)
C
C     Written by  S. Knecht         - March 13 2008
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: NVEC = NVEC
C           IADD = IADD
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUOUTLIST(*)
      DIMENSION SCRRED(*), LUOUT2LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT, IOFFSET_LUOUT2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
      LOGICAL STORE_F
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_LUOUT2    = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
      IOFFSET_INT_LUOUT2= 0
      STORE_F = .FALSE.
      IF( NVEC + IADD - 1 - NROOT_INFO .gt. 0 ) STORE_F = .TRUE.
C
      REDSCRVAR = 0.0D0
      CALL DZERO(SCRRED,NVEC+IADD)
      CALL DZERO(SUBSPH,NVEC+IADD)
CSK   WRITE(LUWRT,*) ' NVEC + IADD - 1',  NVEC + IADD - 1
CSK   WRITE(LUWRT,*) ' LUIN1,LUIN2,LUOUT,LUOUT2', 
CSK  &                 LUIN1,LUIN2,LUOUT,LUOUT2
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IVEC - 1
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK 
C
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN2 in P1..._3 100',
CSK  &                  IOFFSET_IN_LUIN2
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN2 in P1..._3 100'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)

C
        DO 100 IVEC = 1, NROOT_INFO
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
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
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN1 in P1..._3 100',
CSK  &                   IOFFSET_IN_LUIN1
CSK       WRITE(LUWRT,*) 'This is my INT_OFFSET for LUIN1 in 
CSK  &                    P1..._3 100',IOFFSET_INT_IN1
CSK       WRITE(LUWRT,*) 'THIS IS MY LU1LIST inside P1_B_PAR_RL_3 100'
CSK       CALL IWRTMAMN(LUIN1LIST,1,IALL_LU1,1,IALL_LU1,LUWRT)
C
          CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                         LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 100'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          SCRRED(IVEC) = SCRRED(IVEC) -
     &                   DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  100   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
C     communicate SCRRED to get full OVERLAP matrix
C
      CALL REDVEC_REL(SCRRED,SUBSPH,NROOT_INFO,2,MPI_SUM,
     &                MPI_COMM_WORLD,-1)
C
C
C
CSK   WRITE(LUWRT,*) ' THIS IS MY SUBSPH in P1..._3'
CSK   CALL WRTMATMN(SUBSPH,1,NVEC+IADD-1,1,NVEC+IADD-1,LUWRT)
C
C     zero scratch offsets
C
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      REDSCRVAR         = 0.0D0
C
C
CSK   WRITE(LUWRT,*) ' THIS IS MY LUIN2LIST in P1..._3'
CSK   CALL IWRTMAMN(LUIN2LIST,1,IALL_LU3,1,IALL_LU3,LUWRT)
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IVEC - 1
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK  &                  IOFFSET_IN_LUIN2
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN2'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)

C
        DO 200 IVEC = 1, NROOT_INFO
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
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
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN1 in P1..._3 200',
CSK  &                   IOFFSET_IN_LUIN1
C
          CALL RDVEC_BATCH_DRV4(LUIN1,VEC1,LBATCH(ISBATCH),
     &                         IBATCH(1,I1BATCH(ISBATCH)),
     &                         IOFFSET_IN_LUIN1,IOFFSET_INT_IN1,
     &                         LUIN1LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 200'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          FACTOR = SUBSPH(IVEC)
C
C         VEC2 == VEC2 + VEC1 * FACTOR
C
          CALL DAXPY(LEBATCH(ISBATCH), FACTOR, VEC1, 1, VEC2, 1) 
C
  200   CONTINUE
C
CSK     WRITE(LUWRT,*) 'final VEC2 to write on LUOUT2 '
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
        IF( STORE_F )THEN
C
C         new offset for writing on LUOUT2 --> ILU5
C
          IOFFSET_LUOUT2     = MY_IOFF_LUOUT2 + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
          CALL WTVEC_BATCH_DRV4(LUOUT2,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                          LUOUT2LIST,NUM_ACTIVE_BATCH)
        ELSE
C
          REDSCRVAR = REDSCRVAR
     &              + DDOT( LEBATCH(ISBATCH), VEC2, 1, VEC2, 1)
C
C         new offset for writing on LUIN2
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
          IOFFSET_INT_IN2  = 1 + NUM_BLK
C
          CALL WTVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
        END IF
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      IF( STORE_F )THEN
C
C        zero scratch offsets
C
         NUM_BLK           = 0
         IOFFSET_SCRATCH   = 0
         CALL DZERO(SCRRED,NVEC+IADD)
         CALL DZERO(SUBSPH,NVEC+IADD)
C
         DO ISBATCH = 1, NBATCH
C
C          offset for batch ISBATCH w.r.t JOFF
C
           CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C          set new offset
C
C          position in file is at the end of vector IVEC - 1
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
           IOFFSET_INT_IN2  = 1 + NUM_BLK 
C
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2 in P1..._3 100',
CSK  &                     IOFFSET_IN_LUIN2
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK        WRITE(LUWRT,*) 'initial VEC2 on LUIN2 in P1..._3 100'
CSK        CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                   LUWRT)

C
           DO 300 IVEC = 1, (NVEC + IADD - 1 -NROOT_INFO)
C
             CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C            set new offset
C
C            position in file is at the end of vector IVEC - 1
C
             IOFFSET_LUOUT     = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                         ( IVEC - 1 )  * MY_VEC1_IOFF
             IOFFSET_INT_LUOUT = 1 + NUM_BLK   +
     &                         ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK          WRITE(LUWRT,*) 'This is my OFFSET for LUOUT in P1..._3 100',
CSK  &                       IOFFSET_IN_LUOUT
CSK          WRITE(LUWRT,*) 'This is my INT_OFFSET for LUOUT in 
CSK  &                       P1..._3 100',IOFFSET_INT_LUOUT
C
             CALL RDVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                             IBATCH(1,I1BATCH(ISBATCH)),
     &                             IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                             LUOUTLIST,NUM_ACTIVE_BATCH)
C
CSK          WRITE(LUWRT,*) 'initial VEC1 on LUOUT in P1..._3 100'
CSK          CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                     LUWRT)
C
             SCRRED(IVEC) = SCRRED(IVEC) -
     &                      DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
  300     CONTINUE
C
C         keep track of correct offset
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
          NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
        END DO
C
C       communicate SCRRED to get full OVERLAP matrix
C
        CALL REDVEC_REL(SCRRED,SUBSPH,NVEC+IADD-1-NROOT_INFO,2,MPI_SUM,
     &                  MPI_COMM_WORLD,-1)
C
C       zero scratch offsets
C
        NUM_BLK           = 0
        IOFFSET_SCRATCH   = 0
        REDSCRVAR         = 0.0D0
C
        DO ISBATCH = 1, NBATCH
C
C         offset for batch ISBATCH w.r.t JOFF
C
          CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C         set new offset
C
C         position in file is at the end of vector IVEC - 1
C
          IOFFSET_LUOUT2     = MY_IOFF_LUOUT2 + IOFFSET_SCRATCH
          IOFFSET_INT_LUOUT2 = 1 + NUM_BLK
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUOUT2',
CSK  &                    IOFFSET_IN_LUOUT2
C
          CALL RDVEC_BATCH_DRV4(LUOUT2,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT2,IOFFSET_INT_LUOUT2,
     &                          LUOUT2LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC2 on LUOUT2'
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
          DO 400 IVEC = 1, (NVEC+IADD-1-NROOT_INFO)
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
C           set new offset
C
C           position in file is at the end of vector IVEC - 1
C
            IOFFSET_LUOUT     = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                        ( IVEC - 1 )  * MY_VEC1_IOFF
            IOFFSET_INT_LUOUT = 1 + NUM_BLK   +
     &                        ( IVEC - 1 )  * MY_ACT_BLK1
C
CSK         WRITE(LUWRT,*) 'This is my OFFSET for LUOUT in P1..._3 100',
CSK  &                      IOFFSET_IN_LUOUT
CSK         WRITE(LUWRT,*) 'This is my INT_OFFSET for LUOUT in
CSK  &                      P1..._3 100',IOFFSET_INT_LUOUT
C
            CALL RDVEC_BATCH_DRV4(LUOUT,VEC1,LBATCH(ISBATCH),
     &                            IBATCH(1,I1BATCH(ISBATCH)),
     &                            IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                            LUOUTLIST,NUM_ACTIVE_BATCH)
CSK         WRITE(LUWRT,*) 'initial VEC1 on LUIN1 in P1..._3 200'
CSK         CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                    LUWRT)
C
            FACTOR = SUBSPH(IVEC)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH), FACTOR, VEC1, 1, VEC2, 1) 
C
  400     CONTINUE
C
CSK       WRITE(LUWRT,*) 'final VEC2 to write on LUIN2 '
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
          REDSCRVAR = REDSCRVAR
     &              + DDOT( LEBATCH(ISBATCH), VEC2, 1, VEC2, 1)
C
C         new offset for writing on LUIN2
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
          IOFFSET_INT_IN2  = 1 + NUM_BLK
C
          CALL WTVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
C         keep track of correct offset
          IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
          NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
        END DO
C
      END IF
C     ^ NVEC + IADD - 1 - NROOT > 0 ( STORE_F == .TRUE. )
C
      SCALEVEC = 0.0D0
C
C     communicate REDSCRVAR to get full scale factor
C       
      CALL REDVEC_REL(REDSCRVAR,SCALEVEC,1,2,MPI_SUM,MPI_COMM_WORLD,-1)
C
C     1.4 normalizing the new vector
C
C     zero scratch offsets
C
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
C
C
      FACTOR = 1.0D0 / SQRT( SCALEVEC )
csk   WRITE(LUWRT,*) 'THIS IS MY SCALING FACTOR',FACTOR
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
C       set new offset
C
C       position in file is at the end of vector IVEC - 1
C
        IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH
        IOFFSET_INT_IN2  = 1 + NUM_BLK
C
C
CSK     WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK  &                  IOFFSET_IN_LUIN2
C
        CALL RDVEC_BATCH_DRV4(LUIN2,VEC2,LBATCH(ISBATCH),
     &                        IBATCH(1,I1BATCH(ISBATCH)),
     &                        IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                        LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK     WRITE(LUWRT,*) 'initial VEC2 on LUIN2'
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)

C
        CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
C
C       set new offset
C
C       position in file is at the end of vector NVEC + IADD - 1 - NROOT
C
        IOFFSET_LUOUT  = MY_IOFF_LUOUT + IOFFSET_SCRATCH +
     &                 ( NVEC + IADD - 1 - NROOT_INFO ) * MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT = 1 + NUM_BLK +
     &                 ( NVEC + IADD - 1 - NROOT_INFO ) * MY_ACT_BLK1
C
csk     WRITE(LUWRT,*) 'This is my OFFSET for LUOUT, IOFFSET_INT_LUOUT',
csk  &                  IOFFSET_LUOUT, IOFFSET_INT_LUOUT
C
csk     WRITE(LUWRT,*) 'absolute final new vec on VEC2 to LUOUT'
csk     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
csk  &                LUWRT)
C
        IDEBUGPRNT = 0
C
        CALL WTVEC_BATCH_DRV4SP(LUOUT,VEC2,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_LUOUT,IOFFSET_INT_LUOUT,
     &                          LUOUTLIST,NUM_ACTIVE_BATCH,IDEBUGPRNT)
C
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE INPROD_B_PAR_RL_LUCI2(LUIN1,LUIN2,LUIN3,VEC1,VEC2,
     &                                 SUBSPH,NBATCH,
     &                                 LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                                 MY_IOFF_LUIN1,MY_IOFF_LUIN2,
     &                                 MY_IOFF_LUIN3,LUIN1LIST,
     &                                 LUIN2LIST,LUIN3LIST,IVEC,
     &                                 NVEC,IMUSTRED,ISTRED)
C
C     Written by  S. Knecht         - March 14 2008
C
C**********************************************************************
C
C     calculating dot product between two vectors on file LUIN1 resp.
C     LUIN2
C
C     NOTE: IVEC = IVEC
C           NVEC = NVEC
C
C     active blocks on the MPI-files are flagged by a nonzero length
C
C     Last revision:     S. Knecht       - June  2007
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUIN1LIST(*), LUIN2LIST(*), LUIN3LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN1
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN1, IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN3
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN1  = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_IN_LUIN3  = 0
      IOFFSET_INT_IN1   = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_IN3   = 0
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
     &                     ( JVEC_SF ) * MY_VEC1_IOFF
        IOFFSET_INT_IN1  = 1 + NUM_BLK  +
     &                     ( JVEC_SF ) * MY_ACT_BLK1
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
        DO 100 JVEC = 1, NROOT_INFO
C
C          set new offset and zero read-in vector
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                        ( JVEC - 1 )  * MY_VEC1_IOFF
           IOFFSET_INT_IN2  = 1 + NUM_BLK  + 
     &                        ( JVEC - 1 ) * MY_ACT_BLK1
C
CSK        WRITE(LUWRT,*) 'This is my OFFSET for LUIN2',
CSK     &                  IOFFSET_IN_LUIN2
C
C
C          read in batch ISBATCH from LUIN1 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                           LUIN2LIST,NUM_ACTIVE_BATCH)
C
csk        WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
csk        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,
csk  &                   LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
C          SUBSPH(IJ) == VEC1 * VEC2
C
           IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
csk        WRITE(LUWRT,*) ' IJ in loop 1', IJ
C
           SUBSPH(IJ) = SUBSPH(IJ) + 
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
C          keep track of memory offset and the 'reduction' counter
C
           IF( ISBATCH .eq. 1 ) THEN
             IMUSTRED = IMUSTRED + 1
             IF( IVEC .eq. 1 .and. JVEC .eq. 1 ) ISTRED = IJ
           END IF
C
C
  100   CONTINUE
C
        JJVEC = 0
C
        DO 200 JVEC = NROOT_INFO+1 , NVEC+IVEC
C
C          set new offset and zero read-in vector
C
           JJVEC = JJVEC + 1
C
           CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
           IOFFSET_IN_LUIN3 = MY_IOFF_LUIN3 + IOFFSET_SCRATCH +
     &                        ( JJVEC - 1 )  * MY_VEC1_IOFF
           IOFFSET_INT_IN3  = 1 + NUM_BLK  + 
     &                        ( JJVEC - 1 )  * MY_ACT_BLK1
C
csk        WRITE(LUWRT,*) 'This is my OFFSET for LUIN3',
csk  &                     IOFFSET_IN_LUIN3, IOFFSET_INT_IN3
C
C
C          read in batch ISBATCH from LUIN3 to VEC1
C
           CALL RDVEC_BATCH_DRV4(LUIN3,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN3,IOFFSET_INT_IN3,
     &                           LUIN3LIST,NUM_ACTIVE_BATCH)
C
csk        WRITE(LUWRT,*) 'initial VEC1 on LUIN3'
csk        CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,
csk  &                   LEBATCH(ISBATCH),LUWRT)
C
C
C          IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
C          SUBSPH(IJ) == VEC1 * VEC2
C
           IJ = (IVEC+NVEC)*(IVEC+NVEC-1)/2 + JVEC
csk        WRITE(LUWRT,*) ' IJ in loop 2', IJ
C
           SUBSPH(IJ) = SUBSPH(IJ) + 
     &                  DDOT(LEBATCH(ISBATCH),VEC1,1,VEC2,1)
C
C          keep track of the 'reduction' counter
C
           IF( ISBATCH .eq. 1 ) THEN
             IMUSTRED = IMUSTRED + 1
           END IF
C
C
  200   CONTINUE
C
C       keep track of correct offset
        IOFFSET_SCRATCH = IOFFSET_SCRATCH + LEBATCH(ISBATCH)
        NUM_BLK         = NUM_BLK + NUM_ACTIVE_BATCH
C
      END DO
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE P3_B_PAR_RL_LUCI1(VEC1,VEC2,SUBSPH,
     &                             LUINLIST,LUIN2LIST,LUOUTLIST,
     &                             NBATCH,LBATCH,LEBATCH,I1BATCH,IBATCH,
     &                             MY_IOFF_LUIN,MY_IOFF_LUIN2,
     &                             MY_IOFF_LUOUT,NVEC,NVEC2,IROOT,
     &                             LUIN,LUIN2,LUOUT)
C
C     Written by  S. Knecht         - March 14 2008
C
C**********************************************************************
C
C     calculating scaled vecsum between two vectors on file LUIN resp.
C     LUIN2; saving on LUOUT
C
C     NOTE: IROOT = IROOT
C           NVEC  = NVEC
C           NVEC2 = NROOT
C
C     active blocks on the MPI-files are flagged by a nonzero list entry
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT DOUBLE PRECISION ( A-H,O-Z)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION VEC1(*), VEC2(*), SUBSPH(*)
      DIMENSION LBATCH(*), LEBATCH(*), I1BATCH(*), IBATCH(8,*)
      DIMENSION LUINLIST(*), LUOUTLIST(*), LUIN2LIST(*)
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) MY_IOFF_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_IN_LUIN2
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_LUOUT
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET_SCRATCH
      INTEGER NUM_BLK
C
C     initialize scratch offsets
      NUM_BLK           = 0
      IOFFSET_SCRATCH   = 0
      IOFFSET_IN_LUIN   = 0
      IOFFSET_IN_LUIN2  = 0
      IOFFSET_LUOUT     = 0
      IOFFSET_INT_IN    = 0
      IOFFSET_INT_IN2   = 0
      IOFFSET_INT_LUOUT = 0
C
C
      DO ISBATCH = 1, NBATCH
C
C       offset for batch ISBATCH w.r.t JOFF
C
        CALL DZERO(VEC2,LEBATCH(ISBATCH))
C
        DO 100 IVEC = 1, NVEC2
C
          IJ = (( IROOT - 1 ) * NVEC ) + 1 + ( IVEC - 1)
C
          FACTOR = SUBSPH( IJ )
C
C         set new offset
C
C         position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN  = MY_IOFF_LUIN + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )   * MY_VEC1_IOFF
          IOFFSET_INT_IN   = 1 + NUM_BLK  +
     &                     ( IVEC - 1 )   * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK  &                    IOFFSET_IN_LUIN
C
          IF( IVEC .eq. 1 ) THEN
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC2,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC2 on LUIN'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                    LUWRT)
CSK         WRITE(LUWRT,*) 'scaling factor for this vector',FACTOR
C
            CALL DSCAL(LEBATCH(ISBATCH),FACTOR,VEC2,1)
CSK         WRITE(LUWRT,*) ' VEC2 after first scaling '
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                    LUWRT)
C
          ELSE
C
            CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
            CALL RDVEC_BATCH_DRV4(LUIN,VEC1,LBATCH(ISBATCH),
     &                           IBATCH(1,I1BATCH(ISBATCH)),
     &                           IOFFSET_IN_LUIN,IOFFSET_INT_IN,
     &                           LUINLIST,NUM_ACTIVE_BATCH)
C
CSK         WRITE(LUWRT,*) 'initial VEC1 on LUIN'
CSK         CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                 LUWRT)
C
C           VEC2 == VEC2 + VEC1 * FACTOR
C
            CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
CSK         WRITE(LUWRT,*) 'final VEC2 after DAXPY in P3...1'
CSK         CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                    LUWRT)
C
          END IF
C         ^ IVEC == 1 ?
C
  100   CONTINUE
C
        DO 200 IVEC = 1, (NVEC - NVEC2)
C
          IJ = (( IROOT - 1 ) * NVEC ) + 1 + ( IVEC - 1) + NVEC2
C
          FACTOR = SUBSPH( IJ )
C
C         set new offset
C
C         position in file is at the end of vector IVEC - 1
C
          IOFFSET_IN_LUIN2 = MY_IOFF_LUIN2 + IOFFSET_SCRATCH +
     &                     ( IVEC - 1 )    * MY_VEC1_IOFF
          IOFFSET_INT_IN2  = 1 + NUM_BLK   +
     &                     ( IVEC - 1 )    * MY_ACT_BLK1
C
CSK       WRITE(LUWRT,*) 'This is my OFFSET for LUIN',
CSK     &                 IOFFSET_IN_LUIN
C
          CALL DZERO(VEC1,LEBATCH(ISBATCH))
C
          CALL RDVEC_BATCH_DRV4(LUIN2,VEC1,LBATCH(ISBATCH),
     &                          IBATCH(1,I1BATCH(ISBATCH)),
     &                          IOFFSET_IN_LUIN2,IOFFSET_INT_IN2,
     &                          LUIN2LIST,NUM_ACTIVE_BATCH)
C
CSK       WRITE(LUWRT,*) 'initial VEC1 on LUIN2'
CSK       CALL WRTMATMN(VEC1,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
CSK       WRITE(LUWRT,*) 'scaling factor for this vector',FACTOR
C
C         VEC2 == VEC2 + VEC1 * FACTOR
C
          CALL DAXPY(LEBATCH(ISBATCH),FACTOR,VEC1,1,VEC2,1)
C
CSK       WRITE(LUWRT,*) 'final VEC2 after 2nd DAXPY in P3...1'
CSK       CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                  LUWRT)
C
  200   CONTINUE
C
CSK     WRITE(LUWRT,*) 'final VEC2 to write on position',IROOT -1
CSK     CALL WRTMATMN(VEC2,1,LEBATCH(ISBATCH),1,LEBATCH(ISBATCH),
CSK  &                LUWRT)
C
C       new offset for writing on LUOUT
C
        IOFFSET_LUOUT      =  MY_IOFF_LUOUT + IOFFSET_SCRATCH + 
     &                        ( IROOT - 1 ) * MY_VEC1_IOFF
C
        IOFFSET_INT_LUOUT  = 1 + NUM_BLK    +
     &                        ( IROOT - 1 ) * MY_ACT_BLK1
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
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE UPDATE_LUC_LIST2(ISCLFAC_GROUP,LUCLIST,RCCTOS,CB,
     &                            NPARBLOCK,IBLOCKL,IGROUPLIST,
     &                            IPROCLIST,IRILP,BLOCKTIME)
C
C     make an update of of grouplist for c-vector file based on 
C     different list gathered from MPI_COMM_WORLD
C
C
C     Written by  S. Knecht         - March 08 2008 
C
C     OUTPUT: ISCLFAC_GROUP and updated file ILUC
C
C**********************************************************************
#include "implicit.h"
      INTEGER*8 KSCALLOC2, KSCALLOC3
!               for addressing of WORK
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
#include "wrkspc.inc"
      DIMENSION ISCLFAC_GROUP(*), LUCLIST(*) 
      DIMENSION CB(*)
      DIMENSION NPARBLOCK(*), IBLOCKL(*)
      DIMENSION IGROUPLIST(*), IPROCLIST(*)
      CHARACTER*12 WALLTID3, SECTID
      INTEGER RCCTOS(*)
C     some scratch
      INTEGER NZERO
C
      NZERO = 0
C
C     set mark for local memory
C
      IDUM = 0
      CALL MEMMAN(KDUM,  IDUM,    'MARK  ',IDUM,'UPLIST')
C
      CALL MEMMAN(KSCALLOC2,NUM_BLOCKS2,'ADDL  ',1,'ICLLC2')
      CALL MEMMAN(KSCALLOC3,NUM_BLOCKS2,'ADDL  ',1,'ICLLC3')
C
C     fill complete local iscalfac arrays with zero's
      CALL IZERO(WORK(KSCALLOC2), NUM_BLOCKS2)
      CALL IZERO(WORK(KSCALLOC3), NUM_BLOCKS2)
      CALL IZERO(ISCLFAC_GROUP  , NUM_BLOCKS2)
C
#ifdef LUCI_DEBUG
      WRITE(luwrt,*) '  start of subroutine UPDATE_LUC_LIST2 speaking'
      WRITE(luwrt,*) 'LUCLIST:'
      CALL IWRTMAMN(LUCLIST,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
      WRITE(luwrt,*) 'ISCLFAC_GROUP:'
      CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
      WRITE(luwrt,*) 'WORK(KSCALLOC2):'
      CALL IWRTMAMN(WORK(KSCALLOC2),1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
      WRITE(luwrt,*) 'WORK(KSCALLOC3):'
      CALL IWRTMAMN(WORK(KSCALLOC3),1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
#endif
C
      starttimer = MPI_WTIME()
C
C     "mpi_allsum" local LUCLIST which then on all
C     nodes will contain the number of non-zero C-blocks in
C     the complete CI-vector
C
      CALL REDVEC_REL(LUCLIST,WORK(KSCALLOC2),NUM_BLOCKS2,1,
     &                MPI_SUM,MPI_COMM_WORLD,-1)
C
C     find all c-blocks connecting to all sigma-blocks on each cpu
C
      CALL ICOPY(NUM_BLOCKS2,RCCTOS,1,WORK(KSCALLOC3),1)
C
C     case 1: number of CPUs in new group not equal to total number
C
C     case 2: number of CPUs in new group equal to total number
C
C
#ifdef LUCI_DEBUG
      write(luwrt,*) 'NEWCOMM_PROC,LUCI_NMPROC',NEWCOMM_PROC,LUCI_NMPROC
      write(luwrt,*) 'ILUC',ILUC
      write(luwrt,*) 'MYNEW_COMM',ILUC
      write(luwrt,*) 'ICOMM',ICOMM
      write(luwrt,*) 'ICOMM_id, ICOMM_SIZE',ICOMM_id,ICOMM_SIZE
#endif
*
      IF( NEWCOMM_PROC .ne. LUCI_NMPROC ) THEN
*
        CALL REDVEC_REL(WORK(KSCALLOC3),ISCLFAC_GROUP,NUM_BLOCKS2,1,
     &                  MPI_SUM,MYNEW_COMM,0)
*
*       all local node-masters call this routine!
*
        IF( MYNEW_ID .eq. 0 ) THEN
           CALL COPVCD_PAR_BDRIV5_REL(ILUC,ILUC,CB,NPARBLOCK,
     &                                WORK(KSCALLOC2),ISCLFAC_GROUP,
     &                                IBLOCKL,NUM_BLOCKS,ICOMM,
     &                                IGROUPLIST,IPROCLIST,IRILP)
C               COPVCD_PAR_BDRIV5_REL(LUIN,LUOUT,SEGMNT,IBLOCKD,
C     &                               ISCALFAC,ISCALFAC_GROUP,
C     &                               IBLOCKL,NBLOCK,JCOMM,
C     &                               IGROUPLIST,IPROCLIST,IRILP)
C

        END IF
           CALL MPI_BCAST(ISCLFAC_GROUP,NUM_BLOCKS2,
     &                    my_MPI_INTEGER,0,MYNEW_COMM,IERR)
*
      ELSE
*
C
         CALL UPDATE_GEN_LIST(WORK(KSCALLOC3),WORK(KSCALLOC2),
     &                        NUM_BLOCKS2)
C
C        to be consistent with output of case 1
C
         CALL IZERO(ISCLFAC_GROUP,NUM_BLOCKS2)
         CALL ICOPY(NUM_BLOCKS2,WORK(KSCALLOC3),1,ISCLFAC_GROUP,1)
C
      END IF
C     ^ NEWCOMM_PROC == LUCI_NMPROC ?
C
#ifdef LUCI_DEBUG
      WRITE(luwrt,*) '  subroutine UPDATE_LUC_LIST2 speaking'
      WRITE(luwrt,*) 'LUCLIST:'
      CALL IWRTMAMN(LUCLIST,1,NUM_BLOCKS2,1,NUM_BLOCKS2,luwrt)
      WRITE(luwrt,*) 'ISCLFAC_GROUP:'
      CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,luwrt)
#endif
C
C     final timing for block distribution
      blocktime = blocktime + MPI_WTIME() - starttimer
C
C     flush local memory
C
      IDUM = 0
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',2,'UPLIST')
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE BLOCK_DISTR_DRV(NBLOCK,IBLOCKL,NBLOCKD,IBLOCKS_FNODE,
     &                           SCALFAC,NVAR,IPROCLIST)
#include "implicit.h"
      INTEGER*8 KICCTOS, KCWEIGHT, KCWEIGHTF, KBLCKWT, KIBTOTW
!               for addressing of WORK
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
#include "parluci.h"
#include "wrkspc.inc"
! from cands.inc: icsm
#include "cands.inc"
      INTEGER*8 IABSOLUTE_WEIGHT
      IABSOLUTE_WEIGHT = 0
C
      IDUM  = 0
      CALL MEMMAN(KDUM,IDUM,'MARK  ',IDUM,'BLKDRV')
C
C     allocate local scratch arrays
C
      CALL MEMMAN(KICCTOS,      NBLOCK**2,'ADDL  ',1,'ICCTOS')
      CALL MEMMAN(KCWEIGHT,        NBLOCK,'ADDL  ',1,'ICWHT ')
      CALL MEMMAN(KCWEIGHTF,       NBLOCK,'ADDL  ',1,'ICWHTF')
      CALL MEMMAN(KBLCKWT,  2*LUCI_NMPROC,'ADDL  ',1,'IBLCKW')
      CALL MEMMAN(KIBTOTW,         NBLOCK,'ADDL  ',2,'IBTOTW')

      CALL IZERO(WORK(KBLCKWT),2*LUCI_NMPROC)
      CALL IZERO(WORK(KCWEIGHT),NBLOCK)
      CALL IZERO(WORK(KCWEIGHTF),NBLOCK)
      CALL IZERO(WORK(KICCTOS),NBLOCK**2)
C
      CALL FIND_IMAT_SC(IBLOCKL,SCALFAC,WORK(KICCTOS),WORK(KCWEIGHT),
     &                  WORK(KIBTOTW),WORK(KCWEIGHTF),NBLOCK,
     &                  IABSOLUTE_WEIGHT)
C
      IF(IDISTROUTE.EQ.1) THEN
        CALL DISTBLKND_1(NBLOCK,WORK(KCWEIGHTF),NBLOCKD,WORK(KIBTOTW),
     &                   WORK(KBLCKWT),NVAR,WORK(KICCTOS),IBLOCKL,
     &                   IPROCLIST,WORK(KCWEIGHT),IABSOLUTE_WEIGHT)
      ELSE
        CALL DISTBLKND_2(NBLOCK,WORK(KCWEIGHTF),NBLOCKD,IBLOCKL,icsm)
      END IF
C
C     find all c-blocks connecting to a given sigma-block,
C     information will be stored in CBLOCKS_FNODE
C
      CALL FIND_CCTOS(IBLOCKS_FNODE,NBLOCKD,WORK(KICCTOS),NBLOCK)
C
C     eliminate local memory
      IDUM = 0
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',IDUM,'BLKDRV')
C
      END 
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS_PAR(LUIN,LBLK,BLK_A,SEGMNT,
     &                                  NBLOCK,IBLOCKD)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*  parallel version
*
      IMPLICIT REAL*8(A-H,O-Z)
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION IBLOCKD(NBLOCK)
*
      CALL REWINE(LUIN,LBLK)
*
      NBLK_A = 0
*. Loop over blocks
      DO 500 IBLK = 1, NBLOCK
*        
        IF(LUCI_MYPROC.NE.IBLOCKD(IBLK))THEN
          BLK_A(IBLK) = 0.0D0
          GOTO 300
        ELSE
          CALL IFRMDS(LBL,1,-1,LUIN)
          IF( LBL .GE. 0 ) THEN
            IF(LBLK .GE.0 ) THEN
              KBLK = LBL
            ELSE
              KBLK = -1
            END IF
            NO_ZEROING = 1
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
            IF(IMZERO.EQ.0) THEN
             NBLK_A = NBLK_A + 1
             BLK_A(IBLK) = 1.0D0
            ELSE
             BLK_A(IBLK) = 0.0D0
            END IF
          END IF
        END IF
 300  CONTINUE
*
 500  CONTINUE
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,*)'myproc',LUCI_MYPROC
        WRITE(6,'(A,I8,I8)')
     &  ' FIND_A.... Number of total and active Blocks',NBLOCK,NBLK_A
      END IF
*
      END
***********************************************************************
      SUBROUTINE FNDMND_PAR(LU,LBLK,SEGMNT,NSUBMX,NSUB,ISCR,
     &                     SCR,ISCAT,SUBVAL,IBLOCKL,IBLOCKD,
     &                     NBLOCK,NTESTG)
*
* FIND NSUB LOWEST ELEMENTS OF VECTOR RESIDING ON FILE
* LU. ENSURE THAT NO DEGENERENCIES ARE SPLIT
*
*
* INPUT
*=======
* LU :    FILE WHERE VECTOR OF INTEREST IS RESIDING, REWOUND
* LBLK :  DEFINES FILE STRUCTURE ON FILE LU
* NSUBMX: LARGEST ALLOWED NUMBER OF SORTED ELEMENTS
*
* OUTPUT
*=======
* NSUB : ACTUAL NUMBER OF ELEMENTS OBTAINED. CAN BE SMALLER
*        THAN NSUBMX IF THE LAST ELEMENT BELONGS TO A DEGENERATE
*        SET
*ISCAT:  SCATTERING ARRAY, ISCAT(I) GIVES FULL ADRESS OF SORTED
*        ELEMENT I
*SUBVAL: VALUE OF SORTED ELEMENTS

      IMPLICIT REAL*8           ( A-H,O-Z)
      INTEGER*8 KGATHERA, KGATHERB, KGATHERC, KGATHERD, KGATHERE
!               for addressing of WORK
#include "maxorb.h"
#include "infpar.h"
#include "mpif.h"
      DIMENSION my_STATUS(MPI_STATUS_SIZE)
#include "parluci.h"
#include "wrkspc.inc"
      
      DIMENSION SEGMNT(*), ISCAT(*),SUBVAL(*),SCR(*),ISCR(*)
      DIMENSION IBLOCKL(NBLOCK), IBLOCKD(NBLOCK)
      INTEGER(KIND=MPI_OFFSET_KIND) IOFF_SCR
C
      NTESTL = 0000
      NTEST = MAX(NTESTG,NTESTL)
      NTEST = 000
C     offset initialization
      IOFF_SCR = 0
      IOFF_SCR = IOFF_SCR + MY_DIA_OFF
C
      IBASE = 1
      LSUB = 0
*     loop over blocks
      DO 1000 II = 1, NBLOCK    
*
        IF( IBLOCKD(II) .eq. LUCI_MYPROC )THEN
          LBL = IBLOCKL(II) 
        ELSE
*         useful to set all other blocks to 0?
          LBL = 0
        ENDIF
*
        IF(NTEST.GE.10) THEN
          WRITE(LUWRT,*) ' Info about block ',II
          WRITE(LUWRT,*) ' Number of elements ',LBL
        END IF
        IF(LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
          IF( IBLOCKD(II) .eq. LUCI_MYPROC )THEN 
             CALL MPI_FILE_READ_AT(IDIA,IOFF_SCR,SEGMNT,LBL,
     &                             MPI_REAL8,my_STATUS,IERR)
          ENDIF
          IF(NTEST.GE.100) THEN
            WRITE(LUWRT,*) ' Elements read in '
            CALL WRTMATMN(SEGMNT,1,LBL,1,LBL,LUWRT)
          END IF
          IF(LBL .GE. 0 ) THEN
*. LOWEST ELEMENTS IN SEGMNT  ( ADD TO PREVIOUS LIST )
            MSUBMX = MIN(NSUBMX,LBL)
            IF(LBL.GE.1) THEN
              CALL SORLOW(SEGMNT,SCR(1+LSUB),ISCR(1+LSUB),LBL,
     &                    MSUBMX,MSUB,NTEST)
            ELSE
              MSUB = 0
            END IF
            DO 10 I = 1, MSUB
   10         ISCR(LSUB+I) = ISCR(LSUB+I) + IBASE - 1
* SORT COMBINED LIST
            MSUBMX = MIN(NSUBMX,LSUB+MSUB)
            IF(MSUBMX.GT.0) THEN
              CALL SORLOW(SCR,SUBVAL,ISCAT,LSUB+MSUB,MSUBMX,LLSUB,
     &                    NTEST)
            ELSE
              LLSUB = 0
            END IF
            LSUB = LLSUB
            DO 20 I = 1, LSUB
              ISCR(I+2*NSUBMX) = ISCR(ISCAT(I))
   20       CONTINUE
*
            CALL ICOPVE(ISCR(1+2*NSUBMX),ISCR(1),LSUB)
            CALL DCOPY(LSUB,SUBVAL,1,SCR,1)

            IF(NTEST .GE. 20 ) THEN
              WRITE(LUWRT,*)' Lowest elements and their original place'
              WRITE(LUWRT,*)' Number of elements obtained ', LSUB
              CALL WRTMATMN(SUBVAL,1,LSUB,1,LSUB,LUWRT)
              CALL IWRTMAMN(ISCR,1,LSUB,1,LSUB,LUWRT)
            END IF
          END IF
*
        END IF
        IOFF_SCR = IOFF_SCR + LBL
C       set to lbl to true value
        LBL = IBLOCKL(II)
        IBASE = IBASE + LBL
C
 1000 CONTINUE
*
      NTEST = 00
      NSUB = LSUB
      CALL ICOPVE(ISCR,ISCAT,NSUB)
      IF(NTEST .GE. 20) THEN
        WRITE(LUWRT,*) ' Lowest elements and their original place '
        WRITE(LUWRT,*) ' Number of elements obtained ', NSUB
        CALL WRTMATMN(SUBVAL,1,NSUB,1,NSUB,LUWRT)
        CALL IWRTMAMN(ISCAT,1,NSUB,1,NSUB,LUWRT)
      END IF
*
      IDUM = 0
      CALL MEMMAN(KDUM,IDUM,'MARK  ',IDUM,'GATHER')
*
      CALL MEMMAN(KGATHERA,LUCI_NMPROC*NSUBMX,'ADDL  ',2,'PARRA1')
      CALL MEMMAN(KGATHERB,LUCI_NMPROC*NSUBMX,'ADDL  ',2,'PARRA2')
      CALL MEMMAN(KGATHERC,LUCI_NMPROC*NSUBMX,'ADDL  ',1,'PARIA1')
      CALL MEMMAN(KGATHERD,LUCI_NMPROC*NSUBMX,'ADDL  ',1,'PARIA2')
      CALL MEMMAN(KGATHERE,NSUBMX,'ADDL  ',1,'PARIA3')
*. We gather all lowest values from each node 
*. and build up a combined list of those
      CALL GATHER_LOW_PAR(NSUB,NSUBMX,SUBVAL,ISCAT,
     &                    WORK(KGATHERA),WORK(KGATHERB),
     &                    WORK(KGATHERC),WORK(KGATHERD),
     &                    WORK(KGATHERE),NTESTG)
*     update SCR1 and ISCR1
      CALL DCOPY(NSUBMX,SUBVAL,1,SCR,1)
      CALL ICOPVE(ISCAT,ISCR,NSUBMX)
     
      IF(NTEST.GE.20)THEN
        WRITE(LUWRT,*)'after search: SUBVAL and ISCAT'
        CALL WRTMATMN(SUBVAL,1,NSUBMX,1,NSUBMX,LUWRT)
        CALL IWRTMAMN(ISCAT,1,NSUBMX,1,NSUBMX,LUWRT)
      END IF
CSK      NTEST = 0
           
      IF(NSUB.NE.NSUBMX.AND.NTEST.GE.20)THEN
        WRITE(LUWRT,*)'Warning! NSUB is lower than NSUBMX'
        WRITE(LUWRT,*)'NSUB is set to be equal to NSUBMX'
        NSUB = NSUBMX
      END IF
*
*. Eliminate local memory
      IDUM = 0
      CALL MEMMAN(KDUM ,IDUM,'FLUSM ',IDUM,'GATHER')
*
      END
*****************************************************
#else   /* VAR_MPI */
* dummy routine for normal compilation
       SUBROUTINE DO_NOTHING2(NOTHING)
       IMPLICIT REAL*8(A-H,O-Z)
*        write(6,*) 'let us do nothing',nothing  
       END 
#endif    /* VAR_MPI */
