C
C...   Copyright (c) 2004 by the authors of Dirac (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   Dirac, a relativistic ab initio electronic structure program,
C...   release DIRAC04.0 (2004),
C...   written by H. J. Aa. Jensen, T. Saue, and L. Visscher
C...   with contributions from V. Bakken, E. Eliav, T. Enevoldsen, T. Fleig,
C...   O. Fossgaard, T. Helgaker, J. Henriksson, J. K. Laerdahl, C. V. Larsen,
C...   P. Norman, J. Olsen, M. Pernpointner, J. K. Pedersen, K. Ruud,
C...   P. Salek, J. N. P. van Stralen, J. Thyssen, O. Visser, and T. Winther
C...   (http://dirac.chem.sdu.dk).
C...
C...   For a suitable BibTEX entry, see:
C...      http://dirac.chem.sdu.dk/doc/reference.shtml
C...
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dirac,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dirac-admin@dirac.chem.sdu.dk
C...
C...   For information on how to get a licence see:
C...      http://dirac.chem.sdu.dk
C
#if defined (VAR_MPI2)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
      WRITE(LUWRT,'(A19,7X,I3,3X,1E18.13,3X,1F19.10)')
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
#include "wrkspc.inc"
      DIMENSION ISCLFAC_GROUP(*), LUCLIST(*) 
      DIMENSION CB(*)
      DIMENSION NPARBLOCK(*), IBLOCKL(*)
      DIMENSION IGROUPLIST(*), IPROCLIST(*)
      CHARACTER*12 WALLTID3, SECTID
      DIMENSION RCCTOS(*)
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
csk   WRITE(LUCIWT,*) '  start of subroutine UPDATE_LUC_LIST2 speaking'
csk   WRITE(LUCIWT,*) 'LUCLIST:'
csk   CALL IWRTMAMN(LUCLIST,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
csk   WRITE(LUCIWT,*) 'ISCLFAC_GROUP:'
csk   CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
csk   WRITE(LUCIWT,*) 'WORK(KSCALLOC2):'
csk   CALL IWRTMAMN(WORK(KSCALLOC2),1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
csk   WRITE(LUCIWT,*) 'WORK(KSCALLOC3):'
csk   CALL IWRTMAMN(WORK(KSCALLOC3),1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
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
csk   WRITE(LUCIWT,*) '  subroutine UPDATE_LUC_LIST2 speaking'
csk   WRITE(LUCIWT,*) 'LUCLIST:'
csk   CALL IWRTMAMN(LUCLIST,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
csk   WRITE(LUCIWT,*) 'ISCLFAC_GROUP:'
csk   CALL IWRTMAMN(ISCLFAC_GROUP,1,NUM_BLOCKS2,1,NUM_BLOCKS2,LUCIWT)
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
      SUBROUTINE CALC_OFF_MPI_FILE2(FILENAME,IGROUPLIST,IBLOCKL,
     &                              IBLOCKD,IVEC_IN)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
      DIMENSION IGROUPLIST(*), IBLOCKL(*), IBLOCKD(*)
      CHARACTER*6 FILENAME
      INTEGER(KIND=MPI_OFFSET_KIND) ND_VEC_IOFF, ND_VEC_IOFF2
      INTEGER(KIND=MPI_OFFSET_KIND) MULT1, MULT2, JVEC
      INTEGER(KIND=MPI_OFFSET_KIND) I_VEC_LEN_SCR
C
      I_VEC_LEN_SCR = 0
      JVEC  = IVEC_IN
      MULT1 = 0
      MULT2 = 0
C
      IF( NPTEST_VAR .ge. 10 ) 
     & WRITE(LUWRT,*) 'Calculating offset for file', FILENAME
C
C     count active blocks for each cpu, calculate individual offset
C
      ND_VEC_IOFF = 0
      ND_VEC_IOFF2 = 0
      NBLK_ACT = 0
C
C     scratch file for c-vector
C
      IF( FILENAME .eq. 'ILUC  ' ) THEN
        DO JBLK = 1, NUM_BLOCKS
          IF( IBLOCKL(JBLK) .gt. 0 ) THEN
            NBLK_ACT = NBLK_ACT + 1
          END IF
        END DO
C
        MY_ACT_BLK_ALL = NUM_BLOCKS
C
        GOTO 101
C
      END IF
C
      NBLK_ACT = 0
C
C
C     we have group files, so only cpus within a group are relevant
C
      MY_VEC1_IOFF = 0
      MY_VEC2_IOFF = 0
C
      DO IPROC = 1, NEWCOMM_PROC
C     
        JPROC = IGROUPLIST(IPROC)
C
        DO JBLK = 1, NUM_BLOCKS
          IF(IBLOCKD(JBLK) .eq. JPROC ) THEN
C?            WRITE(LUWRT,*) 'assigning a block to JPROC',JBLK,JPROC
            ND_VEC_IOFF = ND_VEC_IOFF + IBLOCKL(JBLK)
            IF( LUCI_MYPROC .eq. JPROC ) THEN
              MY_VEC1_IOFF = MY_VEC1_IOFF + IBLOCKL(JBLK) 
              MY_VEC2_IOFF = MY_VEC2_IOFF + IBLOCKL(JBLK) 
              NBLK_ACT = NBLK_ACT + 1
            END IF
          END IF
        END DO
C
        IF( LUCI_MYPROC .eq. JPROC ) THEN
          MULT1 = 1
C
C         MY_VEC1_IOFF = length of one vector
C
          MY_VEC1_IOFF = MY_VEC1_IOFF
C
C         MY_VEC1_IOFF = length of one vector x complex multiplier ( x 1 = real!)
C
          I_VEC_LEN_SCR = MY_VEC2_IOFF
          MY_VEC2_IOFF = I_VEC_LEN_SCR * MULT1
C
          MY_ACT_BLK1 = NBLK_ACT
          MY_ACT_BLK2 = NBLK_ACT * MULT1

C
          IF( FILENAME .eq. 'LUDIA ' ) THEN
C
C           initialize MY_DIA_OFF
C
            MY_DIA_OFF = 0
C
            MY_DIA_OFF = ND_VEC_IOFF2 * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU1  ') THEN
C
C           initialize MY_LU1_OFF
C
            MY_LU1_OFF = 0
C
            MY_LU1_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU2  ') THEN
C
C           initialize MY_LU2_OFF
C
            MY_LU2_OFF = 0
C
            MY_LU2_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU3  ') THEN
C
C           initialize MY_LU3_OFF
C
            MY_LU3_OFF = 0
C
            MY_LU3_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU4  ') THEN
C
C           initialize MY_LU4_OFF
C
            MY_LU4_OFF = 0
C
            MY_LU4_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU5  ') THEN
C
C           initialize MY_LU5_OFF
C
            MY_LU5_OFF = 0
C
            MY_LU5_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU6  ') THEN
C
C           initialize MY_LU6_OFF
C
            MY_LU6_OFF = 0
C
            MY_LU6_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          ELSE IF( FILENAME .eq. 'ILU7  ') THEN
C
C           initialize MY_LU5_OFF
C
            MY_LU7_OFF = 0
C
            MY_LU7_OFF = ND_VEC_IOFF2 * JVEC * MULT1
C
          END IF
C         ^ filenames
        END IF
C
C       reset ND_VEC_IOFF2 and ND_VEC_IOFF values
C
        ND_VEC_IOFF2 = ND_VEC_IOFF
C
        NBLK_ACT = 0
C
      END DO
C
      IF( NPTESTVAR .ge. 10 ) THEN
        IF( FILENAME .eq. 'LUDIA ' ) THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'LUDIA offset :', MY_DIA_OFF
        ELSE IF( FILENAME .eq. 'ILU1  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU1  offset :', MY_LU1_OFF
        ELSE IF( FILENAME .eq. 'ILU2  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU2  offset :', MY_LU2_OFF
        ELSE IF( FILENAME .eq. 'ILU3  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU3  offset :', MY_LU3_OFF
        ELSE IF( FILENAME .eq. 'ILU4  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU4  offset :', MY_LU4_OFF
        ELSE IF( FILENAME .eq. 'ILU5  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU5  offset :', MY_LU5_OFF
        ELSE IF( FILENAME .eq. 'ILU6  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU6  offset :', MY_LU6_OFF
        ELSE IF( FILENAME .eq. 'ILU7  ') THEN
          WRITE(LUWRT,'(2X,A,1X,I18)') 'ILU7  offset :', MY_LU7_OFF
        END IF
        WRITE(LUWRT,*) '  MY_VEC1_IOFF :', MY_VEC1_IOFF
        WRITE(LUWRT,*) '  MY_VEC2_IOFF :', MY_VEC2_IOFF
        WRITE(LUWRT,*) '  MY_ACT_BLK1  :', MY_ACT_BLK1
        WRITE(LUWRT,*) '  MY_ACT_BLK2  :', MY_ACT_BLK2
      END IF
C
 101  CONTINUE     
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE GROUP_CONSTRUCTOR(IGROUPLIST,IPROCLIST)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
      DIMENSION IPROCLIST(LUCI_NMPROC),IGROUPLIST(LUCI_NMPROC)
      INTEGER IKEY, ICOLOR, JKEY, JCOLOR
      ITEST = 00
csk   ITEST = 10
C
C     'intra-node' communicator MYNEW_COMM ( I/O communicator )
C
      IKEY = LUCI_MYPROC
      ICOLOR = IPROCLIST(LUCI_MYPROC+1)
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,ICOLOR,IKEY,MYNEW_COMM,IERR)
C
C     collect useful information about each group,
C     store on common block
C
      NEWCOMM_PROC = 0
      MYNEW_ID = 0
      CALL MPI_COMM_SIZE(MYNEW_COMM,NEWCOMM_PROC,IERR)
      CALL MPI_COMM_RANK(MYNEW_COMM,MYNEW_ID,IERR)
C
C     inter-node communicator ICOMM
C
      IF( MYNEW_ID .eq. 0 ) THEN
        JKEY = IPROCLIST(LUCI_MYPROC+1)
        JCOLOR = 2
      ELSE
        JKEY = IPROCLIST(LUCI_MYPROC+1)
        JCOLOR = 3
      END IF
C
      CALL MPI_COMM_SPLIT(MPI_COMM_WORLD,JCOLOR,JKEY,ICOMM,IERR)
C
C     collect again ...
C
      ICOMM_ID = 0
      ICOMM_SIZE = 0
      CALL MPI_COMM_RANK(ICOMM,ICOMM_ID,IERR)
      CALL MPI_COMM_SIZE(ICOMM,ICOMM_SIZE,IERR)
C
C     set-up and store group information
C
      CALL SET_GROUP_TABLE(IGROUPLIST,IPROCLIST,ICOLOR)
C
C     store personal group number
C
      MY_GROUPN = ICOLOR
      IF ( ITEST .ge. 10 )THEN
        WRITE(LUOUT,*) ' '
        WRITE(LUOUT,*) ' OUTPUT FROM GROUP_CONSTRUCTOR'
        WRITE(LUOUT,*) ' '
        WRITE(LUOUT,*) ' size of MYNEW_COMM     :',NEWCOMM_PROC
        WRITE(LUOUT,*) ' size of ICOMM          :',ICOMM_SIZE
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE GROUP_DESTRUCTOR(JCOMM1,JCOMM2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
      CALL MPI_COMM_FREE(JCOMM1,IERR)
      CALL MPI_COMM_FREE(JCOMM2,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SET_GROUP_TABLE(IGROUPLIST,IPROCLIST,ICOLOR)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
C     INPUT
      DIMENSION IPROCLIST(LUCI_NMPROC)
C     OUTPUT
      DIMENSION IGROUPLIST(LUCI_NMPROC)
C
      INUMB = 1
      DO IPROC = 1, LUCI_NMPROC
C
        IF( IPROCLIST(IPROC) .eq. ICOLOR ) THEN
C         put in the process tag ( master is 0 )
          IGROUPLIST(INUMB) = IPROC - 1
          INUMB = INUMB + 1
        END IF
C
      END DO
C
      IF( INUMB-1 .gt. NEWCOMM_PROC ) THEN
        WRITE(LUOUT,*) 'Error in SET_GROUP_TABLE: more CPUs assigned
     & as included in this group!'
        WRITE(LUOUT,*) 'assigned CPUs, group size',INUMB-1,NEWCOMM_PROC
        CALL Abend2('Error detected in SET_GROUP_TABLE')
      END IF
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_OPEN(NFILE_ID)
C
C     OUTPUT
C     ======
C
C     open files ready for MPI-I/O
C     files handles stored on common block LUCIPAR
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      CHARACTER*10 PLU1,PLU2,PLU3,PLU4,PLU5,PLU6,PLU7,PLU8,PLUC,PDIA
      CHARACTER*6 PLU1BL,PLU2BL,PLU3BL,PLU4BL,PLU6BL,PLU7BL,PLU5BL
      CHARACTER*6 PLUDIA,PLUCBL
      INTEGER PFILELENGTH
      INTEGER(KIND=MPI_OFFSET_KIND) IDISP
      IDISP = 0
      PFILELENGTH = 0
C
      PLUDIA = 'DIAPAR'
      PLU1BL = 'LU1PAR'
      PLU2BL = 'LU2PAR'
      PLU3BL = 'LU3PAR'
      PLU4BL = 'LU4PAR'
      PLU5BL = 'LU5PAR'
      PLU6BL = 'LU6PAR'
      PLU7BL = 'LU7PAR'
      PLUCBL = 'LUCPAR'
C
C     set correct names
C     -----------------
      IF( NFILE_ID .lt. 10 ) THEN    ! MPI ID has one digit
        WRITE (PDIA,'(A6,A1,I1)') PLUDIA,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I1)') PLUCBL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I1)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I1)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I1)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I1)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I1)') PLU5BL,'.',NFILE_ID
        WRITE (PLU6,'(A6,A1,I1)') PLU6BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I1)') PLU7BL,'.',NFILE_ID
        PFILELENGTH = 8
      ELSE IF( NFILE_ID .lt. 100)THEN  ! MPI ID has two digits
        WRITE (PDIA,'(A6,A1,I2)') PLUDIA,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I2)') PLUCBL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I2)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I2)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I2)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I2)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I2)') PLU5BL,'.',NFILE_ID
        WRITE (PLU6,'(A6,A1,I2)') PLU6BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I2)') PLU7BL,'.',NFILE_ID
        PFILELENGTH = 9
      ELSE                        ! MPI ID has three digits
        WRITE (PDIA,'(A6,A1,I3)') PLUDIA,'.',NFILE_ID
        WRITE (PLUC,'(A6,A1,I3)') PLUCBL,'.',NFILE_ID
        WRITE (PLU1,'(A6,A1,I3)') PLU1BL,'.',NFILE_ID
        WRITE (PLU2,'(A6,A1,I3)') PLU2BL,'.',NFILE_ID
        WRITE (PLU3,'(A6,A1,I3)') PLU3BL,'.',NFILE_ID
        WRITE (PLU4,'(A6,A1,I3)') PLU4BL,'.',NFILE_ID
        WRITE (PLU5,'(A6,A1,I3)') PLU5BL,'.',NFILE_ID
        WRITE (PLU6,'(A6,A1,I3)') PLU6BL,'.',NFILE_ID
        WRITE (PLU7,'(A6,A1,I3)') PLU7BL,'.',NFILE_ID
        PFILELENGTH = 10
      END IF
C
C     ... open
C
      CALL MPI_FILE_OPEN(MYNEW_COMM,PDIA(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,IDIA,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU1(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU1,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU2(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU2,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU3(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU3,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU4(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU4,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU5(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU5,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU6(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU6,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLU7(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILU7,IERR)
      CALL MPI_FILE_OPEN(MYNEW_COMM,PLUC(1:PFILELENGTH),
     &     MPI_MODE_CREATE + MPI_MODE_RDWR,MPI_INFO_NULL,ILUC,IERR)
C
C     ... set fileview
C
      CALL MPI_FILE_SET_VIEW(IDIA,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU1,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU2,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU3,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU4,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU5,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU6,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILU7,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
      CALL MPI_FILE_SET_VIEW(ILUC,IDISP,MPI_REAL8,MPI_REAL8,
     &        "native",MPI_INFO_NULL,IERR)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE SETUNITS_PAR_CLOSE
C
C     OUTPUT
C     ======
C
C     close files
C     files handles stored on common block LUCIPAR are set to MPI_NULL
C
C     Last revision:     S. Knecht       - March  2008
C
************************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
C
C     ... close
C
      CALL MPI_FILE_CLOSE(IDIA,IERR)
      CALL MPI_FILE_CLOSE(ILU1,IERR)
      CALL MPI_FILE_CLOSE(ILU2,IERR)
      CALL MPI_FILE_CLOSE(ILU3,IERR)
      CALL MPI_FILE_CLOSE(ILU4,IERR)
      CALL MPI_FILE_CLOSE(ILU5,IERR)
      CALL MPI_FILE_CLOSE(ILU6,IERR)
      CALL MPI_FILE_CLOSE(ILU7,IERR)
      CALL MPI_FILE_CLOSE(ILUC,IERR)
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
#include "infpar.h"
#include "mpif.h"
#include "parluci.h"
#include "wrkspc.inc"
      INTEGER*8 IABSOLUTE_WEIGHT
      IABSOLUTE_WEIGHT = 0
C
      IDUM  = 0
      CALL MEMMAN(KDUM,IDUM,'MARK  ',IDUM,'BLKDRV')
C
C     allocate local scratch arrays
C
      CALL MEMMAN(KICCTOS,  NBLOCK**2,'ADDL  ',1,'ICCTOS')
      CALL MEMMAN(KCWEIGHT, NBLOCK   ,'ADDL  ',1,'ICWHT ')
      CALL MEMMAN(KCWEIGHTF,NBLOCK   ,'ADDL  ',1,'ICWHTF')
      CALL MEMMAN(KBLCKWT,2*LUCI_NMPROC,'ADDL  ',1,'IBLCKW')
      CALL MEMMAN(KIBTOTW,  NBLOCK   ,'ADDL  ',2,'IBTOTW')
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
        CALL DISTBLKND_2(NBLOCK,WORK(KCWEIGHTF),NBLOCKD,IBLOCKL)
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_BLOCK_DISTR(JCTCMP,MNSBLOCK,JISBLOCK,
     &           JICOCOC,JICSMOS,JICBLTP,MNSSOA,MNSSOB,
     &           MNOCTPA,MNOCTPB,
     &           MNSMST,JLC,JIDC,JIDOH2,JISTRFL,
     &           PPS,LLCBLOCK,LLECBLOCK,JI1CBLOCK,JICBLOCK,
     &           JIRESTRICT,JICONSPA,JICONSPB,TSCLFAC,
     &           JICBAT_RES,LCBAT_INI,LCBAT_END,JCBLINF,
     &           MNJJJBLK,MTEMPP,JMPINFA,LLUC,CCB,CC2,CC,NYTOTB,
     &           ICSBTCHB)
*
*
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
*. Specific input
      INTEGER JISBLOCK(8,*)
*.General input
      INTEGER JICOCOC(MNOCTPA,MNOCTPB)
      INTEGER JICSMOS(MNSMST)
      INTEGER JICBLTP(*)
      INTEGER MNSSOA(MNSMST ,*), MNSSOB(MNSMST ,*)
      INTEGER JICONSPA(MNOCTPA,MNOCTPA), JICONSPB(MNOCTPB,MNOCTPB)
*.Scratch
      INTEGER   LLCBLOCK(*),JI1CBLOCK(*),JICBLOCK(8,*),LLECBLOCK(*)
*
      DIMENSION LASM(4),LBSM(4),LATP(4),LBTP(4),LSGN(5),LTRP(5)
      DIMENSION TSCLFAC(*)
* new scratch
      INTEGER MNJJJBLK,ILOOPSB,MTEMPP,ISCRLNGTH
      DIMENSION JCBLINF(MNJJJBLK)
      DIMENSION CCB(*),CC2(*)
* switch between 'send' and 'broadcast' of the CB-array
*      INTEGER I_USE_BCAST
*  ==========
*.   output 
*  ==========
*
      DIMENSION JCTCMP(MNJJJBLK,4), JMPINFA(MTEMPP,2), NYTOTB(MNJJJBLK)
      DIMENSION ICSBTCHB(MNJJJBLK)
*   
* JCTCMP: array containing information about a c-block, whether it is
* contributing, together with its length; 
* third number specifies how often it interacts with a given s-block    
*
*.
      ISCRLNGTH = 0
      ILOOPSB = 0
*      I_USE_BCAST = 1
CCC!      IF(LUCI_MYPROC.EQ.LUCI_MASTER)THEN
CCC!        REWIND LLUC
CCC!      END IF
*      WRITE(6,*)'JIDC', JIDC
* ===========================
* 1 : Arrays for accessing C
* ============================
*. Find batches of C - strings
*. Should be changed to REPART or ??
*. Not neccesary, all zero blocks detected otherwise
      ITTSS_ORD = 2
      CALL PART_CIV2(JIDC,JICBLTP,MNSSOA,MNSSOB,MNOCTPA,MNOCTPB,MNSMST,
     &              JLC,JICOCOC,JICSMOS,MCBATCH,LLCBLOCK,LLECBLOCK,
     &              JI1CBLOCK,JICBLOCK,0,ITTSS_ORD)
*
*
* Loop over batches over C blocks
      IF(JIDOH2.EQ.1) THEN
        MXEXC  = 2
      ELSE
        MXEXC = 1
      END IF
      IF(JICBAT_RES.EQ.1) THEN
        WRITE(6,*) ' Restricted set of C batches '
        WRITE(6,*) ' LCBAT_INI LCBAT_END', LCBAT_INI,LCBAT_END
        JCBAT_INI = LCBAT_INI
        JCBAT_END = LCBAT_END
      ELSE
        JCBAT_INI = 1
        JCBAT_END = MCBATCH
      END IF
*
      DO 20000 JCBATCH = JCBAT_INI,JCBAT_END

*. Read C blocks into core
*
        ICOFF = 1
        NJBLOCK = LLCBLOCK(JCBATCH)
*        WRITE(6,*)'MCBATCH = ',MCBATCH
*        WRITE(6,*)'NJBLOCK =',NJBLOCK
        DO JJCBLOCK = 1, NJBLOCK
          JBLOCK = JI1CBLOCK(JCBATCH)-1+JJCBLOCK
*. Will this block be needed ??
          INTERACT = 0
*          WRITE(6,*)'TSCLFAC(JBLOCK) =',TSCLFAC(JBLOCK)
          IF(TSCLFAC(JBLOCK).EQ. 1.0D0) THEN
C?          WRITE(6,*) ' Nonvanishing block ', JBLOCK, myproc
          JATP = JICBLOCK(1,JBLOCK)
          JBTP = JICBLOCK(2,JBLOCK)
          JASM = JICBLOCK(3,JBLOCK)
          JBSM = JICBLOCK(4,JBLOCK)
          JOFF = JICBLOCK(5,JBLOCK)
          CALL PRMBLK(JIDC,JISTRFL,JASM,JBSM,JATP,JBTP,PPS,PL,
     &                LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
          DO IPERM = 1, NPERM
            LLASM = LASM(IPERM)
            LLBSM = LBSM(IPERM)
            LLATP = LATP(IPERM)
            LLBTP = LBTP(IPERM)
*.Loop over Sigma blocks in batch
            DO JSBLOCK = 1, MNSBLOCK
            IDENT = 0
            IF(JISBLOCK(1,JSBLOCK).GT.0) THEN
              IATP = JISBLOCK(1,JSBLOCK)
              IBTP = JISBLOCK(2,JSBLOCK)
              IASM = JISBLOCK(3,JSBLOCK)
              IBSM = JISBLOCK(4,JSBLOCK)
*. Are the two blocks connected by allowed excitation
              IF(MXEXC.EQ.1) THEN
                 IF((JICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (JICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((JICONSPA(IATP,LLATP).LE.1.AND.
     &               JICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (JICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (JICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
*. Are they identical ?
              IDENT = 0
              IF(IASM.EQ.JASM.AND.IATP.EQ.JATP.AND.
     &           IBSM.EQ.JBSM.AND.IBTP.EQ.JBTP) IDENT = 1
*
            END IF
            END DO
          END DO
*.        ^ End of checking whether C-block is needed
          END IF
*         ^ Checking was only done for nonvanishing blocks
*
          ISCALE = 0
          IF(INTERACT.EQ.1) THEN
*
            JCTCMP(JBLOCK,1) = 1
*          write(6,*)'active: JCTCMP(JBLOCK,1) =',JCTCMP(JBLOCK,1),myproc
            JCTCMP(JBLOCK,2) = JCBLINF(JBLOCK)
*          write(6,*)'active length: JCTCMP(JBLOCK,2) =',JCTCMP(JBLOCK,2)
          ELSE
*
            JCTCMP(JBLOCK,1) = 0
*            write(6,*)'inact: JCTCMP(JBLOCK,1) =',JCTCMP(JBLOCK,1)
            JCTCMP(JBLOCK,2) = JCBLINF(JBLOCK)
            TSCLFAC(JBLOCK) = 0.0D0
*
          END IF
*
        END DO
*
*
       DO 9000 ICBLK = JI1CBLOCK(JCBATCH), JI1CBLOCK(JCBATCH)-1+NJBLOCK
*
*            WRITE(6,*)'ICBLK 9000loop =',ICBLK
            JATP = JICBLOCK(1,ICBLK)
            JBTP = JICBLOCK(2,ICBLK)
            JASM = JICBLOCK(3,ICBLK)
            JBSM = JICBLOCK(4,ICBLK)
            ICOFF = JICBLOCK(5,ICBLK)
            NJA = MNSSOA(JASM,JATP)
*            WRITE(6,*)'NJA =',NJA
            NJB = MNSSOB(JBSM,JBTP)
*
            IF(TSCLFAC(ICBLK).NE.0.0D0) THEN
*. Other symmetry blocks that can be obtained from this block
              CALL PRMBLK(JIDC,JISTRFL,JASM,JBSM,JATP,JBTP,PPS,PL,
     &                    LATP,LBTP,LASM,LBSM,LSGN,LTRP,NPERM)
*. Start with transposed block
                ILOOPSB = 0
*                WRITE(6,*)'NPERM 9000loop=',NPERM
*
              DO 8765 IPERM = NPERM,1, -1
                LLASM = LASM(IPERM)
                LLBSM = LBSM(IPERM)
                LLATP = LATP(IPERM)
                LLBTP = LBTP(IPERM)
                NLLA = MNSSOA(LLASM,LLATP)
                NLLB = MNSSOB(LLBSM,LLBTP)
*
                DO 10000 ISBLK = 1, MNSBLOCK
                IF(JISBLOCK(1,ISBLK) .GT. 0 ) THEN
                  IATP = JISBLOCK(1,ISBLK)
                  IBTP = JISBLOCK(2,ISBLK)
                  IASM = JISBLOCK(3,ISBLK)
                  IBSM = JISBLOCK(4,ISBLK)
                  ISOFF = JISBLOCK(5,ISBLK)
                  NIA = MNSSOA(IASM,IATP)
                  NIB = MNSSOB(IBSM,IBTP)
*
                  IF(NIA*NIB.EQ.0) GOTO 9999
                  IF(JIRESTRICT.EQ.1.AND.
     &               (JASM.GT.IASM.OR.
     &               JASM.EQ.IASM.AND.JATP.GT.IATP.OR.
     &               JASM.EQ.IASM.AND.JATP.EQ.IATP.AND.JBTP.GT.IBTP))
     &            GOTO 9999
*. Are the two blocks connected by allowed excitation
              INTERACT = 0
              IF(MXEXC.EQ.1) THEN
                 IF((JICONSPA(IATP,LLATP).LE.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM  ) .OR.
     &              (JICONSPB(IBTP,LLBTP).LE.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM  )     )THEN
                       INTERACT = 1
                 ENDIF
              ELSE IF(MXEXC.EQ.2) THEN
                 IF((JICONSPA(IATP,LLATP).LE.1.AND.
     &               JICONSPB(IBTP,LLBTP).LE.1     )   .OR.
     &              (JICONSPA(IATP,LLATP).EQ.MXEXC.AND.
     &               IBTP.EQ.LLBTP.AND.IBSM.EQ.LLBSM) .OR.
     &              (JICONSPB(IBTP,LLBTP).EQ.MXEXC.AND.
     &               IATP.EQ.LLATP.AND.IASM.EQ.LLASM)     )THEN
                        INTERACT = 1
                 END IF
              END IF
C
                 IF(INTERACT.EQ.0) GOTO 9999
*
                 ILOOPSB = ILOOPSB + 1
*              INSERT SOMETHING HERE !!!
*
*
              END IF
 9999         CONTINUE
              ILOOPSB = ILOOPSB + 0
10000         CONTINUE
*./\  End of loop over sigma blocks
 8765         CONTINUE
            END IF
*
            ILOOPSB = ILOOPSB * NPERM
*
            JCTCMP(ICBLK,3) = ILOOPSB
*           if c-block is interacting and at least one s-block then 
*           NYTOTB(ICBLK) > 0, else 0
            NYTOTB(ICBLK) = 
     &      JCTCMP(ICBLK,3)*JCTCMP(ICBLK,2)*JCTCMP(ICBLK,1) 
*
*
            ICSBTCHB(ICBLK) = ICSBTCHB(ICBLK) + NYTOTB(ICBLK)
*
*            WRITE(6,*)'NYTOTB(ICBLK) =',NYTOTB(ICBLK), myproc
*            WRITE(6,*)'ICSBTCHB(ICBLK) =',ICSBTCHB(ICBLK), myproc
            
 9000       CONTINUE
*. /\ End of loop over C blocks in Batch
*. End of loop over S blocks in batch
20000   CONTINUE
*. End of loop over batches of C blocks

*. we wanna find a distribution of the c-blocks among the nodes (incl. master)
*. as even as possible, here we go ! 
*
      IF(ICOMPUTE.GT.0.AND.IUSE_VER.EQ.2) THEN
        IF(NEXCL_MR.EQ.0) THEN
*
          CALL FIND_EVEN_DISTR(JCTCMP,MNJJJBLK,JCBLINF,MTEMPP,
     &         JMPINFA,NYTOTB)
*
        ELSE IF(NEXCL_MR.GT.0) THEN
*
          CALL FIND_EVEN_DISTR_EXM(JCTCMP,MNJJJBLK,JCBLINF,MTEMPP,
     &         JMPINFA,NYTOTB)
*
        END IF
      END IF
*
      RETURN 
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_SBATCH_DISTR(JCSBTCHB,NNJJJBLK,NTEMPP,IDISTLIST)
*

      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
*
*. input
* NNJJJBLK: total number of c-blocks
* NTEMPP: number of processes available for the distribution
* ICSBTCHB(NNJJJBLK): total length of each block
*. output
* IDISTLIST(NNJJJBLK): list containing the number of a node  
* that computes the block
*
      INTEGER NNJJJBLK, ICOUNTABLK
      INTEGER NTEST
      INTEGER*8 ITOTBLCKL,MMPINFA
      DIMENSION MMPINFA(NTEMPP,2)
      DIMENSION IDISTLIST(NNJJJBLK), JCSBTCHB(NNJJJBLK)
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: actual number of the block showing the LONGEST! :)
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*      DATA MXSIZE,ITEMPW1 /20000000000000000,20000000000000000/
*
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 1900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1900000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO I = 1, NTEMPP
       MMPINFA(I,1) = 0
       MMPINFA(I,2) = 0
      END DO
*. end of initializtion
      DO ICBL = 1, NNJJJBLK
*
        IF(JCSBTCHB(ICBL).NE.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
*           IF(MYNEW_ID.EQ.0)THEN
*             WRITE(6,*) 'This is a active block with length',
*     &ICBL,MCTCMP(ICBL,2)
*           END IF
           ITOTBLCKL = ITOTBLCKL + JCSBTCHB(ICBL)
        END IF
*        IF(MYNEW_ID.EQ.0)THEN
*           WRITE(6,*)'JCSBTCHB(ICBL) is',JCSBTCHB(ICBL)
*        END IF
*
      END DO
*
      IF(NTEST.GT.0) THEN
       IF(MYNEW_ID.EQ.0) THEN
         WRITE(6,*) '  total number of processes to distribute on:',
     &                 NTEMPP -1
         WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
         WRITE(6,*) '  overall active block length:', ITOTBLCKL
       END IF
      END IF
*
      IF(ICOUNTABLK.LT.NTEMPP-1) THEN
        IMINNP = MIN(ICOUNTABLK,NTEMPP)
        IF(NTEST.GT.0) THEN
          IF(MYNEW_ID.EQ.0) THEN
            WRITE(6,*) 'subroutine FIND_SBATCH_DISTR speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP+1,'procs'
          END IF
        END IF
        NTEMPP = IMINNP+1
      END IF
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
        DO 3000 II = 1, NNJJJBLK
*
          ITEMPL = JCSBTCHB(II)
          ITEMPD = IDISTLIST(II)
          ITEMPN = II

*          WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
*          WRITE(6,*) 'MXSIZE:',MXSIZE
*          WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.LT.0) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
*                  WRITE(6,*) 'MXNUMB',MXNUMB
*                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
*        WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 2, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR,2)
*          WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> MYNEW_ID = 0 !
          IAMTPROC = IPR
*          WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.2) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
*             WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
       IF(NTEST.GT.0) THEN
       WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
       IDISTLIST(MXNUMB)  = IAMGPROC - 1
*
*.      raising MMPINFA(proc,1) by 1
*.      adding on MMPINFA(proc,2) the weight of the new block
        MXSIZE = MXSZTMP
        MMPINFA(IAMGPROC,1) = MMPINFA(IAMGPROC,1) + 1
        MMPINFA(IAMGPROC,2) = MMPINFA(IAMGPROC,2) + MXSIZE
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
      IF(NTEST.GT.0) THEN
        IF(MYNEW_ID.EQ.0) THEN
      WRITE(6,'(3X,A,I4)')
     &'SUMMATION OF EVEN DISTRIBUTION OF LUCI_MYPROC:',MYNEW_ID
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MMPINFA(ISTI,1),' blocks with a
     & total length of',MMPINFA(ISTI,2)
      END DO
*
*      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
      END IF
*
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_EVEN_DISTR(MCTCMP,NNJJJBLK,MCBLINF,NTEMPP,
     &           MMPINFA,JYTOTB)
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
*
*. input
* NNJJJBLK: total number of c-blocks 
* MCBLINF(NNJJJBLK): length of c-block
* NTEMPP: number of processes available for the distribution
* MCTCMP: Blocks info array 
*   MCTCMP(*,1) : block contributes 1 or not 0
*   MCTCMP(*,2) : block length
*   MCTCMP(*,3) : count of interacts with s-blocks
*   MCTCMP(*,4) : node number that will get the block; default is -1
*   JYTOTB(*) : length of ACTIVE c-block * count of interacts with s-blocks 
*
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks 
      INTEGER NNJJJBLK, ICOUNTABLK 
      INTEGER NTEST 
      INTEGER*8 ITOTBLCKL
      DIMENSION MCTCMP(NNJJJBLK,4), MCBLINF(NNJJJBLK)
      DIMENSION MMPINFA(NTEMPP,2), JYTOTB(NNJJJBLK) 
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: actual number of the block showing the LONGEST! :)
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 1900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1900000000 
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO ICBL = 1, NNJJJBLK
*
        IF(MCTCMP(ICBL,1).NE.0) THEN 
           ICOUNTABLK = ICOUNTABLK + 1
*           IF(LUCI_MYPROC.EQ.LUCI_MASTER)THEN
*             WRITE(6,*) 'This is a active block with length',
*     &ICBL,MCTCMP(ICBL,2)
*           END IF
           ITOTBLCKL = ITOTBLCKL + JYTOTB(ICBL)
        END IF
*        WRITE(6,*)'MCTCMP(ICBL,4) is',MCTCMP(ICBL,4) 
*
      END DO 
*
      IF(NTEST.GT.0) THEN
       WRITE(6,*) '  total number of processes to distribute on:',NTEMPP
       WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
       WRITE(6,*) '  overall active block length:', ITOTBLCKL
      END IF 
*
      IF(ICOUNTABLK.LT.NTEMPP) THEN
        IMINNP = MIN(ICOUNTABLK,NTEMPP)
        IF(NTEST.GT.0) THEN
          WRITE(6,*) 'subroutine FIND_EVEN... speaking!'
          WRITE(6,*) 'number of active blocks lower than number of 
     & procs --> I only use',IMINNP,'procs'
        END IF
        NTEMPP = IMINNP
      END IF
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
        DO 3000 II = 1, NNJJJBLK
*
          
          ITEMPL = JYTOTB(II)
          ITEMPD = MCTCMP(II,4)
          ITEMPN = II
*          WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
*          WRITE(6,*) 'MXSIZE:',MXSIZE
*          WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.LT.0) THEN
              IF(ITEMPL.LE.MXSIZE)THEN 
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
*                  WRITE(6,*) 'MXNUMB',MXNUMB
*                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
*        WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 1, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR,2)
*          WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> LUCI_MASTER !
          IAMTPROC = IPR 
*          WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.1) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC 
*             WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF   
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M incl.)
       IF(NTEST.GT.0) THEN
       WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
        MCTCMP(MXNUMB,4) = IAMGPROC - 1
*              
*.      raising MMPINFA(proc,1) by 1
*.      adding on MMPINFA(proc,2) the weight of the new block      
        MXSIZE = MXSZTMP
        MMPINFA(IAMGPROC,1) = MMPINFA(IAMGPROC,1) + 1
        MMPINFA(IAMGPROC,2) = MMPINFA(IAMGPROC,2) + MXSIZE
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
      IF(NTEST.GT.0) THEN
        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
      WRITE(6,'(3X,A,I4)')
     &'SUMMATION OF EVEN DISTRIBUTION OF LUCI_MYPROC:',LUCI_MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MMPINFA(ISTI,1),' blocks with a  
     & total length of',MMPINFA(ISTI,2)
      END DO
*
      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_EVEN_DISTR_EXM(MCTCMP,NNJJJBLK,MCBLINF,NTEMPP,
     &           MMPINFA,JYTOTB)

*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
*
*. input
* NNJJJBLK: total number of c-blocks
* MCBLINF(NNJJJBLK): length of c-block
* NTEMPP: number of processes available for the distribution
* MCTCMP: Blocks info array
*   MCTCMP(*,1) : block contributes 1 or not 0
*   MCTCMP(*,2) : block length
*   MCTCMP(*,3) : count of interacts with s-blocks
*   MCTCMP(*,4) : length of ACTIVE c-block * count of interacts with s-blocks
*   MCTCMP(*,5) : node number that will get the block; default is -1
*
*   MMPINFA(*,1): total number of blocks for the n-th proc
*   MMPINFA(*,2): total length ('weight') of the blocks
      INTEGER NNJJJBLK, ICOUNTABLK
      INTEGER NTEST
      INTEGER*8 ITOTBLCKL
      DIMENSION MCTCMP(NNJJJBLK,4), MCBLINF(NNJJJBLK)
      DIMENSION MMPINFA(NTEMPP,2), JYTOTB(NNJJJBLK)
*.scratch
* MXSIZE: actual maximum size for one of all blocks
* MXNUMB: actual number of the block showing the LONGEST! :)
      INTEGER*8 MXSIZE, MXSZTMP, ITEMPL, ITEMPW1, ITEMPW2
      INTEGER ITEMPN, IAMGPROC, IAMTPROC, MXNUMB, IRUN
*
      NTEST = 0
      ICOUNTABLK = 0
      MXSIZE = 1900000000
      MXSZTMP = 0
      MXNUMB = 0
      IRUN = 0
      ITEMPW1 = 1900000000
      ITEMPW2 = 0
      ITEMPN = 0
      ITEMPL = 0
      IAMGPROC = 0
      IAMTPROC = 0
      ITOTBLCKL = 0
      DO ICBL = 1, NNJJJBLK
*
        IF(MCTCMP(ICBL,1).NE.0) THEN
           ICOUNTABLK = ICOUNTABLK + 1
*           IF(LUCI_MYPROC.EQ.LUCI_MASTER)THEN
*             WRITE(6,*) 'This is a active block with length',
*     &ICBL,MCTCMP(ICBL,2)
*           END IF
           ITOTBLCKL = ITOTBLCKL + JYTOTB(ICBL)
        END IF
*        WRITE(6,*)'MCTCMP(ICBL,5) is',MCTCMP(ICBL,5)
*
      END DO
*
*      IF(NTEST.GT.0) THEN
       IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
        WRITE(6,'(3X,A,1X,I4)') 'total number of processes to 
     &distribute on:',NTEMPP-1
        WRITE(6,*) '  total number of active blocks:', ICOUNTABLK
        WRITE(6,*) '  overall active block length:', ITOTBLCKL
       END IF
*      END IF
*
      IF(ICOUNTABLK.LT.NTEMPP-1) THEN
        IMINNP = MIN(ICOUNTABLK,NTEMPP)
*        IF(NTEST.GT.0) THEN
          IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
            WRITE(6,*) 'subroutine FIND_EVEN...EXM speaking!'
            WRITE(6,*) 'number of active blocks lower than number of
     & procs --> I only use',IMINNP+1,'procs'
         END IF
*        END IF
        NTEMPP = IMINNP+1
      END IF
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
        DO 3000 II = 1, NNJJJBLK
*
          ITEMPL = JYTOTB(II)
          ITEMPD = MCTCMP(II,4)
          ITEMPN = II

*          WRITE(6,*) 'ITEMPL,ITEMPD,ITEMPN:',ITEMPL,ITEMPD,ITEMPN
*
*          WRITE(6,*) 'MXSIZE:',MXSIZE
*          WRITE(6,*) 'MXSZTMP:',MXSZTMP
          IF(ITEMPL.GT.0)THEN
            IF(ITEMPD.LT.0) THEN
              IF(ITEMPL.LE.MXSIZE)THEN
                IF(ITEMPL.GT.MXSZTMP) THEN
                   MXSZTMP = ITEMPL
                   MXNUMB = ITEMPN
*                  WRITE(6,*) 'MXNUMB',MXNUMB
*                  WRITE(6,*) 'MXSZTMP (2):',MXSZTMP
                END IF
              END IF
            END IF
          END IF
*
3000    CONTINUE
*
*        WRITE(6,*)'NTEMPP:',NTEMPP
        DO 4000 IPR = 2, NTEMPP
*
          ITEMPW2 = MMPINFA(IPR,2)
*          WRITE(6,*) 'ITEMPW2 and ITEMPW1:',ITEMPW2,ITEMPW1
* attention: IAMTPROC = 1 --> LUCI_MASTER !
          IAMTPROC = IPR
*          WRITE(6,*) 'IAMTPROC',IAMTPROC
*
          IF(IPR.EQ.2) ITEMPW1 = ITEMPW2
*
          IF(ITEMPW2.LE.ITEMPW1) THEN
             ITEMPW1 = ITEMPW2
             IAMGPROC = IAMTPROC
*             WRITE(6,*)'ITEMPW1,IAMGPROC',ITEMPW1,IAMGPROC
          END IF
*
4000    CONTINUE
*. now we should have found a proc and a block --> put both together !
*
*.      calculation of block MXNUMB by proc IAMGPROC -1 (M excl.)
       IF(NTEST.GT.0) THEN
       WRITE(6,*)'calculation of block MXNUMB by proc',MXNUMB,IAMGPROC-1
       ENDIF
        MCTCMP(MXNUMB,4) = IAMGPROC - 1
*
*.      raising MMPINFA(proc,1) by 1
*.      adding on MMPINFA(proc,2) the weight of the new block
        MXSIZE = MXSZTMP
        MMPINFA(IAMGPROC,1) = MMPINFA(IAMGPROC,1) + 1
        MMPINFA(IAMGPROC,2) = MMPINFA(IAMGPROC,2) + MXSIZE
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
*      IF(NTEST.GT.0) THEN
        IF(LUCI_MYPROC.EQ.LUCI_MASTER) THEN
      WRITE(6,'(3X,A,I4)')
     &'SUMMATION OF EVEN DISTRIBUTION OF LUCI_MYPROC:',LUCI_MYPROC
      DO ISTI = 1, NTEMPP
        WRITE(6,'(3X,A,I2,A,I5,A,1X,I17)')
     & 'process',ISTI -1,' calculates',MMPINFA(ISTI,1),' blocks with a
     & total length of',MMPINFA(ISTI,2)
      END DO
*
      WRITE(6,'(3X,A)') 'even distribution finished!'
*
        END IF
*      END IF
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
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
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
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
     &                             MPI_REAL8,ISTAT,IERR)
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
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_SPAR_CCOEFF(LCFILE,LUOUT,SCR,SCALFAC,
     &                              IBLOCKL,IBLOCKD,NDIM,IGROUPLIST)
C
C     copy vector from serial I/O file LCFILE to 
C     parallel I/O file MPILUSCR
C
      IMPLICIT REAL*8(A-H,O-Z)
*
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      INTEGER(KIND=MPI_OFFSET_KIND) IOFFSET
C     INPUT 
      DIMENSION SCR(*), IBLOCKL(NDIM), IBLOCKD(NDIM)
      DIMENSION IGROUPLIST(LUCI_NMPROC)
C     OUTPUT
      DIMENSION SCALFAC(NDIM)
C     SCRATCH
      INTEGER MIN2
      MIN2 = -2
C
C     Initialize SCALFAC to 0.0D0
C
      ZERO = 0.0D0
      CALL SETVEC(SCALFAC,ZERO,NDIM) 
CSKT      write(6,*) 'I arrived at the copy routine',myproc
C
C     file offset is at 0
      IOFFSET = 0
C
C     Loop over blocks
C
      DO IBLK = 1, NDIM
C
        LBL = IBLOCKL(IBLK)
        NDBLOCK = IBLOCKD(IBLK)
C
        IF( NDBLOCK .eq. MIN2 ) GOTO 500
        IF( NDBLOCK .eq. LUCI_MYPROC ) GOTO 100
        IGROUPBLK = 0
        IF(NEWCOMM_PROC .ne. LUCI_NMPROC ) THEN
          DO IPROC = 1, NEWCOMM_PROC
            IF( NDBLOCK .eq. IGROUPLIST(IPROC) ) IGROUPBLK = 1 
          END DO
          IF( IGROUPBLK .eq. 0) GOTO 500
        END IF
C       ... well this nonzero block is at least inside this file-group
 100    CONTINUE
        IF( LUCI_MYPROC .eq. NDBLOCK ) THEN
          CALL IFRMDS(LBL,1,-1,LCFILE)
          CALL ITODS (LBL,1,-1,LUOUT)
          KBLK = -1
          NO_ZEROING = 1
          CALL FRMDSC2(SCR,LBL,KBLK,LCFILE,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          IF( IMZERO .eq. 0 ) THEN
            SCALFAC(IBLK) = 1.0D0
            CALL MPI_FILE_WRITE_AT(MPILUSCR2,IOFFSET,SCR,LBL,MPI_REAL8,
     &           MPI_STATUS_IGNORE,IERR)
CSKT           WRITE(6,*) 'THIS IS WHAT I HAVE WRITTEN',LUCI_MYPROC
CSKT           WRITE(6,*) 'IOFFSET',IOFFSET
CSKT            CALL WRTMAT(SCR,1,LBL,1,LBL)
            IF(IAMPACK.EQ.0) THEN
              CALL TODSC (SCR,LBL,KBLK,LUOUT)
            ELSE
              CALL TODSCP(SCR,LBL,KBLK,LUOUT)
            END IF
          ELSE
            CALL ZERORC(LBL,LUOUT,IAMPACK)
          END IF
CSK          IOFFSET = IOFFSET + LBL
        ELSE
CSK        IOFFSET = IOFFSET + LBL
        END IF
 500    CONTINUE
        IOFFSET = IOFFSET + LBL
      END DO
C
C     find EOF-mark
C
      IF( NROOT_INFO .gt. 0 ) THEN
        CALL IFRMDS(LBL,1,-1,LCFILE)
CSK        IF(LBL.LT.0) WRITE(LUDUMPF,*)'I am fine'
        IF(LBL.GE.0) THEN
          WRITE(6,*) 'Error in COPVCD_SPAR_CCOEFF detected'
          CALL Abend2('I assume all blocks to be distributed among the
     & nodes but there are still some left')
        ENDIF
      END IF
C
C     write EOF-mark
C
      CALL ITODS(-1,1,-1,LUOUT)
C
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE FIND_ACTIVE_BLOCKS_MR(LUIN,LBLK,BLK_A,SEGMNT,IFHBLK)
*
*. Find the active (nonvanishing blocks) on LUIN
*. Non vanishing block is flagged by a 1.0 ( note : real)
*  in BLK_A
*  master version
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
*. Output
      DIMENSION BLK_A(*)
*. Scratch
      DIMENSION SEGMNT(*)
*      INTEGER ISTBLK
*
      CALL REWINE(LUIN,LBLK)
*
      IBLK = 0
*      IBLK = ISTBLK
      NBLK_A = 0
      NBLK = 0
**      write(6,*)'ISTBLK',ISTBLK
*. Loop over blocks
 1000 CONTINUE
        IBLK = IBLK + 1
*        write(6,*)'IBLK =',IBLK
        IF(IBLK.GT.IFHBLK) GOTO 1001
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
          READ(LUIN) LBL
        ELSE IF  (LBLK .LT. 0 ) THEN
          CALL IFRMDS(LBL,1,-1,LUIN)
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
C         CALL FRMDSC(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK)
          IF(IMZERO.EQ.0) THEN
           NBLK_A = NBLK_A + 1
           BLK_A(IBLK) = 1.0D0
          ELSE
           BLK_A(IBLK) = 0.0D0
          END IF
        END IF
*      write(6,*)'lbl of IBLK =',lbl, iblk
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
1001  CONTINUE
      NBLK =  IBLK - 1
*
      NTEST = 0
      IF(NTEST.GE.1) THEN
        WRITE(6,'(A,I8,I8)')
     &  ' FIND_A.... Number of total and active Blocks',NBLK,NBLK_A
      END IF
      IF(NTEST.GE.100) THEN
        WRITE(6,*) ' Active blocks (master)'
        CALL WRTMAT(BLK_A,1,NBLK,1,NBLK)
      END IF
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE ADDSB(A,B,NDIM)
*
* add recieved SB vector to current SB vector
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      DIMENSION A(*),B(*)
*
      DO 100 I = 1,NDIM
        A(I) = A(I) + B(I)
  100 CONTINUE
*
      RETURN
      END
*
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE DIST_NODES(NTOTALND,NACTROOT,JDISTROOT)
*
* find an even distribution for the sum NTOTALND of the nodes 
* w.r.t. the number of roots NACTROOT that have to be calculated
* distribution is stored in array JDISTROOT 
*
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
*
      INTEGER NTOTALND,NACTROOT,ISTNODE,ISTIADD 
      DIMENSION JDISTROOT(NTOTALND)
*
      ISTNODE = 1
      ISTIADD = 1
 1    CONTINUE
      IF(ISTIADD.LE.NACTROOT)THEN
*
        JDISTROOT(ISTNODE) = ISTIADD
        ISTIADD = ISTIADD + 1
        ISTNODE = ISTNODE + 1
        GOTO 1
*        
      ELSE IF(ISTIADD.GT.NACTROOT) THEN
*
        ISTIADD = NACTROOT
        GOTO 2
*
      END IF
 2    CONTINUE
*
      ISTNODE = NTOTALND
 3    CONTINUE
      ISTIADD = NACTROOT
 4    CONTINUE
      IF(ISTIADD.NE.0.AND.ISTNODE.NE.NACTROOT) THEN
        JDISTROOT(ISTNODE) = ISTIADD
        ISTIADD = ISTIADD - 1
        ISTNODE = ISTNODE - 1
        GOTO 4
      ELSE IF(ISTIADD.EQ.0.AND.ISTNODE.GT.NACTROOT) THEN
        GOTO 3 
      ELSE IF(ISTIADD.NE.0.AND.ISTNODE.EQ.NACTROOT) THEN
        GOTO 5
      ELSE IF(ISTIADD.EQ.0.AND.ISTNODE.EQ.NACTROOT) THEN
        GOTO 5
      END IF
 5    CONTINUE
*    
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaption by Timo Fleig                *
* parallelization by Stefan Knecht                                    *
*                                                                     *
***********************************************************************
      SUBROUTINE MSSNDBLCK(LUIN,LUOUT,SEGMNT,JCOMM,LBLK,IDISTLIST,
     &           NMBLOCK)
*
* Read vector (master) from disc file luin to array segmnt
* and send to node IDISTLIST(block) the vector where it is 
* written in file luout
*
* LBLK DEFINES STRUCTURE OF FILE
* COMM IS THE COMMUNICATOR USED FOR MPI-COMMUNICATION
* --> DEFAULT: MPI_COMM_WORLD
* YOUR CALL SHOULD LOOK LIKE: CALL MSS...(LUIN,...,MPI_COMM_WORLD,LBLK)
*
*
* Type of file LUOUT is inherited from LUIN
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*),IDISTLIST(NMBLOCK)
      INTEGER ISENDN,NROUND
C
*
      IREW = 1
*
      IF( IREW .NE. 0 ) THEN
        IF (MYNEW_ID.EQ.0) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (MYNEW_ID.NE.0) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      ISENDN = 0
      NROUND = 0
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
      NROUND = NROUND + 1
      IF(NROUND.LE.NMBLOCK+1) THEN
        IF(NROUND.EQ.NMBLOCK+1) THEN
          GOTO 9000
        ELSE
          ISENDN = IDISTLIST(NROUND)
*          if (myproc.eq.master) WRITE(6,*)'NROUND=',NROUND
*          if (myproc.eq.master) WRITE(6,*)'I..D(NROUND)=ISENDN =',ISENDN
*          if (myproc.eq.master) WRITE(6,*)'LBLK =',LBLK
*          if (myproc.eq.master) WRITE(6,*)'LUIN =',LUIN
        END IF
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (MYNEW_ID.EQ.0.and.ISENDN.GT.0) then
              call mpi_send(LBL,1,my_MPI_INTEGER,ISENDN,60,
     &                      JCOMM,ISTAT,IERR)
            end if
            if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
              call mpi_recv(LBL,1,my_MPI_INTEGER,0,60,
     &                      JCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (MYNEW_ID.EQ.0) then
            READ(LUIN) LBL
          if (MYNEW_ID.EQ.0.and.ISENDN.GT.0) then
            call mpi_send(LBL,1,my_MPI_INTEGER,ISENDN,60,
     &                    JCOMM,ISTAT,IERR)
          end if
          end if


          if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,0,60,
     &                    JCOMM,ISTAT,IERR)
            WRITE(LUOUT) LBL
        end if
        if (MYNEW_ID.EQ.0) then
        write(6,*) 'SR MSSNDBLCK:  LBLK = ', LBLK
        end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (MYNEW_ID.EQ.0) then
*            write(6,*)'master reads in block',NROUND
            CALL IFRMDS(LBL,1,-1,LUIN)
          if (MYNEW_ID.EQ.0.and.ISENDN.GT.0) then
            call mpi_send(LBL,1,my_MPI_INTEGER,ISENDN,60,
     &                    JCOMM,ISTAT,IERR)
*            write(6,*)'master has sent block',NROUND
          end if
          end if
*
*
          if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,0,60,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS(LBL,1,-1,LUOUT)
C?        write(6,*) ' COPVCD LBL : ', LBL
          end if
*
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (MYNEW_ID.EQ.0) then
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
            if (ISENDN.GT.0) then
              call mpi_send(IMZERO,1,my_MPI_INTEGER,ISENDN,61,
     &             JCOMM,ISTAT,IERR)
              call mpi_send(IAMPACK,1,my_MPI_INTEGER,ISENDN,62,
     &             JCOMM,ISTAT,IERR)
              call mpi_send(SEGMNT,LBL,MPI_REAL8,ISENDN,63,
     &             JCOMM,ISTAT,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
             call mpi_recv(IMZERO,1,my_MPI_INTEGER,0,61,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,my_MPI_INTEGER,0,62,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,0,63,
     &                     JCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
                CALL TODSC(SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (MYNEW_ID.NE.0.AND.MYNEW_ID.EQ.ISENDN) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
9000  CONTINUE
      IF (MYNEW_ID.NE.0) THEN
        CALL ITODS(-1,1,-1,LUOUT)
      END IF
*
      RETURN
      END
*
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MSBCVC_MPI(LUIN,LUOUT,SEGMNT,COMM,LBLK)
*
* Read vector (master) from disc file luin to array segmnt
* and broadcast the vector to all nodes
* where vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
*
* LBLK DEFINES STRUCTURE OF FILE
* COMM IS THE COMMUNICATOR USED FOR MPI-COMMUNICATION
* --> DEFAULT: MPI_COMM_WORLD
* YOUR CALL SHOULD LOOK LIKE: CALL MSB...(LUIN,...,MPI_COMM_WORLD,LBLK)
*
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*)
C
*
      IREW = 0
*
      IF( IREW .NE. 0 ) THEN
        IF (LUCI_MYPROC.EQ.LUCI_MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (LUCI_MYPROC.NE.LUCI_MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      ihere = 0
C
C LOOP OVER BLOCKS
C
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
              call mpi_bcast(LBL,1,my_MPI_INTEGER,LUCI_MASTER,
     &                      COMM,IERR)
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (myproc.eq.master) then
            READ(LUIN) LBL
          end if
            call mpi_bcast(LBL,1,my_MPI_INTEGER,LUCI_MASTER,
     &                    COMM,IERR)
           if (myproc.ne.master) then
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
           end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (myproc.eq.master) then
            CALL IFRMDS(LBL,1,-1,LUIN)
          end if
            call mpi_bcast(LBL,1,my_MPI_INTEGER,LUCI_MASTER,
     &                    COMM,IERR)
*
*
          if (myproc.ne.master) then
            CALL ITODS (LBL,1,-1,LUOUT)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (myproc.eq.master) then
            ihere = ihere + 1
**            write(6,*)'i am here for the x time', ihere
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                   NO_ZEROING)
          end if
          call mpi_bcast(IMZERO,1,my_MPI_INTEGER,LUCI_MASTER,
     &                  COMM,IERR)
          call mpi_bcast(IAMPACK,1,my_MPI_INTEGER,LUCI_MASTER,
     &                  COMM,IERR)
          call mpi_bcast(SEGMNT,LBL,MPI_REAL8,LUCI_MASTER,
     &                  COMM,IERR)
*
**          if (myproc.eq.master) then
**            write(6,*)'lbl,IMZERO,IAMPACK = ',lbl,IMZERO,IAMPACK
**            write(6,*)'this is segment(LBL)', SEGMNT(LBL)
**          end if
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (myproc.ne.master) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (myproc.ne.master) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE MSCPCFF(SEGMNT,LUCFF,NJBATCH,JRESTC,JLOOPC,JLOOPCS,
     &                   NYCOMM,JCBLOCK,NNBLOCKS,LBLK)
*
* Read coff. blocks (master) from disc file lucff to array segmnt
* and send to node X
* where vector is written in file lucff
*
* MPI adaption         Stefan Knecht, Aug. 17., 2006
* MSDTAG: send-tag; NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*), JCBLOCK(*), NNBLOCKS(NEWCOMM_PROC)
*.scratch
      INTEGER MTAG,NBLCK,JMBLOCK,JSBLOCK,JNCBLOCK,IGETRNR,ICHNBLCK
C
*
      IREW = 1
      JMBLOCK = 0
      JSBLOCK = 0
      NBLCK = 0
      IGETRNR = 0
      NODEX = 0
      LBLK = -1
*
      IF( IREW .NE. 0 ) THEN
        IF (MYNEW_ID.EQ.0) THEN
          CALL REWINE( LUCFF ,LBLK)
        ELSE IF (MYNEW_ID.NE.0) THEN
          CALL REWINE( LUCFF ,LBLK)
        END IF
      END IF
      MTAG = 100
*
      IF(MYNEW_ID.EQ.0) THEN
*
        DO ISK = 1, JLOOPC
*
           JCBATCH = ISK
           JMBLOCK = JMBLOCK + JCBLOCK(JCBATCH)
*
        END DO
*

*
        IF(NEWCOMM_PROC.LT.NJBATCH) THEN
**
          DO IGTB = 1, NEWCOMM_PROC -1
             JNCBLOCK = 0
             IF(IGTB.GT.(NEWCOMM_PROC-JRESTC)) THEN
               IGETRNR = IABS(JRESTC - NEWCOMM_PROC + IGTB)
               ISKIPR = (IGTB*JLOOPCS) + IGETRNR
               DO ISK = 1, JLOOPC+1
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
               END DO
               NNBLOCKS(IGTB) = JNCBLOCK
             ELSE IF(IGTB.EQ.(NEWCOMM_PROC-JRESTC)) THEN
               ISKIPR = (IGTB*JLOOPCS)
               DO ISK = 1, JLOOPC+1
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
               END DO
               NNBLOCKS(IGTB) = JNCBLOCK
             ELSE IF(IGTB.LT.(NEWCOMM_PROC-JRESTC)) THEN
               ISKIPR = (IGTB*JLOOPCS)
               DO ISK = 1, JLOOPC
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
               END DO
               NNBLOCKS(IGTB) = JNCBLOCK
             END IF
          END DO
**
        ELSE IF(NEWCOMM_PROC.EQ.NJBATCH) THEN
*
          DO IGTB = 1, NEWCOMM_PROC -1
             JNCBLOCK = 0
             ISKIPR = (IGTB*JLOOPCS)
             DO ISK = 1, JLOOPC
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
             END DO
               NNBLOCKS(IGTB) = JNCBLOCK
          END DO
**
        ELSE IF(NEWCOMM_PROC.GT.NJBATCH) THEN
            DO IGTB = 1, (NJBATCH - 1)
            JNCBLOCK = 0
            ISKIPR = (IGTB*JLOOPCS)
            DO ISK = 1, JLOOPC
                  JCBATCH = ISKIPR + ISK
                  JNCBLOCK = JNCBLOCK + JCBLOCK(JCBATCH)
             END DO
               NNBLOCKS(IGTB) = JNCBLOCK
          END DO
        END IF
**
*
      ELSE IF (MYNEW_ID.NE.0) THEN
*
        ISKIPR = 0
* how many records do i have to skip ?
*
        IF(NEWCOMM_PROC.LT.NJBATCH) THEN
*
          IF(MYNEW_ID.GT.(NEWCOMM_PROC-JRESTC)) THEN
*
            IGETRNR = IABS(JRESTC - NEWCOMM_PROC + MYNEW_ID)
            ISKIPR = (MYNEW_ID*JLOOPCS) + IGETRNR
*
          ELSE IF(MYNEW_ID.LE.(NEWCOMM_PROC-JRESTC)) THEN
*
            ISKIPR = (MYNEW_ID*JLOOPCS)
*
          END IF
        ELSE IF(NEWCOMM_PROC.GE.NJBATCH) THEN
*
            ISKIPR = (MYNEW_ID*JLOOPCS)
*
        END IF
*
        DO ISK = 1, JLOOPC
          JCBATCH = ISKIPR + ISK
          JSBLOCK = JSBLOCK + JCBLOCK(JCBATCH)
        END DO
*
      END IF
C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (MYNEW_ID.eq.0.and.NEWCOMM_PROC.gt.1) then
              if(nodex.gt.0) then
                call mpi_send(LBL,1,my_MPI_INTEGER,NODEX,MTAG-2,
     &                        NYCOMM,IERR)
              end if
            end if
            if (MYNEW_ID.ne.0) then
              call mpi_recv(LBL,1,my_MPI_INTEGER,0,MTAG-2,
     &                      NYCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (MYNEW_ID.eq.0) then
            READ(LUCFF) LBL
            if(nodex.gt.0.and.NEWCOMM_PROC.gt.1) then
              call mpi_send(LBL,1,my_MPI_INTEGER,NODEX,MTAG-2,
     &                      NYCOMM,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.ne.0) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,0,MTAG-2,
     &                    NYCOMM,ISTAT,IERR)
            WRITE(LUCFF) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
            end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (MYNEW_ID.eq.0) then
            CALL IFRMDS(LBL,1,-1,LUCFF)
            if(nodex.gt.0.and.NEWCOMM_PROC.gt.1) then
              call mpi_send(LBL,1,my_MPI_INTEGER,NODEX,MTAG-2,
     &                      NYCOMM,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.ne.0) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,LUCI_MASTER,MTAG-2,
     &                    NYCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUCFF)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (MYNEW_ID.eq.0) then
            CALL FRMDSC2(SEGMNT,LBL,KBLK,LUCFF,IMZERO,IAMPACK,
     &                 NO_ZEROING)
            if(nodex.gt.0.and.NEWCOMM_PROC.gt.1) then
              call mpi_send(IMZERO,1,my_MPI_INTEGER,NODEX,MTAG+1,
     &                      NYCOMM,IERR)
              call mpi_send(IAMPACK,1,my_MPI_INTEGER,NODEX,MTAG+2,
     &                      NYCOMM,IERR)
              call mpi_send(SEGMNT,LBL,MPI_REAL8,NODEX,MTAG+3,
     &                      NYCOMM,IERR)
            end if
          end if
*
*
          if (MYNEW_ID.ne.0) then
             call mpi_recv(IMZERO,1,my_MPI_INTEGER,0,MTAG+1,
     &                     NYCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,my_MPI_INTEGER,0,MTAG+2,
     &                     NYCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,0,MTAG+3,
     &                     NYCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (MYNEW_ID.ne.0) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUCFF)
              end if
            ELSE
              if (MYNEW_ID.ne.0) then
                CALL TODSCP(SEGMNT,LBL,KBLK,LUCFF)
              end if
            END IF
          ELSE
            if (MYNEW_ID.ne.0) then
              CALL ZERORC(LBL,LUCFF,IAMPACK)
            end if
          END IF
        END IF
*
      NBLCK = NBLCK + 1
**      if(myproc.eq.master) then
**        write(6,*)'NBLCK (myproc) =', NBLCK,myproc
**      end if
      IF(MYNEW_ID.EQ.0.and.NEWCOMM_PROC.gt.1) THEN
*
        IF(NBLCK.LT.JMBLOCK) THEN
*           write(6,*)'i am here?', myproc
           GOTO 999
        END IF
*
        IF(NBLCK.EQ.JMBLOCK) THEN
*           write(6,*)'start with node1 at NBLCK=',NBLCK
           NODEX = 1
           GOTO 999
        END IF
*
        IF(NODEX.GT.0) THEN
          ICHNBLCK = 0
          ICHNBLCK = JMBLOCK
          DO IFCH = 1, NODEX
             ICHNBLCK = ICHNBLCK + NNBLOCKS(IFCH)
          END DO
**          write(6,*)'ICHNBLCK=',ICHNBLCK
          IF(NBLCK.EQ.ICHNBLCK) THEN
*.? is the complete number of nodes relevant? I assume, it is NOT!
            NODEX = NODEX + 1
*            write(6,*)'master comm. now with nodex',NODEX
*.?
            IF(NEWCOMM_PROC.LE.NJBATCH) THEN
              IF(NODEX.LE.(NEWCOMM_PROC-1)) THEN
**                write(6,*)'master comm. now with nodex',NODEX
                GOTO 999
              ELSE IF(NODEX.GT.(NEWCOMM_PROC-1)) THEN
**                write(6,*)'master comm. NOT with nodex',NODEX
                GOTO 1001
              END IF
*
            ELSE IF(NEWCOMM_PROC.GT.NJBATCH) THEN
              IF(NODEX.LE.(NJBATCH-1)) THEN
**                write(6,*)'master comm. now with nodex',NODEX
                GOTO 999
              ELSE IF(NODEX.GT.(NJBATCH-1)) THEN
**                write(6,*)'master comm. NOT with nodex',NODEX
                GOTO 1001
              END IF
            END IF
          END IF
*.        /\ end if of (NBLCK.EQ.(JMBLOCK + NNBLOCKS(NODEX))
        END IF
      ELSE IF (MYNEW_ID.NE.0) THEN
        IF(NBLCK.EQ.JSBLOCK) THEN
          CALL ITODS (-1,1,-1,LUCFF)
          GOTO 1001
        END IF
      END IF
*
999   CONTINUE
*      write(6,*)'after 999 (myproc)', myproc
*      write(6,*)'LBL, LBLK',LBL,LBLK
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
1001  CONTINUE
**      write(6,*)'everything done', myproc

      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_MS(LUIN,LUOUT,SEGMNT,NOFSLV,JCOMM,LBLK)
*
* Read vector (master) from disc file luin to array segmnt
* and send to node nofslv using communicator jcomm
* where vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
* NOFSLV: process-# of the slave
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*)
C
*
      IREW = 0
      IF(LBLK.NE.-1) THEN
        LBLK = -1
      END IF
*
      IF( IREW .NE. 0 ) THEN
        IF (LUCI_MYPROC.EQ.LUCI_MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (LUCI_MYPROC.NE.LUCI_MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      MYNSLV = NOFSLV

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (myproc.eq.master) then
              call mpi_send(LBL,1,my_MPI_INTEGER,MYNSLV,80,
     &                      JCOMM,IERR)
            end if
            if (myproc.ne.master) then
              call mpi_recv(LBL,1,my_MPI_INTEGER,LUCI_MASTER,80,
     &                      JCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (myproc.eq.master) then
            READ(LUIN) LBL
            call mpi_send(LBL,1,my_MPI_INTEGER,MYNSLV,80,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.ne.master) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,LUCI_MASTER,80,
     &                    JCOMM,ISTAT,IERR)
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
            end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (myproc.eq.master) then
            CALL IFRMDS(LBL,1,-1,LUIN)
            call mpi_send(LBL,1,my_MPI_INTEGER,MYNSLV,80,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.ne.master) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,LUCI_MASTER,80,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUOUT)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (myproc.eq.master) then
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          call mpi_send(IMZERO,1,my_MPI_INTEGER,MYNSLV,81,
     &                  JCOMM,IERR)
          call mpi_send(IAMPACK,1,my_MPI_INTEGER,MYNSLV,82,
     &                  JCOMM,IERR)
          call mpi_send(SEGMNT,LBL,MPI_REAL8,MYNSLV,83,
     &                  JCOMM,IERR)
          end if
*
*
          if (myproc.ne.master) then
             call mpi_recv(IMZERO,1,my_MPI_INTEGER,LUCI_MASTER,81,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,my_MPI_INTEGER,LUCI_MASTER,82,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,LUCI_MASTER,83,
     &                     JCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (myproc.ne.master) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (myproc.ne.master) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_SPECIAL(LUIN,LUOUT,SEGMNT,LBLK,ACT_BLK,NUMBLK)
*
* Read vector from disc file luin to array segmnt
* and send to master using communicator jcomm.
* vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
* NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*),ACT_BLK(*)
C
*
      IREW = 0
      IBLOCK = 0
*
      IF( IREW .NE. 0 ) THEN
          CALL REWINE( LUIN ,LBLK)
          CALL REWINE( LUOUT ,LBLK)
      END IF

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IBLOCK = IBLOCK + 1
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
        ELSE IF ( LBLK .EQ. 0 ) THEN
            READ(LUIN) LBL
*
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
*
        ELSE IF  (LBLK .LT. 0 ) THEN
            CALL IFRMDS(LBL,1,-1,LUIN)
*
            CALL ITODS (LBL,1,-1,LUOUT)
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          IF(IMZERO.EQ.0) THEN
            ACT_BLK(IBLOCK) = 1.0D0 
            IF(IAMPACK.EQ.0) THEN
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
              ACT_BLK(IBLOCK) = 0.0D0
              CALL ZERORC(LBL,LUOUT,IAMPACK)
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
*
      NUMBLK = IBLOCK - 1
*
      END
***********************************************************************
*                                                                     *
* LUCITA, by Jeppe Olsen, DIRAC adaptation by Timo Fleig              *
*                                                                     *
***********************************************************************
      SUBROUTINE COPVCD_PAR_SM(LUIN,LUOUT,SEGMNT,NOFSLV,JCOMM,LBLK)
*
* Read vector (node # mynslv) from disc file luin to array segmnt
* and send to master using communicator jcomm.
* vector is written in file luout
*
* MPI adaption         Stefan Knecht, Feb 21, 2006
* NOFSLV: process-# of the node
*
* LBLK DEFINES STRUCTURE OF FILE
*
* Type of file LUOUT is inherited from LUIN
      IMPLICIT REAL*8(A-H,O-Z)
#include "infpar.h"
#include "mpif.h"
      DIMENSION ISTAT(MPI_STATUS_SIZE)
#include "parluci.h"
      DIMENSION SEGMNT(*)
C
*
      IREW = 0
*
      IF( IREW .NE. 0 ) THEN
        IF (LUCI_MYPROC.NE.LUCI_MASTER) THEN
          CALL REWINE( LUIN ,LBLK)
        ELSE IF (LUCI_MYPROC.EQ.LUCI_MASTER) THEN
          CALL REWINE( LUOUT ,LBLK)
        END IF
      END IF
      MYNSLV = NOFSLV

C
C LOOP OVER BLOCKS
C
C?      write(6,*) ' COPVCD LBLK : ', LBLK
 1000 CONTINUE
        IF(LBLK .GT. 0 ) THEN
          LBL = LBLK
            if (myproc.ne.master) then
              call mpi_send(LBL,1,my_MPI_INTEGER,LUCI_MASTER,90,
     &                      JCOMM,IERR)
            end if
            if (myproc.eq.master) then
              call mpi_recv(LBL,1,my_MPI_INTEGER,MYNSLV,90,
     &                      JCOMM,ISTAT,IERR)
            end if
        ELSE IF ( LBLK .EQ. 0 ) THEN
          if (myproc.ne.master) then
            READ(LUIN) LBL
            call mpi_send(LBL,1,my_MPI_INTEGER,LUCI_MASTER,90,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.eq.master) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,MYNSLV,90,
     &                    JCOMM,ISTAT,IERR)
            WRITE(LUOUT) LBL
C?        write(6,*) ' COPVCD LBL : ', LBL
            end if
*
        ELSE IF  (LBLK .LT. 0 ) THEN
          if (myproc.ne.master) then
            CALL IFRMDS(LBL,1,-1,LUIN)
            call mpi_send(LBL,1,my_MPI_INTEGER,LUCI_MASTER,90,
     &                    JCOMM,IERR)
          end if
*
*
          if (myproc.eq.master) then
            call mpi_recv(LBL,1,my_MPI_INTEGER,MYNSLV,90,
     &                    JCOMM,ISTAT,IERR)
            CALL ITODS (LBL,1,-1,LUOUT)
          end if
        END IF
*
        IF( LBL .GE. 0 ) THEN
          IF(LBLK .GE.0 ) THEN
            KBLK = LBL
          ELSE
            KBLK = -1
          END IF
C?        write(6,*) ' LBL and KBLK ', LBL,KBLK
          NO_ZEROING = 1
*
          if (myproc.ne.master) then
          CALL FRMDSC2(SEGMNT,LBL,KBLK,LUIN,IMZERO,IAMPACK,
     &                 NO_ZEROING)
          call mpi_send(IMZERO,1,my_MPI_INTEGER,LUCI_MASTER,91,
     &                  JCOMM,IERR)
          call mpi_send(IAMPACK,1,my_MPI_INTEGER,LUCI_MASTER,92,
     &                  JCOMM,IERR)
          call mpi_send(SEGMNT,LBL,MPI_REAL8,LUCI_MASTER,93,
     &                  JCOMM,IERR)
          end if
*
*
          if (myproc.eq.master) then
             call mpi_recv(IMZERO,1,my_MPI_INTEGER,MYNSLV,91,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(IAMPACK,1,my_MPI_INTEGER,MYNSLV,92,
     &                     JCOMM,ISTAT,IERR)
             call mpi_recv(SEGMNT,LBL,MPI_REAL8,MYNSLV,93,
     &                     JCOMM,ISTAT,IERR)
          end if
*
*
          IF(IMZERO.EQ.0) THEN
            IF(IAMPACK.EQ.0) THEN
              if (myproc.eq.master) then
                CALL TODSC (SEGMNT,LBL,KBLK,LUOUT)
              end if
            ELSE
              CALL TODSCP(SEGMNT,LBL,KBLK,LUOUT)
            END IF
          ELSE
            if (myproc.eq.master) then
              CALL ZERORC(LBL,LUOUT,IAMPACK)
            end if
          END IF
        END IF
      IF( LBL .GE. 0 .AND. LBLK .LE. 0 ) GOTO 1000
C
      RETURN
      END
*****************************************************
* dummy routine for normal compilation
       SUBROUTINE DO_NOTHING(NOTHING)
       IMPLICIT REAL*8(A-H,O-Z)
*        write(6,*) 'let us do nothing',nothing
       END
*****************************************************
#else
* dummy routine for normal compilation
       SUBROUTINE DO_NOTHING2(NOTHING)
       IMPLICIT REAL*8(A-H,O-Z)
*        write(6,*) 'let us do nothing',nothing  
       END 
#endif 
