C
C  /* Deck so_rspdrv */
      SUBROUTINE SO_RSPDRV(WORK,LWORK)
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Rasmus Faber, 2017
C
C     PURPOSE: A response wrapper for the SO_RSPLRS AOSOPPA linear
C     response driver
C
      use so_info, only: sop_dp
      implicit none
C From infrsp we need
C FREQ <- Array of requecies
C NFREQ <- Number of frequencies
#include "infrsp.h"
C inflr.h below requires MAXLBL from rspprp.h
#include "rspprp.h"
C We need
C LBLLR <- Label of operators
C NGPLR <- Number of operators of each symmetry
C THCLR <- Convergence threshold
#include "inflr.h"
C  nsym
#include "inforb.h"
      REAL(SOP_DP), INTENT(INOUT) :: WORK(LWORK)
      INTEGER, INTENT(IN) :: LWORK

      character(len=9), parameter :: myname = 'SO_RSPDRV'
      integer :: labcount, klab, kend, lwrk, laboffset(8)

      CALL QENTER(MYNAME)
C
C  Copy labels from response
      labcount = SUM(NGPLR(1:NSYM))
      klab = 1
      kend = klab + labcount
      lwrk = lwork - kend
C
      CALL SO_RSP_LABELCOPY(WORK(KLAB), LABOFFSET, LABCOUNT)
C
C  Call AO-SOPPA Driver routine
      CALL SO_RSPLRS(WORK(KLAB), NGPLR, LABOFFSET,
     &               TRPLET, THCLR, FREQ, NFREQ, 0,
     &               WORK(kend), LWRK)
      CALL QEXIT(MYNAME)

      END SUBROUTINE

      SUBROUTINE SO_RSP_LABELCOPY(SOPPA_LABELS,
     &                            LABEL_OFFSETS,
     &                            LABCOUNT)
C
C  Translate the rsp labels on LBLLR to the list used by aosoppa
C
C
      IMPLICIT NONE
C  nsym
#include "inforb.h"
C inflr.h below requires MAXLBL from rspprp.h
#include "rspprp.h"
C We need
C LBLLR <- Label of operators
C NGPLR <- Number of operators of each symmetry
#include "inflr.h"
      CHARACTER(LEN=8), INTENT(OUT) :: SOPPA_LABELS(LABCOUNT)
      INTEGER,INTENT(OUT) :: LABEL_OFFSETS(8)
      INTEGER,INTENT(IN) :: LABCOUNT

      INTEGER :: ISYM, IOPER, ICOUNT

      ICOUNT = 0
      DO ISYM = 1, NSYM
         LABEL_OFFSETS(ISYM) = ICOUNT
         DO IOPER = 1, NGPLR(ISYM)
            ICOUNT = ICOUNT + 1
            SOPPA_LABELS(ICOUNT) = LBLLR(ISYM,IOPER)
         END DO
      END DO

      END SUBROUTINE

