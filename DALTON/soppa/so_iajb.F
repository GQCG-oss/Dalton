C  /* Deck so_iajb */
      SUBROUTINE SO_IAJB(XAIBJ,WORK,LWORK)
C
C     The Original subroutine CCSD_IAJB is written by Henrik Koch 27-Mar-1990.
C
C     This Routine has been modified such that it can calculate and
C      transfrom the intregrals in the SOPPA program.
C
C     It still contains the MO integrals needed for gradients and frozen
C       core FOP.
C     NB: Frozen core does NOT work when doing SOPPA calculations
C           without the CC input, i.e. not using the CC program.
C
C     Writen by Lilli Irene Ør Kristensen Januar 2017
C
#include "implicit.h"
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0)
#include "priunit.h"
#include "dummy.h"
#include "maxorb.h"
#include "maxash.h"
#include "mxcent.h"
#include "aovec.h"
#include "iratdef.h"
#include "ccorb.h"
#include "ccisao.h"
#include "blocks.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "distcl.h"
#include "eritap.h"
#include "ccfro.h"
#include "ccfop.h"
      DIMENSION XAIBJ(*),WORK(LWORK),INDEXA(MXCORB)
      INTEGER KEND0
C
      LOGICAL OPENED
C
      CHARACTER*8 FNTOC
C
      CALL QENTER('SO_IAJB')
C
C-------------------------------------------------------
C     Initialize the XAIBJ integral array.
C-------------------------------------------------------
C
      CALL DZERO(XAIBJ,NT2AM(ISYMOP))
C
C-------------------------------------------------------
C     Check whether SIRIFC is still open and close it.
C-------------------------------------------------------
C
      CALL GPINQ('SIRIFC','OPENE',OPENED)
      IF (OPENED) THEN
         INQUIRE (FILE='SIRIFC',NUMBER=LUSIFC)
         WRITE(LUPRI,'(2A,I3)') ' SO_IAJB: file SIRIFC is already ',
     &                          'opened with unit number ',LUSIFC
         CALL GPCLOSE (LUSIFC,'KEEP')
         WRITE(LUPRI,'(2A,I3)') ' RP_IAJB: file SIRIFC was closed again'
      ENDIF
C
C-------------------------------------------------------
C     Dynamic allocation of space.
C-------------------------------------------------------
C
      KCMO   = 1
      KEND1 = KCMO   + NLAMDS
C
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space in SO_IAJB')
      ENDIF
C
C---------------------------------------------------------------------
C     Initialize CMO vector:
C---------------------------------------------------------------------
C
      CALL SO_GETMO(WORK(KCMO),NLAMDS,WORK(KEND1),LWRK1)
C
C====================================================
C     Start the loop over distributions of integrals.
C====================================================
C
      IF (DEBUG) THEN
C        IPRERI = 5
         WRITE(LUPRI,'(1X,A,I10)') 'LWORK = ',LWORK
      END IF
C
      IF (DIRECT) THEN
         DTIME  = SECOND()
         IF (HERDIR) THEN
            CALL HERDI1(WORK(KEND1),LWRK1,IPRERI)
         ELSE
            KCCFB1 = KEND1
            KINDXB = KCCFB1 + MXPRIM*MXCONT
            KEND1  = KINDXB + (8*MXSHEL*MXCONT + 1)/IRAT
            LWRK1  = LWORK  - KEND1
            CALL ERIDI1(KODCL1,KODCL2,KODBC1,KODBC2,KRDBC1,KRDBC2,
     &                  KODPP1,KODPP2,KRDPP1,KRDPP2,
     &                  KFREE,LFREE,KEND1,WORK(KCCFB1),WORK(KINDXB),
     &                  WORK(KEND1),LWRK1,IPRERI)
            KEND1 = KFREE
            LWRK1 = LFREE
         ENDIF
         NTOSYM = 1
      ELSE
         NTOSYM = NSYM
      ENDIF
C
      THRDIS = 1.0D-8
      ICOUNT1 = 0
      ICOUNT2 = 0
C
      KENDSV = KEND1
      LWRKSV = LWRK1
C
      DO 100 ISYMD1 = 1,NTOSYM
C
         IF (DIRECT) THEN
            IF (HERDIR) THEN
               NTOT = MAXSHL
            ELSE
               NTOT = MXCALL
            ENDIF
         ELSE
            NTOT = NBAS(ISYMD1)
         ENDIF
C
         DO 110 ILLL = 1,NTOT
C
C-------------------------------------------------------
C           If direct calculate the integrals.
C-------------------------------------------------------
C
            IF (DIRECT) THEN
C
               KEND1 = KENDSV
               LWRK1 = LWRKSV
C
               IF (HERDIR) THEN
                  CALL HERDI2(WORK(KEND1),LWRK1,INDEXA,ILLL,NUMDIS,
     &                        IPRERI)
               ELSE
                  CALL ERIDI2(ILLL,INDEXA,NUMDIS,0,0,
     &                        WORK(KODCL1),WORK(KODCL2),WORK(KODBC1),
     &                        WORK(KODBC2),WORK(KRDBC1),WORK(KRDBC2),
     &                        WORK(KODPP1),WORK(KODPP2),WORK(KRDPP1),
     &                        WORK(KRDPP2),WORK(KCCFB1),WORK(KINDXB),
     &                        WORK(KEND1), LWRK1,IPRERI)
               ENDIF
C
               KRECNR = KEND1
               KEND1  = KRECNR + (NBUFX(0) - 1)/IRAT + 1
               LWRK1  = LWORK  - KEND1
               IF (LWRK1 .LT. 0) THEN
                  CALL QUIT('Insufficient core in CCRHSN')
               END IF
C
            ELSE
               KRECNR = KEND1
               NUMDIS = 1
            ENDIF
C
C-------------------------------------------------------
C           Loop over number of distributions in disk.
C-------------------------------------------------------
C
            DO 120 IDEL2 = 1,NUMDIS
C
               IF (DIRECT) THEN
                  IDEL  = INDEXA(IDEL2)
                  ISYMD = ISAO(IDEL)
               ELSE
                  IDEL  = IBAS(ISYMD1) + ILLL
                  ISYMD = ISYMD1
               ENDIF
C
               ISYMB  = ISYMD
               ISYDIS = MULD2H(ISYMD,ISYMOP)
C
C-------------------------------------------------------
C              Dynamic allocation of work space.
C-------------------------------------------------------
C
               KXINT = KEND1
               KSCR1 = KXINT + NDISAO(ISYDIS)
               KSCR2 = KSCR1 + NBAST*NBAST
               KEND2 = KSCR2 + NBAST*NRHFT
               LWRK2 = LWORK - KEND2
C
               IF (LWRK2 .LT. 0) THEN
                  CALL QUIT('Insufficient space in SO_IAJB')
               ENDIF
C
C-------------------------------------------------------
C              Read in batch of integrals.
C-------------------------------------------------------
C
               IOFFU21 = NDISAO(ISYDIS)
               CALL DZERO(WORK(KXINT),2*NDISAO(ISYDIS))
               CALL CCRDAO(WORK(KXINT),IDEL,IDEL2,WORK(KEND2),LWRK2,
     *                     WORK(KRECNR),DIRECT)
C
C-------------------------------------------------------
C              compute the AO-Fock matrix:
C-------------------------------------------------------
C
C              CALL CC_AOFOCK(WORK(KXINT),WORK(KDNSHF),WORK(KFCKHF),
C     *                        WORK(KEND2),LWRK2,IDEL,ISYMD,.FALSE.,
C     *                        DUMMY,1)
C
C--------------------------------------------------------------------
C        Transform (alpha beta|gamma delta) to (alpha beta| j delta).
C--------------------------------------------------------------------
C
               DO 130 ISYMG = 1,NSYM
C
                  ISYMAB = MULD2H(ISYMG,ISYDIS)
                  ISYMJ  = ISYMG
                  ISYMBJ = MULD2H(ISYMB,ISYMJ)
                  ISYMAI = MULD2H(ISYMBJ,ISYMOP)
C
                  IF (ISYMAI .GT. ISYMBJ) GOTO 130
C
                  KOFF1 = KXINT  + IDSAOG(ISYMG,ISYDIS)
                  KOFF2 = KCMO + ILMRHF(ISYMJ)
                  KOFF6 = KCMO + ILMRHF(ISYMJ)
C
                  IF (LWRK2 .LT. 2*NNBST(ISYMAB)*NRHF(ISYMJ)) THEN
                     CALL QUIT('Insufficient core in SO_IAJB')
                  ENDIF
C
C--------------------------------------------------------
C                 Analyse size of integral distributions.
C--------------------------------------------------------
C
                  DO 140 G = 1,NBAS(ISYMG)
C
                     KOFFG = KXINT + IDSAOG(ISYMG,ISYDIS)
     *                             + NNBST(ISYMAB)*(G - 1)
                     NAB   = NNBST(ISYMAB)
C
                     DO 150 IAB = 1,NAB
  150                CONTINUE
C
                     ICOUNT1 = ICOUNT1 + 1
C
  158                CONTINUE
C
                     ICOUNT2 = ICOUNT2 + 1
C
  140             CONTINUE
C
C-------------------------------------------------------------------
C                 Transform the gamma index in the integral (AB|GD).
C-------------------------------------------------------------------
C
                  NNBSAB = MAX(NNBST(ISYMAB),1)
                  NBASG  = MAX(NBAS(ISYMG),1)
                  CALL DGEMM('N','N',NNBST(ISYMAB),NRHF(ISYMJ),
     *                       NBAS(ISYMG),ONE,WORK(KOFF1),NNBSAB,
     *                       WORK(KOFF2),NBASG,ZERO,WORK(KEND2),
     *                       NNBSAB)
C      (alpha beta|j delta)=(alpha beta| gamma delta)*Lamda(gamma j).
C
C--------------------------------------------------------------------
C                 Transform integrals and add to the result vector.
C--------------------------------------------------------------------
C
                  KOFF4  = IT2AM(ISYMAI,ISYMBJ) + 1
C
                  CALL SO_AIBJ2(WORK(KEND2),XAIBJ(KOFF4),WORK(KCMO),
     *                            WORK(KCMO),WORK(KSCR1),WORK(KSCR2),
     *                            IDEL,ISYMD,ISYMJ,ISYMAB,LUTOC,FNTOC,
     *                            .FALSE.)
C
  130          CONTINUE
C
  120       CONTINUE
C
  110    CONTINUE
C
  100 CONTINUE
C
      KEND1 = KENDSV
      LWRK1 = LWRKSV
C
      IF (IPRINT .GT. 40) THEN
        CALL AROUND('(ia|jb) integral vector')
        DO 252 ISYMBJ = 1,NSYM
             ISYMAI = ISYMBJ
             KOFF   = IT2AM(ISYMAI,ISYMBJ) + 1
             NTOTAI = NT1AM(ISYMAI)
             CALL OUTPUT(XAIBJ(KOFF),1,NTOTAI,1,1,NTOTAI,1,1,LUPRI)
  252   CONTINUE
      END IF
C
      CALL QEXIT('SO_IAJB')

      RETURN
C
      END
