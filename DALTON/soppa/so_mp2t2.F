C  /* Deck so_mp2t2 */
      SUBROUTINE SO_MP2T2(FOCKD,LFOCKD,T2AM,LT2AMH)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, November 1995
C     Stephan P. A. Sauer: 10.11.2003: merge with Dalton 2.0
C
C     PURPOSE: Gets (ia|jb) integrals in T2AM and
C              calculates the T² MP2 amplitudes.
C
C     Written by Lilli Irene Ør Kristensen
C
#include "implicit.h"
#include "priunit.h"
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION FOCKD(LFOCKD),T2AM(LT2AMH)
C
#include "ccorb.h"
#include "ccsdinp.h"
#include "ccsdsym.h"
#include "soppinf.h"
C
C------------------------------
C     Statement function INDEX.
C------------------------------
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_MP2T2')
C
C
C-----------------------------------
C    Creating NAIBJ
C-----------------------------------
C
      DO 100 ISYMBJ = 1,NSYM
         ISYMAI = ISYMBJ !At his piont AI and BJ is the same
         DO 110 ISYMJ = 1,NSYM
            ISYMB = MULD2H(ISYMJ,ISYMBJ)
            DO 120 ISYMI = 1,NSYM
               ISYMA = MULD2H(ISYMI,ISYMAI)
               DO 130 J = 1,NRHF(ISYMJ)
                  KOFFJ = IRHF(ISYMJ) + J  !Position for KOFFJ
                  DO 140 B = 1, NVIR(ISYMB)
                     KOFFB = IVIR(ISYMB) + B
                     NBJ = IT1AM(ISYMB,ISYMJ)+ NVIR(ISYMB)*(J-1)+B
C                     NBJ is a loop over 1,NT1AM(ISYMBJ), but placed right
                     DO 150 I = 1,NRHF(ISYMI)
                        KOFFI = IRHF(ISYMI) + I
                        DO 160 A = 1,NVIR(ISYMA)
                           KOFFA = IVIR(ISYMA) + A
                           NAI = IT1AM(ISYMA,ISYMI)+NVIR(ISYMA)*(I-1)+A
C
                           IF (NAI .GT. NBJ) GOTO 160
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)+INDEX(NAI,NBJ)
C
                           WRITE(LUPRI,*) 'A',A,'I',I,'B',B,'J',J
                           WRITE(LUPRI,*) 'NAI',NAI,'NBJ',NBJ
                           WRITE(LUPRI,*) 'NAIBJ',NAIBJ
                           T2AM(NAIBJ) = T2AM(NAIBJ) /
     *                                     (FOCKD(KOFFI) + FOCKD(KOFFJ)
     *                                    - FOCKD(KOFFA) - FOCKD(KOFFB))
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C Printing the amplitudes
      IF (IPRINT .GT. 10) THEN
         CALL AROUND('MP2 amplitudes')
         DO 250 ISYMBJ = 1,NSYM
            ISYMAI = ISYMBJ
            KOFF   = IT2AM(ISYMAI,ISYMBJ) + 1
            NTOTAI = NT1AM(ISYMAI)
            CALL OUTPUT(T2AM(KOFF),1,NTOTAI,1,1,NTOTAI,1,1,LUPRI)
250      CONTINUE
      END IF
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL FLSHFO(LUPRI)
C
      CALL QEXIT('SO_MP2T2')
C
      RETURN
      END
