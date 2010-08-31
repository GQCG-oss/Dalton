C
C  /* Deck so_tmltr */
      SUBROUTINE SO_TMLTR(T2AM,SCAL,ISYOPE)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Andrea Ligabue, January 2004
C
C     Develepode starting frm CCSD_TCMEPKX 
C     Henrik Koch and Alfredo Sanchez.                Dec 1994
C     Made workable for non-symmetric T2AM, Keld Bak, Dec 1996
C
C     Purpose: calculate the left T_aibj starting from the right ones 
C     (I need to call that subroutine with SCAL = HALF)
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0, THREE = 3.0D0)
      PARAMETER (FOUR = 4.0D0)
C
      DIMENSION T2AM(*)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J) - 3)/2 + I + J
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_TMLTR')
C
      FAC = TWO/THREE
C
      DO 100 ISYMIJ = 1,NSYM
C
         ISYMAB = MULD2H(ISYMIJ,ISYOPE)
C
         DO 110 ISYMJ = 1,NSYM
C
            ISYMI = MULD2H(ISYMJ,ISYMIJ)
C
            IF (ISYMI .GT. ISYMJ) GOTO 110
C
            DO 120 ISYMB = 1,NSYM
C
               ISYMA = MULD2H(ISYMB,ISYMAB)
C
               IF (ISYMA .GT. ISYMB) GOTO 120
C
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMBJ = MULD2H(ISYMB,ISYMJ)
               ISYMBI = MULD2H(ISYMB,ISYMI)
               ISYMAJ = MULD2H(ISYMA,ISYMJ)
C
               DO 130 J = 1,NRHF(ISYMJ)
C
                  IF (ISYMI .EQ. ISYMJ) THEN
                     NRHFI =  J
                  ELSE
                     NRHFI = NRHF(ISYMI)
                  ENDIF
C
               IF ( ISYMAI .EQ. ISYMBJ ) THEN
C
                  DO 140 I = 1,NRHFI
C
                     DO 150 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 160 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + INDEX(NAI,NBJ)
C
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + INDEX(NAJ,NBI)
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ)+T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI)+T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  160                   CONTINUE
  150                CONTINUE
  140             CONTINUE
C
               ELSE IF ((ISYMAI.LT.ISYMBJ).AND.(ISYMAJ.LT.ISYMBI)) THEN
C
                  DO 240 I = 1,NRHFI
C
                     DO 250 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 260 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI) * (NBJ - 1) + NAI
C
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ) * (NBI - 1) + NAJ
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ)+T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI)+T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  260                   CONTINUE
  250                CONTINUE
  240             CONTINUE
C
               ELSE IF ((ISYMAI.GT.ISYMBJ).AND.(ISYMAJ.GT.ISYMBI)) THEN
C
                  DO 340 I = 1,NRHFI
C
                     DO 350 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 360 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                           + NT1AM(ISYMBJ) * (NAI - 1) + NBJ
C
                           NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                           + NT1AM(ISYMBI) * (NAJ - 1) + NBI
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ) + T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI) + T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  360                   CONTINUE
  350                CONTINUE
  340             CONTINUE
C
               ELSE IF ((ISYMAI.LT.ISYMBJ).AND.(ISYMAJ.GT.ISYMBI)) THEN
C
                  DO 440 I = 1,NRHFI
C
                     DO 450 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 460 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMAI,ISYMBJ)
     *                           + NT1AM(ISYMAI) * (NBJ - 1) + NAI
C
                           NAJBI = IT2AM(ISYMBI,ISYMAJ)
     *                           + NT1AM(ISYMBI) * (NAJ - 1) + NBI
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ) + T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI) + T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  460                   CONTINUE
  450                CONTINUE
  440             CONTINUE
C
               ELSE IF ((ISYMAI.GT.ISYMBJ).AND.(ISYMAJ.LT.ISYMBI)) THEN
C
                  DO 540 I = 1,NRHFI
C
                     DO 550 B = 1,NVIR(ISYMB)
C
                        IF (ISYMB .EQ. ISYMA) THEN
                           NVIRA = B
                        ELSE
                           NVIRA = NVIR(ISYMA)
                        ENDIF
C
                        NBI = IT1AM(ISYMB,ISYMI)
     *                      + NVIR(ISYMB)*(I - 1) + B
                        NBJ = IT1AM(ISYMB,ISYMJ)
     *                      + NVIR(ISYMB)*(J - 1) + B
C
                        DO 560 A = 1,NVIRA
C
                           NAI = IT1AM(ISYMA,ISYMI)
     *                         + NVIR(ISYMA)*(I - 1) + A
                           NAJ = IT1AM(ISYMA,ISYMJ)
     *                         + NVIR(ISYMA)*(J - 1) + A
C
                           NAIBJ = IT2AM(ISYMBJ,ISYMAI)
     *                           + NT1AM(ISYMBJ) * (NAI - 1) + NBJ
C
                           NAJBI = IT2AM(ISYMAJ,ISYMBI)
     *                           + NT1AM(ISYMAJ) * (NBI - 1) + NAJ
C
                           XAIBJ = FAC*(TWO*T2AM(NAIBJ) + T2AM(NAJBI))
                           XAJBI = FAC*(TWO*T2AM(NAJBI) + T2AM(NAIBJ))
C
                           T2AM(NAIBJ) = XAIBJ
                           T2AM(NAJBI) = XAJBI
C
  560                   CONTINUE
  550                CONTINUE
  540             CONTINUE
C
               END IF
C
  130          CONTINUE
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C---------------------------------------
C     Scale diagonal elements of result.
C---------------------------------------
C
      IF (ISYOPE .NE. 1) GOTO 1000
C
      DO 600 ISYMAI = 1,NSYM
         DO 610 NAI = 1,NT1AM(ISYMAI)
            NAIAI = IT2AM(ISYMAI,ISYMAI) + INDEX(NAI,NAI)
            T2AM(NAIAI) = SCAL*T2AM(NAIAI)
  610    CONTINUE
  600 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
 1000 CALL QEXIT('SO_TMLTR')
C
      RETURN
      END
