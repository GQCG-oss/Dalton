C  /* Deck fijk01 */
      SUBROUTINE FIJK01(D2IJK,CMO,IDEL,ISYMD)
C
C     Written by Asger Halkier 11/6 - 1998
C
C     Purpose: To calculate the first contribution to the
C              CCS 2-electron density including frozen-core
C              contributions - based on DIJK01.
C
#include <implicit.h>
      PARAMETER(FOUR = 4.0D0)
      DIMENSION D2IJK(*), CMO(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccfro.h>
C
      ISYDEN = ISYMD
      ISYMK  = ISYMD
      ISYMIJ = MULD2H(ISYDEN,ISYMK)
C
C------------------------------------------
C     Calculate adresses and add to result.
C------------------------------------------
C
      DO 100 K = 1,NRHFS(ISYMK)
C
         DO 110 ISYMJ = 1,NSYM
C
            DO 120 J = 1,NRHFS(ISYMJ)
C
               NIJK = IFRIJK(ISYMIJ,ISYMK) + NFROIJ(ISYMIJ)*(K - 1)
     *              + IFROIJ(ISYMJ,ISYMJ) + NRHFS(ISYMJ)*(J - 1) + J
               NDEL = ILRHSI(ISYMK) + NBAS(ISYMD)*(K - 1)
     *              + IDEL - IBAS(ISYMD)
C
               D2IJK(NIJK) = D2IJK(NIJK) + FOUR*CMO(NDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck fijk02 */
      SUBROUTINE FIJK02(D2IJK,CMO,IDEL,ISYMD)
C
C     Written by Asger Halkier 11/6 - 1998
C
C     Purpose: To calculate the second contribution to the
C              CCS 2-electron density including frozen-core
C              contributions - based on DIJK02.
C
#include <implicit.h>
      PARAMETER(TWO = 2.0D0)
      DIMENSION D2IJK(*), CMO(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccfro.h>
C
      ISYDEN = ISYMD
      ISYMI  = ISYMD
C
C------------------------------------------
C     Calculate adresses and add to result.
C------------------------------------------
C
      DO 100 ISYMK = 1,NSYM
C
         ISYMIJ = MULD2H(ISYMI,ISYMK)
C
         DO 110 K = 1,NRHFS(ISYMK)
C
            DO 120 I = 1,NRHFS(ISYMI)
C
               NIJK = IFRIJK(ISYMIJ,ISYMK) + NFROIJ(ISYMIJ)*(K - 1)
     *              + IFROIJ(ISYMI,ISYMK) + NRHFS(ISYMI)*(K - 1) + I
               NDEL = ILRHSI(ISYMI) + NBAS(ISYMD)*(I - 1)
     *              + IDEL - IBAS(ISYMD)
C
               D2IJK(NIJK) = D2IJK(NIJK) - TWO*CMO(NDEL)
C
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
C  /* Deck f2_pqao */
      SUBROUTINE F2_PQAO(D2IJG,D2IJK,CMO,ISYDEN)
C
C     Written by Asger Halkier 11/6 - 1998
C
C     Purpose: To calculate the backtransformation of the CCS
C              two-electron density d(ij,k;del) to d(ij,gam;del)
C              including the contributions from frozen core orbitals.
C
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION D2IJG(*), D2IJK(*), CMO(*)
#include <ccorb.h>
#include <ccsdsym.h>
#include <cclr.h>
#include <ccfro.h>
C
C----------------------------------
C     Calculate the transformation.
C----------------------------------
C
         DO 100 ISYMIJ = 1,NSYM
C
            ISYMK  = MULD2H(ISYMIJ,ISYDEN)
            ISYMG  = ISYMK
C
            KOFF1  = IFRIJK(ISYMIJ,ISYMK) + 1
            KOFF2  = ILRHSI(ISYMK) + 1
            KOFF3  = IF2IJG(ISYMIJ,ISYMG) + 1
C
            NTOTIJ = MAX(NFROIJ(ISYMIJ),1)
            NTOTG  = MAX(NBAS(ISYMG),1)
C
            CALL DGEMM('N','T',NFROIJ(ISYMIJ),NBAS(ISYMG),NRHFS(ISYMK),
     *                 ONE,D2IJK(KOFF1),NTOTIJ,CMO(KOFF2),NTOTG,
     *                 ONE,D2IJG(KOFF3),NTOTIJ)
C
  100    CONTINUE
C
      RETURN
      END
