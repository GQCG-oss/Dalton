C
C  /* Deck so_init */
      SUBROUTINE SO_INIT
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Set the pointer-arrays IIJDEN, IABDEN, and IAIDEN
C              and the dimensions NIJDEN, NABDEN, and NAIDEN and
C              NIJDET, NABDET, and NAIDET.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <soppinf.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_INIT')
C
      CALL DZERO(SOTIME,LSOTIM)
      CALL DZERO(SOWTIM,LSOWTI)
      CALL DZERO(SOORWC,LORWCI)
C
      DO 100 ISYMIJ = 1,NSYM
C
         ICOUN1 = 0
         ICOUN2 = 0
         ICOUN3 = 0
C
         DO 200 ISYMI = 1,NSYM
C
            ISYMJ = MULD2H(ISYMI,ISYMIJ)
C
            IIJDEN(ISYMI,ISYMJ) = ICOUN1
            IABDEN(ISYMI,ISYMJ) = ICOUN2
            IAIDEN(ISYMI,ISYMJ) = ICOUN3
C
            ICOUN1 = ICOUN1 + NRHF(ISYMI) * NRHF(ISYMJ)
            ICOUN2 = ICOUN2 + NVIR(ISYMI) * NVIR(ISYMJ)
            ICOUN3 = ICOUN3 + NVIR(ISYMI) * NRHF(ISYMJ)
C
  200    CONTINUE
C
         NIJDEN(ISYMIJ) = ICOUN1
         NABDEN(ISYMIJ) = ICOUN2
         NAIDEN(ISYMIJ) = ICOUN3
C
  100 CONTINUE
C
      NIJDET = 0
      NABDET = 0
      NAIDET = 0
C
      DO 300 ISYMIJ = 1,NSYM
C
         NIJDET = NIJDET + NIJDEN(ISYMIJ)
         NABDET = NABDET + NABDEN(ISYMIJ)
         NAIDET = NAIDET + NAIDEN(ISYMIJ)
C
  300 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_INIT')
C
      RETURN
      END
