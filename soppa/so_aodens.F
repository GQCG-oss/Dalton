#include <single.h>
C
C  /* Deck so_aodens */
      SUBROUTINE SO_AODENS(DENS,LDENS,CMO,LCMO,TR1E,LTR1E,
     &                     TR1D,LTR1D,ISYMTR,WORK,LWORK)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, December 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate RPA AO-density matrix DENS to be used for
C              SOPPA.
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION DENS(LDENS)
      DIMENSION CMO(LCMO),   TR1E(LTR1E),    TR1D(LTR1D)
      DIMENSION WORK(LWORK)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_AODENS')
C
C-------------------------------
C     Initialize density matrix.
C-------------------------------
C
      CALL DZERO(DENS,LDENS)
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA = MULD2H(ISYMI,ISYMTR)
C
         LSCR  = NRHF(ISYMI)*NBAS(ISYMA)
C
         KSCR    = 1
         KEND    = KSCR  + LSCR
         LWORK   = LWORK - KEND
C
         CALL SO_MEMMAX ('SO_AODENS',LWORK1)
         IF (LWORK .LT. 0) CALL STOPIT('SO_AODENS',' ',KEND,LWORK)
C
C-----------------------------------------
C        Calculate density matrix DI1. See
C        Eq. (26) in C.P. 172, 13 (1993).
C-----------------------------------------
C
         ISYMAL = ISYMI
         ISYMBE = ISYMA
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
         NTOTI  = MAX(NRHF(ISYMI),1)
C
         KOFF1  = IT1AM(ISYMA,ISYMI)    + 1
         KOFF2  = ILMVIR(ISYMA)         + 1
         KOFF3  = ILMRHF(ISYMI)         + 1
         KOFF4  = IAODIS(ISYMAL,ISYMBE) + 1
C
         CALL DGEMM('T','T',NRHF(ISYMI),NBAS(ISYMBE),
     &              NVIR(ISYMA),ONE,TR1E(KOFF1),NTOTA,
     &              CMO(KOFF2),NTOTBE,ZERO,WORK(KSCR),NTOTI)
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NBAS(ISYMBE),
     &              NRHF(ISYMI),ONE,CMO(KOFF3),NTOTAL,
     &              WORK(KSCR),NTOTI,ONE,DENS(KOFF4),NTOTAL)
C
C-----------------------------------------------------------
C        Calculate density matrix DI2 and subtract from DI1.
C        See Eqs. (27) and (25) in C.P. 172, 13 (1993).
C-----------------------------------------------------------
C
         ISYMAL = ISYMA
         ISYMBE = ISYMI
C
         NTOTAL = MAX(NBAS(ISYMAL),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTBE = MAX(NBAS(ISYMBE),1)
C
         KOFF1  = ILMVIR(ISYMA)         + 1
         KOFF2  = IT1AM(ISYMA,ISYMI)    + 1
         KOFF3  = ILMRHF(ISYMI)         + 1
         KOFF4  = IAODIS(ISYMAL,ISYMBE) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),
     &              NVIR(ISYMA),ONE,CMO(KOFF1),NTOTAL,
     &              TR1D(KOFF2),NTOTA,ZERO,WORK(KSCR),NTOTAL)
C
         CALL DGEMM('N','T',NBAS(ISYMAL),NBAS(ISYMBE),
     &              NRHF(ISYMI),-ONE,WORK(KSCR),NTOTAL,
     &              CMO(KOFF3),NTOTBE,ONE,DENS(KOFF4),NTOTAL)
C
  100 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_AODENS')
C
      RETURN
      END
