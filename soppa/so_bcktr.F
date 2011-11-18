#include <single.h>
C
C  /* Deck so_bcktr */
      SUBROUTINE SO_BCKTR(TR1E,LTR1E,TR1D,LTR1D,BTR1E,LBTR1E,BTR1D,
     &                    LBTR1D,BTJ1E,LBTJ1E,BTJ1D,LBTJ1D,CMO,
     &                    LCMO,ISYMTR)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Henrik Koch, Stephan Sauer and Keld Bak, September 1995
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate one-index backtransformed trial-vectors
C
#include <implicit.h>
#include <priunit.h>
C
      PARAMETER (ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION TR1E(LTR1E),   TR1D(LTR1D)
      DIMENSION BTR1E(LBTR1E), BTR1D(LBTR1D), BTJ1E(LBTJ1E)
      DIMENSION BTJ1D(LBTJ1D), CMO(LCMO)
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_BCKTR')
C
C------------------------------------------------------------------
C     Calculate trialvectors with the virtual index backtransformed
C     to an ao index. The results are written in BTR1E and BTR1D.
C------------------------------------------------------------------
C
      DO 100 ISYMI = 1,NSYM
C
         ISYMA  = MULD2H(ISYMI,ISYMTR)
C
         NTOTA  = MAX(NVIR(ISYMA),1)
         NTOTAL = MAX(NBAS(ISYMA),1)
C
         KOFF1  = ILMVIR(ISYMA) + 1
         KOFF2  = IT1AM(ISYMA,ISYMI) + 1
         KOFF3  = IT1AO(ISYMA,ISYMI) + 1
C
         CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NVIR(ISYMA),ONE,
     &              CMO(KOFF1),NTOTAL,TR1E(KOFF2),NTOTA,ZERO,
     &              BTR1E(KOFF3),NTOTAL)
C
         CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NVIR(ISYMA),ONE,
     &              CMO(KOFF1),NTOTAL,TR1D(KOFF2),NTOTA,ZERO,
     &              BTR1D(KOFF3),NTOTAL)
C
  100 CONTINUE
C
C-------------------------------------------------------------------
C     Calculate trialvectors with the occupied index backtransformed
C     to an ao index. The results are written in BTJ1E and BTJ1D.
C-------------------------------------------------------------------
C
      DO 200 ISYMA = 1,NSYM
C
         ISYMI  = MULD2H(ISYMA,ISYMTR)
C
         NTOTI  = MAX(NBAS(ISYMI),1)
         NTOTA  = MAX(NVIR(ISYMA),1)
C
         KOFF1  = ILMRHF(ISYMI) + 1
         KOFF2  = IT1AM(ISYMA,ISYMI) + 1
         KOFF3  = IMATAV(ISYMI,ISYMA) + 1
C
         CALL DGEMM('N','T',NBAS(ISYMI),NVIR(ISYMA),NRHF(ISYMI),ONE,
     &              CMO(KOFF1),NTOTI,TR1E(KOFF2),NTOTA,ZERO,
     &              BTJ1E(KOFF3),NTOTI)
C
         CALL DGEMM('N','T',NBAS(ISYMI),NVIR(ISYMA),NRHF(ISYMI),ONE,
     &              CMO(KOFF1),NTOTI,TR1D(KOFF2),NTOTA,ZERO,
     &              BTJ1D(KOFF3),NTOTI)
C
  200 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_BCKTR')
C
      RETURN
      END
