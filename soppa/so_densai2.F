C
C  /* Deck so_densai2 */
      SUBROUTINE SO_DENSAI2(DENSAI,LDENSAI,FOCKD,LFOCKD)
C
C     This routine is part of the atomic integral direct SOPPA program.
C
C     Keld Bak, July 1997
C     Stephan P. A. Sauer, November 2003: merge with DALTON 2.0
C
C     PURPOSE: Calculate the final SOPPA second order one particle
C              density matrix my multipling the sums added into
C              DENSAI in SO_DENSAI1 with inverse MO-energy differences.
C
#include <implicit.h>
#include <priunit.h>
C
#include <ccorb.h>
#include <ccsdinp.h>
#include <ccsdsym.h>
C
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
C
      DIMENSION DENSAI(LDENSAI), FOCKD(LFOCKD)
C
C------------------
C     Add to trace.
C------------------
C
      CALL QENTER('SO_DENSAI2')
C
C---------------------------
C     Divide by E(i) - E(a).
C---------------------------
C
      DO 300 ISYMI = 1,NSYM
C
         ISYMA = ISYMI
C
         DO 200 I = 1,NRHF(ISYMI)
C
            NI    = IRHF(ISYMI) + I
            KOFFA = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1)
C
            DO 100 A = 1,NVIR(ISYMA)      
C
               NA  = IVIR(ISYMA) + A
               NAI = KOFFA       + A
C
               DENSAI(NAI) = DENSAI(NAI) / ( FOCKD(NI) - FOCKD(NA) )
C
  100       CONTINUE
C
  200    CONTINUE
C
  300 CONTINUE
C
C-----------------------
C     Remove from trace.
C-----------------------
C
      CALL QEXIT('SO_DENSAI2')
C
      RETURN
      END
