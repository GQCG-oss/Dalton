C  /* Deck so_t2ao */
      SUBROUTINE SO_T2AO(T2AM,XLAMDH,ISYMLH,SCRM,WORK,LWORK,
     *                   IDEL,ISYMD,ISYMTR)
C
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C     The Orgininal routine: CC_T2AO, is written by Henrik Koch 22-dec-1993.
C 
C
C    Symmetry due to Alfredo Sanchez and Henrik Koch 11-July 1994
C     Nontotal symmetric amplitudes Ove Christiansen 14-2-1995.
C     LAMDH is still assumed tot. sym.
C     Asger Halkier 13/2-1996: Generalised to handle "non-direct"
C     AO-index gamma in lampda matrix (IOPT = 2), as well as the
C     usual "direct" delta AO-index (IOPT = 1).
C     Ove Christiansen 16-6-1996:
C     Generalised to non-total symmetric Lamdba matrices.
C     PURPOSE:
C             Tdjci -> Tci,j (delta)
C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
C
#include "implicit.h"
#include "priunit.h"
#include "iratdef.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
      DIMENSION T2AM(*),XLAMDH(*)
      DIMENSION SCRM(*),WORK(LWORK)
#include "ccorb.h"
#include "ccsdsym.h"
C
C-----------------------------------------------------
C     Calculate the transformed t2-amplitude and save.
C-----------------------------------------------------
C
      ISYDVI = MULD2H(ISYMLH,ISYMD)
      ISYMM = MULD2H(ISYMTR,ISYDVI)
C
      CALL DZERO(SCRM,NT2BCD(ISYMM))
C
      IF ( LWORK .LT. NVIR(ISYDVI)) THEN
         CALL QUIT('Insufficient core in SO_T2AO')
      ENDIF
C
      CALL DZERO(WORK,NVIR(ISYDVI))
C
      KOFF1 = IGLMVI(ISYMD,ISYDVI) + IDEL - IBAS(ISYMD)
C
      CALL DCOPY(NVIR(ISYDVI),XLAMDH(KOFF1),NBAS(ISYMD),WORK,1)
C       Copies XLAMDDH(KOFF1) to WORK

      DO 100 ISYMJ = 1,NSYM
C
         ISYMDJ = MULD2H(ISYMJ,ISYDVI)
         ISYMCI = MULD2H(ISYMTR,ISYMDJ)
C
         NTOTCI = MAX(NT1AM(ISYMCI),1)
C
         DO 110 J = 1,NRHF(ISYMJ)
C
            KDJ   = IT1AM(ISYDVI,ISYMJ) + NVIR(ISYDVI)*(J-1) + 1
            KOFF2 = IT2SQ(ISYMCI,ISYMDJ)
     *            + NT1AM(ISYMCI)*(KDJ - 1) + 1
            KOFF3 = IT2BCD(ISYMCI,ISYMJ)
     *            + NT1AM(ISYMCI)*(J-1) + 1
C
            CALL DGEMV('N',NT1AM(ISYMCI),NVIR(ISYDVI),ONE,
     *                 T2AM(KOFF2),NTOTCI,WORK,1,ZERO,
     *                 SCRM(KOFF3),1)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
