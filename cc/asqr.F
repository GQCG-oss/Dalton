C  /* Deck ccsd_asymsq */
      SUBROUTINE CCSD_ASYMSQ(DISTAB,ISYMAB,SCR,ISYMG,ISYMD)
C
C     Antisymmetric Squareup of the integral distribution,
C     for orbit-orbit Breit-Pauli correction
C     S. Coriani, April 2003. Based on CCSD_SYMSQO.
C
#include <implicit.h>
      DIMENSION DISTAB(*), SCR(*)
      PARAMETER (ONE = 1.0D0)
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>
C
      CALL QENTER('CCSD_ASYMSQ')
C
      IF (ISYMAB .EQ. 1) THEN
C
         KOFF1 = 1
         KOFF2 = 1
         DO 100 ISYMB = 1,NSYM
            CALL ASQMATR(NBAS(ISYMB),DISTAB(KOFF1),SCR(KOFF2))
            KOFF1 = KOFF1 + NBAS(ISYMB)*(NBAS(ISYMB)+1)/2
            KOFF2 = KOFF2 + NBAS(ISYMB)*NBAS(ISYMB)
  100    CONTINUE
C
      ELSE
         KOFF1 = 1
         KOFF2 = 1
         DO 200 ISYMB = 1,NSYM
C
            ISYMA = MULD2H(ISYMB,ISYMAB)

            IF (ISYMB .GT. ISYMA) THEN
C
               NTOT  = NBAS(ISYMA)*NBAS(ISYMB)
C
               KOFF2 = KOFF1
               KOFF3 = IAODIS(ISYMB,ISYMA) + 1
               DO 210 B = 1,NBAS(ISYMB)
                  CALL DCOPY(NBAS(ISYMA),DISTAB(KOFF2),1,SCR(KOFF3),
     *                       NBAS(ISYMB))
                  KOFF2 = KOFF2 + NBAS(ISYMA)
                  KOFF3 = KOFF3 + 1
  210          CONTINUE
C
               KOFF4 = IAODIS(ISYMA,ISYMB) + 1
               CALL DCOPY(NTOT,DISTAB(KOFF1),1,SCR(KOFF4),1)
               CALL DSCAL(NTOT,-ONE,SCR(KOFF4),1)
C
               KOFF1 = KOFF1 + NTOT
C
            ENDIF
C
  200    CONTINUE
C
      ENDIF
C
      CALL QEXIT('CCSD_ASYMSQ')
C
      RETURN
      END
C  /* Deck asqmatr */
      SUBROUTINE ASQMATR(NDIM,PKMAT,SQMAT)
C
C     PURPOSE:
C      Antisymmetric Square up packed matrix.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PKMAT(*),SQMAT(NDIM,NDIM)
C
      DO 100 I = 1,NDIM
         DO 110 J = 1,I
C
            IJ = I*(I-1)/2 + J
            SQMAT(I,J) = PKMAT(IJ)
            SQMAT(J,I) = -PKMAT(IJ)
C
  110    CONTINUE
  100 CONTINUE
C
      RETURN
      END
