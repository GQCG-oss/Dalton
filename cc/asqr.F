!
!...   Copyright (c) 2011 by the authors of Dalton (see below).
!...   All Rights Reserved.
!...
!...   The source code in this file is part of
!...   "Dalton, a molecular electronic structure program,
!...    Release DALTON2011 (2011), see http://daltonprogram.org"
!...
!...   This source code is provided under a written licence and may be
!...   used, copied, transmitted, or stored only in accord with that
!...   written licence.
!...
!...   In particular, no part of the source code or compiled modules may
!...   be distributed outside the research group of the licence holder.
!...   This means also that persons (e.g. post-docs) leaving the research
!...   group of the licence holder may not take any part of Dalton,
!...   including modified files, with him/her, unless that person has
!...   obtained his/her own licence.
!...
!...   For further information, including how to get a licence, see:
!...      http://daltonprogram.org
!
!
      SUBROUTINE CCSD_ASYMSQ(DISTAB,ISYMAB,SCR,ISYMG,ISYMD)
C
C     Antisymmetric Squareup of the integral distribution,
C     for orbit-orbit Breit-Pauli correction
C     S. Coriani, April 2003. Based on CCSD_SYMSQO.
C     Modified for [T1+T2,r12]-Integrals (Elena Vollmer, September 2003)
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
            CALL ASQMATR(NBAS(ISYMB),DISTAB(KOFF1),SCR(KOFF2),ISYMG)
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
                 IF (ISYMG .EQ. 0) THEN
                  CALL DCOPY(NBAS(ISYMA),DISTAB(KOFF2),1,SCR(KOFF3),
     *                       NBAS(ISYMB))
                 ELSE
                 CALL DSCAL(NBAS(ISYMA),-ONE,SCR(KOFF3),NBAS(ISYMB))
                 CALL DAXPY(NBAS(ISYMA),ONE,DISTAB(KOFF2),1,SCR(KOFF3),
     *                      NBAS(ISYMB))
                 END IF
                  KOFF2 = KOFF2 + NBAS(ISYMA)
                  KOFF3 = KOFF3 + 1
  210          CONTINUE
C
               KOFF4 = IAODIS(ISYMA,ISYMB) + 1
               IF (ISYMG .EQ. 0) THEN
                CALL DCOPY(NTOT,DISTAB(KOFF1),1,SCR(KOFF4),1)
                CALL DSCAL(NTOT,-ONE,SCR(KOFF4),1)
               ELSE
                CALL DSCAL(NTOT,-ONE,SCR(KOFF4),1)
                CALL DAXPY(NTOT,-ONE,DISTAB(KOFF1),1,SCR(KOFF4),1)
               END IF
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
      SUBROUTINE ASQMATR(NDIM,PKMAT,SQMAT,ISYMG)
C
C     PURPOSE:
C      Antisymmetric Square up packed matrix.
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION PKMAT(*),SQMAT(NDIM,NDIM)
 
C
      IF (ISYMG .EQ. 0) THEN

       IJ = 0
       DO  I = 1,NDIM
          DO  J = 1,I
C
             IJ = IJ + 1
             SQMAT(I,J) =   PKMAT(IJ) 
             SQMAT(J,I) = - PKMAT(IJ) 
          ENDDO
       ENDDO

      ELSE
      IJ = 0
      DO 100 I = 1,NDIM
          DO 110 J = 1,I
           IJ = IJ + 1
           IF (I.EQ.J) THEN
             SQMAT(J,I) =  - PKMAT(IJ)  -  SQMAT(J,I)
           ELSE
             SQMAT(I,J) =    PKMAT(IJ)  -  SQMAT(I,J)
             SQMAT(J,I) =  - PKMAT(IJ)  -  SQMAT(J,I)
           END IF
  110     CONTINUE
  100 CONTINUE
      ENDIF

      RETURN
      END
