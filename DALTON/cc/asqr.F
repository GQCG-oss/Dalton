!
!  Dalton, a molecular electronic structure program
!  Copyright (C) 2018 by the authors of Dalton.
!
!  This program is free software; you can redistribute it and/or
!  modify it under the terms of the GNU Lesser General Public
!  License version 2.1 as published by the Free Software Foundation.
!
!  This program is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
!  Lesser General Public License for more details.
!
!  If a copy of the GNU LGPL v2.1 was not distributed with this
!  code, you can obtain one at https://www.gnu.org/licenses/old-licenses/lgpl-2.1.en.html.
!
!
      SUBROUTINE CCSD_ASYMSQ(DISTAB,ISYMAB,SCR,ISYMG,ISYMD)
C
C     Antisymmetric Squareup of the integral distribution,
C     for orbit-orbit Breit-Pauli correction
C     S. Coriani, April 2003. Based on CCSD_SYMSQO.
C     Modified for [T1+T2,r12]-Integrals (Elena Vollmer, September 2003)
#include "implicit.h"
      DIMENSION DISTAB(*), SCR(*)
      PARAMETER (ONE = 1.0D0)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
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
