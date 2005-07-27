C
C...   Copyright (c) 2005 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 2.0
C...   (2005), written by C. Angeli, K. L. Bak,  V. Bakken, 
C...   O. Christiansen, R. Cimiraglia, S. Coriani, P. Dahle,
C...   E. K. Dalskov, T. Enevoldsen, B. Fernandez, C. Haettig,
C...   K. Hald, A. Halkier, H. Heiberg, T. Helgaker, H. Hettema, 
C...   H. J. Aa. Jensen, D. Jonsson, P. Joergensen, S. Kirpekar, 
C...   W. Klopper, R.Kobayashi, H. Koch, O. B. Lutnaes, K. V. Mikkelsen, 
C...   P. Norman, J.Olsen, M. J. Packer, T. B. Pedersen, Z. Rinkevicius,
C...   E. Rudberg, T. A. Ruden, K. Ruud, P. Salek, A. Sanchez de Meras,
C...   T. Saue, S. P. A. Sauer, B. Schimmelpfennig, K. O. Sylvester-Hvid, 
C...   P. R. Taylor, O. Vahtras, D. J. Wilson, H. Agren.
C...   This source code is provided under a written licence and may be
C...   used, copied, transmitted, or stored only in accord with that
C...   written licence.
C...
C...   In particular, no part of the source code or compiled modules may
C...   be distributed outside the research group of the licence holder.
C...   This means also that persons (e.g. post-docs) leaving the research
C...   group of the licence holder may not take any part of Dalton,
C...   including modified files, with him/her, unless that person has
C...   obtained his/her own licence.
C...
C...   For questions concerning this copyright write to:
C...      dalton-admin@kjemi.uio.no
C...
C...   For information on how to get a licence see:
C...      http://www.kjemi.uio.no/software/dalton/dalton.html
C
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
      DIMENSION PKMAT(*),SQMAT(NDIM,NDIM),SKMAT(NDIM,NDIM)
 
C
      IF (ISYMG .EQ. 0) THEN
       DO  I = 1,NDIM
          DO  J = 1,I
C
             IJ = I*(I-1)/2 + J
             SQMAT(I,J) =  PKMAT(IJ) 
             SQMAT(J,I) = - PKMAT(IJ) 
          ENDDO
       ENDDO

      ELSE
      DO I = 1,NDIM
         DO J = 1,I
             SKMAT(I,J) = SQMAT(I,J)
             SKMAT(J,I) = SQMAT(J,I)  
         ENDDO
      ENDDO  
      DO 100 I = 1,NDIM
          DO 110 J = 1,I

             IJ = I*(I-1)/2 + J
             SQMAT(I,J) =  PKMAT(IJ)  -  SKMAT(I,J)
             SQMAT(J,I) =  -PKMAT(IJ)  -  SKMAT(J,I)

  110     CONTINUE
  100 CONTINUE
      ENDIF

      RETURN
      END
