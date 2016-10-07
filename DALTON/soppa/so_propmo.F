
      SUBROUTINE SO_PROPMO(ISYMPROP,PROPVAL,DOAB,IMAGPROP,
     &                     PRPMO,LPRPMO,
     &                     PDENSIJ,LPDENSIJ,
     &                     PDENSAB,LPDENSAB,
     &                     PDENSAI,LPDENSAI)
C
C     Rasmus Faber, Nov. 2015
C
C     PURPOSE: Contract the SOPPA perturbed density matrices with
C              the relevant property integrals
C              First implementation use simply densities and properties
C              in MO basis.
C
C     INPUT:
C        ISYMPROP    Symmetry of the perturbed density
C        DOAB        Calculate also contribution from VV part
C        IMAGPROP    Are the property operator imaginary
C        PRPMO(LPRPMO) Property integrals (MO basis)
C                 Perturbed density matrices
C        PDENSIJ(LPDENSIJ), PDENSAB(LPDENSAB), PDENSAI(LPDENSAI)
C
C     OUTPUT:
C        PROPVAL     Result: Value of the property
C
C     :
c#include "implicit.h"
      implicit none

C Symmetry offsets
#include "ccorb.h"
#include "ccsdsym.h"
#include "soppinf.h"
C Arguments
      CHARACTER*8       LABEL
      INTEGER,INTENT(IN) ::  ISYMPROP,
     &             LPRPMO, LPDENSIJ, LPDENSAB, LPDENSAI
      LOGICAL,INTENT(IN) ::  DOAB, IMAGPROP
      REAL(8), INTENT(OUT)::  PROPVAL
      REAL(8), INTENT(IN) ::  PDENSIJ(LPDENSIJ), PDENSAI(LPDENSAI),
     6                       PDENSAB(LPDENSAB),
     &                       PRPMO(LPRPMO)
C External routines
      EXTERNAL DDOT
      DOUBLE PRECISION DDOT
C Local variables
C      INTEGER           KPRP1
c      INTEGER           LPRP1
C      INTEGER           KEND1
C      INTEGER           LWORK1
      INTEGER           ISYMA, ISYMB, ISYMI, ISYMJ
      INTEGER           IOFFP, IOFFD, IOFFDOO, IOFFDVO
      DOUBLE PRECISION  PROPVAL1, PROPVAL2

      CHARACTER*8       RTNLBL(2)

      CALL QENTER('SO_PROPMO')
C
C
      PROPVAL1 = 0.0D0
      PROPVAL2 = 0.0D0
C
C-----------------------------------------------
C     Do P_{ai} * ^pD_{ai}* AND P_{ji} * ^pD_{ji}
C-----------------------------------------------
C
      IOFFP = 1
      DO ISYMI = 1, NSYM
         ISYMA = MULD2H(ISYMI,ISYMPROP)
         ISYMJ = ISYMA
C
C         IOFFDOO = IIJDEN(ISYMJ,ISYMI) + 1
C  Should for some reason be like this
         IOFFDOO = IIJDEN(ISYMI,ISYMJ) + 1
         IOFFDVO = IAIDEN(ISYMA,ISYMI) + 1
C         print *, ioffdvo, nrhf(isymi),it1am(isyma,isymi)
         DO I = 1, NRHF(ISYMI)
C           HANDLE THE OO BLOCK
            PROPVAL1 = PROPVAL1 + DDOT(NRHF(ISYMJ),PRPMO(IOFFP),1,
     &                               PDENSIJ(IOFFDOO),1)
            IOFFP = IOFFP + NRHF(ISYMJ)
            IOFFDOO = IOFFDOO + NRHF(ISYMJ)
C
C           HANDLE THE VO BLOCK

            PROPVAL2 = PROPVAL2 +
     &                     DDOT(NVIR(ISYMA),PRPMO(IOFFP),1,
     &                          PDENSAI(IOFFDVO),1)
            IOFFP = IOFFP + NVIR(ISYMA)
            IOFFDVO = IOFFDVO + NVIR(ISYMA)
         END DO
      END DO
C
C --> this would be the end for RPA
      IF ( DOAB ) THEN
C
C---------------------------------
C     Do P_{ab} * ^pD_{ab}
C---------------------------------
C
! IOFFP should already have this value
         IOFFP = NLRHFR(ISYMPROP) + 1
         DO ISYMB = 1, NSYM
            ISYMA = MULD2H(ISYMB,ISYMPROP)
C
C         IOFFD = IABDEN(ISYMA,ISYMB) + 1
C         print *, ioffd, iabden(isymb,isyma)+1
C     Skip elements in the OV block
            IOFFD = IABDEN(ISYMB,ISYMA) + 1
            DO I = 1, NVIR(ISYMB)
               IOFFP = IOFFP + NRHF(ISYMA)
               PROPVAL1 = PROPVAL1 + DDOT(NVIR(ISYMA),PRPMO(IOFFP),1,
     &                               PDENSAB(IOFFD),1)
               IOFFP = IOFFP + NVIR(ISYMA)
               IOFFD = IOFFD + NVIR(ISYMA)
            END DO
         END DO
C
      END IF
C
C
C  For imaginary perturbations, there is an
C  overall factor of i*i = -1,
C  but VO cont. is complex conjugated, canceling this factor
      IF (IMAGPROP) THEN
         PROPVAL = -PROPVAL1 + PROPVAL2
      ELSE
         PROPVAL =  PROPVAL1 + PROPVAL2
      ENDIF
C
C
      CALL QEXIT('SO_PROPMO')

      END


