
      SUBROUTINE SO_PROPMO(ISYMPROP,PROPVAL,DOAB,
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
C        LABEL       Label of property integrals to be read         
C        ISYMPROP    Symmetry of the perturbed density
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
      LOGICAL,INTENT(IN) ::  DOAB
      REAL*8, INTENT(OUT)::  PROPVAL 
      REAL*8, INTENT(IN) ::  PDENSIJ(LPDENSIJ), PDENSAI(LPDENSAI),
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

      CHARACTER*8       RTNLBL(2)

      CALL QENTER('SO_PROPMO')
C
      PROPVAL = 0.0D0
C
C-----------------------------------------------
C     Do P_{ai} * ^pD_{ai} AND P_{ji} * ^pD_{ji} 
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
            PROPVAL = PROPVAL + DDOT(NRHF(ISYMJ),PRPMO(IOFFP),1,
     &                               PDENSIJ(IOFFDOO),1)
            IOFFP = IOFFP + NRHF(ISYMJ)          
            IOFFDOO = IOFFDOO + NRHF(ISYMJ)
C
C           HANDLE THE VO BLOCK
            
            PROPVAL = PROPVAL + DDOT(NVIR(ISYMA),PRPMO(IOFFP),1,
     &                               PDENSAI(IOFFDVO),1)
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
               PROPVAL = PROPVAL + DDOT(NVIR(ISYMA),PRPMO(IOFFP),1,
     &                               PDENSAB(IOFFD),1)
               IOFFP = IOFFP + NVIR(ISYMA)
               IOFFD = IOFFD + NVIR(ISYMA)
            END DO
         END DO
C         
      END IF         
C
      CALL QEXIT('SO_PROPMO')

      END


