*=====================================================================*
      SUBROUTINE CCKAPPASQ(KAPPASQ,KAPPA,ISYKAP,TRANS)
*---------------------------------------------------------------------*
*
*     Purpose: resort kappa vector to full matrix scheme
*
*              TRANS = 'N' :  KAPPASQ <-- KAPPA
*              TRANS = 'T' :  KAPPASQ <-- KAPPA^T
*
*     Christof Haettig 8-2-1999
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <priunit_dec.h>
#include <priunit.h>
#include <ccorb_dec.h>
#include <ccorb.h>
#include <ccsdsym_dec.h>
#include <ccsdsym.h>
#include <ccfro.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      CHARACTER*(1) TRANS
      INTEGER ISYKAP

#if defined (SYS_CRAY)
      REAL KAPPASQ(*), KAPPA(*)
#else
      DOUBLE PRECISION KAPPASQ(*), KAPPA(*)
#endif

      INTEGER ISYMA, ISYMI, IORBA, IORBI, KKIA, KKAI, KSAI, KSIA

*---------------------------------------------------------------------*
*     resort kappa vector:
*---------------------------------------------------------------------*
      CALL DZERO(KAPPASQ,N2BST(ISYKAP))

      DO ISYMI = 1, NSYM
         ISYMA = MULD2H(ISYMI,ISYKAP)

         DO I = 1, NRHFS(ISYMI)
         DO A = 1, NVIRS(ISYMA)

            IORBI = I
            IORBA = NRHFS(ISYMA) + A

            KKAI = IALLAI(ISYMA,ISYMI) + (I-1)*NVIRS(ISYMA) + A
            KSAI = IAODIS(ISYMA,ISYMI) + (IORBI-1)*NORBS(ISYMA) + IORBA

            KKIA = NALLAI(ISYKAP) + KKAI
            KSIA = IAODIS(ISYMI,ISYMA) + (IORBA-1)*NORBS(ISYMI) + IORBI

            IF      (TRANS.EQ.'N' .OR. TRANS.EQ.'n') THEN
C              KAPPASQ(KSAI) = -KAPPA(KKIA)
               KAPPASQ(KSAI) = KAPPA(KKAI)
               KAPPASQ(KSIA) = KAPPA(KKIA)
C              KAPPASQ(KSIA) = - KAPPA(KKAI)
            ELSE IF (TRANS.EQ.'T' .OR. TRANS.EQ.'t') THEN
C              KAPPASQ(KSIA) = -KAPPA(KKIA)
               KAPPASQ(KSIA) = KAPPA(KKAI)
               KAPPASQ(KSAI) = KAPPA(KKIA)
C              KAPPASQ(KSAI) = - KAPPA(KKAI)
            ELSE
               CALL QUIT('Illegal value of TRANS in CCKAPPASQ.')
            END IF

         END DO
         END DO

      END DO

*---------------------------------------------------------------------*
*     print to output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CCKAPPASQ> input kappa vector:'
         WRITE (LUPRI,'(5X,I5,F12.8)') (I,KAPPA(I),I=1,2*NALLAI(ISYKAP))
         WRITE (LUPRI,*) 'CCKAPPASQ> resorted orbital '//
     &        'relaxation matrix:'
         CALL CC_PRONELAO(KAPPASQ,ISYKAP)
      END IF

      RETURN
      END
*======================================================================*
