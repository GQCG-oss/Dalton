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
#include <single.h>
C
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
#include <priunit.h>
#include <ccorb.h>
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
