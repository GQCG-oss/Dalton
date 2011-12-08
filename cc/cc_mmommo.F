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
#include <single.h>
C
*=====================================================================*
      SUBROUTINE CC_MMOMMO(TRANSA,TRANSB,ALPHA,AMAT,ISYMA,BMAT,ISYMB,
     &                     BETA,CMAT,ISYMC)
*---------------------------------------------------------------------*
*
*     Purpose: DGEMM like multiplication of two symmetry blocked
*              MO matrices
*  
*               CMAT := alpha AMAT x BMAT + beta CMAT
*              
*               TRANSA: 'N'/'T' transpose / do not transpose A
*               TRANSB: 'N'/'T' transpose / do not transpose B
*
*     Christof Haettig, March 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccsdsym.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYMC, ISYMA, ISYMB
      CHARACTER*(*) TRANSA, TRANSB

#if defined (SYS_CRAY)
      REAL AMAT(*), BMAT(*), CMAT(*), ALPHA, BETA
#else
      DOUBLE PRECISION AMAT(*), BMAT(*), CMAT(*), ALPHA, BETA
#endif

      INTEGER ISYA1, ISYA2, ISYB1, ISYB2, ISYC1, ISYC2, LENK
      INTEGER LDA, LDB, LDC, KOFF1, KOFF2, KOFF3

*---------------------------------------------------------------------*
*     check if the symmetries match:
*---------------------------------------------------------------------*
      IF ( ISYMC .NE. MULD2H(ISYMA,ISYMB) ) THEN
         WRITE (LUPRI,*) 'Symmetry mismatch in CC_MAOMAO.'
         CALL QUIT('Symmetry mismatch in CC_MAOMAO.')
      END IF

*---------------------------------------------------------------------*
*     do the matrix multiplication using DGEMM
*---------------------------------------------------------------------*
      DO ISYA1 = 1, NSYM

         ISYA2 = MULD2H(ISYA1,ISYMA)

         IF      (TRANSA(1:1).EQ.'N' .OR. TRANSA(1:1).EQ.'n') THEN
            ISYC1  = ISYA1
            LENK   = NBAS(ISYA2)
         ELSE IF (TRANSA(1:1).EQ.'T' .OR. TRANSA(1:1).EQ.'t') THEN
            ISYC1  = ISYA2
            LENK   = NBAS(ISYA1)
         END IF

         ISYC2  = MULD2H(ISYC1,ISYMC)

         IF      (TRANSB(1:1).EQ.'N' .OR. TRANSB(1:1).EQ.'n') THEN
            ISYB2  = ISYC2
            ISYB1  = MULD2H(ISYB2,ISYMB)
         ELSE IF (TRANSB(1:1).EQ.'T' .OR. TRANSB(1:1).EQ.'t') THEN
            ISYB1  = ISYC1
            ISYB2  = MULD2H(ISYB1,ISYMB)
         END IF

         LDA = MAX(NORBS(ISYA1),1)
         LDB = MAX(NORBS(ISYB1),1)
         LDC = MAX(NORBS(ISYC1),1)

         KOFF1 = IAODIS(ISYA1,ISYA2) + 1
         KOFF2 = IAODIS(ISYB1,ISYB2) + 1
         KOFF3 = IAODIS(ISYC1,ISYC2) + 1

         CALL DGEMM(TRANSA,TRANSB,NORBS(ISYC1),NORBS(ISYC2),LENK,
     &              ALPHA,AMAT(KOFF1),LDA,BMAT(KOFF2),LDB,
     &              BETA, CMAT(KOFF3),LDC)

      END DO

      RETURN
      END
*=====================================================================*
