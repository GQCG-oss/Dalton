C
C...   Copyright (c) 2001 by the authors of Dalton (see below).
C...   All Rights Reserved.
C...
C...   The source code in this file is part of
C...   "Dalton, a molecular electronic structure program, Release 1.2
C...   (2001), written by T. Helgaker, H. J. Aa. Jensen, P. Joergensen,
C...   J. Olsen, K. Ruud, H. Aagren, A.A. Auer, K.L. Bak, V. Bakken,
C...   O. Christiansen, S. Coriani, P. Dahle, E. K. Dalskov,
C...   T. Enevoldsen, B. Fernandez, C. Haettig, K. Hald, A. Halkier,
C...   H. Heiberg, H. Hettema, D. Jonsson, S. Kirpekar, R. Kobayashi,
C...   H. Koch, K. V. Mikkelsen, P. Norman, M. J. Packer,
C...   T. B. Pedersen, T. A. Ruden, A. Sanchez, T. Saue, S. P. A. Sauer,
C...   B. Schimmelpfennig, K. O. Sylvester-Hvid, P. R. Taylor,
C...   and O. Vahtras"
C...
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
*=====================================================================*
      SUBROUTINE CC_MAOMAO(TRANSA,TRANSB,ALPHA,AMAT,ISYMA,BMAT,ISYMB,
     &                     BETA,CMAT,ISYMC)
*---------------------------------------------------------------------*
*
*     Purpose: DGEMM like multiplication of two symmetry blocked
*              AO matrices
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
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
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
         WRITE (LUPRI,*) 'ISYMA,ISYMB,ISYMC:',ISYMA,ISYMB,ISYMC
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

         LDA = MAX(NBAS(ISYA1),1)
         LDB = MAX(NBAS(ISYB1),1)
         LDC = MAX(NBAS(ISYC1),1)

         KOFF1 = IAODIS(ISYA1,ISYA2) + 1
         KOFF2 = IAODIS(ISYB1,ISYB2) + 1
         KOFF3 = IAODIS(ISYC1,ISYC2) + 1

         CALL DGEMM(TRANSA,TRANSB,NBAS(ISYC1),NBAS(ISYC2),LENK,
     &              ALPHA,AMAT(KOFF1),LDA,BMAT(KOFF2),LDB,
     &              BETA, CMAT(KOFF3),LDC)

      END DO

      RETURN
      END
*=====================================================================*
