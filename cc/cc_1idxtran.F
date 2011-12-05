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
      SUBROUTINE CC_1IDXTRAN(AMAT,ISYMA,BMAT,ISYMB,CMAT,ISYMC)
*---------------------------------------------------------------------*
*
*     Purpose: evaluate the 1-index transformation
*
*           C = A B + B A^T
*
*     Christof Haettig 7-2-1999
*
*---------------------------------------------------------------------*
      IMPLICIT NONE
#include <ccorb.h>
#include <ccsdsym.h>
#include <priunit.h>

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)
 
      INTEGER ISYMA, ISYMB, ISYMC

#if defined (SYS_CRAY)
      REAL AMAT(*), BMAT(*), CMAT(*), ONE
#else
      DOUBLE PRECISION AMAT(*), BMAT(*), CMAT(*), ONE
#endif
      PARAMETER( ONE = 1.0D0 )

      INTEGER ISYMP,ISYMQ,ISYMR,NBASP,NBASQ,NBASR,KOFF1,KOFF2,KOFF3

*---------------------------------------------------------------------*
*     check symmetries and initialize output matrix:
*---------------------------------------------------------------------*
      IF (ISYMC .NE. MULD2H(ISYMA,ISYMB)) THEN
         CALL QUIT('Symmetry mismatch in CC_1IDXTRAN.')
      END IF

      CALL DZERO(CMAT,N2BST(ISYMC))

*---------------------------------------------------------------------*
*     Calculate A x B  and add to output matrix:
*---------------------------------------------------------------------*
      DO ISYMP = 1, NSYM

         ISYMQ = MULD2H(ISYMP,ISYMA)
         ISYMR = MULD2H(ISYMQ,ISYMB)

         KOFF1 = IAODIS(ISYMP,ISYMQ) + 1
         KOFF2 = IAODIS(ISYMQ,ISYMR) + 1
         KOFF3 = IAODIS(ISYMP,ISYMR) + 1

         NBASP = MAX(1,NBAS(ISYMP))
         NBASQ = MAX(1,NBAS(ISYMQ))

         CALL DGEMM('N','N',NBAS(ISYMP),NBAS(ISYMR),NBAS(ISYMQ),
     *              ONE,AMAT(KOFF1),NBASP,BMAT(KOFF2),NBASQ,
     *              ONE,CMAT(KOFF3),NBASP)
      END DO

 
*---------------------------------------------------------------------*
*     Calculate B x A^T  and add to output matrix:
*---------------------------------------------------------------------*
      DO ISYMP = 1, NSYM

         ISYMQ = MULD2H(ISYMP,ISYMB)
         ISYMR = MULD2H(ISYMQ,ISYMA)

         KOFF1 = IAODIS(ISYMP,ISYMQ) + 1
         KOFF2 = IAODIS(ISYMR,ISYMQ) + 1
         KOFF3 = IAODIS(ISYMP,ISYMR) + 1

         NBASP = MAX(1,NBAS(ISYMP))
         NBASR = MAX(1,NBAS(ISYMR))

         CALL DGEMM('N','T',NBAS(ISYMP),NBAS(ISYMR),NBAS(ISYMQ),
     *              ONE,BMAT(KOFF1),NBASP,AMAT(KOFF2),NBASR,
     *              ONE,CMAT(KOFF3),NBASP)
      END DO

 
*---------------------------------------------------------------------*
*     print to output & return:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
         WRITE (LUPRI,*) 'CC_1IDXTRAN> result of one-index '//
     &        'transformation:'
         CALL CC_PRONELAO(CMAT,ISYMC)
      END IF

      RETURN
      END
*======================================================================*
