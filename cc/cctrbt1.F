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
C  /* Deck cctrbt1 */
      SUBROUTINE CCTRBT1(XINT,DSRHF,XLAMDP,ISYMLP,WORK,LWORK,
     &                   ISYDIS,IOPT,SQRINT)
*---------------------------------------------------------------------*
*
*     Purpose: Transform three-index integral batch.
*
*     Written by Henrik Koch 3-Jan-1994
*     Symmetry by Henrik Koch and Alfredo Sanchez. 12-July-1994
*     Ove Christiansen 14-6-1996: General sym. lambda matrix ISYMLP
*
*     if IOPT = 0 overwrite result matrix
*     if IOPT = 1 add to previous
*     Sonia 19.8.99
*
*     SQRINT = .FALSE.  (alpha,beta) triangular packed integrals
*     SQRINT = .TRUE.   (alpha,beta) squared integrals
*     Christof Haettig 7.11.99
*
*=====================================================================*
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION XINT(*),DSRHF(*),XLAMDP(*),WORK(LWORK)
      LOGICAL SQRINT
C
#include <ccorb.h>
#include <ccsdsym.h>
C
      IF (IOPT.EQ.0) THEN
        OPTION = ZERO
      ELSE IF (IOPT.EQ.1) THEN
        OPTION = ONE
      ELSE
        CALL QUIT('Unknown option in CCTRBT1')
      ENDIF

      DO ISYMG = 1, NSYM
C
         ISYMJ  = MULD2H(ISYMLP,ISYMG)
         ISYMAB = MULD2H(ISYMG,ISYDIS)
C
         IF (SQRINT) THEN
           NDIMAB = N2BST(ISYMAB)
           KOFF1  = IDSAOGSQ(ISYMG,ISYDIS) + 1
           KOFF2  = IGLMRH(ISYMG,ISYMJ)    + 1
           KOFF3  = IDSRHFSQ(ISYMAB,ISYMJ) + 1
         ELSE
           NDIMAB = NNBST(ISYMAB)
           KOFF1  = IDSAOG(ISYMG,ISYDIS) + 1
           KOFF2  = IGLMRH(ISYMG,ISYMJ)  + 1
           KOFF3  = IDSRHF(ISYMAB,ISYMJ) + 1
         END IF
C
         NBASG  = MAX(NBAS(ISYMG),1)
C
         IF (NDIMAB.GT.0) THEN
           CALL DGEMM('N','N',NDIMAB,NRHF(ISYMJ),NBAS(ISYMG),
     *                ONE,XINT(KOFF1),NDIMAB,XLAMDP(KOFF2),NBASG,
     *                OPTION,DSRHF(KOFF3),NDIMAB)
         END IF
C
      END DO
C
      RETURN
      END
*=====================================================================*
