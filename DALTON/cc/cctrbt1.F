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
#include "implicit.h"
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
C
      DIMENSION XINT(*),DSRHF(*),XLAMDP(*),WORK(LWORK)
      LOGICAL SQRINT
C
#include "ccorb.h"
#include "ccsdsym.h"
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
