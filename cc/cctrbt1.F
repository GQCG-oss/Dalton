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
