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
C /* Deck cctrbt2 */
      SUBROUTINE CCTRBT2(XINT,DSRHF,XLAMDP,ISYMLP,WORK,LWORK,
     &                        ISYDIS,IOPT,LSQRINT,LSQRUP,SGNINT)
*---------------------------------------------------------------------*
*
*     Purpose: Transform gamma index of integral batch 
*              I_{al be, gamma}^del to occupied.
*
*     XLAMDP,ISYMLP = lambda matrix and its symmetry
*     XINT, ISYDIS  = I_{al be, gamma} batch and its symmetry
*     Options:
*       if IOPT = 0 overwrite result matrix
*       if IOPT = 1 add to previous
*       LSQRINT = TRUE, (alpha beta|* *) is full matrix (not packed) 
*       LSQRUP  = TRUE, square up (a b| after transformation of gamma 
*                       to k 
*       SGNINT  = sign of integral distribution 
*
*   Written by Sonia Coriani 19-11-99, based on CCTRBT
*
*=====================================================================*
#include <implicit.h>
      PARAMETER(ZERO = 0.0D0, ONE = 1.0D0)
*
      DIMENSION XINT(*),DSRHF(*),XLAMDP(*),WORK(LWORK)
      LOGICAL LSQRUP,LSQRINT
*
#include <ccorb.h>
#include <ccsdsym.h>
*
      IF (IOPT.EQ.0) THEN
         FAC = ZERO
      ELSE IF (IOPT.EQ.1) THEN
         FAC = ONE
      ELSE
        CALL QUIT('Unknown option in CCTRBT2')
      ENDIF
*
* memory check when squaring
*
      IF (LSQRUP) THEN
         DO ISYMJ = 1, NSYM
            ISYMG    = MULD2H(ISYMLP,ISYMJ)
            ISYMAB   = MULD2H(ISYMG,ISYDIS)
            ISYDSRHF = MULD2H(ISYMAB,ISYMJ)
            IF (LWORK.LT.NDSRHF(ISYDSRHF)) THEN
              CALL QUIT('Insufficient memory in CCTRBT2')
            END IF
         END DO
      END IF
*
* Calculate (al be|j)^del = sum_gam I^del_{al be, gam} Lambda_{gam j}
*
      DO ISYMJ = 1,NSYM
*
         ISYMG  = MULD2H(ISYMLP,ISYMJ)
         ISYMAB = MULD2H(ISYMG,ISYDIS)
         NBASG  = MAX(NBAS(ISYMG),1)

         KOFF2  = 1 + IGLMRH(ISYMG,ISYMJ)

         IF (LSQRINT) THEN
            KOFF1  = 1 + IDSAOGSQ(ISYMG,ISYDIS)
            KOFF3  = 1 + IDSRHFSQ(ISYMAB,ISYMJ)  
            NDIMAB = N2BST(ISYMAB) 
         ELSE
            KOFF1  = 1 + IDSAOG(ISYMG,ISYDIS)
            KOFF3  = 1 + IDSRHF(ISYMAB,ISYMJ)
            NDIMAB = NNBST(ISYMAB)
         END IF

         NALBEM = MAX(NDIMAB,1)

         IF (LSQRUP) THEN

            CALL DGEMM('N','N',NDIMAB,NRHF(ISYMJ),NBAS(ISYMG),
     *                 ONE,XINT(KOFF1),NALBEM,XLAMDP(KOFF2),NBASG,
     *                 ZERO,WORK,NALBEM)            

            ! Resort (al>=be, k) to (al be| k) 
            ! Put in DSRHF which is dimensioned full (a b| from input
            DO J = 1, NRHF(ISYMJ)
              KOFF4 = NNBST(ISYMAB)*(J-1) + 1
              KOFF5 = IDSRHFSQ(ISYMAB,ISYMJ) + N2BST(ISYMAB)*(J-1) + 1
              CALL CCSD_SYMSQ(WORK(KOFF4),ISYMAB,DSRHF(KOFF5))
            END DO

         ELSE
            CALL DGEMM('N','N',NDIMAB,NRHF(ISYMJ),NBAS(ISYMG),
     *              SGNINT,XINT(KOFF1),NALBEM,XLAMDP(KOFF2),NBASG,
     *                 FAC,DSRHF(KOFF3),NALBEM)
         END IF
 
      END DO
 
      RETURN
      END
*---------------------------------------------------------------------*
