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
      SUBROUTINE CC_GIM1(DSRHF,ISYRHF,XMGD,ISYMGD,GIM,IOPT,LSQRAB,
     &                                                 WORK,LWORK)
*---------------------------------------------------------------------*
*
* Purpose: Calculate G intermediate from effective BF density and
*          (**|k del) integrals
*
*    IOPT = 0  :  compute left-hand-side GZeta intermediate
*
*    IOPT = 1  :  compute right-hand-side G intermediate, using
*                 2 Cou - Exc combination of XMGD 
*                 Note: this overwrites XMGD!
*
*    symmetries:     ISYRHF  --  DSRHF
*                    ISYMGD  --  XMGD
*
*    LSQRAB =.TRUE. Distribution DSRHF is already squared in (al bet|
*                   otherwise it is packed and it gets squared up here
*
* Christof Haettig, November 1998
* LSQRAB option, Sonia Coriani, November 1999
*
*---------------------------------------------------------------------*
#include <implicit.h>
#include <ccorb.h>
#include <ccsdsym.h>
      PARAMETER(ZERO=0.0D0,ONE=1.0D0)
      DIMENSION DSRHF(*),XMGD(*),GIM(*),WORK(LWORK)
      LOGICAL LSQRAB
C
C     ----------------------------------------------------------
C     if right-hand-side G requested replace XMGD by 2Cou - Exc:
C     ----------------------------------------------------------
      IF (IOPT.EQ.1) THEN
         CALL CC_MGDTCME(XMGD,ISYMGD,WORK,LWORK)
      END IF
C
      DO ISYMK = 1,NSYM
C
         ISYMAG = MULD2H(ISYMK,ISYRHF)
         ISYMGI = MULD2H(ISYMK,ISYMGD)
C
         IF ((.NOT.LSQRAB).AND.(LWORK .LT. N2BST(ISYMAG))) THEN
            CALL QUIT('Insufficient space in CC_GIM')
         ENDIF
C
         DO K = 1,NRHF(ISYMK)
C
            IF (LSQRAB) THEN
              KOFF1 = IDSRHFSQ(ISYMAG,ISYMK) + N2BST(ISYMAG)*(K-1) + 1
            ELSE
              KOFF1 = IDSRHF(ISYMAG,ISYMK) + NNBST(ISYMAG)*(K - 1) + 1
              CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYMAG,WORK)
            END IF
C
            DO ISYMI = 1,NSYM
C
               ISYMG = MULD2H(ISYMI,ISYMGI)
               ISYMA = MULD2H(ISYMG,ISYMAG)
C
               KOFF5 = IT2BGD(ISYMGI,ISYMK) + NT1AO(ISYMGI)*(K - 1) 
     *               + IT1AO(ISYMG,ISYMI) + 1
               KOFF6 = IT1AO(ISYMA,ISYMI) + 1
C
               NBASG = MAX(NBAS(ISYMG),1)
               NBASA = MAX(NBAS(ISYMA),1)
C
               IF (LSQRAB) THEN
                 KOFF4 = IAODIS(ISYMA,ISYMG) + KOFF1
                 CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     *                       ONE,DSRHF(KOFF4),NBASA,XMGD(KOFF5),NBASG,
     *                       ONE,GIM(KOFF6), NBASA)
               ELSE
                 KOFF4 = IAODIS(ISYMA,ISYMG) + 1
                 CALL DGEMM('N','N',NBAS(ISYMA),NRHF(ISYMI),NBAS(ISYMG),
     *                       ONE,WORK(KOFF4),NBASA,XMGD(KOFF5),NBASG,
     *                       ONE,GIM(KOFF6), NBASA)
               END IF
C
            END DO
         END DO
      END DO
 
      RETURN
      END
*=====================================================================*
*                   END OF SUBROUTINE CC_GIM1                         *
*=====================================================================*
