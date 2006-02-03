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
C  /* Deck cc_sort4o */
      SUBROUTINE CC_SORT4O(XKJIL,ISYM4O,XJKLI,IOPT)
*---------------------------------------------------------------------*
*     Purpose: resort I(kj,i;l) to I(jk,l;i)
*
*     IOPT = 1 : XJKLI area is initialized here
*     IOPT = 2 : XJKLI is added to gamma intermediate
*                       already stored as I(jk,l;i)
*
*     Sonia Coriani, 14/09-1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
      INTEGER ISYM4O, IOPT
#if defined (SYS_CRAY)
      REAL XKJIL(*), XJKLI(*)
      REAL ZERO, ONE, HALF, DDOT, XNORM, FAC
#else
      DOUBLE PRECISION XKJIL(*), XJKLI(*)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM, FAC
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*
      INTEGER ISYML, ISYKJI, ISYMI, ISYMKJ, ISYJKL, ISYMJK
      INTEGER ISYMK, ISYMJ,  KJIL,  JKLI
*
*     --------------------------------------
*     Initialize result area with zero's or add to previous
*     --------------------------------------
      IF (IOPT.EQ.1) THEN
        FAC = ZERO
      ELSE
        FAC = ONE 
      END IF
*     --------------------------------------
*     Reorder thru loops on all 4 indices
*     --------------------------------------
      DO ISYML = 1,NSYM
         ISYKJI = MULD2H(ISYM4O,ISYML)
         DO L = 1, NRHF(ISYML)
            DO ISYMI = 1, NSYM
               ISYMKJ = MULD2H(ISYKJI,ISYMI)
               ISYJKL = MULD2H(ISYM4O,ISYMI)
               ISYMJK = MULD2H(ISYJKL,ISYML)
               DO I = 1, NRHF(ISYMI)
                  DO ISYMK = 1, NSYM
                     ISYMJ = MULD2H(ISYMJK,ISYMK)
                     DO K = 1, NRHF(ISYMK)
                     DO J = 1, NRHF(ISYMJ)

                KJIL = I3ORHF(ISYKJI,ISYML) + NMAIJK(ISYKJI)*(L-1)+
     &                 IMAIJK(ISYMKJ,ISYMI) + NMATIJ(ISYMKJ)*(I-1)+
     &                 IMATIJ(ISYMK,ISYMJ)  + NRHF(ISYMK)*(J-1) + K

                JKLI = I3ORHF(ISYJKL,ISYMI) + NMAIJK(ISYJKL)*(I-1)+
     &                 IMAIJK(ISYMJK,ISYML) + NMATIJ(ISYMJK)*(L-1)+
     &                 IMATIJ(ISYMJ,ISYMK)  + NRHF(ISYMJ)*(K-1) + J

               XJKLI(JKLI) = FAC*XJKLI(JKLI) + XKJIL(KJIL)

                     END DO         !J
                     END DO         !K
                  END DO            !ISYMK
               END DO               !I
            END DO                  !ISYMI
         END DO                     !L
      END DO                        !ISYML
*     ---------------------------------------
*     Finished, return
*     ---------------------------------------
      RETURN
      END
*=====================================================================*
C
*=====================================================================*
C  /* Deck cc_sort4o2 */
      SUBROUTINE CC_SORT4O2(XKJIL,ISYM4O,XJLKI,IOPT)
*---------------------------------------------------------------------*
*     Purpose: resort I(kj,i;l) to I(jl;ki)
*
*     based on CC_SORT4O
*
*     Christian Neiss, 20/10-2005
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>
      INTEGER ISYM4O, IOPT
#if defined (SYS_CRAY)
      REAL XKJIL(*), XJLKI(*)
      REAL ZERO, ONE, HALF, DDOT, XNORM, FAC
#else
      DOUBLE PRECISION XKJIL(*), XJLKI(*)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM, FAC
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*
      INTEGER ISYML, ISYKJI, ISYMI, ISYMKJ, ISYJKL, ISYMJK
      INTEGER ISYMK, ISYMJ,  KJIL,  JLKI, IDXKI, ISYMJL, ISYMKI
*
*     --------------------------------------
*     Initialize result area with zero's or add to previous
*     --------------------------------------
      IF (IOPT.EQ.1) THEN
        FAC = ZERO
      ELSE
        FAC = ONE 
      END IF
*     --------------------------------------
*     Reorder thru loops on all 4 indices
*     --------------------------------------
      DO ISYML = 1,NSYM
         ISYKJI = MULD2H(ISYM4O,ISYML)
         DO L = 1, NRHF(ISYML)
            DO ISYMI = 1, NSYM
               ISYMKJ = MULD2H(ISYKJI,ISYMI)
               ISYJKL = MULD2H(ISYM4O,ISYMI)
               ISYMJK = MULD2H(ISYJKL,ISYML)
               DO I = 1, NRHF(ISYMI)
                  DO ISYMK = 1, NSYM
                     ISYMJ = MULD2H(ISYMJK,ISYMK)
                     ISYMJL = MULD2H(ISYMJ,ISYML)
                     ISYMKI = MULD2H(ISYMK,ISYMI)
                     DO K = 1, NRHF(ISYMK)
                     DO J = 1, NRHF(ISYMJ)

                KJIL = I3ORHF(ISYKJI,ISYML) + NMAIJK(ISYKJI)*(L-1)+
     &                 IMAIJK(ISYMKJ,ISYMI) + NMATIJ(ISYMKJ)*(I-1)+
     &                 IMATIJ(ISYMK,ISYMJ)  + NRHF(ISYMK)*(J-1) + K

                IDXKI = IMATIJ(ISYMK,ISYMI) + NRHF(ISYMK)*(I-1) + K

                JLKI = IGAMSQ(ISYMJL,ISYMKI) + 
     &                 NMATIJ(ISYMJL)*(IDXKI-1) +
     &                 IMATIJ(ISYMJ,ISYML) + NRHF(ISYMJ)*(L-1) + J

               XJLKI(JLKI) = FAC*XJLKI(JLKI) + XKJIL(KJIL)

                     END DO         !J
                     END DO         !K
                  END DO            !ISYMK
               END DO               !I
            END DO                  !ISYMI
         END DO                     !L
      END DO                        !ISYML
*     ---------------------------------------
*     Finished, return
*     ---------------------------------------
      RETURN
      END
*=====================================================================*
