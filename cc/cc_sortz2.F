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
C  /* Deck cc_sortz2 */
      SUBROUTINE CC_SORTZ2(ZKJAI,ZJKIA,ISYZA2,IOPT)
*---------------------------------------------------------------------*
*     Purpose: resort ZA2(kj,ai) to ZA2(jk,i;a)
*
*     IOPT = 1 : ZJKIA area is initialized here
*     IOPT = 2 : ZJKIA is added to intermediate
*                       already stored as ZA2(jk,l;a)
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
      INTEGER ISYZA2, IOPT
#if defined (SYS_CRAY)
      REAL ZKJAI(*), ZJKIA(*)
      REAL ZERO, ONE, HALF, DDOT, XNORM, FAC
#else
      DOUBLE PRECISION ZKJAI(*), ZJKIA(*)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM, FAC
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
*
      INTEGER ISYKJI, ISYMI, ISYMKJ, ISYJKI, ISYMJK, ISYM
      INTEGER ISYMK, ISYMJ,  KJIA,  JKLI, ISYMA, AI, JKIA,ISYMAI
      INTEGER ICOUNT, ISKJAI(8,8), NSKJAI(8), ISYKJA, KJAI
*     --------------------------------------
*     precalculate symmetry array for ZKJAI:
*     --------------------------------------
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAI = 1, NSYM
           ISYMKJ = MULD2H(ISYMAI,ISYM)
           ISKJAI(ISYMKJ,ISYMAI) = ICOUNT
           ICOUNT = ICOUNT + NMATIJ(ISYMKJ)*NT1AM(ISYMAI)
        END DO
        NSKJAI(ISYM) = ICOUNT
      END DO
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
      DO ISYMI = 1,NSYM
         ISYKJA = MULD2H(ISYZA2,ISYMI)
         DO I = 1, NRHF(ISYMI)
            DO ISYMA = 1, NSYM
               ISYMAI = MULD2H(ISYMA,ISYMI)
               ISYMKJ = MULD2H(ISYKJA,ISYMA)
               ISYJKI = MULD2H(ISYZA2,ISYMA)
               ISYMJK = MULD2H(ISYJKI,ISYMI)
               DO A = 1, NVIR(ISYMA)
                  DO ISYMJ = 1, NSYM
                     ISYMK = MULD2H(ISYMKJ,ISYMJ)
                     DO J = 1, NRHF(ISYMJ)
                     DO K = 1, NRHF(ISYMK)

                AI   = IT1AM(ISYMA,ISYMI) + NVIR(ISYMA)*(I-1) + A
                KJAI = ISKJAI(ISYMKJ,ISYMAI) + NMATIJ(ISYMKJ)*(AI-1)+
     &                 IMATIJ(ISYMK,ISYMJ)   + NRHF(ISYMK)*(J-1) + K

                JKIA = I3OVIR(ISYJKI,ISYMA) + NMAIJK(ISYJKI)*(A-1)+
     &                 IMAIJK(ISYMJK,ISYMI) + NMATIJ(ISYMJK)*(I-1)+
     &                 IMATIJ(ISYMJ,ISYMK)  + NRHF(ISYMJ)*(K-1) + J

               ZJKIA(JKIA) = FAC*ZJKIA(JKIA) + ZKJAI(KJAI)

                     END DO         !K
                     END DO         !J
                  END DO            !ISYMJ
               END DO               !A
            END DO                  !ISYMA
         END DO                     !I
      END DO                        !ISYMI
*     ---------------------------------------
*     Finished, return
*     ---------------------------------------
      RETURN
      END
*=====================================================================*
