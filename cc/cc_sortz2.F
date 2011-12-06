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
