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
C  /* Deck cc_zkjai */
      SUBROUTINE CC_ZKJAM(CTR2,ISYCTR,TAMP,ISYTAM,ZKJAI)
*---------------------------------------------------------------------*
*
*     Purpose: transform the b index of a Zeta_bj,ai (CTR2)
*              to k by contraction with a T_bk (TAMP) type of 
*              vector --> Zeta_kj,ai
*
*     Sonia Coriani, 10/09-1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include <implicit.h>
#endif
#include <ccorb.h>
#include <maxorb.h>
#include <ccsdsym.h>

      INTEGER ISYCTR,ISYTAM

#if defined (SYS_CRAY)
      REAL CTR2(*), TAMP(*), ZKJAM(*)
      REAL ZERO, ONE, HALF, DDOT, XNORM
#else
      DOUBLE PRECISION CTR2(*), TAMP(*), ZKJAI(*)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
C
      INTEGER AI, ISYMAI, ISYMBJ, ISYMB, ISYMJ, ISYMK, ISYMKJ
      INTEGER ICOUNT,ISYM,KOFFZ,KOFFT,KOFFR,ISYRES
      INTEGER ISKJAI(8,8), NSKJAI(8), NVIRB, NRHFK
C
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
*
      ISYRES = MULD2H(ISYCTR,ISYTAM)

*     ---------------------------------------
*     Calculate Z_kj,ai = sum_b t_bk Z_bj,ai
*     ---------------------------------------

      DO ISYMAI = 1, NSYM
         ISYMBJ = MULD2H(ISYCTR,ISYMAI)
         DO ISYMB = 1, NSYM
            ISYMJ = MULD2H(ISYMBJ,ISYMB)
            ISYMK = MULD2H(ISYMB,ISYTAM)
            ISYMKJ = MULD2H(ISYMK,ISYMJ)
            DO AI = 1, NT1AM(ISYMAI)
               KOFFZ = IT2SQ(ISYMBJ,ISYMAI) + NT1AM(ISYMBJ)*(AI-1) +
     &                 IT1AM(ISYMB,ISYMJ)   + 1
               KOFFT = IT1AM(ISYMB,ISYMK) + 1
               KOFFR = ISKJAI(ISYMKJ,ISYMAI) + NMATIJ(ISYMKJ)*(AI-1) +
     &                 IMATIJ(ISYMK,ISYMJ)  + 1
 
               NVIRB = MAX(NVIR(ISYMB),1)
               NRHFK = MAX(NRHF(ISYMK),1)

               CALL DGEMM('T','N',NRHF(ISYMK),NRHF(ISYMJ),NVIR(ISYMB),
     &                     ONE,TAMP(KOFFT),NVIRB,CTR2(KOFFZ),NVIRB,
     &                     ZERO,ZKJAI(KOFFR),NRHFK)
            END DO  !AI
         END DO     !ISYMB
      END DO        !ISYMAI
C
      RETURN
      END
*=====================================================================*
