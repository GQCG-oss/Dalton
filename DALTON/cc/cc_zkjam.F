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
#  include "implicit.h"
#endif
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"

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
