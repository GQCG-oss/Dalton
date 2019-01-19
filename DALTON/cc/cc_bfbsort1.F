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
C  /* Deck cc_bfbsort1 */
      SUBROUTINE CC_BFBSORT1(DSRHF,BSRHF,ISYRHF,LSQRAB)
*---------------------------------------------------------------------*
*
*     Purpose: presort DSRHF integral array for the BF intermediate
*              calculation in the B matrix transformation
*
*     DSRHF  : (alp bet|k delta) integrals for a fixed delta
*     BSRHF  : integrals sorted as I(alp k;bet)^del
*     ISYRHF : symmetry of the integral arrays DSRHF,BSRHF
*
*     Written by Christof Haettig July/October 1998
*     Updated by Sonia Coriani November 1999 to handle full (a b|
*---------------------------------------------------------------------*

      use dyn_iadrpk

#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"
#include "symsq.h"

      LOGICAL LSQRAB
      INTEGER ISYRHF, ISYM, ISYMAK, ISYBET, ISYMK, ISYMAB, ISYALP
      INTEGER ICOUNT, NBSRHF(8), IBSRHF(8,8)
      INTEGER NABK, NAKB, NAK, KOFF1, IJSQ

#if defined (SYS_CRAY)
      REAL DSRHF(*), BSRHF(*)
#else
      DOUBLE PRECISION DSRHF(*), BSRHF(*)
#endif
C
C     --------------------------------------
C     precalculate symmetry array for BSRHF:
C     --------------------------------------
C
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBSRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBSRHF(ISYM) = ICOUNT
      END DO
C
C     -------------------
C     sort the integrals:
C     -------------------
C
      DO ISYMAK = 1, NSYM
      DO ISYMK  = 1, NSYM
C 
         ISYBET = MULD2H(ISYMAK,ISYRHF)
         ISYALP = MULD2H(ISYMK,ISYMAK)
         ISYMAB = MULD2H(ISYALP,ISYBET)
C
C        --------------------------------------------------------
C        get (alp k;bet) blocks out of (alp bet|k del) integrals:
C        --------------------------------------------------------
C
         DO K = 1, NRHF(ISYMK)
C
            IF (LSQRAB) THEN
              KOFF1  = IDSRHFSQ(ISYMAB,ISYMK) + N2BST(ISYMAB)*(K-1)
            ELSE
              KOFF1  = IDSRHF(ISYMAB,ISYMK) + NNBST(ISYMAB)*(K-1) 
            END IF
C
            DO A = 1, NBAS(ISYALP)
            DO B = 1, NBAS(ISYBET)
C
               IJSQ = IAODIS(ISYALP,ISYBET) + NBAS(ISYALP)*(B-1) + A
               IF (LSQRAB) THEN
                 NABK = KOFF1  + IJSQ                        !not quite sure
               ELSE
                 NABK = KOFF1  + IADRPK( I2BST(ISYMAB) + IJSQ )
               END IF
               NAK  = IT1AO(ISYALP,ISYMK)   + NBAS(ISYALP)*(K-1) + A
               NAKB = IBSRHF(ISYMAK,ISYBET) +NT1AO(ISYMAK)*(B-1) + NAK
C
               BSRHF(NAKB) = DSRHF(NABK)
C
            END DO
            END DO
C
         END DO
C
      END DO
      END DO
C
      RETURN
      END
*=====================================================================*
