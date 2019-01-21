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
      SUBROUTINE CCRHS_IPM1(XINT,XINTP,XINTM,SCRAB,INDV1,INDV2,
     *                     ISYMAB,ISYMG,NUMG,IG1,IG2,IOPT)
C
C     Written by Henrik Koch 17-aug-1994.
C
C     Purpose: Making plus and minus combination of integrals.
C              (a>=g|bd) -> K+ and K- where
C                           K+- = (ag|bd) +- (bg|ad) a<=b,g<=d
C
C     Modified by Sonia Coriani 26-oct-1999 to 
C              handle XINT with squared (ag|bd) (ag part): 
C              If IOPT = 0, XINT in input is packed  (a>=g|bd)
C              If IOPT = 1, XINT in input is squared (ag|bd)
C        
C
#include "implicit.h"
#include "maxorb.h"
      PARAMETER(ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION XINT(*),XINTP(*),XINTM(*),SCRAB(*)
      DIMENSION INDV1(*), INDV2(*)
#include "ccorb.h"
#include "ccsdsym.h"
#include "symsq.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYDIS = MULD2H(ISYMAB,ISYMG)
C
C
      DO 100 G = IG1,IG2
C
         IG = G - IG1 + 1
C
         DO 110 ISYMB = 1,NSYM
C
            ISYMA  = MULD2H(ISYMB,ISYMAB)
            ISYMAG = MULD2H(ISYMA,ISYMG)
C
            NTOTA  = MAX(NBAS(ISYMA),1)
cs
            IF (IOPT.EQ.0) THEN
               NTOTAG = MAX(NNBST(ISYMAG),1)
            ELSE IF (IOPT.EQ.1) THEN
               NTOTAG = MAX(N2BST(ISYMAG),1)
            ELSE
               CALL QUIT('Unknown option in CCRHS_IPM1')
            END IF
C
            DO 120 A = 1,NBAS(ISYMA)
C
               IF (IOPT.EQ.0) THEN
                  IF (ISYMA .EQ. ISYMG) THEN
                     KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                     + INDEX(G,A)
                  ELSE IF (ISYMA .LT. ISYMG) THEN
                     KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                     + NBAS(ISYMA)*(G - 1) + A
                  ELSE
                     KOFF1 = IDSAOG(ISYMB,ISYDIS) + IAODPK(ISYMA,ISYMG)
     *                     + NBAS(ISYMG)*(A - 1) + G
                  ENDIF
               ELSE IF (IOPT.EQ.1) THEN
                  KOFF1 = IDSAOGSQ(ISYMB,ISYDIS) + IAODIS(ISYMA,ISYMG)
     &                  + NBAS(ISYMA)*(G - 1) + A
               END IF
C
               KOFF2 = IAODIS(ISYMA,ISYMB) + A
C
               CALL DCOPY(NBAS(ISYMB),XINT(KOFF1),NTOTAG,
     *                    SCRAB(KOFF2),NTOTA)
C
  120       CONTINUE
C
  110    CONTINUE
C
         KOFF = NNBST(ISYMAB)*(IG - 1)
C
#if !defined (SYS_CRAY)
         DO 130 I = 1,NNBST(ISYMAB)
C
            XINTP(KOFF + I) = SCRAB(INDV1(I))
            XINTM(KOFF + I) = SCRAB(INDV2(I))
C
  130    CONTINUE
#else
         CALL GATHER(NNBST(ISYMAB),XINTP(KOFF + 1),SCRAB,INDV1)
         CALL GATHER(NNBST(ISYMAB),XINTM(KOFF + 1),SCRAB,INDV2)
#endif
C
  100 CONTINUE
C
C
      NTOT = NNBST(ISYMAB)*NUMG
C
      CALL DAXPY(NTOT,ONE,XINTM,1,XINTP,1)
      CALL DSCAL(NTOT,-TWO,XINTM,1)
      CALL DAXPY(NTOT,ONE,XINTP,1,XINTM,1)
C
      RETURN
      END
