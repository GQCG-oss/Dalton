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
C  /* Deck cc_bfif1 */
      SUBROUTINE CC_BFIF1(BFRHF0,ISYRHF0,XMGD1,ISYMGD1,
     &                   BFRHF1,ISYRHF1,XMGD2,ISYMGD2,
     &                   LUBFI,FNBFI,IADRBF,IADR,IDEL,WORK,LWORK)
*---------------------------------------------------------------------*
*
*     Purpose: contract effective density with (**|k delta) integrals
*              to the BF intermediate in the F matrix transformation
*              (special version of the CC_BFI routine for F matrix)
*
*     BFRHF(0 or 1)  : (**|k delta) integral distribution, sym ISYRHF(0,1)
*     XMGD(1 or 2)   : effective density for BF term, sym ISYMGD(1,2)
*
*     Generalization of CC_BFIF to calculate the BZ(QA) intermediate
*     Sonia and Poul, August 1999
*---------------------------------------------------------------------*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "dummy.h"
#include "ccorb.h"
#include "maxorb.h"
#include "ccsdsym.h"

      CHARACTER*(*) FNBFI
      INTEGER LWORK, ISYMGD1, ISYMGD2, ISYRHF0, ISYRHF1
      INTEGER IADR, IDEL, LUBFI, IADRBF(*)

#if defined (SYS_CRAY)
      REAL BFRHF0(*), BFRHF1(*), XMGD1(*), XMGD2(*), WORK(LWORK)
      REAL ZERO, ONE, HALF, DDOT, XNORM
#else
      DOUBLE PRECISION BFRHF0(*), BFRHF1(*)
      DOUBLE PRECISION XMGD1(*), XMGD2(*), WORK(LWORK)
      DOUBLE PRECISION ZERO, ONE, HALF, DDOT, XNORM
#endif
      PARAMETER(ZERO = 0.0D0, HALF = 0.5D0, ONE = 1.0D0)
C
      INTEGER NBFRHF(8), IBFRHF(8,8), ISYM, ICOUNT, ISYBET
      INTEGER ISYMIJ, ISYMIJB, KOFF1, KOFF2, KOFF3, NBASB, NTOTAK
      INTEGER ISYMAK, ISYMAK1, ISYMAK2
C
C     --------------------------------------
C     precalculate symmetry array for BFRHF:
C     --------------------------------------
      DO ISYM = 1, NSYM
        ICOUNT = 0
        DO ISYMAK = 1, NSYM
           ISYBET = MULD2H(ISYMAK,ISYM)
           IBFRHF(ISYMAK,ISYBET) = ICOUNT
           ICOUNT = ICOUNT + NT1AO(ISYMAK)*NBAS(ISYBET)
        END DO
        NBFRHF(ISYM) = ICOUNT
      END DO
C
      ISYMIJB = MULD2H(ISYRHF1,ISYMGD1)
      IF (ISYMIJB.NE.MULD2H(ISYRHF0,ISYMGD2)) THEN
        CALL QUIT('Symmetry mismatch 1 in CC_BFIF1.')
      END IF
      
C
      IF (LWORK .LT. ND2IJG(ISYMIJB)) THEN
         WRITE (LUPRI,*) 'LWORK =',LWORK
         WRITE (LUPRI,*) 'need ',ND2IJG(ISYMIJB)
         CALL QUIT('Insufficient memory in CC_BFIF1.')
      END IF
C
      DO ISYMIJ = 1, NSYM
C
         ISYMAK1 = MULD2H(ISYMGD1,ISYMIJ)
         ISYMAK2 = MULD2H(ISYMGD2,ISYMIJ)
         ISYBET  = MULD2H(ISYMAK1,ISYRHF1)
         IF (ISYBET.NE.MULD2H(ISYMAK2,ISYRHF0)) THEN
           CALL QUIT('Symmetry mismatch 2 in CC_BFIF1.')
         END IF
C
         KOFF1  = IT2AOIJ(ISYMAK1,ISYMIJ) + 1
         KOFF3  = IBFRHF(ISYMAK1,ISYBET)  + 1
         KOFF2  = ID2IJG(ISYMIJ,ISYBET)   + 1
         NTOTAK = MAX(NT1AO(ISYMAK1),1)
         NBASB  = MAX(NBAS(ISYBET),1)
C
         CALL DGEMM('T','N',NBAS(ISYBET),NMATIJ(ISYMIJ),NT1AO(ISYMAK1),
     &              ONE, BFRHF1(KOFF3),NTOTAK, XMGD1(KOFF1),NTOTAK,
     &              ZERO,WORK(KOFF2),NBASB)
C

         KOFF1  = IT2AOIJ(ISYMAK2,ISYMIJ) + 1
         KOFF3  = IBFRHF(ISYMAK2,ISYBET)  + 1
         NTOTAK = MAX(NT1AO(ISYMAK2),1)

         CALL DGEMM('T','N',NBAS(ISYBET),NMATIJ(ISYMIJ),NT1AO(ISYMAK2),
     &              ONE, BFRHF0(KOFF3),NTOTAK, XMGD2(KOFF1),NTOTAK,
     &              ONE,WORK(KOFF2),NBASB)

      END DO
C
c      WRITE(LUPRI,*) 'Test BZQAeta intermediate in CC_BFIF1'
c      XNORM = DDOT(ND2IJG(ISYMIJB), WORK, 1, WORK,1)
c      WRITE(LUPRI,*) 'Norm of BZQAeta', XNORM,' for Delta ', IDEL

      IADRBF(IDEL) = IADR
C
c     IF (LUBFI .LE. 0) THEN
c        CALL GPOPEN(LUBFI,FNBFI,'UNKNOWN',' ','UNFORMATTED',
c    *               IDUMMY,.FALSE.)
c     ENDIF
CCH
      IF (LUBFI.LE.0) THEN
        CALL QUIT('LUBFI=0 in CC_BFIF1!')
      END IF
CCH

      CALL PUTWA2(LUBFI,FNBFI,WORK,IADR,ND2IJG(ISYMIJB))
C
      IADR = IADR + ND2IJG(ISYMIJB)
C
c     CALL GPCLOSE(LUBFI,'KEEP')
C
      RETURN
      END
*=====================================================================*
