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
      SUBROUTINE CC_LAMBDAQ(XLAMDPQ,XLAMDHQ,CMOPQ,CMOHQ,ISYMQ,T1AMP0, 
     &                      XKAPPA, RMAT, IREAL, IOPT, WORK, LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: Calculate the derivative LambdaQ matrices which are
*              build from a derivative MO vector, accounting for
*              orbital relaxation and reorthogonalization (the 
*              contribution from the connection matrix)
*
*     Unfinished!!! : needs some fixes for frozen/deleted orbitals
*
*              IOPT = 0 :
*                 CMOPQ   = C (R - kappa)^*            (Sirius order)
*                 CMOHQ   = C (R - kappa)              (Sirius order)
*
*              IOPT = 1 :
*                 CMOQ    = C (R - kappa)                (CC order)
*                 XLAMDPQ = C (R - kappa)^* (1 - T1^T)   (CC order)
*                 XLAMDHQ = C (R - kappa) (1 + T1  )     (CC order)
*
*
*              IREAL = +1  :  R and kappa are real
*                      -1  :  R and kappa are pure imaginary
*              
*              Symmetries:     ISYMQ  -- LambdaQ, R, kappa 
*                              1      -- C, T1 
*
*
*     Christof Haettig, spring 1999
*     generalized for CMOPQ different from CMOHQ, november 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "dummy.h"      
#include "ccorb.h"
#include "ccfro.h"
#include "ccsdsym.h"
#include "inftap.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER IOPT, ISYMQ, IREAL, LWORK

#if defined (SYS_CRAY)
      REAL ONE, ZERO
      REAL XLAMDPQ(*), XLAMDHQ(*), CMOPQ(*), CMOHQ(*), T1AMP0(*)
      REAL XKAPPA(*), RMAT(*), WORK(LWORK) 
#else
      DOUBLE PRECISION ONE, ZERO
      DOUBLE PRECISION XLAMDPQ(*),XLAMDHQ(*),CMOPQ(*),CMOHQ(*),T1AMP0(*)
      DOUBLE PRECISION XKAPPA(*), RMAT(*), WORK(LWORK) 
#endif
      PARAMETER(ONE=1.0D0, ZERO=0.0D0)

      LOGICAL NOKAPPA
      INTEGER NCMO(8), ICMO(8,8)
      INTEGER KCMO, KQMATP, KQMATH, KCMOPQ, KCMOHQ
      INTEGER KSCR, KEND1, LWRK1, NBASA, NBASB
      INTEGER ISYM, ISYM1, ISYM2, ICOUNT, ISYALP, ISYBET, NORBSA
      INTEGER KOFF1, KOFF2, KOFF3
 
*---------------------------------------------------------------------*
*     print some debug output:
*---------------------------------------------------------------------*
      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'entered CC_LAMBDAQ...'
        WRITE (LUPRI,*) 'connection matrix in AO:'
        CALL CC_PRONELAO(RMAT,ISYMQ)
      END IF

*---------------------------------------------------------------------*
*     memory allocation and some setup:
*---------------------------------------------------------------------*
      KCMO    = 1
      KQMATP  = KCMO   + N2BST(ISYM0)
      KQMATH  = KQMATP + N2BST(ISYMQ)
      KCMOPQ  = KQMATH + N2BST(ISYMQ)
      KCMOHQ  = KCMOPQ + N2BST(ISYMQ)
      KSCR    = KCMOHQ + N2BST(ISYMQ)
      KEND1   = KSCR   + N2BST(ISYMQ)
      LWRK1   = LWORK  - KEND1

      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_LAMBDAQ.')
      END IF
 
      DO ISYM = 1, NSYM
         ICOUNT = 0
         DO ISYM2 = 1, NSYM
            ISYM1 = MULD2H(ISYM,ISYM2)
            ICMO(ISYM1,ISYM2) = ICOUNT
            ICOUNT = ICOUNT + NBAS(ISYM1)*NORBS(ISYM2)
         END DO
         NCMO(ISYM) = ICOUNT
      END DO

*---------------------------------------------------------------------*
*     read (undifferentiated) MO coefficients from file:
*---------------------------------------------------------------------*
      CALL GPOPEN(LUSIFC,'SIRIFC','OLD',' ','UNFORMATTED',IDUMMY,
     &            .FALSE.)
      REWIND(LUSIFC)
      CALL MOLLAB('TRCCINT ',LUSIFC,LUPRI)
      READ(LUSIFC)
      READ(LUSIFC)
      READ(LUSIFC) (WORK(KCMO+I-1),I=1,NLAMDS)
      CALL GPCLOSE(LUSIFC,'KEEP')

*---------------------------------------------------------------------*
*     get Q matrix in MO representation:
*---------------------------------------------------------------------*
      NOKAPPA = .FALSE.
      CALL CC_QMAT(WORK(KQMATP),WORK(KQMATH),RMAT,XKAPPA,
     &             IREAL,ISYMQ,NOKAPPA,WORK(KCMO),WORK(KEND1),LWRK1)

      IF (LOCDBG) THEN
        WRITE (LUPRI,*) 'Q^h matrix in MO:'
        CALL CC_PRONELAO(WORK(KQMATH),ISYMQ)
      END IF

*---------------------------------------------------------------------*
*     transform leading index to contravariant AO basis:
*             CMOQ^h = CMO x Q^h;    CMOQ^p = CMO x Q^p
*---------------------------------------------------------------------*
      DO ISYALP = 1, NSYM
         ISYBET = MULD2H(ISYALP,ISYMQ)

         NBASA  = MAX(NBAS(ISYALP),1)
         NORBSA = MAX(NORBS(ISYALP),1)

         KOFF1 = KCMO   + ICMO(ISYALP,ISYALP)
         KOFF2 = KQMATH + IAODIS(ISYALP,ISYBET)
         KOFF3 = KCMOHQ + ICMO(ISYALP,ISYBET)

         CALL DGEMM('N','N',NBAS(ISYALP),NORBS(ISYBET),NORBS(ISYALP),
     &              ONE,WORK(KOFF1),NBASA,WORK(KOFF2),NORBSA,
     &              ZERO,WORK(KOFF3),NBASA)

         KOFF2 = KQMATP + IAODIS(ISYALP,ISYBET)
         KOFF3 = KCMOPQ + ICMO(ISYALP,ISYBET)

         CALL DGEMM('N','N',NBAS(ISYALP),NORBS(ISYBET),NORBS(ISYALP),
     &              ONE,WORK(KOFF1),NBASA,WORK(KOFF2),NORBSA,
     &              ZERO,WORK(KOFF3),NBASA)

      END DO

*---------------------------------------------------------------------*
*     reorder to CC standard ordering and calculate XLAMBDAQ matrices:
*---------------------------------------------------------------------*
      IF      (IOPT.EQ.0) THEN

         CALL DCOPY(NGLMDT(ISYMQ),WORK(KCMOPQ),1,CMOPQ,1)
         CALL DCOPY(NGLMDT(ISYMQ),WORK(KCMOHQ),1,CMOHQ,1)

      ELSE IF (IOPT.EQ.1) THEN

         CALL CMO_REORDER2(WORK(KCMOPQ),CMOPQ,ISYMQ)
         CALL CMO_REORDER2(WORK(KCMOHQ),CMOHQ,ISYMQ)

         CALL DCOPY(NGLMDT(ISYMQ),CMOPQ,1,XLAMDPQ,1)
         CALL DCOPY(NGLMDT(ISYMQ),CMOHQ,1,XLAMDHQ,1)

         CALL LAMDA2(XLAMDPQ,XLAMDHQ,ISYMQ,T1AMP0,ISYM0,
     &               CMOPQ,CMOHQ,ISYMQ)

      ELSE 
         CALL QUIT('Illegal option in CC_LAMBDAQ.')
      END IF

      RETURN
      END
*=====================================================================*
