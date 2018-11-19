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
      SUBROUTINE CC_GET_LAMBDA1(IRELAX, IOPER, IOPT, XLAMDPQ, XLAMDHQ, 
     &                          CMOPQ,CMOHQ,ISYMQ,T1AMP0,WORK,LWORK )
*---------------------------------------------------------------------*
*
*     Purpose: Calculate the derivative LambdaQ matrices 
*              (for more information see CC_LAMBDAQ routine)
*
*     Christof Haettig, summer 1999
*
*=====================================================================*
#if defined (IMPLICIT_NONE)
      IMPLICIT NONE
#else
#  include "implicit.h"
#endif
#include "priunit.h"
#include "ccorb.h"
#include "ccfro.h"
#include "ccsdsym.h"
#include "ccroper.h"

      LOGICAL LOCDBG
      PARAMETER (LOCDBG = .FALSE.)

      INTEGER ISYM0
      PARAMETER (ISYM0 = 1)

      INTEGER ISYMQ, LWORK, IOPER, IRELAX, IOPT

#if defined (SYS_CRAY)
      REAL XLAMDPQ(*), XLAMDHQ(*), CMOPQ(*), CMOHQ(*), T1AMP0(*)
      REAL WORK(LWORK) 
#else
      DOUBLE PRECISION XLAMDPQ(*),XLAMDHQ(*),CMOPQ(*),CMOHQ(*),T1AMP0(*)
      DOUBLE PRECISION WORK(LWORK) 
#endif

      CHARACTER*(10) MODEL
      INTEGER KAPPA, KRMAT, KEND, LWRK, IREAL

*---------------------------------------------------------------------*
*     memory allocation:
*---------------------------------------------------------------------*
      KAPPA = 1
      KRMAT = KAPPA + 2*NALLAI(ISYMQ)
      KEND  = KRMAT + N2BST(ISYMQ)
      LWRK  = LWORK - KEND

      IF (LWRK .LT. 0) THEN
         CALL QUIT('Insufficient work space in CC_GET_LAMBDA1.')
      END IF

      IF (LBLOPR(IOPER).EQ.'HAM0    ' .AND. IRELAX.GE.1) THEN
         WRITE (LUPRI,*) 'Test case "HAM0"... no relaxation '//
     &        'vector used.'
         CALL DZERO(WORK(KAPPA),2*NALLAI(ISYMQ))
      ELSE IF (IRELAX.GE.1) THEN
         CALL CC_RDHFRSP('R1 ',IRELAX,ISYMQ,WORK(KAPPA))
      ELSE
         CALL DZERO(WORK(KAPPA),2*NALLAI(ISYMQ))
      END IF

      CALL CC_GET_RMAT(WORK(KRMAT),IOPER,1,ISYMQ,WORK(KEND),LWRK)

      IOPT  = 1
      IREAL = ISYMAT(IOPER)
      CALL CC_LAMBDAQ(XLAMDPQ,XLAMDHQ,CMOPQ,CMOHQ,ISYMQ,
     &                T1AMP0,WORK(KAPPA),WORK(KRMAT),
     &                IREAL,IOPT,WORK(KEND),LWRK)
      

      RETURN
      END
*=====================================================================*
