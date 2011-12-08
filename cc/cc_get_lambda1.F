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
#  include <implicit.h>
#endif
#include <priunit.h>
#include <ccorb.h>
#include <ccfro.h>
#include <ccsdsym.h>
#include <ccroper.h>

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
