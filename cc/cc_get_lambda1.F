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
#  include <priunit_dec.h>
#  include <ccorb_dec.h>
#  include <ccsdsym_dec.h>
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
         IOPT = 1
         CALL CC_RDHFRSP('R1',IRELAX,ISYMQ,IOPT,MODEL,WORK(KAPPA))
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
