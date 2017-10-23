C  /* Deck so_aikdelta */
      SUBROUTINE SO_AIKDELTA(DSRHF,ISYMT2,XLAMDH,XLAMPC,ISYMPC,ISYMHC,
     *                   WORK,LWORK,IDEL,ISYMD,FACTD,LUD,DFIL)
C
C     The original routine CCRHS_D is written by Henrik Koch 9-Jan-1994.
C
C     Purpose: Set up the preliminary memory allocation needed to
C              calculate (ai|k delta) (is done in SO_AIKDELTA1),
C              where the loop over delta is outside this routine.
C
C     List of variables:
C       DSRHF:   one-index transfromed-integrals (alpha gamma|j delta)
C       ISYMT2: Is ISYMTR (sym. of trialvector = 1)
C       XLAMDP: 
C       XLAMDH: 
C       XLAMPC: 
C       ISYMPC: 
C       IDEL:   
C       ISYMD:  
C       FACTC:  HALF
C       LUD:
C       DFIL:
C     Written by Lilli Irene Ør Kristensen, spring 2017.
C
#include "implicit.h"
      PARAMETER (ONE = 1.0D0, TWO = 2.0D0)
      DIMENSION DSRHF(*),WORK(LWORK)
      DIMENSION XLAMDH(*)
      DIMENSION XLAMPC(*)
      CHARACTER DFIL*(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
      IF (ISYMT2 .NE. ISYMPC) CALL QUIT('Symmetry Mismatch in CCRHS_D' )
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KSCR1  = 1
      KSCR2  = KSCR1  + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
C      KSCR3  = KSCR2  + NT2BGD(ISYDIS)
C      KEND1  = KSCR3  + NT2BGD(ISYMD)
      KEND1  =KSCR2 + NT2BGD(ISYDIS)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in CCRHS_D')
      ENDIF
C
C--------------------------------
C     Calculate (ai|k delta)
C--------------------------------
C
      CALL SO_AIKDELTA1(DSRHF,WORK(KSCR1),WORK(KSCR2),
     *              XLAMDH,XLAMPC,ISYMPC,ISYMHC,
     *              WORK(KEND1),LWRK1,ISYDIS,IDEL,
     *              ISYMD,LUD,DFIL)
C
      RETURN
      END
