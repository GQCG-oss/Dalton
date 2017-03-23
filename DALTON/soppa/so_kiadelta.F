C  /* Deck so_kiadelta */
      SUBROUTINE SO_KIADELTA(XINT,ISYMT2,XLAMDP,XLAMDH,XLAMPC,
     *                   ISYMPC,WORK,LWORK,IDEL,ISYMD,FACTC,LUC,CFIL)
C
C     Original routine CCRHS_C was written by Henrik Koch 3-Jan-1994.
C
C     Purpose: Set up the preliminary memory allocation needed to
C               calculate (ki|a delta) (is done in SO_KIADELTA1),
C               where the loop over delta is outsidei this routine.
C
C     List of variables:
C       XINT:   AO-integrals
C       ISYMT2: Is ISYMTR (sym. of trialvector = 1)
C       XLAMDP: 
C       XLAMDH: 
C       XLAMPC: 
C       ISYMPC: 
C       IDEL:   
C       ISYMD:  
C       FACTC:  HALF
C       LUC:
C       CFIL:
C
C     Written by Lilli Irene Ør Kristensen, spring 2017.
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
      DIMENSION XINT(*),WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMPC(*),XLAMDH(*)
      CHARACTER CFIL*(*),CPFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
#include "ccsdinp.h"
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
C--------------------------------------
C     Dynamic allocation of work space.
C--------------------------------------
C
      KSCR1 = 1
      KSCR2 = KSCR1 + MAX(NT2BCD(ISYAIK),NT2BCD(ISYDIS))
      KSCR3 = KSCR2 + NT2BGD(ISYDIS)
      IF (ICON .EQ. 2) THEN
         KEND1  = KSCR3  + NT2BGD(ISYMD)
      ELSE
         KEND1  = KSCR3  + NT2BGD(ISYAIK)
      ENDIF

      LWRK1 = LWORK - KEND1
      IF (LWRK1 .LT. 0) THEN
         CALL QUIT('Insufficient space for allocation in CCRHS_C')
      END IF
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
         CALL SO_KIADELTA1(XINT,ISYMT2,WORK(KSCR1),WORK(KSCR2),
     *                 XLAMDP,XLAMDH,XLAMPC,ISYMPC
     *                 WORK(KEND1),LWRK1,
     *                 ISYDIS,IDEL,FACTC,LUC,CFIL)
C
C
      RETURN
      END
