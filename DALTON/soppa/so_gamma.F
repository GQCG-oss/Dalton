C  /* Deck so_gamma */
      SUBROUTINE SO_GAMMA(DSRHF,GAMMA,XLAMDP,XLAMDH,SCRM,
     *                     WORK,LWORK,IDEL,ISYMD)
C
C     Original routine CCRHS_GAM is written by Henrik Koch 3-Jan-1994
C     Symmetry by Henrik Koch and Alfredo Sanchez. 21-July-1994
C
C     Purpose: Calculate the gamma intermediate.
C
C     List of variables:
C       DSRHF:  (alpha beta | j delta)-integral
C       GAMMA:  Result array
C       XLAMDP: 
C       XLAMDH: 
C       SCRM:   T_(ij)^(c delta) 
C       WORK:   
C       LWORK:  
C       IDEL:   
C       ISYMD:  
C
C     Written by Lilli Irene Ør Kristensen, Spring 2017
C
#include "implicit.h"
      DIMENSION DSRHF(*),GAMMA(*),SCRM(*)
      DIMENSION WORK(LWORK)
      DIMENSION XLAMDP(*),XLAMDH(*)
#include "priunit.h"
#include "ccorb.h"
#include "ccsdsym.h"
C
C------------------------
C     Dynamic allocation.
C------------------------
C
      KLAMDA = 1
      KEND1  = KLAMDA + NRHF(ISYMD)
      LWRK1  = LWORK  - KEND1
C
      IF (LWRK1 .LT. 0) THEN
         WRITE(LUPRI,*) 'Need : ',KEND1,'Available : ',LWORK
         CALL QUIT('Insufficient space in SO_GAMMA (1) ')
      ENDIF
C
C---------------------------------------
C     Copy XLAMDH vector for given IDEL.
C---------------------------------------
C
      KOFF1 = ILMRHF(ISYMD) + IDEL - IBAS(ISYMD)
      CALL DCOPY(NRHF(ISYMD),XLAMDH(KOFF1),NBAS(ISYMD),WORK(KLAMDA),1)
C
C--------------------------------
C     Calculate the contribution.
C--------------------------------
C
      ISYDIS = MULD2H(ISYMD,ISYMOP)
C
      DO 100 ISYML = 1,NSYM
C
         ISYMAG = MULD2H(ISYML,ISYDIS)
C
C---------------------------
C        Dynamic allocation.
C---------------------------
C
         KSCR1  = KEND1
         KSCR2  = KSCR1  + N2BST(ISYMAG)
C         KSCR3  = KSCR2  + NT1AO(ISYMAG)
C         KSCR4  = KSCR3  + NT1AM(ISYMAG)
         KSCR3  = KSCR2  + NT1AO(ISYMAG)
C         KEND2  = KSCR4  + NMATIJ(ISYMAG)
         KEND2  = KSCR3  + NMATIJ(ISYMAG)
         LWRK2  = LWORK  - KEND2
C
         IF (LWRK2 .LT. 0) THEN
            WRITE(LUPRI,*) 'Need : ',KEND2,'Available : ',LWORK
            CALL QUIT('Insufficient space in SO_GAMMA (2)')
         ENDIF
C
         CALL SO_GAMMA1(DSRHF,GAMMA,WORK(KLAMDA),XLAMDP,XLAMDH,
     *                   WORK(KSCR1),WORK(KSCR2),WORK(KSCR3),
     *                  WORK(KEND2),LWRK2,ISYMD,ISYML,
     *                   ISYMAG)
C
  100 CONTINUE
C
      RETURN
      END
