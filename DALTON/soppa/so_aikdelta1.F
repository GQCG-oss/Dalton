C    /* Deck so_aikdelta1 */
      SUBROUTINE SO_AIKDELTA1(DSRHF,SCR1,SCR2,XLAMDH,XLAMPC,ISYMPC,
     *                    ISYMHC,WORK,LWORK,ISYDIS,IDEL,ISYDEL,FACTD,
     *                    LUD,DFIL)
C
C     The orignal routine CCRHS_D1 is written by Henrik Koch 3-Jan-1994.
C
C     This routine transfrom the AO-integral (alpha gamma|beta delta) to
C       (a i | k delta).
C     In ROUTINE NAME (ai|k delta) is multiplied with the 2nd order
C       doubles amplitudes.
C
C     LIST OF VARIALBES:
C      DSRHF:  Contains (alpha beta| j delta)
C      SCR1:   Here is (ai|k delta) is calculated
C      SCR2:   Used during calculations
C      XLAMDH: 
C      XLAMPC: XLAMDP 
C      ISYMPC: ISYMTR
C      ISYMHC: ISYMTR
C      ISYDIS: In soppa this is ISYMD
C      IDEL:   ?
C      ISYDEL: Is ISYMD, symmetry of delta 
C      FACTD:  HALF
C      LUD:    LUNIT?
C      DFIL:   The File name of where we save the integral.
C
C      Written by Lilli Irene Ør Kristensen, spring 2017
C
#include "implicit.h"
#include "priunit.h"
#include "maxorb.h"
#include "ccsdinp.h"
      PARAMETER(ZERO=0.0D0,ONE=1.0D0,HALF=0.5D0,XMHALF=-0.5D0)
      PARAMETER(TWO=2.0D0)
      DIMENSION DSRHF(*)
      DIMENSION SCR1(*),SCR2(*)
      DIMENSION XLAMDH(*)
      DIMENSION XLAMPC(*)
      DIMENSION WORK(LWORK)
      CHARACTER DFIL*(*)
#include "ccorb.h"
#include "symsq.h"
#include "ccsdsym.h"
#include "ccsdio.h"
C
      INDEX(I,J) = MAX(I,J)*(MAX(I,J)-3)/2 + I + J
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      DO 340 ISYMK = 1,NSYM
C
         ISYALG = MULD2H(ISYMK,ISYDIS)
         ISYALI = MULD2H(ISYMHC,ISYALG)
         NT1AOM = MAX(NT1AO(ISYALG),NT1AO(ISYALI))
C
         KSCR10 = 1
         KSCR11 = KSCR10 + N2BST(ISYALG)
         KEND1  = KSCR11 + NT1AOM
         LWRK1  = LWORK  - KEND1
         IF (LWRK1 .LT. 0) THEN
            CALL QUIT('Insufficient space for allocation in CCRHS_D1')
         END IF
C
         DO 350 K = 1,NRHF(ISYMK)
C
            KOFF1 = IDSRHF(ISYALG,ISYMK) + NNBST(ISYALG)*(K - 1) + 1
            CALL CCSD_SYMSQ(DSRHF(KOFF1),ISYALG,WORK(KSCR10))
C
            ISYALI = ISYALG
            CALL DZERO(WORK(KSCR11),NT1AO(ISYALI))
C
C-------------------------------------------------------------------
C           Transform the gamma index in (alpha gamma| j delta) to i.
C--------------------------------------------------------------------
C
            DO 360 ISYMI = 1,NSYM
C
               ISYMAL = MULD2H(ISYMI,ISYALI)
               ISYMG  = ISYMI
C
               NBASAL = MAX(NBAS(ISYMAL),1)
               NBASG = MAX(NBAS(ISYMG),1)
C
               KOFF2 = KSCR10 + IAODIS(ISYMAL,ISYMG)
               KOFF3 = ILMRHF(ISYMI) + 1
               KOFF4 = KSCR11 + IT1AO(ISYMAL,ISYMI)
C
               CALL DGEMM('N','N',NBAS(ISYMAL),NRHF(ISYMI),NBAS(ISYMG),
     *                    ONE,WORK(KOFF2),NBASAL,XLAMDH(KOFF3),NBASG,
     *                    ZERO,WORK(KOFF4),NBASAL)
C
  360       CONTINUE
C
            NRHFK = MAX(NRHF(ISYMK),1)
            KOFF5 = IT2BGT(ISYMK,ISYALI) + K
C            CALL DAXPY(NT1AO(ISYALI),TWO,WORK(KSCR11),1,SCR2(KOFF5),
C     *                 NRHFK)
            CALL DCOPY(NT1AO(ISYALI),WORK(KSCR11),1,SCR2(KOFF5),NRHFK)
C               CHANGE THIS TO DCOPY
C
  350    CONTINUE
C
  340 CONTINUE
C
C--------------------------------------------------------------------------
C     Transform the alpha index of (alpha i| k delta) to a.
C--------------------------------------------------------------------------
C
      ISYAIK = MULD2H(ISYDIS,ISYMPC)
C
      DO 710 ISYMAI = 1,NSYM
C
         ISYMK = MULD2H(ISYMAI,ISYAIK)
         NRHFK = MAX(NRHF(ISYMK),1)
C
         DO 720 ISYMI = 1,NSYM
C
            ISYMA  = MULD2H(ISYMI,ISYMAI)
            ISYMAL = MULD2H(ISYMPC,ISYMA)
            ISYALI = MULD2H(ISYMAL,ISYMI)
            NBASAL = MAX(NBAS(ISYMAL),1)
C
            DO 730 I = 1,NRHF(ISYMI)
C
               NAI   = IT1AM(ISYMA,ISYMI)   + NVIR(ISYMA)*(I - 1) + 1
               MALI  = IT1AO(ISYMAL,ISYMI)  + NBAS(ISYMAL)*(I - 1) + 1
C
               KOFF1 = IT2BGT(ISYMK,ISYALI) + NRHF(ISYMK)*(MALI - 1) + 1
               KOFF2 = IGLMVI(ISYMAL,ISYMA) + 1
               KOFF3 = IT2BCT(ISYMK,ISYMAI) + NRHF(ISYMK)*(NAI - 1) + 1
C
               CALL DGEMM('N','N',NRHF(ISYMK),NVIR(ISYMA),NBAS(ISYMAL),
     *                    ONE,SCR2(KOFF1),NRHFK,XLAMPC(KOFF2),NBASAL,
     *                    FACTD ,SCR1(KOFF3),NRHFK)
C
  730       CONTINUE
  720    CONTINUE
  710 CONTINUE
C
C---------------------------------------
C     Dump to Disk.
C---------------------------------------
C
      IOFF = IT2DEL(IDEL) + 1
C
      CALL PUTWA2(LUD,DFIL,SCR1,IOFF,NT2BCD(ISYAIK))
C
      RETURN
      END
