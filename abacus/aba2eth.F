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
C  /* Deck u2driv */
      SUBROUTINE U2DRIV(AOINT,HCINT,HCRNT,
     &                  COEF34,INDHER,CONT3,CONT4,WORK,
     &                  LWORK,NPCO3,NPCO4,NUCS34,IPRINT,LMNV12,LMNV34,
     &                  IODD12,IODD34,NPNT34,NRED34,R12EIN)
C     Copy of C2DRIV for [r12,T1+T2] integrals (WK/UniKA/15-11-2002).
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
      LOGICAL R12EIN
      DIMENSION AOINT(*), HCINT(*), WORK(LWORK), COEF34(*), CONT3(*),
     &          CONT4(*), NPCO3(*), NPCO4(*), NUCS34(*), INDHER(*),
     &          LMNV12(*), LMNV34(*), IODD12(*), IODD34(*),
     &          NPNT34(*), NRED34(*), HCRNT(*)
#include <twoao.h>
#include <twosta.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from U2DRIV','*',103)
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 1010) NHKT3, NHKT4
         WRITE (LUPRI, 1020) KCKT3, KCKT4
         WRITE (LUPRI, 1025) KHKT3, KHKT4
         WRITE (LUPRI, 1030) KCKT34
         WRITE (LUPRI, 1040) DIAG34
         WRITE (LUPRI, 1060) NUC34
         WRITE (LUPRI, 1070) I340X, I340Y, I340Z
         WRITE (LUPRI, 1005) DC2H
         WRITE (LUPRI, 1006) DC2E
      END IF
C
C     Allocations
C
      KNHCC  = 1
      IF (DC2E) THEN
         KLAST  = KNHCC  + (18*KCKT34 + 1)/IRAT
      ELSE
         KLAST  = KNHCC  + ( 2*KCKT34 + 1)/IRAT
      END IF
      IF (KLAST .GT. LWORK) CALL STOPIT('U2DRIV',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      MWC2DR = MAX(MWC2DR,KLAST)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
      IF (TKTIME) TIMSTR = SECOND()
      CALL U2HINT(AOINT,IODD12,IODD34,HCINT,HCRNT,COEF34,CONT3,
     &            CONT4,WORK(KLAST),LWRK,LMNV34,NPCO3,NPCO4,NUCS34,
     &            IPRINT,WORK(KNHCC),INDHER,NPNT34,NRED34,R12EIN)
      IF (TKTIME) TC2HIN = TC2HIN + SECOND() - TIMSTR
      LWTOT  = LWTOT - KLAST
C
      RETURN
 1005 FORMAT (' DC2H    ',L7)
 1006 FORMAT (' DC2E    ',L7)
 1010 FORMAT (' NHKT3/4 ',2I7)
 1020 FORMAT (' KCKT3/4 ',2I7)
 1025 FORMAT (' KHKT3/4 ',2I7)
 1030 FORMAT (' KCKT34  ',I7)
 1040 FORMAT (' DIAG34  ',L7)
 1060 FORMAT (' NUC34   ',I7)
 1070 FORMAT (' I340:   ',3I7)
      END
C  /* Deck u2hint */
      SUBROUTINE U2HINT(AOINT,IODD12,IODD34,HCINT,HCRNT,
     &                  COEF34,CONT3,CONT4,
     &                  WORK,LWORK,LMNV34,NPCO3,NPCO4,NUCS34,IPRINT,
     &                  NHCC,INDHER,NPNT34,NRED34,R12EIN)
C     Copy of C2HINT for [r12,T1+T2] integrals (WK/UniKA/15-11-2002).
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.00D0)
#include <maxaqn.h>
      LOGICAL DORDER(2), DHOVEC, R12EIN
      DIMENSION IXDER(6), IYDER(6), IZDER(6), IDHO(2,2), LMNV34(*)
      DIMENSION AOINT(*), IODD12(KCKT12,2), IODD34(KCKT34,2),
     &          HCINT(*), WORK(LWORK),
     &          COEF34(*),NHCC(KCKT34,2), CONT3(*), CONT4(*),
     &          NPCO3(*), NPCO4(*), NUCS34(*), INDHER(*), NPNT34(*),
     &          NRED34(*), HCRNT(*)
#include <twoao.h>
#include <twosta.h>
#include <derzer.h>
      COMMON /COMC2H/ JODDIF(27), ITADD(27), IUADD(27), IVADD(27),
     &                FACTOR(27), IOFFHC(27), IOFFCC(27)
      COMMON /DHCINF/ IDHC(10)
      COMMON /DHODIR/ DHOVEC(18)
      COMMON /DHOADR/ IHOVEC(18)
      COMMON /DHOFAC/ FHOVEC(18)
      DATA IXDER /0,1,1,2,2,2/
     &     IYDER /0,1,0,2,1,0/
     &     IZDER /0,0,1,0,1,2/
C
C                *********************
C                ***** HOVEC(18) *****
C                *********************
C
      DATA IDHO /3,0,  12,6/
C
C     ARRANGEMENT OF VECTORS HOVEC(18)
C
C     1   XP00 YP00 ZP00
C     4   XQ00 YQ00 ZQ00
C     7   XXPP XYPP XZPP YYPP YZPP ZZPP
C     13  XXQQ XYQQ XZQQ YYQQ YZQQ ZZQQ
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from U2HINT','*',103)
C
C     **************************
C     ***** Initialization *****
C     **************************
C
C     NCOOR, ITADD, IUADD, IVADD, JODDIF, IOFFCC & FACTOR
C
      ICOOR = 1
      ITADD(ICOOR)  = 0
      IUADD(ICOOR)  = 0
      IVADD(ICOOR)  = 0
      IOFFHC(ICOOR) = IDHC(1)
      JODDIF(ICOOR) = 0
      IOFFCC(ICOOR) = IZERO
      FACTOR(ICOOR) = FZERO
      DORDER(1) = IDERIV .GE. 1
      DORDER(2) = IDERIV .EQ. 2
      NCOOR = ICOOR
C
      NHCMAX = 0
      DO 200 ITYPE = 1, 2
         IF (ITYPE .EQ. 1) THEN
            ICMPMX = KCKT34
         ELSE
            ICMPMX = KHKT34
         END IF
         DO 210 ICMP34 = 1, ICMPMX
            IF (IHRSYM .EQ. 0) THEN
               NHCCMP = NCOOR*KHKT12
            ELSE
               NHCCMP = 0
               JODD34 = IODD34(ICMP34,ITYPE)
               DO 220 ICMP12 = 1, KHKT12
                  IF (IODD12(ICMP12,2) .EQ. JODD34)
     &            NHCCMP = NHCCMP + 1
  220          CONTINUE
            END IF
            NHCC(ICMP34,ITYPE) = NHCCMP
            NHCMAX = MAX(NHCMAX,NHCCMP)
  210    CONTINUE
  200 CONTINUE
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 2010)
         WRITE (LUPRI, 2020) DHOVEC(1)
         WRITE (LUPRI, 2030) (DHOVEC(I), I = 1, 3)
         WRITE (LUPRI, 2031) (DHOVEC(I), I = 4, 6)
         WRITE (LUPRI, 2032) (DHOVEC(I), I = 7, 12)
         WRITE (LUPRI, 2033) (DHOVEC(I), I = 13, 18)
         WRITE (LUPRI, 2050) PATH1
         WRITE (LUPRI, 2060) KCKT12, KCKT34
         WRITE (LUPRI, 2065) KHKT12, KHKT34
         WRITE (LUPRI, 2070) NORB12, NORB34
         WRITE (LUPRI, 2080) NUC34
         WRITE (LUPRI, 2120) IHRSYM
         WRITE (LUPRI, 2140) NCOOR
         WRITE (LUPRI, 2150) (JODDIF(I), I = 1, NCOOR)
         WRITE (LUPRI, 2160) (ITADD(I), I = 1, NCOOR)
         WRITE (LUPRI, 2170) (IUADD(I), I = 1, NCOOR)
         WRITE (LUPRI, 2180) (IVADD(I), I = 1, NCOOR)
         WRITE (LUPRI, 2190) (IOFFCC(I), I = 1, NCOOR)
         WRITE (LUPRI, 2200) (FACTOR(I), I = 1, NCOOR)
         WRITE (LUPRI, '(1X,A,10I7)') 'NHCC 1',(NHCC(I,1),I=1,KCKT34)
         WRITE (LUPRI, '(1X,A,10I7)') 'NHCC 2',(NHCC(I,2),I=1,KHKT34)
      END IF
 2010 FORMAT (1X,   '           INITIALIZATION ',/)
 2020 FORMAT(1X,'DH0000',L7)
 2030 FORMAT(1X,'DHXP00',3L7)
 2031 FORMAT(1X,'DHXQ00',3L7)
 2032 FORMAT(1X,'DHXXPP',6L7)
 2033 FORMAT(1X,'DHXXQQ',6L7)
 2050 FORMAT(1X,'PATH1 ',L7)
 2060 FORMAT(1X,'KCKT  ',3I7)
 2065 FORMAT(1X,'KHKT  ',3I7)
 2070 FORMAT(1X,'NORB  ',2I7)
 2080 FORMAT(1X,'NUC34 ',I7)
 2120 FORMAT(1X,'IHRSYM',I7)
 2140 FORMAT(1X,'NCOOR ',I7)
 2150 FORMAT(1X,'JODDIF',(10I7))
 2160 FORMAT(1X,'ITADD ',(10I7))
 2170 FORMAT(1X,'IUADD ',(10I7))
 2180 FORMAT(1X,'IVADD ',(10I7))
 2190 FORMAT(1X,'IOFFCC',(10I7))
 2200 FORMAT(1X,'FACTOR',(10F5.2))
C
      IF (NHCMAX .GT. 0) THEN
C
C        Work space allocations in U2HIN1
C
C        SSINT | CCONT | CCPRIM | ETUV
C                               | SCR1   | SCR2
C                      | CSINT
C
         LSSINT = NO1234*NHCMAX*KHKT34
         LCCPRM = NCCPP*NHCMAX
         IF (SPHR3 .AND. SPHR4) THEN
            LCSINT = NO1234*NHCMAX*KHKT4
         ELSE
            LCSINT = 0
         END IF
         IF (SPHR34) THEN
            LCCONT = NO1234*NHCMAX*KCKT34
         ELSE
            LCCONT = 0
         END IF
         IF (GEN34) THEN
            LSCR1= NUCR3*NUCR4*NHCMAX*NORB12
            LSCR2= NORR3*NUCR4*NHCMAX*NORB12
         ELSE
            LSCR1 = 0
            LSCR2 = 0
         END IF
C
         KSSINT = 1
         KCCONT = KSSINT + LSSINT
C
         KCCPRM = KCCONT + LCCONT
         KETUV  = KCCPRM + LCCPRM
         KLAST1 = KETUV  + NCCPP
C
         KSCR1  = KCCPRM + LCCPRM
         KSCR2  = KSCR1  + LSCR1
         KLAST2 = KSCR2  + LSCR2
C
         KCSINT = KCCONT + LCCONT
         KLAST3 = KCSINT + LCSINT
C
         KLAST = MAX(KLAST1,KLAST2,KLAST3)
         IF (KLAST.GT.LWORK) CALL STOPIT('U2HINT',' ',KLAST,LWORK)
         MWC2HI = MAX(MWC2HI,KLAST)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
         CALL U2HIN1(AOINT,IODD12,IODD34,WORK(KCCPRM),NHCC,HCINT,
     &               HCRNT,COEF34,CONT3,CONT4,WORK(KETUV),WORK(KCCONT),
     &               WORK(KSCR1),WORK(KSCR2),NCOOR,LMNV34,NPCO3,NPCO4,
     &               NUCS34,IPRINT,INDHER,NPNT34,NRED34,
     &               WORK(KSSINT),WORK(KCSINT),NHCMAX,R12EIN)
         LWTOT  = LWTOT - KLAST
      END IF
      RETURN
      END
C  /* Deck u2hin1 */
      SUBROUTINE U2HIN1(AOINT,IODD12,IODD34,CCPRIM,NHCC,HCINT,
     &                  HCRNT,COEF34,CONT3,CONT4,ETUV,
     &                  CCONT,SCR1,SCR2,NCOOR,LMNV34,
     &                  NPCO3,NPCO4,NUCS34,IPRINT,INDHER,NPNT34,NRED34,
     &                  SSINT,CSINT,NHCMAX,R12EIN)
C
C     [r12,T1+T2] integrals are computed according to:
C     W.Klopper, R.Roehse, Theor Chim Acta (1992) 83:441, Eq.(21).
C     ------------------------------------------------------------
C     (WK/UniKA/15-11-2002).
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DP5 = 0.5D0)
#include <maxaqn.h>
      INTEGER T, U, V
      COMMON /COMC2H/ JODDIF(27), ITADD(27), IUADD(27), IVADD(27),
     &                FACTOR(27), IOFFHC(27), IOFFCC(27)
      DIMENSION AOINT(NOABCD,KHKTCD*KHKTAB,*),
     &          HCRNT(NCCPP,NTUV34,KHKT12,*),
     &          IODD12(KCKT12,2), IODD34(KCKT34,2),
     &          CCPRIM(NCCPP,*), HCINT(NCCPP,NTUV34,KHKT12,*),
     &          SSINT(NO1234*NHCMAX,KHKT34), CSINT(NO1234*NHCMAX,KHKT4),
     &          CCONT(NO1234*NHCMAX,KCKT34),
     &          ETUV(NCCPP), SCR1(*), SCR2(*),
     &          COEF34(MXUC34,0:JMAX3+JMAX4,0:JMAX3,0:JMAX4,3,*),
     &          LMNV34(KCKMX,5,2), NHCC(KCKT34,2),
     &          CONT3(*), CONT4(*), NPCO3(*), NPCO4(*), NUCS34(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), NPNT34(*), NRED34(*)
      LOGICAL R12EIN
#include <hertop.h>
#include <twoao.h>
#include <twocom.h>
#include <sphtrm.h>
#include <codata.h>
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from U2HIN1','*',103)
C
      IF (PATH1) THEN
         SGN = - D1
      ELSE
         SGN = D1
      END IF
      INCRMT = I340X + 1
      INCRMU = I340Y + 1
      INCRMV = I340Z + 1
      ICMP34 = 0
      DO 100 ICOMP3 = 1,KCKT3
         L3 = LMNV34(ICOMP3,1,1)
         M3 = LMNV34(ICOMP3,2,1)
         N3 = LMNV34(ICOMP3,3,1)
         MAX4 = KCKT4
         IF (DIAC34) MAX4 = ICOMP3
         DO 150 ICOMP4 = 1, MAX4
            ICMP34 = ICMP34 + 1
            NHCCMP = NHCC(ICMP34,1)
            IF (NHCCMP .GT. 0) THEN
C
C              *******************************
C              ***** Primitive Integrals *****
C              *******************************
C
               JODD34 = IODD34(ICMP34,1)
               CALL DZERO(CCPRIM,NHCCMP*NCCPP)
C
               IF (KHKT34 .EQ. 1) THEN
C
C                 Expansion coefficients := (a-b)/(a+b)
C                 =====================================
C
                  DO 220 I = 1, NUC34
                     ECOEFI = COEF34(I,0,0,0,1,3)
                     IF (R12EIN) ECOEFI = ECOEFI * DP5
                     IJ = I
                     DO 230 J = 1, NORB12
                        ETUV(IJ) = ECOEFI
                        IJ = IJ + NUC34
  230                CONTINUE
  220             CONTINUE
C
C                 Cartesian integrals := (a-b)/(a+b) (ab|1/r12|cd)
C                 ================================================
C
                  IF (R12EIN) THEN
                     INTHCT = INDHER(2,0,0)
                     INTHCU = INDHER(0,2,0)
                     INTHCV = INDHER(0,0,2)
                     ICCTYP = 1
                     DO 215 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. 0) THEN
                        DO 245 J = 1, NCCPP
                           HHCRNT = HCRNT(J,INTHCT,I,1)
     &                            + HCRNT(J,INTHCU,I,1)
     &                            + HCRNT(J,INTHCV,I,1)
                           CCPRIM(J,ICCTYP) = ETUV(J)*HHCRNT
  245                   CONTINUE
                        ICCTYP = ICCTYP + 1
                     END IF
  215                CONTINUE
                  ELSE
                     INTHC = INDHER(0,0,0)
                     ICCTYP = 1
                     DO 210 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. 0) THEN
                        DO 240 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = ETUV(J)*HCINT(J,INTHC,I,1)
  240                   CONTINUE
                        ICCTYP = ICCTYP + 1
                     END IF
  210                CONTINUE
                  END IF
C
                  IF (I340X .EQ. 0) THEN
C
C                    Expansion coefficients := Rx*Ey*Ez
C                    ==================================
C
                     DO 221 I = 1, NUC34
                        ECOEFI = SGN*COEF34(I,0,0,0,1,2)
                        IJ = I
                        DO 231 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  231                   CONTINUE
  221                CONTINUE
C
C                    Cartesian integrals := dRx dPx (ab|r12|cd)
C                    ==========================================
C
                     INTHC = INDHER(1,0,0)
                     ICCTYP = 1
                     DO 211 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. 0) THEN
                        DO 241 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                            + ETUV(J)*HCRNT(J,INTHC,I,1)
  241                   CONTINUE
                       ICCTYP = ICCTYP + 1
                     END IF
  211                CONTINUE
                  END IF
C
                  IF (I340Y .EQ. 0) THEN
C
C                    Expansion coefficients := Ex*Ry*Ez
C                    ==================================
C
                     DO 222 I = 1, NUC34
                        ECOEFI = SGN*COEF34(I,0,0,0,2,2)
                        IJ = I
                        DO 232 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  232                   CONTINUE
  222                CONTINUE
C
C                    Cartesian integrals := dRy dPy (ab|r12|cd)
C                    ==========================================
C
                     INTHC = INDHER(0,1,0)
                     ICCTYP = 1
                     DO 212 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. 0) THEN
                        DO 242 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                            + ETUV(J)*HCRNT(J,INTHC,I,1)
  242                   CONTINUE
                       ICCTYP = ICCTYP + 1
                     END IF
  212                CONTINUE
                  END IF
C
                  IF (I340Z .EQ. 0) THEN
C
C                    Expansion coefficients := Ex*Ey*Rz
C                    ==================================
C
                     DO 223 I = 1, NUC34
                        ECOEFI = SGN*COEF34(I,0,0,0,3,2)
                        IJ = I
                        DO 233 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  233                   CONTINUE
  223                CONTINUE
C
C                    Cartesian integrals := dRz dPz (ab|r12|cd)
C                    ==========================================
C
                     INTHC = INDHER(0,0,1)
                     ICCTYP = 1
                     DO 213 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. 0) THEN
                        DO 243 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                            + ETUV(J)*HCRNT(J,INTHC,I,1)
  243                   CONTINUE
                       ICCTYP = ICCTYP + 1
                     END IF
  213                CONTINUE
                  END IF
               ELSE
                  L4 = LMNV34(ICOMP4,1,2)
                  M4 = LMNV34(ICOMP4,2,2)
                  N4 = LMNV34(ICOMP4,3,2)
                  MAXT = L3 + L4
                  MAXU = M3 + M4
                  MAXV = N3 + N4
                  MINT = IBTAND(MAXT,INCRMT - 1)
                  MINU = IBTAND(MAXU,INCRMU - 1)
                  MINV = IBTAND(MAXV,INCRMV - 1)
                  DO 300 V = MINV, MAXV, INCRMV
                  DO 300 U = MINU, MAXU, INCRMU
                  DO 300 T = MINT, MAXT, INCRMT
C
C                    Expansion coefficients := (a-b)/(a+b)
C                    =====================================
C
                     DO 310 I = 1, NUC34
                        ECOEFI = COEF34(I,T,L3,L4,1,1)
     &                         * COEF34(I,U,M3,M4,2,1)
     &                         * COEF34(I,V,N3,N4,3,1)
     &                         * COEF34(I,0,0 ,0 ,1,3)
                        IF (R12EIN) ECOEFI = ECOEFI * DP5
                        IJ = I
                        DO 320 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  320                   CONTINUE
  310                CONTINUE
C
C                    Cartesian integrals := (a-b)/(a+b) (ab|1/r12|cd)
C                    ================================================
C
                     IF (R12EIN) THEN
                        INTHCT = INDHER(T+2,U,V)
                        INTHCU = INDHER(T,U+2,V)
                        INTHCV = INDHER(T,U,V+2)
                        ICCTYP = 1
                        DO 335 I = 1, KHKT12
                        IF (IODD12(I,2) .EQ. JODD34) THEN
                           DO 345 J = 1, NCCPP
                              HHCRNT = HCRNT(J,INTHCT,I,1)
     &                               + HCRNT(J,INTHCU,I,1)
     &                               + HCRNT(J,INTHCV,I,1)
                              CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                                         + ETUV(J) * HHCRNT
  345                      CONTINUE
                           ICCTYP = ICCTYP + 1
                        END IF
  335                   CONTINUE
                     ELSE
                        INTHC = INDHER(T,U,V)
                        ICCTYP = 1
                        DO 330 I = 1, KHKT12
                        IF (IODD12(I,2) .EQ. JODD34) THEN
                           DO 340 J = 1, NCCPP
                              CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                               + ETUV(J)*HCINT(J,INTHC,I,1)
  340                      CONTINUE
                           ICCTYP = ICCTYP + 1
                        END IF
  330                   CONTINUE
                     END IF
  300             CONTINUE
C
                  MINT = IBTAND(MAXT + 1,INCRMT - 1)
                  MINU = IBTAND(MAXU    ,INCRMU - 1)
                  MINV = IBTAND(MAXV    ,INCRMV - 1)
                  DO 301 V = MINV, MAXV, INCRMV
                  DO 301 U = MINU, MAXU, INCRMU
                  DO 301 T = MINT, MAXT, INCRMT
C
C                    Expansion coefficients := Rx*Ey*Ez
C                    ==================================
C
                     DO 311 I = 1, NUC34
                        ECOEFI = COEF34(I,T,L3,L4,1,2)
     &                         * COEF34(I,U,M3,M4,2,1)
     &                         * COEF34(I,V,N3,N4,3,1)
     &                         * SGN
                        IJ = I
                        DO 321 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  321                   CONTINUE
  311                CONTINUE
C
C                    Cartesian integrals := dRx dPx (ab|r12|cd)
C                    ==========================================
C
                     INTHC = INDHER(T+1,U,V)
                     ICCTYP = 1
                     DO 331 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. JODD34) THEN
                        DO 341 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                          + ETUV(J) * HCRNT(J,INTHC,I,1)
  341                   CONTINUE
                        ICCTYP = ICCTYP + 1
                     END IF
  331                CONTINUE
  301             CONTINUE
                  MINT = IBTAND(MAXT    ,INCRMT - 1)
                  MINU = IBTAND(MAXU + 1,INCRMU - 1)
                  MINV = IBTAND(MAXV    ,INCRMV - 1)
                  DO 302 V = MINV, MAXV, INCRMV
                  DO 302 U = MINU, MAXU, INCRMU
                  DO 302 T = MINT, MAXT, INCRMT
C
C                    Expansion coefficients := Ex*Ry*Ez
C                    ==================================
C
                     DO 312 I = 1, NUC34
                        ECOEFI = COEF34(I,T,L3,L4,1,1)
     &                         * COEF34(I,U,M3,M4,2,2)
     &                         * COEF34(I,V,N3,N4,3,1)
     &                         * SGN
                        IJ = I
                        DO 322 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  322                   CONTINUE
  312                CONTINUE
C
C                    Cartesian integrals := dRy dRz (ab|r12|cd)
C                    ==========================================
C
                     INTHC = INDHER(T,U+1,V)
                     ICCTYP = 1
                     DO 332 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. JODD34) THEN
                        DO 342 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                          + ETUV(J) * HCRNT(J,INTHC,I,1)
  342                   CONTINUE
                        ICCTYP = ICCTYP + 1
                     END IF
  332                CONTINUE
  302             CONTINUE
                  MINT = IBTAND(MAXT    ,INCRMT - 1)
                  MINU = IBTAND(MAXU    ,INCRMU - 1)
                  MINV = IBTAND(MAXV + 1,INCRMV - 1)
                  DO 303 V = MINV, MAXV, INCRMV
                  DO 303 U = MINU, MAXU, INCRMU
                  DO 303 T = MINT, MAXT, INCRMT
C
C                    Expansion coefficients := Ex*Ey*Rz
C                    ==================================
C
                     DO 313 I = 1, NUC34
                        ECOEFI = COEF34(I,T,L3,L4,1,1)
     &                         * COEF34(I,U,M3,M4,2,1)
     &                         * COEF34(I,V,N3,N4,3,2)
     &                         * SGN
                        IJ = I
                        DO 323 J = 1, NORB12
                           ETUV(IJ) = ECOEFI
                           IJ = IJ + NUC34
  323                   CONTINUE
  313                CONTINUE
C
C                    Cartesian integrals := dRz dPz (ab|r12|cd)
C                    ==========================================
C
                     INTHC = INDHER(T,U,V+1)
                     ICCTYP = 1
                     DO 333 I = 1, KHKT12
                     IF (IODD12(I,2) .EQ. JODD34) THEN
                        DO 343 J = 1, NCCPP
                           CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                          + ETUV(J) * HCRNT(J,INTHC,I,1)
  343                   CONTINUE
                        ICCTYP = ICCTYP + 1
                     END IF
  333                CONTINUE
  303             CONTINUE
               END IF
C
C              ********************************
C              ***** Contracted Integrals *****
C              ********************************
C
               IF (SPHR34) THEN
                  CALL C2CONT(CCPRIM,CCONT(1,ICMP34),CONT3,CONT4,SCR1,
     &                        SCR2,NPCO3,NPCO4,NUCS34,NHCCMP,NPNT34,
     &                        NRED34,IPRINT,.TRUE.)
               ELSE
                  CALL C2CONT(CCPRIM,SSINT(1,ICMP34),CONT3,CONT4,SCR1,
     &                        SCR2,NPCO3,NPCO4,NUCS34,NHCCMP,NPNT34,
     &                        NRED34,IPRINT,.TRUE.)
               END IF
            END IF
  150    CONTINUE
  100 CONTINUE
C
C     Spherical integrals
C     ===================
C
      IF (SPHR34) THEN
         CALL C2SPHR(CCONT,CSINT,SSINT,NHCC,CSP(ISPADR(NHKT3)),
     &               CSP(ISPADR(NHKT4)),NHCMAX,IPRINT)
      END IF
C
C     ***********************************************
C     ***** Multiply by factors and distribute ******
C     ***********************************************
C
      DO 700 ICMP34 = 1, KHKT34
         JODD34 = IODD34(ICMP34,2)
         IWORK = 0
         IF (PATH1) THEN
            DO 800 ICOOR = 1, NCOOR
               ITYPE  = IOFFCC(ICOOR)
               FAC    = FACTOR(ICOOR)
               IODDIF = IBTXOR(JODDIF(ICOOR),JODD34)
               ICMP   = ICMP34
               DO 810 ICMP12 = 1, KHKT12
                  IF (IODD12(ICMP12,2) .EQ. IODDIF) THEN
                     DO 820 I = 1, NO1234
                        AOINT(I,ICMP,ITYPE) = AOINT(I,ICMP,ITYPE)
     &                              + FAC*SSINT(IWORK + I,ICMP34)
  820                CONTINUE
                     IWORK = IWORK + NO1234
                  END IF
                  ICMP = ICMP + KHKT34
  810          CONTINUE
  800       CONTINUE
         ELSE
            ICMP0 = (ICMP34 - 1)*KHKT12
            DO 830 ICOOR = 1, NCOOR
               ITYPE  = IOFFCC(ICOOR)
               FAC    = FACTOR(ICOOR)
               IODDIF = IBTXOR(JODDIF(ICOOR),JODD34)
               DO 840 ICMP12 = 1, KHKT12
               IF (IODD12(ICMP12,2) .EQ. IODDIF) THEN
                  ICMP = ICMP0 + ICMP12
                  DO 850 I = 1, NORB12
                     IJ = I
                     DO 860 J = 1, NORB34
                        AOINT(IJ,ICMP,ITYPE) = AOINT(IJ,ICMP,ITYPE)
     &                                  + FAC*SSINT(IWORK+J,ICMP34)
                        IJ = IJ + NORB12
  860                CONTINUE
                     IWORK = IWORK  + NORB34
  850             CONTINUE
               END IF
  840          CONTINUE
  830       CONTINUE
         END IF
  700 CONTINUE
C
C     *************************
C     ***** Print Section *****
C     *************************
C
      IF (IPRINT .GE. 15) THEN
         CALL HEADER('Final spherical integrals - U2HIN1',-1)
         DO 900 ICOOR = 1, NCOOR
            ICMP = 0
            K = IOFFCC(ICOOR)
            DO 910 ICMPAB = 1, KHKTAB
               DO 920 ICMPCD = 1, KHKTCD
                  WRITE (LUPRI, '(/1X,A,3I3)')' ICOOR, ICMPAB, ICMPCD ',
     &                                          ICOOR, ICMPAB, ICMPCD
                  ICMP = ICMP + 1
                  WRITE(LUPRI,'(1P,6D12.4)')(AOINT(I,ICMP,K),I=1,NO1234)
  920          CONTINUE
  910       CONTINUE
  900    CONTINUE
      END IF
      RETURN
      END
C  /* Deck lab64u */
      SUBROUTINE LAB64U(SO,N,IA,IB,IC,ID,XIJKL,THRESH,GTTHRS,IPRINT)
C
C     Make label and change sign of [r12,T1], [r12,T2] integrals
C     according to permutation of ia,ib,ic,id.
C
C     WK-02/25/1995-IPS/ETHZ
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.00D0)
#include <r12int.h>
      LOGICAL GTTHRS
      INTEGER*4 M(2)
      DIMENSION SO(N,*)
      EQUIVALENCE(X,M(1))
      IF (NOPP12 .GT. 1) THEN
         IF (U12INT) THEN
            IF (IA .GT. IB) THEN
               M(1) = IA*65536 + IB
            ELSE
               M(1) = IB*65536 + IA
               SO(1,IADU12) = - SO(1,IADU12)
            END IF
            IF (IC .GT. ID) THEN
               M(2) = IC*65536 + ID
            ELSE
               M(2) = ID*65536 + IC
               SO(1,IADU21) = - SO(1,IADU21)
            END IF
            IF (M(2) .GT. M(1)) THEN
               MMMM = M(1)
               M(1) = M(2)
               M(2) = MMMM
               XXXX = SO(1,IADU12)
               SO(1,IADU12) = SO(1,IADU21)
               SO(1,IADU21) = XXXX
            END IF
         ELSE
            M(1) = MAX(IA,IB)*65535 + IA + IB
            M(2) = MAX(IC,ID)*65535 + IC + ID
            IF (M(2) .GT. M(1)) THEN
               MMMM = M(1)
               M(1) = M(2)
               M(2) = MMMM
            END IF
         END IF
         GTTHRS = .FALSE.
         DO 100 I = 1, NOPP12
            IF (ABS(SO(1,I)) .LE. THRESH) THEN
               SO(1,I) = D0
            ELSE
               GTTHRS = .TRUE.
            END IF
  100    CONTINUE
         IF (GTTHRS .AND. IPRINT .GE. 10)
     &      WRITE(LUPRI,'(4I4,4(1P,E15.6))')
     &      IA, IB, IC, ID, (SO(1,I),I=1,NOPP12)
      ELSE
         M(1) = MAX(IA,IB)*65535 + IA + IB
         M(2) = MAX(IC,ID)*65535 + IC + ID
         IF (M(2) .GT. M(1)) THEN
            MMMM = M(1)
            M(1) = M(2)
            M(2) = MMMM
         END IF
         GTTHRS = ABS(SO(1,1)) .GT. THRESH
CWMK     IF (GTTHRS .AND. IPRINT .GE. 10)
         IF (GTTHRS)
     &      WRITE(LUPRI,'(4I4,1P,E15.6,2I10)')
     &      IA, IB, IC, ID, SO(1,1), M(1), M(2)
      END IF
      XIJKL = X
      RETURN
      END
