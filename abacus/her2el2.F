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
C  /* Deck c2driv */
      SUBROUTINE C2DRIV(AOINT,HCINT,COEF34,INDHER,CONT3,CONT4,WORK,
     &                  LWORK,NPCO3,NPCO4,NUCS34,IPRINT,LMNV12,LMNV34,
     &                  IODD12,IODD34,NPNT34,NRED34)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
#include <iratdef.h>
#include <maxaqn.h>
      DIMENSION AOINT(*), HCINT(*), WORK(LWORK), COEF34(*), CONT3(*),
     &          CONT4(*), NPCO3(*), NPCO4(*), NUCS34(*), INDHER(*),
     &          LMNV12(*), LMNV34(*), IODD12(*), IODD34(*),
     &          NPNT34(*), NRED34(*)
#include <twoao.h>
#include <twosta.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from C2DRIV','*',103)
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
      IF (KLAST .GT. LWORK) CALL STOPIT('C2DRIV',' ',KLAST,LWORK)
      LWRK   = LWORK - KLAST + 1
      MWC2DR = MAX(MWC2DR,KLAST)
      LWTOT  = LWTOT + KLAST
      MWTOT  = MAX(MWTOT,LWTOT)
C
C     Derivatives on common nucleus
C     and undifferentiated integrals
C
      IF (DC2H) THEN
         IF (TKTIME) TIMSTR = SECOND()
         CALL C2HINT(AOINT,IODD12,IODD34,HCINT,COEF34,CONT3,
     &               CONT4,WORK(KLAST),LWRK,LMNV34,NPCO3,NPCO4,NUCS34,
     &               IPRINT,WORK(KNHCC),INDHER,NPNT34,NRED34)
         IF (TKTIME) TC2HIN = TC2HIN + SECOND() - TIMSTR
      END IF
C
C     Derivatives on each separate nucleus
C
      IF (DC2E) THEN
         IF (TKTIME) TIMSTR = SECOND()
         CALL C2EINT(AOINT,HCINT,COEF34,CONT3,CONT4,WORK(KLAST),LWRK,
     &               LMNV34,NPCO3,NPCO4,NUCS34,IODD12,IODD34,IPRINT,
     &               WORK(KNHCC),INDHER,NPNT34,NRED34)
         IF (TKTIME) TC2EIN = TC2EIN + SECOND() - TIMSTR
      END IF
      LWTOT  = LWTOT - KLAST
C
      RETURN
 1005 FORMAT (1X,'DC2H    ',L7)
 1006 FORMAT (1X,'DC2E    ',L7)
 1010 FORMAT (1X,'NHKT3/4 ',2I7)
 1020 FORMAT (1X,'KCKT3/4 ',2I7)
 1025 FORMAT (1X,'KHKT3/4 ',2I7)
 1030 FORMAT (1X,'KCKT34  ',I7)
 1040 FORMAT (1X,'DIAG34  ',L7)
 1060 FORMAT (1X,'NUC34   ',I7)
 1070 FORMAT (1X,'I340:   ',3I7)
      END
C  /* Deck c2hint */
      SUBROUTINE C2HINT(AOINT,IODD12,IODD34,HCINT,COEF34,CONT3,CONT4,
     &                  WORK,LWORK,LMNV34,NPCO3,NPCO4,NUCS34,IPRINT,
     &                  NHCC,INDHER,NPNT34,NRED34)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D1 = 1.00 D00)
#include <maxaqn.h>
#include <maxorb.h>
      LOGICAL DORDER(2), DHOVEC
      DIMENSION IXDER(6), IYDER(6), IZDER(6), IDHO(2,2), LMNV34(*)
      DIMENSION AOINT(*), IODD12(KCKT12,2), IODD34(KCKT34,2),
     &          HCINT(*), WORK(LWORK),
     &          COEF34(*),NHCC(KCKT34,2), CONT3(*), CONT4(*),
     &          NPCO3(*), NPCO4(*), NUCS34(*), INDHER(*), NPNT34(*),
     &          NRED34(*)
#include <twoao.h>
#include <twosta.h>
#include <derzer.h>
#include <drw2el.h>
#include <expopt.h>
      COMMON /COMC2H/ JODDIF(27), ITADD(27), IUADD(27), IVADD(27),
     &                FACTOR(27), IOFFHC(27), IOFFCC(27)
      COMMON /DHCINF/ IDHC(10), IHCSYM(10)
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
      IF (IPRINT .GT. 5) CALL TITLER('Output from C2HINT','*',103)
C
C     **************************
C     ***** Initialization *****
C     **************************
C
C     NCOOR, ITADD, IUADD, IVADD, JODDIF, IOFFCC & FACTOR
C
      ICOOR = 0
      IF (DZER) THEN
         ICOOR = ICOOR + 1
         ITADD(ICOOR)  = 0
         IUADD(ICOOR)  = 0
         IVADD(ICOOR)  = 0
         IOFFHC(ICOOR) = IDHC(1)
         JODDIF(ICOOR) = 0
         IOFFCC(ICOOR) = IZERO
         FACTOR(ICOOR) = FZERO
      END IF
      IF (EXPGRA) THEN
         DO I = 1, 2
            ICOOR = ICOOR + 1
            IOFFHC(ICOOR) = IDHC(1)
            JODDIF(ICOOR) = 0
            IF (PATH1) THEN
               IOFFCC(ICOOR) = I + 2
            ELSE
               IOFFCC(ICOOR) = I 
            END IF
            FACTOR(ICOOR) = D1 
         END DO
      ELSE IF (BPH2OO) THEN
         ICOOR = ICOOR + 1
         IOFFHC(ICOOR) = IDHC(1)
         JODDIF(ICOOR) = 0
         IOFFCC(ICOOR) = 1 
         FACTOR(ICOOR) = D1 
      ELSE
         DORDER(1) = IDERIV .GE. 1
         DORDER(2) = IDERIV .EQ. 2
         IF (IDERIV .GT. 0) THEN
            DO 100 IORDER = 1, IDERIV
               IF (DORDER(IORDER)) THEN
                  IF (PATH1 .AND. (IORDER .EQ. 1)) THEN
                     SGN = - D1
                  ELSE
                     SGN = D1
                  END IF
                  JDHO = IDHO(IPATH,IORDER)
                  DO 110 ICOOR2 = 1, (IORDER + 1)*(IORDER + 2)/2
                     IX = IORDER - IXDER(ICOOR2)
                     IY = IYDER(ICOOR2)
                     IZ = IZDER(ICOOR2)
                     JODD2 = IBTAND(
     &                       IHRSYM,MOD(IX,2)+2*MOD(IY,2)+4*MOD(IZ,2))
                     IADRHO = JDHO + ICOOR2
C
C                    Differentiation on a single electron
C
                     IADRHC = 1
                     IF (DHOVEC(IADRHO)) THEN
                        ICOOR = ICOOR + 1
                        ITADD(ICOOR) = IX
                        IUADD(ICOOR) = IY
                        IVADD(ICOOR) = IZ
                        IOFFHC(ICOOR) = IDHC(IADRHC)
                        JODDIF(ICOOR) = JODD2
                        IOFFCC(ICOOR) = IHOVEC(IADRHO)
                        FACTOR(ICOOR) = SGN*FHOVEC(IADRHO)
                     END IF
  110             CONTINUE
               END IF
  100       CONTINUE
         END IF
      END IF
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
                  IO1234 = IBTXOR(IODD12(ICMP12,2),JODD34)
                  DO 230 ICOOR = 1, NCOOR
                     IF (IO1234.EQ.JODDIF(ICOOR)) NHCCMP = NHCCMP + 1
  230             CONTINUE
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
C        Work space allocations in C2HIN1
C
C        SSINT | CCONT | CCPRIM | ETUV
C                               | SCR1   | SCR2
C                      | CSINT
C
         LSSINT = NO1234*NHCMAX*KHKT34
         LCCPRM = NCCPP*NHCMAX
         IF ((KHKT34 .EQ. 1) .AND. 
     &       (.NOT.(FINDPT.OR.DPTINT.OR.BPH2OO))) THEN
            LETUV = 0
         ELSE
            LETUV  = NCCPP
         END IF
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
         KLAST1 = KETUV  + LETUV
C
         KSCR1  = KCCPRM + LCCPRM
         KSCR2  = KSCR1  + LSCR1
         KLAST2 = KSCR2  + LSCR2
C
         KCSINT = KCCONT + LCCONT
         KLAST3 = KCSINT + LCSINT
C
         KLAST = MAX(KLAST1,KLAST2,KLAST3)
         IF (KLAST.GT.LWORK) CALL STOPIT('C2HINT',' ',KLAST,LWORK)
         MWC2HI = MAX(MWC2HI,KLAST)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
         CALL C2HIN1(AOINT,IODD12,IODD34,WORK(KCCPRM),NHCC,HCINT,
     &               COEF34,CONT3,CONT4,WORK(KETUV),WORK(KCCONT),
     &               WORK(KSCR1),WORK(KSCR2),NCOOR,LMNV34,NPCO3,NPCO4,
     &               NUCS34,IPRINT,INDHER,NPNT34,NRED34,
     &               WORK(KSSINT),WORK(KCSINT),NHCMAX)
         LWTOT  = LWTOT - KLAST
      END IF
      RETURN
      END
C  /* Deck c2hin1 */
      SUBROUTINE C2HIN1(AOINT,IODD12,IODD34,CCPRIM,NHCC,HCINT,COEF34,
     &                  CONT3,CONT4,ETUV,CCONT,SCR1,SCR2,NCOOR,LMNV34,
     &                  NPCO3,NPCO4,NUCS34,IPRINT,INDHER,NPNT34,NRED34,
     &                  SSINT,CSINT,NHCMAX)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0, DP25 = 0.25D0)
#include <maxaqn.h>
#include <maxorb.h>
      INTEGER T, U, V
      COMMON /COMC2H/ JODDIF(27), ITADD(27), IUADD(27), IVADD(27),
     &                FACTOR(27), IOFFHC(27), IOFFCC(27)
      DIMENSION AOINT(NOABCD,KHKTCD*KHKTAB,*),
     &          IODD12(KCKT12,2), IODD34(KCKT34,2),
     &          CCPRIM(NCCPP,*), HCINT(NCCPP,NTUV34,KHKT12,*),
     &          SSINT(NO1234*NHCMAX,KHKT34), CSINT(NO1234*NHCMAX,KHKT4),
     &          CCONT(NO1234*NHCMAX,KCKT34),
     &          ETUV(NCCPP), SCR1(*), SCR2(*),
     &          COEF34(MXUC34,0:JMAX3+JMAX4,0:JMAX3,0:JMAX4,3,*),
     &          LMNV34(KCKMX,5,2), NHCC(KCKT34,2),
     &          CONT3(*), CONT4(*), NPCO3(*), NPCO4(*), NUCS34(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), NPNT34(*), NRED34(*)
#include <hertop.h>
#include <twoao.h>
#include <twocom.h>
#include <sphtrm.h>
#include <codata.h>
#include <drw2el.h>
#include <expopt.h>
      DIMENSION INCDPT(3,3), INCHOO(3,3), ITYDPT(3,3), ITYHOO(3,3), 
     &          INCEXP(3,0:3)
      DATA INCDPT/2, 0, 0, 0, 2, 0, 0, 0, 2/
      DATA INCHOO/1, 0, 0, 0, 1, 0, 0, 0, 1/
      DATA ITYDPT/5, 1, 1, 1, 5, 1, 1, 1, 5/
      DATA ITYHOO/3, 1, 1, 1, 3, 1, 1, 1, 3/
      DATA INCEXP/0, 0, 0, 2, 0, 0, 0, 2, 0, 0, 0, 2/
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from C2HIN1','*',103)
      IF (FINDPT) ECOFAC = DPTFAC * DP25 * ALPHAC ** 2
      IF (DPTINT) ECOFAC = DP25 * ALPHAC ** 2
      IF (NO2DPT) ECOFAC = D0
C
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
               IF ((KHKT34 .EQ. 1).AND.
     &             (.NOT.(FINDPT.OR.DPTINT.OR.EXPGRA.OR.BPH2OO))) THEN
                  ICCTYP = 1
                  DO 200 ICOOR = 1, NCOOR
                     IODDIF= JODDIF(ICOOR)
                     INTHC = INDHER(ITADD(ICOOR),IUADD(ICOOR),
     &                                           IVADD(ICOOR))
                     DO 210 ICMP12 = 1, KHKT12
                     IF (IODD12(ICMP12,2) .EQ. IODDIF) THEN
                        DO 220 J = 1, NCCPP
                          CCPRIM(J,ICCTYP) =
     &                        HCINT(J,INTHC,ICMP12,IOFFHC(ICOOR))
  220                   CONTINUE
                        ICCTYP = ICCTYP + 1
                     END IF
  210                CONTINUE
  200             CONTINUE
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
                  IF (IDERIV .EQ. 0 .OR. FINDPT .OR. DPTINT) THEN
                     IF (DPTINT) THEN
                        ITYPE = 2
                     ELSE
                        ITYPE  = IOFFHC(1)
                     END IF
                     DO 300 V = MINV, MAXV, INCRMV
                     DO 300 U = MINU, MAXU, INCRMU
                     DO 300 T = MINT, MAXT, INCRMT
C
C                       Expansion coefficients
C                       ======================
C
                        DO 310 I = 1, NUC34
                           ECOEFI = COEF34(I,T,L3,L4,1,1)
     &                            * COEF34(I,U,M3,M4,2,1)
     &                            * COEF34(I,V,N3,N4,3,1)
                           IJ = I
                           DO 320 J = 1, NORB12
                              ETUV(IJ) = ECOEFI
                              IJ = IJ + NUC34
  320                      CONTINUE
  310                   CONTINUE
C
C                       Cartesian integrals
C                       ===================
C
                        ITUV = INDHER(T,U,V)
                        ICCTYP = 1
                        DO 330 I = 1, KHKT12
                        IF (IODD12(I,2) .EQ. JODD34) THEN
                           DO 340 J = 1, NCCPP
                              CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                           + ETUV(J)*HCINT(J,ITUV,I,ITYPE)
  340                      CONTINUE
                           ICCTYP = ICCTYP + 1
                        END IF
  330                   CONTINUE
  300                CONTINUE
                     IF (FINDPT.OR.DPTINT) THEN
C
C     =============== Add or do first-order DPT perturbation ==============
C     Wim Klopper & Asger Halkier, University of Utrecht, November 25, 1999
C
                     DO 333 ICOOR = 1, 3
                        DO 301 V = MINV, MAXV + INCDPT(3,ICOOR), INCRMV
                        DO 301 U = MINU, MAXU + INCDPT(2,ICOOR), INCRMU
                        DO 301 T = MINT, MAXT + INCDPT(1,ICOOR), INCRMT
C
C                          Expansion coefficients
C                          ======================
C
                           DO 311 I = 1, NUC34
                            ECOEFI = COEF34(I,T,L3,L4,1,ITYDPT(1,ICOOR))
     &                             * COEF34(I,U,M3,M4,2,ITYDPT(2,ICOOR))
     &                             * COEF34(I,V,N3,N4,3,ITYDPT(3,ICOOR))
     &                             * ECOFAC
                            IJ = I
                            DO 321 J = 1, NORB12
                              ETUV(IJ) = ECOEFI
                              IJ = IJ + NUC34
  321                       CONTINUE
  311                      CONTINUE
C
C                          Cartesian integrals
C                          ===================

                           ITUV = INDHER(T,U,V)
                           ICCTYP = 1
                           DO 331 I = 1, KHKT12
                            IF (IODD12(I,2) .EQ. JODD34) THEN
                              DO 341 J = 1, NCCPP
                                 CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                              + ETUV(J)*HCINT(J,ITUV,I,1)
  341                         CONTINUE
                              ICCTYP = ICCTYP + 1
                            END IF
  331                      CONTINUE
  301                   CONTINUE
  333                CONTINUE
                     END IF
C
C                 Orbital-orbital interaction 
C                 ===========================
C
C                 tuh & wk 07.04.03
C
                  ELSE IF (BPH2OO) THEN
                     DO ICOOR = 1, 3
                        MAXT = L3 + L4 + INCHOO(1,ICOOR)
                        MAXU = M3 + M4 + INCHOO(2,ICOOR)
                        MAXV = N3 + N4 + INCHOO(3,ICOOR)
                        MINT = IBTAND(MAXT,INCRMT - 1)
                        MINU = IBTAND(MAXU,INCRMU - 1)
                        MINV = IBTAND(MAXV,INCRMV - 1)
                        DO V = MINV, MAXV, INCRMV
                        DO U = MINU, MAXU, INCRMU
                        DO T = MINT, MAXT, INCRMT
                           DO I = 1, NUC34
                              ECOEFI=COEF34(I,T,L3,L4,1,ITYHOO(1,ICOOR))
     &                              *COEF34(I,U,M3,M4,2,ITYHOO(2,ICOOR))
     &                              *COEF34(I,V,N3,N4,3,ITYHOO(3,ICOOR))
                              IJ = I
                              DO J = 1, NORB12
                                 ETUV(IJ) = ECOEFI
                                 IJ = IJ + NUC34
                              END DO
                           END DO
                           ITUV = INDHER(T,U,V)
                           ICCTYP = 1
                           DO I = 1, KHKT12
                              IF (IODD12(I,2) .EQ. JODD34) THEN
                                 DO J = 1, NCCPP
                                    CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                                 + ETUV(J)*HCINT(J,ITUV,I,ICOOR)
                                 END DO
                                 ICCTYP = ICCTYP + 1
                              END IF
                           END DO
                        END DO
                        END DO
                        END DO
                     END DO
C
C                 Exponent gradient
C                 =================
C
                  ELSE IF (EXPGRA) THEN
                     DO K = 0, 3
                     DO V = MINV, MAXV + INCEXP(3,K), INCRMV
                     DO U = MINU, MAXU + INCEXP(2,K), INCRMU
                     DO T = MINT, MAXT + INCEXP(1,K), INCRMT
                        ICCTYP = 1
                        DO ICOOR = 1, 2
C
C                          Expansion coefficients
C                          ======================
C
                           DO I = 1, NUC34
                              IF (K.EQ.0) THEN
                                 ECX = COEF34(I,T,L3,L4,1,1+ICOOR)
     &                               * COEF34(I,U,M3,M4,2,1)
     &                               * COEF34(I,V,N3,N4,3,1)
                                 ECY = COEF34(I,T,L3,L4,1,1)
     &                               * COEF34(I,U,M3,M4,2,1+ICOOR)
     &                               * COEF34(I,V,N3,N4,3,1)
                                 ECZ = COEF34(I,T,L3,L4,1,1)
     &                               * COEF34(I,U,M3,M4,2,1)
     &                               * COEF34(I,V,N3,N4,3,1+ICOOR)
                                 ECF = ECX + ECY + ECZ
                              ELSE IF (ICOOR.EQ.1) THEN
                                 ECF =-COEF34(I,T,L3+INCEXP(1,K),L4,1,1)
     &                                *COEF34(I,U,M3+INCEXP(2,K),M4,2,1)
     &                                *COEF34(I,V,N3+INCEXP(3,K),N4,3,1)
                              ELSE
                                 ECF =-COEF34(I,T,L3,L4+INCEXP(1,K),1,1)
     &                                *COEF34(I,U,M3,M4+INCEXP(2,K),2,1)
     &                                *COEF34(I,V,N3,N4+INCEXP(3,K),3,1)
                              END IF
                              IJ = I
                              DO J = 1, NORB12
                                 ETUV(IJ) = ECF
                                 IJ = IJ + NUC34
                              END DO
                           END DO
C
C                          Cartesian integrals
C                          ===================
C
                           ITUV = INDHER(T,U,V)
                           DO I = 1, KHKT12
                           IF (IODD12(I,2) .EQ. JODD34) THEN
                              DO J = 1, NCCPP
                                 CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                              + ETUV(J)*HCINT(J,ITUV,I,1)
                              END DO
                              ICCTYP = ICCTYP + 1
                           END IF
                           END DO
                        END DO
                     END DO
                     END DO
                     END DO
                     END DO
                  ELSE
                     DO 400 V = MINV, MAXV, INCRMV
                     DO 400 U = MINU, MAXU, INCRMU
                     DO 400 T = MINT, MAXT, INCRMT
C
C                       Expansion coefficients
C                       ======================
C
                        DO 410 I = 1, NUC34
                           ECOEFI = COEF34(I,T,L3,L4,1,1)
     &                            * COEF34(I,U,M3,M4,2,1)
     &                            * COEF34(I,V,N3,N4,3,1)
                           IJ = I
                           DO 420 J = 1, NORB12
                              ETUV(IJ) = ECOEFI
                              IJ = IJ + NUC34
  420                      CONTINUE
  410                   CONTINUE
C
C                       Cartesian integrals
C                       ===================
C
                        ICCTYP = 1
                        DO 430 ICOOR = 1, NCOOR
                           JT = T + ITADD(ICOOR)
                           JU = U + IUADD(ICOOR)
                           JV = V + IVADD(ICOOR)
                           INTHC  = INDHER(JT,JU,JV)
                           ITYPE  = IOFFHC(ICOOR)
                           IODDIF = IBTXOR(JODDIF(ICOOR),JODD34)
                           DO 440 I = 1, KHKT12
                           IF (IODD12(I,2) .EQ. IODDIF) THEN
                              DO 450 J = 1, NCCPP
                                 CCPRIM(J,ICCTYP) = CCPRIM(J,ICCTYP)
     &                              + ETUV(J)*HCINT(J,INTHC,I,ITYPE)
  450                         CONTINUE
                              ICCTYP = ICCTYP + 1
                           END IF
  440                      CONTINUE
  430                   CONTINUE
  400                CONTINUE
                  END IF
               END IF
C
C              ********************************
C              ***** Contracted Integrals *****
C              ********************************
C
               IF (SPHR34) THEN
                  CALL C2CONT(CCPRIM,CCONT(1,ICMP34),CONT3,CONT4,SCR1,
     &                        SCR2,NPCO3,NPCO4,NUCS34,NHCCMP,NPNT34,
     &                        NRED34,IPRINT,.FALSE.)
               ELSE
                  CALL C2CONT(CCPRIM,SSINT(1,ICMP34),CONT3,CONT4,SCR1,
     &                        SCR2,NPCO3,NPCO4,NUCS34,NHCCMP,NPNT34,
     &                        NRED34,IPRINT,.FALSE.)
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
         CALL HEADER('Final spherical integrals - C2HIN1',-1)
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
C  /* Deck c2eint */
      SUBROUTINE C2EINT(AOINT,HCINT,COEF34,CONT3,CONT4,WORK,LWORK,
     &                  LMNV34,NPCO3,NPCO4,NUCS34,IODD12,IODD34,IPRINT,
     &                  NHCC,INDHER,NPNT34,NRED34)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      LOGICAL DCOOR, DCOOR0, DCOORH, DCOORE, DATOM,
     &        DATOM0, DATOMH, DATOME, DCORDH(2), DCORDE(2), DEOVEC,
     &        DETVEC, DMXVEC, DEO, DET, DMX
      DIMENSION IEOFF(0:2), IXDER(6), IYDER(6), IZDER(6), JODD1(3),
     &          IEND(6), IDEO(6,2,2), IDET(3,3,2,2,2), IDMX(3,3,2,2),
     &          NHCC(KCKT34,2,9), LMNV34(*)
      DIMENSION AOINT(*), HCINT(*), COEF34(*),
     &          IODD12(KCKT12,2), IODD34(KCKT34,2),
     &          CONT3(*), CONT4(*), WORK(LWORK), NPCO3(*), NPCO4(*),
     &          NUCS34(*), INDHER(*), NPNT34(*), NRED34(*)
#include <twoao.h>
#include <twosta.h>
      COMMON /DEODIR/ DEOVEC(54)
      COMMON /DEOFAC/ FEOVEC(54)
      COMMON /DEOADR/ IEOVEC(54)
      COMMON /DETDIR/ DETVEC(36)
      COMMON /DETFAC/ FETVEC(36)
      COMMON /DETADR/ IETVEC(36)
      COMMON /DMXDIR/ DMXVEC(36)
      COMMON /DMXFAC/ FMXVEC(36)
      COMMON /DMXADR/ JMXVEC(36)
      COMMON /DHCINF/ IHCVEC(10), IHCSYM(10)
      COMMON /COMC2E/ JSTRHC(27,9), JENDHC(27,9),
     &                INCMXT(9), INCMXU(9), INCMXV(9),
     &                IOFFEX(27), IOFFEY(27), IOFFEZ(27),
     &                IOFF2X(27), IOFF2Y(27), IOFF2Z(27),
     &                IOFFHC(81), IOFFCC(81), JODDIF(81),
     &                FACTOR(81), JSTRAT(9), JENDAT(9)
      DATA IXDER /0,1,1,2,2,2/
     &     IYDER /0,1,0,2,1,0/
     &     IZDER /0,0,1,0,1,2/
      DATA IEND  /3,4,4,3,4,3/
      DATA IEOFF /1,2,4/
C
C            ************************
C            ***** DEO(6,4,2,2) *****
C            ************************
C
      DEO(ICOOR2,IATOM2,IDER,IPATH)
     &     = DEOVEC(IDEO(ICOOR2,IDER,IPATH) + IATOM2)
C
      DATA IDEO / 2,  6, 10,  0,  0,  0,
     &           15, 22, 30, 37, 44, 51,
     &            0,  4,  8,  0,  0,  0,
     &           12, 18, 26, 34, 40, 48/
C
C     Arrangement of vector DEOVEC(54)
C
C     1   XA00 XB00             3   XC00 XD00
C     5   YA00 YB00             7   YC00 YD00
C     9   ZA00 ZB00             11  ZC00 ZD00
C     13  XXAA XXAB XXBB        16  XXCC XXCD XXDD
C     19  XYAA XYAB XYBA XYBB   23  XYCC XYCD XYDC XYDD
C     27  XZAA XZAB XZBA XZBB   31  XZCC XZCD XZDC XZDD
C     35  YYAA YYAB YYBB        38  YYCC YYCD YYDD
C     41  YZAA YZAB YZBA YZBB   45  YZCC YZCD YZDC YZDD
C     49  ZZAA ZZAB ZZBB        52  ZZCC ZZCD ZZDD
C
C            **************************
C            ***** DET(3,3,2,2,2) *****
C            **************************
C
      DET(ICOOR1,ICOOR2,IATOM1,IATOM2,IPATH)
     &     = DETVEC(IDET(ICOOR2,ICOOR1,IATOM2,IATOM1,IPATH))
C
      DATA IDET /  1,  2,  3,  4,  5,  6,  7,  8,  9,
     &            10, 11, 12, 13, 14, 15, 16, 17, 18,
     &            19, 20, 21, 22, 23, 24, 25, 26, 27,
     &            28, 29, 30, 31, 32, 33, 34, 35, 36,
     &             1,  4,  7,  2,  5,  8,  3,  6,  9,
     &            19, 22, 25, 20, 23, 26, 21, 24, 27,
     &            10, 13, 16, 11, 14, 17, 12, 15, 18,
     &            28, 31, 34, 29, 32, 35, 30, 33, 36/
C
C     Arrangement of vector DETVEC(36)
C
C     1   XXAC XYAC XZAC
C     4   XYCA YYAC YZAC
C     7   XZCA YZCA ZZAC
C
C     10  XXAD XYAD XZAD
C     13  XYDA YYAD YZAD
C     16  XZDA YZDA ZZAD
C
C     19  XXBC XYBC XZBC
C     22  XYCB YYBC YZBC
C     25  XZCB YZCB ZZBC
C
C     28  XXBD XYBD XZBD
C     31  XYDB YYBD YZBD
C     34  XZDB YZDB ZZBD
C
C            ************************
C            ***** DMX(3,3,2,2) *****
C            ************************
C
      DMX(ICOOR1,ICOOR2,IATOM2,IPATH)
     *     = DMXVEC(IDMX(ICOOR2,ICOOR1,IATOM2,IPATH))
C
      DATA IDMX / 1,  2,  3,  4,  5,  6,  7,  8,  9,
     *           10, 11, 12, 13, 14, 15, 16, 17, 18,
     *           19, 22, 25, 20, 23, 26, 21, 24, 27,
     *           28, 31, 34, 29, 32, 35, 30, 33, 36/
     *
C
C     Arrangement of vector DMXVEC(36)
C
C     1   XXPC XYPC XZPC
C     4   XYCP YYPC YZPC
C     7   XZCP YZCP ZZPC
C
C     10  XXPD XYPD XZPD
C     13  XYDP YYPD YZPD
C     16  XZDP YZDP ZZPD
C
C     19  XXAQ XYAQ XZAQ
C     22  XYQA YYAQ YZAQ
C     25  XZQA YZQA ZZAQ
C
C     28  XXBQ XYBQ XZBQ
C     31  XYQB YYBQ YZBQ
C     34  XZQB YZQB ZZBQ
C
C            *********************
C            ***** HCVEC(10) *****
C            *********************
C
C     Arrangement of elements in COMMON /DHCINF/
C
C     HC00
C     HCHX HCEX1 HCEX2
C     HCHY HCEY1 HCEY2
C     HCHZ HCEZ1 HCEZ2
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from C2EINT','*',103)
C
C     ***** Initialization *****
C
      IF (CROSS) THEN
         JODD1(1) = IBTAND(IHRSYM,1)
         JODD1(2) = IBTAND(IHRSYM,2)
         JODD1(3) = IBTAND(IHRSYM,4)
      END IF
C
C     MAX(ICOOR) = 3 + 6 = 9
C     MAX(IATOM) = 3*2 + (3 + 4 + 4 + 3 + 4 + 3) = 27
C     MAX(IHCINT) = 3*2*(1 + 3 + 3*2) + (3 + 4 + 4 + 3 + 4 + 3) = 81
C
      ICOOR = 0
      IATOM = 0
      IHCINT = 0
C
C     ***** First Derivatives and Cross Derivatives *****
C
      DO 10 ICOOR2 = 1, 3
         JDEO = IDEO(ICOOR2,1,IPATH)
         DCOOR0 = DEOVEC(JDEO + 1) .OR. DEOVEC(JDEO + 2)
         DCOOR = DCOOR0
         IF (CROSS) THEN
            DCORDH(1) = DMX(1,ICOOR2,1,IPATH)
     &             .OR. DMX(2,ICOOR2,1,IPATH)
     &             .OR. DMX(3,ICOOR2,1,IPATH)
            DCORDH(2) = DMX(1,ICOOR2,2,IPATH)
     &             .OR. DMX(2,ICOOR2,2,IPATH)
     &             .OR. DMX(3,ICOOR2,2,IPATH)
            DCORDE(1) = DET(1,ICOOR2,1,1,IPATH)
     &             .OR. DET(1,ICOOR2,2,1,IPATH)
     &             .OR. DET(2,ICOOR2,1,1,IPATH)
     &             .OR. DET(2,ICOOR2,2,1,IPATH)
     &             .OR. DET(3,ICOOR2,1,1,IPATH)
     &             .OR. DET(3,ICOOR2,2,1,IPATH)
            DCORDE(2) = DET(1,ICOOR2,1,2,IPATH)
     &             .OR. DET(1,ICOOR2,2,2,IPATH)
     &             .OR. DET(2,ICOOR2,1,2,IPATH)
     &             .OR. DET(2,ICOOR2,2,2,IPATH)
     &             .OR. DET(3,ICOOR2,1,2,IPATH)
     &             .OR. DET(3,ICOOR2,2,2,IPATH)
            DCOORH = DCORDH(1) .OR. DCORDH(2)
            DCOORE = DCORDE(1) .OR. DCORDE(2)
            DCOOR = DCOOR .OR. DCOORH .OR. DCOORE
         END IF
         IF (DCOOR) THEN
            ICOOR = ICOOR + 1
            JSTRAT(ICOOR) = IATOM + 1
            IWORK = IHCINT
            IX = 1 - IXDER(ICOOR2)
            IY = IYDER(ICOOR2)
            IZ = IZDER(ICOOR2)
            INCMXT(ICOOR) = IX
            INCMXU(ICOOR) = IY
            INCMXV(ICOOR) = IZ
            IEXOFF = IEOFF(IX)
            IEYOFF = IEOFF(IY)
            IEZOFF = IEOFF(IZ)
            JODD2  = IBTAND(IHRSYM,MOD(IX,2)+2*MOD(IY,2)+4*MOD(IZ,2))
            IATOM2 = 0
            DO 20 IATOMX = 0, IX
            DO 20 IATOMY = 0, IY
            DO 20 IATOMZ = 0, IZ
               IATOM2 = IATOM2 + 1
               IADREO = JDEO + IATOM2
               DATOM0 = DEOVEC(IADREO)
               DATOM = DATOM0
               IF (CROSS) THEN
                  DATOMH = DCORDH(IATOM2)
                  DATOME = DCORDE(IATOM2)
                  DATOM = DATOM .OR. DATOMH .OR. DATOME
               END IF
               IF (DATOM) THEN
                  IATOM = IATOM + 1
                  JSTRHC(IATOM,ICOOR) = IHCINT + 1
                  IOFFEX(IATOM) = IEXOFF + IATOMX
                  IOFFEY(IATOM) = IEYOFF + IATOMY
                  IOFFEZ(IATOM) = IEZOFF + IATOMZ
                  IF (DATOM0) THEN
                     IHCINT = IHCINT + 1
                     JODDIF(IHCINT) = JODD2
                     FACTOR(IHCINT) = FEOVEC(IADREO)
                     IOFFCC(IHCINT) = IEOVEC(IADREO)
                     IOFFHC(IHCINT) = IHCVEC(1)
                  END IF
                  IF (CROSS) THEN
                     DO 30 ICOOR1 = 1, 3
                        JODD = IBTXOR(JODD1(ICOOR1),JODD2)
                        IADRHC = 3*ICOOR1 - 1
                        IF (DATOMH) THEN
                           IADRMX = IDMX(ICOOR2,ICOOR1,IATOM2,IPATH)
                           IF (DMXVEC(IADRMX)) THEN
                              IHCINT = IHCINT + 1
                              JODDIF(IHCINT) = JODD
                              FACTOR(IHCINT) = FMXVEC(IADRMX)
                              IOFFCC(IHCINT) = JMXVEC(IADRMX)
                              IOFFHC(IHCINT) = IHCVEC(IADRHC)
                           END IF
                        END IF
                        IF (DATOME) THEN
                           DO 40 IATOM1 = 1, 2
                              IADRET = IDET(ICOOR2,ICOOR1,IATOM2,
     &                                      IATOM1,IPATH)
                              IF (DETVEC(IADRET)) THEN
                                 IHCINT = IHCINT + 1
                                 JODDIF(IHCINT)=JODD
                                 FACTOR(IHCINT)=FETVEC(IADRET)
                                 IOFFCC(IHCINT)=IETVEC(IADRET)
                                 IOFFHC(IHCINT)=IHCVEC(IADRHC+IATOM1)
                              END IF
   40                      CONTINUE
                        END IF
   30                CONTINUE
                  END IF
                 JENDHC(IATOM,ICOOR) = IHCINT
               END IF
   20       CONTINUE
            JENDAT(ICOOR) = IATOM
         END IF
   10 CONTINUE
C
C     ***** Second Derivatives *****
C
      IF (IDERIV .EQ. 2 .OR. SPIORB) THEN
         DO 50 ICOOR2 = 1, 6
            JDEO = IDEO(ICOOR2,2,IPATH)
            DCOOR = .FALSE.
            DO 55 IATOM2 = 1, IEND(ICOOR2)
               DCOOR = DCOOR .OR. DEOVEC(JDEO + IATOM2)
   55       CONTINUE
            IF (DCOOR) THEN
               ICOOR = ICOOR + 1
               JSTRAT(ICOOR) = IATOM + 1
               IWORK = IHCINT
               IX = 2 - IXDER(ICOOR2)
               IY = IYDER(ICOOR2)
               IZ = IZDER(ICOOR2)
               INCMXT(ICOOR) = IX
               INCMXU(ICOOR) = IY
               INCMXV(ICOOR) = IZ
               IEXOFF = IEOFF(IX)
               IEYOFF = IEOFF(IY)
               IEZOFF = IEOFF(IZ)
               JODD2  = IBTAND(IHRSYM,MOD(IX,2)+2*MOD(IY,2)+4*MOD(IZ,2))
               IATOM2 = 0
               DO 60 IATOMX = 0, IX
               DO 60 IATOMY = 0, IY
               DO 60 IATOMZ = 0, IZ
                  IATOM2 = IATOM2 + 1
                  IADREO = JDEO + IATOM2
                  IF (DEOVEC(IADREO)) THEN
                     IATOM = IATOM + 1
                     IHCINT = IHCINT + 1
                     JSTRHC(IATOM,ICOOR) = IHCINT
                     JENDHC(IATOM,ICOOR) = IHCINT
                     IOFFEX(IATOM) = IEXOFF + IATOMX
                     IOFFEY(IATOM) = IEYOFF + IATOMY
                     IOFFEZ(IATOM) = IEZOFF + IATOMZ
                     IF (SPIORB) THEN
                        IOFF2X(IATOM) = IEXOFF + IX - IATOMX
                        IOFF2Y(IATOM) = IEYOFF + IY - IATOMY
                        IOFF2Z(IATOM) = IEZOFF + IZ - IATOMZ
                     END IF
                     JODDIF(IHCINT) = JODD2
                     FACTOR(IHCINT) = FEOVEC(IADREO)
                     IOFFCC(IHCINT) = IEOVEC(IADREO)
                     IOFFHC(IHCINT) = IHCVEC(1)
                  END IF
   60          CONTINUE
               JENDAT(ICOOR) = IATOM
            END IF
   50    CONTINUE
      END IF
      NCOOR = ICOOR
C
      NHCMAX = 0
      DO 100 ITYPE = 1, 2
         IF (ITYPE .EQ. 1) THEN
            ICMPMX = KCKT34
         ELSE
            ICMPMX = KHKT34
         END IF
         DO 110 ICMP34 = 1, ICMPMX
            JODD34 = IODD34(ICMP34,ITYPE)
            DO 120 ICOOR = 1, NCOOR
               ISTART = JSTRHC(JSTRAT(ICOOR),ICOOR)
               ILAST  = JENDHC(JENDAT(ICOOR),ICOOR)
               NHCCMP = 0
               DO 130 ICMP12 = 1, KHKT12
                 IO1234 = IBTXOR(IODD12(ICMP12,2),JODD34)
                 DO 140 I = ISTART, ILAST
                    IF (IO1234 .EQ. JODDIF(I)) NHCCMP = NHCCMP + 1
  140            CONTINUE
  130          CONTINUE
               NHCC(ICMP34,ITYPE,ICOOR) = NHCCMP
               NHCMAX = MAX(NHCMAX,NHCCMP)
  120       CONTINUE
  110    CONTINUE
  100 CONTINUE
C
C     *************************
C     ***** PRINT SECTION *****
C     *************************
C
      IF (IPRINT .GE. 10) THEN
         WRITE (LUPRI, 2010)
         WRITE (LUPRI, 2020) PATH1
         WRITE (LUPRI, 2030) CROSS
         WRITE (LUPRI, 2040) IDERIV
         WRITE (LUPRI, 2050) KCKT12, KCKT34
         WRITE (LUPRI, 2060) NORB12, NORB34
         WRITE (LUPRI, 2070) NCOOR
         WRITE (LUPRI, 2080) IATOM
         WRITE (LUPRI, 2090) IHCINT
         WRITE (LUPRI, 2100) (INCMXT(I), I = 1, NCOOR)
         WRITE (LUPRI, 2110) (INCMXU(I), I = 1, NCOOR)
         WRITE (LUPRI, 2120) (INCMXV(I), I = 1, NCOOR)
         WRITE (LUPRI, 2130) (JSTRAT(I), I = 1, NCOOR)
         WRITE (LUPRI, 2140) (JENDAT(I), I = 1, NCOOR)
         WRITE (LUPRI, 2150)
         DO 80 I = 1, NCOOR
            WRITE (LUPRI, 2151) (JSTRHC(J,I),J=JSTRAT(I),JENDAT(I))
   80    CONTINUE
         WRITE (LUPRI, 2160)
         DO 90 I = 1, NCOOR
            WRITE (LUPRI, 2161) (JENDHC(J,I),J=JSTRAT(I),JENDAT(I))
   90    CONTINUE
         WRITE (LUPRI, 2170) (IOFFEX(I), I = 1, IATOM)
         WRITE (LUPRI, 2180) (IOFFEY(I), I = 1, IATOM)
         WRITE (LUPRI, 2190) (IOFFEZ(I), I = 1, IATOM)
         IF (SPIORB) THEN
            WRITE (LUPRI, 2175) (IOFF2X(I), I = 1, IATOM)
            WRITE (LUPRI, 2185) (IOFF2Y(I), I = 1, IATOM)
            WRITE (LUPRI, 2195) (IOFF2Z(I), I = 1, IATOM)
         END IF
         WRITE (LUPRI, 2200) (IOFFHC(I), I = 1, IHCINT)
         WRITE (LUPRI, 2210) (IOFFCC(I), I = 1, IHCINT)
         WRITE (LUPRI, 2220) (JODDIF(I), I = 1, IHCINT)
         WRITE (LUPRI, 2230) (FACTOR(I), I = 1, IHCINT)
         WRITE (LUPRI, 2260) IHRSYM
         WRITE (LUPRI, 2290) NUC34
         DO 95 I = 1, NCOOR
            WRITE (LUPRI,'(1X,A,I1,A,6I7/,(20X,6I7))')
     &         'NHCC 1 for ICOOR = ',I,':',(NHCC(J,1,I),J=1,KCKT34)
            WRITE (LUPRI,'(1X,A,I1,A,6I7/,(20X,6I7))')
     &         'NHCC 2 for ICOOR = ',I,':',(NHCC(J,2,I),J=1,KHKT34)
   95    CONTINUE
         WRITE (LUPRI, '(1X,A,I7)') 'NHCMAX ', NHCMAX
      END IF
 2010 FORMAT (   1X,'           INITIALIZATION ',/)
 2020 FORMAT(1X,'PATH1 ',L7)
 2030 FORMAT(1X,'CROSS ',L7)
 2040 FORMAT(1X,'IDERIV',I7)
 2050 FORMAT(1X,'KCKT  ',2I7)
 2060 FORMAT(1X,'NORB  ',2I7)
 2070 FORMAT(1X,'NCOOR ',I7)
 2080 FORMAT(1X,'IATOM ',I7)
 2090 FORMAT(1X,'IHCINT',I7)
 2100 FORMAT(1X,'INCMXT',9(I7))
 2110 FORMAT(1X,'INCMXU',9(I7))
 2120 FORMAT(1X,'INCMXV',9(I7))
 2130 FORMAT(1X,'JSTRAT',9(I7))
 2140 FORMAT(1X,'JENDAT',9(I7))
 2150 FORMAT(/,1X,'JSTRHC')
 2151 FORMAT(7X,10I7)
 2160 FORMAT(/,1X,'JENDHC')
 2161 FORMAT(7X,10I7)
 2170 FORMAT(/,1X,'IOFFEX',(10I6))
 2180 FORMAT(1X,'IOFFEY',(10I6))
 2190 FORMAT(1X,'IOFFEZ',(10I6))
 2175 FORMAT(/,1X,'IOFF2X',(10I6))
 2185 FORMAT(1X,'IOFF2Y',(10I6))
 2195 FORMAT(1X,'IOFF2Z',(10I6))
 2200 FORMAT(1X,'IOFFHC',(10I6))
 2210 FORMAT(1X,'IOFFCC',(10I6))
 2220 FORMAT(1X,'JODDIF',(10I6))
 2230 FORMAT(1X,'FACTOR',(10F6.2))
 2260 FORMAT(1X,'IHRSYM',I7)
 2290 FORMAT(1X,'NUC34 ',I7)
C
      IF (NHCMAX .GT. 0) THEN
C
C        Work space allocations in C2EIN1
C
C        SSINT | CCONT | CCPRIM | ETUV2 | ETUV1
C                               | SCR1  | SCR2
C                      | CSINT
C
         LSSINT = NO1234*NHCMAX*KHKT34
         LCCPRM = NCCPP*NHCMAX
         LETUV1 = NUC34
         LETUV2 = NCCPP
         IF (SPHR34) THEN
            LCCONT = NO1234*NHCMAX*KCKT34
         ELSE
            LCCONT = 0
         END IF
         IF (SPHR3 .AND. SPHR4) THEN
            LCSINT = NO1234*NHCMAX*KHKT4
         ELSE
            LCSINT = 0
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
            KCCPRM = KCCONT + LCCONT
C
               KETUV2 = KCCPRM + LCCPRM
               KETUV1 = KETUV2 + LETUV2
               KLAST1 = KETUV1 + LETUV1
C
               KSCR1  = KCCPRM + LCCPRM
               KSCR2  = KSCR1  + LSCR1
               KLAST2 = KSCR2  + LSCR2
C
            KCSINT = KCCONT + LCCONT
            KLAST3 = KCSINT + LCSINT
C
         KLAST = MAX(KLAST1,KLAST2,KLAST3)
C
         IF (KLAST .GT. LWORK) CALL STOPIT('C2EINT',' ',KLAST,LWORK)
         MWC2EI = MAX(MWC2EI,KLAST)
         LWTOT  = LWTOT + KLAST
         MWTOT  = MAX(MWTOT,LWTOT)
         CALL C2EIN1(AOINT,HCINT,COEF34,CONT3,CONT4,WORK(KCCPRM),NHCC,
     &               WORK(KCCONT),WORK(KSSINT),WORK(KCSINT),
     &               WORK(KSCR1),WORK(KSCR2),WORK(KETUV1),
     &               WORK(KETUV2),NCOOR,LMNV34,NPCO3,NPCO4,NUCS34,
     &               IODD12,IODD34,INDHER,NPNT34,NRED34,IPRINT,NHCMAX)
         LWTOT = LWTOT - KLAST
      END IF
      RETURN
      END
C  /* Deck c2ein1 */
      SUBROUTINE C2EIN1(AOINT,HCINT,COEF34,CONT3,CONT4,CCPRIM,NHCC,
     &                  CCONT,SSINT,CSINT,SCR1,SCR2,ETUV1,ETUV2,NCOOR,
     &                  LMNV34,NPCO3,NPCO4,NUCS34,IODD12,IODD34,INDHER,
     &                  NPNT34,NRED34,IPRINT,NHCMAX)
C
C     TUH
C
#include <implicit.h>
#include <priunit.h>
#include <maxaqn.h>
      PARAMETER (D0 = 0.0D0)
      INTEGER T, U, V
      COMMON /COMC2E/ JSTRHC(27,9), JENDHC(27,9),
     &                INCMXT(9), INCMXU(9), INCMXV(9),
     &                IOFFEX(27), IOFFEY(27), IOFFEZ(27),
     &                IOFF2X(27), IOFF2Y(27), IOFF2Z(27),
     &                IOFFHC(81), IOFFCC(81), JODDIF(81),
     &                FACTOR(81), JSTRAT(9), JENDAT(9)
      DIMENSION NHCC(KCKT34,2,9), LMNV34(KCKMX,5,2), CCPRIM(*),
     &          ETUV1(NUC34), ETUV2(NCCPP),HCINT(NCCPP,NTUV34,KHKT12,*),
     &          SSINT(NO1234*NHCMAX,KHKT34), CSINT(NO1234*NHCMAX,KHKT4),
     &          CCONT(NO1234*NHCMAX,KCKT34),
     &          AOINT(NOABCD,KHKTCD*KHKTAB,*),
     &          CONT3(*), CONT4(*),
     &          COEF34(MXUC34,0:JMAX3+JMAX4,0:JMAX3,0:JMAX4,3,*),
     &          IODD12(KCKT12,2), IODD34(KCKT34,2), SCR1(*), SCR2(*),
     &          NPCO3(*), NPCO4(*), NUCS34(*), NRED34(*),
     &          INDHER(0:JTOP,0:JTOP,0:JTOP), NPNT34(*)
#include <hertop.h>
#include <twoao.h>
#include <twocom.h>
#include <sphtrm.h>
C
#include <ibtfun.h>
C
      IF (IPRINT .GT. 5) CALL TITLER('Output from C2EIN1','*',103)
C
      INCRMT = I340X + 1
      INCRMU = I340Y + 1
      INCRMV = I340Z + 1
C
      DO 100 ICOOR = 1, NCOOR
         ISTRAT = JSTRAT(ICOOR)
         IENDAT = JENDAT(ICOOR)
C
C        *******************************
C        ***** Cartesian integrals *****
C        *******************************
C
         ICMP34 = 0
         DO 200 ICOMP3 = 1, KCKT3
            L3 = LMNV34(ICOMP3,1,1)
            M3 = LMNV34(ICOMP3,2,1)
            N3 = LMNV34(ICOMP3,3,1)
            MAX4 = KCKT4
            IF (DIAC34) MAX4 = ICOMP3
            DO 210 ICOMP4 = 1, MAX4
               ICMP34 = ICMP34 + 1
               NHCCMP = NHCC(ICMP34,1,ICOOR)
               IF (NHCCMP .GT. 0) THEN
                  L4 = LMNV34(ICOMP4,1,2)
                  M4 = LMNV34(ICOMP4,2,2)
                  N4 = LMNV34(ICOMP4,3,2)
                  JODD34 = IODD34(ICMP34,1)
C
C                 Primitive Cartesian Integrals
C                 =============================
C
                  CALL DZERO(CCPRIM,NCCPP*NHCCMP)
                  MAXT   = L3 + L4 + INCMXT(ICOOR)
                  MAXU   = M3 + M4 + INCMXU(ICOOR)
                  MAXV   = N3 + N4 + INCMXV(ICOOR)
                  MAXTUV = MAXT + MAXU + MAXV - 1
                  MINT   = IBTAND(MAXT,INCRMT - 1)
                  MINU   = IBTAND(MAXU,INCRMU - 1)
                  MINV   = IBTAND(MAXV,INCRMV - 1)
                  DO 300 V = MINV, MAXV, INCRMV
                  DO 300 U = MINU, MAXU, INCRMU
                  DO 300 T = MINT, MAXT, INCRMT
                     ITUV  = INDHER(T,U,V)
                     IWORK = 0
                     DO 400 IATOM = ISTRAT, IENDAT
                        ISTRHC = JSTRHC(IATOM,ICOOR)
                        IENDHC = JENDHC(IATOM,ICOOR)
                        IODD   = IBTXOR(JODDIF(ISTRHC),JODD34)
C
C                       Expansion coefficients
C                       ======================
C
                        IF (SPIORB) THEN
                           IF (T + U + V .GT. MAXTUV) GO TO 300
                           I1X = IOFFEX(IATOM)
                           I1Y = IOFFEY(IATOM)
                           I1Z = IOFFEZ(IATOM)
                           I2X = IOFF2X(IATOM)
                           I2Y = IOFF2Y(IATOM)
                           I2Z = IOFF2Z(IATOM)
                           DO 500 I = 1, NUC34
                              ETUV1(I) = ( COEF34(I,T,L3,L4,1,I1X)
     &                                    *COEF34(I,U,M3,M4,2,I1Y)
     &                                    *COEF34(I,V,N3,N4,3,I1Z))
     &                                  -( COEF34(I,T,L3,L4,1,I2X)
     &                                    *COEF34(I,U,M3,M4,2,I2Y)
     &                                    *COEF34(I,V,N3,N4,3,I2Z))
  500                      CONTINUE
                        ELSE
                           IEX = IOFFEX(IATOM)
                           IEY = IOFFEY(IATOM)
                           IEZ = IOFFEZ(IATOM)
                           DO 510 I = 1, NUC34
                              ETUV1(I) = COEF34(I,T,L3,L4,1,IEX)
     &                                 * COEF34(I,U,M3,M4,2,IEY)
     &                                 * COEF34(I,V,N3,N4,3,IEZ)
  510                      CONTINUE
                        END IF
C
C                       Primitive Cartesian integrals for KHKT12 = 1
C                       ============================================
C
                        IF ((KHKT12.EQ.1) .AND. (ISTRHC.EQ.IENDHC)) THEN
                           IF (IODD12(1,2) .EQ. IODD) THEN
                              INTHC = IOFFHC(ISTRHC)
                              IJ = 0
                              DO 520 I = 1, NORB12
                                 DO 530 J = 1, NUC34
                                    CCPRIM(IWORK+J) = CCPRIM(IWORK+J)
     &                               + ETUV1(J)*HCINT(IJ+J,ITUV,1,INTHC)
  530                            CONTINUE
                                 IJ = IJ + NUC34
                                 IWORK = IWORK + NUC34
  520                         CONTINUE
                           END IF
C
C                       Primitive Cartesian integrals for KHKT12 > 1
C                       ============================================
C
                        ELSE
                           IJ = 0
                           DO 540 I = 1, NORB12
                           DO 540 J = 1, NUC34
                              IJ = IJ + 1
                              ETUV2(IJ) = ETUV1(J)
  540                      CONTINUE
                           DO 550 IHCINT = ISTRHC, IENDHC
                              INTHC  = IOFFHC(IHCINT)
                              IODDIF = IBTXOR(JODDIF(IHCINT),JODD34)
                              DO 560 I = 1,KHKT12
                              IF (IODD12(I,2) .EQ. IODDIF) THEN
                                 DO 570 J = 1, NCCPP
                                    CCPRIM(IWORK+J) = CCPRIM(IWORK+J)
     &                               + ETUV2(J)*HCINT(J,ITUV,I,INTHC)
  570                            CONTINUE
                                 IWORK = IWORK + NCCPP
                              END IF
  560                         CONTINUE
  550                      CONTINUE
                        END IF
  400                CONTINUE
  300             CONTINUE
C
C                 Contracted Integrals
C                 ====================
C
                  IF (SPHR34) THEN
                     CALL C2CONT(CCPRIM,CCONT(1,ICMP34),CONT3,CONT4,
     &                           SCR1,SCR2,NPCO3,NPCO4,NUCS34,NHCCMP,
     &                           NPNT34,NRED34,IPRINT,.FALSE.)
                  ELSE
                     CALL C2CONT(CCPRIM,SSINT(1,ICMP34),CONT3,CONT4,
     &                           SCR1,SCR2,NPCO3,NPCO4,NUCS34,NHCCMP,
     &                           NPNT34,NRED34,IPRINT,.FALSE.)
                  END IF
               END IF
  210       CONTINUE
  200    CONTINUE
C
C        *******************************
C        ***** Spherical integrals *****
C        *******************************
C
         IF (SPHR34) THEN
            CALL C2SPHR(CCONT,CSINT,SSINT,NHCC(1,1,ICOOR),
     &                  CSP(ISPADR(NHKT3)),CSP(ISPADR(NHKT4)),
     &                  NHCMAX,IPRINT)
         END IF
C
C        *****************************************
C        ***** Multiply by factors and order *****
C        *****************************************
C
C        Swap indices 12 and 34 for path 2
C
         DO 700 ICMP34 = 1, KHKT34
            JODD34 = IODD34(ICMP34,2)
            IWORK = 0
            IF (PATH1) THEN
               DO 800 IATOM = ISTRAT, IENDAT
                  ISTRHC = JSTRHC(IATOM,ICOOR)
                  IENDHC = JENDHC(IATOM,ICOOR)
                  DO 810 IHCINT = ISTRHC, IENDHC
                     ITYPE  = IOFFCC(IHCINT)
                     ICMP   = ICMP34
                     FAC    = FACTOR(IHCINT)
                     IODDIF = IBTXOR(JODDIF(IHCINT),JODD34)
                     DO 820 ICMP12 = 1, KHKT12
                        IF (IODD12(ICMP12,2) .EQ. IODDIF) THEN
                           DO 830 I = 1, NO1234
                             AOINT(I,ICMP,ITYPE) = AOINT(I,ICMP,ITYPE)
     &                                   + FAC*SSINT(IWORK + I,ICMP34)
  830                      CONTINUE
                           IWORK = IWORK + NO1234
                        END IF
                        ICMP = ICMP + KHKT34
  820                CONTINUE
  810             CONTINUE
  800          CONTINUE
            ELSE
               ICMP0 = (ICMP34 - 1)*KHKT12
               DO 900 IATOM = ISTRAT, IENDAT
                  ISTRHC = JSTRHC(IATOM,ICOOR)
                  IENDHC = JENDHC(IATOM,ICOOR)
                  DO 910 IHCINT = ISTRHC, IENDHC
                     ITYPE  = IOFFCC(IHCINT)
                     FAC    = FACTOR(IHCINT)
                     IODDIF = IBTXOR(JODDIF(IHCINT),JODD34)
                     DO 920 ICMP12 = 1, KHKT12
                     IF (IODD12(ICMP12,2) .EQ. IODDIF) THEN
                        ICMP = ICMP0 + ICMP12
                        DO 930 I = 1, NORB12
                           IJ = I
                           DO 940 J = 1, NORB34
                              AOINT(IJ,ICMP,ITYPE)=AOINT(IJ,ICMP,ITYPE)
     &                                      + FAC*SSINT(IWORK+J,ICMP34)
                              IJ = IJ + NORB12
  940                      CONTINUE
                           IWORK  = IWORK  + NORB34
  930                   CONTINUE
                     END IF
  920                CONTINUE
  910             CONTINUE
  900          CONTINUE
            END IF
            IF (IPRINT .GE. 10) THEN
               WRITE (LUPRI,'(1X,A, I7 )') 'ICMP34',ICMP34
               WRITE (LUPRI,'(1X,A,3I7)') ' MINT  ',MINT, MINU, MINV
               WRITE (LUPRI,'(1X,A,3I7)') ' MAXT  ',MAXT, MAXU, MAXV
               WRITE (LUPRI,'(1X,A,3I7)') ' INCRMT',INCRMT,INCRMU,INCRMV
               IF (IPRINT .GE. 15) THEN
                  CALL HEADER('Final Cartesian Integrals - C2EINT',-1)
                  WRITE (LUPRI,'(/1X,A,I3)') ' ICMP34 ', ICMP34
                  DO 1000 IATOM = ISTRAT, IENDAT
                     WRITE (LUPRI,'(/1X,A,I3)') '  IATOM ', IATOM
                     DO 1010 I = JSTRHC(IATOM,ICOOR),JENDHC(IATOM,ICOOR)
                        WRITE (LUPRI,'(/1X,A,I3)') ' IHCINT ', I
                        ITYPE = IOFFCC(I)
                        IF (PATH1) THEN
                           IADR = ICMP34
                           IADD = KHKT34
                        ELSE
                           IADR = (ICMP34-1)*KHKT12 + 1
                           IADD = 1
                        END IF
                        IODDIF = IBTXOR(JODDIF(ICOOR),JODD34)
                        DO 1020 ICMP12 = 1, KHKT12
                           WRITE (LUPRI,'(/1X,A,I3)') ' ICMP12 ', ICMP12
                           IF (IODD12(ICMP12,2).EQ. IODDIF) THEN
                              WRITE (LUPRI,'(1P,6D12.4)')
     &                              (AOINT(K,IADR,ITYPE),K=1,NO1234)
                           ELSE
                              WRITE (LUPRI,'(A)')
     &                           ' Integrals zero by symmetry'
                           END IF
                           IADR = IADR + IADD
 1020                   CONTINUE
 1010                CONTINUE
 1000             CONTINUE
               END IF
            END IF
  700    CONTINUE
  100 CONTINUE
      RETURN
      END
C  /* Deck c2cont */
      SUBROUTINE C2CONT(CCPRIM,CCONT,CONT3,CONT4,SCR1,SCR2,NPCO3,NPCO4,
     &                  NUCS34,NHCCMP,NPNT34,NRED34,IPRINT,ANTICO)
#include <implicit.h>
#include <priunit.h>
C
C     tuh Mar 1988
C     Modified for no transformations (NOCNT) tuh Apr 1989
C     Modified for segmented contraction tuh Feb 1992
C     Modified for antihermitean integrals (WK/UniKA/18-11-2002).
C
C     Purpose: Transformation of two innermost indices
C              Third index is next innermost
C              Fourth index is innermost
C
C     Note: After transformation order of two innermost indices
C           is reversed
C
C     In:  CCPRIM(NINT12*NUC4*NUC3)
C          COEF3(NORR4*NUC4)
C          COEF4(NORR3*NUC3)
C
C     Out: CCONT(NINT12*NORR3*NORR4)
C
C     Scratch: SCR1(NUCR3*NUCR4*NINT12)
C              SCR2(NORR3*NUCR4*NINT12)
C
C     DM1 and THRESH are added for antihermitean 
C     integrals (WK/UniKA/18-11-2002).
      PARAMETER (D0 = 0.0D0, D1 = 1.0D0, DM1 = -1.0D0, THRESH = 1.0D-13)
#include <twoao.h>
      LOGICAL ANTICO
      DIMENSION CCPRIM(*), CCONT(*), SCR1(*), SCR2(*),
     &          CONT3(*), CONT4(*), NPCO3(NSET3,2), NPCO4(NSET4,2),
     &          NUCS34(*), NPNT34(NUC3*NUC4,2), NRED34(*)
C
      IF (IPRINT .GT. 20) THEN
         CALL HEADER('Uncontracted Cartesian integrals in C2CONT',-1)
         WRITE (LUPRI,'(/7X,A,I3)')'Number of columns: NHCCMP =',NHCCMP
         WRITE (LUPRI,'(7X,A,I3/)')'Number of rows:    NCCPP  =',NCCPP
         CALL OUTPUT(CCPRIM,1,NHCCMP,1,NCCPP,NHCCMP,NCCPP,1,LUPRI)
      END IF
C
      NINT12 = NHCCMP*NORB12
C
C     Case (i): No contraction
C     ========================
C
      IF (NOCNT) THEN
#ifndef OLD_NOCONT
C980827-hjaaj: now the NOCNT flag is controlled in CCDRIV
C              and it checks if a segmented contr. is uncontracted
         CALL DCOPY(NINT12*NORB34,CCPRIM,1,CCONT,1)
#else
C980827-hjaaj: old .NOCONT did not work, I have made
C              TWOINT stop if NOCONT is set.
         IF (TPRI34) THEN
            N34D = NUC4*(NUC4 + 1)/2
            CALL DCOPY(NINT12*N34D,CCPRIM,1,CCONT,1)
         ELSE
            IJ = 1
            N34Q = NUC4*NUC3
            DO 10 I = 1, NUC4
               DO 20 J = 1, NUC3
                  JI = (J - 1)*NUC4 + I
                  CALL DCOPY(NINT12,CCPRIM(IJ),N34Q,CCONT(JI),N34Q)
                  IJ = IJ + 1
   20          CONTINUE
   10       CONTINUE
         END IF
#endif
C
C     Case (ii): General contraction
C     ==============================
C
      ELSE IF (GEN34) THEN
C
C
C        Transform fourth index
C
         IF (TPRI34) THEN
            N34Q = NUCR3*NUCR4
            IF (RPRI34) THEN
               CALL DZERO(SCR1,N34Q*NINT12)
               DO 100 I = 1, NUC34
                  IJ = NPNT34(I,1)
                  JI = NPNT34(I,2)
                  CALL DCOPY(NINT12,CCPRIM(I),NUC34,SCR1(IJ),N34Q)
                  IF (IJ .NE. JI) THEN
C                    Antihermitean integrals (WK/UniKA/18-11-2002).
                     IF (ANTICO) THEN
                        CALL VSMUL(CCPRIM(I),NUC34,DM1,
     &                             SCR1(JI),N34Q,NINT12)
                     ELSE
                        CALL DCOPY(NINT12,CCPRIM(I),
     &                             NUC34,SCR1(JI),N34Q)
                     END IF
                  ELSE IF (ANTICO) THEN
                     IADS = I
                     DO 110 INT12 = 1, NINT12
                        IF (ABS(CCPRIM(IADS)) .GT. THRESH) THEN
                           WRITE (LUPRI,'(A,2(1P,E14.6))')
     &                     ' MATRIX NOT ANTI-SYMM. (C2CONT) ',
     &                     CCPRIM(IADS), THRESH
                           CALL QUIT('*** ERROR IN C2CONT ***')
                        END IF
                        IADS = IADS + NUC34
  110                CONTINUE
                  END IF
  100          CONTINUE
            ELSE
               N34D = NUC3*(NUC3 + 1)/2
               IADR = 1
               DO 200 I = 1, NUC4
               DO 200 J = 1, I
                  IJ = (I - 1)*NUC4 + J
                  JI = (J - 1)*NUC4 + I
                  CALL DCOPY(NINT12,CCPRIM(IADR),N34D,SCR1(IJ),N34Q)
                  IF (IJ .NE. JI) THEN
                     IF (ANTICO) THEN
                        CALL VSMUL(CCPRIM(IADR),N34D,DM1,
     &                             SCR1(JI),N34Q,NINT12)
                     ELSE
                        CALL DCOPY(NINT12,CCPRIM(IADR),
     &                             N34D,SCR1(JI),N34Q)
                     END IF
                  ELSE IF (ANTICO) THEN
                     IADS = IADR
                     DO 210 INT12 = 1, NINT12
                        IF (ABS(CCPRIM(IADS)) .GT. THRESH) THEN
                           WRITE (LUPRI,'(A,2(1P,E14.6))')
     &                     ' MATRIX NOT ANTI-SYMM. (C2CONT) ',
     &                     CCPRIM(IADS), THRESH
                           CALL QUIT('*** ERROR IN C2CONT ***')
                        END IF
                        IADS = IADS + N34D
  210                CONTINUE
                  END IF
                  IADR = IADR + 1
  200          CONTINUE
            END IF
            CALL DGEMM('N','N',NORR3,NUCR4*NINT12,NUCR3,D1,CONT3,NORR3,
     &           SCR1,NUCR3,D0,SCR2,NORR3)
         ELSE
            IF (RPRI34) THEN
               N34Q = NUCR3*NUCR4
               CALL DZERO(SCR1,N34Q*NINT12)
               DO 300 I = 1, NUC34
                  CALL DCOPY(NINT12,CCPRIM(I),NUC34,
     &                              SCR1(NPNT34(I,1)),N34Q)
  300          CONTINUE
               CALL DGEMM('N','N',NORR3,NUCR4*NINT12,NUCR3,D1,CONT3,
     &              NORR3,SCR1,NUCR3,D0,SCR2,NORR3)
            ELSE
               CALL DGEMM('N','N',NORR3,NUCR4*NINT12,NUCR3,D1,CONT3,
     &              NORR3,CCPRIM,NUCR3,D0,SCR2,NORR3)
            END IF
         END IF
C
C        Change order of third and fourth indices
C
         NPR34 = NUCR4*NORR3
         IADR10 = 0
         IADR20 = 0
         DO 400 I = 1, NUCR4
            DO 410 J = 1, NORR3
               IADR1 = IADR10 + J - 1
               IADR2 = IADR20 + (J - 1)*NUCR4
               DO 420 K = 1, (NINT12 - 1)*NPR34 + 1, NPR34
                  SCR1(IADR2 + K) = SCR2(IADR1 + K)
  420          CONTINUE
  410       CONTINUE
            IADR10 = IADR10 + NORR3
            IADR20 = IADR20 + 1
  400    CONTINUE
C
C        Transform third index
C
         NCT412 = NORR3*NINT12
         IF (TCON34 .OR. RCNT34) THEN
            CALL DGEMM('N','N',NORR4,NCT412,NUCR4,D1,CONT4,NORR4,SCR1,
     &           NUCR4,D0,SCR2,NORR4)
            NCTSQ = NORR3*NORR4
            DO 800 I = 1, NORB34
               CALL DCOPY(NINT12,SCR2(NRED34(I)),NCTSQ,CCONT(I),NORB34)
  800       CONTINUE
         ELSE
            CALL DGEMM('N','N',NORR4,NCT412,NUCR4,D1,CONT4,NORR4,SCR1,
     &           NUCR4,D0,CCONT,NORR4)
         END IF
C
C     Case (iii): Segmented contraction
C     =================================
C
      ELSE
         IF (ANTICO .AND. TPRI34 .AND. NORB34 .GT. 1) THEN
            CALL QUIT('*** NOT IMPLEMENTED (C2CONT) ***')
         ELSE IF (ANTICO .AND. TPRI34) THEN
            DO 930 I = 1, NINT12
               CCONT(I) = D0
  930       CONTINUE
         ELSE
            IOFFCN = 1
            IOFFPR = 0
            DO 900 I = 1, NINT12
               DO 910 IORB34 = 1, NORB34
                  NPRIM = NUCS34(IORB34)
                  IF (NPRIM .GT. 0) THEN
                     CNT = CCPRIM(IOFFPR+1)
                     DO 920 J = 2, NPRIM
                        CNT = CNT + CCPRIM(IOFFPR+J)
  920                CONTINUE
                     CCONT(IOFFCN) = CNT
                     IOFFPR = IOFFPR + NPRIM
                  ELSE
                     CCONT(IOFFCN) = D0
                  END IF
                  IOFFCN = IOFFCN + 1
  910          CONTINUE
  900       CONTINUE
         END IF
      END IF
      RETURN
      END
C  /* Deck c2sphr */
      SUBROUTINE C2SPHR(CCONT,CSINT,SSINT,NHCC,CSP3,CSP4,NHCMAX,IPRINT)
C
#include <implicit.h>
#include <priunit.h>
      PARAMETER (D0 = 0.0D0)
      DIMENSION CCONT(NO1234*NHCMAX,KCKT34),
     &          CSINT(NO1234*NHCMAX,KHKT4),
     &          SSINT(NO1234*NHCMAX,KHKT34),
     &          CSP3(KHKT3,KCKT3), CSP4(KHKT4,KCKT4), NHCC(KCKT34,2)
#include <twoao.h>
C
      CALL DZERO(SSINT,NO1234*NHCMAX*KHKT34)
C
C     Transformation of both indices
C     ==============================
C
      IF (SPHR3 .AND. SPHR4) THEN
         ICMP34 = 0
         DO 100 ICOMP3 = 1,KCKT3
C
C           First half transformation
C
            CALL DZERO(CSINT,NO1234*NHCMAX*KHKT4)
            DO 200 ICOMP4 = 1, KCKT4
               ICMP34 = ICMP34 + 1
               NHCCMP = NHCC(ICMP34,1)
               IF (NHCCMP .GT. 0) THEN
                  DO 210 IKOMP4 = 1, KHKT4
                     SPHFAC = CSP4(IKOMP4,ICOMP4)
                     IF (ABS(SPHFAC).GT.D0) THEN
                        DO 220 I = 1, NHCCMP*NO1234
                           CSINT(I,IKOMP4) = CSINT(I,IKOMP4)
     &                                     + SPHFAC*CCONT(I,ICMP34)
  220                   CONTINUE
                     END IF
  210             CONTINUE
               END IF
  200       CONTINUE
C
C           Second half transformation
C
            IKMP34 = 0
            DO 300 IKOMP3 = 1, KHKT3
               SPHFAC = CSP3(IKOMP3,ICOMP3)
               IF (ABS(SPHFAC) .GT. D0) THEN
                  MAX4 = KHKT4
                  IF (DIAG34) MAX4 = IKOMP3
                  DO 310 IKOMP4 = 1, MAX4
                     IKMP34 = IKMP34 + 1
                     NHCCMP = NHCC(IKMP34,2)
                     DO 320 I = 1, NHCCMP*NO1234
                        SSINT(I,IKMP34) = SSINT(I,IKMP34)
     &                           + SPHFAC*CSINT(I,IKOMP4)
  320                CONTINUE
  310             CONTINUE
               ELSE IF (DIAG34) THEN
                  IKMP34 = IKMP34 + IKOMP3
               ELSE
                  IKMP34 = IKMP34 + KHKT4
               END IF
  300       CONTINUE
  100    CONTINUE
C
C     Transformation of first index only
C     ==================================
C
      ELSE IF (SPHR3) THEN
         DO 400 ICOMP3 = 1, KCKT3
         DO 400 IKOMP3 = 1, KHKT3
            SPHFAC = CSP3(IKOMP3,ICOMP3)
            IF (ABS(SPHFAC) .GT. D0) THEN
               DO 420 IKOMP4 = 1, KHKT4
                  ICMP34 = (ICOMP3 - 1)*KHKT4 + IKOMP4
                  IKMP34 = (IKOMP3 - 1)*KHKT4 + IKOMP4
                  NHCCMP = NHCC(IKMP34,2)
                  DO 430 I = 1, NHCCMP*NO1234
                     SSINT(I,IKMP34) = SSINT(I,IKMP34)
     &                        + SPHFAC*CCONT(I,ICMP34)
  430             CONTINUE
  420          CONTINUE
            END IF
  400    CONTINUE
C
C     Transformation of second index only
C     ===================================
C
      ELSE
         DO 500 IKOMP3 = 1, KHKT3
         DO 500 ICOMP4 = 1, KCKT4
            ICMP34 = (IKOMP3 - 1)*KCKT4 + ICOMP4
            NHCCMP = NHCC(ICMP34,1)
            IF (NHCCMP .GT. 0) THEN
               DO 510 IKOMP4 = 1, KHKT4
                  SPHFAC = CSP4(IKOMP4,ICOMP4)
                  IF (ABS(SPHFAC).GT.D0) THEN
                     IKMP34 = (IKOMP3 - 1)*KHKT4 + IKOMP4
                     DO 520 I = 1, NHCCMP*NO1234
                        SSINT(I,IKMP34) = SSINT(I,IKMP34)
     &                           + SPHFAC*CCONT(I,ICMP34)
  520                CONTINUE
                  END IF
  510          CONTINUE
            END IF
  500    CONTINUE
      END IF
      RETURN
      END
